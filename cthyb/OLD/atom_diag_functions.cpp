#include "./atom_diag_functions.hpp"

namespace cthyb {

  double const threshold = 1.e-11;

  double partition_function(atom_diag const &atom, double beta) {
    double z = 0;
    for (auto const &es : atom.get_eigensystem())
      for (auto e : es.eigenvalues) z += std::exp(-beta * e);
    return z;
  }

  // -----------------------------------------------------------------

  block_matrix_t atomic_density_matrix(atom_diag const &atom, double beta) {
    double z     = partition_function(atom, beta);
    int n_blocks = atom.n_blocks();
    block_matrix_t dm(n_blocks);
    for (int bl = 0; bl < n_blocks; ++bl) {
      int bl_size = atom.get_block_dim(bl);
      dm[bl]      = matrix<h_scalar_t>(bl_size, bl_size);
      for (int u = 0; u < bl_size; ++u) {
        for (int v = 0; v < bl_size; ++v) { dm[bl](u, v) = (u == v) ? std::exp(-beta * atom.get_eigenvalue(bl, u)) / z : 0; }
      }
    }
    return dm;
  }

  // -----------------------------------------------------------------

  block_gf<imtime> atomic_gf(atom_diag const &atom, double beta, std::map<std::string, indices_t> const &gf_struct, int n_tau,
                             std::vector<std::pair<int, int>> const &excluded_states) {

    double z = partition_function(atom, beta);

    std::vector<std::string> block_names;
    std::vector<gf<imtime>> gf_blocks;
    auto const &fops = atom.get_fops();

    auto is_excluded = [&excluded_states](int A, int ia) {
      return std::find(excluded_states.begin(), excluded_states.end(), std::make_pair(A, ia)) != excluded_states.end();
    };

    for (auto const &block : gf_struct) {
      block_names.push_back(block.first);
      int bl_size = block.second.size();
      auto g      = gf<imtime>{{beta, Fermion, n_tau}, {bl_size, bl_size}};

      for (int inner_index1 = 0; inner_index1 < bl_size; ++inner_index1)
        for (int inner_index2 = 0; inner_index2 < bl_size; ++inner_index2) {
          int n1 = fops[{block.first, block.second[inner_index1]}]; // linear_index of c
          int n2 = fops[{block.first, block.second[inner_index2]}]; // linear_index of c_dag

          for (int A = 0; A < atom.n_blocks(); ++A) {    // index of the A block. sum over all
            int B = atom.cdag_connection(n2, A);         // index of the block connected to A by operator c_n
            if (B == -1) continue;                       // no matrix element
            if (atom.c_connection(n1, B) != A) continue; //
            for (int ia = 0; ia < atom.get_block_dim(A); ++ia)
              for (int ib = 0; ib < atom.get_block_dim(B); ++ib) {
                auto Ea = atom.get_eigenvalue(A, ia);
                auto Eb = atom.get_eigenvalue(B, ib);
                if (is_excluded(A, ia) || is_excluded(B, ib)) continue;
                for (auto tau : g.mesh())
                  g[tau](inner_index1, inner_index2) +=
                     -atom.cdag_matrix(n2, A)(ib, ia) * atom.c_matrix(n1, B)(ia, ib) * std::exp(-(Eb - Ea) * tau - beta * Ea) / z;
              }
          }
        }
      g.singularity()(1) = 1.0;
      gf_blocks.push_back(std::move(g));
    }

    return make_block_gf(block_names, gf_blocks);
  }

  // -----------------------------------------------------------------

  double trace_rho_op(block_matrix_t const &density_matrix, many_body_op_t const &op, atom_diag const &atom) {
    h_scalar_t result = 0;
    if (atom.n_blocks() != density_matrix.size()) TRIQS_RUNTIME_ERROR << "trace_rho_op : size mismatch : number of blocks differ";
    for (int bl = 0; bl < atom.n_blocks(); ++bl) {
      if (atom.get_block_dim(bl) != first_dim(density_matrix[bl]))
        TRIQS_RUNTIME_ERROR << "trace_rho_op : size mismatch : size of block " << bl << " differ";
      for (auto const &x : op) {
        auto b_m = atom.matrix_element_of_monomial(x.monomial, bl);
        if (b_m.first == bl) result += x.coef * dot_product(b_m.second, density_matrix[bl]);
      }
    }
    if (imag(result) > threshold) TRIQS_RUNTIME_ERROR << "trace_rho_op: the result is not real.";
    return real(result);
  }

  // -----------------------------------------------------------------

  full_hilbert_space_state_t act(many_body_op_t const &op, full_hilbert_space_state_t const &st, atom_diag const &atom) {
    full_hilbert_space_state_t result(st.size());
    result() = 0;
    for (auto const &x : op) {
      for (int bl = 0; bl < atom.n_blocks(); ++bl) {
        auto b_m = atom.matrix_element_of_monomial(x.monomial, bl);
        if (b_m.first == -1) continue;
        result(atom.index_range_of_block(b_m.first)) += x.coef * b_m.second * st(atom.index_range_of_block(bl));
      }
    }
    return result;
  }

  //---------------------

  std::vector<std::vector<quantum_number_t>> quantum_number_eigenvalues(many_body_op_t const &op, atom_diag const &atom) {

    auto commutator = op * atom.get_h_atomic() - atom.get_h_atomic() * op;
    if (!commutator.is_zero()) TRIQS_RUNTIME_ERROR << "The operator is not a quantum number";

    std::vector<std::vector<quantum_number_t>> result;

    for (int bl = 0; bl < atom.n_blocks(); ++bl) {
      auto dim = atom.get_block_dim(bl);
      result.push_back(std::vector<quantum_number_t>(dim, 0));
      for (auto const &x : op) {
        auto b_m = atom.matrix_element_of_monomial(x.monomial, bl);
        if (b_m.first != bl) continue;
        for (int i = 0; i < dim; ++i)
          result.back()[i] += real(x.coef * b_m.second(i, i)); //FIXME leave general, cast into quantum_number_t before returning -- is this possible?
      }
    }
    return result;
  }

  //---------------------

  template <typename M>
  // require (ImmutableMatrix<M>)
  bool is_diagonal(M const &m) {
    return ((sum(abs(m)) - trace(abs(m))) < threshold);
  }

  std::vector<std::vector<quantum_number_t>> quantum_number_eigenvalues2(many_body_op_t const &op, atom_diag const &atom) {

    auto commutator = op * atom.get_h_atomic() - atom.get_h_atomic() * op;
    if (!commutator.is_zero()) TRIQS_RUNTIME_ERROR << "The operator is not a quantum number";

    auto d = atom.get_full_hilbert_space_dim();
    matrix<quantum_number_t> M(d, d);
    std::vector<std::vector<quantum_number_t>> result;

    for (int bl = 0; bl < atom.n_blocks(); ++bl) {
      auto dim = atom.get_block_dim(bl);
      for (auto const &x : op) {
        auto b_m = atom.matrix_element_of_monomial(x.monomial, bl);
        if (b_m.first == -1) continue;
        M(atom.index_range_of_block(b_m.first), atom.index_range_of_block(bl)) += real(x.coef * b_m.second);
      }
    }
    if (!is_diagonal(M)) TRIQS_RUNTIME_ERROR << "The Matrix of the operator is not diagonal !!!";

    for (int bl = 0; bl < atom.n_blocks(); ++bl) {
      auto dim = atom.get_block_dim(bl);
      result.push_back(std::vector<quantum_number_t>(dim, 0));
      for (int i = 0; i < dim; ++i) result.back()[i] = M(atom.flatten_block_index(bl, i), atom.flatten_block_index(bl, i));
    }

    return result;
  }
}
