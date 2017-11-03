/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2016, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#include "../functions.hpp"
#include <cmath>
#include <algorithm>
#include <limits>
#include <triqs/arrays.hpp>
#include <triqs/utility/numeric_ops.hpp>
#include <triqs/utility/legendre.hpp>

namespace triqs {
namespace atom_diag {

using namespace triqs::arrays;

#define ATOM_DIAG atom_diag<Complex>
#define ATOM_DIAG_R atom_diag<false>
#define ATOM_DIAG_C atom_diag<true>
#define ATOM_DIAG_T typename atom_diag<Complex>

template<bool Complex>
double partition_function(ATOM_DIAG const& atom, double beta) {
 double z = 0;
 for (auto const& es : atom.get_eigensystems()) z += sum(exp(-beta * es.eigenvalues));
 return z;
}
template double partition_function(ATOM_DIAG_R const&, double);
template double partition_function(ATOM_DIAG_C const&, double);

// -----------------------------------------------------------------
template<bool Complex>
ATOM_DIAG_T::block_matrix_t atomic_density_matrix(ATOM_DIAG const& atom, double beta) {
 double z = partition_function(atom, beta);
 int n_blocks = atom.n_subspaces();
 ATOM_DIAG_T::block_matrix_t dm(n_blocks);
 for (int bl = 0; bl < n_blocks; ++bl) {
  int bl_size = atom.get_subspace_dim(bl);
  dm[bl] = ATOM_DIAG_T::matrix_t(bl_size, bl_size);
  for (int u = 0; u < bl_size; ++u) {
   for (int v = 0; v < bl_size; ++v) {
    dm[bl](u, v) = (u == v) ? std::exp(-beta * atom.get_eigenvalue(bl, u)) / z : 0;
   }
  }
 }
 return dm;
}
template typename ATOM_DIAG_R::block_matrix_t atomic_density_matrix(ATOM_DIAG_R const&, double);
template typename ATOM_DIAG_C::block_matrix_t atomic_density_matrix(ATOM_DIAG_C const&, double);

// -----------------------------------------------------------------

// Given a monomial (ccccc), and a subspace B, returns
//   - the subspace connected by ccccc from B
//   - the corresponding matrix (not necessarily square)
template<bool Complex>
auto matrix_element_of_monomial(ATOM_DIAG const& atom, operators::monomial_t const& op_vec, int B)
 -> std::pair<int, ATOM_DIAG_T::matrix_t> {

 auto m = triqs::arrays::make_unit_matrix<ATOM_DIAG_T::scalar_t>(atom.get_subspace_dim(B));
 auto const& fops = atom.get_fops();
 for (int i = op_vec.size() - 1; i >= 0; --i) {
  int ind = fops[op_vec[i].indices];
  int Bp = (op_vec[i].dagger ? atom.cdag_connection(ind, B) : atom.c_connection(ind, B));
  if (Bp == -1) return {-1, std::move(m)};
  m = (op_vec[i].dagger ? atom.cdag_matrix(ind, B) : atom.c_matrix(ind, B)) * m;
  B = Bp;
 }
 return {B, std::move(m)};
}

// -----------------------------------------------------------------

template<bool Complex>
ATOM_DIAG_T::scalar_t trace_rho_op(ATOM_DIAG_T::block_matrix_t const& density_matrix,
                                   ATOM_DIAG_T::many_body_op_t const& op,
                                   ATOM_DIAG const& atom) {
 ATOM_DIAG_T::scalar_t result = 0;
 if (atom.n_subspaces() != density_matrix.size())
  TRIQS_RUNTIME_ERROR << "trace_rho_op : size mismatch : number of blocks differ";
 for (int sp = 0; sp < atom.n_subspaces(); ++sp) {
  if (atom.get_subspace_dim(sp) != first_dim(density_matrix[sp]))
   TRIQS_RUNTIME_ERROR << "trace_rho_op : size mismatch : size of block " << sp << " differ";
  for (auto const& x : op) {
   auto b_m = matrix_element_of_monomial(atom, x.monomial, sp);
   if (b_m.first == sp) result += x.coef * trace(b_m.second * density_matrix[sp]);
  }
 }
 return result;
}
template ATOM_DIAG_R::scalar_t trace_rho_op(ATOM_DIAG_R::block_matrix_t const&,
                                            ATOM_DIAG_R::many_body_op_t const&,
                                            ATOM_DIAG_R const&);
template ATOM_DIAG_C::scalar_t trace_rho_op(ATOM_DIAG_C::block_matrix_t const&,
                                            ATOM_DIAG_C::many_body_op_t const&,
                                            ATOM_DIAG_C const&);

// -----------------------------------------------------------------
template<bool Complex>
auto act(ATOM_DIAG_T::many_body_op_t const& op,
         ATOM_DIAG_T::full_hilbert_space_state_t const& st,
         ATOM_DIAG const& atom) -> ATOM_DIAG_T::full_hilbert_space_state_t {
 ATOM_DIAG_T::full_hilbert_space_state_t result(st.size());
 result() = 0;
 for (auto const& x : op) {
  for (int bl = 0; bl < atom.n_subspaces(); ++bl) {
   auto b_m = matrix_element_of_monomial(atom, x.monomial, bl);
   if (b_m.first == -1) continue;
   result(atom.index_range_of_subspace(b_m.first)) += x.coef * b_m.second * st(atom.index_range_of_subspace(bl));
  }
 }
 return result;
}
template ATOM_DIAG_R::full_hilbert_space_state_t act(ATOM_DIAG_R::many_body_op_t const&,
                                                     ATOM_DIAG_R::full_hilbert_space_state_t const&,
                                                     ATOM_DIAG_R const&);
template ATOM_DIAG_C::full_hilbert_space_state_t act(ATOM_DIAG_C::many_body_op_t const&,
                                                     ATOM_DIAG_C::full_hilbert_space_state_t const&,
                                                     ATOM_DIAG_C const&);

// -----------------------------------------------------------------
template<bool Complex>
auto quantum_number_eigenvalues(ATOM_DIAG_T::many_body_op_t const& op,
                                ATOM_DIAG const& atom) -> std::vector<std::vector<quantum_number_t>> {

 auto commutator = op * atom.get_h_atomic() - atom.get_h_atomic() * op;
 if (!commutator.is_zero())
  TRIQS_RUNTIME_ERROR << "The operator is not a quantum number";

 std::vector<std::vector<quantum_number_t>> result;

 for (int sp = 0; sp < atom.n_subspaces(); ++sp) {
  auto dim = atom.get_subspace_dim(sp);
  result.push_back(std::vector<quantum_number_t>(dim, 0));
  for (auto const& x : op) {
   auto b_m = matrix_element_of_monomial(atom, x.monomial, sp);
   if (b_m.first != sp) continue;
   for (int i = 0; i < dim; ++i) result.back()[i] += real(x.coef * b_m.second(i, i));
  }
 }
 return result;
}
template std::vector<std::vector<quantum_number_t>> quantum_number_eigenvalues(ATOM_DIAG_R::many_body_op_t const&,
                                                                               ATOM_DIAG_R const&);
template std::vector<std::vector<quantum_number_t>> quantum_number_eigenvalues(ATOM_DIAG_C::many_body_op_t const&,
                                                                               ATOM_DIAG_C const&);

// -----------------------------------------------------------------
template <typename M>
// require (ImmutableMatrix<M>)
inline bool is_diagonal(M const& m) {
 double static const threshold = 1.e-11;
 return triqs::utility::is_zero(sum(abs(m)) - trace(abs(m)),threshold);
}

template<bool Complex>
auto quantum_number_eigenvalues_checked(ATOM_DIAG_T::many_body_op_t const& op,
                                        ATOM_DIAG const& atom) -> std::vector<std::vector<quantum_number_t>> {

 auto commutator = op * atom.get_h_atomic() - atom.get_h_atomic() * op;
 if (!commutator.is_zero()) TRIQS_RUNTIME_ERROR << "The operator is not a quantum number";

 auto d = atom.get_full_hilbert_space_dim();
 matrix<quantum_number_t> M(d, d);
 M() = 0;
 std::vector<std::vector<quantum_number_t>> result;

 for (int sp = 0; sp < atom.n_subspaces(); ++sp) {
  auto dim = atom.get_subspace_dim(sp);
  for (auto const& x : op) {
   auto b_m = matrix_element_of_monomial(atom, x.monomial, sp);
   if (b_m.first == -1) continue;
   M(atom.index_range_of_subspace(b_m.first), atom.index_range_of_subspace(sp)) += real(x.coef * b_m.second);
  }
 }
 if (!is_diagonal(M)) TRIQS_RUNTIME_ERROR << "The matrix of the operator is not diagonal !!!";

 for (int sp = 0; sp < atom.n_subspaces(); ++sp) {
   auto dim = atom.get_subspace_dim(sp);
   result.push_back(std::vector<quantum_number_t>(dim, 0));
   for (int i = 0; i < dim; ++i) result.back()[i] = M(atom.flatten_subspace_index(sp, i),
                                                      atom.flatten_subspace_index(sp, i));
 }

 return result;
}
template std::vector<std::vector<quantum_number_t>> quantum_number_eigenvalues_checked(ATOM_DIAG_R::many_body_op_t const&,
                                                                                       ATOM_DIAG_R const&);
template std::vector<std::vector<quantum_number_t>> quantum_number_eigenvalues_checked(ATOM_DIAG_C::many_body_op_t const&,
                                                                                       ATOM_DIAG_C const&);
}}
