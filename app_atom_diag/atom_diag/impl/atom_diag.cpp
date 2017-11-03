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
#include "../atom_diag.hpp"
#include "./worker.hpp"

#include <triqs/arrays.hpp>

using namespace triqs::arrays;

namespace triqs {
namespace atom_diag {

// Methods of atom_diag
#define ATOM_DIAG_CONSTRUCTOR(ARGS) template<bool Complex> atom_diag<Complex>::atom_diag ARGS
#define ATOM_DIAG_METHOD(RET,F) template<bool Complex> auto atom_diag<Complex>::F -> RET

ATOM_DIAG_CONSTRUCTOR((many_body_op_t const& h, fundamental_operator_set const& fops, std::vector<many_body_op_t> const& qn_vector)) :
 h_atomic(h), fops(fops), full_hs(fops), vacuum(full_hs.size()) {
 atom_diag_worker<Complex>{this}.partition_with_qn(qn_vector);
 fill_first_eigenstate_of_subspace();
 compute_vacuum();
}

// -----------------------------------------------------------------

ATOM_DIAG_CONSTRUCTOR((many_body_op_t const& h, fundamental_operator_set const& fops)) :
 h_atomic(h), fops(fops), full_hs(fops), vacuum(full_hs.size()) {
 atom_diag_worker<Complex>{this}.autopartition();
 fill_first_eigenstate_of_subspace();
 compute_vacuum();
}

// -----------------------------------------------------------------

ATOM_DIAG_METHOD(void,fill_first_eigenstate_of_subspace()) {
 // Calculate the index of the first eigenstate of each block
 first_eigenstate_of_subspace.resize(n_subspaces());
 first_eigenstate_of_subspace[0] = 0;
 for (int sp = 1; sp < n_subspaces(); ++sp)
  first_eigenstate_of_subspace[sp] = first_eigenstate_of_subspace[sp - 1] + get_subspace_dim(sp - 1);
}

// -----------------------------------------------------------------

ATOM_DIAG_METHOD(void,compute_vacuum()) {
 // Compute vacuum vector in the eigenbasis
 vacuum() = 0;
 for(int sp : range(sub_hilbert_spaces.size())) {
  if(sub_hilbert_spaces[sp].has_state(fock_state_t(0))) {
   vacuum_subspace_index = sp;
   vacuum(index_range_of_subspace(sp)) = dagger(eigensystems[sp].unitary_matrix)(range(),0);
   break;
  }
 }
}

// -----------------------------------------------------------------

ATOM_DIAG_METHOD(std::vector<std::vector<double>>,get_energies() const) {
 std::vector<std::vector<double>> R;
 for (auto const& es : eigensystems)
  R.emplace_back(es.eigenvalues.begin(), es.eigenvalues.end());
 return R;
}

#undef ATOM_DIAG_METHOD

// -----------------------------------------------------------------

// Free functions

template<bool Complex>
void h5_write(triqs::h5::group fg, std::string const& name, atom_diag<Complex> const& ad) {
 using matrix_t = typename atom_diag<Complex>::matrix_t;
 auto gr = fg.create_group(name);
 gr.write_triqs_hdf5_data_scheme(ad);
 h5_write(gr, "creation_connection", ad.creation_connection);
 h5_write(gr, "annihilation_connection", ad.annihilation_connection);

 auto write_sparse = [&](std::string na, std::vector<std::vector<matrix_t>> const& Mvv) {
 auto gr2 = gr.create_group(na);
 for (int i = 0; i < Mvv.size(); ++i)
  for (int j = 0; j < Mvv[i].size(); ++j)
   if (!Mvv[i][j].is_empty()) h5_write(gr2, std::to_string(i) + ' ' + std::to_string(j), Mvv[i][j]);
  };

 write_sparse("c_matrices", ad.c_matrices);
 write_sparse("cdag_matrices", ad.cdag_matrices);

 h5_write(gr, "full_hs", ad.full_hs);
 h5_write(gr, "sub_hilbert_spaces", ad.sub_hilbert_spaces);
 h5_write(gr, "eigensystems", ad.eigensystems);
 h5_write(gr, "gs_energy", ad.gs_energy);
 h5_write(gr, "vacuum_subspace_index", ad.vacuum_subspace_index);
 h5_write(gr, "vacuum", ad.vacuum);
 h5_write_attribute(gr, "fops", ad.fops);
}

// -----------------------------------------------------------------
template<bool Complex>
void h5_read(triqs::h5::group fg, std::string const& name, atom_diag<Complex> & ad) {
 using matrix_t = typename atom_diag<Complex>::matrix_t;
 auto gr = fg.open_group(name);

 h5_read(gr, "creation_connection", ad.creation_connection);
 h5_read(gr, "annihilation_connection", ad.annihilation_connection);

 auto read_sparse = [&](std::string na, std::vector<std::vector<matrix_t>>& Mvv) {
  Mvv.resize(first_dim(ad.creation_connection), std::vector<matrix_t>(second_dim(ad.creation_connection)));
  auto gr2 = gr.open_group(na);
  for (auto s : gr2.get_all_dataset_names()) {
   std::stringstream ss(s);
   std::string item1, item2;
   std::getline(ss, item1, ' ');
   std::getline(ss, item2, ' ');
   int i = std::stoi(item1), j = std::stoi(item2);
   h5_read(gr2, s, Mvv[i][j]);
  }
 };

 read_sparse("c_matrices", ad.c_matrices);
 read_sparse("cdag_matrices", ad.cdag_matrices);

 h5_read(gr, "full_hs", ad.full_hs);
 h5_read(gr, "sub_hilbert_spaces", ad.sub_hilbert_spaces);
 h5_read(gr, "eigensystems", ad.eigensystems);
 h5_read(gr, "gs_energy", ad.gs_energy);
 h5_read(gr, "vacuum_subspace_index", ad.vacuum_subspace_index);
 h5_read(gr, "vacuum", ad.vacuum);
 h5_read_attribute(gr, "fops", ad.fops);
 ad.fill_first_eigenstate_of_subspace();
}

// -----------------------------------------------------------------

template<bool Complex>
std::ostream& operator<<(std::ostream& os, atom_diag<Complex> const& ad) {
 os << "Dimension of full Hilbert space: " << ad.get_full_hilbert_space_dim() << std::endl;
 os << "Number of invariant subspaces: " << ad.n_subspaces() << std::endl;
 for (int n_sp = 0; n_sp < ad.n_subspaces(); ++n_sp) {
  os << "Subspace " << n_sp << ", ";
  os << "lowest energy level : " << ad.eigensystems[n_sp].eigenvalues[0] << std::endl;
  os << "Subspace dimension = " << ad.eigensystems[n_sp].eigenvalues.size() << std::endl;
  os << "-------------------------" << std::endl;
 }
 return os;
}

// -----------------------------------------------------------------

// Explicit instantiations for real and complex atom_diag
template class atom_diag<false>;
template class atom_diag<true>;

template void h5_write(triqs::h5::group fg, std::string const&, atom_diag<false> const&);
template void h5_write(triqs::h5::group fg, std::string const&, atom_diag<true> const&);
template void h5_read(triqs::h5::group fg, std::string const&, atom_diag<false> &);
template void h5_read(triqs::h5::group fg, std::string const&, atom_diag<true> &);

template std::ostream& operator<<(std::ostream&, atom_diag<true> const&);
template std::ostream& operator<<(std::ostream&, atom_diag<false> const&);

}}
