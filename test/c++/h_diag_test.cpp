#include <triqs/test_tools/arrays.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>

#include "./atom_diag.hpp"

using namespace cthyb;
using triqs::operators::c;
using triqs::operators::c_dag;

// define operators
auto n_up = c_dag("up") * c("up");
auto n_down = c_dag("down") * c("down");

// Hubbard atom in magnetic field
TEST(h_diag, Hubbard) {
 // choose the fundamental operator set
 fundamental_operator_set fops;
 fops.insert("up");
 fops.insert("down");

 double U = 2.0;
 double h = 0.5;
 double d = 0.5;

 auto H = U * n_up * n_down + h*(n_up - n_down);

 // put quantum numbers in a vector
 std::vector<many_body_op_t> qn_list{n_up, n_down};

 // Divide the full Hilbert space
 atom_diag h_diag(H, fops, qn_list);
 EXPECT_EQ(4, h_diag.n_blocks());

 double block_size_ref[] = {1, 1, 1, 1};
 double gs_energy_ref[] = {0, 0.5, 1, 2.5};
 for(int i = 0; i < h_diag.n_blocks(); ++i) {
  EXPECT_EQ(block_size_ref[i], h_diag.get_eigensystem()[i].eigenvalues.size());
  EXPECT_CLOSE(gs_energy_ref[i], h_diag.get_eigensystem()[i].eigenvalues[0]);
 }
}

// Half-filled Hubbard atom + spin flips + anomalous terms
TEST(h_diag, Hubbard_anomalous) {
 // choose the fundamental operator set
 fundamental_operator_set fops;
 fops.insert("up");
 fops.insert("down");

 double U = 2.0;
 double J = 0.3;
 double d = 0.1;

 auto H = U * (n_up - 0.5) * (n_down - 0.5);
 H += J*(c_dag("up") * c("down") + c_dag("down") * c("up"));
 H += d * (c_dag("up") * c_dag("down") - c("up") * c("down"));

 // Divide the full Hilbert space
 atom_diag h_diag(H, fops);
 EXPECT_EQ(2, h_diag.n_blocks());

 double block_size_ref[] = {2, 2};
 double gs_energy_ref[] = {0, 1.2};
 for(int i = 0; i < h_diag.n_blocks(); ++i) {
  EXPECT_EQ(block_size_ref[i], h_diag.get_eigensystem()[i].eigenvalues.size());
  EXPECT_CLOSE(gs_energy_ref[i], h_diag.get_eigensystem()[i].eigenvalues[0]);
 }
}

MAKE_MAIN;
