/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2014, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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
#pragma once
#include <triqs/gfs.hpp>
#include "triqs/statistics/histograms.hpp" //FIXME
#include "./qmc_data.hpp"
#include "./impurity_trace.hpp"

namespace cthyb {

using namespace triqs::gfs;

// Measure the four body correlator C_abcd(tau - tau') = < (c+_a c_b) (tau) (c+_c c_d) (tau') >
struct measure_four_body_corr {
 using mc_sign_type = double;
 using node = impurity_trace::node;
 using trace_t = impurity_trace::trace_t;
 using block_and_matrix = impurity_trace::block_and_matrix;

 qmc_data const& data;
 gf_view<imtime, scalar_valued> correlator;
 mc_sign_type z;
 qmc_data::trace_t new_atomic_weight, new_atomic_reweighting;
 arrays::array<dcomplex, 4> coefficients;          // Coefficients of op*op, where op is a quadratic operator
 arrays::array<dcomplex, 2> coefficients_one_pair; // Max abs coefficient of op
 bool anticommute;                                 // Do the cdag and c operators anticommute?
 impurity_trace::rb_tree_t& tree = data.imp_trace.tree;

//FIXME
 statistics::histogram_segment_bin binned_taus = {0, data.config.beta(), 100, "histo_binned_taus.dat"};

 measure_four_body_corr(qmc_data const& data, gf_view<imtime, scalar_valued> correlator, fundamental_operator_set const & fops, many_body_operator const & A, bool anticommute);
 void accumulate(mc_sign_type s);
 void collect_results(triqs::mpi::communicator const& c);

 impurity_trace::block_and_matrix integral(int b_i, time_pt tau_i, time_pt tau_f, op_desc const& op1, op_desc const& op2);
 impurity_trace::block_and_matrix integral(int b_i, time_pt tau_i, time_pt tau_f, op_desc const& op1, op_desc const& op2,
                                           op_desc const& op3, op_desc const& op4);
 double compute_evolution_integral(double lamb1, double lamb2, double lamb3, double lamb4, double lamb5);
 double compute_evolution_integral(double lamb1, double lamb2, double lamb3);
 std::pair<trace_t, trace_t> compute_sliding_trace_integral(std::vector<node> const& flat_config, int index_node_l,
                                                            int index_node_r, std::vector<int> const& blocks);
 std::pair<trace_t, trace_t> compute_sliding_trace_integral_one_pair(std::vector<node> const& flat_config, int index_node,
                                                                     std::vector<int> const& blocks);
};
}