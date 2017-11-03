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
#pragma once

#include <vector>
#include <triqs/gfs.hpp>
#include "./atom_diag.hpp"

namespace triqs {
namespace atom_diag {

 using namespace triqs::gfs;

#define ATOM_DIAG atom_diag<Complex>
#define ATOM_DIAG_T typename atom_diag<Complex>

/// Lehmann representation of one matrix element of a GF: list of (pole, residue) pairs
template<bool Complex>
using gf_scalar_lehmann_t = std::vector<std::pair<double, ATOM_DIAG_T::scalar_t>>;

/// Lehmann representation of a matrix-valued GF
template<bool Complex>
using gf_lehmann_t = std::vector<matrix<gf_scalar_lehmann_t<Complex>>>;

/// List of excluded eigenstates, (subspace index, inner index) pairs
using excluded_states_t = std::vector<std::pair<int, int>>;

/// The atomic Green's function, Lehmann representation
/**
 * @tparam Complex Do we have a diagonalization problem with a complex-valued Hamiltonian?
 * @param atom Solved diagonalization problem.
 * @param beta Inverse temperature.
 * @param gf_struct Block structure of the Green's function, block name -> list of inner indices.
 * @param excluded_states Excluded eigenstates as pairs (subspace index, inner index).
 * @return Atomic Green's function in the Lehmann representation
 * @include triqs/atom_diag/gf.hpp
 */
template<bool Complex>
gf_lehmann_t<Complex> atomic_g_lehmann(ATOM_DIAG const& atom,
                                       double beta,
                                       gf_struct_t const& gf_struct,
                                       excluded_states_t excluded_states = {});

/// The atomic imaginary time Green's function, constructed from precomputed Lehmann representation
/**
 * @tparam Complex Do we have a diagonalization problem with a complex-valued Hamiltonian?
 * @param lehmann Lehmann representation.
 * @param gf_struct Block structure of the Green's function, block name -> list of inner indices.
 * @param mesh Imaginary time mesh used in construction.
 * @return Atomic Green's function :math:`G_{at}(\tau)`
 * @include triqs/atom_diag/gf.hpp
 */
template<bool Complex>
block_gf<imtime> atomic_g_tau(gf_lehmann_t<Complex> const& lehmann,
                              gf_struct_t const& gf_struct,
                              gf_mesh<imtime> const& mesh);

/// The atomic imaginary time Green's function, possibly with excluded states (none by default)
/**
 * @tparam Complex Are we using Lehmann representation with complex matrix elements?
 * @param atom Solved diagonalization problem.
 * @param beta Inverse temperature.
 * @param gf_struct Block structure of the Green's function, block name -> list of inner indices.
 * @param n_tau Number of imaginary time points.
 * @param excluded_states Excluded eigenstates as pairs (subspace index, inner index).
 * @return Atomic Green's function :math:`G_{at}(\tau)`
 * @include triqs/atom_diag/gf.hpp
 */
template<bool Complex>
block_gf<imtime> atomic_g_tau(ATOM_DIAG const& atom,
                              double beta,
                              gf_struct_t const& gf_struct,
                              int n_tau,
                              excluded_states_t const& excluded_states = {});

/// The atomic Matsubara Green's function, constructed from precomputed Lehmann representation
/**
 * @tparam Complex Are we using Lehmann representation with complex matrix elements?
 * @param lehmann Lehmann representation.
 * @param gf_struct Block structure of the Green's function, block name -> list of inner indices.
 * @param mesh Matsubara mesh used in construction.
 * @return Atomic Green's function :math:`G_{at}(i\omega)`.
 * @include triqs/atom_diag/gf.hpp
 */
template<bool Complex>
block_gf<imfreq> atomic_g_iw(gf_lehmann_t<Complex> const& lehmann,
                             gf_struct_t const& gf_struct,
                             gf_mesh<imfreq> const& mesh);

/// The atomic Matsubara Green's function, possibly with excluded states (none by default)
/**
 * @tparam Complex Do we have a diagonalization problem with a complex-valued Hamiltonian?
 * @param atom Solved diagonalization problem.
 * @param beta Inverse temperature.
 * @param gf_struct Block structure of the Green's function, block name -> list of inner indices.
 * @param n_iw Number of Matsubara frequencies.
 * @param excluded_states Excluded eigenstates as pairs (subspace index, inner index).
 * @return Atomic Green's function :math:`G_{at}(i\omega)`.
 * @include triqs/atom_diag/gf.hpp
 */
template<bool Complex>
block_gf<imfreq> atomic_g_iw(ATOM_DIAG const& atom,
                             double beta,
                             gf_struct_t const& gf_struct,
                             int n_iw,
                             excluded_states_t const& excluded_states = {});

/// The atomic Green's function in Legendre basis, constructed from precomputed Lehmann representation
/**
 * @tparam Complex Are we using Lehmann representation with complex matrix elements?
 * @param lehmann Lehmann representation.
 * @param gf_struct Block structure of the Green's function, block name -> list of inner indices.
 * @param mesh Legendre mesh used in construction.
 * @return Atomic Green's function :math:`G_{at}(\ell)`.
 * @include triqs/atom_diag/gf.hpp
 */
template<bool Complex>
block_gf<legendre> atomic_g_l(gf_lehmann_t<Complex> const& lehmann,
                              gf_struct_t const& gf_struct,
                              gf_mesh<legendre> const& mesh);

/// The atomic Green's function in Legendre basis, possibly with excluded states (none by default)
/**
 * @tparam Complex Do we have a diagonalization problem with a complex-valued Hamiltonian?
 * @param atom Solved diagonalization problem.
 * @param beta Inverse temperature.
 * @param gf_struct Block structure of the Green's function, block name -> list of inner indices.
 * @param n_l Number of Legendre coefficients.
 * @param excluded_states Excluded eigenstates as pairs (subspace index, inner index).
 * @return Atomic Green's function :math:`G_{at}(\ell)`.
 * @include triqs/atom_diag/gf.hpp
 */
template<bool Complex>
block_gf<legendre> atomic_g_l(ATOM_DIAG const& atom,
                              double beta,
                              gf_struct_t const& gf_struct,
                              int n_l,
                              excluded_states_t const& excluded_states = {});

/// The atomic retarded Green's function, constructed from precomputed Lehmann representation
/**
 * @tparam Complex Do we have a diagonalization problem with a complex-valued Hamiltonian?
 * @param lehmann Lehmann representation.
 * @param gf_struct Block structure of the Green's function, block name -> list of inner indices.
 * @param mesh Real frequency mesh used in construction.
 * @param broadening Lorentian broadening of the spectrum (imaginary frequency shift).
 * @return Atomic Green's function :math:`G_{at}(\omega)`.
 * @include triqs/atom_diag/gf.hpp
 */
template<bool Complex>
block_gf<refreq> atomic_g_w(gf_lehmann_t<Complex> const& lehmann,
                              gf_struct_t const& gf_struct,
                              gf_mesh<refreq> const& mesh,
                              double broadening = 0);

/// The atomic retarded Green's function, possibly with excluded states (none by default)
/**
 * @tparam Complex Do we have a diagonalization problem with a complex-valued Hamiltonian?
 * @param atom Solved diagonalization problem.
 * @param beta Inverse temperature.
 * @param gf_struct Block structure of the Green's function, block name -> list of inner indices.
 * @param energy_window Energy window :math:`(\omega_{min}, \omega_{max})`.
 * @param n_w Number of frequency points.
 * @param broadening Lorentian broadening of the spectrum (imaginary frequency shift).
 * @param excluded_states Excluded eigenstates as pairs (subspace index, inner index).
 * @return Atomic Green's function :math:`G_{at}(\omega)`.
 * @include triqs/atom_diag/gf.hpp
 */
template<bool Complex>
block_gf<refreq> atomic_g_w(ATOM_DIAG const& atom,
                            double beta, gf_struct_t const& gf_struct,
                            std::pair<double, double> const& energy_window,
                            int n_w,
                            double broadening = 0,
                            excluded_states_t const& excluded_states = {});

#undef ATOM_DIAG
#undef ATOM_DIAG_T

}}
