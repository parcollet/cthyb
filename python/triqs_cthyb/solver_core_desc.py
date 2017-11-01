# Generated automatically using the command :
# c++2py.py ../c++/solver_core.hpp -I../cbuild/c++ -I../c++ -p -m pytriqs.applications.impurity_solvers.cthyb -o cthyb --moduledoc "The cthyb solver"
from cpp2py.wrap_generator import *

# The module
#module = module_(full_name = "pytriqs.applications.impurity_solvers.cthyb", doc = "The cthyb solver", app_name = "pytriqs.applications.impurity_solvers.cthyb")
module = module_(full_name = "solver_core", doc = "", app_name = "solver_core")

# Imports
import pytriqs.gf
import pytriqs.operators
import pytriqs.statistics

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("cthyb/solver_core.hpp")
module.add_include("cthyb/atom_diag.hpp")

module.add_enum("block_order", ["block_order::AABB","block_order::ABBA"], "cthyb", "Order of block indices for Block2Gf objects")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/map.hpp>
#include <cpp2py/converters/set.hpp>
#include <cpp2py/converters/optional.hpp>
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/vector.hpp>
#include <cpp2py/converters/tuple.hpp>

#include <triqs/cpp2py_converters/arrays.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
#include <triqs/cpp2py_converters/operators_real_complex.hpp>
#include <triqs/cpp2py_converters/fundamental_operator_set.hpp>
#include <triqs/cpp2py_converters/variant.hpp>

using namespace triqs::gfs;
using triqs::operators::many_body_operator;
using namespace cthyb;
#include "solver_core_converters.hxx"
""")

# The class solver_core
c = class_(
        py_type = "SolverCore",  # name of the python class
        c_type = "solver_core",   # name of the C++ class
        doc = r"Core class of the cthyb solver",   # doc of the C++ class
)

c.add_constructor("""(double beta, std::map<std::string,indices_type> gf_struct, int n_iw = 1025, int n_tau = 10001, int n_l = 50)""",
                  doc = """ """)

c.add_member(c_name = "G_tau",
             c_type = "std::optional<G_tau_t>",
             read_only= False,
             doc = """Single-particle Green\'s function :math:`G(\\tau)` in imaginary time. """)

c.add_member(c_name = "G_tau_accum",
             c_type = "std::optional<G_tau_G_target_t>",
             read_only= False,
             doc = """Intermediate Green\'s function to accumulate g(tau), either real or complex """)

c.add_member(c_name = "G_l",
             c_type = "std::optional<G_l_t>",
             read_only= False,
             doc = """Single-particle Green\'s function :math:`G_l` in Legendre polynomial representation. """)

c.add_member(c_name = "G2_tau",
             c_type = "std::optional<G2_tau_t>",
             read_only= False,
             doc = """Two-particle Green\'s function :math:`G^{(2)}(\\tau_1,\\tau_2,\\tau_3)` (three Fermionic imaginary times) """)

c.add_member(c_name = "G2_iw",
             c_type = "std::optional<G2_iw_t>",
             read_only= False,
             doc = """Two-particle Green\'s function :math:`G^{(2)}(i\\nu,i\\nu\',i\\nu\'\')` (three Fermionic frequencies) """)

c.add_member(c_name = "G2_iw_pp",
             c_type = "std::optional<G2_iw_t>",
             read_only= False,
             doc = """Two-particle Green\'s function :math:`G^{(2)}(i\\omega,i\\nu,i\\nu\')` in the pp-channel (one bosonic matsubara and two fermionic) """)

c.add_member(c_name = "G2_iw_ph",
             c_type = "std::optional<G2_iw_t>",
             read_only= False,
             doc = """Two-particle Green\'s function :math:`G^{(2)}(i\\omega,i\\nu,i\\nu\')` in the ph-channel (one bosonic matsubara and two fermionic) """)

c.add_member(c_name = "G2_iwll_pp",
             c_type = "std::optional<G2_iwll_t>",
             read_only= False,
             doc = """Two-particle Green\'s function :math:`G^{(2)}(i\\omega,l,l\')` in the pp-channel (one bosonic matsubara and two legendre) """)

c.add_member(c_name = "G2_iwll_ph",
             c_type = "std::optional<G2_iwll_t>",
             read_only= False,
             doc = """Two-particle Green\'s function :math:`G^{(2)}(i\\omega,l,l\')` in the ph-channel (one bosonic matsubara and two legendre) """)

c.add_method("""void solve (**cthyb::solve_parameters_t)""",
             doc = """+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| Parameter Name                | Type                                           | Default                       | Documentation                                                                                           |
+===============================+================================================+===============================+=========================================================================================================+
| h_int                         | Operator                                       |                               | Interacting part of the atomic Hamiltonian                                                              |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| n_cycles                      | int                                            |                               | Number of QMC cycles                                                                                    |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| partition_method              | str                                            | "autopartition"               | Partition method                                                                                        |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| quantum_numbers               | list(Operator)                                 | []                            | Quantum numbers                                                                                         |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| length_cycle                  | int                                            | 50                            | Length of a single QMC cycle                                                                            |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| n_warmup_cycles               | int                                            | 5000                          | Number of cycles for thermalization                                                                     |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| random_seed                   | int                                            | 34788 + 928374 * MPI.rank     | Seed for random number generator                                                                        |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| random_name                   | str                                            | ""                            | Name of random number generator                                                                         |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| max_time                      | int                                            | -1 = infinite                 | Maximum runtime in seconds, use -1 to set infinite                                                      |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| verbosity                     | int                                            | 3 on MPI rank 0, 0 otherwise. | Verbosity level                                                                                         |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| move_shift                    | bool                                           | true                          | Add shifting an operator as a move?                                                                     |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| move_double                   | bool                                           | false                         | Add double insertions as a move?                                                                        |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| use_trace_estimator           | bool                                           | false                         | Calculate the full trace or use an estimate?                                                            |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_G_tau                 | bool                                           | true                          | Measure G(tau)?                                                                                         |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_G_l                   | bool                                           | false                         | Measure G_l (Legendre)?                                                                                 |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_G2_tau                | bool                                           | false                         | Measure G^4(tau,tau',tau'') with three fermionic times.                                                 |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_G2_iw                 | bool                                           | false                         | Measure G^4(inu,inu',inu'') with three fermionic frequencies.                                           |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_G2_iw_pp              | bool                                           | false                         | Measure G^4(iomega,inu,inu') within the particle-particle channel.                                      |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_G2_iw_ph              | bool                                           | false                         | Measure G^4(iomega,inu,inu') within the particle-hole channel.                                          |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_G2_iwll_pp            | bool                                           | false                         | Measure G^2(iomega,l,l') within the particle-particle channel.                                          |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_G2_iwll_ph            | bool                                           | false                         | Measure G^2(iomega,l,l') within the particle-hole channel.                                              |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_G2_block_order        | cthyb::block_order                             | block_order::AABB             | Order of block indices in the definition of G^2.                                                        |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_G2_blocks             | std::set<std::pair<std::string, std::string> > | measure all blocks            | List of block index pairs of G^2 to measure.                                                            |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_G2_n_tau              | int                                            | 10                            | Number of imaginary time slices for G^4 measurement.                                                    |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_G2_n_bosonic          | int                                            | 30                            | Number of bosonic Matsubara frequencies for G^4 measurement.                                            |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_G2_n_fermionic        | int                                            | 30                            | Number of fermionic Matsubara frequencies for G^4 measurement.                                          |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_G2_n_l                | int                                            | 20                            | Number of Legendre coefficients for G^4(iomega,l,l') measurement.                                       |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_G2_iwll_nfft_buf_size | int                                            | 100                           | NFFT buffer size for G^4(iomega,l,l') measurement.                                                      |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| nfft_buf_sizes                | std::map<std::string, int>                     | 100 for every block           | NFFT buffer sizes for different blocks                                                                  |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_pert_order            | bool                                           | false                         | Measure perturbation order?                                                                             |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_density_matrix        | bool                                           | false                         | Measure the reduced impurity density matrix?                                                            |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| use_norm_as_weight            | bool                                           | false                         | Use the norm of the density matrix in the weight if true, otherwise use Trace                           |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| performance_analysis          | bool                                           | false                         | Analyse performance of trace computation with histograms (developers only)?                             |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| proposal_prob                 | dict(str:float)                                | {}                            | Operator insertion/removal probabilities for different blocks                                           |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| move_global                   | dict(str : dict(indices : indices))            | {}                            | List of global moves (with their names). Each move is specified with an index substitution dictionary.  |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| move_global_prob              | double                                         | 0.05                          | Overall probability of the global moves                                                                 |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| imag_threshold                | double                                         | 1.e-15                        | Threshold below which imaginary components of Delta and h_loc are set to zero                           |
+-------------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+ """)

c.add_property(name = "h_loc",
               getter = cfunction("many_body_op_t h_loc ()"),
               doc = """The local Hamiltonian of the problem: :math:`H_{loc}` used in the last call to ``solve()``. """)

c.add_property(name = "last_solve_parameters",
               getter = cfunction("cthyb::solve_parameters_t last_solve_parameters ()"),
               doc = """Set of parameters used in the last call to ``solve()``. """)

c.add_property(name = "Delta_tau",
               getter = cfunction("block_gf_view<triqs::gfs::imtime> Delta_tau ()"),
               doc = """:math:`\\Delta(\\tau)` in imaginary time. """)

c.add_property(name = "G0_iw",
               getter = cfunction("block_gf_view<triqs::gfs::imfreq> G0_iw ()"),
               doc = """:math:`G_0(i\\omega)` in imaginary frequencies. """)

c.add_property(name = "atomic_gf",
               getter = cfunction("block_gf_view<triqs::gfs::imtime> atomic_gf ()"),
               doc = """Atomic :math:`G(\\tau)` in imaginary time. """)

c.add_property(name = "density_matrix",
               getter = cfunction("std::vector<matrix_t> density_matrix ()"),
               doc = """Accumulated density matrix. """)

c.add_property(name = "h_loc_diagonalization",
               getter = cfunction("cthyb::atom_diag h_loc_diagonalization ()"),
               doc = """Diagonalization of :math:`H_{loc}`. """)

c.add_property(name = "perturbation_order_total",
               getter = cfunction("triqs::statistics::histogram get_perturbation_order_total ()"),
               doc = """Histogram of the total perturbation order. """)

c.add_property(name = "perturbation_order",
               getter = cfunction("histo_map_t get_perturbation_order ()"),
               doc = """Histograms of the perturbation order for each block. """)

c.add_property(name = "performance_analysis",
               getter = cfunction("histo_map_t get_performance_analysis ()"),
               doc = """Histograms related to the performance analysis. """)

c.add_property(name = "average_sign",
               getter = cfunction("mc_weight_t average_sign ()"),
               doc = """Monte Carlo average sign. """)

c.add_property(name = "solve_status",
               getter = cfunction("int solve_status ()"),
               doc = """Status of the ``solve()`` on exit. """)

module.add_class(c)


# The class atom_diag
c = class_(
        py_type = "AtomDiag",  # name of the python class
        c_type = "atom_diag",   # name of the C++ class
        hdf5 = True,
        doc = r"",   # doc of the C++ class
)

c.add_constructor("""(many_body_op_t h_, triqs::hilbert_space::fundamental_operator_set fops)""",
                  doc = """ """)

c.add_constructor("""(many_body_op_t h_, triqs::hilbert_space::fundamental_operator_set fops, std::vector<many_body_op_t> qn_vector)""",
                  doc = """ """)

c.add_method("""int get_block_dim (int b)""",
             doc = """The dimension of block b """)

c.add_method("""int flatten_block_index (int block_index, int i)""",
             doc = """Returns the index in the full hilbert space for block_index and i, the index within the block. """)

c.add_method("""double get_eigenvalue (int block_index, int i)""",
             doc = """Get the i-th eigenvalue of block bl """)

c.add_method("""long c_connection (int op_linear_index, int block_index)""",
             doc = """Connections for fundamental operators C\n\n op_linear_index : the linear index (i.e. number) of the c operator, as defined by the fundamental_operator_set fops\n block_number : the number of the initial block\n @return : the number of the final block """)

c.add_method("""long cdag_connection (int op_linear_index, int block_index)""",
             doc = """Connections for fundamental operators C^\\dagger\n\n op_linear_index : the linear index (i.e. number) of the c operator, as defined by the fundamental_operator_set fops\n block_number : the number of the initial block\n @return : the number of the final block """)

c.add_method("""matrix<h_scalar_t> c_matrix (int op_linear_index, int block_index)""",
             doc = """Matrix for fundamental operators C\n\n op_linear_index : the linear index (i.e. number) of the c operator, as defined by the fundamental_operator_set fops\n block_number : the number of the initial block\n @return : the number of the final block """)

c.add_method("""matrix<h_scalar_t> cdag_matrix (int op_linear_index, int block_index)""",
             doc = """Matrix for fundamental operators C^\\dagger\n\n op_linear_index : the linear index (i.e. number) of the c operator, as defined by the fundamental_operator_set fops\n block_number : the number of the initial block\n @return : the number of the final block """)

c.add_property(name = "h_atomic",
               getter = cfunction("many_body_op_t get_h_atomic ()"),
               doc = """The Hamiltonian """)

c.add_property(name = "fops",
               getter = cfunction("triqs::hilbert_space::fundamental_operator_set get_fops ()"),
               doc = """The fundamental operator set used at construction """)

c.add_property(name = "full_hilbert_space_dim",
               getter = cfunction("int get_full_hilbert_space_dim ()"),
               doc = """Dimension of the full Hilbert space """)

c.add_property(name = "n_blocks",
               getter = cfunction("int n_blocks ()"),
               doc = """Number of Blocks """)

c.add_property(name = "fock_states",
               getter = cfunction("std::vector<std::vector<fock_state_t>> get_fock_states ()"),
               doc = """The list of fock states for each block """)

c.add_property(name = "unitary_matrices",
               getter = cfunction("std::vector<matrix<h_scalar_t>> get_unitary_matrices ()"),
               doc = """Unitary matrices that transform from Fock states to atomic eigenstates """)

c.add_property(name = "energies",
               getter = cfunction("std::vector<std::vector<double>> get_energies ()"),
               doc = """A vector of all the energies, by blocks. result[block_number][i] is the energy """)

c.add_property(name = "quantum_numbers",
               getter = cfunction("std::vector<std::vector<double>> get_quantum_numbers ()"),
               doc = """A vector of all the QNs, by blocks : result[block_number][qn_index] is the ..... """)

c.add_property(name = "gs_energy",
               getter = cfunction("double get_gs_energy ()"),
               doc = """Ground state energy (i.e. min of all subspaces). """)

c.add_property(name = "vacuum_block_index",
               getter = cfunction("int get_vacuum_block_index ()"),
               doc = """Returns the block index of the vacuum state. """)

c.add_property(name = "vacuum_inner_index",
               getter = cfunction("int get_vacuum_inner_index ()"),
               doc = """Returns the inner index of the vacuum state. """)

c.add_property(name = "vacuum_state",
               getter = cfunction("full_hilbert_space_state_t get_vacuum_state ()"),
               doc = """Returns the vacuum state as a long vector in the full Hilbert space. """)

module.add_class(c)

#  Free functions

module.add_function ("double partition_function (cthyb::atom_diag atom, double beta)", doc = """The atomic partition function""")

module.add_function ("block_matrix_t atomic_density_matrix (cthyb::atom_diag atom, double beta)", doc = """The atomic density matrix""")

module.add_function ("block_gf<imtime> atomic_gf (cthyb::atom_diag atom, double beta, std::map<std::string,indices_t> indices_list, int n_tau, std::vector<std::pair<int,int>> excluded_states = {})", doc = """The atomic green function, possibly with excluded states (default none)""")

module.add_function ("double trace_rho_op (block_matrix_t density_matrix, many_body_op_t op, cthyb::atom_diag atom)", doc = "Trace (op * density_matrix)")

module.add_function ("full_hilbert_space_state_t act (many_body_op_t op, full_hilbert_space_state_t st, cthyb::atom_diag atom)", doc = """Act with operator op on state st""")

module.add_function ("std::vector<std::vector<double>> quantum_number_eigenvalues (many_body_op_t op, cthyb::atom_diag atom)", doc = """The operator op is supposed to be a quantum number (if not -> exception)\n @return the eigenvalue by block""")

module.add_function ("std::vector<std::vector<double>> quantum_number_eigenvalues2 (many_body_op_t op, cthyb::atom_diag atom)", doc = """The operator op is supposed to be a quantum number (if not -> exception)\n @return the eigenvalue by block""")

module.generate_code()
