from pytriqs.archive import HDFArchive
import pytriqs.utility.mpi as mpi
from pytriqs.gf.local import *
from pytriqs.operators import *
from pytriqs.applications.impurity_solvers.cthyb import *

#  Example of DMFT single site solution with CTQMC

# set up a few parameters
half_bandwidth = 0.05
U = 10.0
mu = U/2.0
beta = 1.00

use_prop = True

gf_struct = {'up':[0],'down':[0]}

# Construct solver    
S = Solver(beta=beta, gf_struct=gf_struct, n_iw=1025, n_tau=3000, n_l=30)

# Local Hamiltonian
H = U*n("up",0)*n("down",0)

# init the Green function
S.G_iw << SemiCircular(half_bandwidth)

# Impose Paramagnetism
gpara = 0.5*(S.G_iw['up']+S.G_iw['down'])
for name, g in S.G_iw: g << gpara

# Compute G0
for name, g0 in S.G0_iw:
  g0 << inverse( iOmega_n + mu - (half_bandwidth/2.0)**2  * S.G_iw[name] )

# Parameters
p = {}
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["length_cycle"] = 50
p["n_warmup_cycles"] = 5000
p["n_cycles"] = 500000
p["make_histograms"] = True
p["use_proposed"] = use_prop

S.solve(h_loc=H, **p)

# Calculation is done. Now save a few things
if use_prop:
    filename = "hight_with.h5"
else:
    filename = "hight_without.h5"

Results = HDFArchive(filename,'w')
Results["Sigma_iw"] = S.Sigma_iw
Results["G_tau"] = S.G_tau
Results["G_iw"] = S.G_iw
