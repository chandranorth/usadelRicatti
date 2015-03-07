"""
Testing current in a single wire
"""
from __future__ import division
from testutils import *

from numpy import *
from usadel1 import *
from geometries import *
import scipy.interpolate
import scipy.integrate

__revision__ = "$Id: test_singlewire.py 3184 2006-10-02 07:05:45Z pauli $"

# 
# Solve a single N-wire between two S-terminals, first treating the wire
# as one part, then splitting it into two. This should test whether
# scaling of currents etc. quantities goes correctly.
# 
# Below, we want to evaluate the supercurrent at one point. It fastest
# to do this in a specialized way, since solutions to the kinetic
# equations are known:
# 
def evaluate_IS(coef, T0):
     # Use Scipy's spline routines -- they seem to be fast enough
     rep = scipy.interpolate.splrep(coef.E,coef.ijE[:,0,0],k=3,s=0)
     def fun(E):
         return scipy.interpolate.splev(E, rep) * tanh(.5*E/T0)
     result = scipy.integrate.quad(fun, 1e-10, max(coef.E),
                                   points=(2*T0,), full_output=True)
     return result[0]

# 
# Then, to the calculation:
# 
def test():

    E = linspace_weighed(1e-7, 1000, 1000, ((10, 1, 20),))
    x = linspace(0, 1, 101)
    phi = 1.27*pi/2
    gs = [ geometry_SNS(x, 0, 1e-6, phi, 1e9),
            geometry_SNS_split(x, 0, 1e-6, phi, 1e9) ]

    T0s = ([1e-6, 1e-5, 1e-4, 1e-3 ]
         + list(arange(0.01, 0.09, 0.01))
         + list(arange(0.1, 0.9, 0.1))
         + list(arange(1, 9, 0.25))
         + list(arange(10, 100, 1)))

    data = [[], []]
    for ig, geometry in enumerate(gs):
         solver = CurrentSolver(geometry, E, chunksize=250)
         #solver.solver.set_solvers(sp_solver=SP_SOLVER_TWPBVP)
         solver.solve_spectral()
         
         for T0 in T0s:
             IS = evaluate_IS(solver.coefficient, T0)
             data[ig].append((phi, T0, IS))
# 
# Then check that the results are equal. They of course should be, but
# this test should catch any problems with scaling of energies etc.
# 
    assert arrays_equal(array(data[0]), array(data[1]),
                         tolerance=1e-5)
# 
# Compare the maximum supercurrent to a result in [dubos01]_:
# 
    IS = scipy.interpolate.interp1d(array(data[0])[:,1],
                                     array(data[0])[:,2])
        
    assert is_zero(IS(1e-3) - 10.82, tolerance=0.02)
    assert is_zero(IS(25) - 0.04, tolerance=0.005)
# 
# 
# .. [dubos01] P. Dubos et al., Phys. Rev. B 63, 064502 (2001).

if __name__ == "__main__":
    run_tests()
