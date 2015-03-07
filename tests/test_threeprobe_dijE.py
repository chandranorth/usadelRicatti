"""
Regression test in a three-probe structure

A problem may occur when solving $j_E$: at least TWPBVP introduces
quite a large numerical noise in the solution, which can prevent the
problem from being solved.
"""
from __future__ import division
from testutils import *

import scipy.interpolate
import scipy.integrate
import usadel1

from usadel1 import *
from usadel1.data import *
from numpy import *

__revision__ = "$Id: test_threeprobe_dijE.py 3184 2006-10-02 07:05:45Z pauli $"

#
# Initialize:
#
def setup():
     global solver, geometry, E, x
     E = linspace_weighed(0, 10, 200, ((.10, 1, 30),))
     x = linspace(0, 1, 101)

     geometry = geometry_NS3(x, 99, 1e-6, 1.27*pi/2, 1e9)
     solver = Solver()
     solver.set_geometry(geometry)
     solver.set_solvers(kin_solver=KIN_SOLVER_TWPBVP)

def teardown():
    global solver, geometry
    solver = None
    geometry = None

#
# Test that COLNEW produces reasonable results:
#
def test_NS3_COLNEW():
     solver.set_solvers(sp_solver=KIN_SOLVER_COLNEW)
     spectral = solver.sp_solve(E, x)

     coefficient = KineticCoefficients(geometry, spectral)
     solver.set_kinetic(coefficient)
     kinetic = solver.kin_solve(E, x)

     assert is_zero(array_normmax(coefficient.dijE)
                    / array_normmax(coefficient.ijE), 1e-3)

#
# Test that TWPBVP (actually TWPBVP) produces reasonable results:
#
def test_NS3_TWPBVP():
     solver.set_solvers(sp_solver=KIN_SOLVER_TWPBVP)
     spectral = solver.sp_solve(E, x)

     coefficient = KineticCoefficients(geometry, spectral)
     solver.set_kinetic(coefficient)
     kinetic = solver.kin_solve(E, x)

     assert is_zero(array_normmax(coefficient.dijE)
                    / array_normmax(coefficient.ijE), 1e-3)

#
# The Geometry to use:
#
def geometry_NS3(x, V0, T0, phi, Delta):
     """N-wire between S terminals."""
     g = Geometry(3, 4, x)

     g.t_type = [ NODE_CLEAN_S_TERMINAL,
                  NODE_CLEAN_S_TERMINAL,
                  NODE_CLEAN_N_TERMINAL,
                  NODE_CLEAN_NODE ]

     g.t_delta = [ Delta, Delta, 0, 0 ]
     g.t_phase = [ 0, phi, 0, 0 ]

     g.t_inelastic = 1e-9
     g.t_spinflip = 0
     g.t_t = T0
     g.t_mu = [0, 0, V0, 0]

     g.w_type = [ WIRE_TYPE_N, WIRE_TYPE_N, WIRE_TYPE_N ]
     g.w_length = 1
     g.w_conductance = 1
     g.w_inelastic = 1e-9
     g.w_spinflip = 0

     g.w_ends[0,:] = [ 0, 3 ]
     g.w_ends[1,:] = [ 1, 3 ]
     g.w_ends[2,:] = [ 2, 3 ]

     return g

if __name__ == "__main__":
    run_tests()
