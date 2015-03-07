"""
Spectral supercurrent in a three-probe structure
"""
from __future__ import division
from testutils import *
from numpy import *
from usadel1 import *
import tables

__revision__ = "$Id: test_threeprobe.py 3261 2007-01-23 13:15:06Z pauli $"

#
# Let's calculate the spectral supercurrent flowing in a three-probe
# structure::
#
#
#        $j_S$
#      -------->
#
#     S----+----S
#          |
#          |
#          |
#          N
#
# As always, this starts by specifying a geometry.
#

def setup():
    global solver, g

    x = linspace(0, 1, 71)
    g = Geometry(3, 4, x)
#
# (Specifying ``x`` is not mandatory, but let's do it nonetheless.)
#
    g.t_type = [NODE_CLEAN_S_TERMINAL,
                 NODE_CLEAN_S_TERMINAL,
                 NODE_CLEAN_N_TERMINAL,
                 NODE_CLEAN_NODE]

    g.t_delta = [ 50, 50, 0, 0 ]
    g.t_phase = array([ -.5, .5, 0, 0 ]) * (pi/2)

    g.t_inelastic = 1e-9
    g.t_spinflip = 0
    g.t_t = 1e-6
    g.t_mu = 0

    g.w_type = [ WIRE_TYPE_N, WIRE_TYPE_N, WIRE_TYPE_N ]
    g.w_length = [ 0.5, 0.5, 1.35 ]
    g.w_conductance = [ 1, 1, 1.2 ]
    g.w_inelastic = array([ 1e-9, 1e-9, 1e-9 ])
    g.w_spinflip = [ 0, 0, 0 ]

    g.w_ends[0,:] = [ 0, 3 ]
    g.w_ends[1,:] = [ 1, 3 ]
    g.w_ends[2,:] = [ 2, 3 ]

    E = linspace_weighed(0, 100, 350, ((0, 1, 50),))
    solver = CurrentSolver(g, E, chunksize=100)

def teardown():
    global solver, g
    solver = None
    g = None

#
# Then solve the problem:
#
def test_1():
    global solver, g
    solver.solve()
#
# Specifying ``E`` is also not mandatory -- but we do it here to
# reduce the number of points in energy discretization.
#
# Now, the data is in ``solver.kinetic``, ``solver.spectral`` and
# ``solver.coefficient``.
#
#
# Checking the results
# --------------------
#
# Check now the results for $j_S$ against known-good data:
#
    jS = load_ascii_array('test_threeprobe.1.dat')
    is_zero(difference(jS), tolerance=5e-2)

def difference(jS):
     global solver, g
     return curve_max_difference(4*jS[:,4], 2*jS[:,8],
         solver.coefficient.E[:], solver.coefficient.ijE[:,0,50])

#
# The multiplication factors are needed since the known-good data is
# scaled differently from this data.


#
# Trying other solvers
# --------------------
#
# The previous run was done using the default solver engine. Let's
# now do runs with other engines:
#
def test_twpbvp_1():
     global solver, g
     solver.set_solvers(sp_solver=SP_SOLVER_TWPBVP)
     solver.solve()
     jS = load_ascii_array('test_threeprobe.1.dat')
     assert is_zero(difference(jS), tolerance=5e-2)

def test_colnew_1():
     global solver, g
     solver.set_solvers(sp_solver=SP_SOLVER_COLNEW)
     solver.solve()
     jS = load_ascii_array('test_threeprobe.1.dat')
     assert is_zero(difference(jS), tolerance=5e-2)

#
# Different geometry
# ------------------
#
# Let's change the geometry slightly:
#
def test_other_geometry():
    global solver, g
    g.t_delta = array([ 20, 20, 0, 0 ])
    g.w_length = [ 0.5, 0.5, 0.55 ]
    g.w_conductance = [ 1, 1, 2.6 ]
    solver.set_solvers(sp_solver=SP_SOLVER_DEFAULT)
    solver.solve()
    jS2 = load_ascii_array('test_threeprobe.2.dat')
    assert is_zero(difference(jS2), tolerance=1e-1)

#
# Quasiparticle current
#

def test_summation():
    global solver, g
    g.t_mu[0:3] = [0, 0, 30]
    g.t_t[0:3] = [50, 50, 50]
    g.t_phase = array([ -.5, .5, 0, 0 ]) * pi * 0.95
    solver.set_solvers(sp_solver=SP_SOLVER_TWPBVP,
                       kin_solver=KIN_SOLVER_BLOCK)
    solver.solve()

    print solver.kinetic.jT[20,:,1]*g.w_conductance
    print solver.coefficient.ijE[20,:,1]

    # FIXME: no asserts

def test_supercurrent():
    global solver, g

    import tables

    phis = r_[0.5*pi, 0.99*pi]
    g.t_mu[0:3] = [0,0,5]

    for i, phi in enumerate(phis):
        g.t_phase[0:2] = [-.5*phi,.5*phi]
        solver.solve()
        Ic, IE = solver.get_currents(ix=0)

    # FIXME: no asserts

if __name__ == "__main__":
    run_tests()
