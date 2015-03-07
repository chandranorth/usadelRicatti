"""
Test that the phase jumps work in simple cases.

"""
from __future__ import division
from testutils import *
from numpy import *
from usadel1 import *
import tables

#
#
#

def get_geometry(phi, Phi):
    x = linspace(0, 1, 71)
    g = Geometry(4, 4, x)

    g.t_type = [NODE_CLEAN_S_TERMINAL,
                NODE_CLEAN_NODE,
                NODE_CLEAN_NODE,
                NODE_CLEAN_S_TERMINAL]
    
    g.t_delta = [ 50, 0, 0, 50]
    g.t_phase = array([ -.5, .5, 0, 0 ]) * phi
    
    g.t_inelastic = 1e-3
    g.t_spinflip = 0
    g.t_t = 1e-6
    g.t_mu = 0
    
    g.w_type = WIRE_TYPE_N
    g.w_length = [ 1./3, 1./3, 1./3, 1./3 ]
    g.w_conductance = [ 1, 1, 1, 1 ]

    g.w_phase_jump[1] =  Phi/2
    g.w_phase_jump[2] = -Phi/2

    g.w_ends[0,:] = [ 0, 1 ]
    g.w_ends[1,:] = [ 1, 2 ]
    g.w_ends[2,:] = [ 1, 2 ]
    g.w_ends[3,:] = [ 2, 3 ]

#
# Then solve the problem:
#

def test_phase_shift():

    g = Geometry(1, 2)

    g.t_type  = [NODE_CLEAN_S_TERMINAL, NODE_CLEAN_S_TERMINAL]
    g.t_delta = [50, 50]
    g.t_inelastic = 3e-3
    g.t_t = 1e-6

    g.w_type = WIRE_TYPE_N
    g.w_length = 1
    g.w_conductance = 1

    g.w_ends[0,:] = [0, 1]

    currents = CurrentSolver(g)
    currents.set_solvers(sp_solver=SP_SOLVER_TWPBVP)
    
    # 1.

    g.t_phase = [0, pi/2]
    g.w_phase_jump = 0

    currents.solve()
    Ic_1, Ie_1 = currents.get_currents(ix=0)

    a_1 = currents.spectral.a.copy()

    # 2.

    g.t_phase = [0, 0]
    g.w_phase_jump = pi/2

    currents.solve()
    Ic_2, Ie_2 = currents.get_currents(ix=0)

    a_2 = currents.spectral.a.copy()

    # Compare

    assert allclose(Ic_1, Ic_2), (Ic_1, Ic_2)

if __name__ == "__main__":
    run_tests()
