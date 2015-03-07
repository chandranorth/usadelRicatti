"""
Test the simplest (and only!) charge imbalance-dependent boundary
conditions that we have available.

:Author: Pauli Virtanen <pauli@ltl.tkk.fi>
"""
from __future__ import division
from testutils import *

from numpy import *
from usadel1 import *
from geometries import *
import sys

__revision__ = "$Id: test_simple_charge_imbalance_bc.py 3262 2007-01-23 14:13:51Z pauli $"

#
# First, specify the geometry
#
def setup_module():
    global g1, g2, s1, s2

    g1 = geometry_SnnNnNnnS(phi=1)
    g2 = geometry_SnnNnNnnS(phi=1)

    g1.t_type = [ NODE_CLEAN_N_TERMINAL,
                  NODE_CLEAN_N_TERMINAL,
                  NODE_CLEAN_S_TERMINAL,
                  NODE_CLEAN_S_TERMINAL,
                  NODE_CLEAN_NODE,
                  NODE_CLEAN_NODE ]
    g1.w_length = [1,1,1,1,1]
    g1.w_conductance = 1
    g1.t_delta = array([ 0, 0, 60, 60, 0, 0 ])

    g2.t_type = [ NODE_CLEAN_N_TERMINAL,
                  NODE_CLEAN_N_TERMINAL,
                  NODE_CLEAN_S_TERMINAL_CIB,
                  NODE_CLEAN_S_TERMINAL_CIB,
                  NODE_CLEAN_NODE,
                  NODE_CLEAN_NODE ]
    g2.w_length = [1,1,1,1,1]
    g2.w_conductance = 1
    g2.t_delta = array([ 0, 0, 60, 60, 0, 0 ])

    s1 = CurrentSolver(g1, maxE=450)
    s2 = CurrentSolver(g2, maxE=450)

    s1.set_solvers(kin_solver=KIN_SOLVER_BLOCK)
    s2.set_solvers(kin_solver=KIN_SOLVER_BLOCK)

    s1.solve_spectral_if_needed(calculate_G=True)
    s2.solve_spectral_if_needed(calculate_G=True)

def teardown_module():
    global g1, g2, s1, s2
    g1 = None
    g2 = None
    s1 = None
    s2 = None

#
# Calculate NN conductance at *zero* temperature without charge imbalance
#
@slow
def test_zero_temp():
    global g1, g2, s1, s2

    g1.t_t = 1e-2
    g2.t_t = 1e-2

    G1_VV, G1_VT, G1_TV, G1_TT = s1.get_linear_response_from_G()
    G2_VV, G2_VT, G2_TV, G2_TT = s2.get_linear_response_from_G()

    R1 = linalg.pinv(G1_VV[0:2,0:2])
    R2 = linalg.pinv(G2_VV[0:2,0:2])

    RNN1 = dot([[1, -1]], dot(R1, [[1], [1]]))
    RNN2 = dot([[1, -1]], dot(R2, [[1], [1]]))

    R_N_without_loop = 1 + 1 + 1
    R_N_with_loop    = 1 + 1 + 1/(1/2 + 1/1)

    print R_N_with_loop, "~=", RNN1, R_N_with_loop, "~=", RNN2

    assert allclose(RNN1, R_N_with_loop, rtol=1e-2), RNN1
    assert allclose(RNN2, R_N_with_loop, rtol=1e-2), RNN2

#
# Calculate NN conductance at *zero* temperature without charge imbalance
#
@slow
def test_large_temp():
    global g1, g2, s1, s2

    g1.t_t = 100
    g2.t_t = 100

    G1_VV, G1_VT, G1_TV, G1_TT = s1.get_linear_response_from_G()
    G2_VV, G2_VT, G2_TV, G2_TT = s2.get_linear_response_from_G()

    R1 = linalg.pinv(G1_VV[0:2,0:2])
    R2 = linalg.pinv(G2_VV[0:2,0:2])

    RNN1 = dot([[1, -1]], dot(R1, [[1], [1]]))
    RNN2 = dot([[1, -1]], dot(R2, [[1], [1]]))

    R_N_without_loop = 1 + 1 + 1
    R_N_with_loop    = 1 + 1 + 1/(1/2 + 1/1)

    print R_N_with_loop, "~=", RNN1, R_N_without_loop, "~=", RNN2

    assert allclose(RNN1, R_N_with_loop, rtol=3e-2), RNN1
    assert allclose(RNN2, R_N_without_loop, rtol=3e-2), RNN2

if __name__ == "__main__":
    run_tests()
