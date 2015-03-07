"""
Test self-consistent iteration near an NS interface
"""
from __future__ import division
from testutils import *

from usadel1 import *
from numpy import *

__revision__ = "$Id: test_selfcons_ns.py 3184 2006-10-02 07:05:45Z pauli $"

def test():
    raise KnownFailure("Needs updating: semantics for self-consistent "
                       "iteration have changed")
#
# Let's aim to calculate the behavior of the order parameter in an
# S-wire connected to an N-terminal, when there is a potential
# difference over the junction:
#
    Delta = 60

    g = Geometry(2, 3)
    g.t_type = [NODE_CLEAN_S_TERMINAL,
                 NODE_CLEAN_NODE,
                 NODE_CLEAN_N_TERMINAL]
    g.t_delta = Delta
    g.t_phase = 0
    g.t_inelastic = 1e-9
    g.t_t = 1e-9
    g.t_mu = [0, 0, 1]
    g.w_type = [WIRE_TYPE_S, WIRE_TYPE_N]
    g.w_length = [0.5, 1]
    g.w_conductance = 1
    g.w_inelastic = 1e-8
    g.w_spinflip = 0
    g.w_ends[0,:] = [0,1]
    g.w_ends[1,:] = [2,1]
#
# We need to fill in these parameters:
#
    g.omega_D = 4000
    g.coupling_lambda = 1/arccosh(g.omega_D/Delta)
#
# Also, put in an initial guess for the wire delta and phase:
#
    g.w_phase = 0
    g.w_delta = Delta
#
# So, solve the problem
#
    #sc = SelfConsistentIteration(g, 'test.h5', maxE=200)
    #sc.set_solvers(kin_solver=KIN_SOLVER_BLOCK)
    #sc.iterate()
#
# Convergence probably indicates that everything went ok...
#
# FIXME: this test should be more rigorous

if __name__ == "__main__":
    run_tests()
