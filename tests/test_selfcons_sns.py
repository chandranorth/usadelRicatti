"""
Test self-consistent iteration for an SNS structure; esp. check current
conservation.

"""
from __future__ import division
from testutils import *

from usadel1 import *
from numpy import *

__revision__ = "$Id: test_selfcons_sns.py 3184 2006-10-02 07:05:45Z pauli $"

@slow
def test():
#
# Let's aim to calculate the behavior of the order parameter in an SNS
# junction:
#
    Delta = 60

    x = linspace(0, 1, 300);
    g = Geometry(1, 2, x=x)
    g.t_type = [NODE_CLEAN_S_TERMINAL,
                 NODE_CLEAN_S_TERMINAL]
    g.t_delta = Delta
    g.t_phase = array([-.5, .5]) * pi/2
    g.t_inelastic = 1e-9
    g.t_t = 1e-9
    g.t_mu = 0
    g.w_type = [WIRE_TYPE_S]
    g.w_length = [1]
    g.w_conductance = 1
    g.w_inelastic = 1e-8
    g.w_spinflip = 0
    g.w_ends[0,:] = [0,1]
#
# We need to fill in these parameters:
#
    LS = 0.9
    x1 = LS / 2
    x2 = 1 - LS / 2

    g.omega_D = 4000
    g.coupling_lambda = 1/arccosh(g.omega_D/Delta)
    g.coupling_lambda[0, (g.x > x1) & (g.x < x2)] = 0
#
# Also, put in some initial guesses for the wire delta and phase:
#
    g.w_phase[0, g.x < x1] = g.t_phase[0]
    g.w_phase[0, g.x > x2] = g.t_phase[1]
    g.w_delta[0, (g.x < x1) | (g.x > x2)] = Delta
#
# So, solve the problem
#
    sol = CurrentSolver(g, maxE=200, ne=300)
    it = self_consistent_realtime_iteration(sol)
    for k, d, v in it:
        print "Iteration %d: residual %g (current-conservation violation %g)"%(
            k, d.residual_norm(), v)
        if d.residual_norm() < 1e-3 and v < 5e-3:
            break
    else:
        raise RuntimeError("Self-consistent iteration did not converge")
#
# Convergence probably indicates that everything went ok...
#
# FIXME: this test should be more rigorous

if __name__ == "__main__":
    run_tests()
