"""
Test self-consistent iteration in simple cases

:Author: Pauli Virtanen <pauli@ltl.tkk.fi>
"""
from __future__ import division
from testutils import *
from usadel1 import *
from numpy import *

__revision__ = "$Id: test_selfcons_simple.py 3184 2006-10-02 07:05:45Z pauli $"

def setup():
    global g, Delta
#
# Bulk
# ----
#
# Let's aim to calculate the following value for the energy gap:
#
    Delta = 60
#
# in a single S wire connected between N and S terminals:
#
    g = Geometry(1, 2)
    g.t_type = [NODE_CLEAN_S_TERMINAL, NODE_CLEAN_N_TERMINAL]
    g.t_delta = Delta
    g.t_inelastic = 1e-9
    g.t_t = 0.4*Delta
    g.t_mu = 0
    g.w_type = [WIRE_TYPE_S]
    g.w_length = 1
    g.w_conductance = 1
    g.w_inelastic = 1e-9
    g.w_spinflip = 0
    g.w_ends[0,:] = [0,1]
#
# We need to fill in these parameters:
#
    g.omega_D = 4000
    g.coupling_lambda = 1/arccosh(g.omega_D/Delta)
#
# Also, put in some nasty initial guess for the wire delta and phase:
#
    g.w_phase = 0
    g.w_delta = 10
#
# So, solve the problem
#

def test():
    global g, Delta

    g.w_phase = 0
    g.w_delta = 10
    sol = CurrentSolver(g, ne=400, chunksize=400)
    sol.set_solvers(sp_solver=SP_SOLVER_COLNEW)
    it = self_consistent_realtime_iteration(sol)
    for k, d, v in it:
        print "Residual:", d.residual_norm()
        if d.residual_norm() < 1e-3:
            break
    else:
        raise RuntimeError("Iteration did not converge")

    realtime_delta = g.w_delta.copy()

    g.w_phase = 0
    g.w_delta = 10
    it = self_consistent_matsubara_iteration(g)
    for k, d, v in it:
        print "Residual:", d.residual_norm()
        if d.residual_norm() < 1e-3:
            break
    else:
        raise RuntimeError("Iteration did not converge")

    matsubara_delta = g.w_delta.copy()

    print around(realtime_delta, 3)
    print around(matsubara_delta, 3)

    assert allclose(realtime_delta, matsubara_delta, atol=0.5, rtol=5e-2)

    assert is_zero(abs(realtime_delta[0,0]) - Delta, tolerance=5)
    assert is_zero(abs(realtime_delta[0,-1]), tolerance=0.5)

    assert is_zero(abs(matsubara_delta[0,0]) - Delta, tolerance=5)
    assert is_zero(abs(matsubara_delta[0,-1]), tolerance=0.5)

if __name__ == "__main__":
    run_tests()
