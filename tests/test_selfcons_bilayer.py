"""
Test self-consistent iteration of a bilayer structure

Verify the averaging of coupling coefficients in an NS bilayer.
"""
from __future__ import division

from testutils import *
from scipy import *
from scipy import optimize, integrate
import usadel1 as u

__revision__ = "$Id: test_selfcons_bilayer.py 3196 2006-10-04 10:31:09Z pauli $"

def setup():
     global g, Delta_0, omega_D, lambda_0, d_A, d_B

     #raise RuntimeError("This test does not work yet.")

     d_A = 0.01
     d_B = d_A / 4

     Delta_0 = 1
     omega_D = Delta_0 * 125 # 400 K
     lambda_0 = 1/log(2*omega_D/Delta_0)

     x = linspace(0, 1, 200)
     g = get_geometry(d_A, d_B, lambda_0, omega_D, 0.0002, x=x)

@slow
def test_bilayer_realtime():
     global g, Delta_0, omega_D, lambda_0, d_A, d_B

     g.w_delta[...] = Delta_0

     Delta_1 = bilayer_delta(lambda_0, omega_D, d_A, d_B)

     print Delta_1

     # Note: need the whole energy range, because the Thouless energy is much
     #       larger than the enery gap and omega_D!
     E = u.linspace_weighed(1e-4, omega_D, 1493, [(0, 10, 0.15),
                                                  (0, 10, 5)])

     currents = u.CurrentSolver(g, E=E, chunksize=1500)
     currents.set_solvers(sp_solver=u.SP_SOLVER_TWPBVP)

     it = u.self_consistent_realtime_iteration(currents)

     for k, d, violation in it:
          print "Iter %d: relative residual %.2g" % (
               k, d.relative_residual_norm())

          if (d.relative_residual_norm() < 1e-4 and violation < 1e-5):
               break
     else:
          raise RuntimeError("Did not converge")

     assert allclose(g.w_delta[0,:], Delta_1, rtol=1e-2), \
            (g.w_delta[0,0], Delta_1)

def test_bilayer_matsubara():
     global g, Delta_0, omega_D, lambda_0, d_A, d_B

     g.w_delta[...] = Delta_0

     Delta_1 = bilayer_delta(lambda_0, omega_D, d_A, d_B)

     # Note: again need the whole energy range up to omega_D,
     #       because the Thouless energy is much larger than the enery gap
     #       and omega_D!
     it = u.self_consistent_matsubara_iteration(g, max_ne=200)

     for k, d, violation in it:
          print "Iter %d: relative residual %.2g" % (
               k, d.relative_residual_norm())

          if (d.relative_residual_norm() < 1e-4 and violation < 1e-5):
               break
     else:
          raise RuntimeError("Did not converge")

     assert allclose(g.w_delta[0,:], Delta_1, rtol=1e-2), g.w_delta[0,0]

def bilayer_delta(lambda_0, omega_D, d_S, d_N):
     """
     Compute the Delta(T=0) in the S part of a extremely thin SN bilayer;
     d_S, d_N << coherence length.

     """

     def delta_avg_fun(delta):
          delta = abs(delta)
          c = d_A/(d_A + d_B)
          ii, ier = integrate.quad(
               lambda ee: 1/sqrt((ee+1e-9j)**2 - c**2 * delta**2).real,
               c*delta, omega_D)
          res = 1 - lambda_0*c*ii
          return res

     Delta_1, info, ier, mesg = optimize.fsolve(delta_avg_fun, Delta_0/2,
                                                full_output=1)
     assert ier == 1
     return abs(Delta_1)

def get_geometry(L_1, L_2, lambda_0, omega_D, T0, x=None):
     g = u.Geometry(2, 3, x=x)

     g.t_type = [u.NODE_FREE_INTERFACE,
                 u.NODE_FREE_INTERFACE,
                 u.NODE_CLEAN_NODE]

     g.t_delta = 0
     g.t_phase = 0

     g.omega_D = omega_D

     g.coupling_lambda = lambda_0
     g.w_phase = 0

     g.t_inelastic = 1e-3
     g.t_spinflip = 0
     g.t_t = T0
     g.t_mu = 0

     g.w_type = [u.WIRE_TYPE_S, u.WIRE_TYPE_N]
     g.w_length = [L_1, L_2]
     g.w_conductance = 1
     g.w_inelastic = 1e-3
     g.w_spinflip = 0

     g.w_ends[0,:] = [ 0, 2 ]
     g.w_ends[1,:] = [ 1, 2 ]
     g.w_delta = 2 * omega_D * exp(-1/lambda_0)

     return g

if __name__ == "__main__":
    run_tests()
