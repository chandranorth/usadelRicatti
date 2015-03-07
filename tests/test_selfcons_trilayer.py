"""
Test self-consistent iteration of a trilayer structure
"""
from __future__ import division
from testutils import *
from scipy import *
import usadel1 as u

__revision__ = "$Id: test_selfcons_trilayer.py 3196 2006-10-04 10:31:09Z pauli $"

def setup():
     global g

     d_A = 90
     d_B = 58
     d_C = 30

     Delta_0 = 3.17 * 10
     omega_D = Delta_0 * 125 # 400 K
     lambda_0 = 1/arccosh(omega_D/Delta_0)

     lambda_1 = d_B*lambda_0/(d_A+d_B+d_C)
     lambda_2 = d_B*lambda_0/(    d_B+d_C)
     lambda_3 = d_B*lambda_0/(    d_B    )

     L1 = 2
     L2 = 1.06
     L3 = 1

     x = linspace(0, 1, 200)
     g = get_geometry(d_A, d_B, d_C, lambda_1, lambda_2, lambda_3,
                      omega_D, 1e-6, x=x)

def test_trilayer():
     global g

     raise KnownFailure("Needs updating: semantics for self-consistent "
                        "iteration have changed. Moreover, this test was "
                        "completely written...")

     #solver = u.SelfConsistentIteration.resume(g, 'test.h5',
     #                                          ne=300, maxE=60)
     #solver.iterate()

def get_geometry(L_1, L_2, L_3, lambda_1, lambda_2, lambda_3,
                  omega_D, T0, x=None):
     g = u.Geometry(1, 2, x=x)

     g.t_type = [ u.NODE_CLEAN_S_TERMINAL,
                  u.NODE_FREE_INTERFACE ]

     Delta_1 = omega_D/cosh(1/lambda_1)
     Delta_2 = omega_D/cosh(1/lambda_2)
     Delta_3 = omega_D/cosh(1/lambda_3)
     g.t_delta = array([ Delta_1, 0 ])
     g.t_phase = array([ 0, 0 ])

     g.omega_D[:,:] = omega_D

     j1 = int(where(g.x > L_1/(L_1+L_2+L_3))[0][0])
     j2 = int(where(g.x > (L_1+L_2)/(L_1+L_2+L_3))[0][0])

     g.coupling_lambda[0,0:j1]  = lambda_1
     g.coupling_lambda[0,j1:j2] = lambda_2
     g.coupling_lambda[0,j2:]   = lambda_3

     g.w_phase[:,:] = 0

     g.t_inelastic[:] = 1e-9
     g.t_spinflip[:] = 0
     g.t_t[:] = T0
     g.t_mu[:] = 0

     g.w_type = [ u.WIRE_TYPE_S ]
     g.w_length = array([ L_1 + L_2 + L_3 ])
     g.w_conductance[:] = 1
     g.w_inelastic[:] = 0
     g.w_spinflip[:] = 0

     g.w_ends[0,:] = [ 0, 1 ]
     g.w_delta[:,:] = Delta_3

     return g

if __name__ == "__main__":
    run_tests()
