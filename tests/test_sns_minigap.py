"""
Test the minigap in an SNS junction
"""
from __future__ import division
from testutils import *
from scipy import *
import usadel1 as u

__revision__ = "$Id: test_sns_minigap.py 3185 2006-10-02 07:07:23Z pauli $"

#
# Approximate analytic expressions for the minigap in an SNS junction
# are discussed in [ivanov02]:
#

C_1 = 0.0921
C_2 = 3.122
C_3 = pi**2 / 4

def approx_minigap_0(phi):
    """Approximate minigap at phi = 0"""
    return C_2 * (1 - C_1 * phi**2)

def approx_minigap_pi(phi):
    """Approximate minigap at phi = pi"""
    return C_3 * (pi - phi)

#
# Let's check that we get the same results from the code:
#
@slow
def test_minigap():
    solver = u.Solver()
    #solver.set_solvers(sp_solver=u.SP_SOLVER_TWPBVP)
    Delta = 1e6

    for phi in [0.01, 0.1, 0.3, 0.5, pi - 0.08, pi - 0.01]:
        g = get_geometry(Delta, phi)

        if phi > pi/2:
            E_g = approx_minigap_pi(phi)
        else:
            E_g = approx_minigap_0(phi)

        Egg = linspace(0.99*E_g, 1.01*E_g, 250)
        E = unique(r_[linspace(1e-4, 0.99*E_g, 60), Egg])
        solver.set_geometry(g)
        solver.set_solvers(sp_solver=u.SP_SOLVER_COLNEW, sp_tol=1e-8)
        s = solver.sp_solve(E, g.x)

        dos = ((1 + s.a*s.b)/(1 - s.a*s.b)).real
        gap = all(dos[:,0,:] < 1e-5, axis=1)
        E_g_num = max(E[gap])

        print '%', E_g, E_g_num, " __ ", min(Egg), '..', max(Egg)

        #import matplotlib.pyplot as plt
        #plt.plot(E, dos[:,0,50])
        #plt.show()

        if phi < pi/2:
            assert allclose(E_g, E_g_num, rtol=5e-4, atol=2*diff(Egg).max())
        else:
            assert allclose(E_g, E_g_num, rtol=1e-1, atol=2*diff(Egg).max())

def get_geometry(Delta, phi):
    g = u.Geometry(1, 2)

    g.t_type = [u.NODE_CLEAN_S_TERMINAL, u.NODE_CLEAN_S_TERMINAL]

    g.t_delta = [Delta, Delta]
    g.t_phase = [-phi/2, phi/2]
    g.t_inelastic = 1e-9

    g.w_type = [ u.WIRE_TYPE_N ]
    g.w_length = 1
    g.w_conductance = 1

    g.w_ends[0,:] = [ 0, 1 ]

    return g

#
# .. [ivanov02]
#    D. A. Ivanov, R. von Roten, G. Blatter.
#    "Minigap in a long disordered SNS junction: Analytical results".
#    Phys. Rev. B. *66*, 052507 (2002).
#

if __name__ == "__main__":
    run_tests()
