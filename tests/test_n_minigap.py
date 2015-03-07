"""
Test the minigap in an N layer

Check that the minigap in an N layer matches analytical results.
"""

from __future__ import division
from testutils import *

from scipy import *
from scipy import optimize, integrate
import usadel1 as u

__revision__ = "$Id: test_n_minigap.py 3184 2006-10-02 07:05:45Z pauli $"

# 
# An implicit analytic expression for the minigap in an N-layer atop a
# superconductor exists:
# 
# .. latex::
#    \begin{align}
#    \sqrt{E_g}/L_N &= \max_{y<z} \frac{1}{2\sqrt{2}}
#                   \int_y^z \frac{d x}{\sqrt{(x^2 + 1/4)(x - y)}}
#    \,, \\
#    z &= -\frac{1}{2} \sinh {\rm Re}\,\theta_0 \,,
#    \end{align}
# 
# Where :m:`\theta_0` is the value of :m:`\theta` on the NS interface.
# If we assume that it is the bulk value
# 
# .. latex::
#    \[
#    \theta_0 = {\rm artanh}\, \frac{\lvert\Delta\rvert}{E_g + 0 i}
#    \]
# 
# then it is straightforward to evaluate it numerically:
# 

def exact_minigap(Delta):
    def integrand(x, y):
        return 1/sqrt((x**2 + .25) * (x - y))
    def func_1(y, z):
        fval,err = integrate.quad(integrand, y, z, args=(y,))
        return -fval / (2*sqrt(2))
    def func_2(E_g):
        theta_0 = arctanh(Delta / (E_g + 0j))
        z = -.5*sinh(theta_0.real)
        yopt, fval, ierr, numfunc = optimize.fminbound(
            func_1, -1 - 5*abs(z), z, args=(z,), full_output=True)
        return E_g - fval**2

    E_g = optimize.fsolve(func_2, min(0.5*Delta, 0.7))
    return E_g

#     
# Let's check that we get the same results from the code
# 

def test_bilayer():
    solver = u.Solver()
    for Delta in [0.1, 1, 10, 50]:
        E_g = exact_minigap(Delta)

        g = get_geometry(Delta)
        E = linspace(0, 2*E_g, 200)
        solver.set_geometry(g)
        solver.set_solvers(sp_solver=u.SP_SOLVER_COLNEW)
        solution = solver.sp_solve(E, g.x)

        dos = ((1 + solution.a*solution.b)/(1 - solution.a*solution.b)).real
        gap = all(dos[:,0,:] < 1e-5, axis=1)
        E_g_num = max(E[gap])

        assert allclose(E_g, E_g_num, rtol=1e-8, atol=2*max(E)/len(E)),\
               (E_g, E_g_num)

def get_geometry(Delta_0, x=None):
    g = u.Geometry(1, 2, x=x)

    g.t_type = [u.NODE_CLEAN_S_TERMINAL, u.NODE_FREE_INTERFACE]

    g.t_delta = [Delta_0, 0]
    g.t_phase = 0
    g.t_inelastic = 1e-9

    g.w_type = [ u.WIRE_TYPE_N ]
    g.w_length = 1
    g.w_conductance = 1

    g.w_ends[0,:] = [ 0, 1 ]

    return g

if __name__ == "__main__":
    run_tests()
