"""
Test spectral thermoelectric matrix

:Author: Pauli Virtanen <pauli@ltl.tkk.fi>
"""
from __future__ import division
from testutils import *
from numpy import *
from usadel1 import *
from geometries import *

__revision__ = "$Id: test_spectral_thermoelectric_matrix.py 3194 2006-10-04 10:08:13Z pauli $"

#
# The idea of a spectral thermoelectric matrix is based on the following
# notion:
#
#   1. Kinetic equations are linear
#   2. Hence, the spectral currents can be expressed in the form
#
#      .. latex::
#         $$
#         j_{\alpha,i}(E) = \sum_{\beta, r} G_{\alpha,i}^{\beta,r}(E)
#                                               f_{\beta,r}(E)
#         $$
#
#      Here, :m:`G_{\alpha,i}^{\beta,r}(E)` is the spectral
#      thermoelectric matrix. The indices :m:`\alpha` and :m:`\beta`
#      label "L" and "T", :m:`r` labels terminals and :m:`i` wires.
#
#      The "diagonal" elements in the "L/T" space correspond to thermal
#      and electrical conduction, and the off-diagonal elements are related
#      to thermopower and Thompson effect.
#
# I have written some code to calculate the elements of this matrix.
# Let's test it here.
#

#
# Three-probe system
# ------------------
#
def test_three_probe():
#
# First, initialize the solver and the geometry
#
    g = geometry_SnnNnS(V0=0, T0=1e-6, phi=1)
    sc = CurrentSolver(g, chunksize=300)
#
# COLNEW still has problems solving the kinetic equations. Hence, use a
# different solver.
#
    sc.set_solvers(sp_solver=SP_SOLVER_COLNEW, kin_solver=KIN_SOLVER_TWPBVP,
                   sp_tol=1e-4)
#
# Then, just calculate the full thermoelectric matrix...
#
    sc.solve_spectral_if_needed(calculate_G=True)
#
# The following calculates the currents first directly, then using the G
# matrix and compares the results:
#
    def cmp_results():
         Ic, IE = sc.get_currents_from_G()
         Ic_0, IE_0 = sc.get_currents(E=sc.kinetic.E, ix=0)
         assert arrays_equal(Ic, Ic_0, tolerance=1e-2), [Ic, Ic_0]
         assert arrays_equal(IE, IE_0, tolerance=1e-2), [IE, IE_0]
         return Ic, IE
#
# Let's test this for some parameters
#
    g.t_t  = [1, 1, 1, 0]
    g.t_mu = [0, 0, 0, 0]
    Ic, IE = cmp_results()
    assert arrays_equal(Ic, [6.8, -6.8, 0.0], tolerance=1e-1)
    assert arrays_equal(IE, [0.0, 0.0, 0.0],  tolerance=1e-4)
#
    g.t_t  = [5, 5, 3, 0]
    g.t_mu = [0, 0, 0, 0]
    Ic, IE = cmp_results()
    assert arrays_equal(Ic, [4.2, -4.2, 0.0], tolerance=1e-1)
    assert arrays_equal(IE, [0.0, 0.0, 0.0],  tolerance=1e-4)
#
    g.t_t  = [1, 1, 1, 0]
    g.t_mu = [0, 0, 5, 0]
    Ic, IE = cmp_results()
    assert arrays_equal(Ic, [3.5, -2.5, -1.0], tolerance=1e-1)
    assert arrays_equal(IE, [0.0, 0.0, 0.0],   tolerance=1e-4)
#
    g.t_t  = [2, 3, 4, 0]
    g.t_mu = [0, 0, 6, 0]
    Ic, IE = cmp_results()
    assert arrays_equal(Ic, [2.4, -1.2, -1.2], tolerance=1e-1)
    assert arrays_equal(IE, [0.0, 0.0, 0.0],   tolerance=1e-4)
#
#
# Four-probe system & thermopower
# -------------------------------
#
@slow
def test_four_probe():
#
# Initialize spectral quantities: use an asymmetric geometry
#
    g = geometry_SnnNnNnnS(1e-6, 1e-6, 1e-6, 1)
    g.w_length = [1.1, 0.8, 0.7, 1.3, 1]
    g.w_length /= sum(g.w_length[3:6])
    g.w_conductance = [0.7, 1.2, 0.8, 0.9, 1.1]
    g.t_delta[0:4] = [0, 0, 50, 50]
    sc = CurrentSolver(g, chunksize=100, maxE=300)
    sc.set_solvers(kin_solver=KIN_SOLVER_BLOCK)
    sc.solve_spectral_if_needed(calculate_G=True)
#
# Also, let's compare to a different kinetic solver:
#
    sc.set_solvers(kin_solver=KIN_SOLVER_TWPBVP)
#
# Thermopower machinery:
#
    def optimize_thermopower(current_func):
        def zero_current():
            Ic, Ie = current_func(w_jT=[0,1], w_jL=[], ix=0)
            return [Ic[0], Ic[1]]
        def adjust_potentials(z):
            g.t_mu[0:2] = z[0:2]
        z0 = array([0, 0], float64)
        z = optimize_parameters_for(z0, zero_current, adjust_potentials,
                                    xtol=1e-3, epsfcn=1e-5)
        adjust_potentials(z)
        return z
#
# Another way to calculate linear-response thermopower is to obtain it
# directly from the G-matrix:
#
    def linear_response_thermopower_from_G_matrix():
        sc.geometry.t_mu = 0
        G_VV, G_VT, G_TV, G_TT = \
              sc.get_linear_response_from_G()

        avg_T = average(g.t_t[0:4])
        dT = asmatrix(g.t_t[0:4] - avg_T).T

        M1 = G_VT[ix_([0,1], [0,1,2,3])]
        M2 = G_VV[ix_([0,1], [0,1])]

        dV = -linalg.solve(M2, M1 * dT)

        g.t_mu[0:2] = ravel(dV)
        g.t_mu[2:] = 0

        return ravel(dV)
#
# Then test it:
#
    g.t_t[0:4] = [1, 1, 1, 1]
    mu1 = optimize_thermopower(sc.get_currents)
    mu2 = optimize_thermopower(sc.get_currents_from_G)
    mu3 = linear_response_thermopower_from_G_matrix()
    assert arrays_equal(mu1, [0.0, 0.0], tolerance=1e-5), locals()
    assert arrays_equal(mu2, [0.0, 0.0], tolerance=1e-5), locals()
    assert arrays_equal(mu3, [0.0, 0.0], tolerance=1e-5), locals()
#
    g.t_t[0:4] = [1.01, 0.99, 1, 1]
    mu1 = optimize_thermopower(sc.get_currents)
    mu2 = optimize_thermopower(sc.get_currents_from_G)
    mu3 = linear_response_thermopower_from_G_matrix()
    assert arrays_equal(mu1, [0.000991, 0.000683], relative=1e-2), locals()
    assert arrays_equal(mu2, [0.000991, 0.000683], relative=1e-2), locals()
    assert arrays_equal(mu3, [0.000991, 0.000683], relative=1e-1), locals()
#
    g.t_t[0:4] = [5.1, 5.9, 5, 5]
    mu1 = optimize_thermopower(sc.get_currents)
    mu2 = optimize_thermopower(sc.get_currents_from_G)
    mu3 = linear_response_thermopower_from_G_matrix()
    assert arrays_equal(mu1, [-0.00710, -0.00278], relative=1e-2), locals()
    assert arrays_equal(mu2, [-0.00710, -0.00278], relative=1e-2), locals()
    assert arrays_equal(mu3, [-0.00710, -0.00278], relative=1e-1), locals()
#
# For testing, change the solver.
#
    sc.set_solvers(kin_solver=KIN_SOLVER_BLOCK)
#
# One test at higher temperatures
#
    g.t_t[0:4] = [26, 24, 25, 25]
    mu1 = optimize_thermopower(sc.get_currents)
    mu2 = optimize_thermopower(sc.get_currents_from_G)
    mu3 = linear_response_thermopower_from_G_matrix()
    assert arrays_equal(mu1, [0.000717, -0.000103], relative=1e-2),locals()
    assert arrays_equal(mu2, [0.000717, -0.000103], relative=1e-2),locals()
    assert arrays_equal(mu3, [0.000717, -0.000103], relative=1e-2),locals()
#
# A test far from linear response
#
    g.t_t[0:4] = [26, 24, 1, 1]
    mu1 = optimize_thermopower(sc.get_currents)
    mu2 = optimize_thermopower(sc.get_currents_from_G)
    mu3 = linear_response_thermopower_from_G_matrix()
    assert arrays_equal(mu1, [0.000774, -0.000148], relative=1e-2),locals()
    assert arrays_equal(mu2, [0.000774, -0.000148], relative=1e-2),locals()
    assert arrays_equal(mu3, [0.000774, -0.000148], relative=1), locals()
#
# The linear response estimate is off by a factor of 2 or so.
#
# Test the above-gap contribution from superconductors:
#
    g.t_t[0:4] = [40, 40, 39, 41]
    mu1 = optimize_thermopower(sc.get_currents)
    mu2 = optimize_thermopower(sc.get_currents_from_G)
    mu3 = linear_response_thermopower_from_G_matrix()
    assert arrays_equal(mu1, [1.391e-5, 1.833e-5], relative=1e-2), locals()
    assert arrays_equal(mu2, [1.391e-5, 1.833e-5], relative=1e-2), locals()
    assert arrays_equal(mu3, [1.391e-5, 1.833e-5], relative=1e-1), locals()

if __name__ == "__main__":
    run_tests()
