"""
Evaluate thermopower in the 4-probe 5-wire structure

:Author: Pauli Virtanen <pauli@ltl.tkk.fi>
"""
from __future__ import division
from testutils import *

from numpy import *
from usadel1 import *
from geometries import *
import sys

__revision__ = "$Id: test_optimization_snns_thermopower.py 3196 2006-10-04 10:31:09Z pauli $"

#
# First, specify the geometry and solve the spectral data
#
def setup_module():
    global sc

    T = [2.3360300000, 2.3832200000, 2.3596300000]

    geometry = geometry_SnnNnNnnS(T[0], T[1], T[2], 1)
    E = linspace_weighed(0, 200, 400, ((0, 1, 20),));
    sc = CurrentSolver(geometry, E=E, chunksize=100)
    sc.set_solvers(kin_solver=KIN_SOLVER_BLOCK)
    sc.solve_spectral_if_needed()

def teardown_module():
    global sc
    sc = None

#
# Then, solve the thermopower for the given geometry and temperature
# values
#

def test():
    global sc
    def zero_current():
         Ic, Ie = sc.get_currents_lazy(w_jT=[0,1], w_jL=[])
         return [Ic[0], Ic[1]]
    def adjust_potentials(z):
         sc.geometry.t_mu[0:2] = z[0:2]

    z0 = array([0, 0], float64)
    z = optimize_parameters_for(z0, zero_current, adjust_potentials,
                                 xtol=1e-3, epsfcn=1e-5)
    mu = array(z)

#
# Compare to a result from an old code of mine:
#

    target_mu = array([-0.22640172484e-2, -0.22631056449e-2])
    assert is_zero(max(abs((mu - target_mu)/mu)), tolerance=1e-1)

#
# Same test, with G
# -----------------
#
# Calculating currents from G is much faster:
#
def test_G():
    global sc
    def zero_current():
         Ic, Ie = sc.get_currents_from_G(w_jT=[0,1], w_jL=[])
         return [Ic[0], Ic[1]]
    def adjust_potentials(z):
         sc.geometry.t_mu[0:2] = z[0:2]

    z0 = array([0, 0], float64)
    z = optimize_parameters_for(z0, zero_current, adjust_potentials,
                                 xtol=1e-3, epsfcn=1e-5)
    mu = array(z)

#
# Compare to a result from an old code of mine:
#

    target_mu = array([-0.22640172484e-2, -0.22631056449e-2])
    assert is_zero(max(abs((mu - target_mu)/mu)), tolerance=1e-1)

if __name__ == "__main__":
    run_tests()
