# -*- encoding: utf-8 -*-
"""
Supress supercurrent with applied voltage in a 3-probe structure

:Author: Pauli Virtanen <pauli@ltl.tkk.fi>
"""
from __future__ import division
from testutils import *

from numpy import *
from usadel1 import *
from geometries import *
import sys

__revision__ = "$Id: test_optimization_snnns.py 3196 2006-10-04 10:31:09Z pauli $"

#
# Non-trivial test for a three-probe setup: suppress supercurrent
# by applying potential in the normal arm.
#
# First, specify the geometry and solve the spectral equations.  We use
# here `CurrentSolver` that saves its data to an output file. If
# the output file is present, and contains spectral data for this
# geometry, we do not recalculate it.
#

def setup_module():
    global g, sc
    g = geometry_SnnNnS(V0=5, T0=1e-6, phi=1)
    sc = CurrentSolver(g, chunksize=300)

    #
    # COLNEW has problems (at present when writing this test...) solving the
    # kinetic equations. Hence, use a different solver.
    #

    sc.set_solvers(kin_solver=KIN_SOLVER_TWPBVP)
    sc.solve_spectral_if_needed()

def teardown_module():
    global g, sc
    g = None
    sc = None

#
# Specify the optimization problem
#

@slow
def test():
    global sc, g
    def zero_current():
         Ic, IE = sc.get_currents_lazy(w_jT=[0,1], w_jL=[])
         sys.stderr.write('%% %g => %g\n' % (sc.geometry.t_mu[2],
                                             Ic[0] - Ic[1]))
         return [Ic[0] - Ic[1]]
    def adjust_potential(z):
         sc.geometry.t_mu[2] = z[0]

#
# First, solve with one solution of the spectral equations,
#

    z0 = [5]
    z1 = optimize_parameters_for(z0, zero_current, adjust_potential,
                                  xtol=1e-3)

#
# And check the results:
#

    assert len(z1) == 1
    assert is_zero(z1[0] - 8, tolerance=.5)

#
# The expected position of the zero was found from [Heikkila2002]_.
#
#
# .. [Heikkila2002] Heikkilä et al., PRB 66, 184513 (2002)
#

#
# Same test, using G
# ------------------
#
def test_G():
    global sc, g
    def zero_current():
#
# Calculating currents from G is much faster than solving the kinetic
# equations!
#
         Ic, IE = sc.get_currents_from_G(w_jT=[0,1], w_jL=[])
         sys.stderr.write('%% %g => %g\n' % (sc.geometry.t_mu[2],
                                             Ic[0] - Ic[1]))
         return [Ic[0] - Ic[1]]
    def adjust_potential(z):
         sc.geometry.t_mu[2] = z[0]
    z0 = [5]
    z1 = optimize_parameters_for(z0, zero_current, adjust_potential,
                                  xtol=1e-3)
    assert len(z1) == 1
    assert is_zero(z1[0] - 8, tolerance=.5)

if __name__ == "__main__":
    run_tests()
