"""
Test a trivial optimization problem in an NnN structure

:Author: Pauli Virtanen <pauli@ltl.tkk.fi>
"""

from __future__ import division
from testutils import *

from geometries import *

from usadel1 import *
from numpy import *

__revision__ = "$Id: test_optimization_nnn.py 3184 2006-10-02 07:05:45Z pauli $"

# 
# A trivial test: optimize charge current to zero between two normal terminals.
# 
def test():
 
    g = geometry_NnN(V0=15, T0=1e-6, phi=1)
    sc = CurrentSolver(g, output_function=lambda s: None)
 
    def zero_current():
         Ic, Ie = sc.get_currents_lazy(w_jT=[0], w_jL=[])
         return [Ic[0]]
        
    def adjust_potential(z):
         sc.geometry.t_mu[0] = z[0]-1
# 
# We put here z = V+1
# 
    z0 = [15+1]
    sc.solve_spectral()
    z = optimize_parameters_for(z0, zero_current, adjust_potential,
                                 xtol=1e-3)
# 
# Check that the result is sensible
# 
    assert len(z) == 1
    assert is_zero(abs(z[0]-1), tolerance=1e-3)

if __name__ == "__main__":
    run_tests()
