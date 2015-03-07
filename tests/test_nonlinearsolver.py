"""
Test nonlinear solver

:Author: Pauli Virtanen <pauli@ltl.tkk.fi>
"""
from __future__ import division
from testutils import *

import numpy as _n
import usadel1.nonlinearsolver as nls
 
__revision__ = "$Id: test_nonlinearsolver.py 3184 2006-10-02 07:05:45Z pauli $"

def test():
    def test_iteration(solver_cls, x0, func, maxiter):
         solver = solver_cls(x0)
         for i in range(maxiter):
             solver.add(func(solver.get()))
             if solver.residual_norm() < 1e-12:
                 break
         return (solver.get(), solver)
# 
# Simple problem
# --------------
# 
# Check that a simple 1D fixed-point iteration converges reasonably fast.
# 
    func = lambda x: _n.sqrt(1 - x)
    x_exact = .5*(_n.sqrt(5.0) - 1)

    x, s = test_iteration(nls.FixedPointSolver, 0.3, func, 6)
    assert is_zero(x - x_exact, tolerance=1e-10)
    
# 
# For this problem, Broyden's Quasi-Newton method is very nice.  A brute
# force fixed-point iteration takes much more work:
#

    x, s = test_iteration(nls.DummyFixedPointSolver, 0.3, func, 150)
    assert is_zero(x - x_exact, tolerance=1e-10)

# 
# Then, solve the same problem, but using arrays -- to check whether
# dimensions work as intended:
#

    shape = (50, 4, 2)
    x, s = test_iteration(nls.FixedPointSolver, _n.zeros(shape)+0.3, func, 9)
    assert is_zero(array_normmax(x - x_exact), tolerance=1e-10)

# 
# Test also now the auxiliary functions
#
    assert s.residual().shape == shape
    assert arrays_equal(s.residual(), 0, tolerance=1e-9)

if __name__ == "__main__":
    run_tests()
