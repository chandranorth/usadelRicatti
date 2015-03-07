"""
:Author:       Pauli Virtanen <pauli@ltl.tkk.fi>
:Organization: Low Temperature Laboratory, Helsinki University of Technology
:Date:         2005-2006

Quasi-Newton solvers (full-rank).

The solvers can be used to solve multi-dimensional nonlinear equations
and fixed-point problems. There are a few classes available:

- `BroydenSolver`: general Broyden solver for ``f(x) = 0``
- `FixedPointSolver`: solver for ``f(x) = x``
  (in these problems scaling and error estimation is easier).
- `DummyFixedPointSolver`: trivial solver for ``f(x) = x``.

The `BroydenSolver` and `FixedPointSolver` have a few choices for
Jacobian updates:

- `BroydenSolver.UPDATE_ICUM`: Inverse column update.
- `BroydenSolver.UPDATE_GOOD_BROYDEN`: "Good" Broyden method.

As usually, the performance of each depends on the problem.
"""

from __future__ import division
import numpy as _n
import util as _util
import warnings
import error

__docformat__ = "restructuredtext en"

class BroydenWarning(error.CoreWarning):
    pass

warnings.simplefilter("ignore", BroydenWarning)

class BroydenSolver(object):
    """Solver a non-linear system of equations, f(x) = 0.

    This solver uses Broyden's Quasi-Newton method.
    
    Usage:
    
        >>> solver = BroydenSolver(x0, initial_scale)
        ... while (solver.relative_change().max() > rel_xtol or
        ...        solver.change().max() > abs_xtol or
        ...        solver.residual().max() > abs_ftol):
        ...     x = solver.get()
        ...     y = f(x)
        ...     solver.add(y)
        ... fixed_point = solver.get()
    """

    UPDATE_ICUM = 1
    UPDATE_GOOD_BROYDEN = 2
    UPDATE_BAD_BROYDEN = 3
    
    def __init__(self, x0, initial_scale, update_type=None,
                 skip_treshold=None):
        """Setup the Broyden non-linear solver.

        :Parameters:
         - `x0`:            Initial guess for the fixed point.
         - `initial_scale`: Initial scaling of the Jacobian.
         - `update_type`: Type of rank-1 Jacobian update to use.
           Possible choices are:
           `BroydenSolver.UPDATE_ICUM` (Inverse column update), and
           `BroydenSolver.UPDATE_GOOD_BROYDEN` ("Good" Broyden).
        """

        if not update_type:
            update_type = BroydenSolver.UPDATE_BAD_BROYDEN

        self.original_shape = _n.shape(x0)
        """Original shape of the iterating vector"""
        self.n = _n.product(self.original_shape)
        """Number of elements in the iterating vector"""

        self.initial_scale = initial_scale
        """Initial scaling of the problem Jacobian"""
        self.update_type = update_type
        """The update type"""

        self.x = _n.reshape(_n.array(x0, copy=True), (self.n, 1))
        """Current point (latest)"""
        self.next_x = self.x
        """Next point"""

        self.H = _n.identity(self.n) / self.initial_scale
        """Approximation to the inverse of the Jacobian (latest)"""
        self.F = _n.zeros((self.n, 1))
        """The target function (latest)"""

        self.iteration = 0
        """Current iteration"""
        self.first_step = True
        """Is the next step the first step?"""

        if skip_treshold is None:
            skip_treshold = 2.0
        self.skip_treshold = skip_treshold
        """Skip steps that change the residual by less than this factor."""

        self.max_step = 1.0
        """Maximum relative step size"""

        self.last_reset = False

    def add(self, next_f):
        """Inform the solver of a new function value."""
        next_f = _n.reshape(_n.array(next_f, copy=True), (self.n, 1))

        # Update the Jacobian, if applicable
        if not self.first_step:
            self.__update_jacobian(self.next_x, next_f)

        # Check if improvement is found
        if (_util.norm2(next_f) < self.skip_treshold*_util.norm2(self.F)
                or self.first_step):
            self.x = self.next_x
            self.F = next_f
            self.last_reset = False
            F = next_f
        else:
            new_H = _n.identity(self.n) / self.initial_scale
            if self.last_reset:
                warnings.warn("%% Broyden reset II", BroydenWarning)
                # also last step failed...
                self.next_x = self.next_x - next_f
                self.F = next_f
            else:
                warnings.warn("%% Broyden reset", BroydenWarning)
                self.next_x = self.x - self.F
            self.last_reset = True
            self.H = new_H
            self.iteration += 1
            return

        if self.first_step:
            self.first_norm = self.residual_norm()
            self.first_step = False


        step = _n.dot(self.H, F)
        factor = min(1., self.max_step*_util.norm2(self.x)/_util.norm2(step))
        step *= factor
        
        self.next_x = self.x - step
        
        self.iteration += 1

    def reset(self):
        """Reset Jacobian."""
        self.H = _n.identity(self.n) / self.initial_scale
        self.F = _n.zeros((self.n, 1))
        self.x = self.next_x
        self.first_step = True

    def __update_jacobian(self, x, F):
        """Update the Jacobian, via a Quasi-Newton method."""
        old_err = _n.seterr(divide='raise')

        try:
            y = F - self.F
            s = x - self.x

            zt = None
            if self.update_type == BroydenSolver.UPDATE_ICUM:
                maxi = abs(_n.ravel(y)).argmax()
                zt = _n.transpose(_n.zeros((1,self.n), _n.float_))
                zt[0, maxi] = 1
            elif self.update_type == BroydenSolver.UPDATE_GOOD_BROYDEN:
                # (Good) Broyden update
                zt = _n.dot(_n.transpose(s), self.H)
            elif self.update_type == BroydenSolver.UPDATE_BAD_BROYDEN:
                # (Bad) Broyden update
                zt = _n.transpose(y)
            else:
                raise ValueError("Unknown update type %s" % (self.update_type))

            self.H = self.H \
                     + _n.dot(s - _n.dot(self.H, y), zt) / _n.dot(zt, y)
        except FloatingPointError:
            warnings.warn("%% Broyden reset: singular", BroydenWarning)
            self.H = _n.identity(self.n) / self.initial_scale

        _n.seterr(**old_err)

    def get(self):
        """Get a point where to evaluate the next function value."""
        return _n.reshape(self.next_x, self.original_shape)

    def residual(self):
        """Return the residual."""
        return _n.reshape(self.F, self.original_shape)

    def residual_norm(self):
        return abs(self.residual()).max()

    def relative_residual_norm(self):
        if self.first_step:
            return _n.inf
        return self.residual_norm() / self.first_norm

    def change(self):
        """Return the change in x."""
        return _n.reshape(self.next_x - self.x, self.original_shape)
        

class FixedPointSolver(BroydenSolver):
    """Solver for a fixed point problem of the form  x = f(x).

    This solver uses Broyden's Quasi-Newton method to accelerate the
    fixed-point iteration.
    
    Usage:
        >>> solver = FixedPointSolver(x0)
        ... while solver.relative_residual.max() > tolerance:
        ...     x = solver.get()
        ...     y = f(x)
        ...     solver.add(y)
        ... fixed_point = solver.get()
    """
    
    def __init__(self, x0, initial_scale=None, update_type=None,
                 skip_treshold=None):
        """Setup the Broyden accelerated fixed-point iteration.

        :Parameters:
         - `x0`:            Initial guess for the fixed point.
         - `initial_scale`: Initial scaling of the Jacobian.
         - `update_type`:   Type of rank-1 Jacobian update to use.
        """
        if not initial_scale:
            initial_scale = 1
        BroydenSolver.__init__(self, x0, initial_scale, update_type,
                               skip_treshold)

    def add(self, next_f): 
        """Inform the solver of a new function value."""
        next_f = _n.reshape(next_f, (self.n, 1))
        BroydenSolver.add(self, self.next_x - next_f)

class DummyFixedPointSolver(object):
    """Simple fixed point iteration."""

    def __init__(self, x0):
        self.x = _n.array(x0, copy=True)
        self.last_x = self.x
        self.iteration = 0

    def get(self):
        return self.x

    def add(self, next_f):
        self.last_x = self.x
        self.x = _n.array(next_f, copy=True)
        if self.iteration == 0:
            self.first_norm = self.residual_norm()
        self.iteration += 1

    def reset(self):
        pass

    def residual(self):
        return self.x - self.last_x

    def residual_norm(self):
        return abs(self.residual()).max()

    def relative_residual_norm(self):
        if self.iteration == 0:
            return _n.inf
        return self.residual_norm() / self.first_norm
