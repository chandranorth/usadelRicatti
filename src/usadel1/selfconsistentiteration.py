"""
:Author:       Pauli Virtanen <pauli@ltl.tkk.fi>
:Organization: Low Temperature Laboratory, Helsinki University of Technology
:Date:         2005-2006

Self-consistent iteration.
"""

from __future__ import division

import numpy as _n
import sys as _sys
import re as _re
import os as _os
import scipy, scipy.integrate, scipy.interpolate, scipy.optimize
import tables as _tbl
import shutil as _shutil
import datetime as _datetime

from nonlinearsolver import *
from util import *
import data
from solver import *
import error

__all__ = ['self_consistent_realtime_iteration',
           'self_consistent_matsubara_iteration',
           'dump_delta', 'bulk_delta']

__docformat__ = "restructuredtext en"

    
#------------------------------------------------------------------------------
# Real-time self-consistent iteration
#------------------------------------------------------------------------------

def _set_self_consistent_delta(g, sp, kin, coef, output_func=lambda x: x):
    """Evaluate Delta self-consistently for given data (real-time)"""
    
    nx = len(sp.x)
    ne = len(sp.E)

    # Evaluate integrand

    for w in range(g.nwire):
        if not (g.coupling_lambda[w,:].max() > 0 and g.w_type[w]&WIRE_TYPE_S):
            continue

        ### Calculate the contribution in the solved energy interval
            
        # Restrict integral below omega_D for each wire
        E = stretch_shape(sp.E, (-1, nx) )[:,:]
        k = (E > g.omega_D[w,:])
        E[k] = g.omega_D[w, k[:,1]]

        a = sp.a[:,w,:]
        b = sp.b[:,w,:]

        ab_1 = a*b - 1

        fL = scipy.interpolate.interp1d(kin.E, kin.fL[:,w,:], axis=0,
                                        fill_value=0, bounds_error=False)
        fT = scipy.interpolate.interp1d(kin.E, kin.fT[:,w,:], axis=0,
                                        fill_value=0, bounds_error=False)

        # Evaluate
        f_K = 2*(  (-a/ab_1 - (b/ab_1).conj()) * fL(sp.E)
                 + (-a/ab_1 + (b/ab_1).conj()) * fT(sp.E))

        # Integrate
        z = scipy.integrate.trapz(f_K, E, axis=0)

        #-- Cutoff elimination

        delta_0 = g.omega_D[w,:] / _n.cosh(1./g.coupling_lambda[w,:])
        z0 = 2 * delta_0 * _n.arccosh(E[-1] / delta_0)

        def darccosh(z):
            q = _n.arccosh(g.omega_D[w]/z) - _n.arccosh(E[-1]/z)
            q[~_n.isfinite(q)] = 0
            return q
        # NB: this remainder term is not very significant
        remainder = delta_0 * g.w_delta[w] * (
            darccosh(g.w_delta[w]) - darccosh(delta_0))
        remainder = remainder * _n.exp(1j*_n.angle(z))
        
        z = 1/z0 * (delta_0 * z + remainder)
        z[~_n.isfinite(z)] = 0

        #-- Go back to polar representation
        g.w_delta[w,:] = abs(z)
        g.w_phase[w,:] = continuous_phase(_n.angle(z), center=True)


def get_current_conservation_violation(g, kin, coef):
    """Evaluates the amount by which current conservation is violated.

    The returned value is the maximum of standard deviations of current
    in the wires. (Exact current conservation at nodes is assumed.)
    """

    nx = len(kin.x)

    stds = []

    for w in range(g.nwire):
        jT = kin.jT[:,w,:]
        jL = kin.jL[:,w,:]

        E = kin.E[:,None]
        Ic = scipy.integrate.trapz(jT, E, axis=0) * g.w_conductance[w]
        IE = scipy.integrate.trapz(jL*E, E, axis=0) * g.w_conductance[w]

        def stdavg(x, y):
            avg = scipy.integrate.trapz(y[:], x[:])
            std = _n.sqrt(scipy.integrate.trapz((y[:] - avg)**2, x[:]))
            return std, avg

        stdIc, avgIc = stdavg(kin.x, Ic)
        stdIE, avgIE = stdavg(kin.x, IE)

        stds.append((stdIc + stdIE) / (abs(avgIc) + abs(avgIE) + 1e-6))

    return max(stds)

def self_consistent_realtime_iteration(currents,
                                       max_iterations=100,
                                       output_func=_sys.stderr.write,
                                       iterator=FixedPointSolver):
    """
    Self-consistent real-time iteration for the order parameter.

    Adjusts the parameters in the given `Geometry` object until convergence
    is met.

    .. note::

       Always check that your results are independent of the energy cutoffs,
       by adjusting the energy range and number of points in the
       :class:`CurrentSolver` you use.

    Parameters
    ----------
    solver : CurrentSolver
        A solver containing the geometry in which to solve for :math:`\\Delta`
    max_iterations : int, optional
        Maximum number of iterations to make
    output_func : func(message), optional
        Function to use for printing messages. If omitted, prints to stderr.
    iterator
        Fixed-point iteration accellerator. Suitable choices are
        `FixedPointSolver` (Broyden solver, default) and `DummyFixedPointSolver`
        (simple fixed-point iteration).

    Yields
    ------
    k : int
        Number of current iteration
    d : fixed-point solver object
        Has the methods `residual_norm()` that returns a max-norm estimating
        the error in :math:`\\Delta`, and `relative_residual_norm()` that returns
        the norm scaled by the initial norm.
    v : float
        Amount by which current conservation is violated.

    Notes
    -----
    In real-time iteration, some care must be taken to choose the energy grid
    so that it is sufficiently dense in the energy range where the peak in
    the :math:`F^K`-function will be.

    The energy cutoff is eliminated by using the formula:

    .. math::
       \\Delta = [\\int_0^{\\epsilon_0}d\\epsilon F^K_0(\\epsilon)]^{-1} \\left(
       \\Delta_0 \\int_0^{\\epsilon_0}d\\epsilon F^K(\\epsilon)
       + \\int_{\\epsilon_0}^{\\omega_c}d\\epsilon
       [\\Delta_0 F^K(\\epsilon) - \\Delta F^K_0(\\epsilon)]
       \\right)

    and in the latter term taking the bulk expressions for :math:`F^K`
    and :math:`F^K_0`, the latter being the bulk function corresponding to
    :math:`\\Delta_0` (the bulk gap). Here, :math:`\\epsilon_0` is the energy
    limit up to which numerical solutions are calculated for :math:`F^K`.

    Examples
    --------
    >>> solver = CurrentSolver(geometry)
    >>> it = self_consistent_realtime_iteration(solver)
    >>> for k, d, v in it:
    >>>     if d.residual_norm() < 0.1 and v < 1e-4:
    >>>         break
    >>> else:
    >>>     raise RuntimeError("Iteration did not converge")

    """

    geometry = currents.geometry

    delta_ex = iterator(geometry.w_delta * _n.exp(1j*geometry.w_phase))

    for k in range(max_iterations):
        zdelta = delta_ex.get()
        geometry.w_delta = abs(zdelta)
        geometry.w_phase = _n.angle(zdelta)

        if output_func is not None:
            dump_delta(geometry, output_func)
        
        currents.solve()

        _set_self_consistent_delta(geometry, currents.spectral,
                                   currents.kinetic, currents.coefficient,
                                   output_func=output_func)

        delta = geometry.w_delta[:,:] * _n.exp(1j*geometry.w_phase[:,:])
        delta_ex.add(delta)

        violation = get_current_conservation_violation(currents.geometry,
                                                       currents.kinetic,
                                                       currents.coefficient)
        yield k, delta_ex, violation

def dump_delta(geometry, _output=_sys.stderr.write):
    """Dump the Delta in all S-wires to output."""

    delta = geometry.w_delta
    phase = geometry.w_phase

    jstep = 1 + delta.shape[1]//5
    j0 = (delta.shape[1] % jstep) // 2

    def fmt_x_array(z):
        return " ".join(["%10.2g" % x for x in z[j0::jstep]])

    for i in range(geometry.nwire):
        if geometry.w_type[i] & WIRE_TYPE_S:
            _output("DELTA[%2d]    = [%s]\n" % (i, fmt_x_array(delta[i])))
    for i in range(geometry.nwire):
        if geometry.w_type[i] & WIRE_TYPE_S:
            _output("PHASE[%2d]/PI = [%s]\n" % (i, fmt_x_array(phase[i]/_n.pi)))

#------------------------------------------------------------------------------
# Matsubara self-consistent iteration
#------------------------------------------------------------------------------

class MatsubaraWarning(error.CoreWarning):
    pass

def self_consistent_matsubara_iteration(geometry,
                                        max_iterations=100,
                                        max_ne=300,
                                        output_func=_sys.stderr.write,
                                        iterator=FixedPointSolver,
                                        E_max=None,
                                        force_integral=False):
    """
    Self-consistent Matsubara iteration for the order parameter.

    Adjusts the parameters in the given `Geometry` object until convergence
    is met. Applicable only to equilibrium situations.

    .. note::

       Always check that your results are independent of the energy cutoffs,
       by adjusting the *max_ne* and *E_max* parameters.

    Parameters
    ----------
    geometry : Geometry
        The geometry in which to solve for :math:`\\Delta`
    max_iterations : int, optional
        Maximum number of iterations to make
    max_ne : int, optional
        Number of Matsubara frequencies to sum over.

        .. note:: If the number of frequencies given is smaller than that
           required to sum to frequencies over :math:`\\Delta`, an energy integral
           is performed instead of a discrete sum.

    output_func : func(message), optional
        Function to use for printing messages. If omitted, prints to stderr.
    iterator
        Fixed-point iteration accellerator. Suitable choices are
        `FixedPointSolver` (Broyden solver, default) and `DummyFixedPointSolver`
        (simple fixed-point iteration).
    E_max : float
        Cutoff energy to use (smaller than Debye, but larger than energy
        scales of the structure).

        If ``None``, determined automatically. The automatic cutoff satisfies:
        - It is no larger than ``omega_D``.
        - The corresponding length scale is smaller than 0.05 of shortest wire
          or 0.1 of coherence length

    force_integral : bool
        Force integration, even if `max_ne` is large enough for summing
        up to E_max.

    Yields
    ------
    k : int
        Number of current iteration
    d : fixed-point solver object
        Has the methods `residual_norm()` that returns a max-norm estimating
        the error in :math:`\\Delta`, and `relative_residual_norm()` that returns
        the norm scaled by the initial norm.
    v : float
        Amount by which current conservation is violated.
        (NB: not currently implemented in Matsubara iteration, always zero.)

    Notes
    -----
    The energy cutoff is eliminated by using the formula:

    .. math::
       \\Delta = [\\sum_{\\omega_k < \\epsilon_0} F_0(\\omega_k)]^{-1} \\left(
       \\Delta_0 \\sum_{\\omega_k < \\epsilon_0} F(\\omega_k)
       + \\sum_{\\epsilon_0 < \\omega_k < \\omega_c} [\\Delta_0 F(\\omega_k) - \\Delta F_0(\\omega_k)]
       \\right)

    and in the latter term taking the bulk expressions for :math:`F`
    and :math:`F_0`, the latter being the bulk function corresponding to
    :math:`\\Delta_0` (the bulk gap). Here, :math:`\\epsilon_0` is the energy
    limit up to which numerical solutions are calculated for `F`.

    Examples
    --------
    >>> it = self_consistent_matsubara_iteration(geometry)
    >>> for k, d, v in it:
    >>>     if d.residual_norm() < 0.1 and v < 1e-4:
    >>>         break
    >>> else:
    >>>     raise RuntimeError("Iteration did not converge")

    """

    if not ((_n.diff(geometry.t_t) == 0).all() and (geometry.t_mu == 0).all()):
        raise ValueError("The geometry does not describe an equilibrium "
                         "situation.")

    T = geometry.t_t[0]

    # XXX: the Matsubara iteration still needs work

    delta_ex = iterator(geometry.w_delta * _n.exp(1j*geometry.w_phase))

    ne = 0

    _E_max = E_max

    for k in range(max_iterations):
        zdelta = delta_ex.get()
        geometry.w_delta = abs(zdelta)
        geometry.w_phase = _n.angle(zdelta)

        if output_func is not None:
            dump_delta(geometry, output_func)

        # Pick an energy cutoff:
        #
        # Criteria are explained in the docstring, except for the technical
        # one:
        #
        # - it should not vary much between different iterations
        #
        if _E_max is None:
            delta_0 = bulk_delta(geometry.coupling_lambda,
                                 geometry.omega_D, T)
            E_T = 1/geometry.w_length**2
            E_T[geometry.w_type != WIRE_TYPE_S] = 0
            E_max = _n.minimum(
                _n.maximum(2500*E_T[:,None],
                           100*_n.clip(geometry.w_delta,
                                       geometry.t_delta.max(), _n.inf)),
                geometry.omega_D
                ).max()
        else:
            E_max = _E_max
        delta_ne = int(round(E_max/(2*_n.pi*T))) + 10

        if delta_ne > ne:
            ne = delta_ne
        elif delta_ne < 3*ne/4:
            ne = delta_ne

        if ne > max_ne or force_integral:
            # Calculate an integral instead:
            # - sum the low-energy part more densely
            # - and the high-energy part less densely, with energy grid
            #   E_n ~ n^2
            #
            ne_x = max_ne//4

            nq = ne_x + (ne-ne_x) * linspace(0, 1, max_ne-ne_x)**2

            # Midpoint scheme
            weights = _n.r_[_n.ones(ne_x), _n.diff(nq)]
            n = _n.r_[_n.arange(ne_x), .5*(nq[1:] + nq[:-1])]
            E = 2j*_n.pi*T*(n + .5)

            # Check energy discretization
            def check_worry(Ex):
                ddE = _n.diff(abs(E)[ne_x:]).min()
                return ((Ex[...,None] > abs(E)[ne_x])
                        & (Ex[...,None] < 10*ddE)).any()
            
            worry_flag = check_worry(geometry.w_delta)
            if worry_flag:
                warnings.warn("The structure has energy scales "
                              "in energy range that may be too sparsely "
                              "sampled. Consider checking that your results "
                              "do not depend on cutoff parameters `max_ne` "
                              "and `E_max`. Currently, max_ne=%d"
                              % (max_ne,),
                              MatsubaraWarning)
        else:
            E = 2j*_n.pi*T*(_n.arange(ne) + .5)
            weights = _n.ones([ne])

        # cut away too high energies
        j = (abs(E) <= geometry.omega_D.max())
        E = E[j]
        weights = weights[j]
        ne = len(E)

        # solve
        solver = Solver()
        solver.set_geometry(geometry)
        sol = solver.sp_solve(E, geometry.x)

        # compute a new Delta
        _set_self_consistent_delta_m(geometry, sol, T, weights,
                                     output_func=output_func)

        delta_ex.add(geometry.w_delta[:,:] * _n.exp(1j*geometry.w_phase[:,:]))

        #violation = get_current_conservation_violation(geometry,
        #                                               kinetic,
        #                                               coefficient)

        violation = 0
        
        yield k, delta_ex, violation

def _set_self_consistent_delta_m(g, sp, T, weights, output_func=lambda x: x):

    nx = len(sp.x)
    ne = len(sp.E)

    for w in range(g.nwire):
        if not (g.coupling_lambda[w,:].max() > 0 and g.w_type[w]&WIRE_TYPE_S):
            continue
        
        ### Calculate the contribution in the solved energy interval

        E = stretch_shape(sp.E, (-1, nx) )[:,:]
        
        a = sp.a[:,w,:]
        b = sp.b[:,w,:]
        ab_1 = a*b - 1
        
        # Evaluate
        f_K_sp = 4j*(-a/ab_1) * weights[:,None]
        f_K_sp[abs(E) > g.omega_D[w]] = 0

        # Sum
        z = f_K_sp.sum(axis=0)

        ### Cutoff elimination
        delta_0 = bulk_delta(g.coupling_lambda[w], g.omega_D[w], T)
        f_K_sp0_per_delta0 = 2/_n.sqrt(abs(delta_0)**2 + abs(E)**2) * weights[:,None]
        f_K_sp0_per_delta0[abs(E) > g.omega_D[w]] = 0
        z0_per_delta0 = f_K_sp0_per_delta0.sum(axis=0)


        # Remainder
        def darccosh(z):
            q = _n.arccosh(g.omega_D[w]/z) - _n.arccosh(abs(E[-1])/z)
            q[~_n.isfinite(q)] = 0
            q[abs(E[-1]) > g.omega_D[w]] = 0
            return q
        # NB: this remainder term is not very significant
        remainder = delta_0 * g.w_delta[w] * (
            darccosh(g.w_delta[w]) - darccosh(delta_0))
        remainder = remainder * _n.exp(1j*_n.angle(z))

        remainder /= z0_per_delta0*delta_0

        z = z/z0_per_delta0 + remainder

        ### Cannot use cutoff elimination where bulk gap vanishes, fall back
        ### to summation
        m = (delta_0 == 0)
        z[m] = g.coupling_lambda[w][m] * f_K_sp.sum(axis=0)[m]

        ### Go back to polar representation
        g.w_delta[w,:] = abs(z)
        g.w_phase[w,:] = continuous_phase(_n.angle(z), center=True)


_bulk_delta_cache = {}

@_n.vectorize
def bulk_delta(lambda_0, omega_D, T):
    """
    Evaluate :math:`Delta` in a bulk superconductor.

    Parameters
    ----------
    lambda_0
        Coupling parameter
    omega_D
        Debye cutoff frequency
    T
        Temperature

    """
    key = (float(lambda_0), float(omega_D), float(T))
    if key in _bulk_delta_cache:
        return _bulk_delta_cache[key]

    if lambda_0 == 0:
        return 0.0

    Delta_0 = 2*omega_D*_n.exp(-1/lambda_0)

    ne = int(omega_D/(2*_n.pi*T))
    if ne < 1000:
        def delta_fun(delta):
            delta = abs(delta)
            w = 2*_n.pi*T*(_n.arange(0, ne+1) + .5)
            ii = (2*_n.pi*T/_n.sqrt(w**2 + delta**2)).sum()
            res = 1 - lambda_0*ii
            return res
    else:
        def delta_fun(delta):
            delta = abs(delta)
            ii, ier = scipy.integrate.quad(
                lambda ee: (1/_n.sqrt((ee+1e-9j)**2 - delta**2)).real * _n.tanh(ee/(2*T)),
                delta, omega_D)
            res = 1 - lambda_0*ii
            return res

    if delta_fun(Delta_0*2)*delta_fun(0) > 0:
        # No superconducting state possible
        result = 0.0 # note: must return float so that vectorize detects otype properly
    else:
        result = scipy.optimize.brenth(delta_fun, 0, Delta_0*2)

    _bulk_delta_cache[key] = result
    return result
