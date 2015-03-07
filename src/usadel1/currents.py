"""
:Author:       Pauli Virtanen <pauli@ltl.tkk.fi>
:Organization: Low Temperature Laboratory, Helsinki University of Technology
:Date:         2005-2006

Evaluating currents.

Basic usage
===========

The simplest interface:

  >>> geometry = get_some_geometry()
  >>> solver = currents.CurrentSolver.resume('data.h5', geometry)
  >>> solver.solve_spectral()
  >>> solver.save('data.h5')
  >>> Ic, IE = solver.get_currents_lazy()

This picks a sensible energy and spatial discretizations, and
stores/resumes the run from the specified file.


Optimization
------------

Another thing to do is optimizing parameter values, for example to get
the thermopower. For example:

  >>> def zero_currents()
  ...     Ic, IE = solver.get_currents_lazy(w_jL=[], w_jT=[0,1])
  ...     return Ic[0:2]
  
  >>> def adjust_potentials(z)
  ...     solver.geometry.t_mu[0:2] = z

  >>> z = currents.optimize_parameters_for(x0=[0,0],
  ...                                      zero_currents,
  ...                                      adjust_potentials,
  ...                                      xtol=1e-3, epsfcn=1e-5)
  >>> solver.geometry.t_mu[0:2] = z

This adjusts potentials in terminals 0, 1 so that no current flows in
wires 0, 1.
"""

import numpy as _n
import scipy
import scipy.integrate
import scipy.optimize
import scipy.linalg

import sys, os, shutil, traceback
from util import *
import nonlinearsolver as _nl
from solver import *
import data, tables
from error import *

import tables

__all__ = ['calculate_currents', 'optimize_parameters_for',
           'calculate_currents_lazy',
           'calculate_currents_from_G_with_f',
           'calculate_linear_response_from_G_with_f',
           'approximate_G',
           'CurrentSolver',
           'twolevel_Egrid',
           'AlreadyConvergedException']
__docformat__ = "restructuredtext en"

class AlreadyConvergedException(Exception):
    """An exception that the user can raise during `optimize_parameters_for`\
    to indicate that the optimum has been already reached."""
    pass


def calculate_currents_lazy(geometry, coefs, solver,
                            w_jL=None, w_jT=None, epsabs=None, epsrel=None):
    """Calculate currents flowing in wires, using adaptive integration \
    in energy.

    The current is evaluated at x-index 0.

    :Parameters:
      - `geometry`: the geometry to use
      - `solver`: `solver.Solver` object to use
      - `coefs`: `solver.KineticCoefficients` to use
      - `w_jL`: a list of wire indices for which to calculate energy current
        ``None`` indicates all wires.
      - `w_jT`: a list of wire indices for which to calculate charge current
        ``None`` indicates all wires
      - `epsabs`: absolute error tolerance
      - `epsrel`: relative error tolerance
    """

    if w_jT is None:
        w_jT = range(geometry.nwire)
    if w_jL is None:
        w_jL = range(geometry.nwire)

    def kin_solve(E):
        E = _n.array([E], _n.float_)
        return solver.kin_solve(E, geometry.x)
        
    def integrand_jT(E, wire):
        rk = kin_solve(E)
        return rk.jT[0,wire,0]

    def integrand_jL(E, wire):
        rk = kin_solve(E)
        return E*rk.jL[0,wire,0]

    Ic = _n.zeros([geometry.nwire]) + _n.nan
    IE = _n.zeros([geometry.nwire]) + _n.nan

    def do_quad(f, a0, b0, *args, **kwargs):
        breakpoints = scipy.unique(
            [a0] + list(geometry.t_mu) + list(2*geometry.t_t)
            + list(geometry.t_delta) + [b0])
        breakpoints = filter(lambda x: b0 >= x and x >= a0 and x < scipy.inf,
                             breakpoints)
        r = scipy.integrate.quad(f, a0, b0, full_output=True,
                                 points=breakpoints, *args, **kwargs)
        if len(r) == 3:
            z, err, info = r
        elif len(r) == 2:
            z, err = r
        elif len(r) == 4:
            z, err, info, warning = r
            warnings.warn(str((warning, info)), NoConvergenceWarning)
        else:
            raise NoConvergenceException("Integration failed", r)
        return z, err

    maxE = 2*max(coefs.E)

    if epsabs is None: epsabs = 1e-4
    if epsrel is None: epsrel = 1e-2

    for w in w_jT:
        Ic[w], abserr = do_quad(
            integrand_jT, 0, maxE, args=(w,), epsabs=epsabs, epsrel=epsrel,
            limit=100)
        Ic[w] *= geometry.w_conductance[w]
    for w in w_jL:
        IE[w], abserr = do_quad(
            integrand_jL, 0, maxE, args=(w,), epsabs=epsabs, epsrel=epsrel,
            limit=100)
        IE[w] *= geometry.w_conductance[w]

    return (Ic, IE)
        

def calculate_currents(geometry, kinetic, ix=None,
                       w_jL=None, w_jT=None):
    """Calculate the charge and energy currents corresponding to the \
    given data.

    :Parameters:
      - `geometry`: The geometry to use
      - `kinetic`: The kinetic data to use
      - `ix`: The x-indices where to evaluate
      - `w_jL`: a list of wire indices for which to calculate energy current
        ``None`` indicates all wires.
      - `w_jT`: a list of wire indices for which to calculate charge current
        ``None`` indicates all wires
    """

    if not ix is None:
        kinetic = kinetic[:,:,slice(ix,ix+1)]

    nx = len(kinetic.x)
    ne = len(kinetic.E)
    nw = geometry.nwire

    if w_jT is None:
        w_jT = range(nw)
    if w_jL is None:
        w_jL = range(nw)

    Ic = _n.zeros([nw, nx], _n.float_) + _n.nan
    IE = _n.zeros([nw, nx], _n.float_) + _n.nan

    E = kinetic.E[:,_n.newaxis]

    for i, w in enumerate(w_jT):
        Ic[w,:] = geometry.w_conductance[w] * scipy.integrate.trapz(
            kinetic.jT[:,w,:], E, axis=0)

    for i, w in enumerate(w_jL):
        IE[w,:] = geometry.w_conductance[w] * scipy.integrate.trapz(
            kinetic.jL[:,w,:]*E, E, axis=0)
    
    if not ix is None:
        Ic = Ic[:,0]
        IE = IE[:,0]

    return (Ic, IE)


def _find_all_terminals(geometry):
    terminals = [NODE_CLEAN_N_TERMINAL, NODE_CLEAN_S_TERMINAL,
                 NODE_TUNNEL_N_TERMINAL, NODE_TUNNEL_S_TERMINAL]
    return [t for t in range(geometry.nnode)
            if (geometry.t_type[t] & BCTYPE_MASK) in terminals]

def calculate_currents_from_G_with_f(geometry, E, G, f_func,
                                     w_jT=None, w_jL=None):
    g = geometry

    ne = len(E)
    nw = g.nwire

    terminals = _find_all_terminals(geometry)

    if w_jT is None:
        w_jT = range(nw)
    if w_jL is None:
        w_jL = range(nw)

    Ic = _n.empty([nw], _n.float_) + _n.nan
    IE = _n.empty([nw], _n.float_) + _n.nan

    jL = _n.zeros([ne, nw], _n.float_)
    jT = _n.zeros([ne, nw], _n.float_)

    fL, fT = f_func(E[:,None], g.t_mu[None,:], g.t_t[None,:])

    for t in terminals:
        jL += G[:,:,0,0,t] * fL[:,None,t]
        jL += G[:,:,0,1,t] * fT[:,None,t]
        jT += G[:,:,1,0,t] * fL[:,None,t]
        jT += G[:,:,1,1,t] * fT[:,None,t]

    for i, w in enumerate(w_jT):
        Ic[w] = g.w_conductance[w] * scipy.integrate.trapz(
            jT[:,w], E, axis=0)


    for i, w in enumerate(w_jL):
        IE[w] = g.w_conductance[w] * scipy.integrate.trapz(
            jL[:,w]*E, E, axis=0)

    return (Ic, IE)

def calculate_linear_response_from_G_with_f(geometry, E, G, df_func,
                                            w_jT=None, w_jL=None):
    g = geometry

    ne = len(E)
    nw = g.nwire

    terminals = _find_all_terminals(geometry)

    if w_jT is None:
        w_jT = range(nw)
    if w_jL is None:
        w_jL = range(nw)

    G_VV = _n.asmatrix(_n.empty([nw, g.nnode], _n.float_)) + _n.nan
    G_VT = _n.asmatrix(_n.empty([nw, g.nnode], _n.float_)) + _n.nan
    G_TV = _n.asmatrix(_n.empty([nw, g.nnode], _n.float_)) + _n.nan
    G_TT = _n.asmatrix(_n.empty([nw, g.nnode], _n.float_)) + _n.nan

    dfL_dT, dfL_dV, dfT_dT, dfT_dV = df_func(E[:,None],
                                             g.t_mu[None,:],
                                             g.t_t[None,:])

    for t in terminals:
        # Temperature response
        jL  = G[:,:,0,0,t] * dfL_dT[:,None,t]
        jL += G[:,:,0,1,t] * dfT_dT[:,None,t]
        jT  = G[:,:,1,0,t] * dfL_dT[:,None,t]
        jT += G[:,:,1,1,t] * dfT_dT[:,None,t]

        for w in w_jT:
            G_VT[w,t] = g.w_conductance[w] * scipy.integrate.trapz(
                jT[:,w], E, axis=0)
        for w in w_jL:
            G_TT[w,t] = g.w_conductance[w] * scipy.integrate.trapz(
                jL[:,w]*E, E, axis=0)

        # Potential response
        jL  = G[:,:,0,0,t] * dfL_dV[:,None,t]
        jL += G[:,:,0,1,t] * dfT_dV[:,None,t]
        jT  = G[:,:,1,0,t] * dfL_dV[:,None,t]
        jT += G[:,:,1,1,t] * dfT_dV[:,None,t]

        for w in w_jT:
            G_VV[w,t] = g.w_conductance[w] * scipy.integrate.trapz(
                jT[:,w], E, axis=0)
        for w in w_jL:
            G_TV[w,t] = g.w_conductance[w] * scipy.integrate.trapz(
                jL[:,w]*E, E, axis=0)

    return G_VV, G_VT, G_TV, G_TT


def approximate_G(geometry, coefficient, use_jS=True, use_T=True):
    """
    Evaluate an approximation for the G-matrix.

    It is solved from equations valid in the first order in j_S and T.
    """
    trapz = scipy.integrate.trapz
    
    g = geometry
    c = coefficient
    ne = len(c.E)

    # Evaluate some coefficients

    R = g.w_length / g.w_conductance

    MT = trapz(1/c.DT, c.x[None,None,:], axis=2) * R[None,:]
    ML = trapz(1/c.DL, c.x[None,None,:], axis=2) * R[None,:]

    if use_T:
        TT = trapz(c.TT/c.DL/c.DT, c.x[None,None,:], axis=2) * R[None,:]
        TT /= MT*ML

    if use_jS:
        sgn = _n.sign(c.x[:, None] - c.x[None, :])[None,None,:,:]
        gamma = trapz(
            1/c.DT * trapz(sgn/c.DL[:,:,None,:], c.x[None,None,None,:],axis=3),
            c.x[None,None,:],
            axis=2) * R[None,:]*R[None,:]
        gamma /= MT*ML

        jS = c.ijE[:,:,0] * g.w_conductance[None,:]

    # Formulate node-to-node conductance matrix,  I_i = G_ij f_j

    L = _n.zeros([ne, 2, g.nwire, 2, g.nnode])

    node_wires = {}
    
    for i, (w_type, w_ends) in enumerate(zip(g.w_type, g.w_ends)):
        
        if w_type != WIRE_TYPE_N:
            raise RuntimeError("Only N-wires can be handled!")

        # fix clean contacts to superconductors

        if g.t_type[w_ends[0]] == NODE_CLEAN_S_TERMINAL:
            subgap = (c.E < g.t_delta[w_ends[0]])
            ML[subgap,i]    = _n.inf
            if use_jS:
                gamma[subgap,i] = 1
            if use_T:
                TT[subgap,i]    = 0
        elif g.t_type[w_ends[1]] == NODE_CLEAN_S_TERMINAL:
            subgap = (c.E < g.t_delta[w_ends[1]])
            ML[subgap,i]    = _n.inf
            if use_jS:
                gamma[subgap,i] = -1
            if use_T:
                TT[subgap,i]    = 0

        # approximation up to first order in j_S, {\cal T}
        
        L[:, 0, i, 0, w_ends[1]] +=  1/ML[:, i]
        L[:, 0, i, 0, w_ends[0]] += -1/ML[:, i]

        L[:, 1, i, 1, w_ends[1]] +=  1/MT[:, i]
        L[:, 1, i, 1, w_ends[0]] += -1/MT[:, i]

        if use_T:
            L[:, 0, i, 1, w_ends[1]] += -TT[:,i]
            L[:, 0, i, 1, w_ends[0]] += +TT[:,i]

            L[:, 1, i, 0, w_ends[1]] += +TT[:,i]
            L[:, 1, i, 0, w_ends[0]] += -TT[:,i]

        if use_jS:
            L[:, 0, i, 1, w_ends[1]] += jS[:,i] * (1 - gamma[:,i])/2
            L[:, 0, i, 1, w_ends[0]] += jS[:,i] * (1 + gamma[:,i])/2

            L[:, 1, i, 0, w_ends[1]] += jS[:,i] * (1 + gamma[:,i])/2
            L[:, 1, i, 0, w_ends[0]] += jS[:,i] * (1 - gamma[:,i])/2

        node_wires.setdefault(w_ends[1], {})[i] = True
        node_wires.setdefault(w_ends[0], {})[i] = True

    # Condense inner nodes out (separately at each energy)
    
    nodes     = (geometry.t_type == NODE_CLEAN_NODE)
    terminals = ~nodes

    ninner = nodes.sum()
    nouter = terminals.sum()

    for ie in xrange(ne):

        # The condensing equation (conservation of currents):
        #
        # K_ij g_jk f_k = 0,  K_ij = { 1, current j associated with node i,
        #                            { 0, otherwise,
        #                     i \in { nodes }
        #
        # Note that orientation of currents is enforced elsewhere.

        Kg = _n.zeros([2, ninner, 2, g.nnode])

        for i, n in enumerate(_n.where(nodes)[0]):
            for j in node_wires[n].keys():
                Kg[:,i,:,:] += L[ie,:,j,:,:]

        lhs = _n.reshape( Kg[:,:,:,nodes], (2*ninner, 2*ninner))
        rhs = _n.reshape(-Kg[:,:,:,terminals], (2*ninner, 2*nouter))

        f_nodes = scipy.linalg.solve(lhs, rhs)
        f_nodes.shape = (2, ninner, 2, nouter)

        # condense!
        L[ie,...][...,terminals] += _n.tensordot(L[ie,...][...,nodes], f_nodes,
                                                 axes=2)
    L[:,:,:,:,nodes] = 0

    # Normalize to the convention
    return _n.swapaxes(L, 1, 2)

class CacheFifo(object):
    """FIFO cache for numeric array (key, value) pairs.

    This is needed, as scipy.optimize.fsolve initially calls the
    function multiple times with the same arguments: this is costly.
    """
    
    def __init__(self, nmax=1):
        self.nmax = nmax
        self.data = []

    def append(self, key, value):
        """Add a key-value pair (a copy of it) in the fifo."""
        # Copying is absolutely important: passing things by reference
        # invites hard-to-find problems. (In fact, I already
        # encountered such problems...)
        key = _n.array(key, copy=True)
        value = _n.array(value, copy=True)
        for i, (key2, value2) in enumerate(self.data):
            if _n.alltrue(key == key2):
                self.data[i] = [key, value]
                return
        self.data.insert(0, [key, value])
        if len(self.data) > self.nmax:
            self.data.pop()

    def get(self, key):
        for key2, value2 in self.data:
            if _n.alltrue(key == key2):
                return _n.array(value2, copy=True)
        return None

def optimize_parameters_for(x0, zero_func, adjust_func, **kwargs):
    """
    Find a configuration where the specified quantity vanishes,
    by adjusting the specified variables.

    Parameters
    ----------
    x0 : array
        Initial guess for all n parameters
    zero_func : func()
        Function whose zero to look for; should return a vector of the same size
        as `x0`.
    adjust_func : function(x)
        Function to call to set current parameters

    Other keyword arguments accepted by scipy.optimize.fsolve can be
    given.

    Notes
    -----

    If `AlreadyConvergedException` is raised by `zero_func`, the
    optimization terminates successfully with the current value.

    """

    cache = CacheFifo()
    ok = [True, None]

    def eval_func(x):
        """Evaluate the given function."""
        if not ok[0]:
            return 0*x
        y = cache.get(x)
        if not y is None:
            return y
        adjust_func(x)
        try:
            y = _n.array(zero_func(), _n.float_)
            cache.append(x, y)
            return y
        except AlreadyConvergedException:
            return 0*x
        except Exception, e:
            ok[0] = False
            ok[1] = ''.join(traceback.format_exception(
                sys.exc_type, sys.exc_value, sys.exc_traceback))
            return 0*x

    x0 = _n.array(x0, _n.float_)
    try:
        adjust_func(x0)
        f0 = zero_func()
    except AlreadyConvergedException:
        return x0

    if 'epsfcn' not in kwargs:
        # Provide a reasonable initial step size for determining the
        # finite-difference Jacobian.
        kwargs['epsfcn'] = (1 + norm2(f0))/(1 + norm2(x0)) * 1e-3

    cache.append(x0, f0)

    z = scipy.optimize.fsolve(eval_func, x0, full_output=False, **kwargs)

    if not ok[0]:
        raise RuntimeError(ok[1])

    return _n.reshape(_n.array(z), x0.shape)


#############################################################################

class InequivalentDataError(RuntimeError):
    """Indicates that two pieces of data were not calculated for the same parameters."""
    pass

class InequivalentGeometryError(InequivalentDataError):
    """Indicates that two geometries were inequivalent."""
    pass

class CurrentSolver(object):
    """
    Solver for currents flowing in a given geometry.

    Parameters
    ----------
    geometry : Geometry
        The geometry to solve for
    E : array of floats, optional
        Energy points to evaluate all quantities at.
        If omitted, the default is to use a sensibly chosen grid.
    maxE : float, optional
        If `E` is omitted: the maximum energy for the automatic energy grid.
    ne : int, optional
        If `E` is omitted: the number of energy points to use.
    automatic_energy : bool, optional
        Whether to choose an energy grid automatically.
        (Yes, if `E` is omitted.)
    chunksize : int
        How often to display the progress of calculation.
    output_function : func(message)
        Function to use to print any output. If omitted, print to stderr

    Attributes
    ----------
    geometry : Geometry
        Current geometry
    spectral : SpResult
        Solution to the spectral equation
    kinetic : KinResult
        Solution to the kinetic equations
    coefficient : KineticCoefficient
        Coefficients in the kinetic equations
    solver : Solver
        The low-level solver.
    G : array of floats, shape (ne, nw, 2, 2, nnode)
        The spectral conductance/thermoelectric matrix. See eg. [1].

    References
    ----------
    .. [1] Applied Physics A, 89, 625-637 (2007)

    """
    
    def __init__(self, geometry, E=None, chunksize=None,
                 output_function=None, automatic_energy=False,
                 maxE=None, ne=None):
        """Initialize the solver.

        :Parameters:
          - `E`: energy discretization. If ``None``, a sensible choice is
            taken
          - `chunksize`: how many subintervals to divide the spectral
            calculation
          - `output_function`: function(text) that prints output
          - `automatic_energy`: Whether to pick the energy scale used
            for kinetic energy automatically.
          - `maxE`: maximum energy if energy discretization is automatic
          - `ne`: number of points in automatic energy discretization
        """

        if maxE is None: maxE = 500
        if ne is None: ne = 750
    
        if output_function is None:
            output_function = (lambda s: sys.stderr.write(s))

        if E is None:
            self.automatic_energy = True
            E = twolevel_Egrid([maxE], geometry, ne)
        else:
            self.automatic_energy = automatic_energy

        if chunksize is None:
            chunksize = 300
        if chunksize <= 0:
            chunksize = 1e99

        self.geometry = geometry
        self.chunksize = chunksize
        self.output_function = output_function

        self.E = E

        self.kinetic = KinResult()
        self.spectral = SpResult()
        self.coefficient = KineticCoefficients()

        self.G = None
        self.G_E = None
        
        self.solver = Solver()
        self.solver.set_geometry(self.geometry)

    def _output(self, s):
        """Output a string. This is the default output function."""
        self.output_function(s)

    def set_solvers(self, *args, **kwargs):
        """
        Set solver types and tolerances.

        Parameters
        ----------
        sp_solver : int, optional
            Spectral solver to use.
            Choices are ``SP_SOLVER_COLNEW`` (default) and ``SP_SOLVER_TWPBVP``.
        sp_tol : float, optional
            Tolerance to use in the spectral solver.
        kin_solver : int, optional
            Kineticsolver to use.
            Choices are ``KIN_SOLVER_COLNEW`` and ``KIN_SOLVER_TWPBVP`` (default),
            and ``KIN_SOLVER_BLOCK``.
        kin_tol : float, optional
            Tolerance to use in the kinetic solver.

        """
        self.solver.set_solvers(*args, **kwargs)

    def solve_spectral(self):
        """Solve the spectral equations.

        Store the results to `self.spectral` and `self.kinetic`.
        """

        E = self.E

        x = self.geometry.x

        nw = self.geometry.nwire
        ne = len(E)
        nx = len(x)

        self.spectral.setshape(ne, nw, nx)
        self.coefficient.setshape(ne, nw, nx)
        
        self.solver.set_geometry(self.geometry)

        continued = False

        # Note: continuation works best when solving from high energies toward
        #       low energies; hence spectral equations are solved in this way.
        
        k = len(E)
        while k > 0:
            maxk = k
            mink = k - min(self.chunksize, k)
            k = mink

            self._output("%% E = %11.5g ... %11.5g (%5d, %d) \t[ %3d %% ]\t"%(
                E[maxk-1], E[mink], ne - mink, (maxk - mink),
                int(float(ne - maxk) / (ne-1) * 100.0 )))

            # Solve spectral equations
            self._output(" sp")
            rs = self.solver.sp_solve(E[mink:maxk], x, continued)
            self.spectral[mink:maxk,:,:] = rs

            # Evaluate coefficients for kinetic equations
            self._output(" coef")
            coefs = KineticCoefficients(self.geometry, rs)
            #coefs.regularize(self.solver.get_subgap_bcs(E[mink:maxk]))
            self.coefficient[mink:maxk,:,:] = coefs

            # Continue
            continued = True
            self._output("\n")

        self.solver.set_kinetic(self.coefficient)
    
    def solve_spectral_if_needed(self, calculate_G=True):
        """
        Solve for spectral data if there is none yet.

        Parameters
        ----------
        calculate_G : bool, optional
            Whether to calculate the spectral conductance/thermoelectric matrix.

        """
        need_save = False
        if self.coefficient.DL is None:
            self.solve_spectral()
            need_save = True
        if self.G is None and calculate_G:
            self.calculate_G()
            need_save = True
        return need_save


    def solve_spectral_and_save_if_needed(self, filename,
                                          calculate_G=True, **kw):
        """
        Solve for spectral data if there is none yet, and save the result to a file.

        Parameters
        ----------
        filename : str
            Name of the file to solve data to.
        calculate_G : bool, optional
            Whether to calculate the spectral conductance/thermoelectric matrix.

        Other parameters as in `save`.

        """
        if self.solve_spectral_if_needed(calculate_G):
            self.save(filename, **kw)

    def solve_kinetic(self, E=None, ne=None):
        """
        Solve the kinetic equations and store the result to `self.kinetic`.

        Parameters
        ----------
        E : array of floats, optional
            The energy discretization to use. None indicates that
            use either the same as for spectral, or, if
            `self.automatic_energy` is ``True``, pick a sensible choice.
        ne : int, optional
            How many points to use in energy discretization, if
            choosing the energy discretization automatically.

        """
        E = self._prepare_kinetic(E, ne)
        self.kinetic = self.solver.kin_solve(E, self.geometry.x)

    def _prepare_kinetic(self, E=None, ne=None):
        self.solver.set_geometry(self.geometry, True)
        if E is None:
            if self.automatic_energy:
                if ne is None: ne = len(self.E) - 1 # -1 is for catching errors
                E = twolevel_Egrid(self.E, self.geometry, ne)
            else:
                E = self.E
        return E

    def solve(self):
        """Solve both spectral and kinetic equations."""
        self.solve_spectral()
        self.solve_kinetic()

    @classmethod
    def load(cls, src, path=None, **kw):
        """Load data from an open HDF5 file.

        :Parameters:
          - `src`: a tables.Group where to load data from,
            or a file name joined with node path.
          - Other arguments (except ``geometry``) are passed on to __init__
        """

        if isinstance(src, str):
            own_file = True
            file, src = data.find_node_in_file(src)
        else:
            own_file = False
            if isinstance(src, tables.File):
                file = src
                if not hasattr(file, 'root'):
                    raise IOError("File does not have root")
                src = file.root
            else:
                file = src._v_file

        if path:
            src = file.getNode('/'.join([src._v_pathname, path]))

        try:
            geometry = data.load(src, 'geometry')

            kw.setdefault('automatic_energy', False)
            self = cls(geometry, **kw)

            self.E = None
            if hasattr(src, 'spectral'):
                self.spectral = data.load(src, 'spectral')
                if self.E is None and self.spectral.E is not None:
                    self.E = self.spectral.E

            if hasattr(src, 'coefficient'):
                self.coefficient = data.load(src, 'coefficient')
                if self.coefficient.E is not None:
                    self.solver.set_kinetic(self.coefficient)
                if self.E is None and self.kinetic.E is not None:
                    self.E = self.coefficient.E
            
            if hasattr(src, 'kinetic'):
                self.kinetic = data.load(src, 'kinetic')
                if self.E is None and self.kinetic.E is not None:
                    self.E = self.kinetic.E

            if (hasattr(src, 'G') and hasattr(src.G, 'G')
                and hasattr(src.G, 'E')):
                
                self.G = data.load(src, 'G/G')
                self.G_E = data.load(src, 'G/E')
        finally:
            if own_file:
                del src
                file.close()
        return self

    @classmethod
    def resume(cls, src, geometry, path=None, compare_E=False, **kw):
        """Load data from a data file if parameters match or return \
        a new solver if not."""

        g_kw = {}

        for name in ('compare_delta', 'compare_kinetic', 'compare_phases'):
            if name in kw:
                g_kw[name] = kw[name]
                del kw[name]

        try:
            solver = cls.load(src, path=path, **kw)
        except (IOError, LookupError):
            solver = None

        if not solver:
            return cls(geometry, **kw)
        else:
            new_solver = cls(geometry, **kw)

        equal = (new_solver.geometry.equal(solver.geometry, **g_kw)
                 and (not compare_E or _n.alltrue(new_solver.E == solver.E))
                 )
        if not equal:
            solver = new_solver

        # Ensure that the solver is a reasonable state.
        # Note: discard the loaded geometry, since a typical usage pattern
        #       is to rely on the fact that the geometry object is the same
        solver.geometry = geometry
        solver.solver.set_geometry(solver.geometry)
        if solver.coefficient.DL is not None:
            solver.solver.set_kinetic(solver.coefficient)
        return solver

    def save(self, to, path=None, save_coefficient=True, save_spectral=True,
             save_kinetic=True):
        """
        Save date to an open HDF5 file.

        Parameters
        ----------
        to: str or tables.File
            A file name or a parent HDF5 node to save to
        path : str, optional
            HDF5 name to save under
        save_coefficient : bool, optional
            Whether to save kinetic coefficients
        save_spectral : bool, optional
            Whether to save spectral data
        save_kinetic : bool, optional
            Whether to save kinetic data.

        Notes
        -----
        Geometry is saved in all cases.

        """

        if isinstance(to, str):
            file = tables.openFile(to, 'w')
            to = file.root
        else:
            file = None

        if path:
            to = data.create_group(to, path)

        try:
            data.save(to, 'geometry', self.geometry)
            if save_spectral:
                data.save(to, 'spectral', self.spectral)
            if save_coefficient:
                data.save(to, 'coefficient', self.coefficient)
            if save_kinetic:
                data.save(to, 'kinetic', self.kinetic)
            if self.G != None:
                data.save(to, 'G/G', self.G)
                data.save(to, 'G/E', self.G_E)
        finally:
            if file is not None:
                del to
                file.close()

    def get_currents(self, E=None, ix=None, w_jL=None, w_jT=None, ne=None):
        """
        Calculate currents, using fixed-grid discretization.

        Parameters
        ----------
        ix : int, optional
            Index of position to evaluate currents at in each wire.
            If omitted, current is evaluated at all positions.
        w_jT : list of int, optional
            For which wires to compute charge current for. If omitted, compute
            for all wires.
        w_jL : list of int, optional
            For which wires to compute energy current for. If omitted, compute
            for all wires.
        ne : int, optional
            How many points to put in energy discretization, if chosen automatically.


        Returns
        -------
        Ic : array of floats, shape (nwire, nx')
            Charge current in each wire. Contains `nan` in entries that were
            not calculated. If `ix` was given, `nx' == 1`, otherwise `nx' == nx`.
        IE : array of floats, shape (nwire, nx')
            Energy current in each wire.

        """
        self.solve_kinetic(ne=ne)
        return calculate_currents(self.geometry, self.kinetic,
                                  ix=ix, w_jL=w_jL, w_jT=w_jT)

    def get_currents_lazy(self, w_jL=None, w_jT=None,
                          epsabs=None, epsrel=None, ix=None):
        """
        Calculate currents, using an adaptive integrator.

        Will solve kinetic equations at energy points where the solutions
        are needed.

        Parameters
        ----------
        epsabs : float, optional
            Absolute tolerance for the currents.
        epsrel : float, optional
            Relative tolerance for the currents.

        Also takes the same parameters as `get_currents`.

        Returns
        -------
        Returns similar output as `get_currents`.

        """
        self.solver.set_geometry(self.geometry, True)
        return calculate_currents_lazy(self.geometry, self.coefficient,
                                       self.solver, w_jL=w_jL, w_jT=w_jT,
                                       epsabs=epsabs, epsrel=epsrel)


    def _get_all_terminals(self):
        return _find_all_terminals(self.geometry)

    def calculate_G(self, superconductors_in_equilibrium=False,
                    only_terminals=None):
        """
        Compute the spectral conductance matrix.

        .. warning::

           Spectral quantities (such as G) are not conserved in
           superconducting wires. Be aware that the G is evaluated
           at x[0] in each wire.

        Parameters
        ----------
        superconductors_in_equilibrium : bool, optional
            Assume superconductors are at equilibrium, and skip calculating
            the conductance matrix entries for them.
        only_terminals : list of int, optional
            Terminals for which to calculate the conductance matrix entries.
            If omitted, calculate for all terminals.

        """
        g = self.geometry
        E = self.E
        nw = g.nwire
        ne = len(E)

        self.solver.set_geometry(self.geometry, True)
        self.solver.set_kinetic(self.coefficient)

        terminals = self._get_all_terminals()

        if only_terminals is not None:
            terminals = filter(lambda x: x in only_terminals, terminals)

        # Check sanity
        if not terminals:
            raise ValueError('Some real terminals must be present')

        # Calculate the characteristic currents
        clear_mask = BCSUB_MASK

        self.G = (_n.zeros([ne, nw, 2, 2, g.nnode], _n.float_)
                       + _n.nan)
        self.G_E = E

        self._output('%% Evaluating G matrix...')

        for index, mask in [(0, BCSUB_CHARACTERISTIC_L),
                            (1, BCSUB_CHARACTERISTIC_T)]:

            for t in terminals:
                # Push in the characteristic specification
                g.t_type &= ~clear_mask
                for t2 in terminals:
                    if t2 == t: continue
                    g.t_type[t2] |= BCSUB_CHARACTERISTIC_ZERO
                g.t_type[t] |= mask

                if index == 0:
                    self._output(" L%d" % t)
                else:
                    self._output(" T%d" % t)

                if superconductors_in_equilibrium:
                    if (g.t_type[t] & BCTYPE_MASK) == NODE_CLEAN_S_TERMINAL \
                           and index == 1:
                        self.G[:,:,:,index,t] = 0
                        continue

                # Solve
                self.solve_kinetic(E=E)

                # Get the results
                self.G[:,:,0,index,t] = self.kinetic.jL[:,:,0]
                self.G[:,:,1,index,t] = self.kinetic.jT[:,:,0]

        self._output(".\n")
        g.t_type &= ~clear_mask # cleanup

    def approximate_G(self, use_jS=True, use_T=True):
        """
        Find an approximation for the G-matrix,
        up to first order in :math:`j_S` and :math:`{\\cal T}`.

        Stores the result to `self.G`.

        Parameters
        ----------
        use_jS : bool, optional
            Use the spectral supercurrent in the approximation.
        use_T : bool, optional
            Use the coefficient :math:`{\\cal T}` in the approximation.

        """
        self.G = approximate_G(self.geometry, self.coefficient,
                               use_jS=use_jS, use_T=use_T)
        self.G_E = self.coefficient.E

    def get_currents_from_G_with_f(self, f_func,
                                   w_jT=None, w_jL=None, ix=None):
        """
        Calculate currents from the G-matrix, for given distribution functions
        at the terminals.

        Parameters
        ----------
        f_func : func(E, mu, T)
            Function that takes an energy grid of shape (ne, 1) and potentials
            and temperatures [shape (1, nnode)] and returns (fL, fT) where each
            entry is an array of shape (ne, nnode) describing the distribution
            function at each energy in each terminal.
        w_jT : list of int, optional
            For which wires to compute charge current for. If omitted, compute
            for all wires.
        w_jL : list of int, optional
            For which wires to compute energy current for. If omitted, compute
            for all wires.

        Returns
        -------
        Ic : array of floats, shape (nwire,)
            Charge current in each wire. Contains `nan` in entries that were
            not calculated.
        IE : array of floats, shape (nwire,)
            Energy current in each wire.

        """
        return calculate_currents_from_G_with_f(self.geometry,
                                                self.G_E, self.G,
                                                f_func, w_jT, w_jL)

    def get_linear_response_from_G_with_f(self, df_func, w_jT=None, w_jL=None):
        """
        Compute the linear response in currents to changes in distribution
        functions in the terminals.

        Parameters
        ----------
        df_func : func(E, mu, T)
            Function that takes an energy grid of shape (ne, 1) and potentials
            and temperatures [shape (1, nnode)] and returns
            `(dfL_dT, dfL_dV, dfT_dT, dfT_dV)` where each entry is an array
            of shape (ne, nnode) describing the change in the distribution function
            with regard to perturbation in each parameter.
        w_jT : list of int, optional
            For which wires to evaluate the charge current related entries.
        w_jL : list of int, optional
            For which wires to evaluate the energy current related entries.

        Returns
        -------
        G_VV : array, shape (nwire, nnode)
            Conductance matrix; contains entries :math:`dI_{c,i}/dV_j`.
            Has `nan` at entries that were not calculated.
        G_VT : array, shape (nwire, nnode)
            Thermoelectric matrix; contains entries :math:`dI_{c,i}/dT_j`.
        G_TV : array, shape (nwire, nnode)
            Thermoelectric matrix; contains entries :math:`dI_{E,i}/dV_j`.
        G_TT : array, shape (nwire, nnode)
            Thermal conductance matrix; contains entries :math:`dI_{E,i}/dT_j`.

        """
        return calculate_linear_response_from_G_with_f(self.geometry,
                                                       self.G_E, self.G,
                                                       df_func,
                                                       w_jT, w_jL)

    def get_currents_from_G(self, *args, **kw):
        """
        Calculate currents, from the spectral conductance matrix.

        Parameters
        ----------
        w_jT : list of int, optional
            For which wires to compute charge current for. If omitted, compute
            for all wires.
        w_jL : list of int, optional
            For which wires to compute energy current for. If omitted, compute
            for all wires.

        Returns
        -------
        Ic : array of floats, shape (nwire,)
            Charge current in each wire. Contains `nan` in entries that were
            not calculated.
        IE : array of floats, shape (nwire,)
            Energy current in each wire.

        """
        def equilibrium_terminal_f(E, mu, T):
            fL = .5*(_n.tanh(.5*(E + mu)/T) + _n.tanh(.5*(E - mu)/T))
            fT = .5*(_n.tanh(.5*(E + mu)/T) - _n.tanh(.5*(E - mu)/T))
            return fL, fT
        return self.get_currents_from_G_with_f(equilibrium_terminal_f,
                                               *args, **kw)

    def get_linear_response_from_G(self, *args, **kw):
        """
        Compute the linear response in currents to changes in potentials and
        temperatures in the terminals.

        This function assumes that the distribution function in each terminal
        is a Fermi function.

        Parameters
        ----------
        w_jT : list of int, optional
            For which wires to evaluate the charge current related entries.
        w_jL : list of int, optional
            For which wires to evaluate the energy current related entries.

        Returns
        -------
        G_VV : array, shape (nwire, nnode)
            Conductance matrix; contains entries :math:`dI_{c,i}/dV_j`.
            Has `nan` at entries that were not calculated.
        G_VT : array, shape (nwire, nnode)
            Thermoelectric matrix; contains entries :math:`dI_{c,i}/dT_j`.
        G_TV : array, shape (nwire, nnode)
            Thermoelectric matrix; contains entries :math:`dI_{E,i}/dV_j`.
        G_TT : array, shape (nwire, nnode)
            Thermal conductance matrix; contains entries :math:`dI_{E,i}/dT_j`.

        """
        def equilibrium_terminal_df(E, mu, T):
            sechp = _n.cosh(.5*(E + mu)/T)**(-2)
            sechm = _n.cosh(.5*(E - mu)/T)**(-2)
            sechpv = sechp * (E + mu)
            sechmv = sechm * (E - mu)
            
            dfL_dT = -.25/T**2 * (sechpv + sechmv)
            dfT_dT = -.25/T**2 * (sechpv - sechmv)
            dfL_dV =  .25/T    * (sechp  - sechm )
            dfT_dV =  .25/T    * (sechp  + sechm )
            return dfL_dT, dfL_dV, dfT_dT, dfT_dV
        return self.get_linear_response_from_G_with_f(equilibrium_terminal_df,
                                                      *args, **kw)

def twolevel_Egrid(spectral_E, geometry, n):
    """
    Generate a grid on the energy scales of both kinetic and spectral
    equations.

    This uses "expert knowledge" to weight energy scales of

      - 1:         assumed Thouless energy
      - t_t, t_mu: temperatures in kinetic equations
      - t_delta:   increase density near terminal energy gaps

    """
    sp_maxE = max(spectral_E)
    kin_maxE = max(abs(geometry.t_mu)) + 10 * max(abs(geometry.t_t))
    maxE = max(sp_maxE, kin_maxE)

    minE = 1e-4

    peaks = []

    if maxE > sp_maxE:
        sys.stderr.write(" Warning:maxe ")
        maxE = sp_maxE

    if kin_maxE < sp_maxE and kin_maxE > 0:
        peaks.append((minE + .5*kin_maxE, 1, kin_maxE*.5))

    jspeak = 2
    if jspeak > max(abs(geometry.t_delta)):
        jspeak = max(abs(geometry.t_delta)) + 1e-5

    ns = 0
    for i in range(geometry.nwire):
        if (geometry.t_type[i] & BCTYPE_MASK) == NODE_CLEAN_S_TERMINAL:
            ns += 1

    for i in range(geometry.nwire):
        if (geometry.t_type[i] & BCTYPE_MASK) == NODE_CLEAN_S_TERMINAL:
            peaks.append((minE + geometry.t_delta[i], 0.5/ns, jspeak))

    peaks.append((minE, 1, 20 * jspeak))
    peaks.append((jspeak, 1, 5 * jspeak))
    peaks.append((minE, 0.3, 1 * jspeak))
    peaks.append((minE, 0.3, maxE))

    return linspace_weighed(minE, maxE, n, peaks)
