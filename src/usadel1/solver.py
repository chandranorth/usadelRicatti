"""
:Author:       Pauli Virtanen <pauli@ltl.tkk.fi>
:Organization: Low Temperature Laboratory, Helsinki University of Technology
:Date:         2005-2006

Spectral and kinetic Usadel equation solvers.

This is the main interface for solving spectral and kinetic Usadel
equations.

Some use cases:

- `currents`: Calculating currents
- `selfconsistentiteration`: Self consistent iteration.


Usage
-----

The general usage is

1.  First initialize a `Geometry` object

       >>> x = linspace(0, 1, 100)
       >>> geometry = Geometry(nwire, nnode, len(x))
       >>> geometry.x = x
       ...

    Also set its other parameters.

2.  Then create a `Solver` object and initialize it

        >>> solver = Solver()
        >>> solver.set_geometry(geometry)

3.  Solve the spectral equations at some energies

        >>> E = linspace(0, 100, 750)
        >>> spectral_solution = solver.sp_solve(E, x)

    The result is a `SpResult` object.

4.  Calculate `KineticCoefficients` for the kinetic equations

        >>> coefficient = KineticCoefficients(geometry, spectral_solution)
        >>> solver.set_kinetic(coefficient)

5.  Solve kinetic equations

        >>> kinetic_solution = solver.kin_solve(E, x)

    The result is a `KinResult` object.

After this, you can

- Use tools in `currents` for integrating observable currents:

      >>> Ic, IE = usadel1.currents.calculate_currents(geometry,
      ...                                              kinetic_solution)

- Solve optimization problems for finding for example the thermopower etc.
  See `currents.optimize_parameters_for`.

- Save data to files using `data`:

      >>> import tables
      >>> file = tables.openFile("datafile.h5", mode="w")
      >>> usadel1.data.save(file.root, 'geometry', geometry)
      >>> usadel1.data.save(file.root, 'spectral', spectral_solution)
      >>> usadel1.data.save(file.root, 'kinetic', kinetic_solution)
      >>> usadel1.data.save(file.root, 'Ic', Ic)
      >>> usadel1.data.save(file.root, 'IE', IE)
      >>> file.close()

  The data files are saved in the standard HDF5_ format.
  You can load these e.g. in Matlab [#]_, most easily by using the
  ``loadhdf5struct.m`` m-file bundled with usadel1.

.. [#] But note that Matlab versions at least prior to R14SP3 have
       bugs in the HDF5 support.
.. _HDF5: http://hdf.ncsa.uiuc.edu/HDF5/

Internal logic
--------------

Part of the boundary condition logic is written in Python, embedded in
`Solver`, `EquationSystem`, and `BoundaryConditionList`.

The main reason is that `_solvercore` expects specially ordered input
variables ``bctype``, ``bcconnect``, and ``nrec``. These are defined as
follows:

============= ===============================================================
``bctype``    (``mstar``) vector of type codes (``solvercore.BCTYPE_*``)
              for boundary conditions.

``bcconnect`` a (``nwire+1``, ``mstar``) array of data used by boundary
              conditions in ``kin_equations.f90`` and ``sp_equations.f90``.
              The contents of ``bcconnect(:, j)`` for the corresponding
              ``bctype(j)`` are:

              - ``BCTYPE_*_TERMINAL``: [ wire, terminal, 0, 0, ... ]
              - ``BCTYPE_CLEAN_NODE``: [ wire_1, wire_2, ..., wire_r, 0,
                0, ... ]

``nrec``      number of boundary conditions in ``bcconnect`` corresponding to
              end-points of wires.
============= ===============================================================

Note that all indices should be in Fortran form, i.e., 1-based.

This code forms the above arrays using `BoundaryConditionList`, which
determines on its initialization the proper boundary conditions. This
data is then used by `EquationSystem.get_bc` which formats it to the form
expected by the Fortran routines.

"""


from __future__ import division

import solvercore as _sc
import numpy as _n
import scipy as _scipy

from util import ArrayProperty, Sliceable, linspace_weighed, \
                 divide_vector_to_chunks


__docformat__ = "restructuredtext en"
__all__ = ['Geometry', 'Solver', 'SpResult', 'KinResult',
           'KineticCoefficients',
           'NODE_CLEAN_N_TERMINAL',
           'NODE_CLEAN_S_TERMINAL',
           'NODE_CLEAN_S_TERMINAL_CIB',
           'NODE_CLEAN_NODE',
           'NODE_FREE_INTERFACE',
           'NODE_FREE_INTERFACE_FIX_PHASE',
           'NODE_TUNNEL_N_TERMINAL',
           'NODE_TUNNEL_S_TERMINAL',
           'NODE_TUNNEL_NODE',
           'WIRE_TYPE_N',
           'WIRE_TYPE_S',
           'WIRE_TYPE_NULL',
           'SP_SOLVER_DEFAULT',
           'SP_SOLVER_COLNEW',
           'SP_SOLVER_TWPBVP',
           'KIN_SOLVER_DEFAULT',
           'KIN_SOLVER_COLNEW',
           'KIN_SOLVER_TWPBVP',
           'KIN_SOLVER_BLOCK',
           'BCTYPE_MASK',
           'BCSUB_MASK',
           'BCSUB_CHARACTERISTIC_ZERO',
           'BCSUB_CHARACTERISTIC_L',
           'BCSUB_CHARACTERISTIC_T',
           ]

WIRE_TYPE_N = 0
WIRE_TYPE_S = 1
WIRE_TYPE_NULL = 2

BCTYPE_CLEAN_N_TERMINAL         = 1
BCTYPE_CLEAN_S_TERMINAL         = 2
BCTYPE_CLEAN_NODE               = 3
BCTYPE_CLEAN_S_TERMINAL_CIB     = 4
BCTYPE_FREE_INTERFACE           = 5
BCTYPE_FREE_INTERFACE_FIX_PHASE = 6
BCTYPE_TUNNEL_N_TERMINAL        = 7
BCTYPE_TUNNEL_S_TERMINAL        = 8
BCTYPE_TUNNEL_NODE              = 9

BCSUB_CHARACTERISTIC_ZERO      =  1 << 8
BCSUB_CHARACTERISTIC_L         =  2 << 8
BCSUB_CHARACTERISTIC_T         =  3 << 8

BCTYPE_MASK = 0x00ff
BCSUB_MASK  = 0xff00

NODE_CLEAN_N_TERMINAL         = BCTYPE_CLEAN_N_TERMINAL
NODE_CLEAN_S_TERMINAL         = BCTYPE_CLEAN_S_TERMINAL
NODE_CLEAN_NODE               = BCTYPE_CLEAN_NODE
NODE_CLEAN_S_TERMINAL_CIB     = BCTYPE_CLEAN_S_TERMINAL_CIB
NODE_FREE_INTERFACE           = BCTYPE_FREE_INTERFACE
NODE_FREE_INTERFACE_FIX_PHASE = BCTYPE_FREE_INTERFACE_FIX_PHASE
NODE_TUNNEL_N_TERMINAL        = BCTYPE_TUNNEL_N_TERMINAL
NODE_TUNNEL_S_TERMINAL        = BCTYPE_TUNNEL_S_TERMINAL
NODE_TUNNEL_NODE              = BCTYPE_TUNNEL_NODE

## Solver types
SP_SOLVER_COLNEW  = 1
SP_SOLVER_TWPBVP  = 2
SP_SOLVER_DEFAULT = SP_SOLVER_COLNEW

KIN_SOLVER_COLNEW  = 1
KIN_SOLVER_TWPBVP  = 2
KIN_SOLVER_BLOCK   = 3
KIN_SOLVER_DEFAULT = KIN_SOLVER_TWPBVP


##############################################################################
class Geometry(object):
    """
    Geometry of the setup.

    Parameters
    ----------
    nwire : int
        Number of wires in the setup
    nnode : int
        Number of nodes in the setup
    x : array of floats, optional
        Array of x-positions at which quantities are specified.
        These must be in the range [0, 1].

    Examples
    --------
    A three-probe structure::

                  3
           0 S----*----S 1
               0  |  1
                  |
                  |2
                  |
                  N
                  2

       >>> g = Geometry(nwire=3, nnode=4)
       >>> g.t_type = [NODE_CLEAN_S_TERMINAL,
       ...             NODE_CLEAN_S_TERMINAL,
       ...             NODE_CLEAN_NODE]
       >>> g.w_type = WIRE_TYPE_N
       >>> g.w_ends[0,:] = [0, 3]
       >>> g.w_ends[1,:] = [1, 3]
       >>> g.w_ends[2,:] = [2, 3]
       >>> g.w_length = [0.5, 0.5, 2]
       >>> g.t_delta = [100, 100, 0]
       >>> g.t_phase = array([-.5, .5, 0]) * pi/2
       >>> g.t_t = 1
       >>> g.t_mu = 0

    """

    t_type = ArrayProperty("t_type", """
        array(nnode) of types of nodes

        Possible choices are:
        NODE_CLEAN_N_TERMINAL (normal terminal; clean interface),
        NODE_CLEAN_S_TERMINAL (superconducting terminal; clean interface),
        NODE_CLEAN_NODE (node connecting several wires; clean interface),
        NODE_FREE_INTERFACE (vacuum interface).
        """)
    t_delta = ArrayProperty("t_delta", "array(nnode) of energy gaps of nodes")
    t_phase = ArrayProperty("t_phase", "array(nnode) of phases of nodes")
    t_inelastic = ArrayProperty("t_inelastic",
                                "array(nnode) of :math:`\Gamma`-s of nodes")
    t_spinflip = ArrayProperty("t_spinflip", """
        FIXME: this is not yet functional
        array(nnode) of :math:`\Gamma_{sf}`-s of node""")
    t_t = ArrayProperty("t_t", "array(nnode) of node temperatures")
    t_mu = ArrayProperty("t_mu", "array(nnode) of node potentials")
    t_resistance = ArrayProperty("t_resistance", """
        array(nnode) of node resistances

        Has effect only for NODE_TUNNEL_* nodes.
        """)
    w_type = ArrayProperty("w_type", """
        array(nwire) of wire types

        Possible choices are: WIRE_TYPE_N (normal wire) and
        WIRE_TYPE_S (superconducting wire).
        """)
    w_length = ArrayProperty("w_length", "array(nwire) of wire lengths")
    w_conductance = ArrayProperty("w_conductance", """
        array(nwire) of conductance*area of wires""")
    w_inelastic = ArrayProperty("w_inelastic",
                                "array(nwire) of :math:`\Gamma`-s of wires")
    w_spinflip = ArrayProperty("w_spinflip", """
        array(nwire) of :math:`\Gamma_{sf}`-s of wires""")
    w_ends = ArrayProperty("w_ends",
                           "array(nwire,2) that maps wire ends -> nodes")
    coupling_lambda = ArrayProperty("coupling_lambda", """
        array(nwire, nx) of :math:`\lambda` in wires""")
    omega_D = ArrayProperty("omega_D",
                            "array(nwire, nx) of Debye frequency in wires")
    x = ArrayProperty("x", """
        array(nx) of discretization points.

        NOTE: These must be in the range [0,1] !

        Note also that this does not affect the actual mesh chosen
        by all solvers -- however, it always affects the discretization
        used for spectral solutions and kinetic coefficients.""")
    w_delta = ArrayProperty("w_delta",
                            "array(nwire, nx) of energy gaps in wires")
    w_phase = ArrayProperty("w_phase", "array(nwire, nx) of phase in wires")
    nnode = ArrayProperty("nnode", "how many nodes in geometry")
    nwire = ArrayProperty("nwire", "how many wires in geometry")

    w_phase_jump = ArrayProperty("w_phase_jump", """
       array(nwire) phase jump between the ends of the wire, due to vector
       potential parallel to the wire.""")

    def __init__(self, nwire, nnode, x=None):
        """Initialize a geometry.

        :Parameters:
          - `nwire`: how many wires in geometry
          - `nnode`: how many nodes in geometry
          - `x`: the x-discretization.
            NOTE: These must be in the range [0,1] !
            If ``None``, a sensible choice is picked.
        """
        self.nwire = nwire
        self.nnode = nnode

        if x is None:
            x = linspace_weighed(0, 1, 101, [(0, 1, 0.1), (1, 1, 0.1)])

        nx = len(x)

        self.t_type = _n.zeros([nnode], _n.int_)
        self.t_delta = _n.zeros([nnode], _n.float_)
        self.t_phase = _n.zeros([nnode], _n.float_)
        self.t_inelastic = _n.zeros([nnode], _n.float_)
        self.t_resistance = _n.zeros([nnode], _n.float_)
        self.t_spinflip = _n.zeros([nnode], _n.float_)
        self.t_t = _n.zeros([nnode], _n.float_)
        self.t_mu = _n.zeros([nnode], _n.float_)

        self.w_type = _n.zeros([nwire], _n.int_)
        self.w_length = _n.zeros([nwire], _n.float_)
        self.w_conductance = _n.zeros([nwire], _n.float_)
        self.w_inelastic = _n.zeros([nwire], _n.float_)
        self.w_spinflip = _n.zeros([nwire], _n.float_)

        self.w_ends = _n.zeros([nwire, 2], _n.int_)

        self.coupling_lambda = _n.zeros([nwire, nx], _n.float_)
        self.omega_D = _n.zeros([nwire, nx], _n.float_)
        self.x = _n.array(x, _n.float_, copy=True)
        self.w_delta = _n.zeros([nwire, nx], _n.float_)
        self.w_phase = _n.zeros([nwire, nx], _n.float_)
        self.w_phase_jump = _n.zeros([nwire], _n.float_)

    def equal(self, other, compare_delta=True, compare_kinetic=True,
              compare_phases=True):
        """
        Compare two geometries.

        .. rubric:: Parameters

        other : Geometry
            The geometry to compare this one to.
        compare_kinetic : bool, optional
            Whether to compare kinetic quantities (temperatures, potentials)
        compare_delta : bool, optional
            Whether to compare :math:`|\\Delta|` between the Geometries.
        compare_phases : bool, optional
            Whether to compare :math:`{\\rm arg}\\Delta` between the Geometries.

        """
        delta_vars = ('w_delta', 'w_phase')
        kinetic_vars = ('t_t', 't_mu')
        phase_vars = ('t_phase', 'w_phase')

        try:
            names = [n for n, v in self.__class__.__dict__.items()
                     if isinstance(v, ArrayProperty)]

            for name in names + ["nwire", "nnode"]:
                if not compare_delta and name in delta_vars:
                    continue
                if not compare_kinetic and name in kinetic_vars:
                    continue
                if not compare_phases and name in phase_vars:
                    continue
                v1 = getattr(self, name)
                v2 = getattr(other, name)
                if not _n.alltrue(_n.ravel(v1 == v2)):
                    return False
        except ValueError:
            return False
        return True

    def __eq__(self, other):
        """Compare geometries."""
        return equal(self, other)

    def get_node_wires(self):
        """Collect lists of wires connected to each node, and label nodes
           by which end of the wires they correspond to."""
        nodewires = [[] for i in range(self.nnode)]
        nodetype = [-1] * self.nnode

        for k in range(self.nwire):
            for j in range(2):
                ni = self.w_ends[k, j]
                nodewires[ni].append(k)
                if nodetype[ni] < 0:
                    nodetype[ni] = j
                elif nodetype[ni] != j:
                    raise RuntimeError, \
                          ("Wire ends for wires %d and %d do not match " +
                           "at node %d. Add dummy wires, if necessary.") \
                           % ( nodewires[ni][0], k, ni )

        return (nodewires, nodetype)

    def get_idstr(self):
        """Get a descriptive string for this geometry."""
        strD = ""
        strP = ""
        strL = ""
        strA = ""
        strV = ""
        strT = ""
        for w in range(self.nwire):
            strL = "%s  %7g" % (strL, self.w_length[w])
            strA = "%s  %7g" % (strA, self.w_conductance[w])
        for w in range(self.nnode):
            if (self.t_type[w] & BCTYPE_MASK) == NODE_CLEAN_NODE:
                continue
            strD = "%s  %7g" % (strD, self.t_delta[w])
            strP = "%s  %7g" % (strP, self.t_phase[w]/_n.pi)
            strV = "%s  %7g" % (strV, self.t_mu[w])
            strT = "%s  %7g" % (strT, self.t_t[w])

        s = "Delta  = %s\nPhase  = %s\nLength = %s\nArea   = %s\nV      = %s\nT      = %s" % (strD, strP, strL, strA, strV, strT)
        return s

    def hash_this(self, h):
        names = [n for n, v in self.__class__.__dict__.items()
                 if isinstance(v, ArrayProperty)]
        for name in names:
            hash_array(h, getattr(self, name))

    def copy(self):
        other = Geometry(nwire=self.nwire, nnode=self.nnode, x=self.x)
        for field in ('t_type', 't_delta', 't_phase', 't_inelastic',
                      't_spinflip', 't_t', 't_mu', 't_resistance',
                      'w_type', 'w_length', 'w_conductance', 'w_inelastic',
                      'w_spinflip', 'w_ends', 'coupling_lambda', 'omega_D',
                      'w_delta', 'w_phase', 'w_phase_jump'):
            setattr(other, field, getattr(self, field))
        return other

def hash_array(h, a):
    a = _n.asarray(a)
    shape = a.shape
    h.update(repr(a.shape))
    for k in _n.ravel(a):
        h.update(repr(k))

##############################################################################
class Solver(object):
    """
    Low-level interface to the Usadel solver.

    """

    __last_user = 0
    __next_id = 0

    def __init__(self):
        """Initialize the solver."""
        self.__geometry = None

        self.__params_set = False
        self.__delta_set = False
        self.__kinetic_set = False

        self.__sp_initialized = False
        self.__kin_initialized = False

        self.__last_E_sp = None
        self.__last_E_kin = None

        self.equations = None

        self.__id = Solver.__next_id
        Solver.__next_id += 1

        self.__solver_settings = dict(
            sp_solver=SP_SOLVER_COLNEW,
            kin_solver=KIN_SOLVER_TWPBVP,
            sp_tol=0.0,
            kin_tol=0.0
        )

    def __check_serial(self):
        """The core solver is a singleton, but we want to provide
           an object-like interface for it. This method checks
           whether we need to reset the core."""
        if Solver.__last_user != self.__id:
            self.__sp_initialized = False
            self.__kin_initialized = False
            self.__delta_set = False
            self.__params_set = False
            self.__kinetic_set = False
            Solver.__last_user = self.__id


    def set_geometry(self, g, preserve=False):
        """
        Set the geometry to use.

        Parameters
        ----------
        g : Geometry
            The geometry to use
        preserve : bool, optional
            Whether to preserve currently set values for :math:`\\Delta`
            and the kinetic coefficients.

        """
        self.__check_serial()

        self.__params_set = True
        self.__geometry = g

        if not preserve:
            self.__delta_set = False
            self.__kinetic_set = False

            self.equations = EquationSystem(self.__geometry, 0)
            self.equations.sanity_check_geometry()

        _sc.set_params(
            g.w_length, g.w_conductance, g.w_inelastic, g.w_spinflip,
            g.w_phase_jump,
            g.t_delta, g.t_phase, g.t_inelastic, g.t_spinflip, g.t_t, g.t_mu)

    def __reset_delta(self):
        """Set the superconducting order parameter"""
        self.__check_serial()

        if not self.__params_set:
            raise RuntimeError, "spectralsolver geometry not set"
        g = self.__geometry

        nx = _n.zeros([g.nwire]) + len(g.x)

        # Make sure that normal wires have no delta,
        # and that no-phase wires have an undefined phase
        for w in range(g.nwire):
            if g.w_type[w] != WIRE_TYPE_S:
                g.w_delta[w,:] = 0

        _sc.set_delta(g.x,
                      _n.transpose(g.w_delta),
                      _n.transpose(g.w_phase))
        self.__delta_set = True

    def set_kinetic(self, coefs):
        """
        Set the coefficients for the kinetic equations.

        Parameters
        ----------
        coefs : KineticCoefficients
            Kinetic coefficients for the equations.

        """
        self.__check_serial()

        if not self.__params_set:
            raise RuntimeError, "spectralsolver geometry not set"

        g = self.__geometry
        nx = _n.zeros([g.nwire]) + len(coefs.x)

        # Set them
        _sc.set_kinetic(coefs.E,
                        coefs.x,
                        _n.transpose(coefs.DL[:,:,:]),
                        _n.transpose(coefs.DT[:,:,:]),
                        _n.transpose(coefs.TT[:,:,:]),
                        _n.transpose(coefs.ijE[:,:,:]),
                        _n.transpose(coefs.dDL[:,:,:]),
                        _n.transpose(coefs.dDT[:,:,:]),
                        _n.transpose(coefs.dTT[:,:,:]),
                        _n.transpose(coefs.dijE[:,:,:]),
                        _n.transpose(coefs.cTL[:,:,:]),
                        _n.transpose(coefs.cTT[:,:,:]))
        self.__kinetic_coefs = coefs
        self.__kinetic_set = True

    def get_equilibrium_T(self):
        """
        Get the equilibrium temperature.

        Returns
        -------
        T_eq : float
            Temperature of the system

        Raises
        ------
        RuntimeError
            If the system parameters do not correspond to equilibrium

        """
        return self.equations.get_equilibrium_T()

    def __sp_initialize(self):
        """Initialize spectral solver."""
        self.__check_serial()

        if not self.__params_set:
            raise RuntimeError, "Parameters not set"

        # Set Delta if not set yet
        if not self.__delta_set:
            self.__reset_delta()

        w_type, bc_type, bc_data, bc_connect, nrec = \
                self.equations.get_bc(False)
        #print "w_type:", w_type
        #print "bc_type:", bc_type
        #print "bc_connect:\n", bc_connect
        #print "nrec:", nrec

        self.set_solvers()
        _sc.sp_initialize(w_type, bc_type, _n.ravel(bc_data),
                          _n.ravel(bc_connect))
        self.__sp_initialized = True

    def __kin_initialize(self, E):
        """Initialize kinetic solver."""
        self.__check_serial()

        if not self.__params_set:
            raise RuntimeError, "Parameters not set"
        self.equations = EquationSystem(self.__geometry, E)
        w_type, bc_type, bc_data, bc_connect, nrec = \
                self.equations.get_bc(True)
        self.set_solvers()
        _sc.kin_initialize(w_type, bc_type, _n.ravel(bc_data),
                           _n.ravel(bc_connect), nrec)
        self.__kin_initialized = True

    def __delta_configuration(self, E):
        z = _n.array(self.__geometry.t_delta, copy=True)
        z.sort()
        return _scipy.searchsorted(z, E)

    def sp_solve(self, E, x, continued=False):
        """
        Solve spectral equations.

        Parameters
        ----------
        E : array of floats
            Energies to solve the equations at
        x : array of floats, optional
            Positions to return quantities at. Defaults to the same as those in
            Geometry.
        continued : bool, optional
            Whether to use old solution as a starting point.

        Returns
        -------
        sol : SpResult
            Solution to the equations.

        """
        self.__check_serial()

        nwire = self.__geometry.nwire
        a = _n.zeros([len(E), nwire, len(x)], _n.complex128)
        b = _n.zeros([len(E), nwire, len(x)], _n.complex128)
        da = _n.zeros([len(E), nwire, len(x)], _n.complex128)
        db = _n.zeros([len(E), nwire, len(x)], _n.complex128)

        # solve by continuation from largest energies to smallest

        mink = len(E)
        energies = divide_vector_to_chunks(E, self.__geometry.t_delta)
        for energy in reversed(energies):
            maxk = mink
            mink = maxk - len(energy)

            if self.__delta_configuration(self.__last_E_sp) \
               != self.__delta_configuration(energy[0]):
                continued = False

            if not continued or not self.__sp_initialized:
                self.__sp_initialize()

            # In-place Fortran array assignment
            ca = a[mink:maxk,:,:]
            cb = b[mink:maxk,:,:]
            cda = da[mink:maxk,:,:]
            cdb = db[mink:maxk,:,:]

            _sc.sp_solve(energy, x,
                         _n.transpose(ca),
                         _n.transpose(cb),
                         _n.transpose(cda),
                         _n.transpose(cdb))

            self.__last_E_sp = energy[0]
            continued = True

        return SpResult(E, x, (a, b, da, db))

    def __kin_solve_equilibrium(self, E, x):
        """'Solve' kinetic equations, at equilibrium."""
        self.__check_serial()

        nwire = self.__geometry.nwire
        fT = _n.zeros([len(E), nwire, len(x)], _n.float_)
        fL = _n.zeros([len(E), nwire, len(x)], _n.float_)
        dfT = _n.zeros([len(E), nwire, len(x)], _n.float_)
        dfL = _n.zeros([len(E), nwire, len(x)], _n.float_)
        jL = _n.zeros([len(E), nwire, len(x)], _n.float_)
        jT = _n.zeros([len(E), nwire, len(x)], _n.float_)

        equilibrium_T = self.equations.get_equilibrium_T()

        fL[...] = _n.tanh(.5*E/(equilibrium_T+1e-20))[:,_n.newaxis,_n.newaxis]

        ijE = _scipy.interpolate.interp1d(self.__kinetic_coefs.E,
                                          self.__kinetic_coefs.ijE,
                                          axis=0,
                                          bounds_error=False,
                                          fill_value=0)
        iE = ijE(E)
        if fL.shape != iE.shape:
            iE = _n.swapaxes(iE, 1, 2)
        jT = fL * iE

        return (fL, fT, dfL, dfT, jL, jT)

    def kin_solve(self, E, x, continued=False):
        """
        Solve kinetic equations.

        Parameters
        ----------
        E : array of floats
            Energies to solve the equations at
        x : array of floats, optional
            Positions to return quantities at. Defaults to the same as those in
            Geometry.
        continued : bool, optional
            Whether to use old solution as a starting point.

        Returns
        -------
        sol : KinResult
            Solution to the equations.

        """
        self.__check_serial()

        assert (_n.asarray(E).imag == 0).all(), \
               "Energy should be real when kinetic equations are solved"

        if not self.__kinetic_set:
            raise RuntimeError, "spectralsolver kinetic coefficients not set"

        if not self.equations.get_equilibrium_T() is None:
            return KinResult(E, x, self.__kin_solve_equilibrium(E, x))

        nwire = self.__geometry.nwire
        fT = _n.zeros([len(E), nwire, len(x)], _n.float_)
        fL = _n.zeros([len(E), nwire, len(x)], _n.float_)
        dfT = _n.zeros([len(E), nwire, len(x)], _n.float_)
        dfL = _n.zeros([len(E), nwire, len(x)], _n.float_)
        jL = _n.zeros([len(E), nwire, len(x)], _n.float_)
        jT = _n.zeros([len(E), nwire, len(x)], _n.float_)

        maxk = 0
        energies = divide_vector_to_chunks(E, self.__geometry.t_delta)
        for energy in energies:
            mink = maxk
            maxk = mink + len(energy)

            if self.__delta_configuration(self.__last_E_kin) \
               != self.__delta_configuration(energy[0]):
                continued = False

            if not continued or not self.__kin_initialized:
                self.__kin_initialize(energy[0])

            # In-place Fortran array assignment
            cfL = fL[mink:maxk,:,:]
            cfT = fT[mink:maxk,:,:]
            cdfL = dfL[mink:maxk,:,:]
            cdfT = dfT[mink:maxk,:,:]
            cjL = jL[mink:maxk,:,:]
            cjT = jT[mink:maxk,:,:]

            _sc.kin_solve(E[mink:maxk], x,
                          _n.transpose(cfL),
                          _n.transpose(cfT),
                          _n.transpose(cdfL),
                          _n.transpose(cdfT),
                          _n.transpose(cjL),
                          _n.transpose(cjT))

            self.__last_E_kin = energy[-1]
            continued = True

        return KinResult(E, x, (fL, fT, dfL, dfT, jL, jT))

    def set_solvers(self, sp_solver=None, kin_solver=None,
                    sp_tol=None, kin_tol=None):
        """
        Set solver types and tolerances.

        Parameters
        ----------
        sp_solver : int, optional
            Spectral solver to use.
            Choices are ``SP_SOLVER_COLNEW`` (default) and ``SP_SOLVER_TWPBVP``.
        sp_tol : float, optional
            Tolerance to use in the spectral solver.
            If 0.0, a solver-dependent builtin default value is used.
        kin_solver : int, optional
            Kineticsolver to use.
            Choices are ``KIN_SOLVER_COLNEW`` and ``KIN_SOLVER_TWPBVP`` (default),
            and ``KIN_SOLVER_BLOCK``.
        kin_tol : float, optional
            Tolerance to use in the kinetic solver.
            If 0.0, a solver-dependent builtin default value is used.

        """
        self.__check_serial()

        if sp_solver is not None:
            self.__solver_settings['sp_solver'] = sp_solver
        if kin_solver is not None:
            self.__solver_settings['kin_solver'] = kin_solver
        if sp_tol is not None:
            self.__solver_settings['sp_tol'] = sp_tol
        if kin_tol is not None:
            self.__solver_settings['kin_tol'] = kin_tol

        _sc.set_solvers(self.__solver_settings['sp_solver'],
                        self.__solver_settings['kin_solver'],
                        self.__solver_settings['sp_tol'],
                        self.__solver_settings['kin_tol'])
        self.__sp_initialized = False
        self.__kin_initialized = False

    def get_subgap_bcs(self, energies):
        """Get a list of boundary conditions corresponding to S-interfaces
           below the gap"""
        bcs = []
        for i, E in enumerate(energies):
            eqs = EquationSystem(self.__geometry, E)
            bcs = bcs + [ [bc[0], i, bc[1]] for bc in eqs.get_subgab_bc() ]

        return bcs

class EquationSystem(object):
    """Data for equations at a given energy.

    This class takes care of formatting input data to the form that
    the Fortran code expects.
    """

    def __init__(self, geometry, energy):
        self.geometry = geometry
        self.energy = energy

    def get_equilibrium_T(self):
        """Check whether the system is in equilibrium.
           Returns None if not, otherwise the temperature."""
        g = self.geometry

        equilibrium_T = g.t_t[0]

        for ni in range(g.nnode):
            if (g.t_type[ni] & BCTYPE_MASK) != NODE_CLEAN_NODE:
                if (g.t_t[ni] != equilibrium_T or g.t_mu[ni] != 0
                        or (g.t_type[ni] & BCSUB_MASK) != 0):
                    return None
        return equilibrium_T

    def sanity_check_geometry(self):
        """Check that the geometry is consistent"""
        g = self.geometry

        has_normal = False

        for ni in range(g.nnode):
            if (g.t_type[ni] & BCTYPE_MASK) == NODE_CLEAN_S_TERMINAL:
                if g.t_mu[ni] != 0:
                    raise ValueError, \
                          ("Superconducting terminal %d has " +
                           "a finite bias voltage %g. Non-stationary " +
                           "effects are unsupported.") % ( ni, g.t_mu[ni] )
            elif (g.t_type[ni] & BCTYPE_MASK) == NODE_CLEAN_N_TERMINAL:
                has_normal = True

        if not has_normal and self.get_equilibrium_T() is None:
            raise ValueError, "No normal terminals in the setup, and " \
                  "a non-equilibrium situation was specified. This code " \
                  "does not solve this class of problems: in these cases " \
                  "inelastic relaxation may be important for sub-gap " \
                  "transport. Please specify an equilibrium situation " \
                  "instead."

    def get_subgab_bc(self):
        """Get boundary conditions that are at sub-gap S-interfaces."""
        g = self.geometry
        nodewires, nodetype = self.geometry.get_node_wires()
        subgap_bcs = []
        for i in range(g.nnode):
            if (g.t_type[i] & BCTYPE_MASK) == NODE_CLEAN_S_TERMINAL \
                   and self.energy < g.t_delta[i]:
                for wire in nodewires[i]:
                    subgap_bcs.append([wire, nodetype[i]])
        return subgap_bcs


    def variable_count(self, iskinetic):
        """Count the number of variables in the equation system."""
        g = self.geometry

        mstar = 0
        for k in range(g.nwire):
            if iskinetic:
                mstar = mstar + 4
            else:
                mstar = mstar + 8

        return mstar


    def get_bc(self, iskinetic):
        """Return boundary condition arrays used by the fortran code:
           (w_type, bctype, bcconnect, nrec)"""

        g = self.geometry

        # Collect wires that are connected to each node
        nodewires, wire_end = self.geometry.get_node_wires()

        # Get boundary conditions
        bcs = BoundaryConditionList(self.energy, self.geometry, iskinetic)

        # Boundary conditions must be sorted s.t. those at the right come
        # after those at the left:
        bcs.sort(lambda bc1, bc2: bc1.wire_end - bc2.wire_end)

        # Push the boundary conditions to the correct place
        # and count those at the right
        nrec = 0
        bctype = _n.zeros([len(bcs)], _n.int_)
        bcdata = _n.zeros([16, len(bcs)], _n.float_)
        bcconnect = _n.zeros([self.geometry.nwire+1, len(bcs)], _n.int_)

        for i, bc in enumerate(bcs):
            if bc.wire_end == 1:
                nrec += 1
            bctype[i] = bc.type
            bcdata[:,i] = bc.data
            bcconnect[0:len(bc.connect), i] = _n.array(bc.connect) + 1
            # note the Fortran index offset

        # Check that there are enough boundary conditions
        mstar = self.variable_count(iskinetic)
        if len(bctype) != mstar:
            raise ValueError(("Wrong number of boundary conditions, %d!=%d. "
                              "Your geometry specification probably contains "
                              "errors; please check.") % (len(bctype), mstar))

        #print "nrec=", nrec
        #print "w_type=", g.w_type
        #print "bctype=", bctype
        #print "bcconnect=", bcconnect

        #_DEBUG("bctype =", bctype)

        return (g.w_type, bctype, bcdata, bcconnect, nrec)

class BoundaryCondition(object):
    def __init__(self, type, connect, wire_end, data=None):
        self.type = type
        self.connect = connect
        self.wire_end = wire_end
        self.data = _n.zeros((16,), _n.float_)
        if data is not None:
            self.data[...] = data

class BoundaryConditionList(list):
    """A list of all boundary conditions."""
    def __init__(self, energy, geometry, iskinetic):
        list.__init__(self)
        self.energy = energy
        self.geometry = geometry
        self.iskinetic = iskinetic

        self.nodewires, self.wire_ends = self.geometry.get_node_wires()

        self.node = None
        self.subgap = None

        self.node_sub = 0

        self.__init_boundary_conditions()

    def __init_boundary_conditions(self):
        """Go through all nodes and append the corresponding boundary
           conditions."""

        # A map of all known boundary conditions: (condition, sp, kin)
        bctypes = (
            ( lambda ntype: ntype == NODE_CLEAN_S_TERMINAL,
              self.append_s_terminal_spectral,
              self.append_s_terminal_kinetic ),
            ( lambda ntype: ntype == NODE_CLEAN_S_TERMINAL_CIB,
              self.append_s_terminal_spectral,
              self.append_s_terminal_cib_kinetic),
            ( lambda ntype: ntype == NODE_CLEAN_N_TERMINAL,
              self.append_n_terminal_spectral,
              self.append_n_terminal_kinetic ),
            ( lambda ntype: ntype == NODE_CLEAN_NODE,
              self.append_node_spectral,
              self.append_node_kinetic ),
            ( lambda ntype: ntype == NODE_FREE_INTERFACE,
              self.append_free_interface_spectral,
              self.append_free_interface_kinetic ),
            ( lambda ntype: ntype == NODE_FREE_INTERFACE_FIX_PHASE,
              self.append_free_interface_fix_phase_spectral,
              self.append_free_interface_kinetic ),
            ( lambda ntype: ntype == NODE_TUNNEL_N_TERMINAL,
              self.append_tunnel_n_terminal_spectral,
              self.append_tunnel_n_terminal_kinetic ),
            ( lambda ntype: ntype == NODE_TUNNEL_S_TERMINAL,
              self.append_tunnel_s_terminal_spectral,
              self.append_tunnel_s_terminal_kinetic ),
            ( lambda ntype: ntype == NODE_TUNNEL_NODE,
              self.append_tunnel_node_spectral,
              self.append_tunnel_node_kinetic ),
            )

        for node in range(self.geometry.nnode):
            self.node = node
            self.wires = self.nodewires[node]

            if self.energy < self.geometry.t_delta[self.node]:
                self.subgap = True
            else:
                self.subgap = False

            self.node_sub = (self.geometry.t_type[node] & BCSUB_MASK)

            found = False
            for bctype in bctypes:
                if bctype[0](self.geometry.t_type[node] & BCTYPE_MASK):
                    found = True
                    if self.iskinetic:
                        bctype[2]()
                    else:
                        bctype[1]()
                    break
            if not found:
                raise ValueError("Unknown boundary condition at node %d" % (
                    node))

    ### Utilities
    def append_bcs(self, base, indices, connection, data=None):
        self += [ BoundaryCondition(((base | self.node_sub) << 8) | j,
                                    connection,
                                    self.wire_ends[self.node], data=data)
                 for j in indices ]

    def append_terminal_bcs(self, base, indices, data=None):
        self.append_bcs(base, indices, (self.wires[0], self.node), data=data)

    ### Superconducting terminals
    def append_s_terminal_kinetic(self):
        """Boundary conditions at an S terminal."""
        if self.subgap:
            self.append_terminal_bcs(BCTYPE_CLEAN_S_TERMINAL, (1,3,))
        else:
            self.append_terminal_bcs(BCTYPE_CLEAN_S_TERMINAL, (1,2))

    def append_s_terminal_spectral(self):
        """Boundary conditions at an S terminal."""
        self.append_terminal_bcs(BCTYPE_CLEAN_S_TERMINAL, (1,2,3,4))

    ### Superconducting "terminals", no e-ph relaxation
    def append_s_terminal_cib_kinetic(self):
        """Boundary conditions at an S "terminal", no e-ph relaxation,
           strong charge imbalance."""
        if self.subgap:
            self.append_terminal_bcs(BCTYPE_CLEAN_S_TERMINAL, (1,3,))
        else:
            self.append_terminal_bcs(BCTYPE_FREE_INTERFACE, (1,2))

    ### Normal terminals
    def append_n_terminal_kinetic(self):
        """Boundary conditions at an N terminal."""
        self.append_terminal_bcs(BCTYPE_CLEAN_N_TERMINAL, (1,2))
    def append_n_terminal_spectral(self):
        """Boundary conditions at an N terminal."""
        self.append_terminal_bcs(BCTYPE_CLEAN_N_TERMINAL, (1,2,3,4))

    ### Nodes
    def append_node_kinetic(self):
        # Current conservation
        self.append_bcs(BCTYPE_CLEAN_NODE, (3,4), self.wires)

        # Continuity
        for w in range(1, len(self.wires)):
            self.append_bcs(BCTYPE_CLEAN_NODE, (1,2),
                            (self.wires[w], self.wires[w-1]))
    def append_node_spectral(self):
        # Current conservation
        self.append_bcs(BCTYPE_CLEAN_NODE, (5,6,7,8), self.wires)

        # Continuity
        for i in range(1, len(self.wires)):
            w = self.wires[i]
            self.append_bcs(BCTYPE_CLEAN_NODE, (1,2,3,4),
                            (self.wires[i], self.wires[i-1]))


    ### Free interfaces
    def append_free_interface_kinetic(self):
        self.append_terminal_bcs(BCTYPE_FREE_INTERFACE, (1,2))
    def append_free_interface_spectral(self):
        self.append_terminal_bcs(BCTYPE_FREE_INTERFACE, (1,2,3,4))

    ### Free interfaces (fixed phase)
    def append_free_interface_fix_phase_spectral(self):
        self.append_terminal_bcs(BCTYPE_FREE_INTERFACE_FIX_PHASE, (1,2,3,4))

    ### Tunnel interface (N)

    def _check_resistance(self):
        r = self.geometry.t_resistance[self.node]
        if r < 0:
            raise ValueError("Node %d has invalid resistance %g"
                             % (self.node, r))

    def append_tunnel_n_terminal_kinetic(self):
        self._check_resistance()
        self.append_terminal_bcs(BCTYPE_TUNNEL_N_TERMINAL, (1,2,3,4),
                                 data=self.geometry.t_resistance[self.node])

        raise NotImplementedError("Kinetic K-L boundary conditions are not implemented yet.")

    def append_tunnel_n_terminal_spectral(self):
        self._check_resistance()
        self.append_terminal_bcs(BCTYPE_TUNNEL_N_TERMINAL, (1,2,3,4),
                                 data=self.geometry.t_resistance[self.node])

    ### Tunnel interface (S)
    def append_tunnel_s_terminal_kinetic(self):
        self._check_resistance()
        self.append_terminal_bcs(BCTYPE_TUNNEL_S_TERMINAL, (1,2,3,4),
                                 data=self.geometry.t_resistance[self.node])

        raise NotImplementedError("Kinetic K-L boundary conditions are not implemented yet.")

    def append_tunnel_s_terminal_spectral(self):
        self._check_resistance()
        self.append_terminal_bcs(BCTYPE_TUNNEL_S_TERMINAL, (1,2,3,4),
                                 data=self.geometry.t_resistance[self.node])

    ### Tunnel node
    def append_tunnel_node_kinetic(self):
        self._check_resistance()

        if len(self.wires) != 2:
            raise ValueError("Tunnel junction boundary condition can connect "
                             "only two wires.")

        self.append_bcs(BCTYPE_TUNNEL_NODE,
                        (1,2,3,4),
                        (self.wires[0], self.wires[1]),
                        data=self.geometry.t_resistance[self.node])
        self.append_bcs(BCTYPE_TUNNEL_NODE,
                        (1,2,3,4),
                        (self.wires[1], self.wires[0]),
                        data=self.geometry.t_resistance[self.node])

        raise NotImplementedError("Kinetic K-L boundary conditions are not implemented yet.")

    def append_tunnel_node_spectral(self):
        self._check_resistance()

        if len(self.wires) != 2:
            raise ValueError("Tunnel junction boundary condition can connect "
                             "only two wires.")

        self.append_bcs(BCTYPE_TUNNEL_NODE,
                        (1,2,3,4),
                        (self.wires[0], self.wires[1]),
                        data=self.geometry.t_resistance[self.node])
        self.append_bcs(BCTYPE_TUNNEL_NODE,
                        (1,2,3,4),
                        (self.wires[1], self.wires[0]),
                        data=self.geometry.t_resistance[self.node])


##############################################################################
class SpResult(Sliceable):
    """
    The result from a spectral calculation.

    The parameters `a` and `b` correspond to the representation

    .. math::

       \\hat{G}^R = \\frac{1}{1 - a b} \\begin{pmatrix}
       1 + a b & 2 a \\\\ -2 b & -(1 + a b) \\end{pmatrix}

    of the retarded Green function.

    Attributes
    ----------
    E : array of floats, shape (ne,)
        Energies the solutions are evaluated at
    x : array of floats, shape (nx,)
        Positions the solutions are evaluated at; scaled to [0, 1]
    a : array of floats, shape (ne, nwire, nx)
        Riccati parameter `a`
    b : array of floats, shape (ne, nwire, nx)
        Riccati parameter `b`
    da : array of floats, shape (ne, nwire, nx)
        Derivative of `a`
    db : array of floats, shape (ne, nwire, nx)
        Derivative of `b`
    ne : int
        Number of energy points
    nx : int
        Number of x-points
    nwire : int
        Number of wires

    """

    def __init__(self, E=None, x=None, r=None):
        self.E = E
        self.x = x

        if r:
            self.ne     = _n.shape(r[0])[0]
            self.nwire  = _n.shape(r[0])[1]
            self.nx     = _n.shape(r[0])[2]

            self.a  = r[0]
            self.b  = r[1]
            self.da = r[2]
            self.db = r[3]
        else:
            self.nwire  = None
            self.ne     = None
            self.nx     = None

            self.a = None
            self.b = None
            self.da = None
            self.db = None


        # Indicate how one can slice this
        axis = dict.fromkeys(["a","b","da","db"],
                             (0, 1, 2))
        axis.update({'E': (0, None, None),
                     'x': (None, None, 0)})
        Sliceable.__init__(self, axis)

    def setshape(self, ne, nw, nx):
        self.ne = ne
        self.nwire = nw
        self.nx = nx
        self.E = _n.zeros((ne,), _n.float_)
        self.x = _n.zeros((nx,), _n.float_)
        for name in ["a","b","da","db"]:
            setattr(self, name, _n.zeros((ne, nw, nx), _n.complex128))

class KinResult(Sliceable):
    """
    The result from a kinetic calculation.

    Attributes
    ----------
    E : array of floats, shape (ne,)
        Energies the solutions are evaluated at
    x : array of floats, shape (nx,)
        Positions the solutions are evaluated at; scaled to [0, 1]
    fL : array of floats, shape (ne, nwire, nx)
        The distribution function :math:`f_L`
    fT : array of floats, shape (ne, nwire, nx)
        The distribution function :math:`f_T`
    dfL : array of floats, shape (ne, nwire, nx)
        Derivative of fL
    dfT : array of floats, shape (ne, nwire, nx)
        Derivative of fT
    jL : array of floats, shape (ne, nwire, nx)
        The spectral current :math:`j_L`
    jT : array of floats, shape (ne, nwire, nx)
        The spectral current :math:`j_T`
    ne : int
        Number of energy points
    nx : int
        Number of x-points
    nwire : int
        Number of wires

    """
    def __init__(self, E=None, x=None, r=None):
        self.E = E
        self.x = x

        if r:
            self.ne     = _n.shape(r[0])[0]
            self.nwire  = _n.shape(r[0])[1]
            self.nx     = _n.shape(r[0])[2]

            self.fL  = r[0]
            self.fT  = r[1]
            self.dfL = r[2]
            self.dfT = r[3]
            self.jL  = r[4]
            self.jT  = r[5]

        # Indicate how one can slice this
        axis = dict.fromkeys(["fL", "fT", "dfL", "dfT", "jL", "jT"],
                             (0, 1, 2))
        axis.update({'E': (0, None, None),
                     'x': (None, None, 0)})
        Sliceable.__init__(self, axis)

    def setshape(self, ne, nw, nx):
        self.ne = ne
        self.nwire = nw
        self.nx = nx
        self.E = _n.zeros((ne,), _n.float_)
        self.x = _n.zeros((nx,), _n.float_)
        for name in ["fL", "fT", "dfL", "dfT", "jL", "jT"]:
            setattr(self, name, _n.zeros((ne,nw,nx), _n.float_))



##############################################################################


class KineticCoefficients(Sliceable):
    """
    Coefficients in the kinetic equations.

    Attributes
    ----------
    x : array of floats, shape (nx,)
        x-positions the coefficients are evaluated at; scaled to range [0, 1]
    E : array of floats, shape (ne,)
        Energies the coefficients are evaluated at
    DL : array, shape (ne, nwire, nx)
        Coefficient :math:`D_L`
    DT : array, shape (ne, nwire, nx)
        Coefficient :math:`D_T`
    TT : array, shape (ne, nwire, nx)
        Coefficient :math:`{\\cal T}`
    rjE : array, shape (ne, nwire, nx)
        Coefficient :math:`\\Re j_E`
    ijE : array, shape (ne, nwire, nx)
        Coefficient :math:`\\Im j_E`
    dDL : array, shape (ne, nwire, nx)
        Derivative of DL
    dDT : array, shape (ne, nwire, nx)
        Derivative of DT
    dTT : array, shape (ne, nwire, nx)
        Derivative of TT
    drjE : array, shape (ne, nwire, nx)
        Derivative of rjE
    dijE : array, shape (ne, nwire, nx)
        Derivative of ijE
    cTL : array, shape (ne, nwire, nx)
        prefactor of :math:`f_L` in the sink term for :math:`j_T`
    cTT : array, shape (ne, nwire, nx)
        prefactor of :math:`f_T` in the sink term for :math:`j_T`

    """

    def __init__(self, g=None, r=None, to_evaluate=None):
        """Calculate coefficients in given geometry for given spectral data."""
        self.x = None
        self.E = None
        self.DL = None
        self.DT = None
        self.TT = None
        self.rjE = None
        self.ijE = None
        self.dDL = None
        self.dDT = None
        self.dTT = None
        self.drjE = None
        self.dijE = None
        self.cTL = None
        self.cTT = None

        if to_evaluate is None:
            self.to_evaluate = ["DL", "DT", "TT", "rjE","ijE",
                                "dDL", "dDT", "dTT", "drjE","dijE",
                                "cTL", "cTT"]
        else:
            self.to_evaluate = to_evaluate

        if not g is None and not r is None:
            self.__calculate_from(g, r)

        # Indicate how one can slice this
        axis = dict.fromkeys(["DL", "DT", "TT", "rjE","ijE",
                              "dDL", "dDT", "dTT", "drjE","dijE",
                              "cTL", "cTT"], (0, 1, 2))
        axis.update({'E': (0, None, None),
                     'x': (None, None, 0)})
        Sliceable.__init__(self, axis)

    def setshape(self, ne, nw, nx):
        self.E = _n.zeros((ne,), _n.float_)
        self.x = _n.zeros((nx,), _n.float_)
        for name in ["DL", "DT", "TT", "rjE","ijE",
                     "dDL", "dDT", "dTT", "drjE","dijE",
                     "cTL", "cTT"]:
            setattr(self, name, _n.zeros((ne,nw,nx), _n.float_))

    def __calculate_from(self, g, r):
        a = r.a
        b = r.b
        da = r.da
        db = r.db

        ca = a.conjugate()
        cb = b.conjugate()

        denom = abs(a*b-1)**2
        ddenom = -2 * (da*b + a*db) / (a*b - 1)

        self.x = r.x
        self.E = r.E
        nwire = r.nwire
        ne = len(self.E)

        if "DL" in self.to_evaluate:
            self.DL = ((a*ca - 1)*(b*cb - 1)/denom).real
        if "DT" in self.to_evaluate:
            self.DT = ((a*ca + 1)*(b*cb + 1)/denom).real
        if "TT" in self.to_evaluate:
            self.TT = ((b*cb - a*ca)/denom).real

        if "dDL" in self.to_evaluate:
            self.dDL = (
                + 2*da*ca*(b*cb - 1)/denom
                + 2*db*cb*(a*ca - 1)/denom
                + self.DL * ddenom
                ).real

        if "dDT" in self.to_evaluate:
            self.dDT = (
                + 2*da*ca*(b*cb + 1)/denom
                + 2*db*cb*(a*ca + 1)/denom
                + self.DT * ddenom
                ).real
        if "dTT" in self.to_evaluate:
            self.dTT = (
                + 2*(db*cb - da*ca)/denom
                + self.TT * ddenom
                ).real

        if "rjE" in self.to_evaluate or "ijE" in self.to_evaluate:
            jE = -2j*(a*db - b*da)/(a*b-1)**2
            self.rjE = jE.real
            self.ijE = jE.imag

        if ("drjE" in self.to_evaluate or "dijE" in self.to_evaluate
            or "cTT" in self.to_evaluate or "cTL" in self.to_evaluate):

            delta = g.w_delta[_n.newaxis,:,:]
            phase = g.w_phase[_n.newaxis,:,:]

            djE = 2*delta*(_n.exp(1j*phase)*b - _n.exp(-1j*phase)*a)/(a*b - 1)

            RR = 1j*delta*(
                + _n.exp(-1j*phase)*a*b*cb
                + _n.exp(+1j*phase)*a*ca*b
                + _n.exp(+1j*phase)*b
                + _n.exp(-1j*phase)*a)/denom

            LL = 1j*delta*(
                + _n.exp(-1j*phase)*a*b*cb
                - _n.exp(+1j*phase)*a*ca*b
                + _n.exp(+1j*phase)*b
                - _n.exp(-1j*phase)*a)/denom

            self.drjE = djE.real
            self.dijE = djE.imag

            self.cTT = 2 * RR.real
            self.cTL = 2 * LL.real
