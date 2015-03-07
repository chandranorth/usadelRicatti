"""Sample geometries for testing."""
from __future__ import division
from usadel1 import *
from numpy import *

__revision__ = "$Id: geometries.py 3260 2007-01-23 13:14:22Z pauli $"

def geometry_NnN(V0=0, T0=1e-4, phi=0):
    """Simple normal wire between normal terminals."""
    g = Geometry(1, 2)

    g.t_type = [ NODE_CLEAN_N_TERMINAL, NODE_CLEAN_N_TERMINAL ]

    g.t_delta = 0
    g.t_phase = 0

    g.t_inelastic = 0
    g.t_spinflip = 0
    g.t_t = [ T0, T0 ]
    g.t_mu = [ V0, 0 ]

    g.w_type = [ WIRE_TYPE_N ]
    g.w_length = 1
    g.w_conductance = 1
    g.w_inelastic = 0
    g.w_spinflip = 0

    g.w_ends[0,:] = [ 0, 1 ]

    return g

def geometry_SnnNnS(V0=0, T0=1e-4, phi=0):
    """The three-probe setup."""
    g = Geometry(3, 4)

    g.t_type = [ NODE_CLEAN_S_TERMINAL,
                 NODE_CLEAN_S_TERMINAL,
                 NODE_CLEAN_N_TERMINAL,
                 NODE_CLEAN_NODE ]

    g.t_delta = array([ 5000, 5000, 0, 0 ])
    g.t_phase = array([ -1, 1, 0, 0 ]) * (pi/2) * 0.5 * phi

    g.t_inelastic = 1e-9
    g.t_spinflip = 0
    g.t_t = T0
    g.t_mu = [ 0, 0, V0, 0 ]

    g.w_type = [ WIRE_TYPE_N, WIRE_TYPE_N, WIRE_TYPE_N ]
    g.w_length = [ 0.5, 0.5, 5 ]
    g.w_conductance = 1
    g.w_inelastic = 1e-9
    g.w_spinflip = 0

    g.w_ends[0,:] = [ 0, 3 ]
    g.w_ends[1,:] = [ 1, 3 ]
    g.w_ends[2,:] = [ 2, 3 ]

    return g


def geometry_SnnNnNnnS(T1=1e-4, T2=1e-4, T0=1e-4, phi=0):
    """The 4-probe 5-wire setup."""
    g = Geometry(5, 6)

    g.t_type = [ NODE_CLEAN_N_TERMINAL,
                 NODE_CLEAN_N_TERMINAL,
                 NODE_CLEAN_S_TERMINAL, 
                 NODE_CLEAN_S_TERMINAL,
                 NODE_CLEAN_NODE,
                 NODE_CLEAN_NODE ]

    g.t_delta = array([ 0, 0, 5000, 5000, 0, 0 ])
    g.t_phase = array([ 0, 0, -1, 1, 0, 0 ]) * (pi/2) * 0.5 * phi

    g.t_inelastic = 1e-8
    g.t_spinflip = 0
    g.t_t = [ T1, T2, T0, T0, T0, T0 ]
    g.t_mu = [ 0, 0, 0, 0, 0, 0 ]

    g.w_type = [ WIRE_TYPE_N, WIRE_TYPE_N, WIRE_TYPE_N,
                 WIRE_TYPE_N, WIRE_TYPE_N ]
    g.w_length = array([ 1, 1, 1, 1, 1  ]) / 3
    g.w_conductance = 1
    g.w_inelastic = 1e-9
    g.w_spinflip = 0

    g.w_ends[0,:] = [ 4, 0 ]
    g.w_ends[1,:] = [ 1, 5 ]
    g.w_ends[2,:] = [ 4, 2 ]
    g.w_ends[3,:] = [ 3, 5 ]
    g.w_ends[4,:] = [ 4, 5 ]

    return g

def geometry_SNS(x, V0, T0, phi, Delta):
    """N-wire between S terminals."""
    g = Geometry(1, 2, x)
    
    g.t_type = [ NODE_CLEAN_S_TERMINAL,
                 NODE_CLEAN_S_TERMINAL ]
    
    g.t_delta = [ Delta, Delta ]
    g.t_phase = [ -.5*phi, .5*phi ]
    
    g.t_inelastic = 1e-9
    g.t_spinflip = 0
    g.t_t = T0
    g.t_mu = 0
    
    g.w_type = WIRE_TYPE_N
    g.w_length = 1
    g.w_conductance = 1
    g.w_inelastic = 1e-9
    g.w_spinflip = 0
    
    g.w_ends[0,:] = [ 0, 1 ]
    
    return g

def geometry_SNS_split(x, V0, T0, phi, Delta):
    """N-wire between S terminals, split into two parts."""
    g = Geometry(2, 3, x)
    
    g.t_type = [ NODE_CLEAN_S_TERMINAL,
                 NODE_CLEAN_NODE,
                 NODE_CLEAN_S_TERMINAL ]
    
    g.t_delta = [ Delta, 0, Delta ]
    g.t_phase = [ -.5*phi, 0, .5*phi ]
    
    g.t_inelastic = 1e-9
    g.t_spinflip = 0
    g.t_t = T0
    g.t_mu = 0
    
    g.w_type = WIRE_TYPE_N
    g.w_length = [ 1.0/3.0, 2.0/3.0 ]
    g.w_conductance = 1
    g.w_inelastic = 1e-9
    g.w_spinflip = 0
    
    g.w_ends[0,:] = [ 0, 1 ]
    g.w_ends[1,:] = [ 2, 1 ]
    
    return g

