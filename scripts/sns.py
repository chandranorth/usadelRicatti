#!/usr/bin/env python
"""
SNS junction, varying thickness
"""
from __future__ import division
import sys, glob, os

# add build/lib* to sys.path. Not needed if you did setup.py install
# or set PYTHONPATH.
paths = [os.path.join(p, "src")
         for p in glob.glob(os.path.join('build','lib.*'))]
sys.path.extend(paths)
os.environ['LD_LIBRARY_PATH'] = os.pathsep.join(
    paths + [os.environ.get('LD_LIBRARY_PATH', '')])
#

from scipy import *
import usadel1 as u


## Geometry

Delta = 40
omega_D = 600

g = u.Geometry(nwire=3, nnode=4)

g.t_type = [u.NODE_CLEAN_S_TERMINAL,
            u.NODE_CLEAN_S_TERMINAL,
            u.NODE_CLEAN_NODE,
            u.NODE_CLEAN_NODE]

g.t_delta[0:2] = [Delta, Delta]
g.t_phase[0:2] = [-.25*pi, .25*pi]
g.t_inelastic  = 1e-9

g.w_type = [u.WIRE_TYPE_S, u.WIRE_TYPE_N, u.WIRE_TYPE_S]
g.w_length = [0.5, 1, 0.5]
g.w_conductance = [1, 0.8, 1]
g.w_inelastic = 1e-9

g.w_spinflip = [0.4, 0.1, 0.4]

g.w_ends[0,:] = [0, 2]  # wire end to node mapping
g.w_ends[1,:] = [3, 2]
g.w_ends[2,:] = [3, 1]

g.omega_D = omega_D
g.coupling_lambda = 1/arccosh(g.omega_D/Delta)

g.w_phase[0,:] = g.t_phase[0]
g.w_phase[2,:] = g.t_phase[1]
g.w_delta[0,:] = Delta
g.w_delta[2,:] = Delta


## Solve

solver = u.SelfConsistentIteration(g, 'sns.h5', ne=500, maxE=Delta*2,
                                   iterator=u.DummyFixedPointSolver,
				   delta_tolerance=1e-4)
solver.iterate()
