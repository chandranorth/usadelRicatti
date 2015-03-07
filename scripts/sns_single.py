#!/usr/bin/env python
"""
SNS wire, constant thickness
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

Delta = 40    # Energies are in units of the Thouless energy correspondign
              # to the length 1

omega_D = 600 # Debye energy
LS = 0.9      # Portion of the superconducting part

x = linspace(0, 1, 300)
g = u.Geometry(nwire=1, nnode=2, x=x)

g.t_type = [u.NODE_CLEAN_S_TERMINAL,
            u.NODE_CLEAN_S_TERMINAL]

g.t_delta[0:2] = [Delta, Delta]
g.t_phase[0:2] = array([-.5, .5]) * pi/2
g.t_inelastic  = 1e-9

g.w_type = [u.WIRE_TYPE_S]
g.w_length = 1
g.w_conductance = 1
g.w_inelastic = 1e-9

g.w_ends[0,:] = [0, 1] # wire end to node mapping

g.omega_D = omega_D
g.coupling_lambda = 1/arccosh(g.omega_D/Delta)

## Split the wire to superconducting and normal parts:

x1 = LS / 2
x2 = 1 - LS / 2

g.coupling_lambda = 1/arccosh(g.omega_D/Delta)
g.coupling_lambda[0, (g.x > x1) & (g.x < x2)] = 0

g.w_phase[0, g.x < x1] = g.t_phase[0]
g.w_phase[0, g.x > x2] = g.t_phase[1]
g.w_delta[0, (g.x < x1) | (g.x > x2)] = Delta

## Solve

solver = u.SelfConsistentIteration(g, 'sns_single.h5', ne=300, maxE=Delta*3)
solver.iterate()
