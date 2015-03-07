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

g.t_type = [u.NODE_CLEAN_S_TERMINAL, u.NODE_CLEAN_S_TERMINAL]

g.t_delta[0:2] = [Delta, Delta]
g.t_phase[0:2] = array([-.5, .5]) * pi/2
g.t_inelastic  = 1e-7

g.w_type = [u.WIRE_TYPE_N]
g.w_length = 1
g.w_conductance = 1

g.w_ends[0,:] = [0, 1] # wire end to node mapping

g.w_phase[0,:] = (g.x - .5) * pi/2

## Solve

solver = u.CurrentSolver.resume('sns_simple.h5', g, ne=300, maxE=Delta*3)
solver.set_solvers(kin_solver=u.KIN_SOLVER_BLOCK)
solver.solve_spectral_and_save_if_needed('sns_simple.h5', calculate_G=False)
