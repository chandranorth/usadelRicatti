#!/usr/bin/env python
"""
Print minigap E_g(phi) in an SNS junction

:author: Pauli Virtanen <pauli@ltl.tkk.fi>
"""
from __future__ import division
from scipy import *
import usadel1 as u
import traceback

__revision__ = "$Id: sns_minigap_vs_phi.py 3275 2007-02-14 11:00:03Z pauli $"

def print_minigap_vs_phi():
    seterr(invalid='ignore')

    solver = u.Solver()
    Delta = 1e6

    f = open('sns_minigap_vs_phi.dat', 'a')

    print >> f, "% phi   E_g"
    for phi in linspace(0, pi, 20):
        g = get_geometry(Delta, phi)
        solver.set_geometry(g)
        solver.set_solvers(sp_tol=1e-6)

        try:
            E = linspace(0.5*3*cos(phi/2), 1.5*3*cos(phi/2), 1500)
            solution = solver.sp_solve(E, g.x)
            gap = all(((1+solution.a*solution.b)/(1-solution.a*solution.b))[:,0,:].real < 1e-5, axis=1)
            if gap.any():
                E_g = max(E[gap])
            else:
                E_g = 0
            print phi, E_g
            print >> f, phi, E_g
        except:
            traceback.print_exc()
            raise

    f.close()

def get_geometry(Delta, phi):
    g = u.Geometry(1, 2)
    
    g.t_type = [u.NODE_CLEAN_S_TERMINAL, u.NODE_CLEAN_S_TERMINAL]
     
    g.t_delta = [Delta, Delta]
    g.t_phase = [-phi/2, phi/2]
    g.t_inelastic = 1e-9
    
    g.w_type = [ u.WIRE_TYPE_N ]
    g.w_length = 1
    g.w_conductance = 1
    
    g.w_ends[0,:] = [ 0, 1 ]
    
    return g

if __name__ == "__main__":
    print_minigap_vs_phi()
