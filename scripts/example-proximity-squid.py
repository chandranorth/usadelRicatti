#!/usr/bin/env python
"""
Supercurrent in a proximity squid.

"""
from __future__ import division
from numpy import *
import usadel1 as u

#
# The proximity SQUID
#
#       0
#       |
#       | 0
#       ~
#       1
#      / \
#     /   \
#  1 |     | 2
#     \   /
#      \ /
#       2
#       |
#       | 3
#       |
#       ~
#       3
#

def get_geometry(phi, phi_loop):
    
    # 4 nodes and 4 wires
    g = u.Geometry(4, 4)

    g.t_type = [u.NODE_CLEAN_S_TERMINAL,
                u.NODE_CLEAN_NODE,
                u.NODE_CLEAN_NODE,
                u.NODE_CLEAN_S_TERMINAL]

    # Energy gaps at nodes (has no effect except at terminals)
    g.t_delta = [ 50, 0, 0, 50]
    # Superconducting phases at nodes (has no effect except at terminals)
    g.t_phase = array([ -.5, 0, 0, .5 ]) * phi
    # Inelastic scattering parameter \Gamma
    g.t_inelastic = 1e-9
    # Spin-flip scattering parameter \Gamma_sf
    g.t_spinflip = 0
    # Temperature
    g.t_t = 10
    # Potentials
    g.t_mu = 0

    # Normal-metal structure
    g.w_type = u.WIRE_TYPE_N
    # Wire lengths
    g.w_length = [ 1./3, 1./3, 1./3, 1./3 ]
    # Wire conductance-area products
    g.w_conductance = [ 1, 1, 1, 1 ]

    # Phase jumps at ends of wires 1 and 2, due to applied field
    g.w_phase_jump[1] = phi_loop
    g.w_phase_jump[2] = 0

    # Specify the connections between wires and nodes
    g.w_ends[0,:] = [ 1, 0 ]
    g.w_ends[1,:] = [ 1, 2 ]
    g.w_ends[2,:] = [ 1, 2 ]
    g.w_ends[3,:] = [ 3, 2 ]

    return g

#
# Solve the supercurrent as a function of the two phases
#

def main():
    #---
    phi = pi/2
    phi_loop = 0
    T = 10
    geometry = get_geometry(phi, phi_loop)
    geometry.t_delta = 50
    currents = u.CurrentSolver(geometry, ne=1600,
                               output_function=lambda x: None)
    currents.set_solvers(sp_solver=u.SP_SOLVER_TWPBVP)
    currents.solve_spectral()
    geometry.t_t = T
    Ic, Ie = currents.get_currents(ix=0)
    print "  %14.5g %14.5g %14.5g %14.5g %14.5g" % (
        phi, phi_loop, T, Ic[0], Ic[1] - Ic[2])
    #---
    return
    
    output = open('proximity-squid.dat', 'w')

    print >> output, "%% %14s %14s %14s %14s %14s" % (
        "phi", "phi_loop", "T", "current", "circulating")
    
    for phi in linspace(-pi, pi, 21):
        for phi_loop in linspace(-pi, pi, 21):
            geometry = get_geometry(phi, phi_loop)
            currents = u.CurrentSolver(geometry, ne=160,
                                       output_function=lambda x: None)
            currents.set_solvers(sp_solver=u.SP_SOLVER_TWPBVP)
            currents.solve_spectral()

            print phi, phi_loop, "..."

            ## Uncomment the following line out to dump the solutions
            ## to the spectral & kinetic equations to a HDF5 file.
            ## Use e.g. Matlab and the "h5load.m" script to inspect the output.
            #currents.save('dump.h5')
            
            for T in linspace(1e-6, 20, 100):
                geometry.t_t = T
                Ic, Ie = currents.get_currents(ix=0)

                print >> output, "  %14.5g %14.5g %14.5g %14.5g %14.5g" % (
                    phi, phi_loop, T, Ic[0], Ic[1] - Ic[2])

if __name__ == "__main__":
    main()
