#!/usr/bin/env python
"""
Nonlocal thermovoltage structure

"""
from __future__ import division
from numpy import *
import usadel1 as u

#
# The structure::
#
#   S0      S1
#     \    /
#      4--5
#     /    \
#    /      \
#   N2       N3
#

def get_geometry(phi):
    
    # 6 nodes and 5 wires
    g = u.Geometry(nnode=6, nwire=5)

    g.t_type = [u.NODE_CLEAN_S_TERMINAL, # 0
                u.NODE_CLEAN_S_TERMINAL, # 1
                u.NODE_CLEAN_N_TERMINAL, # 2
                u.NODE_CLEAN_N_TERMINAL, # 3
                u.NODE_CLEAN_NODE,       # 4
                u.NODE_CLEAN_NODE,       # 5
    ]

    g.t_delta[0:2] = 100
    g.t_phase[0] = -.5*phi
    g.t_phase[1] = +.5*phi
    g.t_inelastic = 1e-9
    
    g.w_type      = u.WIRE_TYPE_N
    g.w_length = 1./3
    g.w_conductance = 1

    g.w_ends[0,:]  = [ 4, 0 ]
    g.w_ends[1,:]  = [ 1, 5 ]
    g.w_ends[2,:]  = [ 4, 2 ]
    g.w_ends[3,:]  = [ 3, 5 ]
    g.w_ends[4,:]  = [ 4, 5 ]
    return g

#
# Solve the thermoelectric voltage between N and S
#

def main():
    g = get_geometry(pi/2)

    file_name = 'thermo_4probe_spectral.h5'

    # Load data from a file, if it exists
    solver = u.CurrentSolver.resume(file_name, g)
    solver.set_solvers(kin_solver=u.KIN_SOLVER_BLOCK,
                       sp_solver=u.SP_SOLVER_TWPBVP)

    # Solve and save spectral quantities to a file.
    # This takes some time, so we'll want to avoid redoing it unnecessarily.
    #
    # calculage_G=True instructs the solver also to calculate a spectral
    # conductance matrix for the circuit, which can be used to conveniently
    # evaluate currents given distribution functions in the terminals.
    #
    solver.solve_spectral_and_save_if_needed(file_name, calculate_G=True)

    # Solve thermovoltage vs. temperature
    dT = 0.001

    # ... and write results to a text file
    output = open('thermo_4probe.dat', 'w')
    print >> output, "%% %14s %14s" % ("T (E_T)", "dV/dT (ueV/K)")

    for T in logspace(log10(dT + 1e-4), log10(20), 100):
        # Compute the thermovoltage:
        g.t_mu = 0
        g.t_t = T
        g.t_t[2] += dT/2
        g.t_t[3] -= dT/2

        # Make the terminals 2, 3 to float
        # Currents entering them flow in wires 2, 3
        
        def zero_currents():
            Ic, Ie = solver.get_currents_from_G(w_jT=[2,3], w_jL=[])
            return [Ic[2], Ic[3]]

        def set_potentials(z):
            g.t_mu[2], g.t_mu[3] = z

        u.optimize_parameters_for([0,0], zero_currents, set_potentials)

        # Print the thermovoltage at terminal 2 in ueV/K
        print >> output, "  %14g %14g" % (T, g.t_mu[2] / dT * 86.17343)

    # NOTE: This is the simplest way to calculate the thermovoltage;
    #       it would also be possible to calculate directly the linear-response
    #       voltage, as the G-matrix is known.

if __name__ == "__main__":
    main()
