#!/usr/bin/env python
"""
Nonlocal thermovoltage structure

"""
from __future__ import division
from numpy import *
import usadel1 as u
import sys

#
# The structure::
#
#                   N0
#                   |0
#                   |
#              1    ~    9   4
#       .....S1====>6<====9--->N4
#       .           ^
#       .           |6
#       .           |    7
#       .   Phi     7---->11
#       .           |
#       .           |8
#       .      2    ~    10  5
#       .....S2====>8<====10-->N5
#                   ^
#                   |
#                  3|
#                   N3
#
#
# N: normal terminals
# S: superconducting terminals
#
# ... wire omitted from the model
# --- normal wire
# === superconducting wire
#

def get_geometry(Phi):
    
    # 12 nodes and 11 wires
    g = u.Geometry(nnode=12, nwire=11)

    g.t_type = [u.NODE_CLEAN_N_TERMINAL, # 0
                u.NODE_CLEAN_S_TERMINAL, # 1
                u.NODE_CLEAN_S_TERMINAL, # 2
                u.NODE_CLEAN_N_TERMINAL, # 3
                u.NODE_CLEAN_N_TERMINAL, # 4
                u.NODE_CLEAN_N_TERMINAL, # 5
                u.NODE_CLEAN_NODE,       # 6
                u.NODE_CLEAN_NODE,       # 7
                u.NODE_CLEAN_NODE,       # 8
                u.NODE_CLEAN_NODE,       # 9
                u.NODE_CLEAN_NODE,       # 10
                u.NODE_CLEAN_N_TERMINAL, # 11
    ]

    # Energy gaps at nodes (has no effect except at terminals)
    g.t_delta[1] = 100
    g.t_delta[2] = 100
    # Superconducting phases at nodes (has no effect except at terminals)
    g.t_phase[1] = -.5*Phi
    g.t_phase[2] = +.5*Phi
    # Inelastic scattering parameter \Gamma
    g.t_inelastic = 1e-9
    # Spin-flip scattering parameter \Gamma_sf
    g.t_spinflip = 0
    # Temperature
    g.t_t = 5
    # Potentials
    g.t_mu = 0
    
    # Normal-metal structure, except for wires 1, 2, 9 and 10
    g.w_type      = u.WIRE_TYPE_N
    g.w_type[1]   = u.WIRE_TYPE_S
    g.w_type[2]   = u.WIRE_TYPE_S
    g.w_type[9]   = u.WIRE_TYPE_S
    g.w_type[10]  = u.WIRE_TYPE_S
    # Guess a non-self-consistent solution for the order parameter
    g.w_delta[1]  = g.t_delta[1]
    g.w_delta[9]  = g.t_delta[1]
    g.w_delta[2]  = g.t_delta[2]
    g.w_delta[10] = g.t_delta[2]
    g.w_phase[1]  = g.t_phase[1]
    g.w_phase[9]  = g.t_phase[1]
    g.w_phase[2]  = g.t_phase[2]
    g.w_phase[10] = g.t_phase[2]
    # Set parameters for self-consistent iteration
    g.omega_D[1:3]          = 1000
    g.omega_D[9:11]         = 1000
    g.coupling_lambda[1:3]  = 1/arccosh(g.omega_D[9:11]/g.w_delta[9:11])
    g.coupling_lambda[9:11] = 1/arccosh(g.omega_D[9:11]/g.w_delta[9:11])
    # Wire lengths
    g.w_length = 1
    g.w_length[1] = 0.2
    g.w_length[2] = 0.2
    g.w_length[6] = 0.5
    g.w_length[8] = 0.5
    g.w_length[9] = 0.2
    g.w_length[10] = 0.2
    # Wire conductance-area products
    g.w_conductance = 1

    # Specify the connections between wires and nodes
    g.w_ends[0,:]  = [ 0, 6 ]
    g.w_ends[1,:]  = [ 1, 6 ]
    g.w_ends[2,:]  = [ 2, 8 ]
    g.w_ends[3,:]  = [ 3, 8 ]
    g.w_ends[4,:]  = [ 9, 4 ]
    g.w_ends[5,:]  = [10, 5 ]
    g.w_ends[6,:]  = [ 7, 6 ]
    g.w_ends[7,:]  = [ 7, 11]
    g.w_ends[8,:]  = [ 7, 8 ]
    g.w_ends[9,:]  = [ 9, 6 ]
    g.w_ends[10,:] = [10, 8 ]

    return g

#
# Solve the thermoelectric voltage at 4 when terminal 11 is heated
#

def main():
    g = get_geometry(pi/2)

    file_name = 'nonlocal_thermovoltage_spectral.h5'

    # It is possible to do a self-consistent iteration, but this is quite
    # slow.
    #
    # Change False to True below, if you want to try that:
    do_selfconsistent_iteration = False

    # Load data from a file, if it exists
    solver = u.CurrentSolver.resume(file_name, g, ne=200, chunksize=50)
    solver.set_solvers(kin_solver=u.KIN_SOLVER_BLOCK,
                       sp_solver=u.SP_SOLVER_TWPBVP)

    # Solve and save spectral quantities to a file.
    # This takes some time, so we'll want to avoid redoing it unnecessarily.
    #
    # calculage_G=True instructs the solver also to calculate a spectral
    # conductance matrix for the circuit, which can be used to conveniently
    # evaluate currents given distribution functions in the terminals.
    #
    if not do_selfconsistent_iteration:
        solver.solve_spectral_and_save_if_needed(file_name, calculate_G=True)
 

    # Solve thermovoltage vs. temperature
    dT = 0.05

    # ... and write results to a text file
    output = open('nonlocal_thermovoltage.dat', 'w')
    print >> output, "%% %14s %14s" % ("T (E_T)", "dV/dT (ueV/K)")

    if do_selfconsistent_iteration:
        Ts = logspace(log10(dT + 1e-4), log10(10), 15)[1:]
    else:
        Ts = logspace(log10(dT + 1e-4), log10(10), 100)

    for T in Ts:
        
        # (Optional) self-consistent iteration
        if do_selfconsistent_iteration:
            g.t_mu = 0
            g.t_t = T
            it = u.self_consistent_matsubara_iteration(g, max_ne=50)
            #it = u.self_consistent_realtime_iteration(solver)
            for k, d, I_error in it:
                print >> sys.stderr, "%% Self-consistent iteration %d (residual %g)" % (k, d.residual_norm())
                if (d.residual_norm() < 1e-3 * 100 and I_error < 1e-5):
                   break
            else:
                raise RuntimeError("Self-cons. iteration didn't converge!")

            solver.solve_spectral()
            solver.calculate_G()
            solver.save("nonlocal_thermovoltage_T_%.2f.h5" % T)

        # Compute the thermovoltage:
        g.t_mu = 0
        g.t_t = T
        g.t_t[11] += dT

        # Make the terminals 0, 3, 4, 5, 11 to float
        # Currents entering them flow in wires 0, 3, 4, 5, 7
        
        def zero_currents():
            Ic, Ie = solver.get_currents_from_G(w_jT=[0,3,4,5,7], w_jL=[])
            #Ic, Ie = solver.get_currents(w_jT=[0,3,4,5,7], w_jL=[], ix=0)
            return [Ic[0], Ic[3], Ic[4], Ic[5], Ic[7]]

        def set_potentials(z):
            g.t_mu[0], g.t_mu[3], g.t_mu[4], g.t_mu[5], g.t_mu[11] = z

        u.optimize_parameters_for([0,0,0,0,0], zero_currents, set_potentials)

        # Print the thermovoltage at terminal 4 in ueV/K
        print >> output, "  %14g %14g" % (T, g.t_mu[4] / dT * 86.17343)

    # Solve the kinetic equations for some temperature and a larger
    # temperature difference, and dump the result for inspection.
    g.t_t = 1e-4
    g.t_t[11] = 8

    solver.solve_kinetic()
    solver.save('dump.h5')

    output.close()

if __name__ == "__main__":
    main()
