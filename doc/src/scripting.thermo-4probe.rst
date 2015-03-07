.. index:: thermoelectricity

.. _thermoelectricity-in-4-probe-structure:

**************************************
Thermoelectricity in 4-probe structure
**************************************

Consider the simple 4-probe structure:

    .. aafig::

       S0       S1
         \     /
          -----
         /     \
        N2      N3

What is the thermovoltage induced between N2, N3 and S0, S1 when N2
and N3 are at different temperatures and their potentials float?

.. seealso:: `example-thermo-4probe.py <_static/scripts/example-thermo-4probe.py>`__

.. index:: geometry

Geometry specification
======================

Geometry is specified exactly as in the other examples:

.. sourcecode:: python

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

        g.w_type = u.WIRE_TYPE_N
        g.w_length = 1./3
        g.w_conductance = 1

        g.w_ends[0,:]  = [ 4, 0 ]
        g.w_ends[1,:]  = [ 1, 5 ]
        g.w_ends[2,:]  = [ 4, 2 ]
        g.w_ends[3,:]  = [ 3, 5 ]
        g.w_ends[4,:]  = [ 4, 5 ]
        return g

Note that wires connected to a node must either all end or start at
the node.

.. index::
   pair: optimization; thermovoltage

Computing the thermovoltage
===========================

Solve the thermoelectric voltage between N2 and S, when there is a
small temperature difference between N2 and N3.

.. sourcecode:: python

    def main():
        g = get_geometry(pi/2)

        file_name = 'thermo_4probe_spectral.h5'

        solver = u.CurrentSolver.resume(file_name, g)
        solver.set_solvers(kin_solver=u.KIN_SOLVER_BLOCK,
                           sp_solver=u.SP_SOLVER_TWPBVP)

        solver.solve_spectral_and_save_if_needed(file_name, calculate_G=True)

Here, because solving the spectral equations takes some time, we want
to save the result to a file. If the file already exists, solving
the equations is skipped, and old data is reused.

We also instruct the solver to compute the energy-dependent
conductances between all terminals, which allows fast calculation of
currents.

The following is a very straightforward calculation for the
thermovoltage: we set the temperatures of the terminals explicitly,
and adjust the potentials of N2 and N3 until no current enters them:

.. sourcecode:: python

        dT = 0.001

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

The differential thermovoltage could have also been calculated
directly from the energy-dependent conductances.

The result for the S-N voltage at phase difference :math:`\phi=\pi/2`
between the S-terminals of course coincides with the published results: [VT04]_

    .. image:: thermo-4probe.png


.. [VT04]
   P\. Virtanen, and T.T. Heikkilä, Physical Review Letters **92**, 177004 (2004)
