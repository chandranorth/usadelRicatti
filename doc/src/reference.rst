Reference
=========

.. automodule:: usadel1
   :members:


.. index:: geometry

Geometry
--------
.. autoclass:: Geometry
   :members:


Solver
------
.. autoclass:: Solver

   .. automethod:: set_geometry

   .. automethod:: set_kinetic

   .. automethod:: sp_solve

   .. automethod:: kin_solve


Results
-------
.. autoclass:: SpResult

.. autoclass:: KinResult

.. autoclass:: KineticCoefficients


Currents
--------

.. autoclass:: CurrentSolver

   .. rubric:: Solving equations
   
   .. automethod:: solve

   .. automethod:: solve_spectral

   .. automethod:: solve_kinetic

   .. automethod:: solve_spectral_if_needed

   .. automethod:: solve_spectral_and_save_if_needed

   .. automethod:: calculate_G
   
   .. automethod:: approximate_G

   .. automethod:: set_solvers

   .. rubric:: Computing currents

   .. automethod:: get_currents

   .. automethod:: get_currents_lazy

   .. automethod:: get_currents_from_G

   .. automethod:: get_currents_from_G_with_f

   .. rubric:: Computing linear response

   .. automethod:: get_linear_response_from_G

   .. automethod:: get_linear_response_from_G_with_f

   .. rubric:: Saving and loading data:

   .. automethod:: load

   .. automethod:: save

   .. automethod:: resume


.. index:: optimization; parameters

Optimizing parameters
---------------------

.. autofunction:: optimize_parameters_for

.. autoexception:: AlreadyConvergedException


.. index:: self-consistency

Self-consistent iteration
-------------------------
.. currentmodule:: usadel1

.. autofunction:: self_consistent_realtime_iteration

.. autofunction:: self_consistent_matsubara_iteration


.. index:: HDF5, data file

HDF5 file layout
----------------

The HDF5 files produced by :meth:`CurrentSolver.save` have the layout::

    /geometry/w_type
    /geometry/t_type
    ...
    /spectral/a
    /spectral/b
    ...
    /kinetic/fL
    /kinetic/fT
    ...
    /coefficient/DL
    /coefficient/DT
    ...

i.e. the attributes of the :class:`Geometry`, :class:`SpResult`, 
:class:`KinResult`, and :class:`KineticCoefficient` objects are simply dumped
to the file to an appropriate place.

You can easily load the HDF5 files by using the supplied ``h5load.m`` script
(probably requires a recent version of Matlab).

.. seealso:: `h5load.m <_static/scripts/h5load.m>`__

