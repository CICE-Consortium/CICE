:tocdepth: 3

.. _adding:

Other things
=============


Reproducible Sums
----------------------

The ‘reproducible’ option (`DITTO`) makes diagnostics bit-for-bit when
varying the number of processors. (The simulation results are
bit-for-bit regardless, because they do not require global sums or
max/mins as do the diagnostics.) This was done mainly by increasing the
precision for the global reduction calculations, except for regular
double-precision (r8) calculations involving MPI; MPI can not handle
MPI\_REAL16 on some architectures. Instead, these cases perform sums or
max/min calculations across the global block structure, so that the
results are bit-for-bit as long as the block distribution is the same
(the number of processors can be different).

A more flexible option is available for double-precision MPI
calculations, using the namelist variable `bfbflag`. When true, this flag
produces bit-for-bit identical diagnostics with different tasks,
threads, blocks and grid decompositions.


.. _addtimer:

Adding Timers
-----------------

Timing any section of code, or multiple sections, consists of defining
the timer and then wrapping the code with start and stop commands for
that timer. Printing of the timer output is done simultaneously for all
timers. To add a timer, first declare it (`timer\_[tmr]`) at the top of
**ice\_timers.F90** (we recommend doing this in both the **mpi/** and
**serial/** directories), then add a call to *get\_ice\_timer* in the
subroutine *init\_ice\_timers*. In the module containing the code to be
timed, `call ice\_timer\_start`(`timer\_[tmr]`) at the beginning of the
section to be timed, and a similar call to `ice\_timer\_stop` at the end.
A use `ice\_timers` statement may need to be added to the subroutine being
modified. Be careful not to have one command outside of a loop and the
other command inside. Timers can be run for individual blocks, if
desired, by including the block ID in the timer calls.

.. _addhist:

Adding History fields
-------------------------

To add a variable to be printed in the history output, search for
‘example’ in **ice\_history\_shared.F90**:

#. add a frequency flag for the new field

#. add the flag to the namelist (here and also in **ice\_in**)

#. add an index number

and in **ice\_history.F90**:

#. broadcast the flag

#. add a call to `define\_hist\_field`

#. add a call to `accum\_hist\_field`

The example is for a standard, two-dimensional (horizontal) field; for
other array sizes, choose another history variable with a similar shape
as an example. Some history variables, especially tracers, are grouped
in other files according to their purpose (bgc, melt ponds, etc.).

To add an output frequency for an existing variable, see
section :ref:`history`.

.. _addtrcr:

Adding Tracers
--------------------- 

We require that any changes made to the code be implemented in such a way that they can
be "turned off" through namelist flags.  In most cases, code run with such changes should 
be bit-for-bit identical with the unmodified code.  Occasionally, non-bit-for-bit changes
are necessary, e.g. associated with an unavoidable change in the order of operations. In
these cases, changes should be made in stages to isolate the non-bit-for-bit changes, 
so that those that should be bit-for-bit can be tested separately.

Tracers added to CICE will also require extensive modifications to the Icepack
driver, including initialization, namelist flags 
and restart capabilities.  Modifications to the Icepack driver should reflect
the modifications needed in CICE but are not expected to match completely.
We recommend that the logical namelist variable
``tr_[tracer]`` be used for all calls involving the new tracer outside of
**ice\_[tracer].F90**, in case other users do not want to use that
tracer.

A number of optional tracers are available in the code, including ice
age, first-year ice area, melt pond area and volume, brine height,
aerosols, and level ice area and volume (from which ridged ice
quantities are derived). Salinity, enthalpies, age, aerosols, level-ice
volume, brine height and most melt pond quantities are volume-weighted
tracers, while first-year area, pond area, and level-ice area are area-weighted 
tracers. Biogeochemistry tracers in the skeletal layer are area-weighted,
and vertical biogeochemistry tracers are volume-weighted.  In
the absence of sources and sinks, the total mass of a volume-weighted
tracer such as aerosol (kg) is conserved under transport in horizontal
and thickness space (the mass in a given grid cell will change), whereas
the aerosol concentration (kg/m) is unchanged following the motion, and
in particular, the concentration is unchanged when there is surface or
basal melting. The proper units for a volume-weighted mass tracer in the
tracer array are kg/m.

In several places in the code, tracer computations must be performed on
the conserved "tracer volume" rather than the tracer itself; for
example, the conserved quantity is :math:`h_{pnd}a_{pnd}a_{lvl}a_{i}`,
not :math:`h_{pnd}`. Conserved quantities are thus computed according to
the tracer dependencies (weights), which are tracked using the arrays
``trcr_depend`` (indicates dependency on area, ice volume or snow volume),
``trcr_base`` (a dependency mask), ``n_trcr_strata`` (the number of
underlying tracer layers), and ``nt_strata`` (indices of underlying layers). 
Additional information about tracers can be found in the
`Icepack documentation <https://cice-consortium-icepack.readthedocs.io/en/master/developer_guide/index.html>`_.

To add a tracer, follow these steps using one of the existing tracers as
a pattern.

#. **icepack\_tracers.F90** and **icepack\_[tracer].F90**: declare tracers,
add flags and indices, and create physics routines as described in the
`Icepack documentation <https://cice-consortium-icepack.readthedocs.io/en/master/developer_guide/dg_adding_tracers.html>`_

#. **ice_arrays_column.F90**: declare arrays

#. **ice_init_column.F90**: initialize arrays

#. **ice\_init.F90**: (some of this may be done in **icepack\_[tracer].F90**
   instead)

   -  declare ``tr_[tracer]``  and ``nt_[tracer]`` as needed

   -  add logical namelist variables ``tr_[tracer]``, ``restart_[tracer]``

   -  initialize and broadcast namelist variables

   -  check for potential conflicts, aborting if any occur

   -  print namelist variables to diagnostic output file

   -  initialize tracer flags etc in icepack (call *icepack_init_tracer_flags* etc)

   -  increment number of tracers in use based on namelist input (``ntrcr``)

   -  define tracer dependencies

#.  **CICE\_InitMod.F90**: initialize tracer (includes reading restart file)

#.  **CICE\_RunMod.F90**, **ice\_step\_mod.F90** (and elsewhere as needed):

   -  call routine to write tracer restart data

   -  call Icepack or other routines to update tracer value 
      (often called from **ice\_step\_mod.F90**)

#.  **ice\_restart.F90**: define restart variables (for binary, netCDF and PIO)

#.  **ice\_restart\_column.F90**: create routines to read, write tracer restart data

#.  **ice\_fileunits.F90**: add new dump and restart file units

#.  **ice\_history\_[tracer].F90**: add history variables
   (Section :ref:`addhist`)

#.  **ice\_in**: add namelist variables to *tracer\_nml* and
   *icefields\_nml*. Best practice is to set the namelist values so that the 
   new capability is turned off, and create an option file with your preferred
   configuration in **configuration/scripts/options**.

#.  If strict conservation is necessary, add diagnostics as noted for
   topo ponds in the `Icepack documentation <https://cice-consortium-icepack.readthedocs.io/en/master/science_guide/index.html>`_.

#. Update documentation, including **cice_index.rst** and **ug_case_settings.rst**
