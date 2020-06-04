:tocdepth: 3

.. _adding:

Other things
=============


.. _debugger:

Running with a Debugger
-------------------------

Availability and usage of interactive debuggers varies across machines.  Contact your 
system administrator for additional information about what’s available on your system.  
To run with an interactive debugger, the following general steps should be taken.

- Setup a case
- Modify the env file and Macros file to add appropriate modules and compiler/ linker flags
- Build the model
- Get interactive hardware resources as needed
- Open a csh shell
- Source the env.${machine} file
- Source cice.settings
- Change directories to the run directory
- Manually launch the executable thru the debugger



Reproducible Sums
----------------------

Reproducible sums in the CICE diagnostics are set with the namelist `bfbflag`.
CICE prognostics results do NOT depend on the global sum implementation.  The
results are bit-for-bit identical with any `bfbflag`.  The `bfbflag` only impacts
the results and performance of the global diagnostics written to the CICE
log file.  For best performance, the off (or lsum8 which is equivalent) setting is recommended.
This will probably not produce bit-for-bit results with different decompositions.
For bit-for-bit results, the reprosum setting is recommended.  This should be
only slightly slower than the lsum8 implementation.

Global sums of real types are not reproducible due to different order of operations of the
sums of the individual data which introduced roundoff errors.  
This is caused when the model data is laid out in different
block decompositions or on different pe counts so the data is stored in memory
in different orders.  Integer data should be bit-for-bit identical regardless of
the order of operation of the sums.

The `bfbflag` namelist is a character string with several valid settings.
The tradeoff in these settings is the likelihood for bit-for-bit results versus
their cost.  The `bfbflag` settings are implemented as follows,

off is the default and equivalent to lsum8.

lsum4 is a local sum computed with single precision (4 byte) data and a scalar mpi allreduce.
This is extremely unlikely to be bit-for-bit for different decompositions.
This should generally not be used as the accuracy is very poor for a model
implemented with double precision (8 byte) variables.

lsum8 is a local sum computed with double precision data and a scalar mpi allreduce.
This is extremely unlikely to be bit-for-bit for different decompositions
but is fast.  For CICE implemented in double precision, the differences in global sums
for different decompositions should be at the roundoff level.

lsum16 is a local sum computed with quadruple precision (16 byte) data and a scalar mpi allreduce.
This is very likely to be bit-for-bit for different decompositions.  However,
it should be noted that this implementation is not available or does not work
properly with some compiler and some MPI implementation.  Support for quad precision 
and consistency between underlying fortran and c datatypes can result in inability to
compile or incorrect results.  The source code associated with this implementation
can be turned off with the cpp, NO_R16.  Otherwise, it is recommended that this
option NOT be used or that results be carefully validated on any platform before
it is used.

reprosum is a fixed point method based on ordered double integer sums
that requires two scalar reductions per global sum.  This is extremely likely to be bfb,
but will be slightly more expensive than the lsum algorithms.  See :cite:`Mirin12`

ddpdd is a parallel double-double algorithm using single scalar reduction.
This is very likely to be bfb, but is not as fast or accurate as the reprosum
implementation.  See :cite:`He01`


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
aerosols, water isotopes, and level ice area and volume (from which ridged ice
quantities are derived). Salinity, enthalpies, age, aerosols, isotopes, level-ice
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
`Icepack documentation <https://cice-consortium-icepack.readthedocs.io/en/master/developer_guide/index.html>`__.

To add a tracer, follow these steps using one of the existing tracers as
a pattern.

  1)  **icepack\_tracers.F90** and **icepack\_[tracer].F90**: declare tracers,
      add flags and indices, and create physics routines as described in the
      `Icepack documentation <https://cice-consortium-icepack.readthedocs.io/en/master/developer_guide/dg_adding_tracers.html>`__

  2)  **ice_arrays_column.F90**: declare arrays

  3)  **ice_init_column.F90**: initialize arrays

  4)  **ice\_init.F90**: (some of this may be done in **icepack\_[tracer].F90**
      instead)

     -  declare ``tr_[tracer]``  and ``nt_[tracer]`` as needed

     -  add logical namelist variables ``tr_[tracer]``, ``restart_[tracer]``

     -  initialize and broadcast namelist variables

     -  check for potential conflicts, aborting if any occur

     -  print namelist variables to diagnostic output file

     -  initialize tracer flags etc in icepack (call *icepack_init_tracer_flags* etc)

     -  increment number of tracers in use based on namelist input (``ntrcr``)

     -  define tracer dependencies

  5)  **CICE\_InitMod.F90**: initialize tracer (includes reading restart file)

  6)  **CICE\_RunMod.F90**, **ice\_step\_mod.F90** (and elsewhere as needed):

     -  call routine to write tracer restart data

     -  call Icepack or other routines to update tracer value 
        (often called from **ice\_step\_mod.F90**)

  7)  **ice\_restart.F90**: define restart variables (for binary, netCDF and PIO)

  8)  **ice\_restart\_column.F90**: create routines to read, write tracer restart data

  9)  **ice\_fileunits.F90**: add new dump and restart file units

  10)  **ice\_history\_[tracer].F90**: add history variables
       (Section :ref:`addhist`)

  11)  **ice\_in**: add namelist variables to *tracer\_nml* and
       *icefields\_nml*. Best practice is to set the namelist values so that the 
       new capability is turned off, and create an option file with your preferred
       configuration in **configuration/scripts/options**.

  12)  If strict conservation is necessary, add diagnostics as noted for
       topo ponds in the `Icepack documentation <https://cice-consortium-icepack.readthedocs.io/en/master/science_guide/index.html>`__.

  13)  Update documentation, including **cice_index.rst** and **ug_case_settings.rst**
