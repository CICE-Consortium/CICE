:tocdepth: 3

.. _adding:

Adding things
=============

.. _addtimer:

~~~~~~
Timers
~~~~~~

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

~~~~~~~~~~~~~~
History fields
~~~~~~~~~~~~~~

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

~~~~~~~
Tracers
~~~~~~~

Each optional tracer has its own module, **ice\_[tracer].F90**, which
also contains as much of the additional tracer code as possible, and for
backward compatibility of binary restart files, each new tracer has its
own binary restart file. We recommend that the logical namelist variable
`tr\_[tracer]` be used for all calls involving the new tracer outside of
**ice\_[tracer].F90**, in case other users do not want to use that
tracer.

A number of optional tracers are available in the code, including ice
age, first-year ice area, melt pond area and volume, brine height,
aerosols, and level ice area and volume (from which ridged ice
quantities are derived). Salinity, enthalpies, age, aerosols, level-ice
volume, brine height and most melt pond quantities are volume-weighted
tracers, while first-year area, pond area, level-ice area and all of the
biogeochemistry tracers in this release are area-weighted tracers. In
the absence of sources and sinks, the total mass of a volume-weighted
tracer such as aerosol (kg) is conserved under transport in horizontal
and thickness space (the mass in a given grid cell will change), whereas
the aerosol concentration (kg/m) is unchanged following the motion, and
in particular, the concentration is unchanged when there is surface or
basal melting. The proper units for a volume-weighted mass tracer in the
tracer array are kg/m.

In several places in the code, tracer computations must be performed on
the conserved “tracer volume" rather than the tracer itself; for
example, the conserved quantity is :math:`h_{pnd}a_{pnd}a_{lvl}a_{i}`,
not :math:`h_{pnd}`. Conserved quantities are thus computed according to
the tracer dependencies, and code must be included to account for new
dependencies (e.g., :math:`a_{lvl}` and :math:`a_{pnd}` in
**ice\_itd.F90** and **ice\_mechred.F90**).

To add a tracer, follow these steps using one of the existing tracers as
a pattern.

#. **ice\_domain\_size.F90**: increase `max\_ntrcr` (can also add option
   to **comp\_ice** and **bld/Macros.\***)

#. **ice\_state.F90**: declare `nt\_[tracer]` and `tr\_[tracer]`

#. **ice\_[tracer].F90**: create initialization, physics, restart
   routines

#. **ice\_fileunits.F90**: add new dump and restart file units

#. **ice\_init.F90**: (some of this may be done in **ice\_[tracer].F90**
   instead)

   -  add new module and `tr\_[tracer]` to list of used modules and
      variables

   -  add logical namelist variable `tr\_[tracer]`

   -  initialize namelist variable

   -  broadcast namelist variable

   -  print namelist variable to diagnostic output file

   -  increment number of tracers in use based on namelist input (`ntrcr`)

   -  define tracer types (`trcr\_depend` = 0 for ice area tracers, 1 for
      ice volume, 2 for snow volume, 2+nt\_[tracer] for dependence on
      other tracers)

#. **ice\_itd.F90**, **ice\_mechred.F90**: Account for new dependencies
   if needed.

#. **CICE\_InitMod.F90**: initialize tracer (includes reading restart
   file)

#. **CICE\_RunMod.F90**, **ice\_step\_mod.F90**:

   -  call routine to write tracer restart data

   -  call physics routines in **ice\_[tracer].F90** (often called from
      **ice\_step\_mod.F90**)

#. **ice\_restart.F90**: define restart variables (for binary,  and PIO)

#. **ice\_history\_[tracer].F90**: add history variables
   (Section :ref:`addhist`)

#. **ice\_in**: add namelist variables to *tracer\_nml* and
   *icefields\_nml*

#. If strict conservation is necessary, add diagnostics as noted for
   topo ponds in Section :ref:`ponds`.
