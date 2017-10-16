.. role:: math(raw)
   :format: html latex
..

Troubleshooting 
================

Check the FAQ: http://oceans11.lanl.gov/drupal/CICE/FAQ.

.. _setup:

Initial setup
-------------

The script **comp\_ice** is configured so that the files **grid**,
**kmt**, **ice\_in**, **run\_ice**, **iced\_gx3\_v5.0** and
**ice.restart\_file** are NOT overwritten after the first setup. If you
wish to make changes to the original files in **input\_templates/**
rather than those in the run directory, either remove the files from the
run directory before executing **comp\_ice** or edit the script.

The code may abort during the setup phase for any number of reasons, and
often the buffer containing the diagnostic output fails to print before
the executable exits. The quickest way to get the diagnostic information
is to run the code in an interactive shell with just the command `cice`
for serial runs or “`mpirun -np N cice`” for MPI runs, where N is the
appropriate number of processors (or a command appropriate for your
computer’s software).

If the code fails to compile or run, or if the model configuration is
changed, try the following:

-  create **Macros.\***. **Makefile.\*** and **run\_ice.\*** files for
   your particular platform, if they do not already exist (type ‘uname
   -s’ at the prompt and compare the result with the file suffixes; we
   rename `UNICOS/mp` as `UNICOS` for simplicity).

-  modify the `INCLUDE` directory path and other settings for your system
   in the scripts, **Macros.\*** and **Makefile.\*** files.

-  alter directory paths, file names and the execution command as needed
   in **run\_ice** and **ice\_in**.

-  ensure that `nprocs` in **ice\_in** is equal to `NTASK` in **comp\_ice**.

-  ensure that the block size `NXBLOCK`, `NYBLOCK` in **comp\_ice** is
   compatible with the processor\_shape and other domain options in
   **ice\_in**

-  if using the rake or space-filling curve algorithms for block
   distribution (`distribution\_type` in **ice\_in**) the code will abort
   if `MXBLCKS` is not large enough. The correct value is provided in the
   diagnostic output.

-  if starting from a restart file, ensure that kcatbound is the same as
   that used to create the file (`kcatbound` = 0 for the files included in
   this code distribution). Other configuration parameters, such as
   `NICELYR`, must also be consistent between runs.

-  for stand-alone runs, check that `-Dcoupled` is *not* set in the
   **Macros.\*** file.

-  for coupled runs, check that `-Dcoupled` and other
   coupled-model-specific (e.g., CESM, popcice or hadgem) preprocessing
   options are set in the **Macros.\*** file.

-  edit the grid size and other parameters in **comp\_ice**.

-  remove the **compile/** directory completely and recompile.

.. _restarttrouble:

Restarts
--------

CICE version 5 introduces a new model configuration that makes
restarting from older simulations difficult. In particular, the number
of ice categories, the category boundaries, and the number of vertical
layers within each category must be the same in the restart file and in
the run restarting from that file. Moreover, significant differences in
the physics, such as the salinity profile, may cause the code to fail
upon restart. Therefore, new model configurations may need to be started
using `runtype` = ‘initial’. Binary restart files that were provided with
CICE v4.1 were made using the BL99 thermodynamics with 4 layers and 5
thickness categories (`kcatbound` = 0) and therefore can not be used for
the default CICE v5 configuration (7 layers). In addition, CICE’s
default restart file format is now  instead of binary.

Restarting a run using `runtype` = ‘continue’ requires restart data for
all tracers used in the new run. If tracer restart data is not
available, use `runtype` = ‘initial’, setting `ice\_ic` to the name of the
core restart file and setting to true the namelist restart flags for
each tracer that is available. The unavailable tracers will be
initialized to their default settings.

On tripole grids, use `restart\_ext` = true when using either binary or
regular (non-PIO) netcdf.

Provided that the same number of ice layers (default: 4) will be used
for the new runs, it is possible to convert v4.1 restart files to the
new file structure and then to  format. If the same physical
parameterizations are used, the code should be able to execute from
these files. However if different physics is used (for instance, mushy
thermo instead of BL99), the code may still fail. To convert a v4.1
restart file:

#. Edit the code **input\_templates/convert\_restarts.f90** for your
   model configuration and path names. Compile and run this code to
   create a binary restart file that can be read using v5. Copy the
   resulting file to the **restart/** subdirectory in your working
   directory.

#. In your working directory, turn off all tracer restart flags in
   **ice\_in** and set the following:

   -  runtype = ‘initial’

   -  ice\_ic = ‘./restart/[your binary file name]’

   -  restart = .true.

   -  use\_restart\_time = .true.

#. In **CICE\_InitMod.F90**, comment out the call to
   restartfile(ice\_ic) and uncomment the call to
   restartfile\_v4(ice\_ic) immediately below it. This will read the
   v4.1 binary file and write a v5  file containing the same
   information.

If restart files are taking a long time to be written serially (i.e.,
not using PIO), see the next section.

Slow execution
--------------

On some architectures, underflows (:math:`10^{-300}` for example) are
not flushed to zero automatically. Usually a compiler flag is available
to do this, but if not, try uncommenting the block of code at the end of
subroutine *stress* in **ice\_dyn\_evp.F90** or **ice\_dyn\_eap.F90**.
You will take a hit for the extra computations, but it will not be as
bad as running with the underflows.

In some configurations, multiple calls to scatter or gather global
variables may overfill MPI’s buffers, causing the code to slow down
(particularly when writing large output files such as restarts). To
remedy this problem, set `BARRIERS yes` in **comp\_ice**. This
synchronizes MPI messages, keeping the buffers in check.

Debugging hints
---------------

Several utilities are available that can be helpful when debugging the
code. Not all of these will work everywhere in the code, due to possible
conflicts in module dependencies.

*debug\_ice* (**CICE.F90**)
    A wrapper for *print\_state* that is easily called from numerous
    points during the timestepping loop (see
    **CICE\_RunMod.F90\_debug**, which can be substituted for
    **CICE\_RunMod.F90**).

*print\_state* (**ice\_diagnostics.F90**)
    Print the ice state and forcing fields for a given grid cell.

`dbug` = true (**ice\_in**)
    Print numerous diagnostic quantities.

`print\_global` (**ice\_in**)
    If true, compute and print numerous global sums for energy and mass
    balance analysis. This option can significantly degrade code
    efficiency.

`print\_points` (**ice\_in**)
    If true, print numerous diagnostic quantities for two grid cells,
    one near the north pole and one in the Weddell Sea. This utility
    also provides the local grid indices and block and processor numbers
    (`ip`, `jp`, `iblkp`, `mtask`) for these points, which can be used in
    conjunction with `check\_step`, to call *print\_state*. These flags
    are set in **ice\_diagnostics.F90**. This option can be fairly slow,
    due to gathering data from processors.

*global\_minval, global\_maxval, global\_sum* (**ice\_global\_reductions.F90**)
    Compute and print the minimum and maximum values for an individual
    real array, or its global sum.

Known bugs
----------

#. Fluxes sent to the CESM coupler may have incorrect values in grid
   cells that change from an ice-free state to having ice during the
   given time step, or vice versa, due to scaling by the ice area. The
   authors of the CESM flux coupler insist on the area scaling so that
   the ice and land models are treated consistently in the coupler (but
   note that the land area does not suddenly become zero in a grid cell,
   as does the ice area).

#. With the old CCSM radiative scheme (`shortwave` = ‘default’ or
   ‘ccsm3’), a sizable fraction (more than 10%) of the total shortwave
   radiation is absorbed at the surface but should be penetrating into
   the ice interior instead. This is due to use of the aggregated,
   effective albedo rather than the bare ice albedo when
   `snowpatch` :math:`< 1`.

#. The date-of-onset diagnostic variables, `melt\_onset` and `frz\_onset`,
   are not included in the core restart file, and therefore may be
   incorrect for the current year if the run is restarted after Jan 1.
   Also, these variables were implemented with the Arctic in mind and
   may be incorrect for the Antarctic.

#. The single-processor *system\_clock* time may give erratic results on
   some architectures.

#. History files that contain time averaged data (`hist\_avg` = true in
   **ice\_in**) will be incorrect if restarting from midway through an
   averaging period.

#. In stand-alone runs, restarts from the end of `ycycle` will not be
   exact.

#. Using the same frequency twice in `histfreq` will have unexpected
   consequences and causes the code to abort.

#. Latitude and longitude fields in the history output may be wrong when
   using padding.

Interpretation of albedos
-------------------------

The snow-and-ice albedo, `albsni`, and diagnostic albedos `albice`, `albsno`,
and `albpnd` are merged over categories but not scaled (divided) by the
total ice area. (This is a change from CICE v4.1 for `albsni`.) The latter
three history variables represent completely bare or completely snow- or
melt-pond-covered ice; that is, they do not take into account the snow
or melt pond fraction (`albsni` does, as does the code itself during
thermodyamic computations). This is to facilitate comparison with
typical values in measurements or other albedo parameterizations. The
melt pond albedo `albpnd` is only computed for the Delta-Eddington
shortwave case.

With the Delta-Eddington parameterization, the albedo depends on the
cosine of the zenith angle (:math:`\cos\varphi`, `coszen`) and is zero if
the sun is below the horizon (:math:`\cos\varphi < 0`). Therefore
time-averaged albedo fields would be low if a diurnal solar cycle is
used, because zero values would be included in the average for half of
each 24-hour period. To rectify this, a separate counter is used for the
averaging that is incremented only when :math:`\cos\varphi > 0`. The
albedos will still be zero in the dark, polar winter hemisphere.

Proliferating subprocess parameterizations
------------------------------------------

With the addition of several alternative parameterizations for sea ice
processes, a number of subprocesses now appear in multiple parts of the
code with differing descriptions. For instance, sea ice porosity and
permeability, along with associated flushing and flooding, are
calculated separately for mushy thermodynamics, topo and level-ice melt
ponds, and for the brine height tracer, each employing its own
equations. Likewise, the BL99 and mushy thermodynamics compute freeboard
and snow–ice formation differently, and the topo and level-ice melt pond
schemes both allow fresh ice to grow atop melt ponds, using slightly
different formulations for Stefan freezing. These various process
parameterizations will be compared and their subprocess descriptions
possibly unified in the future.
