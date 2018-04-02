:tocdepth: 3

.. _running:

Execution procedures
====================

To compile and execute the code: in the source directory,

#. Download the forcing data used for testing from the CICE-Consortium github page,
   https://github.com/CICE-Consortium .

#. Create **Macros.\*** and **run\_ice.\*** files for your particular
   platform, if they do not already exist (type ‘uname -s’ at the prompt
   to get :math:`\langle`\ OS\ :math:`\rangle`).

#. Alter directories in the script **comp\_ice**.

#. Run **comp\_ice** to set up the run directory and make the executable
   ‘**cice**’.

#. | To clean the compile directory and start fresh, simply execute
     ‘/bin/rm -rf compile’ from the run directory.

In the run directory,

#. Alter `atm\_data\_dir` and `ocn\_data\_dir` in the namelist file
   **ice\_in**.

#. Alter the script **run\_ice** for your system.

#. Execute **run\_ice**.

If this fails, see Section :ref:`setup`.

This procedure creates the output log file **ice.log.[ID]**, and if
`npt` is long enough compared with `dumpfreq` and `histfreq`, dump files
**iced.[timeID]** and   (or binary) history output files
**iceh\_[timeID].nc (.da)**. Using the :math:`\left<3^\circ\right>`
grid, the log file should be similar to
**ice.log.\ :math:`\langle`\ OS\ :math:`\rangle`**, provided for the
user’s convenience. These log files were created using MPI on 4
processors on the :math:`\left<3^\circ\right>` grid.

Several options are available in **comp\_ice** for configuring the run,
shown in :ref:`comp-ice`. If `NTASK` = 1, then the **serial/**
code is used, otherwise the code in **mpi/** is used. Loops over blocks
have been threaded throughout the code, so that their work will be
divided among `OMP\_NUM\_THREADS` if `THRD` is ‘yes.’ Note that the value of
`NTASK` in **comp\_ice** must equal the value of `nprocs` in **ice\_in**.
Generally the value of `MXBLCKS` computed by **comp\_ice** is sufficient,
but sometimes it will need to be set explicitly, as discussed in
Section :ref:`performance`. To conserve memory, match the tracer requests
in **comp\_ice** with those in **ice\_in**. CESM uses 3 aerosol tracers;
the number given in **comp\_ice** must be less than or equal to the
maximum allowed in **ice\_domain\_size.F90**.

The scripts define a number of environment variables, mostly as
directories that you will need to edit for your own environment.
`$SYSTEM\_USERDIR`, which on machines at Oak Ridge National Laboratory
points automatically to scratch space, is intended to be a disk where
the run directory resides. `SHRDIR` is a path to the CESM shared code.

:ref:`comp-ice` : Configuration options available in **comp_ice**.

.. _comp-ice:

.. table:: Table 6

   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+
   | variable            | options                              | description                                                                        |
   +=====================+======================================+====================================================================================+
   |RES                  | col, gx3, gx1                        | grid resolution                                                                    |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+
   |NTASK                | (integer)                            | total number of processors                                                         |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+
   |BLCKX                | (integer)                            | number of grid cells on each block in the x-direction :math:`^\dagger`             |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+
   |BLCKY                | (integer)                            | number of grid cells on each block in the y-direction :math:`^\dagger`             |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+
   |MXBLCKS              | (integer)                            | maximum number of blocks per processor                                             |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+
   |NICELYR              | (integer)                            | number of vertical layers in the ice                                               |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+
   |NSNWLYR              | (integer)                            | number of vertical layers in the snow                                              |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+
   |NICECAT              | (integer)                            | number of ice thickness categories                                                 |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+ 
   |TRAGE                | 0 or 1                               | set to 1 for ice age tracer                                                        |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+
   |TRFY                 | 0 or 1                               | set to 1 for first-year ice age tracer                                             |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+
   |TRLVL                | 0 or 1                               | set to 1 for level and deformed ice tracers                                        |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+
   |TRPND                | 0 or 1                               | set to 1 for melt pond tracers                                                     |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+
   |NTRAERO              | 0 or 1                               | number of aerosol tracers                                                          |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+
   |TRBRINE              | set to 1 for brine height tracer     |                                                                                    |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+
   |NBGCLYR              | (integer)                            | number of vertical layers for biogeochemical transport                             |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+
   |IO_TYPE              | none/netcdf/pio                      | use ‘none’ if  library is unavailable,‘pio’ for PIO                                |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+ 
   |DITTO                | yes/no                               | for reproducible diagnostics                                                       |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+
   |BARRIERS             | yes/no                               | flushes MPI buffers during global scatters and gathers                             |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+
   |THRD                 | yes/no                               | set to yes for OpenMP threaded parallelism                                         |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+
   |OMP_NUM_THREADS      | (integer)                            | the number of OpenMP threads requested                                             |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+
   |NUMIN                | (integer)                            | smallest unit number assigned to CICE files                                        |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+
   |NUMAX                | (integer)                            | largest unit number assigned to CICE files                                         |
   +---------------------+--------------------------------------+------------------------------------------------------------------------------------+

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

CICE namelist variables available for changes after compile time appear
in **ice.log.\*** with values read from the file **ice\_in**; their
definitions are given in Section :ref:`index`. For example, to run for a
different length of time, say three days, set `npt` = 72 in **ice\_in**.
At present, the user supplies the time step `dt`, the number of
dynamics/advection/ridging subcycles `ndtd`, and for classic EVP, the
number of EVP subcycles `ndte`; `dte` is then calculated in subroutine
*init\_evp*. The primary reason for doing it this way is to ensure that
`ndte` is an integer. (This is done differently for `revised\_evp` = true.;
see Section :ref:`dynam`).

To restart from a previous run, set restart = true in **ice\_in**. There
are two ways of restarting from a given file. The restart pointer file
**ice.restart\_file** (created by the previous run) contains the name of
the last written data file (**iced.[timeID]**). Alternatively, a
filename can be assigned to ice\_ic in **ice\_in**. Consult
Section :ref:`init` for more details. Restarts are exact for MPI or
single processor runs.

~~~~~~~
Scripts
~~~~~~~

~~~~~~~~~~~
Directories
~~~~~~~~~~~

~~~~~~~~~~~~~~~~~~~
Local modifications
~~~~~~~~~~~~~~~~~~~

~~~~~~~~~~~~
Forcing data
~~~~~~~~~~~~
