:tocdepth: 3

.. _running_cice:

Running CICE
====================

Quick-start instructions are provided in the :ref:`quickstart` section.

.. _scripts:

Scripts
-------

The CICE scripts are written to allow quick setup of cases and tests.  Once a case is 
generated, users can manually modify the namelist and other files to custom configure
the case.  Several settings are available via scripts as well.

Overview
~~~~~~~~

Most of the scripts that configure, build and run CICE are contained in 
the directory **configuration/scripts/**, except for **cice.setup**, which is
in the main directory.  **cice.setup** is the main script that generates a case. 

Users may need to port the scripts to their local machine.
Specific instructions for porting are provided in :ref:`porting`.

``cice.setup -h`` will provide the latest information about how to use the tool.
``cice.setup --help`` will provide an extended version of the help.
There are three usage modes,

* ``--case`` or ``-c`` creates individual stand alone cases.
* ``--test`` creates individual tests.  Tests are just cases that have some extra automation in order to carry out particular tests such as exact restart.
* ``--suite`` creates a test suite.  Test suites are predefined sets of tests and ``--suite`` provides the ability to quickly setup, build, and run a full suite of tests.

All modes will require use of ``--mach`` or ``-m`` to specify the machine and case and test modes 
can use ``--set`` or ``-s`` to define specific options.  ``--test`` and ``--suite`` will require ``--testid`` to be set 
and both of the test modes can use ``--bdir``, ``--bgen``, ``--bcmp``, and ``--diff`` to generate (save) results and compare results with prior results.
Testing will be described in greater detail in the :ref:`testing` section.

Again, ``cice.setup --help`` will show the latest usage information including 
the available ``--set`` options, the current ported machines, and the test choices.

To create a case, run **cice.setup**::

  cice.setup -c mycase -m machine
  cd mycase

Once a case/test is created, several files are placed in the case directory

- **env.[machine]** defines the environment
- **cice.settings** defines many variables associated with building and running the model
- **makdep.c** is a tool that will automatically generate the make dependencies
- **Macros.[machine]** defines the Makefile macros
- **Makefile** is the makefile used to build the model
- **cice.build** is a script that builds and compiles the model
- **ice\_in** is the namelist input file
- **setup\_run\_dirs.csh** is a script that will create the run directories.  This will be called automatically from the **cice.run** script if the user does not invoke it.
- **cice.run** is a batch run script
- **cice.submit** is a simple script that submits the cice.run script

All scripts and namelist are fully resolved in the case.  Users can edit any
of the files in the case directory manually to change the model configuration,
build options, or batch settings.  The file
dependency is indicated in the above list.  For instance, if any of the files before
**cice.build** in the list are edited, **cice.build** should be rerun.

The **casescripts/** directory holds scripts used to create the case and can 
largely be ignored.  Once a case is created, the **cice.build** script should be run
interactively and then the case should be submitted by executing the 
**cice.submit** script interactively.  The **cice.submit** script
simply submits the **cice.run script**.  
You can also submit the **cice.run** script on the command line.

Some hints:

- To change the tracer numbers or block sizes required at build time, edit the **cice.settings** file.
- To change namelist, manually edit the **ice_in** file
- To change batch settings, manually edit the top of the **cice.run** or **cice.test** (if running a test) file
- To turn on the debug compiler flags, set ``ICE_BLDDEBUG`` in **cice.setttings** to true
- To change compiler options, manually edit the Macros file
- To clean the build before each compile, set ``ICE_CLEANBUILD`` in **cice.settings** to true.  To not clean before the build, set ``ICE_CLEANBUILD`` in **cice.settings** to false

To build and run::

  ./cice.build
  ./cice.submit

The build and run log files will be copied into the logs directory in the case directory.
Other model output will be in the run directory.  The run directory is set in **cice.settings**
via the ``ICE_RUNDIR`` variable.  To modify the case setup, changes should be made in the
case directory, NOT the run directory.

.. _case_options:

Command Line Options
~~~~~~~~~~~~~~~~~~~~

``cice.setup -h`` provides a summary of the command line options.  There are three different modes, ``--case``, ``--test``, and ``--suite``.  This section provides details about the relevant options for setting up cases with examples.
Testing will be described in greater detail in the :ref:`testing` section.

``--help``, ``-h`` 
  prints ``cice.setup`` help information to the terminal and exits.

``--version``
  prints the CICE version to the terminal and exits.

``--setvers VERSION``
  updates the CICE version in your sandbox.  The version should be some like n.m.p.string.

``--case``, ``-c`` CASE
  specifies the case name.  This can be either a relative path of an absolute path.  This cannot be used with --test or --suite.  Either ``--case``, ``--test``, or ``--suite`` is required.

``--mach``, ``-m`` MACHINE
  specifies the machine name.  This should be consistent with the name defined in the Macros and env files in **configurations/scripts/machines**.  This is required in all modes.

``--env``,  ``-e`` ENVIRONMENT1,ENVIRONMENT2,ENVIRONMENT3
  specifies the environment or compiler associated with the machine.  This should be consistent with the name defined in the Macros and env files in **configurations/scripts/machines**.  Each machine can have multiple supported environments including support for different compilers or other system setups.  When used with ``--suite`` or ``--test``, the ENVIRONMENT can be a set of comma deliminated values with no spaces and the tests will then be run for all of those environments.  With ``--case``, only one ENVIRONMENT should be specified. (default is intel)
  
``--pes``,  ``-p`` MxN[[xBXxBY[xMB]
  specifies the number of tasks and threads the case should be run on.  This only works with ``--case``.  The format is tasks x threads or "M"x"N" where M is tasks and N is threads and both are integers. BX, BY, and MB can also be set via this option where BX is the x-direction blocksize, BY is the y-direction blocksize, and MB is the max-blocks setting.  If BX, BY, and MB are not set, they will be computed automatically based on the grid size and the task/thread count.  More specifically, this option has three modes, --pes MxN, --pes MxNxBXxBY, and --pes MxNxBXxBYxMB.  (default is 4x1)

``--acct``  ACCOUNT
  specifies a batch account number.  This is optional.  See :ref:`account` for more information.

``--grid``, ``-g`` GRID
  specifies the grid.  This is a string and for the current CICE driver, gx1 and gx3 are supported. (default = gx3)

``--set``,  ``-s`` SET1,SET2,SET3
  specifies the optional settings for the case.  The settings for ``--suite`` are defined in the suite file.  Multiple settings can be specified by providing a comma deliminated set of values without spaces between settings.  The available settings are in **configurations/scripts/options** and ``cice.setup --help`` will also list them.  These settings files can change either the namelist values or overall case settings (such as the debug flag).

For CICE, when setting up cases, the ``--case`` and ``--mach`` must be specified.  
It's also recommended that ``--env`` be set explicitly as well.  
``--pes`` and ``--grid`` can be very useful.
``--acct`` is not normally used.  A more convenient method 
is to use the **~/cice\_proj** file, see :ref:`account`.  The ``--set`` option can be 
extremely handy.  The ``--set`` options are documented in :ref:`settings`.

.. _settings:

Preset Options
~~~~~~~~~~~~~~

There are several preset options.  These are hardwired in 
**configurations/scripts/options** and are specfied for a case or test by 
the ``--set`` command line option.  You can see the full list of settings 
by doing ``cice.setup --help``.  

The default CICE namelist and CICE settings are specified in the 
files **configuration/scripts/ice_in** and 
**configuration/scripts/cice.settings** respectively.  When picking a 
preset setting (option), the set_env.setting and set_nml.setting will be used to 
change the defaults.  This is done as part of the ``cice.setup`` and the
modifications are resolved in the **cice.settings** and **ice_in** file placed in 
the case directory.  If multiple options are chosen and then conflict, then the last
option chosen takes precedent.  Not all options are compatible with each other.

Some of the options are

``debug`` which turns on the compiler debug flags

``short``, ``medium``, ``long`` which change the batch time limit

``gx3`` and ``gx1`` are associate with grid specific settings

``diag1`` which turns on diagnostics each timestep

``leap`` which turns on the leap year

``pondcesm``, ``pondlvl``, ``pondtopo`` which turn on the various pond schemes

``run10day``, ``run1year``, etc which specifies a run length

``dslenderX1``, ``droundrobin``, ``dspacecurve``, etc specify decomposition options

``dynevp``, ``dyneap``, ``dynoff``, and ``dynrevp`` specify dynamics choices

``therm0``, ``thermBL``, and ``thermmushy`` are thermodynamics settings

``swccsm3`` which turns on the ccsm3 shortwave and albedo computation

``bgc*`` which turns of various bgc configurations

and there are others.  These may change as needed.  Use ``cice.setup --help`` to see the latest.  To add a new option, just add the appropriate file in **configuration/scripts/options**.  For more information, see :ref:`dev_options`

Examples
~~~~~~~~~

The simplest case is just to setup a default configurations specifying the
case name, machine, and environment::

  cice.setup --case mycase1 --mach spirit --env intel

To add some optional settings, one might do::

  cice.setup --case mycase2 --mach spirit --env intel --set debug,diag1,run1year,pondtopo

Once the cases are created, users are free to modify the cice.settings and ice_in namelist to further modify their setup.

.. _porting:

Porting
-------

To port, an **env.[machine]_[environment]** and **Macros.[machine]_[environment}** file have to be added to the
**configuration/scripts/machines/** directory and the 
**configuration/scripts/cice.batch.csh** file needs to be modified.
In general, the machine is specified in ``cice.setup`` with ``--mach``
and the environment (compiler) is specified with ``--env``.
 
- cd to **configuration/scripts/machines/**

- Copy an existing env and a Macros file to new names for your new machine

- Edit your env and Macros files

- cd .. to **configuration/scripts/**

- Edit the **cice.batch.csh** script to add a section for your machine 
  with batch settings and job launch settings

- Download and untar a forcing dataset to the location defined by 
  ``ICE_MACHINE_INPUTDATA`` in the env file

In fact, this process almost certainly will require some iteration.  The easiest way 
to carry this out is to create an initial set of changes as described above, then 
create a case and manually modify the **env.[machine]** file and **Macros.[machine]** 
file until the case can build and run.  Then copy the files from the case 
directory back to **configuration/scripts/machines/** and update 
the **configuration/scripts/cice.batch.csh** file, retest, 
and then add and commit the updated machine files to the repository.

.. _account:

Machine Account Settings
~~~~~~~~~~~~~~~~~~~~~~~~

The machine account default is specified by the variable ``ICE_MACHINE_ACCT`` in 
the **env.[machine]** file.  The easiest way to change a user's default is to 
create a file in your home directory called **.cice\_proj** and add your 
preferred account name to the first line.  
There is also an option (``--acct``) in **cice.setup** to define the account number.  
The order of precedent is **cice.setup** command line option, 
**.cice\_proj** setting, and then value in the **env.[machine]** file.

.. _force:

Forcing data
------------

The input data space is defined on a per machine basis by the ``ICE_MACHINE_INPUTDATA`` 
variable in the **env.[machine]** file.  That file space is often shared among multiple 
users, and it can be desirable to consider using a common file space with group read 
and write permissions such that a set of users can update the inputdata area as 
new datasets are available.

CICE input datasets are stored on an anonymous ftp server.  More information about
how to download the input data can be found at https://github.com/CICE-Consortium/CICE/wiki.
Test forcing datasets are available for various grids at the ftp site.  
These data files are designed only for testing the code, not for use in production runs 
or as observational data. Please do not publish results based on these data sets.


Run Directories
---------------

The **cice.setup** script creates a case directory.  However, the model 
is actually built and run under the ``ICE_OBJDIR`` and ``ICE_RUNDIR`` directories
as defined in the **cice.settings** file.

Build and run logs will be copied from the run directory into the case **logs/** 
directory when complete.


Local modifications
-------------------

Scripts and other case settings can be changed manually in the case directory and
used.  Source code can be modified in the main sandbox.  When changes are made, the code
should be rebuilt before being resubmitted.  It is always recommended that users
modify the scripts and input settings in the case directory, NOT the run directory.
In general, files in the run directory are overwritten by versions in the case
directory when the model is built, submitted, and run.



Old Notes
------------

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

