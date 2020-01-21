:tocdepth: 3

.. _running_cice:

Running CICE
====================

Quick-start instructions are provided in the :ref:`quickstart` section.

.. _software:

Software Requirements
----------------------

To run stand-alone, CICE requires

- gmake (GNU Make)
- Fortran and C	compilers (Intel, PGI, GNU, Cray, and NAG have been tested)
- NetCDF
- MPI (this is actually	optional but without it	you can	only run on 1 processor)

Below are lists of software versions that the Consortium has tested at some point.  There is no
guarantee that all compiler versions work with all CICE model versions.  At any given
point, the Consortium is regularly testing on several different compilers, but not 
necessarily on all possible versions or combinations.  A CICE goal is to be relatively portable
across different hardware, compilers, and other software.  As a result, the coding
implementation tends to be on the conservative side at times.  If there are problems 
porting to a particular system, please let the Consortium know.

The Consortium has tested the following compilers at some point,

- Intel 15.0.3.187
- Intel 16.0.1.150
- Intel 17.0.1.132
- Intel 17.0.2.174
- Intel 17.0.5.239
- Intel 18.0.1.163
- Intel 19.0.2
- Intel 19.0.3.199
- PGI 16.10.0
- GNU 6.3.0
- GNU 7.2.0
- GNU 7.3.0
- Cray 8.5.8
- Cray 8.6.4
- NAG 6.2

The Consortium has tested the following mpi versions,

- MPICH 7.3.2
- MPICH 7.5.3
- MPICH 7.6.2
- MPICH 7.6.3
- MPICH 7.7.6
- Intel MPI 18.0.1
- MPT 2.14
- MPT 2.17
- MPT 2.18
- MPT 2.19
- OpenMPI 1.6.5

The NetCDF implementation is relatively general and should work with any version of NetCDF 3 or 4.  The Consortium has tested

- NetCDF 4.3.0
- NetCDF 4.3.2
- NetCDF 4.4.0
- NetCDF 4.4.1.1.32
- NetCDF 4.4.1.1
- NetCDF 4.4.2
- NetCDF 4.5.0
- NetCDF 4.6.1.3

Please email the Consortium if this list can be extended.

.. _scripts:

Scripts
-------

The CICE scripts are written to allow quick setup of cases and tests.  Once a case is 
generated, users can manually modify the namelist and other files to custom configure
the case.  Several settings are available via scripts as well.

.. _overview:

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
and both of the test modes can use ``--bdir``, ``--bgen``, ``--bcmp``, and ``--diff`` to generate (save) results and compare results with prior results as well as ``--tdir`` to specify the location of the test directory.
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
- **cice.build** is a script that calls the Makefile and compiles the model
- **ice\_in** is the namelist input file
- **setup\_run\_dirs.csh** is a script that will create the run directories.  This will be called automatically from the **cice.run** script if the user does not invoke it.
- **cice.run** is a batch run script
- **cice.submit** is a simple script that submits the cice.run script

Once the case is created, all scripts and namelist are fully resolved. Users can edit any
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

- To change the block sizes required at build time, edit the **cice.settings** file.
- To change namelist, manually edit the **ice_in** file
- To change batch settings, manually edit the top of the **cice.run** or **cice.test** (if running a test) file
- To turn on the debug compiler flags, set ``ICE_BLDDEBUG`` in **cice.setttings** to true.  It is also possible to use the ``debug`` option  (``-s debug``) when creating the case with **cice.setup** to set this option automatically.
- To change compiler options, manually edit the Macros file
- To clean the build before each compile, set ``ICE_CLEANBUILD`` in **cice.settings** to true (this is the default value), or use the ``buildclean`` option (``-s buildclean``)  when creating the case with **cice.setup**.  To not clean before the build, set ``ICE_CLEANBUILD`` in **cice.settings** to false, or use the ``buildincremental`` option  (``-s buildincremental``) when creating the case with **cice.setup**.  It is recommended that the ``ICE_CLEANBUILD`` be set to true if there are any questions about whether the build is proceeding properly.

To build and run::

  ./cice.build
  ./cice.submit

The build and run log files will be copied into the logs directory in the case directory.
Other model output will be in the run directory.  The run directory is set in **cice.settings**
via the ``ICE_RUNDIR`` variable.  To modify the case setup, changes should be made in the
case directory, NOT the run directory.

.. _case_options:

**cice.setup** Command Line Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``cice.setup -h`` provides a summary of the command line options.  There are three different modes, ``--case``, ``--test``, and ``--suite``.  This section provides details about the relevant options for setting up cases with examples.
Testing will be described in greater detail in the :ref:`testing` section.

``--help``, ``-h`` 
  prints ``cice.setup`` help information to the terminal and exits.

``--version``
  prints the CICE version to the terminal and exits.

``--setvers VERSION``
  internally updates the CICE version in your sandbox. Those changes can then be commited (or not)
  to the repository. --version will show the updated value. The argument VERSION is typically a
  string like "5.1.2" but could be any alphanumeric string.

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
  specifies the grid.  This is a string and for the current CICE driver, gx1, gx3, and tx1 are supported. (default = gx3)

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
**configuration/scripts/cice.settings** respectively.  When picking
settings (options), the set_env.setting and set_nml.setting will be used to 
change the defaults.  This is done as part of the ``cice.setup`` and the
modifications are resolved in the **cice.settings** and **ice_in** file placed in 
the case directory.  If multiple options are chosen and then conflict, then the last
option chosen takes precedent.  Not all options are compatible with each other.

Some of the options are

``debug`` which turns on the compiler debug flags

``buildclean`` which turns on the option to clean the build before each compile

``buildincremental`` which turns off the option to clean the build before each compile

``short``, ``medium``, ``long`` which change the batch time limit

``gx3``, ``gx1``, ``tx1`` are associate with grid specific settings

``diag1`` which turns on diagnostics each timestep

``run10day``, ``run1year``, etc which specifies a run length

``dslenderX1``, ``droundrobin``, ``dspacecurve``, etc specify decomposition options

``bgcISPOL`` and ``bgcNICE`` specify bgc options

``boxadv``, ``boxdyn``, and ``boxrestore`` are simple box configurations

``alt*`` which turns on various combinations of dynamics and physics options for testing

and there are others.  These may change as needed.  Use ``cice.setup --help`` to see the latest.  
To add a new option, just add the appropriate file in **configuration/scripts/options**.  
For more information, see :ref:`dev_test_options`

Examples
~~~~~~~~~

The simplest case is just to setup a default configuration specifying the
case name, machine, and environment::

  cice.setup --case mycase1 --mach spirit --env intel

To add some optional settings, one might do::

  cice.setup --case mycase2 --mach spirit --env intel --set debug,diag1,run1year

Once the cases are created, users are free to modify the cice.settings and ice_in namelist to further modify their setup.

.. _cicebuild:

More about **cice.build**
~~~~~~~~~~~~~~~~~~~~~~~~~~

**cice.build** is copied into the case directory and should be run interactively from the
case directory to build the model.  CICE is built with make and there is a generic
Makefile and a machine specific Macros file in the case directory.  **cice.build**
is a wrapper for a call to make that includes several other features.  

CICE is built as follows.  First, the makdep binary is created by compiling a small
C program.  The makdep binary is then run and dependency files are created.  The dependency
files are included into the Makefile automatically.  As a result, make dependencies do not 
need to be explicitly defined by the user.  In the next step, make compiles the CICE
code and generates the cice binary.

The standard and recommended way to run is with 
no arguments
::

  cice.build

However, **cice.build** does support a couple other use modes.
::

  cice.build [-h|--help] 

provides a summary of the usage.
::

  cice.build [make arguments] [target]

turns off most of the features of the cice.build script and turns it into a wrapper
for the make call.  The arguments and/or target are passed to make and invoked more
or less like  make [make arguments] [target].  This will be the case if either or 
both the arguments or target are passed to cice.build.  Some examples of that are
::

  cice.build --version

which will pass --version to make.
::

  cice.build targets

is a valid target of the CICE Makefile and simply echos all the valid
targets of the Makefile.
::

  cice.build cice

or ::

  cice.build all

are largely equivalent to running **cice.build** without an argument,
although as noted earlier, many of the extra features of the cice.build script
are turned off when calling cice.build with a target or an argument.  Any of the
full builds will compile makdep, generate the source code dependencies, and
compile the source code.
::

  cice.build [clean|realclean]
  cice.build [db_files|db_flags]
  cice.build [makdep|depends]

are other valid options for cleaning the build, writing out information about
the Makefile setup, and building just the makdep tool or the dependency file.
It is also possible to target a particular CICE object file.

Finally, there is one important parameter in **cice.settings**.  The ``ICE_CLEANBUILD``
variable defines whether the model is cleaned before a build is carried out.  By
default, this variable is true which means each invokation of **cice.build** will
automatically clean the prior build.  If incremental builds are desired to save
time during development, the ``ICE_CLEANBUILD`` setting in **cice.settings** should
be modified.


.. _porting:

Porting
-------

There are four basic issues that need to be addressed when porting, and these are addressed in four separate files in the script system,

- setup of the environment such as compilers, environment variables, and other support software (in **env.[machine]_[environment]**)

- setup of the Macros file to support the model build (in **Macros.[machine]_[environment]**)

- setup of the batch submission scripts (in **cice.batch.csh**)

- setup of the model launch command (in **cice.launch.csh**)

To port, an **env.[machine]_[environment]** and **Macros.[machine]_[environment]** file have to be added to the
**configuration/scripts/machines/** directory and the 
**configuration/scripts/cice.batch.csh** and **configuration/scripts/cice.launch.csh** files need to be modified.
In general, the machine is specified in ``cice.setup`` with ``--mach``
and the environment (compiler) is specified with ``--env``.
 
- cd to **configuration/scripts/machines/**

- Copy an existing env and a Macros file to new names for your new machine

- Edit your env and Macros files

- cd .. to **configuration/scripts/**

- Edit the **cice.batch.csh** script to add a section for your machine 
  with batch settings

- Edit the **cice.batch.csh** script to add a section for your machine 
  with job launch settings

- Download and untar a forcing dataset to the location defined by 
  ``ICE_MACHINE_INPUTDATA`` in the env file

In fact, this process almost certainly will require some iteration.  The easiest way 
to carry this out is to create an initial set of changes as described above, then 
create a case and manually modify the **env.[machine]** file and **Macros.[machine]** 
file until the case can build and run.  Then copy the files from the case 
directory back to **configuration/scripts/machines/** and update 
the **configuration/scripts/cice.batch.csh** and **configuratin/scripts/cice.launch.csh** files, retest, 
and then add and commit the updated machine files to the repository.

.. _machvars: 

Machine variables
~~~~~~~~~~~~~~~~~~~~~~

There are several machine specific variables defined in the **env.$[machine]**.  These
variables are used to generate working cases for a given machine, compiler, and batch
system.  Some variables are optional.

.. csv-table:: *Machine Settings*
   :header: "variable", "format", "description"
   :widths: 15, 15, 25

   "ICE_MACHINE_ENVNAME", "string", "machine name"
   "ICE_MACHINE_COMPILER", "string", "compiler"
   "ICE_MACHINE_MAKE", "string", "make command"
   "ICE_MACHINE_WKDIR", "string", "root work directory"
   "ICE_MACHINE_INPUTDATA", "string", "root input data directory"
   "ICE_MACHINE_BASELINE", "string", "root regression baseline directory"
   "ICE_MACHINE_SUBMIT", "string", "batch job submission command"
   "ICE_MACHINE_TPNODE", "integer", "machine maximum MPI tasks per node"
   "ICE_MACHINE_MAXPES", "integer", "machine maximum total processors per job (optional)"
   "ICE_MACHINE_MAXTHREADS", "integer", "machine maximum threads per mpi task (optional)"
   "ICE_MACHINE_MAXRUNLENGTH", "integer", "batch wall time limit in hours (optional)"
   "ICE_MACHINE_ACCT", "string", "batch default account"
   "ICE_MACHINE_QUEUE", "string", "batch default queue"
   "ICE_MACHINE_BLDTHRDS", "integer", "number of threads used during build"
   "ICE_MACHINE_QSTAT", "string", "batch job status command (optional)"
   "ICE_MACHINE_QUIETMODE", "true/false", "flag to reduce build output (optional)"

.. _cross_compiling:

Cross-compiling
~~~~~~~~~~~~~~~
It can happen that the model must be built on a platform and run on another, for example when the run environment is only available in a batch queue. The program **makdep** (see :ref:`overview`), however, is both compiled and run as part of the build process.

In order to support this, the Makefile uses a variable ``CFLAGS_HOST`` that can hold compiler flags specfic to the build machine for the compilation of makdep. If this feature is needed, add the variable ``CFLAGS_HOST`` to the **Macros.[machine]_[environment]** file. For example : ::

  CFLAGS_HOST = -xHost

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

.. _queue:

Machine Queue Settings
~~~~~~~~~~~~~~~~~~~~~~~~

Supported machines will have a default queue specified by the variable ``ICE_MACHINE_QUEUE``
in the **env.[machine]** file.  This can also be manually changed in the **cice.run** or
**cice.test** scripts or even better, use the ``--queue`` option in **cice.setup**.

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

.. _timeseries:

Timeseries Plotting
-------------------

The CICE scripts include two scripts that will generate timeseries figures from a 
diagnostic output file, a Python version (``timeseries.py``) and a csh version 
(``timeseries.csh``).  Both scripts create the same set of plots, but the Python 
script has more capabilities, and it's likely that the csh
script will be removed in the future.  

To use the ``timeseries.py`` script, the following requirements must be met:

* Python v2.7 or later
* numpy Python package
* matplotlib Python package
* datetime Python package

See :ref:`CodeCompliance` for additional information about how to setup the Python 
environment, but we recommend using ``pip`` as follows: ::

  pip install --user numpy
  pip install --user matplotlib
  pip install --user datetime

When creating a case or test via ``cice.setup``, the ``timeseries.csh`` and 
``timeseries.py`` scripts are automatically copied to the case directory.  
Alternatively, the plotting scripts can be found in ``./configuration/scripts``, and can be
run from any directory.

The Python script can be passed a directory, a specific log file, or no directory at all:

  - If a directory is passed, the script will look either in that directory or in 
    directory/logs for a filename like cice.run*.  As such, users can point the script
    to either a case directory or the ``logs`` directory directly.  The script will use 
    the file with the most recent creation time.
  - If a specific file is passed the script parses that file, assuming that the file
    matches the same form of cice.run* files.
  - If nothing is passed, the script will look for log files or a ``logs`` directory in the 
    directory from where the script was run.

For example:

Run the timeseries script on the desired case. ::

$ python timeseries.py /p/work1/turner/CICE_RUNS/conrad_intel_smoke_col_1x1_diag1_run1year.t00/

or ::

$ python timeseries.py /p/work1/turner/CICE_RUNS/conrad_intel_smoke_col_1x1_diag1_run1year.t00/logs
    
The output figures are placed in the directory where the ``timeseries.py`` script is run.

The plotting script will plot the following variables by default, but you can also select 
specific plots to create via the optional command line arguments.

  - total ice area (:math:`km^2`)
  - total ice extent (:math:`km^2`)
  - total ice volume (:math:`m^3`)
  - total snow volume (:math:`m^3`)
  - RMS ice speed (:math:`m/s`)

For example, to plot only total ice volume and total snow volume ::

$ python timeseries.py /p/work1/turner/CICE_RUNS/conrad_intel_smoke_col_1x1_diag1_run1year.t00/ --volume --snw_vol

To generate plots for all of the cases within a suite with a testid, create and run a script such as  ::

     #!/bin/csh
     foreach dir (`ls -1  | grep testid`)
       echo $dir
       python timeseries.py $dir
     end

Plots are only made for a single output file at a time.  The ability to plot output from 
a series of cice.run* files is not currently possible, but may be added in the future.
However, using the ``--bdir`` option will plot two datasets (from log files) on the
same figure.

For the latest help information for the script, run ::

$ python timeseries.py -h

The ``timeseries.csh`` script works basically the same way as the Python version, however it
does not include all of the capabilities present in the Python version.  

To use the C-Shell version of the script, ::

$ ./timeseries.csh /p/work1/turner/CICE_RUNS/conrad_intel_smoke_col_1x1_diag1_run1year.t00/

