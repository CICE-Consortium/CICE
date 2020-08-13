:tocdepth: 3

.. _running_cice:

Running CICE
====================

Quick-start instructions are provided in the :ref:`quickstart` section.

.. _software:

Software Requirements
----------------------

To run stand-alone, CICE requires

- bash and csh
- gmake (GNU Make)
- Fortran and C	compilers (Intel, PGI, GNU, Cray, and NAG have been tested)
- NetCDF (this is actually optional but required to test out of the box configurations)
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

All modes will require use of ``--mach`` or ``-m`` to specify the machine.  Use of ``--env`` is also recommended to specify the compilation environment.  ``--case`` and ``--test`` modes can use ``--set`` or ``-s`` which will turn on various model options.  ``--test`` and ``--suite`` will require ``--testid`` to be set and can use ``--bdir``, ``--bgen``, ``--bcmp``, and ``--diff`` to generate (save) results for regression testing (comparison with prior results). ``--tdir`` will specify the location of the test directory.
Testing will be described in greater detail in the :ref:`testing` section.

Again, ``cice.setup --help`` will show the latest usage information including 
the available ``--set`` options, the current ported machines, and the test choices.

To create a case, run **cice.setup**::

  cice.setup -c mycase -m machine -e intel
  cd mycase

Once a case/test is created, several files are placed in the case directory

- **env.[machine]_[env]** defines the environment
- **cice.settings** defines many variables associated with building and running the model
- **makdep.c** is a tool that will automatically generate the make dependencies
- **Macros.[machine]_[env]** defines the Makefile macros
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
submits the **cice.run script** or **cice.test** script.  These scripts can
also be run interactively or submitted manually without the **cice.submit** script.

Some hints:

- To change namelist, manually edit the **ice_in** file
- To change batch settings, manually edit the top of the **cice.run** or **cice.test** (if running a test) file
- When the run scripts are submitted, the current **ice_in**, **cice.settings**, and **env.[machine]** files are copied from the case directory into the run directory.  Users should generally not edit files in the run directory as these are overwritten when following the standard workflow.  **cice.settings** can be sourced to establish the case values in the login shell.  An alias like the following can be established to quickly switch between case and run directories::

    alias  cdrun 'cd `\grep "setenv ICE_RUNDIR"  cice.settings | awk "{print "\$"NF}"`'
    alias cdcase 'cd `\grep "setenv ICE_CASEDIR" cice.settings | awk "{print "\$"NF}"`'

- To turn on the debug compiler flags, set ``ICE_BLDDEBUG`` in **cice.setttings** to true.  It is also possible to use the ``debug`` option  (``-s debug``) when creating the case with **cice.setup** to set this option automatically.
- To change compiler options, manually edit the Macros file. To add user defined preprocessor macros, modify ``ICE_CPPDEFS`` in **cice.settings** using the syntax ``-DCICE_MACRO``.
- To clean the build before each compile, set ``ICE_CLEANBUILD`` in **cice.settings** to true (this is the default value), or use the ``buildclean`` option (``-s buildclean``)  when creating the case with **cice.setup**.  To not clean before the build, set ``ICE_CLEANBUILD`` in **cice.settings** to false, or use the ``buildincremental`` option  (``-s buildincremental``) when creating the case with **cice.setup**.  It is recommended that the ``ICE_CLEANBUILD`` be set to true if there are any questions about whether the build is proceeding properly.

To build and run::

  ./cice.build
  ./cice.submit

The build and run log files will be copied into the logs subdirectory in the case directory.
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
  specifies the machine name.  This should be consistent with the name defined in the Macros and env files in **configurations/scripts/machines**.  This is required in all modes and is paired with ``--env`` to define the compilation environment.

``--env``,  ``-e`` ENVIRONMENT1,ENVIRONMENT2,ENVIRONMENT3
specifies the compilation environment associated with the machine.  This should be consistent with the name defined in the Macros and env files in **configurations/scripts/machines**.  Each machine can have multiple supported environments including support for different compilers, different compiler versions, different mpi libraries, or other system settigs.  When used with ``--suite`` or ``--test``, the ENVIRONMENT can be a set of comma deliminated values with no spaces and the tests will then be run for all of those environments.  With ``--case``, only one ENVIRONMENT should be specified. (default is intel)
  
``--pes``,  ``-p`` MxN[[xBXxBY[xMB]
  specifies the number of tasks and threads the case should be run on.  This only works with ``--case``.  The format is tasks x threads or "M"x"N" where M is tasks and N is threads and both are integers. BX, BY, and MB can also be set via this option where BX is the x-direction blocksize, BY is the y-direction blocksize, and MB is the max-blocks setting.  If BX, BY, and MB are not set, they will be computed automatically based on the grid size and the task/thread count.  More specifically, this option has three modes, --pes MxN, --pes MxNxBXxBY, and --pes MxNxBXxBYxMB.  (default is 4x1)

``--acct``  ACCOUNT
  specifies a batch account number.  This is optional.  See :ref:`account` for more information.

``--queue`` QUEUE
  specifies a batch queue name.  This is optional.  See :ref:`queue` for more information.

``--grid``, ``-g`` GRID
  specifies the grid.  This is a string and for the current CICE driver, gx1, gx3, and tx1 are supported. (default = gx3)

``--set``,  ``-s`` SET1,SET2,SET3
  specifies the optional settings for the case.  The settings for ``--suite`` are defined in the suite file.  Multiple settings can be specified by providing a comma deliminated set of values without spaces between settings.  The available settings are in **configurations/scripts/options** and ``cice.setup --help`` will also list them.  These settings files can change either the namelist values or overall case settings (such as the debug flag).  For cases and tests (not suites), settings defined in **~/.cice_set** (if it exists) will be included in the --set options.  This behaviour can be overridden with the `--ignore-user-set`` command line option.

``--ignore-user-set``
  ignores settings defined in **~/.cice.set** (if it exists) for cases and tests.  **~/.cice_set** is always ignored for test suites.

For CICE, when setting up cases, the ``--case`` and ``--mach`` must be specified.  
It's also recommended that ``--env`` be set explicitly as well.  
``--pes`` and ``--grid`` can be very useful.
``--acct`` and ``--queue`` are not normally used.  A more convenient method 
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
the case directory.  If multiple options are chosen that conflict, then the last
option chosen takes precedence.  Not all options are compatible with each other.

Settings defined in **~/.cice_set** (if it exists) will be included in the ``--set`` 
options.  This behaviour can be overridden with the `--ignore-user-set`` command 
line option.  The format of the **~/.cice_set** file is a identical to the
``--set`` option, a single comma-delimited line of options.  Settings on the 
command line will take precedence over settings defined in **~/.cice_set**.

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

Once the cases are created, users are free to modify the **cice.settings** and 
**ice_in** namelist to further modify their setup.

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

.. _cicecpps:

C Preprocessor (CPP) Macros
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are a number of C Preprocessing Macros supported in the CICE model.  These
allow certain coding features like NetCDF, MPI, or specific Fortran features to be 
excluded or included during the compile.  

The CPPs are defined by the `CPPDEFS` variable in the Makefile.  They are defined
by passing the -D[CPP] to the C and Fortran compilers (ie. -DUSE_NETCDF) and this
is what needs to be set in the `CPPDEFS` variable.  The value of `ICE_CPPDEFS` in
**cice.settings** is copied into the Makefile `CPPDEFS` variable as are settings
hardwired into the **Macros.[machine]_[environment]** file.

In general, ``-DFORTRANUNDERSCORE`` should always be set to support the Fortran/C
interfaces in **ice_shr_reprosum.c**.  In addition, if NetCDF is used, ``-DUSE_NETCDF``
should also be defined.  A list of available CPPs can be found in
:ref:`tabcpps`.

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
and the environment (compiler) is specified with ``--env``.  mach and env 
in combination define the compiler, compiler version, supporting libaries,
and batch information.  Multiple compilation environments can be created for
a single machine by choosing unique env names.
 
- cd to **configuration/scripts/machines/**

- Copy an existing env and a Macros file to new names for your new machine

- Edit your env and Macros files, update as needed

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

   "ICE_MACHINE_MACHNAME", "string", "machine name"
   "ICE_MACHINE_MACHINFO", "string", "machine information"
   "ICE_MACHINE_ENVNAME", "string", "env/compiler name"
   "ICE_MACHINE_ENVINFO", "string", "env/compiler information"
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
The order of precedence is **cice.setup** command line option, 
**.cice\_proj** setting, and then value in the **env.[machine]** file.

.. _queue:

Machine Queue Settings
~~~~~~~~~~~~~~~~~~~~~~~~

Supported machines will have a default queue specified by the variable ``ICE_MACHINE_QUEUE``
in the **env.[machine]** file.  This can also be manually changed in the **cice.run** or
**cice.test** scripts or even better, use the ``--queue`` option in **cice.setup**.

.. _laptops:

Porting to Laptop or Personal Computers
-----------------------------------------
To get the required software necessary to build and run CICE, and use the plotting and quality control scripts included in the repository, a `conda <https://docs.conda.io/en/latest/>`_ environment file is available at :

``configuration/scripts/machines/environment.yml``.

This configuration is supported by the Consortium on a best-effort basis on macOS and GNU/Linux. It is untested under Windows, but might work using the `Windows Subsystem for Linux <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_.

Once you have installed Miniconda and created the ``cice`` conda environment by following the procedures in this section, CICE should run on your machine without having to go through the formal :ref:`porting` process outlined above.

.. _install_miniconda:

Installing Miniconda
~~~~~~~~~~~~~~~~~~~~

We recommend the use of the `Miniconda distribution <https://docs.conda.io/en/latest/miniconda.html>`_ to create a self-contained conda environment from the ``environment.yml`` file.
This process has to be done only once.
If you do not have Miniconda or Anaconda installed, you can install Miniconda by following the `official instructions  <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_, or with these steps:

On macOS:

.. code-block:: bash

  # Download the Miniconda installer to ~/Downloads/miniconda.sh
  curl -L https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o ~/Downloads/miniconda.sh
  # Install Miniconda
  bash ~/Downloads/miniconda.sh
  
  # Follow the prompts
  
  # Close and reopen your shell


On GNU/Linux:

.. code-block:: bash

  # Download the Miniconda installer to ~/miniconda.sh
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
  # Install Miniconda
  bash ~/miniconda.sh
  
  # Follow the prompts
  
  # Close and reopen your shell

Note: on some Linux distributions (including Ubuntu and its derivatives), the csh shell that comes with the system is not compatible with conda.
You will need to install the tcsh shell (which is backwards compatible with csh), and configure your system to use tcsh as csh:
 
.. code-block:: bash
 
  # Install tcsh
  sudo apt-get install tcsh
  # Configure your system to use tcsh as csh
  sudo update-alternatives --set csh /bin/tcsh
 
  

.. _init_shell:

Initializing your shell for use with conda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We recommend initializing your default shell to use conda.
This process has to be done only once.

The Miniconda installer should ask you if you want to do that as part of the installation procedure.
If you did not answer "yes", you can use one of the following procedures depending on your default shell.
Bash should be your default shell if you are on macOS (10.14 and older) or GNU/Linux.

Note: answering "yes" during the Miniconda installation procedure will only initialize the Bash shell for use with conda.

If your Mac has macOS 10.15 or higher, your default shell is Zsh. 

These instructions make sure that the ``conda`` command is available when you start your shell by modifying your shell's startup file.
Also, they make sure not to activate the "base" conda environment when you start your shell.
This conda environment is created during the Miniconda installation but is not used for CICE. 

For Bash:

.. code-block:: bash

  # Install miniconda as indicated above, then initialize your shell to use conda:
  source $HOME/miniconda3/bin/activate
  conda init bash
  
  # Don't activate the "base" conda environment on shell startup
  conda config --set auto_activate_base false
  
  # Close and reopen your shell

For Zsh (Z shell):

.. code-block:: bash

  # Initialize Zsh to use conda
  source $HOME/miniconda3/bin/activate
  conda init zsh
  
  # Don't activate the "base" conda environment on shell startup
  conda config --set auto_activate_base false
  
  # Close and reopen your shell

For tcsh:

.. code-block:: bash
  
  # Install miniconda as indicated above, then initialize your shell to use conda:
  source $HOME/miniconda3/etc/profile.d/conda.csh
  conda init tcsh
  
  # Don't activate the "base" conda environment on shell startup
  conda config --set auto_activate_base false
  
  # Close and reopen your shell

For fish:

.. code-block:: bash
  
  # Install miniconda as indicated above, then initialize your shell to use conda:
  source $HOME/miniconda3/etc/fish/conf.d/conda.fish
  conda init fish
  
  # Don't activate the "base" conda environment on shell startup
  conda config --set auto_activate_base false
  
  # Close and reopen your shell

For xonsh:

.. code-block:: bash

  # Install miniconda as indicated above, then initialize your shell to use conda:
  source-bash $HOME/miniconda3/bin/activate
  conda init xonsh
  
  # Don't activate the "base" conda environment on shell startup
  conda config --set auto_activate_base false
  
  # Close and reopen your shell

.. _init_shell_manually:

Initializing your shell for conda manually
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you prefer not to modify your shell startup files, you will need to run the appropriate ``source`` command below (depending on your default shell) before using any conda command, and before compiling and running CICE.
These instructions make sure the ``conda`` command is available for the duration of your shell session.

For Bash and Zsh:

.. code-block:: bash

  # Initialize your shell session to use conda:
  source $HOME/miniconda3/bin/activate

For tcsh:

.. code-block:: bash
  
  # Initialize your shell session to use conda:
  source $HOME/miniconda3/etc/profile.d/conda.csh


For fish:

.. code-block:: bash
  
  # Initialize your shell session to use conda:
  source $HOME/miniconda3/etc/fish/conf.d/conda.fish

For xonsh:

.. code-block:: bash

  # Initialize your shell session to use conda:
  source-bash $HOME/miniconda3/bin/activate


.. _create_conda_env:

Creating CICE directories and the conda environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The conda configuration expects some directories and files to be present at ``$HOME/cice-dirs``:

.. code-block:: bash

  cd $HOME
  mkdir -p cice-dirs/runs cice-dirs/baseline cice-dirs/input
  # Download the required forcing from https://github.com/CICE-Consortium/CICE/wiki/CICE-Input-Data
  # and untar it at $HOME/cice-dirs/input

This step needs to be done only once.

If you prefer that some or all of the CICE directories be located somewhere else, you can create a symlink from your home to another location:

.. code-block:: bash

  
  # Create the CICE directories at your preferred location
  cd ${somewhere}
  mkdir -p cice-dirs/runs cice-dirs/baseline cice-dirs/input
  # Download the required forcing from https://github.com/CICE-Consortium/CICE/wiki/CICE-Input-Data
  # and untar it at cice-dirs/input
  
  # Create a symlink to cice-dirs in your $HOME
  cd $HOME
  ln -s ${somewhere}/cice-dirs cice-dirs

Note: if you wish, you can also create a complete machine port for your computer by leveraging the conda configuration as a starting point. See :ref:`porting`.

Next, create the "cice" conda environment from the ``environment.yml`` file in the CICE source code repository.  You will need to clone CICE to run the following command:

.. code-block:: bash

  conda env create -f configuration/scripts/machines/environment.yml

This step needs to be done only once.

.. _using_conda_env:

Using the conda configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Follow the general instructions in :ref:`overview`, using the ``conda`` machine name and ``macos`` or ``linux`` as compiler names.

On macOS:

.. code-block:: bash

  ./cice.setup -m conda -e macos -c ~/cice-dirs/cases/case1
  cd ~/cice-dirs/cases/case1
  ./cice.build
  ./cice.run

On GNU/Linux:

.. code-block:: bash

  ./cice.setup -m conda -e linux -c ~/cice-dirs/cases/case1
  cd ~/cice-dirs/cases/case1
  ./cice.build
  ./cice.run

A few notes about the conda configuration:

- This configuration always runs the model interactively, such that ``./cice.run`` and ``./cice.submit`` are the same.
- You should not update the packages in the ``cice`` conda environment, nor install additional packages.
- Depending on the numbers of CPUs in your machine, you might not be able to run with the default MPI configuration (``-p 4x1``). You likely will get an OpenMPI error such as:

    There are not enough slots available in the system to satisfy the 4 slots that were requested by the application:  ./cice
    
  You can run CICE in serial mode by specifically requesting only one process:
  
  .. code-block:: bash
  
    ./cice.setup -m conda -e linux -p 1x1 ...
  
  If you do want to run with more MPI processes than the number of available CPUs in your machine, you can add the ``--oversubscribe`` flag to the ``mpirun`` call in ``cice.run``:
  
  .. code-block:: bash
  
    # For a specific case:
    # Open cice.run and replace the line
    mpirun -np <num> ./cice >&! $ICE_RUNLOG_FILE
    # with
    mpirun -np <num> --oversubscribe ./cice >&! $ICE_RUNLOG_FILE
  
    # For all future cases:
    # Open configuration/scripts/cice.launch.csh and replace the line
    mpirun -np ${ntasks} ./cice >&! \$ICE_RUNLOG_FILE
    # with
    mpirun -np ${ntasks} --oversubscribe ./cice >&! \$ICE_RUNLOG_FILE
  
- It is not recommeded to run other test suites than ``quick_suite`` or ``travis_suite`` on a personal computer.
- The conda environment is automatically activated when compiling or running the model using the ``./cice.build`` and ``./cice.run`` scripts in the case directory. These scripts source the file ``env.conda_{linux.macos}``, which calls ``conda activate cice``.
- To use the "cice" conda environment with the Python plotting (see :ref:`timeseries`) and quality control scripts (see :ref:`CodeCompliance`), you must manually activate the environment:

  .. code-block:: bash
  
    cd ~/cice-dirs/cases/case1
    conda activate cice
    python timeseries.py ~/cice-dirs/cases/case1/logs
    conda deactivate  # to deactivate the environment
  
- The environment also contains the Sphinx package necessesary to build the HTML documentation :

  .. code-block:: bash
  
    cd doc
    conda activate cice
    make html
    # Open build/html/index.html in your browser
    conda deactivate  # to deactivate the environment


.. _force:

Forcing data
------------

The input data space is defined on a per machine basis by the ``ICE_MACHINE_INPUTDATA`` 
variable in the **env.[machine]** file.  That file space is often shared among multiple 
users, and it can be desirable to consider using a common file space with group read 
and write permissions such that a set of users can update the inputdata area as 
new datasets are available.

CICE input datasets are stored on an anonymous ftp server.  More information about
how to download the input data can be found at https://github.com/CICE-Consortium/CICE/wiki/CICE-Input-Data.
Test forcing datasets are available for various grids at the ftp site.  
These data files are designed only for testing the code, not for use in production runs 
or as observational data. Please do not publish results based on these data sets.


Run Directories
---------------

The **cice.setup** script creates a case directory.  However, the model 
is actually built and run under the ``ICE_OBJDIR`` and ``ICE_RUNDIR`` directories
as defined in the **cice.settings** file.  It's important to note that when the
run scripts are submitted, the current **ice_in**, **cice.settings**, and **env.[machine]**
files are copied from the case directory into the run directory.  Users should 
generally not edit files in the run directory as these are overwritten when following
the standard workflow.

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
