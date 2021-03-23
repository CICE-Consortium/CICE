:tocdepth: 3

.. _dev_scripts:

Scripts 
========

The scripts are the third part of the cice package.  They support setting up
cases, building, and running the cice stand-alone model.

File List
--------------

The directory structure under configure/scripts is as follows.

| **configuration/scripts/**
|        **Makefile**              primary makefile
|        **cice.batch.csh**        creates batch scripts for particular machines
|        **cice.build**            compiles the code
|        **cice.decomp.csh**       computes a decomposition given a grid and task/thread count
|        **cice.launch.csh**       creates script logic that runs the executable
|        **cice.run.setup.csh**    sets up the run scripts
|        **cice.settings**         defines environment, model configuration and run settings
|        **cice.test.setup.csh**   creates configurations for testing the model
|        **ice_in**                namelist input data
|        **machines/**             machine specific files to set env and Macros
|        **makdep.c**              determines module dependencies
|        **options/**              other namelist configurations available from the cice.setup command line
|        **parse_namelist.sh**     replaces namelist with command-line configuration
|        **parse_namelist_from_settings.sh**   replaces namelist with values from cice.settings
|        **parse_settings.sh**     replaces settings with command-line configuration
|        **setup_run_dirs.csh**    creates the case run directories
|        **set_version_number.csh** updates the model version number from the **cice.setup** command line
|        **timeseries.csh**        generates PNG timeseries plots from output files, using GNUPLOT
|        **timeseries.py**         generates PNG timeseries plots from output files, using Python
|        **tests/**                scripts for configuring and running basic tests

.. _dev_strategy:

Strategy
-----------

The cice scripts are implemented such that everything is resolved after
**cice.setup** is called.  This is done by both copying specific files
into the case directory and running scripts as part of the **cice.setup**
command line to setup various files.

**cice.setup** drives the case setup.  It is written in csh.  All supporting
scripts are relatively simple csh or sh scripts.  See :ref:`scripts` for additional
details.

The file **cice.settings** specifies a set of env defaults for the case.  The file
**ice_in** defines the namelist input for the cice driver.


.. _dev_preset_options:

Preset Case Options
---------------------

The ``cice.setup --set`` option allows the user to choose some predetermined cice
settings and namelist.  Those options are defined in **configurations/scripts/options/**
and the files are prefixed by either set_env or set_nml.  When **cice.setup**
is executed, the appropriate files are read from **configurations/scripts/options/**
and the **cice.settings** and/or **ice_in** files are updated in the case directory
based on the values in those files.

The filename suffix determines the name of the -s option.  So, for instance, 

  ``cice.setup -s diag1,debug,bgcISPOL``

will search for option files with suffixes of diag1, debug, and bgcISPOL and then
apply those settings.  

**parse_namelist.sh**, **parse_settings.sh**, and **parse_namelist_from_settings.sh** 
are the three scripts that modify **ice_in** and **cice.settings**.

To add new options, just add new files to the **configurations/scripts/options/** directory
with appropriate names and syntax.  The set_nml file syntax is the same as namelist
syntax and the set_env files are consistent with csh setenv syntax.  See other files for
examples of the syntax.

.. _build:

Build Scripts
--------------

CICE uses GNU Make to build the model.  There is a common **Makefile** for all machines.  
Each machine provides a Macros file to define some Makefile variables
and and an env file to specify the modules/software stack for each compiler.
The machine is built by the cice.build script which invokes Make.
There is a special trap for circular dependencies in the cice.build script to
highlight this error when it occurs.

The **cice.build** script has some additional features including the ability to 
pass a Makefile target.  This is documented in :ref:`cicebuild`.  In addition, there
is a hidden feature in the **cice.build** script that allows for reuse of 
executables.  This is used by the test suites to significantly reduce cost of
building the model.  It is invoked with the ``--exe`` argument to **cice.build**
and should not be invoked by users interactively.

.. _dev_machines:

Machines
-----------

Machine specific information is contained in **configuration/scripts/machines**.  That
directory contains a Macros file and an env file for each supported machine.
One other files will need to be
changed to support a port, that is **configuration/scripts/cice.batch.csh**.
To port to a new machine, see :ref:`porting`.  

.. _dev_test_options:

Test Options
---------------

Values that are associated with the `--sets` cice.setup are defined in 
**configuration/scripts/options**.  Those files are text files and cice.setup
uses the values in those files to modify the `cice.settings` and `ice_in` files
in the case as the case is created.  Files name `set_env.$option` are associated
with values in the `cice.settings` file.  Files named `set_nml.$option` are associated
with values in `ice.in`.  These files contain simple keyword pair values one line
at a time.  A line starting with # is a comment.  Files names that start with `test_`
are used specifically for tests.

That directory also contains files named `set_files.$option`.  This provides an
extra layer on top of the individual setting files that allows settings to be
defined based on groups of other settings.  The `set_files.$option` files
contain a list of `--sets` options to be applied.  

The $option part of the filename is the argument to `--sets` argument in `cice.setup`.
Multiple options can be specified by creating a comma delimited list.  In the case
where settings contradict each other, the last defined is used.

.. _dev_testing:

Test scripts
-------------

Under **configuration/scripts/tests** are several files including the scripts to 
setup the various tests, such as smoke and restart tests (**test_smoke.script**, **test_restart.script**)
and the files that describe with options files are needed for each test (ie. **test_smoke.files**, **test_restart.files**).
A baseline test script (**baseline.script**) is also there to setup the general regression
and comparison testing.  That directory also contains the preset test suites 
(ie. **base_suite.ts**) and a script (**report_results.csh**) that pushes results from 
test suites back to the CICE-Consortium test results wiki page.

To add a new test (for example newtest), several files may be needed,

- **configuration/scripts/tests/test_newtest.script** defines how to run the test.  This chunk
  of script will be incorporated into the case test script
- **configuration/scripts/tests/test_newtest.files** list the set of options files found in
  **configuration/scripts/options/** needed to
  run this test.  Those files will be copied into the test directory when the test is invoked
  so they are available for the **test_newtest.script** to use.
- some new files may be needed in **configuration/scripts/options/**.  These could be relatively
  generic **set_nml** or **set_env** files, or they could be test specific files typically carrying
  a prefix of **test_nml**.

Generating a new test, particularly the **test_newtest.script** usually takes some iteration before
it's working properly.

.. _dev_validation:

Code Validation Script
----------------------

The code validation (aka QC or quality control) test validates non bit-for-bit model changes.  The directory 
**configuration/scripts/tests/QC** contains scripts related to the validation testing,
and this process is described in :ref:`validation`.  This section will describe a set
of scripts that test and validate the QC process.  This should be done 
when the QC test or QC test scripts (i.e., ``cice.t-test.py``) are modified.  
Again, this section **documents a validation process for the QC scripts**; it does not
describe to how run the validation test itself.  

Two scripts have been created to automatically validate the QC script.  
These scripts are:

* ``gen_qc_cases.csh``, which creates the 4 test cases required for validation,
  builds the executable, and submits to the queue.
* ``compare_qc_cases.csh``, which runs the QC script on three combinations
  of the 4 test cases and outputs whether or not the correct response was received.

The ``gen_qc_cases.csh`` script allows users to pass some arguments similar
to the ``cice.setup`` script.  These options include:

* ``--mach, -m``: Machine (REQUIRED)
* ``--env,  -e``: Compiler
* ``--pes,  -p``: tasks x threads
* ``--acct``    : Account number for batch submission
* ``--grid, -g``: Grid
* ``--queue``   : Queue for the batch submission
* ``--testid``  : test ID, user-defined id for testing

The script creates 4 test cases, with testIDs ``qc_base``, ``qc_bfb``, ``qc_nonbfb``,
and ``qc_fail``.  ``qc_base`` is the base test case with the default QC namelist.
``qc_bfb`` is identical to ``qc_base``.  ``qc_nonbfb`` is a test that is not bit-for-bit
when compared to ``qc_base``, but not climate changing.  ``qc_fail`` is a test that is not
bit-for-bit and also climate changing.

In order to run the ``compare_qc_cases.csh`` script, the following requirements must be met:

* Python v2.7 or later
* netcdf Python package
* numpy Python package

To install the necessary Python packages, the ``pip`` Python utility can be used.

.. code-block:: bash

  pip install --user netCDF4
  pip install --user numpy

**Note:** Some machines might report ``pip: Command not found.``  If you encounter this error,
check to see if there is any Python module (``module avail python``) that you might need
to load prior to using ``pip``.

To perform the QC validation, execute the following commands.

.. code-block:: bash

  # From the CICE base directory
  cp configuration/scripts/tests/QC/gen_qc_cases.csh .
  cp configuration/scripts/tests/QC/compare_qc_cases.csh
  
  # Create the required test cases
  ./gen_qc_cases.csh -m <machine> --acct <acct>

  # Wait for all 4 jobs to complete

  # Perform the comparisons
  ./compare_qc_cases.csh

The ``compare_qc_cases.csh`` script will run the QC script on the following combinations:

* ``qc_base`` vs. ``qc_bfb``
* ``qc_base`` vs. ``qc_nonbfb``
* ``qc_base`` vs. ``qc_fail``

An example of the output from ``compare_qc_cases.csh`` is shown below.::

  ===== Running QC tests and writing output to validate_qc.log =====
  Running QC test on base and bfb directories.
  Expected result: PASSED
  Result: PASSED
  -----------------------------------------------
  Running QC test on base and non-bfb directories.
  Expected result: PASSED
  Result: PASSED
  -----------------------------------------------
  Running QC test on base and climate-changing directories.
  Expected result: FAILED
  Result: FAILED
  
  
  QC Test has validated

