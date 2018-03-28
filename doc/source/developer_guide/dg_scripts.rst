:tocdepth: 3

.. _dev_scripts:

Scripts Implementation
========================

The scripts are the third part of the icepack package.  They support setting up
cases, building, and running the icepack stand-alone model.

File List
--------------

The directory structure under configure/scripts is as follows.

| **configuration/scripts/**
|        **Makefile**              primary makefile
|        **icepack.batch.csh**     creates batch scripts for particular machines
|        **icepack.build**         compiles the code
|        **icepack.launch.csh**    creates script logic that runs the executable
|        **icepack.run.setup.csh** sets up the run scripts
|        **icepack.run.suite.csh** sets up the test suite
|        **icepack.settings**      defines environment, model configuration and run settings
|        **icepack.test.setup.csh**   creates configurations for testing the model
|        **icepack_decomp.csh**    defines the grid size
|        **icepack_in**            namelist input data
|        **machines/**             machine specific files to set env and Macros
|        **makdep.c**              determines module dependencies
|        **options/**              other namelist configurations available from the icepack.setup command line
|        **parse_namelist.sh**     replaces namelist with command-line configuration
|        **parse_namelist_from_settings.sh**   replaces namelist with values from icepack.settings
|        **parse_settings.sh**     replaces settings with command-line configuration
|        **tests/**                scripts for configuring and running basic tests

.. _dev_strategy:

Strategy
-----------

The icepack scripts are implemented such that everything is resolved after
**icepack.setup** is called.  This is done by both copying specific files
into the case directory and running scripts as part of the **icepack.setup**
command line to setup various files.

**icepack.setup** drives the case setup.  It is written in csh.  All supporting
scripts are relatively simple csh or sh scripts.

The file **icepack.settings** specifies a set of env defaults for the case.  The file
**icepack_in** defines the namelist input for the icepack driver.

.. _dev_options:

Preset Case Options
---------------------


``icepack.setup -s`` option allows the user to choose some predetermined icepack
settings and namelist.  Those options are defined in **configurations/scripts/options/**
and the files are prefixed by either set_env, set_nml, or test_nml.  When **icepack.setup**
is executed, the appropriate files are read from **configurations/scripts/options/**
and the **icepack.settings** and/or **icepack_in** files are updated in the case directory
based on the values in those files.

The filename suffix determines the name of the -s option.  So, for instance, 

  ``icepack.setup -s diag1,debug,bgcISPOL``

will search for option files with suffixes of diag1, debug, and bgcISPOL and then
apply those settings.  

**parse_namelist.sh**, **parse_settings.sh**, and **parse_namelist_from_settings.sh** 
are the three scripts that modify **icepack_in** and **icepack.settings**.

To add new options, just add new files to the **configurations/scripts/options/** directory
with appropriate names and syntax.  The set_nml file syntax is the same as namelist
syntax and the set_env files are consistent with csh setenv syntax.  See other files for
examples of the syntax.

.. _dev_machines:

Machines
-----------

Machine specific information is contained in **configuration/scripts/machines**.  That
directory contains a Macros file and an env file for each supported machine.
One other files will need to be
changed to support a port, that is **configuration/scripts/icepack.batch.csh**.
To port to a new machine, see :ref:`porting`.  

.. _dev_testing:

Test scripts
-------------

Under **configuration/scripts/tests** are several files including the scripts to 
setup the smoke and restart tests (**test_smoke.script**, **test_restart.script*).
A baseline test script (**baseline.script**) is also there to setup the regression
and comparison testing.  That directory also contains the preset test suites 
(ie. **base_suite.ts**) and a file that supports post-processing on the model
output (**timeseries.csh**).  

There is a subdirectory, **configuration/scripts/tests/CTest**, that supports the
CTest scripts.  These scripts allow test reporting to CDash.

To add a new test, a file associated with that test will need to be added to the
**configuration/scripts/tests** directory similar to **test_smoke.script** 
and **test_restart.script**.  In addition, some new options files in 
**configuration/scripts/options** may need to be added similar to **test_nml.restart1**,
**test_nml.restart2**, and **set_nml.restart**.  
