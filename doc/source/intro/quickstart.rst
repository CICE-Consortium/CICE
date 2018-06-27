:tocdepth: 3


.. _quickstart:

Quick Start
===========

Download the model from the CICE-Consortium repository, 
    https://github.com/CICE-Consortium/CICE

Instructions for working in github with CICE (and Icepack) can be
found in the `CICE Git and Workflow Guide <https://docs.google.com/document/d/1rR6WAvZQT9iAMUp-m_HZ06AUCCI19mguFialsMCYs9o>`_.

You will probably have to download some inputdata, see the `CICE wiki <https://github.com/cice-consortium/CICE/wiki>`_ or :ref:`force`.

From your main CICE directory, execute::

  ./cice.setup -c ~/mycase1 -g gx3 -m testmachine -s diag1,thread -p 8x1
  cd ~/mycase1
  ./cice.build
  ./cice.submit


``testmachine`` is a generic machine name included with the cice scripts.
The local machine name will have to be substituted for ``testmachine`` and
there are working ports for several different machines.  However, it may be necessary
to port the model to a new machine.  See :ref:`porting` for 
more information about how to port and :ref:`scripts` for more information about 
how to use the cice.setup script.

Please cite any use of the CICE code. More information can be found at :ref:`citing`.

~~~~~~~~~~~~
More Details
~~~~~~~~~~~~

``cice.setup -h`` will provide the latest information about how to use the tool.
``cice.setup --help`` will provide an extended version of the help.
There are three usage modes,

* ``--case`` or ``-c`` creates individual stand alone cases.
* ``--test`` creates individual tests.  Tests are just cases that have some extra automation in order to carry out particular tests such as exact restart.
* ``--suite`` creates a test suite.  Test suites are predefined sets of tests and ``--suite`` provides the ability to quick setup, build, and run a full suite of tests.

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

 - env.${machine} defines the environment

 - cice.settings defines many variables associated with building and running the model

 - makdep.c is a tool that will automatically generate the make dependencies

 - Macros.${machine} defines the Makefile Macros

 - Makefile is the makefile used to build the model

 - cice.build is a script that build the model

 - ice_in is the namelist file

 - cice.run is a batch run script

 - cice.submit is a simple script that submits the cice.run script

Once the case is created, all scripts and namelist are fully resolved.  Users can edit any
of the files in the case directory manually to change the model configuration.  The file
dependency is indicated in the above list.  For instance, if any of the files before
cice.build in the list are edited, cice.build should be rerun.

The casescripts directory holds scripts used to create the case and can largely be ignored.  

In general, when cice.build is executed, the model will build from scratch due to the large
dependence on cpps.  To change this behavior, edit the env variable ``ICE_CLEANBUILD`` in
cice.settings.  

The cice.submit script just submits the cice.run script.  You can use cice.submit or just
submit the cice.run script on the command line.

The model will run in the directory defined by the env variable ``ICE_RUNDIR`` in cice.settings.  
Build and run logs will be copied into the case logs directory when complete.

To port, an env.machine and Macros.machine file have to be added to scripts/machines and the cice.run.setup.csh file needs to be modified.
 - cd to consortium/scripts/machines
 - Copy an existing env and Macros file to new names for your new machine
 - Edit the env and Macros file
 - cd to consortium/scripts
 - Edit the cice.run.setup.csh script to add a section for your machine for the batch settings and for the job launch settings
 - Download and untar the 1997 dataset to the location defined by ``ICE_MACHINE_INPUTDATA`` in the env file
 - Create a file in your home directory called .cice_proj and add your preferred account name to the first line.
 - You can now create a case and test.  If there are problems, you can manually edit the env, Macros, and cice.run files in the case directory until things are working properly.  Then you can copy the env and Macros files back to consortium/scripts/machines.  You will have to manually modify the cice.run.setup.csh script if there any changes needed there.

~~~~~~~~~~~~
Forcing data
~~~~~~~~~~~~

The code is currently configured to run in standalone mode on a 3 degree grid using 
atmospheric data from 1997, available as detailed on the `wiki <https://github.com/CICE-Consortium/CICE/wiki/Testing-CICE>`_.
These data files are designed only for testing the code, not for use in production 
runs or as observational data.  Please do not publish results based on these data
sets.  Module cicecore/dynamics/cicedynB/ice_forcing.F90 can be modified to change the 
forcing data. 

As currently configured, the model runs on 4 processors.  MPI is used for message passing 
between processors, and OpenMP threading is available.  The grid provided here is too 
small for the code to scale well beyond about 8 processors. A 1 degree grid is provided also, 
and details about this grid can be found on the `wiki <https://github.com/CICE-Consortium/CICE/wiki/Testing-CICE>`_.

