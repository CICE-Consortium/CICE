Quick Start guide
============================================

Get the model
-------------

Checkout the model from the CICE-Consortium repository,

  github.com/CICE-Consortium

For more details about how to work in github with CICE, a document can be 
found `here <https://docs.google.com/document/d/1rR6WAvZQT9iAMUp-m_HZ06AUCCI19mguFialsMCYs9o>`_.


Running the model
-----------------

> cd consortium

> ./create.case -c ~/mycase1 -g gx3 -m thunder -s diag1,thread -p 8x1

> cd ~/mycase1

> ./cice.build

> ./cice.submit


More Details:
-------------

create.case generates a case, use "create.case -h" for help with the tool.
  -c is the case name and location (required)

  -m is the machine name (required). Currently, there are working ports for NCAR yellowstone and cheyenne, AFRL thunder, NavyDSRC gordon and conrad, and LANL’s wolf machines.

  -g is the resolution (default is gx3)

  -p is the task x thread/task values (default is 4x1)

  -s are comma separated optional env or namelist settings (default is "null")

  -t is the test name and location (cannot be used with -c).

  -bd is used to specify the location of the baseline datasets (only used with -t)

  -bg is used to specify the cice version name for generating baseline datasets (only used with -t)

  -bc is used to specify the cice versoin name for comparison. I.e., the version name for the baseline dataset (only used with -t)

  -testid is used to specify a test ID (used only with -t or -ts)

  -ts is used to generate all test cases for a given test suite.


Several files are placed in the case directory

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
dependence on cpps.  To change this behavior, edit the env variable ICE_CLEANBUILD in
cice.settings.  

The cice.submit script just submits the cice.run script.  You can use cice.submit or just
submit the cice.run script on the command line.

The model will run in the directory defined by the env variable CICE_RUNDIR in cice.settings.  
Build and run logs will be copied into the case logs directory when complete.

To port, an env.machine and Macros.machine file have to be added to scripts/machines and the cice.run.setup.csh file needs to be modified.
 - cd to consortium/scripts/machines
 - Copy an existing env and Macros file to new names for your new machine
 - Edit the env and Macros file
 - cd to consortium/scripts
 - Edit the cice.run.setup.csh script to add a section for your machine for the batch settings and for the job launch settings
 - Download and untar the 1997 dataset to the location defined by ICE_MACHINE_INPUTDATA in the env file
 - Create a file in your home directory called .cice_proj and add your preferred account name to the first line.
 - You can now create a case and test.  If there are problems, you can manually edit the env, Macros, and cice.run files in the case directory until things are working properly.  Then you can copy the env and Macros files back to consortium/scripts/machines.  You will have to manually modify the cice.run.setup.csh script if there any changes needed there.


Forcing data
------------

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


Online resources
----------------

**DO WE WANT TO KEEP THESE?**

primary wiki page:


FAQ:


instructions for code developers:


ongoing or planned development projects:


list of users and publications:


Please send references to your publications using the CICE model to ...


Authors
-------

Please report any bugs to 
Elizabeth Hunke (eclare@lanl.gov)
 
Good luck!
