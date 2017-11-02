:tocdepth: 3

Introduction - CICE5
============================================

The Los Alamos sea ice model (CICE) is the result of an effort to
develop a computationally efficient sea ice component for a fully
coupled atmosphere--land global climate model. It was
designed to be compatible with the Parallel Ocean Program
(POP), an ocean circulation model developed at 
Los Alamos National Laboratory for use on massively parallel computers
:cite:`SDM92,DSM93,DSM94`. The current version of the
model has been enhanced greatly through collaborations with members of
the community.

CICE has several interacting components: a thermodynamic model that
computes local growth rates of snow and ice due to vertical conductive, 
radiative and turbulent fluxes, along with snowfall; a model of ice 
dynamics, which predicts the velocity field of the ice pack based on 
a model of the material strength of the ice; a transport model that 
describes advection of the areal concentration, ice volumes and other 
state variables; and a ridging parameterization that transfers ice among
thickness categories based on energetic balances and 
rates of strain.External routines would prepare and execute data exchanges with an
external "flux coupler," which then passes the data to other climate
model components such as POP.

This model release is CICE version 5.1, available from http://oceans11.lanl.gov/trac/CICE/wiki.
It updates CICE5.0, which was released in September 2013. With so many new parameterizations,
we must admit that all combinations have not been tested.  Also, different parameterizations for 
various sub-processes (e.g., snow infiltration by melt or sea water) have been introduced as part 
of the new code options, which need to be unified or pulled out as separate options.  

This document uses the following text conventions:
Variable names used in the code are ``typewritten``.
Subroutine names are given in *italic*.
File and directory names are in **boldface**.
A comprehensive :ref:`index`, including glossary of symbols with many of their values, appears
at the end of this guide.

=================
Quick Start guide
=================

~~~~~~~~~~~~~
Get the model
~~~~~~~~~~~~~

Checkout the model from the CICE-Consortium repository,

  github.com/CICE-Consortium

For more details about how to work in github with CICE, a document can be 
found `here <https://docs.google.com/document/d/1rR6WAvZQT9iAMUp-m_HZ06AUCCI19mguFialsMCYs9o>`_.

~~~~~~~~~~~~~~~~~
Running the model
~~~~~~~~~~~~~~~~~

> cd consortium

> ./create.case -c ~/mycase1 -g gx3 -m thunder -s diag1,thread -p 8x1

> cd ~/mycase1

> ./cice.build

> ./cice.submit/Users/duvivier/Documents/Research/github/CICE-Consortium/CICE/doc/source/all_orig/cice_2_quick_start.rst

~~~~~~~~~~~~
More Details
~~~~~~~~~~~~

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

~~~~~~~~~~~~~~~~
Online resources
~~~~~~~~~~~~~~~~

**DO WE WANT TO KEEP THESE?**

primary wiki page:


FAQ:


instructions for code developers:


ongoing or planned development projects:


list of users and publications:


Please send references to your publications using the CICE model to ...


Please report any bugs to 
Elizabeth Hunke (eclare@lanl.gov)
 
Good luck!


=============
Major updates
=============

~~~~~~~~~
CICE V5.1
~~~~~~~~~

- include ice velocity in atm-ice coupling updates (e.g. stress) for high-frequency coupling
- allow a variable coefficient for the ice-ocean heat flux
- several new namelist options improve flexibility, especially for coupled model configurations:
   - ice-ocean heat flux
   - 'virtual' or 'real' topo melt pond water
   - ocean freezing temperature
   - high-frequency coupling
   - coupling and computational grids may be different
   - and more
- additional software enhancements improve flexibility and compatibility with CESM, Hadley Centre, and U.S. Navy coupled models
- new diagnostics and updated documentation
- various bug fixes 

~~~~~~~~~
CICE V5.0
~~~~~~~~~

- A method for prognosing sea ice salinity, including improved snow-ice formation
- Two new explicit melt pond parameterizations (topo and level-ice)
- Sea ice biogeochemistry
- Elastic-Anisotropic-Plastic rheology
- Improved parameterization of form drag 
- The "revised EVP" under-damping approach
- Gracefully handles the case when an internal layer melts completely
- Gregorian calendar with leap years
- Reduced memory and floating-point operations for tracer calculations
- Ice and snow enthalpy defined as tracers
- New history variables for melt ponds, ridging diagnostics, biogeochemistry and more
- Read/write variables on the extended grid, including ghost (halo) cells
- Parallelism option via OpenMP threads
- Improved parallel efficiency through halo masks and new block distributions
- Parallel I/O option via the PIO library
- Restarts in binary or netCDF formats
- CPP options for categories, layers and tracers 
- Corrected bugs, particularly for nonstandard configurations.

======================
Acknowledgements
======================
This work has been supported under the Department of Energy’s Climate,
Ocean and Sea Ice Modeling project through the Computer Hardware Applied
Mathematics and Model Physics (CHAMMP) program, Climate Change
Prediction Program (CCPP), Improving the Characterization of Clouds,
Aerosols and the Cryosphere in Climate Models (Cloud-Cryo) program and
Scientific Discovery through Advanced Computing (SCIDAC) program, with
additional support from the T-3 Fluid Dynamics and Solid Mechanics Group
at Los Alamos National Laboratory. Special thanks are due to the
following people:

-  members of the CESM Polar Climate Working Group, including David
   Bailey, Alice DuVivier, Cecilia Bitz, Bruce Briegleb, Tony Craig, 
   Marika Holland, John Dennis, Julie Schramm, Bonnie Light and Phil Jones.

-  Andrew Roberts of the Naval Postgraduate School,

-  David Hebert and Olivier Lecomte for their melt pond work,

-  Jonathan Gregory of the University of Reading and the U.K. MetOffice
   for supplying tripole T-fold code and documentation,

-  Alison McLaren, Ann Keen and others working with the Hadley Centre
   GCM for testing non-standard model configurations and providing their
   code to us,

-  Daniel Feltham and his research group for several new
   parameterizations and documentation,

-  Sylvain Bouillon for the revised EVP approach,

-  the many researchers who tested beta versions of CICE 5 and waited
   patiently for the official release.

======================
Copyright
======================
© Copyright 2013, LANS LLC. All rights reserved. Unless otherwise
indicated, this information has been authored by an employee or
employees of the Los Alamos National Security, LLC (LANS), operator of
the Los Alamos National Laboratory under Contract No. DE-AC52-06NA25396
with the U.S. Department of Energy. The U.S. Government has rights to
use, reproduce, and distribute this information. The public may copy and
use this information without charge, provided that this Notice and any
statement of authorship are reproduced on all copies. Neither the
Government nor LANS makes any warranty, express or implied, or assumes
any liability or responsibility for the use of this information.
Beginning with version 4.0, the CICE code carries Los Alamos Software
Release number LA-CC-06-012.


