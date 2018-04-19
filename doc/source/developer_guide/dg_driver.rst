:tocdepth: 3

.. _dev_driver:


Driver and Coupling Implementation
====================================

The driver and coupling layer is found in **cicecore/drivers/**.  The standalone driver is found
under **cicecore/drivers/cice/** and other high level coupling layers are found in other directories.
In general, CICE will build with only one of these drivers, depending how the model is run and
coupled.  Within the **cicecore/drivers/cice/** directory, the following files are found,

**CICE.F90** is the top level program file and that calls CICE_Initialize, CICE_Run, and CICE_Finalize methods.
**CICE_InitMod.F90** contains the CICE_Initialize method and other next level source code.
**CICE_RunMod.F90** contains the CICE_Run method and other next level source code.
**CICE_FinalMod.F90 ** contains the CICE_Finalize method and other next level source code.

Other **cicecore/drivers/** directories are similarly implemented with a top level coupling layer,
that is largely specified by an external coupled system and then some version of the CICE_InitMod.F90,
CICE_RunMod.F90, and CICE_FinalMod.F90 files.


Calling Sequence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~