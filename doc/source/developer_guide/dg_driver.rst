:tocdepth: 3

.. _dev_driver:

Driver Implementation
========================

The icepack driver is Fortran source code and exists to test the column physics
in a stand-alone mode for some simple column configurations.

File List
-------------------

The icepack driver consists of the following files 

|  **configuration/driver/**       driver for testing Icepack in stand-alone mode
|        **icedrv_MAIN.F90**        main program
|        **icedrv_InitMod.F90**     routines for initializing a run
|        **icedrv_RunMod.F90**      main driver routines for time stepping
|        **icedrv_arrays_column.F90**    essential arrays to describe the state of the ice
|        **icedrv_calendar.F90**    keeps track of what time it is
|        **icedrv_constants.F90**   physical and numerical constants and parameters
|        **icedrv_diagnostics.F90** miscellaneous diagnostic and debugging routines
|        **icedrv_diagnostics_bgc.F90**  diagnostic routines for biogeochemistry
|        **icedrv_domain_size.F90** domain sizes
|        **icedrv_flux.F90**        fluxes needed/produced by the model
|        **icedrv_forcing.F90**     routines to read and interpolate forcing data for stand-alone model runs
|        **icedrv_forcing_bgc.F90** routines to read and interpolate forcing data for bgc stand-alone model runs
|        **icedrv_init.F90**        general initialization routines
|        **icedrv_init_column.F90** initialization routines specific to the column physics
|        **icedrv_restart.F90**     driver for reading/writing restart files
|        **icedrv_restart_bgc.F90**  restart routines specific to the column physics
|        **icedrv_restart_shared.F90**  code shared by all restart options
|        **icedrv_state.F90**       essential arrays to describe the state of the ice
|        **icedrv_step.F90**        routines for time stepping the major code components
|        **icedrv_system.F90**      overall system management calls

Overview
------------

The icepack driver exists to test the column physics.  At the present time, it is hardwired
to run 4 different gridcells on one processor with the same forcing used for all gridcells.  
There is no MPI and no threading built into the icepack driver.  There is limited IO capabilities,
no history files, and no netcdf restart files.  The model generally runs very quickly.

Forcing data and details on these data are available in :ref:`force`.
