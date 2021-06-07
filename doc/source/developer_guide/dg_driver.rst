:tocdepth: 3

.. _dev_driver:


Driver and Coupling 
====================

The driver and coupling layer is found in **cicecore/drivers/**.  The standalone driver is found
under **cicecore/drivers/standalone/cice/** and other high level coupling layers are found in other directories.
CICE is designed to build with only one of these drivers at a time, depending how the model is run and coupled.  Within the **cicecore/drivers/standalone/cice/** directory, the following files are found,

**CICE.F90** is the top level program file and that calls CICE_Initialize, CICE_Run, and CICE_Finalize methods.
**CICE_InitMod.F90** contains the CICE_Initialize method and other next level source code.
**CICE_RunMod.F90** contains the CICE_Run method and other next level source code.
**CICE_FinalMod.F90** contains the CICE_Finalize method and other next level source code.

The files provide the top level sequencing for calling the standalone CICE model.

Adding a New Driver
------------------------

The drivers directory contains two levels of subdirectories.  The first layer indicates the coupling infrastructure or strategy and the second later indicates the application or coupler the driver is written for.  At the present time, the directory structures is::

  drivers/direct/hadgem3
  drivers/mct/cesm1
  drivers/nuopc/cmeps
  drivers/standalone/cice

The standalone driver is **drivers/standalone/cice**, and this is the driver used when running with the CICE scripts in standalone mode.  New drivers can be added as needed when coupling to new infrastructure or in new applications.  We encourage the community to use the drivers directory to facilitate reuse with the understanding that the driver code could also reside in the application.  Users should follow the naming strategy as best as possible. Drivers should be added under the appropriate subdirectory indicative of the coupling infrastructure.  New subdirectories (such as oasis or esmf) can be added in the future as needed.  The community will have to decide when it's appropriate to share drivers between different applications, when to update drivers, and when to create new drivers.  There are a number of trade-offs to consider including backwards compatibility with earlier versions of applications, code reuse, and independence.  As a general rule, driver directories should not be deleted and names should not be reused to avoid confusion with prior versions that were fundamentally different.  The number of drivers will likely increase over time as new infrastructure and applications are added and as versions evolve in time.

The current drivers subdirectories are mct, nuopc, standalone, and direct.  The standalone subdirectory contains drivers to run the model in standalone mode as a standalone program.  The direct subdirectory contains coupling interfaces that supporting calling the ice model directory from other models as subroutines.  The subdirectory mct contains subdirectories for applications/couplers that provide coupling via mct interfaces.  And the subdirectory nuopc contains subdirectories for applications/couplers that provide coupling via nuopc interfaces.

The varied **cicecore/drivers/** directories are generally implemented similar to the standalone cice case with versions of **CICE_InitMod.F90**, **CICE_RunMod.F90**, and **CICE_FinalMod.F90** files in addition to files consistent with the coupling layer.

As features are added to the CICE model over time that require changes in the calling sequence, it's possible that all drivers will need to be updated.  These kinds of changes are impactful and not taken lightly.  It will be up to the community as a whole to work together to maintain the various drivers in these situations.


Calling Sequence
------------------------

The initialize calling sequence looks something like::

      call init_communicate     ! initial setup for message passing
      call init_fileunits       ! unit numbers
      call icepack_configure()  ! initialize icepack
      call input_data           ! namelist variables
      call init_zbgc            ! vertical biogeochemistry namelist
      call count_tracers        ! count tracers
      call init_domain_blocks   ! set up block decomposition
      call init_grid1           ! domain distribution
      call alloc_*              ! allocate arrays
      call init_ice_timers      ! initialize all timers
      call init_grid2           ! grid variables
      call init_zbgc            ! vertical biogeochemistry initialization
      call init_calendar        ! initialize some calendar stuff
      call init_hist (dt)       ! initialize output history file
      call init_dyn (dt_dyn)    ! define dynamics parameters, variables
      if (kdyn == 2) then
         call init_eap          ! define eap dynamics parameters, variables
      else if (kdyn == 3) then
         call init_vp           ! define vp dynamics parameters, variables
      endif
      call init_coupler_flux    ! initialize fluxes exchanged with coupler
      call init_thermo_vertical ! initialize vertical thermodynamics
      call icepack_init_itd(ncat, hin_max)  ! ice thickness distribution
      if (tr_fsd) call icepack_init_fsd_bounds  ! floe size distribution
      call calendar(time)       ! determine the initial date
      call init_forcing_ocn(dt) ! initialize sss and sst from data
      call init_state           ! initialize the ice state
      call init_transport       ! initialize horizontal transport
      call ice_HaloRestore_init ! restored boundary conditions
      call init_restart         ! initialize restart variables
      call init_diags           ! initialize diagnostic output points
      call init_history_therm   ! initialize thermo history variables
      call init_history_dyn     ! initialize dynamic history variables
      call init_shortwave       ! initialize radiative transfer
      call init_forcing_atmo    ! initialize atmospheric forcing (standalone)
      if (tr_fsd .and. wave_spec) call get_wave_spec ! wave spectrum in ice
      call get_forcing*         ! read forcing data (standalone)

See a **CICE_InitMod.F90** file for the latest.

The run sequence within a time loop looks something like::

         call init_mass_diags   ! diagnostics per timestep
         call init_history_therm
         call init_history_bgc

         do iblk = 1, nblocks
            if (calc_Tsfc) call prep_radiation (dt, iblk)
            call step_therm1     (dt, iblk) ! vertical thermodynamics
            call biogeochemistry (dt, iblk) ! biogeochemistry
            call step_therm2     (dt, iblk) ! ice thickness distribution thermo
         enddo ! iblk

         call update_state (dt, daidtt, dvidtt, dagedtt, offset)

         if (tr_fsd .and. wave_spec) call step_dyn_wave(dt)
         do k = 1, ndtd
            call step_dyn_horiz (dt_dyn)
            do iblk = 1, nblocks
               call step_dyn_ridge (dt_dyn, ndtd, iblk)
            enddo
            call update_state (dt_dyn, daidtd, dvidtd, dagedtd, offset)
         enddo

         do iblk = 1, nblocks
            call step_radiation (dt, iblk)
            call coupling_prep (iblk)
         enddo ! iblk

         ! write data
         ! update forcing

See a **CICE_RunMod.F90** file for the latest.
