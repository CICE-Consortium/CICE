:tocdepth: 3

.. _dev_driver:


Driver and Coupling 
====================

The driver and coupling layer is found in **cicecore/drivers/**.  The standalone driver is found
under **cicecore/drivers/cice/** and other high level coupling layers are found in other directories.
In general, CICE will build with only one of these drivers, depending how the model is run and
coupled.  Within the **cicecore/drivers/cice/** directory, the following files are found,

**CICE.F90** is the top level program file and that calls CICE_Initialize, CICE_Run, and CICE_Finalize methods.
**CICE_InitMod.F90** contains the CICE_Initialize method and other next level source code.
**CICE_RunMod.F90** contains the CICE_Run method and other next level source code.
**CICE_FinalMod.F90** contains the CICE_Finalize method and other next level source code.

Other **cicecore/drivers/** directories are similarly implemented with a top level coupling layer,
that is largely specified by an external coupled system and then some version of the **CICE_InitMod.F90**,
**CICE_RunMod.F90**, and **CICE_FinalMod.F90** files.


Calling Sequence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The initialize calling sequence looks something like::

      call init_communicate     ! initial setup for message passing
      call init_fileunits       ! unit numbers
      call icepack_configure()  ! initialize icepack
      call input_data           ! namelist variables
      call init_zbgc            ! vertical biogeochemistry namelist
      call init_domain_blocks   ! set up block decomposition
      call init_grid1           ! domain distribution
      call alloc_*              ! allocate arrays
      call init_ice_timers      ! initialize all timers
      call init_grid2           ! grid variables
      call init_calendar        ! initialize some calendar stuff
      call init_hist (dt)       ! initialize output history file
      if (kdyn == 2) then
         call init_eap (dt_dyn) ! define eap dynamics parameters, variables
      else                      ! for both kdyn = 0 or 1
         call init_evp (dt_dyn) ! define evp dynamics parameters, variables
      endif
      call init_coupler_flux    ! initialize fluxes exchanged with coupler
      call init_thermo_vertical ! initialize vertical thermodynamics
      call icepack_init_itd(ncat, hin_max)  ! ice thickness distribution
      call calendar(time)       ! determine the initial date
      call init_forcing_ocn(dt) ! initialize sss and sst from data
      call init_state           ! initialize the ice state
      call init_transport       ! initialize horizontal transport
      call ice_HaloRestore_init ! restored boundary conditions
      call init_restart         ! initialize restart variables
      call init_diags           ! initialize diagnostic output points
      call init_history_therm   ! initialize thermo history variables
      call init_history_dyn     ! initialize dynamic history variables
      call init_shortwave    ! initialize radiative transfer
      call init_forcing_atmo    ! initialize atmospheric forcing (standalone)

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

See a **CICE_RunMod.F90** file for the latest.
