!  SVN:$Id: CICE_RunMod.F90 746 2013-09-28 22:47:56Z eclare $
!=======================================================================
!
!  Main driver for time stepping of CICE.
!
!  authors Elizabeth C. Hunke, LANL
!          Philip W. Jones, LANL
!          William H. Lipscomb, LANL
!
! 2006 ECH: moved exit timeLoop to prevent execution of unnecessary timestep
! 2006 ECH: Streamlined for efficiency
! 2006 ECH: Converted to free source form (F90)
! 2007 BPB: Modified Delta-Eddington shortwave interface
! 2008 ECH: moved ESMF code to its own driver

      module CICE_RunMod

      use ice_kinds_mod
      use perf_mod,        only : t_startf, t_stopf, t_barrierf
      use ice_fileunits, only: nu_diag

      implicit none
      private
      public :: CICE_Run, ice_step

!=======================================================================

      contains

!=======================================================================
!
!  This is the main driver routine for advancing CICE forward in time.
!
!  author Elizabeth C. Hunke, LANL
!         Philip W. Jones, LANL
!         William H. Lipscomb, LANL

      subroutine CICE_Run

      use ice_aerosol, only: faero_default
      use ice_algae, only: get_forcing_bgc
      use ice_calendar, only: istep, istep1, time, dt, stop_now, calendar
      use ice_forcing, only: get_forcing_atmo, get_forcing_ocn
      use ice_flux, only: init_flux_atm, init_flux_ocn
      use ice_state, only: tr_aero
      use ice_timers, only: ice_timer_start, ice_timer_stop, &
          timer_couple, timer_step
      use ice_zbgc_shared, only: skl_bgc

   !--------------------------------------------------------------------
   !  initialize error code and step timer
   !--------------------------------------------------------------------

      call ice_timer_start(timer_step)   ! start timing entire run

   !--------------------------------------------------------------------
   ! timestep loop
   !--------------------------------------------------------------------

   !      timeLoop: do

         istep  = istep  + 1    ! update time step counters
         istep1 = istep1 + 1
         time = time + dt       ! determine the time and date

         call init_flux_atm     ! initialize atmosphere fluxes sent to coupler
         call init_flux_ocn     ! initialize ocean fluxes sent to coupler

         call calendar(time)    ! at the end of the timestep

         call ice_step

!         if (stop_now >= 1) exit timeLoop

!      enddo timeLoop

   !--------------------------------------------------------------------
   ! end of timestep loop
   !--------------------------------------------------------------------

      call ice_timer_stop(timer_step)   ! end timestepping loop timer

      end subroutine CICE_Run

!=======================================================================
!
!  Calls drivers for physics components, some initialization, and output
!
!  author Elizabeth C. Hunke, LANL
!         William H. Lipscomb, LANL

      subroutine ice_step

      use ice_age, only: write_restart_age
      use ice_aerosol, only: write_restart_aero
      use ice_boundary, only: ice_HaloUpdate
      use ice_brine, only: hbrine_diags, write_restart_hbrine
      use ice_calendar, only: dt, dt_dyn, ndtd, diagfreq, write_restart, istep, idate, sec
      use ice_constants, only: field_loc_center, field_type_scalar
      use ice_diagnostics, only: init_mass_diags, runtime_diags, print_points_state
      use ice_domain, only: halo_info, nblocks
      use ice_domain_size, only: nslyr
      use ice_dyn_eap, only: write_restart_eap
      use ice_dyn_shared, only: kdyn
      use ice_firstyear, only: write_restart_FY
      use ice_flux, only: scale_factor, init_history_therm
      use ice_history, only: accum_hist
      use ice_lvl, only: write_restart_lvl
      use ice_restart, only: final_restart
      use ice_restart_driver, only: dumpfile
      use ice_meltpond_cesm, only: write_restart_pond_cesm
      use ice_meltpond_lvl, only: write_restart_pond_lvl
      use ice_meltpond_topo, only: write_restart_pond_topo
      use ice_restoring, only: restore_ice, ice_HaloRestore
      use ice_state, only: nt_qsno, trcrn, tr_iage, tr_FY, tr_lvl, &
          tr_pond_cesm, tr_pond_lvl, tr_pond_topo, tr_brine, tr_aero
      use ice_step_mod, only: prep_radiation, step_therm1, step_therm2, &
          post_thermo, step_dynamics, step_radiation
      use ice_therm_shared, only: calc_Tsfc
      use ice_timers, only: ice_timer_start, ice_timer_stop, &
          timer_diags, timer_column, timer_thermo, timer_bound, &
          timer_hist, timer_readwrite
      use ice_algae, only: bgc_diags, write_restart_bgc
      use ice_zbgc, only: init_history_bgc, biogeochemistry
      use ice_zbgc_shared, only: skl_bgc
      use ice_communicate, only: MPI_COMM_ICE
      use ice_prescribed_mod

      integer (kind=int_kind) :: &
         iblk        , & ! block index
         k               ! dynamics supercycling index

      !-----------------------------------------------------------------
      ! restoring on grid boundaries
      !-----------------------------------------------------------------

         if (restore_ice) call ice_HaloRestore

      !-----------------------------------------------------------------
      ! initialize diagnostics
      !-----------------------------------------------------------------

         call ice_timer_start(timer_diags)  ! diagnostics/history
         call init_mass_diags   ! diagnostics per timestep
         call init_history_therm
         call init_history_bgc
         call ice_timer_stop(timer_diags)   ! diagnostics/history

         if(prescribed_ice) then  ! read prescribed ice
            call t_barrierf('cice_run_presc_BARRIER',MPI_COMM_ICE)
            call t_startf ('cice_run_presc')
            call ice_prescribed_run(idate, sec)
            call t_stopf ('cice_run_presc')
         endif

         call ice_timer_start(timer_column)  ! column physics
         call ice_timer_start(timer_thermo)  ! thermodynamics

         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! Scale radiation fields
      !-----------------------------------------------------------------

            if (calc_Tsfc) call prep_radiation (dt, iblk)

      !-----------------------------------------------------------------
      ! thermodynamics
      !-----------------------------------------------------------------

            call step_therm1     (dt, iblk) ! vertical thermodynamics
            call biogeochemistry (dt, iblk) ! biogeochemistry
            if (.not.prescribed_ice) &
              call step_therm2   (dt, iblk) ! ice thickness distribution thermo

         enddo ! iblk
         !$OMP END PARALLEL DO

         call post_thermo (dt)             ! finalize thermo update

         call ice_timer_stop(timer_thermo) ! thermodynamics
         call ice_timer_stop(timer_column) ! column physics

      !-----------------------------------------------------------------
      ! dynamics, transport, ridging
      !-----------------------------------------------------------------

         if (.not.prescribed_ice .and. kdyn>0) then
         do k = 1, ndtd
            call step_dynamics (dt_dyn, ndtd)
         enddo
         endif

      !-----------------------------------------------------------------
      ! albedo, shortwave radiation
      !-----------------------------------------------------------------

         call ice_timer_start(timer_column)  ! column physics
         call ice_timer_start(timer_thermo)  ! thermodynamics

         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks

            call step_radiation (dt, iblk)

      !-----------------------------------------------------------------
      ! get ready for coupling and the next time step
      !-----------------------------------------------------------------

            call coupling_prep (iblk)

         enddo ! iblk
         !$OMP END PARALLEL DO

         call ice_timer_start(timer_bound)
         call ice_HaloUpdate (scale_factor,     halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_timer_stop(timer_bound)

         call ice_timer_stop(timer_thermo) ! thermodynamics
         call ice_timer_stop(timer_column) ! column physics

      !-----------------------------------------------------------------
      ! write data
      !-----------------------------------------------------------------

         call ice_timer_start(timer_diags)  ! diagnostics
         if (mod(istep,diagfreq) == 0) then
            call runtime_diags(dt)          ! log file
            if (skl_bgc)  call bgc_diags (dt)
            if (tr_brine) call hbrine_diags (dt)
         endif
         call ice_timer_stop(timer_diags)   ! diagnostics

         call ice_timer_start(timer_hist)   ! history
         call accum_hist (dt)               ! history file
         call ice_timer_stop(timer_hist)    ! history

         call ice_timer_start(timer_readwrite)  ! reading/writing
         if (write_restart == 1) then
            call dumpfile     ! core variables for restarting
            if (tr_iage)      call write_restart_age
            if (tr_FY)        call write_restart_FY
            if (tr_lvl)       call write_restart_lvl
            if (tr_pond_cesm) call write_restart_pond_cesm
            if (tr_pond_lvl)  call write_restart_pond_lvl
            if (tr_pond_topo) call write_restart_pond_topo
            if (tr_aero)      call write_restart_aero
            if (skl_bgc)      call write_restart_bgc
            if (tr_brine)     call write_restart_hbrine
            if (kdyn == 2)    call write_restart_eap
            call final_restart
         endif

         call ice_timer_stop(timer_readwrite)  ! reading/writing

      end subroutine ice_step

!=======================================================================
!
! Prepare for coupling
!
! authors: Elizabeth C. Hunke, LANL

      subroutine coupling_prep (iblk)

      use ice_blocks, only: block, nx_block, ny_block
      use ice_calendar, only: dt, nstreams
      use ice_constants, only: c0, c1, puny, rhofresh
      use ice_domain_size, only: ncat
      use ice_flux, only: alvdf, alidf, alvdr, alidr, albice, albsno, &
          albpnd, albcnt, apeff_ai, coszen, fpond, fresh, l_mpond_fresh, &
          alvdf_ai, alidf_ai, alvdr_ai, alidr_ai, fhocn_ai, &
          fresh_ai, fsalt_ai, fsalt, &
          fswthru_ai, fhocn, fswthru, scale_factor, &
          fswthruvdr, fswthruvdf, fswthruidr, fswthruidf, &
          swvdr, swidr, swvdf, swidf, Tf, Tair, Qa, strairxT, strairyt, &
          fsens, flat, fswabs, flwout, evap, Tref, Qref, Uref, faero_ocn, &
          fsurfn_f, flatn_f, scale_fluxes, frzmlt_init, frzmlt, wind, &
          snowfrac
      use ice_grid, only: tmask
      use ice_ocean, only: oceanmixed_ice, ocean_mixed_layer
      use ice_shortwave, only: alvdfn, alidfn, alvdrn, alidrn, &
                               albicen, albsnon, albpndn, apeffn, snowfracn
      use ice_state, only: aicen, aice, aice_init, nbtrcr
      use ice_therm_shared, only: calc_Tsfc
      use ice_timers, only: timer_couple, ice_timer_start, ice_timer_stop
      use ice_zbgc_shared, only: flux_bio, flux_bio_ai

      integer (kind=int_kind), intent(in) :: &
         iblk            ! block index

      ! local variables

      integer (kind=int_kind) :: &
         n           , & ! thickness category index
         i,j         , & ! horizontal indices
         k               ! tracer index

      real (kind=dbl_kind) :: cszn ! counter for history averaging

      real (kind=dbl_kind) :: netsw

      !-----------------------------------------------------------------
      ! Save current value of frzmlt for diagnostics.
      ! Update mixed layer with heat and radiation from ice.
      !-----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            frzmlt_init  (i,j,iblk) = frzmlt(i,j,iblk)
         enddo
         enddo

         call ice_timer_start(timer_couple,iblk)   ! atm/ocn coupling

         if (oceanmixed_ice) &
         call ocean_mixed_layer (dt,iblk) ! ocean surface fluxes and sst

      !-----------------------------------------------------------------
      ! Aggregate albedos
      !-----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            alvdf(i,j,iblk) = c0
            alidf(i,j,iblk) = c0
            alvdr(i,j,iblk) = c0
            alidr(i,j,iblk) = c0

            albice(i,j,iblk) = c0
            albsno(i,j,iblk) = c0
            albpnd(i,j,iblk) = c0
            apeff_ai(i,j,iblk) = c0
            snowfrac(i,j,iblk) = c0

            ! for history averaging
            cszn = c0
            netsw = swvdr(i,j,iblk)+swidr(i,j,iblk)+swvdf(i,j,iblk)+swidf(i,j,iblk)
            if (netsw > puny) cszn = c1
            do n = 1, nstreams
               albcnt(i,j,iblk,n) = albcnt(i,j,iblk,n) + cszn
            enddo
         enddo
         enddo
         do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            alvdf(i,j,iblk) = alvdf(i,j,iblk) &
               + alvdfn(i,j,n,iblk)*aicen(i,j,n,iblk)
            alidf(i,j,iblk) = alidf(i,j,iblk) &
               + alidfn(i,j,n,iblk)*aicen(i,j,n,iblk)
            alvdr(i,j,iblk) = alvdr(i,j,iblk) &
               + alvdrn(i,j,n,iblk)*aicen(i,j,n,iblk)
            alidr(i,j,iblk) = alidr(i,j,iblk) &
               + alidrn(i,j,n,iblk)*aicen(i,j,n,iblk)

            if (coszen(i,j,iblk) > puny) then ! sun above horizon
            albice(i,j,iblk) = albice(i,j,iblk) &
               + albicen(i,j,n,iblk)*aicen(i,j,n,iblk)
            albsno(i,j,iblk) = albsno(i,j,iblk) &
               + albsnon(i,j,n,iblk)*aicen(i,j,n,iblk)
            albpnd(i,j,iblk) = albpnd(i,j,iblk) &
               + albpndn(i,j,n,iblk)*aicen(i,j,n,iblk)
            endif

            apeff_ai(i,j,iblk) = apeff_ai(i,j,iblk) &       ! for history
               + apeffn(i,j,n,iblk)*aicen(i,j,n,iblk)
            snowfrac(i,j,iblk) = snowfrac(i,j,iblk) &       ! for history
               + snowfracn(i,j,n,iblk)*aicen(i,j,n,iblk)
         enddo
         enddo
         enddo

         do j = 1, ny_block
         do i = 1, nx_block

      !-----------------------------------------------------------------
      ! reduce fresh by fpond for coupling
      !-----------------------------------------------------------------

            if (l_mpond_fresh) then
               fpond(i,j,iblk) = fpond(i,j,iblk) * rhofresh/dt
               fresh(i,j,iblk) = fresh(i,j,iblk) - fpond(i,j,iblk)
            endif

      !----------------------------------------------------------------
      ! Store grid box mean albedos and fluxes before scaling by aice
      !----------------------------------------------------------------

            alvdf_ai  (i,j,iblk) = alvdf  (i,j,iblk)
            alidf_ai  (i,j,iblk) = alidf  (i,j,iblk)
            alvdr_ai  (i,j,iblk) = alvdr  (i,j,iblk)
            alidr_ai  (i,j,iblk) = alidr  (i,j,iblk)
            fresh_ai  (i,j,iblk) = fresh  (i,j,iblk)
            fsalt_ai  (i,j,iblk) = fsalt  (i,j,iblk)
            fhocn_ai  (i,j,iblk) = fhocn  (i,j,iblk)
            fswthru_ai(i,j,iblk) = fswthru(i,j,iblk)

            if (nbtrcr > 0) then
            do k = 1, nbtrcr
              flux_bio_ai  (i,j,k,iblk) = flux_bio  (i,j,k,iblk)
            enddo
            endif

      !-----------------------------------------------------------------
      ! Save net shortwave for scaling factor in scale_factor
      !-----------------------------------------------------------------
            scale_factor(i,j,iblk) = &
                       swvdr(i,j,iblk)*(c1 - alvdr_ai(i,j,iblk)) &
                     + swvdf(i,j,iblk)*(c1 - alvdf_ai(i,j,iblk)) &
                     + swidr(i,j,iblk)*(c1 - alidr_ai(i,j,iblk)) &
                     + swidf(i,j,iblk)*(c1 - alidf_ai(i,j,iblk))

         enddo
         enddo

      !-----------------------------------------------------------------
      ! Divide fluxes by ice area
      !  - the CCSM coupler assumes fluxes are per unit ice area
      !  - also needed for global budget in diagnostics
      !-----------------------------------------------------------------

         call scale_fluxes (nx_block,              ny_block,             &
                            tmask     (:,:,iblk)  , nbtrcr,              &
                            aice      (:,:,iblk)  , Tf      (:,:,iblk),  &
                            Tair      (:,:,iblk)  , Qa      (:,:,iblk),  &
                            strairxT  (:,:,iblk)  , strairyT(:,:,iblk),  &
                            fsens     (:,:,iblk)  , flat    (:,:,iblk),  &
                            fswabs    (:,:,iblk)  , flwout  (:,:,iblk),  &
                            evap      (:,:,iblk)  ,                      &
                            Tref      (:,:,iblk)  , Qref    (:,:,iblk),  &
                            fresh     (:,:,iblk)  , fsalt   (:,:,iblk),  &
                            fhocn     (:,:,iblk)  , fswthru (:,:,iblk),  &
                            fswthruvdr(:,:,iblk),  fswthruvdf(:,:,iblk), &
                            fswthruidr(:,:,iblk) , fswthruidf(:,:,iblk), &
                            faero_ocn (:,:,:,iblk),                      &
                            alvdr     (:,:,iblk)  , alidr   (:,:,iblk),  &
                            alvdf     (:,:,iblk)  , alidf   (:,:,iblk),  &
                            flux_bio  (:,:,1:nbtrcr,iblk),               &
                            Uref=Uref (:,:,iblk), wind=wind(:,:,iblk) )

!echmod - comment this out for efficiency, if .not. calc_Tsfc
         if (.not. calc_Tsfc) then

       !---------------------------------------------------------------
       ! If surface fluxes were provided, conserve these fluxes at ice
       ! free points by passing to ocean.
       !---------------------------------------------------------------

            call sfcflux_to_ocn &
                         (nx_block,              ny_block,             &
                          tmask   (:,:,iblk),    aice_init(:,:,iblk),  &
                          fsurfn_f (:,:,:,iblk), flatn_f(:,:,:,iblk),  &
                          fresh    (:,:,iblk),   fhocn    (:,:,iblk))
         endif
!echmod

         call ice_timer_stop(timer_couple,iblk)   ! atm/ocn coupling

      end subroutine coupling_prep

!=======================================================================
!
! If surface heat fluxes are provided to CICE instead of CICE calculating
! them internally (i.e. .not. calc_Tsfc), then these heat fluxes can
! be provided at points which do not have ice.  (This is could be due to
! the heat fluxes being calculated on a lower resolution grid or the
! heat fluxes not recalculated at every CICE timestep.)  At ice free points,
! conserve energy and water by passing these fluxes to the ocean.
!
! author: A. McLaren, Met Office

      subroutine sfcflux_to_ocn(nx_block,   ny_block,     &
                                tmask,      aice,         &
                                fsurfn_f,   flatn_f,      &
                                fresh,      fhocn)

      use ice_domain_size, only: ncat

      integer (kind=int_kind), intent(in) :: &
          nx_block, ny_block  ! block dimensions

      logical (kind=log_kind), dimension (nx_block,ny_block), &
          intent(in) :: &
          tmask       ! land/boundary mask, thickness (T-cell)

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(in):: &
          aice        ! initial ice concentration

      real (kind=dbl_kind), dimension(nx_block,ny_block,ncat), &
          intent(in) :: &
          fsurfn_f, & ! net surface heat flux (provided as forcing)
          flatn_f     ! latent heat flux (provided as forcing)

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(inout):: &
          fresh        , & ! fresh water flux to ocean         (kg/m2/s)
          fhocn            ! actual ocn/ice heat flx           (W/m**2)

#ifdef CICE_IN_NEMO

      ! local variables
      integer (kind=int_kind) :: &
          i, j, n    ! horizontal indices

      real (kind=dbl_kind)    :: &
          rLsub            ! 1/Lsub

      rLsub = c1 / Lsub

      do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            if (tmask(i,j) .and. aice(i,j) <= puny) then
               fhocn(i,j)      = fhocn(i,j)              &
                            + fsurfn_f(i,j,n) + flatn_f(i,j,n)
               fresh(i,j)      = fresh(i,j)              &
                                 + flatn_f(i,j,n) * rLsub
            endif
         enddo   ! i
         enddo   ! j
      enddo      ! n

#endif

      end subroutine sfcflux_to_ocn

!=======================================================================

      end module CICE_RunMod

!=======================================================================
