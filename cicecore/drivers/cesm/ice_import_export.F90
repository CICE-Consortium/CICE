module ice_import_export

  use shr_kind_mod      , only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use shr_sys_mod       , only: shr_sys_abort, shr_sys_flush
  use shr_mpi_mod       , only: shr_mpi_max, shr_mpi_sum
  use ice_kinds_mod     , only: int_kind, dbl_kind, char_len, char_len_long, log_kind
  use ice_constants     , only: c0, c1, spval_dbl
  use ice_constants     , only: field_loc_center, field_type_scalar
  use ice_constants     , only: field_type_vector, c100
  use ice_constants     , only: p001, p5
  use ice_blocks        , only: block, get_block, nx_block, ny_block
  use ice_flux          , only: strairxt, strairyt, strocnxt, strocnyt           
  use ice_flux          , only: alvdr, alidr, alvdf, alidf, Tref, Qref, Uref
  use ice_flux          , only: flat, fsens, flwout, evap, fswabs, fhocn, fswthru
  use ice_flux          , only: fresh, fsalt, zlvl, uatm, vatm, potT, Tair, Qa
  use ice_flux          , only: rhoa, swvdr, swvdf, swidr, swidf, flw, frain
  use ice_flux          , only: fsnow, uocn, vocn, sst, ss_tltx, ss_tlty, frzmlt
  use ice_flux          , only: sss, tf, wind, fsw, init_flux_atm, init_flux_ocn
  use ice_flux_bgc      , only: faero_atm, flux_bio_atm, nit, amm, sil, dmsp, &
                          dms, hum, algalN, doc, don, dic, fed, fep, zaeros, &
                          falgalN, fdoc, fdon, fnit, famm, fsil, fdmsp, fdms, &
                          fhum, ffed, ffep, fdust, cpl_bgc
  use ice_exit, only: abort_ice
  use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
  use icepack_intfc, only: icepack_query_parameters, icepack_query_tracer_flags
  use ice_arrays_column , only: Cdn_atm, R_C2N
  use ice_state         , only: vice, vsno, aice, trcr

  use ice_domain        , only: nblocks, blocks_ice, halo_info, distrb_info
  use ice_domain_size   , only: nx_global, ny_global, block_size_x, block_size_y, max_blocks
  use ice_grid          , only: tlon, tlat, tarea, tmask, anglet, hm
  use ice_grid          , only: grid_type, t2ugrid_vector
  use ice_boundary      , only: ice_HaloUpdate 
  use ice_communicate   , only: my_task, master_task, MPI_COMM_ICE, get_num_procs
  use ice_calendar      , only: istep, istep1, diagfreq
  use ice_fileunits     , only: nu_diag
  use ice_prescribed_mod
  use ice_cpl_indices
  use perf_mod          , only: t_startf, t_stopf, t_barrierf

  implicit none
  public

#ifdef RASM_MODS
! (1)  Andrew Roberts:  Added artificial correction to snow and rain division
!      This is to be consistent with VIC in the Regional Arctic System Model
  logical, parameter :: rasm_snowrain_split = .true.
#else
  logical, parameter :: rasm_snowrain_split = .false.
#endif

  !==============================================================================
contains
  !==============================================================================

  subroutine ice_import( x2i )

    !-----------------------------------------------------
    ! Arguments
    real(r8), intent(inout) :: x2i(:,:)
    !
    ! Local variables
    integer     :: i, j, iblk, n
    integer     :: ilo, ihi, jlo, jhi !beginning and end of physical domain
    type(block) :: this_block         ! block information for current block
    integer,parameter       :: nflds=17,nfldv=6,nfldb=27
    real (kind=dbl_kind),allocatable :: aflds(:,:,:,:)
    real (kind=dbl_kind)    :: workx, worky
    real (kind=dbl_kind)    :: MIN_RAIN_TEMP, MAX_SNOW_TEMP 
    character(len=char_len) :: tfrz_option
    logical (kind=log_kind) :: modal_aero, z_tracers, skl_bgc
    logical (kind=log_kind) :: tr_aero, tr_iage, tr_FY, tr_pond
    logical (kind=log_kind) :: tr_lvl, tr_zaero, tr_bgc_Nit 
    real (kind=dbl_kind)    :: tffresh
    logical (kind=log_kind) :: first_call = .true.
    character(len=*), parameter :: subname = '(ice_import)'
    !-----------------------------------------------------

    call icepack_query_parameters(tfrz_option_out=tfrz_option, &
       modal_aero_out=modal_aero, z_tracers_out=z_tracers, skl_bgc_out=skl_bgc, &
       Tffresh_out=Tffresh)
    call icepack_query_tracer_flags(tr_aero_out=tr_aero, tr_iage_out=tr_iage, &
       tr_FY_out=tr_FY, tr_pond_out=tr_pond, tr_lvl_out=tr_lvl, &
       tr_zaero_out=tr_zaero, tr_bgc_Nit_out=tr_bgc_Nit)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
        file=__FILE__, line=__LINE__)

    ! Note that the precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2 which is what CICE requires.
    ! Note also that the read in below includes only values needed
    ! by the thermodynamic component of CICE.  Variables uocn, vocn,
    ! ss_tltx, and ss_tlty are excluded. Also, because the SOM and
    ! DOM don't  compute SSS.   SSS is not read in and is left at
    ! the initilized value (see ice_flux.F init_coupler_flux) of
    ! 34 ppt

    ! Use aflds to gather the halo updates of multiple fields
    ! Need to separate the scalar from the vector halo updates

    allocate(aflds(nx_block,ny_block,nflds,nblocks))
    aflds = c0

    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo, jhi
          do i = ilo, ihi

             n = n+1
             aflds(i,j, 1,iblk)   = x2i(index_x2i_So_t,n)
             aflds(i,j, 2,iblk)   = x2i(index_x2i_So_s,n)
             aflds(i,j, 3,iblk)   = x2i(index_x2i_Sa_z,n)
             aflds(i,j, 4,iblk)   = x2i(index_x2i_Sa_ptem,n)
             aflds(i,j, 5,iblk)   = x2i(index_x2i_Sa_tbot,n)
             aflds(i,j, 6,iblk)   = x2i(index_x2i_Sa_shum,n)
             aflds(i,j, 7,iblk)   = x2i(index_x2i_Sa_dens,n)
             aflds(i,j, 8,iblk)   = x2i(index_x2i_Fioo_q,n)
             aflds(i,j, 9,iblk)   = x2i(index_x2i_Faxa_swvdr,n)
             aflds(i,j,10,iblk)   = x2i(index_x2i_Faxa_swndr,n)
             aflds(i,j,11,iblk)   = x2i(index_x2i_Faxa_swvdf,n)
             aflds(i,j,12,iblk)   = x2i(index_x2i_Faxa_swndf,n)
             aflds(i,j,13,iblk)   = x2i(index_x2i_Faxa_lwdn,n)
             aflds(i,j,14,iblk)   = x2i(index_x2i_Faxa_rain,n)
             aflds(i,j,15,iblk)   = x2i(index_x2i_Faxa_snow,n)
             if (index_x2i_So_z > 0) &
                aflds(i,j,16,iblk)   = x2i(index_x2i_So_z,n)
             if (index_x2i_So_logz0 > 0) &
                aflds(i,j,17,iblk)   = x2i(index_x2i_So_logz0,n)

          enddo    !i
       enddo    !j

    enddo        !iblk

    if (.not.prescribed_ice) then
       call t_startf ('cice_imp_halo')
       call ice_HaloUpdate(aflds, halo_info, field_loc_center, &
            field_type_scalar)
       call t_stopf ('cice_imp_halo')
    endif

    if (rasm_snowrain_split) then
       MIN_RAIN_TEMP = Tffresh-c1
       MAX_SNOW_TEMP = Tffresh+c0
    endif

    !$OMP PARALLEL DO PRIVATE(iblk,i,j)
    do iblk = 1, nblocks
       do j = 1,ny_block
          do i = 1,nx_block
             sst  (i,j,iblk)   = aflds(i,j, 1,iblk)
             sss  (i,j,iblk)   = aflds(i,j, 2,iblk)
             zlvl (i,j,iblk)   = aflds(i,j, 3,iblk)
             potT (i,j,iblk)   = aflds(i,j, 4,iblk)
             Tair (i,j,iblk)   = aflds(i,j, 5,iblk)
             Qa   (i,j,iblk)   = aflds(i,j, 6,iblk)
             rhoa (i,j,iblk)   = aflds(i,j, 7,iblk)
             frzmlt (i,j,iblk) = aflds(i,j, 8,iblk)
             swvdr(i,j,iblk)   = aflds(i,j, 9,iblk)
             swidr(i,j,iblk)   = aflds(i,j,10,iblk)
             swvdf(i,j,iblk)   = aflds(i,j,11,iblk)
             swidf(i,j,iblk)   = aflds(i,j,12,iblk)
             flw  (i,j,iblk)   = aflds(i,j,13,iblk)
             frain(i,j,iblk)   = aflds(i,j,14,iblk)
             fsnow(i,j,iblk)   = aflds(i,j,15,iblk)
!             if (index_x2i_So_z > 0) &
!                zlvlo(i,j,iblk) = aflds(i,j,16,iblk)
!
!             ! Convert log(z0) from the ocean to z0. log(z0) is transmitted
!             ! through the coupler in case any area waiting is done inside
!             ! the coupler.
!             if (index_x2i_So_logz0 > 0) &
!                z0ocn(i,j,iblk)   = exp(aflds(i,j,17,iblk))

             if (rasm_snowrain_split) then
                !--- Artificial correction to snow and rain for RASM
                if (Tair(i,j,iblk)<MIN_RAIN_TEMP) then
                   fsnow(i,j,iblk)=fsnow(i,j,iblk)+frain(i,j,iblk)
                   frain(i,j,iblk)=0
                elseif (Tair(i,j,iblk)>MAX_SNOW_TEMP) then
                   frain(i,j,iblk)=fsnow(i,j,iblk)+frain(i,j,iblk)
                   fsnow(i,j,iblk)=0
                else
                   frain(i,j,iblk)=fsnow(i,j,iblk)+frain(i,j,iblk)
                   fsnow(i,j,iblk)=frain(i,j,iblk)
                   frain(i,j,iblk)=frain(i,j,iblk)*(Tair(i,j,iblk)-MIN_RAIN_TEMP) / &
                                                    (MAX_SNOW_TEMP-MIN_RAIN_TEMP)
                   fsnow(i,j,iblk)=fsnow(i,j,iblk)-frain(i,j,iblk)
                endif
                !--- end artificial RASM correction
             endif  ! rasm_snowrain_split

          enddo    !i
       enddo    !j
    enddo        !iblk
    !$OMP END PARALLEL DO

    deallocate(aflds)
    allocate(aflds(nx_block,ny_block,nfldv,nblocks))
    aflds = c0

    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo, jhi
          do i = ilo, ihi
             n = n+1
             aflds(i,j, 1,iblk)   = x2i(index_x2i_So_u,n)
             aflds(i,j, 2,iblk)   = x2i(index_x2i_So_v,n)
             aflds(i,j, 3,iblk)   = x2i(index_x2i_Sa_u,n)
             aflds(i,j, 4,iblk)   = x2i(index_x2i_Sa_v,n)
             aflds(i,j, 5,iblk)   = x2i(index_x2i_So_dhdx,n)
             aflds(i,j, 6,iblk)   = x2i(index_x2i_So_dhdy,n)
          enddo
       enddo
    enddo

    if (.not.prescribed_ice) then
       call t_startf ('cice_imp_halo')
       call ice_HaloUpdate(aflds, halo_info, field_loc_center, &
            field_type_vector)
       call t_stopf ('cice_imp_halo')
    endif

    !$OMP PARALLEL DO PRIVATE(iblk,i,j)
    do iblk = 1, nblocks
       do j = 1,ny_block
          do i = 1,nx_block
             uocn (i,j,iblk)   = aflds(i,j, 1,iblk)
             vocn (i,j,iblk)   = aflds(i,j, 2,iblk)
             uatm (i,j,iblk)   = aflds(i,j, 3,iblk)
             vatm (i,j,iblk)   = aflds(i,j, 4,iblk)
             ss_tltx(i,j,iblk) = aflds(i,j, 5,iblk)
             ss_tlty(i,j,iblk) = aflds(i,j, 6,iblk)
          enddo    !i
       enddo    !j
    enddo        !iblk
    !$OMP END PARALLEL DO

    deallocate(aflds)

    !-------------------------------------------------------
    ! Set aerosols from coupler 
    !-------------------------------------------------------

    allocate(aflds(nx_block,ny_block,nfldb,nblocks))
    aflds = c0

    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo, jhi
          do i = ilo, ihi

             n = n+1
             if (tr_aero .or. tr_zaero) then
                if (modal_aero) then

                   ! BC species 1 (=intersitial/external BC)
                   aflds(i,j,1,iblk) = x2i(index_x2i_Faxa_bcphodry,n) &
                              + x2i(index_x2i_Faxa_bcphidry,n)

                   ! BC species 2 (=cloud_water/within-ice BC)
                   aflds(i,j,2,iblk) = x2i(index_x2i_Faxa_bcphiwet,n)

                   ! Combine all of the dust into one category
                   aflds(i,j,3,iblk) = x2i(index_x2i_Faxa_dstwet1,n) &
                              + x2i(index_x2i_Faxa_dstdry1,n) &
                              + x2i(index_x2i_Faxa_dstwet2,n) &
                              + x2i(index_x2i_Faxa_dstdry2,n) &
                              + x2i(index_x2i_Faxa_dstwet3,n) &
                              + x2i(index_x2i_Faxa_dstdry3,n) &
                              + x2i(index_x2i_Faxa_dstwet4,n) &
                              + x2i(index_x2i_Faxa_dstdry4,n)
                else
                   aflds(i,j,1,iblk) = x2i(index_x2i_Faxa_bcphodry,n)

                   aflds(i,j,2,iblk) = x2i(index_x2i_Faxa_bcphidry,n) &
                              + x2i(index_x2i_Faxa_bcphiwet,n)
                   ! Combine all of the dust into one category
                   aflds(i,j,3,iblk) = x2i(index_x2i_Faxa_dstwet1,n) &
                              + x2i(index_x2i_Faxa_dstdry1,n) &
                              + x2i(index_x2i_Faxa_dstwet2,n) &
                              + x2i(index_x2i_Faxa_dstdry2,n) &
                              + x2i(index_x2i_Faxa_dstwet3,n) &
                              + x2i(index_x2i_Faxa_dstdry3,n) &
                              + x2i(index_x2i_Faxa_dstwet4,n) &
                              + x2i(index_x2i_Faxa_dstdry4,n)
                endif
             endif
             if (cpl_bgc .or. (z_tracers .and. tr_bgc_Nit) .or. skl_bgc) then
                aflds(i,j,4,iblk)      = x2i(index_x2i_So_diat, n)
                aflds(i,j,5,iblk)      = x2i(index_x2i_So_sp, n)
                aflds(i,j,6,iblk)      = x2i(index_x2i_So_phaeo, n)
                aflds(i,j,7,iblk)      = x2i(index_x2i_So_doc, n) * p5 ! split evenly for now
                aflds(i,j,8,iblk)      = x2i(index_x2i_So_doc, n) * p5 !x2i(index_x2i_So_doc2, n)
                aflds(i,j,9,iblk)      = c0
                aflds(i,j,10,iblk)     = c0  !x2i(index_x2i_So_dic, n) 
                aflds(i,j,11,iblk)     = x2i(index_x2i_So_don, n)
                aflds(i,j,12,iblk)     = x2i(index_x2i_So_no3, n)
                aflds(i,j,13,iblk)     = x2i(index_x2i_So_sio3, n)
                aflds(i,j,14,iblk)     = x2i(index_x2i_So_nh4, n)
                aflds(i,j,15,iblk)     = x2i(index_x2i_So_dms, n)
                aflds(i,j,16,iblk)     = x2i(index_x2i_So_dmsp, n)
                aflds(i,j,17,iblk)     = x2i(index_x2i_So_donr, n)
                aflds(i,j,18,iblk)     = c0 !x2i(index_x2i_So_fep1, n)
                aflds(i,j,19,iblk)     = c0 !x2i(index_x2i_So_fep2, n)
                aflds(i,j,20,iblk)     = x2i(index_x2i_So_fed, n)
                aflds(i,j,21,iblk)     = c0 !x2i(index_x2i_So_fed2, n)
                aflds(i,j,22,iblk)     = c0 !x2i(index_x2i_So_zaer1, n) 
                aflds(i,j,23,iblk)     = c0 !x2i(index_x2i_So_zaer2, n) 
                aflds(i,j,24,iblk)     = c0 !x2i(index_x2i_So_zaer3, n) 
                aflds(i,j,25,iblk)     = c0 !x2i(index_x2i_So_zaer4, n) 
                aflds(i,j,26,iblk)     = c0 !x2i(index_x2i_So_zaer5, n) 
                aflds(i,j,27,iblk)     = c0 !x2i(index_x2i_So_zaer6, n) 
             endif
          enddo
       enddo
    enddo

    if (.not.prescribed_ice) then
       call t_startf ('cice_imp_halo')
       call ice_HaloUpdate(aflds, halo_info, field_loc_center, &
            field_type_scalar)
       call t_stopf ('cice_imp_halo')
    endif

    !$OMP PARALLEL DO PRIVATE(iblk,i,j)
    do iblk = 1, nblocks
       do j = 1,ny_block
          do i = 1,nx_block
             faero_atm(i,j,1,iblk) = aflds(i,j,1,iblk)
             faero_atm(i,j,2,iblk) = aflds(i,j,2,iblk)
             faero_atm(i,j,3,iblk) = aflds(i,j,3,iblk)    
          enddo    !i
       enddo    !j
    enddo        !iblk
    !$OMP END PARALLEL DO

    if (cpl_bgc) then
       !$OMP PARALLEL DO PRIVATE(iblk,i,j)
       do iblk = 1, nblocks
          do j = 1,ny_block
             do i = 1,nx_block
                algalN(i,j,1,iblk)    = aflds(i,j,4,iblk)    
                algalN(i,j,2,iblk)    = aflds(i,j,5,iblk)
                algalN(i,j,3,iblk)    = aflds(i,j,6,iblk)
                doc(i,j,1,iblk)       = aflds(i,j,7,iblk)
                doc(i,j,2,iblk)       = aflds(i,j,8,iblk)
                doc(i,j,3,iblk)       = aflds(i,j,9,iblk)
                dic(i,j,1,iblk)       = aflds(i,j,10,iblk)
                don(i,j,1,iblk)       = aflds(i,j,11,iblk)
                nit(i,j,iblk)         = aflds(i,j,12,iblk)
                sil(i,j,iblk)         = aflds(i,j,13,iblk)
                amm(i,j,iblk)         = aflds(i,j,14,iblk)
                dms(i,j,iblk)         = aflds(i,j,15,iblk)
                dmsp(i,j,iblk)        = aflds(i,j,16,iblk)
                hum(i,j,iblk)         = aflds(i,j,17,iblk)
                fep(i,j,1,iblk)       = aflds(i,j,18,iblk)
                fep(i,j,2,iblk)       = aflds(i,j,19,iblk)
                fed(i,j,1,iblk)       = aflds(i,j,20,iblk)
                fed(i,j,2,iblk)       = aflds(i,j,21,iblk)
                zaeros(i,j,1,iblk)    = aflds(i,j,22,iblk)
                zaeros(i,j,2,iblk)    = aflds(i,j,23,iblk)
                zaeros(i,j,3,iblk)    = aflds(i,j,24,iblk)
                zaeros(i,j,4,iblk)    = aflds(i,j,25,iblk)
                zaeros(i,j,5,iblk)    = aflds(i,j,26,iblk)
                zaeros(i,j,6,iblk)    = aflds(i,j,27,iblk)
             enddo    !i
          enddo      !j
       enddo        !iblk
       !$OMP END PARALLEL DO
    endif

    deallocate(aflds)

    !-----------------------------------------------------------------
    ! rotate zonal/meridional vectors to local coordinates
    ! compute data derived quantities
    ! unit conversions
    !-----------------------------------------------------------------

    ! Vector fields come in on T grid, but are oriented geographically
    ! need to rotate to pop-grid FIRST using ANGLET
    ! then interpolate to the U-cell centers  (otherwise we
    ! interpolate across the pole)
    ! use ANGLET which is on the T grid !

    call t_startf ('cice_imp_ocn')
    !$OMP PARALLEL DO PRIVATE(iblk,i,j,workx,worky)
    do iblk = 1, nblocks

       do j = 1,ny_block
          do i = 1,nx_block

             ! ocean
             workx      = uocn  (i,j,iblk) ! currents, m/s 
             worky      = vocn  (i,j,iblk)
             uocn(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid 
                  + worky*sin(ANGLET(i,j,iblk))
             vocn(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                  - workx*sin(ANGLET(i,j,iblk))

             workx      = ss_tltx  (i,j,iblk)           ! sea sfc tilt, m/m
             worky      = ss_tlty  (i,j,iblk)
             ss_tltx(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid 
                  + worky*sin(ANGLET(i,j,iblk))
             ss_tlty(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                  - workx*sin(ANGLET(i,j,iblk))

             sst(i,j,iblk) = sst(i,j,iblk) - Tffresh       ! sea sfc temp (C)

#ifdef RASM_MODS
             sss(i,j,iblk)=max(sss(i,j,iblk),c0)
#endif

             if (tfrz_option == 'minus1p8') then
                Tf (i,j,iblk) = -1.8_dbl_kind 
             elseif (tfrz_option == 'linear_salt') then
                Tf (i,j,iblk) = -0.0544_r8*sss(i,j,iblk)   ! THIS IS THE ORIGINAL POP FORMULA
             elseif (tfrz_option == 'mushy') then
                if (sss(i,j,iblk) > c0) then
                   Tf (i,j,iblk) = sss(i,j,iblk) / (-18.48_dbl_kind &
                                   + ((18.48_dbl_kind*p001)*sss(i,j,iblk)))
                else
                   Tf (i,j,iblk) = c0
                endif
             else
                write(nu_diag,*) subname,' ERROR: unknown tfrz_option = ',trim(tfrz_option)
                call shr_sys_abort(subname//' ERROR: unknown tfrz_option = '//trim(tfrz_option))
             endif

             if (cpl_bgc) then
                ! convert from mmol C/m^3 to mmol N/m^3
                algalN(i,j,1,iblk)    = algalN(i,j,1,iblk)/R_C2N(1)
                algalN(i,j,2,iblk)    = algalN(i,j,2,iblk)/R_C2N(2)
                algalN(i,j,3,iblk)    = algalN(i,j,3,iblk)/R_C2N(3)

                ! convert from mmol Fe/m^3 to umol Fe/m^3
                fep(i,j,1,iblk)       = fep(i,j,1,iblk) * 1000.0_dbl_kind
                fep(i,j,2,iblk)       = fep(i,j,2,iblk) * 1000.0_dbl_kind
                fed(i,j,1,iblk)       = fed(i,j,1,iblk) * 1000.0_dbl_kind
                fed(i,j,2,iblk)       = fed(i,j,2,iblk) * 1000.0_dbl_kind
             endif
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call t_stopf ('cice_imp_ocn')

    ! Interpolate ocean dynamics variables from T-cell centers to 
    ! U-cell centers.

    if (.not.prescribed_ice) then
       call t_startf ('cice_imp_t2u')
       call t2ugrid_vector(uocn)
       call t2ugrid_vector(vocn)
       call t2ugrid_vector(ss_tltx)
       call t2ugrid_vector(ss_tlty)
       call t_stopf ('cice_imp_t2u')
    end if

    ! Atmosphere variables are needed in T cell centers in
    ! subroutine stability and are interpolated to the U grid
    ! later as necessary.

    call t_startf ('cice_imp_atm')
    !$OMP PARALLEL DO PRIVATE(iblk,i,j,workx,worky)
    do iblk = 1, nblocks
       do j = 1, ny_block
          do i = 1, nx_block

             ! atmosphere
             workx      = uatm(i,j,iblk) ! wind velocity, m/s
             worky      = vatm(i,j,iblk) 
             uatm (i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid
                             + worky*sin(ANGLET(i,j,iblk))   ! note uatm, vatm, wind
             vatm (i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) & ! are on the T-grid here
                  - workx*sin(ANGLET(i,j,iblk))

             wind (i,j,iblk) = sqrt(uatm(i,j,iblk)**2 + vatm(i,j,iblk)**2)
             fsw  (i,j,iblk) = swvdr(i,j,iblk) + swvdf(i,j,iblk) &
                             + swidr(i,j,iblk) + swidf(i,j,iblk)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call t_stopf ('cice_imp_atm')

  end subroutine ice_import

  !===============================================================================

  subroutine ice_export( i2x )   

    !-----------------------------------------------------
    !
    ! Arguments
    real(r8), intent(inout) :: i2x(:,:)
    !
    ! Local Variables
    integer :: i, j, iblk, n, ij 
    integer :: ilo, ihi, jlo, jhi !beginning and end of physical domain
    integer (kind=int_kind)                                :: icells ! number of ocean/ice cells
    integer (kind=int_kind), dimension (nx_block*ny_block) :: indxi  ! compressed indices in i
    integer (kind=int_kind), dimension (nx_block*ny_block) :: indxj  ! compressed indices in i

    real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
         Tsrf  &      ! surface temperature
         ,  tauxa &      ! atmo/ice stress
         ,  tauya &
         ,  tauxo &      ! ice/ocean stress
         ,  tauyo &
         ,  ailohi       ! fractional ice area

    real (kind=dbl_kind) :: &
         workx, worky           ! tmps for converting grid

    real (kind=dbl_kind) :: &
         vonkar, zref, iceruf, tffresh 

    type(block)        :: this_block       ! block information for current block
    integer :: icnt,icnt1,iblk1,icnt1sum,icnt1max  ! gridcell and block counters
    logical :: flag
    character(len=*), parameter :: subname = '(ice_export)'
    !-----------------------------------------------------

    call icepack_query_parameters(Tffresh_out=Tffresh, vonkar_out=vonkar, zref_out=zref, &
       iceruf_out=iceruf)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
        file=__FILE__, line=__LINE__)

    flag=.false.

    !calculate ice thickness from aice and vice. Also
    !create Tsrf from the first tracer (trcr) in ice_state.F

    !$OMP PARALLEL DO PRIVATE(iblk,i,j,workx,worky)
    do iblk = 1, nblocks
       do j = 1, ny_block
          do i = 1, nx_block

             ! ice fraction
             ailohi(i,j,iblk) = min(aice(i,j,iblk), c1)

             ! surface temperature
             Tsrf(i,j,iblk)  = Tffresh + trcr(i,j,1,iblk)     !Kelvin (original ???)

             ! wind stress  (on POP T-grid:  convert to lat-lon)
             workx = strairxT(i,j,iblk)                             ! N/m^2
             worky = strairyT(i,j,iblk)                             ! N/m^2
             tauxa(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) &
                             - worky*sin(ANGLET(i,j,iblk))
             tauya(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                             + workx*sin(ANGLET(i,j,iblk))

             ! ice/ocean stress (on POP T-grid:  convert to lat-lon)
             workx = -strocnxT(i,j,iblk)                            ! N/m^2
             worky = -strocnyT(i,j,iblk)                            ! N/m^2
             tauxo(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) &
                             - worky*sin(ANGLET(i,j,iblk))
             tauyo(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                             + workx*sin(ANGLET(i,j,iblk))

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    do iblk = 1, nblocks
       do j = 1, ny_block
          do i = 1, nx_block
             if (tmask(i,j,iblk) .and. ailohi(i,j,iblk) < c0 ) then
                flag = .true.
             endif
          end do
       end do
    end do
    if (flag) then
       do iblk = 1, nblocks
          do j = 1, ny_block
             do i = 1, nx_block
                if (tmask(i,j,iblk) .and. ailohi(i,j,iblk) < c0 ) then
                   write(nu_diag,*) &
                        ' (ice) send: ERROR ailohi < 0.0 ',i,j,ailohi(i,j,iblk)
                   call shr_sys_flush(nu_diag)
                endif
             end do
          end do
       end do
    endif

    ! Fill export state i2x_i

    i2x(:,:) = spval_dbl

    n=0
    icnt1 = 0
    iblk1 = 0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       icnt = 0
       do j = jlo, jhi
          do i = ilo, ihi

             n = n+1

             !--- zero out fields with tmask for proper coupler accumulation in ice free areas

             if ( tmask(i,j,iblk)) i2x(:,n) = c0

             !-------states-------------------- 
             i2x(index_i2x_Si_ifrac ,n)    = ailohi(i,j,iblk)   

             if ( tmask(i,j,iblk) .and. ailohi(i,j,iblk) > c0 ) then
                icnt = icnt + 1
                !-------states-------------------- 
                i2x(index_i2x_Si_t     ,n)    = Tsrf(i,j,iblk)
                i2x(index_i2x_Si_avsdr ,n)    = alvdr(i,j,iblk)
                i2x(index_i2x_Si_anidr ,n)    = alidr(i,j,iblk)
                i2x(index_i2x_Si_avsdf ,n)    = alvdf(i,j,iblk)
                i2x(index_i2x_Si_anidf ,n)    = alidf(i,j,iblk)
                i2x(index_i2x_Si_u10  ,n)     = Uref(i,j,iblk)
                i2x(index_i2x_Si_tref  ,n)    = Tref(i,j,iblk)
                i2x(index_i2x_Si_qref  ,n)    = Qref(i,j,iblk)
                i2x(index_i2x_Si_snowh ,n)    = vsno(i,j,iblk) &
                     / ailohi(i,j,iblk)

                if (index_i2x_Si_logz0 > 0) then
                if (Cdn_atm(i,j,iblk) > c0) then
                 i2x(index_i2x_Si_logz0 ,n) = log(zref)-(vonkar/sqrt(Cdn_atm(i,j,iblk)))
                else
                 !--- tcraig, this should not happen but if it does, continue gracefully
                 write(nu_diag,*) trim(subname),&
                      ' WARNING: Cdn_atm error ',Cdn_atm(i,j,iblk),i,j,iblk
                 i2x(index_i2x_Si_logz0 ,n)    = log(iceruf)
                endif
                endif

                !--- a/i fluxes computed by ice
                i2x(index_i2x_Faii_taux ,n)   = tauxa(i,j,iblk)    
                i2x(index_i2x_Faii_tauy ,n)   = tauya(i,j,iblk)    
                i2x(index_i2x_Faii_lat  ,n)   = flat(i,j,iblk)     
                i2x(index_i2x_Faii_sen  ,n)   = fsens(i,j,iblk)    
                i2x(index_i2x_Faii_lwup ,n)   = flwout(i,j,iblk)   
                i2x(index_i2x_Faii_evap ,n)   = evap(i,j,iblk)     
                i2x(index_i2x_Faii_swnet,n)   = fswabs(i,j,iblk)

                !--- i/o fluxes computed by ice
                i2x(index_i2x_Fioi_melth,n)   = fhocn(i,j,iblk)
                i2x(index_i2x_Fioi_swpen,n)   = fswthru(i,j,iblk) ! hf from melting          
                i2x(index_i2x_Fioi_meltw,n)   = fresh(i,j,iblk)   ! h2o flux from melting    ???
                i2x(index_i2x_Fioi_salt ,n)   = fsalt(i,j,iblk)   ! salt flux from melting   ???
                i2x(index_i2x_Fioi_taux ,n)   = tauxo(i,j,iblk)   ! stress : i/o zonal       ???
                i2x(index_i2x_Fioi_tauy ,n)   = tauyo(i,j,iblk)   ! stress : i/o meridional  ???

                ! export biogeochemistry fields, if configured
                ! convert from mmol N/m^3 to mmol C/m^3
                if (index_i2x_Fioi_diat  > 0) i2x(index_i2x_Fioi_diat  ,n) = falgalN(i,j,1,iblk) * R_C2N(1)
                if (index_i2x_Fioi_sp    > 0) i2x(index_i2x_Fioi_sp    ,n) = falgalN(i,j,2,iblk) * R_C2N(2)
                if (index_i2x_Fioi_phaeo > 0) i2x(index_i2x_Fioi_phaeo ,n) = falgalN(i,j,3,iblk) * R_C2N(3)
                if (index_i2x_Fioi_doc   > 0) i2x(index_i2x_Fioi_doc   ,n) = fdoc(i,j,1,iblk) + fdoc(i,j,2,iblk)  
                if (index_i2x_Fioi_doc2  > 0) i2x(index_i2x_Fioi_doc2  ,n) = c0 !fdoc(i,j,2,iblk) 
                if (index_i2x_Fioi_doc3  > 0) i2x(index_i2x_Fioi_doc3  ,n) = c0 !fdoc(i,j,3,iblk)
                if (index_i2x_Fioi_dic   > 0) i2x(index_i2x_Fioi_dic   ,n) = c0 !fdic(i,j,1,iblk)
                if (index_i2x_Fioi_don   > 0) i2x(index_i2x_Fioi_don   ,n) = fdon(i,j,1,iblk) 
                if (index_i2x_Fioi_no3   > 0) i2x(index_i2x_Fioi_no3   ,n) = fnit(i,j,iblk)      
                if (index_i2x_Fioi_sio3  > 0) i2x(index_i2x_Fioi_sio3  ,n) = fsil(i,j,iblk)      
                if (index_i2x_Fioi_nh4   > 0) i2x(index_i2x_Fioi_nh4   ,n) = famm(i,j,iblk)       
                if (index_i2x_Fioi_dms   > 0) i2x(index_i2x_Fioi_dms   ,n) = fdms(i,j,iblk)  
                if (index_i2x_Fioi_dmspp > 0) i2x(index_i2x_Fioi_dmspp ,n) = c0                 
                if (index_i2x_Fioi_dmsp  > 0) i2x(index_i2x_Fioi_dmsp  ,n) = fdmsp(i,j,iblk)   
                if (index_i2x_Fioi_donr  > 0) i2x(index_i2x_Fioi_donr  ,n) = fhum(i,j,iblk)    
                ! convert from umol Fe/m^3 to mmol Fe/m^3
                if (index_i2x_Fioi_fep1  > 0) i2x(index_i2x_Fioi_fep1  ,n) = c0 !ffep(i,j,1,iblk) / 1000.0_dbl_kind
                if (index_i2x_Fioi_fep2  > 0) i2x(index_i2x_Fioi_fep2  ,n) = c0 !ffep(i,j,2,iblk) / 1000.0_dbl_kind
                if (index_i2x_Fioi_fed   > 0) i2x(index_i2x_Fioi_fed   ,n) = ffed(i,j,1,iblk) / 1000.0_dbl_kind
                if (index_i2x_Fioi_fed2  > 0) i2x(index_i2x_Fioi_fed2  ,n) = c0 !ffed(i,j,2,iblk) / 1000.0_dbl_kind
                if (index_i2x_Fioi_dust  > 0) i2x(index_i2x_Fioi_dust  ,n) = fdust(i,j,iblk)
            endif
          enddo    !i
       enddo    !j
!tcx       write(1000+my_task,'(a,5i8)') 'ice count0: ',my_task,istep1,nblocks,iblk,icnt
       if (icnt > 0) iblk1 = iblk1 + 1
       icnt1 = icnt1 + icnt
    enddo        !iblk

    if (mod(istep,diagfreq) == 0) then
!tcx    write(1000+my_task,'(a,5i8)') 'ice count1: ',my_task,istep1,nblocks,iblk1,icnt1
       call shr_mpi_max(icnt1,icnt1max,MPI_COMM_ICE,'icnt1max',all=.false.)
       call shr_mpi_sum(icnt1,icnt1sum,MPI_COMM_ICE,'icnt1sum',all=.false.)
       if (my_task == master_task) then
          write(nu_diag,'(a,f12.2,i8)') 'ice gridcell counts mean, max = ',float(icnt1sum)/float(get_num_procs()),icnt1max
       endif
    endif

  end subroutine ice_export

end module ice_import_export
