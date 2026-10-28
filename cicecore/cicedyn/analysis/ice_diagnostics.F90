!=======================================================================

! Diagnostic information output during run
!
! authors: Elizabeth C. Hunke, LANL
!          Bruce P. Briegleb, NCAR
!
! 2004: Block structure added by William Lipscomb
! 2006: Converted to free source form (F90) by Elizabeth Hunke

      module ice_diagnostics

      use ice_kinds_mod
      use ice_blocks, only: nx_block, ny_block
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0, c1
      use ice_calendar, only: istep1
      use ice_domain_size, only: nslyr
      use ice_fileunits, only: nu_diag
      use ice_fileunits, only: flush_fileunit
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_max_aero, icepack_max_iso
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_query_tracer_flags
      use icepack_intfc, only: icepack_query_tracer_indices

      implicit none
      private
      public :: runtime_diags, init_mass_diags, init_diags, debug_ice, &
                print_state, diagnostic_abort, &
                init_hdiag, hdiag_cell_values, hdiag_reduce, hdiag_derived, hdiag_global

      ! diagnostic output file
      character (len=char_len), public :: diag_file

      ! point print data

      logical (kind=log_kind), public :: &
         debug_model      , & ! if true, debug model at high level
         print_points     , & ! if true, print point data
         print_global         ! if true, print global data

      integer (kind=int_kind), public :: &
         debug_model_step = 0 ! begin printing at istep1=debug_model_step

      integer (kind=int_kind), parameter, public :: &
         npnt = 2             ! total number of points to be printed

      ! Set to true to identify unstable fast-moving ice.
      logical (kind=log_kind), parameter ::  &
         check_umax = .false. ! if true, check for speed > umax_stab

      real (kind=dbl_kind), parameter :: &
         umax_stab   = 1.0_dbl_kind     ! ice speed threshold for instability (m/s)

      real (kind=dbl_kind), parameter, public :: &
         aice_extmin = 0.15_dbl_kind    ! min aice value for ice extent calc

      real (kind=dbl_kind), parameter, public :: &
         maxval_spval = -0.9_dbl_kind*1.0e36_dbl_kind  ! spval to detect
                 ! undefined values returned from global_maxval.  if global_maxval
                 ! is applied to a region that does not exist (for instance
                 ! southern hemisphere in box cases), global_maxval
                 ! returns -1.0e36 which we want to avoid writing.  The
                 ! return value is checked against maxval_spval before writing.

      !-----------------------------------------------------------------
      ! Hemispheric diagnostics, for the log and the netcdf diagnostics
      ! (ice_diagnostics_netcdf)
      !
      ! Each diagnostic is a per grid cell value that is reduced to a value
      ! for each hemisphere, either an area weighted global sum or a global
      ! maximum (hdiag_reduce).  The same per grid cell values are used for the
      ! log, instantaneous netcdf output and, accumulated every timestep, for
      ! netcdf time averages.
      !
      ! The conservation diagnostics (and pond volume) are added when they are
      ! used, for the log (print_global) or netcdf, and the conservation budgets
      ! are derived in hdiag_global.
      !
      ! To add a diagnostic:
      !   1. add an index hd_* below, and for an optional diagnostic add it in
      !      init_hdiag with add_hdiag
      !   2. set its names, units and reduction in init_hdiag
      !   3. compute its per grid cell value in hdiag_cell_values
      !-----------------------------------------------------------------

      integer (kind=int_kind), parameter, public :: &
         hd_aice   = 1   , & ! sea-ice area
         hd_extent = 2   , & ! sea-ice extent
         hd_vice   = 3   , & ! sea-ice volume
         hd_vsno    = 4   , & ! snow volume, netcdf snow mass
         hd_ke     = 5   , & ! kinetic energy
         hd_albn   = 6   , & ! albedo numerator, for derived albedo
         hd_albd   = 7   , & ! albedo denominator, for derived albedo
         hd_vmax   = 8   , & ! max ice volume per unit area
         hd_umax   = 9   , & ! max ice speed
         hd_pmax   = 10  , & ! max ice strength
         n_hdiag_core = 10   ! number of hemispheric diagnostics always computed

      ! optional hemispheric diagnostics, index 0 if not used
      integer (kind=int_kind), public :: &
         hd_ptot   = 0   , & ! pond volume (tr_pond_topo)
         hd_etot   = 0   , & ! ice and snow energy (conservation)
         hd_stot   = 0   , & ! salt (conservation)
         hd_evp    = 0   , & ! evaporation (conservation)
         hd_brine  = 0   , & ! brine tracer volume (conservation and tr_brine)
         hd_sfsalt = 0   , & ! salt flux (conservation)
         hd_sfresh = 0   , & ! fresh water flux (conservation)
         hd_spond  = 0   , & ! pond water flux (conservation and tr_pond_topo)
         hd_fhocn  = 0   , & ! ocean heat flux (conservation)
         hd_fhatm  = 0   , & ! atmosphere heat flux (conservation)
         hd_fswnet = 0   , & ! absorbed shortwave (conservation)
         hd_fswdn  = 0   , & ! incoming shortwave (conservation)
         hd_fhfrz  = 0   , & ! frazil heat flux (conservation)
         hd_rn     = 0   , & ! rain (conservation)
         hd_sn     = 0   , & ! snow (conservation)
         hd_frz    = 0   , & ! frazil water flux (conservation)
         n_hdiag   = n_hdiag_core ! number of hemispheric diagnostics

      integer (kind=int_kind), dimension(icepack_max_iso), public :: &
         hd_fisoa  = 0   , & ! isotope flux from atmosphere (conservation and tr_iso)
         hd_fisoo  = 0   , & ! isotope flux to ocean
         hd_isotot = 0   , & ! isotope total
         hd_isomx  = 0       ! isotope maximum

      integer (kind=int_kind), dimension(icepack_max_aero), public :: &
         hd_faeroa = 0   , & ! aerosol flux from atmosphere (conservation and tr_aero)
         hd_faeroo = 0   , & ! aerosol flux to ocean
         hd_aerotot= 0   , & ! aerosol total
         hd_aeromx = 0       ! aerosol maximum

      type, public :: hdiag_field
         character (len=char_len) :: &
            vname_n     , & ! netcdf variable name, northern hemisphere
            vname_s     , & ! netcdf variable name, southern hemisphere
            lname       , & ! long name
            units       , & ! units after scale is applied
            sname           ! standard_name, blank for none
         character (len=4) :: &
            reduce          ! 'sum' area weighted global sum or 'max' global maximum
         real (kind=dbl_kind) :: &
            scale           ! unit conversion applied to the hemispheric value
         logical (kind=log_kind) :: &
            lwrite          ! written to file, false if only used for derived diagnostics
      end type hdiag_field

      type (hdiag_field), dimension(:), allocatable, public :: &
         hdiag               ! per grid cell diagnostics

      ! conservation diagnostics for one hemisphere, see hdiag_global.
      ! Water, salt and isotope/aerosol fluxes are amounts in one timestep.
      type, public :: hdiag_global_values
         real (kind=dbl_kind) :: &
            rn, sn, evp, frz, spond, sfresh, & ! water fluxes (kg in dt)
            mice, msnw, mpnd, mtot,          & ! ice, snow, pond and total mass (kg)
            delmi, delmx, flux, werr,        & ! mass change, ice mass change, water flux, error
            etot, delei,                     & ! energy, energy change (J)
            fhatm, fhocn, fhfrz, hnet, herr, & ! heat fluxes (W), net heat (J), error
            fswnet, fswdn, swerr,            & ! absorbed and incoming shortwave (W), error
            shmaxt,                          & ! brine tracer volume (m^3)
            stot, mslt, delmslt, sfsalt, serr  ! salt, salt mass, change, flux (kg in dt), error
         real (kind=dbl_kind), dimension(icepack_max_iso) :: &
            fisoa, fisoo, isotot, isor, isomx  ! isotope fluxes (kg in dt), total, error, max
         real (kind=dbl_kind), dimension(icepack_max_aero) :: &
            faeroa, faeroo, aerotot, aerr, aeromx ! aerosol fluxes (kg in dt), total, error, max
      end type hdiag_global_values

      logical (kind=log_kind), public :: &
         hdiag_initialized = .false., & ! hdiag has been initialized
         hdiag_global_on = .false.      ! conservation diagnostics are computed, for
                                        ! the log (print_global) or netcdf

      logical (kind=log_kind), dimension(2), public :: &
         hdiag_hemi          ! hemisphere has grid cells, otherwise values are missing


      real (kind=dbl_kind), dimension(npnt), public :: &
         latpnt           , & !  latitude of diagnostic points
         lonpnt               ! longitude of diagnostic points

      integer (kind=int_kind) :: &
         iindx            , & ! i index for points
         jindx            , & ! j index for points
         bindx                ! block index for points

      ! for water and heat budgets
      real (kind=dbl_kind), dimension(npnt) :: &
         pdhi             , & ! change in mean ice thickness (m)
         pdhs             , & ! change in mean snow thickness (m)
         pde                  ! change in ice and snow energy (W m-2)

      real (kind=dbl_kind), dimension(npnt), public :: &
         plat, plon           ! latitude, longitude of points

      integer (kind=int_kind), dimension(npnt), public :: &
         piloc, pjloc, pbloc, pmloc  ! location of diagnostic points

      integer (kind=int_kind), public :: &
         debug_model_i = -1,    &  ! location of debug_model point, local i index
         debug_model_j = -1,    &  ! location of debug_model point, local j index
         debug_model_iblk = -1, &  ! location of debug_model point, local block number
         debug_model_task = -1   ! location of debug_model point, local task number

      ! for hemispheric water and heat budgets
      real (kind=dbl_kind) :: &
         totmn            , & ! total ice/snow water mass (nh)
         totms            , & ! total ice/snow water mass (sh)
         totmin           , & ! total ice water mass (nh)
         totmis           , & ! total ice water mass (sh)
         totsn            , & ! total salt mass (nh)
         totss            , & ! total salt mass (sh)
         toten            , & ! total ice/snow energy (J)
         totes                ! total ice/snow energy (J)

      real (kind=dbl_kind), dimension(icepack_max_iso) :: &
         totison         , & ! total isotope mass
         totisos             ! total isotope mass

      real (kind=dbl_kind), dimension(icepack_max_aero) :: &
         totaeron         , & ! total aerosol mass
         totaeros             ! total aerosol mass



!=======================================================================

      contains

!=======================================================================

! Writes diagnostic info (max, min, global sums, etc) to standard out
!
! authors: Elizabeth C. Hunke, LANL
!          Bruce P. Briegleb, NCAR
!          Cecilia M. Bitz, UW

      subroutine runtime_diags (dt)

      use ice_arrays_column, only: floe_rad_c
      use ice_broadcast, only: broadcast_scalar
      use ice_constants, only: c1, c1000, c2, p001, p5, &
          field_loc_center, m2_to_km2
      use ice_domain, only: distrb_info, nblocks
      use ice_domain_size, only: ncat, n_iso, n_aero, max_blocks, nfsd
      use ice_flux, only: alvdr, alidr, alvdf, alidf, evap, fsnow, frazil, &
          fswabs, fswthru, flw, flwout, fsens, fsurf, flat, frzmlt_init, frain, fpond, &
          fhocn_ai, fsalt_ai, fresh_ai, frazil_diag, &
          update_ocn_f, cpl_frazil, Tair, Qa, fsw, fcondtop, meltt, meltb, meltl, snoice, &
          dsnow, congel, sst, sss, Tf, fhocn, &
          swvdr, swvdf, swidr, swidf, &
          alvdr_init, alvdf_init, alidr_init, alidf_init
      use ice_flux_bgc, only: faero_atm, faero_ocn, fiso_atm, fiso_ocn
      use ice_global_reductions, only: global_sum, global_sum_prod, global_maxval
      use ice_grid, only: lmask_n, lmask_s, tarean, tareas, grid_ice, grid_average_X2Y
      use ice_state   ! everything
! tcraig, this is likely to cause circular dependency because ice_prescribed_mod is high level routine
#ifdef CESMCOUPLED
      use ice_prescribed_mod, only: prescribed_ice
#endif

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         i, j, k, n, iblk, nc, &
         ktherm, &
         nt_tsfc, nt_aero, nt_fbri, nt_apnd, nt_hpnd, nt_fsd, &
         nt_isosno, nt_isoice, nt_rsnw, nt_rhos, nt_smice, nt_smliq

      logical (kind=log_kind) :: &
         tr_pond_topo, tr_brine, tr_iso, tr_aero, calc_Tsfc, tr_fsd, &
         tr_snow, snwgrain

      real (kind=dbl_kind) :: &
         rhow, rhos, rhoi, puny, awtvdr, awtidr, awtvdf, awtidf, &
         rhofresh, lfresh, lvap, ice_ref_salinity, Tffresh

      character (len=char_len) :: &
         snwredist, saltflux_option

      ! hemispheric state quantities
      real (kind=dbl_kind) :: &
         umaxn,   hmaxn,   shmaxn,    arean,   snwmxn, extentn, shmaxnt, &
         umaxs,   hmaxs,   shmaxs,    areas,   snwmxs, extents, shmaxst, &
         etotn,   mtotn,   micen,     msnwn,   pmaxn,  ketotn, &
         etots,   mtots,   mices,     msnws,   pmaxs,  ketots, &
         stotn, &
         stots, &
         urmsn,   albtotn,            mpndn,   ptotn,  spondn, &
         urmss,   albtots,            mpnds,   ptots,  sponds

      ! hemispheric flux quantities
      real (kind=dbl_kind) :: &
         rnn, snn, frzn,  hnetn, fhocnn, fhatmn,  fhfrzn, &
         rns, sns, frzs,  hnets, fhocns, fhatms,  fhfrzs, &
         fswnetn, fswnets, fswdnn, fswdns, swerrn, swerrs, &
         sfsaltn, sfreshn, evpn, fluxn , delmxn,  delmin, &
         sfsalts, sfreshs, evps, fluxs , delmxs,  delmis, &
         delein, werrn, herrn, msltn, delmsltn, serrn, &
         deleis, werrs, herrs, mslts, delmslts, serrs

      ! isotope diagnostics
      real (kind=dbl_kind), dimension(icepack_max_aero) :: &
         fisoan, fisoon, isorn, &
         fisoas, fisoos, isors, &
         isomx1n, isomx1s, &
         isototn, isotots

      ! aerosol diagnostics
      real (kind=dbl_kind), dimension(icepack_max_aero) :: &
         faeran, faeron, aerrn, &
         faeras, faeros, aerrs, &
         aeromx1n, aeromx1s, &
         aerototn, aerotots

      ! fields at diagnostic points
      real (kind=dbl_kind), dimension(npnt) :: &
         paice, pTair, pQa, pfsnow, pfrain, pfsw, pflw, &
         pTsfc, pevap, pfswabs, pflwout, pflat, pfsens, &
         pfsurf, pfcondtop, psst, psss, pTf, hiavg, hsavg, hbravg, &
         pfhocn, psalt, fsdavg, &
         pmeltt, pmeltb, pmeltl, psnoice, pdsnow, pfrazil, pcongel, &
         prsnwavg, prhosavg, psmicetot, psmliqtot, psmtot

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1, work2      ! temporary

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable :: &
         hfld              ! per grid cell hemispheric diagnostics

      real (kind=dbl_kind), dimension(n_hdiag) :: &
         hvaln, hvals      ! hemispheric diagnostics, northern and southern hemisphere

      type (hdiag_global_values) :: &
         gn, gs            ! print_global diagnostics, northern and southern hemisphere

      character(len=*), parameter :: subname = '(runtime_diags)'

      call icepack_query_parameters(ktherm_out=ktherm, calc_Tsfc_out=calc_Tsfc)
      call icepack_query_tracer_flags(tr_brine_out=tr_brine, tr_aero_out=tr_aero, &
           tr_pond_topo_out=tr_pond_topo, tr_fsd_out=tr_fsd, tr_iso_out=tr_iso, &
           tr_snow_out=tr_snow)
      call icepack_query_tracer_indices(nt_fbri_out=nt_fbri, nt_Tsfc_out=nt_Tsfc, &
           nt_aero_out=nt_aero, nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
           nt_fsd_out=nt_fsd,nt_isosno_out=nt_isosno, nt_isoice_out=nt_isoice, &
           nt_rsnw_out=nt_rsnw, nt_rhos_out=nt_rhos, &
           nt_smice_out=nt_smice, nt_smliq_out=nt_smliq)
      call icepack_query_parameters(Tffresh_out=Tffresh, rhos_out=rhos, &
           rhow_out=rhow, rhoi_out=rhoi, puny_out=puny, &
           awtvdr_out=awtvdr, awtidr_out=awtidr, awtvdf_out=awtvdf, awtidf_out=awtidf, &
           rhofresh_out=rhofresh, lfresh_out=lfresh, lvap_out=lvap, &
           ice_ref_salinity_out=ice_ref_salinity,snwredist_out=snwredist, &
           snwgrain_out=snwgrain, saltflux_option_out=saltflux_option)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! state of the ice
      !-----------------------------------------------------------------
      ! hemispheric quantities

      ! area, extent, volume, kinetic energy, albedo and maxima are the
      ! hemispheric diagnostics also written to netcdf (hdiag_reduce)
      if (.not. hdiag_initialized) call init_hdiag()
      allocate(hfld(nx_block,ny_block,max_blocks,n_hdiag))
      call hdiag_cell_values(hfld)
      call hdiag_reduce(hfld, c1, hvaln, hvals)
      deallocate(hfld)

      ! total ice area
      arean = hvaln(hd_aice)
      areas = hvals(hd_aice)
      arean = arean * m2_to_km2
      areas = areas * m2_to_km2

      ! ice extent (= area of grid cells with aice > aice_extmin)
      extentn = hvaln(hd_extent)
      extents = hvals(hd_extent)
      extentn = extentn * m2_to_km2
      extents = extents * m2_to_km2

      ! total ice volume
      shmaxn = hvaln(hd_vice)
      shmaxs = hvals(hd_vice)

      ! total snow volume
      snwmxn = hvaln(hd_vsno)
      snwmxs = hvals(hd_vsno)

      ! total pond volume
      ptotn = c0
      ptots = c0
      if (hd_ptot > 0) then
         ptotn = hvaln(hd_ptot)
         ptots = hvals(hd_ptot)
      endif

      ! total ice-snow kinetic energy
      ketotn = hvaln(hd_ke)
      ketots = hvals(hd_ke)

      ! rms ice speed, and average ice albedo over sunlit points
      call hdiag_derived(hvaln, urmsn, albtotn)
      call hdiag_derived(hvals, urmss, albtots)

      ! maximum ice volume (= mean thickness including open water)
      hmaxn = hvaln(hd_vmax)
      hmaxs = hvals(hd_vmax)

      ! maximum ice speed
      umaxn = hvaln(hd_umax)
      umaxs = hvals(hd_umax)

      ! Write warning message if ice speed is too big
      ! (Ice speeds of ~1 m/s or more usually indicate instability)

      if (check_umax) then
         call diag_ice_speed(work1)
         if (umaxn > umax_stab) then
            !$OMP PARALLEL DO PRIVATE(iblk,i,j)
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  if (abs(work1(i,j,iblk) - umaxn) < puny) then
                     write(nu_diag,*) ' '
                     write(nu_diag,*) 'Warning, large ice speed'
                     write(nu_diag,*) 'my_task, iblk, i, j, umaxn:', &
                                       my_task, iblk, i, j, umaxn
                  endif
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
         elseif (umaxs > umax_stab) then
            !$OMP PARALLEL DO PRIVATE(iblk,i,j)
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  if (abs(work1(i,j,iblk) - umaxs) < puny) then
                     write(nu_diag,*) ' '
                     write(nu_diag,*) 'Warning, large ice speed'
                     write(nu_diag,*) 'my_task, iblk, i, j, umaxs:', &
                                       my_task, iblk, i, j, umaxs
                  endif
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
         endif   ! umax
      endif      ! check_umax

      ! maximum ice strength

      pmaxn = hvaln(hd_pmax)
      pmaxs = hvals(hd_pmax)

      pmaxn = pmaxn / c1000   ! convert to kN/m
      pmaxs = pmaxs / c1000

      if (print_global) then
         ! conservation diagnostics, the same as the netcdf diagnostics
         call hdiag_global(hvaln, 1, dt, gn)
         call hdiag_global(hvals, 2, dt, gs)

         etotn = gn%etot    ; etots = gs%etot
         stotn = gn%stot    ; stots = gs%stot
         evpn  = gn%evp     ; evps  = gs%evp
         shmaxnt = gn%shmaxt; shmaxst = gs%shmaxt
         sfsaltn = gn%sfsalt; sfsalts = gs%sfsalt
         sfreshn = gn%sfresh; sfreshs = gs%sfresh
         spondn = gn%spond  ; sponds = gs%spond
         fhocnn = gn%fhocn  ; fhocns = gs%fhocn
         fhatmn = gn%fhatm  ; fhatms = gs%fhatm
         fswnetn = gn%fswnet; fswnets = gs%fswnet
         fswdnn = gn%fswdn  ; fswdns = gs%fswdn
         fhfrzn = gn%fhfrz  ; fhfrzs = gs%fhfrz
         rnn = gn%rn        ; rns = gs%rn
         snn = gn%sn        ; sns = gs%sn
         frzn = gn%frz      ; frzs = gs%frz
         micen = gn%mice    ; mices = gs%mice
         msnwn = gn%msnw    ; msnws = gs%msnw
         mpndn = gn%mpnd    ; mpnds = gs%mpnd
         mtotn = gn%mtot    ; mtots = gs%mtot
         delmin = gn%delmi  ; delmis = gs%delmi
         delmxn = gn%delmx  ; delmxs = gs%delmx
         fluxn = gn%flux    ; fluxs = gs%flux
         werrn = gn%werr    ; werrs = gs%werr
         delein = gn%delei  ; deleis = gs%delei
         hnetn = gn%hnet    ; hnets = gs%hnet
         herrn = gn%herr    ; herrs = gs%herr
         swerrn = gn%swerr  ; swerrs = gs%swerr
         msltn = gn%mslt    ; mslts = gs%mslt
         delmsltn = gn%delmslt; delmslts = gs%delmslt
         serrn = gn%serr    ; serrs = gs%serr

         do n = 1, n_iso
            fisoan(n)  = gn%fisoa(n) ; fisoas(n)  = gs%fisoa(n)
            fisoon(n)  = gn%fisoo(n) ; fisoos(n)  = gs%fisoo(n)
            isototn(n) = gn%isotot(n); isotots(n) = gs%isotot(n)
            isomx1n(n) = gn%isomx(n) ; isomx1s(n) = gs%isomx(n)
            isorn(n)   = gn%isor(n)  ; isors(n)   = gs%isor(n)
         enddo

         do n = 1, n_aero
            faeran(n)   = gn%faeroa(n) ; faeras(n)   = gs%faeroa(n)
            faeron(n)   = gn%faeroo(n) ; faeros(n)   = gs%faeroo(n)
            aerototn(n) = gn%aerotot(n); aerotots(n) = gs%aerotot(n)
            aeromx1n(n) = gn%aeromx(n) ; aeromx1s(n) = gs%aeromx(n)
            aerrn(n)    = gn%aerr(n)   ; aerrs(n)    = gs%aerr(n)
         enddo
      endif                     ! print_global

      if (print_points) then

      !-----------------------------------------------------------------
      ! state of the ice and associated fluxes for 2 defined points
      ! NOTE these are computed for the last timestep only (not avg)
      !-----------------------------------------------------------------

         call total_energy (work1)
         call total_salt   (work2)

         do n = 1, npnt
            if (my_task == pmloc(n)) then
               i = piloc(n)
               j = pjloc(n)
               iblk = pbloc(n)

               pTair(n) = Tair(i,j,iblk) - Tffresh ! air temperature
               pQa(n) = Qa(i,j,iblk)               ! specific humidity
               pfsnow(n) = fsnow(i,j,iblk)*dt/rhos ! snowfall
               pfrain(n) = frain(i,j,iblk)*dt/rhow ! rainfall
               pfsw(n) = fsw(i,j,iblk)             ! shortwave radiation
               pflw(n) = flw(i,j,iblk)             ! longwave radiation
               paice(n) = aice(i,j,iblk)           ! ice area

               fsdavg(n) = c0                      ! avg floe effective radius
               hiavg(n) = c0                       ! avg snow/ice thickness
               hsavg(n) = c0
               hbravg(n) = c0                      ! avg brine thickness
               if (paice(n) /= c0) then
                  hiavg(n) = vice(i,j,iblk)/paice(n)
                  hsavg(n) = vsno(i,j,iblk)/paice(n)
                  if (tr_brine) hbravg(n) = trcr(i,j,nt_fbri,iblk)* hiavg(n)
                  if (tr_fsd) then
! not sure why this does not work
!                     do k = 1, nfsd
!                        fsdavg(n) = fsdavg(n) &
!                                  + trcr(i,j,nt_fsd+k-1,iblk) * floe_rad_c(k) &
!                                  / trcr(i,j,nt_fsd+k-1,iblk)
!                     enddo
! this works
                     do nc = 1, ncat
                     do k = 1, nfsd
                        fsdavg(n) = fsdavg(n) &
                                  + trcrn(i,j,nt_fsd+k-1,nc,iblk) * floe_rad_c(k) &
                                  * aicen(i,j,nc,iblk) / paice(n)
                     enddo
                     enddo
                  endif
               endif
               if (tr_snow) then      ! snow tracer quantities
                  prsnwavg (n) = c0   ! avg snow grain radius
                  prhosavg (n) = c0   ! avg snow density
                  psmicetot(n) = c0   ! total mass of ice in snow (kg/m2)
                  psmliqtot(n) = c0   ! total mass of liquid in snow (kg/m2)
                  psmtot   (n) = c0   ! total mass of snow volume (kg/m2)
                  if (vsno(i,j,iblk) > c0) then
                     do k = 1, nslyr
                        prsnwavg (n) = prsnwavg (n) + trcr(i,j,nt_rsnw +k-1,iblk) ! snow grain radius
                        prhosavg (n) = prhosavg (n) + trcr(i,j,nt_rhos +k-1,iblk) ! compacted snow density
                        psmicetot(n) = psmicetot(n) + trcr(i,j,nt_smice+k-1,iblk) * vsno(i,j,iblk)
                        psmliqtot(n) = psmliqtot(n) + trcr(i,j,nt_smliq+k-1,iblk) * vsno(i,j,iblk)
                     end do
                  endif
                  psmtot   (n) = rhos * vsno(i,j,iblk) ! mass of ice in standard density snow
                  prsnwavg (n) = prsnwavg (n) / real(nslyr,kind=dbl_kind) ! snow grain radius
                  prhosavg (n) = prhosavg (n) / real(nslyr,kind=dbl_kind) ! compacted snow density
                  psmicetot(n) = psmicetot(n) / real(nslyr,kind=dbl_kind) ! mass of ice in snow
                  psmliqtot(n) = psmliqtot(n) / real(nslyr,kind=dbl_kind) ! mass of liquid in snow
               end if
               psalt(n) = c0
               if (vice(i,j,iblk) /= c0) psalt(n) = work2(i,j,iblk)/vice(i,j,iblk)
               pTsfc(n) = trcr(i,j,nt_Tsfc,iblk)   ! ice/snow sfc temperature
               pevap(n) = evap(i,j,iblk)*dt/rhoi   ! sublimation/condensation
               pfswabs(n) = fswabs(i,j,iblk)       ! absorbed solar flux
               pflwout(n) = flwout(i,j,iblk)       ! outward longwave flux
               pflat(n) = flat(i,j,iblk)           ! latent heat flux
               pfsens(n) = fsens(i,j,iblk)         ! sensible heat flux
               pfsurf(n) = fsurf(i,j,iblk)         ! total sfc heat flux
               pfcondtop(n) = fcondtop(i,j,iblk)   ! top sfc cond flux
               pmeltt(n) = meltt(i,j,iblk)         ! top melt
               pmeltb(n) = meltb(i,j,iblk)         ! bottom melt
               pmeltl(n) = meltl(i,j,iblk)         ! lateral melt
               psnoice(n) = snoice(i,j,iblk)       ! snow ice
               pdsnow(n) = dsnow(i,j,iblk)         ! snow change
               pfrazil(n) = frazil(i,j,iblk)       ! frazil ice
               pcongel(n) = congel(i,j,iblk)       ! congelation ice
               pdhi(n) = vice(i,j,iblk) - pdhi(n)  ! ice thickness change
               pdhs(n) = vsno(i,j,iblk) - pdhs(n)  ! snow thickness change
               pde(n) =-(work1(i,j,iblk)- pde(n))/dt ! ice/snow energy change
               psst(n) = sst(i,j,iblk)             ! sea surface temperature
               psss(n) = sss(i,j,iblk)             ! sea surface salinity
               pTf(n) = Tf(i,j,iblk)               ! freezing temperature
               pfhocn(n) = -fhocn(i,j,iblk)        ! ocean heat used by ice

            endif  ! my_task = pmloc

            call broadcast_scalar(pTair    (n), pmloc(n))
            call broadcast_scalar(pQa      (n), pmloc(n))
            call broadcast_scalar(pfsnow   (n), pmloc(n))
            call broadcast_scalar(pfrain   (n), pmloc(n))
            call broadcast_scalar(pfsw     (n), pmloc(n))
            call broadcast_scalar(pflw     (n), pmloc(n))
            call broadcast_scalar(paice    (n), pmloc(n))
            call broadcast_scalar(hsavg    (n), pmloc(n))
            call broadcast_scalar(hiavg    (n), pmloc(n))
            call broadcast_scalar(fsdavg   (n), pmloc(n))
            call broadcast_scalar(psalt    (n), pmloc(n))
            call broadcast_scalar(hbravg   (n), pmloc(n))
            call broadcast_scalar(pTsfc    (n), pmloc(n))
            call broadcast_scalar(pevap    (n), pmloc(n))
            call broadcast_scalar(pfswabs  (n), pmloc(n))
            call broadcast_scalar(pflwout  (n), pmloc(n))
            call broadcast_scalar(pflat    (n), pmloc(n))
            call broadcast_scalar(pfsens   (n), pmloc(n))
            call broadcast_scalar(pfsurf   (n), pmloc(n))
            call broadcast_scalar(pfcondtop(n), pmloc(n))
            call broadcast_scalar(pmeltt   (n), pmloc(n))
            call broadcast_scalar(pmeltb   (n), pmloc(n))
            call broadcast_scalar(pmeltl   (n), pmloc(n))
            call broadcast_scalar(psnoice  (n), pmloc(n))
            call broadcast_scalar(pdsnow   (n), pmloc(n))
            call broadcast_scalar(psmtot   (n), pmloc(n))
            call broadcast_scalar(prsnwavg (n), pmloc(n))
            call broadcast_scalar(prhosavg (n), pmloc(n))
            call broadcast_scalar(psmicetot(n), pmloc(n))
            call broadcast_scalar(psmliqtot(n), pmloc(n))
            call broadcast_scalar(pfrazil  (n), pmloc(n))
            call broadcast_scalar(pcongel  (n), pmloc(n))
            call broadcast_scalar(pdhi     (n), pmloc(n))
            call broadcast_scalar(pdhs     (n), pmloc(n))
            call broadcast_scalar(pde      (n), pmloc(n))
            call broadcast_scalar(psst     (n), pmloc(n))
            call broadcast_scalar(psss     (n), pmloc(n))
            call broadcast_scalar(pTf      (n), pmloc(n))
            call broadcast_scalar(pfhocn   (n), pmloc(n))

         enddo                  ! npnt
      endif                     ! print_points

      !-----------------------------------------------------------------
      ! start spewing
      !-----------------------------------------------------------------

      if (my_task == master_task) then

        write(nu_diag,899) 'Arctic','Antarctic'

        write(nu_diag,901) 'total ice area  (km^2) = ',arean,  areas
        write(nu_diag,901) 'total ice extent(km^2) = ',extentn,extents
        write(nu_diag,901) 'total ice volume (m^3) = ',shmaxn, shmaxs
        write(nu_diag,901) 'total snw volume (m^3) = ',snwmxn, snwmxs
        write(nu_diag,901) 'tot kinetic energy (J) = ',ketotn, ketots
        write(nu_diag,900) 'rms ice speed    (m/s) = ',urmsn,  urmss
        write(nu_diag,900) 'average albedo         = ',albtotn,albtots
        write(nu_diag,900) 'max ice volume     (m) = ',hmaxn,  hmaxs
        write(nu_diag,900) 'max ice speed    (m/s) = ',umaxn,  umaxs
        write(nu_diag,900) 'max strength    (kN/m) = ',pmaxn,  pmaxs

        if (print_global) then  ! global diags for conservations checks

! tcraig, this is likely to cause circular dependency because ice_prescribed_mod is high level routine
#ifdef CESMCOUPLED
         if (prescribed_ice) then
          write (nu_diag,*) '----------------------------'
          write (nu_diag,*)   'This is the prescribed ice option.'
          write (nu_diag,*)   'Heat and water will not be conserved.'
          write (nu_diag,*) '----------------------------'
         endif
#endif

         write(nu_diag,*) '----------------------------'
         write(nu_diag,901) 'arwt rain h2o kg in dt = ',rnn,rns
         write(nu_diag,901) 'arwt snow h2o kg in dt = ',snn,sns
         write(nu_diag,901) 'arwt evap h2o kg in dt = ',evpn,evps
         write(nu_diag,901) 'arwt frzl h2o kg in dt = ',frzn,frzs
         if (tr_pond_topo) &
         write(nu_diag,901) 'arwt fpnd h2o kg in dt = ',spondn,sponds
         write(nu_diag,901) 'arwt frsh h2o kg in dt = ',sfreshn,sfreshs

         write(nu_diag,901) 'arwt ice mass (kg)     = ',micen,mices
         write(nu_diag,901) 'arwt snw mass (kg)     = ',msnwn,msnws
         if (tr_pond_topo) &
         write(nu_diag,901) 'arwt pnd mass (kg)     = ',mpndn,mpnds

         write(nu_diag,901) 'arwt tot mass (kg)     = ',mtotn,mtots
         write(nu_diag,901) 'arwt tot mass chng(kg) = ',delmin,delmis
         write(nu_diag,901) 'arwt water flux        = ',fluxn,fluxs
         if (update_ocn_f) then
           write (nu_diag,*) '(=rain+snow+evap-fresh)  '
         else
           write (nu_diag,*) '(=rain+snow+evap+frzl-fresh)  '
         endif
         write(nu_diag,901) 'water flux error       = ',werrn,werrs

         write(nu_diag,*) '----------------------------'
         write(nu_diag,901) 'arwt atm heat flux (W) = ',fhatmn,fhatms
         write(nu_diag,901) 'arwt ocn heat flux (W) = ',fhocnn,fhocns
         write(nu_diag,901) 'arwt frzl heat flux(W) = ',fhfrzn,fhfrzs
         write(nu_diag,901) 'arwt tot energy    (J) = ',etotn,etots
         write(nu_diag,901) 'arwt net heat      (J) = ',hnetn,hnets
         write(nu_diag,901) 'arwt tot energy chng(J)= ',delein,deleis
         write(nu_diag,901) 'arwt heat error        = ',herrn,herrs
         write(nu_diag,*) '----------------------------'
         write(nu_diag,901) 'arwt incoming sw (W)   = ',fswdnn,fswdns
         write(nu_diag,901) 'arwt absorbed sw (W)   = ',fswnetn,fswnets
         write(nu_diag,901) 'arwt swdn error        = ',swerrn,swerrs

         write(nu_diag,*) '----------------------------'
         write(nu_diag,901) 'total brine tr (m^3)   = ',shmaxnt, shmaxst
         write(nu_diag,901) 'arwt salt mass (kg)    = ',msltn,mslts
         write(nu_diag,901) 'arwt salt mass chng(kg)= ',delmsltn, &
                                                        delmslts
         write(nu_diag,901) 'arwt salt flx in dt(kg)= ',sfsaltn, &
                                                        sfsalts
         write(nu_diag,901) 'arwt salt flx error    = ',serrn,serrs

         write(nu_diag,*) '----------------------------'
         if (tr_iso) then
         do n = 1, n_iso
         write(nu_diag,*)   '  isotope ',n
         write(nu_diag,901) 'fiso_atm (kg/m2)       = ', fisoan(n), fisoas(n)
         write(nu_diag,901) 'fiso_ocn (kg/m2)       = ', fisoon(n), fisoos(n)
         write(nu_diag,901) 'total iso (kg/m2)      = ', isototn(n), isotots(n)
         write(nu_diag,901) 'iso error              = ', isorn(n), isors(n)
         write(nu_diag,901) 'maximum iso (kg/m2)    = ', isomx1n(n),isomx1s(n)
         enddo
         write(nu_diag,*) '----------------------------'
         endif ! tr_iso
         if (tr_aero) then
         do n = 1, n_aero
         write(nu_diag,*)   '  aerosol ',n
         write(nu_diag,901) 'faero_atm (kg/m2)      = ', faeran(n), faeras(n)
         write(nu_diag,901) 'faero_ocn (kg/m2)      = ', faeron(n), faeros(n)
         write(nu_diag,901) 'total aero (kg/m2)     = ', aerototn(n), aerotots(n)
         write(nu_diag,901) 'aero error             = ', aerrn(n), aerrs(n)
         write(nu_diag,901) 'maximum aero (kg/m2)   = ', aeromx1n(n),aeromx1s(n)
         enddo
         write(nu_diag,*) '----------------------------'
         endif ! tr_aero

        endif                    ! print_global

       call flush_fileunit(nu_diag)

      !-----------------------------------------------------------------
      ! diagnostics for Arctic and Antarctic points
      !-----------------------------------------------------------------

       if (print_points) then

        write(nu_diag,*) '                         '
        write(nu_diag,902) '       Lat, Long         ',plat(1),plon(1), &
                                                       plat(2),plon(2)
        write(nu_diag,903) '  my_task, iblk, i, j     ', &
                              pmloc(1),pbloc(1),piloc(1),pjloc(1), &
                              pmloc(2),pbloc(2),piloc(2),pjloc(2)
        write(nu_diag,*) '----------atm----------'
        write(nu_diag,900) 'air temperature (C)    = ',pTair(1),pTair(2)
        write(nu_diag,900) 'specific humidity      = ',pQa(1),pQa(2)
        write(nu_diag,900) 'snowfall (m)           = ',pfsnow(1), &
                                                       pfsnow(2)
        write(nu_diag,900) 'rainfall (m)           = ',pfrain(1), &
                                                       pfrain(2)
        if (.not.calc_Tsfc) then
           write(nu_diag,900) 'total surface heat flux= ',pfsurf(1),pfsurf(2)
           write(nu_diag,900) 'top sfc conductive flux= ',pfcondtop(1), &
                                                          pfcondtop(2)
           write(nu_diag,900) 'latent heat flx        = ',pflat(1),pflat(2)
        else
           write(nu_diag,900) 'shortwave radiation sum= ',pfsw(1),pfsw(2)
           write(nu_diag,900) 'longwave radiation     = ',pflw(1),pflw(2)
        endif
        write(nu_diag,*) '----------ice----------'
        write(nu_diag,900) 'area fraction          = ',paice(1),paice(2)
        write(nu_diag,900) 'avg ice thickness (m)  = ',hiavg(1),hiavg(2)
        write(nu_diag,900) 'avg snow depth (m)     = ',hsavg(1),hsavg(2)
        write(nu_diag,900) 'avg salinity (ppt)     = ',psalt(1),psalt(2)
        write(nu_diag,900) 'avg brine thickness (m)= ',hbravg(1),hbravg(2)
        if (tr_fsd) &
        write(nu_diag,900) 'avg fsd rep radius (m) = ',fsdavg(1),fsdavg(2)

        if (calc_Tsfc) then
           write(nu_diag,900) 'surface temperature(C) = ',pTsfc(1),pTsfc(2)
           write(nu_diag,900) 'absorbed shortwave flx = ',pfswabs(1), &
                                                          pfswabs(2)
           write(nu_diag,900) 'outward longwave flx   = ',pflwout(1), &
                                                          pflwout(2)
           write(nu_diag,900) 'sensible heat flx      = ',pfsens(1), &
                                                          pfsens(2)
           write(nu_diag,900) 'latent heat flx        = ',pflat(1),pflat(2)
        endif
        write(nu_diag,900) 'subl/cond (m ice)      = ',pevap(1),pevap(2)
        write(nu_diag,900) 'top melt (m)           = ',pmeltt(1) &
                                                      ,pmeltt(2)
        write(nu_diag,900) 'bottom melt (m)        = ',pmeltb(1) &
                                                      ,pmeltb(2)
        write(nu_diag,900) 'lateral melt (m)       = ',pmeltl(1) &
                                                      ,pmeltl(2)
        write(nu_diag,900) 'new ice (m)            = ',pfrazil(1), &
                                                       pfrazil(2)
        write(nu_diag,900) 'congelation (m)        = ',pcongel(1), &
                                                       pcongel(2)
        write(nu_diag,900) 'snow-ice (m)           = ',psnoice(1), &
                                                       psnoice(2)
        write(nu_diag,900) 'snow change (m)        = ',pdsnow(1), &
                                                       pdsnow(2)
        write(nu_diag,900) 'effective dhi (m)      = ',pdhi(1),pdhi(2)
        write(nu_diag,900) 'effective dhs (m)      = ',pdhs(1),pdhs(2)
        write(nu_diag,900) 'intnl enrgy chng(W/m^2)= ',pde (1),pde (2)

        if (tr_snow) then
           if (trim(snwredist) /= 'none') then
              write(nu_diag,900) 'avg snow density(kg/m3)= ',prhosavg(1) &
                                                            ,prhosavg(2)
           endif
           if (snwgrain) then
              write(nu_diag,900) 'avg snow grain radius  = ',prsnwavg(1) &
                                                            ,prsnwavg(2)
              write(nu_diag,900) 'mass ice in snow(kg/m2)= ',psmicetot(1) &
                                                            ,psmicetot(2)
              write(nu_diag,900) 'mass liq in snow(kg/m2)= ',psmliqtot(1) &
                                                            ,psmliqtot(2)
              write(nu_diag,900) 'mass std snow   (kg/m2)= ',psmtot(1) &
                                                            ,psmtot(2)
              write(nu_diag,900) 'max  ice+liq    (kg/m2)= ',rhow * hsavg(1) &
                                                            ,rhow * hsavg(2)
           endif
        endif

        write(nu_diag,*) '----------ocn----------'
        write(nu_diag,900) 'sst (C)                = ',psst(1),psst(2)
        write(nu_diag,900) 'sss (ppt)              = ',psss(1),psss(2)
        write(nu_diag,900) 'freezing temp (C)      = ',pTf(1),pTf(2)
        write(nu_diag,900) 'heat used (W/m^2)      = ',pfhocn(1), &
                                                       pfhocn(2)

       endif                    ! print_points
      endif                     ! my_task = master_task

  899 format (27x,a24,2x,a24)
  900 format (a25,2x,f24.17,2x,f24.17)
  901 format (a25,2x,1pe24.17,2x,1pe24.17)
  902 format (a25,10x,f6.1,1x,f6.1,9x,f6.1,1x,f6.1)
  903 format (a25,5x,i4,1x,i4,1x,i4,1x,i4,7x,i4,1x,i4,1x,i4,1x,i4)

      end subroutine runtime_diags

!=======================================================================

! Computes global combined ice and snow mass sum
!
! author: Elizabeth C. Hunke, LANL

      subroutine init_mass_diags

      use ice_constants, only: field_loc_center
      use ice_domain, only: distrb_info, nblocks
      use ice_domain_size, only: n_iso, n_aero, ncat, max_blocks
      use ice_global_reductions, only: global_sum
      use ice_grid, only: tareas, tarean
      use ice_state, only: aicen, vice, vsno, trcrn, trcr

      integer (kind=int_kind) :: n, i, j, k, iblk, &
         nt_hpnd, nt_apnd, nt_aero, nt_isosno, nt_isoice

      logical (kind=log_kind) :: &
         tr_iso, tr_aero, tr_pond_topo

      real (kind=dbl_kind) :: &
         shmaxn, snwmxn,  shmaxs, snwmxs, totpn, totps, &
         rhoi, rhos, rhofresh

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1, work2

      character(len=*), parameter :: subname = '(init_mass_diags)'

      call icepack_query_tracer_flags(tr_aero_out=tr_aero, tr_pond_topo_out=tr_pond_topo)
      call icepack_query_tracer_flags(tr_iso_out=tr_iso)
      call icepack_query_tracer_indices(nt_isosno_out=nt_isosno, nt_isoice_out=nt_isoice, &
         nt_hpnd_out=nt_hpnd, nt_apnd_out=nt_apnd, nt_aero_out=nt_aero)
      call icepack_query_parameters( &
         rhoi_out=rhoi, rhos_out=rhos, rhofresh_out=rhofresh)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      ! total ice volume
      shmaxn = global_sum(vice, distrb_info, field_loc_center, tarean)
      shmaxs = global_sum(vice, distrb_info, field_loc_center, tareas)

      ! total snow volume
      snwmxn = global_sum(vsno, distrb_info, field_loc_center, tarean)
      snwmxs = global_sum(vsno, distrb_info, field_loc_center, tareas)

      ! north/south ice mass
      totmin = rhoi*shmaxn
      totmis = rhoi*shmaxs

      ! north/south ice+snow mass
      totmn = totmin + rhos*snwmxn
      totms = totmis + rhos*snwmxs

      ! north/south ice+snow energy
      call total_energy (work1)
      toten = global_sum(work1, distrb_info, field_loc_center, tarean)
      totes = global_sum(work1, distrb_info, field_loc_center, tareas)

      ! north/south salt
      call total_salt (work2)
      totsn = global_sum(work2, distrb_info, field_loc_center, tarean)
      totss = global_sum(work2, distrb_info, field_loc_center, tareas)


      if (print_points) then
         do n = 1, npnt
            if (my_task == pmloc(n)) then
               i = piloc(n)
               j = pjloc(n)
               iblk = pbloc(n)

               pdhi(n) = vice(i,j,iblk)
               pdhs(n) = vsno(i,j,iblk)
               pde(n) = work1(i,j,iblk)
            endif
         enddo  ! npnt
      endif                     ! print_points

      if (tr_iso) then
         do n=1,n_iso
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,k)
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  work1(i,j,iblk) = c0
                  do k = 1, n_iso
                     work1(i,j,iblk) = work1(i,j,iblk) &
                                     + vsno(i,j,iblk)*trcr(i,j,nt_isosno+k-1,iblk) &
                                     + vice(i,j,iblk)*trcr(i,j,nt_isoice+k-1,iblk)
                  enddo
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
            totison(n)= global_sum(work1, distrb_info, field_loc_center, tarean)
            totisos(n)= global_sum(work1, distrb_info, field_loc_center, tareas)
         enddo
      endif

      if (tr_aero) then
         do n=1,n_aero
            !$OMP PARALLEL DO PRIVATE(iblk,i,j)
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  work1(i,j,iblk) = trcr(i,j,nt_aero  +4*(n-1),iblk)*vsno(i,j,iblk) &
                                  + trcr(i,j,nt_aero+1+4*(n-1),iblk)*vsno(i,j,iblk) &
                                  + trcr(i,j,nt_aero+2+4*(n-1),iblk)*vice(i,j,iblk) &
                                  + trcr(i,j,nt_aero+3+4*(n-1),iblk)*vice(i,j,iblk)
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
            totaeron(n)= global_sum(work1, distrb_info, field_loc_center, tarean)
            totaeros(n)= global_sum(work1, distrb_info, field_loc_center, tareas)
         enddo
      endif

      if (tr_pond_topo) then
         totpn = c0
         totps = c0
         !$OMP PARALLEL DO PRIVATE(iblk,i,j,n)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               work1(i,j,iblk) = c0
               do n = 1, ncat
                  work1(i,j,iblk) = work1(i,j,iblk)  &
                                  + aicen(i,j,n,iblk) &
                                  * trcrn(i,j,nt_apnd,n,iblk) &
                                  * trcrn(i,j,nt_hpnd,n,iblk)
               enddo
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
         totpn = global_sum(work1, distrb_info, field_loc_center, tarean)
         totps = global_sum(work1, distrb_info, field_loc_center, tareas)

         ! north/south ice+snow+pond mass
         totmn = totmn + totpn*rhofresh
         totms = totms + totps*rhofresh
      endif

      end subroutine init_mass_diags

!=======================================================================

! Computes total energy of ice and snow in a grid cell.
!
! authors: E. C. Hunke, LANL

      subroutine total_energy (work)

      use ice_domain, only: nblocks
      use ice_domain_size, only: ncat, nilyr, nslyr, max_blocks
      use ice_grid, only: tmask
      use ice_state, only: vicen, vsnon, trcrn

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(out) :: &
         work      ! total energy

      ! local variables

      integer (kind=int_kind) :: &
        icells                ! number of ocean/ice cells

      integer (kind=int_kind), dimension (nx_block*ny_block) :: &
        indxi, &              ! compressed indices in i/j directions
        indxj

      integer (kind=int_kind) :: &
        i, j, k, n, iblk, ij, &
        nt_qice, nt_qsno

      character(len=*), parameter :: subname = '(total_energy)'

      call icepack_query_tracer_indices(nt_qice_out=nt_qice, nt_qsno_out=nt_qsno)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,n,k,ij,icells,indxi,indxj)
      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

         icells = 0
         do j = 1, ny_block
         do i = 1, nx_block
            if (tmask(i,j,iblk)) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif                  ! tmask
         enddo
         enddo

         work(:,:,iblk) = c0

      !-----------------------------------------------------------------
      ! Aggregate
      !-----------------------------------------------------------------

         do n = 1, ncat
            do k = 1, nilyr
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  work(i,j,iblk) = work(i,j,iblk) &
                                 + trcrn(i,j,nt_qice+k-1,n,iblk) &
                                 * vicen(i,j,n,iblk) / real(nilyr,kind=dbl_kind)
               enddo            ! ij
            enddo               ! k

            do k = 1, nslyr
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  work(i,j,iblk) = work(i,j,iblk) &
                                 + trcrn(i,j,nt_qsno+k-1,n,iblk) &
                                 * vsnon(i,j,n,iblk) / real(nslyr,kind=dbl_kind)
               enddo            ! ij
            enddo               ! k
         enddo                  ! n

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      end subroutine total_energy

!=======================================================================

! Computes bulk salinity of ice and snow in a grid cell.
! author: E. C. Hunke, LANL

      subroutine total_salt (work)

      use ice_domain, only: nblocks
      use ice_domain_size, only: ncat, nilyr, max_blocks
      use ice_grid, only: tmask
      use ice_state, only: vicen, trcrn

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(out) :: &
         work      ! total salt

      ! local variables

      integer (kind=int_kind) :: &
        icells                ! number of ocean/ice cells

      integer (kind=int_kind), dimension (nx_block*ny_block) :: &
        indxi, &              ! compressed indices in i/j directions
        indxj

      integer (kind=int_kind) :: &
        i, j, k, n, iblk, ij, &
        nt_sice

      character(len=*), parameter :: subname = '(total_salt)'

      call icepack_query_tracer_indices(nt_sice_out=nt_sice)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,n,k,ij,icells,indxi,indxj)
      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

         icells = 0
         do j = 1, ny_block
         do i = 1, nx_block
            if (tmask(i,j,iblk)) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif                  ! tmask
         enddo
         enddo

         work(:,:,iblk) = c0

      !-----------------------------------------------------------------
      ! Aggregate
      !-----------------------------------------------------------------

         do n = 1, ncat
            do k = 1, nilyr
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  work(i,j,iblk) = work(i,j,iblk) &
                                 + trcrn(i,j,nt_sice+k-1,n,iblk) &
                                 * vicen(i,j,n,iblk) / real(nilyr,kind=dbl_kind)
               enddo            ! ij
            enddo               ! k
         enddo                  ! n

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      end subroutine total_salt

!=======================================================================

!  Find tasks for diagnostic points.
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL

      subroutine init_diags

      use ice_grid, only: hm, TLAT, TLON
      use ice_blocks, only: block, get_block
      use ice_constants, only: c180, c360, p5
      use ice_domain, only: blocks_ice, distrb_info, nblocks
      use ice_global_reductions, only: global_minval, global_maxval

      real (kind=dbl_kind) :: &
         rad_to_deg, &
         puny, &
         latdis  , & ! latitude distance
         londis  , & ! longitude distance
         totdis  , & ! total distance
         mindis  , & ! minimum distance from desired location
         mindis_g    ! global minimum distance from desired location

      integer (kind=int_kind) :: &
         n           , & ! index for point search
         i,j         , & ! grid indices
         iblk        , & ! block index
         ilo,ihi,jlo,jhi ! beginning and end of physical domain

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(init_diags)'

!tcraig, do this all the time now for print_points_state usage
!      if (print_points) then

         call icepack_query_parameters(puny_out=puny, rad_to_deg_out=rad_to_deg)
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
            file=__FILE__, line=__LINE__)

         if (my_task==master_task) then
            write(nu_diag,*) ' '
            write(nu_diag,*) ' Find indices of diagnostic points '
         endif

         piloc(:) = -1
         pjloc(:) = -1
         pbloc(:) = -1
         pmloc(:) = -999
         plat(:)  = -999._dbl_kind
         plon(:)  = -999._dbl_kind

         ! find minimum distance to diagnostic points on this processor
         do n = 1, npnt
            if (lonpnt(n) > c180) lonpnt(n) = lonpnt(n) - c360

            iindx = 0
            jindx = 0
            bindx = 0
            mindis = 540.0_dbl_kind !  360. + 180.

            if (abs(latpnt(n)) < c360 .and. abs(lonpnt(n)) < c360) then

            ! MDT, 09/2017: Comment out OpenMP directives since loop is not thread-safe
            ! This is computing closest point, Could add a CRITICAL but it's just initialization
            !!$XXXOMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,latdis,londis,totdis)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi

               do j = jlo, jhi
               do i = ilo, ihi
                  if (hm(i,j,iblk) > p5) then
                     latdis = abs(latpnt(n)-TLAT(i,j,iblk)*rad_to_deg)
                     londis = abs(lonpnt(n)-TLON(i,j,iblk)*rad_to_deg) &
                            * cos(TLAT(i,j,iblk))
                     totdis = sqrt(latdis**2 + londis**2)
                     if (totdis < mindis) then
                        mindis = totdis
                        jindx = j
                        iindx = i
                        bindx = iblk
                     endif      ! totdis < mindis
                  endif         ! hm > p5
               enddo            ! i
               enddo            ! j
            enddo               ! iblk
            !!$XXXOMP END PARALLEL DO

            endif

            ! find global minimum distance to diagnostic points
            mindis_g = global_minval(mindis, distrb_info)

            ! save indices of minimum-distance grid cell
            if (mindis <= 180.0 .and. abs(mindis_g - mindis) < puny) then
               piloc(n) = iindx
               pjloc(n) = jindx
               pbloc(n) = bindx
               pmloc(n) = my_task
               plat(n) = TLAT(iindx,jindx,bindx)*rad_to_deg
               plon(n) = TLON(iindx,jindx,bindx)*rad_to_deg
            endif

            ! communicate to all processors
            piloc(n) = global_maxval(piloc(n), distrb_info)
            pjloc(n) = global_maxval(pjloc(n), distrb_info)
            pbloc(n) = global_maxval(pbloc(n), distrb_info)
            pmloc(n) = global_maxval(pmloc(n), distrb_info)
            plat(n)  = global_maxval(plat(n), distrb_info)
            plon(n)  = global_maxval(plon(n), distrb_info)

            ! write to log file
            if (my_task==master_task) then
               write(nu_diag,*) ' '
               write(nu_diag,100) n,latpnt(n),lonpnt(n),plat(n),plon(n), &
                    piloc(n), pjloc(n), pbloc(n), pmloc(n)
            endif
 100        format(' found point',i4/ &
               '   lat    lon   TLAT   TLON     i     j   block  task'/ &
                4(f6.1,1x),1x,4(i4,2x) )

         enddo                  ! npnt
!      endif                     ! print_points

      ! hemispheric diagnostics for the log and netcdf
      call init_hdiag()

      end subroutine init_diags

!=======================================================================

! Kinetic energy of ice and snow per unit area on T points (J/m^2)

      subroutine diag_kinetic_energy(work)

      use ice_constants, only: p5
      use ice_domain, only: nblocks
      use ice_domain_size, only: max_blocks
      use ice_grid, only: grid_ice, grid_average_X2Y
      use ice_state, only: uvel, vvel, uvelE, vvelE, uvelN, vvelN, vice, vsno

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(out) :: &
         work              ! kinetic energy per unit area

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk

      real (kind=dbl_kind) :: &
         rhos, rhoi

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         uvelT, vvelT      ! u,v on T points

      character(len=*), parameter :: subname = '(diag_kinetic_energy)'

      call icepack_query_parameters(rhos_out=rhos, rhoi_out=rhoi)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (grid_ice == 'B') then
         call grid_average_X2Y('A',uvel ,'U',uvelT,'T')
         call grid_average_X2Y('A',vvel ,'U',vvelT,'T')
      elseif (grid_ice == 'C') then
         call grid_average_X2Y('A',uvelE,'E',uvelT,'T')
         call grid_average_X2Y('A',vvelN,'N',vvelT,'T')
      elseif (grid_ice == 'CD') then
         call grid_average_X2Y('A',uvelE,'E',uvelN,'N',uvelT,'T')
         call grid_average_X2Y('A',vvelE,'E',vvelN,'N',vvelT,'T')
      endif

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            work(i,j,iblk) = p5 * (rhos*vsno(i,j,iblk) + rhoi*vice(i,j,iblk)) &
                                * (uvelT(i,j,iblk)**2 + vvelT(i,j,iblk)**2)
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine diag_kinetic_energy

!=======================================================================

! Ice speed (m/s), the maximum over E and N points for C and CD grids

      subroutine diag_ice_speed(work)

      use ice_domain, only: nblocks
      use ice_domain_size, only: max_blocks
      use ice_grid, only: grid_ice, grid_average_X2Y
      use ice_state, only: uvel, vvel, uvelE, vvelE, uvelN, vvelN

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(out) :: &
         work              ! ice speed

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work2             ! temporary

      character(len=*), parameter :: subname = '(diag_ice_speed)'

      if (grid_ice == 'CD') then
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               work(i,j,iblk) = max(sqrt(uvelE(i,j,iblk)**2 + vvelE(i,j,iblk)**2), &
                                    sqrt(uvelN(i,j,iblk)**2 + vvelN(i,j,iblk)**2))
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
      elseif (grid_ice == 'C') then
         ! map uvelE to N and vvelN to E then compute max on E and N
         call grid_average_X2Y('A',uvelE,'E',work ,'N')  ! work  =~ uvelN
         call grid_average_X2Y('A',vvelN,'N',work2,'E')  ! work2 =~ vvelE
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               work(i,j,iblk) = max(sqrt(uvelE(i,j,iblk)**2 + work2(i,j,iblk)**2), &
                                    sqrt(work(i,j,iblk)**2 + vvelN(i,j,iblk)**2))
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
      else
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               work(i,j,iblk) = sqrt(uvel(i,j,iblk)**2 + vvel(i,j,iblk)**2)
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
      endif

      end subroutine diag_ice_speed

!=======================================================================

! This routine is useful for debugging
! author Elizabeth C. Hunke, LANL

      subroutine debug_ice(iblk, plabeld)

      character (len=*), intent(in) :: plabeld
      integer (kind=int_kind), intent(in) :: iblk

      ! local
      character(len=*), parameter :: subname='(debug_ice)'

      if (istep1 >= debug_model_step) then

         ! set debug point to 1st global point if not set as local values
         if (debug_model_i < 0 .and. debug_model_j < 0 .and. &
             debug_model_iblk < 0 .and. debug_model_task < 0) then
            debug_model_i    = piloc(1)
            debug_model_j    = pjloc(1)
            debug_model_task = pmloc(1)
            debug_model_iblk = pbloc(1)
         endif

         ! if debug point is messed up, abort
         if (debug_model_i < 0 .or. debug_model_j < 0 .or. &
             debug_model_iblk < 0 .or. debug_model_task < 0) then
            call abort_ice (subname//'ERROR: debug_model_[i,j,iblk,mytask] not set correctly')
         endif

         ! write out debug info
         if (debug_model_iblk == iblk .and. debug_model_task == my_task) then
            call print_state(plabeld,debug_model_i,debug_model_j,debug_model_iblk)
         endif

      endif

      end subroutine debug_ice

!=======================================================================

! This routine is useful for debugging.
! author: Elizabeth C. Hunke, LANL

      subroutine print_state(plabel,i,j,iblk)

      use ice_grid, only: grid_ice
      use ice_blocks, only: block, get_block
      use ice_domain, only: blocks_ice
      use ice_domain_size, only: ncat, nilyr, nslyr, nfsd
      use ice_grid, only: TLAT, TLON
      use ice_state, only: aice, aice0, aicen, vicen, vsnon, uvel, vvel, &
                           uvelE, vvelE, uvelN, vvelN, trcrn
      use ice_flux, only: uatm, vatm, potT, Tair, Qa, flw, frain, fsnow, &
          fsens, flat, evap, flwout, swvdr, swvdf, swidr, swidf, rhoa, &
          frzmlt, sst, sss, Tf, Tref, Qref, Uref, uocn, vocn, strtltxU, strtltyU

      character (len=*), intent(in) :: plabel

      integer (kind=int_kind), intent(in) :: &
          i, j       , & ! horizontal indices
          iblk           ! block index

      ! local variables

      real (kind=dbl_kind) :: &
           eidebug, esdebug, &
           qi, qs, Tsnow, si, &
           rad_to_deg, puny, rhoi, lfresh, rhos, cp_ice

      integer (kind=int_kind) :: n, k, nt_Tsfc, nt_qice, nt_qsno, nt_fsd, &
           nt_isosno, nt_isoice, nt_sice, nt_smice, nt_smliq

      logical (kind=log_kind) :: tr_fsd, tr_iso, tr_snow

      type (block) :: &
           this_block    ! block information for current block

      character(len=*), parameter :: subname = '(print_state)'

      call icepack_query_tracer_flags(tr_fsd_out=tr_fsd, tr_iso_out=tr_iso, &
           tr_snow_out=tr_snow)
      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice, &
           nt_qsno_out=nt_qsno, nt_sice_out=nt_sice, nt_fsd_out=nt_fsd, &
           nt_isosno_out=nt_isosno, nt_isoice_out=nt_isoice, &
           nt_smice_out=nt_smice, nt_smliq_out=nt_smliq)
      call icepack_query_parameters( &
           rad_to_deg_out=rad_to_deg, puny_out=puny, rhoi_out=rhoi, lfresh_out=lfresh, &
           rhos_out=rhos, cp_ice_out=cp_ice)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      this_block = get_block(blocks_ice(iblk),iblk)

      write(nu_diag,*) subname,' ',trim(plabel)
      write(nu_diag,*) subname,' istep1, my_task, i, j, iblk:', &
                        istep1, my_task, i, j, iblk
      write(nu_diag,*) subname,' Global block:', this_block%block_id
      write(nu_diag,*) subname,' Global i and j:', &
                        this_block%i_glob(i), &
                        this_block%j_glob(j)
      write (nu_diag,*) subname,' Lat, Lon (degrees):', &
                        TLAT(i,j,iblk)*rad_to_deg, &
                        TLON(i,j,iblk)*rad_to_deg
      write(nu_diag,*) ' '
      write(nu_diag,*) 'aice ', aice(i,j,iblk)
      write(nu_diag,*) 'aice0', aice0(i,j,iblk)
      do n = 1, ncat
         write(nu_diag,*) ' '
         write(nu_diag,*) 'n =',n
         write(nu_diag,*) 'aicen', aicen(i,j,n,iblk)
         write(nu_diag,*) 'vicen', vicen(i,j,n,iblk)
         write(nu_diag,*) 'vsnon', vsnon(i,j,n,iblk)
         if (aicen(i,j,n,iblk) > puny) then
            write(nu_diag,*) 'hin', vicen(i,j,n,iblk)/aicen(i,j,n,iblk)
            write(nu_diag,*) 'hsn', vsnon(i,j,n,iblk)/aicen(i,j,n,iblk)
         endif
         write(nu_diag,*) 'Tsfcn',trcrn(i,j,nt_Tsfc,n,iblk)
         if (tr_fsd) write(nu_diag,*) 'afsdn',trcrn(i,j,nt_fsd,n,iblk) ! fsd cat 1
! layer 1 diagnostics
!         if (tr_iso)  write(nu_diag,*) 'isosno',trcrn(i,j,nt_isosno,n,iblk) ! isotopes in snow
!         if (tr_iso)  write(nu_diag,*) 'isoice',trcrn(i,j,nt_isoice,n,iblk) ! isotopes in ice
!         if (tr_snow) write(nu_diag,*) 'smice', trcrn(i,j,nt_smice, n,iblk) ! ice mass in snow
!         if (tr_snow) write(nu_diag,*) 'smliq', trcrn(i,j,nt_smliq, n,iblk) ! liquid mass in snow
         write(nu_diag,*) ' '

! dynamics (transport and/or ridging) causes the floe size distribution to become non-normal
!         if (tr_fsd) then
!         if (abs(sum(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n,iblk))-c1) > puny) &
!            write(nu_diag,*) 'afsdn not normal', &
!                 sum(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n,iblk)), &
!                     trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n,iblk)
!         endif

      enddo                     ! n

      eidebug = c0
      do n = 1,ncat
         do k = 1,nilyr
            qi = trcrn(i,j,nt_qice+k-1,n,iblk)
            write(nu_diag,*) 'qice, cat ',n,' layer ',k, qi
            eidebug = eidebug + qi
            if (aicen(i,j,n,iblk) > puny) then
               write(nu_diag,*)  'qi/rhoi', qi/rhoi
            endif
         enddo
         write(nu_diag,*) ' '
      enddo
      write(nu_diag,*) 'qice(i,j)',eidebug
      write(nu_diag,*) ' '

      esdebug = c0
      do n = 1,ncat
         if (vsnon(i,j,n,iblk) > puny) then
            do k = 1,nslyr
               qs = trcrn(i,j,nt_qsno+k-1,n,iblk)
               write(nu_diag,*) 'qsnow, cat ',n,' layer ',k, qs
               esdebug = esdebug + qs
               Tsnow = (Lfresh + qs/rhos) / cp_ice
               write(nu_diag,*) 'qs/rhos', qs/rhos
               write(nu_diag,*) 'Tsnow', Tsnow
            enddo
            write(nu_diag,*) ' '
         endif
      enddo
      write(nu_diag,*) 'qsnow(i,j)',esdebug
      write(nu_diag,*) ' '

      do n = 1,ncat
         do k = 1,nilyr
            si = trcrn(i,j,nt_sice+k-1,n,iblk)
            write(nu_diag,*) 'sice, cat ',n,' layer ',k, si
         enddo
      enddo
      write(nu_diag,*) ' '

      write(nu_diag,*) 'uvel(i,j)',uvel(i,j,iblk)
      write(nu_diag,*) 'vvel(i,j)',vvel(i,j,iblk)
      if (grid_ice == 'C') then
         write(nu_diag,*) 'uvelE(i,j)',uvelE(i,j,iblk)
         write(nu_diag,*) 'uvelN(i,j)',uvelN(i,j,iblk)
      elseif (grid_ice == 'CD') then
         write(nu_diag,*) 'uvelE(i,j)',uvelE(i,j,iblk)
         write(nu_diag,*) 'vvelE(i,j)',vvelE(i,j,iblk)
         write(nu_diag,*) 'uvelN(i,j)',uvelN(i,j,iblk)
         write(nu_diag,*) 'vvelN(i,j)',vvelN(i,j,iblk)
      endif

      write(nu_diag,*) ' '
      write(nu_diag,*) 'atm states and fluxes'
      write(nu_diag,*) '            uatm    = ',uatm (i,j,iblk)
      write(nu_diag,*) '            vatm    = ',vatm (i,j,iblk)
      write(nu_diag,*) '            potT    = ',potT (i,j,iblk)
      write(nu_diag,*) '            Tair    = ',Tair (i,j,iblk)
      write(nu_diag,*) '            Qa      = ',Qa   (i,j,iblk)
      write(nu_diag,*) '            rhoa    = ',rhoa (i,j,iblk)
      write(nu_diag,*) '            swvdr   = ',swvdr(i,j,iblk)
      write(nu_diag,*) '            swvdf   = ',swvdf(i,j,iblk)
      write(nu_diag,*) '            swidr   = ',swidr(i,j,iblk)
      write(nu_diag,*) '            swidf   = ',swidf(i,j,iblk)
      write(nu_diag,*) '            flw     = ',flw  (i,j,iblk)
      write(nu_diag,*) '            frain   = ',frain(i,j,iblk)
      write(nu_diag,*) '            fsnow   = ',fsnow(i,j,iblk)
      write(nu_diag,*) ' '
      write(nu_diag,*) 'ocn states and fluxes'
      write(nu_diag,*) '            frzmlt  = ',frzmlt  (i,j,iblk)
      write(nu_diag,*) '            sst     = ',sst     (i,j,iblk)
      write(nu_diag,*) '            sss     = ',sss     (i,j,iblk)
      write(nu_diag,*) '            Tf      = ',Tf      (i,j,iblk)
      write(nu_diag,*) '            uocn    = ',uocn    (i,j,iblk)
      write(nu_diag,*) '            vocn    = ',vocn    (i,j,iblk)
      write(nu_diag,*) '            strtltxU= ',strtltxU(i,j,iblk)
      write(nu_diag,*) '            strtltyU= ',strtltyU(i,j,iblk)
      write(nu_diag,*) ' '
      write(nu_diag,*) 'srf states and fluxes'
      write(nu_diag,*) '            Tref    = ',Tref  (i,j,iblk)
      write(nu_diag,*) '            Qref    = ',Qref  (i,j,iblk)
      write(nu_diag,*) '            Uref    = ',Uref  (i,j,iblk)
      write(nu_diag,*) '            fsens   = ',fsens (i,j,iblk)
      write(nu_diag,*) '            flat    = ',flat  (i,j,iblk)
      write(nu_diag,*) '            evap    = ',evap  (i,j,iblk)
      write(nu_diag,*) '            flwout  = ',flwout(i,j,iblk)
      write(nu_diag,*) ' '
      call flush_fileunit(nu_diag)

      end subroutine print_state

!=======================================================================

! prints error information prior to aborting

      subroutine diagnostic_abort(istop, jstop, iblk, stop_label)

      use ice_blocks, only: block, get_block
      use ice_domain, only: blocks_ice

      integer (kind=int_kind), intent(in) :: &
         istop, jstop, & ! indices of grid cell where model aborts
         iblk            ! block index

      character (len=*), intent(in) :: stop_label

      ! local variables

      real (kind=dbl_kind) :: rad_to_deg

      type (block) :: &
         this_block      ! block information for current block

      character(len=*), parameter :: subname = '(diagnostic_abort)'

      call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      this_block = get_block(blocks_ice(iblk),iblk)

      call flush_fileunit(nu_diag)
      if (istop > 0 .and. jstop > 0) then
         call print_state(trim(stop_label),istop,jstop,iblk)
      else
         write (nu_diag,*) subname,' istep1, my_task, iblk =', &
                            istep1, my_task, iblk
         write (nu_diag,*) subname,' Global block:', this_block%block_id
      endif
      call flush_fileunit(nu_diag)
      call abort_ice (subname//'ERROR: '//trim(stop_label))

      end subroutine diagnostic_abort

!=======================================================================

! Define the per grid cell hemispheric netcdf diagnostics, see the
! comments with hd_* in the module header

      subroutine init_hdiag

      use ice_constants, only: field_loc_center
      use ice_domain, only: distrb_info
      use ice_domain_size, only: max_blocks, n_iso, n_aero
      use ice_global_reductions, only: global_sum
      use ice_grid, only: tarean, tareas
      use ice_calendar, only: diag_netcdf_freq

      integer (kind=int_kind) :: &
         n                 ! tracer index

      logical (kind=log_kind) :: &
         tr_pond_topo, tr_brine, tr_iso, tr_aero

      character (len=char_len) :: &
         cn                ! tracer index string

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1             ! temporary

      real (kind=dbl_kind) :: &
         rhos              ! snow density

      character(len=*), parameter :: subname = '(init_hdiag)'

      if (hdiag_initialized) return

      call icepack_query_parameters(rhos_out=rhos)
      call icepack_query_tracer_flags(tr_pond_topo_out=tr_pond_topo, tr_brine_out=tr_brine, &
           tr_iso_out=tr_iso, tr_aero_out=tr_aero)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      allocate(hdiag(n_hdiag_core + 16 + 4*n_iso + 4*n_aero))
      n_hdiag = n_hdiag_core

      ! conservation diagnostics for the log (print_global) and always for netcdf
      hdiag_global_on = print_global .or. diag_netcdf_freq /= 'x'

      ! a hemisphere without grid cells (e.g. regional domains) is written as missing
      work1(:,:,:) = c1
      hdiag_hemi(1) = global_sum(work1, distrb_info, field_loc_center, tarean) > c0
      hdiag_hemi(2) = global_sum(work1, distrb_info, field_loc_center, tareas) > c0

      ! CMIP sea-ice diagnostics, in CMIP units
      call def_hdiag(hd_aice, 'siarean'  , 'siareas'  , 'sea-ice area'  , '1e6 km2', 'sum', &
                      1.0e-12_dbl_kind, sname='sea_ice_area')
      call def_hdiag(hd_extent, 'siextentn', 'siextents', 'sea-ice extent', '1e6 km2', 'sum', &
                      1.0e-12_dbl_kind, sname='sea_ice_extent')
      call def_hdiag(hd_vice, 'sivoln'   , 'sivols'   , 'sea-ice volume', '1e3 km3', 'sum', &
                      1.0e-12_dbl_kind, sname='sea_ice_volume')
      call def_hdiag(hd_vsno, 'sisnmassn', 'sisnmasss', 'snow mass on sea ice', 'kg', 'sum', &
                      rhos)

      ! other diagnostics, the same quantities as in the log
      call def_hdiag(hd_ke  , 'kinetic_energy_nh'  , 'kinetic_energy_sh'  , &
                      'total kinetic energy', 'J', 'sum', c1)
      call def_hdiag(hd_albn, 'albedo_num_nh'      , 'albedo_num_sh'      , &
                      'albedo numerator', '1', 'sum', c1, lwrite=.false.)
      call def_hdiag(hd_albd, 'albedo_den_nh'      , 'albedo_den_sh'      , &
                      'albedo denominator', '1', 'sum', c1, lwrite=.false.)
      call def_hdiag(hd_vmax, 'ice_volume_max_nh'  , 'ice_volume_max_sh'  , &
                      'max ice volume per unit area', 'm', 'max', c1)
      call def_hdiag(hd_umax, 'ice_speed_max_nh'   , 'ice_speed_max_sh'   , &
                      'max ice speed', 'm/s', 'max', c1)
      call def_hdiag(hd_pmax, 'ice_strength_max_nh', 'ice_strength_max_sh', &
                      'max ice strength', 'kN/m', 'max', 1.0e-3_dbl_kind)

      ! optional diagnostics, used for the conservation budgets (hdiag_global)
      if (tr_pond_topo) call add_hdiag(hd_ptot, 'pond_volume', 'total pond volume', 'm^3', 'sum')
      if (hdiag_global_on) then
         call add_hdiag(hd_etot  , 'total_energy'    , 'total ice and snow energy', 'J'  , 'sum')
         call add_hdiag(hd_stot  , 'total_salt'      , 'total salt'               , 'ppt m^3', 'sum')
         call add_hdiag(hd_evp   , 'evap_flux'       , 'evaporation'              , 'kg s-1', 'sum')
         if (tr_brine) &
         call add_hdiag(hd_brine , 'brine_volume'    , 'total brine tracer volume', 'm^3', 'sum')
         call add_hdiag(hd_sfsalt, 'salt_flux'       , 'salt flux'                , 'kg s-1', 'sum')
         call add_hdiag(hd_sfresh, 'fresh_flux'      , 'fresh water flux'         , 'kg s-1', 'sum')
         if (tr_pond_topo) &
         call add_hdiag(hd_spond , 'pond_flux'       , 'pond water flux'          , 'kg s-1', 'sum')
         call add_hdiag(hd_fhocn , 'ocn_heat_flux'   , 'ocean heat flux'          , 'W'  , 'sum')
         call add_hdiag(hd_fhatm , 'atm_heat_flux'   , 'atmosphere heat flux'     , 'W'  , 'sum')
         call add_hdiag(hd_fswnet, 'sw_absorbed_flux', 'absorbed shortwave'       , 'W'  , 'sum')
         call add_hdiag(hd_fswdn , 'sw_down_flux'    , 'incoming shortwave'       , 'W'  , 'sum')
         call add_hdiag(hd_fhfrz , 'frazil_heat_flux', 'frazil heat flux'         , 'W'  , 'sum')
         call add_hdiag(hd_rn    , 'rain_flux'       , 'rain'                     , 'kg s-1', 'sum')
         call add_hdiag(hd_sn    , 'snow_flux'       , 'snow'                     , 'kg s-1', 'sum')
         call add_hdiag(hd_frz   , 'frazil_flux'     , 'frazil ice formation'     , 'kg s-1', 'sum')
         if (tr_iso) then
            do n = 1, n_iso
               write(cn,'(i0)') n
               call add_hdiag(hd_fisoa(n) , 'isotope'//trim(cn)//'_atm_flux', 'isotope flux from atmosphere', 'kg s-1', 'sum')
               call add_hdiag(hd_fisoo(n) , 'isotope'//trim(cn)//'_ocn_flux', 'isotope flux to ocean'       , 'kg s-1', 'sum')
               call add_hdiag(hd_isotot(n), 'isotope'//trim(cn)//'_total'   , 'total isotope'               , 'kg'    , 'sum')
               call add_hdiag(hd_isomx(n) , 'isotope'//trim(cn)//'_max'     , 'maximum isotope'             , 'kg m-2', 'max')
            enddo
         endif
         if (tr_aero) then
            do n = 1, n_aero
               write(cn,'(i0)') n
               call add_hdiag(hd_faeroa(n) , 'aerosol'//trim(cn)//'_atm_flux', 'aerosol flux from atmosphere', 'kg s-1', 'sum')
               call add_hdiag(hd_faeroo(n) , 'aerosol'//trim(cn)//'_ocn_flux', 'aerosol flux to ocean'       , 'kg s-1', 'sum')
               call add_hdiag(hd_aerotot(n), 'aerosol'//trim(cn)//'_total'   , 'total aerosol'               , 'kg'    , 'sum')
               call add_hdiag(hd_aeromx(n) , 'aerosol'//trim(cn)//'_max'     , 'maximum aerosol'             , 'kg m-2', 'max')
            enddo
         endif
      endif

      hdiag_initialized = .true.

      end subroutine init_hdiag

!=======================================================================

! Add an optional per grid cell hemispheric diagnostic, used for the
! conservation budgets.  The netcdf output is from hdiag_global, not directly.

      subroutine add_hdiag(k, vname, lname, units, reduce)

      integer (kind=int_kind), intent(out) :: &
         k               ! index in hdiag

      character (len=*), intent(in) :: &
         vname       , & ! variable name, _nh and _sh are appended
         lname       , & ! long name
         units       , & ! units
         reduce          ! 'sum' or 'max'

      character(len=*), parameter :: subname = '(add_hdiag)'

      if (n_hdiag >= size(hdiag)) then
         call abort_ice(subname//' ERROR: too many hemispheric diagnostics', &
            file=__FILE__, line=__LINE__)
      endif
      n_hdiag = n_hdiag + 1
      k = n_hdiag
      call def_hdiag(k, trim(vname)//'_nh', trim(vname)//'_sh', lname, units, reduce, c1, &
                     lwrite=.false.)

      end subroutine add_hdiag

!=======================================================================

! Set the definition of one per grid cell hemispheric diagnostic

      subroutine def_hdiag(k, vname_n, vname_s, lname, units, reduce, scale, sname, lwrite)

      integer (kind=int_kind), intent(in) :: &
         k               ! index in hdiag

      character (len=*), intent(in) :: &
         vname_n     , & ! netcdf variable name, northern hemisphere
         vname_s     , & ! netcdf variable name, southern hemisphere
         lname       , & ! long name
         units       , & ! units after scale is applied
         reduce          ! 'sum' area weighted global sum or 'max' global maximum

      real (kind=dbl_kind), intent(in) :: &
         scale           ! unit conversion applied to the hemispheric value

      character (len=*), intent(in), optional :: &
         sname           ! standard_name

      logical (kind=log_kind), intent(in), optional :: &
         lwrite          ! written to file, default true

      character(len=*), parameter :: subname = '(def_hdiag)'

      if (reduce /= 'sum' .and. reduce /= 'max') then
         call abort_ice(subname//' ERROR: reduce must be sum or max for '//trim(vname_n), &
            file=__FILE__, line=__LINE__)
      endif

      hdiag(k)%vname_n = vname_n
      hdiag(k)%vname_s = vname_s
      hdiag(k)%lname   = lname
      hdiag(k)%units   = units
      hdiag(k)%reduce  = reduce
      hdiag(k)%scale   = scale
      hdiag(k)%sname   = ''
      if (present(sname)) hdiag(k)%sname = sname
      hdiag(k)%lwrite  = .true.
      if (present(lwrite)) hdiag(k)%lwrite = lwrite

      end subroutine def_hdiag

!=======================================================================

! Compute the per grid cell values of the hemispheric diagnostics.
! Values to be summed are per unit area; they are multiplied by the cell
! area when reduced.

      subroutine hdiag_cell_values(fld)

      use ice_calendar, only: dt
      use ice_domain, only: nblocks
      use ice_domain_size, only: max_blocks, ncat, n_iso, n_aero
      use ice_flux, only: alvdr, alidr, alvdf, alidf, evap, fsnow, frazil, &
          fswabs, fswthru, flw, flwout, fsens, fsurf, flat, frzmlt_init, frain, fpond, &
          fhocn_ai, fsalt_ai, fresh_ai, frazil_diag, update_ocn_f, cpl_frazil, &
          swvdr, swvdf, swidr, swidf, alvdr_init, alvdf_init, alidr_init, alidf_init
      use ice_flux_bgc, only: faero_atm, faero_ocn, fiso_atm, fiso_ocn
      use ice_state, only: aice, vice, vsno, strength, aice_init, aicen, trcr, trcrn

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks,n_hdiag), intent(out) :: &
         fld               ! per grid cell values

      ! local variables

      integer (kind=int_kind) :: &
         i, j, k, n, iblk, ktherm, &
         nt_fbri, nt_apnd, nt_hpnd, nt_isosno, nt_isoice, nt_aero

      logical (kind=log_kind) :: &
         calc_Tsfc

      real (kind=dbl_kind) :: &
         puny, rhoi, awtvdr, awtidr, awtvdf, awtidf, alb

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1, work2      ! temporary

      character(len=*), parameter :: subname = '(hdiag_cell_values)'

      call icepack_query_parameters(puny_out=puny, rhoi_out=rhoi, ktherm_out=ktherm, &
           calc_Tsfc_out=calc_Tsfc, &
           awtvdr_out=awtvdr, awtidr_out=awtidr, awtvdf_out=awtvdf, awtidf_out=awtidf)
      call icepack_query_tracer_indices(nt_fbri_out=nt_fbri, nt_apnd_out=nt_apnd, &
           nt_hpnd_out=nt_hpnd, nt_isosno_out=nt_isosno, nt_isoice_out=nt_isoice, &
           nt_aero_out=nt_aero)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call diag_kinetic_energy(work1)
      call diag_ice_speed(work2)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,alb)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            fld(i,j,iblk,hd_aice) = aice(i,j,iblk)
            fld(i,j,iblk,hd_extent) = c0
            if (aice(i,j,iblk) >= aice_extmin) fld(i,j,iblk,hd_extent) = c1
            fld(i,j,iblk,hd_vice) = vice(i,j,iblk)
            fld(i,j,iblk,hd_vsno) = vsno(i,j,iblk)
            fld(i,j,iblk,hd_ke  ) = work1(i,j,iblk)
            ! ice area weighted albedo, masked where the sun is below the horizon
            alb = alvdr(i,j,iblk)*awtvdr + alidr(i,j,iblk)*awtidr &
                + alvdf(i,j,iblk)*awtvdf + alidf(i,j,iblk)*awtidf
            fld(i,j,iblk,hd_albn) = c0
            fld(i,j,iblk,hd_albd) = c0
            if (alb > puny) then
               fld(i,j,iblk,hd_albn) = aice(i,j,iblk)*alb
               fld(i,j,iblk,hd_albd) = aice(i,j,iblk)
            endif
            fld(i,j,iblk,hd_vmax) = vice(i,j,iblk)
            fld(i,j,iblk,hd_umax) = work2(i,j,iblk)
            fld(i,j,iblk,hd_pmax) = strength(i,j,iblk)
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      ! total pond volume
      if (hd_ptot > 0) then
         !$OMP PARALLEL DO PRIVATE(iblk,i,j,n)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               fld(i,j,iblk,hd_ptot) = c0
               do n = 1, ncat
                  fld(i,j,iblk,hd_ptot) = fld(i,j,iblk,hd_ptot)  &
                                        + aicen(i,j,n,iblk) &
                                        * trcrn(i,j,nt_apnd,n,iblk) &
                                        * trcrn(i,j,nt_hpnd,n,iblk)
               enddo
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
      endif

      if (.not. hdiag_global_on) return

      !-----------------------------------------------------------------
      ! conservation diagnostics, for the log (print_global) or netcdf
      ! evap, fsens, and flwout need to be multiplied by aice because
      ! regrettably they have been divided by aice for the coupler
      !-----------------------------------------------------------------

      ! total ice/snow internal energy and salt
      call total_energy (fld(:,:,:,hd_etot))
      call total_salt   (fld(:,:,:,hd_stot))

      fld(:,:,:,hd_evp) = evap(:,:,:)*aice(:,:,:)
      if (hd_brine > 0) fld(:,:,:,hd_brine) = vice(:,:,:)*trcr(:,:,nt_fbri,:)
      fld(:,:,:,hd_sfsalt) = fsalt_ai(:,:,:)
      fld(:,:,:,hd_sfresh) = fresh_ai(:,:,:)
      if (hd_spond > 0) fld(:,:,:,hd_spond) = fpond(:,:,:)

      ! ocean heat
      ! Note: fswthru not included because it does not heat ice
      fld(:,:,:,hd_fhocn) = fhocn_ai(:,:,:)

      ! atmo heat flux
      ! Note: flwout includes the reflected longwave down, needed by the
      !  atmosphere as an upwards radiative boundary condition.
      ! Also note: fswabs includes solar radiation absorbed in ocean,
      !  which must be subtracted here.
      if (calc_Tsfc) then
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               fld(i,j,iblk,hd_fhatm) = &
                  (fswabs(i,j,iblk) - fswthru(i,j,iblk) &
                 + fsens (i,j,iblk) + flwout (i,j,iblk)) &
                                    * aice      (i,j,iblk) &
                 + flw   (i,j,iblk) * aice_init (i,j,iblk)
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
      else   ! fsurf is computed by atmosphere model
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               fld(i,j,iblk,hd_fhatm) = &
                          (fsurf(i,j,iblk) - flat(i,j,iblk)) &
                           * aice(i,j,iblk)
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
      endif     ! calc_Tsfc

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            fld(i,j,iblk,hd_fswnet) = fswabs(i,j,iblk) * aice(i,j,iblk)
            fld(i,j,iblk,hd_fswdn) = (aice_init(i,j,iblk)-alvdr_init(i,j,iblk))*swvdr(i,j,iblk) &
                                   + (aice_init(i,j,iblk)-alidr_init(i,j,iblk))*swidr(i,j,iblk) &
                                   + (aice_init(i,j,iblk)-alvdf_init(i,j,iblk))*swvdf(i,j,iblk) &
                                   + (aice_init(i,j,iblk)-alidf_init(i,j,iblk))*swidf(i,j,iblk)
            ! freezing potential
            fld(i,j,iblk,hd_fhfrz) = max(c0,frzmlt_init(i,j,iblk))
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      fld(:,:,:,hd_rn) = frain(:,:,:)*aice_init(:,:,:)
      fld(:,:,:,hd_sn) = fsnow(:,:,:)*aice_init(:,:,:)

      ! frazil ice growth !! should not be multiplied by aice
      ! m/step->kg/m^2/s
      fld(:,:,:,hd_frz) = frazil(:,:,:)*rhoi/dt
      if (.not. update_ocn_f .and. ktherm == 2 .and. cpl_frazil == 'fresh_ice_correction') then
         fld(:,:,:,hd_frz) = (frazil(:,:,:)-frazil_diag(:,:,:))*rhoi/dt
      endif

      ! isotopes
      do n = 1, n_iso
         if (hd_fisoa(n) == 0) cycle
         fld(:,:,:,hd_fisoa(n)) = fiso_atm(:,:,n,:)*aice_init(:,:,:)
         fld(:,:,:,hd_fisoo(n)) = fiso_ocn(:,:,n,:)*aice(:,:,:)
         !$OMP PARALLEL DO PRIVATE(iblk,i,j,k)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               work1(i,j,iblk) = c0
               do k = 1, n_iso
                  work1(i,j,iblk) = work1(i,j,iblk) &
                                  + vsno(i,j,iblk)*trcr(i,j,nt_isosno+k-1,iblk) &
                                  + vice(i,j,iblk)*trcr(i,j,nt_isoice+k-1,iblk)
               enddo
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
         fld(:,:,:,hd_isotot(n)) = work1(:,:,:)
         fld(:,:,:,hd_isomx(n)) = work1(:,:,:)
      enddo

      ! aerosols
      do n = 1, n_aero
         if (hd_faeroa(n) == 0) cycle
         fld(:,:,:,hd_faeroa(n)) = faero_atm(:,:,n,:)*aice_init(:,:,:)
         fld(:,:,:,hd_faeroo(n)) = faero_ocn(:,:,n,:)*aice(:,:,:)
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               work1(i,j,iblk) = &
                trcr(i,j,nt_aero  +4*(n-1),iblk)*vsno(i,j,iblk) &
              + trcr(i,j,nt_aero+1+4*(n-1),iblk)*vsno(i,j,iblk) &
              + trcr(i,j,nt_aero+2+4*(n-1),iblk)*vice(i,j,iblk) &
              + trcr(i,j,nt_aero+3+4*(n-1),iblk)*vice(i,j,iblk)
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
         fld(:,:,:,hd_aerotot(n)) = work1(:,:,:)
         fld(:,:,:,hd_aeromx(n)) = work1(:,:,:)
      enddo

      end subroutine hdiag_cell_values

!=======================================================================
! Reduce the per grid cell hemispheric diagnostics to a value for each
! hemisphere, before the scale is applied.  Sums are multiplied by rn
! (eg 1/number of accumulated steps).  Called on all tasks.

      subroutine hdiag_reduce(fld, rn, valn, vals)

      use ice_constants, only: field_loc_center
      use ice_domain, only: distrb_info
      use ice_domain_size, only: max_blocks
      use ice_global_reductions, only: global_sum, global_maxval
      use ice_grid, only: lmask_n, lmask_s, tarean, tareas

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks,n_hdiag), intent(in) :: &
         fld               ! per grid cell values

      real (kind=dbl_kind), intent(in) :: &
         rn                ! factor for sums

      real (kind=dbl_kind), dimension(n_hdiag), intent(out) :: &
         valn, vals        ! northern and southern hemisphere values

      ! local variables

      integer (kind=int_kind) :: &
         k

      character(len=*), parameter :: subname = '(hdiag_reduce)'

      do k = 1, n_hdiag
         if (hdiag(k)%reduce == 'max') then
            valn(k) = global_maxval(fld(:,:,:,k), distrb_info, lmask_n)
            vals(k) = global_maxval(fld(:,:,:,k), distrb_info, lmask_s)
            if (valn(k) < maxval_spval) valn(k) = c0
            if (vals(k) < maxval_spval) vals(k) = c0
         else
            valn(k) = global_sum(fld(:,:,:,k), distrb_info, field_loc_center, tarean)*rn
            vals(k) = global_sum(fld(:,:,:,k), distrb_info, field_loc_center, tareas)*rn
         endif
      enddo

      end subroutine hdiag_reduce

!=======================================================================

! Hemispheric diagnostics derived from the reduced values for one hemisphere:
! rms ice speed from the kinetic energy and ice and snow mass, and the ice area
! weighted albedo over sunlit points

      subroutine hdiag_derived(val, urms, alb)

      use ice_constants, only: c2

      real (kind=dbl_kind), dimension(n_hdiag), intent(in) :: &
         val               ! reduced values for one hemisphere, before scale

      real (kind=dbl_kind), intent(out) :: &
         urms          , & ! rms ice speed (m/s)
         alb               ! average albedo

      ! local variables

      real (kind=dbl_kind) :: &
         puny, rhoi, rhos

      character(len=*), parameter :: subname = '(hdiag_derived)'

      call icepack_query_parameters(puny_out=puny, rhoi_out=rhoi, rhos_out=rhos)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      urms = c2*val(hd_ke)/(rhoi*val(hd_vice) + rhos*val(hd_vsno) + puny)
      if (urms > puny) then
         urms = sqrt(urms)
      else
         urms = c0
      endif

      if (val(hd_albd) > c0) then
         alb = val(hd_albn) / val(hd_albd)
      else
         alb = c0
      endif

      end subroutine hdiag_derived

!=======================================================================

! The conservation diagnostics for one hemisphere, from the reduced values.
! The changes and errors are over the last timestep, using the totals
! from the start of the timestep (init_mass_diags), so are only meaningful
! for instantaneous values.

      subroutine hdiag_global(val, hemi, dt, g)

      use ice_constants, only: p001
      use ice_domain_size, only: n_iso, n_aero
      use ice_flux, only: update_ocn_f

      real (kind=dbl_kind), dimension(n_hdiag), intent(in) :: &
         val               ! reduced values for one hemisphere, before scale

      integer (kind=int_kind), intent(in) :: &
         hemi              ! hemisphere, 1 north, 2 south

      real (kind=dbl_kind), intent(in) :: &
         dt                ! time step

      type (hdiag_global_values), intent(out) :: &
         g                 ! conservation diagnostics

      ! local variables

      integer (kind=int_kind) :: &
         n

      real (kind=dbl_kind) :: &
         rhoi, rhos, rhofresh, Lfresh, Lvap, ice_ref_salinity, &
         ptot, totm, totmi, tote, tots

      real (kind=dbl_kind), dimension(icepack_max_iso) :: &
         totiso            ! total isotope at the start of the timestep

      real (kind=dbl_kind), dimension(icepack_max_aero) :: &
         totaero           ! total aerosol at the start of the timestep

      character (len=char_len) :: &
         saltflux_option

      character(len=*), parameter :: subname = '(hdiag_global)'

      call icepack_query_parameters(rhoi_out=rhoi, rhos_out=rhos, rhofresh_out=rhofresh, &
           lfresh_out=Lfresh, lvap_out=Lvap, ice_ref_salinity_out=ice_ref_salinity, &
           saltflux_option_out=saltflux_option)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (hemi == 1) then
         totm  = totmn
         totmi = totmin
         tote  = toten
         tots  = totsn
         totiso(:)  = totison(:)
         totaero(:) = totaeron(:)
      else
         totm  = totms
         totmi = totmis
         tote  = totes
         tots  = totss
         totiso(:)  = totisos(:)
         totaero(:) = totaeros(:)
      endif

      g%etot = val(hd_etot)
      g%stot = val(hd_stot)

      g%evp = val(hd_evp)
      g%evp = g%evp*dt

      g%shmaxt = c0
      if (hd_brine > 0) g%shmaxt = val(hd_brine)

      g%sfsalt = val(hd_sfsalt)
      g%sfsalt = g%sfsalt*dt

      g%sfresh = val(hd_sfresh)
      g%sfresh = g%sfresh*dt

      g%spond = c0
      if (hd_spond > 0) then
         g%spond = val(hd_spond)
         g%spond = g%spond*dt
      endif

      g%fhocn  = val(hd_fhocn)
      g%fhatm  = val(hd_fhatm)
      g%fswnet = val(hd_fswnet)
      g%fswdn  = val(hd_fswdn)
      g%fhfrz  = val(hd_fhfrz)

      g%rn = val(hd_rn)
      g%rn = g%rn*dt
      g%sn = val(hd_sn)
      g%sn = g%sn*dt
      g%frz = val(hd_frz)
      g%frz = g%frz*dt

      ! ice, snow, pond mass
      ptot = c0
      if (hd_ptot > 0) ptot = val(hd_ptot)
      g%mice = rhoi*val(hd_vice)
      g%msnw = rhos*val(hd_vsno)
      g%mpnd = rhofresh*ptot

      ! total ice, snow and pond mass
      g%mtot = g%mice + g%msnw + g%mpnd

      ! mass change since beginning of time step
      g%delmi = g%mtot - totm

      ! ice mass change including frazil ice formation
      g%delmx = g%mice - totmi
      if (.not. update_ocn_f) then
        ! ice mass change excluding frazil ice formation
        g%delmx = g%delmx - g%frz
      endif

      ! total water flux
      g%flux = c0
      if (val(hd_aice) > c0) then
        ! water associated with frazil ice included in fresh
        g%flux = g%rn + g%sn + g%evp - g%sfresh
        if (.not. update_ocn_f) then
          g%flux = g%flux + g%frz
        endif
      endif

      g%werr = (g%flux-g%delmi)/(g%mtot + c1)

      ! energy change
      g%delei = g%etot - tote

      g%fhatm = g%fhatm + ( - g%sn * Lfresh + g%evp * Lvap ) / dt

      g%hnet = (g%fhatm - g%fhocn - g%fhfrz) * dt

      g%herr = (g%hnet - g%delei) / (g%etot - c1)

      g%swerr = (g%fswnet - g%fswdn) / (g%fswnet - c1)

      ! salt mass
      if (saltflux_option == 'prognostic') then
         ! compute the total salt mass
         g%mslt = g%stot*rhoi*p001

         ! change in salt mass
         g%delmslt = rhoi*(g%stot-tots)*p001
      else
         g%mslt = g%mice*ice_ref_salinity*p001

         ! change in salt mass
         g%delmslt = g%delmx*ice_ref_salinity*p001
      endif

      ! salt error
      g%serr = (g%sfsalt + g%delmslt) / (g%mslt + c1)

      ! isotopes
      g%fisoa  = c0
      g%fisoo  = c0
      g%isotot = c0
      g%isomx  = c0
      g%isor   = c0
      do n = 1, n_iso
         if (hd_fisoa(n) == 0) cycle
         g%fisoa(n) = val(hd_fisoa(n))
         g%fisoa(n) = g%fisoa(n)*dt
         g%fisoo(n) = val(hd_fisoo(n))
         g%fisoo(n) = g%fisoo(n)*dt
         g%isotot(n) = val(hd_isotot(n))
         g%isomx(n) = val(hd_isomx(n))
         g%isor(n) = (totiso(n)-g%isotot(n)+g%fisoa(n)-g%fisoo(n))/(g%isotot(n)+c1)
      enddo

      ! aerosols
      g%faeroa  = c0
      g%faeroo  = c0
      g%aerotot = c0
      g%aeromx  = c0
      g%aerr    = c0
      do n = 1, n_aero
         if (hd_faeroa(n) == 0) cycle
         g%faeroa(n) = val(hd_faeroa(n))
         g%faeroa(n) = g%faeroa(n)*dt
         g%faeroo(n) = val(hd_faeroo(n))
         g%faeroo(n) = g%faeroo(n)*dt
         g%aerotot(n) = val(hd_aerotot(n))
         g%aeromx(n) = val(hd_aeromx(n))
         g%aerr(n) = (totaero(n)-g%aerotot(n)+g%faeroa(n)-g%faeroo(n)) &
                   / (g%aerotot(n) + c1)
      enddo

      end subroutine hdiag_global

!=======================================================================

      end module ice_diagnostics

!=======================================================================
