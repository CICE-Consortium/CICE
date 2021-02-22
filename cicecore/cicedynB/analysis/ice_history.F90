!=======================================================================
! Driver for core history output
!
! The following variables are currently hard-wired as snapshots 
!   (instantaneous rather than time-averages):
!   divu, shear, sig1, sig2, sigP, trsig, mlt_onset, frz_onset, hisnap, aisnap
!
! Options for histfreq: '1','h','d','m','y','x', where x means that
!   output stream will not be used (recommended for efficiency).  
! histfreq_n can be any nonnegative integer, where 0 means that the 
!   corresponding histfreq frequency will not be used.
! The flags (f_<field>) can be set to '1','h','d','m','y' or 'x', where
!   n means the field will not be written.  To output the same field at
!   more than one frequency, for instance monthy and daily, set 
!   f_<field> = 'md'.
!
! authors Tony Craig and Bruce Briegleb, NCAR
!         Elizabeth C. Hunke and William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! 2004 WHL: Block structure added 
! 2006 ECH: Accepted some CESM code into mainstream CICE
!           Added ice_present, aicen, vicen; removed aice1...10, vice1...1.
!           Added histfreq_n and histfreq='h' options, removed histfreq='w'
!           Converted to free source form (F90)
!           Added option for binary output instead of netCDF
! 2009 D Bailey and ECH: Generalized for multiple frequency output
! 2010 Alison McLaren and ECH: Added 3D capability

      module ice_history

      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0, c1, c2, c100, c360, c180, &
          p001, p25, p5, mps_to_cmpdy, kg_to_g, spval
      use ice_fileunits, only: nu_nml, nml_filename, nu_diag, &
          get_fileunit, release_fileunit, flush_fileunit
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_snow_temperature, icepack_ice_temperature
      use icepack_intfc, only: icepack_query_parameters, &
          icepack_query_tracer_flags, icepack_query_tracer_indices

      implicit none
      private
      public :: init_hist, accum_hist
      
!=======================================================================

      contains

!=======================================================================

! Initialize history files
!
! authors Tony Craig, NCAR
!         Elizabeth C. Hunke, LANL
!         C.M. Bitz, UW
!         Bruce P. Briegleb, NCAR
!         William H. Lipscomb, LANL

      subroutine init_hist (dt)

      use ice_blocks, only: nx_block, ny_block
      use ice_broadcast, only: broadcast_scalar, broadcast_array
      use ice_calendar, only: yday, days_per_year, histfreq, &
          histfreq_n, nstreams
      use ice_domain_size, only: max_blocks, max_nstrm, nilyr, nslyr, nblyr, ncat, nfsd
      use ice_dyn_shared, only: kdyn
      use ice_flux, only: mlt_onset, frz_onset, albcnt
      use ice_history_shared ! everything
      use ice_history_mechred, only: init_hist_mechred_2D, init_hist_mechred_3Dc
      use ice_history_pond, only: init_hist_pond_2D, init_hist_pond_3Dc
      use ice_history_bgc, only:init_hist_bgc_2D, init_hist_bgc_3Dc, &
          init_hist_bgc_3Db, init_hist_bgc_3Da
      use ice_history_drag, only: init_hist_drag_2D
      use ice_history_fsd, only: init_hist_fsd_2D, init_hist_fsd_3Df, &
          init_hist_fsd_4Df, f_afsd, f_afsdn
      use ice_restart_shared, only: restart

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      real (kind=dbl_kind) :: rhofresh, Tffresh, secday, rad_to_deg
      logical (kind=log_kind) :: formdrag
      logical (kind=log_kind) :: tr_iage, tr_FY, tr_lvl, tr_pond, tr_aero, tr_brine
      logical (kind=log_kind) :: tr_fsd
      logical (kind=log_kind) :: skl_bgc, solve_zsal, solve_zbgc, z_tracers
      integer (kind=int_kind) :: n, ns, ns1, ns2
      integer (kind=int_kind), dimension(max_nstrm) :: &
         ntmp
      integer (kind=int_kind) :: nml_error ! namelist i/o error flag
      character(len=*), parameter :: subname = '(init_hist)'

      !-----------------------------------------------------------------
      ! set history dimensions
      !-----------------------------------------------------------------

      ncat_hist = ncat    ! number of thickness categories written <= ncat
      nfsd_hist = nfsd    ! number of floe size categories written <= nfsd
      nzilyr = nilyr      ! vertical dimension (allows alternative grids)
      nzslyr = nslyr      ! snow
      nzblyr = nblyr+2    ! bio grid
      nzalyr = nblyr+4    ! aerosols (2 snow & nblyr+2 bio)

      !-----------------------------------------------------------------
      ! read namelist
      !-----------------------------------------------------------------

      call icepack_query_parameters(rhofresh_out=rhofresh, Tffresh_out=Tffresh, &
         secday_out=secday, rad_to_deg_out=rad_to_deg)
      call icepack_query_parameters(formdrag_out=formdrag, skl_bgc_out=skl_bgc, &
         solve_zsal_out=solve_zsal, solve_zbgc_out=solve_zbgc, z_tracers_out=z_tracers)
      call icepack_query_tracer_flags(tr_iage_out=tr_iage, tr_FY_out=tr_FY, &
         tr_lvl_out=tr_lvl, tr_pond_out=tr_pond, tr_aero_out=tr_aero, &
         tr_brine_out=tr_brine, tr_fsd_out=tr_fsd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call get_fileunit(nu_nml)
      if (my_task == master_task) then
         open (nu_nml, file=nml_filename, status='old',iostat=nml_error)
         if (nml_error /= 0) then
            nml_error = -1
         else
            nml_error =  1
         endif
         do while (nml_error > 0)
            read(nu_nml, nml=icefields_nml,iostat=nml_error)
         end do
         if (nml_error == 0) close(nu_nml)
      endif
      call release_fileunit(nu_nml)

      call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
         close (nu_nml)
         call abort_ice(subname//'ERROR: reading icefields_nml')
      endif

      ! histfreq options ('1','h','d','m','y')
      nstreams = 0
      do ns = 1, max_nstrm
         if (histfreq(ns) == '1' .or. histfreq(ns) == 'h' .or. &
             histfreq(ns) == 'd' .or. histfreq(ns) == 'm' .or. &
             histfreq(ns) == 'y') then
                nstreams = nstreams + 1
                if (ns >= 2) then
                   if (histfreq(ns-1) == 'x') then
                      call abort_ice(subname//'ERROR: histfreq all non x must be at start of array')
                   endif
                endif
         else if (histfreq(ns) /= 'x') then
             call abort_ice(subname//'ERROR: histfreq contains illegal element')
         endif
      enddo
      if (nstreams == 0) write (nu_diag,*) 'WARNING: No history output'
      do ns1 = 1, nstreams
         do ns2 = 1, nstreams
            if (histfreq(ns1) == histfreq(ns2) .and. ns1/=ns2 &
               .and. my_task == master_task) then
               call abort_ice(subname//'ERROR: histfreq elements must be unique')
            endif
         enddo
      enddo

      if (.not. tr_iage) then
         f_iage = 'x'
         f_dagedtt = 'x'
         f_dagedtd = 'x'
      endif
      if (.not. tr_FY)   f_FY   = 'x'
      if (kdyn /= 2) then
           f_a11       = 'x'
           f_a12       = 'x'
           f_e11       = 'x'
           f_e12       = 'x'
           f_e22       = 'x'
           f_s11       = 'x'
           f_s12       = 'x'
           f_s22       = 'x'
           f_yieldstress11 = 'x'
           f_yieldstress12 = 'x'
           f_yieldstress22 = 'x'
      endif

      ! these must be output at the same frequency because of 
      ! cos(zenith angle) averaging
      if (f_albice(1:1) /= 'x' .and. f_albsni(1:1) /= 'x') f_albice = f_albsni
      if (f_albsno(1:1) /= 'x') f_albsno = f_albice
      if (f_albpnd(1:1) /= 'x') f_albpnd = f_albice
      if (f_coszen(1:1) /= 'x' .and. f_albice(1:1) /= 'x') f_coszen = f_albice
      if (f_coszen(1:1) /= 'x' .and. f_albsni(1:1) /= 'x') f_coszen = f_albsni

      ! to prevent array-out-of-bounds when aggregating
      if (f_fmeltt_ai(1:1) /= 'x') f_fmelttn_ai = f_fmeltt_ai

     ! Turn on all CMIP fields in one go.

      if (f_CMIP(1:1) /= 'x') then
         f_sithick = 'mxxxx'
         f_sisnthick = 'mxxxx'
         f_siage = 'mxxxx'
         f_sitemptop = 'mxxxx'
         f_sitempsnic = 'mxxxx'
         f_sitempbot = 'mxxxx'
         f_sispeed = 'mxxxx'
         f_siu = 'mxxxx'
         f_siv = 'mxxxx'
         f_sidmasstranx = 'mxxxx'
         f_sidmasstrany = 'mxxxx'
         f_sistrxdtop = 'mxxxx'
         f_sistrydtop = 'mxxxx'
         f_sistrxubot = 'mxxxx'
         f_sistryubot = 'mxxxx'
         f_sicompstren = 'mxxxx'
         f_sialb = 'mxxxx'
         f_sihc = 'mxxxx'
         f_sisnhc = 'mxxxx'
         f_sidconcth = 'mxxxx'
         f_sidconcdyn = 'mxxxx'
         f_sidmassth = 'mxxxx'
         f_sidmassdyn = 'mxxxx'
         f_sidmassgrowthwat = 'mxxxx'
         f_sidmassgrowthbot = 'mxxxx'
         f_sidmasssi = 'mxxxx'
         f_sidmassevapsubl = 'mxxxx'
         f_sndmasssubl = 'mxxxx'
         f_sidmassmelttop = 'mxxxx'
         f_sidmassmeltbot = 'mxxxx'
         f_sidmasslat = 'mxxxx'
         f_sndmasssnf = 'mxxxx'
         f_sndmassmelt = 'mxxxx'
         f_siflswdtop = 'mxxxx'
         f_siflswutop = 'mxxxx'
         f_siflswdbot = 'mxxxx'
         f_sifllwdtop = 'mxxxx'
         f_sifllwutop = 'mxxxx'
         f_siflsenstop = 'mxxxx'
         f_siflsensupbot = 'mxxxx'
         f_sifllatstop = 'mxxxx'
         f_siflcondtop = 'mxxxx'
         f_siflcondbot = 'mxxxx'
         f_sipr = 'mxxxx'
         f_sifb = 'mxxxx'
         f_siflsaltbot = 'mxxxx'
         f_siflfwbot = 'mxxxx'
         f_siflfwdrain = 'mxxxx'
         f_siforcetiltx = 'mxxxx'
         f_siforcetilty = 'mxxxx'
         f_siforcecoriolx = 'mxxxx'
         f_siforcecorioly = 'mxxxx'
         f_siforceintstrx = 'mxxxx'
         f_siforceintstry = 'mxxxx'
         f_sidragtop = 'mxxxx'
         f_sistreave = 'mxxxx'
         f_sistremax = 'mxxxx'
         f_sirdgthick = 'mxxxx'
         f_siitdconc = 'mxxxx'
         f_siitdthick = 'mxxxx'
         f_siitdsnthick = 'mxxxx'
         f_aicen = 'mxxxx'
      endif

      if (f_CMIP(2:2) == 'd') then
         f_icepresent = f_CMIP
         f_aice = f_CMIP
         f_sithick = f_CMIP
         f_sisnthick = f_CMIP
         f_sitemptop = f_CMIP
         f_siu = f_CMIP
         f_siv = f_CMIP
         f_sispeed = f_CMIP
      endif

#ifndef ncdf
      f_bounds = .false.
#endif

      ! write dimensions for 3D or 4D history variables
      ! note: list of variables checked here is incomplete
      if (f_aicen(1:1) /= 'x' .or. f_vicen(1:1) /= 'x' .or. &
          f_Tinz (1:1) /= 'x' .or. f_Sinz (1:1) /= 'x') f_NCAT  = .true.
      if (f_Tinz (1:1) /= 'x' .or. f_Sinz (1:1) /= 'x') f_VGRDi = .true.
      if (f_Tsnz (1:1) /= 'x')                          f_VGRDs = .true.
      if (tr_fsd)                                       f_NFSD  = .true.

      call broadcast_scalar (f_tmask, master_task)
      call broadcast_scalar (f_blkmask, master_task)
      call broadcast_scalar (f_tarea, master_task)
      call broadcast_scalar (f_uarea, master_task)
      call broadcast_scalar (f_dxt, master_task)
      call broadcast_scalar (f_dyt, master_task)
      call broadcast_scalar (f_dxu, master_task)
      call broadcast_scalar (f_dyu, master_task)
      call broadcast_scalar (f_HTN, master_task)
      call broadcast_scalar (f_HTE, master_task)
      call broadcast_scalar (f_ANGLE, master_task)
      call broadcast_scalar (f_ANGLET, master_task)
      call broadcast_scalar (f_bounds, master_task)
      call broadcast_scalar (f_NCAT, master_task)
      call broadcast_scalar (f_VGRDi, master_task)
      call broadcast_scalar (f_VGRDs, master_task)
      call broadcast_scalar (f_VGRDb, master_task)
      call broadcast_scalar (f_VGRDa, master_task)
      call broadcast_scalar (f_NFSD, master_task)

!     call broadcast_scalar (f_example, master_task)
      call broadcast_scalar (f_hi, master_task)
      call broadcast_scalar (f_hs, master_task)
      call broadcast_scalar (f_snowfrac, master_task)
      call broadcast_scalar (f_snowfracn, master_task)
      call broadcast_scalar (f_Tsfc, master_task)
      call broadcast_scalar (f_aice, master_task)
      call broadcast_scalar (f_uvel, master_task)
      call broadcast_scalar (f_vvel, master_task)
      call broadcast_scalar (f_uatm, master_task)
      call broadcast_scalar (f_vatm, master_task)
      call broadcast_scalar (f_atmspd, master_task)
      call broadcast_scalar (f_atmdir, master_task)
      call broadcast_scalar (f_sice, master_task)
      call broadcast_scalar (f_fswup, master_task)
      call broadcast_scalar (f_fswdn, master_task)
      call broadcast_scalar (f_flwdn, master_task)
      call broadcast_scalar (f_snow, master_task)
      call broadcast_scalar (f_snow_ai, master_task)
      call broadcast_scalar (f_rain, master_task)
      call broadcast_scalar (f_rain_ai, master_task)
      call broadcast_scalar (f_sst, master_task)
      call broadcast_scalar (f_sss, master_task)
      call broadcast_scalar (f_uocn, master_task)
      call broadcast_scalar (f_vocn, master_task)
      call broadcast_scalar (f_ocnspd, master_task)
      call broadcast_scalar (f_ocndir, master_task)
      call broadcast_scalar (f_frzmlt, master_task)
      call broadcast_scalar (f_fswfac, master_task)
      call broadcast_scalar (f_fswint_ai, master_task)
      call broadcast_scalar (f_fswabs, master_task)
      call broadcast_scalar (f_fswabs_ai, master_task)
      call broadcast_scalar (f_albsni, master_task)
      call broadcast_scalar (f_alvdr, master_task)
      call broadcast_scalar (f_alidr, master_task)
      call broadcast_scalar (f_alvdf, master_task)
      call broadcast_scalar (f_alidf, master_task)
      call broadcast_scalar (f_alvdr_ai, master_task)
      call broadcast_scalar (f_alidr_ai, master_task)
      call broadcast_scalar (f_alvdf_ai, master_task)
      call broadcast_scalar (f_alidf_ai, master_task)
      call broadcast_scalar (f_albice, master_task)
      call broadcast_scalar (f_albsno, master_task)
      call broadcast_scalar (f_albpnd, master_task)
      call broadcast_scalar (f_coszen, master_task)
      call broadcast_scalar (f_flat, master_task)
      call broadcast_scalar (f_flat_ai, master_task)
      call broadcast_scalar (f_fsens, master_task)
      call broadcast_scalar (f_fsens_ai, master_task)
      call broadcast_scalar (f_flwup, master_task)
      call broadcast_scalar (f_flwup_ai, master_task)
      call broadcast_scalar (f_evap, master_task)
      call broadcast_scalar (f_evap_ai, master_task)
      call broadcast_scalar (f_Tair, master_task)
      call broadcast_scalar (f_Tref, master_task)
      call broadcast_scalar (f_Qref, master_task)
      call broadcast_scalar (f_congel, master_task)
      call broadcast_scalar (f_frazil, master_task)
      call broadcast_scalar (f_snoice, master_task)
      call broadcast_scalar (f_dsnow, master_task)
      call broadcast_scalar (f_meltt, master_task)
      call broadcast_scalar (f_melts, master_task)
      call broadcast_scalar (f_meltb, master_task)
      call broadcast_scalar (f_meltl, master_task)
      call broadcast_scalar (f_fresh, master_task)
      call broadcast_scalar (f_fresh_ai, master_task)
      call broadcast_scalar (f_fsalt, master_task)
      call broadcast_scalar (f_fsalt_ai, master_task)
      call broadcast_scalar (f_fbot, master_task)
      call broadcast_scalar (f_fhocn, master_task)
      call broadcast_scalar (f_fhocn_ai, master_task)
      call broadcast_scalar (f_fswthru, master_task)
      call broadcast_scalar (f_fswthru_ai, master_task)
      call broadcast_scalar (f_strairx, master_task)
      call broadcast_scalar (f_strairy, master_task)
      call broadcast_scalar (f_strtltx, master_task)
      call broadcast_scalar (f_strtlty, master_task)
      call broadcast_scalar (f_strcorx, master_task)
      call broadcast_scalar (f_strcory, master_task)
      call broadcast_scalar (f_strocnx, master_task)
      call broadcast_scalar (f_strocny, master_task)
      call broadcast_scalar (f_strintx, master_task)
      call broadcast_scalar (f_strinty, master_task)
      call broadcast_scalar (f_taubx, master_task)
      call broadcast_scalar (f_tauby, master_task)
      call broadcast_scalar (f_strength, master_task)
      call broadcast_scalar (f_divu, master_task)
      call broadcast_scalar (f_shear, master_task)
      call broadcast_scalar (f_sig1, master_task)
      call broadcast_scalar (f_sig2, master_task)
      call broadcast_scalar (f_sigP, master_task)
      call broadcast_scalar (f_dvidtt, master_task)
      call broadcast_scalar (f_dvidtd, master_task)
      call broadcast_scalar (f_daidtt, master_task)
      call broadcast_scalar (f_daidtd, master_task)
      call broadcast_scalar (f_dagedtt, master_task)
      call broadcast_scalar (f_dagedtd, master_task)
      call broadcast_scalar (f_mlt_onset, master_task)
      call broadcast_scalar (f_frz_onset, master_task)
      call broadcast_scalar (f_aisnap, master_task)
      call broadcast_scalar (f_hisnap, master_task)
      call broadcast_scalar (f_sithick, master_task)
      call broadcast_scalar (f_siage, master_task)
      call broadcast_scalar (f_sisnthick, master_task)
      call broadcast_scalar (f_sitemptop, master_task)
      call broadcast_scalar (f_sitempsnic, master_task)
      call broadcast_scalar (f_sitempbot, master_task)
      call broadcast_scalar (f_siu, master_task)
      call broadcast_scalar (f_siv, master_task)
      call broadcast_scalar (f_sidmasstranx, master_task)
      call broadcast_scalar (f_sidmasstrany, master_task)
      call broadcast_scalar (f_sistrxdtop, master_task)
      call broadcast_scalar (f_sistrydtop, master_task)
      call broadcast_scalar (f_sistrxubot, master_task)
      call broadcast_scalar (f_sistryubot, master_task)
      call broadcast_scalar (f_sicompstren, master_task)
      call broadcast_scalar (f_sispeed, master_task)
      call broadcast_scalar (f_sidir, master_task)
      call broadcast_scalar (f_sialb, master_task)
      call broadcast_scalar (f_sihc, master_task)
      call broadcast_scalar (f_sisnhc, master_task)
      call broadcast_scalar (f_sidconcth, master_task)
      call broadcast_scalar (f_sidconcdyn, master_task)
      call broadcast_scalar (f_sidmassth, master_task)
      call broadcast_scalar (f_sidmassdyn, master_task)
      call broadcast_scalar (f_sidmassgrowthwat, master_task)
      call broadcast_scalar (f_sidmassgrowthbot, master_task)
      call broadcast_scalar (f_sidmasssi, master_task)
      call broadcast_scalar (f_sidmassevapsubl, master_task)
      call broadcast_scalar (f_sndmasssubl, master_task)
      call broadcast_scalar (f_sidmassmelttop, master_task)
      call broadcast_scalar (f_sidmassmeltbot, master_task)
      call broadcast_scalar (f_sidmasslat, master_task)
      call broadcast_scalar (f_sndmasssnf, master_task)
      call broadcast_scalar (f_sndmassmelt, master_task)
      call broadcast_scalar (f_siflswdtop, master_task)
      call broadcast_scalar (f_siflswutop, master_task)
      call broadcast_scalar (f_siflswdbot, master_task)
      call broadcast_scalar (f_sifllwdtop, master_task)
      call broadcast_scalar (f_sifllwutop, master_task)
      call broadcast_scalar (f_siflsenstop, master_task)
      call broadcast_scalar (f_siflsensupbot, master_task)
      call broadcast_scalar (f_sifllatstop, master_task)
      call broadcast_scalar (f_siflcondtop, master_task)
      call broadcast_scalar (f_siflcondbot, master_task)
      call broadcast_scalar (f_sipr, master_task)
      call broadcast_scalar (f_sifb, master_task)
      call broadcast_scalar (f_siflsaltbot, master_task)
      call broadcast_scalar (f_siflfwbot, master_task)
      call broadcast_scalar (f_siflfwdrain, master_task)
      call broadcast_scalar (f_siforcetiltx, master_task)
      call broadcast_scalar (f_siforcetilty, master_task)
      call broadcast_scalar (f_siforcecoriolx, master_task)
      call broadcast_scalar (f_siforcecorioly, master_task)
      call broadcast_scalar (f_siforceintstrx, master_task)
      call broadcast_scalar (f_siforceintstry, master_task)
      call broadcast_scalar (f_siitdconc, master_task)
      call broadcast_scalar (f_siitdthick, master_task)
      call broadcast_scalar (f_siitdsnthick, master_task)
      call broadcast_scalar (f_sidragtop, master_task)
      call broadcast_scalar (f_sistreave, master_task)
      call broadcast_scalar (f_sistremax, master_task)
      call broadcast_scalar (f_sirdgthick, master_task)

      call broadcast_scalar (f_aicen, master_task)
      call broadcast_scalar (f_vicen, master_task)
      call broadcast_scalar (f_vsnon, master_task)
      call broadcast_scalar (f_trsig, master_task)
      call broadcast_scalar (f_icepresent, master_task)
      call broadcast_scalar (f_fsurf_ai, master_task)
      call broadcast_scalar (f_fcondtop_ai, master_task)
      call broadcast_scalar (f_fmeltt_ai, master_task)
      call broadcast_scalar (f_fsurfn_ai, master_task)
      call broadcast_scalar (f_fcondtopn_ai, master_task)
      call broadcast_scalar (f_fmelttn_ai, master_task)
      call broadcast_scalar (f_flatn_ai, master_task)
      call broadcast_scalar (f_fsensn_ai, master_task)

!      call broadcast_scalar (f_field3dz, master_task)
      call broadcast_scalar (f_keffn_top, master_task)
      call broadcast_scalar (f_Tinz, master_task)
      call broadcast_scalar (f_Sinz, master_task)
      call broadcast_scalar (f_Tsnz, master_task)

      call broadcast_scalar (f_iage, master_task)
      call broadcast_scalar (f_FY, master_task)

      call broadcast_scalar (f_a11, master_task)
      call broadcast_scalar (f_a12, master_task)
      call broadcast_scalar (f_e11, master_task)
      call broadcast_scalar (f_e12, master_task)
      call broadcast_scalar (f_e22, master_task)
      call broadcast_scalar (f_s11, master_task)
      call broadcast_scalar (f_s12, master_task) 
      call broadcast_scalar (f_s22, master_task)
      call broadcast_scalar (f_yieldstress11, master_task)
      call broadcast_scalar (f_yieldstress12, master_task)
      call broadcast_scalar (f_yieldstress22, master_task)

      ! 2D variables
      do ns1 = 1, nstreams
      if (histfreq(ns1) /= 'x') then

!!!!! begin example
!         call define_hist_field(n_example,"example","m",tstr2D, tcstr, & 
!            "example: mean ice thickness",                           &
!            "ice volume per unit grid cell area", c1, c0,            &
!            ns1, f_example)
!!!!! end example

         call define_hist_field(n_hi,"hi","m",tstr2D, tcstr,        & 
            "grid cell mean ice thickness",                       &
            "ice volume per unit grid cell area", c1, c0,         &
            ns1, f_hi)

         call define_hist_field(n_hs,"hs","m",tstr2D, tcstr,        &
             "grid cell mean snow thickness",                     &
             "snow volume per unit grid cell area", c1, c0,       &
             ns1, f_hs)

         call define_hist_field(n_snowfrac,"snowfrac","1",tstr2D, tcstr, &
             "grid cell mean snow fraction",                     &
             "snow fraction per unit grid cell area", c1, c0,       &
             ns1, f_snowfrac)

         call define_hist_field(n_Tsfc,"Tsfc","C",tstr2D, tcstr,    &
             "snow/ice surface temperature",                      &
             "averaged with Tf if no ice is present", c1, c0,     &
             ns1, f_Tsfc)
      
         call define_hist_field(n_aice,"aice","1",tstr2D, tcstr,    &
             "ice area  (aggregate)",                             &
             "none", c1, c0,                                      &
             ns1, f_aice)
      
         call define_hist_field(n_uvel,"uvel","m/s",ustr2D, ucstr,  &
             "ice velocity (x)",                                  &
             "positive is x direction on U grid", c1, c0,         &
             ns1, f_uvel)
      
         call define_hist_field(n_vvel,"vvel","m/s",ustr2D, ucstr,  &
             "ice velocity (y)",                                  &
             "positive is y direction on U grid", c1, c0,         &
             ns1, f_vvel)
      
         call define_hist_field(n_uatm,"uatm","m/s",ustr2D, ucstr,  &
             "atm velocity (x)",                                  &
             "positive is x direction on U grid", c1, c0,         &
             ns1, f_uatm)
      
         call define_hist_field(n_vatm,"vatm","m/s",ustr2D, ucstr,  &
             "atm velocity (y)",                                  &
             "positive is y direction on U grid", c1, c0,         &
             ns1, f_vatm)

         call define_hist_field(n_atmspd,"atmspd","m/s",ustr2D, ucstr, &
             "atmosphere wind speed",                                  &
             "vector magnitude", c1, c0,                               &
             ns1, f_atmspd)
      
         call define_hist_field(n_atmdir,"atmdir","deg",ustr2D, ucstr, &
             "atmosphere wind direction",                              &
             "vector direction - coming from", c1, c0,                 &
             ns1, f_atmdir)
      
         call define_hist_field(n_sice,"sice","ppt",tstr2D, tcstr,  &
             "bulk ice salinity",                                 &
             "none", c1, c0,                                      &
             ns1, f_sice)
      
         call define_hist_field(n_fswup,"fswup","W/m^2",tstr2D, tcstr, &
             "upward solar flux",                                      &
             "positive upward", c1, c0,                            &
             ns1, f_fswup)
      
         call define_hist_field(n_fswdn,"fswdn","W/m^2",tstr2D, tcstr, &
             "down solar flux",                                      &
             "positive downward", c1, c0,                            &
             ns1, f_fswdn)
      
         call define_hist_field(n_flwdn,"flwdn","W/m^2",tstr2D, tcstr, &
             "down longwave flux",                                   &
             "positive downward", c1, c0,                            &
             ns1, f_flwdn)
      
         call define_hist_field(n_snow,"snow","cm/day",tstr2D, tcstr, &
             "snowfall rate (cpl)",                                 &
             "none", mps_to_cmpdy/rhofresh, c0,                     &
             ns1, f_snow)
      
         call define_hist_field(n_snow_ai,"snow_ai","cm/day",tstr2D, tcstr, &
             "snowfall rate",                                             &
             "weighted by ice area", mps_to_cmpdy/rhofresh, c0,           &
             ns1, f_snow_ai)
      
         call define_hist_field(n_rain,"rain","cm/day",tstr2D, tcstr, &
             "rainfall rate (cpl)",                                 &
             "none", mps_to_cmpdy/rhofresh, c0,                     &
             ns1, f_rain)
      
         call define_hist_field(n_rain_ai,"rain_ai","cm/day",tstr2D, tcstr, &
             "rainfall rate",                                             &
             "weighted by ice area", mps_to_cmpdy/rhofresh, c0,           &
             ns1, f_rain_ai)
      
         call define_hist_field(n_sst,"sst","C",tstr2D, tcstr, &
             "sea surface temperature",                      &
             "none", c1, c0,                                 &
             ns1, f_sst)
      
         call define_hist_field(n_sss,"sss","ppt",tstr2D, tcstr, &
             "sea surface salinity",                           &
             "none", c1, c0,                                   &
             ns1, f_sss)
      
         call define_hist_field(n_uocn,"uocn","m/s",ustr2D, ucstr, &
             "ocean current (x)",                                &
             "positive is x direction on U grid", c1, c0,        &
             ns1, f_uocn)
      
         call define_hist_field(n_vocn,"vocn","m/s",ustr2D, ucstr, &
             "ocean current (y)",                                &
             "positive is y direction on U grid", c1, c0,        &
             ns1, f_vocn)

         call define_hist_field(n_ocnspd,"ocnspd","m/s",ustr2D, ucstr, &
             "ocean current speed",                                    &
             "vector magnitude", c1, c0,                               &
             ns1, f_ocnspd)
      
         call define_hist_field(n_ocndir,"ocndir","deg",ustr2D, ucstr, &
             "ocean current direction",                                &
             "vector direction - going to", c1, c0,                    &
             ns1, f_ocndir)
      
         call define_hist_field(n_frzmlt,"frzmlt","W/m^2",tstr2D, tcstr, &
             "freeze/melt potential",                                  &
             "if >0, new ice forms; if <0, ice melts", c1, c0,         &
             ns1, f_frzmlt)
      
         call define_hist_field(n_fswfac,"scale_factor","1",tstr2D, tcstr, &
             "shortwave scaling factor",                           &
             "ratio of netsw new:old", c1, c0,                     &
             ns1, f_fswfac)

         call define_hist_field(n_fswint_ai,"fswint_ai","W/m^2",tstr2D, tcstr, &
             "shortwave absorbed in ice interior",                             &
             "does not include surface", c1, c0,                               &
             ns1, f_fswint_ai)

         call define_hist_field(n_fswabs,"fswabs","W/m^2",tstr2D, tcstr, &
             "snow/ice/ocn absorbed solar flux (cpl)",                 &
             "positive downward", c1, c0,                              &
             ns1, f_fswabs)
      
         call define_hist_field(n_fswabs_ai,"fswabs_ai","W/m^2",tstr2D, tcstr, &
             "snow/ice/ocn absorbed solar flux",                             &
             "weighted by ice area", c1, c0,                                 &
             ns1, f_fswabs_ai)
      
         call define_hist_field(n_albsni,"albsni","%",tstr2D, tcstr, &
             "snow/ice broad band albedo",                         &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_albsni)
      
         call define_hist_field(n_alvdr,"alvdr","%",tstr2D, tcstr, &
             "visible direct albedo",                            &
             "scaled (divided) by aice", c100, c0,               &
             ns1, f_alvdr)
      
         call define_hist_field(n_alidr,"alidr","%",tstr2D, tcstr, &
             "near IR direct albedo",                            &
             "scaled (divided) by aice", c100, c0,               &
             ns1, f_alidr)

         call define_hist_field(n_alvdf,"alvdf","%",tstr2D, tcstr, &
             "visible diffuse albedo",                            &
             "scaled (divided) by aice", c100, c0,               &
             ns1, f_alvdf)
      
         call define_hist_field(n_alidf,"alidf","%",tstr2D, tcstr, &
             "near IR diffuse albedo",                            &
             "scaled (divided) by aice", c100, c0,               &
             ns1, f_alidf)

         call define_hist_field(n_alvdr_ai,"alvdr_ai","%",tstr2D, tcstr, &
             "visible direct albedo",                            &
             " ", c100, c0,               &
             ns1, f_alvdr_ai)
      
         call define_hist_field(n_alidr_ai,"alidr_ai","%",tstr2D, tcstr, &
             "near IR direct albedo",                            &
             " ", c100, c0,               &
             ns1, f_alidr_ai)

         call define_hist_field(n_alvdf_ai,"alvdf_ai","%",tstr2D, tcstr, &
             "visible diffuse albedo",                            &
             " ", c100, c0,               &
             ns1, f_alvdf_ai)
      
         call define_hist_field(n_alidf_ai,"alidf_ai","%",tstr2D, tcstr, &
             "near IR diffuse albedo",                            &
             " ", c100, c0,               &
             ns1, f_alidf_ai)

         call define_hist_field(n_albice,"albice","%",tstr2D, tcstr, &
             "bare ice albedo",                                    &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_albice)
      
         call define_hist_field(n_albsno,"albsno","%",tstr2D, tcstr, &
             "snow albedo",                                        &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_albsno)
      
         call define_hist_field(n_albpnd,"albpnd","%",tstr2D, tcstr, &
             "melt pond albedo",                                   &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_albpnd)
      
         call define_hist_field(n_coszen,"coszen","radian",tstr2D, tcstr, &
             "cosine of the zenith angle",                              &
             "negative below horizon", c1, c0,                          &
             ns1, f_coszen)

         call define_hist_field(n_flat,"flat","W/m^2",tstr2D, tcstr, &
             "latent heat flux (cpl)",                             &
             "positive downward", c1, c0,                          &
             ns1, f_flat)
      
         call define_hist_field(n_flat_ai,"flat_ai","W/m^2",tstr2D, tcstr, &
             "latent heat flux",                                         &
             "weighted by ice area", c1, c0,                             &
             ns1, f_flat_ai)
      
         call define_hist_field(n_fsens,"fsens","W/m^2",tstr2D, tcstr, &
             "sensible heat flux (cpl)",                             &
             "positive downward", c1, c0,                            &
             ns1, f_fsens)
      
         call define_hist_field(n_fsens_ai,"fsens_ai","W/m^2",tstr2D, tcstr, &
             "sensible heat flux",                                         &
             "weighted by ice area", c1, c0,                               &
             ns1, f_fsens_ai)
      
         call define_hist_field(n_flwup,"flwup","W/m^2",tstr2D, tcstr, &
             "upward longwave flux (cpl)",                           &
             "positive downward", c1, c0,                            &
             ns1, f_flwup)
      
         call define_hist_field(n_flwup_ai,"flwup_ai","W/m^2",tstr2D, tcstr, &
             "upward longwave flux",                                       &
             "weighted by ice area", c1, c0,                               &
             ns1, f_flwup_ai)
      
         call define_hist_field(n_evap,"evap","cm/day",tstr2D, tcstr, &
             "evaporative water flux (cpl)",                        &
             "none", mps_to_cmpdy/rhofresh, c0,                     &
             ns1, f_evap)
      
         call define_hist_field(n_evap_ai,"evap_ai","cm/day",tstr2D, tcstr, &
             "evaporative water flux",                                    &
             "weighted by ice area", mps_to_cmpdy/rhofresh, c0,           &
             ns1, f_evap_ai)
      
         call define_hist_field(n_Tair,"Tair","C",tstr2D, tcstr, &
             "air temperature",                                &
             "none", c1, -Tffresh,                             &
             ns1, f_Tair)
      
         call define_hist_field(n_Tref,"Tref","C",tstr2D, tcstr, &
             "2m reference temperature",                       &
             "none", c1, -Tffresh,                             &
             ns1, f_Tref)
      
         call define_hist_field(n_Qref,"Qref","g/kg",tstr2D, tcstr, &
             "2m reference specific humidity",                    &
             "none", kg_to_g, c0,                                 &
             ns1, f_Qref)
      
         call define_hist_field(n_congel,"congel","cm/day",tstr2D, tcstr, &
             "congelation ice growth",                                  &
             "none", mps_to_cmpdy/dt, c0,                               &
             ns1, f_congel)
      
         call define_hist_field(n_frazil,"frazil","cm/day",tstr2D, tcstr, &
             "frazil ice growth",                                       &
             "none", mps_to_cmpdy/dt, c0,                               &
             ns1, f_frazil)
      
         call define_hist_field(n_snoice,"snoice","cm/day",tstr2D, tcstr, &
             "snow-ice formation",                                      &
             "none", mps_to_cmpdy/dt, c0,                               &
             ns1, f_snoice)
           
         call define_hist_field(n_dsnow,"dsnow","cm/day",tstr2D, tcstr, &
             "snow formation",                                      &
             "none", mps_to_cmpdy/dt, c0,                               &
             ns1, f_dsnow)
      
         call define_hist_field(n_meltt,"meltt","cm/day",tstr2D, tcstr, &
             "top ice melt",                                          &
             "none", mps_to_cmpdy/dt, c0,                             &
             ns1, f_meltt)
      
         call define_hist_field(n_melts,"melts","cm/day",tstr2D, tcstr, &
             "top snow melt",                                          &
             "none", mps_to_cmpdy/dt, c0,                             &
             ns1, f_melts)
      
         call define_hist_field(n_meltb,"meltb","cm/day",tstr2D, tcstr, &
             "basal ice melt",                                        &
             "none", mps_to_cmpdy/dt, c0,                             &
             ns1, f_meltb)
      
         call define_hist_field(n_meltl,"meltl","cm/day",tstr2D, tcstr, &
             "lateral ice melt",                                      &
             "none", mps_to_cmpdy/dt, c0,                             &
             ns1, f_meltl)
      
         call define_hist_field(n_fresh,"fresh","cm/day",tstr2D, tcstr,   &
             "freshwtr flx ice to ocn (cpl)",                           &
             "if positive, ocean gains fresh water",                    &
             mps_to_cmpdy/rhofresh, c0,                                 &
             ns1, f_fresh)
      
         call define_hist_field(n_fresh_ai,"fresh_ai","cm/day",tstr2D, tcstr, &
             "freshwtr flx ice to ocn",                                     &
             "weighted by ice area", mps_to_cmpdy/rhofresh, c0,             &
             ns1, f_fresh_ai)
      
         call define_hist_field(n_fsalt,"fsalt","kg/m^2/s",tstr2D, tcstr, &
             "salt flux ice to ocn (cpl)",                              &
             "if positive, ocean gains salt", c1, c0,                   &
             ns1, f_fsalt)
      
         call define_hist_field(n_fsalt_ai,"fsalt_ai","kg/m^2/s",tstr2D, tcstr, &
             "salt flux ice to ocean",                                        &
             "weighted by ice area", c1, c0,                                  &
             ns1, f_fsalt_ai)
      
         call define_hist_field(n_fbot,"fbot","W/m^2",tstr2D, tcstr, &
             "heat flux ice to ocean (fbot)",                        &
             "if positive, ocean gains heat", c1, c0,                &
             ns1, f_fbot)
      
         call define_hist_field(n_fhocn,"fhocn","W/m^2",tstr2D, tcstr, &
             "heat flux ice to ocn (cpl)",                           &
             "if positive, ocean gains heat", c1, c0,                &
             ns1, f_fhocn)
      
         call define_hist_field(n_fhocn_ai,"fhocn_ai","W/m^2",tstr2D, tcstr, &
             "heat flux ice to ocean (fhocn_ai)",                          &
             "weighted by ice area", c1, c0,                               &
             ns1, f_fhocn_ai)
      
         call define_hist_field(n_fswthru,"fswthru","W/m^2",tstr2D, tcstr, &
             "SW thru ice to ocean (cpl)",                               &
             "if positive, ocean gains heat", c1, c0,                    &
             ns1, f_fswthru)
      
         call define_hist_field(n_fswthru_ai,"fswthru_ai","W/m^2",tstr2D, tcstr,&
             "SW flux thru ice to ocean",                                     &
             "weighted by ice area", c1, c0,                                  &
             ns1, f_fswthru_ai)
      
         call define_hist_field(n_strairx,"strairx","N/m^2",ustr2D, ucstr, &
             "atm/ice stress (x)",                                       &
             "positive is x direction on U grid", c1, c0,                &
             ns1, f_strairx)
      
         call define_hist_field(n_strairy,"strairy","N/m^2",ustr2D, ucstr, &
             "atm/ice stress (y)",                                       &
             "positive is y direction on U grid", c1, c0,                &
             ns1, f_strairy)
      
         call define_hist_field(n_strtltx,"strtltx","N/m^2",ustr2D, ucstr, &
             "sea sfc tilt stress (x)",                                  &
             "none", c1, c0,                                             &
             ns1, f_strtltx)
      
         call define_hist_field(n_strtlty,"strtlty","N/m^2",ustr2D, ucstr, &
             "sea sfc tilt stress (y)",                                  &
             "none", c1, c0,                                             &
             ns1, f_strtlty)
      
         call define_hist_field(n_strcorx,"strcorx","N/m^2",ustr2D, ucstr, &
             "coriolis stress (x)",                                      &
             "positive is x direction on U grid", c1, c0,                &
             ns1, f_strcorx)
      
         call define_hist_field(n_strcory,"strcory","N/m^2",ustr2D, ucstr, &
             "coriolis stress (y)",                                      &
             "positive is y direction on U grid", c1, c0,                &
             ns1, f_strcory)
      
         call define_hist_field(n_strocnx,"strocnx","N/m^2",ustr2D, ucstr, &
             "ocean/ice stress (x)",                                     &
             "positive is x direction on U grid", c1, c0,                &
             ns1, f_strocnx)
      
         call define_hist_field(n_strocny,"strocny","N/m^2",ustr2D, ucstr, &
             "ocean/ice stress (y)",                                     &
             "positive is y direction on U grid", c1, c0,                &
             ns1, f_strocny)
      
         call define_hist_field(n_strintx,"strintx","N/m^2",ustr2D, ucstr, &
             "internal ice stress (x)",                                  &
             "positive is x direction on U grid", c1, c0,                &
             ns1, f_strintx)
      
         call define_hist_field(n_strinty,"strinty","N/m^2",ustr2D, ucstr, &
             "internal ice stress (y)",                                  &
             "positive is y direction on U grid", c1, c0,                &
             ns1, f_strinty)

         call define_hist_field(n_taubx,"taubx","N/m^2",ustr2D, ucstr,   &
             "seabed (basal) stress (x)",                                &
             "positive is x direction on U grid", c1, c0,                &
             ns1, f_taubx)

         call define_hist_field(n_tauby,"tauby","N/m^2",ustr2D, ucstr,   &
             "seabed (basal) stress (y)",                                &
             "positive is y direction on U grid", c1, c0,                &
             ns1, f_tauby)
      
         call define_hist_field(n_strength,"strength","N/m",tstr2D, tcstr, &
             "compressive ice strength",                                 &
             "none", c1, c0,                                             &
             ns1, f_strength)
      
         call define_hist_field(n_divu,"divu","%/day",tstr2D, tcstr, &
             "strain rate (divergence)",                           &
             "none", secday*c100, c0,                              &
             ns1, f_divu)
      
         call define_hist_field(n_shear,"shear","%/day",tstr2D, tcstr, &
             "strain rate (shear)",                                  &
             "none", secday*c100, c0,                                &
             ns1, f_shear)
      
         call define_hist_field(n_sig1,"sig1","1",ustr2D, ucstr, &
             "norm. principal stress 1",                       &
             "sig1 is instantaneous", c1, c0,                  &
             ns1, f_sig1)
      
         call define_hist_field(n_sig2,"sig2","1",ustr2D, ucstr, &
             "norm. principal stress 2",                       &
             "sig2 is instantaneous", c1, c0,                  &
             ns1, f_sig2)
             
         call define_hist_field(n_sigP,"sigP","1",ustr2D, ucstr, &
             "ice pressure",                       &
             "sigP is instantaneous", c1, c0,                  &
             ns1, f_sigP)
      
         call define_hist_field(n_dvidtt,"dvidtt","cm/day",tstr2D, tcstr, &
             "volume tendency thermo",                                  &
             "none", mps_to_cmpdy, c0,                                  &
             ns1, f_dvidtt)
      
         call define_hist_field(n_dvidtd,"dvidtd","cm/day",tstr2D, tcstr, &
             "volume tendency dynamics",                                &
             "none", mps_to_cmpdy, c0,                                  &
             ns1, f_dvidtd)
      
         call define_hist_field(n_daidtt,"daidtt","%/day",tstr2D, tcstr, &
             "area tendency thermo",                                   &
             "none", secday*c100, c0,                                  &
             ns1, f_daidtt)
      
         call define_hist_field(n_daidtd,"daidtd","%/day",tstr2D, tcstr, &
             "area tendency dynamics",                                 &
             "none", secday*c100, c0,                                  &
             ns1, f_daidtd)
      
         call define_hist_field(n_dagedtt,"dagedtt","day/day",tstr2D, tcstr, &
             "age tendency thermo",                                   &
             "excludes time step increment", c1, c0,                  &
             ns1, f_dagedtt)
      
         call define_hist_field(n_dagedtd,"dagedtd","day/day",tstr2D, tcstr, &
             "age tendency dynamics",                                 &
             "excludes time step increment", c1, c0,                  &
             ns1, f_dagedtd)

         call define_hist_field(n_mlt_onset,"mlt_onset","day of year", &
             tstr2D, tcstr,"melt onset date",                            &
             "midyear restart gives erroneous dates", c1, c0,          &
             ns1, f_mlt_onset)

         call define_hist_field(n_frz_onset,"frz_onset","day of year", &
             tstr2D, tcstr,"freeze onset date",                          &
             "midyear restart gives erroneous dates", c1, c0,          &
             ns1, f_frz_onset)

         call define_hist_field(n_hisnap,"hisnap","m",tstr2D, tcstr, &
             "ice volume snapshot",                                &
             "none", c1, c0,                              &
             ns1, f_hisnap)
      
         call define_hist_field(n_aisnap,"aisnap","1",tstr2D, tcstr, &
             "ice area snapshot",                                  &
             "none", c1, c0,                              &
             ns1, f_aisnap)
      
         call define_hist_field(n_trsig,"trsig","N/m",tstr2D, tcstr, &
             "internal stress tensor trace",                         &
             "ice strength approximation", c1, c0,                   &
             ns1, f_trsig)
      
         call define_hist_field(n_icepresent,"ice_present","1",tstr2D, tcstr, &
             "fraction of time-avg interval that ice is present",           &
             "ice extent flag", c1, c0,                                     &
             ns1, f_icepresent)
      
         call define_hist_field(n_fsurf_ai,"fsurf_ai","W/m^2",tstr2D, tcstr, &
             "net surface heat flux",                                      &
             "positive downward, excludes conductive flux, weighted by ice area", &
             c1, c0, &
             ns1, f_fsurf_ai)

         call define_hist_field(n_fcondtop_ai,"fcondtop_ai","W/m^2", &
             tstr2D, tcstr,"top surface conductive heat flux",         &
             "positive downward, weighted by ice area", c1, c0,      &
             ns1, f_fcondtop_ai)

         call define_hist_field(n_fmeltt_ai,"fmeltt_ai","W/m^2",tstr2D, tcstr, &
             "net surface heat flux causing melt",                           &
             "always >= 0, weighted by ice area", c1, c0,                    &
             ns1, f_fmeltt_ai)

         call define_hist_field(n_a11,"a11"," ",tstr2D, tcstr, &
            "a11: component a11 of the structure tensor",      &
            "none", c1, c0,            &
            ns1, f_a11)

         call define_hist_field(n_a12,"a12"," ",tstr2D, tcstr, &
            "a12: component a12 of the structure tensor",      &
            "none", c1, c0,            &
            ns1, f_a12)

         call define_hist_field(n_e11,"e11","1/s",tstr2D, tcstr, &
            "e11: component e11 of the strain rate tensor",      &
            "none", c1, c0,            &
            ns1, f_e11)

         call define_hist_field(n_e12,"e12","1/s",tstr2D, tcstr, &
            "e12: component e12 of the strain rate tensor",      &
            "none", c1, c0,            &
            ns1, f_e12)

         call define_hist_field(n_e22,"e22","1/s",tstr2D, tcstr, &
            "e22: component e22 of the strain rate tensor",      &
            "none", c1, c0,            &
            ns1, f_e22)

         call define_hist_field(n_s11,"s11","kg/s^2",tstr2D, tcstr, &
            "s11: component s11 of the stress tensor",              &
            "none", c1, c0,            &
            ns1, f_s11)

         call define_hist_field(n_s12,"s12","kg/s^2",tstr2D, tcstr, &
            "s12: component s12 of the stress tensor",              &
            "none", c1, c0,            &
            ns1, f_s12)

         call define_hist_field(n_s22,"s22","kg/s^2",tstr2D, tcstr, &
            "s22: component s12 of the stress tensor",              &
            "none", c1, c0,            &
            ns1, f_s22)

         call define_hist_field(n_yieldstress11,"yieldstress11","kg/s^2",tstr2D, tcstr, &
            "yieldstress11: component 11 of the yieldstress tensor",                    &
            "none", c1, c0,            &
            ns1, f_yieldstress11)

         call define_hist_field(n_yieldstress12,"yieldstress12","kg/s^2",tstr2D, tcstr, &
            "yieldstress12: component 12 of the yieldstress tensor",                    &
            "none", c1, c0,            &
            ns1, f_yieldstress12)

         call define_hist_field(n_yieldstress22,"yieldstress22","kg/s^2",tstr2D, tcstr, &
            "yieldstress22: component 12 of the yieldstress tensor",                    &
            "none", c1, c0,            &
            ns1, f_yieldstress22)

      ! Tracers

      ! Ice Age
         call define_hist_field(n_iage,"iage","years",tstr2D, tcstr, &
             "sea ice age",                                        &
             "none", c1/(secday*days_per_year), c0,                &
             ns1, f_iage)

      ! First Year Ice Area
         call define_hist_field(n_FY,"FYarea"," ",tstr2D, tcstr, &
             "first-year ice area",                            &
             "weighted by ice area", c1, c0,                   &
              ns1, f_FY)

      ! CMIP 2D variables

         call define_hist_field(n_sithick,"sithick","m",tstr2D, tcstr, &
             "sea ice thickness",                             &
             "volume divided by area", c1, c0, &
             ns1, f_sithick)
      
         call define_hist_field(n_siage,"siage","s",tstr2D, tcstr,    &
             "sea ice age",                             &
             "none", c1, c0,                                      &
             ns1, f_siage)
      
         call define_hist_field(n_sisnthick,"sisnthick","m",tstr2D, tcstr,    &
             "sea ice snow thickness",                             &
             "snow volume divided by area", c1, c0, &
             ns1, f_sisnthick)
      
         call define_hist_field(n_sitemptop,"sitemptop","K",tstr2D, tcstr,    &
             "sea ice surface temperature", &
             "none", c1, c0,           &
             ns1, f_sitemptop)
      
         call define_hist_field(n_sitempsnic,"sitempsnic","K",tstr2D, tcstr,    &
             "snow ice interface temperature", &
             "surface temperature when no snow present", c1, c0, &
             ns1, f_sitempsnic)
      
         call define_hist_field(n_sitempbot,"sitempbot","K",tstr2D, tcstr,    &
             "sea ice bottom temperature",                             &
             "none", c1, c0,           &
             ns1, f_sitempbot)

         call define_hist_field(n_siu,"siu","m/s",ustr2D, ucstr,  &
             "ice x velocity component", &
             "none", c1, c0,         &
             ns1, f_siu)

         call define_hist_field(n_siv,"siv","m/s",ustr2D, ucstr,  &
             "ice y velocity component", &
             "none", c1, c0,         &
             ns1, f_siv)
      
         call define_hist_field(n_sidmasstranx,"sidmasstranx","kg/s",ustr2D, ucstr,  &
             "x component of snow and sea ice mass transport", &
             "none", c1, c0,         &
             ns1, f_sidmasstranx)
      
         call define_hist_field(n_sidmasstrany,"sidmasstrany","kg/s",ustr2D, ucstr,  &
             "y component of snow and sea ice mass transport", &
             "none", c1, c0,         &
             ns1, f_sidmasstrany)
      
         call define_hist_field(n_sistrxdtop,"sistrxdtop","N m-2",ustr2D, ucstr,  &
             "x component of atmospheric stress on sea ice", &
             "none", c1, c0,         &
             ns1, f_sistrxdtop)
      
         call define_hist_field(n_sistrydtop,"sistrydtop","N m-2",ustr2D, ucstr,  &
             "y component of atmospheric stress on sea ice", &
             "none", c1, c0,         &
             ns1, f_sistrydtop)
      
         call define_hist_field(n_sistrxubot,"sistrxubot","N m-2",ustr2D, ucstr,  &
             "x component of ocean stress on sea ice", &
             "none", c1, c0,         &
             ns1, f_sistrxubot)
      
         call define_hist_field(n_sistryubot,"sistryubot","N m-2",ustr2D, ucstr,  &
             "y component of ocean stress on sea ice", &
             "none", c1, c0,         &
             ns1, f_sistryubot)
      
         call define_hist_field(n_sicompstren,"sicompstren","N m-1",tstr2D, tcstr,  &
             "compressive sea ice strength",                      &
             "none", c1, c0,         &
             ns1, f_sicompstren)

         call define_hist_field(n_sispeed,"sispeed","m/s",ustr2D, ucstr, &
             "ice speed",                                  &
             "none", c1, c0,         &
             ns1, f_sispeed)

         call define_hist_field(n_sidir,"sidir","deg",ustr2D, ucstr,  &
             "ice direction",                                         &
             "vector direction - going to", c1, c0,                   &
             ns1, f_sidir)
      
         call define_hist_field(n_sialb,"sialb","1",tstr2D, tcstr,  &
             "sea ice albedo",                                  &
             "none", c1, c0,         &
             ns1, f_sialb)
      
         call define_hist_field(n_sihc,"sihc","J m-2",tstr2D, tcstr,  &
             "sea ice heat content",                                  &
             "none", c1, c0,         &
             ns1, f_sihc)
      
         call define_hist_field(n_sisnhc,"sisnhc","J m-2",tstr2D, tcstr, &
             "snow heat content",                                  &
             "none", c1, c0,         &
             ns1, f_sisnhc)
      
         call define_hist_field(n_sidconcth,"sidconcth","1/s",tstr2D, tcstr,  &
             "sea ice area change from thermodynamics",              &
             "none", c1, c0,         &
             ns1, f_sidconcth)
      
         call define_hist_field(n_sidconcdyn,"sidconcdyn","1/s",tstr2D, tcstr,  &
             "sea ice area change from dynamics",                      &
             "none", c1, c0,         &
             ns1, f_sidconcdyn)
      
         call define_hist_field(n_sidmassth,"sidmassth","kg m-2 s-1",tstr2D, tcstr,  &
             "sea ice mass change from thermodynamics",              &
             "none", c1, c0,         &
             ns1, f_sidmassth)
      
         call define_hist_field(n_sidmassdyn,"sidmassdyn","kg m-2 s-1",tstr2D, tcstr,  &
             "sea ice mass change from dynamics",                      &
             "none", c1, c0,         &
             ns1, f_sidmassdyn)

         call define_hist_field(n_sidmassgrowthwat,"sidmassgrowthwat","kg m-2 s-1",tstr2D, tcstr,  &
             "sea ice mass change from frazil",                      &
             "none", c1, c0,         &
             ns1, f_sidmassgrowthwat)
      
         call define_hist_field(n_sidmassgrowthbot,"sidmassgrowthbot","kg m-2 s-1",tstr2D, tcstr,  &
             "sea ice mass change from basal growth", &
             "none", c1, c0,         &
             ns1, f_sidmassgrowthbot)
      
         call define_hist_field(n_sidmasssi,"sidmasssi","kg m-2 s-1",tstr2D, tcstr,  &
             "sea ice mass change from snow-ice formation", &
             "none", c1, c0,         &
             ns1, f_sidmasssi)
      
         call define_hist_field(n_sidmassevapsubl,"sidmassevapsubl","kg m-2 s-1",tstr2D, tcstr,  &
             "sea ice mass change from evaporation and sublimation", &
             "none", c1, c0,         &
             ns1, f_sidmassevapsubl)
      
         call define_hist_field(n_sndmasssubl,"sndmassubl","kg m-2 s-1",tstr2D, tcstr,  &
             "snow mass change from evaporation and sublimation", &
             "none", c1, c0,         &
             ns1, f_sndmasssubl)
      
         call define_hist_field(n_sidmassmelttop,"sidmassmelttop","kg m-2 s-1",tstr2D, tcstr,  &
             "sea ice mass change top melt",                      &
             "none", c1, c0,         &
             ns1, f_sidmassmelttop)
      
         call define_hist_field(n_sidmassmeltbot,"sidmassmeltbot","kg m-2 s-1",tstr2D, tcstr,  &
             "sea ice mass change bottom melt",                      &
             "none", c1, c0,         &
             ns1, f_sidmassmeltbot)
      
         call define_hist_field(n_sidmasslat,"sidmasslat","kg m-2 s-1",tstr2D, tcstr,  &
             "sea ice mass change lateral melt",                      &
             "none", c1, c0,         &
             ns1, f_sidmasslat)

        call define_hist_field(n_sndmasssnf,"sndmasssnf","kg m-2 s-1",tstr2D, tcstr,  &
             "snow mass change from snow fall",                      &
             "none", c1, c0,         &
             ns1, f_sndmasssnf)
      
         call define_hist_field(n_sndmassmelt,"sndmassmelt","kg m-2 s-1",tstr2D, tcstr,  &
             "snow mass change from snow melt",                      &
             "none", c1, c0,         &
             ns1, f_sndmassmelt)
      
         call define_hist_field(n_siflswdtop,"siflswdtop","W/m2",tstr2D, tcstr, &
             "down shortwave flux over sea ice", &
             "positive downward", c1, c0,                            &
             ns1, f_siflswdtop)
      
         call define_hist_field(n_siflswutop,"siflswutop","W/m2",tstr2D, tcstr, &
             "upward shortwave flux over sea ice", &
             "positive downward", c1, c0,                            &
             ns1, f_siflswutop)
      
         call define_hist_field(n_siflswdbot,"siflswdbot","W/m2",tstr2D, tcstr, &
             "down shortwave flux at bottom of ice", &
             "positive downward", c1, c0,                            &
             ns1, f_siflswdbot)
      
         call define_hist_field(n_sifllwdtop,"sifllwdtop","W/m2",tstr2D, tcstr, &
             "down longwave flux over sea ice", &
             "positive downward", c1, c0,                            &
             ns1, f_sifllwdtop)
      
         call define_hist_field(n_sifllwutop,"sifllwutop","W/m2",tstr2D, tcstr, &
             "upward longwave flux over sea ice", &
             "positive downward", c1, c0,                            &
             ns1, f_sifllwutop)
      
         call define_hist_field(n_siflsenstop,"siflsenstop","W/m2",tstr2D, tcstr, &
             "sensible heat flux over sea ice", &
             "positive downward", c1, c0,                            &
             ns1, f_siflsenstop)

         call define_hist_field(n_siflsensupbot,"siflsensupbot","W/m2",tstr2D, tcstr, &
             "sensible heat flux at bottom of sea ice", &
             "positive downward", c1, c0,                            &
             ns1, f_siflsensupbot)
      
         call define_hist_field(n_sifllatstop,"sifllatstop","W/m2",tstr2D, tcstr, &
             "latent heat flux over sea ice", &
             "positive downward", c1, c0,                            &
             ns1, f_sifllatstop)
      
         call define_hist_field(n_siflcondtop,"siflcondtop","W/m2",tstr2D, tcstr, &
             "conductive heat flux at top of sea ice", &
             "positive downward", c1, c0,                            &
             ns1, f_siflcondtop)
      
         call define_hist_field(n_siflcondbot,"siflcondbot","W/m2",tstr2D, tcstr, &
             "conductive heat flux at bottom of sea ice", &
             "positive downward", c1, c0,                            &
             ns1, f_siflcondbot)
      
         call define_hist_field(n_sipr,"sipr","kg m-2 s-1",tstr2D, tcstr, &
             "rainfall over sea ice", &
             "none", c1, c0,                            &
             ns1, f_sipr)
      
         call define_hist_field(n_sifb,"sifb","m",tstr2D, tcstr, &
             "sea ice freeboard above sea level",                &
             "none", c1, c0,                            &
             ns1, f_sifb)
      
         call define_hist_field(n_siflsaltbot,"siflsaltbot","kg m-2 s-1",tstr2D, tcstr, &
             "salt flux from sea ice", &
             "positive downward", c1, c0,                            &
             ns1, f_siflsaltbot)
      
         call define_hist_field(n_siflfwbot,"siflfwbot","kg m-2 s-1",tstr2D, tcstr, &
             "fresh water flux from sea ice", &
             "positive downward", c1, c0,                            &
             ns1, f_siflfwbot)

         call define_hist_field(n_siflfwdrain,"siflfwdrain","kg m-2 s-1",tstr2D, tcstr, &
             "fresh water drainage through sea ice", &
             "positive downward", c1, c0,                            &
             ns1, f_siflfwdrain)
      
         call define_hist_field(n_sidragtop,"sidragtop","1",tstr2D, tcstr, &
             "atmospheric drag over sea ice", &
             "none", c1, c0,                            &
             ns1, f_sidragtop)
      
         call define_hist_field(n_sirdgthick,"sirdgthick","m",tstr2D, tcstr, &
             "sea ice ridge thickness", &
             "vrdg divided by ardg", c1, c0, &
             ns1, f_sirdgthick)
      
         call define_hist_field(n_siforcetiltx,"siforcetiltx","N m-2",tstr2D, tcstr, &
             "sea surface tilt term", &
             "none", c1, c0,                            &
             ns1, f_siforcetiltx)
      
         call define_hist_field(n_siforcetilty,"siforcetilty","N m-2",tstr2D, tcstr, &
             "sea surface tile term", &
             "none", c1, c0,                            &
             ns1, f_siforcetilty)
      
         call define_hist_field(n_siforcecoriolx,"siforcecoriolx","N m-2",tstr2D, tcstr, &
             "coriolis term",                                      &
             "none", c1, c0,                            &
             ns1, f_siforcecoriolx)
      
         call define_hist_field(n_siforcecorioly,"siforcecorioly","N m-2",tstr2D, tcstr, &
             "coriolis term",                                      &
             "none", c1, c0,                            &
             ns1, f_siforcecorioly)
      
         call define_hist_field(n_siforceintstrx,"siforceintstrx","N m-2",tstr2D, tcstr, &
             "internal stress term", &
             "none", c1, c0,                            &
             ns1, f_siforceintstrx)

         call define_hist_field(n_siforceintstry,"siforceintstry","N m-2",tstr2D, tcstr, &
             "internal stress term", &
             "none", c1, c0,                            &
             ns1, f_siforceintstry)

         call define_hist_field(n_sistreave,"sistreave","N m-1",ustr2D, ucstr, &
             "average normal stress",                       &
             "sistreave is instantaneous", c1, c0,                  &
             ns1, f_sistreave)
      
         call define_hist_field(n_sistremax,"sistremax","N m-1",ustr2D, ucstr, &
             "maximum shear stress",                       &
             "sistremax is instantaneous", c1, c0,                  &
             ns1, f_sistremax)

      endif ! if (histfreq(ns1) /= 'x') then
      enddo ! ns1

      ! other 2D history variables

      ! mechanical redistribution
      call init_hist_mechred_2D

      ! melt ponds
      call init_hist_pond_2D

      ! biogeochemistry
      call init_hist_bgc_2D

      ! form drag
      call init_hist_drag_2D

      ! floe size distribution
      call init_hist_fsd_2D

      !-----------------------------------------------------------------
      ! 3D (category) variables looped separately for ordering
      !-----------------------------------------------------------------
      do ns1 = 1, nstreams
      if (histfreq(ns1) /= 'x') then

           call define_hist_field(n_aicen,"aicen","1",tstr3Dc, tcstr, & 
              "ice area, categories","none", c1, c0,                  &            
              ns1, f_aicen)

           call define_hist_field(n_vicen,"vicen","m",tstr3Dc, tcstr, & 
              "ice volume, categories","none", c1, c0,                &            
              ns1, f_vicen)

           call define_hist_field(n_vsnon,"vsnon","m",tstr3Dc, tcstr, &
              "snow depth on ice, categories","volume per unit area of snow", c1, c0, &
              ns1, f_vsnon)

           call define_hist_field(n_snowfracn,"snowfracn","1",tstr3Dc, tcstr, &
             "category mean snow fraction",                     &
             "snow fraction per unit grid cell area", c1, c0,       &
              ns1, f_snowfracn)

           call define_hist_field(n_fsurfn_ai,"fsurfn_ai","W/m^2",tstr3Dc, tcstr, & 
              "net surface heat flux, categories","weighted by ice area", c1, c0, &
              ns1, f_fsurfn_ai)
   
           call define_hist_field(n_fcondtopn_ai,"fcondtopn_ai","W/m^2",tstr3Dc, tcstr, &
              "top sfc conductive heat flux, cat","weighted by ice area", c1, c0,       &
              ns1, f_fcondtopn_ai)

           call define_hist_field(n_fmelttn_ai,"fmelttn_ai","W/m^2",tstr3Dc, tcstr, & 
              "net sfc heat flux causing melt, cat","weighted by ice area", c1, c0, &            
              ns1, f_fmelttn_ai)

           call define_hist_field(n_flatn_ai,"flatn_ai","W/m^2",tstr3Dc, tcstr, & 
              "latent heat flux, category","weighted by ice area", c1, c0,      &            
              ns1, f_flatn_ai)

           call define_hist_field(n_fsensn_ai,"fsensn_ai","W/m^2",tstr3Dc, tcstr, & 
              "sensible heat flux, category","weighted by ice area", c1, c0,      &            
              ns1, f_fsensn_ai)

           call define_hist_field(n_keffn_top,"keffn_top","W/m^2/K",tstr3Dc, tcstr, &
              "effective thermal conductivity of the top ice layer, categories", &
              "multilayer scheme", c1, c0,      &           
              ns1, f_keffn_top)

           ! CMIP 3D
           call define_hist_field(n_siitdconc,"siitdconc","1",tstr3Dc, tcstr, &
              "ice area, categories","none", c1, c0,                  &
              ns1, f_siitdconc)

           call define_hist_field(n_siitdthick,"siitdthick","m",tstr3Dc, tcstr, &
              "ice thickness, categories","none", c1, c0, &
              ns1, f_siitdthick)

           call define_hist_field(n_siitdsnthick,"siitdsnthick","m",tstr3Dc, tcstr, &
              "snow thickness, categories","none", c1, c0, &
              ns1, f_siitdsnthick)

      endif ! if (histfreq(ns1) /= 'x') then
      enddo ! ns1

      ! other 3D (category) history variables

      ! mechanical redistribution
      call init_hist_mechred_3Dc

      ! melt ponds
      call init_hist_pond_3Dc

      ! biogeochemistry
      call init_hist_bgc_3Dc

      !-----------------------------------------------------------------
      ! 3D (vertical) variables must be looped separately
      !-----------------------------------------------------------------

!      do ns1 = 1, nstreams
!      if (histfreq(ns1) /= 'x') then

!         call define_hist_field(n_field3dz,"field3dz","1",tstr3Dz, tcstr, & 
!            "example 3dz field",                    &
!            "vertical profile", c1, c0,                  &
!            ns1, f_field3dz)

!      endif ! if (histfreq(ns1) /= 'x') then
!      enddo ! ns1 

      ! biogeochemistry
      call init_hist_bgc_3Db  
      call init_hist_bgc_3Da

      !-----------------------------------------------------------------
      ! 3D (floe size) variables must be looped separately
      !-----------------------------------------------------------------

      ! floe size distribution
      call init_hist_fsd_3Df

      !-----------------------------------------------------------------
      ! 4D (categories, vertical) variables must be looped separately
      !-----------------------------------------------------------------

      do ns1 = 1, nstreams
      if (histfreq(ns1) /= 'x') then

         call define_hist_field(n_Tinz,"Tinz","C",tstr4Di, tcstr, & 
            "ice internal temperatures on CICE grid",          &
            "vertical profile", c1, c0,                    &
            ns1, f_Tinz)

         call define_hist_field(n_Sinz,"Sinz","ppt",tstr4Di, tcstr, & 
            "ice internal bulk salinity",          &
            "vertical profile", c1, c0,                    &
            ns1, f_Sinz)

      endif ! if (histfreq(ns1) /= 'x') then
      enddo ! ns1

      do ns1 = 1, nstreams
      if (histfreq(ns1) /= 'x') then

         call define_hist_field(n_Tsnz,"Tsnz","C",tstr4Ds, tcstr, & 
            "snow internal temperatures",          &
            "vertical profile", c1, c0,                    &
            ns1, f_Tsnz)

      endif ! if (histfreq(ns1) /= 'x') then
      enddo

       if (f_Tinz   (1:1) /= 'x') then
            if (allocated(Tinz4d)) deallocate(Tinz4d)
            allocate(Tinz4d(nx_block,ny_block,nzilyr,ncat_hist))
       endif
       if (f_Sinz   (1:1) /= 'x')  then
            if (allocated(Sinz4d)) deallocate(Sinz4d)
            allocate(Sinz4d(nx_block,ny_block,nzilyr,ncat_hist))
       endif
       if (f_Tsnz   (1:1) /= 'x') then
            if (allocated(Tsnz4d)) deallocate(Tsnz4d)
            allocate(Tsnz4d(nx_block,ny_block,nzslyr,ncat_hist))
       endif

      !-----------------------------------------------------------------
      ! 4D (floe size, thickness categories) variables looped separately
      !-----------------------------------------------------------------

      ! floe size distribution
       call init_hist_fsd_4Df

      !-----------------------------------------------------------------
      ! fill igrd array with namelist values
      !-----------------------------------------------------------------

      igrd=.true.

      igrd(n_tmask     ) = f_tmask
      igrd(n_blkmask   ) = f_blkmask
      igrd(n_tarea     ) = f_tarea
      igrd(n_uarea     ) = f_uarea
      igrd(n_dxt       ) = f_dxt
      igrd(n_dyt       ) = f_dyt
      igrd(n_dxu       ) = f_dxu
      igrd(n_dyu       ) = f_dyu
      igrd(n_HTN       ) = f_HTN
      igrd(n_HTE       ) = f_HTE
      igrd(n_ANGLE     ) = f_ANGLE
      igrd(n_ANGLET    ) = f_ANGLET

      igrdz=.true.
      igrdz(n_NCAT     ) = f_NCAT
      igrdz(n_VGRDi    ) = f_VGRDi
      igrdz(n_VGRDs    ) = f_VGRDs
      igrdz(n_VGRDb    ) = f_VGRDb
      igrdz(n_VGRDa    ) = f_VGRDa
      igrdz(n_NFSD     ) = f_NFSD

      !-----------------------------------------------------------------
      ! diagnostic output
      !-----------------------------------------------------------------

      ntmp(:) = 0
      if (my_task == master_task) then
        write(nu_diag,*) ' '
        write(nu_diag,*) 'total number of history fields = ',num_avail_hist_fields_tot
        write(nu_diag,*) 'max number of history fields   = ',max_avail_hist_fields
        write(nu_diag,*) 'The following variables will be ', &
                         'written to the history tape: '
        write(nu_diag,101) 'description','units','variable','frequency','x'
        do n=1,num_avail_hist_fields_tot
           if (avail_hist_fields(n)%vhistfreq_n /= 0) &
           write(nu_diag,100) avail_hist_fields(n)%vdesc, &
              avail_hist_fields(n)%vunit, avail_hist_fields(n)%vname, &
              avail_hist_fields(n)%vhistfreq,avail_hist_fields(n)%vhistfreq_n
           do ns = 1, nstreams
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) &
                 ntmp(ns)=ntmp(ns)+1
           enddo
        enddo ! num_avail_hist_fields_tot
        write(nu_diag,*) ' '
      endif
  100 format (1x,a40,2x,a16,2x,a12,1x,a1,2x,i6)
  101 format (2x,a19,10x,a16,9x,a12,2x,a,3x,a1)

      call broadcast_array(ntmp, master_task)
      do ns = 1, nstreams
         if (ntmp(ns)==0) histfreq_n(ns) = 0
      enddo

      !-----------------------------------------------------------------
      ! initialize the history arrays
      !-----------------------------------------------------------------

      if (allocated(a2D)) deallocate(a2D)
      if (num_avail_hist_fields_2D > 0) &
      allocate(a2D(nx_block,ny_block,num_avail_hist_fields_2D,max_blocks))

      if (allocated(a3Dc)) deallocate(a3Dc)
      if (num_avail_hist_fields_3Dc > 0) &
      allocate(a3Dc(nx_block,ny_block,ncat_hist,num_avail_hist_fields_3Dc,max_blocks))

      if (allocated(a3Dz)) deallocate(a3Dz)
      if (num_avail_hist_fields_3Dz > 0) &
      allocate(a3Dz(nx_block,ny_block,nzilyr,num_avail_hist_fields_3Dz,max_blocks))

      if (allocated(a3Db)) deallocate(a3Db)
      if (num_avail_hist_fields_3Db > 0) &
      allocate(a3Db(nx_block,ny_block,nzblyr,num_avail_hist_fields_3Db,max_blocks))

      if (allocated(a3Da)) deallocate(a3Da)
      if (num_avail_hist_fields_3Da > 0) &
      allocate(a3Da(nx_block,ny_block,nzalyr,num_avail_hist_fields_3Da,max_blocks))

      if (allocated(a3Df)) deallocate(a3Df)
      if (num_avail_hist_fields_3Df > 0) &
      allocate(a3Df(nx_block,ny_block,nfsd_hist,num_avail_hist_fields_3Df,max_blocks))

      if (allocated(a4Di)) deallocate(a4Di)
      if (num_avail_hist_fields_4Di > 0) &
      allocate(a4Di(nx_block,ny_block,nzilyr,ncat_hist,num_avail_hist_fields_4Di,max_blocks))

      if (allocated(a4Ds)) deallocate(a4Ds)
      if (num_avail_hist_fields_4Ds > 0) &
      allocate(a4Ds(nx_block,ny_block,nzslyr,ncat_hist,num_avail_hist_fields_4Ds,max_blocks))

      if (allocated(a4Df)) deallocate(a4Df)
      if (num_avail_hist_fields_4Df > 0) &
      allocate(a4Df(nx_block,ny_block,nfsd_hist,ncat_hist,num_avail_hist_fields_4Df,max_blocks))

      if (allocated(a2D))  a2D (:,:,:,:)     = c0
      if (allocated(a3Dc)) a3Dc(:,:,:,:,:)   = c0
      if (allocated(a3Dz)) a3Dz(:,:,:,:,:)   = c0
      if (allocated(a3Db)) a3Db(:,:,:,:,:)   = c0
      if (allocated(a3Da)) a3Da(:,:,:,:,:)   = c0
      if (allocated(a3Df)) a3Df(:,:,:,:,:)   = c0
      if (allocated(a4Di)) a4Di(:,:,:,:,:,:) = c0
      if (allocated(a4Ds)) a4Ds(:,:,:,:,:,:) = c0
      if (allocated(a4Df)) a4Df(:,:,:,:,:,:) = c0
      avgct(:) = c0
      albcnt(:,:,:,:) = c0

      if (restart .and. yday >= c2) then
! restarting midyear gives erroneous onset dates
         mlt_onset = 999._dbl_kind 
         frz_onset = 999._dbl_kind 
      else
         mlt_onset = c0
         frz_onset = c0
      endif

      end subroutine init_hist

!=======================================================================

! accumulate average ice quantities or snapshots
!
! author:   Elizabeth C. Hunke, LANL

      subroutine accum_hist (dt)

      use ice_blocks, only: block, get_block, nx_block, ny_block
      use ice_domain, only: blocks_ice, nblocks
      use ice_domain_size, only: nfsd
      use ice_grid, only: tmask, lmask_n, lmask_s, dxu, dyu
      use ice_calendar, only: new_year, write_history, &
                              write_ic, time, histfreq, nstreams, month, &
                              new_month
      use ice_dyn_eap, only: a11, a12, e11, e12, e22, s11, s12, s22, &
          yieldstress11, yieldstress12, yieldstress22
      use ice_dyn_shared, only: kdyn, principal_stress
      use ice_flux, only: fsw, flw, fsnow, frain, sst, sss, uocn, vocn, &
          frzmlt_init, scale_factor, fswabs, fswthru, alvdr, alvdf, alidr, alidf, &
          albice, albsno, albpnd, coszen, flat, fsens, flwout, evap, evaps, evapi, &
          Tair, Tref, Qref, congel, frazil, frazil_diag, snoice, dsnow, &
          melts, meltb, meltt, meltl, fresh, fsalt, fresh_ai, fsalt_ai, &
          fhocn, fhocn_ai, uatm, vatm, fbot, Tbot, Tsnice, &
          fswthru_ai, strairx, strairy, strtltx, strtlty, strintx, strinty, &
          taubx, tauby, strocnx, strocny, fm, daidtt, dvidtt, daidtd, dvidtd, fsurf, &
          fcondtop, fcondbot, fsurfn, fcondtopn, flatn, fsensn, albcnt, &
          stressp_1, stressm_1, stress12_1, &
          stressp_2, &
          stressp_3, &
          stressp_4, sig1, sig2, sigP, &
          mlt_onset, frz_onset, dagedtt, dagedtd, fswint_ai, keffn_top, &
          snowfrac, alvdr_ai, alvdf_ai, alidr_ai, alidf_ai, update_ocn_f
      use ice_arrays_column, only: snowfracn, Cdn_atm
      use ice_history_shared ! almost everything
      use ice_history_write, only: ice_write_hist
      use ice_history_bgc, only: accum_hist_bgc
      use ice_history_mechred, only: accum_hist_mechred
      use ice_history_pond, only: accum_hist_pond
      use ice_history_drag, only: accum_hist_drag
      use icepack_intfc, only: icepack_mushy_density_brine, icepack_mushy_liquid_fraction
      use icepack_intfc, only: icepack_mushy_temperature_mush
      use ice_history_fsd, only: accum_hist_fsd
      use ice_state ! almost everything
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_readwrite

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
           i,j,k,ic,n,ns,nn, &
           iblk             , & ! block index
           ilo,ihi,jlo,jhi  , & ! beginning and end of physical domain
           nstrm                ! nstreams (1 if writing initial condition)

      real (kind=dbl_kind) :: &
           ravgct           , & ! 1/avgct
           ravgctz              ! 1/avgct

      real (kind=dbl_kind) :: & 
           qn                , & ! temporary variable for enthalpy
           sn                    ! temporary variable for salinity

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         worka, workb, ravgip

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat_hist) :: &
         ravgipn, worka3

      real (kind=dbl_kind) :: awtvdr, awtidr, awtvdf, awtidf, puny, secday, rad_to_deg
      real (kind=dbl_kind) :: Tffresh, rhoi, rhos, rhow, ice_ref_salinity
      real (kind=dbl_kind) :: rho_ice, rho_ocn, Tice, Sbr, phi, rhob, dfresh, dfsalt
      logical (kind=log_kind) :: formdrag, skl_bgc
      logical (kind=log_kind) :: tr_pond, tr_aero, tr_brine
      integer (kind=int_kind) :: ktherm
      integer (kind=int_kind) :: nt_sice, nt_qice, nt_qsno, nt_iage, nt_FY, nt_Tsfc, &
                                 nt_alvl, nt_vlvl

      type (block) :: &
         this_block           ! block information for current block
      character(len=*), parameter :: subname = '(accum_hist)'

      call icepack_query_parameters(awtvdr_out=awtvdr, awtidr_out=awtidr, &
           awtvdf_out=awtvdf, awtidf_out=awtidf, puny_out=puny, secday_out=secday, &
           rad_to_deg_out=rad_to_deg)
      call icepack_query_parameters(Tffresh_out=Tffresh, rhoi_out=rhoi, rhos_out=rhos, &
           rhow_out=rhow, ice_ref_salinity_out=ice_ref_salinity)
      call icepack_query_parameters(formdrag_out=formdrag, skl_bgc_out=skl_bgc, ktherm_out=ktherm)
      call icepack_query_tracer_flags(tr_pond_out=tr_pond, tr_aero_out=tr_aero, &
           tr_brine_out=tr_brine)
      call icepack_query_tracer_indices(nt_sice_out=nt_sice, nt_qice_out=nt_qice, &
           nt_qsno_out=nt_qsno, nt_iage_out=nt_iage, nt_FY_out=nt_FY, nt_Tsfc_out=nt_Tsfc, &
           nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !---------------------------------------------------------------
      ! increment step counter
      !---------------------------------------------------------------

      n2D     = num_avail_hist_fields_2D
      n3Dccum = n2D     + num_avail_hist_fields_3Dc
      n3Dzcum = n3Dccum + num_avail_hist_fields_3Dz
      n3Dbcum = n3Dzcum + num_avail_hist_fields_3Db
      n3Dacum = n3Dbcum + num_avail_hist_fields_3Da
      n3Dfcum = n3Dacum + num_avail_hist_fields_3Df
      n4Dicum = n3Dfcum + num_avail_hist_fields_4Di
      n4Dscum = n4Dicum + num_avail_hist_fields_4Ds
      n4Dfcum = n4Dscum + num_avail_hist_fields_4Df ! should equal num_avail_hist_fields_tot

      do ns = 1,nstreams
         if (.not. hist_avg .or. histfreq(ns) == '1') then  ! write snapshots
           do n = 1,n2D
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) &
                  a2D(:,:,n,:) = c0
           enddo
           do n = n2D + 1, n3Dccum                  
              nn = n - n2D
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) &
                  a3Dc(:,:,:,nn,:) = c0
           enddo
           do n = n3Dccum + 1, n3Dzcum
              nn = n - n3Dccum
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) &
                  a3Dz(:,:,:,nn,:) = c0
           enddo
           do n = n3Dzcum + 1, n3Dbcum
              nn = n - n3Dzcum
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) &
                  a3Db(:,:,:,nn,:) = c0
           enddo
           do n = n3Dbcum + 1, n3Dacum
              nn = n - n3Dbcum
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) &
                  a3Da(:,:,:,nn,:) = c0
           enddo
           do n = n3Dacum + 1, n3Dfcum
              nn = n - n3Dacum
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) &
                  a3Df(:,:,:,nn,:) = c0
           enddo
           do n = n3Dfcum + 1, n4Dicum
              nn = n - n3Dfcum
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) &
                  a4Di(:,:,:,:,nn,:) = c0
           enddo
           do n = n4Dicum + 1, n4Dscum
              nn = n - n4Dicum
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) &
                  a4Ds(:,:,:,:,nn,:) = c0
           enddo
           do n = n4Dscum + 1, n4Dfcum
              nn = n - n4Dscum
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) &
                  a4Df(:,:,:,:,nn,:) = c0
           enddo
           avgct(ns) = c1
         else                      ! write averages over time histfreq
           avgct(ns) = avgct(ns) + c1
!           if (avgct(ns) == c1) time_beg(ns) = (time-dt)/int(secday)
           if (avgct(ns) == c1) then
              time_beg(ns) = (time-dt)/int(secday)
              time_beg(ns) = real(time_beg(ns),kind=real_kind)
           endif
         endif
      enddo

      !---------------------------------------------------------------
      ! increment field
      !---------------------------------------------------------------

! MHRI: CHECK THIS OMP ... Maybe ok after "dfresh,dfsalt" added
      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block, &
      !$OMP             k,n,qn,ns,sn,rho_ocn,rho_ice,Tice,Sbr,phi,rhob,dfresh,dfsalt, &
      !$OMP             worka,workb,worka3,Tinz4d,Sinz4d,Tsnz4d)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         workb(:,:) = aice_init(:,:,iblk)

!        if (f_example(1:1) /= 'x') &
!            call accum_hist_field(n_example,iblk, vice(:,:,iblk), a2D)
         if (f_hi     (1:1) /= 'x') &
             call accum_hist_field(n_hi,     iblk, vice(:,:,iblk), a2D)
         if (f_hs     (1:1) /= 'x') &
             call accum_hist_field(n_hs,     iblk, vsno(:,:,iblk), a2D)
         if (f_snowfrac(1:1) /= 'x') &
             call accum_hist_field(n_snowfrac, iblk, snowfrac(:,:,iblk), a2D)
         if (f_Tsfc   (1:1) /= 'x') &
             call accum_hist_field(n_Tsfc,   iblk, trcr(:,:,nt_Tsfc,iblk), a2D)
         if (f_aice   (1:1) /= 'x') &
             call accum_hist_field(n_aice,   iblk, aice(:,:,iblk), a2D)
         if (f_uvel   (1:1) /= 'x') &
             call accum_hist_field(n_uvel,   iblk, uvel(:,:,iblk), a2D)
         if (f_vvel   (1:1) /= 'x') &
             call accum_hist_field(n_vvel,   iblk, vvel(:,:,iblk), a2D)
         if (f_uatm   (1:1) /= 'x') &
             call accum_hist_field(n_uatm,   iblk, uatm(:,:,iblk), a2D)
         if (f_vatm   (1:1) /= 'x') &
             call accum_hist_field(n_vatm,   iblk, vatm(:,:,iblk), a2D)
         if (f_atmspd   (1:1) /= 'x') &
             call accum_hist_field(n_atmspd,   iblk, sqrt( &
                                  (uatm(:,:,iblk)*uatm(:,:,iblk)) + &
                                  (vatm(:,:,iblk)*vatm(:,:,iblk))), a2D)
         if (f_atmdir(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (abs(uatm(i,j,iblk)) > puny .or. abs(vatm(i,j,iblk)) > puny) &
                 worka(i,j) = atan2(uatm(i,j,iblk),vatm(i,j,iblk))*rad_to_deg
                 worka(i,j) = worka(i,j) + c180
           enddo
           enddo
           call accum_hist_field(n_atmdir, iblk, worka(:,:), a2D)
         endif
         if (f_sice   (1:1) /= 'x') then
             do j = jlo, jhi
             do i = ilo, ihi
                worka(i,j) = c0
                do k = 1, nzilyr
                   worka(i,j) = worka(i,j) + trcr(i,j,nt_sice+k-1,iblk)
                enddo
                worka(i,j) = worka(i,j) / nzilyr
             enddo
             enddo
             call accum_hist_field(n_sice,   iblk, worka(:,:), a2D)
         endif

         if (f_fswup(1:1) /= 'x') &
            call accum_hist_field(n_fswup, iblk, &
                 (fsw(:,:,iblk)-fswabs(:,:,iblk)*workb(:,:)), a2D)
         if (f_fswdn  (1:1) /= 'x') &
             call accum_hist_field(n_fswdn,  iblk, fsw(:,:,iblk), a2D)
         if (f_flwdn  (1:1) /= 'x') &
             call accum_hist_field(n_flwdn,  iblk, flw(:,:,iblk), a2D)
         if (f_snow   (1:1) /= 'x') &
             call accum_hist_field(n_snow,   iblk, fsnow(:,:,iblk), a2D)
         if (f_snow_ai(1:1) /= 'x') &
             call accum_hist_field(n_snow_ai,iblk, fsnow(:,:,iblk)*workb(:,:), a2D)
         if (f_rain   (1:1) /= 'x') &
             call accum_hist_field(n_rain,   iblk, frain(:,:,iblk), a2D)
         if (f_rain_ai(1:1) /= 'x') &
             call accum_hist_field(n_rain_ai,iblk, frain(:,:,iblk)*workb(:,:), a2D)

         if (f_sst    (1:1) /= 'x') &
             call accum_hist_field(n_sst,    iblk, sst(:,:,iblk), a2D)
         if (f_sss    (1:1) /= 'x') &
             call accum_hist_field(n_sss,    iblk, sss(:,:,iblk), a2D)
         if (f_uocn   (1:1) /= 'x') &
             call accum_hist_field(n_uocn,   iblk, uocn(:,:,iblk), a2D)
         if (f_vocn   (1:1) /= 'x') &
             call accum_hist_field(n_vocn,   iblk, vocn(:,:,iblk), a2D)
         if (f_ocnspd   (1:1) /= 'x') &
             call accum_hist_field(n_ocnspd,   iblk, sqrt( &
                                  (uocn(:,:,iblk)*uocn(:,:,iblk)) + &
                                  (vocn(:,:,iblk)*vocn(:,:,iblk))), a2D)
         if (f_ocndir(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (abs(uocn(i,j,iblk)) > puny .or. abs(vocn(i,j,iblk)) > puny) &
                 worka(i,j) = atan2(uocn(i,j,iblk),vocn(i,j,iblk))*rad_to_deg
              if (worka(i,j) < 0.0 ) then
                 worka(i,j) = worka(i,j) + c360
              else
                 worka(i,j) = worka(i,j) * c1
              endif
           enddo
           enddo
           call accum_hist_field(n_ocndir, iblk, worka(:,:), a2D)
         endif
         if (f_frzmlt (1:1) /= 'x') &
             call accum_hist_field(n_frzmlt, iblk, frzmlt_init(:,:,iblk), a2D)

         if (f_fswfac (1:1) /= 'x') &
             call accum_hist_field(n_fswfac, iblk, scale_factor(:,:,iblk), a2D)
         if (f_fswabs (1:1) /= 'x') &
             call accum_hist_field(n_fswabs, iblk, fswabs(:,:,iblk), a2D)

         if (f_fswint_ai (1:1) /= 'x') &
             call accum_hist_field(n_fswint_ai, iblk, fswint_ai(:,:,iblk), a2D)

         if (f_fswabs_ai(1:1)/= 'x') &
             call accum_hist_field(n_fswabs_ai, iblk, fswabs(:,:,iblk)*workb(:,:), a2D)

         if (f_albsni (1:1) /= 'x') &
             call accum_hist_field(n_albsni, iblk, &
                                  (awtvdr*alvdr(:,:,iblk) &
                                 + awtidr*alidr(:,:,iblk) &
                                 + awtvdf*alvdf(:,:,iblk) &
                                 + awtidf*alidf(:,:,iblk))*workb(:,:), a2D)
         if (f_alvdr  (1:1) /= 'x') &
             call accum_hist_field(n_alvdr,  iblk, alvdr(:,:,iblk), a2D)
         if (f_alidr  (1:1) /= 'x') &
             call accum_hist_field(n_alidr,  iblk, alidr(:,:,iblk), a2D)
         if (f_alvdf  (1:1) /= 'x') &
             call accum_hist_field(n_alvdf,  iblk, alvdf(:,:,iblk), a2D)
         if (f_alidf  (1:1) /= 'x') &
             call accum_hist_field(n_alidf,  iblk, alidf(:,:,iblk), a2D)
         if (f_alvdr_ai  (1:1) /= 'x') &
             call accum_hist_field(n_alvdr_ai,  iblk, alvdr_ai(:,:,iblk), a2D)
         if (f_alidr_ai  (1:1) /= 'x') &
             call accum_hist_field(n_alidr_ai,  iblk, alidr_ai(:,:,iblk), a2D)
         if (f_alvdf_ai  (1:1) /= 'x') &
             call accum_hist_field(n_alvdf_ai,  iblk, alvdf_ai(:,:,iblk), a2D)
         if (f_alidf_ai  (1:1) /= 'x') &
             call accum_hist_field(n_alidf_ai,  iblk, alidf_ai(:,:,iblk), a2D)

         if (f_albice (1:1) /= 'x') &
             call accum_hist_field(n_albice, iblk, albice(:,:,iblk), a2D)
         if (f_albsno (1:1) /= 'x') &
             call accum_hist_field(n_albsno, iblk, albsno(:,:,iblk), a2D)
         if (f_albpnd (1:1) /= 'x') &
             call accum_hist_field(n_albpnd, iblk, albpnd(:,:,iblk), a2D)
         if (f_coszen (1:1) /= 'x') &
             call accum_hist_field(n_coszen, iblk, coszen(:,:,iblk), a2D)

         if (f_flat   (1:1) /= 'x') &
             call accum_hist_field(n_flat,   iblk, flat(:,:,iblk), a2D)
         if (f_flat_ai(1:1) /= 'x') &
             call accum_hist_field(n_flat_ai,iblk, flat(:,:,iblk)*workb(:,:), a2D)
         if (f_fsens  (1:1) /= 'x') &
             call accum_hist_field(n_fsens,   iblk, fsens(:,:,iblk), a2D)
         if (f_fsens_ai(1:1)/= 'x') &
             call accum_hist_field(n_fsens_ai,iblk, fsens(:,:,iblk)*workb(:,:), a2D)
         if (f_flwup  (1:1) /= 'x') &
             call accum_hist_field(n_flwup,   iblk, flwout(:,:,iblk), a2D)
         if (f_flwup_ai(1:1)/= 'x') &
             call accum_hist_field(n_flwup_ai,iblk, flwout(:,:,iblk)*workb(:,:), a2D)
         if (f_evap   (1:1) /= 'x') &
             call accum_hist_field(n_evap,   iblk, evap(:,:,iblk), a2D)
         if (f_evap_ai(1:1) /= 'x') &
             call accum_hist_field(n_evap_ai,iblk, evap(:,:,iblk)*workb(:,:), a2D)

         if (f_Tair   (1:1) /= 'x') &
             call accum_hist_field(n_Tair,   iblk, Tair(:,:,iblk), a2D)
         if (f_Tref   (1:1) /= 'x') &
             call accum_hist_field(n_Tref,   iblk, Tref(:,:,iblk), a2D)
         if (f_Qref   (1:1) /= 'x') &
             call accum_hist_field(n_Qref,   iblk, Qref(:,:,iblk), a2D)
         if (f_congel (1:1) /= 'x') &
             call accum_hist_field(n_congel, iblk, congel(:,:,iblk), a2D)
         if (f_frazil (1:1) /= 'x') &
             call accum_hist_field(n_frazil, iblk, frazil(:,:,iblk), a2D)
         if (f_snoice (1:1) /= 'x') &
             call accum_hist_field(n_snoice, iblk, snoice(:,:,iblk), a2D)
         if (f_dsnow (1:1) /= 'x') &
             call accum_hist_field(n_dsnow, iblk, dsnow(:,:,iblk), a2D)
         if (f_meltt  (1:1) /= 'x') &
             call accum_hist_field(n_meltt,  iblk, meltt(:,:,iblk), a2D)
         if (f_melts  (1:1) /= 'x') &
              call accum_hist_field(n_melts,  iblk, melts(:,:,iblk), a2D)
         if (f_meltb  (1:1) /= 'x') &
             call accum_hist_field(n_meltb,  iblk, meltb(:,:,iblk), a2D)
         if (f_meltl  (1:1) /= 'x') &
             call accum_hist_field(n_meltl,  iblk, meltl(:,:,iblk), a2D)

         if (f_fresh  (1:1) /= 'x') &
             call accum_hist_field(n_fresh,   iblk, fresh(:,:,iblk), a2D)
         if (f_fresh_ai(1:1)/= 'x') &
             call accum_hist_field(n_fresh_ai,iblk, fresh_ai(:,:,iblk), a2D)
         if (f_fsalt  (1:1) /= 'x') &
             call accum_hist_field(n_fsalt,   iblk, fsalt(:,:,iblk), a2D)
         if (f_fsalt_ai(1:1)/= 'x') &
             call accum_hist_field(n_fsalt_ai,iblk, fsalt_ai(:,:,iblk), a2D)

         if (f_fbot(1:1)/= 'x') &
             call accum_hist_field(n_fbot,iblk, fbot(:,:,iblk), a2D)
         if (f_fhocn  (1:1) /= 'x') &
             call accum_hist_field(n_fhocn,   iblk, fhocn(:,:,iblk), a2D)
         if (f_fhocn_ai(1:1)/= 'x') &
             call accum_hist_field(n_fhocn_ai,iblk, fhocn_ai(:,:,iblk), a2D)
         if (f_fswthru(1:1) /= 'x') &
             call accum_hist_field(n_fswthru, iblk, fswthru(:,:,iblk), a2D)
         if (f_fswthru_ai(1:1)/= 'x') &
             call accum_hist_field(n_fswthru_ai,iblk, fswthru_ai(:,:,iblk), a2D)
               
         if (f_strairx(1:1) /= 'x') &
             call accum_hist_field(n_strairx, iblk, strairx(:,:,iblk), a2D)
         if (f_strairy(1:1) /= 'x') &
             call accum_hist_field(n_strairy, iblk, strairy(:,:,iblk), a2D)
         if (f_strtltx(1:1) /= 'x') &
             call accum_hist_field(n_strtltx, iblk, strtltx(:,:,iblk), a2D)
         if (f_strtlty(1:1) /= 'x') &
             call accum_hist_field(n_strtlty, iblk, strtlty(:,:,iblk), a2D)
         if (f_strcorx(1:1) /= 'x') &
             call accum_hist_field(n_strcorx, iblk, fm(:,:,iblk)*vvel(:,:,iblk), a2D)
         if (f_strcory(1:1) /= 'x') &
             call accum_hist_field(n_strcory, iblk,-fm(:,:,iblk)*uvel(:,:,iblk), a2D)
         if (f_strocnx(1:1) /= 'x') &
             call accum_hist_field(n_strocnx, iblk, strocnx(:,:,iblk), a2D)
         if (f_strocny(1:1) /= 'x') &
             call accum_hist_field(n_strocny, iblk, strocny(:,:,iblk), a2D)
         if (f_strintx(1:1) /= 'x') &
             call accum_hist_field(n_strintx, iblk, strintx(:,:,iblk), a2D)
         if (f_strinty(1:1) /= 'x') &
             call accum_hist_field(n_strinty, iblk, strinty(:,:,iblk), a2D)
         if (f_taubx(1:1) /= 'x') &
             call accum_hist_field(n_taubx, iblk, taubx(:,:,iblk), a2D)
         if (f_tauby(1:1) /= 'x') &
             call accum_hist_field(n_tauby, iblk, tauby(:,:,iblk), a2D)
         if (f_strength(1:1)/= 'x') &
             call accum_hist_field(n_strength,iblk, strength(:,:,iblk), a2D)

! The following fields (divu, shear, sig1, and sig2) will be smeared
!  if averaged over more than a few days.
! Snapshots may be more useful (see below).

!        if (f_divu   (1:1) /= 'x') &
!             call accum_hist_field(n_divu,    iblk, divu(:,:,iblk), a2D)
!        if (f_shear  (1:1) /= 'x') &
!             call accum_hist_field(n_shear,   iblk, shear(:,:,iblk), a2D)
!        if (f_sig1   (1:1) /= 'x') &
!             call accum_hist_field(n_sig1,    iblk, sig1(:,:,iblk), a2D)
!        if (f_sig2   (1:1) /= 'x') &
!             call accum_hist_field(n_sig2,    iblk, sig2(:,:,iblk), a2D)
!        if (f_trsig  (1:1) /= 'x') &
!             call accum_hist_field(n_trsig,   iblk, trsig(:,:,iblk), a2D)

         if (f_dvidtt (1:1) /= 'x') &
             call accum_hist_field(n_dvidtt,  iblk, dvidtt(:,:,iblk), a2D)
         if (f_dvidtd (1:1) /= 'x') &
             call accum_hist_field(n_dvidtd,  iblk, dvidtd(:,:,iblk), a2D)
         if (f_daidtt (1:1) /= 'x') &
             call accum_hist_field(n_daidtt,  iblk, daidtt(:,:,iblk), a2D)
         if (f_daidtd (1:1) /= 'x') &
             call accum_hist_field(n_daidtd,  iblk, daidtd(:,:,iblk), a2D)
         if (f_dagedtt (1:1) /= 'x') &
             call accum_hist_field(n_dagedtt, iblk, dagedtt(:,:,iblk), a2D)
         if (f_dagedtd (1:1) /= 'x') &
             call accum_hist_field(n_dagedtd, iblk, dagedtd(:,:,iblk), a2D)

         if (f_fsurf_ai(1:1)/= 'x') &
             call accum_hist_field(n_fsurf_ai,iblk, fsurf(:,:,iblk)*workb(:,:), a2D)
         if (f_fcondtop_ai(1:1)/= 'x') &
             call accum_hist_field(n_fcondtop_ai, iblk, &
                                                 fcondtop(:,:,iblk)*workb(:,:), a2D)

         if (f_icepresent(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) worka(i,j) = c1
           enddo
           enddo
           call accum_hist_field(n_icepresent, iblk, worka(:,:), a2D)
         endif

         ! 2D CMIP fields

         if (f_sithick(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) worka(i,j) = vice(i,j,iblk)
           enddo
           enddo
           call accum_hist_field(n_sithick, iblk, worka(:,:), a2D)
         endif

         if (f_siage(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) worka(i,j) = aice(i,j,iblk)*trcr(i,j,nt_iage,iblk)
           enddo
           enddo
           call accum_hist_field(n_siage, iblk, worka(:,:), a2D)
         endif

         if (f_sisnthick(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (vsno(i,j,iblk) > puny) &
                 worka(i,j) = vsno(i,j,iblk)
           enddo
           enddo
           call accum_hist_field(n_sisnthick, iblk, worka(:,:), a2D)
         endif

         if (f_sitemptop(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) &
              worka(i,j) = aice(i,j,iblk)*(trcr(i,j,nt_Tsfc,iblk)+Tffresh)
           enddo
           enddo
           call accum_hist_field(n_sitemptop, iblk, worka(:,:), a2D)
         endif

         if (f_sitempsnic(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (vsno(i,j,iblk) > puny .and. aice_init(i,j,iblk) > puny) then
                 worka(i,j) = aice(i,j,iblk)*(Tsnice(i,j,iblk)/aice_init(i,j,iblk)+Tffresh)
              else
                 worka(i,j) = aice(i,j,iblk)*(trcr(i,j,nt_Tsfc,iblk)+Tffresh)
              endif
           enddo
           enddo
           call accum_hist_field(n_sitempsnic, iblk, worka(:,:), a2D)
         endif

         if (f_sitempbot(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice_init(i,j,iblk) > puny) &
                 worka(i,j) = aice(i,j,iblk)*(Tbot(i,j,iblk)/aice_init(i,j,iblk)+Tffresh)
           enddo
           enddo
           call accum_hist_field(n_sitempbot, iblk, worka(:,:), a2D)
         endif

         if (f_siu(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) worka(i,j) = aice(i,j,iblk)*uvel(i,j,iblk)
           enddo
           enddo
           call accum_hist_field(n_siu, iblk, worka(:,:), a2D)
         endif

         if (f_siv(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) worka(i,j) = aice(i,j,iblk)*vvel(i,j,iblk)
           enddo
           enddo
           call accum_hist_field(n_siv, iblk, worka(:,:), a2D)
         endif

         if (f_sispeed(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) worka(i,j) = aice(i,j,iblk) &
                 * sqrt(uvel(i,j,iblk)*uvel(i,j,iblk)+vvel(i,j,iblk)*vvel(i,j,iblk))
           enddo
           enddo
           call accum_hist_field(n_sispeed, iblk, worka(:,:), a2D)
         endif
         if (f_sidir(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (abs(uvel(i,j,iblk)) > puny .or. abs(vvel(i,j,iblk)) > puny) &
                 worka(i,j) = atan2(uvel(i,j,iblk),vvel(i,j,iblk))*rad_to_deg
              if (worka(i,j) < 0.0 ) then
                 worka(i,j) = worka(i,j) + c360
              else
                 worka(i,j) = worka(i,j) * c1
              endif
           enddo
           enddo
           call accum_hist_field(n_sidir, iblk, worka(:,:), a2D)
         endif
         if (f_sidmasstranx(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) &
                 worka(i,j) = (rhoi*p5*(vice(i+1,j,iblk)+vice(i,j,iblk))*dyu(i,j,iblk) &
                            + rhos*p5*(vsno(i+1,j,iblk)+vsno(i,j,iblk))*dyu(i,j,iblk)) &
                            *  p5*(uvel(i,j-1,iblk)+uvel(i,j,iblk))
           enddo
           enddo
           call accum_hist_field(n_sidmasstranx, iblk, worka(:,:), a2D)
         endif

         if (f_sidmasstrany(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) &
                 worka(i,j) = (rhoi*p5*(vice(i,j+1,iblk)+vice(i,j,iblk))*dxu(i,j,iblk) &
                            + rhos*p5*(vsno(i,j+1,iblk)+vsno(i,j,iblk))*dxu(i,j,iblk)) &
                            *  p5*(vvel(i-1,j,iblk)+vvel(i,j,iblk))
           enddo
           enddo
           call accum_hist_field(n_sidmasstrany, iblk, worka(:,:), a2D)
         endif

         if (f_sistrxdtop(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice_init(i,j,iblk) > puny) &
                 worka(i,j) = aice(i,j,iblk)*(aice(i,j,iblk)*strairx(i,j,iblk)/aice_init(i,j,iblk))
           enddo
           enddo
           call accum_hist_field(n_sistrxdtop, iblk, worka(:,:), a2D)
         endif

         if (f_sistrydtop(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice_init(i,j,iblk) > puny) &
                 worka(i,j) = aice(i,j,iblk)*(aice(i,j,iblk)*strairy(i,j,iblk)/aice_init(i,j,iblk))
           enddo
           enddo
           call accum_hist_field(n_sistrydtop, iblk, worka(:,:), a2D)
         endif

         if (f_sistrxubot(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) &
                 worka(i,j) = aice(i,j,iblk)*strocnx(i,j,iblk)
           enddo
           enddo
           call accum_hist_field(n_sistrxubot, iblk, worka(:,:), a2D)
         endif

         if (f_sistryubot(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) &
                 worka(i,j) = aice(i,j,iblk)*strocny(i,j,iblk)
           enddo
           enddo
           call accum_hist_field(n_sistryubot, iblk, worka(:,:), a2D)
         endif

         if (f_sicompstren(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) &
                 worka(i,j) = aice(i,j,iblk)*strength(i,j,iblk)
           enddo
           enddo
           call accum_hist_field(n_sicompstren, iblk, worka(:,:), a2D)
         endif

         if (f_sialb(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (fsw(i,j,iblk) > puny .and. aice_init(i,j,iblk) > puny) then
                 worka(i,j) = aice(i,j,iblk)*(fsw(i,j,iblk)-fswabs(i,j,iblk) &
                            * aice(i,j,iblk)/aice_init(i,j,iblk)) / fsw(i,j,iblk)
              endif
           enddo
           enddo
           call accum_hist_field(n_sialb, iblk, worka(:,:), a2D)
         endif

         if (f_sihc(1:1) /= 'x') then
           worka(:,:) = c0
           do k = 1,nzilyr
           do j = jlo, jhi
           do i = ilo, ihi
              worka(i,j) = worka(i,j) + trcr(i,j,nt_qice+k-1,iblk)*vice(i,j,iblk)/real(nzilyr,kind=dbl_kind)
           enddo
           enddo
           enddo
           call accum_hist_field(n_sihc, iblk, worka(:,:), a2D)
         endif

         if (f_sisnhc(1:1) /= 'x') then
           worka(:,:) = c0
           do k = 1,nzslyr
           do j = jlo, jhi
           do i = ilo, ihi
              worka(i,j) = worka(i,j) + trcr(i,j,nt_qsno+k-1,iblk)*vsno(i,j,iblk)/real(nzslyr,kind=dbl_kind)
           enddo
           enddo
           enddo
           call accum_hist_field(n_sisnhc, iblk, worka(:,:), a2D)
         endif

         if (f_sidconcth(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = daidtt(i,j,iblk)
              endif
           enddo
           enddo
           call accum_hist_field(n_sidconcth, iblk, worka(:,:), a2D)
         endif

         if (f_sidconcdyn(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = daidtd(i,j,iblk)
              endif
           enddo
           enddo
           call accum_hist_field(n_sidconcdyn, iblk, worka(:,:), a2D)
         endif

         if (f_sidmassth(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = dvidtt(i,j,iblk) * rhoi
              endif
           enddo
           enddo
           call accum_hist_field(n_sidmassth, iblk, worka(:,:), a2D)
         endif

         if (f_sidmassdyn(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = dvidtd(i,j,iblk) * rhoi
              endif
           enddo
           enddo
           call accum_hist_field(n_sidmassdyn, iblk, worka(:,:), a2D)
         endif

         if (f_sidmassgrowthwat(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice_init(i,j,iblk) > puny) then
                 worka(i,j) = frazil(i,j,iblk)*rhoi/dt
              endif
           enddo
           enddo
           call accum_hist_field(n_sidmassgrowthwat, iblk, worka(:,:), a2D)
         endif

         if (f_sidmassgrowthbot(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = congel(i,j,iblk)*rhoi/dt
              endif
           enddo
           enddo
           call accum_hist_field(n_sidmassgrowthbot, iblk, worka(:,:), a2D)
         endif

         if (f_sidmasssi(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = snoice(i,j,iblk)*rhoi/dt
              endif
           enddo
           enddo
           call accum_hist_field(n_sidmasssi, iblk, worka(:,:), a2D)
         endif

         if (f_sidmassevapsubl(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = evapi(i,j,iblk)*rhoi
              endif
           enddo
           enddo
           call accum_hist_field(n_sidmassevapsubl, iblk, worka(:,:), a2D)
         endif

         if (f_sidmassmelttop(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = meltt(i,j,iblk)*rhoi/dt
              endif
           enddo
           enddo
           call accum_hist_field(n_sidmassmelttop, iblk, worka(:,:), a2D)
         endif

         if (f_sidmassmeltbot(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = meltb(i,j,iblk)*rhoi/dt
              endif
           enddo
           enddo
           call accum_hist_field(n_sidmassmeltbot, iblk, worka(:,:), a2D)
         endif

         if (f_sidmasslat(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = meltl(i,j,iblk)*rhoi/dt
              endif
           enddo
           enddo
           call accum_hist_field(n_sidmasslat, iblk, worka(:,:), a2D)
         endif

         if (f_sndmasssubl(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = evaps(i,j,iblk)*rhos
              endif
           enddo
           enddo
           call accum_hist_field(n_sndmasssubl, iblk, worka(:,:), a2D)
         endif

         if (f_sndmasssnf(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = aice(i,j,iblk)*fsnow(i,j,iblk)*rhos
              endif
           enddo
           enddo
           call accum_hist_field(n_sndmasssnf, iblk, worka(:,:), a2D)
         endif

         if (f_sndmassmelt(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = melts(i,j,iblk)*rhos/dt
              endif
           enddo
           enddo
           call accum_hist_field(n_sndmassmelt, iblk, worka(:,:), a2D)
         endif

         if (f_siflswdtop(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (fsw(i,j,iblk) > puny .and. aice(i,j,iblk) > puny) then
                 worka(i,j) = aice(i,j,iblk)*fsw(i,j,iblk)
              endif
           enddo
           enddo
           call accum_hist_field(n_siflswdtop, iblk, worka(:,:), a2D)
         endif

         if (f_siflswutop(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (fsw(i,j,iblk) > puny .and. aice_init(i,j,iblk) > puny) then
                 worka(i,j) = aice(i,j,iblk)*(fsw(i,j,iblk)-fswabs(i,j,iblk) &
                            * aice(i,j,iblk)/aice_init(i,j,iblk))
              endif
           enddo
           enddo
           call accum_hist_field(n_siflswutop, iblk, worka(:,:), a2D)
         endif

         if (f_siflswdbot(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = aice(i,j,iblk)*fswthru(i,j,iblk)
              endif
           enddo
           enddo
           call accum_hist_field(n_siflswdbot, iblk, worka(:,:), a2D)
         endif

         if (f_sifllwdtop(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = aice(i,j,iblk)*flw(i,j,iblk)
              endif
           enddo
           enddo
           call accum_hist_field(n_sifllwdtop, iblk, worka(:,:), a2D)
         endif

         if (f_sifllwutop(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = aice(i,j,iblk)*flwout(i,j,iblk)
              endif
           enddo
           enddo
           call accum_hist_field(n_sifllwutop, iblk, worka(:,:), a2D)
         endif

         if (f_siflsenstop(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = aice(i,j,iblk)*fsens(i,j,iblk)
              endif
           enddo
           enddo
           call accum_hist_field(n_siflsenstop, iblk, worka(:,:), a2D)
         endif

         if (f_siflsensupbot(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = aice(i,j,iblk)*fhocn(i,j,iblk)
              endif
           enddo
           enddo
           call accum_hist_field(n_siflsensupbot, iblk, worka(:,:), a2D)
         endif

         if (f_sifllatstop(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = aice(i,j,iblk)*flat(i,j,iblk)
              endif
           enddo
           enddo
           call accum_hist_field(n_sifllatstop, iblk, worka(:,:), a2D)
         endif

         if (f_siflcondtop(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = aice(i,j,iblk)*fcondtop(i,j,iblk)
              endif
           enddo
           enddo
           call accum_hist_field(n_siflcondtop, iblk, worka(:,:), a2D)
         endif

         if (f_siflcondbot(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice_init(i,j,iblk) > puny) then
                 worka(i,j) = aice(i,j,iblk)*fcondbot(i,j,iblk)/aice_init(i,j,iblk)
              endif
           enddo
           enddo
           call accum_hist_field(n_siflcondbot, iblk, worka(:,:), a2D)
         endif

         if (f_sipr(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = aice(i,j,iblk)*frain(i,j,iblk)
              endif
           enddo
           enddo
           call accum_hist_field(n_sipr, iblk, worka(:,:), a2D)
         endif

         if (f_sifb(1:1) /= 'x') then
           worka(:,:) = c0
           rho_ice = rhoi
           rho_ocn = rhow
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 if (ktherm == 2) then
                    rho_ocn = icepack_mushy_density_brine(sss(i,j,iblk))
                    rho_ice = c0
                    do k = 1, nzilyr
                       Tice = icepack_mushy_temperature_mush(trcr(i,j,nt_qice+k-1,iblk),trcr(i,j,nt_sice+k-1,iblk))
                       Sbr = trcr(i,j,nt_sice+k-1,iblk)
                       phi = icepack_mushy_liquid_fraction(Tice,Sbr)
                       rhob = icepack_mushy_density_brine(Sbr)
                       rho_ice = rho_ice + min(phi*rhob+(c1-phi)*rhoi,rho_ocn)
                    enddo
                    rho_ice = rho_ice / real(nzilyr,kind=dbl_kind)
                 endif
                 worka(i,j) = ((rho_ocn-rho_ice)*vice(i,j,iblk) - rhos*vsno(i,j,iblk))/rho_ocn
!                if (worka(i,j) < c0) then
!                   write(nu_diag,*) 'negative fb',rho_ocn,rho_ice,rhos
!                   write(nu_diag,*) vice(i,j,iblk),vsno(i,j,iblk)
!                endif
              endif
           enddo
           enddo
           call accum_hist_field(n_sifb, iblk, worka(:,:), a2D)
         endif

         if (f_siflsaltbot(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
!                Add in frazil flux
                 if (.not. update_ocn_f) then
                 if ( ktherm == 2) then
                    dfresh = -rhoi*(frazil(i,j,iblk)-frazil_diag(i,j,iblk))/dt
                 else
                    dfresh = -rhoi*frazil(i,j,iblk)/dt 
                 endif
                 endif
                 dfsalt = ice_ref_salinity*p001*dfresh
                 worka(i,j) = aice(i,j,iblk)*(fsalt(i,j,iblk)+dfsalt)
              endif
           enddo
           enddo
           call accum_hist_field(n_siflsaltbot, iblk, worka(:,:), a2D)
         endif

         if (f_siflfwbot(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
!                Add in frazil flux
!                Add in frazil flux
                 if (.not. update_ocn_f) then
                 if ( ktherm == 2) then
                    dfresh = -rhoi*(frazil(i,j,iblk)-frazil_diag(i,j,iblk))/dt
                 else
                    dfresh = -rhoi*frazil(i,j,iblk)/dt 
                 endif
                 endif
                 worka(i,j) = aice(i,j,iblk)*(fresh(i,j,iblk)+dfresh)
              endif
           enddo
           enddo
           call accum_hist_field(n_siflfwbot, iblk, worka(:,:), a2D)
         endif

         if (f_siflfwdrain(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = aice(i,j,iblk)*(frain(i,j,iblk)+melts(i,j,iblk)+meltt(i,j,iblk))
              endif
           enddo
           enddo
           call accum_hist_field(n_siflfwdrain, iblk, worka(:,:), a2D)
         endif

         if (f_sidragtop(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = aice(i,j,iblk)*Cdn_atm(i,j,iblk)
              endif
           enddo
           enddo
           call accum_hist_field(n_sidragtop, iblk, worka(:,:), a2D)
         endif

         if (f_sirdgthick(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk)*(c1 - trcr(i,j,nt_alvl,iblk)) > puny) then
                 worka(i,j) = vice(i,j,iblk) * (c1 - trcr(i,j,nt_vlvl,iblk)) &
                           / (aice(i,j,iblk) * (c1 - trcr(i,j,nt_alvl,iblk)))
              endif
           enddo
           enddo
           call accum_hist_field(n_sirdgthick, iblk, worka(:,:), a2D)
         endif

         if (f_siforcetiltx(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = aice(i,j,iblk)*strtltx(i,j,iblk)
              endif
           enddo
           enddo
           call accum_hist_field(n_siforcetiltx, iblk, worka(:,:), a2D)
         endif

         if (f_siforcetilty(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = aice(i,j,iblk)*strtlty(i,j,iblk)
              endif
           enddo
           enddo
           call accum_hist_field(n_siforcetilty, iblk, worka(:,:), a2D)
         endif

         if (f_siforcecoriolx(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = aice(i,j,iblk)*fm(i,j,iblk)*vvel(i,j,iblk)
              endif
           enddo
           enddo
           call accum_hist_field(n_siforcecoriolx, iblk, worka(:,:), a2D)
         endif

         if (f_siforcecorioly(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = -aice(i,j,iblk)*fm(i,j,iblk)*uvel(i,j,iblk)
              endif
           enddo
           enddo
           call accum_hist_field(n_siforcecorioly, iblk, worka(:,:), a2D)
         endif

         if (f_siforceintstrx(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = aice(i,j,iblk)*strintx(i,j,iblk)
              endif
           enddo
           enddo
           call accum_hist_field(n_siforceintstrx, iblk, worka(:,:), a2D)
         endif

         if (f_siforceintstry(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) then
                 worka(i,j) = aice(i,j,iblk)*strinty(i,j,iblk)
              endif
           enddo
           enddo
           call accum_hist_field(n_siforceintstry, iblk, worka(:,:), a2D)
         endif

         ! 3D category fields
         if (f_aicen   (1:1) /= 'x') &
             call accum_hist_field(n_aicen-n2D, iblk, ncat_hist, &
                                   aicen(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_vicen   (1:1) /= 'x') &
             call accum_hist_field(n_vicen-n2D, iblk, ncat_hist, &
                                   vicen(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_vsnon   (1:1) /= 'x') &
             call accum_hist_field(n_vsnon-n2D, iblk, ncat_hist, &
                                   vsnon(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_snowfracn(1:1) /= 'x') &
             call accum_hist_field(n_snowfracn-n2D, iblk, ncat_hist, &
                                   snowfracn(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_snowfracn(1:1) /= 'x') &
             call accum_hist_field(n_snowfracn-n2D, iblk, ncat_hist, &
                                   snowfracn(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_keffn_top (1:1) /= 'x') &
             call accum_hist_field(n_keffn_top-n2D, iblk, ncat_hist, &
                                   keffn_top(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_fsurfn_ai   (1:1) /= 'x') &
             call accum_hist_field(n_fsurfn_ai-n2D, iblk, ncat_hist, &
                  fsurfn(:,:,1:ncat_hist,iblk)*aicen_init(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_fcondtopn_ai   (1:1) /= 'x') &
             call accum_hist_field(n_fcondtopn_ai-n2D, iblk, ncat_hist, &
                  fcondtopn(:,:,1:ncat_hist,iblk)*aicen_init(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_flatn_ai   (1:1) /= 'x') &
             call accum_hist_field(n_flatn_ai-n2D, iblk, ncat_hist, &
                  flatn(:,:,1:ncat_hist,iblk)*aicen_init(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_fsensn_ai   (1:1) /= 'x') &
             call accum_hist_field(n_fsensn_ai-n2D, iblk, ncat_hist, &
                  fsensn(:,:,1:ncat_hist,iblk)*aicen_init(:,:,1:ncat_hist,iblk), a3Dc)
         ! Calculate surface heat flux that causes melt (calculated by the 
         ! atmos in HadGEM3 so needed for checking purposes)
         if (f_fmelttn_ai   (1:1) /= 'x') &
             call accum_hist_field(n_fmelttn_ai-n2D, iblk, ncat_hist, &
                  max(fsurfn(:,:,1:ncat_hist,iblk) - fcondtopn(:,:,1:ncat_hist,iblk),c0) &
                      *aicen_init(:,:,1:ncat_hist,iblk), a3Dc)

         if (f_siitdconc   (1:1) /= 'x') then
           worka3(:,:,:) = c0
           do n = 1,ncat_hist
           do j = jlo, jhi
           do i = ilo, ihi
              if (aicen(i,j,n,iblk) > puny) then
                 worka3(i,j,n) = aicen(i,j,n,iblk)
              endif
           enddo
           enddo
           enddo
           call accum_hist_field(n_siitdconc-n2D, iblk, ncat_hist, worka3(:,:,:), a3Dc)
         endif

         if (f_siitdthick   (1:1) /= 'x') then
           worka3(:,:,:) = c0
           do n = 1,ncat_hist
           do j = jlo, jhi
           do i = ilo, ihi
              if (aicen(i,j,n,iblk) > puny) then
                 worka3(i,j,n) = vicen(i,j,n,iblk)
              endif
           enddo
           enddo
           enddo
           call accum_hist_field(n_siitdthick-n2D, iblk, ncat_hist, worka3(:,:,:), a3Dc)
         endif

         if (f_siitdsnthick   (1:1) /= 'x') then
           worka3(:,:,:) = c0
           do n = 1,ncat_hist
           do j = jlo, jhi
           do i = ilo, ihi
              if (aicen(i,j,n,iblk) > puny) then
                 worka3(i,j,n) = vsnon(i,j,n,iblk)
              endif
           enddo
           enddo
           enddo
           call accum_hist_field(n_siitdsnthick-n2D, iblk, ncat_hist, worka3(:,:,:), a3Dc)
         endif

! example for 3D field (x,y,z)
!         if (f_field3dz   (1:1) /= 'x') &
!             call accum_hist_field(n_field3dz-n3Dccum, iblk, nzilyr, &
!                                   field3dz(:,:,1:nzilyr,iblk), a3Dz)

         ! 4D category fields
         if (f_Tinz   (1:1) /= 'x') then
            Tinz4d(:,:,:,:) = c0
            do n = 1, ncat_hist
               do j = jlo, jhi
               do i = ilo, ihi
                  do k = 1, nzilyr
                     qn = trcrn(i,j,nt_qice+k-1,n,iblk)
                     sn = trcrn(i,j,nt_sice+k-1,n,iblk)
                     Tinz4d(i,j,k,n) = icepack_ice_temperature(qn,sn)
                  enddo
               enddo
               enddo
            enddo
            call accum_hist_field(n_Tinz-n3Dfcum, iblk, nzilyr, ncat_hist, &
                                  Tinz4d(:,:,1:nzilyr,1:ncat_hist), a4Di)
         endif
         if (f_Sinz   (1:1) /= 'x') then
            Sinz4d(:,:,:,:) = c0
            do n = 1, ncat_hist
               do j = jlo, jhi
               do i = ilo, ihi
                  if (vicen(i,j,n,iblk) > puny) then
                     Sinz4d(i,j,1:nzilyr,n) = trcrn(i,j,nt_sice:nt_sice+nzilyr-1,n,iblk)
                  endif
               enddo
               enddo
            enddo
            call accum_hist_field(n_Sinz-n3Dfcum, iblk, nzilyr, ncat_hist, &
                                  Sinz4d(:,:,1:nzilyr,1:ncat_hist), a4Di)
         endif
         
         if (f_Tsnz   (1:1) /= 'x') then
            Tsnz4d(:,:,:,:) = c0
            do n = 1, ncat_hist
               do j = jlo, jhi
               do i = ilo, ihi
                  do k = 1, nzslyr
                     qn = trcrn(i,j,nt_qsno+k-1,n,iblk)
                     Tsnz4d(i,j,k,n) = icepack_snow_temperature(qn)
                  enddo
               enddo
               enddo
            enddo
            call accum_hist_field(n_Tsnz-n4Dicum, iblk, nzslyr, ncat_hist, &
                                  Tsnz4d(:,:,1:nzslyr,1:ncat_hist), a4Ds)
         endif
         
        ! Calculate aggregate surface melt flux by summing category values
        if (f_fmeltt_ai(1:1) /= 'x') then
         do ns = 1, nstreams
           if (n_fmeltt_ai(ns) /= 0) then
              worka(:,:) = c0
              do j = jlo, jhi
              do i = ilo, ihi
               if (tmask(i,j,iblk)) then
                 do n=1,ncat_hist
                    worka(i,j)  = worka(i,j) + a3Dc(i,j,n,n_fmelttn_ai(ns)-n2D,iblk)
                 enddo            ! n
               endif              ! tmask
              enddo                ! i
              enddo                ! j
              a2D(:,:,n_fmeltt_ai(ns),iblk) = worka(:,:)
           endif
         enddo
        endif

      !---------------------------------------------------------------
      ! accumulate other history output
      !---------------------------------------------------------------

         ! mechanical redistribution
         call accum_hist_mechred (iblk)

         ! melt ponds
         call accum_hist_pond (iblk)

         ! biogeochemistry
         call accum_hist_bgc (iblk)

         ! form drag
         call accum_hist_drag (iblk)

         ! floe size distribution
         call accum_hist_fsd (iblk)

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !---------------------------------------------------------------
      ! Write output files at prescribed intervals
      !---------------------------------------------------------------

      nstrm = nstreams
      if (write_ic) nstrm = 1

      do ns = 1, nstrm
      if (write_history(ns) .or. write_ic) then

      !---------------------------------------------------------------
      ! Mask out land points and convert units 
      !---------------------------------------------------------------

        ravgct = c1/avgct(ns)
        !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block, &
        !$OMP                     n,nn,ravgctz,ravgip,ravgipn)
        do iblk = 1, nblocks
           this_block = get_block(blocks_ice(iblk),iblk)         
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi

           ! Ice fraction really needs to be on one of the history
           ! streams, but in case it is not.

           if (n_aice(ns) > 0) then
           do j = jlo, jhi
           do i = ilo, ihi
              if (a2D(i,j,n_aice(ns),iblk) > puny) then
                 ravgip(i,j) = c1/(a2D(i,j,n_aice(ns),iblk))
              else
                 ravgip(i,j) = c0
              endif
           enddo             ! i
           enddo             ! j
           endif
           if (n_aicen(ns) > n2D) then
           do k=1,ncat_hist
           do j = jlo, jhi
           do i = ilo, ihi
              if (a3Dc(i,j,k,n_aicen(ns)-n2D,iblk) > puny) then
                 ravgipn(i,j,k) = c1/(a3Dc(i,j,k,n_aicen(ns)-n2D,iblk))
              else
                 ravgipn(i,j,k) = c0
              endif
           enddo             ! i
           enddo             ! j
           enddo             ! k
           endif

           do n = 1, num_avail_hist_fields_2D
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) then 

              do j = jlo, jhi
              do i = ilo, ihi
                 if (.not. tmask(i,j,iblk)) then ! mask out land points
                    a2D(i,j,n,iblk) = spval
                 else                            ! convert units
                    a2D(i,j,n,iblk) = avail_hist_fields(n)%cona*a2D(i,j,n,iblk) &
                                   * ravgct + avail_hist_fields(n)%conb
                 endif
              enddo             ! i
              enddo             ! j

              ! Only average for timesteps when ice present
              if (index(avail_hist_fields(n)%vname,'sithick') /= 0) then
                 if (f_sithick(1:1) /= 'x' .and. n_sithick(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_sithick(ns),iblk) = &
                             a2D(i,j,n_sithick(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_sithick(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'siage') /= 0) then
                 if (f_siage(1:1) /= 'x' .and. n_siage(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_siage(ns),iblk) = &
                             a2D(i,j,n_siage(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_siage(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'sisnthick') /= 0) then
                 if (f_sisnthick(1:1) /= 'x' .and. n_sisnthick(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_sisnthick(ns),iblk) = &
                             a2D(i,j,n_sisnthick(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_sisnthick(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'sitemptop') /= 0) then
                 if (f_sitemptop(1:1) /= 'x' .and. n_sitemptop(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_sitemptop(ns),iblk) = &
                             a2D(i,j,n_sitemptop(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_sitemptop(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'sitempsnic') /= 0) then
                 if (f_sitempsnic(1:1) /= 'x' .and. n_sitempsnic(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_sitempsnic(ns),iblk) = &
                             a2D(i,j,n_sitempsnic(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_sitempsnic(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'sitempbot') /= 0) then
                 if (f_sitempbot(1:1) /= 'x' .and. n_sitempbot(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_sitempbot(ns),iblk) = &
                             a2D(i,j,n_sitempbot(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_sitempbot(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'siu') /= 0) then
                 if (f_siu(1:1) /= 'x' .and. n_siu(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_siu(ns),iblk) = &
                             a2D(i,j,n_siu(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_siu(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'siv') /= 0) then
                 if (f_siv(1:1) /= 'x' .and. n_siv(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_siv(ns),iblk) = &
                             a2D(i,j,n_siv(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_siv(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'sistrxdtop') /= 0) then
                 if (f_sistrxdtop(1:1) /= 'x' .and. n_sistrxdtop(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_sistrxdtop(ns),iblk) = &
                             a2D(i,j,n_sistrxdtop(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_sistrxdtop(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'sistrydtop') /= 0) then
                 if (f_sistrydtop(1:1) /= 'x' .and. n_sistrydtop(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_sistrydtop(ns),iblk) = &
                             a2D(i,j,n_sistrydtop(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_sistrydtop(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'sistrxubot') /= 0) then
                 if (f_sistrxubot(1:1) /= 'x' .and. n_sistrxubot(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_sistrxubot(ns),iblk) = &
                             a2D(i,j,n_sistrxubot(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_sistrxubot(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'sistryubot') /= 0) then
                 if (f_sistryubot(1:1) /= 'x' .and. n_sistryubot(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_sistryubot(ns),iblk) = &
                             a2D(i,j,n_sistryubot(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_sistryubot(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'sicompstren') /= 0) then
                 if (f_sicompstren(1:1) /= 'x' .and. n_sicompstren(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_sicompstren(ns),iblk) = &
                             a2D(i,j,n_sicompstren(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_sicompstren(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'sispeed') /= 0) then
                 if (f_sispeed(1:1) /= 'x' .and. n_sispeed(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_sispeed(ns),iblk) = &
                             a2D(i,j,n_sispeed(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_sispeed(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'sialb') /= 0) then
                 if (f_sialb(1:1) /= 'x' .and. n_sialb(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_sialb(ns),iblk) = &
                             a2D(i,j,n_sialb(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_sialb(ns),iblk) = spval
                       if (albcnt(i,j,iblk,ns) <= puny) a2D(i,j,n_sialb(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'siflswdtop') /= 0) then
                 if (f_siflswdtop(1:1) /= 'x' .and. n_siflswdtop(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_siflswdtop(ns),iblk) = &
                             a2D(i,j,n_siflswdtop(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_siflswdtop(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'siflswutop') /= 0) then
                 if (f_siflswutop(1:1) /= 'x' .and. n_siflswutop(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_siflswutop(ns),iblk) = &
                             a2D(i,j,n_siflswutop(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_siflswutop(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'siflswdbot') /= 0) then
                 if (f_siflswdbot(1:1) /= 'x' .and. n_siflswdbot(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_siflswdbot(ns),iblk) = &
                             a2D(i,j,n_siflswdbot(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_siflswdbot(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'sifllwdtop') /= 0) then
                 if (f_sifllwdtop(1:1) /= 'x' .and. n_sifllwdtop(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_sifllwdtop(ns),iblk) = &
                             a2D(i,j,n_sifllwdtop(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_sifllwdtop(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'sifllwutop') /= 0) then
                 if (f_sifllwutop(1:1) /= 'x' .and. n_sifllwutop(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_sifllwutop(ns),iblk) = &
                             a2D(i,j,n_sifllwutop(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_sifllwutop(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'siflsenstop') /= 0) then
                 if (f_siflsenstop(1:1) /= 'x' .and. n_siflsenstop(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_siflsenstop(ns),iblk) = &
                             a2D(i,j,n_siflsenstop(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_siflsenstop(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'siflsensupbot') /= 0) then
                 if (f_siflsensupbot(1:1) /= 'x' .and.  n_siflsensupbot(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_siflsensupbot(ns),iblk) = &
                             a2D(i,j,n_siflsensupbot(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_siflsensupbot(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'sifllatstop') /= 0) then
                 if (f_sifllatstop(1:1) /= 'x' .and. n_sifllatstop(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_sifllatstop(ns),iblk) = &
                             a2D(i,j,n_sifllatstop(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_sifllatstop(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'sipr') /= 0) then
                 if (f_sipr(1:1) /= 'x' .and. n_sipr(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_sipr(ns),iblk) = &
                             a2D(i,j,n_sipr(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_sipr(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'sifb') /= 0) then
                 if (f_sifb(1:1) /= 'x' .and. n_sifb(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_sifb(ns),iblk) = &
                             a2D(i,j,n_sifb(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_sifb(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'siflcondtop') /= 0) then
                 if (f_siflcondtop(1:1) /= 'x' .and. n_siflcondtop(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_siflcondtop(ns),iblk) = &
                             a2D(i,j,n_siflcondtop(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_siflcondtop(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'siflcondbot') /= 0) then
                 if (f_siflcondbot(1:1) /= 'x' .and. n_siflcondbot(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_siflcondbot(ns),iblk) = &
                             a2D(i,j,n_siflcondbot(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_siflcondbot(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'siflsaltbot') /= 0) then
                 if (f_siflsaltbot(1:1) /= 'x' .and. n_siflsaltbot(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_siflsaltbot(ns),iblk) = &
                             a2D(i,j,n_siflsaltbot(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_siflsaltbot(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'siflfwbot') /= 0) then
                 if (f_siflfwbot(1:1) /= 'x' .and. n_siflfwbot(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_siflfwbot(ns),iblk) = &
                             a2D(i,j,n_siflfwbot(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_siflfwbot(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'siflfwdrain') /= 0) then
                 if (f_siflfwdrain(1:1) /= 'x' .and. n_siflfwdrain(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_siflfwdrain(ns),iblk) = &
                             a2D(i,j,n_siflfwdrain(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_siflfwdrain(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
             if (index(avail_hist_fields(n)%vname,'sidragtop') /= 0) then
                 if (f_sidragtop(1:1) /= 'x' .and. n_sidragtop(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_sidragtop(ns),iblk) = &
                             a2D(i,j,n_sidragtop(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_sidragtop(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'sirdgthick') /= 0) then
                 if (f_sirdgthick(1:1) /= 'x' .and. n_sirdgthick(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_sirdgthick(ns),iblk) = &
                             a2D(i,j,n_sirdgthick(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_sirdgthick(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'siforcetiltx') /= 0) then
                 if (f_siforcetiltx(1:1) /= 'x' .and. n_siforcetiltx(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_siforcetiltx(ns),iblk) = &
                             a2D(i,j,n_siforcetiltx(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_siforcetiltx(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'siforcetilty') /= 0) then
                 if (f_siforcetilty(1:1) /= 'x' .and. n_siforcetilty(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_siforcetilty(ns),iblk) = &
                             a2D(i,j,n_siforcetilty(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_siforcetilty(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'siforcecoriolx') /= 0) then
                 if (f_siforcecoriolx(1:1) /= 'x' .and.  n_siforcecoriolx(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_siforcecoriolx(ns),iblk) = &
                             a2D(i,j,n_siforcecoriolx(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_siforcecoriolx(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'siforcecorioly') /= 0) then
                 if (f_siforcecorioly(1:1) /= 'x' .and.  n_siforcecorioly(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_siforcecorioly(ns),iblk) = &
                             a2D(i,j,n_siforcecorioly(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_siforcecorioly(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'siforceintstrx') /= 0) then
                 if (f_siforceintstrx(1:1) /= 'x' .and.  n_siforceintstrx(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_siforceintstrx(ns),iblk) = &
                             a2D(i,j,n_siforceintstrx(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_siforceintstrx(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif
              if (index(avail_hist_fields(n)%vname,'siforceintstry') /= 0) then
                 if (f_siforceintstry(1:1) /= 'x' .and.  n_siforceintstry(ns) /= 0) then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a2D(i,j,n_siforceintstry(ns),iblk) = &
                             a2D(i,j,n_siforceintstry(ns),iblk)*avgct(ns)*ravgip(i,j)
                       endif
                       if (ravgip(i,j) == c0) a2D(i,j,n_siforceintstry(ns),iblk) = spval
                    enddo             ! i
                    enddo             ! j
                 endif
              endif

              ! back out albedo/zenith angle dependence
              if (avail_hist_fields(n)%vname(1:6) == 'albice') then
              do j = jlo, jhi
              do i = ilo, ihi
                 if (tmask(i,j,iblk)) then 
                    ravgctz = c0
                    if (albcnt(i,j,iblk,ns) > puny) &
                        ravgctz = c1/albcnt(i,j,iblk,ns)
                    if (f_albice (1:1) /= 'x' .and. n_albice(ns) /= 0) &
                       a2D(i,j,n_albice(ns),iblk) = &
                       a2D(i,j,n_albice(ns),iblk)*avgct(ns)*ravgctz
                    if (f_albsno (1:1) /= 'x' .and. n_albsno(ns) /= 0) &
                       a2D(i,j,n_albsno(ns),iblk) = &
                       a2D(i,j,n_albsno(ns),iblk)*avgct(ns)*ravgctz
                    if (f_albpnd (1:1) /= 'x' .and. n_albpnd(ns) /= 0) &
                       a2D(i,j,n_albpnd(ns),iblk) = &
                       a2D(i,j,n_albpnd(ns),iblk)*avgct(ns)*ravgctz
                 endif
              enddo             ! i
              enddo             ! j
              endif
              if (avail_hist_fields(n)%vname(1:6) == 'albsni') then
              do j = jlo, jhi
              do i = ilo, ihi
                 if (tmask(i,j,iblk)) then 
                    ravgctz = c0
                    if (albcnt(i,j,iblk,ns) > puny) &
                        ravgctz = c1/albcnt(i,j,iblk,ns)
                    if (f_albsni (1:1) /= 'x' .and. n_albsni(ns) /= 0) &
                       a2D(i,j,n_albsni(ns),iblk) = &
                       a2D(i,j,n_albsni(ns),iblk)*avgct(ns)*ravgctz
                 endif
              enddo             ! i
              enddo             ! j
              endif
              if (avail_hist_fields(n)%vname(1:8) == 'alvdr_ai') then
              do j = jlo, jhi
              do i = ilo, ihi
                 if (tmask(i,j,iblk)) then 
                    ravgctz = c0
                    if (albcnt(i,j,iblk,ns) > puny) &
                        ravgctz = c1/albcnt(i,j,iblk,ns)
                    if (f_alvdr_ai (1:1) /= 'x' .and. n_alvdr_ai(ns) /= 0) &
                       a2D(i,j,n_alvdr_ai(ns),iblk) = &
                       a2D(i,j,n_alvdr_ai(ns),iblk)*avgct(ns)*ravgctz
                    if (f_alvdf_ai (1:1) /= 'x' .and. n_alvdf_ai(ns) /= 0) &
                       a2D(i,j,n_alvdf_ai(ns),iblk) = &
                       a2D(i,j,n_alvdf_ai(ns),iblk)*avgct(ns)*ravgctz
                    if (f_alidr_ai (1:1) /= 'x' .and. n_alidr_ai(ns) /= 0) &
                       a2D(i,j,n_alidr_ai(ns),iblk) = &
                       a2D(i,j,n_alidr_ai(ns),iblk)*avgct(ns)*ravgctz
                    if (f_alidf_ai (1:1) /= 'x' .and. n_alidf_ai(ns) /= 0) &
                       a2D(i,j,n_alidf_ai(ns),iblk) = &
                       a2D(i,j,n_alidf_ai(ns),iblk)*avgct(ns)*ravgctz
                 endif
              enddo             ! i
              enddo             ! j
              endif

              endif
           enddo                ! n

           do n = 1, num_avail_hist_fields_3Dc
              nn = n2D + n
              if (avail_hist_fields(nn)%vhistfreq == histfreq(ns)) then 

              do k = 1, ncat_hist
              do j = jlo, jhi
              do i = ilo, ihi
                 if (.not. tmask(i,j,iblk)) then ! mask out land points
                    a3Dc(i,j,k,n,iblk) = spval
                 else                            ! convert units
                    a3Dc(i,j,k,n,iblk) = avail_hist_fields(nn)%cona*a3Dc(i,j,k,n,iblk) &
                                   * ravgct + avail_hist_fields(nn)%conb
                 endif
              enddo             ! i
              enddo             ! j
              enddo             ! k
              if (index(avail_hist_fields(nn)%vname,'siitdthick') /= 0) then
                 if (f_siitdthick(1:1) /= 'x' .and. n_siitdthick(ns)-n2D /= 0) then
                    do k = 1, ncat_hist
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a3Dc(i,j,k,n_siitdthick(ns)-n2D,iblk) = &
                             a3Dc(i,j,k,n_siitdthick(ns)-n2D,iblk)*avgct(ns)*ravgipn(i,j,k)
                       endif
                    enddo             ! i
                    enddo             ! j
                    enddo             ! k
                 endif
              endif
              if (index(avail_hist_fields(nn)%vname,'siitdsnthick') /= 0) then
                 if (f_siitdsnthick(1:1) /= 'x' .and.  n_siitdsnthick(ns)-n2D /= 0) then
                    do k = 1, ncat_hist
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then
                             a3Dc(i,j,k,n_siitdsnthick(ns)-n2D,iblk) = &
                             a3Dc(i,j,k,n_siitdsnthick(ns)-n2D,iblk)*avgct(ns)*ravgipn(i,j,k)
                       endif
                    enddo             ! i
                    enddo             ! j
                    enddo             ! k
                 endif
              endif

              endif

           enddo                ! n

           do n = 1, num_avail_hist_fields_3Dz
              nn = n3Dccum + n
              if (avail_hist_fields(nn)%vhistfreq == histfreq(ns)) then 

              do k = 1, nzilyr
              do j = jlo, jhi
              do i = ilo, ihi
                 if (.not. tmask(i,j,iblk)) then ! mask out land points
                    a3Dz(i,j,k,n,iblk) = spval
                 else                            ! convert units
                    a3Dz(i,j,k,n,iblk) = avail_hist_fields(nn)%cona*a3Dz(i,j,k,n,iblk) &
                                   * ravgct + avail_hist_fields(nn)%conb
                 endif
              enddo             ! i
              enddo             ! j
              enddo             ! k
              endif
           enddo                ! n
           do n = 1, num_avail_hist_fields_3Db
              nn = n3Dzcum + n
              if (avail_hist_fields(nn)%vhistfreq == histfreq(ns)) then 
              do k = 1, nzblyr
              do j = jlo, jhi
              do i = ilo, ihi
                 if (.not. tmask(i,j,iblk)) then ! mask out land points
                    a3Db(i,j,k,n,iblk) = spval
                 else                            ! convert units
                    a3Db(i,j,k,n,iblk) = avail_hist_fields(nn)%cona*a3Db(i,j,k,n,iblk) &
                                   * ravgct + avail_hist_fields(nn)%conb
                 endif
              enddo             ! i
              enddo             ! j
              enddo             ! k
              endif
           enddo                ! n

           do n = 1, num_avail_hist_fields_3Da
              nn = n3Dbcum + n
              if (avail_hist_fields(nn)%vhistfreq == histfreq(ns)) then 
              do k = 1, nzalyr
              do j = jlo, jhi
              do i = ilo, ihi
                 if (.not. tmask(i,j,iblk)) then ! mask out land points
                    a3Da(i,j,k,n,iblk) = spval
                 else                            ! convert units
                    a3Da(i,j,k,n,iblk) = avail_hist_fields(nn)%cona*a3Da(i,j,k,n,iblk) &
                                   * ravgct + avail_hist_fields(nn)%conb
                 endif
              enddo             ! i
              enddo             ! j
              enddo             ! k
              endif
           enddo                ! n

           do n = 1, num_avail_hist_fields_3Df
              nn = n3Dacum + n
              if (avail_hist_fields(nn)%vhistfreq == histfreq(ns)) then
              do k = 1, nfsd
              do j = jlo, jhi
              do i = ilo, ihi
                 if (.not. tmask(i,j,iblk)) then ! mask out land points
                    a3Df(i,j,k,n,iblk) = spval
                 else                            ! convert units
                    a3Df(i,j,k,n,iblk) = avail_hist_fields(nn)%cona*a3Df(i,j,k,n,iblk) &
                                   * ravgct + avail_hist_fields(nn)%conb
                 endif
              enddo             ! i
              enddo             ! j
              enddo             ! k
              endif
           enddo                ! n

           do n = 1, num_avail_hist_fields_4Di
              nn = n3Dfcum + n
              if (avail_hist_fields(nn)%vhistfreq == histfreq(ns)) then 
              do k = 1, nzilyr
              do ic = 1, ncat_hist
              do j = jlo, jhi
              do i = ilo, ihi
                 if (.not. tmask(i,j,iblk)) then ! mask out land points
                    a4Di(i,j,k,ic,n,iblk) = spval
                 else                            ! convert units
                    a4Di(i,j,k,ic,n,iblk) = avail_hist_fields(nn)%cona*a4Di(i,j,k,ic,n,iblk) &
                                   * ravgct + avail_hist_fields(nn)%conb
                 endif
              enddo             ! i
              enddo             ! j
              enddo             ! ic
              enddo             ! k
              endif
           enddo                ! n

           do n = 1, num_avail_hist_fields_4Ds
              nn = n4Dicum + n
              if (avail_hist_fields(nn)%vhistfreq == histfreq(ns)) then 
              do k = 1, nzslyr
              do ic = 1, ncat_hist
              do j = jlo, jhi
              do i = ilo, ihi
                 if (.not. tmask(i,j,iblk)) then ! mask out land points
                    a4Ds(i,j,k,ic,n,iblk) = spval
                 else                            ! convert units
                    a4Ds(i,j,k,ic,n,iblk) = avail_hist_fields(nn)%cona*a4Ds(i,j,k,ic,n,iblk) &
                                   * ravgct + avail_hist_fields(nn)%conb
                 endif
              enddo             ! i
              enddo             ! j
              enddo             ! ic
              enddo             ! k
              endif
           enddo                ! n

           do n = 1, num_avail_hist_fields_4Df
              nn = n4Dscum + n
              if (avail_hist_fields(nn)%vhistfreq == histfreq(ns)) then
              do k = 1, nfsd
              do ic = 1, ncat_hist
              do j = jlo, jhi
              do i = ilo, ihi
                 if (.not. tmask(i,j,iblk)) then ! mask out land points
                    a4Df(i,j,k,ic,n,iblk) = spval
                 else                            ! convert units
                    a4Df(i,j,k,ic,n,iblk) = avail_hist_fields(nn)%cona*a4Df(i,j,k,ic,n,iblk) &
                                   * ravgct + avail_hist_fields(nn)%conb
                 endif
              enddo             ! i
              enddo             ! j
              enddo             ! ic
              enddo             ! k
              endif
           enddo                ! n

      !---------------------------------------------------------------
      ! snapshots
      !---------------------------------------------------------------

          ! compute sig1 and sig2
        
           call principal_stress (nx_block,  ny_block,  &
                                  stressp_1 (:,:,iblk), &
                                  stressm_1 (:,:,iblk), &
                                  stress12_1(:,:,iblk), &
                                  strength  (:,:,iblk), &
                                  sig1      (:,:,iblk), &
                                  sig2      (:,:,iblk), &
                                  sigP      (:,:,iblk))
 
           do j = jlo, jhi
           do i = ilo, ihi
              if (.not. tmask(i,j,iblk)) then ! mask out land points
                 if (n_divu     (ns) /= 0) a2D(i,j,n_divu(ns),     iblk) = spval
                 if (n_shear    (ns) /= 0) a2D(i,j,n_shear(ns),    iblk) = spval
                 if (n_sig1     (ns) /= 0) a2D(i,j,n_sig1(ns),     iblk) = spval
                 if (n_sig2     (ns) /= 0) a2D(i,j,n_sig2(ns),     iblk) = spval
                 if (n_sigP     (ns) /= 0) a2D(i,j,n_sigP(ns),     iblk) = spval
                 if (n_sistreave(ns) /= 0) a2D(i,j,n_sistreave(ns),iblk) = spval
                 if (n_sistremax(ns) /= 0) a2D(i,j,n_sistremax(ns),iblk) = spval
                 if (n_mlt_onset(ns) /= 0) a2D(i,j,n_mlt_onset(ns),iblk) = spval
                 if (n_frz_onset(ns) /= 0) a2D(i,j,n_frz_onset(ns),iblk) = spval
                 if (n_hisnap   (ns) /= 0) a2D(i,j,n_hisnap(ns),   iblk) = spval
                 if (n_aisnap   (ns) /= 0) a2D(i,j,n_aisnap(ns),   iblk) = spval
                 if (n_trsig    (ns) /= 0) a2D(i,j,n_trsig(ns),    iblk) = spval
                 if (n_iage     (ns) /= 0) a2D(i,j,n_iage(ns),     iblk) = spval
                 if (n_FY       (ns) /= 0) a2D(i,j,n_FY(ns),       iblk) = spval

                 if (n_a11      (ns) /= 0) a2D(i,j,n_a11(ns),      iblk) = spval
                 if (n_a12      (ns) /= 0) a2D(i,j,n_a12(ns),      iblk) = spval
                 if (n_e11      (ns) /= 0) a2D(i,j,n_e11(ns),      iblk) = spval
                 if (n_e12      (ns) /= 0) a2D(i,j,n_e12(ns),      iblk) = spval
                 if (n_e22      (ns) /= 0) a2D(i,j,n_e22(ns),      iblk) = spval
                 if (n_s11      (ns) /= 0) a2D(i,j,n_s11(ns),      iblk) = spval
                 if (n_s12      (ns) /= 0) a2D(i,j,n_s12(ns),      iblk) = spval
                 if (n_s22      (ns) /= 0) a2D(i,j,n_s22(ns),      iblk) = spval
                 if (n_yieldstress11 (ns) /= 0) a2D(i,j,n_yieldstress11(ns),iblk) = spval
                 if (n_yieldstress12 (ns) /= 0) a2D(i,j,n_yieldstress12(ns),iblk) = spval
                 if (n_yieldstress22 (ns) /= 0) a2D(i,j,n_yieldstress22(ns),iblk) = spval
              else
                 if (n_divu     (ns) /= 0) a2D(i,j,n_divu(ns),iblk)      = &
                       divu (i,j,iblk)*avail_hist_fields(n_divu(ns))%cona
                 if (n_shear    (ns) /= 0) a2D(i,j,n_shear(ns),iblk)     = &
                       shear(i,j,iblk)*avail_hist_fields(n_shear(ns))%cona
                 if (n_sig1     (ns) /= 0) a2D(i,j,n_sig1(ns),iblk)      = &
                       sig1 (i,j,iblk)*avail_hist_fields(n_sig1(ns))%cona
                 if (n_sig2     (ns) /= 0) a2D(i,j,n_sig2(ns),iblk)      = &
                       sig2 (i,j,iblk)*avail_hist_fields(n_sig2(ns))%cona
                 if (n_sigP     (ns) /= 0) a2D(i,j,n_sigP(ns),iblk)      = &
                       sigP (i,j,iblk)*avail_hist_fields(n_sigP(ns))%cona      
                 if (n_sistreave(ns) /= 0) a2D(i,j,n_sistreave(ns),iblk) = &
                       p5*(sig1(i,j,iblk)+sig2(i,j,iblk))*avail_hist_fields(n_sistreave(ns))%cona
                 if (n_sistremax(ns) /= 0) a2D(i,j,n_sistremax(ns),iblk) = &
                       p5*(sig1(i,j,iblk)-sig2(i,j,iblk))*avail_hist_fields(n_sistremax(ns))%cona
                 if (n_mlt_onset(ns) /= 0) a2D(i,j,n_mlt_onset(ns),iblk) = &
                       mlt_onset(i,j,iblk)
                 if (n_frz_onset(ns) /= 0) a2D(i,j,n_frz_onset(ns),iblk) = &
                       frz_onset(i,j,iblk)
                 if (n_hisnap   (ns) /= 0) a2D(i,j,n_hisnap(ns),iblk)    = &
                       vice(i,j,iblk)
                 if (n_aisnap   (ns) /= 0) a2D(i,j,n_aisnap(ns),iblk)    = &
                       aice(i,j,iblk)

                 if (kdyn == 2) then  ! for EAP dynamics different time of output
                    if (n_trsig    (ns) /= 0) a2D(i,j,n_trsig(ns),iblk ) = &
                                        strength(i,j,iblk)
                 else
                    if (n_trsig    (ns) /= 0) a2D(i,j,n_trsig(ns),iblk ) = &
                                       p25*(stressp_1(i,j,iblk) &
                                          + stressp_2(i,j,iblk) &
                                          + stressp_3(i,j,iblk) &
                                          + stressp_4(i,j,iblk))
                 endif

                 if (n_iage     (ns) /= 0) a2D(i,j,n_iage(ns),iblk)  = &
                       trcr(i,j,nt_iage,iblk)*avail_hist_fields(n_iage(ns))%cona
                 if (n_FY       (ns) /= 0) a2D(i,j,n_FY(ns),iblk)  = &
                       trcr(i,j,nt_FY,iblk)*avail_hist_fields(n_FY(ns))%cona

                 if (n_a11     (ns) /= 0) a2D(i,j,n_a11(ns),iblk)      = &
                       a11 (i,j,iblk)*avail_hist_fields(n_a11(ns))%cona
                 if (n_a12     (ns) /= 0) a2D(i,j,n_a12(ns),iblk)      = &
                       a12 (i,j,iblk)*avail_hist_fields(n_a12(ns))%cona
                 if (n_e11     (ns) /= 0) a2D(i,j,n_e11(ns),iblk)      = &
                       e11 (i,j,iblk)*avail_hist_fields(n_e11(ns))%cona
                 if (n_e12     (ns) /= 0) a2D(i,j,n_e12(ns),iblk)      = &
                       e12 (i,j,iblk)*avail_hist_fields(n_e12(ns))%cona
                 if (n_e22     (ns) /= 0) a2D(i,j,n_e22(ns),iblk)      = &
                       e22 (i,j,iblk)*avail_hist_fields(n_e22(ns))%cona
                 if (n_s11     (ns) /= 0) a2D(i,j,n_s11(ns),iblk)      = &
                       s11 (i,j,iblk)*avail_hist_fields(n_s11(ns))%cona
                 if (n_s12     (ns) /= 0) a2D(i,j,n_s12(ns),iblk)      = &
                       s12 (i,j,iblk)*avail_hist_fields(n_s12(ns))%cona
                 if (n_s22     (ns) /= 0) a2D(i,j,n_s22(ns),iblk)      = &
                       s22 (i,j,iblk)*avail_hist_fields(n_s22(ns))%cona
                 if (n_yieldstress11     (ns) /= 0) a2D(i,j,n_yieldstress11(ns),iblk)      = &
                       yieldstress11 (i,j,iblk)*avail_hist_fields(n_yieldstress11(ns))%cona
                 if (n_yieldstress12     (ns) /= 0) a2D(i,j,n_yieldstress12(ns),iblk)      = &
                       yieldstress12 (i,j,iblk)*avail_hist_fields(n_yieldstress12(ns))%cona
                 if (n_yieldstress22     (ns) /= 0) a2D(i,j,n_yieldstress22(ns),iblk)      = &
                       yieldstress22 (i,j,iblk)*avail_hist_fields(n_yieldstress22(ns))%cona
              endif
           enddo                ! i
           enddo                ! j

        enddo                   ! iblk
        !$OMP END PARALLEL DO

        time_end(ns) = time/int(secday)
        time_end(ns) = real(time_end(ns),kind=real_kind)

      !---------------------------------------------------------------
      ! write file
      !---------------------------------------------------------------

        call ice_timer_start(timer_readwrite)  ! reading/writing
        call ice_write_hist (ns)
        call ice_timer_stop(timer_readwrite)  ! reading/writing

      !---------------------------------------------------------------
      ! reset to zero
      !------------------------------------------------------------
        if (write_ic) then
           if (allocated(a2D))  a2D (:,:,:,:)     = c0
           if (allocated(a3Dc)) a3Dc(:,:,:,:,:)   = c0
           if (allocated(a3Dz)) a3Dz(:,:,:,:,:)   = c0
           if (allocated(a3Db)) a3Db(:,:,:,:,:)   = c0
           if (allocated(a3Da)) a3Da(:,:,:,:,:)   = c0
           if (allocated(a3Df)) a3Df(:,:,:,:,:)   = c0
           if (allocated(a4Di)) a4Di(:,:,:,:,:,:) = c0
           if (allocated(a4Ds)) a4Ds(:,:,:,:,:,:) = c0
           if (allocated(a4Df)) a4Df(:,:,:,:,:,:) = c0
           avgct(:) = c0
           albcnt(:,:,:,:) = c0
           write_ic = .false.        ! write initial condition once at most
        else
           avgct(ns) = c0
           albcnt(:,:,:,ns) = c0
        endif
!        if (write_history(ns)) albcnt(:,:,:,ns) = c0

        do n = 1,n2D
           if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) a2D(:,:,n,:) = c0
        enddo
        do n = n2D + 1, n3Dccum   
           nn = n - n2D               
           if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) a3Dc(:,:,:,nn,:) = c0
        enddo
        do n = n3Dccum + 1, n3Dzcum
           nn = n - n3Dccum
           if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) a3Dz(:,:,:,nn,:) = c0
        enddo
        do n = n3Dzcum + 1, n3Dbcum
           nn = n - n3Dzcum
           if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) a3Db(:,:,:,nn,:) = c0
        enddo
        do n = n3Dbcum + 1, n3Dacum
           nn = n - n3Dbcum
           if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) a3Da(:,:,:,nn,:) = c0
        enddo
        do n = n3Dacum + 1, n3Dfcum
           nn = n - n3Dacum
           if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) a3Df(:,:,:,nn,:) = c0
        enddo
        do n = n3Dfcum + 1, n4Dicum
           nn = n - n3Dfcum
           if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) a4Di(:,:,:,:,nn,:) = c0
        enddo
        do n = n4Dicum + 1, n4Dscum
           nn = n - n4Dicum
           if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) a4Ds(:,:,:,:,nn,:) = c0
        enddo
        do n = n4Dscum + 1, n4Dfcum
           nn = n - n4Dscum
           if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) a4Df(:,:,:,:,nn,:) = c0
        enddo

      endif  ! write_history or write_ic
      enddo  ! nstreams

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         if (new_year) then
            do j=jlo,jhi
            do i=ilo,ihi
               ! reset NH Jan 1
               if (lmask_n(i,j,iblk)) mlt_onset(i,j,iblk) = c0
               ! reset SH Jan 1 
               if (lmask_s(i,j,iblk)) frz_onset(i,j,iblk) = c0
            enddo
            enddo
         endif                  ! new_year

         if ( (month .eq. 7) .and. new_month ) then 
            do j=jlo,jhi
            do i=ilo,ihi
               ! reset SH Jul 1
               if (lmask_s(i,j,iblk)) mlt_onset(i,j,iblk) = c0
               ! reset NH Jul 1
               if (lmask_n(i,j,iblk)) frz_onset(i,j,iblk) = c0
            enddo
            enddo
         endif                  ! 1st of July
      enddo                     ! iblk
      !$OMP END PARALLEL DO

      end subroutine accum_hist

!=======================================================================

      end module ice_history

!=======================================================================
