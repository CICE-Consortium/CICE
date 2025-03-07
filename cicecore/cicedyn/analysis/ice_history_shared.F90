!=======================================================================
!
! Output files: netCDF or binary data, Fortran unformatted dumps
!
! The following variables are currently hard-wired as snapshots
!   (instantaneous rather than time-averages):
!   divu, shear, vort, sig1, sig2, sigP, trsig, mlt_onset,
!   frz_onset, hisnap, aisnap
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
! 2012 Elizabeth Hunke split code from ice_history.F90

      module ice_history_shared

      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_domain_size, only: max_nstrm
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_diag
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted

      implicit none

      private
      public :: define_hist_field, accum_hist_field, icefields_nml, construct_filename

      integer (kind=int_kind), public :: history_precision

      logical (kind=log_kind), public :: &
         hist_avg(max_nstrm)  ! if true, write averaged data instead of snapshots

      character (len=char_len_long), public :: &
         history_file  , & ! output file for history
         incond_file       ! output file for snapshot initial conditions

      character (len=char_len_long), public :: &
         history_dir   , & ! directory name for history file
         incond_dir        ! directory for snapshot initial conditions

      character (len=char_len), public :: &
         version_name

      character (len=char_len), public :: &
         history_format      , & ! history format, cdf1, cdf2, cdf5, etc
         history_rearranger      ! history file rearranger, box or subset for pio

      character (len=char_len), public :: &
         hist_suffix(max_nstrm)  ! appended to history_file in filename

      integer (kind=int_kind), public :: &
         history_iotasks     , & ! iotasks, root, stride defines io pes for pio
         history_root        , & ! iotasks, root, stride defines io pes for pio
         history_stride      , & ! iotasks, root, stride defines io pes for pio
         history_deflate     , & ! compression level for hdf5/netcdf4
         history_chunksize(2) ! chunksize for hdf5/netcdf4

      !---------------------------------------------------------------
      ! Instructions for adding a field: (search for 'example')
      !     Here or in ice_history_[process].F90:
      ! (1) Add to frequency flags (f_<field>)
      ! (2) Add to namelist (here and also in ice_in)
      ! (3) Add to index list
      !     In init_hist (in ice_history.F90):
      ! (4) Add define_hist_field call with vname, vdesc, vunit,
      !     and vcomment, vcellmeas, and conversion factor if necessary.
      ! (5) Add flag to broadcast list
      ! (6) Add accum_hist_field call with appropriate variable
      !---------------------------------------------------------------

      type, public :: ice_hist_field
          character (len=16) :: vname     ! variable name
          character (len=16) :: vunit     ! variable units
          character (len=25) :: vcoord    ! variable coordinates
          character (len=16) :: vcellmeas ! variable cell measures
          character (len=55) :: vdesc     ! variable description
          character (len=55) :: vcomment  ! variable description
          real (kind=dbl_kind) :: cona    ! multiplicative conversion factor
          real (kind=dbl_kind) :: conb    ! additive conversion factor
          character (len=1) :: vhistfreq  ! frequency of history output
          integer (kind=int_kind) :: vhistfreq_n ! number of vhistfreq intervals
          logical (kind=log_kind) :: avg_ice_present ! only average where ice is present
          logical (kind=log_kind) :: mask_ice_free_points ! mask ice-free points
      end type

      integer (kind=int_kind), parameter, public :: &
         max_avail_hist_fields = 800      ! Max number of history fields

      integer (kind=int_kind), public :: &
         num_avail_hist_fields_tot  = 0, & ! Current, total number of defined fields
         num_avail_hist_fields_2D   = 0, & ! Number of 2D fields
         num_avail_hist_fields_3Dz  = 0, & ! Number of 3D fields (vertical)
         num_avail_hist_fields_3Dc  = 0, & ! Number of 3D fields (thickness categories)
         num_avail_hist_fields_3Db  = 0, & ! Number of 3D fields (vertical biology)
         num_avail_hist_fields_3Da  = 0, & ! Number of 3D fields (vertical), snow-biology
         num_avail_hist_fields_3Df  = 0, & ! Number of 3D fields (floe size categories)
         num_avail_hist_fields_4Di  = 0, & ! Number of 4D fields (categories,vertical), ice
         num_avail_hist_fields_4Ds  = 0, & ! Number of 4D fields (categories,vertical), snow
         num_avail_hist_fields_4Df  = 0    ! Number of 4D fields (floe size, thickness categories)

      integer (kind=int_kind), public :: &        ! cumulative counts
         n2D     , & ! num_avail_hist_fields_2D
         n3Dccum , & ! n2D     + num_avail_hist_fields_3Dc
         n3Dzcum , & ! n3Dccum + num_avail_hist_fields_3Dz
         n3Dbcum , & ! n3Dzcum + num_avail_hist_fields_3Db
         n3Dacum , & ! n3Dbcum + num_avail_hist_fields_3Da
         n3Dfcum , & ! n3Dacum + num_avail_hist_fields_3Df
         n4Dicum , & ! n3Dfcum + num_avail_hist_fields_4Di
         n4Dscum , & ! n4Dicum + num_avail_hist_fields_4Ds
         n4Dfcum     ! n4Dscum + num_avail_hist_fields_4Df

      ! could set nzilyr = nilyr + nslyr and write Tin+Tsn together into Tinz
      integer (kind=int_kind), public :: &
         nzilyr , & ! vertical dimension (allows alternative grids)
         nzslyr , & ! snow
         nzblyr , & ! bio grid
         nzalyr     ! aerosols (2 snow & nblyr+2 bio)

      type (ice_hist_field), public :: &
         avail_hist_fields(max_avail_hist_fields)

      integer (kind=int_kind), parameter, public :: &
         ncoord   = 8           , & ! number of coordinate variables: TLON, TLAT, ULON, ULAT, NLON, NLAT, ELON, ELAT
         nvar_grd = 21          , & ! number of grid fields that can be written
                                    !   excluding grid vertices
         nvar_grdz = 6              ! number of category/vertical grid fields written

      integer (kind=int_kind), public :: &
         ncat_hist              , & ! number of thickness categories written <= ncat
         nfsd_hist                  ! number of floe size categories written <= nfsd

      real (kind=real_kind), public :: time_beg(max_nstrm), & ! bounds for averaging
                                       time_end(max_nstrm), &
                                       time_bounds(2)

      character (len=char_len), public :: hist_time_axis

      real (kind=dbl_kind), allocatable, public :: &
         a2D (:,:,:,:)    , & ! field accumulations/averages, 2D
         a3Dz(:,:,:,:,:)  , & ! field accumulations/averages, 3D vertical
         a3Db(:,:,:,:,:)  , & ! field accumulations/averages, 3D vertical biology
         a3Dc(:,:,:,:,:)  , & ! field accumulations/averages, 3D thickness categories
         a3Da(:,:,:,:,:)  , & ! field accumulations/averages, 3D snow+bio
         a3Df(:,:,:,:,:)  , & ! field accumulations/averages, 3D floe size categories
         a4Di(:,:,:,:,:,:), & ! field accumulations/averages, 4D categories,vertical, ice
         a4Ds(:,:,:,:,:,:), & ! field accumulations/averages, 4D categories,vertical, snow
         a4Df(:,:,:,:,:,:)    ! field accumulations/averages, 4D floe size, thickness categories

      real (kind=dbl_kind), allocatable, public :: &
         Tinz4d (:,:,:,:)    , & ! array for Tin
         Tsnz4d (:,:,:,:)    , & ! array for Tsn
         Sinz4d (:,:,:,:)        ! array for Sin

      real (kind=dbl_kind), public :: &
         avgct(max_nstrm)   ! average sample counter

      logical (kind=log_kind), public :: &
         icoord(ncoord) , &    ! true if coord field is written to output file
         igrd (nvar_grd), &    ! true if grid field is written to output file
         igrdz(nvar_grdz)      ! true if category/vertical grid field is written

      character (len=25), public, parameter :: &
         ! T grids
         tcstr   = 'area: tarea'         , & ! vcellmeas for T cell quantities
         tstr2D  = 'TLON TLAT time'      , & ! vcoord for T cell, 2D
         tstr3Dc = 'TLON TLAT NCAT  time', & ! vcoord for T cell, 3D, ncat
         tstr3Da = 'TLON TLAT VGRDa time', & ! vcoord for T cell, 3D, ice-snow-bio
         tstr3Db = 'TLON TLAT VGRDb time', & ! vcoord for T cell, 3D, ice-bio
         tstr3Df = 'TLON TLAT NFSD  time', & ! vcoord for T cell, 3D, fsd
         tstr4Di = 'TLON TLAT VGRDi NCAT', & ! vcoord for T cell, 4D, ice
         tstr4Ds = 'TLON TLAT VGRDs NCAT', & ! vcoord for T cell, 4D, snow
         tstr4Df = 'TLON TLAT NFSD  NCAT', & ! vcoord for T cell, 4D, fsd
         ! U grids
         ucstr   = 'area: uarea'         , & ! vcellmeas for U cell quantities
         ustr2D  = 'ULON ULAT time'      , & ! vcoord for U cell, 2D
         ! N grids
         ncstr   = 'area: narea'         , & ! vcellmeas for N cell quantities
         nstr2D  = 'NLON NLAT time'      , & ! vcoord for N cell, 2D
         ! E grids
         ecstr   = 'area: earea'         , & ! vcellmeas for E cell quantities
         estr2D  = 'ELON ELAT time'          ! vcoord for E cell, 2D

      !---------------------------------------------------------------
      ! flags: write to output file if true or histfreq value
      !---------------------------------------------------------------

      logical (kind=log_kind), public :: &
           f_tlon      = .true., f_tlat       = .true., &
           f_ulon      = .true., f_ulat       = .true., &
           f_nlon      = .true., f_nlat       = .true., &
           f_elon      = .true., f_elat       = .true., &
           f_tmask     = .true., f_umask      = .true., &
           f_nmask     = .true., f_emask      = .true., &
           f_blkmask   = .true., &
           f_tarea     = .true., f_uarea      = .true., &
           f_narea     = .true., f_earea      = .true., &
           f_dxt       = .true., f_dyt        = .true., &
           f_dxu       = .true., f_dyu        = .true., &
           f_dxn       = .true., f_dyn        = .true., &
           f_dxe       = .true., f_dye        = .true., &
           f_HTN       = .true., f_HTE        = .true., &
           f_ANGLE     = .true., f_ANGLET     = .true., &
           f_bounds    = .true., f_NCAT       = .true., &
           f_VGRDi     = .true., f_VGRDs      = .true., &
           f_VGRDb     = .true., f_VGRDa      = .true., &
           f_NFSD      = .false.

      character (len=max_nstrm), public :: &
!          f_example   = 'md', &
           f_hi        = 'm', f_hs         = 'm', &
           f_snowfrac  = 'x', f_snowfracn  = 'x', &
           f_Tsfc      = 'm', f_aice       = 'm', &
           f_uvel      = 'm', f_vvel       = 'm', &
           f_icespd    = 'm', f_icedir     = 'm', &
           f_uvelE     = 'x', f_vvelE      = 'x', &
           f_icespdE   = 'x', f_icedirE    = 'x', &
           f_uvelN     = 'x', f_vvelN      = 'x', &
           f_icespdN   = 'x', f_icedirN    = 'x', &
           f_uatm      = 'm', f_vatm       = 'm', &
           f_atmspd    = 'm', f_atmdir     = 'm', &
           f_fswup     = 'm', &
           f_fswdn     = 'm', f_flwdn      = 'm', &
           f_snow      = 'm', f_snow_ai    = 'm', &
           f_rain      = 'm', f_rain_ai    = 'm', &
           f_sst       = 'm', f_sss        = 'm', &
           f_uocn      = 'm', f_vocn       = 'm', &
           f_ocnspd    = 'm', f_ocndir     = 'm', &
           f_sice      = 'm', f_frzmlt     = 'm', &
           f_fswfac    = 'm', f_fswint_ai  = 'x', &
           f_fswabs    = 'm', f_fswabs_ai  = 'm', &
           f_albsni    = 'm', &
           f_alvdr     = 'x', f_alidr      = 'x', &
           f_alvdf     = 'x', f_alidf      = 'x', &
           f_alvdr_ai  = 'm', f_alidr_ai   = 'm', &
           f_alvdf_ai  = 'm', f_alidf_ai   = 'm', &
           f_albice    = 'm', f_albsno     = 'm', &
           f_albpnd    = 'm', f_coszen     = 'm', &
           f_flat      = 'm', f_flat_ai    = 'm', &
           f_fsens     = 'm', f_fsens_ai   = 'm', &
           f_flwup     = 'm', f_flwup_ai   = 'm', &
           f_evap      = 'm', f_evap_ai    = 'm', &
           f_Tair      = 'm', &
           f_Tref      = 'm', f_Qref       = 'm', &
           f_congel    = 'm', f_frazil     = 'm', &
           f_snoice    = 'm', f_dsnow      = 'm', &
           f_meltt     = 'm', f_melts      = 'm', &
           f_meltb     = 'm', f_meltl      = 'm', &
           f_fresh     = 'm', f_fresh_ai   = 'm', &
           f_fsalt     = 'm', f_fsalt_ai   = 'm', &
           f_fbot      = 'm', &
           f_fhocn     = 'm', f_fhocn_ai   = 'm', &
           f_fswthru   = 'm', f_fswthru_ai = 'm', &
           f_strairx   = 'm', f_strairy    = 'm', &
           f_strtltx   = 'm', f_strtlty    = 'm', &
           f_strcorx   = 'm', f_strcory    = 'm', &
           f_strocnx   = 'm', f_strocny    = 'm', &
           f_strintx   = 'm', f_strinty    = 'm', &
           f_taubx     = 'm', f_tauby      = 'm', &
           f_strairxN  = 'x', f_strairyN   = 'x', &
           f_strtltxN  = 'x', f_strtltyN   = 'x', &
           f_strcorxN  = 'x', f_strcoryN   = 'x', &
           f_strocnxN  = 'x', f_strocnyN   = 'x', &
           f_strintxN  = 'x', f_strintyN   = 'x', &
           f_taubxN    = 'x', f_taubyN     = 'x', &
           f_strairxE  = 'x', f_strairyE   = 'x', &
           f_strtltxE  = 'x', f_strtltyE   = 'x', &
           f_strcorxE  = 'x', f_strcoryE   = 'x', &
           f_strocnxE  = 'x', f_strocnyE   = 'x', &
           f_strintxE  = 'x', f_strintyE   = 'x', &
           f_taubxE    = 'x', f_taubyE     = 'x', &
           f_strength  = 'm', f_vort       = 'm', &
           f_divu      = 'm', f_shear      = 'm', &
           f_sig1      = 'm', f_sig2       = 'm', &
           f_sigP      = 'm', &
           f_dvidtt    = 'm', f_dvidtd     = 'm', &
           f_daidtt    = 'm', f_daidtd     = 'm', &
           f_dagedtt   = 'm', f_dagedtd    = 'm', &
           f_mlt_onset = 'm', f_frz_onset  = 'm', &
           f_iage      = 'm', f_FY         = 'm', &
           f_hisnap    = 'm', f_aisnap     = 'm', &
           f_CMIP = 'x'     , &
           f_sithick   = 'x', f_sisnthick  = 'x', &
           f_siage     = 'x', &
           f_sitemptop = 'x', f_sitempsnic = 'x', &
           f_sitempbot = 'x', &
           f_sispeed   = 'x', f_sidir      = 'x', &
           f_siu       = 'x', f_siv        = 'x', &
           f_sidmasstranx = 'x', f_sidmasstrany = 'x', &
           f_sistrxdtop = 'x', f_sistrydtop = 'x', &
           f_sistrxubot = 'x', f_sistryubot = 'x', &
           f_sicompstren = 'x', &
           f_sialb     = 'x', &
           f_sihc      = 'x', f_sisnhc     = 'x', &
           f_sidconcth = 'x', f_sidconcdyn = 'x', &
           f_sidmassth = 'x', f_sidmassdyn = 'x', &
           f_sidmassgrowthwat = 'x', &
           f_sidmassgrowthbot = 'x', &
           f_sidmasssi = 'x', &
           f_sidmassevapsubl = 'x', &
           f_sndmasssubl = 'x', &
           f_sidmassmelttop = 'x', &
           f_sidmassmeltbot = 'x', &
           f_sidmasslat = 'x', &
           f_sndmasssnf = 'x', &
           f_sndmassmelt = 'x', &
           f_sndmassdyn = 'x', &
           f_siflswdtop = 'x', &
           f_siflswutop = 'x', &
           f_siflswdbot = 'x', &
           f_sifllwdtop = 'x', &
           f_sifllwutop = 'x', &
           f_siflsenstop = 'x', &
           f_siflsensupbot = 'x', &
           f_sifllatstop = 'x', &
           f_siflcondtop = 'x', &
           f_siflcondbot = 'x', &
           f_sipr = 'x', &
           f_sifb = 'x', &
           f_siflsaltbot = 'x', &
           f_siflfwbot = 'x', &
           f_siflfwdrain = 'x', &
           f_siforcetiltx = 'x', &
           f_siforcetilty = 'x', &
           f_siforcecoriolx = 'x', &
           f_siforcecorioly = 'x', &
           f_siforceintstrx = 'x', &
           f_siforceintstry = 'x', &
           f_siitdconc = 'x', &
           f_siitdthick = 'x', &
           f_siitdsnthick = 'x', &
           f_sidragtop = 'x', &
           f_sirdgthick = 'x', &
           f_sistreave = 'x', &
           f_sistremax = 'x', &
           f_aicen     = 'x', f_vicen      = 'x', &
           f_vsnon     = 'x', &
           f_trsig     = 'm', f_icepresent = 'm', &
           f_fsurf_ai  = 'm', f_fcondtop_ai= 'm', &
           f_fmeltt_ai = 'm',                     &
           f_fsurfn_ai = 'x' ,f_fcondtopn_ai='x', &
           f_fmelttn_ai= 'x', f_flatn_ai   = 'x', &
           f_fsensn_ai = 'x', &
!          f_field3dz  = 'x', &
           f_keffn_top = 'x', &
           f_Tinz      = 'x', f_Sinz       = 'x', &
           f_Tsnz      = 'x', &
           f_a11       = 'x', f_a12        = 'x', &
           f_e11       = 'x', f_e12        = 'x', &
           f_e22       = 'x', &
           f_s11       = 'x', f_s12        = 'x', &
           f_s22       = 'x', &
           f_yieldstress11  = 'x', &
           f_yieldstress12  = 'x', &
           f_yieldstress22  = 'x'

      !---------------------------------------------------------------
      ! namelist variables
      !---------------------------------------------------------------

      namelist / icefields_nml /     &
           f_tlon     , f_tlat     , &
           f_ulon     , f_ulat     , &
           f_nlon     , f_nlat     , &
           f_elon     , f_elat     , &
           f_tmask    , f_umask    , &
           f_nmask    , f_emask    , &
           f_blkmask  , &
           f_tarea    , f_uarea    , &
           f_narea    , f_earea    , &
           f_dxt      , f_dyt      , &
           f_dxu      , f_dyu      , &
           f_dxn      , f_dyn      , &
           f_dxe      , f_dye      , &
           f_HTN      , f_HTE      , &
           f_ANGLE    , f_ANGLET   , &
           f_bounds   , f_NCAT     , &
           f_VGRDi    , f_VGRDs    , &
           f_VGRDb    , f_VGRDa    , &
           f_NFSD     , &
!          f_example  , &
           f_hi,        f_hs       , &
           f_snowfrac,  f_snowfracn, &
           f_Tsfc,      f_aice     , &
           f_uvel,      f_vvel     , &
           f_icespd,    f_icedir   , &
!          For now, C and CD grid quantities are controlled by the generic (originally B-grid) namelist flag
!          f_uvelE,     f_vvelE    , &
!          f_icespdE,   f_icedirE  , &
!          f_uvelN,     f_vvelN    , &
!          f_icespdN,   f_icedirN  , &
           f_uatm,      f_vatm     , &
           f_atmspd,    f_atmdir   , &
           f_fswup,     &
           f_fswdn,     f_flwdn    , &
           f_snow,      f_snow_ai  , &
           f_rain,      f_rain_ai  , &
           f_sst,       f_sss      , &
           f_uocn,      f_vocn     , &
           f_ocnspd,    f_ocndir   , &
           f_sice,      f_frzmlt   , &
           f_fswfac,    f_fswint_ai, &
           f_fswabs,    f_fswabs_ai, &
           f_albsni,    &
           f_alvdr,     f_alidr    , &
           f_alvdf,     f_alidf    , &
           f_alvdr_ai,  f_alidr_ai , &
           f_alvdf_ai,  f_alidf_ai , &
           f_albice,    f_albsno   , &
           f_albpnd,    f_coszen   , &
           f_flat,      f_flat_ai  , &
           f_fsens,     f_fsens_ai , &
           f_flwup,     f_flwup_ai , &
           f_evap,      f_evap_ai  , &
           f_Tair,      &
           f_Tref,      f_Qref     , &
           f_congel,    f_frazil   , &
           f_snoice,    f_dsnow    , &
           f_meltt,     f_melts    , &
           f_meltb,     f_meltl    , &
           f_fresh,     f_fresh_ai , &
           f_fsalt,     f_fsalt_ai , &
           f_fbot,      &
           f_fhocn,     f_fhocn_ai , &
           f_fswthru,   f_fswthru_ai,&
           f_strairx,   f_strairy  , &
           f_strtltx,   f_strtlty  , &
           f_strcorx,   f_strcory  , &
           f_strocnx,   f_strocny  , &
           f_strintx,   f_strinty  , &
           f_taubx,     f_tauby    , &
!          f_strairxN,  f_strairyN , &
!          f_strtltxN,  f_strtltyN , &
!          f_strcorxN,  f_strcoryN , &
!          f_strocnxN,  f_strocnyN , &
!          f_strintxN,  f_strintyN , &
!          f_taubxN,    f_taubyN   , &
!          f_strairxE,  f_strairyE , &
!          f_strtltxE,  f_strtltyE , &
!          f_strcorxE,  f_strcoryE , &
!          f_strocnxE,  f_strocnyE , &
!          f_strintxE,  f_strintyE , &
!          f_taubxE,    f_taubyE   , &
           f_strength,  f_vort     , &
           f_divu,      f_shear    , &
           f_sig1,      f_sig2     , &
           f_sigP,      &
           f_dvidtt,    f_dvidtd   , &
           f_daidtt,    f_daidtd   , &
           f_dagedtt,   f_dagedtd  , &
           f_mlt_onset, f_frz_onset, &
           f_iage,      f_FY       , &
           f_hisnap,    f_aisnap   , &
           f_CMIP, &
           f_sithick,   f_sisnthick, &
           f_siage,     &
           f_sitemptop, f_sitempsnic,&
           f_sitempbot, &
           f_sispeed,   f_sidir,     &
           f_siu,       f_siv,       &
           f_sidmasstranx, f_sidmasstrany, &
           f_sistrxdtop, f_sistrydtop, &
           f_sistrxubot, f_sistryubot, &
           f_sicompstren, &
           f_sialb, &
           f_sihc,      f_sisnhc,    &
           f_sidconcth, f_sidconcdyn,&
           f_sidmassth, f_sidmassdyn,&
           f_sidmassgrowthwat, &
           f_sidmassgrowthbot, &
           f_sidmasssi, &
           f_sidmassevapsubl, &
           f_sndmasssubl, &
           f_sidmassmelttop, &
           f_sidmassmeltbot, &
           f_sidmasslat, &
           f_sndmasssnf, &
           f_sndmassmelt, &
           f_sndmassdyn, &
           f_siflswdtop, &
           f_siflswutop, &
           f_siflswdbot, &
           f_sifllwdtop, &
           f_sifllwutop, &
           f_siflsenstop, &
           f_siflsensupbot, &
           f_sifllatstop, &
           f_siflcondtop, &
           f_siflcondbot, &
           f_sipr, &
           f_sifb, &
           f_siflsaltbot, &
           f_siflfwbot, &
           f_siflfwdrain, &
           f_siforcetiltx, &
           f_siforcetilty, &
           f_siforcecoriolx, &
           f_siforcecorioly, &
           f_siforceintstrx, &
           f_siforceintstry, &
           f_siitdconc, &
           f_siitdthick, &
           f_siitdsnthick, &
           f_sidragtop, &
           f_sirdgthick, &
           f_sistreave, &
           f_sistremax, &
           f_aicen,     f_vicen    , &
           f_vsnon,     &
           f_trsig,     f_icepresent,&
           f_fsurf_ai,  f_fcondtop_ai,&
           f_fmeltt_ai, &
           f_fsurfn_ai,f_fcondtopn_ai,&
           f_fmelttn_ai,f_flatn_ai,  &
           f_fsensn_ai, &
!          f_field3dz,  &
           f_keffn_top, &
           f_Tinz,      f_Sinz,      &
           f_Tsnz,      &
           f_a11,       f_a12,       &
           f_e11,       f_e12,       &
           f_e22,       &
           f_s11,       f_s12,       &
           f_s22,       &
           f_yieldstress11, &
           f_yieldstress12, &
           f_yieldstress22

      !---------------------------------------------------------------
      ! field indices
      !---------------------------------------------------------------

      integer (kind=int_kind), parameter, public :: &
           n_tlon       = 1,  &
           n_tlat       = 2,  &
           n_ulon       = 3,  &
           n_ulat       = 4,  &
           n_nlon       = 5,  &
           n_nlat       = 6,  &
           n_elon       = 7,  &
           n_elat       = 8,  &

           n_tmask      = 1,  &
           n_umask      = 2,  &
           n_nmask      = 3,  &
           n_emask      = 4,  &
           n_blkmask    = 5,  &
           n_tarea      = 6,  &
           n_uarea      = 7,  &
           n_narea      = 8,  &
           n_earea      = 9,  &
           n_dxt        = 10, &
           n_dyt        = 11, &
           n_dxu        = 12, &
           n_dyu        = 13, &
           n_dxn        = 14, &
           n_dyn        = 15, &
           n_dxe        = 16, &
           n_dye        = 17, &
           n_HTN        = 18, &
           n_HTE        = 19, &
           n_ANGLE      = 20, &
           n_ANGLET     = 21, &

           n_NCAT       = 1, &
           n_VGRDi      = 2, &
           n_VGRDs      = 3, &
           n_VGRDb      = 4, &
           n_VGRDa      = 5, &
           n_NFSD       = 6, &

           n_lont_bnds  = 1, &
           n_latt_bnds  = 2, &
           n_lonu_bnds  = 3, &
           n_latu_bnds  = 4, &
           n_lonn_bnds  = 5, &
           n_latn_bnds  = 6, &
           n_lone_bnds  = 7, &
           n_late_bnds  = 8

      integer (kind=int_kind), dimension(max_nstrm), public :: &
!          n_example    , &
           n_hi         , n_hs         , &
           n_snowfrac   , n_snowfracn  , &
           n_Tsfc       , n_aice       , &
           n_uvel       , n_vvel       , &
           n_icespd     , n_icedir     , &
           n_uvelE      , n_vvelE      , &
           n_icespdE    , n_icedirE    , &
           n_uvelN      , n_vvelN      , &
           n_icespdN    , n_icedirN    , &
           n_uatm       , n_vatm       , &
           n_atmspd     , n_atmdir     , &
           n_sice       , &
           n_fswup      , &
           n_fswdn      , n_flwdn      , &
           n_snow       , n_snow_ai    , &
           n_rain       , n_rain_ai    , &
           n_sst        , n_sss        , &
           n_uocn       , n_vocn       , &
           n_ocnspd     , n_ocndir     , &
           n_frzmlt     , n_fswfac     , &
           n_fswint_ai  , &
           n_fswabs     , n_fswabs_ai  , &
           n_albsni     , &
           n_alvdr      , n_alidr      , &
           n_alvdf      , n_alidf      , &
           n_alvdr_ai   , n_alidr_ai   , &
           n_alvdf_ai   , n_alidf_ai   , &
           n_albice     , n_albsno     , &
           n_albpnd     , n_coszen     , &
           n_flat       , n_flat_ai    , &
           n_fsens      , n_fsens_ai   , &
           n_flwup      , n_flwup_ai   , &
           n_evap       , n_evap_ai    , &
           n_Tair       , &
           n_Tref       , n_Qref       , &
           n_congel     , n_frazil     , &
           n_snoice     , n_dsnow      , &
           n_meltt      , n_melts      , &
           n_meltb      , n_meltl      , &
           n_fresh      , n_fresh_ai   , &
           n_fsalt      , n_fsalt_ai   , &
           n_vsnon      , &
           n_fbot       , &
           n_fhocn      , n_fhocn_ai   , &
           n_fswthru    , n_fswthru_ai , &
           n_strairx    , n_strairy    , &
           n_strtltx    , n_strtlty    , &
           n_strcorx    , n_strcory    , &
           n_strocnx    , n_strocny    , &
           n_strintx    , n_strinty    , &
           n_taubx      , n_tauby      , &
           n_strairxN   , n_strairyN   , &
           n_strtltxN   , n_strtltyN   , &
           n_strcorxN   , n_strcoryN   , &
           n_strocnxN   , n_strocnyN   , &
           n_strintxN   , n_strintyN   , &
           n_taubxN     , n_taubyN     , &
           n_strairxE   , n_strairyE   , &
           n_strtltxE   , n_strtltyE   , &
           n_strcorxE   , n_strcoryE   , &
           n_strocnxE   , n_strocnyE   , &
           n_strintxE   , n_strintyE   , &
           n_taubxE     , n_taubyE     , &
           n_strength   , n_vort       , &
           n_divu       , n_shear      , &
           n_sig1       , n_sig2       , &
           n_sigP       , &
           n_dvidtt     , n_dvidtd     , &
           n_daidtt     , n_daidtd     , &
           n_dagedtt    , n_dagedtd    , &
           n_mlt_onset  , n_frz_onset  , &
           n_hisnap     , n_aisnap     , &
           n_sithick    , n_sisnthick  , &
           n_siage,       &
           n_sitemptop  , n_sitempsnic , &
           n_sitempbot  , &
           n_sispeed    , n_sidir      , &
           n_siu,         n_siv,         &
           n_sidmasstranx, n_sidmasstrany, &
           n_sistrxdtop,  n_sistrydtop,  &
           n_sistrxubot,  n_sistryubot,  &
           n_sicompstren, &
           n_sialb, &
           n_sihc       , n_sisnhc,      &
           n_sidconcth  , n_sidconcdyn,  &
           n_sidmassth  , n_sidmassdyn,  &
           n_sidmassgrowthwat,  &
           n_sidmassgrowthbot,  &
           n_sidmasssi,  &
           n_sidmassevapsubl,  &
           n_sndmasssubl,  &
           n_sidmassmelttop,  &
           n_sidmassmeltbot,  &
           n_sidmasslat,  &
           n_sndmasssnf,  &
           n_sndmassmelt,  &
           n_sndmassdyn,  &
           n_siflswdtop,  &
           n_siflswutop,  &
           n_siflswdbot,  &
           n_sifllwdtop,  &
           n_sifllwutop,  &
           n_siflsenstop,  &
           n_siflsensupbot,  &
           n_sifllatstop,  &
           n_siflcondtop,  &
           n_siflcondbot,  &
           n_sipr,  &
           n_sifb,  &
           n_siflsaltbot,  &
           n_siflfwbot,  &
           n_siflfwdrain,  &
           n_siforcetiltx,  &
           n_siforcetilty,  &
           n_siforcecoriolx,  &
           n_siforcecorioly,  &
           n_siforceintstrx,  &
           n_siforceintstry,  &
           n_siitdconc, &
           n_siitdthick, &
           n_siitdsnthick, &
           n_sidragtop, &
           n_sirdgthick, &
           n_sistreave, &
           n_sistremax, &
           n_trsig      , n_icepresent , &
           n_iage       , n_FY         , &
           n_fsurf_ai   , &
           n_fcondtop_ai, n_fmeltt_ai  , &
           n_aicen      , n_vicen      , &
           n_fsurfn_ai   , &
           n_fcondtopn_ai, &
           n_fmelttn_ai  , &
           n_flatn_ai    , &
           n_fsensn_ai   , &
!          n_field3dz    , &
           n_keffn_top   , &
           n_Tinz        , n_Sinz      , &
           n_Tsnz        , &
           n_a11         , n_a12       , &
           n_e11         , n_e12       , &
           n_e22         , &
           n_s11         , n_s12       , &
           n_s22         , &
           n_yieldstress11, n_yieldstress12, &
           n_yieldstress22

      interface accum_hist_field ! generic interface
           module procedure accum_hist_field_2D, &
                            accum_hist_field_3D, &
                            accum_hist_field_4D
      end interface

!=======================================================================

      contains

!=======================================================================

      subroutine construct_filename(ncfile,suffix,ns)

      use ice_calendar, only: msec, myear, mmonth, daymo,  &
                              mday, write_ic, histfreq, histfreq_n, &
                              new_year, new_month, new_day, &
                              dt
      use ice_restart_shared, only: lenstr

      character (len=*), intent(inout) :: ncfile
      character (len=*), intent(in) :: suffix
      integer (kind=int_kind), intent(in) :: ns

      integer (kind=int_kind) :: iyear, imonth, iday, isec
      integer (kind=int_kind) :: n
      character (len=char_len) :: cstream
      character (len=char_len_long), save :: ncfile_last(max_nstrm) = 'UnDefineD'
      character(len=*), parameter :: subname = '(construct_filename)'

        iyear = myear
        imonth = mmonth
        iday = mday
        isec = int(msec - dt,int_kind)
        cstream = ''
        if (hist_suffix(ns) /= 'x') cstream = hist_suffix(ns)

        ! construct filename
        if (write_ic) then
           isec = msec
           write(ncfile,'(a,a,i4.4,a,i2.2,a,i2.2,a,i5.5,a,a)')  &
              incond_file(1:lenstr(incond_file)),'.',iyear,'-', &
              imonth,'-',iday,'-',isec,'.',trim(suffix)
        else

           if (hist_avg(ns)) then
              if (histfreq(ns) == '1' .or. histfreq(ns) == 'h'.or.histfreq(ns) == 'H') then
                 ! do nothing
              elseif (new_year) then
                 iyear = iyear - 1
                 imonth = 12
                 iday = daymo(imonth)
              elseif (new_month) then
                 imonth = mmonth - 1
                 iday = daymo(imonth)
              elseif (new_day) then
                 iday = iday - 1
              endif
           endif

           if (hist_avg(ns)) then    ! write averaged data
              if (histfreq(ns) == '1' .and. histfreq_n(ns) == 1)  then ! timestep
                 write(ncfile,'(a,a,i4.4,a,i2.2,a,i2.2,a,i5.5,a,a)')  &
                       history_file(1:lenstr(history_file))//trim(cstream),'_inst.', &
                       iyear,'-',imonth,'-',iday,'-',msec,'.',trim(suffix)
              elseif (histfreq(ns) == '1' .and. histfreq_n(ns) > 1)  then ! timestep
                 write(ncfile,'(a,a,i4.4,a,i2.2,a,i2.2,a,i5.5,a,a)')  &
                       history_file(1:lenstr(history_file))//trim(cstream),'.', &
                       iyear,'-',imonth,'-',iday,'-',msec,'.',trim(suffix)
              elseif (histfreq(ns) == 'h'.or.histfreq(ns) == 'H') then ! hourly
                 write(ncfile,'(a,a,i2.2,a,i4.4,a,i2.2,a,i2.2,a,i5.5,a,a)')  &
                       history_file(1:lenstr(history_file))//trim(cstream),'_', &
                       histfreq_n(ns),'h.',iyear,'-',imonth,'-',iday,'-',msec,'.',trim(suffix)
              elseif (histfreq(ns) == 'd'.or.histfreq(ns) == 'D') then ! daily
                 write(ncfile,'(a,a,i4.4,a,i2.2,a,i2.2,a,a)')  &
                       history_file(1:lenstr(history_file))//trim(cstream),'.', &
                       iyear,'-',imonth,'-',iday,'.',trim(suffix)
              elseif (histfreq(ns) == 'm'.or.histfreq(ns) == 'M') then ! monthly
                 write(ncfile,'(a,a,i4.4,a,i2.2,a,a)')  &
                       history_file(1:lenstr(history_file))//trim(cstream),'.', &
                       iyear,'-',imonth,'.',trim(suffix)
              elseif (histfreq(ns) == 'y'.or.histfreq(ns) == 'Y') then ! yearly
                 write(ncfile,'(a,a,i4.4,a,a)') &
                       history_file(1:lenstr(history_file))//trim(cstream),'.', &
                       iyear,'.',trim(suffix)
              elseif (histfreq(ns) == 'g') then
                 write(ncfile,'(a,a,a,a)')  &
                       history_file(1:lenstr(history_file)),'_grid', &
                       '.',trim(suffix)
              endif

           else                     ! instantaneous
              if (histfreq(ns) == 'g') then
                 write(ncfile,'(a,a,a,a)')  &
                       history_file(1:lenstr(history_file)),'_grid', &
                       '.',trim(suffix)
              else
                 write(ncfile,'(a,a,i4.4,a,i2.2,a,i2.2,a,i5.5,a,a)')  &
                       history_file(1:lenstr(history_file))//trim(cstream),'_inst.', &
                       iyear,'-',imonth,'-',iday,'-',msec,'.',trim(suffix)
              endif
           endif

        endif

        ! Check whether the filename is already in use.
        ! Same filename in multiple history streams leads to files being overwritten (not good).
        ! The current filename convention means we just have to check latest filename,
        ! not all filenames ever generated because of use of current model date/time in filename.

        ! write(nu_diag,'(2a,i2,1x,a)') subname, 'debug ncfile= ',ns,trim(ncfile)
        do n = 1,max_nstrm
           ! write(nu_diag,'(2a,i2,1x,a)') subname, 'debug nfile_last= ',n,trim(ncfile_last(n))
           if (ncfile == ncfile_last(n)) then
              write(nu_diag,*) subname,' history stream = ',ns
              write(nu_diag,*) subname,' history filename = ',trim(ncfile)
              write(nu_diag,*) subname,' filename in use for stream ',n
              write(nu_diag,*) subname,' filename for stream ',trim(ncfile_last(n))
              write(nu_diag,*) subname,' Use namelist hist_suffix so history filenames are unique'
              call abort_ice(subname//' ERROR: history filename already used for another history stream '//trim(ncfile))
           endif
        enddo
        ncfile_last(ns) = ncfile

      end subroutine construct_filename

!=======================================================================

!     Initializes description of an available field and returns location
!     in the available fields array for use in later calls.
!
!     2009 Created by D. Bailey following POP

      subroutine define_hist_field(id, vname, vunit, vcoord, vcellmeas, &
                                   vdesc, vcomment, cona, conb, &
                                   ns, vhistfreq, avg_ice_present, mask_ice_free_points)

      use ice_calendar, only: histfreq, histfreq_n

      integer (int_kind), dimension(:), intent(out) :: &  ! max_nstrm
         id                ! location in avail_fields array for use in
                           ! later routines

      character (len=*), intent(in) :: &
         vname      , & ! variable names
         vunit      , & ! variable units
         vcoord     , & ! variable coordinates
         vcellmeas  , & ! variables cell measures
         vdesc      , & ! variable descriptions
         vcomment       ! variable comments

      real (kind=dbl_kind), intent(in) :: &
         cona       , & ! multiplicative conversion factor
         conb           ! additive conversion factor

      character (len=*), intent(in) :: &
         vhistfreq      ! history frequency

      integer (kind=int_kind), intent(in) :: &
         ns             ! history file stream index

      logical (kind=log_kind), optional, intent(in) :: &
         avg_ice_present       , & ! compute average only when ice is present
         mask_ice_free_points      ! mask ice-free points

      integer (kind=int_kind) :: &
         ns1        , & ! variable stream loop index
         lenf           ! length of namelist string

      character (len=40) :: stmp

      logical (kind=log_kind) :: &
         l_avg_ice_present       , & ! compute average only when ice is present
         l_mask_ice_free_points      ! mask ice-free points

      character(len=*), parameter :: subname = '(define_hist_field)'

      l_avg_ice_present = .false.
      l_mask_ice_free_points = .false.

      if(present(avg_ice_present)) l_avg_ice_present = avg_ice_present
      if(present(mask_ice_free_points)) l_mask_ice_free_points = mask_ice_free_points

      if (histfreq(ns) == 'x') then
         call abort_ice(subname//' ERROR: define_hist_fields has histfreq x')
      endif

      if (ns == 1) id(:) = 0
      lenf = len(trim(vhistfreq))

      do ns1 = 1, lenf
         if (vhistfreq(ns1:ns1) == histfreq(ns)) then

            if (ns1 > 1 .and. index(vhistfreq(1:ns1-1),'x') /= 0) then
               call abort_ice(subname//' ERROR: history frequency variable f_' // vname // ' can''t contain ''x'' along with active frequencies')
            endif

            num_avail_hist_fields_tot = num_avail_hist_fields_tot + 1

            if (vcoord(11:14) == 'time') then
               num_avail_hist_fields_2D  = num_avail_hist_fields_2D + 1
            elseif (vcoord(11:14) == 'NCAT' .and. vcoord(17:20) == 'time') then
               num_avail_hist_fields_3Dc = num_avail_hist_fields_3Dc + 1
            elseif (vcoord(11:14) == 'NFSD' .and. vcoord(17:20) == 'time') then
               num_avail_hist_fields_3Df = num_avail_hist_fields_3Df + 1
            elseif (vcoord(11:15) == 'VGRDi' .and. vcoord(17:20) == 'time') then
               num_avail_hist_fields_3Dz = num_avail_hist_fields_3Dz + 1
            elseif (vcoord(11:15) == 'VGRDb' .and. vcoord(17:20) == 'time') then
               num_avail_hist_fields_3Db = num_avail_hist_fields_3Db + 1
            elseif (vcoord(11:15) == 'VGRDa' .and. vcoord(17:20) == 'time') then
               num_avail_hist_fields_3Da = num_avail_hist_fields_3Da + 1
            elseif (vcoord(11:15) == 'VGRDi' .and. vcoord(17:20) == 'NCAT') then
               num_avail_hist_fields_4Di = num_avail_hist_fields_4Di + 1
            elseif (vcoord(11:15) == 'VGRDs' .and. vcoord(17:20) == 'NCAT') then
               num_avail_hist_fields_4Ds = num_avail_hist_fields_4Ds + 1
            elseif (vcoord(11:14) == 'NFSD' .and. vcoord(17:20) == 'NCAT') then
               num_avail_hist_fields_4Df = num_avail_hist_fields_4Df + 1
            endif

            if (num_avail_hist_fields_tot > max_avail_hist_fields) then
               if (my_task == master_task) then
                  write(nu_diag,*) subname,' num_avail_hist_fields_tot = ',num_avail_hist_fields_tot
                  write(nu_diag,*) subname,' max_avail_hist_fields     = ',max_avail_hist_fields
               endif
               call abort_ice(subname//' ERROR: Need in computation of max_avail_hist_fields')
            endif

            if (num_avail_hist_fields_tot /= &
                num_avail_hist_fields_2D  + &
                num_avail_hist_fields_3Dc + &
                num_avail_hist_fields_3Dz + &
                num_avail_hist_fields_3Db + &
                num_avail_hist_fields_3Da + &
                num_avail_hist_fields_3Df + &
                num_avail_hist_fields_4Di + &
                num_avail_hist_fields_4Ds + &
                num_avail_hist_fields_4Df) then
               if (my_task == master_task) then
                  write(nu_diag,*) subname,' num_avail_hist_fields_tot = ',num_avail_hist_fields_tot
               endif
               call abort_ice(subname//' ERROR: in num_avail_hist_fields')
            endif

            id(ns) = num_avail_hist_fields_tot

            stmp = vname
            if (ns > 1) &
               write(stmp,'(a,a1,a1)') trim(stmp),'_',vhistfreq(ns1:ns1)

            avail_hist_fields(id(ns))%vname = trim(stmp)
            avail_hist_fields(id(ns))%vunit = trim(vunit)
            avail_hist_fields(id(ns))%vcoord = trim(vcoord)
            avail_hist_fields(id(ns))%vcellmeas = trim(vcellmeas)
            avail_hist_fields(id(ns))%vdesc = trim(vdesc)
            avail_hist_fields(id(ns))%vcomment = trim(vcomment)
            avail_hist_fields(id(ns))%cona = cona
            avail_hist_fields(id(ns))%conb = conb
            avail_hist_fields(id(ns))%vhistfreq = vhistfreq(ns1:ns1)
            avail_hist_fields(id(ns))%vhistfreq_n = histfreq_n(ns)
            avail_hist_fields(id(ns))%avg_ice_present = l_avg_ice_present
            avail_hist_fields(id(ns))%mask_ice_free_points = l_mask_ice_free_points

         endif
      enddo

      end subroutine define_hist_field

!=======================================================================

!     Accumulates a history field
!
!     2009 Created by D. Bailey following POP
!     2010 Generalized dimension of variables by N. Jeffery, E. Hunke

      subroutine accum_hist_field_2D(id, iblk, field_accum, field)

      use ice_blocks, only: block, get_block
      use ice_calendar, only: nstreams
      use ice_domain, only: blocks_ice
      use ice_grid, only: tmask

      integer (int_kind), dimension(:), intent(in) :: &  ! max_nstrm
         id                ! location in avail_fields array for use in
                           ! later routines

      integer (kind=int_kind), intent(in) :: iblk

      real (kind=dbl_kind), intent(in) :: &
         field_accum(:,:)

      real (kind=dbl_kind), intent(inout) :: &
         field(:,:,:,:)

      type (block) :: &
         this_block           ! block information for current block

      integer (kind=int_kind) :: i,j, ilo, ihi, jlo, jhi, ns, idns

      character(len=*), parameter :: subname = '(accum_hist_field_2D)'

      !---------------------------------------------------------------
      ! increment field
      !---------------------------------------------------------------

       do ns = 1, nstreams
       idns = id(ns)
       if (idns > 0) then

       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo, jhi
       do i = ilo, ihi
          if (tmask(i,j,iblk)) then
             field(i,j,idns, iblk) = field(i,j,idns, iblk) + field_accum(i,j)
          endif
       enddo
       enddo

       endif
       enddo

      end subroutine accum_hist_field_2D

!=======================================================================

!     Accumulates a history field
!
!     2009 Created by D. Bailey following POP
!     2010 Generalized dimension of variables by N. Jeffery, E. Hunke

      subroutine accum_hist_field_3D(id, iblk, ndim, field_accum, field)

      use ice_blocks, only: block, get_block
      use ice_calendar, only: nstreams
      use ice_domain, only: blocks_ice
      use ice_grid, only: tmask

      integer (int_kind), dimension(:), intent(in) :: &  ! max_nstrm
         id                ! location in avail_fields array for use in
                           ! later routines

      integer (kind=int_kind), intent(in) :: iblk

      integer (kind=int_kind), intent(in) :: &
         ndim              ! third dimension size

      real (kind=dbl_kind), intent(in) :: &
         field_accum(:,:,:)

      real (kind=dbl_kind), intent(inout) :: &
         field(:,:,:,:,:)

      type (block) :: &
         this_block           ! block information for current block

      integer (kind=int_kind) :: i,j,k, ilo, ihi, jlo, jhi, ns, idns

      character(len=*), parameter :: subname = '(accum_hist_field_3D)'

      !---------------------------------------------------------------
      ! increment field
      !---------------------------------------------------------------

       do ns = 1, nstreams
       idns = id(ns)
       if (idns > 0) then

       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do k = 1, ndim
       do j = jlo, jhi
       do i = ilo, ihi
          if (tmask(i,j,iblk)) then
             field(i,j,k,idns,iblk) = field(i,j,k,idns,iblk) + field_accum(i,j,k)
          endif
       enddo
       enddo
       enddo

       endif
       enddo

      end subroutine accum_hist_field_3D

!=======================================================================

!     Accumulates a history field
!
!     2009 Created by D. Bailey following POP
!     2010 Generalized dimension of variables by N. Jeffery, E. Hunke

      subroutine accum_hist_field_4D(id, iblk, ndim3, ndim4, field_accum, field)

      use ice_blocks, only: block, get_block
      use ice_calendar, only: nstreams
      use ice_domain, only: blocks_ice
      use ice_grid, only: tmask

      integer (int_kind), dimension(:), intent(in) :: &  ! max_nstrm
         id                ! location in avail_fields array for use in
                           ! later routines

      integer (kind=int_kind), intent(in) :: iblk

      integer (kind=int_kind), intent(in) :: &
         ndim3  , &        ! third dimension size
         ndim4             ! fourth dimension size

      real (kind=dbl_kind), intent(in) :: &
         field_accum(:,:,:,:)

      real (kind=dbl_kind), intent(inout) :: &
         field(:,:,:,:,:,:)

      type (block) :: &
         this_block           ! block information for current block

      integer (kind=int_kind) :: i,j,k,n,ilo, ihi, jlo, jhi, ns, idns

      character(len=*), parameter :: subname = '(accum_hist_field_4D)'

      !---------------------------------------------------------------
      ! increment field
      !---------------------------------------------------------------

       do ns = 1, nstreams
       idns = id(ns)
       if (idns > 0) then

       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do k = 1, ndim4
       do n = 1, ndim3
       do j = jlo, jhi
       do i = ilo, ihi
          if (tmask(i,j,iblk)) then
             field(i,j,n,k,idns,iblk) = field(i,j,n,k,idns,iblk) + field_accum(i,j,n,k)
          endif
       enddo
       enddo
       enddo
       enddo

       endif
       enddo

      end subroutine accum_hist_field_4D

!=======================================================================

      end module ice_history_shared

!=======================================================================
