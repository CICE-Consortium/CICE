!=======================================================================
! Biogeochemistry history output
!
! authors Elizabeth C. Hunke and Nicole Jeffery, LANL
!
! 2012 Elizabeth Hunke split code from ice_history.F90

      module ice_history_bgc

      use ice_kinds_mod
      use ice_constants
      use ice_fileunits, only: nu_nml, nml_filename, &
          get_fileunit, release_fileunit
      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_max_aero, icepack_max_dic, &
          icepack_max_doc, icepack_max_don, &
          icepack_max_algae, icepack_max_fe
      use icepack_intfc, only: icepack_query_tracer_flags, &
          icepack_query_tracer_indices, icepack_query_parameters, &
          icepack_query_parameters
      use ice_domain_size, only: max_nstrm, n_aero, &
          n_algae, n_dic, n_doc, n_don, n_zaero, n_fed, n_fep 

      implicit none
      private
      public :: init_hist_bgc_2D, init_hist_bgc_3Dc, &
                init_hist_bgc_3Db, init_hist_bgc_3Da,&
                accum_hist_bgc, init_history_bgc
      
      !---------------------------------------------------------------
      ! flags: write to output file if true or histfreq value
      !---------------------------------------------------------------
      ! specified in input_templates
      !--------------------------------------------------------------
      character (len=max_nstrm), public :: & 
           f_faero_atm    = 'x', f_faero_ocn    = 'x', &
           f_aero         = 'x', &
           f_fzsal        = 'm', f_fzsal_ai     = 'm', & 
           f_fzsal_g      = 'm', f_fzsal_g_ai   = 'm', &
           f_zsal         = 'x', &
           f_fbio         = 'x', f_fbio_ai      = 'x', &    
           f_zaero        = 'x', f_bgc_S        = 'x', & 
           f_bgc_N        = 'x', f_bgc_C        = 'x', &
           f_bgc_DOC      = 'x', f_bgc_DIC      = 'x', &
           f_bgc_chl      = 'x', f_bgc_Nit      = 'x', &
           f_bgc_Am       = 'x', f_bgc_Sil      = 'x', &
           f_bgc_DMSPp    = 'x', f_bgc_DMSPd    = 'x', &
           f_bgc_DMS      = 'x', f_bgc_DON      = 'x', & 
           f_bgc_Fe       = 'x', f_bgc_hum      = 'x', &
           f_bgc_PON      = 'x', f_bgc_ml       = 'x', &
           f_upNO         = 'x', f_upNH         = 'x', &  
           f_bTin         = 'x', f_bphi         = 'x', &  
           f_iDi          = 'x', f_iki          = 'x', &  
           f_fbri         = 'x', f_hbri         = 'x', &
           f_zfswin       = 'x', f_grownet      = 'x', &             
           f_bionet       = 'x', f_biosnow      = 'x', &  
           f_PPnet        = 'x', f_algalpeak    = 'x', &
           f_zbgc_frac    = 'x', &
         !------------------------------------------------
         ! specified by combinations of above values
         !-------------------------------------------------
           f_bgc_Fed      = 'x', f_bgc_Fep      = 'x', &
           f_DONnet       = 'x',                       &  
           f_DICnet       = 'x', f_DOCnet       = 'x', &  
           f_chlnet       = 'x', f_Nitnet       = 'x', &   
           f_Amnet        = 'x', f_Cnet         = 'x', &  
           f_Nnet         = 'x', f_DMSPpnet     = 'x', &   
           f_DMSPdnet     = 'x', f_DMSnet       = 'x', &   
           f_Fednet       = 'x', f_Fepnet       = 'x', &  
           f_Silnet       = 'x', f_PONnet       = 'x', &
           f_zaeronet     = 'x', f_humnet       = 'x', &  
           f_chlsnow      = 'x', f_Nitsnow      = 'x', &   
           f_Amsnow       = 'x', f_Csnow        = 'x', &  
           f_Nsnow        = 'x', f_DMSPpsnow    = 'x', &   
           f_DMSPdsnow    = 'x', f_DMSsnow      = 'x', &    
           f_Fedsnow      = 'x', f_Fepsnow      = 'x', &  
           f_Silsnow      = 'x', f_PONsnow      = 'x', & 
           f_humsnow      = 'x',                       &
           f_DICsnow      = 'x', f_DOCsnow      = 'x', & 
           f_DONsnow      = 'x', f_zaerosnow    = 'x', &
           f_chlfrac      = 'x', f_Nitfrac      = 'x', &   
           f_Amfrac       = 'x',                       &  
           f_Nfrac        = 'x', f_DMSPpfrac    = 'x', &   
           f_DMSPdfrac    = 'x', f_DMSfrac      = 'x', &  
           f_Silfrac      = 'x', f_PONfrac      = 'x', & 
           f_humfrac      = 'x',                       &
           f_DICfrac      = 'x', f_DOCfrac      = 'x', & 
           f_DONfrac      = 'x', f_zaerofrac    = 'x', &
           f_Fedfrac      = 'x', f_Fepfrac      = 'x', &
           f_fNit         = 'x', f_fNit_ai      = 'x', &
           f_fAm          = 'x', f_fAm_ai       = 'x', &
           f_fN           = 'x', f_fN_ai        = 'x', &
           f_fDOC         = 'x', f_fDOC_ai      = 'x', &
           f_fDIC         = 'x', f_fDIC_ai      = 'x', &
           f_fDON         = 'x', f_fDON_ai      = 'x', &
           f_fFed         = 'x', f_fFed_ai      = 'x', &
           f_fFep         = 'x', f_fFep_ai      = 'x', &
           f_fSil         = 'x', f_fSil_ai      = 'x', &  
           f_fPON         = 'x', f_fPON_ai      = 'x', &   
           f_fhum         = 'x', f_fhum_ai      = 'x', &   
           f_fDMSPp       = 'x', f_fDMSPp_ai    = 'x', &   
           f_fDMSPd       = 'x', f_fDMSPd_ai    = 'x', &    
           f_fDMS         = 'x', f_fDMS_ai      = 'x', &   
           f_fzaero       = 'x', f_fzaero_ai    = 'x', & 
           f_bgc_Sil_ml   = 'x', &
           f_bgc_Nit_ml   = 'x', f_bgc_Am_ml    = 'x', &
           f_bgc_DMSP_ml  = 'x', f_bgc_DMS_ml   = 'x', &
           f_bgc_DOC_ml   = 'x', f_bgc_DIC_ml   = 'x', &
           f_bgc_N_ml     = 'x', f_bgc_DON_ml   = 'x', &
           f_peakval      = 'x', f_bgc_Fed_ml   = 'x', &
           f_bgc_Fep_ml   = 'x', f_bgc_hum_ml   = 'x', &
           f_bgc_N_cat1   = 'x', f_bgc_DOC_cat1 = 'x', &
	   f_bgc_DIC_cat1 = 'x', f_bgc_Nit_cat1 = 'x', &
           f_bgc_Am_cat1  = 'x', f_bgc_Sil_cat1 = 'x', &
           f_bgc_DMSPd_cat1= 'x', f_bgc_DMS_cat1 = 'x', &
	   f_bgc_DON_cat1 = 'x', f_bgc_Fed_cat1 = 'x', &
           f_bgc_hum_cat1 = 'x', f_bgc_Fep_cat1 = 'x', &
           f_bgc_PON_cat1 = 'x'

      !---------------------------------------------------------------
      ! namelist variables
      !---------------------------------------------------------------

      namelist / icefields_bgc_nml /     &
           f_faero_atm   , f_faero_ocn   , &
           f_aero        , &
           f_fbio        , f_fbio_ai     , &
           f_zaero       , f_bgc_S       , &
           f_bgc_N       , f_bgc_C       , &
           f_bgc_DOC     , f_bgc_DIC     , &
           f_bgc_chl     , f_bgc_Nit     , &
           f_bgc_Am      , f_bgc_Sil     , &
           f_bgc_DMSPp   , f_bgc_DMSPd   , &
           f_bgc_DMS     , f_bgc_DON     , &
           f_bgc_Fe      , f_bgc_hum     , &
           f_bgc_PON     , f_bgc_ml      , &
           f_upNO        , f_upNH        , &     
           f_bTin        , f_bphi        , &
           f_iDi         , f_iki         , & 
           f_fbri        , f_hbri        , &
           f_zfswin      , f_grownet     , & 
           f_bionet      , f_biosnow     , & 
           f_PPnet       , f_algalpeak   , &
           f_zbgc_frac

      !---------------------------------------------------------------
      ! field indices
      !---------------------------------------------------------------

      integer (kind=int_kind), dimension(max_nstrm), public :: &
           n_fzsal      , n_fzsal_ai   , &  
           n_fzsal_g    , n_fzsal_g_ai , & 
           n_zsal       

      integer(kind=int_kind), dimension(icepack_max_aero,max_nstrm) :: &
           n_faero_atm    , &
           n_faero_ocn    , &
           n_aerosn1      , &
           n_aerosn2      , &
           n_aeroic1      , &
           n_aeroic2      , &
           n_zaeros       , & !  using z tracers
           n_fzaero       , &
           n_fzaero_ai    , &
           n_zaeronet     , &
           n_zaerosnow    , &
           n_zaerofrac

      integer(kind=int_kind), dimension(icepack_max_algae, max_nstrm) :: &
           n_bgc_N  , n_bgc_C   , &
           n_bgc_chl, n_bgc_N_ml, &
           n_fN     , n_fN_ai   , &
           n_Nnet   , n_Nsnow   , n_Nfrac,   &
           n_Cnet   , n_Csnow   ,            &
           n_chlnet , n_chlsnow , n_chlfrac, &
           n_algalpeak          , &
           n_peakval, n_bgc_N_cat1

      integer(kind=int_kind), dimension(icepack_max_doc, max_nstrm) :: &
           n_bgc_DOC,  n_bgc_DOC_ml, &
           n_fDOC   ,  n_fDOC_ai   , &
           n_DOCnet ,  n_DOCsnow   , n_DOCfrac, &
           n_bgc_DOC_cat1

      integer(kind=int_kind), dimension(icepack_max_dic, max_nstrm) :: &
           n_bgc_DIC,  n_bgc_DIC_ml, &
           n_fDIC   ,  n_fDIC_ai   , &
           n_DICnet ,  n_DICsnow   , n_DICfrac, &
           n_bgc_DIC_cat1

      integer(kind=int_kind), dimension(icepack_max_don, max_nstrm) :: &
           n_bgc_DON,  n_bgc_DON_ml, &
           n_fDON   ,  n_fDON_ai   , &
           n_DONnet ,  n_DONsnow   , n_DONfrac, &
           n_bgc_DON_cat1

      integer(kind=int_kind), dimension(icepack_max_fe,  max_nstrm) :: &
           n_bgc_Fed ,  n_bgc_Fed_ml , &
           n_fFed    ,  n_fFed_ai    , &
           n_Fednet  ,  n_Fedsnow    , n_Fedfrac, &
           n_bgc_Fep ,  n_bgc_Fep_ml , &
           n_fFep    ,  n_fFep_ai    , &
           n_Fepnet  ,  n_Fepsnow    , n_Fepfrac, &
           n_bgc_Fed_cat1, n_bgc_Fep_cat1

      integer(kind=int_kind), dimension(max_nstrm) :: &
           n_bgc_S       , &  
           n_fNit        , n_fNit_ai     , &
           n_fAm         , n_fAm_ai      , &
           n_fSil        , n_fSil_ai     , &
           n_fDMS        , n_fDMS_ai     , &
           n_fDMSPp      , n_fDMSPp_ai   , &
           n_fDMSPd      , n_fDMSPd_ai   , &
           n_fPON        , n_fPON_ai     , &
           n_fhum        , n_fhum_ai     , &
           n_bgc_Nit                     , &
           n_bgc_Am      , n_bgc_Sil     , &
           n_bgc_hum     ,                 &
           n_bgc_DMSPp   , n_bgc_DMSPd   , &
           n_bgc_DMS     , n_bgc_Sil_ml  , &
           n_bgc_hum_ml  ,                 &
           n_bgc_Nit_ml  , n_bgc_Am_ml   , &
           n_bgc_DMSP_ml , n_bgc_DMS_ml  , &
           n_upNO        , n_upNH        , & 
           n_bTin        , n_bphi        , &
           n_iDi         , n_iki         , &
           n_bgc_PON     , n_bgc_PON_ml  , &
           n_fbri        , n_hbri        , &
           n_zfswin      , n_Nitnet      , & 
           n_Amnet       , n_Silnet      , &  
           n_humnet      , &
           n_DMSPpnet    , n_DMSPdnet    , &  
           n_DMSnet      , n_PONnet      , &  
           n_Nitsnow     , n_Amsnow      , &
           n_Silsnow     , n_humsnow     , &
           n_DMSPpsnow   , n_DMSPdsnow   , &  
           n_DMSsnow     , n_PONsnow     , &   
           n_Nitfrac     , n_Amfrac      , &
           n_Silfrac     , n_zbgc_frac   , &
           n_humfrac     , &
           n_DMSPpfrac   , n_DMSPdfrac   , &  
           n_DMSfrac     , n_PONfrac     , &  
           n_grownet     , n_PPnet       , &
           n_bgc_Nit_cat1, n_bgc_Am_cat1 , &
           n_bgc_Sil_cat1, n_bgc_DMSPd_cat1,&
           n_bgc_DMS_cat1, n_bgc_PON_cat1, &
           n_bgc_hum_cat1

!=======================================================================

      contains

!=======================================================================

      subroutine init_hist_bgc_2D

      use ice_broadcast, only: broadcast_scalar
      use ice_calendar, only: nstreams
      use ice_communicate, only: my_task, master_task
      use ice_history_shared, only: tstr2D, tcstr, define_hist_field, &
          f_fsalt, f_fsalt_ai, f_sice

      integer (kind=int_kind) :: n, ns
      integer (kind=int_kind) :: nml_error ! namelist i/o error flag
      character (len=3) :: nchar
      character (len=16) :: vname_in     ! variable name
      logical (kind=log_kind) :: tr_zaero, tr_aero, tr_brine, &
          tr_bgc_Nit,    tr_bgc_Am,    tr_bgc_Sil,   &
          tr_bgc_DMS,    tr_bgc_PON,                 &
          tr_bgc_N,      tr_bgc_C,     tr_bgc_chl,   &
          tr_bgc_DON,    tr_bgc_Fe,    tr_bgc_hum,   &
          skl_bgc, solve_zsal, z_tracers
      character(len=*), parameter :: subname = '(init_hist_bgc_2D)'

      call icepack_query_parameters(skl_bgc_out=skl_bgc, &
          solve_zsal_out=solve_zsal, z_tracers_out=z_tracers)
      call icepack_query_tracer_flags(tr_zaero_out =tr_zaero, &
          tr_aero_out   =tr_aero,    tr_brine_out  =tr_brine, &
          tr_bgc_Nit_out=tr_bgc_Nit, tr_bgc_Am_out =tr_bgc_Am, &
          tr_bgc_Sil_out=tr_bgc_Sil, tr_bgc_DMS_out=tr_bgc_DMS, &
          tr_bgc_PON_out=tr_bgc_PON, &
          tr_bgc_N_out  =tr_bgc_N,   tr_bgc_C_out  =tr_bgc_C, &
          tr_bgc_chl_out=tr_bgc_chl, tr_bgc_DON_out=tr_bgc_DON, &
          tr_bgc_Fe_out =tr_bgc_Fe,  tr_bgc_hum_out=tr_bgc_hum ) 
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! read namelist
      !-----------------------------------------------------------------

      call get_fileunit(nu_nml)
      if (my_task == master_task) then
         open (nu_nml, file=nml_filename, status='old',iostat=nml_error)
         if (nml_error /= 0) then
            nml_error = -1
         else
            nml_error =  1
         endif
         do while (nml_error > 0)
            read(nu_nml, nml=icefields_bgc_nml,iostat=nml_error)
         end do
         if (nml_error == 0) close(nu_nml)
      endif
      call release_fileunit(nu_nml)

      call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
         close (nu_nml)
         call abort_ice(subname//'ERROR: reading icefields_bgc_nml')
      endif

      if (.not. tr_aero) then
         f_faero_atm = 'x'
         f_faero_ocn = 'x'
         f_aero      = 'x' 
      endif
 
      if (.not. tr_brine)  then
         f_fbri  = 'x'
         f_hbri  = 'x'
      endif
      
      f_zaeronet  = f_bionet
      f_zaerosnow = f_biosnow
      f_zaerofrac = f_zbgc_frac
      f_fzaero    = f_fbio
      f_fzaero_ai = f_fbio_ai

      if (.not. tr_zaero) then
         f_zaero     = 'x'
         f_fzaero    = 'x' 
         f_fzaero_ai = 'x'
         f_zaeronet  = 'x'
         f_zaerosnow = 'x'
         f_zaerofrac = 'x'
      endif

     if (skl_bgc) then
       f_bionet  = 'x'
       f_biosnow = 'x'
       f_zfswin  = 'x'
       f_zbgc_frac = 'x'
       f_algalpeak = 'x'
     endif

      f_bgc_Nit_ml    = f_bgc_ml
      f_bgc_N_ml      = f_bgc_ml
      f_bgc_DOC_ml    = f_bgc_ml
      f_bgc_DIC_ml    = f_bgc_ml
      f_bgc_Sil_ml    = f_bgc_ml
      f_bgc_hum_ml    = f_bgc_ml
      f_bgc_Am_ml     = f_bgc_ml
      f_bgc_DMSP_ml   = f_bgc_ml
      f_bgc_DMS_ml    = f_bgc_ml
      f_bgc_DON_ml    = f_bgc_ml
      f_bgc_Fed_ml    = f_bgc_ml
      f_bgc_Fep_ml    = f_bgc_ml

      f_Nitnet    = f_bionet
      f_Amnet     = f_bionet
      f_Nnet      = f_bionet
      f_chlnet    = f_bionet
      f_Cnet      = f_bionet
      f_DOCnet    = f_bionet
      f_DICnet    = f_bionet
      f_DONnet    = f_bionet
      f_Fednet    = f_bionet
      f_Fepnet    = f_bionet
      f_Silnet    = f_bionet
      f_humnet    = f_bionet
      f_DMSPpnet  = f_bionet
      f_DMSPdnet  = f_bionet
      f_DMSnet    = f_bionet
      f_PONnet    = f_bionet
     
      f_Nitsnow    = f_biosnow
      f_Amsnow     = f_biosnow
      f_Nsnow      = f_biosnow
      f_chlsnow    = f_biosnow
      f_Csnow      = f_biosnow
      f_DOCsnow    = f_biosnow
      f_DICsnow    = f_biosnow
      f_DONsnow    = f_biosnow
      f_Fedsnow    = f_biosnow
      f_Fepsnow    = f_biosnow
      f_Silsnow    = f_biosnow
      f_humsnow    = f_biosnow
      f_DMSPpsnow  = f_biosnow
      f_DMSPdsnow  = f_biosnow
      f_DMSsnow    = f_biosnow
      f_PONsnow    = f_biosnow

      f_Nitfrac    = f_zbgc_frac
      f_Amfrac     = f_zbgc_frac
      f_Nfrac      = f_zbgc_frac
      f_chlfrac    = f_zbgc_frac
      f_DOCfrac    = f_zbgc_frac
      f_DICfrac    = f_zbgc_frac
      f_DONfrac    = f_zbgc_frac
      f_Fedfrac    = f_zbgc_frac
      f_Fepfrac    = f_zbgc_frac
      f_Silfrac    = f_zbgc_frac
      f_humfrac    = f_zbgc_frac
      f_DMSPpfrac  = f_zbgc_frac
      f_DMSPdfrac  = f_zbgc_frac
      f_DMSfrac    = f_zbgc_frac
      f_PONfrac    = f_zbgc_frac

      f_bgc_DMSPp  = f_bgc_DMS
      f_bgc_DMSPd  = f_bgc_DMS

      f_bgc_Fed = f_bgc_Fe
      f_bgc_Fep = f_bgc_Fe

      f_fDOC    = f_fbio
      f_fDIC    = f_fbio
      f_fDON    = f_fbio
      f_fFed    = f_fbio
      f_fFep    = f_fbio
      f_fNit    = f_fbio
      f_fAm     = f_fbio
      f_fN      = f_fbio
      f_fSil    = f_fbio
      f_fhum    = f_fbio
      f_fPON    = f_fbio
      f_fDMSPp  = f_fbio
      f_fDMSPd  = f_fbio
      f_fDMS    = f_fbio

      f_fDOC_ai    = f_fbio_ai
      f_fDIC_ai    = f_fbio_ai
      f_fDON_ai    = f_fbio_ai
      f_fFed_ai    = f_fbio_ai
      f_fFep_ai    = f_fbio_ai
      f_fNit_ai    = f_fbio_ai
      f_fAm_ai     = f_fbio_ai
      f_fN_ai      = f_fbio_ai
      f_fSil_ai    = f_fbio_ai
      f_fhum_ai    = f_fbio_ai
      f_fPON_ai    = f_fbio_ai
      f_fDMSPp_ai  = f_fbio_ai
      f_fDMSPd_ai  = f_fbio_ai
      f_fDMS_ai    = f_fbio_ai

     if (.not. tr_bgc_N) then  
        f_bgc_N     = 'x'
        f_bgc_N_ml  = 'x'
        f_fN        = 'x'
        f_fN_ai     = 'x'
        f_Nsnow     = 'x'
        f_Nfrac     = 'x'
        f_Nnet      = 'x'
        f_algalpeak = 'x'
     endif

     f_peakval = f_algalpeak
     if (.not. tr_bgc_Nit) then  
        f_upNO      = 'x'  
        f_bgc_Nit   = 'x'
        f_bgc_Nit_ml= 'x'
        f_fNit      = 'x'
        f_fNit_ai   = 'x'
        f_Nitsnow   = 'x'
        f_Nitfrac   = 'x'
        f_Nitnet    = 'x'
     endif
     if (.not. tr_bgc_C) then
        f_bgc_C   = 'x'
        f_bgc_DOC = 'x'
        f_bgc_DIC = 'x'
        f_fDOC    = 'x'
        f_fDOC_ai = 'x'
        f_bgc_DOC_ml = 'x'
        f_bgc_DIC_ml = 'x'
        f_Csnow      = 'x'
        f_DOCsnow    = 'x'
        f_DICsnow    = 'x'
        f_DOCfrac    = 'x'
        f_DICfrac    = 'x'
        f_Cnet       = 'x'
        f_DOCnet     = 'x'
        f_DICnet     = 'x'
     endif
     if (.not. tr_bgc_chl) then
        f_bgc_chl = 'x'
        f_chlnet  = 'x'
        f_chlsnow = 'x'
        f_chlfrac = 'x'
     endif
     if (.not. tr_bgc_Am) then  
        f_upNH      = 'x'  
        f_bgc_Am    = 'x'
        f_bgc_Am_ml = 'x'
        f_fAm       = 'x'
        f_fAm_ai    = 'x'
        f_Amsnow    = 'x'
        f_Amfrac    = 'x'
        f_Amnet     = 'x'
     endif
     if (.not. tr_bgc_Sil) then
        f_bgc_Sil    = 'x'
        f_bgc_Sil_ml = 'x'
        f_fSil       = 'x'
        f_fSil_ai    = 'x'
        f_Silnet     = 'x'
        f_Silsnow    = 'x'
        f_Silfrac    = 'x'
     endif
     if (.not. tr_bgc_hum) then
        f_bgc_hum    = 'x'
        f_bgc_hum_ml = 'x'
        f_fhum       = 'x'
        f_fhum_ai    = 'x'
        f_humnet     = 'x'
        f_humsnow    = 'x'
        f_humfrac    = 'x'
     endif
     if (.not. tr_bgc_DMS) then
        f_bgc_DMS   = 'x'
        f_bgc_DMSPp = 'x'
        f_bgc_DMSPd = 'x'
        f_bgc_DMSP_ml= 'x'
        f_bgc_DMS_ml = 'x'
        f_fDMS       = 'x'
        f_fDMSPp     = 'x'
        f_fDMSPd     = 'x'
        f_fDMS_ai    = 'x'
        f_fDMSPp_ai  = 'x'
        f_fDMSPd_ai  = 'x'
        f_DMSnet     = 'x'
        f_DMSPpnet   = 'x'
        f_DMSPdnet   = 'x'
        f_DMSsnow    = 'x'
        f_DMSPpsnow  = 'x'
        f_DMSPdsnow  = 'x'
        f_DMSfrac    = 'x'
        f_DMSPpfrac  = 'x'
        f_DMSPdfrac  = 'x'
     endif   
     if (.not. tr_bgc_DON) then 
        f_bgc_DON    = 'x'
        f_bgc_DON_ml = 'x'
        f_DONsnow    = 'x'
        f_DONfrac    = 'x'
        f_DONnet     = 'x'
        f_fDON       = 'x'
        f_fDON_ai    = 'x'
     endif 
     if (.not. tr_bgc_Fe ) then 
        f_bgc_Fe     = 'x'
        f_bgc_Fed    = 'x'
        f_bgc_Fed_ml = 'x'
        f_Fedsnow    = 'x'
        f_Fedfrac    = 'x'
        f_Fednet     = 'x'
        f_fFed       = 'x'
        f_fFed_ai    = 'x'
        f_bgc_Fep    = 'x'
        f_bgc_Fep_ml = 'x'
        f_Fepsnow    = 'x'
        f_Fepfrac    = 'x'
        f_Fepnet     = 'x'
        f_fFep       = 'x'
        f_fFep_ai    = 'x'
     endif
     if (.not. tr_bgc_PON .or. skl_bgc) then 
        f_bgc_PON   = 'x'
        f_PONsnow    = 'x'
        f_PONfrac    = 'x'
        f_PONnet     = 'x'
        f_fPON       = 'x'
        f_fPON_ai    = 'x'
     endif
       
      f_bgc_Nit_cat1    = f_bgc_Nit 
      f_bgc_Am_cat1     = f_bgc_Am   
      f_bgc_N_cat1      = f_bgc_N
      f_bgc_DOC_cat1    = f_bgc_DOC
      f_bgc_DIC_cat1    = f_bgc_DIC
      f_bgc_DON_cat1    = f_bgc_DON
      f_bgc_Fed_cat1    = f_bgc_Fe 
      f_bgc_Fep_cat1    = f_bgc_Fe  
      f_bgc_Sil_cat1    = f_bgc_Sil 
      f_bgc_hum_cat1    = f_bgc_hum 
      f_bgc_DMSPd_cat1  = f_bgc_DMSPd
      f_bgc_DMS_cat1    = f_bgc_DMS 
      f_bgc_PON_cat1    = f_bgc_PON

      if (solve_zsal) then
         f_fzsal = f_fsalt
         f_fzsal_g = f_fsalt
         f_fzsal_ai = f_fsalt_ai
         f_fzsal_g_ai = f_fsalt_ai
         f_zsal = f_sice
         f_fsalt = 'x'
         f_fsalt_ai = 'x'
         f_sice = 'x'
      else
         f_fzsal      = 'x'
         f_fzsal_g    = 'x'
         f_fzsal_ai   = 'x'
         f_fzsal_g_ai = 'x'
         f_zsal       = 'x'
         f_bgc_S  = 'x'
         f_iki    = 'x'
      endif

      call broadcast_scalar (f_faero_atm,    master_task)
      call broadcast_scalar (f_faero_ocn,    master_task)
      call broadcast_scalar (f_aero,         master_task)
      call broadcast_scalar (f_fbri,         master_task)
      call broadcast_scalar (f_hbri,         master_task)
      call broadcast_scalar (f_fzsal,        master_task)
      call broadcast_scalar (f_fzsal_ai,     master_task)
      call broadcast_scalar (f_fzsal_g,      master_task)
      call broadcast_scalar (f_fzsal_g_ai,   master_task)
      call broadcast_scalar (f_zsal,         master_task)
      call broadcast_scalar (f_fNit,         master_task)
      call broadcast_scalar (f_fNit_ai,      master_task)
      call broadcast_scalar (f_fDOC,         master_task)
      call broadcast_scalar (f_fDOC_ai,      master_task)
      call broadcast_scalar (f_fDIC,         master_task)
      call broadcast_scalar (f_fDIC_ai,      master_task)
      call broadcast_scalar (f_fDON,         master_task)
      call broadcast_scalar (f_fDON_ai,      master_task)
      call broadcast_scalar (f_fFed,         master_task)
      call broadcast_scalar (f_fFed_ai,      master_task)
      call broadcast_scalar (f_fFep,         master_task)
      call broadcast_scalar (f_fFep_ai,      master_task)
      call broadcast_scalar (f_fAm,          master_task)
      call broadcast_scalar (f_fAm_ai,       master_task)
      call broadcast_scalar (f_fN,           master_task)
      call broadcast_scalar (f_fN_ai,        master_task)
      call broadcast_scalar (f_fSil,         master_task)
      call broadcast_scalar (f_fSil_ai,      master_task)
      call broadcast_scalar (f_fhum,         master_task)
      call broadcast_scalar (f_fhum_ai,      master_task)
      call broadcast_scalar (f_fPON,         master_task)
      call broadcast_scalar (f_fPON_ai,      master_task)
      call broadcast_scalar (f_fDMS,         master_task)
      call broadcast_scalar (f_fDMS_ai,      master_task)
      call broadcast_scalar (f_fDMSPp,       master_task)
      call broadcast_scalar (f_fDMSPp_ai,    master_task)
      call broadcast_scalar (f_fDMSPd,       master_task)
      call broadcast_scalar (f_fDMSPd_ai,    master_task)
      call broadcast_scalar (f_fzaero,       master_task)
      call broadcast_scalar (f_fzaero_ai,    master_task)
      call broadcast_scalar (f_zaero,        master_task)
      call broadcast_scalar (f_bgc_N,        master_task)
      call broadcast_scalar (f_bgc_C,        master_task)
      call broadcast_scalar (f_bgc_DOC,      master_task)
      call broadcast_scalar (f_bgc_DIC,      master_task)
      call broadcast_scalar (f_bgc_chl,      master_task)
      call broadcast_scalar (f_bgc_Nit,      master_task)
      call broadcast_scalar (f_bgc_Am,       master_task)
      call broadcast_scalar (f_bgc_Sil,      master_task)
      call broadcast_scalar (f_bgc_hum,      master_task)
      call broadcast_scalar (f_bgc_DMSPp,    master_task)
      call broadcast_scalar (f_bgc_DMSPd,    master_task)
      call broadcast_scalar (f_bgc_DMS,      master_task)
      call broadcast_scalar (f_bgc_PON,      master_task)
      call broadcast_scalar (f_bgc_DON,      master_task)
      call broadcast_scalar (f_bgc_Fe,       master_task)
      call broadcast_scalar (f_bgc_Fed,      master_task)
      call broadcast_scalar (f_bgc_Fep,      master_task)
      call broadcast_scalar (f_bgc_N_cat1,   master_task)
      call broadcast_scalar (f_bgc_DOC_cat1, master_task)
      call broadcast_scalar (f_bgc_DIC_cat1, master_task)
      call broadcast_scalar (f_bgc_Nit_cat1, master_task)
      call broadcast_scalar (f_bgc_Am_cat1,  master_task)
      call broadcast_scalar (f_bgc_Sil_cat1, master_task)
      call broadcast_scalar (f_bgc_hum_cat1, master_task)
      call broadcast_scalar (f_bgc_DMSPd_cat1, master_task)
      call broadcast_scalar (f_bgc_DMS_cat1, master_task)
      call broadcast_scalar (f_bgc_PON_cat1, master_task)
      call broadcast_scalar (f_bgc_DON_cat1, master_task)
      call broadcast_scalar (f_bgc_Fed_cat1, master_task)
      call broadcast_scalar (f_bgc_Fep_cat1, master_task)
      call broadcast_scalar (f_bgc_Nit_ml,   master_task)
      call broadcast_scalar (f_bgc_DOC_ml,   master_task)
      call broadcast_scalar (f_bgc_DIC_ml,   master_task)
      call broadcast_scalar (f_bgc_N_ml,     master_task)
      call broadcast_scalar (f_bgc_Am_ml,    master_task)
      call broadcast_scalar (f_bgc_Sil_ml,   master_task)
      call broadcast_scalar (f_bgc_hum_ml,   master_task)
      call broadcast_scalar (f_bgc_DMSP_ml,  master_task)
      call broadcast_scalar (f_bgc_DMS_ml,   master_task)  
      call broadcast_scalar (f_bgc_DON_ml,   master_task)  
      call broadcast_scalar (f_bgc_Fed_ml,   master_task)  
      call broadcast_scalar (f_bgc_Fep_ml,   master_task)  
      call broadcast_scalar (f_upNO,         master_task)  
      call broadcast_scalar (f_upNH,         master_task) 
      call broadcast_scalar (f_bTin,         master_task)
      call broadcast_scalar (f_bphi,         master_task)
      call broadcast_scalar (f_iDi,          master_task)  
      call broadcast_scalar (f_iki,          master_task) 
      call broadcast_scalar (f_bgc_S,        master_task)  
      call broadcast_scalar (f_zfswin,       master_task) 
      call broadcast_scalar (f_PPnet,        master_task) 
      call broadcast_scalar (f_algalpeak,    master_task)  
      call broadcast_scalar (f_zbgc_frac,    master_task)  
      call broadcast_scalar (f_peakval,      master_task) 
      call broadcast_scalar (f_grownet,      master_task) 
      call broadcast_scalar (f_chlnet,       master_task) 
      call broadcast_scalar (f_Nitnet,       master_task)  
      call broadcast_scalar (f_Nnet,         master_task)  
      call broadcast_scalar (f_Cnet,         master_task)   
      call broadcast_scalar (f_DOCnet,       master_task)   
      call broadcast_scalar (f_DICnet,       master_task)  
      call broadcast_scalar (f_Amnet,        master_task)  
      call broadcast_scalar (f_Silnet,       master_task)   
      call broadcast_scalar (f_humnet,       master_task)  
      call broadcast_scalar (f_DMSPpnet,     master_task)  
      call broadcast_scalar (f_DMSPdnet,     master_task)  
      call broadcast_scalar (f_DMSnet,       master_task)  
      call broadcast_scalar (f_PONnet,       master_task)    
      call broadcast_scalar (f_DONnet,       master_task)    
      call broadcast_scalar (f_Fednet,       master_task)     
      call broadcast_scalar (f_Fepnet,       master_task)      
      call broadcast_scalar (f_zaeronet,     master_task)    
      call broadcast_scalar (f_chlsnow,      master_task) 
      call broadcast_scalar (f_Nitsnow,      master_task)  
      call broadcast_scalar (f_Nsnow,        master_task)  
      call broadcast_scalar (f_Csnow,        master_task)   
      call broadcast_scalar (f_DOCsnow,      master_task)   
      call broadcast_scalar (f_DICsnow,      master_task)  
      call broadcast_scalar (f_Amsnow,       master_task)  
      call broadcast_scalar (f_Silsnow,      master_task)   
      call broadcast_scalar (f_humsnow,      master_task)  
      call broadcast_scalar (f_DMSPpsnow,    master_task)  
      call broadcast_scalar (f_DMSPdsnow,    master_task)  
      call broadcast_scalar (f_DMSsnow,      master_task)   
      call broadcast_scalar (f_PONsnow,      master_task)    
      call broadcast_scalar (f_DONsnow,      master_task)    
      call broadcast_scalar (f_Fedsnow,      master_task)    
      call broadcast_scalar (f_Fepsnow,      master_task)      
      call broadcast_scalar (f_zaerosnow,    master_task)    
      call broadcast_scalar (f_chlfrac,      master_task) 
      call broadcast_scalar (f_Nitfrac,      master_task)  
      call broadcast_scalar (f_Nfrac,        master_task)  
      call broadcast_scalar (f_DOCfrac,      master_task)   
      call broadcast_scalar (f_DICfrac,      master_task)  
      call broadcast_scalar (f_Amfrac,       master_task)  
      call broadcast_scalar (f_Silfrac,      master_task)  
      call broadcast_scalar (f_humfrac,      master_task)  
      call broadcast_scalar (f_DMSPpfrac,    master_task)  
      call broadcast_scalar (f_DMSPdfrac,    master_task)  
      call broadcast_scalar (f_DMSfrac,      master_task)   
      call broadcast_scalar (f_PONfrac,      master_task)    
      call broadcast_scalar (f_DONfrac,      master_task)    
      call broadcast_scalar (f_Fedfrac,      master_task)    
      call broadcast_scalar (f_Fepfrac,      master_task)      
      call broadcast_scalar (f_zaerofrac,    master_task) 

      ! 2D variables

    if (tr_aero .or. tr_brine .or. solve_zsal .or. skl_bgc) then

    do ns = 1, nstreams

      ! zsalinity 
     
      call define_hist_field(n_fzsal,"fzsal","kg/m^2/s",tstr2D, tcstr, &
          "prognostic salt flux ice to ocn (cpl)",                     &
          "if positive, ocean gains salt", c1, c0,                     &
          ns, f_fzsal)
      
      call define_hist_field(n_fzsal_ai,"fzsal_ai","kg/m^2/s",tstr2D, tcstr, &
          "prognostic salt flux ice to ocean",                         &
          "weighted by ice area", c1, c0,                              &
          ns, f_fzsal_ai)
      
      call define_hist_field(n_fzsal_g,"fzsal_g","kg/m^2/s",tstr2D, tcstr, &
          "Gravity drainage salt flux ice to ocn (cpl)",               &
          "if positive, ocean gains salt", c1, c0,                     &
          ns, f_fzsal_g)
      
      call define_hist_field(n_fzsal_g_ai,"fzsal_g_ai","kg/m^2/s",tstr2D, tcstr, &
          "Gravity drainage salt flux ice to ocean",                   &
          "weighted by ice area", c1, c0,                              &
          ns, f_fzsal_g_ai)
      
      call define_hist_field(n_zsal,"zsal_tot","g/m^2",tstr2D, tcstr,  &
          "Total Salt content",                                        &
          "In ice volume*fbri", c1, c0,                                &
          ns, f_zsal)

      ! Aerosols
      if (f_aero(1:1) /= 'x') then
         do n=1,n_aero
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'aerosnossl', trim(nchar)
            call define_hist_field(n_aerosn1(n,:),vname_in,"kg/kg",   &
                tstr2D, tcstr,"snow ssl aerosol mass","none", c1, c0, &
                ns, f_aero)
            write(vname_in,'(a,a)') 'aerosnoint', trim(nchar)
            call define_hist_field(n_aerosn2(n,:),vname_in,"kg/kg",   &
                tstr2D, tcstr,"snow int aerosol mass","none", c1, c0, &
                ns, f_aero)
            write(vname_in,'(a,a)') 'aeroicessl', trim(nchar)
            call define_hist_field(n_aeroic1(n,:),vname_in,"kg/kg",  &
                tstr2D, tcstr,"ice ssl aerosol mass","none", c1, c0, &
                ns, f_aero)
            write(vname_in,'(a,a)') 'aeroiceint', trim(nchar)
            call define_hist_field(n_aeroic2(n,:),vname_in,"kg/kg",  &
                tstr2D, tcstr,"ice int aerosol mass","none", c1, c0, &
                ns, f_aero)
         enddo
      endif

      if (f_faero_atm(1:1) /= 'x') then
         do n=1,n_aero
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'faero_atm', trim(nchar)
            call define_hist_field(n_faero_atm(n,:),vname_in,"kg/m^2 s", &
                tstr2D, tcstr,"aerosol deposition rate","none", c1, c0,  &
                ns, f_faero_atm)
         enddo
      endif

      if (f_faero_ocn(1:1) /= 'x') then
         do n=1,n_aero
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'faero_ocn', trim(nchar)
            call define_hist_field(n_faero_ocn(n,:),vname_in,"kg/m^2 s", &
                tstr2D, tcstr,"aerosol flux to ocean","none", c1, c0,    &
                ns, f_faero_ocn)
         enddo
      endif

      if (skl_bgc) then
        ! skeletal layer tracers


        if (f_bgc_N(1:1) /= 'x') then
          do n = 1, n_algae
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'algal_N', trim(nchar)
            call define_hist_field(n_bgc_N(n,:),vname_in,"mmol/m^2",tstr2D, tcstr, &
             "Bulk ice bottom algae (nitrogen)",                           &
             "skeletal layer: bottom 3 cm", c1, c0,                &
             ns, f_bgc_N)
         enddo
        endif
        if (f_bgc_chl(1:1) /= 'x') then
          do n = 1, n_algae
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'algal_chl', trim(nchar)
            call define_hist_field(n_bgc_chl(n,:),vname_in,"mg chl/m^2",tstr2D, tcstr, &
             "Bulk ice bottom algae (chlorophyll)",                             &
             "skeletal layer: bottom 3 cm", c1, c0,                &
             ns, f_bgc_chl)
          enddo
        endif !f_bgc_chl
        if (f_bgc_C(1:1) /= 'x') then
          do n = 1, n_algae
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'algal_C', trim(nchar)
            call define_hist_field(n_bgc_C(n,:),vname_in,"mmol C/m^2",tstr2D, tcstr, &
             "Bulk ice bottom diatoms (carbon)",                             &
             "skeletal layer: bottom 3 cm", c1, c0,                &
             ns, f_bgc_C)
          enddo
        endif
        if (f_bgc_DOC(1:1) /= 'x') then
          do n = 1, n_doc
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'DOC', trim(nchar)
            call define_hist_field(n_bgc_DOC(n,:),vname_in,"mmol/m^2",tstr2D, tcstr, &
             "Bulk DOC (carbon)",                                      &
             "skeletal layer: bottom 3 cm", c1, c0,                &
             ns, f_bgc_DOC)
          enddo
        endif
        if (f_bgc_DIC(1:1) /= 'x') then
          do n = 1, n_dic
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'DIC', trim(nchar)
            call define_hist_field(n_bgc_DIC(n,:),vname_in,"mmol/m^2",tstr2D, tcstr, &
             "Bulk DIC (carbon)",                                      &
             "skeletal layer: bottom 3 cm", c1, c0,                &
             ns, f_bgc_DIC)
          enddo
        endif
        if (f_bgc_DON(1:1) /= 'x') then
          do n = 1, n_don
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'DON', trim(nchar)
            call define_hist_field(n_bgc_DON(n,:),vname_in,"mmol/m^2",tstr2D, tcstr, &
             "Bulk ice bottom DON (nitrogen)",                             &
             "skeletal layer: bottom 3 cm", c1, c0,                &
             ns, f_bgc_DON)
          enddo
        endif !f_bgc_DON
        if (f_bgc_Fe (1:1) /= 'x') then
          do n = 1, n_fed
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'dFe', trim(nchar)
            call define_hist_field(n_bgc_Fed (n,:),vname_in,"umol/m^2",tstr2D, tcstr, &
             "Bulk ice bottom dissolved Fe (iron)",                             &
             "skeletal layer: bottom 3 cm", c1, c0,                &
             ns, f_bgc_Fed )
          enddo
          do n = 1, n_fep
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'pFe', trim(nchar)
            call define_hist_field(n_bgc_Fep (n,:),vname_in,"umol/m^2",tstr2D, tcstr, &
             "Bulk ice bottom particulate Fe (iron)",                           &
             "skeletal layer: bottom 3 cm", c1, c0,                &
             ns, f_bgc_Fep )
          enddo
        endif !f_bgc_Fe 
 
        if (f_bgc_Nit(1:1) /= 'x') &
          call define_hist_field(n_bgc_Nit,"Nit","mmol/m^2",tstr2D, tcstr, &
             "Bulk skeletal nutrient (nitrate)",                             &
             "skeletal layer: bottom 3 cm", c1, c0,                &
             ns, f_bgc_Nit)
        if (f_bgc_Am(1:1) /= 'x') &
          call define_hist_field(n_bgc_Am,"Am","mmol/m^2",tstr2D, tcstr, &
             "Bulk skeletal nutrient (ammonia/um)",                        &
             "skeletal layer: bottom 3 cm", c1, c0,                &
             ns, f_bgc_Am)
        if (f_bgc_Sil(1:1) /= 'x') &
          call define_hist_field(n_bgc_Sil,"Sil","mmol/m^2",tstr2D, tcstr, &
             "Bulk skeletal nutrient (silicate)",                            &
             "skelelal layer: bottom 3 cm", c1, c0,                &
             ns, f_bgc_Sil)
        if (f_bgc_hum(1:1) /= 'x') &
          call define_hist_field(n_bgc_hum,"hum","mmol/m^2",tstr2D, tcstr, &
             "Bulk skeletal humic material (carbon)",              &
             "skelelal layer: bottom 3 cm", c1, c0,                &
             ns, f_bgc_hum)
        if (f_bgc_PON(1:1) /= 'x') &
          call define_hist_field(n_bgc_PON,"PON","mmol/m^2",tstr2D, tcstr, &
             "Bulk skeletal nutrient (silicate)",                            &
             "skelelal layer: bottom 3 cm", c1, c0,                &
             ns, f_bgc_PON)
        if (f_bgc_DMSPp(1:1) /= 'x') &
          call define_hist_field(n_bgc_DMSPp,"DMSPp","mmol/m^2",tstr2D, tcstr, &
             "Bulk particulate S in algae (DMSPp)",                              &
             "skeletal layer: bottom 3 cm", c1, c0,                &
             ns, f_bgc_DMSPp)
        if (f_bgc_DMSPd(1:1) /= 'x') &
          call define_hist_field(n_bgc_DMSPd,"DMSPd","mmol/m^2",tstr2D, tcstr, &
             "Bulk dissolved skl precursor (DSMPd)",                             &
             "skeletal layer: bottom 3 cm", c1, c0,                &
             ns, f_bgc_DMSPd)
        if (f_bgc_DMS(1:1) /= 'x') &
          call define_hist_field(n_bgc_DMS,"DMS","mmol/m^2",tstr2D, tcstr, &
             "Bulk dissolved skl trace gas (DMS)",                           &
             "skeletal layer: bottom 3 cm", c1, c0,                &
             ns, f_bgc_DMS)
 
      endif  !skl_bgc

      ! vertical and skeletal layer biogeochemistry

        if (f_bgc_DOC_ml(1:1) /= 'x') then
          do n = 1, n_doc
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'ml_DOC', trim(nchar)
            call define_hist_field(n_bgc_DOC_ml(n,:),vname_in,"mmol/m^3",tstr2D, tcstr, &
             "mixed layer DOC (carbon)",                             &
             "upper ocean", c1, c0,                &
             ns, f_bgc_DOC_ml)
          enddo
        endif
        if (f_bgc_DIC_ml(1:1) /= 'x') then
          do n = 1, n_dic
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'ml_DIC', trim(nchar)
            call define_hist_field(n_bgc_DIC_ml(n,:),vname_in,"mmol/m^3",tstr2D, tcstr, &
             "mixed layer DIC (carbon)",                             &
             "upper ocean", c1, c0,                &
             ns, f_bgc_DIC_ml)
          enddo
        endif
        if (f_bgc_DON_ml(1:1) /= 'x') then
          do n = 1, n_don
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'ml_DON', trim(nchar)
            call define_hist_field(n_bgc_DON_ml(n,:),vname_in,"mmol/m^3",tstr2D, tcstr, &
             "mixed layer DON (nitrogen)",                           &
             "upper ocean", c1, c0,                &
             ns, f_bgc_DON_ml)
          enddo
        endif
        if (f_bgc_Fed_ml (1:1) /= 'x') then
          do n = 1, n_fed 
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'ml_dFe', trim(nchar)
            call define_hist_field(n_bgc_Fed_ml (n,:),vname_in,"nM",tstr2D, tcstr, &
             "mixed layer dissolved Fe (iron)",                           &
             "upper ocean", c1, c0,                &
             ns, f_bgc_Fed_ml )
          enddo
        endif
        if (f_bgc_Fep_ml (1:1) /= 'x') then
          do n = 1, n_fep 
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'ml_pFe', trim(nchar)
            call define_hist_field(n_bgc_Fep_ml (n,:),vname_in,"nM",tstr2D, tcstr, &
             "mixed layer particulate Fe (iron)",                           &
             "upper ocean", c1, c0,                &
             ns, f_bgc_Fep_ml )
          enddo
        endif
        if (f_bgc_N_ml(1:1) /= 'x') then
          do n = 1, n_algae
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'ml_N', trim(nchar)
            call define_hist_field(n_bgc_N_ml(n,:),vname_in,"mmol/m^3",tstr2D, tcstr, &
             "mixed layer nitrogen",                             &
             "upper ocean", c1, c0,                &
             ns, f_bgc_N_ml)
          enddo
        endif
        if (f_bgc_Nit_ml(1:1) /= 'x') &
          call define_hist_field(n_bgc_Nit_ml,"ml_Nit","mmol/m^3",tstr2D, tcstr, &
             "mixed layer nutrient (nitrate)",                         &
             "upper ocean", c1, c0,                                  &
             ns, f_bgc_Nit_ml)
        if (f_bgc_Am_ml(1:1) /= 'x') &
          call define_hist_field(n_bgc_Am_ml,"ml_Am","mmol/m^3",tstr2D, tcstr, &
             "mixed layer nutrient (ammonia/um)",                    &
             "upper ocean", c1, c0,                                  &
             ns, f_bgc_Am_ml)
        if (f_bgc_Sil_ml(1:1) /= 'x') &
          call define_hist_field(n_bgc_Sil_ml,"ml_Sil","mmol/m^3",tstr2D, tcstr, &
             "mixed layer nutrient (silicate)",                        &
             "upper ocean", c1, c0,                                  &
             ns, f_bgc_Sil_ml)
        if (f_bgc_hum_ml(1:1) /= 'x') &
          call define_hist_field(n_bgc_hum_ml,"ml_hum","mmol/m^3",tstr2D, tcstr, &
             "mixed layer humic material (carbon)",                             &
             "upper ocean", c1, c0,                                  &
             ns, f_bgc_hum_ml) 
        if (f_bgc_DMSP_ml(1:1) /= 'x') &
          call define_hist_field(n_bgc_DMSP_ml,"ml_DMSP","mmol/m^3",tstr2D, tcstr, &
             "mixed layer precursor (DMSP)",                             &
             "upper ocean", c1, c0,                                  &
             ns, f_bgc_DMSP_ml)
        if (f_bgc_DMS_ml(1:1) /= 'x') &
          call define_hist_field(n_bgc_DMS_ml,"ml_DMS","mmol/m^3",tstr2D, tcstr, &
             "mixed layer trace gas (DMS)",                            &
             "upper ocean", c1, c0,                                  &
             ns, f_bgc_DMS_ml)
            
        if (f_fNit(1:1) /= 'x') &
          call define_hist_field(n_fNit,"fNit","mmol/m^2/s",tstr2D, tcstr, &
             "nitrate flux ice to ocn (cpl)",                           &
             "if positive, ocean gains nitrate", c1, c0,                &
             ns, f_fNit)
      
        if (f_fNit_ai(1:1) /= 'x') &
          call define_hist_field(n_fNit_ai,"fNit_ai","mmol/m^2/s",tstr2D, tcstr, &
             "nitrate flux ice to ocean",                            &
             "weighted by ice area", c1, c0,                         &
             ns, f_fNit_ai)
           
        if (f_fAm(1:1) /= 'x') &
          call define_hist_field(n_fAm,"fAm","mmol/m^2/s",tstr2D, tcstr, &
             "ammonium flux ice to ocn (cpl)",                     &
             "if positive, ocean gains ammonium", c1, c0,          &
             ns, f_fAm)
 
        if (f_fAm_ai(1:1) /= 'x') &
          call define_hist_field(n_fAm_ai,"fAm_ai","mmol/m^2/s",tstr2D, tcstr, &
             "ammonium flux ice to ocean",                           &
             "weighted by ice area", c1, c0,                         &
             ns, f_fAm_ai)             
        if (f_fN(1:1) /= 'x') then
          do n = 1, n_algae
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'fN', trim(nchar)
            call define_hist_field(n_fN(n,:),vname_in,"mmol/m^2/s",tstr2D, tcstr, &
             "algal N flux ice to ocn (cpl)",                     &
             "if positive, ocean gains algal N", c1, c0,          &
             ns, f_fN)
          enddo
        endif
        if (f_fN_ai(1:1) /= 'x') then
          do n = 1, n_algae
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'fN_ai', trim(nchar)
            call define_hist_field(n_fN_ai(n,:),vname_in,"mmol/m^2/s",tstr2D, tcstr, &
             "algal N flux ice to ocean",                            &
             "weighted by ice area", c1, c0,                         &
             ns, f_fN_ai)
          enddo
        endif
        if (f_fDOC(1:1) /= 'x') then
          do n = 1, n_doc
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'fDOC', trim(nchar)
            call define_hist_field(n_fDOC(n,:),vname_in,"mmol/m^2/s",tstr2D, tcstr, &
             "DOC flux ice to ocn (cpl)",                     &
             "positive to ocean", c1, c0,                   &
             ns, f_fDOC)
          enddo
        endif
        if (f_fDOC_ai(1:1) /= 'x') then
          do n = 1, n_doc
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'fDOC_ai', trim(nchar)
            call define_hist_field(n_fDOC_ai(n,:),vname_in,"mmol/m^2/s",tstr2D, tcstr, &
             "DOC flux ice to ocn",                              &
             "weighted by ice area", c1, c0,                         &
             ns, f_fDOC_ai)
          enddo
        endif        
        if (f_fDIC(1:1) /= 'x') then
          do n = 1, n_dic
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'fDIC', trim(nchar)
            call define_hist_field(n_fDIC(n,:),vname_in,"mmol/m^2/s",tstr2D, tcstr, &
             "DIC flux ice to ocn (cpl)",                     &
             "positive to ocean", c1, c0,                   &
             ns, f_fDIC)
          enddo
        endif
        if (f_fDIC_ai(1:1) /= 'x') then
          do n = 1, n_dic
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'fDIC_ai', trim(nchar)
            call define_hist_field(n_fDIC_ai(n,:),vname_in,"mmol/m^2/s",tstr2D, tcstr, &
             "DIC flux ice to ocn",                              &
             "weighted by ice area", c1, c0,                         &
             ns, f_fDIC_ai)
          enddo
        endif               
        if (f_fDON(1:1) /= 'x') then
          do n = 1, n_don
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'fDON', trim(nchar)
            call define_hist_field(n_fDON(n,:),vname_in,"mmol/m^2/s",tstr2D, tcstr, &
             "DON flux ice to ocn (cpl)",                     &
             "positive to ocean", c1, c0,                   &
             ns, f_fDON)
          enddo
        endif
        if (f_fDON_ai(1:1) /= 'x') then
          do n = 1, n_don
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'fDON_ai', trim(nchar)
            call define_hist_field(n_fDON_ai(n,:),vname_in,"mmol/m^2/s",tstr2D, tcstr, &
             "DON flux ice to ocn",                              &
             "weighted by ice area", c1, c0,                         &
             ns, f_fDON_ai)
          enddo
        endif               
        if (f_fFed(1:1) /= 'x') then
          do n = 1, n_fed
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'fdFe', trim(nchar)
            call define_hist_field(n_fFed (n,:),vname_in,"umol/m^2/s",tstr2D, tcstr, &
             "dFe flux ice to ocn (cpl)",                     &
             "positive to ocean", c1, c0,                   &
             ns, f_fFed )
          enddo
        endif              
        if (f_fFed_ai (1:1) /= 'x') then
          do n = 1, n_fed 
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'fdFe_ai', trim(nchar)
            call define_hist_field(n_fFed_ai (n,:),vname_in,"umol/m^2/s",tstr2D, tcstr, &
             "dFe flux ice to ocn",                              &
             "weighted by ice area", c1, c0,                         &
             ns, f_fFed_ai )
          enddo
        endif       
        if (f_fFep(1:1) /= 'x') then
          do n = 1, n_fep
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'fpFe', trim(nchar)
            call define_hist_field(n_fFep (n,:),vname_in,"umol/m^2/s",tstr2D, tcstr, &
             "pFe flux ice to ocn (cpl)",                     &
             "positive to ocean", c1, c0,                   &
             ns, f_fFep )
          enddo
        endif
        if (f_fFep_ai (1:1) /= 'x') then
          do n = 1, n_fep 
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'fpFe_ai', trim(nchar)
            call define_hist_field(n_fFep_ai (n,:),vname_in,"umol/m^2/s",tstr2D, tcstr, &
             "pFe flux ice to ocn",                              &
             "weighted by ice area", c1, c0,                         &
             ns, f_fFep_ai )
          enddo
        endif        
        if (f_fSil(1:1) /= 'x') &
          call define_hist_field(n_fSil,"fSil","mmol/m^2/s",tstr2D, tcstr, &
             "silicate flux ice to ocn (cpl)",                     &
             "positive into ocean", c1, c0,          &
             ns, f_fSil)
      
        if (f_fSil_ai(1:1) /= 'x') &
          call define_hist_field(n_fSil_ai,"fSil_ai","mmol/m^2/s",tstr2D, tcstr, &
             "silicate flux ice to ocean",                               &
             "weighted by ice area", c1, c0,                         &
             ns, f_fSil_ai)
     
        if (f_fhum(1:1) /= 'x') &
          call define_hist_field(n_fhum,"fhum","mmol/m^2/s",tstr2D, tcstr, &
             "humic matter (carbon) flux ice to ocn (cpl)",          &
             "positive into ocean", c1, c0,                   &
             ns, f_fhum)
      
        if (f_fhum_ai(1:1) /= 'x') &
          call define_hist_field(n_fhum_ai,"fhum_ai","mmol/m^2/s",tstr2D, tcstr, &
             "humic matter (carbon) flux ice to ocean",              &
             "weighted by ice area", c1, c0,                         &
             ns, f_fhum_ai)

        if (f_fPON(1:1) /= 'x') &
          call define_hist_field(n_fPON,"fPON","mmol/m^2/s",tstr2D, tcstr, &
             "PON flux ice to ocean",                               &
             "positive into ocean", c1, c0,                         &
             ns, f_fPON)

        if (f_fPON_ai(1:1) /= 'x') &
          call define_hist_field(n_fPON_ai,"fPON_ai","mmol/m^2/s",tstr2D, tcstr, &
             "PON flux ice to ocean",                               &
             "weighted by ice area", c1, c0,                         &
             ns, f_fPON_ai)

        if (f_fDMS(1:1) /= 'x') &
          call define_hist_field(n_fDMS,"fDMS","mmol/m^2/s",tstr2D, tcstr, &
             "DMS flux ice to ocean",                               &
             "positive into ocean", c1, c0,                         &
             ns, f_fDMS)

        if (f_fDMS_ai(1:1) /= 'x') &
          call define_hist_field(n_fDMS_ai,"fDMS_ai","mmol/m^2/s",tstr2D, tcstr, &
             "DMS flux ice to ocean",                               &
             "weighted by ice area", c1, c0,                         &
             ns, f_fDMS_ai)

        if (f_fDMSPd(1:1) /= 'x') &
          call define_hist_field(n_fDMSPd,"fDMSPd","mmol/m^2/s",tstr2D, tcstr, &
             "DMSPd flux ice to ocean",                               &
             "positive into ocean", c1, c0,                         &
             ns, f_fDMSPd)

        if (f_fDMSPd_ai(1:1) /= 'x') &
          call define_hist_field(n_fDMSPd_ai,"fDMSPd_ai","mmol/m^2/s",tstr2D, tcstr, &
             "DMSPd flux ice to ocean",                               &
             "weighted by ice area", c1, c0,                         &
             ns, f_fDMSPd_ai)

        if (f_fDMSPp(1:1) /= 'x') &
          call define_hist_field(n_fDMSPp,"fDMSPp","mmol/m^2/s",tstr2D, tcstr, &
             "DMSPp flux ice to ocean",                               &
             "positive into ocean", c1, c0,                         &
             ns, f_fDMSPp)

        if (f_fDMSPp_ai(1:1) /= 'x') &
          call define_hist_field(n_fDMSPp_ai,"fDMSPp_ai","mmol/m^2/s",tstr2D, tcstr, &
             "DMSPp flux ice to ocean",                               &
             "weighted by ice area", c1, c0,                         &
             ns, f_fDMSPp_ai)

        if (f_PPnet(1:1) /= 'x') &
          call define_hist_field(n_PPnet,"PP_net","mg C/m^2/d",tstr2D, tcstr, &
             "Net Primary Production",            &
             "weighted by ice area", c1, c0,       &
             ns, f_PPnet)

        if (f_grownet(1:1) /= 'x') &
          call define_hist_field(n_grownet,"grow_net","m/d",tstr2D, tcstr, &
             "Net specific growth",                     &
             "weighted by brine or skl volume ", c1, c0,       &
             ns, f_grownet)

        if (f_upNO(1:1) /= 'x') & 
          call define_hist_field(n_upNO,"upNO","mmol/m^2/d",tstr2D, tcstr, &
             "Tot algal Nit uptake rate",                  &
             "weighted by ice area", c1, c0, &
             ns, f_upNO)

        if (f_upNH(1:1) /= 'x') &  
          call define_hist_field(n_upNH,"upNH","mmol/m^2/d",tstr2D, tcstr, &
             "Tot algal Am uptake rate",                  &
             "weighted by ice area", c1, c0,&
             ns, f_upNH)

        ! vertical biogeochemistry  
      if (z_tracers) then

        if (f_fzaero(1:1) /= 'x') then
          do n = 1, n_zaero
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'fzaero', trim(nchar)
            call define_hist_field(n_fzaero(n,:),vname_in,"kg/m^2/s",tstr2D, tcstr, &
             "z aerosol flux ice to ocn (cpl)",                     &
             "positive to ocean", c1, c0,                   &
             ns, f_fzaero)
          enddo
        endif
        if (f_fzaero_ai(1:1) /= 'x') then
          do n = 1, n_zaero
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'fzaero_ai', trim(nchar)
            call define_hist_field(n_fzaero_ai(n,:),vname_in,"kg/m^2/s",tstr2D, tcstr, &
             "z aerosol flux ice to ocn",                     &
             "weighted by ice area", c1, c0,                   &
             ns, f_fzaero_ai)
          enddo
        endif
        if (f_algalpeak(1:1) /= 'x') then
          do n = 1, n_algae
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'peak_loc', trim(nchar)
            call define_hist_field(n_algalpeak(n,:),vname_in,"1",tstr2D, tcstr, &
             "zgrid level of peak chla",            &
             "0 if no distinct peak", c1, c0,       &
             ns, f_algalpeak)
          enddo
        endif
        if (f_peakval(1:1) /= 'x') then
          do n = 1, n_algae
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'peak_val', trim(nchar)
            call define_hist_field(n_peakval(n,:),vname_in,"mg/m^3",tstr2D, tcstr, &
             "Bulk concentration of peak chla",            &
             "at peak_loc", c1, c0,       &
             ns, f_peakval)
          enddo
        endif
        if (f_zaeronet(1:1) /= 'x')  then
          do n = 1, n_zaero
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'zaero_net', trim(nchar)
             call define_hist_field(n_zaeronet(n,:),vname_in,"kg/m^2",tstr2D, tcstr, &
             "Net z aerosol",                     &
             "weighted by ice area ", c1, c0,       &
             ns, f_zaeronet)
          enddo
        endif  !f_zaeronet
        if (f_chlnet(1:1) /= 'x')  then
          do n = 1, n_algae
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'chl_net', trim(nchar)
             call define_hist_field(n_chlnet(n,:),vname_in,"mg chl/m^2",tstr2D, tcstr, &
             "Net chlorophyll",                     &
             "weighted by ice area ", c1, c0,       &
             ns, f_chlnet)
          enddo
        endif  !f_chlnet
        if (f_Nnet(1:1) /= 'x') then
          do n = 1, n_algae
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'algalN_net', trim(nchar)
            call define_hist_field(n_Nnet(n,:),vname_in,"mmol/m^2",tstr2D, tcstr, &
             "Net algal nitrogen",                     &
             "weighted by ice area ", c1, c0,       &
             ns, f_Nnet)
          enddo
        endif !f_Nnet

        if (f_Cnet(1:1) /= 'x') then
          do n = 1, n_algae
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'algalC_net', trim(nchar)
            call define_hist_field(n_Cnet(n,:),vname_in,"mmol/m^2",tstr2D, tcstr, &
             "Net algal carbon",                     &
             "weighted by ice area ", c1, c0,       &
             ns, f_Cnet)
          enddo
        endif !f_Cnet
        if (f_DOCnet(1:1) /= 'x') then
          do n = 1, n_doc
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'DOC_net', trim(nchar)
            call define_hist_field(n_DOCnet(n,:),vname_in,"mmol/m^2",tstr2D, tcstr, &
             "Net DOC",                     &
             "weighted by ice area ", c1, c0,       &
             ns, f_DOCnet)
          enddo
        endif !f_DOCnet
        if (f_DICnet(1:1) /= 'x') then
          do n = 1, n_dic
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'DIC_net', trim(nchar)
            call define_hist_field(n_DICnet(n,:),vname_in,"mmol/m^2",tstr2D, tcstr, &
             "Net DIC",                     &
             "weighted by ice area ", c1, c0,       &
             ns, f_DICnet)
          enddo
        endif !f_DICnet
        if (f_DONnet(1:1) /= 'x') then
          do n = 1, n_don
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'DON_net', trim(nchar)
            call define_hist_field(n_DONnet(n,:),vname_in,"mmol/m^2",tstr2D, tcstr, &
             "Net DON",                     &
             "weighted by ice area ", c1, c0,       &
             ns, f_DONnet)
          enddo
        endif !f_DONnet
        if (f_Fednet (1:1) /= 'x') then
          do n = 1, n_fed 
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'dFe_net', trim(nchar)
            call define_hist_field(n_Fednet (n,:),vname_in,"umol/m^2",tstr2D, tcstr, &
             "Net dFe",                     &
             "weighted by ice area ", c1, c0,       &
             ns, f_Fednet )
          enddo
        endif !f_Fednet 
        if (f_Fepnet (1:1) /= 'x') then
          do n = 1, n_fep 
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'pFe_net', trim(nchar)
            call define_hist_field(n_Fepnet (n,:),vname_in,"umol/m^2",tstr2D, tcstr, &
             "Net pFe",                     &
             "weighted by ice area ", c1, c0,       &
             ns, f_Fepnet )
          enddo
        endif !f_Fepnet 
        if (f_Nitnet(1:1) /= 'x') &
          call define_hist_field(n_Nitnet,"Nit_net","mmol/m^2",tstr2D, tcstr, &
             "Net Nitrate",                     &
             "weighted by ice area", c1, c0,       &
             ns, f_Nitnet)
        if (f_Amnet(1:1) /= 'x') &
          call define_hist_field(n_Amnet,"Am_net","mmol/m^2",tstr2D, tcstr, &
             "Net Ammonium",                     &
             "weighted by ice area", c1, c0,       &
             ns, f_Amnet)
        if (f_Silnet(1:1) /= 'x') &
          call define_hist_field(n_Silnet,"Sil_net","mmol/m^2",tstr2D, tcstr, &
             "Net Silicate",                     &
             "weighted by ice area", c1, c0,       &
             ns, f_Silnet)
        if (f_humnet(1:1) /= 'x') &
          call define_hist_field(n_humnet,"hum_net","mmol/m^2",tstr2D, tcstr, &
             "Net humic material (carbon)",               &
             "weighted by ice area", c1, c0,       &
             ns, f_humnet) 
        if (f_DMSPpnet(1:1) /= 'x') &
          call define_hist_field(n_DMSPpnet,"DMSPp_net","mmol/m^2",tstr2D, tcstr, &
             "Net DMSPp",                     &
             "weighted by ice area", c1, c0,       &
             ns, f_DMSPpnet)
        if (f_DMSPdnet(1:1) /= 'x') &
          call define_hist_field(n_DMSPdnet,"DMSPd_net","mmol/m^2",tstr2D, tcstr, &
             "Net DMSPd",                     &
             "weighted by ice area", c1, c0,       &
             ns, f_DMSPdnet)
        if (f_DMSnet(1:1) /= 'x') &
          call define_hist_field(n_DMSnet,"DMS_net","mmol/m^2",tstr2D, tcstr, &
             "Net DMS",                     &
             "weighted by ice area", c1, c0,       &
             ns, f_DMSnet)
        if (f_PONnet(1:1) /= 'x') &
          call define_hist_field(n_PONnet,"PON_net","mmol/m^2",tstr2D, tcstr, &
             "Net Nitrate if no reactions",                     &
             "weighted by ice area", c1, c0,       &
             ns, f_PONnet)

        if (f_zaerosnow(1:1) /= 'x')  then
          do n = 1, n_zaero 
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'zaero_snow', trim(nchar)
             call define_hist_field(n_zaerosnow(n,:),vname_in,"kg/m^2",tstr2D, tcstr, &
             "Snow z aerosol",                     &
             "weighted by ice area ", c1, c0,       &
             ns, f_zaerosnow)
          enddo
        endif  !f_chlnet
        if (f_chlsnow(1:1) /= 'x')  then
          do n = 1, n_algae
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'chl_snow', trim(nchar)
             call define_hist_field(n_chlsnow(n,:),vname_in,"mg chl/m^2",tstr2D, tcstr, &
             "Snow chlorophyll",                     &
             "weighted by ice area ", c1, c0,       &
             ns, f_chlsnow)
          enddo
        endif  !f_chlnet
        if (f_Nsnow(1:1) /= 'x') then
          do n = 1, n_algae
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'algalN_snow', trim(nchar)
            call define_hist_field(n_Nsnow(n,:),vname_in,"mmol/m^2",tstr2D, tcstr, &
             "Snow algal nitrogen",                     &
             "weighted by ice area ", c1, c0,       &
             ns, f_Nsnow)
          enddo
        endif !f_Nsnow
        if (f_Csnow(1:1) /= 'x') then
          do n = 1, n_algae
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'algalC_snow', trim(nchar)
            call define_hist_field(n_Csnow(n,:),vname_in,"mmol/m^2",tstr2D, tcstr, &
             "Snow algal carbon",                     &
             "weighted by ice area ", c1, c0,       &
             ns, f_Csnow)
          enddo
        endif !f_Csnow
        if (f_DOCsnow(1:1) /= 'x') then
          do n = 1, n_doc
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'DOC_snow', trim(nchar)
            call define_hist_field(n_DOCsnow(n,:),vname_in,"mmol/m^2",tstr2D, tcstr, &
             "Snow DOC",                     &
             "weighted by ice area ", c1, c0,       &
             ns, f_DOCsnow)
          enddo
        endif !f_DOCsnow
        if (f_DICsnow(1:1) /= 'x') then
          do n = 1, n_dic
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'DIC_snow', trim(nchar)
            call define_hist_field(n_DICsnow(n,:),vname_in,"mmol/m^2",tstr2D, tcstr, &
             "Snow DIC",                     &
             "weighted by ice area ", c1, c0,       &
             ns, f_DICsnow)
          enddo
        endif !f_DICsnow
        if (f_DONsnow(1:1) /= 'x') then
          do n = 1, n_don
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'DON_snow', trim(nchar)
            call define_hist_field(n_DONsnow(n,:),vname_in,"mmol/m^2",tstr2D, tcstr, &
             "Snow DON",                     &
             "weighted by ice area ", c1, c0,       &
             ns, f_DONsnow)
          enddo
        endif !f_DONsnow
        if (f_Fedsnow (1:1) /= 'x') then
          do n = 1, n_fed 
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'dFe_snow', trim(nchar)
            call define_hist_field(n_Fedsnow (n,:),vname_in,"umol/m^2",tstr2D, tcstr, &
             "Snow dFe",                     &
             "weighted by ice area ", c1, c0,       &
             ns, f_Fedsnow )
          enddo
        endif !f_Fedsnow 
        if (f_Fepsnow (1:1) /= 'x') then
          do n = 1, n_fed 
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'pFe_snow', trim(nchar)
            call define_hist_field(n_Fepsnow (n,:),vname_in,"umol/m^2",tstr2D, tcstr, &
             "Snow pFe",                     &
             "weighted by ice area ", c1, c0,       &
             ns, f_Fepsnow )
          enddo
        endif !f_Fepsnow 
        if (f_Nitsnow(1:1) /= 'x') &
          call define_hist_field(n_Nitsnow,"Nit_snow","mmol/m^2",tstr2D, tcstr, &
             "Snow Nitrate",                     &
             "weighted by ice area", c1, c0,       &
             ns, f_Nitsnow)
        if (f_Amsnow(1:1) /= 'x') &
          call define_hist_field(n_Amsnow,"Am_snow","mmol/m^2",tstr2D, tcstr, &
             "Snow Ammonium",                     &
             "weighted by ice area", c1, c0,       &
             ns, f_Amsnow)
        if (f_Silsnow(1:1) /= 'x') &
          call define_hist_field(n_Silsnow,"Sil_snow","mmol/m^2",tstr2D, tcstr, &
             "Snow Silicate",                     &
             "weighted by ice area", c1, c0,       &
             ns, f_Silsnow)
        if (f_humsnow(1:1) /= 'x') &
          call define_hist_field(n_humsnow,"hum_snow","mmol/m^2",tstr2D, tcstr, &
             "Snow humic material (carbon)",               &
             "weighted by ice area", c1, c0,       &
             ns, f_humsnow) 
        if (f_DMSPpsnow(1:1) /= 'x') &
          call define_hist_field(n_DMSPpsnow,"DMSPp_snow","mmol/m^2",tstr2D, tcstr, &
             "Snow DMSPp",                     &
             "weighted by ice area", c1, c0,       &
             ns, f_DMSPpsnow)
        if (f_DMSPdsnow(1:1) /= 'x') &
          call define_hist_field(n_DMSPdsnow,"DMSPd_snow","mmol/m^2",tstr2D, tcstr, &
             "Snow DMSPd",                     &
             "weighted by ice area", c1, c0,       &
             ns, f_DMSPdsnow)
        if (f_DMSsnow(1:1) /= 'x') &
          call define_hist_field(n_DMSsnow,"DMS_snow","mmol/m^2",tstr2D, tcstr, &
             "Snow DMS",                     &
             "weighted by ice area", c1, c0,       &
             ns, f_DMSsnow)
        if (f_PONsnow(1:1) /= 'x') &
          call define_hist_field(n_PONsnow,"PON_snow","mmol/m^2",tstr2D, tcstr, &
             "Snow Nitrate if no reactions",                     &
             "weighted by ice area", c1, c0,       &
             ns, f_PONsnow)

        if (f_zaerofrac(1:1) /= 'x')  then
          do n = 1, n_zaero 
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'zaero_frac', trim(nchar)
             call define_hist_field(n_zaerofrac(n,:),vname_in,"1",tstr2D, tcstr, &
             "Mobile frac z aerosol",                     &
             "averaged over depth ", c1, c0,       &
             ns, f_zaerofrac)
          enddo
        endif  !f_chlnet
        if (f_chlfrac(1:1) /= 'x')  then
          do n = 1, n_algae
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'chl_frac', trim(nchar)
             call define_hist_field(n_chlfrac(n,:),vname_in,"mg chl/m^2",tstr2D, tcstr, &
             "Mobile frac chlorophyll",                     &
             "averaged over depth ", c1, c0,       &
             ns, f_chlfrac)
          enddo
        endif  !f_chlnet
        if (f_Nfrac(1:1) /= 'x') then
          do n = 1, n_algae
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'algalN_frac', trim(nchar)
            call define_hist_field(n_Nfrac(n,:),vname_in,"1",tstr2D, tcstr, &
             "Mobile frac algal nitrogen",                     &
             "averaged over depth ", c1, c0,       &
             ns, f_Nfrac)
          enddo
        endif !f_Nfrac
        if (f_DOCfrac(1:1) /= 'x') then
          do n = 1, n_doc
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'DOC_frac', trim(nchar)
            call define_hist_field(n_DOCfrac(n,:),vname_in,"1",tstr2D, tcstr, &
             "Mobile frac DOC",                     &
             "averaged over depth ", c1, c0,       &
             ns, f_DOCfrac)
          enddo
        endif !f_DOCfrac
        if (f_DICfrac(1:1) /= 'x') then
          do n = 1, n_dic
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'DIC_frac', trim(nchar)
            call define_hist_field(n_DICfrac(n,:),vname_in,"1",tstr2D, tcstr, &
             "Mobile frac DIC",                     &
             "averaged over depth ", c1, c0,       &
             ns, f_DICfrac)
          enddo
        endif !f_DICfrac
        if (f_DONfrac(1:1) /= 'x') then
          do n = 1, n_don
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'DON_frac', trim(nchar)
            call define_hist_field(n_DONfrac(n,:),vname_in,"1",tstr2D, tcstr, &
             "Mobile frac DON",                     &
             "averaged over depth ", c1, c0,       &
             ns, f_DONfrac)
         enddo
        endif !f_DONfrac
        if (f_Fedfrac (1:1) /= 'x') then
          do n = 1, n_fed 
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'dFe_frac', trim(nchar)
            call define_hist_field(n_Fedfrac (n,:),vname_in,"1",tstr2D, tcstr, &
             "Mobile frac dFe",                     &
             "averaged over depth ", c1, c0,       &
             ns, f_Fedfrac )
          enddo
        endif !f_Fedfrac 
        if (f_Fepfrac (1:1) /= 'x') then
          do n = 1, n_fep 
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'pFe_frac', trim(nchar)
            call define_hist_field(n_Fepfrac (n,:),vname_in,"1",tstr2D, tcstr, &
             "Mobile frac pFe",                     &
             "averaged over depth ", c1, c0,       &
             ns, f_Fepfrac )
          enddo
        endif !f_Fepfrac 
        if (f_Nitfrac(1:1) /= 'x') &
          call define_hist_field(n_Nitfrac,"Nit_frac","1",tstr2D, tcstr, &
             "Mobile frac Nitrate",                     &
             "averaged over depth", c1, c0,       &
             ns, f_Nitfrac)
        if (f_Amfrac(1:1) /= 'x') &
          call define_hist_field(n_Amfrac,"Am_frac","1",tstr2D, tcstr, &
             "Mobile frac Ammonium",                     &
             "averaged over depth", c1, c0,       &
             ns, f_Amfrac)
        if (f_Silfrac(1:1) /= 'x') &
          call define_hist_field(n_Silfrac,"Sil_frac","1",tstr2D, tcstr, &
             "Mobile frac Silicate",                     &
             "averaged over depth", c1, c0,       &
             ns, f_Silfrac)
        if (f_humfrac(1:1) /= 'x') &
          call define_hist_field(n_humfrac,"hum_frac","1",tstr2D, tcstr, &
             "Mobile frac humic material",               &
             "averaged over depth", c1, c0,       &
             ns, f_humfrac) 
        if (f_DMSPpfrac(1:1) /= 'x') &
          call define_hist_field(n_DMSPpfrac,"DMSPp_frac","1",tstr2D, tcstr, &
             "Mobile frac DMSPp",                     &
             "averaged over depth", c1, c0,       &
             ns, f_DMSPpfrac)
        if (f_DMSPdfrac(1:1) /= 'x') &
          call define_hist_field(n_DMSPdfrac,"DMSPd_frac","1",tstr2D, tcstr, &
             "Mobile frac DMSPd",                     &
             "averaged over depth", c1, c0,       &
             ns, f_DMSPdfrac)
        if (f_DMSfrac(1:1) /= 'x') &
          call define_hist_field(n_DMSfrac,"DMS_frac","1",tstr2D, tcstr, &
             "Mobile frac DMS",                     &
             "averaged over depth", c1, c0,       &
             ns, f_DMSfrac)
        if (f_PONfrac(1:1) /= 'x') &
          call define_hist_field(n_PONfrac,"PON_frac","1",tstr2D, tcstr, &
             "Mobile frac Nitrate if no reactions",                     &
             "averaged over depth", c1, c0,       &
             ns, f_PONfrac)

      endif ! z_tracers

      ! brine
      if (f_hbri(1:1) /= 'x') &
         call define_hist_field(n_hbri,"hbrine","m",tstr2D, tcstr, &
             "Brine height",                     &
             "distance from ice bottom to brine surface", c1, c0,       &
             ns, f_hbri)

    enddo ! nstreams

    endif ! tr_aero, etc 
      
      end subroutine init_hist_bgc_2D

!=======================================================================

      subroutine init_hist_bgc_3Dc

      use ice_calendar, only: nstreams
      use ice_history_shared, only: tstr3Dc, tcstr, define_hist_field

      integer (kind=int_kind) :: ns
      logical (kind=log_kind) :: tr_brine
      character(len=*), parameter :: subname = '(init_hist_bgc_3Dc)'

      call icepack_query_tracer_flags(tr_brine_out=tr_brine)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

    if (tr_brine) then

      ! 3D (category) variables must be looped separately
      do ns = 1, nstreams
        if (f_fbri(1:1) /= 'x') &
         call define_hist_field(n_fbri,"fbrine","1",tstr3Dc, tcstr, &
             "brine tracer fraction of ice volume, cat",             &
             "none", c1, c0,       &
             ns, f_fbri)
      enddo ! ns

    endif

      end subroutine init_hist_bgc_3Dc

!=======================================================================

      subroutine init_hist_bgc_3Db

      use ice_calendar, only: nstreams
      use ice_history_shared, only: tstr3Db, tcstr, define_hist_field

      integer (kind=int_kind) :: ns
      real (kind=dbl_kind) :: secday
      logical (kind=log_kind) :: solve_zsal, z_tracers
      character(len=*), parameter :: subname = '(init_hist_bgc_3Db)'
      
      ! biology vertical grid

      call icepack_query_parameters(secday_out=secday)
      call icepack_query_parameters( &
          solve_zsal_out=solve_zsal, z_tracers_out=z_tracers)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

    if (z_tracers .or. solve_zsal) then

      do ns = 1, nstreams
 
         if (f_bTin(1:1) /= 'x') &
            call define_hist_field(n_bTin,"bTizn","C",tstr3Db, tcstr, &
                "ice internal temperatures on bio grid", &
                "interpolated to bio grid", c1, c0,      &
                ns, f_bTin)

         if (f_bphi(1:1) /= 'x') &
            call define_hist_field(n_bphi,"bphizn","%",tstr3Db, tcstr, &
                "porosity", "brine volume fraction", c100, c0, &
                ns, f_bphi)
         
         if (f_iDi(1:1) /= 'x') &  
             call define_hist_field(n_iDi,"iDin","m^2/d",tstr3Db, tcstr, &
                "interface diffusivity", "on bio interface grid", secday, c0, &
                ns, f_iDi)
      
         if (f_iki(1:1) /= 'x') & 
            call define_hist_field(n_iki,"ikin","mm^2",tstr3Db, tcstr, &
                "permeability", "on bio interface grid", 1.0e6_dbl_kind, c0, &
                ns, f_iki)
 
         if (f_bgc_S(1:1) /= 'x') &
            call define_hist_field(n_bgc_S,"bgc_S","ppt",tstr3Db, tcstr, &
                "bulk salinity", "on bio grid", c1, c0, &
                ns, f_bgc_S)
      
         if (f_zfswin(1:1) /= 'x') &
            call define_hist_field(n_zfswin,"zfswin","W/m^2",tstr3Db, tcstr, &
                "internal ice PAR", "on bio interface grid", c1, c0, &
                ns, f_zfswin)
    
      enddo  ! ns

    endif  ! z_tracers or solve_zsal

      end subroutine init_hist_bgc_3Db

!=======================================================================
!
! write average ice quantities or snapshots

      subroutine accum_hist_bgc (iblk)

      use ice_arrays_column, only: ocean_bio, &
          grow_net, PP_net, upNO, upNH, ice_bio_net, snow_bio_net, &
          hbri, bTiz, bphi, zfswin, iDi, iki, zsal_tot, fzsal, fzsal_g, &
          R_C2N, R_chl2N
      use ice_blocks, only: block, get_block, nx_block, ny_block
      use ice_domain, only: blocks_ice
      use ice_domain_size, only: nblyr
      use ice_flux, only: sss
      use ice_flux_bgc, only: faero_atm, faero_ocn, flux_bio, flux_bio_ai, &
          fzsal_ai, fzsal_g_ai
      use ice_history_shared, only: n2D, a2D, a3Dc, &         
          n3Dzcum, n3Dbcum, a3Db, a3Da, &    
          ncat_hist, accum_hist_field, nzblyr, nzalyr
      use ice_state, only: trcrn, trcr, aicen, aice, vicen

      integer (kind=int_kind), intent(in) :: &
           iblk                 ! block index

      ! local variables

      integer (kind=int_kind) :: &
         i, j, n, k,        & ! loop indices  
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+4) :: & 
         workz, workz2

      real (kind=dbl_kind) :: & 
         maxv, rhos, rhoi, rhow, puny, sk_l

      real (kind=dbl_kind), dimension (nblyr+1) :: & 
         workv

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: & 
         workni, worknj

      integer (kind=int_kind), dimension (1) :: & 
         worki
      integer (kind=int_kind) :: & 
         workii

      logical (kind=log_kind) :: &
         skl_bgc, z_tracers, tr_aero, tr_brine, solve_zsal

      integer(kind=int_kind) :: nt_aero, nt_fbri, &
         nt_bgc_Nit,   nt_bgc_Am,   nt_bgc_Sil, nt_bgc_DMSPp, &
         nt_bgc_DMSPd,  nt_bgc_DMS, nt_bgc_PON, nt_bgc_S, &
         nt_zbgc_frac,  nlt_chl_sw, &
         nlt_bgc_Nit,   nlt_bgc_Am, nlt_bgc_Sil,   &
         nlt_bgc_DMS,   nlt_bgc_PON,  &
         nlt_bgc_DMSPp, nlt_bgc_DMSPd, &
         nt_bgc_hum,    nlt_bgc_hum

      integer (kind=int_kind), dimension(icepack_max_aero) :: &
         nlt_zaero_sw       ! points to aerosol in trcrn_sw
  
      integer (kind=int_kind), dimension(icepack_max_algae) :: &
         nt_bgc_N, nlt_bgc_N, &     ! algae 
         nt_bgc_C, nlt_bgc_C, &     !
         nt_bgc_chl, nlt_bgc_chl    !

      integer (kind=int_kind), dimension(icepack_max_doc) :: &
         nt_bgc_DOC, nlt_bgc_DOC    ! disolved organic carbon

      integer (kind=int_kind), dimension(icepack_max_don) :: &
         nt_bgc_DON, nlt_bgc_DON    !

      integer (kind=int_kind), dimension(icepack_max_dic) :: &
         nt_bgc_DIC, nlt_bgc_DIC    ! disolved inorganic carbon

      integer (kind=int_kind), dimension(icepack_max_fe) :: &
         nt_bgc_Fed, nlt_bgc_Fed, & !
         nt_bgc_Fep, nlt_bgc_Fep    !

      integer (kind=int_kind), dimension(icepack_max_aero) :: &
         nt_zaero, nlt_zaero        ! non-reacting layer aerosols

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(accum_hist_bgc)'

      call icepack_query_parameters(rhos_out=rhos, rhoi_out=rhoi, &
         rhow_out=rhow, puny_out=puny, sk_l_out=sk_l)
      call icepack_query_tracer_flags( &
         tr_aero_out=tr_aero, tr_brine_out=tr_brine)
      call icepack_query_parameters(skl_bgc_out=skl_bgc, &
         solve_zsal_out=solve_zsal, z_tracers_out=z_tracers)
      call icepack_query_tracer_indices( nt_aero_out=nt_aero, &
         nt_fbri_out=nt_fbri, nt_bgc_DOC_out=nt_bgc_DOC, &
         nt_zaero_out=nt_zaero, nt_bgc_DIC_out=nt_bgc_DIC, &
         nt_bgc_DON_out=nt_bgc_DON, nt_bgc_N_out=nt_bgc_N, &
         nt_bgc_C_out=nt_bgc_C, nt_bgc_Am_out=nt_bgc_Am, &
         nt_bgc_chl_out=nt_bgc_chl, nt_bgc_Nit_out=nt_bgc_Nit, &
         nt_bgc_Sil_out=nt_bgc_Sil, nt_bgc_DMSPp_out=nt_bgc_DMSPp, &
         nt_bgc_DMSPd_out=nt_bgc_DMSPd, nt_bgc_DMS_out=nt_bgc_DMS, &
         nt_bgc_PON_out=nt_bgc_PON, &
         nt_bgc_S_out=nt_bgc_S, nt_bgc_Fed_out=nt_bgc_Fed, &
         nt_bgc_Fep_out=nt_bgc_Fep, nt_zbgc_frac_out=nt_zbgc_frac, &
         nlt_zaero_sw_out=nlt_zaero_sw, nlt_chl_sw_out=nlt_chl_sw, &
         nlt_bgc_Nit_out=nlt_bgc_Nit,   nlt_bgc_Am_out=nlt_bgc_Am,   &
         nlt_bgc_Sil_out=nlt_bgc_Sil,   &
         nlt_bgc_DMS_out=nlt_bgc_DMS,   nlt_bgc_PON_out=nlt_bgc_PON,  &
         nlt_bgc_N_out=nlt_bgc_N,     nlt_bgc_C_out=nlt_bgc_C,    &
         nlt_bgc_chl_out=nlt_bgc_chl, nlt_bgc_DIC_out=nlt_bgc_DIC,   &
         nlt_bgc_DOC_out=nlt_bgc_DOC,   nlt_bgc_DON_out=nlt_bgc_DON,  &
         nlt_zaero_out=nlt_zaero,   nlt_bgc_DMSPp_out=nlt_bgc_DMSPp, &
         nlt_bgc_DMSPd_out=nlt_bgc_DMSPd, &
         nlt_bgc_Fed_out=nlt_bgc_Fed,   nlt_bgc_Fep_out=nlt_bgc_Fep,  &
         nt_bgc_hum_out=nt_bgc_hum,    nlt_bgc_hum_out=nlt_bgc_hum)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)
          
      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi

      !---------------------------------------------------------------
      ! increment field
      !---------------------------------------------------------------

    if (tr_aero .or. tr_brine .or. solve_zsal .or. skl_bgc) then
      ! 2d bgc fields


      ! zsalinity
      if (f_fzsal  (1:1) /= 'x') &  
         call accum_hist_field(n_fzsal,     iblk, fzsal(:,:,iblk), a2D)
      if (f_fzsal_ai(1:1)/= 'x') & 
         call accum_hist_field(n_fzsal_ai,  iblk, fzsal_ai(:,:,iblk), a2D)
      if (f_fzsal_g  (1:1) /= 'x') & 
         call accum_hist_field(n_fzsal_g,   iblk, fzsal_g(:,:,iblk), a2D)
      if (f_fzsal_g_ai(1:1)/= 'x') & 
         call accum_hist_field(n_fzsal_g_ai,iblk, fzsal_g_ai(:,:,iblk), a2D)
      if (f_zsal  (1:1) /= 'x') &   
         call accum_hist_field(n_zsal,      iblk, zsal_tot(:,:,iblk), a2D)

      ! Aerosols
      if (f_faero_atm(1:1) /= 'x') then
         do n=1,n_aero
            call accum_hist_field(n_faero_atm(n,:),iblk, &
                                    faero_atm(:,:,n,iblk), a2D)
         enddo
      endif
      if (f_faero_ocn(1:1) /= 'x') then
         do n=1,n_aero
            call accum_hist_field(n_faero_ocn(n,:),iblk, &
                                    faero_ocn(:,:,n,iblk), a2D)
         enddo
      endif
      if (f_aero(1:1) /= 'x') then
         do n=1,n_aero
            call accum_hist_field(n_aerosn1(n,:), iblk, &
                               trcr(:,:,nt_aero  +4*(n-1),iblk)/rhos, a2D)
            call accum_hist_field(n_aerosn2(n,:), iblk, &
                               trcr(:,:,nt_aero+1+4*(n-1),iblk)/rhos, a2D)
            call accum_hist_field(n_aeroic1(n,:), iblk, &
                               trcr(:,:,nt_aero+2+4*(n-1),iblk)/rhoi, a2D)
            call accum_hist_field(n_aeroic2(n,:), iblk, &
                               trcr(:,:,nt_aero+3+4*(n-1),iblk)/rhoi, a2D)
         enddo
      endif

    if (skl_bgc) then

      ! skeletal layer bgc

      if (f_bgc_N(1:1)/= 'x') then
         do n=1,n_algae
            call accum_hist_field(n_bgc_N(n,:), iblk, &
                    sk_l*trcr(:,:,nt_bgc_N(n),   iblk), a2D)
         enddo
      endif
      if (f_bgc_C(1:1)/= 'x') then
         do n=1,n_algae
            call accum_hist_field(n_bgc_C(n,:),        iblk, &
                    sk_l*R_C2N(n)*trcr(:,:,nt_bgc_N(n), iblk), a2D)
         enddo
      endif
      if (f_bgc_DOC(1:1)/= 'x') then
         do n=1,n_doc
            call accum_hist_field(n_bgc_DOC(n,:), iblk, &
                    sk_l*trcr(:,:,nt_bgc_DOC(n),   iblk), a2D)
         enddo
      endif
      if (f_bgc_DIC(1:1)/= 'x') then
         do n=1,n_dic
            call accum_hist_field(n_bgc_DIC(n,:), iblk, &
                    sk_l*trcr(:,:,nt_bgc_DIC(n),   iblk), a2D)
         enddo
      endif
      if (f_bgc_DON(1:1)/= 'x') then
         do n=1,n_don
            call accum_hist_field(n_bgc_DON(n,:), iblk, &
                    sk_l*trcr(:,:,nt_bgc_DON(n),   iblk), a2D)
         enddo
      endif
      if (f_bgc_Fed (1:1)/= 'x') then
         do n=1,n_fed 
            call accum_hist_field(n_bgc_Fed (n,:), iblk, &
                    sk_l*trcr(:,:,nt_bgc_Fed (n),   iblk), a2D)
         enddo
      endif
      if (f_bgc_Fep (1:1)/= 'x') then
         do n=1,n_fep 
            call accum_hist_field(n_bgc_Fep (n,:), iblk, &
                    sk_l*trcr(:,:,nt_bgc_Fep (n),   iblk), a2D)
         enddo
      endif
      if (f_bgc_chl(1:1)/= 'x') then
         do n=1,n_algae
            call accum_hist_field(n_bgc_chl(n,:), iblk, &
                    sk_l*trcr(:,:,nt_bgc_chl(n),   iblk), a2D)
         enddo
      endif
      if (f_bgc_Nit(1:1)/= 'x') &
         call accum_hist_field(n_bgc_Nit,  iblk, &
                sk_l*trcr(:,:,nt_bgc_Nit,  iblk), a2D)  
      if (f_bgc_Am(1:1)/= 'x') &
         call accum_hist_field(n_bgc_Am,   iblk, &
                sk_l*trcr(:,:,nt_bgc_Am,   iblk), a2D)  
      if (f_bgc_Sil(1:1)/= 'x') &
         call accum_hist_field(n_bgc_Sil,  iblk, &
                sk_l*trcr(:,:,nt_bgc_Sil,  iblk), a2D)   
      if (f_bgc_hum(1:1)/= 'x') &
         call accum_hist_field(n_bgc_hum,  iblk, &
                sk_l*trcr(:,:,nt_bgc_hum,  iblk), a2D)  
      if (f_bgc_PON(1:1)/= 'x') &
         call accum_hist_field(n_bgc_PON,  iblk, &
                sk_l*trcr(:,:,nt_bgc_PON,  iblk), a2D)  
      if (f_bgc_DMSPp(1:1)/= 'x') &
         call accum_hist_field(n_bgc_DMSPp,iblk, &
                sk_l*trcr(:,:,nt_bgc_DMSPp,iblk), a2D)  
      if (f_bgc_DMSPd(1:1)/= 'x') &
         call accum_hist_field(n_bgc_DMSPd,iblk, &
                sk_l*trcr(:,:,nt_bgc_DMSPd,iblk), a2D)  
      if (f_bgc_DMS(1:1)/= 'x') &
         call accum_hist_field(n_bgc_DMS,  iblk, &
                sk_l*trcr(:,:,nt_bgc_DMS,  iblk), a2D)   

    endif  !skl_bgc 

      ! skeletal layer and vertical bgc 

      if (f_bgc_DOC_ml(1:1)/= 'x') then
         do n=1,n_doc
            call accum_hist_field(n_bgc_DOC_ml(n,:), iblk, &
                   ocean_bio(:,:,nlt_bgc_DOC(n),      iblk), a2D)
         enddo
      endif
      if (f_bgc_DIC_ml(1:1)/= 'x') then
         do n=1,n_dic
            call accum_hist_field(n_bgc_DIC_ml(n,:), iblk, &
                   ocean_bio(:,:,nlt_bgc_DIC(n),      iblk), a2D)
         enddo
      endif
      if (f_bgc_DON_ml(1:1)/= 'x') then
         do n=1,n_don
            call accum_hist_field(n_bgc_DON_ml(n,:), iblk, &
                   ocean_bio(:,:,nlt_bgc_DON(n),      iblk), a2D)
         enddo
      endif
      if (f_bgc_Fed_ml (1:1)/= 'x') then
         do n=1,n_fed 
            call accum_hist_field(n_bgc_Fed_ml (n,:), iblk, &
                   ocean_bio(:,:,nlt_bgc_Fed (n),      iblk), a2D)
         enddo
      endif
      if (f_bgc_Fep_ml (1:1)/= 'x') then
         do n=1,n_fep 
            call accum_hist_field(n_bgc_Fep_ml (n,:), iblk, &
                   ocean_bio(:,:,nlt_bgc_Fep (n),      iblk), a2D)
         enddo
      endif
      if (f_bgc_N_ml(1:1)/= 'x') then
         do n=1,n_algae
            call accum_hist_field(n_bgc_N_ml(n,:), iblk, &
                   ocean_bio(:,:,nlt_bgc_N(n),      iblk), a2D)
         enddo
      endif
      if (f_bgc_Nit_ml(1:1)/= 'x') &
         call accum_hist_field(n_bgc_Nit_ml,  iblk, &
               ocean_bio(:,:,nlt_bgc_Nit,     iblk), a2D)  
      if (f_bgc_Am_ml(1:1)/= 'x') &
         call accum_hist_field(n_bgc_Am_ml,   iblk, &
               ocean_bio(:,:,nlt_bgc_Am,      iblk), a2D)  
      if (f_bgc_Sil_ml(1:1)/= 'x') &
         call accum_hist_field(n_bgc_Sil_ml,  iblk, &
               ocean_bio(:,:,nlt_bgc_Sil,     iblk), a2D)  
      if (f_bgc_hum_ml(1:1)/= 'x') &
         call accum_hist_field(n_bgc_hum_ml,  iblk, &
               ocean_bio(:,:,nlt_bgc_hum,     iblk), a2D)  
      if (f_bgc_DMSP_ml(1:1)/= 'x') &
         call accum_hist_field(n_bgc_DMSP_ml, iblk, &
               ocean_bio(:,:,nlt_bgc_DMSPd,   iblk), a2D)  
      if (f_bgc_DMS_ml(1:1)/= 'x') &
         call accum_hist_field(n_bgc_DMS_ml,  iblk, &
               ocean_bio(:,:,nlt_bgc_DMS,     iblk), a2D) 

      if (f_fNit  (1:1) /= 'x') &
         call accum_hist_field(n_fNit,     iblk, &
                  flux_bio(:,:,nlt_bgc_Nit,iblk), a2D)
      if (f_fNit_ai(1:1)/= 'x') &
         call accum_hist_field(n_fNit_ai,  iblk, &
               flux_bio_ai(:,:,nlt_bgc_Nit,iblk), a2D)

      if (f_fAm  (1:1) /= 'x') &
         call accum_hist_field(n_fAm,     iblk, &
                  flux_bio(:,:,nlt_bgc_Am,iblk), a2D)
      if (f_fAm_ai(1:1)/= 'x') &
         call accum_hist_field(n_fAm_ai,  iblk, &
               flux_bio_ai(:,:,nlt_bgc_Am,iblk), a2D)
      if (f_fN(1:1)/= 'x') then
         do n=1,n_algae
            call accum_hist_field(n_fN(n,:),    iblk, &
                  flux_bio(:,:,nlt_bgc_N(n),iblk), a2D)
         enddo
      endif
      if (f_fN_ai(1:1)/= 'x') then
         do n=1,n_algae
            call accum_hist_field(n_fN_ai(n,:),    iblk, &
                  flux_bio_ai(:,:,nlt_bgc_N(n),iblk), a2D)
         enddo
      endif
      if (f_fDOC(1:1)/= 'x') then
         do n=1,n_doc
            call accum_hist_field(n_fDOC(n,:),    iblk, &
                  flux_bio(:,:,nlt_bgc_DOC(n),iblk), a2D)
         enddo
      endif
      if (f_fDOC_ai(1:1)/= 'x') then
         do n=1,n_doc
            call accum_hist_field(n_fDOC_ai(n,:),    iblk, &
                  flux_bio_ai(:,:,nlt_bgc_DOC(n),iblk), a2D)
         enddo
      endif
      if (f_fDIC(1:1)/= 'x') then
         do n=1,n_dic
            call accum_hist_field(n_fDIC(n,:),    iblk, &
                  flux_bio(:,:,nlt_bgc_DIC(n),iblk), a2D)
         enddo
      endif
      if (f_fDIC_ai(1:1)/= 'x') then
         do n=1,n_dic
            call accum_hist_field(n_fDIC_ai(n,:),    iblk, &
                  flux_bio_ai(:,:,nlt_bgc_DIC(n),iblk), a2D)
         enddo
      endif
      if (f_fDON(1:1)/= 'x') then
         do n=1,n_don
            call accum_hist_field(n_fDON(n,:),    iblk, &
                  flux_bio(:,:,nlt_bgc_DON(n),iblk), a2D)
         enddo
      endif
      if (f_fDON_ai(1:1)/= 'x') then
         do n=1,n_don
            call accum_hist_field(n_fDON_ai(n,:),    iblk, &
                  flux_bio_ai(:,:,nlt_bgc_DON(n),iblk), a2D)
         enddo
      endif
      if (f_fFed (1:1)/= 'x') then
         do n=1,n_fed 
            call accum_hist_field(n_fFed (n,:),    iblk, &
                  flux_bio(:,:,nlt_bgc_Fed (n),iblk), a2D)
         enddo
      endif
      if (f_fFed_ai (1:1)/= 'x') then
         do n=1,n_fed 
            call accum_hist_field(n_fFed_ai (n,:),    iblk, &
                  flux_bio_ai(:,:,nlt_bgc_Fed (n),iblk), a2D)
         enddo
      endif
      if (f_fFep (1:1)/= 'x') then
         do n=1,n_fep 
            call accum_hist_field(n_fFep (n,:),    iblk, &
                  flux_bio(:,:,nlt_bgc_Fep (n),iblk), a2D)
         enddo
      endif
      if (f_fFep_ai (1:1)/= 'x') then
         do n=1,n_fep 
            call accum_hist_field(n_fFep_ai (n,:),    iblk, &
                  flux_bio_ai(:,:,nlt_bgc_Fep (n),iblk), a2D)
         enddo
      endif
      if (f_fSil  (1:1) /= 'x') &
         call accum_hist_field(n_fSil,    iblk, &
                 flux_bio(:,:,nlt_bgc_Sil,iblk), a2D)
      if (f_fSil_ai(1:1)/= 'x') &
         call accum_hist_field(n_fSil_ai, iblk, &
              flux_bio_ai(:,:,nlt_bgc_Sil,iblk), a2D)
      if (f_fhum  (1:1) /= 'x') &
         call accum_hist_field(n_fhum,    iblk, &
                 flux_bio(:,:,nlt_bgc_hum,iblk), a2D)
      if (f_fhum_ai(1:1)/= 'x') &
         call accum_hist_field(n_fhum_ai, iblk, &
              flux_bio_ai(:,:,nlt_bgc_hum,iblk), a2D)
      if (f_fPON  (1:1) /= 'x') &
         call accum_hist_field(n_fPON,    iblk, &
                 flux_bio(:,:,nlt_bgc_PON,iblk), a2D)
      if (f_fPON_ai(1:1)/= 'x') &
         call accum_hist_field(n_fPON_ai, iblk, &
              flux_bio_ai(:,:,nlt_bgc_PON,iblk), a2D)
      if (f_fDMS  (1:1) /= 'x') &
         call accum_hist_field(n_fDMS,    iblk, &
                 flux_bio(:,:,nlt_bgc_DMS,iblk), a2D)
      if (f_fDMS_ai(1:1)/= 'x') &
         call accum_hist_field(n_fDMS_ai, iblk, &
              flux_bio_ai(:,:,nlt_bgc_DMS,iblk), a2D)
      if (f_fDMSPd  (1:1) /= 'x') &
         call accum_hist_field(n_fDMSPd,    iblk, &
                 flux_bio(:,:,nlt_bgc_DMSPd,iblk), a2D)
      if (f_fDMSPd_ai(1:1)/= 'x') &
         call accum_hist_field(n_fDMSPd_ai, iblk, &
              flux_bio_ai(:,:,nlt_bgc_DMSPd,iblk), a2D)
      if (f_fDMSPp  (1:1) /= 'x') &
         call accum_hist_field(n_fDMSPp,    iblk, &
                 flux_bio(:,:,nlt_bgc_DMSPp,iblk), a2D)
      if (f_fDMSPp_ai(1:1)/= 'x') &
         call accum_hist_field(n_fDMSPp_ai, iblk, &
              flux_bio_ai(:,:,nlt_bgc_DMSPp,iblk), a2D)
      if (f_PPnet  (1:1) /= 'x') &
         call accum_hist_field(n_PPnet,   iblk, &
                               PP_net(:,:,iblk), a2D)
      if (f_grownet  (1:1) /= 'x') &
         call accum_hist_field(n_grownet, iblk, &
                             grow_net(:,:,iblk), a2D)   
      if (f_upNO (1:1) /= 'x') &
         call accum_hist_field(n_upNO,   iblk, &
                   upNO(:,:,iblk), a2D)
      if (f_upNH (1:1) /= 'x') &
         call accum_hist_field(n_upNH,   iblk, &
                   upNH(:,:,iblk), a2D)

      ! vertical biogeochemistry  

    if (z_tracers) then

      if (f_fzaero(1:1)/= 'x') then
         do n=1,n_zaero
            call accum_hist_field(n_fzaero(n,:),    iblk, &
                  flux_bio(:,:,nlt_zaero(n),iblk), a2D)
         enddo
      endif
      if (f_fzaero_ai(1:1)/= 'x') then
         do n=1,n_zaero
            call accum_hist_field(n_fzaero_ai(n,:),    iblk, &
                  flux_bio_ai(:,:,nlt_zaero(n),iblk), a2D)
         enddo
      endif
      if (f_algalpeak  (1:1) /= 'x') then
         do n=1,n_algae
            do j = jlo, jhi
               do i = ilo, ihi
                workii = 0
                maxv = c0
                if (aice(i,j,iblk) > c0) then
                  workv(:) = trcr(i,j,nt_bgc_N(n):nt_bgc_N(n)+nblyr,iblk)
                  worki = maxloc(workv)  ! integer for first location of max
                  maxv = workv(worki(1)) ! value of max
                  workv(worki(1)) = c0   ! remove maximum at first location
                  if (maxval(workv) - maxv < -puny) &
                    workii =  worki(1)
                endif
                workni(i,j) = real(workii,kind=dbl_kind)
                worknj(i,j) = maxv
               enddo  ! i
             enddo    ! j
             call accum_hist_field(n_algalpeak(n,:),    iblk, &
                   workni(:,:), a2D)
             call accum_hist_field(n_peakval(n,:), iblk, &
                   worknj(:,:), a2D)
         enddo      ! n
      endif !f_algalpeak

      !  
      ! ice_bio_net
      !
      if (f_zaeronet  (1:1) /= 'x') then
         do n=1,n_zaero
            call accum_hist_field(n_zaeronet(n,:),    iblk, &
                   ice_bio_net(:,:,nlt_zaero(n), iblk), a2D)
         enddo
      endif !f_zaeronet
      if (f_chlnet  (1:1) /= 'x') then
         do n=1,n_algae
            call accum_hist_field(n_chlnet(n,:),    iblk, &
                   R_chl2N(n)*ice_bio_net(:,:,nlt_bgc_N(n), iblk), a2D)
         enddo
      endif !f_chlnet
      if (f_Nnet  (1:1) /= 'x') then
         do n=1,n_algae
            call accum_hist_field(n_Nnet(n,:),    iblk, &
                   ice_bio_net(:,:,nlt_bgc_N(n), iblk), a2D)
         enddo
      endif !f_Nnet
      if (f_Cnet  (1:1) /= 'x') then
         do n=1,n_algae
            call accum_hist_field(n_Cnet(n,:),    iblk, &
                   R_C2N(n) *ice_bio_net(:,:,nlt_bgc_N(n), iblk), a2D)
         enddo
      endif !f_Cnet
      if (f_DOCnet  (1:1) /= 'x') then
         do n=1,n_doc  
            call accum_hist_field(n_DOCnet(n,:),    iblk, &
                   ice_bio_net(:,:,nlt_bgc_DOC(n), iblk), a2D)
         enddo
      endif !f_DOCnet
      if (f_DICnet  (1:1) /= 'x') then
         do n=1,n_dic  
            call accum_hist_field(n_DICnet(n,:),    iblk, &
                   ice_bio_net(:,:,nlt_bgc_DIC(n), iblk), a2D)
         enddo
      endif !f_DICnet
      if (f_DONnet  (1:1) /= 'x') then
         do n=1,n_don  
            call accum_hist_field(n_DONnet(n,:),    iblk, &
                   ice_bio_net(:,:,nlt_bgc_DON(n), iblk), a2D)
         enddo
      endif !f_DONnet
      if (f_Fednet   (1:1) /= 'x') then
         do n=1,n_fed   
            call accum_hist_field(n_Fednet (n,:),    iblk, &
                   ice_bio_net(:,:,nlt_bgc_Fed (n), iblk), a2D)
         enddo
      endif !f_Fednet 
      if (f_Fepnet   (1:1) /= 'x') then
         do n=1,n_fep   
            call accum_hist_field(n_Fepnet (n,:),    iblk, &
                   ice_bio_net(:,:,nlt_bgc_Fep (n), iblk), a2D)
         enddo
      endif !f_Fepnet 

      if (f_Nitnet  (1:1) /= 'x') &
         call accum_hist_field(n_Nitnet,   iblk, &
                   ice_bio_net(:,:,nlt_bgc_Nit, iblk), a2D)
      if (f_Amnet  (1:1) /= 'x') &
         call accum_hist_field(n_Amnet,   iblk, &
                   ice_bio_net(:,:,nlt_bgc_Am, iblk), a2D)
      if (f_Silnet  (1:1) /= 'x') &
         call accum_hist_field(n_Silnet,   iblk, &
                   ice_bio_net(:,:,nlt_bgc_Sil, iblk), a2D)
      if (f_humnet  (1:1) /= 'x') &
         call accum_hist_field(n_humnet,   iblk, &
                   ice_bio_net(:,:,nlt_bgc_hum, iblk), a2D)
      if (f_DMSPpnet  (1:1) /= 'x') &
         call accum_hist_field(n_DMSPpnet,   iblk, &
                  ice_bio_net(:,:,nlt_bgc_DMSPp, iblk), a2D)
      if (f_DMSPdnet  (1:1) /= 'x') &
         call accum_hist_field(n_DMSPdnet,   iblk, &
                   ice_bio_net(:,:,nlt_bgc_DMSPd, iblk), a2D)
      if (f_DMSnet  (1:1) /= 'x') &
         call accum_hist_field(n_DMSnet,   iblk, &
                   ice_bio_net(:,:,nlt_bgc_DMS, iblk), a2D)
      if (f_PONnet  (1:1) /= 'x') &
         call accum_hist_field(n_PONnet,   iblk, &
                   ice_bio_net(:,:,nlt_bgc_PON, iblk), a2D)
      !
      !  snow_bio_net
      ! 
      if (f_zaerosnow  (1:1) /= 'x') then
         do n=1,n_zaero
            call accum_hist_field(n_zaerosnow(n,:),    iblk, &
                   snow_bio_net(:,:,nlt_zaero(n), iblk), a2D)
         enddo
      endif !f_zaerosnow
      if (f_chlsnow  (1:1) /= 'x') then
         do n=1,n_algae
            call accum_hist_field(n_chlsnow(n,:),    iblk, &
                   R_chl2N(n)*snow_bio_net(:,:,nlt_bgc_N(n), iblk), a2D)
         enddo
      endif !f_chlsnow
      if (f_Nsnow  (1:1) /= 'x') then
         do n=1,n_algae
            call accum_hist_field(n_Nsnow(n,:),    iblk, &
                   snow_bio_net(:,:,nlt_bgc_N(n), iblk), a2D)
         enddo
      endif !f_Nsnow
      if (f_Csnow  (1:1) /= 'x') then
         do n=1,n_algae
            call accum_hist_field(n_Csnow(n,:),    iblk, &
                   R_C2N(n)*snow_bio_net(:,:,nlt_bgc_N(n), iblk), a2D)
         enddo
      endif !f_Csnow
      if (f_DOCsnow  (1:1) /= 'x') then
         do n=1,n_doc  
            call accum_hist_field(n_DOCsnow(n,:),    iblk, &
                   snow_bio_net(:,:,nlt_bgc_DOC(n), iblk), a2D)
         enddo
      endif !f_DOCsnow
      if (f_DICsnow  (1:1) /= 'x') then
         do n=1,n_dic  
            call accum_hist_field(n_DICsnow(n,:),    iblk, &
                   snow_bio_net(:,:,nlt_bgc_DIC(n), iblk), a2D)
         enddo
      endif !f_DICsnow
      if (f_DONsnow  (1:1) /= 'x') then
         do n=1,n_don  
            call accum_hist_field(n_DONsnow(n,:),    iblk, &
                   snow_bio_net(:,:,nlt_bgc_DON(n), iblk), a2D)
         enddo
      endif !f_DONsnow
      if (f_Fedsnow   (1:1) /= 'x') then
         do n=1,n_fed   
            call accum_hist_field(n_Fedsnow (n,:),    iblk, &
                   snow_bio_net(:,:,nlt_bgc_Fed (n), iblk), a2D)
         enddo
      endif !f_Fedsnow 
      if (f_Fepsnow   (1:1) /= 'x') then
         do n=1,n_fep   
            call accum_hist_field(n_Fepsnow (n,:),    iblk, &
                   snow_bio_net(:,:,nlt_bgc_Fep (n), iblk), a2D)
         enddo
      endif !f_Fepsnow 

      if (f_Nitsnow  (1:1) /= 'x') &
         call accum_hist_field(n_Nitsnow,   iblk, &
                   snow_bio_net(:,:,nlt_bgc_Nit, iblk), a2D)
      if (f_Amsnow  (1:1) /= 'x') &
         call accum_hist_field(n_Amsnow,   iblk, &
                   snow_bio_net(:,:,nlt_bgc_Am, iblk), a2D)
      if (f_Silsnow  (1:1) /= 'x') &
         call accum_hist_field(n_Silsnow,   iblk, &
                   snow_bio_net(:,:,nlt_bgc_Sil, iblk), a2D)
      if (f_humsnow  (1:1) /= 'x') &
         call accum_hist_field(n_humsnow,   iblk, &
                   snow_bio_net(:,:,nlt_bgc_hum, iblk), a2D)
      if (f_DMSPpsnow  (1:1) /= 'x') &
         call accum_hist_field(n_DMSPpsnow,   iblk, &
                   snow_bio_net(:,:,nlt_bgc_DMSPp, iblk), a2D)
      if (f_DMSPdsnow  (1:1) /= 'x') &
         call accum_hist_field(n_DMSPdsnow,   iblk, &
                   snow_bio_net(:,:,nlt_bgc_DMSPd, iblk), a2D)
      if (f_DMSsnow  (1:1) /= 'x') &
         call accum_hist_field(n_DMSsnow,   iblk, &
                   snow_bio_net(:,:,nlt_bgc_DMS, iblk), a2D)
      if (f_PONsnow  (1:1) /= 'x') &
         call accum_hist_field(n_PONsnow,   iblk, &
                   snow_bio_net(:,:,nlt_bgc_PON, iblk), a2D)
      !
      ! mobile frac
      ! 
      if (f_zaerofrac  (1:1) /= 'x') then
         do n=1,n_zaero
            call accum_hist_field(n_zaerofrac(n,:),    iblk, &
                   trcr(:,:,nt_zbgc_frac - 1 + nlt_zaero(n), iblk), a2D)
         enddo
      endif !f_zaerofrac
      if (f_chlfrac  (1:1) /= 'x') then
         do n=1,n_algae
            call accum_hist_field(n_chlfrac(n,:),    iblk, &
                   trcr(:,:,nlt_bgc_N(n), iblk), a2D)
         enddo
      endif !f_chlfrac
      if (f_Nfrac  (1:1) /= 'x') then
         do n=1,n_algae
            call accum_hist_field(n_Nfrac(n,:),    iblk, &
                   trcr(:,:,nt_zbgc_frac - 1 + nlt_bgc_N(n), iblk), a2D)
         enddo
      endif !f_Nfrac
      if (f_DOCfrac  (1:1) /= 'x') then
         do n=1,n_doc  
            call accum_hist_field(n_DOCfrac(n,:),    iblk, &
                   trcr(:,:,nt_zbgc_frac - 1 + nlt_bgc_DOC(n), iblk), a2D)
         enddo
      endif !f_DOCfrac
      if (f_DICfrac  (1:1) /= 'x') then
         do n=1,n_dic  
            call accum_hist_field(n_DICfrac(n,:),    iblk, &
                   trcr(:,:,nt_zbgc_frac - 1 + nlt_bgc_DIC(n), iblk), a2D)
         enddo
      endif !f_DICfrac
      if (f_DONfrac  (1:1) /= 'x') then
         do n=1,n_don  
            call accum_hist_field(n_DONfrac(n,:),    iblk, &
                   trcr(:,:,nt_zbgc_frac - 1 + nlt_bgc_DON(n), iblk), a2D)
         enddo
      endif !f_DONfrac
      if (f_Fedfrac   (1:1) /= 'x') then
         do n=1,n_fed   
            call accum_hist_field(n_Fedfrac (n,:),    iblk, &
                   trcr(:,:,nt_zbgc_frac - 1 + nlt_bgc_Fed (n), iblk), a2D)
         enddo
      endif !f_Fedfrac 
      if (f_Fepfrac   (1:1) /= 'x') then
         do n=1,n_fep   
            call accum_hist_field(n_Fepfrac (n,:),    iblk, &
                   trcr(:,:,nt_zbgc_frac - 1 + nlt_bgc_Fep (n), iblk), a2D)
         enddo
      endif !f_Fepfrac 

      if (f_Nitfrac  (1:1) /= 'x') &
         call accum_hist_field(n_Nitfrac,   iblk, &
                   trcr(:,:,nt_zbgc_frac - 1 + nlt_bgc_Nit, iblk), a2D)
      if (f_Amfrac  (1:1) /= 'x') &
         call accum_hist_field(n_Amfrac,   iblk, &
                   trcr(:,:,nt_zbgc_frac - 1 + nlt_bgc_Am, iblk), a2D)
      if (f_Silfrac  (1:1) /= 'x') &
         call accum_hist_field(n_Silfrac,   iblk, &
                   trcr(:,:,nt_zbgc_frac - 1 + nlt_bgc_Sil, iblk), a2D)
      if (f_humfrac  (1:1) /= 'x') &
         call accum_hist_field(n_humfrac,   iblk, &
                   trcr(:,:,nt_zbgc_frac - 1 + nlt_bgc_hum, iblk), a2D) 
      if (f_DMSPpfrac  (1:1) /= 'x') &
         call accum_hist_field(n_DMSPpfrac,   iblk, &
                   trcr(:,:,nt_zbgc_frac - 1 + nlt_bgc_DMSPp, iblk), a2D)
      if (f_DMSPdfrac  (1:1) /= 'x') &
         call accum_hist_field(n_DMSPdfrac,   iblk, &
                   trcr(:,:,nt_zbgc_frac - 1 + nlt_bgc_DMSPd, iblk), a2D)
      if (f_DMSfrac  (1:1) /= 'x') &
         call accum_hist_field(n_DMSfrac,   iblk, &
                   trcr(:,:,nt_zbgc_frac - 1 + nlt_bgc_DMS, iblk), a2D)
      if (f_PONfrac  (1:1) /= 'x') &
         call accum_hist_field(n_PONfrac,   iblk, &
                   trcr(:,:,nt_zbgc_frac - 1 + nlt_bgc_PON, iblk), a2D)

    endif  ! z_tracers

      ! brine
      if (f_hbri  (1:1) /= 'x') &
         call accum_hist_field(n_hbri,     iblk, &
                        hbri(:,:,iblk), a2D)

    endif  ! 2d bgc tracers, tr_aero, tr_brine, solve_zsal, skl_bgc


      ! 3D category fields

    if (tr_brine) then
      ! 3Dc bgc category fields

      if (f_fbri   (1:1) /= 'x') &
         call accum_hist_field(n_fbri-n2D, iblk, ncat_hist, &
                               trcrn(:,:,nt_fbri,1:ncat_hist,iblk), a3Dc)
    endif

    if (z_tracers .or. solve_zsal) then
      ! 3Db category fields

      if (f_bTin  (1:1) /= 'x')  then
         workz(:,:,:) = c0
         do n = 1, ncat_hist
            do j = jlo, jhi
               do i = ilo, ihi
                  if (aicen(i,j,n,iblk) > c0) then
                      workz(i,j,1:nzblyr) =  workz(i,j,1:nzblyr) + bTiz(i,j,1:nzblyr,n,iblk)*&
                               aicen(i,j,n,iblk)
                  endif
               enddo  ! i
            enddo     ! j
         enddo        ! n
         call accum_hist_field(n_bTin-n3Dzcum, iblk, nzblyr, &
                               workz(:,:,1:nzblyr), a3Db)
      endif

      if (f_bphi  (1:1) /= 'x') then
         workz(:,:,:) = c0
         do n = 1, ncat_hist
            do j = jlo, jhi
               do i = ilo, ihi
                  if (aicen(i,j,n,iblk) > c0) then
                      workz(i,j,1:nzblyr) =  workz(i,j,1:nzblyr) + bphi(i,j,1:nzblyr,n,iblk)*&
                               aicen(i,j,n,iblk)
                  endif
               enddo ! i
            enddo    ! j
         enddo       ! n
         call accum_hist_field(n_bphi-n3Dzcum, iblk, nzblyr, &
                               workz(:,:,1:nzblyr), a3Db)
      endif

      if (f_bgc_S   (1:1) /= 'x') then
         workz(:,:,:) = c0
            do j = jlo, jhi
               do i = ilo, ihi
                  if (aice(i,j,iblk) > c0) then
                    workz(i,j,1) = trcr(i,j,nt_bgc_S,iblk)
                    workz(i,j,2:nblyr+1) = trcr(i,j,nt_bgc_S:nt_bgc_S+nblyr-1,iblk)
                    workz(i,j,nblyr+2) = sss(i,j,iblk)
                  endif
                enddo ! i
            enddo     ! j
         call accum_hist_field(n_bgc_S-n3Dzcum, iblk, nzblyr, &
                                  workz(:,:,1:nzblyr), a3Db)
      endif

      if (f_zfswin   (1:1) /= 'x') then
         workz(:,:,:) = c0
         do n = 1, ncat_hist
            do j = jlo, jhi
               do i = ilo, ihi
                if (aicen(i,j,n,iblk) > c0) then
                   workz(i,j,1:nblyr+1) = workz(i,j,1:nblyr+1)+ zfswin(i,j,1:nblyr+1,n,iblk)*&
                       aicen(i,j,n,iblk)
                   workz(i,j,nzblyr) = workz(i,j,nzblyr) + zfswin(i,j,nblyr+1,n,iblk)*&
                       aicen(i,j,n,iblk)
                endif
               enddo ! i
            enddo    ! j
         enddo       ! n
         call accum_hist_field(n_zfswin-n3Dzcum, iblk, nzblyr, &
                                  workz(:,:,1:nzblyr), a3Db)
      endif

      if (f_iDi   (1:1) /= 'x') then
         workz(:,:,:) = c0
         do n = 1, ncat_hist
            do k = 1,nzblyr-1
              do j = jlo, jhi
               do i = ilo, ihi
                if (aicen(i,j,n,iblk) > c0) then
                   workz(i,j,k) = workz(i,j,k) + iDi(i,j,k,n,iblk)*vicen(i,j,n,iblk)**2/aicen(i,j,n,iblk)
                   workz(i,j,nzblyr)   = workz(i,j,nzblyr-1)     
                endif
               enddo ! i
              enddo  ! j
            enddo    ! k        
          enddo      ! n
          call accum_hist_field(n_iDi-n3Dzcum, iblk, nzblyr, &
                                  workz(:,:,1:nzblyr), a3Db)
      endif

      if (f_iki   (1:1) /= 'x') then
         workz(:,:,:) = c0
         do n = 1, ncat_hist
            do k = 1,nzblyr-1
            do j = jlo, jhi
               do i = ilo, ihi
                if (aicen(i,j,n,iblk) > c0) then
                   workz(i,j,k) = workz(i,j,k) + iki(i,j,k,n,iblk)*&
                       aicen(i,j,n,iblk)
                   workz(i,j,nzblyr)   = workz(i,j,nzblyr) + iki(i,j,nblyr+1,n,iblk)*&
                       aicen(i,j,n,iblk)
                endif
               enddo ! i
            enddo    ! j
            enddo    ! k  
         enddo       ! n
         call accum_hist_field(n_iki-n3Dzcum, iblk, nzblyr,  &
                                  workz(:,:,1:nzblyr), a3Db)
      endif

    endif  ! 3Db fields

    if (z_tracers) then
      ! 3Da category fields

      if (f_zaero   (1:1) /= 'x') then
         do k = 1,n_zaero
           workz(:,:,:) = c0
             do j = jlo, jhi
               do i = ilo, ihi
                  if (aice(i,j,iblk) > puny) then  
                    workz(i,j,1:2) =  &   !snow
                         trcr(i,j,nt_zaero(k)+nblyr+1:nt_zaero(k)+nblyr+2,iblk)/rhos
                    workz(i,j,3:nblyr+3) = &  !ice
                         trcr(i,j,nt_zaero(k):nt_zaero(k)+nblyr,iblk)/rhoi
                    workz(i,j,nblyr+4) =  ocean_bio(i,j,nlt_zaero(k),iblk)/rhow !ocean
                  endif
                enddo ! i
             enddo    ! j 
           call accum_hist_field(n_zaeros(k,:)-n3Dbcum, iblk, nzalyr, &
                                  workz(:,:,1:nzalyr), a3Da)
         enddo !k
      endif  !f_zaero
      if (f_bgc_N   (1:1) /= 'x') then
         do k = 1,n_algae
           workz(:,:,:) = c0
           workz2(:,:,:) = c0
             do j = jlo, jhi
               do i = ilo, ihi
                  if (aice(i,j,iblk) > puny) then  
                    workz(i,j,1:2) =  &   !snow
                         trcr(i,j,nt_bgc_N(k)+nblyr+1:nt_bgc_N(k)+nblyr+2,iblk)
                    workz(i,j,3:nblyr+3) = &  !ice
                         trcr(i,j,nt_bgc_N(k):nt_bgc_N(k)+nblyr,iblk)
                    workz(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_N(k),iblk) !ocean
                  endif
                  if (aicen(i,j,1,iblk) > puny) then  
                    workz2(i,j,1:2) =  &   !snow
                         trcrn(i,j,nt_bgc_N(k)+nblyr+1:nt_bgc_N(k)+nblyr+2,1,iblk)
                    workz2(i,j,3:nblyr+3) = &  !ice
                         trcrn(i,j,nt_bgc_N(k):nt_bgc_N(k)+nblyr,1,iblk)
                    workz2(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_N(k),iblk) !ocean
                  endif
                enddo ! i
             enddo    ! j 
           call accum_hist_field(n_bgc_N(k,:)-n3Dbcum, iblk, nzalyr,  &
                                  workz(:,:,1:nzalyr), a3Da)
           call accum_hist_field(n_bgc_N_cat1(k,:)-n3Dbcum, iblk, nzalyr,  &
                                  workz2(:,:,1:nzalyr), a3Da)
         enddo !k
      endif  !f_bgc_N
      if (f_bgc_C   (1:1) /= 'x') then
         do k = 1,n_algae
           workz(:,:,:) = c0
             do j = jlo, jhi
               do i = ilo, ihi
                  if (aice(i,j,iblk) > puny) then  
                    workz(i,j,1:2) =  &   !snow
                         R_C2N(k)*trcr(i,j,nt_bgc_N(k)+nblyr+1:nt_bgc_N(k)+nblyr+2,iblk)
                    workz(i,j,3:nblyr+3) = &  !ice
                         R_C2N(k)*trcr(i,j,nt_bgc_N(k):nt_bgc_N(k)+nblyr,iblk)
                    workz(i,j,nblyr+4) =  R_C2N(k)*ocean_bio(i,j,nlt_bgc_N(k),iblk) !ocean
                  endif
                enddo ! i
             enddo    ! j 
           call accum_hist_field(n_bgc_C(k,:)-n3Dbcum, iblk, nzalyr, &
                                  workz(:,:,1:nzalyr), a3Da)
         enddo !k
      endif  !f_bgc_C
      if (f_bgc_DOC   (1:1) /= 'x') then
         do k = 1,n_doc   
           workz(:,:,:) = c0
           workz2(:,:,:) = c0
             do j = jlo, jhi
               do i = ilo, ihi
                  if (aice(i,j,iblk) > puny) then  
                    workz(i,j,1:2) =  &   !snow
                         trcr(i,j,nt_bgc_DOC(k)+nblyr+1:nt_bgc_DOC(k)+nblyr+2,iblk)
                    workz(i,j,3:nblyr+3) = &  !ice
                         trcr(i,j,nt_bgc_DOC(k):nt_bgc_DOC(k)+nblyr,iblk)
                    workz(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_DOC(k),iblk) !ocean
                  endif
                  if (aicen(i,j,1,iblk) > puny) then  
                    workz2(i,j,1:2) =  &   !snow
                         trcrn(i,j,nt_bgc_DOC(k)+nblyr+1:nt_bgc_DOC(k)+nblyr+2,1,iblk)
                    workz2(i,j,3:nblyr+3) = &  !ice
                         trcrn(i,j,nt_bgc_DOC(k):nt_bgc_DOC(k)+nblyr,1,iblk)
                    workz2(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_DOC(k),iblk) !ocean
                  endif
                enddo ! i
             enddo    ! j 
           call accum_hist_field(n_bgc_DOC(k,:)-n3Dbcum, iblk, nzalyr,  &
                                  workz(:,:,1:nzalyr), a3Da)
           call accum_hist_field(n_bgc_DOC_cat1(k,:)-n3Dbcum, iblk, nzalyr,  &
                                  workz2(:,:,1:nzalyr), a3Da)
         enddo !k
      endif  !f_bgc_DOC
      if (f_bgc_DIC   (1:1) /= 'x') then
         do k = 1,n_dic  
           workz(:,:,:) = c0
           workz2(:,:,:) = c0
             do j = jlo, jhi
               do i = ilo, ihi
                  if (aice(i,j,iblk) > puny) then  
                    workz(i,j,1:2) =  &   !snow
                         trcr(i,j,nt_bgc_DIC(k)+nblyr+1:nt_bgc_DIC(k)+nblyr+2,iblk)
                    workz(i,j,3:nblyr+3) = &  !ice
                         trcr(i,j,nt_bgc_DIC(k):nt_bgc_DIC(k)+nblyr,iblk)
                    workz(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_DIC(k),iblk) !ocean
                  endif
                  if (aicen(i,j,1,iblk) > puny) then  
                    workz2(i,j,1:2) =  &   !snow
                         trcrn(i,j,nt_bgc_DIC(k)+nblyr+1:nt_bgc_DIC(k)+nblyr+2,1,iblk)
                    workz2(i,j,3:nblyr+3) = &  !ice
                         trcrn(i,j,nt_bgc_DIC(k):nt_bgc_DIC(k)+nblyr,1,iblk)
                    workz2(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_DIC(k),iblk) !ocean
                  endif
                enddo ! i
             enddo    ! j 
           call accum_hist_field(n_bgc_DIC(k,:)-n3Dbcum, iblk, nzalyr,  &
                                  workz(:,:,1:nzalyr), a3Da)
          call accum_hist_field(n_bgc_DIC_cat1(k,:)-n3Dbcum, iblk, nzalyr,  &
                                  workz2(:,:,1:nzalyr), a3Da)
         enddo !k
      endif  !f_bgc_DIC
      if (f_bgc_DON   (1:1) /= 'x') then
         do k = 1,n_don  
           workz(:,:,:) = c0
           workz2(:,:,:) = c0
             do j = jlo, jhi
               do i = ilo, ihi
                  if (aice(i,j,iblk) > puny) then  
                    workz(i,j,1:2) =  &   !snow
                         trcr(i,j,nt_bgc_DON(k)+nblyr+1:nt_bgc_DON(k)+nblyr+2,iblk)
                    workz(i,j,3:nblyr+3) = &  !ice
                         trcr(i,j,nt_bgc_DON(k):nt_bgc_DON(k)+nblyr,iblk)
                    workz(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_DON(k),iblk) !ocean
                  endif
                  if (aicen(i,j,1,iblk) > puny) then  
                    workz2(i,j,1:2) =  &   !snow
                         trcrn(i,j,nt_bgc_DON(k)+nblyr+1:nt_bgc_DON(k)+nblyr+2,1,iblk)
                    workz2(i,j,3:nblyr+3) = &  !ice
                         trcrn(i,j,nt_bgc_DON(k):nt_bgc_DON(k)+nblyr,1,iblk)
                    workz2(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_DON(k),iblk) !ocean
                  endif
                enddo ! i
             enddo    ! j 
           call accum_hist_field(n_bgc_DON(k,:)-n3Dbcum, iblk, nzalyr,  &
                                  workz(:,:,1:nzalyr), a3Da)
           call accum_hist_field(n_bgc_DON_cat1(k,:)-n3Dbcum, iblk, nzalyr,  &
                                  workz2(:,:,1:nzalyr), a3Da)
         enddo !k
      endif  !f_bgc_DON
      if (f_bgc_Fed    (1:1) /= 'x') then
         do k = 1,n_fed   
           workz(:,:,:) = c0
           workz2(:,:,:) = c0
             do j = jlo, jhi
               do i = ilo, ihi
                  if (aice(i,j,iblk) > puny) then  
                    workz(i,j,1:2) =  &   !snow
                         trcr(i,j,nt_bgc_Fed (k)+nblyr+1:nt_bgc_Fed (k)+nblyr+2,iblk)
                    workz(i,j,3:nblyr+3) = &  !ice
                         trcr(i,j,nt_bgc_Fed (k):nt_bgc_Fed (k)+nblyr,iblk)
                    workz(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_Fed (k),iblk) !ocean
                  endif
                  if (aicen(i,j,1,iblk) > puny) then  
                    workz2(i,j,1:2) =  &   !snow
                         trcrn(i,j,nt_bgc_Fed (k)+nblyr+1:nt_bgc_Fed (k)+nblyr+2,1,iblk)
                    workz2(i,j,3:nblyr+3) = &  !ice
                         trcrn(i,j,nt_bgc_Fed (k):nt_bgc_Fed (k)+nblyr,1,iblk)
                    workz2(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_Fed (k),iblk) !ocean
                  endif
                enddo ! i
             enddo    ! j 
           call accum_hist_field(n_bgc_Fed (k,:)-n3Dbcum, iblk, nzalyr,  &
                                  workz(:,:,1:nzalyr), a3Da)
           call accum_hist_field(n_bgc_Fed_cat1 (k,:)-n3Dbcum, iblk, nzalyr,  &
                                  workz2(:,:,1:nzalyr), a3Da)
         enddo !k
      endif  !f_bgc_Fed 
      if (f_bgc_Fep    (1:1) /= 'x') then
         do k = 1,n_fep   
           workz(:,:,:) = c0
           workz2(:,:,:) = c0
             do j = jlo, jhi
               do i = ilo, ihi
                  if (aice(i,j,iblk) > puny) then  
                    workz(i,j,1:2) =  &   !snow
                         trcr(i,j,nt_bgc_Fep (k)+nblyr+1:nt_bgc_Fep (k)+nblyr+2,iblk)
                    workz(i,j,3:nblyr+3) = &  !ice
                         trcr(i,j,nt_bgc_Fep (k):nt_bgc_Fep (k)+nblyr,iblk)
                    workz(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_Fep (k),iblk) !ocean
                  endif
                  if (aicen(i,j,1,iblk) > puny) then  
                    workz2(i,j,1:2) =  &   !snow
                         trcrn(i,j,nt_bgc_Fep (k)+nblyr+1:nt_bgc_Fep (k)+nblyr+2,1,iblk)
                    workz2(i,j,3:nblyr+3) = &  !ice
                         trcrn(i,j,nt_bgc_Fep (k):nt_bgc_Fep (k)+nblyr,1,iblk)
                    workz2(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_Fep (k),iblk) !ocean
                  endif
                enddo ! i
             enddo    ! j 
           call accum_hist_field(n_bgc_Fep (k,:)-n3Dbcum, iblk, nzalyr,  &
                                  workz(:,:,1:nzalyr), a3Da)
           call accum_hist_field(n_bgc_Fep_cat1 (k,:)-n3Dbcum, iblk, nzalyr,  &
                                  workz2(:,:,1:nzalyr), a3Da)
         enddo !k
      endif  !f_bgc_Fep 
      if (f_bgc_chl   (1:1) /= 'x') then
         do k = 1,n_algae
           workz(:,:,:) = c0
             do j = jlo, jhi
               do i = ilo, ihi
                  if (aice(i,j,iblk) > puny) then  
                    workz(i,j,1:2) =  &   !snow
                         trcr(i,j,nt_bgc_chl(k)+nblyr+1:nt_bgc_chl(k)+nblyr+2,iblk)
                    workz(i,j,3:nblyr+3) = &  !ice
                         trcr(i,j,nt_bgc_chl(k):nt_bgc_chl(k)+nblyr,iblk)
                    workz(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_chl(k),iblk) !ocean
                  endif
                enddo ! i
             enddo    ! j 
           call accum_hist_field(n_bgc_chl(k,:)-n3Dbcum, iblk, nzalyr,  &
                                  workz(:,:,1:nzalyr), a3Da)
         enddo !k
      endif  !f_bgc_chl
         
      if (f_bgc_Nit   (1:1) /= 'x') then
         workz(:,:,:) = c0
         workz2(:,:,:) = c0
            do j = jlo, jhi
               do i = ilo, ihi
                  if (aice(i,j,iblk) > puny) then
                    workz(i,j,1:2) =  &   !snow
                         trcr(i,j,nt_bgc_Nit+nblyr+2:nt_bgc_Nit+nblyr+3,iblk)
                    workz(i,j,3:nblyr+3) = &  !ice
                         trcr(i,j,nt_bgc_Nit:nt_bgc_Nit+nblyr,iblk)  
                    workz(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_Nit,iblk) !ocean
                  endif
                  if (aicen(i,j,1,iblk) > puny) then
                    workz2(i,j,1:2) =  &   !snow
                         trcrn(i,j,nt_bgc_Nit+nblyr+2:nt_bgc_Nit+nblyr+3,1,iblk)
                    workz2(i,j,3:nblyr+3) = &  !ice
                         trcrn(i,j,nt_bgc_Nit:nt_bgc_Nit+nblyr,1,iblk)  
                    workz2(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_Nit,iblk) !ocean
                  endif
                enddo ! i
             enddo    ! j 
         call accum_hist_field(n_bgc_Nit-n3Dbcum, iblk, nzalyr,  &
                                  workz(:,:,1:nzalyr), a3Da)
         call accum_hist_field(n_bgc_Nit_cat1-n3Dbcum, iblk, nzalyr,  &
                                  workz2(:,:,1:nzalyr), a3Da)
      endif

      if (f_bgc_Am   (1:1) /= 'x') then
         workz(:,:,:) = c0
         workz2(:,:,:) = c0
            do j = jlo, jhi
               do i = ilo, ihi
                  if (aice(i,j,iblk) > puny) then
                    workz(i,j,1:2) =  &   !snow
                         trcr(i,j,nt_bgc_Am+nblyr+1:nt_bgc_Am+nblyr+2,iblk)
                    workz(i,j,3:nblyr+3) = &  !ice
                         trcr(i,j,nt_bgc_Am:nt_bgc_Am+nblyr,iblk)  
                    workz(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_Am,iblk) !ocean
                  endif
                  if (aicen(i,j,1,iblk) > puny) then
                    workz2(i,j,1:2) =  &   !snow
                         trcrn(i,j,nt_bgc_Am+nblyr+1:nt_bgc_Am+nblyr+2,1,iblk)
                    workz2(i,j,3:nblyr+3) = &  !ice
                         trcrn(i,j,nt_bgc_Am:nt_bgc_Am+nblyr,1,iblk)  
                    workz2(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_Am,iblk) !ocean
                  endif
                enddo ! i
             enddo    ! j 
         call accum_hist_field(n_bgc_Am-n3Dbcum, iblk, nzalyr,  &
                                  workz(:,:,1:nzalyr), a3Da)
         call accum_hist_field(n_bgc_Am_cat1-n3Dbcum, iblk, nzalyr,  &
                                  workz2(:,:,1:nzalyr), a3Da)
      endif

      if (f_bgc_Sil   (1:1) /= 'x') then
         workz(:,:,:) = c0
         workz2(:,:,:) = c0
            do j = jlo, jhi
               do i = ilo, ihi
                  if (aice(i,j,iblk) > puny) then
                    workz(i,j,1:2) =  &   !snow
                         trcr(i,j,nt_bgc_Sil+nblyr+1:nt_bgc_Sil+nblyr+2,iblk)
                    workz(i,j,3:nblyr+3) = &  !ice
                         trcr(i,j,nt_bgc_Sil:nt_bgc_Sil+nblyr,iblk) 
                    workz(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_Sil,iblk) !ocean
                  endif
                  if (aicen(i,j,1,iblk) > puny) then
                    workz2(i,j,1:2) =  &   !snow
                         trcrn(i,j,nt_bgc_Sil+nblyr+1:nt_bgc_Sil+nblyr+2,1,iblk)
                    workz2(i,j,3:nblyr+3) = &  !ice
                         trcrn(i,j,nt_bgc_Sil:nt_bgc_Sil+nblyr,1,iblk) 
                    workz2(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_Sil,iblk) !ocean
                  endif
                enddo ! i
             enddo    ! j 
         call accum_hist_field(n_bgc_Sil-n3Dbcum, iblk, nzalyr,  &
                                  workz(:,:,1:nzalyr), a3Da)
         call accum_hist_field(n_bgc_Sil_cat1-n3Dbcum, iblk, nzalyr,  &
                                  workz2(:,:,1:nzalyr), a3Da)
      endif
  

      if (f_bgc_hum   (1:1) /= 'x') then
         workz(:,:,:) = c0
         workz2(:,:,:) = c0
            do j = jlo, jhi
               do i = ilo, ihi
                  if (aice(i,j,iblk) > puny) then
                    workz(i,j,1:2) =  &   !snow
                         trcr(i,j,nt_bgc_hum+nblyr+1:nt_bgc_hum+nblyr+2,iblk)
                    workz(i,j,3:nblyr+3) = &  !ice
                         trcr(i,j,nt_bgc_hum:nt_bgc_hum+nblyr,iblk) 
                    workz(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_hum,iblk) !ocean
                  endif
                  if (aicen(i,j,1,iblk) > puny) then
                    workz2(i,j,1:2) =  &   !snow
                         trcrn(i,j,nt_bgc_hum+nblyr+1:nt_bgc_hum+nblyr+2,1,iblk)
                    workz2(i,j,3:nblyr+3) = &  !ice
                         trcrn(i,j,nt_bgc_hum:nt_bgc_hum+nblyr,1,iblk) 
                    workz2(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_hum,iblk) !ocean
                  endif
                enddo ! i
             enddo    ! j 
         call accum_hist_field(n_bgc_hum-n3Dbcum, iblk, nzalyr,  &
                                  workz(:,:,1:nzalyr), a3Da)
         call accum_hist_field(n_bgc_hum_cat1-n3Dbcum, iblk, nzalyr,  &
                                  workz2(:,:,1:nzalyr), a3Da)
      endif
 
      if (f_bgc_DMSPd   (1:1) /= 'x') then
         workz(:,:,:) = c0
         workz2(:,:,:) = c0
            do j = jlo, jhi
               do i = ilo, ihi
                  if (aice(i,j,iblk) > puny) then
                    workz(i,j,1:2) =  &   !snow
                         trcr(i,j,nt_bgc_DMSPd+nblyr+1:nt_bgc_DMSPd+nblyr+2,iblk)
                    workz(i,j,3:nblyr+3) = &  !ice
                         trcr(i,j,nt_bgc_DMSPd:nt_bgc_DMSPd+nblyr,iblk) 
                    workz(i,j,nblyr+4) = ocean_bio(i,j,nlt_bgc_DMSPd,iblk) !ocean
                  endif
                  if (aicen(i,j,1,iblk) > puny) then
                    workz2(i,j,1:2) =  &   !snow
                         trcrn(i,j,nt_bgc_DMSPd+nblyr+1:nt_bgc_DMSPd+nblyr+2,1,iblk)
                    workz2(i,j,3:nblyr+3) = &  !ice
                         trcrn(i,j,nt_bgc_DMSPd:nt_bgc_DMSPd+nblyr,1,iblk) 
                    workz2(i,j,nblyr+4) = ocean_bio(i,j,nlt_bgc_DMSPd,iblk) !ocean
                  endif
                enddo ! i
             enddo    ! j 
         call accum_hist_field(n_bgc_DMSPd-n3Dbcum, iblk, nzalyr,  &
                                  workz(:,:,1:nzalyr), a3Da)
         call accum_hist_field(n_bgc_DMSPd_cat1-n3Dbcum, iblk, nzalyr,  &
                                  workz2(:,:,1:nzalyr), a3Da)
      endif       

      if (f_bgc_DMSPp   (1:1) /= 'x') then
         workz(:,:,:) = c0
            do j = jlo, jhi
               do i = ilo, ihi
                  if (aice(i,j,iblk) > puny) then
                    workz(i,j,1:2) =  &   !snow
                         trcr(i,j,nt_bgc_DMSPp+nblyr+1:nt_bgc_DMSPp+nblyr+2,iblk)
                    workz(i,j,3:nblyr+3) = &  !ice
                         trcr(i,j,nt_bgc_DMSPp:nt_bgc_DMSPp+nblyr,iblk) 
                    workz(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_DMSPp,iblk) !ocean
                  endif
                enddo ! i
             enddo    ! j 
         call accum_hist_field(n_bgc_DMSPp-n3Dbcum, iblk, nzalyr,  &
                                  workz(:,:,1:nzalyr), a3Da)
      endif

      if (f_bgc_DMS   (1:1) /= 'x') then
         workz(:,:,:) = c0
         workz2(:,:,:) = c0
            do j = jlo, jhi
               do i = ilo, ihi
                  if (aice(i,j,iblk) > puny) then
                    workz(i,j,1:2) =  &   !snow
                         trcr(i,j,nt_bgc_DMS+nblyr+1:nt_bgc_DMS+nblyr+2,iblk)
                    workz(i,j,3:nblyr+3) = &  !ice
                         trcr(i,j,nt_bgc_DMS:nt_bgc_DMS+nblyr,iblk) 
                    workz(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_DMS,iblk) !ocean
                  endif
                  if (aicen(i,j,1,iblk) > puny) then
                    workz2(i,j,1:2) =  &   !snow
                         trcrn(i,j,nt_bgc_DMS+nblyr+1:nt_bgc_DMS+nblyr+2,1,iblk)
                    workz2(i,j,3:nblyr+3) = &  !ice
                         trcrn(i,j,nt_bgc_DMS:nt_bgc_DMS+nblyr,1,iblk) 
                    workz2(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_DMS,iblk) !ocean
                  endif
                enddo ! i
             enddo    ! j 
         call accum_hist_field(n_bgc_DMS-n3Dbcum, iblk, nzalyr,  &
                                  workz(:,:,1:nzalyr), a3Da)
         call accum_hist_field(n_bgc_DMS_cat1-n3Dbcum, iblk, nzalyr,  &
                                  workz2(:,:,1:nzalyr), a3Da)
      endif

      if (f_bgc_PON   (1:1) /= 'x') then
         workz(:,:,:) = c0
         workz2(:,:,:) = c0
            do j = jlo, jhi
               do i = ilo, ihi
                  if (aice(i,j,iblk) > puny) then
                    workz(i,j,1:2) =  &   !snow
                         trcr(i,j,nt_bgc_PON+nblyr+1:nt_bgc_PON+nblyr+2,iblk)
                    workz(i,j,3:nblyr+3) = &  !ice
                         trcr(i,j,nt_bgc_PON:nt_bgc_PON+nblyr,iblk) 
                    workz(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_PON,iblk) !ocean
                  endif
                  if (aicen(i,j,1,iblk) > puny) then
                    workz2(i,j,1:2) =  &   !snow
                         trcrn(i,j,nt_bgc_PON+nblyr+1:nt_bgc_PON+nblyr+2,1,iblk)
                    workz2(i,j,3:nblyr+3) = &  !ice
                         trcrn(i,j,nt_bgc_PON:nt_bgc_PON+nblyr,1,iblk) 
                    workz2(i,j,nblyr+4) =  ocean_bio(i,j,nlt_bgc_PON,iblk) !ocean
                  endif
                enddo ! i
             enddo    ! j 
         call accum_hist_field(n_bgc_PON-n3Dbcum, iblk, nzalyr,  &
                                  workz(:,:,1:nzalyr), a3Da)
         call accum_hist_field(n_bgc_PON_cat1-n3Dbcum, iblk, nzalyr,  &
                                  workz2(:,:,1:nzalyr), a3Da)
      endif

    endif   ! z_tracers, 3Da tracers

      end subroutine accum_hist_bgc

!=======================================================================

      subroutine init_hist_bgc_3Da

      use ice_calendar, only: nstreams
      use ice_history_shared, only: tstr3Da, tcstr, define_hist_field

      integer (kind=int_kind) :: ns, n
      logical (kind=log_kind) :: z_tracers
      character (len=3) :: nchar
      character (len=16):: vname_in     ! variable name
      character(len=*), parameter :: subname = '(init_hist_bgc_3Da)'
      
      call icepack_query_parameters(z_tracers_out=z_tracers)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

    ! snow+bio grid

    if (z_tracers) then

    do ns = 1, nstreams
 
!----------------------------------------------------------------------------
! snow+bio grid ==>
! 1:2 snow (surface layer +interior), 3:nblyr+2 ice (bio grid), nblyr+3 ocean
!----------------------------------------------------------------------------

       if (f_zaero(1:1) /= 'x') then
         do n=1,n_zaero
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'zaero', trim(nchar)
            call define_hist_field(n_zaeros(n,:),vname_in,"kg/kg",tstr3Da, tcstr, &
                "bulk z aerosol mass fraction", "snow+bio grid", c1, c0, &
                ns, f_zaero)
         enddo
       endif
     
       if (f_bgc_Nit(1:1) /= 'x') & 
            call define_hist_field(n_bgc_Nit,"bgc_Nit","mmol/m^3",tstr3Da, tcstr, &
                "bulk nitrate ", "snow+bio grid", c1, c0,     &
                ns, f_bgc_Nit)
 
       if (f_bgc_Am(1:1) /= 'x') &
            call define_hist_field(n_bgc_Am,"bgc_Am","mmol/m^3",tstr3Da, tcstr, &
                "bulk ammonia/um ", "snow+bio grid", c1, c0,  &
                ns, f_bgc_Am)

       if (f_bgc_N(1:1) /= 'x') then
         do n=1,n_algae
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'bgc_N', trim(nchar)
            call define_hist_field(n_bgc_N(n,:),vname_in,"mmol/m^3",tstr3Da, tcstr, &
                "bulk algal N conc. ", "snow+bio grid", c1, c0, &
                ns, f_bgc_N)
         enddo
       endif
       if (f_bgc_C(1:1) /= 'x') then
         do n=1,n_algae
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'bgc_C', trim(nchar)
            call define_hist_field(n_bgc_C(n,:),vname_in,"mmol C/m^3",tstr3Da, tcstr, &
                "bulk algal C conc. ", "snow+bio grid", c1, c0, &
                ns, f_bgc_C)
         enddo
       endif
       if (f_bgc_chl(1:1) /= 'x') then
         do n=1,n_algae
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'bgc_chl', trim(nchar)
            call define_hist_field(n_bgc_chl(n,:),vname_in,"mg chl/m^3",tstr3Da, tcstr, &
                "bulk algal chl conc. ", "snow+bio grid", c1, c0, &
                ns, f_bgc_chl)
         enddo
       endif
       if (f_bgc_DOC(1:1) /= 'x') then
         do n=1,n_doc
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'bgc_DOC', trim(nchar)
            call define_hist_field(n_bgc_DOC(n,:),vname_in,"mmol/m^3",tstr3Da, tcstr, &
                "bulk DOC conc. ", "snow+bio grid", c1, c0, &
                ns, f_bgc_DOC)
         enddo
       endif
       if (f_bgc_DIC(1:1) /= 'x') then
         do n=1,n_dic
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'bgc_DIC', trim(nchar)
            call define_hist_field(n_bgc_DIC(n,:),vname_in,"mmol/m^3",tstr3Da, tcstr, &
                "bulk DIC conc. ", "snow+bio grid", c1, c0, &
                ns, f_bgc_DIC)
         enddo
       endif
       if (f_bgc_DON(1:1) /= 'x') then
         do n=1,n_don
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'bgc_DON', trim(nchar)
            call define_hist_field(n_bgc_DON(n,:),vname_in,"mmol/m^3",tstr3Da, tcstr, &
                "bulk DON conc. ", "snow+bio grid", c1, c0, &
                ns, f_bgc_DON)
         enddo
       endif
       if (f_bgc_Fed (1:1) /= 'x') then
         do n=1,n_fed 
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'bgc_Fed', trim(nchar)
            call define_hist_field(n_bgc_Fed (n,:),vname_in,"umol/m^3",tstr3Da, tcstr, &
                "bulk dFe conc. ", "snow+bio grid", c1, c0, &
                ns, f_bgc_Fed )
         enddo
       endif
       if (f_bgc_Fep (1:1) /= 'x') then
         do n=1,n_fep 
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'bgc_Fep', trim(nchar)
            call define_hist_field(n_bgc_Fep (n,:),vname_in,"umol/m^3",tstr3Da, tcstr, &
                "bulk pFe conc. ", "snow+bio grid", c1, c0, &
                ns, f_bgc_Fep )
         enddo
       endif
     
       if (f_bgc_Sil(1:1) /= 'x') &
            call define_hist_field(n_bgc_Sil,"bgc_Sil","mmol/m^3",tstr3Da, tcstr, &
                "bulk silicate ", "snow+bio grid", c1, c0, &
                ns, f_bgc_Sil)
     
       if (f_bgc_hum(1:1) /= 'x') &
            call define_hist_field(n_bgc_hum,"bgc_hum","mmol/m^3",tstr3Da, tcstr, &
                "bulk humic (carbon) material ", "snow+bio grid", c1, c0, &
                ns, f_bgc_hum)
      
       if (f_bgc_DMSPp(1:1) /= 'x') &
            call define_hist_field(n_bgc_DMSPp,"bgc_DMSPp","mmol/m^3",tstr3Da, tcstr, &
                "bulk algal DMSP ", "snow+bio grid", c1, c0,&
                ns, f_bgc_DMSPp)
      
       if (f_bgc_DMSPd(1:1) /= 'x') &
            call define_hist_field(n_bgc_DMSPd,"bgc_DMSPd","mmol/m^3",tstr3Da, tcstr, &
                "bulk dissolved DMSP ", "snow+bio grid", c1, c0, &
                ns, f_bgc_DMSPd)
  
       if (f_bgc_DMS(1:1) /= 'x') &
            call define_hist_field(n_bgc_DMS,"bgc_DMS","mmol/m^3",tstr3Da, tcstr, &
                "bulk DMS gas ", "snow+bio grid", c1, c0, &
                ns, f_bgc_DMS)
     
       if (f_bgc_PON(1:1) /= 'x') &
            call define_hist_field(n_bgc_PON,"bgc_PON","mmol/m^3",tstr3Da, tcstr, &
                "other bulk nitrogen pool ", "snow+bio grid", c1, c0, &
                ns, f_bgc_PON)

       !--------------------------------------------
       !   Category 1 BGC
       !----------------------------------------------

       if (f_bgc_Nit_cat1(1:1) /= 'x') & 
            call define_hist_field(n_bgc_Nit_cat1,"bgc_Nit_cat1","mmol/m^3",tstr3Da, tcstr, &
                "bulk nitrate in cat 1 ", "snow+bio grid", c1, c0,     &
                ns, f_bgc_Nit_cat1)
 
       if (f_bgc_Am_cat1(1:1) /= 'x') &
            call define_hist_field(n_bgc_Am_cat1,"bgc_Am_cat1","mmol/m^3",tstr3Da, tcstr, &
                "bulk ammonia/um in cat 1", "snow+bio grid", c1, c0,  &
                ns, f_bgc_Am_cat1)

       if (f_bgc_N_cat1(1:1) /= 'x') then
         do n=1,n_algae
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'bgc_N_cat1', trim(nchar)
            call define_hist_field(n_bgc_N_cat1(n,:),vname_in,"mmol/m^3",tstr3Da, tcstr, &
                "bulk algal N conc. in cat 1", "snow+bio grid", c1, c0, &
                ns, f_bgc_N_cat1)
         enddo
       endif
       if (f_bgc_DOC_cat1(1:1) /= 'x') then
         do n=1,n_doc
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'bgc_DOC_cat1', trim(nchar)
            call define_hist_field(n_bgc_DOC_cat1(n,:),vname_in,"mmol/m^3",tstr3Da, tcstr, &
                "bulk DOC conc. in cat 1 ", "snow+bio grid", c1, c0, &
                ns, f_bgc_DOC_cat1)
         enddo
       endif
       if (f_bgc_DIC_cat1(1:1) /= 'x') then
         do n=1,n_dic
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'bgc_DIC_cat1', trim(nchar)
            call define_hist_field(n_bgc_DIC_cat1(n,:),vname_in,"mmol/m^3",tstr3Da, tcstr, &
                "bulk DIC conc. in cat 1", "snow+bio grid", c1, c0, &
                ns, f_bgc_DIC_cat1)
         enddo
       endif
       if (f_bgc_DON_cat1(1:1) /= 'x') then
         do n=1,n_don
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'bgc_DON_cat1', trim(nchar)
            call define_hist_field(n_bgc_DON_cat1(n,:),vname_in,"mmol/m^3",tstr3Da, tcstr, &
                "bulk DON conc. in cat 1", "snow+bio grid", c1, c0, &
                ns, f_bgc_DON_cat1)
         enddo
       endif
       if (f_bgc_Fed_cat1 (1:1) /= 'x') then
         do n=1,n_fed 
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'bgc_Fed_cat1', trim(nchar)
            call define_hist_field(n_bgc_Fed_cat1 (n,:),vname_in,"umol/m^3",tstr3Da, tcstr, &
                "bulk dFe conc. in cat 1 ", "snow+bio grid", c1, c0, &
                ns, f_bgc_Fed_cat1 )
         enddo
       endif
       if (f_bgc_Fep_cat1 (1:1) /= 'x') then
         do n=1,n_fep 
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'bgc_Fep_cat1', trim(nchar)
            call define_hist_field(n_bgc_Fep_cat1 (n,:),vname_in,"umol/m^3",tstr3Da, tcstr, &
                "bulk pFe conc. in cat 1", "snow+bio grid", c1, c0, &
                ns, f_bgc_Fep_cat1 )
         enddo
       endif
     
       if (f_bgc_Sil_cat1(1:1) /= 'x') &
            call define_hist_field(n_bgc_Sil_cat1,"bgc_Sil_cat1","mmol/m^3",tstr3Da, tcstr, &
                "bulk silicate in cat 1", "snow+bio grid", c1, c0, &
                ns, f_bgc_Sil_cat1)
     
       if (f_bgc_hum_cat1(1:1) /= 'x') &
            call define_hist_field(n_bgc_hum,"bgc_hum_cat1","mmol/m^3",tstr3Da, tcstr, &
                "bulk humic (carbon) material in cat 1", "snow+bio grid", c1, c0, &
                ns, f_bgc_hum_cat1)
      
       if (f_bgc_DMSPd_cat1(1:1) /= 'x') &
            call define_hist_field(n_bgc_DMSPd_cat1,"bgc_DMSPd_cat1","mmol/m^3",tstr3Da, tcstr, &
                "bulk dissolved DMSP in cat 1", "snow+bio grid", c1, c0, &
                ns, f_bgc_DMSPd_cat1)
  
       if (f_bgc_DMS_cat1(1:1) /= 'x') &
            call define_hist_field(n_bgc_DMS_cat1,"bgc_DMS_cat1","mmol/m^3",tstr3Da, tcstr, &
                "bulk DMS gas in cat 1", "snow+bio grid", c1, c0, &
                ns, f_bgc_DMS_cat1)
     
       if (f_bgc_PON_cat1(1:1) /= 'x') &
            call define_hist_field(n_bgc_PON_cat1,"bgc_PON_cat1","mmol/m^3",tstr3Da, tcstr, &
                "other bulk nitrogen pool in cat 1", "snow+bio grid", c1, c0, &
                ns, f_bgc_PON_cat1)
    
    enddo  !ns

    endif ! z_tracers

      end subroutine init_hist_bgc_3Da

!=======================================================================

! Initialize bgc fields written to history files
!
! authors: Nicole Jeffery, LANL
!          Elizabeth C. Hunke, LANL

      subroutine init_history_bgc

      use ice_arrays_column, only: PP_net, grow_net, hbri, &
          ice_bio_net, snow_bio_net, fbio_snoice, fbio_atmice, &
          fzsal, fzsal_g, zfswin 
      use ice_flux_bgc, only: flux_bio, flux_bio_ai, fnit, fsil, &
          famm, fdmsp, fdms, fhum, fdust, falgalN, fdoc, fdic, &
          fdon, ffep, ffed

      character(len=*), parameter :: subname = '(init_history_bgc)'

      PP_net        (:,:,:) = c0
      grow_net      (:,:,:) = c0
      hbri          (:,:,:) = c0
      flux_bio    (:,:,:,:) = c0
      flux_bio_ai (:,:,:,:) = c0
      ice_bio_net (:,:,:,:) = c0
      snow_bio_net(:,:,:,:) = c0
      fbio_snoice (:,:,:,:) = c0
      fbio_atmice (:,:,:,:) = c0
      fzsal         (:,:,:) = c0
      fzsal_g       (:,:,:) = c0
      zfswin    (:,:,:,:,:) = c0
      fnit          (:,:,:) = c0
      fsil          (:,:,:) = c0
      famm          (:,:,:) = c0
      fdmsp         (:,:,:) = c0
      fdms          (:,:,:) = c0
      fhum          (:,:,:) = c0
      fdust         (:,:,:) = c0
      falgalN     (:,:,:,:) = c0
      fdoc        (:,:,:,:) = c0
      fdic        (:,:,:,:) = c0
      fdon        (:,:,:,:) = c0
      ffep        (:,:,:,:) = c0
      ffed        (:,:,:,:) = c0

      end subroutine init_history_bgc

!=======================================================================

      end module ice_history_bgc

!=======================================================================
