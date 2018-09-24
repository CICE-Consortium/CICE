!=======================================================================

! Defines the global domain size and number of categories and layers.
! Code originally based on domain_size.F in POP
!
! author Elizabeth C. Hunke, LANL
! 2004: Block structure and snow parameters added by William Lipscomb
!       Renamed (used to be ice_model_size)
! 2006: Converted to free source form (F90) by Elizabeth Hunke
!       Removed hardwired sizes (NX...can now be set in compile scripts)

      module ice_domain_size

      use ice_kinds_mod

!=======================================================================

      implicit none
      private

      integer (kind=int_kind), public :: &
        max_blocks  , & ! max number of blocks per processor
        block_size_x, & ! size of block in first horiz dimension
        block_size_y, & ! size of block in second horiz dimension
        nx_global   , & ! i-axis size
        ny_global       ! j-axis size

      integer (kind=int_kind), parameter, public :: &
        ncat      = NICECAT   , & ! number of categories
        nilyr     = NICELYR   , & ! number of ice layers per category
        nslyr     = NSNWLYR   , & ! number of snow layers per category
        n_aero    = NTRAERO   , & ! number of aerosols in use
        n_zaero   = TRZAERO   , & ! number of z aerosols in use 
        n_algae   = TRALG     , & ! number of algae in use 
        n_doc     = TRDOC     , & ! number of DOC pools in use
        n_dic     = TRDIC     , & ! number of DIC pools in use
        n_don     = TRDON     , & ! number of DON pools in use
        n_fed     = TRFED     , & ! number of Fe  pools in use dissolved Fe
        n_fep     = TRFEP     , & ! number of Fe  pools in use particulate Fe
        nblyr     = NBGCLYR   , & ! number of bio/brine layers per category 
                                  ! maximum number of biology tracers + aerosols
                                  ! *** add to kscavz in icepack_zbgc_shared.F90 
        n_bgc     = (n_algae*2 + n_doc + n_dic + n_don + n_fed + n_fep + n_zaero &
                  + 8)        , & ! nit, am, sil, dmspp, dmspd, dms, pon, humic 
        nltrcr    = (n_bgc*TRBGCZ+TRZS)*TRBRI, & ! number of zbgc (includes zaero)
                                                 ! and zsalinity tracers 
        max_nsw   = (nilyr+nslyr+2) & ! total chlorophyll plus aerosols
                  * (1+TRZAERO),& ! number of tracers active in shortwave calculation
        max_ntrcr =   1         & ! 1 = surface temperature              
                  + nilyr       & ! ice salinity
                  + nilyr       & ! ice enthalpy
                  + nslyr       & ! snow enthalpy
                              !!!!! optional tracers:
                  + TRAGE       & ! age
                  + TRFY        & ! first-year area
                  + TRLVL*2     & ! level/deformed ice
                  + TRPND*3     & ! ponds
                  + n_aero*4    & ! number of aerosols * 4 aero layers
                  + TRBRI       & ! brine height 
                  + TRBGCS*n_bgc           & ! skeletal layer BGC 
                  + TRZS  *TRBRI* nblyr    & ! zsalinity  (off if TRBRI=0)
                  + n_bgc*TRBGCZ*TRBRI*(nblyr+3) & ! zbgc (off if TRBRI=0) 
                  + n_bgc*TRBGCZ           & ! mobile/stationary phase tracer 
                  + 1         , & ! for unused tracer flags
        max_nstrm =   5           ! max number of history output streams

   !*** The model will inform the user of the correct
   !*** values for the parameter below.  A value higher than
   !*** necessary will not cause the code to fail, but will
   !*** allocate more memory than is necessary.  A value that
   !*** is too low will cause the code to exit.  
   !*** A good initial guess is found using
   !*** max_blocks = (nx_global/block_size_x)*(ny_global/block_size_y)/
   !***               num_procs
 
!=======================================================================

      end module ice_domain_size

!=======================================================================
