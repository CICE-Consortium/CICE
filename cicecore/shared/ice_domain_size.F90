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

      ! namelist

      integer (kind=int_kind), public :: &
        max_blocks  , & ! max number of blocks per processor
        block_size_x, & ! size of block in first horiz dimension
        block_size_y, & ! size of block in second horiz dimension
        nx_global   , & ! i-axis size
        ny_global       ! j-axis size

      integer (kind=int_kind), public :: &
        ncat      , & ! number of categories
        nilyr     , & ! number of ice layers per category
        nslyr     , & ! number of snow layers per category
        nblyr     , & ! number of bio/brine layers per category 
        n_aero    , & ! number of aerosols in use
        n_zaero   , & ! number of z aerosols in use 
        n_algae   , & ! number of algae in use 
        n_doc     , & ! number of DOC pools in use
        n_dic     , & ! number of DIC pools in use
        n_don     , & ! number of DON pools in use
        n_fed     , & ! number of Fe  pools in use dissolved Fe
        n_fep         ! number of Fe  pools in use particulate Fe

      integer (kind=int_kind), public, parameter :: &
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
