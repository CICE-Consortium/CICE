!=======================================================================
!
! Primary state variables in various configurations
! Note: other state variables are at the end of this...
! The primary state variable names are:
!-------------------------------------------------------------------
! for each category   aggregated over     units
!                       categories
!-------------------------------------------------------------------
! aicen(i,j,n)         aice(i,j)           ---
! vicen(i,j,n)         vice(i,j)           m
! vsnon(i,j,n)         vsno(i,j)           m
! trcrn(i,j,it,n)      trcr(i,j,it)
!
! Area is dimensionless because aice is the fractional area
! (normalized so that the sum over all categories, including open
! water, is 1.0).  That is why vice/vsno have units of m instead of m^3.
!
! Variable names follow these rules:
!
! (1) For 3D variables (indices i,j,n), write 'ice' or 'sno' or
!     'sfc' and put an 'n' at the end.
! (2) For 2D variables (indices i,j) aggregated over all categories,
!     write 'ice' or 'sno' or 'sfc' without the 'n'.
! (3) For 2D variables (indices i,j) associated with an individual
!     category, write 'i' or 's' instead of 'ice' or 'sno' and put an 'n'
!     at the end: e.g. hin, hsn.  These are not declared here
!     but in individual modules (e.g., icepack_therm_vertical).
!
! authors C. M. Bitz, UW
!         Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! 2004: Block structure added by William Lipscomb
! 2006: Converted to free form source (F90) by Elizabeth Hunke

      module ice_bound_state

      use ice_kinds_mod
      use ice_constants, only: field_loc_center, field_type_scalar, c0
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_diag
!      use icepack_intfc, only: icepack_query_tracer_sizes
!      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted

      implicit none
      private
      public :: bound_state

!=======================================================================

      contains

!=======================================================================
!
! Get ghost cell values for ice state variables in each thickness category.
! NOTE: This subroutine cannot be called from inside a block loop!
!
! author: William H. Lipscomb, LANL

      subroutine bound_state (aicen,        &
                              vicen, vsnon, &
                              ntrcr, trcrn, &
                              restore)

      use ice_boundary, only: ice_halo, ice_HaloMask, ice_HaloUpdate, &
          ice_HaloDestroy
      use ice_domain_size, only: max_blocks, ncat
      use ice_blocks, only: nx_block, ny_block
      use ice_domain, only: halo_info, maskhalo_bound, nblocks
      use ice_restoring, only: restore_ice, ice_HaloRestore

      integer (kind=int_kind), intent(in) :: &
         ntrcr     ! number of tracers in use

      real (kind=dbl_kind), dimension(nx_block,ny_block,ncat,max_blocks), intent(inout) :: &
         aicen , & ! fractional ice area
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), intent(inout), dimension(:,:,:,:,:) :: &  ! (nx_block,ny_block,ntrcr,ncat,max_blocks)
         trcrn     ! ice tracers

      logical (kind=log_kind), intent(in), optional :: &
         restore

      ! local variables

      integer (kind=int_kind) :: i, j, n, iblk

      integer (kind=int_kind), &
         dimension(nx_block,ny_block,max_blocks) :: halomask

      type (ice_halo) :: halo_info_aicemask  ! halo mask

      logical (kind=log_kind) :: lrestore  ! local restart flag

      character(len=*), parameter :: subname = '(bound_state)'

      if (present(restore)) then
         lrestore = restore
      else
         lrestore = .true.
      endif

      call ice_HaloUpdate (aicen,            halo_info, &
                           field_loc_center, field_type_scalar)

      if (maskhalo_bound) then
         halomask(:,:,:) = 0

         !$OMP PARALLEL DO PRIVATE(iblk,n,i,j)
         do iblk = 1, nblocks
         do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            if (aicen(i,j,n,iblk) > c0) halomask(i,j,iblk) = 1
         enddo
         enddo
         enddo
         enddo
         !$OMP END PARALLEL DO

         call ice_HaloMask(halo_info_aicemask, halo_info, halomask)

         call ice_HaloUpdate (trcrn(:,:,:,:,:), halo_info_aicemask, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vicen,            halo_info_aicemask, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vsnon,            halo_info_aicemask, &
                              field_loc_center, field_type_scalar)
         call ice_HaloDestroy(halo_info_aicemask)

      else
         call ice_HaloUpdate (trcrn(:,:,:,:,:), halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vicen,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vsnon,            halo_info, &
                              field_loc_center, field_type_scalar)
      endif

      if (lrestore .and. restore_ice) then
         call ice_HaloRestore(fields='state')
      endif

      end subroutine bound_state

!=======================================================================

      end module ice_bound_state

!=======================================================================
