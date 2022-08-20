! Provides methods for querying memory use

MODULE ice_memusage

!-------------------------------------------------------------------------------
! PURPOSE: memory use query methods
!    Should call ice_memusage_init once before calling other interfaces
!-------------------------------------------------------------------------------

   use ice_kinds_mod, only : dbl_kind, log_kind

   implicit none
   private

! PUBLIC: Public interfaces

   public ::  ice_memusage_getusage, &
              ice_memusage_init, &
              ice_memusage_print

   logical(log_kind), public :: memory_stats

! PRIVATE DATA:

   real(dbl_kind) :: mb_blk = 1.0_dbl_kind
   logical        :: initset = .false.

!===============================================================================

contains

!===============================================================================
! Initialize memory conversion to MB

subroutine ice_memusage_init(iunit)

   implicit none

   !----- arguments -----

   integer, optional :: iunit   !< output unit number for optional writes

   !----- local -----

   ! --- Memory stats ---
   integer :: msize                   ! memory size (high water)
   integer :: mrss0,mrss1,mrss2       ! temporary rss
   integer :: mshare,mtext,mdatastack
   integer :: ierr

   integer :: ice_memusage_gptl

   real(dbl_kind),allocatable :: mem_tmp(:)
   character(*),parameter  :: subname = '(ice_memusage_init)'

   !---------------------------------------------------

   ! return if memory_stats are off
   if (.not. memory_stats) return

   ierr = ice_memusage_gptl (msize, mrss0, mshare, mtext, mdatastack)
   allocate(mem_tmp(1024*1024))    ! 1 MWord, 8 MB
   mem_tmp = -1.0
   ierr = ice_memusage_gptl (msize, mrss1, mshare, mtext, mdatastack)
   deallocate(mem_tmp)
   ierr = ice_memusage_gptl (msize, mrss2, mshare, mtext, mdatastack)
   mb_blk = 1.0_dbl_kind
   if (mrss1 - mrss0 > 0) then
      mb_blk = (8.0_dbl_kind)/((mrss1-mrss0)*1.0_dbl_kind)
      initset = .true.
   endif

   if (present(iunit)) then
      write(iunit,'(A,l4)')    subname//' Initset conversion flag is ',initset
      write(iunit,'(A,f16.2)') subname//' 8 MB memory   alloc in MB is ',(mrss1-mrss0)*mb_blk
      write(iunit,'(A,f16.2)') subname//' 8 MB memory dealloc in MB is ',(mrss1-mrss2)*mb_blk
      write(iunit,'(A,f16.2)') subname//' Memory block size conversion in bytes is ',mb_blk*1024_dbl_kind*1024.0_dbl_kind
   endif

end subroutine ice_memusage_init

!===============================================================================
! Determine memory use

subroutine ice_memusage_getusage(r_msize,r_mrss)

   implicit none

   !----- arguments ---
   real(dbl_kind),intent(out) :: r_msize  !< memory usage value
   real(dbl_kind),intent(out) :: r_mrss   !< memory usage value

   !----- local ---
   integer :: msize,mrss
   integer :: mshare,mtext,mdatastack
   integer :: ierr
   integer :: ice_memusage_gptl
   character(*),parameter  :: subname = '(ice_memusage_getusage)'

   !---------------------------------------------------

   ! return if memory_stats are off
   if (.not. memory_stats) return

   ierr = ice_memusage_gptl (msize, mrss, mshare, mtext, mdatastack)
   r_msize = msize*mb_blk
   r_mrss  = mrss*mb_blk

end subroutine ice_memusage_getusage

!===============================================================================
! Print memory use

subroutine ice_memusage_print(iunit,string)

   implicit none

   !----- arguments ---
   integer, intent(in) :: iunit    !< unit number to write to
   character(len=*),optional, intent(in) :: string  !< optional string

   !----- local ---
   real(dbl_kind)     :: msize,mrss
   character(len=128) :: lstring
   character(*),parameter  :: subname = '(ice_memusage_print)'

   !---------------------------------------------------

   ! return if memory_stats are off
   if (.not. memory_stats) return

   lstring = ' '
   if (present(string)) then
      lstring = string
   endif

   call ice_memusage_getusage(msize,mrss)

   if (initset) then
      write(iunit,'(2a,2f14.4,1x,a)') subname,' memory use (MB) = ',msize,mrss,trim(lstring)
   else
      write(iunit,'(2a,2f14.4,1x,a)') subname,' memory use (??) = ',msize,mrss,trim(lstring)
   endif

end subroutine ice_memusage_print

!===============================================================================

END MODULE ice_memusage
