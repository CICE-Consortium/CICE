!  SVN:$Id: ice_restart_shared.F90 607 2013-03-29 15:49:42Z eclare $
!=======================================================================

      module ice_restart_shared

      use ice_kinds_mod
      implicit none
      private
      public :: lenstr

      logical (kind=log_kind), public :: &
         restart    , &   ! if true, initialize using restart file instead of defaults
         restart_ext, &   ! if true, read/write extended grid (with ghost cells)
         use_restart_time ! if true, use time written in core restart file

      character (len=char_len), public :: &
         runtype           ! initial, continue, hybrid, branch

      character (len=char_len_long), public :: &
         restart_file  , & ! output file for restart dump
         restart_dir   , & ! directory name for restart dump
         runid             ! identifier for CESM coupled run or bering

      character (len=char_len_long), public :: &
         pointer_file      ! input pointer file for restarts

      character (len=char_len), public :: &
         restart_format    ! format of restart files 'nc'

      logical (kind=log_kind), public :: lcdf64

!=======================================================================

      contains

!=======================================================================

! Compute length of string by finding first non-blank
! character from the right.

      integer function lenstr(label)

      character*(*) label

      ! local variables

      integer (kind=int_kind) :: &
         length, & ! length of character string
         n         ! loop index

      length = len(label)
      do n=length,1,-1
        if( label(n:n) /= ' ' ) exit
      enddo
      lenstr = n

      end function lenstr

!=======================================================================

      end module ice_restart_shared

!=======================================================================
