!=======================================================================

      module ice_restart_shared

      use ice_kinds_mod

      implicit none
      private
      public :: lenstr

      logical (kind=log_kind), public :: &
         restart    , &   ! if true, initialize using restart file instead of defaults
         restart_ext, &   ! if true, read/write extended grid (with ghost cells)
         restart_coszen, &   ! if true, read/write coszen
         use_restart_time ! if true, use time written in core restart file

      character (len=char_len), public :: &
         runtype           ! initial, continue, hybrid, branch

      character (len=char_len_long), public :: &
         restart_file  , & ! output file for restart dump
         restart_dir   , & ! directory name for restart dump
         runid             ! identifier for CESM coupled run or bering

      character (len=char_len_long), public :: &
         pointer_file      ! input pointer file for restarts

      logical (kind=log_kind), public :: &
         pointer_date =  .false.   ! if true, append datestamp to pointer file

      character (len=char_len), public :: &
         restart_format      , & ! format of restart files 'nc'
         restart_mod         , & ! restart modification option, "none", "adjust_aice"
         restart_rearranger      ! restart file rearranger, box or subset for pio

      integer (kind=int_kind), public :: &
         restart_iotasks     , & ! iotasks, root, stride defines io pes for pio
         restart_root        , & ! iotasks, root, stride defines io pes for pio
         restart_stride      , & ! iotasks, root, stride defines io pes for pio
         restart_deflate     , & ! compression level for hdf5/netcdf4
         restart_chunksize(2)    ! chunksize for hdf5/netcdf4


!=======================================================================

      contains

!=======================================================================

! Compute length of string by finding first non-blank
! character from the right.

      integer function lenstr(label)

      character(len=*) :: label

      character(len=*),parameter :: subname='(lenstr)'

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
