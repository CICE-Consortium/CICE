
      module optargs_subs

      implicit none
      private

      integer, public, parameter :: dp = kind(1.d0)

      logical, public :: computeA = .false., &
                         computeB = .false., &
                         computeC = .false., &
                         computeD = .false., &
                         computeE = .false.

      integer, public :: oa_error = -99, &
                         oa_OK    =   0, &
                         oa_A     =   1, &
                         oa_B     =   2, &
                         oa_C     =   4, &
                         oa_D     =   8, &
                         oa_E     =  16

      public :: oa_layer1, oa_count1

!-----------------------------------
CONTAINS
!-----------------------------------

      subroutine oa_count1(Ai1,Ao,B,Ci1,Co,Di1,Di2,Do,Ei,Eo,ierr)

      real(dp), intent(in)   , optional :: Ai1, Di1, Di2
      real(dp), intent(out)  , optional :: Ao, Do
      real(dp), intent(inout), optional :: B
      real(dp), intent(in)              :: Ci1
      real(dp), intent(out)             :: Co
      real(dp), intent(in)   , optional, dimension(:) :: Ei
      real(dp), intent(out)  , optional, dimension(:) :: Eo
      integer , intent(inout)           :: ierr

      call oa_count2(Ai1,Ao,B,Ci1,Co,Di1,Di2,Do,Ei,Eo,ierr)

!      write(6,*) 'debug oa_count1 ',ierr

      end subroutine oa_count1

!-----------------------------------

      subroutine oa_count2(Ai1,Ao,B,Ci1,Co,Di1,Di2,Do,Ei,Eo,ierr)

      real(dp), intent(in)   , optional :: Ai1, Di1, Di2
      real(dp), intent(out)  , optional :: Ao, Do
      real(dp), intent(inout), optional :: B
      real(dp), intent(in)              :: Ci1
      real(dp), intent(out)             :: Co
      real(dp), intent(in)   , optional, dimension(:) :: Ei
      real(dp), intent(out)  , optional, dimension(:) :: Eo
      integer , intent(inout)           :: ierr

      ierr = 3 ! Ci1, Co, ierr have to be passed
      if (present(Ai1)) ierr = ierr + 1
      if (present(Ao) ) ierr = ierr + 1
      if (present(B)  ) ierr = ierr + 1
      if (present(Di1)) ierr = ierr + 1
      if (present(Di2)) ierr = ierr + 1
      if (present(Do) ) ierr = ierr + 1
      if (present(Ei) ) ierr = ierr + 1
      if (present(Eo) ) ierr = ierr + 1

!      write(6,*) 'debug oa_count2 ',ierr

      end subroutine oa_count2

!-----------------------------------

      subroutine oa_layer1(Ai1,Ao,B,Ci1,Co,Di1,Di2,Do,Ei,Eo,ierr)

      real(dp), intent(in)   , optional :: Ai1, Di1, Di2
      real(dp), intent(out)  , optional :: Ao, Do
      real(dp), intent(inout), optional :: B
      real(dp), intent(in)              :: Ci1
      real(dp), intent(out)             :: Co
      real(dp), intent(in)   , optional, dimension(:) :: Ei
      real(dp), intent(out)  , optional, dimension(:) :: Eo
      integer , intent(inout)           :: ierr

      ierr = oa_OK
      if (computeA) then
         if (.not.(present(Ai1).and.present(Ao))) then
            ierr = oa_error
         endif
      endif
      if (computeB) then
         if (.not.(present(B))) then
            ierr = oa_error
         endif
      endif
      if (computeD) then
         if (.not.(present(Di1).and.present(Di2).and.present(Do))) then
            ierr = oa_error
         endif
      endif
      if (computeE) then
         if (.not.(present(Ei).and.present(Eo))) then
            ierr = oa_error
         endif
      endif

      if (ierr == oa_OK) then
         call oa_layer2(Ai1,Ao,B,Ci1,Co,Di1,Di2,Do,Ei,Eo,ierr)
      endif

      end subroutine oa_layer1

!-----------------------------------

      subroutine oa_layer2(Ai1,Ao,B,Ci1,Co,Di1,Di2,Do,Ei,Eo,ierr)

! Note: optional arrays must have an optional attribute, otherwise they seg fault
! Scalars do not seem to have this problem

      real(dp), intent(in)              :: Ai1, Di1, Di2
      real(dp), intent(out)             :: Ao, Do
      real(dp), intent(inout)           :: B
      real(dp), intent(in)              :: Ci1
      real(dp), intent(out)             :: Co
      real(dp), intent(in)   , optional, dimension(:) :: Ei
      real(dp), intent(out)  , optional, dimension(:) :: Eo
      integer , intent(inout)           :: ierr

      call oa_compute(Ai1,Ao,B,Ci1,Co,Di1,Di2,Do,Ei,Eo,ierr)

      end subroutine oa_layer2

!-----------------------------------

      subroutine oa_compute(Ai1,Ao,B,Ci1,Co,Di1,Di2,Do,Ei,Eo,ierr)

! Note: optional arrays must have an optional attribute, otherwise they seg fault
! Scalars do not seem to have this problem

      real(dp), intent(in)              :: Ai1, Di1, Di2
      real(dp), intent(out)             :: Ao, Do
      real(dp), intent(inout)           :: B
      real(dp), intent(in)              :: Ci1
      real(dp), intent(out)             :: Co
      real(dp), intent(in)   , optional, dimension(:) :: Ei
      real(dp), intent(out)  , optional, dimension(:) :: Eo
      integer , intent(inout)           :: ierr

      integer :: n

      if (computeA) then
         Ao = Ai1 - 1.
         ierr = ierr + oa_A
      endif

      if (computeB) then
         B = B + 5.
         ierr = ierr + oa_B
      endif

      if (computeC) then
         Co = Ci1 * (2.)
         ierr = ierr + oa_C
      endif

      if (computeD) then
         Do = Di1 + Di2
         ierr = ierr + oa_D
      endif

      if (computeE) then
         ierr = ierr + oa_E
         do n = 1,size(Eo)
           Eo(n) = Ei(n) + n
         enddo
      endif

      return
      end subroutine oa_compute

!-----------------------------------

      end module optargs_subs
