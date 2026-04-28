
      program mpif08

      ! This tests the availability of mpi_f08 on porting platforms

#define USE_MPI_F08
#ifdef USE_MPI
      use mpi
#else
      use mpi_f08
#endif

      implicit none

      integer, parameter :: &
                            int_kind  = selected_int_kind(6), &
                            real_kind = selected_real_kind(6), &
                            dbl_kind  = selected_real_kind(13)

      integer(int_kind) :: i,j,k,n

      integer(int_kind), parameter :: dsize = 8
      integer(int_kind) :: arrayI(dsize,dsize)
      real(real_kind)   :: array4(dsize,dsize)
      real(dbl_kind)    :: array8(dsize,dsize)
      real(dbl_kind)    :: carray8(dsize,dsize),rarray8(dsize,dsize)

      integer(int_kind) :: my_task,master_task,npes
      integer(int_kind) :: ierr,tag
      integer(int_kind) :: sro ! send recv pe offset
#ifdef USE_MPI
      integer :: MPI_COMM_ICE
      integer :: rrequest
      integer :: rstatus(MPI_STATUS_SIZE)
#else
      type(MPI_COMM)    :: MPI_COMM_ICE
      type(MPI_REQUEST) :: rrequest
      type(MPI_STATUS)  :: rstatus
#endif

      integer(int_kind), parameter :: ntests1 = 4
      integer(int_kind) :: iflag, gflag
      integer(int_kind) :: cval,ival,sumi
      real(real_kind)   :: sum4
      real(dbl_kind)    :: sum8
      character(len=8)  :: errorflag0
      character(len=8)  :: errorflag1(ntests1)
      character(len=32) :: stringflag1(ntests1)

      character(len=*), parameter ::  &
         passflag = 'PASS', &
         failflag = 'FAIL'

      character(len=*), parameter :: subname = '(mpif08)'

      ! ---------------------------

      call MPI_INIT(ierr)
      if (ierr /= MPI_SUCCESS) write(6,*) ' ERROR MPI_INIT'
      call MPI_COMM_DUP(MPI_COMM_WORLD, MPI_COMM_ICE, ierr)
      if (ierr /= MPI_SUCCESS) write(6,*) ' ERROR MPI_COMM_DUP'

      master_task = 0
      call MPI_COMM_SIZE(MPI_COMM_ICE, npes, ierr)
      if (ierr /= MPI_SUCCESS) write(6,*) ' ERROR MPI_COMM_SIZE'
      call MPI_COMM_RANK(MPI_COMM_ICE, my_task, ierr)
      if (ierr /= MPI_SUCCESS) write(6,*) ' ERROR MPI_COMM_RANK'
      sro = max(npes/2,1)

      if (my_task == master_task) then
         write(6,*) ' '
         write(6,*) '=========================================================='
         write(6,*) ' '
         write(6,*) 'RunningUnitTest MPIF08'
#ifdef USE_MPI
         write(6,*) ' with    use mpi'
#else
         write(6,*) ' with    use mpi_f08'
#endif
         write(6,*) ' '
         write(6,*) ' npes         = ',npes
         write(6,*) ' my_task      = ',my_task
         write(6,*) ' sr pe offset = ',sro
         write(6,*) ' '
         flush(6)
      endif

      ! ---------------------------

      errorflag0 = passflag
      errorflag1 = passflag
      stringflag1 = ' '

      ival = 0
      do k = 0,npes-1
      do j = 1,dsize
      do i = 1,dsize
         cval = (k*100 + j*10 + i)
         ival = ival + cval
      enddo
      enddo
      enddo

      do j = 1,dsize
      do i = 1,dsize
         cval = (my_task*100 + j*10 + i)
         arrayI(i,j) = cval
         array4(i,j) = real(cval,real_kind)
         array8(i,j) = real(cval,dbl_kind)
      enddo
      enddo

      ! for test 4, communicate array to neighbor task
      do j = 1,dsize
      do i = 1,dsize
         cval = (mod(my_task+sro,npes)*100 + j*10 + i)
         carray8(i,j) = real(cval,dbl_kind)
      enddo
      enddo

!      if (my_task == master_task) write(6,*) 'ival = ',ival

      do n = 1,ntests1
         iflag = 0
         gflag = 0
         if (n == 1) then
            stringflag1(n) = 'allreduce_int'
            call MPI_ALLREDUCE(sum(arrayI), sumi, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_ICE, ierr)
            if (ierr /= MPI_SUCCESS) write(6,*) ' ERROR MPI_ALLREDUCE arrayI'
            if (sumi /= ival) iflag=1
            if (my_task == master_task) write(6,*) ' Test ',n,ival,sumi
            flush(6)
         elseif (n == 2) then
            stringflag1(n) = 'allreduce_r4'
            call MPI_ALLREDUCE(sum(array4), sum4, 1, MPI_REAL, MPI_SUM, MPI_COMM_ICE, ierr)
            if (ierr /= MPI_SUCCESS) write(6,*) ' ERROR MPI_ALLREDUCE array4'
            if (sum4 /= real(ival,real_kind)) iflag=1
            if (my_task == master_task) write(6,*) ' Test ',n,ival,sum4
            flush(6)
         elseif (n == 3) then
            stringflag1(n) = 'allreduce_r8'
            call MPI_ALLREDUCE(sum(array8), sum8, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_ICE, ierr)
            if (ierr /= MPI_SUCCESS) write(6,*) ' ERROR MPI_ALLREDUCE array8'
            if (sum8 /= real(ival,dbl_kind)) iflag=1
            if (my_task == master_task) write(6,*) ' Test ',n,ival,sum8
            flush(6)
         elseif (n == 4) then
            stringflag1(n) = 'irecv_send_r8'
            tag = 1001
            call MPI_IRECV(rarray8, dsize*dsize, MPI_DOUBLE, mod(my_task+sro,npes), tag, MPI_COMM_ICE, rrequest, ierr)
            if (ierr /= MPI_SUCCESS) write(6,*) ' ERROR MPI_IRECV'
            call MPI_SEND(array8, dsize*dsize, MPI_DOUBLE, mod(my_task-sro+npes,npes), tag, MPI_COMM_ICE, ierr)
            if (ierr /= MPI_SUCCESS) write(6,*) ' ERROR MPI_SEND'
            call MPI_WAIT(rrequest, rstatus, ierr)
            if (ierr /= MPI_SUCCESS) write(6,*) ' ERROR MPI_WAIT'
            do j = 1,dsize
            do i = 1,dsize
               if (rarray8(i,j) /= carray8(i,j)) iflag=1
            enddo
            enddo
            if (my_task == master_task) write(6,*) ' Test ',n,sum(carray8),sum(rarray8)
            flush(6)
         endif
         call MPI_REDUCE(iflag, gflag, 1, MPI_INTEGER, MPI_MAX, master_task, MPI_COMM_ICE, ierr)
         if (ierr /= MPI_SUCCESS) write(6,*) ' ERROR MPI_REDUCE iflag'
         if (my_task == master_task) then
            if (gflag /= 0) then
               write(6,*) '   **** ERROR test ',n
               flush(6)
               errorflag1(n) = failflag
               errorflag0 = failflag
            endif
         endif
      enddo

      call MPI_BARRIER(MPI_COMM_ICE, ierr)
      if (ierr /= MPI_SUCCESS) write(6,*) ' ERROR MPI_BARRIER 2'

      ! ---------------------------

      if (my_task == master_task) then
         write(6,*) ' '
         do k = 1,ntests1
            write(6,*) errorflag1(k),stringflag1(k)
         enddo
         write(6,*) ' '
         write(6,*) 'MPIF08 COMPLETED SUCCESSFULLY'
         if (errorflag0 == passflag) then
            write(6,*) 'MPIF08 TEST COMPLETED SUCCESSFULLY'
         else
            write(6,*) 'MPIF08 TEST FAILED'
         endif
      endif

      flush(6)

      ! ---------------------------
      ! exit gracefully

      call MPI_FINALIZE(ierr)
      if (ierr /= MPI_SUCCESS) write(6,*) ' ERROR MPI_FINALIZE'

      end program mpif08

!=======================================================================
