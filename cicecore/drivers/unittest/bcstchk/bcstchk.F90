
      program bcstchk

      ! This tests the CICE ice_broadcast infrastructure by calling the
      ! methods with hardwired input and known outputs and verifying the
      ! results.

      use ice_kinds_mod, only: int_kind, dbl_kind, real_kind, log_kind
      use ice_communicate, only: my_task, master_task, get_num_procs, get_rank, MPI_COMM_ICE
      use ice_communicate, only: init_communicate, get_num_procs, ice_barrier
      use ice_global_reductions, only: global_maxval
      use ice_fileunits, only: flush_fileunit
      use ice_exit, only: abort_ice, end_run
      use ice_broadcast

      implicit none

      integer(int_kind) :: n, k, k1, k2, k3

      integer(int_kind), parameter :: dsize = 10
      integer(int_kind) :: ival, i0, i1(dsize), i2(dsize,dsize), i3(dsize,dsize,dsize)
      logical(log_kind) :: lval, l0, l1(dsize), l2(dsize,dsize), l3(dsize,dsize,dsize)
      real(real_kind)   :: rval, r0, r1(dsize), r2(dsize,dsize), r3(dsize,dsize,dsize)
      real(dbl_kind)    :: dval, d0, d1(dsize), d2(dsize,dsize), d3(dsize,dsize,dsize)
      character(len=32) :: cval, c0

      real(dbl_kind)    :: xval

      integer(int_kind), parameter :: ntests1 = 17
      character(len=8)  :: errorflag1(ntests1)
      character(len=32) :: stringflag1(ntests1)

      integer(int_kind) :: ierr, npes, bcst_pe
      integer(int_kind) :: iflag, gflag
      character(len=8)  :: errorflag0
      character(len=16)  :: teststr
      character(len=*), parameter ::  &
         passflag = 'PASS', &
         failflag = 'FAIL'

      character(len=*), parameter :: subname = '(bcstchk)'

      ! ---------------------------

      call init_communicate()
      npes = get_num_procs()
      my_task = get_rank()
      master_task = 0

      if (my_task == master_task) then
         write(6,*) ' '
         write(6,*) '=========================================================='
         write(6,*) ' '
         write(6,*) 'RunningUnitTest BCSTCHK'
         write(6,*) ' '
         write(6,*) ' npes         = ',npes
         write(6,*) ' my_task      = ',my_task
         write(6,*) ' '
      endif

      errorflag0 = passflag
      errorflag1 = passflag
      stringflag1 = ' '

      ! ---------------------------
      ! Test ice_broadcast methods
      ! Test broadcast from root and from npes
      ! ---------------------------

      do k = 1,2
         if (k == 1) then
            bcst_pe = 0
         else
            bcst_pe = max(npes,1) - 1
         endif
         if (my_task == master_task) then
            write(6,*) ' '
            write(6,*) ' bcst_pe      = ',bcst_pe
         endif

         xval = -999._dbl_kind
         rval = 21.5_real_kind + real(bcst_pe,kind=real_kind)
         dval = 17.3_dbl_kind + real(bcst_pe,kind=dbl_kind)
         ival = 223 + bcst_pe
         write(cval,'(a,i4.4)') 'string is passed from ',bcst_pe
         lval = (k == 1)

         do n = 1,ntests1
            i0 = xval
            i1 = xval
            i2 = xval
            i3 = xval
            r0 = xval
            r1 = xval
            r2 = xval
            r3 = xval
            d0 = xval
            d1 = xval
            d2 = xval
            d3 = xval
            l0 = .not.lval
            l1 = .not.lval
            l2 = .not.lval
            l3 = .not.lval
            c0 = 'nothing to see here'

            if (my_task == bcst_pe) then
               i0 = ival
               i1 = ival
               i2 = ival
               i3 = ival
               r0 = rval
               r1 = rval
               r2 = rval
               r3 = rval
               d0 = dval
               d1 = dval
               d2 = dval
               d3 = dval
               l0 = lval
               l1 = lval
               l2 = lval
               l3 = lval
               c0 = cval
            endif

            iflag = 0
            gflag = -1
            write(teststr,'(a,1x,i2.2)') '   test',n
            if (n == 1) then
               stringflag1(n) = ' bcst_scalar_dbl'
               call broadcast_scalar(d0, bcst_pe)
               if (my_task == master_task) write(6,*) trim(teststr),trim(stringflag1(n)),d0,dval
               if (d0 /= dval) iflag=1
            elseif (n == 2) then
               stringflag1(n) = ' bcst_array_dbl_1d'
               call broadcast_array(d1, bcst_pe)
               if (my_task == master_task) write(6,*) trim(teststr),trim(stringflag1(n)),minval(d1),maxval(d1),dval
               if (minval(d1) /= dval) iflag=1
               if (maxval(d1) /= dval) iflag=1
            elseif (n == 3) then
               stringflag1(n) = ' bcst_array_dbl_2d'
               call broadcast_array(d2, bcst_pe)
               if (my_task == master_task) write(6,*) trim(teststr),trim(stringflag1(n)),minval(d2),maxval(d2),dval
               if (minval(d2) /= dval) iflag=1
               if (maxval(d2) /= dval) iflag=1
            elseif (n == 4) then
               stringflag1(n) = ' bcst_array_dbl_3d'
               call broadcast_array(d3, bcst_pe)
               if (my_task == master_task) write(6,*) trim(teststr),trim(stringflag1(n)),minval(d3),maxval(d3),dval
               if (minval(d3) /= dval) iflag=1
               if (maxval(d3) /= dval) iflag=1
            elseif (n == 5) then
               stringflag1(n) = ' bcst_scalar_real'
               call broadcast_scalar(r0, bcst_pe)
               if (my_task == master_task) write(6,*) trim(teststr),trim(stringflag1(n)),r0,rval
               if (r0 /= rval) iflag=1
            elseif (n == 6) then
               stringflag1(n) = ' bcst_array_real_1d'
               call broadcast_array(r1, bcst_pe)
               if (my_task == master_task) write(6,*) trim(teststr),trim(stringflag1(n)),minval(r1),maxval(r1),rval
               if (minval(r1) /= rval) iflag=1
               if (maxval(r1) /= rval) iflag=1
            elseif (n == 7) then
               stringflag1(n) = ' bcst_array_real_2d'
               call broadcast_array(r2, bcst_pe)
               if (my_task == master_task) write(6,*) trim(teststr),trim(stringflag1(n)),minval(r2),maxval(r2),rval
               if (minval(r2) /= rval) iflag=1
               if (maxval(r2) /= rval) iflag=1
            elseif (n == 8) then
               stringflag1(n) = ' bcst_array_real_3d'
               call broadcast_array(r3, bcst_pe)
               if (my_task == master_task) write(6,*) trim(teststr),trim(stringflag1(n)),minval(r3),maxval(r3),rval
               if (minval(r3) /= rval) iflag=1
               if (maxval(r3) /= rval) iflag=1
            elseif (n == 9) then
               stringflag1(n) = ' bcst_scalar_int'
               call broadcast_scalar(i0, bcst_pe)
               if (my_task == master_task) write(6,*) trim(teststr),trim(stringflag1(n)),i0,ival
               if (i0 /= ival) iflag=1
            elseif (n == 10) then
               stringflag1(n) = ' bcst_array_int_1d'
               call broadcast_array(i1, bcst_pe)
               if (my_task == master_task) write(6,*) trim(teststr),trim(stringflag1(n)),minval(i1),maxval(i1),ival
               if (minval(i1) /= ival) iflag=1
               if (maxval(i1) /= ival) iflag=1
            elseif (n == 11) then
               stringflag1(n) = ' bcst_array_int_2d'
               call broadcast_array(i2, bcst_pe)
               if (my_task == master_task) write(6,*) trim(teststr),trim(stringflag1(n)),minval(i2),maxval(i2),ival
               if (minval(i2) /= ival) iflag=1
               if (maxval(i2) /= ival) iflag=1
            elseif (n == 12) then
               stringflag1(n) = ' bcst_array_int_3d'
               call broadcast_array(i3, bcst_pe)
               if (my_task == master_task) write(6,*) trim(teststr),trim(stringflag1(n)),minval(i3),maxval(i3),ival
               if (minval(i3) /= ival) iflag=1
               if (maxval(i3) /= ival) iflag=1
            elseif (n == 13) then
               stringflag1(n) = ' bcst_scalar_logical'
               call broadcast_scalar(l0, bcst_pe)
               if (my_task == master_task) write(6,*) trim(teststr),trim(stringflag1(n)),l0,lval
               if (l0 .neqv. lval) iflag=1
            elseif (n == 14) then
               stringflag1(n) = ' bcst_array_logical_1d'
               call broadcast_array(l1, bcst_pe)
               if (my_task == master_task) write(6,*) trim(teststr),trim(stringflag1(n)),l1(1),lval
               do k1 = 1,dsize
                  if (l1(k1) .neqv. lval) iflag=1
               enddo
            elseif (n == 15) then
               stringflag1(n) = ' bcst_array_logical_2d'
               call broadcast_array(l2, bcst_pe)
               if (my_task == master_task) write(6,*) trim(teststr),trim(stringflag1(n)),l2(1,1),lval
               do k2 = 1,dsize
               do k1 = 1,dsize
                  if (l2(k1,k2) .neqv. lval) iflag=1
               enddo
               enddo
            elseif (n == 16) then
               stringflag1(n) = ' bcst_array_logical_3d'
               call broadcast_array(l3, bcst_pe)
               if (my_task == master_task) write(6,*) trim(teststr),trim(stringflag1(n)),l3(1,1,1),lval
               do k3 = 1,dsize
               do k2 = 1,dsize
               do k1 = 1,dsize
                  if (l3(k1,k2,k3) .neqv. lval) iflag=1
               enddo
               enddo
               enddo
            elseif (n == 17) then
               stringflag1(n) = ' bcst_scalar_char'
               call broadcast_scalar(c0, bcst_pe)
               if (my_task == master_task) write(6,*) trim(teststr),trim(stringflag1(n)),' ',trim(c0),' : ',trim(cval)
               if (c0 /= cval) iflag=1
            else
               call abort_ice(subname//' illegal k bcst',file=__FILE__,line=__LINE__)
            endif

            gflag = global_maxval(iflag, MPI_COMM_ICE)
            if (gflag /= 0) then
               if (my_task == master_task) write(6,*) '   **** ERROR test ',n
               errorflag1(n) = failflag
               errorflag0 = failflag
            endif
         enddo  ! n
      enddo  ! k

      call flush_fileunit(6)
      call ice_barrier()

      ! ---------------------------

      if (my_task == master_task) then
         write(6,*) ' '
         do k = 1,ntests1
            write(6,*) errorflag1(k),stringflag1(k)
         enddo
         write(6,*) ' '
         write(6,*) 'BCSTCHK COMPLETED SUCCESSFULLY'
         if (errorflag0 == passflag) then
            write(6,*) 'BCSTCHK TEST COMPLETED SUCCESSFULLY'
         else
            write(6,*) 'BCSTCHK TEST FAILED'
         endif
      endif

      ! ---------------------------
      ! exit gracefully

      call end_run()

      end program bcstchk

!=======================================================================
