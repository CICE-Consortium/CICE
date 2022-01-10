
      program sumchk

      ! This tests the CICE ice_global_reductions infrastructure by 
      ! using CICE_InitMod (from the standalone model) to read/initialize
      ! a CICE grid/configuration.  Then methods in ice_global_reductions
      ! are verified using hardwired inputs with known outputs.
      ! A grid needs to be initialized because most of the global reduction
      ! infrastructure assumes haloed and distributed arrays are passed
      ! possibly with a tripole seam.  These interfaces are more than just
      ! layers on top of MPI.  They have the CICE grid/decomposition
      ! infrastructure built-in.

      use CICE_InitMod
      use ice_kinds_mod, only: int_kind, dbl_kind, real_kind
      use ice_communicate, only: my_task, master_task, get_num_procs
      use ice_domain_size, only: nx_global, ny_global
      use ice_domain_size, only: block_size_x, block_size_y, max_blocks
      use ice_domain, only: distrb_info
      use ice_blocks, only: block, get_block, nx_block, ny_block, nblocks_tot
      use ice_distribution, only: ice_distributionGetBlockID, ice_distributionGet
      use ice_constants, only: field_loc_center, field_loc_Nface
      use ice_fileunits, only: bfbflag
      use ice_global_reductions
      use ice_exit, only: abort_ice, end_run

      implicit none

      integer(int_kind) :: i, j, k, l, m, n, iblock, ib, ie, jb, je
      integer(int_kind) :: blockID, numBlocks
      type (block) :: this_block  

      real(dbl_kind)   ,allocatable :: arrayA(:,:,:),arrayB(:,:,:),arrayC(:,:,:)
      integer(int_kind),allocatable :: arrayiA(:,:,:),arrayiB(:,:,:)
      real(dbl_kind)   ,allocatable :: array8(:,:,:),array82(:,:,:)
      real(real_kind)  ,allocatable :: array4(:,:,:),array42(:,:,:)
      integer(int_kind),allocatable :: arrayi1(:,:,:),arrayi2(:,:,:)
      real(dbl_kind)   ,allocatable :: mmask8(:,:,:)
      real(real_kind)  ,allocatable :: mmask4(:,:,:)
      integer(int_kind),allocatable :: mmaski(:,:,:)
      logical          ,allocatable :: lmask (:,:,:)
      real(dbl_kind)   ,allocatable :: vec8(:),sum8(:)
      real(dbl_kind)    :: locval, corval, minval, maxval  ! local, correct, min, max values
      real(dbl_kind)    :: locval8, sumval8, minval8, maxval8
      real(real_kind)   :: locval4, sumval4, minval4, maxval4
      integer(int_kind) :: iocval, locvali, sumvali, corvali, minvali, maxvali
      real(dbl_kind)    :: reldig,reldigchk_now
      real(dbl_kind) ,allocatable :: reldigchk(:,:)

      character(len=8)  :: errorflag0
      character(len=32) :: string
      integer(int_kind),parameter :: ntests1 = 19
      character(len=8)  :: errorflag1(ntests1)
      character(len=32) :: stringflag1(ntests1)
      integer(int_kind),parameter :: ntests2 = 6
      character(len=8)  :: errorflag2(ntests2)
      character(len=32) :: stringflag2(ntests2)
      integer(int_kind),parameter :: ntests3 = 3
      character(len=8)  :: errorflag3(ntests3)
      character(len=32) :: stringflag3(ntests3)
      integer(int_kind),parameter :: ntests4 = 1
      character(len=8)  :: errorflag4(ntests4)
      character(len=32) :: stringflag4(ntests4)

      integer(int_kind) :: npes, ierr, ntask


      integer(int_kind), parameter :: mfld_loc = 2
      integer(int_kind), parameter :: field_loc(mfld_loc) = &
         (/ field_loc_center, field_loc_Nface /)
      character(len=16), parameter :: field_loc_string(mfld_loc) = &
         (/ 'field_loc_center', 'field_loc_Nface ' /)

      integer(int_kind), parameter :: nscale = 4
      real(dbl_kind), parameter    :: lscale(nscale) = &
         (/ 1.0_dbl_kind, &
            1.0e8_dbl_kind, &
            1.0e16_dbl_kind, &
            1.0e32_dbl_kind /)

      integer(int_kind), parameter :: nbflags = 6
      character(len=8), parameter  :: bflags(1:nbflags) = &
         (/ 'off     ','lsum8   ','lsum16  ','lsum4   ','ddpdd   ','reprosum' /)
      character(len=*), parameter ::  &
         passflag = 'PASS', &
         failflag = 'FAIL'
      character(len=*), parameter :: subname='(sumchk)'

      !-----------------------------------------------------------------
      ! Initialize CICE
      !-----------------------------------------------------------------

      call CICE_Initialize

      !-----------------------------------------------------------------
      ! Testing
      !-----------------------------------------------------------------

      errorflag0 = passflag
      errorflag1 = passflag
      errorflag2 = passflag
      errorflag3 = passflag
      errorflag4 = passflag
      npes = get_num_procs()

      if (my_task == master_task) then
         write(6,*) ' '
         write(6,*) '=========================================================='
         write(6,*) ' '
         write(6,*) 'RunningUnitTest SUMCHK'
         write(6,*) ' '
         write(6,*) ' npes         = ',npes
         write(6,*) ' my_task      = ',my_task
         write(6,*) ' nx_global    = ',nx_global
         write(6,*) ' ny_global    = ',ny_global
         write(6,*) ' block_size_x = ',block_size_x
         write(6,*) ' block_size_y = ',block_size_y
         write(6,*) ' nblocks_tot  = ',nblocks_tot
         write(6,*) ' '
      endif

      ! ---------------------------
      ! TEST GLOBAL SUMS
      ! ---------------------------
      ! test difficult sum
      ! fill array with constant value that sums to corval when 2 gridcells per block are excluded
      ! fill those two gridcells per block with very large and opposite signed values
      ! arrayA should sum to corval, arrayB should sum to corval when mask is applied on 2 gridcells
      ! fill 2 extra gridcells with special values
      ! lscale defines relative size of large values
      ! arrayA has large and opposite values in upper right hand corner of block
      ! arrayB has large and same size values in upper right hand corner to check masks
      ! arrayC has large and opposite values in first two values of block
      ! arrayA should add large values at end of a local sum (bad)
      ! arrayC should add large values first then rest of values (not so bad)

      if (my_task == master_task) write(6,*) ' '

      allocate(arrayA (nx_block,ny_block,max_blocks))
      allocate(arrayB (nx_block,ny_block,max_blocks))
      allocate(arrayC (nx_block,ny_block,max_blocks))
      allocate(arrayiA(nx_block,ny_block,max_blocks))
      allocate(arrayiB(nx_block,ny_block,max_blocks))
      allocate(array4 (nx_block,ny_block,max_blocks))
      allocate(array8 (nx_block,ny_block,max_blocks))
      allocate(array42(nx_block,ny_block,max_blocks))
      allocate(array82(nx_block,ny_block,max_blocks))
      allocate(arrayi1(nx_block,ny_block,max_blocks))
      allocate(arrayi2(nx_block,ny_block,max_blocks))
      allocate(mmask4 (nx_block,ny_block,max_blocks))
      allocate(mmask8 (nx_block,ny_block,max_blocks))
      allocate(mmaski (nx_block,ny_block,max_blocks))
      allocate(lmask  (nx_block,ny_block,max_blocks))

      call ice_distributionGet(distrb_info, numLocalBlocks = numBlocks)

      ! correct results for relative digits check in sum
      allocate(reldigchk(nbflags,nscale))
#ifdef NO_R16
      ! lsum16 will revert to a double precision calc like lsum8
      reldigchk(:,:) = 15.7
      reldigchk(1:3,1) = 14.
      reldigchk(4,1) = 3.9
      reldigchk(1:3,2) = 9.
      reldigchk(4,2) = 1.
      reldigchk(1:3,3) = 1.
      reldigchk(4,3) = 0.
      reldigchk(1:3,4) = 0.
      reldigchk(4,4) = 0.
      reldigchk(5,4) = 15.
      if (nx_global == 360 .and. ny_global == 240) then
         reldigchk(1:3,1) = 13.
         reldigchk(5,4) = 14.
      endif
#else
      reldigchk(:,:) = 15.7
      reldigchk(1:2,1) = 14.
      reldigchk(4,1) = 3.9
      reldigchk(1:2,2) = 9.
      reldigchk(4,2) = 1.
      reldigchk(1:2,3) = 1.
      reldigchk(4,3) = 0.
      reldigchk(1:2,4) = 0.
      reldigchk(3,4) = 3.
      reldigchk(4,4) = 0.
      reldigchk(5,4) = 15.
      if (nx_global == 360 .and. ny_global == 240) then
         reldigchk(1:2,1) = 13.
         reldigchk(5,4) = 14.
      endif
#endif

      ! test list
      n = 1    ; stringflag1(n) = 'dble sum easy'
      n = n + 1; stringflag1(n) = 'dble sum'
      n = n + 1; stringflag1(n) = 'real sum'
      n = n + 1; stringflag1(n) = 'intg sum'
      n = n + 1; stringflag1(n) = 'dble sum + dble mask'
      n = n + 1; stringflag1(n) = 'real sum + real mask'
      n = n + 1; stringflag1(n) = 'intg sum + intg mask'
      n = n + 1; stringflag1(n) = 'dble sum + logical mask'
      n = n + 1; stringflag1(n) = 'real sum + logical mask'
      n = n + 1; stringflag1(n) = 'intg sum + logical mask'
      n = n + 1; stringflag1(n) = 'dble prod sum'
      n = n + 1; stringflag1(n) = 'real prod sum'
      n = n + 1; stringflag1(n) = 'intg prod sum'
      n = n + 1; stringflag1(n) = 'dble prod sum + dble mask'
      n = n + 1; stringflag1(n) = 'real prod sum + real mask'
      n = n + 1; stringflag1(n) = 'intg prod sum + intg mask'
      n = n + 1; stringflag1(n) = 'dble prod sum + logical mask'
      n = n + 1; stringflag1(n) = 'real prod sum + logical mask'
      n = n + 1; stringflag1(n) = 'intg prod sum + logical mask'

      do m = 1, mfld_loc

         ! set corval to something a little interesting (not 1.0 for instance which gives atypical results)
         corval = 4.0_dbl_kind/3.0_dbl_kind
         iocval = 8
         ! tuned for gx3 and tx1 only
         if ((nx_global == 100 .and. ny_global == 116) .or. &
             (nx_global == 360 .and. ny_global == 240)) then
            if (field_loc(m) == field_loc_Nface .and. nx_global == 360 .and. ny_global == 240) then
               ! tx1 tripole face, need to adjust local value to remove half of row at ny_global
               ! or modify corval to account for different sum
               locval = corval / real((nblocks_tot*(block_size_x*block_size_y-2)-nx_global/2),dbl_kind)
               corvali = (nblocks_tot*(block_size_x*block_size_y-2)-nx_global/2)*iocval
            else
               locval = corval / real(nblocks_tot*(block_size_x*block_size_y-2),dbl_kind)
               corvali = nblocks_tot*(block_size_x*block_size_y-2)*iocval
            endif
         else
            call abort_ice(subname//' ERROR not set for this grid ')
         endif

      do l = 1, nscale
         if (my_task == master_task) then
            write(6,*) ' '
            write(6,'(a,i4,a,i4)') 'test: m = ',m,': l = ', l
            write(6,'(a,a    )') 'field_loc = ',trim(field_loc_string(m))
            write(6,'(a,e11.4)') 'lscale    = ',lscale(l)
            write(6,*) 'local array value = ',locval
            write(6,*) 'correct value     = ',corval
            write(6,*) 'correct value int = ',corvali
            write(6,*) ' '
            write(6,'(6x,a,26x,a,10x,a,10x,a)') 'test','bfbflag','sum','digits of precision (max is 16)'
         endif

         arrayA(:,:,:) = locval
         arrayB(:,:,:) = locval
         arrayC(:,:,:) = locval
         lmask(:,:,:) = .true.
         do iblock = 1,numBlocks
            call ice_distributionGetBlockID(distrb_info, iblock, blockID)
            this_block = get_block(blockID, blockID)
            ib = this_block%ilo
            ie = this_block%ihi
            jb = this_block%jlo
            je = this_block%jhi

            lmask(ie,je-1,iblock)   = .false.
            lmask(ie,je-2,iblock) = .false.
            arrayA(ie,je-1,iblock)   = locval * lscale(l)
            arrayA(ie,je-2,iblock) = -arrayA(ie,je-1,iblock)
            arrayB(ie,je-1,iblock)   = locval * lscale(l)
            arrayB(ie,je-2,iblock) =  arrayB(ie,je-1,iblock)
            arrayC(ib,jb,iblock)   = locval * lscale(l)
            arrayC(ib+1,jb,iblock) = -arrayC(ib,jb,iblock)
            arrayiA(:,:,iblock) = iocval
            arrayiB(:,:,iblock) = iocval
            arrayiA(ie,je-1,iblock)   = 13 * iocval
            arrayiA(ie,je-2,iblock) = -arrayiA(ie,je-1,iblock)
         enddo

         do k = 1,ntests1
         do n = 1,nbflags
            bfbflag = bflags(n)
            string = stringflag1(k)
            sumval8 = 888.0e12
            sumvali = 8888888

            if (k == 1) then
               array8(:,:,:) = arrayC(:,:,:)
               sumval8 = global_sum(array8, distrb_info, field_loc(m))
            elseif (k == 2) then
               array8(:,:,:) = arrayA(:,:,:)
               sumval8 = global_sum(array8, distrb_info, field_loc(m))
            elseif (k == 3) then
               array4(:,:,:) = arrayA(:,:,:)
               sumval4 = global_sum(array4, distrb_info, field_loc(m))
               sumval8 = sumval4
            elseif (k == 4) then
               arrayi1 = arrayiA
               sumvali = global_sum(arrayi1, distrb_info, field_loc(m))
            elseif (k == 5) then
               mmask8(:,:,:) = 6.0_dbl_kind
               array8(:,:,:) = arrayA(:,:,:)/mmask8(:,:,:)
               sumval8 = global_sum(array8, distrb_info, field_loc(m), mmask=mmask8)
            elseif (k == 6) then
               mmask4(:,:,:) = 6.0_real_kind
               array4(:,:,:) = arrayA(:,:,:)/mmask4(:,:,:)
               sumval4 = global_sum(array4, distrb_info, field_loc(m), mmask=mmask4)
               sumval8 = sumval4
            elseif (k == 7) then
               mmaski(:,:,:) = 2
               arrayi1(:,:,:) = arrayiA(:,:,:)/mmaski(:,:,:)
               sumvali = global_sum(arrayi1, distrb_info, field_loc(m), mmask=mmaski)
            elseif (k == 8) then
               array8(:,:,:) = arrayB(:,:,:)
               sumval8 = global_sum(array8, distrb_info, field_loc(m), lmask=lmask)
            elseif (k == 9) then
               array4(:,:,:) = arrayB(:,:,:)
               sumval4 = global_sum(array4, distrb_info, field_loc(m), lmask=lmask)
               sumval8 = sumval4
            elseif (k == 10) then
               arrayi1(:,:,:) = arrayiB(:,:,:)
               sumvali = global_sum(arrayi1, distrb_info, field_loc(m), lmask=lmask)
            elseif (k == 11) then
               array82(:,:,:) = 7.0_dbl_kind
               array8(:,:,:) = arrayA(:,:,:)/array82(:,:,:)
               sumval8 = global_sum_prod(array8, array82, distrb_info, field_loc(m))
            elseif (k == 12) then
               array42(:,:,:) = 7.0_real_kind
               array4(:,:,:) = arrayA(:,:,:)/array42(:,:,:)
               sumval4 = global_sum_prod(array4, array42, distrb_info, field_loc(m))
               sumval8 = sumval4
            elseif (k == 13) then
               arrayi2(:,:,:) = 4
               arrayi1(:,:,:) = arrayiA(:,:,:)/arrayi2(:,:,:)
               sumvali = global_sum_prod(arrayi1, arrayi2, distrb_info, field_loc(m))
            elseif (k == 14) then
               array82(:,:,:) = 7.0_dbl_kind
               mmask8(:,:,:) = 6.0_dbl_kind
               array8(:,:,:) = arrayA(:,:,:)/(mmask8(:,:,:)*array82(:,:,:))
               sumval8 = global_sum_prod(array8, array82, distrb_info, field_loc(m), mmask=mmask8)
            elseif (k == 15) then
               array42(:,:,:) = 7.0_real_kind
               mmask4(:,:,:) = 6.0_real_kind
               array4(:,:,:) = arrayA(:,:,:)/(mmask4(:,:,:)*array42(:,:,:))
               sumval4 = global_sum_prod(array4, array42, distrb_info, field_loc(m), mmask=mmask4)
               sumval8 = sumval4
            elseif (k == 16) then
               arrayi2(:,:,:) = 2
               mmaski(:,:,:) = 2
               arrayi1(:,:,:) = arrayiA(:,:,:)/(arrayi2(:,:,:)*mmaski(:,:,:))
               sumvali = global_sum_prod(arrayi1, arrayi2, distrb_info, field_loc(m), mmask=mmaski)
            elseif (k == 17) then
               array82(:,:,:) = 7.0_dbl_kind
               array8(:,:,:) = arrayB(:,:,:)/array82(:,:,:)
               sumval8 = global_sum_prod(array8, array82, distrb_info, field_loc(m), lmask=lmask)
            elseif (k == 18) then
               array42(:,:,:) = 7.0_real_kind
               array4(:,:,:) = arrayB(:,:,:)/array42(:,:,:)
               sumval4 = global_sum_prod(array4, array42, distrb_info, field_loc(m), lmask=lmask)
               sumval8 = sumval4
            elseif (k == 19) then
               arrayi2(:,:,:) = 4
               arrayi1(:,:,:) = arrayiB(:,:,:)/(arrayi2(:,:,:))
               sumvali = global_sum_prod(arrayi1, arrayi2, distrb_info, field_loc(m), lmask=lmask)
            else
               call abort_ice(subname//' illegal k sum',file=__FILE__,line=__LINE__)
            endif

            if (string(1:4) == 'intg') then
               ! integer
               if (my_task == master_task) then
                  write(6,'(1x,a,a10,i12)') string,trim(bfbflag), sumvali
               endif
               if (sumvali /= corvali) then
                  errorflag1(k) = failflag
                  errorflag0 = failflag
                  if (my_task == master_task) then
                     write(6,*) '**** ERROR ',sumvali,corvali
                  endif
               endif
            else
               ! real/dbl
               if (sumval8 == corval) then
                  reldig = 16.0_dbl_kind
               elseif (sumval8 == 0._dbl_kind) then
                  reldig = 0
               else
                  reldig = -log10(abs(corval-sumval8)/corval)
               endif
               if (my_task == master_task) then
                  write(6,'(1x,a,a10,g25.17,f8.2)') string,trim(bfbflag), sumval8, reldig
               endif

               ! (real*4) can't have more than 8 digits of precision
               reldigchk_now = reldigchk(n,l)
               if (string(1:4) == 'real') reldigchk_now = min(reldigchk(n,l),7.0)
               if (reldig < reldigchk_now) then
                  errorflag1(k) = failflag
                  errorflag0 = failflag
                  if (my_task == master_task) then
                     write(6,*) '**** ERROR ',reldig,reldigchk_now
                  endif
               endif
            endif
         enddo  ! n
         enddo  ! k
      enddo ! l
      enddo ! m

      ! ---------------------------
      ! Test Global Min/Max
      ! ---------------------------

      if (my_task == master_task) write(6,*) ' '

      n = 1    ; stringflag2(n) = 'dble min/max'
      n = n + 1; stringflag2(n) = 'real min/max'
      n = n + 1; stringflag2(n) = 'intg min/max'
      n = n + 1; stringflag2(n) = 'dble min/max + logical mask'
      n = n + 1; stringflag2(n) = 'real min/max + logical mask'
      n = n + 1; stringflag2(n) = 'intg min/max + logical mask'

      minval = -17.
      maxval = 37.

      ! fill arrays with large values as default
      array8 = 999.0e10_dbl_kind
      array4 = 999.0e10_real_kind
      arrayi1 = 9999999

      n = 1
      ! fill active part of arrays with values between 0 and 10
      do iblock = 1,numBlocks
         call ice_distributionGetBlockID(distrb_info, iblock, blockID)
         this_block = get_block(blockID, blockID)
         ib = this_block%ilo
         ie = this_block%ihi
         jb = this_block%jlo
         je = this_block%jhi
         do j = jb,je
         do i = ib,ie
            n = n + 1
            array8(i,j,iblock) = real(mod(n,10),dbl_kind)
            array4(i,j,iblock) = real(mod(n,8),real_kind)
            arrayi1(i,j,iblock) = mod(n,9)
         enddo
         enddo
      enddo

      ! fill one gridcell with a min and max value
      ntask = max(npes-1,1)-1
      iblock = max(numBlocks-1,1)
      call ice_distributionGetBlockID(distrb_info, iblock, blockID)
      this_block = get_block(blockID, blockID)
      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi
      i = max(ie-3,ib)
      j = max(je-4,jb)
      if (my_task == ntask) then
         array8(i,j,iblock) = minval
         array4(i,j,iblock) = minval
         arrayi1(i,j,iblock) = minval
      endif

      ntask = min(npes,2)-1
      iblock = min(numBlocks,2)
      call ice_distributionGetBlockID(distrb_info, iblock, blockID)
      this_block = get_block(blockID, blockID)
      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi
      i = min(ib+1,ie)
      j = min(jb+2,je)
      if (my_task == ntask) then
         array8(i,j,iblock) = maxval
         array4(i,j,iblock) = maxval
         arrayi1(i,j,iblock) = maxval
      endif

      do k = 1,ntests2
         string = stringflag2(k)
         minval8 =  888e12
         maxval8 = -888e12
         if (k == 1) then
            minval8 = global_minval(array8, distrb_info)
            maxval8 = global_maxval(array8, distrb_info)
         elseif (k == 2) then
            minval4 = global_minval(array4, distrb_info)
            maxval4 = global_maxval(array4, distrb_info)
            minval8 = minval4
            maxval8 = maxval4
         elseif (k == 3) then
            minvali = global_minval(arrayi1, distrb_info)
            maxvali = global_maxval(arrayi1, distrb_info)
            minval8 = minvali
            maxval8 = maxvali
         elseif (k == 4) then
            minval8 = global_minval(array8, distrb_info, lmask=lmask)
            maxval8 = global_maxval(array8, distrb_info, lmask=lmask)
         elseif (k == 5) then
            minval4 = global_minval(array4, distrb_info, lmask=lmask)
            maxval4 = global_maxval(array4, distrb_info, lmask=lmask)
            minval8 = minval4
            maxval8 = maxval4
         elseif (k == 6) then
            minvali = global_minval(arrayi1, distrb_info, lmask=lmask)
            maxvali = global_maxval(arrayi1, distrb_info, lmask=lmask)
            minval8 = minvali
            maxval8 = maxvali
         else
            call abort_ice(subname//' illegal k minmax',file=__FILE__,line=__LINE__)
         endif

         if (my_task == master_task) then
            write(6,'(1x,a,2g16.8)') string, minval8, maxval8
         endif

         if (minval8 /= minval .or. maxval8 /= maxval) then
           errorflag2(k) = failflag
           errorflag0 = failflag
           if (my_task == master_task) then
              write(6,*) '**** ERROR ', minval8, minval, maxval8, maxval
           endif
         endif
      enddo

      ! ---------------------------
      ! Test Scalar Reductions
      ! ---------------------------

      if (my_task == master_task) write(6,*) ' '

      n = 1    ; stringflag3(n) = 'dble scalar min/max/sum'
      n = n + 1; stringflag3(n) = 'real scalar min/max/sum'
      n = n + 1; stringflag3(n) = 'intg scalar min/max/sum'

      minval = -5.
      maxval =  8.

      locval8 = 1.
      locval4 = 1.
      locvali = 1.

      ! fill one gridcell with a min and max value
      ntask = max(npes-1,1)-1
      if (my_task == ntask) then
         locval8 = minval
         locval4 = minval
         locvali = minval
      endif
      ntask = min(npes,2)-1
      if (my_task == ntask) then
         locval8 = maxval
         locval4 = maxval
         locvali = maxval
      endif

      ! compute correct results
      if (npes == 1) then
         minval = maxval
         corval = maxval
      else
         corval = (npes - 2) * 1.0 + minval + maxval
      endif

      do k = 1,ntests3
         string = stringflag3(k)
         minval8 =  888e12
         maxval8 = -888e12
         sumval8 = -888e12
         if (k == 1) then
            minval8 = global_minval(locval8, distrb_info)
            maxval8 = global_maxval(locval8, distrb_info)
            sumval8 = global_sum   (locval8, distrb_info)
         elseif (k == 2) then
            minval4 = global_minval(locval4, distrb_info)
            maxval4 = global_maxval(locval4, distrb_info)
            sumval4 = global_sum   (locval4, distrb_info)
            minval8 = minval4
            maxval8 = maxval4
            sumval8 = sumval4
         elseif (k == 3) then
            minvali = global_minval(locvali, distrb_info)
            maxvali = global_maxval(locvali, distrb_info)
            sumvali = global_sum   (locvali, distrb_info)
            minval8 = minvali
            maxval8 = maxvali
            sumval8 = sumvali
         else
            call abort_ice(subname//' illegal k scalar',file=__FILE__,line=__LINE__)
         endif

         if (my_task == master_task) then
            write(6,'(1x,a,3g16.8)') string, minval8, maxval8, sumval8
         endif

         if (minval8 /= minval .or. maxval8 /= maxval .or. sumval8 /= corval) then
           errorflag3(k) = failflag
           errorflag0 = failflag
           if (my_task == master_task) then
              write(6,*) '**** ERROR ', minval8, minval, maxval8, maxval, sumval8, corval
           endif
         endif
      enddo

      ! ---------------------------
      ! Test Vector Reductions
      ! ---------------------------

      if (my_task == master_task) write(6,*) ' '

      n = 1    ; stringflag4(n) = 'dble sum vector'
      allocate(vec8(3))
      allocate(sum8(3))

      minval = -5.
      maxval =  8.

      vec8(1) = 1.

      ! fill one gridcell with a min and max value
      ntask = max(npes-1,1)-1
      if (my_task == ntask) then
         vec8(1) = minval
      endif
      ntask = min(npes,2)-1
      if (my_task == ntask) then
         vec8(1) = maxval
      endif
      vec8(2) = 2. * vec8(1)
      vec8(3) = 3. * vec8(1)

      ! compute correct results
      if (npes == 1) then
         minval = maxval
         corval = maxval
      else
         corval = (npes - 2) * 1.0 + minval + maxval
      endif

      do k = 1,ntests4
         string = stringflag4(k)
         sum8 = -888e12
         if (k == 1) then
            sum8 = global_allreduce_sum(vec8, distrb_info)
         else
            call abort_ice(subname//' illegal k vector',file=__FILE__,line=__LINE__)
         endif

         if (my_task == master_task) then
            write(6,'(1x,a,3g16.8)') string, sum8(1),sum8(2),sum8(3)
         endif

         if (sum8(1) /= corval .or. sum8(2) /= 2.*corval .or. sum8(3) /= 3.*corval) then
           errorflag4(k) = failflag
           errorflag0 = failflag
           if (my_task == master_task) then
              write(6,*) '**** ERROR ', sum8(1),sum8(2),sum8(3),corval
           endif
         endif
      enddo

      ! ---------------------------

      if (my_task == master_task) then
         write(6,*) ' '
         do k = 1,ntests1
            write(6,*) errorflag1(k),stringflag1(k)
         enddo
         do k = 1,ntests2
            write(6,*) errorflag2(k),stringflag2(k)
         enddo
         do k = 1,ntests3
            write(6,*) errorflag3(k),stringflag3(k)
         enddo
         do k = 1,ntests4
            write(6,*) errorflag4(k),stringflag4(k)
         enddo
         write(6,*) ' '
         write(6,*) 'SUMCHK COMPLETED SUCCESSFULLY'
         if (errorflag0 == passflag) then
            write(6,*) 'SUMCHK TEST COMPLETED SUCCESSFULLY'
         else
            write(6,*) 'SUMCHK TEST FAILED'
         endif
         write(6,*) ' '
         write(6,*) '=========================================================='
         write(6,*) ' '
      endif


      !-----------------------------------------------------------------
      ! Gracefully end
      !-----------------------------------------------------------------

      call end_run()

      end program sumchk

!=======================================================================
