
      program gridavgchk

      ! This tests the CICE grid_average_X2Y methods by
      ! using CICE_InitMod (from the standalone model) to read/initialize
      ! a CICE grid/configuration.  Then methods in grid_average_X2Y
      ! are verified using hardwired inputs with known outputs.
      ! There are lots of issues here
      !   areas (T, U, N, E) are not locally conservative, affect X2YF
      !   X2YF is unmasked which can create havoc in U2T type directions
      !   X2YS is masked but there can be no active cells to average (for instance,
      !        single gridcell wide channels U2T where resuilt is zero)
      !   land block elimination can lead to missing data on halo
      ! This test tries to deal with all these things....

      use CICE_InitMod
      use ice_kinds_mod, only: int_kind, dbl_kind
      use ice_blocks, only: block, get_block, nx_block, ny_block, nblocks_tot
      use ice_constants, only: c0, c1, c2, p25, &
          field_loc_center, field_loc_NEcorner, &
          field_loc_Nface, field_loc_Eface, field_type_scalar
      use ice_communicate, only: my_task, master_task, get_num_procs, MPI_COMM_ICE
      use ice_distribution, only: ice_distributionGetBlockID, ice_distributionGet
      use ice_domain_size, only: nx_global, ny_global, &
          block_size_x, block_size_y, max_blocks
      use ice_domain, only: distrb_info, halo_info, landblockelim
      use ice_fileunits, only: bfbflag
      use ice_exit, only: abort_ice, end_run
      use ice_global_reductions, only: global_minval, global_maxval
      use ice_grid, only: grid_average_X2Y,tarea,uarea,narea,earea,tmask,umask,nmask,emask, &
          hm,uvm,epm,npm

      implicit none

      integer(int_kind) :: i, j, n, ib, ie, jb, je, iblock
      integer(int_kind) :: iglob, jglob
      integer(int_kind) :: blockID, numBlocks
      type (block) :: this_block

      real(dbl_kind) ,allocatable :: array1x(:,:,:), array1y(:,:,:)  ! input
      real(dbl_kind) ,allocatable :: arraysx(:,:,:), arraysy(:,:,:)  ! extra input for NE2T, NE2U
      real(dbl_kind) ,allocatable :: array2x(:,:,:), array2y(:,:,:)  ! output
      real(dbl_kind) ,allocatable :: array3x(:,:,:), array3y(:,:,:)  ! error
      real(dbl_kind) ,allocatable :: wght1(:,:,:), mask1(:,:,:), array2z(:,:,:)  ! extra for explicit
      real(dbl_kind) :: amin, amax, fmax, errtol, errx, erry
      real(dbl_kind) :: deltax0, deltay0, deltax, deltay

      integer(int_kind), parameter :: maxtests = 3
      integer(int_kind), parameter :: maxgroups = 4
      integer(int_kind) :: numtests_cnt, numgroups_cnt
      character(len=16) :: numtests_name(maxtests)
      integer(int_kind) :: nbase(maxgroups)
      character(len=16) :: numgroups_name(maxgroups)
      real(dbl_kind) :: errmax(maxgroups,maxtests)
      integer(int_kind) :: npes, ierr, ntask, testcnt, tottest, mtests, cnt, ng
      integer(int_kind) :: errorflag0,gflag
      integer(int_kind), allocatable :: errorflag(:)
      character(len=32), allocatable :: stringflag(:)
      integer(int_kind), parameter ::  &
         passflag = 0, &
         failflag = 1
      character(len=8), allocatable  :: avgname(:)
      logical, allocatable :: dmask(:,:,:,:)
      real(dbl_kind), allocatable :: errtolconst(:),errtolijind(:),errtolarea(:)

      real(dbl_kind), parameter :: fillval = -1.0e36_dbl_kind
      real(dbl_kind), parameter :: testconst = 100._dbl_kind
      character(len=*), parameter :: subname='(gridavgchk)'

      !-----------------------------------------------------------------
      ! Initialize CICE
      !-----------------------------------------------------------------

      call CICE_Initialize
      npes = get_num_procs()

      numtests_name(1) = 'constant'
      numtests_name(2) = 'ijindex'
      numtests_name(3) = 'area'
      numgroups_name(1) = 'X2YA'
      numgroups_name(2) = 'X2YS'
      numgroups_name(3) = 'X2YF'
      numgroups_name(4) = 'NE2YS'
      nbase(1) = 16
      nbase(2) = 16
      nbase(3) = 0
      nbase(4) = 4
      errmax = c0

      if (.not. landblockelim) nbase(3) = nbase(2)  ! no land block elimination, can test F mappings
      mtests = nbase(1) + nbase(2) + nbase(3) + nbase(4)

      allocate(avgname(mtests))
      allocate(errtolconst(mtests))
      allocate(errtolijind(mtests))
      allocate(errtolarea(mtests))
      errtolconst = c0
      errtolijind = c0
      errtolarea  = c0
      tottest = maxtests * mtests
      allocate(errorflag(tottest))
      allocate(stringflag(tottest))
      allocate(dmask(nx_block,ny_block,max_blocks,mtests))

      n = 0
      errtolconst(n+1:n+nbase(1)) = 0.00001_dbl_kind
      errtolijind(n+1:n+nbase(1)) = 0.10_dbl_kind
      errtolarea (n+1:n+nbase(1)) = 0.04_dbl_kind
      if (nx_global > 200 .and. ny_global > 200) then
         errtolijind(n+1:n+nbase(1)) = 0.03_dbl_kind
         errtolarea (n+1:n+nbase(1)) = 0.003_dbl_kind
      endif
      n=n+1; avgname(n) = 'T2TA' ; dmask(:,:,:,n) = tmask(:,:,:)
      n=n+1; avgname(n) = 'T2UA' ; dmask(:,:,:,n) = umask(:,:,:)
      n=n+1; avgname(n) = 'T2NA' ; dmask(:,:,:,n) = nmask(:,:,:)
      n=n+1; avgname(n) = 'T2EA' ; dmask(:,:,:,n) = emask(:,:,:)
      n=n+1; avgname(n) = 'U2TA' ; dmask(:,:,:,n) = tmask(:,:,:)
      n=n+1; avgname(n) = 'U2UA' ; dmask(:,:,:,n) = umask(:,:,:)
      n=n+1; avgname(n) = 'U2NA' ; dmask(:,:,:,n) = nmask(:,:,:)
      n=n+1; avgname(n) = 'U2EA' ; dmask(:,:,:,n) = emask(:,:,:)
      n=n+1; avgname(n) = 'N2TA' ; dmask(:,:,:,n) = tmask(:,:,:)
      n=n+1; avgname(n) = 'N2UA' ; dmask(:,:,:,n) = umask(:,:,:)
      n=n+1; avgname(n) = 'N2NA' ; dmask(:,:,:,n) = nmask(:,:,:)
      n=n+1; avgname(n) = 'N2EA' ; dmask(:,:,:,n) = emask(:,:,:)
      n=n+1; avgname(n) = 'E2TA' ; dmask(:,:,:,n) = tmask(:,:,:)
      n=n+1; avgname(n) = 'E2UA' ; dmask(:,:,:,n) = umask(:,:,:)
      n=n+1; avgname(n) = 'E2NA' ; dmask(:,:,:,n) = nmask(:,:,:)
      n=n+1; avgname(n) = 'E2EA' ; dmask(:,:,:,n) = emask(:,:,:)

      errtolconst(n+1:n+nbase(2)) = 0.00001_dbl_kind
      errtolijind(n+1:n+nbase(2)) = 0.51_dbl_kind
      errtolarea (n+1:n+nbase(2)) = 0.19_dbl_kind
      if (nx_global > 200 .and. ny_global > 200) then
         errtolarea (n+1:n+nbase(2)) = 0.06_dbl_kind
      endif
      n=n+1; avgname(n) = 'T2TS' ; dmask(:,:,:,n) = tmask(:,:,:)
      n=n+1; avgname(n) = 'T2US' ; dmask(:,:,:,n) = umask(:,:,:)
      n=n+1; avgname(n) = 'T2NS' ; dmask(:,:,:,n) = nmask(:,:,:)
      n=n+1; avgname(n) = 'T2ES' ; dmask(:,:,:,n) = emask(:,:,:)
      n=n+1; avgname(n) = 'U2TS' ; dmask(:,:,:,n) = tmask(:,:,:)
      n=n+1; avgname(n) = 'U2US' ; dmask(:,:,:,n) = umask(:,:,:)
      n=n+1; avgname(n) = 'U2NS' ; dmask(:,:,:,n) = nmask(:,:,:)
      n=n+1; avgname(n) = 'U2ES' ; dmask(:,:,:,n) = emask(:,:,:)
      n=n+1; avgname(n) = 'N2TS' ; dmask(:,:,:,n) = tmask(:,:,:)
      n=n+1; avgname(n) = 'N2US' ; dmask(:,:,:,n) = umask(:,:,:)
      n=n+1; avgname(n) = 'N2NS' ; dmask(:,:,:,n) = nmask(:,:,:)
      n=n+1; avgname(n) = 'N2ES' ; dmask(:,:,:,n) = emask(:,:,:)
      n=n+1; avgname(n) = 'E2TS' ; dmask(:,:,:,n) = tmask(:,:,:)
      n=n+1; avgname(n) = 'E2US' ; dmask(:,:,:,n) = umask(:,:,:)
      n=n+1; avgname(n) = 'E2NS' ; dmask(:,:,:,n) = nmask(:,:,:)
      n=n+1; avgname(n) = 'E2ES' ; dmask(:,:,:,n) = emask(:,:,:)

     if (nbase(3) > 0) then
      errtolconst(n+1:n+nbase(3)) = 0.0065_dbl_kind
      errtolijind(n+1:n+nbase(3)) = 0.65_dbl_kind
      errtolarea (n+1:n+nbase(3)) = 0.04_dbl_kind
      if (nx_global > 200 .and. ny_global > 200) then
         errtolijind(n+1:n+nbase(3)) = 0.22_dbl_kind
         errtolarea (n+1:n+nbase(3)) = 0.004_dbl_kind
      endif
      n=n+1; avgname(n) = 'T2TF' ; dmask(:,:,:,n) = tmask(:,:,:)
      n=n+1; avgname(n) = 'T2UF' ; dmask(:,:,:,n) = umask(:,:,:)
      n=n+1; avgname(n) = 'T2NF' ; dmask(:,:,:,n) = nmask(:,:,:)
      n=n+1; avgname(n) = 'T2EF' ; dmask(:,:,:,n) = emask(:,:,:)
      n=n+1; avgname(n) = 'U2TF' ; dmask(:,:,:,n) = tmask(:,:,:)
      n=n+1; avgname(n) = 'U2UF' ; dmask(:,:,:,n) = umask(:,:,:)
      n=n+1; avgname(n) = 'U2NF' ; dmask(:,:,:,n) = nmask(:,:,:)
      n=n+1; avgname(n) = 'U2EF' ; dmask(:,:,:,n) = emask(:,:,:)
      n=n+1; avgname(n) = 'N2TF' ; dmask(:,:,:,n) = tmask(:,:,:)
      n=n+1; avgname(n) = 'N2UF' ; dmask(:,:,:,n) = umask(:,:,:)
      n=n+1; avgname(n) = 'N2NF' ; dmask(:,:,:,n) = nmask(:,:,:)
      n=n+1; avgname(n) = 'N2EF' ; dmask(:,:,:,n) = emask(:,:,:)
      n=n+1; avgname(n) = 'E2TF' ; dmask(:,:,:,n) = tmask(:,:,:)
      n=n+1; avgname(n) = 'E2UF' ; dmask(:,:,:,n) = umask(:,:,:)
      n=n+1; avgname(n) = 'E2NF' ; dmask(:,:,:,n) = nmask(:,:,:)
      n=n+1; avgname(n) = 'E2EF' ; dmask(:,:,:,n) = emask(:,:,:)
     endif

      errtolconst(n+1:n+nbase(4)) = 0.00001_dbl_kind
      errtolijind(n+1:n+nbase(4)) = 0.51_dbl_kind
      errtolarea (n+1:n+nbase(4)) = 0.12_dbl_kind
      if (nx_global > 200 .and. ny_global > 200) then
         errtolijind(n+1:n+nbase(4)) = 0.26_dbl_kind
         errtolarea (n+1:n+nbase(4)) = 0.03_dbl_kind
      endif
      n=n+1; avgname(n) = 'NE2TS'; dmask(:,:,:,n) = tmask(:,:,:)
      n=n+1; avgname(n) = 'EN2TS'; dmask(:,:,:,n) = tmask(:,:,:)
      n=n+1; avgname(n) = 'NE2US'; dmask(:,:,:,n) = umask(:,:,:)
      n=n+1; avgname(n) = 'EN2US'; dmask(:,:,:,n) = umask(:,:,:)

      if (n /= mtests) then
        call abort_ice(subname//' n ne mtests')
      endif

      !-----------------------------------------------------------------
      ! Testing
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         write(6,*) ' '
         write(6,*) '=========================================================='
         write(6,*) ' '
         write(6,*) 'RunningUnitTest GRIDAVGCHK'
         write(6,*) ' '
         write(6,*) ' npes         = ',npes
         write(6,*) ' my_task      = ',my_task
         write(6,*) ' nx_global    = ',nx_global
         write(6,*) ' ny_global    = ',ny_global
         write(6,*) ' block_size_x = ',block_size_x
         write(6,*) ' block_size_y = ',block_size_y
         write(6,*) ' nblocks_tot  = ',nblocks_tot
         write(6,*) ' tottest      = ',tottest
         write(6,*) ' '
      endif

      errorflag0 = passflag
      errorflag  = passflag
      stringflag = ' '

      ! ---------------------------
      ! TEST GRID AVERAGES
      ! ---------------------------

      if (my_task == master_task) write(6,*) ' '

      allocate(array1x(nx_block,ny_block,max_blocks))
      allocate(array1y(nx_block,ny_block,max_blocks))
      allocate(arraysx(nx_block,ny_block,max_blocks))
      allocate(arraysy(nx_block,ny_block,max_blocks))
      allocate(array2x(nx_block,ny_block,max_blocks))
      allocate(array2y(nx_block,ny_block,max_blocks))
      allocate(array3x(nx_block,ny_block,max_blocks))
      allocate(array3y(nx_block,ny_block,max_blocks))
      allocate(wght1  (nx_block,ny_block,max_blocks))
      allocate(mask1  (nx_block,ny_block,max_blocks))
      allocate(array2z(nx_block,ny_block,max_blocks))

      call ice_distributionGet(distrb_info, numLocalBlocks = numBlocks)

      testcnt = 0

      !----------------
      ! Test constant field
      !----------------

      numtests_cnt = 1
      if (my_task == master_task) then
         write(6,*) ''
         write(6,*) 'TEST constant field, test ',numtests_cnt
      endif

      array1x = testconst
      arraysx = testconst

      do n = 1,mtests
         testcnt = testcnt + 1

         cnt = 0
         do ng = 1,maxgroups
            if (n > cnt) numgroups_cnt = ng
            cnt = cnt + nbase(ng)
         enddo

         errtol = errtolconst(n)
         stringflag(testcnt) = trim(avgname(n))//' const'
         if (my_task == master_task) then
            write(6,*) ''
            write(6,110) trim(stringflag(testcnt))//' test ',testcnt,errtol,numtests_cnt,numgroups_cnt
         endif

         array2x = c0
         if (len_trim(avgname(n)) == 4) then
            call grid_average_X2Y(avgname(n)(4:4),array1x,avgname(n)(1:1),array2x,avgname(n)(3:3))
         else   ! len_trim(avgname(n)) == 5
            if (avgname(n)(1:2) == 'NE') then
               call grid_average_X2Y(avgname(n)(5:5),array1x,avgname(n)(1:1),arraysx,avgname(n)(2:2),array2x,avgname(n)(4:4))
            else  ! EN, swap needed
               call grid_average_X2Y(avgname(n)(5:5),arraysx,avgname(n)(1:1),array1x,avgname(n)(2:2),array2x,avgname(n)(4:4))
            endif
         endif

         array3x = c0
         do iblock = 1,numBlocks
            call ice_distributionGetBlockID(distrb_info, iblock, blockID)
            this_block = get_block(blockID, blockID)
            ib = this_block%ilo
            ie = this_block%ihi
            jb = this_block%jlo
            je = this_block%jhi
            do j = jb,je
             jglob = this_block%j_glob(j)
             do i = ib,ie
               iglob = this_block%i_glob(i)
               array3x(i,j,iblock) = (array2x(i,j,iblock) - testconst)/testconst
               ! if array2 is c0, then there are no valid surrounding points and ignore it
               if (array2x(i,j,iblock) == c0) array3x(i,j,iblock) = c0
               errx = abs(array3x(i,j,iblock))
               ! flag points that are active and error numerically
               if (dmask(i,j,iblock,n) .and. errx > errtol .and. array2x(i,j,iblock) /= c0) then
                  errorflag(testcnt) = failflag
                  errorflag0         = failflag
                  write(100+my_task,*) ''
                  write(100+my_task,100) 'error const  '//trim(avgname(n)),my_task,iblock,i,j,iglob,jglob
                  write(100+my_task,101) 'value, error ',array2x(i,j,iblock),errx
               endif
             enddo
            enddo
         enddo
         gflag = global_maxval(errorflag(testcnt), MPI_COMM_ICE)
         if (my_task == master_task .and. gflag == failflag) write(6,*) ' *** FAIL ***'
         amin = global_minval(array1x, distrb_info)
         amax = global_maxval(array1x, distrb_info)
         if (my_task == master_task) write(6,102) 'input  min/max = ',amin,amax
         amin = global_minval(array2x, distrb_info, dmask(:,:,:,n))
         amax = global_maxval(array2x, distrb_info, dmask(:,:,:,n))
         if (my_task == master_task) write(6,102) 'result min/max = ',amin,amax
         amin = global_minval(array3x, distrb_info, dmask(:,:,:,n))
         amax = global_maxval(array3x, distrb_info, dmask(:,:,:,n))
         if (my_task == master_task) write(6,102) 'error  min/max = ',amin,amax
         amax = global_maxval(abs(array3x), distrb_info, dmask(:,:,:,n))
         errmax(numgroups_cnt,numtests_cnt) = max(errmax(numgroups_cnt,numtests_cnt), amax)
      enddo

      !----------------
      ! Test global i, j fields
      ! for NE2T, NE2U, inputs should result in exact calcs
      !----------------

      numtests_cnt = 2
      if (my_task == master_task) then
         write(6,*) ''
         write(6,*) 'TEST global i, j fields, test ',numtests_cnt
      endif

      array1x = -999.
      arraysx = -999.

      do iblock = 1,numBlocks
         call ice_distributionGetBlockID(distrb_info, iblock, blockID)
         this_block = get_block(blockID, blockID)
         ib = this_block%ilo
         ie = this_block%ihi
         jb = this_block%jlo
         je = this_block%jhi
         do j = jb,je
         do i = ib,ie
            array1x(i,j,iblock) = real(this_block%i_glob(i),kind=dbl_kind)
            array1y(i,j,iblock) = real(this_block%j_glob(j),kind=dbl_kind)
         enddo
         enddo
      enddo

      ! Fill in ghost cells with locally appropriate value
      do iblock = 1,numBlocks
         call ice_distributionGetBlockID(distrb_info, iblock, blockID)
         this_block = get_block(blockID, blockID)
         ib = this_block%ilo
         ie = this_block%ihi
         jb = this_block%jlo
         je = this_block%jhi
         ! skip corners
         do i = ib,ie
            array1x(i,jb-1,iblock) = array1x(i,jb,iblock)
            array1y(i,jb-1,iblock) = array1y(i,jb,iblock) - 1.0_dbl_kind
            array1x(i,je+1,iblock) = array1x(i,je,iblock)
            array1y(i,je+1,iblock) = array1y(i,je,iblock) + 1.0_dbl_kind
         enddo
         ! set corners
         do j = 1,ny_block
            array1x(ib-1,j,iblock) = array1x(ib,j,iblock) - 1.0_dbl_kind
            array1y(ib-1,j,iblock) = array1y(ib,j,iblock)
            array1x(ie+1,j,iblock) = array1x(ie,j,iblock) + 1.0_dbl_kind
            array1y(ie+1,j,iblock) = array1y(ie,j,iblock)
         enddo
      enddo

      arraysx = array1x + 0.5_dbl_kind
      arraysy = array1y - 0.5_dbl_kind

      do n = 1,mtests
         testcnt = testcnt + 1

         cnt = 0
         do ng = 1,maxgroups
            if (n > cnt) numgroups_cnt = ng
            cnt = cnt + nbase(ng)
         enddo

         stringflag(testcnt) = trim(avgname(n))//' ijind'
         if (my_task == master_task) then
            write(6,*) ''
            write(6,110) trim(stringflag(testcnt))//' test ',testcnt,errtolijind(n),numtests_cnt,numgroups_cnt
         endif

         deltax0 = 0.0_dbl_kind
         deltay0 = 0.0_dbl_kind
         if (avgname(n)(1:3) == 'T2U' .or. &
             avgname(n)(1:3) == 'T2E' .or. &
             avgname(n)(1:3) == 'N2U' .or. &
             avgname(n)(1:3) == 'N2E' .or. &
             avgname(n)(1:4) == 'NE2U'.or. &
             avgname(n)(1:4) == 'EN2U') then
            deltax0 = 0.5_dbl_kind
         elseif (avgname(n)(1:3) == 'U2T' .or. &
             avgname(n)(1:3) == 'U2N' .or. &
             avgname(n)(1:3) == 'E2T' .or. &
             avgname(n)(1:3) == 'E2N' ) then
            deltax0 = -0.5_dbl_kind
         endif
         if (avgname(n)(1:3) == 'T2U' .or. &
             avgname(n)(1:3) == 'T2N' .or. &
             avgname(n)(1:3) == 'E2U' .or. &
             avgname(n)(1:3) == 'E2N' ) then
            deltay0 = 0.5_dbl_kind
         elseif (avgname(n)(1:3) == 'U2T' .or. &
             avgname(n)(1:3) == 'U2E' .or. &
             avgname(n)(1:3) == 'N2T' .or. &
             avgname(n)(1:3) == 'N2E' .or. &
             avgname(n)(1:4) == 'NE2T'.or. &
             avgname(n)(1:4) == 'EN2T') then
            deltay0 = -0.5_dbl_kind
         endif

         array2x = c0
         array2y = c0
         if (len_trim(avgname(n)) == 4) then
            call grid_average_X2Y(avgname(n)(4:4),array1x,avgname(n)(1:1),array2x,avgname(n)(3:3))
            call grid_average_X2Y(avgname(n)(4:4),array1y,avgname(n)(1:1),array2y,avgname(n)(3:3))
         else   ! len_trim(avgname(n)) == 5
            if (avgname(n)(1:2) == 'NE') then
               call grid_average_X2Y(avgname(n)(5:5),array1x,avgname(n)(1:1),arraysx,avgname(n)(2:2),array2x,avgname(n)(4:4))
               call grid_average_X2Y(avgname(n)(5:5),array1y,avgname(n)(1:1),arraysy,avgname(n)(2:2),array2y,avgname(n)(4:4))
            else  ! EN, swap needed array1 is N, arrays is E
               call grid_average_X2Y(avgname(n)(5:5),arraysx,avgname(n)(1:1),array1x,avgname(n)(2:2),array2x,avgname(n)(4:4))
               call grid_average_X2Y(avgname(n)(5:5),arraysy,avgname(n)(1:1),array1y,avgname(n)(2:2),array2y,avgname(n)(4:4))
            endif
         endif

         array3x = c0
         errtol = errtolijind(n)
         do iblock = 1,numBlocks
            call ice_distributionGetBlockID(distrb_info, iblock, blockID)
            this_block = get_block(blockID, blockID)
            ib = this_block%ilo
            ie = this_block%ihi
            jb = this_block%jlo
            je = this_block%jhi
            do j = jb,je
             jglob = this_block%j_glob(j)
             do i = ib,ie
               iglob = this_block%i_glob(i)
               deltax = deltax0
               deltay = deltay0
               array3x(i,j,iblock) = array2x(i,j,iblock)-array1x(i,j,iblock)-deltax
               ! if array2 is c0, then there are no valid surrounding points and ignore it
               if (array2x(i,j,iblock) == c0) array3x(i,j,iblock) = c0
               errx = abs(array3x(i,j,iblock))
               array3y(i,j,iblock) = array2y(i,j,iblock)-array1y(i,j,iblock)-deltay
               ! if array2 is c0, then there are no valid surrounding points and ignore it
               if (array2y(i,j,iblock) == c0) array3y(i,j,iblock) = c0
               erry = abs(array3y(i,j,iblock))
               ! flag points that are active and error numerically
               if (dmask(i,j,iblock,n) .and. (errx > errtol .or. erry > errtol)) then
                  errorflag(testcnt) = failflag
                  errorflag0         = failflag
                  write(100+my_task,*) ''
                  write(100+my_task,100) 'error ijind '//trim(avgname(n)),my_task,iblock,i,j,iglob,jglob
                  write(100+my_task,101) 'array2x, err',array2x(i,j,iblock),errx
                  write(100+my_task,101) 'array1x j+1 ',array1x(i-1,j+1,iblock),array1x(i,j+1,iblock),array1x(i+1,j+1,iblock)
                  write(100+my_task,101) 'array1x j   ',array1x(i-1,j  ,iblock),array1x(i,j  ,iblock),array1x(i+1,j  ,iblock)
                  write(100+my_task,101) 'array1x j-1 ',array1x(i-1,j-1,iblock),array1x(i,j-1,iblock),array1x(i+1,j-1,iblock)
                  write(100+my_task,101) 'array2y, err',array2y(i,j,iblock),erry
                  write(100+my_task,101) 'array1y j+1 ',array1y(i-1,j+1,iblock),array1y(i,j+1,iblock),array1y(i+1,j+1,iblock)
                  write(100+my_task,101) 'array1y j   ',array1y(i-1,j  ,iblock),array1y(i,j  ,iblock),array1y(i+1,j  ,iblock)
                  write(100+my_task,101) 'array1y j-1 ',array1y(i-1,j-1,iblock),array1y(i,j-1,iblock),array1y(i+1,j-1,iblock)
                  write(100+my_task,101) 'tarea       ',tarea(i,j,iblock)
                  write(100+my_task,101) 'uarea   j+1 ',uarea  (i-1,j+1,iblock),uarea  (i,j+1,iblock),uarea  (i+1,j+1,iblock)
                  write(100+my_task,101) 'uarea   j   ',uarea  (i-1,j  ,iblock),uarea  (i,j  ,iblock),uarea  (i+1,j  ,iblock)
                  write(100+my_task,101) 'uarea   j-1 ',uarea  (i-1,j-1,iblock),uarea  (i,j-1,iblock),uarea  (i+1,j-1,iblock)
                  write(100+my_task,101) 'hm      j+1 ',hm     (i-1,j+1,iblock),hm     (i,j+1,iblock),hm     (i+1,j+1,iblock)
                  write(100+my_task,101) 'hm      j   ',hm     (i-1,j  ,iblock),hm     (i,j  ,iblock),hm     (i+1,j  ,iblock)
                  write(100+my_task,101) 'hm      j-1 ',hm     (i-1,j-1,iblock),hm     (i,j-1,iblock),hm     (i+1,j-1,iblock)
                  write(100+my_task,101) 'uvm     j+1 ',uvm    (i-1,j+1,iblock),uvm    (i,j+1,iblock),uvm    (i+1,j+1,iblock)
                  write(100+my_task,101) 'uvm     j   ',uvm    (i-1,j  ,iblock),uvm    (i,j  ,iblock),uvm    (i+1,j  ,iblock)
                  write(100+my_task,101) 'uvm     j-1 ',uvm    (i-1,j-1,iblock),uvm    (i,j-1,iblock),uvm    (i+1,j-1,iblock)
               endif
             enddo
            enddo
         enddo

         gflag = global_maxval(errorflag(testcnt), MPI_COMM_ICE)
         if (my_task == master_task .and. gflag == failflag) write(6,*) ' *** FAIL ***'
         amin = global_minval(array1x, distrb_info)
         amax = global_maxval(array1x, distrb_info)
         if (my_task == master_task) write(6,102) 'i_glob   min/max = ',amin,amax
         amin = global_minval(array1y, distrb_info)
         amax = global_maxval(array1y, distrb_info)
         if (my_task == master_task) write(6,102) 'j_glob   min/max = ',amin,amax
         amin = global_minval(array2x, distrb_info, dmask(:,:,:,n))
         amax = global_maxval(array2x, distrb_info, dmask(:,:,:,n))
         if (my_task == master_task) write(6,102) 'i result min/max = ',amin,amax
         amin = global_minval(array2y, distrb_info, dmask(:,:,:,n))
         amax = global_maxval(array2y, distrb_info, dmask(:,:,:,n))
         if (my_task == master_task) write(6,102) 'j result min/max = ',amin,amax
         amin = global_minval(array3x, distrb_info, dmask(:,:,:,n))
         amax = global_maxval(array3x, distrb_info, dmask(:,:,:,n))
         if (my_task == master_task) write(6,102) 'i error  min/max = ',amin,amax
         amin = global_minval(array3y, distrb_info, dmask(:,:,:,n))
         amax = global_maxval(array3y, distrb_info, dmask(:,:,:,n))
         if (my_task == master_task) write(6,102) 'j error  min/max = ',amin,amax
         amax = global_maxval(abs(array3x), distrb_info, dmask(:,:,:,n))
         errmax(numgroups_cnt,numtests_cnt) = max(errmax(numgroups_cnt,numtests_cnt), amax)
         amax = global_maxval(abs(array3y), distrb_info, dmask(:,:,:,n))
         errmax(numgroups_cnt,numtests_cnt) = max(errmax(numgroups_cnt,numtests_cnt), amax)

      enddo

      !----------------
      ! Test area fields
      !----------------

      numtests_cnt = 3
      if (my_task == master_task) then
         write(6,*) ''
         write(6,*) 'TEST area fields, test ',numtests_cnt
      endif

      do n = 1,mtests
         testcnt = testcnt + 1

         cnt = 0
         do ng = 1,maxgroups
            if (n > cnt) numgroups_cnt = ng
            cnt = cnt + nbase(ng)
         enddo

         stringflag(testcnt) = trim(avgname(n))//' area'
         if (my_task == master_task) then
            write(6,*) ''
            write(6,110) trim(stringflag(testcnt))//' test ',testcnt,errtolarea(n),numtests_cnt,numgroups_cnt
         endif

         array1x = -999.
         arraysx = -999.
         mask1   = -999.
         wght1   = -999.
         if     (avgname(n)(1:2) == 'T2') then
            array1x = tarea
            wght1 = tarea
            mask1 = hm
         elseif (avgname(n)(1:2) == 'U2') then
            array1x = uarea
            wght1 = uarea
            mask1 = uvm
         elseif (avgname(n)(1:2) == 'E2') then
            array1x = earea
            wght1 = earea
            mask1 = epm
         elseif (avgname(n)(1:2) == 'N2') then
            array1x = narea
            wght1 = narea
            mask1 = npm
         elseif (avgname(n)(1:3) == 'NE2') then
            array1x = narea
            arraysx = earea
         elseif (avgname(n)(1:3) == 'EN2') then
            array1x = earea
            arraysx = narea
         else
            call abort_ice(subname//' avgname not supported 1x = '//trim(avgname(n)))
         endif

         array2y = -999.
         if     (avgname(n)(2:3) == '2T' .or. &
                 avgname(n)(3:4) == '2T') then
            array2y = tarea
         elseif (avgname(n)(2:3) == '2U' .or. &
                 avgname(n)(3:4) == '2U') then
            array2y = uarea
         elseif (avgname(n)(2:3) == '2E') then
            array2y = earea
         elseif (avgname(n)(2:3) == '2N') then
            array2y = narea
         else
            call abort_ice(subname//' avgname not supported 2y = '//trim(avgname(n)))
         endif

         array2x = c0
         if (len_trim(avgname(n)) == 4) then
!            call grid_average_X2Y(trim(avgname(n)),array1x,array2x)
            call grid_average_X2Y(avgname(n)(4:4),array1x,avgname(n)(1:1),array2x,avgname(n)(3:3))
            ! ------
            ! Extra Explicit Calc Test
            ! ------
            if (avgname(n)(2:2) == '2' .and. (avgname(n)(4:4) == 'S' .or. avgname(n)(4:4) == 'A')) then
               stringflag(testcnt) = trim(stringflag(testcnt))//' + explicit'
               if (avgname(n)(4:4) == 'S') then
                  ! test direct mapping compared to S, array1x*wght1*mask1 where wght1=area and mask1=mask
                  call grid_average_X2Y(avgname(n)(4:4),array1x,avgname(n)(1:1),wght1,mask1,array2z,avgname(n)(3:3))
               elseif (avgname(n)(4:4) == 'A') then
                  ! test direct mapping compared to A, array1x*wght1 where wght1=area and mask1=1.0
                  mask1 = c1
                  call grid_average_X2Y(avgname(n)(4:4),array1x,avgname(n)(1:1),wght1,mask1,array2z,avgname(n)(3:3))
               endif
               fmax = global_maxval(abs(array1x), distrb_info)
               amax = global_maxval(abs(array2z-array2x), distrb_info)
!              tcraig, errtol=c0 doesn't work here, diff seems smaller than roundoff? - interesting
!               errtol = c0
               errtol = 1.0e-20_dbl_kind
               if (amax < fmax * errtol) then
                  if (my_task == master_task) write(6,103) 'PASS explicit avg vs implicit avg ',errtol
               else
                  errorflag(testcnt) = failflag
                  errorflag0         = failflag
                  if (my_task == master_task) write(6,103) 'FAIL explicit avg vs implicit avg ',amax,fmax*errtol
                  amin = global_minval(array2x, distrb_info)
                  amax = global_maxval(array2x, distrb_info)
                  if (my_task == master_task) write(6,103) 'output min/max = ',amin,amax
                  amin = global_minval(array2z, distrb_info)
                  amax = global_maxval(array2z, distrb_info)
                  if (my_task == master_task) write(6,103) 'expout min/max = ',amin,amax
               endif
            endif

         else   ! len_trim(avgname(n)) == 5
             ! no swap needed 1x and sx set based on NE or EN
             call grid_average_X2Y(avgname(n)(5:5),array1x,avgname(n)(1:1),arraysx,avgname(n)(2:2),array2x,avgname(n)(4:4))
         endif

         array3x = c1
         array3y = c1

         do iblock = 1,numBlocks
            call ice_distributionGetBlockID(distrb_info, iblock, blockID)
            this_block = get_block(blockID, blockID)
            ib = this_block%ilo
            ie = this_block%ihi
            jb = this_block%jlo
            je = this_block%jhi
            do j = jb,je
            jglob = this_block%j_glob(j)
            do i = ib,ie
               iglob = this_block%i_glob(i)
               array3x(i,j,iblock) = array2x(i,j,iblock)/array2y(i,j,iblock) - c1
               ! if array2 is c0, then there are no valid surrounding points and ignore it
               if (array2x(i,j,iblock) == c0) array3x(i,j,iblock) = c0
               errx = abs(array3x(i,j,iblock))
               ! flag points that are active and error numerically
               if (dmask(i,j,iblock,n) .and. errx > errtolarea(n)) then
                  errorflag(testcnt) = failflag
                  errorflag0         = failflag
                  write(100+my_task,*) ''
                  write(100+my_task,100) 'error area '//trim(avgname(n)),my_task,iblock,i,j,iglob,jglob
                  write(100+my_task,101) 'out,exact,err',array2x(i,j,iblock),array2y(i,j,iblock),array3x(i,j,iblock)
               endif
            enddo
            enddo
         enddo
         gflag = global_maxval(errorflag(testcnt), MPI_COMM_ICE)
         if (my_task == master_task .and. gflag == failflag) write(6,*) ' *** FAIL ***'
         amin = global_minval(array1x, distrb_info)
         amax = global_maxval(array1x, distrb_info)
         if (my_task == master_task) write(6,103) 'input  min/max = ',amin,amax
         amin = global_minval(array2x, distrb_info, dmask(:,:,:,n))
         amax = global_maxval(array2x, distrb_info, dmask(:,:,:,n))
         if (my_task == master_task) write(6,103) 'output min/max = ',amin,amax
         amin = global_minval(array2y, distrb_info, dmask(:,:,:,n))
         amax = global_maxval(array2y, distrb_info, dmask(:,:,:,n))
         if (my_task == master_task) write(6,103) 'exact  min/max = ',amin,amax
         amin = global_minval(array3x, distrb_info, dmask(:,:,:,n))
         amax = global_maxval(array3x, distrb_info, dmask(:,:,:,n))
         if (my_task == master_task) write(6,102) 'error  min/max = ',amin,amax
         amax = global_maxval(abs(array3x), distrb_info, dmask(:,:,:,n))
         errmax(numgroups_cnt,numtests_cnt) = max(errmax(numgroups_cnt,numtests_cnt), amax)
      enddo

100   format(a,10i8)
101   format(a,3g16.7)
102   format(a,3f16.7)
103   format(a,2g16.7,f16.7)
110   format(a,i8,g16.7,6i8)

      if (my_task == master_task) then
         write(6,*) ' '
         write(6,*) 'Max Errors:'
         do i = 1,maxgroups
         do j = 1,maxtests
            write(6,'(2x,a16,2x,a16,2x,f23.16)') trim(numgroups_name(i)),trim(numtests_name(j)),errmax(i,j)
         enddo
         enddo
      endif

      gflag = global_maxval(errorflag0, MPI_COMM_ICE)
      errorflag0 = gflag
      do n = 1,tottest
         gflag = global_maxval(errorflag(n), MPI_COMM_ICE)
         errorflag(n) = gflag
      enddo

      if (my_task == master_task) then
         write(6,*) ' '
         write(6,*) 'GRIDAVGCHK COMPLETED SUCCESSFULLY'
         do n = 1,tottest
            if (errorflag(n) == passflag) then
               write(6,*) 'PASS ',trim(stringflag(n))
            else
               write(6,*) 'FAIL ',trim(stringflag(n))
            endif
         enddo
         if (errorflag0 == passflag) then
            write(6,*) 'GRIDAVGCHK TEST COMPLETED SUCCESSFULLY'
         else
            write(6,*) 'GRIDAVGCHK TEST FAILED'
         endif
         write(6,*) ' '
         write(6,*) '=========================================================='
         write(6,*) ' '
      endif


      !-----------------------------------------------------------------
      ! Gracefully end
      !-----------------------------------------------------------------

      call end_run()

      end program gridavgchk

!=======================================================================
