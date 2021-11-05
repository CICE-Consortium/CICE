
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
      use ice_boundary, only: ice_HaloUpdate
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
          hm,uvm

      implicit none

      integer(int_kind) :: i, j, n, ib, ie, jb, je, iblock
      integer(int_kind) :: iglob, jglob
      integer(int_kind) :: blockID, numBlocks
      type (block) :: this_block  

      real(dbl_kind) ,allocatable :: array1x(:,:,:), array1y(:,:,:)
      real(dbl_kind) ,allocatable :: array2x(:,:,:), array2y(:,:,:)
      real(dbl_kind) ,allocatable :: array3x(:,:,:), array3y(:,:,:)
      real(dbl_kind) :: amin, amax, errtol, errx, erry
      real(dbl_kind) :: deltax0, deltay0, deltax, deltay

      integer(int_kind) :: npes, ierr, ntask, ntest, maxtest, navg
      integer(int_kind), parameter :: maxgroup = 3
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

      navg = 12
      if (.not. landblockelim) navg=24  ! no land block elimination, can test F mappings
      allocate(avgname(navg))
      allocate(errtolconst(navg))
      allocate(errtolijind(navg))
      allocate(errtolarea(navg))
      maxtest = maxgroup * navg
      allocate(errorflag(maxtest))
      allocate(stringflag(maxtest))
      allocate(dmask(nx_block,ny_block,max_blocks,navg))

      errtolconst(1:12) = 0.0001_dbl_kind
      errtolijind(1:12) = 0.51_dbl_kind
      errtolarea (1:12) = 0.75_dbl_kind
      if (nx_global > 200 .and. ny_global > 200) then
         errtolarea (1:12) = 0.20_dbl_kind
      endif
      avgname(1)  = 'T2US'; dmask(:,:,:,1)  = umask(:,:,:)
      avgname(2)  = 'T2NS'; dmask(:,:,:,2)  = nmask(:,:,:)
      avgname(3)  = 'T2ES'; dmask(:,:,:,3)  = emask(:,:,:)
      avgname(4)  = 'U2TS'; dmask(:,:,:,4)  = tmask(:,:,:)
      avgname(5)  = 'U2NS'; dmask(:,:,:,5)  = nmask(:,:,:)
      avgname(6)  = 'U2ES'; dmask(:,:,:,6)  = emask(:,:,:)
      avgname(7)  = 'N2TS'; dmask(:,:,:,7)  = tmask(:,:,:)
      avgname(8)  = 'N2US'; dmask(:,:,:,8)  = umask(:,:,:)
      avgname(9)  = 'N2ES'; dmask(:,:,:,9)  = emask(:,:,:)
      avgname(10) = 'E2TS'; dmask(:,:,:,10) = tmask(:,:,:)
      avgname(11) = 'E2US'; dmask(:,:,:,11) = umask(:,:,:)
      avgname(12) = 'E2NS'; dmask(:,:,:,12) = nmask(:,:,:)
      if (navg > 12) then
      errtolconst(13:24) = 0.008_dbl_kind
      errtolijind(13:24) = 0.65_dbl_kind
      errtolarea (13:24) = 0.55_dbl_kind
      if (nx_global > 200 .and. ny_global > 200) then
         errtolijind(13:24) = 0.25_dbl_kind
         errtolarea (13:24) = 0.15_dbl_kind
      endif
      avgname(13) = 'T2UF'; dmask(:,:,:,13) = umask(:,:,:)
      avgname(14) = 'T2NF'; dmask(:,:,:,14) = nmask(:,:,:)
      avgname(15) = 'T2EF'; dmask(:,:,:,15) = emask(:,:,:)
      avgname(16) = 'U2TF'; dmask(:,:,:,16) = tmask(:,:,:)
      avgname(17) = 'U2NF'; dmask(:,:,:,17) = nmask(:,:,:)
      avgname(18) = 'U2EF'; dmask(:,:,:,18) = emask(:,:,:)
      avgname(19) = 'N2TF'; dmask(:,:,:,19) = tmask(:,:,:)
      avgname(20) = 'N2UF'; dmask(:,:,:,20) = umask(:,:,:)
      avgname(21) = 'N2EF'; dmask(:,:,:,21) = emask(:,:,:)
      avgname(22) = 'E2TF'; dmask(:,:,:,22) = tmask(:,:,:)
      avgname(23) = 'E2UF'; dmask(:,:,:,23) = umask(:,:,:)
      avgname(24) = 'E2NF'; dmask(:,:,:,24) = nmask(:,:,:)
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
         write(6,*) ' maxtest      = ',maxtest
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
      allocate(array2x(nx_block,ny_block,max_blocks))
      allocate(array2y(nx_block,ny_block,max_blocks))
      allocate(array3x(nx_block,ny_block,max_blocks))
      allocate(array3y(nx_block,ny_block,max_blocks))

      call ice_distributionGet(distrb_info, numLocalBlocks = numBlocks)

      ntest = 0

      !----------------
      ! Test constant field
      !----------------

      if (my_task == master_task) then
         write(6,*) ''
         write(6,*) 'TEST constant field'
      endif

      array1x = testconst

      do n = 1,navg
         ntest = ntest + 1

         stringflag(ntest) = trim(avgname(n))//' const'
         if (my_task == master_task) then
            write(6,*) ''
            write(6,*) trim(stringflag(ntest)),' test ',ntest,errtolconst(n)
         endif

         array2x = c0
         call grid_average_X2Y(trim(avgname(n)),array1x,array2x)

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
             errtol = errtolconst(n) * testconst
             do i = ib,ie
               iglob = this_block%i_glob(i)
               array3x(i,j,iblock) = array2x(i,j,iblock) - array1x(i,j,iblock)
               ! if array2 is c0, then there are no valid surrounding points and ignore it
               if (array2x(i,j,iblock) == c0) array3x(i,j,iblock) = c0
               errx = abs(array3x(i,j,iblock))
               ! flag points that are active and error numerically
               if (dmask(i,j,iblock,n) .and. errx > errtol .and. array2x(i,j,iblock) /= c0) then
                  errorflag(ntest) = failflag
                  errorflag0       = failflag
                  write(100+my_task,*) ''
                  write(100+my_task,100) 'error const  '//trim(avgname(n)),my_task,iblock,i,j,iglob,jglob
                  write(100+my_task,101) 'value, error ',array2x(i,j,iblock),errx
               endif
             enddo
            enddo
         enddo
         amin = global_minval(array1x, distrb_info)
         amax = global_maxval(array1x, distrb_info)
         if (my_task == master_task) write(6,102) 'input  min/max = ',amin,amax
         amin = global_minval(array2x, distrb_info, dmask(:,:,:,n))
         amax = global_maxval(array2x, distrb_info, dmask(:,:,:,n))
         if (my_task == master_task) write(6,102) 'result min/max = ',amin,amax
         amin = global_minval(array3x, distrb_info, dmask(:,:,:,n))
         amax = global_maxval(array3x, distrb_info, dmask(:,:,:,n))
         if (my_task == master_task) write(6,102) 'error  min/max = ',amin,amax
      enddo

      !----------------
      ! Test global i, j fields
      !----------------

      if (my_task == master_task) then
         write(6,*) ''
         write(6,*) 'TEST global i, j fields'
      endif

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

      call ice_HaloUpdate(array1x, halo_info, field_loc_center, field_type_scalar, fillval)
      call ice_HaloUpdate(array1y, halo_info, field_loc_center, field_type_scalar, fillval)

      ! Overwrite the i wraparound points to deal with i/j index average on wraparound
      do iblock = 1,numBlocks
         call ice_distributionGetBlockID(distrb_info, iblock, blockID)
         this_block = get_block(blockID, blockID)
         ib = this_block%ilo
         ie = this_block%ihi
         jb = this_block%jlo
         je = this_block%jhi
         do j = 1,ny_block
         do i = ib,ie
            if (this_block%i_glob(i) == 1        ) array1x(i-1,j,iblock) = 0
            if (this_block%i_glob(i) == nx_global) array1x(i+1,j,iblock) = nx_global+1
         enddo
         enddo
      enddo

      do n = 1,navg
         ntest = ntest + 1

         stringflag(ntest) = trim(avgname(n))//' ijind'
         if (my_task == master_task) then
            write(6,*) ''
            write(6,*) trim(stringflag(ntest)),' test ',ntest,errtolijind(n)
         endif

         deltax0 = 0.0_dbl_kind
         deltay0 = 0.0_dbl_kind
         if (avgname(n)(1:3) == 'T2U' .or. &
             avgname(n)(1:3) == 'T2E' .or. &
             avgname(n)(1:3) == 'N2U' .or. &
             avgname(n)(1:3) == 'N2E') then
            deltax0 = 0.5_dbl_kind
         elseif (avgname(n)(1:3) == 'U2T' .or. &
             avgname(n)(1:3) == 'U2N' .or. &
             avgname(n)(1:3) == 'E2T' .or. &
             avgname(n)(1:3) == 'E2N') then
            deltax0 = -0.5_dbl_kind
         endif
         if (avgname(n)(1:3) == 'T2U' .or. &
             avgname(n)(1:3) == 'T2N' .or. &
             avgname(n)(1:3) == 'E2U' .or. &
             avgname(n)(1:3) == 'E2N') then
            deltay0 = 0.5_dbl_kind
         elseif (avgname(n)(1:3) == 'U2T' .or. &
             avgname(n)(1:3) == 'U2E' .or. &
             avgname(n)(1:3) == 'N2T' .or. &
             avgname(n)(1:3) == 'N2E') then
            deltay0 = -0.5_dbl_kind
         endif

         array2x = c0
         array2y = c0
         call grid_average_X2Y(trim(avgname(n)),array1x,array2x)
         call grid_average_X2Y(trim(avgname(n)),array1y,array2y)

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
                  errorflag(ntest) = failflag
                  errorflag0       = failflag
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

      enddo

      !----------------
      ! Test area fields
      !----------------

      if (my_task == master_task) then
         write(6,*) ''
         write(6,*) 'TEST area fields'
      endif

      do n = 1,navg
         ntest = ntest + 1

         stringflag(ntest) = trim(avgname(n))//' area'
         if (my_task == master_task) then
            write(6,*) ''
            write(6,*) trim(stringflag(ntest)),' test ',ntest,errtolarea(n)
         endif

         array1x = tarea   ! input
         array2y = uarea   ! result
         call ice_HaloUpdate(array1x, halo_info, field_loc_center, field_type_scalar, fillval)
         array2x = c0
         call grid_average_X2Y(trim(avgname(n)),array1x,array2x)

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
                  errorflag(ntest) = failflag
                  errorflag0       = failflag
                  write(100+my_task,*) ''
                  write(100+my_task,100) 'error area '//trim(avgname(n)),my_task,iblock,i,j,iglob,jglob
                  write(100+my_task,101) 'out,exact,err',array2x(i,j,iblock),array2y(i,j,iblock),array3x(i,j,iblock)
               endif
            enddo
            enddo
         enddo
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
      enddo

100   format(a,10i8)
101   format(a,3g16.7)
102   format(a,3f16.7)
103   format(a,2g16.7,f16.7)

      gflag = global_maxval(errorflag0, MPI_COMM_ICE)
      errorflag0 = gflag
      do n = 1,maxtest
         gflag = global_maxval(errorflag(n), MPI_COMM_ICE)
         errorflag(n) = gflag
      enddo

      if (my_task == master_task) then
         write(6,*) ' '
         write(6,*) 'GRIDAVGCHK COMPLETED SUCCESSFULLY'
         do n = 1,maxtest
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
