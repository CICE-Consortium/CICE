
      program gridavgchk

      ! This tests the CICE grid_average_X2Y methods by
      ! using CICE_InitMod (from the standalone model) to read/initialize
      ! a CICE grid/configuration.  Then methods in grid_average_X2Y
      ! are verified using hardwired inputs with known outputs.

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
      use ice_domain, only: distrb_info, halo_info
      use ice_fileunits, only: bfbflag
      use ice_exit, only: abort_ice, end_run
      use ice_global_reductions, only: global_minval, global_maxval
      use ice_grid, only: grid_average_X2Y,tarea,uarea,narea,earea,tmask,umask,nmask,emask

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

      integer(int_kind) :: npes, ierr, ntask, ntest
      integer(int_kind), parameter :: maxtest = 36
      integer(int_kind) :: errorflag0,errorflag(maxtest),gflag
      character(len=32) :: stringflag(maxtest)
      integer(int_kind), parameter ::  &
         passflag = 0, &
         failflag = 1
      integer(int_kind), parameter :: navg = 12
      character(len=8)  :: avgname(navg)
      logical, allocatable :: dmask(:,:,:,:)

      real(dbl_kind), parameter :: fillval = -1.0e36_dbl_kind
      real(dbl_kind), parameter :: testconst = 100._dbl_kind
      real(dbl_kind), parameter :: errtolconst = 0.01_dbl_kind  ! error tolerance relative to const field
      real(dbl_kind), parameter :: errtolijind = 0.65_dbl_kind  ! absolute error tolerance for ij index field
      real(dbl_kind), parameter :: errtolarea  = 0.06_dbl_kind  ! relative error tolerance for area field ratio
      character(len=*), parameter :: subname='(gridavgchk)'

      !-----------------------------------------------------------------
      ! Initialize CICE
      !-----------------------------------------------------------------

      call CICE_Initialize
      npes = get_num_procs()

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

      allocate(dmask(nx_block,ny_block,max_blocks,navg))
      avgname(1)  = 'T2U'; dmask(:,:,:,1)  = umask(:,:,:)
      avgname(2)  = 'T2N'; dmask(:,:,:,2)  = nmask(:,:,:)
      avgname(3)  = 'T2E'; dmask(:,:,:,3)  = emask(:,:,:)
      avgname(4)  = 'U2T'; dmask(:,:,:,4)  = tmask(:,:,:)
      avgname(5)  = 'U2N'; dmask(:,:,:,5)  = nmask(:,:,:)
      avgname(6)  = 'U2E'; dmask(:,:,:,6)  = emask(:,:,:)
      avgname(7)  = 'N2T'; dmask(:,:,:,7)  = tmask(:,:,:)
      avgname(8)  = 'N2U'; dmask(:,:,:,8)  = umask(:,:,:)
      avgname(9)  = 'N2E'; dmask(:,:,:,9)  = emask(:,:,:)
      avgname(10) = 'E2T'; dmask(:,:,:,10) = tmask(:,:,:)
      avgname(11) = 'E2U'; dmask(:,:,:,11) = umask(:,:,:)
      avgname(12) = 'E2N'; dmask(:,:,:,12) = nmask(:,:,:)

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
            write(6,*) trim(stringflag(ntest)),' test ',ntest
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
             errtol = errtolconst * testconst
             do i = ib,ie
               iglob = this_block%i_glob(i)
               array3x(i,j,iblock) = array2x(i,j,iblock) - array1x(i,j,iblock)
               errx = abs(array3x(i,j,iblock))
               if (dmask(i,j,iblock,n) .and. errx > errtol) then
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
!         write(6,*) 'tcx1x ',my_task,iblock,minval(this_block%i_glob(ib:ie)),maxval(this_block%i_glob(ib:ie))
!         write(6,*) 'tcx1y ',my_task,iblock,minval(this_block%j_glob(jb:je)),maxval(this_block%j_glob(jb:je))
         do j = jb,je
         do i = ib,ie
            array1x(i,j,iblock) = real(this_block%i_glob(i),kind=dbl_kind)
            array1y(i,j,iblock) = real(this_block%j_glob(j),kind=dbl_kind)
         enddo
         enddo
      enddo

      call ice_HaloUpdate(array1x, halo_info, field_loc_center, field_type_scalar, fillval)
      call ice_HaloUpdate(array1y, halo_info, field_loc_center, field_type_scalar, fillval)

      do n = 1,navg
         ntest = ntest + 1

         stringflag(ntest) = trim(avgname(n))//' ijind'
         if (my_task == master_task) then
            write(6,*) ''
            write(6,*) trim(stringflag(ntest)),' test ',ntest
         endif

         deltax0 = 0.0_dbl_kind
         deltay0 = 0.0_dbl_kind
         if (avgname(n) == 'T2U' .or. &
             avgname(n) == 'T2E' .or. &
             avgname(n) == 'N2U' .or. &
             avgname(n) == 'N2E') then
            deltax0 = 0.5_dbl_kind
         elseif (avgname(n) == 'U2T' .or. &
             avgname(n) == 'U2N' .or. &
             avgname(n) == 'E2T' .or. &
             avgname(n) == 'E2N') then
            deltax0 = -0.5_dbl_kind
         endif
         if (avgname(n) == 'T2U' .or. &
             avgname(n) == 'T2N' .or. &
             avgname(n) == 'E2U' .or. &
             avgname(n) == 'E2N') then
            deltay0 = 0.5_dbl_kind
         elseif (avgname(n) == 'U2T' .or. &
             avgname(n) == 'U2E' .or. &
             avgname(n) == 'N2T' .or. &
             avgname(n) == 'N2E') then
            deltay0 = -0.5_dbl_kind
         endif

         array2x = c0
         array2y = c0
         call grid_average_X2Y(trim(avgname(n)),array1x,array2x)
         call grid_average_X2Y(trim(avgname(n)),array1y,array2y)

         array3x = c0
         errtol = errtolijind
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
               ! adjust deltax at wraparound
               if (iglob == 1         .and. deltax < -0.01_dbl_kind) deltax = deltax0 + nx_global/c2
               if (iglob == nx_global .and. deltax >  0.01_dbl_kind) deltax = deltax0 - nx_global/c2
               array3x(i,j,iblock) = array2x(i,j,iblock)-array1x(i,j,iblock)-deltax
               errx = abs(array3x(i,j,iblock))
               array3y(i,j,iblock) = array2y(i,j,iblock)-array1y(i,j,iblock)-deltay
               erry = abs(array3y(i,j,iblock))
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
!                  write(100+my_task,101) 'uarea       ',uarea(i,j,iblock)
!                  write(100+my_task,101) 'tarea   j+1 ',tarea  (i-1,j+1,iblock),tarea  (i,j+1,iblock),tarea  (i+1,j+1,iblock)
!                  write(100+my_task,101) 'tarea   j   ',tarea  (i-1,j  ,iblock),tarea  (i,j  ,iblock),tarea  (i+1,j  ,iblock)
!                  write(100+my_task,101) 'tarea   j-1 ',tarea  (i-1,j-1,iblock),tarea  (i,j-1,iblock),tarea  (i+1,j-1,iblock)
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
            write(6,*) trim(stringflag(ntest)),' test ',ntest
         endif

         array1x = tarea   ! input
         array2y = uarea   ! result
         call ice_HaloUpdate(array1x, halo_info, field_loc_center, field_type_scalar, fillval)
         array2x = c0
         call grid_average_X2Y('T2U',array1x,array2x)

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
               errx = abs(array3x(i,j,iblock))
               if (dmask(i,j,iblock,n) .and. errx > errtolarea) then
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
