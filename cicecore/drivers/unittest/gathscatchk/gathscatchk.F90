
      module gathscatchk_data

      use CICE_InitMod
      use ice_kinds_mod, only: int_kind, dbl_kind, real_kind, log_kind
      use ice_blocks, only: block, get_block, nx_block, ny_block, nblocks_tot, nghost, &
          i_global, j_global, nblocks_x, nblocks_y, &
          ew_boundary_type, ns_boundary_type
      use ice_boundary, only: ice_HaloUpdate, ice_HaloUpdate_stress
      use ice_gather_scatter
      use ice_constants, only: c0, c1, c2, p5, spval_dbl, spval, spval_int, &
          field_loc_unknown, field_loc_noupdate, &
          field_loc_center, field_loc_NEcorner, &
          field_loc_Nface, field_loc_Eface, &
          field_type_unknown, field_type_noupdate, &
          field_type_scalar, field_type_vector, field_type_angle
      use ice_communicate, only: my_task, master_task, get_num_procs, MPI_COMM_ICE
      use ice_distribution, only: ice_distributionGetBlockID, ice_distributionGet, &
          ice_distributionGetBlockLoc
      use ice_domain_size, only: nx_global, ny_global, &
          block_size_x, block_size_y, max_blocks
      use ice_domain, only: distrb_info, halo_info, nblocks
      use ice_exit, only: abort_ice, end_run
      use ice_global_reductions, only: global_minval, global_maxval, global_sum

      implicit none

      integer(int_kind), parameter ::  &
         passflag = 0, &
         failflag = 1

      end module gathscatchk_data

!=======================================================================

      program gathscatchk

      ! This tests the CICE halo update methods by
      ! using CICE_InitMod (from the standalone model) to read/initialize
      ! a CICE grid/configuration.

      use gathscatchk_data

      implicit none

      integer(int_kind) :: nl, nt, nf, n, nd
      integer(int_kind) :: i, j, ii, jj, isrc, jsrc
      integer(int_kind) :: ilo, ihi, jlo, jhi
      integer(int_kind) :: nxg, nyg, nxgx, nygx, nxb, nyb
      integer(int_kind) :: iblock, ioffset, joffset
      integer(int_kind) :: blockID, numBlocks, processor
      type (block) :: this_block

      ! temporary fields sent for computation
      real(dbl_kind)   , allocatable :: dg1(:,:), dx1(:,:), da1(:,:,:), da2(:,:,:)
      real(real_kind)  , allocatable :: rg1(:,:), rx1(:,:), ra1(:,:,:), ra2(:,:,:)
      integer(int_kind), allocatable :: ig1(:,:), ix1(:,:), ia1(:,:,:), ia2(:,:,:)
      logical(log_kind), allocatable :: lg1(:,:), lx1(:,:), la1(:,:,:), la2(:,:,:)

      ! ba0 is index based array
      ! ba1 is halo update of that
      ! e is with land block elimination
      real(dbl_kind), allocatable :: bx0 (:,:),ba0 (:,:,:)
      real(dbl_kind), allocatable :: bx1 (:,:),ba1 (:,:,:)
      real(dbl_kind), allocatable :: bx0e(:,:),ba0e(:,:,:)
      real(dbl_kind), allocatable :: bx1e(:,:),ba1e(:,:,:)

      ! store the baseline arrays between tests
      real(dbl_kind), allocatable :: bx1h (:,:),ba1h (:,:,:)
      real(dbl_kind), allocatable :: bx1eh(:,:),ba1eh(:,:,:)

      integer(int_kind), parameter :: maxtests1 = 4  ! n  types of tests
      integer(int_kind), parameter :: maxtests2 = 4  ! nd datatype (real, int, etc)
      integer(int_kind), parameter :: maxtypes = 4   ! nt field_type
      integer(int_kind), parameter :: maxlocs = 5    ! nl field_loc
      integer(int_kind), parameter :: maxfills = 2   ! nf fill
      integer(int_kind), parameter :: nz1 = 0
      integer(int_kind), parameter :: nz2 = 0
      real(dbl_kind)    :: chkval,rsign,efac
      character(len=16) :: locs_name(maxlocs), types_name(maxtypes), fill_name(maxfills), data_name(maxtests2)
      integer(int_kind) :: field_loc(maxlocs), field_type(maxtypes)
      logical(log_kind) :: callfill(maxfills)
      integer(int_kind) :: npes, testcnt, tottest, tpcnt, tfcnt
      integer(int_kind) :: errorflag0, gflag, ptcntsum, failcntsum
      integer(int_kind), allocatable :: errorflag(:)
      integer(int_kind), allocatable :: ptcnt(:), failcnt(:)
      character(len=128), allocatable :: teststring(:)
      character(len=32) :: halofld
      character(len=128) :: header0, header1, header2
      logical :: first_call = .true.

      ! debug points
      logical(log_kind), parameter :: debugpts = .false.
      integer(int_kind) :: icg=0, jcg=240         ! point to check, global index
      integer(int_kind) :: ic=1, jc=40, ibc=8     ! point to check, local index

      real(dbl_kind)   , parameter :: dnothing  = -90001.0_dbl_kind
      real(dbl_kind)   , parameter :: dlbepoint = -54545.0_dbl_kind
      real(dbl_kind)   , parameter :: dlbepthal = -53335.0_dbl_kind
      real(dbl_kind)   , parameter :: dfillval  = -87654.0_dbl_kind
      real(real_kind)  , parameter :: rfillval  = -87654.0_real_kind
      integer(int_kind), parameter :: ifillval  = -87654
      logical(log_kind), parameter :: lfillval  = .false.
      real(dbl_kind)   , parameter :: dfillval2 = -45678.0_dbl_kind
      real(real_kind)  , parameter :: rfillval2 =-45678.0_real_kind
      integer(int_kind), parameter :: ifillval2 =-45678
      logical(log_kind), parameter :: lfillval2 =.false.
      real(dbl_kind)               :: hupd_dfillval = c0
      real(real_kind)              :: hupd_rfillval = 0._real_kind
      integer(int_kind)            :: hupd_ifillval = 0
      real(dbl_kind)               :: gath_dfillval
      real(real_kind)              :: gath_rfillval
      integer(int_kind)            :: gath_ifillval
      real(dbl_kind)   , parameter :: scat_dfillval = c0
      character(len=*) , parameter :: subname='(gathscatchk)'

      !-----------------------------------------------------------------
      ! Initialize CICE
      !-----------------------------------------------------------------

      call CICE_Initialize
      npes = get_num_procs()
      call ice_distributionGet(distrb_info, numLocalBlocks = numBlocks)

      gath_dfillval = spval_dbl
      gath_rfillval = spval
      gath_ifillval = spval_int

      locs_name (:) = 'unknown'
      types_name(:) = 'unknown'
      fill_name (:) = 'unknown'
      data_name (:) = 'unknown'
      field_type(:) = field_type_unknown
      field_loc (:) = field_loc_unknown

      types_name(1) = 'scalar'
      field_type(1) = field_type_scalar
      types_name(2) = 'vector'
      field_type(2) = field_type_vector
      types_name(3) = 'angle'
      field_type(3) = field_type_angle
      types_name(4) = 'none'
      field_type(4) = field_type_noupdate
!      types_name(5) = 'unknown'
!      field_type(5) = field_type_unknown  ! aborts in CICE, as expected

      locs_name (1) = 'center'
      field_loc (1)  = field_loc_center
      locs_name (2) = 'NEcorn'
      field_loc (2)  = field_loc_NEcorner
      locs_name (3) = 'Nface'
      field_loc (3)  = field_loc_Nface
      locs_name (4) = 'Eface'
      field_loc (4)  = field_loc_Eface
      locs_name (5) = 'none'
      field_loc (5)  = field_loc_noupdate
!      locs_name (6) = 'unknown'
!      field_loc (6)  = field_loc_unknown  ! aborts in CICE, as expected

      fill_name (1) = 'fill'
      callfill(1) = .true.
      fill_name (2) = 'nofill'
      callfill(2) = .false.

      data_name (1) = 'R8'
      data_name (2) = 'R4'
      data_name (3) = 'I4'
      data_name (4) = 'L1'

      tottest = (maxtests1+1) * maxtests2 * maxlocs * maxtypes * maxfills
      allocate(errorflag(tottest))
      allocate(teststring(tottest))
      allocate(ptcnt(tottest))
      allocate(failcnt(tottest))
      ptcnt(:) = 0
      failcnt(:) = 0

      !-----------------------------------------------------------------
      ! Testing
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         write(6,*) ' '
         write(6,*) '=========================================================='
         write(6,*) ' '
         write(6,*) 'RunningUnitTest GATHSCATCHK'
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

      teststring(:) = ' '

      !-----------------------------------------------------------------
      ! Test gathscat
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         nxg = nx_global
         nyg = ny_global
         nxgx = nx_global+2*nghost
         nygx = ny_global+2*nghost
      else
         nxg = 1
         nyg = 1
         nxgx = 1
         nygx = 1
      endif
      nxb = nx_block
      nyb = ny_block
      allocate(dg1(nxg,nyg), dx1(nxgx,nygx), da1(nx_block,ny_block,max_blocks), da2(nx_block,ny_block,max_blocks))
      allocate(rg1(nxg,nyg), rx1(nxgx,nygx), ra1(nx_block,ny_block,max_blocks), ra2(nx_block,ny_block,max_blocks))
      allocate(ig1(nxg,nyg), ix1(nxgx,nygx), ia1(nx_block,ny_block,max_blocks), ia2(nx_block,ny_block,max_blocks))
      allocate(lg1(nxg,nyg), lx1(nxgx,nygx), la1(nx_block,ny_block,max_blocks), la2(nx_block,ny_block,max_blocks))

      ! on all tasks, bx0 is index based
      !               bx0e is bx0 with lbe
      !               bx1 is halo updated global array
      !               bx1e is bx1 with lbe
      allocate(bx0 (1-nghost:nx_global+nghost,1-nghost:ny_global+nghost), ba0 (nx_block,ny_block,max_blocks))
      allocate(bx1 (1-nghost:nx_global+nghost,1-nghost:ny_global+nghost), ba1 (nx_block,ny_block,max_blocks))
      allocate(bx0e(1-nghost:nx_global+nghost,1-nghost:ny_global+nghost), ba0e(nx_block,ny_block,max_blocks))
      allocate(bx1e(1-nghost:nx_global+nghost,1-nghost:ny_global+nghost), ba1e(nx_block,ny_block,max_blocks))

      ! these are the original updates arrays, store and reset bx1 and bx1e for each test due to tripole
      allocate(bx1h (1-nghost:nx_global+nghost,1-nghost:ny_global+nghost), ba1h (nx_block,ny_block,max_blocks))
      allocate(bx1eh(1-nghost:nx_global+nghost,1-nghost:ny_global+nghost), ba1eh(nx_block,ny_block,max_blocks))

      ! i_global is 1-nghost:n*_global+nghost
      ! local etxended grids are 1:n*_global+2*nghost

      ! fill active points
      ! fill halo for bx0 from indices
      ! then do basic halo update on bx1
      bx0(:,:) = dnothing
      bx1(:,:) = dnothing
      do j = 1-nghost,ny_global+nghost
      do i = 1-nghost,nx_global+nghost
         bx0(i,j) = real(i*1000+j,dbl_kind)
         if (i >= 1 .and. i <= nx_global .and. j >= 1 .and. j <= ny_global) then
            bx1(i,j) = bx0(i,j)
         endif
         if (debugpts .and. (i==icg.or.i==icg+1) .and. (j==jcg.or.j==jcg+1)) then
            write(100+my_task,*) 'tcxa11 ',i,j,bx0(i,j)
            write(100+my_task,*) 'tcxa12 ',i,j,bx1(i,j)
         endif
      enddo
      enddo

      ! mark eliminated land blocks, distinguish interior and outer halo
      bx0e(:,:) = bx0(:,:)
      bx1e(:,:) = bx1(:,:)
      do blockID = 1,nblocks_tot
         call ice_distributionGetBlockLoc(distrb_info, blockID, processor, iblock)
         if (iblock == 0) then
            this_block = get_block(blockID, blockID)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            if (debugpts) write(100+my_task,'(a,i4,a,i6,a,i6,a,i6,a,i6)') 'Eliminate Block ',blockID, &
               ' iglob=',this_block%i_glob(ilo),':',this_block%i_glob(ihi), &
               ' jglob=',this_block%j_glob(jlo),':',this_block%j_glob(jhi)
            do j = jlo, jhi
            do i = ilo, ihi
               ! mark interior gridpoints
               ii = this_block%i_glob(i)
               jj = this_block%j_glob(j)
               !if (debugpts) write(100+my_task,*) 'Eliminate point',blockID,ii,jj
               bx0e(ii,jj) = dlbepoint
               bx1e(ii,jj) = dlbepoint
               ! mark outer halo
               if (ii == 1) then
                  bx0e(1-nghost:0,jj) = dlbepthal
                  bx1e(1-nghost:0,jj) = dlbepthal
                  if (jj == 1) then
                     bx0e(1-nghost:0,1-nghost:0) = dlbepthal
                     bx1e(1-nghost:0,1-nghost:0) = dlbepthal
                  endif
                  if (jj == ny_global) then
                     bx0e(1-nghost:0,ny_global+1:ny_global+nghost) = dlbepthal
                     bx1e(1-nghost:0,ny_global+1:ny_global+nghost) = dlbepthal
                  endif
               endif
               if (ii == nx_global) then
                  bx0e(nx_global+1:nx_global+nghost,jj) = dlbepthal
                  bx1e(nx_global+1:nx_global+nghost,jj) = dlbepthal
                  if (jj == 1) then
                     bx0e(nx_global+1:nx_global+nghost,1-nghost:0) = dlbepthal
                     bx1e(nx_global+1:nx_global+nghost,1-nghost:0) = dlbepthal
                  endif
                  if (jj == ny_global) then
                     bx0e(nx_global+1:nx_global+nghost,ny_global+1:ny_global+nghost) = dlbepthal
                     bx1e(nx_global+1:nx_global+nghost,ny_global+1:ny_global+nghost) = dlbepthal
                  endif
               endif
               if (jj == 1) then
                  bx0e(ii,1-nghost:0) = dlbepthal
                  bx1e(ii,1-nghost:0) = dlbepthal
               endif
               if (jj == ny_global) then
                  bx0e(ii,ny_global+1:ny_global+nghost) = dlbepthal
                  bx1e(ii,ny_global+1:ny_global+nghost) = dlbepthal
               endif
            enddo
            enddo
         endif
      enddo

      do j = 1-nghost,ny_global+nghost
      do i = 1-nghost,nx_global+nghost
         if (debugpts .and. (i==icg.or.i==icg+1) .and. (j==jcg.or.j==jcg+1)) then
            write(100+my_task,*) 'tcxa21 ',i,j,bx0(i,j),bx0e(i,j)
            write(100+my_task,*) 'tcxa22 ',i,j,bx1(i,j),bx1e(i,j)
         endif
      enddo
      enddo

      !--- bottom edge
      do j = 1-nghost, 0
      do i = 1, nx_global
         if (ns_boundary_type == 'cyclic') then
            bx1(i,j) = bx1(i,j+ny_global)
            if (bx1e(i,j) /= dlbepoint .and. bx1e(i,j) /= dlbepthal) bx1e(i,j) = bx1e(i,j+ny_global)
         elseif (ns_boundary_type == 'zero_gradient') then
            bx1(i,j) = bx1(i,1)
            if (bx1e(i,j) /= dlbepoint .and. bx1e(i,j) /= dlbepthal) bx1e(i,j) = bx1e(i,1)
         elseif (ns_boundary_type == 'linear_extrap') then
            efac = real(2-j,dbl_kind)
            bx1(i,j) = efac*bx1(i,1) - (efac-c1)*bx1(i,2)
            if (bx1e(i,j) /= dlbepoint .and. bx1e(i,j) /= dlbepthal) bx1e(i,j) = efac*bx1e(i,1) - (efac-c1)*bx1e(i,2)
         endif
      enddo
      enddo

      !--- top edge
      do j = ny_global+1, ny_global+nghost
      do i = 1, nx_global
         if (ns_boundary_type == 'cyclic') then
            bx1(i,j) = bx1(i,j-ny_global)
            if (bx1e(i,j) /= dlbepoint .and. bx1e(i,j) /= dlbepthal) bx1e(i,j) = bx1e(i,j-ny_global)
         elseif (ns_boundary_type == 'zero_gradient') then
            bx1(i,j) = bx1(i,ny_global)
            if (bx1e(i,j) /= dlbepoint .and. bx1e(i,j) /= dlbepthal) bx1e(i,j) = bx1e(i,ny_global)
         elseif (ns_boundary_type == 'linear_extrap') then
            efac = real(j-ny_global+1,dbl_kind)
            bx1(i,j) = efac*bx1(i,ny_global) - (efac-c1)*bx1(i,ny_global-1)
            if (bx1e(i,j) /= dlbepoint .and. bx1e(i,j) /= dlbepthal) bx1e(i,j) = efac*bx1e(i,ny_global) - (efac-c1)*bx1e(i,ny_global-1)
         endif
      enddo
      enddo

      do j = 1-nghost, ny_global+nghost
      !--- left edge
      do i = 1-nghost, 0
         if (ew_boundary_type == 'cyclic') then
            bx1(i,j) = bx1(i+nx_global,j)
            if (bx1e(i,j) /= dlbepoint .and. bx1e(i,j) /= dlbepthal) bx1e(i,j) = bx1e(i+nx_global,j)
         elseif (ew_boundary_type == 'zero_gradient') then
            bx1(i,j) = bx1(1,j)
            if (bx1e(i,j) /= dlbepoint .and. bx1e(i,j) /= dlbepthal) bx1e(i,j) = bx1e(1,j)
         elseif (ew_boundary_type == 'linear_extrap') then
            efac = real(2-i,dbl_kind)
            bx1(i,j) = efac*bx1(1,j) - (efac-c1)*bx1(2,j)
            if (bx1e(i,j) /= dlbepoint .and. bx1e(i,j) /= dlbepthal) bx1e(i,j) = efac*bx1e(1,j) - (efac-c1)*bx1e(2,j)
         endif
      enddo
      !--- right edge
      do i = nx_global+1, nx_global+nghost
         if (ew_boundary_type == 'cyclic') then
            bx1(i,j) = bx1(i-nx_global,j)
            if (bx1e(i,j) /= dlbepoint .and. bx1e(i,j) /= dlbepthal) bx1e(i,j) = bx1e(i-nx_global,j)
         elseif (ew_boundary_type == 'zero_gradient') then
            bx1(i,j) = bx1(nx_global,j)
            if (bx1e(i,j) /= dlbepoint .and. bx1e(i,j) /= dlbepthal) bx1e(i,j) = bx1e(nx_global,j)
         elseif (ew_boundary_type == 'linear_extrap') then
            efac = real(i-nx_global+1,dbl_kind)
            bx1(i,j) = efac*bx1(nx_global,j) - (efac-c1)*bx1(nx_global-1,j)
            if (bx1e(i,j) /= dlbepoint .and. bx1e(i,j) /= dlbepthal) bx1e(i,j) = efac*bx1e(nx_global,j) - (efac-c1)*bx1e(nx_global-1,j)
         endif
      enddo
      enddo

      do j = 1-nghost,ny_global+nghost
      do i = 1-nghost,nx_global+nghost
         if (debugpts .and. (i==icg.or.i==icg+1) .and. (j==jcg.or.j==jcg+1)) then
            write(100+my_task,*) 'tcxa31 ',i,j,bx0(i,j),bx0e(i,j)
            write(100+my_task,*) 'tcxa32 ',i,j,bx1(i,j),bx1e(i,j)
         endif
      enddo
      enddo

      ! fill ba0, ba1 with global values from bx0 and bx1
      ba0 (:,:,:) = dnothing
      ba1 (:,:,:) = dnothing
      ba0e(:,:,:) = dnothing
      ba1e(:,:,:) = dnothing
      do iblock = 1,numBlocks
         call ice_distributionGetBlockID(distrb_info, iblock, blockID)
         this_block = get_block(blockID, blockID)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         do j = 1,ny_block
         do i = 1,nx_block
            ii = this_block%i_glob(ilo) + i - nghost - 1
            jj = this_block%j_glob(jlo) + j - nghost - 1
            if (ii >= 1-nghost .and. ii <= nx_global+nghost .and. &
                jj >= 1-nghost .and. jj <= ny_global+nghost) then
               ba0 (i,j,iblock) = bx0 (ii,jj)
               ba1 (i,j,iblock) = bx1 (ii,jj)
               ba0e(i,j,iblock) = bx0e(ii,jj)
               ba1e(i,j,iblock) = bx1e(ii,jj)
            endif
         enddo
         enddo
      enddo

      !---------------------------------------------
      ! Test
      !---------------------------------------------

      rsign = c1   ! set the default here

      ! save original baselines, tripole updates them
      bx1h(:,:) = bx1(:,:)
      bx1eh(:,:) = bx1e(:,:)
      ba1h(:,:,:) = ba1(:,:,:)
      ba1eh(:,:,:) = ba1e(:,:,:)

      header1 = 'test   n     i     j       gathscat      expected       ig    jg'
      header2 = 'test   n     i     j  iblock     gathscat      expected       ig    jg'
 1001 format(a,3i6,2g16.6,2i6)
 1002 format(a,4i6,2g16.6,2i6)
      testcnt = 0
      errorflag(:)  = passflag
      do nl = 1, maxlocs
      do nt = 1, maxtypes
      do nf = 1, maxfills
      do n  = 1, maxtests1
      do nd = 1, maxtests2

         bx1(:,:) = bx1h(:,:)
         bx1e(:,:) = bx1eh(:,:)
         ba1(:,:,:) = ba1h(:,:,:)
         ba1e(:,:,:) = ba1eh(:,:,:)

         if (debugpts .and. nd==1) then
            write(100+my_task,*) 'tcx001 ',icg,jcg,bx1(icg,jcg),bx1(icg,jcg+1)
            write(100+my_task,*) 'tcx002 ',icg+1,jcg,bx1(icg+1,jcg),bx1(icg+1,jcg+1)
         endif

         ! Update tripole except with noupdate
         if (field_loc (nl) /= field_loc_noupdate .and. &
             field_type(nt) /= field_type_noupdate .and. &
            (ns_boundary_type == 'tripole' .or. ns_boundary_type == 'tripoleT')) then

            if (debugpts .and. nd==1) then
               write(100+my_task,*) 'tcx011 ',icg,jcg,bx1(icg,jcg),bx1(icg,jcg+1)
               write(100+my_task,*) 'tcx012 ',icg+1,jcg,bx1(icg+1,jcg),bx1(icg+1,jcg+1)
            endif

            ! flip sign for vector/angle except for logicals (nd = 4)
            rsign = c1
            if ((field_type(nt) == field_type_vector .or. field_type(nt) == field_type_angle) .and. &
                 .not. (nd == 4)) then
               rsign = -c1
            endif

            ioffset = -999
            joffset = -999

            ! center offset
            if (ns_boundary_type == 'tripole') then
               ioffset = 0
               joffset = 0
            else ! tripoleT fold
               ioffset = -1
               joffset = 1
            endif

            ! adjust
            ! joffset == 1 is a redundant j line at j=ny_global
            if (field_loc(nl) == field_loc_Eface .or. field_loc(nl) == field_loc_NEcorner) then
               ioffset = ioffset + 1
            endif
            if (field_loc(nl) == field_loc_Nface .or. field_loc(nl) == field_loc_NEcorner) then
               joffset = joffset + 1
            endif

!            if (debugpts .and. nd==1) write(100+my_task,*) 'tcx015',ioffset,joffset

            ! north of active cells
            do j = ny_global, ny_global+nghost
            do i = 1, nx_global
               isrc = nx_global - i + 1 - ioffset   ! ioffset = 0 for tripole center, ioffset = -1 for tripoleT center
               jsrc = ny_global - (j-ny_global) - joffset + 1  ! joffset = 0 for tripole center, joffset = 1 for tripoleT center
               if (isrc < 1        ) isrc = isrc + nx_global
               if (isrc > nx_global) isrc = isrc - nx_global

               !*** for center and Eface on u-fold, and NE corner and Nface
               !*** on T-fold, do not need to replace
               !*** top row of physical domain, so jsrc should be greater than j

!               if (debugpts .and. nd==1) write(100+my_task,*) 'tcx016',i,j,isrc,jsrc
!               if (debugpts .and. nd==1) write(100+my_task,*) 'tcx017',i,j,bx1(i,j)
               if (jsrc > j) then
                  ! do nothing
               elseif (jsrc == j .and. .not.(index(halofld,'STRESS') > 0)) then
                  ! average, but don't average corner points or point if it's redundant with itself (i == isrc)
                  ! GATH/SCAT does NOT do averaging step, skip this relative to haloUpdate
                  ! bx1(i,j) = 0.5_dbl_kind * (bx0(i,jsrc) + rsign*bx0(isrc,jsrc))
               else
                  ! copy
                  bx1(i,j) = rsign * bx0(isrc,jsrc)
                  if (bx0e(isrc,jsrc) == dlbepoint .or. bx0e(isrc,jsrc) == dlbepthal) then
                     bx1e(i,j) = bx0e(isrc,jsrc)
                  else
                     bx1e(i,j) = rsign * bx0e(isrc,jsrc)
                  endif
               endif
!               if (debugpts .and. nd==1) write(100+my_task,*) 'tcx018',i,j,bx1(i,j)
            enddo
            enddo

            ! Update tripole corners
            do j = 1-nghost, ny_global+nghost
            !--- left edge
            do i = 1-nghost, 0
               if (ew_boundary_type == 'cyclic') then
                  bx1(i,j) = bx1(i+nx_global,j)
                  if (bx1e(i,j) /= dlbepoint .and. bx1e(i,j) /= dlbepthal) bx1e(i,j) = bx1e(i+nx_global,j)
               elseif (ew_boundary_type == 'zero_gradient') then
                  bx1(i,j) = bx1(1,j)
                  if (bx1e(i,j) /= dlbepoint .and. bx1e(i,j) /= dlbepthal) bx1e(i,j) = bx1e(1,j)
               elseif (ew_boundary_type == 'linear_extrap') then
                  efac = real(2-i,dbl_kind)
                  bx1(i,j) = efac*bx1(1,j) - (efac-c1)*bx1(2,j)
                  if (bx1e(i,j) /= dlbepoint .and. bx1e(i,j) /= dlbepthal) bx1e(i,j) = efac*bx1e(1,j) - (efac-c1)*bx1e(2,j)
               endif
            enddo
            !--- right edge
            do i = nx_global+1, nx_global+nghost
               if (ew_boundary_type == 'cyclic') then
                  bx1(i,j) = bx1(i-nx_global,j)
                  if (bx1e(i,j) /= dlbepoint .and. bx1e(i,j) /= dlbepthal) bx1e(i,j) = bx1e(i-nx_global,j)
               elseif (ew_boundary_type == 'zero_gradient') then
                  bx1(i,j) = bx1(nx_global,j)
                  if (bx1e(i,j) /= dlbepoint .and. bx1e(i,j) /= dlbepthal) bx1e(i,j) = bx1e(nx_global,j)
               elseif (ew_boundary_type == 'linear_extrap') then
                  efac = real(i-nx_global+1,dbl_kind)
                  bx1(i,j) = efac*bx1(nx_global,j) - (efac-c1)*bx1(nx_global-1,j)
                  if (bx1e(i,j) /= dlbepoint .and. bx1e(i,j) /= dlbepthal) bx1e(i,j) = efac*bx1e(nx_global,j) - (efac-c1)*bx1e(nx_global-1,j)
               endif
            enddo
            enddo

            if (debugpts .and. nd==1) then
               write(100+my_task,*) 'tcx021 ',icg,jcg,bx1(icg,jcg),bx1(icg,jcg+1)
               write(100+my_task,*) 'tcx022 ',icg+1,jcg,bx1(icg+1,jcg),bx1(icg+1,jcg+1)
            endif

            ! update ba1 and ba1e on the tripole
            do iblock = 1,numBlocks
               call ice_distributionGetBlockID(distrb_info, iblock, blockID)
               this_block = get_block(blockID, blockID)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = 1,ny_block
               do i = 1,nx_block
                  ii = this_block%i_glob(ilo) + i - nghost - 1
                  jj = this_block%j_glob(jlo) + j - nghost - 1
                  if (ii >= 1-nghost .and. ii <= nx_global+nghost .and. &
                      jj >= ny_global .and. jj <= ny_global+nghost) then
                     ba1 (i,j,iblock) = bx1 (ii,jj)
                     ba1e(i,j,iblock) = bx1e(ii,jj)
                  endif
               enddo
               enddo
            enddo

         endif  ! tripole

         write(header0,'(a2,1x,i4,1x,3a8,2a16)') data_name(nd),n,locs_name(nl),types_name(nt),fill_name(nf),trim(ew_boundary_type),trim(ns_boundary_type)
         first_call = .true.
         if (n == 1 .and. nl == 1 .and. nt == 1) then
            testcnt = testcnt + 1
            teststring(testcnt) = 'GATH_EXT '//trim(header0)
            da1(:,:,:) = ba0(:,:,:)
            dx1(:,:) = dnothing

            if (nd == 1) then
               if (my_task == master_task) write(6,*) 'Testing ',nl,nt,nf,n,nd
               if (callfill(nf)) then
                  call gather_global(dx1,da1,master_task,distrb_info,fillValue=dfillval,grid_ext=.true.)
               else
                  call gather_global(dx1,da1,master_task,distrb_info,grid_ext=.true.)
               endif
            elseif (nd == 2) then
               if (my_task == master_task) write(6,*) 'Testing ',nl,nt,nf,n,nd
               ra1 = real(da1,real_kind)
               rx1 = real(dx1,real_kind)
               if (callfill(nf)) then
                  call gather_global(rx1,ra1,master_task,distrb_info,fillValue=rfillval,grid_ext=.true.)
               else
                  call gather_global(rx1,ra1,master_task,distrb_info,grid_ext=.true.)
               endif
               dx1 = real(rx1,dbl_kind)
            elseif (nd == 3) then
               if (my_task == master_task) write(6,*) 'Testing ',nl,nt,nf,n,nd
               ia1 = nint(da1)
               ix1 = nint(dx1)
               if (callfill(nf)) then
                  call gather_global(ix1,ia1,master_task,distrb_info,fillValue=ifillval,grid_ext=.true.)
               else
                  call gather_global(ix1,ia1,master_task,distrb_info,grid_ext=.true.)
               endif
               dx1 = real(ix1,dbl_kind)
            elseif (nd == 4) then
               if (my_task == master_task) write(6,*) 'Testing ',nl,nt,nf,n,nd
               la1 = .false.
               lx1 = .false.
               where (mod(nint(da1),2) == 0) la1 = .true.  ! even numbers are true
               if (callfill(nf)) then
                  call gather_global(lx1,la1,master_task,distrb_info,fillValue=lfillval,grid_ext=.true.)
               else
                  call gather_global(lx1,la1,master_task,distrb_info,grid_ext=.true.)
               endif
               dx1 = c0               ! c0 false
               where (lx1) dx1 = c1   ! c1 true
            endif

            if (my_task == master_task) then
               do j = 1-nghost,ny_global+nghost
               do i = 1-nghost,nx_global+nghost
                  ptcnt(testcnt) = ptcnt(testcnt) + 1
                  chkval = bx0e(i,j)
                  if (debugpts .and. nd==1 .and. (i==icg.or.i==icg+1) .and. (j==jcg.or.j==jcg+1)) then
                     write(100+my_task,*) 'tcx111 ',i,j,bx0e(i,j)
                     write(100+my_task,*) 'tcx112 ',i,j,chkval
                  endif
                  if (chkval == dlbepoint) then    ! lbe interior points get filled
                     if (callfill(nf)) then
                        if (nd == 1) chkval = dfillval
                        if (nd == 2) chkval = rfillval
                        if (nd == 3) chkval = ifillval
                        if (nd == 4) chkval = c1   ! odd is false here
                     else
                        if (nd == 1) chkval = gath_dfillval
                        if (nd == 2) chkval = gath_rfillval
                        if (nd == 3) chkval = gath_ifillval
                        if (nd == 4) chkval = c1   ! odd is false here
                     endif
                  elseif (chkval == dlbepthal) then
                     chkval = dnothing             ! lbe halo points are unchanged
                     if (nd == 4) chkval = c1      ! which is false for logicals
                  endif
                  if (debugpts .and. nd==1 .and. (i==icg.or.i==icg+1) .and. (j==jcg.or.j==jcg+1)) then
                     write(100+my_task,*) 'tcx122 ',i,j,chkval
                  endif
                  if (nd == 2) then
                     chkval = real(real(chkval,real_kind),dbl_kind)
                  elseif (nd == 3) then
                     chkval = real(nint(chkval),dbl_kind)
                  elseif (nd == 4) then
                     chkval = c1 - real(mod(nint(abs(chkval)),2),dbl_kind)  ! even number are true = c1
                  endif
                  if (debugpts .and. nd==1 .and. (i==icg.or.i==icg+1) .and. (j==jcg.or.j==jcg+1)) then
                     write(100+my_task,*) 'tcx132 ',i,j,chkval
                  endif
                  if (dx1(i+nghost,j+nghost) /= chkval) then
                     errorflag(testcnt) = failflag
                     failcnt(testcnt) = failcnt(testcnt) + 1
                     if (first_call) then
                        write(100+my_task,*) ' '
                        write(100+my_task,'(a,i4,1x,2a)') '------ TEST = ',testcnt,trim(teststring(testcnt))
                        write(100+my_task,*) ' '
                        write(100+my_task,'(a)') trim(header1)
                        first_call = .false.
                     endif
                     write(100+my_task,1001) 'e1 ',n,i,j,dx1(i+nghost,j+nghost),chkval,i,j
                  endif
               enddo
               enddo
            endif

         elseif (n == 2 .and. nl == 1 .and. nt == 1) then
            testcnt = testcnt + 1
            teststring(testcnt) = 'GATH_STD '//trim(header0)
            da1(:,:,:) = ba0(:,:,:)
            dg1(:,:) = dnothing

            if (nd == 1) then
               if (my_task == master_task) write(6,*) 'Testing ',nl,nt,nf,n,nd
               if (callfill(nf)) then
                  call gather_global(dg1,da1,master_task,distrb_info,fillValue=dfillval)
               else
                  call gather_global(dg1,da1,master_task,distrb_info)
               endif
            elseif (nd == 2) then
               if (my_task == master_task) write(6,*) 'Testing ',nl,nt,nf,n,nd
               rg1 = real(dg1,real_kind)
               ra1 = real(da1,real_kind)
               if (callfill(nf)) then
                  call gather_global(rg1,ra1,master_task,distrb_info,fillValue=rfillval)
               else
                  call gather_global(rg1,ra1,master_task,distrb_info)
               endif
               dg1 = real(rg1,dbl_kind)
            elseif (nd == 3) then
               if (my_task == master_task) write(6,*) 'Testing ',nl,nt,nf,n,nd
               ig1 = nint(dg1)
               ia1 = nint(da1)
               if (callfill(nf)) then
                  call gather_global(ig1,ia1,master_task,distrb_info,fillValue=ifillval)
               else
                  call gather_global(ig1,ia1,master_task,distrb_info)
               endif
               dg1 = real(ig1,dbl_kind)
            elseif (nd == 4) then
               if (my_task == master_task) write(6,*) 'Testing ',nl,nt,nf,n,nd
               lg1 = .false.
               la1 = .false.
               where (mod(nint(da1),2) == 0) la1 = .true.   ! even numbers are true
               if (callfill(nf)) then
                  call gather_global(lg1,la1,master_task,distrb_info,fillValue=lfillval)
               else
                  call gather_global(lg1,la1,master_task,distrb_info)
               endif
               dg1 = c0               ! c0 false
               where (lg1) dg1 = c1   ! c1 true
            endif

            if (my_task == master_task) then
               do j = 1,ny_global
               do i = 1,nx_global
                  ptcnt(testcnt) = ptcnt(testcnt) + 1
                  chkval = bx0e(i,j)
                  if (chkval == dlbepoint) then
                     if (callfill(nf)) then
                        if (nd == 1) chkval = dfillval
                        if (nd == 2) chkval = rfillval
                        if (nd == 3) chkval = ifillval
                        if (nd == 4) chkval = c1   ! odd is false here
                     else
                        if (nd == 1) chkval = gath_dfillval
                        if (nd == 2) chkval = gath_rfillval
                        if (nd == 3) chkval = gath_ifillval
                        if (nd == 4) chkval = c1   ! odd is false here
                     endif
                  endif
                  if (nd == 2) then
                     chkval = real(real(chkval,real_kind),dbl_kind)
                  elseif (nd == 3) then
                     chkval = real(nint(chkval),dbl_kind)
                  elseif (nd == 4) then
                     chkval = c1 - real(mod(nint(abs(chkval)),2),dbl_kind)  ! even number are true = c1
                  endif
                  if (dg1(i,j) /= chkval) then
                     errorflag(testcnt) = failflag
                     failcnt(testcnt) = failcnt(testcnt) + 1
                     if (first_call) then
                        write(100+my_task,*) ' '
                        write(100+my_task,'(a,i4,1x,2a)') '------ TEST = ',testcnt,trim(teststring(testcnt))
                        write(100+my_task,*) ' '
                        write(100+my_task,'(a)') trim(header1)
                        first_call = .false.
                     endif
                     write(100+my_task,1001) 'e2 ',n,i,j,dg1(i,j),chkval,i,j
                  endif
               enddo
               enddo
            endif

         elseif (n == 3 .and. nl == 1 .and. nt == 1) then
            testcnt = testcnt + 1
            teststring(testcnt) = 'SCAT_EXT '//trim(header0)
            if (my_task == master_task) dx1(1:nxgx,1:nygx) = bx1(1-nghost:nx_global+nghost,1-nghost:ny_global+nghost)
            da1(:,:,:) = spval_dbl

            if (nd == 1) then
               if (my_task == master_task) write(6,*) 'Testing ',nl,nt,nf,n,nd
               if (callfill(nf)) then
                  call scatter_global(da1,dx1,master_task,distrb_info,fillValue=dfillval,grid_ext=.true.)
               else
                  call scatter_global(da1,dx1,master_task,distrb_info,grid_ext=.true.)
               endif
            elseif (nd == 2) then
               if (my_task == master_task) write(6,*) 'Testing ',nl,nt,nf,n,nd
               ra1 = spval
               rx1 = real(dx1,real_kind)
               if (callfill(nf)) then
                  call scatter_global(ra1,rx1,master_task,distrb_info,fillValue=rfillval,grid_ext=.true.)
               else
                  call scatter_global(ra1,rx1,master_task,distrb_info,grid_ext=.true.)
               endif
               da1 = real(ra1,dbl_kind)
            elseif (nd == 3) then
               if (my_task == master_task) write(6,*) 'Testing ',nl,nt,nf,n,nd
               ia1 = spval_int
               ix1 = nint(dx1)
               if (callfill(nf)) then
                  call scatter_global(ia1,ix1,master_task,distrb_info,fillValue=ifillval,grid_ext=.true.)
               else
                  call scatter_global(ia1,ix1,master_task,distrb_info,grid_ext=.true.)
               endif
               da1 = real(ia1,dbl_kind)
            elseif (nd == 4) then
               if (my_task == master_task) write(6,*) 'Testing ',nl,nt,nf,n,nd
               la1 = .false.
               lx1 = .false.
               where (mod(nint(dx1),2) == 0) lx1 = .true.   ! even numbers are true
               if (callfill(nf)) then
                  call scatter_global(la1,lx1,master_task,distrb_info,fillValue=lfillval,grid_ext=.true.)
               else
                  call scatter_global(la1,lx1,master_task,distrb_info,grid_ext=.true.)
               endif
               da1 = c0               ! c0 false
               where (la1) da1 = c1   ! c1 true
            endif

            do iblock = 1,numBlocks
               call ice_distributionGetBlockID(distrb_info, iblock, blockID)
               this_block = get_block(blockID, blockID)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo-nghost,jhi+nghost
               do i = ilo-nghost,ihi+nghost
                  ptcnt(testcnt) = ptcnt(testcnt) + 1
                  chkval = ba1(i,j,iblock)
                  if (nd == 2) then
                     chkval = real(real(chkval,real_kind),dbl_kind)
                  elseif (nd == 3) then
                     chkval = real(nint(chkval),dbl_kind)
                  elseif (nd == 4) then
                     chkval = c1 - real(mod(nint(abs(chkval)),2),dbl_kind)  ! even number are true = c1
                  endif
                  if (da1(i,j,iblock) /= chkval) then
                     errorflag(testcnt) = failflag
                     failcnt(testcnt) = failcnt(testcnt) + 1
                     if (first_call) then
                        write(100+my_task,*) ' '
                        write(100+my_task,'(a,i4,1x,2a)') '------ TEST = ',testcnt,trim(teststring(testcnt))
                        write(100+my_task,*) ' '
                        write(100+my_task,'(a)')  trim(header2)
                        first_call = .false.
                     endif
                     write(100+my_task,1002) 'e3 ',n,i,j,iblock,da1(i,j,iblock),chkval,this_block%i_glob(i),this_block%j_glob(j)
                  endif
               enddo
               enddo
            enddo

         ! keep these two tests together, gathhalo test depends on da1 from the scat_std test
         elseif (n == 4) then
            testcnt = testcnt + 1
            teststring(testcnt) = 'SCAT_STD '//trim(header0)
            if (my_task == master_task) then
               dg1(1:nx_global,1:ny_global) = bx0(1:nx_global,1:ny_global)
            endif
            da1(:,:,:) = dnothing

            if (nd == 1) then
               if (my_task == master_task) write(6,*) 'Testing ',nl,nt,nf,n,nd
               if (callfill(nf)) then
                  call scatter_global(da1,dg1,master_task,distrb_info,field_loc(nl),field_type(nt),fillValue=dfillval)
               else
                  call scatter_global(da1,dg1,master_task,distrb_info,field_loc(nl),field_type(nt))
               endif
            elseif (nd == 2) then
               if (my_task == master_task) write(6,*) 'Testing ',nl,nt,nf,n,nd
               ra1 = real(da1,real_kind)
               rg1 = real(dg1,real_kind)
               if (callfill(nf)) then
                  call scatter_global(ra1,rg1,master_task,distrb_info,field_loc(nl),field_type(nt),fillValue=rfillval)
               else
                  call scatter_global(ra1,rg1,master_task,distrb_info,field_loc(nl),field_type(nt))
               endif
               da1 = real(ra1,dbl_kind)
            elseif (nd == 3) then
               if (my_task == master_task) write(6,*) 'Testing ',nl,nt,nf,n,nd
               ia1 = nint(da1)
               ig1 = nint(dg1)
               if (callfill(nf)) then
                  call scatter_global(ia1,ig1,master_task,distrb_info,field_loc(nl),field_type(nt),fillValue=ifillval)
               else
                  call scatter_global(ia1,ig1,master_task,distrb_info,field_loc(nl),field_type(nt))
               endif
               da1 = real(ia1,dbl_kind)
            elseif (nd == 4) then
               if (my_task == master_task) write(6,*) 'Testing ',nl,nt,nf,n,nd
               la1 = .false.
               lg1 = .false.
               where (mod(nint(dg1),2) == 0) lg1 = .true.    ! even numbers are true
               if (callfill(nf)) then
                  call scatter_global(la1,lg1,master_task,distrb_info,field_loc(nl),field_type(nt),fillValue=lfillval)
               else
                  call scatter_global(la1,lg1,master_task,distrb_info,field_loc(nl),field_type(nt))
               endif
               da1 = c0               ! c0 false
               where (la1) da1 = c1   ! c1 true
            endif

            do iblock = 1,numBlocks
               call ice_distributionGetBlockID(distrb_info, iblock, blockID)
               this_block = get_block(blockID, blockID)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo-nghost,jhi+nghost
               do i = ilo-nghost,ihi+nghost
                  ptcnt(testcnt) = ptcnt(testcnt) + 1
                  chkval = ba1(i,j,iblock)
                  if (debugpts .and. nd==1 .and. (i==ic.or.i==ic+1) .and. (j==jc.or.j==jc+1) .and. iblock==ibc) then
                     write(100+my_task,*) 'tcx441',i,j,iblock,da1(i,j,iblock)
                     write(100+my_task,*) 'tcx442',i,j,iblock,ba1(i,j,iblock)
                     write(100+my_task,*) 'tcx443',i,j,iblock,chkval
                  endif
                  ! if chkval is set to dnothing or it's a no update call, then there is no outside halo update done
                  ! and values returned are the fill values, set chkval to fill value in either case
                  if ((chkval == dnothing) .or. &
                      ((field_loc(nl) == field_loc_noupdate .or. field_type(nt) == field_type_noupdate) .and. &
                       ((this_block%i_glob(ilo) == 1 .and. i < ilo) .or. &
                        (this_block%i_glob(ihi) == nx_global .and. i > ihi) .or. &
                        (this_block%j_glob(jlo) == 1 .and. j < jlo) .or. &
                        (this_block%j_glob(jhi) == ny_global .and. j > jhi)))) then
                     if (nd == 4) then
                        chkval = c1  ! odd = false, false is lfillval and scat_dfillval for logicals
                     elseif (callfill(nf)) then
                        chkval = dfillval
                     else
                        chkval = scat_dfillval
                     endif
                  endif
                  if (nd == 2) then
                     chkval = real(real(chkval,real_kind),dbl_kind)
                  elseif (nd == 3) then
                     chkval = real(nint(chkval),dbl_kind)
                  elseif (nd == 4) then
                     chkval = c1 - real(mod(nint(abs(chkval)),2),dbl_kind)  ! even number are true = c1
                  endif
                  if (debugpts .and. nd==1 .and. (i==ic.or.i==ic+1) .and. (j==jc.or.j==jc+1) .and. iblock==ibc) then
                     write(100+my_task,*) 'tcx446',i,j,iblock,da1(i,j,iblock)
                     write(100+my_task,*) 'tcx447',i,j,iblock,ba1(i,j,iblock)
                     write(100+my_task,*) 'tcx448',i,j,iblock,chkval
                  endif
                  if (da1(i,j,iblock) /= chkval) then
                     errorflag(testcnt) = failflag
                     failcnt(testcnt) = failcnt(testcnt) + 1
                     if (first_call) then
                        write(100+my_task,*) ' '
                        write(100+my_task,'(a,i4,1x,2a)') '------ TEST = ',testcnt,trim(teststring(testcnt))
                        write(100+my_task,*) ' '
                        write(100+my_task,'(a)')  trim(header2)
                        first_call = .false.
                     endif
                     write(100+my_task,1002) 'e4 ',n,i,j,iblock,da1(i,j,iblock),chkval,this_block%i_glob(i),this_block%j_glob(j)
                  endif
               enddo
               enddo
            enddo

            testcnt = testcnt + 1
            teststring(testcnt) = 'GATHHALO '//trim(header0)
            ! these ifs are strictly correct, but aren't really needed, do they slow the code down?
            if (nd == 1) da2(:,:,:) = dfillval2
            if (nd == 2) ra2(:,:,:) = rfillval2
            if (nd == 3) ia2(:,:,:) = ifillval2
            if (nd == 4) la2(:,:,:) = .false.
            do iblock = 1,numBlocks
               call ice_distributionGetBlockID(distrb_info, iblock, blockID)
               this_block = get_block(blockID, blockID)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo,jhi
               do i = ilo,ihi
                  if (nd == 1) da2(i,j,iblock) = da1(i,j,iblock)
                  if (nd == 2) ra2(i,j,iblock) = ra1(i,j,iblock)
                  if (nd == 3) ia2(i,j,iblock) = ia1(i,j,iblock)
                  if (nd == 4) la2(i,j,iblock) = la1(i,j,iblock)
               enddo
               enddo
            enddo

            if (nd == 1) then
               if (callfill(nf)) then
                  call ice_haloUpdate(da2, halo_info, field_loc(nl), field_type(nt), fillvalue=dfillval)
               else
                  call ice_haloUpdate(da2, halo_info, field_loc(nl), field_type(nt))
               endif
            elseif (nd == 2) then
               if (callfill(nf)) then
                  call ice_haloUpdate(ra2, halo_info, field_loc(nl), field_type(nt), fillvalue=rfillval)
               else
                  call ice_haloUpdate(ra2, halo_info, field_loc(nl), field_type(nt))
               endif
               da2 = real(ra2,dbl_kind)
            elseif (nd == 3) then
               if (callfill(nf)) then
                  call ice_haloUpdate(ia2, halo_info, field_loc(nl), field_type(nt), fillvalue=ifillval)
               else
                  call ice_haloUpdate(ia2, halo_info, field_loc(nl), field_type(nt))
               endif
               da2 = real(ia2,dbl_kind)
            elseif (nd == 4) then
               if (callfill(nf)) then
                  call ice_haloUpdate(la2, halo_info, field_loc(nl), field_type(nt), fillvalue=lfillval)
               else
                  call ice_haloUpdate(la2, halo_info, field_loc(nl), field_type(nt))
               endif
               da2 = c0               ! c0 false
               where (la2) da2 = c1   ! c1 true
            endif

!            if (debugpts .and. nd==1) then
!               write(100+my_task,*) 'tcx521',da1(1-nghost:nx_global+nghost,ny_global,1)
!               write(100+my_task,*) 'tcx522',da1(1-nghost:nx_global+nghost,ny_global+1,1)
!               write(100+my_task,*) 'tcx531',da2(1-nghost:nx_global+nghost,ny_global,1)
!               write(100+my_task,*) 'tcx532',da2(1-nghost:nx_global+nghost,ny_global+1,1)
!            endif

            do iblock = 1,numBlocks
               call ice_distributionGetBlockID(distrb_info, iblock, blockID)
               this_block = get_block(blockID, blockID)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo-nghost,jhi+nghost
               do i = ilo-nghost,ihi+nghost
                  ii = this_block%i_glob(ilo) + i - nghost - 1
                  jj = this_block%j_glob(jlo) + j - nghost - 1
                  ptcnt(testcnt) = ptcnt(testcnt) + 1
                  if (nd == 1) chkval = da1(i,j,iblock)
                  if (nd == 2) chkval = real(ra1(i,j,iblock),dbl_kind)
                  if (nd == 3) chkval = real(ia1(i,j,iblock),dbl_kind)
                  if (nd == 4) then
                     chkval = c0
                     if (la1(i,j,iblock)) chkval = c1
                  endif

                  if (debugpts .and. nd==1 .and. (i==ic.or.i==ic+1) .and. (j==jc.or.j==jc+1) .and. iblock==ibc) then
                     write(100+my_task,*) 'tcx531',i,j,iblock,da1(i,j,iblock)
                     write(100+my_task,*) 'tcx532',i,j,iblock,da2(i,j,iblock)
                     write(100+my_task,*) 'tcx533',i,j,iblock,chkval
                  endif

                  ! no fill for open or closed leaves halo original values, find points with ba1e=donothing
                  if (.not.callfill(nf) .and. ba1e(i,j,iblock) == dnothing) then
                     if (nd == 1) chkval = dfillval2
                     if (nd == 2) chkval = rfillval2
                     if (nd == 3) chkval = ifillval2
                     if (nd == 4) chkval = c0  ! c0 = false
                  endif

                  ! no fill for tripole on southern boundary is filled with c0
                  if (.not.callfill(nf) .and. j < jlo .and. this_block%j_glob(jlo) == 1 .and. &
                      (ns_boundary_type == 'tripole' .or. ns_boundary_type == 'tripoleT')) then
                     chkval = c0
                  endif

                  if (debugpts .and. nd==1 .and. (i==ic.or.i==ic+1) .and. (j==jc.or.j==jc+1) .and. iblock==ibc) then
                     write(100+my_task,*) 'tcx541',i,j,iblock,da1(i,j,iblock)
                     write(100+my_task,*) 'tcx542',i,j,iblock,da2(i,j,iblock)
                     write(100+my_task,*) 'tcx543',i,j,iblock,chkval
                  endif

                  ! eliminated blocks get fill value
                  if (ba1e(i,j,iblock) == dlbepoint .or. ba1e(i,j,iblock) == dlbepthal) then
                     if (callfill(nf)) then
                        if (nd == 1) chkval = dfillval
                        if (nd == 2) chkval = rfillval
                        if (nd == 3) chkval = ifillval
                        if (nd == 4) chkval = c0  ! c0 = false
                     else
                        if (nd == 1) chkval = hupd_dfillval
                        if (nd == 2) chkval = hupd_rfillval
                        if (nd == 3) chkval = hupd_ifillval
                        if (nd == 4) chkval = c0  ! c0 = false
                     endif
                     ! this is really strange, the halo updates the corner point of an eliminated block
                     ! on the outer edge of the domain. it probably shouldn't.  with tripole and tripoleT
                     ! it does it on the reflected line too so add an extra check for tripole at j=ny_global.
                     ! force chkval to match da2 for these points.
                     ! again, this is just if an eliminated block is on the outer edge
                     if (((ii < 1 .or. ii > nx_global) .and. (j < jlo .or. j > jhi)) .or. &
                         ((jj < 1 .or. jj > nx_global) .and. (i < ilo .or. i > ihi))) then
                        chkval = da2(i,j,iblock)
                     endif
                     if ((ns_boundary_type == 'tripole' .or. ns_boundary_type == 'tripoleT') .and. &
                         ((jj == ny_global) .and. (i < ilo .or. i > ihi))) then
                        chkval = da2(i,j,iblock)
                     endif
                     ! this is wrong, the halo update is switching the sign of the fill values
                     ! when reflecting across the tripole if rsign = -1., need to adjust chkval
                     ! accordingly but ultimately, need to fix the haloupdate.  to limit this fix,
                     ! only apply it when rsign=-1 (could apply it when rsign=1 but da2 is never -chkval).
                     ! we are setting chkval to the incorrect value to match the haloupdate.
                     if ((ns_boundary_type == 'tripole' .or. ns_boundary_type == 'tripoleT') .and. &
                         rsign == -c1 .and. da2(i,j,iblock) == -chkval) chkval = -chkval
                  endif

                  if (debugpts .and. nd==1 .and. (i==ic.or.i==ic+1) .and. (j==jc.or.j==jc+1) .and. iblock==ibc) then
                     write(100+my_task,*) 'tcx551',i,j,iblock,da1(i,j,iblock)
                     write(100+my_task,*) 'tcx552',i,j,iblock,da2(i,j,iblock)
                     write(100+my_task,*) 'tcx553',i,j,iblock,chkval
                     write(100+my_task,*) 'tcx554',i,j,iblock,ba0(i,j,iblock)
                     write(100+my_task,*) 'tcx555',i,j,iblock,ba0e(i,j,iblock)
                     write(100+my_task,*) 'tcx556',i,j,iblock,ba1(i,j,iblock)
                     write(100+my_task,*) 'tcx557',i,j,iblock,ba1e(i,j,iblock)
                  endif

                  ! halo update and scatter do NOT match on the redundant line, so ignore check and force a pass there
                  if ((this_block%j_glob(j) == ny_global .and. ns_boundary_type == 'tripole' .and. &
                      (field_loc(nl) == field_loc_NEcorner .or. field_loc(nl) == field_loc_Nface)) .or. &
                      (this_block%j_glob(j) == ny_global .and. ns_boundary_type == 'tripoleT' .and. &
                      (field_loc(nl) == field_loc_center .or. field_loc(nl) == field_loc_Eface))) then
                     chkval = da2(i,j,iblock)
                  endif

                  if (debugpts .and. nd==1 .and. (i==ic.or.i==ic+1) .and. (j==jc.or.j==jc+1) .and. iblock==ibc) then
                     write(100+my_task,*) 'tcx561',i,j,iblock,da1(i,j,iblock)
                     write(100+my_task,*) 'tcx562',i,j,iblock,da2(i,j,iblock)
                     write(100+my_task,*) 'tcx563',i,j,iblock,chkval
                  endif

                  ! no update leaves halo original values
                  if ((field_loc(nl) == field_loc_noupdate .or. field_type(nt) == field_type_noupdate) .and. &
                      (i < ilo .or. i > ihi .or. j < jlo .or. j > jhi)) then
                     if (nd == 1) chkval = dfillval2
                     if (nd == 2) chkval = rfillval2
                     if (nd == 3) chkval = ifillval2
                     if (nd == 4) chkval = c0  ! c0 = false
                  endif

                  if (debugpts .and. nd==1 .and. (i==ic.or.i==ic+1) .and. (j==jc.or.j==jc+1) .and. iblock==ibc) then
                     write(100+my_task,*) 'tcx571',i,j,iblock,da1(i,j,iblock)
                     write(100+my_task,*) 'tcx572',i,j,iblock,da2(i,j,iblock)
                     write(100+my_task,*) 'tcx573',i,j,iblock,chkval
                  endif

                  if (da2(i,j,iblock) /= chkval) then
                     errorflag(testcnt) = failflag
                     failcnt(testcnt) = failcnt(testcnt) + 1
                     if (first_call) then
                        write(100+my_task,*) ' '
                        write(100+my_task,'(a,i4,1x,2a)') '------ TEST = ',testcnt,trim(teststring(testcnt))
                        write(100+my_task,*) ' '
                        write(100+my_task,'(a)')  trim(header2)
                        first_call = .false.
                     endif
                     write(100+my_task,1002) 'e5 ',n,i,j,iblock,da2(i,j,iblock),chkval,this_block%i_glob(i),this_block%j_glob(j)
                  endif
               enddo
               enddo
            enddo

         endif

      enddo
      enddo
      enddo
      enddo
      enddo

      ! ---------------------------
      ! SUMMARY
      ! ---------------------------

      errorflag0 = passflag
      do n = 1,testcnt
         gflag = global_maxval(errorflag(n), MPI_COMM_ICE)
         errorflag(n) = gflag
         ptcntsum = global_sum(ptcnt(n),distrb_info)
         ptcnt(n) = ptcntsum
         failcntsum = global_sum(failcnt(n),distrb_info)
         failcnt(n) = failcntsum
         if (errorflag(n) == failflag) errorflag0 = failflag
      enddo

      if (my_task == master_task) then
         write(6,*) ' '
         write(6,*) 'GATHSCATCHK COMPLETED SUCCESSFULLY'
         write(6,*) ' '
         tpcnt = 0
         tfcnt = 0
         do n = 1,testcnt
            if (errorflag(n) == passflag) then
               tpcnt = tpcnt + 1
               write(6,'(2a,2i9)') 'PASS ',trim(teststring(n)),ptcnt(n),failcnt(n)
            else
               tfcnt = tfcnt + 1
               write(6,'(2a,2i9)') 'FAIL ',trim(teststring(n)),ptcnt(n),failcnt(n)
            endif
         enddo
         write(6,*) ' '
         write(6,*) ' total pass = ',tpcnt
         write(6,*) ' total fail = ',tfcnt
         write(6,*) ' '
         if (errorflag0 == passflag) then
            write(6,*) 'GATHSCATCHK TEST COMPLETED SUCCESSFULLY'
         else
            write(6,*) 'GATHSCATCHK TEST FAILED'
         endif
         write(6,*) ' '
         write(6,*) '=========================================================='
         write(6,*) ' '
      endif

      !-----------------------------------------------------------------
      ! Gracefully end
      !-----------------------------------------------------------------

      call end_run()

      end program gathscatchk

!=======================================================================
