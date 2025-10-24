
      module halochk_data

      use CICE_InitMod
      use ice_kinds_mod, only: int_kind, dbl_kind, real_kind, log_kind
      use ice_blocks, only: block, get_block, nx_block, ny_block, nblocks_tot, nghost
      use ice_boundary, only: ice_HaloUpdate, ice_HaloUpdate_stress
      use ice_constants, only: c0, c1, p5, &
          field_loc_unknown, field_loc_noupdate, &
          field_loc_center, field_loc_NEcorner, &
          field_loc_Nface, field_loc_Eface, &
          field_type_unknown, field_type_noupdate, &
          field_type_scalar, field_type_vector, field_type_angle
      use ice_communicate, only: my_task, master_task, get_num_procs, MPI_COMM_ICE
      use ice_distribution, only: ice_distributionGetBlockID, ice_distributionGet
      use ice_domain_size, only: nx_global, ny_global, &
          block_size_x, block_size_y, max_blocks
      use ice_domain, only: distrb_info, halo_info, &
          ew_boundary_type, ns_boundary_type
      use ice_exit, only: abort_ice, end_run
      use ice_global_reductions, only: global_minval, global_maxval, global_sum

      implicit none

      integer(int_kind), parameter ::  &
         passflag = 0, &
         failflag = 1

      end module halochk_data

!=======================================================================

      program halochk

      ! This tests the CICE halo update methods by
      ! using CICE_InitMod (from the standalone model) to read/initialize
      ! a CICE grid/configuration.

      use halochk_data

      implicit none

      integer(int_kind) :: nn, nl, nt, nf, i, j, k1, k2, n, ib, ie, jb, je
      integer(int_kind) :: iblock, itrip, ioffset, joffset
      integer(int_kind) :: blockID, numBlocks, jtrip
      type (block) :: this_block

      ! fields sent to the haloupdate
      real(dbl_kind)   , allocatable :: darrayi1(:,:,:)    , darrayj1(:,:,:)
      real(dbl_kind)   , allocatable :: darrayi2(:,:,:,:)  , darrayj2(:,:,:,:)
      real(dbl_kind)   , allocatable :: darrayi3(:,:,:,:,:), darrayj3(:,:,:,:,:)
      real(real_kind)  , allocatable :: rarrayi1(:,:,:)    , rarrayj1(:,:,:)
      real(real_kind)  , allocatable :: rarrayi2(:,:,:,:)  , rarrayj2(:,:,:,:)
      real(real_kind)  , allocatable :: rarrayi3(:,:,:,:,:), rarrayj3(:,:,:,:,:)
      integer(int_kind), allocatable :: iarrayi1(:,:,:)    , iarrayj1(:,:,:)
      integer(int_kind), allocatable :: iarrayi2(:,:,:,:)  , iarrayj2(:,:,:,:)
      integer(int_kind), allocatable :: iarrayi3(:,:,:,:,:), iarrayj3(:,:,:,:,:)
      logical(log_kind), allocatable :: larrayi1(:,:,:)    , larrayj1(:,:,:)
      real(dbl_kind)   , allocatable :: darrayi1str(:,:,:) , darrayj1str(:,:,:)
      real(dbl_kind)   , allocatable :: darrayi10(:,:,:)   , darrayj10(:,:,:)

      ! expected results
      real(dbl_kind), allocatable :: cidata_bas(:,:,:,:,:),cjdata_bas(:,:,:,:,:)  ! baseline

      integer(int_kind), parameter :: maxtests = 11
      integer(int_kind), parameter :: maxtypes = 4
      integer(int_kind), parameter :: maxlocs = 5
      integer(int_kind), parameter :: maxfills = 2
      integer(int_kind), parameter :: nz1 = 3
      integer(int_kind), parameter :: nz2 = 4
      real(dbl_kind)    :: aichk,ajchk,cichk,cjchk,rival,rjval,rsign
      real(dbl_kind)    :: fillexpected
      character(len=16) :: locs_name(maxlocs), types_name(maxtypes), fill_name(maxfills)
      integer(int_kind) :: field_loc(maxlocs), field_type(maxtypes)
      logical :: halofill
      integer(int_kind) :: npes, ierr, ntask, testcnt, tottest, tpcnt, tfcnt
      integer(int_kind) :: errorflag0, gflag, k1m, k2m, ptcntsum, failcntsum
      integer(int_kind), allocatable :: errorflag(:)
      integer(int_kind), allocatable :: ptcnt(:), failcnt(:)
      character(len=128), allocatable :: teststring(:)
      character(len=32) :: halofld
      logical :: tripole_average, tripole_pole
      logical :: first_call = .true.

      real(dbl_kind)   , parameter :: fillval = -88888.0_dbl_kind
      real(dbl_kind)   , parameter :: dhalofillval = -999.0_dbl_kind
      real(real_kind)  , parameter :: rhalofillval = -999.0_real_kind
      integer(int_kind), parameter :: ihalofillval = -999
      character(len=*) , parameter :: subname='(halochk)'

      !-----------------------------------------------------------------
      ! Initialize CICE
      !-----------------------------------------------------------------

      call CICE_Initialize
      npes = get_num_procs()

      locs_name (:) = 'unknown'
      types_name(:) = 'unknown'
      fill_name (:) = 'unknown'
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
      fill_name (2) = 'nofill'

      tottest = maxtests * maxlocs * maxtypes * maxfills
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
         write(6,*) 'RunningUnitTest HALOCHK'
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

      errorflag0    = passflag
      errorflag(:)  = passflag
      teststring(:) = ' '

      ! ---------------------------
      ! TEST HALO UPDATE
      ! ---------------------------

      if (my_task == master_task) write(6,*) ' '

      allocate(darrayi1   (nx_block,ny_block,max_blocks))
      allocate(darrayj1   (nx_block,ny_block,max_blocks))
      allocate(darrayi2   (nx_block,ny_block,nz1,max_blocks))
      allocate(darrayj2   (nx_block,ny_block,nz1,max_blocks))
      allocate(darrayi3   (nx_block,ny_block,nz1,nz2,max_blocks))
      allocate(darrayj3   (nx_block,ny_block,nz1,nz2,max_blocks))
      allocate(rarrayi1   (nx_block,ny_block,max_blocks))
      allocate(rarrayj1   (nx_block,ny_block,max_blocks))
      allocate(rarrayi2   (nx_block,ny_block,nz1,max_blocks))
      allocate(rarrayj2   (nx_block,ny_block,nz1,max_blocks))
      allocate(rarrayi3   (nx_block,ny_block,nz1,nz2,max_blocks))
      allocate(rarrayj3   (nx_block,ny_block,nz1,nz2,max_blocks))
      allocate(iarrayi1   (nx_block,ny_block,max_blocks))
      allocate(iarrayj1   (nx_block,ny_block,max_blocks))
      allocate(iarrayi2   (nx_block,ny_block,nz1,max_blocks))
      allocate(iarrayj2   (nx_block,ny_block,nz1,max_blocks))
      allocate(iarrayi3   (nx_block,ny_block,nz1,nz2,max_blocks))
      allocate(iarrayj3   (nx_block,ny_block,nz1,nz2,max_blocks))
      allocate(larrayi1   (nx_block,ny_block,max_blocks))
      allocate(larrayj1   (nx_block,ny_block,max_blocks))
      allocate(darrayi1str(nx_block,ny_block,max_blocks))
      allocate(darrayj1str(nx_block,ny_block,max_blocks))
      allocate(darrayi10  (nx_block,ny_block,max_blocks))
      allocate(darrayj10  (nx_block,ny_block,max_blocks))

      allocate(cidata_bas(nx_block,ny_block,nz1,nz2,max_blocks))
      allocate(cjdata_bas(nx_block,ny_block,nz1,nz2,max_blocks))

      darrayi1 = fillval
      darrayj1 = fillval
      darrayi2 = fillval
      darrayj2 = fillval
      darrayi3 = fillval
      darrayj3 = fillval
      rarrayi1 = fillval
      rarrayj1 = fillval
      rarrayi2 = fillval
      rarrayj2 = fillval
      rarrayi3 = fillval
      rarrayj3 = fillval
      iarrayi1 = fillval
      iarrayj1 = fillval
      iarrayi2 = fillval
      iarrayj2 = fillval
      iarrayi3 = fillval
      iarrayj3 = fillval
      larrayi1 = .false.
      larrayj1 = .true.
      darrayi1str = fillval
      darrayj1str = fillval
      darrayi10  = fillval
      darrayj10  = fillval
      cidata_bas = fillval
      cjdata_bas = fillval

      call ice_distributionGet(distrb_info, numLocalBlocks = numBlocks)

      !--- baseline data ---
      ! set to the global index
      ! i/j valid everywhere for "cyclic"
      ! i/j valid for "open" with extrapolation on outer boundary
      ! i/j zero on outer boundary for "closed"

      do iblock = 1,numBlocks
         call ice_distributionGetBlockID(distrb_info, iblock, blockID)
         this_block = get_block(blockID, blockID)
         do k2 = 1,nz2
         do k1 = 1,nz1
         do j = 1,ny_block
         do i = 1,nx_block
            cidata_bas(i,j,k1,k2,iblock) = real(this_block%i_glob(i),kind=dbl_kind) + &
                   real(k1,kind=dbl_kind)*1000._dbl_kind + real(k2,kind=dbl_kind)*10000._dbl_kind
            cjdata_bas(i,j,k1,k2,iblock) = real(this_block%j_glob(j),kind=dbl_kind) + &
                   real(k1,kind=dbl_kind)*1000._dbl_kind + real(k2,kind=dbl_kind)*10000._dbl_kind
         enddo
         enddo
         enddo
         enddo
      enddo

      !---------------------------------------------------------------

      testcnt = 0
      do nn = 1, maxtests
      do nl = 1, maxlocs
      do nt = 1, maxtypes
      do nf = 1, maxfills

         !--- setup test ---
         first_call = .true.
         testcnt = testcnt + 1
         if (nf == 1) then
            halofill = .true.
            fillexpected = dhalofillval
         elseif (nf == 2) then
            halofill = .false.
            fillexpected = fillval
         else
            write(6,*) subname,' nf = ',nf
            if (my_task == master_task) then
               write(6,*) ' '
               write(6,*) 'HALOCHK FAILED'
               write(6,*) ' '
            endif
            call abort_ice(subname//' invalid value of nf',file=__FILE__,line=__LINE__)
         endif
         if (testcnt > tottest) then
            if (my_task == master_task) then
               write(6,*) ' '
               write(6,*) 'HALOCHK FAILED'
               write(6,*) ' '
            endif
            call abort_ice(subname//' testcnt > tottest',file=__FILE__,line=__LINE__)
         endif

         !--- fill arrays ---
         darrayi1(:,:,:) = fillval
         darrayj1(:,:,:) = fillval
         darrayi2(:,:,:,:) = fillval
         darrayj2(:,:,:,:) = fillval
         darrayi3(:,:,:,:,:) = fillval
         darrayj3(:,:,:,:,:) = fillval
         darrayi1str(:,:,:) = fillval
         darrayj1str(:,:,:) = fillval
         do iblock = 1,numBlocks
            call ice_distributionGetBlockID(distrb_info, iblock, blockID)
            this_block = get_block(blockID, blockID)
            ib = this_block%ilo
            ie = this_block%ihi
            jb = this_block%jlo
            je = this_block%jhi
            do j = jb,je
               do i = ib,ie
                  k1 = 1
                  k2 = 1
                  darrayi1(i,j,iblock) = cidata_bas(i,j,k1,k2,iblock)
                  darrayj1(i,j,iblock) = cjdata_bas(i,j,k1,k2,iblock)
                  do k1 = 1,nz1
                     k2 = 1
                     darrayi2(i,j,k1,iblock) = cidata_bas(i,j,k1,k2,iblock)
                     darrayj2(i,j,k1,iblock) = cjdata_bas(i,j,k1,k2,iblock)
                     do k2 = 1,nz2
                        darrayi3(i,j,k1,k2,iblock) = cidata_bas(i,j,k1,k2,iblock)
                        darrayj3(i,j,k1,k2,iblock) = cjdata_bas(i,j,k1,k2,iblock)
                     enddo
                  enddo
               enddo
            enddo
         enddo

         ! copy original darray1 for "stress" compare
         darrayi10 = darrayi1
         darrayj10 = darrayj1

         !--- halo update ---

         if (nn == 1) then
            k1m = 1
            k2m = 1
            halofld = '2DR8'
            if (halofill) then
               call ice_haloUpdate(darrayi1, halo_info, field_loc(nl), field_type(nt), fillvalue=dhalofillval)
               call ice_haloUpdate(darrayj1, halo_info, field_loc(nl), field_type(nt), fillvalue=dhalofillval)
            else
               call ice_haloUpdate(darrayi1, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate(darrayj1, halo_info, field_loc(nl), field_type(nt))
            endif
         elseif (nn == 2) then
            k1m = nz1
            k2m = 1
            halofld = '3DR8'
            if (halofill) then
               call ice_haloUpdate(darrayi2, halo_info, field_loc(nl), field_type(nt), fillvalue=dhalofillval)
               call ice_haloUpdate(darrayj2, halo_info, field_loc(nl), field_type(nt), fillvalue=dhalofillval)
            else
               call ice_haloUpdate(darrayi2, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate(darrayj2, halo_info, field_loc(nl), field_type(nt))
            endif
         elseif (nn == 3) then
            k1m = nz1
            k2m = nz2
            halofld = '4DR8'
            if (halofill) then
               call ice_haloUpdate(darrayi3, halo_info, field_loc(nl), field_type(nt), fillvalue=dhalofillval)
               call ice_haloUpdate(darrayj3, halo_info, field_loc(nl), field_type(nt), fillvalue=dhalofillval)
            else
               call ice_haloUpdate(darrayi3, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate(darrayj3, halo_info, field_loc(nl), field_type(nt))
            endif
         elseif (nn == 4) then
            k1m = 1
            k2m = 1
            halofld = '2DR4'
            rarrayi1 = real(darrayi1,kind=real_kind)
            rarrayj1 = real(darrayj1,kind=real_kind)
            if (halofill) then
               call ice_haloUpdate(rarrayi1, halo_info, field_loc(nl), field_type(nt), fillvalue=rhalofillval)
               call ice_haloUpdate(rarrayj1, halo_info, field_loc(nl), field_type(nt), fillvalue=rhalofillval)
            else
               call ice_haloUpdate(rarrayi1, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate(rarrayj1, halo_info, field_loc(nl), field_type(nt))
            endif
            darrayi1 = real(rarrayi1,kind=dbl_kind)
            darrayj1 = real(rarrayj1,kind=dbl_kind)
         elseif (nn == 5) then
            k1m = nz1
            k2m = 1
            halofld = '3DR4'
            rarrayi2 = real(darrayi2,kind=real_kind)
            rarrayj2 = real(darrayj2,kind=real_kind)
            if (halofill) then
               call ice_haloUpdate(rarrayi2, halo_info, field_loc(nl), field_type(nt), fillvalue=rhalofillval)
               call ice_haloUpdate(rarrayj2, halo_info, field_loc(nl), field_type(nt), fillvalue=rhalofillval)
            else
               call ice_haloUpdate(rarrayi2, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate(rarrayj2, halo_info, field_loc(nl), field_type(nt))
            endif
            darrayi2 = real(rarrayi2,kind=dbl_kind)
            darrayj2 = real(rarrayj2,kind=dbl_kind)
         elseif (nn == 6) then
            k1m = nz1
            k2m = nz2
            halofld = '4DR4'
            rarrayi3 = real(darrayi3,kind=real_kind)
            rarrayj3 = real(darrayj3,kind=real_kind)
            if (halofill) then
               call ice_haloUpdate(rarrayi3, halo_info, field_loc(nl), field_type(nt), fillvalue=rhalofillval)
               call ice_haloUpdate(rarrayj3, halo_info, field_loc(nl), field_type(nt), fillvalue=rhalofillval)
            else
               call ice_haloUpdate(rarrayi3, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate(rarrayj3, halo_info, field_loc(nl), field_type(nt))
            endif
            darrayi3 = real(rarrayi3,kind=dbl_kind)
            darrayj3 = real(rarrayj3,kind=dbl_kind)
         elseif (nn == 7) then
            k1m = 1
            k2m = 1
            halofld = '2DI4'
            iarrayi1 = nint(darrayi1)
            iarrayj1 = nint(darrayj1)
            if (halofill) then
               call ice_haloUpdate(iarrayi1, halo_info, field_loc(nl), field_type(nt), fillvalue=ihalofillval)
               call ice_haloUpdate(iarrayj1, halo_info, field_loc(nl), field_type(nt), fillvalue=ihalofillval)
            else
               call ice_haloUpdate(iarrayi1, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate(iarrayj1, halo_info, field_loc(nl), field_type(nt))
            endif
            darrayi1 = real(iarrayi1,kind=dbl_kind)
            darrayj1 = real(iarrayj1,kind=dbl_kind)
         elseif (nn == 8) then
            k1m = nz1
            k2m = 1
            halofld = '3DI4'
            iarrayi2 = nint(darrayi2)
            iarrayj2 = nint(darrayj2)
            if (halofill) then
               call ice_haloUpdate(iarrayi2, halo_info, field_loc(nl), field_type(nt), fillvalue=ihalofillval)
               call ice_haloUpdate(iarrayj2, halo_info, field_loc(nl), field_type(nt), fillvalue=ihalofillval)
            else
               call ice_haloUpdate(iarrayi2, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate(iarrayj2, halo_info, field_loc(nl), field_type(nt))
            endif
            darrayi2 = real(iarrayi2,kind=dbl_kind)
            darrayj2 = real(iarrayj2,kind=dbl_kind)
         elseif (nn == 9) then
            k1m = nz1
            k2m = nz2
            halofld = '4DI4'
            iarrayi3 = nint(darrayi3)
            iarrayj3 = nint(darrayj3)
            if (halofill) then
               call ice_haloUpdate(iarrayi3, halo_info, field_loc(nl), field_type(nt), fillvalue=ihalofillval)
               call ice_haloUpdate(iarrayj3, halo_info, field_loc(nl), field_type(nt), fillvalue=ihalofillval)
            else
               call ice_haloUpdate(iarrayi3, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate(iarrayj3, halo_info, field_loc(nl), field_type(nt))
            endif
            darrayi3 = real(iarrayi3,kind=dbl_kind)
            darrayj3 = real(iarrayj3,kind=dbl_kind)
         elseif (nn == 10) then
            k1m = 1
            k2m = 1
            halofld = '2DL1'
            where (darrayi1 == fillval) 
               larrayi1 = .false.
            elsewhere
               larrayi1 = (mod(nint(darrayi1),2) == 1) 
            endwhere
            where (darrayj1 == fillval) 
               larrayj1 = .true.
            elsewhere
               larrayj1 = (mod(nint(darrayj1),2) == 1) 
            endwhere
            if (halofill) then
               call ice_haloUpdate(larrayi1, halo_info, field_loc(nl), field_type(nt), fillvalue=0)
               call ice_haloUpdate(larrayj1, halo_info, field_loc(nl), field_type(nt), fillvalue=1)
            else
               call ice_haloUpdate(larrayi1, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate(larrayj1, halo_info, field_loc(nl), field_type(nt))
            endif
            darrayi1 = c0
            where (larrayi1) darrayi1 = c1
            darrayj1 = c0
            where (larrayj1) darrayj1 = c1
         elseif (nn == 11) then
            k1m = 1
            k2m = 1
            halofld = 'STRESS'
            darrayi1str = -darrayi1  ! flip sign for testing
            darrayj1str = -darrayj1
            if (halofill) then
               call ice_haloUpdate_stress(darrayi1, darrayi1str, halo_info, field_loc(nl), field_type(nt), fillvalue=dhalofillval)
               call ice_haloUpdate_stress(darrayj1, darrayj1str, halo_info, field_loc(nl), field_type(nt), fillvalue=dhalofillval)
            else
               call ice_haloUpdate_stress(darrayi1, darrayi1str, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate_stress(darrayj1, darrayj1str, halo_info, field_loc(nl), field_type(nt))
            endif
         endif

         write(teststring(testcnt),'(6a8)') trim(halofld),trim(locs_name(nl)),trim(types_name(nt)),trim(fill_name(nf)), &
                      trim(ew_boundary_type),trim(ns_boundary_type)

         do iblock = 1,numBlocks
            call ice_distributionGetBlockID(distrb_info, iblock, blockID)
            this_block = get_block(blockID, blockID)
            ib = this_block%ilo
            ie = this_block%ihi
            jb = this_block%jlo
            je = this_block%jhi
            ! just check non-padded gridcells
            do j = jb-nghost, je+nghost
            do i = ib-nghost, ie+nghost
            do k1 = 1,k1m
            do k2 = 1,k2m
               tripole_average = .false.
               tripole_pole = .false.
               if (index(halofld,'2D') > 0) then
                  aichk = darrayi1(i,j,iblock)
                  ajchk = darrayj1(i,j,iblock)
               elseif (index(halofld,'STRESS') > 0) then
                  aichk = darrayi1(i,j,iblock)
                  ajchk = darrayj1(i,j,iblock)
               elseif (index(halofld,'3D') > 0) then
                  aichk = darrayi2(i,j,k1,iblock)
                  ajchk = darrayj2(i,j,k1,iblock)
               elseif (index(halofld,'4D') > 0) then
                  aichk = darrayi3(i,j,k1,k2,iblock)
                  ajchk = darrayj3(i,j,k1,k2,iblock)
               else
                  if (my_task == master_task) then
                     write(6,*) ' '
                     write(6,*) 'HALOCHK FAILED'
                     write(6,*) ' '
                  endif
                  call abort_ice(subname//' halofld not matched '//trim(halofld),file=__FILE__,line=__LINE__)
               endif

               cichk = cidata_bas(i,j,k1,k2,iblock)
               cjchk = cjdata_bas(i,j,k1,k2,iblock)

               ! halo special cases

               if (field_loc (nl) == field_loc_noupdate .or. &
                   field_type(nt) == field_type_noupdate) then
                  if (i < ib .or. j < jb .or. i > ie .or. j > je) then
                     ! no halo update anywhere, doesn't even see fillvalue passed in
                     cichk = fillval
                     cjchk = fillval
                  endif

               else
                  ! if ew_boundary_type is not cyclic we expect just fill values on outer boundary
                  if (ew_boundary_type /= 'cyclic' .and. &
                      ((this_block%i_glob(ib) == 1         .and. i < ib) .or. &  ! west outer face
                       (this_block%i_glob(ie) == nx_global .and. i > ie))) then  ! east outer face
                     cichk = fillexpected
                     cjchk = fillexpected
                  endif

                  ! if ns_boundary_type is not cyclic we expect just fill values on outer boundary except
                  ! - tripole north edge will be haloed and is updated below, default to fill value for now
                  ! - tripole south edge will be set to the fillvalue or to haloupdate internal default (c0)
                  !   tripole basically assumes south edge is land or always ice free in CICE
                  if (ns_boundary_type /= 'cyclic' .and. &
                      ((this_block%j_glob(jb) == 1         .and. j < jb) .or. &  ! south outer face
                       (this_block%j_glob(je) == ny_global .and. j > je))) then  ! north outer face
                     ! ns_boundary_type is not cyclic and on outer boundary
                     if ((ns_boundary_type == 'tripole' .or. &
                          ns_boundary_type == 'tripoleT') .and. &
                         .not. halofill) then
                        cichk = c0
                        cjchk = c0
                     else
                        cichk = fillexpected
                        cjchk = fillexpected
                     endif
                  endif

                  if (index(halofld,'STRESS') > 0) then
                     ! only updates on tripole zipper for tripole grids
                     ! darrayi10 is copy of darrayi1 before halo call
                     cichk = darrayi10(i,j,iblock)
                     cjchk = darrayj10(i,j,iblock)
                  endif

                  !--- tripole on north boundary, need to hardcode ---
                  !--- tripole and tripoleT slightly different     ---
                  !--- establish special set of points here        ---
                  if ((this_block%j_glob(je) == ny_global) .and. &
                     ((ns_boundary_type == 'tripole'  .and. &
                       (j > je .or. &
                        (j == je .and. (field_loc(nl) == field_loc_Nface .or. field_loc(nl) == field_loc_NEcorner)))) .or. &
                      (ns_boundary_type == 'tripoleT' .and. &
                       (j >= je)))) then

                     ! flip sign for vector/angle except for logical halo updates
                     rsign = c1
                     if ((field_type(nt) == field_type_vector .or. field_type(nt) == field_type_angle) .and. &
                         .not. (index(halofld,'L1') > 0)) then
                        rsign = -c1
                     endif

                     ! for tripole
                     if (ns_boundary_type == 'tripole') then

                        ! compute itrip and jtrip, these are the location where the halo values are defined for i,j
                        ! for j=je averaging, itrip and jtrip are the 2nd gridpoint associated with averaging

                        ! standard center tripole u-fold
                        itrip = nx_global-this_block%i_glob(i)+1
                        jtrip = max(je - (j-je) + 1 , je)
                        ioffset = 0
                        joffset = 0

                        if (field_loc(nl) == field_loc_NEcorner .or. field_loc(nl) == field_loc_Nface) then
                           ! need j offset
                           joffset = -1
                           if (j == je) then
                              tripole_average = .true.
                           endif
                        endif

                        if (field_loc(nl) == field_loc_NEcorner .or. field_loc(nl) == field_loc_Eface) then
                           ! fold plus cell offset
                           ioffset = -1
                           ! CICE treats j=ny_global tripole edge points incorrectly
                           ! should do edge wraparound and average
                           ! CICE does not update those points, assumes it's "land"
                           if (j == je) then
                              if (this_block%i_glob(i) == nx_global/2) then
                                  tripole_pole = .true.
                              elseif (this_block%i_glob(i) == nx_global  ) then
                                  tripole_pole = .true.
                              endif
                           endif
                        endif

                     ! for tripoleT
                     elseif (ns_boundary_type == 'tripoleT') then

                        ! compute itrip and jtrip, these are the location where the halo values are defined for i,j
                        ! for j=je averaging, itrip and jtrip are the 2nd gridpoint associated with averaging

                        ! standard center tripoleT t-fold
                        itrip = nx_global-this_block%i_glob(i)+2
                        jtrip = je - (j-je)
                        ioffset = 0
                        joffset = 0

                        if (field_loc(nl) == field_loc_NEcorner .or. field_loc(nl) == field_loc_Eface) then
                           ! fold plus cell offset
                           ioffset = -1
                        endif

                        if (field_loc(nl) == field_loc_NEcorner .or. field_loc(nl) == field_loc_Nface) then
                           ! need j offset
                           joffset = -1
                        endif

                        if (field_loc(nl) == field_loc_Center .or. field_loc(nl) == field_loc_Eface) then
                           if (j == je) then
                              tripole_average = .true.
                           endif
                        endif

                        ! center point poles need to be treated special
                        if (field_loc(nl) == field_loc_Center) then
                           if (j == je .and. &
                              (this_block%i_glob(i) == 1 .or. this_block%i_glob(i) == nx_global/2+1)) then
                              tripole_pole = .true.
                           endif
                        endif

                     endif

                     itrip = mod(itrip + ioffset + nx_global-1,nx_global)+1
                     jtrip = jtrip + joffset

                     rival = (real(itrip,kind=dbl_kind) + &
                              real(k1,kind=dbl_kind)*1000._dbl_kind + real(k2,kind=dbl_kind)*10000._dbl_kind)
                     rjval = (real(this_block%j_glob(jtrip),kind=dbl_kind) + &
                              real(k1,kind=dbl_kind)*1000._dbl_kind + real(k2,kind=dbl_kind)*10000._dbl_kind)

                     if (index(halofld,'STRESS') > 0) then
                        ! only updates on tripole zipper for tripole grids, not tripoleT
                        ! note: L1 and STRESS never overlap so don't worry about L1 here
                        if (tripole_pole) then
                           ! flip sign due to sign of darrayi1str
                           ! ends of tripole seam not averaged in CICE
                           cichk = -rsign * cidata_bas(i,j,k1,k2,iblock)
                           cjchk = -rsign * cjdata_bas(i,j,k1,k2,iblock)
                        else
                           cichk = -rsign * rival
                           cjchk = -rsign * rjval
                        endif

                     elseif (tripole_pole) then
                        ! ends of tripole seam not averaged in CICE
                        cichk = rsign * cidata_bas(i,j,k1,k2,iblock)
                        cjchk = rsign * cjdata_bas(i,j,k1,k2,iblock)

                     elseif (tripole_average) then
                        if (index(halofld,'L1') > 0) then
                           ! logical math doesn't work this way, force to correct answer
                           cichk = aichk ! p5 * (mod(nint(cidata_bas(i,j,k1,k2,iblock)),2) + rsign * mod(nint(rival),2))
                           cjchk = ajchk ! p5 * (mod(nint(cidata_bas(i,j,k1,k2,iblock)),2) + rsign * mod(nint(rjval),2))
                        else
                           cichk = p5 * (cidata_bas(i,j,k1,k2,iblock) + rsign * rival)
                           cjchk = p5 * (cjdata_bas(i,j,k1,k2,iblock) + rsign * rjval)
                        endif

                     else
                        ! standard tripole fold
                        cichk = rsign * rival
                        cjchk = rsign * rjval
                     endif

                  endif  ! tripole or tripoleT

               endif

               if (index(halofld,'I4') > 0) then
                  cichk = real(nint(cichk),kind=dbl_kind)
                  cjchk = real(nint(cjchk),kind=dbl_kind)
               endif

               if (index(halofld,'L1') > 0) then
                  if (cichk == dhalofillval .or. cichk == fillval) then
                     cichk = c0
                  else
                     cichk = mod(nint(cichk),2)
                  endif
                  if (cjchk == dhalofillval .or. cjchk == fillval) then
                     cjchk = c1
                  else
                     cjchk = mod(nint(cjchk),2)
                  endif
               endif

               ptcnt(testcnt) = ptcnt(testcnt) + 1
               call chkresults(aichk,cichk,errorflag(testcnt),testcnt,failcnt(testcnt), &
                    i,j,k1,k2,iblock,first_call,teststring(testcnt),trim(halofld)//'_I')
               call chkresults(ajchk,cjchk,errorflag(testcnt),testcnt,failcnt(testcnt), &
                    i,j,k1,k2,iblock,first_call,teststring(testcnt),trim(halofld)//'_J')
            enddo  ! k2
            enddo  ! k1
            enddo  ! i
            enddo  ! j
         enddo  ! iblock

      enddo  ! maxfills
      enddo  ! maxtypes
      enddo  ! maxlocs
      enddo  ! maxtests

      ! ---------------------------
      ! SUMMARY
      ! ---------------------------

      do n = 1,tottest
         gflag = global_maxval(errorflag(n), MPI_COMM_ICE)
         errorflag(n) = gflag
         ptcntsum = global_sum(ptcnt(n),distrb_info)
         ptcnt(n) = ptcntsum
         failcntsum = global_sum(failcnt(n),distrb_info)
         failcnt(n) = failcntsum
      enddo
      errorflag0 = maxval(errorflag(:))

      if (my_task == master_task) then
         write(6,*) ' '
         write(6,*) 'HALOCHK COMPLETED SUCCESSFULLY'
         write(6,*) ' '
         tpcnt = 0
         tfcnt = 0
         do n = 1,tottest
            if (errorflag(n) == passflag) then
               tpcnt = tpcnt + 1
               write(6,'(2a,2i8)') 'PASS ',trim(teststring(n)),ptcnt(n),failcnt(n)
            else
               tfcnt = tfcnt + 1
               write(6,'(2a,2i8)') 'FAIL ',trim(teststring(n)),ptcnt(n),failcnt(n)
            endif
         enddo
         write(6,*) ' '
         write(6,*) ' total pass = ',tpcnt
         write(6,*) ' total fail = ',tfcnt
         write(6,*) ' '
         if (errorflag0 == passflag) then
            write(6,*) 'HALOCHK TEST COMPLETED SUCCESSFULLY'
         else
            write(6,*) 'HALOCHK TEST FAILED'
         endif
         write(6,*) ' '
         write(6,*) '=========================================================='
         write(6,*) ' '
      endif


      !-----------------------------------------------------------------
      ! Gracefully end
      !-----------------------------------------------------------------

      call end_run()

      end program halochk

!=======================================================================

      subroutine chkresults(a1,r1,errorflag,testcnt,failcnt,i,j,k1,k2,iblock,first_call,teststring,halofld)

      use halochk_data

      implicit none

      real(dbl_kind)   , intent(in)    :: a1,r1
      integer(int_kind), intent(inout) :: errorflag, failcnt
      integer(int_kind), intent(in)    :: i,j,k1,k2,iblock,testcnt
      logical          , intent(inout) :: first_call
      character(len=*) , intent(in)    :: teststring,halofld

      logical,parameter :: print_always = .false.
      character(len=*) , parameter :: subname='(chkresults)'

      if (a1 /= r1 .or. print_always) then
         if (a1 /= r1) then
            errorflag = failflag
            failcnt = failcnt + 1
         endif
         if (first_call) then
            write(100+my_task,*) ' '
            write(100+my_task,'(a,i4,2a)') '------- TEST = ',testcnt,' ',trim(teststring)
            write(100+my_task,*) ' '
            write(100+my_task,'(a)') '           test  task    i     j    k1    k2  iblock  expected   halocomp       diff'
            first_call = .false.
         endif
         write(100+my_task,1001) trim(halofld),testcnt,my_task,i,j,k1,k2,iblock,r1,a1,r1-a1
      endif

 1001 format(a8,7i6,3f12.3)

      end subroutine chkresults
!=======================================================================
