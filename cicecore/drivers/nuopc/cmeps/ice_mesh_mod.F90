module ice_mesh_mod

  use ESMF
  use NUOPC            , only : NUOPC_CompAttributeGet
  use ice_kinds_mod    , only : dbl_kind, int_kind, char_len, char_len_long
  use ice_domain_size  , only : nx_global, ny_global, max_blocks
  use ice_domain       , only : nblocks, blocks_ice, distrb_info
  use ice_blocks       , only : block, get_block, nx_block, ny_block, nblocks_x, nblocks_y
  use ice_shr_methods  , only : chkerr
  use ice_fileunits    , only : nu_diag
  use ice_communicate  , only : my_task, master_task
  use ice_exit         , only : abort_ice
  use icepack_intfc    , only : icepack_query_parameters
  use icepack_intfc    , only : icepack_warnings_flush, icepack_warnings_aborted
  implicit none
  private

  public  :: ice_mesh_set_distgrid
  public  :: ice_mesh_setmask_from_maskfile
  public  :: ice_mesh_check

  private :: ice_mesh_create_mask

  ! Only relevant for lat-lon grids gridcell value of [1 - (land fraction)] (T-cell)
  real (dbl_kind), allocatable, public :: ocn_gridcell_frac(:,:,:)

  character(*), parameter :: u_FILE_u = &
       __FILE__

!=======================================================================
contains
!=======================================================================

  subroutine ice_mesh_set_distgrid(localpet, npes, distgrid, rc)

    ! Determine the global index space needed for the distgrid

    ! input/output variables
    integer             , intent(in)    :: localpet
    integer             , intent(in)    :: npes
    type(ESMF_DistGrid) , intent(inout) :: distgrid
    integer             , intent(out)   :: rc

    ! local variables
    integer               :: n,c,g,i,j,m        ! indices
    integer               :: iblk, jblk         ! indices
    integer               :: ig, jg             ! indices
    integer               :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    integer               :: lsize              ! local size of coupling array
    type(block)           :: this_block         ! block information for current block
    integer               :: num_elim_global
    integer               :: num_elim_local
    integer               :: num_elim
    integer               :: num_ice
    integer               :: num_elim_gcells    ! local number of eliminated gridcells
    integer               :: num_elim_blocks    ! local number of eliminated blocks
    integer               :: num_total_blocks
    integer               :: my_elim_start, my_elim_end
    integer , allocatable :: gindex(:)
    integer , allocatable :: gindex_ice(:)
    integer , allocatable :: gindex_elim(:)
    integer               :: globalID
    character(len=*), parameter :: subname = ' ice_mesh_set_distgrid: '
    !----------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! number the local grid to get allocation size for gindex_ice
    lsize = 0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
          do i = ilo, ihi
             lsize = lsize + 1
          enddo
       enddo
    enddo

    ! set global index array
    allocate(gindex_ice(lsize))
    n = 0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
          do i = ilo, ihi
             n = n+1
             ig = this_block%i_glob(i)
             jg = this_block%j_glob(j)
             gindex_ice(n) = (jg-1)*nx_global + ig
          enddo
       enddo
    enddo

    ! Determine total number of eliminated blocks globally
    globalID = 0
    num_elim_global = 0  ! number of eliminated blocks
    num_total_blocks = 0
    do jblk=1,nblocks_y
       do iblk=1,nblocks_x
          globalID = globalID + 1
          num_total_blocks = num_total_blocks + 1
          if (distrb_info%blockLocation(globalID) == 0) then
             num_elim_global = num_elim_global + 1
          end if
       end do
    end do

    if (num_elim_global > 0) then

       ! Distribute the eliminated blocks in a round robin fashion amoung processors
       num_elim_local = num_elim_global / npes
       my_elim_start = num_elim_local*localPet + min(localPet, mod(num_elim_global, npes)) + 1
       if (localPet < mod(num_elim_global, npes)) then
          num_elim_local = num_elim_local + 1
       end if
       my_elim_end = my_elim_start + num_elim_local - 1

       ! Determine the number of eliminated gridcells locally
       globalID = 0
       num_elim_blocks = 0  ! local number of eliminated blocks
       num_elim_gcells = 0
       do jblk=1,nblocks_y
          do iblk=1,nblocks_x
             globalID = globalID + 1
             if (distrb_info%blockLocation(globalID) == 0) then
                num_elim_blocks = num_elim_blocks + 1
                if (num_elim_blocks >= my_elim_start .and. num_elim_blocks <= my_elim_end) then
                   this_block = get_block(globalID, globalID)
                   num_elim_gcells = num_elim_gcells + &
                        (this_block%jhi-this_block%jlo+1) * (this_block%ihi-this_block%ilo+1)
                end if
             end if
          end do
       end do

       ! Determine the global index space of the eliminated gridcells
       allocate(gindex_elim(num_elim_gcells))
       globalID = 0
       num_elim_gcells = 0  ! local number of eliminated gridcells
       num_elim_blocks = 0  ! local number of eliminated blocks
       do jblk=1,nblocks_y
          do iblk=1,nblocks_x
             globalID = globalID + 1
             if (distrb_info%blockLocation(globalID) == 0) then
                this_block = get_block(globalID, globalID)
                num_elim_blocks = num_elim_blocks + 1
                if (num_elim_blocks >= my_elim_start .and. num_elim_blocks <= my_elim_end) then
                   do j=this_block%jlo,this_block%jhi
                      do i=this_block%ilo,this_block%ihi
                         num_elim_gcells = num_elim_gcells + 1
                         ig = this_block%i_glob(i)
                         jg = this_block%j_glob(j)
                         gindex_elim(num_elim_gcells) = (jg-1)*nx_global + ig
                      end do
                   end do
                end if
             end if
          end do
       end do

       ! create a global index that includes both active and eliminated gridcells
       num_ice  = size(gindex_ice)
       num_elim = size(gindex_elim)
       allocate(gindex(num_elim + num_ice))
       do n = 1,num_ice
          gindex(n) = gindex_ice(n)
       end do
       do n = num_ice+1,num_ice+num_elim
          gindex(n) = gindex_elim(n-num_ice)
       end do

       deallocate(gindex_elim)

    else

       ! No eliminated land blocks
       num_ice = size(gindex_ice)
       allocate(gindex(num_ice))
       do n = 1,num_ice
          gindex(n) = gindex_ice(n)
       end do

    end if

    !---------------------------------------------------------------------------
    ! Create distGrid from global index array
    !---------------------------------------------------------------------------

    DistGrid = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    deallocate(gindex_ice)
    deallocate(gindex)

  end subroutine ice_mesh_set_distgrid
  
  !=======================================================================
  subroutine ice_mesh_setmask_from_maskfile(gcomp, ice_mesh, rc)

    use ice_scam           , only : scmlat, scmlon, single_column
    use ice_grid           , only : tlon, tlat, hm, tarea, ULON, ULAT, HTN, HTE, ANGLE, ANGLET
    use ice_grid           , only : uarea, uarear, tarear, tinyarea
    use ice_grid           , only : dxt, dyt, dxu, dyu, dyhx, dxhy, cyp, cxp, cym, cxm
    use ice_grid           , only : kmt_file,  makemask, tmask
    use ice_boundary       , only : ice_HaloUpdate
    use ice_domain         , only : blocks_ice, nblocks, halo_info, distrb_info
    use ice_constants      , only : c0, c1, c2, p25, radius
    use ice_constants      , only : field_loc_center, field_type_scalar
    use ice_read_write     , only : ice_open_nc, ice_close_nc
    use netcdf

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    type(ESMF_Mesh)     , intent(in)    :: ice_mesh
    integer             , intent(out)   :: rc

    ! local variables
    character(len=char_len_long) :: ice_maskfile
    integer                      :: i,j,n
    integer                      :: iblk, jblk           ! indices
    integer                      :: ilo, ihi, jlo, jhi   ! beginning and end of physical domain
    type(ESMF_Field)             :: areaField
    real(dbl_kind) , pointer     :: mesh_areas(:)
    integer                      :: numownedelements
    real(dbl_kind) , pointer     :: ownedElemCoords(:)
    integer                      :: spatialDim
    integer (int_kind)           :: ni, nj, ncid
    integer (int_kind)           :: dimid, varid, ier
    type (block)                 :: this_block           ! block information for current block
    real (dbl_kind)              :: closelat             ! Single-column latitude value
    real (dbl_kind)              :: closelon             ! Single-column longitude value
    real (dbl_kind)              :: closelatidx          ! Single-column latitude index to retrieve
    real (dbl_kind)              :: closelonidx          ! Single-column longitude index to retrieve
    integer (int_kind)           :: start(2)             ! Start index to read in
    integer (int_kind)           :: count(2)             ! Number of points to read in
    integer (int_kind)           :: start3(3)            ! Start index to read in
    integer (int_kind)           :: count3(3)            ! Number of points to read in
    integer (int_kind)           :: status               ! status flag
    real (dbl_kind), allocatable :: lats(:)              ! temporary
    real (dbl_kind), allocatable :: lons(:)              ! temporary
    real (dbl_kind), allocatable :: pos_lons(:)          ! temporary
    real (dbl_kind), allocatable :: glob_grid(:,:)       ! temporary
    real (dbl_kind)              :: pos_scmlon           ! temporary
    real (dbl_kind)              :: scamdata             ! temporary
    integer (int_kind), pointer  :: ice_mask(:)
    real(dbl_kind)    , pointer  :: ice_frac(:)
    real(dbl_kind)               :: pi
    real(dbl_kind)               :: c180
    real(dbl_kind)               :: puny
    real(dbl_kind)               :: deg_to_rad
    logical                      :: isPresent, isSet
    character(len=char_len_long) :: cvalue
    character(len=*), parameter  :: subname = ' ice_mesh_setmask_from_maskfile'
    !---------------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine mask input file and create the mask
    call NUOPC_CompAttributeGet(gcomp, name='mesh_mask', value=ice_maskfile, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (my_task == master_task) then
       write(nu_diag,*)'mask file for cice domain is ',trim(ice_maskfile)
    end if

    ! Determine if single column
    call NUOPC_CompAttributeGet(gcomp, name='single_column', value=cvalue, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) single_column
       if (single_column) then
          call NUOPC_CompAttributeGet(gcomp, name='scmlon', value=cvalue, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          read(cvalue,*) scmlon
          call NUOPC_CompAttributeGet(gcomp, name='scmlat', value=cvalue, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          read(cvalue,*) scmlat
       end if
    else
       single_column = .false.
    end if

    ! Determine start/count to read in for either single column or global lat-lon grid
    ! If single_column, then assume that only master_task is used since there is only one task

    if (.not. single_column) then

       ! Obtain the model mask and model frac by mapping the mesh created by reading
       ! in the model_maskfile to the model mesh and then resetting the model mesh mask
       call ice_mesh_create_mask(ice_mesh, ice_maskfile, ice_mask, ice_frac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Obtain mesh areas in radians^2
       areaField = ESMF_FieldCreate(ice_mesh, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldRegridGetArea(areaField, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(areaField, farrayPtr=mesh_areas, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Obtain mesh lons and lats in degrees
       call ESMF_MeshGet(ice_mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       allocate(ownedElemCoords(spatialDim*numownedelements))
       call ESMF_MeshGet(ice_mesh, ownedElemCoords=ownedElemCoords)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_MeshGet(ice_mesh, ownedElemCoords=ownedElemCoords, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Allocate module variable ocn_gridcell_frac
       allocate(ocn_gridcell_frac(nx_block,ny_block,max_blocks))

       ! Get required constants
       call icepack_query_parameters(pi_out=pi, puny_out=puny, c180_out=c180)
       if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
            file=__FILE__, line=__LINE__)
       deg_to_rad = pi/c180

       ! Set tlon, tlat, tarea, hm and ocn_gridcell_frac
       ! Convert mesh areas from radians^2 to m^2 (tarea is in m^2)
       ! Convert lons and lats from degrees to radians
       n = 0
       do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi
          do j = jlo, jhi
             do i = ilo, ihi
                n = n + 1
                tlon(i,j,iblk) = ownedElemCoords(2*n-1) * deg_to_rad
                tlat(i,j,iblk) = ownedElemCoords(2*n) * deg_to_rad
                tarea(i,j,iblk) = mesh_areas(n) * (radius*radius)
                hm(i,j,iblk) = real(ice_mask(n),kind=dbl_kind)
                ocn_gridcell_frac(i,j,iblk) = ice_frac(n)
             enddo
          enddo
       enddo

       ! Dealocate memory
       deallocate(ownedElemCoords)
       call ESMF_FieldDestroy(areaField)

    else ! single column mode

       if (my_task == master_task) then
          call ice_open_nc(kmt_file, ncid)
          status = nf90_inq_dimid (ncid, 'ni', dimid)
          status = nf90_inquire_dimension(ncid, dimid, len=ni)
          status = nf90_inq_dimid (ncid, 'nj', dimid)
          status = nf90_inquire_dimension(ncid, dimid, len=nj)
       end if

       ! Check for consistency
       if (my_task == master_task) then
          if ((nx_global /= 1).or. (ny_global /= 1)) then
             write(nu_diag,*) 'Because you have selected the column model flag'
             write(nu_diag,*) 'Please set nx_global=ny_global=1 in file'
             write(nu_diag,*) 'ice_domain_size.F and recompile'
             call abort_ice ('ice_mesh_setmask_from_maskfile: check nx_global, ny_global')
          endif
       end if

       ! Read in domain file for single column
       allocate(lats(nj))
       allocate(lons(ni))
       allocate(pos_lons(ni))
       allocate(glob_grid(ni,nj))

       start3=(/1,1,1/)
       count3=(/ni,nj,1/)
       status = nf90_inq_varid(ncid, 'xc' , varid)
       if (status /= nf90_noerr) call abort_ice (subname//' inq_varid xc')
       status = nf90_get_var(ncid, varid, glob_grid, start3, count3)
       if (status /= nf90_noerr) call abort_ice (subname//' get_var xc')
       do i = 1,ni
          lons(i) = glob_grid(i,1)
       end do

       status = nf90_inq_varid(ncid, 'yc' , varid)
       if (status /= nf90_noerr) call abort_ice (subname//' inq_varid yc')
       status = nf90_get_var(ncid, varid, glob_grid, start3, count3)
       if (status /= nf90_noerr) call abort_ice (subname//' get_var yc')
       do j = 1,nj
          lats(j) = glob_grid(1,j)
       end do

       ! convert lons array and scmlon to 0,360 and find index of value closest to 0
       ! and obtain single-column longitude/latitude indices to retrieve

       pos_lons(:)= mod(lons(:) + 360._dbl_kind,360._dbl_kind)
       pos_scmlon = mod(scmlon  + 360._dbl_kind,360._dbl_kind)
       start(1) = (MINLOC(abs(pos_lons-pos_scmlon),dim=1))
       start(2) = (MINLOC(abs(lats    -scmlat    ),dim=1))

       deallocate(lats)
       deallocate(lons)
       deallocate(pos_lons)
       deallocate(glob_grid)

       status = nf90_inq_varid(ncid, 'xc' , varid)
       if (status /= nf90_noerr) call abort_ice (subname//' inq_varid xc')
       status = nf90_get_var(ncid, varid, scamdata, start)
       if (status /= nf90_noerr) call abort_ice (subname//' get_var xc')
       TLON = scamdata
       status = nf90_inq_varid(ncid, 'yc' , varid)
       if (status /= nf90_noerr) call abort_ice (subname//' inq_varid yc')
       status = nf90_get_var(ncid, varid, scamdata, start)
       if (status /= nf90_noerr) call abort_ice (subname//' get_var yc')
       TLAT = scamdata
       status = nf90_inq_varid(ncid, 'area' , varid)
       if (status /= nf90_noerr) call abort_ice (subname//' inq_varid area')
       status = nf90_get_var(ncid, varid, scamdata, start)
       if (status /= nf90_noerr) call abort_ice (subname//' get_var are')
       tarea = scamdata
       status = nf90_inq_varid(ncid, 'mask' , varid)
       if (status /= nf90_noerr) call abort_ice (subname//' inq_varid mask')
       status = nf90_get_var(ncid, varid, scamdata, start)
       if (status /= nf90_noerr) call abort_ice (subname//' get_var mask')
       hm = scamdata
       status = nf90_inq_varid(ncid, 'frac' , varid)
       if (status /= nf90_noerr) call abort_ice (subname//' inq_varid frac')
       status = nf90_get_var(ncid, varid, scamdata, start)
       if (status /= nf90_noerr) call abort_ice (subname//' get_var frac')
       ocn_gridcell_frac = scamdata

       if (my_task == master_task) then
          call ice_close_nc(ncid)
       end if

       do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi

          do j = jlo, jhi
             do i = ilo, ihi
                ! Convert from degrees to radians
                TLON(i,j,iblk) = pi*TLON(i,j,iblk)/180._dbl_kind

                ! Convert from degrees to radians
                TLAT(i,j,iblk) = pi*TLAT(i,j,iblk)/180._dbl_kind

                ! Convert from radians^2 to m^2
                ! (area in domain file is in radians^2 and tarea is in m^2)
                tarea(i,j,iblk) = tarea(i,j,iblk) * (radius*radius)
             end do
          end do
       end do

    end if

    call ice_HaloUpdate (TLON  , halo_info, field_loc_center, field_type_scalar, fillValue=c1)
    call ice_HaloUpdate (TLAT  , halo_info, field_loc_center, field_type_scalar, fillValue=c1)
    call ice_HaloUpdate (tarea , halo_info, field_loc_center, field_type_scalar, fillValue=c1)
    call ice_HaloUpdate (hm    , halo_info, field_loc_center, field_type_scalar, fillValue=c1)

    !-----------------------------------------------------------------
    ! CALCULATE various geometric 2d arrays
    ! The U grid (velocity) is not used when run with sequential CAM
    ! because we only use thermodynamic sea ice.  However, ULAT is used
    ! in the default initialization of CICE so we calculate it here as
    ! a "dummy" so that CICE will initialize with ice.  If a no ice
    ! initialization is OK (or desired) this can be commented out and
    ! ULAT will remain 0 as specified above.  ULAT is located at the
    ! NE corner of the grid cell, TLAT at the center, so here ULAT is
    ! hacked by adding half the latitudinal spacing (in radians) to TLAT.
    !-----------------------------------------------------------------

    ANGLET(:,:,:) = c0

    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo, jhi
          do i = ilo, ihi

             if (ny_global == 1) then
                uarea(i,j,iblk)  = tarea(i,j,  iblk)
             else
                uarea(i,j,iblk)  = p25*  &
                     (tarea(i,j,  iblk) + tarea(i+1,j,  iblk) &
                     + tarea(i,j+1,iblk) + tarea(i+1,j+1,iblk))
             endif
             tarear(i,j,iblk)   = c1/tarea(i,j,iblk)
             uarear(i,j,iblk)   = c1/uarea(i,j,iblk)
             tinyarea(i,j,iblk) = puny*tarea(i,j,iblk)

             if (single_column) then
                ULAT  (i,j,iblk) = TLAT(i,j,iblk)+(pi/nj)
             else
                if (ny_global == 1) then
                   ULAT  (i,j,iblk) = TLAT(i,j,iblk)
                else
                   ULAT  (i,j,iblk) = TLAT(i,j,iblk)+(pi/ny_global)
                endif
             endif
             ULON  (i,j,iblk) = c0
             ANGLE (i,j,iblk) = c0

             HTN   (i,j,iblk) = 1.e36_dbl_kind
             HTE   (i,j,iblk) = 1.e36_dbl_kind
             dxt   (i,j,iblk) = 1.e36_dbl_kind
             dyt   (i,j,iblk) = 1.e36_dbl_kind
             dxu   (i,j,iblk) = 1.e36_dbl_kind
             dyu   (i,j,iblk) = 1.e36_dbl_kind
             dxhy  (i,j,iblk) = 1.e36_dbl_kind
             dyhx  (i,j,iblk) = 1.e36_dbl_kind
             cyp   (i,j,iblk) = 1.e36_dbl_kind
             cxp   (i,j,iblk) = 1.e36_dbl_kind
             cym   (i,j,iblk) = 1.e36_dbl_kind
             cxm   (i,j,iblk) = 1.e36_dbl_kind
          enddo
       enddo
    enddo

    call ice_HaloUpdate (ULAT, halo_info, field_loc_center, field_type_scalar, fillValue=c1)

    ! Set the boundary values for the T cell land mask (hm) and
    ! make the logical land masks for T and U cells (tmask, umask).
    ! Also create hemisphere masks (mask-n northern, mask-s southern)
    call makemask()

  end subroutine ice_mesh_setmask_from_maskfile

  !===============================================================================
  subroutine ice_mesh_create_mask(ice_mesh, ice_maskfile, ice_mask, ice_frac, rc)

    use ice_constants, only : c0, c1

    ! input/out variables
    type(ESMF_Mesh)          , intent(in)  :: ice_mesh
    character(len=*)         , intent(in)  :: ice_maskfile
    integer        , pointer , intent(out) :: ice_mask(:)
    real(dbl_kind) , pointer , intent(out) :: ice_frac(:)
    integer                  , intent(out) :: rc

    ! local variables:
    type(ESMF_Mesh)          :: mesh_mask
    type(ESMF_Field)         :: field_mask
    type(ESMF_Field)         :: field_dst
    type(ESMF_RouteHandle)   :: rhandle
    integer                  :: srcMaskValue = 0
    integer                  :: dstMaskValue = -987987 ! spval for RH mask values
    integer                  :: srcTermProcessing_Value = 0
    logical                  :: checkflag = .false.
    real(dbl_kind) , pointer :: mask_src(:) ! on mesh created from ice_maskfile
    real(dbl_kind) , pointer :: dataptr1d(:)
    type(ESMF_DistGrid)      :: distgrid_mask
    type(ESMF_Array)         :: elemMaskArray
    integer                  :: lsize_mask, lsize_dst
    integer                  :: n, spatialDim
    real(dbl_kind)           :: fminval = 0.001_dbl_kind ! TODO: make this a share constant
    real(dbl_kind)           :: fmaxval = 1._dbl_kind
    real(dbl_kind)           :: lfrac
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    mesh_mask = ESMF_MeshCreate(trim(ice_maskfile), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MeshGet(ice_mesh, spatialDim=spatialDim, numOwnedElements=lsize_dst, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ice_mask(lsize_dst))
    allocate(ice_frac(lsize_dst))

    ! create fields on source and destination meshes
    field_mask = ESMF_FieldCreate(mesh_mask, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    field_dst = ESMF_FieldCreate(ice_mesh, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create route handle to map source mask (assume ocean) to destination mesh (assume atm/lnd)
    call ESMF_FieldRegridStore(field_mask, field_dst, routehandle=rhandle, &
         srcMaskValues=(/srcMaskValue/), dstMaskValues=(/dstMaskValue/), &
         regridmethod=ESMF_REGRIDMETHOD_CONSERVE, normType=ESMF_NORMTYPE_DSTAREA, &
         srcTermProcessing=srcTermProcessing_Value, &
         ignoreDegenerate=.true., unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! fill in values for field_mask with mask on source mesh
    call ESMF_MeshGet(mesh_mask, elementdistGrid=distgrid_mask, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_DistGridGet(distgrid_mask, localDe=0, elementCount=lsize_mask, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(mask_src(lsize_mask))
    elemMaskArray = ESMF_ArrayCreate(distgrid_mask, mask_src, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! The following call fills in the values of mask_src
    call ESMF_MeshGet(mesh_mask, elemMaskArray=elemMaskArray, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! The following call fills in the values of field_mask
    call ESMF_FieldGet(field_mask, farrayptr=dataptr1d, rc=rc)
    dataptr1d(:) = mask_src(:)

    ! map source mask to destination mesh - to obtain destination mask and frac
    call ESMF_FieldRegrid(field_mask, field_dst, routehandle=rhandle, &
         termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_TOTAL, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MeshGet(ice_mesh, spatialDim=spatialDim, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_dst, farrayptr=dataptr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! now determine ice_mask and ice_frac
    do n = 1,lsize_dst
       lfrac = c1 - dataptr1d(n)
       if (lfrac > fmaxval) lfrac = c1
       if (lfrac < fminval) lfrac = c0
       ice_frac(n) = c1 - lfrac
       if (ice_frac(n) == c0) then
          ice_mask(n) = 0
       else
          ice_mask(n) = 1
       end if
    enddo

    ! reset the model mesh mask
    call ESMF_MeshSet(ice_mesh, elementMask=ice_mask, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! deallocate memory
    call ESMF_RouteHandleDestroy(rhandle, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldDestroy(field_mask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldDestroy(field_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    deallocate(mask_src)

  end subroutine ice_mesh_create_mask

  !===============================================================================
  subroutine ice_mesh_check(gcomp, ice_mesh, rc)

    ! Check CICE mesh

    use ice_constants, only : c1,c0,c360
    use ice_grid     , only : tlon, tlat

    ! input/output parameters
    type(ESMF_GridComp) , intent(inout) :: gcomp
    type(ESMF_Mesh)     , intent(inout) :: ice_mesh
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_DistGrid)          :: distGrid
    integer                      :: n,c,g,i,j,m        ! indices
    integer                      :: iblk, jblk         ! indices
    integer                      :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    type(block)                  :: this_block         ! block information for current block
    integer                      :: spatialDim
    integer                      :: numOwnedElements
    real(dbl_kind), pointer      :: ownedElemCoords(:)
    real(dbl_kind), pointer      :: lat(:), latMesh(:)
    real(dbl_kind), pointer      :: lon(:), lonMesh(:)
    real(dbl_kind)               :: diff_lon
    real(dbl_kind)               :: diff_lat
    real(dbl_kind)               :: rad_to_deg
    real(dbl_kind)               :: tmplon, eps_imesh
    logical                      :: isPresent, isSet
    character(len=char_len_long) :: cvalue
    character(len=char_len_long) :: logmsg
    character(len=*), parameter  :: subname = ' ice_mesh_check: '
    !---------------------------------------------------

    ! Determine allowed mesh error
    call NUOPC_CompAttributeGet(gcomp, name='eps_imesh', value=cvalue, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) eps_imesh
    else
       eps_imesh = 1.0e-1_dbl_kind
    end if
    write(logmsg,*) eps_imesh
    call ESMF_LogWrite(trim(subname)//' eps_imesh = '//trim(logmsg), ESMF_LOGMSG_INFO)

    ! error check differences between internally generated lons and those read in
    call ESMF_MeshGet(ice_mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numownedelements))
    allocate(lonmesh(numOwnedElements))
    allocate(latmesh(numOwnedElements))
    call ESMF_MeshGet(ice_mesh, ownedElemCoords=ownedElemCoords)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 1,numOwnedElements
       lonMesh(n) = ownedElemCoords(2*n-1)
       latMesh(n) = ownedElemCoords(2*n)
    end do

    ! obtain internally generated cice lats and lons for error checks
    call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
    allocate(lon(numOwnedElements))
    allocate(lat(numOwnedElements))
    lon(:) = 0.
    lat(:) = 0.
    n = 0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
          do i = ilo, ihi
             n = n + 1
             lon(n) = tlon(i,j,iblk)*rad_to_deg
             lat(n) = tlat(i,j,iblk)*rad_to_deg

             tmplon = lon(n)
             if(tmplon < c0)tmplon = tmplon + c360

             ! error check differences between internally generated lons and those read in
             diff_lon = abs(mod(lonMesh(n) - tmplon,360.0))
             if (diff_lon > eps_imesh ) then
                write(6,100)n,lonMesh(n),tmplon, diff_lon
                call abort_ice(error_message=subname, &
                     file=__FILE__, line=__LINE__)
             end if
             diff_lat = abs(latMesh(n) - lat(n))
             if (diff_lat > eps_imesh) then
                write(6,101)n,latMesh(n),lat(n), diff_lat
                call abort_ice(error_message=subname, &
                     file=__FILE__, line=__LINE__)
             end if

          enddo
       enddo
    enddo

100 format('ERROR: CICE n, lonmesh, lon, diff_lon = ',i6,2(f21.13,3x),d21.5)
101 format('ERROR: CICE n, latmesh, lat, diff_lat = ',i6,2(f21.13,3x),d21.5)

    ! deallocate memory
    deallocate(ownedElemCoords)
    deallocate(lon, lonMesh)
    deallocate(lat, latMesh)

  end subroutine ice_mesh_check

end module ice_mesh_mod
