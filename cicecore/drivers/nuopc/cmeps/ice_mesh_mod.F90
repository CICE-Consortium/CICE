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
  public  :: ice_mesh_create_scolumn
  public  :: ice_mesh_init_tlon_tlat_area_hm
  public  :: ice_mesh_check

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
  subroutine ice_mesh_setmask_from_maskfile(ice_maskfile, ice_mesh, rc)

    use ice_grid      , only : tlon, tlat, hm, tarea
    use ice_constants , only : c0, c1, c2, p25, radius

    ! input/output variables
    character(len=*) , intent(in)    :: ice_maskfile
    type(ESMF_Mesh)  , intent(inout) :: ice_mesh
    integer          , intent(out)   :: rc

    ! local variables
    integer                     :: i, j, n
    integer (int_kind)          :: ni, nj
    integer                     :: iblk, jblk           ! indices
    integer                     :: ilo, ihi, jlo, jhi   ! beginning and end of physical domain
    type (block)                :: this_block           ! block information for current block
    real(dbl_kind) , pointer    :: ice_frac(:)
    type(ESMF_Field)            :: areaField
    type(ESMF_Mesh)             :: mesh_mask
    type(ESMF_Field)            :: field_mask
    type(ESMF_Field)            :: field_dst
    type(ESMF_RouteHandle)      :: rhandle
    integer                     :: srcMaskValue = 0
    integer                     :: dstMaskValue = -987987 ! spval for RH mask values
    integer                     :: srcTermProcessing_Value = 0
    logical                     :: checkflag = .false.
    integer, pointer            :: ice_mask(:)
    real(dbl_kind) , pointer    :: mask_src(:) ! on mesh created from ice_maskfile
    real(dbl_kind) , pointer    :: dataptr1d(:)
    type(ESMF_DistGrid)         :: distgrid_mask
    type(ESMF_Array)            :: elemMaskArray
    integer                     :: lsize_mask, lsize_dst
    integer                     :: spatialDim
    real(dbl_kind)              :: fminval = 0.001_dbl_kind ! TODO: make this a share constant
    real(dbl_kind)              :: fmaxval = 1._dbl_kind
    real(dbl_kind)              :: lfrac
    real(dbl_kind) , pointer    :: mesh_areas(:)
    integer                     :: numownedelements
    real(dbl_kind) , pointer    :: ownedElemCoords(:)
    real(dbl_kind)              :: pi
    real(dbl_kind)              :: c180
    real(dbl_kind)              :: puny
    real(dbl_kind)              :: deg_to_rad
    character(len=*), parameter :: subname = ' ice_mesh_setmask_from_maskfile'
    !---------------------------------------------------

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
    do n = 1,size(dataptr1d)
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

    ! Allocate module variable ocn_gridcell_frac
    allocate(ocn_gridcell_frac(nx_block,ny_block,max_blocks))

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

    ! Get required constants
    call icepack_query_parameters(pi_out=pi, c180_out=c180)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)
    deg_to_rad = pi/c180

    ! Set tlon, tlat, tarea, hm
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

  end subroutine ice_mesh_setmask_from_maskfile

  !===============================================================================
  subroutine ice_mesh_create_scolumn(scol_lon, scol_lat, ice_mesh, rc)

    use ice_constants , only : c0, c1
    use ice_scam      , only : scmlat, scmlon, scol_area, scol_mask, scol_frac, scol_nj
    use netcdf

    ! Create the model mesh from the domain file - for either single column mode
    ! or for a regional grid

    ! input/output variables
    real(dbl_kind)  , intent(in)    :: scol_lon
    real(dbl_kind)  , intent(in)    :: scol_lat
    type(ESMF_Mesh) , intent(inout) :: ice_mesh
    integer         , intent(out)   :: rc

    ! local variables
    type(ESMF_Grid) :: lgrid
    integer         :: maxIndex(2)
    real(dbl_kind)  :: mincornerCoord(2)
    real(dbl_kind)  :: maxcornerCoord(2)
    integer         :: i, j,iblk, jblk      ! indices
    integer         :: ilo, ihi, jlo, jhi   ! beginning and end of physical domain
    type (block)    :: this_block           ! block information for current block
    character(len=*), parameter  :: subname = ' ice_mesh_create_scolumn'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Use center and come up with arbitrary area delta lon and lat = .1 degree
    maxIndex(1)       = 1                ! number of lons
    maxIndex(2)       = 1                ! number of lats
    mincornerCoord(1) = scol_lon - .1_dbl_kind ! min lon
    mincornerCoord(2) = scol_lat - .1_dbl_kind ! min lat
    maxcornerCoord(1) = scol_lon + .1_dbl_kind ! max lon
    maxcornerCoord(2) = scol_lat + .1_dbl_kind ! max lat

    ! create the ESMF grid
    lgrid = ESMF_GridCreateNoPeriDimUfrm (maxindex=maxindex, &
         mincornercoord=mincornercoord, maxcornercoord= maxcornercoord, &
         staggerloclist=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create the mesh from the lgrid
    ice_mesh = ESMF_MeshCreate(lgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Allocate module variable ocn_gridcell_frac
    allocate(ocn_gridcell_frac(nx_block,ny_block,max_blocks))
    ocn_gridcell_frac(:,:,:) = scol_frac
        
  end subroutine ice_mesh_create_scolumn

  !===============================================================================
  subroutine ice_mesh_init_tlon_tlat_area_hm()

    use ice_grid      , only : tlon, tlat, hm, tarea, ULON, ULAT, HTN, HTE, ANGLE, ANGLET
    use ice_grid      , only : uarea, uarear, tarear, tinyarea
    use ice_grid      , only : dxt, dyt, dxu, dyu, dyhx, dxhy, cyp, cxp, cym, cxm
    use ice_grid      , only : makemask
    use ice_boundary  , only : ice_HaloUpdate
    use ice_domain    , only : blocks_ice, nblocks, halo_info, distrb_info
    use ice_constants , only : c0, c1, p25
    use ice_constants , only : field_loc_center, field_type_scalar
    use ice_scam      , only : scmlat, scmlon, scol_area, scol_mask, scol_frac, scol_nj, single_column

    ! local variables
    integer        :: i,j,n
    integer        :: iblk, jblk           ! indices
    integer        :: ilo, ihi, jlo, jhi   ! beginning and end of physical domain
    type (block)   :: this_block           ! block information for current block
    real(dbl_kind) :: puny
    real(dbl_kind) :: pi
    character(len=*), parameter  :: subname = ' ice_mesh_init_tlon_tlat_area_hm'
    ! ----------------------------------------------

    ! Get required constants
    call icepack_query_parameters(pi_out=pi, puny_out=puny)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

    ! Check for consistency
    if (single_column) then
       if ((nx_global /= 1).or. (ny_global /= 1)) then
          write(nu_diag,*) 'nx_global = ',nx_global
          write(nu_diag,*) 'ny_global = ',ny_global
          write(nu_diag,*) 'Because you have selected the column model flag'
          write(nu_diag,*) 'then require nx_global=ny_global=1 in file ice_domain_size.F'
          call abort_ice(' ice_mesh_init_tlon_tlat_area_hm: nx_global and ny_global need to be 1 for single column')
       else
          write(nu_diag,'(a,f10.5)')' single column mode lon/lat does contain ocn with ocn fraction ',scol_frac
       end if

       TLON  = scmlon
       TLAT  = scmlat
       tarea = scol_area
       hm    = scol_mask
       ULAT  = TLAT + pi/scol_nj
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

             if (.not. single_column) then
                if (ny_global == 1) then
                   ULAT(i,j,iblk) = TLAT(i,j,iblk)
                else
                   ULAT(i,j,iblk) = TLAT(i,j,iblk)+(pi/ny_global)
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

  end subroutine ice_mesh_init_tlon_tlat_area_hm

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
                !call abort_ice(error_message=subname, file=__FILE__, line=__LINE__)
             end if
             diff_lat = abs(latMesh(n) - lat(n))
             if (diff_lat > eps_imesh) then
                write(6,101)n,latMesh(n),lat(n), diff_lat
                !call abort_ice(error_message=subname, file=__FILE__, line=__LINE__)
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
