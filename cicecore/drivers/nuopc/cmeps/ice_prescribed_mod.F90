module ice_prescribed_mod

  ! !DESCRIPTION:
  ! The prescribed ice model reads in ice concentration data from a netCDF
  ! file.  Ice thickness, temperature, the ice temperature profile are
  ! prescribed.  Air/ice fluxes are computed to get surface temperature,
  ! Ice/ocean fluxes are set to zero, and ice dynamics are not calculated.
  ! Regridding and data cycling capabilities are included.

#ifndef CESMCOUPLED

  use ice_kinds_mod

  implicit none
  private ! except

  public  :: ice_prescribed_init      ! initialize input data stream
  logical(kind=log_kind), parameter, public :: prescribed_ice = .false.     ! true if prescribed ice

contains
  ! This is a stub routine for now
  subroutine ice_prescribed_init(mpicom, compid, gindex)
    integer(kind=int_kind), intent(in) :: mpicom
    integer(kind=int_kind), intent(in) :: compid
    integer(kind=int_kind), intent(in) :: gindex(:)
    ! do nothing
  end subroutine ice_prescribed_init

#else 

  use shr_nl_mod        , only : shr_nl_find_group_name
  use shr_strdata_mod
  use shr_dmodel_mod
  use shr_string_mod
  use shr_ncread_mod
  use shr_sys_mod
  use shr_mct_mod
  use mct_mod
  use pio
  use ice_broadcast
  use ice_communicate   , only : my_task, master_task, MPI_COMM_ICE
  use ice_kinds_mod
  use ice_fileunits
  use ice_exit          , only : abort_ice
  use ice_domain_size   , only : nx_global, ny_global, ncat, nilyr, nslyr, max_blocks
  use ice_constants
  use ice_blocks        , only : nx_block, ny_block, block, get_block
  use ice_domain        , only : nblocks, distrb_info, blocks_ice
  use ice_grid          , only : TLAT, TLON, hm, tmask, tarea, grid_type, ocn_gridcell_frac
  use ice_calendar      , only : idate, calendar_type
  use ice_arrays_column , only : hin_max
  use ice_read_write
  use ice_exit          , only: abort_ice
  use icepack_intfc     , only: icepack_warnings_flush, icepack_warnings_aborted
  use icepack_intfc     , only: icepack_query_tracer_indices, icepack_query_tracer_sizes
  use icepack_intfc     , only: icepack_query_parameters

  implicit none
  private ! except

  ! MEMBER FUNCTIONS:
  public  :: ice_prescribed_init      ! initialize input data stream
  public  :: ice_prescribed_run       ! get time slices and time interp
  public  :: ice_prescribed_phys      ! set prescribed ice state and fluxes

  ! !PUBLIC DATA MEMBERS:
  logical(kind=log_kind), public   :: prescribed_ice      ! true if prescribed ice
  integer(kind=int_kind),parameter :: nFilesMaximum = 400 ! max number of files
  integer(kind=int_kind)           :: stream_year_first   ! first year in stream to use
  integer(kind=int_kind)           :: stream_year_last    ! last year in stream to use
  integer(kind=int_kind)           :: model_year_align    ! align stream_year_first with this model year
  character(len=char_len_long)     :: stream_fldVarName
  character(len=char_len_long)     :: stream_fldFileName(nFilesMaximum)
  character(len=char_len_long)     :: stream_domTvarName
  character(len=char_len_long)     :: stream_domXvarName
  character(len=char_len_long)     :: stream_domYvarName
  character(len=char_len_long)     :: stream_domAreaName
  character(len=char_len_long)     :: stream_domMaskName
  character(len=char_len_long)     :: stream_domFileName
  character(len=char_len_long)     :: stream_mapread
  logical(kind=log_kind)           :: prescribed_ice_fill ! true if data fill required
  type(shr_strdata_type)           :: sdat                ! prescribed data stream
  character(len=char_len_long)     :: fldList             ! list of fields in data stream
  real(kind=dbl_kind),allocatable  :: ice_cov(:,:,:)      ! ice cover

contains

  subroutine ice_prescribed_init(mpicom, compid, gindex)

    !    Prescribed ice initialization - needed to
    !    work with new shr_strdata module derived type

    use shr_pio_mod, only : shr_pio_getiotype, shr_pio_getiosys, shr_pio_getioformat

    implicit none
    include 'mpif.h'

    ! !nput/output parameters:
    integer(kind=int_kind), intent(in) :: mpicom
    integer(kind=int_kind), intent(in) :: compid
    integer(kind=int_kind), intent(in) :: gindex(:)

    !----- Local ------
    type(mct_gsMap)        :: gsmap_ice
    type(mct_gGrid)        :: dom_ice
    integer(kind=int_kind) :: lsize
    integer(kind=int_kind) :: gsize
    integer(kind=int_kind) :: nml_error ! namelist i/o error flag
    integer(kind=int_kind) :: n, nFile, ierr
    character(len=8)       :: fillalgo
    character(*),parameter :: subName = '(ice_prescribed_init)'

    namelist /ice_prescribed_nml/  &
         prescribed_ice,      &
         model_year_align,    &
         stream_year_first ,  &
         stream_year_last  ,  &
         stream_fldVarName ,  &
         stream_fldFileName,  &
         stream_domTvarName,  &
         stream_domXvarName,  &
         stream_domYvarName,  &
         stream_domAreaName,  &
         stream_domMaskName,  &
         stream_domFileName,  &
         stream_mapread,      &
         prescribed_ice_fill

    ! default values for namelist
    prescribed_ice         = .false.          ! if true, prescribe ice
    stream_year_first      = 1                ! first year in  pice stream to use
    stream_year_last       = 1                ! last  year in  pice stream to use
    model_year_align       = 1                ! align stream_year_first with this model year
    stream_fldVarName      = 'ice_cov'
    stream_fldFileName(:)  = ' '
    stream_domTvarName     = 'time'
    stream_domXvarName     = 'lon'
    stream_domYvarName     = 'lat'
    stream_domAreaName     = 'area'
    stream_domMaskName     = 'mask'
    stream_domFileName     = ' '
    stream_mapread         = 'NOT_SET'
    prescribed_ice_fill    = .false.          ! true if pice data fill required

    ! read from input file
    call get_fileunit(nu_nml)
    if (my_task == master_task) then
       open (nu_nml, file=nml_filename, status='old',iostat=nml_error)
       call shr_nl_find_group_name(nu_nml, 'ice_prescribed_nml', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, ice_prescribed_nml, iostat=nml_error)
          if (nml_error > 0) then
             call shr_sys_abort( 'problem on read of ice_prescribed namelist in ice_prescribed_mod' )
          endif
       endif
    end if
    call release_fileunit(nu_nml)
    call broadcast_scalar(prescribed_ice, master_task)

    ! *** If not prescribed ice then return ***
    if (.not. prescribed_ice) RETURN

    call broadcast_scalar(model_year_align,master_task)
    call broadcast_scalar(stream_year_first,master_task)
    call broadcast_scalar(stream_year_last,master_task)
    call broadcast_scalar(stream_fldVarName,master_task)
    call broadcast_scalar(stream_domTvarName,master_task)
    call broadcast_scalar(stream_domXvarName,master_task)
    call broadcast_scalar(stream_domYvarName,master_task)
    call broadcast_scalar(stream_domAreaName,master_task)
    call broadcast_scalar(stream_domMaskName,master_task)
    call broadcast_scalar(stream_domFileName,master_task)
    call broadcast_scalar(stream_mapread,master_task)
    call broadcast_scalar(prescribed_ice_fill,master_task)
    call mpi_bcast(stream_fldFileName, len(stream_fldFileName(1))*NFilesMaximum, &
         MPI_CHARACTER, 0, MPI_COMM_ICE, ierr)

    nFile = 0
    do n=1,nFilesMaximum
       if (stream_fldFileName(n) /= ' ') nFile = nFile + 1
    end do

    ! Read shr_strdata_nml namelist
    if (prescribed_ice_fill) then
       fillalgo='nn'
    else
       fillalgo='none'
    endif

    if (my_task == master_task) then
       write(nu_diag,*) ' '
       write(nu_diag,*) 'This is the prescribed ice coverage option.'
       write(nu_diag,*) '  stream_year_first  = ',stream_year_first
       write(nu_diag,*) '  stream_year_last   = ',stream_year_last
       write(nu_diag,*) '  model_year_align   = ',model_year_align
       write(nu_diag,*) '  stream_fldVarName  = ',trim(stream_fldVarName)
       do n = 1,nFile
          write(nu_diag,*) '  stream_fldFileName = ',trim(stream_fldFileName(n)),n
       end do
       write(nu_diag,*) '  stream_domTvarName = ',trim(stream_domTvarName)
       write(nu_diag,*) '  stream_domXvarName = ',trim(stream_domXvarName)
       write(nu_diag,*) '  stream_domYvarName = ',trim(stream_domYvarName)
       write(nu_diag,*) '  stream_domFileName = ',trim(stream_domFileName)
       write(nu_diag,*) '  stream_mapread     = ',trim(stream_mapread)
       write(nu_diag,*) '  stream_fillalgo    = ',trim(fillalgo)
       write(nu_diag,*) ' '
    endif

    gsize = nx_global*ny_global
    lsize = size(gindex)
    call mct_gsMap_init( gsmap_ice, gindex, MPI_COMM_ICE, compid, lsize, gsize)
    call ice_prescribed_set_domain( lsize, MPI_COMM_ICE, gsmap_ice, dom_ice )

    call shr_strdata_create(sdat,name="prescribed_ice", &
         mpicom=MPI_COMM_ICE, compid=compid,   &
         gsmap=gsmap_ice, ggrid=dom_ice,       &
         nxg=nx_global,nyg=ny_global,          &
         yearFirst=stream_year_first,          &
         yearLast=stream_year_last,            &
         yearAlign=model_year_align,           &
         offset=0,                             &
         domFilePath='',                       &
         domFileName=trim(stream_domFileName), &
         domTvarName=stream_domTvarName,       &
         domXvarName=stream_domXvarName,       &
         domYvarName=stream_domYvarName,       &
         domAreaName=stream_domAreaName,       &
         domMaskName=stream_domMaskName,       &
         filePath='',                          &
         filename=stream_fldFileName(1:nFile), &
         fldListFile=stream_fldVarName,        &
         fldListModel=stream_fldVarName,       &
         fillalgo=trim(fillalgo),              &
         calendar=trim(calendar_type),         &
         mapread=trim(stream_mapread))

    if (my_task == master_task) then
       call shr_strdata_print(sdat,'SPRESICE data')
    endif

    !-----------------------------------------------------------------
    ! For one ice category, set hin_max(1) to something big
    !-----------------------------------------------------------------
    if (ncat == 1) then
       hin_max(1) = 999._dbl_kind
    end if
  end subroutine ice_prescribed_init

  !=======================================================================
  subroutine ice_prescribed_run(mDateIn, secIn)

    ! !DESCRIPTION:
    !  Finds two time slices bounding current model time, remaps if necessary

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    integer(kind=int_kind), intent(in) :: mDateIn  ! Current model date (yyyymmdd)
    integer(kind=int_kind), intent(in) :: secIn    ! Elapsed seconds on model date

    ! local varaibles
    integer(kind=int_kind) :: i,j,n,iblk       ! loop indices and counter
    integer(kind=int_kind) :: ilo,ihi,jlo,jhi  ! beginning and end of physical domain
    type (block)           :: this_block
    real(kind=dbl_kind)    :: aice_max         ! maximun ice concentration
    logical, save          :: first_time = .true.
    character(*),parameter :: subName = '(ice_prescribed_run)'
    character(*),parameter :: F00 = "(a,2g20.13)"

    !------------------------------------------------------------------------
    ! Interpolate to new ice coverage
    !------------------------------------------------------------------------

    call shr_strdata_advance(sdat,mDateIn,SecIn,MPI_COMM_ICE,'cice_pice')

    if (first_time) then
       allocate(ice_cov(nx_block,ny_block,max_blocks))
    endif

    ice_cov(:,:,:) = c0  ! This initializes ghost cells as well

    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo, jhi
          do i = ilo, ihi
             n = n+1
             ice_cov(i,j,iblk) = sdat%avs(1)%rAttr(1,n)
          end do
       end do
    end do

    !--------------------------------------------------------------------
    ! Check to see that ice concentration is in fraction, not percent
    !--------------------------------------------------------------------
    if (first_time) then
       aice_max = maxval(ice_cov)

       if (aice_max > c10) then
          write(nu_diag,F00) subname//" ERROR: Ice conc data must be in fraction, aice_max= ",&
               aice_max
          call abort_ice(subName)
       end if
       first_time = .false.
    end if

    !-----------------------------------------------------------------
    ! Set prescribed ice state and fluxes
    !-----------------------------------------------------------------

    call ice_prescribed_phys()

  end subroutine ice_prescribed_run

  !===============================================================================
  subroutine ice_prescribed_phys

    ! Set prescribed ice state using input ice concentration;
    ! set surface ice temperature to atmospheric value; use
    ! linear temperature gradient in ice to ocean temperature.

    ! !USES:
    use ice_flux
    use ice_state
    use icepack_intfc, only : icepack_aggregate
    use ice_dyn_evp
    implicit none

    !----- Local ------
    integer(kind=int_kind) :: layer    ! level index
    integer(kind=int_kind) :: nc       ! ice category index
    integer(kind=int_kind) :: i,j,k    ! longitude, latitude and level indices
    integer(kind=int_kind) :: iblk
    integer(kind=int_kind) :: nt_Tsfc, nt_sice, nt_qice, nt_qsno, ntrcr

    real(kind=dbl_kind) :: slope     ! diff in underlying ocean tmp and ice surface tmp
    real(kind=dbl_kind) :: Ti        ! ice level temperature
    real(kind=dbl_kind) :: Tmlt      ! ice level melt temperature
    real(kind=dbl_kind) :: qin_save(nilyr)
    real(kind=dbl_kind) :: qsn_save(nslyr)
    real(kind=dbl_kind) :: hi        ! ice prescribed (hemispheric) ice thickness
    real(kind=dbl_kind) :: hs        ! snow thickness
    real(kind=dbl_kind) :: zn        ! normalized ice thickness
    real(kind=dbl_kind) :: salin(nilyr)  ! salinity (ppt)
    real(kind=dbl_kind) :: rad_to_deg, pi, puny
    real(kind=dbl_kind) :: rhoi, rhos, cp_ice, cp_ocn, lfresh, depressT

    real(kind=dbl_kind), parameter :: nsal    = 0.407_dbl_kind
    real(kind=dbl_kind), parameter :: msal    = 0.573_dbl_kind
    real(kind=dbl_kind), parameter :: saltmax = 3.2_dbl_kind   ! max salinity at ice base (ppm)
    character(*),parameter :: subName = '(ice_prescribed_phys)'

    call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_sice_out=nt_sice, &
       nt_qice_out=nt_qice, nt_qsno_out=nt_qsno)
    call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
    call icepack_query_parameters(rad_to_deg_out=rad_to_deg, pi_out=pi, &
       puny_out=puny, rhoi_out=rhoi, rhos_out=rhos, cp_ice_out=cp_ice, cp_ocn_out=cp_ocn, &
       lfresh_out=lfresh, depressT_out=depressT)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
       file=__FILE__, line=__LINE__)

    !-----------------------------------------------------------------
    ! Set ice cover over land to zero, not sure if this should be
    ! be done earier, before time/spatial interp??????
    !-----------------------------------------------------------------
    do iblk = 1,nblocks
       do j = 1,ny_block
          do i = 1,nx_block
             if (tmask(i,j,iblk)) then
                if (ice_cov(i,j,iblk) .lt. eps04) ice_cov(i,j,iblk) = c0
                if (ice_cov(i,j,iblk) .gt. c1)    ice_cov(i,j,iblk) = c1
             else
                ice_cov(i,j,iblk) = c0
             end if
          enddo
       enddo
    enddo

    do iblk = 1,nblocks
       do j = 1,ny_block
          do i = 1,nx_block

             if (tmask(i,j,iblk)) then   ! Over ocean points

                !--------------------------------------------------------------
                ! Place ice where ice concentration > .0001
                !--------------------------------------------------------------

                if (ice_cov(i,j,iblk) >= eps04) then

                   hi = 0.0_dbl_kind
                   !----------------------------------------------------------
                   ! Set ice thickness in each hemisphere
                   !----------------------------------------------------------
                   if(TLAT(i,j,iblk)*rad_to_deg > 40.0_dbl_kind) then
                      hi  = 2.0_dbl_kind
                   else if(TLAT(i,j,iblk)*rad_to_deg < -40.0_dbl_kind) then
                      hi  = 1.0_dbl_kind
                   end if

                   !----------------------------------------------------------
                   ! All ice in appropriate thickness category
                   !----------------------------------------------------------
                   do nc = 1,ncat

                      if(hin_max(nc-1) < hi .and. hi < hin_max(nc)) then

                         if (aicen(i,j,nc,iblk) > c0) then
                            hs = vsnon(i,j,nc,iblk) / aicen(i,j,nc,iblk)
                         else
                            hs = c0
                         endif

                         aicen(i,j,nc,iblk) = ice_cov(i,j,iblk)
                         vicen(i,j,nc,iblk) = hi*aicen(i,j,nc,iblk)
                         vsnon(i,j,nc,iblk) = hs*aicen(i,j,nc,iblk)

                         !---------------------------------------------------------
                         ! make linear temp profile and compute enthalpy
                         !---------------------------------------------------------

                         if (abs(trcrn(i,j,nt_qice,nc,iblk)) < puny) then

                            if (aice(i,j,iblk) < puny) &
                                 trcrn(i,j,nt_Tsfc,nc,iblk) = Tf(i,j,iblk)

                            slope = Tf(i,j,iblk) - trcrn(i,j,nt_Tsfc,nc,iblk)
                            do k = 1, nilyr
                               zn = (real(k,kind=dbl_kind)-p5) / real(nilyr,kind=dbl_kind)
                               Ti = trcrn(i,j,nt_Tsfc,nc,iblk) + slope*zn
                               salin(k) = (saltmax/c2)*(c1-cos(pi*zn**(nsal/(msal+zn))))
                               Tmlt = -salin(k)*depressT
                               trcrn(i,j,nt_sice+k-1,nc,iblk) = salin(k)
                               trcrn(i,j,nt_qice+k-1,nc,iblk) =                      &
                                    -(rhoi * (cp_ice*(Tmlt-Ti) &
                                    + Lfresh*(c1-Tmlt/Ti) - cp_ocn*Tmlt))
                            enddo

                            do k=1,nslyr
                               trcrn(i,j,nt_qsno+k-1,nc,iblk) =                      &
                                    -rhos*(Lfresh - cp_ice*trcrn(i,j,nt_Tsfc,nc,iblk))
                            enddo

                         endif  ! aice < puny
                      end if    ! hin_max
                   enddo        ! ncat
                else
                   trcrn(i,j,nt_Tsfc,:,iblk) = Tf(i,j,iblk)
                   aicen(i,j,:,iblk) = c0
                   vicen(i,j,:,iblk) = c0
                   vsnon(i,j,:,iblk) = c0
                   trcrn(i,j,nt_sice:nt_sice+nilyr-1,:,iblk) = c0
                   trcrn(i,j,nt_qice:nt_qice+nilyr-1,:,iblk) = c0
                   trcrn(i,j,nt_qsno:nt_qsno+nslyr-1,:,iblk) = c0
                end if          ! ice_cov >= eps04

                !--------------------------------------------------------------------
                ! compute aggregate ice state and open water area
                !--------------------------------------------------------------------
                call icepack_aggregate(ncat  = ncat,                  &
                                       aicen = aicen(i,j,:,iblk),     &
                                       trcrn = trcrn(i,j,1:ntrcr,:,iblk), &
                                       vicen = vicen(i,j,:,iblk),     &
                                       vsnon = vsnon(i,j,:,iblk),     &
                                       aice  = aice (i,j,  iblk),     &
                                       trcr  = trcr (i,j,1:ntrcr,iblk), &
                                       vice  = vice (i,j,  iblk),     &
                                       vsno  = vsno (i,j,  iblk),     &
                                       aice0 = aice0(i,j,  iblk),     &
                                       ntrcr = ntrcr,                 &
                                       trcr_depend   = trcr_depend(1:ntrcr),   &
                                       trcr_base     = trcr_base(1:ntrcr,:),   &
                                       n_trcr_strata = n_trcr_strata(1:ntrcr), &
                                       nt_strata     = nt_strata(1:ntrcr,:))
             end if             ! tmask
          enddo                 ! i
       enddo                 ! j
    enddo                 ! iblk

    do iblk = 1, nblocks
       do j = 1, ny_block
          do i = 1, nx_block
             aice_init(i,j,iblk) = aice(i,j,iblk)
          enddo
       enddo
    enddo

    !--------------------------------------------------------------------
    ! set non-computed fluxes, ice velocities, ice-ocn stresses to zero
    !--------------------------------------------------------------------

    frzmlt    (:,:,:) = c0
    uvel      (:,:,:) = c0
    vvel      (:,:,:) = c0
    strocnxT  (:,:,:) = c0
    strocnyT  (:,:,:) = c0

    !-----------------------------------------------------------------
    ! other atm and ocn fluxes
    !-----------------------------------------------------------------
    call init_flux_atm
    call init_flux_ocn

  end subroutine ice_prescribed_phys

  !===============================================================================
  subroutine ice_prescribed_set_domain( lsize, mpicom, gsmap_i, dom_i )

    ! Arguments
    integer        , intent(in)    :: lsize
    integer        , intent(in)    :: mpicom
    type(mct_gsMap), intent(in)    :: gsMap_i
    type(mct_ggrid), intent(inout) :: dom_i

    ! Local Variables
    integer                 :: i, j, iblk, n      ! indices
    integer                 :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    real(dbl_kind), pointer :: data1(:)           ! temporary
    real(dbl_kind), pointer :: data2(:)           ! temporary
    real(dbl_kind), pointer :: data3(:)           ! temporary
    real(dbl_kind), pointer :: data4(:)           ! temporary
    real(dbl_kind), pointer :: data5(:)           ! temporary
    real(dbl_kind), pointer :: data6(:)           ! temporary
    integer       , pointer :: idata(:)           ! temporary
    real(kind=dbl_kind)     :: rad_to_deg
    type(block)             :: this_block         ! block information for current block
    character(*),parameter :: subName = '(ice_prescribed_set_domain)'
    !--------------------------------

    call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
       file=__FILE__, line=__LINE__)

    ! Initialize mct domain type
    call mct_gGrid_init(GGrid=dom_i, &
         CoordChars='lat:lon:hgt', OtherChars='area:aream:mask:frac', lsize=lsize )
    call mct_aVect_zero(dom_i%data)

    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    call mct_gsMap_orderedPoints(gsMap_i, my_task, idata)
    call mct_gGrid_importIAttr(dom_i,'GlobGridNum',idata,lsize)
    deallocate(idata)

    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value

    allocate(data1(lsize))
    allocate(data2(lsize))
    allocate(data3(lsize))
    allocate(data4(lsize))
    allocate(data5(lsize))
    allocate(data6(lsize))

    data1(:) = -9999.0_dbl_kind
    data2(:) = -9999.0_dbl_kind
    data3(:) = -9999.0_dbl_kind
    data4(:) = -9999.0_dbl_kind
    call mct_gGrid_importRAttr(dom_i,"lat"  ,data1,lsize)
    call mct_gGrid_importRAttr(dom_i,"lon"  ,data2,lsize)
    call mct_gGrid_importRAttr(dom_i,"area" ,data3,lsize)
    call mct_gGrid_importRAttr(dom_i,"aream",data4,lsize)
    data5(:) = 0.0_dbl_kind
    data6(:) = 0.0_dbl_kind
    call mct_gGrid_importRAttr(dom_i,"mask" ,data5,lsize)
    call mct_gGrid_importRAttr(dom_i,"frac" ,data6,lsize)

    ! Fill in correct values for domain components
    ! lat/lon in degrees,  area in radians^2, mask is 1 (ocean), 0 (non-ocean)
    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
       do i = ilo, ihi
          n = n+1

          data1(n) = TLON(i,j,iblk)*rad_to_deg
          data2(n) = TLAT(i,j,iblk)*rad_to_deg
          data3(n) = tarea(i,j,iblk)/(radius*radius)

          data5(n) = real(nint(hm(i,j,iblk)),kind=dbl_kind)
          if (trim(grid_type) == 'latlon') then
             data6(n) = ocn_gridcell_frac(i,j,iblk)
          else
             data6(n) = real(nint(hm(i,j,iblk)),kind=dbl_kind)
          end if

       enddo   !i
       enddo   !j
    enddo      !iblk
    call mct_gGrid_importRattr(dom_i,"lon" ,data1,lsize)
    call mct_gGrid_importRattr(dom_i,"lat" ,data2,lsize)
    call mct_gGrid_importRattr(dom_i,"area",data3,lsize)
    call mct_gGrid_importRattr(dom_i,"mask",data5,lsize)
    call mct_gGrid_importRattr(dom_i,"frac",data6,lsize)

    deallocate(data1, data2, data3, data4, data5, data6)

  end subroutine ice_prescribed_set_domain

#endif

end module ice_prescribed_mod
