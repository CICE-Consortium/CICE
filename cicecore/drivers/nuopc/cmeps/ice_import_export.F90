module ice_import_export

  use ESMF
  use NUOPC
  use NUOPC_Model
  use ice_kinds_mod      , only : int_kind, dbl_kind, char_len, char_len_long, log_kind
  use ice_constants      , only : c0, c1, spval_dbl, radius
  use ice_constants      , only : field_loc_center, field_type_scalar, field_type_vector
  use ice_blocks         , only : block, get_block, nx_block, ny_block
  use ice_domain         , only : nblocks, blocks_ice, halo_info, distrb_info
  use ice_domain_size    , only : nx_global, ny_global, block_size_x, block_size_y, max_blocks, ncat
  use ice_domain_size    , only : nfreq, nfsd
  use ice_exit           , only : abort_ice
  use ice_flux           , only : strairxT, strairyT, strocnxT_iavg, strocnyT_iavg
  use ice_flux           , only : alvdr, alidr, alvdf, alidf, Tref, Qref, Uref
  use ice_flux           , only : flat, fsens, flwout, evap, fswabs, fhocn, fswthru
  use ice_flux           , only : fswthru_vdr, fswthru_vdf, fswthru_idr, fswthru_idf
  use ice_flux           , only : send_i2x_per_cat, fswthrun_ai
  use ice_flux_bgc       , only : faero_atm, faero_ocn
  use ice_flux_bgc       , only : fiso_atm, fiso_ocn, fiso_evap
  use ice_flux_bgc       , only : Qa_iso, Qref_iso, HDO_ocn, H2_18O_ocn, H2_16O_ocn
  use ice_flux           , only : fresh, fsalt, zlvl, uatm, vatm, potT, Tair, Qa
  use ice_flux           , only : rhoa, swvdr, swvdf, swidr, swidf, flw, frain
  use ice_flux           , only : fsnow, uocn, vocn, sst, ss_tltx, ss_tlty, frzmlt
  use ice_flux           , only : send_i2x_per_cat
  use ice_flux           , only : sss, Tf, wind, fsw
  use ice_arrays_column  , only : floe_rad_c, wave_spectrum
  use ice_state          , only : vice, vsno, aice, aicen_init, trcr, trcrn
  use ice_grid           , only : tlon, tlat, tarea, tmask, anglet, hm
  use ice_grid           , only : grid_format
  use ice_mesh_mod       , only : ocn_gridcell_frac
  use ice_boundary       , only : ice_HaloUpdate
  use ice_fileunits      , only : nu_diag, flush_fileunit
  use ice_communicate    , only : my_task, master_task, MPI_COMM_ICE
  use ice_prescribed_mod , only : prescribed_ice
  use ice_shr_methods    , only : chkerr, state_reset
  use icepack_intfc      , only : icepack_warnings_flush, icepack_warnings_aborted
  use icepack_intfc      , only : icepack_query_parameters, icepack_query_tracer_flags
  use icepack_intfc      , only : icepack_liquidus_temperature
  use icepack_intfc      , only : icepack_sea_freezing_temperature
  use icepack_intfc      , only : icepack_query_tracer_indices
  use icepack_parameters , only : puny, c2
  use cice_wrapper_mod   , only : t_startf, t_stopf, t_barrierf
#ifdef CESMCOUPLED
  use shr_frz_mod        , only : shr_frz_freezetemp
  use shr_mpi_mod        , only : shr_mpi_min, shr_mpi_max
#endif

  implicit none
  public

  public  :: ice_advertise_fields
  public  :: ice_realize_fields
  public  :: ice_import
  public  :: ice_export

  private :: fldlist_add
  private :: fldlist_realize
  private :: state_FldChk

  interface state_getfldptr
     module procedure state_getfldptr_1d
     module procedure state_getfldptr_2d
  end interface state_getfldptr
  private :: state_getfldptr

  interface state_getimport
     module procedure state_getimport_4d
     module procedure state_getimport_3d
  end interface state_getimport
  private :: state_getimport

  interface state_setexport
     module procedure state_setexport_4d
     module procedure state_setexport_3d
  end interface state_setexport
  private :: state_setexport

  ! Private module data

  type fld_list_type
    character(char_len) :: stdname
    integer :: ungridded_lbound = 0
    integer :: ungridded_ubound = 0
  end type fld_list_type

  ! area correction factors for fluxes send and received from mediator
  real(dbl_kind), allocatable :: mod2med_areacor(:) ! ratios of model areas to input mesh areas
  real(dbl_kind), allocatable :: med2mod_areacor(:) ! ratios of input mesh areas to model areas

  integer, parameter       :: fldsMax = 100
  integer                  :: fldsToIce_num = 0
  integer                  :: fldsFrIce_num = 0
  type (fld_list_type)     :: fldsToIce(fldsMax)
  type (fld_list_type)     :: fldsFrIce(fldsMax)

  logical                  :: flds_wave           ! wave ice coupling
  integer     , parameter  :: io_dbug = 10        ! i/o debug messages
  character(*), parameter  :: u_FILE_u = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine ice_advertise_fields(gcomp, importState, exportState, flds_scalar_name, rc)

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    type(ESMF_State)               :: importState
    type(ESMF_State)               :: exportState
    character(len=*) , intent(in)  :: flds_scalar_name
    integer          , intent(out) :: rc

    ! local variables
    integer             :: n
    character(char_len) :: stdname
    character(char_len) :: cvalue
    logical             :: flds_wiso         ! use case
    logical             :: isPresent, isSet
    character(len=*), parameter :: subname='(ice_import_export:ice_advertise_fields)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (io_dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Determine if ice sends multiple ice category info back to mediator
    send_i2x_per_cat = .false.
    call NUOPC_CompAttributeGet(gcomp, name='flds_i2o_per_cat', value=cvalue, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) send_i2x_per_cat
    end if
    if (my_task == master_task) then
       write(nu_diag,*)'send_i2x_per_cat = ',send_i2x_per_cat
    end if
    if (.not.send_i2x_per_cat) then
       if (allocated(fswthrun_ai)) then
          deallocate(fswthrun_ai)
       end if
    end if

    ! Determine if the following attributes are sent by the driver and if so read them in
    flds_wiso = .false.
    call NUOPC_CompAttributeGet(gcomp, name='flds_wiso', value=cvalue, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_wiso
    end if
    if (my_task == master_task) then
       write(nu_diag,*)'flds_wiso = ',flds_wiso
    end if

    flds_wave = .false.
    call NUOPC_CompAttributeGet(gcomp, name='wav_coupling_to_cice', value=cvalue, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_wave
    end if
    if (my_task == master_task) then
       write(nu_diag,*)'flds_wave = ',flds_wave
    end if

    !-----------------
    ! advertise import fields
    !-----------------

    call fldlist_add(fldsToIce_num, fldsToIce, trim(flds_scalar_name))

    ! from ocean
    call fldlist_add(fldsToIce_num, fldsToIce, 'So_dhdx' )
    call fldlist_add(fldsToIce_num, fldsToIce, 'So_dhdy' )
    call fldlist_add(fldsToIce_num, fldsToIce, 'So_t'    )
    call fldlist_add(fldsToIce_num, fldsToIce, 'So_s'    )
    call fldlist_add(fldsToIce_num, fldsToIce, 'So_u'    )
    call fldlist_add(fldsToIce_num, fldsToIce, 'So_v'    )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Fioo_q'  )
    if (flds_wiso) then
       call fldlist_add(fldsToIce_num, fldsToIce, 'So_roce_wiso', ungridded_lbound=1, ungridded_ubound=3)
    end if

    ! from atmosphere
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_z'       )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_u'       )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_v'       )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_shum'    )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_tbot'    )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_pbot'    )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_swvdr' )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_swvdf' )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_swndr' )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_swndf' )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_lwdn'  )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_rain'  )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_snow'  )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_ptem'    ) !cesm
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_dens'    ) !cesm

    ! the following are advertised but might not be connected if they are not present
    ! in the cmeps esmFldsExchange_xxx_mod.F90 that is model specific
    ! from atm - black carbon deposition fluxes (3)
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_bcph',  ungridded_lbound=1, ungridded_ubound=3)
    ! from atm - wet dust deposition fluxes (4 sizes)
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_dstwet', ungridded_lbound=1, ungridded_ubound=4)
    ! from atm - dry dust deposition fluxes (4 sizes)
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_dstdry', ungridded_lbound=1, ungridded_ubound=4)

    ! the following are advertised but might not be connected if they are not advertised in the
    ! in the cmeps esmFldsExchange_xxx_mod.F90 that is model specific
    ! from wave
    if (flds_wave) then
       call fldlist_add(fldsToIce_num, fldsToIce, 'Sw_elevation_spectrum', ungridded_lbound=1, &
            ungridded_ubound=25)
    end if

    do n = 1,fldsToIce_num
       call NUOPC_Advertise(importState, standardName=fldsToIce(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    !-----------------
    ! advertise export fields
    !-----------------

    call fldlist_add(fldsFrIce_num, fldsFrIce, trim(flds_scalar_name))

    ! ice states
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_imask' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_ifrac' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_t'     )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_vice'  )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_vsno'  )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_tref'  )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_qref'  )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_snowh' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_u10'   )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_avsdr' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_avsdf' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_anidr' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_anidf' )

    ! the following are advertised but might not be connected if they are not present
    ! in the cmeps esmFldsExchange_xxx_mod.F90 that is model specific
    if (send_i2x_per_cat) then
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_ifrac_n', &
            ungridded_lbound=1, ungridded_ubound=ncat)
    end if
    if (flds_wave) then
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_thick'    )
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_floediam' )
    end if

    ! ice/atm fluxes computed by ice
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_taux'      )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_tauy'      )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_lat'       )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_sen'       )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_lwup'      )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_evap'      )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_swnet'     )

    ! ice/ocn fluxes computed by ice
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_melth'     )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_swpen'     )

    if (.not.prescribed_ice) then
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_swpen_vdr' )
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_swpen_vdf' )
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_swpen_idr' )
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_swpen_idf' )
    endif

    if (send_i2x_per_cat) then
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_swpen_ifrac_n', &
            ungridded_lbound=1, ungridded_ubound=ncat)
    end if
    call fldlist_add(fldsFrIce_num , fldsFrIce, 'Fioi_meltw' )
    call fldlist_add(fldsFrIce_num , fldsFrIce, 'Fioi_salt'  )
    call fldlist_add(fldsFrIce_num , fldsFrIce, 'Fioi_taux'  )
    call fldlist_add(fldsFrIce_num , fldsFrIce, 'Fioi_tauy'  )

    ! the following are advertised but might not be connected if they are not present
    ! in the cmeps esmFldsExchange_xxx_mod.F90 that is model specific
    call fldlist_add(fldsFrIce_num , fldsFrIce, 'Fioi_bcpho'  )
    call fldlist_add(fldsFrIce_num , fldsFrIce, 'Fioi_bcphi'  )
    call fldlist_add(fldsFrIce_num , fldsFrIce, 'Fioi_flxdst' )

    if (flds_wiso) then
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_meltw_wiso', &
            ungridded_lbound=1, ungridded_ubound=3)
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_evap_wiso', &
            ungridded_lbound=1, ungridded_ubound=3)
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_qref_wiso', &
            ungridded_lbound=1, ungridded_ubound=3)
    end if

    do n = 1,fldsFrIce_num
       call NUOPC_Advertise(exportState, standardName=fldsFrIce(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    if (io_dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ice_advertise_fields

  !==============================================================================
  subroutine ice_realize_fields(gcomp, mesh, flds_scalar_name, flds_scalar_num, rc)
    use ice_scam, only : single_column

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    type(ESMF_Mesh)  , intent(in)  :: mesh
    character(len=*) , intent(in)  :: flds_scalar_name
    integer          , intent(in)  :: flds_scalar_num
    integer          , intent(out) :: rc

    ! local variables
    type(ESMF_State)            :: importState
    type(ESMF_State)            :: exportState
    type(ESMF_Field)            :: lfield
    integer                     :: numOwnedElements
    integer                     :: i, j, iblk, n
    integer                     :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    type(block)                 :: this_block         ! block information for current block
    real(dbl_kind), allocatable :: mesh_areas(:)
    real(dbl_kind), allocatable :: model_areas(:)
    real(dbl_kind), pointer     :: dataptr(:)
    real(dbl_kind)              :: max_mod2med_areacor
    real(dbl_kind)              :: max_med2mod_areacor
    real(dbl_kind)              :: min_mod2med_areacor
    real(dbl_kind)              :: min_med2mod_areacor
    real(dbl_kind)              :: max_mod2med_areacor_glob
    real(dbl_kind)              :: max_med2mod_areacor_glob
    real(dbl_kind)              :: min_mod2med_areacor_glob
    real(dbl_kind)              :: min_med2mod_areacor_glob
    character(len=*), parameter :: subname='(ice_import_export:realize_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=ExportState, &
         fldList=fldsFrIce, &
         numflds=fldsFrIce_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':CICE_Export',&
         mesh=mesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=importState, &
         fldList=fldsToIce, &
         numflds=fldsToIce_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':CICE_Import',&
         mesh=mesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
#ifdef CESMCOUPLED

    ! allocate area correction factors
    call ESMF_MeshGet(mesh, numOwnedElements=numOwnedElements, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate (mod2med_areacor(numOwnedElements))
    allocate (med2mod_areacor(numOwnedElements))

    if (single_column) then

       mod2med_areacor(:) = 1._dbl_kind
       med2mod_areacor(:) = 1._dbl_kind

    else

       ! Get mesh areas from second field - using second field since the
       ! first field is the scalar field

       call ESMF_StateGet(exportState, itemName=trim(fldsFrIce(2)%stdname), field=lfield, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldRegridGetArea(lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayPtr=dataptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       allocate(mesh_areas(numOwnedElements))
       mesh_areas(:) = dataptr(:)

       ! Determine flux correction factors (module variables)
       allocate(model_areas(numOwnedElements))
       mod2med_areacor(:) = 1._dbl_kind
       med2mod_areacor(:) = 1._dbl_kind
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
                model_areas(n) = tarea(i,j,iblk)/(radius*radius)
                mod2med_areacor(n) = model_areas(n) / mesh_areas(n)
                med2mod_areacor(n) = mesh_areas(n) / model_areas(n)
             enddo
          enddo
       enddo
       deallocate(model_areas)
       deallocate(mesh_areas)
    end if

    min_mod2med_areacor = minval(mod2med_areacor)
    max_mod2med_areacor = maxval(mod2med_areacor)
    min_med2mod_areacor = minval(med2mod_areacor)
    max_med2mod_areacor = maxval(med2mod_areacor)
    call shr_mpi_max(max_mod2med_areacor, max_mod2med_areacor_glob, mpi_comm_ice)
    call shr_mpi_min(min_mod2med_areacor, min_mod2med_areacor_glob, mpi_comm_ice)
    call shr_mpi_max(max_med2mod_areacor, max_med2mod_areacor_glob, mpi_comm_ice)
    call shr_mpi_min(min_med2mod_areacor, min_med2mod_areacor_glob, mpi_comm_ice)

    if (my_task == master_task) then
       write(nu_diag,'(2A,2g23.15,A )') trim(subname),' :  min_mod2med_areacor, max_mod2med_areacor ',&
            min_mod2med_areacor_glob, max_mod2med_areacor_glob, 'CICE6'
       write(nu_diag,'(2A,2g23.15,A )') trim(subname),' :  min_med2mod_areacor, max_med2mod_areacor ',&
            min_med2mod_areacor_glob, max_med2mod_areacor_glob, 'CICE6'
    end if
#endif

  end subroutine ice_realize_fields

  !==============================================================================
  subroutine ice_import( importState, rc )

    ! input/output variables
    type(ESMF_State) , intent(in)  :: importState
    integer          , intent(out) :: rc

    ! local variables
    integer,parameter                :: nflds=16
    integer,parameter                :: nfldv=6
    integer                          :: i, j, iblk, n, k
    integer                          :: ilo, ihi, jlo, jhi !beginning and end of physical domain
    type(block)                      :: this_block         ! block information for current block
    real (kind=dbl_kind),allocatable :: aflds(:,:,:,:)
    real (kind=dbl_kind)             :: workx, worky
    real (kind=dbl_kind)             :: MIN_RAIN_TEMP, MAX_SNOW_TEMP
    real (kind=dbl_kind)             :: Tffresh
    real (kind=dbl_kind)             :: inst_pres_height_lowest
    real (kind=dbl_kind), pointer    :: dataptr2d(:,:)
    real (kind=dbl_kind), pointer    :: dataptr1d(:)
    real (kind=dbl_kind), pointer    :: dataptr2d_dstwet(:,:)
    real (kind=dbl_kind), pointer    :: dataptr2d_dstdry(:,:)
    character(len=char_len)          :: tfrz_option
    integer(int_kind)                :: ktherm
    character(len=*),   parameter    :: subname = 'ice_import'
    character(len=1024)              :: msgString
    !-----------------------------------------------------

    call icepack_query_parameters(Tffresh_out=Tffresh)
    call icepack_query_parameters(tfrz_option_out=tfrz_option)
    call icepack_query_parameters(ktherm_out=ktherm)

    if (io_dbug > 5) then
       write(msgString,'(A,i8)')trim(subname)//' tfrz_option = ' &
            // trim(tfrz_option)//', ktherm = ',ktherm
       call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
    end if

    !    call icepack_query_parameters(tfrz_option_out=tfrz_option, &
    !       modal_aero_out=modal_aero, z_tracers_out=z_tracers, skl_bgc_out=skl_bgc, &
    !       Tffresh_out=Tffresh)
    !    call icepack_query_tracer_flags(tr_aero_out=tr_aero, tr_iage_out=tr_iage, &
    !       tr_FY_out=tr_FY, tr_pond_out=tr_pond, tr_lvl_out=tr_lvl, &
    !       tr_zaero_out=tr_zaero, tr_bgc_Nit_out=tr_bgc_Nit)

    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
        file=u_FILE_u, line=__LINE__)

    ! Note that the precipitation fluxes received from the mediator
    ! are in units of kg/s/m^2 which is what CICE requires.
    ! Note also that the read in below includes only values needed
    ! by the thermodynamic component of CICE.  Variables uocn, vocn,
    ! ss_tltx, and ss_tlty are excluded. Also, because the SOM and
    ! DOM don't  compute SSS.   SSS is not read in and is left at
    ! the initilized value (see ice_flux.F init_coupler_flux) of
    ! 34 ppt

    ! Use aflds to gather the halo updates of multiple fields
    ! Need to separate the scalar from the vector halo updates

    allocate(aflds(nx_block,ny_block,nflds,nblocks))
    aflds = c0

    ! import ocean states

    call state_getimport(importState, 'So_t', output=aflds, index=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'So_s', output=aflds, index=2, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! import atm states

    call state_getimport(importState, 'Sa_z', output=aflds, index=3, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (State_FldChk(importState, 'Sa_ptem') .and. State_fldchk(importState, 'Sa_dens')) then
       call state_getimport(importState, 'Sa_ptem', output=aflds, index=4, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'Sa_dens', output=aflds, index=5, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else if (State_FldChk(importState, 'Sa_pbot')) then
       call state_getimport(importState, 'Sa_pbot', output=aflds, index=6, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call abort_ice(trim(subname)//&
            ": ERROR either Sa_ptem and Sa_dens OR Sa_pbot must be in import state")
    end if

    call state_getimport(importState, 'Sa_tbot', output=aflds, index=7, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Sa_shum', output=aflds, index=8, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! import ocn/ice fluxes

    call state_getimport(importState, 'Fioo_q', output=aflds, index=9, &
         areacor=med2mod_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! import atm fluxes

    call state_getimport(importState, 'Faxa_swvdr', output=aflds, index=10, &
         areacor=med2mod_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_swndr', output=aflds, index=11, &
         areacor=med2mod_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_swvdf', output=aflds, index=12, &
         areacor=med2mod_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_swndf', output=aflds, index=13, &
         areacor=med2mod_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_lwdn', output=aflds, index=14, &
         areacor=med2mod_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_rain', output=aflds, index=15, &
         areacor=med2mod_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_snow', output=aflds, index=16, &
         areacor=med2mod_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! perform a halo update

    if (.not.prescribed_ice) then
       call t_startf ('cice_imp_halo')
       call ice_HaloUpdate(aflds, halo_info, field_loc_center, field_type_scalar)
       call t_stopf ('cice_imp_halo')
    endif

    ! now fill in the ice internal data types

    !$OMP PARALLEL DO PRIVATE(iblk,i,j)
    do iblk = 1, nblocks
       do j = 1,ny_block
          do i = 1,nx_block
             sst  (i,j,iblk)         = aflds(i,j, 1,iblk)
             sss  (i,j,iblk)         = aflds(i,j, 2,iblk)
             zlvl (i,j,iblk)         = aflds(i,j, 3,iblk)
             ! see below for 4,5,6
             Tair (i,j,iblk)         = aflds(i,j, 7,iblk)
             Qa   (i,j,iblk)         = aflds(i,j, 8,iblk)
             frzmlt (i,j,iblk)       = aflds(i,j, 9,iblk)
             swvdr(i,j,iblk)         = aflds(i,j,10,iblk)
             swidr(i,j,iblk)         = aflds(i,j,11,iblk)
             swvdf(i,j,iblk)         = aflds(i,j,12,iblk)
             swidf(i,j,iblk)         = aflds(i,j,13,iblk)
             flw  (i,j,iblk)         = aflds(i,j,14,iblk)
             frain(i,j,iblk)         = aflds(i,j,15,iblk)
             fsnow(i,j,iblk)         = aflds(i,j,16,iblk)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    ! import wave elevation spectrum from wave  (frequencies 1-25, assume that nfreq is 25)
    if (State_FldChk(importState, 'Sw_elevation_spectrum')) then
       if (nfreq /= 25) then
          call abort_ice(trim(subname)//": ERROR nfreq not equal to 25 ")
       end if
       call state_getfldptr(importState, 'Sw_elevation_spectrum', fldptr=dataPtr2d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do k = 1,nfreq
          n = 0
          do iblk = 1, nblocks
             this_block = get_block(blocks_ice(iblk),iblk)
             ilo = this_block%ilo; ihi = this_block%ihi
             jlo = this_block%jlo; jhi = this_block%jhi
             do j = jlo, jhi
                do i = ilo, ihi
                   n = n+1
                   wave_spectrum(i,j,k,iblk) = dataPtr2d(k,n)
                end do
             end do
          end do
       end do
    end if

    if ( State_fldChk(importState, 'Sa_ptem') .and. State_fldchk(importState,'Sa_dens')) then
       !$OMP PARALLEL DO PRIVATE(iblk,i,j)
       do iblk = 1, nblocks
          do j = 1,ny_block
             do i = 1,nx_block
                potT (i,j,iblk) = aflds(i,j, 4,iblk)
                rhoa (i,j,iblk) = aflds(i,j, 5,iblk)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else if (State_fldChk(importState, 'Sa_pbot')) then
       !$OMP PARALLEL DO PRIVATE(iblk,i,j)
       do iblk = 1, nblocks
          do j = 1,ny_block
             do i = 1,nx_block
                inst_pres_height_lowest = aflds(i,j,6,iblk)
                if (inst_pres_height_lowest > 0.0_ESMF_KIND_R8) then
                   potT (i,j,iblk) = Tair(i,j,iblk) * (100000._ESMF_KIND_R8/inst_pres_height_lowest)**0.286_ESMF_KIND_R8
                else
                   potT (i,j,iblk) = 0.0_ESMF_KIND_R8
                end if
                if (Tair(i,j,iblk) /= 0._ESMF_KIND_R8) then
                   rhoa(i,j,iblk) = inst_pres_height_lowest / &
                        (287.058_ESMF_KIND_R8*(1._ESMF_KIND_R8+0.608_ESMF_KIND_R8*Qa(i,j,iblk))*Tair(i,j,iblk))
                else
                   rhoa(i,j,iblk) = 1.2_ESMF_KIND_R8
                endif
             end do !i
          end do !j
       end do !iblk
       !$OMP END PARALLEL DO
    end if

    deallocate(aflds)
    allocate(aflds(nx_block,ny_block,nfldv,nblocks))
    aflds = c0

    ! Get velocity fields from ocean and atm and slope fields from ocean

    call state_getimport(importState, 'So_u', output=aflds, index=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'So_v', output=aflds, index=2, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Sa_u', output=aflds, index=3, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'Sa_v', output=aflds, index=4, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'So_dhdx', output=aflds, index=5, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'So_dhdy', output=aflds, index=6, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return


    if (.not.prescribed_ice) then
       call t_startf ('cice_imp_halo')
       call ice_HaloUpdate(aflds, halo_info, field_loc_center, field_type_vector)
       call t_stopf ('cice_imp_halo')
    endif

    !$OMP PARALLEL DO PRIVATE(iblk,i,j)
    do iblk = 1, nblocks
       do j = 1,ny_block
          do i = 1,nx_block
             uocn (i,j,iblk)   = aflds(i,j, 1,iblk)
             vocn (i,j,iblk)   = aflds(i,j, 2,iblk)
             uatm (i,j,iblk)   = aflds(i,j, 3,iblk)
             vatm (i,j,iblk)   = aflds(i,j, 4,iblk)
             ss_tltx(i,j,iblk) = aflds(i,j, 5,iblk)
             ss_tlty(i,j,iblk) = aflds(i,j, 6,iblk)
          enddo  !i
       enddo     !j
    enddo        !iblk
    !$OMP END PARALLEL DO

    deallocate(aflds)

    !-------------------------------------------------------
    ! Get aerosols from mediator
    !-------------------------------------------------------

    if (State_FldChk(importState, 'Faxa_bcph')) then
       ! the following indices are based on what the atmosphere is sending
       ! bcphidry  ungridded_index=1
       ! bcphodry  ungridded_index=2
       ! bcphiwet  ungridded_index=3

       call state_getfldptr(importState, 'Faxa_bcph', fldptr=dataPtr2d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       n = 0
       do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo; ihi = this_block%ihi
          jlo = this_block%jlo; jhi = this_block%jhi
          do j = jlo, jhi
             do i = ilo, ihi
                n = n+1
                faero_atm(i,j,1,iblk)  = dataPtr2d(2,n) * med2mod_areacor(n) ! bcphodry
                faero_atm(i,j,2,iblk)  = (dataptr2d(1,n) + dataPtr2d(3,n)) * med2mod_areacor(n) ! bcphidry + bcphiwet
             end do
          end do
       end do
    end if

    ! Sum over all dry and wet dust fluxes from ath atmosphere
    if (State_FldChk(importState, 'Faxa_dstwet') .and. State_FldChk(importState, 'Faxa_dstdry')) then
       call state_getfldptr(importState, 'Faxa_dstwet', fldptr=dataPtr2d_dstwet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, 'Faxa_dstdry', fldptr=dataPtr2d_dstdry, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       n = 0
       do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo; ihi = this_block%ihi
          jlo = this_block%jlo; jhi = this_block%jhi
          do j = jlo, jhi
             do i = ilo, ihi
                n = n+1
                faero_atm(i,j,3,iblk)  = dataPtr2d_dstwet(1,n) + dataptr2d_dstdry(1,n) + &
                                         dataPtr2d_dstwet(2,n) + dataptr2d_dstdry(2,n) + &
                                         dataPtr2d_dstwet(3,n) + dataptr2d_dstdry(3,n) + &
                                         dataPtr2d_dstwet(4,n) + dataptr2d_dstdry(4,n)
                faero_atm(i,j,3,iblk) = faero_atm(i,j,3,iblk) * med2mod_areacor(n)
             end do
          end do
       end do
    end if

    !-------------------------------------------------------
    ! Water isotopes from the mediator
    !-------------------------------------------------------

    ! 16O => ungridded_index=1
    ! 18O => ungridded_index=2
    ! HDO => ungridded_index=3

    if (State_FldChk(importState, 'shum_wiso')) then
       call state_getimport(importState, 'Sa_shum_wiso', output=Qa_iso, index=1, ungridded_index=3, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'Sa_shum_wiso', output=Qa_iso, index=2, ungridded_index=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'Sa_shum_wiso', output=Qa_iso, index=3, ungridded_index=2, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! call state_getimport(importState, 'mean_prec_rate_wiso', output=fiso_rain, index=1, ungridded_index=3, &
       !      areacor=med2mod_areacor, rc=rc)
       ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! call state_getimport(importState, 'mean_prec_rate_wiso', output=fiso_rain, index=2, ungridded_index=1, &
       !      areacor=med2mod_areacor, rc=rc)
       ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! call state_getimport(importState, 'mean_prec_rate_wiso', output=fiso_rain, index=3, ungridded_index=2, &
       !      areacor=med2mod_areacor, rc=rc)
       ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call state_getimport(importState, 'Faxa_snow_wiso', output=fiso_atm, index=1, ungridded_index=3, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'Faxa_snow_wiso', output=fiso_atm, index=2, ungridded_index=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'Faxa_snow_wiso', output=fiso_atm, index=3, ungridded_index=2, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call state_getimport(importState, 'So_roce_wiso', output=HDO_ocn   , ungridded_index=3, &
            areacor=med2mod_areacor, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'So_roce_wiso', output=H2_16O_ocn, ungridded_index=1, &
            areacor=med2mod_areacor, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'So_roce_wiso', output=H2_18O_ocn, ungridded_index=2, &
            areacor=med2mod_areacor, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !-----------------------------------------------------------------
    ! rotate zonal/meridional vectors to local coordinates
    ! compute data derived quantities
    !-----------------------------------------------------------------

    ! Vector fields come in on T grid, but are oriented geographically
    ! need to rotate to pop-grid FIRST using ANGLET
    ! then interpolate to the U-cell centers  (otherwise we
    ! interpolate across the pole)
    ! use ANGLET which is on the T grid !

    call t_startf ('cice_imp_ocn')

    !$OMP PARALLEL DO PRIVATE(iblk,i,j,workx,worky)
    do iblk = 1, nblocks

       do j = 1,ny_block
          do i = 1,nx_block
             ! ocean
             workx      = uocn  (i,j,iblk) ! currents, m/s
             worky      = vocn  (i,j,iblk)

             uocn(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! rotate to align with model i,j
                            + worky*sin(ANGLET(i,j,iblk))
             vocn(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                            - workx*sin(ANGLET(i,j,iblk))

             workx      = ss_tltx  (i,j,iblk)           ! sea sfc tilt, m/m
             worky      = ss_tlty  (i,j,iblk)

             ss_tltx(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! rotate to align with model i,j
                               + worky*sin(ANGLET(i,j,iblk))
             ss_tlty(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                               - workx*sin(ANGLET(i,j,iblk))

             sst(i,j,iblk) = sst(i,j,iblk) - Tffresh       ! sea sfc temp (C)

             sss(i,j,iblk) = max(sss(i,j,iblk),c0)

          enddo
       enddo
    end do

#ifdef CESMCOUPLED
    ! Use shr_frz_mod for this
    do iblk = 1, nblocks
       Tf(:,:,iblk) = shr_frz_freezetemp(sss(:,:,iblk))
    end do
#else
    !$OMP PARALLEL DO PRIVATE(iblk,i,j)
    do iblk = 1, nblocks
       do j = 1,ny_block
          do i = 1,nx_block
            Tf(i,j,iblk) = icepack_sea_freezing_temperature(sss(i,j,iblk))
          end do
       end do
    end do
    !$OMP END PARALLEL DO
#endif

    call t_stopf ('cice_imp_ocn')

    ! Interpolate ocean dynamics variables from T-cell centers to
    ! U-cell centers.

    if (.not.prescribed_ice) then
       call t_startf ('cice_imp_t2u')
       call ice_HaloUpdate(uocn, halo_info, field_loc_center, field_type_vector)
       call ice_HaloUpdate(vocn, halo_info, field_loc_center, field_type_vector)
       call ice_HaloUpdate(ss_tltx, halo_info, field_loc_center, field_type_vector)
       call ice_HaloUpdate(ss_tlty, halo_info, field_loc_center, field_type_vector)
       call t_stopf ('cice_imp_t2u')
    end if

    ! Atmosphere variables are needed in T cell centers in
    ! subroutine stability and are interpolated to the U grid
    ! later as necessary.

    call t_startf ('cice_imp_atm')
    !$OMP PARALLEL DO PRIVATE(iblk,i,j,workx,worky)
    do iblk = 1, nblocks
       do j = 1, ny_block
          do i = 1, nx_block

             ! atmosphere
             workx      = uatm(i,j,iblk) ! wind velocity, m/s
             worky      = vatm(i,j,iblk)
             uatm (i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid
                             + worky*sin(ANGLET(i,j,iblk))   ! note uatm, vatm, wind
             vatm (i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) & ! are on the T-grid here
                             - workx*sin(ANGLET(i,j,iblk))

             wind (i,j,iblk) = sqrt(uatm(i,j,iblk)**2 + vatm(i,j,iblk)**2)
             fsw  (i,j,iblk) = swvdr(i,j,iblk) + swvdf(i,j,iblk) &
                             + swidr(i,j,iblk) + swidf(i,j,iblk)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call t_stopf ('cice_imp_atm')

  end subroutine ice_import

  !===============================================================================
  subroutine ice_export( exportState, rc )

    use ice_scam, only : single_column

    ! input/output variables
    type(ESMF_State), intent(inout) :: exportState
    integer         , intent(out)   :: rc

    ! local variables
    type(block)             :: this_block                           ! block information for current block
    integer                 :: i, j, iblk, n, k                     ! indices
    integer                 :: n2                                   ! thickness category index
    integer                 :: ilo, ihi, jlo, jhi                   ! beginning and end of physical domain
    real    (kind=dbl_kind) :: workx, worky                         ! tmps for converting grid
    integer (kind=int_kind) :: icells                               ! number of ocean/ice cells
    logical                 :: flag
    integer (kind=int_kind) :: indxi (nx_block*ny_block)            ! compressed indices in i
    integer (kind=int_kind) :: indxj (nx_block*ny_block)            ! compressed indices in i
    real    (kind=dbl_kind) :: Tsrf  (nx_block,ny_block,max_blocks) ! surface temperature
    real    (kind=dbl_kind) :: tauxa (nx_block,ny_block,max_blocks) ! atmo/ice stress
    real    (kind=dbl_kind) :: tauya (nx_block,ny_block,max_blocks) ! atm/ice stress
    real    (kind=dbl_kind) :: tauxo (nx_block,ny_block,max_blocks) ! ice/ocean stress
    real    (kind=dbl_kind) :: tauyo (nx_block,ny_block,max_blocks) ! ice/ocean stress
    real    (kind=dbl_kind) :: ailohi(nx_block,ny_block,max_blocks) ! fractional ice area
    real    (kind=dbl_kind) :: floediam(nx_block,ny_block,max_blocks)
    real    (kind=dbl_kind) :: floethick(nx_block,ny_block,max_blocks) ! ice thickness
    logical (kind=log_kind) :: tr_fsd
    integer (kind=int_kind) :: nt_fsd
    real    (kind=dbl_kind) :: Tffresh, stefan_boltzmann
    real    (kind=dbl_kind), allocatable :: tempfld(:,:,:)
    real    (kind=dbl_kind), pointer :: dataptr_ifrac_n(:,:)
    real    (kind=dbl_kind), pointer :: dataptr_swpen_n(:,:)
    logical (kind=log_kind), save :: first_call = .true.
    character(len=*),parameter :: subname = 'ice_export'
    !-----------------------------------------------------

    rc = ESMF_SUCCESS
    if (io_dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call icepack_query_parameters(Tffresh_out=Tffresh)
    call icepack_query_parameters(stefan_boltzmann_out=stefan_boltzmann)
    !    call icepack_query_parameters(tfrz_option_out=tfrz_option, &
    !       modal_aero_out=modal_aero, z_tracers_out=z_tracers, skl_bgc_out=skl_bgc, &
    !       Tffresh_out=Tffresh)
    !    call icepack_query_tracer_flags(tr_aero_out=tr_aero, tr_iage_out=tr_iage, &
    !       tr_FY_out=tr_FY, tr_pond_out=tr_pond, tr_lvl_out=tr_lvl, &
    !       tr_zaero_out=tr_zaero, tr_bgc_Nit_out=tr_bgc_Nit)

    call icepack_query_tracer_indices(nt_fsd_out=nt_fsd)
    call icepack_query_tracer_flags(tr_fsd_out=tr_fsd)

    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
        file=u_FILE_u, line=__LINE__)

    !calculate ice thickness from aice and vice. Also
    !create Tsrf from the first tracer (trcr) in ice_state.F

    ailohi(:,:,:) = c0
    Tsrf(:,:,:)  = c0
    tauxa(:,:,:) = c0
    tauya(:,:,:) = c0
    tauxo(:,:,:) = c0
    tauyo(:,:,:) = c0
    floediam(:,:,:) = c0
    floethick(:,:,:) = c0

    !$OMP PARALLEL DO PRIVATE(iblk,i,j,k,workx,worky, this_block, ilo, ihi, jlo, jhi)
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo,jhi
          do i = ilo,ihi
             ! ice fraction
             ailohi(i,j,iblk) = min(aice(i,j,iblk), c1)

             ! surface temperature
             Tsrf(i,j,iblk)  = Tffresh + trcr(i,j,1,iblk)     !Kelvin (original ???)

             if (flds_wave) then
                ! floe thickness (m)
                if (aice(i,j,iblk) > puny) then
                   floethick(i,j,iblk) = vice(i,j,iblk) / aice(i,j,iblk)
                else
                   floethick(i,j,iblk) = c0
                end if

                if (tr_fsd) then
                   ! floe diameter (m)
                   workx = c0
                   worky = c0
                   do n = 1, ncat
                      do k = 1, nfsd
                         workx = workx + floe_rad_c(k) * aicen_init(i,j,n,iblk) * trcrn(i,j,nt_fsd+k-1,n,iblk)
                         worky = worky + aicen_init(i,j,n,iblk) * trcrn(i,j,nt_fsd+k-1,n,iblk)
                      end do
                   end do
                   if (worky > c0) workx = c2*workx / worky
                   floediam(i,j,iblk) = MAX(c2*floe_rad_c(1),workx)
                else ! with FSD off
                   ! floe diameter (m)
                   floediam(i,j,iblk) = 50.0_dbl_kind
                endif
             endif

             ! wind stress  (on POP T-grid:  convert to lat-lon)
             workx = strairxT(i,j,iblk)                             ! N/m^2
             worky = strairyT(i,j,iblk)                             ! N/m^2
             tauxa(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) - worky*sin(ANGLET(i,j,iblk))
             tauya(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) + workx*sin(ANGLET(i,j,iblk))

             ! ice/ocean stress (on POP T-grid:  convert to lat-lon)
             workx = -strocnxT_iavg(i,j,iblk)                       ! N/m^2
             worky = -strocnyT_iavg(i,j,iblk)                       ! N/m^2
             tauxo(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) - worky*sin(ANGLET(i,j,iblk))
             tauyo(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) + workx*sin(ANGLET(i,j,iblk))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    flag=.false.
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo,jhi
          do i = ilo,ihi
             if (tmask(i,j,iblk) .and. ailohi(i,j,iblk) < c0 ) then
                flag = .true.
             endif
          end do
       end do
    end do
    if (flag) then
       do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi
          do j = jlo,jhi
             do i = ilo,ihi
                if (tmask(i,j,iblk) .and. ailohi(i,j,iblk) < c0 ) then
                   write(nu_diag,*) &
                        ' (ice) send: ERROR ailohi < 0.0 ',i,j,ailohi(i,j,iblk)
                   call flush_fileunit(nu_diag)
                endif
             end do
          end do
       end do
    endif

    !---------------------------------
    ! Create the export state
    !---------------------------------

    ! Zero out fields with tmask for proper coupler accumulation in ice free areas
    if (first_call .or. .not.single_column) then
       call state_reset(exportState, c0, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       first_call = .false.
    endif

    ! Create a temporary field
    allocate(tempfld(nx_block,ny_block,nblocks))

    ! Fractions and mask
    call state_setexport(exportState, 'Si_ifrac', input=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (trim(grid_format) == 'meshnc') then
       call state_setexport(exportState, 'Si_imask', input=ocn_gridcell_frac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       tempfld(:,:,:) = c0
       do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi
          do j = jlo, jhi
             do i = ilo, ihi
                tempfld(i,j,iblk) = real(nint(hm(i,j,iblk)),kind=dbl_kind)
             end do
          end do
       end do
       call state_setexport(exportState, 'Si_imask', input=tempfld, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! ----
    ! States from ice
    ! ----

    ! surface temperature of ice covered portion (degK)
    call state_setexport(exportState, 'Si_t', input=Tsrf , lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! albedo vis dir
    call state_setexport(exportState, 'Si_avsdr', input=alvdr, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! albedo nir dir
    call state_setexport(exportState, 'Si_anidr', input=alidr, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! albedo vis dif
    call state_setexport(exportState, 'Si_avsdf', input=alvdf, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! albedo nir dif
    call state_setexport(exportState, 'Si_anidf', input=alidf, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! 10m atm reference wind speed (m/s)
    call state_setexport(exportState, 'Si_u10'  , input=Uref , lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! 2m atm reference temperature (K)
    call state_setexport(exportState, 'Si_tref' , input=Tref , lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! 2m atm reference spec humidity (kg/kg)
    call state_setexport(exportState, 'Si_qref' , input=Qref , lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Snow volume
    call state_setexport(exportState, 'Si_vsno' , input=vsno , lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Ice volume
    call state_setexport(exportState, 'Si_vice' , input=vice , lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Snow height
    tempfld(:,:,:) = c0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
          do i = ilo, ihi
             if ( tmask(i,j,iblk) .and. ailohi(i,j,iblk) > c0 ) then
                tempfld(i,j,iblk) = vsno(i,j,iblk)/ailohi(i,j,iblk)
             end if
          end do
       end do
    end do
    call state_setexport(exportState, 'Si_snowh' , input=tempfld , lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ------
    ! optional floe diameter and ice thickness to wave
    ! ------

    ! Sea ice thickness (m)
    if (State_FldChk(exportState, 'Si_thick')) then
       call state_setexport(exportState, 'Si_thick' , input=floethick , lmask=tmask, ifrac=ailohi, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Sea ice floe diameter (m)
    if (State_FldChk(exportState, 'Si_floediam')) then
       call state_setexport(exportState, 'Si_floediam' , input=floediam , lmask=tmask, ifrac=ailohi, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! ------
    ! ice/atm fluxes computed by ice
    ! ------

    ! Zonal air/ice stress
    call state_setexport(exportState, 'Faii_taux' , input=tauxa, lmask=tmask, ifrac=ailohi, &
         areacor=mod2med_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Meridional air/ice stress
    call state_setexport(exportState, 'Faii_tauy' , input=tauya, lmask=tmask, ifrac=ailohi, &
         areacor=mod2med_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Latent heat flux (atm into ice)
    call state_setexport(exportState, 'Faii_lat' , input=flat, lmask=tmask, ifrac=ailohi, &
         areacor=mod2med_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Sensible heat flux (atm into ice)
    call state_setexport(exportState, 'Faii_sen' , input=fsens, lmask=tmask, ifrac=ailohi, &
         areacor=mod2med_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Fix outgoing longwave if aice_init = 0, but aice > 0.
    tempfld(:,:,:) = flwout(:,:,:)
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
          do i = ilo, ihi
             if ( tmask(i,j,iblk) .and. ailohi(i,j,iblk) > c0 .and. flwout(i,j,iblk) > -puny) then
                 tempfld(i,j,iblk) = (-stefan_boltzmann *(Tf(i,j,iblk) + Tffresh)**4) / ailohi(i,j,iblk)
             end if
          end do
       end do
    end do
    ! longwave outgoing (upward), average over ice fraction only
    call state_setexport(exportState, 'Faii_lwup' , input=tempfld, lmask=tmask, ifrac=ailohi, &
         areacor=mod2med_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    deallocate(tempfld)

    ! Evaporative water flux (kg/m^2/s)
    call state_setexport(exportState, 'Faii_evap' , input=evap, lmask=tmask, ifrac=ailohi, &
         areacor=mod2med_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Shortwave flux absorbed in ice and ocean (W/m^2)
    call state_setexport(exportState, 'Faii_swnet' , input=fswabs, lmask=tmask, ifrac=ailohi, &
         areacor=mod2med_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ------
    ! ice/ocn fluxes computed by ice
    ! ------

    ! flux of shortwave through ice to ocean
    call state_setexport(exportState, 'Fioi_swpen' , input=fswthru, lmask=tmask, ifrac=ailohi, &
         areacor=mod2med_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (.not.prescribed_ice) then

    ! flux of vis dir shortwave through ice to ocean
    call state_setexport(exportState, 'Fioi_swpen_vdr' , input=fswthru_vdr, lmask=tmask, ifrac=ailohi, &
         areacor=mod2med_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! flux of vis dif shortwave through ice to ocean
    call state_setexport(exportState, 'Fioi_swpen_vdf' , input=fswthru_vdf, lmask=tmask, ifrac=ailohi, &
         areacor=mod2med_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! flux of ir dir shortwave through ice to ocean
    call state_setexport(exportState, 'Fioi_swpen_idr' , input=fswthru_idr, lmask=tmask, ifrac=ailohi, &
         areacor=mod2med_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! flux of ir dif shortwave through ice to ocean
    call state_setexport(exportState, 'Fioi_swpen_idf' , input=fswthru_idf, lmask=tmask, ifrac=ailohi, &
         areacor=mod2med_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    endif

    ! flux of heat exchange with ocean
    call state_setexport(exportState, 'Fioi_melth' , input=fhocn, lmask=tmask, ifrac=ailohi, &
         areacor=mod2med_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! flux fresh water to ocean (h2o flux from melting)
    call state_setexport(exportState, 'Fioi_meltw' , input=fresh, lmask=tmask, ifrac=ailohi, &
         areacor=mod2med_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! flux of salt to ocean (salt flux from melting)
    call state_setexport(exportState, 'Fioi_salt' , input=fsalt, lmask=tmask, ifrac=ailohi, &
         areacor=mod2med_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! stress n i/o zonal
    call state_setexport(exportState, 'Fioi_taux' , input=tauxo, lmask=tmask, ifrac=ailohi, &
         areacor=mod2med_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! stress n i/o meridional
    call state_setexport(exportState, 'Fioi_tauy' , input=tauyo, lmask=tmask, ifrac=ailohi, &
         areacor=mod2med_areacor, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ------
    ! optional aerosol fluxes to ocean
    ! ------

    ! hydrophobic bc
    if (State_FldChk(exportState, 'Fioi_bcpho')) then
       call state_setexport(exportState, 'Fioi_bcpho' , input=faero_ocn, index=1, lmask=tmask, ifrac=ailohi, &
            areacor=mod2med_areacor, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! hydrophilic bc
    if (State_FldChk(exportState, 'Fioi_bcphi')) then
       call state_setexport(exportState, 'Fioi_bcphi' , input=faero_ocn, index=2, lmask=tmask, ifrac=ailohi, &
            areacor=mod2med_areacor, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! dust
    if (State_FldChk(exportState, 'Fioi_flxdst')) then
       call state_setexport(exportState, 'Fioi_flxdst' , input=faero_ocn, index=3, lmask=tmask, ifrac=ailohi, &
            areacor=mod2med_areacor, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! ------
    ! optional water isotope fluxes to ocean
    ! ------

    if (State_FldChk(exportState, 'Fioi_meltw_wiso')) then
       ! 16O => ungridded_index=1
       ! 18O => ungridded_index=2
       ! HDO => ungridded_index=3

       call state_setexport(exportState, 'Fioi_meltw_wiso' , input=fiso_ocn, index=1, &
            lmask=tmask, ifrac=ailohi, ungridded_index=3, areacor=mod2med_areacor, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport(exportState, 'Fioi_meltw_wiso' , input=fiso_ocn, index=2, &
            lmask=tmask, ifrac=ailohi, ungridded_index=1, areacor=mod2med_areacor, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport(exportState, 'Fioi_meltw_wiso' , input=fiso_ocn, index=3, &
            lmask=tmask, ifrac=ailohi, ungridded_index=2, areacor=mod2med_areacor, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! ------
    ! optional water isotope fluxes to atmospehre
    ! ------

    if (State_FldChk(exportState, 'Faii_evap_wiso')) then
       !  Isotope evap to atm
       call state_setexport(exportState, 'Faii_evap_wiso' , input=fiso_evap, index=1, &
            lmask=tmask, ifrac=ailohi, ungridded_index=3, areacor=mod2med_areacor, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport(exportState, 'Faii_evap_wiso' , input=fiso_evap, index=2, &
            lmask=tmask, ifrac=ailohi, ungridded_index=1, areacor=mod2med_areacor, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport(exportState, 'Faii_evap_wiso' , input=fiso_evap, index=3, &
            lmask=tmask, ifrac=ailohi, ungridded_index=2, areacor=mod2med_areacor, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !  qref to atm
       call state_setexport(exportState, 'Si_qref_wiso' , input=Qref_iso, index=1, &
            lmask=tmask, ifrac=ailohi, ungridded_index=3, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport(exportState, 'Si_qref_wiso' , input=Qref_iso, index=2, &
            lmask=tmask, ifrac=ailohi, ungridded_index=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport(exportState, 'Si_qref_wiso' , input=Qref_iso, index=3, &
            lmask=tmask, ifrac=ailohi, ungridded_index=2, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    ! ------
    ! optional short wave penetration to ocean ice category
    ! ------

    ! ice fraction by category
    if ( State_FldChk(exportState, 'Si_ifrac_n') .and. &
         State_FldChk(exportState, 'Fioi_swpen_ifrac_n')) then
       do n = 1,ncat
          call state_setexport(exportState, 'Si_ifrac_n', input=aicen_init, index=n, &
               ungridded_index=n, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! penetrative shortwave by category
          ! Note: no need zero out pass-through fields over land for benefit of x2oacc fields in cpl hist files since
          ! the export state has been zeroed out at the beginning
          call state_setexport(exportState, 'Fioi_swpen_ifrac_n', input=fswthrun_ai, index=n, &
               lmask=tmask, ifrac=ailohi, ungridded_index=n, areacor=mod2med_areacor, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end do
    end if

  end subroutine ice_export

  !===============================================================================
  subroutine fldlist_add(num, fldlist, stdname, ungridded_lbound, ungridded_ubound)

    ! input/output variables
    integer             , intent(inout) :: num
    type(fld_list_type) , intent(inout) :: fldlist(:)
    character(len=*)    , intent(in)    :: stdname
    integer, optional   , intent(in)    :: ungridded_lbound
    integer, optional   , intent(in)    :: ungridded_ubound

    ! local variables
    character(len=*), parameter :: subname='(fldlist_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
       call abort_ice(trim(subname)//": ERROR num > fldsMax "//trim(stdname))
    endif
    fldlist(num)%stdname = trim(stdname)

    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
       fldlist(num)%ungridded_lbound = ungridded_lbound
       fldlist(num)%ungridded_ubound = ungridded_ubound
    end if

  end subroutine fldlist_add

  !===============================================================================
  subroutine fldlist_realize(state, fldList, numflds, flds_scalar_name, flds_scalar_num, mesh, grid, tag, rc)

    use NUOPC, only : NUOPC_IsConnected, NUOPC_Realize
    use ESMF , only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF , only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Mesh, ESMF_StateRemove
    use ESMF , only : ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LOGERR_PASSTHRU
    use ESMF , only : ESMF_VM

    ! input/output variables
    type(ESMF_State)          , intent(inout) :: state
    type(fld_list_type)       , intent(in)    :: fldList(:)
    integer                   , intent(in)    :: numflds
    character(len=*)          , intent(in)    :: flds_scalar_name
    integer                   , intent(in)    :: flds_scalar_num
    character(len=*)          , intent(in)    :: tag
    type(ESMF_Mesh), optional , intent(in)    :: mesh
    type(ESMF_Grid), optional , intent(in)    :: grid
    integer                   , intent(inout) :: rc

    ! local variables
    integer                :: n
    type(ESMF_Field)       :: field
    character(char_len)    :: stdname
    character(ESMF_MAXSTR) :: msg
    character(len=*),parameter  :: subname='(ice_import_export:fld_list_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO)
             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             if (present(mesh)) then
                ! Create the field
                if (fldlist(n)%ungridded_lbound > 0 .and. fldlist(n)%ungridded_ubound > 0) then
                   field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                        ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                        ungriddedUbound=(/fldlist(n)%ungridded_ubound/), &
                        gridToFieldMap=(/2/), rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   write(msg, '(a,i4,2x,i4)') trim(subname)//trim(tag)//" Field = "//trim(stdname)//&
                        " is connected using mesh with lbound, ubound = ",&
                        fldlist(n)%ungridded_lbound,fldlist(n)%ungridded_ubound
                   call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO)
                else
                   field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   write(msg, '(a,i4,a,i4)') trim(subname)//trim(tag)//" Field = "//trim(stdname)//&
                        " is connected using mesh without ungridded dimension"
                   call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO)
                end if
             else if (present(grid)) then
                call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using grid", &
                     ESMF_LOGMSG_INFO)
                if (fldlist(n)%ungridded_lbound > 0 .and. fldlist(n)%ungridded_ubound > 0) then
                   field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, name=stdname, indexflag=ESMF_INDEX_DELOCAL, &
                        ungriddedLBound=(/1,1/), ungriddedUBound=(/max_blocks,fldlist(n)%ungridded_ubound/), rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                else
                   field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, name=stdname, indexflag=ESMF_INDEX_DELOCAL, &
                        ungriddedLBound=(/1/), ungriddedUBound=(/max_blocks/), rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                end if
             else
                call ESMF_LogWrite(subname // 'input must be grid or mesh', ESMF_LOGMSG_INFO)
                rc = ESMF_FAILURE
                return
             end if
          end if ! if not scalar field

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO)
             call ESMF_StateRemove(state, (/stdname/), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if
    end do

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)
      ! ----------------------------------------------
      ! create a field with scalar data on the root pe
      ! ----------------------------------------------
      use ESMF, only : ESMF_Field, ESMF_DistGrid, ESMF_Grid
      use ESMF, only : ESMF_DistGridCreate, ESMF_GridCreate, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
      use ESMF, only : ESMF_FieldCreate, ESMF_GridCreate, ESMF_TYPEKIND_R8

      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      real(ESMF_KIND_R8), pointer :: fldptr2d(:,:)
      character(len=*), parameter :: subname='(ice_import_export:SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), gridToFieldMap=(/2/), rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      ! initialize fldptr to zero
      call ESMF_FieldGet(field, farrayPtr=fldptr2d, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      fldptr2d(:,:) = 0.0

    end subroutine SetScalarField

  end subroutine fldlist_realize

  !===============================================================================
  logical function State_FldChk(State, fldname)
    ! ----------------------------------------------
    ! Determine if field is in state
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State) , intent(in)  :: State
    character(len=*) , intent(in)  :: fldname

    ! local variables
    type(ESMF_StateItem_Flag) :: itemType
    ! ----------------------------------------------

    call ESMF_StateGet(State, trim(fldname), itemType)
    State_FldChk = (itemType /= ESMF_STATEITEM_NOTFOUND)

  end function State_FldChk

  !===============================================================================
  subroutine state_getimport_4d(state, fldname, output, index, ungridded_index, areacor, rc)

    ! ----------------------------------------------
    ! Map import state field to output array
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)              , intent(in)    :: state
    character(len=*)              , intent(in)    :: fldname
    real (kind=dbl_kind)          , intent(inout) :: output(:,:,:,:)
    integer                       , intent(in)    :: index
    integer, optional             , intent(in)    :: ungridded_index
    real(kind=dbl_kind), optional , intent(in)    :: areacor(:)
    integer                       , intent(out)   :: rc

    ! local variables
    type(block)                  :: this_block            ! block information for current block
    integer                      :: ilo, ihi, jlo, jhi    ! beginning and end of physical domain
    integer                      :: i, j, iblk, n, i1, j1 ! incides
    real(kind=dbl_kind), pointer :: dataPtr1d(:)          ! mesh
    real(kind=dbl_kind), pointer :: dataPtr2d(:,:)        ! mesh
    character(len=*), parameter  :: subname='(ice_import_export:state_getimport_4d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! check that fieldname exists
    if (.not. State_FldChk(state, trim(fldname))) return

    ! get field pointer
    if (present(ungridded_index)) then
       call state_getfldptr(state, trim(fldname), dataPtr2d, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call state_getfldptr(state, trim(fldname), dataPtr1d, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! set values of output array
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
             if (present(ungridded_index)) then
                output(i,j,index,iblk)  = dataPtr2d(ungridded_index,n)
             else
                output(i,j,index,iblk)  = dataPtr1d(n)
             end if
          end do
       end do
    end do
    if (present(areacor)) then
       n = 0
       do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo; ihi = this_block%ihi
          jlo = this_block%jlo; jhi = this_block%jhi
          do j = jlo, jhi
             do i = ilo, ihi
                n = n + 1
                output(i,j,index,iblk) = output(i,j,index,iblk) * areacor(n)
             end do
          end do
       end do
    end if

  end subroutine state_getimport_4d

  !===============================================================================
  subroutine state_getimport_3d(state, fldname, output, ungridded_index, areacor, rc)

    ! ----------------------------------------------
    ! Map import state field to output array
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)              , intent(in)    :: state
    character(len=*)              , intent(in)    :: fldname
    real (kind=dbl_kind)          , intent(inout) :: output(:,:,:)
    integer, optional             , intent(in)    :: ungridded_index
    real(kind=dbl_kind), optional , intent(in)    :: areacor(:)
    integer                       , intent(out)   :: rc

    ! local variables
    type(block)                  :: this_block            ! block information for current block
    integer                      :: ilo, ihi, jlo, jhi    ! beginning and end of physical domain
    integer                      :: i, j, iblk, n, i1, j1 ! incides
    real(kind=dbl_kind), pointer :: dataPtr1d(:)          ! mesh
    real(kind=dbl_kind), pointer :: dataPtr2d(:,:)        ! mesh
    character(len=*) , parameter :: subname='(ice_import_export:state_getimport_3d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! check that fieldname exists
    if (.not. State_FldChk(state, trim(fldname))) return

    ! get field pointer
    if (present(ungridded_index)) then
       call state_getfldptr(state, trim(fldname), dataPtr2d, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call state_getfldptr(state, trim(fldname), dataPtr1d, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! determine output array
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
             if (present(ungridded_index)) then
                output(i,j,iblk)  = dataPtr2d(ungridded_index,n)
             else
                output(i,j,iblk) = dataPtr1d(n)
             end if
          end do
       end do
    end do
    if (present(areacor)) then
       n = 0
       do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo; ihi = this_block%ihi
          jlo = this_block%jlo; jhi = this_block%jhi
          do j = jlo, jhi
             do i = ilo, ihi
                n = n + 1
                output(i,j,iblk) = output(i,j,iblk) * areacor(n)
             end do
          end do
       end do
    end if

  end subroutine state_getimport_3d

  !===============================================================================
  subroutine state_setexport_4d(state, fldname, input, index, lmask, ifrac, ungridded_index, areacor, rc)

    ! ----------------------------------------------
    ! Map 4d input array to export state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    ,           intent(inout) :: state
    character(len=*)    ,           intent(in)    :: fldname
    real(kind=dbl_kind) ,           intent(in)    :: input(:,:,:,:)
    integer             ,           intent(in)    :: index
    logical             , optional, intent(in)    :: lmask(:,:,:)
    real(kind=dbl_kind) , optional, intent(in)    :: ifrac(:,:,:)
    integer             , optional, intent(in)    :: ungridded_index
    real(kind=dbl_kind) , optional, intent(in)    :: areacor(:)
    integer             ,           intent(out)   :: rc

    ! local variables
    type(block)                  :: this_block            ! block information for current block
    integer                      :: ilo, ihi, jlo, jhi    ! beginning and end of physical domain
    integer                      :: i, j, iblk, n, i1, j1 ! indices
    real(kind=dbl_kind), pointer :: dataPtr1d(:)          ! mesh
    real(kind=dbl_kind), pointer :: dataPtr2d(:,:)        ! mesh
    integer                      :: ice_num
    character(len=*), parameter  :: subname='(ice_import_export:state_setexport_4d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! check that fieldname exists
    if (.not. State_FldChk(state, trim(fldname))) return

    ! get field pointer
    if (present(ungridded_index)) then
       call state_getfldptr(state, trim(fldname), dataPtr2d, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (ungridded_index == 1) then
          dataptr2d(:,:) = c0
       end if
       n = 0
       do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo; ihi = this_block%ihi
          jlo = this_block%jlo; jhi = this_block%jhi
          if (present(lmask) .and. present(ifrac)) then
             do j = jlo, jhi
                do i = ilo, ihi
                   n = n+1
                   if ( lmask(i,j,iblk) .and. ifrac(i,j,iblk) > c0 ) then
                      dataPtr2d(ungridded_index,n) = input(i,j,index,iblk)
                   else
                      dataPtr2d(ungridded_index,n) = c0
                   end if
                end do
             end do
          else
             do j = jlo, jhi
                do i = ilo, ihi
                   n = n+1
                   dataPtr2d(ungridded_index,n) = input(i,j,index,iblk)
                end do
             end do
          end if
       end do
       ice_num = n
       if (present(areacor)) then
          do n = 1,ice_num
             dataPtr2d(ungridded_index,n) = dataPtr2d(ungridded_index,n) * areacor(n)
          end do
       end if
    else
       call state_getfldptr(state, trim(fldname), dataPtr1d, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       dataptr1d(:) = c0
       n = 0
       do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo; ihi = this_block%ihi
          jlo = this_block%jlo; jhi = this_block%jhi
          if (present(lmask) .and. present(ifrac)) then
             do j = jlo, jhi
                do i = ilo, ihi
                   n = n+1
                   if ( lmask(i,j,iblk) .and. ifrac(i,j,iblk) > c0 ) then
                      dataPtr1d(n) = input(i,j,index,iblk)
                   end if
                end do
             end do
          else
             do j = jlo, jhi
                do i = ilo, ihi
                   n = n+1
                   dataPtr1d(n) = input(i,j,index,iblk)
                end do
             end do
          end if
       end do
       ice_num = n
       if (present(areacor)) then
          do n = 1,ice_num
             dataPtr1d(n) = dataPtr1d(n) * areacor(n)
          end do
       end if
    end if

  end subroutine state_setexport_4d

  !===============================================================================
  subroutine state_setexport_3d(state, fldname, input, lmask, ifrac, ungridded_index, areacor, rc)

    ! ----------------------------------------------
    ! Map 3d input array to export state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)               , intent(inout) :: state
    character(len=*)               , intent(in)    :: fldname
    real(kind=dbl_kind)            , intent(in)    :: input(:,:,:)
    logical             , optional , intent(in)    :: lmask(:,:,:)
    real(kind=dbl_kind) , optional , intent(in)    :: ifrac(:,:,:)
    integer             , optional , intent(in)    :: ungridded_index
    real(kind=dbl_kind) , optional , intent(in)    :: areacor(:)
    integer                        , intent(out)   :: rc

    ! local variables
    type(block)                  :: this_block            ! block information for current block
    integer                      :: ilo, ihi, jlo, jhi    ! beginning and end of physical domain
    integer                      :: i, j, iblk, n, i1, j1 ! incides
    real(kind=dbl_kind), pointer :: dataPtr1d(:)          ! mesh
    real(kind=dbl_kind), pointer :: dataPtr2d(:,:)        ! mesh
    integer                      :: num_ice
    character(len=*), parameter  :: subname='(ice_import_export:state_setexport_3d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! check that fieldname exists
    if (.not. State_FldChk(state, trim(fldname))) return

    ! get field pointer
    if (present(ungridded_index)) then
       call state_getfldptr(state, trim(fldname), dataPtr2d, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call state_getfldptr(state, trim(fldname), dataPtr1d, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

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
             if (present(lmask) .and. present(ifrac)) then
                if ( lmask(i,j,iblk) .and. ifrac(i,j,iblk) > c0 ) then
                   if (present(ungridded_index)) then
                      dataPtr2d(ungridded_index,n) = input(i,j,iblk)
                   else
                      dataPtr1d(n) = input(i,j,iblk)
                   end if
                end if
             else
                if (present(ungridded_index)) then
                   dataPtr2d(ungridded_index,n) = input(i,j,iblk)
                else
                   dataPtr1d(n) = input(i,j,iblk)
                end if
             end if
          end do
       end do
    end do
    num_ice = n
    if (present(areacor)) then
       if (present(ungridded_index)) then
          do n = 1,num_ice
             dataPtr2d(:,n) = dataPtr2d(:,n) * areacor(n)
          end do
       else
          do n = 1,num_ice
             dataPtr1d(n) = dataPtr1d(n) * areacor(n)
          end do
       end if
    end if

  end subroutine state_setexport_3d

  !===============================================================================
  subroutine State_GetFldPtr_1d(State, fldname, fldptr, rc)
    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)              , intent(in)     :: State
    character(len=*)              , intent(in)     :: fldname
    real(kind=dbl_kind) , pointer , intent(inout)  :: fldptr(:)
    integer, optional             , intent(out)    :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    character(len=*),parameter :: subname='(ice_import_export:State_GetFldPtr_1d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine State_GetFldPtr_1d

  !===============================================================================
  subroutine State_GetFldPtr_2d(State, fldname, fldptr, rc)
    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    ,            intent(in)     :: State
    character(len=*)    ,            intent(in)     :: fldname
    real(kind=dbl_kind) , pointer ,  intent(inout)  :: fldptr(:,:)
    integer             , optional , intent(out)    :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    character(len=*),parameter :: subname='(ice_import_export:State_GetFldPtr_2d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine State_GetFldPtr_2d

end module ice_import_export
