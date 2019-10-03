module ice_import_export

  use ESMF
  use NUOPC
  use NUOPC_Model
  use shr_sys_mod        , only : shr_sys_abort, shr_sys_flush
  use shr_frz_mod        , only : shr_frz_freezetemp
  use shr_kind_mod       , only : r8 => shr_kind_r8, cl=>shr_kind_cl, cs=>shr_kind_cs 
  use ice_kinds_mod      , only : int_kind, dbl_kind, char_len_long, log_kind
  use ice_constants      , only : c0, c1, tffresh, spval_dbl
  use ice_constants      , only : field_loc_center, field_type_scalar, field_type_vector
  use ice_blocks         , only : block, get_block, nx_block, ny_block
  use ice_domain         , only : nblocks, blocks_ice, halo_info, distrb_info
  use ice_domain_size    , only : nx_global, ny_global, block_size_x, block_size_y, max_blocks, ncat
  use ice_flux           , only : strairxt, strairyt, strocnxt, strocnyt
  use ice_flux           , only : alvdr, alidr, alvdf, alidf, Tref, Qref, Uref
  use ice_flux           , only : flat, fsens, flwout, evap, fswabs, fhocn, fswthru
  use ice_flux           , only : fswthruvdr, fswthruvdf, fswthruidr, fswthruidf
  use ice_flux           , only : send_i2x_per_cat, fswthrun_ai
  use ice_flux           , only : fresh, fsalt, zlvl, uatm, vatm, potT, Tair, Qa
  use ice_flux           , only : rhoa, swvdr, swvdf, swidr, swidf, flw, frain
  use ice_flux           , only : fsnow, uocn, vocn, sst, ss_tltx, ss_tlty, frzmlt
  use ice_flux           , only : sss, tf, wind, fsw
  use ice_flux           , only : faero_atm, faero_ocn
  use ice_flux           , only : fiso_atm, fiso_ocn, fiso_rain, fiso_evap
  use ice_flux           , only : Qa_iso, Qref_iso, HDO_ocn, H2_18O_ocn, H2_16O_ocn
  use ice_state          , only : vice, vsno, aice, aicen_init, trcr
  use ice_grid           , only : tlon, tlat, tarea, tmask, anglet, hm, ocn_gridcell_frac
  use ice_grid           , only : grid_type, t2ugrid_vector
  use ice_boundary       , only : ice_HaloUpdate
  use ice_fileunits      , only : nu_diag
  use ice_communicate    , only : my_task, master_task, MPI_COMM_ICE
  use ice_prescribed_mod , only : prescribed_ice
  use ice_shr_methods    , only : chkerr, state_reset
  use perf_mod           , only : t_startf, t_stopf, t_barrierf

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
     module procedure state_getfldptr_3d
     module procedure state_getfldptr_4d
  end interface state_getfldptr
  private :: state_getfldptr

  interface state_getimport
     module procedure state_getimport_4d_output
     module procedure state_getimport_3d_output
  end interface state_getimport
  private :: state_getimport

  interface state_setexport
     module procedure state_setexport_4d_input
     module procedure state_setexport_3d_input
  end interface state_setexport
  private :: state_setexport

  ! Private module data

  type fld_list_type
    character(len=128) :: stdname
    integer :: ungridded_lbound = 0
    integer :: ungridded_ubound = 0
  end type fld_list_type

  integer, parameter       :: fldsMax = 100
  integer                  :: fldsToIce_num = 0
  integer                  :: fldsFrIce_num = 0
  type (fld_list_type)     :: fldsToIce(fldsMax)
  type (fld_list_type)     :: fldsFrIce(fldsMax)
  type(ESMF_GeomType_Flag) :: geomtype

  integer     , parameter  :: dbug = 10        ! i/o debug messages
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
    integer       :: n
    character(CS) :: stdname
    character(CS) :: cvalue
    logical       :: flds_wiso         ! use case
    logical       :: flds_i2o_per_cat  ! .true. => select per ice thickness category
    character(len=*), parameter :: subname='(ice_import_export:ice_advertise_fields)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call NUOPC_CompAttributeGet(gcomp, name='flds_wiso', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_wiso
    call ESMF_LogWrite('flds_wiso = '// trim(cvalue), ESMF_LOGMSG_INFO)

    call NUOPC_CompAttributeGet(gcomp, name='flds_i2o_per_cat', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) send_i2x_per_cat
    call ESMF_LogWrite('flds_i2o_per_cat = '// trim(cvalue), ESMF_LOGMSG_INFO)

    !-----------------
    ! advertise import fields
    !-----------------

    call fldlist_add(fldsToIce_num, fldsToIce, trim(flds_scalar_name))

    ! from ocean
    call fldlist_add(fldsToIce_num, fldsToIce, 'sea_surface_slope_zonal' )
    call fldlist_add(fldsToIce_num, fldsToIce, 'sea_surface_slope_merid' )
    call fldlist_add(fldsToIce_num, fldsToIce, 'sea_surface_temperature' )
    call fldlist_add(fldsToIce_num, fldsToIce, 's_surf'                  )
    call fldlist_add(fldsToIce_num, fldsToIce, 'ocn_current_zonal'       )
    call fldlist_add(fldsToIce_num, fldsToIce, 'ocn_current_merid'       )
    call fldlist_add(fldsToIce_num, fldsToIce, 'freezing_melting_potential'                  )
    if (flds_wiso) then
       call fldlist_add(fldsToIce_num, fldsToIce, 'So_roce_wiso', ungridded_lbound=1, ungridded_ubound=3)
    end if

    ! from atmosphere
    call fldlist_add(fldsToIce_num, fldsToIce, 'inst_height_lowest'            )
    call fldlist_add(fldsToIce_num, fldsToIce, 'inst_zonal_wind_height_lowest' )
    call fldlist_add(fldsToIce_num, fldsToIce, 'inst_merid_wind_height_lowest' )
    call fldlist_add(fldsToIce_num, fldsToIce, 'inst_spec_humid_height_lowest' )
    call fldlist_add(fldsToIce_num, fldsToIce, 'inst_temp_height_lowest'       )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_ptem'                       )
    call fldlist_add(fldsToIce_num, fldsToIce, 'air_density_height_lowest'                       )
    call fldlist_add(fldsToIce_num, fldsToIce, 'mean_down_sw_vis_dir_flx'      )
    call fldlist_add(fldsToIce_num, fldsToIce, 'mean_down_sw_ir_dir_flx'       )
    call fldlist_add(fldsToIce_num, fldsToIce, 'mean_down_sw_vis_dif_flx'      )
    call fldlist_add(fldsToIce_num, fldsToIce, 'mean_down_sw_ir_dif_flx'       )
    call fldlist_add(fldsToIce_num, fldsToIce, 'mean_down_lw_flx'              )
    call fldlist_add(fldsToIce_num, fldsToIce, 'mean_prec_rate'                )
    call fldlist_add(fldsToIce_num, fldsToIce, 'mean_fprec_rate'               )

    ! from atm - black carbon deposition fluxes (3)
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_bcph',  ungridded_lbound=1, ungridded_ubound=3)

    ! from atm - wet dust deposition frluxes (4 sizes)
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_dstwet', ungridded_lbound=1, ungridded_ubound=4)

    ! from - atm dry dust deposition frluxes (4 sizes)
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_dstdry', ungridded_lbound=1, ungridded_ubound=4)

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
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'ice_mask'                )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'ice_fraction'            )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'sea_ice_temperature'     )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'mean_ice_volume'         )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'mean_snow_volume'        )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_tref'                 )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_qref'                 )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_snowh'                )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_u10'                  )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'inst_ice_vis_dir_albedo' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'inst_ice_ir_dir_albedo'  )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'inst_ice_vis_dif_albedo' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'inst_ice_ir_dif_albedo'  )
    if (send_i2x_per_cat) then
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'ice_fraction_n', &
            ungridded_lbound=1, ungridded_ubound=ncat)
    end if

    ! ice/atm fluxes computed by ice
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'stress_on_air_ice_zonal'          )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'stress_on_air_ice_merid'          )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'mean_laten_heat_flx_atm_into_ice' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'mean_sensi_heat_flx_atm_into_ice' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'mean_up_lw_flx_ice'               )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'mean_evap_rate_atm_into_ice'      )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_swnet'                       )

    ! ice/ocn fluxes computed by ice
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'net_heat_flx_to_ocn'     )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'mean_sw_pen_to_ocn'      )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'mean_sw_pen_to_ocn_vis_dir_flx' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'mean_sw_pen_to_ocn_vis_dif_flx' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'mean_sw_pen_to_ocn_ir_dir_flx'  )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'mean_sw_pen_to_ocn_ir_dif_flx'  )
    if (send_i2x_per_cat) then
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'mean_sw_pen_to_ocn_ifrac_n', &
            ungridded_lbound=1, ungridded_ubound=ncat)
    end if
    call fldlist_add(fldsFrIce_num , fldsFrIce, 'mean_fresh_water_to_ocean_rate' )
    call fldlist_add(fldsFrIce_num , fldsFrIce, 'mean_salt_rate'                 )
    call fldlist_add(fldsFrIce_num , fldsFrIce, 'stress_on_ocn_ice_zonal'        )
    call fldlist_add(fldsFrIce_num , fldsFrIce, 'stress_on_ocn_ice_merid'        )
    call fldlist_add(fldsFrIce_num , fldsFrIce, 'Fioi_bcpho'                     )
    call fldlist_add(fldsFrIce_num , fldsFrIce, 'Fioi_bcphi'                     )
    call fldlist_add(fldsFrIce_num , fldsFrIce, 'Fioi_flxdst'                    )
    if (flds_wiso) then
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'mean_fresh_water_to_ocean_rate_wiso', &
            ungridded_lbound=1, ungridded_ubound=3)
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'mean_evap_rate_atm_into_ice_wiso', &
            ungridded_lbound=1, ungridded_ubound=3)
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_qref_wiso', &
            ungridded_lbound=1, ungridded_ubound=3)
    end if

    do n = 1,fldsFrIce_num
       call NUOPC_Advertise(exportState, standardName=fldsFrIce(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ice_advertise_fields

!==============================================================================

  subroutine ice_realize_fields(gcomp, mesh, grid, flds_scalar_name, flds_scalar_num, rc)

    ! input/output variables
    type(ESMF_GridComp)                      :: gcomp
    type(ESMF_Mesh) , optional , intent(in)  :: mesh
    type(ESMF_Grid) , optional , intent(in)  :: grid
    character(len=*)           , intent(in)  :: flds_scalar_name
    integer                    , intent(in)  :: flds_scalar_num 
    integer                    , intent(out) :: rc

    ! local variables
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    character(len=*), parameter :: subname='(ice_import_export:realize_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (present(mesh)) then

       geomtype = ESMF_GEOMTYPE_MESH

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

    else if (present(grid)) then

       geomtype = ESMF_GEOMTYPE_GRID

       call fldlist_realize( &
            state=ExportState, &
            fldList=fldsFrIce, &
            numflds=fldsFrIce_num, &
            flds_scalar_name=flds_scalar_name, &
            flds_scalar_num=flds_scalar_num, &
            tag=subname//':CICE_Export',&
            grid=grid, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call fldlist_realize( &
            state=importState, &
            fldList=fldsToIce, &
            numflds=fldsToIce_num, &
            flds_scalar_name=flds_scalar_name, &
            flds_scalar_num=flds_scalar_num, &
            tag=subname//':CICE_Import',&
            grid=grid, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

  end subroutine ice_realize_fields

  !==============================================================================

  subroutine ice_import( importState, rc )

    ! input/output variables
    type(ESMF_State) , intent(in)  :: importState
    integer          , intent(out) :: rc

    ! local variables
    integer,parameter                :: nflds=15
    integer,parameter                :: nfldv=6
    integer                          :: i, j, iblk, n
    integer                          :: ilo, ihi, jlo, jhi !beginning and end of physical domain
    type(block)                      :: this_block         ! block information for current block
    real (kind=dbl_kind),allocatable :: aflds(:,:,:,:)
    real (kind=dbl_kind)             :: workx, worky
    real (kind=dbl_kind)             :: MIN_RAIN_TEMP, MAX_SNOW_TEMP
    character(len=*),   parameter    :: subname = 'ice_import'
    !-----------------------------------------------------

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

    call state_getimport(importState, 'sea_surface_temperature', output=aflds, index=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 's_surf', output=aflds, index=2, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! import ocean states

    call state_getimport(importState, 'inst_height_lowest', output=aflds, index=3, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Sa_ptem', output=aflds, index=4, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'inst_temp_height_lowest', output=aflds, index=5, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'inst_spec_humid_height_lowest', output=aflds, index=6, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'air_density_height_lowest', output=aflds, index=7, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! import ocn/ice fluxes

    call state_getimport(importState, 'freezing_melting_potential', output=aflds, index=8, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! import atm fluxes

    call state_getimport(importState, 'mean_down_sw_vis_dir_flx', output=aflds, index=9, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'mean_down_sw_ir_dir_flx', output=aflds, index=10, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'mean_down_sw_vis_dif_flx', output=aflds, index=11, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'mean_down_sw_ir_dif_flx', output=aflds, index=12, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'mean_down_lw_flx', output=aflds, index=13, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'mean_prec_rate', output=aflds, index=14, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'mean_fprec_rate', output=aflds, index=15, rc=rc)
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
             sst  (i,j,iblk)   = aflds(i,j, 1,iblk)
             sss  (i,j,iblk)   = aflds(i,j, 2,iblk)
             zlvl (i,j,iblk)   = aflds(i,j, 3,iblk)
             potT (i,j,iblk)   = aflds(i,j, 4,iblk)
             Tair (i,j,iblk)   = aflds(i,j, 5,iblk)
             Qa   (i,j,iblk)   = aflds(i,j, 6,iblk)
             rhoa (i,j,iblk)   = aflds(i,j, 7,iblk)
             frzmlt (i,j,iblk) = aflds(i,j, 8,iblk)
             swvdr(i,j,iblk)   = aflds(i,j, 9,iblk)
             swidr(i,j,iblk)   = aflds(i,j,10,iblk)
             swvdf(i,j,iblk)   = aflds(i,j,11,iblk)
             swidf(i,j,iblk)   = aflds(i,j,12,iblk)
             flw  (i,j,iblk)   = aflds(i,j,13,iblk)
             frain(i,j,iblk)   = aflds(i,j,14,iblk)
             fsnow(i,j,iblk)   = aflds(i,j,15,iblk)
          enddo    !i
       enddo    !j
    enddo        !iblk
    !$OMP END PARALLEL DO

    deallocate(aflds)
    allocate(aflds(nx_block,ny_block,nfldv,nblocks))
    aflds = c0

    ! Get velocity fields from ocean and atm and slope fields from ocean

    call state_getimport(importState, 'ocn_current_zonal', output=aflds, index=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'ocn_current_merid', output=aflds, index=2, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'inst_zonal_wind_height_lowest', output=aflds, index=3, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'inst_merid_wind_height_lowest', output=aflds, index=4, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'sea_surface_slope_zonal', output=aflds, index=5, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'sea_surface_slope_merid', output=aflds, index=6, rc=rc)
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

       ! bcphodry
       call state_getimport(importState, 'Faxa_bcph', output=faero_atm, index=1,                ungridded_index=2, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! bcphidry + bcphiwet
       call state_getimport(importState, 'Faxa_bcph', output=faero_atm, index=2,                ungridded_index=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'Faxa_bcph', output=faero_atm, index=2, do_sum=.true., ungridded_index=3, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Sum over all dry and wet dust fluxes from ath atmosphere
    if (State_FldChk(importState, 'Faxa_dstwet') .and. State_FldChk(importState, 'Faxa_dstdry')) then
       call state_getimport(importState, 'Faxa_dstwet', output=faero_atm,  index=3,                ungridded_index=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'Faxa_dstdry', output=faero_atm,  index=3, do_sum=.true., ungridded_index=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'Faxa_dstwet', output=faero_atm,  index=3, do_sum=.true., ungridded_index=2, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'Faxa_dstdry', output=faero_atm,  index=3, do_sum=.true., ungridded_index=2, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'Faxa_dstwet', output=faero_atm,  index=3, do_sum=.true., ungridded_index=3, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'Faxa_dstdry', output=faero_atm,  index=3, do_sum=.true., ungridded_index=3, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'Faxa_dstwet', output=faero_atm,  index=3, do_sum=.true., ungridded_index=4, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'Faxa_dstdry', output=faero_atm,  index=3, do_sum=.true., ungridded_index=4, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !-------------------------------------------------------
    ! Water isotopes from the mediator
    !-------------------------------------------------------

    ! 16O => ungridded_index=1
    ! 18O => ungridded_index=2
    ! HDO => ungridded_index=3

    if (State_FldChk(importState, 'shum_wiso')) then
       call state_getimport(importState, 'inst_spec_humid_height_lowest_wiso', output=Qa_iso, index=1, ungridded_index=3, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'inst_spec_humid_height_lowest_wiso', output=Qa_iso, index=2, ungridded_index=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'inst_spec_humid_height_lowest_wiso', output=Qa_iso, index=3, ungridded_index=2, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call state_getimport(importState, 'mean_prec_rate_wiso', output=fiso_rain, index=1, ungridded_index=3, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'mean_prec_rate_wiso', output=fiso_rain, index=2, ungridded_index=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'mean_prec_rate_wiso', output=fiso_rain, index=3, ungridded_index=2, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call state_getimport(importState, 'mean_fprec_rate_wiso', output=fiso_atm, index=1, ungridded_index=3, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'mean_fprec_rate_wiso', output=fiso_atm, index=2, ungridded_index=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'mean_fprec_rate_wiso', output=fiso_atm, index=3, ungridded_index=2, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call state_getimport(importState, 'So_roce_wiso', output=HDO_ocn   , ungridded_index=3, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'So_roce_wiso', output=H2_16O_ocn, ungridded_index=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getimport(importState, 'So_roce_wiso', output=H2_18O_ocn, ungridded_index=2, rc=rc)
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
             uocn(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid
                            + worky*sin(ANGLET(i,j,iblk))
             vocn(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                            - workx*sin(ANGLET(i,j,iblk))

             workx      = ss_tltx  (i,j,iblk)           ! sea sfc tilt, m/m
             worky      = ss_tlty  (i,j,iblk)
             ss_tltx(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid
                               + worky*sin(ANGLET(i,j,iblk))
             ss_tlty(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                               - workx*sin(ANGLET(i,j,iblk))

             sst(i,j,iblk) = sst(i,j,iblk) - Tffresh       ! sea sfc temp (C)

             sss(i,j,iblk) = max(sss(i,j,iblk),c0)
          enddo
       enddo

       ! Use shr_frz_mod for this
       Tf(:,:,iblk) = shr_frz_freezetemp(sss(:,:,iblk))

    enddo
    !$OMP END PARALLEL DO
    call t_stopf ('cice_imp_ocn')

    ! Interpolate ocean dynamics variables from T-cell centers to
    ! U-cell centers.

    if (.not.prescribed_ice) then
       call t_startf ('cice_imp_t2u')
       call t2ugrid_vector(uocn)
       call t2ugrid_vector(vocn)
       call t2ugrid_vector(ss_tltx)
       call t2ugrid_vector(ss_tlty)
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

    ! input/output variables
    type(ESMF_State), intent(inout) :: exportState
    integer         , intent(out)   :: rc

    ! local variables
    type(block)             :: this_block                           ! block information for current block
    integer                 :: i, j, iblk, n                        ! incides
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
    real    (kind=dbl_kind), allocatable :: tempfld(:,:,:)
    character(len=*),parameter :: subname = 'ice_export'
    !-----------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !calculate ice thickness from aice and vice. Also
    !create Tsrf from the first tracer (trcr) in ice_state.F

    ailohi(:,:,:) = c0
    Tsrf(:,:,:)  = c0
    tauxa(:,:,:) = c0
    tauya(:,:,:) = c0
    tauxo(:,:,:) = c0
    tauyo(:,:,:) = c0

    !$OMP PARALLEL DO PRIVATE(iblk,i,j,workx,worky, this_block, ilo, ihi, jlo, jhi)
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

             ! wind stress  (on POP T-grid:  convert to lat-lon)
             workx = strairxT(i,j,iblk)                             ! N/m^2
             worky = strairyT(i,j,iblk)                             ! N/m^2
             tauxa(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) - worky*sin(ANGLET(i,j,iblk))
             tauya(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) + workx*sin(ANGLET(i,j,iblk))

             ! ice/ocean stress (on POP T-grid:  convert to lat-lon)
             workx = -strocnxT(i,j,iblk)                            ! N/m^2
             worky = -strocnyT(i,j,iblk)                            ! N/m^2
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
                   call shr_sys_flush(nu_diag)
                endif
             end do
          end do
       end do
    endif

    !---------------------------------
    ! Create the export state
    !---------------------------------

    ! Zero out fields with tmask for proper coupler accumulation in ice free areas
    call state_reset(exportState, c0, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Create a temporary field
    allocate(tempfld(nx_block,ny_block,nblocks))

    ! Fractions and mask
    call state_setexport(exportState, 'ice_fraction', input=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (trim(grid_type) == 'latlon') then
       call state_setexport(exportState, 'ice_mask', input=ocn_gridcell_frac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
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
       call state_setexport(exportState, 'ice_mask', input=tempfld, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! ----
    ! States from ice
    ! ----

    ! surface temperature of ice covered portion (degK)
    call state_setexport(exportState, 'sea_ice_temperature', input=Tsrf , lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! albedo vis dir
    call state_setexport(exportState, 'inst_ice_vis_dir_albedo', input=alvdr, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! albedo nir dir
    call state_setexport(exportState, 'inst_ice_ir_dir_albedo', input=alidr, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! albedo vis dif
    call state_setexport(exportState, 'inst_ice_vis_dif_albedo', input=alvdf, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! albedo nir dif
    call state_setexport(exportState, 'inst_ice_ir_dif_albedo', input=alidf, lmask=tmask, ifrac=ailohi, rc=rc)
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
    call state_setexport(exportState, 'mean_snow_volume' , input=vsno , lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Ice volume
    call state_setexport(exportState, 'mean_ice_volume' , input=vice , lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Snow height
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
    ! ice/atm fluxes computed by ice
    ! ------

    ! Zonal air/ice stress
    call state_setexport(exportState, 'stress_on_air_ice_zonal' , input=tauxa, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Meridional air/ice stress
    call state_setexport(exportState, 'stress_on_air_ice_merid' , input=tauya, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Latent heat flux (atm into ice)
    call state_setexport(exportState, 'mean_laten_heat_flx_atm_into_ice' , input=flat, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Sensible heat flux (atm into ice)
    call state_setexport(exportState, 'mean_sensi_heat_flx_atm_into_ice' , input=fsens, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! longwave outgoing (upward), average over ice fraction only
    call state_setexport(exportState, 'mean_up_lw_flx_ice' , input=flwout, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Evaporative water flux (kg/m^2/s)
    call state_setexport(exportState, 'mean_evap_rate_atm_into_ice' , input=evap, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Shortwave flux absorbed in ice and ocean (W/m^2)
    call state_setexport(exportState, 'Faii_swnet' , input=fswabs, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ------
    ! ice/ocn fluxes computed by ice
    ! ------

    ! flux of shortwave through ice to ocean
    call state_setexport(exportState, 'mean_sw_pen_to_ocn' , input=fswthru, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! flux of vis dir shortwave through ice to ocean
    call state_setexport(exportState, 'mean_sw_pen_to_ocn_vis_dir_flx' , input=fswthruvdr, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! flux of vis dif shortwave through ice to ocean
    call state_setexport(exportState, 'mean_sw_pen_to_ocn_vis_dif_flx' , input=fswthruvdf, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! flux of ir dir shortwave through ice to ocean
    call state_setexport(exportState, 'mean_sw_pen_to_ocn_ir_dir_flx' , input=fswthruidr, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! flux of ir dif shortwave through ice to ocean
    call state_setexport(exportState, 'mean_sw_pen_to_ocn_ir_dif_flx' , input=fswthruidf, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! heat exchange with ocean
    call state_setexport(exportState, 'net_heat_flx_to_ocn' , input=fhocn, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! fresh water to ocean (h2o flux from melting)
    call state_setexport(exportState, 'mean_fresh_water_to_ocean_rate' , input=fresh, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! salt to ocean (salt flux from melting)
    call state_setexport(exportState, 'mean_salt_rate' , input=fsalt, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! stress n i/o zonal
    call state_setexport(exportState, 'stress_on_ocn_ice_zonal' , input=tauxo, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! stress n i/o meridional
    call state_setexport(exportState, 'stress_on_ocn_ice_merid' , input=tauyo, lmask=tmask, ifrac=ailohi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ------
    ! optional aerosol fluxes to ocean
    ! ------

    ! hydrophobic bc
    if (State_FldChk(exportState, 'Fioi_bcpho')) then
       call state_setexport(exportState, 'Fioi_bcpho' , input=faero_ocn, index=1, lmask=tmask, ifrac=ailohi, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! hydrophilic bc
    if (State_FldChk(exportState, 'Fioi_bcphi')) then
       call state_setexport(exportState, 'Fioi_bcphi' , input=faero_ocn, index=2, lmask=tmask, ifrac=ailohi, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! dust
    if (State_FldChk(exportState, 'Fioi_flxdst')) then
       call state_setexport(exportState, 'Fioi_flxdst' , input=faero_ocn, index=3, lmask=tmask, ifrac=ailohi, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! ------
    ! optional water isotope fluxes to ocean
    ! ------

    if (State_FldChk(exportState, 'mean_fresh_water_to_ocean_rate_wiso')) then
       ! 16O => ungridded_index=1
       ! 18O => ungridded_index=2
       ! HDO => ungridded_index=3

       call state_setexport(exportState, 'mean_fresh_water_to_ocean_rate_wiso' , input=fiso_ocn, index=1, &
            lmask=tmask, ifrac=ailohi, ungridded_index=3, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport(exportState, 'mean_fresh_water_to_ocean_rate_wiso' , input=fiso_ocn, index=2, &
            lmask=tmask, ifrac=ailohi, ungridded_index=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport(exportState, 'mean_fresh_water_to_ocean_rate_wiso' , input=fiso_ocn, index=3, &
            lmask=tmask, ifrac=ailohi, ungridded_index=2, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! ------
    ! optional water isotope fluxes to atmospehre
    ! ------

    if (State_FldChk(exportState, 'mean_evap_rate_atm_into_ice_wiso')) then
       !  Isotope evap to atm
       call state_setexport(exportState, 'mean_evap_rate_atm_into_ice_wiso' , input=fiso_evap, index=1, &
            lmask=tmask, ifrac=ailohi, ungridded_index=3, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport(exportState, 'mean_evap_rate_atm_into_ice_wiso' , input=fiso_evap, index=2, &
            lmask=tmask, ifrac=ailohi, ungridded_index=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport(exportState, 'mean_evap_rate_atm_into_ice_wiso' , input=fiso_evap, index=3, &
            lmask=tmask, ifrac=ailohi, ungridded_index=2, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !  Isotope evap to atm
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
    if ( State_FldChk(exportState, 'ice_fraction_n') .and. &
         State_FldChk(exportState, 'mean_sw_pen_to_ocn_ifrac_n')) then
       do n = 1,ncat
          call state_setexport(exportState, 'ice_fraction_n', input=aicen_init, index=n, &
               ungridded_index=n, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! penetrative shortwave by category
          ! Note: no need zero out pass-through fields over land for benefit of x2oacc fields in cpl hist files since
          ! the export state has been zeroed out at the beginning
          call state_setexport(exportState, 'mean_sw_pen_to_ocn_ifrac_n', input=aicen_init, index=n, &
               lmask=tmask, ifrac=ailohi, ungridded_index=n, rc=rc)
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
       call shr_sys_abort(trim(subname)//": ERROR num > fldsMax "//trim(stdname))
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
    character(len=80)      :: stdname
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
                call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                     ESMF_LOGMSG_INFO)
                ! Create the field
                if (fldlist(n)%ungridded_lbound > 0 .and. fldlist(n)%ungridded_ubound > 0) then
                   field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                        ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                        ungriddedUbound=(/fldlist(n)%ungridded_ubound/), &
                        gridToFieldMap=(/2/), rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                else
                   field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
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

  subroutine state_getimport_4d_output(state, fldname, output, index, do_sum, ungridded_index, rc)

    ! ----------------------------------------------
    ! Map import state field to output array
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)           , intent(in)    :: state
    character(len=*)           , intent(in)    :: fldname
    real (kind=dbl_kind)       , intent(inout) :: output(:,:,:,:)
    integer                    , intent(in)    :: index
    logical, optional          , intent(in)    :: do_sum
    integer, optional          , intent(in)    :: ungridded_index
    integer                    , intent(out)   :: rc

    ! local variables
    type(block)                  :: this_block            ! block information for current block
    integer                      :: ilo, ihi, jlo, jhi    ! beginning and end of physical domain
    integer                      :: i, j, iblk, n, i1, j1 ! incides
    real(kind=dbl_kind), pointer :: dataPtr1d(:)          ! mesh
    real(kind=dbl_kind), pointer :: dataPtr2d(:,:)        ! mesh
    real(kind=dbl_kind), pointer :: dataPtr3d(:,:,:)      ! grid
    real(kind=dbl_kind), pointer :: dataPtr4d(:,:,:,:)    ! grid
    character(len=*), parameter  :: subname='(ice_import_export:state_getimport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (geomtype == ESMF_GEOMTYPE_MESH) then

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
                if (present(do_sum)) then ! do sum
                   if (present(ungridded_index)) then
                      output(i,j,index,iblk)  = output(i,j,index, iblk) + dataPtr2d(ungridded_index,n)
                   else
                      output(i,j,index,iblk)  = output(i,j,index, iblk) + dataPtr1d(n)
                   end if
                else ! do not do sum
                   if (present(ungridded_index)) then
                      output(i,j,index,iblk)  = dataPtr2d(ungridded_index,n)
                   else
                      output(i,j,index,iblk)  = dataPtr1d(n)
                   end if
                end if
             end do
          end do
       end do

    else if (geomtype == ESMF_GEOMTYPE_GRID) then

       ! get field pointer
       if (present(ungridded_index)) then
          call state_getfldptr(state, trim(fldname), dataptr4d, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call state_getfldptr(state, trim(fldname), dataptr3d, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! set values of output array
       do iblk = 1,nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi
          do j = jlo,jhi
             do i = ilo,ihi
                i1 = i - ilo + 1
                j1 = j - jlo + 1
                if (present(do_sum)) then
                   if (present(ungridded_index)) then
                      output(i,j,index,iblk) = output(i,j,index,iblk) + dataPtr4d(i1,j1,iblk,ungridded_index)
                   else
                      output(i,j,index,iblk) = output(i,j,index,iblk) + dataPtr3d(i1,j1,iblk)
                   end if
                else
                   if (present(ungridded_index)) then
                      output(i,j,index,iblk) = dataPtr4d(i1,j1,iblk,ungridded_index)
                   else
                      output(i,j,index,iblk) = dataPtr3d(i1,j1,iblk)
                   end if
                end if
             end do
          end do
       end do

    end if

  end subroutine state_getimport_4d_output

  !===============================================================================

  subroutine state_getimport_3d_output(state, fldname, output, do_sum, ungridded_index, rc)

    ! ----------------------------------------------
    ! Map import state field to output array
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)           , intent(in)    :: state
    character(len=*)           , intent(in)    :: fldname
    real (kind=dbl_kind)       , intent(inout) :: output(:,:,:)
    logical, optional          , intent(in)    :: do_sum
    integer, optional          , intent(in)    :: ungridded_index
    integer                    , intent(out)   :: rc

    ! local variables
    type(block)                  :: this_block            ! block information for current block
    integer                      :: ilo, ihi, jlo, jhi    ! beginning and end of physical domain
    integer                      :: i, j, iblk, n, i1, j1 ! incides
    real(kind=dbl_kind), pointer :: dataPtr1d(:)          ! mesh
    real(kind=dbl_kind), pointer :: dataPtr2d(:,:)        ! mesh
    real(kind=dbl_kind), pointer :: dataPtr3d(:,:,:)      ! grid
    real(kind=dbl_kind), pointer :: dataPtr4d(:,:,:,:)    ! grid
    character(len=*) , parameter :: subname='(ice_import_export:state_getimport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (geomtype == ESMF_GEOMTYPE_MESH) then

       ! get field pointer
       if (present(ungridded_index)) then
          call state_getfldptr(state, trim(fldname), dataPtr2d, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call state_getfldptr(state, trim(fldname), dataPtr1d, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! determine output array
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
                if (present(do_sum) .and. present(ungridded_index)) then
                   output(i,j,iblk)  = output(i,j,iblk) + dataPtr2d(ungridded_index,n)
                else if (present(do_sum)) then
                   output(i,j,iblk)  = output(i,j,iblk) + dataPtr1d(n)
                else if (present(ungridded_index)) then
                   output(i,j,iblk)  = dataPtr2d(ungridded_index,n)
                else
                   output(i,j,iblk) = dataPtr1d(n)
                end if
             end do
          end do
       end do

    else if (geomtype == ESMF_GEOMTYPE_GRID) then

       ! get field pointer
       if (present(ungridded_index)) then
          call state_getfldptr(state, trim(fldname), dataptr4d, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call state_getfldptr(state, trim(fldname), dataptr3d, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! set values of output array
       do iblk = 1,nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi
          do j = jlo,jhi
             do i = ilo,ihi
                i1 = i - ilo + 1
                j1 = j - jlo + 1
                if (present(do_sum) .and. present(ungridded_index)) then
                   output(i,j,iblk) = output(i,j,iblk) + dataPtr4d(i1,j1,iblk,ungridded_index)
                else if (present(do_sum)) then
                   output(i,j,iblk) = output(i,j,iblk) + dataPtr3d(i1,j1,iblk)
                else if (present(ungridded_index)) then
                   output(i,j,iblk) = dataPtr4d(i1,j1,iblk, ungridded_index)
                else
                   output(i,j,iblk) = dataPtr3d(i1,j1,iblk)
                end if
             end do
          end do
       end do

    end if

  end subroutine state_getimport_3d_output

  !===============================================================================

  subroutine state_setexport_4d_input(state, fldname, input, index, lmask, ifrac, ungridded_index, rc)

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
    integer             ,           intent(out)   :: rc

    ! local variables
    type(block)                  :: this_block            ! block information for current block
    integer                      :: ilo, ihi, jlo, jhi    ! beginning and end of physical domain
    integer                      :: i, j, iblk, n, i1, j1 ! indices
    real(kind=dbl_kind), pointer :: dataPtr1d(:)          ! mesh
    real(kind=dbl_kind), pointer :: dataPtr2d(:,:)        ! mesh
    real(kind=dbl_kind), pointer :: dataPtr3d(:,:,:)      ! grid
    real(kind=dbl_kind), pointer :: dataPtr4d(:,:,:,:)    ! grid
    character(len=*), parameter  :: subname='(ice_import_export:state_setexport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (geomtype == ESMF_GEOMTYPE_MESH) then

       if (present(ungridded_index)) then
          write(6,*)'DEBUG: fldname = ',trim(fldname),' has ungridded index= ',ungridded_index
       end if

       ! get field pointer
       if (present(ungridded_index)) then
          call state_getfldptr(state, trim(fldname), dataPtr2d, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call state_getfldptr(state, trim(fldname), dataPtr1d, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! set values of field pointer
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
                         dataPtr2d(ungridded_index,n) = input(i,j,index,iblk)
                      else
                         dataPtr1d(n) = input(i,j,index,iblk)
                      end if
                   end if
                else
                   if (present(ungridded_index)) then
                      dataPtr2d(ungridded_index,n) = input(i,j,index,iblk)
                   else
                      dataPtr1d(n) = input(i,j,index,iblk)
                   end if
                end if
             end do
          end do
       end do

    else if (geomtype == ESMF_GEOMTYPE_GRID) then

       if (present(ungridded_index)) then
          call state_getfldptr(state, trim(fldname), dataptr4d, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call state_getfldptr(state, trim(fldname), dataptr3d, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       do iblk = 1,nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi
          do j = jlo,jhi
             do i = ilo,ihi
                i1 = i - ilo + 1
                j1 = j - jlo + 1
                if (present(lmask) .and. present(ifrac)) then
                   if ( lmask(i,j,iblk) .and. ifrac(i,j,iblk) > c0 ) then
                      if (present(ungridded_index)) then
                         dataPtr4d(i1,j1,iblk,ungridded_index) = input(i,j,index,iblk)
                      end if
                   else
                      dataPtr3d(i1,j1,iblk) = input(i,j,index,iblk)
                   end if
                else
                   if (present(ungridded_index)) then
                      dataPtr4d(i1,j1,iblk,ungridded_index) = input(i,j,index,iblk)
                   else
                      dataPtr3d(i1,j1,iblk) = input(i,j,index,iblk)
                   end if
                end if
             end do
          end do
       end do

    end if

  end subroutine state_setexport_4d_input

  !===============================================================================

  subroutine state_setexport_3d_input(state, fldname, input, lmask, ifrac, ungridded_index, rc)

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
    integer                        , intent(out)   :: rc

    ! local variables
    type(block)                  :: this_block            ! block information for current block
    integer                      :: ilo, ihi, jlo, jhi    ! beginning and end of physical domain
    integer                      :: i, j, iblk, n, i1, j1 ! incides
    real(kind=dbl_kind), pointer :: dataPtr1d(:)          ! mesh
    real(kind=dbl_kind), pointer :: dataPtr2d(:,:)        ! mesh
    real(kind=dbl_kind), pointer :: dataPtr3d(:,:,:)      ! grid
    real(kind=dbl_kind), pointer :: dataPtr4d(:,:,:,:)    ! grid
    character(len=*), parameter  :: subname='(ice_import_export:state_setexport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (geomtype == ESMF_GEOMTYPE_MESH) then

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

    else if (geomtype == ESMF_GEOMTYPE_GRID) then

       ! get field pointer
       if (present(ungridded_index)) then
          call state_getfldptr(state, trim(fldname), dataptr4d, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call state_getfldptr(state, trim(fldname), dataptr3d, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       do iblk = 1,nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi
          do j = jlo,jhi
             do i = ilo,ihi
                i1 = i - ilo + 1
                j1 = j - jlo + 1
                if (present(lmask) .and. present(ifrac)) then
                   if ( lmask(i,j,iblk) .and. ifrac(i,j,iblk) > c0 ) then
                      if (present(ungridded_index)) then
                         dataPtr4d(i1,j1,iblk,ungridded_index) = input(i,j,iblk)
                      else
                         dataPtr3d(i1,j1,iblk) = input(i,j,iblk)
                      end if
                   end if
                else
                   if (present(ungridded_index)) then
                      dataPtr4d(i1,j1,iblk,ungridded_index) = input(i,j,iblk)
                   else
                      dataPtr3d(i1,j1,iblk) = input(i,j,iblk)
                   end if
                end if
             end do
          end do
       end do

    end if

  end subroutine state_setexport_3d_input

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

  !===============================================================================

  subroutine State_GetFldPtr_3d(State, fldname, fldptr, rc)
    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    ,            intent(in)     :: State
    character(len=*)    ,            intent(in)     :: fldname
    real(kind=dbl_kind) , pointer ,  intent(inout)  :: fldptr(:,:,:)
    integer             , optional , intent(out)    :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    character(len=*),parameter :: subname='(ice_import_export:State_GetFldPtr_3d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  end subroutine State_GetFldPtr_3d

  !===============================================================================

  subroutine State_GetFldPtr_4d(State, fldname, fldptr, rc)
    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    ,            intent(in)     :: State
    character(len=*)    ,            intent(in)     :: fldname
    real(kind=dbl_kind) , pointer ,  intent(inout)  :: fldptr(:,:,:,:)
    integer             , optional , intent(out)    :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    character(len=*),parameter :: subname='(ice_import_export:State_GetFldPtr_3d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  end subroutine State_GetFldPtr_4d

end module ice_import_export
