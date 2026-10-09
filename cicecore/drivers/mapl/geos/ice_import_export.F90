module ice_import_export

  use ESMF
  use ice_kinds_mod      , only : int_kind, dbl_kind, real_kind, char_len, log_kind
  use ice_constants      , only : c0, c1, c2, p5, p25, spval_dbl, radius
  use ice_constants      , only : field_loc_center, field_type_scalar, field_type_vector
  use ice_blocks         , only : block, get_block, nx_block, ny_block, nghost
  use ice_domain         , only : nblocks, blocks_ice, halo_info, distrb_info
  use ice_domain_size    , only : nx_global, ny_global, block_size_x, block_size_y, max_blocks, ncat
  use ice_exit           , only : abort_ice
  use ice_flux           , only : strairxt, strairyt, strocnxT_iavg, strocnyT_iavg
  use ice_flux           , only : alvdr, alidr, alvdf, alidf, Tref, Qref, Uref
  use ice_flux           , only : flat, fsens, flwout, evap, fswabs, fhocn, fswthru
  use ice_flux           , only : evapn_f, fsurfn_f, dfsurfndTsfc_f, dflatndTsfc_f
  use ice_flux           , only : flatn_f, coszen
  use ice_flux           , only : fswthru_uvrdr, fswthru_uvrdf, fswthru_pardr, fswthru_pardf
  use ice_flux           , only : send_i2x_per_cat, fswthrun_ai
  use ice_flux_bgc       , only : faero_atm, faero_ocn
  use ice_flux_bgc       , only : fiso_atm, fiso_ocn, fiso_evap
  use ice_flux_bgc       , only : Qa_iso, Qref_iso, HDO_ocn, H2_18O_ocn, H2_16O_ocn
  use ice_flux           , only : fresh, fsalt, zlvl, uatm, vatm, potT, Tair, Qa
  use ice_flux           , only : fresh_ai, fsalt_ai, fhocn_ai
  use ice_flux           , only : rhoa, swvdr, swvdf, swidr, swidf, flw, frain
  use ice_flux           , only : swuvrdr, swuvrdf, swpardr, swpardf
  use ice_flux           , only : fcondtop
  use ice_flux           , only : fsnow, uocn, vocn, sst, ss_tltx, ss_tlty, frzmlt
  use ice_flux           , only : send_i2x_per_cat
  use ice_flux           , only : sss, Tf, wind, fsw
  use ice_state          , only : vice, vsno, aice, aicen, trcr, trcrn
  use ice_state          , only : Tsfcn_init, aice_init, uvel, vvel
  use ice_grid           , only : tlon, tlat, tarea, tmask, umask, anglet, ocn_gridcell_frac, hm
  use ice_grid           , only : dxu, dyu, dxE, dyE, dxN, dyN, nmask, emask
  use ice_grid           , only : dxT, dyT
  use ice_grid           , only : grid_ice, grid_ocn
  use ice_boundary       , only : ice_HaloUpdate
  use ice_shr_methods    , only : chkerr
  use ice_fileunits      , only : nu_diag, flush_fileunit
  use ice_communicate    , only : my_task, master_task, MPI_COMM_ICE
  use ice_prescribed_mod , only : prescribed_ice
  use icepack_intfc      , only : icepack_warnings_flush, icepack_warnings_aborted
  use icepack_intfc      , only : icepack_query_parameters, icepack_query_tracer_flags
  use icepack_intfc      , only : icepack_liquidus_temperature
  use icepack_intfc      , only : icepack_sea_freezing_temperature

  implicit none
  public

  public  :: ice_import_thermo1
  public  :: ice_export_thermo1
  public  :: ice_import_grid
  !public  :: ice_import_thermo2
  !public  :: ice_export_thermo2
  public  :: ice_import_dyna
  public  :: ice_export_dyna
  public  :: ice_export_field

  interface ice_export_field
     module procedure ice_export_field_2d
     module procedure ice_export_field_3d
  end interface ice_export_field

  private :: state_FldChk

  interface state_getfldptr
     module procedure state_getfldptr_1d
     module procedure state_getfldptr_2d
     module procedure state_getfldptr_3d
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

  interface arr_setexport
     module procedure arr_setexport_4d
     module procedure arr_setexport_3d
  end interface arr_setexport
  private :: arr_setexport

  ! Private module data

  type fld_list_type
    character(len=128) :: stdname
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

  integer     , parameter  :: io_dbug = 10        ! i/o debug messages
  character(*), parameter  :: u_FILE_u = &
       __FILE__

!==============================================================================
contains
!==============================================================================
  subroutine ice_import_grid( fro, rc )

    ! input/output variables
    real(kind=real_kind), dimension(:,:), intent(in)  :: fro
    integer                             , intent(out) :: rc

    integer                          :: i, j, i1, j1, k, iblk, n
    integer                          :: ilo, ihi, jlo, jhi !beginning and end of physical domain
    type(block)                      :: this_block         ! block information for current block

    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
          j1 = j - nghost
          do i = ilo, ihi
             i1 = i - nghost
             ocn_gridcell_frac(i,j,iblk)  = real(fro(i1, j1), kind=dbl_kind)
          enddo
       enddo
    enddo

    rc = ESMF_SUCCESS

  end subroutine ice_import_grid
  !===============================================================================

  subroutine ice_import_thermo1( importState, rc )

    ! input/output variables
    type(ESMF_State) , intent(in)  :: importState
    integer          , intent(out) :: rc

    ! local variables
    integer,parameter                :: nfldu=5
    integer,parameter                :: nfld=12
    integer                          :: i, j, k, iblk, n
    integer                          :: ilo, ihi, jlo, jhi !beginning and end of physical domain
    type(block)                      :: this_block         ! block information for current block
    real (kind=dbl_kind)             :: Tffresh
    real (kind=dbl_kind),allocatable :: afldu(:,:,:,:,:)
    real (kind=dbl_kind),allocatable :: afld(:,:,:,:)
    character(len=*),   parameter    :: subname = 'ice_import_thermo1'
    character(len=1024)              :: msgString
    !-----------------------------------------------------

    call icepack_query_parameters(Tffresh_out=Tffresh)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
        file=u_FILE_u, line=__LINE__)


    allocate(afldu(nx_block,ny_block,ncat,nfldu,nblocks))
    afldu = c0
    allocate( afld(nx_block,ny_block,      nfld,nblocks))
    afld = c0

    !call state_getimport(importState,  'TSKINICE', output=afldu, index=1, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState,      'EVAP', output=afldu, index=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState,     'FSURF', output=afldu, index=2, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'DFSURFDTS', output=afldu, index=3, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState,   'DLHFDTS', output=afldu, index=4, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState,       'LHF', output=afldu, index=5, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return


    call state_getimport(importState,      'SNOW', output=afld,  index=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState,      'RAIN', output=afld,  index=2, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState,     'DRPAR', output=afld,  index=3, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState,     'DFPAR', output=afld,  index=4, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState,     'DRNIR', output=afld,  index=5, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState,     'DFNIR', output=afld,  index=6, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState,     'DRUVR', output=afld,  index=7, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState,     'DFUVR', output=afld,  index=8, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState,      'COSZ', output=afld,  index=9, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState,      'SST', output=afld,  index=10, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState,      'SSS', output=afld,  index=11, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState,   'FRZMLT', output=afld,  index=12, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! now fill in the ice internal data types
    do k = 1, ncat
       !$OMP PARALLEL DO PRIVATE(iblk,i,j)
       do iblk = 1, nblocks
          do j = 1,ny_block
             do i = 1,nx_block
                !trcrn  (i,j,1,k,iblk)      = afldu(i,j,k,1,iblk) - Tffresh
                evapn_f  (i,j,k,iblk)      = afldu(i,j,k,1,iblk)
                fsurfn_f (i,j,k,iblk)      = afldu(i,j,k,2,iblk)
                dfsurfndTsfc_f(i,j,k,iblk) = afldu(i,j,k,3,iblk)
                dflatndTsfc_f(i,j,k,iblk)  = afldu(i,j,k,4,iblk)
                flatn_f  (i,j,k,iblk)      = afldu(i,j,k,5,iblk)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end do

    !$OMP PARALLEL DO PRIVATE(iblk,i,j)
    do iblk = 1, nblocks
       do j = 1,ny_block
          do i = 1,nx_block
             fsnow (i,j,iblk)         = afld(i,j,1,iblk)
             frain (i,j,iblk)         = afld(i,j,2,iblk)
             swvdr (i,j,iblk)         = afld(i,j,3,iblk) + afld(i,j,7,iblk)
             swidr (i,j,iblk)         = afld(i,j,5,iblk)
             swvdf (i,j,iblk)         = afld(i,j,4,iblk) + afld(i,j,8,iblk)
             swidf (i,j,iblk)         = afld(i,j,6,iblk)
             swuvrdr(i,j,iblk)        = afld(i,j,7,iblk)
             swuvrdf(i,j,iblk)        = afld(i,j,8,iblk)
             swpardr(i,j,iblk)        = afld(i,j,3,iblk)
             swpardf(i,j,iblk)        = afld(i,j,4,iblk)
             coszen(i,j,iblk)         = afld(i,j,9,iblk)
             sst   (i,j,iblk)         = afld(i,j,10,iblk) - Tffresh
             sss   (i,j,iblk)         = afld(i,j,11,iblk)
             frzmlt(i,j,iblk)         = afld(i,j,12,iblk)
             fsw   (i,j,iblk)         = swvdr(i,j,iblk) + swvdf(i,j,iblk) &
                                      + swidr(i,j,iblk) + swidf(i,j,iblk)
          end do
       end do
    end do
    !$OMP END PARALLEL DO



    !== will change to read in from coupler once Tf from MOM is ready
    !$OMP PARALLEL DO PRIVATE(iblk,i,j)
    do iblk = 1, nblocks
       do j = 1,ny_block
          do i = 1,nx_block
            Tf(i,j,iblk) = icepack_sea_freezing_temperature(sss(i,j,iblk))
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    deallocate(afldu)
    deallocate(afld)

    rc = ESMF_SUCCESS

  end subroutine ice_import_thermo1

  subroutine ice_import_radiation( importState, rc )

    ! input/output variables
    type(ESMF_State) , intent(in)  :: importState
    integer          , intent(out) :: rc

    ! local variables
    integer,parameter                :: nfld=1
    integer                          :: i, j, k, iblk
    integer                          :: ilo, ihi, jlo, jhi !beginning and end of physical domain
    type(block)                      :: this_block         ! block information for current block
    real (kind=dbl_kind),allocatable :: afld(:,:,:,:)
    character(len=*),   parameter    :: subname = 'ice_import_radiation'
    character(len=1024)              :: msgString
    !-----------------------------------------------------


    allocate( afld(nx_block,ny_block,      nfld,nblocks))
    afld = c0

    call state_getimport(importState,      'COSZ', output=afld,  index=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !$OMP PARALLEL DO PRIVATE(iblk,i,j)
    do iblk = 1, nblocks
       do j = 1,ny_block
          do i = 1,nx_block
             coszen(i,j,iblk)         = afld(i,j,1,iblk)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    deallocate(afld)

    rc = ESMF_SUCCESS

  end subroutine ice_import_radiation

  subroutine ice_import_dyna( taux, tauy, slv, uoa, voa, uob, vob, uoc, voc, rc )

    ! input/output variables
    real(kind=real_kind) ,     intent(in) :: taux(:,:)
    real(kind=real_kind) ,     intent(in) :: tauy(:,:)
    real(kind=real_kind) ,     intent(in) :: slv(:,:)
    real(kind=real_kind) ,     intent(in) :: uoa(:,:)
    real(kind=real_kind) ,     intent(in) :: voa(:,:)
    real(kind=real_kind) ,     intent(in) :: uob(:,:)
    real(kind=real_kind) ,     intent(in) :: vob(:,:)
    real(kind=real_kind) ,     intent(in) :: uoc(:,:)
    real(kind=real_kind) ,     intent(in) :: voc(:,:)
    integer              ,     intent(out):: rc

    ! local variables
    integer                          :: i, j, k, iblk
    integer                          :: i1, j1
    integer                          :: ilo, ihi, jlo, jhi !beginning and end of physical domain
    type(block)                      :: this_block         ! block information for current block
    real(kind=dbl_kind)              :: workx, worky
    real(kind=dbl_kind)              :: slp_L, slp_R, slp_C
    real(kind=dbl_kind)              :: u_min, u_max, slope
    real(kind=dbl_kind)              :: ssh(nx_block,ny_block,max_blocks)
    character(len=*),   parameter    :: subname = 'ice_import_dyna'
    character(len=1024)              :: msgString
    !-----------------------------------------------------


    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
          j1 = j - nghost
          do i = ilo, ihi
             i1 = i - nghost
                if(tmask(i,j,iblk)) then
                   workx  = real(taux(i1, j1), kind=dbl_kind)
                   worky  = real(tauy(i1, j1), kind=dbl_kind)
                   strairxT(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid
                                      + worky*sin(ANGLET(i,j,iblk))   ! note strax, stray, wind
                   strairyT(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) & !  are on the T-grid here
                                      - workx*sin(ANGLET(i,j,iblk))
                   ! multiply by aice to properly treat free drift
                   strairxT(i,j,iblk) = strairxT(i,j,iblk) * aice_init(i,j,iblk)
                   strairyT(i,j,iblk) = strairyT(i,j,iblk) * aice_init(i,j,iblk)
                   ssh(i,j,iblk)      = real(slv(i1,j1), kind=dbl_kind)
                else
                   strairxT(i,j,iblk) = c0
                   strairyT(i,j,iblk) = c0
                   ssh(i,j,iblk)      = c0
                endif
                if(trim(grid_ocn) == 'A') then
                   if(tmask(i,j,iblk)) then
                     workx = real(uoa(i1,j1), kind=dbl_kind)
                     worky = real(voa(i1,j1), kind=dbl_kind)
                     uocn(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid
                                    + worky*sin(ANGLET(i,j,iblk))   ! note strax, stray, wind
                     vocn(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) & !  are on the T-grid here
                                    - workx*sin(ANGLET(i,j,iblk))
                   else
                     uocn(i,j,iblk) = c0
                     vocn(i,j,iblk) = c0
                   endif
                elseif(trim(grid_ice) == 'B') then
                   if(umask(i,j,iblk)) then
                     uocn(i,j,iblk) = real(uob(i1,j1), kind=dbl_kind)
                     vocn(i,j,iblk) = real(vob(i1,j1), kind=dbl_kind)
                   else
                     uocn(i,j,iblk) = c0
                     vocn(i,j,iblk) = c0
                   endif
                elseif(trim(grid_ice) == 'C') then
                   if(emask(i,j,iblk)) then
                     uocn(i,j,iblk) = real(uoc(i1,j1), kind=dbl_kind)
                   else
                     uocn(i,j,iblk) = c0
                   endif
                   if(nmask(i,j,iblk)) then
                     vocn(i,j,iblk) = real(voc(i1,j1), kind=dbl_kind)
                   else
                     vocn(i,j,iblk) = c0
                   endif
                else
                   call abort_ice(error_message='unknown grid_ice', &
                       file=u_FILE_u, line=__LINE__)
                endif
          enddo
       enddo
    enddo

    call ice_HaloUpdate (ssh,     halo_info, &
                         field_loc_center, field_type_scalar)

    !*** if  C-grid ice dynamics is on, the following needs to be revised
    !*** so a query of which dynamics (B- or C-) should be made and branch into
    !*** different computation of ss_tlt* terms
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
       do i = ilo, ihi
           if(trim(grid_ocn) == 'A') then
             if(tmask(i,j,iblk)) then
                slp_L = ssh(I,j,iblk) - ssh(I-1,j,iblk)
                if(.not. emask(i-1,j,iblk)) slp_L = c0
                slp_R = ssh(I+1,j,iblk) - ssh(I,j,iblk)
                if(.not. emask(i,j,iblk)) slp_R = c0
                slp_C = p5 * (slp_L + slp_R)
                if ( (slp_L * slp_R) > c0 ) then
                  ! This limits the slope so that the edge values are bounded by the
                  ! two cell averages spanning the edge.
                  u_min = min( ssh(i-1,j,iblk), ssh(i,j,iblk), ssh(i+1,j,iblk) )
                  u_max = max( ssh(i-1,j,iblk), ssh(i,j,iblk), ssh(i+1,j,iblk) )
                  slope = sign( min( abs(slp_C), c2*min( ssh(i,j,iblk) - u_min, u_max - ssh(i,j,iblk) ) ), slp_C )
                else
                  ! Extrema in the mean values require a PCM reconstruction avoid generating
                  ! larger extreme values.
                  slope = c0
                endif
                ss_tltx(i,j,iblk) = slope / dxT(i,j,iblk)

                slp_L = ssh(I,j,iblk) - ssh(I,j-1,iblk)
                if(.not. nmask(i,j-1,iblk)) slp_L = c0
                slp_R = ssh(I,j+1,iblk) - ssh(I,j,iblk)
                if(.not. nmask(i,j,iblk)) slp_R = c0
                slp_C = p5 * (slp_L + slp_R)
                if ( (slp_L * slp_R) > c0 ) then
                  ! This limits the slope so that the edge values are bounded by the
                  ! two cell averages spanning the edge.
                  u_min = min( ssh(i,j-1,iblk), ssh(i,j,iblk), ssh(i,j+1,iblk) )
                  u_max = max( ssh(i,j-1,iblk), ssh(i,j,iblk), ssh(i,j+1,iblk) )
                  slope = sign( min( abs(slp_C), c2*min( ssh(i,j,iblk) - u_min, u_max - ssh(i,j,iblk) ) ), slp_C )
                else
                  ! Extrema in the mean values require a PCM reconstruction avoid generating
                  ! larger extreme values.
                  slope = c0
                endif
                ss_tlty(i,j,iblk) = slope / dyT(i,j,iblk)
             else
                ss_tltx(i,j,iblk) = c0
                ss_tlty(i,j,iblk) = c0
             endif
           elseif(trim(grid_ice) == 'B') then
             if(umask(i,j,iblk)) then
                ss_tltx(i,j,iblk) = p5*(ssh(i+1,j+1,iblk)-ssh(i,j+1,iblk)  &
                                       +ssh(i+1,j  ,iblk)-ssh(i,j  ,iblk)) &
                                       /dxu(i,j,iblk)
                ss_tlty(i,j,iblk) = p5*(ssh(i+1,j+1,iblk)+ssh(i,j+1,iblk)  &
                                       -ssh(i+1,j  ,iblk)-ssh(i,j  ,iblk)) &
                                       /dyu(i,j,iblk)
             else
                ss_tltx(i,j,iblk) = c0
                ss_tlty(i,j,iblk) = c0
             endif
           elseif(trim(grid_ice) == 'C') then
             if(emask(i,j,iblk)) then
                ss_tltx(i,j,iblk) = (ssh(i+1,j  ,iblk)-ssh(i,j  ,iblk)) &
                                     /dxE(i,j,iblk)
             else
                ss_tltx(i,j,iblk) = c0
             endif
             if(nmask(i,j,iblk)) then
                ss_tlty(i,j,iblk) = (ssh(i,j+1,iblk) - ssh(i,j  ,iblk)) &
                                     /dyN(i,j,iblk)
             else
                ss_tlty(i,j,iblk) = c0
             endif
           else
             call abort_ice(error_message='unknown grid_ice', &
                     file=u_FILE_u, line=__LINE__)
           endif
       enddo
       enddo
    enddo

    rc = ESMF_SUCCESS

  end subroutine ice_import_dyna

  subroutine ice_export_dyna( tauxo, tauyo, ui, vi, rc )

    ! input/output variables
    real(kind=real_kind) ,     intent(out) :: tauxo(:,:)
    real(kind=real_kind) ,     intent(out) :: tauyo(:,:)
    real(kind=real_kind) ,     intent(out) :: ui   (:,:)
    real(kind=real_kind) ,     intent(out) :: vi   (:,:)
    integer              ,     intent(out):: rc

    ! local variables
    integer                          :: i, j, k, iblk
    integer                          :: i1, j1
    integer                          :: ilo, ihi, jlo, jhi !beginning and end of physical domain
    type(block)                      :: this_block         ! block information for current block
    real(kind=dbl_kind)              :: workx, worky
    character(len=*),   parameter    :: subname = 'ice_export_dyna'
    character(len=1024)              :: msgString
    !-----------------------------------------------------


    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
          j1 = j - nghost
          do i = ilo, ihi
             i1 = i - nghost
             ! ice/ocean stress (on POP T-grid:  convert to lat-lon)
             workx        = strocnxT_iavg(i,j,iblk)                            ! N/m^2
             worky        = strocnyT_iavg(i,j,iblk)                            ! N/m^2
             tauxo(i1,j1) = real(workx*cos(ANGLET(i,j,iblk)) - &
                                 worky*sin(ANGLET(i,j,iblk)), kind=real_kind)
             tauyo(i1,j1) = real(worky*cos(ANGLET(i,j,iblk)) + &
                                 workx*sin(ANGLET(i,j,iblk)), kind=real_kind)
             workx        = p25*(uvel(i,j  ,iblk) + uvel(i-1,j  ,iblk) & ! cell-centered velocity
                               + uvel(i,j-1,iblk) + uvel(i-1,j-1,iblk))  ! assumes wind components
             worky        = p25*(vvel(i,j  ,iblk) + vvel(i-1,j  ,iblk) & ! are also cell-centered
                               + vvel(i,j-1,iblk) + vvel(i-1,j-1,iblk))
             ui(i1,j1)    = real(workx*cos(ANGLET(i,j,iblk)) - &
                                 worky*sin(ANGLET(i,j,iblk)), kind=real_kind)
             vi(i1,j1)    = real(worky*cos(ANGLET(i,j,iblk)) + &
                                 workx*sin(ANGLET(i,j,iblk)), kind=real_kind)
          enddo
       enddo
    enddo

    rc = ESMF_SUCCESS

  end subroutine ice_export_dyna

  !===============================================================================
  subroutine ice_export_thermo1( exportState, rc )

    ! input/output variables
    type(ESMF_State), intent(inout) :: exportState
    integer         , intent(out)   :: rc

    ! local variables
    type(block)             :: this_block                           ! block information for current block
    integer                 :: i, j, iblk, n                        ! incides
    integer                 :: n2                                   ! thickness category index
    integer                 :: ilo, ihi, jlo, jhi                   ! beginning and end of physical domain
    real    (kind=dbl_kind) :: workx, worky                         ! tmps for converting grid
    logical                 :: flag
    integer (kind=int_kind) :: indxi (nx_block*ny_block)            ! compressed indices in i
    integer (kind=int_kind) :: indxj (nx_block*ny_block)            ! compressed indices in i
    real    (kind=dbl_kind) :: dTsrf  (nx_block,ny_block,ncat,max_blocks) ! surface temperature
    real    (kind=dbl_kind) :: Tffresh
    !integer,               parameter :: nfldu=1
    !real (kind=dbl_kind),allocatable :: afldu(:,:,:,:,:)
    character(len=*),parameter :: subname = 'ice_export_thermo1'
    !-----------------------------------------------------

    rc = ESMF_SUCCESS
    if (io_dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call icepack_query_parameters(Tffresh_out=Tffresh)
    !    call icepack_query_parameters(tfrz_option_out=tfrz_option, &
    !       modal_aero_out=modal_aero, z_tracers_out=z_tracers, skl_bgc_out=skl_bgc, &
    !       Tffresh_out=Tffresh)
    !    call icepack_query_tracer_flags(tr_aero_out=tr_aero, tr_iage_out=tr_iage, &
    !       tr_FY_out=tr_FY, tr_pond_out=tr_pond, tr_lvl_out=tr_lvl, &
    !       tr_zaero_out=tr_zaero, tr_bgc_Nit_out=tr_bgc_Nit)

    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
        file=u_FILE_u, line=__LINE__)

    !---------------------------------
    ! Create the export state
    !---------------------------------
    !allocate(afldu(nx_block,ny_block,ncat,nfldu,nblocks))
    !afldu = c0
    !call state_getimport(exportState,  'TSKINICE', output=afldu, index=1, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return

    dTsrf = trcrn(:,:,1,:,:) - Tsfcn_init

    call state_setexport(exportState, 'DTS', input=dTsrf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'ALBVR', input=alvdr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport(exportState, 'ALBVF', input=alvdf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport(exportState, 'ALBNR', input=alidr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport(exportState, 'ALBNF', input=alidf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'PENUVR', input=fswthru_uvrdr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport(exportState, 'PENUVF', input=fswthru_uvrdf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport(exportState, 'PENPAR', input=fswthru_pardr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport(exportState, 'PENPAF', input=fswthru_pardf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'GHTSKIN', input=fcondtop, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !deallocate(afldu)

  end subroutine ice_export_thermo1

  subroutine ice_export_radiation( exportState, rc )

    ! input/output variables
    type(ESMF_State), intent(inout) :: exportState
    integer         , intent(out)   :: rc

    ! local variables
    character(len=*),parameter :: subname = 'ice_export_radiation'
    !-----------------------------------------------------

    rc = ESMF_SUCCESS
    if (io_dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call state_setexport(exportState, 'ALBVR', input=alvdr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport(exportState, 'ALBVF', input=alvdf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport(exportState, 'ALBNR', input=alidr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport(exportState, 'ALBNF', input=alidf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine ice_export_radiation

  subroutine ice_export_field_3d(fldname, fld, rc)

    ! input/output variables
    character(len=*)     ,     intent(in)  :: fldname
    real(kind=real_kind) ,     intent(out) :: fld(:,:,:)
    integer              ,     intent(out) :: rc

    ! local variables
    character(len=*),parameter :: subname = 'ice_export_field_3d'
    real(kind=dbl_kind)        :: Tffresh

    rc = ESMF_SUCCESS

    call icepack_query_parameters(Tffresh_out=Tffresh)
    !    call icepack_query_parameters(tfrz_option_out=tfrz_option, &
    !       modal_aero_out=modal_aero, z_tracers_out=z_tracers, skl_bgc_out=skl_bgc, &
    !       Tffresh_out=Tffresh)
    !    call icepack_query_tracer_flags(tr_aero_out=tr_aero, tr_iage_out=tr_iage, &
    !       tr_FY_out=tr_FY, tr_pond_out=tr_pond, tr_lvl_out=tr_lvl, &
    !       tr_zaero_out=tr_zaero, tr_bgc_Nit_out=tr_bgc_Nit)

    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
        file=u_FILE_u, line=__LINE__)

    if (trim(fldname) == 'TI') then
        call arr_setexport_4d(fld,  trcrn(:,:,1,:,:), rc)
        fld(:,:,:)  = real(Tffresh, kind=real_kind) + fld(:,:,:)    !Kelvin (original ???)
    elseif (trim(fldname) == 'FRSEAICE') then
        call arr_setexport_4d(fld,  aicen,            rc)
    else
        call ESMF_LogWrite(trim(subname)//": "//trim(fldname)//" not available for export", ESMF_LOGMSG_ERROR)
        rc = ESMF_FAILURE
        return
    endif

  end subroutine ice_export_field_3d

  subroutine ice_export_field_2d(fldname, fld, rc)

    ! input/output variables
    character(len=*)     ,     intent(in)  :: fldname
    real(kind=real_kind) ,     intent(out) :: fld(:,:)
    integer              ,     intent(out) :: rc

    ! local variables
    character(len=*),parameter :: subname = 'ice_export_field_2d'


    rc = ESMF_SUCCESS

    if (trim(fldname) == 'FRACICE') then
        call arr_setexport_3d(fld,  aice,        rc)
    elseif (trim(fldname) == 'FHOCN') then
        call arr_setexport_3d(fld,  fhocn_ai,    rc)
    elseif (trim(fldname) == 'FRESH') then
        call arr_setexport_3d(fld,  fresh_ai,    rc)
    elseif (trim(fldname) == 'FSALT') then
        call arr_setexport_3d(fld,  fsalt_ai,    rc)
    else
        call ESMF_LogWrite(trim(subname)//": "//trim(fldname)//" not available for export", ESMF_LOGMSG_ERROR)
        rc = ESMF_FAILURE
        return
    endif

  end subroutine ice_export_field_2d

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
  subroutine state_getimport_4d(state, fldname, output, index, rc)

    ! ----------------------------------------------
    ! Map import state field to output array
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)              , intent(in)    :: state
    character(len=*)              , intent(in)    :: fldname
    real (kind=dbl_kind)          , intent(inout) :: output(:,:,:,:,:)
    integer                       , intent(in)    :: index
    integer                       , intent(out)   :: rc

    ! local variables
    type(block)                   :: this_block            ! block information for current block
    integer                       :: ilo, ihi, jlo, jhi    ! beginning and end of physical domain
    integer                       :: i, j, k, iblk, n, i1, j1 ! incides
    !real(kind=real_kind), pointer :: dataPtr1d(:,:)          ! mesh
    real(kind=real_kind), pointer :: dataPtr3d(:,:,:)        ! mesh

    character(len=*), parameter  :: subname='(ice_import_export:state_getimport_4d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! check that fieldname exists
    if (.not. State_FldChk(state, trim(fldname))) return

    ! get field pointer
    call state_getfldptr(state, trim(fldname), dataPtr3d, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set values of output array
    do k = 1, ncat
       do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi
          do j = jlo, jhi
             j1 = j - nghost
             do i = ilo, ihi
                i1 = i - nghost
                output(i,j,k,index,iblk)  = real(dataPtr3d(i1,j1,k), kind=dbl_kind)
             end do
          end do
       end do
    end do

  end subroutine state_getimport_4d

  !===============================================================================
  subroutine state_getimport_3d(state, fldname, output, index, rc)

    ! ----------------------------------------------
    ! Map import state field to output array
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)              , intent(in)    :: state
    character(len=*)              , intent(in)    :: fldname
    real (kind=dbl_kind)          , intent(inout) :: output(:,:,:,:)
    integer                       , intent(in)    :: index
    integer                       , intent(out)   :: rc

    ! local variables
    type(block)                  :: this_block            ! block information for current block
    integer                      :: ilo, ihi, jlo, jhi    ! beginning and end of physical domain
    integer                      :: i, j, iblk, n, i1, j1 ! incides
    real(kind=real_kind),pointer :: dataPtr2d(:,:)        ! mesh
    character(len=*) , parameter :: subname='(ice_import_export:state_getimport_3d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! check that fieldname exists
    if (.not. State_FldChk(state, trim(fldname))) return

    ! get field pointer
    call state_getfldptr(state, trim(fldname), dataPtr2d, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine output array
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
          j1 = j - nghost
          do i = ilo, ihi
             i1 = i - nghost
             output(i,j,index,iblk)  = real(dataPtr2d(i1, j1), kind=dbl_kind)
          end do
       end do
    end do

  end subroutine state_getimport_3d

  !===============================================================================
  subroutine arr_setexport_4d(output, input, rc)

    ! ----------------------------------------------
    ! Map 4d input array to export state field
    ! ----------------------------------------------

    ! input/output variables
    real(kind=real_kind) ,         intent(out)    :: output(:,:,:)
    real(kind=dbl_kind) ,           intent(in)    :: input(:,:,:,:)
    integer             ,           intent(out)   :: rc

    ! local variables
    type(block)                  :: this_block            ! block information for current block
    integer                      :: ilo, ihi, jlo, jhi    ! beginning and end of physical domain
    integer                      :: i, j, k, iblk, n, i1, j1 ! indices
    character(len=*), parameter  :: subname='(ice_import_export:arr_setexport_4d)'
    ! ----------------------------------------------


    do k = 1, ncat
       do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo; ihi = this_block%ihi
          jlo = this_block%jlo; jhi = this_block%jhi
          do j = jlo, jhi
            j1 = j - nghost
            do i = ilo, ihi
                i1 = i - nghost
                output(i1,j1,k) = real(input(i,j,k,iblk), kind=real_kind)
            end do
          end do
      end do
    end do

    rc = ESMF_SUCCESS

  end subroutine arr_setexport_4d

  !===============================================================================
  subroutine arr_setexport_3d(output, input, rc)

    ! ----------------------------------------------
    ! Map 3d input array to export array
    ! ----------------------------------------------

    ! input/output variables
    real(kind=real_kind) ,         intent(out)    :: output(:,:)
    real(kind=dbl_kind) ,           intent(in)    :: input(:,:,:)
    integer             ,           intent(out)   :: rc

    ! local variables
    type(block)                  :: this_block            ! block information for current block
    integer                      :: ilo, ihi, jlo, jhi    ! beginning and end of physical domain
    integer                      :: i, j, k, iblk, n, i1, j1 ! indices
    character(len=*), parameter  :: subname='(ice_import_export:arr_setexport_3d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do iblk = 1, nblocks
        this_block = get_block(blocks_ice(iblk),iblk)
        ilo = this_block%ilo; ihi = this_block%ihi
        jlo = this_block%jlo; jhi = this_block%jhi
        do j = jlo, jhi
          j1 = j - nghost
          do i = ilo, ihi
             i1 = i - nghost
             output(i1,j1) = real(input(i,j,iblk), kind=real_kind)
          end do
        end do
    end do

  end subroutine arr_setexport_3d

  !===============================================================================
  subroutine state_setexport_4d(state, fldname, input, rc)

    ! ----------------------------------------------
    ! Map 4d input array to export state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    ,           intent(inout) :: state
    character(len=*)    ,           intent(in)    :: fldname
    real(kind=dbl_kind) ,           intent(in)    :: input(:,:,:,:)
    integer             ,           intent(out)   :: rc

    ! local variables
    type(block)                  :: this_block            ! block information for current block
    integer                      :: ilo, ihi, jlo, jhi    ! beginning and end of physical domain
    integer                      :: i, j, k, iblk, n, i1, j1 ! indices
    real(kind=real_kind), pointer :: dataPtr3d(:,:,:)        ! mesh
    character(len=*), parameter  :: subname='(ice_import_export:state_setexport_4d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! check that fieldname exists
    if (.not. State_FldChk(state, trim(fldname))) return

    ! get field pointer
    call state_getfldptr(state, trim(fldname), dataPtr3d, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do k = 1, ncat
       do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo; ihi = this_block%ihi
          jlo = this_block%jlo; jhi = this_block%jhi
          do j = jlo, jhi
            j1 = j - nghost
            do i = ilo, ihi
                i1 = i - nghost
                dataPtr3d(i1,j1,k) = real(input(i,j,k,iblk), kind=real_kind)
            end do
          end do
      end do
    end do

  end subroutine state_setexport_4d

  !===============================================================================
  subroutine state_setexport_3d(state, fldname, input, rc)

    ! ----------------------------------------------
    ! Map 3d input array to export state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)               , intent(inout) :: state
    character(len=*)               , intent(in)    :: fldname
    real(kind=dbl_kind)            , intent(in)    :: input(:,:,:)
    integer                        , intent(out)   :: rc

    ! local variables
    type(block)                  :: this_block            ! block information for current block
    integer                      :: ilo, ihi, jlo, jhi    ! beginning and end of physical domain
    integer                      :: i, j, iblk, n, i1, j1 ! incides
    real(kind=real_kind), pointer :: dataPtr2d(:,:)        ! mesh
    integer                      :: num_ice
    character(len=*), parameter  :: subname='(ice_import_export:state_setexport_3d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! check that fieldname exists
    if (.not. State_FldChk(state, trim(fldname))) return

    ! get field pointer
    call state_getfldptr(state, trim(fldname), dataPtr2d, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    n = 0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
          j1 = j - nghost
          do i = ilo, ihi
             i1 = i - nghost
             dataPtr2d(i1, j1) = input(i,j,iblk)
          end do
       end do
    end do

  end subroutine state_setexport_3d

  !===============================================================================
  subroutine State_GetFldPtr_1d(State, fldname, fldptr, rc)
    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)              , intent(in)     :: State
    character(len=*)              , intent(in)     :: fldname
    real(kind=real_kind), pointer , intent(inout)  :: fldptr(:)
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
  subroutine state_getfldptr_2d(state, fldname, fldptr, rc)
    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    ,            intent(in)     :: state
    character(len=*)    ,            intent(in)     :: fldname
    real(kind=real_kind),  pointer , intent(inout)  :: fldptr(:,:)
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

  end subroutine state_getfldptr_2d

  subroutine State_GetFldPtr_3d(State, fldname, fldptr, rc)
    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    ,            intent(in)     :: State
    character(len=*)    ,            intent(in)     :: fldname
    real(kind=real_kind),  pointer , intent(inout)  :: fldptr(:,:,:)
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

end module ice_import_export
