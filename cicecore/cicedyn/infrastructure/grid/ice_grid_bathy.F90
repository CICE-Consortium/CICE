    module ice_grid_bathy

    use ice_kinds_mod
    use ice_communicate, only: my_task, master_task
    use ice_fileunits, only: nu_diag
    use icepack_intfc, only: icepack_warnings_flush, icepack_query_parameters, &
        icepack_init_parameters, icepack_warnings_aborted
    use ice_exit, only: abort_ice

    implicit none
    private
    public :: get_bathymetry

    contains
!=======================================================================
! ocean bathymetry for grounded sea ice (seabed stress) or icebergs
! currently hardwired for 40 levels (gx3, gx1 grids)
! should be read from a file instead (see subroutine read_seabedstress_bathy)

    subroutine get_bathymetry(use_bathymetry, bathymetry_file, bathymetry, kmt)

        use ice_constants, only: c0
        use ice_domain, only: nblocks
        use ice_blocks, only: nx_block, ny_block

        ! arguments
        character (len=char_len_long), intent(in) :: bathymetry_file
        real (kind=dbl_kind), dimension (:,:,:), intent(out) :: bathymetry
        logical (kind=log_kind) :: use_bathymetry ! flag for reading in bathymetry_file
        real (kind=dbl_kind), dimension (:,:,:) :: kmt        ! ocean topography mask for bathymetry (T-cell)


        integer (kind=int_kind) :: &
            i, j, k, iblk      ! loop indices

        integer (kind=int_kind), parameter :: &
            nlevel = 40        ! number of layers (gx3 grid)

        real (kind=dbl_kind), dimension(nlevel) :: &
            depth              ! total depth, m

        real (kind=dbl_kind) :: &
            puny

        logical (kind=log_kind) :: &
            calc_dragio

        real (kind=dbl_kind), dimension(nlevel), parameter :: &
            thick  = (/ &                        ! ocean layer thickness, m
                10.01244_dbl_kind,  10.11258_dbl_kind,  10.31682_dbl_kind, &
                10.63330_dbl_kind,  11.07512_dbl_kind,  11.66145_dbl_kind, &
                12.41928_dbl_kind,  13.38612_dbl_kind,  14.61401_dbl_kind, &
                16.17561_dbl_kind,  18.17368_dbl_kind,  20.75558_dbl_kind, &
                24.13680_dbl_kind,  28.63821_dbl_kind,  34.74644_dbl_kind, &
                43.20857_dbl_kind,  55.16812_dbl_kind,  72.30458_dbl_kind, &
                96.74901_dbl_kind,  130.0392_dbl_kind,  170.0489_dbl_kind, &
                207.9933_dbl_kind,  233.5694_dbl_kind,  245.2719_dbl_kind, &
                248.9804_dbl_kind,  249.8322_dbl_kind,  249.9787_dbl_kind, &
                249.9979_dbl_kind,  249.9998_dbl_kind,  250.0000_dbl_kind, &
                250.0000_dbl_kind,  250.0000_dbl_kind,  250.0000_dbl_kind, &
                250.0000_dbl_kind,  250.0000_dbl_kind,  250.0000_dbl_kind, &
                250.0000_dbl_kind,  250.0000_dbl_kind,  250.0000_dbl_kind, &
                250.0000_dbl_kind   /)

        character(len=*), parameter :: subname = '(get_bathymetry)'

        call icepack_query_parameters(puny_out=puny, calc_dragio_out=calc_dragio)
        call icepack_warnings_flush(nu_diag)
        if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
            file=__FILE__, line=__LINE__)

        if (use_bathymetry) then

            call read_seabedstress_bathy(bathymetry_file, bathymetry)

        else

            ! convert to total depth
            depth(1) = thick(1)
            do k = 2, nlevel
                depth(k) = depth(k-1) + thick(k)
            enddo

            bathymetry = c0
            do iblk = 1, nblocks
                do j = 1, ny_block
                do i = 1, nx_block
                    k = min(nint(kmt(i,j,iblk)),nlevel)
                    if (k > nlevel) call abort_ice(subname//' kmt gt nlevel error', &
                                        file=__FILE__, line=__LINE__)
                    if (k > 0) bathymetry(i,j,iblk) = depth(k)
                enddo
                enddo
            enddo

            ! For consistency, set thickness_ocn_layer1 in Icepack if 'calc_dragio' is active
            if (calc_dragio) then
                call icepack_init_parameters(thickness_ocn_layer1_in=thick(1))
                call icepack_warnings_flush(nu_diag)
                if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
                    file=__FILE__, line=__LINE__)
            endif

        endif ! bathymetry_file

    end subroutine get_bathymetry
        
    !=======================================================================

    ! Read bathymetry data for seabed stress calculation (grounding scheme for
    ! landfast ice) in CICE stand-alone mode. When CICE is in coupled mode
    ! (e.g. CICE-NEMO), hwater should be uptated at each time level so that
    ! it varies with ocean dynamics.
    !
    ! author: Fred Dupont, CMC

    subroutine read_seabedstress_bathy(bathymetry_file, bathymetry)

        use ice_read_write, only: ice_open_nc, ice_read_nc, ice_close_nc
        use ice_constants, only: field_loc_center , field_type_scalar
  
        ! arguments
        character (len=char_len_long), intent(in) :: bathymetry_file
        real (kind=dbl_kind), dimension (:,:,:), intent(out) :: bathymetry

        ! local variables
        integer (kind=int_kind) :: &
           fid_init        ! file id for netCDF init file
  
        character (char_len_long) :: &        ! input data file names
           fieldname
  
        logical (kind=log_kind) :: diag=.true.
  
        character(len=*), parameter :: subname = '(read_seabedstress_bathy)'
  
        if (my_task == master_task) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'Bathymetry file: ', trim(bathymetry_file)
            call icepack_warnings_flush(nu_diag)
        endif
  
        call ice_open_nc(bathymetry_file,fid_init)
  
        fieldname='Bathymetry'
  
        if (my_task == master_task) then
           write(nu_diag,*) subname,' reading ',TRIM(fieldname)
           call icepack_warnings_flush(nu_diag)
        endif
        call ice_read_nc(fid_init,1,fieldname,bathymetry,diag, &
                      field_loc=field_loc_center, &
                      field_type=field_type_scalar)
  
        call ice_close_nc(fid_init)
  
        if (my_task == master_task) then
           write(nu_diag,*) subname,' closing file ',TRIM(bathymetry_file)
           call icepack_warnings_flush(nu_diag)
        endif
  
    end subroutine read_seabedstress_bathy

    end module ice_grid_bathy