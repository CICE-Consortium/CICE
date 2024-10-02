    module ice_grid_bathy

    use ice_kinds_mod
    use ice_communicate, only: my_task, master_task
    use ice_fileunits, only: nu_diag
    use icepack_intfc, only: icepack_warnings_flush

    implicit none
    private
    public :: read_seabedstress_bathy

    contains

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