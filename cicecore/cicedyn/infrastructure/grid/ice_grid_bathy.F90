    module ice_grid_bathy

    use ice_kinds_mod
    use ice_communicate, only: my_task, master_task
    use ice_fileunits, only: nu_diag
    use icepack_intfc, only: icepack_warnings_flush, icepack_query_parameters, &
        icepack_init_parameters, icepack_warnings_aborted
    use ice_exit, only: abort_ice
    use ice_domain, only: nblocks
    use ice_blocks, only: nx_block, ny_block

    implicit none
    private
    public :: get_bathymetry, get_bathymetry_popfile

    contains
!=======================================================================
! ocean bathymetry for grounded sea ice (seabed stress) or icebergs
! currently hardwired for 40 levels (gx3, gx1 grids)
! should be read from a file instead (see subroutine read_seabedstress_bathy)

    subroutine get_bathymetry(use_bathymetry, bathymetry_file, bathymetry, kmt)

        use ice_constants, only: c0

        character (len=char_len_long), intent(in) :: bathymetry_file

        real (kind=dbl_kind), dimension (:,:,:), intent(out) :: bathymetry

        logical (kind=log_kind), intent(in) :: use_bathymetry ! flag for reading in bathymetry_file

        real (kind=dbl_kind), dimension (:,:,:), intent(in) :: kmt 

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
  
        character (len=char_len_long), intent(in) :: bathymetry_file

        real (kind=dbl_kind), dimension (:,:,:), intent(out) :: bathymetry

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

!=======================================================================
! with use_bathymetry = false, vertical depth profile generated for max KMT
! with use_bathymetry = true, expects to read in pop vert_grid file

    subroutine get_bathymetry_popfile(use_bathymetry, bathymetry_file, bathymetry, kmt)

        use ice_domain, only: distrb_info
        use ice_global_reductions, only: global_maxval
        use ice_fileunits, only: get_fileunit, release_fileunit
        use ice_broadcast, only: broadcast_array

        character (len=char_len_long), intent(in) :: bathymetry_file

        real (kind=dbl_kind), dimension (:,:,:), intent(out) :: bathymetry

        logical (kind=log_kind), intent(in) :: use_bathymetry ! flag for reading in bathymetry_file

        real (kind=dbl_kind), dimension (:,:,:), intent(in) :: kmt 

        integer (kind=int_kind) :: &
           i, j, k, iblk      ! loop indices
  
        integer (kind=int_kind) :: &
           ntmp, nlevel   , & ! number of levels (max KMT)
           k1             , & ! levels
           ierr           , & ! error tag
           fid                ! fid unit number
  
        real (kind=dbl_kind), dimension(:),allocatable :: &
           depth          , & ! total depth, m
           thick              ! layer thickness, cm -> m
  
        logical (kind=log_kind) :: &
           calc_dragio
  
        character(len=*), parameter :: subname = '(get_bathymetry_popfile)'
  
        ntmp = maxval(nint(KMT))
        nlevel = global_maxval(ntmp,distrb_info)
  
        if (my_task==master_task) then
           write(nu_diag,*) subname,' KMT max = ',nlevel
        endif
  
        allocate(depth(nlevel),thick(nlevel))
        thick = -999999.
        depth = -999999.
  
        if (use_bathymetry) then
  
           write (nu_diag,*) subname,' Bathymetry file = ', trim(bathymetry_file)
           if (my_task == master_task) then
              call get_fileunit(fid)
              open(fid,file=bathymetry_file,form='formatted',iostat=ierr)
              if (ierr/=0) call abort_ice(subname//' open error', file=__FILE__, line=__LINE__)
              do k = 1,nlevel
                 read(fid,*,iostat=ierr) thick(k)
                 if (ierr/=0) call abort_ice(subname//' read error', file=__FILE__, line=__LINE__)
              enddo
              call release_fileunit(fid)
           endif
  
           call broadcast_array(thick,master_task)
  
        else
  
           ! create thickness profile
           k1 = min(5,nlevel)
           do k = 1,k1
              thick(k) = max(10000._dbl_kind/float(nlevel),500._dbl_kind)
           enddo
           do k = k1+1,nlevel
              thick(k) = min(thick(k-1)*1.2_dbl_kind,20000._dbl_kind)
           enddo
  
        endif
  
        ! convert thick from cm to m
        thick = thick / 100._dbl_kind
  
        ! convert to total depth
        depth(1) = thick(1)
        do k = 2, nlevel
           depth(k) = depth(k-1) + thick(k)
           if (depth(k) < 0.) call abort_ice(subname//' negative depth error', file=__FILE__, line=__LINE__)
        enddo
  
        if (my_task==master_task) then
           do k = 1,nlevel
             write(nu_diag,'(2a,i6,2f13.7)') subname,'   k, thick(m), depth(m) = ',k,thick(k),depth(k)
           enddo
        endif
  
        bathymetry = 0._dbl_kind
        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block
              k = nint(kmt(i,j,iblk))
              if (k > nlevel) call abort_ice(subname//' kmt gt nlevel error', file=__FILE__, line=__LINE__)
              if (k > 0) bathymetry(i,j,iblk) = depth(k)
           enddo
           enddo
        enddo
  
        ! For consistency, set thickness_ocn_layer1 in Icepack if 'calc_dragio' is active
        call icepack_query_parameters(calc_dragio_out=calc_dragio)
        if (calc_dragio) then
           call icepack_init_parameters(thickness_ocn_layer1_in=thick(1))
        endif
        call icepack_warnings_flush(nu_diag)
        if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
           file=__FILE__, line=__LINE__)
  
        deallocate(depth,thick)
  
    end subroutine get_bathymetry_popfile  

    end module ice_grid_bathy