!=======================================================================
!
! Reads and interpolates forcing data for biogeochemistry
!
! authors:  Nicole Jeffery, LANL
!          Elizabeth C. Hunke, LANL
!
      module ice_forcing_bgc

      use ice_kinds_mod
      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size, only: max_blocks
      use ice_communicate, only: my_task, master_task
      use ice_calendar, only: dt, istep, sec, mday, month
      use ice_fileunits, only: nu_diag
      use ice_arrays_column, only: restore_bgc, &
         bgc_data_dir, fe_data_type
      use ice_constants, only: c0, p1
      use ice_constants, only: field_loc_center, field_type_scalar
      use ice_exit, only: abort_ice
      use ice_forcing, only: bgc_data_type
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_nspint, icepack_max_aero, &
          icepack_max_algae, icepack_max_doc, icepack_max_dic
      use icepack_intfc, only: icepack_query_tracer_flags, &
          icepack_query_parameters, icepack_query_parameters, &
          icepack_query_tracer_indices

      implicit none
      private
      public :: get_forcing_bgc, get_atm_bgc, fzaero_data, alloc_forcing_bgc, &
                init_bgc_data, faero_data, faero_default, faero_optics

      integer (kind=int_kind) :: &
         bgcrecnum = 0   ! old record number (save between steps)

      real (kind=dbl_kind), dimension(:,:,:), allocatable :: &
          nitdat      , & ! data value toward which nitrate is restored
          sildat          ! data value toward which silicate is restored

      real (kind=dbl_kind), dimension(:,:,:,:), allocatable, save :: &
         nit_data, & ! field values at 2 temporal data points
         sil_data

!=======================================================================

      contains

!=======================================================================
!
! Allocate space for forcing_bgc variables
!
      subroutine alloc_forcing_bgc

      integer (int_kind) :: ierr

      allocate( &
          nitdat  (nx_block,ny_block,max_blocks), & ! data value toward which nitrate is restored
          sildat  (nx_block,ny_block,max_blocks), & ! data value toward which silicate is restored
          nit_data(nx_block,ny_block,2,max_blocks), & ! field values at 2 temporal data points
          sil_data(nx_block,ny_block,2,max_blocks), &
          stat=ierr)
      if (ierr/=0) call abort_ice('(alloc_forcing_bgc): Out of memory')

      end subroutine alloc_forcing_bgc

!=======================================================================
!
! Read and interpolate annual climatologies of silicate and nitrate.
! Restore model quantities to data if desired.
!
! author: Elizabeth C. Hunke, LANL

      subroutine get_forcing_bgc

      use ice_blocks, only: block, get_block
      use ice_domain, only: nblocks, blocks_ice
      use ice_arrays_column, only: ocean_bio_all
      use ice_calendar, only:  yday
      use ice_flux, only: sss
      use ice_flux_bgc, only: sil, nit
      use ice_forcing, only: trestore, trest, fyear, &
          read_clim_data_nc, interpolate_data, &
          interp_coeff_monthly, interp_coeff,  &
          read_data_nc_point, c1intp, c2intp

      integer (kind=int_kind) :: &
         i, j, iblk     , & ! horizontal indices
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         ixm,ixp, ixx   , & ! record numbers for neighboring months
         maxrec         , & ! maximum record number
         recslot        , & ! spline slot for current record
         midmonth       , & ! middle day of month
         recnum         , & ! record number
         dataloc        , & ! = 1 for data located in middle of time interval
                            ! = 2 for date located at end of time interval
         ks                 ! bgc tracer index (bio_index_o)

      character (char_len_long) :: &
         met_file,   &    ! netcdf filename
         fieldname        ! field name in netcdf file

      real (kind=dbl_kind), dimension(2), save :: &
         sil_data_p      , &  ! field values at 2 temporal data points
         nit_data_p           ! field values at 2 temporal data points

      real (kind=dbl_kind) :: &
          secday      , &     ! number of seconds in day
          sec1hr              ! number of seconds in 1 hour

      logical (kind=log_kind) :: readm, read1, tr_bgc_Nit, tr_bgc_Sil

      character (char_len_long) :: &        ! input data file names
         nit_file   , & ! nitrate input file
         sil_file       ! silicate input file

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(get_forcing_bgc)'

      call icepack_query_parameters(secday_out=secday)
      call icepack_query_tracer_flags(tr_bgc_Nit_out=tr_bgc_Nit, &
           tr_bgc_Sil_out=tr_bgc_Sil)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (trim(bgc_data_type) == 'clim') then

         nit_file = 'nitrate_climatologyWOA_gx1v6f_20150107.nc'
                             !'nitrate_WOA2005_surface_monthly'  ! gx1 only
         sil_file = 'silicate_climatologyWOA_gx1v6f_20150107.nc'
                             !'silicate_WOA2005_surface_monthly' ! gx1 only

         nit_file = trim(bgc_data_dir)//'/'//trim(nit_file)
         sil_file = trim(bgc_data_dir)//'/'//trim(sil_file)

         if (my_task == master_task .and. istep == 1) then
         if (tr_bgc_Sil) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'silicate data interpolated to timestep:'
            write (nu_diag,*) trim(sil_file)
         endif
         if (tr_bgc_Nit) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'nitrate data interpolated to timestep:'
            write (nu_diag,*) trim(nit_file)
            if (restore_bgc) write (nu_diag,*) &
              'bgc restoring timescale (days) =', trestore
         endif
         endif                     ! my_task, istep

    !-------------------------------------------------------------------
    ! monthly data
    !
    ! Assume that monthly data values are located in the middle of the
    ! month.
    !-------------------------------------------------------------------

         midmonth = 15          ! data is given on 15th of every month
!!!      midmonth = fix(p5 * real(daymo(month)))  ! exact middle

         ! Compute record numbers for surrounding months
         maxrec = 12
         ixm  = mod(month+maxrec-2,maxrec) + 1
         ixp  = mod(month,         maxrec) + 1
         if (mday >= midmonth) ixm = -99 ! other two points will be used
         if (mday <  midmonth) ixp = -99

         ! Determine whether interpolation will use values 1:2 or 2:3
         ! recslot = 2 means we use values 1:2, with the current value (2)
         !  in the second slot
         ! recslot = 1 means we use values 2:3, with the current value (2)
         !  in the first slot
         recslot = 1            ! latter half of month
         if (mday < midmonth) recslot = 2 ! first half of month

         ! Find interpolation coefficients
         call interp_coeff_monthly (recslot)

         readm = .false.
         if (istep==1 .or. (mday==midmonth .and. sec==0)) readm = .true.

      endif   ! 'clim prep'

    !-------------------------------------------------------------------
    ! Read two monthly silicate values and interpolate.
    ! Restore toward interpolated value.
    !-------------------------------------------------------------------
    
      if (trim(bgc_data_type)=='clim'  .AND. tr_bgc_Sil) then
        ! call read_clim_data (readm, 0, ixm, month, ixp, &
        !                      sil_file,  sil_data, &
        !                      field_loc_center, field_type_scalar)
         fieldname = 'silicate'
         call read_clim_data_nc (readm, 0, ixm, month, ixp, &
                              sil_file, fieldname, sil_data, &
                              field_loc_center, field_type_scalar)
         call interpolate_data (sil_data, sildat)
         
         if (istep == 1 .or. .NOT. restore_bgc) then

            !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block)
      	    do iblk = 1, nblocks

               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi

	       do j = jlo, jhi
      	       do i = ilo, ihi

                  sil(i,j,iblk) = sildat(i,j,iblk)
                  ks = 2*icepack_max_algae + icepack_max_doc + 3 + icepack_max_dic
                  ocean_bio_all(i,j,ks,iblk) = sil(i,j,iblk)                       !Sil
               enddo
               enddo
         enddo

         !$OMP END PARALLEL DO
         elseif (restore_bgc) then

            !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block)
      	    do iblk = 1, nblocks

               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi

	       do j = jlo, jhi
      	       do i = ilo, ihi

                  sil(i,j,iblk) = sil(i,j,iblk)  &
                         + (sildat(i,j,iblk)-sil(i,j,iblk))*dt/trest
                  ks = 2*icepack_max_algae + icepack_max_doc + 3 + icepack_max_dic
                  ocean_bio_all(i,j,ks,iblk) = sil(i,j,iblk)                       !Sil
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
         endif  !restore
        elseif (tr_bgc_Sil) then ! bgc_data_type /= 'clim'
            !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block)
      	    do iblk = 1, nblocks

               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi

	       do j = jlo, jhi
      	       do i = ilo, ihi

	          sil(i,j,iblk) = 25.0_dbl_kind
                  ks = 2*icepack_max_algae + icepack_max_doc + 3 + icepack_max_dic
                  ocean_bio_all(i,j,ks,iblk) = sil(i,j,iblk)                       !Sil
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

      endif  !tr_bgc_Sil
    !-------------------------------------------------------------------
    ! Read two monthly nitrate values and interpolate.
    ! Restore toward interpolated value.
    !-------------------------------------------------------------------

      if (trim(bgc_data_type)=='clim' .AND. tr_bgc_Nit) then 
        ! call read_clim_data (readm, 0, ixm, month, ixp, &
        !                      nit_file, nit_data, &
        !                      field_loc_center, field_type_scalar)
         fieldname = 'nitrate'
         call read_clim_data_nc (readm, 0, ixm, month, ixp, &
                              nit_file, fieldname, nit_data, &
                              field_loc_center, field_type_scalar)
         call interpolate_data (nit_data, nitdat)

         if (istep == 1 .or. .NOT. restore_bgc) then
            !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block)
      	    do iblk = 1, nblocks

               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi

	       do j = jlo, jhi
      	       do i = ilo, ihi

	          nit(i,j,iblk) = nitdat(i,j,iblk)
                  ks = icepack_max_algae + 1
                  ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !nit
                  ks =  2*icepack_max_algae + icepack_max_doc + 7 + icepack_max_dic
                  ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !PON
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
         elseif (restore_bgc ) then
            !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block)
      	    do iblk = 1, nblocks

               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi

	       do j = jlo, jhi
      	       do i = ilo, ihi

                  nit(i,j,iblk) = nit(i,j,iblk)  &
                         + (nitdat(i,j,iblk)-nit(i,j,iblk))*dt/trest     
                  ks = icepack_max_algae + 1
                  ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !nit
                  ks =  2*icepack_max_algae + icepack_max_doc + 7 + icepack_max_dic
                  ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !PON
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
        endif  !restore_bgc

!      elseif (trim(nit_data_type) == 'sss'  .AND.  tr_bgc_Nit) then 
!           !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block)
!      	    do iblk = 1, nblocks

!               this_block = get_block(blocks_ice(iblk),iblk)
!               ilo = this_block%ilo
!               ihi = this_block%ihi
!               jlo = this_block%jlo
!               jhi = this_block%jhi

!               do j = jlo, jhi
!               do i = ilo, ihi

!                  nit(i,j,iblk) =  sss(i,j,iblk)      
!                  ks = icepack_max_algae + 1
!                  ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !nit 
!                  ks =  2*icepack_max_algae + icepack_max_doc + 7 + icepack_max_dic
!                  ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !PON      
!               enddo
!               enddo
!            enddo
!            !$OMP END PARALLEL DO

      elseif (tr_bgc_Nit) then ! bgc_data_type /= 'clim'
            !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block)
      	    do iblk = 1, nblocks

               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi

	       do j = jlo, jhi
      	       do i = ilo, ihi
       
                  nit(i,j,iblk) = 12.0_dbl_kind
                  ks = icepack_max_algae + 1
                  ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !nit 
                  ks =  2*icepack_max_algae + icepack_max_doc + 7 + icepack_max_dic
                  ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !PON      
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

      endif   !tr_bgc_Nit

    !-------------------------------------------------------------------
    ! Data from Papdimitrious et al., 2007, Limnol. Oceanogr. 
    ! and WOA at 68oS, 304.5oE : 
    ! daily data located at the end of the 24-hour period. 
    !-------------------------------------------------------------------

      if (trim(bgc_data_type) == 'ISPOL') then

         nit_file = trim(bgc_data_dir)//'nutrients_daily_ISPOL_WOA_field3.nc'
         sil_file = trim(bgc_data_dir)//'nutrients_daily_ISPOL_WOA_field3.nc' 

         if (my_task == master_task .and. istep == 1) then
         if (tr_bgc_Sil) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'silicate data interpolated to timestep:'
            write (nu_diag,*) trim(sil_file)
         endif
         if (tr_bgc_Nit) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'nitrate data interpolated to timestep:'
            write (nu_diag,*) trim(nit_file)
            if (restore_bgc) write (nu_diag,*) &
              'bgc restoring timescale (days) =', trestore
         endif
         endif                     ! my_task, istep

        dataloc = 2                          ! data located at end of interval
        sec1hr = secday                      ! seconds in day
        maxrec = 365                         ! 

        ! current record number
        recnum = int(yday)   

        ! Compute record numbers for surrounding data (2 on each side)
        ixm = mod(recnum+maxrec-2,maxrec) + 1
        ixx = mod(recnum-1,       maxrec) + 1
       
        recslot = 2
        ixp = -99
        call interp_coeff (recnum, recslot, sec1hr, dataloc)

        read1 = .false.
        if (istep==1 .or. bgcrecnum .ne. recnum) read1 = .true.
 
 
        if (tr_bgc_Sil) then
          met_file = sil_file
          fieldname= 'silicate' 
          call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
                    maxrec, met_file, fieldname, sil_data_p, &
                    field_loc_center, field_type_scalar)
      
          sil(:,:,:) = c1intp * sil_data_p(1) &
                     + c2intp * sil_data_p(2)
         endif

         if (tr_bgc_Nit) then
           met_file = nit_file
           fieldname= 'nitrate' 
           call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
                    maxrec, met_file, fieldname, nit_data_p, &
                    field_loc_center, field_type_scalar)
      
           nit(:,:,:) = c1intp * nit_data_p(1) &
                      + c2intp * nit_data_p(2)
         endif
         
            !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block)
      	    do iblk = 1, nblocks

               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi

	       do j = jlo, jhi
      	       do i = ilo, ihi
       
                  ks = 2*icepack_max_algae + icepack_max_doc + 3 + icepack_max_dic
                  ocean_bio_all(i,j,ks,iblk) = sil(i,j,iblk)                       !Sil  
                  ks = icepack_max_algae + 1
                  ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !nit
                  ks =  2*icepack_max_algae + icepack_max_doc + 7 + icepack_max_dic
                  ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !PON    
               enddo
               enddo
            enddo
          !$OMP END PARALLEL DO

       ! Save record number for next time step
       bgcrecnum = recnum
      endif

      end subroutine get_forcing_bgc

!=======================================================================
!
! author: Nicole Jeffery, LANL

      subroutine get_atm_bgc 

      use ice_blocks, only: block, get_block
      use ice_domain, only: nblocks, blocks_ice
      use ice_domain_size, only: n_zaero 
      use ice_flux_bgc, only: flux_bio_atm, faero_atm

      !  local variables

      integer (kind=int_kind) :: &
         i, j, nn       , & ! horizontal indices
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         iblk               ! block index 

      logical (kind=log_kind) :: &
         tr_zaero

      integer (kind=int_kind), dimension(icepack_max_aero) :: &
         nlt_zaero

      type (block) :: &
         this_block      ! block information for current block

      character(len=*), parameter :: subname = '(get_atm_bgc)'

      !-----------------------------------------------------------------
      ! initialize
      !-----------------------------------------------------------------

      call icepack_query_tracer_flags(tr_zaero_out=tr_zaero)
      call icepack_query_tracer_indices(nlt_zaero_out=nlt_zaero)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      flux_bio_atm(:,:,:,:) = c0
      if (tr_zaero) then
      !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,nn)
      do iblk = 1, nblocks

      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi
             
      do nn = 1, n_zaero 
         do j = jlo, jhi
         do i = ilo, ihi 
            flux_bio_atm(i,j,nlt_zaero(nn),iblk) = faero_atm(i,j,nn,iblk)
         enddo
         enddo
      enddo

      enddo ! iblk
      !$OMP END PARALLEL DO
      endif

      end subroutine get_atm_bgc

!=======================================================================

! constant values for atmospheric aerosols
!
! authors: Elizabeth Hunke, LANL

      subroutine faero_default

      use ice_flux_bgc, only: faero_atm

      character(len=*), parameter :: subname = '(faero_default)'

        faero_atm(:,:,1,:) = 1.e-12_dbl_kind ! kg/m^2 s
        faero_atm(:,:,2,:) = 1.e-13_dbl_kind
        faero_atm(:,:,3,:) = 1.e-14_dbl_kind 
        faero_atm(:,:,4,:) = 1.e-14_dbl_kind 
        faero_atm(:,:,5,:) = 1.e-14_dbl_kind 
        faero_atm(:,:,6,:) = 1.e-14_dbl_kind 

      end subroutine faero_default

!=======================================================================

! read atmospheric aerosols
!
! authors: Elizabeth Hunke, LANL

      subroutine faero_data

      use ice_calendar, only: month, mday, istep, sec
      use ice_domain_size, only: max_blocks
      use ice_blocks, only: nx_block, ny_block
      use ice_flux_bgc, only: faero_atm
      use ice_forcing, only: interp_coeff_monthly, read_clim_data_nc, interpolate_data

#ifdef ncdf 
      ! local parameters

      real (kind=dbl_kind), dimension(:,:,:,:), allocatable, &
         save :: &
         aero1_data    , & ! field values at 2 temporal data points
         aero2_data    , & ! field values at 2 temporal data points
         aero3_data        ! field values at 2 temporal data points

      character (char_len_long) :: & 
         aero_file,   &   ! netcdf filename
         fieldname        ! field name in netcdf file

      integer (kind=int_kind) :: & 
         ixm,ixp     , & ! record numbers for neighboring months
         maxrec      , & ! maximum record number
         recslot     , & ! spline slot for current record
         midmonth        ! middle day of month

      logical (kind=log_kind) :: readm

      character(len=*), parameter :: subname = '(faero_data)'

      allocate( aero1_data(nx_block,ny_block,2,max_blocks), &
                aero2_data(nx_block,ny_block,2,max_blocks), &
                aero3_data(nx_block,ny_block,2,max_blocks)  )


    !-------------------------------------------------------------------
    ! monthly data 
    !
    ! Assume that monthly data values are located in the middle of the 
    ! month.
    !-------------------------------------------------------------------

      midmonth = 15  ! data is given on 15th of every month
!      midmonth = fix(p5 * real(daymo(month)))  ! exact middle

      ! Compute record numbers for surrounding months
      maxrec = 12
      ixm  = mod(month+maxrec-2,maxrec) + 1
      ixp  = mod(month,         maxrec) + 1
      if (mday >= midmonth) ixm = 99  ! other two points will be used
      if (mday <  midmonth) ixp = 99

      ! Determine whether interpolation will use values 1:2 or 2:3
      ! recslot = 2 means we use values 1:2, with the current value (2)
      !  in the second slot
      ! recslot = 1 means we use values 2:3, with the current value (2)
      !  in the first slot
      recslot = 1                             ! latter half of month
      if (mday < midmonth) recslot = 2        ! first half of month

      ! Find interpolation coefficients
      call interp_coeff_monthly (recslot)

      ! Read 2 monthly values 
      readm = .false.
      if (istep==1 .or. (mday==midmonth .and. sec==0)) readm = .true.

!      aero_file = trim(atm_data_dir)//'faero.nc'   
      aero_file = '/usr/projects/climate/eclare/DATA/gx1v3/faero.nc'   

      fieldname='faero_atm001'
      call read_clim_data_nc (readm, 0,  ixm, month, ixp, &
                              aero_file, fieldname, aero1_data, &
                              field_loc_center, field_type_scalar)

      fieldname='faero_atm002'
      call read_clim_data_nc (readm, 0,  ixm, month, ixp, &
                              aero_file, fieldname, aero2_data, &
                              field_loc_center, field_type_scalar)

      fieldname='faero_atm003'
      call read_clim_data_nc (readm, 0,  ixm, month, ixp, &
                              aero_file, fieldname, aero3_data, &
                              field_loc_center, field_type_scalar)

      call interpolate_data (aero1_data, faero_atm(:,:,1,:)) ! W/m^2 s
      call interpolate_data (aero2_data, faero_atm(:,:,2,:))
      call interpolate_data (aero3_data, faero_atm(:,:,3,:))

      where (faero_atm(:,:,:,:) > 1.e20) faero_atm(:,:,:,:) = c0

      deallocate( aero1_data, aero2_data, aero3_data )
#endif

      end subroutine faero_data

!=======================================================================

! read atmospheric aerosols
!
! authors: Elizabeth Hunke, LANL

      subroutine fzaero_data

      use ice_blocks, only: nx_block, ny_block
      use ice_flux_bgc, only: faero_atm
      use ice_forcing, only: interp_coeff_monthly, read_clim_data_nc, interpolate_data

#ifdef ncdf 
      ! local parameters

      real (kind=dbl_kind), dimension(:,:,:,:), allocatable, &
         save :: &
         aero_data    ! field values at 2 temporal data points

      character (char_len_long) :: & 
         aero_file,   &   ! netcdf filename
         fieldname        ! field name in netcdf file

      integer (kind=int_kind) :: & 
         ixm,ixp     , & ! record numbers for neighboring months
         maxrec      , & ! maximum record number
         recslot     , & ! spline slot for current record
         midmonth        ! middle day of month

      integer (kind=int_kind), dimension(icepack_max_aero) :: &
         nlt_zaero

      logical (kind=log_kind) :: readm

      character(len=*), parameter :: subname = '(fzaero_data)'

      call icepack_query_tracer_indices(nlt_zaero_out=nlt_zaero)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      allocate( aero_data(nx_block,ny_block,2,max_blocks) )

    !-------------------------------------------------------------------
    ! monthly data 
    !
    ! Assume that monthly data values are located in the middle of the 
    ! month.
    !-------------------------------------------------------------------

      midmonth = 15  ! data is given on 15th of every month
!      midmonth = fix(p5 * real(daymo(month)))  ! exact middle

      ! Compute record numbers for surrounding months
      maxrec = 12
      ixm  = mod(month+maxrec-2,maxrec) + 1
      ixp  = mod(month,         maxrec) + 1
      if (mday >= midmonth) ixm = -99  ! other two points will be used
      if (mday <  midmonth) ixp = -99

      ! Determine whether interpolation will use values 1:2 or 2:3
      ! recslot = 2 means we use values 1:2, with the current value (2)
      !  in the second slot
      ! recslot = 1 means we use values 2:3, with the current value (2)
      !  in the first slot
      recslot = 1                             ! latter half of month
      if (mday < midmonth) recslot = 2        ! first half of month

      ! Find interpolation coefficients
      call interp_coeff_monthly (recslot)

      ! Read 2 monthly values 
      readm = .false.
      if (istep==1 .or. (mday==midmonth .and. sec==0)) readm = .true.

!      aero_file = trim(atm_data_dir)//'faero.nc'   
      ! Cam5 monthly total black carbon deposition on the gx1 grid"
      aero_file = '/usr/projects/climate/njeffery/DATA/CAM/Hailong_Wang/Cam5_bc_monthly_popgrid.nc'   

      fieldname='bcd'
      call read_clim_data_nc (readm, 0,  ixm, month, ixp, &
                              aero_file, fieldname, aero_data, &
                              field_loc_center, field_type_scalar)


      call interpolate_data (aero_data, faero_atm(:,:,nlt_zaero(1),:)) ! kg/m^2/s

      where (faero_atm(:,:,nlt_zaero(1),:) > 1.e20) faero_atm(:,:,nlt_zaero(1),:) = c0

      deallocate( aero_data )
#endif

      end subroutine fzaero_data

!=======================================================================

! Initialize ocean iron from file
!
! authors: Nicole Jeffery, LANL

      subroutine init_bgc_data (fed1,fep1)

      use ice_read_write, only: ice_open_nc, ice_read_nc, ice_close_nc

#ifdef ncdf
      use netcdf
#endif
           
      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(inout) :: &
           fed1, &  ! first dissolved iron pool (nM)
           fep1    ! first particulate iron pool (nM)

      ! local parameters

      integer (kind=int_kind) :: &
         fid              ! file id for netCDF file 

      logical (kind=log_kind) :: diag

      character (char_len_long) :: & 
         iron_file,   &   ! netcdf filename
         fieldname        ! field name in netcdf file

      character(len=*), parameter :: subname = '(init_bgc_data)'

    !-------------------------------------------------------------------
    ! Annual average data from Tagliabue, 2012 (top 50 m average
    ! poisson grid filled on gx1v6
    !-------------------------------------------------------------------

      if (trim(fe_data_type) == 'clim') then
       	diag = .true.   ! write diagnostic information 
        iron_file = trim(bgc_data_dir)//'dFe_50m_annual_Tagliabue_gx1.nc'

        if (my_task == master_task) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'Dissolved iron ocean concentrations from:'
            write (nu_diag,*) trim(iron_file)
            call ice_open_nc(iron_file,fid)
        endif

        fieldname='dFe'
        ! Currently only first fed  value is read
        call ice_read_nc(fid,1,fieldname,fed1,diag) 
        where ( fed1(:,:,:) > 1.e20) fed1(:,:,:) = p1  

        if (my_task == master_task) call ice_close_nc(fid)  

       	diag = .true.   ! write diagnostic information 
        iron_file = trim(bgc_data_dir)//'pFe_bathy_gx1.nc'

        if (my_task == master_task) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'Particulate iron ocean concentrations from:'
            write (nu_diag,*) trim(iron_file)
            call ice_open_nc(iron_file,fid)
        endif

        fieldname='pFe'
        ! Currently only first fep value is read
        call ice_read_nc(fid,1,fieldname,fep1,diag) 
        where ( fep1(:,:,:) > 1.e20) fep1(:,:,:) = p1  

        if (my_task == master_task) call ice_close_nc(fid)  
    
      endif
    
      end subroutine init_bgc_data

!=======================================================================
!
! Aerosol optical properties for bulk and modal aerosol formulation
! X_bc_tab properties are from snicar_optics_5bnd_mam_c140303 (Mark Flanner 2009)
! ==> "Mie optical parameters for CLM snowpack treatment" Includes
! ice (effective radii from 30-1500um), black carbon, organic carbon and dust
!
! authors: Elizabeth Hunke, LANL

      subroutine faero_optics

      use ice_broadcast, only: broadcast_array
      use ice_read_write, only: ice_open_nc, ice_close_nc
      use ice_communicate, only: my_task, master_task
      use ice_arrays_column, only: &
         kaer_tab, & ! aerosol mass extinction cross section (m2/kg)
         waer_tab, & ! aerosol single scatter albedo (fraction)
         gaer_tab, & ! aerosol asymmetry parameter (cos(theta))
         kaer_bc_tab, & ! BC mass extinction cross section (m2/kg)
         waer_bc_tab, & ! BC single scatter albedo (fraction)
         gaer_bc_tab, & ! BC aerosol asymmetry parameter (cos(theta))
         bcenh          ! BC absorption enhancement facto

#ifdef ncdf
      use netcdf
#endif

      ! local parameters

      integer (kind=int_kind) :: & 
         varid          , & ! variable id
         status         , & ! status output from netcdf routines
         n,  k              ! index

      real (kind=dbl_kind) :: &
         amin, amax, asum   ! min, max values and sum of input array

      integer (kind=int_kind) :: &
         fid                ! file id for netCDF file 

      logical (kind=log_kind) :: modal_aero

      character (char_len_long) :: & 
         optics_file,   &   ! netcdf filename
         fieldname          ! field name in netcdf file

      character(len=*), parameter :: subname = '(faero_optics)'

      ! this data is used in bulk aerosol treatment in dEdd radiation
      kaer_tab = reshape((/ &      ! aerosol mass extinction cross section (m2/kg)
          11580.61872,   5535.41835,   2793.79690, &
          25798.96479,  11536.03871,   4688.24207, &
            196.49772,    204.14078,    214.42287, &
           2665.85867,   2256.71027,    820.36024, &
            840.78295,   1028.24656,   1163.03298, &
            387.51211,    414.68808,    450.29814/), &
            (/icepack_nspint,icepack_max_aero/))
      waer_tab = reshape((/ &      ! aerosol single scatter albedo (fraction)
              0.29003,      0.17349,      0.06613, &
              0.51731,      0.41609,      0.21324, &
              0.84467,      0.94216,      0.95666, &
              0.97764,      0.99402,      0.98552, &
              0.94146,      0.98527,      0.99093, &
              0.90034,      0.96543,      0.97678/), &
              (/icepack_nspint,icepack_max_aero/))
      gaer_tab = reshape((/ &      ! aerosol asymmetry parameter (cos(theta))
              0.35445,      0.19838,      0.08857, &
              0.52581,      0.32384,      0.14970, &
              0.83162,      0.78306,      0.74375, &
              0.68861,      0.70836,      0.54171, &
              0.70239,      0.66115,      0.71983, &
              0.78734,      0.73580,      0.64411/), &
              (/icepack_nspint,icepack_max_aero/))

      ! this data is used in MODAL AEROSOL treatment in dEdd radiation
      kaer_bc_tab = reshape((/ &      ! aerosol mass extinction cross section (m2/kg)
             12955.44732,   5946.89461, 2772.33366, &
             12085.30664,   7438.83131, 3657.13084, &
              9753.99698,   7342.87139, 4187.79304, &
              7815.74879,   6659.65096, 4337.98863, &
              6381.28194,   5876.78408, 4254.65054, &
              5326.93163,   5156.74532, 4053.66581, &
              4538.09763,   4538.60875, 3804.10884, &
              3934.17604,   4020.20799, 3543.27199, &
              3461.20656,   3587.80962, 3289.98060, &
              3083.03396,   3226.27231, 3052.91441/), &
              (/icepack_nspint,10/))

      waer_bc_tab = reshape((/ &      ! aerosol single scatter albedo (fraction)
              0.26107,      0.15861,    0.06535, &
              0.37559,      0.30318,    0.19483, &
              0.42224,      0.36913,    0.27875, &
              0.44777,      0.40503,    0.33026, &
              0.46444,      0.42744,    0.36426, &
              0.47667,      0.44285,    0.38827, &
              0.48635,      0.45428,    0.40617, &
              0.49440,      0.46328,    0.42008, &
              0.50131,      0.47070,    0.43128, &
              0.50736,      0.47704,    0.44056/), &
              (/icepack_nspint,10/))

      gaer_bc_tab = reshape((/ &      ! aerosol asymmetry parameter (cos(theta))
              0.28328,      0.19644,      0.10498, &
              0.44488,      0.32615,      0.19612, &
              0.54724,      0.41611,      0.26390, &
              0.61711,      0.48475,      0.31922, &
              0.66673,      0.53923,      0.36632, &
              0.70296,      0.58337,      0.40732, &
              0.73002,      0.61960,      0.44344, &
              0.75064,      0.64959,      0.47551, &
              0.76663,      0.67461,      0.50415, &
              0.77926,      0.69561,      0.52981/),&
              (/icepack_nspint,10/))

      bcenh(:,:,:)     = c0

    call icepack_query_parameters(modal_aero_out=modal_aero)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
       file=__FILE__, line=__LINE__)

#ifdef ncdf
    if (modal_aero) then
       optics_file =  &
        '/usr/projects/climate/njeffery/DATA/CAM/snicar/snicar_optics_5bnd_mam_c140303.nc'

        if (my_task == master_task) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'Read optics for modal aerosol treament in'
            write (nu_diag,*) trim(optics_file)
            call ice_open_nc(optics_file,fid)
        endif

        fieldname='bcint_enh_mam_cice'
        if (my_task == master_task) then

          status = nf90_inq_varid(fid, trim(fieldname), varid)
 
           if (status /= nf90_noerr) then
             call abort_ice (subname//'ERROR: Cannot find variable '//trim(fieldname))
           endif
           status = nf90_get_var( fid, varid, bcenh, &
               start=(/1,1,1,1/), & 
               count=(/3,10,8,1/) )
           do n=1,10
            amin = minval(bcenh(:,n,:))
            amax = maxval(bcenh(:,n,:))
            asum = sum   (bcenh(:,n,:))
            write(nu_diag,*) ' min, max, sum =', amin, amax, asum
           enddo
           call ice_close_nc(fid)      
         endif  !master_task
         do n=1,3
            do k=1,8
                call broadcast_array(bcenh(n,:,k),      master_task)
            enddo
         enddo          
      endif      ! modal_aero
#else
    if (modal_aero) then
      call abort_ice(subname//'ERROR: netcdf required for modal_aero')
    endif
#endif

      end subroutine faero_optics

!=======================================================================

      end module ice_forcing_bgc

!=======================================================================
