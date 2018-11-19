!=======================================================================
!
! Drivers for remapping and upwind ice transport
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL 
!
! 2004: Revised by William Lipscomb from ice_transport_mpdata.
!       Stripped out mpdata, retained upwind, and added block structure.
! 2006: Incorporated remap transport driver and renamed from
!       ice_transport_upwind.  
! 2011: ECH moved edgearea arrays into ice_transport_remap.F90

      module ice_transport_driver

      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0, c1, p5, &
          field_loc_center, &
          field_type_scalar, field_type_vector, &
          field_loc_Nface, field_loc_Eface
      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_compute_tracers
      use icepack_intfc, only: icepack_query_tracer_flags, &
          icepack_query_tracer_numbers, icepack_query_tracer_indices, &
          icepack_query_parameters

      implicit none
      private
      public :: init_transport, transport_remap, transport_upwind

      character (len=char_len), public ::     &
         advection   ! type of advection scheme used
                     ! 'upwind' => 1st order donor cell scheme
                     ! 'remap' => remapping scheme

      logical, parameter :: & ! if true, prescribe area flux across each edge  
         l_fixed_area = .false.

! NOTE: For remapping, hice and hsno are considered tracers.
!       ntrace is not equal to ntrcr!

      integer (kind=int_kind) ::                      &
         ntrace              ! number of tracers in use
                          
      integer (kind=int_kind), dimension(:), allocatable ::             &
         tracer_type       ,&! = 1, 2, or 3 (depends on 0, 1 or 2 other tracers)
         depend              ! tracer dependencies (see below)

      logical (kind=log_kind), dimension (:), allocatable ::             &
         has_dependents      ! true if a tracer has dependent tracers

      integer (kind=int_kind), parameter ::                      &
         integral_order = 3   ! polynomial order of quadrature integrals
                              ! linear=1, quadratic=2, cubic=3

      logical (kind=log_kind), parameter ::     &
         l_dp_midpt = .true.  ! if true, find departure points using
                              ! corrected midpoint velocity
                          
!=======================================================================

      contains

!=======================================================================
!
! This subroutine is a wrapper for init_remap, which initializes the
! remapping transport scheme.  If the model is run with upwind
! transport, no initializations are necessary.
!
! authors William H. Lipscomb, LANL

      subroutine init_transport

      use ice_state, only: trcr_depend
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_advect
      use ice_transport_remap, only: init_remap

      integer (kind=int_kind) ::       &
         k, nt, nt1     ! tracer indices

      integer (kind=int_kind) :: ntrcr, nt_Tsfc, nt_qice, nt_qsno, &
          nt_sice, nt_fbri, nt_iage, nt_FY, nt_alvl, nt_vlvl, &
          nt_apnd, nt_hpnd, nt_ipnd, nt_bgc_Nit, nt_bgc_S

      character(len=*), parameter :: subname = '(init_transport)'

      call ice_timer_start(timer_advect)  ! advection 

      call icepack_query_tracer_numbers(ntrcr_out=ntrcr)
      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice, &
          nt_qsno_out=nt_qsno, nt_sice_out=nt_sice, nt_fbri_out=nt_fbri, &
          nt_iage_out=nt_iage, nt_FY_out=nt_FY, nt_alvl_out=nt_alvl, &
          nt_vlvl_out=nt_vlvl, nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
          nt_ipnd_out=nt_ipnd, nt_bgc_Nit_out=nt_bgc_Nit, nt_bgc_S_out=nt_bgc_S)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      ntrace = 2 + ntrcr ! hice,hsno,trcr

      if (allocated(tracer_type)) deallocate(tracer_type)
      if (allocated(depend)) deallocate(depend)
      if (allocated(has_dependents)) deallocate(has_dependents)

      allocate (tracer_type   (ntrace), &
                depend        (ntrace), &
                has_dependents(ntrace))

         ! define tracer dependency arrays
         ! see comments in remapping routine

          depend(1:2)         = 0 ! hice, hsno
          tracer_type(1:2)    = 1 ! no dependency
      
          k = 2

          do nt = 1, ntrcr
             depend(k+nt) = trcr_depend(nt) ! 0 for ice area tracers
                                            ! 1 for ice volume tracers
                                            ! 2 for snow volume tracers
             tracer_type(k+nt) = 2          ! depends on 1 other tracer
             if (trcr_depend(nt) == 0) then
                tracer_type(k+nt) = 1       ! depends on no other tracers
             elseif (trcr_depend(nt) > 2) then
                if (trcr_depend(trcr_depend(nt)-2) > 0) then
                   tracer_type(k+nt) = 3    ! depends on 2 other tracers
                endif
             endif
          enddo

          has_dependents = .false.
          do nt = 1, ntrace
             if (depend(nt) > 0) then
                nt1 = depend(nt)
                has_dependents(nt1) = .true.
                if (nt1 > nt) then
                   write(nu_diag,*)     &
                      'Tracer nt2 =',nt,' depends on tracer nt1 =',nt1
                   call abort_ice(subname//       &
                      'ERROR: remap transport: Must have nt2 > nt1')
                endif
             endif
          enddo                 ! ntrace

          ! diagnostic output
          if (my_task == master_task) then
          write (nu_diag, *) 'tracer        index      depend        type has_dependents'
             nt = 1
                write(nu_diag,*) '   hi  ',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             nt = 2
                write(nu_diag,*) '   hs  ',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
          k=2
          do nt = k+1, k+ntrcr
             if (nt-k==nt_Tsfc) &
                write(nu_diag,*) 'nt_Tsfc',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_qice) &
                write(nu_diag,*) 'nt_qice',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_qsno) &
                write(nu_diag,*) 'nt_qsno',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_sice) &
                write(nu_diag,*) 'nt_sice',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_fbri) &
                write(nu_diag,*) 'nt_fbri',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_iage) &
                write(nu_diag,*) 'nt_iage',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_FY) &
                write(nu_diag,*) 'nt_FY  ',  nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_alvl) &
                write(nu_diag,*) 'nt_alvl',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_vlvl) &
                write(nu_diag,*) 'nt_vlvl',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_apnd) &
                write(nu_diag,*) 'nt_apnd',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_hpnd) &
                write(nu_diag,*) 'nt_hpnd',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_ipnd) &
                write(nu_diag,*) 'nt_ipnd',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_bgc_Nit) &
                write(nu_diag,*) 'nt_bgc_Nit',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_bgc_S) &
                write(nu_diag,*) 'nt_bgc_S',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
          enddo
          endif ! master_task

          if (trim(advection)=='remap') call init_remap    ! grid quantities

      call ice_timer_stop(timer_advect)  ! advection 

      end subroutine init_transport

!=======================================================================
!
! This subroutine solves the transport equations for one timestep
! using the conservative remapping scheme developed by John Dukowicz
! and John Baumgardner (DB) and modified for sea ice by William
! Lipscomb and Elizabeth Hunke.
!
! This scheme preserves monotonicity of ice area and tracers.  That is,
! it does not produce new extrema.  It is second-order accurate in space,
! except where gradients are limited to preserve monotonicity. 
!
! authors William H. Lipscomb, LANL

      subroutine transport_remap (dt)

      use ice_blocks, only: nx_block, ny_block
      use ice_boundary, only: ice_HaloUpdate
      use ice_global_reductions, only: global_sum, global_sum_prod
      use ice_domain, only: nblocks, distrb_info, blocks_ice, halo_info
      use ice_domain_size, only: ncat, max_blocks
      use ice_blocks, only: nx_block, ny_block, block, get_block, nghost
      use ice_state, only: aice0, aicen, vicen, vsnon, trcrn, &
          uvel, vvel, bound_state
      use ice_grid, only: tarea
      use ice_calendar, only: istep1
      use ice_timers, only: ice_timer_start, ice_timer_stop, &
          timer_advect, timer_bound
      use ice_transport_remap, only: horizontal_remap, make_masks

      real (kind=dbl_kind), intent(in) ::     &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) ::     &
         iblk           ,&! block index
         ilo,ihi,jlo,jhi,&! beginning and end of physical domain
         n              ,&! ice category index
         nt, nt1, nt2     ! tracer indices

      real (kind=dbl_kind),      &
         dimension (nx_block,ny_block,0:ncat,max_blocks) ::     &
         aim            ,&! mean ice category areas in each grid cell
         aimask           ! = 1. if ice is present, = 0. otherwise

      real (kind=dbl_kind),      &
         dimension (nx_block,ny_block,ntrace,ncat,max_blocks) ::     &
         trm            ,&! mean tracer values in each grid cell
         trmask           ! = 1. if tracer is present, = 0. otherwise

      logical (kind=log_kind) ::     &
         l_stop           ! if true, abort the model

      integer (kind=int_kind) ::     &
         istop, jstop     ! indices of grid cell where model aborts 

      integer (kind=int_kind), dimension(0:ncat,max_blocks) ::     &
         icellsnc         ! number of cells with ice

      integer (kind=int_kind),      &
         dimension(nx_block*ny_block,0:ncat,max_blocks) ::     &
         indxinc, indxjnc   ! compressed i/j indices

      integer (kind=int_kind) :: &
         ntrcr

      type (block) :: &
         this_block           ! block information for current block
      
      ! variables related to optional bug checks

      logical (kind=log_kind), parameter ::     &
         l_conservation_check = .false. ,&! if true, check conservation
         l_monotonicity_check = .false.   ! if true, check monotonicity

      real (kind=dbl_kind), dimension(0:ncat) ::     &
         asum_init      ,&! initial global ice area
         asum_final       ! final global ice area

      real (kind=dbl_kind), dimension(ntrace,ncat) ::     &
         atsum_init     ,&! initial global ice area*tracer
         atsum_final      ! final global ice area*tracer

      real (kind=dbl_kind), dimension (:,:,:,:,:), allocatable ::     &
         tmin         ,&! local min tracer
         tmax           ! local max tracer

      integer (kind=int_kind) :: alloc_error

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      character(len=*), parameter :: subname = '(transport_remap)'

      call ice_timer_start(timer_advect)  ! advection 
      call icepack_query_tracer_numbers(ntrcr_out=ntrcr)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

!---!-------------------------------------------------------------------
!---! Prepare for remapping.
!---! Initialize, update ghost cells, fill tracer arrays.
!---!-------------------------------------------------------------------

      l_stop = .false.
      istop = 0
      jstop = 0

    !-------------------------------------------------------------------
    ! Compute open water area in each grid cell.
    ! Note: An aggregate_area call is needed only if the open
    !       water area has changed since the previous call.
    !       Here we assume that aice0 is up to date.
    !-------------------------------------------------------------------

!      !$OMP PARALLEL DO PRIVATE(i,j,iblk)
!      do iblk = 1, nblocks
!      do j = 1, ny_block
!      do i = 1, nx_block
!         call aggregate_area (ncat,
!                              aicen(i,j,:,iblk),     &
!                              aice (i,j,  iblk),     &
!                              aice0(i,j,  iblk)) 
!      enddo
!      enddo
!      enddo
!      !$OMP END PARALLEL DO

    !-------------------------------------------------------------------
    ! Ghost cell updates for state variables.
    ! Commented out because ghost cells are updated after cleanup_itd.
    !-------------------------------------------------------------------
!      call ice_timer_start(timer_bound)

!      call ice_HaloUpdate (aice0,            halo_info,     &
!                           field_loc_center, field_type_scalar)

!      call bound_state (aicen,        &
!                        vicen, vsnon, &
!                        ntrcr, trcrn)

!      call ice_timer_stop(timer_bound)

    !-------------------------------------------------------------------
    ! Ghost cell updates for ice velocity.
    ! Commented out because ghost cell velocities are computed
    !  in ice_dyn_evp.
    !-------------------------------------------------------------------

!      call ice_timer_start(timer_bound)
!      call ice_HaloUpdate (uvel,               halo_info,     &
!                           field_loc_NEcorner, field_type_vector)
!      call ice_HaloUpdate (vvel,               halo_info,     &
!                           field_loc_NEcorner, field_type_vector)
!      call ice_timer_stop(timer_bound)


      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks

    !-------------------------------------------------------------------
    ! Fill arrays with fields to be remapped.
    !-------------------------------------------------------------------

         call state_to_tracers(nx_block,          ny_block,             &
                               ntrcr,             ntrace,               &
                               aice0(:,:,  iblk), aicen(:,:,:,iblk),    &
                               trcrn(:,:,:,:,iblk),                     &
                               vicen(:,:,:,iblk), vsnon(:,:,  :,iblk),  &
                               aim  (:,:,:,iblk), trm  (:,:,:,:,iblk))

      enddo
      !$OMP END PARALLEL DO

!---!-------------------------------------------------------------------
!---! Optional conservation and monotonicity checks.
!---!-------------------------------------------------------------------

      if (l_conservation_check) then

    !-------------------------------------------------------------------
    ! Compute initial values of globally conserved quantities.
    !-------------------------------------------------------------------

         do n = 0, ncat
            asum_init(n) = global_sum(aim(:,:,n,:),     distrb_info,       &
                                      field_loc_center, tarea)
         enddo

         do n = 1, ncat
            do nt = 1, ntrace
               if (tracer_type(nt)==1) then ! does not depend on another tracer
                  atsum_init(nt,n) =      &
                      global_sum_prod(trm(:,:,nt,n,:), aim(:,:,n,:),       &
                                      distrb_info,     field_loc_center,   &
                                      tarea)
               elseif (tracer_type(nt)==2) then ! depends on another tracer
                  nt1 = depend(nt)
                  work1(:,:,:) = trm(:,:,nt,n,:)*trm(:,:,nt1,n,:)
                  atsum_init(nt,n) =     &
                      global_sum_prod(work1(:,:,:), aim(:,:,n,:),          &
                                      distrb_info,  field_loc_center,      &
                                      tarea)
               elseif (tracer_type(nt)==3) then ! depends on two tracers
                  nt1 = depend(nt)
                  nt2 = depend(nt1)
                  work1(:,:,:) = trm(:,:,nt,n,:)*trm(:,:,nt1,n,:)          &
                                                *trm(:,:,nt2,n,:)
                  atsum_init(nt,n) =     &
                      global_sum_prod(work1(:,:,:), aim(:,:,n,:),          &
                                      distrb_info,  field_loc_center,      &
                                      tarea)
               endif            ! tracer_type
            enddo               ! nt
         enddo                  ! n

      endif                     ! l_conservation_check
      
      if (l_monotonicity_check) then

         allocate(tmin(nx_block,ny_block,ntrace,ncat,max_blocks),     &
                  tmax(nx_block,ny_block,ntrace,ncat,max_blocks),     &
                  STAT=alloc_error)

         if (alloc_error /= 0)      &
              call abort_ice (subname//'ERROR: allocation error')

         tmin(:,:,:,:,:) = c0
         tmax(:,:,:,:,:) = c0

         !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block,n)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

    !------------------------------------------------------------------- 
    ! Compute masks.
    ! Masks are used to prevent tracer values in cells without ice
    !  from being used in the monotonicity check.
    !------------------------------------------------------------------- 

            call make_masks (nx_block,          ny_block,              &
                             ilo, ihi,          jlo, jhi,              &
                             nghost,            ntrace,                &
                             has_dependents,                           &
                             icellsnc(:,iblk),                         &
                             indxinc(:,:,iblk), indxjnc(:,:,iblk),     &
                             aim(:,:,:,iblk),   aimask(:,:,:,iblk),    &
                             trm(:,:,:,:,iblk), trmask(:,:,:,:,iblk))

    !-------------------------------------------------------------------
    ! Compute local max and min of tracer fields.
    !-------------------------------------------------------------------

            do n = 1, ncat
               call local_max_min                                      &  
                            (nx_block,           ny_block,             &
                             ilo, ihi,           jlo, jhi,             &
                             trm (:,:,:,n,iblk),                       &
                             tmin(:,:,:,n,iblk), tmax  (:,:,:,n,iblk), &
                             aimask(:,:,n,iblk), trmask(:,:,:,n,iblk))
            enddo
         enddo
         !$OMP END PARALLEL DO

         call ice_timer_start(timer_bound)
         call ice_HaloUpdate (tmin,             halo_info,     &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (tmax,             halo_info,     &
                              field_loc_center, field_type_scalar)
         call ice_timer_stop(timer_bound)

         !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block,n)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do n = 1, ncat
               call quasilocal_max_min (nx_block, ny_block,     &
                                        ilo, ihi, jlo, jhi,     &
                                        tmin(:,:,:,n,iblk),      &
                                        tmax(:,:,:,n,iblk))
            enddo
         enddo
         !$OMP END PARALLEL DO

      endif                     ! l_monotonicity_check

    !-------------------------------------------------------------------
    ! Main remapping routine: Step ice area and tracers forward in time.
    !-------------------------------------------------------------------
   
         call horizontal_remap (dt,                ntrace,             &
                                uvel      (:,:,:), vvel      (:,:,:),  &
                                aim     (:,:,:,:), trm   (:,:,:,:,:),  &
                                l_fixed_area,                          &
                                tracer_type,       depend,             &
                                has_dependents,    integral_order,     &
                                l_dp_midpt)
         
    !-------------------------------------------------------------------
    ! Given new fields, recompute state variables.
    !-------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks

         call tracers_to_state (nx_block,          ny_block,            &
                                ntrcr,             ntrace,              &
                                aim  (:,:,:,iblk), trm  (:,:,:,:,iblk), &
                                aice0(:,:,  iblk), aicen(:,:,:,iblk),   &
                                trcrn(:,:,:,:,iblk),                    &
                                vicen(:,:,:,iblk), vsnon(:,:,  :,iblk))

      enddo                     ! iblk
      !$OMP END PARALLEL DO

    !-------------------------------------------------------------------
    ! Ghost cell updates for state variables.
    !-------------------------------------------------------------------

      call ice_timer_start(timer_bound)

      call bound_state (aicen,        &
                        vicen, vsnon, &
                        ntrcr, trcrn)

      call ice_timer_stop(timer_bound)

!---!-------------------------------------------------------------------
!---! Optional conservation and monotonicity checks
!---!-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! Compute final values of globally conserved quantities.
    ! Check global conservation of area and area*tracers.  (Optional)
    !-------------------------------------------------------------------

      if (l_conservation_check) then

         do n = 0, ncat
            asum_final(n) = global_sum(aim(:,:,n,:),     distrb_info,      &
                                       field_loc_center, tarea)
         enddo

         do n = 1, ncat
            do nt = 1, ntrace
               if (tracer_type(nt)==1) then ! does not depend on another tracer
                  atsum_final(nt,n) =      &
                      global_sum_prod(trm(:,:,nt,n,:), aim(:,:,n,:),       &
                                      distrb_info,     field_loc_center,   &
                                      tarea)
               elseif (tracer_type(nt)==2) then ! depends on another tracer
                  nt1 = depend(nt)
                  work1(:,:,:) = trm(:,:,nt,n,:)*trm(:,:,nt1,n,:)
                  atsum_final(nt,n) =     &
                      global_sum_prod(work1(:,:,:), aim(:,:,n,:),          &
                                      distrb_info,  field_loc_center,      &
                                      tarea)
               elseif (tracer_type(nt)==3) then ! depends on two tracers
                  nt1 = depend(nt)
                  nt2 = depend(nt1)
                  work1(:,:,:) = trm(:,:,nt,n,:)*trm(:,:,nt1,n,:)          &
                                                *trm(:,:,nt2,n,:)
                  atsum_final(nt,n) =     &
                      global_sum_prod(work1(:,:,:), aim(:,:,n,:),          &
                                      distrb_info,  field_loc_center,      &
                                      tarea)
               endif            ! tracer_type
            enddo               ! nt
         enddo                  ! n

         if (my_task == master_task) then
            call global_conservation (l_stop,     &
                                      asum_init(0), asum_final(0))

            if (l_stop) then
               write (nu_diag,*) 'istep1, my_task, iblk =',     &
                                  istep1, my_task, iblk
               write (nu_diag,*) 'transport: conservation error, cat 0'
               call abort_ice(subname//'ERROR: conservation error1')
            endif

            do n = 1, ncat               
               call global_conservation                                 &
                                     (l_stop,                           &
                                      asum_init(n),    asum_final(n),   &
                                      atsum_init(:,n), atsum_final(:,n))

               if (l_stop) then
                  write (nu_diag,*) 'istep1, my_task, iblk, cat =',     &
                                     istep1, my_task, iblk, n
                  write (nu_diag,*) 'transport: conservation error, cat ',n
                  call abort_ice(subname//'ERROR: conservation error2')
               endif
            enddo               ! n

         endif                  ! my_task = master_task

      endif                     ! l_conservation_check

    !-------------------------------------------------------------------
    ! Check tracer monotonicity.  (Optional)
    !-------------------------------------------------------------------

      if (l_monotonicity_check) then
         !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block,n,l_stop,istop,jstop)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            l_stop = .false.
            istop = 0
            jstop = 0

            do n = 1, ncat
               call check_monotonicity      &
                               (nx_block,           ny_block,     &
                                ilo, ihi, jlo, jhi,     &
                                tmin(:,:,:,n,iblk), tmax(:,:,:,n,iblk),  &
                                aim (:,:,  n,iblk), trm (:,:,:,n,iblk),  &
                                l_stop,     &
                                istop,              jstop)

               if (l_stop) then
                  write (nu_diag,*) 'istep1, my_task, iblk, cat =',     &
                                     istep1, my_task, iblk, n
                  call abort_ice(subname//'ERROR: monotonicity error')
               endif
            enddo               ! n

         enddo                  ! iblk
         !$OMP END PARALLEL DO

         deallocate(tmin, tmax, STAT=alloc_error)
         if (alloc_error /= 0) call abort_ice (subname//'ERROR: deallocation error')

      endif                     ! l_monotonicity_check

      call ice_timer_stop(timer_advect)  ! advection 
           
      end subroutine transport_remap

!=======================================================================
!
! Computes the transport equations for one timestep using upwind. Sets
! several fields into a work array and passes it to upwind routine.

      subroutine transport_upwind (dt)

      use ice_boundary, only: ice_HaloUpdate
      use ice_blocks, only: nx_block, ny_block, block, get_block, nx_block, ny_block
      use ice_domain, only: blocks_ice, halo_info, nblocks
      use ice_domain_size, only: ncat, max_blocks
      use ice_state, only: aice0, aicen, vicen, vsnon, trcrn, &
          uvel, vvel, trcr_depend, bound_state, trcr_base, &
          n_trcr_strata, nt_strata
      use ice_grid, only: HTE, HTN, tarea
      use ice_timers, only: ice_timer_start, ice_timer_stop, &
          timer_bound, timer_advect

      real (kind=dbl_kind), intent(in) ::     &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) ::     &
         ntrcr, &           !
         narr               ! max number of state variable arrays

      integer (kind=int_kind) ::     &
         i, j, iblk       ,&! horizontal indices
         ilo,ihi,jlo,jhi    ! beginning and end of physical domain

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblocks) ::     &
         uee, vnn           ! cell edge velocities

      real (kind=dbl_kind),     &
         dimension (:,:,:,:), allocatable :: &
         works              ! work array

      type (block) ::     &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(transport_upwind)'

      call ice_timer_start(timer_advect)  ! advection 

      call icepack_query_tracer_numbers(ntrcr_out=ntrcr)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      narr = 1 + ncat*(3+ntrcr) ! max number of state variable arrays

      allocate (works(nx_block,ny_block,narr,max_blocks))

    !-------------------------------------------------------------------
    ! Get ghost cell values of state variables.
    ! (Assume velocities are already known for ghost cells, also.)
    !-------------------------------------------------------------------
!      call bound_state (aicen,        &
!                        vicen, vsnon, &
!                        ntrcr, trcrn)

    !-------------------------------------------------------------------
    ! Average corner velocities to edges.
    !-------------------------------------------------------------------
      
      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            uee(i,j,iblk) = p5*(uvel(i,j,iblk) + uvel(i,j-1,iblk))
            vnn(i,j,iblk) = p5*(vvel(i,j,iblk) + vvel(i-1,j,iblk))
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (uee,             halo_info,     &
                           field_loc_Eface, field_type_vector)
      call ice_HaloUpdate (vnn,             halo_info,     &
                           field_loc_Nface, field_type_vector)
      call ice_timer_stop(timer_bound)

      !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi


      !-----------------------------------------------------------------
      ! fill work arrays with fields to be advected
      !-----------------------------------------------------------------

         call state_to_work (nx_block,             ny_block,             &
                             ntrcr,                                      &
                             narr,                 trcr_depend,          &
                             aicen (:,:,  :,iblk), trcrn (:,:,:,:,iblk), &
                             vicen (:,:,  :,iblk), vsnon (:,:,  :,iblk), &
                             aice0 (:,:,    iblk), works (:,:,  :,iblk))

      !-----------------------------------------------------------------
      ! advect
      !-----------------------------------------------------------------

         call upwind_field (nx_block,       ny_block,               &
                            ilo, ihi,       jlo, jhi,               &
                            dt,                                     &
                            narr,           works(:,:,:,iblk),      &
                            uee(:,:,iblk),  vnn    (:,:,iblk),      &
                            HTE(:,:,iblk),  HTN    (:,:,iblk),      &
                            tarea(:,:,iblk))

      !-----------------------------------------------------------------
      ! convert work arrays back to state variables
      !-----------------------------------------------------------------

         call work_to_state (nx_block,            ny_block,             &
                             ntrcr,               narr,                 &
                             trcr_depend(:),      trcr_base(:,:),       &
                             n_trcr_strata(:),    nt_strata(:,:),       &
                             aicen(:,:,  :,iblk), trcrn (:,:,:,:,iblk), &
                             vicen(:,:,  :,iblk), vsnon (:,:,  :,iblk), &
                             aice0(:,:,    iblk), works (:,:,  :,iblk)) 

      enddo                     ! iblk
      !$OMP END PARALLEL DO
 
      deallocate (works)

    !-------------------------------------------------------------------
    ! Ghost cell updates for state variables.
    !-------------------------------------------------------------------

      call ice_timer_start(timer_bound)

      call bound_state (aicen,        &
                        vicen, vsnon, &
                        ntrcr, trcrn)

      call ice_timer_stop(timer_bound)

      call ice_timer_stop(timer_advect)  ! advection 

      end subroutine transport_upwind

!=======================================================================
! The next few subroutines (through check_monotonicity) are called
! by transport_remap.
!=======================================================================
!
! Fill ice area and tracer arrays.
! Assume that the advected tracers are hicen, hsnon, trcrn, 
!  qicen(1:nilyr), and qsnon(1:nslyr).
! This subroutine must be modified if a different set of tracers
!   is to be transported.  The rule for ordering tracers
!   is that a dependent tracer (such as qice) must have a larger
!   tracer index than the tracer it depends on (i.e., hice).
!
! author William H. Lipscomb, LANL

      subroutine state_to_tracers (nx_block, ny_block,   &
                                   ntrcr,    ntrace,     &
                                   aice0,    aicen,      &
                                   trcrn,                &
                                   vicen,    vsnon,      &
                                   aim,      trm)

      use ice_domain_size, only: ncat, nslyr

      integer (kind=int_kind), intent(in) ::     &
           nx_block, ny_block, & ! block dimensions
           ntrcr             , & ! number of tracers in use
           ntrace                ! number of tracers in use incl. hi, hs

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
           aice0     ! fractional open water area

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), intent(in) :: &
           aicen   ,&! fractional ice area
           vicen   ,&! volume per unit area of ice          (m)
           vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat), intent(in) :: &
           trcrn     ! ice area tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block,0:ncat), intent(out) :: &
           aim       ! mean ice area in each grid cell

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrace,ncat), intent(out) :: &
           trm       ! mean tracer values in each grid cell

      ! local variables

      integer (kind=int_kind) ::     &
           nt_qsno      ,&!
           i, j, n      ,&! standard indices
           it, kt       ,&! tracer indices
           ij             ! combined i/j index

      real (kind=dbl_kind) ::     &
           puny         ,&!
           rhos         ,&!
           Lfresh       ,&!
           w1             ! work variable

      integer (kind=int_kind), dimension(nx_block*ny_block,0:ncat) ::  &
           indxi        ,&! compressed i/j indices
           indxj

      integer (kind=int_kind), dimension(0:ncat) ::     &
           icells         ! number of cells with ice

      character(len=*), parameter :: subname = '(state_to_tracers)'

      call icepack_query_parameters(puny_out=puny, rhos_out=rhos, &
           Lfresh_out=Lfresh)
      call icepack_query_tracer_indices(nt_qsno_out=nt_qsno)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      aim(:,:,0) = aice0(:,:)

      do n = 1, ncat

         trm(:,:,:,n) = c0

    !-------------------------------------------------------------------
    ! Find grid cells where ice is present and fill area array.
    !-------------------------------------------------------------------

         icells(n) = 0
         do j = 1, ny_block
         do i = 1, nx_block
            aim(i,j,n) = aicen(i,j,n)
            if (aim(i,j,n) > puny) then
               icells(n) = icells(n) + 1
               ij = icells(n)
               indxi(ij,n) = i
               indxj(ij,n) = j
            endif               ! aim > puny
         enddo
         enddo
      
    !-------------------------------------------------------------------
    ! Fill tracer array
    ! Note: If aice > 0, then hice > 0, but we can have hsno = 0.
    ! Alse note: We transport qice*nilyr rather than qice, so as to
    !  avoid extra operations here and in tracers_to_state.
    !-------------------------------------------------------------------

         do ij = 1, icells(n)
            i = indxi(ij,n)
            j = indxj(ij,n)
            w1 = c1 / aim(i,j,n)
            trm(i,j,1,n) = vicen(i,j,n) * w1 ! hice
            trm(i,j,2,n) = vsnon(i,j,n) * w1 ! hsno
         enddo
         kt = 2

         do it = 1, ntrcr
            if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
               do ij = 1, icells(n)
                  i = indxi(ij,n)
                  j = indxj(ij,n)
                  trm(i,j,kt+it,n) = trcrn(i,j,it,n) + rhos*Lfresh ! snow enthalpy
               enddo
            else
               do ij = 1, icells(n)
                  i = indxi(ij,n)
                  j = indxj(ij,n)
                  trm(i,j,kt+it,n) = trcrn(i,j,it,n) ! other tracers
               enddo
            endif
         enddo
      enddo                     ! ncat
 
      end subroutine state_to_tracers

!=======================================================================
!
! Convert area and tracer arrays back to state variables.
!
! author William H. Lipscomb, LANL

      subroutine tracers_to_state (nx_block, ny_block,   &
                                   ntrcr,    ntrace,     &
                                   aim,      trm,        &
                                   aice0,    aicen,      &
                                   trcrn,                &
                                   vicen,    vsnon)

      use ice_domain_size, only: ncat, nslyr

      integer (kind=int_kind), intent(in) ::     &
           nx_block, ny_block, & ! block dimensions
           ntrcr             , & ! number of tracers in use
           ntrace                ! number of tracers in use incl. hi, hs

      real (kind=dbl_kind), dimension (nx_block,ny_block,0:ncat), intent(in) :: &
           aim       ! fractional ice area

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrace,ncat), intent(in) :: &
           trm       ! mean tracer values in each grid cell

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
           aice0     ! fractional ice area

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), intent(inout) :: &
           aicen   ,&! fractional ice area
           vicen   ,&! volume per unit area of ice          (m)
           vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat), intent(inout) :: &
           trcrn     ! tracers

      ! local variables

      integer (kind=int_kind) ::     &
           nt_qsno         ,&!
           i, j, n      ,&! standard indices
           it, kt          ,&! tracer indices
           icells          ,&! number of cells with ice
           ij

      real (kind=dbl_kind) :: &
           rhos, &
           Lfresh

      integer (kind=int_kind), dimension (nx_block*ny_block) ::     &
           indxi, indxj      ! compressed indices

      character(len=*), parameter :: subname = '(tracers_to_state)'

      call icepack_query_parameters(rhos_out=rhos, Lfresh_out=Lfresh)
      call icepack_query_tracer_indices(nt_qsno_out=nt_qsno)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      aice0(:,:) = aim(:,:,0)

      do n = 1, ncat

      icells = 0
      do j = 1, ny_block
      do i = 1, nx_block
         if (aim(i,j,n) > c0) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif
      enddo
      enddo

    !-------------------------------------------------------------------
    ! Compute state variables.
    !-------------------------------------------------------------------

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            aicen(i,j,n) = aim(i,j,n)
            vicen(i,j,n) = aim(i,j,n)*trm(i,j,1,n) ! aice*hice
            vsnon(i,j,n) = aim(i,j,n)*trm(i,j,2,n) ! aice*hsno
         enddo                  ! ij
         kt = 2

         do it = 1, ntrcr
            if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  trcrn(i,j,it,n) = trm(i,j,kt+it,n) - rhos*Lfresh ! snow enthalpy
               enddo
               else
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  trcrn(i,j,it,n) = trm(i,j,kt+it,n)  ! other tracers
               enddo
            endif
         enddo
      enddo                     ! ncat

      end subroutine tracers_to_state

!=======================================================================
!
! Check whether values of conserved quantities have changed.
! An error probably means that ghost cells are treated incorrectly.
!
! author William H. Lipscomb, LANL

      subroutine global_conservation (l_stop,                     &
                                      asum_init,  asum_final,     &
                                      atsum_init, atsum_final)

      real (kind=dbl_kind), intent(in) ::     &
         asum_init   ,&! initial global ice area
         asum_final    ! final global ice area

      real (kind=dbl_kind), dimension(ntrace), intent(in), optional :: &
         atsum_init  ,&! initial global ice area*tracer
         atsum_final   ! final global ice area*tracer

      logical (kind=log_kind), intent(inout) ::     &
         l_stop    ! if true, abort on return

      ! local variables

      integer (kind=int_kind) ::     &
           nt            ! tracer index

      real (kind=dbl_kind) ::     &
           puny        ,&!
           diff          ! difference between initial and final values

      character(len=*), parameter :: subname = '(global_conservation)'

      call icepack_query_parameters(puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (asum_init > puny) then
         diff = asum_final - asum_init
         if (abs(diff/asum_init) > puny) then
            l_stop = .true.
            write (nu_diag,*)
            write (nu_diag,*) 'Ice area conserv error'
            write (nu_diag,*) 'Initial global area =', asum_init
            write (nu_diag,*) 'Final global area =', asum_final
            write (nu_diag,*) 'Fractional error =', abs(diff)/asum_init
            write (nu_diag,*) 'asum_final-asum_init =', diff
         endif
      endif

      if (present(atsum_init)) then
       do nt = 1, ntrace
         if (abs(atsum_init(nt)) > puny) then
            diff = atsum_final(nt) - atsum_init(nt)
            if (abs(diff/atsum_init(nt)) > puny) then
               l_stop = .true.
               write (nu_diag,*)
               write (nu_diag,*) 'area*tracer conserv error'
               write (nu_diag,*) 'tracer index =', nt
               write (nu_diag,*) 'Initial global area*tracer =',   &
                                  atsum_init(nt)
               write (nu_diag,*) 'Final global area*tracer =',     &
                                  atsum_final(nt)
               write (nu_diag,*) 'Fractional error =',             &
                                  abs(diff)/atsum_init(nt)
               write (nu_diag,*) 'atsum_final-atsum_init =', diff
            endif
         endif
       enddo
      endif                     ! present(atsum_init)

      end subroutine global_conservation

!=======================================================================
!
! At each grid point, compute the local max and min of a scalar
! field phi: i.e., the max and min values in the nine-cell region
! consisting of the home cell and its eight neighbors.
! 
! To extend to the neighbors of the neighbors (25 cells in all),
! follow this call with a call to quasilocal_max_min.
!
! author William H. Lipscomb, LANL

      subroutine local_max_min (nx_block, ny_block,     &
                                ilo, ihi, jlo, jhi,     &
                                trm,                    &
                                tmin,     tmax,         &
                                aimask,   trmask)

      integer (kind=int_kind), intent(in) ::     &
           nx_block, ny_block,&! block dimensions
           ilo,ihi,jlo,jhi     ! beginning and end of physical domain

      real (kind=dbl_kind), intent(in), dimension(nx_block,ny_block) :: &
           aimask         ! ice area mask

      real (kind=dbl_kind), intent(in), dimension (nx_block,ny_block,ntrace) :: &
           trm          ,&! tracer fields
           trmask         ! tracer mask

      real (kind=dbl_kind), intent(out), dimension (nx_block,ny_block,ntrace) :: &
           tmin         ,&! local min tracer
           tmax           ! local max tracer

      ! local variables

      integer (kind=int_kind) ::     &
           i, j         ,&! horizontal indices
           nt, nt1        ! tracer indices

      real (kind=dbl_kind), dimension(nx_block,ny_block) ::     &
           phimask        ! aimask or trmask, as appropriate

      real (kind=dbl_kind) ::     &
           phi_nw, phi_n, phi_ne ,&! field values in 8 neighbor cells
           phi_w, phi_e          ,&
           phi_sw, phi_s, phi_se

      character(len=*), parameter :: subname = '(local_max_min)'

      do nt = 1, ntrace

         if (tracer_type(nt)==1) then  ! does not depend on another tracer

            do j = 1, ny_block
            do i = 1, nx_block
               phimask(i,j) = aimask(i,j)
            enddo
            enddo

         else   ! depends on another tracer

            nt1 = depend(nt)
            do j = 1, ny_block
            do i = 1, nx_block
               phimask(i,j) = trmask(i,j,nt1)
            enddo
            enddo

         endif

!-----------------------------------------------------------------------
!  Store values of trm in the 8 neighbor cells.
!  If aimask = 1, use the true value; otherwise use the home cell value
!  so that non-physical values of phi do not contribute to the gradient.
!-----------------------------------------------------------------------

         do j = jlo, jhi
            do i = ilo, ihi

               phi_nw = phimask(i-1,j+1) * trm(i-1,j+1,nt)     &
                  + (c1-phimask(i-1,j+1))* trm(i,  j,  nt)
               phi_n  = phimask(i,  j+1) * trm(i,  j+1,nt)     &
                  + (c1-phimask(i,  j+1))* trm(i,  j,  nt)
               phi_ne = phimask(i+1,j+1) * trm(i+1,j+1,nt)     &
                  + (c1-phimask(i+1,j+1))* trm(i,  j,  nt)
               phi_w  = phimask(i-1,j)   * trm(i-1,j,  nt)     &
                  + (c1-phimask(i-1,j))  * trm(i,  j,  nt)
               phi_e  = phimask(i+1,j)   * trm(i+1,j,  nt)     &
                  + (c1-phimask(i+1,j))  * trm(i,  j,  nt)
               phi_sw = phimask(i-1,j-1) * trm(i-1,j-1,nt)     &
                  + (c1-phimask(i-1,j-1))* trm(i,  j,  nt)
               phi_s  = phimask(i,  j-1) * trm(i,  j-1,nt)     &
                  + (c1-phimask(i,  j-1))* trm(i,  j,  nt)
               phi_se = phimask(i+1,j-1) * trm(i+1,j-1,nt)     &
                  + (c1-phimask(i+1,j-1))* trm(i,  j,  nt)

!-----------------------------------------------------------------------
!     Compute the minimum and maximum among the nine local cells.
!-----------------------------------------------------------------------

               tmax(i,j,nt) = max (phi_nw, phi_n,  phi_ne, phi_w,     &
                      trm(i,j,nt), phi_e,  phi_sw, phi_s,  phi_se)

               tmin(i,j,nt) = min (phi_nw, phi_n,  phi_ne, phi_w,     &
                      trm(i,j,nt), phi_e,  phi_sw, phi_s,  phi_se)

            enddo               ! i
         enddo                  ! j

      enddo                     ! nt

      end subroutine local_max_min

!=======================================================================
!
! Extend the local max and min by one grid cell in each direction.
! Incremental remapping is monotone for the "quasilocal" max and min,
! but in rare cases may violate monotonicity for the local max and min.
!
! author William H. Lipscomb, LANL

      subroutine quasilocal_max_min (nx_block, ny_block,     &
                                     ilo, ihi, jlo, jhi,     &
                                     tmin,     tmax)

      integer (kind=int_kind), intent(in) ::     &
         nx_block, ny_block,&! block dimensions
         ilo,ihi,jlo,jhi     ! beginning and end of physical domain

      real (kind=dbl_kind), intent(inout), dimension (nx_block,ny_block,ntrace) :: &
           tmin         ,&! local min tracer
           tmax           ! local max tracer

      ! local variables

      integer (kind=int_kind) ::     &
           i, j          ,&! horizontal indices
           nt              ! tracer index

      character(len=*), parameter :: subname = '(quasilocal_max_min)'

      do nt = 1, ntrace

         do j = jlo, jhi
         do i = ilo, ihi

            tmax(i,j,nt) =     &
              max (tmax(i-1,j+1,nt), tmax(i,j+1,nt), tmax(i+1,j+1,nt),     &
                   tmax(i-1,j,  nt), tmax(i,j,  nt), tmax(i+1,j,  nt),     &
                   tmax(i-1,j-1,nt), tmax(i,j-1,nt), tmax(i+1,j-1,nt))

            tmin(i,j,nt) =     &
              min (tmin(i-1,j+1,nt), tmin(i,j+1,nt), tmin(i+1,j+1,nt),     &
                   tmin(i-1,j,  nt), tmin(i,j,  nt), tmin(i+1,j,  nt),     &
                   tmin(i-1,j-1,nt), tmin(i,j-1,nt), tmin(i+1,j-1,nt))

         enddo                  ! i
         enddo                  ! j

      enddo

      end subroutine quasilocal_max_min

!======================================================================
!
! At each grid point, make sure that the new tracer values
! fall between the local max and min values before transport.
!
! author William H. Lipscomb, LANL

      subroutine check_monotonicity (nx_block, ny_block,     &
                                     ilo, ihi, jlo, jhi,     &
                                     tmin,     tmax,         &
                                     aim,      trm,          &
                                     l_stop,                 &
                                     istop,    jstop)

      integer (kind=int_kind), intent(in) ::     &
           nx_block, ny_block,&! block dimensions
           ilo,ihi,jlo,jhi     ! beginning and end of physical domain

      real (kind=dbl_kind), intent(in), dimension (nx_block,ny_block) ::     &
           aim            ! new ice area

      real (kind=dbl_kind), intent(in), dimension (nx_block,ny_block,ntrace) ::     &
           trm            ! new tracers

      real (kind=dbl_kind), intent(in), dimension (nx_block,ny_block,ntrace) ::     &
           tmin         ,&! local min tracer
           tmax           ! local max tracer

      logical (kind=log_kind), intent(inout) ::     &
         l_stop    ! if true, abort on return

      integer (kind=int_kind), intent(inout) ::     &
         istop, jstop     ! indices of grid cell where model aborts 

      ! local variables

      integer (kind=int_kind) ::     &
           i, j           ,&! horizontal indices
           nt, nt1, nt2     ! tracer indices

      real (kind=dbl_kind) ::     &
           puny         ,&!
           w1, w2         ! work variables

      logical (kind=log_kind), dimension (nx_block, ny_block) ::   &
           l_check        ! if true, check monotonicity

      character(len=*), parameter :: subname = '(check_monotonicity)'

      call icepack_query_parameters(puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do nt = 1, ntrace

    !-------------------------------------------------------------------
    ! Load logical array to identify tracers that need checking.
    !-------------------------------------------------------------------

         if (tracer_type(nt)==1) then ! does not depend on another tracer

            do j = jlo, jhi
            do i = ilo, ihi
               if (aim(i,j) > puny) then 
                  l_check(i,j) = .true.
               else
                  l_check(i,j) = .false.
               endif
            enddo
            enddo

         elseif (tracer_type(nt)==2) then ! depends on another tracer

            nt1 = depend(nt)
            do j = jlo, jhi
            do i = ilo, ihi
               if (abs(trm(i,j,nt1)) > puny) then
                  l_check(i,j) = .true.
               else
                  l_check(i,j) = .false.
               endif
            enddo
            enddo

         elseif (tracer_type(nt)==3) then ! depends on two tracers

            nt1 = depend(nt)
            nt2 = depend(nt1)
            do j = jlo, jhi
            do i = ilo, ihi
               if (abs(trm(i,j,nt1)) > puny .and.     &
                   abs(trm(i,j,nt2)) > puny) then
                  l_check(i,j) = .true.
               else
                  l_check(i,j) = .false.
               endif
            enddo
            enddo
         endif

    !-------------------------------------------------------------------
    ! Make sure new values lie between tmin and tmax
    !-------------------------------------------------------------------

         do j = jlo, jhi
         do i = ilo, ihi

            if (l_check(i,j)) then
               ! w1 and w2 allow for roundoff error when abs(trm) is big
               w1 = max(c1, abs(tmin(i,j,nt)))
               w2 = max(c1, abs(tmax(i,j,nt)))
               if (trm(i,j,nt) < tmin(i,j,nt)-w1*puny) then
                  l_stop = .true.
                  istop = i
                  jstop = j
                  write (nu_diag,*) ' '
                  write (nu_diag,*) 'new tracer < tmin'
                  write (nu_diag,*) 'i, j, nt =', i, j, nt
                  write (nu_diag,*) 'new tracer =', trm (i,j,nt)
                  write (nu_diag,*) 'tmin ='      , tmin(i,j,nt)
                  write (nu_diag,*) 'ice area ='  , aim(i,j)
               elseif (trm(i,j,nt) > tmax(i,j,nt)+w2*puny) then
                  l_stop = .true.
                  istop = i
                  jstop = j
                  write (nu_diag,*) ' '
                  write (nu_diag,*) 'new tracer > tmax'
                  write (nu_diag,*) 'i, j, nt =', i, j, nt
                  write (nu_diag,*) 'new tracer =', trm (i,j,nt)
                  write (nu_diag,*) 'tmax ='      , tmax(i,j,nt)
                  write (nu_diag,*) 'ice area ='  , aim(i,j)
               endif
            endif

         enddo                  ! i
         enddo                  ! j

      enddo                     ! nt

      end subroutine check_monotonicity

!=======================================================================
! The remaining subroutines are called by transport_upwind.
!=======================================================================
!
! Fill work array with state variables in preparation for upwind transport

      subroutine state_to_work (nx_block, ny_block,        &
                                ntrcr,                     &
                                narr,     trcr_depend,     &
                                aicen,    trcrn,           &
                                vicen,    vsnon,           &
                                aice0,    works)

      use ice_domain_size, only: ncat

      integer (kind=int_kind), intent(in) ::     &
         nx_block, ny_block, & ! block dimensions
         ntrcr             , & ! number of tracers in use
         narr        ! number of 2D state variable arrays in works array

      integer (kind=int_kind), dimension (ntrcr), intent(in) ::     &
         trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), intent(in) ::     &
         aicen   ,&! concentration of ice
         vicen   ,&! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat), intent(in) ::     &
         trcrn     ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::        &
         aice0     ! concentration of open water

      real (kind=dbl_kind), dimension(nx_block,ny_block,narr), intent (out) ::      &
         works     ! work array

      ! local variables

      integer (kind=int_kind) :: &
         nt_alvl, nt_apnd, nt_fbri

      logical (kind=log_kind) :: &
         tr_pond_cesm, tr_pond_lvl, tr_pond_topo

      integer (kind=int_kind) ::      &
         i, j, n, it    ,&! counting indices
         narrays          ! counter for number of state variable arrays

      character(len=*), parameter :: subname = '(state_to_work)'

      call icepack_query_tracer_flags(tr_pond_cesm_out=tr_pond_cesm, &
           tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
      call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_apnd_out=nt_apnd, &
           nt_fbri_out=nt_fbri)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! This array is used for performance (balance memory/cache vs
      ! number of bound calls);  a different number of arrays may perform
      ! better depending on the machine used, number of processors, etc.
      ! --tested on SGI R2000, using 4 pes for the ice model under MPI
      !-----------------------------------------------------------------

      do j = 1, ny_block
      do i = 1, nx_block
         works(i,j,1) = aice0(i,j)
      enddo
      enddo
      narrays = 1

      do n=1, ncat

         do j = 1, ny_block
         do i = 1, nx_block
            works(i,j,narrays+1) = aicen(i,j,n)
            works(i,j,narrays+2) = vicen(i,j,n)
            works(i,j,narrays+3) = vsnon(i,j,n)
         enddo                  ! i
         enddo                  ! j
         narrays = narrays + 3

         do it = 1, ntrcr
            if (trcr_depend(it) == 0) then
               do j = 1, ny_block
               do i = 1, nx_block
                  works(i,j,narrays+it) = aicen(i,j,n)*trcrn(i,j,it,n)
               enddo
               enddo
            elseif (trcr_depend(it) == 1) then
               do j = 1, ny_block
               do i = 1, nx_block
                  works(i,j,narrays+it) = vicen(i,j,n)*trcrn(i,j,it,n)
               enddo
               enddo
            elseif (trcr_depend(it) == 2) then
               do j = 1, ny_block
               do i = 1, nx_block
                  works(i,j,narrays+it) = vsnon(i,j,n)*trcrn(i,j,it,n)
               enddo
               enddo
            elseif (trcr_depend(it) == 2+nt_alvl) then
               do j = 1, ny_block
               do i = 1, nx_block
                  works(i,j,narrays+it) = aicen(i,j,n) &
                                        * trcrn(i,j,nt_alvl,n) &
                                        * trcrn(i,j,it,n)
               enddo
               enddo
            elseif (trcr_depend(it) == 2+nt_apnd .and. &
                    tr_pond_cesm .or. tr_pond_topo) then
               do j = 1, ny_block
               do i = 1, nx_block
                  works(i,j,narrays+it) = aicen(i,j,n) &
                                        * trcrn(i,j,nt_apnd,n) &
                                        * trcrn(i,j,it,n)
               enddo
               enddo
            elseif (trcr_depend(it) == 2+nt_apnd .and. &
                    tr_pond_lvl) then
               do j = 1, ny_block
               do i = 1, nx_block
                  works(i,j,narrays+it) = aicen(i,j,n) &
                                        * trcrn(i,j,nt_alvl,n) &
                                        * trcrn(i,j,nt_apnd,n) &
                                        * trcrn(i,j,it,n)
               enddo
               enddo
            elseif (trcr_depend(it) == 2+nt_fbri) then
               do j = 1, ny_block
               do i = 1, nx_block
                  works(i,j,narrays+it) = vicen(i,j,n) &
                                        * trcrn(i,j,nt_fbri,n) &
                                        * trcrn(i,j,it,n)
               enddo
               enddo
            endif
         enddo
         narrays = narrays + ntrcr

      enddo                     ! n

      if (narr /= narrays) write(nu_diag,*)      &
           "Wrong number of arrays in transport bound call"

      end subroutine state_to_work

!=======================================================================
!
! Convert work array back to state variables

      subroutine work_to_state (nx_block, ny_block, &
                                ntrcr,    narr,     &
                                trcr_depend,        &
                                trcr_base,          &
                                n_trcr_strata,      &
                                nt_strata,          &
                                aicen,    trcrn,    &
                                vicen,    vsnon,    &
                                aice0,    works)

      use ice_domain_size, only: ncat

      integer (kind=int_kind), intent (in) ::                       &
         nx_block, ny_block, & ! block dimensions
         ntrcr             , & ! number of tracers in use
         narr        ! number of 2D state variable arrays in works array

      integer (kind=int_kind), dimension (ntrcr), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (ntrcr,3), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (ntrcr,2), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      real (kind=dbl_kind), intent (in) ::                          &
         works (nx_block,ny_block,narr)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), intent(out) :: &
         aicen   ,&! concentration of ice
         vicen   ,&! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat),intent(out) :: &
         trcrn     ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) :: &
         aice0     ! concentration of open water

      ! local variables

      integer (kind=int_kind) ::      &
         i, j, ij, n    ,&! counting indices
         narrays        ,&! counter for number of state variable arrays
         icells           ! number of ocean/ice cells

      integer (kind=int_kind), dimension (nx_block*ny_block) ::        &
        indxi, indxj

      real (kind=dbl_kind), dimension (nx_block*ny_block,narr) ::      &
         work 

      character(len=*), parameter :: subname = '(work_to_state)'

      ! for call to compute_tracers
      icells = 0
      do j = 1, ny_block
      do i = 1, nx_block
         icells = icells + 1
         indxi(icells) = i
         indxj(icells) = j
         work (icells,:) = works(i,j,:)
      enddo
      enddo

      do j=1,ny_block
      do i=1,nx_block
         aice0(i,j) = works(i,j,1)
      enddo
      enddo
      narrays = 1               ! aice0 is first array

      do n=1,ncat

         do j=1,ny_block
         do i=1,nx_block
            aicen(i,j,n) = works(i,j,narrays+1)
            vicen(i,j,n) = works(i,j,narrays+2)
            vsnon(i,j,n) = works(i,j,narrays+3)
         enddo
         enddo
         narrays = narrays + 3

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            call icepack_compute_tracers (ntrcr,        trcr_depend(:),     &
                                         work (ij,narrays+1:narrays+ntrcr), &
                                         aicen(i,j,n),                     &
                                         vicen(i,j,n), vsnon(i,j,n),       &
                                         trcr_base(:,:), n_trcr_strata(:), &
                                         nt_strata(:,:), &
                                         trcrn(i,j,:,n))
         enddo
         narrays = narrays + ntrcr

      enddo                     ! ncat

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      end subroutine work_to_state

!=======================================================================
!
! upwind transport algorithm

      subroutine upwind_field (nx_block, ny_block,   &
                               ilo, ihi, jlo, jhi,   &
                               dt,                   &
                               narrays,  phi,        &
                               uee,      vnn,        &
                               HTE,      HTN,        &
                               tarea)

      integer (kind=int_kind), intent (in) ::     &
         nx_block, ny_block ,&! block dimensions
         ilo,ihi,jlo,jhi    ,&! beginning and end of physical domain
         narrays              ! number of 2D arrays to be transported

      real (kind=dbl_kind), intent(in) ::         &
         dt                   ! time step

      real (kind=dbl_kind), dimension(nx_block,ny_block,narrays), intent(inout) :: &
         phi                  ! scalar field

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
         uee, vnn             ! cell edge velocities

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
         HTE                ,&! length of east cell edge 
         HTN                ,&! length of north cell edge
         tarea                ! grid cell area

      ! local variables

      integer (kind=int_kind) :: &
         i, j, n              ! standard indices

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         worka, workb

      character(len=*), parameter :: subname = '(upwind_field)'

    !-------------------------------------------------------------------
    ! upwind transport
    !-------------------------------------------------------------------

      do n = 1, narrays

         do j = 1, jhi
         do i = 1, ihi
            worka(i,j)=     &
               upwind(phi(i,j,n),phi(i+1,j,n),uee(i,j),HTE(i,j),dt)
            workb(i,j)=     &
               upwind(phi(i,j,n),phi(i,j+1,n),vnn(i,j),HTN(i,j),dt)
         enddo
         enddo

         do j = jlo, jhi
         do i = ilo, ihi
            phi(i,j,n) = phi(i,j,n) - ( worka(i,j)-worka(i-1,j)      &
                                      + workb(i,j)-workb(i,j-1) )    &
                                      / tarea(i,j)
         enddo
         enddo

      enddo                     ! narrays

      end subroutine upwind_field

!=======================================================================

    !-------------------------------------------------------------------
    ! Define upwind function
    !-------------------------------------------------------------------

      real(kind=dbl_kind) function upwind(y1,y2,a,h,dt)

      real(kind=dbl_kind), intent(in) :: y1,y2,a,h,dt

      upwind = p5*dt*h*((a+abs(a))*y1+(a-abs(a))*y2)

      end function upwind

!=======================================================================

      end module ice_transport_driver

!=======================================================================
