!=======================================================================
!
! Elastic-anisotropic sea ice dynamics model
! Computes ice velocity and deformation
!
! See:
!
! Wilchinsky, A.V. and D.L. Feltham (2006). Modelling the rheology of 
! sea ice as a collection of diamond-shaped floes. 
! Journal of Non-Newtonian Fluid Mechanics, 138(1), 22-32.
!
! Tsamados, M., D.L. Feltham, and A.V. Wilchinsky (2013). Impact on new
! anisotropic rheology on simulations of Arctic sea ice. JGR, 118, 91-107.
!
! authors: Michel Tsamados, CPOM 
!          David Schroeder, CPOM

      module ice_dyn_eap

      use ice_kinds_mod
      use ice_blocks, only: nx_block, ny_block
      use ice_communicate, only: my_task, master_task
      use ice_domain_size, only: max_blocks, ncat
      use ice_constants, only: c0, c1, c2, c3, c4, c12, p1, p2, p5, &
          p001, p027, p055, p111, p166, p222, p25, p333
      use ice_fileunits, only: nu_diag, nu_dump_eap, nu_restart_eap
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_ice_strength

      implicit none
      private
      public :: eap, init_eap, write_restart_eap, read_restart_eap, &
                alloc_dyn_eap

      ! Look-up table needed for calculating structure tensor
      integer (int_kind), parameter :: & 
        nx_yield            =  41, &
        ny_yield            =  41, &
        na_yield            =  21

      real (kind=dbl_kind), dimension (nx_yield,ny_yield,na_yield) :: & 
        s11r, s12r, s22r, s11s, s12s, s22s           

      real (kind=dbl_kind), dimension (:,:,:), allocatable :: &
         a11_1, a11_2, a11_3, a11_4,                  & ! components of 
         a12_1, a12_2, a12_3, a12_4                     ! structure tensor

      ! history
      real (kind=dbl_kind), dimension(:,:,:), allocatable, public :: &
         e11      , & ! components of strain rate tensor (1/s)
         e12      , & 
         e22      , &
         yieldstress11, & ! components of yield stress tensor (kg/s^2)
         yieldstress12, &
         yieldstress22, &
         s11      , & ! components of stress tensor (kg/s^2)
         s12      , &
         s22      , &
         a11      , & ! components of structure tensor ()
         a12

!=======================================================================

      contains

!=======================================================================
!
! Allocate space for all variables 
!
      subroutine alloc_dyn_eap

      integer (int_kind) :: ierr

      allocate( a11_1        (nx_block,ny_block,max_blocks), &
                a11_2        (nx_block,ny_block,max_blocks), &
                a11_3        (nx_block,ny_block,max_blocks), &
                a11_4        (nx_block,ny_block,max_blocks), &
                a12_1        (nx_block,ny_block,max_blocks), &
                a12_2        (nx_block,ny_block,max_blocks), &
                a12_3        (nx_block,ny_block,max_blocks), &
                a12_4        (nx_block,ny_block,max_blocks), &
                e11          (nx_block,ny_block,max_blocks), &
                e12          (nx_block,ny_block,max_blocks), &
                e22          (nx_block,ny_block,max_blocks), &
                yieldstress11(nx_block,ny_block,max_blocks), &
                yieldstress12(nx_block,ny_block,max_blocks), &
                yieldstress22(nx_block,ny_block,max_blocks), &
                s11          (nx_block,ny_block,max_blocks), &
                s12          (nx_block,ny_block,max_blocks), &
                s22          (nx_block,ny_block,max_blocks), &
                a11          (nx_block,ny_block,max_blocks), &
                a12          (nx_block,ny_block,max_blocks), &
                stat=ierr)
      if (ierr/=0) call abort_ice('(alloc_dyn_eap): Out of memory')

      end subroutine alloc_dyn_eap

!=======================================================================
!
! Elastic-anisotropic-plastic dynamics driver
! based on subroutine evp

      subroutine eap (dt)

#ifdef CICE_IN_NEMO
! Wind stress is set during this routine from the values supplied
! via NEMO (unless calc_strair is true).  These values are supplied  
! rotated on u grid and multiplied by aice.  strairxT = 0 in this  
! case so operations in dyn_prep1 are pointless but carried out to  
! minimise code changes.
#endif

      use ice_arrays_column, only: Cdn_ocn
      use ice_boundary, only: ice_halo, ice_HaloMask, ice_HaloUpdate, &
          ice_HaloDestroy
      use ice_blocks, only: block, get_block
      use ice_constants, only: field_loc_center, field_loc_NEcorner, &
          field_type_scalar, field_type_vector
      use ice_domain, only: nblocks, blocks_ice, halo_info, maskhalo_dyn
      use ice_dyn_shared, only: fcor_blk, ndte, dtei, &
          denom1, uvel_init, vvel_init, arlx1i, &
          dyn_prep1, dyn_prep2, stepu, dyn_finish, &
          seabed_stress_factor_LKD, seabed_stress_factor_prob, &
          seabed_stress_method, seabed_stress, &
          stack_velocity_field, unstack_velocity_field
      use ice_flux, only: rdg_conv, strairxT, strairyT, &
          strairx, strairy, uocn, vocn, ss_tltx, ss_tlty, iceumask, fm, &
          strtltx, strtlty, strocnx, strocny, strintx, strinty, taubx, tauby, &
          strocnxT, strocnyT, strax, stray, &
          Tbu, hwater, &
          stressp_1, stressp_2, stressp_3, stressp_4, &
          stressm_1, stressm_2, stressm_3, stressm_4, &
          stress12_1, stress12_2, stress12_3, stress12_4
      use ice_grid, only: tmask, umask, dxt, dyt, dxhy, dyhx, cxp, cyp, cxm, cym, &
          tarear, uarear, to_ugrid, t2ugrid_vector, u2tgrid_vector
      use ice_state, only: aice, vice, vsno, uvel, vvel, divu, shear, &
          aice_init, aice0, aicen, vicen, strength
!      use ice_timers, only: timer_dynamics, timer_bound, &
!          ice_timer_start, ice_timer_stop, &
!          timer_tmp1, timer_tmp2, timer_tmp3
      use ice_timers, only: timer_dynamics, timer_bound, &
          ice_timer_start, ice_timer_stop

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: & 
         ksub           , & ! subcycle step
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, ij

      integer (kind=int_kind), dimension(max_blocks) :: & 
         icellt   , & ! no. of cells where icetmask = 1
         icellu       ! no. of cells where iceumask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block, max_blocks) :: &
         indxti   , & ! compressed index in i-direction
         indxtj   , & ! compressed index in j-direction
         indxui   , & ! compressed index in i-direction
         indxuj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         tmass    , & ! total mass of ice and snow (kg/m^2)
         waterx   , & ! for ocean stress calculation, x (m/s)
         watery   , & ! for ocean stress calculation, y (m/s)
         forcex   , & ! work array: combined atm stress and ocn tilt, x
         forcey   , & ! work array: combined atm stress and ocn tilt, y
         aiu      , & ! ice fraction on u-grid
         umass    , & ! total mass of ice and snow (u grid)
         umassdti     ! mass of U-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), allocatable :: fld2(:,:,:,:)

      real (kind=dbl_kind), dimension(nx_block,ny_block,8):: &
         strtmp       ! stress combinations for momentum equation

      logical (kind=log_kind) :: calc_strair

      integer (kind=int_kind), dimension (nx_block,ny_block,max_blocks) :: &
         icetmask, &  ! ice extent mask (T-cell)
         halomask     ! ice mask for halo update

      type (ice_halo) :: &
         halo_info_mask !  ghost cell update info for masked halo

      type (block) :: &
         this_block           ! block information for current block
      
      character(len=*), parameter :: subname = '(eap)'

      call ice_timer_start(timer_dynamics) ! dynamics

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      allocate(fld2(nx_block,ny_block,2,max_blocks))

       ! This call is needed only if dt changes during runtime.
!      call set_evp_parameters (dt)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         do j = 1, ny_block 
         do i = 1, nx_block 
            rdg_conv (i,j,iblk) = c0 
!           rdg_shear(i,j,iblk) = c0 
            divu (i,j,iblk) = c0 
            shear(i,j,iblk) = c0 
            e11(i,j,iblk) = c0
            e12(i,j,iblk) = c0
            e22(i,j,iblk) = c0
            s11(i,j,iblk) = c0
            s12(i,j,iblk) = c0
            s22(i,j,iblk) = c0
            yieldstress11(i,j,iblk) = c0
            yieldstress12(i,j,iblk) = c0
            yieldstress22(i,j,iblk) = c0
         enddo
         enddo

      !-----------------------------------------------------------------
      ! preparation for dynamics
      !-----------------------------------------------------------------

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         call dyn_prep1 (nx_block,           ny_block,           & 
                         ilo, ihi,           jlo, jhi,           &
                         aice    (:,:,iblk), vice    (:,:,iblk), & 
                         vsno    (:,:,iblk), tmask   (:,:,iblk), & 
                         strairxT(:,:,iblk), strairyT(:,:,iblk), & 
                         strairx (:,:,iblk), strairy (:,:,iblk), & 
                         tmass   (:,:,iblk), icetmask(:,:,iblk))

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (icetmask,          halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! convert fields from T to U grid
      !-----------------------------------------------------------------

      call to_ugrid(tmass,umass)
      call to_ugrid(aice_init, aiu)

      !----------------------------------------------------------------
      ! Set wind stress to values supplied via NEMO or other forcing
      ! This wind stress is rotated on u grid and multiplied by aice
      !----------------------------------------------------------------
      call icepack_query_parameters(calc_strair_out=calc_strair)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (.not. calc_strair) then       
         strairx(:,:,:) = strax(:,:,:)
         strairy(:,:,:) = stray(:,:,:)
      else
         call t2ugrid_vector(strairx)
         call t2ugrid_vector(strairy)
      endif

! tcraig, tcx, turned off this threaded region, in evp, this block and 
! the icepack_ice_strength call seems to not be thread safe.  more
! debugging needed
      !$TCXOMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! more preparation for dynamics
      !-----------------------------------------------------------------

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         call dyn_prep2 (nx_block,             ny_block,             & 
                         ilo, ihi,             jlo, jhi,             &
                         icellt(iblk),         icellu(iblk),         & 
                         indxti      (:,iblk), indxtj      (:,iblk), & 
                         indxui      (:,iblk), indxuj      (:,iblk), & 
                         aiu       (:,:,iblk), umass     (:,:,iblk), & 
                         umassdti  (:,:,iblk), fcor_blk  (:,:,iblk), & 
                         umask     (:,:,iblk),                       & 
                         uocn      (:,:,iblk), vocn      (:,:,iblk), & 
                         strairx   (:,:,iblk), strairy   (:,:,iblk), & 
                         ss_tltx   (:,:,iblk), ss_tlty   (:,:,iblk), &  
                         icetmask  (:,:,iblk), iceumask  (:,:,iblk), & 
                         fm        (:,:,iblk), dt,                   & 
                         strtltx   (:,:,iblk), strtlty   (:,:,iblk), & 
                         strocnx   (:,:,iblk), strocny   (:,:,iblk), & 
                         strintx   (:,:,iblk), strinty   (:,:,iblk), & 
                         taubx     (:,:,iblk), tauby     (:,:,iblk), & 
                         waterx    (:,:,iblk), watery    (:,:,iblk), & 
                         forcex    (:,:,iblk), forcey    (:,:,iblk), & 
                         stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), & 
                         stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), & 
                         stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), & 
                         stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), & 
                         stress12_1(:,:,iblk), stress12_2(:,:,iblk), & 
                         stress12_3(:,:,iblk), stress12_4(:,:,iblk), & 
                         uvel_init (:,:,iblk), vvel_init (:,:,iblk), &
                         uvel      (:,:,iblk), vvel      (:,:,iblk), &
                         Tbu       (:,:,iblk))

      !-----------------------------------------------------------------
      ! Initialize structure tensor
      !-----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            if (icetmask(i,j,iblk)==0) then
            ! structure tensor
               a11_1(i,j,iblk) = p5
               a11_2(i,j,iblk) = p5
               a11_3(i,j,iblk) = p5
               a11_4(i,j,iblk) = p5
               a12_1(i,j,iblk) = c0
               a12_2(i,j,iblk) = c0
               a12_3(i,j,iblk) = c0
               a12_4(i,j,iblk) = c0
            endif                  ! icetmask
         enddo                     ! i
         enddo                     ! j

      !-----------------------------------------------------------------
      ! ice strength
      ! New strength used in Ukita Moritz rheology         
      !-----------------------------------------------------------------

         strength(:,:,iblk) = c0  ! initialize
         do ij = 1, icellt(iblk)
            i = indxti(ij, iblk)
            j = indxtj(ij, iblk)
            call icepack_ice_strength(ncat=ncat,                 &
                                      aice     = aice    (i,j,  iblk), & 
                                      vice     = vice    (i,j,  iblk), & 
                                      aice0    = aice0   (i,j,  iblk), & 
                                      aicen    = aicen   (i,j,:,iblk), &  
                                      vicen    = vicen   (i,j,:,iblk), & 
                                      strength = strength(i,j,  iblk) )
         enddo  ! ij
      enddo  ! iblk
      !$TCXOMP END PARALLEL DO

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (strength,           halo_info, &
                           field_loc_center,   field_type_scalar)
      ! velocities may have changed in dyn_prep2
      call stack_velocity_field(uvel, vvel, fld2)
      call ice_HaloUpdate (fld2,               halo_info, &
                           field_loc_NEcorner, field_type_vector)
      call unstack_velocity_field(fld2, uvel, vvel)
      call ice_timer_stop(timer_bound)

      if (maskhalo_dyn) then
         call ice_timer_start(timer_bound)
         halomask = 0
         where (iceumask) halomask = 1
         call ice_HaloUpdate (halomask,          halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_timer_stop(timer_bound)
         call ice_HaloMask(halo_info_mask, halo_info, halomask)
      endif

      !-----------------------------------------------------------------
      ! seabed stress factor Tbu (Tbu is part of Cb coefficient)  
      !-----------------------------------------------------------------
      
      if (seabed_stress) then

       !$OMP PARALLEL DO PRIVATE(iblk)
       do iblk = 1, nblocks
          
          if ( seabed_stress_method == 'LKD' ) then
             
             call seabed_stress_factor_LKD (nx_block,         ny_block,       &
                                            icellu  (iblk),                   &
                                            indxui(:,iblk),   indxuj(:,iblk), &
                                            vice(:,:,iblk),   aice(:,:,iblk), &
                                            hwater(:,:,iblk), Tbu(:,:,iblk))

          elseif ( seabed_stress_method == 'probabilistic' ) then
             
             call seabed_stress_factor_prob (nx_block,         ny_block,                   &
                                             icellt(iblk), indxti(:,iblk), indxtj(:,iblk), &
                                             icellu(iblk), indxui(:,iblk), indxuj(:,iblk), &
                                             aicen(:,:,:,iblk), vicen(:,:,:,iblk),         &
                                             hwater(:,:,iblk), Tbu(:,:,iblk))
          endif

       enddo
       !$OMP END PARALLEL DO 
      endif
      
      do ksub = 1,ndte        ! subcycling

      !-----------------------------------------------------------------
      ! stress tensor equation, total surface stress
      !-----------------------------------------------------------------

         !$TCXOMP PARALLEL DO PRIVATE(iblk,strtmp)
         do iblk = 1, nblocks

!      call ice_timer_start(timer_tmp1) ! dynamics
            call stress_eap  (nx_block,             ny_block,             &
                              ksub,                 ndte,                 &
                              icellt(iblk),                               &
                              indxti      (:,iblk), indxtj      (:,iblk), &
                              arlx1i,               denom1,         &
                              uvel      (:,:,iblk), vvel      (:,:,iblk), &
                              dxt       (:,:,iblk), dyt       (:,:,iblk), &
                              dxhy      (:,:,iblk), dyhx      (:,:,iblk), &
                              cxp       (:,:,iblk), cyp       (:,:,iblk), &
                              cxm       (:,:,iblk), cym       (:,:,iblk), &
                              tarear    (:,:,iblk), strength  (:,:,iblk), &
                              a11_1     (:,:,iblk), a11_2   (:,:,iblk),   &
                              a11_3     (:,:,iblk), a11_4   (:,:,iblk),   &
                              a12_1     (:,:,iblk), a12_2   (:,:,iblk),   &
                              a12_3     (:,:,iblk), a12_4   (:,:,iblk),   &
                              stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), &
                              stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), &
                              stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), &
                              stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), &
                              stress12_1(:,:,iblk), stress12_2(:,:,iblk), &
                              stress12_3(:,:,iblk), stress12_4(:,:,iblk), &
                              shear     (:,:,iblk), divu      (:,:,iblk), &
                              e11       (:,:,iblk), e12       (:,:,iblk), &
                              e22       (:,:,iblk),                       &
                              s11       (:,:,iblk), s12       (:,:,iblk), &
                              s22       (:,:,iblk),                       &
                              yieldstress11 (:,:,iblk),                   &
                              yieldstress12 (:,:,iblk),                   &
                              yieldstress22 (:,:,iblk),                   &
!                             rdg_conv  (:,:,iblk), rdg_shear (:,:,iblk), &
                              rdg_conv  (:,:,iblk), &
                              strtmp    (:,:,:))
!      call ice_timer_stop(timer_tmp1) ! dynamics

      !-----------------------------------------------------------------
      ! momentum equation
      !-----------------------------------------------------------------

            call stepu (nx_block,            ny_block,           & 
                        icellu       (iblk), Cdn_ocn (:,:,iblk), & 
                        indxui     (:,iblk), indxuj    (:,iblk), & 
                        ksub,                                    &
                        aiu      (:,:,iblk), strtmp  (:,:,:),    & 
                        uocn     (:,:,iblk), vocn    (:,:,iblk), &     
                        waterx   (:,:,iblk), watery  (:,:,iblk), & 
                        forcex   (:,:,iblk), forcey  (:,:,iblk), & 
                        umassdti (:,:,iblk), fm      (:,:,iblk), & 
                        uarear   (:,:,iblk),                     & 
                        strintx  (:,:,iblk), strinty (:,:,iblk), &
                        taubx    (:,:,iblk), tauby   (:,:,iblk), & 
                        uvel_init(:,:,iblk), vvel_init(:,:,iblk),&
                        uvel     (:,:,iblk), vvel    (:,:,iblk), &
                        Tbu      (:,:,iblk))

      !-----------------------------------------------------------------
      ! evolution of structure tensor A
      !-----------------------------------------------------------------

!      call ice_timer_start(timer_tmp3) ! dynamics
            if (mod(ksub,10) == 1) then ! only called every 10th timestep
            call stepa (nx_block,          ny_block,                &
                        dtei,              icellt     (iblk),       &
                        indxti   (:,iblk), indxtj    (:,iblk),      &
                        a11    (:,:,iblk), a12  (:,:,iblk),         &
                        a11_1  (:,:,iblk), a11_2   (:,:,iblk),      &
                        a11_3  (:,:,iblk), a11_4   (:,:,iblk),      &
                        a12_1  (:,:,iblk), a12_2   (:,:,iblk),      &
                        a12_3  (:,:,iblk), a12_4   (:,:,iblk),      &
                        stressp_1(:,:,iblk), stressp_2(:,:,iblk),   &
                        stressp_3(:,:,iblk), stressp_4(:,:,iblk),   &
                        stressm_1(:,:,iblk), stressm_2(:,:,iblk),   &
                        stressm_3(:,:,iblk), stressm_4(:,:,iblk),   &
                        stress12_1(:,:,iblk), stress12_2(:,:,iblk), &
                        stress12_3(:,:,iblk), stress12_4(:,:,iblk))
            endif
!      call ice_timer_stop(timer_tmp3) ! dynamics
         enddo
         !$TCXOMP END PARALLEL DO

         call stack_velocity_field(uvel, vvel, fld2)
         call ice_timer_start(timer_bound)
         if (maskhalo_dyn) then
            call ice_HaloUpdate (fld2,               halo_info_mask, &
                                 field_loc_NEcorner, field_type_vector)
         else
            call ice_HaloUpdate (fld2,               halo_info, &
                                 field_loc_NEcorner, field_type_vector)
         endif
         call ice_timer_stop(timer_bound)
         call unstack_velocity_field(fld2, uvel, vvel)

      enddo                     ! subcycling

      deallocate(fld2)
      if (maskhalo_dyn) call ice_HaloDestroy(halo_info_mask)

      !-----------------------------------------------------------------
      ! ice-ocean stress
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks

         call dyn_finish                               & 
              (nx_block,           ny_block,           & 
               icellu      (iblk), Cdn_ocn (:,:,iblk), & 
               indxui    (:,iblk), indxuj    (:,iblk), & 
               uvel    (:,:,iblk), vvel    (:,:,iblk), & 
               uocn    (:,:,iblk), vocn    (:,:,iblk), & 
               aiu     (:,:,iblk), fm      (:,:,iblk), &
               strintx (:,:,iblk), strinty (:,:,iblk), &
               strairx (:,:,iblk), strairy (:,:,iblk), & 
               strocnx (:,:,iblk), strocny (:,:,iblk), & 
               strocnxT(:,:,iblk), strocnyT(:,:,iblk))

      enddo
      !$OMP END PARALLEL DO

      call u2tgrid_vector(strocnxT)    ! shift
      call u2tgrid_vector(strocnyT)

      call ice_timer_stop(timer_dynamics)    ! dynamics

      end subroutine eap

!=======================================================================

! Initialize parameters and variables needed for the eap dynamics
! (based on init_dyn)

      subroutine init_eap

      use ice_blocks, only: nx_block, ny_block
      use ice_domain, only: nblocks

      ! local variables

      integer (kind=int_kind) :: &
         i, j, &
         iblk          ! block index

      real (kind=dbl_kind), parameter :: & 
         eps6 = 1.0e-6_dbl_kind

      integer (kind=int_kind) :: & 
         ix, iy, iz, ia

      integer (kind=int_kind), parameter :: & 
         nz = 100

      real (kind=dbl_kind) :: & 
         ainit, xinit, yinit, zinit, &
         da, dx, dy, dz, &
         pi, pih, piq, phi

      character(len=*), parameter :: subname = '(init_eap)'

      call icepack_query_parameters(pi_out=pi, pih_out=pih, piq_out=piq)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)
      phi = pi/c12 ! diamond shaped floe smaller angle (default phi = 30 deg)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
      do j = 1, ny_block
      do i = 1, nx_block
         e11(i,j,iblk) = c0
         e12(i,j,iblk) = c0
         e22(i,j,iblk) = c0
         s11(i,j,iblk) = c0
         s12(i,j,iblk) = c0
         s22(i,j,iblk) = c0
         yieldstress11(i,j,iblk) = c0
         yieldstress12(i,j,iblk) = c0
         yieldstress22(i,j,iblk) = c0
         a11_1 (i,j,iblk) = p5
         a11_2 (i,j,iblk) = p5
         a11_3 (i,j,iblk) = p5
         a11_4 (i,j,iblk) = p5
         a12_1 (i,j,iblk) = c0
         a12_2 (i,j,iblk) = c0
         a12_3 (i,j,iblk) = c0
         a12_4 (i,j,iblk) = c0
      enddo                     ! i
      enddo                     ! j
      enddo                     ! iblk
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! create lookup table for eap dynamics (see Appendix A1)
      !-----------------------------------------------------------------

      da = p5/real(na_yield-1,kind=dbl_kind)
      ainit = p5 - da
      dx = pi/real(nx_yield-1,kind=dbl_kind)
      xinit = pi + piq - dx
      dz = pi/real(nz,kind=dbl_kind)
      zinit = -pih
      dy = pi/real(ny_yield-1,kind=dbl_kind)
      yinit = -dy

      do ia=1,na_yield
       do ix=1,nx_yield
        do iy=1,ny_yield
         s11r(ix,iy,ia) = c0
         s12r(ix,iy,ia) = c0
         s22r(ix,iy,ia) = c0
         s11s(ix,iy,ia) = c0
         s12s(ix,iy,ia) = c0
         s22s(ix,iy,ia) = c0
         if (ia <= na_yield-1) then
          do iz=1,nz
          s11r(ix,iy,ia) = s11r(ix,iy,ia) + 1*w1(ainit+ia*da)* &
           exp(-w2(ainit+ia*da)*(zinit+iz*dz)*(zinit+iz*dz))* &
           s11kr(xinit+ix*dx,yinit+iy*dy,zinit+iz*dz,phi)*dz/sin(c2*phi)
          s12r(ix,iy,ia) = s12r(ix,iy,ia) + 1*w1(ainit+ia*da)* &
           exp(-w2(ainit+ia*da)*(zinit+iz*dz)*(zinit+iz*dz))* &
           s12kr(xinit+ix*dx,yinit+iy*dy,zinit+iz*dz,phi)*dz/sin(c2*phi)
          s22r(ix,iy,ia) = s22r(ix,iy,ia) + 1*w1(ainit+ia*da)* &
           exp(-w2(ainit+ia*da)*(zinit+iz*dz)*(zinit+iz*dz))* &
           s22kr(xinit+ix*dx,yinit+iy*dy,zinit+iz*dz,phi)*dz/sin(c2*phi) 
          s11s(ix,iy,ia) = s11s(ix,iy,ia) + 1*w1(ainit+ia*da)* &
           exp(-w2(ainit+ia*da)*(zinit+iz*dz)*(zinit+iz*dz))* &
           s11ks(xinit+ix*dx,yinit+iy*dy,zinit+iz*dz,phi)*dz/sin(c2*phi)
          s12s(ix,iy,ia) = s12s(ix,iy,ia) + 1*w1(ainit+ia*da)* &
           exp(-w2(ainit+ia*da)*(zinit+iz*dz)*(zinit+iz*dz))* &
           s12ks(xinit+ix*dx,yinit+iy*dy,zinit+iz*dz,phi)*dz/sin(c2*phi)
          s22s(ix,iy,ia) = s22s(ix,iy,ia) + 1*w1(ainit+ia*da)* &
           exp(-w2(ainit+ia*da)*(zinit+iz*dz)*(zinit+iz*dz))* &
           s22ks(xinit+ix*dx,yinit+iy*dy,zinit+iz*dz,phi)*dz/sin(c2*phi)
          enddo
          if (abs(s11r(ix,iy,ia)) < eps6) s11r(ix,iy,ia) = c0
          if (abs(s12r(ix,iy,ia)) < eps6) s12r(ix,iy,ia) = c0
          if (abs(s22r(ix,iy,ia)) < eps6) s22r(ix,iy,ia) = c0
          if (abs(s11s(ix,iy,ia)) < eps6) s11s(ix,iy,ia) = c0
          if (abs(s12s(ix,iy,ia)) < eps6) s12s(ix,iy,ia) = c0
          if (abs(s22s(ix,iy,ia)) < eps6) s22s(ix,iy,ia) = c0
         else
          s11r(ix,iy,ia) = p5*s11kr(xinit+ix*dx,yinit+iy*dy,c0,phi)/sin(c2*phi)
          s12r(ix,iy,ia) = p5*s12kr(xinit+ix*dx,yinit+iy*dy,c0,phi)/sin(c2*phi)
          s22r(ix,iy,ia) = p5*s22kr(xinit+ix*dx,yinit+iy*dy,c0,phi)/sin(c2*phi)
          s11s(ix,iy,ia) = p5*s11ks(xinit+ix*dx,yinit+iy*dy,c0,phi)/sin(c2*phi)
          s12s(ix,iy,ia) = p5*s12ks(xinit+ix*dx,yinit+iy*dy,c0,phi)/sin(c2*phi)
          s22s(ix,iy,ia) = p5*s22ks(xinit+ix*dx,yinit+iy*dy,c0,phi)/sin(c2*phi)
          if (abs(s11r(ix,iy,ia)) < eps6) s11r(ix,iy,ia) = c0
          if (abs(s12r(ix,iy,ia)) < eps6) s12r(ix,iy,ia) = c0
          if (abs(s22r(ix,iy,ia)) < eps6) s22r(ix,iy,ia) = c0
          if (abs(s11s(ix,iy,ia)) < eps6) s11s(ix,iy,ia) = c0
          if (abs(s12s(ix,iy,ia)) < eps6) s12s(ix,iy,ia) = c0
          if (abs(s22s(ix,iy,ia)) < eps6) s22s(ix,iy,ia) = c0
         endif
        enddo
       enddo
      enddo

      end subroutine init_eap

!=======================================================================
! Function : w1 (see Gaussian function psi in Tsamados et al 2013)

      FUNCTION w1(a)
      double precision, intent(in) :: a

      real (kind=dbl_kind) :: w1
      character(len=*), parameter :: subname = '(w1)'

      w1 = - 223.87569446_dbl_kind &
           + 2361.2198663_dbl_kind*a &
           - 10606.56079975_dbl_kind*a*a &
           + 26315.50025642_dbl_kind*a*a*a &
           - 38948.30444297_dbl_kind*a*a*a*a &
           + 34397.72407466_dbl_kind*a*a*a*a*a &
           - 16789.98003081_dbl_kind*a*a*a*a*a*a &
           + 3495.82839237_dbl_kind*a*a*a*a*a*a*a

      end FUNCTION w1

!=======================================================================
! Function : w2 (see Gaussian function psi in Tsamados et al 2013)

      FUNCTION w2(a)
      double precision, intent(in) :: a

      real (kind=dbl_kind) :: w2
      character(len=*), parameter :: subname = '(w2)'

      w2 = - 6670.68911883_dbl_kind &
           + 70222.33061536_dbl_kind*a &
           - 314871.71525448_dbl_kind*a*a &
           + 779570.02793492_dbl_kind*a*a*a &
           - 1151098.82436864_dbl_kind*a*a*a*a &
           + 1013896.59464498_dbl_kind*a*a*a*a*a &
           - 493379.44906738_dbl_kind*a*a*a*a*a*a &
           + 102356.551518_dbl_kind*a*a*a*a*a*a*a

      end FUNCTION w2

!=======================================================================
! Function : s11kr

      FUNCTION s11kr(x,y,z,phi) 

      real (kind=dbl_kind), intent(in) :: &
        x,y,z,phi

      real (kind=dbl_kind) :: &
      s11kr, p

      real (kind=dbl_kind) :: &
      n1t2i11, n1t2i12, n1t2i21, n1t2i22, &
      n2t1i11, n2t1i12, n2t1i21, n2t1i22, &
!     t1t2i11, t1t2i12, t1t2i21, t1t2i22, &
!     t2t1i11, t2t1i12, t2t1i21, t2t1i22, &
      d11, d12, d22, &
      IIn1t2, IIn2t1, &
!     IIt1t2, &
      Hen1t2, Hen2t1, &
      pih, puny
      character(len=*), parameter :: subname = '(s11kr)'

      call icepack_query_parameters(pih_out=pih, puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      p = phi

      n1t2i11 = cos(z+pih-p) * cos(z+p)
      n1t2i12 = cos(z+pih-p) * sin(z+p)
      n1t2i21 = sin(z+pih-p) * cos(z+p)
      n1t2i22 = sin(z+pih-p) * sin(z+p)
      n2t1i11 = cos(z-pih+p) * cos(z-p)
      n2t1i12 = cos(z-pih+p) * sin(z-p)
      n2t1i21 = sin(z-pih+p) * cos(z-p)
      n2t1i22 = sin(z-pih+p) * sin(z-p)
!     t1t2i11 = cos(z-p) * cos(z+p)
!     t1t2i12 = cos(z-p) * sin(z+p)
!     t1t2i21 = sin(z-p) * cos(z+p)
!     t1t2i22 = sin(z-p) * sin(z+p)
!     t2t1i11 = cos(z+p) * cos(z-p)
!     t2t1i12 = cos(z+p) * sin(z-p)
!     t2t1i21 = sin(z+p) * cos(z-p)
!     t2t1i22 = sin(z+p) * sin(z-p)
! In expression of tensor d, with this formulatin d(x)=-d(x+pi)
! Solution, when diagonalizing always check sgn(a11-a22) if > then keep x else x=x-pi/2
      d11 = cos(y)*cos(y)*(cos(x)+sin(x)*tan(y)*tan(y))
      d12 = cos(y)*cos(y)*tan(y)*(-cos(x)+sin(x))
      d22 = cos(y)*cos(y)*(sin(x)+cos(x)*tan(y)*tan(y))
      IIn1t2 = n1t2i11 * d11 + (n1t2i12 + n1t2i21) * d12 + n1t2i22 * d22
      IIn2t1 = n2t1i11 * d11 + (n2t1i12 + n2t1i21) * d12 + n2t1i22 * d22
!     IIt1t2 = t1t2i11 * d11 + (t1t2i12 + t1t2i21) * d12 + t1t2i22 * d22

      if (-IIn1t2>=puny) then
      Hen1t2 = c1
      else
      Hen1t2 = c0
      endif

      if (-IIn2t1>=puny) then
      Hen2t1 = c1
      else
      Hen2t1 = c0
      endif

      s11kr = (- Hen1t2 * n1t2i11 - Hen2t1 * n2t1i11)

      end FUNCTION s11kr

!=======================================================================
! Function : s12kr

      FUNCTION s12kr(x,y,z,phi)

      real (kind=dbl_kind), intent(in) :: &
        x,y,z,phi

      real (kind=dbl_kind) :: &
      s12kr, s12r0, s21r0, p

      real (kind=dbl_kind) :: &
      n1t2i11, n1t2i12, n1t2i21, n1t2i22, &
      n2t1i11, n2t1i12, n2t1i21, n2t1i22, &
!     t1t2i11, t1t2i12, t1t2i21, t1t2i22, &
!     t2t1i11, t2t1i12, t2t1i21, t2t1i22, &
      d11, d12, d22, &
      IIn1t2, IIn2t1, &
!     IIt1t2, &
      Hen1t2, Hen2t1, &
      pih, puny
      character(len=*), parameter :: subname = '(s12kr)'

      call icepack_query_parameters(pih_out=pih, puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      p = phi

      n1t2i11 = cos(z+pih-p) * cos(z+p)
      n1t2i12 = cos(z+pih-p) * sin(z+p)
      n1t2i21 = sin(z+pih-p) * cos(z+p)
      n1t2i22 = sin(z+pih-p) * sin(z+p)
      n2t1i11 = cos(z-pih+p) * cos(z-p)
      n2t1i12 = cos(z-pih+p) * sin(z-p)
      n2t1i21 = sin(z-pih+p) * cos(z-p)
      n2t1i22 = sin(z-pih+p) * sin(z-p)
!     t1t2i11 = cos(z-p) * cos(z+p)
!     t1t2i12 = cos(z-p) * sin(z+p)
!     t1t2i21 = sin(z-p) * cos(z+p)
!     t1t2i22 = sin(z-p) * sin(z+p)
!     t2t1i11 = cos(z+p) * cos(z-p)
!     t2t1i12 = cos(z+p) * sin(z-p)
!     t2t1i21 = sin(z+p) * cos(z-p)
!     t2t1i22 = sin(z+p) * sin(z-p)
      d11 = cos(y)*cos(y)*(cos(x)+sin(x)*tan(y)*tan(y))
      d12 = cos(y)*cos(y)*tan(y)*(-cos(x)+sin(x))
      d22 = cos(y)*cos(y)*(sin(x)+cos(x)*tan(y)*tan(y))
      IIn1t2 = n1t2i11 * d11 + (n1t2i12 + n1t2i21) * d12 + n1t2i22 * d22
      IIn2t1 = n2t1i11 * d11 + (n2t1i12 + n2t1i21) * d12 + n2t1i22 * d22
!     IIt1t2 = t1t2i11 * d11 + (t1t2i12 + t1t2i21) * d12 + t1t2i22 * d22

      if (-IIn1t2>=puny) then
      Hen1t2 = c1
      else
      Hen1t2 = c0
      endif

      if (-IIn2t1>=puny) then
      Hen2t1 = c1
      else
      Hen2t1 = c0
      endif

      s12r0 = (- Hen1t2 * n1t2i12 - Hen2t1 * n2t1i12)
      s21r0 = (- Hen1t2 * n1t2i21 - Hen2t1 * n2t1i21)
      s12kr=p5*(s12r0+s21r0)

      end FUNCTION s12kr

!=======================================================================
! Function : s22r

      FUNCTION s22kr(x,y,z,phi)

      real (kind=dbl_kind), intent(in) :: &
        x,y,z,phi

      real (kind=dbl_kind) :: &
      s22kr, p

      real (kind=dbl_kind) :: &
      n1t2i11, n1t2i12, n1t2i21, n1t2i22, &
      n2t1i11, n2t1i12, n2t1i21, n2t1i22, &
!     t1t2i11, t1t2i12, t1t2i21, t1t2i22, &
!     t2t1i11, t2t1i12, t2t1i21, t2t1i22, &
      d11, d12, d22, &
      IIn1t2, IIn2t1, &
!     IIt1t2, &
      Hen1t2, Hen2t1, &
      pih, puny
      character(len=*), parameter :: subname = '(s22kr)'

      call icepack_query_parameters(pih_out=pih, puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      p = phi

      n1t2i11 = cos(z+pih-p) * cos(z+p)
      n1t2i12 = cos(z+pih-p) * sin(z+p)
      n1t2i21 = sin(z+pih-p) * cos(z+p)
      n1t2i22 = sin(z+pih-p) * sin(z+p)
      n2t1i11 = cos(z-pih+p) * cos(z-p)
      n2t1i12 = cos(z-pih+p) * sin(z-p)
      n2t1i21 = sin(z-pih+p) * cos(z-p)
      n2t1i22 = sin(z-pih+p) * sin(z-p)
!     t1t2i11 = cos(z-p) * cos(z+p)
!     t1t2i12 = cos(z-p) * sin(z+p)
!     t1t2i21 = sin(z-p) * cos(z+p)
!     t1t2i22 = sin(z-p) * sin(z+p)
!     t2t1i11 = cos(z+p) * cos(z-p)
!     t2t1i12 = cos(z+p) * sin(z-p)
!     t2t1i21 = sin(z+p) * cos(z-p)
!     t2t1i22 = sin(z+p) * sin(z-p)
      d11 = cos(y)*cos(y)*(cos(x)+sin(x)*tan(y)*tan(y))
      d12 = cos(y)*cos(y)*tan(y)*(-cos(x)+sin(x))
      d22 = cos(y)*cos(y)*(sin(x)+cos(x)*tan(y)*tan(y))
      IIn1t2 = n1t2i11 * d11 + (n1t2i12 + n1t2i21) * d12 + n1t2i22 * d22
      IIn2t1 = n2t1i11 * d11 + (n2t1i12 + n2t1i21) * d12 + n2t1i22 * d22
!     IIt1t2 = t1t2i11 * d11 + (t1t2i12 + t1t2i21) * d12 + t1t2i22 * d22

      if (-IIn1t2>=puny) then
      Hen1t2 = c1
      else
      Hen1t2 = c0
      endif

      if (-IIn2t1>=puny) then
      Hen2t1 = c1
      else
      Hen2t1 = c0
      endif

      s22kr = (- Hen1t2 * n1t2i22 - Hen2t1 * n2t1i22)

      end FUNCTION s22kr

!=======================================================================
! Function : s11ks

      FUNCTION s11ks(x,y,z,phi)

      real (kind=dbl_kind), intent(in):: &
        x,y,z,phi

      real (kind=dbl_kind) :: &
      s11ks, p

      real (kind=dbl_kind) :: &
      n1t2i11, n1t2i12, n1t2i21, n1t2i22, &
      n2t1i11, n2t1i12, n2t1i21, n2t1i22, &
      t1t2i11, &
      t1t2i12, t1t2i21, t1t2i22, &
      t2t1i11, &
!     t2t1i12, t2t1i21, t2t1i22, &
      d11, d12, d22, &
      IIn1t2, IIn2t1, IIt1t2, &
      Hen1t2, Hen2t1, &
      pih, puny
      character(len=*), parameter :: subname = '(s11ks)'

      call icepack_query_parameters(pih_out=pih, puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      p = phi

      n1t2i11 = cos(z+pih-p) * cos(z+p)
      n1t2i12 = cos(z+pih-p) * sin(z+p)
      n1t2i21 = sin(z+pih-p) * cos(z+p)
      n1t2i22 = sin(z+pih-p) * sin(z+p)
      n2t1i11 = cos(z-pih+p) * cos(z-p)
      n2t1i12 = cos(z-pih+p) * sin(z-p)
      n2t1i21 = sin(z-pih+p) * cos(z-p)
      n2t1i22 = sin(z-pih+p) * sin(z-p)
      t1t2i11 = cos(z-p) * cos(z+p)
      t1t2i12 = cos(z-p) * sin(z+p)
      t1t2i21 = sin(z-p) * cos(z+p)
      t1t2i22 = sin(z-p) * sin(z+p)
      t2t1i11 = cos(z+p) * cos(z-p)
!     t2t1i12 = cos(z+p) * sin(z-p)
!     t2t1i21 = sin(z+p) * cos(z-p)
!     t2t1i22 = sin(z+p) * sin(z-p)
      d11 = cos(y)*cos(y)*(cos(x)+sin(x)*tan(y)*tan(y))
      d12 = cos(y)*cos(y)*tan(y)*(-cos(x)+sin(x))
      d22 = cos(y)*cos(y)*(sin(x)+cos(x)*tan(y)*tan(y))
      IIn1t2 = n1t2i11 * d11 + (n1t2i12 + n1t2i21) * d12 + n1t2i22 * d22
      IIn2t1 = n2t1i11 * d11 + (n2t1i12 + n2t1i21) * d12 + n2t1i22 * d22
      IIt1t2 = t1t2i11 * d11 + (t1t2i12 + t1t2i21) * d12 + t1t2i22 * d22

      if (-IIn1t2>=puny) then
      Hen1t2 = c1
      else
      Hen1t2 = c0
      endif

      if (-IIn2t1>=puny) then
      Hen2t1 = c1
      else
      Hen2t1 = c0
      endif

      s11ks = sign(c1,IIt1t2+puny)*(Hen1t2 * t1t2i11 + Hen2t1 * t2t1i11)

      end FUNCTION s11ks

!=======================================================================
! Function : s12ks

      FUNCTION s12ks(x,y,z,phi)

      real (kind=dbl_kind), intent(in) :: &
        x,y,z,phi

      real (kind=dbl_kind) :: &
      s12ks,s12s0,s21s0,p

      real (kind=dbl_kind) :: &
      n1t2i11, n1t2i12, n1t2i21, n1t2i22, &
      n2t1i11, n2t1i12, n2t1i21, n2t1i22, &
      t1t2i11, t1t2i12, t1t2i21, t1t2i22, &
!     t2t1i11, t2t1i22, &
      t2t1i12, t2t1i21, &
      d11, d12, d22, &
      IIn1t2, IIn2t1, IIt1t2, &
      Hen1t2, Hen2t1, &
      pih, puny
      character(len=*), parameter :: subname = '(s12ks)'

      call icepack_query_parameters(pih_out=pih, puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      p =phi

      n1t2i11 = cos(z+pih-p) * cos(z+p)
      n1t2i12 = cos(z+pih-p) * sin(z+p)
      n1t2i21 = sin(z+pih-p) * cos(z+p)
      n1t2i22 = sin(z+pih-p) * sin(z+p)
      n2t1i11 = cos(z-pih+p) * cos(z-p)
      n2t1i12 = cos(z-pih+p) * sin(z-p)
      n2t1i21 = sin(z-pih+p) * cos(z-p)
      n2t1i22 = sin(z-pih+p) * sin(z-p)
      t1t2i11 = cos(z-p) * cos(z+p)
      t1t2i12 = cos(z-p) * sin(z+p)
      t1t2i21 = sin(z-p) * cos(z+p)
      t1t2i22 = sin(z-p) * sin(z+p)
!     t2t1i11 = cos(z+p) * cos(z-p)
      t2t1i12 = cos(z+p) * sin(z-p)
      t2t1i21 = sin(z+p) * cos(z-p)
!     t2t1i22 = sin(z+p) * sin(z-p)
      d11 = cos(y)*cos(y)*(cos(x)+sin(x)*tan(y)*tan(y))
      d12 = cos(y)*cos(y)*tan(y)*(-cos(x)+sin(x))
      d22 = cos(y)*cos(y)*(sin(x)+cos(x)*tan(y)*tan(y))
      IIn1t2 = n1t2i11 * d11 + (n1t2i12 + n1t2i21) * d12 + n1t2i22 * d22
      IIn2t1 = n2t1i11 * d11 + (n2t1i12 + n2t1i21) * d12 + n2t1i22 * d22
      IIt1t2 = t1t2i11 * d11 + (t1t2i12 + t1t2i21) * d12 + t1t2i22 * d22

      if (-IIn1t2>=puny) then
      Hen1t2 = c1
      else
      Hen1t2 = c0
      endif

      if (-IIn2t1>=puny) then
      Hen2t1 = c1
      else
      Hen2t1 = c0
      endif

      s12s0 = sign(c1,IIt1t2+puny)*(Hen1t2 * t1t2i12 + Hen2t1 * t2t1i12)
      s21s0 = sign(c1,IIt1t2+puny)*(Hen1t2 * t1t2i21 + Hen2t1 * t2t1i21)
      s12ks=p5*(s12s0+s21s0)

      end FUNCTION s12ks

!=======================================================================
! Function : s22ks

      FUNCTION s22ks(x,y,z,phi) 

      real (kind=dbl_kind), intent(in) :: &
        x,y,z,phi

      real (kind=dbl_kind) :: &
      s22ks,p

      real (kind=dbl_kind) :: &
      n1t2i11, n1t2i12, n1t2i21, n1t2i22, &
      n2t1i11, n2t1i12, n2t1i21, n2t1i22, &
      t1t2i11, t1t2i12, t1t2i21, t1t2i22, &
!     t2t1i11, t2t1i12, t2t1i21, &
      t2t1i22, &
      d11, d12, d22, &
      IIn1t2, IIn2t1, IIt1t2, &
      Hen1t2, Hen2t1, &
      pih, puny
      character(len=*), parameter :: subname = '(s22ks)'

      call icepack_query_parameters(pih_out=pih, puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      p = phi

      n1t2i11 = cos(z+pih-p) * cos(z+p)
      n1t2i12 = cos(z+pih-p) * sin(z+p)
      n1t2i21 = sin(z+pih-p) * cos(z+p)
      n1t2i22 = sin(z+pih-p) * sin(z+p)
      n2t1i11 = cos(z-pih+p) * cos(z-p)
      n2t1i12 = cos(z-pih+p) * sin(z-p)
      n2t1i21 = sin(z-pih+p) * cos(z-p)
      n2t1i22 = sin(z-pih+p) * sin(z-p)
      t1t2i11 = cos(z-p) * cos(z+p)
      t1t2i12 = cos(z-p) * sin(z+p)
      t1t2i21 = sin(z-p) * cos(z+p)
      t1t2i22 = sin(z-p) * sin(z+p)
!     t2t1i11 = cos(z+p) * cos(z-p)
!     t2t1i12 = cos(z+p) * sin(z-p)
!     t2t1i21 = sin(z+p) * cos(z-p)
      t2t1i22 = sin(z+p) * sin(z-p)
      d11 = cos(y)*cos(y)*(cos(x)+sin(x)*tan(y)*tan(y))
      d12 = cos(y)*cos(y)*tan(y)*(-cos(x)+sin(x))
      d22 = cos(y)*cos(y)*(sin(x)+cos(x)*tan(y)*tan(y))
      IIn1t2 = n1t2i11 * d11 + (n1t2i12 + n1t2i21) * d12 + n1t2i22 * d22
      IIn2t1 = n2t1i11 * d11 + (n2t1i12 + n2t1i21) * d12 + n2t1i22 * d22
      IIt1t2 = t1t2i11 * d11 + (t1t2i12 + t1t2i21) * d12 + t1t2i22 * d22

      if (-IIn1t2>=puny) then
      Hen1t2 = c1
      else
      Hen1t2 = c0
      endif

      if (-IIn2t1>=puny) then
      Hen2t1 = c1
      else
      Hen2t1 = c0
      endif

      s22ks = sign(c1,IIt1t2+puny)*(Hen1t2 * t1t2i22 + Hen2t1 * t2t1i22)

      end FUNCTION s22ks


!=======================================================================

! Computes the rates of strain and internal stress components for
! each of the four corners on each T-grid cell.
! Computes stress terms for the momentum equation
! (based on subroutine stress)

      subroutine stress_eap  (nx_block,   ny_block,       &
                              ksub,       ndte,           &
                              icellt,                     &
                              indxti,     indxtj,         &
                              arlx1i,     denom1,         &
                              uvel,       vvel,           &
                              dxt,        dyt,            &
                              dxhy,       dyhx,           &
                              cxp,        cyp,            &
                              cxm,        cym,            &
                              tarear,     strength,       &
                              a11_1, a11_2, a11_3, a11_4, &
                              a12_1, a12_2, a12_3, a12_4, &
                              stressp_1,  stressp_2,      &
                              stressp_3,  stressp_4,      &
                              stressm_1,  stressm_2,      &
                              stressm_3,  stressm_4,      &
                              stress12_1, stress12_2,     &
                              stress12_3, stress12_4,     &
                              shear,      divu,           &
                              e11,        e12,            &
                              e22,                        &
                              s11,        s12,            &
                              s22,                        &
                              yieldstress11,              &
                              yieldstress12,              &
                              yieldstress22,              &
!                             rdg_conv,   rdg_shear,      &
                              rdg_conv, &
                              strtmp)

!echmod tmp
!      use ice_timers, only:  &
!          ice_timer_start, ice_timer_stop, &
!          timer_tmp1, timer_tmp2, timer_tmp3

      integer (kind=int_kind), intent(in) :: & 
         nx_block, ny_block, & ! block dimensions
         ksub              , & ! subcycling step
         ndte              , & ! number of subcycles
         icellt                ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj       ! compressed index in j-direction

      real (kind=dbl_kind), intent(in) :: &
         arlx1i   , & ! dte/2T (original) or 1/alpha1 (revised)
         denom1       ! constant for stress equation

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         strength , & ! ice strength (N/m)
         uvel     , & ! x-component of velocity (m/s)
         vvel     , & ! y-component of velocity (m/s)
         dxt      , & ! width of T-cell through the middle (m)
         dyt      , & ! height of T-cell through the middle (m)
         dxhy     , & ! 0.5*(HTE - HTE)
         dyhx     , & ! 0.5*(HTN - HTN)
         cyp      , & ! 1.5*HTE - 0.5*HTE
         cxp      , & ! 1.5*HTN - 0.5*HTN
         cym      , & ! 0.5*HTE - 1.5*HTE
         cxm      , & ! 0.5*HTN - 1.5*HTN
         tarear       ! 1/tarea

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         stressp_1, stressp_2, stressp_3, stressp_4, & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4, & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4   ! sigma12

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         a11_1, a11_2, a11_3, a11_4, & ! structure tensor
         a12_1, a12_2, a12_3, a12_4              ! structure tensor

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         shear    , & ! strain rate II component (1/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         e11      , & ! components of strain rate tensor (1/s)
         e12      , & ! 
         e22      , & ! 
         s11      , & ! components of stress tensor (kg/s^2)
         s12      , & ! 
         s22      , & ! 
         yieldstress11, & ! components of yield stress tensor (kg/s^2)
         yieldstress12, &
         yieldstress22, &
         rdg_conv     ! convergence term for ridging (1/s)
!        rdg_shear    ! shear term for ridging (1/s)

      real (kind=dbl_kind), dimension(nx_block,ny_block,8), intent(out) :: &
         strtmp       ! stress combinations

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
         stressptmp_1, stressptmp_2, stressptmp_3, stressptmp_4, & ! sigma11+sigma22
         stressmtmp_1, stressmtmp_2, stressmtmp_3, stressmtmp_4, & ! sigma11-sigma22
         stress12tmp_1,stress12tmp_2,stress12tmp_3,stress12tmp_4   ! sigma12

      real (kind=dbl_kind) :: &
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        ssigpn, ssigps, ssigpe, ssigpw            , &
        ssigmn, ssigms, ssigme, ssigmw            , &
        ssig12n, ssig12s, ssig12e, ssig12w        , &
        ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122, &
        csigpne, csigpnw, csigpse, csigpsw        , &
        csigmne, csigmnw, csigmse, csigmsw        , &
        csig12ne, csig12nw, csig12se, csig12sw    , &
        str12ew, str12we, str12ns, str12sn        , &
        strp_tmp, strm_tmp, puny

      real (kind=dbl_kind) :: &
        alpharne, alpharnw, alpharsw, alpharse,     &
        alphasne, alphasnw, alphassw, alphasse

      character(len=*), parameter :: subname = '(stress_eap)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      call icepack_query_parameters(puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      strtmp(:,:,:) = c0

      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
         ! divergence  =  e_11 + e_22
         divune    = cyp(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  ) &
                   + cxp(i,j)*vvel(i  ,j  ) - dxt(i,j)*vvel(i  ,j-1)
         divunw    = cym(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  ) &
                   + cxp(i,j)*vvel(i-1,j  ) - dxt(i,j)*vvel(i-1,j-1)
         divusw    = cym(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1) &
                   + cxm(i,j)*vvel(i-1,j-1) + dxt(i,j)*vvel(i-1,j  )
         divuse    = cyp(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1) &
                   + cxm(i,j)*vvel(i  ,j-1) + dxt(i,j)*vvel(i  ,j  )

         ! tension strain rate  =  e_11 - e_22
         tensionne = -cym(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  ) &
                   +  cxm(i,j)*vvel(i  ,j  ) + dxt(i,j)*vvel(i  ,j-1)
         tensionnw = -cyp(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  ) &
                   +  cxm(i,j)*vvel(i-1,j  ) + dxt(i,j)*vvel(i-1,j-1)
         tensionsw = -cyp(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1) &
                   +  cxp(i,j)*vvel(i-1,j-1) - dxt(i,j)*vvel(i-1,j  )
         tensionse = -cym(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1) &
                   +  cxp(i,j)*vvel(i  ,j-1) - dxt(i,j)*vvel(i  ,j  )

         ! shearing strain rate  =  2*e_12
         shearne = -cym(i,j)*vvel(i  ,j  ) - dyt(i,j)*vvel(i-1,j  ) &
                 -  cxm(i,j)*uvel(i  ,j  ) - dxt(i,j)*uvel(i  ,j-1)
         shearnw = -cyp(i,j)*vvel(i-1,j  ) + dyt(i,j)*vvel(i  ,j  ) &
                 -  cxm(i,j)*uvel(i-1,j  ) - dxt(i,j)*uvel(i-1,j-1)
         shearsw = -cyp(i,j)*vvel(i-1,j-1) + dyt(i,j)*vvel(i  ,j-1) &
                 -  cxp(i,j)*uvel(i-1,j-1) + dxt(i,j)*uvel(i-1,j  )
         shearse = -cym(i,j)*vvel(i  ,j-1) - dyt(i,j)*vvel(i-1,j-1) &
                 -  cxp(i,j)*uvel(i  ,j-1) + dxt(i,j)*uvel(i  ,j  )

      !-----------------------------------------------------------------
      ! Stress updated depending on strain rate and structure tensor
      !-----------------------------------------------------------------
!      call ice_timer_start(timer_tmp2) ! dynamics

         ! ne
         call update_stress_rdg (ksub, ndte, divune, tensionne, &
                                 shearne, a11_1(i,j), a12_1(i,j), &
                                 stressptmp_1, stressmtmp_1, &
                                 stress12tmp_1, strength(i,j), &
                                 alpharne, alphasne)
         ! nw
         call update_stress_rdg (ksub, ndte, divunw, tensionnw, &
                                 shearnw, a11_2(i,j), a12_2(i,j), &
                                 stressptmp_2, stressmtmp_2, &
                                 stress12tmp_2, strength(i,j), &
                                 alpharnw, alphasnw)
         ! sw
         call update_stress_rdg (ksub, ndte, divusw, tensionsw, &
                                 shearsw, a11_3(i,j), a12_3(i,j), &
                                 stressptmp_3, stressmtmp_3, &
                                 stress12tmp_3, strength(i,j), &
                                 alpharsw, alphassw)
         ! se
         call update_stress_rdg (ksub, ndte, divuse, tensionse, &
                                 shearse, a11_4(i,j), a12_4(i,j), &
                                 stressptmp_4, stressmtmp_4, &
                                 stress12tmp_4, strength(i,j), &
                                 alpharse, alphasse)

!      call ice_timer_stop(timer_tmp2) ! dynamics
      !-----------------------------------------------------------------
      ! on last subcycle, save quantities for mechanical redistribution
      !-----------------------------------------------------------------
         if (ksub == ndte) then

            ! diagnostic only
            ! shear = sqrt(tension**2 + shearing**2)
            shear(i,j) = p25*tarear(i,j)*sqrt( &
                 (tensionne + tensionnw + tensionse + tensionsw)**2 &
                +  (shearne +   shearnw +   shearse +   shearsw)**2)

            divu(i,j) = p25*(divune + divunw + divuse + divusw) * tarear(i,j)
            rdg_conv(i,j)  = -min(p25*(alpharne + alpharnw &
                                     + alpharsw + alpharse),c0) * tarear(i,j)
            !rdg_shear=0 for computing closing_net in ridge_prep
            !rdg_shear(i,j) = p25*(alphasne + alphasnw &
            !                    + alphassw + alphasse) * tarear(i,j)
         endif

         e11(i,j) = p5*p25*(divune + divunw + divuse + divusw + &
                    tensionne + tensionnw + tensionse + tensionsw) * tarear(i,j)

         e12(i,j) = p5*p25*(shearne + shearnw + shearse + shearsw) * tarear(i,j)

         e22(i,j) = p5*p25*(divune + divunw + divuse + divusw - &
                    tensionne - tensionnw - tensionse - tensionsw) * tarear(i,j)

      !-----------------------------------------------------------------
      ! elastic relaxation, see Eq. A12-A14
      !-----------------------------------------------------------------

         stressp_1(i,j) = (stressp_1(i,j) + stressptmp_1*arlx1i) &
                          * denom1
         stressp_2(i,j) = (stressp_2(i,j) + stressptmp_2*arlx1i) &
                          * denom1
         stressp_3(i,j) = (stressp_3(i,j) + stressptmp_3*arlx1i) &
                          * denom1
         stressp_4(i,j) = (stressp_4(i,j) + stressptmp_4*arlx1i) &
                          * denom1

         stressm_1(i,j) = (stressm_1(i,j) + stressmtmp_1*arlx1i) &
                          * denom1
         stressm_2(i,j) = (stressm_2(i,j) + stressmtmp_2*arlx1i) &
                          * denom1
         stressm_3(i,j) = (stressm_3(i,j) + stressmtmp_3*arlx1i) &
                          * denom1
         stressm_4(i,j) = (stressm_4(i,j) + stressmtmp_4*arlx1i) &
                          * denom1

         stress12_1(i,j) = (stress12_1(i,j) + stress12tmp_1*arlx1i) &
                          * denom1
         stress12_2(i,j) = (stress12_2(i,j) + stress12tmp_2*arlx1i) &
                          * denom1
         stress12_3(i,j) = (stress12_3(i,j) + stress12tmp_3*arlx1i) &
                          * denom1
         stress12_4(i,j) = (stress12_4(i,j) + stress12tmp_4*arlx1i) &
                          * denom1

         s11(i,j) = p5 * p25 * (stressp_1(i,j) + stressp_2(i,j) &
                              + stressp_3(i,j) + stressp_4(i,j) &
                              + stressm_1(i,j) + stressm_2(i,j) &
                              + stressm_3(i,j) + stressm_4(i,j))
         s22(i,j) = p5 * p25 * (stressp_1(i,j) + stressp_2(i,j) &
                              + stressp_3(i,j) + stressp_4(i,j) &
                              - stressm_1(i,j) - stressm_2(i,j) &
                              - stressm_3(i,j) - stressm_4(i,j))
         s12(i,j) =      p25 * (stress12_1(i,j) + stress12_2(i,j) &
                              + stress12_3(i,j) + stress12_4(i,j))

         yieldstress11(i,j) = p5 * p25 * (stressptmp_1 + stressptmp_2 &
                                        + stressptmp_3 + stressptmp_4 &
                                        + stressmtmp_1 + stressmtmp_2 &
                                        + stressmtmp_3 + stressmtmp_4)
         yieldstress22(i,j) = p5 * p25 * (stressptmp_1 + stressptmp_2 &
                                        + stressptmp_3 + stressptmp_4 &
                                        - stressmtmp_1 - stressmtmp_2 &
                                        - stressmtmp_3 - stressmtmp_4)
         yieldstress12(i,j) =      p25 * (stress12tmp_1 + stress12tmp_2 &
                                        + stress12tmp_3 + stress12tmp_4)

      !-----------------------------------------------------------------
      ! Eliminate underflows.
      ! The following code is commented out because it is relatively
      ! expensive and most compilers include a flag that accomplishes
      ! the same thing more efficiently.  This code is cheaper than
      ! handling underflows if the compiler lacks a flag; uncomment
      ! it in that case.  The compiler flag is often described with the
      ! phrase "flush to zero".
      !-----------------------------------------------------------------

!      stressp_1(i,j) = sign(max(abs(stressp_1(i,j)),puny),stressp_1(i,j))
!      stressp_2(i,j) = sign(max(abs(stressp_2(i,j)),puny),stressp_2(i,j))
!      stressp_3(i,j) = sign(max(abs(stressp_3(i,j)),puny),stressp_3(i,j))
!      stressp_4(i,j) = sign(max(abs(stressp_4(i,j)),puny),stressp_4(i,j))

!      stressm_1(i,j) = sign(max(abs(stressm_1(i,j)),puny),stressm_1(i,j))
!      stressm_2(i,j) = sign(max(abs(stressm_2(i,j)),puny),stressm_2(i,j))
!      stressm_3(i,j) = sign(max(abs(stressm_3(i,j)),puny),stressm_3(i,j))
!      stressm_4(i,j) = sign(max(abs(stressm_4(i,j)),puny),stressm_4(i,j))

!      stress12_1(i,j) = sign(max(abs(stress12_1(i,j)),puny),stress12_1(i,j))
!      stress12_2(i,j) = sign(max(abs(stress12_2(i,j)),puny),stress12_2(i,j))
!      stress12_3(i,j) = sign(max(abs(stress12_3(i,j)),puny),stress12_3(i,j))
!      stress12_4(i,j) = sign(max(abs(stress12_4(i,j)),puny),stress12_4(i,j))

      !-----------------------------------------------------------------
      ! combinations of the stresses for the momentum equation ! kg/s^2
      !-----------------------------------------------------------------

         ssigpn  = stressp_1(i,j) + stressp_2(i,j)
         ssigps  = stressp_3(i,j) + stressp_4(i,j)
         ssigpe  = stressp_1(i,j) + stressp_4(i,j)
         ssigpw  = stressp_2(i,j) + stressp_3(i,j)
         ssigp1  =(stressp_1(i,j) + stressp_3(i,j))*p055
         ssigp2  =(stressp_2(i,j) + stressp_4(i,j))*p055

         ssigmn  = stressm_1(i,j) + stressm_2(i,j)
         ssigms  = stressm_3(i,j) + stressm_4(i,j)
         ssigme  = stressm_1(i,j) + stressm_4(i,j)
         ssigmw  = stressm_2(i,j) + stressm_3(i,j)
         ssigm1  =(stressm_1(i,j) + stressm_3(i,j))*p055
         ssigm2  =(stressm_2(i,j) + stressm_4(i,j))*p055

         ssig12n = stress12_1(i,j) + stress12_2(i,j)
         ssig12s = stress12_3(i,j) + stress12_4(i,j)
         ssig12e = stress12_1(i,j) + stress12_4(i,j)
         ssig12w = stress12_2(i,j) + stress12_3(i,j)
         ssig121 =(stress12_1(i,j) + stress12_3(i,j))*p111
         ssig122 =(stress12_2(i,j) + stress12_4(i,j))*p111

         csigpne = p111*stressp_1(i,j) + ssigp2 + p027*stressp_3(i,j)
         csigpnw = p111*stressp_2(i,j) + ssigp1 + p027*stressp_4(i,j)
         csigpsw = p111*stressp_3(i,j) + ssigp2 + p027*stressp_1(i,j)
         csigpse = p111*stressp_4(i,j) + ssigp1 + p027*stressp_2(i,j)
         
         csigmne = p111*stressm_1(i,j) + ssigm2 + p027*stressm_3(i,j)
         csigmnw = p111*stressm_2(i,j) + ssigm1 + p027*stressm_4(i,j)
         csigmsw = p111*stressm_3(i,j) + ssigm2 + p027*stressm_1(i,j)
         csigmse = p111*stressm_4(i,j) + ssigm1 + p027*stressm_2(i,j)
         
         csig12ne = p222*stress12_1(i,j) + ssig122 &
                  + p055*stress12_3(i,j)
         csig12nw = p222*stress12_2(i,j) + ssig121 &
                  + p055*stress12_4(i,j)
         csig12sw = p222*stress12_3(i,j) + ssig122 &
                  + p055*stress12_1(i,j)
         csig12se = p222*stress12_4(i,j) + ssig121 &
                  + p055*stress12_2(i,j)

         str12ew = p5*dxt(i,j)*(p333*ssig12e + p166*ssig12w)
         str12we = p5*dxt(i,j)*(p333*ssig12w + p166*ssig12e)
         str12ns = p5*dyt(i,j)*(p333*ssig12n + p166*ssig12s)
         str12sn = p5*dyt(i,j)*(p333*ssig12s + p166*ssig12n)

      !-----------------------------------------------------------------
      ! for dF/dx (u momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dyt(i,j)*(p333*ssigpn  + p166*ssigps)
         strm_tmp  = p25*dyt(i,j)*(p333*ssigmn  + p166*ssigms)

         ! northeast (i,j)
         strtmp(i,j,1) = -strp_tmp - strm_tmp - str12ew &
              + dxhy(i,j)*(-csigpne + csigmne) + dyhx(i,j)*csig12ne

         ! northwest (i+1,j)
         strtmp(i,j,2) = strp_tmp + strm_tmp - str12we &
              + dxhy(i,j)*(-csigpnw + csigmnw) + dyhx(i,j)*csig12nw

         strp_tmp  = p25*dyt(i,j)*(p333*ssigps  + p166*ssigpn)
         strm_tmp  = p25*dyt(i,j)*(p333*ssigms  + p166*ssigmn)

         ! southeast (i,j+1)
         strtmp(i,j,3) = -strp_tmp - strm_tmp + str12ew &
              + dxhy(i,j)*(-csigpse + csigmse) + dyhx(i,j)*csig12se

         ! southwest (i+1,j+1)
         strtmp(i,j,4) = strp_tmp + strm_tmp + str12we &
              + dxhy(i,j)*(-csigpsw + csigmsw) + dyhx(i,j)*csig12sw

      !-----------------------------------------------------------------
      ! for dF/dy (v momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dxt(i,j)*(p333*ssigpe  + p166*ssigpw)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigme  + p166*ssigmw)

         ! northeast (i,j)
         strtmp(i,j,5) = -strp_tmp + strm_tmp - str12ns &
              - dyhx(i,j)*(csigpne + csigmne) + dxhy(i,j)*csig12ne

         ! southeast (i,j+1)
         strtmp(i,j,6) = strp_tmp - strm_tmp - str12sn &
              - dyhx(i,j)*(csigpse + csigmse) + dxhy(i,j)*csig12se

         strp_tmp  = p25*dxt(i,j)*(p333*ssigpw  + p166*ssigpe)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigmw  + p166*ssigme)

         ! northwest (i+1,j)
         strtmp(i,j,7) = -strp_tmp + strm_tmp + str12ns &
              - dyhx(i,j)*(csigpnw + csigmnw) + dxhy(i,j)*csig12nw

         ! southwest (i+1,j+1)
         strtmp(i,j,8) = strp_tmp - strm_tmp + str12sn &
              - dyhx(i,j)*(csigpsw + csigmsw) + dxhy(i,j)*csig12sw

      enddo                     ! ij

      end subroutine stress_eap

!=======================================================================

! Updates the stress depending on values of strain rate and structure
! tensor and for ksub=ndte it computes closing and sliding rate

      subroutine update_stress_rdg (ksub, ndte, divu, tension, &
                                   shear, a11, a12, &
                                   stressp,  stressm, &
                                   stress12, strength, &
                                   alphar, alphas)

      integer (kind=int_kind), intent(in) :: &
         ksub, &
         ndte

      real (kind=dbl_kind), intent(in) :: &
         a11,     a12,                   &
         divu,    tension, shear,        &
         strength

      real (kind=dbl_kind), intent(out) :: &
         stressp, stressm, stress12, &
         alphar, alphas     

      ! local variables

      integer (kind=int_kind) :: &
         kx ,ky, ka

      real (kind=dbl_kind) :: &
         stemp11r, stemp12r, stemp22r,   &
         stemp11s, stemp12s, stemp22s, &
         a22, Qd11Qd11, Qd11Qd12, Qd12Qd12, &
         Q11Q11, Q11Q12, Q12Q12, &
         dtemp11, dtemp12, dtemp22, &
         rotstemp11r, rotstemp12r, rotstemp22r,   &
         rotstemp11s, rotstemp12s, rotstemp22s, &
         sig11, sig12, sig22, &
         sgprm11, sgprm12, sgprm22, &
         invstressconviso, &
         Angle_denom_gamma,  Angle_denom_alpha, &
         Tany_1, Tany_2, &
         x, y, dx, dy, da, &
         invdx, invdy, invda, invsin, &
         dtemp1, dtemp2, atempprime, &
         kxw, kyw, kaw, &
         puny, pi, pi2, piq, pih

      real (kind=dbl_kind), parameter :: &
         kfriction = 0.45_dbl_kind

      ! tcraig, temporary, should be moved to namelist
      ! turns on interpolation in stress_rdg
      logical(kind=log_kind), parameter :: &
         interpolate_stress_rdg = .false.

      character(len=*), parameter :: subname = '(update_stress_rdg)'

         call icepack_query_parameters(puny_out=puny, &
            pi_out=pi, pi2_out=pi2, piq_out=piq, pih_out=pih)
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
            file=__FILE__, line=__LINE__)

! Factor to maintain the same stress as in EVP (see Section 3)
! Can be set to 1 otherwise

         invstressconviso = c1/(c1+kfriction*kfriction)
         invsin = c1/sin(pi2/c12) * invstressconviso

! compute eigenvalues, eigenvectors and angles for structure tensor, strain rates

! 1) structure tensor

         a22 = c1-a11

! gamma: angle between general coordiantes and principal axis of A
! here Tan2gamma = 2 a12 / (a11 - a22) 

         Q11Q11 = c1
         Q12Q12 = puny
         Q11Q12 = puny

         if((ABS(a11 - a22) > puny).or.(ABS(a12) > puny)) then
           Angle_denom_gamma = sqrt( ( a11 - a22 )*( a11 - a22) + &
                                     c4*a12*a12 )

           Q11Q11 = p5 + ( a11 - a22 )*p5/Angle_denom_gamma   !Cos^2
           Q12Q12 = p5 - ( a11 - a22 )*p5/Angle_denom_gamma   !Sin^2
           Q11Q12 = a12/Angle_denom_gamma                             !CosSin
         endif

! rotation Q*atemp*Q^T
         atempprime = Q11Q11*a11 + c2*Q11Q12*a12 + Q12Q12*a22 

! make first principal value the largest
         atempprime = max(atempprime, c1 - atempprime)

! 2) strain rate

         dtemp11 = p5*(divu + tension)
         dtemp12 = shear*p5
         dtemp22 = p5*(divu - tension)

! here Tan2alpha = 2 dtemp12 / (dtemp11 - dtemp22) 

         Qd11Qd11 = c1
         Qd12Qd12 = puny
         Qd11Qd12 = puny

         if((ABS( dtemp11 - dtemp22) > puny).or.(ABS(dtemp12) > puny)) then
           Angle_denom_alpha = sqrt( ( dtemp11 - dtemp22 )* &
                        ( dtemp11 - dtemp22 ) + c4*dtemp12*dtemp12)

           Qd11Qd11 = p5 + ( dtemp11 - dtemp22 )*p5/Angle_denom_alpha   !Cos^2 
           Qd12Qd12 = p5 - ( dtemp11 - dtemp22 )*p5/Angle_denom_alpha   !Sin^2
           Qd11Qd12 = dtemp12/Angle_denom_alpha                          !CosSin
         endif 

         dtemp1 = Qd11Qd11*dtemp11 + c2*Qd11Qd12*dtemp12 + Qd12Qd12*dtemp22
         dtemp2 = Qd12Qd12*dtemp11 - c2*Qd11Qd12*dtemp12 + Qd11Qd11*dtemp22

! In cos and sin values
         x = c0

         if ((ABS(dtemp1) > puny).or.(ABS(dtemp2) > puny)) then
!           invleng = c1/sqrt(dtemp1*dtemp1 + dtemp2*dtemp2) ! not sure if this is neccessary
!           dtemp1 = dtemp1*invleng
!           dtemp2 = dtemp2*invleng
           if (dtemp1 == c0) then
             x = pih
           else
             x = atan2(dtemp2,dtemp1)
           endif
         endif

!echmod to ensure the angle lies between pi/4 and 9 pi/4
         if (x < piq) x = x + pi2
!echmod require 0 <= x < (nx_yield-1)*dx = 2 pi
!         x = mod(x+pi2, pi2)

! y: angle between major principal axis of strain rate and structure tensor
! y = gamma - alpha
! Expressesed componently with
! Tany = (Singamma*Cosgamma - Sinalpha*Cosgamma)/(Cos^2gamma - Sin^alpha)

         Tany_1 = Q11Q12 - Qd11Qd12
         Tany_2 = Q11Q11 - Qd12Qd12

         y = c0

         if ((ABS(Tany_1) > puny).or.(ABS(Tany_2) > puny)) then
!           invleng = c1/sqrt(Tany_1*Tany_1 + Tany_2*Tany_2) ! not sure if this is neccessary
!           Tany_1 = Tany_1*invleng
!           Tany_2 = Tany_2*invleng
           if (Tany_2 == c0) then
             y = pih
           else
             y = atan2(Tany_1,Tany_2)
           endif
         endif

! to make sure y is between 0 and pi
         if (y > pi) y = y - pi
         if (y < 0)  y = y + pi

! Now calculate updated stress tensor
         dx   = pi/real(nx_yield-1,kind=dbl_kind)
         dy   = pi/real(ny_yield-1,kind=dbl_kind)
         da   = p5/real(na_yield-1,kind=dbl_kind)
         invdx = c1/dx
         invdy = c1/dy
         invda = c1/da

         if (interpolate_stress_rdg) then

! Interpolated lookup

           ! if (x>=9*pi/4) x=9*pi/4-puny; end
           ! if (y>=pi/2)   y=pi/2-puny; end
           ! if (atempprime>=1.0), atempprime=1.0-puny; end

           ! % need 8 coords and 8 weights
           ! % range in kx

           kx  = int((x-piq-pi)*invdx) + 1
           kxw = c1 - ((x-piq-pi)*invdx - (kx-1))

           ky  = int(y*invdy) + 1
           kyw = c1 - (y*invdy - (ky-1))

           ka  = int((atempprime-p5)*invda) + 1
           kaw = c1 - ((atempprime-p5)*invda - (ka-1))

! % Determine sigma_r(A1,Zeta,y) and sigma_s (see Section A1)

           stemp11r =  kxw* kyw      * kaw      * s11r(kx  ,ky  ,ka  ) &
               + (c1-kxw) * kyw      * kaw      * s11r(kx+1,ky  ,ka  ) &
               + kxw      * (c1-kyw) * kaw      * s11r(kx  ,ky+1,ka  ) &
               + kxw      * kyw      * (c1-kaw) * s11r(kx  ,ky  ,ka+1) &
               + (c1-kxw) * (c1-kyw) * kaw      * s11r(kx+1,ky+1,ka  ) &
               + (c1-kxw) * kyw      * (c1-kaw) * s11r(kx+1,ky  ,ka+1) &
               + kxw      * (c1-kyw) * (c1-kaw) * s11r(kx  ,ky+1,ka+1) &
               + (c1-kxw) * (c1-kyw) * (c1-kaw) * s11r(kx+1,ky+1,ka+1)

           stemp12r =  kxw* kyw      * kaw      * s12r(kx  ,ky  ,ka  ) &
               + (c1-kxw) * kyw      * kaw      * s12r(kx+1,ky  ,ka  ) &
               + kxw      * (c1-kyw) * kaw      * s12r(kx  ,ky+1,ka  ) &
               + kxw      * kyw      * (c1-kaw) * s12r(kx  ,ky  ,ka+1) &
               + (c1-kxw) * (c1-kyw) * kaw      * s12r(kx+1,ky+1,ka  ) &
               + (c1-kxw) * kyw      * (c1-kaw) * s12r(kx+1,ky  ,ka+1) &
               + kxw      * (c1-kyw) * (c1-kaw) * s12r(kx  ,ky+1,ka+1) &
               + (c1-kxw) * (c1-kyw) * (c1-kaw) * s12r(kx+1,ky+1,ka+1)

           stemp22r = kxw * kyw      * kaw      * s22r(kx  ,ky  ,ka  ) &
               + (c1-kxw) * kyw      * kaw      * s22r(kx+1,ky  ,ka  ) &
               + kxw      * (c1-kyw) * kaw      * s22r(kx  ,ky+1,ka  ) &
               + kxw      * kyw      * (c1-kaw) * s22r(kx  ,ky  ,ka+1) &
               + (c1-kxw) * (c1-kyw) * kaw      * s22r(kx+1,ky+1,ka  ) &
               + (c1-kxw) * kyw      * (c1-kaw) * s22r(kx+1,ky  ,ka+1) &
               + kxw      * (c1-kyw) * (c1-kaw) * s22r(kx  ,ky+1,ka+1) &
               + (c1-kxw) * (c1-kyw) * (c1-kaw) * s22r(kx+1,ky+1,ka+1)

           stemp11s =  kxw* kyw      * kaw      * s11s(kx  ,ky  ,ka  ) &
               + (c1-kxw) * kyw      * kaw      * s11s(kx+1,ky  ,ka  ) &
               + kxw      * (c1-kyw) * kaw      * s11s(kx  ,ky+1,ka  ) &
               + kxw      * kyw      * (c1-kaw) * s11s(kx  ,ky  ,ka+1) &
               + (c1-kxw) * (c1-kyw) * kaw      * s11s(kx+1,ky+1,ka  ) &
               + (c1-kxw) * kyw      * (c1-kaw) * s11s(kx+1,ky  ,ka+1) &
               + kxw      * (c1-kyw) * (c1-kaw) * s11s(kx  ,ky+1,ka+1) &
               + (c1-kxw) * (c1-kyw) * (c1-kaw) * s11s(kx+1,ky+1,ka+1)

           stemp12s =  kxw* kyw      * kaw      * s12s(kx  ,ky  ,ka  ) &
               + (c1-kxw) * kyw      * kaw      * s12s(kx+1,ky  ,ka  ) &
               + kxw      * (c1-kyw) * kaw      * s12s(kx  ,ky+1,ka  ) &
               + kxw      * kyw      * (c1-kaw) * s12s(kx  ,ky  ,ka+1) &
               + (c1-kxw) * (c1-kyw) * kaw      * s12s(kx+1,ky+1,ka  ) &
               + (c1-kxw) * kyw      * (c1-kaw) * s12s(kx+1,ky  ,ka+1) &
               + kxw      * (c1-kyw) * (c1-kaw) * s12s(kx  ,ky+1,ka+1) &
               + (c1-kxw) * (c1-kyw) * (c1-kaw) * s12s(kx+1,ky+1,ka+1)

           stemp22s =  kxw* kyw      * kaw      * s22s(kx  ,ky  ,ka  ) &
               + (c1-kxw) * kyw      * kaw      * s22s(kx+1,ky  ,ka  ) &
               + kxw      * (c1-kyw) * kaw      * s22s(kx  ,ky+1,ka  ) &
               + kxw      * kyw      * (c1-kaw) * s22s(kx  ,ky  ,ka+1) &
               + (c1-kxw) * (c1-kyw) * kaw      * s22s(kx+1,ky+1,ka  ) &
               + (c1-kxw) * kyw      * (c1-kaw) * s22s(kx+1,ky  ,ka+1) &
               + kxw      * (c1-kyw) * (c1-kaw) * s22s(kx  ,ky+1,ka+1) &
               + (c1-kxw) * (c1-kyw) * (c1-kaw) * s22s(kx+1,ky+1,ka+1)

         else
           kx = int((x-piq-pi)*invdx) + 1
           ky = int(y*invdy) + 1
           ka = int((atempprime-p5)*invda) + 1

! Determine sigma_r(A1,Zeta,y) and sigma_s (see Section A1) 
           stemp11r = s11r(kx,ky,ka)     
           stemp12r = s12r(kx,ky,ka)
           stemp22r = s22r(kx,ky,ka)

           stemp11s = s11s(kx,ky,ka)
           stemp12s = s12s(kx,ky,ka)
           stemp22s = s22s(kx,ky,ka)
         endif

! Calculate mean ice stress over a collection of floes (Equation 3)

         stressp  = strength*(stemp11r + kfriction*stemp11s &
                            + stemp22r + kfriction*stemp22s) * invsin
         stress12 = strength*(stemp12r + kfriction*stemp12s) * invsin 
         stressm  = strength*(stemp11r + kfriction*stemp11s &
                            - stemp22r - kfriction*stemp22s) * invsin

! Back - rotation of the stress from principal axes into general coordinates

! Update stress
         sig11 = p5*(stressp + stressm)
         sig12 = stress12
         sig22 = p5*(stressp - stressm)

         sgprm11 = Q11Q11*sig11 + Q12Q12*sig22 -        c2*Q11Q12 *sig12 
         sgprm12 = Q11Q12*sig11 - Q11Q12*sig22 + (Q11Q11 - Q12Q12)*sig12
         sgprm22 = Q12Q12*sig11 + Q11Q11*sig22 +        c2*Q11Q12 *sig12 

         stressp  = sgprm11 + sgprm22
         stress12 = sgprm12
         stressm  = sgprm11 - sgprm22

! Compute ridging and sliding functions in general coordinates (Equation 11)
         if (ksub == ndte) then
           rotstemp11r = Q11Q11*stemp11r - c2*Q11Q12* stemp12r &
                       + Q12Q12*stemp22r
           rotstemp12r = Q11Q11*stemp12r +    Q11Q12*(stemp11r-stemp22r) &
                       - Q12Q12*stemp12r
           rotstemp22r = Q12Q12*stemp11r + c2*Q11Q12* stemp12r & 
                       + Q11Q11*stemp22r

           rotstemp11s = Q11Q11*stemp11s - c2*Q11Q12* stemp12s &
                       + Q12Q12*stemp22s
           rotstemp12s = Q11Q11*stemp12s +    Q11Q12*(stemp11s-stemp22s) &
                       - Q12Q12*stemp12s
           rotstemp22s = Q12Q12*stemp11s + c2*Q11Q12* stemp12s & 
                       + Q11Q11*stemp22s

           alphar =  rotstemp11r*dtemp11 + c2*rotstemp12r*dtemp12 &
                   + rotstemp22r*dtemp22
           alphas =  rotstemp11s*dtemp11 + c2*rotstemp12s*dtemp12 &
                   + rotstemp22s*dtemp22
         endif

      end subroutine update_stress_rdg

!=======================================================================

! Solves evolution equation for structure tensor (A19, A20) 

      subroutine stepa  (nx_block,   ny_block,       &
                         dtei,       icellt,         &
                         indxti,     indxtj,         &
                         a11, a12,                   &
                         a11_1, a11_2, a11_3, a11_4, &
                         a12_1, a12_2, a12_3, a12_4, &
                         stressp_1,  stressp_2,      &
                         stressp_3,  stressp_4,      &
                         stressm_1,  stressm_2,      &
                         stressm_3,  stressm_4,      &
                         stress12_1, stress12_2,     &
                         stress12_3, stress12_4)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellt                ! no. of cells where icetmask = 1

      real (kind=dbl_kind), intent(in) :: &
         dtei        ! 1/dte, where dte is subcycling timestep (1/s)

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         ! ice stress tensor (kg/s^2) in each corner of T cell
         stressp_1, stressp_2, stressp_3, stressp_4, & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4, & ! sigma11-sigma22
         stress12_1, stress12_2, stress12_3, stress12_4    ! sigma12

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         ! structure tensor () in each corner of T cell
         a11, a12, a11_1, a11_2, a11_3, a11_4, & ! components of 
         a12_1, a12_2, a12_3, a12_4              ! structure tensor ()

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
        mresult11, mresult12, &
        dteikth, p5kth

      real (kind=dbl_kind), parameter :: &
        kth  = p2*p001             

      character(len=*), parameter :: subname = '(stepa)'

      dteikth = c1 / (dtei + kth)
      p5kth = p5 * kth

      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

! ne 
         call calc_ffrac(1, stressp_1(i,j), stressm_1(i,j),  &
                           stress12_1(i,j),                  &
                           a11_1(i,j),                       &
                           mresult11)

         call calc_ffrac(2, stressp_1(i,j), stressm_1(i,j),  &
                           stress12_1(i,j),                  &
                           a12_1(i,j),                       &
                           mresult12)

         a11_1(i,j) = (a11_1(i,j)*dtei + p5kth - mresult11) * dteikth ! implicit
         a12_1(i,j) = (a12_1(i,j)*dtei - mresult12) * dteikth ! implicit

 
! nw
         call calc_ffrac(1, stressp_2(i,j), stressm_2(i,j),  &
                           stress12_2(i,j),                  &
                           a11_2(i,j),                       &
                           mresult11)

         call calc_ffrac(2, stressp_2(i,j), stressm_2(i,j),  &
                           stress12_2(i,j),                  &
                           a12_2(i,j),                       &
                           mresult12)

         a11_2(i,j) = (a11_2(i,j)*dtei + p5kth - mresult11) * dteikth ! implicit
         a12_2(i,j) = (a12_2(i,j)*dtei - mresult12) * dteikth ! implicit

! sw
         call calc_ffrac(1, stressp_3(i,j), stressm_3(i,j),  &
                           stress12_3(i,j),                  &
                           a11_3(i,j),                       &
                           mresult11)

         call calc_ffrac(2, stressp_3(i,j), stressm_3(i,j),  &
                           stress12_3(i,j),                  &
                           a12_3(i,j),                       &
                           mresult12)

         a11_3(i,j) = (a11_3(i,j)*dtei + p5kth - mresult11) * dteikth ! implicit
         a12_3(i,j) = (a12_3(i,j)*dtei - mresult12) * dteikth ! implicit
                                
! se
         call calc_ffrac(1, stressp_4(i,j), stressm_4(i,j),  &
                           stress12_4(i,j),                  &
                           a11_4(i,j),                       &
                           mresult11)

         call calc_ffrac(2, stressp_4(i,j), stressm_4(i,j),  &
                           stress12_4(i,j),                  &
                           a12_4(i,j),                       &
                           mresult12)

         a11_4(i,j) = (a11_4(i,j)*dtei + p5kth - mresult11) * dteikth ! implicit
         a12_4(i,j) = (a12_4(i,j)*dtei - mresult12) * dteikth ! implicit

! average structure tensor

         a11(i,j) = p25*(a11_1(i,j) + a11_2(i,j) + a11_3(i,j) + a11_4(i,j))
         a12(i,j) = p25*(a12_1(i,j) + a12_2(i,j) + a12_3(i,j) + a12_4(i,j))
               
      enddo                     ! ij
      
      end subroutine stepa

!=======================================================================

! computes term in evolution equation for structure tensor which determines
! the ice floe re-orientation due to fracture
! Eq. 7: Ffrac = -kf(A-S) or = 0 depending on sigma_1 and sigma_2

      subroutine calc_ffrac (blockno, stressp, stressm, &
                             stress12,                  &
                             a1x,                       &
                             mresult)

      integer(kind=int_kind), intent(in) :: &
         blockno

      real (kind=dbl_kind), intent(in) :: &
         stressp, stressm, stress12, a1x

      real (kind=dbl_kind), intent(out) :: &
         mresult

      ! local variables

      real (kind=dbl_kind) :: &
         sigma11, sigma12, sigma22, &
         gamma, sigma_1, sigma_2, pih, &
         Q11, Q12, Q11Q11, Q11Q12, Q12Q12

      real (kind=dbl_kind), parameter :: &
         kfrac = p001, &
         threshold = c3*p1

      character(len=*), parameter :: subname = '(calc_ffrac)'

       call icepack_query_parameters(pih_out=pih)
       call icepack_warnings_flush(nu_diag)
       if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

       sigma11 = p5*(stressp+stressm) 
       sigma12 = stress12 
       sigma22 = p5*(stressp-stressm) 

       if ((sigma11-sigma22) == c0) then
         gamma = p5*(pih)
       else
         gamma = p5*atan2((c2*sigma12),(sigma11-sigma22))
       endif

! rotate tensor to get into sigma principal axis

       Q11 = cos(gamma)
       Q12 = sin(gamma)

       Q11Q11 = Q11*Q11
       Q11Q12 = Q11*Q12
       Q12Q12 = Q12*Q12

       sigma_1 = Q11Q11*sigma11 + c2*Q11Q12*sigma12  &
               + Q12Q12*sigma22 ! S(1,1)
       sigma_2 = Q12Q12*sigma11 - c2*Q11Q12*sigma12  &
               + Q11Q11*sigma22 ! S(2,2)

! Pure divergence
       if ((sigma_1 >= c0).and.(sigma_2 >= c0))  then
         mresult = c0

! Unconfined compression: cracking of blocks not along the axial splitting direction
! which leads to the loss of their shape, so we again model it through diffusion
       elseif ((sigma_1 >= c0).and.(sigma_2 < c0))  then
         if (blockno == 1) mresult = kfrac * (a1x - Q12Q12)
         if (blockno == 2) mresult = kfrac * (a1x + Q11Q12)

! Shear faulting
       elseif (sigma_2 == c0) then
         mresult = c0
       elseif ((sigma_1 <= c0).and.(sigma_1/sigma_2 <= threshold)) then
         if (blockno == 1) mresult = kfrac * (a1x - Q12Q12)
         if (blockno == 2) mresult = kfrac * (a1x + Q11Q12)

! Horizontal spalling
       else 
         mresult = c0
       endif

      end subroutine calc_ffrac

!=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
!=======================================================================

! Dumps all values needed for a restart

      subroutine write_restart_eap ()

      use ice_restart, only: write_restart_field

      ! local variables

      logical (kind=log_kind) :: diag

      character(len=*), parameter :: subname = '(write_restart_eap)'

      diag = .true.

      !-----------------------------------------------------------------
      ! structure tensor 
      !-----------------------------------------------------------------

      call write_restart_field(nu_dump_eap,0,a11_1,'ruf8','a11_1',1,diag)
      call write_restart_field(nu_dump_eap,0,a11_3,'ruf8','a11_3',1,diag)
      call write_restart_field(nu_dump_eap,0,a11_2,'ruf8','a11_2',1,diag)
      call write_restart_field(nu_dump_eap,0,a11_4,'ruf8','a11_4',1,diag)

      call write_restart_field(nu_dump_eap,0,a12_1,'ruf8','a12_1',1,diag)
      call write_restart_field(nu_dump_eap,0,a12_3,'ruf8','a12_3',1,diag)
      call write_restart_field(nu_dump_eap,0,a12_2,'ruf8','a12_2',1,diag)
      call write_restart_field(nu_dump_eap,0,a12_4,'ruf8','a12_4',1,diag)

      end subroutine write_restart_eap

!=======================================================================

! Reads all values needed for elastic anisotropic plastic dynamics restart

      subroutine read_restart_eap()

      use ice_blocks, only: nghost
      use ice_boundary, only: ice_HaloUpdate_stress
      use ice_constants, only:  &
          field_loc_center, field_type_scalar
      use ice_domain, only: nblocks, halo_info
      use ice_grid, only: grid_type
      use ice_restart, only: read_restart_field

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk

      logical (kind=log_kind) :: &
         diag

      character(len=*), parameter :: subname = '(read_restart_eap)'

      diag = .true.

      !-----------------------------------------------------------------
      ! Structure tensor must be read and scattered in pairs in order
      ! to properly match corner values across a tripole grid cut.
      !-----------------------------------------------------------------
         if (my_task == master_task) write(nu_diag,*) &
           'structure tensor restart data'
    
      call read_restart_field(nu_restart_eap,0,a11_1,'ruf8', &
           'a11_1',1,diag,field_loc_center,field_type_scalar) ! a11_1
      call read_restart_field(nu_restart_eap,0,a11_3,'ruf8', &
           'a11_3',1,diag,field_loc_center,field_type_scalar) ! a11_3
      call read_restart_field(nu_restart_eap,0,a11_2,'ruf8', &
           'a11_2',1,diag,field_loc_center,field_type_scalar) ! a11_2
      call read_restart_field(nu_restart_eap,0,a11_4,'ruf8', &
           'a11_4',1,diag,field_loc_center,field_type_scalar) ! a11_4

      call read_restart_field(nu_restart_eap,0,a12_1,'ruf8', &
           'a12_1',1,diag,field_loc_center,field_type_scalar) ! a12_1
      call read_restart_field(nu_restart_eap,0,a12_3,'ruf8', &
           'a12_3',1,diag,field_loc_center,field_type_scalar) ! a12_3
      call read_restart_field(nu_restart_eap,0,a12_2,'ruf8', &
           'a12_2',1,diag,field_loc_center,field_type_scalar) ! a12_2
      call read_restart_field(nu_restart_eap,0,a12_4,'ruf8', &
           'a12_4',1,diag,field_loc_center,field_type_scalar) ! a12_4

      if (trim(grid_type) == 'tripole') then

      call ice_HaloUpdate_stress(a11_1, a11_3, halo_info, &
                                 field_loc_center,  field_type_scalar)
      call ice_HaloUpdate_stress(a11_3, a11_1, halo_info, &
                                 field_loc_center,  field_type_scalar)
      call ice_HaloUpdate_stress(a11_2, a11_4, halo_info, &
                                 field_loc_center,  field_type_scalar)
      call ice_HaloUpdate_stress(a11_4, a11_2, halo_info, &
                                 field_loc_center,  field_type_scalar)
      call ice_HaloUpdate_stress(a12_1, a12_3, halo_info, &
                                 field_loc_center,  field_type_scalar)
      call ice_HaloUpdate_stress(a12_3, a12_1, halo_info, &
                                 field_loc_center,  field_type_scalar)
      call ice_HaloUpdate_stress(a12_2, a12_4, halo_info, &
                                 field_loc_center,  field_type_scalar)
      call ice_HaloUpdate_stress(a12_4, a12_2, halo_info, &
                                 field_loc_center,  field_type_scalar)

      endif

      !-----------------------------------------------------------------
      ! Ensure unused values in west and south ghost cells are 0
      !-----------------------------------------------------------------

         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, nghost
            do i = 1, nx_block
               a11_1 (i,j,iblk) = c0
               a11_2 (i,j,iblk) = c0
               a11_3 (i,j,iblk) = c0
               a11_4 (i,j,iblk) = c0
               a12_1 (i,j,iblk) = c0
               a12_2 (i,j,iblk) = c0
               a12_3 (i,j,iblk) = c0
               a12_4 (i,j,iblk) = c0
            enddo
            enddo
            do j = 1, ny_block
            do i = 1, nghost
               a11_1 (i,j,iblk) = c0
               a11_2 (i,j,iblk) = c0
               a11_3 (i,j,iblk) = c0
               a11_4 (i,j,iblk) = c0
               a12_1 (i,j,iblk) = c0
               a12_2 (i,j,iblk) = c0
               a12_3 (i,j,iblk) = c0
               a12_4 (i,j,iblk) = c0
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO

      end subroutine read_restart_eap

!=======================================================================

      end module ice_dyn_eap

!=======================================================================
