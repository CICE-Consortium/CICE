!  SVN:$Id: ice_dyn_evp.F90 1228 2017-05-23 21:33:34Z tcraig $
!=======================================================================
!
! Elastic-viscous-plastic sea ice dynamics model
! Computes ice velocity and deformation
!
! See:
!
! Hunke, E. C., and J. K. Dukowicz (1997). An elastic-viscous-plastic model
! for sea ice dynamics. {\em J. Phys. Oceanogr.}, {\bf 27}, 1849--1867.
!
! Hunke, E. C. (2001).  Viscous-Plastic Sea Ice Dynamics with the EVP Model:
! Linearization Issues. {\em Journal of Computational Physics}, {\bf 170},
! 18--38.
!
! Hunke, E. C., and J. K. Dukowicz (2002).  The Elastic-Viscous-Plastic
! Sea Ice Dynamics Model in General Orthogonal Curvilinear Coordinates
! on a Sphere---Incorporation of Metric Terms. {\em Monthly Weather Review},
! {\bf 130}, 1848--1865.
!
! Hunke, E. C., and J. K. Dukowicz (2003).  The sea ice momentum
! equation in the free drift regime.  Los Alamos Tech. Rep. LA-UR-03-2219.
!
! Bouillon, S., T. Fichefet, V. Legat and G. Madec (submitted 2013).  The 
! revised elastic-viscous-plastic method.  Ocean Modelling.
!
! author: Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb (LANL)
! 2004: Block structure added by William Lipscomb
! 2005: Removed boundary calls for stress arrays (WHL)
! 2006: Streamlined for efficiency by Elizabeth Hunke
!       Converted to free source form (F90)

      module ice_dyn_evp

      use ice_kinds_mod
      use ice_constants, only: field_loc_center, field_loc_NEcorner, &
          field_type_scalar, field_type_vector
      use ice_constants, only: c0, c4, p027, p055, p111, p166, &
          p2, p222, p25, p333, p5, c1
      use ice_dyn_shared, only: stepu, evp_prep1, evp_prep2, evp_finish, &
          ndte, yield_curve, ecci, denom1, arlx1i, fcor_blk, uvel_init,  &
          vvel_init, basal_stress_coeff, basalstress, Ktens
      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_ice_strength
#ifdef CICE_IN_NEMO
      use icepack_intfc, only: calc_strair
#endif

      implicit none
      private
      public :: evp

!=======================================================================

      contains

!=======================================================================

! Elastic-viscous-plastic dynamics driver
!
#ifdef CICE_IN_NEMO
! Wind stress is set during this routine from the values supplied
! via NEMO (unless calc_strair is true).  These values are supplied 
! rotated on u grid and multiplied by aice.  strairxT = 0 in this 
! case so operations in evp_prep1 are pointless but carried out to 
! minimise code changes.
#endif
!
! author: Elizabeth C. Hunke, LANL

      subroutine evp (dt)

      use ice_arrays_column, only: Cdn_ocn
      use ice_boundary, only: ice_halo, ice_HaloMask, ice_HaloUpdate, &
          ice_HaloDestroy, ice_HaloUpdate_stress
      use ice_blocks, only: block, get_block, nx_block, ny_block
      use ice_domain, only: nblocks, blocks_ice, halo_info, maskhalo_dyn
      use ice_domain_size, only: max_blocks, ncat
      use ice_flux, only: rdg_conv, rdg_shear, strairxT, strairyT, &
          strairx, strairy, uocn, vocn, ss_tltx, ss_tlty, iceumask, fm, &
          strtltx, strtlty, strocnx, strocny, strintx, strinty, &
          strocnxT, strocnyT, strax, stray, &
          Cbu, taubx, tauby, hwater, &
          stressp_1, stressp_2, stressp_3, stressp_4, &
          stressm_1, stressm_2, stressm_3, stressm_4, &
          stress12_1, stress12_2, stress12_3, stress12_4
      use ice_grid, only: tmask, umask, dxt, dyt, dxhy, dyhx, cxp, cyp, cxm, cym, &
          tarear, uarear, tinyarea, to_ugrid, t2ugrid_vector, u2tgrid_vector, &
          grid_type
      use ice_state, only: aice, vice, vsno, uvel, vvel, divu, shear, &
          aice_init, aice0, aicen, vicen, strength
      use ice_timers, only: timer_dynamics, timer_bound, &
          ice_timer_start, ice_timer_stop
#ifdef DMI_EVP
      use ice_grid, only: HTE,HTN
      use evp_kernel1d
      use ice_communicate, only: my_task, master_task
!      use ice_domain_size, only: nx_global, ny_global 
#else
#ifdef DMI_TEST_EVP
      use ice_grid, only: HTE,HTN
#endif
#endif
#ifdef DMI_EVP_TEST
      use ice_calendar, only: istep
#endif

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

      integer (kind=int_kind), dimension (nx_block,ny_block,max_blocks) :: &
         icetmask, &  ! ice extent mask (T-cell)
         halomask     ! generic halo mask

      type (ice_halo) :: &
         halo_info_mask !  ghost cell update info for masked halo

      type (block) :: &
         this_block           ! block information for current block
      
      call ice_timer_start(timer_dynamics) ! dynamics

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      allocate(fld2(nx_block,ny_block,2,max_blocks))

       ! This call is needed only if dt changes during runtime.
!      call set_evp_parameters (dt)

      !-----------------------------------------------------------------
      ! boundary updates
      ! commented out because the ghost cells are freshly 
      ! updated after cleanup_itd
      !-----------------------------------------------------------------

!      call ice_timer_start(timer_bound)
!      call ice_HaloUpdate (aice,              halo_info, &
!                           field_loc_center,  field_type_scalar)
!      call ice_HaloUpdate (vice,              halo_info, &
!                           field_loc_center,  field_type_scalar)
!      call ice_HaloUpdate (vsno,              halo_info, &
!                           field_loc_center,  field_type_scalar)
!      call ice_timer_stop(timer_bound)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks

         do j = 1, ny_block 
         do i = 1, nx_block 
            rdg_conv (i,j,iblk) = c0 
            rdg_shear(i,j,iblk) = c0 
            divu (i,j,iblk) = c0 
            shear(i,j,iblk) = c0 
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

         call evp_prep1 (nx_block,           ny_block,           & 
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

#ifdef CICE_IN_NEMO
      !----------------------------------------------------------------
      ! Set wind stress to values supplied via NEMO
      ! This wind stress is rotated on u grid and multiplied by aice
      !----------------------------------------------------------------
      if (.not. calc_strair) then       
         strairx(:,:,:) = strax(:,:,:)
         strairy(:,:,:) = stray(:,:,:)
      else
#endif
      call t2ugrid_vector(strairx)
      call t2ugrid_vector(strairy)
#ifdef CICE_IN_NEMO
      endif      
#endif

! tcraig, tcx, threading here leads to some non-reproducbile results and failures in icepack_ice_strength
! need to do more debugging
      !$TCXOMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! more preparation for dynamics
      !-----------------------------------------------------------------

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         call evp_prep2 (nx_block,             ny_block,             & 
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
                         Cbu       (:,:,iblk))

      !-----------------------------------------------------------------
      ! ice strength
      !-----------------------------------------------------------------

         strength(:,:,iblk) = c0  ! initialize
         do ij = 1, icellt(iblk)
            i = indxti(ij, iblk)
            j = indxtj(ij, iblk)
            call icepack_ice_strength (ncat,                 &
                                      aice    (i,j,  iblk), & 
                                      vice    (i,j,  iblk), & 
                                      aice0   (i,j,  iblk), & 
                                      aicen   (i,j,:,iblk), &  
                                      vicen   (i,j,:,iblk), & 
                                      strength(i,j,  iblk) )
         enddo  ! ij

         ! load velocity into array for boundary updates
         fld2(:,:,1,iblk) = uvel(:,:,iblk)
         fld2(:,:,2,iblk) = vvel(:,:,iblk)

      enddo  ! iblk
      !$TCXOMP END PARALLEL DO

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) &
         call abort_ice(error_message="subname", &
         file=__FILE__, line=__LINE__)

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (strength,           halo_info, &
                           field_loc_center,   field_type_scalar)
      ! velocities may have changed in evp_prep2
      call ice_HaloUpdate (fld2,               halo_info, &
                           field_loc_NEcorner, field_type_vector)
      call ice_timer_stop(timer_bound)

      ! unload
      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
         uvel(:,:,iblk) = fld2(:,:,1,iblk)
         vvel(:,:,iblk) = fld2(:,:,2,iblk)
      enddo
      !$OMP END PARALLEL DO

      if (maskhalo_dyn) then
         call ice_timer_start(timer_bound)
         halomask = 0
         where (iceumask) halomask = 1
         call ice_HaloUpdate (halomask,          halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_timer_stop(timer_bound)
         call ice_HaloMask(halo_info_mask, halo_info, halomask)
      endif

#ifdef DMI_EVP
!      na=icellt(1)    ! MHRI/DMI: assume iblk==1
!      nx=nx_global+2  ! MHRI/DMI: assume iblk==1  or nx_block
!      ny=ny_global+2  ! MHRI/DMI: assume iblk==1  or ny_block
!      call evp_copyin(nx_global+2,ny_global+2,icellt(1),               &
!      write(*,*)'Entering evp_kernel'
      if (my_task == master_task) then
        call evp_copyin(nx_block,ny_block,icellt(1),                    &
          icetmask,iceumask, HTE,HTN,                                   &
          cdn_ocn,aiu,uocn,vocn,waterx,watery,forcex,forcey,            &
          umassdti,fm,uarear,tarear,strintx,strinty,uvel_init,vvel_init,&
          strength,uvel,vvel,dxt,dyt,dxhy,dyhx,cyp,cxp,cym,cxm,tinyarea,&
          stressp_1 ,stressp_2, stressp_3, stressp_4,                   &
          stressm_1 ,stressm_2, stressm_3, stressm_4,                   &
          stress12_1,stress12_2,stress12_3,stress12_4                   )
        call evp_kernel()
        call evp_copyout(nx_block,ny_block,                             &
          uvel,vvel, strintx,strinty, strocnx,strocny,                  &
          stressp_1, stressp_2, stressp_3, stressp_4,                   &
          stressm_1, stressm_2, stressm_3, stressm_4,                   &
          stress12_1,stress12_2,stress12_3,stress12_4,                  &
          divu,rdg_conv,rdg_shear,shear                                 )
      endif

#else

#ifdef DMI_TEST_EVP
      !------------------------
      write(*,*) 'DMI: Save input files for testing evp kernels'
      write(*,*) 'na = icellt = ',icellt
      ! Use j as temporary input-unit
      j=999

      !------------------------
      open(j, file='EVP_v0input.bin', form='unformatted', access='stream',&
          action='write', status='replace', iostat=i)
      write(j,iostat=i)      &
                    icellt,     &
        indxti    , indxtj,     &
        uvel      , vvel,       &
        dxt       , dyt,        &
        dxhy      , dyhx,       &
        cxp       , cyp,        &
        cxm       , cym,        &
        tarear    , tinyarea,   &
        strength  ,             &
        stressp_1 , stressp_2,  &
        stressp_3 , stressp_4,  &
        stressm_1 , stressm_2,  &
        stressm_3 , stressm_4,  &
        stress12_1, stress12_2, &
        stress12_3, stress12_4, &
        shear     , divu,       &
!MHRI        prs_sig   ,             &
        rdg_conv  , rdg_shear
      close(j)

      !------------------------
      open(j, file='EVP_v0input2.bin', form='unformatted', access='stream',&
          action='write', status='replace', iostat=i)
      write(j,iostat=i)      &
       icellu   , Cdn_ocn  , &
       indxui   , indxuj   , &
       aiu      ,            &
       uocn     , vocn     , &
       waterx   , watery   , &
       forcex   , forcey   , &
       umassdti , fm       , &
       uarear   ,            &
       strocnx  , strocny  , &
       strintx  , strinty  , &
       uvel_init, vvel_init
      close(j)

      !------------------------
      open(j,file='EVP_input2D.bin', form='unformatted', access='stream', &
             action='write', iostat=i)
      if (i/=0) stop 'Problem opening file input2D.bin'
      write(j,iostat=i)                                                &
        icetmask, iceumask,                                            &
        uvel_init, vvel_init, uvel, vvel,                              &
        dxt, dyt, dxhy, dyhx, cxp, cyp, cxm, cym,                      &
        tarear, tinyarea, strength,                                    &
        stressp_1, stressp_2, stressp_3, stressp_4,                    &
        stressm_1, stressm_2, stressm_3, stressm_4,                    &
        stress12_1, stress12_2, stress12_3, stress12_4,                &
        uarear,cdn_ocn,aiu, uocn, vocn, waterx, watery, forcex, forcey,&
        umassdti, fm, strintx, strinty
      if (i/=0) stop 'Problem writing file EVP_input2D.bin'
      close(j)

      !------------------------
      open(j,file='EVP_input2D_Grid.bin',form='unformatted', &
             access='stream', action='write', iostat=i)
      if (i/=0) stop 'Problem opening file EVP_input2D_Grid.bin'
      write(j,iostat=i) HTE,HTN
      if (i/=0) stop 'Problem reading file EVP_input2D_Grid.bin'
      close(j)

#endif

      do ksub = 1,ndte        ! subcycling

      !-----------------------------------------------------------------
      ! stress tensor equation, total surface stress
      !-----------------------------------------------------------------

         !$OMP PARALLEL DO PRIVATE(iblk,strtmp)
         do iblk = 1, nblocks

!            if (trim(yield_curve) == 'ellipse') then
               call stress (nx_block,             ny_block,             & 
                            ksub,                 icellt(iblk),         & 
                            indxti      (:,iblk), indxtj      (:,iblk), & 
                            uvel      (:,:,iblk), vvel      (:,:,iblk), &     
                            dxt       (:,:,iblk), dyt       (:,:,iblk), & 
                            dxhy      (:,:,iblk), dyhx      (:,:,iblk), & 
                            cxp       (:,:,iblk), cyp       (:,:,iblk), & 
                            cxm       (:,:,iblk), cym       (:,:,iblk), & 
                            tarear    (:,:,iblk), tinyarea  (:,:,iblk), & 
                            strength  (:,:,iblk),                       & 
                            stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), & 
                            stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), & 
                            stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), & 
                            stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), & 
                            stress12_1(:,:,iblk), stress12_2(:,:,iblk), & 
                            stress12_3(:,:,iblk), stress12_4(:,:,iblk), & 
                            shear     (:,:,iblk), divu      (:,:,iblk), & 
                            rdg_conv  (:,:,iblk), rdg_shear (:,:,iblk), & 
                            strtmp    (:,:,:) )
!            endif               ! yield_curve

#ifdef DMI_TEST_EVP
      write(*,*)'-----------------------------------------------'
      write(*,*)'Statistics for stepu relevant 2D variables'
      write(*,*)'-----------------------------------------------'
      write(*,'(a11,3i25)') '      indxui:',   minval(indxui(:,1)), maxval(indxui(:,1)),sum(indxui(:,1))
      write(*,'(a11,3i25)') '      indxuj:',   minval(indxuj(:,1)), maxval(indxuj(:,1)),sum(indxuj(:,1))
      write(*,'(a11,3f25.8)') ' uvel_init:',minval( uvel_init(:,:,1)),maxval( uvel_init(:,:,1)) ,sum( uvel_init(:,:,1))
      write(*,'(a11,3f25.8)') ' vvel_init:',minval( vvel_init(:,:,1)),maxval( vvel_init(:,:,1)) ,sum( vvel_init(:,:,1))
      write(*,'(a11,3f25.8)') '   strintx:',   minval(strintx(:,:,1)), maxval(strintx(:,:,1)),sum(strintx(:,:,1))
      write(*,'(a11,3f25.8)') '   strinty:',   minval(strinty(:,:,1)), maxval(strinty(:,:,1)),sum(strinty(:,:,1))
      write(*,'(a11,3f25.8)') '   cdn_ocn:',   minval( cdn_ocn(:,:,1)),maxval( cdn_ocn(:,:,1)),sum( cdn_ocn(:,:,1))
      write(*,'(a11,3f25.8)') '       aiu:',      minval( aiu(:,:,1)),      maxval( aiu(:,:,1)),sum( aiu(:,:,1))
      write(*,'(a11,3f25.8)') '      uocn:',     minval( uocn(:,:,1)),     maxval( uocn(:,:,1)),sum( uocn(:,:,1))
      write(*,'(a11,3f25.8)') '      vocn:',     minval( vocn(:,:,1)),     maxval( vocn(:,:,1)),sum( vocn(:,:,1))
      write(*,'(a11,3f25.8)') '    waterx:',   minval( waterx(:,:,1)),   maxval( waterx(:,:,1)),sum( waterx(:,:,1))
      write(*,'(a11,3f25.8)') '    watery:',   minval( watery(:,:,1)),   maxval( watery(:,:,1)),sum( watery(:,:,1))
      write(*,'(a11,3f25.8)') '    forcex:',   minval( forcex(:,:,1)),   maxval( forcex(:,:,1)),sum( forcex(:,:,1))
      write(*,'(a11,3f25.8)') '    forcey:',   minval( forcey(:,:,1)),   maxval( forcey(:,:,1)),sum( forcey(:,:,1))
      write(*,'(a11,3f25.8)') '  umassdti:', minval( umassdti(:,:,1)),maxval( umassdti(:,:,1)),sum( umassdti(:,:,1))
      write(*,'(a11,3f25.8)') '        fm:',       minval( fm(:,:,1)),       maxval( fm(:,:,1)),sum( fm(:,:,1))
      write(*,'(a11,3f25.8)') '    uarear:',   minval( uarear(:,:,1)),maxval( uarear(:,:,1)),sum( uarear(:,:,1))
      write(*,'(a11,3f25.8)') '      Cbu:',   minval( Cbu(:,:,1)),maxval( Cbu(:,:,1)),sum( Cbu(:,:,1))
      write(*,*)'-----------------------------------------------'
      write(*,*)'END OF Statistics for stepu relevant 2D variables'
      write(*,*)'-----------------------------------------------'
#endif

      !-----------------------------------------------------------------
      ! basal stress calculation (landfast ice)
      !-----------------------------------------------------------------
      
            if (basalstress) then        
               call basal_stress_coeff (nx_block,       ny_block,       &
                                        icellu  (iblk),                 &
                                        indxui(:,iblk), indxuj(:,iblk), &
                                        vice(:,:,iblk), aice(:,:,iblk), &
                                        hwater(:,:,iblk),               &
                                        uvel(:,:,iblk), vvel(:,:,iblk), &
                                        Cbu(:,:,iblk)) 
            endif

      !-----------------------------------------------------------------
      ! momentum equation
      !-----------------------------------------------------------------

            call stepu (nx_block,            ny_block,           &
                        icellu       (iblk), Cdn_ocn (:,:,iblk), & 
                        indxui     (:,iblk), indxuj    (:,iblk), & 
                        aiu      (:,:,iblk), strtmp  (:,:,:),    & 
                        uocn     (:,:,iblk), vocn    (:,:,iblk), &     
                        waterx   (:,:,iblk), watery  (:,:,iblk), & 
                        forcex   (:,:,iblk), forcey  (:,:,iblk), & 
                        umassdti (:,:,iblk), fm      (:,:,iblk), & 
                        uarear   (:,:,iblk),                     & 
                        strocnx  (:,:,iblk), strocny (:,:,iblk), & 
                        strintx  (:,:,iblk), strinty (:,:,iblk), & 
                        uvel_init(:,:,iblk), vvel_init(:,:,iblk),&
                        uvel     (:,:,iblk), vvel    (:,:,iblk), &
                        Cbu      (:,:,iblk))

            ! load velocity into array for boundary updates
            fld2(:,:,1,iblk) = uvel(:,:,iblk)
            fld2(:,:,2,iblk) = vvel(:,:,iblk)
         enddo
         !$OMP END PARALLEL DO

#ifdef DMI_TEST_EVP
     !------------------------------------
     ! stress output
     ! Use j as temporary input-unit
     j=999

     !------------------------------------
     open(j, file='EVP_output.bin', form='unformatted', access='stream',&
            action='write', status='replace', iostat=i)
     write(j,iostat=i)                           &
       stressp_1, stressp_2, stressp_3, stressp_4 , &
       stressm_1, stressm_2, stressm_3, stressm_4 , &
       stress12_1,stress12_2,stress12_3,stress12_4, &
       Cbu
!       prs_sig  , &
!       shear    , &
!       divu     , &
!       rdg_conv , &
!       rdg_shear, &
!       strtmp
     close(j)
     !------------------------
     ! stepu output
     open(j, file='EVP_output2.bin', form='unformatted', access='stream',&
            action='write', status='replace', iostat=i)
     write(j,iostat=i)       &
       Cdn_ocn  , aiu     , &
       uocn     , vocn    , &
       waterx   , watery  , &
       forcex   , forcey  , &
       umassdti , fm      , &
       uarear   ,           &
       strocnx  , strocny , &
       strintx  , strinty , &
       uvel_init, vvel_init, Cbu
     close(j)
     !--------------------
     write(*,*)'AFTER CALLING stress+stepu --'
     write(*,'(a11,3f25.8)') ' stressp_1:',   minval(stressp_1(:,:,1)), maxval(stressp_1(:,:,1)),sum(stressp_1(:,:,1))
     write(*,'(a11,3f25.8)') ' stressp_2:',   minval(stressp_2(:,:,1)), maxval(stressp_2(:,:,1)),sum(stressp_2(:,:,1))
     write(*,'(a11,3f25.8)') ' stressp_3:',   minval(stressp_3(:,:,1)), maxval(stressp_3(:,:,1)),sum(stressp_3(:,:,1))
     write(*,'(a11,3f25.8)') ' stressp_4:',   minval(stressp_4(:,:,1)), maxval(stressp_4(:,:,1)),sum(stressp_4(:,:,1))
     write(*,'(a11,3f25.8)') ' stressm_1:',   minval(stressm_1(:,:,1)), maxval(stressm_1(:,:,1)),sum(stressm_1(:,:,1))
     write(*,'(a11,3f25.8)') ' stressm_2:',   minval(stressm_2(:,:,1)), maxval(stressm_2(:,:,1)),sum(stressm_2(:,:,1))
     write(*,'(a11,3f25.8)') ' stressm_3:',   minval(stressm_3(:,:,1)), maxval(stressm_3(:,:,1)),sum(stressm_3(:,:,1))
     write(*,'(a11,3f25.8)') ' stressm_4:',   minval(stressm_4(:,:,1)), maxval(stressm_4(:,:,1)),sum(stressm_4(:,:,1))
     write(*,'(a11,3f25.8)') 'stress12_1:',    minval(stress12_1(:,:,1)), maxval(stress12_1(:,:,1)),sum(stress12_1(:,:,1))
     write(*,'(a11,3f25.8)') 'stress12_2:',    minval(stress12_2(:,:,1)), maxval(stress12_2(:,:,1)),sum(stress12_2(:,:,1))
     write(*,'(a11,3f25.8)') 'stress12_3:',    minval(stress12_3(:,:,1)), maxval(stress12_3(:,:,1)),sum(stress12_3(:,:,1))
     write(*,'(a11,3f25.8)') 'stress12_4:',    minval(stress12_4(:,:,1)), maxval(stress12_4(:,:,1)),sum(stress12_4(:,:,1))
!MHRI     write(*,'(a11,3f25.8)') '   prs_sig:',   minval(prs_sig(:,:,1)), maxval(prs_sig(:,:,1)),sum(prs_sig(:,:,1))
     write(*,'(a11,3f25.8)') '      str1:',   minval(strtmp(:,:,1)), maxval(strtmp(:,:,1)),sum(strtmp(:,:,1))
     write(*,'(a11,3f25.8)') '      str2:',   minval(strtmp(:,:,2)), maxval(strtmp(:,:,2)),sum(strtmp(:,:,2))
     write(*,'(a11,3f25.8)') '      str3:',   minval(strtmp(:,:,3)), maxval(strtmp(:,:,3)),sum(strtmp(:,:,3))
     write(*,'(a11,3f25.8)') '      str4:',   minval(strtmp(:,:,4)), maxval(strtmp(:,:,4)),sum(strtmp(:,:,4))
     write(*,'(a11,3f25.8)') '      str5:',   minval(strtmp(:,:,5)), maxval(strtmp(:,:,5)),sum(strtmp(:,:,5))
     write(*,'(a11,3f25.8)') '      str6:',   minval(strtmp(:,:,6)), maxval(strtmp(:,:,6)),sum(strtmp(:,:,6))
     write(*,'(a11,3f25.8)') '      str7:',   minval(strtmp(:,:,7)), maxval(strtmp(:,:,7)),sum(strtmp(:,:,7))
     write(*,'(a11,3f25.8)') '      str8:',   minval(strtmp(:,:,8)), maxval(strtmp(:,:,8)),sum(strtmp(:,:,8))
     write(*,'(a11,3f25.8)') '   strocnx:',   minval(strocnx(:,:,1)), maxval(strocnx(:,:,1)),sum(strocnx(:,:,1))
     write(*,'(a11,3f25.8)') '   strocny:',   minval(strocny(:,:,1)), maxval(strocny(:,:,1)),sum(strocny(:,:,1))
     write(*,'(a11,3f25.8)') '   strintx:',   minval(strintx(:,:,1)), maxval(strintx(:,:,1)),sum(strintx(:,:,1))
     write(*,'(a11,3f25.8)') '   strinty:',   minval(strinty(:,:,1)), maxval(strinty(:,:,1)),sum(strinty(:,:,1))
     write(*,'(a11,3f25.8)') '      uvel:',   minval(uvel(:,:,1)), maxval(uvel(:,:,1)),sum(uvel(:,:,1))
     write(*,'(a11,3f25.8)') '      vvel:',   minval(vvel(:,:,1)), maxval(vvel(:,:,1)),sum(vvel(:,:,1))
     write(*,'(a11,3f25.8)') '       Cbu:',   minval( Cbu(:,:,1)), maxval( Cbu(:,:,1)),sum( Cbu(:,:,1))
     write(*,*)'AFTER CALLING --'
     !--------------------
     stop
#endif

         call ice_timer_start(timer_bound)
         if (maskhalo_dyn) then
            call ice_HaloUpdate (fld2,               halo_info_mask, &
                                 field_loc_NEcorner, field_type_vector)
         else
            call ice_HaloUpdate (fld2,               halo_info, &
                                 field_loc_NEcorner, field_type_vector)
         endif
         call ice_timer_stop(timer_bound)

         ! unload
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            uvel(:,:,iblk) = fld2(:,:,1,iblk)
            vvel(:,:,iblk) = fld2(:,:,2,iblk)
         enddo
         !$OMP END PARALLEL DO
         
      enddo                     ! subcycling
#endif
      
      ! calculate basal stress component for outputs
      if ( basalstress ) then
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            taubx(:,:,iblk) = -Cbu(:,:,iblk)*uvel(:,:,iblk)
            tauby(:,:,iblk) = -Cbu(:,:,iblk)*vvel(:,:,iblk)
         enddo
         !$OMP END PARALLEL DO
      endif

      deallocate(fld2)
      if (maskhalo_dyn) call ice_HaloDestroy(halo_info_mask)

      ! Force symmetry across the tripole seam
      if (trim(grid_type) == 'tripole') then
      if (maskhalo_dyn) then
         !-------------------------------------------------------
         ! set halomask to zero because ice_HaloMask always keeps
         ! local copies AND tripole zipper communication
         !-------------------------------------------------------
         halomask = 0
         call ice_HaloMask(halo_info_mask, halo_info, halomask)

         call ice_HaloUpdate_stress(stressp_1, stressp_3, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_3, stressp_1, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_2, stressp_4, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_4, stressp_2, halo_info_mask, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stressm_1, stressm_3, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_3, stressm_1, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_2, stressm_4, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_4, stressm_2, halo_info_mask, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stress12_1, stress12_3, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_3, stress12_1, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_2, stress12_4, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_4, stress12_2, halo_info_mask, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloDestroy(halo_info_mask)
      else
         call ice_HaloUpdate_stress(stressp_1, stressp_3, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_3, stressp_1, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_2, stressp_4, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_4, stressp_2, halo_info, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stressm_1, stressm_3, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_3, stressm_1, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_2, stressm_4, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_4, stressm_2, halo_info, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stress12_1, stress12_3, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_3, stress12_1, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_2, stress12_4, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_4, stress12_2, halo_info, &
                              field_loc_center,  field_type_scalar)
      endif   ! maskhalo
      endif   ! tripole

      !-----------------------------------------------------------------
      ! ice-ocean stress
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks

         call evp_finish                               & 
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

#ifdef DMI_EVP_TEST
     !------------------------------------
     ! uvel,vvel output
     ! Use j as temporary input-unit
     j=999

     !------------------------------------
     open(j, file='EVP_output_uvel_vvel.bin', form='unformatted', access='stream',&
            action='write', status='replace', iostat=i)
     write(j,iostat=i)                           &
         uvel,vvel
     close(j)
     write(*,*)'-----------------------------------------------'
     write(*,*)'Statistics for uvel,vvel 2D variables: istep=',istep
     write(*,*)'-----------------------------------------------'
     write(*,'(a11,3f25.10)') '       uvel:',   minval(uvel), maxval(uvel),sum(uvel)
     write(*,'(a11,3f25.10)') '       vvel:',   minval(vvel), maxval(vvel),sum(vvel)
#endif

      end subroutine evp

!=======================================================================

! Computes the rates of strain and internal stress components for
! each of the four corners on each T-grid cell.
! Computes stress terms for the momentum equation
!
! author: Elizabeth C. Hunke, LANL

      subroutine stress (nx_block,   ny_block,   & 
                         ksub,       icellt,     & 
                         indxti,     indxtj,     & 
                         uvel,       vvel,       & 
                         dxt,        dyt,        & 
                         dxhy,       dyhx,       & 
                         cxp,        cyp,        & 
                         cxm,        cym,        & 
                         tarear,     tinyarea,   & 
                         strength,               & 
                         stressp_1,  stressp_2,  & 
                         stressp_3,  stressp_4,  & 
                         stressm_1,  stressm_2,  & 
                         stressm_3,  stressm_4,  & 
                         stress12_1, stress12_2, & 
                         stress12_3, stress12_4, & 
                         shear,      divu,       & 
                         rdg_conv,   rdg_shear,  & 
                         str )

      integer (kind=int_kind), intent(in) :: & 
         nx_block, ny_block, & ! block dimensions
         ksub              , & ! subcycling step
         icellt                ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), & 
         intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj       ! compressed index in j-direction

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
         tarear   , & ! 1/tarea
         tinyarea     ! puny*tarea

      real (kind=dbl_kind), dimension (nx_block,ny_block), & 
         intent(inout) :: &
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4    ! sigma12

      real (kind=dbl_kind), dimension (nx_block,ny_block), & 
         intent(inout) :: &
         shear    , & ! strain rate II component (1/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         rdg_conv , & ! convergence term for ridging (1/s)
         rdg_shear    ! shear term for ridging (1/s)

      real (kind=dbl_kind), dimension(nx_block,ny_block,8), & 
         intent(out) :: &
         str          ! stress combinations

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        Deltane, Deltanw, Deltase, Deltasw        , & ! Delt
        c0ne, c0nw, c0se, c0sw                    , & ! useful combinations
        c1ne, c1nw, c1se, c1sw                    , &
        ssigpn, ssigps, ssigpe, ssigpw            , &
        ssigmn, ssigms, ssigme, ssigmw            , &
        ssig12n, ssig12s, ssig12e, ssig12w        , &
        ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122, &
        csigpne, csigpnw, csigpse, csigpsw        , &
        csigmne, csigmnw, csigmse, csigmsw        , &
        csig12ne, csig12nw, csig12se, csig12sw    , &
        str12ew, str12we, str12ns, str12sn        , &
        strp_tmp, strm_tmp, tmp

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      str(:,:,:) = c0

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
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

         ! shearing strain rate  =  e_12
         shearne = -cym(i,j)*vvel(i  ,j  ) - dyt(i,j)*vvel(i-1,j  ) &
                 -  cxm(i,j)*uvel(i  ,j  ) - dxt(i,j)*uvel(i  ,j-1)
         shearnw = -cyp(i,j)*vvel(i-1,j  ) + dyt(i,j)*vvel(i  ,j  ) &
                 -  cxm(i,j)*uvel(i-1,j  ) - dxt(i,j)*uvel(i-1,j-1)
         shearsw = -cyp(i,j)*vvel(i-1,j-1) + dyt(i,j)*vvel(i  ,j-1) &
                 -  cxp(i,j)*uvel(i-1,j-1) + dxt(i,j)*uvel(i-1,j  )
         shearse = -cym(i,j)*vvel(i  ,j-1) - dyt(i,j)*vvel(i-1,j-1) &
                 -  cxp(i,j)*uvel(i  ,j-1) + dxt(i,j)*uvel(i  ,j  )
         
         ! Delta (in the denominator of zeta, eta)
         Deltane = sqrt(divune**2 + ecci*(tensionne**2 + shearne**2))
         Deltanw = sqrt(divunw**2 + ecci*(tensionnw**2 + shearnw**2))
         Deltase = sqrt(divuse**2 + ecci*(tensionse**2 + shearse**2))
         Deltasw = sqrt(divusw**2 + ecci*(tensionsw**2 + shearsw**2))

#ifndef DMI_TEST_EVP
      !-----------------------------------------------------------------
      ! on last subcycle, save quantities for mechanical redistribution
      !-----------------------------------------------------------------
         if (ksub == ndte) then
            divu(i,j) = p25*(divune + divunw + divuse + divusw) * tarear(i,j)
            tmp = p25*(Deltane + Deltanw + Deltase + Deltasw)   * tarear(i,j)
            rdg_conv(i,j)  = -min(divu(i,j),c0)
            rdg_shear(i,j) = p5*(tmp-abs(divu(i,j))) 

            ! diagnostic only
            ! shear = sqrt(tension**2 + shearing**2)
            shear(i,j) = p25*tarear(i,j)*sqrt( &
                 (tensionne + tensionnw + tensionse + tensionsw)**2 &
                +  (shearne +   shearnw +   shearse +   shearsw)**2)

         endif
#endif

      !-----------------------------------------------------------------
      ! replacement pressure/Delta                   ! kg/s
      ! save replacement pressure for principal stress calculation
      !-----------------------------------------------------------------
         c0ne = strength(i,j)/max(Deltane,tinyarea(i,j))
         c0nw = strength(i,j)/max(Deltanw,tinyarea(i,j))
         c0sw = strength(i,j)/max(Deltasw,tinyarea(i,j))
         c0se = strength(i,j)/max(Deltase,tinyarea(i,j))

         c1ne = c0ne*arlx1i
         c1nw = c0nw*arlx1i
         c1sw = c0sw*arlx1i
         c1se = c0se*arlx1i

         c0ne = c1ne*ecci
         c0nw = c1nw*ecci
         c0sw = c1sw*ecci
         c0se = c1se*ecci

      !-----------------------------------------------------------------
      ! the stresses                            ! kg/s^2
      ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
      !-----------------------------------------------------------------

         stressp_1(i,j) = (stressp_1(i,j) + c1ne*(divune*(c1+Ktens) - Deltane*(c1-Ktens))) &
                          * denom1
         stressp_2(i,j) = (stressp_2(i,j) + c1nw*(divunw*(c1+Ktens) - Deltanw*(c1-Ktens))) &
                          * denom1
         stressp_3(i,j) = (stressp_3(i,j) + c1sw*(divusw*(c1+Ktens) - Deltasw*(c1-Ktens))) &
                          * denom1
         stressp_4(i,j) = (stressp_4(i,j) + c1se*(divuse*(c1+Ktens) - Deltase*(c1-Ktens))) &
                          * denom1

         stressm_1(i,j) = (stressm_1(i,j) + c0ne*tensionne*(c1+Ktens)) * denom1
         stressm_2(i,j) = (stressm_2(i,j) + c0nw*tensionnw*(c1+Ktens)) * denom1
         stressm_3(i,j) = (stressm_3(i,j) + c0sw*tensionsw*(c1+Ktens)) * denom1
         stressm_4(i,j) = (stressm_4(i,j) + c0se*tensionse*(c1+Ktens)) * denom1

         stress12_1(i,j) = (stress12_1(i,j) + c0ne*shearne*p5*(c1+Ktens)) * denom1
         stress12_2(i,j) = (stress12_2(i,j) + c0nw*shearnw*p5*(c1+Ktens)) * denom1
         stress12_3(i,j) = (stress12_3(i,j) + c0sw*shearsw*p5*(c1+Ktens)) * denom1
         stress12_4(i,j) = (stress12_4(i,j) + c0se*shearse*p5*(c1+Ktens)) * denom1

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
         str(i,j,1) = -strp_tmp - strm_tmp - str12ew &
              + dxhy(i,j)*(-csigpne + csigmne) + dyhx(i,j)*csig12ne

         ! northwest (i+1,j)
         str(i,j,2) = strp_tmp + strm_tmp - str12we &
              + dxhy(i,j)*(-csigpnw + csigmnw) + dyhx(i,j)*csig12nw

         strp_tmp  = p25*dyt(i,j)*(p333*ssigps  + p166*ssigpn)
         strm_tmp  = p25*dyt(i,j)*(p333*ssigms  + p166*ssigmn)

         ! southeast (i,j+1)
         str(i,j,3) = -strp_tmp - strm_tmp + str12ew &
              + dxhy(i,j)*(-csigpse + csigmse) + dyhx(i,j)*csig12se

         ! southwest (i+1,j+1)
         str(i,j,4) = strp_tmp + strm_tmp + str12we &
              + dxhy(i,j)*(-csigpsw + csigmsw) + dyhx(i,j)*csig12sw

      !-----------------------------------------------------------------
      ! for dF/dy (v momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dxt(i,j)*(p333*ssigpe  + p166*ssigpw)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigme  + p166*ssigmw)

         ! northeast (i,j)
         str(i,j,5) = -strp_tmp + strm_tmp - str12ns &
              - dyhx(i,j)*(csigpne + csigmne) + dxhy(i,j)*csig12ne

         ! southeast (i,j+1)
         str(i,j,6) = strp_tmp - strm_tmp - str12sn &
              - dyhx(i,j)*(csigpse + csigmse) + dxhy(i,j)*csig12se

         strp_tmp  = p25*dxt(i,j)*(p333*ssigpw  + p166*ssigpe)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigmw  + p166*ssigme)

         ! northwest (i+1,j)
         str(i,j,7) = -strp_tmp + strm_tmp + str12ns &
              - dyhx(i,j)*(csigpnw + csigmnw) + dxhy(i,j)*csig12nw

         ! southwest (i+1,j+1)
         str(i,j,8) = strp_tmp - strm_tmp + str12sn &
              - dyhx(i,j)*(csigpsw + csigmsw) + dxhy(i,j)*csig12sw

      enddo                     ! ij

      end subroutine stress

!=======================================================================

      end module ice_dyn_evp

!=======================================================================
