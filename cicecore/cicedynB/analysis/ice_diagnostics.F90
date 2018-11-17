!=======================================================================

! Diagnostic information output during run
!
! authors: Elizabeth C. Hunke, LANL
!          Bruce P. Briegleb, NCAR
!
! 2004: Block structure added by William Lipscomb
! 2006: Converted to free source form (F90) by Elizabeth Hunke

      module ice_diagnostics

      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0
      use ice_calendar, only: istep1
      use ice_fileunits, only: nu_diag
      use ice_fileunits, only: flush_fileunit
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_max_aero
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_query_tracer_flags
      use icepack_intfc, only: icepack_query_tracer_indices

      implicit none
      private
      public :: runtime_diags, init_mass_diags, init_diags, &
                print_state, print_points_state, diagnostic_abort


      ! diagnostic output file
      character (len=char_len), public :: diag_file

      ! point print data

      logical (kind=log_kind), public :: &
         print_points     , & ! if true, print point data
         print_global         ! if true, print global data

      integer (kind=int_kind), parameter, public :: &
         npnt = 2             ! total number of points to be printed

      ! Set to true to identify unstable fast-moving ice.
      logical (kind=log_kind), parameter ::  &
         check_umax = .false. ! if true, check for speed > umax_stab

      real (kind=dbl_kind), parameter :: &
         umax_stab   = 1.0_dbl_kind , & ! ice speed threshold for instability (m/s)
         aice_extmin = 0.15_dbl_kind    ! min aice value for ice extent calc
 
      real (kind=dbl_kind), dimension(npnt), public :: &
         latpnt           , & !  latitude of diagnostic points
         lonpnt               ! longitude of diagnostic points

      integer (kind=int_kind) :: &
         iindx            , & ! i index for points
         jindx            , & ! j index for points
         bindx                ! block index for points

      ! for water and heat budgets
      real (kind=dbl_kind), dimension(npnt) :: &
         pdhi             , & ! change in mean ice thickness (m)
         pdhs             , & ! change in mean snow thickness (m)
         pde                  ! change in ice and snow energy (W m-2)

      real (kind=dbl_kind), dimension(npnt), public :: &
         plat, plon           ! latitude, longitude of points

      integer (kind=int_kind), dimension(npnt), public :: &
         piloc, pjloc, pbloc, pmloc  ! location of diagnostic points

      ! for hemispheric water and heat budgets
      real (kind=dbl_kind) :: &
         totmn            , & ! total ice/snow water mass (nh)
         totms            , & ! total ice/snow water mass (sh)
         totmin           , & ! total ice water mass (nh)
         totmis           , & ! total ice water mass (sh)
         toten            , & ! total ice/snow energy (J)
         totes                ! total ice/snow energy (J)

      real (kind=dbl_kind), dimension(icepack_max_aero) :: &
         totaeron         , & ! total aerosol mass
         totaeros             ! total aerosol mass

      ! printing info for routine print_state
      ! iblkp, ip, jp, mtask identify the grid cell to print
!     character (char_len) :: plabel
      integer (kind=int_kind), parameter, public :: &
         check_step = 999999999, & ! begin printing at istep1=check_step
         iblkp = 1, &      ! block number 
         ip = 3, &         ! i index
         jp = 5, &         ! j index
         mtask = 0         ! my_task

!=======================================================================

      contains

!=======================================================================

! Writes diagnostic info (max, min, global sums, etc) to standard out
!
! authors: Elizabeth C. Hunke, LANL
!          Bruce P. Briegleb, NCAR
!          Cecilia M. Bitz, UW

      subroutine runtime_diags (dt)

      use ice_blocks, only: nx_block, ny_block
      use ice_broadcast, only: broadcast_scalar
      use ice_constants, only: c1, c1000, c2, p001, p5, &
          field_loc_center, m2_to_km2
      use ice_domain, only: distrb_info, nblocks
      use ice_domain_size, only: ncat, n_aero, max_blocks
      use ice_flux, only: alvdr, alidr, alvdf, alidf, evap, fsnow, frazil, &
          fswabs, fswthru, flw, flwout, fsens, fsurf, flat, frzmlt_init, frain, fpond, &
          fhocn_ai, fsalt_ai, fresh_ai, frazil_diag, &
          update_ocn_f, Tair, Qa, fsw, fcondtop, meltt, meltb, meltl, snoice, &
          dsnow, congel, sst, sss, Tf, fhocn, &
          swvdr, swvdf, swidr, swidf, &
          alvdr_init, alvdf_init, alidr_init, alidf_init
      use ice_flux_bgc, only: faero_atm, faero_ocn
      use ice_global_reductions, only: global_sum, global_sum_prod, global_maxval
      use ice_grid, only: lmask_n, lmask_s, tarean, tareas
      use ice_state   ! everything
#ifdef CESMCOUPLED
      use ice_prescribed_mod, only: prescribed_ice
#endif

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         i, j, n, iblk, &
         ktherm, &
         nt_tsfc, nt_aero, nt_fbri, nt_apnd, nt_hpnd

      logical (kind=log_kind) :: &
         tr_pond_topo, tr_brine, tr_aero, calc_Tsfc

      real (kind=dbl_kind) :: &
         rhow, rhos, rhoi, puny, awtvdr, awtidr, awtvdf, awtidf, &
         rhofresh, lfresh, lvap, ice_ref_salinity, Tffresh

      ! hemispheric state quantities
      real (kind=dbl_kind) :: &
         umaxn,   hmaxn,   shmaxn,    arean,   snwmxn, extentn, shmaxnt, &
         umaxs,   hmaxs,   shmaxs,    areas,   snwmxs, extents, shmaxst, &
         etotn,   mtotn,   micen,     msnwn,   pmaxn,  ketotn, &
         etots,   mtots,   mices,     msnws,   pmaxs,  ketots, &
         urmsn,   albtotn, arean_alb, mpndn,   ptotn,  spondn, &
         urmss,   albtots, areas_alb, mpnds,   ptots,  sponds 

      ! hemispheric flux quantities
      real (kind=dbl_kind) :: &
         rnn, snn, frzn,  hnetn, fhocnn, fhatmn,  fhfrzn, &
         rns, sns, frzs,  hnets, fhocns, fhatms,  fhfrzs, &
         fswnetn, fswnets, fswdnn, fswdns, swerrn, swerrs, &
         sfsaltn, sfreshn, evpn, fluxn , delmxn,  delmin, &
         sfsalts, sfreshs, evps, fluxs , delmxs,  delmis, &
         delein, werrn, herrn, msltn, delmsltn, serrn, &
         deleis, werrs, herrs, mslts, delmslts, serrs

      ! aerosol diagnostics
      real (kind=dbl_kind), dimension(icepack_max_aero) :: &
         faeran, faeron, aerrn, &
         faeras, faeros, aerrs, &
         aeromx1n, aeromx1s, &
         aerototn, aerotots

      ! fields at diagnostic points
      real (kind=dbl_kind), dimension(npnt) :: &
         paice, pTair, pQa, pfsnow, pfrain, pfsw, pflw, & 
         pTsfc, pevap, pfswabs, pflwout, pflat, pfsens, &
         pfsurf, pfcondtop, psst, psss, pTf, hiavg, hsavg, hbravg, &
         pfhocn, psalt, &
         pmeltt, pmeltb, pmeltl, psnoice, pdsnow, pfrazil, pcongel

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1, work2

      character(len=*), parameter :: subname = '(runtime_diags)'

      call icepack_query_parameters(ktherm_out=ktherm, calc_Tsfc_out=calc_Tsfc)
      call icepack_query_tracer_flags(tr_brine_out=tr_brine, tr_aero_out=tr_aero, &
           tr_pond_topo_out=tr_pond_topo)
      call icepack_query_tracer_indices(nt_fbri_out=nt_fbri, nt_Tsfc_out=nt_Tsfc, &
           nt_aero_out=nt_aero, nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd)
      call icepack_query_parameters(Tffresh_out=Tffresh, rhos_out=rhos, &
           rhow_out=rhow, rhoi_out=rhoi, puny_out=puny, &
           awtvdr_out=awtvdr, awtidr_out=awtidr, awtvdf_out=awtvdf, awtidf_out=awtidf, &
           rhofresh_out=rhofresh, lfresh_out=lfresh, lvap_out=lvap, &
           ice_ref_salinity_out=ice_ref_salinity)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! state of the ice
      !-----------------------------------------------------------------
      ! hemispheric quantities

      ! total ice area
      arean = global_sum(aice, distrb_info, field_loc_center, tarean)
      areas = global_sum(aice, distrb_info, field_loc_center, tareas)
      arean = arean * m2_to_km2
      areas = areas * m2_to_km2

      ! ice extent (= area of grid cells with aice > aice_extmin)
      work1(:,:,:) = c0
      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            if (aice(i,j,iblk) >= aice_extmin) work1(i,j,iblk) = c1
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      extentn = global_sum(work1, distrb_info, field_loc_center, &
                           tarean)
      extents = global_sum(work1, distrb_info, field_loc_center, &
                           tareas)
      extentn = extentn * m2_to_km2
      extents = extents * m2_to_km2

      ! total ice volume
      shmaxn = global_sum(vice, distrb_info, field_loc_center, tarean)
      shmaxs = global_sum(vice, distrb_info, field_loc_center, tareas)

      ! total snow volume
      snwmxn = global_sum(vsno, distrb_info, field_loc_center, tarean)
      snwmxs = global_sum(vsno, distrb_info, field_loc_center, tareas)

      ! total pond volume
      ptotn = c0
      ptots = c0
      if (tr_pond_topo) then
         !$OMP PARALLEL DO PRIVATE(iblk,i,j,n)
         do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            work1(i,j,iblk) = c0
            do n = 1, ncat
               work1(i,j,iblk) = work1(i,j,iblk)  &
                               + aicen(i,j,n,iblk) &
                               * trcrn(i,j,nt_apnd,n,iblk) & 
                               * trcrn(i,j,nt_hpnd,n,iblk)
            enddo
         enddo
         enddo
         enddo
         !$OMP END PARALLEL DO
         ptotn = global_sum(work1, distrb_info, field_loc_center, tarean)
         ptots = global_sum(work1, distrb_info, field_loc_center, tareas)
      endif

      ! total ice-snow kinetic energy
      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            work1(i,j,iblk) = p5 &
                           * (rhos*vsno(i,j,iblk) + rhoi*vice(i,j,iblk)) &
                           * (uvel(i,j,iblk)**2 + vvel(i,j,iblk)**2)
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      ketotn = global_sum(work1, distrb_info, field_loc_center, tarean)
      ketots = global_sum(work1, distrb_info, field_loc_center, tareas)

      ! rms ice speed
      urmsn = c2*ketotn/(rhoi*shmaxn + rhos*snwmxn + puny)
      if (urmsn > puny) then
         urmsn = sqrt(urmsn)
      else
         urmsn = c0
      endif

      urmss = c2*ketots/(rhoi*shmaxs + rhos*snwmxs + puny)
      if (urmss > puny) then
         urmss = sqrt(urmss)
      else
         urmss = c0
      endif

      ! average ice albedo
      ! mask out cells where sun is below horizon (for delta-Eddington)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            work1(i,j,iblk) = alvdr(i,j,iblk)*awtvdr &
                            + alidr(i,j,iblk)*awtidr &
                            + alvdf(i,j,iblk)*awtvdf &
                            + alidf(i,j,iblk)*awtidf
            if (work1(i,j,iblk) > puny) then
               work2(i,j,iblk) = tarean(i,j,iblk)
            else
               work2(i,j,iblk) = c0
            endif
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      
      arean_alb = global_sum(aice, distrb_info, field_loc_center, work2)      

      albtotn = global_sum_prod(aice, work1, distrb_info, &
                                field_loc_center, work2)

      if (arean_alb > c0) then
         albtotn = albtotn / arean_alb
      else
         albtotn = c0
      endif

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            if (work1(i,j,iblk) > puny) then
               work2(i,j,iblk) = tareas(i,j,iblk)
            else
               work2(i,j,iblk) = c0
            endif
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      areas_alb = global_sum(aice, distrb_info, field_loc_center, work2)      

      albtots = global_sum_prod(aice, work1, distrb_info, &
                                field_loc_center, work2)

      if (areas_alb > c0) then
         albtots = albtots / areas_alb
      else
         albtots = c0
      endif

      ! maximum ice volume (= mean thickness including open water)
      hmaxn = global_maxval(vice, distrb_info, lmask_n)
      hmaxs = global_maxval(vice, distrb_info, lmask_s)

      ! maximum ice speed
      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            work1(i,j,iblk) = sqrt(uvel(i,j,iblk)**2 &
                                 + vvel(i,j,iblk)**2)
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      umaxn = global_maxval(work1, distrb_info, lmask_n)
      umaxs = global_maxval(work1, distrb_info, lmask_s)

      ! Write warning message if ice speed is too big
      ! (Ice speeds of ~1 m/s or more usually indicate instability)

      if (check_umax) then
         if (umaxn > umax_stab) then
            !$OMP PARALLEL DO PRIVATE(iblk,i,j)
            do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               if (abs(work1(i,j,iblk) - umaxn) < puny) then
                  write(nu_diag,*) ' '
                  write(nu_diag,*) 'Warning, large ice speed'
                  write(nu_diag,*) 'my_task, iblk, i, j, umaxn:', &
                                    my_task, iblk, i, j, umaxn
               endif
            enddo
            enddo
            enddo
            !$OMP END PARALLEL DO
         elseif (umaxs > umax_stab) then
            !$OMP PARALLEL DO PRIVATE(iblk,i,j)
            do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               if (abs(work1(i,j,iblk) - umaxs) < puny) then
                  write(nu_diag,*) ' '
                  write(nu_diag,*) 'Warning, large ice speed'
                  write(nu_diag,*) 'my_task, iblk, i, j, umaxs:', &
                                    my_task, iblk, i, j, umaxs
               endif
            enddo
            enddo
            enddo
            !$OMP END PARALLEL DO
         endif   ! umax
      endif      ! check_umax

      ! maximum ice strength

      pmaxn = global_maxval(strength, distrb_info, lmask_n)
      pmaxs = global_maxval(strength, distrb_info, lmask_s)

      pmaxn = pmaxn / c1000   ! convert to kN/m
      pmaxs = pmaxs / c1000 

      if (print_global) then

         ! total ice/snow internal energy
         call total_energy (work1)
      
         etotn = global_sum(work1, distrb_info, &
                            field_loc_center, tarean)
         etots = global_sum(work1, distrb_info, &
                            field_loc_center, tareas)

      !-----------------------------------------------------------------
      ! various fluxes
      !-----------------------------------------------------------------
      ! evap, fsens, and flwout need to be multiplied by aice because
      ! regrettably they have been divided by aice for the coupler
      !-----------------------------------------------------------------

         ! evaporation

         evpn = global_sum_prod(evap, aice, distrb_info, &
                                field_loc_center, tarean)
         evps = global_sum_prod(evap, aice, distrb_info, &
                                field_loc_center, tareas)
         evpn = evpn*dt
         evps = evps*dt

         ! total brine tracer
         shmaxnt = c0
         shmaxst = c0
         if (tr_brine) then
         shmaxnt = global_sum(vice(:,:,:)*trcr(:,:,nt_fbri,:), distrb_info, &
                                   field_loc_center, tarean)
         shmaxst = global_sum(vice(:,:,:)*trcr(:,:,nt_fbri,:), distrb_info, &
                                   field_loc_center, tareas)
         endif

         ! salt flux
         sfsaltn = global_sum(fsalt_ai, distrb_info, &
                                   field_loc_center, tarean)
         sfsalts = global_sum(fsalt_ai, distrb_info, &
                                   field_loc_center, tareas)
         sfsaltn = sfsaltn*dt
         sfsalts = sfsalts*dt

         ! fresh water flux
         sfreshn = global_sum(fresh_ai, distrb_info, &
                                   field_loc_center, tarean)
         sfreshs = global_sum(fresh_ai, distrb_info, &
                                   field_loc_center, tareas)
         sfreshn = sfreshn*dt
         sfreshs = sfreshs*dt

         ! pond water flux
         spondn = c0
         sponds = c0
         if (tr_pond_topo) then
         spondn = global_sum(fpond, distrb_info, &
                                   field_loc_center, tarean)
         sponds = global_sum(fpond, distrb_info, &
                                   field_loc_center, tareas)
         spondn = spondn*dt
         sponds = sponds*dt
         endif

         ! ocean heat
         ! Note: fswthru not included because it does not heat ice
         fhocnn = global_sum(fhocn_ai, distrb_info, &
                                  field_loc_center, tarean)
         fhocns = global_sum(fhocn_ai, distrb_info, &
                                  field_loc_center, tareas)

         ! latent heat
         ! You may be wondering, where is the latent heat flux?
         ! It is not included here because it cancels with
         ! the evaporative flux times the enthalpy of the
         ! ice/snow that evaporated.

         ! atmo heat flux
         ! Note: flwout includes the reflected longwave down, needed by the
         !  atmosphere as an upwards radiative boundary condition.
         ! Also note: fswabs includes solar radiation absorbed in ocean,
         !  which must be subtracted here.

         if (calc_Tsfc) then

            !$OMP PARALLEL DO PRIVATE(iblk,i,j)
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  work1(i,j,iblk) = &
                     (fswabs(i,j,iblk) - fswthru(i,j,iblk) &
                    + fsens (i,j,iblk) + flwout (i,j,iblk)) &
                                                  * aice      (i,j,iblk) &
                    + flw   (i,j,iblk) * aice_init (i,j,iblk)
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         else   ! fsurf is computed by atmosphere model 

            !$OMP PARALLEL DO PRIVATE(iblk,i,j)
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  work1(i,j,iblk) = &
                             (fsurf(i,j,iblk) - flat(i,j,iblk)) & 
                              * aice(i,j,iblk)
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         endif     ! calc_Tsfc

         fhatmn = global_sum(work1, distrb_info, &
                             field_loc_center, tarean)
         fhatms = global_sum(work1, distrb_info, &
                             field_loc_center, tareas)
  
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               work1(i,j,iblk) = &
                  fswabs(i,j,iblk) * aice(i,j,iblk)
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO

         fswnetn = global_sum(work1, distrb_info, &
                             field_loc_center, tarean)
         fswnets = global_sum(work1, distrb_info, &
                             field_loc_center, tareas)

         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               work1(i,j,iblk) = (aice_init(i,j,iblk)-alvdr_init(i,j,iblk))*swvdr(i,j,iblk) &
                               + (aice_init(i,j,iblk)-alidr_init(i,j,iblk))*swidr(i,j,iblk) &
                               + (aice_init(i,j,iblk)-alvdf_init(i,j,iblk))*swvdf(i,j,iblk) &
                               + (aice_init(i,j,iblk)-alidf_init(i,j,iblk))*swidf(i,j,iblk)
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO

         fswdnn = global_sum(work1, distrb_info, &
                             field_loc_center, tarean)
         fswdns = global_sum(work1, distrb_info, &
                             field_loc_center, tareas)

         ! freezing potential
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               work1(i,j,iblk) = max(c0,frzmlt_init(i,j,iblk))
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
         fhfrzn = global_sum(work1, distrb_info, &
                             field_loc_center, tarean)
         fhfrzs = global_sum(work1, distrb_info, &
                             field_loc_center, tareas)

         ! rain
         rnn = global_sum_prod(frain, aice_init, distrb_info, &
                               field_loc_center, tarean)
         rns = global_sum_prod(frain, aice_init, distrb_info, &
                               field_loc_center, tareas)
         rnn = rnn*dt
         rns = rns*dt

         ! snow
         snn = global_sum_prod(fsnow, aice_init, distrb_info, &
                               field_loc_center, tarean)
         sns = global_sum_prod(fsnow, aice_init, distrb_info, &
                               field_loc_center, tareas)
         snn = snn*dt
         sns = sns*dt

         ! frazil ice growth !! should not be multiplied by aice
         ! m/step->kg/m^2/s
         work1(:,:,:) = frazil(:,:,:)*rhoi/dt
         if (ktherm == 2 .and. .not.update_ocn_f) &
            work1(:,:,:) = (frazil(:,:,:)-frazil_diag(:,:,:))*rhoi/dt
         frzn = global_sum(work1, distrb_info, &
                           field_loc_center, tarean)
         frzs = global_sum(work1, distrb_info, &
                           field_loc_center, tareas)
         frzn = frzn*dt
         frzs = frzs*dt

         ! ice, snow, pond mass
         micen = rhoi*shmaxn
         msnwn = rhos*snwmxn
         mices = rhoi*shmaxs
         msnws = rhos*snwmxs
         mpndn = rhofresh*ptotn
         mpnds = rhofresh*ptots

         ! total ice, snow and pond mass
         mtotn = micen + msnwn + mpndn
         mtots = mices + msnws + mpnds
  
         ! mass change since beginning of time step
         delmin = mtotn - totmn
         delmis = mtots - totms

         ! ice mass change including frazil ice formation
         delmxn = micen - totmin
         delmxs = mices - totmis
         if (.not. update_ocn_f) then
           ! ice mass change excluding frazil ice formation
           delmxn = delmxn - frzn
           delmxs = delmxs - frzs
         endif

         ! total water flux
         fluxn  = c0
         fluxs  = c0
         if( arean > c0) then
           ! water associated with frazil ice included in fresh
           fluxn = rnn + snn + evpn - sfreshn 
           if (.not. update_ocn_f) then
             fluxn = fluxn + frzn
           endif
         endif
         if( areas > c0) then
           ! water associated with frazil ice included in fresh
           fluxs = rns + sns + evps - sfreshs 
           if (.not. update_ocn_f) then
             fluxs = fluxs + frzs
           endif
         endif

         werrn = (fluxn-delmin)/(mtotn + c1)
         werrs = (fluxs-delmis)/(mtots + c1)

         ! energy change
         delein = etotn - toten
         deleis = etots - totes

         fhatmn = fhatmn + ( - snn * Lfresh + evpn * Lvap ) / dt
         fhatms = fhatms + ( - sns * Lfresh + evps * Lvap ) / dt

         hnetn = (fhatmn - fhocnn - fhfrzn) * dt
         hnets = (fhatms - fhocns - fhfrzs) * dt

         herrn = (hnetn - delein) / (etotn - c1)
         herrs = (hnets - deleis) / (etots - c1)

         swerrn = (fswnetn - fswdnn) / (fswnetn - c1)
         swerrs = (fswnets - fswdns) / (fswnets - c1)

         ! salt mass
         msltn = micen*ice_ref_salinity*p001
         mslts = mices*ice_ref_salinity*p001

         ! change in salt mass
         delmsltn = delmxn*ice_ref_salinity*p001
         delmslts = delmxs*ice_ref_salinity*p001

         ! salt error
         serrn = (sfsaltn + delmsltn) / (msltn + c1)
         serrs = (sfsalts + delmslts) / (mslts + c1)

         ! aerosols
         if (tr_aero) then
         do n = 1, n_aero
            faeran(n) = global_sum_prod(faero_atm(:,:,n,:), aice_init, &
                                        distrb_info, field_loc_center, tarean)
            faeras(n) = global_sum_prod(faero_atm(:,:,n,:), aice_init, &
                                        distrb_info, field_loc_center, tareas)
            faeran(n) = faeran(n)*dt
            faeras(n) = faeras(n)*dt
            faeron(n) = global_sum_prod(faero_ocn(:,:,n,:), aice, &
                                        distrb_info, field_loc_center, tarean)
            faeros(n) = global_sum_prod(faero_ocn(:,:,n,:), aice, &
                                        distrb_info, field_loc_center, tareas)
            faeron(n) = faeron(n)*dt
            faeros(n) = faeros(n)*dt

            !$OMP PARALLEL DO PRIVATE(iblk,i,j)
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  work1(i,j,iblk) = &
                   trcr(i,j,nt_aero  +4*(n-1),iblk)*vsno(i,j,iblk) &
                 + trcr(i,j,nt_aero+1+4*(n-1),iblk)*vsno(i,j,iblk) &
                 + trcr(i,j,nt_aero+2+4*(n-1),iblk)*vice(i,j,iblk) &
                 + trcr(i,j,nt_aero+3+4*(n-1),iblk)*vice(i,j,iblk)
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
            aerototn(n) = global_sum(work1, distrb_info, field_loc_center, tarean)
            aerotots(n) = global_sum(work1, distrb_info, field_loc_center, tareas)
            aeromx1n(n) = global_maxval(work1, distrb_info, lmask_n)
            aeromx1s(n) = global_maxval(work1, distrb_info, lmask_s)

            aerrn(n) = (totaeron(n)-aerototn(n)+faeran(n)-faeron(n)) &
                                 / (aerototn(n) + c1)
            aerrs(n) = (totaeros(n)-aerotots(n)+faeras(n)-faeros(n)) &
                                 / (aerotots(n) + c1)
         enddo ! n_aero
         endif ! tr_aero

      endif                     ! print_global

      if (print_points) then

      !-----------------------------------------------------------------
      ! state of the ice and associated fluxes for 2 defined points
      ! NOTE these are computed for the last timestep only (not avg)
      !-----------------------------------------------------------------

         call total_energy (work1)
         call total_salt   (work2)

         do n = 1, npnt
            if (my_task == pmloc(n)) then
               i = piloc(n)
               j = pjloc(n)
               iblk = pbloc(n)

               pTair(n) = Tair(i,j,iblk) - Tffresh ! air temperature
               pQa(n) = Qa(i,j,iblk)               ! specific humidity
               pfsnow(n) = fsnow(i,j,iblk)*dt/rhos ! snowfall
               pfrain(n) = frain(i,j,iblk)*dt/rhow ! rainfall
               pfsw(n) = fsw(i,j,iblk)             ! shortwave radiation
               pflw(n) = flw(i,j,iblk)             ! longwave radiation
               paice(n) = aice(i,j,iblk)           ! ice area
               
               hiavg(n) = c0                       ! avg snow/ice thickness
               hsavg(n) = c0
               hbravg(n) = c0                      ! avg brine thickness
               if (paice(n) /= c0) then
                  hiavg(n) = vice(i,j,iblk)/paice(n)
                  hsavg(n) = vsno(i,j,iblk)/paice(n)
                  if (tr_brine) hbravg(n) = trcr(i,j,nt_fbri,iblk)* hiavg(n)
               endif
               if (vice(i,j,iblk) /= c0) psalt(n) = work2(i,j,iblk)/vice(i,j,iblk)
               pTsfc(n) = trcr(i,j,nt_Tsfc,iblk)   ! ice/snow sfc temperature
               pevap(n) = evap(i,j,iblk)*dt/rhoi   ! sublimation/condensation
               pfswabs(n) = fswabs(i,j,iblk)       ! absorbed solar flux
               pflwout(n) = flwout(i,j,iblk)       ! outward longwave flux
               pflat(n) = flat(i,j,iblk)           ! latent heat flux
               pfsens(n) = fsens(i,j,iblk)         ! sensible heat flux
               pfsurf(n) = fsurf(i,j,iblk)         ! total sfc heat flux
               pfcondtop(n) = fcondtop(i,j,iblk)   ! top sfc cond flux
               pmeltt(n) = meltt(i,j,iblk)         ! top melt
               pmeltb(n) = meltb(i,j,iblk)         ! bottom melt
               pmeltl(n) = meltl(i,j,iblk)         ! lateral melt
               psnoice(n) = snoice(i,j,iblk)       ! snow ice
               pdsnow(n) = dsnow(i,j,iblk)         ! snow change
               pfrazil(n) = frazil(i,j,iblk)       ! frazil ice
               pcongel(n) = congel(i,j,iblk)       ! congelation ice
               pdhi(n) = vice(i,j,iblk) - pdhi(n)  ! ice thickness change
               pdhs(n) = vsno(i,j,iblk) - pdhs(n)  ! snow thickness change
               pde(n) =-(work1(i,j,iblk)- pde(n))/dt ! ice/snow energy change 
               psst(n) = sst(i,j,iblk)             ! sea surface temperature
               psss(n) = sss(i,j,iblk)             ! sea surface salinity
               pTf(n) = Tf(i,j,iblk)               ! freezing temperature
               pfhocn(n) = -fhocn(i,j,iblk)        ! ocean heat used by ice

            endif  ! my_task = pmloc

            call broadcast_scalar(pTair    (n), pmloc(n))             
            call broadcast_scalar(pQa      (n), pmloc(n))             
            call broadcast_scalar(pfsnow   (n), pmloc(n))             
            call broadcast_scalar(pfrain   (n), pmloc(n))             
            call broadcast_scalar(pfsw     (n), pmloc(n))             
            call broadcast_scalar(pflw     (n), pmloc(n))             
            call broadcast_scalar(paice    (n), pmloc(n))             
            call broadcast_scalar(hsavg    (n), pmloc(n))             
            call broadcast_scalar(hiavg    (n), pmloc(n))              
            call broadcast_scalar(psalt    (n), pmloc(n))
            call broadcast_scalar(hbravg   (n), pmloc(n))
            call broadcast_scalar(pTsfc    (n), pmloc(n))             
            call broadcast_scalar(pevap    (n), pmloc(n))             
            call broadcast_scalar(pfswabs  (n), pmloc(n)) 
            call broadcast_scalar(pflwout  (n), pmloc(n)) 
            call broadcast_scalar(pflat    (n), pmloc(n)) 
            call broadcast_scalar(pfsens   (n), pmloc(n)) 
            call broadcast_scalar(pfsurf   (n), pmloc(n)) 
            call broadcast_scalar(pfcondtop(n), pmloc(n)) 
            call broadcast_scalar(pmeltt   (n), pmloc(n)) 
            call broadcast_scalar(pmeltb   (n), pmloc(n)) 
            call broadcast_scalar(pmeltl   (n), pmloc(n)) 
            call broadcast_scalar(psnoice  (n), pmloc(n)) 
            call broadcast_scalar(pdsnow   (n), pmloc(n)) 
            call broadcast_scalar(pfrazil  (n), pmloc(n)) 
            call broadcast_scalar(pcongel  (n), pmloc(n)) 
            call broadcast_scalar(pdhi     (n), pmloc(n)) 
            call broadcast_scalar(pdhs     (n), pmloc(n)) 
            call broadcast_scalar(pde      (n), pmloc(n)) 
            call broadcast_scalar(psst     (n), pmloc(n)) 
            call broadcast_scalar(psss     (n), pmloc(n)) 
            call broadcast_scalar(pTf      (n), pmloc(n)) 
            call broadcast_scalar(pfhocn   (n), pmloc(n))
            
         enddo                  ! npnt
      endif                     ! print_points

      !-----------------------------------------------------------------
      ! start spewing
      !-----------------------------------------------------------------

      if (my_task == master_task) then

        write(nu_diag,899) 'Arctic','Antarctic'

        write(nu_diag,901) 'total ice area  (km^2) = ',arean,  areas
        write(nu_diag,901) 'total ice extent(km^2) = ',extentn,extents
        write(nu_diag,901) 'total ice volume (m^3) = ',shmaxn, shmaxs
        write(nu_diag,901) 'total snw volume (m^3) = ',snwmxn, snwmxs
        write(nu_diag,901) 'tot kinetic energy (J) = ',ketotn, ketots
        write(nu_diag,900) 'rms ice speed    (m/s) = ',urmsn,  urmss
        write(nu_diag,900) 'average albedo         = ',albtotn,albtots
        write(nu_diag,900) 'max ice volume     (m) = ',hmaxn,  hmaxs
        write(nu_diag,900) 'max ice speed    (m/s) = ',umaxn,  umaxs
        write(nu_diag,900) 'max strength    (kN/m) = ',pmaxn,  pmaxs

        if (print_global) then  ! global diags for conservations checks

#ifdef CESMCOUPLED
         if (prescribed_ice) then
          write (nu_diag,*) '----------------------------'
          write (nu_diag,*)   'This is the prescribed ice option.'
          write (nu_diag,*)   'Heat and water will not be conserved.'
          write (nu_diag,*) '----------------------------'
         endif
#endif

         write(nu_diag,*) '----------------------------'
         write(nu_diag,901) 'arwt rain h2o kg in dt = ',rnn,rns
         write(nu_diag,901) 'arwt snow h2o kg in dt = ',snn,sns
         write(nu_diag,901) 'arwt evap h2o kg in dt = ',evpn,evps
         write(nu_diag,901) 'arwt frzl h2o kg in dt = ',frzn,frzs
         if (tr_pond_topo) &
         write(nu_diag,901) 'arwt fpnd h2o kg in dt = ',spondn,sponds
         write(nu_diag,901) 'arwt frsh h2o kg in dt = ',sfreshn,sfreshs

         write(nu_diag,901) 'arwt ice mass (kg)     = ',micen,mices
         write(nu_diag,901) 'arwt snw mass (kg)     = ',msnwn,msnws
         if (tr_pond_topo) &
         write(nu_diag,901) 'arwt pnd mass (kg)     = ',mpndn,mpnds
 
         write(nu_diag,901) 'arwt tot mass (kg)     = ',mtotn,mtots
         write(nu_diag,901) 'arwt tot mass chng(kg) = ',delmin,delmis
         write(nu_diag,901) 'arwt water flux        = ',fluxn,fluxs
         if (update_ocn_f) then
           write (nu_diag,*) '(=rain+snow+evap-fresh)  '
         else
           write (nu_diag,*) '(=rain+snow+evap+frzl-fresh)  '
         endif
         write(nu_diag,901) 'water flux error       = ',werrn,werrs

         write(nu_diag,*) '----------------------------'
         write(nu_diag,901) 'arwt atm heat flux (W) = ',fhatmn,fhatms
         write(nu_diag,901) 'arwt ocn heat flux (W) = ',fhocnn,fhocns
         write(nu_diag,901) 'arwt frzl heat flux(W) = ',fhfrzn,fhfrzs
         write(nu_diag,901) 'arwt tot energy    (J) = ',etotn,etots
         write(nu_diag,901) 'arwt net heat      (J) = ',hnetn,hnets
         write(nu_diag,901) 'arwt tot energy chng(J)= ',delein,deleis
         write(nu_diag,901) 'arwt heat error        = ',herrn,herrs
         write(nu_diag,*) '----------------------------'
         write(nu_diag,901) 'arwt incoming sw (W)   = ',fswdnn,fswdns
         write(nu_diag,901) 'arwt absorbed sw (W)   = ',fswnetn,fswnets
         write(nu_diag,901) 'arwt swdn error        = ',swerrn,swerrs

         write(nu_diag,*) '----------------------------'
         write(nu_diag,901) 'total brine tr (m^3)   = ',shmaxnt, shmaxst
         write(nu_diag,901) 'arwt salt mass (kg)    = ',msltn,mslts
         write(nu_diag,901) 'arwt salt mass chng(kg)= ',delmsltn, &
                                                        delmslts
         write(nu_diag,901) 'arwt salt flx in dt(kg)= ',sfsaltn, &
                                                        sfsalts
         write(nu_diag,901) 'arwt salt flx error    = ',serrn,serrs

         write(nu_diag,*) '----------------------------'
         if (tr_aero) then
         do n = 1, n_aero
         write(nu_diag,*)   '  aerosol ',n
         write(nu_diag,901) 'faero_atm (kg/m2)      = ', faeran(n), faeras(n)
         write(nu_diag,901) 'faero_ocn (kg/m2)      = ', faeron(n), faeros(n)
         write(nu_diag,901) 'total aero (kg/m2)     = ', aerototn(n), aerotots(n)
         write(nu_diag,901) 'aero error             = ', aerrn(n), aerrs(n)
         write(nu_diag,901) 'maximum aero (kg/m2)   = ', aeromx1n(n),aeromx1s(n)
         enddo
         write(nu_diag,*) '----------------------------'
         endif ! tr_aero

        endif                    ! print_global

       call flush_fileunit(nu_diag)

      !-----------------------------------------------------------------
      ! diagnostics for Arctic and Antarctic points
      !-----------------------------------------------------------------

       if (print_points) then

        write(nu_diag,*) '                         '
        write(nu_diag,902) '       Lat, Long         ',plat(1),plon(1), &
                                                       plat(2),plon(2)
        write(nu_diag,903) '  my_task, iblk, i, j     ', &
                              pmloc(1),pbloc(1),piloc(1),pjloc(1), &
                              pmloc(2),pbloc(2),piloc(2),pjloc(2)
        write(nu_diag,*) '----------atm----------'
        write(nu_diag,900) 'air temperature (C)    = ',pTair(1),pTair(2)
        write(nu_diag,900) 'specific humidity      = ',pQa(1),pQa(2)
        write(nu_diag,900) 'snowfall (m)           = ',pfsnow(1), &
                                                       pfsnow(2)
        write(nu_diag,900) 'rainfall (m)           = ',pfrain(1), &
                                                       pfrain(2)
        if (.not.calc_Tsfc) then
           write(nu_diag,900) 'total surface heat flux= ',pfsurf(1),pfsurf(2)
           write(nu_diag,900) 'top sfc conductive flux= ',pfcondtop(1), &
                                                          pfcondtop(2)
           write(nu_diag,900) 'latent heat flx        = ',pflat(1),pflat(2)
        else
           write(nu_diag,900) 'shortwave radiation sum= ',pfsw(1),pfsw(2)
           write(nu_diag,900) 'longwave radiation     = ',pflw(1),pflw(2)
        endif
        write(nu_diag,*) '----------ice----------'
        write(nu_diag,900) 'area fraction          = ',paice(1),paice(2)
        write(nu_diag,900) 'avg ice thickness (m)  = ',hiavg(1),hiavg(2)
        write(nu_diag,900) 'avg snow depth (m)     = ',hsavg(1),hsavg(2)
        write(nu_diag,900) 'avg salinity (ppt)     = ',psalt(1),psalt(2)
        write(nu_diag,900) 'avg brine thickness (m)= ',hbravg(1),hbravg(2)

        if (calc_Tsfc) then
           write(nu_diag,900) 'surface temperature(C) = ',pTsfc(1),pTsfc(2)
           write(nu_diag,900) 'absorbed shortwave flx = ',pfswabs(1), &
                                                          pfswabs(2)
           write(nu_diag,900) 'outward longwave flx   = ',pflwout(1), &
                                                          pflwout(2)
           write(nu_diag,900) 'sensible heat flx      = ',pfsens(1), &
                                                          pfsens(2)
           write(nu_diag,900) 'latent heat flx        = ',pflat(1),pflat(2)
        endif
        write(nu_diag,900) 'subl/cond (m ice)      = ',pevap(1),pevap(2)
        write(nu_diag,900) 'top melt (m)           = ',pmeltt(1) &
                                                      ,pmeltt(2)
        write(nu_diag,900) 'bottom melt (m)        = ',pmeltb(1) &
                                                      ,pmeltb(2)
        write(nu_diag,900) 'lateral melt (m)       = ',pmeltl(1) &
                                                      ,pmeltl(2)
        write(nu_diag,900) 'new ice (m)            = ',pfrazil(1), &
                                                       pfrazil(2)
        write(nu_diag,900) 'congelation (m)        = ',pcongel(1), &
                                                       pcongel(2)
        write(nu_diag,900) 'snow-ice (m)           = ',psnoice(1), &
                                                       psnoice(2)
        write(nu_diag,900) 'snow change (m)        = ',pdsnow(1), &
                                                       pdsnow(2)
        write(nu_diag,900) 'effective dhi (m)      = ',pdhi(1),pdhi(2)
        write(nu_diag,900) 'effective dhs (m)      = ',pdhs(1),pdhs(2)
        write(nu_diag,900) 'intnl enrgy chng(W/m^2)= ',pde (1),pde (2)
        write(nu_diag,*) '----------ocn----------'
        write(nu_diag,900) 'sst (C)                = ',psst(1),psst(2)
        write(nu_diag,900) 'sss (ppt)              = ',psss(1),psss(2)
        write(nu_diag,900) 'freezing temp (C)      = ',pTf(1),pTf(2)
        write(nu_diag,900) 'heat used (W/m^2)      = ',pfhocn(1), &
                                                       pfhocn(2)

       endif                    ! print_points
      endif                     ! my_task = master_task

  799 format (27x,a24)
  800 format (a25,2x,f24.17)
  801 format (a25,2x,1pe24.17)
  899 format (27x,a24,2x,a24)
  900 format (a25,2x,f24.17,2x,f24.17)
  901 format (a25,2x,1pe24.17,2x,1pe24.17)
  902 format (a25,10x,f6.1,1x,f6.1,9x,f6.1,1x,f6.1)
  903 format (a25,5x,i4,1x,i4,1x,i4,1x,i4,7x,i4,1x,i4,1x,i4,1x,i4)

      end subroutine runtime_diags

!=======================================================================

! Computes global combined ice and snow mass sum
!
! author: Elizabeth C. Hunke, LANL

      subroutine init_mass_diags

      use ice_blocks, only: nx_block, ny_block
      use ice_constants, only: field_loc_center
      use ice_domain, only: distrb_info, nblocks
      use ice_domain_size, only: n_aero, ncat, max_blocks
      use ice_global_reductions, only: global_sum
      use ice_grid, only: tareas, tarean
      use ice_state, only: aicen, vice, vsno, trcrn, trcr

      integer (kind=int_kind) :: n, i, j, iblk, &
         nt_hpnd, nt_apnd, nt_aero

      logical (kind=log_kind) :: &
         tr_aero, tr_pond_topo

      real (kind=dbl_kind) :: &
         shmaxn, snwmxn,  shmaxs, snwmxs, totpn, totps, &
         rhoi, rhos, rhofresh

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      character(len=*), parameter :: subname = '(init_mass_diags)'

      call icepack_query_tracer_flags(tr_aero_out=tr_aero, tr_pond_topo_out=tr_pond_topo)
      call icepack_query_tracer_indices( &
         nt_hpnd_out=nt_hpnd, nt_apnd_out=nt_apnd, nt_aero_out=nt_aero)
      call icepack_query_parameters( &
         rhoi_out=rhoi, rhos_out=rhos, rhofresh_out=rhofresh)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      ! total ice volume
      shmaxn = global_sum(vice, distrb_info, field_loc_center, tarean)
      shmaxs = global_sum(vice, distrb_info, field_loc_center, tareas)

      ! total snow volume
      snwmxn = global_sum(vsno, distrb_info, field_loc_center, tarean)
      snwmxs = global_sum(vsno, distrb_info, field_loc_center, tareas)

      ! north/south ice mass
      totmin = rhoi*shmaxn
      totmis = rhoi*shmaxs

      ! north/south ice+snow mass
      totmn = totmin + rhos*snwmxn
      totms = totmis + rhos*snwmxs

      ! north/south ice+snow energy
      call total_energy (work1)
      toten = global_sum(work1, distrb_info, field_loc_center, tarean)
      totes = global_sum(work1, distrb_info, field_loc_center, tareas)

      if (print_points) then
         do n = 1, npnt
            if (my_task == pmloc(n)) then
               i = piloc(n)
               j = pjloc(n)
               iblk = pbloc(n)

               pdhi(n) = vice(i,j,iblk)
               pdhs(n) = vsno(i,j,iblk)
               pde(n) = work1(i,j,iblk)
            endif
         enddo  ! npnt
      endif                     ! print_points

      if (tr_aero) then
         do n=1,n_aero
            !$OMP PARALLEL DO PRIVATE(iblk,i,j)
            do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               work1(i,j,iblk) = trcr(i,j,nt_aero  +4*(n-1),iblk)*vsno(i,j,iblk) &
                               + trcr(i,j,nt_aero+1+4*(n-1),iblk)*vsno(i,j,iblk) &
                               + trcr(i,j,nt_aero+2+4*(n-1),iblk)*vice(i,j,iblk) &
                               + trcr(i,j,nt_aero+3+4*(n-1),iblk)*vice(i,j,iblk)
            enddo
            enddo
            enddo
            !$OMP END PARALLEL DO
            totaeron(n)= global_sum(work1, distrb_info, field_loc_center, tarean)
            totaeros(n)= global_sum(work1, distrb_info, field_loc_center, tareas)
         enddo
      endif

      if (tr_pond_topo) then
         totpn = c0
         totps = c0
         !$OMP PARALLEL DO PRIVATE(iblk,i,j,n)
         do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            work1(i,j,iblk) = c0
            do n = 1, ncat
               work1(i,j,iblk) = work1(i,j,iblk)  &
                               + aicen(i,j,n,iblk) &
                               * trcrn(i,j,nt_apnd,n,iblk) & 
                               * trcrn(i,j,nt_hpnd,n,iblk)
            enddo
         enddo
         enddo
         enddo
         !$OMP END PARALLEL DO
         totpn = global_sum(work1, distrb_info, field_loc_center, tarean)
         totps = global_sum(work1, distrb_info, field_loc_center, tareas)

         ! north/south ice+snow+pond mass
         totmn = totmn + totpn*rhofresh
         totms = totms + totps*rhofresh
      endif

      end subroutine init_mass_diags

!=======================================================================

! Computes total energy of ice and snow in a grid cell.
!
! authors: E. C. Hunke, LANL

      subroutine total_energy (work)

      use ice_blocks, only: nx_block, ny_block
      use ice_domain, only: nblocks
      use ice_domain_size, only: ncat, nilyr, nslyr, max_blocks
      use ice_grid, only: tmask
      use ice_state, only: vicen, vsnon, trcrn

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(out) :: &
         work      ! total energy

      ! local variables

      integer (kind=int_kind) :: &
        icells                ! number of ocean/ice cells

      integer (kind=int_kind), dimension (nx_block*ny_block) :: &
        indxi, &              ! compressed indices in i/j directions
        indxj

      integer (kind=int_kind) :: &
        i, j, k, n, iblk, ij, &
        nt_qice, nt_qsno

      character(len=*), parameter :: subname = '(total_energy)'

      call icepack_query_tracer_indices(nt_qice_out=nt_qice, nt_qsno_out=nt_qsno)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,n,k,ij,icells,indxi,indxj)
      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

         icells = 0
         do j = 1, ny_block
         do i = 1, nx_block
            if (tmask(i,j,iblk)) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif                  ! tmask
         enddo
         enddo

         work(:,:,iblk) = c0

      !-----------------------------------------------------------------
      ! Aggregate
      !-----------------------------------------------------------------

         do n = 1, ncat
            do k = 1, nilyr
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  work(i,j,iblk) = work(i,j,iblk) &
                                 + trcrn(i,j,nt_qice+k-1,n,iblk) &
                                 * vicen(i,j,n,iblk) / real(nilyr,kind=dbl_kind)
               enddo            ! ij
            enddo               ! k

            do k = 1, nslyr
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  work(i,j,iblk) = work(i,j,iblk) &
                                 + trcrn(i,j,nt_qsno+k-1,n,iblk) &
                                 * vsnon(i,j,n,iblk) / real(nslyr,kind=dbl_kind)
               enddo            ! ij
            enddo               ! k
         enddo                  ! n

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      end subroutine total_energy

!=======================================================================

! Computes bulk salinity of ice and snow in a grid cell.
! author: E. C. Hunke, LANL

      subroutine total_salt (work)

      use ice_blocks, only: nx_block, ny_block
      use ice_domain, only: nblocks
      use ice_domain_size, only: ncat, nilyr, max_blocks
      use ice_grid, only: tmask
      use ice_state, only: vicen, trcrn

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(out) :: &
         work      ! total salt

      ! local variables

      integer (kind=int_kind) :: &
        icells                ! number of ocean/ice cells

      integer (kind=int_kind), dimension (nx_block*ny_block) :: &
        indxi, &              ! compressed indices in i/j directions
        indxj

      integer (kind=int_kind) :: &
        i, j, k, n, iblk, ij, &
        nt_sice

      character(len=*), parameter :: subname = '(total_salt)'

      call icepack_query_tracer_indices(nt_sice_out=nt_sice)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,n,k,ij,icells,indxi,indxj)
      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

         icells = 0
         do j = 1, ny_block
         do i = 1, nx_block
            if (tmask(i,j,iblk)) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif                  ! tmask
         enddo
         enddo

         work(:,:,iblk) = c0

      !-----------------------------------------------------------------
      ! Aggregate
      !-----------------------------------------------------------------

         do n = 1, ncat
            do k = 1, nilyr
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  work(i,j,iblk) = work(i,j,iblk) &
                                 + trcrn(i,j,nt_sice+k-1,n,iblk) &
                                 * vicen(i,j,n,iblk) / real(nilyr,kind=dbl_kind)
               enddo            ! ij
            enddo               ! k
         enddo                  ! n

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      end subroutine total_salt

!=======================================================================

!  Find tasks for diagnostic points.
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL

      subroutine init_diags

      use ice_grid, only: hm, TLAT, TLON
      use ice_blocks, only: block, get_block
      use ice_constants, only: c180, c360, p5
      use ice_domain, only: blocks_ice, distrb_info, nblocks
      use ice_global_reductions, only: global_minval, global_maxval

      real (kind=dbl_kind) :: &
         rad_to_deg, &
         puny, &
         latdis  , & ! latitude distance
         londis  , & ! longitude distance
         totdis  , & ! total distance
         mindis  , & ! minimum distance from desired location
         mindis_g    ! global minimum distance from desired location

      integer (kind=int_kind) :: &
         n           , & ! index for point search
         i,j         , & ! grid indices
         iblk        , & ! block index
         ilo,ihi,jlo,jhi ! beginning and end of physical domain

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(init_diags)'

!tcraig, do this all the time now for print_points_state usage
!      if (print_points) then

         call icepack_query_parameters(puny_out=puny, rad_to_deg_out=rad_to_deg)
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
            file=__FILE__, line=__LINE__)

         if (my_task==master_task) then
            write(nu_diag,*) ' '
            write(nu_diag,*) ' Find indices of diagnostic points '
         endif

         piloc(:) = 0
         pjloc(:) = 0
         pbloc(:) = 0
         pmloc(:) = -999
         plat(:)  = -999._dbl_kind
         plon(:)  = -999._dbl_kind

         ! find minimum distance to diagnostic points on this processor 
         do n = 1, npnt
            if (lonpnt(n) > c180) lonpnt(n) = lonpnt(n) - c360

            iindx = 0
            jindx = 0
            bindx = 0
            mindis = 540.0_dbl_kind !  360. + 180.

            if (abs(latpnt(n)) < c360 .and. abs(lonpnt(n)) < c360) then

            ! MDT, 09/2017: Comment out OpenMP directives since loop is not thread-safe
            !!$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,latdis,londis,totdis)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)         
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi

               do j = jlo, jhi
               do i = ilo, ihi
                  if (hm(i,j,iblk) > p5) then
                     latdis = abs(latpnt(n)-TLAT(i,j,iblk)*rad_to_deg)
                     londis = abs(lonpnt(n)-TLON(i,j,iblk)*rad_to_deg) &
                            * cos(TLAT(i,j,iblk))
                     totdis = sqrt(latdis**2 + londis**2)
                     if (totdis < mindis) then
                        mindis = totdis
                        jindx = j
                        iindx = i
                        bindx = iblk
                     endif      ! totdis < mindis
                  endif         ! hm > p5
               enddo            ! i
               enddo            ! j
            enddo               ! iblk
            !!$OMP END PARALLEL DO

            endif

            ! find global minimum distance to diagnostic points 
            mindis_g = global_minval(mindis, distrb_info)

            ! save indices of minimum-distance grid cell
            if (mindis <= 180.0 .and. abs(mindis_g - mindis) < puny) then
               piloc(n) = iindx
               pjloc(n) = jindx
               pbloc(n) = bindx
               pmloc(n) = my_task
               plat(n) = TLAT(iindx,jindx,bindx)*rad_to_deg
               plon(n) = TLON(iindx,jindx,bindx)*rad_to_deg
            endif

            ! communicate to all processors
            piloc(n) = global_maxval(piloc(n), distrb_info)
            pjloc(n) = global_maxval(pjloc(n), distrb_info)
            pbloc(n) = global_maxval(pbloc(n), distrb_info)
            pmloc(n) = global_maxval(pmloc(n), distrb_info)
            plat(n)  = global_maxval(plat(n), distrb_info)
            plon(n)  = global_maxval(plon(n), distrb_info)

            ! write to log file
            if (my_task==master_task) then
               write(nu_diag,*) ' '
               write(nu_diag,100) n,latpnt(n),lonpnt(n),plat(n),plon(n), &
                    piloc(n), pjloc(n), pbloc(n), pmloc(n)
            endif
 100        format(' found point',i4/ &
               '   lat    lon   TLAT   TLON     i     j   block  task'/ &
                4(f6.1,1x),1x,4(i4,2x) )

         enddo                  ! npnt
!      endif                     ! print_points

      end subroutine init_diags

!=======================================================================

! This routine is useful for debugging.
! Calls to it should be inserted in the form (after thermo, for example)
!      do iblk = 1, nblocks
!      do j=jlo,jhi
!      do i=ilo,ihi
!         plabel = 'post thermo'
!         if (istep1 >= check_step .and. iblk==iblkp .and i==ip &
!             .and. j==jp .and. my_task == mtask) &
!         call print_state(plabel,i,j,iblk)
!      enddo
!      enddo
!      enddo
!
! 'use ice_diagnostics' may need to be inserted also
! author: Elizabeth C. Hunke, LANL

      subroutine print_state(plabel,i,j,iblk)

      use ice_blocks, only: block, get_block
      use ice_domain, only: blocks_ice
      use ice_domain_size, only: ncat, nilyr, nslyr
      use ice_state, only: aice0, aicen, vicen, vsnon, uvel, vvel, trcrn
      use ice_flux, only: uatm, vatm, potT, Tair, Qa, flw, frain, fsnow, &
          fsens, flat, evap, flwout, swvdr, swvdf, swidr, swidf, rhoa, &
          frzmlt, sst, sss, Tf, Tref, Qref, Uref, uocn, vocn, strtltx, strtlty

      character (len=20), intent(in) :: plabel

      integer (kind=int_kind), intent(in) :: & 
          i, j       , & ! horizontal indices
          iblk           ! block index

      ! local variables

      real (kind=dbl_kind) :: &
           eidebug, esdebug, &
           qi, qs, Tsnow, &
           rad_to_deg, puny, rhoi, lfresh, rhos, cp_ice

      integer (kind=int_kind) :: n, k, nt_Tsfc, nt_qice, nt_qsno

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(print_state)'

      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice, &
           nt_qsno_out=nt_qsno)
      call icepack_query_parameters( &
           rad_to_deg_out=rad_to_deg, puny_out=puny, rhoi_out=rhoi, lfresh_out=lfresh, &
           rhos_out=rhos, cp_ice_out=cp_ice)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      this_block = get_block(blocks_ice(iblk),iblk)         

      write(nu_diag,*) plabel
      write(nu_diag,*) 'istep1, my_task, i, j, iblk:', &
                        istep1, my_task, i, j, iblk
      write(nu_diag,*) 'Global i and j:', &
                        this_block%i_glob(i), &
                        this_block%j_glob(j) 
      write(nu_diag,*) ' '
      write(nu_diag,*) 'aice0', aice0(i,j,iblk)
      do n = 1, ncat
         write(nu_diag,*) ' '
         write(nu_diag,*) 'n =',n
         write(nu_diag,*) 'aicen', aicen(i,j,n,iblk)
         write(nu_diag,*) 'vicen', vicen(i,j,n,iblk)
         write(nu_diag,*) 'vsnon', vsnon(i,j,n,iblk)
         if (aicen(i,j,n,iblk) > puny) then
            write(nu_diag,*) 'hin', vicen(i,j,n,iblk)/aicen(i,j,n,iblk)
            write(nu_diag,*) 'hsn', vsnon(i,j,n,iblk)/aicen(i,j,n,iblk)
         endif
         write(nu_diag,*) 'Tsfcn',trcrn(i,j,nt_Tsfc,n,iblk)
         write(nu_diag,*) ' '
      enddo                     ! n

      eidebug = c0
      do n = 1,ncat
         do k = 1,nilyr
            qi = trcrn(i,j,nt_qice+k-1,n,iblk)
            write(nu_diag,*) 'qice, cat ',n,' layer ',k, qi
            eidebug = eidebug + qi
            if (aicen(i,j,n,iblk) > puny) then
               write(nu_diag,*)  'qi/rhoi', qi/rhoi
            endif
         enddo
         write(nu_diag,*) ' '
      enddo
      write(nu_diag,*) 'qice(i,j)',eidebug
      write(nu_diag,*) ' '

      esdebug = c0
      do n = 1,ncat
         if (vsnon(i,j,n,iblk) > puny) then
            do k = 1,nslyr
               qs = trcrn(i,j,nt_qsno+k-1,n,iblk)
               write(nu_diag,*) 'qsnow, cat ',n,' layer ',k, qs
               esdebug = esdebug + qs
               Tsnow = (Lfresh + qs/rhos) / cp_ice
               write(nu_diag,*) 'qs/rhos', qs/rhos
               write(nu_diag,*) 'Tsnow', Tsnow
            enddo
            write(nu_diag,*) ' '
         endif
      enddo
      write(nu_diag,*) 'qsnow(i,j)',esdebug
      write(nu_diag,*) ' '

      write(nu_diag,*) 'uvel(i,j)',uvel(i,j,iblk)
      write(nu_diag,*) 'vvel(i,j)',vvel(i,j,iblk)

      write(nu_diag,*) ' '
      write(nu_diag,*) 'atm states and fluxes'
      write(nu_diag,*) '            uatm    = ',uatm (i,j,iblk)
      write(nu_diag,*) '            vatm    = ',vatm (i,j,iblk)
      write(nu_diag,*) '            potT    = ',potT (i,j,iblk)
      write(nu_diag,*) '            Tair    = ',Tair (i,j,iblk)
      write(nu_diag,*) '            Qa      = ',Qa   (i,j,iblk)
      write(nu_diag,*) '            rhoa    = ',rhoa (i,j,iblk)
      write(nu_diag,*) '            swvdr   = ',swvdr(i,j,iblk)
      write(nu_diag,*) '            swvdf   = ',swvdf(i,j,iblk)
      write(nu_diag,*) '            swidr   = ',swidr(i,j,iblk)
      write(nu_diag,*) '            swidf   = ',swidf(i,j,iblk)
      write(nu_diag,*) '            flw     = ',flw  (i,j,iblk)
      write(nu_diag,*) '            frain   = ',frain(i,j,iblk)
      write(nu_diag,*) '            fsnow   = ',fsnow(i,j,iblk)
      write(nu_diag,*) ' '
      write(nu_diag,*) 'ocn states and fluxes'
      write(nu_diag,*) '            frzmlt  = ',frzmlt (i,j,iblk)
      write(nu_diag,*) '            sst     = ',sst    (i,j,iblk)
      write(nu_diag,*) '            sss     = ',sss    (i,j,iblk)
      write(nu_diag,*) '            Tf      = ',Tf     (i,j,iblk)
      write(nu_diag,*) '            uocn    = ',uocn   (i,j,iblk)
      write(nu_diag,*) '            vocn    = ',vocn   (i,j,iblk)
      write(nu_diag,*) '            strtltx = ',strtltx(i,j,iblk)
      write(nu_diag,*) '            strtlty = ',strtlty(i,j,iblk)
      write(nu_diag,*) ' '
      write(nu_diag,*) 'srf states and fluxes'
      write(nu_diag,*) '            Tref    = ',Tref  (i,j,iblk)
      write(nu_diag,*) '            Qref    = ',Qref  (i,j,iblk)
      write(nu_diag,*) '            Uref    = ',Uref  (i,j,iblk)
      write(nu_diag,*) '            fsens   = ',fsens (i,j,iblk)
      write(nu_diag,*) '            flat    = ',flat  (i,j,iblk)
      write(nu_diag,*) '            evap    = ',evap  (i,j,iblk)
      write(nu_diag,*) '            flwout  = ',flwout(i,j,iblk)
      write(nu_diag,*) ' '

      end subroutine print_state

!=======================================================================

! This routine is useful for debugging.
! Calls can be inserted anywhere and it will print info on print_points points
!      call print_points_state(plabel)
!
! 'use ice_diagnostics' may need to be inserted also

      subroutine print_points_state(plabel,ilabel)

      use ice_blocks, only: block, get_block
      use ice_domain, only: blocks_ice
      use ice_domain_size, only: ncat, nilyr, nslyr
      use ice_state, only: aice0, aicen, vicen, vsnon, uvel, vvel, trcrn
      use ice_flux, only: uatm, vatm, potT, Tair, Qa, flw, frain, fsnow, &
          fsens, flat, evap, flwout, swvdr, swvdf, swidr, swidf, rhoa, &
          frzmlt, sst, sss, Tf, Tref, Qref, Uref, uocn, vocn, strtltx, strtlty

      character (len=*), intent(in),optional :: plabel
      integer          , intent(in),optional :: ilabel

      ! local variables

      real (kind=dbl_kind) :: &
           eidebug, esdebug, &
           qi, qs, &
           puny

      integer (kind=int_kind) :: m, n, k, i, j, iblk, nt_Tsfc, nt_qice, nt_qsno
      character(len=256) :: llabel

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(print_points_state)'
      ! ----------------------

      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice, &
           nt_qsno_out=nt_qsno)
      call icepack_query_parameters(puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do m = 1, npnt
      if (my_task == pmloc(m)) then
         i = piloc(m)
         j = pjloc(m)
         iblk = pbloc(m)
         this_block = get_block(blocks_ice(iblk),iblk)         

         if (present(ilabel)) then
            write(llabel,'(i6,a1,i3,a1)') ilabel,':',m,':'
         else
            write(llabel,'(i3,a1)') m,':'
         endif
         if (present(plabel)) then
            write(llabel,'(a)') 'pps:'//trim(plabel)//':'//trim(llabel)
         else
            write(llabel,'(a)') 'pps:'//trim(llabel)
         endif

         write(nu_diag,*) trim(llabel),'istep1, my_task, i, j, iblk=', &
                           istep1, my_task, i, j, iblk
         write(nu_diag,*) trim(llabel),'Global i and j=', &
                           this_block%i_glob(i), &
                           this_block%j_glob(j) 
         write(nu_diag,*) trim(llabel),'aice0=', aice0(i,j,iblk)

      do n = 1, ncat
         write(nu_diag,*) trim(llabel),'aicen=', n,aicen(i,j,n,iblk)
         write(nu_diag,*) trim(llabel),'vicen=', n,vicen(i,j,n,iblk)
         write(nu_diag,*) trim(llabel),'vsnon=', n,vsnon(i,j,n,iblk)
         if (aicen(i,j,n,iblk) > puny) then
            write(nu_diag,*) trim(llabel),'hin=', n,vicen(i,j,n,iblk)/aicen(i,j,n,iblk)
            write(nu_diag,*) trim(llabel),'hsn=', n,vsnon(i,j,n,iblk)/aicen(i,j,n,iblk)
         endif
         write(nu_diag,*) trim(llabel),'Tsfcn=',n,trcrn(i,j,nt_Tsfc,n,iblk)
      enddo

      eidebug = c0
      do n = 1,ncat
         do k = 1,nilyr
            qi = trcrn(i,j,nt_qice+k-1,n,iblk)
            write(nu_diag,*) trim(llabel),'qice= ',n,k, qi
            eidebug = eidebug + qi
         enddo
      enddo
      write(nu_diag,*) trim(llabel),'qice=',eidebug

      esdebug = c0
      do n = 1,ncat
         if (vsnon(i,j,n,iblk) > puny) then
            do k = 1,nslyr
               qs = trcrn(i,j,nt_qsno+k-1,n,iblk)
               write(nu_diag,*) trim(llabel),'qsnow=',n,k, qs
               esdebug = esdebug + qs
            enddo
         endif
      enddo
      write(nu_diag,*) trim(llabel),'qsnow=',esdebug

      write(nu_diag,*) trim(llabel),'uvel=',uvel(i,j,iblk)
      write(nu_diag,*) trim(llabel),'vvel=',vvel(i,j,iblk)

      write(nu_diag,*) ' '
      write(nu_diag,*) 'atm states and fluxes'
      write(nu_diag,*) '            uatm    = ',uatm (i,j,iblk)
      write(nu_diag,*) '            vatm    = ',vatm (i,j,iblk)
      write(nu_diag,*) '            potT    = ',potT (i,j,iblk)
      write(nu_diag,*) '            Tair    = ',Tair (i,j,iblk)
      write(nu_diag,*) '            Qa      = ',Qa   (i,j,iblk)
      write(nu_diag,*) '            rhoa    = ',rhoa (i,j,iblk)
      write(nu_diag,*) '            swvdr   = ',swvdr(i,j,iblk)
      write(nu_diag,*) '            swvdf   = ',swvdf(i,j,iblk)
      write(nu_diag,*) '            swidr   = ',swidr(i,j,iblk)
      write(nu_diag,*) '            swidf   = ',swidf(i,j,iblk)
      write(nu_diag,*) '            flw     = ',flw  (i,j,iblk)
      write(nu_diag,*) '            frain   = ',frain(i,j,iblk)
      write(nu_diag,*) '            fsnow   = ',fsnow(i,j,iblk)
      write(nu_diag,*) ' '
      write(nu_diag,*) 'ocn states and fluxes'
      write(nu_diag,*) '            frzmlt  = ',frzmlt (i,j,iblk)
      write(nu_diag,*) '            sst     = ',sst    (i,j,iblk)
      write(nu_diag,*) '            sss     = ',sss    (i,j,iblk)
      write(nu_diag,*) '            Tf      = ',Tf     (i,j,iblk)
      write(nu_diag,*) '            uocn    = ',uocn   (i,j,iblk)
      write(nu_diag,*) '            vocn    = ',vocn   (i,j,iblk)
      write(nu_diag,*) '            strtltx = ',strtltx(i,j,iblk)
      write(nu_diag,*) '            strtlty = ',strtlty(i,j,iblk)
      write(nu_diag,*) ' '
      write(nu_diag,*) 'srf states and fluxes'
      write(nu_diag,*) '            Tref    = ',Tref  (i,j,iblk)
      write(nu_diag,*) '            Qref    = ',Qref  (i,j,iblk)
      write(nu_diag,*) '            Uref    = ',Uref  (i,j,iblk)
      write(nu_diag,*) '            fsens   = ',fsens (i,j,iblk)
      write(nu_diag,*) '            flat    = ',flat  (i,j,iblk)
      write(nu_diag,*) '            evap    = ',evap  (i,j,iblk)
      write(nu_diag,*) '            flwout  = ',flwout(i,j,iblk)
      write(nu_diag,*) ' '

      endif   ! my_task
      enddo   ! ncnt

      end subroutine print_points_state

!=======================================================================

! prints error information prior to aborting

      subroutine diagnostic_abort(istop, jstop, iblk, istep1, stop_label)

      use ice_blocks, only: block, get_block
      use ice_communicate, only: my_task
      use ice_domain, only: blocks_ice
      use ice_grid, only: TLAT, TLON
      use ice_state, only: aice

      integer (kind=int_kind), intent(in) :: &
         istop, jstop, & ! indices of grid cell where model aborts
         iblk        , & ! block index
         istep1          ! time step number

      character (char_len), intent(in) :: stop_label

      ! local variables

      real (kind=dbl_kind) :: rad_to_deg

      type (block) :: &
         this_block      ! block information for current block

      character(len=*), parameter :: subname = '(diagnostic_abort)'

      call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      this_block = get_block(blocks_ice(iblk),iblk)         

      write (nu_diag,*) 'istep1, my_task, iblk =', &
                         istep1, my_task, iblk
      write (nu_diag,*) 'Global block:', this_block%block_id
      if (istop > 0 .and. jstop > 0) &
      write (nu_diag,*) 'Global i and j:', &
                          this_block%i_glob(istop), &
                          this_block%j_glob(jstop) 
      write (nu_diag,*) 'Lat, Lon:', &
                         TLAT(istop,jstop,iblk)*rad_to_deg, &
                         TLON(istop,jstop,iblk)*rad_to_deg
      write (nu_diag,*) 'aice:', &
                         aice(istop,jstop,iblk)
      call abort_ice (subname//'ERROR: '//trim(stop_label))

      end subroutine diagnostic_abort

!=======================================================================

      end module ice_diagnostics

!=======================================================================
