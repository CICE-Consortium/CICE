!  SVN:$Id: ice_init_column.F90 1228 2017-05-23 21:33:34Z tcraig $
!=========================================================================
!
! Initialization routines for the column package.
!
! author: Elizabeth C. Hunke, LANL
!
! 2014: Moved subroutines from column package modules

      module ice_init_column

      use ice_kinds_mod
      use ice_constants
      use ice_domain_size, only: ncat, nilyr, nslyr, max_blocks 

      implicit none
      save

      private
      public :: init_thermo_vertical, init_shortwave, &
                init_age, init_FY, init_lvl, &
                init_meltponds_cesm, init_meltponds_lvl, init_meltponds_topo, &
                init_aerosol, init_bgc, init_hbrine, init_zbgc

!=======================================================================

      contains

!=======================================================================
!
! Initialize the vertical profile of ice salinity and melting temperature.
!
! authors: C. M. Bitz, UW
!          William H. Lipscomb, LANL

      subroutine init_thermo_vertical

      use ice_blocks, only: nx_block, ny_block
      use icepack_intfc, only: icepack_init_thermo
      use ice_flux, only: salinz, Tmltz

      integer (kind=int_kind) :: &
         i, j, iblk, &  ! horizontal indices
         k              ! ice layer index

      real (kind=dbl_kind), dimension(nilyr+1) :: &
         sprofile                         ! vertical salinity profile

      !-----------------------------------------------------------------
      ! initialize heat_capacity, l_brine, and salinity profile
      !-----------------------------------------------------------------

      call icepack_init_thermo(nilyr, sprofile)

      !-----------------------------------------------------------------
      ! Prescibe vertical profile of salinity and melting temperature.
      ! Note this profile is only used for BL99 thermodynamics.
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,k)
      do iblk = 1,max_blocks
      do j = 1, ny_block
      do i = 1, nx_block
         do k = 1, nilyr+1
            salinz(i,j,k,iblk) = sprofile(k)
            Tmltz (i,j,k,iblk) = -salinz(i,j,k,iblk)*depressT
         enddo ! k
      enddo    ! i
      enddo    ! j
      enddo    ! iblk
      !$OMP END PARALLEL DO

      end subroutine init_thermo_vertical

!=======================================================================
!
!  Initialize shortwave

      subroutine init_shortwave

      use ice_arrays_column, only: fswpenln, Iswabsn, Sswabsn, albicen, &
          albsnon, alvdrn, alidrn, alvdfn, alidfn, fswsfcn, fswthrun, &
          fswintn, albpndn, apeffn, trcrn_sw, dhsn, ffracn, snowfracn, &
          kaer_tab, waer_tab, gaer_tab, kaer_bc_tab, waer_bc_tab, gaer_bc_tab, bcenh, &
          swgrid, igrid
      use ice_blocks, only: block, get_block, nx_block, ny_block
      use ice_calendar, only: nstreams, istep1, dt, calendar_type, &
          days_per_year, nextsw_cday, yday, sec
      use ice_communicate, only: my_task, master_task
      use ice_diagnostics, only: npnt, print_points, pmloc, piloc, pjloc, &
          diagnostic_abort
      use ice_domain, only: nblocks, blocks_ice
      use ice_domain_size, only: n_aero, n_zaero, ncat, nilyr, nslyr, n_algae, nblyr
      use ice_fileunits, only: nu_diag
      use ice_flux, only: alvdf, alidf, alvdr, alidr, &
                          alvdr_ai, alidr_ai, alvdf_ai, alidf_ai, &
                          swvdr, swvdf, swidr, swidf, scale_factor, snowfrac, &
                          albice, albsno, albpnd, apeff_ai, albcnt, coszen, fsnow
      use ice_grid, only: tlat, tlon, tmask
      use ice_restart_shared, only: restart, runtype
      use ice_state, only: aicen, vicen, vsnon, trcrn

      ! column package includes
      use icepack_intfc, only: icepack_step_radiation, icepack_init_orbit, &
           icepack_clear_warnings, icepack_print_warnings
      use icepack_intfc_shared, only: shortwave, dEdd_algae, modal_aero
      use icepack_intfc_tracers, only: nt_Tsfc, &
          nt_alvl, nt_apnd, nt_hpnd, nt_ipnd, nt_aero, tr_bgc_N, &
          tr_zaero, nlt_chl_sw, nlt_zaero_sw, ntrcr, nbtrcr, nbtrcr_sw, nt_fbri, tr_brine, &
          nt_zaero

      integer (kind=int_kind) :: &
         i, j , k    , & ! horizontal indices
         iblk        , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         n                  ! thickness category index

      real (kind=dbl_kind) :: &
         cszn        , & ! counter for history averaging
         netsw           ! flag for shortwave radiation presence

      type (block) :: &
         this_block      ! block information for current block

      logical (kind=log_kind) :: &
         l_stop       , & ! if true, abort the model
         l_print_point, & ! flag to print designated grid point diagnostics
         debug            ! if true, print diagnostics

      character (char_len) :: stop_label

      integer (kind=int_kind) :: &
         ipoint

      real (kind=dbl_kind), dimension(ncat) :: &
         fbri                 ! brine height to ice thickness

      real(kind= dbl_kind), dimension(ntrcr, ncat) :: &
         ztrcr

      real(kind= dbl_kind), dimension(ntrcr, ncat) :: &
         ztrcr_sw

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,n,ilo,ihi,jlo,jhi,this_block, &
      !$OMP                     cszn,l_print_point,debug,ipoint)
      do iblk=1,nblocks

         ! Initialize
         fswpenln(:,:,:,:,iblk) = c0
         Iswabsn(:,:,:,:,iblk) = c0
         Sswabsn(:,:,:,:,iblk) = c0

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         
         do j = 1, ny_block ! can be jlo, jhi
         do i = 1, nx_block ! can be ilo, ihi

            l_print_point = .false.
            debug = .false.
            if (debug .and. print_points) then
               do ipoint = 1, npnt
                  if (my_task == pmloc(ipoint) .and. &
                       i == piloc(ipoint) .and. &
                       j == pjloc(ipoint)) &
                       l_print_point = .true.
                       write (nu_diag, *) 'my_task = ',my_task
               enddo ! n
            endif

            alvdf(i,j,iblk) = c0
            alidf(i,j,iblk) = c0
            alvdr(i,j,iblk) = c0
            alidr(i,j,iblk) = c0
            alvdr_ai(i,j,iblk) = c0
            alidr_ai(i,j,iblk) = c0
            alvdf_ai(i,j,iblk) = c0
            alidf_ai(i,j,iblk) = c0
            albice(i,j,iblk) = c0
            albsno(i,j,iblk) = c0
            albpnd(i,j,iblk) = c0
            snowfrac(i,j,iblk) = c0
            apeff_ai(i,j,iblk) = c0

            do n = 1, ncat
               alvdrn(i,j,n,iblk) = c0
               alidrn(i,j,n,iblk) = c0
               alvdfn(i,j,n,iblk) = c0
               alidfn(i,j,n,iblk) = c0
               fswsfcn(i,j,n,iblk) = c0
               fswintn(i,j,n,iblk) = c0
               fswthrun(i,j,n,iblk) = c0
            enddo   ! ncat

         enddo
         enddo
         do j = jlo, jhi
         do i = ilo, ihi

            if (trim(shortwave) == 'dEdd') then ! delta Eddington

#ifndef CCSMCOUPLED
               ! initialize orbital parameters
               ! These come from the driver in the coupled model.
               call icepack_clear_warnings()
               call icepack_init_orbit(l_stop, stop_label)
               call icepack_print_warnings(nu_diag)

               if (l_stop) then
                  call diagnostic_abort(i, j, iblk, istep1, stop_label)
               endif
#endif
            endif

         fbri(:) = c0
         ztrcr_sw(:,:) = c0
         do n = 1, ncat
           do k = 1, ntrcr
             ztrcr(k,n) = trcrn(i,j,k,n,iblk)
           enddo
           if (tr_brine)  fbri(n) = trcrn(i,j,nt_fbri,n,iblk)
         enddo

         if (tmask(i,j,iblk)) then
            call icepack_clear_warnings()
            call icepack_step_radiation (dt,         ncat,                    &
                          n_algae,   tr_zaero, nblyr,                     &
                          ntrcr,     nbtrcr,   nbtrcr_sw,                 &
                          nilyr,    nslyr,       n_aero,                  &
                          n_zaero,  dEdd_algae,  nlt_chl_sw,              &
                          nlt_zaero_sw(:),                                &
                          swgrid(:),           igrid(:),                  &
                          fbri(:),                                        &
                          aicen(i,j,:,iblk),     vicen(i,j,:,iblk),       &
                          vsnon(i,j,:,iblk),                              &
                          trcrn(i,j,nt_Tsfc,:,iblk),                      &
                          trcrn(i,j,nt_alvl,:,iblk),                      &
                          trcrn(i,j,nt_apnd,:,iblk),                      &
                          trcrn(i,j,nt_hpnd,:,iblk),                      &
                          trcrn(i,j,nt_ipnd,:,iblk),                      &
                          trcrn(i,j,nt_aero:nt_aero+4*n_aero-1,:,iblk),   &
                          ztrcr_sw,                                       &
                          ztrcr,                                          &
                          TLAT(i,j,iblk),        TLON(i,j,iblk),          &
                          calendar_type,         days_per_year,           &
                          nextsw_cday,           yday,                    &
                          sec,                                            &
                          kaer_tab, waer_tab,                             &
                          gaer_tab,                                       &
                          kaer_bc_tab(:,:),      waer_bc_tab(:,:),        &
                          gaer_bc_tab(:,:),      bcenh(:,:,:),            &
                          modal_aero,                                     &
                          swvdr(i,j,iblk),       swvdf(i,j,iblk),         &
                          swidr(i,j,iblk),       swidf(i,j,iblk),         &
                          coszen(i,j,iblk),      fsnow(i,j,iblk),         &
                          alvdrn(i,j,:,iblk),    alvdfn(i,j,:,iblk),      &
                          alidrn(i,j,:,iblk),    alidfn(i,j,:,iblk),      &
                          fswsfcn(i,j,:,iblk),   fswintn(i,j,:,iblk),     &
                          fswthrun(i,j,:,iblk),  fswpenln(i,j,:,:,iblk),  &
                          Sswabsn(i,j,:,:,iblk), Iswabsn(i,j,:,:,iblk),   &
                          albicen(i,j,:,iblk),   albsnon(i,j,:,iblk),     &
                          albpndn(i,j,:,iblk),   apeffn(i,j,:,iblk),      &
                          snowfracn(i,j,:,iblk),                          &
                          dhsn(i,j,:,iblk),      ffracn(i,j,:,iblk),      &
                          l_print_point,                                  &
                          initonly = .true.)
            call icepack_print_warnings(nu_diag)
         endif
         
      !-----------------------------------------------------------------
      ! Define aerosol tracer on shortwave grid
      !-----------------------------------------------------------------

      if (dEdd_algae .and. (tr_zaero .or. tr_bgc_N)) then
        do n = 1, ncat
           do k = 1, nbtrcr_sw
              trcrn_sw(i,j,k,n,iblk) = ztrcr_sw(k,n)
           enddo
        enddo
      endif

      !-----------------------------------------------------------------
      ! Aggregate albedos 
      !-----------------------------------------------------------------

            do n = 1, ncat
               
               if (aicen(i,j,n,iblk) > puny) then
                  
                  alvdf(i,j,iblk) = alvdf(i,j,iblk) &
                       + alvdfn(i,j,n,iblk)*aicen(i,j,n,iblk)
                  alidf(i,j,iblk) = alidf(i,j,iblk) &
                       + alidfn(i,j,n,iblk)*aicen(i,j,n,iblk)
                  alvdr(i,j,iblk) = alvdr(i,j,iblk) &
                       + alvdrn(i,j,n,iblk)*aicen(i,j,n,iblk)
                  alidr(i,j,iblk) = alidr(i,j,iblk) &
                       + alidrn(i,j,n,iblk)*aicen(i,j,n,iblk)
                  
                  netsw = swvdr(i,j,iblk) + swidr(i,j,iblk) &
                        + swvdf(i,j,iblk) + swidf(i,j,iblk)
                  if (netsw > puny) then ! sun above horizon
                     albice(i,j,iblk) = albice(i,j,iblk) &
                          + albicen(i,j,n,iblk)*aicen(i,j,n,iblk)
                     albsno(i,j,iblk) = albsno(i,j,iblk) &
                          + albsnon(i,j,n,iblk)*aicen(i,j,n,iblk)
                     albpnd(i,j,iblk) = albpnd(i,j,iblk) &
                          + albpndn(i,j,n,iblk)*aicen(i,j,n,iblk)
                  endif
                  
                  apeff_ai(i,j,iblk) = apeff_ai(i,j,iblk) &
                       + apeffn(i,j,n,iblk)*aicen(i,j,n,iblk)
                  snowfrac(i,j,iblk) = snowfrac(i,j,iblk) &
                       + snowfracn(i,j,n,iblk)*aicen(i,j,n,iblk)
               
               endif ! aicen > puny

            enddo  ! ncat

      !----------------------------------------------------------------
      ! Store grid box mean albedos and fluxes before scaling by aice
      !----------------------------------------------------------------

            alvdf_ai  (i,j,iblk) = alvdf  (i,j,iblk)
            alidf_ai  (i,j,iblk) = alidf  (i,j,iblk)
            alvdr_ai  (i,j,iblk) = alvdr  (i,j,iblk)
            alidr_ai  (i,j,iblk) = alidr  (i,j,iblk)
            
            ! for history averaging
!echmod?            cszn = c0
!echmod            if (coszen(i,j,iblk) > puny) cszn = c1
!echmod            do n = 1, nstreams
!echmod               albcnt(i,j,iblk,n) = albcnt(i,j,iblk,n) + cszn
!echmod            enddo
            
      !----------------------------------------------------------------
      ! Save net shortwave for scaling factor in scale_factor
      !----------------------------------------------------------------
            if (runtype == 'initial' .and. .not. restart) then
               scale_factor(i,j,iblk) = &
	 	      swvdr(i,j,iblk)*(c1 - alvdr_ai(i,j,iblk)) &
	 	    + swvdf(i,j,iblk)*(c1 - alvdf_ai(i,j,iblk)) &
 	            + swidr(i,j,iblk)*(c1 - alidr_ai(i,j,iblk)) &
	 	    + swidf(i,j,iblk)*(c1 - alidf_ai(i,j,iblk))
	    endif

         enddo ! i
         enddo ! j
      enddo ! iblk

      end subroutine init_shortwave

!=======================================================================

!  Initialize ice age tracer (call prior to reading restart data)

      subroutine init_age(iage)

      real(kind=dbl_kind), dimension(:,:,:), intent(out) :: iage

      iage(:,:,:) = c0

      end subroutine init_age

!=======================================================================

!  Initialize ice FY tracer (call prior to reading restart data)

      subroutine init_FY(firstyear)

      real(kind=dbl_kind), dimension(:,:,:), intent(out) :: firstyear

      firstyear(:,:,:) = c0

      end subroutine init_FY

!=======================================================================

!  Initialize ice lvl tracers (call prior to reading restart data)

      subroutine init_lvl(alvl, vlvl) 

      use ice_constants, only: c1

      real(kind=dbl_kind), dimension(:,:,:), intent(out) :: &
         alvl , & ! level ice area fraction
         vlvl     ! level ice volume

      alvl(:,:,:) = c1 ! level ice area fraction
      vlvl(:,:,:) = c1 ! level ice volume

      end subroutine init_lvl

!=======================================================================

!  Initialize melt ponds.

      subroutine init_meltponds_cesm(apnd, hpnd)

      real(kind=dbl_kind), dimension(:,:,:), intent(out) :: &
         apnd , & ! melt pond area fraction
         hpnd     ! melt pond depth

      apnd(:,:,:) = c0
      hpnd(:,:,:) = c0

      end subroutine init_meltponds_cesm

!=======================================================================

!  Initialize melt ponds.

      subroutine init_meltponds_lvl(apnd, hpnd, ipnd, dhsn)

      real(kind=dbl_kind), dimension(:,:,:), intent(out) :: &
         apnd , & ! melt pond area fraction
         hpnd , & ! melt pond depth
         ipnd , & ! melt pond refrozen lid thickness
         dhsn     ! depth difference for snow on sea ice and pond ice

      apnd(:,:,:) = c0
      hpnd(:,:,:) = c0
      ipnd(:,:,:) = c0
      dhsn(:,:,:) = c0

      end subroutine init_meltponds_lvl

!=======================================================================

!  Initialize melt ponds.

      subroutine init_meltponds_topo(apnd, hpnd, ipnd)

      real(kind=dbl_kind), dimension(:,:,:), intent(out) :: &
         apnd , & ! melt pond area fraction
         hpnd , & ! melt pond depth
         ipnd     ! melt pond refrozen lid thickness

      apnd(:,:,:) = c0
      hpnd(:,:,:) = c0
      ipnd(:,:,:) = c0
        
      end subroutine init_meltponds_topo

!=======================================================================

!  Initialize ice aerosol tracer (call prior to reading restart data)

      subroutine init_aerosol(aero)

      real(kind=dbl_kind), dimension(:,:,:,:), intent(out) :: &
         aero ! aerosol tracers

      aero(:,:,:,:) = c0

      end subroutine init_aerosol

!=======================================================================

!  Initialize vertical profile for biogeochemistry

      subroutine init_bgc() 

      use ice_arrays_column, only: zfswin, trcrn_sw, &
          ocean_bio_all, ice_bio_net, snow_bio_net, &
          cgrid, igrid, bphi, iDi, bTiz, iki, &
          Rayleigh_criteria, Rayleigh_real
      use ice_blocks, only: block, get_block, nx_block, ny_block
      use ice_calendar, only: dt, istep1
      use ice_communicate, only: my_task
      use ice_diagnostics, only: npnt, print_points, pmloc, piloc, pjloc, &
          diagnostic_abort
      use ice_domain, only: nblocks, blocks_ice
      use ice_domain_size, only: nblyr, nilyr
      use ice_fileunits, only: nu_diag
      use ice_flux, only: sss
      use ice_flux_bgc, only: nit, amm, sil, dmsp, dms, algalN, &
          doc, don, dic, fed, fep, zaeros, hum
      use ice_forcing_bgc, only: init_bgc_data, get_forcing_bgc
      use ice_restart_column, only: restart_zsal, &
          read_restart_bgc, restart_bgc
      use ice_state, only: trcrn, aicen, vicen, vsnon
      use icepack_intfc_shared, only: solve_zsal, &
         max_algae, max_don, max_doc, max_dic, max_aero, max_fe, &
         max_nbtrcr

      ! column package includes
      use icepack_intfc, only: icepack_init_bgc, icepack_init_zsalinity, &
          icepack_init_ocean_conc, icepack_init_OceanConcArray
      use icepack_intfc_tracers, only: nbtrcr, ntrcr, nt_bgc_S, ntrcr_o, &
          nt_sice, nt_fbri

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk       , & ! horizontal indices
         ilo,ihi,jlo,jhi  , & ! beginning and end of physical domain
         k,m              , & ! vertical index
         n                , & ! category index
         ipoint

      logical (kind=log_kind) :: &
         l_print_point, & ! flag to print designated grid point diagnostics
         debug        , & ! prints debugging output if true
         l_stop       , & ! if true, print diagnostics and abort on return
         RayleighC
        
      character (char_len) :: stop_label

      type (block) :: &
         this_block      ! block information for current block

      real(kind=dbl_kind), dimension(ntrcr,ncat) :: &
         trcrn_bgc 
      
      real(kind=dbl_kind), dimension(nilyr,ncat) :: &
         sicen    

      real(kind=dbl_kind) :: &
         RayleighR

      ! Initialize

      l_stop = .false.

      bphi(:,:,:,:,:) = c0   ! initial porosity for no ice 
      iDi (:,:,:,:,:) = c0   ! interface diffusivity
      bTiz(:,:,:,:,:) = c0   ! initial bio grid ice temperature
      iki (:,:,:,:,:) = c0   ! permeability

      ocean_bio_all(:,:,:,:)   = c0
      ice_bio_net  (:,:,:,:)   = c0 ! integrated ice tracer conc (mmol/m^2 or mg/m^2) 
      snow_bio_net (:,:,:,:)   = c0 ! integrated snow tracer conc (mmol/m^2 or mg/m^2)
      zfswin       (:,:,:,:,:) = c0 ! shortwave flux on bio grid
      trcrn_sw     (:,:,:,:,:) = c0 ! tracers active in the shortwave calculation
      trcrn_bgc    (:,:)       = c0
      RayleighR                = c0
      RayleighC                = .false.

      !-----------------------------------------------------------------
      ! zsalinity initialization
      !-----------------------------------------------------------------
      
      if (solve_zsal) then  ! default values 

         !$OMP PARALLEL DO PRIVATE(iblk,i,j,n,ilo,ihi,jlo,jhi,this_block)
         do iblk = 1, nblocks

            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j = jlo, jhi
            do i = ilo, ihi  
               call icepack_init_zsalinity(nblyr, ntrcr_o, RayleighC, &
                      RayleighR, trcrn_bgc, nt_bgc_S, ncat, sss(i,j,iblk))
               if (.not. restart_zsal) then
                  Rayleigh_real    (i,j,iblk) = RayleighR
                  Rayleigh_criteria(i,j,iblk) = RayleighC
                  do n = 1,ncat
                     do k  = 1, nblyr
                        trcrn    (i,j,nt_bgc_S+k-1,        n,iblk) = &
                        trcrn_bgc(    nt_bgc_S+k-1-ntrcr_o,n)
                     enddo
                  enddo
               endif
            enddo      ! i
            enddo      ! j  
         enddo         ! iblk
      endif ! solve_zsal

      if (.not. solve_zsal) restart_zsal = .false.

      !-----------------------------------------------------------------
      ! biogeochemistry initialization
      !-----------------------------------------------------------------

      if (.not. restart_bgc) then       
     
      !-----------------------------------------------------------------
      ! Initial Ocean Values if not coupled to the ocean bgc
      !-----------------------------------------------------------------
         !$OMP PARALLEL DO PRIVATE(iblk,i,j,n,ilo,ihi,jlo,jhi,this_block)
         do iblk = 1, nblocks

            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j = jlo, jhi
            do i = ilo, ihi
               call icepack_init_ocean_conc ( &
                    amm   (i,j,  iblk), dmsp(i,j,  iblk), dms(i,j,  iblk), &
                    algalN(i,j,:,iblk), doc (i,j,:,iblk), dic(i,j,:,iblk), &
                    don   (i,j,:,iblk), fed (i,j,:,iblk), fep(i,j,:,iblk), &
                    hum   (i,j,  iblk), nit (i,j,  iblk), sil(i,j,  iblk), &
                    zaeros(i,j,:,iblk), max_dic, max_don, max_fe, max_aero)
            enddo  ! i
            enddo  ! j

         enddo     ! iblk

         call init_bgc_data(fed(:,:,1,:),fep(:,:,1,:)) ! input dFe from file
         call get_forcing_bgc                          ! defines nit and sil

      endif     ! .not. restart

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,n,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi  

            do n = 1, ncat
            do k = 1, nilyr
               sicen(k,n) = trcrn(i,j,nt_sice+k-1,n,iblk)
            enddo
            do k = ntrcr_o+1, ntrcr
               trcrn_bgc(k-ntrcr_o,n) = trcrn(i,j,k,n,iblk)
            enddo
            enddo

            call icepack_init_OceanConcArray(max_nbtrcr,         &
                                 max_algae, max_don,  max_doc,   &
                                 max_dic,   max_aero, max_fe,    &
                                 nit (i,j,  iblk), amm(i,j,  iblk), sil   (i,j,  iblk), &
                                 dmsp(i,j,  iblk), dms(i,j,  iblk), algalN(i,j,:,iblk), &
                                 doc (i,j,:,iblk), don(i,j,:,iblk), dic   (i,j,:,iblk), &  
                                 fed (i,j,:,iblk), fep(i,j,:,iblk), zaeros(i,j,:,iblk), &
                                 ocean_bio_all(i,j,:,iblk),         hum   (i,j,  iblk))

            if (l_stop) then
               call diagnostic_abort(i, j, iblk, istep1, stop_label)
            endif

         enddo  ! i
         enddo  ! j

      enddo     ! iblk

      if (.not. restart_bgc) then       
     
         !$OMP PARALLEL DO PRIVATE(iblk,i,j,n,ilo,ihi,jlo,jhi,this_block)
         do iblk = 1, nblocks

            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j = jlo, jhi
            do i = ilo, ihi  

            call icepack_init_bgc(dt, ncat, nblyr, nilyr, ntrcr_o, &
               cgrid, igrid, ntrcr, nbtrcr, &
               sicen(:,:), trcrn_bgc(:,:), sss(i,j, iblk), &
               ocean_bio_all(i,j,:,iblk), &
               l_stop, stop_label)
            enddo  ! i
            enddo  ! j
         enddo     ! iblk

      endif ! .not. restart

      !-----------------------------------------------------------------
      ! read restart to complete BGC initialization
      !-----------------------------------------------------------------

      if (restart_zsal .or. restart_bgc) call read_restart_bgc  

      end subroutine init_bgc

!=======================================================================

!  Initialize brine height tracer

      subroutine init_hbrine()

      use ice_arrays_column, only: first_ice, bgrid, igrid, cgrid, &
          icgrid, swgrid
      use ice_domain_size, only: nblyr
      use ice_state, only: trcrn
      use icepack_intfc, only: icepack_init_hbrine
      use icepack_intfc_tracers, only: nt_fbri, tr_brine
      use icepack_intfc_shared, only: phi_snow

      call icepack_init_hbrine(bgrid, igrid, cgrid, icgrid, &
            swgrid, nblyr, nilyr, phi_snow)

      first_ice(:,:,:,:) = .true.            
      if (tr_brine) trcrn(:,:,nt_fbri,:,:) = c1

      end subroutine init_hbrine

!=======================================================================

! Namelist variables, set to default values; may be altered at run time
! 
! author Elizabeth C. Hunke, LANL
!        Nicole Jeffery, LANL

      subroutine init_zbgc

      use ice_broadcast, only: broadcast_scalar
      use ice_communicate, only: my_task, master_task
      use icepack_constants, only: c1, p5, c0, c5, rhos, rhoi, p1
      use ice_domain_size, only: max_ntrcr, nblyr, nilyr, nslyr, &
                           n_algae, n_zaero, n_doc, n_dic, n_don, &
                           n_fed, n_fep, max_nsw, n_bgc
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_nml, nml_filename, get_fileunit, &
                               release_fileunit, nu_diag
      use ice_restart_column, only: restart_bgc, restart_zsal, &
          restart_hbrine
      use ice_state, only: trcr_base, trcr_depend, n_trcr_strata, &
          nt_strata      

      use icepack_intfc, only:icepack_init_bgc_trcr,  icepack_init_zbgc

      use icepack_intfc_tracers, only: tr_brine, &
          ntrcr,         nbtrcr,       nbtrcr_sw,    &
          ntrcr_o, &
          tr_bgc_Nit,    tr_bgc_Am,    tr_bgc_Sil,   &
          tr_bgc_DMS,    tr_bgc_PON,   tr_bgc_S,     &
          tr_bgc_N,      tr_bgc_C,     tr_bgc_chl,   &
          tr_bgc_DON,    tr_bgc_Fe,    tr_zaero,     &
          tr_bgc_hum,    tr_aero,      nt_fbri,      &  
          nt_bgc_Nit,    nt_bgc_Am,    nt_bgc_Sil,   &
          nt_bgc_DMS,    nt_bgc_PON,   nt_bgc_S,     &
          nt_bgc_N,      nt_bgc_C,     nt_bgc_chl,   &
          nt_bgc_DOC,    nt_bgc_DON,   nt_bgc_DIC,   &
          nt_zaero,      nt_bgc_DMSPp, nt_bgc_DMSPd, &
          nt_bgc_Fed,    nt_bgc_Fep,  nt_zbgc_frac, &
          nlt_zaero_sw,  nlt_chl_sw, &
          nlt_bgc_N, nlt_bgc_Nit, nlt_bgc_Am, nlt_bgc_Sil, &
          nlt_bgc_DMS, nlt_bgc_DMSPp, nlt_bgc_DMSPd, nlt_bgc_C, nlt_bgc_chl, &
          nlt_bgc_DIC, nlt_bgc_DOC, nlt_bgc_PON, &
          nlt_bgc_DON, nlt_bgc_Fed, nlt_bgc_Fep, nlt_zaero, &
          nt_bgc_hum,  nlt_bgc_hum, bio_index_o, bio_index
 
      use icepack_intfc_shared, only: ktherm, shortwave, solve_zsal, &
          skl_bgc, z_tracers, scale_bgc, dEdd_algae, solve_zbgc, &
          bgc_data_dir, sil_data_type, nit_data_type, fe_data_type, &
          bgc_flux_type, grid_o, l_sk, grid_o_t, initbio_frac, &
          frazil_scav, grid_oS, l_skS, max_nbtrcr, max_algae, max_aero, &
          max_doc, max_dic, max_don, max_fe, restore_bgc, phi_snow, &
          modal_aero, &
          ratio_Si2N_diatoms , ratio_Si2N_sp      , ratio_Si2N_phaeo   ,  &
          ratio_S2N_diatoms  , ratio_S2N_sp       , ratio_S2N_phaeo    ,  &
          ratio_Fe2C_diatoms , ratio_Fe2C_sp      , ratio_Fe2C_phaeo   ,  &
          ratio_Fe2N_diatoms , ratio_Fe2N_sp      , ratio_Fe2N_phaeo   ,  &
          ratio_Fe2DON       , ratio_Fe2DOC_s     , ratio_Fe2DOC_l     ,  &
          fr_resp            , tau_min            , tau_max            ,  &
          algal_vel          , R_dFe2dust         , dustFe_sol         ,  &
          chlabs_diatoms     , chlabs_sp          , chlabs_phaeo       ,  &
          alpha2max_low_diatoms,alpha2max_low_sp  , alpha2max_low_phaeo,  &
          beta2max_diatoms   , beta2max_sp        , beta2max_phaeo     ,  &
          mu_max_diatoms     , mu_max_sp          , mu_max_phaeo       ,  &
          grow_Tdep_diatoms  , grow_Tdep_sp       , grow_Tdep_phaeo    ,  &
          fr_graze_diatoms   , fr_graze_sp        , fr_graze_phaeo     ,  &
          mort_pre_diatoms   , mort_pre_sp        , mort_pre_phaeo     ,  &
          mort_Tdep_diatoms  , mort_Tdep_sp       , mort_Tdep_phaeo    ,  &
          k_exude_diatoms    , k_exude_sp         , k_exude_phaeo      ,  &
          K_Nit_diatoms      , K_Nit_sp           , K_Nit_phaeo        ,  &
          K_Am_diatoms       , K_Am_sp            , K_Am_phaeo         ,  &
          K_Sil_diatoms      , K_Sil_sp           , K_Sil_phaeo        ,  &
          K_Fe_diatoms       , K_Fe_sp            , K_Fe_phaeo         ,  &
          f_don_protein      , kn_bac_protein     , f_don_Am_protein   ,  &
          f_doc_s            , f_doc_l            , f_exude_s          ,  &
          f_exude_l          , k_bac_s            , k_bac_l            ,  &
          T_max              , fsal               , op_dep_min         ,  &
          fr_graze_s         , fr_graze_e         , fr_mort2min        ,  &
          fr_dFe             , k_nitrif           , t_iron_conv        ,  &
          max_loss           , max_dfe_doc1       , fr_resp_s          ,  &
          y_sk_DMS           , t_sk_conv          , t_sk_ox            ,  &
          algaltype_diatoms  , algaltype_sp       , algaltype_phaeo    ,  &
          nitratetype        , ammoniumtype       , silicatetype       ,  &
          dmspptype          , dmspdtype          , humtype            ,  &
          doctype_s          , doctype_l          , dontype_protein    ,  &
          fedtype_1          , feptype_1          , zaerotype_bc1      ,  &
          zaerotype_bc2      , zaerotype_dust1    , zaerotype_dust2    ,  &
          zaerotype_dust3    , zaerotype_dust4    , ratio_C2N_diatoms  ,  &
          ratio_C2N_sp       , ratio_C2N_phaeo    , ratio_chl2N_diatoms,  & 
          ratio_chl2N_sp     , ratio_chl2N_phaeo  , F_abs_chl_diatoms  ,  &
          F_abs_chl_sp       , F_abs_chl_phaeo    , ratio_C2N_proteins 

      integer (kind=int_kind) :: &
        nml_error, & ! namelist i/o error flag
        k            ! loop index  

      !------------------------------------------------------------
      !        Tracers have mobile and  stationary phases. 
      ! ice growth allows for retention, ice melt facilitates mobility
      ! bgc_tracer_type defines the exchange timescales between these phases
      ! -1 : entirely in the mobile phase, no exchange  (this is the default)
      !  0 : retention time scale is tau_min, release time scale is tau_max
      !  1 : retention time scale is tau_max, release time scale is tau_min
      ! 0.5: retention time scale is tau_min, release time scale is tau_min
      !  2 : retention time scale is tau_max, release time scale is tau_max
      ! tau_min and tau_max are defined in icepack_intfc_shared.f90
      !------------------------------------------------------------

      !-----------------------------------------------------------------
      ! namelist variables
      !-----------------------------------------------------------------

      namelist /zbgc_nml/  &
        tr_brine, restart_hbrine, tr_zaero, modal_aero, skl_bgc, &
        z_tracers, dEdd_algae, solve_zbgc, bgc_flux_type, &
        restore_bgc, restart_bgc, scale_bgc, solve_zsal, restart_zsal, &
        bgc_data_dir, sil_data_type, nit_data_type,  fe_data_type, &
        tr_bgc_Nit, tr_bgc_C, tr_bgc_chl, tr_bgc_Am, tr_bgc_Sil, &
        tr_bgc_DMS, tr_bgc_PON, tr_bgc_hum, tr_bgc_DON, tr_bgc_Fe, &
        grid_o, grid_o_t, l_sk, grid_oS, &   
        l_skS, phi_snow,  initbio_frac, frazil_scav, &
        ratio_Si2N_diatoms , ratio_Si2N_sp      , ratio_Si2N_phaeo   ,  &
        ratio_S2N_diatoms  , ratio_S2N_sp       , ratio_S2N_phaeo    ,  &
        ratio_Fe2C_diatoms , ratio_Fe2C_sp      , ratio_Fe2C_phaeo   ,  &
        ratio_Fe2N_diatoms , ratio_Fe2N_sp      , ratio_Fe2N_phaeo   ,  &
        ratio_Fe2DON       , ratio_Fe2DOC_s     , ratio_Fe2DOC_l     ,  &
        fr_resp            , tau_min            , tau_max            ,  &
        algal_vel          , R_dFe2dust         , dustFe_sol         ,  &
        chlabs_diatoms     , chlabs_sp          , chlabs_phaeo       ,  &
        alpha2max_low_diatoms,alpha2max_low_sp  , alpha2max_low_phaeo,  &
        beta2max_diatoms   , beta2max_sp        , beta2max_phaeo     ,  &
        mu_max_diatoms     , mu_max_sp          , mu_max_phaeo       ,  &
        grow_Tdep_diatoms  , grow_Tdep_sp       , grow_Tdep_phaeo    ,  &
        fr_graze_diatoms   , fr_graze_sp        , fr_graze_phaeo     ,  &
        mort_pre_diatoms   , mort_pre_sp        , mort_pre_phaeo     ,  &
        mort_Tdep_diatoms  , mort_Tdep_sp       , mort_Tdep_phaeo    ,  &
        k_exude_diatoms    , k_exude_sp         , k_exude_phaeo      ,  &
        K_Nit_diatoms      , K_Nit_sp           , K_Nit_phaeo        ,  &
        K_Am_diatoms       , K_Am_sp            , K_Am_phaeo         ,  &
        K_Sil_diatoms      , K_Sil_sp           , K_Sil_phaeo        ,  &
        K_Fe_diatoms       , K_Fe_sp            , K_Fe_phaeo         ,  &
        f_don_protein      , kn_bac_protein     , f_don_Am_protein   ,  &
        f_doc_s            , f_doc_l            , f_exude_s          ,  &
        f_exude_l          , k_bac_s            , k_bac_l            ,  &
        T_max              , fsal               , op_dep_min         ,  &
        fr_graze_s         , fr_graze_e         , fr_mort2min        ,  &
        fr_dFe             , k_nitrif           , t_iron_conv        ,  &
        max_loss           , max_dfe_doc1       , fr_resp_s          ,  &
        y_sk_DMS           , t_sk_conv          , t_sk_ox            ,  &
        algaltype_diatoms  , algaltype_sp       , algaltype_phaeo    ,  &
        nitratetype        , ammoniumtype       , silicatetype       ,  &
        dmspptype          , dmspdtype          , humtype            ,  &
        doctype_s          , doctype_l          , dontype_protein    ,  &
        fedtype_1          , feptype_1          , zaerotype_bc1      ,  &
        zaerotype_bc2      , zaerotype_dust1    , zaerotype_dust2    ,  &
        zaerotype_dust3    , zaerotype_dust4    , ratio_C2N_diatoms  ,  &
        ratio_C2N_sp       , ratio_C2N_phaeo    , ratio_chl2N_diatoms,  & 
        ratio_chl2N_sp     , ratio_chl2N_phaeo  , F_abs_chl_diatoms  ,  &
        F_abs_chl_sp       , F_abs_chl_phaeo      , ratio_C2N_proteins 

      !-----------------------------------------------------------------
      ! default values
      !-----------------------------------------------------------------
      tr_brine        = .false.  ! brine height differs from ice height
      tr_zaero        = .false.  ! z aerosol tracers
      modal_aero      = .false.  ! use modal aerosol treatment of aerosols
      restore_bgc     = .false.  ! restore bgc if true
      solve_zsal      = .false.  ! update salinity tracer profile from solve_S_dt
      bgc_data_dir    = 'unknown_bgc_data_dir'
      sil_data_type   = 'default'
      nit_data_type   = 'default'
      fe_data_type    = 'default'
      restart_bgc     = .false.  ! biogeochemistry restart
      restart_zsal    = .false.  ! salinity restart
      restart_hbrine  = .false.  ! hbrine restart
      scale_bgc       = .false.  ! initial bgc tracers proportional to S  
      skl_bgc         = .false.  ! solve skeletal biochemistry 
      z_tracers       = .false.  ! solve vertically resolved tracers
      dEdd_algae      = .false.  ! dynamic algae contributes to shortwave absorption
                                 ! in delta-Eddington calculation
      solve_zbgc      = .false.  ! turn on z layer biochemistry 
      tr_bgc_PON      = .false.  !---------------------------------------------   
      tr_bgc_Nit      = .false.  ! biogeochemistry (skl or zbgc)
      tr_bgc_C        = .false.  ! if skl_bgc = .true. then skl
      tr_bgc_chl      = .false.  ! if z_tracers = .true. then vertically resolved
      tr_bgc_Sil      = .false.  ! if z_tracers + solve_zbgc = .true. then
      tr_bgc_Am       = .false.  ! vertically resolved with reactions  
      tr_bgc_DMS      = .false.  !------------------------------------------------
      tr_bgc_DON      = .false.  ! 
      tr_bgc_hum      = .false.  !
      tr_bgc_Fe       = .false.  ! 
      tr_bgc_N        = .true.   !

      ! brine height parameter
      phi_snow        = p5       ! snow porosity

      ! skl biology parameters
      bgc_flux_type   = 'Jin2006'! type of ocean-ice poston velocity ('constant')

      ! z biology parameters  
      grid_o          = c5           ! for bottom flux        
      grid_o_t        = c5           ! for top flux        
      l_sk            = 7.0_dbl_kind ! characteristic diffusive scale (m)   
      initbio_frac    = c1           ! fraction of ocean trcr concentration in bio trcrs
      frazil_scav     = c1           ! increase in initial bio tracer from ocean scavenging 
      ratio_Si2N_diatoms = 1.8_dbl_kind    ! algal Si to N (mol/mol)                       
      ratio_Si2N_sp      = c0              ! diatoms, small plankton, phaeocystis
      ratio_Si2N_phaeo   = c0
      ratio_S2N_diatoms  = 0.03_dbl_kind   ! algal S  to N (mol/mol)
      ratio_S2N_sp       = 0.03_dbl_kind 
      ratio_S2N_phaeo    = 0.03_dbl_kind
      ratio_Fe2C_diatoms = 0.0033_dbl_kind ! algal Fe to C  (umol/mol)
      ratio_Fe2C_sp      = 0.0033_dbl_kind
      ratio_Fe2C_phaeo   = p1
      ratio_Fe2N_diatoms = 0.023_dbl_kind  ! algal Fe to N  (umol/mol)
      ratio_Fe2N_sp      = 0.023_dbl_kind
      ratio_Fe2N_phaeo   = 0.7_dbl_kind
      ratio_Fe2DON       = 0.023_dbl_kind  ! Fe to N of DON (nmol/umol)
      ratio_Fe2DOC_s     = p1              ! Fe to C of DOC (nmol/umol) saccharids
      ratio_Fe2DOC_l     = 0.033_dbl_kind  ! Fe to C of DOC (nmol/umol) lipids
      fr_resp            = 0.05_dbl_kind   ! frac of algal growth lost due to respiration      
      tau_min            = 5200.0_dbl_kind ! rapid mobile to stationary exchanges (s)
      tau_max            = 1.73e5_dbl_kind ! long time mobile to stationary exchanges (s)
      algal_vel          = 1.11e-8_dbl_kind! 0.5 cm/d(m/s) Lavoie 2005  1.5 cm/day
      R_dFe2dust         = 0.035_dbl_kind  !  g/g (3.5% content) Tagliabue 2009
      dustFe_sol         = 0.005_dbl_kind  ! solubility fraction
      chlabs_diatoms     = 0.03_dbl_kind   ! chl absorption (1/m/(mg/m^3))
      chlabs_sp          = 0.01_dbl_kind
      chlabs_phaeo       = 0.05_dbl_kind
      alpha2max_low_diatoms = 0.8_dbl_kind ! light limitation (1/(W/m^2))  
      alpha2max_low_sp      = 0.67_dbl_kind
      alpha2max_low_phaeo   = 0.67_dbl_kind
      beta2max_diatoms   = 0.018_dbl_kind  ! light inhibition (1/(W/m^2))  
      beta2max_sp        = 0.0025_dbl_kind
      beta2max_phaeo     = 0.01_dbl_kind
      mu_max_diatoms     = 1.2_dbl_kind    ! maximum growth rate (1/day) 
      mu_max_sp          = 0.851_dbl_kind
      mu_max_phaeo       = 0.851_dbl_kind
      grow_Tdep_diatoms  = 0.06_dbl_kind ! Temperature dependence of growth (1/C)
      grow_Tdep_sp       = 0.06_dbl_kind
      grow_Tdep_phaeo    = 0.06_dbl_kind
      fr_graze_diatoms   = 0.01_dbl_kind ! Fraction grazed
      fr_graze_sp        = p1
      fr_graze_phaeo     = p1
      mort_pre_diatoms   = 0.007_dbl_kind! Mortality (1/day)
      mort_pre_sp        = 0.007_dbl_kind
      mort_pre_phaeo     = 0.007_dbl_kind
      mort_Tdep_diatoms  = 0.03_dbl_kind ! T dependence of mortality (1/C)
      mort_Tdep_sp       = 0.03_dbl_kind
      mort_Tdep_phaeo    = 0.03_dbl_kind
      k_exude_diatoms    = c0            ! algal exudation (1/d)
      k_exude_sp         = c0
      k_exude_phaeo      = c0
      K_Nit_diatoms      = c1            ! nitrate half saturation (mmol/m^3)
      K_Nit_sp           = c1
      K_Nit_phaeo        = c1
      K_Am_diatoms       = 0.3_dbl_kind  ! ammonium half saturation (mmol/m^3)
      K_Am_sp            = 0.3_dbl_kind
      K_Am_phaeo         = 0.3_dbl_kind
      K_Sil_diatoms      = 4.0_dbl_kind  ! silicate half saturation (mmol/m^3)
      K_Sil_sp           = c0
      K_Sil_phaeo        = c0
      K_Fe_diatoms       = c1            ! iron half saturation (nM)
      K_Fe_sp            = 0.2_dbl_kind
      K_Fe_phaeo         = p1
      f_don_protein      = 0.6_dbl_kind  ! fraction of spilled grazing to proteins           
      kn_bac_protein     = 0.03_dbl_kind ! Bacterial degredation of DON (1/d)                
      f_don_Am_protein   = 0.25_dbl_kind ! fraction of remineralized DON to ammonium         
      f_doc_s            = 0.4_dbl_kind  ! fraction of mortality to DOC 
      f_doc_l            = 0.4_dbl_kind
      f_exude_s          = c1            ! fraction of exudation to DOC
      f_exude_l          = c1
      k_bac_s            = 0.03_dbl_kind ! Bacterial degredation of DOC (1/d)
      k_bac_l            = 0.03_dbl_kind
      T_max              = c0            ! maximum temperature (C)
      fsal               = c1            ! Salinity limitation (ppt)
      op_dep_min         = p1            ! Light attenuates for optical depths exceeding min
      fr_graze_s         = p5            ! fraction of grazing spilled or slopped
      fr_graze_e         = p5            ! fraction of assimilation excreted 
      fr_mort2min        = p5            ! fractionation of mortality to Am
      fr_dFe             = 0.3_dbl_kind  ! fraction of remineralized nitrogen
                                         ! (in units of algal iron)
      k_nitrif           = c0            ! nitrification rate (1/day)           
      t_iron_conv        = 3065.0_dbl_kind ! desorption loss pFe to dFe (day)
      max_loss           = 0.9_dbl_kind ! restrict uptake to % of remaining value 
      max_dfe_doc1       = 0.2_dbl_kind ! max ratio of dFe to saccharides in the ice 
                                         !(nM Fe/muM C)    
      fr_resp_s          = 0.75_dbl_kind ! DMSPd fraction of respiration loss as DMSPd
      y_sk_DMS           = p5            ! fraction conversion given high yield
      t_sk_conv          = 3.0_dbl_kind  ! Stefels conversion time (d)
      t_sk_ox            = 10.0_dbl_kind ! DMS oxidation time (d)
      algaltype_diatoms  = c0            ! ------------------
      algaltype_sp       = p5            !
      algaltype_phaeo    = p5            !
      nitratetype        = -c1           ! mobility type between
      ammoniumtype       = c1            ! stationary <-->  mobile
      silicatetype       = -c1           !
      dmspptype          = p5            !
      dmspdtype          = -c1           !
      humtype            = c1            !
      doctype_s          = p5            !
      doctype_l          = p5            !
      dontype_protein    = p5            !
      fedtype_1          = p5            !
      feptype_1          = p5            !
      zaerotype_bc1      = c1            !
      zaerotype_bc2      = c1            !
      zaerotype_dust1    = c1            !
      zaerotype_dust2    = c1            !
      zaerotype_dust3    = c1            !
      zaerotype_dust4    = c1            !--------------------
      ratio_C2N_diatoms  = 7.0_dbl_kind  ! algal C to N ratio (mol/mol)
      ratio_C2N_sp       = 7.0_dbl_kind
      ratio_C2N_phaeo    = 7.0_dbl_kind
      ratio_chl2N_diatoms= 2.1_dbl_kind  ! algal chlorophyll to N ratio (mg/mmol)
      ratio_chl2N_sp     = 1.1_dbl_kind
      ratio_chl2N_phaeo  = 0.84_dbl_kind
      F_abs_chl_diatoms  = 2.0_dbl_kind  ! scales absorbed radiation for dEdd
      F_abs_chl_sp       = 4.0_dbl_kind
      F_abs_chl_phaeo    = 5.0
      ratio_C2N_proteins = 7.0_dbl_kind  ! ratio of C to N in proteins (mol/mol)       

      ! z salinity  parameters
      grid_oS         = c5            ! for bottom flux         
      l_skS           = 7.0_dbl_kind  ! characteristic diffusive scale (m)  

      !-----------------------------------------------------------------
      ! read from input file
      !-----------------------------------------------------------------

      call get_fileunit(nu_nml)

      if (my_task == master_task) then
         open (nu_nml, file=trim(nml_filename), status='old',iostat=nml_error)
         if (nml_error /= 0) then
            nml_error = -1
         else
            nml_error =  1
         endif 

         print*,'Reading zbgc_nml'
         do while (nml_error > 0)
            read(nu_nml, nml=zbgc_nml,iostat=nml_error)
         end do
         if (nml_error == 0) close(nu_nml)
      endif
      call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
         call abort_ice('error reading zbgc namelist')
      endif
      call release_fileunit(nu_nml)

      !-----------------------------------------------------------------
      ! zsalinity and brine
      !-----------------------------------------------------------------
      if (solve_zsal .and. TRZS == 0) then
         write(nu_diag,*) 'WARNING: solve_zsal=T but 0 zsalinity tracers'
         write(nu_diag,*) 'WARNING: setting solve_zsal = F'
         solve_zsal = .false.      
      elseif (solve_zsal .and. nblyr < 1)  then
         write(nu_diag,*) 'WARNING: solve_zsal=T but 0 zsalinity tracers'
         write(nu_diag,*) 'WARNING: setting solve_zsal = F'
         solve_zsal = .false.     
      endif 

      if (solve_zsal .and. ((.not. tr_brine) .or. (ktherm /= 1))) then
         write(nu_diag,*) 'WARNING: solve_zsal needs tr_brine=T and ktherm=1'
         write(nu_diag,*) 'WARNING: setting tr_brine=T and ktherm=1'
         tr_brine = .true.
         ktherm = 1
      endif

      if (tr_brine .and. TRBRI == 0 ) then
         write(nu_diag,*) &
            'WARNING: tr_brine=T but no brine height compiled'
         write(nu_diag,*) &
            'WARNING: setting solve_zsal and tr_brine = F'
         solve_zsal = .false.
         tr_brine  = .false.
      elseif (tr_brine .and. nblyr < 1 ) then
         write(nu_diag,*) &
            'WARNING: tr_brine=T but no biology layers compiled'
         write(nu_diag,*) &
            'WARNING: setting solve_zsal and tr_brine = F'
         solve_zsal = .false.
         tr_brine  = .false.
      endif 

      call broadcast_scalar(solve_zsal,         master_task)  
      call broadcast_scalar(restart_zsal,       master_task)  
      call broadcast_scalar(tr_brine,           master_task)
      call broadcast_scalar(restart_hbrine,     master_task) 

      call broadcast_scalar(phi_snow,           master_task)
      call broadcast_scalar(grid_oS,            master_task)
      call broadcast_scalar(l_skS,              master_task)

      if (my_task == master_task) then
         write(nu_diag,1010) ' tr_brine                  = ', tr_brine
         if (tr_brine) then
         write(nu_diag,1010) ' restart_hbrine            = ', restart_hbrine
         write(nu_diag,1005) ' phi_snow                  = ', phi_snow
         endif
         if (solve_zsal) then
         write(nu_diag,1010) ' solve_zsal                = ', solve_zsal
         write(nu_diag,1010) ' restart_zsal              = ', restart_zsal
         write(nu_diag,1000) ' grid_oS                   = ', grid_oS
         write(nu_diag,1005) ' l_skS                     = ', l_skS
         endif
      endif

      !-----------------------------------------------------------------
      ! biogeochemistry
      !-----------------------------------------------------------------

      if (.not. tr_brine) then
         if (solve_zbgc) then
            write(nu_diag,*) 'WARNING: tr_brine = F and solve_zbgc = T'
            write(nu_diag,*) 'WARNING: setting solve_zbgc = F'
            solve_zbgc = .false.
         endif
         if (tr_zaero) then
            write(nu_diag,*) 'WARNING: tr_brine = F and tr_zaero = T'
            write(nu_diag,*) 'WARNING: setting tr_zaero = F'
            tr_zaero = .false.
         endif
      endif

      if ((skl_bgc .AND. solve_zbgc) .or. (skl_bgc .AND. z_tracers)) &
              call abort_ice('error:skl_bgc &
              and solve_zbgc or z_tracers are both true')

      if (skl_bgc .AND. tr_zaero) then
         write(nu_diag,*) 'WARNING: skl bgc does not use vertical tracers'
         write(nu_diag,*) 'WARNING: setting tr_zaero = F'
         tr_zaero = .false.
      endif

      if (dEdd_algae .AND. trim(shortwave) /= 'dEdd') then 
         write(nu_diag,*) 'WARNING: dEdd_algae = T but shortwave /= dEdd'
         write(nu_diag,*) 'WARNING: setting dEdd_algae = F'
         dEdd_algae = .false.
      endif

      if (dEdd_algae .AND. (.NOT. tr_bgc_N) .AND. (.NOT. tr_zaero)) then 
         write(nu_diag,*) 'WARNING: need tr_bgc_N or tr_zaero for dEdd_algae'
         write(nu_diag,*) 'WARNING: setting dEdd_algae = F'
         dEdd_algae = .false.
      endif

      if (modal_aero .AND. (.NOT. tr_zaero) .AND. (.NOT. tr_aero)) then
         modal_aero = .false.
      endif
         
      if (modal_aero .AND. trim(shortwave) /= 'dEdd') then 
         write(nu_diag,*) 'WARNING: modal_aero = T but shortwave /= dEdd'
         write(nu_diag,*) 'WARNING: setting modal_aero = F'
         modal_aero = .false.
      endif
      if (n_algae > max_algae) call abort_ice('error:number of algal &
            types exceeds max_algae')
      if (n_doc > max_doc) call abort_ice('error:number of doc &
            types exceeds max_doc')
      if (n_dic > max_doc) call abort_ice('error:number of dic &
            types exceeds max_dic')
      if (n_don > max_don) call abort_ice('error:number of don &
            types exceeds max_don')
      if (n_fed  > max_fe ) call abort_ice('error:number of dissolved fe &
            types exceeds max_fe ')
      if (n_fep  > max_fe ) call abort_ice('error:number of particulate fe &
            types exceeds max_fe ')
      if ((TRBGCS == 0 .and. skl_bgc) .or. (TRALG == 0 .and. skl_bgc)) then
         write(nu_diag,*) &
            'WARNING: skl_bgc=T but 0 bgc or algal tracers compiled'
         write(nu_diag,*) &
            'WARNING: setting skl_bgc = F'
         skl_bgc = .false.
      endif

      if ((TRBGCZ == 0 .and. solve_zbgc) .or. (TRALG == 0 .and. solve_zbgc)) then
         write(nu_diag,*) &
            'WARNING: solve_zbgc=T but 0 zbgc or algal tracers compiled'
         write(nu_diag,*) &
            'WARNING: setting solve_zbgc = F'
         solve_zbgc = .false.
      endif

      if (solve_zbgc .and. .not. z_tracers) z_tracers = .true.
      if (skl_bgc .or. solve_zbgc) then
         tr_bgc_N         = .true.   ! minimum NP biogeochemistry
         tr_bgc_Nit       = .true.
      else
         tr_bgc_N         = .false.
         tr_bgc_C         = .false.
         tr_bgc_chl       = .false.
         tr_bgc_Nit       = .false.
         tr_bgc_Am        = .false.
         tr_bgc_Sil       = .false.
         tr_bgc_hum       = .false.
         tr_bgc_DMS       = .false.
         tr_bgc_PON       = .false.
         tr_bgc_DON       = .false.
         tr_bgc_Fe        = .false.
      endif

      call broadcast_scalar(solve_zbgc,         master_task)
      call broadcast_scalar(skl_bgc,            master_task)
      call broadcast_scalar(restart_bgc,        master_task)
      call broadcast_scalar(bgc_flux_type,      master_task)
      call broadcast_scalar(restore_bgc,        master_task)
      call broadcast_scalar(bgc_data_dir,       master_task)
      call broadcast_scalar(sil_data_type,      master_task)
      call broadcast_scalar(nit_data_type,      master_task)
      call broadcast_scalar(fe_data_type,       master_task)
      call broadcast_scalar(tr_bgc_N,           master_task)
      call broadcast_scalar(tr_bgc_C,           master_task)
      call broadcast_scalar(tr_bgc_chl,         master_task)
      call broadcast_scalar(tr_bgc_Nit,         master_task)
      call broadcast_scalar(tr_bgc_Am,          master_task)
      call broadcast_scalar(tr_bgc_Sil,         master_task)
      call broadcast_scalar(tr_bgc_hum,         master_task)
      call broadcast_scalar(tr_bgc_DMS,         master_task) 
      call broadcast_scalar(tr_bgc_PON,         master_task) 
      call broadcast_scalar(tr_bgc_DON,         master_task) 
      call broadcast_scalar(tr_bgc_Fe,          master_task) 

      !-----------------------------------------------------------------
      ! z layer aerosols
      !-----------------------------------------------------------------
      if (tr_zaero .and. .not. z_tracers) z_tracers = .true.

      if (n_zaero > max_aero) call abort_ice('error:number of z aerosols &
            exceeds max_aero')
         
      call broadcast_scalar(z_tracers,          master_task)
      call broadcast_scalar(tr_zaero,           master_task)
      call broadcast_scalar(dEdd_algae,         master_task) 
      call broadcast_scalar(modal_aero,         master_task)
      call broadcast_scalar(grid_o,             master_task)
      call broadcast_scalar(grid_o_t,           master_task)
      call broadcast_scalar(l_sk,               master_task)
      call broadcast_scalar(scale_bgc,          master_task)
      call broadcast_scalar(initbio_frac,       master_task)
      call broadcast_scalar(frazil_scav,        master_task)
      call broadcast_scalar(ratio_Si2N_diatoms, master_task)
      call broadcast_scalar(ratio_Si2N_sp,      master_task)
      call broadcast_scalar(ratio_Si2N_phaeo,   master_task)
      call broadcast_scalar(ratio_S2N_diatoms,  master_task)
      call broadcast_scalar(ratio_S2N_sp,       master_task)
      call broadcast_scalar(ratio_S2N_phaeo,    master_task)
      call broadcast_scalar(ratio_Fe2C_diatoms, master_task)
      call broadcast_scalar(ratio_Fe2C_sp,      master_task)
      call broadcast_scalar(ratio_Fe2C_phaeo,   master_task)
      call broadcast_scalar(ratio_Fe2N_diatoms, master_task)
      call broadcast_scalar(ratio_Fe2N_sp,      master_task)
      call broadcast_scalar(ratio_Fe2N_phaeo,   master_task)
      call broadcast_scalar(ratio_Fe2DON   ,    master_task)
      call broadcast_scalar(ratio_Fe2DOC_s ,    master_task)
      call broadcast_scalar(ratio_Fe2DOC_l ,    master_task)
      call broadcast_scalar(fr_resp     ,       master_task)
      call broadcast_scalar(tau_min  ,          master_task)
      call broadcast_scalar(tau_max  ,          master_task)
      call broadcast_scalar(algal_vel  ,        master_task)
      call broadcast_scalar(R_dFe2dust ,        master_task)
      call broadcast_scalar(dustFe_sol      ,   master_task)
      call broadcast_scalar(chlabs_diatoms ,  master_task)
      call broadcast_scalar(chlabs_sp      ,  master_task)
      call broadcast_scalar(chlabs_phaeo     ,  master_task)
      call broadcast_scalar(alpha2max_low_diatoms ,  master_task) 
      call broadcast_scalar(alpha2max_low_sp      ,  master_task)
      call broadcast_scalar(alpha2max_low_phaeo   ,  master_task)
      call broadcast_scalar(beta2max_diatoms ,  master_task) 
      call broadcast_scalar(beta2max_sp      ,  master_task) 
      call broadcast_scalar(beta2max_phaeo   ,  master_task) 
      call broadcast_scalar(mu_max_diatoms   ,  master_task) 
      call broadcast_scalar(mu_max_sp        ,  master_task) 
      call broadcast_scalar(mu_max_phaeo     ,  master_task) 
      call broadcast_scalar(grow_Tdep_diatoms,  master_task) 
      call broadcast_scalar(grow_Tdep_sp     ,  master_task) 
      call broadcast_scalar(grow_Tdep_phaeo  ,  master_task) 
      call broadcast_scalar(fr_graze_diatoms ,  master_task) 
      call broadcast_scalar(fr_graze_sp      ,  master_task) 
      call broadcast_scalar(fr_graze_phaeo   ,  master_task)  
      call broadcast_scalar(mort_pre_diatoms ,  master_task) 
      call broadcast_scalar(mort_pre_sp      ,  master_task)  
      call broadcast_scalar(mort_pre_phaeo   ,  master_task)  
      call broadcast_scalar(mort_Tdep_diatoms,  master_task) 
      call broadcast_scalar(mort_Tdep_sp     ,  master_task)   
      call broadcast_scalar(mort_Tdep_phaeo  ,  master_task)   
      call broadcast_scalar(k_exude_diatoms  ,  master_task) 
      call broadcast_scalar(k_exude_sp       ,  master_task)   
      call broadcast_scalar(k_exude_phaeo    ,  master_task)   
      call broadcast_scalar(K_Nit_diatoms    ,  master_task) 
      call broadcast_scalar(K_Nit_sp         ,  master_task)
      call broadcast_scalar(K_Nit_phaeo      ,  master_task)
      call broadcast_scalar(K_Am_diatoms     ,  master_task)
      call broadcast_scalar(K_Am_sp          ,  master_task)
      call broadcast_scalar(K_Am_phaeo       ,  master_task)
      call broadcast_scalar(K_Sil_diatoms    ,  master_task)
      call broadcast_scalar(K_Sil_sp         ,  master_task)
      call broadcast_scalar(K_Sil_phaeo      ,  master_task)
      call broadcast_scalar(K_Fe_diatoms     ,  master_task)
      call broadcast_scalar(K_Fe_sp          ,  master_task)
      call broadcast_scalar(K_Fe_phaeo       ,  master_task)
      call broadcast_scalar(f_don_protein    ,  master_task)
      call broadcast_scalar(kn_bac_protein   ,  master_task)
      call broadcast_scalar(f_don_Am_protein ,  master_task)
      call broadcast_scalar(f_doc_s          ,  master_task)
      call broadcast_scalar(f_doc_l          ,  master_task)  
      call broadcast_scalar(f_exude_s        ,  master_task)
      call broadcast_scalar(f_exude_l        ,  master_task)
      call broadcast_scalar(k_bac_s          ,  master_task) 
      call broadcast_scalar(k_bac_l          ,  master_task)
      call broadcast_scalar(T_max            ,  master_task)
      call broadcast_scalar(fsal             ,  master_task)
      call broadcast_scalar(op_dep_min       ,  master_task)
      call broadcast_scalar(fr_graze_s       ,  master_task) 
      call broadcast_scalar(fr_graze_e       ,  master_task) 
      call broadcast_scalar(fr_mort2min      ,  master_task) 
      call broadcast_scalar(fr_dFe           ,  master_task)
      call broadcast_scalar(k_nitrif         ,  master_task)
      call broadcast_scalar(t_iron_conv      ,  master_task)
      call broadcast_scalar(max_loss         ,  master_task)
      call broadcast_scalar(max_dfe_doc1     ,  master_task)
      call broadcast_scalar(fr_resp_s        ,  master_task)
      call broadcast_scalar(y_sk_DMS         ,  master_task)
      call broadcast_scalar(t_sk_conv        ,  master_task) 
      call broadcast_scalar(t_sk_ox          ,  master_task)
      call broadcast_scalar(algaltype_diatoms,  master_task)
      call broadcast_scalar(algaltype_sp       ,  master_task) 
      call broadcast_scalar(algaltype_phaeo    ,  master_task) 
      call broadcast_scalar(nitratetype        ,  master_task)
      call broadcast_scalar(ammoniumtype       ,  master_task)
      call broadcast_scalar(silicatetype       ,  master_task)
      call broadcast_scalar(dmspptype          ,  master_task)  
      call broadcast_scalar(dmspdtype          ,  master_task) 
      call broadcast_scalar(humtype            ,  master_task)
      call broadcast_scalar(doctype_s          ,  master_task) 
      call broadcast_scalar(doctype_l          ,  master_task)
      call broadcast_scalar(dontype_protein    ,  master_task)
      call broadcast_scalar(fedtype_1          ,  master_task)
      call broadcast_scalar(feptype_1          ,  master_task)
      call broadcast_scalar(zaerotype_bc1      ,  master_task)
      call broadcast_scalar(zaerotype_bc2      ,  master_task)
      call broadcast_scalar(zaerotype_dust1    ,  master_task)
      call broadcast_scalar(zaerotype_dust2    ,  master_task)
      call broadcast_scalar(zaerotype_dust3    ,  master_task)
      call broadcast_scalar(zaerotype_dust4    ,  master_task)
      call broadcast_scalar(ratio_C2N_diatoms  ,  master_task)
      call broadcast_scalar(ratio_C2N_sp       ,  master_task)
      call broadcast_scalar(ratio_C2N_phaeo    ,  master_task)
      call broadcast_scalar(ratio_chl2N_diatoms,  master_task)
      call broadcast_scalar(ratio_chl2N_sp     ,  master_task)
      call broadcast_scalar(ratio_chl2N_phaeo  ,  master_task)
      call broadcast_scalar(F_abs_chl_diatoms  ,  master_task)
      call broadcast_scalar(F_abs_chl_sp       ,  master_task)
      call broadcast_scalar(F_abs_chl_phaeo    ,  master_task)
      call broadcast_scalar(ratio_C2N_proteins ,  master_task) 

      if (skl_bgc .and. n_bgc < 2) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'comp_ice must have number of bgc tracers >= 2'
         write (nu_diag,*) 'number of bgc tracers compiled:',n_bgc
         call abort_ice ('init_zbgc error: skl_bgc and n_bgc < 2')
      endif

      if (solve_zbgc .and. n_bgc < 2) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'comp_ice must have number of zbgc tracers >= 2'
         write (nu_diag,*) 'number of bgc tracers compiled:',n_bgc
         call abort_ice ('init_zbgc error: solve_zbgc and n_bgc < 2')
      endif

      if (tr_zaero .and. TRZAERO <  1) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'comp_ice must have number of TRZAERO > 0'
         write (nu_diag,*) 'in order to solve z aerosols:',TRZAERO
         call abort_ice ('init_zbgc error: tr_zaero and tr zaero < 1')
      endif

      !-----------------------------------------------------------------
      ! initialize tracers etc in the column package
      !-----------------------------------------------------------------
      call icepack_init_zbgc (nblyr, nilyr, nslyr, &
                 n_algae, n_zaero, n_doc, n_dic, n_don, n_fed, n_fep, &
                 trcr_base, trcr_depend, n_trcr_strata, nt_strata, nbtrcr_sw, &
                 tr_brine, nt_fbri, ntrcr, nbtrcr, nt_bgc_Nit, nt_bgc_Am, &
                 nt_bgc_Sil, nt_bgc_DMS, nt_bgc_PON, nt_bgc_S, nt_bgc_N, &
                 nt_bgc_C, nt_bgc_chl, nt_bgc_DOC, nt_bgc_DON, nt_bgc_DIC, & 
                 nt_zaero, nt_bgc_DMSPp, nt_bgc_DMSPd, nt_bgc_Fed, nt_bgc_Fep, &
                 nt_zbgc_frac, tr_bgc_Nit, tr_bgc_Am, tr_bgc_Sil, tr_bgc_DMS, &
                 tr_bgc_PON, tr_bgc_S, tr_bgc_N, tr_bgc_C, tr_bgc_chl, &
                 tr_bgc_DON, tr_bgc_Fe, tr_zaero, nlt_zaero_sw, nlt_chl_sw, &
                 nlt_bgc_N, nlt_bgc_Nit, nlt_bgc_Am, nlt_bgc_Sil, &
                 nlt_bgc_DMS, nlt_bgc_DMSPp, nlt_bgc_DMSPd, &
                 nlt_bgc_C, nlt_bgc_chl, nlt_bgc_DIC, nlt_bgc_DOC, &
                 nlt_bgc_PON, nlt_bgc_DON, nlt_bgc_Fed, nlt_bgc_Fep, &
                 nlt_zaero, &
                 nt_bgc_hum, nlt_bgc_hum, tr_bgc_hum, solve_zsal, &
                 skl_bgc, z_tracers, dEdd_algae, solve_zbgc, &
                 frazil_scav, initbio_frac, bio_index_o, bio_index, ntrcr_o, &
                 max_algae, max_doc, max_dic, max_don, max_fe, &
                 ratio_Si2N_diatoms, ratio_Si2N_sp, ratio_Si2N_phaeo, &
                 ratio_S2N_diatoms, ratio_S2N_sp, ratio_S2N_phaeo, &
                 ratio_Fe2C_diatoms, ratio_Fe2C_sp, ratio_Fe2C_phaeo, &
                 ratio_Fe2N_diatoms, ratio_Fe2N_sp, ratio_Fe2N_phaeo, &
                 ratio_Fe2DON, ratio_Fe2DOC_s,  ratio_Fe2DOC_l, & 
                 chlabs_diatoms, chlabs_sp, chlabs_phaeo, &    
                 alpha2max_low_diatoms, alpha2max_low_sp, alpha2max_low_phaeo, &  
                 beta2max_diatoms, beta2max_sp, beta2max_phaeo, &    
                 mu_max_diatoms, mu_max_sp, mu_max_phaeo, &      
                 grow_Tdep_diatoms, grow_Tdep_sp, grow_Tdep_phaeo, &      
                 fr_graze_diatoms, fr_graze_sp, fr_graze_phaeo, &    
                 mort_pre_diatoms, mort_pre_sp, mort_pre_phaeo, &        
                 mort_Tdep_diatoms, mort_Tdep_sp, mort_Tdep_phaeo, &
                 k_exude_diatoms, k_exude_sp, k_exude_phaeo, &   
                 K_Nit_diatoms, K_Nit_sp, K_Nit_phaeo, &     
                 K_Am_diatoms, K_Am_sp, K_Am_phaeo, &     
                 K_Sil_diatoms, K_Sil_sp, K_Sil_phaeo, &     
                 K_Fe_diatoms, K_Fe_sp, K_Fe_phaeo, & 
                 f_don_protein, kn_bac_protein, &   
                 f_don_Am_protein ,f_doc_s, f_doc_l, &
                 f_exude_s, f_exude_l, k_bac_s,  k_bac_l, &
                 algaltype_diatoms, algaltype_sp, algaltype_phaeo, &
                 doctype_s, doctype_l, dontype_protein, &
                 fedtype_1, feptype_1, zaerotype_bc1, zaerotype_bc2, &
                 zaerotype_dust1, zaerotype_dust2, zaerotype_dust3, &
                 zaerotype_dust4, &
                 ratio_C2N_diatoms, ratio_C2N_sp, ratio_C2N_phaeo, &
                 ratio_chl2N_diatoms, ratio_chl2N_sp, ratio_chl2N_phaeo, &
                 F_abs_chl_diatoms, F_abs_chl_sp, F_abs_chl_phaeo, &
                 ratio_C2N_proteins, &
                 nitratetype, ammoniumtype, dmspptype, dmspdtype, &
                 silicatetype, humtype, tau_min, tau_max)

      !-----------------------------------------------------------------
      ! final consistency checks
      !----------------------------------------------------------------- 
      if (nbtrcr > max_nbtrcr) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'nbtrcr > max_nbtrcr'
         write (nu_diag,*) 'nbtrcr, max_nbtrcr:',nbtrcr, max_nbtrcr
         call abort_ice ('init_zbgc error: nbtrcr > max_nbtrcr')
      endif
      if (.NOT. dEdd_algae) nbtrcr_sw = 1

      if (nbtrcr_sw > max_nsw) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'nbtrcr_sw > max_nsw'
         write (nu_diag,*) 'nbtrcr_sw, max_nsw:',nbtrcr_sw, max_nsw
         call abort_ice ('init_zbgc error: nbtrcr_sw > max_nsw')
      endif

      if (ntrcr > max_ntrcr) then
         write(nu_diag,*) 'max_ntrcr < number of namelist tracers'
         write(nu_diag,*) 'max_ntrcr = ',max_ntrcr,' ntrcr = ',ntrcr
         call abort_ice('max_ntrcr < number of namelist tracers')
      endif                               

      !-----------------------------------------------------------------
      ! spew
      !-----------------------------------------------------------------
      if (my_task == master_task) then
      if (skl_bgc) then

         write(nu_diag,1010) ' skl_bgc                   = ', skl_bgc
         write(nu_diag,1030) ' bgc_flux_type             = ', bgc_flux_type
         write(nu_diag,1010) ' restart_bgc               = ', restart_bgc
         write(nu_diag,1010) ' restore_bgc               = ', restore_bgc
         write(nu_diag,*)    ' bgc_data_dir              = ', &
                               trim(bgc_data_dir)
         write(nu_diag,*)    ' sil_data_type             = ', &
                               trim(sil_data_type)
         write(nu_diag,*)    ' nit_data_type             = ', &
                               trim(nit_data_type)
         write(nu_diag,*)    ' fe_data_type              = ', &
                               trim(fe_data_type)
         write(nu_diag,1020) ' number of bio tracers     = ', nbtrcr
         write(nu_diag,1020) ' number of Isw tracers     = ', nbtrcr_sw
         write(nu_diag,1020) ' number of autotrophs      = ', n_algae
         write(nu_diag,1020) ' number of doc          = ', n_doc
         write(nu_diag,1020) ' number of dic          = ', n_dic
         write(nu_diag,1020) ' number of don          = ', n_don
         write(nu_diag,1020) ' number of fed          = ', n_fed
         write(nu_diag,1020) ' number of fep          = ', n_fep
         write(nu_diag,1010) ' tr_bgc_N               = ', tr_bgc_N
         write(nu_diag,1010) ' tr_bgc_C               = ', tr_bgc_C
         write(nu_diag,1010) ' tr_bgc_chl             = ', tr_bgc_chl
         write(nu_diag,1010) ' tr_bgc_Nit             = ', tr_bgc_Nit
         write(nu_diag,1010) ' tr_bgc_Am              = ', tr_bgc_Am
         write(nu_diag,1010) ' tr_bgc_Sil             = ', tr_bgc_Sil
         write(nu_diag,1010) ' tr_bgc_hum             = ', tr_bgc_hum
         write(nu_diag,1010) ' tr_bgc_DMS             = ', tr_bgc_DMS
         write(nu_diag,1010) ' tr_bgc_PON             = ', tr_bgc_PON
         write(nu_diag,1010) ' tr_bgc_DON             = ', tr_bgc_DON
         write(nu_diag,1010) ' tr_bgc_Fe              = ', tr_bgc_Fe 
        
      elseif (z_tracers) then

         write(nu_diag,*)    ' sil_data_type             = ', &
                               trim(sil_data_type)
         write(nu_diag,*)    ' nit_data_type             = ', &
                               trim(nit_data_type)
         write(nu_diag,*)    ' fe_data_type              = ', &
                               trim(fe_data_type)
         write(nu_diag,*)    ' bgc_data_dir              = ', &
                               trim(bgc_data_dir)
         write(nu_diag,1010) ' restart_bgc               = ', restart_bgc
         write(nu_diag,1010) ' dEdd_algae                = ', dEdd_algae  
         write(nu_diag,1010) ' modal_aero                = ', modal_aero  
         write(nu_diag,1010) ' scale_bgc                 = ', scale_bgc
         write(nu_diag,1010) ' solve_zbgc                = ', solve_zbgc
         write(nu_diag,1020) ' number of ztracers        = ', nbtrcr
         write(nu_diag,1020) ' number of Isw tracers     = ', nbtrcr_sw
         write(nu_diag,1020) ' number of autotrophs      = ', n_algae
         write(nu_diag,1020) ' number of doc             = ', n_doc
         write(nu_diag,1020) ' number of dic             = ', n_dic
         write(nu_diag,1020) ' number of fed             = ', n_fed
         write(nu_diag,1020) ' number of fep             = ', n_fep
         write(nu_diag,1020) ' number of aerosols        = ', n_zaero
         write(nu_diag,1010) ' tr_zaero                  = ', tr_zaero
         write(nu_diag,1010) ' tr_bgc_Nit                = ', tr_bgc_Nit
         write(nu_diag,1010) ' tr_bgc_N                  = ', tr_bgc_N
         write(nu_diag,1010) ' tr_bgc_Am                 = ', tr_bgc_Am
         write(nu_diag,1010) ' tr_bgc_C                  = ', tr_bgc_C
         write(nu_diag,1010) ' tr_bgc_Sil                = ', tr_bgc_Sil
         write(nu_diag,1010) ' tr_bgc_hum                = ', tr_bgc_hum
         write(nu_diag,1010) ' tr_bgc_chl                = ', tr_bgc_chl
         write(nu_diag,1010) ' tr_bgc_DMS                = ', tr_bgc_DMS
         write(nu_diag,1010) ' tr_bgc_PON                = ', tr_bgc_PON
         write(nu_diag,1010) ' tr_bgc_DON                = ', tr_bgc_DON
         write(nu_diag,1010) ' tr_bgc_Fe                 = ', tr_bgc_Fe 
         !bio parameters
         write(nu_diag,1000) ' grid_o                    = ', grid_o
         write(nu_diag,1000) ' grid_o_t                  = ', grid_o_t
         write(nu_diag,1005) ' l_sk                      = ', l_sk
         write(nu_diag,1000) ' initbio_frac              = ', initbio_frac
         write(nu_diag,1000) ' frazil_scav               = ', frazil_scav  

      endif  ! skl_bgc or solve_bgc
      endif  ! master_task

 1000    format (a30,2x,f9.2)  ! a30 to align formatted, unformatted statements
 1005    format (a30,2x,f9.6)  ! float
 1010    format (a30,2x,l6)    ! logical
 1020    format (a30,2x,i6)    ! integer
 1030    format (a30,   a8)    ! character

      end subroutine init_zbgc

!=======================================================================

      end module ice_init_column

!=======================================================================

