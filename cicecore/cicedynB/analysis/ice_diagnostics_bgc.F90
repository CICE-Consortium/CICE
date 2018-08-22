!=======================================================================

! Diagnostic information output during run
!
! authors: Elizabeth C. Hunke, LANL
!          Bruce P. Briegleb, NCAR
!
! 2004: Block structure added by William Lipscomb
! 2006: Converted to free source form (F90) by Elizabeth Hunke

      module ice_diagnostics_bgc

      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0, mps_to_cmpdy, c100, p5, c1
      use ice_fileunits, only: nu_diag
      use ice_fileunits, only: flush_fileunit
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_max_algae, icepack_max_aero
      use icepack_intfc, only: icepack_max_doc, icepack_max_don, icepack_max_fe
      use icepack_intfc, only: icepack_query_parameters, icepack_query_tracer_flags
      use icepack_intfc, only: icepack_query_tracer_indices

      implicit none
      private
      public :: hbrine_diags, bgc_diags, zsal_diags

!=======================================================================

      contains

!=======================================================================

!
! Writes diagnostic info (max, min, global sums, etc) to standard out
!
! authors: Elizabeth C. Hunke, LANL
!          Bruce P. Briegleb, NCAR
!          Cecilia M. Bitz, UW
!          Nicole Jeffery, LANL

      subroutine hbrine_diags
              
      use ice_arrays_column, only: darcy_V
      use ice_broadcast, only: broadcast_scalar, broadcast_array
      use ice_diagnostics, only: npnt, print_points, pmloc, piloc, pjloc, pbloc
      use ice_domain_size, only: nilyr
      use ice_state, only: aice, aicen, vicen, vice, trcr, trcrn

      ! local variables

      integer (kind=int_kind) :: &
         i, j, k, n, iblk, &
         nt_sice, nt_fbri, ktherm

      ! fields at diagnostic points
      real (kind=dbl_kind), dimension(npnt) :: &
         phinS, phinS1, pdarcy_V, pfbri

      real (kind=dbl_kind), dimension(npnt,nilyr) :: &
         pSin, pSin1

      character(len=*), parameter :: subname = '(hbrine_diags)'

      call icepack_query_parameters(ktherm_out=ktherm)
      call icepack_query_tracer_indices(nt_sice_out=nt_sice, nt_fbri_out=nt_fbri)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! Dynamic brine height
      !-----------------------------------------------------------------

      if (print_points) then

      !-----------------------------------------------------------------
      ! state of the ice and associated fluxes for 2 defined points
      ! NOTE these are computed for the last timestep only (not avg)
      !-----------------------------------------------------------------

         do n = 1, npnt
            if (my_task == pmloc(n)) then
               i = piloc(n)
               j = pjloc(n)
               iblk = pbloc(n)            
               phinS1(n) = c0             
               phinS(n) = c0            
               pfbri(n) = trcrn(i,j,nt_fbri,1,iblk) 
               pdarcy_V(n) = darcy_V(i,j,1,iblk)
               if (aice(i,j,iblk) > c0) &
                       phinS(n) = trcr(i,j,nt_fbri,iblk)*vice(i,j,iblk)/aice(i,j,iblk)
               if (aicen(i,j,1,iblk)> c0) &
                       phinS1(n) = trcrn(i,j,nt_fbri,1,iblk)*vicen(i,j,1,iblk)/&
                                                aicen(i,j,1,iblk)                    
               do k = 1,nilyr
                  pSin1(n,k) = trcrn(i,j,nt_sice+k-1,1,iblk)
                  pSin(n,k)  = trcr(i,j,nt_sice+k-1,iblk)
               enddo
            endif                 ! my_task = pmloc
           
            call broadcast_array (pSin    (n,:), pmloc(n))   
            call broadcast_array (pSin1   (n,:), pmloc(n))   
            call broadcast_scalar(pfbri   (n),   pmloc(n))  
            call broadcast_scalar(phinS1  (n),   pmloc(n))  
            call broadcast_scalar(phinS   (n),   pmloc(n)) 
            call broadcast_scalar(pdarcy_V(n),   pmloc(n))
         enddo                  ! npnt
      endif                     ! print_points

      !-----------------------------------------------------------------
      ! start spewing
      !-----------------------------------------------------------------

      if (my_task == master_task) then

      call flush_fileunit(nu_diag)

      !-----------------------------------------------------------------
      ! diagnostics for Arctic and Antarctic points
      !-----------------------------------------------------------------

      if (print_points) then
         write(nu_diag,*) '------ hbrine ------'
         write(nu_diag,900) 'hbrine, (m)        = ',phinS(1),phinS(2)
         write(nu_diag,900) 'fbri, cat1 (m)     = ',pfbri(1),pfbri(2)
         write(nu_diag,900) 'hbrine cat1, (m)   = ',phinS1(1),phinS1(2)  
         write(nu_diag,900) 'darcy_V cat1, (m/s)= ',pdarcy_V(1),pdarcy_V(2)   
         if (ktherm == 2) then          
            write(nu_diag,*) '                         '
            write(nu_diag,*) '------ Thermosaline Salinity ------'
            write(nu_diag,803) 'Sice1(1) cat1 S (ppt)','Sice1(2) cat1 S'
            write(nu_diag,*) '---------------------------------------------------'
            write(nu_diag,802) ((pSin1(n,k),n=1,2), k = 1,nilyr)              
            write(nu_diag,*) '                         '
            write(nu_diag,803) 'Sice(1) bulk S (ppt) ','Sice(2) bulk S'
            write(nu_diag,*) '---------------------------------------------------'
            write(nu_diag,802) ((pSin(n,k),n=1,2), k = 1,nilyr)              
            write(nu_diag,*) '                         '
         endif
      endif                   ! print_points
      endif                   ! my_task = master_task 

  802 format (f24.17,2x,f24.17)
  803 format (a25,2x,a25)
  900 format (a25,2x,f24.17,2x,f24.17)

      end subroutine hbrine_diags

!=======================================================================
!
! Writes diagnostic info (max, min, global sums, etc) to standard out
!
! authors: Elizabeth C. Hunke, LANL
!          Nicole Jeffery, LANL

      subroutine bgc_diags

      use ice_arrays_column, only: ocean_bio, zfswin, fbio_atmice, fbio_snoice, &
          Zoo, grow_net, ice_bio_net, trcrn_sw
      use ice_broadcast, only: broadcast_scalar, broadcast_array
      use ice_diagnostics, only: npnt, print_points, pmloc, piloc, pjloc, pbloc
      use ice_domain_size, only: ncat, nblyr, n_algae, n_zaero, &
          n_doc, n_don, n_fed, n_fep, nilyr, nslyr
      use ice_flux_bgc, only: flux_bio, flux_bio_atm
      use ice_state, only: vicen, vice, trcr

      ! local variables

      integer (kind=int_kind) :: &
         i, j, k, n, nn, iblk,kk, klev
      ! fields at diagnostic points
      real (kind=dbl_kind), dimension(npnt) :: &
         pNit_sk, pAm_sk, pSil_sk, phum_sk, &
         pDMSPp_sk, pDMSPd_sk, pDMS_sk, &
         pNit_ac, pAm_ac, pSil_ac, pDMSP_ac, pDMS_ac, &
         pflux_NO, pflux_Am,  phum_ac, &
         pflux_snow_NO, pflux_snow_Am,  &
         pflux_atm_NO, pflux_atm_Am,  pgrow_net, &
         pflux_hum

      logical (kind=log_kind) :: &
         skl_bgc, z_tracers, dEdd_algae

      logical (kind=log_kind) :: &
         tr_bgc_DMS, tr_bgc_PON, &
         tr_bgc_N, tr_bgc_C, tr_bgc_DON, tr_zaero, tr_bgc_hum, &
         tr_bgc_Nit, tr_bgc_Am, tr_bgc_Sil, tr_bgc_Fe

      integer (kind=int_kind) :: &
         nt_fbri, nt_sice, nt_bgc_nit, nt_bgc_am, nt_bgc_sil, &
         nt_bgc_hum, nt_bgc_pon, nt_bgc_dmspp, nt_bgc_dmspd, nt_bgc_dms, &
         nlt_bgc_hum, nlt_bgc_Nit, nlt_bgc_Am, nlt_bgc_Sil, nlt_chl_sw, &
         nlt_bgc_DMSPp, nlt_bgc_DMS
      integer (kind=int_kind), dimension(icepack_max_algae) :: &
         nt_bgc_n, nlt_bgc_N
      integer (kind=int_kind), dimension(icepack_max_doc) :: &
         nt_bgc_doc, nlt_bgc_DOC
      integer (kind=int_kind), dimension(icepack_max_don) :: &
         nt_bgc_don, nlt_bgc_DON 
      integer (kind=int_kind), dimension(icepack_max_aero) :: &
         nt_zaero, nlt_zaero, nlt_zaero_sw
      integer (kind=int_kind), dimension(icepack_max_fe) :: &
         nt_bgc_fed, nt_bgc_fep, nlt_bgc_Fed, nlt_bgc_Fep

      real (kind=dbl_kind), dimension(npnt,icepack_max_algae) :: &
         pN_ac, pN_tot, pN_sk, pflux_N
      real (kind=dbl_kind), dimension(npnt,icepack_max_doc) :: &
         pDOC_ac, pDOC_sk
      real (kind=dbl_kind), dimension(npnt,icepack_max_don) :: &
         pDON_ac, pDON_sk
      real (kind=dbl_kind), dimension(npnt,icepack_max_fe ) :: &
         pFed_ac,  pFed_sk, pFep_ac, pFep_sk 
      real (kind=dbl_kind), dimension(npnt,icepack_max_aero) :: &
        pflux_zaero, pflux_snow_zaero, pflux_atm_zaero, &
        pflux_atm_zaero_s

      ! vertical  fields of category 1 at diagnostic points for bgc layer model
      real (kind=dbl_kind), dimension(npnt,2) :: &
         pNOs, pAms, pPONs, phums
      real (kind=dbl_kind), dimension(npnt,2,icepack_max_algae) :: &
         pNs
      real (kind=dbl_kind), dimension(npnt,2,icepack_max_doc) :: &
         pDOCs
      real (kind=dbl_kind), dimension(npnt,2,icepack_max_don) :: &
         pDONs
      real (kind=dbl_kind), dimension(npnt,2,icepack_max_fe ) :: &
         pFeds, pFeps 
      real (kind=dbl_kind), dimension(npnt,2,icepack_max_aero) :: &
         pzaeros
      real (kind=dbl_kind), dimension(npnt,nblyr+1) :: &
         pNO, pAm, pPON, pzfswin, pZoo, phum
      real (kind=dbl_kind), dimension(npnt,nblyr+1,icepack_max_algae) :: &
         pN
      real (kind=dbl_kind), dimension(npnt,nblyr+1,icepack_max_aero) :: &
         pzaero
      real (kind=dbl_kind), dimension(npnt,nblyr+1,icepack_max_doc) :: &
         pDOC
      real (kind=dbl_kind), dimension(npnt,nblyr+1,icepack_max_don) :: &
         pDON
      real (kind=dbl_kind), dimension(npnt,nblyr+1,icepack_max_fe ) :: &
         pFed, pFep 
      real (kind=dbl_kind), dimension (nblyr+1) :: & 
         zspace
      real (kind=dbl_kind), dimension (npnt,nslyr+nilyr+2) :: & 
         pchlsw
      real (kind=dbl_kind), dimension(npnt,nslyr+nilyr+2,icepack_max_aero) :: &
         pzaerosw
      character(len=*), parameter :: subname = '(bgc_diags)'

      call icepack_query_parameters(skl_bgc_out=skl_bgc, z_tracers_out=z_tracers, dEdd_algae_out=dEdd_algae)
      call icepack_query_tracer_flags( &
         tr_bgc_DMS_out=tr_bgc_DMS, tr_bgc_PON_out=tr_bgc_PON, &
         tr_bgc_N_out=tr_bgc_N, tr_bgc_C_out=tr_bgc_C, &
         tr_bgc_DON_out=tr_bgc_DON, tr_zaero_out=tr_zaero, tr_bgc_hum_out=tr_bgc_hum, &
         tr_bgc_Nit_out=tr_bgc_Nit, tr_bgc_Am_out=tr_bgc_Am, tr_bgc_Sil_out=tr_bgc_Sil, &
         tr_bgc_Fe_out=tr_bgc_Fe)
      call icepack_query_tracer_indices( &
         nt_fbri_out=nt_fbri, nt_sice_out=nt_sice, nt_zaero_out=nt_zaero, &
         nt_bgc_n_out=nt_bgc_n, nt_bgc_doc_out=nt_bgc_doc, nt_bgc_don_out=nt_bgc_don, &
         nt_bgc_fed_out=nt_bgc_fed, nt_bgc_fep_out=nt_bgc_fep, &
         nt_bgc_nit_out=nt_bgc_nit, nt_bgc_am_out=nt_bgc_am, nt_bgc_sil_out=nt_bgc_sil, &
         nt_bgc_hum_out=nt_bgc_hum, nt_bgc_pon_out=nt_bgc_pon, &
         nt_bgc_dmspp_out=nt_bgc_dmspp, nt_bgc_dmspd_out=nt_bgc_dmspd, nt_bgc_dms_out=nt_bgc_dms, &
         nlt_bgc_N_out=nlt_bgc_N, nlt_zaero_out=nlt_zaero, nlt_chl_sw_out=nlt_chl_sw, &
         nlt_zaero_sw_out=nlt_zaero_sw, nlt_bgc_Fed_out=nlt_bgc_Fed, nlt_bgc_Fep_out=nlt_bgc_Fep, &
         nlt_bgc_hum_out=nlt_bgc_hum, nlt_bgc_Nit_out=nlt_bgc_Nit, nlt_bgc_Am_out=nlt_bgc_Am, &
         nlt_bgc_Sil_out=nlt_bgc_Sil, &
         nlt_bgc_DOC_out=nlt_bgc_DOC, nlt_bgc_DON_out=nlt_bgc_DON, nlt_bgc_DMSPp_out=nlt_bgc_DMSPp, &
         nlt_bgc_DMS_out=nlt_bgc_DMS)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      zspace(:) = c1/real(nblyr,kind=dbl_kind)
      zspace(1) = zspace(1)*p5
      zspace(nblyr+1) = zspace(nblyr+1)*p5    

      klev = 1+nilyr+nslyr
      !-----------------------------------------------------------------
      ! biogeochemical state of the ice
      !-----------------------------------------------------------------

      if (print_points) then

      !-----------------------------------------------------------------
      ! state of the ice and associated fluxes for 2 defined points
      ! NOTE these are computed for the last timestep only (not avg)
      !-----------------------------------------------------------------

         do n = 1, npnt
            if (my_task == pmloc(n)) then
               i = piloc(n)
               j = pjloc(n)
               iblk = pbloc(n)
               pAm_ac(n)   = c0
               pSil_ac(n)  = c0
               phum_ac(n)  = c0
               pDMSP_ac(n) = c0
               pDMS_ac(n)  = c0
               pN_ac(n,:)  = c0
               pDOC_ac(n,:)= c0
               pDON_ac(n,:)= c0
               pFed_ac(n,:)= c0
               pFep_ac(n,:)= c0
               pNit_ac(n) = c0
               if (tr_bgc_N) then
                 do k = 1,n_algae
                    pN_ac(n,k)    = ocean_bio(i,j,nlt_bgc_N(k),iblk) 
                 enddo  !n_algae
               endif    !tr_bgc_N
               if (tr_bgc_C) then
                 do k = 1,n_doc
                    pDOC_ac(n,k)    = ocean_bio(i,j,nlt_bgc_DOC(k),iblk) 
                 enddo  !n_algae
               endif    !tr_bgc_N
               if (tr_bgc_DON) then
                 do k = 1,n_don
                    pDON_ac(n,k)    = ocean_bio(i,j,nlt_bgc_DON(k),iblk) 
                 enddo 
               endif
               if (tr_bgc_Fe ) then
                 do k = 1,n_fed 
                    pFed_ac (n,k)   = ocean_bio(i,j,nlt_bgc_Fed (k),iblk) 
                 enddo 
                 do k = 1,n_fep 
                    pFep_ac (n,k)   = ocean_bio(i,j,nlt_bgc_Fep (k),iblk) 
                 enddo 
               endif
               if (tr_bgc_Nit) &
               pNit_ac(n)  = ocean_bio(i,j,nlt_bgc_Nit,iblk)  ! nit(i,j,iblk)
               if (tr_bgc_Am) &
               pAm_ac(n)   = ocean_bio(i,j,nlt_bgc_Am,iblk)   ! amm(i,j,iblk)
               if (tr_bgc_Sil) &
               pSil_ac(n)  = ocean_bio(i,j,nlt_bgc_Sil,iblk)  ! sil(i,j,iblk)
               if (tr_bgc_hum) &
               phum_ac(n)  = ocean_bio(i,j,nlt_bgc_hum,iblk)  ! hum(i,j,iblk)
               if (tr_bgc_DMS) then
               pDMSP_ac(n) = ocean_bio(i,j,nlt_bgc_DMSPp,iblk)! dmsp(i,j,iblk)
               pDMS_ac(n)  = ocean_bio(i,j,nlt_bgc_DMS,iblk)  ! dms(i,j,iblk)
               endif

               ! fluxes in mmol/m^2/d
               ! concentrations are bulk in mmol/m^3
               ! iron is in 10^-3 mmol/m^3  (nM)

               if (skl_bgc) then
                  pNit_sk(n)   = c0
                  pAm_sk(n)    = c0
                  pSil_sk(n)   = c0
                  phum_sk(n)   = c0
                  pDMSPp_sk(n) = c0
                  pDMSPd_sk(n) = c0
                  pDMS_sk(n)   = c0
                  pN_sk(n,:)   = c0
                  pflux_N(n,:) = c0
                  pDOC_sk(n,:) = c0
                  pDON_sk(n,:) = c0
                  pFed_sk(n,:) = c0
                  pFep_sk(n,:) = c0
                  
                  do k = 1,n_algae            
                    pN_sk(n,k)       = trcr    (i,j,nt_bgc_N(k),   iblk)
                    pflux_N(n,k)     = flux_bio(i,j,nlt_bgc_N(k),  iblk)*mps_to_cmpdy/c100 
                  enddo
                  if (tr_bgc_C) then
                     do k = 1,n_doc
                        pDOC_sk(n,k) = trcr    (i,j,nt_bgc_DOC(k),   iblk)
                     enddo
                  endif
                  if (tr_bgc_DON) then
                     do k = 1,n_don
                        pDON_sk(n,k) = trcr    (i,j,nt_bgc_DON(k),   iblk)
                     enddo
                  endif
                  if (tr_bgc_Fe ) then
                     do k = 1,n_fed 
                        pFed_sk (n,k)= trcr    (i,j,nt_bgc_Fed(k),   iblk)
                     enddo
                     do k = 1,n_fep 
                        pFep_sk (n,k)= trcr    (i,j,nt_bgc_Fep(k),   iblk)
                     enddo
                  endif
                  if (tr_bgc_Nit) then
                     pNit_sk(n)  = trcr    (i,j, nt_bgc_Nit, iblk) 
                     pflux_NO(n) = flux_bio(i,j,nlt_bgc_Nit, iblk)*mps_to_cmpdy/c100 
                  endif
                  if (tr_bgc_Am) then
                     pAm_sk(n)   = trcr    (i,j, nt_bgc_Am,  iblk)
                     pflux_Am(n) = flux_bio(i,j,nlt_bgc_Am,  iblk)*mps_to_cmpdy/c100 
                  endif
                  if (tr_bgc_Sil) then
                     pSil_sk(n)  = trcr    (i,j, nt_bgc_Sil, iblk) 
                  endif
                  if (tr_bgc_hum) then
                     phum_sk(n)  = trcr    (i,j, nt_bgc_hum, iblk) 
                     pflux_hum(n)= flux_bio(i,j,nlt_bgc_hum,  iblk)*mps_to_cmpdy/c100 
                  endif
                  if (tr_bgc_DMS) then
                     pDMSPp_sk(n) = trcr   (i,j,nt_bgc_DMSPp,iblk)
                     pDMSPd_sk(n) = trcr   (i,j,nt_bgc_DMSPd,iblk)
                     pDMS_sk  (n) = trcr   (i,j,nt_bgc_DMS,  iblk)
                  endif

               elseif (z_tracers) then   ! zbgc
                  pflux_NO(n) = c0
                  pN_tot(n,:) = c0
                  pflux_Am(n) = c0
                  pflux_hum(n) = c0
                  pflux_atm_Am(n) = c0
                  pflux_snow_Am(n) = c0
                  pflux_N(n,:) = c0
                  pflux_NO(n) = c0
                  pflux_atm_NO(n) = c0
                  pflux_snow_NO(n) = c0
                  pflux_zaero(n,:) = c0
                  pflux_atm_zaero_s(n,:) = c0
                  pflux_atm_zaero(n,:) = c0
                  pflux_snow_zaero(n,:) = c0
                  if (tr_bgc_Nit) then
                    pflux_NO(n)       = flux_bio(i,j,nlt_bgc_Nit,iblk)*mps_to_cmpdy/c100 
                    pflux_atm_NO(n)   = fbio_atmice(i,j,nlt_bgc_Nit,iblk)*mps_to_cmpdy/c100 
                    pflux_snow_NO(n)  = fbio_snoice(i,j,nlt_bgc_Nit,iblk)*mps_to_cmpdy/c100
                  endif
                  if (tr_bgc_Am) then
                    pflux_Am(n)       = flux_bio(i,j,nlt_bgc_Am,iblk)*mps_to_cmpdy/c100 
                    pflux_atm_Am(n)   = fbio_atmice(i,j,nlt_bgc_Am,iblk)*mps_to_cmpdy/c100 
                    pflux_snow_Am(n)  = fbio_snoice(i,j,nlt_bgc_Am,iblk)*mps_to_cmpdy/c100
                  endif 
                  if (tr_bgc_hum) then
                    pflux_hum(n)       = flux_bio(i,j,nlt_bgc_hum,iblk)*mps_to_cmpdy/c100 
                  endif
                  if (tr_bgc_N)  then
                     do k = 1,n_algae
                     pflux_N(n,k)     = flux_bio(i,j,nlt_bgc_N(k),iblk)*mps_to_cmpdy/c100 
                     enddo
                  endif
                  if (tr_zaero)  then
                     do k = 1,n_zaero
                     pflux_zaero(n,k)      = flux_bio(i,j,nlt_zaero(k),iblk)*mps_to_cmpdy/c100 
                     pflux_atm_zaero_s(n,k)= flux_bio_atm(i,j,nlt_zaero(k),iblk)*mps_to_cmpdy/c100 !*aice
                     pflux_atm_zaero(n,k)  = fbio_atmice(i,j,nlt_zaero(k),iblk)*mps_to_cmpdy/c100
                     pflux_snow_zaero(n,k) = fbio_snoice(i,j,nlt_zaero(k),iblk)*mps_to_cmpdy/c100
                     enddo
                  endif
                  do k = 1, nblyr+1
                   pzfswin(n,k) = c0
                   pZoo(n,k) = c0
                   do nn = 1,ncat
                     pzfswin(n,k) = pzfswin(n,k) + zfswin(i,j,k,nn,iblk)*vicen(i,j,nn,iblk)
                     pZoo(n,k) = pZoo(n,k) + Zoo(i,j,k,nn,iblk)*vicen(i,j,nn,iblk)
                   enddo !nn
                   if (vice(i,j,iblk) > c0) then
                     pzfswin(n,k) = pzfswin(n,k)/vice(i,j,iblk)
                     pZoo(n,k)    = pZoo(n,k)/vice(i,j,iblk)
                   endif !vice
                   pAm(n,k) = c0
                   pN(n,k,:) = c0
                   pDOC(n,k,:) = c0
                   pDON(n,k,:) = c0
                   pFed(n,k,:) = c0
                   pFep(n,k,:) = c0
                   pzaero(n,k,:) = c0
                   pPON(n,k) = c0
                   phum(n,k) = c0
                   pNO(n,k) = c0
                   if (tr_bgc_Nit) pNO(n,k) =  trcr(i,j,nt_bgc_Nit+k-1,iblk)                   
                   if (tr_bgc_Am) pAm(n,k) =  trcr(i,j,nt_bgc_Am+k-1,iblk)     
                   if (tr_bgc_N) then
                     do nn = 1, n_algae
                        pN(n,k,nn)   =  trcr(i,j,nt_bgc_N(nn)+k-1,iblk)
                     enddo   
                   endif     
                   if (tr_bgc_C) then
                     do nn = 1, n_doc
                        pDOC(n,k,nn)   =  trcr(i,j,nt_bgc_DOC(nn)+k-1,iblk)
                     enddo   
                   endif   
                   if (tr_bgc_DON) then
                     do nn = 1, n_don
                        pDON(n,k,nn)   =  trcr(i,j,nt_bgc_DON(nn)+k-1,iblk)
                     enddo   
                   endif    
                   if (tr_bgc_Fe)  then
                     do nn = 1, n_fed
                        pFed(n,k,nn)   =  trcr(i,j,nt_bgc_Fed(nn)+k-1,iblk)
                     enddo   
                     do nn = 1, n_fep
                        pFep(n,k,nn)   =  trcr(i,j,nt_bgc_Fep(nn)+k-1,iblk)
                     enddo   
                   endif   
                   if (tr_zaero) then
                     do nn = 1, n_zaero
                        pzaero(n,k,nn)   =  trcr(i,j,nt_zaero(nn)+k-1,iblk)
                     enddo   
                   endif
                   if (tr_bgc_PON) pPON(n,k) =  trcr(i,j,nt_bgc_PON+k-1,iblk)
                   if (tr_bgc_hum) phum(n,k) =  trcr(i,j,nt_bgc_hum+k-1,iblk)
                 enddo  !k
                 if (tr_bgc_N) then
                   do nn = 1,n_algae
                      pN_tot(n,nn) = ice_bio_net(i,j,nlt_bgc_N(nn),iblk)
                   enddo
                   pgrow_net(n) =  grow_net(i,j,iblk)
                 endif !tr_bgc_N
                 do k = 1,2  !snow concentration
                   pAms(n,k) = c0
                   pNs(n,k,:) = c0
                   pDOCs(n,k,:) = c0
                   pDONs(n,k,:) = c0
                   pFeds (n,k,:)= c0
                   pFeps (n,k,:)= c0
                   pzaeros(n,k,:) = c0
                   pPONs(n,k) = c0
                   phums(n,k) = c0
                   pNOs(n,k) = c0
                   if (tr_bgc_Nit) pNOs(n,k)  = trcr(i,j,nt_bgc_Nit+nblyr+k,iblk)  
                   if (tr_bgc_Am) pAms(n,k) = trcr(i,j,nt_bgc_Am+nblyr+k,iblk)
                   if (tr_bgc_N) then
                     do nn = 1, n_algae
                       pNs(n,k,nn) =  trcr(i,j,nt_bgc_N(nn)+nblyr+k,iblk)
                     enddo
                   endif
                   if (tr_bgc_C) then
                     do nn = 1, n_doc
                       pDOCs(n,k,nn) =  trcr(i,j,nt_bgc_DOC(nn)+nblyr+k,iblk)
                     enddo
                   endif
                   if (tr_bgc_DON) then
                     do nn = 1, n_don
                       pDONs(n,k,nn) =  trcr(i,j,nt_bgc_DON(nn)+nblyr+k,iblk)
                     enddo
                   endif
                   if (tr_bgc_Fe ) then
                     do nn = 1, n_fed 
                       pFeds(n,k,nn) =  trcr(i,j,nt_bgc_Fed(nn)+nblyr+k,iblk)
                     enddo
                     do nn = 1, n_fep 
                       pFeps(n,k,nn) =  trcr(i,j,nt_bgc_Fep(nn)+nblyr+k,iblk)
                     enddo
                   endif
                   if (tr_zaero) then
                     do nn = 1, n_zaero
                       pzaeros(n,k,nn) =  trcr(i,j,nt_zaero(nn)+nblyr+k,iblk)
                     enddo
                   endif
                   if (tr_bgc_PON)pPONs(n,k) =trcr(i,j,nt_bgc_PON+nblyr+k,iblk)
                   if (tr_bgc_hum)phums(n,k) =trcr(i,j,nt_bgc_hum+nblyr+k,iblk)
                 enddo   !k 
               endif
               pchlsw(n,:) = c0
               pzaerosw(n,:,:) = c0
               if (dEdd_algae) then
                 do k = 0, klev
                    if (tr_bgc_N) pchlsw(n,k+1) = trcrn_sw(i,j,nlt_chl_sw+k,1,iblk)
                    if (tr_zaero) then
                    do nn = 1,n_zaero
                       pzaerosw(n,k+1,nn) =  trcrn_sw(i,j,nlt_zaero_sw(nn) + k,1,iblk)
                    enddo
                    endif
                  enddo
               endif              ! dEdd_algae            
            endif                 ! my_task = pmloc
            
            call broadcast_scalar   (pNit_ac          (n),     pmloc(n))             
            call broadcast_scalar   (pAm_ac           (n),     pmloc(n))             
            call broadcast_scalar   (pSil_ac          (n),     pmloc(n))             
            call broadcast_scalar   (phum_ac          (n),     pmloc(n))             
            call broadcast_scalar   (pDMSP_ac         (n),     pmloc(n))             
            call broadcast_scalar   (pDMS_ac          (n),     pmloc(n))  
            call broadcast_scalar   (pflux_NO         (n),     pmloc(n))             
            call broadcast_scalar   (pflux_Am         (n),     pmloc(n))            
            call broadcast_scalar   (pflux_hum        (n),     pmloc(n))
            call broadcast_array    (pN_ac            (n,:),   pmloc(n))
            call broadcast_array    (pflux_N          (n,:),   pmloc(n))
            call broadcast_array    (pDOC_ac          (n,:),   pmloc(n))
            call broadcast_array    (pDON_ac          (n,:),   pmloc(n))
            call broadcast_array    (pFed_ac          (n,:),   pmloc(n))
            call broadcast_array    (pFep_ac          (n,:),   pmloc(n))
            call broadcast_array    (pchlsw           (n,:),   pmloc(n)) 
            call broadcast_array    (pzaerosw         (n,:,:), pmloc(n)) 
            if (skl_bgc) then              ! skl_bgc
               call broadcast_array (pN_sk            (n,:),   pmloc(n))
               call broadcast_array (pDOC_sk          (n,:),   pmloc(n))
               call broadcast_array (pDON_sk          (n,:),   pmloc(n))
               call broadcast_array (pFed_sk          (n,:),   pmloc(n))
               call broadcast_array (pFep_sk          (n,:),   pmloc(n))

               call broadcast_scalar(pNit_sk          (n),     pmloc(n))             
               call broadcast_scalar(pAm_sk           (n),     pmloc(n))             
               call broadcast_scalar(pSil_sk          (n),     pmloc(n))             
               call broadcast_scalar(phum_sk          (n),     pmloc(n))             
               call broadcast_scalar(pDMSPp_sk        (n),     pmloc(n))             
               call broadcast_scalar(pDMSPd_sk        (n),     pmloc(n))             
               call broadcast_scalar(pDMS_sk          (n),     pmloc(n))        
            endif   !tr_bgc_sk

            if (z_tracers) then                   !  z_bgc
               call broadcast_array (pN_tot           (n,:),   pmloc(n))    
               call broadcast_array (pflux_zaero      (n,:),   pmloc(n)) 
               call broadcast_array (pflux_atm_zaero_s(n,:),   pmloc(n)) 
               call broadcast_array (pflux_atm_zaero  (n,:),   pmloc(n)) 
               call broadcast_array (pflux_snow_zaero (n,:),   pmloc(n))
               call broadcast_scalar(pflux_atm_NO     (n),     pmloc(n))            
               call broadcast_scalar(pflux_atm_Am     (n),     pmloc(n))    
               call broadcast_scalar(pflux_snow_NO    (n),     pmloc(n))             
               call broadcast_scalar(pflux_snow_Am    (n),     pmloc(n))
               call broadcast_scalar(pgrow_net        (n),     pmloc(n))
               call broadcast_array (pzfswin          (n,:),   pmloc(n))
               call broadcast_array (pZoo             (n,:),   pmloc(n))
               call broadcast_array (pNO              (n,:),   pmloc(n))
               call broadcast_array (pAm              (n,:),   pmloc(n))
               call broadcast_array (pPON             (n,:),   pmloc(n))
               call broadcast_array (phum             (n,:),   pmloc(n))
               call broadcast_array (pN               (n,:,:), pmloc(n))
               call broadcast_array (pDOC             (n,:,:), pmloc(n))
               call broadcast_array (pDON             (n,:,:), pmloc(n))
               call broadcast_array (pFed             (n,:,:), pmloc(n))
               call broadcast_array (pFep             (n,:,:), pmloc(n))
               call broadcast_array (pzaero           (n,:,:), pmloc(n))
               call broadcast_array (pNOs             (n,:),   pmloc(n))
               call broadcast_array (pAms             (n,:),   pmloc(n))
               call broadcast_array (pPONs            (n,:),   pmloc(n))
               call broadcast_array (phums            (n,:),   pmloc(n))
               call broadcast_array (pNs              (n,:,:), pmloc(n))   
               call broadcast_array (pDOCs            (n,:,:), pmloc(n))   
               call broadcast_array (pDONs            (n,:,:), pmloc(n))   
               call broadcast_array (pFeds            (n,:,:), pmloc(n))   
               call broadcast_array (pFeps            (n,:,:), pmloc(n))   
               call broadcast_array (pzaeros          (n,:,:), pmloc(n))   
            endif    ! z_tracers
         enddo                  ! npnt
      endif                     ! print_points

      !-----------------------------------------------------------------
      ! start spewing
      !-----------------------------------------------------------------

      if (my_task == master_task) then

      call flush_fileunit(nu_diag)

      !-----------------------------------------------------------------
      ! diagnostics for Arctic and Antarctic points
      !-----------------------------------------------------------------

     if (print_points) then
       if (z_tracers) then
         write(nu_diag,803) 'zfswin(1) PAR  ','zfswin(2) PAR '
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,802) ((pzfswin(n,k),n=1,2), k = 1,nblyr+1)              
         write(nu_diag,*) '      '          
         write(nu_diag,803) 'Losses: Zoo(1)(mmol/m^3)  ','Zoo(2)'
         write(nu_diag,803) '        Brine Conc.       ',' Brine Conc'
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,802) ((pZoo(n,k),n=1,2), k = 1,nblyr+1)              
         write(nu_diag,*) '      '     
       endif     
       if (tr_bgc_Nit) then
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,*) '    nitrate conc. (mmol/m^3) or flux (mmol/m^2/d)'
         write(nu_diag,900) 'Ocean conc       = ',pNit_ac(1),pNit_ac(2)
         write(nu_diag,900) 'ice-ocean flux   = ',pflux_NO(1),pflux_NO(2)
         if (skl_bgc) then
           write(nu_diag,900) 'Bulk ice conc.   = ',pNit_sk(1),pNit_sk(2)
         elseif (z_tracers) then
           write(nu_diag,900) 'atm-ice flux     = ',pflux_atm_NO(1),pflux_atm_NO(2)
           write(nu_diag,900) 'snow-ice flux    = ',pflux_snow_NO(1),pflux_snow_NO(2)
           write(nu_diag,*) '             snow + ice conc'
           write(nu_diag,803) '    nitrate(1)','   nitrate(2)'
           write(nu_diag,802) ((pNOs(n,k),n=1,2), k = 1,2)              
           write(nu_diag,802) ((pNO(n,k),n=1,2), k = 1,nblyr+1)              
           write(nu_diag,*) '    '        
         endif
      endif
      if (tr_bgc_PON .and. z_tracers) then
           write(nu_diag,*) '---------------------------------------------------'
           write(nu_diag,*) '    PON snow + ice conc. (mmol/m^3)'
           write(nu_diag,803) '    PON(1)','    PON(2)'
           write(nu_diag,802) ((pPONs(n,k),n=1,2), k = 1,2)              
           write(nu_diag,802) ((pPON(n,k),n=1,2), k = 1,nblyr+1)              
           write(nu_diag,*) ' '
      endif
      if (tr_bgc_hum) then
           write(nu_diag,*) '---------------------------------------------------'
           write(nu_diag,*) '    hum snow + ice conc. (mmolC/m^3)'
           write(nu_diag,900) 'Ocean conc       = ',phum_ac(1),phum_ac(2)
           write(nu_diag,900) 'ice-ocean flux   = ',pflux_hum(1),pflux_hum(2)
         if (skl_bgc) then
           write(nu_diag,900) 'Bulk ice conc.   = ',phum_sk(1),phum_sk(2)
         elseif (z_tracers) then
           write(nu_diag,803) '    hum(1)','    hum(2)'
           write(nu_diag,802) ((phums(n,k),n=1,2), k = 1,2)              
           write(nu_diag,802) ((phum(n,k),n=1,2), k = 1,nblyr+1)              
           write(nu_diag,*) ' '
         endif
      endif
      if (tr_bgc_Am) then
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,*) '    ammonium conc. (mmol/m^3) or flux (mmol/m^2/d)'
         write(nu_diag,900) 'Ocean conc       = ',pAm_ac(1),pAm_ac(2)
         write(nu_diag,900) 'ice-ocean flux   = ',pflux_Am(1),pflux_Am(2)
         if (skl_bgc) then
           write(nu_diag,900) 'Bulk ice conc.   = ',pAm_sk(1),pAm_sk(2)
         elseif (z_tracers) then
           write(nu_diag,900) 'atm-ice flux     = ',pflux_atm_Am(1),pflux_atm_Am(2)
           write(nu_diag,900) 'snow-ice flux    = ',pflux_snow_Am(1),pflux_snow_Am(2)
           write(nu_diag,*) '             snow + ice conc.'
           write(nu_diag,803) '  ammonium(1)','  ammonium (2)'
           write(nu_diag,802) ((pAms(n,k),n=1,2), k = 1,2)              
           write(nu_diag,802) ((pAm(n,k),n=1,2), k = 1,nblyr+1)              
           write(nu_diag,*) '       '     
         endif
      endif
      if (tr_bgc_N) then
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,900) 'tot algal growth (1/d) = ',pgrow_net(1),pgrow_net(2)
       do kk = 1,n_algae
         write(nu_diag,*) '  algal conc. (mmol N/m^3) or flux (mmol N/m^2/d)'
         write(nu_diag,1020) '  type:', kk
         write(nu_diag,900) 'Ocean conc           = ',pN_ac(1,kk),pN_ac(2,kk)
         write(nu_diag,900) 'ice-ocean flux       = ',pflux_N(1,kk),pflux_N(2,kk)
         if (skl_bgc) then
           write(nu_diag,900) 'Bulk ice conc.   = ',pN_sk(1,kk),pN_sk(2,kk)
         elseif (z_tracers) then
           write(nu_diag,900) 'Tot ice (mmolN/m^2) = ',pN_tot(1,kk),pN_tot(2,kk)
           write(nu_diag,*) '             snow + ice conc.'
           write(nu_diag,803) '  algal N(1)','  algal N(2) '
           write(nu_diag,802) ((pNs(n,k,kk),n=1,2), k = 1,2)              
           write(nu_diag,802) ((pN(n,k,kk),n=1,2), k = 1,nblyr+1)              
           write(nu_diag,*) '         '   
         endif
       enddo
      endif
      if (tr_bgc_C) then
       do kk = 1,1 !n_doc
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,*) '  DOC conc. (mmol C/m^3)'
         write(nu_diag,1020) '  type:', kk
         write(nu_diag,900)  'Ocean conc       = ',(pDOC_ac(n,kk),n=1,2)
         if (skl_bgc) then
           write(nu_diag,900)'Bulk ice conc.   = ',pDOC_sk(1,kk),pDOC_sk(2,kk)
         elseif (z_tracers) then
           write(nu_diag,*) '             snow + ice conc.'
           write(nu_diag,803) '  DOC(1)','  DOC(2) '
           write(nu_diag,802) ((pDOCs(n,k,kk),n=1,2), k = 1,2)              
           write(nu_diag,802) ((pDOC(n,k,kk),n=1,2), k = 1,nblyr+1)              
           write(nu_diag,*) '      '      
         endif
       enddo
      endif
      if (tr_bgc_DON) then
       do kk = 1,n_don
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,*) '  DON conc. (mmol N/m^3)'
         write(nu_diag,1020) '  type:', kk
         write(nu_diag,900)  'Ocean conc       = ',(pDON_ac(n,kk),n=1,2)
         if (skl_bgc) then
           write(nu_diag,900)'Bulk ice conc.   = ',pDON_sk(1,kk),pDON_sk(2,kk)
         elseif (z_tracers) then
           write(nu_diag,*) '             snow + ice conc.'
           write(nu_diag,803) '  DON(1)','  DON(2) '
           write(nu_diag,802) ((pDONs(n,k,kk),n=1,2), k = 1,2)              
           write(nu_diag,802) ((pDON(n,k,kk),n=1,2), k = 1,nblyr+1)              
           write(nu_diag,*) '      '      
         endif
       enddo
      endif
      if (tr_bgc_Fe ) then
       do kk = 1,n_fed
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,*) ' dFe  conc. (nM)'
         write(nu_diag,1020) '  type:', kk
         write(nu_diag,900)  'Ocean conc       = ',(pFed_ac (n,kk),n=1,2)
         if (skl_bgc) then
           write(nu_diag,900)'Bulk ice conc.   = ',pFed_sk (1,kk),pFed_sk (2,kk)
         elseif (z_tracers) then
           write(nu_diag,*) '             snow + ice conc.'
           write(nu_diag,803) '  Fed (1)','  Fed (2) '
           write(nu_diag,802) ((pFeds (n,k,kk),n=1,2), k = 1,2)              
           write(nu_diag,802) ((pFed (n,k,kk),n=1,2), k = 1,nblyr+1)              
           write(nu_diag,*) '      '      
         endif
       enddo
       do kk = 1,n_fep
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,*) ' pFe  conc. (nM)'
         write(nu_diag,1020) '  type:', kk
         write(nu_diag,900)  'Ocean conc       = ',(pFep_ac (n,kk),n=1,2)
         if (skl_bgc) then
           write(nu_diag,900)'Bulk ice conc.   = ',pFep_sk (1,kk),pFep_sk (2,kk)
         elseif (z_tracers) then
           write(nu_diag,*) '             snow + ice conc.'
           write(nu_diag,803) '  Fep (1)','  Fep (2) '
           write(nu_diag,802) ((pFeps (n,k,kk),n=1,2), k = 1,2)              
           write(nu_diag,802) ((pFep (n,k,kk),n=1,2), k = 1,nblyr+1)              
           write(nu_diag,*) '      '      
         endif
       enddo
      endif
      if (tr_bgc_DMS) then
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,*) '    DMS (mmol/m^3)      '
         write(nu_diag,900) 'Ocean DMSP    = ',pDMSP_ac(1),pDMSP_ac(2)
         write(nu_diag,900) 'Ocean DMS     = ',pDMS_ac(1),pDMS_ac(2)
         if (skl_bgc) then
          write(nu_diag,900) 'Ice DMSPp    = ',pDMSPp_sk(1),pDMSPp_sk(2)
          write(nu_diag,900) 'Ice DMSPd    = ',pDMSPd_sk(1),pDMSPd_sk(2)
          write(nu_diag,900) 'Ice DMS      = ',pDMS_sk(1),pDMS_sk(2)    
         endif
      endif
      if (tr_zaero .and. z_tracers) then
       do kk = 1,n_zaero
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,*) '  aerosol conc. (kg/m^3) or flux (kg/m^2/d)'
         write(nu_diag,1020) '  type: ',kk
         write(nu_diag,900) 'Atm source flux     = ',pflux_atm_zaero_s(1,kk),pflux_atm_zaero_s(2,kk)
         write(nu_diag,900) 'ice-ocean flux*aice = ',pflux_zaero(1,kk),pflux_zaero(2,kk)
         write(nu_diag,900) 'atm-ice flux*aice   = ',pflux_atm_zaero(1,kk),pflux_atm_zaero(2,kk)
         write(nu_diag,900) 'snow-ice flux*aice  = ',pflux_snow_zaero(1,kk),pflux_snow_zaero(2,kk)
         write(nu_diag,*) '             snow + ice conc.'
         write(nu_diag,803) ' aerosol(1)','    aerosol(2) '
         write(nu_diag,802) ((pzaeros(n,k,kk),n=1,2), k = 1,2)              
         write(nu_diag,802) ((pzaero(n,k,kk),n=1,2), k = 1,nblyr+1)   
         write(nu_diag,*) '            '
       enddo
      endif
      if (dEdd_algae) then
          if (tr_zaero) then
          do kk = 1,n_zaero
          write(nu_diag,*) '---------------------------------------------------'
          write(nu_diag,*) '  Cat 1 aerosol conc. (kg/m^3) on delta-Eddington grid  '       
          write(nu_diag,802) ((pzaerosw(n,k,kk),n=1,2), k = 1,klev +1)   
          enddo
          endif
         if (tr_bgc_N) then
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,*) '  Cat 1 chl (mg/m^3) on delta-Eddington grid  '       
         write(nu_diag,802) ((pchlsw(n,k),n=1,2), k = 1,klev +1)   
         endif
      endif
      endif                   ! print_points
      endif                   ! my_task = master_task 

  802 format (f24.17,2x,f24.17)
  803 format (a25,2x,a25)
  900 format (a25,2x,f24.17,2x,f24.17)
  902 format (a25,10x,f6.1,1x,f6.1,9x,f6.1,1x,f6.1)
 1020 format (a30,2x,i6)    ! integer

      end subroutine bgc_diags

!=======================================================================
!
! Writes diagnostic info (max, min, global sums, etc) to standard out
!
! authors: Elizabeth C. Hunke, LANL
!          Bruce P. Briegleb, NCAR
!          Cecilia M. Bitz, UW
!          Nicole Jeffery, LANL

      subroutine zsal_diags

      use ice_arrays_column, only: fzsal, fzsal_g, sice_rho, bTiz, &
          iDi, bphi, dhbr_top, dhbr_bot, darcy_V
      use ice_broadcast, only: broadcast_scalar, broadcast_array
      use ice_diagnostics, only: npnt, print_points, pmloc, piloc, pjloc, &
          pbloc
      use ice_domain_size, only: nblyr, ncat, nilyr
      use ice_state, only: aicen, aice, vice, trcr, trcrn, vicen, vsno

      ! local variables

      integer (kind=int_kind) :: &
         i, j, k, n, nn, iblk

      ! fields at diagnostic points
      real (kind=dbl_kind), dimension(npnt) :: &
         phinS, phinS1,&
         phbrn,pdh_top1,pdh_bot1, psice_rho, pfzsal, & 
         pfzsal_g, pdarcy_V1 

      ! vertical  fields of category 1 at diagnostic points for bgc layer model
      real (kind=dbl_kind), dimension(npnt,nblyr+2) :: &
         pphin, pphin1
      real (kind=dbl_kind), dimension(npnt,nblyr) :: &
         pSin, pSice, pSin1

      real (kind=dbl_kind), dimension(npnt,nblyr+1) :: &
         pbTiz, piDin

      real (kind=dbl_kind) :: &
         rhosi, rhow, rhos

      logical (kind=log_kind) :: tr_brine

      integer (kind=int_kind) :: nt_fbri, nt_bgc_S, nt_sice
      character(len=*), parameter :: subname = '(zsal_diags)'

      call icepack_query_parameters(rhosi_out=rhosi, rhow_out=rhow, rhos_out=rhos)
      call icepack_query_tracer_flags(tr_brine_out=tr_brine)
      call icepack_query_tracer_indices(nt_fbri_out=nt_fbri, nt_bgc_S_out=nt_bgc_S, &
           nt_sice_out=nt_sice)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! salinity and microstructure  of the ice
      !-----------------------------------------------------------------

      if (print_points) then

      !-----------------------------------------------------------------
      ! state of the ice and associated fluxes for 2 defined points
      ! NOTE these are computed for the last timestep only (not avg)
      !-----------------------------------------------------------------

         do n = 1, npnt
            if (my_task == pmloc(n)) then
               i = piloc(n)
               j = pjloc(n)
               iblk = pbloc(n)

               pfzsal(n) = fzsal(i,j,iblk)   
               pfzsal_g(n) = fzsal_g(i,j,iblk)            
               phinS(n) = c0            
               phinS1(n) = c0 
               phbrn(n) = c0
               psice_rho(n) = c0
               pdh_top1(n) = c0
               pdh_bot1(n) = c0
               pdarcy_V1(n) = c0
               do nn = 1,ncat
                  psice_rho(n) = psice_rho(n) + sice_rho(i,j,nn,iblk)*aicen(i,j,nn,iblk)
               enddo
               if (aice(i,j,iblk) > c0) &
                  psice_rho(n) = psice_rho(n)/aice(i,j,iblk)
               if (tr_brine .and. aice(i,j,iblk) > c0) &
                  phinS(n) = trcr(i,j,nt_fbri,iblk)*vice(i,j,iblk)/aice(i,j,iblk)

               if (aicen(i,j,1,iblk)> c0) then
                  if (tr_brine) phinS1(n) = trcrn(i,j,nt_fbri,1,iblk) &
                                          * vicen(i,j,1,iblk)/aicen(i,j,1,iblk)
                  pdh_top1(n) = dhbr_top(i,j,1,iblk)
                  pdh_bot1(n) = dhbr_bot(i,j,1,iblk)
                  pdarcy_V1(n) = darcy_V(i,j,1,iblk)
               endif  
               if (tr_brine .AND. aice(i,j,iblk) > c0) &
                  phbrn(n) = (c1 - rhosi/rhow)*vice(i,j,iblk)/aice(i,j,iblk) &
                                 - rhos/rhow  *vsno(i,j,iblk)/aice(i,j,iblk)
               do k = 1, nblyr+1
                  pbTiz(n,k) = c0
                  piDin(n,k) = c0
                  do nn = 1,ncat
                     pbTiz(n,k) = pbTiz(n,k) + bTiz(i,j,k,nn,iblk)*vicen(i,j,nn,iblk)
                     piDin(n,k) = piDin(n,k) +  iDi(i,j,k,nn,iblk)*vicen(i,j,nn,iblk)
                  enddo
                  if (vice(i,j,iblk) > c0) then
                     pbTiz(n,k) = pbTiz(n,k)/vice(i,j,iblk)
                     piDin(n,k) = piDin(n,k)/vice(i,j,iblk) 
                  endif
               enddo                 ! k
               do k = 1, nblyr+2
                  pphin(n,k) = c0
                  pphin1(n,k) = c0
                  if (aicen(i,j,1,iblk) > c0) pphin1(n,k) = bphi(i,j,k,1,iblk)
                  do nn = 1,ncat
                     pphin(n,k) = pphin(n,k) + bphi(i,j,k,nn,iblk)*vicen(i,j,nn,iblk)
                  enddo
                  if (vice(i,j,iblk) > c0) then
                     pphin(n,k) = pphin(n,k)/vice(i,j,iblk)
                  endif
               enddo
               do k = 1,nblyr
                  pSin(n,k) = c0
                  pSin1(n,k) = c0                    
                  pSin(n,k)= trcr(i,j,nt_bgc_S+k-1,iblk)     
                  if (aicen(i,j,1,iblk) > c0) pSin1(n,k) = trcrn(i,j,nt_bgc_S+k-1,1,iblk)
               enddo 
               do k = 1,nilyr
                  pSice(n,k) = trcr(i,j,nt_sice+k-1,iblk)
               enddo
            endif                 ! my_task = pmloc

            call broadcast_scalar(phinS    (n),   pmloc(n)) 
            call broadcast_scalar(phinS1   (n),   pmloc(n)) 
            call broadcast_scalar(phbrn    (n),   pmloc(n)) 
            call broadcast_scalar(pdh_top1 (n),   pmloc(n)) 
            call broadcast_scalar(pdh_bot1 (n),   pmloc(n)) 
            call broadcast_scalar(psice_rho(n),   pmloc(n))  
            call broadcast_scalar(pfzsal_g (n),   pmloc(n))  
            call broadcast_scalar(pdarcy_V1(n),   pmloc(n)) 
            call broadcast_scalar(pfzsal   (n),   pmloc(n)) 
            call broadcast_array (pbTiz    (n,:), pmloc(n))
            call broadcast_array (piDin    (n,:), pmloc(n))
            call broadcast_array (pphin    (n,:), pmloc(n))
            call broadcast_array (pphin1   (n,:), pmloc(n))
            call broadcast_array (pSin     (n,:), pmloc(n))
            call broadcast_array (pSin1    (n,:), pmloc(n))
            call broadcast_array (pSice    (n,:), pmloc(n))
         enddo                  ! npnt
      endif                     ! print_points

      !-----------------------------------------------------------------
      ! start spewing
      !-----------------------------------------------------------------

      if (my_task == master_task) then

      call flush_fileunit(nu_diag)

      !-----------------------------------------------------------------
      ! diagnostics for Arctic and Antarctic points
      !-----------------------------------------------------------------

      if (print_points) then

        write(nu_diag,*) '                         '
        write(nu_diag,902) '      Brine height       '
        write(nu_diag,900) 'hbrin                   = ',phinS(1),phinS(2)
        write(nu_diag,900) 'hbrin cat 1             = ',phinS1(1),phinS1(2)
        write(nu_diag,900) 'Freeboard               = ',phbrn(1),phbrn(2)
        write(nu_diag,900) 'dhbrin cat 1 top        = ',pdh_top1(1),pdh_top1(2)
        write(nu_diag,900) 'dhbrin cat 1 bottom     = ',pdh_bot1(1),pdh_bot1(2)
        write(nu_diag,*) '                         '
        write(nu_diag,902) '     zSalinity         '
        write(nu_diag,900) 'Avg density (kg/m^3)   = ',psice_rho(1),psice_rho(2)
        write(nu_diag,900) 'Salt flux (kg/m^2/s)   = ',pfzsal(1),pfzsal(2)
        write(nu_diag,900) 'Grav. Drain. Salt flux = ',pfzsal_g(1),pfzsal_g(2)
        write(nu_diag,900) 'Darcy V cat 1 (m/s)    = ',pdarcy_V1(1),pdarcy_V1(2)
        write(nu_diag,*) '                         '
        write(nu_diag,*) ' Top down bgc Layer Model'
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'bTiz(1) ice temp',' bTiz(2) ice temp  '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((pbTiz(n,k),n = 1,2), k = 1,nblyr+1)
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'iDi(1) diffusivity  ','iDi(2) diffusivity  '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((piDin(n,k),n=1,2), k = 1,nblyr+1)
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'bphi(1) porosity   ','bphi(2) porosity  '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((pphin(n,k),n=1,2), k = 1,nblyr+1)
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'phi1(1) porosity   ','phi1(2) porosity  '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((pphin1(n,k),n=1,2), k = 1,nblyr+1)
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'zsal(1) cat 1 ','zsal(2) cat 1 '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((pSin1(n,k),n=1,2), k = 1,nblyr)                         
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'zsal(1) Avg S ','zsal(2) Avg S '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((pSin(n,k),n=1,2), k = 1,nblyr)            
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'Sice(1) Ice S ','Sice(2) Ice S '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((pSice(n,k),n=1,2), k = 1,nilyr)            
        write(nu_diag,*) '                         '

      endif                     ! print_points
      endif                     ! my_task = master_task

  802 format (f24.17,2x,f24.17)
  803 format (a25,2x,a25)
  900 format (a25,2x,f24.17,2x,f24.17)
  902 format (a25,10x,f6.1,1x,f6.1,9x,f6.1,1x,f6.1)
  903 format (a25,5x,i4,1x,i4,1x,i4,1x,i4,7x,i4,1x,i4,1x,i4,1x,i4)

      end subroutine zsal_diags

!=======================================================================

      end module ice_diagnostics_bgc

!=======================================================================
