! converts a restart file from CICE 4.1 to 5.0 format
! single processor (not parallelized)
! core restarts only (not including tracers except Tsfc)

! search for "modify as needed" to find potential configuration changes.
! NB: restarts will likely work only for exactly the same configuration 
!     (BL99 thermo in particular).

! intel compiler:
! ifort convert_restarts.f90 -o convert_restarts -i4 -r8 -convert big_endian -assume byterecl

      program convert_restarts

      implicit none

      integer, parameter :: char_len_long  = 256, &
                            log_kind  = kind(.true.), &
                            int_kind  = selected_int_kind(6), &
                            dbl_kind  = selected_real_kind(13)

!!!!!!! modify as needed (begin)
      ! these values are for the standard gx1 configuration
      integer (kind=int_kind), parameter :: &
                            ncat = 5, &       ! number of thickness cats
                            nilyr = 4, &      ! number of ice layers
                            nslyr = 1, &      ! number of snow layers
                            nx_block = 320, & ! global grid size
                            ny_block = 384

      ! these values are for the standard gx3 configuration
!      integer (kind=int_kind), parameter :: &
!                            ncat = 5, &       ! number of thickness cats
!                            nilyr = 4, &      ! number of ice layers
!                            nslyr = 1, &      ! number of snow layers
!                            nx_block = 100, & ! global grid size
!                            ny_block = 116

      ! flags
      logical (kind=log_kind), parameter :: &
         oceanmixed_ice = .true., & ! if true, read/write ocean mixed layer fields
         heat_capacity  = .true., & ! if true, ice has nonzero heat capacity
                                    ! if false, use zero-layer thermodynamics
         diag = .true.              ! write min/max diagnostics for fields

      ! file names
      character (len=char_len_long), parameter :: &
         iced_4_1 = '/scratch/eclare/tmp/restarts/iced_gx1_v4.0_kcatbound0', &  ! gx1
         iced_5_0 = '/scratch/eclare/tmp/restarts/iced_gx1_v4.0_kcatbound0_converted'
!         iced_4_1 = 'iced_gx3_v4.0_kcatbound0', &  ! gx3
!         iced_5_0 = 'iced_gx3_v4.0_kcatbound0_converted''

!!!!!!! modify as needed (end)

      ! array sizes
      integer (kind=int_kind), parameter :: &
         ntilyr    = ncat*nilyr, & ! number of ice layers in all categories
         ntslyr    = ncat*nslyr, & ! number of snow layers in all categories
         max_ntrcr =   1         & ! 1 = surface temperature              
                   + nilyr       & ! ice salinity
                   + nilyr       & ! ice enthalpy
                   + nslyr         ! snow enthalpy

      integer (kind=int_kind), dimension(ncat) :: &
         ilyr1          , & ! starting ice layer number for each category
         slyr1              ! starting snow layer number for each category

      integer (kind=int_kind) :: &
         ntrcr          , & ! number of required tracers
         i,j,n,k,m          ! indices

      ! tracer indices
      integer (kind=int_kind) :: &
         nt_Tsfc  , & ! ice/snow temperature
         nt_qice  , & ! volume-weighted ice enthalpy (in layers)
         nt_qsno  , & ! volume-weighted snow enthalpy (in layers)
         nt_sice      ! volume-weighted ice bulk salinity (CICE grid layers)

      ! time info
      integer (kind=int_kind) :: &
         istep1    ! counter, number of steps at current timestep

      real (kind=dbl_kind) :: &
         time           , & ! total elapsed time (s)
         time_forc          ! time of last forcing update (s)

      ! restart fields
      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ncat) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,max_ntrcr,ncat) :: &
         trcrn     ! tracers
                   ! 1: surface temperature of ice/snow (C)

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ntilyr) :: &
         eicen     ! energy of melting for each ice layer  (J/m^2)

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ntslyr) :: &
         esnon     ! energy of melting for each snow layer (J/m^2)

      real (kind=dbl_kind), dimension(nx_block,ny_block) :: &
         uvel    , & ! x-component of velocity (m/s)
         vvel    , & ! y-component of velocity (m/s)
         swvdr   , & ! sw down, visible, direct  (W/m^2)
         swvdf   , & ! sw down, visible, diffuse (W/m^2)
         swidr   , & ! sw down, near IR, direct  (W/m^2)
         swidf   , & ! sw down, near IR, diffuse (W/m^2)
         scale_factor, &! scaling factor for shortwave components
         strocnxT, & ! ice-ocean stress, x-direction
         strocnyT, & ! ice-ocean stress, y-direction
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4, & ! sigma12
         sst     , & ! sea surface temperature (C)
         frzmlt  , & ! freezing/melting potential (W/m^2)
         iceumask    ! ice extent mask (U-cell)

      ! flags
      logical (kind=log_kind) :: &
         l_brine     ! if true, treat brine pocket effects

      ! numbers
      real (kind=dbl_kind), parameter :: &
         spval_dbl = 1.0e30_dbl_kind, & ! special value (double precision)
         puny = 1.0e-11_dbl_kind, &
         pi   = 3.14159265358979323846_dbl_kind,&! pi
         c0   = 0.0_dbl_kind, &
         c1   = 1.0_dbl_kind, &
         c2   = 2.0_dbl_kind, &
         p5   = 0.5_dbl_kind

      ! physical parameters
      real (kind=dbl_kind), parameter :: &
         Lsub      = 2.835e6_dbl_kind ,&! latent heat, sublimation freshwater (J/kg)
         Lvap      = 2.501e6_dbl_kind ,&! latent heat, vaporization freshwater (J/kg)
         Lfresh    = Lsub-Lvap        ,&! latent heat of melting of fresh ice (J/kg)
         cp_ice    = 2106._dbl_kind   ,&! specific heat of fresh ice (J/kg/K)
         rhos      = 330.0_dbl_kind   ,&! density of snow (kg/m^3)
         hs_min  = 1.e-4_dbl_kind, &  ! min snow thickness for computing Tsno (m)
         nsal    = 0.407_dbl_kind, &
         msal    = 0.573_dbl_kind, &
         min_salin = 0.1_dbl_kind, &  ! threshold for brine pocket treatment 
         saltmax = 3.2_dbl_kind       ! max salinity at ice base (ppt)

      ! useful temps
      real (kind=dbl_kind), dimension(nilyr+1) :: &
         salin       ! initial salinity  profile (ppt)   

      real (kind=dbl_kind)    :: &
         Tmin, Tmax,    & ! min and max snow temperature
         zTsn,          & ! snow temperature
         zn,            & ! thickness
         rnslyr,        & ! real(nslyr)
         rnilyr           ! real(nilyr)

      ! count tracers
         nt_Tsfc = 1           ! index tracers, starting with Tsfc = 1
         ntrcr = 1             ! count tracers, starting with Tsfc = 1
         nt_qice = ntrcr + 1
         ntrcr = ntrcr + nilyr ! qice in nilyr layers
         nt_qsno = ntrcr + 1
         ntrcr = ntrcr + nslyr ! qsno in nslyr layers
         nt_sice = ntrcr + 1
         ntrcr = ntrcr + nilyr ! sice in nilyr layers
         if (ntrcr /= max_ntrcr) write (*,*) 'Tracer number mismatch'

      ! ice salinity profile
      if (saltmax > min_salin .and. heat_capacity) then
         l_brine = .true.
      else
         l_brine = .false.
      endif

      if (l_brine) then
         do k = 1, nilyr
            zn = (real(k,kind=dbl_kind)-p5) /  &
                  real(nilyr,kind=dbl_kind)
            salin(k)=(saltmax/c2)*(c1-cos(pi*zn**(nsal/(msal+zn))))
!            salin(k)=saltmax ! for isosaline ice ! modify as needed
         enddo
         salin(nilyr+1) = saltmax
      else
         do k = 1, nilyr+1
            salin(k) = c0
         enddo
      endif
      do k = 1, nilyr
         trcrn(:,:,nt_sice+k-1,:) = salin(k)
      enddo

      ! read cice 4.1 restart data
      call restartfile_4_1 (iced_4_1)

      ! convert eicen, esnon to qicen, qsnon
      ilyr1(1) = 1
      slyr1(1) = 1
      do n = 2, ncat
         ilyr1(n) = ilyr1(n-1) + nilyr
         slyr1(n) = slyr1(n-1) + nslyr
      enddo

      rnslyr = real(nslyr, kind=dbl_kind)
      rnilyr = real(nilyr, kind=dbl_kind)
      Tmin = -100.  ! minimum allowed snow temperature
      do j = 1, ny_block
      do i = 1, nx_block
      do n = 1, ncat
      if (aicen(i,j,n) > puny) then
      if (vsnon(i,j,n)/aicen(i,j,n) > hs_min) then
         do k = 1, nslyr
            ! qsn, esnon < 0              
            trcrn(i,j,nt_qsno+k-1,n) = &
               min(esnon(i,j,slyr1(n)+k-1)*rnslyr/vsnon(i,j,n),-rhos*Lfresh) 
            Tmax = -trcrn(i,j,nt_qsno+k-1,n)*puny*rnslyr/(rhos*cp_ice*vsnon(i,j,n))
            if (.not. heat_capacity) then
               trcrn(i,j,nt_qsno+k-1,n) = -rhos * Lfresh
               Tmax = puny
            endif
!!!!!!! modify as needed (begin)
! if your restarts do not work, try uncommenting this section
!            ! snow temperature
!            zTsn = (Lfresh + trcrn(i,j,nt_qsno+k-1,n)/rhos)/cp_ice
!            ! zap entire snow volume if temperature is out of bounds
!            if (zTsn < Tmin .or. zTsn > Tmax) then
!               print*, 'zapping snow volume ', i,j,n,vsnon(i,j,n)
!               vsnon(i,j,n) = c0
!               do m = 1, nslyr
!                  trcrn(i,j,nt_qsno+m-1,n) = c0
!               enddo
!            endif
!!!!!!! modify as needed (begin)
         enddo
      else
         vsnon(i,j,n) = c0
         do k = 1, nslyr
            trcrn(i,j,nt_qsno+k-1,n) = c0
         enddo
      endif
      endif
      do k = 1, nilyr
!         if (vicen(i,j,n) > puny) then
         if (aicen(i,j,n) > puny) then  ! matches v4.1
            trcrn(i,j,nt_qice+k-1,n) = eicen(i,j,ilyr1(n)+k-1)*rnilyr/vicen(i,j,n)
         else
            trcrn(i,j,nt_qice+k-1,n) = c0
         endif
      enddo
      enddo
      enddo
      enddo

      ! write cice 5.0 restart data
      call dumpfile_5_0 (iced_5_0)

!=======================================================================

      contains

!=======================================================================
!---! these subroutines write/read Fortran unformatted data files
!=======================================================================

      subroutine restartfile_4_1 (ice_ic)

      character (*) :: ice_ic

      integer (kind=int_kind) :: i, j, k, n

      character(len=char_len_long) :: &
         filename

      integer (kind=int_kind), parameter :: nu_restart = 20

         filename = ice_ic
         open(nu_restart,file=filename,form='unformatted')

         write(*,*) 'Using restart dump=', trim(filename)
         read (nu_restart) istep1,time,time_forc
         write(*,*) 'Restart read at istep=',istep1,time,time_forc

      !-----------------------------------------------------------------
      ! state variables
      ! Tsfc is the only tracer read in this file.  All other
      ! tracers are in their own dump/restart files.
      !-----------------------------------------------------------------
      do n=1,ncat
              write(*,*) 'cat ',n, &
                               ' min/max area, vol ice, vol snow, Tsfc'

         call ice_read(nu_restart,aicen(:,:,n),diag)
         call ice_read(nu_restart,vicen(:,:,n),diag)
         call ice_read(nu_restart,vsnon(:,:,n),diag)
         call ice_read(nu_restart,trcrn(:,:,nt_Tsfc,n),diag)
      enddo

           write(*,*) 'min/max eicen for each layer'
      do k=1,ntilyr
         call ice_read(nu_restart,eicen(:,:,k),diag)
      enddo

           write(*,*) 'min/max esnon for each layer'
      do k=1,ntslyr
         call ice_read(nu_restart,esnon(:,:,k),diag)
      enddo

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
           write(*,*) 'min/max velocity components'

      call ice_read(nu_restart,uvel,diag)
      call ice_read(nu_restart,vvel,diag)

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------
         write(*,*) 'radiation fields'

      call ice_read(nu_restart,scale_factor,diag)
      call ice_read(nu_restart,swvdr,diag)
      call ice_read(nu_restart,swvdf,diag)
      call ice_read(nu_restart,swidr,diag)
      call ice_read(nu_restart,swidf,diag)

      !-----------------------------------------------------------------
      ! ocean stress
      !-----------------------------------------------------------------
           write(*,*) 'min/max ocean stress components'

      call ice_read(nu_restart,strocnxT,diag)
      call ice_read(nu_restart,strocnyT,diag)

      !-----------------------------------------------------------------
      ! internal stress
      !-----------------------------------------------------------------
           write(*,*) 'internal stress components'
      
      call ice_read(nu_restart,stressp_1,diag) ! stressp_1
      call ice_read(nu_restart,stressp_3,diag) ! stressp_3

      call ice_read(nu_restart,stressp_2,diag) ! stressp_2
      call ice_read(nu_restart,stressp_4,diag) ! stressp_4

      call ice_read(nu_restart,stressm_1,diag) ! stressm_1
      call ice_read(nu_restart,stressm_3,diag) ! stressm_3

      call ice_read(nu_restart,stressm_2,diag) ! stressm_2
      call ice_read(nu_restart,stressm_4,diag) ! stressm_4

      call ice_read(nu_restart,stress12_1,diag) ! stress12_1
      call ice_read(nu_restart,stress12_3,diag) ! stress12_3

      call ice_read(nu_restart,stress12_2,diag) ! stress12_2
      call ice_read(nu_restart,stress12_4,diag) ! stress12_4

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
           write(*,*) 'ice mask for dynamics'

      call ice_read(nu_restart,iceumask,diag)

      ! for mixed layer model
      if (oceanmixed_ice) then
              write(*,*) 'min/max sst, frzmlt'

         call ice_read(nu_restart,sst,diag)
         call ice_read(nu_restart,frzmlt,diag)
      endif

      !-----------------------------------------------------------------
      ! mask out all-land blocks
      !-----------------------------------------------------------------
      where (aicen        > 0.5*spval_dbl) aicen        = c0
      where (vicen        > 0.5*spval_dbl) vicen        = c0
      where (vsnon        > 0.5*spval_dbl) vsnon        = c0
      where (trcrn        > 0.5*spval_dbl) trcrn        = c0
      where (eicen        > 0.5*spval_dbl) eicen        = c0
      where (esnon        > 0.5*spval_dbl) esnon        = c0
      where (uvel         > 0.5*spval_dbl) uvel         = c0
      where (vvel         > 0.5*spval_dbl) vvel         = c0
      where (scale_factor > 0.5*spval_dbl) scale_factor = c0
      where (swvdr        > 0.5*spval_dbl) swvdr        = c0
      where (swvdf        > 0.5*spval_dbl) swvdf        = c0
      where (swidr        > 0.5*spval_dbl) swidr        = c0
      where (swidf        > 0.5*spval_dbl) swidf        = c0
      where (strocnxT     > 0.5*spval_dbl) strocnxT     = c0
      where (strocnyT     > 0.5*spval_dbl) strocnyT     = c0
      where (stressp_1    > 0.5*spval_dbl) stressp_1    = c0
      where (stressp_2    > 0.5*spval_dbl) stressp_2    = c0
      where (stressp_3    > 0.5*spval_dbl) stressp_3    = c0
      where (stressp_4    > 0.5*spval_dbl) stressp_4    = c0
      where (stressm_1    > 0.5*spval_dbl) stressm_1    = c0
      where (stressm_2    > 0.5*spval_dbl) stressm_2    = c0
      where (stressm_3    > 0.5*spval_dbl) stressm_3    = c0
      where (stressm_4    > 0.5*spval_dbl) stressm_4    = c0
      where (stress12_1   > 0.5*spval_dbl) stress12_1   = c0
      where (stress12_2   > 0.5*spval_dbl) stress12_2   = c0
      where (stress12_3   > 0.5*spval_dbl) stress12_3   = c0
      where (stress12_4   > 0.5*spval_dbl) stress12_4   = c0
      where (iceumask     > 0.5*spval_dbl) iceumask     = c0
      if (oceanmixed_ice) then
         where (sst       > 0.5*spval_dbl) sst          = c0
         where (frzmlt    > 0.5*spval_dbl) frzmlt       = c0
      endif

      close(nu_restart)

      end subroutine restartfile_4_1

!=======================================================================

      subroutine dumpfile_5_0(filename_spec)

      character(len=char_len_long), intent(in) :: filename_spec

      integer (kind=int_kind) :: i, j, k, n

      character(len=char_len_long) :: filename

      integer (kind=int_kind), parameter :: nu_dump = 21

        filename = trim(filename_spec)
        open(nu_dump,file=filename,form='unformatted')

        write(nu_dump) istep1,time,time_forc
        write(*,*) 'Writing ',trim(filename)
        write(*,*) 'Restart written ',istep1,time,time_forc

      !-----------------------------------------------------------------
      ! state variables
      ! Tsfc is the only tracer written to this file.  All other
      ! tracers are written to their own dump/restart files.
      !-----------------------------------------------------------------

      do n=1,ncat
              write(*,*) 'cat ',n, &
                               ' min/max area, vol ice, vol snow, Tsfc'

         call ice_write(nu_dump,aicen(:,:,n),diag)
         call ice_write(nu_dump,vicen(:,:,n),diag)
         call ice_write(nu_dump,vsnon(:,:,n),diag)
         call ice_write(nu_dump,trcrn(:,:,nt_Tsfc,n),diag)

           write(*,*) 'min/max sicen for each layer'
         do k=1,nilyr
            call ice_write(nu_dump,trcrn(:,:,nt_sice+k-1,n),diag)
         enddo

           write(*,*) 'min/max qicen for each layer'
         do k=1,nilyr
            call ice_write(nu_dump,trcrn(:,:,nt_qice+k-1,n),diag)
         enddo

           write(*,*) 'min/max qsnon for each layer'
         do k=1,nslyr
            call ice_write(nu_dump,trcrn(:,:,nt_qsno+k-1,n),diag)
         enddo
      enddo ! ncat

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
           write(*,*) 'min/max velocity components'

      call ice_write(nu_dump,uvel,diag)
      call ice_write(nu_dump,vvel,diag)

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------
         write(*,*) 'radiation fields'

      call ice_write(nu_dump,scale_factor,diag)
      call ice_write(nu_dump,swvdr,diag)
      call ice_write(nu_dump,swvdf,diag)
      call ice_write(nu_dump,swidr,diag)
      call ice_write(nu_dump,swidf,diag)

      !-----------------------------------------------------------------
      ! ocean stress (for bottom heat flux in thermo)
      !-----------------------------------------------------------------
           write(*,*) 'min/max ocean stress components'

      call ice_write(nu_dump,strocnxT,diag)
      call ice_write(nu_dump,strocnyT,diag)

      !-----------------------------------------------------------------
      ! internal stress
      !-----------------------------------------------------------------
           write(*,*) 'internal stress components'

      call ice_write(nu_dump,stressp_1,diag)
      call ice_write(nu_dump,stressp_3,diag)
      call ice_write(nu_dump,stressp_2,diag)
      call ice_write(nu_dump,stressp_4,diag)

      call ice_write(nu_dump,stressm_1,diag)
      call ice_write(nu_dump,stressm_3,diag)
      call ice_write(nu_dump,stressm_2,diag)
      call ice_write(nu_dump,stressm_4,diag)

      call ice_write(nu_dump,stress12_1,diag)
      call ice_write(nu_dump,stress12_3,diag)
      call ice_write(nu_dump,stress12_2,diag)
      call ice_write(nu_dump,stress12_4,diag)

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
           write(*,*) 'ice mask for dynamics'
      
      call ice_write(nu_dump,iceumask,diag)

      ! for mixed layer model
      if (oceanmixed_ice) then
              write(*,*) 'min/max sst, frzmlt'

         call ice_write(nu_dump,sst,diag)
         call ice_write(nu_dump,frzmlt,diag)
      endif

      close(nu_dump)

      end subroutine dumpfile_5_0

!=======================================================================

      subroutine ice_read(nu, work_g1, diag)

      integer (kind=int_kind), intent(in) :: &
           nu                ! unit number

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
           intent(out) :: &
           work_g1           ! output array (real, 8-byte)

      logical (kind=log_kind), intent(in) :: &
           diag              ! if true, write diagnostic output

      integer (kind=int_kind) :: i, j

      real (kind=dbl_kind) :: &
         amin, amax         ! min and max values of input array

               read(nu) ((work_g1(i,j),i=1,nx_block),j=1,ny_block)

      if (diag) then
         amin = minval(work_g1)
         amax = maxval(work_g1, mask = work_g1 /= spval_dbl)
         write(*,*) ' read_global ',nu, amin, amax
      endif

      end subroutine ice_read

!=======================================================================

      subroutine ice_write(nu, work_g1, diag)

      integer (kind=int_kind), intent(in) :: &
           nu                ! unit number

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
           intent(in) :: &
           work_g1           ! input array (real, 8-byte)

      logical (kind=log_kind), intent(in) :: &
           diag              ! if true, write diagnostic output
!
      integer (kind=int_kind) :: i, j

      real (kind=dbl_kind) :: &
         amin, amax     ! min and max values of ouput array

            write(nu) ((work_g1(i,j),i=1,nx_block),j=1,ny_block)

         if (diag) then
            amin = minval(work_g1)
            amax = maxval(work_g1)
            write(*,*) ' write_global ', nu, amin, amax
         endif

      end subroutine ice_write

!=======================================================================

      end program convert_restarts

!=======================================================================
