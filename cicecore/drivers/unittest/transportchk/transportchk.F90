
      program transportchk

      ! This tests state_to_tracers and tracers_to_state in ice_transport_driver.
      ! CICE_InitMod (from the standalone model) is used to initialize the
      ! tracer configuration.  The routines are then called on one block of
      ! hardwired inputs and the results are checked against known values:
      !   state_to_tracers fills aim and trm, trm is zero where aicen <= puny
      !   round trip state -> tracers -> state recovers the state where aicen > puny
      !   round trip where 0 < aicen <= puny keeps aicen, zeroes volumes
      !   tracers_to_state where ice is removed (aim = 0) zeroes aicen, vicen, vsnon
      !     and leaves trcrn unchanged
      !   tracers_to_state where ice is added to an ice-free cell
      !   tracers_to_state where a cell is ice-free before and after
      !   round trip for a block with no ice

      use CICE_InitMod
      use ice_kinds_mod, only: int_kind, dbl_kind
      use ice_blocks, only: nx_block, ny_block
      use ice_constants, only: c0, c1, p5
      use ice_communicate, only: my_task, master_task, MPI_COMM_ICE
      use ice_domain_size, only: ncat, nslyr
      use ice_exit, only: abort_ice, end_run
      use ice_fileunits, only: nu_diag
      use ice_global_reductions, only: global_maxval
      use ice_transport_driver, only: state_to_tracers, tracers_to_state
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_query_tracer_sizes, icepack_query_tracer_indices

      implicit none

      integer(int_kind) :: i, j, n, it, ntrcr, ntrace, nt_qsno
      integer(int_kind) :: cnt(3)
      real(dbl_kind) :: puny, rhos, Lfresh, errmax

      real(dbl_kind), allocatable :: &
         aice0(:,:), aicen(:,:,:), vicen(:,:,:), vsnon(:,:,:), trcrn(:,:,:,:), &
         aice0r(:,:), aicenr(:,:,:), vicenr(:,:,:), vsnonr(:,:,:), trcrnr(:,:,:,:), &
         aim(:,:,:), trm(:,:,:,:)

      integer(int_kind), parameter :: ntests = 7
      integer(int_kind) :: errorflag(ntests), errorflag0, gflag
      character(len=64) :: stringflag(ntests)
      integer(int_kind), parameter ::  &
         passflag = 0, &
         failflag = 1

      real(dbl_kind), parameter :: fillval = -1.0e36_dbl_kind
      real(dbl_kind), parameter :: errtol = 1.0e-14_dbl_kind
      real(dbl_kind), parameter :: aimadd = 0.2_dbl_kind   ! area added to ice-free cells
      character(len=*), parameter :: subname='(transportchk)'

      !-----------------------------------------------------------------
      ! Initialize CICE
      !-----------------------------------------------------------------

      call CICE_Initialize

      call icepack_query_parameters(puny_out=puny, rhos_out=rhos, Lfresh_out=Lfresh)
      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_query_tracer_indices(nt_qsno_out=nt_qsno)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)
      ntrace = 2 + ntrcr

      errorflag = passflag
      stringflag(1) = 'state_to_tracers aim and trm'
      stringflag(2) = 'round trip where aicen > puny'
      stringflag(3) = 'round trip where 0 < aicen <= puny'
      stringflag(4) = 'tracers_to_state ice removed'
      stringflag(5) = 'tracers_to_state ice added to ice-free cell'
      stringflag(6) = 'tracers_to_state ice-free before and after'
      stringflag(7) = 'round trip with no ice in block'

      if (my_task == master_task) then
         write(6,*) ' '
         write(6,*) '=========================================================='
         write(6,*) ' '
         write(6,*) 'RunningUnitTest TRANSPORTCHK'
         write(6,*) ' '
         write(6,*) ' nx_block, ny_block = ', nx_block, ny_block
         write(6,*) ' ncat, ntrcr, nslyr = ', ncat, ntrcr, nslyr
         write(6,*) ' nt_qsno, puny      = ', nt_qsno, puny
         write(6,*) ' '
      endif

      allocate(aice0 (nx_block,ny_block), aicen (nx_block,ny_block,ncat), &
               vicen (nx_block,ny_block,ncat), vsnon (nx_block,ny_block,ncat), &
               trcrn (nx_block,ny_block,ntrcr,ncat), &
               aice0r(nx_block,ny_block), aicenr(nx_block,ny_block,ncat), &
               vicenr(nx_block,ny_block,ncat), vsnonr(nx_block,ny_block,ncat), &
               trcrnr(nx_block,ny_block,ntrcr,ncat), &
               aim(nx_block,ny_block,0:ncat), trm(nx_block,ny_block,ntrace,ncat))

      !-----------------------------------------------------------------
      ! Hardwired state: cycle through ice-free, aicen < puny, aicen = puny,
      ! and normal ice, with some snow-free cells
      !-----------------------------------------------------------------

      call init_state(aice0, aicen, vicen, vsnon, trcrn)

      !-----------------------------------------------------------------
      ! Test 1: state_to_tracers
      ! aim is the area, trm is zero where aicen <= puny, otherwise
      ! hice, hsno, snow enthalpy + rhos*Lfresh, other tracers
      !-----------------------------------------------------------------

      aim = fillval
      trm = fillval
      call state_to_tracers (nx_block, ny_block, ntrcr, ntrace, &
                             aice0, aicen, trcrn, vicen, vsnon, aim, trm)

      if (any(aim(:,:,0) /= aice0)) errorflag(1) = failflag
      if (any(aim(:,:,1:ncat) /= aicen)) errorflag(1) = failflag
      errmax = c0
      do n = 1, ncat
      do j = 1, ny_block
      do i = 1, nx_block
         if (aicen(i,j,n) > puny) then
            errmax = max(errmax, relerr(trm(i,j,1,n), vicen(i,j,n)/aicen(i,j,n)))
            errmax = max(errmax, relerr(trm(i,j,2,n), vsnon(i,j,n)/aicen(i,j,n)))
            do it = 1, ntrcr
               if (trm(i,j,2+it,n) /= trcrn(i,j,it,n) + snwoff_of(it)) errorflag(1) = failflag
            enddo
         else
            if (any(trm(i,j,:,n) /= c0)) errorflag(1) = failflag
         endif
      enddo
      enddo
      enddo
      if (errmax > errtol) errorflag(1) = failflag
      if (my_task == master_task) write(6,*) ' test 1 max relative error hice, hsno = ', errmax

      !-----------------------------------------------------------------
      ! Tests 2 and 3: round trip state -> tracers -> state
      ! The state arrays start from fill values, so every cell must be set
      ! (trcrn is not set where aim = 0, so it is not checked there).
      !-----------------------------------------------------------------

      aice0r = fillval
      aicenr = fillval
      vicenr = fillval
      vsnonr = fillval
      trcrnr = fillval
      call tracers_to_state (nx_block, ny_block, ntrcr, ntrace, &
                             aim, trm, aice0r, aicenr, trcrnr, vicenr, vsnonr)

      if (any(aice0r /= aice0)) errorflag(2) = failflag
      if (any(aicenr /= aicen)) errorflag(2) = failflag
      errmax = c0
      cnt = 0
      do n = 1, ncat
      do j = 1, ny_block
      do i = 1, nx_block
         if (aicen(i,j,n) > puny) then
            cnt(1) = cnt(1) + 1
            errmax = max(errmax, relerr(vicenr(i,j,n), vicen(i,j,n)))
            errmax = max(errmax, relerr(vsnonr(i,j,n), vsnon(i,j,n)))
            do it = 1, ntrcr
               errmax = max(errmax, relerr(trcrnr(i,j,it,n), trcrn(i,j,it,n)))
            enddo
         elseif (aicen(i,j,n) > c0) then
            ! tiny ice keeps its area, volumes and tracers are zeroed
            cnt(2) = cnt(2) + 1
            if (vicenr(i,j,n) /= c0 .or. vsnonr(i,j,n) /= c0) errorflag(3) = failflag
            do it = 1, ntrcr
               if (trcrnr(i,j,it,n) /= c0 - snwoff_of(it)) errorflag(3) = failflag
            enddo
         else
            cnt(3) = cnt(3) + 1
            if (vicenr(i,j,n) /= c0 .or. vsnonr(i,j,n) /= c0) errorflag(2) = failflag
         endif
      enddo
      enddo
      enddo
      if (errmax > errtol) errorflag(2) = failflag
      if (cnt(1) == 0) errorflag(2) = failflag
      if (cnt(2) == 0) errorflag(3) = failflag
      if (my_task == master_task) then
         write(6,*) ' test 2 cells aicen > puny, max relative error = ', cnt(1), errmax
         write(6,*) ' test 3 cells 0 < aicen <= puny = ', cnt(2)
         write(6,*) ' test 2 cells aicen = 0 = ', cnt(3)
      endif

      !-----------------------------------------------------------------
      ! Tests 4, 5 and 6: tracers_to_state after a change in aim similar
      ! to transport.  aim, trm are from state_to_tracers above, then
      !   ice removed (aim = 0) from some cells with aicen > puny
      !   ice added (aim = aimadd) to some cells with aicen = 0
      !   other cells with aicen = 0 stay ice-free
      ! The state arrays start from the original state.
      !-----------------------------------------------------------------

      do n = 1, ncat
      do j = 1, ny_block
      do i = 1, nx_block
         if (aicen(i,j,n) > puny .and. mod(i + 2*j + n, 5) == 0) then
            aim(i,j,n) = c0
         elseif (aicen(i,j,n) == c0 .and. mod(i+j,2) == 0) then
            aim(i,j,n) = aimadd
            trm(i,j,1,n) = 1.1_dbl_kind
            trm(i,j,2,n) = 0.05_dbl_kind
            do it = 1, ntrcr
               trm(i,j,2+it,n) = -0.5_dbl_kind*real(it,dbl_kind) + snwoff_of(it)
            enddo
         endif
      enddo
      enddo
      enddo

      aice0r = aice0
      aicenr = aicen
      vicenr = vicen
      vsnonr = vsnon
      trcrnr = trcrn
      call tracers_to_state (nx_block, ny_block, ntrcr, ntrace, &
                             aim, trm, aice0r, aicenr, trcrnr, vicenr, vsnonr)

      cnt = 0
      do n = 1, ncat
      do j = 1, ny_block
      do i = 1, nx_block
         if (aicen(i,j,n) > c0 .and. aim(i,j,n) == c0) then
            ! ice removed: area and volumes zeroed, tracers unchanged
            cnt(1) = cnt(1) + 1
            if (aicenr(i,j,n) /= c0 .or. vicenr(i,j,n) /= c0 .or. vsnonr(i,j,n) /= c0) &
               errorflag(4) = failflag
            if (any(trcrnr(i,j,:,n) /= trcrn(i,j,:,n))) errorflag(4) = failflag
         elseif (aicen(i,j,n) == c0 .and. aim(i,j,n) > c0) then
            ! ice added to an ice-free cell
            cnt(2) = cnt(2) + 1
            if (aicenr(i,j,n) /= aimadd) errorflag(5) = failflag
            if (vicenr(i,j,n) /= aimadd*1.1_dbl_kind)  errorflag(5) = failflag
            if (vsnonr(i,j,n) /= aimadd*0.05_dbl_kind) errorflag(5) = failflag
            do it = 1, ntrcr
               if (trcrnr(i,j,it,n) /= trm(i,j,2+it,n) - snwoff_of(it)) errorflag(5) = failflag
            enddo
         elseif (aicen(i,j,n) == c0 .and. aim(i,j,n) == c0) then
            ! ice-free before and after: area and volumes zero, tracers unchanged
            cnt(3) = cnt(3) + 1
            if (aicenr(i,j,n) /= c0 .or. vicenr(i,j,n) /= c0 .or. vsnonr(i,j,n) /= c0) &
               errorflag(6) = failflag
            if (any(trcrnr(i,j,:,n) /= trcrn(i,j,:,n))) errorflag(6) = failflag
         endif
      enddo
      enddo
      enddo
      if (any(aice0r /= aim(:,:,0))) errorflag(4) = failflag
      if (cnt(1) == 0) errorflag(4) = failflag
      if (cnt(2) == 0) errorflag(5) = failflag
      if (cnt(3) == 0) errorflag(6) = failflag
      if (my_task == master_task) then
         write(6,*) ' test 4 cells ice removed = ', cnt(1)
         write(6,*) ' test 5 cells ice added = ', cnt(2)
         write(6,*) ' test 6 cells ice-free before and after = ', cnt(3)
      endif

      !-----------------------------------------------------------------
      ! Test 7: round trip for a block with no ice, tracers keep their
      ! (arbitrary) values, everything else is zero
      !-----------------------------------------------------------------

      aicen = c0
      vicen = c0
      vsnon = c0
      aice0 = c1
      aim = fillval
      trm = fillval
      call state_to_tracers (nx_block, ny_block, ntrcr, ntrace, &
                             aice0, aicen, trcrn, vicen, vsnon, aim, trm)
      if (any(aim(:,:,0) /= c1) .or. any(aim(:,:,1:ncat) /= c0)) errorflag(7) = failflag
      if (any(trm /= c0)) errorflag(7) = failflag

      aice0r = fillval
      aicenr = fillval
      vicenr = fillval
      vsnonr = fillval
      trcrnr = trcrn
      call tracers_to_state (nx_block, ny_block, ntrcr, ntrace, &
                             aim, trm, aice0r, aicenr, trcrnr, vicenr, vsnonr)
      if (any(aice0r /= c1)) errorflag(7) = failflag
      if (any(aicenr /= c0) .or. any(vicenr /= c0) .or. any(vsnonr /= c0)) errorflag(7) = failflag
      if (any(trcrnr /= trcrn)) errorflag(7) = failflag

      !-----------------------------------------------------------------
      ! Report
      !-----------------------------------------------------------------

      errorflag0 = maxval(errorflag)
      gflag = global_maxval(errorflag0, MPI_COMM_ICE)
      errorflag0 = gflag
      do n = 1, ntests
         gflag = global_maxval(errorflag(n), MPI_COMM_ICE)
         errorflag(n) = gflag
      enddo

      if (my_task == master_task) then
         write(6,*) ' '
         write(6,*) 'TRANSPORTCHK COMPLETED SUCCESSFULLY'
         do n = 1, ntests
            if (errorflag(n) == passflag) then
               write(6,*) 'PASS ',trim(stringflag(n))
            else
               write(6,*) 'FAIL ',trim(stringflag(n))
            endif
         enddo
         if (errorflag0 == passflag) then
            write(6,*) 'TRANSPORTCHK TEST COMPLETED SUCCESSFULLY'
         else
            write(6,*) 'TRANSPORTCHK TEST FAILED'
         endif
         write(6,*) ' '
         write(6,*) '=========================================================='
         write(6,*) ' '
      endif

      !-----------------------------------------------------------------
      ! Gracefully end
      !-----------------------------------------------------------------

      call end_run()

      contains

!=======================================================================
! Hardwired state

      subroutine init_state(aice0, aicen, vicen, vsnon, trcrn)

      real(dbl_kind), intent(out) :: aice0(:,:), aicen(:,:,:), vicen(:,:,:), &
                                     vsnon(:,:,:), trcrn(:,:,:,:)

      integer(int_kind) :: i, j, n, it

      do n = 1, ncat
      do j = 1, ny_block
      do i = 1, nx_block
         select case (mod(i + 3*j + 5*n, 8))
         case (0)
            aicen(i,j,n) = c0
         case (1)
            aicen(i,j,n) = p5*puny
         case (2)
            aicen(i,j,n) = puny
         case (3)
            aicen(i,j,n) = c0
         case default
            aicen(i,j,n) = 0.01_dbl_kind + 0.1_dbl_kind*real(mod(i*j+n,9),dbl_kind)
         end select
         vicen(i,j,n) = aicen(i,j,n)*(0.3_dbl_kind + 0.7_dbl_kind*real(n,dbl_kind) &
                                    + 0.01_dbl_kind*real(i,dbl_kind))
         if (mod(i+j,3) == 0) then
            vsnon(i,j,n) = c0
         else
            vsnon(i,j,n) = aicen(i,j,n)*0.03_dbl_kind*real(mod(j,7)+1,dbl_kind)
         endif
         do it = 1, ntrcr
            if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
               trcrn(i,j,it,n) = -rhos*Lfresh - 1.0e5_dbl_kind*real(it+n,dbl_kind)
            else
               trcrn(i,j,it,n) = -1.7_dbl_kind + 0.37_dbl_kind*real(it,dbl_kind) &
                               + 0.013_dbl_kind*real(i-j,dbl_kind) + 0.11_dbl_kind*real(n,dbl_kind)
            endif
         enddo
      enddo
      enddo
      enddo
      aice0(:,:) = c1 - sum(aicen,dim=3)

      end subroutine init_state

!=======================================================================
! Offset added to snow enthalpy tracers by state_to_tracers

      real(dbl_kind) function snwoff_of(it)

      integer(int_kind), intent(in) :: it

      if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
         snwoff_of = rhos*Lfresh
      else
         snwoff_of = c0
      endif

      end function snwoff_of

!=======================================================================

      real(dbl_kind) function relerr(a, b)

      real(dbl_kind), intent(in) :: a, b

      if (b == c0) then
         relerr = abs(a)
      else
         relerr = abs(a-b)/abs(b)
      endif

      end function relerr

!=======================================================================

      end program transportchk

!=======================================================================
