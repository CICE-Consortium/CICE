!=======================================================================

! Calendar routines for managing time
!
! authors: Elizabeth C. Hunke, LANL
!          Tony Craig, NCAR
!          Craig MacLachlan, UK Met Office 
!
! 2006 ECH: Removed 'w' option for history; added 'h' and histfreq_n.
!           Converted to free form source (F90).
! 2010 CM : Fixed support for Gregorian calendar: subroutines
!           sec2time, time2sec and set_calendar added.

      module ice_calendar

      use ice_kinds_mod
      use ice_constants, only: c0, c1, c100, c30, c360, c365, c3600, &
          c4, c400
      use ice_domain_size, only: max_nstrm
      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters

      implicit none
      private

      public :: init_calendar, calendar, time2sec, sec2time, hc_jday

      integer (kind=int_kind), public :: &
         days_per_year        , & ! number of days in one year
         daymo(12)            , & ! number of days in each month
         daycal(13)               ! day number at end of month

      ! 360-day year data
      integer (kind=int_kind) :: &
         daymo360(12)         , & ! number of days in each month
         daycal360(13)            ! day number at end of month
      data daymo360 /   30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30/
      data daycal360/ 0,30, 60, 90,120,150,180,210,240,270,300,330,360/

      ! 365-day year data
      integer (kind=int_kind) :: &
         daymo365(12)         , & ! number of days in each month
         daycal365(13)            ! day number at end of month
      data daymo365 /   31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      data daycal365/ 0,31, 59, 90,120,151,181,212,243,273,304,334,365/

      ! 366-day year data (leap year)
      integer (kind=int_kind) :: &
         daymo366(12)         , & ! number of days in each month
         daycal366(13)            ! day number at end of month
      data daymo366 /   31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      data daycal366/ 0,31, 60, 91,121,152,182,213,244,274,305,335,366/

      real (kind=dbl_kind), parameter :: &
	days_per_4c = 146097.0_dbl_kind, &
	days_per_c  = 36524.0_dbl_kind,  &
	days_per_4y = 1461.0_dbl_kind,   &
	days_per_y  = 365.0_dbl_kind

      integer (kind=int_kind), public :: &
         istep    , & ! local step counter for time loop
         istep0   , & ! counter, number of steps taken in previous run
         istep1   , & ! counter, number of steps at current timestep
         mday     , & ! day of the month
         hour     , & ! hour of the day
         month    , & ! month number, 1 to 12
         monthp   , & ! last month
         year_init, & ! initial year
         nyr      , & ! year number
         idate    , & ! date (yyyymmdd)
         idate0   , & ! initial date (yyyymmdd)
         sec      , & ! elapsed seconds into date
         npt      , & ! total number of time steps (dt)
         ndtd     , & ! number of dynamics subcycles: dt_dyn=dt/ndtd
         stop_now     , & ! if 1, end program execution
         write_restart, & ! if 1, write restart now
         diagfreq     , & ! diagnostic output frequency (10 = once per 10 dt)
         dumpfreq_n   , & ! restart output frequency (10 = once per 10 d,m,y)
         nstreams     , & ! number of history output streams
         histfreq_n(max_nstrm) ! history output frequency 

      real (kind=dbl_kind), public :: &
         dt             , & ! thermodynamics timestep (s)
         dt_dyn         , & ! dynamics/transport/ridging timestep (s)
         time           , & ! total elapsed time (s)
         time_forc      , & ! time of last forcing update (s)
         yday           , & ! day of the year
         tday           , & ! absolute day number
         dayyr          , & ! number of days per year
         nextsw_cday    , & ! julian day of next shortwave calculation
         basis_seconds      ! Seconds since calendar zero

      logical (kind=log_kind), public :: &
         new_year       , & ! new year = .true.
         new_month      , & ! new month = .true.
         new_day        , & ! new day = .true.
         new_hour       , & ! new hour = .true.
         use_leap_years , & ! use leap year functionality if true
         write_ic       , & ! write initial condition now
         dump_last      , & ! write restart file on last time step
         force_restart_now, & ! force a restart now
         write_history(max_nstrm) ! write history now

      character (len=1), public :: &
         histfreq(max_nstrm), & ! history output frequency, 'y','m','d','h','1'
         dumpfreq               ! restart frequency, 'y','m','d'

      character (len=char_len), public :: &
         calendar_type       ! differentiates Gregorian from other calendars
                             ! default = ' '

!=======================================================================

      contains

!=======================================================================

! Initialize calendar variables
!
! authors: Elizabeth C. Hunke, LANL
!          Tony Craig, NCAR
!          Craig MacLachlan, UK Met Office

      subroutine init_calendar

      real    (kind=dbl_kind) :: secday           ! seconds per day

      character(len=*),parameter :: subname='(init_calendar)'

      call icepack_query_parameters(secday_out=secday)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      istep = 0         ! local timestep number
      time=istep0*dt    ! s
      yday=c0           ! absolute day number
      mday=0            ! day of the month
      month=0           ! month
      nyr=0             ! year
      idate=00000101    ! date
      sec=0             ! seconds into date
      istep1 = istep0   ! number of steps at current timestep
                        ! real (dumped) or imagined (use to set calendar)
      stop_now = 0      ! end program execution if stop_now=1
      dt_dyn = dt/real(ndtd,kind=dbl_kind) ! dynamics et al timestep
      force_restart_now = .false.

      ! Check that the number of days per year is set correctly when using
      ! leap years. If not, set days_per_year correctly and warn the user.
      if (use_leap_years .and. days_per_year /= 365) then
         days_per_year = 365
         write(nu_diag,*) 'Warning: days_per_year has been set to 365', &
              ' because use_leap_years = .true.'
      end if

#ifdef CESMCOUPLED
      ! calendar_type set by coupling
#else
      calendar_type = ' '
      if (use_leap_years .and. days_per_year == 365) calendar_type = 'Gregorian'
#endif
      
      dayyr = real(days_per_year, kind=dbl_kind)
      if (days_per_year == 360) then
         daymo  = daymo360
         daycal = daycal360
      elseif (days_per_year == 365) then
         daymo  = daymo365
         daycal = daycal365
      else 
         call abort_ice(subname//'ERROR: days_per_year must be 360 or 365')
      endif

      ! Get the time in seconds from calendar zero to start of initial year
      call time2sec(year_init,1,1,basis_seconds)

      ! determine initial date (assumes namelist year_init, istep0 unchanged)     
      sec = mod(time,secday)            ! elapsed seconds into date at
                                        ! end of dt
      tday = (time-sec)/secday + c1     ! absolute day number

      ! Convert the current timestep into a calendar date
      call sec2time(nyr,month,mday,basis_seconds+sec)

      yday = mday + daycal(month)  ! day of the year
      nyr = nyr - year_init + 1    ! year number

      idate0 = (nyr+year_init-1)*10000 + month*100 + mday ! date (yyyymmdd) 

      end subroutine init_calendar

!=======================================================================

! Determine the date at the end of the time step
!
! authors: Elizabeth C. Hunke, LANL
!          Tony Craig, NCAR
!          Craig MacLachlan, UK Met Office

      subroutine calendar(ttime)

      use ice_communicate, only: my_task, master_task

      real (kind=dbl_kind), intent(in) :: &
         ttime                          ! time variable

      ! local variables

      integer (kind=int_kind) :: &
         ns                         , & ! loop index
         nyrp,mdayp,hourp           , & ! previous year, day, hour
         elapsed_days               , & ! since beginning this run
         elapsed_months             , & ! since beginning this run
         elapsed_hours              , & ! since beginning this run
         month0
      real    (kind=dbl_kind) :: secday ! seconds per day
      character(len=*),parameter :: subname='(calendar)'

      call icepack_query_parameters(secday_out=secday)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      nyrp=nyr
      monthp=month
      mdayp=mday
      hourp=hour
      new_year=.false.
      new_month=.false.
      new_day=.false.
      new_hour=.false.
      write_history(:)=.false.
      write_restart=0

      sec = mod(ttime,secday)           ! elapsed seconds into date at
                                        ! end of dt
      tday = (ttime-sec)/secday + c1    ! absolute day number

      ! Deterime the current date from the timestep
      call sec2time(nyr,month,mday,basis_seconds+ttime)

      yday = mday + daycal(month)   ! day of the year
      nyr = nyr - year_init + 1     ! year number
      
      hour = int((ttime)/c3600) + c1 ! hour

      month0 = int((idate0 - int(idate0 / 10000) * 10000) / 100)

      elapsed_months = (nyr - 1)*12 + (month - month0)
      elapsed_days = int((istep * dt) / secday)
      elapsed_hours = int(ttime/3600)

      idate = (nyr+year_init-1)*10000 + month*100 + mday ! date (yyyymmdd) 

#ifndef CESMCOUPLED
      if (istep >= npt+1)  stop_now = 1
      if (istep == npt .and. dump_last) write_restart = 1 ! last timestep
#endif
      if (nyr   /= nyrp)   new_year = .true.
      if (month /= monthp) new_month = .true.
      if (mday  /= mdayp)  new_day = .true.
      if (hour  /= hourp)  new_hour = .true.


      do ns = 1, nstreams
         if (histfreq(ns)=='1' .and. histfreq_n(ns)/=0) then
             if (mod(istep1, histfreq_n(ns))==0) &
                write_history(ns)=.true.
         endif
      enddo

      if (istep > 1) then

        do ns = 1, nstreams

           select case (histfreq(ns))
           case ("y", "Y")
             if (new_year  .and. histfreq_n(ns)/=0) then
                if (mod(nyr, histfreq_n(ns))==0) &
                   write_history(ns) = .true.
             endif
           case ("m", "M")
             if (new_month .and. histfreq_n(ns)/=0) then
                if (mod(elapsed_months,histfreq_n(ns))==0) &
                   write_history(ns) = .true.
             endif
           case ("d", "D")
             if (new_day  .and. histfreq_n(ns)/=0) then
                if (mod(elapsed_days,histfreq_n(ns))==0) &
                   write_history(ns) = .true.
             endif
           case ("h", "H")
             if (new_hour  .and. histfreq_n(ns)/=0) then
                if (mod(elapsed_hours,histfreq_n(ns))==0) &
                   write_history(ns) = .true.
             endif
           end select

        enddo ! nstreams

        select case (dumpfreq)
        case ("y", "Y")
          if (new_year  .and. mod(nyr, dumpfreq_n)==0) &
                write_restart = 1
        case ("m", "M")
          if (new_month .and. mod(elapsed_months,dumpfreq_n)==0) &
                write_restart = 1
        case ("d", "D")
          if (new_day   .and. mod(elapsed_days, dumpfreq_n)==0) &
                write_restart = 1
        case ("h", "H")
          if (new_hour  .and. mod(elapsed_hours, dumpfreq_n)==0) &
                write_restart = 1
        end select

        if (force_restart_now) write_restart = 1
      
      endif !  istep > 1

      if (my_task == master_task .and. mod(istep,diagfreq) == 0 &
                                 .and. stop_now /= 1) then
        write(nu_diag,*) ' '
        write(nu_diag,'(a7,i10,4x,a6,i10,4x,a4,i10)') &
             'istep1:', istep1, 'idate:', idate, 'sec:', sec
      endif

      end subroutine calendar

!=======================================================================

! Convert the date to seconds since calendar zero.
!  ** This is based on the UM routine TIME2SEC **
!
! authors: Craig MacLachlan, UK Met Office

      subroutine time2sec(year,month,day,tsec)

      integer (kind=int_kind), intent(in) :: year  ! year
      integer (kind=int_kind), intent(in) :: month ! month
      integer (kind=int_kind), intent(in) :: day   ! year
      real (kind=dbl_kind),   intent(out) :: tsec  ! seconds since calendar zero

      ! local variables

      real    (kind=dbl_kind) :: days_since_calz   ! days since calendar zero
      real    (kind=dbl_kind) :: secday            ! seconds per day
      integer (kind=int_kind) :: years_since_calz  ! days since calendar zero
      character(len=*),parameter :: subname='(time2sec)'

      call icepack_query_parameters(secday_out=secday)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (dayyr == 360) then
         days_since_calz = c360*year + c30*(month-1) + day - c1
         tsec = secday * days_since_calz

      else
         
         if (use_leap_years) then

            call set_calendar(year)

            ! Add on the days from this year
            days_since_calz = day + daycal(month) - c1

            ! Subtract a year because we only want to count whole years
            years_since_calz = year - 1
            
            ! Add days from preceeding years
            days_since_calz  = days_since_calz &
                             + int(years_since_calz/c400)*days_per_4c
            years_since_calz = years_since_calz &
                             - int(years_since_calz/c400)*400

            days_since_calz  = days_since_calz &
                             + int(years_since_calz/c100)*days_per_c
            years_since_calz = years_since_calz &
                             - int(years_since_calz/c100)*100

            days_since_calz  = days_since_calz &
                             + int(years_since_calz/c4)*days_per_4y
            years_since_calz = years_since_calz &
                             - int(years_since_calz/c4)*4

            days_since_calz  = days_since_calz &
                             + years_since_calz*days_per_y

            tsec = secday * days_since_calz
            
         else ! Using fixed 365-day calendar
            
            days_since_calz = c365*year + daycal365(month) + day - c1
            tsec = secday * days_since_calz

         end if

      end if

      end subroutine time2sec

!=======================================================================

! Convert the time in seconds since calendar zero to a date.
!
! authors: Craig MacLachlan, UK Met Office

      subroutine sec2time(year,month,day,tsec)

      integer (kind=int_kind), intent(out) :: year     ! year
      integer (kind=int_kind), intent(out) :: month    ! month
      integer (kind=int_kind), intent(out) :: day      ! year
      real (kind=dbl_kind),    intent(in)  :: tsec     ! seconds since calendar zero      

      ! local variables

      real    (kind=dbl_kind) :: days_since_calz  ! days since calendar zero
      real    (kind=dbl_kind) :: secday           ! seconds per day
      integer (kind=int_kind) :: k                ! counter
      character(len=*),parameter :: subname='(sec2time)'

      call icepack_query_parameters(secday_out=secday)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      days_since_calz = int(tsec/secday)

      if (dayyr == 360) then

         year = int(days_since_calz/c360)
         month = mod(int(days_since_calz/c30),12) + 1
         day = mod(int(days_since_calz),30) + 1

      else

         if (use_leap_years) then

            year = int(days_since_calz/days_per_4c)*400
            days_since_calz = days_since_calz &
                            - int(days_since_calz/days_per_4c)*days_per_4c
            
            if (days_since_calz == 4*days_per_c) then
               year = year + 400
               days_since_calz = days_per_y + 1
            else
               year = year + int(days_since_calz/days_per_c)*100
               days_since_calz = days_since_calz &
                               - int(days_since_calz/days_per_c)*days_per_c
               
               year = year + int(days_since_calz/days_per_4y)*4
               days_since_calz = days_since_calz &
                               - int(days_since_calz/days_per_4y)*days_per_4y
               
               if (days_since_calz == 4*days_per_y) then
                  year = year + 4
                  days_since_calz = days_per_y + 1
               else
                  year = year + int(days_since_calz/days_per_y) + 1
                  days_since_calz = days_since_calz &
                                  - int(days_since_calz/days_per_y)*days_per_y + c1
               endif
            endif

            ! Ensure the calendar variables are correct for this year.
            call set_calendar(year)

            ! Calculate the month
            month = 1
            do k = 1, 12
               if (days_since_calz > daycal(k)) month = k
            enddo

            ! Calculate the day of the month
            day = days_since_calz - daycal(month)

         else ! Using fixed 365-day calendar
            
            year = int(days_since_calz/c365)
            days_since_calz = days_since_calz - year*365 + 1
      
            ! Calculate the month
            month = 1
            do k = 1, 12
               if (days_since_calz > daycal365(k)) month = k
            enddo

            ! Calculate the day of the month
            day = days_since_calz - daycal365(month)

         end if

      end if

      end subroutine sec2time

!=======================================================================

! Set the "days per month", "days per year", etc variables for the 
! current year.
!
! authors: Craig MacLachlan, UK Met Office

      subroutine set_calendar(year)

      integer (kind=int_kind), intent(in) :: year   ! current year

      ! Internal variable
      logical (kind=log_kind) :: isleap   ! Leap year logical
      character(len=*),parameter :: subname='(set_calendar)'

      isleap = .false. ! not a leap year
      if (mod(year,  4) == 0) isleap = .true.
      if (mod(year,100) == 0) isleap = .false.
      if (mod(year,400) == 0) isleap = .true.
      
      ! Ensure the calendar is set correctly
      if (isleap) then
         daycal = daycal366
         daymo = daymo366
         dayyr=real(daycal(13), kind=dbl_kind)
         days_per_year=int(dayyr)
      else
         daycal = daycal365
         daymo = daymo365
         dayyr=real(daycal(13), kind=dbl_kind)
         days_per_year=int(dayyr)
      endif

    end subroutine set_calendar

!=======================================================================

      real(kind=dbl_kind) function hc_jday(iyear,imm,idd,ihour)
!--------------------------------------------------------------------
! converts "calendar" date to HYCOM julian day:
!   1) year,month,day,hour  (4 arguments)
!   2) year,doy,hour        (3 arguments)
!
! HYCOM model day is calendar days since 31/12/1900
!--------------------------------------------------------------------
        real(kind=dbl_kind)     :: dtime
        integer(kind=int_kind)  :: iyear,iyr,imm,idd,idoy,ihr
        integer(kind=int_kind), optional :: ihour

        if (present(ihour)) then
          !-----------------
          ! yyyy mm dd HH
          !-----------------
          iyr=iyear-1901
          if (mod(iyr,4)==3) then
            dtime = floor(365.25_dbl_kind*iyr)*c1 + daycal366(imm)*c1 + idd*c1 + ihour/24._dbl_kind
          else
            dtime = floor(365.25_dbl_kind*iyr)*c1 + daycal365(imm)*c1 + idd*c1 + ihour/24._dbl_kind
          endif

        else
          !-----------------
          ! yyyy DOY HH
          !-----------------
          ihr   = idd   ! redefine input
          idoy  = imm   ! redefine input
          iyr   = iyear - 1901
          dtime = floor(365.25_dbl_kind*iyr)*c1 + idoy*c1 + ihr/24._dbl_kind

        endif

        hc_jday=dtime

        return
      end function hc_jday

!=======================================================================

      end module ice_calendar

!=======================================================================
