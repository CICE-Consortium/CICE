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
! 2020 TC : Significant refactor to move away from time as prognostic

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

      ! PUBLIC

      public :: init_calendar, advance_time, calendar, hc_jday
      public :: calendar_compute_elapsed_days, calendar_compute_days_between
      public :: calendar_set_date_from_timesecs
!      public :: date2rdate, rdate2date

      character(len=*), public, parameter :: &
         ice_calendar_gregorian = 'Gregorian', &  ! calendar name, actually proleptic gregorian here
         ice_calendar_noleap    = 'NO_LEAP', &    ! 365 day per year calendar
         ice_calendar_360day    = '360day'        ! 360 day calendar with 30 days per month

      integer (kind=int_kind), public, parameter :: &
         months_per_year = 12, &     ! months per year
         hours_per_day   = 24        ! hours per day

      integer (kind=int_kind), public :: &
         days_per_year         , & ! number of days in one year
         daymo(months_per_year), & ! number of days in each month
         daycal(months_per_year+1) ! accumulated days in year to end of prior month

      integer (kind=int_kind), public :: &
         ! step counters
         istep    , & ! local step counter for current run in time loop
         istep0   , & ! counter, number of steps at start of run
         istep1   , & ! counter, number of steps at current timestep
         ! basic time variables
         nyr      , & ! year number
         month    , & ! month number, 1 to months_per_year
         mday     , & ! day of the month
         sec      , & ! elapsed seconds into date
         ! initial time
         year_init, & ! initial year
         month_init,& ! initial month
         day_init, & ! initial day of month
         sec_init , & ! initial seconds
         ! other stuff
         idate    , & ! date (yyyymmdd)
         idate0   , & ! initial date (yyyymmdd), associated with year_init, month_init, day_init
         dayyr    , & ! number of days in the current year
         npt      , & ! total number of time steps (dt)
         ndtd     , & ! number of dynamics subcycles: dt_dyn=dt/ndtd
         stop_now     , & ! if 1, end program execution
         write_restart, & ! if 1, write restart now
         diagfreq     , & ! diagnostic output frequency (10 = once per 10 dt)
         dumpfreq_n   , & ! restart output frequency (10 = once per 10 d,m,y)
         nstreams     , & ! number of history output streams
         histfreq_n(max_nstrm) ! history output frequency

      logical (kind=log_kind), public :: &
         new_year       , & ! new year = .true.
         new_month      , & ! new month = .true.
         new_day        , & ! new day = .true.
         new_hour           ! new hour = .true.

      real (kind=dbl_kind), public :: &
         dt             , & ! thermodynamics timestep (s)
         dt_dyn         , & ! dynamics/transport/ridging timestep (s)
         timesecs       , & ! total elapsed time (s)
         rdate_forc     , & ! time of last forcing update (s)
         yday           , & ! day of the year
!tcx1         tday           , & ! absolute day number
!tcx1         dayyr          , & ! number of days per year
         nextsw_cday        ! julian day of next shortwave calculation
!tcx1         basis_seconds      ! Seconds since calendar zero

      logical (kind=log_kind), public :: &
         use_leap_years , & ! use leap year functionality if true
         write_ic       , & ! write initial condition now
         dump_last      , & ! write restart file on last time step
         force_restart_now, & ! force a restart now
         write_history(max_nstrm) ! write history now

      character (len=1), public :: &
         npt_unit,            & ! run length unit, 'y', 'm', 'd', 'h', 's', '1'
         histfreq(max_nstrm), & ! history output frequency, 'y','m','d','h','1'
         dumpfreq               ! restart frequency, 'y','m','d'

      character (len=char_len), public :: &
         calendar_type       ! differentiates Gregorian from other calendars
                             ! default = ' '

      ! PRIVATE

      integer (kind=int_kind) :: &
         hour         ! hour of the day

      ! 360-day year data
      integer (kind=int_kind) :: &
         daymo360(months_per_year)    ! number of days in each month
      data daymo360 /   30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30/

      ! 365-day year data
      integer (kind=int_kind) :: &
         daymo365(months_per_year)    ! number of days in each month
      data daymo365 /   31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

      ! 366-day year data (leap year)
      integer (kind=int_kind) :: &
         daymo366(months_per_year)    ! number of days in each month
      data daymo366 /   31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/


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
      nyr=year_init     ! year
      month=month_init  ! month
      mday=day_init     ! day of the month
      sec=sec_init      ! seconds into date
      istep1 = istep0   ! number of steps at current timestep
                        ! real (dumped) or imagined (use to set calendar)
      idate0 = (nyr)*10000 + month*100 + mday ! date (yyyymmdd) 
      stop_now = 0      ! end program execution if stop_now=1
      dt_dyn = dt/real(ndtd,kind=dbl_kind) ! dynamics et al timestep
      force_restart_now = .false.

#ifdef CESMCOUPLED
      ! calendar_type set by coupling
#else
      calendar_type = ''
      if (use_leap_years) then
         if (days_per_year == 365) then
            calendar_type = trim(ice_calendar_gregorian)
         else
            call abort_ice(subname//'ERROR: use_leap_years is true, must set days_per_year to 365')
         endif
      else
         if (days_per_year == 365) then
            calendar_type = trim(ice_calendar_noleap)
         elseif (days_per_year == 360) then
            calendar_type = trim(ice_calendar_360day)
         else
            call abort_ice(subname//'ERROR: days_per_year only 365 or 360 supported')
         endif
      endif
#endif

      call set_calendar(nyr)
      call calendar()

      end subroutine init_calendar

!=======================================================================

! Determine the date at the end of the time step
!
! authors: Elizabeth C. Hunke, LANL
!          Tony Craig, NCAR
!          Craig MacLachlan, UK Met Office

      subroutine advance_time()

      use ice_communicate, only: my_task, master_task

      ! local variables

      integer(kind=int_kind) :: &
         idt       ! integer dt
      character(len=*),parameter :: subname='(advance_time)'

      istep = istep + 1
      istep1 = istep1 + 1
      idt = nint(dt)
      sec = sec + idt

      call calendar()

      end subroutine advance_time

!=======================================================================

! Determine the date at the end of the time step
!
! authors: Elizabeth C. Hunke, LANL
!          Tony Craig, NCAR
!          Craig MacLachlan, UK Met Office

      subroutine calendar()

      use ice_communicate, only: my_task, master_task

!      real (kind=dbl_kind), intent(in), optional :: &
!         ttime                          ! time variable

      ! local variables

      integer (kind=int_kind) :: &
         ns                         , & ! loop index
         nyrp,monthp,mdayp,hourp    , & ! previous year, day, hour
         elapsed_days               , & ! since beginning this run
         elapsed_months             , & ! since beginning this run
         elapsed_hours                  ! since beginning this run
      real    (kind=dbl_kind) :: secday ! seconds per day
      integer (kind=int_kind) :: isecday  ! seconds per day
      character(len=*),parameter :: subname='(calendar)'

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

      !--- recompute yr/mon/day/sec

      call icepack_query_parameters(secday_out=secday)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      isecday = nint(secday)
      if (sec >= isecday) then
         mday = mday + int(sec/isecday)
         sec = mod(sec,isecday)
      endif
      do while (mday > daymo(month))
         mday = mday - daymo(month)
         month = month + 1
         do while (month > months_per_year)
            month = month - months_per_year
            nyr = nyr + 1
            call set_calendar(nyr)
         enddo
      enddo

      idate = (nyr)*10000 + month*100 + mday ! date (yyyymmdd) 
      yday = daycal(month) + mday            ! day of the year
      dayyr = daycal(months_per_year+1)
      elapsed_months = (nyr - year_init)*months_per_year + month - month_init
      elapsed_days = calendar_compute_days_between(year_init,month_init,day_init,nyr,month,mday)
      elapsed_hours = elapsed_days * hours_per_day
      timesecs = real(elapsed_days,kind=dbl_kind)*secday + &
                 real(sec,kind=dbl_kind)

      !--- compute other stuff

!tcx      hour = int((sec)/c3600) + c1 ! hour

#ifndef CESMCOUPLED
      if (istep >= npt+1)  stop_now = 1
      if (istep == npt .and. dump_last) write_restart = 1 ! last timestep
#endif
      if (nyr   /= nyrp)   new_year = .true.
      if (month /= monthp) new_month = .true.
      if (mday  /= mdayp)  new_day = .true.
      if (hour  /= hourp)  new_hour = .true.

      ! History writing flags

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
         case ("1")
            if (histfreq_n(ns)/=0) then
               if (mod(istep1, histfreq_n(ns))==0) &
                  write_history(ns)=.true.
            endif
         end select

      enddo

      ! Restart writing flag

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
      case ("1")
         if (mod(istep1, dumpfreq_n)==0) &
            write_restart = 1
      end select

      if (force_restart_now) write_restart = 1
      
      if (my_task == master_task .and. mod(istep1,diagfreq) == 0 &
                                 .and. stop_now /= 1) then
        write(nu_diag,*) ' '
        write(nu_diag,'(a7,i10,4x,a6,i10,4x,a4,i10)') &
             'istep1:', istep1, 'idate:', idate, 'sec:', sec
      endif

      end subroutine calendar

!=======================================================================
#if (1 == 0)

! Determine the date at the end of the time step
!
! authors: Elizabeth C. Hunke, LANL
!          Tony Craig, NCAR
!          Craig MacLachlan, UK Met Office

      subroutine calendar_old(ttime)

      use ice_communicate, only: my_task, master_task

      real (kind=dbl_kind), intent(in) :: &
         ttime                          ! time variable

      ! local variables

      integer (kind=int_kind) :: &
         ns                         , & ! loop index
         nyrp,mdayp,hourp           , & ! previous year, month, day, hour
         elapsed_days               , & ! since beginning this run
         elapsed_months             , & ! since beginning this run
         elapsed_hours              , & ! since beginning this run
         month0
      real    (kind=dbl_kind) :: secday ! seconds per day
      character(len=*),parameter :: subname='(calendar_old)'

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
!tcx1      tday = (ttime-sec)/secday + c1    ! absolute day number

      ! Deterime the current date from the timestep
      call sec2time(nyr,month,mday,basis_seconds+ttime)

      yday = mday + daycal(month)   ! day of the year
      nyr = nyr - year_init + 1     ! year number
      
      hour = int((ttime)/c3600) + c1 ! hour

      month0 = int((idate0 - int(idate0 / 10000) * 10000) / 100)

      elapsed_months = (nyr - 1)*months_per_year + (month - month0)
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

      if (dumpfreq == '1') then
         if (mod(istep1, dumpfreq_n)==0) &
            write_restart = 1
      endif

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

      end subroutine calendar_old
#endif
!=======================================================================
! Convert the date to a real value time

      subroutine date2rdate(year,month,day,sec,time)

      integer (kind=int_kind), intent(in) :: year  ! year
      integer (kind=int_kind), intent(in) :: month ! month
      integer (kind=int_kind), intent(in) :: day   ! day
      integer (kind=int_kind), intent(in) :: sec   ! seconds
      real (kind=dbl_kind),   intent(out) :: time  ! ymds expressed in real

      ! local variables

      character(len=*),parameter :: subname='(date2rdate)'

      ! Needs to be consistent with rdate2date

      time = real((year*10000 + month*100 + day),kind=dbl_kind) + &
             real(sec,kind=dbl_kind)/100000._dbl_kind

      end subroutine date2rdate

!=======================================================================
! Convert the real time value to a date

      subroutine rdate2date(time,year,month,day,sec)

      real (kind=dbl_kind),    intent(in)  :: time     ! ymds expressed in real
      integer (kind=int_kind), intent(out) :: year     ! year
      integer (kind=int_kind), intent(out) :: month    ! month
      integer (kind=int_kind), intent(out) :: day      ! day
      integer (kind=int_kind), intent(out) :: sec      ! seconds

      ! local variables

      integer (kind=int_kind) :: itime
      character(len=*),parameter :: subname='(rdate2date)'

      ! Needs to be consistent with date2rdate

      ! Since time is a real and the decimal should >= 0 and <= .86400, 
      ! add 0.001 in case decimal is zero but real has underflowed by roundoff to make
      ! sure "int" generates correct result

      itime = int(time+0.001_dbl_kind)
      year  = itime/100000
      month = mod(itime,10000) / 100
      day   = mod(itime,100)
      sec   = nint(mod(time,1.0_dbl_kind)*100000._dbl_kind)

      end subroutine rdate2date

!=======================================================================

! Set the "days per month", "days per year", etc variables for the 
! current year.
!
! authors: Craig MacLachlan, UK Met Office

      subroutine set_calendar(year)

      integer (kind=int_kind), intent(in) :: year   ! current year

      ! Internal variable
      logical (kind=log_kind) :: isleap   ! Leap year logical
      integer (kind=int_kind) :: n
      character(len=*),parameter :: subname='(set_calendar)'

      if (trim(calendar_type) == trim(ice_calendar_gregorian)) then

         isleap = .false. ! not a leap year
         if (mod(year,  4) == 0) isleap = .true.
         if (mod(year,100) == 0) isleap = .false.
         if (mod(year,400) == 0) isleap = .true.      

         if (isleap) then
            daymo = daymo366
         else
            daymo = daymo365
         endif

      elseif (trim(calendar_type) == trim(ice_calendar_360day)) then
         daymo = daymo360
      else
         daymo = daymo365
      endif

      daycal(1) = 0
      do n = 1, months_per_year
         daycal(n+1) = daycal(n) + daymo(n)
      enddo
      days_per_year=daycal(months_per_year+1)

      end subroutine set_calendar

!=======================================================================

! Compute elapsed days from year0,month0,day0 to year1,month1,day1
! Same day results in 0 elapsed days

      integer function calendar_compute_days_between(year0,month0,day0,year1,month1,day1)

      integer (kind=int_kind), intent(in) :: year0   ! start year
      integer (kind=int_kind), intent(in) :: month0  ! start month
      integer (kind=int_kind), intent(in) :: day0    ! start day
      integer (kind=int_kind), intent(in) :: year1   ! end year
      integer (kind=int_kind), intent(in) :: month1  ! end month
      integer (kind=int_kind), intent(in) :: day1    ! end day

      ! Internal variable
      logical (kind=log_kind) :: isleap   ! Leap year logical
      integer (kind=int_kind) :: nday0, nday1
      character(len=*),parameter :: subname='(calendar_compute_days_between)'

      nday0 = calendar_compute_elapsed_days(year0,month0,day0)
      nday1 = calendar_compute_elapsed_days(year1,month1,day1)

      calendar_compute_days_between = nday1 - nday0

      end function calendar_compute_days_between

!=======================================================================

! Compute elapsed days from 0000-01-01 to year1,month1,day1
! 0000-01-01 is 0 elapsed days

      integer function calendar_compute_elapsed_days(year,month,day)

      integer (kind=int_kind), intent(in) :: year   ! year
      integer (kind=int_kind), intent(in) :: month  ! month
      integer (kind=int_kind), intent(in) :: day    ! day

      ! Internal variable
      logical (kind=log_kind) :: isleap   ! Leap year logical
      integer (kind=int_kind) :: nday, n
      integer (kind=int_kind) :: daymo_local(months_per_year)
      character(len=*),parameter :: subname='(calendar_compute_elapsed_days)'

      ! use 0000-01-01 as base, year 0 is a leap year
      ! this must be implemented consistent with set_calendar

      if (year < 0 .or. month <= 0 .or. month > months_per_year .or. day <= 0) then
         write(nu_diag,*) trim(subname),' ERROR for year,month,day = ',year,month,day
         call abort_ice(subname//'ERROR: illegal date')
      endif

      ! set calendar for year
      daymo_local = daymo
      if (trim(calendar_type) == trim(ice_calendar_gregorian)) then
         isleap = .false. ! not a leap year
         if (mod(year,  4) == 0) isleap = .true.
         if (mod(year,100) == 0) isleap = .false.
         if (mod(year,400) == 0) isleap = .true.      
         if (isleap) then
            daymo_local = daymo366
         else
            daymo_local = daymo365
         endif
      endif

      ! compute days from year 0000-01-01 to year-01-01
      ! don't loop thru years for performance reasons
      if (trim(calendar_type) == trim(ice_calendar_gregorian)) then
         if (year == 0) then
            nday = 0
         else
            nday = year * 365 + 1 + (year-1)/4 - (year-1)/100 + (year-1)/400
         endif
      else
         nday = year * daycal(months_per_year+1)
      endif

      do n = 1, month-1
         nday = nday + daymo_local(n)
      enddo

      if (day <= daymo_local(month)) then
         nday = nday + day - 1
      else
         write(nu_diag,*) trim(subname),' ERROR for year,month,day = ',year,month,day
         call abort_ice(subname//'ERROR: illegal day in month')
      endif

      calendar_compute_elapsed_days = nday

      end function calendar_compute_elapsed_days

!=======================================================================

! Compute elapsed days from 0000-01-01 to year1,month1,day1
! 0000-01-01 is 0 elapsed days

      subroutine calendar_set_date_from_timesecs(ttimesecs)

      real (kind=dbl_kind), intent(in) :: ttimesecs   ! seconds since init date

      ! Internal variable
      integer (kind=int_kind) :: nday0, nday1, n
      integer (kind=int_kind) :: isecday  ! seconds per day
      real (kind=dbl_kind) :: secday
      character(len=*),parameter :: subname='(calendar_set_date_from_timesecs)'

      call icepack_query_parameters(secday_out=secday)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      timesecs = ttimesecs

      nyr = year_init
      month = month_init
      mday = day_init
      sec = ttimesecs
      call set_calendar(nyr)

      isecday = nint(secday)
      if (sec >= isecday) then
         mday = mday + int(sec/isecday)
         sec = mod(sec,isecday)
      endif
      do while (mday > daymo(month))
         mday = mday - daymo(month)
         month = month + 1
         do while (month > months_per_year)
            month = month - months_per_year
            nyr = nyr + 1
            call set_calendar(nyr)
         enddo
      enddo

      end subroutine calendar_set_date_from_timesecs

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
        integer (kind=int_kind) :: n

        if (present(ihour)) then
          !-----------------
          ! yyyy mm dd HH
          !-----------------
          iyr=iyear-1901
          dtime = floor(365.25_dbl_kind*iyr)*c1 + idd*c1 + ihour/24._dbl_kind
          if (mod(iyr,4)==3) then
            do n = 1,imm-1
               dtime = dtime + daymo366(n)
            enddo
          else
            do n = 1,imm-1
               dtime = dtime + daymo365(n)
            enddo
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
