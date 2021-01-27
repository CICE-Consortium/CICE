!=======================================================================

! Calendar routines for managing time
!
! Authors: Elizabeth C. Hunke, LANL
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

      ! INTERFACES

      public :: init_calendar     ! initialize calendar
      public :: advance_timestep  ! advance model 1 timestep and update calendar
      public :: calendar ! update model internal calendar/time information
      public :: set_date_from_timesecs ! set model date from time in seconds
                                       ! (relative to init date)
                                       ! needed for binary restarts
      public :: hc_jday  ! converts "calendar" date to HYCOM julian day

      ! semi-private, only used directly by unit tester
      public :: compute_elapsed_days ! compute elapsed days since 0000-01-01
      public :: compute_days_between ! compute elapsed days between two dates
      public :: update_date          ! input date and delta date, compute new date
      public :: calendar_date2time   ! convert date to time relative to init date
      public :: calendar_time2date   ! convert time to date relative to init date

      ! private functions
      private :: set_calendar          ! sets model calendar type (noleap, etc)
      private :: compute_calendar_data ! compute info about calendar for a given year

      ! PUBLIC

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

      subroutine advance_timestep()

      use ice_communicate, only: my_task, master_task

      ! local variables

      integer(kind=int_kind) :: &
         idt       ! integer dt
      character(len=*),parameter :: subname='(advance_timestep)'

      istep = istep + 1
      istep1 = istep1 + 1
      idt = nint(dt)
      ! dt is historically a real but it should be an integer
      ! make sure dt is very close to an integer
      if (abs(real(idt,kind=dbl_kind)-dt) > 1.0e-7) then
         write(nu_diag,*) trim(subname),' ERROR dt error, needs to be integer number of seconds, dt=',dt
         call abort_ice(subname//'ERROR: improper dt')
      endif
      sec = sec + idt
      call calendar()

      end subroutine advance_timestep

!=======================================================================
! Update the calendar and time manager info

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
#if (1 == 0)
      real    (kind=dbl_kind) :: secday ! seconds per day
      integer (kind=int_kind) :: isecday  ! seconds per day
#endif
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

#if (1 == 0)
      call icepack_query_parameters(secday_out=secday)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)
#endif

      call update_date(nyr,month,mday,sec)
      call set_calendar(nyr)
#if (1 == 0)
      isecday = nint(secday)
      if (sec >= isecday) then
         mday = mday + int(sec/isecday)
         sec = mod(sec,isecday)
      endif
      do while (month > months_per_year)
         month = month - months_per_year
         nyr = nyr + 1
         call set_calendar(nyr)
      enddo
      do while (mday > daymo(month))
         mday = mday - daymo(month)
         month = month + 1
         do while (month > months_per_year)
            month = month - months_per_year
            nyr = nyr + 1
            call set_calendar(nyr)
         enddo
      enddo
#endif

      idate = (nyr)*10000 + month*100 + mday ! date (yyyymmdd) 
      yday = daycal(month) + mday            ! day of the year
      elapsed_months = (nyr - year_init)*months_per_year + month - month_init
      elapsed_days = compute_days_between(year_init,month_init,day_init,nyr,month,mday)
      elapsed_hours = elapsed_days * hours_per_day
#if (1 == 0)
      timesecs = real(elapsed_days,kind=dbl_kind)*secday + &
                 real(sec,kind=dbl_kind)
#endif
      call calendar_date2time(nyr,month,mday,sec,timesecs)

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
! Set the model calendar data for year

      subroutine set_calendar(year)

      integer (kind=int_kind), intent(in) :: year   ! current year

      ! Internal variable
#if (1 == 0)
      logical (kind=log_kind) :: isleap   ! Leap year logical
      integer (kind=int_kind) :: n
#endif
      character(len=*),parameter :: subname='(set_calendar)'

      call compute_calendar_data(year,daymo,daycal,dayyr)

#if (1 == 0)
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
      dayyr=daycal(months_per_year+1)
#endif

      end subroutine set_calendar

!=======================================================================
! Add and reconcile date
! delta time arguments are optional

      subroutine update_date(anyr,amon,aday,asec,dnyr,dmon,dday,dsec)

      use ice_communicate, only: my_task, master_task

      integer (kind=int_kind), intent(inout) :: anyr, amon, aday, asec  ! year, month, day, sec
      integer (kind=int_kind), intent(in), optional :: dnyr, dmon, dday, dsec  ! delta year, month, day, sec

      ! local variables
      integer (kind=int_kind) :: tdaymo (months_per_year)   ! days per month
      integer (kind=int_kind) :: tdaycal(months_per_year+1) ! day count per month
      integer (kind=int_kind) :: tdayyr                     ! days in year
      real    (kind=dbl_kind) :: secday ! seconds per day
      integer (kind=int_kind) :: isecday  ! seconds per day
      integer (kind=int_kind) :: delta
      character(len=*),parameter :: subname='(update_date)'

      call icepack_query_parameters(secday_out=secday)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)
      isecday = nint(secday)

      ! order matters.  think about adding 1 month and 10 days to the 25th of a month
      ! what is the right order?
      ! will add all deltas then reconcile years then months then days then seconds

      if (present(dnyr)) anyr = anyr + dnyr
      if (present(dmon)) amon = amon + dmon
      if (present(dday)) aday = aday + dday
      if (present(dsec)) asec = asec + dsec

      ! adjust negative data first
      ! reconcile months - years
      do while (amon <= 0)
         delta = int((abs(amon))/months_per_year) + 1
         anyr = anyr - delta
         amon = amon + delta*months_per_year
      enddo
      call compute_calendar_data(anyr,tdaymo,tdaycal,tdayyr)

      ! reconcile days - months - years
      do while (aday <= 0)
         amon = amon - 1
         do while (amon <= 0)
            delta = int((abs(amon))/months_per_year) + 1
            anyr = anyr - delta
            amon = amon + delta*months_per_year
            call compute_calendar_data(anyr,tdaymo,tdaycal,tdayyr)
         enddo
         aday = aday + tdaymo(amon)
      enddo

      ! reconcile seconds - days - months - years
      if (asec < 0) then
         delta = int(abs(asec)/isecday) + 1
         aday = aday - delta
         asec = asec + delta*isecday
      endif
      do while (aday <= 0)
         amon = amon - 1
         do while (amon <= 0)
            delta = int((abs(amon))/months_per_year) + 1
            anyr = anyr - delta
            amon = amon + delta*months_per_year
            call compute_calendar_data(anyr,tdaymo,tdaycal,tdayyr)
         enddo
         aday = aday + tdaymo(amon)
      enddo

      ! check for negative data
      if (anyr < 0 .or. amon <= 0 .or. aday <= 0 .or. asec < 0) then
         write(nu_diag,*) trim(subname),' ERROR in dateA, ',anyr,amon,aday,asec
         call abort_ice(subname//'ERROR: in date')
      endif

      ! reconcile months - years
      do while (amon > months_per_year)
         delta = int((amon-1)/months_per_year)
         anyr = anyr + delta
         amon = amon - delta*months_per_year
      enddo
      call compute_calendar_data(anyr,tdaymo,tdaycal,tdayyr)

      ! reconcile days - months - years
      do while (aday > tdaymo(amon))
         aday = aday - tdaymo(amon)
         amon = amon + 1
         do while (amon > months_per_year)
            delta = int((amon-1)/months_per_year)
            anyr = anyr + delta
            amon = amon - delta*months_per_year
            call compute_calendar_data(anyr,tdaymo,tdaycal,tdayyr)
         enddo
      enddo

      ! reconcile seconds - days - months - years
      if (asec >= isecday) then
         delta = int(asec/isecday)
         aday = aday + delta
         asec = asec - delta*isecday
      endif
      do while (aday > tdaymo(amon))
         aday = aday - tdaymo(amon)
         amon = amon + 1
         do while (amon > months_per_year)
            delta = int((amon-1)/months_per_year)
            anyr = anyr + delta
            amon = amon - delta*months_per_year
            call compute_calendar_data(anyr,tdaymo,tdaycal,tdayyr)
         enddo
      enddo

      ! check for negative data, just in case
      if (anyr < 0 .or. amon <= 0 .or. aday <= 0 .or. asec < 0) then
         write(nu_diag,*) trim(subname),' ERROR in dateB, ',anyr,amon,aday,asec
         call abort_ice(subname//'ERROR: in date')
      endif

      end subroutine update_date

!=======================================================================

! Set internal calendar date from timesecs input
! Needed for binary restarts where only timesecs is on the restart file

      subroutine set_date_from_timesecs(ttimesecs)

      real (kind=dbl_kind), intent(in) :: ttimesecs   ! seconds since init date

      ! Internal variable
#if (1 == 0)
      integer (kind=int_kind) :: ndays
      real (kind=dbl_kind) :: secday, rdays
#endif
      character(len=*),parameter :: subname='(set_date_from_timesecs)'

      call calendar_time2date(ttimesecs,nyr,month,mday,sec,year_init,month_init,day_init,sec_init)

#if (1 == 0)
      call icepack_query_parameters(secday_out=secday)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      timesecs = ttimesecs

      nyr = year_init
      month = 1
      mday = 1
      sec = 0
      call set_calendar(nyr)

      ! first estimate of nyr
      call set_calendar(nyr)
      rdays = ttimesecs/secday
      nyr = nyr + int(rdays)/dayyr

      ! reduce estimate of nyr if ndays > rdays
      ndays = compute_days_between(year_init,month_init,day_init,nyr,month,mday)
      if (ndays > int(rdays)) then
         nyr = nyr - (ndays - int(rdays))/dayyr - 1
         ndays = compute_days_between(year_init,month_init,day_init,nyr,month,mday)
      endif
      call set_calendar(nyr)

      ! compute residiual, switch to integers, advance calendar
      rdays = ttimesecs/secday
      mday = int(rdays) - ndays + 1

      do while (mday > daymo(month))
         mday = mday - daymo(month)
         month = month + 1
         do while (month > months_per_year)
            month = month - months_per_year
            nyr = nyr + 1
            call set_calendar(nyr)
         enddo
      enddo

      ndays = compute_days_between(year_init,month_init,day_init,nyr,month,mday)
      sec = int(ttimesecs - real(ndays,kind=dbl_kind)*secday)
      if (sec > secday) then
         write(nu_diag,*) trim(subname),' ERROR in seconds, ',nyr,month,mday,sec
         call abort_ice(subname//'ERROR: in seconds')
      endif
#endif

      end subroutine set_date_from_timesecs

!=======================================================================
! Compute elapsed days from year0,month0,day0 to year1,month1,day1
! Same day results in 0 elapsed days

      integer function compute_days_between(year0,month0,day0,year1,month1,day1)

      integer (kind=int_kind), intent(in) :: year0   ! start year
      integer (kind=int_kind), intent(in) :: month0  ! start month
      integer (kind=int_kind), intent(in) :: day0    ! start day
      integer (kind=int_kind), intent(in) :: year1   ! end year
      integer (kind=int_kind), intent(in) :: month1  ! end month
      integer (kind=int_kind), intent(in) :: day1    ! end day

      ! Internal variable
      logical (kind=log_kind) :: isleap   ! Leap year logical
      integer (kind=int_kind) :: nday0, nday1
      character(len=*),parameter :: subname='(compute_days_between)'

      nday0 = compute_elapsed_days(year0,month0,day0)
      nday1 = compute_elapsed_days(year1,month1,day1)

      compute_days_between = nday1 - nday0

      end function compute_days_between

!=======================================================================
! compute calendar data based on year

      subroutine compute_calendar_data(ayear,adaymo,adaycal,adayyr)

      integer (kind=int_kind), intent(in)  :: ayear   ! year
      integer (kind=int_kind), intent(out) :: adaymo(:)  ! days per month
      integer (kind=int_kind), intent(out) :: adaycal(:) ! day count per month
      integer (kind=int_kind), intent(out) :: adayyr  ! days per year

      ! Internal variable
      logical (kind=log_kind) :: isleap   ! Leap year logical
      integer (kind=int_kind) :: n
      character(len=*),parameter :: subname='(compute_calendar_data)'

      if (ayear < 0) then
         write(nu_diag,*) trim(subname),' ERROR in ayear = ',ayear
         call abort_ice(subname//'ERROR: in ayear')
      endif

      if (size(adaymo)  /= months_per_year .or. &
          size(adaycal) /= months_per_year+1 ) then
         call abort_ice(subname//'ERROR: in argument sizes')
      endif

      if (trim(calendar_type) == trim(ice_calendar_gregorian)) then

         isleap = .false. ! not a leap year
         if (mod(ayear,  4) == 0) isleap = .true.
         if (mod(ayear,100) == 0) isleap = .false.
         if (mod(ayear,400) == 0) isleap = .true.      

         if (isleap) then
            adaymo = daymo366
         else
            adaymo = daymo365
         endif

      elseif (trim(calendar_type) == trim(ice_calendar_360day)) then
         adaymo = daymo360
      else
         adaymo = daymo365
      endif

      adaycal(1) = 0
      do n = 1, months_per_year
         adaycal(n+1) = adaycal(n) + adaymo(n)
      enddo
      adayyr=adaycal(months_per_year+1)

      end subroutine compute_calendar_data

!=======================================================================
! Compute elapsed days from 0000-01-01 to year1,month1,day1
! 0000-01-01 is 0 elapsed days

      integer function compute_elapsed_days(ayear,amonth,aday)

      integer (kind=int_kind), intent(in) :: ayear   ! year
      integer (kind=int_kind), intent(in) :: amonth  ! month
      integer (kind=int_kind), intent(in) :: aday    ! day

      ! Internal variable
      integer (kind=int_kind) :: ced_nday, n
      integer (kind=int_kind) :: lyear,lmonth,lday,lsec
      integer (kind=int_kind) :: tdaymo (months_per_year)   ! days per month
      integer (kind=int_kind) :: tdaycal(months_per_year+1) ! day count per month
      integer (kind=int_kind) :: tdayyr                     ! days in year
      character(len=*),parameter :: subname='(compute_elapsed_days)'

      ! use 0000-01-01 as base, year 0 is a leap year
      ! this must be implemented consistent with set_calendar

      lyear = ayear
      lmonth = amonth
      lday = aday
      lsec = 0

      if (lyear < 0 .or. lmonth <= 0 .or. lday <= 0) then
         write(nu_diag,*) trim(subname),' ERROR for year,month,day = ',lyear,lmonth,lday
         call abort_ice(subname//'ERROR: illegal date')
      elseif (lmonth > months_per_year) then
         call update_date(lyear,lmonth,lday,lsec)
      endif

      ! compute days from year 0000-01-01 to year-01-01
      ! don't loop thru years for performance reasons
      if (trim(calendar_type) == trim(ice_calendar_gregorian)) then
         if (lyear == 0) then
            ced_nday = 0
         else
            ced_nday = lyear * 365 + 1 + (lyear-1)/4 - (lyear-1)/100 + (lyear-1)/400
         endif
      else
         ced_nday = lyear * daycal(months_per_year+1)
      endif

      ! now compute days in this year
      call compute_calendar_data(lyear,tdaymo,tdaycal,tdayyr)

      do n = 1, lmonth-1
         ced_nday = ced_nday + tdaymo(n)
      enddo

      if (lday <= tdaymo(lmonth)) then
         ced_nday = ced_nday + lday - 1
      else
         write(nu_diag,*) trim(subname),' ERROR for year,month,day = ',ayear,amonth,aday
         call abort_ice(subname//'ERROR: illegal day in month')
      endif

      compute_elapsed_days = ced_nday

      end function compute_elapsed_days

!=======================================================================
! Compute time in seconds from input calendar date
! relative to year_init, month_init, day_init, sec_init unless _ref values passed in
! For santity, must pass all four ref values or none

      subroutine calendar_date2time(anyr,amon,aday,asec,atimesecs,nyr_ref,mon_ref,day_ref,sec_ref)

      integer(kind=int_kind), intent(in)  :: &
        anyr,amon,aday,asec                              ! year, month, day, sec of ttimesecs
      real   (kind=dbl_kind), intent(out) :: atimesecs   ! seconds since init date
      integer(kind=int_kind), intent(in), optional  :: &
        nyr_ref,mon_ref,day_ref,sec_ref                  ! year, month, day, sec reference time

      ! Internal variable
      real    (kind=dbl_kind) :: secday
      integer (kind=int_kind) :: elapsed_days ! since beginning this run
      integer (kind=int_kind) :: lnyr_ref,lmon_ref,lday_ref,lsec_ref  ! local reference year, month, day, sec
      integer (kind=int_kind) :: cnt
      character(len=*),parameter :: subname='(calendar_date2time)'

      ! set reference date and check that 0 or 4 optional arguments are passed
      cnt = 0
      if (present(nyr_ref)) then
         lnyr_ref = nyr_ref
         cnt = cnt + 1
      else
         lnyr_ref = year_init
      endif
      if (present(mon_ref)) then
         lmon_ref = mon_ref
         cnt = cnt + 1
      else
         lmon_ref = month_init
      endif
      if (present(day_ref)) then
         lday_ref = day_ref
         cnt = cnt + 1
      else
         lday_ref = day_init
      endif
      if (present(sec_ref)) then
         lsec_ref = sec_ref
         cnt = cnt + 1
      else
         lsec_ref = sec_init
      endif
      if (cnt /= 0 .and. cnt /= 4) then
         write(nu_diag,*) trim(subname),' ERROR in ref args, must pass 0 or 4 '
         call abort_ice(subname//'ERROR: in ref args, must pass 0 or 4')
      endif

      call icepack_query_parameters(secday_out=secday)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      elapsed_days = compute_days_between(lnyr_ref,lmon_ref,lday_ref,anyr,amon,aday)
      atimesecs = real(elapsed_days,kind=dbl_kind)*secday + &
                  real(asec,kind=dbl_kind) - real(lsec_ref,kind=dbl_kind)

      end subroutine calendar_date2time

!=======================================================================
! Compute calendar date from input time in seconds
! relative to year_init, month_init, day_init, sec_init or ref data if passed.
! For sanity, require all four or no ref values.
! Implemented to minimize accumulating errors and avoid overflows
! and perform well.

      subroutine calendar_time2date(atimesecs,anyr,amon,aday,asec,nyr_ref,mon_ref,day_ref,sec_ref)

      real   (kind=dbl_kind), intent(in)  :: atimesecs            ! seconds since init date
      integer(kind=int_kind), intent(out) :: &
        anyr,amon,aday,asec              ! year, month, day, sec of timesecs
      integer(kind=int_kind), intent(in), optional  :: &
        nyr_ref,mon_ref,day_ref,sec_ref  ! year, month, day, sec reference time

      ! Internal variable
      integer (kind=int_kind) :: ndays
      integer (kind=int_kind) :: tnyr, tmon, tday, tsec     ! temporaries
      integer (kind=int_kind) :: tdaymo (months_per_year)   ! days per month
      integer (kind=int_kind) :: tdaycal(months_per_year+1) ! day count per month
      integer (kind=int_kind) :: tdayyr                     ! days in year
      real (kind=dbl_kind) :: secday, rdays, ltimesecs
      integer (kind=int_kind) :: lnyr_ref,lmon_ref,lday_ref,lsec_ref  ! local reference year, month, day, sec
      integer (kind=int_kind) :: cnt
      character(len=*),parameter :: subname='(calendar_time2date)'

      call icepack_query_parameters(secday_out=secday)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      ! we could allow negative atimesecs, but this shouldn't be needed
      if (atimesecs < 0._dbl_kind) then
         write(nu_diag,*) trim(subname),' ERROR in atimesecs ',atimesecs
         call abort_ice(subname//'ERROR: in atimesecs')
      endif

      ! set reference date and check that 0 or 4 optional arguments are passed
      cnt = 0
      if (present(nyr_ref)) then
         lnyr_ref = nyr_ref
         cnt = cnt + 1
      else
         lnyr_ref = year_init
      endif
      if (present(mon_ref)) then
         lmon_ref = mon_ref
         cnt = cnt + 1
      else
         lmon_ref = month_init
      endif
      if (present(day_ref)) then
         lday_ref = day_ref
         cnt = cnt + 1
      else
         lday_ref = day_init
      endif
      if (present(sec_ref)) then
         lsec_ref = sec_ref
         cnt = cnt + 1
      else
         lsec_ref = sec_init
      endif
      if (cnt /= 0 .and. cnt /= 4) then
         write(nu_diag,*) trim(subname),' ERROR in ref args, must pass 0 or 4 '
         call abort_ice(subname//'ERROR: in ref args, must pass 0 or 4')
      endif

! -------------------------------------------------------------------
! tcraig, this is risky because atimesecs is real and could be very large
!      anyr = lnyr_ref
!      amon = lmon_ref
!      aday = lday_ref
!      asec = lsec_ref
!
!      call update_date(anyr,amon,aday,asec,dsec=nint(atimesecs))
!      return
! -------------------------------------------------------------------

      ! initial guess
      tnyr = lnyr_ref
      tmon = 1
      tday = 1
      tsec = 0

      ! add initial seconds to timesecs and treat lsec_ref as zero 
      ltimesecs = atimesecs + real(lsec_ref,kind=dbl_kind)

      ! first estimate of tnyr
      call compute_calendar_data(tnyr,tdaymo,tdaycal,tdayyr)
      rdays = ltimesecs/secday
      tnyr = tnyr + int(rdays)/tdayyr

      ! reduce estimate of tnyr if ndays > rdays
      ndays = compute_days_between(lnyr_ref,lmon_ref,lday_ref,tnyr,tmon,tday)
      if (ndays > int(rdays)) then
         tnyr = tnyr - (ndays - int(rdays))/tdayyr - 1
         ndays = compute_days_between(lnyr_ref,lmon_ref,lday_ref,tnyr,tmon,tday)
      endif
      call compute_calendar_data(tnyr,tdaymo,tdaycal,tdayyr)

      ! compute residual days, switch to integers, compute date
      rdays = ltimesecs/secday
      tday = int(rdays) - ndays + 1

      do while (tday > tdaymo(tmon))
         tday = tday - tdaymo(tmon)
         tmon = tmon + 1
         do while (tmon > months_per_year)
            tmon = tmon - months_per_year
            tnyr = tnyr + 1
            call compute_calendar_data(tnyr,tdaymo,tdaycal,tdayyr)
         enddo
      enddo

      ndays = compute_days_between(lnyr_ref,lmon_ref,lday_ref,tnyr,tmon,tday)
      tsec = int(ltimesecs - real(ndays,kind=dbl_kind)*secday)
      if (tsec > secday) then
         write(nu_diag,*) trim(subname),' ERROR in seconds, ',tnyr,tmon,tday,tsec
         call abort_ice(subname//'ERROR: in seconds')
      endif

      anyr = tnyr
      amon = tmon
      aday = tday
      asec = tsec

      end subroutine calendar_time2date

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
