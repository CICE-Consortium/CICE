
      program calchk

      use ice_kinds_mod, only: int_kind, dbl_kind
      use ice_calendar, only: myear, mmonth, mday, msec
      use ice_calendar, only: year_init, month_init, day_init, sec_init
      use ice_calendar, only: dt, ndtd, istep0, diagfreq, npt, npt_unit
      use ice_calendar, only: months_per_year, daymo, timesecs, seconds_per_day
      use ice_calendar, only: use_leap_years, days_per_year
      use ice_calendar, only: compute_elapsed_days
      use ice_calendar, only: update_date, calc_timesteps
      use ice_calendar, only: init_calendar, calendar
      use ice_calendar, only: set_date_from_timesecs
      use ice_calendar, only: calendar_date2time, calendar_time2date
      use ice_calendar, only: compute_calendar_data
      implicit none

      integer(kind=int_kind) :: yearmax
      integer(kind=int_kind) :: nday,nptc
      integer(kind=int_kind) :: n,m,ny,nm,nd,nf1,nf2,xadd,nfa,nfb,nfc,ns1,ns2
      integer(kind=int_kind) :: yi,mi,di,si
      integer(kind=int_kind) :: dyear,dmon,dday,dsec
      integer(kind=int_kind) :: fyear,fmon,fday,fsec
      character(len=32) :: calstr,unitstr,signstr
      integer (kind=int_kind) :: tdaymo (months_per_year)   ! days per month
      integer (kind=int_kind) :: tdaycal(months_per_year+1) ! day count per month
      integer (kind=int_kind) :: tdayyr                     ! days in year

      integer(kind=int_kind), parameter :: ntests = 8
      character(len=8)  :: errorflag0,errorflag(1:ntests),errorflagtmp
      character(len=32) :: testname(ntests)
      integer(kind=int_kind) :: yearv(ntests),monv(ntests),dayv(ntests),secv(ntests),ndayv(ntests) ! computed values
      integer(kind=int_kind) :: yearc(ntests),monc(ntests),dayc(ntests),secc(ntests),ndayc(ntests) ! correct results
      real(kind=dbl_kind) :: timesecsv(ntests),timesecsc(ntests)
      character(len=*), parameter ::  &
         passflag = 'PASS', &
         failflag = 'FAIL'

      write(6,*) ' '
      write(6,*) 'Running CALCHK'
      write(6,*) ' '

      errorflag0   = passflag
      errorflag(:) = passflag
      testname(:) = ''
      testname(1) = 'compute_elapsed_days'
      testname(2) = 'set_date_from_timesecs'
      testname(3) = 'calendar advance'
      testname(4) = 'date2time time2date'
      testname(5) = 'big add/sub update_date'
      testname(6) = 'small add/sub update_date'
      testname(7) = 'special checks'
      testname(8) = 'calc_timesteps'

      ndtd = 1

      ! test yearmax years from year 0
!      yearmax = 1000
      yearmax = 100000

   ! test 3 calendars
   do n = 1,3

      errorflag(:) = passflag

      if (n == 1) then
         use_leap_years = .false.
         days_per_year = 365
         calstr = 'noleap'
      elseif (n == 2) then
         use_leap_years = .false.
         days_per_year = 360
         calstr = '360day'
      elseif (n == 3) then
         use_leap_years = .true.
         days_per_year = 365
         calstr = 'gregorian'
      endif

      istep0 = 1000
      year_init = 0
      month_init = 1
      day_init = 1
      sec_init = 0
      myear = -1
      mmonth = -1
      mday = -1
      dt = 86400._dbl_kind
      diagfreq = 99999999
      call init_calendar()

      !-----------------
      ! This test makes sure compute_elapsed_days works for different calendars
      ! and multiple years.  This also checks that the timesecs value computed
      ! in calendar and passed into set_date_from_timesecs returns the correct date.
      ! In test1, nday should increment 1 day each loop and the final number
      !   of days is known for 1000 and 100000 years (precomputed)
      ! In test2, set_date_from_timesecs will reset myear, mmonth, mday, msec
      !-----------------

      ndayc(1) = -1   ! prior day
      do ny = 0,yearmax
         do nm = 1,months_per_year
            do nd = 1,daymo(nm)

               errorflagtmp = passflag
               yearv(1) = ny
               monv(1) = nm
               dayv(1) = nd
               secv(1) = 0

               ! check days increment by 1
               ndayv(1) = compute_elapsed_days(yearv(1),monv(1),dayv(1))
               if (ndayv(1) - ndayc(1) /= 1) then
                  errorflagtmp  = failflag
                  errorflag(1) = failflag
                  write(6,*) 'ERROR1: did not increment one day',yearv(1),monv(1),dayv(1),ndayv(1)
               endif

               ! establish internal date and update internal calendar including timesecs
               myear = yearv(1)
               mmonth = monv(1)
               mday = dayv(1)
               msec = secv(1)
               call calendar()
               timesecsv(1) = timesecs

               ! check set_date_from_timesecs
               yearc(2) = myear
               monc(2) = mmonth
               dayc(2) = mday
               secc(2) = msec
               timesecsc(2) = timesecs
               ndayc(2) = ndayv(1)
               myear = -1
               mmonth = -1
               mday = -1
               msec = -1
               timesecs = -1
               call set_date_from_timesecs(timesecsc(2))
               if (myear /= yearc(2) .or. mmonth /= monc(2) .or. mday /= dayc(2) .or. msec /= secc(2) .or. timesecs /= timesecsc(2)) then
                  errorflagtmp  = failflag
                  errorflag(2) = failflag
                  write(6,*) 'ERROR2: timesecs error'
                  write(6,1001) 'e2',ndayc(2),yearc(2),'-',monc(2),'-',dayc(2),':',secc(2),' timesecs = ',timesecsc(2)
               endif
               if (errorflagtmp /= passflag .or. &
                   ndayv(1) <= 10 .or. mod(ndayv(1),yearmax*10) == 0 .or. &
                   (yearv(1) == yearmax .and. monv(1) == months_per_year)) then
                  write(6,1001) ' CHECK1: ',ndayv(1),yearv(1) ,'-',monv(1),'-',dayv(1),':',secv(1) ,' timesecs = ',timesecsv(1)
               endif
               ndayc(1) = ndayv(1)
            enddo
         enddo
      enddo

      ! check total number of days run in yearmax years
      if (yearmax == 1000) then
         if (n == 1) then
            ndayc(1) = 365364
         elseif (n == 2) then
            ndayc(1) = 360359
         elseif (n == 3) then
            ndayc(1) = 365607
         endif
         if (ndayv(1) /= ndayc(1)) then
            errorflag(1) = failflag
            write(6,*) 'ERROR1a: final nday incorrect', ndayv(1), ndayc(1)
         endif
      endif

      ! check total number of days run in yearmax years
      if (yearmax == 100000) then
         if (n == 1) then
            ndayc(1) = 36500364
         elseif (n == 2) then
            ndayc(1) = 36000359
         elseif (n == 3) then
            ndayc(1) = 36524615
         endif
         if (ndayv(1) /= ndayc(1)) then
            errorflag(1) = failflag
            write(6,*) 'ERROR1a: final nday incorrect', ndayv(1), ndayc(1)
         endif
      endif

      !-----------------
      ! check adding arbitrary amounts to each date unit and see if calendar reconciles properly
      ! then subtract same arbitrary amounts in reverse order and make sure it ends at original value
      !-----------------

      yearv(1) = 1000
      monv(1) = 1
      dayv(1) = 1
      secv(1) = 0
      myear = yearv(1)
      mmonth = monv(1)
      mday = dayv(1)
      msec = secv(1)
      call calendar()
      nday = compute_elapsed_days(myear,mmonth,mday)
      dyear = 0
      dmon = 0
      dday = 0
      dsec = 0
      do nfa = 1,-1,-2
         write(6,*) ' '
         write(6,1001) ' CHECK3: ',nday,myear ,'-',mmonth ,'-',mday ,':',msec ,' timesecs = ',timesecs
         do nfb = 1,10
            do nfc = 1,4
               if (nfa == 1) then
                  nf1 = nfb
                  nf2 = nfc
                  signstr = 'Add'
               elseif (nfa == -1) then
                  nf1 = 11-nfb
                  nf2 = 5-nfc
                  signstr = 'Sub'
               endif
               fyear = 0
               fmon = 0
               fday = 0
               fsec = 0
               if (nf2 == 1) then
                  xadd = nf1*nf1
                  unitstr = 'years'
                  myear = myear + nfa*xadd
                  if (nfa == 1) dyear = dyear + nfa*xadd
                  fyear = nfa*xadd
               elseif (nf2 == 2) then
                  xadd = nf1*nf1
                  unitstr = 'months'
                  mmonth = mmonth + nfa*xadd
                  if (nfa == 1) dmon = dmon + nfa*xadd
                  fmon = nfa*xadd
               elseif (nf2 == 3) then
                  xadd = nf1*nf1*nf1*nf1
                  unitstr = 'days'
                  mday = mday + nfa*xadd
                  if (nfa == 1) dday = dday + nfa*xadd
                  fday = nfa*xadd
               elseif (nf2 == 4) then
                  xadd = nf1*nf1*nf1*nf1*nf1*nf1*nf1
                  unitstr = 'seconds'
                  msec = msec + nfa*xadd
                  if (nfa == 1) dsec = dsec + nfa*xadd
                  fsec = nfa*xadd
               endif
               call calendar()
               nday = compute_elapsed_days(myear,mmonth,mday)
               write(6,1002) ' CHECK3: '//trim(signstr)//' ',xadd,trim(unitstr)
               write(6,1001) ' CHECK3: ',nday,myear ,'-',mmonth ,'-',mday ,':',msec ,' timesecs = ',timesecs

               !-----------------
               ! This checks update_date add and subtract to make sure the original value is returned
               !-----------------

               yearc(6) = myear
               monc(6) = mmonth
               dayc(6) = mday
               secc(6) = msec
               timesecsc(6) = timesecs
               yearv(6) = yearc(6)
               monv(6) = monc(6)
               dayv(6) = dayc(6)
               secv(6) = secc(6)
               call update_date(yearv(6),monv(6),dayv(6),secv(6),fyear,fmon,fday,fsec)
               write(6,1001) ' CHECK6: ',-1,yearv(6),'-',monv(6),'-',dayv(6),':',secv(6)
               if (yearc(6) == yearv(6) .and. monc(6) == monv(6) .and. dayc(6) == dayv(6) .and. secc(6) == secv(6) .and. timesecsc(6) == timesecsv(6)) then
                  errorflag(6) = failflag
                  write(6,*) ' '
                  write(6,*) 'ERROR6a: update date error'
                  write(6,1001) 'e6',nday,yearv(6),'-',monv(6),'-',dayv(6),':',secv(6),' timesecs = ',timesecsv(6)
                  write(6,1001) '  ',nday,yearc(6),'-',monc(6),'-',dayc(6),':',secc(6),' timesecs = ',timesecsc(6)
                  write(6,*)    '  ',fyear,fmon,fday,fsec
                  write(6,*) ' '
               endif
               call update_date(yearv(6),monv(6),dayv(6),secv(6),-fyear,-fmon,-fday,-fsec)
               call calendar_date2time(yearc(6),monc(6),dayc(6),secc(6),timesecsv(6))
               if (yearc(6) /= yearv(6) .or. monc(6) /= monv(6) .or. dayc(6) /= dayv(6) .or. secc(6) /= secv(6) .or. timesecsc(6) /= timesecsv(6)) then
                  errorflag(6) = failflag
                  write(6,*) ' '
                  write(6,*) 'ERROR6b: update date error'
                  write(6,1001) 'e6',nday,yearv(6),'-',monv(6),'-',dayv(6),':',secv(6),' timesecs = ',timesecsv(6)
                  write(6,1001) '  ',nday,yearc(6),'-',monc(6),'-',dayc(6),':',secc(6),' timesecs = ',timesecsc(6)
                  write(6,*)    '  ',fyear,fmon,fday,fsec
                  write(6,*) ' '
               endif

               !-----------------
               ! This checks date2time and time2date leveraging the pseudo random dates
               ! plus various reference settings.  Different reference dates means
               ! timesecs won't match, so don't check them.
               !-----------------

               yi = myear/2
               mi = max(mmonth/2,1)
               di = max(mday*7/8,1)
               si = max(msec*7/8,1)
               yearc(4) = myear
               monc(4) = mmonth
               dayc(4) = mday
               secc(4) = msec
               timesecsc(4) = timesecs
               yearv(4) = -1
               monv(4) = -1
               dayv(4) = -1
               secv(4) = -1
               timesecsv(4) = -1
               call calendar_date2time(yearc(4),monc(4),dayc(4),secc(4),timesecsv(4),yi,mi,di,si)
               call calendar_time2date(timesecsv(4),yearv(4),monv(4),dayv(4),secv(4),yi,mi,di,si)
               write(6,*) 'CHECK4: ',timesecsv(4)
               if (yearc(4) /= yearv(4) .or. monc(4) /= monv(4) .or. dayc(4) /= dayv(4) .or. secc(4) /= secv(4)) then
                  errorflag(4) = failflag
                  write(6,*) ' '
                  write(6,*) 'ERROR4: date2time time2date error'
                  write(6,1001) 'e4',nday,yearv(4),'-',monv(4),'-',dayv(4),':',secv(4),' timesecs = ',timesecsv(4)
                  write(6,1001) '  ',nday,yearc(4),'-',monc(4),'-',dayc(4),':',secc(4),' timesecs = ',timesecsc(4)
                  write(6,*) ' '
               endif

            enddo
         enddo

         yearv(3) = myear
         monv(3) = mmonth
         dayv(3) = mday
         secv(3) = msec
         timesecsv(3) = timesecs
         if (nfa == 1) then
            if (n == 1) then
               yearc(3) = 1487
               monc(3) = 1
               dayc(3) = 21
               secc(3) = 22825
               ndayc(3) = 542775
            elseif (n == 2) then
               yearc(3) = 1488
               monc(3) = 1
               dayc(3) = 13
               secc(3) = 22825
               ndayc(3) = 535692
            elseif (n == 3) then
               yearc(3) = 1487
               monc(3) = 1
               dayc(3) = 5
               secc(3) = 22825
               ndayc(3) = 543120
            endif
         elseif (nfa == -1) then
            yearc(3) = yearv(1)
            monc(3) = monv(1)
            dayc(3) = dayv(1)
            secc(3) = secv(1)
            if (n == 1) then
               ndayc(3) = 365000
            elseif (n == 2) then
               ndayc(3) = 360000
            elseif (n == 3) then
               ndayc(3) = 365243
            endif
         endif

         ! check answers
         if (yearv(3) /= yearc(3) .or. monv(3) /= monc(3) .or. dayv(3) /= dayc(3) .or. secv(3) /= secc(3)) then
            errorflag(3) = failflag
            write(6,*) ' '
            write(6,*) 'ERROR3: calendar advance error'
            write(6,1001) 'e3',nday,yearc(3),'-',monc(3),'-',dayc(3),':',secc(3),' timesecs = ',timesecsc(3)
            write(6,1001) '  ',nday,yearv(3),'-',monv(3),'-',dayv(3),':',secv(3),' timesecs = ',timesecsv(3)
            write(6,*) ' '
         endif
      enddo

      write(6,*) ' '
      yearv(1) = 1000
      monv(1) = 1
      dayv(1) = 1
      secv(1) = 0
      yearv(5) = yearv(1)
      monv(5) = monv(1)
      dayv(5) = dayv(1)
      secv(5) = secv(1)
      write(6,1001) ' CHECK5a: ',-1,yearv(5) ,'-',monv(5) ,'-',dayv(5) ,':',secv(5)
      write(6,1002) '  Add ',dyear,'years'
      write(6,1002) '  Add ',dmon,'months'
      write(6,1002) '  Add ',dday,'days'
      write(6,1002) '  Add ',dsec,'seconds'
      call update_date(yearv(5),monv(5),dayv(5),secv(5),dyear,dmon,dday,dsec)
      write(6,1001) ' CHECK5a: ',-1,yearv(5) ,'-',monv(5) ,'-',dayv(5) ,':',secv(5)
      write(6,*) ' '

      ! correct answers
      if (n == 1) then
         yearc(5) = 1487
         monc(5) = 1
         dayc(5) = 24
         secc(5) = 22825
         ndayc(5) = 542775
      elseif (n == 2) then
         yearc(5) = 1488
         monc(5) = 1
         dayc(5) = 13
         secc(5) = 22825
         ndayc(5) = 535692
      elseif (n == 3) then
         yearc(5) = 1487
         monc(5) = 1
         dayc(5) = 7
         secc(5) = 22825
         ndayc(5) = 543120
      endif

      ! check answers
      if (yearv(5) /= yearc(5) .or. monv(5) /= monc(5) .or. dayv(5) /= dayc(5) .or. secv(5) /= secc(5)) then
         errorflag(5) = failflag
         write(6,*) ' '
         write(6,*) 'ERROR5a: calendar advance error'
         write(6,1001) 'e5',nday,yearc(5),'-',monc(5),'-',dayc(5),':',secc(5),' timesecs = ',timesecs
         write(6,1001) '  ',nday,yearv(5),'-',monv(5),'-',dayv(5),':',secv(5),' timesecs = ',timesecs
         write(6,*) ' '
      endif

      write(6,1001) ' CHECK5b: ',-1,yearv(5) ,'-',monv(5) ,'-',dayv(5) ,':',secv(5)
      write(6,1002) '  Sub ',dyear,'years'
      write(6,1002) '  Sub ',dmon,'months'
      write(6,1002) '  Sub ',dday,'days'
      write(6,1002) '  Sub ',dsec,'seconds'
      call update_date(yearv(5),monv(5),dayv(5),secv(5),-dyear,-dmon,-dday,-dsec)
      write(6,1001) ' CHECK5b: ',-1,yearv(5) ,'-',monv(5) ,'-',dayv(5) ,':',secv(5)

      ! correct answers
      yearc(5) = yearv(1)
      monc(5) = monv(1)
      dayc(5) = dayv(1)
      secc(5) = secv(1)
      if (yearv(5) /= yearc(5) .or. monv(5) /= monc(5) .or. dayv(5) /= dayc(5) .or. secv(5) /= secc(5)) then
         errorflag(5) = failflag
         write(6,*) ' '
         write(6,*) 'ERROR5b: calendar advance error'
         write(6,1001) 'e5',nday,yearc(5),'-',monc(5),'-',dayc(5),':',secc(5),' timesecs = ',timesecs
         write(6,1001) '  ',nday,yearv(5),'-',monv(5),'-',dayv(5),':',secv(5),' timesecs = ',timesecs
         write(6,*) ' '
      endif

      !-------------------------
      ! Special checks:
      ! Add a month to the last day of each month
      ! Check date2time for seconds
      !-------------------------

      write(6,*) ' '
      do ny = 1,5
      do nm = 1, months_per_year
         if (ny == 1) yearv(7) = 1900
         if (ny == 2) yearv(7) = 1999
         if (ny == 3) yearv(7) = 2000
         if (ny == 4) yearv(7) = 2004
         if (ny == 5) yearv(7) = 2005
         call compute_calendar_data(yearv(7),tdaymo,tdaycal,tdayyr)
         monv(7) = nm
         dayv(7) = tdaymo(nm)
         secv(7) = 0
         if (tdaymo(mod(nm,months_per_year)+1) >= tdaymo(nm)) then
            monc(7) = mod(nm,months_per_year)+1
            dayc(7) = dayv(7)
         else
            monc(7) = mod(nm+1,months_per_year)+1
            dayc(7) = tdaymo(nm) - tdaymo(mod(nm,months_per_year)+1)
         endif
         yearc(7) = yearv(7)
         if (monc(7) < monv(7)) yearc(7) = yearv(7) + 1
         secc(7) = secv(7)
         call update_date(yearv(7),monv(7),dayv(7),secv(7),dmon=1)
         write(6,1001) ' CHECK7a:',1,yearv(7),'-',monv(7),'-',dayv(7),':',secv(7)
         if (yearv(7) /= yearc(7) .or. monv(7) /= monc(7) .or. dayv(7) /= dayc(7) .or. secv(7) /= secc(7)) then
            errorflag(7) = failflag
            write(6,*) ' '
            write(6,*) 'ERROR7a: add 1 month to end of month error'
            write(6,1001) 'e7',-1,yearc(7),'-',monc(7),'-',dayc(7),':',secc(7)
            write(6,1001) '  ',-1,yearv(7),'-',monv(7),'-',dayv(7),':',secv(7)
            write(6,*) ' '
         endif
      enddo
      enddo

      do ns1 = 0,seconds_per_day,seconds_per_day/4
      do ns2 = 0,seconds_per_day,seconds_per_day/4
         yearv(7) = 2002
         monv(7) = 3
         call compute_calendar_data(yearv(7),tdaymo,tdaycal,tdayyr)
         dayv(7) = tdaymo(monv(7))
         call calendar_date2time(yearv(7),monv(7),dayv(7),ns2,timesecsv(7),yearv(7),monv(7),dayv(7),ns1)
         write(6,*) 'CHECK7b:',ns1,ns2,timesecsv(7)
         if (timesecsv(7) /= ns2-ns1) then
            errorflag(7) = failflag
            write(6,*) ' '
            write(6,*) 'ERROR7b: sec diff same date error'
            write(6,*) '  ',ns1,ns2,timesecsv(7),ns2-ns1
            write(6,*) ' '
         endif
         call calendar_date2time(yearv(7),monv(7)+1,1,ns2,timesecsv(7),yearv(7),monv(7),dayv(7),ns1)
         write(6,*) 'CHECK7c:',ns1,ns2,timesecsv(7)
         if (timesecsv(7) /= ns2-ns1+seconds_per_day) then
            errorflag(7) = failflag
            write(6,*) ' '
            write(6,*) 'ERROR7c: sec diff next day error'
            write(6,*) '  ',ns1,ns2,timesecsv(7),ns2-ns1+seconds_per_day
            write(6,*) ' '
         endif
      enddo
      enddo

      !-------------------------
      ! calc_timesteps
      !-------------------------

      myear = 2000
      mmonth = 2
      mday = 1
      msec = 0
      do nf1 = 1,6
         npt = 10
         dt = 3600._dbl_kind

         if (nf1 == 1) then
            npt_unit = '1'
            nptc = 10
         endif
         if (nf1 == 2) then
            npt_unit = 's'
            npt = 36000.
            nptc = 10
         endif
         if (nf1 == 3) then
            npt_unit = 'h'
            nptc = 10
         endif
         if (nf1 == 4) then
            npt_unit = 'd'
            nptc = 240
         endif
         if (nf1 == 5) then
            npt_unit = 'm'
            if (n == 1) nptc = 7272
            if (n == 2) nptc = 7200
            if (n == 3) nptc = 7296
         endif
         if (nf1 == 6) then
            npt_unit = 'y'
            if (n == 1) nptc = 87600
            if (n == 2) nptc = 86400
            if (n == 3) nptc = 87672
         endif
         call calc_timesteps()
         write(6,*) 'CHECK8:',npt
         if (npt /= nptc) then
            errorflag(8) = failflag
            write(6,*) 'ERROR8: npt error',npt,nptc
         endif
      enddo

      !-------------------------
      ! write test results
      !-------------------------

      write(6,*) ' '
      write(6,*) 'Test Results: ',yearmax,' years'
      do m = 1,ntests
         write(6,*) trim(errorflag(m))," ... ",trim(calstr)," ",trim(testname(m))
         if (errorflag(m) == failflag) errorflag0=failflag
      enddo
      write(6,*) ' '

   enddo  ! do n

 1001 format(a,i10,1x,i7.4,a,i2.2,a,i2.2,a,i5.5,a,e23.16)
 1002 format(a,i10,1x,a)

      write(6,*) ' '
      if (errorflag0 == passflag) then
         write(6,*) 'CALCHK COMPLETED SUCCESSFULLY'
      else
         write(6,*) 'CALCHK FAILED'
      endif

      end program

