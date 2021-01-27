
      program calchk
      use ice_kinds_mod, only: int_kind, dbl_kind
      use ice_calendar
      implicit none

      integer(kind=int_kind) :: nyr1,month1,mday1,sec1
      integer(kind=int_kind) :: nyr3,month3,mday3,sec3,nday3
      integer(kind=int_kind) :: nyr4,month4,mday4,sec4,nday4
      integer(kind=int_kind) :: nyr5a,month5a,mday5a,sec5a,nday5a
      integer(kind=int_kind) :: nyr5b,month5b,mday5b,sec5b,nday5b
      integer(kind=int_kind) :: nyrmax
      integer(kind=int_kind) :: nday,ndayp
      integer(kind=int_kind) :: n,nd,nf,xadd
      integer(kind=int_kind) :: yi,mi,di,si
      integer(kind=int_kind) :: dnyr,dmon,dday,dsec,mday5
      real(kind=dbl_kind) :: timesecs4
      character(len=128) :: str1,str2,xunits
      character(len=8) :: errorflagA
      character(len=8) :: errorflag1 ,errorflag2 ,errorflag3
      character(len=8) :: errorflag1A,errorflag2A,errorflag3A,errorflag4A,errorflag5A,errorflag5B
      character(len=*), parameter ::  &
         passflag = 'PASS', &
         failflag = 'FAIL'

      write(6,*) ' '
      write(6,*) 'Running CALCHK'
      write(6,*) ' '

      errorflagA = passflag
      errorflag1A = passflag
      errorflag2A = passflag
      errorflag3A = passflag
      errorflag4A = passflag
      errorflag5A = passflag
      errorflag5B = passflag

      ndtd = 1

      ! test nyrmax years from year 0
!      nyrmax = 1000
      nyrmax = 100000

    ! test 3 calendars
    do n = 1,3

      if (n == 1) then
         use_leap_years = .false.
         days_per_year = 365
         str1 = 'noleap'
      elseif (n == 2) then
         use_leap_years = .false.
         days_per_year = 360
         str1 = '360day'
      elseif (n == 3) then
         use_leap_years = .true.
         days_per_year = 365
         str1 = 'gregorian'
      endif
         
      istep0 = 1000
      year_init = 0
      month_init = 1
      day_init = 1
      sec_init = 0
      nyr = -1
      month = -1
      mday = -1
      dt = 86400._dbl_kind
      diagfreq = 99999999
      call init_calendar()

      !-----------------
      ! This test makes sure compute_elapsed_days works for different calendars
      ! and multiple years.  This also checks that the timesecs value computed
      ! in calendar and passed into set_date_from_timesecs returns the correct date.
      !-----------------

      ndayp = -1   ! prior day
      sec = 0
      do nyr = 0,nyrmax
         do month = 1,months_per_year
            do mday = 1,daymo(month)
               errorflag1 = passflag
               errorflag2 = passflag
               call calendar()
               nday = compute_elapsed_days(nyr,month,mday)
               if (nday - ndayp /= 1) then
                  errorflagA  = failflag
                  errorflag1  = failflag
                  errorflag1A = failflag
                  write(6,*) 'ERROR1: did not increment one day'
               endif
               nyr1 = nyr
               month1 = month
               mday1 = mday
               sec1 = sec
               call set_date_from_timesecs(timesecs)
               if (nyr /= nyr1 .or. month /= month1 .or. mday /= mday1 .or. sec /= sec1) then
                  errorflagA  = failflag
                  errorflag2  = failflag
                  errorflag2A = failflag
                  write(6,*) 'ERROR2: timesecs error'
                  write(6,1001) 'e2',nday,nyr1,'-',month1,'-',mday1,':',sec1,' timesecs = ',timesecs
               endif
               if (errorflag1 /= passflag .or. errorflag2 /= passflag .or. &
                   nday <= 10 .or. mod(nday,nyrmax*10) == 0 .or. &
                   (nyr == nyrmax .and. month == months_per_year)) then
                  write(6,1001) '  ',nday,nyr ,'-',month ,'-',mday ,':',sec ,' timesecs = ',timesecs
               endif
               ndayp = nday
            enddo
         enddo
      enddo

      ! check total number of days run in nyrmax years
      if (nyrmax == 1000) then
         if (n == 1) then
            nday = 365364
         elseif (n == 2) then
            nday = 360359
         elseif (n == 3) then
            nday = 365607
         endif
         if (ndayp /= nday) then
            errorflagA  = failflag
            errorflag1  = failflag
            errorflag1A = failflag
            write(6,*) 'ERROR1a: final nday incorrect', ndayp, nday
         endif
      endif

      ! check total number of days run in nyrmax years
      if (nyrmax == 100000) then
         if (n == 1) then
            nday = 36500364
         elseif (n == 2) then
            nday = 36000359
         elseif (n == 3) then
            nday = 36524615
         endif
         if (ndayp /= nday) then
            errorflagA  = failflag
            errorflag1  = failflag
            errorflag1A = failflag
            write(6,*) 'ERROR1a: final nday incorrect', ndayp, nday
         endif
      endif

      !-----------------
      ! check adding arbitrary amounts to each date unit and see if calendar reconciles properly
      !-----------------

      nyr = 1000
      month = 1
      mday = 1
      sec = 0
      call calendar()
      nday = compute_elapsed_days(nyr,month,mday)
      write(6,*) ' '
      write(6,1001) '      ',nday,nyr ,'-',month ,'-',mday ,':',sec ,' timesecs = ',timesecs
      nyr1 = nyr
      month1 = month
      mday1 = mday
      sec1 = sec
      dnyr = 0
      dmon = 0
      dday = 0
      dsec = 0
      do nd = 1,10
         do nf = 1,4
            if (nf == 1) then
               xadd = nd*nd
               xunits = 'years'
               nyr = nyr + xadd
               dnyr = dnyr + xadd
            elseif (nf == 2) then
               xadd = nd*nd
               xunits = 'months'
               month = month + xadd
               dmon = dmon + xadd
            elseif (nf == 3) then
               xadd = nd*nd*nd*nd
               xunits = 'days'
               mday = mday + xadd
               dday = dday + xadd
            elseif (nf == 4) then
               xadd = nd*nd*nd*nd*nd*nd*nd
               xunits = 'seconds'
               sec = sec + xadd
               dsec = dsec + xadd
            endif
            call calendar()
            nday = compute_elapsed_days(nyr,month,mday)
            write(6,1002) '  Add ',xadd,trim(xunits)
            write(6,1001) '      ',nday,nyr ,'-',month ,'-',mday ,':',sec ,' timesecs = ',timesecs

            ! This checks date2time and time2date leveraging the pseudo random dates
            ! plus various _init settings.  Reset the _init setting to the original
            ! values at the end of the test
            yi = year_init
            mi = month_init
            di = day_init
            si = sec_init
            year_init = nyr/2
            month_init = max(month/2,1)
            day_init = max(mday*7/8,1)
            sec_init = max(sec*7/8,1)
            nyr4 = -1
            month4 = -1
            mday4 = -1
            sec4 = -1
            call calendar_date2time(nyr,month,mday,sec,timesecs4)
            call calendar_time2date(timesecs4,nyr4,month4,mday4,sec4)
!            write(6,*) 'd2t2d ',timesecs,timesecs4
            if (nyr /= nyr4 .or. month /= month4 .or. mday /= mday4 .or. sec /= sec4) then
               errorflagA  = failflag
               errorflag4A = failflag
               write(6,*) 
               write(6,*) 'ERROR4A: date2time time2date error'
               write(6,1001) 'e4',nday,nyr4,'-',month4,'-',mday4,':',sec4,' timesecs = ',timesecs4
               write(6,1001) '  ',nday,nyr ,'-',month ,'-',mday ,':',sec ,' timesecs = ',timesecs
               write(6,*) 
            endif
            year_init = yi
            month_init = mi
            day_init = di
            sec_init = si
         enddo
      enddo

      write(6,*) ' '
      nyr5a = nyr1
      month5a = month1
      mday5a = mday1
      sec5a = sec1
      write(6,1001) '      ',-1,nyr5a ,'-',month5a ,'-',mday5a ,':',sec5a 
      write(6,1002) '  Add ',dnyr,'years'
      write(6,1002) '  Add ',dmon,'months'
      write(6,1002) '  Add ',dday,'days'
      write(6,1002) '  Add ',dsec,'seconds'
      call update_date(nyr5a,month5a,mday5a,sec5a,dnyr,dmon,dday,dsec)
      write(6,1001) '      ',-1,nyr5a ,'-',month5a ,'-',mday5a ,':',sec5a 
      nyr5b = nyr5a
      month5b = month5a
      mday5b = mday5a
      sec5b = sec5a
      write(6,1001) '      ',-1,nyr5b ,'-',month5b ,'-',mday5b ,':',sec5b 
      write(6,1002) '  Sub ',dnyr,'years'
      write(6,1002) '  Sub ',dmon,'months'
      write(6,1002) '  Sub ',dday,'days'
      write(6,1002) '  Sub ',dsec,'seconds'
      call update_date(nyr5b,month5b,mday5b,sec5b,-dnyr,-dmon,-dday,-dsec)
      write(6,1001) '      ',-1,nyr5b ,'-',month5b ,'-',mday5b ,':',sec5b 

      errorflag3A = passflag
      errorflag5A = passflag
      errorflag5B = passflag
      if (n == 1) then
         nyr3 = 1487
         month3 = 1
         mday3 = 21
         mday5 = 24
         sec3 = 22825
         nday3 = 542775
      elseif (n == 2) then
         nyr3 = 1488
         month3 = 1
         mday3 = 13
         mday5 = 13
         sec3 = 22825
         nday3 = 535692
      elseif (n == 3) then
         nyr3 = 1487
         month3 = 1
         mday3 = 5
         mday5 = 7
         sec3 = 22825
         nday3 = 543120
      endif
      if (nyr /= nyr3 .or. month /= month3 .or. mday /= mday3 .or. sec /= sec3) then
         errorflagA  = failflag
         errorflag3A = failflag
         write(6,*) 
         write(6,*) 'ERROR3A: calendar advance error'
         write(6,1001) 'e3',nday,nyr3,'-',month3,'-',mday3,':',sec3,' timesecs = ',timesecs
         write(6,1001) '  ',nday,nyr ,'-',month ,'-',mday ,':',sec ,' timesecs = ',timesecs
         write(6,*) 
      endif
      if (nyr5a /= nyr3 .or. month5a /= month3 .or. mday5a /= mday5 .or. sec5a /= sec3) then
         errorflagA  = failflag
         errorflag5A = failflag
         write(6,*) 
         write(6,*) 'ERROR5A: calendar advance error'
         write(6,1001) 'e5',nday,nyr3 ,'-',month3 ,'-',mday5 ,':',sec3 ,' timesecs = ',timesecs
         write(6,1001) '  ',nday,nyr5a,'-',month5a,'-',mday5a,':',sec5a,' timesecs = ',timesecs
         write(6,*) 
      endif
      if (nyr5b /= nyr1 .or. month5b /= month1 .or. mday5b /= mday1 .or. sec5b /= sec1) then
         errorflagA  = failflag
         errorflag5B = failflag
         write(6,*) 
         write(6,*) 'ERROR5B: calendar advance error'
         write(6,1001) 'e5',nday,nyr1 ,'-',month1 ,'-',mday1 ,':',sec1 ,' timesecs = ',timesecs
         write(6,1001) '  ',nday,nyr5b,'-',month5b,'-',mday5b,':',sec5b,' timesecs = ',timesecs
         write(6,*) 
      endif

      write(6,*) ' '
      write(6,*) trim(errorflag1A)," ... ",trim(str1)," compute_elapsed_days: ",nyrmax," years "
      write(6,*) trim(errorflag2A)," ... ",trim(str1)," set_date_from_timesecs: ",nyrmax," years "
      write(6,*) trim(errorflag3A)," ... ",trim(str1)," calendar advancing check: "
      write(6,*) trim(errorflag4A)," ... ",trim(str1)," date2time time2date check: "
      write(6,*) trim(errorflag5A)," ... ",trim(str1)," calendar update_date check: "
      write(6,*) trim(errorflag5B)," ... ",trim(str1)," calendar neg update_date check: "
      write(6,*) ' '

    enddo  ! do n

 1001 format(a,i10,1x,i7.4,a,i2.2,a,i2.2,a,i5.5,a,e23.16)
 1002 format(a,i10,1x,a)

      write(6,*) ' '
      if (errorflagA == passflag) then
         write(6,*) 'CALCHK COMPLETED SUCCESSFULLY'
      else
         write(6,*) 'CALCHK FAILED'
      endif

      end program

