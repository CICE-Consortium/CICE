!  SVN:$Id: ice_history_write.F90 567 2013-01-07 02:57:36Z eclare $
!=======================================================================
!
! Writes history in binary format
!
! authors Tony Craig and Bruce Briegleb, NCAR
!         Elizabeth C. Hunke and William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! 2004 WHL: Block structure added 
! 2006 ECH: Accepted some CCSM code into mainstream CICE
!           Added ice_present, aicen, vicen; removed aice1...10, vice1...1.
!           Added histfreq_n and histfreq='h' options, removed histfreq='w'
!           Converted to free source form (F90)
!           Added option for binary output instead of netCDF
! 2009 D Bailey and ECH: Generalized for multiple frequency output
! 2010 Alison McLaren and ECH: Added 3D capability
! 2013 ECH split from ice_history.F90

      module ice_history_write

      implicit none
      private
      public :: ice_write_hist
      save
      
!=======================================================================

      contains

!=======================================================================
!
! write binary history file
!
! This routine writes fewer grid variables compared with the netcdf
! version, to reduce file size.  Grid variables can be obtained from
! the original grid input files.
!
! authors:   E.C.Hunke, LANL

      subroutine ice_write_hist(ns)

      use ice_kinds_mod
      use ice_calendar, only: write_ic, dayyr, histfreq, use_leap_years
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: spval
      use ice_domain_size, only: nx_global, ny_global, max_nstrm
      use ice_read_write, only: ice_open, ice_write
      use ice_fileunits, only: nu_history, nu_hdr, nu_diag
      use ice_grid, only: tarea
      use ice_history_shared
      use ice_restart_shared, only: lenstr, runid

      integer (kind=int_kind), intent(in) :: ns

      ! local variables

      integer (kind=int_kind) :: i,k,n,nn,nrec,nbits
      character (char_len) :: title
      character (char_len_long) :: ncfile(max_nstrm), hdrfile

      integer (kind=int_kind) :: icategory,i_aice

      character (char_len) :: current_date,current_time
      character (len=16) :: c_aice
      logical (kind=log_kind) :: diag

      diag = .false.

      if (my_task == master_task) then

        call construct_filename(ncfile(ns),'da',ns)

        ! add local directory path name to ncfile
        if (write_ic) then
          ncfile(ns) = trim(incond_dir)//ncfile(ns)
        else
          ncfile(ns) = trim(history_dir)//ncfile(ns)
        endif
        hdrfile = trim(ncfile(ns))//'.hdr'

        !-----------------------------------------------------------------
        ! create history files
        !-----------------------------------------------------------------
        nbits = 32 ! single precision
        call ice_open(nu_history, ncfile(ns), nbits) ! direct access
        open(nu_hdr,file=hdrfile,form='formatted',status='unknown') ! ascii

        title  = 'sea ice model: CICE'
        write (nu_hdr, 999) 'source',title,' '

        write (nu_hdr, 999) 'file name contains model date',trim(ncfile(ns)),' '
#ifdef CCSMCOUPLED
        write (nu_hdr, 999) 'runid',runid,' '
#endif
        if (use_leap_years) then
           write (nu_hdr, 999) 'calendar','Gregorian',' '
           write (title,'(a,i3,a)') 'This year has ',int(dayyr),' days'
        else
           write (nu_hdr, 999) 'calendar','noleap',' '
           write (title,'(a,i3,a)') 'All years have exactly ',int(dayyr),' days'
        end if
        write (nu_hdr, 999) 'comment',title,' '
        write (nu_hdr, 999) 'conventions','CICE',' '
        write (nu_hdr, 997) 'missing_value',spval
        write (nu_hdr, 997) '_FillValue',spval

        call date_and_time(date=current_date, time=current_time)
        write (nu_hdr,1000) current_date(1:4), current_date(5:6), &
                            current_date(7:8), current_time(1:2), &
                            current_time(3:4), current_time(5:8)
        write (nu_hdr, *  ) ' '
        write (nu_hdr, *  ) 'Grid size:'
        write (nu_hdr, 998) '  ni',nx_global
        write (nu_hdr, 998) '  nj',ny_global
        write (nu_hdr, 998) '  nk',nzilyr
        write (nu_hdr, 998) '  nc',ncat_hist

        write (nu_hdr, *  ) 'Grid variables: (left column = nrec)'
        nrec = 1
        write (nu_hdr, 996) nrec,'tarea','area of T grid cells','m^2'
        write (nu_hdr, *  ) 'History variables: (left column = nrec)'
      endif  ! my_task = master_task
      call ice_write(nu_history, nrec, tarea, 'rda4', diag)

      do n=1,num_avail_hist_fields_2D
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) then

          nrec = nrec + 1
          if (my_task == master_task) then
            write (nu_hdr, 996) nrec,trim(avail_hist_fields(n)%vname), &
               trim(avail_hist_fields(n)%vdesc),trim(avail_hist_fields(n)%vunit)

            ! Append ice thickness range to aicen comments
            c_aice = TRIM(avail_hist_fields(n)%vname)
            i_aice = lenstr(c_aice)
            if (i_aice > 4 .and. c_aice(1:5) == 'aicen') then
              read(c_aice(6:9), '(i3)') icategory
!             avail_hist_fields(n)%vcomment = &
!                'Ice range: '//c_hi_range(icategory)
            endif
            write (nu_hdr, 995) nrec,trim(avail_hist_fields(n)%vname), &
               trim(avail_hist_fields(n)%vcomment)

            if (histfreq(ns) == '1' .or. .not. hist_avg         &
                .or. n==n_divu(ns)      .or. n==n_shear(ns)     &  ! snapshots
                .or. n==n_sig1(ns)      .or. n==n_sig2(ns)      &
                .or. n==n_trsig(ns)                             &
                .or. n==n_mlt_onset(ns) .or. n==n_frz_onset(ns) &
                .or. n==n_hisnap(ns)    .or. n==n_aisnap(ns)) then
               write (nu_hdr, 996) nrec,trim(avail_hist_fields(n)%vname), &
                  'time_rep','instantaneous'
            else
               write (nu_hdr, 996) nrec,trim(avail_hist_fields(n)%vname), &
                  'time_rep','averaged'
            endif
          endif

          call ice_write(nu_history, nrec, a2D(:,:,n,:), 'rda4', diag)

        endif
      enddo ! num_avail_hist_fields_2D

      do n = n2D + 1, n3Dccum
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) then

          do nn = 1, ncat_hist
          nrec = nrec + 1
          if (my_task == master_task) then
            write (nu_hdr, 994) nrec,trim(avail_hist_fields(n)%vname), &
               trim(avail_hist_fields(n)%vdesc),trim(avail_hist_fields(n)%vunit),nn

            if (histfreq(ns) == '1' .or. .not. hist_avg) then
               write (nu_hdr, 996) nrec,trim(avail_hist_fields(n)%vname), &
                  'time_rep','instantaneous'
            else
               write (nu_hdr, 996) nrec,trim(avail_hist_fields(n)%vname), &
                  'time_rep','averaged'
            endif
          endif

          call ice_write(nu_history, nrec, a3Dc(:,:,nn,n-n2D,:), 'rda4', diag)
          enddo ! ncat

        endif
      enddo ! num_avail_hist_fields_3Dc

      do n = n3Dccum + 1, n3Dzcum
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) then

          do k = 1, nzilyr
          nrec = nrec + 1
          if (my_task == master_task) then
            write (nu_hdr, 993) nrec,trim(avail_hist_fields(n)%vname), &
               trim(avail_hist_fields(n)%vdesc),trim(avail_hist_fields(n)%vunit),k

            if (histfreq(ns) == '1' .or. .not. hist_avg) then
               write (nu_hdr, 996) nrec,trim(avail_hist_fields(n)%vname), &
                  'time_rep','instantaneous'
            else
               write (nu_hdr, 996) nrec,trim(avail_hist_fields(n)%vname), &
                  'time_rep','averaged'
            endif
          endif

          call ice_write(nu_history, nrec, a3Dz(:,:,k,n-n3Dccum,:), 'rda4', diag)
          enddo ! nzilyr

        endif
      enddo ! num_avail_hist_fields_3Dz

      do n = n3Dzcum + 1, n4Dicum
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) then

          do nn = 1, ncat_hist
          do k = 1, nzilyr
          nrec = nrec + 1
          if (my_task == master_task) then
            write (nu_hdr, 993) nrec,trim(avail_hist_fields(n)%vname), &
               trim(avail_hist_fields(n)%vdesc),trim(avail_hist_fields(n)%vunit),nn,k

            if (histfreq(ns) == '1' .or. .not. hist_avg) then
               write (nu_hdr, 996) nrec,trim(avail_hist_fields(n)%vname), &
                  'time_rep','instantaneous'
            else
               write (nu_hdr, 996) nrec,trim(avail_hist_fields(n)%vname), &
                  'time_rep','averaged'
            endif
          endif

          call ice_write(nu_history, nrec, a4Di(:,:,k,nn,n-n3Dzcum,:), 'rda4', diag)
          enddo ! nzilyr
          enddo ! ncat_hist

        endif
      enddo ! num_avail_hist_fields_4Di

995     format(i3,2x,a,' comment: ',a)
996     format(i3,2x,a,': ',a,',',2x,a)
993     format(i3,2x,a,': ',a,',',2x,a,2x,' cat ',i3,2x,'zlvl ',i3)
994     format(i3,2x,a,': ',a,',',2x,a,2x,' cat ',i3)
997     format(a,': ',es13.6)
998     format(a,': ',i6)
999     format(a,': ',a,2x,a)
1000    format('This dataset was created on ', &
                a,'-',a,'-',a,' at ',a,':',a,':',a)

      if (my_task == master_task) then
        close (nu_hdr)     ! header file
        close (nu_history) ! data file
        write (nu_diag,*) ' '
        write (nu_diag,*) 'Finished writing ',trim(ncfile(ns))
      endif

      end subroutine ice_write_hist

!=======================================================================

      end module ice_history_write

!=======================================================================
