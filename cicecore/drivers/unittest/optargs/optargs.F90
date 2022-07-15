
      program optargs

      use optargs_subs, only: computeA, computeB, computeC, computeD
      use optargs_subs, only: oa_error, oa_OK, oa_A, oa_B, oa_C, oa_D
      use optargs_subs, only: oa_layer1, oa_count1

      implicit none

      real*8  :: Ai1, Ao
      real*8  :: B
      real*8  :: Ci1, Co
      real*8  :: Di1, Di2, Do
      integer :: ierr, ierrV

      integer :: n
      integer, parameter :: ntests = 100
      integer :: iresult
      real*8  :: result, resultV
      real*8, parameter :: errtol = 1.0e-12

      !----------------------

      write(6,*) 'RunningUnitTest optargs'
      write(6,*) ' '

      iresult = 0
      do n = 1,ntests

        Ai1 = -99.; Ao  = -99.
        B   = -99.
        Ci1 = -99.; Co  = -99.
        Di1 = -99.; Di2 = -99.; Do = -99.

        ierr = oa_error
        result = -888.
        resultV = -999.

        computeA = .false.
        computeB = .false.
        computeC = .false.
        computeD = .false.

        select case (n)

! fails to compile as it should
!          case(0)
!            ierrV = oa_OK
!            call oa_layer1()

          ! test counts of present optional arguments at 2nd level
          ! result should be number of arguments
          case(1)
            result = -777.; resultV = -777.
            ierrV = 9
            call oa_count1(Ai1=Ai1,Ao=Ao,B=B,Ci1=Ci1,Co=Co,Di1=Di1,Di2=Di2,Do=Do,ierr=ierr)
          case(2)
            result = -777.; resultV = -777.
            ierrV = 9
            call oa_count1(Ai1,Ao,B,Ci1,Co,Di1,Di2,Do,ierr)
          case(3)
            result = -777.; resultV = -777.
            ierrV = 3
            call oa_count1(Ci1=Ci1,Co=Co,ierr=ierr)
          case(4)
            result = -777.; resultV = -777.
            ierrV = 5
            call oa_count1(Ci1=Ci1,Co=Co,ierr=ierr,Ao=Ao,Di1=Di1)

          ! test optional order
          case(11)
            result = -777.; resultV = -777.
            ierrV = oa_OK
            call oa_layer1(Ai1=Ai1,Ao=Ao,B=B,Ci1=Ci1,Co=Co,Di1=Di1,Di2=Di2,Do=Do,ierr=ierr)
          case(12)
            result = -777.; resultV = -777.
            ierrV = oa_OK
            call oa_layer1(Ci1=Ci1,Co=Co,ierr=ierr)
          case(13)
            result = -777.; resultV = -777.
            ierrV = oa_OK
            call oa_layer1(Ci1=Ci1,Co=Co,ierr=ierr,Ao=Ao,Di1=Di1)

          ! test optional argument checking
          case(21)
            computeA = .true.
            computeB = .true.
            computeC = .true.
            computeD = .true.
            result = -777.; resultV = -777.
            ierrV = oa_error
            ! B missing
            call oa_layer1(Ai1=Ai1,Ao=Ao,Ci1=Ci1,Co=Co,Di1=Di1,Di2=Di2,Do=Do,ierr=ierr)
          case(22)
            computeA = .true.
            computeB = .true.
            computeC = .true.
            computeD = .true.
            result = -777.; resultV = -777.
            ierrV = oa_error
            ! all optional missing
            call oa_layer1(Ci1=Ci1,Co=Co,ierr=ierr)
          case(23)
            computeA = .true.
            computeB = .true.
            computeC = .true.
            computeD = .true.
            result = -777.; resultV = -777.
            ierrV = oa_error
            ! some optional missing
            call oa_layer1(Ci1=Ci1,Co=Co,ierr=ierr,B=B,Ao=Ao,Di1=Di1)
          case(24)
            computeA = .true.
            computeB = .true.
            computeC = .true.
            computeD = .true.
            result = -777.; resultV = -777.
            ierrV = oa_error
            ! one optional missing
            call oa_layer1(Ai1=Ai1,Ao=Ao,B=B,Ci1=Ci1,Co=Co,Di1=Di1,Do=Do,ierr=ierr)

          ! test computations individually
          case(31)
            computeA = .true.
            ierrV = oa_A
            Ai1 = 5.
            resultV = 4.
            call oa_layer1(Ai1=Ai1,Ao=Ao,B=B,Ci1=Ci1,Co=Co,Di1=Di1,Di2=Di2,Do=Do,ierr=ierr)
            result = Ao
          case(32)
            computeB = .true.
            ierrV = oa_B
            B = 15.
            resultV = 20.
            call oa_layer1(ierr=ierr,Ai1=Ai1,Ao=Ao,B=B,Ci1=Ci1,Co=Co,Di1=Di1,Di2=Di2,Do=Do)
            result = B
          case(33)
            computeC = .true.
            ierrV = oa_C
            Ci1 = 7.
            resultV = 14.
            call oa_layer1(B=B,Ci1=Ci1,Co=Co,Di1=Di1,Ai1=Ai1,Ao=Ao,Di2=Di2,Do=Do,ierr=ierr)
            result = Co
          case(34)
            computeD = .true.
            ierrV = oa_D
            Di1 = 19; Di2=11.
            resultV = 30.
            call oa_layer1(Ai1=Ai1,Ao=Ao,Ci1=Ci1,Co=Co,Di1=Di1,Di2=Di2,Do=Do,B=B,ierr=ierr)
            result = Do

          ! test computations individually
          case(41)
            computeA = .true.
            computeC = .true.
            ierrV = oa_A + oa_C
            Ai1 = 6.
            Ci1 = 8.
            resultV = 21.
            call oa_layer1(Ai1=Ai1,Ao=Ao,B=B,Ci1=Ci1,Co=Co,Di1=Di1,Di2=Di2,Do=Do,ierr=ierr)
            result = Ao + Co
          case(42)
            computeB = .true.
            computeC = .true.
            ierrV = oa_B + oa_C
            B = -20.
            Ci1 = 2.
            resultV = -11.
            call oa_layer1(ierr=ierr,Ai1=Ai1,Ao=Ao,B=B,Ci1=Ci1,Co=Co,Di1=Di1,Di2=Di2,Do=Do)
            result = B + Co
          case(43)
            computeB = .true.
            computeD = .true.
            ierrV = oa_B + oa_D
            B = 4.
            Di1 = 3; Di2=19.
            resultV = 31.
            call oa_layer1(B=B,Ci1=Ci1,Co=Co,Di1=Di1,Ai1=Ai1,Ao=Ao,Di2=Di2,Do=Do,ierr=ierr)
            result = B + Do
          case(44)
            computeC = .true.
            computeD = .true.
            ierrV = oa_C + oa_D
            Ci1 = 7.
            Di1 = 6; Di2=7.
            resultV = 27.
            call oa_layer1(Ai1=Ai1,Ao=Ao,Ci1=Ci1,Co=Co,Di1=Di1,Di2=Di2,Do=Do,B=B,ierr=ierr)
            result = Co + Do
          case(45)
            computeA = .true.
            computeB = .true.
            computeC = .true.
            computeD = .true.
            ierrV = oa_A + oa_B + oa_C + oa_D
            Ai1 = 7.
            B   = 9. 
            Ci1 = 7.
            Di1 = 12; Di2=3.
            resultV = 49.
            call oa_layer1(Ao=Ao,B=B,Co=Co,Do=Do,Ai1=Ai1,Ci1=Ci1,Di1=Di1,Di2=Di2,ierr=ierr)
            result = Ao + B + Co + Do
          case(46)
            computeA = .true.
            computeB = .true.
            computeD = .true.
            ierrV = oa_A + oa_B + oa_D
            Ai1 = 10.
            B   = 11. 
            Di1 = 12; Di2=3.
            resultV = 40.
            call oa_layer1(Ao=Ao,B=B,Co=Co,Do=Do,Ai1=Ai1,Ci1=Ci1,Di1=Di1,Di2=Di2,ierr=ierr)
            result = Ao + B + Do

          case DEFAULT
            ierr = -1234

        end select

        ! skip -1234
        if (ierr /= -1234) then
          if (ierr == ierrV .and. abs(result-resultV) < errtol ) then
            write(6,101) 'PASS','optarg test',n,ierr,ierrV,result,resultV,Ao,B,Co,Do
!            write(6,101) 'PASS','optarg test',n,ierr,ierrV,result,resultV
          else
            write(6,101) 'FAIL','optarg test',n,ierr,ierrV,result,resultV,Ao,B,Co,Do
!            write(6,101) 'FAIL','optarg test',n,ierr,ierrV,result,resultV
            iresult = 1
          endif
        endif

      enddo

 101  format(1x,a,1x,a,1x,i2.2,2i6,3x,6g11.4)

      write(6,*) ' '
      write(6,*) 'optargs COMPLETED SUCCESSFULLY'
      if (iresult == 1) then
        write(6,*) 'optargs TEST FAILED'
      else
        write(6,*) 'optargs TEST COMPLETED SUCCESSFULLY'
      endif

      !----------------------

      end program

