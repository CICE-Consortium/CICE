
      program optargs

      use optargs_subs, only: dp
      use optargs_subs, only: computeA, computeB, computeC, computeD, computeE
      use optargs_subs, only: oa_error, oa_OK, oa_A, oa_B, oa_C, oa_D, oa_E
      use optargs_subs, only: oa_layer1, oa_count1

      implicit none

      real(dp):: Ai1, Ao
      real(dp):: B
      real(dp):: Ci1, Co
      real(dp):: Di1, Di2, Do
      real(dp), allocatable :: Ei(:),Eo(:)
      integer :: ierr, ierrV

      integer :: n
      integer, parameter :: ntests = 100
      integer :: iresult
      real(dp):: result, resultV
      real(dp), parameter :: dpic = -99._dp
      real(dp), parameter :: errtol = 1.0e-12

      !----------------------

      write(6,*) 'RunningUnitTest optargs'
      write(6,*) ' '

      allocate(Ei(3),Eo(3))

      iresult = 0
      do n = 1,ntests

        Ai1 = dpic; Ao  = dpic
        B   = dpic
        Ci1 = dpic; Co  = dpic
        Di1 = dpic; Di2 = dpic; Do = dpic
        Ei  = dpic; Eo  = dpic

        ierr = oa_error
        result = -888._dp
        resultV = -999._dp

        computeA = .false.
        computeB = .false.
        computeC = .false.
        computeD = .false.
        computeE = .false.

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
            ierrV = 11
            call oa_count1(Ai1,Ao,B,Ci1,Co,Di1,Di2,Do,Ei,Eo,ierr)
          case(3)
            result = -777.; resultV = -777.
            ierrV = 3
            call oa_count1(Ci1=Ci1,Co=Co,ierr=ierr)
          case(4)
            result = -777.; resultV = -777.
            ierrV = 5
            call oa_count1(Ci1=Ci1,Co=Co,ierr=ierr,Ao=Ao,Di1=Di1)
          case(5)
            result = -777.; resultV = -777.
            ierrV = 8
            call oa_count1(Ai1,Ao,B,Ci1,Co,Di1,Di2,ierr=ierr)

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
          case(14)
            result = -777.; resultV = -777.
            ierrV = oa_OK
            call oa_layer1(Eo=Eo,Ei=Ei,Ci1=Ci1,Co=Co,ierr=ierr,Ao=Ao,Di1=Di1)

          ! test optional argument checking
          case(21)
            computeA = .true.
            computeB = .true.
            computeC = .true.
            computeD = .true.
            computeE = .true.
            result = -777.; resultV = -777.
            ierrV = oa_error
            ! B missing
            call oa_layer1(Ai1=Ai1,Ao=Ao,Ci1=Ci1,Co=Co,Di1=Di1,Di2=Di2,Do=Do,Ei=Ei,Eo=Eo,ierr=ierr)
          case(22)
            computeA = .true.
            computeB = .true.
            computeC = .true.
            computeD = .true.
            computeE = .true.
            result = -777.; resultV = -777.
            ierrV = oa_error
            ! all optional missing
            call oa_layer1(Ci1=Ci1,Co=Co,ierr=ierr)
          case(23)
            computeA = .true.
            computeB = .true.
            computeC = .true.
            computeD = .true.
            computeE = .true.
            result = -777.; resultV = -777.
            ierrV = oa_error
            ! some optional missing
            call oa_layer1(Ci1=Ci1,Co=Co,Eo=Eo,ierr=ierr,B=B,Ao=Ao,Di1=Di1)
          case(24)
            computeA = .true.
            computeB = .true.
            computeC = .true.
            computeD = .true.
            computeE = .true.
            result = -777.; resultV = -777.
            ierrV = oa_error
            ! one optional missing
            call oa_layer1(Ai1=Ai1,Ao=Ao,B=B,Ci1=Ci1,Co=Co,Di1=Di1,Do=Do,Ei=Ei,Eo=Eo,ierr=ierr)
          case(25)
            computeA = .true.
            computeB = .true.
            computeC = .true.
            computeD = .true.
            computeE = .true.
            result = -777.; resultV = -777.
            ierrV = oa_error
            ! Ei missing
            call oa_layer1(Ai1=Ai1,Ao=Ao,B=B,Ci1=Ci1,Co=Co,Di1=Di1,Di2=Di2,Do=Do,Eo=Eo,ierr=ierr)

          ! test computations individually, all args
          case(31)
            computeA = .true.
            ierrV = oa_A
            Ai1 = 5.
            resultV = 4.
            call oa_layer1(Ai1=Ai1,Ao=Ao,B=B,Ci1=Ci1,Co=Co,Di1=Di1,Di2=Di2,Do=Do,Ei=Ei,Eo=Eo,ierr=ierr)
            result = Ao
          case(32)
            computeB = .true.
            ierrV = oa_B
            B = 15.
            resultV = 20.
            call oa_layer1(ierr=ierr,Ai1=Ai1,Ao=Ao,B=B,Ci1=Ci1,Co=Co,Di1=Di1,Di2=Di2,Do=Do,Ei=Ei,Eo=Eo)
            result = B
          case(33)
            computeC = .true.
            ierrV = oa_C
            Ci1 = 7.
            resultV = 14.
            call oa_layer1(B=B,Ci1=Ci1,Co=Co,Di1=Di1,Ai1=Ai1,Ao=Ao,Di2=Di2,Do=Do,ierr=ierr,Ei=Ei,Eo=Eo)
            result = Co
          case(34)
            computeD = .true.
            ierrV = oa_D
            Di1 = 19; Di2=11.
            resultV = 30.
            call oa_layer1(Ai1=Ai1,Ao=Ao,Ci1=Ci1,Co=Co,Di1=Di1,Ei=Ei,Eo=Eo,Di2=Di2,Do=Do,B=B,ierr=ierr)
            result = Do
          case(35)
            computeE = .true.
            ierrV = oa_E
            Ei = 25.
            resultV = 81.
            call oa_layer1(Ai1=Ai1,Ao=Ao,Ci1=Ci1,Co=Co,Ei=Ei,Eo=Eo,Di1=Di1,Di2=Di2,Do=Do,B=B,ierr=ierr)
            result = sum(Eo)

          ! test computations individually, min args
          case(41)
            computeA = .true.
            ierrV = oa_A
            Ai1 = 5.
            resultV = 4.
            call oa_layer1(Ao=Ao,Co=Co,Ai1=Ai1,Ci1=Ci1,ierr=ierr)
            result = Ao
          case(42)
            computeB = .true.
            ierrV = oa_B
            B = 15.
            resultV = 20.
            call oa_layer1(ierr=ierr,Ci1=Ci1,Co=Co,B=B)
            result = B
          case(43)
            computeC = .true.
            ierrV = oa_C
            Ci1 = 7.
            resultV = 14.
            call oa_layer1(Ci1=Ci1,Co=Co,ierr=ierr)
            result = Co
          case(44)
            computeD = .true.
            ierrV = oa_D
            Di1 = 19; Di2=11.
            resultV = 30.
            call oa_layer1(Ci1=Ci1,Di1=Di1,Di2=Di2,Co=Co,Do=Do,ierr=ierr)
            result = Do
          case(45)
            computeE = .true.
            ierrV = oa_E
            Ei = 25.
            resultV = 81.
            call oa_layer1(Ci1=Ci1,Co=Co,Ei=Ei,Eo=Eo,ierr=ierr)
            result = sum(Eo)

          ! test computations in groups, mix of passed arguments
          case(51)
            computeA = .true.
            computeC = .true.
            ierrV = oa_A + oa_C
            Ai1 = 6.
            Ci1 = 8.
            resultV = 21.
            call oa_layer1(Ai1=Ai1,Ao=Ao,B=B,Ci1=Ci1,Co=Co,Eo=Eo,ierr=ierr)
            result = Ao + Co
          case(52)
            computeB = .true.
            computeC = .true.
            ierrV = oa_B + oa_C
            B = -20.
            Ci1 = 2.
            resultV = -11.
            call oa_layer1(ierr=ierr,Ai1=Ai1,Ao=Ao,B=B,Ci1=Ci1,Co=Co,Di1=Di1,Di2=Di2,Do=Do)
            result = B + Co
          case(53)
            computeB = .true.
            computeD = .true.
            ierrV = oa_B + oa_D
            B = 4.
            Di1 = 3; Di2=19.
            resultV = 31.
            call oa_layer1(B=B,Ci1=Ci1,Co=Co,Di1=Di1,Ai1=Ai1,Ao=Ao,Di2=Di2,Do=Do,ierr=ierr)
            result = B + Do
          case(54)
            computeC = .true.
            computeD = .true.
            ierrV = oa_C + oa_D
            Ci1 = 7.
            Di1 = 6; Di2=7.
            resultV = 27.
            call oa_layer1(Ai1=Ai1,Ao=Ao,Ci1=Ci1,Co=Co,Di1=Di1,Di2=Di2,Do=Do,B=B,ierr=ierr)
            result = Co + Do
          case(55)
            computeA = .true.
            computeB = .true.
            computeC = .true.
            computeD = .true.
            computeE = .true.
            ierrV = oa_A + oa_B + oa_C + oa_D + oa_E
            Ai1 = 7.
            B   = 9.
            Ci1 = 7.
            Di1 = 12; Di2=3.
            Ei  = 5
            resultV = 70.
            call oa_layer1(Ao=Ao,B=B,Co=Co,Do=Do,Ai1=Ai1,Ci1=Ci1,Di1=Di1,Di2=Di2,Ei=Ei,Eo=Eo,ierr=ierr)
            result = Ao + B + Co + Do + sum(Eo)
          case(56)
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
          case(57)
            computeB = .true.
            computeE = .true.
            ierrV = oa_B + oa_E
            B = 4.
            Ei = 8.
            resultV = 39.
            call oa_layer1(B=B,Ci1=Ci1,Co=Co,Di1=Di1,Ai1=Ai1,Ao=Ao,Di2=Di2,Do=Do,Ei=Ei,Eo=Eo,ierr=ierr)
            result = B + sum(Eo)

          case DEFAULT
            ierr = -1234

        end select

        ! skip -1234
        if (ierr /= -1234) then
          if (ierr == ierrV .and. abs(result-resultV) < errtol ) then
            write(6,101) 'PASS','optarg test',n,ierr,ierrV,result,resultV,Ao,B,Co,Do,sum(Eo)
!            write(6,101) 'PASS','optarg test',n,ierr,ierrV,result,resultV
          else
            write(6,101) 'FAIL','optarg test',n,ierr,ierrV,result,resultV,Ao,B,Co,Do,sum(Eo)
!            write(6,101) 'FAIL','optarg test',n,ierr,ierrV,result,resultV
            iresult = 1
          endif
        endif

      enddo

 101  format(1x,a,1x,a,1x,i2.2,2i6,3x,8g11.4)

      write(6,*) ' '
      write(6,*) 'optargs COMPLETED SUCCESSFULLY'
      if (iresult == 1) then
        write(6,*) 'optargs TEST FAILED'
      else
        write(6,*) 'optargs TEST COMPLETED SUCCESSFULLY'
      endif

      !----------------------

      end program

