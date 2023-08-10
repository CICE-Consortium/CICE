!===============================================================================
! Copyright (C) 2023, Intel Corporation
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!===============================================================================
module myomp
  use ice_kinds_mod
  !- directives ----------------------------------------------------------------
  implicit none
  private

!- interfaces ----------------------------------------------------------------
  interface domp_get_domain
    module procedure domp_get_domain_rlu
  end interface

  integer(int_kind), private :: domp_iam, domp_nt

  private :: domp_get_domain_rlu

  !- public vars & methods -----------------------------------------------------
  public  :: domp_init, domp_get_domain, domp_get_thread_no

#if defined (_OPENMP)
  ! Please note, this constant will create a compiler info for a constant
  ! expression in IF statements:
  real(kind=dbl_kind), private :: rdomp_iam, rdomp_nt
!$OMP THREADPRIVATE(domp_iam,domp_nt,rdomp_iam,rdomp_nt)
#endif

contains

! ----------------------------------------------------------------------------

  subroutine domp_init(nt_out)

#if defined (_OPENMP)
    use omp_lib, only : omp_get_thread_num, omp_get_num_threads
#endif

    !- argument(s) -------------------------------------------------------------
    integer(int_kind), intent(out) :: nt_out

!$OMP PARALLEL DEFAULT(none)
#if defined (_OPENMP)
    domp_iam  = omp_get_thread_num()
    rdomp_iam = real(domp_iam,8)
    domp_nt   = omp_get_num_threads()
    rdomp_nt  = real(domp_nt,8)
#else
    domp_iam  = 0
    domp_nt   = 1
#endif
!$OMP END PARALLEL

#if defined (_OPENMP)
    write(*,'(a26)') 'Build with openMP support'
#endif

#ifdef _OPENMP_TARGET
    write(*,'(a26)') 'Build with openMP OFFLOAD support'
#endif

    !- echo #threads:
    if (domp_nt > 1) then
      write(*,'(a20,i5,a8)') 'Running openMP with ', domp_nt, ' threads'
    else
#if defined (_OPENMP)
      write(*,'(a35)') 'Running openMP with a single thread'
#else
      write(*,'(a22)') 'Running without openMP'
#endif
    endif

    !- return value of #threads:
    nt_out = domp_nt

  end subroutine domp_init

  ! ----------------------------------------------------------------------------

  subroutine domp_get_domain_rlu(lower,upper,d_lower,d_upper)
  use ice_constants, only : p5
#if defined (_OPENMP)
    use omp_lib,   only : omp_in_parallel
#endif

    !- arguments ---------------------------------------------------------------
    integer(int_kind), intent(in)  :: lower,upper
    integer(int_kind), intent(out) :: d_lower,d_upper

#if defined (_OPENMP)
    !  local variables ---------------------------------------------------------
    real(kind=dbl_kind)    :: dlen
#endif

    ! proper action in "null" cases:
    if (upper <= 0 .or. upper < lower) then
      d_lower = 0
      d_upper = -1
      return
    endif

    ! proper action in serial sections
    d_lower = lower
    d_upper = upper

#if defined (_OPENMP)
    if (omp_in_parallel()) then
      dlen    = real(upper-lower+1, 8)
      d_lower = lower    + floor((rdomp_iam*dlen+p5)/rdomp_nt, 4)
      d_upper = lower -1 + floor((rdomp_iam*dlen+dlen+p5)/rdomp_nt, 4)
    endif
#endif

!  if (.false.) then
!    write(*,'(a14,i3,a24,i10,i10)') 'openMP thread ', domp_iam,               &
!         ' handles range: ', d_lower, d_upper
!  endif

  end subroutine domp_get_domain_rlu

  subroutine domp_get_thread_no (tnum)
    implicit none

    integer(int_kind), intent(out) :: tnum

    tnum = domp_iam + 1

  end subroutine domp_get_thread_no
  ! ----------------------------------------------------------------------------
end module myomp
!===============================================================================

