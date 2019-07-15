!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------

!** fgmres - Flexible generalized minimum residual method (with restarts).
!
      subroutine sol_fgmres2d (solution, matvec, rhs_b, tolerance, maxinner, maxouter, nbiter, conv, level)
      use dyn_fisl_options
      use glb_ld
      use ldnh
      use prec
      use sol
      use HORgrid_options, only: Grd_yinyang_L

      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: level

      ! Initial guess on input, approximate solution on output
      real*8, dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy), intent(inout) :: solution

      ! A matrix-vector product routine (A.*v).
      interface
         subroutine matvec(v, prod, level)
            use ldnh, only: ldnh_minx, ldnh_maxx, ldnh_miny, ldnh_maxy
            implicit none
            integer, intent(in) :: level
            real*8, dimension (ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy), intent(in) :: v
            real*8, dimension (ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy), intent(out) :: prod
         end subroutine
      end interface

      ! Right hand side of the linear system.
      real*8, dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy), intent(in) :: rhs_b

      ! Tolerance to achieve. The algorithm terminates when the relative
      ! residual is below tolerance.
      real*8, intent(in) :: tolerance

      ! Restarts the method every maxinner inner iterations.
      integer, intent(in) :: maxinner

      ! Specifies the maximum number of outer iterations.
      ! Iteration will stop after maxinner*maxouter steps
      ! even if the specified tolerance has not been achieved.
      integer, intent(in) :: maxouter

      ! Total number of iteration performed
      integer, intent(out) :: nbiter

      real*8, intent(out) :: conv

      ! Author
      !     Stéphane Gaudreault, Abdessamad Qaddouri -- March 2018
      !
      ! References
      ! C. T. Kelley, Iterative Methods for Linear and Nonlinear Equations, SIAM, 1995
      ! (https://www.siam.org/books/textbooks/fr16_book.pdf)
      !
      ! Y. Saad, Iterative Methods for Sparse Linear Systems. SIAM, 2003.
      ! (http://www-users.cs.umn.edu/~saad/IterMethBook_2ndEd.pdf)
      !
      integer :: i, j, k, k1, ii, jj, ierr

      integer :: initer, outiter, nextit, it
      real*8 :: relative_tolerance, r0
      real*8 :: norm_residual, nu, t
      real*8, dimension(maxinner+1, maxinner) :: hessenberg

      real*8, dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy) :: work_space
      real*8, dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, maxinner) :: ww
      real*8, dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, maxinner+1) :: vv

      real*8 :: local_dot, dotprod
      real*8, dimension(maxinner) :: dotprod_local

      real*8, dimension(maxinner+1) :: rot_cos, rot_sin, gg
      logical almost_zero

      integer i0, in, j0, jn
      integer niloc,njloc
      character(len=9) :: communicate_S

      communicate_S = "GRID"
      if (Grd_yinyang_L) communicate_S = "MULTIGRID"


      niloc = (l_ni-pil_e)-(1+pil_w)+1
      njloc = (l_nj-pil_n)-(1+pil_s)+1


      ! Here we go !

      i0 = 1  + sol_pil_w
      in = l_ni - sol_pil_e
      j0 = 1  + sol_pil_s
      jn = l_nj - sol_pil_n

      outiter = 0
      nbiter = 0

      conv = 1.d0

      ! Residual of the initial iterate
      call matvec(solution, work_space, level)

      do j=j0,jn
         do i=i0,in
            vv(i,j,1) = rhs_b(i,j) - work_space(i,j)
         end do
      end do

      do

         local_dot = 0.0d0
         do j=j0,jn
!DIR$ SIMD
            do i=i0,in
               local_dot = local_dot + (vv(i, j, 1) * vv(i, j, 1))
            end do
         end do

         call RPN_COMM_allreduce(local_dot, norm_residual, 1, "MPI_double_precision", "MPI_sum", communicate_S, ierr)
         norm_residual = sqrt(norm_residual)

         ! Current guess is a good enough solution
         if (norm_residual < tolerance) then
            return
         end if

         nu = 1.0d0 / norm_residual
         do j=j0,jn
            do i=i0,in
               vv(i,j,1) = vv(i,j,1) * nu
            end do
         end do

         if (outiter == 0) then
            relative_tolerance = tolerance * norm_residual
            r0 = norm_residual
         end if

         conv = norm_residual / r0

         ! initialize 1-st term of rhs of hessenberg system.
         gg(1) = norm_residual
         gg(2:) = 0.d0

         initer = 0

         do

            nbiter = nbiter + 1
            initer = initer + 1
            nextit = initer + 1

            if (sol2D_precond_S == 'JACOBI')   then
               call pre_jacobi2D ( work_space(i0:in,j0:jn), &
                                   vv(i0:in,j0:jn, initer), &
                                   Prec_xevec_8, niloc, njloc,&
                                   Prec_ai_8, Prec_bi_8, Prec_ci_8 )
            else
               work_space(i0:in,j0:jn) = vv(i0:in,j0:jn, initer)
            endif

            ww(i0:in,j0:jn, initer) = work_space(i0:in,j0:jn)

            call matvec ( work_space, vv(:,:,nextit), level )

            ! Classical Gram-Schmidt orthogonalisation process
            dotprod_local = 0.d0
            do it=1,initer
                local_dot = 0.0d0
                do j=j0,jn
!DIR$ SIMD
                   do i=i0,in
                      local_dot = local_dot + (vv(i, j, it) * vv(i, j, nextit))
                   end do
                end do
                dotprod_local(it) = local_dot
            end do

            call RPN_COMM_allreduce(dotprod_local(:), hessenberg(1,initer), initer, "MPI_double_precision", "MPI_sum", communicate_S, ierr)

            do it=1,initer
               do j=j0,jn
                  do i=i0,in
                     vv(i, j, nextit) = vv(i, j, nextit) - hessenberg(it,initer) * vv(i, j, it)
                  end do
               end do
            end do

            local_dot = 0.d0
            do j=j0,jn
!DIR$ SIMD
               do i=i0,in
                  local_dot = local_dot + (vv(i, j, nextit) * vv(i, j, nextit))
               end do
            end do

            call RPN_COMM_allreduce(local_dot,dotprod,1,"MPI_double_precision","MPI_sum",communicate_S,ierr)

            hessenberg(nextit,initer) = sqrt(dotprod)

            ! Watch out for happy breakdown
            if (.not. almost_zero( hessenberg(nextit,initer) ) ) then
               nu = 1.d0 / hessenberg(nextit,initer)
               do j=j0,jn
                  do i=i0,in
                     vv(i, j, nextit) = vv(i, j, nextit) * nu
                  end do
               end do
            end if

            ! Form and store the information for the new Givens rotation
            if (initer > 1) then
               do k=2,initer
                  k1 = k-1
                  t = hessenberg(k1,initer)
                  hessenberg(k1,initer) = rot_cos(k1)*t + rot_sin(k1)*hessenberg(k,initer)
                  hessenberg(k,initer) = -rot_sin(k1)*t + rot_cos(k1)*hessenberg(k,initer)
               end do

            end if

            nu = sqrt(hessenberg(initer,initer)**2 + hessenberg(nextit,initer)**2)
            if (.not. almost_zero(nu)) then
               rot_cos(initer) = hessenberg(initer,initer) / nu
               rot_sin(initer) = hessenberg(nextit,initer) / nu

               gg(nextit) = -rot_sin(initer) * gg(initer)
               gg(initer) =  rot_cos(initer) * gg(initer)

               hessenberg(initer,initer) = rot_cos(initer) * hessenberg(initer,initer) + rot_sin(initer) * hessenberg(nextit,initer)
            end if

            norm_residual = abs(gg(nextit))

            conv = norm_residual / r0

            if ((initer >= maxinner) .or. (norm_residual <= relative_tolerance)) then
               exit
            endif

         end do

         ! At this point either the maximum number of inner iterations
         ! was reached or the absolute residual is below the scaled tolerance.

         ! Solve upper triangular system
         gg(initer) = gg(initer) / hessenberg(initer,initer)
         do ii=2,initer
            k  = initer - ii + 1
            k1 = k + 1
            t  = gg(k)
            do j=k1,initer
               t = t - hessenberg(k,j) * gg(j)
            end do
            gg(k) = t / hessenberg(k,k)
         end do

         ! Form linear combination to get solution.
         do it=1,initer
            t = gg(it)

            do j=j0,jn
               do i=i0,in
                  solution(i, j) = solution(i, j) + t * ww(i, j, it)
               end do
            end do

         end do

         outiter = outiter + 1

         if (norm_residual <= relative_tolerance .or. outiter > maxouter) then
            return
         end if

         ! Solution is not convergent : compute residual vector and continue.
         do it=1,initer
            jj = nextit - it + 1
            gg(jj-1) = -rot_sin(jj-1) * gg(jj)
            gg(jj)   =  rot_cos(jj-1) * gg(jj)
         end do

         do it=1,nextit
            t = gg(it)
            if (it == 1) then
               t = t - 1.d0
            end if

            do j=j0,jn
!DIR$ SIMD
               do i=i0,in
                  vv(i, j, 1) = vv(i, j, 1) + t * vv(i, j, it)
               end do
            end do

         end do

      end do

      return
      end
