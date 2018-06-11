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

!** fgmres - flexible GMRES routine to allow a variable preconditioner
!
      subroutine fgmres2(n,im,rhs,sol,i,vv,w,wk1, wk2,eps,maxits, &
                                                  its,conv,icode)
      implicit none
!#include <arch_specific.hf>

      include 'mpif.h'   ! MPI Fortran include file

      integer n, im, i, maxits, its, icode
      real*8 rhs(*), sol(*), vv(n,im+1),w(n,im), wk1(n), wk2(n), eps,conv
!
!author Y. Saad (modified by A. Malevsky Feb1, 1995)
!
!revision
! v3_30 - Abdessamad Qaddouri ; add RPN_comm calls
!-----------------------------------------------------------------------
! flexible GMRES routine. This is a version of GMRES which allows a
! a variable preconditioner. Implemented with a reverse communication
! protocole for flexibility -
! DISTRIBUTED VERSION (USES DISTDOT FOR DDOT)
! explicit (exact) residual norms for restarts
! written by Y. Saad, modified by A. Malevsky, version February 1, 1995
!  revision :
!            Abdessamad Qaddouri ; add RPN_comm calls
!-----------------------------------------------------------------------
! This Is A Reverse Communication Implementation.
!-------------------------------------------------
! USAGE: (see also comments for icode below). CGMRES
! should be put in a loop and the loop should be active for as
! long as icode is not equal to 0. On return fgmres will
!    1) either be requesting the new preconditioned vector applied
!       to wk1 in case icode == 1 (result should be put in wk2)
!    2) or be requesting the product of A applied to the vector wk1
!       in case icode == 2 (result should be put in wk2)
!    3) or be terminated in case icode == 0.
! on entry always set icode = 0. So icode should be set back to zero
! upon convergence.
!-----------------------------------------------------------------------
! Here is a typical way of running fgmres:
!
!      icode = 0
! 1    continue
!      call fgmres2(n,im,rhs,sol,i,vv,w,wk1, wk2,eps,maxits,its,conv,icode)
!
!      if (icode == 1) then
!         call  precon(n, wk1, wk2)    <--- user's variable preconditioning
!         goto 1
!      else if (icode >= 2) then
!         call  matvec (n,wk1, wk2)    <--- user's matrix vector product.
!         goto 1
!      else
!         ----- done ----
!         .........
!-----------------------------------------------------------------------
! list of parameters
!-------------------
!
! n     == integer. the dimension of the problem
! im    == size of Krylov subspace:  should not exceed 100 in this
!          version (can be reset in code. looking at comment below)
! rhs   == vector of length n containing the right hand side
! sol   == initial guess on input, approximate solution on output
! vv    == work space of size n x (im+1)
! w     == work space of length n x im
! wk1,
! wk2,  == two work vectors of length n each used for the reverse
!          communication protocole. When on return (icode /= 1)
!          the user should call fgmres again with wk2 = precon * wk1
!          and icode untouched. When icode == 1 then it means that
!          convergence has taken place.
!
! eps   == tolerance for stopping criterion. process is stopped
!          as soon as ( ||.|| is the euclidean norm):
!          || current residual||/||initial residual|| <= eps
!
! maxits== maximum number of iterations allowed
!
!
! icode = integer. indicator for the reverse communication protocole.
!         ON ENTRY : icode should be set to icode = 0.
!         ON RETURN:
!       * icode == 1 value means that fgmres has not finished
!         and that it is requesting a preconditioned vector before
!         continuing. The user must compute M**(-1) wk1, where M is
!         the preconditioing  matrix (may vary at each call) and wk1 is
!         the vector as provided by fgmres upun return, and put the
!         result in wk2. Then fgmres must be called again without
!         changing any other argument.
!       * icode == 2 value means that fgmres has not finished
!         and that it is requesting a matrix vector product before
!         continuing. The user must compute  A * wk1, where A is the
!         coefficient  matrix and wk1 is the vector provided by
!         upon return. The result of the operation is to be put in
!         the vector wk2. Then fgmres must be called again without
!         changing any other argument.
!       * icode == 0 means that fgmres has finished and sol contains
!         the approximate solution.
!         comment: typically fgmres must be implemented in a loop
!         with fgmres being called as long icode is returned with
!         a value /= 0.
!-----------------------------------------------------------------------
!     local variables --
      real*8 hh(101,100),hhloc(101,100),c(100),s(100), &
           rs(101),t,tloc,ro,dsqrt,epsmac, ddot, eps1, gam ,r0
      integer n1, j, i1, k, k1, ii, jj, ierr
!      integer mproc, myproc ! use when printing out ..
!-------------------------------------------------------------
!     arnoldi size should not exceed 100 in this version..
!     to reset modify sizes of hh, c, s, rs
!-------------------------------------------------------------
!
      save
! ##
!     used for printing out only -- ignore
!      call MPI_Comm_rank(MPI_COMM_WORLD,mproc,ierr)
!      myproc = mproc+1
!
      data epsmac/1.d-16/
!
!     computed goto
!
      goto (100,200,300,11) icode +1
 100  continue
      n1 = n + 1
      its = 0
!-------------------------------------------------------------
!     **  outer loop starts here..
!--------------compute initial residual vector --------------
 10   continue
      call dcopy (n, sol, 1, wk1, 1)
      icode = 3
      return
 11   continue
      do 21 j=1,n
         vv(j,1) = rhs(j) - wk2(j)
 21   continue
 20   continue
       tloc=ddot(n, vv, 1, vv,1)
      call MPI_allreduce(tloc,ro,1,MPI_double_precision, &
           MPI_sum,MPI_COMM_WORLD,ierr)
!       call RPN_COMM_allreduce(tloc,ro,1,"MPI_double_precision", &
!           "MPI_sum","grid",ierr)
       ro = dsqrt(ro)
!##
!      if (mproc == 0) write (19,*) ro
      if (ro == 0.0d0) goto 999
      t = 1.0d0/ ro
      call dscal(n,t,vv(1,1),1)
      if (its == 0) eps1=eps*ro
      if (its == 0) r0 = ro

      conv = ro/r0
!
!     initialize 1-st term  of rhs of hessenberg system..
!
      rs(1) = ro
      i = 0
 4    i=i+1
      its = its + 1
      i1 = i + 1
      call dcopy(n,vv(1,i),1,wk1,1)
!
!     return
!
      icode = 1
      return
 200  continue
      call dcopy(n, wk2, 1, w(1,i), 1)
!
!     call matvec operation
!
      icode = 2
      call dcopy(n, wk2, 1, wk1, 1)
!
!     return
!
      return
 300  continue
!
!     first call to ope corresponds to intialization goto back to 11.
!
!      if (icode == 3) goto 11
      call  dcopy (n, wk2, 1, vv(1,i1), 1)
!
!     classical gram - schmidt...
!
      do 55 j=1, i
         hhloc(j,i) = ddot(n, vv(1,j), 1, vv(1,i1), 1)
 55   continue
      call MPI_allreduce(hhloc(1,i),hh(1,i),i,MPI_double_precision, &
                         MPI_sum,MPI_COMM_WORLD,ierr)
!       call RPN_COMM_allreduce(hhloc(1,i),hh(1,i),i,"MPI_double_precision", &
!           "MPI_sum","grid",ierr)

      do 56 j=1, i
         call daxpy(n, -hh(j,i), vv(1,j), 1, vv(1,i1), 1)
 56   continue
      tloc = ddot(n, vv(1,i1), 1, vv(1,i1), 1)
!
      call MPI_allreduce(tloc,t,1,MPI_double_precision, &
                         MPI_sum,MPI_COMM_WORLD,ierr)
!       call RPN_COMM_allreduce(tloc,t,1,"MPI_double_precision", &
!           "MPI_sum","grid",ierr)

      t = sqrt(t)
      hh(i1,i) = t
      if (t == 0.0d0) goto 58
      t = 1.0d0 / t
      call dscal(n,t,vv(1,i1),1)
!
!     done with classical gram schimd and arnoldi step.
!     now  update factorization of hh
!
 58   if (i == 1) goto 121
!
!     perfrom previous transformations  on i-th column of h
!
      do 66 k=2,i
         k1 = k-1
         t = hh(k1,i)
         hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
         hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
 66   continue
 121  gam = sqrt(hh(i,i)**2 + hh(i1,i)**2)
      if (gam == 0.0d0) gam = epsmac
!-----------#determinenextplane rotation  #-------------------
      c(i) = hh(i,i)/gam
      s(i) = hh(i1,i)/gam
      rs(i1) = -s(i)*rs(i)
      rs(i) =  c(i)*rs(i)
!
!     determine res. norm. and test for convergence-
!
      hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
      ro = abs(rs(i1))
!##
!      if (mproc == 0) write (19,*) ro
!
      conv = ro/r0

      if ((i < im) .and. (ro > eps1)) then

      goto 4
      endif
!
!     now compute solution. first solve upper triangular system.
!
      rs(i) = rs(i)/hh(i,i)
      do 30 ii=2,i
         k=i-ii+1
         k1 = k+1
         t=rs(k)
         do 40 j=k1,i
            t = t-hh(k,j)*rs(j)
 40      continue
         rs(k) = t/hh(k,k)
 30   continue
!
!     done with back substitution..
!     now form linear combination to get solution
!
      do 16 j=1, i
         t = rs(j)
         call daxpy(n, t, w(1,j), 1, sol,1)
 16   continue
!
!     test for return
      if (ro <= eps1 .or. its >= maxits) goto 999
!
!     else compute residual vector and continue..
!
!      goto 10
      do 24 j=1,i
         jj = i1-j+1
         rs(jj-1) = -s(jj-1)*rs(jj)
         rs(jj) = c(jj-1)*rs(jj)
 24   continue
      do 25  j=1,i1
         t = rs(j)
         if (j == 1)  t = t-1.0d0
         call daxpy (n, t, vv(1,j), 1,  vv, 1)
 25   continue
!
!     restart outer loop.
!
      goto 20
 999  icode = 0
!
      return
!-----end-of-fgmres-----------------------------------------------------
!-----------------------------------------------------------------------
      end
