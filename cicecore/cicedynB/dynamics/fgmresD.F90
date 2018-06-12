      subroutine fgmres (n,im,rhs,sol,i,vv,w,wk1, wk2, &
                  eps,maxits,iout,icode,its) 

!-----------------------------------------------------------------------
! jfl Dec 1st 2006. We modified the routine so that it is double precison.
! Here are the modifications:
! 1) implicit real (a-h,o-z) becomes implicit real*8 (a-h,o-z) 
! 2) real bocomes real*8
! 3) subroutine scopy.f has been changed for dcopy.f
! 4) subroutine saxpy.f has been changed for daxpy.f
! 5) function sdot.f has been changed for ddot.f
! 6) 1e-08 becomes 1d-08
!
! Be careful with the dcopy, daxpy and ddot code...there is a slight 
! difference with the single precision versions (scopy, saxpy and sdot).
! In the single precision versions, the array are declared sightly differently.
! It is written for single precision:
!
! modified 12/3/93, array(1) declarations changed to array(*)
!-----------------------------------------------------------------------

      implicit double precision (a-h,o-z) !jfl modification
      integer n, im, maxits, iout, icode
      double precision rhs(*), sol(*), vv(n,im+1),w(n,im)
      double precision wk1(n), wk2(n), eps
!-----------------------------------------------------------------------
! flexible GMRES routine. This is a version of GMRES which allows a 
! a variable preconditioner. Implemented with a reverse communication 
! protocole for flexibility -
! DISTRIBUTED VERSION (USES DISTDOT FOR DDOT) 
! explicit (exact) residual norms for restarts  
! written by Y. Saad, modified by A. Malevsky, version February 1, 1995
!-----------------------------------------------------------------------
! This Is A Reverse Communication Implementation. 
!------------------------------------------------- 
! USAGE: (see also comments for icode below). FGMRES
! should be put in a loop and the loop should be active for as
! long as icode is not equal to 0. On return fgmres will
!    1) either be requesting the new preconditioned vector applied
!       to wk1 in case icode.eq.1 (result should be put in wk2) 
!    2) or be requesting the product of A applied to the vector wk1
!       in case icode.eq.2 (result should be put in wk2) 
!    3) or be terminated in case icode .eq. 0. 
! on entry always set icode = 0. So icode should be set back to zero
! upon convergence.
!-----------------------------------------------------------------------
! Here is a typical way of running fgmres: 
!
!      icode = 0
! 1    continue
!      call fgmres (n,im,rhs,sol,i,vv,w,wk1, wk2,eps,maxits,iout,icode)
!
!      if (icode .eq. 1) then
!         call  precon(n, wk1, wk2)    <--- user's variable preconditioning
!         goto 1
!      else if (icode .ge. 2) then
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
! im    == size of Krylov subspace:  should not exceed 50 in this
!          version (can be reset in code. looking at comment below)
! rhs   == vector of length n containing the right hand side
! sol   == initial guess on input, approximate solution on output
! vv    == work space of size n x (im+1)
! w     == work space of length n x im 
! wk1,
! wk2,  == two work vectors of length n each used for the reverse
!          communication protocole. When on return (icode .ne. 1)
!          the user should call fgmres again with wk2 = precon * wk1
!          and icode untouched. When icode.eq.1 then it means that
!          convergence has taken place.
!          
! eps   == tolerance for stopping criterion. process is stopped
!          as soon as ( ||.|| is the euclidean norm):
!          || current residual||/||initial residual|| <= eps
!
! maxits== maximum number of iterations allowed
!
! iout  == output unit number number for printing intermediate results
!          if (iout .le. 0) no statistics are printed.
! 
! icode = integer. indicator for the reverse communication protocole.
!         ON ENTRY : icode should be set to icode = 0.
!         ON RETURN: 
!       * icode .eq. 1 value means that fgmres has not finished
!         and that it is requesting a preconditioned vector before
!         continuing. The user must compute M**(-1) wk1, where M is
!         the preconditioing  matrix (may vary at each call) and wk1 is
!         the vector as provided by fgmres upun return, and put the 
!         result in wk2. Then fgmres must be called again without
!         changing any other argument. 
!       * icode .eq. 2 value means that fgmres has not finished
!         and that it is requesting a matrix vector product before
!         continuing. The user must compute  A * wk1, where A is the
!         coefficient  matrix and wk1 is the vector provided by 
!         upon return. The result of the operation is to be put in
!         the vector wk2. Then fgmres must be called again without
!         changing any other argument. 
!       * icode .eq. 0 means that fgmres has finished and sol contains 
!         the approximate solution.
!         comment: typically fgmres must be implemented in a loop
!         with fgmres being called as long icode is returned with 
!         a value .ne. 0. 
!-----------------------------------------------------------------------
!     local variables -- !jfl modif
      double precision hh(201,200),c(200),s(200),rs(201),t,ro,ddot,sqrt 
!
!-------------------------------------------------------------
!     arnoldi size should not exceed 50 in this version..
!     to reset modify sizes of hh, c, s, rs       
!-------------------------------------------------------------

      save
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
! 10   continue
      call dcopy (n, sol, 1, wk1, 1) !jfl modification
      icode = 3
      return
 11   continue
      do j=1,n
         vv(j,1) = rhs(j) - wk2(j) 
      enddo
 20   ro = ddot(n, vv, 1, vv,1) !jfl modification
      ro = sqrt(ro)
      if (ro .eq. 0.0d0) goto 999 
      t = 1.0d0/ ro 
      do j=1, n
         vv(j,1) = vv(j,1)*t 
      enddo
      if (its .eq. 0) eps1=eps
      if (its .eq. 0) r0 = ro
      if (iout .gt. 0) write(*, 199) its, ro!&
!           print *,'chau',its, ro !write(iout, 199) its, ro
!     
!     initialize 1-st term  of rhs of hessenberg system..
!     
      rs(1) = ro
      i = 0
 4    i=i+1
      its = its + 1
      i1 = i + 1
      do k=1, n
         wk1(k) = vv(k,i) 
      enddo
!     
!     return
!     
      icode = 1

      return
 200  continue
      do k=1, n
         w(k,i) = wk2(k) 
      enddo
!     
!     call matvec operation
!     
      icode = 2
      call dcopy(n, wk2, 1, wk1, 1) !jfl modification
!
!     return
!     
      return
 300  continue
!     
!     first call to ope corresponds to intialization goto back to 11.
!     
!      if (icode .eq. 3) goto 11
      call  dcopy (n, wk2, 1, vv(1,i1), 1) !jfl modification
!     
!     modified gram - schmidt...
!     
      do j=1, i
         t = ddot(n, vv(1,j), 1, vv(1,i1), 1) !jfl modification
         hh(j,i) = t
         call daxpy(n, -t, vv(1,j), 1, vv(1,i1), 1) !jfl modification
      enddo
      t = sqrt(ddot(n, vv(1,i1), 1, vv(1,i1), 1)) !jfl modification
      hh(i1,i) = t
      if (t .eq. 0.0d0) goto 58
      t = 1.0d0 / t
      do k=1,n
         vv(k,i1) = vv(k,i1)*t
      enddo
!     
!     done with modified gram schimd and arnoldi step. 
!     now  update factorization of hh
!     
 58   if (i .eq. 1) goto 121
!     
!     perfrom previous transformations  on i-th column of h
!     
      do k=2,i
         k1 = k-1
         t = hh(k1,i)
         hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
         hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
      enddo
 121  gam = sqrt(hh(i,i)**2 + hh(i1,i)**2)
      if (gam .eq. 0.0d0) gam = epsmac
!-----------#determine next plane rotation  #-------------------
      c(i) = hh(i,i)/gam
      s(i) = hh(i1,i)/gam
      rs(i1) = -s(i)*rs(i)
      rs(i) =  c(i)*rs(i)
!     
!     determine res. norm. and test for convergence-
!     
      hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
      ro = abs(rs(i1))
      if (iout .gt. 0) &
           write(*, 199) its, ro
      if (i .lt. im .and. (ro .gt. eps1))  goto 4
!     
!     now compute solution. first solve upper triangular system.
!     
      rs(i) = rs(i)/hh(i,i)
      do ii=2,i
         k=i-ii+1
         k1 = k+1
         t=rs(k)
         do j=k1,i
            t = t-hh(k,j)*rs(j)
         enddo
         rs(k) = t/hh(k,k)
      enddo
!     
!     done with back substitution..
!     now form linear combination to get solution
!     
      do j=1, i
         t = rs(j)
         call daxpy(n, t, w(1,j), 1, sol,1) !jfl modification
      enddo
!     
!     test for return 
!     
      if (ro .le. eps1 .or. its .ge. maxits) goto 999
!     
!     else compute residual vector and continue..
!     
!       goto 10

     do j=1,i
        jj = i1-j+1
        rs(jj-1) = -s(jj-1)*rs(jj)
        rs(jj) = c(jj-1)*rs(jj)
     enddo
     do j=1,i1
        t = rs(j)
        if (j .eq. 1)  t = t-1.0d0
        call daxpy (n, t, vv(1,j), 1,  vv, 1)
     enddo
!     
!     restart outer loop.
!     
     goto 20
 999  icode = 0

 199  format('   -- fmgres its =', i4, ' res. norm =', d26.16)
!     
      return 
!-----end-of-fgmres----------------------------------------------------- 
!-----------------------------------------------------------------------
      end

