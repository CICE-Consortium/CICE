
!**s/r pgmres -  preconditionner for GEM_H : PGmres
!

       subroutine pgmres(nx_block, ny_block, nblocks, &
                         max_blocks, icellu,          &
                         indxui,     indxuj,          &
                         icellt,                      &
                         indxti,     indxtj,          &
                         dxt,        dyt,             & 
                         dxhy,       dyhx,            & 
                         cxp,        cyp,             & 
                         cxm,        cym,             & 
                         tarear,     tinyarea,        & 
                         vrel,       Cb,              &
                         zetaD,      aiu,             & 
                         umassdti,   fm,              &
                         uarear,     diagvec,         &
                         sol,        rhs,             &
                         n,          im,              &
                         eps,        maxits,          &
                         iout,       ierr)
                         
!-----------------------------------------------------------------------

!        use grid_options
!        use prec
       use ice_kinds_mod
       use ice_dyn_vp, only: matvec, arrays_to_vec, vec_to_arrays, precond_diag
       use ice_fileunits, only: nu_diag

       implicit none

!#include <arch_specific.hf>

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         nblocks,            & ! nb of blocks
         max_blocks            ! max nb of blocks

         
      integer (kind=int_kind), dimension (max_blocks), intent(in) :: &
         icellu  , &
         icellt                ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block, max_blocks), &
         intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj  , & ! compressed index in j-direction
         indxti  , & ! compressed index in i-direction
         indxtj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(in) :: &
         dxt      , & ! width of T-cell through the middle (m)
         dyt      , & ! height of T-cell through the middle (m)
         dxhy     , & ! 0.5*(HTE - HTE)
         dyhx     , & ! 0.5*(HTN - HTN)
         cyp      , & ! 1.5*HTE - 0.5*HTE
         cxp      , & ! 1.5*HTN - 0.5*HTN
         cym      , & ! 0.5*HTE - 1.5*HTE
         cxm      , & ! 0.5*HTN - 1.5*HTN
         tarear   , & ! 1/tarea
         tinyarea     ! puny*tarea   
         
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), &
         intent(in) :: &
         vrel    , & ! coefficient for tauw
         Cb      , & ! coefficient for basal stress
         aiu     , & ! ice fraction on u-grid
         umassdti, & ! mass of U-cell/dt (kg/m^2 s)
         fm      , & ! Coriolis param. * mass in U-cell (kg/s)
         uarear      ! 1/uarea

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks,4), & 
         intent(in) :: &
         zetaD          ! 2*zeta   
         
      real (kind=dbl_kind), dimension (n), intent(in) :: &
         diagvec
      
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         utp     , & ! x-component of velocity (m/s)
         vtp     , & ! y-component of velocity (m/s)         
         Au      , & ! matvec, Fx = Au - bx (N/m^2)! jfl
         Av          ! matvec, Fy = Av - by (N/m^2)! jfl       

       integer n, im, maxits, iout, ierr, iblk
       real*8 rhs(n), sol(n) ,eps ! wk11, wk22, eps
!       Abdessamad Qaddouri -  2018
!
!revision
! v5.0 - Qaddouri A.       - initial version

       real*8 vv(n,im+1), gam,eps1
       real*8  wk(n),r0
!----------------------------------------------------------------------*
       integer kmax,ii,i,j,n1,its,k1,i1,jj,k
       parameter (kmax=50)

       real*8 hh(kmax+1,kmax), c(kmax), s(kmax), rs(kmax+1),t
       real*8 hhloc(kmax+1,kmax)
!-------------------------------------------------------------
! arnoldi size should not exceed kmax=50 in this version..
! to reset modify paramter kmax accordingly.
!-------------------------------------------------------------
       real*8 epsmac ,ro,ddot,dnrm2
       parameter (epsmac=1.d-16)
       integer l
!       character(len= 9) communicate_S
!       communicate_S = "GRID"
!       if (Grd_yinyang_L) communicate_S = "MULTIGRID"


       n1 = n + 1
       its = 0
       sol=0.0 !JFL ...veut-on vraiment mettre sol = 0 ici??????
!-------------------------------------------------------------
! outer loop starts here..
!-------------- compute initial residual vector --------------
       do 21 j=1,n
         vv(j,1) = rhs(j) 
 21    continue
       
!-------------------------------------------------------------
  20   continue
       ro = ddot(n, vv,1,vv,1)
       ro = dsqrt(ro)

        r0=ro
 
       if (ro .eq. 0.0d0) goto 999
       t = 1.0d0/ ro
       do 210 j=1, n
          vv(j,1) = vv(j,1)*t
 210   continue
       if (its .eq. 0) eps1=eps*ro
       if (iout .gt. 0 .and. its .eq. 0)&
          write(nu_diag, 199) its, ro ,eps1
!     ** initialize 1-st term  of rhs of hessenberg system..
       rs(1) = ro
       i = 0
 4     i=i+1
       its = its + 1
       i1 = i + 1

       do l=1,n
       rhs(l)= 0.0 
       wk(l)=  vv(l,i)
       enddo
! precond
       call precond_diag (n,            & 
                          diagvec (:),     &
                          wk (:), rhs (:) )

!       rhs = wk !!! JFL       

!  matrix-vector
!        call sol_matvec_H JFL

       call vec_to_arrays (nx_block, ny_block, nblocks,      &
                           max_blocks, icellu (:), n,        & 
                           indxui    (:,:), indxuj(:,:),     &
                           rhs (:),                         &
                           utp (:,:,:), vtp (:,:,:))    
                         

       !$OMP PARALLEL DO PRIVATE(iblk)
       do iblk = 1, nblocks                          
        call matvec (nx_block             , ny_block,            &
                     icellu   (iblk)      , icellt   (iblk)    , & 
                     indxui   (:,iblk)    , indxuj   (:,iblk)  , &
                     indxti   (:,iblk)    , indxtj   (:,iblk)  , &
                     dxt      (:,:,iblk)  , dyt      (:,:,iblk), & 
                     dxhy     (:,:,iblk)  , dyhx     (:,:,iblk), & 
                     cxp      (:,:,iblk)  , cyp      (:,:,iblk), & 
                     cxm      (:,:,iblk)  , cym      (:,:,iblk), & 
                     tarear   (:,:,iblk)  , tinyarea (:,:,iblk), &
                     utp      (:,:,iblk)  , vtp      (:,:,iblk), &      
                     vrel     (:,:,iblk)  , Cb       (:,:,iblk), &  
                     zetaD    (:,:,iblk,:), aiu      (:,:,iblk), &
                     umassdti (:,:,iblk)  , fm       (:,:,iblk), & 
                     uarear   (:,:,iblk)  ,                      & 
                     Au       (:,:,iblk)  , Av       (:,:,iblk)) 
                    
       enddo
       !$OMP END PARALLEL DO                     
                    
               ! form wk2 from Au and Av arrays        
       call arrays_to_vec (nx_block, ny_block, nblocks,      &
                           max_blocks, icellu (:), n,        & 
                           indxui    (:,:), indxuj(:,:),     &
                           Au      (:,:,:), Av    (:,:,:),   &
                           vv(1,i1))                  
                         
!     classical gram - schmidt...
!
      do 55 j=1, i
         hhloc(j,i) = ddot(n, vv(1,j), 1, vv(1,i1), 1)
         hh(j,i) =  hhloc(j,i)
 55   continue

      do 56 j=1, i
         call daxpy(n, -hh(j,i), vv(1,j), 1, vv(1,i1), 1)
 56   continue
      t = ddot(n, vv(1,i1), 1, vv(1,i1), 1)
!
       t=dsqrt(t)
!

       hh(i1,i) = t
       if ( t .eq. 0.0d0) goto 58
       t = 1.0d0/t
       do 57  k=1,n
          vv(k,i1) = vv(k,i1)*t
 57    continue
!
!     done with modified gram schimd and arnoldi step..
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

!     if gamma is zero then any small value will do...
!     will affect only residual estimate
!
      if (gam == 0.0d0) gam = epsmac
!-----------#determinenextplane rotation  #-------------------
      c(i) = hh(i,i)/gam
      s(i) = hh(i1,i)/gam
      rs(i1) = -s(i)*rs(i)
      rs(i) =  c(i)*rs(i)

!
!     detrermine residual norm and test for convergence-
!
       hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
       ro = abs(rs(i1))
       if (iout .gt. 0) &
           write(nu_diag, 199) its, ro , eps1
       if (i .lt. im .and. (ro .gt. eps1))  goto 4
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
 40       continue
          rs(k) = t/hh(k,k)
 30    continue
!
!     form linear combination of
!,i)'s to get solution
!
       t = rs(1)
       do 15 k=1, n
          rhs(k) = vv(k,1)*t
 15    continue
       do 16 j=2, i
          t = rs(j)
          do 161 k=1, n
             rhs(k) = rhs(k)+t*vv(k,j)
 161      continue
 16    continue
!
!     call preconditioner.
!

       do l=1,n
       wk(l)=  rhs(l)
       rhs(l)=0.0
       enddo
! precond
       call precond_diag (n,            & 
                          diagvec (:),     &
                          wk (:), rhs (:) )
!       rhs = wk !!! JFL       

       do 17 k=1, n
          sol(k) = sol(k) + rhs(k)
 17    continue
!
!     restart outer loop  when necessary
!
       if (ro .le. eps1) goto 990
       if (its .ge. maxits) goto 991
!
!     else compute residual vector and continue..
!
       do 24 j=1,i
          jj = i1-j+1
          rs(jj-1) = -s(jj-1)*rs(jj)
          rs(jj) = c(jj-1)*rs(jj)
 24    continue
       do 25  j=1,i1
          t = rs(j)
          if (j .eq. 1)  t = t-1.0d0
          call daxpy (n, t, vv(1,j), 1,  vv, 1)
 25    continue
 199   format('monitor_pgmres: iter_pmgres=', i4, ' L2norm=', d26.16, ' epsprecond*initial_L2norm=', d26.6)
!     restart outer loop.
       goto 20
 990   ierr = 0
!           write(iout, 198) its, ro/r0 
 198   format('   its =', i4, ' conv =', d20.6)
       return
 991   ierr = 1
!           write(iout, 198) its, ro/r0

       return
 999   continue
       ierr = -1
!           write(iout, 198) its, ro/r0

       return
!-----------------end of pgmres ---------------------------------------
!-----------------------------------------------------------------------
       end
