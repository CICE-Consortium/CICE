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
!=======================================================================
!
! Krylov subspace methods for sea-ice dynamics
!
! See:
!
! C. T. Kelley, Iterative Methods for Linear and Nonlinear Equations, SIAM, 1995
! (https://www.siam.org/books/textbooks/fr16_book.pdf)
!
! Y. Saad, Iterative Methods for Sparse Linear Systems. SIAM, 2003.
! (http://www-users.cs.umn.edu/~saad/IterMethBook_2ndEd.pdf)
!
! author: Philippe Blain, ECCC

      module ice_krylov

      use ice_blocks, only: nx_block, ny_block
      use ice_constants, only: field_loc_NEcorner, c0, c1
      use ice_domain, only: distrb_info, nblocks
      use ice_domain_size, only: max_blocks
      use ice_dyn_vp, only: matvec, residual_vec, calc_L2norm_squared
      use ice_flux, only:  fm, iceumask
      use ice_global_reductions, only: global_sum
      use ice_grid, only: dxt, dyt, dxhy, dyhx, cxp, cyp, cxm, cym, &
          uarear, tinyarea
      use ice_kinds_mod

      implicit none
      private
      public :: fgmres, pgmres

!=======================================================================

      contains

!=======================================================================

! FGMRES: Flexible generalized minimum residual method (with restarts). 
! Solves A x = b using GMRES with a varying (right) preconditioner
!
! authors: Stéphane Gaudreault, Abdessamad Qaddouri, Philippe Blain, ECCC

      subroutine fgmres (icellt,   icellu,  &
                         indxti,   indxtj,  &
                         indxui,   indxuj,  &
                         zetaD,             &
                         Cb,       vrel,    &
                         umassdti,          &
                         solx,     soly,    &
                         bx,       by,      &
                         tolerance, maxinner, maxouter, nbiter, conv)

      integer (kind=int_kind), dimension(max_blocks), intent(in) :: & 
         icellt   , & ! no. of cells where icetmask = 1
         icellu       ! no. of cells where iceumask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block, max_blocks), intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj   , & ! compressed index in j-direction
         indxui   , & ! compressed index in i-direction
         indxuj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks,4), intent(in) :: &
         zetaD      ! zetaD = 2*zeta (viscous coefficient)

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(in) :: &
         vrel     , & ! coefficient for tauw 
         Cb       , & ! seabed stress coefficient
         umassdti     ! mass of U-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(inout) :: &
         solx     , & ! Initial guess on input, approximate solution on output (x components)
         soly         ! Initial guess on input, approximate solution on output (y components)

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(in) :: &
         bx       , & ! Right hand side of the linear system (x components)
         by           ! Right hand side of the linear system (y components)

      real (kind=dbl_kind), intent(in) :: &
         tolerance   ! Tolerance to achieve. The algorithm terminates when the relative
                     ! residual is below tolerance

      integer (kind=int_kind), intent(in) :: &
         maxinner    ! Restart the method every maxinner inner iterations

      integer (kind=int_kind), intent(in) :: &
         maxouter    ! Maximum number of outer iterations
                     ! Iteration will stop after maxinner*maxouter steps
                     ! even if the specified tolerance has not been achieved

      integer (kind=int_kind), intent(out) :: &
         nbiter      ! Total number of iteration performed

      real (kind=dbl_kind), intent(out) :: &
         conv        ! !phb DESCRIBE IF WE KEEP

      ! local variables

      integer (kind=int_kind) :: &
         iblk    , & ! block index
         ij      , & ! compressed index
         i, j        ! grid indices
         
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         Au       , & ! matvec result, Au = Ax * u (N/m^2)
         Av       , & ! matvec result, Av = Ay * v (N/m^2)
         Fx       , & ! residual vector (x components), Fx = Au - bx (N/m^2)
         Fy           ! residual vector (y components), Fy = Av - by (N/m^2)
! real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks) :: work_space
      ! integer (kind=int_kind) :: i, j, k, k1, ii, jj, ierr

      real (kind=dbl_kind), dimension (max_blocks) :: &
         norm_squared   ! array to accumulate squared norm of grid function over blocks

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks, maxinner+1) :: &
         arnoldi_basis_x , & ! arnoldi basis (x components) !phb == vv
         arnoldi_basis_y     ! arnoldi basis (y components)

      real (kind=dbl_kind) :: &
         norm_residual   , & ! current L^2 norm of residual vector
         inverse_norm    , & ! inverse of the residual norm
         t                   ! 

      integer (kind=int_kind) :: &
         initer          , & ! inner (Arnoldi) loop counter
         outiter         , & ! outer (restarts) loop counter
         nextit          , & ! nextit == initer+1
         it                  ! reusable loop counter

      real (kind=dbl_kind), dimension(maxinner+1) :: &
         rot_cos         , & ! 
         rot_sin         , & ! 
         rhs_hess            ! right hand side vector of the Hessenberg (least square) system

      real (kind=dbl_kind) :: relative_tolerance, r0

      real (kind=dbl_kind), dimension(maxinner+1, maxinner) :: hessenberg

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks, maxinner) :: ww


      real (kind=dbl_kind) :: local_dot, dotprod
      real (kind=dbl_kind), dimension(maxinner) :: dotprod_local

      
      logical (kind=log_kind) almost_zero

      ! integer (kind=int_kind) i0, in, j0, jn
      ! integer (kind=int_kind) niloc,njloc

      ! Here we go !

      outiter = 0
      nbiter = 0

      conv = c1

      ! Residual of the initial iterate
      
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
                      solx     (:,:,iblk)  , soly     (:,:,iblk), &
                      vrel     (:,:,iblk)  , Cb       (:,:,iblk), &
                      zetaD    (:,:,iblk,:),                      &
                      umassdti (:,:,iblk)  , fm       (:,:,iblk), &
                      uarear   (:,:,iblk)  ,                      &
                      Au       (:,:,iblk)  , Av       (:,:,iblk))
         call residual_vec (nx_block           , ny_block,           &
                            icellu       (iblk),                     & 
                            indxui     (:,iblk), indxuj    (:,iblk), &
                            bx       (:,:,iblk), by      (:,:,iblk), &
                            Au       (:,:,iblk), Av      (:,:,iblk), &
                            arnoldi_basis_x (:,:,iblk, 1),           &
                            arnoldi_basis_y (:,:,iblk, 1))
      enddo
      !$OMP END PARALLEL DO
!       do j=j0,jn
!          do i=i0,in
!             vv(i,j,1) = rhs_b(i,j) - work_space(i,j)
!          end do
!       end do
! 
      ! Start outer (restarts) loop
      do
         ! Compute norm of initial residual
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            call calc_L2norm_squared(nx_block,        ny_block,        &
                                     icellu   (iblk),                  &
                                     indxui (:,iblk), indxuj(:, iblk), &
                                     arnoldi_basis_x(:,:,iblk, 1),     &
                                     arnoldi_basis_y(:,:,iblk, 1),     &
                                     norm_squared(iblk))

         enddo
         !$OMP END PARALLEL DO 
         norm_residual = sqrt(global_sum(sum(norm_squared), distrb_info))
! 
!          local_dot = 0.0d0
!          do j=j0,jn
! !DIR$ SIMD
!             do i=i0,in
!                local_dot = local_dot + (vv(i, j, 1) * vv(i, j, 1))
!             end do
!          end do
! 
!          call RPN_COMM_allreduce(local_dot, norm_residual, 1, "MPI_double_precision", "MPI_sum", communicate_S, ierr)
!          norm_residual = sqrt(norm_residual)

      ! Current guess is a good enough solution
      if (norm_residual < tolerance) then
         return
      end if
      
      ! Normalize the first Arnoldi vector
      inverse_norm = c1 / norm_residual
      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
         do ij =1, icellu(iblk)
            i = indxui(ij, iblk)
            j = indxuj(ij, iblk)

            arnoldi_basis_x(i, j, iblk, 1) = arnoldi_basis_x(i, j, iblk, 1) * inverse_norm
            arnoldi_basis_y(i, j, iblk, 1) = arnoldi_basis_y(i, j, iblk, 1) * inverse_norm
         enddo ! ij
      enddo
      !$OMP END PARALLEL DO 
! 
!          nu = 1.0d0 / norm_residual
!          do j=j0,jn
!             do i=i0,in
!                vv(i,j,1) = vv(i,j,1) * nu
!             end do
!          end do
! 
      if (outiter == 0) then
         relative_tolerance = tolerance * norm_residual
         r0 = norm_residual
      end if
      
      conv = norm_residual / r0
      
      ! Initialize 1-st term of RHS of hessenberg system
      gg(1) = norm_residual
      gg(2:) = 0.d0
! 
!          initer = 0
! 
!          do
! 
!             nbiter = nbiter + 1
!             initer = initer + 1
!             nextit = initer + 1
! 
!             ! here call precond_diag OR PGMRES
!             ! if (sol2D_precond_S == 'JACOBI')   then
!             !    call pre_jacobi2D ( work_space(i0:in,j0:jn), &
!             !                        vv(i0:in,j0:jn, initer), &
!             !                        Prec_xevec_8, niloc, njloc,&
!             !                        Prec_ai_8, Prec_bi_8, Prec_ci_8 )
!             ! else
!             !    work_space(i0:in,j0:jn) = vv(i0:in,j0:jn, initer)
!             ! endif
! 
!             ww(i0:in,j0:jn, initer) = work_space(i0:in,j0:jn)
! 
!             call matvec ( work_space, vv(:,:,nextit), level )
! 
!             ! Classical Gram-Schmidt orthogonalisation process
!             dotprod_local = 0.d0
!             do it=1,initer
!                 local_dot = 0.0d0
!                 do j=j0,jn
! !DIR$ SIMD
!                    do i=i0,in
!                       local_dot = local_dot + (vv(i, j, it) * vv(i, j, nextit))
!                    end do
!                 end do
!                 dotprod_local(it) = local_dot
!             end do
! 
!             call RPN_COMM_allreduce(dotprod_local(:), hessenberg(1,initer), initer, "MPI_double_precision", "MPI_sum", communicate_S, ierr)
! 
!             do it=1,initer
!                do j=j0,jn
!                   do i=i0,in
!                      vv(i, j, nextit) = vv(i, j, nextit) - hessenberg(it,initer) * vv(i, j, it)
!                   end do
!                end do
!             end do
! 
!             local_dot = 0.d0
!             do j=j0,jn
! !DIR$ SIMD
!                do i=i0,in
!                   local_dot = local_dot + (vv(i, j, nextit) * vv(i, j, nextit))
!                end do
!             end do
! 
!             call RPN_COMM_allreduce(local_dot,dotprod,1,"MPI_double_precision","MPI_sum",communicate_S,ierr)
! 
!             hessenberg(nextit,initer) = sqrt(dotprod)
! 
!             ! Watch out for happy breakdown
!             if (.not. almost_zero( hessenberg(nextit,initer) ) ) then
!                nu = 1.d0 / hessenberg(nextit,initer)
!                do j=j0,jn
!                   do i=i0,in
!                      vv(i, j, nextit) = vv(i, j, nextit) * nu
!                   end do
!                end do
!             end if
! 
!             ! Form and store the information for the new Givens rotation
!             if (initer > 1) then
!                do k=2,initer
!                   k1 = k-1
!                   t = hessenberg(k1,initer)
!                   hessenberg(k1,initer) = rot_cos(k1)*t + rot_sin(k1)*hessenberg(k,initer)
!                   hessenberg(k,initer) = -rot_sin(k1)*t + rot_cos(k1)*hessenberg(k,initer)
!                end do
! 
!             end if
! 
!             nu = sqrt(hessenberg(initer,initer)**2 + hessenberg(nextit,initer)**2)
!             if (.not. almost_zero(nu)) then
!                rot_cos(initer) = hessenberg(initer,initer) / nu
!                rot_sin(initer) = hessenberg(nextit,initer) / nu
! 
!                gg(nextit) = -rot_sin(initer) * gg(initer)
!                gg(initer) =  rot_cos(initer) * gg(initer)
! 
!                hessenberg(initer,initer) = rot_cos(initer) * hessenberg(initer,initer) + rot_sin(initer) * hessenberg(nextit,initer)
!             end if
! 
!             norm_residual = abs(gg(nextit))
! 
!             conv = norm_residual / r0
! 
!             if ((initer >= maxinner) .or. (norm_residual <= relative_tolerance)) then
!                exit
!             endif
! 
!          end do
! 
!          ! At this point either the maximum number of inner iterations
!          ! was reached or the absolute residual is below the scaled tolerance.
! 
!          ! Solve upper triangular system
!          gg(initer) = gg(initer) / hessenberg(initer,initer)
!          do ii=2,initer
!             k  = initer - ii + 1
!             k1 = k + 1
!             t  = gg(k)
!             do j=k1,initer
!                t = t - hessenberg(k,j) * gg(j)
!             end do
!             gg(k) = t / hessenberg(k,k)
!          end do
! 
!          ! Form linear combination to get solution.
!          do it=1,initer
!             t = gg(it)
! 
!             do j=j0,jn
!                do i=i0,in
!                   solution(i, j) = solution(i, j) + t * ww(i, j, it)
!                end do
!             end do
! 
!          end do
! 
!          outiter = outiter + 1
! 
!          if (norm_residual <= relative_tolerance .or. outiter > maxouter) then
!             return
!          end if
! 
!          ! Solution is not convergent : compute residual vector and continue.
!          do it=1,initer
!             jj = nextit - it + 1
!             gg(jj-1) = -rot_sin(jj-1) * gg(jj)
!             gg(jj)   =  rot_cos(jj-1) * gg(jj)
!          end do
! 
!          do it=1,nextit
!             t = gg(it)
!             if (it == 1) then
!                t = t - 1.d0
!             end if
! 
!             do j=j0,jn
! !DIR$ SIMD
!                do i=i0,in
!                   vv(i, j, 1) = vv(i, j, 1) + t * vv(i, j, it)
!                end do
!             end do
! 
!          end do
! 
      end do ! end of outer (restarts) loop
! 
!       return
      end subroutine fgmres

!=======================================================================

! PGMRES: Right-preconditioned generalized minimum residual method (with restarts). 
! Solves A x = b using GMRES with a right preconditioner
!
! authors: Stéphane Gaudreault, Abdessamad Qaddouri, Philippe Blain, ECCC

      subroutine pgmres()
         
      end subroutine pgmres

!=======================================================================

end module ice_krylov

!=======================================================================
