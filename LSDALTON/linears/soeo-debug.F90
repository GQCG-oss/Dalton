!> @file
!> Contains the soeo_debug module
!
!> brief Routines used in the debugging of SOEO (mainly getting the hessian)
!> author C. Nygaard
!> date 2010
module soeo_debug

use precision
use soeo_typedef
use matrix_util
use soeo_matop
use soeo_util
use soeo_redspace, only: soeo_projectout_gradplus
use Fock_evaluator
use files

implicit none

!contains
! soeo_get_full_hessian      (gets the full hessian from linear transformation)
! soeo_finite_difference     (gets the full hessian from finite difference)
! soeo_fd_get_energy         (for use in finite difference)
! soeo_fd_add_delta_to_K     (for use in finite difference)
! soeo_fd_add_delta_to_theta (for use in finite difference)

contains

!> \brief Finds the minimum eigenvalue of a matrix Ain
!> \author C. Nygaard
!> \date 2010-06-30
!> \param iter Dimension of Ain
!> \param Ain The matrix to which we want to find mineval
!> \param mineval The minimum eigenvalue of Ain
!=======================================================================
subroutine soeo_debug_find_mineval (lupri, iter, Ain, mineval, minevec)

implicit none

!I/O
integer, intent(in)      :: iter, lupri
real(realk), intent(in)  :: Ain(:,:)
real(realk), intent(out) :: mineval, minevec(:)
!Other
integer                  :: i, j, IERR
integer                  :: IV1(iter)
real(realk)              :: er(iter), ei(iter), X(iter,iter), FV1(iter), FV2(iter)
real(realk)              :: A(iter,iter), tmp
logical                  :: symmetric

A = Ain(1:iter,1:iter)

symmetric = .true.
do i=1,iter
  do j=1,i-1
    tmp = A(i,j) - A(j,i)
    if (tmp > 1.0E-5_realk .or. tmp < -1.0E-05_realk) then
      symmetric = .false.
    endif
  enddo
enddo

ei = 0.0E0_realk
if (symmetric) then
  call RS (iter, iter, A, er, 1, X, FV1, FV2, IERR)
else
  call lsquit ('redspace is not symmetric, it should be!', lupri)
  call RG (iter, iter, A, er, ei, 1, X, IV1, FV1, IERR)
endif

do i=1,iter
  if (ei(i) > 1.0E-3_realk .or. ei(i) < -1.0E-3_realk) then
    print *, "Redspace has complex eigenvalues!"
    call lsquit ("Redspace has complex eigenvalues!", lupri)
  endif
enddo
if (IERR /= 0) then
  print *, "Error in RG / RS (soeo_find_mineval)"
  call lsquit ("Error in RG / RS",lupri)
endif

mineval = 1.0E5_realk
do i=1,iter
  if (er(i) < mineval) then
    mineval = er(i)
    minevec = X(:,i)
  endif
enddo

!debug
write (lupri, *) 'Eigenvalues for A:'
call ls_output (er, 1, iter, 1, 1, iter, 1, 1, lupri)
write (lupri, *) 'Corresponding eigenvectors for A:'
call ls_output (X, 1, iter, 1, iter, iter, iter, 1, lupri)
!end debug

end subroutine soeo_debug_find_mineval
!=======================================================================



!> \brief Gets full Hessian
!> \author C. Nygaard
!> \date 2010
!> \param soeo Contains all matrices and information needed for soeo
!> \param hessian The full Hessian
!=======================================================================
subroutine soeo_get_full_hessian (soeo, mu, hessian)

implicit none

!I/O
type(soeoItem), intent(in)   :: soeo
real(realk), intent(in)      :: mu
real(realk), intent(inout)   :: hessian(:,:)
!Other
integer                      :: i, j, ib, Nbast, hessdim, Kdim
type(matrix)                 :: e_mat, e_mat_vec, tmpmat, hessvec
type(matrix)                 :: e_vec, tmpvec
type(matrix)                 :: gradp
real(realk)                  :: tmp, diffnorm, afac
real(realk), pointer         :: elm(:), t(:)
integer                      :: ndim, counter

!Initialization
!-----------------------------------------------------------------------
if (soeo%cfg_unres) then
  call mem_alloc (elm,2)
  call mem_alloc (t,2)
else
  call mem_alloc (elm,1)
  call mem_alloc (t,1)
endif

hessdim = size(hessian(:,1))
if (soeo%cfg_unres) then
  hessdim = hessdim / 2
endif
Kdim = hessdim - soeo%space%Nact
ndim = size(soeo%mats%nfirst)

call mat_init (e_mat, soeo%space%Nbast, soeo%space%Nbast)
call mat_init (e_mat_vec, Kdim, 1)
call mat_init (tmpmat, soeo%space%Nbast, soeo%space%Nbast)
call mat_init (hessvec, Kdim, 1)
call mat_init (e_vec, soeo%space%Nact, 1)
call mat_init (tmpvec, soeo%space%Nact, 1)
call mat_init (gradp, soeo%space%Nact, 1)

call mat_zero (e_mat)
call mat_zero (e_mat_vec)
call mat_zero (e_vec)

hessian = 0.0E0_realk

call soeo_gradplus (soeo%space, soeo%cfg_unres, soeo%mats%oldtheta, gradp)
call soeo_get_gradient (soeo, tmpmat, tmpvec)
afac = mat_dotproduct (tmpvec, gradp)
afac = afac / mat_dotproduct (gradp, gradp)

afac = 0.0E0_realk
!-----------------------------------------------------------------------

!write (soeo%lupri, *) 'soeo_get_full_hessian started'
!write (soeo%lupri, *) 'hessdim =', hessdim, 'Kdim =', Kdim

!Setting up the matrix-part
!-----------------------------------------------------------------------
do i=1,2*Kdim
  if (i>Kdim) then
    if (.not. soeo%cfg_unres) exit
    ib = i-Kdim + hessdim
  else
    ib = i
  endif

  if (soeo%cfg_unres) then
    if (i > Kdim) then
      call mat_create_ab_elms (i-Kdim, 1, (/ 0.0E0_realk, 1.0E0_realk /), e_mat_vec)
    else
      call mat_create_ab_elms (i, 1, (/ 1.0E0_realk, 0.0E0_realk /), e_mat_vec)
    endif
  else
    call mat_create_elm (i, 1, 1.0E0_realk, e_mat_vec)
  endif
  call mat_vec_to_mat ('a', e_mat_vec, e_mat)

  call soeo_linear_transform (soeo, e_mat, e_vec, tmpmat, tmpvec)
!  call soeo_projectout_gradplus (soeo%cfg_grandcan, soeo%lupri, gradp, tmpvec)
  call soeo_daxpy (-mu, e_mat, e_vec, tmpmat, tmpvec)

  call mat_to_vec ('a', tmpmat, hessvec)
  do j=1,Kdim
    if (soeo%cfg_unres) then
      call mat_get_ab_elms (hessvec, j, 1, elm)
      hessian(j,ib) = elm(1)
      hessian(hessdim+j,ib) = elm(2)
    else
      call mat_get_elm (hessvec, j, 1, elm(1))
      hessian(j,ib) = elm(1)
    endif
  enddo
  do j=1,soeo%space%Nact
    if (soeo%cfg_unres) then
      call mat_get_ab_elms (tmpvec, j, 1, elm)
      hessian(Kdim+j,ib) = elm(1)
      hessian(hessdim+Kdim+j,ib) = elm(2)
    else
      call mat_get_elm (tmpvec, j, 1, elm(1))
      hessian(Kdim+j,ib) = elm(1)
    endif
  enddo
  call mat_zero (e_mat)
  call mat_zero (e_mat_vec)
enddo
!-----------------------------------------------------------------------

!Setting up the vector-part
!-----------------------------------------------------------------------
do i=1,ndim
  if (i>soeo%space%Nact) then
    if (.not. soeo%cfg_unres) exit
    ib = i - soeo%space%Nact + hessdim + Kdim
  else
    ib = i + Kdim
  endif

  call mat_zero (e_vec)
  if (i > soeo%space%Nact) then
    call mat_create_ab_elms (i-soeo%space%Nact, 1, (/ 0.0E0_realk, 1.0E0_realk /), e_vec)
  else
    call mat_create_ab_elms (i, 1, (/ 1.0E0_realk, 0.0E0_realk /), e_vec)
  endif
!  call soeo_projectout_gradplus (soeo%cfg_grandcan, soeo%lupri, gradp, e_vec)

  call soeo_linear_transform (soeo, e_mat, e_vec, tmpmat, tmpvec)
!  call soeo_projectout_gradplus (soeo%cfg_grandcan, soeo%lupri, gradp, tmpvec)
  call soeo_daxpy (-mu, e_mat, e_vec, tmpmat, tmpvec)

  call mat_to_vec ('a', tmpmat, hessvec)
  do j=1,Kdim
    if (soeo%cfg_unres) then
      call mat_get_ab_elms (hessvec, j, 1, elm)
      hessian(j,ib) = elm(1)
      hessian(hessdim+j,ib) = elm(2)
    else
      call mat_get_elm (hessvec, j, 1, elm(1))
      hessian(j,ib) = elm(1)
    endif
  enddo
  do j=1,soeo%space%Nact
    if (soeo%cfg_unres) then
      call mat_get_ab_elms (tmpvec, j, 1, elm)
      hessian(Kdim+j,ib) = elm(1)
      if (Kdim+j==ib) then
        call mat_get_ab_elms (soeo%mats%oldtheta, soeo%space%Nocc+j, 1, t)
        hessian(Kdim+j,ib) = hessian(Kdim+j,ib) + afac * 2.0E0_realk * &
                            &  (dcos(t(1))*dcos(t(1))&
                            & - dsin(t(1))*dsin(t(1)))
      endif
      hessian(hessdim+Kdim+j,ib) = elm(2)
      if (hessdim+Kdim+j==ib) then
        call mat_get_ab_elms (soeo%mats%oldtheta, soeo%space%Nocc+j, 1, t)
        hessian(Kdim+j,ib) = hessian(Kdim+j,ib) + afac * 2.0E0_realk * &
                            &  (dcos(t(2))*dcos(t(2))&
                            & - dsin(t(2))*dsin(t(2)))
      endif
    else
      call mat_get_elm (tmpvec, j, 1, elm(1))
      hessian(Kdim+j,ib) = elm(1)
      if (Kdim+j==ib) then
        call mat_get_elm (soeo%mats%oldtheta, soeo%space%Nocc+j, 1, t(1))
        hessian(Kdim+j,ib) = hessian(Kdim+j,ib) + afac * 2.0E0_realk * &
                            &  (dcos(t(1))*dcos(t(1))&
                            & - dsin(t(1))*dsin(t(1)))
      endif
    endif
  enddo
enddo
!-----------------------------------------------------------------------

!Checking for symmetry
!-----------------------------------------------------------------------
write (soeo%lupri, *)
write (soeo%lupri, *) "Is the full hessian symmetric?"
if (soeo%cfg_unres) then
  write (soeo%lupri, *) "(alpha-part only)"
endif
write (soeo%lupri, *) "======================================"

write (soeo%lupri, *) "Full hessian ="
call ls_output (hessian, 1, hessdim, 1, hessdim,&
           & hessdim, hessdim, 1, soeo%lupri)

        
counter = 0
diffnorm = 0.0E0_realk
do i=1,hessdim
  do j=1,hessdim
    tmp = hessian(i,j) - hessian(j,i)
    if (tmp > 1.0E-10_realk .or. tmp < -1.0E-10_realk) then
      diffnorm = diffnorm + tmp*tmp
      counter = counter + 1
    endif
  enddo
enddo
if (counter /= 0) then
  write (soeo%lupri, *) "No, the Hessian is NOT SYMMETRIC!!!"
  write (soeo%lupri, *) "Number of elements in hessian (not diagonal):",&
                  & hessdim*hessdim-hessdim
  write (soeo%lupri, *) "Number of different elements =", counter
  write (soeo%lupri, *) "Accumulated difference =", diffnorm
  write (soeo%lupri, *) "Average difference =", diffnorm/counter
  write (soeo%lupri, *) "======================================"
else
  write (soeo%lupri, *) "Yes, the Hessian is SYMMETRIC!"
endif

counter = 0
diffnorm = 0.0E0_realk
do i=1,Kdim
  do j=1,Kdim
    tmp = hessian(i,j) - hessian(j,i)
    if (tmp < -1.0E-10_realk .or. tmp > 1.0E-10_realk) then
      diffnorm = diffnorm + tmp*tmp
      counter = counter + 1
    endif
  enddo
enddo
if (counter /= 0) then
  write (soeo%lupri, *) "KK-part of the hessian is NOT SYMMETRIC!!!"
  write (soeo%lupri, *) "Number of elements in KK-part (not diagonal):",&
                  & Kdim*Kdim-Kdim
  write (soeo%lupri, *) "Number of different elements in KK-part:",&
                  & counter
  write (soeo%lupri, *) "Accumulated difference =", diffnorm
  write (soeo%lupri, *) "Average difference =", diffnorm/counter
endif

counter = 0
diffnorm = 0.0E0_realk
do i=Kdim+1,hessdim
  do j=Kdim+1,hessdim
    tmp = hessian(i,j) - hessian(j,i)
    if (tmp > 1.0E-10_realk .or. tmp < -1.0E-10_realk) then
      diffnorm = diffnorm + tmp*tmp
      counter = counter + 1
    endif
  enddo
enddo
if (counter /= 0) then
  write (soeo%lupri, *) "tt-part of the hessian is NOT SYMMETRIC!!!"
  write (soeo%lupri, *) "Number of elements in tt-part (not diagonal):",&
                  & soeo%space%Nact*soeo%space%Nact-soeo%space%Nact
  write (soeo%lupri, *) "Number of different elements in tt-part:",&
                  & counter
  write (soeo%lupri, *) "Accumulated difference =", diffnorm
  write (soeo%lupri, *) "Average difference =", diffnorm/counter
endif

counter = 0
diffnorm = 0.0E0_realk
do i=1,Kdim
  do j=Kdim+1,hessdim
    tmp = hessian(i,j) - hessian(j,i)
    if (tmp < -1.0E-10_realk .or. tmp > 1.0E-10_realk) then
      diffnorm = diffnorm + tmp*tmp
      counter = counter + 1
    endif
  enddo
enddo
if (counter /= 0) then
  write (soeo%lupri, *) "Kt-part of the hessian is NOT SYMMETRIC!!!"
  write (soeo%lupri, *) "Number of elements in Kt-part:",&
                  & Kdim*soeo%space%Nact
  write (soeo%lupri, *) "Number of different elements in Kt-part:",&
                  & counter
  write (soeo%lupri, *) "Accumulated difference =", diffnorm
  write (soeo%lupri, *) "Average difference =", diffnorm/counter
endif
write (soeo%lupri, *) "======================================"
!-----------------------------------------------------------------------

!Closing
!-----------------------------------------------------------------------
call mat_free (e_mat)
call mat_free (e_mat_vec)
call mat_free (tmpmat)
call mat_free (hessvec)
call mat_free (e_vec)
call mat_free (tmpvec)
call mat_free (gradp)

call mem_dealloc (elm)
call mem_dealloc (t)
!-----------------------------------------------------------------------

end subroutine soeo_get_full_hessian
!=======================================================================

!> \brief Gets full projected tt-hessian P(H + H+)P
!> \author C. Nygaard
!> \date Mar 20 2012
!> \param soeo Contains all matrices and information needed
!=======================================================================
subroutine soeo_get_full_tt_hessian (soeo)

implicit none

type(soeoItem), intent(in) :: soeo

type(matrix)               :: ei, ej, zeromat
type(matrix)               :: sj, smat, spj
type(matrix)               :: gmat, gvec, gp
integer                    :: i, j
real(realk)                :: sofac, fomine, somine, mine
real(realk), pointer       :: fominevec(:), sominevec(:), minevec(:)
real(realk), pointer       :: hess(:,:), fohess(:,:), sohess(:,:)
logical                    :: sym
real(realk)                :: tmp

if (soeo%cfg_unres) then
  !nothing
else
  call mem_alloc (hess,soeo%space%Nact,soeo%space%Nact)
  call mem_alloc (fohess,soeo%space%Nact,soeo%space%Nact)
  call mem_alloc (sohess,soeo%space%Nact,soeo%space%Nact)
  call mem_alloc (fominevec,soeo%space%Nact)
  call mem_alloc (sominevec,soeo%space%Nact)
  call mem_alloc (minevec,soeo%space%Nact)

  call mat_init (zeromat, soeo%space%Nbast, soeo%space%Nbast)
  call mat_zero (zeromat)
  call mat_init (ei, soeo%space%Nact, 1)
  call mat_init (ej, soeo%space%Nact, 1)
  call mat_init (smat, soeo%space%Nbast, soeo%space%Nbast)
  call mat_init (sj, soeo%space%Nact, 1)
  call mat_init (spj, soeo%space%Nact, 1)
  call mat_init (gmat, soeo%space%Nbast, soeo%space%Nbast)
  call mat_init (gvec, soeo%space%Nact, 1)
  call mat_init (gp, soeo%space%Nact, 1)

  call soeo_get_gradient (soeo, gmat, gvec)
  call soeo_gradplus (soeo%space, soeo%cfg_unres, soeo%mats%oldtheta, gp)
  sofac = -1.0E0_realk * mat_dotproduct (gp, gvec) / mat_sqnorm2 (gp)

  write (soeo%lupri, *)
  write (soeo%lupri, *) 'g =', mat_sqnorm2 (gvec)
  call mat_print (gvec, 1, soeo%space%Nact, 1, 1, soeo%lupri)
  write (soeo%lupri, *) 'gp =', mat_sqnorm2 (gp)
  call mat_print (gp, 1, soeo%space%Nact, 1, 1, soeo%lupri)
  call soeo_projectout_gradplus (soeo%cfg_grandcan, soeo%lupri, gp, gvec)
  write (soeo%lupri, *) 'Pg =', mat_sqnorm2 (gvec)
  call mat_print (gvec, 1, soeo%space%Nact, 1, 1, soeo%lupri)
  write (soeo%lupri, *)

  do j=1,soeo%space%Nact
    call mat_zero (ej)
    call mat_create_elm (j, 1, 1.0E0_realk, ej)
    call soeo_projectout_gradplus (soeo%cfg_grandcan, soeo%lupri, gp, ej)
    call soeo_linear_transform (soeo, zeromat, ej, smat, sj)
    call soeo_projectout_gradplus (soeo%cfg_grandcan, soeo%lupri, gp, sj)
    call soeo_transformplus (soeo%space, soeo%cfg_unres, soeo%cfg_grandcan,&
                           & soeo%mats%oldtheta, ej, spj)
    call soeo_projectout_gradplus (soeo%cfg_grandcan, soeo%lupri, gp, spj)
    do i=1,soeo%space%Nact
      call mat_zero (ei)
      call mat_create_elm (i, 1, 1.0E0_realk, ei)
      call soeo_projectout_gradplus (soeo%cfg_grandcan, soeo%lupri, gp, ei)
      fohess(i,j) = mat_dotproduct (ei, sj)
      sohess(i,j) = sofac * mat_dotproduct (ei, spj)
    enddo
  enddo

  write (soeo%lupri, *) '------------------------------------------------'

  sym = .true.
  do i=1,soeo%space%Nact
    do j=1,soeo%space%Nact
      tmp = fohess(i,j) - fohess(j,i)
      tmp = dsqrt(tmp*tmp)
      if (tmp > 1.0E-5_realk) then
        sym = .false.
      endif
    enddo
  enddo
  call soeo_debug_find_mineval (soeo%lupri, soeo%space%Nact, fohess, fomine, fominevec)
  write (soeo%lupri, *) 'First order tt-hessian:'
  write (soeo%lupri, *) 'Symmetric?', sym
  call ls_output (fohess, 1, soeo%space%Nact, 1, soeo%space%Nact, soeo%space%Nact, soeo%space%Nact, 1, soeo%lupri)
  write (soeo%lupri, *) 'Min eval =', fomine
  write (soeo%lupri, *) 'Crsp evec ='
  call ls_output (fominevec, 1, soeo%space%Nact, 1, 1, soeo%space%Nact, 1, 1, soeo%lupri)

  sym = .true.
  do i=1,soeo%space%Nact
    do j=1,soeo%space%Nact
      tmp = sohess(i,j) - sohess(j,i)
      tmp = dsqrt(tmp*tmp)
      if (tmp > 1.0E-5_realk) then
        sym = .false.
      endif
    enddo
  enddo
  call soeo_debug_find_mineval (soeo%lupri, soeo%space%Nact, sohess, somine, sominevec)
  write (soeo%lupri, *) 'Second order tt-hessian:'
  write (soeo%lupri, *) 'Symmetric?', sym
  call ls_output (sohess, 1, soeo%space%Nact, 1, soeo%space%Nact, soeo%space%Nact, soeo%space%Nact, 1, soeo%lupri)
  write (soeo%lupri, *) 'Min eval =', somine
  write (soeo%lupri, *) 'Crsp evec ='
  call ls_output (sominevec, 1, soeo%space%Nact, 1, 1, soeo%space%Nact, 1, 1, soeo%lupri)

  hess = fohess + sohess
  sym = .true.
  do i=1,soeo%space%Nact
    do j=1,soeo%space%Nact
      tmp = hess(i,j) - hess(j,i)
      tmp = dsqrt(tmp*tmp)
      if (tmp > 1.0E-5_realk) then
        sym = .false.
      endif
    enddo
  enddo
  call soeo_debug_find_mineval (soeo%lupri, soeo%space%Nact, hess, mine, minevec)
  write (soeo%lupri, *) 'Total order tt-hessian:'
  write (soeo%lupri, *) 'Symmetric?', sym
  call ls_output (hess, 1, soeo%space%Nact, 1, soeo%space%Nact, soeo%space%Nact, soeo%space%Nact, 1, soeo%lupri)
  write (soeo%lupri, *) 'Min eval =', mine
  write (soeo%lupri, *) 'Crsp evec ='
  call ls_output (minevec, 1, soeo%space%Nact, 1, 1, soeo%space%Nact, 1, 1, soeo%lupri)

  write (soeo%lupri, *) '------------------------------------------------'
  
  call mat_free (zeromat)
  call mat_free (ei)
  call mat_free (ej)
  call mat_free (smat)
  call mat_free (sj)
  call mat_free (spj)
  call mat_free (gmat)
  call mat_free (gvec)
  call mat_free (gp)

  call mem_dealloc (hess)
  call mem_dealloc (fohess)
  call mem_dealloc (sohess)
  call mem_dealloc (fominevec)
  call mem_dealloc (sominevec)
  call mem_dealloc (minevec)
endif

end subroutine soeo_get_full_tt_hessian
!=======================================================================


!> \brief Gets full gradient and Hessian by Finite Difference
!> \author C. Nygaard
!> \date 2010-10-05
!> \param soeo Contains all matrices and information needed for soeo
!> \param gradient The gradient
!> \param hessian The full Hessian
!=======================================================================
subroutine soeo_finite_difference (soeo, mu, gradient, hessian)

implicit none

type(soeoItem), intent(in)   :: soeo
real(realk), intent(in)      :: mu
real(realk), intent(inout)   :: gradient(:), hessian(:,:)

type(matrix)                 :: theta, K
real(realk)                  :: delta, xnorm
real(realk)                  :: Ep1, Em1, Ep2, Em2
real(realk)                  :: Ep1p1, Ep1m1, Em1p1, Em1m1
real(realk)                  :: Ep2p2, Ep2m2, Em2p2, Em2m2
real(realk)                  :: G1, G2
integer                      :: hessdim, Kdim, ndim
integer                      :: i, j, row, col

hessdim = size(gradient)
if (soeo%cfg_unres) then
  hessdim = hessdim / 2
endif
Kdim = hessdim - soeo%space%Nact
ndim = size(soeo%mats%nfirst)

call mat_init (theta, soeo%space%Nbast, 1)
call mat_init (K, soeo%space%Nbast, soeo%space%Nbast)

delta = 0.005E0_realk

!Matrix-part
!-----------------------------------------------------------------------
do i=1,2*Kdim
  if (i>Kdim) then
    if (.not. soeo%cfg_unres) exit
    row = hessdim + i - Kdim
  else
    row = i
  endif

  !Energies with change in i'th variable in K
  !---------------------------------------------------------------------
  call mat_assign(theta,soeo%mats%oldtheta)

  call mat_zero (K)
  call soeo_fd_add_delta_to_K (K, delta, i, soeo%cfg_unres)
  xnorm = delta
  call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                         & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Ep1)

  call mat_zero (K)
  call soeo_fd_add_delta_to_K (K, -delta, i, soeo%cfg_unres)
  xnorm = delta
  call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                         & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Em1)

  call mat_zero (K)
  call soeo_fd_add_delta_to_K (K, 2.0E0_realk*delta, i, soeo%cfg_unres)
  xnorm = 2.0E0_realk*delta
  call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                         & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Ep2)

  call mat_zero (K)
  call soeo_fd_add_delta_to_K (K, -2.0E0_realk*delta, i, soeo%cfg_unres)
  xnorm = 2.0E0_realk*delta
  call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                         & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Em2)
  !---------------------------------------------------------------------

  !Matrix-part of the gradient
  ! gradient = (8Ep1 - 8Em1 - Ep2 + Em2)/(12delta)
  !---------------------------------------------------------------------
  gradient(row) = (8.0E0_realk*Ep1 - 8.0E0_realk*Em1 - Ep2 + Em2) / (12.0E0_realk*delta)
  !---------------------------------------------------------------------

  !Matrix-matrix-part
  !---------------------------------------------------------------------
  do j=1,2*Kdim
    if (j>Kdim) then
      if (.not. soeo%cfg_unres) exit
      col = hessdim + j - Kdim
    else
      col = j
    endif

    !Energies with change in i'th and j'th variable in K
    !-------------------------------------------------------------------
    call mat_assign(theta,soeo%mats%oldtheta)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, delta, i, soeo%cfg_unres)
    call soeo_fd_add_delta_to_K (K, delta, j, soeo%cfg_unres)
    if (i==j) then
      xnorm = 2.0E0_realk*delta
    else
      xnorm = sqrt(2.0E0_realk)*delta
    endif
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                         & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Ep1p1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, delta, i, soeo%cfg_unres)
    call soeo_fd_add_delta_to_K (K, -delta, j, soeo%cfg_unres)
    if (i==j) then
      xnorm = 0.0E0_realk
    else
      xnorm = sqrt(2.0E0_realk)*delta
    endif
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                         & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Ep1m1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -delta, i, soeo%cfg_unres)
    call soeo_fd_add_delta_to_K (K, delta, j, soeo%cfg_unres)
    if (i==j) then
      xnorm = 0.0E0_realk
    else
      xnorm = sqrt(2.0E0_realk)*delta
    endif
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                         & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Em1p1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -delta, i, soeo%cfg_unres)
    call soeo_fd_add_delta_to_K (K, -delta, j, soeo%cfg_unres)
    if (i==j) then
      xnorm = 2.0E0_realk*delta
    else
      xnorm = sqrt(2.0E0_realk)*delta
    endif
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                         & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Em1m1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, 2.0E0_realk*delta, i, soeo%cfg_unres)
    call soeo_fd_add_delta_to_K (K, 2.0E0_realk*delta, j, soeo%cfg_unres)
    if (i==j) then
      xnorm = 4.0E0_realk*delta
    else
      xnorm = sqrt(8.0E0_realk)*delta
    endif
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                         & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Ep2p2)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, 2.0E0_realk*delta, i, soeo%cfg_unres)
    call soeo_fd_add_delta_to_K (K, -2.0E0_realk*delta, j, soeo%cfg_unres)
    if (i==j) then
      xnorm = 0.0E0_realk
    else
      xnorm = sqrt(8.0E0_realk)*delta
    endif
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                         & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Ep2m2)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -2.0E0_realk*delta, i, soeo%cfg_unres)
    call soeo_fd_add_delta_to_K (K, 2.0E0_realk*delta, j, soeo%cfg_unres)
    if (i==j) then
      xnorm = 0.0E0_realk
    else
      xnorm = sqrt(8.0E0_realk)*delta
    endif
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                         & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Em2p2)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -2.0E0_realk*delta, i, soeo%cfg_unres)
    call soeo_fd_add_delta_to_K (K, -2.0E0_realk*delta, j, soeo%cfg_unres)
    if (i==j) then
      xnorm = 4.0E0_realk*delta
    else
      xnorm = sqrt(8.0E0_realk)*delta
    endif
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                         & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Em2m2)
    !-------------------------------------------------------------------

    !Matrix-matrix-part of the Hessian
    ! hessian = (16G1 - G2)/(48delta^2)
    !  G1 = Ep1p1 - Ep1m1 - Em1p1 + Em1m1
    !  G2 = Ep2p2 - Ep2m2 - Em2p2 + Em2m2
    !-------------------------------------------------------------------
    G1 = Ep1p1 - Ep1m1 - Em1p1 + Em1m1
    G2 = Ep2p2 - Ep2m2 - Em2p2 + Em2m2
    hessian(row,col) = (16.0E0_realk*G1 - G2) / (48.0E0_realk*delta*delta)
    !-------------------------------------------------------------------
  enddo
  !---------------------------------------------------------------------

  !Matrix-vector-part
  do j=1,ndim
    if (j>soeo%space%Nact) then
      if (.not. soeo%cfg_unres) exit
      col = j - soeo%space%Nact + hessdim + Kdim
    else
      col = j + Kdim
    endif

    !Energies with change in i'th variable in K
    !                    and j'th variable in theta
    !-------------------------------------------------------------------
    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, delta, i, soeo%cfg_unres)
    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, delta, j,&
                                   & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    xnorm = sqrt(2.0E0_realk)*delta
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Ep1p1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, delta, i, soeo%cfg_unres)
    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, -delta, j,&
                                   & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    xnorm = sqrt(2.0E0_realk)*delta
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Ep1m1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -delta, i, soeo%cfg_unres)
    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, delta, j,&
                                   & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    xnorm = sqrt(2.0E0_realk)*delta
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Em1p1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -delta, i, soeo%cfg_unres)
    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, -delta, j,&
                                   & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    xnorm = sqrt(2.0E0_realk)*delta
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Em1m1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, 2.0E0_realk*delta, i, soeo%cfg_unres)
    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, 2.0E0_realk*delta, j,&
                                   & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    xnorm = sqrt(8.0E0_realk)*delta
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Ep2p2)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, 2.0E0_realk*delta, i, soeo%cfg_unres)
    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, -2.0E0_realk*delta, j,&
                                   & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    xnorm = sqrt(8.0E0_realk)*delta
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Ep2m2)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -2.0E0_realk*delta, i, soeo%cfg_unres)
    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, 2.0E0_realk*delta, j,&
                                   & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    xnorm = sqrt(8.0E0_realk)*delta
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Em2p2)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -2.0E0_realk*delta, i, soeo%cfg_unres)
    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, -2.0E0_realk*delta, j,&
                                   & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    xnorm = sqrt(8.0E0_realk)*delta
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Em2m2)
    !-------------------------------------------------------------------

    !Matrix-vector-part of the Hessian
    ! hessian = (16G1 - G2)/(48delta^2)
    !  G1 = Ep1p1 - Ep1m1 - Em1p1 + Em1m1
    !  G2 = Ep2p2 - Ep2m2 - Em2p2 + Em2m2
    !-------------------------------------------------------------------
    G1 = Ep1p1 - Ep1m1 - Em1p1 + Em1m1
    G2 = Ep2p2 - Ep2m2 - Em2p2 + Em2m2
    hessian(row,col) = (16.0E0_realk*G1 - G2) / (48.0E0_realk*delta*delta)
    !-------------------------------------------------------------------
  enddo
enddo
!-----------------------------------------------------------------------

!Vector-part
!-----------------------------------------------------------------------
call mat_zero (K)
do i=1,ndim
  if (i>soeo%space%Nact) then
    if (.not. soeo%cfg_unres) exit
    row = i - soeo%space%Nact + hessdim + Kdim
  else
    row = i + Kdim
  endif

  !Energies with change in i'th variable in theta
  !---------------------------------------------------------------------
  call mat_zero (K)

  call mat_assign(theta,soeo%mats%oldtheta)
  call soeo_fd_add_delta_to_theta (theta, delta, i,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
  xnorm = delta
  call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                         & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Ep1)

  call mat_assign(theta,soeo%mats%oldtheta)
  call soeo_fd_add_delta_to_theta (theta, -delta, i,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
  xnorm = delta
  call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                         & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Em1)

  call mat_assign(theta,soeo%mats%oldtheta)
  call soeo_fd_add_delta_to_theta (theta, 2.0E0_realk*delta, i,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
  xnorm = 2.0E0_realk*delta
  call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                         & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Ep2)

  call mat_assign(theta,soeo%mats%oldtheta)
  call soeo_fd_add_delta_to_theta (theta, -2.0E0_realk*delta, i,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
  xnorm = 2.0E0_realk*delta
  call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                         & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Em2)
  !---------------------------------------------------------------------

  !Vector-part of the gradient
  ! gradient = (8Ep1 - 8Em1 - Ep2 + Em2)/(12delta)
  !---------------------------------------------------------------------
  gradient(row) = (8.0E0_realk*Ep1 - 8.0E0_realk*Em1 - Ep2 + Em2) / (12.0E0_realk*delta)
  !---------------------------------------------------------------------

  !Vector-matrix-part
  !---------------------------------------------------------------------
  do j=1,2*Kdim
    if (j>Kdim) then
      if (.not. soeo%cfg_unres) exit
      col = hessdim + j - Kdim
    else
      col = j
    endif

    !Energies with change in i'th variable in theta
    !                    and j'th variable in K
    !-------------------------------------------------------------------
    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, delta, j, soeo%cfg_unres)
    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, delta, i,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    xnorm = sqrt(2.0E0_realk)*delta
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Ep1p1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, delta, j, soeo%cfg_unres)
    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, -delta, i,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    xnorm = sqrt(2.0E0_realk)*delta
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Ep1m1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -delta, j, soeo%cfg_unres)
    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, delta, i,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    xnorm = sqrt(2.0E0_realk)*delta
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Em1p1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -delta, j, soeo%cfg_unres)
    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, -delta, i,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    xnorm = sqrt(2.0E0_realk)*delta
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Em1m1)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, 2.0E0_realk*delta, j, soeo%cfg_unres)
    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, 2.0E0_realk*delta, i,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    xnorm = sqrt(8.0E0_realk)*delta
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Ep2p2)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, 2.0E0_realk*delta, j, soeo%cfg_unres)
    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, -2.0E0_realk*delta, i,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    xnorm = sqrt(8.0E0_realk)*delta
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Ep2m2)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -2.0E0_realk*delta, j, soeo%cfg_unres)
    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, 2.0E0_realk*delta, i,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    xnorm = sqrt(8.0E0_realk)*delta
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Em2p2)

    call mat_zero (K)
    call soeo_fd_add_delta_to_K (K, -2.0E0_realk*delta, j, soeo%cfg_unres)
    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, -2.0E0_realk*delta, i,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    xnorm = sqrt(8.0E0_realk)*delta
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Em2m2)
    !-------------------------------------------------------------------

    !Vector-matrix-part of the Hessian
    ! hessian = (16G1 - G2)/(48delta^2)
    !  G1 = Ep1p1 - Ep1m1 - Em1p1 + Em1m1
    !  G2 = Ep2p2 - Ep2m2 - Em2p2 + Em2m2
    !-------------------------------------------------------------------
    G1 = Ep1p1 - Ep1m1 - Em1p1 + Em1m1
    G2 = Ep2p2 - Ep2m2 - Em2p2 + Em2m2
    hessian(row,col) = (16.0E0_realk*G1 - G2) / (48.0E0_realk*delta*delta)
    !-------------------------------------------------------------------
  enddo
  !---------------------------------------------------------------------

  !Vector-vector-part
  !---------------------------------------------------------------------
  do j=1,ndim
    if (j>soeo%space%Nact) then
      if (.not. soeo%cfg_unres) exit
      col = j - soeo%space%Nact + hessdim + Kdim
    else
      col = j + Kdim
    endif

    !Energies with change in i'th and j'th variable in theta
    !-------------------------------------------------------------------
    call mat_zero (K)

    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, delta, i,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    call soeo_fd_add_delta_to_theta (theta, delta, j,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    if (i == j) then
      xnorm = 2.0E0_realk*delta
    else
      xnorm = sqrt(2.0E0_realk)*delta
    endif
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Ep1p1)

    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, delta, i,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    call soeo_fd_add_delta_to_theta (theta, -delta, j,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    if (i == j) then
      xnorm = 0.0E0_realk
    else
      xnorm = sqrt(2.0E0_realk)*delta
    endif
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Ep1m1)

    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, -delta, i,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    call soeo_fd_add_delta_to_theta (theta, delta, j,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    if (i == j) then
      xnorm = 0.0E0_realk
    else
      xnorm = sqrt(2.0E0_realk)*delta
    endif
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Em1p1)

    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, -delta, i,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    call soeo_fd_add_delta_to_theta (theta, -delta, j,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    if (i == j) then
      xnorm = 2.0E0_realk*delta
    else
      xnorm = sqrt(2.0E0_realk)*delta
    endif
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Em1m1)

    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, 2.0E0_realk*delta, i,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    call soeo_fd_add_delta_to_theta (theta, 2.0E0_realk*delta, j,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    if (i == j) then
      xnorm = 4.0E0_realk*delta
    else
      xnorm = sqrt(8.0E0_realk)*delta
    endif
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Ep2p2)

    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, 2.0E0_realk*delta, i,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    call soeo_fd_add_delta_to_theta (theta, -2.0E0_realk*delta, j,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    if (i == j) then
      xnorm = 0.0E0_realk
    else
      xnorm = sqrt(8.0E0_realk)*delta
    endif
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Ep2m2)

    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, -2.0E0_realk*delta, i,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    call soeo_fd_add_delta_to_theta (theta, 2.0E0_realk*delta, j,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    if (i == j) then
      xnorm = 0.0E0_realk
    else
      xnorm = sqrt(8.0E0_realk)*delta
    endif
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Em2p2)

    call mat_assign(theta,soeo%mats%oldtheta)
    call soeo_fd_add_delta_to_theta (theta, -2.0E0_realk*delta, i,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    call soeo_fd_add_delta_to_theta (theta, -2.0E0_realk*delta, j,&
                                 & soeo%space%Nocc, soeo%space%Nact, soeo%cfg_unres)
    if (i == j) then
      xnorm = 4.0E0_realk*delta
    else
      xnorm = sqrt(8.0E0_realk)*delta
    endif
    call soeo_fd_get_energy (soeo%cfg_unres, K, theta, soeo%mats%C, soeo%mats%S,&
                           & xnorm, mu, soeo%settings%trust,&
                         & soeo%space%Nelec, Em2m2)
    !-------------------------------------------------------------------

    !Vector-vector-part of the Hessian
    ! hessian = (16G1 - G2)/(48delta^2)
    !  G1 = Ep1p1 - Ep1m1 - Em1p1 + Em1m1
    !  G2 = Ep2p2 - Ep2m2 - Em2p2 + Em2m2
    !-------------------------------------------------------------------
    G1 = Ep1p1 - Ep1m1 - Em1p1 + Em1m1
    G2 = Ep2p2 - Ep2m2 - Em2p2 + Em2m2
    hessian(row,col) = (16.0E0_realk*G1 - G2) / (48.0E0_realk*delta*delta)
    !-------------------------------------------------------------------
  enddo
  !---------------------------------------------------------------------  
enddo
!-----------------------------------------------------------------------

call mat_free (theta)
call mat_free (K)

end subroutine soeo_finite_difference
!=======================================================================


!> \brief Subroutine that gets the energy as function of K and theta
!> \author C. Nygaard
!> \date 2010-10-06
!> \param unres True if calculation is unrestricted
!> \param K The antisymmetric rotation matrix
!> \param theta The vector of occupation-angles
!> \param Cin The MO coefficient matrix
!> \param S The AO overlap matrix
!> \param E The energy for given K and theta
!=======================================================================
subroutine soeo_fd_get_energy (unres, K, theta, Cin, S, xnorm, mu, trust, Ne, E)

implicit none

logical, intent(in)      :: unres
type(matrix), intent(in) :: K, theta, Cin, S
real(realk), intent(in)  :: xnorm, mu, trust
integer, intent(in)      :: Ne
real(realk), intent(out) :: E

integer                  :: Nbast, i
real(realk), pointer     :: tmp(:)
type(matrix)             :: C, Dao, Dmo, Fao, tmptheta

if (unres) then
  call mem_alloc (tmp,2)
else
  call mem_alloc (tmp,1)
endif

Nbast = Cin%nrow

call mat_init (C, Nbast, Nbast)
call mat_init (Dao, Nbast, Nbast)
call mat_init (Dmo, Nbast, Nbast)
call mat_init (Fao, Nbast, Nbast)
call mat_init (tmptheta, Nbast, 1)

call mat_assign(tmptheta,theta)

call mat_assign(C,Cin)
call soeo_getnew_C (K, C)
call soeo_getnew_Dmo (tmptheta, Dmo, unres)
call util_MO_to_AO_2 (S, C, Dmo, Dao, .false.)
call FCK_get_Fock (Dao, Fao, E)

!E = E - 0.5E0_realk*mu*(xnorm*xnorm - trust*trust)
E = E - 0.5E0_realk*mu*xnorm*xnorm

call mat_free (C)
call mat_free (Dao)
call mat_free (Dmo)
call mat_free (Fao)
call mat_free (tmptheta)

call mem_dealloc (tmp)

end subroutine soeo_fd_get_energy
!=======================================================================

!> /brief Adds delta to the i'th variable in the antisymmetric matrix K
!> /author C. Nygaard
!> /date 2010-10-06
!> \param K The antisymmetric rotation matrix
!> \param delta A small perturbation to K
!> \param i Defines what element in K should be perturbed
!> \param unres Logical, true if calculation is unrestricted
!=======================================================================
subroutine soeo_fd_add_delta_to_K (K, delta, i, unres)

implicit none

type(matrix), intent(inout) :: K
real(realk), intent(in)     :: delta
integer, intent(in)         :: i
logical, intent(in)         :: unres

integer                     :: N, Kdim
type(matrix)                :: Kvec
real(realk), pointer        :: Kelm(:)

if (unres) then
  call mem_alloc (Kelm,2)
else
  call mem_alloc (Kelm,1)
endif
N = K%nrow
Kdim = (N*(N-1))/2
call mat_init (Kvec, Kdim, 1)

call mat_to_vec ('a', K, Kvec)
if (i <= Kdim) then
  call mat_get_ab_elms (Kvec, i, 1, Kelm)
  Kelm(1) = Kelm(1) + delta
  call mat_create_ab_elms (i, 1, Kelm, Kvec)
elseif (i > Kdim .and. i <= 2*Kdim .and. unres) then
  call mat_get_ab_elms (Kvec, i-Kdim, 1, Kelm)
  Kelm(2) = Kelm(2) + delta
  call mat_create_ab_elms (i-Kdim, 1, Kelm, Kvec)
else
  print *, 'i =', i, 'Kdim =', Kdim, 'unres =', unres
  stop "K does not contain an element in the i'th place"
endif
call mat_vec_to_mat ('a', Kvec, K)

call mat_free (Kvec)
call mem_dealloc (Kelm)

end subroutine soeo_fd_add_delta_to_K
!=======================================================================

!> \brief Adds delta to the i'th variable in the vector theta
!> \author C. Nygaard
!> \date 2010-10-06
!> \param theta The vector of occupation angles
!> \param delta A small perturbation in theta
!> \param i Defines what angle should be perturbed
!> \param Nocc Size of occupied space
!> \param Nact Size of active space
!> \param unres Logical, true if calculation is unrestricted
!=======================================================================
subroutine soeo_fd_add_delta_to_theta (theta, delta, i, Nocc, Nact, unres)

implicit none

type(matrix), intent(inout) :: theta
real(realk), intent(in)     :: delta
integer, intent(in)         :: i, Nocc, Nact
logical, intent(in)         :: unres

real(realk), pointer        :: thetaelm(:)

if (unres) then
  call mem_alloc (thetaelm,2)
else
  call mem_alloc (thetaelm,1)
endif

if (i <= Nact) then
  call mat_get_ab_elms (theta, Nocc+i, 1, thetaelm)
  thetaelm(1) = thetaelm(1) + delta
  call mat_create_ab_elms (Nocc+i, 1, thetaelm, theta)
elseif (i > Nact .and. i <= 2*Nact .and. unres) then
  call mat_get_ab_elms (theta, Nocc+i-Nact, 1, thetaelm)
  thetaelm(2) = thetaelm(2) + delta
  call mat_create_ab_elms (Nocc+i-Nact, 1, thetaelm, theta)
else
  print *, 'i =', i, 'Nact =', Nact, 'unres =', unres
  stop "theta does not contain an element in the i'th place"
endif

call mem_dealloc (thetaelm)

end subroutine soeo_fd_add_delta_to_theta
!=======================================================================


!> \brief Compare FD Hessian with full Hessian from linear transform
!> \author C. Nygaard
!> \date Nov 21. 2012
!> \param
!=======================================================================
subroutine soeo_compare_hessians (soeo)

implicit none

type(soeoItem), intent(in) :: soeo
real(realk), pointer       :: hessian(:,:), hessian2(:,:), gradient(:), gradient2(:)
real(realk)                :: mu

integer                    :: tmpint, i, j
real(realk)                :: tmp

mu = 0.0E0_realk
!tmpint = size of hessian
tmpint = (soeo%space%Nbast * (soeo%space%Nbast - 1) / 2) + soeo%space%Nact
if (soeo%cfg_unres) then
  tmpint = 2 * tmpint
endif
call mem_alloc(hessian,tmpint,tmpint)
call mem_alloc(gradient,tmpint)
call mem_alloc(hessian2,tmpint,tmpint)
call mem_alloc(gradient2,tmpint)

!!tmpdim = number of variables in kappa
!if (soeo%cfg_unres) then
!  tmpdim = tmpint/2
!else
!  tmpdim = tmpint
!endif
!tmpdim = tmpdim - soeo%space%Nact
!call mat_init (tmpvec, tmpdim, 1) !for storing kappas

!get full hessian from linear transform
call soeo_get_full_hessian (soeo, mu, hessian)
!write (soeo%lupri, *) 'hessian by linear transform'
!call ls_output (hessian, 1, tmpint, 1, tmpint, &
!           & tmpint, tmpint, 1, soeo%lupri)
!
!lueval = -1 ; luevecs = -1
!call lsopen (lueval, "evals", "UNKNOWN", "UNFORMATTED")
!call lsopen (luevecs, "evecs", "UNKNOWN", "UNFORMATTED")
!call lsclose (lueval, "KEEP")
!call lsclose (luevecs, "KEEP")
!call soeo_find_mineval (soeo%lupri, tmpint, hessian, tmp)
!  
!  !get gradient from linear transform
!  call soeo_get_gradient (soeo, grad_m, grad_v)
!  call mat_to_vec ('a', grad_m, tmpvec)
!  if (soeo%cfg_unres) then
!    allocate (step(2))
!  else
!    allocate (step(1))
!  endif
!  do i=1,tmpdim
!    call mat_get_ab_elms (tmpvec, i, 1, step)
!    gradient(i) = step(1)
!    if (soeo%cfg_unres) then
!      gradient(i+tmpint/2) = step(2)
!    endif
!  enddo
!  do i=1,soeo%space%Nact
!    call mat_get_ab_elms (grad_v, i, 1, step)
!    gradient(tmpdim+i) = step(1)
!    if (soeo%cfg_unres) then
!      gradient(tmpint/2+tmpdim+i) = step(2)
!    endif
!  enddo
!  write (soeo%lupri, *) 'gradient in iteration', iter
!  call ls_output (gradient, 1, tmpint, 1, 1, tmpint, 1, 1, soeo%lupri)
!  deallocate (step)
!  
!Finite difference Hessian and gradient
call soeo_finite_difference (soeo, mu, gradient2, hessian2)
!write (soeo%lupri, *) 'hessian from finite difference'
!call ls_output (hessian2, 1, tmpint, 1, tmpint, &
!           & tmpint, tmpint, 1, soeo%lupri)
!write (soeo%lupri, *) 'and the gradient'
!call ls_output (gradient2, 1, tmpint, 1, 1, tmpint, 1, 1, soeo%lupri)

write (soeo%lupri, *) 'hessian from lt:'
call ls_output (hessian, 1, tmpint, 1, tmpint, &
           & tmpint, tmpint, 1, soeo%lupri)

write (soeo%lupri, *) 'hessians from fd:'
call ls_output (hessian2, 1, tmpint, 1, tmpint, &
           & tmpint, tmpint, 1, soeo%lupri)

write (soeo%lupri, *) 'difference between hessians:'
hessian = hessian-hessian2
call ls_output (hessian, 1, tmpint, 1, tmpint, &
           & tmpint, tmpint, 1, soeo%lupri)
tmp = 0.0E0_realk
do i=1,tmpint
  do j=1,tmpint
    tmp = tmp + hessian(i,j)*hessian(i,j)
  enddo
enddo
tmp = sqrt(tmp)
write (soeo%lupri, *) 'diffnorm hessian', tmp
!  write (soeo%lupri, *) 'difference between gradients:'
!  gradient = gradient-gradient2
!  call ls_output (gradient, 1, tmpint, 1, 1, tmpint, 1, 1, soeo%lupri)
!  tmp=0.0E0_realk
!  do i=1,tmpint
!    tmp = tmp+gradient(i)*gradient(i)
!  enddo
!  tmp = sqrt(tmp)
!  write (soeo%lupri, *) 'diffnorm gradient', tmp
!  
!  write (soeo%lupri, *) 'mu =', mu
!  
!  
!  call mat_free (tmpvec)
  
call mem_dealloc(hessian)
call mem_dealloc(gradient)
call mem_dealloc(hessian2)
call mem_dealloc(gradient2)

end subroutine soeo_compare_hessians
!=======================================================================

end module soeo_debug
