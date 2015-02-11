module soeo_redspace

use soeo_typedef
use soeo_matop
use files
use soeo_util

use Fock_evaluator, only: FCK_get_fock

implicit none

!Contains:
!  soeo_solver
!  soeo_projectout_gradplus
!  soeo_solve_linear_system
!  soeo_solver_getX
!  soeo_solver_getres
!  soeo_find_mineval
!  soeo_binary_search
!  soeo_correct_deltatheta
!  soeo_label_orbitals
!  soeo_precond
!  soeo_cleanup_trialvector

contains

!> \brief Solves the SOEO-Newton equations in a reduced space
!> \author C. Nygaard
!> \date 2010
!> \param soeo Contains all information and matrices used in soeo
!> \param K The MO orbital rotation matrix
!> \param dt The change in the occupation angles
!> \param mu Level shift, output
!> \param gm Matrix-part of gradient
!> \param gv Vector-part of gradient
!> \param resp True if used for calculation of properties
!=======================================================================
subroutine soeo_solver (soeo, K, dt, mu, gm, gv, resp)

implicit none

!Input/output
type(soeoItem), intent(inout):: soeo
type(matrix), intent(inout)  :: K, dt
real(realk), intent(out)     :: mu
type(matrix), intent(in)     :: gm, gv
logical, intent(in)          :: resp

!Other:
integer                      :: iter
integer                      :: lub_m, lub_v, lus_m, lus_v, lusp
integer                      :: Nbast, Nact
integer                      :: i, j
logical                      :: err
!Gradient, trial-vectors and transformed trial-vectors:
type(matrix)                 :: grad_m, bi_m, biter_m, si_m, siter_m
type(matrix)                 :: grad_v, bi_v, biter_v, si_v, siter_v
type(matrix)                 :: gradp, spi, spiter
!Useful numbers:
real(realk)                  :: graddot, Apfac, alpha
real(realk)                  :: mine, minevec(soeo%settings%micromaxiter)
real(realk)                  :: thresh
!The reduced space to solve:
real(realk)                  :: Ared(soeo%settings%micromaxiter,soeo%settings%micromaxiter),&
                               & gradred(soeo%settings%micromaxiter),&
                               & xred(soeo%settings%micromaxiter)
!Solution and residual:
type(matrix)                 :: X_m, X_v, res_m, res_v, dX_m, dX_v
real(realk)                  :: resnorm
!For sorting out redundant orbital rotations:
character, pointer           :: orbs(:)
!Tmps:
real(realk), pointer         :: tmp1(:), tmp2(:), tmpAp(:,:)
real(realk)                  :: tmp, val1, val2
type(matrix)                 :: tmpmat1, tmpvec1, tmpmat2, tmpvec2, tmpmat3, tmpvec3,&
                              & tmpoldt, tmpDmo, tmpDao, tmpFao, tmpC
integer                      :: int1, int2, int3
logical :: OnMaster
OnMaster=.TRUE.
!Initializations:
Nbast = soeo%space%Nbast
Nact  = soeo%space%Nact
call mat_init (grad_m , Nbast, Nbast)
call mat_init (grad_v , Nact , 1    )
call mat_init (gradp  , Nact , 1    )
call mat_init (bi_m   , Nbast, Nbast)
call mat_init (bi_v   , Nact , 1    )
call mat_init (biter_m, Nbast, Nbast)
call mat_init (biter_v, Nact , 1    )
call mat_init (si_m   , Nbast, Nbast)
call mat_init (si_v   , Nact , 1    )
call mat_init (siter_m, Nbast, Nbast)
call mat_init (siter_v, Nact , 1    )
call mat_init (spi    , Nact , 1    )
call mat_init (spiter , Nact , 1    )
call mat_init (X_m    , Nbast, Nbast)
call mat_init (X_v    , Nact , 1    )
call mat_init (res_m  , Nbast, Nbast)
call mat_init (res_v  , Nact , 1    )
call mat_init (dX_m   , Nbast, Nbast)
call mat_init (dX_v   , Nact , 1    )

lub_m = -1 ; lub_v = -1 ; lus_m = -1 ; lus_v = -1 ; lusp = -1
call lsopen (lub_m, "bmats", "unknown", "UNFORMATTED")
call lsopen (lub_v, "bvecs", "unknown", "UNFORMATTED")
call lsopen (lus_m, "sigmamats", "unknown", "UNFORMATTED")
call lsopen (lus_v, "sigmavecs", "unknown", "UNFORMATTED")
call lsopen (lusp, "sigmaps", "unknown", "UNFORMATTED")

!Sort out the redundant orbital rotations
!-----------------------------------------------------------------------
!Create a vector with elements 'i' for inactive, 'a' for active and 
! 'v' for virtual orbitals
if (soeo%cfg_unres) then
  call mem_alloc (orbs,2*Nbast)
else
  call mem_alloc (orbs,Nbast)
endif
call soeo_label_orbitals (soeo%mats%Dmo, orbs, soeo%cfg_unres, soeo%lupri, soeo%prnt)

if (soeo%prnt>1) then
  int1=0 ; int2=0 ; int3=0
  do i=1,size(orbs)
    if (orbs(i) == 'i' .or. orbs(i) == 'I') int1=int1+1
    if (orbs(i) == 'a' .or. orbs(i) == 'A') int2=int2+1
    if (orbs(i) == 'v' .or. orbs(i) == 'V') int3=int3+1
  enddo
  write (soeo%lupri, *) 
  write (soeo%lupri, '("Orbital spaces in this Newton-iteration:")')
  write (soeo%lupri, '("Inactive: ", i6)') int1
  write (soeo%lupri, '("Active:   ", i6)') int2
  write (soeo%lupri, '("Virtual:  ", i6)') int3
  write (soeo%lupri, *) 
endif
!-----------------------------------------------------------------------

!Get the gradient (both the matrix-part and the vector-part)
!-----------------------------------------------------------------------
if (resp) then
  call mat_assign (grad_m, gm) ; call mat_assign (grad_v, gv)
else
  call soeo_get_gradient (soeo, grad_m, grad_v)
endif
call soeo_precond (soeo%space, grad_m, grad_v, orbs, soeo%cfg_unres)
call soeo_gradplus (soeo%space, soeo%cfg_unres, soeo%mats%oldtheta, gradp)
graddot = mat_dotproduct (gradp, grad_v)
if (mat_sqnorm2 (gradp) > 1.0E-5_realk) then
  Apfac = -1.0E0_realk * graddot / mat_sqnorm2 (gradp)
else
  Apfac = 0.0E0_realk
endif
call soeo_projectout_gradplus (soeo%cfg_grandcan, soeo%lupri, gradp, grad_v)
!-----------------------------------------------------------------------

!Setup the first trial-vector b1 = norm(pre(norm(grad)))
!-----------------------------------------------------------------------
call mat_assign(biter_m,grad_m) ; call mat_assign(biter_v,grad_v)
call soeo_precond (soeo%space, biter_m, biter_v, orbs, soeo%cfg_unres)
call soeo_normalize (biter_m, biter_v)
call soeo_projectout_gradplus (soeo%cfg_grandcan, soeo%lupri, gradp, biter_v)
call soeo_normalize (biter_m, biter_v)

thresh = soeo_norm (grad_m, grad_v)
if (thresh > 1.0E0_realk) then
  thresh = 0.1E0_realk * thresh
  if (thresh > 1.0E0_realk) then
    thresh = 1.0E0_realk
  endif
else
  thresh = 0.1E0_realk*thresh*thresh
endif
if (thresh < soeo%settings%microthresh .or. resp) then
  thresh = soeo%settings%microthresh
endif

rewind (lub_m) ; rewind (lub_v)
call mat_write_to_disk (lub_m, biter_m,OnMaster)
call mat_write_to_disk (lub_v, biter_v,OnMaster)
!-----------------------------------------------------------------------

!Building the reduced space (first entries):
!-----------------------------------------------------------------------
iter = 1
Ared = 0.0E0_realk ; gradred = 0.0E0_realk

call soeo_linear_transform (soeo, biter_m, biter_v, siter_m, siter_v)
call soeo_precond (soeo%space, siter_m, siter_v, orbs, soeo%cfg_unres)
call soeo_projectout_gradplus (soeo%cfg_grandcan, soeo%lupri, gradp, siter_v)

call soeo_transformplus (soeo%space, soeo%cfg_unres, soeo%cfg_grandcan,&
                       & soeo%mats%oldtheta, biter_v, spiter)
call soeo_projectout_gradplus (soeo%cfg_grandcan, soeo%lupri, gradp, spiter)

rewind (lus_m) ; rewind (lus_v) ; rewind (lusp)
call mat_write_to_disk (lus_m, siter_m,OnMaster)
call mat_write_to_disk (lus_v, siter_v,OnMaster)
call mat_write_to_disk (lusp, spiter,OnMaster)

Ared(1,1) = soeo_dotproduct (biter_m, biter_v, siter_m, siter_v) ! Ared
Ared(1,1) = Ared(1,1) + Apfac*mat_dotproduct(biter_v, spiter) ! + Aredplus

gradred(1) = -soeo_dotproduct (biter_m, biter_v, grad_m, grad_v) ! gradred
!-----------------------------------------------------------------------

mu = 0.0E0_realk
call mat_zero (X_m)
call mat_zero (X_v)
do

  write (soeo%lupri, *) 
  if (soeo%prnt>2) then
    write (soeo%lupri, '("Microiteration ", i3)') iter
    write (soeo%lupri, '("==============================")')
  endif

  !Find mu (so that step is inside trust region)
  !---------------------------------------------------------------------
  mu = 0.0E0_realk
  call soeo_find_mineval (soeo%lupri, iter, Ared, mine, minevec)
  if (soeo%prnt>2) then
    write (soeo%lupri, '("Step control:")')
    write (soeo%lupri, '("------------------------------")')
    write (soeo%lupri, '("Minimum eigenvalue of reduced Hessian:", f15.7)') mine
    write (soeo%lupri, *)
  endif

!  if (mine < 1.0E-3_realk) then
!
!    print *, 'Negative eigenvalue!!!'
!
!    !debug: examine the eigenvector corresponding to this eigenvalue
!    call mat_init (tmpmat1, soeo%space%Nbast, soeo%space%Nbast)
!    call mat_init (tmpvec1, soeo%space%Nact, 1)
!    call mat_init (tmpmat2, soeo%space%Nbast, soeo%space%Nbast)
!    call mat_init (tmpvec2, soeo%space%Nact, 1)
!    call mat_init (tmpmat3, soeo%space%Nbast, soeo%space%Nbast)
!    call mat_init (tmpvec3, soeo%space%Nact, 1)
!    call mat_init (tmpoldt, soeo%space%Nbast, 1)
!    call mat_init (tmpDmo, soeo%space%Nbast, soeo%space%Nbast)
!    call mat_init (tmpDao, soeo%space%Nbast, soeo%space%Nbast)
!    call mat_init (tmpFao, soeo%space%Nbast, soeo%space%Nbast)
!    call mat_init (tmpC, soeo%space%Nbast, soeo%space%Nbast)
!
!    !create eigenvector-X
!    call soeo_solver_getX (iter, lub_m, lub_v, lusp, minevec, gradp, tmpmat1, tmpvec1)
!
!    !create eigenvector-X without second-order contribution
!    call mat_zero (tmpmat3)
!    call mat_zero (tmpvec3)
!    rewind (lub_m) ; rewind (lub_v)
!    do i=1,iter
!      call mat_read_from_disk (lub_m, bi_m)
!      call mat_read_from_disk (lub_v, bi_v)
!
!      call soeo_daxpy (minevec(i), bi_m, bi_v, tmpmat3, tmpvec3)
!    enddo
!
!    write (soeo%lupri, *) 'Evec * gradp =', mat_dotproduct (tmpvec1, gradp)
!    write (soeo%lupri, *) 'Evec with so-correction ='
!    call mat_print (tmpvec1, 1, soeo%space%Nact, 1, 1, soeo%lupri)
!    write (soeo%lupri, *) 'Evec without so-correction ='
!    call mat_print (tmpvec3, 1, soeo%space%Nact, 1, 1, soeo%lupri)
!
!    write (soeo%lupri, *) '1DELTAN#    scal         deltaN'
!    write (soeo%lupri, *) '1EEE#       scal         E     '
!    write (soeo%lupri, *) '1PPP#       scal         dEpred'
!    write (soeo%lupri, *) '2DELTAN#    scal         deltaN'
!    write (soeo%lupri, *) '2EEE#       scal         E     '
!    write (soeo%lupri, *) '2PPP#       scal         dEpred'
!    do i=-200,200
!      !For X
!      !scale eigenvector
!      tmpmat2 = tmpmat3
!      tmpvec2 = tmpvec3
!      call soeo_scal (i/100.0E0_realk, tmpmat2, tmpvec2)
!      call soeo_transformplus (soeo%space, soeo%cfg_unres, soeo%cfg_grandcan,&
!                             & soeo%mats%oldtheta, tmpvec2, tmpvec1)
!      tmp = -0.5E0_realk * mat_dotproduct(tmpvec2,tmpvec1) / mat_sqnorm2(gradp)
!      call mat_daxpy (tmp, gradp, tmpvec2)
!
!      !calculate number of electrons
!      tmp = soeo_numberof_electrons (soeo%space%Nocc, soeo%mats%oldtheta, tmpvec2, soeo%cfg_unres)
!      tmp = tmp - mat_Tr(soeo%mats%Dmo)
!      write (soeo%lupri, *) '1DELTAN  ', i/100.0E0_realk, tmp
!
!      !calculate energy
!      tmpoldt = soeo%mats%oldtheta
!      call soeo_getnew_theta (tmpvec2, soeo%space, tmpoldt, soeo%cfg_unres)
!      tmpDmo = soeo%mats%Dmo
!      call soeo_getnew_Dmo (tmpoldt, tmpDmo, soeo%cfg_unres)
!      tmpC = soeo%mats%C
!      call soeo_getnew_C (tmpmat2, tmpC)
!      call util_MO_to_AO_2 (soeo%mats%S, tmpC, tmpDmo, tmpDao, .false.)
!      call FCK_get_Fock (tmpDao, tmpFao, tmp)
!      write (soeo%lupri, *) '1EEE  ', i/100.0E0_realk, tmp-soeo%Etotal
!
!      !calculate predicted energy
!      val1 = soeo%dEpred
!      call soeo_get_Epred (soeo, tmpmat2, tmpvec2)
!      val2 = soeo%dEpred
!      soeo%dEpred = val1
!      write (soeo%lupri, *) '1PPP  ', i/100.0E0_realk, val2
!
!      !For X without second-order correction
!      !scale eigenvector
!      tmpmat2 = tmpmat3
!      tmpvec2 = tmpvec3
!      call soeo_scal (i/100.0E0_realk, tmpmat2, tmpvec2)
!
!      !calculate number of electrons
!      tmp = soeo_numberof_electrons (soeo%space%Nocc, soeo%mats%oldtheta, tmpvec2, soeo%cfg_unres)
!      tmp = tmp - mat_Tr(soeo%mats%Dmo)
!      write (soeo%lupri, *) '2DELTAN  ', i/100.0E0_realk, tmp
!
!      !calculate energy
!      tmpoldt = soeo%mats%oldtheta
!      call soeo_getnew_theta (tmpvec2, soeo%space, tmpoldt, soeo%cfg_unres)
!      tmpDmo = soeo%mats%Dmo
!      call soeo_getnew_Dmo (tmpoldt, tmpDmo, soeo%cfg_unres)
!      tmpC = soeo%mats%C
!      call soeo_getnew_C (tmpmat2, tmpC)
!      call util_MO_to_AO_2 (soeo%mats%S, tmpC, tmpDmo, tmpDao, .false.)
!      call FCK_get_Fock (tmpDao, tmpFao, tmp)
!      write (soeo%lupri, *) '2EEE  ', i/100.0E0_realk, tmp-soeo%Etotal
!
!      !calculate predicted energy
!      val1 = soeo%dEpred
!      call soeo_get_Epred (soeo, tmpmat2, tmpvec2)
!      val2 = soeo%dEpred
!      soeo%dEpred = val1
!      write (soeo%lupri, *) '2PPP  ', i/100.0E0_realk, val2
!    enddo
!
!    call mat_free (tmpmat1)
!    call mat_free (tmpvec1)
!    call mat_free (tmpmat2)
!    call mat_free (tmpvec2)
!    call mat_free (tmpmat3)
!    call mat_free (tmpvec3)
!    call mat_free (tmpoldt)
!    call mat_free (tmpDmo)
!    call mat_free (tmpDao)
!    call mat_free (tmpFao)
!    call mat_free (tmpC)
!
!    stop 'velociraptor'
!    !end debug
!  endif
  if (resp) then
    mu = 0.0E0_realk
  else
    call soeo_binary_search (iter, soeo%space, soeo%lupri, soeo%prnt, lub_m, lub_v, lusp,&
                           & Ared, gradred, mine, soeo%settings%trust, mu)
    if (soeo%prnt>2) then
      write (soeo%lupri, '("------------------------------")')
      write (soeo%lupri, *)
    endif
  endif
  !---------------------------------------------------------------------

  !Solving the reduced space:
  !---------------------------------------------------------------------
  call soeo_solve_linear_system (iter, Ared, mu, gradred, xred, soeo%lupri)

  call mat_assign(dX_m,X_m)
  call mat_assign(dX_v,X_v)
  call soeo_solver_getX (iter, lub_m, lub_v, lusp, xred, gradp, X_m, X_v)
  call soeo_daxpy (-1.0E0_realk, X_m, X_v, dX_m, dX_v)

!  if (mine < 1.0E-5_realk .and. mine > -1.0E-5_realk) then
!
!    print *, 'Zero eigenvalue!!!'
!    write (soeo%lupri, *) 'Zero eigenvalue:'
!    write (soeo%lupri, *) '  Level-shift found:', mu
!    write (soeo%lupri, *) '  Solution vector found in the direction'
!    write (soeo%lupri, *) '    Norm(X) = ', soeo_norm (X_m, X_v)
!    write (soeo%lupri, *) '    Norm(X_m) = ', 0.5E0_realk*dsqrt(mat_sqnorm2(X_m))
!    write (soeo%lupri, *) '    Norm(X_v) = ', dsqrt(mat_sqnorm2(X_v))
!    write (soeo%lupri, *) '  Linesearch this direction for lowest predicted energy'
!
!    call soeo_linesearch_X (soeo, lub_m, lub_v, xred, iter, gradp, X_m, X_v)
!    exit
!
!  endif

  if (soeo%prnt>2) then
    write (soeo%lupri, '("Solution of linear equations:")')
    write (soeo%lupri, '("------------------------------")')
    write (soeo%lupri, '("Change in X:  ", es15.5)') soeo_dotproduct (dX_m, dX_v, dX_m, dX_v)
    write (soeo%lupri, '("Change in X_m:", es15.5)') 0.5E0_realk*mat_dotproduct (dX_m, dX_m)
    write (soeo%lupri, '("Change in X_v:", es15.5)') mat_dotproduct (dX_v, dX_v)
    write (soeo%lupri, '("------------------------------")')
    write (soeo%lupri, *)
  endif

!  if (soeo_dotproduct (dX_m, dX_v, dX_m, dX_v) < 1.0E-5_realk) then
!    write (soeo%lupri, *) 'WARNING: Change in X in this microiteration is very small', iter
!  endif

!  !debug: calculate alpha to see how large it is
!  alpha = 0.0E0_realk
!  rewind (lub_v)
!  do i=1,iter
!    call mat_read_from_disk (lub_v, bi_v,OnMaster)
!    rewind (lusp)
!    do j=1,iter
!      call mat_read_from_disk (lusp, spi,OnMaster)
!      alpha = alpha + xred(i)*xred(j)*mat_dotproduct(bi_v, spi)
!    enddo
!  enddo
!  alpha = alpha / (2.0E0_realk * mat_sqnorm2(gradp))
!
!  write (soeo%lupri, *) 'alpha in microiter:', iter, alpha
!  !end debug

!  !debug: Calculate change in occupation numbers compared to last microiteration
!  write (soeo%lupri, *) 'Microiter: Change in occs', iter
!  write (soeo%lupri, *) 'i, change, new, old'
!  if (soeo%cfg_unres) then
!    allocate (tmp1(2), tmp2(2))
!  else
!    allocate (tmp1(1), tmp2(1))
!  endif
!  do i=1,soeo%space%Nact
!    call mat_get_ab_elms (soeo%mats%oldtheta, soeo%space%Nocc+i, 1, tmp1)
!    call mat_get_ab_elms (X_v, i, 1, tmp2)
!    tmp2 = tmp1 + tmp2
!    tmp1(1) = cos(tmp1(1))*cos(tmp1(1))
!    tmp2(1) = cos(tmp2(1))*cos(tmp2(1))
!    if (soeo%cfg_unres) then
!      tmp1(2) = cos(tmp1(2))*cos(tmp1(2))
!      tmp2(2) = cos(tmp2(2))*cos(tmp2(2))
!    endif
!    write (soeo%lupri, *) i, tmp2-tmp1, tmp2, tmp1
!  enddo
!  deallocate (tmp1, tmp2)
!  !end debug

  call soeo_solver_getres (iter, lus_m, lus_v, lusp, xred, Apfac,&
                         & grad_m, grad_v, res_m, res_v)
!  call soeo_daxpy (-mu, X_m, X_v, res_m, res_v)
  rewind (lub_m) ; rewind (lub_v)
  do i=1,iter
    call mat_read_from_disk (lub_m, bi_m,OnMaster)
    call mat_read_from_disk (lub_v, bi_v,OnMaster)
    call soeo_daxpy (-mu*xred(i), bi_m, bi_v, res_m, res_v)
  enddo

!  !debug: calculate res as (PHP + Apfac*PH+P -mu*P)X + Pg
!  call mat_init (tmpmat1, soeo%space%Nbast, soeo%space%Nbast)
!  call mat_init (tmpvec1, soeo%space%Nact, 1)
!  call mat_init (tmpmat2, soeo%space%Nbast, soeo%space%Nbast)
!  call mat_init (tmpvec2, soeo%space%Nact, 1)
!  call mat_init (tmpvec3, soeo%space%Nact, 1)
!
!  tmpmat1 = X_m
!  tmpvec1 = X_v
!  call soeo_precond (soeo%space, tmpmat1, tmpvec1, orbs, soeo%cfg_unres)
!  call soeo_projectout_gradplus (soeo%cfg_grandcan, soeo%lupri, gradp, tmpvec1)
!
!  call soeo_linear_transform (soeo, tmpmat1, tmpvec1, tmpmat2, tmpvec2)
!  call soeo_precond (soeo%space, tmpmat2, tmpvec2, orbs, soeo%cfg_unres)
!  call soeo_projectout_gradplus (soeo%cfg_grandcan, soeo%lupri, gradp, tmpvec2)
!
!  call soeo_transformplus (soeo%space, soeo%cfg_unres, soeo%cfg_grandcan,&
!                         & soeo%mats%oldtheta, tmpvec1, tmpvec3)
!  call soeo_projectout_gradplus (soeo%cfg_grandcan, soeo%lupri, gradp, tmpvec3)
!
!  call mat_daxpy (Apfac, tmpvec3, tmpvec2)
!  call soeo_daxpy (-mu, tmpmat1, tmpvec1, tmpmat2, tmpvec2)
!  call soeo_daxpy (1.0E0_realk, grad_m, grad_v, tmpmat2, tmpvec2)
!
!  write (soeo%lupri, *) 'alternative resnorm = ', soeo_norm (tmpmat2, tmpvec2)
!
!  call soeo_daxpy (-1.0E0_realk, res_m, res_v, tmpmat2, tmpvec2)
!  write (soeo%lupri, *) 'difference in resnorms =', soeo_norm (tmpmat2, tmpvec2)
!
!
!  call mat_free (tmpmat1)
!  call mat_free (tmpvec1)
!  call mat_free (tmpmat2)
!  call mat_free (tmpvec2)
!  call mat_free (tmpvec3)
!  !end debug
  !---------------------------------------------------------------------

  !Test for convergence and stats:
  !---------------------------------------------------------------------
  resnorm = soeo_norm (res_m, res_v)

  if (iter == 1 .and. soeo%prnt>1) then
    write (soeo%lupri, *) "RAPTOR "
    write (soeo%lupri, *) "RAPTOR thresh =", thresh
    write (soeo%lupri, *) "RAPTOR trust =", soeo%settings%trust
    write (soeo%lupri, '(a97)') " RAPTOR iter  &
                                  & Nelec  &
                                  & x*gp      &
                                  & resnorm   &
                                  & resnormK  &
                                  & resnormt  &
                                  & xnorm     &
                                  & xnormK    &
                                  & xnormt  "
    write (soeo%lupri, '(a97)') " RAPTOR-------&
                                  &--------&
                                  &-----------&
                                  &-----------&
                                  &-----------&
                                  &-----------&
                                  &-----------&
                                  &-----------&
                                  &---------"
  endif
  tmp = soeo_numberof_electrons (soeo%space%Nocc, soeo%mats%oldtheta, X_v, soeo%cfg_unres)
  if (soeo%prnt>1) then
    write (soeo%lupri, '(" RAPTOR ", i4, f8.3, 7D11.3)')&
             & iter, tmp, mat_dotproduct(gradp, X_v), &
             & resnorm, sqrt(mat_dotproduct(res_m, res_m)), sqrt(mat_dotproduct(res_v, res_v)), &
             & soeo_norm (X_m, X_v), sqrt(mat_dotproduct(X_m, X_m)), sqrt(mat_dotproduct(X_v, X_v))
  endif

  if (resnorm < thresh) then
    write (soeo%lupri, '("Reduced space calculation converged in ", &
                   & i4, " iterations")') iter
    exit
  elseif (iter >= soeo%settings%micromaxiter) then
    write (soeo%lupri, '("Reduced space calculation NOT converged &
                   &after ", i4, " iterations!!!")') iter
    exit
!    call lsquit ("Error in soeo reduced space", soeo%lupri)
  endif
  write (soeo%lupri, *)

  iter = iter + 1
  !---------------------------------------------------------------------

  !The new trial_vectors: biter = norm(pre(norm(res)))
  !---------------------------------------------------------------------
  call mat_assign(biter_m,res_m) ; call mat_assign(biter_v,res_v)
  call soeo_cleanup_trialvector (soeo%cfg_unres, biter_m)
  call soeo_precond (soeo%space, biter_m, biter_v, orbs, soeo%cfg_unres)
  call soeo_normalize (biter_m, biter_v)
  call soeo_projectout_gradplus (soeo%cfg_grandcan, soeo%lupri, gradp, biter_v)
  call soeo_cleanup_trialvector (soeo%cfg_unres, biter_m)
  call soeo_normalize (biter_m, biter_v)
  rewind (lub_m) ; rewind (lub_v)
  do i=1,iter-1
    call mat_read_from_disk (lub_m, bi_m,OnMaster)
    call mat_read_from_disk (lub_v, bi_v,OnMaster)
    call soeo_orthonormalize (bi_m, bi_v, biter_m, biter_v, err)

    if (err) then
      write (soeo%lupri, '("Linear dependencies in soeo_solver!!!")')
      write (soeo%lupri, '("  Microiteration: ", i4)') iter
      write (soeo%lupri, '("  Norm(biter) =   ", es15.5)') soeo_norm (biter_m, biter_v)
      write (soeo%lupri, '("biter * b", i3, " =", es15.5)') i, soeo_dotproduct (biter_m, biter_v, bi_m, bi_v)
      call lsquit ("Linear dependencies in soeo_solver!!!",soeo%lupri)
    endif
  enddo

  call mat_write_to_disk (lub_m, biter_m,OnMaster)
  call mat_write_to_disk (lub_v, biter_v,OnMaster)

!  !debug: Check orthogonality between different parts of trial-vectors
!  write (soeo%lupri, *) '---------------------------------------------------'
!  write (soeo%lupri, *) 'biter_v * gradp =', mat_dotproduct (biter_v, gradp)
!  rewind (lub_m) ; rewind (lub_v)
!  do i=1,iter
!    call mat_read_from_disk (lub_m, bi_m,OnMaster)
!    call mat_read_from_disk (lub_v, bi_v,OnMaster)
!    write (soeo%lupri, *) 'biter_v * bi_v =', i, mat_dotproduct (biter_v, bi_v)
!    write (soeo%lupri, *) 'biter_m * bi_m =', i, mat_dotproduct (biter_m, bi_m)
!  enddo
!  write (soeo%lupri, *) '---------------------------------------------------'
!  !end debug

  !---------------------------------------------------------------------

  !Expanding the reduced space:
  !---------------------------------------------------------------------
  call soeo_linear_transform (soeo, biter_m, biter_v, siter_m, siter_v)
  call soeo_cleanup_trialvector (soeo%cfg_unres, siter_m)
  call soeo_precond (soeo%space, siter_m, siter_v, orbs, soeo%cfg_unres)
  call soeo_projectout_gradplus (soeo%cfg_grandcan, soeo%lupri, gradp, siter_v)

  call soeo_transformplus (soeo%space, soeo%cfg_unres, soeo%cfg_grandcan,&
                         & soeo%mats%oldtheta, biter_v, spiter)
  call soeo_projectout_gradplus (soeo%cfg_grandcan, soeo%lupri, gradp, spiter)

  call mat_write_to_disk (lus_m, siter_m,OnMaster)
  call mat_write_to_disk (lus_v, siter_v,OnMaster)
  call mat_write_to_disk (lusp, spiter,OnMaster)

  rewind (lub_m) ; rewind (lub_v)
  rewind (lus_m) ; rewind (lus_v) ; rewind (lusp)
  do i=1,iter

    call mat_read_from_disk (lub_m, bi_m,OnMaster)
    call mat_read_from_disk (lub_v, bi_v,OnMaster)

    Ared(i,iter) = soeo_dotproduct (bi_m, bi_v, siter_m, siter_v)
    Ared(i,iter) = Ared(i,iter) + Apfac * mat_dotproduct (bi_v, spiter)

    call mat_read_from_disk (lus_m, si_m,OnMaster)
    call mat_read_from_disk (lus_v, si_v,OnMaster)
    call mat_read_from_disk (lusp, spi,OnMaster)

    Ared(iter,i) = soeo_dotproduct (biter_m, biter_v, si_m, si_v)
    Ared(iter,i) = Ared(iter,i) + Apfac * mat_dotproduct (biter_v, spi)

  enddo

  !check if Ared is symmetric
  do i=1,iter
    do j=i,iter
      tmp = Ared(i,j) - Ared(j,i)
      if (tmp < -1.0E-5_realk .or. tmp > 1.0E-5_realk) then
        write (soeo%lupri, '("------------------------------")')
        write (soeo%lupri, '("Ared(", i3, ",", i3, ") - Ared(", i3, ",", i3, ") =", es15.5)') i, j, j, i, tmp
        write (soeo%lupri, '("------------------------------")')
        call lsquit ('Ared is not symmetric!', soeo%lupri)
      endif
    enddo
  enddo

  gradred(iter) = -soeo_dotproduct (biter_m, biter_v, grad_m, grad_v)
  !---------------------------------------------------------------------

  if (soeo%prnt>2) then
    write (soeo%lupri, '("==============================")')
  endif
enddo

if (soeo%prnt>2) then
  write (soeo%lupri, '("Change in total number of electrons:", es9.2)') &
          & soeo_numberof_electrons (soeo%space%Nocc, soeo%mats%oldtheta, X_v, soeo%cfg_unres) &
          & - soeo%space%Nelec
  write (soeo%lupri, '("==============================")')
  write (soeo%lupri, *)
endif

!Finally setting the output:
!-----------------------------------------------------------------------
call mat_assign(K,X_m)
call mat_assign(dt,X_v)
!-----------------------------------------------------------------------

!Finalizing
!-----------------------------------------------------------------------
call mem_dealloc (orbs)

call lsclose (lub_m, "DELETE")
call lsclose (lub_v, "DELETE")
call lsclose (lus_m, "DELETE")
call lsclose (lus_v, "DELETE")
call lsclose (lusp, "DELETE")

call mat_free (grad_m )
call mat_free (grad_v )
call mat_free (gradp  )
call mat_free (bi_m   )
call mat_free (bi_v   )
call mat_free (biter_m)
call mat_free (biter_v)
call mat_free (si_m   )
call mat_free (si_v   )
call mat_free (spi    )
call mat_free (siter_m)
call mat_free (siter_v)
call mat_free (spiter )
call mat_free (X_m    )
call mat_free (X_v    )
call mat_free (res_m  )
call mat_free (res_v  )
call mat_free (dX_m   )
call mat_free (dX_v   )
!-----------------------------------------------------------------------

end subroutine soeo_solver
!=======================================================================



!> \brief Projects out the direction where the number of electrons is changed
!> \author C. Nygaard
!> \date apr 12 2011
!> \param grandcan Should the number of electrons be preserved
!> \param lupri LU for output file
!> \param gradplus The direction to remove from the vector
!> \param vec The vector to remove the direction from
!=======================================================================
subroutine soeo_projectout_gradplus (grandcan, lupri, gradplus, vec)

implicit none

logical, intent(in)         :: grandcan
integer, intent(in)         :: lupri
type(matrix), intent(in)    :: gradplus
type(matrix), intent(inout) :: vec

integer                     :: N
real(realk)                 :: num, denom, proj

N = vec%nrow

if (.not. grandcan) then
  num = mat_dotproduct (gradplus, vec)
  denom = mat_dotproduct (gradplus, gradplus)
  if (denom > 1.0E-5_realk) then
    proj = num/denom
  else
    proj = 0.0E0_realk
  endif
  call mat_daxpy (-proj, gradplus, vec)
else
  !Nothing, the variations are allowed to change total number of electrons
  ! i.e. they are allowed to vary in the direction of gradplus
endif

end subroutine soeo_projectout_gradplus
!=======================================================================



!> \brief Solves linear system (A-mu*I)x = b
!> \author C. Nygaard
!> \date 2010
!> \param iter Size of the linear system
!> \param Ain The matrix A (in)
!> \param mu A level shift
!> \param bin The matrix b (in)
!> \param xout The matrix x (out)
!> \param lupri Logical unit number for output file (DALTON.OUT)
!=======================================================================
subroutine soeo_solve_linear_system (iter, Ain, mu, bin, xout, lupri)
implicit none

!I/O:
integer, intent(in)          :: iter
real(realk), intent(in)      :: Ain(:,:), bin(:), mu
integer, intent(in)          :: lupri
real(realk), intent(out)     :: xout(:)
!Other:
integer                      :: INFO, i
real(realk)                  :: A(iter,iter), b(iter)
integer                      :: IPIV(iter,iter)
INFO=0

xout = 0.0E0_realk
A = Ain(1:iter,1:iter)
do i=1,iter
  A(i,i) = A(i,i) - mu
enddo
b = bin(1:iter)

!See man dgesv
call dgesv (iter, 1, A, iter, IPIV, b, iter, INFO)

if (INFO < 0) then
  !INFO = -i : the i'th argument had an illegal value
  write (lupri, '(/A, i4)') "Illegal argument in dgesv, INFO =", INFO
  call lsquit(" Problem in dgesv", lupri)
elseif (INFO > 0) then
  !INFO = i : U(i,i) is exactly zero
  !           A is singular and solution cannot be computed
  write (lupri, '(/A)') "Problem in dgesv"
  write (lupri, '(/A, i4)') "A is singular, INFO =", INFO
  call lsquit(" Problem in dgesv", lupri)
else
  !INFO = 0 : successful exit of dgesv
  xout(1:iter) = b
endif

end subroutine soeo_solve_linear_system
!=======================================================================

!> \brief Finds X from the redspace solutions and trial-vectors
!> \author C. Nygaard
!> \date 2010-06-30
!> \param iter Size of the reduced space
!> \param lub_m LU for file containing matrix-part of trial vectors
!> \param lub_v LU for file containing vector-part of trial vectors
!> \param lusp LU for file containing the plus-transformed trial-vectors
!> \param xred Solution to the reduced linear system
!> \param gradp The gradient of the total number of electrons
!> \param X_m Matrix part of the solution to the linear system in the reduced space
!> \param X_v Vector part of the solution to the linear system in the reduced space
!
!> X = sum_{i}{xred(i) * bi}
!>     - 1/2*sum_{ij}{xred(i)xred(j)*(bi*sigmapj)/normsq(gradp) * gradp}
!=======================================================================
subroutine soeo_solver_getX (iter, lub_m, lub_v, lusp, xred, gradp, X_m, X_v)

implicit none

integer, intent(in)           :: iter, lub_m, lub_v, lusp
real(realk), intent(in)       :: xred(:)
type(matrix), intent(in)      :: gradp
type(matrix), intent(inout)   :: X_m, X_v

integer                       :: i, j, Nbast, Nact
type(matrix)                  :: bi_m, bi_v, spj
real(realk)                   :: sofac, gradpsqnorm
logical                       :: OnMaster
OnMaster = .TRUE.
Nbast = X_m%nrow
Nact = X_v%nrow

call mat_init (bi_m, Nbast, Nbast)
call mat_init (bi_v, Nact, 1)
call mat_init (spj, Nact, 1)

gradpsqnorm = mat_sqnorm2 (gradp)

call mat_zero (X_m)
call mat_zero (X_v)
sofac = 0.0E0_realk
rewind (lub_m) ; rewind (lub_v)
do i=1,iter
  !First order: xred(i) * bi:
  call mat_read_from_disk (lub_m, bi_m,OnMaster)
  call mat_read_from_disk (lub_v, bi_v,OnMaster)
  call mat_daxpy (xred(i), bi_m, X_m)
  call mat_daxpy (xred(i), bi_v, X_v)

  rewind (lusp)
  do j=1,iter
    !Second order: xred(i)*xred(j)*sofac * gradp
    !              sofac = -(bi*sigmapj)/(normsq(gradp))
    call mat_read_from_disk (lusp, spj,OnMaster)
    sofac = sofac + xred(i) * xred(j) * mat_dotproduct (bi_v, spj)
  enddo
enddo
if (gradpsqnorm > 1.0E-5_realk) then
  sofac = -0.5E0_realk*sofac / gradpsqnorm
  call mat_daxpy (sofac, gradp, X_v)
endif

call mat_free (bi_m)
call mat_free (bi_v)
call mat_free (spj)

end subroutine soeo_solver_getX
!=======================================================================

!> \brief Finds res from the redspace solutions and trial-vectors
!> \author C. Nygaard
!> \date 2010-06-30
!> \param iter Size of the reduced space
!> \param lus_m LU for file containing matrix-part of sigma
!> \param lus_v LU for file containing vector-part of sigma
!> \param lusp LU for file containing plus-transformed trial-vectors
!> \param xred Solution to the reduced linear system
!> \param sofac Factor in front of second-order term
!> \param grad_m Matrix part of the gradient
!> \param grad_v Vector part of the gradient
!> \param red_m Matrix part of the residual
!> \param red_v Vector part of the residual
!
!> res = sum_{i}{xred(i) * sigmai}
!>       + sum_{i}{xred(i)*sofac * sigmapi}
!>       + grad
!> sofac = -dot(grad,gradp)/sqnorm(gradp)
!>         where grad here is before gradp is projected out
!=======================================================================
subroutine soeo_solver_getres (iter, lus_m, lus_v, lusp, xred, sofac,&
                             & grad_m, grad_v, res_m, res_v)

implicit none

integer, intent(in)         :: iter, lus_m, lus_v, lusp
real(realk), intent(in)     :: xred(:), sofac
type(matrix), intent(in)    :: grad_m, grad_v
type(matrix), intent(inout) :: res_m, res_v

integer                     :: i, Nbast, Nact
type(matrix)                :: si_m, si_v, spi
logical                     :: OnMaster
OnMaster = .TRUE.
Nbast = res_m%nrow
Nact = res_v%nrow
call mat_init (si_m, Nbast, Nbast)
call mat_init (si_v, Nact, 1)
call mat_init (spi, Nact, 1)

call mat_zero (res_m)
call mat_zero (res_v)
rewind (lus_m) ; rewind (lus_v) ; rewind (lusp)

do i=1,iter
  call mat_read_from_disk (lus_m, si_m,OnMaster)
  call mat_read_from_disk (lus_v, si_v,OnMaster)
  call mat_daxpy (xred(i), si_m, res_m)
  call mat_daxpy (xred(i), si_v, res_v)

  call mat_read_from_disk (lusp, spi,OnMaster)
  call mat_daxpy (xred(i)*sofac, spi, res_v)
enddo
call mat_daxpy (1.0E0_realk, grad_m, res_m)
call mat_daxpy (1.0E0_realk, grad_v, res_v)

call mat_free (si_m)
call mat_free (si_v)
call mat_free (spi)

end subroutine soeo_solver_getres
!=======================================================================



!> \brief Finds the minimum eigenvalue of a matrix Ain
!> \author C. Nygaard
!> \date 2010-06-30
!> \param iter Dimension of Ain
!> \param Ain The matrix to which we want to find mineval
!> \param mineval The minimum eigenvalue of Ain
!=======================================================================
subroutine soeo_find_mineval (lupri, iter, Ain, mineval, minevec)

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

integer :: lueval, luevecs
logical :: evalsexist

A = Ain(1:iter,1:iter)
minevec = 0.0E0_realk

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
    do j=1,iter
      minevec(j) = X(j,i)
    enddo
  endif
enddo

!!debug
!write (lupri, *) 'Eigenvalues for A:'
!call ls_output (er, 1, iter, 1, 1, iter, 1, 1, lupri)
!write (lupri, *) 'Corresponding eigenvectors for A:'
!call ls_output (X, 1, iter, 1, iter, iter, iter, 1, lupri)

lueval = -1 ; luevecs = -1
inquire (file='evals',EXIST=evalsexist)
if (evalsexist) then
  call lsopen (lueval, "evals", "OLD", "FORMATTED")
  call lsopen (luevecs, "evecs", "OLD", "FORMATTED")
  write (lueval, *) iter
  write (luevecs, *) iter
  do i=1,iter
    write (lueval, *) er(i)
    do j=1,iter
      write (luevecs, *) X(j,i)
    enddo
  enddo
  call lsclose (lueval, "KEEP")
  call lsclose (luevecs, "KEEP")
endif
!!end debug

end subroutine soeo_find_mineval
!=======================================================================

!> \brief Binary search to find mu (the levelshift)
!> \author C. Nygaard
!> \date 2010-06-30
!> \param iter Dimension of the reduced space
!> \param space Structure containing information about the active space
!> \param lupri LU for standard output
!> \param prnt print level (if prnt > 2 search is printet to output)
!> \param lub_m LU for file containing matrix-part of trial vectors
!> \param lub_v LU for file containing vector-part of trial vectors
!> \param Ared The reduced space hessian
!> \param gradred The reduced space gradient
!> \param mine Minimum eigenvalue of Ared
!> \param h The trust-radius
!> \param mu The levelshift to control steplength
!> \param lupri LU for dalton output (DALTON.OUT)
!=======================================================================
subroutine soeo_binary_search (iter, space, lupri, prnt, lub_m, lub_v, lusp,&
                             & Ared, gradred, mine, h, mu)

implicit none

type(soeoItem_space), intent(inout):: space
integer, intent(in)          :: iter, lub_m, lub_v, lusp, lupri, prnt
real(realk), intent(in)      :: Ared(:,:), gradred(:), mine, h
real(realk), intent(out)     :: mu

integer                      :: Nbast, Nact, prntlimit
real(realk)                  :: highin, lowin
integer                      :: N, ct
real(realk), pointer         :: xred(:)
real(realk)                  :: high, middle, low, normh, normm, norml
type(matrix)                 :: X_m, X_v, zero

Nbast = space%Nbast
Nact = space%Nact

prntlimit = 2

N = size(gradred)
call mem_alloc (xred,N)

call mat_init (X_m, Nbast, Nbast)
call mat_init (X_v, Nact, 1)
call mat_init (zero, Nact, 1)
call mat_zero (zero)


!Determining where to search:
!-----------------------------------------------------------------------
if (mine >= 0.0E0_realk) then
  !Search from 0 and down
  if (prnt>prntlimit) then
    write (lupri, '("Minimum eigenvalue is larger than or equal to 0")')
    write (lupri, '(" -- search from 0 and down")')
    write (lupri, *)
  endif
  highin = 0.0E0_realk
  lowin = -16.0E0_realk
elseif (mine < 0.0E0_realk) then
  !Search from mine and down
  if (prnt>prntlimit) then
    write (lupri, '("Minimum eigenvalue is smaller than 0")')
    write (lupri, '(" -- search from mineval and down")')
    write (lupri, *)
  endif
  highin = (1.0E0_realk+1.0E-10_realk)*mine
  lowin = 10.0E0_realk*mine
endif
!-----------------------------------------------------------------------

!The search:
!-----------------------------------------------------------------------
high = highin ; low = lowin

call soeo_solve_linear_system (iter, Ared, 0.0E0_realk, gradred, xred, lupri)
call soeo_solver_getX (iter, lub_m, lub_v, lusp, xred, zero, X_m, X_v)

normh = soeo_norm (X_m, X_v)

if (prnt>prntlimit) then
  write (lupri, '("Trust-radius =", es12.5)') h
  write (lupri, '("Before the binary search (at mu=0):")')
  write (lupri, '("  Norm(X) =", es12.5)') normh
endif

call soeo_solve_linear_system (iter, Ared, high, gradred, xred, lupri)
call soeo_solver_getX (iter, lub_m, lub_v, lusp, xred, zero, X_m, X_v)

normh = soeo_norm (X_m, X_v)

if (prnt>prntlimit) then
  write (lupri, '("Before the binary search (at maximum mu):")')
  write (lupri, '("  Norm(X) =", es12.5)') normh
endif

if (normh <= h) then
  !Step will allways fall inside trust-region
  mu = high
else
  !check if wanted mu is in the given interval
  ct = 0
  do
    ct = ct + 1
    call soeo_solve_linear_system (iter, Ared, low, gradred, xred, lupri)
    call soeo_solver_getX (iter, lub_m, lub_v, lusp, xred, zero, X_m, X_v)
    norml = soeo_norm (X_m, X_v)
    if (norml > h) then
      high = low
      normh = norml
      low = low * 2.0E0_realk
      !write (lupri, *) 'mu is not in given interval - interval expanded'
    else
      exit
    endif
    if (ct > 100) then
      write (lupri, '("WARNING: mu is not in given interval!")')
      exit
    endif
  enddo

  !performing the search
  ct = 0
  do
    ct = ct + 1

    middle = (high + low)*0.5E0_realk

    call soeo_solve_linear_system (iter, Ared, middle, gradred, xred, lupri)

    call soeo_solver_getX (iter, lub_m, lub_v, lusp, xred, zero, X_m, X_v)
    normm = soeo_norm (X_m, X_v)

    if (normm - h < 1.0E-8_realk .and. normm - h > -1.0E-8_realk) then
      mu = middle
      normh = normm
      exit
    elseif (normm - h < 0.0E0_realk) then
      low = middle
      norml = normm
    elseif (normm - h > 0.0E0_realk) then
      high = middle
      normh = normm
    endif

    if (high-low < 1.0E-8_realk) then
      mu = low
      normh = norml
      exit
    endif

    if (ct > 200) then
      print '("WARNING: mu not found in 200 iterations in soeo_binary_search")'
      stop 'error when finding mu'
      exit
    endif
  enddo
endif

if (prnt>prntlimit) then
  write (lupri, '("After the binary search:")')
  write (lupri, '("  mu      =", es12.5)') mu
  write (lupri, '("  Norm(X) =", es12.5)') normh
endif

if (prnt==4) then
  !make a plot of |x| as a function of mu
  write (lupri, *)
  write (lupri, '("FINDMU")')
  write (lupri, '("FINDMU #       mu             norm    ")')
  middle = 5.0E0_realk
  do
    call soeo_solve_linear_system (iter, Ared, middle, gradred, xred, lupri)
    call soeo_solver_getX (iter, lub_m, lub_v, lusp, xred, zero, X_m, X_v)
    normm = soeo_norm (X_m, X_v)
  
    write (lupri, '("FINDMU", 2es15.5)') middle, normm
  
    middle = middle - 1.0E-1_realk
    if (middle < -20.0E0_realk) exit
  enddo
  write (lupri, *)
  write (lupri, *) 'grep FINDMU to see |X| as function of mu'
endif

call mat_free (X_m)
call mat_free (X_v)
call mat_free (zero)

call mem_dealloc (xred)

end subroutine soeo_binary_search
!=======================================================================

!> \brief Search for minimum predicted energy in direction of given X
!> \author C. Nygaard
!> \date Mar 14 2012
!> \param soeo Structure containing information and matrices
!> \param lub_m LUN for file containg matrix-part of trial-vectors
!> \param lub_v LUN for file containg vector-part of trial-vectors
!> \param xred Solution vector from reduced space
!> \param iter Number of trial-vectors in reduced space
!> \param gradp Gradient of the change of total number of electrons
!> \param X_m Matrix part of solution vector (out)
!> \param X_v Vector part of solution vector (out)
subroutine soeo_linesearch_X (soeo, lub_m, lub_v, xred, iter, gradp, X_m, X_v)

implicit none

type(soeoItem), intent(inout) :: soeo
integer, intent(in)           :: lub_m, lub_v, iter
real(realk), intent(in)       :: xred(:)
type(matrix), intent(in)      :: gradp
type(matrix), intent(inout)   :: X_m, X_v

type(matrix)                  :: mat, mato, matn !unscaled, old, new
type(matrix)                  :: vec, veco, vecn
real(realk)                   :: dEo, dEn
logical                       :: minus

type(matrix)                  :: bi_m, bi_v, sp
real(realk)                   :: delta, tmp
integer                       :: i
logical                     :: OnMaster
OnMaster = .TRUE.

call mat_init (mat, soeo%space%Nbast, soeo%space%Nbast)
call mat_init (vec, soeo%space%Nact, 1)
call mat_init (mato, soeo%space%Nbast, soeo%space%Nbast)
call mat_init (veco, soeo%space%Nact, 1)
call mat_init (matn, soeo%space%Nbast, soeo%space%Nbast)
call mat_init (vecn, soeo%space%Nact, 1)

call mat_init (bi_m, soeo%space%Nbast, soeo%space%Nbast)
call mat_init (bi_v, soeo%space%Nact, 1)
call mat_init (sp, soeo%space%Nact, 1)

!creating PX (mat,vec)
call mat_zero (mat)
call mat_zero (vec)
rewind (lub_m) ; rewind (lub_v)
do i=1,iter
  call mat_read_from_disk (lub_m, bi_m,OnMaster)
  call mat_read_from_disk (lub_v, bi_v,OnMaster)

  call soeo_daxpy (xred(i), bi_m, bi_v, mat, vec)
enddo

!zero
call mat_zero (mato)
call mat_zero (veco)
dEo = 0.0E0_realk

minus = .false.
i=0
do
  if (minus) then
    i = i-1
  else
    i = i+1
  endif
  delta = i/100.0E0_realk

  !First order
  matn = mat
  vecn = vec
  call soeo_scal (delta, matn, vecn)
  !Second order
  call soeo_transformplus (soeo%space, soeo%cfg_unres, soeo%cfg_grandcan,&
                         & soeo%mats%oldtheta, vecn, sp)
  tmp = -0.5E0_realk * mat_dotproduct(vecn, sp) / mat_sqnorm2(gradp)
  call mat_daxpy (tmp, gradp, vecn)
  !dEpred
  tmp = soeo%dEpred
  call soeo_get_Epred (soeo, matn, vecn)
  dEn = soeo%dEpred
  soeo%dEpred = tmp

  if (dEn < dEo) then
    !step gives smaller energy than last step
    mato = matn
    veco = vecn
    dEo = dEn
  else
    if (i == 1) then
      if (minus) then
        !We are already in a minimum when X=0
        exit
      else
        !plus-X-direction gives larger energy, try minus-X
        minus = .true.
        i = 0
        cycle
      endif
    else
      !last step was minimum
      exit
    endif
  endif
  if (i == 100 .or. i == -100) exit !maximum steplength reached
enddo

call mat_assign(X_m,mato)
call mat_assign(X_v,veco)

call mat_free (mat)
call mat_free (vec)
call mat_free (mato)
call mat_free (veco)
call mat_free (matn)
call mat_free (vecn)

call mat_free (bi_m)
call mat_free (bi_v)
call mat_free (sp)

end subroutine soeo_linesearch_X


!> \brief Corrects deltatheta to eliminate higher-order errors
!> \author C. Nygaard
!> \date Apr 13 2011
!> \param soeo Contains information about the calculation
!> \param alpha How much of gradplus is added back to deltatheta
!> \param deltatheta The vector to be corrected
!> \param err True if alpha could not be found, false otherwise
!>
!> Makes deltatheta --> deltatheta + alpha*gradplus
!>  so that the total number of electrons becomes correct
!=======================================================================
subroutine soeo_correct_deltatheta (soeo, alpha, deltatheta, err)

implicit none

type(soeoItem), intent(in)  :: soeo
type(matrix), intent(inout) :: deltatheta
real(realk), intent(out)    :: alpha
logical, intent(out)        :: err

integer                     :: N, i
type(matrix)                :: dt, gp
real(realk)                 :: alphap, alpham, delta
real(realk)                 :: slope0, slope1
integer                     :: limitdir, limitslope, limitN
real(realk)                 :: Nelec, Nelecp, Nelecm
character, pointer          :: orbs(:)
real(realk), pointer        :: tmp(:), tmpt(:)
real(realk)                 :: val, fact
real(realk), parameter      :: pi = 3.1415926553897932

err = .false.

if (soeo%cfg_unres) then
  call mem_alloc (orbs,2*soeo%space%Nbast)
  call mem_alloc (tmp,2)
  call mem_alloc (tmpt,2)
else
  call mem_alloc (orbs,soeo%space%Nbast)
  call mem_alloc (tmp,1)
  call mem_alloc (tmpt,1)
endif

N = soeo%space%Nact

call mat_init (dt, N, 1)
call mat_init (gp, N, 1)

call mat_assign(dt,deltatheta)

call soeo_gradplus (soeo%space, soeo%cfg_unres, soeo%mats%oldtheta, gp)
!call soeo_label_orbitals (soeo%mats%Dmo, orbs, soeo%cfg_unres, soeo%lupri, soeo%prnt)
!do i=1,N
!  if (orbs(soeo%space%Nocc+i) == 'I' .or. orbs(soeo%space%Nocc+i) == 'V') then
!    call mat_get_ab_elms (gp, soeo%space%Nocc+i, 1, tmp)
!    tmp(1) = 0.0E0_realk
!    call mat_create_ab_elms (soeo%space%Nocc+i, 1, tmp, gp)
!  endif
!  if (soeo%cfg_unres) then
!    if (orbs(soeo%space%Nbast+soeo%space%Nocc+i) == 'I' &
!      & .or. orbs(soeo%space%Nbast+soeo%space%Nocc+i) == 'V') then
!      call mat_get_ab_elms (gp, soeo%space%Nocc+i, 1, tmp)
!      tmp(2) = 0.0E0_realk
!      call mat_create_ab_elms (soeo%space%Nocc+i, 1, tmp, gp)
!    endif
!  endif
!enddo


!!debug
!write (soeo%lupri, *) ' TOOL #     alpha              Nelec'
!alpha = -1.0E0_realk
!do
!  dt = deltatheta
!  call mat_daxpy (alpha, gp, dt)
!  Nelec = soeo_numberof_electrons (soeo%space%Nocc, soeo%mats%oldtheta, &
!                                 & dt, soeo%cfg_unres)
!  write (soeo%lupri, *) 'TOOL', alpha, Nelec
!  alpha = alpha + 1.0E-2_realk
!  if (alpha > 1.0E0_realk) exit
!enddo
!!end debug

!finding limits for the binary search
!------------------------------------------------------------------------
alpha = 0.0E0_realk
delta = 1.0E-5_realk

!first find the slope of N(alpha) at alpha=0.0E0_realk
call mat_assign(dt,deltatheta)
call mat_daxpy (alpha, gp, dt)
Nelec = soeo_numberof_electrons (soeo%space%Nocc, soeo%mats%oldtheta, &
                               & dt, soeo%cfg_unres)

if (Nelec-soeo%space%Nelec < 1.0E-10_realk .and.&
            & Nelec-soeo%space%Nelec > -1.0E-10_realk) then
  !alpha = 0.0E0_realk is the right choice!
  alpha = 0.0E0_realk
else
  alphap = alpha+delta
  call mat_assign(dt,deltatheta)
  call mat_daxpy (alphap, gp, dt)
  Nelecp = soeo_numberof_electrons (soeo%space%Nocc, soeo%mats%oldtheta, &
                                  & dt, soeo%cfg_unres)
  
  alpham = alpha-delta
  call mat_assign(dt,deltatheta)
  call mat_daxpy (alpham, gp, dt)
  Nelecm = soeo_numberof_electrons (soeo%space%Nocc, soeo%mats%oldtheta, &
                                  & dt, soeo%cfg_unres)
  
  slope0 = (Nelecp - Nelecm)/(2.0E0_realk*delta)
  
  !setting restrictions for the other limit of the binary search
  if (slope0 < 0.0E0_realk) then
    if (Nelec < soeo%space%Nelec) then
      limitdir = -1    !search direction is negative
      limitslope = -1  !slope in limit should be negative
      limitN = 1       !Nelec in limit should be larger than desired Nelec
    else
      limitdir = 1
      limitslope = -1
      limitN = -1
    endif
  else
    if (Nelec < soeo%space%Nelec) then
      limitdir = 1
      limitslope = 1
      limitN = 1
    else
      limitdir = -1
      limitslope = 1
      limitN = -1
    endif
  endif
  
!  !debug
!  write (soeo%lupri, *)
!  write (soeo%lupri, *) 'DEBUGGING ALPHA'
!  write (soeo%lupri, *) 'dir, slope, N:', limitdir, limitslope, limitN
!  write (soeo%lupri, *) 'slope0 =', slope0
!  write (soeo%lupri, *) 'Nelec(0) =', Nelec
!  !end debug
  
  !then find a limit in the right direction so that the right alpha is between
  ! zero and that limit
  i = 0
  alpha = limitdir * sqrt(mat_dotproduct(deltatheta,deltatheta))
  do
    call mat_assign(dt,deltatheta)
    call mat_daxpy (alpha, gp, dt)
    Nelec = soeo_numberof_electrons (soeo%space%Nocc, soeo%mats%oldtheta, &
                                   & dt, soeo%cfg_unres)
    
    alphap = alpha+delta
    call mat_assign(dt,deltatheta)
    call mat_daxpy (alphap, gp, dt)
    Nelecp = soeo_numberof_electrons (soeo%space%Nocc, soeo%mats%oldtheta, &
                                    & dt, soeo%cfg_unres)
    
    alpham = alpha-delta
    call mat_assign(dt,deltatheta)
    call mat_daxpy (alpham, gp, dt)
    Nelecm = soeo_numberof_electrons (soeo%space%Nocc, soeo%mats%oldtheta, &
                                    & dt, soeo%cfg_unres)
    
    slope1 = (Nelecp - Nelecm)/(2.0E0_realk*delta)
 
    !debug
    write (soeo%lupri, *) 'limit =', alpha
    write (soeo%lupri, *) 'slope1 =', slope1
    write (soeo%lupri, *) 'Nelec(limit) =', Nelec
    !end debug
 
    if (limitslope*slope1 > 0.0E0_realk) then
      if (limitN*(Nelec-soeo%space%Nelec) > 0.0E0_realk) then
        !everything is as it should be
        exit
      else
        !slope has correct sign, but Nelec is on wrong side of desired Nelec
        alpha = 1.2E0_realk*alpha
      endif
    else
      if (limitN*(Nelec-soeo%space%Nelec) > 0.0E0_realk) then
        !slope has wrong sign but Nelec is on right side of desired Nelec
        alpha = 0.7*alpha
      else
        !both slope and Nelec are wrong
        alpha = 0.7*alpha
      endif
    endif
  
    !security
    i = i + 1
    if (i > 200) then
      write (soeo%lupri, *) 'WARNING: soeo_correct_deltatheta could not &
                           &find limit for alpha'
      write (soeo%lupri, *) 'rejection will be performed'
      err = .true.
      exit
!      call lsquit ('Error in soeo_correct_deltatheta: &
!                 & Limit not found', soeo%lupri)
    endif
  enddo
  
  !setting limits for binary search
  if (alpha > 0.0E0_realk) then
    alpham = 0.0E0_realk
    alphap = alpha
  else
    alpham = alpha
    alphap = 0.0E0_realk
  endif
  alpha = 0.5E0_realk*(alphap+alpham)
  !--------------------------------------------------------------------- 

  if (.not. err) then
    !doing binary search
    !---------------------------------------------------------------------
    i=0
    do
      call mat_assign(dt,deltatheta)
      call mat_daxpy (alpha, gp, dt)
      Nelec = soeo_numberof_electrons (soeo%space%Nocc, soeo%mats%oldtheta, &
                                     & dt, soeo%cfg_unres)
    
      !converged?
      if (Nelec-soeo%space%Nelec < 1.0E-10_realk .and. &
                  & Nelec-soeo%space%Nelec > -1.0E-10_realk) then
        exit
      endif
    
      if (limitslope*(Nelec-soeo%space%Nelec) < 0.0E0_realk) then
        alpham = alpha
      else
        alphap = alpha
      endif
      alpha = 0.5E0_realk*(alpham+alphap)
    
      !security
      i=i+1
      if (i>2000) then
        write (soeo%lupri, *) 'WARNING: soeo_correct_deltatheta could not find alpha'
        write (soeo%lupri, *) 'backstep will be performed'
        err = .true.
        exit
      endif
    enddo
    !---------------------------------------------------------------------
  endif
endif

call mat_daxpy (alpha, gp, deltatheta)


!!debug, see if it works
!!adjusting so that 0 <= theta(i) <= pi/2
!do i=1,soeo%space%Nact
!  call mat_get_ab_elms (deltatheta, i, 1, tmp)
!  call mat_get_ab_elms (soeo%mats%oldtheta, soeo%space%Nocc+i, 1, tmpt)
!  write(soeo%lupri, *) 'i =', i
!  write(soeo%lupri, *) 'deltatheta(i) =', tmp
!  write(soeo%lupri, *) 'oldtheta(Nocc+i) =', tmpt
!  do
!    val = tmp(1)+tmpt(1)
!    if (val/pi > 1.0E2_realk) then
!      fact=real(int(log10(val/pi)),realk)
!      fact=10.0E0_realk**fact
!    else
!        fact=real(int(val/pi)+1,realk)
!    endif
!    if (val > pi/2.0E0_realk) then
!      tmp(1) = tmp(1)-fact*pi
!    elseif (val < 0.0E0_realk) then
!      tmp(1) = -2.0E0_realk*val + tmp(1)
!    else
!      exit
!    endif
!  enddo
!  if (soeo%cfg_unres) then
!    do
!      val = tmp(2)+tmpt(2)
!      if (val/pi > 1.0E2_realk) then
!        fact=real(int(log10(val/pi)),realk)
!        fact=10.0E0_realk**fact
!      else
!        fact=real(int(val/pi)+1,realk)
!      endif
!      if (val > pi/2.0E0_realk) then
!        tmp(2) = tmp(2)-fact*pi
!      elseif (val < 0.0E0_realk) then
!        tmp(2) = -2.0E0_realk*val + tmp(2)
!      else
!        exit
!      endif
!    enddo
!  endif
!  call mat_create_ab_elms (i, 1, tmp, deltatheta)
!enddo
!!end debug

call mat_free (dt)
call mat_free (gp)

call mem_dealloc (orbs)
call mem_dealloc (tmp)
call mem_dealloc (tmpt)

end subroutine soeo_correct_deltatheta
!=======================================================================



!> \brief Labels the orbitals inactive, active and virtual
!> \author C. Nygaard
!> \date Mar. 9. 2011
!> \param Dmo The MO density matrix
!> \param orbs Contains labels 'ii', 'i', 'a', 'v' or 'vv'
!> \param unres True if calculation is unrestricted
!> \param lupri Logical unit number for DALTON.OUT
!> \param prnt Print orbital occupations if printlvl > 2
!=======================================================================
subroutine soeo_label_orbitals (Dmo, orbs, unres, lupri, prnt)

implicit none

type(matrix), intent(in) :: Dmo
character, intent(out)   :: orbs(:)
logical, intent(in)      :: unres
integer, intent(in)      :: lupri, prnt

real(realk), pointer     :: occs(:)
real(realk)              :: tmp
integer                  :: Nbast, i, prntlimit
real(realk)              :: limit1, limit2
real(realk)              :: suma, sumb, sumtot


Nbast = Dmo%nrow
limit1 = 1.0E-1_realk
limit2 = 1.0E-8_realk

prntlimit=2

!debug
if (unres) then
  call mem_alloc (occs,2)
else
  call mem_alloc (occs,1)
endif
suma = 0.0E0_realk ; sumb = 0.0E0_realk
if (prnt>prntlimit) then
  write (lupri, *)
  write (lupri, '("orbital   occupation (a and b)")')
endif
do i=1,Nbast
  call mat_get_ab_elms (Dmo, i, i, occs)
  if (prnt>prntlimit) then
    write (lupri, '(i4, f20.5)') i, occs
  endif
  suma = suma + occs(1)
  if (unres) then
    sumb = sumb + occs(2)
  endif
enddo
sumtot = suma + sumb
!write (lupri, *) 'ZZZ Total           a                b'
!write (lupri, *) 'ZZZ ', sumtot, suma, sumb
write (lupri, *)

call mem_dealloc (occs)
!end debug

do i=1,Nbast
  tmp = soeo_diag_elma (Dmo, i)
  if (tmp-1.0E0_realk < limit1 .and. tmp-1.0E0_realk > -limit1) then !Dmo(i,i) = 1
    if (tmp-1.0E0_realk < limit2 .and. tmp-1.0E0_realk > -limit2) then
      orbs(i) = 'I'
    else
      orbs(i) = 'i'
    endif
  elseif (tmp < limit1 .and. tmp > -limit1) then !Dmo(i,i) = 0
    if (tmp < limit2 .and. tmp > -limit2) then
      orbs(i) = 'V'
    else
      orbs(i) = 'v'
    endif
  elseif (tmp <= 1.0E0_realk-limit1 .and. tmp >= limit1) then !0 < Dmo(i,i) < 1
    orbs(i) = 'a'
  else !Something is wrong!
    write (lupri, *) 'alpha orbital:', i, 'occupation:', tmp
    call lsquit ('Occupation number not in interval [0;1]', lupri)
  endif
  if (unres) then
    tmp = soeo_diag_elmb (Dmo, i)
    if (tmp-1.0E0_realk < limit1 .and. tmp-1.0E0_realk > -limit1) then !Dmo(i,i) = 1
      if (tmp-1.0E0_realk < limit2 .and. tmp-1.0E0_realk > -limit2) then
        orbs(Nbast+i) = 'I'
      else
        orbs(Nbast+i) = 'i'
      endif
    elseif (tmp < limit1 .and. tmp > -limit1) then !Dmo(i,i) = 0
      if (tmp < limit2 .and. tmp > -limit2) then
        orbs(Nbast+i) = 'V'
      else
        orbs(Nbast+i) = 'v'
      endif
    elseif (tmp <= 1.0E0_realk-limit1 .and. tmp >= limit1) then !0 < Dmo(i,i) < 1
      orbs(Nbast+i) = 'a'
    else !Something is wrong!
      write (lupri, *) 'beta orbital:', i, 'occupation:', tmp
      call lsquit ('Occupation number not in interval [0;1]', lupri)
    endif
  endif
enddo

end subroutine soeo_label_orbitals
!=======================================================================



!> \brief Projects out the redundant rotations of the orbitals
!> \author C. Nygaard
!> \date Mar. 1. 2011
!> \param space Contains information about the active space
!> \param b_m The trial vector to be preconditioned, matrix part
!> \param b_v The trial vector to be preconditioned, vector part
!> \param orbs Array containing information about the orbital space
!> \param unres True if calculation is unrestricted
!=======================================================================
subroutine soeo_precond (space, b_m, b_v, orbs, unres)

implicit none

type(soeoItem_space), intent(in) :: space
type(matrix), intent(inout)      :: b_m, b_v
character, intent(in)            :: orbs(:)
logical, intent(in)              :: unres

integer                     :: Nbast, Nocc, Nact, i, j
real(realk), pointer        :: tmp(:)

if (unres) then
  call mem_alloc (tmp,2)
else
  call mem_alloc (tmp,1)
endif
Nbast = space%Nbast
Nocc = space%Nocc
Nact = space%Nact

!!first the vector-part
!do i=1,Nact
!  if (orbs(Nocc+i) == 'I') then
!    call mat_get_ab_elms (b_v, i, 1, tmp)
!    tmp(1) = 0.0E0_realk
!    call mat_create_ab_elms (i, 1, tmp, b_v)
!  elseif (orbs(Nocc+i) == 'V') then
!    call mat_get_ab_elms (b_v, i, 1, tmp)
!    tmp(1) = 0.0E0_realk
!    call mat_create_ab_elms (i, 1, tmp, b_v)
!  endif
!
!  if (unres) then
!    if (orbs(Nbast+Nocc+i) == 'I') then
!    call mat_get_ab_elms (b_v, i, 1, tmp)
!    tmp(2) = 0.0E0_realk
!    call mat_create_ab_elms (i, 1, tmp, b_v)
!  elseif (orbs(Nbast+Nocc+i) == 'V') then
!    call mat_get_ab_elms (b_v, i, 1, tmp)
!    tmp(2) = 0.0E0_realk
!    call mat_create_ab_elms (i, 1, tmp, b_v)
!  endif
!  endif
!enddo

!then the matrix-part
do i=1,Nbast
  if (orbs(i) == 'I') then
    !inactive orbitals do not rotate with other inactive orbitals
    do j=1,Nbast
      if (orbs(j) == 'I') then
        call mat_get_ab_elms (b_m, i, j, tmp)
        tmp(1) = 0.0E0_realk
        call mat_create_ab_elms (i, j, tmp, b_m)
        call mat_get_ab_elms (b_m, j, i, tmp)
        tmp(1) = 0.0E0_realk
        call mat_create_ab_elms (j, i, tmp, b_m)
      endif
    enddo
  elseif (orbs(i) == 'V') then
    !virtual orbitals do not rotate with other virtual orbitals
    do j=1,Nbast
      if (orbs(j) == 'V') then
        call mat_get_ab_elms (b_m, i, j, tmp)
        tmp(1) = 0.0E0_realk
        call mat_create_ab_elms (i, j, tmp, b_m)
        call mat_get_ab_elms (b_m, j, i, tmp)
        tmp(1) = 0.0E0_realk
        call mat_create_ab_elms (j, i, tmp, b_m)
      endif
    enddo
  elseif (orbs(i) == 'a' .or. orbs(i) == 'i' .or. orbs(i) == 'v') then
    !nothing
!    !try instead to not rotate active with other active
!    do j=1,Nbast
!      if (orbs(j) == 'a' .or. orbs(i) == 'i' .or. orbs(i) == 'v') then
!        call mat_get_ab_elms (b_m, i, j, tmp)
!        tmp(1) = 0.0E0_realk
!        call mat_create_ab_elms (i, j, tmp, b_m)
!        call mat_get_ab_elms (b_m, j, i, tmp)
!        tmp(1) = 0.0E0_realk
!        call mat_create_ab_elms (j, i, tmp, b_m)
!      endif
!    enddo
  else
    call lsquit ('Something is wrong with the orbital categorization', -1)
  endif

  if (unres) then !and for the beta-part
    if (orbs(Nbast+i) == 'i' .or. orbs(Nbast+i) == 'I') then
      !inactive orbitals do not rotate with other inactive orbitals
      do j=1,Nbast
        if (orbs(Nbast+j) == 'i' .or. orbs(Nbast+j) == 'I') then
          call mat_get_ab_elms (b_m, i, j, tmp)
          tmp(2) = 0.0E0_realk
          call mat_create_ab_elms (i, j, tmp, b_m)
          call mat_get_ab_elms (b_m, j, i, tmp)
          tmp(2) = 0.0E0_realk
          call mat_create_ab_elms (j, i, tmp, b_m)
        endif
      enddo
    elseif (orbs(Nbast+i) == 'v' .or. orbs(Nbast+i) == 'V') then
      !virtual orbitals do not rotate with other virtual orbitals
      do j=1,Nbast
        if (orbs(Nbast+j) == 'v' .or. orbs(Nbast+j) == 'V') then
          call mat_get_ab_elms (b_m, i, j, tmp)
          tmp(2) = 0.0E0_realk
          call mat_create_ab_elms (i, j, tmp, b_m)
          call mat_get_ab_elms (b_m, j, i, tmp)
          tmp(2) = 0.0E0_realk
          call mat_create_ab_elms (j, i, tmp, b_m)
        endif
      enddo
    elseif (orbs(Nbast+i) == 'a') then
      !nothing
!      !try instead to not rotate active with other active
!      do j=1,Nbast
!        if (orbs(Nbast+j) == 'a') then
!          call mat_get_ab_elms (b_m, i, j, tmp)
!          tmp(2) = 0.0E0_realk
!          call mat_create_ab_elms (i, j, tmp, b_m)
!          call mat_get_ab_elms (b_m, j, i, tmp)
!          tmp(2) = 0.0E0_realk
!          call mat_create_ab_elms (j, i, tmp, b_m)
!        endif
!      enddo
    else
      call lsquit ('Something is wrong with the orbital categorization', -1)
    endif
  endif
enddo

call mem_dealloc (tmp)

end subroutine soeo_precond
!=======================================================================


!> \brief Routine that makes sure that the trial-vector is antisymmetric
!> \author C. Nygaard
!> \date Jul 25 2011
!> \param unres True if calculation is unrestricted
!> \param b The trial vector to be cleaned
!=======================================================================
subroutine soeo_cleanup_trialvector (unres, b)

implicit none

type(matrix), intent(inout) :: b
logical, intent(in)         :: unres

integer                     :: i, j, Nbast
real(realk), pointer        :: bij(:), bji(:), avg(:)

Nbast = b%nrow

if (unres) then
  call mem_alloc (bij,2)
  call mem_alloc (bji,2)
  call mem_alloc (avg,2)  
else
  call mem_alloc (bij,1)
  call mem_alloc (bji,1)
  call mem_alloc (avg,1)
endif

do i=1,Nbast
  do j=i,Nbast
    call mat_get_ab_elms (b, i, j, bij)
    call mat_get_ab_elms (b, j, i, bji)

    !if the difference is too large: make an error message
    if (bij(1)+bji(1) < -1.0E-5_realk .or. bij(1)+bji(1) > 1.0E-5_realk) then
      print *, 'i, j =', i, j
      print *, 'bij =', bij(1)
      print *, 'bji =', bji(1)
      print *, 'diff =', bij(1) + bji(1)
      call lsquit ('a: Trial vector is too asymmetric, something is wrong', -1)
    endif
    if (unres) then
      if (bij(2)+bji(2) < -1.0E-5_realk .or. bij(2)+bji(2) > 1.0E-5_realk) then
        print *, 'i, j =', i, j
        print *, 'bij =', bij(2)
        print *, 'bji =', bji(2)
        print *, 'diff =', bij(2) + bji(2)
        call lsquit ('b: Trial vector is too asymmetric, something is wrong', -1)
      endif
    endif

    avg = 0.5E0_realk*(bij-bji)

    call mat_create_ab_elms (i, j, avg, b)
    call mat_create_ab_elms (j, i, -1.0E0_realk*avg, b)
  enddo
enddo

call mem_dealloc (bij)
call mem_dealloc (bji)
call mem_dealloc (avg)

end subroutine soeo_cleanup_trialvector
!=======================================================================

end module soeo_redspace
