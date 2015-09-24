!> @file
!> Contains the soeo_loop module
!
!> brief The main loop of SOEO
!> author C. Nygaard
!> date 2010
module soeo_loop
use files

use precision
use soeo_typedef
use Matrix_module
use soeo_util
use soeo_redspace
use matrix_util
use Fock_evaluator, only: FCK_get_fock
use Matrix_Operations
use Matrix_Operations_aux
use soeo_debug
use extra_output
use dal_interface, only: di_get_dipole
private
public :: soeoloop, soeo_restart

!Contains:
!  soeoloop
!  soeo_init
!  soeo_setold
!  soeo_stepforth
!  soeo_stepback
!  soeo_restart


contains

!> \brief Makes Second Order Ensemble Optimization
!> \author C. Nygaard
!> \date 2010
!> \param Fao The Fock matrix in AO basis
!> \param S The overlap matrix
!> \param Dao The density matrix in the AO basis
!> \param soeoinp Information from input
!
!> Ordinary RH-optimization (arh or similar) should be done first
!>  FC = SCe
!> such that a coefficient matrix C is available
!> together with an AO-density-matrix and an overlap matrix
!=======================================================================
subroutine soeoloop (Fao, S, Dao, soeoinp)

implicit none

!Needed matrices
type(soeoItem_input), intent(in) :: soeoinp
type(soeoItem)                   :: soeo
type(matrix), intent(in)         :: S, Dao, Fao
type(matrix)                     :: K, deltatheta
type(matrix)                     :: Dmo_diff
type(matrix)                     :: grad_m, grad_v, gradp
type(matrix)                     :: V(3)
real(realk), pointer             :: occs(:), orbE(:)
real(realk)                      :: ratio, Dmo_diffnorm, gradnorm
real(realk)                      :: thresh, trust, mu, lambda, alpha !
real(realk)                      :: dip(3), pol(3,3), hyp(3,3,3)
integer                          :: maxiter       ! Those are for the red-space
integer                          :: Nelec, Nocc, Nact, Nvirt
integer                          :: i, j, m, l, u, iter, ncycles
logical                          :: stable, loopdone, tmplog, aufbau, grandcan
type(matrix)                     :: tmpmat, tmpvec, tmps, tmpt, tmpdmo, tmpc, tmpdao, tmpfao

!Other
logical                          :: converged, err
real(realk)                      :: tstart, tloop1, tloop2, tend, dt
real(realk)                      :: tmp, val1, val2, x
real(realk), pointer             :: tmp1(:), tmp2(:)
real(realk), pointer             :: hessian(:,:), gradient(:), step(:), smallhessian(:,:)
real(realk), pointer             :: hessian2(:,:), gradient2(:)
integer                          :: tmpint, tmpdim
real(realk)                      :: suma, sumb, sumtot
integer                          :: lusave, luint, lueval, luevecs
logical                          :: integrals_exist,OnMaster
character, pointer               :: orbs(:)
logical                          :: debug
OnMaster=.TRUE.
!print *, 'soeo%lupri =', soeo%lupri
!print *, 'soeoinp%lupri =', soeoinp%lupri

        write (soeoinp%lupri, *) "SOEOLOOP STARTED!"
        write (soeoinp%lupri, *) "============================================="

call cpu_time (tstart)

!Initializations:
!-----------------------------------------------------------------------
Nelec = nint(mat_TrAB (Dao, S))
!sanity check sum(given orbital occupations) = Nelec ?
if (.not. soeoinp%cfg_grandcan) then
  if (soeoinp%occsinput) then
    tmp = soeoinp%Nfullocc + sum(soeoinp%fracoccs(1:soeoinp%Nfracocc))
    if (soeoinp%cfg_unres) then
      tmp = tmp + soeoinp%Nfulloccb + sum(soeoinp%fracoccsb(1:soeoinp%Nfracoccb))
    endif
  else
    tmp = soeoinp%Nocc
    if (soeoinp%cfg_unres) then
      tmp = tmp * 2
    endif
  endif
  if (tmp-Nelec < -1.0E-5_realk .or. tmp-Nelec > 1.0E-5_realk) then
    print '("Nelec     =", i10)', Nelec
    print '("Tr(SD)    =", f10.2)', mat_TrAB(Dao,S)
    print '("Nfullocc  =", i10)', soeoinp%Nfullocc
    print '("Nfracocc  =", i10)', soeoinp%Nfracocc
    do i=1,soeoinp%Nfracocc,5
      print '("fracoccs =", 5f8.5)', soeoinp%fracoccs(i:soeoinp%Nfracocc)
    enddo
    print '("Nfulloccb =", i10)', soeoinp%Nfulloccb
    print '("Nfracoccb =", i10)', soeoinp%Nfracoccb
    do i=1,soeoinp%Nfracoccb,5
      print '("fracoccsb =", 5f8.5)', soeoinp%fracoccsb(i:soeoinp%Nfracoccb)
    enddo
    print '("sum(fracoccs) = ", f15.5)', sum(soeoinp%fracoccs(1:soeoinp%Nfracocc))
    print '("sum(fracoccsb) =", f15.5)', sum(soeoinp%fracoccsb(1:soeoinp%Nfracoccb))
    print '("Your given number of electrons =", f15.5)', tmp
    call lsquit ('wrong number of electrons in LSDALTON.INP', soeoinp%lupri)
  endif
endif

if (soeoinp%cfg_unres) then
  call mem_alloc (orbE,2*S%nrow)
  call mem_alloc (occs,2*S%nrow)
  call mem_alloc (orbs,2*S%nrow)
else
  call mem_alloc (orbE,S%nrow)
  call mem_alloc (occs,S%nrow)
  call mem_alloc (orbs,S%nrow)
endif
call mat_init (K        , S%nrow, S%nrow)
call mat_init (Dmo_diff , S%nrow, S%nrow)
call mat_init (grad_m   , S%nrow, S%nrow)

call soeoItem_init (S%nrow, Nelec, soeo, soeoinp)
call soeo_init (Dao, Fao, S, soeo)

call mat_init (deltatheta, soeo%space%Nact, 1)
call mat_init (grad_v, soeo%space%Nact, 1)
call mat_init (gradp, soeo%space%Nact, 1)
aufbau = .false.

!DEBUGGING PRINTS
!
!write (soeo%lupri, *) 'Old Dao ='
!call mat_print (Dao, 1, soeo%Nbast, 1, soeo%Nbast, soeo%lupri)
!write (soeo%lupri, *) 'Old Fao ='
!call mat_print (Fao, 1, soeo%Nbast, 1, soeo%Nbast, soeo%lupri)
!write (soeo%lupri, *) 'Old S ='
!call mat_print (S, 1, soeo%Nbast, 1, soeo%Nbast, soeo%lupri)
!write (soeo%lupri, *) 'My Dao ='
!call mat_print (soeo%Dao, 1, soeo%Nbast, 1, soeo%Nbast, soeo%lupri)
!write (soeo%lupri, *) 'My Fao ='
!call mat_print (soeo%Fao, 1, soeo%Nbast, 1, soeo%Nbast, soeo%lupri)
!write (soeo%lupri, *) 'My Dmo ='
!call mat_print (soeo%Dmo, 1, soeo%Nbast, 1, soeo%Nbast, soeo%lupri)
!write (soeo%lupri, *) 'My Fmo ='
!call mat_print (soeo%Fmo, 1, soeo%Nbast, 1, soeo%Nbast, soeo%lupri)
!write (soeo%lupri, *) 'My C ='
!call mat_print (soeo%C, 1, soeo%Nbast, 1, soeo%Nbast, soeo%lupri)
!
!END DEBUGGING PRINTS
!-----------------------------------------------------------------------


!Printing active space:
!-----------------------------------------------------------------------
write (soeo%lupri, *)
if (soeo%cfg_unres) then
  write (soeo%lupri, '("Number of electrons:      ", i7)') Nelec
else
  write (soeo%lupri, '("Number of electron pairs: ", i7)') Nelec
endif
write (soeo%lupri, '("Number of basis-functions:", i7)') soeo%space%Nbast
write (soeo%lupri, '("Size of occupied space:   ", i7)') soeo%space%Nocc
write (soeo%lupri, '("Size of active space:     ", i7)') soeo%space%Nact
write (soeo%lupri, '("Size of virtual space:    ", i7)') &
                       & soeo%space%Nbast - soeo%space%Nocc - soeo%space%Nact
write (soeo%lupri, *)
!-----------------------------------------------------------------------

!Reading input for starting occupations
!------------------------------------------------------------------------
if (soeoinp%occsinput) then
  call soeo_setoccs (soeoinp, soeo%mats%Dmo, soeo%mats%oldtheta)
  aufbau = .false.
else
  write (soeo%lupri, '("No starting occupations given, no permutations")')
  write (soeo%lupri, '(" - NO SOEO CALCULATION!")')
  aufbau = .true.
endif
!------------------------------------------------------------------------

call mat_zero (K) ; call mat_zero (deltatheta)
loopdone = .false.
iter = 0
ncycles = 0
mu = 0.0E0_realk
lambda = 0.0E0_realk
alpha = 0.0E0_realk

if (.not. aufbau) then
  !Start the soeo-loop
  !---------------------------------------------------------------------
  loopdone = .true.
  do

    call cpu_time (tloop1)
    iter = iter + 1
    soeo%iter = iter
    write (soeo%lupri, '("Macroiteration ", i4)') iter
    write (soeo%lupri, '("******************************")')

!write (soeo%lupri, *) 'Nelec =', soeo_numberof_electrons (soeo%space%Nocc, soeo%mats%oldtheta,&
!                                                        & deltatheta, soeo%cfg_unres)

    !Get the new matrices:
    !-------------------------------------------------------------------
    call soeo_stepforth (K, deltatheta, soeo)
    !-------------------------------------------------------------------

!!DEBUGGING
!write (soeo%lupri, *) 'dE and dEpred', soeo%dE, soeo%dEpred
!write (soeo%lupri, *) "norm(K,deltatheta), trust =", soeo_norm (K, deltatheta), soeo%settings%trust
!!END DEBUGGING

!    !debug: (??) Set up tt-hessian and diagonalize it to optimize the 
!    !       occupation numbers in the orbitals found in last iteration
!    !-------------------------------------------------------------------
!    call soeo_get_full_tt_hessian (soeo)
!    !-------------------------------------------------------------------


    !Test for convergence:
    !-------------------------------------------------------------------
    call soeo_label_orbitals (soeo%mats%Dmo, orbs, soeo%cfg_unres, soeo%lupri, 0)
    call soeo_get_gradient (soeo, grad_m, grad_v)
    call soeo_precond (soeo%space, grad_m, grad_v, orbs, soeo%cfg_unres)
    call soeo_gradplus (soeo%space, soeo%cfg_unres, soeo%mats%oldtheta, gradp)
!write (soeo%lupri, *) 'grad_v ='
!call mat_print (grad_v, 1, soeo%space%Nact, 1, 1, soeo%lupri)
!write (soeo%lupri, *) 'gradp and lambda =', lambda
!call mat_print (gradp, 1, soeo%space%Nact, 1, 1, soeo%lupri)
    call soeo_projectout_gradplus (soeo%cfg_grandcan, soeo%lupri, gradp, grad_v)
!    call mat_daxpy (alpha, gradp, grad_v)
!write (soeo%lupri, *) 'new grad_v ='
!call mat_print (grad_v, 1, soeo%space%Nact, 1, 1, soeo%lupri)
    gradnorm = soeo_norm (grad_m, grad_v)
!    write (soeo%lupri, *) "gradnorm =", gradnorm
!    write (soeo%lupri, *) 'grad-sqnorm K =', 0.5E0_realk * mat_sqnorm2(grad_m)
!    write (soeo%lupri, *) 'grad-sqnorm t =', mat_sqnorm2(grad_v)
!    write (soeo%lupri, *) 'grad-sqnorm =',&
!                         & 0.5E0_realk * mat_sqnorm2(grad_m) + mat_sqnorm2(grad_v)

    if (soeo%test .and. iter==1) then
      write (soeo%lupri, *) 'TEST: First gradnorm =', gradnorm
    endif

!    call mat_add (1.0E0_realk, soeo%old_Dmo, -1.0E0_realk, soeo%Dmo, Dmo_diff)
!    Dmo_diffnorm = mat_sqnorm2 (Dmo_diff)
!    Dmo_diffnorm = sqrt(Dmo_diffnorm)
!    write (soeo%lupri, *) "Dmo_diffnorm =", Dmo_diffnorm
    !Stats
    !-------------------------------------------------------------------
    if (iter == 1) then
      write (soeo%lupri, '(a104)') " OOO  iter    Etotal        &
                                          & dEpred        &
                                          & dE         &
                                          & ratio  &
                                          & mu    &
                                          & alpha   &
                                          & xnorm    &
                                          & gradnorm      "
      write (soeo%lupri, '(a104)') " OOO------------------------&
                                          &---------------&
                                          &---------------&
                                          &-------&
                                          &-------&
                                          &-------&
                                          &----------&
                                          &---------------"
      print ('(a101)'), "   iter    Etotal        &
                                          & dEpred        &
                                          & dE         &
                                          & ratio  &
                                          & mu    &
                                          & alpha   &
                                          & xnorm    &
                                          & gradnorm      "
      print ('(a101)'), " ------------------------&
                                          &---------------&
                                          &---------------&
                                          &-------&
                                          &-------&
                                          &-------&
                                          &----------&
                                          &---------------"
    endif
!    if (soeo%dE > 1.0E-7_realk .or. soeo%dE < -1.0E-7_realk) then
    if (soeo%dE /= 0.0E0_realk) then
      ratio = soeo%dEpred/soeo%dE
      if (ratio > 1.0E0_realk .and. soeo%dEpred /= 0.0E0_realk) then
        ratio = 1.0E0_realk / ratio
      endif
    else
      ! ratio is not changed from last iteration (because last iteration was rejected)
    endif
    write (soeo%lupri, '(" OOO ", i4, f15.8, 2D15.4, 3f7.2, d10.2, d15.4)')&
           & iter, soeo%Etotal, soeo%dEpred, soeo%dE, ratio, mu, alpha,&
           & soeo_norm(K,deltatheta), gradnorm
    print ('("  ", i4, f15.8, 2D15.4, 3f7.2, d10.2, d15.4)'),&
           & iter, soeo%Etotal, soeo%dEpred, soeo%dE, ratio, mu, alpha,&
           & soeo_norm(K,deltatheta), gradnorm
    !-------------------------------------------------------------------


!    write (soeo%lupri, *) "Etotal and dE:", soeo%Etotal, soeo%dE

!!debug
!!call soeo_isit_aufbau (soeo, aufbau)
!    call mat_init (tmpmat, soeo%space%Nbast, soeo%space%Nbast)
!    call mat_init (tmpvec, soeo%space%Nact, 1)
!    tmpmat = grad_m ; call mat_scal (-1.0E-8_realk,tmpmat)
!    tmpvec = grad_v ; call mat_scal (-1.0E-8_realk,tmpvec)
!    val1 = soeo%dE ; val2 = soeo%dEpred
!    call soeo_stepforth (tmpmat, tmpvec, soeo)
!    write (soeo%lupri, *) 'dE when step is -1.0E-8_realk*grad =', iter, soeo%dE
!    call soeo_stepback (soeo)
!    soeo%dE = val1 ; soeo%dEpred = val2
!    call mat_free (tmpmat)
!    call mat_free (tmpvec)
!!end debug

    converged = .false.

    !convergence criteria: gradnorm < thresh
    if (gradnorm < soeo%settings%macrothresh .and. iter > 1 .and. soeo%dE < 0.0E0_realk) then
      converged = .true.
    endif

    !convergence criteria: |dE/E| < thresh
    if ((soeo%dE/soeo%Etotal > -soeo%settings%macrothresh**2 .and.&
        & soeo%dE/soeo%Etotal < soeo%settings%macrothresh**2) .and. soeo%dE < 0.0E0_realk) then
      converged = .true.
    endif

    if (converged) then
      if (soeo%space%Nact /= 0) then
        call soeo_isit_aufbau2 (soeo%mats%Dao, soeo%mats%Fao, soeo%mats%S, soeo%cfg_unres, soeo%lupri, aufbau)
      else
        aufbau = .true.
      endif
      if (aufbau) then
        exit
      else
        !perturb
        exit !exit anyway
!        write (soeo%lupri, *) 'Transistion state found, occupation numbers perturbed'
!        print *, 'Transistion state found, occupation numbers perturbed'
!        write (soeo%lupri, *) 'Results for transistion state:'
!        write (soeo%lupri, *) '========================================='
!        call soeo_output_results (soeo)
!        write (soeo%lupri, *) '========================================='
!        write (soeo%lupri, *)
        iter = 0
!stop 'velociraptor'
        ncycles = ncycles + 1
        soeo%settings%trust = soeoinp%trust
        if (ncycles > 2) then
          call lsquit ('soeo keeps finding non-aufbau solutions', soeo%lupri)
        endif
        call mat_init (tmpmat, soeo%space%Nbast, soeo%space%Nbast)
        call mat_init (tmpvec, soeo%space%Nact, 1)
        call mat_zero (tmpmat) ; call mat_zero (tmpvec)
        call soeo_stepforth (tmpmat, tmpvec, soeo)
        call soeo_setold (soeo%mats, soeo%old)
        call mat_free (tmpmat)
        call mat_free (tmpvec)
      endif
    endif

    if (iter > soeo%settings%macromaxiter) then
      print *, "Maximum soeo-iterations reached - sorry"
      exit
    endif
    !-------------------------------------------------------------------

    !Update the trust-radius
    !-------------------------------------------------------------------
    if (iter > 1 .and. (.not. err)) then
      if (soeo%dE > 0.0E0_realk) then
        !rejection
        print *, 'Rejection: dE > 0.0E0_realk'
        write (soeo%lupri, *) ' OOO Rejection: dE > 0.0E0_realk'
        soeo%settings%trust = min(0.7E0_realk*soeo_norm(K, deltatheta),0.7E0_realk*soeo%settings%trust)

        call soeo_stepback (soeo)


      elseif (ratio > 0.75E0_realk .and. mu < -1.0E-8_realk) then
        !increase trust radius
        soeo%settings%trust = 1.2E0_realk * soeo%settings%trust
      elseif (ratio < 0.25E0_realk) then
        !decrease trust radius
        if (soeo%dE<-1.0E-7_realk) then
          soeo%settings%trust = 0.7E0_realk * soeo%settings%trust
        else
          !nothing, energy change is too small to determine ratio in a proper way
        endif
        if (ratio < 0.0E0_realk) then
          !rejection
          if (soeo%dE > 0.0E0_realk) then
            print *, 'Rejection: dE > 0.0E0_realk'
            write (soeo%lupri, *) ' OOO Rejection: dE > 0.0E0_realk'
            call soeo_stepback (soeo)
          elseif (soeo%dEpred > 0.0E0_realk) then
            print *, 'WARNING: dEpred > 0.0E0_realk but dE < 0.0E0_realk'
            write (soeo%lupri, *) ' OOO WARNING: dEpred > 0.0E0_realk but dE < 0.0E0_realk'
          endif
        endif
      endif
    endif
    if (soeo%settings%trust < soeo%settings%macrothresh) then
      !debug
      print *, "WARNING: Trust-radius too small, a step in the gradient direction will be taken"
      write (soeo%lupri, *) "WARNING: Trust-radius too small, a step in the gradient-direction will be taken"
      
      call mat_assign(K,grad_m)
      call mat_assign(deltatheta,grad_v)
      
      call mat_scal (-1.0E-8_realk, K)
      call mat_scal (-1.0E-8_realk, deltatheta)
      
      call soeo_stepforth (K, deltatheta, soeo)
      soeo%settings%trust = 0.1E0_realk
      
      print *, 'A step in the -1.0e-8_realk grad direction gives dE =', soeo%dE
      write (soeo%lupri, *) 'A step in the -1.0e-8_realk grad direction gives dE =', soeo%dE
      
      !exit
      !end debug
      call lsquit ('trust radius is now smaller than macrothresh', soeo%lupri)
    endif
!    if (soeo%settings%trust > 0.9E0_realk) then
!      soeo%settings%trust = 0.9E0_realk
!    endif
    !-------------------------------------------------------------------

    !Solve the soeo-equations in reduced space
    !-------------------------------------------------------------------
    call mat_zero (K) ; call mat_zero (deltatheta)

    call soeo_solver (soeo, K, deltatheta, mu, grad_m, grad_v, .false.)

    !in order to get the right number of electrons, we must eliminate the
    ! second and higher order errors in deltatheta
    if (.not. soeo%cfg_grandcan) then
      call soeo_correct_deltatheta (soeo, alpha, deltatheta, err)
      if (err) then
        call mat_zero (K)
        call mat_zero (deltatheta)
        soeo%settings%trust = 0.7E0_realk*soeo%settings%trust
        alpha = 0.0E0_realk
        print *, 'Rejection: Could not find alpha'
        write (soeo%lupri, *) ' OOO Rejection: Could not find alpha'
      endif
!      write (soeo%lupri, *) 'Correction factor, alpha =', iter+1, alpha, &
!              & soeo_numberof_electrons(soeo%space%Nocc, soeo%mats%oldtheta, deltatheta, soeo%cfg_unres)
    endif

!    !debug
!    write (soeo%lupri, *) 'Change in occupation numbers and actual occupation numbers:'
!    write (soeo%lupri, *) 'i, change, new, old'
!    if (soeo%cfg_unres) then
!      call mem_alloc (tmp1, 2)
!      call mem_alloc (tmp2, 2)
!    else
!      call mem_alloc (tmp1, 1)
!      call mem_alloc (tmp2, 1)
!    endif
!    do i=1,soeo%space%Nact
!      call mat_get_ab_elms (soeo%mats%oldtheta, soeo%space%Nocc+i, 1, tmp1)
!      call mat_get_ab_elms (deltatheta, i, 1, tmp2)
!      tmp2 = tmp1 + tmp2
!      tmp1(1) = cos(tmp1(1))*cos(tmp1(1))
!      tmp2(1) = cos(tmp2(1))*cos(tmp2(1))
!      if (soeo%cfg_unres) then
!        tmp1(2) = cos(tmp1(2))*cos(tmp1(2))
!        tmp2(2) = cos(tmp2(2))*cos(tmp2(2))
!      endif
!      write (soeo%lupri, *) i, tmp2-tmp1, tmp2, tmp1
!    enddo
!    call mem_dealloc (tmp1)
!    call mem_dealloc (tmp2)
!    !end debug
    !-------------------------------------------------------------------

    !Finding the predicted energy change
    !-------------------------------------------------------------------
    call soeo_get_Epred (soeo, K, deltatheta)
!    if (soeo%dEpred > 0.0E0_realk) then
!      write (soeo%lupri, *) " OOO Step ruined by alpha, change in occupations set to zero"
!      print *, "Step ruined by alpha, change in occupations set to zero"
!      call mat_zero (deltatheta)
!      call soeo_get_Epred (soeo, K, deltatheta)
!      write (soeo%lupri, *) " OOO New dEpred =", soeo%dEpred
!      print *, "New dEpred =", soeo%dEpred
!    endif
    !-------------------------------------------------------------------

    !Maybe some averaging here?

    !Energy?
    !-------------------------------------------------------------------
    call cpu_time (tloop2)
    dt = tloop2 - tloop1

    write (soeo%lupri, '("Iteration ", i3, " done in ", f15.10, " seconds")') &
                  & iter, dt
!    print '("Iteration ", i3, " done in ", f15.10, " seconds")', iter, dt

    write (soeo%lupri, '("******************************")')
    write (soeo%lupri, *)
  enddo
  !---------------------------------------------------------------------
else
  call soeo_stepforth (K, deltatheta, soeo)
  write (soeo%lupri, *) "Your system is stable with the previously calculated&
                   & occupation numbers - "
  write (soeo%lupri, *) "No need for SOEO-calculation."
endif

write (soeo%lupri, *) 'Final results from SOEO:'
write (soeo%lupri, *) '================================================='
call soeo_output_results (soeo)
write (soeo%lupri, *) '================================================='

!call soeo_compare_hessians (soeo)

if (soeoinp%cfg_dipole) then

  do i=1,3
    call mat_init (V(i), soeo%space%Nbast, soeo%space%Nbast)
    call mat_zero (V(i))
  enddo
  call di_get_dipole (V)

  call soeo_dipole (V, soeo%mats%Dao, soeo%cfg_unres, soeo%lupri)

  call soeo_polarizability (soeo, V)

  do i=1,3
    call mat_free (V(i))
  enddo

endif

if (soeoinp%cfg_save) then
  !print Fao and Dao in file soeosave.out
  lusave = -1
  call lsopen (lusave, "soeosave.out", "unknown", "UNFORMATTED")
  rewind (lusave)
  call mat_write_to_disk (lusave, soeo%mats%Fao,OnMaster)
  call mat_write_to_disk (lusave, soeo%mats%Dmo,OnMaster)
  call lsclose (lusave, "KEEP")
endif

!Deallocations:
!-----------------------------------------------------------------------
call mat_free (K)
call mat_free (Dmo_diff)
call mat_free (grad_m)

call mem_dealloc (orbE)
call mem_dealloc (occs)
call mem_dealloc (orbs)

call mat_free (deltatheta)
call mat_free (grad_v)
call mat_free (gradp)

call soeoItem_free (soeo)

inquire (file='ERI.Integrals',EXIST=integrals_exist)
if (integrals_exist) then
  luint = -1
  call lsopen (luint,'ERI.Integrals','UNKNOWN','UNFORMATTED')
  call lsclose (luint, 'DELETE')
endif

!-----------------------------------------------------------------------

        print *, "SOEOLOOP FINISHED!"
        write (soeo%lupri, *) "SOEOLOOP FINISHED!"
        write (soeo%lupri, *) "Number of soeo-iterations =", iter
        write (soeo%lupri, *) "============================================="

call cpu_time (tend)
dt = tend - tstart
write (soeo%lupri, '("Time spend in soeoloop:", f15.10)') dt

end subroutine soeoloop
!=======================================================================


!> \brief Creates starting point for soeo calculation from input
!> \author C. Nygaard
!> \date Feb 10. 2011
!> \param inp Contains input read from LSDALTON.INP
!> \param Dmo The MO density matrix
!=======================================================================
subroutine soeo_setoccs (inp, Dmo, theta)

implicit none

type(soeoItem_input), intent(in) :: inp
type(matrix), intent(inout)      :: Dmo, theta

real(realk), pointer             :: tmpd(:), tmpt(:)
integer                          :: i

if (inp%cfg_unres) then
  call mem_alloc (tmpd, 2)
  call mem_alloc (tmpt, 2)
else
  call mem_alloc (tmpd, 1)
  call mem_alloc (tmpt, 1)
endif

!info:
if (inp%cfg_unres) then
  write (inp%lupri, *) 'Fractional occupations in alpha-part of the active space:'
  do i=1,inp%Nfracocc,5
  write (inp%lupri, '(5f8.5)') inp%fracoccs(i:inp%Nfracocc)
  enddo
  write (inp%lupri, *) 'Fractional occupations in beta-part of the active space:'
  do i=1,inp%Nfracoccb,5
  write (inp%lupri, '(5f8.5)') inp%fracoccsb(i:inp%Nfracoccb)
  enddo
else
  write (inp%lupri, *) 'Fractional occupations in the active space:'
  do i=1,inp%Nfracocc,5
  write (inp%lupri, '(5f8.5)') inp%fracoccs(i:inp%Nfracocc)
  enddo
endif
write (inp%lupri, *)

!sanity check:
do i=1,inp%Nfracocc
  if (inp%fracoccs(i) > 1.0E0_realk) then
    call lsquit ('Fractional occupation larger than 1',inp%lupri)
  elseif (inp%fracoccs(i) < 0.0E0_realk) then
    call lsquit ('Fractional occupation smaller than 0',inp%lupri)
  endif
enddo
do i=1,inp%Nfracoccb
  if (inp%fracoccsb(i) > 1.0E0_realk) then
    call lsquit ('Fractional occupation larger than 1',inp%lupri)
  elseif (inp%fracoccsb(i) < 0.0E0_realk) then
    call lsquit ('Fractional occupation smaller than 0',inp%lupri)
  endif
enddo

!perturbing Dmo and theta:
call mat_zero (Dmo)
do i=1,Dmo%nrow
  if (i<=inp%Nfullocc) then
    tmpd(1) = 1.0E0_realk
  elseif (i>inp%Nfullocc .and. i<=inp%Nfullocc+inp%Nfracocc) then
    tmpd(1) = inp%fracoccs(i-inp%Nfullocc)
  else
    tmpd(1) = 0.0E0_realk
  endif
  tmpt(1) = sqrt(tmpd(1)) ; tmpt(1) = acos(tmpt(1))
  if (inp%cfg_unres) then
    if (i<=inp%Nfulloccb) then
      tmpd(2) = 1.0E0_realk
    elseif (i>inp%Nfulloccb .and. i<=inp%Nfulloccb+inp%Nfracoccb) then
      tmpd(2) = inp%fracoccsb(i-inp%Nfulloccb)
    else
      tmpd(2) = 0.0E0_realk
    endif
     tmpt(2) = sqrt(tmpd(2)) ; tmpt(2) = acos(tmpt(2))
  endif
  call mat_create_ab_elms (i, i, tmpd, Dmo)
  call mat_create_ab_elms (i, 1, tmpt, theta)
enddo

call mem_dealloc (tmpd)
call mem_dealloc (tmpt)

end subroutine soeo_setoccs
!=======================================================================



!> \brief Creates the matrices used in soeo
!> \author C. Nygaard
!> \date Feb 3. 2011
!> \param Dao The input ao density matrix
!> \param Fao The input ao Fock matrix
!> \param S The input ao overlap matrix
!> \param soeo The structure containing all the matrices
!=======================================================================
subroutine soeo_init (Dao, Fao, S, soeo)

implicit none

type(matrix), intent(in)      :: Dao, Fao, S
type(soeoItem), intent(inout) :: soeo

integer                       :: Nact, Nbast, i, j, ndim
real(realk), pointer          :: orbE(:), tmp(:)

Nact  = soeo%space%Nact
Nbast = soeo%space%Nbast

if (soeo%cfg_unres) then
  call mem_alloc (orbE,2*Nbast)
  call mem_alloc (tmp,2)
  ndim = 2*Nact
else
  call mem_alloc (orbE,Nbast)
  call mem_alloc (tmp,1)
  ndim = Nact
endif

call mat_assign(soeo%mats%Dao,Dao) !Dao
call mat_assign(soeo%mats%Fao,Fao) !Fao
call mat_assign(soeo%mats%S,S)   !S

!C:
call mat_diag_f (Fao, S, orbE, soeo%mats%C)

!Fmo:
call mat_zero (soeo%mats%Fmo) 
do i=1,Nbast
  tmp(1) = orbE(i)
  if (soeo%cfg_unres) then
    tmp(2) = orbE(Nbast+i)
  endif
  call mat_create_ab_elms (i, i, tmp, soeo%mats%Fmo)
enddo

!Dmo:
call soeo_Dmo_init (S, soeo%mats%C, Dao, soeo%mats%Dmo, soeo%cfg_unres)
!oldtheta:
call soeo_oldtheta_init (soeo%mats%Dmo, soeo%mats%oldtheta, soeo%cfg_unres)

do i=1,ndim
  call soeo_get_nfirst2(soeo%space, soeo%mats%oldtheta, i, soeo%mats%nfirst(i), soeo%cfg_unres)
  do j=1,ndim
    call soeo_get_nsecond2(soeo%space, soeo%mats%oldtheta, i, j, soeo%mats%nsecond(i,j), soeo%cfg_unres)
  enddo
enddo

call mem_dealloc (orbE)
call mem_dealloc (tmp)

end subroutine soeo_init
!=======================================================================



!> \brief
!> \author C. Nygaard
!> \date Feb 3. 2011
!> \param mats Structure containing the matrices to be stored
!> \param old Structure containing the stored matrices
!=======================================================================
subroutine soeo_setold (mats, old)

implicit none

type(soeoItem_mats), intent(inout) :: mats
type(soeoItem_old), intent(inout)  :: old

integer                            :: ndim, i, j

call mat_assign(old%Dao,mats%Dao)
call mat_assign(old%Fao,mats%Fao)
call mat_assign(old%C,mats%C)
call mat_assign(old%Dmo,mats%Dmo)
call mat_assign(old%Fmo,mats%Fmo)
call mat_assign(old%oldtheta,mats%oldtheta)

ndim = size(mats%nfirst)
do i=1,ndim
  call mat_assign(old%nfirst(i),mats%nfirst(i))
  do j=1,ndim
    call mat_assign(old%nsecond(i,j),mats%nsecond(i,j))
  enddo
enddo

end subroutine soeo_setold
!=======================================================================



!> \brief Updates the matrices used in soeo
!> \author C. Nygaard
!> \date Feb 4. 2011
!> \param K The orbital rotation matrix
!> \param deltatheta The change in occupations
!> \param soeo Contains the matrices to be updated
!=======================================================================
subroutine soeo_stepforth (K, deltatheta, soeo)

implicit none

type(matrix), intent(in)      :: K, deltatheta
type(soeoItem), intent(inout) :: soeo

integer                       :: Nact, Nbast, ndim, i, j
real(realk), pointer          :: orbE(:), tmp(:)

Nact = soeo%space%Nact
Nbast = soeo%space%Nbast

call soeo_setold (soeo%mats, soeo%old)

if (soeo%cfg_unres) then
  call mem_alloc (orbE,2*Nbast)
  call mem_alloc (tmp,2)
  ndim = 2*Nact
else
  call mem_alloc (orbE,Nbast)
  call mem_alloc (tmp,1)
  ndim = Nact
endif

if (K%nrow /= Nbast .or. K%ncol /= Nbast) then
  call lsquit ('Wrong dimensions of K in soeo_stepforth', soeo%lupri)
endif
if (deltatheta%nrow /= Nact .or. deltatheta%ncol /= 1) then
  call lsquit ('Wrong dimensions of deltatheta in soeo_stepforth', soeo%lupri)
endif

!oldtheta:
call soeo_getnew_theta (deltatheta, soeo%space, soeo%mats%oldtheta, soeo%cfg_unres)
!Dmo:
call soeo_getnew_Dmo (soeo%mats%oldtheta, soeo%mats%Dmo, soeo%cfg_unres)
!C:
call soeo_getnew_C (K, soeo%mats%C)
!Dao:
call util_MO_to_AO_2 (soeo%mats%S, soeo%mats%C, soeo%mats%Dmo, soeo%mats%Dao, .false.)
!Fao, Etotal and dE:
soeo%dE = soeo%Etotal
call FCK_get_Fock (soeo%mats%Dao, soeo%mats%Fao, soeo%Etotal)
soeo%dE = soeo%Etotal - soeo%dE
!Fmo:
call util_AO_to_MO_2 (soeo%mats%S, soeo%mats%C, soeo%mats%Fao, soeo%mats%Fmo, .true.)

!nfirst and nsecond:
do i=1,ndim
  call soeo_get_nfirst2(soeo%space, soeo%mats%oldtheta, i, soeo%mats%nfirst(i), soeo%cfg_unres)
  do j=1,ndim
    call soeo_get_nsecond2(soeo%space, soeo%mats%oldtheta, i, j, soeo%mats%nsecond(i,j), soeo%cfg_unres)
  enddo
enddo

call mem_dealloc (orbE)
call mem_dealloc (tmp)

end subroutine soeo_stepforth
!=======================================================================



!> \brief Rejects the new matrices and replaces them by old ones
!> \author C. Nygaard
!> \date Feb 4. 2011
!> \param soeo The structure containing the matrices
!=======================================================================
subroutine soeo_stepback (soeo)

implicit none

type(soeoItem), intent(inout) :: soeo

integer                       :: ndim, i, j

call mat_assign(soeo%mats%Dao,soeo%old%Dao)
call mat_assign(soeo%mats%Fao,soeo%old%Fao)
call mat_assign(soeo%mats%C,soeo%old%C)
call mat_assign(soeo%mats%Dmo,soeo%old%Dmo)
call mat_assign(soeo%mats%Fmo,soeo%old%Fmo)
call mat_assign(soeo%mats%oldtheta,soeo%old%oldtheta)

ndim = size(soeo%mats%nfirst)
do i=1,ndim
  call mat_assign(soeo%mats%nfirst(i),soeo%old%nfirst(i))
  do j=1,ndim
    call mat_assign(soeo%mats%nsecond(i,j),soeo%old%nsecond(i,j))
  enddo
enddo

soeo%Etotal = soeo%Etotal - soeo%dE

end subroutine soeo_stepback
!=======================================================================

!> \brief Reads F and D from soeosave.out, and prepares them for restart
!> \author C. Nygaard
!> \date Nov. 19 2012
!> \param unres True if calculation is unrestricted, false otherwise
!> \param lupri Output lun
!> \param S Orbital overlap matrix
!> \param F Fock matrix, output
!> \param D Density matrix, output
!=======================================================================
subroutine soeo_restart (unres, lupri, S, F, D)
implicit none
logical, intent(in)         :: unres
integer, intent(in)         :: lupri
type(matrix), intent(in)    :: S
type(matrix), intent(inout) :: F, D

real(realk), pointer        :: eival(:)
type(matrix)                :: Cmo, Dmo
integer                     :: lusoeo, nbast
logical                     :: OnMaster
OnMaster = .true.

nbast = S%nrow

if (unres) then
   call mem_alloc (eival,2*nbast)
else
   call mem_alloc (eival,nbast)
endif
call mat_init (Cmo, nbast, nbast)
call mat_init (Dmo, nbast, nbast)

lusoeo = -1
call lsopen (lusoeo, "soeosave.out", "unknown", "UNFORMATTED")
rewind (lusoeo)
call mat_read_from_disk (lusoeo, F,OnMaster)
call mat_read_from_disk (lusoeo, Dmo,OnMaster)
print *, '.SOEORST specified, Fao and Dao read from soeosave.out'
write (lupri, *) '.SOEORST specified, Fao and Dao read from soeosave.out'
call lsclose (lusoeo, "KEEP")

call mat_diag_f (F, S, eival, Cmo)
call util_MO_to_AO_2 (S, Cmo, Dmo, D, .false.)

call mat_free (Cmo)
call mat_free (Dmo)
call mem_dealloc (eival)

end subroutine soeo_restart
!=======================================================================

!> \brief 
!> \author
!> \date
!> \param
!=======================================================================
subroutine soeo_isit_aufbau2 (Dao, Fao, S, unres, lupri, aufbau)

implicit none

type(matrix), intent(in)     :: Dao, Fao, S
logical, intent(in)          :: unres
integer, intent(in)          :: lupri
logical, intent(out)         :: aufbau

integer                      :: i, j, luc
integer                      :: Nbast, Nocc, Nact, Nvirt, Ndim
integer                      :: Noccb, Nactb, Nvirtb
type(matrix)                 :: SDS, Cmo, Cocc, Cvirt, Cact, Focc, Fact, Fvirt, SDSact
real(realk), pointer         :: occ(:), orbE(:)
type(matrix)                 :: tmp, eivec
real(realk), pointer         :: eival(:)
real(realk)                  :: efermi, ehomo, elumo, val
real(realk)                  :: check

logical :: OnMaster
OnMaster=.TRUE.

!initialisations
Nocc = nint(mat_TrAB (Dao, S))
Nbast = S%nrow
if (unres) then
  Nocc = Nocc/2
  Ndim = 2*Nbast
else
  Ndim = Nbast
endif
Nvirt = Nbast - Nocc
call mem_alloc (occ,Ndim)
call mem_alloc (orbE,Ndim)
call mat_init (Cmo, Nbast, Nbast)
call mat_init (SDS, Nbast, Nbast)

call mat_init (tmp, Nbast, Nbast)
call mat_mul (S, Dao, 'n', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
call mat_mul (tmp, S, 'n', 'n', -1.0E0_realk, 0.0E0_realk, SDS)
           !minus so the occupied orbitals to comes first
call mat_free (tmp)

call mat_diag_f (SDS, S, occ, Cmo)
occ = -occ !minus to get positive occupation numbers

!Create spaces
Nocc = 0 ; Nact = 0 ; Nvirt = 0
Noccb = 0 ; Nactb = 0 ; Nvirtb = 0
do i=1,Nbast
  if (occ(i)-1.0E0_realk < 1.0E-8_realk .and. occ(i)-1.0E0_realk > -1.0E-8_realk) then
    !orbital is occupied
    Nocc = Nocc + 1
  elseif (occ(i) < 1.0E-8_realk .and. occ(i) > -1.0E-8_realk) then
    !orbital is active
    Nvirt = Nvirt + 1
  else
    !orbital is active
    Nact = Nact + 1
  endif
  if (unres) then
    if (occ(i+Nbast) < 1.0E-8_realk .and.  occ(i+Nbast) > -1.0E-8_realk) then
      !orbital is occupied
      Noccb = Nocc + 1
    elseif (occ(i+Nbast)-1.0E0_realk < 1.0E-8_realk .and. occ(i+Nbast)-1.0E0_realk > -1.0E-8_realk) then
      !orbital is virtual
      Nvirt = Nvirt + 1
    else
      !orbital is active
      Nact = Nact + 1
    endif
  endif
enddo
if (Nocc+Nact+Nvirt /= Nbast) then
  call lsquit ("Something is wrong in print_orbital_info_es", lupri)
endif
if (unres) then
  Nocc = min(Nocc, Noccb)
  Nvirt = min(Nvirt, Nvirtb)
  Nact = Nbast - Nocc - Nvirt
endif

!if (Nact == 0) then
!  call print_orbital_info_ps (Dao, Fao, S, filename, unres, lupri)
!else

call mat_init (Cocc, Nbast, Nocc)
call mat_init (Focc, Nocc, Nocc)
call mat_init (Cact, Nbast, Nact)
call mat_init (Fact, Nact, Nact)
call mat_init (Cvirt, Nbast, Nvirt)
call mat_init (Fvirt, Nvirt, Nvirt)
call mat_init (SDSact, Nact, Nact)

!Make the occupied part of Fock matrix: Focc = Cocc^T Fao Cocc
if (Nocc > 0) then
  call mat_init (tmp, Nocc, Nbast)
  call mat_section (Cmo, 1, Nbast, 1, Nocc, Cocc)
  call mat_mul (Cocc, Fao, 't', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
  call mat_mul (tmp, Cocc, 'n', 'n', 1.0E0_realk, 0.0E0_realk, Focc)
  call mat_free (tmp)
  !Diagonalise Focc
  call mat_init (tmp, Nbast, Nocc)
  call mat_init (eivec, Nocc, Nocc)
  if (unres) then
    call mem_alloc (eival,2*Nocc)
  else
    call mem_alloc (eival,Nocc)
  endif
  call mat_assign(eivec,Focc)
  call mat_dsyev (eivec, eival, Nocc)
  call mat_mul (Cocc, eivec, 'n', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
  call mat_assign(Cocc,tmp)
  call mat_insert_section (Cocc, 1, Nbast, 1, Nocc, Cmo)
  orbE(1:Nocc) = eival(1:Nocc)
  if (unres) orbE(Nbast+1:Nbast+Nocc) = eival(Nocc+1:2*Nocc)
  call mem_dealloc (eival)
  call mat_free (eivec)
  call mat_free (tmp)
endif

if (Nvirt > 0) then
  !Make the unoccupied part of Fock matrix: Fvirt = Cvirt^T Fao Cvirt
  call mat_init (tmp, Nvirt, Nbast)
  call mat_section (Cmo, 1, Nbast, Nocc+Nact+1, Nbast, Cvirt)
  call mat_mul (Cvirt, Fao, 't', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
  call mat_mul (tmp, Cvirt, 'n', 'n', 1.0E0_realk, 0.0E0_realk, Fvirt)
  call mat_free (tmp)
  !Diagonalise Fvirt
  call mat_init (tmp, Nbast, Nvirt)
  call mat_init (eivec, Nvirt, Nvirt)
  if (unres) then
    call mem_alloc (eival,2*Nvirt)
  else
    call mem_alloc (eival,Nvirt)
  endif
  call mat_assign(eivec,Fvirt)
  call mat_dsyev (eivec, eival, Nvirt)
  call mat_mul (Cvirt, eivec, 'n', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
  call mat_assign(Cvirt,tmp)
  call mat_insert_section (Cvirt, 1, Nbast, Nocc+Nact+1, Nbast, Cmo)
  orbE(Nocc+Nact+1:Nbast) = eival(1:Nvirt)
  if (unres) orbE(Nbast+Nocc+Nact+1:2*Nbast) = eival(Nvirt+1:2*Nvirt)
  call mem_dealloc (eival)
  call mat_free (eivec)
  call mat_free (tmp)
endif

if (Nact > 0) then
  !Make the active part of Fock matrix: Fact = Cact^T Fao Cact
  call mat_init (tmp, Nact, Nbast)
  call mat_section (Cmo, 1, Nbast, Nocc+1, Nocc+Nact, Cact)
  call mat_mul (Cact, Fao, 't', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
  call mat_mul (tmp, Cact, 'n', 'n', 1.0E0_realk, 0.0E0_realk, Fact)
  call mat_free (tmp)
  !Diagonalise Fact
  call mat_init (tmp, Nbast, Nact)
  call mat_init (eivec, Nact, Nact)
  if (unres) then
    call mem_alloc (eival,2*Nact)
  else
    call mem_alloc (eival,Nact)
  endif
  call mat_assign(eivec,Fact)
  call mat_dsyev (eivec, eival, Nact)
  call mat_mul (Cact, eivec, 'n', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
  call mat_assign(Cact,tmp)
  call mat_free (eivec)
  call mat_free (tmp)
  !Check if all active orbitals have same energy
  efermi = eival(1)
  do i=2,Nact
    if (efermi-eival(i) < -1.0E-5_realk .or. efermi-eival(i) > 1.0E-5_realk) then
      write (lupri, *) 'WARNING: Active orbitals do not have the same energy!'
    endif
  enddo
  orbE(Nocc+1:Nocc+Nact) = eival(1:Nact)
  if (unres) orbE(Nbast+Nocc+1:Nbast+Nocc+Nact) = eival(Nact+1:2*Nact)
  call mem_dealloc (eival)
  !Make the active part of SDS
  call mat_init (tmp, Nact, Nbast)
  call mat_mul (Cact, SDS, 't', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
  call mat_mul (tmp, Cact, 'n', 'n', 1.0E0_realk, 0.0E0_realk, SDSact)
  call mat_free (tmp)
  !Diagonalise SDSact
  call mat_init (tmp, Nbast, Nact)
  call mat_init (eivec, Nact, Nact)
  if (unres) then
    call mem_alloc (eival,2*Nact)
  else
    call mem_alloc (eival,Nact)
  endif
  call mat_assign(eivec,SDSact)
  call mat_dsyev (eivec, eival, Nact)
  call mat_mul (Cact, eivec, 'n', 'n', 1.0E0_realk, 0.0E0_realk, tmp)
  call mat_assign(Cact,tmp)
  call mat_free (eivec)
  call mat_free (tmp)
  !Create active part of Cmo and occupations
  call mat_insert_section (Cact, 1, Nbast, Nocc+1, Nocc+Nact, Cmo)
  occ(Nocc+1:Nocc+Nact) = -eival(1:Nact)
  if (unres) occ(Nbast+Nocc+1:Nbast+Nocc+Nact) = eival(Nact+1:2*Nact)
  call mem_dealloc (eival)
endif

aufbau=.true.
efermi=orbE(Nocc+1)
if (Nocc > 0) then
  ehomo=orbE(Nocc)
else
  ehomo=efermi-1.0E0_realk !should just be smaller than efermi
endif
elumo=orbE(Nocc+Nact+1)
if (Nact > 0) then
  if (ehomo>efermi .or. ehomo>elumo .or. efermi>elumo) aufbau=.false.
  do i=Nocc+1,Nocc+Nact
    val = efermi - orbE(i)
    if (val < -1.0E-5_realk .or. val > 1.0E-5_realk) aufbau=.false.
  enddo
else
  if (ehomo>elumo) aufbau=.false.
endif

if (.not. aufbau) then
  write (lupri, *)
  write (lupri, *) 'Non-aufbau solution found'
  write (lupri, *) '----------------------------------------------------'

  !output orbital energies and occupations to standard output
  if (unres) then
    write (lupri, *) 'Orbital energies of occupied alpha-orbitals:'
    write (lupri, *) orbE(1:Nocc)
    write (lupri, *) 'Orbital energies of occupied beta-orbitals:'
    write (lupri, *) orbE(Nbast+1:Nbast+Nocc)
    write (lupri, *) 'Orbital occupations and energies of active alpha-orbitals:'
    do i=Nocc+1, Nocc+Nact
      write (lupri, *) i, occ(i), orbE(i)
    enddo
    write (lupri, *) 'Orbital occupations and energies of active beta-orbitals:'
    do i=Nbast+Nocc+1, Nbast+Nocc+Nact
      write (lupri, *) i, occ(i), orbE(i)
    enddo
    write (lupri, *) 'Orbital energies of virtual alpha-orbitals:'
    write (lupri, *) orbE(Nocc+Nact+1:Nbast)
    write (lupri, *) 'Orbital energies of virtual beta-orbitals:'
    write (lupri, *) orbE(Nbast+Nocc+Nact+1:2*Nbast)
  else
    write (lupri, *) 'Orbital energies of occupied orbitals:'
    write (lupri, *) orbE(1:Nocc)
    write (lupri, *) 'Orbital occupations and energies of active orbitals:'
    do i=Nocc+1, Nocc+Nact
      write (lupri, *) i, occ(i), orbE(i)
    enddo
    write (lupri, *) 'Orbital energies of virtual orbitals:'
    write (lupri, *) orbE(Nocc+Nact+1:Nbast)
  endif

  write (lupri, *) '----------------------------------------------------'
  write (lupri, *)
endif

  !closing down
  call mat_free (Cocc)
  call mat_free (Focc)
  call mat_free (Cact)
  call mat_free (Fact)
  call mat_free (Cvirt)
  call mat_free (Fvirt)
  call mat_free (SDSact)

!endif

call mat_free (Cmo)
call mat_free (SDS)
call mem_dealloc (occ)
call mem_dealloc (orbE)

end subroutine soeo_isit_aufbau2
!=======================================================================


!> \brief Testst if a given solution is an aufbau solution
!> \author C. Nygaard
!> \date May 22. 2011
!> \param soeo Structure containing all information about the matrices
!> \param aufbau True if solution is aufbau, false otherwise
!=======================================================================
subroutine soeo_isit_aufbau (soeo, aufbau)

implicit none

type(soeoItem), intent(inout) :: soeo
logical, intent(out)          :: aufbau

real(realk), pointer          :: occs(:), orbE(:), theta(:)
character, pointer            :: labels(:)
real(realk), pointer          :: ehomo(:), elumo(:), ehfmo(:), elfmo(:)
integer, pointer              :: homo(:), lumo(:), hfmo(:), lfmo(:)
integer, pointer              :: Ni(:), Na(:), Nv(:)
integer, pointer              :: high(:), low(:)
real(realk)                   :: fermi
integer                       :: u, i, j, k
logical                       :: strong

if (soeo%cfg_unres) then
  u=2
else
  u=1
endif

call mem_alloc (occs,u)
call mem_alloc (orbE,u*soeo%space%Nbast)
call mem_alloc (theta,u)
call mem_alloc (labels,u*soeo%space%Nbast)
call mem_alloc (homo,u)
call mem_alloc (ehomo,u)
call mem_alloc (lfmo,u)
call mem_alloc (elfmo,u)
call mem_alloc (hfmo,u)
call mem_alloc (ehfmo,u)
call mem_alloc (lumo,u)
call mem_alloc (elumo,u)
call mem_alloc (Ni,u)
call mem_alloc (Na,u)
call mem_alloc (Nv,u)
call mem_alloc (high,u)
call mem_alloc (low,u)

!label the orbitals inactive 'i', active 'a' and virtual 'v'
Ni = 0 ; Na = 0 ; Nv = 0
do i=1,soeo%space%Nbast
  call mat_get_ab_elms (soeo%mats%Dmo, i, i, occs)
  do j=1,u
    if (j==1) then
      k = i !alpha part
    else
      k = soeo%space%Nbast + i !beta part
    endif

    if (occs(j)-1.0E0_realk<1.0E-5_realk .and. occs(j)-1.0E0_realk>-1.0E-5_realk) then !occ=1
      labels(k) = 'i'
      Ni(j) = Ni(j) + 1
    elseif (occs(j)<1.0E-5_realk .and. occs(j)>-1.0E-5_realk) then !occ=0
      labels(k) = 'v'
      Nv(j) = Nv(j) + 1
    else !occ=frac
      labels(k) = 'a'
      Na(j) = Na(j) + 1
    endif
  enddo
enddo

write(soeo%lupri, *) 'alpha:', labels(1:soeo%space%Nbast)
if (soeo%cfg_unres) then
  write(soeo%lupri, *) 'beta:', labels(soeo%space%Nbast+1:2*soeo%space%Nbast)
endif

!finding orbital energies
call mat_diag_f (soeo%mats%Fao, soeo%mats%S, orbE, soeo%mats%C)

!finding the homo and the lumo (only fully- or non-occupied orbitals)
homo = 0 ; ehomo = -1.0E9_realk
lfmo = 0 ; elfmo = 1.0E9_realk
hfmo = 0 ; ehfmo = -1.0E9_realk
lumo = 0 ; elumo = 1.0E9_realk
do i=1,soeo%space%Nbast
  do j=1,u
    if (j==1) then
      k=i !alpha
    else
      k=soeo%space%Nbast+i !beta
    endif

    if (labels(k) == 'i' .and. orbE(k) > ehomo(j)) then
      homo(j) = i
      ehomo(j) = orbE(k)
    endif
    if (labels(k) == 'v' .and. orbE(k) < elumo(j)) then
      lumo(j) = i
      elumo(j) = orbE(k)
    endif
    if (labels(k) == 'a') then
      if (orbE(k) < elfmo(j)) then
        lfmo(j) = i
        elfmo(j) = orbE(k)
      endif
      if (orbE(k) > ehfmo(j)) then
        hfmo(j) = i
        ehfmo(j) = orbE(k)
      endif
    endif
  enddo
enddo

do j=1,u
  if (Ni(j)==0) then
    homo(j) = 0
  endif
  if (Nv(j)==0) then
    lumo(j) = soeo%space%Nbast+1
  endif
enddo

!finding the space of orbitals to be changed
aufbau = .true.
strong = .false.
if (Na(1) /= 0) then
  fermi = orbE(Ni(1)+1)
else
  if (soeo%cfg_unres) then
    if (Na(2) /= 0) then
      fermi = orbE(soeo%space%Nbast+Ni(2)+1)
    endif
  endif
endif
do j=1,u
write(soeo%lupri, *) 'u, j =', u, j
  !are there fractionally occupied orbitals?
  if (lfmo(j)/=0 .and. hfmo(j)/=0) then !yes, fon's are present
write(soeo%lupri, *) 'fons present'
    !are inactive nested?
    if (Ni(j)==homo(j) .or. Ni(j)==0) then !yes, inactive are nested
write(soeo%lupri, *) 'inactive are nested'
      !are virtual nested?
      if (Nv(j)==soeo%space%Nbast-lumo(j)+1 .or. Nv(j)==0) then !yes, virtual are nested
write(soeo%lupri, *) 'virtuals are nested'
        !do active have the same energy?
        do i=Ni(j)+1,Ni(j)+Na(j)
          if (orbE(i)-fermi < 1.0E-4_realk .and. orbE(i)-fermi > -1.0E-4_realk) then !yes, same energy
            low(j) = lumo(j)
            high(j) = homo(j)
            !aufbau solution for the j'th part
          else !no, not same energy
write(soeo%lupri, *) 'fons do NOT have same energy'
            aufbau = .false.
            low(j) = lfmo(j)
            high(j) = hfmo(j)
          endif
        enddo
      else !no, virtual are not nested
write(soeo%lupri, *) 'virtuals are NOT nested'
        aufbau = .false.
        strong = .true.
        low(j) = min(lfmo(j), lumo(j))
        high(j) = hfmo(j)
      endif
    else !no, inactive are not nested
write(soeo%lupri, *) 'inactive are NOT nested'
      aufbau = .false.
      strong = .true.
      low(j) = min(lfmo(j), lumo(j))
      high(j) = max(hfmo(j), homo(j))
    endif
  else !no, only integral occupations are present
write(soeo%lupri, *) 'only integral occupations'
    !is homo and lumo in right order?
    if (homo(j) < lumo(j)) then !yes
write(soeo%lupri, *) 'homo lumo has the right order'
      low(j) = lumo(j)
      high(j) = homo(j)
      !aufbau solution for the j'th part
    else !no
write(soeo%lupri, *) 'homo lumo has the WRONG order'
      aufbau = .false.
      strong = .true.
      low(j) = lumo(j)
      high(j) = homo(j)
    endif
  endif
enddo

write(soeo%lupri, *) 'first low =', low
write(soeo%lupri, *) 'first high =', high


!if unres: connect alpha and beta spaces
if (u==2) then
  do
write (soeo%lupri, *)
    if (orbE(high(1))>=orbE(soeo%space%Nbast+high(2)+1)) then
write (soeo%lupri, *) 'high1 >= high2+1'
write (soeo%lupri, *) 'high1 =', orbE(high(1))
write (soeo%lupri, *) 'high2+1 =', orbE(soeo%space%Nbast+high(2)+1)
      high(2) = high(2) + 1
    elseif (orbE(soeo%space%Nbast+high(2))>=orbE(high(1)+1)) then
write (soeo%lupri, *) 'high2 >= high1+1'
write (soeo%lupri, *) 'high2 =', orbE(soeo%space%Nbast+high(2))
write (soeo%lupri, *) 'high1+1 =', orbE(high(1)+1)
      high(1) = high(1) + 1
    elseif (orbE(low(1))<=orbE(soeo%space%Nbast+low(2)-1)) then
write (soeo%lupri, *) 'low1 <= low2-1'
write (soeo%lupri, *) 'low1 =', orbE(low(1))
write (soeo%lupri, *) 'high2-1 =', orbE(soeo%space%Nbast+low(2)-1)
      low(2) = low(2) - 1
    elseif (orbE(soeo%space%Nbast+low(2))<=orbE(low(1)-1)) then
write (soeo%lupri, *) 'low2 <= low1-1'
write (soeo%lupri, *) 'low2 =', orbE(soeo%space%Nbast+low(2))
write (soeo%lupri, *) 'high1-1 =', orbE(low(1)-1)
      low(1) = low(1) - 1
    else
      exit
    endif
write (soeo%lupri, *) 'new low =', low
write (soeo%lupri, *) 'new high =', high
  enddo
endif

write (soeo%lupri, *) 'final low =', low
write (soeo%lupri, *) 'final high =', high

if (.not. aufbau) then
  if (strong) then
    write (soeo%lupri, *) 'strong non-aufbau'
    call soeo_perturboccs_strong (soeo, orbE, low, high, labels)
  else
    write (soeo%lupri, *) 'weak non-aufbau'
    call soeo_perturboccs_weak (soeo, orbE, low, high, labels)
  endif
endif

call mem_dealloc (occs)
call mem_dealloc (orbE)
call mem_dealloc (theta)
call mem_dealloc (labels)
call mem_dealloc (homo)
call mem_dealloc (ehomo)
call mem_dealloc (lfmo)
call mem_dealloc (elfmo)
call mem_dealloc (hfmo)
call mem_dealloc (ehfmo)
call mem_dealloc (lumo)
call mem_dealloc (elumo)
call mem_dealloc (Ni)
call mem_dealloc (Na)
call mem_dealloc (Nv)
call mem_dealloc (low)
call mem_dealloc (high)

end subroutine soeo_isit_aufbau
!=======================================================================



!> \brief perturbs the occupations away from strong non-aufbau solution
!> \author C. Nygaard
!> \date May 26 2011
!> \param soeo Structure containing the matrices to be changed
!> \param low The lowest orbital to be changed
!> \param high The highest orbital to be changed
!> \param labels contains labels that mark the status of the orbitals
!=======================================================================
subroutine soeo_perturboccs_strong (soeo, orbE, low, high, labels)

implicit none

type(soeoItem), intent(inout) :: soeo
real(realk), intent(in)       :: orbE(:)
integer, intent(inout)        :: low(:), high(:)
character, intent(in)         :: labels(:)

real(realk), pointer          :: occ(:), theta(:)
real(realk), pointer          :: occh(:), occl(:)
real(realk)                   :: Nelec, energy, elow, ehigh
integer                       :: u, i, j, k, N
real(realk)                   :: diffh, diffl

if (soeo%cfg_unres) then
  u=2
else
  u=1
endif

call mem_alloc (occ,u)
call mem_alloc (theta,u)
call mem_alloc (occh,u)
call mem_alloc (occl,u)

write (soeo%lupri, *) 'perturboccs: low =', low
write (soeo%lupri, *) 'perturboccs: high =', high
write (soeo%lupri, *) 'alpha labels =', labels(low(1):high(1))
if (soeo%cfg_unres) then
  write (soeo%lupri, *) 'beta labels  =', labels(soeo%space%Nbast+low(2):soeo%space%Nbast+high(2))
endif

Nelec = 0.0E0_realk ; energy = 0.0E0_realk
N = 0
do j=1,u
!  if (low(j) /= high(j)) then
    do i=low(j),high(j)
      if (j==1) then
        k=i
      else
        k=soeo%space%Nbast+i
      endif
      call mat_get_ab_elms (soeo%mats%Dmo, i, i, occ)
      Nelec = Nelec + occ(j)
      N = N + 1
      energy = energy + occ(j)*orbE(k)
    enddo
!  endif
enddo

energy = energy/Nelec

write (soeo%lupri, *)
write (soeo%lupri, *) 'energy =', energy
write (soeo%lupri, *) 'Nelec =', Nelec

do j=1,u
  if (j==1) then
    ehigh = orbE(high(j))
    elow = orbE(low(j))
  else
    ehigh = orbE(soeo%space%Nbast+high(j))
    elow = orbE(soeo%space%Nbast+low(j))
  endif

write (soeo%lupri, *) 'j =', j

    do
      call mat_get_ab_elms (soeo%mats%Dmo, high(j), high(j), occh)
      call mat_get_ab_elms (soeo%mats%Dmo, low(j), low(j), occl)

write (soeo%lupri, *) 'high, old occh =', high(j), occh(j)
write (soeo%lupri, *) 'low, old occl =', low(j), occl(j)

      occh(j) = Nelec/N
      occl(j) = Nelec/N
!  
!write (soeo%lupri, *) 'high, intermediate occh =', high(j), occh(j)
!write (soeo%lupri, *) 'low, intermediate occl =', low(j), occl(j)

!      !difference between occh and 0
!      diffh = occh(j)
!      !difference between occl and 1
!      diffl = 1-occl(j)
!      !minimum difference / 10 is the perturbation in the occupation numbers
!      occh(j) = occh(j)-min(diffh,diffl)/10.0E0_realk
!      occl(j) = occl(j)+min(diffh,diffl)/10.0E0_realk

      if (occl(j) > 0.5E0_realk) then
        occh(j) = occh(j) - 0.5E0_realk*(1.0E0_realk-occl(j))
        occl(j) = occl(j) + 0.5E0_realk*(1.0E0_realk-occl(j))
      else
        occl(j) = occl(j) + 0.5E0_realk*occh(j)
        occh(j) = occh(j) - 0.5E0_realk*occh(j)
      endif

!      !distribute electrons so 60% are equally distributed in first half of active orbitals
!      ! and 40% in second half of active orbitals
!      if (i<=(low(j)+high(j)-1)/2) then
!        occ(j) = (0.6E0_realk*Nelec/u) / ((high(j)-low(j)+1)/2)
!      else
!        occ(j) = (0.4E0_realk*Nelec/u) / ((high(j)-low(j)+1) - (high(j)-low(j)+1)/2)
!      endif

!      occ(j) = Nelec / N !electrons gets equally distributed

!      occ(j) = occ(j) + ((energy - orbE(k))/(ehigh-elow))*occ(j)
  
!      if (orbE(k) > energy) then
!        occ(j) = 0.5E0_realk*occ(j)
!      elseif (orbE(k) < energy) then
!        occ(j) = occ(j) + 0.5E0_realk*(1.0E0_realk-occ(j))
!      endif

!write (soeo%lupri, *) 'orbE =', orbE(k)
write (soeo%lupri, *) 'high, new occh =', high(j), occh(j)
write (soeo%lupri, *) 'low, new occl =', low(j), occl(j)
  
      call mat_get_ab_elms (soeo%mats%oldtheta, high(j), 1, theta)
      theta(j) = sqrt(occh(j))
      theta(j) = acos(theta(j))
      call mat_create_ab_elms (high(j), 1, theta, soeo%mats%oldtheta)
      call mat_get_ab_elms (soeo%mats%oldtheta, low(j), 1, theta)
      theta(j) = sqrt(occl(j))
      theta(j) = acos(theta(j))
      call mat_create_ab_elms (low(j), 1, theta, soeo%mats%oldtheta)

      high(j) = high(j) - 1
      low(j) = low(j) + 1
      if (high(j) < low(j)) then
        write (soeo%lupri, *) 'high < low :', high(j) - low(j)
        exit
      elseif (high(j) == low(j)) then
        write (soeo%lupri, *) 'high = low =', high(j)
        call mat_get_ab_elms (soeo%mats%Dmo, high(j), high(j), occ)

        occ(j) = Nelec/N
        call mat_get_ab_elms (soeo%mats%oldtheta, high(j), 1, theta)
        theta(j) = sqrt(occ(j))
        theta(j) = acos(theta(j))
        call mat_create_ab_elms (high(j), 1, theta, soeo%mats%oldtheta)

!        if (occ(j) < occh(j)) then
!          occl(j) = occ(j)
!          high(j) = high(j) + 1
!        else
!          occh(j) = occ(j)
!          low(j) = low(j) - 1
!        endif
!
!write (soeo%lupri, *) 'high, old occh =', high(j), occh(j)
!write (soeo%lupri, *) 'low, old occl =', low(j), occl(j)
!
!        !difference between occh and 0
!        diffh = occh(j)
!        !difference between occl and 1
!        diffl = 1-occl(j)
!        !minimum difference / 10 is the perturbation in the occupation numbers
!        occh(j) = occh(j)-min(diffh,diffl)/10.0E0_realk
!        occl(j) = occl(j)+min(diffh,diffl)/10.0E0_realk
!
!write (soeo%lupri, *) 'high, new occh =', high(j), occh(j)
!write (soeo%lupri, *) 'low, new occl =', low(j), occl(j)
!
!        call mat_get_ab_elms (soeo%mats%oldtheta, high(j), 1, theta)
!        theta(j) = sqrt(occh(j))
!        theta(j) = acos(theta(j))
!        call mat_create_ab_elms (high(j), 1, theta, soeo%mats%oldtheta)
!        call mat_get_ab_elms (soeo%mats%oldtheta, low(j), 1, theta)
!        theta(j) = sqrt(occl(j))
!        theta(j) = acos(theta(j))
!        call mat_create_ab_elms (low(j), 1, theta, soeo%mats%oldtheta)
        exit
      endif

!    call mat_create_ab_elms (i, i, occ, soeo%mats%Dmo)
    enddo
!  endif
enddo

call mem_dealloc (occ)
call mem_dealloc (theta)
call mem_dealloc (occh)
call mem_dealloc (occl)

end subroutine soeo_perturboccs_strong
!=======================================================================

!> \brief perturbs the occupations away from weak non-aufbau solution
!> \author C. Nygaard
!> \date May 26 2011
!> \param soeo Structure containing the matrices to be changed
!> \param low The lowest orbital to be changed
!> \param high The highest orbital to be changed
!> \param labels contains labels that mark the status of the orbitals
!=======================================================================
subroutine soeo_perturboccs_weak (soeo, orbE, low, high, labels)

implicit none

type(soeoItem), intent(inout) :: soeo
real(realk), intent(in)       :: orbE(:)
integer, intent(inout)        :: low(:), high(:)
character, intent(in)         :: labels(:)

real(realk), pointer          :: occ(:), theta(:)
real(realk), pointer          :: occh(:), occl(:)
real(realk)                   :: Nelec, energy, elow, ehigh
integer                       :: u, i, j, k, N
real(realk)                   :: diffh, diffl

if (soeo%cfg_unres) then
  u=2
else
  u=1
endif

call mem_alloc (occ,u)
call mem_alloc (theta,u)
call mem_alloc (occh,u)
call mem_alloc (occl,u)

write (soeo%lupri, *) 'perturboccs: low =', low
write (soeo%lupri, *) 'perturboccs: high =', high
write (soeo%lupri, *) 'alpha labels =', labels(low(1):high(1))
if (soeo%cfg_unres) then
  write (soeo%lupri, *) 'beta labels  =', labels(soeo%space%Nbast+low(2):soeo%space%Nbast+high(2))
endif

Nelec = 0.0E0_realk ; energy = 0.0E0_realk
N = 0
do j=1,u
  do i=low(j),high(j)
    if (j==1) then
      k=i
    else
      k=soeo%space%Nbast+i
    endif
    call mat_get_ab_elms (soeo%mats%Dmo, i, i, occ)
    Nelec = Nelec + occ(j)
    N = N + 1
    energy = energy + occ(j)*orbE(k)
  enddo
enddo

energy = energy/Nelec

write (soeo%lupri, *)
write (soeo%lupri, *) 'energy =', energy
write (soeo%lupri, *) 'Nelec =', Nelec

do j=1,u
  if (j==1) then
    ehigh = orbE(high(j))
    elow = orbE(low(j))
  else
    ehigh = orbE(soeo%space%Nbast+high(j))
    elow = orbE(soeo%space%Nbast+low(j))
  endif

write (soeo%lupri, *) 'j =', j

  do
    call mat_get_ab_elms (soeo%mats%Dmo, high(j), high(j), occh)
    call mat_get_ab_elms (soeo%mats%Dmo, low(j), low(j), occl)

write (soeo%lupri, *) 'high, old occh =', high(j), occh(j)
write (soeo%lupri, *) 'low, old occl =', low(j), occl(j)

    !difference between occh and 0
    diffh = occh(j)
    !difference between occl and 1
    diffl = 1-occl(j)
    !minimum difference / 10 is the perturbation in the occupation numbers
    occh(j) = occh(j)-min(diffh,diffl)/10.0E0_realk
    occl(j) = occl(j)+min(diffh,diffl)/10.0E0_realk

write (soeo%lupri, *) 'high, new occh =', high(j), occh(j)
write (soeo%lupri, *) 'low, new occl =', low(j), occl(j)
  
    call mat_get_ab_elms (soeo%mats%oldtheta, high(j), 1, theta)
    theta(j) = sqrt(occh(j))
    theta(j) = acos(theta(j))
    call mat_create_ab_elms (high(j), 1, theta, soeo%mats%oldtheta)
    call mat_get_ab_elms (soeo%mats%oldtheta, low(j), 1, theta)
    theta(j) = sqrt(occl(j))
    theta(j) = acos(theta(j))
    call mat_create_ab_elms (low(j), 1, theta, soeo%mats%oldtheta)

    high(j) = high(j) - 1
    low(j) = low(j) + 1
    if (high(j) < low(j)) then
      write (soeo%lupri, *) 'high < low :', high(j) - low(j)
      exit
    elseif (high(j) == low(j)) then
      write (soeo%lupri, *) 'high = low =', high(j)
      call mat_get_ab_elms (soeo%mats%Dmo, high(j), high(j), occ)

      if (occ(j) < occh(j)) then
        occl(j) = occ(j)
        high(j) = high(j) + 1
      else
        occh(j) = occ(j)
        low(j) = low(j) - 1
      endif

write (soeo%lupri, *) 'high, old occh =', high(j), occh(j)
write (soeo%lupri, *) 'low, old occl =', low(j), occl(j)

      !difference between occh and 0
      diffh = occh(j)
      !difference between occl and 1
      diffl = 1-occl(j)
      !minimum difference / 10 is the perturbation in the occupation numbers
      occh(j) = occh(j)-min(diffh,diffl)/10.0E0_realk
      occl(j) = occl(j)+min(diffh,diffl)/10.0E0_realk

write (soeo%lupri, *) 'high, new occh =', high(j), occh(j)
write (soeo%lupri, *) 'low, new occl =', low(j), occl(j)

      call mat_get_ab_elms (soeo%mats%oldtheta, high(j), 1, theta)
      theta(j) = sqrt(occh(j))
      theta(j) = acos(theta(j))
      call mat_create_ab_elms (high(j), 1, theta, soeo%mats%oldtheta)
      call mat_get_ab_elms (soeo%mats%oldtheta, low(j), 1, theta)
      theta(j) = sqrt(occl(j))
      theta(j) = acos(theta(j))
      call mat_create_ab_elms (low(j), 1, theta, soeo%mats%oldtheta)
      exit
    endif
  enddo
enddo

call mem_dealloc (occ)
call mem_dealloc (theta)
call mem_dealloc (occh)
call mem_dealloc (occl)

end subroutine soeo_perturboccs_weak
!=======================================================================


!> \brief Output results from a SOEO calculation
!> \author C. Nygaard
!> \date May 27 2011
!> \param soeo Structure containing the information to be outputtet
!=======================================================================
subroutine soeo_output_results (soeo)

implicit none

type(soeoItem), intent(inout) :: soeo
real(realk) :: t1, t2

write (soeo%lupri, *) 'Final SOEO energy', soeo%Etotal
call print_orbital_info2(soeo%mats%Dao, soeo%mats%Fao, soeo%mats%S, &
    &'cmo.out', soeo%cfg_unres, .false., soeo%lupri)

if (soeo%test) then
  write (soeo%lupri, *) 'TEST: Final SOEO energy', soeo%Etotal
  call mat_get_elm (soeo%mats%oldtheta, 1, 1, t1)
  call mat_get_elm (soeo%mats%oldtheta, 2, 1, t2)
  write (soeo%lupri, *) 'TEST: Final occupation 1', dcos(t1)**2
  write (soeo%lupri, *) 'TEST: Final occupation 2', dcos(t2)**2
endif

!real(realk), allocatable   :: orbE(:), occs(:)
!integer                    :: i, j, m, l, u
!
!if (soeo%cfg_unres) then
!  u = 2
!else
!  u = 1
!endif
!allocate (orbE(u*soeo%space%Nbast), occs(u*soeo%space%Nbast))
!
!call fck_get_fock (soeo%mats%Dao, soeo%mats%Fao, soeo%Etotal)
!call mat_diag_f (soeo%mats%Fao, soeo%mats%S, orbE, soeo%mats%C)
!write (soeo%lupri, *) "SOEO-energy (fck_get_fock)=", soeo%Etotal
!
!j=0 ; m=0
!do i=1,soeo%space%Nbast
!  occs(i) = soeo_diag_elma (soeo%mats%Dmo, i)
!  if (occs(i) > 1.0E-8_realk) then
!    j = j + 1
!    l = i
!  endif
!  if(soeo%cfg_unres) then
!    occs(soeo%space%Nbast+i) = soeo_diag_elmb (soeo%mats%Dmo, i)
!    if (occs(soeo%space%Nbast+i) > 1.0E-8_realk) then
!      m = m + 1
!      u = i
!    endif
!  endif
!enddo
!
!if (soeo%cfg_unres) then
!  write (soeo%lupri, *) ' (a)Orb energy     (a)Orb occ  (b)Orb energy     (b)Orb occ  '
!else
!   write (soeo%lupri, *) 'Orb energy     Orb occ'
!endif
!do i=1,soeo%space%Nbast
!  if (soeo%cfg_unres) then
!    write (soeo%lupri, '(4f15.5)') orbE(i), occs(i),&
!                                 & orbE(soeo%space%Nbast+i), occs(soeo%space%Nbast+i)
!  else
!    write (soeo%lupri, *) orbE(i), occs(i)
!  endif
!enddo
!
!write (soeo%lupri, *)
!write (soeo%lupri, *) 'MO-coefficient matrix for the occupied orbitals'
!write (soeo%lupri, *) ' (occupied => occupation > 1.0E-8_realk)'
!!if (l > u) then
!!  call mat_print (soeo%mats%C, 1, soeo%space%Nbast, 1, l, soeo%lupri)
!!else
!!  call mat_print (soeo%mats%C, 1, soeo%space%Nbast, 1, u, soeo%lupri)
!!endif
!!u = soeo%space%Nocc + soeo%space%Nact
!u = soeo%space%Nbast
!call mat_print (soeo%mats%C, 1, soeo%space%Nbast, 1, u, soeo%lupri)

!deallocate (orbE, occs)

end subroutine soeo_output_results
!=======================================================================

!> \brief Calculates electronic part of dipole moment
!> \author C. Nygaard
!> \date Nov. 20 2012
!> \param V Matrices containing dipole integrals in x, y and z
!> \param D Density matrix in AO basis
!> \param unres True if calculation is unrestricted
!> \lupri Output lun
!=======================================================================
subroutine soeo_dipole (V, D, unres, lupri)

implicit none

type(matrix), intent(in) :: V(:), D
logical, intent(in)      :: unres
integer, intent(in)      :: lupri

real(realk)              :: dip(3)
integer                  :: i

if (size(V)/=3) call lsquit ('Wrong dimensions of V in soeo_dipole!', lupri)

dip = 0.0E0_realk
do i=1,3
  dip(i) = -mat_TrAB (V(i), D)
enddo
if (.not. unres) then
  dip = 2.0E0_realk * dip
endif

write (lupri, '("Electronic contribution to dipole moment in a.u.")')
write (lupri, '("x ", e15.6)') dip(1)
write (lupri, '("y ", e15.6)') dip(2)
write (lupri, '("z ", e15.6)') dip(3)

end subroutine soeo_dipole
!=======================================================================

!> \brief Calculates polarizability
!> \author C. Nygaard
!> \date Nov. 20 2012
!> \param 
!=======================================================================
subroutine soeo_polarizability (soeo, V)

implicit none

type(soeoItem), intent(inout) :: soeo
type(matrix), intent(in)      :: V(:)

real(realk)                   :: pol(3,3), mu
integer                       :: i, j
type(matrix)                  :: gpm(3), gpv(3)
type(matrix)                  :: vm1(3), vv1(3)
type(matrix)                  :: tmpm, tmpv

if (size(V)/=3) call lsquit ('Wrong dimensions of V in soeo_polarizability!',&
                           & soeo%lupri)

do i=1,3
  call mat_init (gpm(i), soeo%space%Nbast, soeo%space%Nbast)
  call mat_init (vm1(i), soeo%space%Nbast, soeo%space%Nbast)
  call mat_init (gpv(i), soeo%space%Nact, 1)
  call mat_init (vv1(i), soeo%space%Nact, 1)
enddo
call mat_init (tmpm, soeo%space%Nbast, soeo%space%Nbast)
call mat_init (tmpv, soeo%space%Nact, 1)

do i=1,3
  call soeo_get_propgrad (soeo, V(i), gpm(i), gpv(i))
  call soeo_solver (soeo, vm1(i), vv1(i), mu, gpm(i), gpv(i), .true.)
enddo

do i=1,3
  call soeo_linear_transform (soeo, vm1(i), vv1(i), tmpm, tmpv)
  do j=1,3
    pol(j,i) = soeo_dotproduct (vm1(j), vv1(j), tmpm, tmpv)
  enddo
enddo

!debug
do i=1,3
  write (soeo%lupri, *) ''
  write (soeo%lupri, *) '-----------------------------------------------'
  write (soeo%lupri, *) 'i =', i
  write (soeo%lupri, *) '-----------------------------------------------'
  write (soeo%lupri, *) 'V(i)'
  call mat_print (V(i), 1, soeo%space%Nbast, 1, soeo%space%Nbast, soeo%lupri)
  write (soeo%lupri, *) '-----------------------------------------------'
  write (soeo%lupri, *) 'gpm(i)'
  call mat_print (gpm(i), 1, soeo%space%Nbast, 1, soeo%space%Nbast, soeo%lupri)
  write (soeo%lupri, *) '-----------------------------------------------'
  write (soeo%lupri, *) 'gpv(i)'
  call mat_print (gpv(i), 1, soeo%space%Nact, 1, 1, soeo%lupri)
  write (soeo%lupri, *) '-----------------------------------------------'
  write (soeo%lupri, *) 'vm1(i)'
  call mat_print (vm1(i), 1, soeo%space%Nbast, 1, soeo%space%Nbast, soeo%lupri)
  write (soeo%lupri, *) '-----------------------------------------------'
  write (soeo%lupri, *) 'vv1(i)'
  call mat_print (vv1(i), 1, soeo%space%Nact, 1, 1, soeo%lupri)
enddo
!end debug

write (soeo%lupri, '("Electronic contribution to polarizability in a.u.")')
write (soeo%lupri, '("  ", 3a15)') "x", "y", "z"
write (soeo%lupri, '("x ", 3e15.6)') pol(1,1), pol(1,2), pol(1,3)
write (soeo%lupri, '("y ", 3e15.6)') pol(2,1), pol(2,2), pol(2,3)
write (soeo%lupri, '("z ", 3e15.6)') pol(3,1), pol(3,2), pol(3,3)

!debug polarizability
!calculate the change in occupation numbers along v1(z) to see if this 
! can explain the difference in polarizability btw e and ps
!end debug

do i=1,3
  call mat_free (gpm(i))
  call mat_free (vm1(i))
  call mat_free (gpv(i))
  call mat_free (vv1(i))
enddo
call mat_free (tmpm)
call mat_free (tmpv)

end subroutine soeo_polarizability
!=======================================================================

end module soeo_loop
