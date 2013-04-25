!> @file 
!> Contains Augmented Roothaan-Hall / Direct Density optimization drivers.

!> \brief Augmented RH/Direct dens driver.
!> \author S. Host
!> \date 2007
!>
!> ARH: S. Host, B. Jansik, J. Olsen et al. PCCP 10, 5344 (2008)   \n
!>      S. Host, J. Olsen, B. Jansik et al. JCP 129, 124106 (2008) \n
!>                                                                 \n
!> Direct density: P. Salek, S. Host, L. Thogersen et al. JCP 126, 114110 (2007)
!>
MODULE ARHmodule
   use direct_dens_util
   use files
   use queue_ops
   use arhDensity
   use dal_interface
   use KS_settings, only: SaveF0andD0
   use davidson_settings
   use driver
   use II_XC_interfaceModule
   use IntegralInterfaceMOD
contains

   !> \brief This routine calls the appropriate ARH solver 
   !> \author S. Host
   !> \date 2007
   !>
   !> Calculate the density that minimizes the Eepsilon = 2Tr(F*D(X)) + ARH terms. \n
   !> D ,the start guess, is updated during the optimization
   !> the optimized density is returned in D. .\n
   !> This routine is a bit obsolete (made more sense when it was possible to do
   !> more than one Newton iteration). Could be merged with driver.
   !>
   subroutine arh_get_density(arh,decomp,F,fifoqueue,Dnew,SCF_iteration,davidCFG,H1,ls)
   use matrix_util
   use precision
    implicit none
    type(lsitem) ::ls
    !> Contains solver info (ARH/TrFD)
    type(solverItem),intent(inout) :: arh
    !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
    type(decompItem),intent(inout) :: decomp
    !> Current Fock/KS matrix 
    type(Matrix), intent(in)       :: F
    !> Contains Fock/KS and density matrices from previous SCF iterations
    type(modFIFO),intent(inout)    :: fifoqueue
    !> Input: Current density matrix. Output: New density matrix.
    type(matrix), intent(inout)    :: Dnew 
    !> One electron Hamiltonian
    type(matrix),intent(in) :: H1
    !> Current SCF iteration
    integer, intent(in)            :: SCF_iteration
    !> solver info for davidson solver
    type(RedSpaceItem),intent(inout)  :: davidCFG
    type(Matrix)                   :: x,wrk,wrk2,hes
    type(Matrix),target            :: P
    real(realk)                    :: xnorm
    integer                        :: ndim, ndens, i, hesdim,j
    real(realk)                    :: t1, t2
!Print variables
      real(realk), pointer        :: weights(:)
      real(realk)                 :: actual_change_norm, expanded_change_norm
      type(Matrix), pointer       :: Dpointer, Fpointer, D0, Di, dummyp, Fi, F0
      type(debugItem)             :: debug
      type(DDitem)                :: DD



    ndim = decomp%S%nrow
    call mat_init(x,ndim,ndim)
    call mat_init(wrk,ndim,ndim)
    call mat_init(wrk2,ndim,ndim)

     
    !write (arh%lupri,*) 'Incoming density, arh_get_density, OAO basis:'
    !call MAT_PRINT(D, 1, D%nrow, 1, D%ncol, arh%LUPRI)

    if (arh%set_arhterms) then
       ndens = fifoqueue%offset
       if (ndens > 0) then 
          allocate(arh%fifometric(ndens,ndens))
          allocate(arh%inv_fifometric(ndens,ndens))
          allocate(arh%fifoM(ndens,ndens))
          call fifo_inverse_metric(arh,fifoqueue)
          call arh_get_M(arh,fifoqueue)
       endif
    endif
 
     !write(arh%lupri,*) 'F, AO:'
     !call mat_print(F,1,ndim,1,ndim,arh%lupri)
     !write(arh%lupri,*) 'D, AO:'
     !call mat_print(Dnew,1,ndim,1,ndim,arh%lupri)

    if (arh%step_accepted) call get_oao_transformed_matrices(decomp,F,Dnew)

    !FIXME: repair debug stuff for new input structure
    if (arh%debug_dd) then
       arh%scfit = SCF_iteration
       call dd_debug_homolumo(decomp,arh%diag_hlgap)
       if (SCF_iteration > arh%cfg_nits_debug) then
          call DD_init(decomp,DD)
          call DD_homolumo_and_heseigen(DD,decomp,debug,.false.,0,fifoqueue=fifoqueue)
          arh%iter_hlgap = debug%iter_hlgap
          arh%heseival   = debug%heseival
          arh%arheival   = debug%arheival
          call DD_shutdown(decomp,DD)
       endif
    else if (arh%debug_dd_homolumo) then
       call dd_debug_homolumo(decomp,debug%diag_hlgap)
    endif

    if (arh%debug_hessian) then
       if (decomp%cfg_unres) call lsquit('Debug routine not tested for unrestricted',decomp%lupri)
       hesdim = ndim*(ndim+1)/2 - ndim
       call mat_init(hes,hesdim,hesdim)
       call debug_get_hessian(arh,decomp,fifoqueue,hes)
       call util_diag(arh%lupri,hes,.false.,0.25E0_realk,'Lowest Hessian eigenvalue:')
       call mat_free(hes)
    endif

    !write(arh%lupri,*) 'FU:'
    !call mat_print(decomp%FU,1,ndim,1,ndim,arh%lupri)
    !write(arh%lupri,*) 'DU:'
    !call mat_print(decomp%DU,1,ndim,1,ndim,arh%lupri)
    call get_OAO_gradient(decomp%FU, decomp%DU, wrk) !wrk = gradient
    call mat_scal(0.25E0_realk,wrk) !To match linear transformation, also divided by 4!

    arh%OAO_gradnrm = sqrt(mat_sqnorm2(wrk)) 
    arh%OAO_gradnrm_exist = .true.
    
    write (arh%lupri,*) 'OAO gradnorm', arh%OAO_gradnrm
    davidCFG%arh_gradnorm=arh%OAO_gradnrm

    if (arh%cfg_nodamp) then
       call arh_PCG(arh, decomp, wrk, x, 0.0E0_realk, antisymmetric,fifoqueue)
    else if (arh%cfg_arh_truncate .or. arh%cfg_arh_crop) then
       if (davidCFG%arh_davidson) then
            davidCFG%symm=antisymmetric
            CALL LSTIMER('START ',t1,t2,decomp%lupri)
	    call arh_davidson_solver(davidCFG,arh,decomp,wrk,X,SCF_iteration,H1,wrk2,ls) 
            CALL LSTIMER('ARH DAVIDS. ',t1,t2,arh%lupri)
       else
            call arh_crop_solver(decomp, arh, debug, wrk, antisymmetric, x, fifoqueue)
       end if
    else
       call lsquit('Something wrong, no ARH solver chosen',decomp%lupri)
    endif

     if (arh%xnorm < 0.2E0_realk .and. SCF_iteration > 4 ) then
        arh%set_local = .true.
        if (arh%cfg_2nd_order_local .and. .not. arh%set_do_2nd_order) then
           write(arh%lupri,*) 'Local region - switching to second order optimization!'
           arh%set_do_2nd_order = .true.
           arh%cfg_arh_truncate = .false.
        endif
     endif
     
     if (.not. davidCFG%arh_davidson) then
           arh%xnorm=dsqrt(mat_sqnorm2(X))
          !Now we construct the new density from X:
          CALL LSTIMER('START ',t1,t2,decomp%lupri)
          call oao_density_param(x,decomp%DU,wrk) !wrk = new D(X)
          call oao_purify(wrk,wrk2) !wrk2 = new purified density in oao basis
          CALL LSTIMER('NEW D ',t1,t2,arh%lupri)
      endif



     if (arh%set_arhterms) then
        !Determine part of the change in density which can be expanded in the subspace:
        arh%D_para_D_tot = 0.0E0_realk !So we don't get undefined output if ndens = 0
        if (ndens > 0) then
           call mem_alloc(weights,ndens)

           call arh_get_weights(arh,wrk2,x,fifoqueue,weights)

           call mat_zero(x) !x is now part of the change in density which can be expanded in the subspace
           call get_from_modFIFO(fifoqueue, 0, dummyp, D0)
           do i = 1, fifoqueue%offset
              call get_from_modFIFO(fifoqueue, i, dummyp, Di) 
              call MAT_ADD(1.0E0_realk, Di, -1.0E0_realk, D0, wrk)
              call mat_daxpy(weights(i),wrk,x)
           enddo

           !Actual change in density is Dchol - DU
           !Part of change which can be expanded is res
           call MAT_ADD(1.0E0_realk, wrk2, -1.0E0_realk, D0, wrk)  
           actual_change_norm = SQRT(mat_sqnorm2(wrk))
           expanded_change_norm = SQRT(mat_sqnorm2(x))

           arh%D_para_D_tot = expanded_change_norm/actual_change_norm
           call mem_dealloc(weights)
        endif
     endif

     !Stinne 24/1-07: Why tranform to AO basis? We already calculated the corresponding Fock matrix...
     !30/1-07: Because otherwise a wrong density is written in dens.restart!!
     if (decomp%cfg_do_in_oao) then
        call mat_assign(Dnew,wrk2)
     else
        call x_from_oao_basis(decomp,wrk2, Dnew) 
     endif

    if (associated(arh%fifometric)) then
       deallocate(arh%fifometric)
       nullify(arh%fifometric)
    endif
    if (associated(arh%inv_fifometric)) then
       deallocate(arh%inv_fifometric)
       nullify(arh%inv_fifometric)
    endif
    if (associated(arh%fifoM)) then
       deallocate(arh%fifoM)
       nullify(arh%fifoM)
    endif
    call mat_free(x)
    call mat_free(wrk)
    call mat_free(wrk2)
   end subroutine arh_get_density


subroutine arh_davidson_solver(CFG,arh,decomp,wrk,X,SCF_iteration,H1,wrk2,ls)
implicit none
type(RedSpaceItem)         :: CFG
type(lsitem) ::ls
type(SolverItem)           :: arh
type(matrix),intent(inout) :: X,wrk,wrk2
!> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
type(decompItem),intent(inout) :: decomp
!> One electron Hamiltonian
type(matrix),intent(in) :: H1
!> Current SCF iteration
integer, intent(in)     :: SCF_iteration
integer, parameter      :: ndensmat=10 
real(realk)             :: t1, t2, gradnorm
real(realk)    :: stepfactor
type(matrix)   :: densmat(ndensmat),ax,G
real(realk)    :: energy_array(ndensmat),scaling,saveTHRESHOLD
real(realk)    :: arh_Ediff_actual,Ecurrent,Eold
integer        :: indx(1),i,j,counter,idebug
logical        :: saveLSDASCREEN,saveSaveF0andD0
logical        :: TightLineSearchEsitmatesThresholds
type(matrix),pointer :: debugFockMatArray(:)
real(realk),pointer  :: EDFTY(:),debugEtotalarray(:)
real(realk)          :: alpha(ndensmat),xnorm

call mat_init(G,wrk%nrow,wrk%ncol)
call mat_copy(-1.0_realk,wrk,G)  !Scale gradient since davidson solver is made for (H-mu)X= -G
call mat_init(ax,x%nrow,x%ncol)
counter = 0
WHILELOOP: do j=1,10
  if (CFG%stepsize < 0.0001_realk) then
      write(ls%lupri,*) ' **************************************************************'
      write(ls%lupri,*) ' * Too many rejections in SCF optimization                    *'
      write(ls%lupri,*) ' * Consider using a looser SCF convergence threshold          *'
      write(ls%lupri,*) ' * or try too tighten integral accuracy for arh line search.  *'
      write(ls%lupri,*) ' **************************************************************'
      call lsquit('Too many rejections in SCF opt', ls%lupri)
   end if
Linesearch: if (CFG%arh_linesearch) then
   counter = counter+1
   !Scale gradient since davdison solver is made for (H-mu)X= -G
   CALL LSTIMER('START ',t1,t2,decomp%lupri)
   call solver(CFG,G,X)
   CALL LSTIMER('DAVID. SOLVER ',t1,t2,arh%lupri)

    !Now we construct the new density from X:
  if (CFG%arh_gradnorm < 1E-2_realk) stepfactor=0.25_realk
  if (CFG%arh_gradnorm .ge. 1E-2_realk) stepfactor=0.5_realk
  if (CFG%arh_gradnorm .ge. 1E-1_realk) stepfactor=0.75_realk
  if (CFG%arh_gradnorm .ge. 1E0_realk) stepfactor=1.0_realk
  scaling = 1.0_realk-stepfactor
  do i=1,ndensmat
      scaling = scaling + stepfactor 
      if (abs(scaling) <1E-5_realk) scaling=scaling+stepfactor 
      alpha(i)=scaling
  end do
  
  CALL LSTIMER('START ',t1,t2,decomp%lupri)
  do i =1, ndensmat
     call mat_init(densmat(i),wrk%nrow,wrk%ncol)
     call mat_copy(alpha(i),x,aX)
     call oao_density_param(aX,decomp%DU,densmat(i)) !wrk = new D(X)
     call oao_purify(densmat(i),wrk2) !wrk2 = new purified density in oao basis
     call x_from_oao_basis(decomp,wrk2, densmat(i)) 
  end do

  call linesearch_thresholds(CFG,ls,SCF_iteration,TightLineSearchEsitmatesThresholds)

  call di_SCF_EnergyCont(densmat,H1,energy_array,ndensmat,&
       & CFG%lupri,CFG%lupri,ls,CFG%LSmodthresh)
  CALL LSTIMER('NEW D ',t1,t2,CFG%lupri)
  CFG%MaxLineSearchEnergyDiff = 10000.0E0_realk
  indx = MINLOC(energy_array) 
  
  if (indx(1) ==  1)THEN
     CFG%MaxLineSearchEnergyDiff=ABS(energy_array(indx(1))-energy_array(indx(1)+1))               
  elseif (indx(1) <  ndensmat)THEN
     CFG%MaxLineSearchEnergyDiff=ABS(energy_array(indx(1))-energy_array(indx(1)+1))
     CFG%MaxLineSearchEnergyDiff=MIN(CFG%MaxLineSearchEnergyDiff,&
     &ABS(energy_array(indx(1))-energy_array(indx(1)-1)))
  elseif (indx(1) == ndensmat)then
     CFG%MaxLineSearchEnergyDiff=ABS(energy_array(indx(1))-energy_array(indx(1)-1))
  endif

  ! If energy obtained is larger than old energy, reduce stepsize and call solver
  if (energy_array(indx(1)) > arh%old_energy) then
     CFG%stepsize = (2.0_realk/3.0_realk)*CFG%stepsize
     do i=1,ndensmat
        call mat_free(densmat(i))
     end do
     cycle
  end if

  write(arh%lupri,*)
  write(arh%lupri,'(a)') '***** LINESEARCH INFORMATION *****'
  write(arh%lupri,*)
  if (CFG%arh_davidson_debug) write(arh%lupri,'(a,ES20.11)') 'Old energy',arh%old_energy
  xnorm=dsqrt(mat_dotproduct(x,x))
  if (CFG%arh_davidson_debug) then
    do i=1,ndensmat
       write(arh%lupri,'(a,ES20.11,a,f10.4)') 'Energy', energy_array(i), &
       &' Tot. stepsize:', alpha(i)*xnorm 
    end do
  endif
  call mat_scal(alpha(indx(1)),X)
  arh%xnorm = sqrt(mat_sqnorm2(x))
  write(arh%lupri,'(a,i4,a,ES14.7)') 'ARHLS: Minimum function value for ', indx(1),&
  &' using step size ', arh%xnorm
  CFG%arh_linesE=energy_array(indx(1))
  if (CFG%arh_davidson_debug) write(arh%lupri,'(a,ES20.11)') &
  &'ARHLS: Energy From Linesearch   ',CFG%arh_linesE 
  if (CFG%arh_davidson_debug) write(arh%lupri,'(a,ES20.11)') &
  &'ARHLS: Maximum Energy difference', CFG%MaxLineSearchEnergyDiff
  write(arh%lupri,*)
  write(arh%lupri,'(a)') '***** END  LINESEARCH INFORMATION *****'
  write(arh%lupri,*)

  call x_to_oao_basis(decomp,densmat(indx(1)),wrk2)
  do i=1,ndensmat
     call mat_free(densmat(i))
  end do
  EXIT
else
  CALL LSTIMER('START ',t1,t2,decomp%lupri)
  call solver(CFG,G,X)
  CALL LSTIMER('DAVID. SOLVER ',t1,t2,arh%lupri)
  arh%xnorm=dsqrt(mat_sqnorm2(X))
  !Now we construct the new density from X:
  CALL LSTIMER('START ',t1,t2,decomp%lupri)
  call oao_density_param(x,decomp%DU,wrk) !wrk = new D(X)
  call oao_purify(wrk,wrk2) !wrk2 = new purified density in oao basis
  CALL LSTIMER('NEW D ',t1,t2,arh%lupri)
  EXIT
end if Linesearch 
end do WHILELOOP

!print info 
call print_info_davidson(CFG,arh%xnorm,SCF_iteration)
call mat_free(G)
call mat_free(ax)



end subroutine arh_davidson_solver



subroutine linesearch_thresholds(CFG,ls,SCF_iteration,TightLineSearchEsitmatesThresholds)
implicit none
type(RedSpaceItem)  :: CFG
type(lsitem)        :: ls
integer, intent(in) :: SCF_iteration
logical             :: TightLineSearchEsitmatesThresholds


  IF(SCF_iteration.EQ.1)THEN
     CFG%EnergyDiffset = .FALSE.
     CFG%LSmodthresh = 1000.0E0_realk     
     ls%setting%scheme%LSDASCREEN_THRLOG = 0
  ENDIF
  IF(SCF_iteration.EQ.2)THEN
     !CFG%LSmodthresh = 1000 means that the thresholds for energy eval is 
     !a factor 1000 more loose than the Fock matrix thresholds
     CFG%LSmodthresh = 1000.0E0_realk     
  ENDIF
  IF(SCF_iteration.EQ.3)THEN
     CFG%LSmodthresh = 100.0E0_realk
  ENDIF
  IF(SCF_iteration.EQ.4)THEN
     CFG%LSmodthresh = 10.0E0_realk
  ENDIF
  IF(SCF_iteration.EQ.5)THEN
     CFG%LSmodthresh = 1.0E0_realk
     ls%setting%scheme%LSDASCREEN_THRLOG = 1
  ENDIF
  IF(SCF_iteration.EQ.6)THEN
     ls%setting%scheme%LSDASCREEN_THRLOG = 2
  ENDIF
  TightLineSearchEsitmatesThresholds = .FALSE.
  IF(.NOT.SaveF0andD0)THEN
     TightLineSearchEsitmatesThresholds = .TRUE.
  ENDIF
  IF(CFG%EnergyDiffset)THEN
     !we tighten if one of the 3 things happen 
     IF(CFG%MaxLineSearchEnergyDiff/CFG%ActualEnergyDiff .LT. 100)THEN
        !the difference between the linesearch estimate and actual energy less than a factor 100
        TightLineSearchEsitmatesThresholds = .TRUE.
!        WRITE(CFG%lupri,*)'We tighten LS thresholds due to bad ration'
        WRITE(CFG%lupri,*)CFG%MaxLineSearchEnergyDiff/CFG%ActualEnergyDiff,' .LT. 100'
     ENDIF
     IF(ls%SETTING%SCHEME%THRESHOLD*CFG%LSmodthresh.GT.0.1E0_realk*CFG%arh%xnorm*CFG%arh%xnorm)THEN
        TightLineSearchEsitmatesThresholds = .TRUE.
!        WRITE(CFG%lupri,*)'We tighten LS thresholds due to small xnorm'
     ENDIF
     IF(CFG%MaxLineSearchEnergyDiff .LT. ABS(CFG%arh%old_energy-CFG%arh_linesE)) then
        TightLineSearchEsitmatesThresholds = .TRUE.
!        WRITE(CFG%lupri,*)'We tighten LS thresholds due possible wrong linesearch'
     ENDIF
  ENDIF
  IF(TightLineSearchEsitmatesThresholds)THEN
     ls%setting%scheme%LSDASCREEN_THRLOG = MIN(ls%setting%scheme%LSDASCREEN_THRLOG+1,2)
     CFG%LSmodthresh = MAX(0.1E0_realk*CFG%LSmodthresh,1.0E0_realk)
     !CFG%LSmodthresh = 1 means that the Fock matrix and energy have same thresholds
  ENDIF


end subroutine linesearch_thresholds





!> \brief Reduced space solver used for solving Augmented Roothaan-Hall equations
!> \author S. Host
!> \date 2007
!>
!> The solver is based on the Conjugate Residual OPtimal vectors (CROP) scheme: \n
!>  M. Ziolkowski, V. Weijo, P. Jorgensen et al. JCP 128, 204105 \n
!> There are two options: 
!> - Run with the full subspace and keep vectors on disk (cfg_arh_crop=.true.) - this is invoked with .ARH FULL under
!> the *LINSCA section. 
!> - Run with a truncated number of vectors in memory (cfg_arh_truncate=.true.) 
!> (this is standard when using .ARH). Number of vectors to be kept may be set with .MICROVECS. The minimum (and default) is 2, corresponding
!> to two previous trial- and sigma-vectors. Plus the current ones, this means that 3 trial- and sigmas in total are stored. 
!>
!FIXME: Make it an option whether vectors should be kept on disk or in memory.
   subroutine arh_crop_solver(decomp, arh, debug, Grad, symm, x, fifoqueue)
   use direct_dens_util
   use levelshift_mod
      implicit none
      !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
      type(decompItem),intent(in)    :: decomp
      !> Contains solver info (ARH/TrFD)
      type(solverItem),intent(inout) :: arh
      !> If requested, some debug output is put in here to be printed later
      type(debugItem),intent(inout)  :: debug
      !> SCF gradient in orthonormal AO (OAO) basis
      type(Matrix), intent(in)    :: Grad
      !> If symm = 1, trial vector is symmetric. If symm = 2, trial vector is antisymmetric.
      integer, intent(in)         :: symm
      !> Contains Fock/KS and density matrices from previous SCF iterations
      TYPE(modFIFO),intent(inout) :: fifoqueue
      !> Output. The X that minimizes E = TrFD(X) + ARH terms(X)
      type(Matrix),intent(inout)  :: x
      type(Matrix)                :: scrmat, res, resP, xsave
      integer                     :: i, j, k, l, rowdim, coldim, max_it, redspacedim_save, xsave_lu, ndens, dampdim, it
      integer                     :: lub, lusigma
      real(realk)                 :: err, errsave, t1, t2, mumax, thresh
      logical                     :: fileexists, fileopened, done, do_LS,OnMaster
      real(realk)                 :: mu
!Truncate subspace:
      TYPE(modFIFO),target        :: vectorsubspace
      integer                     :: maxvec
!New levelshift
      type(Matrix)                :: xF, sigmaF
      real(realk)                 :: cutoff
!Queue on disk:
!      integer                     :: queue_lu
!      logical                     :: exppoint
      type(lshiftItem)             :: lshift
!Levelshift by homo-lumo gap
       real(realk)                 :: hlgap, tstart, tend
       OnMaster=.TRUE.
   if (arh%set_arhterms) then
      ndens = fifoqueue%offset
      write(arh%lupri,*) 'Number of densities in queue:', ndens
   endif
   done = .false.

   !The SCF energy from the previous iteration is found in the queue:
   rowdim = Grad%nrow
   coldim = Grad%ncol
   !matdim = Grad%nrow
   max_it = 200  !Stinne's first guess as to how many iterations are needed
!   max_it = 60 Brano???

   arh%set_optxelm = .false.
   call set_levelshift_config(arh,lshift)

   !Threshold for convergence:
   thresh = arh%cfg_micro_thresh*sqrt(mat_sqnorm2(Grad))
   if (sqrt(mat_sqnorm2(Grad)) < 1.0E-3_realk) thresh=thresh*0.1_realk
   !if (sqrt(mat_sqnorm2(Grad)) < 1.0E-5_realk) thresh=thresh*0.1_realk

   if (arh%cfg_arh_truncate) then
      maxvec = arh%cfg_arh_microvecs
   else
      maxvec = max_it
   endif
   if (arh%cfg_arh_newdamp) then
      call mat_init(xF,rowdim,coldim)
      call mat_init(sigmaF,rowdim,coldim)
      call mem_alloc(arh%Ared,maxvec+1,maxvec+1)
      call mem_alloc(arh%Gred,maxvec+1)
      call mem_alloc(arh%Sred,maxvec+1,maxvec+1)
      call mem_alloc(arh%CROPmat,maxvec,maxvec)
      dampdim = maxvec+1
   else
      call mem_alloc(arh%Ared,maxvec,maxvec)
      call mem_alloc(arh%Gred,maxvec)
      call mem_alloc(arh%Sred,maxvec,maxvec)
      call mem_alloc(arh%CROPmat,maxvec,maxvec)
      dampdim = maxvec
   endif
      
   call mat_init(scrmat,rowdim,coldim)
   call mat_init(res,rowdim,coldim)
   call mat_init(resP,rowdim,coldim)

   mu = 0.0E0_realk 
   if (arh%cfg_fixed_shift) mu = -arh%cfg_fixed_shift_param

   CALL LSTIMER('START ',t1,t2,arh%lupri)
   if (arh%info_lineq) then 
      if (arh%cfg_arh_truncate) then
         WRITE(arh%LUPRI, "('CROP scheme, truncated subspace (CORE) - max. no. of vectors is', i3)") arh%cfg_arh_microvecs
      else if (arh%cfg_arh_crop) then
         WRITE(arh%LUPRI, "('CROP scheme with full subspace (DISK)')")
      endif
      if (.not. arh%cfg_arh_crop_safe) then
         WRITE(arh%LUPRI, "('WARNING: Double preconditioning has been removed!')")
      endif
   endif
   arh%trustradius_decreased = .false.

   !The maximum number of rejections is currently set to six
   ! - if there are that many rejections, the calculation is definitely
   ! unhealthy and should be stopped. This is almost always caused by lack of
   ! integral accuracy.
   if (arh%Nrejections > 6) then
     WRITE(arh%LUPRI, "('Too many rejections - probably related to lack of integral accuracy!')")
     CALL lsQUIT('Too many rejections - probably related to lack of integral accuracy!',decomp%lupri)
   endif

   arh%Ared = 0.0E0_realk ; arh%Gred = 0.0E0_realk
   arh%Sred = 0.0E0_realk ; arh%CROPmat = 0.0E0_realk

   if (arh%cfg_arh_truncate) then
      call modfifo_init(vectorsubspace,maxvec+1,rowdim,coldim,arh%cfg_arh_disk_micro) 
   else
      lusigma = -1 ; lub = -1
      CALL LSOPEN(lusigma,'sigmavecs','unknown','UNFORMATTED')
      CALL LSOPEN(lub,'bvecs','unknown','UNFORMATTED')
   endif
   !gradnorm = SQRT(mat_sqnorm2(Grad))

   ! Use starting guess, if available
   if (arh%starting_guess_defined) then
      !write(arh%lupri,*) 'starting guess:'
      !call mat_print(arh%starting_guess, 1, arh%starting_guess%nrow, 1, arh%starting_guess%ncol, arh%lupri)
      call mat_copy(1.0E0_realk, arh%starting_guess, x)
   else
      !Save first trial b vector = gradient:
      !call mat_copy(-1.0E0_realk, Grad, scrmat)
      call project_oao_basis(decomp, Grad, symm, x)
      call mat_scal(-1.0E0_realk,x)
   endif

   !First linear transformation:
   call arh_lintrans(arh,decomp,x,symm,0.0E0_realk,scrmat,fifoqueue)

   if (arh%cfg_arh_truncate) then
      call add_to_modFIFO(vectorsubspace, x, scrmat, 0.0E0_realk, .false.)
   else
      call mat_write_to_disk(lub,x,OnMaster)
      call mat_write_to_disk(lusigma,scrmat,OnMaster)
   endif
   arh%Gred(1) = mat_dotproduct(x,Grad)
   arh%Sred(1,1) = mat_dotproduct(x,x)
   arh%Ared(1,1) = mat_dotproduct(x,scrmat)

   !We need to call levelshift here, otherwise mu is wrong when constructing 1st residual
   if (arh%lshift_by_hlgap) then
      arh%cfg_fixed_shift = .false.
      arh%set_local = .true.
      arh%nits_check = 3
      arh%error_decrease = 0.5E0_realk 
      CALL LSTIMER('START ',tstart,tend,arh%lupri)
      call dd_debug_homolumo(decomp,hlgap)
      CALL LSTIMER('HOMO-LUMO gap',tstart,tend,arh%lupri)
      if (hlgap < 0.0E0_realk) then
         mu = 1.10E0_realk*hlgap
      else
         mu = 0.0E0_realk
      endif
      if (.not. arh%step_accepted) then
         !We tried damping from this HL gap before without succes, so now we increase it:
         if (abs(arh%current_mu) > 1.0E-2_realk) then
            mu = arh%current_mu*2.0E0_realk
         else
            if (abs(arh%maxelm) > 0.1) then
               mu = -0.5E0_realk
            else
               mu = -0.1E0_realk
            endif
         endif
         write(arh%lupri,*) 'Previous step was rejected, and mu is increased to:', mu
      endif
   else
      call main_levelshift(lshift,arh%Ared,arh%Sred,arh%Gred,maxvec,lub,1,rowdim,coldim,mu,vectorsubspace=vectorsubspace)
   endif
   write(arh%lupri,*) 'First mu:', mu 
   call mat_add(-1E0_realk,Grad,-1E0_realk, scrmat, res) !res = 1st residual
   call mat_daxpy(mu,x,res)

   !Preconditioning:
   call arh_precond(arh,decomp,res,symm,mu,resP)

   arh%CROPmat(1,1) = mat_dotproduct(res,resP)
   if (arh%info_crop) then
      write(arh%lupri,*) 'First elements in reduced space:'
      write(arh%lupri,*) 'Sred(1,1) =', arh%Sred(1,1) 
      write(arh%lupri,*) 'Ared(1,1) =', arh%Ared(1,1) 
      write(arh%lupri,*) 'CROPmat(1,1) =', arh%CROPmat(1,1) 
   endif

   !if (arh%debug_diag_redspace) then
   !   write(arh%lupri,*) 'First element in reduced space:', arh%Ared(2,2)
   !   debug%final_redspace_eival = arh%Ared(2,2)
   !endif

   j = 0 ; l = 0 ; i = 0
   if (arh%lshift_by_hlgap) then
      k = 1
   else
      k = 5
   endif
   do i = 1, max_it
      if (i == max_it) then
         if (arh%cfg_orbspread) then
             write (arh%lupri,*) 'Convergence NOT obtained, but continue anyway'
             write (arh%lupri,*) '(using solver for localization of orbitals)'
             exit
         else if (arh%set_optxelm) then
            if (arh%INFO_LINEQ) write (arh%lupri,*) 'MAXIT reached: Optimization of xmax failed, I will use old x'
            call mat_init(xsave,rowdim,coldim)
            rewind(xsave_lu)
            call mat_read_from_disk(xsave_lu,xsave,OnMaster)
            call mat_max_elm(xsave, arh%maxelm)
            arh%final_redspacedim = redspacedim_save
            mu = mumax
            call mat_assign(x,xsave)
            call mat_free(xsave)
            exit
         else
            WRITE(arh%LUPRI,'(/A)') &
            &     'ARH: Linear equations not converged'
            CALL lsQUIT('ARH: Linear equations not converged',decomp%lupri)
         endif
      endif
      if (i > 1) call main_levelshift(lshift,arh%Ared,arh%Sred,arh%Gred,dampdim,lub,&
                             & i,rowdim,coldim,mu,xF,sigmaF,vectorsubspace)
      !call arh_crop_x_and_res(fifoqueue,Grad,symm,mu,b_current,scrmat,resP,res)
      call arh_crop_x_and_res(arh,decomp,fifoqueue,Grad,symm,mu,x,scrmat,resP,res)
      !res = resP Stinne bugfix February 2010
      !x = b_current
      !Convergence?
      err = SQRT(mat_sqnorm2(res))
      it = i !"A do-variable within a DO body shall not appear in a variable definition context" pointed out by Thomas
      call arh_test_convergence(arh,err,errsave,x,mu,mumax,thresh,it,j,k,xsave_lu,redspacedim_save,maxvec,do_LS,done)
      call set_levelshift_config(arh,lshift)
      if (do_LS) call main_levelshift(lshift,arh%Ared,arh%Sred,arh%Gred,dampdim,lub,i,rowdim,coldim,mu,xF,sigmaF,vectorsubspace)
      if (done) exit
      !call arh_crop_intermed_sub(vectorsubspace,i,symm,Grad,mu,resP,b_current,scrmat,xF,sigmaF)
      call arh_crop_intermed_sub(arh,decomp,lub,lusigma,vectorsubspace,i,symm,Grad,mu,res,x,scrmat,xF,sigmaF) !resP -> res, Stinne bugfix February 2010
      !Save current on disc:
      if (arh%cfg_arh_truncate) then
         call add_to_modFIFO(vectorsubspace,x,scrmat,0.0E0_realk,.false.)
      else
         call mat_write_to_disk(lub,x,OnMaster)
         call mat_write_to_disk(lusigma,scrmat,OnMaster)
      endif
      !call arh_crop_setup_redsp(vectorsubspace,i,symm,Grad,mu,b_current,scrmat,resP,xF,sigmaF)
      call arh_crop_setup_redsp(arh,decomp,lub,lusigma,vectorsubspace,i,symm,Grad,mu,x,scrmat,resP,res,xF,sigmaF)
   enddo

   arh%current_mu = mu

   !if (arh%debug_diag_redspace .and. .not. arh%cfg_arh_crop) then
   !   write (arh%lupri,'("Final lowest redspace eigenvalue: ",F12.6, "      ")') debug%final_redspace_eival
   !endif

   CALL LSTIMER('CROP solver',t1,t2,arh%lupri)

   if (arh%cfg_arh_crop) call mat_scal(-1.0E0_realk,x) !HACK!!!!

   call arh_get_TR_denom(arh, Grad, x, scrmat, decomp%cfg_unres, mu)

   if (arh%cfg_arh_truncate) then
      call modfifo_free(vectorsubspace)
   else
      CALL LSCLOSE(lusigma,'DELETE')
      CALL LSCLOSE(lub,'DELETE')
   endif

   INQUIRE(file='xsave',EXIST=fileexists,OPENED=fileopened)
   if (fileexists.AND.fileopened) call LSCLOSE(xsave_lu,'DELETE')

   if (arh%cfg_arh_newdamp) then
      call mat_free(xF)
      call mat_free(sigmaF)
   endif
   call mem_dealloc(arh%Ared)
   call mem_dealloc(arh%Gred)
   call mem_dealloc(arh%Sred)
   call mem_dealloc(arh%CROPmat)
   call mat_free(scrmat)
   call mat_free(res)
   call mat_free(resP)
   end subroutine arh_crop_solver

!> \brief Set the parameters that must be passed to level shift module
!> \author S. Host
!> \date March 2010
   subroutine set_levelshift_config(arh,lshift)
   use levelshift_mod
   implicit none
      !> Contains solver info (ARH/TrFD)
      type(solverItem),intent(in) :: arh
      !> Contains info used for level shifting
      type(lshiftItem),intent(inout) :: lshift

      lshift%lupri              = arh%lupri
      lshift%fixed_shift        = arh%cfg_fixed_shift       
      lshift%fixed_shift_param  = arh%cfg_fixed_shift_param 
      lshift%arh_crop           = arh%cfg_arh_crop          
      lshift%arh_truncate       = arh%cfg_arh_truncate      
      lshift%arh_newdamp        = arh%cfg_arh_newdamp       
      lshift%min_lshift         = arh%cfg_min_lshift        
      lshift%max_element        = arh%set_max_element
      lshift%max_step           = arh%set_max_step
      lshift%optxelm            = arh%set_optxelm
      lshift%info_levelshift    = arh%info_levelshift
      lshift%current_mu         = arh%current_mu
      lshift%lshift_by_hlgap    = arh%lshift_by_hlgap
   end subroutine set_levelshift_config

! ------------ Debug section! ---------------!

!> \brief Print debug info.
!> \author S. Host
!> \date 2007
   subroutine arh_debug_print(arh)
   use arhDensity
   implicit none
        !> Contains solver info (ARH/TrFD)
        type(solverItem),intent(in) :: arh

   if (.not. arh%step_accepted) then
      write (arh%lupri, "(i5, F36.6, F12.6, '                ///')") &
           & arh%scfit, arh%final_redspace_eival, arh%arheival
   else          
      write (arh%lupri, "(i5, F12.6, F12.6, F12.6, F12.6, F12.6, '    ///')") &
           & arh%scfit, arh%diag_hlgap, arh%iter_hlgap, arh%final_redspace_eival, &
           & arh%arheival, arh%heseival
   endif
   end subroutine arh_debug_print


  subroutine print_info_davidson(davidCFG,step,it)
  type(RedSpaceItem) :: davidCFG
  real(realk)  :: step 
  integer :: it
  

  write(davidCFG%lupri,'(f10.4,f12.1,5X,f10.4,f11.2,f8.1,f7.1,f14.1,6X,i5,a)') davidCFG%stepsize, 0.0,step, davidCFG%mu,&
  &0.0,0.0,0.0,it+1,'  %%'
 

  end subroutine print_info_davidson
  
  subroutine debug_arh_LineSearch(CFG,arh,ls,H1,densmat)
    implicit none
    type(RedSpaceItem)         :: CFG
    type(lsitem) ::ls
    type(SolverItem)           :: arh
    !> One electron Hamiltonian
    type(matrix),intent(in) :: H1
    !> Current SCF iteration
    integer, parameter      :: ndensmat=10 
    type(matrix)   :: densmat(ndensmat)
    real(realk)    :: energy_array(ndensmat),scaling,saveTHRESHOLD
    real(realk)    :: arh_Ediff_actual,modthresh
    integer        :: indx(1),i,j,counter,idebug
    logical        :: saveLSDASCREEN,saveSaveF0andD0 
    type(matrix),pointer :: debugFockMatArray(:)
    real(realk),pointer  :: EDFTY(:),debugEtotalarray(:)
    modthresh = 1.0E0_realk
    write(arh%lupri,*)
    write(arh%lupri,'(a)') '***** DEBUG LINESEARCH INFORMATION *****'
    write(arh%lupri,*)
    call mem_alloc(debugFockMatArray,ndensmat)
    call mem_alloc(EDFTY,ndensmat)
    call mem_alloc(debugEtotalarray,ndensmat)
    saveTHRESHOLD = ls%SETTING%SCHEME%THRESHOLD
    do i=1,ndensmat
       call mat_init(debugFockMatArray(i),densmat(1)%nrow,densmat(1)%ncol)
    enddo
    do idebug=0,2
       ls%SETTING%SCHEME%THRESHOLD = saveTHRESHOLD*10.0E0_realk**(-idebug) 
       WRITE(arh%lupri,*)'ls%SETTING%SCHEME%THRESHOLD',ls%SETTING%SCHEME%THRESHOLD
       call II_get_Fock_mat(arh%lupri,6,ls%setting,Densmat,.TRUE.,debugFockmatArray,ndensmat,.FALSE.)
       do i=1,ndensmat
          debugEtotalarray(i) = fockenergy_f(debugFockmatArray(i),Densmat(i),H1,&
               & ls%input%dalton%unres,ls%input%potnuc,arh%lupri)
       enddo
       IF(ls%setting%do_dft) THEN
          call II_get_xc_fock_mat(arh%lupri,6,ls%setting,Densmat(1)%nrow,Densmat,&
               & debugFockmatArray,Edfty,ndensmat)
          do i=1,ndensmat
             debugEtotalarray(i) = debugEtotalarray(i) + Edfty(i)
          enddo
       ENDIF
       do i=1,ndensmat
          write(arh%lupri,*) 'Actual Fock Matrix Energy NF:', debugEtotalarray(i)
       enddo
       do i=1,ndensmat
          write(arh%lupri,'(a,ES20.11)') 'Actual Fock Matrix Energy', debugEtotalarray(i)
       enddo

       !=================================================================================
       write(arh%lupri,'(a,ES20.11)') 'Linesearch With full approximations'
       !=================================================================================
       call di_SCF_EnergyCont(densmat,H1,energy_array,ndensmat,&
            & arh%lupri,arh%lupri,ls,modthresh)
       CFG%MaxLineSearchEnergyDiff = 10000.0E0_realk
       indx = MINLOC(energy_array) 
       if (indx(1) ==  1)THEN
          CFG%MaxLineSearchEnergyDiff=ABS(energy_array(indx(1))-energy_array(indx(1)+1))               
          arh_Ediff_actual=ABS(debugEtotalarray(indx(1))-debugEtotalarray(indx(1)+1))               
       elseif (indx(1) <  ndensmat)THEN
          CFG%MaxLineSearchEnergyDiff=ABS(energy_array(indx(1))-energy_array(indx(1)+1))
          CFG%MaxLineSearchEnergyDiff=MIN(CFG%MaxLineSearchEnergyDiff,ABS(energy_array(indx(1))-energy_array(indx(1)-1)))
          arh_Ediff_actual = ABS(debugEtotalarray(indx(1))-debugEtotalarray(indx(1)+1))
          arh_Ediff_actual = MIN(arh_Ediff_actual,ABS(debugEtotalarray(indx(1))-debugEtotalarray(indx(1)-1)))
       elseif (indx(1) == ndensmat)then
          CFG%MaxLineSearchEnergyDiff=ABS(energy_array(indx(1))-energy_array(indx(1)-1))
          arh_Ediff_actual=ABS(debugEtotalarray(indx(1))-debugEtotalarray(indx(1)-1))               
       endif
       write(arh%lupri,'(a,ES20.11)') 'Old energy',arh%old_energy
       !if Old energy + change is not equal to Energy_array it is the first time the
       !linesearch is performed and it would be cunfusing to have the energy change printed
       do i=1,ndensmat
          write(arh%lupri,'(a,ES20.11)') 'Energy', energy_array(i)
       end do
       do i=1,ndensmat
          write(arh%lupri,*) 'Energy NF:', energy_array(i)
       end do
       write(arh%lupri,'(a,i4)') 'Minimum function value for ', indx(1)
       CFG%arh_linesE=energy_array(indx(1))
       write(arh%lupri,'(a,ES20.11)') 'Energy From Linesearch   ', CFG%arh_linesE 
       write(arh%lupri,'(a,ES20.11)') 'Maximum Energy difference', CFG%MaxLineSearchEnergyDiff

       if (CFG%MaxLineSearchEnergyDiff > ABS(debugEtotalarray(indx(1))-CFG%arh_linesE)) then
          WRITE(arh%lupri,'(A)') 'default no error in linesearch'
          WRITE(arh%lupri,'(A,ES22.13)')'actual Etotal from Fock Matrix  ',debugEtotalarray(indx(1))
          WRITE(arh%lupri,'(A,ES22.13)')'linesearch Energy Estimate      ',CFG%arh_linesE
          WRITE(arh%lupri,'(A,ES22.13)')'Energy difference between estimate and actual Energy',&
               & ABS(debugEtotalarray(indx(1))-CFG%arh_linesE)
          WRITE(arh%lupri,'(A,ES22.13)')'Maximum difference in linesearch energies',&
               & CFG%MaxLineSearchEnergyDiff
          WRITE(arh%lupri,'(A,ES22.13)')'Maximum difference in actual energies',&
               & arh_Ediff_actual
!          IF(ABS(CFG%MaxLineSearchEnergyDiff) < 1.0E-10_realk)THEN
!             !turn off incremental build in linesearch, to improve accuracy
!             !              IF(config%opt%cfg_saveF0andD0)THEN
!             WRITE(arh%lupri,*)'Turn off incremental linesearch energy estimates to improve accuracy if cfg_saveF0andD0'
!             !              ENDIF
!          ENDIF
          IF(CFG%MaxLineSearchEnergyDiff < & 
               & 5.0E+2_realk*LS%SETTING%SCHEME%THRESHOLD*LS%SETTING%SCHEME%J_THR.OR.&
               & CFG%MaxLineSearchEnergyDiff < &
               & 1.5E+2_realk*ABS(debugEtotalarray(indx(1))-CFG%arh_linesE))THEN
             !improve accuracy on Fock matrix as the energy differences in the linesearch are small
             !compared to Fock matrix accuracy
             WRITE(arh%lupri,'(A)')'The SCF convergence have reached a point where we should'
             WRITE(arh%lupri,'(A)')'increase the accuracy of the Fock matrix build or increase the'
             WRITE(arh%lupri,'(A)')'accuracy of the linesearch energy estimates. We do both!'
          ENDIF
       endif

       !=================================================================================
       write(arh%lupri,'(a,ES20.11)') 'Linesearch With no approximations'
       !=================================================================================
       saveSaveF0andD0 = SaveF0andD0 
       saveLSDASCREEN = ls%setting%scheme%LSDASCREEN
       SaveF0andD0 = .FALSE.
       ls%setting%scheme%LSDASCREEN = .FALSE.
       call di_SCF_EnergyCont(densmat,H1,energy_array,ndensmat,&
            & arh%lupri,arh%lupri,ls,modthresh)
       CFG%MaxLineSearchEnergyDiff = 10000.0E0_realk
       indx = MINLOC(energy_array) 
       if (indx(1) ==  1)THEN
          CFG%MaxLineSearchEnergyDiff=ABS(energy_array(indx(1))-energy_array(indx(1)+1))               
          arh_Ediff_actual=ABS(debugEtotalarray(indx(1))-debugEtotalarray(indx(1)+1))               
       elseif (indx(1) <  ndensmat)THEN
          CFG%MaxLineSearchEnergyDiff=ABS(energy_array(indx(1))-energy_array(indx(1)+1))
          CFG%MaxLineSearchEnergyDiff=MIN(CFG%MaxLineSearchEnergyDiff,ABS(energy_array(indx(1))-energy_array(indx(1)-1)))
          arh_Ediff_actual = ABS(debugEtotalarray(indx(1))-debugEtotalarray(indx(1)+1))
          arh_Ediff_actual = MIN(arh_Ediff_actual,ABS(debugEtotalarray(indx(1))-debugEtotalarray(indx(1)-1)))
       elseif (indx(1) == ndensmat)then
          CFG%MaxLineSearchEnergyDiff=ABS(energy_array(indx(1))-energy_array(indx(1)-1))
          arh_Ediff_actual=ABS(debugEtotalarray(indx(1))-debugEtotalarray(indx(1)-1))               
       endif
       write(arh%lupri,'(a,ES20.11)') 'Old energy',arh%old_energy
       !if Old energy + change is not equal to Energy_array it is the first time the
       !linesearch is performed and it would be cunfusing to have the energy change printed
       do i=1,ndensmat
          write(arh%lupri,'(a,ES20.11)') 'Energy', energy_array(i)
       end do
       do i=1,ndensmat
          write(arh%lupri,*) 'Energy NF:', energy_array(i)
       end do
       write(arh%lupri,'(a,i4)') 'Minimum function value for ', indx(1)
       CFG%arh_linesE=energy_array(indx(1))
       write(arh%lupri,'(a,ES20.11)') 'Energy From Linesearch   ', CFG%arh_linesE 
       write(arh%lupri,'(a,ES20.11)') 'Maximum Energy difference', CFG%MaxLineSearchEnergyDiff
       ls%setting%scheme%LSDASCREEN = saveLSDASCREEN 
       SaveF0andD0 = saveSaveF0andD0 
       if (CFG%MaxLineSearchEnergyDiff > ABS(debugEtotalarray(indx(1))-CFG%arh_linesE)) then
          WRITE(arh%lupri,'(A)') 'default no error in linesearch'
          WRITE(arh%lupri,'(A,ES22.13)')'actual Etotal from Fock Matrix  ',debugEtotalarray(indx(1))
          WRITE(arh%lupri,'(A,ES22.13)')'linesearch Energy Estimate      ',CFG%arh_linesE
          WRITE(arh%lupri,'(A,ES22.13)')'Energy difference between estimate and actual Energy',&
               & ABS(debugEtotalarray(indx(1))-CFG%arh_linesE)
          WRITE(arh%lupri,'(A,ES22.13)')'Maximum difference in linesearch energies',&
               & CFG%MaxLineSearchEnergyDiff
          WRITE(arh%lupri,'(A,ES22.13)')'Maximum difference in actual energies',&
               & arh_Ediff_actual
!          IF(ABS(CFG%MaxLineSearchEnergyDiff) < 1.0E-10_realk)THEN
!             !turn off incremental build in linesearch, to improve accuracy
!             !              IF(config%opt%cfg_saveF0andD0)THEN
!             WRITE(arh%lupri,*)'Turn off incremental linesearch energy estimates to improve accuracy if cfg_saveF0andD0'
!             !              ENDIF
!          ENDIF
          IF(CFG%MaxLineSearchEnergyDiff < & 
               & 5.0E+2_realk*LS%SETTING%SCHEME%THRESHOLD*LS%SETTING%SCHEME%J_THR.OR.&
               & CFG%MaxLineSearchEnergyDiff < &
               & 1.5E+2_realk*ABS(debugEtotalarray(indx(1))-CFG%arh_linesE))THEN
             !improve accuracy on Fock matrix as the energy differences in the linesearch are small
             !compared to Fock matrix accuracy
             WRITE(arh%lupri,'(A)')'The SCF convergence have reached a point where we should'
             WRITE(arh%lupri,'(A)')'increase the accuracy of the Fock matrix build or increase the'
             WRITE(arh%lupri,'(A)')'accuracy of the linesearch energy estimates. We do both!'
          ENDIF
       endif
    enddo
    ls%SETTING%SCHEME%THRESHOLD = saveTHRESHOLD

    do i=1,ndensmat
       call mat_free(debugFockMatArray(i))
    enddo
    call mem_dealloc(debugFockMatArray)
    call mem_dealloc(EDFTY)
    call mem_dealloc(debugEtotalarray)
    write(arh%lupri,*)
    write(arh%lupri,'(a)') '***** END DEBUG LINESEARCH INFORMATION *****'
    write(arh%lupri,*)
  end subroutine debug_arh_LineSearch

  END MODULE ARHmodule
