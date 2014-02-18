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
   use davidson_solv_mod
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
    call mat_init(wrk,ndim,ndim)
    call get_OAO_gradient(decomp%FU, decomp%DU, wrk) !wrk = gradient
    call mat_scal(0.25E0_realk,wrk) !To match linear transformation, also divided by 4!


    arh%OAO_gradnrm = sqrt(mat_sqnorm2(wrk)) 
    arh%OAO_gradnrm_exist = .true.
    
    write (arh%lupri,*) 'OAO gradnorm', arh%OAO_gradnrm
    davidCFG%arh_gradnorm=arh%OAO_gradnrm

    call mat_init(wrk2,ndim,ndim)
    call mat_init(x,ndim,ndim)
    if (arh%cfg_nodamp) then
       call arh_PCG(arh, decomp, wrk, x, 0.0E0_realk, arh_antisymmetric,fifoqueue)
    else if (arh%cfg_arh_truncate .or. arh%cfg_arh_crop) then
       if (davidCFG%arh_davidson) then
            davidCFG%symm=arh_antisymmetric
            CALL LSTIMER('START ',t1,t2,decomp%lupri)
	    call arh_davidson_solver(davidCFG,arh,decomp,wrk,X,SCF_iteration,H1,wrk2,ls) 
            CALL LSTIMER('ARH DAVIDS. ',t1,t2,arh%lupri)
       else
            call arh_crop_solver(decomp, arh, debug, wrk, arh_antisymmetric, x, fifoqueue)
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
     call mat_free(x)
     call mat_free(wrk)

     !Stinne 24/1-07: Why tranform to AO basis? We already calculated the corresponding Fock matrix...
     !30/1-07: Because otherwise a wrong density is written in dens.restart!!
     call x_from_oao_basis(decomp,wrk2, Dnew) 

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
  real(realk)             :: t1, t2, gradnorm
  real(realk)    :: stepfactor
  type(matrix)   :: ax,G,scratchMat
  real(realk)    :: Eold,scaling,saveTHRESHOLD
  real(realk)    :: arh_Ediff_actual,Ecurrent
  integer        :: indx(1),i,j,counter,idebug
  logical        :: saveLSDASCREEN,saveSaveF0andD0
  logical        :: TightLineSearchEsitmatesThresholds
  type(matrix),pointer :: debugFockMatArray(:)
  real(realk),pointer  :: EDFTY(:),debugEtotalarray(:)
  real(realk)          :: xnorm,tmp,max_alpha
  real(realk)   :: coef(2)
  !For Energy Calculation
  type(matrix),pointer   :: densmat(:)
  integer :: ndensmat
  real(realk),pointer  :: alpha(:),energy_array(:)
  !For Storage
  integer :: nEnergies
  integer,pointer :: tmpindex(:)
  real(realk),pointer  :: fullalpha(:),fullenergy_array(:)
  real(realk) :: Optimal_alpha,dd2,dd3,EstVal
  logical :: CalcEnergyForFit,nondebug
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
        nondebug = .TRUE.
        debugLinesearch: IF(nondebug)THEN
           counter = counter+1
           !Scale gradient since davdison solver is made for (H-mu)X= -G
           CALL LSTIMER('START ',t1,t2,decomp%lupri)
           call davidson_solver(CFG,G,X)
           CALL LSTIMER('DAVID. SOLVER ',t1,t2,arh%lupri)
           IF(ABS(CFG%mu).GT.0.0E-14_realk)THEN
              ndensmat = 10
              call mem_alloc(alpha,ndensmat)
              call mem_alloc(densmat,ndensmat)
              call mem_alloc(energy_array,ndensmat)
              
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
                 call mem_dealloc(alpha)
                 call mem_dealloc(densmat)
                 call mem_dealloc(energy_array)              
                 cycle WHILELOOP
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
              write(arh%lupri,'(a,ES14.7)') 'Alpha Value ',alpha(indx(1))
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
              call mem_dealloc(alpha)
              call mem_dealloc(densmat)
              call mem_dealloc(energy_array)              
              EXIT WHILELOOP             
           ELSE
              !no linesearch
              arh%xnorm = sqrt(mat_sqnorm2(x))
              write(arh%lupri,'(a,ES14.7)') 'No linesearch due to no levelshift, step size= ', arh%xnorm
              !Now we construct the new density from X:
              CALL LSTIMER('START ',t1,t2,decomp%lupri)
              call oao_density_param(x,decomp%DU,wrk) !wrk = new D(X)
              call oao_purify(wrk,wrk2) !wrk2 = new purified density in oao basis
              CALL LSTIMER('NEW D ',t1,t2,arh%lupri)
              EXIT WHILELOOP
           ENDIF
        ELSE !debugLinesearch
           !this is the possible way new way with a several Fits
           counter = counter+1
           !Scale gradient since davdison solver is made for (H-mu)X= -G
           CALL LSTIMER('START ',t1,t2,decomp%lupri)
           call davidson_solver(CFG,G,X)
           CALL LSTIMER('DAVID. SOLVER ',t1,t2,arh%lupri)
           CALL LSTIMER('START ',t1,t2,decomp%lupri)

           IF(ABS(CFG%mu).GT.0.0E-14_realk)THEN
              nEnergies = 1
              call mem_alloc(fullalpha,nEnergies)
              call mem_alloc(fullenergy_array,nEnergies)
              nEnergies = 0
              ! step 1 we calculate Energy of D(X)
              ndensmat = 1
              call mem_alloc(alpha,ndensmat)
              call mem_alloc(densmat,ndensmat)
              do i=1,ndensmat
                 call mat_init(densmat(i),X%nrow,X%ncol)
              end do
              call mem_alloc(energy_array,ndensmat)
              alpha(1) = 1.0E0_realk
              call build_Dmat_from_alpha_X(alpha(1),X,Densmat(1),decomp)
              !set the linesearch_thresholds
              call linesearch_thresholds(CFG,ls,SCF_iteration,TightLineSearchEsitmatesThresholds)
              IF(CFG%arh_davidson_debug)THEN
                 call debug_arh_LineSearch(Densmat,H1,energy_array,ndensmat,&
                      & CFG%lupri,CFG%lupri,ls,CFG%LSmodthresh)
              ELSE
                 call di_SCF_EnergyCont(Densmat,H1,energy_array,ndensmat,&
                      & CFG%lupri,CFG%lupri,ls,CFG%LSmodthresh)
              ENDIF
              call add_energy_array_to_full(nEnergies,fullenergy_array,fullalpha,&
                   & ndensmat,energy_array,alpha)
              call mem_dealloc(alpha)
              call mem_dealloc(energy_array)
              do i=1,ndensmat
                 call mat_free(densmat(i))
              end do
              call mem_dealloc(densmat)
              do i=1,nEnergies
                 WRITE(ls%lupri,*)'energy_array(',i,')',fullenergy_array(i)
              enddo
              
              ! If found energy is larger than old energy, reduce stepsize and call solver
              if (fullenergy_array(1) > arh%old_energy) then
                 WRITE(ls%lupri,*)'found energy  ',fullenergy_array(1)
                 WRITE(ls%lupri,*)'arh%old_energy',arh%old_energy
                 WRITE(ls%lupri,*)'step rejection'
                 CFG%stepsize = (2.0_realk/3.0_realk)*CFG%stepsize
                 cycle WHILELOOP
              end if
              ! We have found a good direction now determine the step.
              
              ! we determine the maximum stepsize max_alpha
              if (CFG%arh_gradnorm < 1E-2_realk) stepfactor=0.25_realk
              if (CFG%arh_gradnorm .ge. 1E-2_realk) stepfactor=0.5_realk
              if (CFG%arh_gradnorm .ge. 1E-1_realk) stepfactor=0.75_realk
              if (CFG%arh_gradnorm .ge. 1E0_realk) stepfactor=1.0_realk
              ndensmat = 10
              scaling = 1.0_realk-stepfactor
              do i=1,ndensmat
                 scaling = scaling + stepfactor 
                 if (abs(scaling) <1E-5_realk) scaling=scaling+stepfactor 
                 !      alpha(i)=scaling
              end do
              !  max_alpha = 6*CFG%stepsize
              max_alpha = scaling
              write(ls%lupri,*)'stepfactor',stepfactor,'max_alpha',max_alpha,'CFG%stepsize',CFG%stepsize
              ! we make ndensmat matrices for a linesearch  
              ndensmat = 3 
              call mem_alloc(alpha,ndensmat)
              call mem_alloc(densmat,ndensmat)
              call mem_alloc(energy_array,ndensmat)
              do i = 1, ndensmat
                 alpha(i) = i*max_alpha/ndensmat
                 call mat_init(densmat(i),X%nrow,X%ncol)
                 call build_Dmat_from_alpha_X(alpha(i),X,Densmat(i),decomp)
              end do
              !set the linesearch_thresholds
              call linesearch_thresholds(CFG,ls,SCF_iteration,TightLineSearchEsitmatesThresholds)
              IF(CFG%arh_davidson_debug)THEN
                 call debug_arh_LineSearch(Densmat,H1,energy_array,ndensmat,&
                      & CFG%lupri,CFG%lupri,ls,CFG%LSmodthresh)
              ELSE
                 call di_SCF_EnergyCont(Densmat,H1,energy_array,ndensmat,&
                      & CFG%lupri,CFG%lupri,ls,CFG%LSmodthresh)
              ENDIF
              call add_energy_array_to_full(nEnergies,fullenergy_array,fullalpha,&
                   & ndensmat,energy_array,alpha)
              call mem_dealloc(alpha)
              call mem_dealloc(energy_array)
              do i=1,ndensmat
                 call mat_free(densmat(i))
              end do
              call mem_dealloc(densmat)
              
              ! we fit the energies to a quadratic fit 
              call fit_arh_LineSearch2(nEnergies,fullenergy_array,&
                   & CFG%lupri,fullalpha,Optimal_alpha,dd2,dd3,EstVal)
              
              WRITE(ls%lupri,*)'Optimal_alpha     :',Optimal_alpha
              WRITE(ls%lupri,*)'Standard Diviation:',dd2
              do i=1,nEnergies
                 WRITE(ls%lupri,*)'fit_energy_array(',i,')',fullenergy_array(i),fullalpha(i)
              enddo
              
              IF(dd2.GT.0.1E0_realk.AND.(Optimal_alpha.LT.fullalpha(1).OR.Optimal_alpha.GT.fullalpha(nEnergies)))THEN
                 Write(ls%lupri,*)'Warning: Clealy not quadratic surface standard diviation too large'
                 !we must chose a conservative alpha
                 !looks like a dobbel well 
                 indx = MINLOC(fullenergy_array)   
                 ndensmat = 0
                 do i=1,nEnergies
                    IF(fullalpha(i).LE.max_alpha+0.0E-5_realk)THEN
                       ndensmat = ndensmat + 1               
                    ENDIF
                 enddo
                 call mem_alloc(tmpindex,ndensmat)
                 call mem_alloc(alpha,ndensmat)
                 call mem_alloc(energy_array,ndensmat)
                 ndensmat = 0
                 do i=1,nEnergies
                    IF(fullalpha(i).LE.max_alpha+0.0E-5_realk)THEN
                       ndensmat = ndensmat + 1               
                       alpha(ndensmat) = fullalpha(i)
                       energy_array(ndensmat) = fullenergy_array(i)
                       tmpindex(ndensmat) = i
                    ENDIF
                 enddo
                 indx = MINLOC(energy_array)   
                 indx(1) = tmpindex(indx(1))
                 call mem_dealloc(tmpindex)
                 call mem_dealloc(alpha)
                 call mem_dealloc(energy_array)
                 Write(ls%lupri,*)'due to non quadratic surface we chose conservative step'
                 Write(ls%lupri,*)'indx  ',indx(1)
                 Write(ls%lupri,*)'Energy',fullenergy_array(indx(1))
                 Write(ls%lupri,*)'alpha ',fullalpha(indx(1))
              ELSE
                 !if the optimal_alpha is far from largest alpha the fit is probably bad
                 !so we add a few more densities in between
                 IF(Optimal_alpha .GT. fullalpha(nEnergies)+fullalpha(nEnergies)/2)THEN
                    !Calculate the energy of the Density from the fit
                    ndensmat = 2
                    call mem_alloc(alpha,ndensmat)
                    alpha(1) = fullalpha(nEnergies) + (Optimal_alpha - fullalpha(nEnergies))/2
                    alpha(2) = Optimal_alpha
                    call mem_alloc(densmat,ndensmat)
                    call mem_alloc(energy_array,ndensmat)
                    do i = 1, ndensmat
                       call mat_init(densmat(i),X%nrow,X%ncol)
                       call build_Dmat_from_alpha_X(alpha(i),X,Densmat(i),decomp)
                    end do
                    !set the linesearch_thresholds
                    call linesearch_thresholds(CFG,ls,SCF_iteration,TightLineSearchEsitmatesThresholds)
                    IF(CFG%arh_davidson_debug)THEN
                       call debug_arh_LineSearch(Densmat,H1,energy_array,ndensmat,&
                            & CFG%lupri,CFG%lupri,ls,CFG%LSmodthresh)
                    ELSE
                       call di_SCF_EnergyCont(Densmat,H1,energy_array,ndensmat,&
                            & CFG%lupri,CFG%lupri,ls,CFG%LSmodthresh)
                    ENDIF
                    call add_energy_array_to_full(nEnergies,fullenergy_array,fullalpha,&
                         & ndensmat,energy_array,alpha)
                    call mem_dealloc(alpha)
                    call mem_dealloc(energy_array)
                    do i = 1, ndensmat
                       call mat_free(densmat(i))
                    enddo
                    call mem_dealloc(densmat)
                    
                    ! we fit the energies to a quadratic fit 
                    call fit_arh_LineSearch2(nEnergies,fullenergy_array,&
                         & CFG%lupri,fullalpha,Optimal_alpha,dd2,dd3,EstVal)
                    
                    WRITE(ls%lupri,*)'Optimal_alpha     :',Optimal_alpha
                    WRITE(ls%lupri,*)'Standard Diviation:',dd2
                    do i=1,nEnergies
                       WRITE(ls%lupri,*)'fit_energy_array(',i,')',fullenergy_array(i),fullalpha(i)
                    enddo
                 ENDIF
                 
                 IF(dd2.GT.0.1E0_realk.AND.(Optimal_alpha.LT.fullalpha(1).OR.Optimal_alpha.GT.fullalpha(nEnergies)))THEN
                    Write(ls%lupri,*)'Warning: Clealy not quadratic surface standard diviation too large'
                    !we must chose a conservative alpha
                    !looks like a dobbel well 
                    indx = MINLOC(fullenergy_array)   
                    ndensmat = 0
                    do i=1,nEnergies
                       IF(fullalpha(i).LE.max_alpha+0.0E-5_realk)THEN
                          ndensmat = ndensmat + 1               
                       ENDIF
                    enddo
                    call mem_alloc(tmpindex,ndensmat)
                    call mem_alloc(alpha,ndensmat)
                    call mem_alloc(energy_array,ndensmat)
                    ndensmat = 0
                    do i=1,nEnergies
                       IF(fullalpha(i).LE.max_alpha+0.0E-5_realk)THEN
                          ndensmat = ndensmat + 1               
                          alpha(ndensmat) = fullalpha(i)
                          energy_array(ndensmat) = fullenergy_array(i)
                          tmpindex(ndensmat) = i
                       ENDIF
                    enddo
                    indx = MINLOC(energy_array)   
                    indx(1) = tmpindex(indx(1))
                    call mem_dealloc(tmpindex)
                    call mem_dealloc(alpha)
                    call mem_dealloc(energy_array)
                    Write(ls%lupri,*)'due to non quadratic surface we chose conservative step'
                    Write(ls%lupri,*)'indx  ',indx(1)
                    Write(ls%lupri,*)'Energy',fullenergy_array(indx(1))
                    Write(ls%lupri,*)'alpha ',fullalpha(indx(1))
                 ELSE
                    CalcEnergyForFit = CFG%arh_davidson_debug
                    IF(CalcEnergyForFit)THEN 
                       !Calculate the energy of the Density from the fit
                       ndensmat = 1
                       call mem_alloc(alpha,ndensmat)
                       alpha(1) = Optimal_alpha
                       call mem_alloc(densmat,ndensmat)
                       call mem_alloc(energy_array,ndensmat)
                       call mat_init(densmat(1),X%nrow,X%ncol)
                       call build_Dmat_from_alpha_X(Optimal_alpha,X,Densmat(1),decomp)
                       !set the linesearch_thresholds
                       call linesearch_thresholds(CFG,ls,SCF_iteration,TightLineSearchEsitmatesThresholds)
                       IF(CFG%arh_davidson_debug)THEN
                          call debug_arh_LineSearch(Densmat,H1,energy_array,ndensmat,&
                               & CFG%lupri,CFG%lupri,ls,CFG%LSmodthresh)
                       ELSE
                          call di_SCF_EnergyCont(Densmat,H1,energy_array,ndensmat,&
                               & CFG%lupri,CFG%lupri,ls,CFG%LSmodthresh)
                       ENDIF
                       WRITE(ls%lupri,*)'the energy from fit',energy_array(1)
                       call add_energy_array_to_full(nEnergies,fullenergy_array,fullalpha,&
                            & ndensmat,energy_array,alpha)
                       call mem_dealloc(alpha)
                       call mem_dealloc(energy_array)
                       call mat_free(densmat(1))
                       call mem_dealloc(densmat)
                       
                       ! we fit the energies to a quadratic fit 
                       call fit_arh_LineSearch2(nEnergies,fullenergy_array,&
                            & CFG%lupri,fullalpha,Optimal_alpha,dd2,dd3,EstVal)
                       
                       WRITE(ls%lupri,*)'Optimal_alpha     :',Optimal_alpha
                       WRITE(ls%lupri,*)'Standard Diviation:',dd2
                       do i=1,nEnergies
                          WRITE(ls%lupri,*)'fit_energy_array(',i,')',fullenergy_array(i),fullalpha(i)
                       enddo
                       ndensmat = 1
                       call mem_alloc(alpha,ndensmat)
                       alpha(1) = Optimal_alpha
                       call mem_alloc(densmat,ndensmat)
                       call mem_alloc(energy_array,ndensmat)
                       call mat_init(densmat(1),X%nrow,X%ncol)
                       call build_Dmat_from_alpha_X(Optimal_alpha,X,Densmat(1),decomp)
                       !set the linesearch_thresholds
                       energy_array(1) = -EstVal
                       call add_energy_array_to_full(nEnergies,fullenergy_array,fullalpha,&
                            & ndensmat,energy_array,alpha)
                       call mem_dealloc(alpha)
                       call mem_dealloc(energy_array)
                       call mat_free(densmat(1))
                       call mem_dealloc(densmat)
                    ELSE
                       ndensmat = 1
                       call mem_alloc(alpha,ndensmat)
                       alpha(1) = Optimal_alpha
                       call mem_alloc(densmat,ndensmat)
                       call mem_alloc(energy_array,ndensmat)
                       call mat_init(densmat(1),X%nrow,X%ncol)
                       call build_Dmat_from_alpha_X(Optimal_alpha,X,Densmat(1),decomp)
                       !set the linesearch_thresholds
                       energy_array(1) = -EstVal
                       call add_energy_array_to_full(nEnergies,fullenergy_array,fullalpha,&
                            & ndensmat,energy_array,alpha)
                       call mem_dealloc(alpha)
                       call mem_dealloc(energy_array)
                       call mat_free(densmat(1))
                       call mem_dealloc(densmat)
                    ENDIF

                    !we choose the Density with the lowest energy
                    indx = MINLOC(fullenergy_array)   
                 ENDIF
              ENDIF
              CFG%MaxLineSearchEnergyDiff = 10000.0E0_realk


              write(arh%lupri,*)
              write(arh%lupri,'(a)') '***** LINESEARCH INFORMATION *****'
              write(arh%lupri,*)
              if (CFG%arh_davidson_debug) write(arh%lupri,'(a,ES20.11)') 'Old energy',arh%old_energy
              xnorm=dsqrt(mat_dotproduct(x,x))
              if (CFG%arh_davidson_debug) then
                 do i=1,nEnergies
                    write(arh%lupri,'(a,ES20.11,a,f10.4)') 'Energy', fullenergy_array(i), &
                         &' Tot. stepsize:', fullalpha(i)*xnorm 
                 end do
              endif
              call mat_scal(fullalpha(indx(1)),X)
              arh%xnorm = sqrt(mat_sqnorm2(x))
              write(arh%lupri,'(a,ES14.7)') 'Alpha Value ',fullalpha(indx(1))
              write(arh%lupri,'(a,i4,a,ES14.7)') 'ARHLS: Minimum function value for ', indx(1),&
                   &' using step size ', arh%xnorm
              CFG%arh_linesE=fullenergy_array(indx(1))
              if (CFG%arh_davidson_debug) write(arh%lupri,'(a,ES20.11)') &
                   &'ARHLS: Energy From Linesearch   ',CFG%arh_linesE 
              write(arh%lupri,*)
              write(arh%lupri,'(a)') '***** END  LINESEARCH INFORMATION *****'
              write(arh%lupri,*)
              call mat_init(scratchMat,X%nrow,X%ncol)
              call oao_density_param(X,decomp%DU,scratchMat) !D = new D(X)
              call oao_purify(scratchMat,wrk2) !wrk2 = new purified density in oao basis
              call mat_free(scratchMat)
              CALL LSTIMER('NEW D ',t1,t2,arh%lupri)
              EXIT WHILELOOP
           ELSE
              !no linesearch
              arh%xnorm = sqrt(mat_sqnorm2(x))
              write(arh%lupri,'(a,ES14.7)') 'No linesearch due to no levelshift, step size= ', arh%xnorm
              !Now we construct the new density from X:
              CALL LSTIMER('START ',t1,t2,decomp%lupri)
              call oao_density_param(x,decomp%DU,wrk) !wrk = new D(X)
              call oao_purify(wrk,wrk2) !wrk2 = new purified density in oao basis
              CALL LSTIMER('NEW D ',t1,t2,arh%lupri)
              EXIT WHILELOOP
           endif
        endif debugLinesearch
     else
        CALL LSTIMER('START ',t1,t2,decomp%lupri)
        call davidson_solver(CFG,G,X)
        CALL LSTIMER('DAVID. SOLVER ',t1,t2,arh%lupri)
        arh%xnorm=dsqrt(mat_sqnorm2(X))
        !Now we construct the new density from X:
        CALL LSTIMER('START ',t1,t2,decomp%lupri)
        call oao_density_param(x,decomp%DU,wrk) !wrk = new D(X)
        call oao_purify(wrk,wrk2) !wrk2 = new purified density in oao basis
        CALL LSTIMER('NEW D ',t1,t2,arh%lupri)
        EXIT WHILELOOP
     end if Linesearch
  end do WHILELOOP

  !print info 
  call print_info_davidson(CFG,arh%xnorm,SCF_iteration)
  call mat_free(G)
  call mat_free(ax)

end subroutine arh_davidson_solver

subroutine add_energy_array_to_full(nEnergies,fullenergy_array,fullalpha,&
     & ndensmat,energy_array,alpha)
  implicit none
  integer :: nEnergies,ndensmat
  real(realk),pointer  :: alpha(:),energy_array(:)
  !For Storage
  real(realk),pointer  :: fullalpha(:),fullenergy_array(:)
  integer :: n1,i
  real(realk),pointer  :: tmp(:)
  n1 = nEnergies + ndensmat
  call mem_alloc(tmp,n1)
  do i = 1,nEnergies
     tmp(i) = fullalpha(i)
  enddo
  do i = 1,nDensmat
     tmp(nEnergies+i) = alpha(i)
  enddo
  call mem_dealloc(fullalpha)
  call mem_alloc(fullalpha,n1)
  do i = 1,n1
     fullalpha(i) = tmp(i)
  enddo
  do i = 1,nEnergies
     tmp(i) = fullenergy_array(i)
  enddo
  do i = 1,nDensmat
     tmp(nEnergies+i) = energy_array(i)
  enddo
  call mem_dealloc(fullenergy_array)
  call mem_alloc(fullenergy_array,n1)
  do i = 1,n1
     fullenergy_array(i) = tmp(i)
  enddo
  nEnergies=n1
end subroutine add_energy_array_to_full

subroutine build_Dmat_from_alpha_X(alpha,X,D,decomp)
  implicit none
  real(realk),intent(in) :: alpha
  type(matrix),intent(in) :: X
  type(matrix),intent(inout) :: D
  !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
  type(decompItem),intent(inout) :: decomp
  !
  type(matrix) :: wrk
  call mat_init(wrk,X%nrow,X%ncol)
  call mat_copy(alpha,X,wrk)
  call oao_density_param(wrk,decomp%DU,D) !D = new D(X)
  call oao_purify(D,wrk) !wrk2 = new purified density in oao basis
  call x_from_oao_basis(decomp,wrk,D) 
  call mat_free(wrk)
end subroutine build_Dmat_from_alpha_X

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
!      call mat_init(xF,rowdim,coldim)
!      call mat_init(sigmaF,rowdim,coldim)
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
   call mat_init(scrmat,rowdim,coldim)
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
   call mat_init(res,rowdim,coldim)
   call mat_add(-1E0_realk,Grad,-1E0_realk, scrmat, res) !res = 1st residual
   call mat_daxpy(mu,x,res)

   !Preconditioning:
   call mat_init(resP,rowdim,coldim)
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

   if (arh%cfg_arh_newdamp) then
      call mat_init(xF,rowdim,coldim)
      call mat_init(sigmaF,rowdim,coldim)
   endif
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
   call mat_free(resP)
   call mat_free(res)
   if (arh%cfg_arh_newdamp) then
      call mat_free(xF)
      call mat_free(sigmaF)
   endif

   arh%current_mu = mu

   !if (arh%debug_diag_redspace .and. .not. arh%cfg_arh_crop) then
   !   write (arh%lupri,'("Final lowest redspace eigenvalue: ",F12.6, "      ")') debug%final_redspace_eival
   !endif

   CALL LSTIMER('CROP solver',t1,t2,arh%lupri)

   if (arh%cfg_arh_crop) call mat_scal(-1.0E0_realk,x) !HACK!!!!

   call arh_get_TR_denom(arh, Grad, x, scrmat, decomp%cfg_unres, mu)
   call mat_free(scrmat)

   if (arh%cfg_arh_truncate) then
      call modfifo_free(vectorsubspace)
   else
      CALL LSCLOSE(lusigma,'DELETE')
      CALL LSCLOSE(lub,'DELETE')
   endif

   INQUIRE(file='xsave',EXIST=fileexists,OPENED=fileopened)
   if (fileexists.AND.fileopened) call LSCLOSE(xsave_lu,'DELETE')

   call mem_dealloc(arh%Ared)
   call mem_dealloc(arh%Gred)
   call mem_dealloc(arh%Sred)
   call mem_dealloc(arh%CROPmat)
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
  


  write(davidCFG%lupri,'(i5,a,f7.3,a,f7.3,a,f6.2,a,i3,a)') it+1,' step norm =',step,' trust radius = '&
  &,davidCFG%stepsize,' mu = ',davidCFG%mu, ' micro it. = ',davidCFG%it, '  %#%' 

 

  end subroutine print_info_davidson
  
  subroutine debug_arh_LineSearch(densmat,H1,energy_array,ndensmat,lupri,luerr,ls,LSmodthresh)
    implicit none
    type(lsitem) ::ls
    !> One electron Hamiltonian
    type(matrix),intent(in) :: H1
    !> Current SCF iteration
    integer,intent(in)      :: ndensmat,lupri,luerr
    type(matrix),intent(inout)   :: densmat(ndensmat)
    real(realk),intent(inout)    :: energy_array(ndensmat),LSmodthresh
    !
    logical :: DaJengine,DaCoulomb

    write(lupri,*)
    write(lupri,'(a)') '***** DEBUG LINESEARCH INFORMATION *****'
    write(lupri,*)
    DaJengine = ls%setting%scheme%LSDaJengine
    DaCoulomb = ls%setting%scheme%LSDaCoulomb
    ls%setting%scheme%LSDaJengine = .TRUE.
    ls%setting%scheme%LSDaCoulomb = .FALSE.
    write(lupri,*)'DaJengine',ls%setting%scheme%LSDaJengine
    write(lupri,*)'DaCoulomb',ls%setting%scheme%LSDaCoulomb
    write(lupri,*)'LSDASCREEN',ls%setting%scheme%LSDASCREEN
    call di_SCF_EnergyCont(densmat,H1,energy_array,ndensmat,&
         & lupri,luerr,ls,LSmodthresh)

    ls%setting%scheme%LSDaJengine = .FALSE.
    ls%setting%scheme%LSDaCoulomb = .TRUE.
    write(lupri,*)'DaJengine',ls%setting%scheme%LSDaJengine
    write(lupri,*)'DaCoulomb',ls%setting%scheme%LSDaCoulomb
    write(lupri,*)'LSDASCREEN',ls%setting%scheme%LSDASCREEN
    call di_SCF_EnergyCont(densmat,H1,energy_array,ndensmat,&
         & lupri,luerr,ls,LSmodthresh)

    ls%setting%scheme%LSDaJengine = DaJengine
    ls%setting%scheme%LSDaCoulomb = DaCoulomb

    write(lupri,*)
    write(lupri,'(a)') '***** END DEBUG LINESEARCH INFORMATION *****'
    write(lupri,*)
  end subroutine debug_arh_LineSearch

  subroutine fit_arh_LineSearch2(ndensmat,energy_array,lupri,alpha,optimal_alpha,dd2,dd3,EstVal)
    implicit none
    !> Current SCF iteration
    integer,intent(in)           :: ndensmat,lupri
    real(realk),intent(inout)    :: energy_array(ndensmat),alpha(ndensmat)
    real(realk),intent(inout)    :: optimal_alpha,dd2,dd3,EstVal
    !
    integer,parameter :: MYSIZE=25
    integer :: morder,lorder,i
    real(realk) :: xarray(MYSIZE),yarray(MYSIZE),carray(0:MYSIZE),e1,ddx
    ! m  : the order of the fit (2)
    ! e1 : the error reduction factor (0)
    ! n  : number of datapoints
    ! l  : the order found (I think)
    ! x  : X coord vector 
    ! y  : Y coord vector 
    ! c  : coeff vector 
    !dd  : standard diviation
    morder = 2 
    e1 = 0
    do i = 1,ndensmat
       xarray(i) = alpha(i)
       yarray(i) = -energy_array(i)
    enddo
    call LS_POLY(morder,e1,ndensmat,lorder,xarray,yarray,carray,dd2)
    write(lupri,*)'Coefficients are: '
    print *,' '
    do i=0, lorder
       write(lupri,*)  i, carray(i)
    end do
!    coef(1) = carray(1); coef(2) = carray(2)
    write(lupri,*)'FIT: standard diviation: ', dd2
    optimal_alpha = -carray(1)/(2*carray(2))    
    write(lupri,*)'FIT: optimal_alpha = -b/(2a) = ',optimal_alpha    
    Estval = carray(2)*optimal_alpha**2+carray(1)*optimal_alpha+carray(0)
    write(lupri,*)'FIT: Estimated value: ',-Estval

    morder = 3
    e1 = 0
    do i = 1,ndensmat
       xarray(i) = alpha(i)
       yarray(i) = -energy_array(i)
    enddo
    call LS_POLY(morder,e1,ndensmat,lorder,xarray,yarray,carray,dd3)
    write(lupri,*)'FIT Coefficients are: '
    print *,' '
    do i=0, lorder
       write(lupri,*)  i, carray(i)
    end do
    write(lupri,*)'FIT: standard diviation3: ', dd3

    IF(ndensmat.EQ.20)THEN      
       morder = 4
       e1 = 0
       do i = 1,ndensmat
          xarray(i) = alpha(i)
          yarray(i) = -energy_array(i)
       enddo
       call LS_POLY(morder,e1,ndensmat,lorder,xarray,yarray,carray,ddx)
       write(lupri,*)'FIT Coefficients are: '
       print *,' '
       do i=0, lorder
          write(lupri,*)  i, carray(i)
       end do       
       write(lupri,*)'FIT: standard diviation: ', ddx

       morder = 5
       e1 = 0
       do i = 1,ndensmat
          xarray(i) = alpha(i)
          yarray(i) = -energy_array(i)
       enddo
       call LS_POLY(morder,e1,ndensmat,lorder,xarray,yarray,carray,ddx)
       write(lupri,*)'FIT Coefficients are: '
       print *,' '
       do i=0, lorder
          write(lupri,*)  i, carray(i)
       end do       
       write(lupri,*)'FIT: standard diviation: ', ddx

    ENDIF

  end subroutine fit_arh_LineSearch2

  !*****************************************************************
  !*         LEAST SQUARES POLYNOMIAL FITTING PROCEDURE            *
  !* ------------------------------------------------------------- *
  !* This program least squares fits a polynomial to input data.   *
  !* forsythe orthogonal polynomials are used in the fitting.      *
  !* The number of data points is n.                               *
  !* The data is input to the subroutine in x[i], y[i] pairs.      *
  !* The coefficients are returned in c[i],                        *
  !* the smoothed data is returned in v[i],                        *
  !* the order of the fit is specified by m.                       *
  !* The standard deviation of the fit is returned in d.           *
  !* There are two options available by use of the parameter e:    *
  !*  1. if e = 0, the fit is to order m,                          *
  !*  2. if e > 0, the order of fit increases towards m, but will  *
  !*     stop if the relative standard deviation does not decrease *
  !*     by more than e between successive fits.                   *
  !* The order of the fit then obtained is l.                      *
  !*****************************************************************
  ! m  : the order of the fit (2)
  ! e1 : the error reduction factor (0)
  ! n  : number of datapoints
  ! l  : the order found (I think)
  ! x  : X coord vector 
  ! y  : Y coord vector 
  ! c  : coeff vector 
  !dd  : standard diviation

  Subroutine LS_POLY(m,e1,n,l,x,y,c,dd)
    implicit none
    integer,parameter :: MYSIZE=25
    !Labels: 10,15,20,30,50
    real*8 x(MYSIZE),y(MYSIZE),v(MYSIZE),a(MYSIZE),b(MYSIZE)
    real*8 c(0:MYSIZE),d(MYSIZE),c2(MYSIZE),e(MYSIZE),f(MYSIZE)
    integer i,l,l2,m,n,n1
    real*8 a1,a2,b1,b2,c1,dd,d1,e1,f1,f2,v1,v2,w,vv
    n1 = m + 1; l=0
    v1 = 1.d7
    ! Initialize the arrays
    do i=1, n1
       a(i) = 0.d0; b(i) = 0.d0; f(i) = 0.d0
    end do
    do i=1, n
       v(i) = 0.d0; d(i) = 0.d0
    end do
    d1 = dsqrt(dfloat(n)); w = d1;
    do i=1, n
       e(i) = 1.d0 / w
    end do
    f1 = d1; a1 = 0.d0
    do i=1, n
       a1 = a1 + x(i) * e(i) * e(i)
    end do
    c1 = 0.d0
    do i=1, n
       c1 = c1 + y(i) * e(i)
    end do
    b(1) = 1.d0 / f1; f(1) = b(1) * c1
    do i=1, n
       v(i) = v(i) + e(i) * c1
    end do
    m = 1
    ! Save latest results
10  do i=1, l
       c2(i) = c(i)
    end do
    l2 = l; v2 = v1; f2 = f1; a2 = a1; f1 = 0.d0
    do i=1, n
       b1 = e(i)
       e(i) = (x(i) - a2) * e(i) - f2 * d(i)
       d(i) = b1
       f1 = f1 + e(i) * e(i)
    end do
    f1 = dsqrt(f1)
    do i=1, n
       e(i) = e(i) / f1
    end do
    a1 = 0.d0
    do i=1, n  
       a1 = a1 + x(i) * e(i) * e(i)
    end do
    c1 = 0.d0
    do i=1, n  
       c1 = c1 + e(i) * y(i)
    end do
    m = m + 1; i = 0
15  l = m - i; b2 = b(l); d1 = 0.d0
    if (l > 1)  d1 = b(l - 1)
    d1 = d1 - a2 * b(l) - f2 * a(l)
    b(l) = d1 / f1; a(l) = b2; i = i + 1
    if (i.ne.m) goto 15
    do i=1, n
       v(i) = v(i) + e(i) * c1 
    end do
    do i=1, n1
       f(i) = f(i) + b(i) * c1
       c(i) = f(i)
    end do
    vv = 0.d0
    do i=1, n
       vv = vv + (v(i) - y(i)) * (v(i) - y(i))
    end do
    !Note the division is by the number of degrees of freedom
    vv = dsqrt(vv / dfloat(n - l - 1)); l = m
    if (e1.eq.0.d0) goto 20
    !Test for minimal improvement
    if (dabs(v1 - vv) / vv < e1) goto 50
    !if error is larger, quit
    if (e1 * vv > e1 * v1) goto 50
    v1 = vv
20  if (m.eq.n1) goto 30
    goto 10
    !Shift the c[i] down, so c(0) is the constant term
30  do i=1, l  
       c(i - 1) = c(i)
    end do
    c(l) = 0.d0
    ! l is the order of the polynomial fitted
    l = l - 1; dd = vv
    return
    ! Aborted sequence, recover last values
50  l = l2; vv = v2
    do i=1, l  
       c(i) = c2(i)
    end do
    goto 30
  end Subroutine LS_POLY

  END MODULE ARHmodule
