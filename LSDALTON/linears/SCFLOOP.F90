!> @file 
!> Contains the scfloop routine and 
!> scfloop module
!> \author Stinne Hoest
!> \date 2008
module scfloop_module
  use precision
  use configurationType
  use TYPEDEFTYPE!, only: lsitem
  use Matrix_Module
  use Matrix_Operations
  use Matrix_Operations_aux!, only: mat_report_sparsity
  use Matrix_util!, only: get_AO_gradient
  use lstiming!, only: lstimer
  use ks_settings!, only: ks_init_linesearch_fock, ks_free_linesearch_fock, SaveF0andD0
  use IntegralInterfaceMOD!, only: II_setIncremental
  use queue_module!, only: modfifo
  use queue_ops!, only: modfifo_init, fifoqueue_on_disk, modfifo_free, fifoqueue_from_disk
  use direct_dens_util!, only: debugitem, dditem, dd_init, DD_homolumo_and_heseigen, dd_shutdown
  use direct_dens_util_unres!,only: dd_homolumo_and_heseigen_unres
  use av_utilities!, only: util_HistoryStore, queue_on_disk, queue_init, queue_free, queue_from_disk
  use diagonalization!, only: dopt_get_density_shutdown, dopt_get_density_init
  USE scf_stats!, only: scf_stats_init, scf_stats_update, &
!       & scf_stats_end_print,scf_stats_shutdown
  use files!, only: lsopen,lsclose
  use Fock_evaluator!, only: fck_unscale_virt
  uSE density_optimization!, only: dopt_get_density
  use decompMod!, only: get_oao_transformed_matrices
  use density_subspace_min!, only: density_subspace_minimization
  use trustradius_mod!, only: update_trustradius
  use arh_debugging!, only: debug_diag_full_hessian
  private 
  public :: scfloop, scf_afterplay
CONTAINS
!> \brief Contains the main SCF loop for LSDALTON
!> \author T. Kjaergaard, S. Reine
!> \date 2008-10-26
!> \param H1 One-electron Hamiltonian
!> \param F Fock/Kohn-Sham matrix
!> \param D Density matrix
!> \param S Overlap matrix
!> \param ls Contains information read from LSDALTON.INP 
!> \param config Contains all info about configuration/settings for SCF calculation
SUBROUTINE scfloop(H1,F,D,S,E,ls,config)
   implicit none   
   type(matrix),intent(in)         :: H1
   type(matrix),intent(inout)      :: D(1), F(1), S
   type(configItem),intent(inout)  :: config
   type(lsitem) :: ls
   real(realk)                     :: E(1)
!
   integer                         :: nbast
   TYPE(util_HistoryStore)         :: queue
   TYPE(modFIFO),target             :: fifoqueue
   TYPE(Matrix)                    :: grad,tempm1,tempm2,tempm3,tempm4,tempm5 
   real(realk)                     :: gradnrm, hessian_eigenval, ehomo, elumo
   integer                         :: iteration, matmul1, matmul2, matmultot
   integer                         :: number_atoms, queuesize,nnz,denslun,ndmat
   integer                         :: istart,ndmat2
   REAL                            :: tstart, tend, t0, MFLOPS, norm
   logical                         :: energy_converged
   logical     :: file_exists 
   real(realk) :: TSTR, TEN, TIMSTR, TIMEND, t1, t2,E2(2),gradnrm2
   type(matrix) :: D2(2), F2(2),grad2,tmpgrad
   type(matrix) :: x, Dchol !For debug
   type(matrix) :: unitmat
   type(matrix), pointer :: Dpointer, Fpointer
   LOGICAL :: dalink, incremental,onmaster,cs00,NotLastSCFLevel,gradalloc,ForcePrint
   real(realk) :: acceptratio, limitratio
   real(realk) :: h
   type(matrix) :: Dtest
   !
   ForcePrint = .TRUE.
   ndmat = 1
   OnMaster=.true.
   NotLastSCFLevel = config%integral%LOW_ACCURACY_START
   gradalloc = .FALSE.
   CALL LSTIMER('START',TSTR,TEN,config%LUPRI)
   !INITIALISE SCF CYCLE THINGS 
   nbast = H1%nrow
   acceptratio = config%av%cfg_settings(config%av%CFG_set_type)%min_density_overlap
   limitratio = config%av%cfg_settings(config%av%CFG_set_type)%max_dorth_ratio
   call DOPT_get_density_init(nbast,acceptratio,limitratio,config%diag)
   !call direct_dens_init
   if (config%opt%cfg_density_method == config%opt%cfg_f2d_arh) then
      queuesize = config%av%cfg_settings(config%av%cfg_set_type)%max_history_size
      !FIXME: this queue cannot be put on disk, get_from_modFIFO_disk won't work! 
      call modfifo_init(fifoqueue,queuesize,nbast,nbast,.false.)
   else
      call queue_init(config%av,queue)
   endif
   !if (DEBUG_DCHANGE) call debug_dchange_init
   
   call scf_stats_init(config%opt)
   
   call mat_report_sparsity(S,'S    ',nnz,config%lupri)
   call mat_report_sparsity(D(1),'AO D ',nnz,config%lupri)
   !    call mat_report_sparsity(F,'AO F ',nnz,config%lupri) not build yet
   write(config%lupri,*)' Relative convergence threshold for solver:', config%solver%cfg_micro_thresh
   WRITE(config%lupri,*)' SCF Convergence criteria for gradient norm:',config%opt%set_convergence_threshold
   
   CALL LSTIMER('INIT SCF',TSTR,TEN,config%LUPRI)
   
   ! attach pointers for RedSpaceItem to make use of davidson solver possible
   if (.not. config%solver%set_arhterms) write(config%lupri,*) ' ARH PART OF LINEAR TRANSFORMATIONS IS TURNED OFF ' 
   nullify(config%davidSCF%fifoqueue)
   config%davidSCF%arh => config%solver
   config%davidSCF%decomp => config%decomp
   config%davidSCF%fifoqueue => fifoqueue
   config%davidSCF%lupri = config%lupri
   config%davidSCF%stepsize=config%davidSCF%max_stepsize
   config%davidSCF%arh_linesearch=config%davidSCF%arh_inp_linesearch
   config%davidSCF%arh_extravec=config%davidSCF%arh_inp_extravec
   
   IF(config%opt%cfg_saveF0andD0)THEN
      call ks_init_linesearch_fock(nbast)
   ENDIF
   config%solver%OAO_gradnrm_exist=.false. ! reset
   gradnrm = 10000.0E0_realk
      
   istart = 1
   if(config%opt%cfg_start_guess=='TRILEVEL') then
      IF(ls%optlevel.EQ.3.AND.config%opt%add_atoms_start)THEN
         CALL LSTIMER('START',TIMSTR,TIMEND,config%LUPRI)
         dalink = ls%setting%scheme%DALINK
         ls%setting%scheme%DALINK = .false.
         cs00 = ls%setting%scheme%DFT%CS00
         ls%setting%scheme%DFT%CS00 = .FALSE.
         call mat_init(D2(1),D(1)%nrow,D(1)%ncol)
         call mat_init(D2(2),D(1)%nrow,D(1)%ncol)
         
         ndmat2 = 2
         call atoms_start(config,D2,H1,S,ls,ndmat2)
         call mat_assign(D2(2),D(1))

         call mat_init(F2(1),F(1)%nrow,F(1)%ncol)
         call mat_init(F2(2),F(1)%nrow,F(1)%ncol)
         iteration = 1
         incremental = .FALSE. !config%opt%cfg_incremental .AND. iteration.NE. 1
         call II_setIncremental(ls%setting%scheme,incremental)
         CALL get_fock(config, fifoqueue, queue, iteration, D2, H1, F2, ndmat2,E2,ls)
         IF(config%opt%cfg_oao_gradnrm)THEN
            call get_oao_transformed_matrices(config%decomp,F2(1),D2(1))
            call mat_init(grad,nbast,nbast)
            gradalloc = .TRUE.
            call get_OAO_gradient(config%decomp%FU, config%decomp%DU, grad) !wrk = gradient
            call mat_scal(0.25E0_realk,grad) !To match linear transformation, also divided by 4!
            gradnrm = sqrt(mat_sqnorm2(grad))

            call get_oao_transformed_matrices(config%decomp,F2(2),D2(2))
            call mat_init(grad2,nbast,nbast)
            call get_OAO_gradient(config%decomp%FU, config%decomp%DU, grad2) !wrk = gradient
            call mat_scal(0.25E0_realk,grad2) !To match linear transformation, also divided by 4!
            gradnrm2 = sqrt(mat_sqnorm2(grad2))


            iteration = 0
            call scf_stats_update(iteration,gradnrm,E2(1),config%opt)
            call scf_stats_debug_mem(config%lupri,iteration)
            iteration = 1
            call scf_stats_update(iteration,gradnrm2,E2(2),config%opt)
            call scf_stats_debug_mem(config%lupri,iteration)
            
            !grad is now AO grad
            CALL get_AO_gradient(F2(1), D2(1), S, grad)
            CALL get_AO_gradient(F2(2), D2(2), S, grad2)

         ELSE
            gradalloc = .TRUE.
            call mat_init(grad,nbast,nbast)
            CALL get_AO_gradient(F2(1), D2(1), S, grad)
            gradnrm = sqrt(mat_sqnorm2(grad))

            call mat_init(grad2,nbast,nbast)
            CALL get_AO_gradient(F2(2), D2(2), S, grad2)
            gradnrm2 = sqrt(mat_sqnorm2(grad2))

            iteration = 0
            call scf_stats_update(iteration,gradnrm,E2(1),config%opt)
            call scf_stats_debug_mem(config%lupri,iteration)
            iteration = 1
            call scf_stats_update(iteration,gradnrm2,E2(2),config%opt)
            call scf_stats_debug_mem(config%lupri,iteration)
         ENDIF
         iteration = 0
         config%solver%step_accepted = .TRUE.
         CALL Density_subspace_minimization(config, fifoqueue, queue, E2(1),S,H1, grad, F2(1), D2(1),iteration)
         call mat_free(grad)
         iteration = 1
         config%solver%step_accepted = .TRUE.
         CALL Density_subspace_minimization(config, fifoqueue, queue, E2(2),S,H1, grad2, F2(2), D2(2),iteration)
         call mat_free(grad2)
         gradalloc = .FALSE.
         call mat_assign(D(1),D2(2))
         call mat_assign(F(1),F2(2))
         gradnrm = gradnrm2
         E(1) = E2(2)
         istart = 2

         call mat_free(D2(1))
         call mat_free(D2(2))
         call mat_free(F2(1))
         call mat_free(F2(2))
         
         call mat_no_of_matmuls(matmul1)
         CALL DOPT_get_density(config, fifoqueue, queue, F(1), H1, D(1), iteration,ls)
         call mat_no_of_matmuls(matmul2)
         WRITE(config%LUPRI,'("No. of matmuls in get_density: ",I5)') matmul2-matmul1
         CALL LSTIMER('G_DENS',TIMSTR,TIMEND,config%LUPRI)
         
         ls%setting%scheme%DALINK = dalink     !Turn DaLink back on, if requested:
         ls%setting%scheme%DFT%CS00 = CS00
      endif
   endif

!
! SCF iterations
!
   DO iteration = istart, config%opt%cfg_max_linscf_iterations
      CALL LSTIMER('START ',t1,t2,config%lupri)
!     Incremental scheme set for the density-fitting gradient contribution. /SR 2010-10-19
      incremental = config%opt%cfg_incremental .AND. iteration.NE. 1
      call II_setIncremental(ls%setting%scheme,incremental)
      WRITE(config%LUPRI,'("** Get Fock matrix number ",i3)') iteration
      CALL LSTIMER('START ',TIMSTR,TIMEND,config%lupri)
      if (iteration == 1) then
         dalink = ls%setting%scheme%DALINK
         ls%setting%scheme%DALINK = .false.
         cs00 = ls%setting%scheme%DFT%CS00
         ls%setting%scheme%DFT%CS00 = .FALSE.
         CALL get_fock(config, fifoqueue, queue, iteration,D,H1,F,ndmat,E,ls)         
         ls%setting%scheme%DALINK = dalink     !Turn DaLink back on, if requested:
         ls%setting%scheme%DFT%CS00 = CS00
      else
         if (iteration == istart+1) then
            IF(config%opt%crashcalc)THEN
               print*,'Calculation was intentionally crashed due to keyword .CRASHCALC'
               print*,'This keyword is only used for debug and testing purposes'
               print*,'We want to be able to test the .RESTART keyword'
               WRITE(config%lupri,*)'Calculation was intentionally crashed due to keyword .CRASHCALC'
               WRITE(config%lupri,*)'This keyword is only used for debug and testing purposes'
               WRITE(config%lupri,*)'We want to be able to test the .RESTART keyword'
               call lsquit('Crash SCF calculation due to keyword .CRASHCALC',config%lupri)
            ENDIF
         endif
 !        ! START DEBUG: PATRICK
 !         call mat_init(Dtest,D(1)%nrow,D(1)%ncol)
 !         call mat_assign(Dtest,D(1))

 !        h=1E-5_realk
 !         dalink = ls%setting%scheme%DALINK
 !         ls%setting%scheme%DALINK = .false.
 !         cs00 = ls%setting%scheme%DFT%CS00
 !         ls%setting%scheme%DFT%CS00 = .FALSE.
 !         D(1)%elms(1)= D(1)%elms(1)-h
 !      WRITE(config%LUPRI,*) "** DEBUG ADMM - PRINT DENSITY MAT for D(1,1)min"
 !        call mat_print(D(1),1,D(1)%nrow,1,D(1)%ncol,config%lupri)
 !         CALL get_fock(config, fifoqueue, queue, iteration,D,H1,F,ndmat,E,ls)         
 !         call mat_assign(D(1),Dtest)
 !         ls%setting%scheme%DALINK = dalink     !Turn DaLink back on, if requested:
 !         ls%setting%scheme%DFT%CS00 = CS00
 !      write(config%lupri,*) "** DEBUG ADMM - PRINT ENERGY for D(1,1)min"
 !      write(config%lupri,*) "Emin", E(1)
 !      WRITE(config%LUPRI,*) "** DEBUG ADMM - PRINT FOCK MAT for D(1,1)min"
 !        call mat_print(F(1),1,F(1)%nrow,1,F(1)%ncol,config%lupri)


 !         dalink = ls%setting%scheme%DALINK
 !         ls%setting%scheme%DALINK = .false.
 !         cs00 = ls%setting%scheme%DFT%CS00
 !         ls%setting%scheme%DFT%CS00 = .FALSE.
 !         D(1)%elms(1)= D(1)%elms(1)+1E0_realk*h
 !         CALL get_fock(config, fifoqueue, queue, iteration,D,H1,F,ndmat,E,ls)         
 !         call mat_assign(D(1),Dtest)
 !         ls%setting%scheme%DALINK = dalink     !Turn DaLink back on, if requested:
 !         ls%setting%scheme%DFT%CS00 = CS00
 !      write(config%lupri,*) "** DEBUG ADMM - PRINT ENERGY for D(1,1)max"
 !      write(config%lupri,*) "Emax",E(1)
 !      WRITE(config%LUPRI,*) "** DEBUG ADMM - PRINT DENSITY MAT for D(1,1)max"
 !        call mat_print(D(1),1,D(1)%nrow,1,D(1)%ncol,config%lupri)
 !      WRITE(config%LUPRI,*) "** DEBUG ADMM - PRINT FOCK MAT for D(1,1)max"
 !        call mat_print(F(1),1,F(1)%nrow,1,F(1)%ncol,config%lupri)

 !  ! back to original density
 !         CALL get_fock(config, fifoqueue, queue, iteration,D,H1,F,ndmat,E,ls)

 !      write(config%lupri,*) "** DEBUG ADMM - PRINT ENERGY for D0"
 !      write(config%lupri,*) "Eexact",E(1)
 !      WRITE(config%LUPRI,*) "** DEBUG ADMM - PRINT DENSITY MAT for D0"
 !        call mat_print(D(1),1,D(1)%nrow,1,D(1)%ncol,config%lupri)
 !      WRITE(config%LUPRI,*) "** DEBUG ADMM - PRINT FOCK MAT for D0"
 !        call mat_print(F(1),1,F(1)%nrow,1,F(1)%ncol,config%lupri)

 !         call mat_free(Dtest)
 ! call lsquit('Stop SCF after first iteration for DEBUGGING ADMM',-1)
 !        ! END  DEBUG: PATRICK

         CALL get_fock(config, fifoqueue, queue, iteration,D,H1,F,ndmat,E,ls)
      endif
      CALL LSTIMER('FCK_FO ',TIMSTR,TIMEND,config%LUPRI)

!      if (config%solver%step_accepted)Then
         IF(config%opt%cfg_oao_gradnrm)THEN
            call get_oao_transformed_matrices(config%decomp,F(1),D(1)) 
            call mat_init(grad,nbast,nbast)
            gradalloc = .TRUE.
            call get_OAO_gradient(config%decomp%FU, config%decomp%DU, grad) !wrk = gradient
            call mat_scal(0.25E0_realk,grad) !To match linear transformation, also divided by 4!
            gradnrm = sqrt(mat_sqnorm2(grad))         !gradnrm in OAO
            CALL get_AO_gradient(F(1), D(1), S, grad) !grad in AO grad
         ELSE
            call mat_init(grad,nbast,nbast)
            gradalloc = .TRUE.
            CALL get_AO_gradient(F(1), D(1), S, grad) !grad in AO grad
            gradnrm = sqrt(mat_sqnorm2(grad))         !gradnrm in OAO
         ENDIF
!!$      else
!!$            !in principel the printet gradient norm is wrong
!!$            !but the calculation of the gradient 
!!$            !(the get_oao_transformed_matrices call) 
!!$            !would modify config%decomp%FU and config%decomp%DU which 
!!$            !is not desirerable when we want to revert back to old D
!!$         IF(.NOT.gradalloc)call mat_init(grad,nbast,nbast)
!!$      endif
      CALL LSTIMER('G_GRAD',TIMSTR,TIMEND,config%LUPRI)

      ! Statistic stuff
      call scf_stats_update(iteration,gradnrm,E(1),config%opt)
      call scf_stats_debug_mem(config%lupri,iteration)
      !if (DEBUG_DCHANGE) call debug_dchange_update(D,S)

      ! Test for convergence
      IF(gradnrm < config%opt%set_convergence_threshold) then 
         ! Write dens.restart in case *someone* is using DEC (which postulate the existence of
	 ! dens.restart) and a molecule whos density is converged already at starting guess 
	 if (iteration==1) then
            denslun = -1
            IF(config%opt%optlevel.EQ.2)THEN
               call lsopen(denslun,'vdens.restart','UNKNOWN','UNFORMATTED')
            ELSE
               call lsopen(denslun,'dens.restart','UNKNOWN','UNFORMATTED')
            ENDIF
            rewind denslun
            call mat_write_to_disk(denslun,D(1),OnMaster)
            call mat_write_info_to_disk(denslun,config%decomp%cfg_gcbasis)
            call lsclose(denslun,'KEEP')
         endif
         if(.NOT.NotLastSCFLevel)THEN
            Write(config%lupri,'(A,i6,A)')'SCF converged in ',iteration,' iterations'
            EXIT
         endif
      else if ((config%opt%cfg_hesonly .or. config%opt%cfg_diaghesonly) .and. &
            & iteration == 2) then
         write(config%lupri,*) 'You have requested only calculation of first SCF energy'
         write(config%lupri,*) '- now moving on to calculate lowest Hessian eigenvalue for this point.'
         exit
      ENDIF
      WRITE(config%LUPRI,'("** Make average of the last F and D matrices")')
      CALL Density_subspace_minimization(config, fifoqueue, queue, E(1), S, H1, grad, F(1), D(1), iteration)
      call mat_free(grad)
      gradalloc = .FALSE.
      CALL LSTIMER('AVERAG',TIMSTR,TIMEND,config%LUPRI)

      WRITE(config%LUPRI,'("** Get new density ")')
      call mat_no_of_matmuls(matmul1)
      CALL DOPT_get_density(config, fifoqueue, queue, F(1), H1, D(1), iteration,ls) 
      call mat_no_of_matmuls(matmul2)
      WRITE(config%LUPRI,'("No. of matmuls in get_density: ",I5)') matmul2-matmul1
      CALL LSTIMER('G_DENS',TIMSTR,TIMEND,config%LUPRI)

      CALL LSTIMER('SCF iteration',t1,t2,config%lupri,ForcePrint)

      IF(NotLastSCFLevel)THEN
         IF(gradnrm < config%opt%set_convergence_threshold) then 
            Write(config%lupri,'(A,i6,A)')'SCF converged in ',iteration,' iterations'
            EXIT
         ENDIF
      ENDIF
   END DO
   CALL LSTIMER('**ITER',TSTR,TEN,config%LUPRI,ForcePrint)
   IF(gradalloc)THEN
      call mat_free(grad)
   ENDIF

   IF(config%opt%cfg_saveF0andD0)THEN
      call ks_free_linesearch_fock()
   ENDIF

   config%solver%set_do_2nd_order = .false.
   call mat_no_of_matmuls(matmultot)
   write(config%lupri,*)
   WRITE(config%LUPRI,'("Total no. of matmuls in SCF optimization: ",I10)') matmultot
   call scf_afterplay(config,H1,S,D(1),E(1),F(1))

   if (config%solver%cfg_2nd_order_all) then
      config%solver%set_do_2nd_order = config%solver%cfg_do_2nd_order
   endif
!
! FREE SCF OPTIMIZATION STUFF (USEUALLY DONE AFTER SCF_AFTERPLAY)
!
  config%solver%OAO_gradnrm_exist=.false. ! reset
 
  call DOPT_get_density_SHUTDOWN(config%diag)
  if (config%opt%cfg_density_method == config%opt%cfg_f2d_arh) then
     call modfifo_free(fifoqueue)
  !endif
  else
     call queue_free(config%av,queue)
  endif
  call scf_stats_shutdown
!!!!Moved outside scfloop, needed also for lcv localization
!#if 0
!   if (cfg_density_method == cfg_f2d_direct_dens .or. cfg_density_method == cfg_f2d_arh .or. &
!       & cfg_check_converged_solution .or. cfg_rsp_nexcit > 0) then   
!     call decomp_shutdown
!     call dd_shutdown
!  endif
!#endif
nullify(config%davidSCF%arh)
nullify(config%davidSCF%decomp)
!should this not be deallocated?
nullify(config%davidSCF%fifoqueue)


END SUBROUTINE scfloop

!> \brief Wrapper for Fock/KS matrix.
!> \author L. Thogersen
!> \date 2002
!>
!> Get the new fock-matrix F(D). If we already evaluated one in densopt because
!> of a configuration-shift test, this one is used.
!>
subroutine get_fock(config,fifoqueue,queue,iteration,D,H1,F,ndmat,Etotal,ls)
   IMPLICIT NONE
   !> Contains all info about configuration/settings for SCF calculation
   type(configItem),intent(inout)         :: config
   !> New queue type: Contains Fock/KS and density matrices from previous SCF iterations (if ARH)
   type(modFIFO), intent(inout)                  :: fifoqueue
   !> Old queue type: Contains Fock/KS and density matrices from previous SCF iterations (if DIIS)
   type(util_HistoryStore), intent(inout) :: queue
   !> Current SCF iteration
   integer, intent(in)                    :: iteration
   !> number of density matrices
   integer, intent(in)                    :: ndmat
   !> Current density matrix
   TYPE(Matrix),intent(inout),target      :: D(ndmat)
   !> One-electron Hamiltonian
   TYPE(Matrix),intent(in)                :: H1
   !> Fock/KS matrix
   TYPE(Matrix),intent(inout)             :: F(ndmat)
   !> SCF energy corresponding to constructed Fock/KS matrix 
   real(realk), INTENT(OUT)               :: Etotal(ndmat)
   !> Contains settings for integral code
   type(lsitem),intent(inout)             :: ls
   type(matrix)                :: wrk
   integer :: queue_lu, ndim
   logical, external :: do_dft
   integer :: previous
   real(realk) :: r  ! ratio for quadraticity of trust region
   ls%setting%scheme%DFT%CS00eHOMO = config%diag%eHOMO
   ndim = F(1)%nrow
   if (config%diag%cfg_no_confs_checked_in_rh) then 
      if (config%opt%cfg_queue_on_disk) then
         !Possibility to dump queue to disk while contruction Fock matrix for saving memory
         if (config%opt%cfg_density_method == config%opt%cfg_f2d_arh) then
            call fifoqueue_on_disk(fifoqueue,queue_lu,ndim)
         else
            call queue_on_disk(queue,queue_lu,ndim)
         endif
      endif
      CALL di_get_fock_LSDALTON(D,h1,F,ndmat,Etotal,config%lupri,config%luerr,ls)
      if (config%opt%cfg_queue_on_disk) then
         !Restore queue if it has been dumped to disk
         if (config%opt%cfg_density_method == config%opt%cfg_f2d_arh) then
            call fifoqueue_from_disk(fifoqueue,queue_lu,ndim)
         else
            call queue_from_disk(queue,queue_lu,ndim)
         endif
      endif
      if ((config%opt%cfg_density_method == config%opt%cfg_f2d_arh .or. &
           config%solver%cfg_do_2nd_order) .and. iteration > 1) then
         if (config%davidSCF%arh_davidson) then
            if (Etotal(ndmat) - config%solver%old_energy < 0) then
               if (config%davidSCF%arh_linesearch) then
                  config%davidSCF%ActualEnergyDiff = ABS(Etotal(ndmat)-config%davidSCF%arh_linesE)
                  config%davidSCF%EnergyDiffset = .TRUE.
                  IF(config%davidSCF%arh%xnorm < 9.0E-5_realk) then 
                     !maybe test on gradnrm < XFACTOR*config%opt%set_convergence_threshold
                     !close to convergence deactivating ARH LINESEARCH
                     !deactivating ARH LINESEARCH due to problems with integral 
                     !accuracy the energy differences between the different densities 
                     !are now soo small that it is basicly within the nummerical accuracy of
                     !machine precision.
                     config%solver%step_accepted=.TRUE.
                     config%davidSCF%arh_linesearch=.false.
                     write(config%lupri,*)'ARHLS: ARH LINESEARCH TURNED OFF!!'
!                     write(config%lupri,*)'ARHLS: xnorm',config%davidSCF%arh%xnorm
!                     write(config%lupri,*)'ARHLS: config%opt%set_convergence_threshold',config%opt%set_convergence_threshold
                  ELSE !default not close to convergence
                     if (config%davidSCF%MaxLineSearchEnergyDiff > ABS(Etotal(ndmat)-config%davidSCF%arh_linesE)) then
                        !default no error in linesearch
                        config%solver%step_accepted=.true.
                        config%davidSCF%stepsize=min(2.5_realk*config%davidSCF%stepsize,&
                             &config%davidSCF%max_stepsize)
!                        write(config%lupri,*)'ARHLS: xnorm',config%davidSCF%arh%xnorm
!                        WRITE(config%lupri,'(A,ES22.13)')'ARHLS: actual Etotal from Fock Matrix  ',Etotal(ndmat)
!                        WRITE(config%lupri,'(A,ES22.13)')'ARHLS: linesearch Energy Estimate      ',config%davidSCF%arh_linesE
                        WRITE(config%lupri,'(A,ES22.13)')'ARHLS: Energy difference between estimate and actual Energy',&
                             & ABS(Etotal(ndmat)-config%davidSCF%arh_linesE)
                        WRITE(config%lupri,'(A,ES22.13)')'ARHLS: Maximum difference in linesearch energies',&
                             & config%davidSCF%MaxLineSearchEnergyDiff
                        IF(ABS(config%davidSCF%MaxLineSearchEnergyDiff) < 1.0E-10_realk)THEN
                           !turn off incremental build in linesearch, to improve accuracy
                           IF(config%opt%cfg_saveF0andD0)THEN
                              WRITE(config%lupri,*)'ARHLS: Turn off incremental linesearch energy estimates to improve accuracy'
                              saveF0andD0=.FALSE.                           
                              call ks_free_linesearch_fock()
                              config%opt%cfg_saveF0andD0 = .FALSE.
                           ENDIF
                        ENDIF
                        !
                     else
                        !possible error in linesearch
                        IF(ABS(config%davidSCF%MaxLineSearchEnergyDiff).LT.ls%setting%scheme%THRESHOLD)THEN
                           config%solver%step_accepted=.TRUE.
                           config%davidSCF%arh_linesearch=.false.
!                           write(config%lupri,'(A,ES24.12)')'ARHLS: actual Etotal   ',Etotal(ndmat)
!                           write(config%lupri,'(A,ES24.12)')'ARHLS: linesearch E    ',config%davidSCF%arh_linesE
                           write(config%lupri,'(A,ES24.12)')'ARHLS: actual Ediff    ',ABS(Etotal(ndmat)-config%davidSCF%arh_linesE)
                           write(config%lupri,'(A,ES24.12)')'ARHLS: Ediff linesearch',config%davidSCF%MaxLineSearchEnergyDiff
                           write(config%lupri,*)'ARHLS: differences between linesearch energies are too small'
                           write(config%lupri,*)'ARHLS: ARH LINESEARCH TURNED OFF!!'
                        ELSE
                           print*,'ARHLS: warning possible error in linesearch'
                           write(config%lupri,'(A)')'ARHLS: warning possible error in linesearch'
                           
                           print*,'ARHLS: actual Etotal   ',Etotal(ndmat)
                           print*,'ARHLS: linesearch E    ',config%davidSCF%arh_linesE
                           print*,'ARHLS: actual Ediff    ',ABS(Etotal(ndmat)-config%davidSCF%arh_linesE)
                           print*,'ARHLS: Ediff linesearch',config%davidSCF%MaxLineSearchEnergyDiff
                           write(config%lupri,'(A,ES24.12)')'ARHLS: actual Etotal   ',Etotal(ndmat)
                           write(config%lupri,'(A,ES24.12)')'ARHLS: linesearch E    ',config%davidSCF%arh_linesE
                           write(config%lupri,'(A,ES24.12)')'ARHLS: actual Ediff    ',ABS(Etotal(ndmat)-config%davidSCF%arh_linesE)
                           write(config%lupri,'(A,ES24.12)')'ARHLS: Ediff linesearch',config%davidSCF%MaxLineSearchEnergyDiff
                           write(config%lupri,'(a,ES13.4,ES13.4)') 'ARHLS: ARH LINESEARCH: Ediff linesearch vs. actual Ediff :', &
                                &config%davidSCF%MaxLineSearchEnergyDiff,  ABS(Etotal(ndmat)-config%davidSCF%arh_linesE)
!                           WRITE(config%lupri,*)'ARHLS: If the energy changes are very small you may face an issue '
!                           WRITE(config%lupri,*)'ARHLS: with accuracy and we suggest to use the .NOECONTINCREM keyword'
!                           WRITE(config%lupri,*)'ARHLS: and possible the .NLSDASCREENO keyword'
!                           WRITE(config%lupri,*)'ARHLS: .NOECONTINCREM keyword deactivates the calculation of energies exploiting'
!                           WRITE(config%lupri,*)'ARHLS: a incremental (recursive) Energy construction'
!                           WRITE(config%lupri,*)'ARHLS: .NOLSDASCREEN keyword deactivates screening directly on the '
!                           WRITE(config%lupri,*)'ARHLS: Energy and reverts back to screening on Fock matrix elements'
!                           print*,'ARHLS: If the energy changes are very small you may face an issue '
!                           print*,'ARHLS: with accuracy and we suggest to use the .NOECONTINCREM keyword'
!                           print*,'ARHLS: and possible the .NOLSDASCREEN keyword'
!                           print*,'ARHLS: .NOECONTINCREM keyword deactivates the calculation of energies exploiting'
!                           print*,'ARHLS: a incremental (recursive) Energy construction'
!                           print*,'ARHLS: .NOLSDASCREEN keyword deactivates screening directly on the '
!                           print*,'ARHLS: Energy and reverts back to screening on Fock matrix elements'

                           !we accept the step because the step does have a lower energy 
                           !the linesearch maybe wrong or non optimal but the result is fine
                           config%solver%step_accepted=.TRUE.
                           config%davidSCF%arh_linesearch=.false.
                           write(config%lupri,*)'ARHLS: ARH LINESEARCH TURNED OFF!!'
                        ENDIF
                     endif 
                  ENDIF !turn off ARH 
               else
                  !
                  config%solver%step_accepted=.true.
                  config%davidSCF%stepsize=min(2.5_realk*config%davidSCF%stepsize,&
                       &config%davidSCF%max_stepsize)
               endif !config%davidSCF%arh_linesearch
            else               
               !not accepted step
               config%solver%step_accepted=.false.
               config%davidSCF%stepsize=0.5_realk*config%davidSCF%stepsize
            endif !Etotal - config%solver%old_energy < 0
         else
            !not davidson
            call update_trustradius(config%solver, ls, iteration, Etotal(ndmat), fifoqueue%offset)
         end if
         if (.not. config%solver%step_accepted)then
            Etotal(ndmat) = config%solver%old_energy
         endif
      endif
   else
      IF(ndmat.GT.1)call lsquit('Error in get_fock: ndmat.gt.1 special case',-1)
      WRITE(config%lupri,'("** Fock matrix was already found in RH-step '// &
           & 'exploring a configuration shift ")')
      call mat_assign(F(1),queue%F(queue%current_position))
      call mat_assign(D(1),queue%D(queue%current_position))
      Etotal(1) = queue%Energy(queue%current_position)
   endif
   !WRITE(LUPRI,*) 'E_SCF right after evaluation: ',Etotal
   !endif
end subroutine get_fock

!> @file
!> Contains main SCF driver and some wrappers to starting guess, Fock/KS matrix etc.



!> \brief After SCF opt., calculate HOMO-LUMO gap, lowest Hes. eigenvalue, print statistic 'n'stuff.
!> \author L. Thogersen
!> \date 2003
subroutine scf_afterplay(config,H1,S,D,E,F)
   implicit none
   !> Contains all info about configuration/settings for SCF calculation
   type(configItem), intent(inout) :: config
   !> One-electron Hamiltonian
   type(matrix), intent(in) :: H1
   !> Overlap matrix
   type(matrix), intent(in) :: S
   !> Converged AO density matrix
   type(matrix), intent(in) :: D
   !> Converged SCF energy
   real(realk), intent(in) :: E
   !> Converged Fock/KS matrix. Output only if config%opt%cfg_scale_virt = .true.
   type(Matrix), intent(inout) :: F
   type(debugItem) :: debug
   type(DDitem)    :: DD

   !write (lupri,*) 'Incoming D, scf_afterplay:'
   !call MAT_PRINT(D, 1, D%nrow, 1, D%ncol, LUPRI)
   !write (lupri,*) 'Incoming F, scf_afterplay:'
   !call MAT_PRINT(F, 1, D%nrow, 1, D%ncol, LUPRI)
   if (config%decomp%cfg_rsp_nexcit > 0) then
      write(config%lupri,*) 'Postponing calculation of HOMO-LUMO gap to response part...'
   else if (config%opt%cfg_density_method == config%opt%cfg_f2d_direct_dens .or. &
          & config%decomp%cfg_check_converged_solution .or. &
          & config%opt%cfg_density_method == config%opt%cfg_f2d_arh .or. &
          & (config%opt%CFG_density_method == config%opt%cfg_f2d_roothaan .AND. .NOT.&
          & (config%opt%cfg_start_guess == 'TRILEVEL' .AND. config%opt%optlevel == 2 ))) then
      !debug_arh_hessian = .false.
      call get_oao_transformed_matrices(config%decomp,F,D)
      call dd_init(config%decomp,DD)
      if (.not.config%diag%nofinalhomolumo) then
         if (config%diag%cfg_unres) then
            call dd_homolumo_and_heseigen_unres(DD,config%decomp,debug,.false.,0)
         else
           call dd_homolumo_and_heseigen(DD,config%decomp,debug,.false.,0)
         endif
         if (config%opt%cfg_diaghesonly .or. config%opt%debug_diag_hessian) then 
            call debug_diag_full_hessian(config%solver,config%decomp)
         endif
      endif
      call dd_shutdown(config%decomp,DD)
   endif
   if (config%opt%cfg_scale_virt) then
     !Fock matrix is modified - unmodify it
     config%opt%cfg_scale_virt = .false.
     call fck_unscale_virt(H1,S,D,F,config%decomp%nocc)
   endif
   CALL scf_stats_end_print(config%opt)

end subroutine scf_afterplay

!> \brief This subroutine calculates number of electron for neutral molecule defined in ls%setting structure
!> \author T. Kjaergaard, S. Reine
!> \date 2008-10-26
!> \param ls Contains information read from LSDALTON.INP 
!> \param nel Number of electrons
SUBROUTINE get_num_electrons_neutral(nel,ls)
  implicit none
  TYPE(lsitem) :: ls
  integer, intent(out)        :: nel
  integer :: i

  nel=0
  do i=1,ls%setting%MOLECULE(1)%p%nAtoms
    nel = nel + ls%setting%MOLECULE(1)%p%Atom(i)%Charge
  enddo
END SUBROUTINE get_num_electrons_neutral

end module scfloop_module
