!> @file
!> Contains density subspace minimization module.

!> \brief Branches out to the chosen method for density subspace minimization.
!> \author L. Thogersen
!> \date 2003
module density_subspace_min
!  use ARHmodule
  use configurationType,only: configitem
  use matrix_module,only: matrix
  use av_utilities!, only: util_HistoryStore, add_to_queue, flush_queue
  use queue_module,only: modFIFO
  use queue_ops, only: add_to_modfifo
  use precision  
  use decompMod
  use matrix_operations
  use matrix_operations_aux
  use matrix_util,only: dumpmats
  USE linsca_dsm
  USE linsca_diis
!  USE linsca_ediis
  use rsp_util, only: util_save_moinfo, util_free_mostuff, util_ao_to_mo
  use diagonalization, only: SCF_iteration, diag_lshift_none
  use ARHmodule, only: arh_debug_print
CONTAINS
!> \brief Branches out to the chosen method for density subspace minimization.
!>
!> Computes an improved new potential (Fock matrix)
!> taking previous iterations into account (or not, depending on the
!> method). \n
!> Both the new (type(modFIFO)) and old (type(util_HistoryStore)) queue types are passed as
!> inputs, because some of the old code uses the old and impractical type(util_HistoryStore)
!> for matrix storage. Preferably, the old code should be rewritten to use the new
!> queue type. This however, involves a huge amount of work and is not likely to
!> happen.
!>  
!> \param config Contains all info about configuration/settings for SCF calculation
!> \param fifoqueue New queue type (holds F,D matrices if ARH is optimization method, otherwise empty)
!> \param queue Old queue type (holds F,D matrices if ARH is not optimization method, otherwise empty)
!> \param E SCF energy corresponding to input F, D 
!> \param S Overlap matrix
!> \param H1 One-electron Hamiltonian
!> \param grad SCF gradient corresponding to input F, D 
!> \param F Fock/KS matrix to be put in queue
!> \param D Density matrix to be put in queue
!> \param iteration Current SCF iteration
   SUBROUTINE Density_subspace_minimization(config,fifoqueue, queue, E, S, H1,grad, F, D, iteration)
      IMPLICIT NONE
      type(configItem),intent(inout) :: config
      TYPE(util_HistoryStore)  :: queue
      type(modFIFO) :: fifoqueue
      TYPE(Matrix), INTENT(IN) :: S, H1,grad
      TYPE(Matrix)                :: F, D
      real(realk), intent(in)  :: E
      type(Matrix) :: diff, Fmo, Fov
      logical :: restart
      integer, intent(in)      :: iteration
      integer                  :: ndim,nnz
      TYPE(Matrix) :: tmp1, tmp2
      real(realk)  :: maxelm

      ndim = F%nrow

      if (config%solver%debug_dd) then
         if (SCF_iteration > config%solver%cfg_nits_debug) call arh_debug_print(config%solver)
      endif
      if (config%solver%cfg_2nd_order_local .and. config%solver%set_do_2nd_order) then
         write(config%lupri,*) 'Switching to second order optimization - turning off averaging!'
         config%av%cfg_averaging = config%av%cfg_avg_none
      endif

       if (config%opt%cfg_density_method == config%opt%cfg_f2d_arh) then
         call mat_init(tmp1,ndim,ndim)
         call mat_init(tmp2,ndim,ndim)
         if (iteration == 1) then    !If initial D, transform to cholesky basis and add
            config%solver%old_energy = E
            call mat_report_sparsity(D,'AO D ',nnz,config%lupri)
            call mat_report_sparsity(F,'AO F ',nnz,config%lupri) 
            
            call x_to_oao_basis(config%decomp, D, tmp1)
            call res_to_oao_basis(config%decomp, F, tmp2)
            
            call mat_report_sparsity(tmp1,'OAO D',nnz,config%lupri)
            call mat_report_sparsity(tmp2,'OAO F',nnz,config%lupri) 
            
            call add_to_modFIFO(fifoqueue, tmp2, tmp1, E, .true.)
            
            if (config%opt%dumpmatrices) call dumpmats(iteration,tmp1,tmp2,config%opt%optlevel)
          else 
            if (config%solver%step_accepted) then !Only add to queue if it is an accepted step!
               config%solver%Nrejections = 0
               config%solver%old_energy = E
               call mat_report_sparsity(D,'AO D ',nnz,config%lupri)
               call mat_report_sparsity(F,'AO F ',nnz,config%lupri) 
               
               call x_to_oao_basis(config%decomp, D, tmp1)
               call res_to_oao_basis(config%decomp, F, tmp2)
               
               call mat_report_sparsity(tmp1,'OAO D',nnz,config%lupri)
               call mat_report_sparsity(tmp2,'OAO F',nnz,config%lupri) 
               
               call add_to_modFIFO(fifoqueue, tmp2, tmp1, E, .true.)
               
               if (config%opt%dumpmatrices) call dumpmats(iteration,tmp1,tmp2,config%opt%optlevel)
            else
               !FIXME: Remove unnecessary ao-oao transformations!
               call res_from_oao_basis(config%decomp,fifoqueue%F_exp,F)
               call x_from_oao_basis(config%decomp,fifoqueue%D_exp,D)
            endif
         endif
         call mat_free(tmp1)
         call mat_free(tmp2)
      else if (config%solver%cfg_do_2nd_order) then
         if (iteration == 1) then    !If initial D, transform to cholesky basis and add
              config%solver%old_energy = E
              CALL add_to_queue(config%av, F, D, S, E, grad, queue)
          else
            if (config%solver%step_accepted) then !Only add to queue if it is an accepted step!
               config%solver%Nrejections = 0
               config%solver%old_energy = E
               CALL add_to_queue(config%av, F, D, S, E, grad, queue)
            else
               call mat_assign(F,queue%F(queue%current_position))
               call mat_assign(D,queue%F(queue%current_position))
            endif
         endif
      else if (config%av%CFG_averaging == config%av%CFG_AVG_van_lenthe) then
         if (iteration > 1 .and. config%av%vanlentheCounter == 0) then
            !Check size of largest occupied-virtual element in Fock matrix:
            !1. Convert F to MO basis:
            call mat_init(Fmo, F%nrow, F%ncol)
            call util_AO_to_MO(S,F,Fmo,.true.)
            write(config%lupri,*) 'config%decomp%nocc, nbast:', config%decomp%nocc, F%nrow
            !write(config%lupri,*) 'Fmo:'
            !call mat_print(Fmo, 1, F%nrow, 1, F%ncol, lupri)
            !2. Get occupied-virtual block of Fock matrix:
            call mat_init(Fov, config%decomp%nocc, F%nrow-config%decomp%nocc)
            call mat_section(Fmo,1,config%decomp%nocc,config%decomp%nocc+1,F%nrow,Fov)
            call mat_free(Fmo)
            !write(config%lupri,*) 'Fov:'
            !call mat_print(Fov, 1, Fov%nrow, 1, Fov%ncol, lupri)
            call mat_abs_max_elm(Fov, maxelm)
            write(config%lupri,*) 'Fov maxelm', maxelm
            if (maxelm < 0.1E0_realk) then
               call flush_queue(config%av,S,queue%current_position,0,queue)
               write(config%lupri,*) 'Queue flushed!'
               config%av%vanlentheCounter = 1
            endif
            call mat_free(Fov)
            call util_free_mostuff()
         endif
         if (config%av%vanlentheCounter == 0) then
            call util_save_MOinfo(F,S,config%decomp%nocc)
         endif
         CALL add_to_queue(config%av, F, D, S, E, grad, queue)
      else
         if (config%diag%cfg_no_confs_checked_in_rh) then
           CALL add_to_queue(config%av, F, D, S, E, grad, queue) 
         endif
      endif

      if (config%av%cfg_averaging == config%av%CFG_AVG_NONE) then
         ! nothing
      else if (config%av%cfg_averaging == config%av%CFG_AVG_DIIS) then
         CALL diis(config%av, queue, D, F)
      else if (config%av%cfg_averaging == config%av%CFG_AVG_EDIIS) then
         call lsquit('EDIIS feature was removed',-1)
!         CALL ediis(config%av, queue, D, F)
      else if (config%av%cfg_averaging == config%av%cfg_avg_van_lenthe) then
         if (config%av%vanlentheCounter > 0) then
            config%av%vanlentheCounter = config%av%vanlentheCounter + 1
         endif
         if (config%av%vanlentheCounter == 4) then
            !Turn off level shift, turn on DIIS
            config%diag%CFG_lshift = diag_lshift_none
            config%av%usequeue = .true.
         endif
         if (config%av%usequeue) then
            CALL diis(config%av, queue, D, F)
            write(config%lupri,*) 'Used DIIS averaging'
         else !simple averaging: Fav = 1/2*(F(n) + F(n-1))
            call simple_averaging(config%av, iteration, D, F)
            write(config%lupri,*) 'Used simple averaging'
         endif
         call mat_assign(config%av%Fprev,F)
         call mat_assign(config%av%Dprev,D)
      else
         stop 'Unknown type of averaging'
      endif
     !write(config%lupri,*) 'F after DSM:'
     !call mat_print(F,1,ndim,1,ndim,lupri)
     !write(config%lupri,*) 'D after DSM:'
     !call mat_print(D,1,ndim,1,ndim,lupri)

      !if (DEBUG_DCHANGE .and. cfg_averaging /= CFG_AVG_NONE) then
      !  call mat_init(diff,S%nrow,S%ncol)
      !  call mat_add(1E0_realk,D,-1E0_realk,queue%D(queue%current_position),diff)
      !  Ddiff(stat_current_iteration,2) = util_Snorm(diff,S)
      !  call mat_free(diff)
      !endif
   END SUBROUTINE Density_subspace_minimization

end module density_subspace_min
