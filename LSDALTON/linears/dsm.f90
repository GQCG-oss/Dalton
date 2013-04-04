!> @file 
!> Contains DSM module.

!> Density Subspace minimization (DSM) module. JCP, 121, 16-27, 2004; JCP 123, 074103, 2005
module LINSCA_DSM
   use precision
   use av_utilities
   use dsm_x_term
   use matrix_util, only: util_Snorm
   use matrix_module
   use matrix_operations
   private
   public :: dsm_get_avg_density
   !logical, parameter :: former = .false.
   !> Eigenvalues smaller than this are considered negative!
   real(realk), save :: dsm_hes_zero = 0.0E0_realk
   !> Size limit on accepted step
   real(realk), save :: changelimit = 1.2E0_realk
   !real(realk), save :: timing1, timing2
   !> Keep track of current position etc. for subspace of densities
   integer, allocatable, save :: his_p(:)

CONTAINS                                      
   SUBROUTINE dsm_init
   END SUBROUTINE dsm_init

!> \brief Construct average F and D from subspace from previous SCF iterations.
!> \author L. Thogersen
!> \date 2003
!> \param av Contains info about SCF averaging
!> \param queue Subspace of Fock/KS and density matrices from previous SCF iterations
!> \param S Overlap matrix
!> \param H1 One-electron Hamiltonian
!> \param F The averaged Fock/KS matrix
!> \param D The averaged density matrix
!> \param restart True is DSM should be restarted, e.g. after removal of one or more vectors
   SUBROUTINE dsm_get_avg_density(av, queue, S, H1, F, D,restart)
!      use fock_evaluator
      IMPLICIT NONE
      type(avItem), intent(inout) :: av
      TYPE(util_HistoryStore)  :: queue
      type(Matrix), intent(in) :: S,H1
      type(Matrix) :: F,D   !output
      logical, intent(out) :: restart
      type(Matrix) :: Dbar,Co
      real(realk),allocatable :: weights(:)
      real(realk) :: EHF
      integer :: i,msize

      changelimit = 1.2E0_realk
      restart = .false.
      !print *, "dsm:queue%used_entries", queue%used_entries
      if (queue%used_entries == 1) then
         call mat_assign(D,queue%D(1))
         call mat_assign(F,queue%F(1))
      else
         !** Find DSM weights
         if (av%CFG_SET_type == av%CFG_THR_dft) then !HACK TO TEST....
           dsm_hes_zero = 0.1E0_realk
         endif
         if (av%debug_dsm_dhistory) then
           call CLEAN_DHISTORY(av,S,queue)
         endif
         msize = MIN(queue%used_entries,av%dsm_history_size)
         WRITE(av%lupri,*) 'DSM - Vectors in queue:', msize
         ALLOCATE(weights(msize))
         ALLOCATE(his_p(msize))
         call set_his_p(av,msize,queue%current_position) 
         call DSM_WEIGHTS(av,queue,msize, S,H1,weights,restart)
         DEALLOCATE(his_p)
         if (restart) then
           DEALLOCATE(weights)
           return
         endif
         if (av%info_weight_final .or. av%info_weights) then
            WRITE(av%LUPRI,*) 'WEIGHTS,POS',av%dsm_pos
            call OUTPUT(weights,1,msize,1,1,msize,1,1,av%lupri)
         endif
         !** Construct average density and fock matrix
         if (av%cfg_safe) then
           call lsquit('cfg_safe diabled! /Stinne',av%lupri)
           !call mat_init(Dbar,S%nrow,S%ncol)
           !call get_AVERAGE_arr(av,'D',queue,av%dsm_history_size,weights,Dbar)
           !!** diagonalize SDS (SDS C = SCe) giving the corresponding
           !! idempotent D = <C|C>
           !call GET_Didem(Dbar,S,D)
           !!** Find the corresponding Fock matrix
           !call fck_get_fock(D,F,EHF)
!          ! !** Overwrite the latest F and D in the queue with these
!          ! call retype_last_in_queue(F,D,S,EHF, queue)
           !call mat_free(Dbar)
         else
           call get_AVERAGE_arr(av,'D',queue,av%dsm_history_size,weights,D)
           !call PURIFY(Ndim,Dav,S) !If new D is not found diagonalizing F
           call get_AVERAGE_arr(av,'F',queue,av%dsm_history_size,weights,F)
           if (av%cfg_dsm_app == av%cfg_dsm_xtra_term) then
             !** + corrections Delta_parr
             call mat_init(Dbar,S%nrow,S%ncol)
             call eval_delta(D,S,Dbar)
             call util_GET_PROJ_PART(av,queue,Dbar,S,weights)
             call get_average_arr(av,'D',queue, av%dsm_history_size, weights, Dbar)
             call mat_daxpy(1E0_realk,Dbar,D)
             call get_average_arr(av,'F',queue, av%dsm_history_size, weights, Dbar)
             call mat_daxpy(1E0_realk,Dbar,F)
             call mat_free(Dbar)
           endif
         endif
         DEALLOCATE(weights)
      endif

   end subroutine DSM_GET_AVG_DENSITY

!> \brief Get weights for averaging Fock/KS and density matrices in DSM.
!> \author L. Thogersen
!> \date 2003
!> \param av Contains info about SCF averaging
!> \param queue Subspace of Fock/KS and density matrices from previous SCF iterations
!> \param msize Size of subspace
!> \param S Overlap matrix
!> \param H1 One-electron Hamiltonian
!> \param weights The weights for averaging Fock/KS and density matrices
!> \param restart True is DSM should be restarted, e.g. after removal of one or more vectors
   subroutine DSM_WEIGHTS(av,queue,msize,S,H1,weights,restart)
      !Find and return the proper DSM weights
      USE scf_stats, only: stat_current_iteration, stat_tab
      implicit none
      type(avItem), intent(inout) :: av
      TYPE(util_HistoryStore)  :: queue
      integer, intent(in)      :: msize 
      type(Matrix), intent(in) :: S, H1
      real(realk), intent(out) :: weights(msize)
      logical, intent(out) :: restart
      integer, parameter :: maxstep = 10
      real(realk) :: &
           & trustr, DSMgrad(msize-1), DSMhes(msize-1,msize-1), mineival(maxstep),&
           & mineivalold,DSMener,coefcor(msize-1),xval,coef(msize-1),&
           & deltac_tot(msize),minabs_eival,DSMgradnorm_0,minEner,S2(msize-1,msize-1),&
           & remove_dir(5,msize),gdtest(msize-1),Newtonstep(maxstep),TrDSDS,TrDavSDavS,&
           & DSMener_start
      type(Matrix) :: Dav,Fav,Dchange,Dtilde,Delta 
      type(Matrix) :: SD,SDS !Poul test
      integer :: i, i1,j, j1,k, nqdsm, naccstep, ierr, minstart, steprej, n_dir
      integer :: pos, ndim, his_start
      real :: t1, t2
      logical :: posdef,outproj,end_dsm,pure_Newton,done

      !** Initializations
      end_dsm = .false.
      nqdsm = 0       !number of dsm iterations
      naccstep = 0
      weights = 0.0E0_realk
      coef = 0.0E0_realk
      trustr = 0.2E0_realk   !normally it works well with 0.2
      steprej = 0
      posdef = .false.
      deltac_tot = 0.0E0_realk
      outproj = .false.
      pure_Newton = .false.
      Newtonstep = 0.0E0_realk
      stat_tab(stat_current_iteration+1,6) = 0
!      pos = queue%current_position
      pos = av%dsm_pos
      ndim = queue%D(queue%current_position)%nrow
      done = .false. !for debug purpose only
      !** Allocations
      call mat_init(Dav,ndim,ndim)
      call mat_init(Fav,ndim,ndim)
      call mat_init(Delta,ndim,ndim)
      call mat_init(SD,S%nrow,S%ncol)
      call mat_init(SDS,S%nrow,S%ncol)
      if (av%cfg_dsm_app == av%cfg_dsm_xtra_term) then
        call xterm_start(av,ndim,msize)
      endif
      !** Make initial weight-guess
      !choosing the one with lowest energy as the starting point
      j = queue%current_position
      minEner = queue%Energy(j)
      i = av%dsm_pos
      minstart = i
      his_start = j
      do k = 1,msize-1
        j = j-1; if(j==0) j=av%cfg_settings(av%cfg_set_type)%max_history_size
        i = i-1; if(i==0) i=av%dsm_history_size
        if (queue%Energy(j) < minEner) then
          minEner = queue%Energy(j)
          minstart = i
          his_start = j
        endif
      enddo
      weights(minstart) = 1.0E0_realk
      if (av%info_weights) then
         WRITE(av%LUPRI,*) 'WEIGHTS, dsm it. ',nqdsm,minEner
         WRITE(av%LUPRI,*) minstart,1E0_realk
      endif
      !create the S2 metric from the TrDSDS array in queue%metric 
      i1 = 0
      do i = 1,msize
        if (i /= minstart) then
          i1 = i1+1
          j1 = 0
          do j = 1,msize
            if (j /= minstart) then
              j1 = j1+1
              S2(i1,j1) = queue%metric(i,j) - queue%metric(i,minstart) &
                     &- queue%metric(minstart,j) + queue%metric(minstart,minstart)
            endif
          enddo
        endif
      enddo
      !Find the initial DSM energy
      call mat_assign(Fav,queue%F(his_start))
      call mat_assign(Dav,queue%D(his_start))
      call eval_delta(Dav,S,Delta)
      if (av%cfg_dsm_app == av%cfg_dsm_xtra_term) then
        call xterm_build_Gpara(av,msize,queue,Delta,S,H1)
      endif
      call DSMenergy(av,msize,his_start,queue,S,H1,Fav,Dav,Delta,DSMener)
      DSMener_start = DSMener   !for use in restart with projections
      if (av%info_dsm_energy) then
         write(av%lupri,*) 'START',' nqdsm',nqdsm,'DSMener',DSMener
      endif
      stat_tab(stat_current_iteration+1,5) = DSMener
      do
         nqdsm = nqdsm + 1
         !** Find the gradient and hessian
         !call CPU_TIME(t1)
         !call DSMgrad_and_hes(msize,his_start,minstart,queue,Delta,S,H1,Dav,Fav,DSMgrad,DSMhes)
         !call CPU_TIME(t2)
         !WRITE(LUPRI,*) 'old routine - ',t2-t1
         !do i = 1,msize-1
         !  WRITE(LUPRI,'(f9.4,4x,8f9.4)') DSMgrad(i),(DSMhes(i,j),j=1,msize-1)
         !enddo
         !WRITE(LUPRI,*)
         !call CPU_TIME(t1)
         call DSMgrad_and_hes_better(msize,his_start,minstart,queue,Delta,S,H1,Dav,Fav,DSMgrad,DSMhes)
         !call CPU_TIME(t2)
         !WRITE(LUPRI,*) 'new routine - ',t2-t1
         !do i = 1,msize-1
         !  WRITE(LUPRI,'(f9.4,4x,8f9.4)') DSMgrad(i),(DSMhes(i,j),j=1,msize-1)
         !enddo
         !WRITE(LUPRI,*)
         !** Analyze the gradient and hessian
         call ANALYZE_G_H(av,nqdsm,msize,queue,S2,naccstep,DSMgrad,DSMhes,&
              &outproj,mineival,minabs_eival,n_dir,remove_dir)
         if (mineival(naccstep+1) > dsm_hes_zero) then
            posdef = .true.
         elseif (posdef .and..not. outproj) then
            !hessian goes from pos.def. to neg.def.
            stat_tab(stat_current_iteration+1,6) = 4E0_realk
            naccstep = naccstep + 1
            posdef = .false.
            call EXIT_OR_RESTART(av,msize,queue,S,maxstep, &
                                &his_start,coefcor,DSMener,DSMener_start,&
                                &end_dsm,steprej,posdef,deltac_tot,&
                                &naccstep,outproj,nqdsm,&
                                &coef,Fav,Dav,mineival,trustr,Newtonstep,restart)
            if (end_dsm) then
               exit
            elseif (restart) then
               call mat_free(Dav)
               call mat_free(Fav)
               call mat_free(SD)
               call mat_free(SDS)
               call mat_free(Delta)
               if (av%cfg_dsm_app == av%cfg_dsm_xtra_term) then
                 call xterm_end()
               endif
               return
            else
               cycle
            endif
         endif
         if (av%info_dsm_cnorm_mu_fig .and. stat_current_iteration+1 < 7) then
            !find the delta c norm for different mu
            call MU_DCNORM(av,msize,queue,S2,DSMgrad,DSMhes,n_dir,remove_dir,trustr) 
         endif
         !** find delta c (coefcor)
         call DSMstep(av,naccstep,msize, queue,trustr, S2, DSMgrad,DSMhes,mineival(naccstep+1),&
              &n_dir,remove_dir,Newtonstep,pure_Newton,coefcor)
         !** checking for convergence and deciding whether the step
         ! should be rejected.
         call STEP_EVAL(av,nqdsm,msize,his_start,n_dir,remove_dir,minstart,queue,DSMgradnorm_0,&
                       &coefcor,S2,DSMgrad,DSMhes,S,H1,&
                       &posdef,pure_Newton,Newtonstep,coef,DSMener,trustr,&
                       &deltac_tot,Fav,Dav,Delta,naccstep,steprej,end_dsm)
         call EXIT_OR_RESTART(av,msize,queue,S,maxstep, &
                                &his_start,coefcor,DSMener,DSMener_start,&
                                &end_dsm,steprej,posdef,deltac_tot,&
                                &naccstep,outproj,nqdsm,&
                                &coef,Fav,Dav,mineival,trustr,Newtonstep,restart)
         if (steprej == 0 .and. (av%cfg_dsm_app == av%cfg_dsm_one .or. &
            & av%cfg_dsm_app == av%cfg_dsm_search .or. av%debug_dsm_linesearch) ) then
           !end after this first real step
           !if cfg_Escf_min a linesearch is made in the direction for the min of Escf
           if (av%cfg_dsm_app == av%cfg_dsm_one .or. av%cfg_dsm_app == av%cfg_dsm_search) then
             call use_first_direction(av,msize,minstart,his_start,queue,S,H1,Delta,coef)
             end_dsm = .true.
           endif
           !if (av%debug_dsm_linesearch .and. stat_current_iteration+1 == 7 .and. .not. done) then
           !  call debug_linesearch(msize,minstart,his_start,queue,S,H1,coef)
           !  done = .true.
           !endif
         endif
         if (end_dsm) then
           exit
         elseif (restart) then
           call mat_free(Dav)
           call mat_free(Fav)
           call mat_free(SD)
           call mat_free(SDS)
           call mat_free(Delta)
           if (av%cfg_dsm_app == av%cfg_dsm_xtra_term) then
             call xterm_end()
           endif
           return
         endif
      enddo
      TrDavSDavS = util_Snorm(Dav,S)
      !call cpu_time(timing1)
      call mat_mul(S,Dav,'n','n',1E0_realk,0E0_realk,SD)
      call mat_mul(SD,S,'n','n',1E0_realk,0E0_realk,SDS)
      !call cpu_time(timing2)
      !cfg_dsmtime = cfg_dsmtime + (timing2 - timing1)
!      TrDSDS = mat_dotproduct(queue%D(queue%current_position),SDS)
!      /stat_tab(stat_current_iteration+1,11) = TrDSDS/sqrt(TrDavSDavS*CFG_NOCC)
      if(av%info_dsm_nit) WRITE(av%LUPRI,*) nqdsm,' iterations made in DSM'
      !create weights from coef
      j = 0
      do i = 1,msize
        if (i == minstart) then
          weights(i) = 1.0E0_realk - SUM(coef)
        else
          j = j + 1
          weights(i) = coef(j)
        endif
      enddo
      !if (av%DEBUG_DSM_DCHANGE) then
      !  call debug_dchange_dsm(msize,minstart,his_start,queue,S,weights,coef,Dav)
      !endif
      !if (av%DEBUG_DSM_EMODEL) then
      !  call debug_emodel(msize,his_start,weights,DSMener_start,queue,S,H1)
      !endif
      call mat_free(Dav)
      call mat_free(Fav)
      call mat_free(Delta)
      call mat_free(SD)
      call mat_free(SDS)
      if (av%cfg_dsm_app == av%cfg_dsm_xtra_term) then
        call xterm_end()
      endif


   end subroutine DSM_WEIGHTS

!> \brief Set up linear equations and solve, finding delta c
   subroutine DSMstep(av,naccstep,msize,queue,trustr, S2, DSMgd,DSMhes,mineival,&
        &             n_dir,remove_dir,Newtonstep,pure_Newton,coefcor)
      implicit none
      type(avitem),intent(in) :: av
      integer, intent(in)     :: naccstep, msize,n_dir
      TYPE(util_HistoryStore) :: queue
      real(realk), intent(in) :: trustr,remove_dir(:,:),mineival
      real(realk), intent(in) :: S2(msize-1,msize-1),DSMgd(msize-1), DSMhes(msize-1,msize-1)
      real(realk), intent(inout) :: Newtonstep(:)
      logical, intent(out)    :: pure_Newton
      real(realk), intent(out) :: coefcor(msize-1)
      !
      real(realk) :: stepminustrustr
      integer :: i
      real(realk) ::  xstep,x1mult,x2mult,xaccmult,fx2,fx1,ds,lowsing
      logical :: STEPFINE

      STEPFINE = .false.   !is the Newton step short enough?
      pure_Newton = .false. 
      if (av%info_dsm_step) then
         !     Newton step - no damping of Hessian
         call SOLVE_LINEQ(av,msize,queue,S2,DSMgd,DSMhes,n_dir,&
              &           remove_dir,0.0E0_realk,trustr,stepminustrustr,coefcor)
         write(av%lupri,*)' DSMstep Newton step',stepminustrustr+trustr
         WRITE(av%LUPRI,*) 'coefcor'
         call output(coefcor,1,msize-1,1,1,msize-1,1,1,av%lupri)
      endif
      !** Construct system of linear equations A [delta c] = b
      !print *,'number of bad directions',n_dir
      if (n_dir == 0 .and. mineival < 0.0E0_realk) then   !lowsing < 0.0E0_realk) then  
         ! Newton step not taken, restricted step should be found
         ! with a damping at least the size of the negative eigenvalue.
         x1mult = mineival - 0.1E0_realk   !lowest singularity
         call SOLVE_LINEQ(av,msize,queue,S2,DSMgd,DSMhes,n_dir,&
              &           remove_dir,x1mult,trustr,fx1,coefcor)
         if (fx1 < 0.0E0_realk) STEPFINE = .true.
      else
         !     Newton step - no damping of Hessian
         call SOLVE_LINEQ(av,msize,queue,S2,DSMgd,DSMhes,n_dir,&
              &           remove_dir,0.0E0_realk,trustr,stepminustrustr,coefcor)
         Newtonstep(naccstep+1) = stepminustrustr+trustr
         if (stepminustrustr < 0.0E0_realk) then
            !Newton step within trust radius
            if (av%info_dsm_step) then
               write(av%lupri,*)' DSMstep Newton step',stepminustrustr+trustr ,&
                    &'less than trustr :',trustr
               WRITE(av%LUPRI,*) 'coefcor'
               call output(coefcor,1,msize-1,1,1,msize-1,1,1,av%lupri)
            endif
            STEPFINE = .true.
            pure_Newton = .true.
         else
            !Newton step too long, step should be found within trustr
            !proper damping of Hessian should be found 
            if (av%info_dsm_step) &
                 &WRITE(av%LUPRI,*) 'Newton step larger than trustr, thus not taken'
         endif
         x1mult = 0.0E0_realk
         fx1 = stepminustrustr
      endif
      !** If Newton step too large or negative hessian
      ! step on border of ** trustradius is found
      if (.not. STEPFINE) then
         !Find bracketing x-values, x1mult is the rigthmost bracket 
         ds = 5.0E0_realk
         x2mult = x1mult - ds 
         call SOLVE_LINEQ(av,msize,queue,S2,DSMgd,DSMhes,n_dir,&
              &           remove_dir,x2mult,trustr,fx2,coefcor) 
         if (fx2 >= fx1) then  
            STOP 'something wrong, curve not descending'
         endif
         !Search for negative bracketing values
         do
            if (fx2 < 0.0E0_realk) then
               !found left bracket
               if (av%info_dsm_step_bracket) then
                  WRITE(av%LUPRI,*) 'bracketing x-values:',x2mult,x1mult,'fx1,fx2',fx2,fx1
               endif
               exit
            endif
            !updating right bracket
            x1mult = x2mult
            ds = ds * 2.0E0_realk
            x2mult = x2mult - ds
            call SOLVE_LINEQ(av,msize,queue, S2,DSMgd,DSMhes,n_dir,&
                 &           remove_dir,x2mult,trustr,fx2,coefcor)
         enddo
         !** find the step where the norm equals the trustr. The correct damping is found
         !   between x2mult and x1mult by bisec-search
         call STEPBISEC(av,msize,queue,S2,DSMgd,DSMhes,n_dir,remove_dir,&
              &         x1mult,x2mult,trustr,xaccmult,coefcor)
      end if

   end subroutine DSMstep

!> \brief Returns the solution vector delta c (coefcor) for a certain damping of the hessian (xlambda).
   subroutine SOLVE_LINEQ(av,msize,queue,S2,DSMgd,DSMhes,n_dir,remove_dir,&
        &                 xlambda,trustr,stepminustrust,coefcor)
      implicit none
      type(avitem),intent(in) :: av
      integer, intent(in)      ::  msize,n_dir
      TYPE(util_HistoryStore)  :: queue
      real(realk), intent(in)  :: S2(msize-1,msize-1),DSMgd(msize-1), DSMhes(msize-1,msize-1)
      real(realk), intent(in)  :: xlambda,trustr,remove_dir(:,:) 
      real(realk), intent(out) :: stepminustrust 
      real(realk), intent(out) :: coefcor(msize-1)
      !
      real(realk) :: Amat(msize-1,msize-1),bvec(msize-1),xstep,&
           &         S2x(msize-1),S2xS2xT(msize-1,msize-1)
      integer :: error, i,j,i1,j1,piv(msize-1)
      !
      !** Construct system of linear equations;
      !   xlambda*S[2] = damping of Hessian
      if (ABS(xlambda) > 1.0E-13_realk) then
         Amat = DSMhes - xlambda*S2
      else
         Amat = DSMhes
      endif
      bvec = -1.0E0_realk*DSMgd
      if (n_dir /= 0) then
         !remove bad directions i: E[2] - alpha*(S[2]x_i)(S[2]x_i)^T
         !alpha = some large number = 1.0E6_realk*E[1]_diagbase_i;
         !x_i = spectral decomposition of the bad direction
         do i=1,n_dir
            if (av%info_dsm_proj) WRITE(av%LUPRI,*)&
                 & 'removing direction, 1.0E6_realk*',ABS(remove_dir(i,msize))
            call DGEMM('N','N',msize-1,1,msize-1,1.0E0_realk,&
               & S2,msize-1,&
               & remove_dir(i,1:msize-1),msize-1,0.0E0_realk,S2x,msize-1)
            call DGEMM('N','T',msize-1,msize-1,1,1.0E0_realk,&
               &S2x,msize-1,&
               &S2x,msize-1,0.0E0_realk,S2xS2xT,msize-1)
            Amat = Amat - 1.0E6_realk*ABS(remove_dir(i,msize))*S2xS2xT
         enddo
      endif
      if (av%info_dsm_eq) then
         write(av%lupri,*)'bvec'
         call output(bvec,1,msize-1,1,1,msize-1,1,1,av%lupri)
         write(av%lupri,*)'Amat '
         call output(Amat,1,msize-1,1,msize-1,msize-1,msize-1,1,av%lupri)
      endif
      !** Solve system of linear equations - solution vector put in bvec
      CALL DGESV(msize-1,1,Amat,msize-1,piv,bvec,msize-1,error)  !LAPACK
      if (error /= 0) then 
         write(av%lupri,*) &
              &'  stepmul: Error', error,&
              &' in gaussian elimination , xlamda:', xlambda
         call output(Amat,1,msize-1,1,msize-1,msize-1,msize-1,1,av%lupri)
         stop 'stepmul: error in gaussian_solve'
      end if
      if (av%info_dsm_eq) then
         write(av%lupri,*)'Solution vector'
         call output(bvec,1,msize-1,1,1,msize-1,1,1,av%lupri)
      endif
      !** Evaluating the stepnorm 
      xstep = ABS(dot_product(bvec,MATMUL(S2,bvec)))
      xstep = SQRT(xstep)
      if (av%info_dsm_step) write(av%lupri,*)'stepnorm:',xstep,'damping:',xlambda
      !** If stepminustrust < 0 then the step is small enough, if > 0 then too big
      stepminustrust = xstep - trustr
      coefcor = bvec

   end subroutine SOLVE_LINEQ

!> \brief By bisection, the proper damping of the hessian is found
   subroutine STEPBISEC(av,msize,queue,S2,DSMgd,DSMhes,n_dir,&
        &               remove_dir,x1mult,x2mult,trustr,xaccmult,coefcor)
      implicit none
      type(avitem),intent(in) :: av
      integer, intent(in) ::  msize,n_dir
      TYPE(util_HistoryStore)  :: queue
      real(realk), intent(in)  :: S2(msize-1,msize-1),DSMgd(msize-1), DSMhes(msize-1,msize-1)
      real(realk), intent(in)  :: x1mult,x2mult,trustr,remove_dir(:,:) 
      real(realk), intent(out) :: xaccmult 
      real(realk), intent(out) :: coefcor(:)
      !
      real(realk) :: fx1, fmid, rtbis, xmid, dx
      integer :: j 
      integer, parameter :: jmax = 40

      call SOLVE_LINEQ(av,msize,queue,S2,DSMgd,DSMhes,n_dir,remove_dir,x1mult,trustr,fx1,coefcor)        
      call SOLVE_LINEQ(av,msize,queue,S2,DSMgd,DSMhes,n_dir,remove_dir,x2mult,trustr,fmid,coefcor)
      if (fx1*fmid.gt. 0.0E0_realk) then
         write(av%lupri,*)' Root must be bracketed for bisection'
         write(av%lupri,*)'x1mult,x2mult,fx1,fmid',x1mult,x2mult,fx1,fmid
         stop 'stepbisec: Root must be bracketed for bisection'
      end if
      if (fx1 .lt. 0.0E0_realk) then
         rtbis = x1mult
         dx    = x2mult - x1mult
      else
         rtbis = x2mult
         dx    = x1mult - x2mult
      end if
      do j = 1, jmax
         dx = dx * 0.5E0_realk
         xmid = rtbis + dx
         call SOLVE_LINEQ(av,msize,queue,S2,DSMgd,DSMhes,n_dir,&
              &           remove_dir,xmid,trustr,fmid,coefcor)
         if (fmid.le. 0.0E0_realk) then
            rtbis = xmid
         end if
         if (ABS(fmid) .lt. trustr*0.1) then 
            xaccmult = xmid
            if (av%info_dsm_step) write(av%lupri,*)' Convergence to trust radius: xmid,fmid,trustr', & 
                 & xmid,fmid,trustr 
            return
         end if
      end do
      write(av%lupri,*)' max number of iterations for bisection exceded'
      stop ' stepbisec: max iterations exceeded'

   end subroutine STEPBISEC

!> Analyze the gradient and hessian.
   subroutine ANALYZE_G_H(av,nqdsm,msize,queue, S2,naccstep,DSMgrad,DSMhes,&
        &                 outproj,mineival,&
        &                 minabs_eival,n_dir,remove_dir)
      implicit none
      type(avitem),intent(in) :: av
      TYPE(util_HistoryStore)  :: queue
      integer, intent(in)      :: nqdsm,msize,naccstep
      real(realk), intent(in)  :: S2(msize-1,msize-1), DSMgrad(msize-1), DSMhes(msize-1,msize-1)
      logical, intent(in)      :: outproj
      real(realk), intent(inout) :: mineival(:)
      real(realk), intent(out) :: minabs_eival,remove_dir(5,msize)
      integer, intent(out)     :: n_dir

      real(realk) :: eivec(msize-1,msize-1), Id(msize-1,msize-1),eival(msize-1),&
           & minS2eival,DSMgrad_diagbase(msize-1)
      integer :: i, ierr

      !** Initialization
      n_dir = 0  !number of bad directions 
      if (av%debug_dsm_metric) then
         !** diagonalize the S metric matrix - is it positive definite?
         eivec = S2
         Id = 0.0E0_realk
         do i = 1,msize-1
            Id(i,i) = 1.0E0_realk
         enddo
         call my_DSYGV(msize-1,eivec,Id,eival,"ANALYZE_G_H:1       ")
         !find the lowest eigenvalue
         minS2eival = eival(1)
         do i = 2,msize-1
            if (eival(i) < minS2eival) then
               minS2eival = eival(i)
            endif
         enddo
         if (minS2eival < 0) then
            WRITE(av%LUPRI,*) 'WARNING, metric matrix in DSM not positive definite'
            WRITE(av%LUPRI,*) 'mineival:',minS2eival
         endif
         !print eigenvalues
         write(av%lupri,*) 'S2 eigenvalues'
         write(av%lupri,*) (eival(i),i=1,msize-1)
      endif
      !** diagonalize the hessian and find eigenvalues
      eivec = DSMhes
      if (av%info_dsm_metric .or. av%debug_dsm_metric) then
         WRITE(av%LUPRI,*) 'S[2]'
         call OUTPUT(S2,1,msize-1,1,msize-1,msize-1,msize-1,1,av%lupri)
         WRITE(av%LUPRI,*) 'E[2]'
         call OUTPUT(eivec,1,msize-1,1,msize-1,msize-1,msize-1,1,av%lupri)
      endif
      !S2 is the metric
      Id = S2
      call my_DSYGV(msize-1,eivec,Id,eival,"ANALYZE_G_H         ")
      !find the lowest eigenvalue
      minabs_eival = ABS(eival(1))
      mineival(naccstep+1) = eival(1)
      do i = 2,msize-1
         if (eival(i) < mineival(naccstep+1)) then
            mineival(naccstep+1) = eival(i)
         endif
         if (ABS(eival(i)) < minabs_eival) then
            minabs_eival = ABS(eival(i))
         endif
      enddo
      if (av%info_dsm_eigenval) then
         !print eigenvalues
         write(av%lupri,*) 'hessian eigenvalues'
         write(av%lupri,*) (eival(i),i=1,msize-1)
      endif
      if (outproj) then
         !** get E[1] in the diagonal base too, that is accomplised by -eivec^T*E[1]
         call DGEMM('T','N',msize-1,1,msize-1,-1.0E0_realk,eivec,msize-1,DSMgrad,msize-1,0.0E0_realk,DSMgrad_diagbase,msize-1)
         do i = 1,msize-1
            if (eival(i) < dsm_hes_zero) then
               !this direction should not be considered
               n_dir = n_dir + 1
               if (n_dir > SIZE(remove_dir,1)) then
                  WRITE(av%LUPRI,*) 'too many bad directions:',SIZE(remove_dir,1)
                  STOP 'too many bad directions'
               endif
               remove_dir(n_dir,1:msize-1) = eivec(:,n_dir)
               remove_dir(n_dir,msize) = DSMgrad_diagbase(i)
               !det i'te element skal ikke medregnes i gradient norm i STEP_EVAL
            else
               !eigenvalues in increasing order, thus no more negative values
               exit
            endif
         enddo
      endif
      if (av%info_dsm_grad) then
         if (.not. outproj) then
           !** get E[1] in the diagonal base too, that is accomplised by -eivec^T*E[1]
           call DGEMM('T','N',msize-1,1,msize-1,-1.0E0_realk,eivec,msize-1,DSMgrad,msize-1,0.0E0_realk,DSMgrad_diagbase,msize-1)
         endif
         write(av%lupri,*) 'grad_diag_norm',sqrt(dot_product(DSMgrad_diagbase,DSMgrad_diagbase))  
         WRITE(av%LUPRI,'("   HESEIVAL    E[1]_DIAGBASE     E[1]")')    
         do i = 1,msize-1
            WRITE(av%LUPRI,'(3f15.8)') eival(i),DSMgrad_diagbase(i),DSMgrad(i)
         enddo
      endif

   end subroutine ANALYZE_G_H

!> Evaluate DSM step.
   subroutine STEP_EVAL(av,nqdsm,msize,his_start,n_dir,remove_dir,minstart,queue,DSMgradnorm_0,&
        &coefcor,S2,DSMgd,DSMhes,S,H1,&
        &posdef,pure_Newton,Newtonstep,coef,DSMener,trustr,deltac_tot,Fav,Dav,Delta,&
        &naccstep,steprej,end_dsm)
      use scf_stats, only: stat_current_iteration, stat_tab
      implicit none
      type(avItem),intent(inout) :: av
      integer, intent(in) :: nqdsm, msize, his_start,n_dir,minstart
      TYPE(util_HistoryStore)    :: queue
      real(realk), intent(in)    :: coefcor(msize-1),S2(msize-1,msize-1),&
                                  & DSMgd(msize-1),DSMhes(msize-1,msize-1),&
                                  & Newtonstep(:),remove_dir(:,:)
      type(Matrix), intent(in)   :: S,H1
      logical, intent(in)        :: posdef,pure_Newton
      real(realk), intent(inout) :: coef(msize-1), DSMener,trustr,deltac_tot(msize-1),DSMgradnorm_0  
      type(Matrix),intent(inout) :: Fav,Dav,Delta  !output
      integer, intent(inout)     :: naccstep, steprej
      logical, intent(out)       :: end_dsm
      real(realk) :: stpnrm, coef_old(msize-1),DSMenerold,totstp,DSMgradnorm,weights(msize),&
                   &S2x(msize-1),DSMgd_outproj(msize-1)
      type(Matrix) :: Davold, Favold , Dchange,Dtilde
      logical :: step_rej
      integer :: i,j

      !** Initializing
      end_dsm = .false.
      call mat_init(Davold,Dav%nrow,Dav%ncol)
      call mat_init(Favold,Fav%nrow,Fav%ncol)
      !** updating coefs 
      coef_old = coef
      coef = coef + coefcor
      !** Construct average density and fock matrix
      call mat_assign(Davold,Dav)
      call mat_assign(Favold,Fav)
      j = 0
      do i = 1,msize
        if (i == minstart) then
          weights(i) = 1.0E0_realk - SUM(coef)
        else
          j = j + 1
          weights(i) = coef(j)
        endif
      enddo
      call get_AVERAGE_arr(av,'D',queue,av%dsm_history_size,weights,Dav)
      call get_AVERAGE_arr(av,'F',queue,av%dsm_history_size,weights,Fav)
      call eval_delta(Dav,S,Delta)
      if (av%cfg_dsm_app == av%cfg_dsm_xtra_term) then
        call xterm_build_Gpara(av,msize,queue,Delta,S,H1)
      endif

      !** updating energy
      DSMenerold = DSMener
      call DSMenergy(av,msize,his_start,queue,S,H1,Fav,Dav,Delta,DSMener)
      if (av%info_dsm_energy) then
         WRITE(av%LUPRI,*) 'DSMener',DSMener,'nqdsm',nqdsm
      endif
      !** Evaluating the stepnorm
      stpnrm = ABS(dot_product(coefcor,MATMUL(S2,coefcor)))
      stpnrm = SQRT(stpnrm)    
      if (av%info_dsm_step) WRITE(av%LUPRI,*) 'stepnorm:',stpnrm
      !** Converged?
      !gradnorm:
      if (n_dir > 0) then
        !project out the direction in the gradient before evaluating
        !(1 - (S[2]x_i)(x_i)^T) DSMgd
        DSMgd_outproj = DSMgd
        DSMgradnorm = sqrt(dot_product(DSMgd,MATMUL(S2,DSMgd)))
        WRITE(av%LUPRI,*) 'DSMgradnorm without projections',DSMgradnorm
        do i = 1,n_dir
          call DGEMM('N','N',msize-1,1,msize-1,1.0E0_realk,&
               & S2,msize-1,remove_dir(i,1:msize-1),msize-1,0.0E0_realk,S2x,msize-1)
          DSMgd_outproj = DSMgd_outproj - dot_product(remove_dir(i,1:msize-1),DSMgd)*S2x
        enddo
        DSMgradnorm = sqrt(dot_product(DSMgd_outproj,MATMUL(S2,DSMgd_outproj)))
      else
        DSMgradnorm = sqrt(dot_product(DSMgd,MATMUL(S2,DSMgd)))
      endif
      if (nqdsm == 1) then
        DSMgradnorm_0 = DSMgradnorm
      endif
      if (av%info_dsm_step) WRITE(av%LUPRI,*) 'DSMgradnorm,DSMgradnorm_0',DSMgradnorm,DSMgradnorm_0
      if ((posdef .or. n_dir > 0) .and. DSMgradnorm < DSMgradnorm_0*1E-3_realk) then
         naccstep = naccstep + 1
         if (av%info_dsm_step) write(av%lupri,*) 'number of small steps taken:',naccstep
         deltac_tot = deltac_tot + coefcor
         totstp = ABS(dot_product(deltac_tot,MATMUL(S2,deltac_tot)))
         totstp = SQRT(totstp)
         if (av%info_dsm_step_total) then
            WRITE(av%LUPRI,*) 'deltac_tot',totstp
            call OUTPUT(deltac_tot,1,msize-1,1,1,msize-1,1,1,av%lupri)
         endif 
         if (pure_Newton) then
           stat_tab(stat_current_iteration+1,6) = 0
         else
           stat_tab(stat_current_iteration+1,6) = 9
         endif
         if (av%DEBUG_DSM_DCHANGE) then
           call mat_init(Dchange,S%nrow,S%ncol)
           call mat_init(Dtilde,S%nrow,S%ncol)
           !||Dbar - Do||
           call mat_add(1E0_realk,Dav,-1E0_realk,queue%D(his_start),Dchange)
           write(av%lupri,*) '||Dbar - Do|| = ',util_Snorm(Dchange,S),'      grepdchange'
           call get_Dtilde(S,Dav,Dtilde)
           call mat_add(1E0_realk,Dtilde,-1E0_realk,Dav,Dchange)
           write(av%lupri,*) '||Dtilde - Dbar|| = ',util_Snorm(Dchange,S),'          grepdchange'
           call mat_add(1E0_realk,Dtilde,-1E0_realk,queue%D(his_start),Dchange)
           write(av%lupri,*) '||Dtilde - Do|| = ',util_Snorm(Dchange,S),'             grepdchange'
           call mat_free(Dchange)
           call mat_free(Dtilde)
         endif
         call mat_free(Davold)
         call mat_free(Favold)
         end_dsm = .true.
         return
      endif
!    if (stat_current_iteration+1 == 10) then
!        !evaluate dE_DSM and dE_SCF and print
!        totstp = ABS(dot_product(coefcor,MATMUL(S2,coefcor)))
!        totstp = SQRT(totstp)
!        call COMPARE_EHF_EDSM(S,totstp,queue,his_start,DSMener,Dav)
!        !make trustr smaller
!        trustr = trustr * 0.5E0_realk
!        if (trustr < 1E-3_realk) STOP 'dadaa'
!        !reject step
!        coef = coef_old
!        DSMener = DSMenerold
!        Dav = Davold
!        Fav = Favold
!    else
      call update_trustr(av,msize,coefcor,DSMener,DSMenerold,DSMgd,DSMhes,trustr,step_rej)
      if (step_rej) then
         !** STEP REJECTED  Too bold a step
         steprej = steprej + 1
         coef = coef_old
         DSMener = DSMenerold
         call mat_assign(Dav,Davold)
         call mat_assign(Fav,Favold)
         call eval_delta(Dav,S,Delta)
         if (av%cfg_dsm_app == av%cfg_dsm_xtra_term) then
           call xterm_build_Gpara(av,msize,queue,Delta,S,H1)
         endif
      else
         !** STEP TAKEN
         steprej = 0
         naccstep = naccstep + 1
         if (av%info_dsm_step) write(av%lupri,*) 'number of small steps taken:',naccstep
         deltac_tot = deltac_tot + coefcor
         totstp = ABS(dot_product(deltac_tot,MATMUL(S2,deltac_tot)))
         totstp = SQRT(totstp)
         if (av%info_dsm_step_total) then
            WRITE(av%LUPRI,*) 'deltac_tot',totstp
            call OUTPUT(deltac_tot,1,msize-1,1,1,msize-1,1,1,av%lupri)
         endif
         !** gone too far?
         if( totstp > changelimit) then 
            if (posdef .or. n_dir > 0) then
              if (Newtonstep(naccstep)>Newtonstep(naccstep-1)) then
                !Newton stalls
                stat_tab(stat_current_iteration+1,6) = 2
              else
                stat_tab(stat_current_iteration+1,6) = 1
              endif
              end_dsm = .true.
            else
              stat_tab(stat_current_iteration+1,6) = 3
            endif
         endif
         !if (DEBUG_DSM_DCHANGE) then
         !  call mat_init(Dchange,S%nrow,S%ncol)
         !  call mat_init(Dtilde,S%nrow,S%ncol)
         !  !||Dbar - Do||
         !  call mat_add(1E0_realk,Dav,-1E0_realk,queue%D(his_start),Dchange)
         !  write(lupri,*) '||Dbar - Do|| = ',util_Snorm(Dchange,S),'      grepdchange'
         !  call get_Dtilde(S,Dav,Dtilde)
         !  call mat_add(1E0_realk,Dtilde,-1E0_realk,Dav,Dchange)
         !  write(lupri,*) '||Dtilde - Dbar|| = ',util_Snorm(Dchange,S),'          grepdchange'
         !  call mat_add(1E0_realk,Dtilde,-1E0_realk,queue%D(his_start),Dchange)
         !  write(lupri,*) '||Dtilde - Do|| = ',util_Snorm(Dchange,S),'             grepdchange'
         !  call mat_free(Dchange)
         !  call mat_free(Dtilde)
         !endif
         if (av%debug_dsm_Ecomp_fig) then
            call lsquit('debug_dsm_Ecomp_fig disabled! /Stinne',av%lupri)
            !if (nqdsm == 1) WRITE(LUPRI,'("# iteration ",i3,"   dsmit,dEDSM,dEHF")') &
            !                                &stat_current_iteration
            !call COMPARE_EHF_EDSM(nqdsm,S,totstp,queue,his_start,DSMener,Dav)
         endif
      endif
!   endif
      call mat_free(Davold)
      call mat_free(Favold)

   end subroutine STEP_EVAL

!> \brief Should the trust-radius be changed and should the step be rejected?
   subroutine update_trustr(av,msize,coefcor,DSMener,DSMenerold,DSMgd,DSMhes,trustr,steprej)
     implicit none
      type(avitem),intent(in) :: av
      integer, intent(in) :: msize
      real(realk), intent(in)    :: coefcor(msize-1),DSMgd(msize-1),DSMhes(msize-1,msize-1)
      real(realk), intent(in) :: DSMener,DSMenerold
      real(realk), intent(inout) :: trustr
      logical, intent(out) :: steprej
      real(realk)  ::  ratio,DSMenerpred

      steprej = .false.
      !** make a prediction of the energy after the new step, based on 
      !   the assumption that we are in the local region
      call DSMenergy_pred(msize,coefcor,DSMenerold,DSMgd,DSMhes,DSMenerpred)
      !** deciding whether the trust radius should be changed and this step rejected
      !compute "ratio", if ratio is close to 1, we are in the local region and our
      !approximation is a good one => large steps could be taken.
      !if ration is far from 1, we should take small careful steps.      
      if (ABS(DSMenerpred - DSMenerold) > 1.0E-20_realk) then
         ratio = (DSMener - DSMenerold) / (DSMenerpred - DSMenerold)
      else
         ratio = -1.0E0_realk
      endif
      !   WRITE(LUPRI,"('Ratio = E-Eo/Epred-Eo',f19.11,' = ',f19.11,'/',f19.11)") &
      !        &ratio,DSMener-DSMenerold,DSMenerpred-DSMenerold
      if ((ratio > 0.0E0_realk .and. ratio < 0.40E0_realk) .or. 1.60E0_realk < ratio) then
         !** STEP REJECTED  Too bold a step
         trustr = trustr * 0.5E0_realk
         if (av%info_dsm_ratio) then
            write(av%lupri,*) 'ratio = ',ratio,' => far from local region => smaller steps should be taken'
            write(av%lupri,*) 'predicted energy, actual energy, old energy',DSMenerpred,DSMener,DSMenerold
            write(av%lupri,*) 'predicted energy change:',DSMenerpred - DSMenerold
            write(av%lupri,*) 'actual energy change:',DSMener - DSMenerold
         endif
         !** reject step  
         if (av%info_dsm_step) write(av%lupri,*) 'step rejected'
         steprej = .true.
      else         
         !** STEP TAKEN
         if (ratio > 0.75E0_realk .and. ratio < 1.25E0_realk) then
            trustr = trustr * 1.3E0_realk
            if (av%info_dsm_trustr) write(av%lupri,*) 'trustr enlarged',trustr
         endif
         if (av%info_dsm_ratio) then
            write(av%lupri,*) 'ratio = ',ratio
            write(av%lupri,*) 'predicted energy, actual energy, old energy',DSMenerpred,DSMener,DSMenerold
            write(av%lupri,*) 'predicted energy change:',DSMenerpred - DSMenerold
            write(av%lupri,*) 'actual energy change:',DSMener - DSMenerold
         endif
      endif
   end subroutine update_trustr

!> \brief DSM is ended or restartet depending on the conditions. In particular, neg. hes. eigenvalues are not accepted as ending condition.
   subroutine EXIT_OR_RESTART(av,msize,queue,S,maxstep, &
        &his_start,coefcor,DSMener,DSMener_start,&
        &end_dsm,steprej,posdef,deltac_tot,&
        &naccstep,outproj,nqdsm,&
        &coef,Fav,Dav,mineival,trustr,Newtonstep,restart)
      USE scf_stats, only: stat_current_iteration, stat_tab
      implicit none
      type(avItem),intent(inout) :: av
      TYPE(util_HistoryStore) :: queue
      type(matrix), intent(in) :: S
      integer, intent(in) :: msize, maxstep, his_start
      real(realk), intent(in) :: coefcor(msize-1),DSMener_start
      real(realk), intent(inout) :: DSMener
      logical, intent(inout) :: end_dsm, outproj,posdef
      logical, intent(out) :: restart
      integer, intent(inout) :: steprej, naccstep,nqdsm
      real(realk), intent(inout) :: coef(msize-1),deltac_tot(msize-1)
      real(realk), intent(inout) :: mineival(maxstep),trustr,Newtonstep(maxstep)
      type(Matrix)  :: Fav, Dav
      integer :: qreduced

      restart = .false.
!
! end_dsm is set to true in STEP_EVAL if the iterations have converged
! or the max-size of the step is reached.
! The hessian is positive def. or the neg. eigenvals are projected out
!
      if (.not. end_dsm) then
!
! Have we reached the maximum number of steps or the maximum number of rejected steps?
!
         if (steprej > 5) stat_tab(stat_current_iteration+1,6) = 5E0_realk
         if (naccstep == maxstep) then
           stat_tab(stat_current_iteration+1,6) = 6E0_realk 
         endif
!
! If the Hessian is pos.def. or we have projected out the negative directions 
! and a code is found in stat_tab(k,6), end_dsm=true
!  
         if (stat_tab(stat_current_iteration+1,6) > 1E-3_realk) then
            if (naccstep.EQ.0) then
               end_dsm = .true.
            else
               if (mineival(naccstep) > dsm_hes_zero .or. outproj) then
                  end_dsm = .true.
               endif
            endif
         endif 
      endif
!
      if (.not. end_dsm .and. .not. outproj) then 
!
! Do we have a problem with the Hessian ??
!
         if (naccstep > 1 .and. mineival(1) < dsm_hes_zero)then
            IF(mineival(naccstep) < dsm_hes_zero) then
               if ((mineival(naccstep) - dsm_hes_zero) < (mineival(1) - dsm_hes_zero)*1.2E0_realk) then
                  ! lowest Hessian eigenval getting more negative
                  restart = .true.
               endif
               if (naccstep > 3 .and. ABS(mineival(naccstep)/mineival(1)) > 0.7E0_realk) then
                  ! lowest Hessian eigenval is negative and stalling
                  restart = .true.
               endif
            endif
         endif
         select case(NINT(stat_tab(stat_current_iteration+1,6)))
         case (3) 
           ! step size limit reached, but neg. hes. eigenval
           ! enlarge limit and hope to see it get pos.def.
           changelimit = changelimit*1.5E0_realk
           stat_tab(stat_current_iteration+1,6) = 0.0E0_realk
         case (4) 
           ! lowest Hessian eigenval goes from positive to negative
           restart = .true.
         case (5)
           ! Rejected more than 5 steps in a row, but hes. is neg.
           restart = .true.
         case (6)
           ! The maximum number of steps have been taken, but hes. is neg.
           restart = .true.
         end select
      endif
!
      if (restart) then
!          
! If more than 2 vectors on stock, remove one and run a new DSM
!                 

!STINNE AND POUL CHANGE, 18/1-2005
        if (queue%used_entries > 2) then
          if (queue%used_entries > 3) then
             qreduced = 1  !A '2' here didn't work well, so for now it's the same as before!
          else
             qreduced = 1
          endif
          if (av%info_dsm_proj) WRITE(av%LUPRI,*) 'chose to restart dsm with',qreduced, 'less vectors'
          call flush_queue(av,S,queue%current_position,queue%used_entries-qreduced,queue)
          WRITE(av%LUPRI,*) "*** Removing vector, now we're down to ",queue%used_entries
        else

!        if (queue%used_entries > 2) then
!          if (info_dsm_proj) WRITE(LUPRI,*) 'chose to restart dsm with one less vector'
!          call flush_queue(S,queue%current_position,queue%used_entries-1,queue)
!          WRITE(LUPRI,*) "*** Removing vector, now we're down to ",queue%used_entries
!        else
!END STINNE AND POUL CHANGE, 18/1-2005

!                  
! If less than or 2 vectors, [project out bad directions -not done!]
!                      
          restart = .false.
          end_dsm = .true.
          stat_tab(stat_current_iteration+1,6) = 7E0_realk
          coef = 0.0E0_realk
          DSMener = DSMener_start
          if (av%info_dsm_proj) WRITE(av%LUPRI,*) '<= 2 vectors - it is crap, but we are ending dsm - no steps taken'
!          if (info_dsm_proj) WRITE(LUPRI,*) '<= 2 vectors, projecting out bad directions'
!             outproj = .true.            
           !reinitializing stuff
 !           nqdsm = 0       !number of dsm iterations
   !          naccstep = 0; steprej = 0
    !         coef = 0.0E0_realk;
     !        trustr = 0.2E0_realk   !normally it works well with 0.2
      !       posdef = .false.
 !            deltac_tot = 0.0E0_realk
  !           Newtonstep = 0.0E0_realk
   !          stat_tab(stat_current_iteration+1,6) = 0E0_realk
    !         Fav = queue%F(his_start)
     !        Dav = queue%D(his_start)
      !       DSMener = DSMener_start
        endif
      endif
      if (end_dsm) then 
         stat_tab(stat_current_iteration+1,5) = &
                 &DSMener - stat_tab(stat_current_iteration+1,5)
         if (NINT(stat_tab(stat_current_iteration+1,6)) == 4) then !pos. to neg. def. hes.
           coef = coef - coefcor !backstep
         endif      
         if (av%info_dsm_exit) then
           select case(NINT(stat_tab(stat_current_iteration+1,6)))
           case (0)
             write(av%lupri,*)&
                    & 'DSM converged in :',nqdsm,' iter.'
             write(av%lupri,*) 'pos.def. hessian and pure Newton step'
           case (1)
             WRITE(av%LUPRI,*) 'exiting DSM, hit limit (pos.def. hessian)'
             write(av%lupri,*) 'Newton steps were getting smaller though, so it looked fine!!'
           case (2)
             WRITE(av%LUPRI,*) 'exiting DSM, hit limit (pos.def. hessian)'
             write(av%lupri,*) 'WARNING!! Newton steps were stalling!!'
           case (3)
             WRITE(av%LUPRI,*) 'exiting DSM, hit limit'
             write(av%lupri,*) 'WARNING!! Negative hessian eigenvalues!!'
           case (4)
             WRITE(av%LUPRI,*) 'Hessian goes from positive definite ',&
                            &'to neg. eigenvals - exiting'
           case (5)
             WRITE(av%LUPRI,*) 'rejected more than 5 steps in a row - thus exiting'
           case (6)
             WRITE(av%LUPRI,*) 'maximum number',naccstep,' of steps has been taken - exiting'
           case (7)
             WRITE(av%LUPRI,*) 'removed vectors - down to 2, but still problems - no steps taken'
           case (9)
             write(av%lupri,*)&
                    & 'DSM converged in :',nqdsm,' iter.'
             write(av%lupri,*) 'pos.def. hessian and damped Newton step'
           case default
             write(av%lupri,*) 'WARNING! something wrong in the program, unknown exit'
           end select
         endif
         if (av%info_dsm_step) WRITE(av%LUPRI,*) 'number of DSMsteps:',naccstep
         if (av%info_dsm_energy) WRITE(av%LUPRI,*) 'FINAL DSMener',DSMener
         if (outproj) stat_tab(stat_current_iteration+1,6) = &
                    &stat_tab(stat_current_iteration+1,6) + 20.0E0_realk
         return  !and exit the loop
      endif
   end subroutine EXIT_OR_RESTART

!> \brief Find the DSM energy Edsm = EHF(Do)+TrFoDav-TrFavDo-TrFoDo-TrFavDav+6TrFavDavSDav-4TrFavDavSDavSDav
   subroutine DSMenergy(av,msize,his_start,queue,S,H1,fav,dav,Delta,DSMener)
      !  Edsm = EHF(Do) +TrFoDav - TrFavDo - TrFoDo - TrFavDav
      !      + 6TrFavDavSDav - 4TrFavDavSDavSDav  
      !FIXME: Make this simpler using the Delta already available
      implicit none
      type(avItem),intent(inout) :: av
      integer, intent(in)      :: msize, his_start
      TYPe(util_HistoryStore)  :: queue
      type(Matrix), intent(in) :: S, H1, Fav, Dav, Delta
      real(realk), intent(out) :: DSMener
      real(realk) :: cie_di, trcifidi, trfavdav, trfavdavsdav,&
           & trfavdavsdavsdav,TrFD(msize), DDOT, TrDGD
      type(Matrix) :: DavS, DavSDav, FavDavS ! , pointer
      integer :: i, j, k

      call mat_init(DavS,S%nrow,S%nrow)
      call mat_init(DavSDav,S%nrow,S%nrow)
      call mat_init(FavDavS,S%nrow,S%nrow)
      DSMener = 0.0E0_realk
      !Edsm = EHF(Do) + TrFoDav - TrFavDo - TrFoDo
      DSMener = DSMener + queue%energy(his_start)                       &
              & + mat_dotproduct(queue%F(his_start),Dav)                &
              & - mat_dotproduct(Fav,queue%D(his_start))                &
              & - mat_dotproduct(queue%F(his_start),queue%D(his_start))
      TrFavDav = mat_dotproduct(Fav,Dav)
      !call cpu_time(timing1)
      call mat_mul(Dav,S,'n','n',1E0_realk,0E0_realk,DavS)
      call mat_mul(DavS,Dav,'n','n',1E0_realk,0E0_realk,DavSDav)
      call mat_mul(Fav,DavS,'n','n',1E0_realk,0E0_realk,FavDavS)
      !call cpu_time(timing2)
      !cfg_dsmtime = cfg_dsmtime + (timing2 - timing1)
      TrFavDavSDav = mat_dotproduct(Fav,DavSDav)
      trfavdavsdavsdav = mat_dotproduct(FavDavS,DavSDav)
      DSMener = DSMener - TrFavDav + &
           &    6E0_realk*TrFavDavSDav - 4E0_realk*TrFavDavSDavSDav
      if (av%cfg_dsm_app == av%cfg_dsm_xtra_term) then
        call xterm_add_Eterm(av,msize,queue,Delta,Dav,S,DSMener)
      endif
      call mat_free(DavS)
      call mat_free(DavSDav)
      call mat_free(FavDavS)

   end subroutine DSMenergy

!> \brief Find the predicted DSM energy. 
   subroutine DSMenergy_pred(msize,coef,DSMenerold,DSMgd,DSMhes,DSMenerpred)
      implicit none
      integer, intent(in) :: msize
      real(realk), intent(in) :: coef(msize-1), DSMenerold,DSMgd(msize-1), &
           & DSMhes(msize-1,msize-1)
      real(realk), intent(out) :: DSMenerpred
      real(realk) :: TMPC
      integer :: i,j

      DSMenerpred  = DSMenerold 
      do i = 1,msize-1
         DSMenerpred  = DSMenerpred + DSMgd(i)*coef(i)
      enddo
      do j = 1,msize-1
         TMPC = 0.5E0_realk*coef(j)
         do i = 1,msize-1
            DSMenerpred = DSMenerpred + coef(i)*DSMhes(i,j)*TMPC
         enddo
      enddo

   end subroutine DSMenergy_pred

!> \brief Find d_tilde = 3dsd - 2dsdsd, where d = d_bar
   subroutine get_Dtilde(S,Dbar,Dtilde)
      implicit none
      type(Matrix), intent(in) :: S,Dbar
      type(Matrix), intent(inout) :: Dtilde !output
      type(Matrix) :: SD, DSD, DSDSD

      call mat_init(SD,S%nrow,S%ncol)
      call mat_init(DSD,S%nrow,S%ncol)
      call mat_init(DSDSD,S%nrow,S%ncol)

      !** Find d_tilde = 3dsd - 2dsdsd  ; where d = d_bar
      !call cpu_time(timing1)
      call mat_mul(S,Dbar,'n','n',1E0_realk,0E0_realk,SD)
      call mat_mul(Dbar,SD,'n','n',1E0_realk,0E0_realk,DSD)
      call mat_mul(DSD,SD,'n','n',1E0_realk,0E0_realk,DSDSD)
      !call cpu_time(timing2)
      !cfg_dsmtime = cfg_dsmtime + (timing2 - timing1)
      call mat_add(3.0E0_realk,dsd, -2E0_realk,dsdsd, dtilde)

      call mat_free(SD)
      call mat_free(DSD)
      call mat_free(DSDSD)

   end subroutine get_Dtilde

!> \brief Get DSM gradient and Hessian, improved routine!
   subroutine DSMgrad_and_hes_better(msize,his_start,minstart,queue,Delta,&
                                     &S,H1,Dav,Fav,grad_out,hes_out)
      implicit none
      integer, intent(in) :: msize,his_start,minstart
      type(util_HistoryStore)  :: queue
      type(Matrix), intent(in) :: Delta,S, H1,Dav, Fav
      real(realk), intent(out) :: grad_out(msize-1), hes_out(msize-1,msize-1)
      real(realk)              :: grad(msize),hes(msize,msize)
      type(Matrix)             :: FavD(msize),FavDS(msize),DavSD(msize),&
                                & FDavS(msize), DavSDav,DavSDavSDav,FavDavS,wrk,wrk2
      integer :: i, j, i1, j1, Ndim

      Ndim = S%nrow
      grad = 0.0E0_realk
      hes = 2.0E0_realk*mat_dotproduct(queue%D(his_start),queue%F(his_start))
      call mat_init(wrk,ndim,ndim)
      call mat_init(wrk2,ndim,ndim)
      call mat_init(DavSDav,ndim,ndim)
      call mat_init(FavDavS,ndim,ndim)
      !call cpu_time(timing1)
      call mat_mul(Dav,S,'n','n',1E0_realk,0E0_realk,wrk)
      call mat_mul(wrk,Dav,'n','n',1E0_realk,0E0_realk,DavSDav)
      call mat_mul(Fav,wrk,'n','n',1E0_realk,0E0_realk,FavDavS)
      !call cpu_time(timing2)
      !cfg_dsmtime = cfg_dsmtime + (timing2 - timing1)
      !build the matrix products FavD_i, SD_i, FavD_iS, DavSD_i 
      do i=1,msize
         call mat_init(FDavS(i),ndim,ndim)
         !call cpu_time(timing1)
         call mat_mul(queue%F(his_p(i)),wrk,'n','n',1E0_realk,0E0_realk,FDavS(i))
         call mat_mul(S,queue%D(his_p(i)),'n','n',1E0_realk,0E0_realk,wrk2)
         !call cpu_time(timing2)
         !cfg_dsmtime = cfg_dsmtime + (timing2 - timing1)
         call mat_init(DavSD(i),ndim,ndim)
         !call cpu_time(timing1)
         call mat_mul(Dav,wrk2,'n','n',1E0_realk,0E0_realk,DavSD(i))
         !call cpu_time(timing2)
         !cfg_dsmtime = cfg_dsmtime + (timing2 - timing1)
         call mat_init(FavDS(i),ndim,ndim)
         !call cpu_time(timing1)
         call mat_mul(Fav,wrk2,'n','t',1E0_realk,0E0_realk,FavDS(i))
         !call cpu_time(timing2)
         !cfg_dsmtime = cfg_dsmtime + (timing2 - timing1)
         call mat_init(FavD(i),ndim,ndim)
         !call cpu_time(timing1)
         call mat_mul(Fav,queue%D(his_p(i)),'n','n',1E0_realk,0E0_realk,FavD(i))
         !call cpu_time(timing2)
         !cfg_dsmtime = cfg_dsmtime + (timing2 - timing1)
      enddo   

      do i = 1,msize
         call mat_trans(DavSD(i),wrk)
         grad(i) = grad(i) + mat_dotproduct(queue%D(his_p(i)),queue%F(his_start)) &
                 &         - mat_dotproduct(queue%D(his_start),queue%F(his_p(i))) &
                 &         - mat_dotproduct(Dav,queue%F(his_p(i)))                &
                 &         - mat_dotproduct(queue%D(his_p(i)),Fav)                &
                 &         - mat_dotproduct(Dav,queue%F(his_start))               &
                 &         - mat_dotproduct(queue%D(his_start),Fav)               &
                 &         + 6.0E0_realk*mat_dotproduct(Fav,DavSD(i))                   &
                 &         + 6.0E0_realk*mat_dotproduct(queue%F(his_p(i)),DavSDav)      &
                 &         + 6.0E0_realk*mat_dotproduct(Fav,wrk)                        &
                 &         - 4.0E0_realk*mat_TrAB(FavDavS,DavSD(i))                     &
                 &         - 4.0E0_realk*mat_dotproduct(FavDS(i),DavSDav)               &
                 &         - 4.0E0_realk*mat_dotproduct(FDavS(i),DavSDav)               &
                 &         - 4.0E0_realk*mat_TrAB(FavDavS,wrk)
         do j = 1,i
           hes(i,j) = hes(i,j) - mat_dotproduct(queue%D(his_p(i)),queue%F(his_p(j)))  &
                    &          - mat_dotproduct(queue%D(his_p(j)),queue%F(his_p(i)))  &
                    &          - mat_dotproduct(queue%D(his_start),queue%F(his_p(i))) &
                    &          - mat_dotproduct(queue%D(his_start),queue%F(his_p(j))) &
                    &          - mat_dotproduct(queue%D(his_p(i)),queue%F(his_start)) &
                    &          - mat_dotproduct(queue%D(his_p(j)),queue%F(his_start)) &
                    &          + 6.0E0_realk*mat_dotproduct(queue%F(his_p(i)),DavSD(j))     &
                    &          + 6.0E0_realk*mat_dotproduct(queue%F(his_p(j)),DavSD(i))     &
                    &          + 6.0E0_realk*mat_dotproduct(FavDS(j),queue%D(his_p(i)))     &
                    &          + 6.0E0_realk*mat_dotproduct(FavDS(i),queue%D(his_p(j)))     &
                    &          - 4.0E0_realk*mat_TrAB(FDavS(i),DavSD(j))                    &
                    &          - 4.0E0_realk*mat_TrAB(FDavS(j),DavSD(i))                    &
                    &          - 4.0E0_realk*mat_TrAB(FavDS(j),DavSD(i))                    &
                    &          - 4.0E0_realk*mat_TrAB(FavDS(i),DavSD(j))                                     
           call mat_trans(DavSD(j),wrk)
           hes(i,j) = hes(i,j) + 6.0E0_realk*mat_dotproduct(queue%F(his_p(i)),wrk)          &
                    &          - 4.0E0_realk*mat_TrAB(FDavS(i),wrk)                         
           call mat_trans(FDavS(i),wrk2)
           hes(i,j) = hes(i,j) - 4.0E0_realk*mat_TrAB(wrk2,wrk)                             &
                    &          - 4.0E0_realk*mat_TrAB(FavDS(i),wrk)
           call mat_trans(DavSD(i),wrk)
           hes(i,j) = hes(i,j) - 4.0E0_realk*mat_TrAB(FDavS(j),wrk)                         &
                    &          - 4.0E0_realk*mat_TrAB(FavDS(j),wrk)                         &
                    &          + 6.0E0_realk*mat_dotproduct(queue%F(his_p(j)),wrk)
           call mat_trans(FDavS(j),wrk2)
           hes(i,j) = hes(i,j) - 4.0E0_realk*mat_TrAB(wrk2,wrk)                             
           call mat_trans(FavDS(i),wrk)
           hes(i,j) = hes(i,j) - 4.0E0_realk*mat_TrAB(wrk,DavSD(j))
           call mat_trans(FavDS(j),wrk)
           hes(i,j) = hes(i,j) - 4.0E0_realk*mat_TrAB(wrk,DavSD(i)) 
         enddo
      enddo

!
! Symmetrize the Hessian
!
      do i = 1,msize-1
        do j = i+1,msize
          hes(i,j) = hes(j,i)
        enddo
      enddo
!
! Since we are in the unconstrained framework......
!
      i1 = 0
      do i = 1,msize
        if (i /= minstart) then
          i1 = i1+1
          grad_out(i1) = grad(i) - grad(minstart) 
          j1 = 0
          do j = 1,msize
            if (j /= minstart) then
              j1 = j1+1
              hes_out(i1,j1) = hes(i,j) - hes(i,minstart) - hes(minstart,j) + hes(minstart,minstart)
            endif
          enddo
        endif
      enddo
      do i = 1,msize
        call mat_free(FDavS(i))
        call mat_free(DavSD(i))
        call mat_free(FavDS(i))
        call mat_free(FavD(i))
      enddo
      call mat_free(wrk)
      call mat_free(wrk2)
      call mat_free(DavSDav)
      call mat_free(FavDavS)

   end subroutine DSMgrad_and_hes_better

!!> \brief Get DSM gradient and Hessian.
!   subroutine DSMgrad_and_hes(msize,his_start,minstart,queue,Delta,&
!        &S,H1,Dav,Fav,grad_out,hes_out)
!      implicit none
!      integer, intent(in) :: msize,his_start,minstart
!      type(util_HistoryStore)  :: queue
!      type(Matrix), intent(in) :: Delta,S, H1,Dav, Fav
!      real(realk), intent(out) :: grad_out(msize-1), hes_out(msize-1,msize-1)
!      real(realk)              :: grad(msize),hes(msize,msize),TrDoF(msize),TrFoD(msize),hesterm
!      type(Matrix) :: SD(msize),dDelta(msize),SDav,wrk1,wrk2,DavSDi,ddDelta
!      ! MatrixPointer above
!      integer :: i, j, i1,j1,Ndim
!
!      Ndim = S%nrow
!
!      call mat_init(SDav,ndim,ndim)
!      call mat_init(wrk1,ndim,ndim)
!      call mat_init(wrk2,ndim,ndim)
!      call mat_init(DavSDi,ndim,ndim)
!      call mat_init(ddDelta,ndim,ndim)
!      do i=1,msize
!         call mat_init(SD(i),ndim,ndim)
!         call mat_init(dDelta(i),ndim,ndim)
!      enddo
!      grad = 0
!      hes = 0
!      !call cpu_time(timing1)
!      call mat_mul(S,Dav,'n','n',1E0_realk,0E0_realk,SDav)
!      do i = 1,msize
!        call mat_mul(S,queue%D(his_p(i)),'n','n',1E0_realk,0E0_realk,SD(i))
!      enddo
!      !call cpu_time(timing2)
!      !cfg_dsmtime = cfg_dsmtime + (timing2 - timing1)
!!
!! Build derivatives of Delta (super matrices) 
!!    - the Hessian submatrices are utilized on the fly, the graident elements are stored
!!      
!      do i = 1,msize
!        dDelta(i) = queue%D(his_p(i))
!        !call mat_scal(-1E0_realk,dDelta(i))  !- Di
!        !call cpu_time(timing1)
!        call mat_mul(Dav,SD(i),'n','n',1E0_realk,0E0_realk,DavSDi)
!        call mat_mul(SDav,DavSDi,'t','n',1E0_realk,0E0_realk,wrk1)
!        !call cpu_time(timing2)
!        !cfg_dsmtime = cfg_dsmtime + (timing2 - timing1)
!        call mat_add(3E0_realk,DavSDi,-2E0_realk,wrk1,wrk2)
!        call util_get_symm_part(wrk2)
!        !call mat_daxpy(2E0_realk,wrk2,dDelta(i)) !2[3DavSDi - 2DavSDavSDi]^S
!        call mat_add(2E0_realk,wrk2, -1E0_realk,queue%D(his_p(i)), dDelta(i))
!        !call cpu_time(timing1)
!        call mat_mul(DavSDi,SDav,'n','n',1E0_realk,0E0_realk,wrk2)
!        !call cpu_time(timing2)
!        !cfg_dsmtime = cfg_dsmtime + (timing2 - timing1)
!        call mat_daxpy(-2E0_realk,wrk2,dDelta(i)) !- 2DavSDiSDav
!!
!! Hessian stuff
!!
!        do j = 1,i
!          !call cpu_time(timing1)
!          call mat_mul(queue%D(his_p(j)),SD(i),'n','n',1E0_realk,0E0_realk,wrk2) !wrk2 = DjSDi
!          !call cpu_time(timing2)
!          !cfg_dsmtime = cfg_dsmtime + (timing2 - timing1)
!          ddDelta = wrk2
!          !call mat_scal(6E0_realk,ddDelta)                  !2*3DjSDi
!          !call cpu_time(timing1)
!          call mat_mul(wrk2,SDav,'n','n',1E0_realk,0E0_realk,wrk1)
!          !call cpu_time(timing2)
!          !cfg_dsmtime = cfg_dsmtime + (timing2 - timing1)
!          !call mat_daxpy(-4E0_realk,wrk1,ddDelta)           !2*- 2DjSDiSDav
!          call mat_add(-4E0_realk,wrk1,6E0_realk,wrk2,ddDelta)
!          !call cpu_time(timing1)
!          call mat_mul(SDav,wrk2,'t','n',1E0_realk,0E0_realk,wrk1)
!          !call cpu_time(timing2)
!          !cfg_dsmtime = cfg_dsmtime + (timing2 - timing1)
!          call mat_daxpy(-4E0_realk,wrk1,ddDelta)           !2*- 2DavSDjSDi
!          !call cpu_time(timing1)
!          call mat_mul(SD(j),DavSDi,'t','n',1E0_realk,0E0_realk,wrk1)
!          !call cpu_time(timing2)
!          !cfg_dsmtime = cfg_dsmtime + (timing2 - timing1)
!          call mat_daxpy(-4E0_realk,wrk1,ddDelta)           !2*- 2DjSDavSDi
!          call util_get_symm_part(ddDelta)             ![]^S
!          !
!          ! d² Delta/d c_i d c_j has now been build
!          !
!          hes(i,j) = hes(i,j) + 2E0_realk*mat_dotproduct(Fav,ddDelta)  !2Tr Fav ddDelta
!          if (cfg_dsm_app == cfg_dsm_xtra_term) then
!            call xterm_build_hes(msize,queue,ddDelta,Delta,S,H1,hesterm) 
!            hes(i,j) = hes(i,j) + hesterm
!          endif
!        enddo
!      enddo
!!
!! Finish gradient and Hessian
!!
!      ! The constant terms added all elements are redundant in the unconstrained optimization - see end of routine
!      grad = grad - mat_dotproduct(Dav,queue%F(his_start)) - mat_dotproduct(queue%D(his_start),Fav) !-TrDavFo - TrDoFav
!      hes = hes + 2E0_realk*mat_dotproduct(queue%D(his_start),queue%F(his_start)) !2TrDoFo
!      do i = 1,msize
!        TrDoF(i) = mat_dotproduct(queue%D(his_start),queue%F(his_p(i)))
!        TrFoD(i) = mat_dotproduct(queue%D(his_p(i)),queue%F(his_start))        
!        grad(i) = grad(i) &
!                   & + 2E0_realk*mat_dotproduct(queue%F(his_p(i)),Delta)               & !2Tr Fi Delta
!                   & + 2E0_realk*mat_dotproduct(Fav,dDelta(i))                         & !2Tr Fav dDeltai 
!                   & +      mat_dotproduct(Dav,queue%F(his_p(i)))                 & !Tr Dav Fi
!                   & +      mat_dotproduct(queue%D(his_p(i)),Fav)                 & !Tr Di Fav
!                   & +    TrFoD(i) - TrDoF(i)                                       !Tr Fo Di - Tr Do Fi
!        do j = 1,i
!          hes(i,j) = hes(i,j) &
!                   & + 2E0_realk*mat_dotproduct(queue%F(his_p(i)),dDelta(j))           & !2Tr Fi dDeltaj
!                   & + 2E0_realk*mat_dotproduct(queue%F(his_p(j)),dDelta(i))           & !2Tr Fj dDeltai
!                   & +      mat_dotproduct(queue%D(his_p(i)),queue%F(his_p(j)))   & !TrDiFj
!                   & +      mat_dotproduct(queue%D(his_p(j)),queue%F(his_p(i)))     !TrDjFi
!        enddo
!      enddo
!      do i = 1,msize
!        do j = 1,i
!          hes(i,j) = hes(i,j) - TrDoF(i) - TrFoD(i) - TrFoD(j) - TrDoF(j) 
!        enddo
!      enddo
!      if (cfg_dsm_app == cfg_dsm_xtra_term) then
!        call xterm_build_deriv(msize,queue,Delta,dDelta,S,H1,grad,hes)
!      endif
!
!!
!! Symmetrize the Hessian
!!
!      do i = 1,msize-1
!        do j = i+1,msize
!          hes(i,j) = hes(j,i)
!        enddo
!      enddo
!!
!! Since we are in the unconstrained framework......
!!
!      i1 = 0
!      do i = 1,msize
!        if (i /= minstart) then
!          i1 = i1+1
!          grad_out(i1) = grad(i) - grad(minstart) 
!          j1 = 0
!          do j = 1,msize
!            if (j /= minstart) then
!              j1 = j1+1
!              hes_out(i1,j1) = hes(i,j) - hes(i,minstart) - hes(minstart,j) + hes(minstart,minstart)
!            endif
!          enddo
!        endif
!      enddo
!
!      call mat_free(SDav)
!      call mat_free(wrk1)
!      call mat_free(wrk2)
!      call mat_free(DavSDi)
!      call mat_free(ddDelta)
!      do i=1,msize
!         call mat_free(SD(i))
!         call mat_free(dDelta(i))
!      enddo
!
!   end subroutine DSMgrad_and_hes
!
!   subroutine DSMgrad_and_hes2(msize,his_start,minstart,queue,coef,S,Dav,Fav,grad_out,hes_out)
!      use linsca_debug
!      implicit none
!      integer, intent(in) :: msize,his_start,minstart
!      type(util_HistoryStore)  :: queue
!      real(realk), intent(in) :: coef(msize-1) 
!      type(Matrix), intent(in) :: S, Dav, Fav
!      real(realk), intent(out) :: grad_out(msize-1), hes_out(msize-1,msize-1)
!      real(realk)              :: grad_x(msize),hes_x(msize,msize),  & !The "EDIIS" terms
!                                & grad(msize),hes(msize,msize)
!      type(Matrix) :: SD(msize),SDav,X,FD,FDSD
!      ! MatrixPointer above
!      real(realk) :: ciTrFiDi, DDOT, TrFD1, TrFDSD, TrFDSDSD, TrFoDo, TrFoDx, TrFxDo
!      integer :: i, j, Ndim,k,i1,i2,i3,i4,in1,in2,i3end,i4end,itest1,nalloc,nmultps,&
!               & nmult, j1, j2, n1, n2,p 
!
!      Ndim = S%nrow
!      call mat_init(SDav,ndim,ndim)
!      call mat_init(X,ndim,ndim)
!      call mat_init(FD,ndim,ndim)
!      call mat_init(FDSD,ndim,ndim)
!      do i=1,msize
!         call mat_init(SD(i),ndim,ndim)
!      enddo
!      if (debug_dsm_derivatives) then
!        call debug_EDbar(msize,coef,his_start,queue,grad,hes)
!        WRITE(LUPRI,*) 'Finite difference derivatives'
!        write(lupri,*) 'DSM gradient of E(Dbar) part'
!        call output(grad,1,msize,1,1,msize,1,1,lupri)
!        write(lupri,*) 'DSM hessian of E(Dbar) part'
!        call output(hes,1,msize,1,msize,msize,msize,1,lupri)
!      endif
!!** Initialize
!      nalloc = 4 + msize   !number of matrices allocated
!      nmult = 0            !number of matrix multiplications
!      nmultps = 0          !number of trace of matrix products
!      ciTrFiDi = 0.0E0_realk
!      hes = 0.0E0_realk
!      hes_x = 0.0E0_realk
!      grad = 0.0E0_realk
!      grad_x = 0.0E0_realk
!
!      j = queue%current_position
!      i = dsm_pos
!      do k = 1,msize
!        !grad(i) = TrFoDi - TrFiDo
!        TrFoDx = mat_dotproduct(queue%F(his_start),queue%D(j))
!        TrFxDo = mat_dotproduct(queue%F(j),queue%D(his_start))
!        grad_x(i) = TrFoDx - TrFxDo
!        !hes(i,j) = - TrFoDi - TrFoDj - TrFiDo - TrFjDo
!        do p = 1,i
!          hes_x(i,p) = hes_x(i,p) - TrFoDx - TrFxDo
!        enddo
!        do p = i+1,msize
!          hes_x(p,i) = hes_x(p,i) - TrFoDx - TrFxDo
!        enddo
!        !call cpu_time(timing1)
!        call mat_mul(S,queue%D(j),'n','n',1E0_realk,0E0_realk,SD(i))
!        !call cpu_time(timing2)
!        !cfg_dsmtime = cfg_dsmtime + (timing2 - timing1)
!        j = j-1; if(j==0) j=cfg_settings(cfg_set_type)%max_history_size
!        i = i-1; if(i==0) i=dsm_history_size
!      enddo
!      !call cpu_time(timing1)
!      call mat_mul(S,Dav,'n','n',1E0_realk,0E0_realk,SDav)
!      !call cpu_time(timing2)
!      !cfg_dsmtime = cfg_dsmtime + (timing2 - timing1)
!      nmult = nmult + msize + 1
!      TrFoDo = mat_dotproduct(queue%F(his_start),queue%D(his_start))
!      !grad = grad + 2TrFoDo - TrFavDo - TrFoDav
!      do i = 1,msize
!        grad_x(i) = grad_x(i) + 2.0E0_realk*TrFoDo - mat_dotproduct(Fav,queue%D(his_start)) &
!                & - mat_dotproduct(queue%F(his_start),Dav)
!!        write(lupri,*) i,m2TrFoDo,' 2TrFoDo'
!!        write(lupri,*) i,- mat_dotproduct(Fav,queue%D(his_start)),' -TrFavDo'
!!        write(lupri,*) i,- mat_dotproduct(queue%F(his_start),Dav),' -TrFoDav'
!      enddo
!      !hes = hes + 2TrFoDo
!      do i = 1,msize
!        do j = 1,i
!          if (i == j) then
!            hes_x(i,j) = hes_x(i,j) + TrFoDo
!          else
!            hes_x(i,j) = hes_x(i,j) + 2.0E0_realk*TrFoDo
!          endif
!        enddo
!      enddo
!  !** chose F
!      j1 = queue%current_position
!      i1 = dsm_pos
!      n1 = 0
!      do
!  !** chose D
!        j2 = queue%current_position
!        i2 = dsm_pos
!        n2 = 0
!        do
!          !call cpu_time(timing1)
!          if (j1 == 0 .and. j2 == 0) then
!            call mat_mul(Fav,Dav,'n','n',1E0_realk,0E0_realk,FD)
!          elseif (j1 == 0) then
!            call mat_mul(Fav,queue%D(j2),'n','n',1E0_realk,0E0_realk,FD)
!          elseif (j2 == 0) then
!            call mat_mul(queue%F(j1),Dav,'n','n',1E0_realk,0E0_realk,FD)
!          else
!            call mat_mul(queue%F(j1),queue%D(j2),'n','n',1E0_realk,0E0_realk,FD)
!          endif
!          !call cpu_time(timing2)
!          !cfg_dsmtime = cfg_dsmtime + (timing2 - timing1)
!          nmult = nmult + 1
!          TrFD1 = mat_Tr(FD)
!          i3end = msize
!          if (MAX(i1,i2) == i1+i2 .and. i1 /= i2) then   !one of them is av         
!            grad_x(i1+i2) = grad_x(i1+i2) + TrFD1  !+TrFavDx + TrFxDav
!            grad(i1+i2) = grad(i1+i2) - 2.0E0_realk*TrFD1  !-2TrFavDx -2TrFxDav
!!            write(lupri,*) i1+i2,-TrFD1,' -TrFavDx - TrFxDav'
!          elseif (i1 /= 0 .and. i2 /= 0) then  !none of them is av
!            if (i1 == i2) then
!              hes_x(i1,i1) = hes_x(i1,i1) + TrFD1
!              hes(i1,i1) = hes(i1,i1) - 2.0E0_realk*TrFD1
!            else
!              hes_x(MAX(i1,i2),MIN(i1,i2)) = hes_x(MAX(i1,i2),MIN(i1,i2)) + TrFD1       !+ TrFxDy + TrFyDx
!              hes(MAX(i1,i2),MIN(i1,i2)) = hes(MAX(i1,i2),MIN(i1,i2)) - 2.0E0_realk*TrFD1 !-2TrFxDy -2TrFyDx
!            endif
!            i3end = 0 !if none of them are av, only contributions from SDav left
!          endif
!  !** chose SD1
!          do i3 = 0,i3end
!            if (i3 == 0) then
!              X = SDav
!            else
!              X = SD(i3)
!            endif
!            !call cpu_time(timing1)
!            call mat_mul(FD,X,'n','n',1E0_realk,0E0_realk,FDSD)
!            !call cpu_time(timing2)
!            !cfg_dsmtime = cfg_dsmtime + (timing2 - timing1)
!            nmult = nmult + 1
!            TrFDSD = mat_Tr(FDSD)
!            i4end = msize
!            if (MAX(i1,i2) == i1+i2 .and. i1 /= i2) then
!              if (i3 == 0) then
!                grad(i1+i2) = grad(i1+i2) + 6E0_realk*TrFDSD  !+6TrFavDxSDav +6TrFxDavSDav
!!                write(lupri,*) i1+i2,6.0E0_realk*TrFDSD,' +6TrFavDxSDav +6TrFxDavSDav'
!              else
!                !+6TrFavDySDx +6TrFavDxSDy +6TrFyDavSDx + 6TrFxDavSDy
!                hes(MAX(i1,i2,i3),MIN(MAX(i1,i2),i3)) = hes(MAX(i1,i2,i3),MIN(MAX(i1,i2),i3)) + 6E0_realk*TrFDSD
!                i4end = 0
!              endif 
!            elseif (i1 /= 0 .and. i2 /= 0) then  !i3 = 0
!              !+6TrFyDxSDav +6TrFxDySDav
!              hes(MAX(i1,i2),MIN(i1,i2)) = hes(MAX(i1,i2),MIN(i1,i2)) + 6E0_realk*TrFDSD
!              i4end = 0
!            else   !both i1 and i2 are zero
!              if (i3 /= 0) then  !6TrFavDavSDx
!                grad(i3) = grad(i3) + 6E0_realk*TrFDSD
!!                write(lupri,*) i3,6E0_realk*TrFDSD,' +6TrFavDavSDx'
!              endif
!            endif    
!  !** choose SD2
!            do i4 = 0,i4end
!              if (i4 == 0) then
!                X = SDav
!              else
!                X = SD(i4)
!              endif
!              TrFDSDSD = mat_dotproduct(FDSD,X)
!              nmultps = nmultps + 1
!              if (MAX(i1,i2) == i1+i2 .and. i1 /= i2) then
!                if (i3 == 0 .and. i4 == 0) then
!                  grad(i1+i2) = grad(i1+i2) - 4E0_realk*TrFDSDSD  !-4TrFavDxSDavSDav -4TrFxDavSDavSDav
!!                  write(lupri,*) i1+i2,-4E0_realk*TrFDSDSD,' -4TrFavDxSDavSDav -4TrFxDavSDavSDav' 
!                else
!                  !-4TrFavDySDavSDx -4TrFavDxSDavSDy -4TrFyDavSDavSDx -4TrFxDavSDavSDy
!                  !-4TrFavDySDxSDav -4TrFavDxSDySDav -4TrFyDavSDxSDav -4TrFxDavSDySDav
!                  in1 = MAX(i1,i2); in2 = MAX(i3,i4)
!                  hes(MAX(in1,in2),MIN(in1,in2)) = hes(MAX(in1,in2),MIN(in1,in2)) - 4E0_realk*TrFDSDSD
!                endif
!              elseif (i1 /= 0 .and. i2 /= 0) then !i3 = 0, i4 = 0
!                !-4TrFyDxSDavSdav -4TrFxDySDavSDav
!                hes(MAX(i1,i2),MIN(i1,i2)) = hes(MAX(i1,i2),MIN(i1,i2)) - 4E0_realk*TrFDSDSD
!              else  !both i1 and i2 are zero
!                if (MAX(i3,i4) == i3+i4 .and. i3 /= i4) then
!                  grad(i3+i4) = grad(i3+i4) - 4E0_realk*TrFDSDSD  !-4TrFavDavSDavSDx -4TrFavDavSDxSDav
!!                  write(lupri,*) i3+i4,- 4E0_realk*TrFDSDSD,' -4TrFavDavSDavSDx -4TrFavDavSDxSDav' 
!                elseif (i3 /= 0 .and. i4 /= 0) then
!                  !-4TrFavDavSDySDx -4TrFavDavSDxSDy
!                  hes(MAX(i3,i4),MIN(i3,i4)) = hes(MAX(i3,i4),MIN(i3,i4)) - 4E0_realk*TrFDSDSD
!                endif
!              endif
!            enddo
!          enddo
!          if (j2 == 0) exit
!          call update_indexes(n2,msize,i2,j2)
!        enddo
!        if (j1 == 0) exit
!        call update_indexes(n1,msize,i1,j1)
!      enddo
!      do i = 1,msize
!        do j = 1,i
!          hes(j,i) = hes(i,j) + hes(j,i)
!          hes_x(j,i) = hes_x(i,j) + hes_x(j,i)
!        enddo
!      enddo
!      
!      if (info_dsm_derivatives .or. debug_dsm_derivatives) then
!        i1 = 0
!        do i = 1,msize
!          if (i /= minstart) then
!            i1 = i1+1
!            grad_out(i1) = grad_x(i) - grad_x(minstart) 
!            j1 = 0
!            do j = 1,msize
!              if (j /= minstart) then
!                j1 = j1+1
!                hes_out(i1,j1) = hes_x(i,j) - hes_x(i,minstart) - hes_x(minstart,j) + hes_x(minstart,minstart)
!              endif
!            enddo
!          endif
!        enddo
!         !print gradient and hessian
!         write(lupri,*) 'DSM gradient E(Dbar) term'
!         call output(grad_out,1,msize-1,1,1,msize-1,1,1,lupri)
!         write(lupri,*) 'DSM hessian E(Dbar) term'
!         call output(hes_out,1,msize-1,1,msize-1,msize-1,msize-1,1,lupri)
!      !   !print gradient and hessian
!      !   write(lupri,*) 'DSM gradient'
!      !   call output(grad,1,msize,1,1,msize,1,1,lupri)
!      !   write(lupri,*) 'DSM hessian'
!      !   call output(hes,1,msize,1,msize,msize,msize,1,lupri)
!      endif
!      grad = grad + grad_x
!      hes = hes + hes_x
!      i1 = 0
!      do i = 1,msize
!        if (i /= minstart) then
!          i1 = i1+1
!          grad_out(i1) = grad(i) - grad(minstart) 
!          j1 = 0
!          do j = 1,msize
!            if (j /= minstart) then
!              j1 = j1+1
!              hes_out(i1,j1) = hes(i,j) - hes(i,minstart) - hes(minstart,j) + hes(minstart,minstart)
!            endif
!          enddo
!        endif
!      enddo
!      if (info_dsm_derivatives .or. debug_dsm_derivatives) then
!         !print gradient and hessian
!         write(lupri,*) 'DSM gradient'
!         call output(grad_out,1,msize-1,1,1,msize-1,1,1,lupri)
!         write(lupri,*) 'DSM hessian'
!         call output(hes_out,1,msize-1,1,msize-1,msize-1,msize-1,1,lupri)
!      endif
!      !WRITE(LUPRI,*)'nalloc=',nalloc,' nmult=',nmult,' nmultps=',nmultps
!      do i=1,msize
!         call mat_free(SD(i))
!      enddo
!      call mat_free(SDav)
!      call mat_free(X)
!      call mat_free(FD)
!      call mat_free(FDSD)
!
!   end subroutine DSMgrad_and_hes2

!> \brief Update DSM indices.
   subroutine update_indexes(av,n,msize,i,j)
     implicit none
     type(avItem), intent(in) :: av
     integer, intent(inout) :: n  !total passes in the loop
     integer, intent(in)    :: msize ! to number of vectors to run through
     integer, intent(inout) :: i,j !pointer to dsm and all vectors respectively

     n = n + 1
     if (n == msize) then
       j = 0  !when indexes are 0, Dav and Fav are chosen
       i = 0
     else
       j = j - 1; if (j == 0) j = av%cfg_settings(av%cfg_set_type)%max_history_size
       i = i - 1; if (i == 0) i = av%dsm_history_size
     endif

   end subroutine update_indexes
   
!> \brief DSM debug routine.
!   subroutine debug_dchange_dsm(msize,minstart,his_start,queue,S,weights,coef,Dav)
!     use scf_stats
!     use scf_utilities
!     implicit none
!     integer, intent(in) :: msize,minstart,his_start
!     type(util_HistoryStore) :: queue
!     type(matrix),intent(in) :: S
!     real(realk), intent(in) :: weights(msize), coef(msize-1)
!     type(Matrix),intent(inout) :: Dav !scratch in a way
!     type(Matrix) :: Dchange,Dtilde
!     real(realk) :: coef_norm
!     integer :: i
!
!     call mat_init(Dchange,S%nrow,S%ncol)
!     call mat_init(Dtilde,S%nrow,S%ncol)
!     if (stat_current_iteration+1 == 3) then
!       write(lupri,*) '||Dbar - Do||^2                                  grepdchangefinalgraf1'
!       write(lupri,*) '||Dtilde - Dbar||^2                                grepdchangefinalgraf2'
!       write(lupri,*) '||c||^2                                   grepdchangefinalgraf3'
!     endif
!     call get_AVERAGE_arr(av,'D',queue,dsm_history_size,weights,Dav)
!
!     !||Dbar - Do||^2
!     call mat_add(1E0_realk,Dav,-1E0_realk,queue%D(his_start),Dchange)
!     write(lupri,'(i5,f20.13,10x,"grepdchangefinalgraf1")') stat_current_iteration+1,&
!                                                          & util_Snorm(Dchange,S)
!     
!     !||Dtilde - Dbar||^2 
!     call get_Dtilde(S,Dav,Dtilde)
!     call mat_add(1E0_realk,Dtilde,-1E0_realk,Dav,Dchange)
!     write(lupri,'(i5,f20.13,10x,"grepdchangefinalgraf2")') stat_current_iteration+1,&
!                                                          & util_Snorm(Dchange,S)
!     
!     !||c||^2
!     coef_norm = 0.0E0_realk
!     do i = 1,msize-1
!       coef_norm = coef_norm + coef(i)**2
!     enddo
!     write(lupri,'(i5,f20.13,10x,"grepdchangefinalgraf3")') stat_current_iteration+1,&
!                                                                  & coef_norm
!
!     if (.true.) then
!       write(lupri,*) ' c_i  ||D_i - D_0||^2 , it. ',stat_current_iteration+1,'    grepdchangespecific'
!       do i = 1,msize
!         if (i /= minstart) then
!           call mat_add(1E0_realk,queue%D(his_p(i)),-1E0_realk,queue%D(his_start),Dchange)
!           write(lupri,'(f10.5,f20.10 ,10x,"grepdchangespecific")') weights(i), util_Snorm(Dchange,S) 
!         endif
!       enddo
!     endif
!     call mat_free(Dchange)
!     call mat_free(Dtilde)
!
!   end subroutine debug_dchange_dsm   

!>  \brief Check the densities one by one, comparing it to Do. If the purification doesn't converge, they are removed from the history
   subroutine CLEAN_DHISTORY(av,S,queue)
      implicit none
      type(avItem), intent(in) :: av
      type(Matrix), intent(in) :: S
      type(util_HistoryStore) :: queue
      type(Matrix) :: Dbar_k,X,DSD,DSDSD
      real(realk) :: minEner,precond_norm
      integer :: k,his_start,ndim
 
      !The density corresponding to the lowest SCF energy is Do
      minEner = queue%Energy(queue%current_position)
      his_start = queue%current_position 
      do k = 1,queue%used_entries
        if (queue%Energy(k) < minEner) then
          minEner = queue%Energy(k)
          his_start = k
        endif
      enddo
      ndim = S%nrow
      call mat_init(Dbar_k,ndim,ndim)
      call mat_init(X,ndim,ndim)
      call mat_init(DSD,ndim,ndim)
      call mat_init(DSDSD,ndim,ndim)
      do k = 1,queue%used_entries
        if (k/=his_start) then
          call mat_add(1E0_realk,queue%D(his_start),0.2E0_realk,queue%D(k),Dbar_k)
          !call cpu_time(timing1)
          call mat_mul(S,Dbar_k,'n','n',1E0_realk,0E0_realk,X)
          call mat_mul(Dbar_k,X,'n','n',1E0_realk,0E0_realk,DSD)
          call mat_mul(DSD,X,'n','n',1E0_realk,0E0_realk,DSDSD)
          !call cpu_time(timing2)
          !cfg_dsmtime = cfg_dsmtime + (timing2 - timing1)
          call mat_add(3.0E0_realk,DSD,-2.0E0_realk,DSDSD,X)
          call mat_daxpy(-1.0E0_realk,Dbar_k,X)   !X = Dtilde_k - Dbar_k
          call mat_daxpy(-1.0E0_realk,queue%D(his_start),Dbar_k)   !Dbar_k = Dbar_k - Do
          precond_norm = util_Snorm(X,S)/util_Snorm(Dbar_k,S)
          write(av%lupri,*) "k's precond_norm:",k,precond_norm,util_Snorm(X,S),util_Snorm(Dbar_k,S)
          if (precond_norm > 0.5E0_realk) then
            !remove from history 
          endif
        endif
      enddo
      call mat_free(Dbar_k)
      call mat_free(X)
      call mat_free(DSD)
      call mat_free(DSDSD)
   end subroutine CLEAN_DHISTORY
 
!> \brief Find Mu DC norm ?
   subroutine MU_DCNORM(av,msize,queue,S2,DSMgd,DSMhes,n_dir,remove_dir,trustr)
      use scf_stats, only: stat_current_iteration, stat_tab
      implicit none
      type(avitem),intent(in) :: av
      integer, intent(in) ::  msize,n_dir
      Type(util_HistoryStore) :: queue
      real(realk), intent(in) :: S2(msize-1,msize-1),DSMgd(msize-1), DSMhes(msize-1,msize-1)
      real(realk), intent(in) :: trustr,remove_dir(:,:)
      real(realk) :: mu, stepminustrust, coefcor(msize-1),stepmtr_uproj

      mu = -10.0E0_realk
      WRITE(av%LUPRI,*) 'iteration ',stat_current_iteration+1,'   grepdsm'
      WRITE(av%LUPRI,"('mu      stepnorm       stepnorm u. proj.   ;   trustr = ',f8.3,'  grepdsm')") trustr
      do 
         if (mu > 10.0E0_realk) exit
         call SOLVE_LINEQ(av,msize,queue,S2,DSMgd,DSMhes,n_dir,remove_dir,mu,trustr,stepminustrust,coefcor)
         call SOLVE_LINEQ(av,msize,queue,S2,DSMgd,DSMhes,0,remove_dir,mu,trustr,stepmtr_uproj,coefcor)
         WRITE(av%LUPRI,"(f10.3,2f30.20,'  grepdsm')") mu,stepminustrust+trustr,stepmtr_uproj+trustr
         mu = mu + 0.02E0_realk
      enddo

   end subroutine MU_DCNORM
      
!> \brief Find the minimum of Escf for alpha*coef
   subroutine use_first_direction(av,msize,minstart,his_start,queue,S,H1,Delta,coef)
      use scf_stats, only: stat_current_iteration, stat_tab
!      use fock_evaluator
      implicit none
      type(avItem), intent(inout) :: av
      !type(diagItem),intent(in) :: diag
      integer, intent(in) :: msize, minstart,his_start
      type(util_historyStore) :: queue
      type(Matrix), intent(in) :: S,H1
      type(matrix), intent(inout) :: Delta
      real(realk), intent(inout) :: coef(msize-1)
      type(Matrix) :: scr1,scr2
      real(realk) :: alpha,cur_coef(msize-1),weights(msize),EHFcur,dESCF, &
                   & DSMener,dESCFold
      integer :: i,j,k,ndim
      real(realk) :: da

      ndim = S%nrow
      call mat_init(scr1,ndim,ndim)
      call mat_init(scr2,ndim,ndim)

      da = 0.2E0_realk !should be smaller than 1 :)
      if (av%debug_dsm_linesearch) then
        !WRITE(LUPRI,*) 'alpha, dESCF, EHFcur    grep_dsmone',stat_current_iteration+1
        !alpha = 0.0E0_realk
        !do k = 1,30
        !  !make the scaling of the direction
        !  cur_coef = alpha*coef
        !  call get_ESCF_from_coefs(av,msize,minstart,queue,S,cur_coef,EHFcur,scr1,scr2)
        !  !** Find the energy decrease comp. to the latest it.
        !  dESCF = EHFcur - stat_tab(stat_current_iteration,1)
        !  !** write to file
        !  WRITE(LUPRI,"(f5.1,f19.15,f18.12,'   grep_dsmone')") alpha,dESCF,EHFcur
        !  !** increase alpha
        !  alpha = alpha + da
        !enddo
      endif
      alpha = 1E0_realk
      if (av%cfg_dsm_app == av%cfg_dsm_search) then
        call lsquit('cfg_dsm_search disabled! /Stinne',av%lupri)
        !!linesearch in the coef direction with alpha as variable
        !alpha = 1E0_realk
        !cur_coef = coef
        !call get_ESCF_from_coefs(av,msize,minstart,queue,S,cur_coef,EHFcur,scr1,scr2)
        !!** Find the energy decrease comp. to the latest it.
        !dESCFold = EHFcur - stat_tab(stat_current_iteration,1)
        !alpha = 1E0_realk + da
        !cur_coef = alpha*coef
        !call get_ESCF_from_coefs(av,msize,minstart,queue,S,cur_coef,EHFcur,scr1,scr2)
        !dESCF = EHFcur - stat_tab(stat_current_iteration,1)
        !write(lupri,*) 'alpha,dESCF',alpha,dESCF
        !if (dESCF < DESCFold) then
        !  !alpha should be > 1
        !  do
        !    dESCFold = dESCF
        !    alpha = alpha + da
        !    cur_coef = alpha*coef
        !    call get_ESCF_from_coefs(av,msize,minstart,queue,S,cur_coef,EHFcur,scr1,scr2)
        !    dESCF = EHFcur - stat_tab(stat_current_iteration,1)
        !    write(lupri,*) 'alpha,dESCF',alpha,dESCF
        !    if (dESCF > dESCFold) then
        !      alpha = alpha - da
        !      WRITE(LUPRI,*) 'chose alpha:',alpha,stat_current_iteration+1
        !      coef = alpha*coef
        !      exit
        !    endif
        !    if (alpha > 500) then
        !      WRITE(LUPRI,*) 'WARNING - some problem, alpha got larger than 500, just stopping it'
        !      alpha = 500
        !      WRITE(LUPRI,*) 'chose alpha:',alpha,stat_current_iteration+1
        !      exit
        !    endif
        !    da = MAX(alpha/5E0_realk,0.2E0_realk)
        !  enddo
        !else
        !  !alpha should be <= 1
        !  alpha = 1E0_realk - da
        !  cur_coef = alpha*coef
        !  call get_ESCF_from_coefs(av,msize,minstart,queue,S,cur_coef,EHFcur,scr1,scr2)
        !  dESCF = EHFcur - stat_tab(stat_current_iteration,1)
        !  write(lupri,*) 'alpha,dESCF',alpha,dESCF
        !  if (dESCF < dESCFold) then
        !    !alpha should be < 1
        !    do
        !      dESCFold = dESCF
        !      alpha = alpha - da
        !      if (alpha < 0.0E0_realk) alpha = 0.0E0_realk
        !      cur_coef = alpha*coef
        !      call get_ESCF_from_coefs(av,msize,minstart,queue,S,cur_coef,EHFcur,scr1,scr2)
        !      dESCF = EHFcur - stat_tab(stat_current_iteration,1)
        !      write(lupri,*) 'alpha,dESCF',alpha,dESCF
        !      if (dESCF > dESCFold) then
        !        alpha = alpha + da
        !        coef = alpha*coef
        !        WRITE(LUPRI,*) 'chose alpha:',alpha,stat_current_iteration+1
        !        exit
        !      endif
        !      if (alpha < 1E-7_realk) then
        !        WRITE(LUPRI,*) 'WARNING - alpha hit zero => no step is taken'
        !        coef = 0.0E0_realk
        !        WRITE(LUPRI,*) 'chose alpha: 0.0E0_realk',stat_current_iteration+1
        !        exit
        !      endif
        !      da = MAX(alpha/5E0_realk,0.2E0_realk)
        !    enddo
        !  else
        !    !alpha = 1
        !    alpha = 1E0_realk
        !    WRITE(LUPRI,*) 'chose alpha: 1E0_realk',stat_current_iteration+1
        !    !coef = coef
        !  endif
        !endif
        !stat_tab(stat_current_iteration+1,6) = 11E0_realk
      else
        stat_tab(stat_current_iteration+1,6) = 10E0_realk
      endif
      stat_tab(stat_current_iteration+1,11) = alpha
      !construct the weights connected to the new coefs
      j = 0
      do i = 1,msize
        if (i == minstart) then
          weights(i) = 1.0E0_realk - SUM(coef)
        else
          j = j + 1
          weights(i) = coef(j)
        endif
      enddo
      !get the corresponding Dbar and Fbar
      call get_AVERAGE_arr(av,'D',queue,av%dsm_history_size,weights,scr1)
      call get_AVERAGE_arr(av,'F',queue,av%dsm_history_size,weights,scr2)
      call eval_delta(scr1,S,Delta)
      if (av%cfg_dsm_app == av%cfg_dsm_xtra_term) then
        call xterm_build_Gpara(av,msize,queue,Delta,S,H1)
      endif
      !get the corresponding DSMener
      call DSMenergy(av,msize,his_start,queue,S,H1,scr2,scr1,Delta,DSMener)
      stat_tab(stat_current_iteration+1,5) = &
              &DSMener - stat_tab(stat_current_iteration+1,5)
      call mat_free(scr1)
      call mat_free(scr2)
   end subroutine use_first_direction

!> \brief DSM debug routine
!   subroutine debug_linesearch(msize,minstart,his_start,queue,S,H1,coef)
!     use dal_interface
!      use scf_stats
!      use fock_evaluator
!      use density_optimization
!!   Find the minimum of Escf for alpha*coef
!      implicit none
!      integer, intent(in) :: msize, minstart, his_start
!      type(util_historyStore) :: queue
!      type(Matrix), intent(in) :: S, H1
!      real(realk), intent(in) :: coef(msize-1)
!      type(Matrix) :: scr1,scr2,Delta, Gorth, Dorth_idem
!      real(realk) :: alpha,cur_coef(msize-1),weights(msize),EHFcur,dESCF, &
!                   & DSMener,DSMener_imp, dESCFold, alpha_opt, dEDSM, dEDSM_imp, &
!                   & coef_parr(msize), TrDparaGDpara, TrDorthGDpara, TrDorthGDorth, junk, TrDeltaGDelta
!      real(realk), allocatable :: Dfull(:), Ffull(:)
!      integer :: i,j,k,ndim
!      real(realk) :: da, dummy1, dummy2
!      logical :: min_found
!!simen
!      real(realk), allocatable :: wrk(:)
!!simen
!
!      ndim = S%nrow
!      call mat_init(scr1,ndim,ndim)
!      call mat_init(scr2,ndim,ndim)
!      call mat_init(Delta,ndim,ndim)
!      call mat_init(Gorth,ndim,ndim)
!      call mat_init(Dorth_idem,ndim,ndim)
!      allocate(Ffull(ndim*ndim),Dfull(ndim*ndim))
!      min_found = .false.
!
!      da = 0.2E0_realk !should be smaller than 1 :)
!      dESCF = 0.0E0_realk
!      WRITE(LUPRI,*) 'alpha, dESCF, dEDSM, dEDSM_imp, TrDpGp, TrDoGp, TrDeltaGDelta    grep_dsmone',stat_current_iteration+1
!      alpha = 0.0E0_realk
!      cur_coef = alpha*coef
!      call get_ESCF_from_coefs(msize,minstart,queue,S,cur_coef,EHFcur,scr1,scr2)
!      dESCF = EHFcur - stat_tab(stat_current_iteration,1)
!      dEDSM = 0.0E0_realk
!      dEDSM_imp = 0.0E0_realk
!      TrDparaGDpara = 0.0E0_realk
!      TrDorthGDpara = 0.0E0_realk
!      TrDeltaGDelta = 0.0E0_realk
!      WRITE(LUPRI,"(f5.1,6f18.10,'   grep_dsmone')") alpha,dESCF,dEDSM,dEDSM_imp,TrDparaGDpara,TrDorthGDpara ,TrDeltaGDelta
!      !** increase alpha
!      alpha = alpha + da
!      do
!        if (alpha > 20.0E0_realk) exit 
!        !make the scaling of the direction
!        cur_coef = alpha*coef
!        call get_ESCF_from_coefs(msize,minstart,queue,S,cur_coef,EHFcur,scr1,scr2)
!        !** Find the energy decrease comp. to the latest it.
!        dESCFold = dESCF
!        dESCF = EHFcur - stat_tab(stat_current_iteration,1)
!        if (dESCF > dESCFold .and. .not. min_found) then
!          alpha_opt = alpha - da
!          min_found = .true.
!        endif
!!** evaluate the DSM energy model           
!        j = 0
!        do i = 1,msize
!          if (i == minstart) then
!            weights(i) = 1.0E0_realk - SUM(cur_coef)
!          else
!            j = j + 1
!            weights(i) = cur_coef(j)
!          endif
!        enddo
!        call get_AVERAGE_arr(av,'D',queue,dsm_history_size,weights,scr1)  !scr1 = Dav
!        call get_AVERAGE_arr(av,'F',queue,dsm_history_size,weights,scr2)  !scr2 = Fav
!        call eval_delta(scr1,S,Delta)
!        call DSMenergy(msize,his_start,queue,S,H1,scr2,scr1,Delta,DSMener)
!        dEDSM = DSMener - stat_tab(stat_current_iteration,1)
!!** evaluate the missing term TrDdelta G(Ddelta)
!         if (cfg_unres) STOP 'debug_linesearch not implemented for open shell'
!         call mat_to_full(Delta,2E0_realk,Dfull)
!
!         CALL di_FCK2AO(.TRUE.,Ffull,Dfull)
!         call mat_set_from_full(Ffull,1E0_realk,Gorth)
!        !call fck_get_fock(Delta,Gorth,junk)
!        !call mat_daxpy(-1E0_realk,H1,Gorth)
!        TrDeltaGDelta = mat_dotproduct(Delta,Gorth)
!!** evaluate the improved DSM energy model with the term Tr (2 Delta - Delta_para)G(Delta_para)
!        call util_GET_PROJ_PART(queue,Delta,S,coef_parr)  !coef for Delta_para 
!        call get_average_arr('F',queue, dsm_history_size, coef_parr, scr1)  !scr1 = Gpara
!        !correction to F if sum_i c_i /= 1
!        !and withdrawel of the one-electron part of F: G=F-h
!        call mat_daxpy(-SUM(coef_parr(1:msize)),H1,scr1)
!        call get_average_arr('D',queue, dsm_history_size, coef_parr, scr2)  !scr2 = Dpara
!        TrDparaGDpara = mat_dotproduct(scr1,scr2)
!        call mat_scal(-1E0_realk,scr2)   ! - Dpara
!        call mat_daxpy(1E0_realk,Delta,scr2)   !Dorth = Delta - Dpara
!        TrDorthGDpara = mat_dotproduct(scr1,scr2)
!!        call GET_Didem(scr2,S,Gorth,Dorth_idem)  !Gorth = scratch
!!        call fck_get_fock(Dorth_idem,Gorth,junk)  
!!        call mat_daxpy(-1E0_realk,H1,Gorth)
!!        TrDorthGDorth = mat_dotproduct(scr2,Gorth)      
!        call mat_daxpy(1E0_realk,Delta,scr2)   !2 Delta - Dpara
!        DSMener_imp = DSMener +  mat_dotproduct(scr2,scr1)
!        dEDSM_imp = DSMener_imp - stat_tab(stat_current_iteration,1)
!        !** write to file
!        WRITE(LUPRI,"(f5.1,6f18.10,'   grep_dsmone')") alpha,dESCF,dEDSM,dEDSM_imp,TrDparaGDpara,TrDorthGDpara ,TrDeltaGDelta
!        !** increase alpha
!        alpha = alpha + da
!      enddo
!      WRITE(LUPRI,*) 'optimal alpha ',alpha_opt,'  grep_dsmone'
!
!      call mat_free(scr1)
!      call mat_free(scr2)
!      call mat_free(Gorth) 
!      call mat_free(Delta)
!      call mat_free(Dorth_idem)
!      deallocate(Dfull,Ffull)
!
!   end subroutine debug_linesearch

!> \brief Find the ESCF for a certain set of coefficients
!   subroutine get_ESCF_from_coefs(av,msize,minstart,queue,S,coef,ESCF,scr1,scr2)
!      use scf_stats
!      use fock_evaluator
!      use density_optimization
!      implicit none
!      type(avItem),intent(inout) :: av
!      integer, intent(in) :: msize, minstart
!      type(util_historyStore) :: queue
!      type(Matrix), intent(in) :: S
!      real(realk), intent(in) :: coef(msize-1)
!      type(Matrix), intent(inout) :: scr1,scr2 !scratchspace preallocated
!      real(realk) :: ESCF
!      type(Matrix) :: Didem,Dav,scr
!      real(realk) :: weights(msize)
!      integer :: i,j,ndim
!
!      ndim = S%nrow
!      call mat_init(Didem,ndim,ndim)
!      
!      j = 0
!      do i = 1,msize
!        if (i == minstart) then
!          weights(i) = 1.0E0_realk - SUM(coef)
!        else
!          j = j + 1
!          weights(i) = coef(j)
!        endif
!      enddo
!      !get the corresponding Dbar and make it idempotent
!      call get_AVERAGE_arr(av,'D',queue,dsm_history_size,weights,scr2)
!      call GET_Didem(scr2,S,Didem)  
!      !** Find new Fock matrix and SCF energy
!      call fck_get_fock(Didem,scr1,ESCF)
!
!      call mat_free(Didem)
!
!   end subroutine get_ESCF_from_coefs

!> \brief DSM debug routine.
!   subroutine debug_emodel(msize,his_start,weights,DSMener_start,queue,S,H1)
!     use scf_stats
!     implicit none
!     integer, intent(in) :: msize, his_start
!     real(realk), intent(in) :: weights(msize),DSMener_start
!     type(util_historyStore) :: queue
!     type(Matrix), intent(in) :: S,H1
!     type(Matrix) :: Fav,Dav,Delta
!     real(realk) :: DSMener, ESCF, ESCF_pur, ESCF_dbar
!     integer :: ndim
!
!     ndim = S%nrow
!     call mat_init(Fav,ndim,ndim)
!     call mat_init(Dav,ndim,ndim)
!     call mat_init(Delta,ndim,ndim)
!
!     if (stat_current_iteration+1 == 3) then
!       WRITE(lupri,*) 'it     EDSM(Do)       ESCF(Do)        EDSM(Dtilde)'&
!            &//'ESCF(Dtilde)      ESCF(Dtilde_idem)     ESCF(Dbar)    '&
!            &//'grep_dsm_emodel'
!     endif
!
!     call get_AVERAGE_arr('D',queue,dsm_history_size,weights,Dav)
!     call get_AVERAGE_arr('F',queue,dsm_history_size,weights,Fav)
!     call eval_delta(Dav,S,Delta)
!     call DSMenergy(msize,his_start,queue,S,H1,fav,dav,Delta,DSMener) 
!
!     call get_ESCF_dtilde(msize,queue,S,weights,ESCF,ESCF_pur,ESCF_dbar)
!     
!     WRITE(LUPRI,'(i4,6f18.10,"    grep_dsm_emodel")') &
!          &stat_current_iteration+1, DSMener_start, queue%energy(his_start),&
!          & DSMener, ESCF,ESCF_pur, ESCF_dbar
!
!     call mat_free(Fav)
!     call mat_free(Dav)
!     call mat_free(Delta)
!
!   end subroutine debug_emodel

!> \brief Find the ESCF for a certain set of coefficients.
!   subroutine get_ESCF_dtilde(av,msize,queue,S,weights,ESCF,ESCF_pur,ESCF_dbar)
!      use scf_stats
!      use fock_evaluator
!      use density_optimization
!      implicit none
!      integer, intent(in) :: msize
!      type(util_historyStore) :: queue
!      type(Matrix), intent(in) :: S
!      real(realk), intent(in) :: weights(msize)
!      real(realk), intent(out) :: ESCF,ESCF_pur, ESCF_dbar
!      type(Matrix) :: Didem,Dbar,Dtilde
!      integer :: ndim
!
!      ndim = S%nrow
!      call mat_init(Dbar,ndim,ndim)
!      call mat_init(Dtilde,ndim,ndim)
!      call mat_init(Didem,ndim,ndim)
!
!      !get Dbar and make it idempotent
!      call get_AVERAGE_arr(av,'D',queue,dsm_history_size,weights,Dbar)
!
!      call GET_Didem(Dbar,S,Didem) 
!      !** Find SCF energy for purified Dbar
!      call fck_get_fock(Didem,Dtilde,ESCF_dbar) !Dtilde is scratch
!
!      !get the corresponding Dtilde
!      call get_Dtilde(S,Dbar,Dtilde)
!      call GET_Didem(Dtilde,S,Didem)  
!      !** Find SCF energy for purified Dtilde 
!      call fck_get_fock(Didem,Dbar,ESCF_pur) !Dbar is scratch
!      !** Find SCF energy for Dtilde 
!      call fck_get_fock(Dtilde,Dbar,ESCF) !Dbar is scratch
!
!      call mat_free(Dbar)
!      call mat_free(Dtilde)
!      call mat_free(Didem)
!
!   end subroutine get_ESCF_dtilde

!> \brief Test routine for comparing the change in HF energy to the DSM energy change
!   subroutine COMPARE_EHF_EDSM(nqdsm,S,totstp,queue,his_start,DSMener,Dav)
!      use scf_stats
!      use fock_evaluator
!      use density_optimization
!      implicit none
!      integer, intent(in) :: nqdsm
!      type(Matrix), intent(in) :: S,Dav
!      real(realk), intent(in) :: totstp, DSMener
!      TYPE(util_HistoryStore)    :: queue
!      integer, intent(in) :: his_start
!      type(Matrix) :: Didem,Fcur
!      real(realk) :: EHFcur
!      integer :: itfck
!
!      call mat_init(Didem,Dav%Nrow,Dav%Ncol)
!      call mat_init(Fcur,Dav%Nrow,Dav%Ncol)
!      itfck = stat_current_iteration
!    !old style: stat_tab(itfck,1) = EHF_start; 
!    !new style: queue%Energy(his_start) = EHF_start, stat_tab(itfck+1,5) = EDSM_start
!      call GET_Didem(Dav,S,Didem)  
!      !** Find new Fock matrix and SCF energy
!      call fck_get_fock(Didem,Fcur,EHFcur)
!      WRITE(LUPRI,'(i4,2f20.10,"    dsmit,dEDSM,dEHF")')&
!            & nqdsm  , DSMener - queue%Energy(his_start), EHFcur - queue%Energy(his_start)
!!           & totstp,DSMener-stat_tab(itfck+1,5),EHFcur-queue%Energy(his_start)
!      call mat_free(Didem)
!      call mat_free(Fcur)
!
!   end subroutine COMPARE_EHF_EDSM

!   subroutine SET_NOCC(D,S,nocc)
!      TYPE(Matrix), INTENT(IN) :: D, S
!      integer, intent(out)     :: nocc
!      !** Nocc is the number of occupied orbitals
!      nocc = NINT(mat_dotproduct(D,S))
!   end subroutine SET_NOCC

!> \brief Delta = 3DavSDav - 2DavSDavSDav - Dav
   subroutine eval_delta(Dav,S,Delta)
     implicit none
     type(matrix), intent(in) :: Dav,S
     type(Matrix), intent(inout) :: Delta
     type(matrix) :: Sd,DSD
     integer :: ndim

     ndim = S%nrow
      call mat_init(SD,ndim,ndim)
      call mat_init(DSD,ndim,ndim)
      !Delta = 3DavSDav - 2DavSDavSDav - Dav
      !call cpu_time(timing1)
      call mat_mul(S,Dav,'n','n',1E0_realk,0E0_realk,SD)
      call mat_mul(Dav,SD,'n','n',1E0_realk,0E0_realk,DSD)
      call mat_mul(DSD,SD,'n','n',-2E0_realk,0E0_realk,Delta)
      !call cpu_time(timing2)
      !cfg_dsmtime = cfg_dsmtime + (timing2 - timing1)
      call mat_daxpy(3E0_realk,DSD,Delta)
      call mat_daxpy(-1E0_realk,Dav,Delta)

      call mat_free(SD)
      call mat_free(DSD)
   end subroutine eval_delta
 
!> \brief Set the global variable his_p, keeping track of current position in queue etc.
   subroutine set_his_p(av,msize,his_pos)

     implicit none
     type(avItem), intent(in) :: av
     integer, intent(in) :: msize,his_pos
     integer :: i,j,k

     j = his_pos
     i = av%dsm_pos
     do k = 1,msize
       his_p(i) = j
       j = j-1 
       if (j==0) then 
         j=av%cfg_settings(av%cfg_set_type)%max_history_size
       endif
       i = i-1
       if(i==0) then 
         i=av%dsm_history_size
       endif
     enddo

   end subroutine set_his_p
!   subroutine get_fd_hes_grad(msize,minstart,his_start,queue,S,H1,coef,diff,grad_out,hes_out)
!      
!      implicit none
!      integer, intent(in) :: msize,minstart,his_start
!      type(util_HistoryStore), intent(in) :: queue
!      type(matrix), intent(in) :: S,H1
!      real(realk), intent(in) :: coef(msize-1),diff
!      real(realk), intent(out) :: grad_out(msize-1),hes_out(msize-1,msize-1)
!
!      type(matrix) :: Dav,Fav,Delta
!      real(realk) :: grad(msize), hes(msize,msize),weights(msize),wrk(msize),grad_ps(msize)
!      real(realk) :: fd(msize),fmd(msize),f0,fdd(msize,msize),fmdmd(msize,msize),&
!                   & f2d(msize),fm2d(msize),f2d2d(msize,msize),fm2dm2d(msize,msize)
!      integer :: i,j,i1,j1,ndim
!      ndim = S%nrow
!
!      call mat_init(Dav,ndim,ndim)
!      call mat_init(Fav,ndim,ndim)
!      call mat_init(Delta,ndim,ndim)
!
!      j = 0
!      do i = 1,msize
!        if (i == minstart) then
!          weights(i) = 1.0E0_realk - SUM(coef)
!        else
!          j = j + 1
!          weights(i) = coef(j)
!        endif
!      enddo
!!** f0
!      call get_AVERAGE_arr('D',queue,dsm_history_size,weights,Dav)
!      call get_AVERAGE_arr('F',queue,dsm_history_size,weights,Fav)
!      call eval_delta(Dav,S,Delta)
!      if (cfg_dsm_app == cfg_dsm_xtra_term) then
!        call xterm_build_Gpara(msize,queue,Delta,S,H1)
!      endif
!      call DSMenergy(msize,his_start,queue,S,H1,fav,dav,Delta,f0) 
!      
!
!      do i = 1,msize
!        wrk = weights
!!** fd
!        wrk(i) = weights(i) + diff
!        call get_AVERAGE_arr('D',queue,dsm_history_size,wrk,Dav)
!        call get_AVERAGE_arr('F',queue,dsm_history_size,wrk,Fav)
!        call eval_delta(Dav,S,Delta)
!        if (cfg_dsm_app == cfg_dsm_xtra_term) then
!          call xterm_build_Gpara(msize,queue,Delta,S,H1)
!          grad_ps(i) = xterm_Eps(ndim,msize,queue,Delta,S,H1)
!        endif
!        call DSMenergy(msize,his_start,queue,S,H1,fav,dav,Delta,fd(i))
!!** fmd
!        wrk(i) = weights(i) - diff
!        call get_AVERAGE_arr('D',queue,dsm_history_size,wrk,Dav)
!        call get_AVERAGE_arr('F',queue,dsm_history_size,wrk,Fav)
!        call eval_delta(Dav,S,Delta)
!        if (cfg_dsm_app == cfg_dsm_xtra_term) then
!          call xterm_build_Gpara(msize,queue,Delta,S,H1)
!          grad_ps(i) = grad_ps(i) - xterm_Eps(ndim,msize,queue,Delta,S,H1)
!          grad_ps(i) = grad_ps(i)/(2.0E0_realk*diff)
!        endif
!        call DSMenergy(msize,his_start,queue,S,H1,fav,dav,Delta,fmd(i))
!!** f2d
!        wrk(i) = weights(i) + 2E0_realk*diff
!        call get_AVERAGE_arr('D',queue,dsm_history_size,wrk,Dav)
!        call get_AVERAGE_arr('F',queue,dsm_history_size,wrk,Fav)
!        call eval_delta(Dav,S,Delta)
!        if (cfg_dsm_app == cfg_dsm_xtra_term) then
!          call xterm_build_Gpara(msize,queue,Delta,S,H1)
!        endif
!        call DSMenergy(msize,his_start,queue,S,H1,fav,dav,Delta,f2d(i))
!!** fm2d
!        wrk(i) = weights(i) - 2E0_realk*diff
!        call get_AVERAGE_arr('D',queue,dsm_history_size,wrk,Dav)
!        call get_AVERAGE_arr('F',queue,dsm_history_size,wrk,Fav)
!        call eval_delta(Dav,S,Delta)
!        if (cfg_dsm_app == cfg_dsm_xtra_term) then
!          call xterm_build_Gpara(msize,queue,Delta,S,H1)
!        endif
!        call DSMenergy(msize,his_start,queue,S,H1,fav,dav,Delta,fm2d(i))
!!** grad
!!        grad(i) = (fd(i) - fmd(i))/(2.0E0_realk*diff)
!        grad(i) = (8.0E0_realk*fd(i)-8.0E0_realk*fmd(i)-f2d(i)+fm2d(i))/(12.0E0_realk*diff)
!!** hes
!!        hes(i,i) = (fd(i) + fmd(i) - 2.0E0_realk*f0)/(diff*diff)
!        hes(i,i) = (-f2d(i) + 16.0E0_realk*fd(i) - 30.0E0_realk*f0 + 16.0E0_realk*fmd(i) - fm2d(i))/(12.0E0_realk*diff*diff)
!        do j = i+1,msize
!!** fdd
!          wrk(i) = weights(i) + diff
!          wrk(j) = weights(j) + diff
!          call get_AVERAGE_arr('D',queue,dsm_history_size,wrk,Dav)
!          call get_AVERAGE_arr('F',queue,dsm_history_size,wrk,Fav)
!          call eval_delta(Dav,S,Delta)
!          if (cfg_dsm_app == cfg_dsm_xtra_term) then
!            call xterm_build_Gpara(msize,queue,Delta,S,H1)
!          endif
!          call DSMenergy(msize,his_start,queue,S,H1,fav,dav,Delta,fdd(i,j)) 
!!** fmdmd
!          wrk(i) = weights(i) - diff
!          wrk(j) = weights(j) - diff
!          call get_AVERAGE_arr('D',queue,dsm_history_size,wrk,Dav)
!          call get_AVERAGE_arr('F',queue,dsm_history_size,wrk,Fav)
!          call eval_delta(Dav,S,Delta)
!          if (cfg_dsm_app == cfg_dsm_xtra_term) then
!            call xterm_build_Gpara(msize,queue,Delta,S,H1)
!          endif
!          call DSMenergy(msize,his_start,queue,S,H1,fav,dav,Delta,fmdmd(i,j))
!!** f2d2d
!          wrk(i) = weights(i) + 2E0_realk*diff
!          wrk(j) = weights(j) + 2E0_realk*diff
!          call get_AVERAGE_arr('D',queue,dsm_history_size,wrk,Dav)
!          call get_AVERAGE_arr('F',queue,dsm_history_size,wrk,Fav)
!          call eval_delta(Dav,S,Delta)
!          if (cfg_dsm_app == cfg_dsm_xtra_term) then
!            call xterm_build_Gpara(msize,queue,Delta,S,H1)
!          endif
!          call DSMenergy(msize,his_start,queue,S,H1,fav,dav,Delta,f2d2d(i,j)) 
!!** fm2dm2d
!          wrk(i) = weights(i) - 2E0_realk*diff
!          wrk(j) = weights(j) - 2E0_realk*diff
!          call get_AVERAGE_arr('D',queue,dsm_history_size,wrk,Dav)
!          call get_AVERAGE_arr('F',queue,dsm_history_size,wrk,Fav)
!          call eval_delta(Dav,S,Delta)
!          if (cfg_dsm_app == cfg_dsm_xtra_term) then
!            call xterm_build_Gpara(msize,queue,Delta,S,H1)
!          endif
!          call DSMenergy(msize,his_start,queue,S,H1,fav,dav,Delta,fm2dm2d(i,j))
!        enddo
!      enddo
!  
!      do i = 1,msize
!        do j = i+1,msize
!          hes(i,j) = (2E0_realk*f0 + fdd(i,j) + fmdmd(i,j) - fd(i) - fmd(i) - fd(j) - fmd(j))/(2E0_realk*diff*diff)
!!          hes(i,j) = (f2d(i) + fm2d(i) + f2d(j) + fm2d(j) - 2E0_realk*f0 - f2d2d(i,j) - fm2dm2d(i,j))/(2E0_realk*diff*diff)
!!          hes(i,j) = hes(i,j) - 4E0_realk*(fd(i) + fmd(i) + fd(j) + fmd(j) - 2E0_realk*f0 - fdd(i,j) - fmdmd(i,j))/(2E0_realk*diff*diff)
!!          hes(i,j) = hes(i,j)/3E0_realk
!          hes(j,i) = hes(i,j)
!        enddo
!      enddo
!!
!! Since we are in the unconstrained framework......
!!
!      WRITE(LUPRI,*) 'grad_ps fd'
!      call OUTPUT(grad_ps,1,msize,1,1,msize,1,1,lupri)
!              WRITE(LUPRI,*) 'fd hes before....'
!         call OUTPUT(hes,1,msize,1,msize,msize,msize,1,lupri) 
!      i1 = 0
!      do i = 1,msize
!        if (i /= minstart) then
!          i1 = i1+1
!          grad_out(i1) = grad(i) - grad(minstart) 
!          j1 = 0
!          do j = 1,msize
!            if (j /= minstart) then
!              j1 = j1+1
!              hes_out(i1,j1) = hes(i,j) - hes(i,minstart) - hes(minstart,j) + hes(minstart,minstart)
!            endif
!          enddo
!        endif
!      enddo
!      call get_AVERAGE_arr('D',queue,dsm_history_size,weights,Dav)
!      call eval_delta(Dav,S,Delta)
!      if (cfg_dsm_app == cfg_dsm_xtra_term) then
!        call xterm_build_Gpara(msize,queue,Delta,S,H1)
!      endif
!
!      call mat_free(Dav)
!      call mat_free(Fav)
!      call mat_free(Delta)
!
!!     f0 = func(x)
!!     fd = func(x+delta)
!!     fmd = func(x-delta)
!!     f2d = func(x+2E0_realk*delta)
!!     fm2d = func(x-2E0_realk*delta)
! 
!!     grad2 = (fd-fmd)/(2.0E0_realk*delta)
!!     hes2 = (fd+fmd-2.0E0_realk*f0)/delta**2
!!     grad = (8.0E0_realk*fd-8.0E0_realk*fmd-f2d+fm2d)/(12.0E0_realk*delta)
!!     hes = (-1.0E0_realk*f2d+16.0E0_realk*fd-30.0E0_realk*f0+16.0E0_realk*fmd-1.0E0_realk*fm2d)/(12.0E0_realk*delta**2)
! 
!   end subroutine get_fd_hes_grad

end module LINSCA_DSM


