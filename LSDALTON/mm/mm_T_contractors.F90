MODULE mm_T_contractors

! Module containing routines to generate contributions to a
! multipole potential on the fly.
! Note that another module allows contraction of T-pairs to generate
! an energy or J-matrix directly (e.g.to exploit symmetry).

   USE mm_stats_mod
   USE mm_global_paras_mod
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_init_T_contractors,          &
             mm_free_T_contractors,          &
             mm_select_T_con,                &
             mm_set_T_con_ptrs,              &
             mm_set_T_contractor_threshold

   ! Public variable to stop the resetting of T_con pointers with open T-buffer
   PUBLIC :: mm_lock_T_con

   REAL(REALK), ALLOCATABLE, SAVE :: T_matrix(:,:)
   INTEGER,            SAVE :: TLDA
   ! when building multiple T matrices together
   REAL(REALK), ALLOCATABLE, SAVE :: T_mats(:,:,:)

   ! Pointers to actual moments and potentials described elsewhere
   REAL(REALK), POINTER, SAVE :: Vff_ptr(:,:)
   REAL(REALK), POINTER, SAVE :: qlm_ptr(:,:)
   REAL(REALK), POINTER, SAVE :: Vff_ptr_local(:,:)

!  Deafult threshold
   REAL(REALK), SAVE :: THRESH_T = 1.0d-16

   ! Diagnostic variables
   CHARACTER(11), SAVE :: T_con_stat 
   LOGICAL,       SAVE :: mm_lock_T_con
!$OMP THREADPRIVATE(T_mats,T_matrix,Vff_ptr_local)

CONTAINS
  
   SUBROUTINE mm_set_T_contractor_threshold(threshold)
   implicit none
   REAL(realk) :: threshold
     THRESH_T = threshold
   END SUBROUTINE mm_set_T_contractor_threshold


!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_T_contractors(scheme,phase,LMAX)
     use mm_mem
     use precision, only: long
     USE memory_handling, only: mem_allocated_mem_real
      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER,      INTENT(IN) :: LMAX, phase
      INTEGER :: lm_max, mode
      integer(kind=long) :: nsize
      lm_max = (LMAX+1)**2
      IF (scheme%dynamic_LMAX_on) CALL LSQUIT('T_cons not all test with DYNLMAX',-1) 
!FIXME: these allocations may fail if using DYNAMIC LMAX!

      mode = scheme%T_con%ID
      IF ( (phase == DO_FQ) .OR. (phase == DO_NN) ) mode = scheme%T_con%NN_ID

      IF ( mode == T_CONTRACTOR_MULTI ) THEN
!FIXME: will need to change TMATM_DF if changed elsewhere
! For some weird reason this does not work on ifort (combi of OMP,THREADPRIVATE,SAVE,POINTER)
!         call mem_alloc_fmm(T_mats,TMATM_DF,lm_max,lm_max)
         ALLOCATE(T_mats(TMATM_DF,lm_max,lm_max))
         nsize = SIZE(T_mats)
         CALL add_fmm_mem(nsize,6)
         call mem_allocated_mem_real(nsize)
         T_mats = zero
      ELSE
!         CALL MEM_ALLOC_FMM(T_matrix,lm_max,lm_max)
         ALLOCATE(T_matrix(lm_max,lm_max))
         nsize = SIZE(T_matrix)
         CALL add_fmm_mem(nsize,2)
         call mem_allocated_mem_real(nsize)
         T_matrix = zero
      END IF
!$OMP MASTER
      TLDA = lm_max
      
      ! code for diagmostic statistics only
      CALL mm_init_contractor_stats
      stat_tvect_builds = zero
!$OMP END MASTER

   END SUBROUTINE mm_init_T_contractors

!-------------------------------------------------------------------------------

   SUBROUTINE mm_free_T_contractors
     use mm_mem
     use precision, only: long
     USE memory_handling, only: mem_deallocated_mem_real
      IMPLICIT NONE
      integer(kind=long) :: nsize

      IF (ALLOCATED(T_matrix))THEN
         nsize = SIZE(T_matrix)
         CALL remove_fmm_mem(nsize,2)
         call mem_deallocated_mem_real(nsize)
         DEALLOCATE(T_matrix)
      ENDIF
      IF (ALLOCATED(T_mats))THEN
         nsize = SIZE(T_mats)
         CALL remove_fmm_mem(nsize,6)
         call mem_deallocated_mem_real(nsize)
         DEALLOCATE(T_mats)
      ENDIF

   END SUBROUTINE mm_free_T_contractors

!-------------------------------------------------------------------------------

   SUBROUTINE mm_select_T_con(T_con_ID)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: T_con_ID
      EXTERNAL mm_store_t_contractor

      SELECT CASE (T_con_ID)
      CASE (T_CONTRACTOR_DIRECT)
         CALL mm_store_t_contractor(mm_T_con_DIRECT)
      CASE (T_CONTRACTOR_TREE)
         CALL mm_store_t_contractor(mm_T_con_TREE)
      CASE (T_CONTRACTOR_DYN)
         CALL mm_store_t_contractor(mm_T_con_DYN)
      CASE (T_CONTRACTOR_SCALE_TREE)
         CALL mm_store_t_contractor(mm_T_con_SCALE_TREE)
      CASE (T_CONTRACTOR_SCALE)
         CALL mm_store_t_contractor(mm_T_con_SCALE)
      CASE (T_CONTRACTOR_MULTI)
         CALL mm_store_t_contractor(mm_T_con_MULTI)
      CASE (T_CONTRACTOR_FULL)
         CALL mm_store_t_contractor(mm_T_con_FULL)
      CASE DEFAULT
         CALL LSQUIT ('invalid T_contractor requested!',-1)
      END SELECT
      ! initialise diagnostics
      T_con_stat = 'initialised'
      mm_lock_T_con = .FALSE.

   END SUBROUTINE mm_select_T_con

!-------------------------------------------------------------------------------

   SUBROUTINE mm_set_T_con_ptrs(Vff,qlm)

      IMPLICIT NONE
      REAL(REALK), TARGET, INTENT(IN) :: Vff(:,:), qlm(:,:)
   
      IF (T_con_stat /= 'initialised') STOP 'no T_contractor preselected!' 
      IF (mm_lock_T_con) STOP 'T_buffer not empty! Cannot reset T_con!'
      NULLIFY (Vff_ptr, qlm_ptr)
      Vff_ptr => Vff(:,:)
      qlm_ptr => qlm(:,:)

   END SUBROUTINE mm_set_T_con_ptrs

!-------------------------------------------------------------------------------
! Built for DIRECT scheme.  Can only take one T_pair at a time.
! Only generates a square T-matrix with dim = MAX(LHS_LMAX,RHS_LMAX)
! Generates local potential to order MAX(LHS_LMAX,RHS_LMAX)
! Adds to Vff only up to LHS_LMAX.
! FIXME: because T_matrix is only built for (l+j)<=LMAX then with dynamic
! contraction we need LEXTRA>0 even for primitive basis sets and FQUAD.

   SUBROUTINE mm_T_con_DIRECT(T_pair)

      USE mm_T_worker, ONLY: mm_get_SPLTSQ_T_matrix,  &
                             mm_contract_Tq

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: T_pair

      REAL(REALK)   :: Vff_tmp(T_pair%lm_max)
      INTEGER :: weight, iLHS,iRHS, hi,l

      iRHS = T_pair%paras%RHS_id !iRHS was not defined here - I hope 
      ! it should be like this (defined further down)
      l = mm_set_T_con_LMAX(T_pair%LMAX,T_pair%r_ab,qlm_ptr(:T_pair%lm_max,iRHS))

      stat_tvect_builds = stat_tvect_builds + one
      CALL mm_get_SPLTSQ_T_matrix(l,T_pair%r_ab,T_matrix)

      iRHS = T_pair%paras%RHS_id
      hi = (l+1)*(l+1)

      Vff_tmp = mm_contract_Tq(T_pair%LMAX,qlm_ptr(:hi,iRHS),T_matrix)
      iLHS = T_pair%paras%LHS_id
      hi = (1+T_pair%paras%LHS_LMAX)**2
      weight = T_pair%paras%weight
!$OMP CRITICAL
      Vff_ptr(:hi,iLHS) = Vff_ptr(:hi,iLHS) + weight*Vff_tmp(:hi)
!$OMP END CRITICAL

   END SUBROUTINE mm_T_con_DIRECT

!-------------------------------------------------------------------------------
   
   SUBROUTINE mm_T_con_TREE(T_pairs)

      USE mm_T_worker, ONLY: mm_get_SPLTSQ_T_matrix,    &
                             mm_contract_Tq

      IMPLICIT NONE
      TYPE(T_pair_list), INTENT(IN) :: T_pairs
 
      REAL(REALK)   :: Vff_tmp(T_pairs%lm_max)
      REAL(REALK)   :: lastlen, r_pq(3), r_pq_mod(3)
      INTEGER :: LMAX, weight, i, p,q, hi,l
 
      r_pq_mod(:) = T_pairs%r_ab(:)
      LMAX = T_pairs%LMAX
      lastlen = zero

      DO i = 1, SIZE(T_pairs%paras) 

         p = T_pairs%paras(i)%LHS_id
         q = T_pairs%paras(i)%RHS_id

         IF(ABS(T_pairs%paras(i)%ratio-lastlen) > DISTINCT_T_TOL) THEN
            stat_tvect_builds = stat_tvect_builds + one
            r_pq = r_pq_mod * T_pairs%paras(i)%ratio
            CALL mm_get_SPLTSQ_T_matrix(LMAX,r_pq,T_matrix)
            lastlen = T_pairs%paras(i)%ratio
         END IF

         hi = (LMAX+1)*(LMAX+1)
         l = mm_set_T_con_LMAX(LMAX,r_pq_mod,qlm_ptr(:hi,q))
         hi = (l+1)*(l+1)

         Vff_tmp = mm_contract_Tq(l,qlm_ptr(:hi,q),T_matrix)
         hi = (1+T_pairs%paras(i)%LHS_LMAX)**2
         weight = T_pairs%paras(i)%weight
!$OMP CRITICAL
         Vff_ptr(:hi,p) = Vff_ptr(:hi,p) + weight*Vff_tmp(:hi)
!$OMP END CRITICAL

      END DO

   END SUBROUTINE mm_T_con_TREE
 
!-------------------------------------------------------------------------------
! Contractor designed for dynamic LMAX, using FULL T-matrix for
! low orders when block terms are important;
! exploit symmetry to only contract the lower triangular part
   
   SUBROUTINE mm_T_con_DYN(T_pairs)

      USE mm_T_worker, ONLY: mm_get_FLTSQ_T_matrix,    &
                             mm_postfac_Vff

      IMPLICIT NONE
      TYPE(T_pair_list), INTENT(IN) :: T_pairs
 
      REAL(REALK)   :: Vff_tmp(T_pairs%lm_max)
      REAL(REALK)   :: lastlen, r_pq(3), r_pq_mod(3)
      ! We pass to BLAS F77, hence set this as 
      DOUBLE PRECISION :: wt
      INTEGER :: LMAX, i, p,q, hi,l
 
      r_pq_mod(:) = T_pairs%r_ab(:)
      LMAX = T_pairs%LMAX
      lastlen = zero

      DO i = 1, SIZE(T_pairs%paras) 

         p = T_pairs%paras(i)%LHS_id
         q = T_pairs%paras(i)%RHS_id

         IF(ABS(T_pairs%paras(i)%ratio-lastlen) > DISTINCT_T_TOL) THEN
            r_pq = r_pq_mod * T_pairs%paras(i)%ratio
            CALL mm_get_FLTSQ_T_matrix(LMAX,r_pq,T_matrix)
            lastlen = T_pairs%paras(i)%ratio
         END IF

         hi = (LMAX+1)*(LMAX+1)
!         l = mm_set_T_con_LMAX(LMAX,r_pq_mod,qlm_ptr(:hi,q))
         l = LMAX
         hi = (l+1)*(l+1)
         wt = T_pairs%paras(i)%weight
         CALL DSYMV('L',hi,wt,T_matrix,TLDA,qlm_ptr(1,q),1,0D0,Vff_tmp,1)
!         hi = (1+T_pairs%paras(i)%LHS_LMAX)**2
         hi = (l+1)*(l+1)
         CALL mm_postfac_Vff(T_pairs%paras(i)%LHS_LMAX,Vff_tmp)
!$OMP CRITICAL
         Vff_ptr(:hi,p) = Vff_ptr(:hi,p) + Vff_tmp(:hi)
!$OMP END CRITICAL

      END DO

   END SUBROUTINE mm_T_con_DYN

!-------------------------------------------------------------------------------
   
   SUBROUTINE mm_T_con_SCALE(T_pairs)

      USE mm_T_worker, ONLY: mm_get_SPLTSQ_T_matrix,    &
                             mm_contract_Tq,            &
                             mm_scale_vec

      IMPLICIT NONE
      TYPE(T_pair_batch), INTENT(IN) :: T_pairs
 
      INTEGER :: LMAX, i, p,q, hi, lastq,l,lq
      REAL(REALK)   :: ratio, lastlen, pref, r_pq(3)
      REAL(REALK)   :: Vff_tmp(T_pairs%items(1)%lm_max)
      REAL(REALK)   :: scaled_qlm(T_pairs%items(1)%lm_max)
      REAL(REALK)   :: scale_vec(T_pairs%items(1)%lm_max)
      LOGICAL       :: new_vec

      pref = one
      lastq = -1
      lastlen = one
      scale_vec(:) = one
      LMAX = T_pairs%items(1)%LMAX
      stat_tvect_builds = stat_tvect_builds + one

      hi = (LMAX+1)*(LMAX+1)
      lq = LMAX
!      lq = 0
!      DO i = 1, T_pairs%ndim
!         q = T_pairs%items(i)%paras%RHS_id
!         lq = max(lq,mm_set_T_con_LMAX(LMAX,T_pairs%items(1)%r_ab,qlm_ptr(:hi,q)))
!      END DO
      CALL mm_get_SPLTSQ_T_matrix(lq,T_pairs%items(1)%r_ab,T_matrix)

      iloop: DO i = 1, T_pairs%ndim 

         p = T_pairs%items(i)%paras%LHS_id
         q = T_pairs%items(i)%paras%RHS_id

         hi = (lq+1)*(lq+1)
         l = LMAX
!         l = mm_set_T_con_LMAX(LMAX,T_pairs%items(i)%r_ab,qlm_ptr(:hi,q))
         hi = (l+1)*(l+1)

         new_vec = .FALSE.
         ratio = T_pairs%items(i)%paras%ratio
         IF (ABS(ratio-lastlen) > DISTINCT_T_TOL) THEN
            CALL mm_scale_vec(LMAX,ratio,scale_vec,pref)
            lastlen = ratio
            new_vec = .TRUE.
         END IF

         hi = T_pairs%items(i)%lm_max
         IF ( new_vec .OR. (q /= lastq) ) THEN
            scaled_qlm = scale_vec(:hi)*qlm_ptr(:hi,q)
            lastq = q
         END IF

         Vff_tmp = mm_contract_Tq(LMAX,scaled_qlm(:hi),T_matrix(:hi,:hi))

         hi = MIN((1+T_pairs%items(i)%paras%LHS_LMAX)**2,hi)
!$OMP CRITICAL
         Vff_ptr(:hi,p) = Vff_ptr(:hi,p) + pref*scale_vec(:hi)*Vff_tmp(:hi)
!$OMP END CRITICAL

      END DO iloop

   END SUBROUTINE mm_T_con_SCALE
 
!-------------------------------------------------------------------------------
   
   SUBROUTINE mm_T_con_SCALE_TREE(T_pairs)

      USE mm_T_worker, ONLY: mm_get_SPLTSQ_T_matrix,    &
                             mm_contract_Tq,            &
                             mm_scale_vec

      IMPLICIT NONE
      TYPE(T_pair_list), INTENT(IN) :: T_pairs
 
      INTEGER :: LMAX, i, p,q, hi, lastq
      REAL(REALK)   :: weight, lastlen, pref, r_pq(3)
      REAL(REALK)   :: Vff_tmp(T_pairs%lm_max)
      REAL(REALK)   :: scaled_qlm(T_pairs%lm_max)
      REAL(REALK)   :: scale_vec(T_pairs%lm_max)
      LOGICAL       :: new_vec
 
      pref = one
      lastq = -1
      lastlen = one
      scale_vec(:) = one
      LMAX = T_pairs%LMAX
      stat_tvect_builds = stat_tvect_builds + one
      CALL mm_get_SPLTSQ_T_matrix(LMAX,T_pairs%r_ab,T_matrix)

      iloop: DO i = 1, SIZE(T_pairs%paras) 

         p = T_pairs%paras(i)%LHS_id
         q = T_pairs%paras(i)%RHS_id

         new_vec = .FALSE.
         IF(ABS(T_pairs%paras(i)%ratio-lastlen) > DISTINCT_T_TOL) THEN
            CALL mm_scale_vec(LMAX,T_pairs%paras(i)%ratio,scale_vec,pref)
            lastlen = T_pairs%paras(i)%ratio
            new_vec = .TRUE.
         END IF

         hi = T_pairs%lm_max
         IF ( new_vec .OR. (q /= lastq) ) THEN
            scaled_qlm = scale_vec(:hi)*qlm_ptr(:hi,q)
            lastq = q
         END IF

         Vff_tmp = mm_contract_Tq(LMAX,scaled_qlm(:hi),T_matrix(:hi,:hi))

         hi = (1+T_pairs%paras(i)%LHS_LMAX)**2
         weight = pref*T_pairs%paras(i)%weight
!$OMP CRITICAL
         Vff_ptr(:hi,p) = Vff_ptr(:hi,p) + weight*scale_vec(:hi)*Vff_tmp(:hi)
!$OMP END CRITICAL

      END DO iloop

   END SUBROUTINE mm_T_con_SCALE_TREE

!-------------------------------------------------------------------------------
! Special contractor designed to take batch of T-pairs ordered in pairs
! such that ( b T1 1; a -T1 2; c T2 1; a -T2 3 ....)
! to halve the number of T matrix builds (using mm_scale...)
! But only the common RHS qlm can be contracted simultaneously which slows
! it down a lot.
!
!   SUBROUTINE mm_T_con_MULTI(T_pairs)
!
!      USE mm_T_worker_multi, ONLY: mm_get_SPLTSQ_T_matrices,    &
!                                   mm_contract_multi_Tq
!      USE mm_T_worker,       ONLY: mm_contract_Tq, mm_scale_vec
!
!      IMPLICIT NONE
!      TYPE(T_pair_batch), INTENT(IN) :: T_pairs
!
!      REAL(REALK), ALLOCATABLE :: Vff_tmp(:,:)
!
!      REAL(REALK), ALLOCATABLE :: scaled_qlm(:)
!      REAL(REALK), ALLOCATABLE :: scale_vec(:)
!
!      REAL(REALK)   :: T_vectors((T_pairs%ndim/2),3)
!      REAL(REALK)   :: pref, weight
!      INTEGER :: i,j, iLHS,iRHS,iRHS_last, hi, LMAX, nT
!
!      ! FIRST DO MULTIPLE T contracted with SAME RHS
!      !---------------------------------------------
!
!      IF (BTEST(T_pairs%ndim,0)) CALL LSQUIT('ndim not EVEN!')
!      nT = T_pairs%ndim/2
!
!      LMAX = 0
!      iRHS_last = T_pairs%items(1)%paras%RHS_id
!      DO i = 1, nT
!         j = 2*i-1 
!         LMAX = MAX(LMAX,T_pairs%items(i)%LMAX)
!         iRHS = T_pairs%items(j)%paras%RHS_id
!         ! get *distinct* T-vectors as every second item in batch
!         T_vectors(i,:) = T_pairs%items(j)%r_ab(:)
!         IF (iRHS /= iRHS_last) CALL LSQUIT('must have same qlm on RHS') 
!         iRHS_last = iRHS
!      END DO
!
!      stat_tvect_builds = stat_tvect_builds + nT
!      CALL mm_get_SPLTSQ_T_matrices(nT,LMAX,T_vectors,T_mats(:nT,:,:))
!
!      hi = (1+LMAX)**2
!      ALLOCATE( Vff_tmp(nT,hi) )
!      ALLOCATE( scaled_qlm(hi) )
!      ALLOCATE( scale_vec(hi) )
!      Vff_tmp(:,:) = zero
!      scale_vec(:) = one
!
!      iRHS = T_pairs%items(1)%paras%RHS_id
!
!      Vff_tmp(:,:hi) = mm_contract_multi_Tq(LMAX,qlm_ptr(:hi,iRHS),    &
!                                            T_mats(:nT,:,:),nT)
!
!      DO i = 1, nT
!         j = 2*i-1 
!         iLHS = T_pairs%items(j)%paras%LHS_id
!         hi = (1+T_pairs%items(j)%paras%LHS_LMAX)**2
!         weight = T_pairs%items(j)%paras%weight
!         Vff_ptr(:hi,iLHS) = Vff_ptr(:hi,iLHS) + weight*Vff_tmp(i,:hi)
!      END DO
!
!      ! NOW DO remaining half of batched T-pairs corresponding to same LHS
!      !-------------------------------------------------------------------
!
!      iLHS = T_pairs%items(2)%paras%LHS_id
!      DO i = 1, nT
!         ! get T_matrix correpxonding to the minus T-vector
!         CALL mm_scale_vec(LMAX,-one,scale_vec,pref)
!         iRHS = T_pairs%items(2*i)%paras%RHS_id
!         scaled_qlm = scale_vec(:hi)*qlm_ptr(:hi,iRHS)
!         Vff_tmp(1,:hi) = mm_contract_Tq(LMAX,scaled_qlm(:hi),T_mats(i,:,:))
!         weight = pref
!         Vff_ptr(:hi,iLHS) = Vff_ptr(:hi,iLHS)  &
!                             + weight*scale_vec(:hi)*Vff_tmp(1,:hi)
!      END DO
!
!      DEALLOCATE(Vff_tmp)
!      DEALLOCATE(scaled_qlm)
!      DEALLOCATE(scale_vec)
!
!   END SUBROUTINE mm_T_con_MULTI

!-------------------------------------------------------------------------------
! Contractor designed to build multiple T matrices simultaneously.
! Exact performance will depend on architecture and choice of NDIM.
! Assumes RHS moments are the same, but T-vectors are different.
 
   SUBROUTINE mm_T_con_MULTI(T_pairs)

      USE mm_T_worker_multi, ONLY: mm_get_SPLTSQ_T_matrices,    &
                                   mm_contract_multi_Tq

      IMPLICIT NONE
      TYPE(T_pair_batch), INTENT(IN) :: T_pairs

      REAL(REALK)   :: Vff_tmp(T_pairs%ndim,T_pairs%items(1)%lm_max)
      REAL(REALK)   :: T_vectors(T_pairs%ndim,3)
      INTEGER :: i, iLHS,iRHS, hi, LMAX, nT, l
      REAL(REALK), POINTER :: qlm(:)


      nT = T_pairs%ndim
      LMAX = T_pairs%items(1)%LMAX

!!$OMP CRITICAL
      DO i = 1, nT
         T_vectors(i,:) = T_pairs%items(i)%r_ab(:)
         iRHS = T_pairs%items(i)%paras%RHS_id
         IF ( iRHS /=  T_pairs%items(MAX(i-1,1))%paras%RHS_id ) THEN
            CALL LSQUIT('RHS moments not sorted in mm_T_con_MULTI',-1)
         END IF
      END DO
!!$OMP END CRITICAL

      stat_tvect_builds = stat_tvect_builds + nT
      CALL mm_get_SPLTSQ_T_matrices(nT,LMAX,T_vectors,T_mats)

      iRHS = T_pairs%items(1)%paras%RHS_id
      hi = (1+LMAX)**2

      Vff_tmp(:,:hi) = mm_contract_multi_Tq(LMAX,qlm_ptr(:hi,iRHS),T_mats,nT)

!$OMP CRITICAL
      DO i = 1, nT
         iLHS = T_pairs%items(i)%paras%LHS_id
         Vff_ptr(:hi,iLHS) = Vff_ptr(:hi,iLHS) + Vff_tmp(i,:hi)
      END DO
!$OMP END CRITICAL

   END SUBROUTINE mm_T_con_MULTI

!-------------------------------------------------------------------------------
! Builds full T-matrix for exact contraction at low orders

   SUBROUTINE mm_T_con_FULL(T_pair)

      USE mm_T_worker, ONLY: mm_get_FLTSQ_T_matrix, mm_postfac_Vff,  &
                             mm_contract_Tq

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: T_pair

      REAL(REALK)   :: Vff_tmp(T_pair%lm_max)
      INTEGER :: iLHS,iRHS, hi

      stat_tvect_builds = stat_tvect_builds + one
      CALL mm_get_FLTSQ_T_matrix(T_pair%LMAX,T_pair%r_ab,T_matrix)

      iRHS = T_pair%paras%RHS_id
      hi = T_pair%lm_max
      CALL DSYMV('L',hi,1d0,T_matrix,TLDA,qlm_ptr(1,iRHS),1,0D0,Vff_tmp,1)

      iLHS = T_pair%paras%LHS_id
      hi = (1+T_pair%paras%LHS_LMAX)**2
      CALL mm_postfac_Vff(T_pair%paras%LHS_LMAX,Vff_tmp)
!$OMP CRITICAL
      Vff_ptr(:hi,iLHS) = Vff_ptr(:hi,iLHS) + Vff_tmp(:hi)
!$OMP END CRITICAL

   END SUBROUTINE mm_T_con_FULL

!-------------------------------------------------------------------------------
   INTEGER FUNCTION mm_set_T_con_LMAX(LMAX,r_ab,qlm)

      IMPLICIT NONE
   
      INTEGER,intent(in) :: LMAX
      REAL(REALK),intent(in)   :: qlm((LMAX+1)**2)
      REAL(REALK),intent(in)   :: r_ab(3)
!
      INTEGER :: hi,l,m
      REAL(REALK)   :: r2i,Ish2(LMAX),qlm_max

      hi = (LMAX+1)**2

      r2i = 1.d0/(r_ab(1)*r_ab(1) + r_ab(2)*r_ab(2) + r_ab(3)*r_ab(3))

      Ish2(1) = r2i
      DO l=2,LMAX
        Ish2(l) = Ish2(l-1)*2*r2i*(2*l-1)
      ENDDO

      mm_set_T_con_LMAX = 1
      DO l=2,LMAX
        qlm_max = 0.d0
        DO m=l*l+1,(l+1)*(l+1)
          qlm_max = max(qlm_max,abs(qlm(m)))
        ENDDO
        qlm_max = qlm_max*qlm_max
        IF ((qlm_max*Ish2(l)*2).GT.THRESH_T*THRESH_T) mm_set_T_con_LMAX = l
      ENDDO
   END FUNCTION mm_set_T_con_LMAX

END MODULE mm_T_contractors

