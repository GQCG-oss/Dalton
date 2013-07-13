! Module containing routines to generate contributions to a
! multipole potential on the fly.
! Note that another module allows contraction of T-pairs to generate
! an energy or J-matrix directly (e.g.to exploit symmetry).


!FIXME:  ************************************************************
!FIXME:  *                                                          *
!FIXME:  *  Array bounds when W-vector is ZERO !!!!!!!!!!           *
!FIXME:  *                                                          *
!FIXME:  ************************************************************


MODULE mm_W_contractors

   USE mm_global_paras_mod
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_init_W_contractors,           &
             mm_free_W_contractors,           &
             mm_select_W_con,                 &
             mm_set_W_con_ptrs

   ! Public variable to stop the resetting of W_con pointers with open T-buffer
   PUBLIC :: mm_lock_W_con
 
   REAL(REALK), POINTER, SAVE :: W_matrix(:,:)
   INTEGER,            SAVE :: WLDA
   ! Pointers to actual moments and potentials described elsewhere
   REAL(REALK), POINTER, SAVE :: old_ptr(:,:)
   REAL(REALK), POINTER, SAVE :: new_ptr(:,:)
   ! Diagnostic variables
   CHARACTER(11), SAVE :: W_con_stat 
   LOGICAL,       SAVE :: mm_lock_W_con

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_W_contractors(LMAX)
     use mm_mem
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LMAX

      call mem_alloc_fmm(W_matrix,(1+LMAX)**2,(1+LMAX)**2)
      WLDA = (1+LMAX)**2
      W_matrix = zero

   END SUBROUTINE mm_init_W_contractors

!-------------------------------------------------------------------------------

   SUBROUTINE mm_free_W_contractors
     use mm_mem

      IMPLICIT NONE
      call mem_dealloc_fmm(W_matrix)

   END SUBROUTINE mm_free_W_contractors

!-------------------------------------------------------------------------------

   SUBROUTINE mm_select_W_con(W_con_ID)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: W_con_ID
      EXTERNAL mm_store_t_contractor

      IF (.NOT.ASSOCIATED(W_matrix)) CALL LSQUIT('W_matrix not allocated!',-1)

      SELECT CASE (W_con_ID)
      CASE (W_CONTRACTOR_DIRECT)
         CALL mm_store_w_contractor(mm_W_contractor_direct)
      CASE (W_CONTRACTOR_B)
         CALL mm_store_w_contractor(mm_W_contractor_B)
      CASE (W_CONTRACTOR_FAST)
         CALL mm_store_w_contractor(mm_W_contractor_FAST)
      CASE DEFAULT
         CALL LSQUIT ('invalid W_contractor requested!',-1)
      END SELECT
      ! initialise diagnostics
      W_con_stat = 'initialised'
      W_con_stat = 'initialised'
      mm_lock_W_con = .FALSE.

   END SUBROUTINE mm_select_W_con

!-------------------------------------------------------------------------------

   SUBROUTINE mm_set_W_con_ptrs(old,new)

      IMPLICIT NONE
      REAL(REALK), TARGET, INTENT(IN) :: old(:,:), new(:,:)
   
      IF (W_con_stat /= 'initialised') STOP 'no W_contractor preselected!' 
      IF (mm_lock_W_con) STOP 'W_buffer not empty! Cannot reset W_con!'
      NULLIFY (old_ptr, new_ptr)
      old_ptr => old(:,:)
      new_ptr => new(:,:)

   END SUBROUTINE mm_set_W_con_ptrs

!-------------------------------------------------------------------------------

   SUBROUTINE mm_check_W_status

      IMPLICIT NONE
      IF ((.NOT.ASSOCIATED(old_ptr)) .OR. (.NOT.ASSOCIATED(new_ptr))) THEN
         CALL LSQUIT('W_contractor pointers not associated as reqd.',-1)
      END IF

   END SUBROUTINE mm_check_W_status

!-------------------------------------------------------------------------------
! Built for DIRECT scheme.  Can only take one W_pair at a time.

   SUBROUTINE mm_W_contractor_direct(W_pair)

      USE mm_W_worker, ONLY: mm_get_ltsqr_W_matrix,       &
                             mm_contract_Wq

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: W_pair

      REAL(REALK)   :: arr_tmp(W_pair%lm_max)
      REAL(REALK)   :: r_mod
      INTEGER :: n,m,w, lm_dim, hi, p,q, LMAX,JMAX
      CHARACTER(1)  :: NT

      CALL mm_check_W_status

      NT = W_pair%N_or_T
      p = W_pair%paras%LHS_id
      q = W_pair%paras%RHS_id
      lm_dim = W_pair%lm_max

      r_mod = SQRT(DOT_PRODUCT(W_pair%r_ab,W_pair%r_ab))
      IF (r_mod > ZERO_VECT_TOL) THEN
         LMAX = W_pair%paras%LHS_LMAX
         JMAX = W_pair%paras%RHS_LMAX
         CALL mm_get_ltsqr_W_matrix(LMAX,JMAX,W_pair%r_ab,W_matrix)
         IF (LMAX /= JMAX) THEN
            n = (1+JMAX)**2
            m = (1+LMAX)**2
            w = W_pair%lm_max 
            CALL mm_contract_Wq(NT,W_matrix,WLDA,old_ptr(:,q),n,new_ptr(:,p),m)
         ELSE
            arr_tmp(:lm_dim) = old_ptr(:lm_dim,q)
            CALL DTRMV('L',NT,'U',lm_dim,W_matrix,WLDA,arr_tmp,1)
            new_ptr(:lm_dim,p) = new_ptr(:lm_dim,p) + arr_tmp(:lm_dim)
         END IF
      ELSE
         hi = (1+W_pair%paras%LHS_LMAX)**2
         new_ptr(:hi,p) = new_ptr(:hi,p) + old_ptr(:hi,q)
      END IF

   END SUBROUTINE mm_W_contractor_direct

!-------------------------------------------------------------------------------
   
   SUBROUTINE mm_W_contractor_B(W_pairs)

      USE mm_W_worker, ONLY: mm_get_ltsqr_W_matrix

      IMPLICIT NONE
      TYPE(T_pair_list), INTENT(IN) :: W_pairs
 
      REAL(REALK)   :: arr_tmp(W_pairs%lm_max)
      REAL(REALK)   :: lastlen, r_pq(3), r_pq_mod(3)
      INTEGER :: LMAX, i,j, p,q, lm_dim, hi
      CHARACTER(1)  :: NT
 
      CALL mm_check_W_status

      LMAX = W_pairs%LMAX
      NT = W_pairs%N_or_T
      r_pq_mod(:) = W_pairs%r_ab(:)
      lm_dim = W_pairs%lm_max 
      lastlen = zero

      DO i = 1, SIZE(W_pairs%paras) 

         p  = W_pairs%paras(i)%LHS_id
         q  = W_pairs%paras(i)%RHS_id

         hi = (1+W_pairs%paras(i)%RHS_LMAX)**2
         IF (hi < SIZE(arr_tmp)) arr_tmp = zero
         arr_tmp(:hi) = old_ptr(:hi,q)

         IF(ABS(W_pairs%paras(i)%ratio) < ZERO_VECT_TOL) THEN
            ! zero translation, so just copy over old to new
            hi = (1+W_pairs%paras(i)%LHS_LMAX)**2
            new_ptr(:hi,p) = new_ptr(:hi,p) + arr_tmp(:hi)
            lastlen = W_pairs%paras(i)%ratio
            CYCLE
         END IF
         ! get matrix if different
         IF(ABS(W_pairs%paras(i)%ratio-lastlen) > DISTINCT_T_TOL) THEN
            r_pq = r_pq_mod * W_pairs%paras(i)%ratio
            lastlen = W_pairs%paras(i)%ratio
            CALL mm_get_ltsqr_W_matrix(LMAX,LMAX,r_pq,W_matrix)
         END IF
         ! now translate
         hi = (1+W_pairs%paras(i)%RHS_LMAX)**2
         IF (hi < lm_dim) arr_tmp = zero
         arr_tmp(:hi) = old_ptr(:hi,q)
         CALL DTRMV('L',NT,'U',lm_dim,W_matrix,WLDA,arr_tmp,1)
         hi = (1+W_pairs%paras(i)%LHS_LMAX)**2
         new_ptr(:hi,p) = new_ptr(:hi,p) + arr_tmp(:hi)

      END DO

   END SUBROUTINE mm_W_contractor_B
 
!-------------------------------------------------------------------------------
   
   SUBROUTINE mm_W_contractor_FAST(W_pairs)

      USE mm_W_worker, ONLY: mm_get_ltsqr_W_matrix,       &
                             mm_contract_Wq

      IMPLICIT NONE
      TYPE(T_pair_list), INTENT(IN) :: W_pairs
 
      REAL(REALK)   :: arr_tmp(W_pairs%lm_max)
      REAL(REALK)   :: lastlen, r_pq(3), r_pq_mod(3)
      INTEGER :: LMAX,JMAX, n,m,lm_dim, i,j, p,q, hi
      CHARACTER(1)  :: NT
 
      CALL mm_check_W_status

      NT = W_pairs%N_or_T
      r_pq_mod(:) = W_pairs%r_ab(:)
      lastlen = zero
      LMAX = W_pairs%LHS_LMAX
      JMAX = W_pairs%RHS_LMAX
      n = (1+JMAX)**2
      m = (1+LMAX)**2
      lm_dim = W_pairs%lm_max 

      DO i = 1, SIZE(W_pairs%paras) 

         p  = W_pairs%paras(i)%LHS_id
         q  = W_pairs%paras(i)%RHS_id
         IF(ABS(W_pairs%paras(i)%ratio) < ZERO_VECT_TOL) THEN
            ! zero translation, so just copy over old to new
            hi = (1+W_pairs%paras(i)%LHS_LMAX)**2
            new_ptr(:hi,p) = new_ptr(:hi,p) + old_ptr(:hi,q)
            lastlen = W_pairs%paras(i)%ratio
            CYCLE
         END IF
         IF(ABS(W_pairs%paras(i)%ratio-lastlen) > DISTINCT_T_TOL) THEN
            r_pq = r_pq_mod * W_pairs%paras(i)%ratio
            lastlen = W_pairs%paras(i)%ratio
            CALL mm_get_ltsqr_W_matrix(LMAX,JMAX,r_pq,W_matrix)
         END IF
         IF (LMAX /= JMAX) THEN
            CALL mm_contract_Wq(NT,W_matrix,WLDA,old_ptr(:,q),n,new_ptr(:,p),m)
         ELSE
            hi = (1+W_pairs%paras(i)%RHS_LMAX)**2
            IF (hi < lm_dim) arr_tmp = zero
            arr_tmp(:hi) = old_ptr(:hi,q)
            CALL DTRMV('L',NT,'U',lm_dim,W_matrix,WLDA,arr_tmp,1)
            hi = (1+W_pairs%paras(i)%LHS_LMAX)**2
            new_ptr(:hi,p) = new_ptr(:hi,p) + arr_tmp(:hi)
         END IF

      END DO

   END SUBROUTINE mm_W_contractor_FAST
 
!-------------------------------------------------------------------------------

END MODULE mm_W_contractors

