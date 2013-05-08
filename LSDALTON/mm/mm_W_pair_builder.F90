MODULE mm_W_pair_builder

   USE mm_global_paras_mod
   USE mm_W_contractors, ONLY: mm_select_W_con

   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_translate_raw_moments,          &
             mm_translate_moments,              &
             mm_translate_parents_Vff,          &
             mm_get_raw_Vff_from_boxed_Vff

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_W_pair_builder(scheme)

      USE mm_W_buffers, ONLY: mm_open_W_buffer
      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      CALL mm_open_W_buffer(scheme)

   END SUBROUTINE mm_init_W_pair_builder

!-------------------------------------------------------------------------------

   SUBROUTINE mm_free_W_pair_builder(scheme)

      USE mm_W_buffers, ONLY: mm_close_W_buffer
      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      CALL mm_close_W_buffer(scheme)

   END SUBROUTINE mm_free_W_pair_builder

!-------------------------------------------------------------------------------

   SUBROUTINE mm_get_W_pair(addr,r_ab,new_LMAX,old_LMAX,object,W_pair)

      IMPLICIT NONE
      TYPE(old_new),       INTENT(IN)  :: addr 
      REAL(REALK),         INTENT(IN)  :: r_ab(3)
      INTEGER,       INTENT(IN)  :: new_LMAX, old_LMAX 
      CHARACTER(3),        INTENT(IN)  :: object 
      TYPE(T_pair_single), INTENT(OUT) :: W_pair

      W_pair%paras%ratio = one ! i.e. W_pair%r_ab is the actual vector

      ! indices to map back to actual moments in separate array
      W_pair%paras%LHS_id = addr%new
      W_pair%paras%RHS_id = addr%old
      ! orders of contraction with T-matrix
      W_pair%paras%LHS_LMAX = new_LMAX
      W_pair%paras%RHS_LMAX = old_LMAX
      
      ! (see mm_W_worker)
      SELECT CASE (object)
      CASE ('qlm')
         W_pair%r_ab(:) = r_ab(:) 
         W_pair%N_or_T  = 'N' 
      CASE ('Vff')
         ! For Vff translations, we require W'(-r_ab)
         W_pair%r_ab(:) = -r_ab(:) 
         ! Use 'T' (transpose) for W_matrix contraction with DTRMV
         W_pair%N_or_T  = 'T' 
      CASE DEFAULT
         CALL LSQUIT('cannot resolve translation object in mm_get_W_pair!',-1)
      END SELECT

      W_pair%LMAX   = MAX(W_pair%paras%LHS_LMAX,W_pair%paras%RHS_LMAX)
      W_pair%lm_max = (1+W_pair%LMAX)**2

   END SUBROUTINE mm_get_W_pair

!-------------------------------------------------------------------------------

   SUBROUTINE mm_translate_raw_moments(scheme,mms_in,mms_out)

      USE mm_W_contractors, ONLY: mm_set_W_con_ptrs
      USE mm_W_buffers,     ONLY: mm_add_to_W_buffer
      USE mm_qlm_processor, ONLY: mm_get_T_sym_qlm

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
!FIXME: do not need to pass all this in! (just the moments)
      TYPE(raw_mm_data),  INTENT(IN)    :: mms_in
      TYPE(box_mm_data),  INTENT(INOUT) :: mms_out

      TYPE(T_pair_single) :: W_pair
      TYPE(old_new)       :: addr 
      INTEGER       :: i, new_LMAX, old_LMAX 
      REAL(REALK)         :: new_centre(3), old_centre(3), trans_vec(3) 

      CALL mm_select_W_con(scheme%W_con%ID)
      old_LMAX = scheme%raw_LMAX
      new_LMAX = scheme%trans_LMAX
      ! set pointers in W_contractor
      CALL mm_set_W_con_ptrs(mms_in%qlm_W,mms_out%qlm_W)

      ! generate W_pairs and pass them to the W_buffer
      CALL mm_init_W_pair_builder(scheme)
      DO i = 1, SIZE(mms_in%paras) 
         ! indices to actual moments to translate from and to
         addr%old = mms_in%paras(i)%id
         addr%new = mms_in%paras(i)%map_up
         IF (addr%new == 0) CALL LSQUIT ('parameter mappings incomplete! 1',-1)
         old_centre = mms_in%paras(i)%cntr
         new_centre = mms_in%paras(i)%box_cntr
         ! optimising LMAX
         IF (scheme%dynamic_LMAX_on) THEN
            old_LMAX = MIN(mms_in%paras(i)%LMIN+scheme%LEXTRA, scheme%raw_LMAX)
         END IF
         ! make W_pair and throw into W-buffer
         trans_vec = new_centre - old_centre
         CALL mm_get_W_pair(addr,trans_vec,new_LMAX,old_LMAX,'qlm',W_pair)
         CALL mm_add_to_W_buffer(W_pair)
      END DO

      ! scale new translated moments for T_matrix symmetry
      ! but first ensure W_buffer is empty
      CALL mm_free_W_pair_builder(scheme)
      CALL mm_get_T_sym_qlm(new_LMAX,mms_out%qlm_W,mms_out%qlm_T)

   END SUBROUTINE mm_translate_raw_moments

!-------------------------------------------------------------------------------

   SUBROUTINE mm_translate_moments(scheme,mms_in,mms_out)

      USE mm_W_contractors, ONLY: mm_set_W_con_ptrs
      USE mm_W_buffers,     ONLY: mm_add_to_W_buffer
      USE mm_qlm_processor, ONLY: mm_get_T_sym_qlm

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      TYPE(box_mm_data),  INTENT(IN)    :: mms_in
      TYPE(box_mm_data),  INTENT(INOUT) :: mms_out

      TYPE(T_pair_single) :: W_pair
      TYPE(old_new)       :: addr 
      INTEGER       :: i, new_LMAX, old_LMAX 
      REAL(REALK)         :: new_centre(3), old_centre(3), trans_vec(3) 

      CALL mm_select_W_con(scheme%W_con%ID)
      old_LMAX = scheme%trans_LMAX
      new_LMAX = scheme%trans_LMAX
      ! set pointers in W_contractor
      CALL mm_set_W_con_ptrs(mms_in%qlm_W,mms_out%qlm_W)

      ! generate W_pairs and pass them to the W_buffer
      CALL mm_init_W_pair_builder(scheme)
      DO i = 1, SIZE(mms_in%RHS_paras) 
         ! indices to actual moments to translate from and to
         addr%old = mms_in%RHS_paras(i)%id
         addr%new = mms_in%RHS_paras(i)%map_up
         IF (addr%new == 0) CALL LSQUIT ('parameter mappings incomplete! 2',-1)
         old_centre = mms_in%RHS_paras(i)%cntr
         new_centre = mms_in%RHS_paras(i)%cntr_up
         ! make W_pair and throw into W-buffer
         trans_vec = new_centre - old_centre
         CALL mm_get_W_pair(addr,trans_vec,new_LMAX,old_LMAX,'qlm',W_pair)
         CALL mm_add_to_W_buffer(W_pair)
      END DO

      ! scale new translated moments for T_matrix symmetry
      ! but first ensure W_buffer is empty
      CALL mm_free_W_pair_builder(scheme)
      CALL mm_get_T_sym_qlm(new_LMAX,mms_out%qlm_W,mms_out%qlm_T)

   END SUBROUTINE mm_translate_moments

!-------------------------------------------------------------------------------

   SUBROUTINE mm_translate_parents_Vff(level,scheme,Vff_p,Vff_c,c_box_paras)

      USE mm_W_contractors, ONLY: mm_set_W_con_ptrs
      USE mm_W_buffers,     ONLY: mm_add_to_W_buffer

      IMPLICIT NONE
      INTEGER,      INTENT(IN)    :: level
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      REAL(REALK),        INTENT(IN)    :: Vff_p(:,:)
      REAL(REALK),        INTENT(INOUT) :: Vff_c(:,:)
      TYPE(box_mm_paras), INTENT(IN)    :: c_box_paras(:)   ! child

      TYPE(T_pair_single) :: W_pair
      TYPE(old_new)       :: addr 
      INTEGER       :: i, new_LMAX, old_LMAX
      REAL(REALK)         :: new_centre(3), old_centre(3), trans_vec(3)

      IF (level <= 2) RETURN

      new_LMAX = scheme%trans_LMAX
      old_LMAX = scheme%trans_LMAX

      CALL mm_select_W_con(scheme%W_con%ID)
      CALL mm_set_W_con_ptrs(Vff_p,Vff_c)
      CALL mm_init_W_pair_builder(scheme)
      DO i = 1, SIZE(Vff_c,2)
         addr%new = c_box_paras(i)%id
         addr%old = c_box_paras(i)%map_up
         IF (addr%old == 0) CALL LSQUIT ('parameter mappings incomplete! 3',-1)
         new_centre = c_box_paras(i)%cntr
         old_centre = c_box_paras(i)%cntr_up
         trans_vec  = new_centre - old_centre
         CALL mm_get_W_pair(addr,trans_vec,new_LMAX,old_LMAX,'Vff',W_pair)
         CALL mm_add_to_W_buffer(W_pair)
      END DO
      CALL mm_free_W_pair_builder(scheme)

   END SUBROUTINE mm_translate_parents_Vff

!-------------------------------------------------------------------------------

   SUBROUTINE mm_get_raw_Vff_from_boxed_Vff(raw_paras,scheme,boxed_Vff,Vff)

      USE mm_W_contractors, ONLY: mm_set_W_con_ptrs
      USE mm_W_buffers,     ONLY: mm_add_to_W_buffer

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(IN)    :: raw_paras(:)
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      REAL(REALK),        INTENT(IN)    :: boxed_Vff(:,:)
      REAL(REALK),        INTENT(INOUT) :: Vff(:,:)

      TYPE(T_pair_single) :: W_pair
      TYPE(old_new)       :: addr
      INTEGER       :: i, new_LMAX, old_LMAX
      REAL(REALK)         :: new_centre(3), old_centre(3), trans_vec(3)

      new_LMAX = scheme%raw_LMAX
      old_LMAX = scheme%trans_LMAX

      CALL mm_select_W_con(scheme%W_con%ID)
      CALL mm_set_W_con_ptrs(boxed_Vff,Vff)
      CALL mm_init_W_pair_builder(scheme)
      DO i = 1, SIZE(raw_paras)
         addr%new = raw_paras(i)%id
         addr%old = raw_paras(i)%map_up
         IF (addr%old == 0) CALL LSQUIT ('parameter mappings incomplete! 4',-1)
         IF (scheme%dynamic_LMAX_on) THEN
            new_LMAX = raw_paras(i)%LMIN + scheme%LEXTRA
            new_LMAX = MIN(scheme%raw_LMAX , raw_paras(i)%LMIN+scheme%LEXTRA)
         END IF
         new_centre = raw_paras(i)%cntr
         old_centre = raw_paras(i)%box_cntr
         trans_vec = new_centre - old_centre
         CALL mm_get_W_pair(addr,trans_vec,new_LMAX,old_LMAX,'Vff',W_pair)
         CALL mm_add_to_W_buffer(W_pair)
      END DO
      CALL mm_free_W_pair_builder(scheme)

   END SUBROUTINE mm_get_raw_Vff_from_boxed_Vff

!-------------------------------------------------------------------------------

END MODULE mm_W_pair_builder
