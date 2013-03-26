MODULE mm_stats_mod

   USE mm_global_paras_mod
   IMPLICIT NONE
   PRIVATE

   ! public procedures
   PUBLIC :: mm_init_stats,             &
             mm_init_buffer_stats,      &
             mm_init_contractor_stats,  &
             mm_init_test_stats,        &
             mm_print_stats
               
   ! public variables
   !-----------------

   REAL(REALK), POINTER, PUBLIC, SAVE :: stat_tpack_unique,      &
                                         stat_tpack_total,       &
                                         stat_tpack_chunks,      &
                                         stat_tvect_builds

   INTEGER, PUBLIC, SAVE :: stat_iteration,             &
                                  stat_1el_2el_or_full,       &
                                  stat_NBAST,                 &
                                  stat_raw_moms,              &
                                  stat_nuc_moms,              &
                                  stat_raw_moms_RHS,          &
                                  stat_pkd_moms_RHS,          &
                                  stat_screened_moms_RHS,     &
                                  stat_raw_moms_LHS,          &
                                  stat_pkd_moms_LHS,          &
                                  stat_min_box_size,          &
                                  stat_deepest_level,         &
                                  stat_common_branch

   REAL(REALK), POINTER, PUBLIC, SAVE :: stat_single_tests,   &
                                         stat_multi_tests

   LOGICAL, PUBLIC, SAVE :: stat_get_E_not_J,  &
                            stat_NN_not_FF,    &
                            stat_do_grad,      & ! do gradient
                            stat_reorder_mom     ! reorder the moments


   ! private variables
   !------------------

   REAL(REALK), TARGET, SAVE :: stat_T_chunks_NN,       &
                                stat_T_chunks_FF,       &
                                stat_W_chunks_NN,       &
                                stat_W_chunks_FF

   REAL(REALK), TARGET, SAVE :: stat_T_direction_NN,       &
                                stat_T_vector_NN,       &
                                stat_T_total_NN,        &
                                stat_T_direction_FF,    &
                                stat_T_vector_FF,       &
                                stat_T_total_FF
   REAL(REALK), TARGET, SAVE :: stat_W_direction_FF,       &
                                stat_W_total_FF

   REAL(REALK), TARGET, SAVE :: stat_single_tests_NN,   &
                                stat_single_tests_FF,   &
                                stat_multi_tests_NN,    &
                                stat_multi_tests_FF

!FIXME: come back to this (requires some detailed structures)
!                    stat_box_moms_RHS(:),       &   ! # boxed moments at level
!                    stat_box_moms_LHS(:),       &
!                    stat_branches(:),           &   ! # branches at each level
               


CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_stats(scheme,mode,gradmode)

      IMPLICIT NONE
      TYPE(scheme_paras), POINTER :: scheme
      LOGICAL,         INTENT(IN) :: mode, gradmode
   
      stat_get_E_not_J = mode
      stat_do_grad = gradmode
      stat_reorder_mom = .false.
   
   END SUBROUTINE mm_init_stats

!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_buffer_stats(mode)
   
      IMPLICIT NONE
      CHARACTER(1), INTENT(IN) :: mode

      SELECT CASE (mode)
         CASE ('T')
            IF (stat_NN_not_FF) THEN
               stat_tpack_chunks => stat_T_chunks_NN
               stat_tpack_unique => stat_T_direction_NN
               stat_tpack_total => stat_T_total_NN
            ELSE
               stat_tpack_chunks => stat_T_chunks_FF
               stat_tpack_unique => stat_T_direction_FF
               stat_tpack_total => stat_T_total_FF
            END IF
         CASE ('W')
            IF (stat_NN_not_FF) THEN
               CALL LSQUIT ('do not expect to translate moments in NN phase',-1)
            ELSE
               stat_tpack_chunks => stat_W_chunks_FF
               stat_tpack_unique => stat_W_direction_FF
               stat_tpack_total => stat_W_total_FF
            END IF
         CASE DEFAULT
            CALL LSQUIT('cannot reconcile buffer statistics requested',-1)
      END SELECT
   
   END SUBROUTINE mm_init_buffer_stats

!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_contractor_stats
   
      IMPLICIT NONE

      IF (stat_NN_not_FF) THEN
         stat_tvect_builds => stat_T_vector_NN
      ELSE
         stat_tvect_builds => stat_T_vector_FF
      END IF

   END SUBROUTINE mm_init_contractor_stats

!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_test_stats(phase)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: phase

      IF (phase == DO_NN) THEN
         stat_single_tests => stat_single_tests_NN
         stat_multi_tests => stat_multi_tests_NN
      ELSE
         stat_single_tests => stat_single_tests_FF
         stat_multi_tests => stat_multi_tests_FF
      END IF
      stat_single_tests = 0
      stat_multi_tests = 0

   END SUBROUTINE mm_init_test_stats

!-------------------------------------------------------------------------------

   SUBROUTINE mm_print_stats
   
      IMPLICIT NONE

      IF (.NOT. mm_stats_printed) THEN
         RETURN
      END IF

      WRITE(LUPRI,'(A,I4)') " deepest level  =", stat_deepest_level
      WRITE(LUPRI,'(/,A)') " MM statistics"
      WRITE(LUPRI,'(A,/)') " -------------"

      IF (stat_get_E_not_J) THEN
         WRITE(LUPRI,'(A)') " building multipole total energy only."
      ELSE
         SELECT CASE (stat_1el_2el_or_full)
            CASE (INTTYPE_ONE_EL)
               WRITE(LUPRI,*) "building 1-el J-matrix"
            CASE (INTTYPE_TWO_EL)
               WRITE(LUPRI,*) "building 2-el J-matrix"
               WRITE(LUPRI,'(A,I4)') " iteration =", stat_iteration
            CASE (INTTYPE_FULL_J)
               WRITE(LUPRI,*) "building full J-matrix"
               WRITE(LUPRI,'(A,I4)') " iteration =", stat_iteration
         END SELECT
      END IF

      WRITE(LUPRI,'(/,A,I9)') " number of nuclear moments  =", stat_nuc_moms
      WRITE(LUPRI,'(A,I9)') " number of raw moments      =", stat_raw_moms
      WRITE(LUPRI,'(A,20X,A,I9)') " RHS raw","=", stat_raw_moms_RHS
      WRITE(LUPRI,'(A,17X,A,I9)') " RHS packed","=", stat_pkd_moms_RHS
      WRITE(LUPRI,'(A,15X,A,I9)') " RHS screened","=", stat_screened_moms_RHS
      WRITE(LUPRI,'(A,20X,A,I9)') " LHS raw","=", stat_raw_moms_LHS
      WRITE(LUPRI,'(A,17X,A,I9)') " LHS packed","=", stat_pkd_moms_LHS

      WRITE(LUPRI,'(/,A)') " NN contraction pairs"
      WRITE(LUPRI,'(A,E11.3)') " total classical pairs  =", stat_T_total_NN
      WRITE(LUPRI,'(A,E11.3)') " T matrices built       =", stat_T_vector_NN
      WRITE(LUPRI,'(A,E11.3)') " total directions       =", stat_T_direction_NN
      WRITE(LUPRI,'(A,E11.3)') " total buffer chunks    =", stat_T_chunks_NN
      WRITE(LUPRI,'(/,A)') " FF contraction pairs"
      WRITE(LUPRI,'(A,E11.3)') " total classical pairs  =", stat_T_total_FF
      WRITE(LUPRI,'(A,E11.3)') " T matrices built       =", stat_T_vector_FF
      WRITE(LUPRI,'(A,E11.3)') " total directions       =", stat_T_direction_FF
      WRITE(LUPRI,'(A,E11.3)') " total buffer chunks    =", stat_T_chunks_FF

      WRITE(LUPRI,'(/,A)') " FF translation vectors"
      WRITE(LUPRI,'(A,E11.3)') " total number         =", stat_W_total_FF
      WRITE(LUPRI,'(A,E11.3)') " total directions     =", stat_W_direction_FF
      WRITE(LUPRI,'(A,E11.3)') " total buffer chunks  =", stat_W_chunks_FF

      WRITE(LUPRI,'(/,A)') " search statistics"
      WRITE(LUPRI,'(A,E11.3)') " NN single tests   = ", stat_single_tests_NN
      WRITE(LUPRI,'(A,E11.3)') " NN multiple tests = ", stat_multi_tests_NN
      WRITE(LUPRI,'(A,E11.3)') " FF single tests   = ", stat_single_tests_FF
      WRITE(LUPRI,'(A,E11.3)') " FF multiple tests = ", stat_multi_tests_FF

      WRITE(LUPRI,'(/,A,I5)') " deepest level  =", stat_deepest_level
      WRITE(LUPRI,'(A,I5,/)') " common branch  =", stat_common_branch
!      WRITE(LUPRI,'(/,A,I5,/)') " TMATM_DF  =", TMATM_DF

   END SUBROUTINE mm_print_stats

!-------------------------------------------------------------------------------

END MODULE mm_stats_mod
