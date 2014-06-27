MODULE mm_scheme_builder

   USE mm_global_paras_mod
   USE mm_stats_mod
   IMPLICIT NONE
   PRIVATE
   ! public procedures
   PUBLIC :: mm_get_scheme,                         &
             mm_set_main_algorithm,                 &
             mm_set_orders,                         &
             mm_set_search_and_test,                &
             mm_set_box_and_branch,                 &
             mm_set_contractors_and_buffers,        &
             mm_verify_scheme_consistency
 
   ! "scheme" contains all the input (and default) parameter data to 
   ! completely define a unique MM run-type;
   ! initialised on input, and saved for use on entry into main MM driver
   !
   TYPE(scheme_paras), TARGET, SAVE :: scheme

CONTAINS

!-------------------------------------------------------------------------------
      
   SUBROUTINE mm_get_scheme(scheme_ptr)

      IMPLICIT NONE
      TYPE(scheme_paras), POINTER :: scheme_ptr
      INTEGER :: iteration = -1

      NULLIFY(scheme_ptr)
      scheme_ptr => scheme 

      iteration = iteration + 1
      stat_iteration = iteration
      !print '(A,I4)', " running multipole moment code; iteration =", iteration
      IF (iteration == 0) CALL mm_print_scheme

   END SUBROUTINE mm_get_scheme

!-------------------------------------------------------------------------------

   SUBROUTINE mm_print_scheme

      IMPLICIT NONE

      WRITE(LUPRI,'(/,A)') " -----------------------------------------------"
      WRITE(LUPRI,'(A)')   " !  (Fast) Multipole Method input parameters   !"
      WRITE(LUPRI,'(A,/)') " -----------------------------------------------"

      SELECT CASE( scheme%algorithm )
         CASE( DO_NN )
            WRITE(LUPRI,*) "Running nearest neighbour MM algorithm."
         CASE( DO_FQ )
            WRITE(LUPRI,*) "Running Full Quadratic MM algorithm."
         CASE( DO_BQ )
            WRITE(LUPRI,*) "Running Box Quadratic MM algorithm."
         CASE( DO_NLOGN )
            WRITE(LUPRI,*) "Running O(NlogN) MM algorithm."
         CASE( DO_FMM )
            WRITE(LUPRI,*) "Running O(N) FMM algorithm."
         CASE DEFAULT
            CALL LSQUIT ('algorithm type not recognised!',-1)
      END SELECT

      IF (.NOT.scheme%inc_NN) WRITE(LUPRI,*) "Skipping all NN evaluations."

!FIXME: get rid of this triangular loop stuff.... OLD and useless!
!      IF ( ALL(scheme%T_searcher(:)%all_square) ) THEN
!         WRITE(LUPRI,*) "Forcing all loops to run square (not triangular)."
!      END IF

      SELECT CASE( scheme%T_con%T_buffer )
         CASE( TREE_T_BUFFER )
            WRITE(LUPRI,*) "Using Tree Buffer for FF T matrices."
         CASE( NULL_T_BUFFER )
            WRITE(LUPRI,*) "Building all FF T matrices on the fly."
         CASE( SKIP_T_BUFFER )
            WRITE(LUPRI,*) "Skipping all FF T matrix contractions."
         CASE( MULTI_T_BUFFER )
            WRITE(LUPRI,*) "Using buffer for multiple FF T matrix build."
         CASE( SCALE_T_BUFFER )
            WRITE(LUPRI,*) "Using buffer for scaled FF T matrix build."
         CASE DEFAULT
            CALL LSQUIT ('T-vector buffer not recognised!',-1)
      END SELECT

      SELECT CASE( scheme%T_con%NN_T_buffer )
         CASE( TREE_T_BUFFER )
            WRITE(LUPRI,*) "Using Tree Buffer for NN T matrices."
         CASE( NULL_T_BUFFER )
            WRITE(LUPRI,*) "Building all NN T matrices on the fly."
         CASE( SKIP_T_BUFFER )
            WRITE(LUPRI,*) "Skipping all NN T matrix contractions."
         CASE( MULTI_T_BUFFER )
            WRITE(LUPRI,*) "Using buffer for multiple NN T matrix build."
         CASE( SCALE_T_BUFFER )
            WRITE(LUPRI,*) "Using buffer for scaled NN T matrix build."
         CASE DEFAULT
            CALL LSQUIT ('T-vector buffer not recognised!',-1)
      END SELECT

      IF ( scheme%T_con%NN_ID == T_CONTRACTOR_FULL ) THEN
         WRITE(LUPRI,*) "Building full T-matrix for NN phase."
      END IF

      SELECT CASE( scheme%W_con%W_buffer )
         CASE( TREE_W_BUFFER )
            WRITE(LUPRI,*) "Using Tree Buffer for W matrices."
         CASE( NULL_W_BUFFER )
            WRITE(LUPRI,*) "Building all W matrices on the fly."
         CASE( SKIP_W_BUFFER )
            WRITE(LUPRI,*) "Skipping all W matrix contractions."
         CASE DEFAULT
            CALL LSQUIT ('W-vector buffer not recognised!',-1)
      END SELECT

      WRITE(LUPRI,'(A,I4)') " LMAX   =", scheme%raw_LMAX
      IF (scheme%algorithm /= DO_FQ) &
         WRITE(LUPRI,'(A,I4)') " TLMAX  =", scheme%trans_LMAX
      IF (scheme%dynamic_LMAX_on) THEN
         WRITE(LUPRI,*) "Running with dynamic LMAX for near field pairs."
         WRITE(LUPRI,'(A,I4)') " LEXTRA =", scheme%LEXTRA
      END IF

      IF (scheme%algorithm == DO_FQ) THEN 
         WRITE(LUPRI,'(/,A,/)') " End of MM input parameters."
         RETURN
      END IF

      IF (scheme%dynamic_levels) THEN
         WRITE(LUPRI,*) "Using dynamic NLEVEL."
         WRITE(LUPRI,'(A,F9.5)')" Smallest box size =", scheme%grain_input
      ELSE
         WRITE(LUPRI,*) "Using static NLEVEL."
         WRITE(LUPRI,'(A,I4)') " NLEVEL =", scheme%NLEVEL
      END IF
      IF (scheme%branch_free) THEN
         WRITE(LUPRI,*) "Running Branch-Free MM algorithm."
         WRITE(LUPRI,'(A,F9.5)')" Cut-off radius =", scheme%cls_radius
      END IF

      IF( ALL(scheme%T_searcher(:)%algorithm == FQ_LOOPS) ) THEN
         WRITE(LUPRI,*) "Using full quadratic search algorithm."
      ELSE IF( ANY(scheme%T_searcher(:)%algorithm == GRID_SEARCH) ) THEN
         WRITE(LUPRI,*) "Using grid search algorithm for NN phase."
         WRITE(LUPRI,*) "Using tree search algorithm for FF phase."
      ELSE IF( ANY(scheme%T_searcher(:)%algorithm == TREE_SEARCH) ) THEN
         WRITE(LUPRI,*) "Using binary tree search algorithm where possible."
      END IF

      IF (scheme%NN_box_pretesting) THEN
         WRITE(LUPRI,*) "Using pre-searching over boxes for near field pairs."
      END IF

!FIXME: put in tests to make ALL_SQUARE if using FULL_J option

      WRITE(LUPRI,'(/,A,/)') " End of MM input parameters."

   END SUBROUTINE mm_print_scheme

!-------------------------------------------------------------------------------

   SUBROUTINE mm_set_main_algorithm(ALGORITHM,INCNN,BRFREE,RPQMIN,  &
                                    SCREEN,CONTRACTED)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ALGORITHM
      LOGICAL,       INTENT(IN) :: BRFREE,INCNN,CONTRACTED
      REAL(REALK),   INTENT(IN) :: RPQMIN,SCREEN

      scheme%algorithm = ALGORITHM
      scheme%INC_NN = INCNN
      ! check for Branch-Free scheme
      scheme%branch_free = BRFREE
      scheme%cls_radius  = RPQMIN
      scheme%dens_screen_thr = SCREEN
      ! flag for use of internal decontraction of basis set
      scheme%contracted = CONTRACTED

      ! miscellaneous scheme parameters
      scheme%pack_LHS = PACK_LHS_DF
      scheme%pack_RHS = PACK_RHS_DF
      scheme%system_size = zero

      ! error checks
      !-------------

      SELECT CASE(ALGORITHM)
         CASE(DO_NN)
         CASE(DO_FQ)
         CASE(DO_BQ)
         CASE(DO_NLOGN)
         CASE(DO_FMM)
         CASE DEFAULT
            CALL LSQUIT ('algorithm type not recognised in mm_input.f90',-1)
      END SELECT
      IF (BRFREE) THEN
         IF (RPQMIN <= zero) THEN
            CALL LSQUIT ('RPQMIN selected for Branch-Free run invalid!',-1) 
         END IF
      END IF

   END SUBROUTINE mm_set_main_algorithm

!-------------------------------------------------------------------------------

   SUBROUTINE mm_set_orders(LMAX,TLMAX,DYNLMAX,LEXTRA)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LMAX, TLMAX, LEXTRA
      LOGICAL,       INTENT(IN) :: DYNLMAX

      ! for near-field contractions
      scheme%raw_LMAX = LMAX
      ! for translations and far-field contractions
      scheme%trans_LMAX = TLMAX
       ! optimise for individual symmetries of overlap distributions 
      scheme%dynamic_LMAX_on = DYNLMAX
      scheme%LEXTRA = LEXTRA

      ! error checks
      !-------------

      IF (LMAX < 0) THEN
          CALL LSQUIT ('LMAX for moment generation must be positive!',-1)
      ELSE IF (LMAX >= LMAX_HIGH) THEN
         WRITE(LUPRI,'(/,A)') "            **** WARNING ****"
         WRITE(LUPRI,'(A,/)') "  LMAX for multipole moments set very high!"
      END IF

      IF ( .NOT. (scheme%algorithm == DO_FQ) ) THEN
         IF (TLMAX < LMAX) THEN
            WRITE(LUPRI,'(/,A)') "TLMAX < LMAX  !!"
            WRITE(LUPRI,'(A,/)') "increase TLMAX or use default"
            CALL LSQUIT ('TLMAX for contractions and translations too low!',-1)
         END IF
         IF (TLMAX < 0) THEN
            CALL LSQUIT ('TLMAX for contractions and translations negative!',-1)
         END IF
         IF (TLMAX > TLMAX_HIGH) THEN
            WRITE(LUPRI,'(/,A)') "            **** WARNING ****"
            WRITE(LUPRI,'(A,/)') "  TLMAX for contractions set very high!"
         END IF
      END IF

      IF (DYNLMAX) THEN
         CALL LSQUIT ('Dynamic LMAX suspected to be broken!',-1)
!FIXME: must fix T_matrix allocation when using dynamic LMAX
         IF (LEXTRA < 0) THEN
            WRITE(LUPRI,'(/,A)') "          **** WARNING ****"
            WRITE(LUPRI,'(A,/)') "  LEXTRA for dynamic LMAX negative! "
         END IF
      END IF
   
   END SUBROUTINE mm_set_orders

!-------------------------------------------------------------------------------

   SUBROUTINE mm_set_search_and_test(ALL_SQUARE,TREE_TEST,NOBOXP,GRID_TEST)

      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: ALL_SQUARE,TREE_TEST,NOBOXP,GRID_TEST
   
      scheme%T_searcher(:)%all_square = ALL_SQUARE
      scheme%NN_box_pretesting = (.NOT.NOBOXP)

      scheme%T_searcher(:)%shape = SHAPE_SQUARE ! (non-square set when allowed)

      ! default quadratic loops for everything
      scheme%T_searcher(:)%algorithm = FQ_LOOPS 

      IF (TREE_TEST) THEN
         ! use local space binary tree searcher where we can
         scheme%T_searcher(DO_NN)%algorithm = TREE_SEARCH
         scheme%T_searcher(DO_FMM)%algorithm = TREE_SEARCH
         scheme%T_searcher(DO_NlogN)%algorithm = TREE_SEARCH
      END IF
      IF (GRID_TEST) THEN
         ! local space direct grid searcher only implemented for NN phase
         scheme%T_searcher(DO_NN)%algorithm = GRID_SEARCH
         scheme%NN_box_pretesting = .FALSE.  ! not strictly needed
      END IF

      IF (NOBOXP .AND. TREE_TEST) THEN
         CALL LSQUIT ('must use box pretesting with binary tree search!',-1)
      END IF

   END SUBROUTINE mm_set_search_and_test

!-------------------------------------------------------------------------------

   SUBROUTINE mm_set_box_and_branch(GRAIN,NLEVEL)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NLEVEL
      REAL(REALK),   INTENT(IN) :: GRAIN
   
      scheme%dynamic_levels = DYNLEV_DF
      scheme%grain_input = GRAIN
      scheme%NLEVEL = NLEVEL

      ! error checks
      !-------------
      ! note an input NLEVEL takes priority over an input GRAIN

      IF (.NOT. DYNLEV_DF) THEN
         CALL LSQUIT('input processing only designed for DYNLEV_DF = T',-1)
      END IF

      IF ( NLEVEL /= NLEVEL_DF ) THEN
         ! NLEVEL has been explicitly input in .dal file
         scheme%dynamic_levels = .FALSE.
         scheme%grain_input = -one   ! arbitrary
      END IF

      IF (scheme%dynamic_levels) THEN
         IF (scheme%grain_input <= zero) THEN
            CALL LSQUIT ('box GRAIN must be positive for dynamic NLEVEL!',-1)
         END IF
      ELSE IF ( scheme%NLEVEL < TOP_LEVEL ) THEN
         CALL LSQUIT ('NLEVEL input too low!',-1)
      END IF

   END SUBROUTINE mm_set_box_and_branch

!-------------------------------------------------------------------------------

   SUBROUTINE mm_set_contractors_and_buffers(DYNLMAX,T_CONTRACTOR)

      IMPLICIT NONE
      LOGICAL,       INTENT(IN) :: DYNLMAX
      INTEGER, INTENT(IN) :: T_CONTRACTOR

      ! all the below assume a symmetric T-matrix
      ! hence auxiliary (prefactor-adapted) moment types 
      scheme%T_con%LHS_mm_type = USE_RAW_QLM
      scheme%T_con%RHS_mm_type = USE_T_SYM_QLM

      ! *** W contractions ***
      scheme%W_con%W_buffer = TREE_W_BUFFER 
      scheme%W_con%ID = W_CONTRACTOR_FAST
      scheme%W_con%sort_para = SORT_BY_SCALE 
!FIXME: cannot currently change W-buffer on input
!      ! *** W contractions ***
!      scheme%W_con%W_buffer = NULL_W_BUFFER 
!      scheme%W_con%ID = W_CONTRACTOR_DIRECT
!      scheme%W_con%sort_para = NO_SORT 

      SELECT CASE( T_CONTRACTOR )
         CASE( 0 )
            ! we do no contractions (for diagnostics)
            scheme%T_con%T_buffer = SKIP_T_BUFFER 
            scheme%W_con%W_buffer = SKIP_W_BUFFER 
            scheme%T_con%NN_T_buffer = SKIP_T_BUFFER
            scheme%T_con%ID = T_CONTRACTOR_DIRECT
            scheme%T_con%NN_ID = T_CONTRACTOR_DIRECT
            ! remainder of the below fields are arbitrary
         CASE( 1 )
            scheme%T_con%T_buffer = NULL_T_BUFFER 
            scheme%T_con%NN_T_buffer = NULL_T_BUFFER
            scheme%T_con%ID = T_CONTRACTOR_DIRECT
            scheme%T_con%NN_ID = T_CONTRACTOR_DIRECT
         CASE( 2 )
            scheme%T_con%T_buffer = SCALE_T_BUFFER
            scheme%T_con%NN_T_buffer = MULTI_T_BUFFER
            scheme%T_con%ID = T_CONTRACTOR_SCALE
            scheme%T_con%NN_ID = T_CONTRACTOR_MULTI
         CASE( 3 )
            scheme%T_con%T_buffer = TREE_T_BUFFER 
            scheme%T_con%NN_T_buffer = TREE_T_BUFFER 
            scheme%T_con%sort_para = SORT_BY_SCALE 
            scheme%T_con%ID = T_CONTRACTOR_TREE
            scheme%T_con%NN_ID = T_CONTRACTOR_TREE
         CASE( 4 )
            scheme%T_con%T_buffer = TREE_T_BUFFER 
            scheme%T_con%NN_T_buffer = TREE_T_BUFFER 
            scheme%T_con%sort_para = SORT_BY_SCALE 
            scheme%T_con%ID = T_CONTRACTOR_SCALE_TREE
            scheme%T_con%NN_ID = T_CONTRACTOR_SCALE_TREE
         CASE( 5 )
            scheme%T_con%T_buffer = SCALE_T_BUFFER 
            scheme%T_con%NN_T_buffer = SCALE_T_BUFFER 
            scheme%T_con%ID = T_CONTRACTOR_SCALE
            scheme%T_con%NN_ID = T_CONTRACTOR_SCALE
         CASE( 6 )
            scheme%T_con%T_buffer = SCALE_T_BUFFER
            scheme%T_con%NN_T_buffer = NULL_T_BUFFER
            scheme%T_con%ID = T_CONTRACTOR_SCALE
            scheme%T_con%NN_ID = T_CONTRACTOR_FULL
            IF (scheme%raw_LMAX > FULL_T_LMAX ) THEN
               CALL LSQUIT('LMAX set very high for FULL T-matrix build',-1)
            END IF
         CASE( 7 )
            scheme%T_con%T_buffer = SCALE_T_BUFFER
            scheme%T_con%NN_T_buffer = NULL_T_BUFFER
            scheme%T_con%ID = T_CONTRACTOR_SCALE
            scheme%T_con%NN_ID = T_CONTRACTOR_DIRECT
         CASE DEFAULT
            CALL LSQUIT('invalid T-contractor specified!',-1)
      END SELECT

      IF ( T_CONTRACTOR == T_CONTRACTOR_DF ) THEN
         ! default T_contractor is used, which depends on
         ! whether we use internal deCONTRACTED.
         IF (.NOT.scheme%contracted) THEN
            ! use full T-matrix to ensure accuracy of low
            ! orders when using decontracted moments
            IF (scheme%raw_LMAX <= FULL_T_LMAX ) THEN
!               scheme%T_con%NN_T_buffer = NULL_T_BUFFER
!               scheme%T_con%NN_ID = T_CONTRACTOR_FULL
            END IF
         END IF
      END IF

      IF (DYNLMAX) THEN
         ! only the low order (NN and FQ) phases are done dynamically
         CALL LSQUIT('.DYNLMAX may be slow or even broken!',-1)
         scheme%T_con%NN_T_buffer = TREE_T_BUFFER 
         scheme%T_con%NN_ID = T_CONTRACTOR_DYN
      END IF

   
   END SUBROUTINE mm_set_contractors_and_buffers

!-------------------------------------------------------------------------------
   
   SUBROUTINE mm_verify_scheme_consistency

      IMPLICIT NONE

!FIXME: insert funky tests here

   END SUBROUTINE mm_verify_scheme_consistency
   
!-------------------------------------------------------------------------------

END MODULE mm_scheme_builder

!===============================================================================

SUBROUTINE mm_init_defs(NLEVEL,GRAIN,TLMAX,ALGORITHM,                         &
                        USEUMAT,BRFREE,RPQMIN,DYNLMAX,LEXTRA,                 &
                        ALLSQR,INCNN,TRSRCH,NOBOXP,GRSRCH,T_CONTRACTOR,SCREEN)
      
   USE mm_global_paras_mod
   IMPLICIT NONE
   INTEGER, INTENT(OUT) :: NLEVEL,TLMAX,LEXTRA,ALGORITHM,T_CONTRACTOR
   REAL(REALK),   INTENT(OUT) :: GRAIN,RPQMIN,SCREEN
   LOGICAL,       INTENT(OUT) :: USEUMAT,BRFREE,GRSRCH,INCNN
   LOGICAL,       INTENT(OUT) :: DYNLMAX,ALLSQR,TRSRCH,NOBOXP

   ! initialise defaults
   NLEVEL       = NLEVEL_DF
   GRAIN        = GRAIN_DF
   TLMAX        = TLMAX_DF
   ALGORITHM    = ALGORITHM_DF
   INCNN        = INC_NN_DF
   USEUMAT      = USEUMAT_DF
   BRFREE       = BRFREE_DF
   RPQMIN       = RPQMIN_DF
   DYNLMAX      = DYNLMAX_DF
   LEXTRA       = LEXTRA_DF
   ALLSQR       = ALLSQR_DF
   TRSRCH       = TRSRCH_DF
   NOBOXP       = NOBOXP_DF
   GRSRCH       = GRSRCH_DF
   T_CONTRACTOR = T_CONTRACTOR_DF
   SCREEN       = DENS_SCREEN_DF

END SUBROUTINE mm_init_defs

!-------------------------------------------------------------------------------
! routine to export named parameters used in MM code for standard input read-in

SUBROUTINE mm_get_named_paras(DO_FQ1,DO_NN1,DO_BQ1,DO_NLOGN1,DO_FMM1)

   USE mm_global_paras_mod
   IMPLICIT NONE
   INTEGER, INTENT(OUT) :: DO_FQ1,DO_NN1,DO_BQ1,DO_NLOGN1,DO_FMM1

   DO_FQ1    = DO_FQ
   DO_NN1    = DO_NN
   DO_BQ1    = DO_BQ
   DO_NLOGN1 = DO_NLOGN
   DO_FMM1   = DO_FMM

END SUBROUTINE mm_get_named_paras

!-------------------------------------------------------------------------------
! routine to initialise all run-type parameters with defaults

SUBROUTINE mm_set_def_scheme(LMAX)

   USE mm_global_paras_mod

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: LMAX

   INTEGER :: NLEVEL,TLMAX,LEXTRA,ALGORITHM,T_CONTRACTOR
   REAL(REALK)   :: GRAIN,RPQMIN,SCREEN
   LOGICAL       :: USEUMAT,BRFREE,GRSRCH,INCNN
   LOGICAL       :: DYNLMAX,ALLSQR,TRSRCH,NOBOXP,CONTRACTED

   CALL mm_init_defs(NLEVEL,GRAIN,TLMAX,ALGORITHM,                            &
                     USEUMAT,BRFREE,RPQMIN,DYNLMAX,LEXTRA,                    &
                     ALLSQR,INCNN,TRSRCH,NOBOXP,GRSRCH,T_CONTRACTOR,SCREEN)
   CONTRACTED = .TRUE. !this is never set for some reason?
   CALL mm_init_scheme(NLEVEL,GRAIN,LMAX,TLMAX,ALGORITHM,                     &
                       USEUMAT,BRFREE,RPQMIN,DYNLMAX,LEXTRA,                  &
                       ALLSQR,INCNN,TRSRCH,NOBOXP,GRSRCH,T_CONTRACTOR,        &
                       SCREEN,CONTRACTED)

END SUBROUTINE mm_set_def_scheme

!-------------------------------------------------------------------------------
! routine to (re)initialise scheme based on input keywords or defaults 

SUBROUTINE mm_init_scheme(NLEVEL,GRAIN,LMAX,TLMAX,ALGORITHM,              &
                          UMAT,BRFREE,RPQMIN,                             &
                          DYNLMAX,LEXTRA,ALL_SQUARE,INCNN,                &
                          TREE_TEST,NOBOXP,GRSRCH,T_CONTRACTOR,           &
                          SCREEN,CONTRACTED)

   USE mm_global_paras_mod
   USE mm_scheme_builder, ONLY: mm_set_main_algorithm,                 &
                                mm_set_orders,                         &
                                mm_set_search_and_test,                &
                                mm_set_box_and_branch,                 &
                                mm_set_contractors_and_buffers,        &
                                mm_verify_scheme_consistency

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: NLEVEL,LMAX,TLMAX,LEXTRA,ALGORITHM
   INTEGER, INTENT(IN) :: T_CONTRACTOR
   LOGICAL,       INTENT(IN) :: UMAT,BRFREE,DYNLMAX
   LOGICAL,       INTENT(IN) :: GRSRCH,INCNN
   LOGICAL,       INTENT(IN) :: ALL_SQUARE,TREE_TEST,NOBOXP,CONTRACTED
   REAL(REALK),   INTENT(IN) :: GRAIN,RPQMIN,SCREEN

   ! first take opportunity to link units with standard Dalton code
   IF ( .NOT. MM_STANDALONE ) LUPRI = 6 ! Stinne 26/10-2010 CALL mm_GTUNIT(LUPRI)

   CALL mm_set_main_algorithm(ALGORITHM,INCNN,BRFREE,RPQMIN,SCREEN,CONTRACTED)
   CALL mm_set_orders(LMAX,TLMAX,DYNLMAX,LEXTRA)
   CALL mm_set_search_and_test(ALL_SQUARE,TREE_TEST,NOBOXP,GRSRCH)
   CALL mm_set_box_and_branch(GRAIN,NLEVEL)
   CALL mm_set_contractors_and_buffers(DYNLMAX,T_CONTRACTOR)

!FIXME: Given this much, the remaining quantities could be generated,
!       in part, automatically??
!FIXME: this testing could be made simpler by using automatic generation above
   CALL mm_verify_scheme_consistency

END SUBROUTINE mm_init_scheme

!-------------------------------------------------------------------------------

SUBROUTINE mm_enable_stats

  USE mm_global_paras_mod
  mm_stats_printed = .TRUE.

END SUBROUTINE mm_enable_stats

!===============================================================================

SUBROUTINE GETMMBUFINFO(LUSEBUFMM,NBUFLEN,NBUFI,NBUFR,NBUFN)
   use cbifmm_module
   INTEGER NBUFLEN, NBUFI, NBUFR, NBUFN
   LOGICAL LUSEBUFMM
   LUSEBUFMM = USEBUFMM
   NBUFLEN   = MMBUFLEN
   NBUFI     = MAXBUFI
   NBUFR     = MAXBUFR
   NBUFN     = MAXBUFN
   RETURN
END

SUBROUTINE SETMMBUFINFO(LUSEBUFMM,NBUFLEN,NBUFI,NBUFR,NBUFN)
      use cbifmm_module
      INTEGER NBUFLEN, NBUFI, NBUFR, NBUFN 
      LOGICAL LUSEBUFMM
      USEBUFMM = LUSEBUFMM
      !      if (dopbc .and. pbcactive) then
      !         USEBUFMM = .FALSE.
      !!         write(lupri,*)'ATTENTION: set USEBUFMM = .FALSE. for pbc'
      !         write(*,*)'ATTENTION: set USEBUFMM = .FALSE. for pbc'
      !      endif
      MMBUFLEN = NBUFLEN
      MAXBUFI  = NBUFI
      MAXBUFR  = NBUFR
      MAXBUFN  = NBUFN
      RETURN
END
