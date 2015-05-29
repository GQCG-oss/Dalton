MODULE CABS_operations
  use memory_handling!,only: mem_alloc, mem_dealloc
  use lsparameters
  use TYPEDEF
  use Matrix_module
  use lowdin_module
  use matrix_operations
  use integralinterfaceMod
  use lstiming
  use dec_typedef_module

SAVE
logical  :: CMO_CABS_save_created
TYPE(Matrix) :: CMO_CABS_save
logical  :: CMO_RI_save_created
TYPE(Matrix) :: CMO_RI_save
logical  :: Save_activated_cabs
logical  :: Save_activated_ri
CONTAINS
  subroutine determine_CABS_nbast(nbast_cabs,nnull,SETTING,lupri)
    implicit none
    integer,intent(inout) :: nbast_cabs,nnull
    integer,intent(in) :: lupri
    TYPE(LSSETTING),intent(inout) :: SETTING
    integer :: nbast
    nbast = getNbasis(AORegular,ContractedintType,SETTING%MOLECULE(1)%p,LUPRI)
    nbast_cabs = getNbasis(AOdfCABS,ContractedintType,SETTING%MOLECULE(1)%p,LUPRI)
    nnull = nbast_cabs-MIN(nbast_cabs,nbast)
  end subroutine determine_CABS_nbast

  subroutine init_cabs(Fulldriver)
    implicit none
    logical,intent(in),optional :: Fulldriver
    CMO_CABS_save_created = .FALSE.
    IF(present(Fulldriver))THEN
       Save_activated_cabs = Fulldriver
    ELSE
       Save_activated_cabs = .FALSE.
    ENDIF
  end subroutine init_cabs

  subroutine init_ri(Fulldriver)
    implicit none
    logical,intent(in),optional :: Fulldriver
    CMO_RI_save_created = .FALSE.
    IF(present(Fulldriver))THEN
       Save_activated_ri = Fulldriver
    ELSE
       Save_activated_ri = .FALSE.
    ENDIF
  end subroutine init_ri

  subroutine free_cabs()
    implicit none
        ! print*,'FREE CABS'
    IF(CMO_CABS_save_created)THEN
       call mat_free(CMO_CABS_save)
    ENDIF
    CMO_CABS_save_created = .FALSE.
    IF(CMO_RI_save_created)THEN
       call mat_free(CMO_RI_save)
    ENDIF
    CMO_RI_save_created = .FALSE.
    Save_activated_ri = .FALSE.
    Save_activated_cabs = .FALSE.
  end subroutine free_cabs

  subroutine build_CABS_MO(CMO_cabs,nbast_cabs,SETTING,lupri)
    implicit none
    integer :: lupri,nbast_cabs
    TYPE(LSSETTING) :: SETTING
    TYPE(MATRIX)    :: CMO_cabs
!
    real(realk)     :: TIMSTR,TIMEND
    type(matrix) :: S,Smix,S_cabs,tmp,S_minus_sqrt,S_minus_sqrt_cabs
    type(matrix) :: tmp_cabs,Vnull,tmp2
    real(realk),pointer :: SV(:),optwrk(:),Sfull(:,:),S_cabsfull(:,:)
    integer,pointer :: IWORK(:)
    integer     :: lwork,nbast,nnull,luerr,IERR,INFO,I,nbastCABO
    logical  :: ODSCREEN,Failed,VerifyCabs

    VerifyCabs = DECinfo%F12debug
#ifdef VAR_LSDEBUG
    VerifyCabs = .TRUE.
#endif
    IERR=0
    IF(CMO_CABS_save_created)THEN
       call mat_assign(CMO_cabs,CMO_CABS_save)
    ELSE
       ODSCREEN = SETTING%SCHEME%OD_SCREEN
       SETTING%SCHEME%OD_SCREEN = .FALSE.
       luerr = 6
       CALL LSTIMER('START ',TIMSTR,TIMEND,lupri)
       nbast = getNbasis(AORegular,ContractedintType,SETTING%MOLECULE(1)%p,LUPRI)
       CALL mat_init(S,nbast,nbast)
       call mat_init(S_minus_sqrt,nbast,nbast)
       
       call mat_init(tmp,nbast,nbast)
       call II_get_mixed_overlap(LUPRI,LUERR,SETTING,S,AORegular,AORegular,.FALSE.,.FALSE.)
       call lowdin_diag(nbast, S%elms,tmp%elms, S_minus_sqrt%elms, lupri)
       call mat_free(tmp)
       call mat_free(S)
       
       CALL mat_init(S_cabs,nbast_cabs,nbast_cabs)
       call mat_init(S_minus_sqrt_cabs,nbast_cabs,nbast_cabs)
       
       call mat_init(tmp_cabs,nbast_cabs,nbast_cabs)
       call II_get_mixed_overlap(LUPRI,LUERR,SETTING,S_cabs,AOdfCABS,AOdfCABS,.FALSE.,.FALSE.)
       If(VerifyCabs) then
          !Test that S_cabs have the structure 
          !  S_cabs =  (  S   |     *      )
          !            (  *   |  cabs_only )
          CALL mat_init(S,nbast,nbast)
          call II_get_mixed_overlap(LUPRI,LUERR,SETTING,S,AORegular,AORegular,.FALSE.,.FALSE.)
          call mem_alloc(Sfull,nbast,nbast)
          CALL mat_to_full(S,1.0E0_realk,Sfull)
          CALL mat_free(S)

          call mem_alloc(S_cabsfull,nbast_cabs,nbast_cabs)
          CALL mat_to_full(S_cabs,1.0E0_realk,S_cabsfull)
          call verifyCABS_HaveRegularFirst(Sfull,nbast,S_cabsfull,nbast_cabs,Failed)
          call mem_dealloc(Sfull)
          IF(Failed)THEN
             CALL LSQUIT('Something Wrong CABSAO do not have regular AO first',-1)
          ELSE
             !check CABSonly is last 
             nbastCABO = nbast_cabs - nbast
             CALL mat_init(S,nbastCABO,nbastCABO)
             call II_get_mixed_overlap(LUPRI,LUERR,SETTING,S,AOdfCABO,AOdfCABO,.FALSE.,.FALSE.)
             call mem_alloc(Sfull,nbastCABO,nbastCABO)
             CALL mat_to_full(S,1.0E0_realk,Sfull)
             CALL mat_free(S)
             call verifyCABS_HaveCABOLast(Sfull,nbastCABO,nbast,S_cabsfull,nbast_cabs,Failed)
             call mem_dealloc(Sfull)
             IF(Failed)THEN
                CALL LSQUIT('Something Wrong CABSAO do not have CABS only last',-1)
             ELSE
                !check last block
                CALL mat_init(S,nbastCABO,nbast)
                call II_get_mixed_overlap(LUPRI,LUERR,SETTING,S,AOdfCABO,AORegular,.FALSE.,.FALSE.)
                call mem_alloc(Sfull,nbastCABO,nbast)
                CALL mat_to_full(S,1.0E0_realk,Sfull)
                CALL mat_free(S)
                call verifyCABS3(Sfull,nbastCABO,nbast,S_cabsfull,nbast_cabs,Failed)
                call mem_dealloc(Sfull)
                IF(Failed)THEN
                   CALL LSQUIT('Something Wrong CABS3',-1)
                ENDIF
             ENDIF
          ENDIF
          call mem_dealloc(S_cabsfull)
          WRITE(lupri,*)'CABS basis test is successful'
       endif

       call lowdin_diag(nbast_cabs, S_cabs%elms,tmp_cabs%elms, S_minus_sqrt_cabs%elms, lupri)
       call mat_free(tmp_cabs)
       CALL mat_free(S_cabs)       
       CALL mat_init(Smix,nbast,nbast_cabs)
       call II_get_mixed_overlap(LUPRI,LUERR,setting,Smix,AORegular,AOdfCABS,.FALSE.,.FALSE.)
       
       call mat_init(tmp,nbast,nbast_cabs)
       call mat_mul(S_minus_sqrt,Smix,'N','N',1.0E0_realk,0.0E0_realk,tmp)
       call mat_mul(tmp,S_minus_sqrt_cabs,'N','N',1.0E0_realk,0.0E0_realk,Smix)

       call mat_free(tmp)
       
       call mat_free(S_minus_sqrt)
       
       call mat_init(tmp,Smix%nrow,Smix%nrow)
       call mat_init(tmp_cabs,Smix%ncol,Smix%ncol)

       call mat_zero(tmp)
       call mat_zero(tmp_cabs)

       call mem_alloc(SV,MIN(Smix%nrow,Smix%ncol))
       SV=0.0E0_realk
       call mem_alloc(optwrk,5) 
       call dgesvd('A','A',Smix%nrow,Smix%ncol,Smix%elms,Smix%nrow, &
            & SV,tmp%elms,tmp%nrow,tmp_cabs%elms,tmp_cabs%nrow,optwrk,-1,INFO)
       lwork = INT(optwrk(1))
       call mem_dealloc(optwrk) 
       call mem_alloc(optwrk,lwork) 
       call dgesvd('A','A',Smix%nrow,Smix%ncol,Smix%elms,Smix%nrow, &
            & SV,tmp%elms,tmp%nrow,tmp_cabs%elms,tmp_cabs%nrow,optwrk,lwork,INFO)
       IF( INFO.GT.0 ) THEN
          WRITE(*,*)'The algorithm computing SVD failed to converge.'
          call lsquit('The algorithm computing SVD failed to converge.',-1)
       ENDIF
       call mem_dealloc(optwrk) 
       nnull = nbast_cabs
       DO I=1,SIZE(SV)
          IF(ABS(SV(I)).GT.1.0E-12_realk)THEN
             nnull=nnull-1
          ENDIF
       ENDDO
       IF(nnull.NE.nbast_cabs-MIN(Smix%nrow,Smix%ncol))THEN       
          print*,'nnull',nnull
          print*,'nbast_cabs-MIN(Smix%nrow,Smix%ncol)',nbast_cabs-MIN(Smix%nrow,Smix%ncol)
          print*,'nbast_cabs,MIN(Smix%nrow,Smix%ncol)',nbast_cabs,MIN(Smix%nrow,Smix%ncol)
          CALL LSQUIT('error in build_CABS_MO',-1)
       ENDIF
       !    nnull = nbast_cabs-MIN(Smix%nrow,Smix%ncol)
       call mem_dealloc(SV)
       call mat_free(tmp)
       call mat_free(Smix)
       
       !Construct CABS MO from V (tmp_cabs)
       call mat_init(tmp,nnull,nbast_cabs)
       call mat_retrieve_block(tmp_cabs,tmp%elms,nnull,nbast_cabs,nbast+1,1)
       call mat_free(tmp_cabs)
       
       call mat_init(Vnull,nbast_cabs,nnull)
       
       call mat_trans(tmp,Vnull)
       call mat_free(tmp)
       
       call mat_mul(S_minus_sqrt_cabs,Vnull,'N','N',1.0E0_realk,0.0E0_realk,CMO_cabs)

       !test of cabs orthonomality
       !If(DECinfo%F12debug) then
        !call test_CABS_MO_orthonomality(CMO_cabs,SETTING,lupri)
       !endif

       call mat_free(S_minus_sqrt_cabs)
       call mat_free(Vnull)       
       
       IF(Save_activated_cabs)THEN
          CMO_CABS_save_created = .TRUE.
          CALL MAT_INIT(CMO_CABS_save,CMO_cabs%nrow,CMO_cabs%ncol)
          call mat_assign(CMO_CABS_save,CMO_cabs)
       ENDIF

       !   print*,'BUILD CABS'

       CALL LSTIMER('build_CABS_MO',TIMSTR,TIMEND,lupri)       
       SETTING%SCHEME%OD_SCREEN = ODSCREEN
    ENDIF
  end subroutine build_CABS_MO

  subroutine verifyCABS_HaveRegularFirst(Sfull,nbast,S_cabsfull,nbast_cabs,Failed)
    implicit none
    integer,intent(in) :: nbast,nbast_cabs
    real(realk),intent(in) :: Sfull(nbast,nbast),S_cabsfull(nbast_cabs,nbast_cabs)
    logical,intent(inout) :: Failed 
    !
    integer :: I,J
    Failed =.FALSE.
    VerifyCabsJloop: DO J=1,nbast
       DO I=1,nbast
          IF(ABS(Sfull(I,J)-S_cabsfull(I,J)).GT.1.0E-12_realk)THEN
             Failed =.TRUE.
             print*,'verifyAOFirst: I=',I,'J=',J
             print*,'verifyAOFirst: Sfull(I,J)',Sfull(I,J)
             print*,'verifyAOFirst: S_cabsfull(I,J)',S_cabsfull(I,J)
             EXIT VerifyCabsJloop
          ENDIF
       ENDDO
    ENDDO VerifyCabsJloop
  end subroutine verifyCABS_HaveRegularFirst

  subroutine verifyCABS_HaveCABOLast(Sfull,nbastCABO,nbast,S_cabsfull,nbast_cabs,Failed)
    implicit none
    integer,intent(in) :: nbast,nbast_cabs,nbastCABO
    real(realk),intent(in) :: Sfull(nbastCABO,nbastCABO),S_cabsfull(nbast_cabs,nbast_cabs)
    logical,intent(inout) :: Failed 
    !
    integer :: I,J
    Failed =.FALSE.
    VerifyCabsJloop2: DO J=1,nbastCABO
       DO I=1,nbastCABO
          IF(ABS(Sfull(I,J)-S_cabsfull(nbast+I,nbast+J)).GT.1.0E-12_realk)THEN
             Failed =.TRUE.
             print*,'verifyCABSLast: I=',I,'J=',J,'nbast=',nbast
             print*,'verifyCABSLast: Sfull(I,J)',Sfull(I,J)
             print*,'verifyCABSLast: S_cabsfull(nbast+I,nbast+J)',S_cabsfull(nbast+I,nbast+J)
             EXIT VerifyCabsJloop2
          ENDIF
       ENDDO
    ENDDO VerifyCabsJloop2
  end subroutine verifyCABS_HaveCABOLast

  subroutine verifyCABS3(Sfull,nbastCABO,nbast,S_cabsfull,nbast_cabs,Failed)
    implicit none
    integer,intent(in) :: nbast,nbast_cabs,nbastCABO
    real(realk),intent(in) :: Sfull(nbastCABO,nbast),S_cabsfull(nbast_cabs,nbast_cabs)
    logical,intent(inout) :: Failed 
    !
    integer :: I,J
    Failed =.FALSE.
    VerifyCabsJloop3: DO J=1,nbast
       DO I=1,nbastCABO
          IF(ABS(Sfull(I,J)-S_cabsfull(nbast+I,J)).GT.1.0E-12_realk)THEN
             Failed =.TRUE.
             print*,'verifyCABS3: I=',I,'J=',J,'nbast=',nbast
             print*,'verifyCABS3: Sfull(I,J)',Sfull(I,J)
             print*,'verifyCABS3: S_cabsfull(nbast+I,nbast+J)',S_cabsfull(nbast+I,nbast+J)
             EXIT VerifyCabsJloop3
          ENDIF
       ENDDO
    ENDDO VerifyCabsJloop3
  end subroutine verifyCABS3


  subroutine build_RI_MO(CMO_RI,nbast_cabs,SETTING,lupri)
    implicit none
    integer :: lupri,nbast_cabs
    TYPE(LSSETTING) :: SETTING
    TYPE(MATRIX)    :: CMO_RI
!
    real(realk)     :: TIMSTR,TIMEND
    type(matrix) :: S,Smix,S_cabs,tmp,S_minus_sqrt,S_minus_sqrt_cabs
    type(matrix) :: tmp_cabs,tmp2
    real(realk),pointer :: SV(:),optwrk(:)
    integer     :: lwork,nbast,nnull,luerr,IERR,INFO,I
    IF(CMO_RI_save_created)THEN
       call mat_assign(CMO_RI,CMO_RI_save)
    ELSE
       luerr = 6
       CALL LSTIMER('START ',TIMSTR,TIMEND,lupri)
       CALL mat_init(S_cabs,nbast_cabs,nbast_cabs)
       call mat_init(tmp_cabs,nbast_cabs,nbast_cabs)
       call II_get_mixed_overlap(LUPRI,LUERR,SETTING,S_cabs,AOdfCABS,AOdfCABS,.FALSE.,.FALSE.)
       call lowdin_diag(nbast_cabs, S_cabs%elms,tmp_cabs%elms, CMO_RI%elms, lupri)
       CALL mat_free(S_cabs)
       call mat_free(tmp_cabs)

       IF(Save_activated_ri)THEN
          CMO_RI_save_created = .TRUE.
          CALL MAT_INIT(CMO_RI_save,CMO_RI%nrow,CMO_RI%ncol)
          call mat_assign(CMO_RI_save,CMO_RI)
       ENDIF
       CALL LSTIMER('build_RI_MO',TIMSTR,TIMEND,lupri)
    ENDIF
  end subroutine build_RI_MO

  subroutine test_CABS_MO_orthonomality(CMO_cabs,SETTING,lupri)
    implicit none
    integer :: lupri
    TYPE(LSSETTING) :: SETTING
    TYPE(MATRIX)    :: CMO_cabs
!
    TYPE(MATRIX)    :: tmp,tmp2,tmp3,S_cabs
    integer ::  nbast_cabs,luerr
    luerr=6
    nbast_cabs = CMO_cabs%nrow
    CALL mat_init(S_cabs,nbast_cabs,nbast_cabs)
    call II_get_mixed_overlap(LUPRI,LUERR,SETTING,S_cabs,AOdfCABS,AOdfCABS,.FALSE.,.FALSE.)

    call mat_init(tmp2,Cmo_cabs%ncol,nbast_cabs)
    call mat_init(tmp,Cmo_cabs%ncol,Cmo_cabs%ncol)
    call mat_mul(Cmo_cabs,S_cabs,'T','N',1.0E0_realk,0.0E0_realk,tmp2)
    call mat_mul(tmp2,Cmo_cabs,'N','N',1.0E0_realk,0.0E0_realk,tmp)  
    
    CALL mat_free(S_cabs)
    call mat_free(tmp2)
    call mat_init(tmp2,tmp%nrow,tmp%ncol)
    call mat_init(tmp3,tmp%nrow,tmp%ncol)
    call mat_identity(tmp2)
    call mat_add(1E0_realk,tmp,-1E0_realk,tmp2,tmp3)
    
    IF(sqrt(mat_sqnorm2(tmp3)/tmp%nrow).GT.1.0E-10_realk)THEN
       write(lupri,*)'sqrt(Ccabs^T*Scabs*Ccabs - I)',sqrt(mat_sqnorm2(tmp3)/tmp%nrow)  
       call mat_print(tmp,1,tmp%nrow,1,tmp%ncol,lupri)
       call lsquit('CABS not Orthonormal',-1)
    ELSE
       write(lupri,*)'sqrt(Ccabs^T*Scabs*Ccabs - I)',sqrt(mat_sqnorm2(tmp3)/tmp%nrow)  
    ENDIF
    call mat_free(tmp)
    call mat_free(tmp2)
    call mat_free(tmp3)

  end subroutine test_CABS_MO_orthonomality

  subroutine test_CABS_MO_orthogonality(CMO,CMO_cabs,setting,lupri)
    implicit none
    TYPE(MATRIX)    :: CMO_cabs,CMO
    TYPE(LSSETTING) :: SETTING
    integer :: lupri
!
    TYPE(MATRIX)    :: tmp,tmp2,Smix
    integer ::  nbast_cabs,nbast,luerr
    luerr=6
    nbast_cabs = CMO_cabs%nrow
    nbast = CMO%nrow
    CALL mat_init(Smix,nbast,nbast_cabs)    
    call II_get_mixed_overlap(LUPRI,LUERR,SETTING,Smix,AORegular,AOdfCABS,.FALSE.,.FALSE.)
    call mat_init (tmp2, nbast, nbast_cabs)
    call mat_init (tmp, nbast, Cmo_cabs%ncol)
    call mat_mul(Cmo,Smix,'T','N',1.0E0_realk,0.0E0_realk,tmp2)
    call mat_mul(tmp2,Cmo_cabs,'N','N',1.0E0_realk,0.0E0_realk,tmp)  
    IF(sqrt(mat_sqnorm2(tmp)/tmp%nrow).GT.1.0E-10_realk)THEN
       write(lupri,*)'Ccabs^T*Scabs*Ccabs = '  
       call mat_print(tmp,1,tmp%nrow,1,tmp%ncol,lupri)
       call lsquit('CABS not Orthogonal to MOs',-1)
    ENDIF
    call mat_free(tmp)
    call mat_free(tmp2)
    call mat_free(Smix)
  end subroutine test_CABS_MO_orthogonality

end MODULE CABS_operations
