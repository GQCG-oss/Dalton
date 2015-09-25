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

private
public :: determine_CABS_nbast, build_CABS_MO, build_RI_MO
CONTAINS
  subroutine determine_CABS_nbast(nbast_cabs,SETTING,lupri)
    implicit none
    integer,intent(inout) :: nbast_cabs
    integer,intent(in) :: lupri
    TYPE(LSSETTING),intent(inout) :: SETTING
    nbast_cabs = getNbasis(AOdfCABS,ContractedintType,SETTING%MOLECULE(1)%p,LUPRI)
  end subroutine determine_CABS_nbast



  subroutine build_CABS_MO(nbast_cabs,lupri,setting,CMO_std,CMO_cabs)
    implicit none
    !> Number of CABS AOs (CABS+)
    integer,intent(in) :: nbast_cabs
    !> Output file unit number
    integer,intent(in) :: lupri
    !> Integral settings
    TYPE(LSSETTING) :: SETTING
    !> Standard MOs
    type(matrix),intent(in) :: CMO_std
    !> CABS MOs - will be initialized inside subroutine!
    TYPE(MATRIX),intent(inout)    :: CMO_cabs
    !
    real(realk)     :: TIMSTR,TIMEND
    type(matrix) :: SAOCABS,S_cabs,tmp,S_minus_sqrt_cabs
    type(matrix) :: tmp_cabs,Vnull,tmp2,SMOCABS
    real(realk),pointer :: SV(:),optwrk(:),Sfull(:,:),S_cabsfull(:,:)
    integer,pointer :: IWORK(:)
    integer     :: lwork,nnull,luerr,IERR,INFO,I,nbastCABO,nov,nAO
    logical  :: ODSCREEN,Failed,doMPI

    CALL LSTIMER('START ',TIMSTR,TIMEND,lupri)

    IERR=0
    ODSCREEN = SETTING%SCHEME%OD_SCREEN
    SETTING%SCHEME%OD_SCREEN = .FALSE.
    luerr = 6
    nov = CMO_std%ncol  ! number of occupied + virtual MOs
    nAO = CMO_std%nrow  ! number of AO basis functions

    CALL mat_init(S_cabs,nbast_cabs,nbast_cabs)
    call mat_init(tmp_cabs,nbast_cabs,nbast_cabs)
    call mat_init(S_minus_sqrt_cabs,nbast_cabs,nbast_cabs)
    doMPI = SETTING%SCHEME%doMPI
    SETTING%SCHEME%doMPI = .FALSE.
    ! CABSAO,CABSAO overlap matrix

    call II_get_mixed_overlap(LUPRI,LUERR,SETTING,S_cabs,AOdfCABS,AOdfCABS,.FALSE.,.FALSE.)
    ! CABS,CABS overlap matrix in orthogonalized basis (S^-1/2)

    call lowdin_diag(nbast_cabs, S_cabs%elms,tmp_cabs%elms, S_minus_sqrt_cabs%elms, lupri)

    call mat_free(tmp_cabs)
    CALL mat_free(S_cabs)
    CALL mat_init(SAOCABS,nAO,nbast_cabs)

    ! AO,CABSAO overlap matrix
    call II_get_mixed_overlap(LUPRI,LUERR,setting,SAOCABS,AORegular,AOdfCABS,.FALSE.,.FALSE.)
    SETTING%SCHEME%doMPI=doMPI

    ! Calculate (MO, ortogonalized CABSAO) overlap: C^T S(AO,CABS) SCABS^(-1/2)
    call mat_init(tmp,nov,nbast_cabs)
    call mat_mul(CMO_std,SAOCABS,'T','N',1.0E0_realk,0.0E0_realk,tmp)
    call mat_init(SMOCABS,nov,nbast_cabs)
    call mat_mul(tmp,S_minus_sqrt_cabs,'N','N',1.0E0_realk,0.0E0_realk,SMOCABS)
    call mat_free(tmp)
    call mat_free(SAOCABS)

    call mat_init(tmp,SMOCABS%nrow,SMOCABS%nrow)
    call mat_init(tmp_cabs,SMOCABS%ncol,SMOCABS%ncol)
    call mat_zero(tmp)
    call mat_zero(tmp_cabs)

    call mem_alloc(SV,MIN(SMOCABS%nrow,SMOCABS%ncol))
    SV=0.0E0_realk
    call mem_alloc(optwrk,5) 
    call dgesvd('A','A',SMOCABS%nrow,SMOCABS%ncol,SMOCABS%elms,SMOCABS%nrow, &
         & SV,tmp%elms,tmp%nrow,tmp_cabs%elms,tmp_cabs%nrow,optwrk,-1,INFO)
    lwork = INT(optwrk(1))
    call mem_dealloc(optwrk) 
    call mem_alloc(optwrk,lwork) 
    call dgesvd('A','A',SMOCABS%nrow,SMOCABS%ncol,SMOCABS%elms,SMOCABS%nrow, &
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

    IF(nnull.NE.nbast_cabs-MIN(SMOCABS%nrow,SMOCABS%ncol))THEN       
       print*,'nnull',nnull
       print*,'nbast_cabs-MIN(SMOCABS%nrow,SMOCABS%ncol)',nbast_cabs-MIN(SMOCABS%nrow,SMOCABS%ncol)
       print*,'nbast_cabs,MIN(SMOCABS%nrow,SMOCABS%ncol)',nbast_cabs,MIN(SMOCABS%nrow,SMOCABS%ncol)
       CALL LSQUIT('error in build_CABS_MO',-1)
    ENDIF
    call mem_dealloc(SV)
    call mat_free(tmp)
    call mat_free(SMOCABS)

    !Construct CABS MO from V (tmp_cabs)
    call mat_init(tmp,nnull,nbast_cabs)
    call mat_retrieve_block(tmp_cabs,tmp%elms,nnull,nbast_cabs,nov+1,1)
    call mat_free(tmp_cabs)

    call mat_init(Vnull,nbast_cabs,nnull)

    call mat_trans(tmp,Vnull)
    call mat_free(tmp)

    call mat_init(CMO_cabs,nbast_cabs,nnull)
    call mat_mul(S_minus_sqrt_cabs,Vnull,'N','N',1.0E0_realk,0.0E0_realk,CMO_cabs)

    !test of cabs orthonomality
    !If(DECinfo%F12debug) then
    !call test_CABS_MO_orthonomality(CMO_cabs,SETTING,lupri)
    !endif

    call mat_free(S_minus_sqrt_cabs)
    call mat_free(Vnull)       

    CALL LSTIMER('build_CABS_MO',TIMSTR,TIMEND,lupri)       
    SETTING%SCHEME%OD_SCREEN = ODSCREEN

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
    !> RI MOs, matrix is initialized here
    TYPE(MATRIX),intent(inout)    :: CMO_RI
    !
    real(realk)     :: TIMSTR,TIMEND
    type(matrix) :: S,Smix,S_cabs,tmp,S_minus_sqrt,S_minus_sqrt_cabs
    type(matrix) :: tmp_cabs,tmp2
    real(realk),pointer :: SV(:),optwrk(:)
    integer     :: lwork,nbast,nnull,luerr,IERR,INFO,I
    logical     :: doMPI

    luerr = 6
    CALL LSTIMER('START ',TIMSTR,TIMEND,lupri)
    CALL mat_init(S_cabs,nbast_cabs,nbast_cabs)
    call mat_init(tmp_cabs,nbast_cabs,nbast_cabs)
    doMPI=SETTING%SCHEME%doMPI
    SETTING%SCHEME%doMPI=.FALSE.
    call II_get_mixed_overlap(LUPRI,LUERR,SETTING,S_cabs,AOdfCABS,AOdfCABS,.FALSE.,.FALSE.)
    SETTING%SCHEME%doMPI=doMPI

    call mat_init(CMO_RI,nbast_cabs,nbast_cabs)
    call lowdin_diag(nbast_cabs, S_cabs%elms,tmp_cabs%elms, CMO_RI%elms, lupri)
    CALL mat_free(S_cabs)
    call mat_free(tmp_cabs)

    CALL LSTIMER('build_RI_MO',TIMSTR,TIMEND,lupri)

  end subroutine build_RI_MO

  subroutine test_CABS_MO_orthonomality(CMO_cabs,SETTING,lupri)
    implicit none
    integer :: lupri
    TYPE(LSSETTING) :: SETTING
    TYPE(MATRIX)    :: CMO_cabs
!
    TYPE(MATRIX)    :: tmp,tmp2,tmp3,S_cabs
    integer ::  nbast_cabs,luerr
    logical     :: doMPI
    luerr=6
    nbast_cabs = CMO_cabs%nrow
    CALL mat_init(S_cabs,nbast_cabs,nbast_cabs)
    doMPI=SETTING%SCHEME%doMPI
    SETTING%SCHEME%doMPI=.FALSE.
    call II_get_mixed_overlap(LUPRI,LUERR,SETTING,S_cabs,AOdfCABS,AOdfCABS,.FALSE.,.FALSE.)
    SETTING%SCHEME%doMPI=doMPI

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
    logical     :: doMPI
    luerr=6
    nbast_cabs = CMO_cabs%nrow
    nbast = CMO%nrow
    CALL mat_init(Smix,nbast,nbast_cabs)    
    doMPI=SETTING%SCHEME%doMPI
    SETTING%SCHEME%doMPI=.FALSE.
    call II_get_mixed_overlap(LUPRI,LUERR,SETTING,Smix,AORegular,AOdfCABS,.FALSE.,.FALSE.)
    SETTING%SCHEME%doMPI=doMPI
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
