!> @file 
!> Contains lstensor structure and associated subroutines
MODULE LSTENSOR_OPERATIONSMOD
  use LSTENSOR_TYPETYPE
  use precision
  use basis_type
  use Matrix_module
  use matrix_operations_csr, only: zeroCSR,mat_csr_allocate
  use matrix_operations_scalapack
  use matrix_operations
  use memory_handling
  use OD_Type
  use AO_Type
  use OD_TypeType
  use AO_TypeType
#ifdef VAR_MPI
  use infpar_module
#endif
  use LSTENSORmem
  INTERFACE Build_mat_from_lst
     MODULE PROCEDURE Build_singlemat_from_lst, &
          &           Build_matarray_from_lst
  END INTERFACE Build_mat_from_lst
  TYPE LSTAUXITEM
     TYPE(MATRIX),pointer :: sMAT
     TYPE(MATRIXP),pointer :: aMAT(:)
     INTEGER,pointer :: COL2(:,:)
     INTEGER,pointer :: ROW2(:,:)
     INTEGER,pointer :: NNZ2(:)
     INTEGER,pointer :: COL1(:)
     INTEGER,pointer :: ROW1(:)
     INTEGER,pointer :: NNZ1
     REAL(REALK),pointer :: VAL1(:)
     REAL(REALK),pointer :: VAL2(:,:)
     REAL(REALK),pointer :: MAT1dim(:)
     REAL(REALK),pointer :: MAT2dim(:,:)
     REAL(REALK),pointer :: MAT3dim(:,:,:)
     REAL(REALK),pointer :: MAT4dim(:,:,:,:)
     REAL(REALK),pointer :: MAT5dim(:,:,:,:,:)
  END TYPE LSTAUXITEM

CONTAINS
  !PRIMITIVE SUBROUTINES

  Subroutine SLSAOTENSOR_nullify(SLSAOTENSORitem)
    implicit none
    TYPE(SLSAOTENSOR) :: SLSAOTENSORitem
    !
    integer :: I
    SLSAOTENSORitem%nelms=0
    NULLIFY(SLSAOTENSORitem%selms)
    NULLIFY(SLSAOTENSORitem%nOrb)
    NULLIFY(SLSAOTENSORitem%startLocalOrb)
    DO I=1,size(SLSAOTENSORitem%nLocal,1)
       SLSAOTENSORitem%nLocal(I) =0
    ENDDO
    DO I=1,size(SLSAOTENSORitem%ATOM,1)
       SLSAOTENSORitem%ATOM(I) =0
    ENDDO
    DO I=1,size(SLSAOTENSORitem%AOBATCH,1)
       SLSAOTENSORitem%AOBATCH(I) =0
    ENDDO
  End Subroutine SLSAOTENSOR_nullify

  Subroutine SLSAOTENSOR_copy(oldSLSAOTENSORitem,newSLSAOTENSORitem,lupri)
    implicit none
    TYPE(SLSAOTENSOR),intent(in) :: oldSLSAOTENSORitem
    TYPE(SLSAOTENSOR),intent(inout) :: newSLSAOTENSORitem
    integer,intent(in) :: lupri
    !
    integer :: I,n1
    newSLSAOTENSORitem%nelms= oldSLSAOTENSORitem%nelms
    newSLSAOTENSORitem%maxBat= oldSLSAOTENSORitem%maxBat
    IF(ASSOCIATED(oldSLSAOTENSORitem%selms))THEN
       n1 = size(oldSLSAOTENSORitem%selms,1)
       call mem_lstpointer_alloc(newSLSAOTENSORitem%selms,n1)
       DO I=1,n1
          newSLSAOTENSORitem%selms(I) = oldSLSAOTENSORitem%selms(I)
       ENDDO
    ELSE
       NULLIFY(newSLSAOTENSORitem%selms)
    ENDIF
    IF(ASSOCIATED(oldSLSAOTENSORitem%nOrb))THEN
       n1 = size(oldSLSAOTENSORitem%nOrb,1)
       call mem_lstpointer_alloc(newSLSAOTENSORitem%nOrb,n1)
       DO I=1,n1
          newSLSAOTENSORitem%nOrb(I) = oldSLSAOTENSORitem%nOrb(I)
       ENDDO
    ELSE
       NULLIFY(newSLSAOTENSORitem%nOrb)
    ENDIF
    IF(ASSOCIATED(oldSLSAOTENSORitem%startLocalOrb))THEN
       n1 = size(oldSLSAOTENSORitem%startLocalOrb,1)
       call mem_lstpointer_alloc(newSLSAOTENSORitem%startLocalOrb,n1)
       DO I=1,n1
          newSLSAOTENSORitem%startLocalOrb(I) = oldSLSAOTENSORitem%startLocalOrb(I)
       ENDDO
    ELSE
       NULLIFY(newSLSAOTENSORitem%startLocalOrb)
    ENDIF
    n1 = size(oldSLSAOTENSORitem%nLocal,1)
    DO I=1,n1
       newSLSAOTENSORitem%nLocal(I) = oldSLSAOTENSORitem%nLocal(I)
    ENDDO
    n1 = size(oldSLSAOTENSORitem%ATOM,1)
    DO I=1,n1
       newSLSAOTENSORitem%ATOM(I) = oldSLSAOTENSORitem%ATOM(I)
    ENDDO
    n1 = size(oldSLSAOTENSORitem%AOBATCH,1)
    DO I=1,n1
       newSLSAOTENSORitem%AOBATCH(I) = oldSLSAOTENSORitem%AOBATCH(I)
    ENDDO
  End Subroutine SLSAOTENSOR_copy

  Subroutine SLSAOTENSOR_p_copy(oldSLSAOTENSORitem,newSLSAOTENSORitem,lupri)
    implicit none
    TYPE(SLSAOTENSOR),intent(in) :: oldSLSAOTENSORitem
    TYPE(SLSAOTENSOR),intent(inout) :: newSLSAOTENSORitem
    integer,intent(in) :: lupri
    !
    integer :: I,n1
    newSLSAOTENSORitem%nelms= oldSLSAOTENSORitem%nelms
    newSLSAOTENSORitem%maxBat= oldSLSAOTENSORitem%maxBat
    IF(ASSOCIATED(oldSLSAOTENSORitem%selms))THEN
       n1 = size(oldSLSAOTENSORitem%selms,1)
       call mem_lstpointer_alloc(newSLSAOTENSORitem%selms,n1)
    ELSE
       NULLIFY(newSLSAOTENSORitem%selms)
    ENDIF
    IF(ASSOCIATED(oldSLSAOTENSORitem%nOrb))THEN
       n1 = size(oldSLSAOTENSORitem%nOrb,1)
       call mem_lstpointer_alloc(newSLSAOTENSORitem%nOrb,n1)
    ELSE
       NULLIFY(newSLSAOTENSORitem%nOrb)
    ENDIF
    IF(ASSOCIATED(oldSLSAOTENSORitem%startLocalOrb))THEN
       n1 = size(oldSLSAOTENSORitem%startLocalOrb,1)
       call mem_lstpointer_alloc(newSLSAOTENSORitem%startLocalOrb,n1)
    ELSE
       NULLIFY(newSLSAOTENSORitem%startLocalOrb)
    ENDIF
    n1 = size(oldSLSAOTENSORitem%nLocal,1)
    DO I=1,n1
       newSLSAOTENSORitem%nLocal(I) = oldSLSAOTENSORitem%nLocal(I)
    ENDDO
    n1 = size(oldSLSAOTENSORitem%ATOM,1)
    DO I=1,n1
       newSLSAOTENSORitem%ATOM(I) = oldSLSAOTENSORitem%ATOM(I)
    ENDDO
    n1 = size(oldSLSAOTENSORitem%AOBATCH,1)
    DO I=1,n1
       newSLSAOTENSORitem%AOBATCH(I) = oldSLSAOTENSORitem%AOBATCH(I)
    ENDDO
  End Subroutine SLSAOTENSOR_p_copy

  ! The Primitive subroutines for TYPE:LSAOTENSOR
  Subroutine LSAOTENSOR_print(LSAOTENSORitem,lupri)
    implicit none
    TYPE(LSAOTENSOR),intent(in) :: LSAOTENSORitem
    integer,intent(in) :: lupri
    !
    integer :: I,J,K,L,M,n1,n2,n3

  End Subroutine LSAOTENSOR_print

  Subroutine LSAOTENSOR_nullify(LSAOTENSORitem)
    implicit none
    TYPE(LSAOTENSOR) :: LSAOTENSORitem
    !
    integer :: I,J,K,L,M
    LSAOTENSORitem%nelms=0
    NULLIFY(LSAOTENSORitem%elms)
    NULLIFY(LSAOTENSORitem%nOrb)
    NULLIFY(LSAOTENSORitem%startLocalOrb)
    NULLIFY(LSAOTENSORitem%startGlobalOrb)
    NULLIFY(LSAOTENSORitem%nAngmom)
    DO I=1,size(LSAOTENSORitem%nLocal,1)
       LSAOTENSORitem%nLocal(I) =0
    ENDDO
    DO I=1,size(LSAOTENSORitem%ATOM,1)
       LSAOTENSORitem%ATOM(I) =0
    ENDDO
    DO I=1,size(LSAOTENSORitem%AOBATCH,1)
       LSAOTENSORitem%AOBATCH(I) =0
    ENDDO
  End Subroutine LSAOTENSOR_nullify

  Subroutine LSAOTENSOR_copy(oldLSAOTENSORitem,newLSAOTENSORitem,lupri)
    implicit none
    TYPE(LSAOTENSOR),intent(in) :: oldLSAOTENSORitem
    TYPE(LSAOTENSOR),intent(inout) :: newLSAOTENSORitem
    integer,intent(in) :: lupri
    !
    integer :: I,n1
    newLSAOTENSORitem%nelms= oldLSAOTENSORitem%nelms
    newLSAOTENSORitem%maxBat= oldLSAOTENSORitem%maxBat
    newLSAOTENSORitem%maxAng= oldLSAOTENSORitem%maxAng
    IF(ASSOCIATED(oldLSAOTENSORitem%elms))THEN
       n1 = size(oldLSAOTENSORitem%elms)
       call mem_lstpointer_alloc(newLSAOTENSORitem%elms,n1)
       DO I=1,n1
          newLSAOTENSORitem%elms(I) = oldLSAOTENSORitem%elms(I)
       ENDDO
    ELSE
       NULLIFY(newLSAOTENSORitem%elms)
    ENDIF
    IF(ASSOCIATED(oldLSAOTENSORitem%nOrb))THEN
       n1 = size(oldLSAOTENSORitem%nOrb)
       call mem_lstpointer_alloc(newLSAOTENSORitem%nOrb,n1)
       DO I=1,n1
          newLSAOTENSORitem%nOrb(I) = oldLSAOTENSORitem%nOrb(I)
       ENDDO
    ELSE
       NULLIFY(newLSAOTENSORitem%nOrb)
    ENDIF
    IF(ASSOCIATED(oldLSAOTENSORitem%startLocalOrb))THEN
       n1 = size(oldLSAOTENSORitem%startLocalOrb)
       call mem_lstpointer_alloc(newLSAOTENSORitem%startLocalOrb,n1)
       DO I=1,n1
          newLSAOTENSORitem%startLocalOrb(I) = oldLSAOTENSORitem%startLocalOrb(I)
       ENDDO
    ELSE
       NULLIFY(newLSAOTENSORitem%startLocalOrb)
    ENDIF
    IF(ASSOCIATED(oldLSAOTENSORitem%startGlobalOrb))THEN
       n1 = size(oldLSAOTENSORitem%startGlobalOrb)
       call mem_lstpointer_alloc(newLSAOTENSORitem%startGlobalOrb,n1)
       DO I=1,n1
          newLSAOTENSORitem%startGlobalOrb(I) = oldLSAOTENSORitem%startGlobalOrb(I)
       ENDDO
    ELSE
       NULLIFY(newLSAOTENSORitem%startGlobalOrb)
    ENDIF
    IF(ASSOCIATED(oldLSAOTENSORitem%nAngmom))THEN
       n1 = size(oldLSAOTENSORitem%nAngmom)
       call mem_lstpointer_alloc(newLSAOTENSORitem%nAngmom,n1)
       DO I=1,n1
          newLSAOTENSORitem%nAngmom(I) = oldLSAOTENSORitem%nAngmom(I)
       ENDDO
    ELSE
       NULLIFY(newLSAOTENSORitem%nAngmom)
    ENDIF
    n1 = size(oldLSAOTENSORitem%nLocal)
    DO I=1,n1
       newLSAOTENSORitem%nLocal(I) = oldLSAOTENSORitem%nLocal(I)
    ENDDO
    n1 = size(oldLSAOTENSORitem%ATOM)
    DO I=1,n1
       newLSAOTENSORitem%ATOM(I) = oldLSAOTENSORitem%ATOM(I)
    ENDDO
    n1 = size(oldLSAOTENSORitem%AOBATCH)
    DO I=1,n1
       newLSAOTENSORitem%AOBATCH(I) = oldLSAOTENSORitem%AOBATCH(I)
    ENDDO
  End Subroutine LSAOTENSOR_copy

  Subroutine LSAOTENSOR_p_copy(oldLSAOTENSORitem,newLSAOTENSORitem,lupri)
    implicit none
    TYPE(LSAOTENSOR),intent(in) :: oldLSAOTENSORitem
    TYPE(LSAOTENSOR),intent(inout) :: newLSAOTENSORitem
    integer,intent(in) :: lupri
    !
    integer :: I,J,K,L,M,n1
    newLSAOTENSORitem%nelms= oldLSAOTENSORitem%nelms
    newLSAOTENSORitem%maxBat= oldLSAOTENSORitem%maxBat
    newLSAOTENSORitem%maxAng= oldLSAOTENSORitem%maxAng
    IF(ASSOCIATED(oldLSAOTENSORitem%elms))THEN
       n1 = size(oldLSAOTENSORitem%elms)
       call mem_lstpointer_alloc(newLSAOTENSORitem%elms,n1)
    ELSE
       NULLIFY(newLSAOTENSORitem%elms)
    ENDIF
    IF(ASSOCIATED(oldLSAOTENSORitem%nOrb))THEN
       n1 = size(oldLSAOTENSORitem%nOrb)
       call mem_lstpointer_alloc(newLSAOTENSORitem%nOrb,n1)
    ELSE
       NULLIFY(newLSAOTENSORitem%nOrb)
    ENDIF
    IF(ASSOCIATED(oldLSAOTENSORitem%startLocalOrb))THEN
       n1 = size(oldLSAOTENSORitem%startLocalOrb)
       call mem_lstpointer_alloc(newLSAOTENSORitem%startLocalOrb,n1)
    ELSE
       NULLIFY(newLSAOTENSORitem%startLocalOrb)
    ENDIF
    IF(ASSOCIATED(oldLSAOTENSORitem%startGlobalOrb))THEN
       n1 = size(oldLSAOTENSORitem%startGlobalOrb)
       call mem_lstpointer_alloc(newLSAOTENSORitem%startGlobalOrb,n1)
    ELSE
       NULLIFY(newLSAOTENSORitem%startGlobalOrb)
    ENDIF
    IF(ASSOCIATED(oldLSAOTENSORitem%nAngmom))THEN
       n1 = size(oldLSAOTENSORitem%nAngmom)
       call mem_lstpointer_alloc(newLSAOTENSORitem%nAngmom,n1)
    ELSE
       NULLIFY(newLSAOTENSORitem%nAngmom)
    ENDIF
    n1 = size(oldLSAOTENSORitem%nLocal)
    DO I=1,n1
       newLSAOTENSORitem%nLocal(I) = oldLSAOTENSORitem%nLocal(I)
    ENDDO
    n1 = size(oldLSAOTENSORitem%ATOM)
    DO I=1,n1
       newLSAOTENSORitem%ATOM(I) = oldLSAOTENSORitem%ATOM(I)
    ENDDO
    n1 = size(oldLSAOTENSORitem%AOBATCH)
    DO I=1,n1
       newLSAOTENSORitem%AOBATCH(I) = oldLSAOTENSORitem%AOBATCH(I)
    ENDDO
  End Subroutine LSAOTENSOR_p_copy

  ! The Primitive subroutines for TYPE:LSTENSOR
  Subroutine LSTENSOR_print(LSTENSORitem,lupri)
    implicit none
    TYPE(LSTENSOR),intent(in) :: LSTENSORitem
    integer,intent(in) :: lupri
    !
    integer :: I,J,K,L,M,n1,IATOM,I2,maxBat,maxAng,n2
    integer :: JATOM,JATOMGL,IATOMGL,GI
    type(slsaotensor),pointer :: SLSAOTENSORitem
    type(lsaotensor),pointer :: LSAOTENSORitem
    type(globallsaotensor),pointer :: G_LSAOTENSORitem
!    type(matrix) :: mat

!!$    IF(LSTENSORitem%screenTensor)THEN
!!$       IF(ASSOCIATED(LSTENSORitem%maxgab))THEN
!!$          n1 = size(LSTENSORitem%maxgab,1)
!!$          n2 = size(LSTENSORitem%maxgab,2)
!!$          Write(lupri,*)'The MaxGab'
!!$          call shortint_output(LSTENSORitem%maxgab,n1,n2,lupri)
!!$       ELSE
!!$          WRITE(lupri,'(A16,A)')' maxgab ' , ' NOT ASSOCIATED ' 
!!$       ENDIF
!!$       IF(ASSOCIATED(LSTENSORitem%maxprimgab))THEN
!!$          n1 = size(LSTENSORitem%maxprimgab,1)
!!$          n2 = size(LSTENSORitem%maxprimgab,2)
!!$          Write(lupri,*)'The MaxPrimGab'
!!$          call shortint_output(LSTENSORitem%maxprimgab,n1,n2,lupri)
!!$       ELSE
!!$          WRITE(lupri,'(A16,A)')' maxprimgab ' , ' NOT ASSOCIATED ' 
!!$       ENDIF
!!$    ENDIF
!    IF(LSTENSORitem%nbast(1).EQ.LSTENSORitem%nbast(2))THEN
!       IF(LSTENSORitem%nbast(3).EQ.1)THEN
!          IF(LSTENSORitem%nbast(4).EQ.1)THEN
!             call mat_init(MAT,LSTENSORitem%nbast(1),LSTENSORitem%nbast(2))
!             call Build_singlemat_from_lst(lupri,LSTENSORitem,MAT)
!             call mat_print(MAT,1,MAT%nrow,1,MAT%nrow,lupri)
!             call mat_free(MAT)
!          ENDIF
!       ENDIF
!    ELSE


    WRITE(lupri,'(A)') ' The LSTENSOR structure'
#ifdef VAR_MPI
    IF(ASSOCIATED(LSTENSORitem%G_LSAO))THEN
       WRITE(lupri,'(A,I9)') ' The Global LSTENSOR structure  G_nLSAO=',LSTENSORitem%G_nLSAO
       DO I=1,LSTENSORitem%G_nLSAO
          WRITE(lupri,'(A,I6,A)')' LSTENSORitem%G_LSAO( ',I,')' 
          G_LSAOTENSORitem => LSTENSORitem%G_LSAO(I)
          WRITE(lupri,'(A)') ' The Global LSAOTENSOR structure'
          WRITE(lupri,'(A16,10I8/,(16X,10I8))')' G_startGlobalOrb' ,&
               & (G_LSAOTENSORitem%G_startGlobalOrb(K),K=1,2)
          WRITE(lupri,'(A16,10I8/,(16X,10I8))')' G_nLocal' ,&
               &(G_LSAOTENSORitem%G_nLocal(I2),I2=1,size(G_LSAOTENSORitem%G_nLocal))
          WRITE(lupri,'(A16,I12)')' MYNUM' ,G_LSAOTENSORitem%mynum
       ENDDO

       WRITE(lupri,'(A16,10I8)')'G_maxnlocal' ,LSTENSORitem%G_maxnlocal
       WRITE(lupri,'(A7,10I8/,(7X,10I8))')'G_nAtom' ,(LSTENSORitem%G_nAtom(I2),I2=1,size(LSTENSORitem%G_nAtom))
       WRITE(lupri,'(A7,10I8/,(7X,10I8))')'G_nbast' ,(LSTENSORitem%G_nbast(I2),I2=1,size(LSTENSORitem%G_nbast))

       IF(ASSOCIATED(LSTENSORitem%G_nAOBATCH))THEN
          DO K=1,2
             n1 = LSTENSORitem%G_nAtom(K) 
             WRITE(lupri,'(A10,I1,A4,10I8/,(15X,10I8))')'G_nAOBATCH',K,':   ',&
                  &(LSTENSORitem%G_nAOBATCH(IATOM,K),IATOM=1,n1)
          ENDDO
       ELSE
          WRITE(lupri,'(A)')'G_nAOBATCH not associated'
       ENDIF
       IF(ASSOCIATED(LSTENSORitem%G_INDEX))THEN
          DO J=1,size(LSTENSORitem%G_INDEX,2)
             WRITE(lupri,'(A7,10I8/,(7X,10I8))')'G_INDEX' ,(LSTENSORitem%G_INDEX(I2,J),I2=1,size(LSTENSORitem%G_INDEX,1))
          ENDDO
       ELSE
          WRITE(lupri,'(A16,A)')' INDEX ' , ' NOT ASSOCIATED ' 
       ENDIF
    ELSE
       WRITE(lupri,'(A)') ' The Global LSTENSOR structure not used'
    ENDIF
#endif

    IF(ASSOCIATED(LSTENSORitem%LSAO))THEN
       WRITE(lupri,'(A)')' LSAOTENSOR structure' 
       DO I=1,LSTENSORitem%nLSAO
          WRITE(lupri,'(A,I6,A)')' LSTENSORitem%LSAO( ',I,')' 
          LSAOTENSORitem => LSTENSORitem%LSAO(I)
          maxBat = LSAOTENSORitem%maxBat
!          maxAng = LSAOTENSORitem%maxAng
          WRITE(lupri,'(A)') ' The LSAOTENSOR structure'
          WRITE(lupri,'(A16,I8)')' nelms ' ,LSAOTENSORitem%nelms
          IF(ASSOCIATED(LSAOTENSORitem%elms))THEN
             call ls_output(LSAOTENSORitem%elms,1,LSAOTENSORitem%nLocal(1),1,LSAOTENSORitem%nLocal(2),&
                  & LSAOTENSORitem%nLocal(1),LSAOTENSORitem%nLocal(2),1,lupri)
          ELSE
             WRITE(lupri,'(A16,A)')' elms ' , ' NOT ASSOCIATED ' 
          ENDIF
          IF(ASSOCIATED(LSAOTENSORitem%nOrb))THEN
!             DO K=1,size(LSAOTENSORitem%nOrb,3)
!                IATOM = LSAOTENSORitem%ATOM(K)
!                n2 = LSTENSORitem%nAOBATCH(IATOM,K)
!                DO J=1,n2
!                   n1 = LSAOTENSORitem%nAngmom(J,K)
!                   WRITE(lupri,'(A16,10I8/,(16X,10I8))')' nOrb ' ,(LSAOTENSORitem%nOrb(I2+(J-1)*maxAng+(K-1)*maxAng*maxBat),I2=1,n1)
!                ENDDO
!             ENDDO
          ELSE
             WRITE(lupri,'(A16,A)')' nOrb ' , ' NOT ASSOCIATED ' 
          ENDIF
          IF(ASSOCIATED(LSAOTENSORitem%startLocalOrb))THEN
!             DO K=1,size(LSAOTENSORitem%startLocalOrb,3)
!                IATOM = LSAOTENSORitem%ATOM(K)
!                n2 = LSTENSORitem%nAOBATCH(IATOM,K)
!                DO J=1,n2
!                   n1 = LSAOTENSORitem%nAngmom(J,K)
!                   WRITE(lupri,'(A16,10I8/,(16X,10I8))')' startLocalOrb ' ,&
!                        & (LSAOTENSORitem%startLocalOrb(I2+(J-1)*maxAng+(K-1)*maxAng*maxBat),I2=1,n1)
!                ENDDO
!             ENDDO
          ELSE
             WRITE(lupri,'(A16,A)')' startLocalOrb ' , ' NOT ASSOCIATED ' 
          ENDIF
          IF(ASSOCIATED(LSAOTENSORitem%startGlobalOrb))THEN
!             DO K=1,size(LSAOTENSORitem%startGlobalOrb,3)
!                IATOM = LSAOTENSORitem%ATOM(K)
!                n2 = LSTENSORitem%nAOBATCH(IATOM,K)
!                DO J=1,n2
!                   n1 = LSAOTENSORitem%nAngmom(J,K)
!                   WRITE(lupri,'(A16,10I8/,(16X,10I8))')' startGlobalOrb ' ,&
!                        & (LSAOTENSORitem%startGlobalOrb(I2+(J-1)*maxAng+(K-1)*maxAng*maxBat),I2=1,n1)
!                ENDDO
!             ENDDO
          ELSE
             WRITE(lupri,'(A16,A)')' startGlobalOrb ' , ' NOT ASSOCIATED ' 
          ENDIF
          IF(ASSOCIATED(LSAOTENSORitem%nAngmom))THEN
!             DO J=1,size(LSAOTENSORitem%nAngmom,2)
!                WRITE(lupri,'(A16,10I8/,(16X,10I8))')' nAngmom ' ,&
!                     &(LSAOTENSORitem%nAngmom(I2+(J-1)*maxAng),I2=1,size(LSAOTENSORitem%nAngmom,1))
!             ENDDO
          ELSE
             WRITE(lupri,'(A16,A)')' nAngmom ' , ' NOT ASSOCIATED ' 
          ENDIF
#ifdef VAR_MPI
          IF(ASSOCIATED(LSTENSORitem%G_fullatoms1).AND.ASSOCIATED(LSTENSORitem%G_fullatoms1))THEN
             IATOM = LSAOTENSORitem%ATOM(1)
             JATOM = LSAOTENSORitem%ATOM(2)
             IATOMGL = LSTENSORitem%G_fullatoms1(IATOM)
             JATOMGL = LSTENSORitem%G_fullatoms2(JATOM)
             WRITE(lupri,'(A,I8,A,I8)')'Local Atoms : ',IATOM,',',JATOM
             WRITE(lupri,'(A,I8,A,I8)')'Global Atoms: ',IATOMGL,',',JATOMGL
             IF(associated(Lstensoritem%g_INDEX))THEN
                gI = Lstensoritem%g_INDEX(IATOMGL,JATOMGL)
                WRITE(lupri,'(A,I8)')'Global Index=',gI
             ENDIF
             IF(associated(Lstensoritem%g_LSAO))THEN
                WRITE(lupri,'(A,I8,A,I8)')'Global start: ',Lstensoritem%g_LSAO(gI)%G_startGlobalOrb(1),',',&
                     &Lstensoritem%g_LSAO(gI)%G_startGlobalOrb(2)
             ENDIF
          ELSE
#endif
             WRITE(lupri,'(A16,10I8/,(16X,10I8))')' nLocal ' ,&
                  &(LSAOTENSORitem%nLocal(I2),I2=1,2)!size(LSAOTENSORitem%nLocal))
             WRITE(lupri,'(A16,10I8/,(16X,10I8))')' ATOM ' ,&
                  &(LSAOTENSORitem%ATOM(I2),I2=1,size(LSAOTENSORitem%ATOM))
             WRITE(lupri,'(A16,10I8/,(16X,10I8))')' AOBATCH ' ,&
                  & (LSAOTENSORitem%AOBATCH(I2),I2=1,size(LSAOTENSORitem%AOBATCH))          
#ifdef VAR_MPI
          ENDIF
#endif
       ENDDO
    ELSE
       WRITE(lupri,'(A16,A)')' LSAO ' , ' NOT ASSOCIATED ' 
    ENDIF
    IF(ASSOCIATED(LSTENSORitem%SLSAO))THEN
       WRITE(lupri,'(A)')' SLSAOTENSOR structure' 
       DO I=1,LSTENSORitem%nSLSAO
          WRITE(lupri,'(A,I6,A)')' LSTENSORitem%SLSAO( ',I,')' 
          SLSAOTENSORitem => LSTENSORitem%SLSAO(I)
          maxBat = SLSAOTENSORitem%maxBat
          WRITE(lupri,'(A)') ' The SLSAOTENSOR structure'
          WRITE(lupri,'(A15,I10)')' nelms ' ,SLSAOTENSORitem%nelms
          IF(ASSOCIATED(SLSAOTENSORitem%selms))THEN
             call shortint_output(SLSAOTENSORitem%selms,SLSAOTENSORitem%nLocal(1),SLSAOTENSORitem%nLocal(2),lupri)
          ELSE
             WRITE(lupri,'(A15,A)')' selms ' , ' NOT ASSOCIATED ' 
          ENDIF
          IF(ASSOCIATED(SLSAOTENSORitem%nOrb))THEN
             DO J=1,2                
                IATOM = SLSAOTENSORitem%ATOM(J)
                n1 = LSTENSORitem%nAOBATCH(IATOM,J)
                WRITE(lupri,'(A15,10I8/,(15X,10I8))')' nOrb ' ,(SLSAOTENSORitem%nOrb(I2+(J-1)*maxBat),I2=1,n1)
             ENDDO
          ELSE
             WRITE(lupri,'(A15,A)')' nOrb ' , ' NOT ASSOCIATED ' 
          ENDIF
          IF(ASSOCIATED(SLSAOTENSORitem%startLocalOrb))THEN
!             DO J=1,size(SLSAOTENSORitem%startLocalOrb,2)
!                IATOM = SLSAOTENSORitem%ATOM(J)
!                n1 = LSTENSORitem%nAOBATCH(IATOM,J)
!                WRITE(lupri,'(A15,10I8/,(15X,10I8))')' startLocalOrb ' ,&
!                     &(SLSAOTENSORitem%startLocalOrb(I2+(J-1)*maxBat),I2=1,n1)
!             ENDDO
          ELSE
             WRITE(lupri,'(A15,A)')' startLocalOrb ' , ' NOT ASSOCIATED ' 
          ENDIF
          WRITE(lupri,'(A15,10I8/,(15X,10I8))')' nLocal ' ,&
               & (SLSAOTENSORitem%nLocal(I2),I2=1,size(SLSAOTENSORitem%nLocal))
          WRITE(lupri,'(A15,10I8/,(15X,10I8))')' ATOM ' ,&
               & (SLSAOTENSORitem%ATOM(I2),I2=1,size(SLSAOTENSORitem%ATOM))
          WRITE(lupri,'(A15,10I8/,(15X,10I8))')' AOBATCH ' ,&
               & (SLSAOTENSORitem%AOBATCH(I2),I2=1,size(SLSAOTENSORitem%AOBATCH))
       ENDDO
    ELSE
       WRITE(lupri,'(A16,A)')' SLSAO ' , ' NOT ASSOCIATED ' 
    ENDIF
    WRITE(lupri,'(A7,10I8/,(7X,10I8))')' nAtom ' ,(LSTENSORitem%nAtom(I2),I2=1,size(LSTENSORitem%nAtom))
    WRITE(lupri,'(A7,10I8/,(7X,10I8))')' nbast ' ,(LSTENSORitem%nbast(I2),I2=1,size(LSTENSORitem%nbast))
    WRITE(lupri,'(A16,I8)')' ndim5 ' ,LSTENSORitem%ndim5
    WRITE(lupri,'(A10,10I8/,(10X,10I8))')' nbatches ' ,(LSTENSORitem%nbatches(I2),I2=1,size(LSTENSORitem%nbatches))
    WRITE(lupri,'(A16,I8)')' nLSAO ' ,LSTENSORitem%nLSAO
    WRITE(lupri,'(A16,I8)')' nSLSAO ' ,LSTENSORitem%nSLSAO
    IF(ASSOCIATED(LSTENSORitem%INDEX))THEN
       DO L=1,size(LSTENSORitem%INDEX,4)
          DO K=1,size(LSTENSORitem%INDEX,3)
             DO J=1,size(LSTENSORitem%INDEX,2)
                WRITE(lupri,'(A7,10I8/,(7X,10I8))')' INDEX ' ,(LSTENSORitem%INDEX(I2,J,K,L),I2=1,size(LSTENSORitem%INDEX,1))
             ENDDO
          ENDDO
       ENDDO
    ELSE
       WRITE(lupri,'(A16,A)')' INDEX ' , ' NOT ASSOCIATED ' 
    ENDIF
    IF(ASSOCIATED(LSTENSORitem%nAOBATCH))THEN
       DO K=1,2
          n1 = LSTENSORitem%nAtom(K) 
          WRITE(lupri,'(A10,I1,A4,10I8/,(15X,10I8))')'  nAOBATCH',K,':   ',&
               &(LSTENSORitem%nAOBATCH(IATOM,K),IATOM=1,n1)
       ENDDO
!       IF(.NOT.LSTENSORitem%Screentensor)THEN
!          DO K=3,4
!             n1 = LSTENSORitem%nAtom(K) 
!             WRITE(lupri,'(A10,I1,A4,10I8/,(15X,10I8))')'  nAOBATCH',K,':   ',&
!                  &(LSTENSORitem%nAOBATCH(IATOM,K),IATOM=1,n1)
!          ENDDO
!       ENDIF
    ENDIF
    WRITE(lupri,'(A16,L1)')' MagGradienttensor' ,LSTENSORitem%MagGradienttensor
    WRITE(lupri,'(A16,L1)')' Gradienttensor ' ,LSTENSORitem%Gradienttensor
    WRITE(lupri,'(A16,L1)')' pChargetensor ' ,LSTENSORitem%pChargetensor
    WRITE(lupri,'(A16,L1)')' EcontribTensor' ,LSTENSORitem%Econtrib
    WRITE(lupri,'(A16,L1)')' Screentensor ' ,LSTENSORitem%Screentensor
    IF(ASSOCIATED(LSTENSORitem%maxgab))THEN
       n1 = size(LSTENSORitem%maxgab,1)
       n2 = size(LSTENSORitem%maxgab,2)
!       DO J=1,size(LSTENSORitem%maxgab,2)
!          WRITE(lupri,'(A16,10I8/,(16X,10I8))')' maxgab ' ,(LSTENSORitem%maxgab(I,J),I=1,size(LSTENSORitem%maxgab,1))
!       ENDDO
       Write(lupri,*)'The MaxGab'
       call shortint_output(LSTENSORitem%maxgab,n1,n2,lupri)
    ELSE
       WRITE(lupri,'(A16,A)')' maxgab ' , ' NOT ASSOCIATED ' 
    ENDIF
    IF(ASSOCIATED(LSTENSORitem%maxprimgab))THEN
       n1 = size(LSTENSORitem%maxprimgab,1)
       n2 = size(LSTENSORitem%maxprimgab,2)
!       DO J=1,size(LSTENSORitem%maxprimgab,2)
!          WRITE(lupri,'(A16,10I8/,(16X,10I8))')' maxprimgab ' ,(LSTENSORitem%maxprimgab(I,J),I=1,size(LSTENSORitem%maxprimgab,1))
!       ENDDO
       Write(lupri,*)'The MaxPrimGab'
       call shortint_output(LSTENSORitem%maxprimgab,n1,n2,lupri)
    ELSE
       WRITE(lupri,'(A16,A)')' maxprimgab ' , ' NOT ASSOCIATED ' 
    ENDIF
    WRITE(lupri,'(A16,I8)')' maxgabelm ' ,LSTENSORitem%maxgabelm
    WRITE(lupri,'(A16,I8)')' maxprimgabelm ' ,LSTENSORitem%maxprimgabelm
    !ENDIF

    Write(lupri,*)'The nMBIE ',LSTENSORitem%nMBIE
    IF(ASSOCIATED(LSTENSORitem%MBIE))THEN
       n1 = size(LSTENSORitem%MBIE,2)
       n2 = size(LSTENSORitem%MBIE,3)
       DO I=1,LSTENSORitem%nMBIE
          call ls_output(LSTENSORitem%MBIE(I,:,:),1,n1,1,n2,n1,n2,1,lupri)
       ENDDO
    ELSE
       WRITE(lupri,'(A16,A)')' MBIE ' , ' NOT ASSOCIATED ' 
    ENDIF

  End Subroutine LSTENSOR_print

  Subroutine LSTENSOR_nullify(LSTENSORitem)
    implicit none
    TYPE(LSTENSOR) :: LSTENSORitem
    !
    integer :: I,J,K,L,M
#ifdef VAR_MPI
    NULLIFY(LSTENSORitem%G_LSAO)
#endif
    NULLIFY(LSTENSORitem%LSAO)
    NULLIFY(LSTENSORitem%SLSAO)
    DO I=1,size(LSTENSORitem%nAtom,1)
       LSTENSORitem%nAtom(I) =0
    ENDDO
    DO I=1,size(LSTENSORitem%nbast,1)
       LSTENSORitem%nbast(I) =0
    ENDDO
    LSTENSORitem%ndim5=0
    DO I=1,size(LSTENSORitem%nbatches,1)
       LSTENSORitem%nbatches(I) =0
    ENDDO
    LSTENSORitem%nLSAO=0
    LSTENSORitem%nSLSAO=0
    NULLIFY(LSTENSORitem%INDEX)
#ifdef VAR_MPI
    NULLIFY(LSTENSORitem%G_nAOBATCH)
    NULLIFY(LSTENSORitem%G_INDEX)
    NULLIFY(LSTENSORitem%G_fullatoms1)
    NULLIFY(LSTENSORitem%G_fullatoms2)
#endif
    LSTENSORitem%MagGradienttensor=.FALSE.
    LSTENSORitem%Gradienttensor=.FALSE.
    LSTENSORitem%pChargetensor=.FALSE.
    LSTENSORitem%Econtrib=.FALSE.
    LSTENSORitem%Screentensor=.FALSE.
    LSTENSORitem%LowerDiagZero = .FALSE.
    LSTENSORitem%PermuteResultTensor = .FALSE.
    NULLIFY(LSTENSORitem%maxgab)
    NULLIFY(LSTENSORitem%maxprimgab)
    NULLIFY(LSTENSORitem%MBIE)
    NULLIFY(LSTENSORitem%nAOBATCH)
    LSTENSORitem%maxgabelm=shortzero
    LSTENSORitem%maxprimgabelm=shortzero
    LSTENSORitem%nMBIE=0
    LSTENSORitem%lstmem_index=0
  End Subroutine LSTENSOR_nullify

  Subroutine LSTENSOR_copy(oldLSTENSORitem,newLSTENSORitem,lupri)
    implicit none
    TYPE(LSTENSOR),intent(in) :: oldLSTENSORitem
    TYPE(LSTENSOR),intent(inout) :: newLSTENSORitem
    integer,intent(in) :: lupri
    !
    integer(kind=long) :: nmemsize,AllocInt,AllocRealk,AllocIntS,nsize
    integer :: I,J,K,L,M,n1,n2,n3,n4,n5
    call lstensor_nullify(newLSTENSORitem)
    IF(ASSOCIATED(oldLSTENSORitem%LSAO).OR.ASSOCIATED(oldLSTENSORitem%SLSAO))THEN
       call retrieve_lstmemval(oldLSTENSORitem%lstmem_index,AllocInt,AllocRealk,AllocIntS)
       call init_lstensorMem(AllocInt,AllocRealk,AllocIntS,newLSTENSORitem%lstmem_index)
       call copy_lstensorMemToCurrent(oldLSTENSORitem%lstmem_index)
    ENDIF
    IF(ASSOCIATED(oldLSTENSORitem%LSAO))THEN
       n1 = size(oldLSTENSORitem%LSAO,1)
       CALL MEM_ALLOC(newLSTENSORitem%LSAO,n1)
       DO I=1,oldLSTENSORitem%nLSAO
          !set pointers and stuff
          call LSAOTENSOR_p_copy(oldLSTENSORitem%LSAO(I),newLSTENSORitem%LSAO(I),lupri)
       ENDDO
    ELSE
       NULLIFY(newLSTENSORitem%LSAO)
    ENDIF
    IF(ASSOCIATED(oldLSTENSORitem%SLSAO))THEN
       n1 = size(oldLSTENSORitem%SLSAO,1)
       CALL MEM_ALLOC(newLSTENSORitem%SLSAO,n1)
       DO I=1,oldLSTENSORitem%nSLSAO
          call SLSAOTENSOR_p_copy(oldLSTENSORitem%SLSAO(I),newLSTENSORitem%SLSAO(I),lupri)
       ENDDO
    ELSE
       NULLIFY(newLSTENSORitem%SLSAO)
    ENDIF
    n1 = size(oldLSTENSORitem%nAtom,1)
    DO I=1,n1
       newLSTENSORitem%nAtom(I) = oldLSTENSORitem%nAtom(I)
    ENDDO
    n1 = size(oldLSTENSORitem%nbast,1)
    DO I=1,n1
       newLSTENSORitem%nbast(I) = oldLSTENSORitem%nbast(I)
    ENDDO
    newLSTENSORitem%ndim5= oldLSTENSORitem%ndim5
    n1 = size(oldLSTENSORitem%nbatches,1)
    DO I=1,n1
       newLSTENSORitem%nbatches(I) = oldLSTENSORitem%nbatches(I)
    ENDDO
    newLSTENSORitem%nLSAO= oldLSTENSORitem%nLSAO
    newLSTENSORitem%nSLSAO= oldLSTENSORitem%nSLSAO
    IF(ASSOCIATED(oldLSTENSORitem%INDEX))THEN
       n1 = size(oldLSTENSORitem%INDEX,1)
       n2 = size(oldLSTENSORitem%INDEX,2)
       n3 = size(oldLSTENSORitem%INDEX,3)
       n4 = size(oldLSTENSORitem%INDEX,4)
       call mem_alloc(newLSTENSORitem%INDEX,n1,n2,n3,n4)
       nsize = n1*n2*n3*n4*mem_intsize
       call mem_allocated_mem_lstensor(nsize)
       DO L=1,n4
          DO K=1,n3
             DO J=1,n2
                DO I=1,n1
                   newLSTENSORitem%INDEX(I,J,K,L) = oldLSTENSORitem%INDEX(I,J,K,L)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ELSE
       NULLIFY(newLSTENSORitem%INDEX)
    ENDIF
    newLSTENSORitem%MagGradienttensor= oldLSTENSORitem%MagGradienttensor
    newLSTENSORitem%Gradienttensor= oldLSTENSORitem%Gradienttensor
    newLSTENSORitem%pChargetensor= oldLSTENSORitem%pChargetensor
    newLSTENSORitem%Econtrib= oldLSTENSORitem%Econtrib
    newLSTENSORitem%Screentensor= oldLSTENSORitem%Screentensor
    IF(ASSOCIATED(oldLSTENSORitem%maxgab))THEN
       n1 = size(oldLSTENSORitem%maxgab,1)
       n2 = size(oldLSTENSORitem%maxgab,2)
       call mem_alloc(newLSTENSORitem%maxgab,n1,n2)
       nsize = n1*n2*mem_shortintsize
       call mem_allocated_mem_lstensor(nsize)

       DO J=1,n2
          DO I=1,n1
             newLSTENSORitem%maxgab(I,J) = oldLSTENSORitem%maxgab(I,J)
          ENDDO
       ENDDO
    ELSE
       NULLIFY(newLSTENSORitem%maxgab)
    ENDIF
    IF(ASSOCIATED(oldLSTENSORitem%maxprimgab))THEN
       n1 = size(oldLSTENSORitem%maxprimgab,1)
       n2 = size(oldLSTENSORitem%maxprimgab,2)
       call mem_alloc(newLSTENSORitem%maxprimgab,n1,n2)
       nsize = n1*n2*mem_shortintsize
       call mem_allocated_mem_lstensor(nsize)

       DO J=1,n2
          DO I=1,n1
             newLSTENSORitem%maxprimgab(I,J) = oldLSTENSORitem%maxprimgab(I,J)
          ENDDO
       ENDDO
    ELSE
       NULLIFY(newLSTENSORitem%maxprimgab)
    ENDIF
    newLSTENSORitem%nMBIE = oldLSTENSORitem%nMBIE
    IF(ASSOCIATED(oldLSTENSORitem%MBIE))THEN
       n1 = size(oldLSTENSORitem%MBIE,2)
       n2 = size(oldLSTENSORitem%MBIE,3)
       call mem_alloc(newLSTENSORitem%MBIE,newLSTENSORitem%nMBIE,n1,n2)
       nsize = SIZE(newLSTENSORitem%MBIE,KIND=long)*mem_realsize
       call mem_allocated_mem_lstensor(nsize)
       DO K=1,n2
        DO J=1,n1
         DO I=1,newLSTENSORitem%nMBIE
          newLSTENSORitem%MBIE(I,J,K)=oldLSTENSORitem%MBIE(I,J,K)
         ENDDO
        ENDDO
       ENDDO
    ELSE
       NULLIFY(newLSTENSORitem%MBIE)
    ENDIF
    newLSTENSORitem%maxgabelm = oldLSTENSORitem%maxgabelm
    newLSTENSORitem%maxprimgabelm = oldLSTENSORitem%maxprimgabelm
    IF(ASSOCIATED(oldLSTENSORitem%nAOBATCH))THEN
       n1 = size(oldLSTENSORitem%nAOBATCH,1)
       n2 = size(oldLSTENSORitem%nAOBATCH,2)
       call mem_alloc(newLSTENSORitem%nAOBATCH,n1,n2)
       nsize = n1*n2*mem_intsize
       call mem_allocated_mem_lstensor(nsize)
       DO J=1,n2
          DO I=1,n1
             newLSTENSORitem%nAOBATCH(I,J) = oldLSTENSORitem%nAOBATCH(I,J)
          ENDDO
       ENDDO
    ELSE
       NULLIFY(newLSTENSORitem%nAOBATCH)
    ENDIF

  End Subroutine LSTENSOR_copy

  Subroutine LSTENSOR_free(LSTENSORitem)
    implicit none
    TYPE(LSTENSOR) :: LSTENSORitem
    integer(kind=long) :: nsize
    !
    integer :: I,J,K,L,M

    call LSTENSOR_local_free(LSTENSORitem)

    IF(ASSOCIATED(LSTENSORitem%maxgab))THEN
       nsize = size(LSTENSORitem%maxgab,KIND=long)*mem_shortintsize
       call mem_deallocated_mem_lstensor(nsize)
       CALL MEM_DEALLOC(LSTENSORitem%maxgab)
    ENDIF
    IF(ASSOCIATED(LSTENSORitem%maxprimgab))THEN
       nsize = size(LSTENSORitem%maxprimgab,KIND=long)*mem_shortintsize
       call mem_deallocated_mem_lstensor(nsize)
       CALL MEM_DEALLOC(LSTENSORitem%maxprimgab)
    ENDIF
    IF(ASSOCIATED(LSTENSORitem%MBIE))THEN
       nsize = size(LSTENSORitem%MBIE,KIND=long)*mem_realsize
       call mem_deallocated_mem_lstensor(nsize)
       CALL MEM_DEALLOC(LSTENSORitem%MBIE)
    ENDIF
#ifdef VAR_MPI
    IF(ASSOCIATED(LSTENSORitem%G_LSAO))THEN
       CALL MEM_DEALLOC(LSTENSORitem%G_LSAO)
    ENDIF
    IF(ASSOCIATED(LSTENSORitem%G_INDEX))THEN
       CALL MEM_DEALLOC(LSTENSORitem%G_INDEX)
    ENDIF
    IF(ASSOCIATED(LSTENSORitem%G_nAOBATCH))THEN
       call mem_dealloc(LSTENSORitem%G_nAOBATCH)
    ENDIF
#endif
  End Subroutine LSTENSOR_free

  Subroutine LSTENSOR_local_free(LSTENSORitem)
    implicit none
    TYPE(LSTENSOR) :: LSTENSORitem
    integer(kind=long) :: nsize
    !
    integer :: I,J,K,L,M

    IF(ASSOCIATED(LSTENSORitem%LSAO).OR.ASSOCIATED(LSTENSORitem%SLSAO))THEN
       call free_lstensorMem(LSTENSORitem%lstmem_index)
    ENDIF
    IF(ASSOCIATED(LSTENSORitem%LSAO))THEN
       CALL MEM_DEALLOC(LSTENSORitem%LSAO)
    ENDIF
    IF(ASSOCIATED(LSTENSORitem%SLSAO))THEN
       CALL MEM_DEALLOC(LSTENSORitem%SLSAO)
    ENDIF
    IF(ASSOCIATED(LSTENSORitem%INDEX))THEN
       nsize = size(LSTENSORitem%INDEX,KIND=long)*mem_intsize
       call mem_deallocated_mem_lstensor(nsize)
       CALL MEM_DEALLOC(LSTENSORitem%INDEX)
    ENDIF
    IF(ASSOCIATED(LSTENSORitem%nAOBATCH))THEN
       nsize = size(LSTENSORitem%nAOBATCH,KIND=long)*mem_intsize
       call mem_deallocated_mem_lstensor(nsize)
       CALL MEM_DEALLOC(LSTENSORitem%nAOBATCH)
    ENDIF
#ifdef VAR_MPI
    IF(ASSOCIATED(LSTENSORitem%G_fullatoms1))THEN
       CALL MEM_DEALLOC(LSTENSORitem%G_fullatoms1)
    ENDIF
    IF(ASSOCIATED(LSTENSORitem%G_fullatoms2))THEN
       CALL MEM_DEALLOC(LSTENSORitem%G_fullatoms2)
    ENDIF
#endif
  End Subroutine LSTENSOR_local_free

  !END PRIMITIVE SUBROUTINE

  !> \brief set the lstensor to zero
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param TENSOR the lstensor
  SUBROUTINE lstensor_zero(TENSOR)
    implicit none
    TYPE(LSTENSOR)     :: TENSOR
    call zero_lstensorMemindex(TENSOR%lstmem_index)
  END SUBROUTINE LSTENSOR_ZERO

  !> \brief initiate the lstensor structure
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param TENSOR the lstensor
  !> \param AO1 the Atomic orbital(AO) item for center 1 
  !> \param AO2 the AO item for center 2 
  !> \param AO3 the AO item for center 3 
  !> \param AO4 the AO item for center 4 
  !> \param nbast1 the size of the 1. dimension
  !> \param nbast2 the size of the 2. dimension
  !> \param nbast3 the size of the 3. dimension
  !> \param nbast4 the size of the 4. dimension
  !> \param nmat the number of density matrices or size of the 5. dimension
  !> \param useAO1 flag to describe if the AO1 item should be used or an empty AO 
  !> \param useAO2 flag to describe if the AO1 item should be used or an empty AO 
  !> \param useAO3 flag to describe if the AO1 item should be used or an empty AO 
  !> \param useAO4 flag to describe if the AO1 item should be used or an empty AO 
  !> \param ODscreen flag: use the overlap distribution screening 
  !> \param lupri the logical unit number for the output
  SUBROUTINE init_lstensor_5dim(TENSOR,AO1,AO2,AO3,AO4,nbast1,nbast2,nbast3,nbast4,&
       & ndim5,useAO1,useAO2,useAO3,useAO4,ODscreen,ODscreen13,lupri)
    implicit none
    INTEGER            :: nbast1,nbast2,nbast3,nbast4,lupri,ndim5
    TYPE(LSTENSOR)     :: TENSOR
    TYPE(AOITEM),target  :: AO1,AO2,AO3,AO4
    logical   :: useAO1,useAO2,useAO3,useAO4,ODscreen,ODscreen13
    !
    TYPE(AOITEMPOINTER) :: AOT(4)
    TYPE(AOITEM),target :: AOTE
    ! 
    INTEGER :: natom(4),nbat(4),AOTbatch(4),batAOT(4),nAng(4)
    !
    INTEGER :: nbastI,Ibat,Iangmom,Iatom,IORB,nOrbA,sA
    INTEGER :: nbastJ,Jbat,Jangmom,Jatom,JORB,nOrbB,sB
    INTEGER :: nbastK,Kbat,Kangmom,Katom,KORB,nOrbC,sC
    INTEGER :: nbastL,Lbat,Langmom,Latom,LORB,nOrbD,sD
    ! COMMON
    INTEGER :: DIM,IELM,IMAT2,I,J,iAO,jAO,maxBat,maxAng,nLocal(4)
    INTEGER :: nElms,startLocalOrb,nAngG,batAOTG,dim2(4),maxnatom
    integer :: offset,lstmem_index
    LOGICAL :: ERROR,Screen,ScreenAtom
    LOGICAL,pointer :: nonscreenAB(:,:),nonscreenAC(:,:),nonscreenCD(:,:)
    integer(kind=long) :: nmemsize,AllocInt,AllocIntS,AllocRealk,nsize
    TYPE(LSAOTENSOR),pointer    :: lsao
!    integer,pointer :: inttmp(:,:)
    integer,pointer :: inttmp2(:)
    real(realk) :: ts,te
    call LSTENSOR_nullify(TENSOR)
    TENSOR%pChargetensor = .FALSE.
    TENSOR%Econtrib = .FALSE.
    TENSOR%Screentensor = .FALSE.
    TENSOR%MagGradienttensor = .FALSE.
    TENSOR%Gradienttensor = .FALSE.
    tensor%screenoutput = .FALSE.
    call SET_EMPTY_AO(AOTE)
    AOT(1)%p => AO1
    AOT(2)%p => AO2
    AOT(3)%p => AO3
    AOT(4)%p => AO4
    IF(.NOT.useAO1) AOT(1)%p => AOTE
    IF(.NOT.useAO2) AOT(2)%p => AOTE
    IF(.NOT.useAO3) AOT(3)%p => AOTE
    IF(.NOT.useAO4) AOT(4)%p => AOTE
    IF(nbast1.NE.AOT(1)%p%nbast)then
       print*,'nbast1',nbast1
       print*,'AOT(1)%p%nbast',AOT(1)%p%nbast
       call lsquit('dim1 mismatch in init_lstensor_5dim',-1)
    ENDIF
    IF(nbast2.NE.AOT(2)%p%nbast)call lsquit('dim2 mismatch in init_lstensor_5dim',-1)
    IF(nbast3.NE.AOT(3)%p%nbast)call lsquit('dim3 mismatch in init_lstensor_5dim',-1)
    IF(nbast4.NE.AOT(4)%p%nbast)call lsquit('dim4 mismatch in init_lstensor_5dim',-1)
    DO IAO = 1,4
       natom(iAO) = AOT(iAO)%p%natoms
    ENDDO

    maxnatom = MAX(natom(1),natom(2),natom(3),natom(4))
    call mem_alloc(TENSOR%nAOBATCH,maxnatom,4)
    nsize = maxnatom*4*mem_intsize
    call mem_allocated_mem_lstensor(nsize)
!    inttmp => TENSOR%nAOBATCH
!    DO IAO = 1,4
!       inttmp2 => AOT(IAO)%p%ATOMICnBatch
!       DO IATOM = 1,natom(iAO)
!          inttmp(IATOM,IAO) = inttmp2(IATOM)
!       ENDDO
!    ENDDO
    DO IAO = 1,4
       inttmp2 => AOT(IAO)%p%ATOMICnBatch
       DO IATOM = 1,natom(iAO)
          TENSOR%nAOBATCH(IATOM,IAO) = inttmp2(IATOM)
       ENDDO
    ENDDO
    DO IAO = 1,4
       TENSOR%natom(iAO) = natom(iAO)
    ENDDO
    TENSOR%nbast(1) = nbast1
    TENSOR%nbast(2) = nbast2
    TENSOR%nbast(3) = nbast3
    TENSOR%nbast(4) = nbast4
    TENSOR%ndim5 = ndim5
    DO IAO = 1,4
       TENSOR%nbatches(iAO) = AOT(iAO)%p%nbatches
    ENDDO

    
    call mem_alloc(nonscreenAB,AOT(1)%p%natoms,AOT(2)%p%natoms)
    IF(ODscreen)THEN
       call build_atomicODScreen(AOT(1)%p,AOT(2)%p,nonscreenAB,AOT(1)%p%natoms,AOT(2)%p%natoms)
    else
       nonscreenAB = .true.
    ENDIF

    call mem_alloc(nonscreenCD,AOT(3)%p%natoms,AOT(4)%p%natoms)
    IF(ODscreen)THEN
       call build_atomicODScreen(AOT(3)%p,AOT(4)%p,nonscreenCD,AOT(3)%p%natoms,AOT(4)%p%natoms)
    else
       nonscreenCD = .true.
    ENDIF

    call mem_alloc(nonscreenAC,AOT(1)%p%natoms,AOT(3)%p%natoms)
    IF(ODscreen13)THEN
       call build_atomicODScreen(AOT(1)%p,AOT(3)%p,nonscreenAC,AOT(1)%p%natoms,AOT(3)%p%natoms)
    else
       nonscreenAC = .true.
    ENDIF

    I = 0
    DO Latom = 1,natom(4)
       DO Katom = 1,natom(3)
          IF(nonscreenCD(Katom,Latom))THEN             
             DO Jatom = 1,natom(2)
                DO Iatom = 1,natom(1)
                   IF(nonscreenAB(Iatom,Jatom))THEN                          
                      IF(nonscreenAC(Iatom,Katom))THEN                          
                         I = I +1
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    TENSOR%nLSAO = I
    CALL MEM_ALLOC(TENSOR%INDEX,natom(1),natom(2),natom(3),natom(4))
    nsize = size(TENSOR%INDEX,KIND=long)*mem_intsize
    call mem_allocated_mem_lstensor(nsize)
    DO Latom = 1,natom(4)
       DO Katom = 1,natom(3)
          DO Jatom = 1,natom(2)
             DO Iatom = 1,natom(1)
                TENSOR%INDEX(Iatom,Jatom,Katom,Latom) = 0 !if 0 lsaotensor not call mem_allocd 
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    AllocInt = 0
    AllocIntS = 0
    AllocRealk = 0
    I = 0
    AOTbatch(4)=0
    DO Latom = 1,natom(4)
       nbat(4) = AOT(4)%p%ATOMICnBatch(LATOM)
       AOTbatch(3)=0
       DO Katom = 1,natom(3)
          nbat(3) = AOT(3)%p%ATOMICnBatch(KATOM)
          IF(nonscreenCD(Katom,Latom))THEN             
             AOTbatch(2)=0
             DO Jatom = 1,natom(2)
                nbat(2) = AOT(2)%p%ATOMICnBatch(JATOM)
                AOTbatch(1)=0
                DO Iatom = 1,natom(1)
                   nbat(1) = AOT(1)%p%ATOMICnBatch(IATOM)
                   IF(nonscreenAB(Iatom,Jatom))THEN                          
                      IF(nonscreenAC(Iatom,Katom))THEN                          
                         I=I+1
                         maxBat = nBat(1)
                         maxBat = MAX(maxBat,nBat(2))
                         maxBat = MAX(maxBat,nBat(3))
                         maxBat = MAX(maxBat,nBat(4))
                         AllocInt = AllocInt + 4*maxBat
                         maxAng = 0
                         do IAO=1,4
                            batAOTG = AOTbatch(IAO)
                            DO Ibat = 1,nbat(iAO)
                               batAOTG = batAOTG+1
                               maxAng = MAX(maxAng,AOT(iAO)%p%BATCH(batAOTG)%nAngmom)
                            ENDDO
                         enddo
                         AllocInt = AllocInt + 12*maxBat*maxAng
                         do IAO=1,4
                            DIM2(iAO) = 0
                            batAOTG = AOTbatch(IAO)
                            DO Ibat = 1,nbat(iAO)
                               batAOTG = batAOTG+1
                               DO Iangmom = 1,AOT(IAO)%p%BATCH(batAOTG)%nAngmom
                                  nOrbA = AOT(iAO)%p%BATCH(batAOTG)%nContracted(Iangmom)*&
                                    & AOT(iAO)%p%BATCH(batAOTG)%nOrbComp(Iangmom)
                                  DIM2(iAO) = DIM2(iAO) + nOrbA
                               ENDDO
                            ENDDO
                            IF(DIM2(iAO).EQ.0)DIM2(iAO)=1
                         ENDDO
                         DIM = DIM2(1)*DIM2(2)*DIM2(3)*DIM2(4)*ndim5
                         AllocRealk = AllocRealk + DIM
                      ENDIF !AC screen
                   ENDIF !AB screen
                   AOTbatch(1) = AOTbatch(1) + nbat(1)
                ENDDO
                AOTbatch(2) = AOTbatch(2) + nbat(2)
             ENDDO
          ENDIF
          AOTbatch(3) = AOTbatch(3) + nbat(3)
       ENDDO
       AOTbatch(4) = AOTbatch(4) + nbat(4)
    ENDDO
    TENSOR%nLSAO = I
    CALL MEM_ALLOC(TENSOR%LSAO,I)
    call nullifyTENSORLSAO(TENSOR%LSAO)
!    call lstimer('START',ts,te,6)
    call init_lstensorMem(AllocInt,AllocRealk,AllocInts,lstmem_index)
!    call lstimer('tensmem',ts,te,6)
    call zero_lstensorMem
    TENSOR%lstmem_index = lstmem_index
    I = 0
    AOTbatch(4)=0
    DO Latom = 1,natom(4)
       nbat(4) = AOT(4)%p%ATOMICnBatch(LATOM)
       AOTbatch(3)=0
       DO Katom = 1,natom(3)
          nbat(3) = AOT(3)%p%ATOMICnBatch(KATOM)
          IF(nonscreenCD(Katom,Latom))THEN             
          AOTbatch(2)=0
          DO Jatom = 1,natom(2)
             nbat(2) = AOT(2)%p%ATOMICnBatch(JATOM)
             AOTbatch(1)=0
             DO Iatom = 1,natom(1)
                nbat(1) = AOT(1)%p%ATOMICnBatch(IATOM)
                IF(nonscreenAB(Iatom,Jatom))THEN                          
                   IF(nonscreenAC(Iatom,Katom))THEN                          
                      I=I+1
                      TENSOR%INDEX(IATOM,JATOM,KATOM,Latom) = I
                      lsao => TENSOR%LSAO(I)
                      lsao%ATOM(1) = Iatom
                      lsao%ATOM(2) = Jatom
                      lsao%ATOM(3) = Katom
                      lsao%ATOM(4) = Latom
                      lsao%AOBATCH(1) = AOTbatch(1)
                      lsao%AOBATCH(2) = AOTbatch(2)
                      lsao%AOBATCH(3) = AOTbatch(3)
                      lsao%AOBATCH(4) = AOTbatch(4)
                      maxBat = nBat(1)
                      maxBat = MAX(maxBat,nBat(2))
                      maxBat = MAX(maxBat,nBat(3))
                      maxBat = MAX(maxBat,nBat(4))
                      lsao%maxBat = maxBat
                      call mem_LSTpointer_alloc(lsao%nAngmom,maxBat*4)
                      maxAng = 0
                      do IAO=1,4
                         batAOTG = AOTbatch(IAO)
                         offset = (IAO-1)*maxBat
                         DO Ibat = 1,nbat(iAO)
                            batAOTG = batAOTG+1
                            lsao%nAngmom(Ibat+offset) = AOT(iAO)%p%BATCH(batAOTG)%nAngmom
                            maxAng = MAX(maxAng,AOT(iAO)%p%BATCH(batAOTG)%nAngmom)
                         ENDDO
                      enddo
                      lsao%maxAng = maxAng
                      call mem_LSTpointer_alloc(lsao%nOrb,maxAng*maxBat*4)
                      call mem_LSTpointer_alloc(lsao%startLocalOrb,maxAng*maxBat*4)
                      call mem_LSTpointer_alloc(lsao%startGlobalOrb,maxAng*maxBat*4)
                      do IAO=1,4
                         nLocal(iAO) = 0
                         batAOTG = AOTbatch(IAO)
                         startLocalOrb = 1
                         DO Ibat = 1,nbat(iAO)
                            batAOTG = batAOTG+1
                            nAngG = AOT(iAO)%p%BATCH(batAOTG)%nAngmom 
                            offset = (Ibat-1)*maxAng+(IAO-1)*maxAng*maxBat
                            DO Iangmom = 1,nAngG
                               lsao%nOrb(Iangmom+offset)=AOT(iAO)%p%BATCH(batAOTG)%nContracted(Iangmom)*&
                                    & AOT(iAO)%p%BATCH(batAOTG)%nOrbComp(Iangmom)
                               nLocal(iAO) = nLocal(iAO) + lsao%nOrb(Iangmom+offset)
                               lsao%startGlobalOrb(Iangmom+offset)=&
                                    & AOT(iAO)%p%BATCH(batAOTG)%startOrbital(Iangmom)
                               lsao%startLocalOrb(Iangmom+offset)=startLocalOrb
                               startLocalOrb = startLocalOrb + lsao%nOrb(Iangmom+offset)
                            ENDDO
                         ENDDO
                         IF(nLocal(iAO).EQ.0)nLocal(iAO)=1
                         lsao%nLocal(iAO) = nLocal(iAO)
                      ENDDO
 !WARNING THIS CAN BE REMOVES AND                     
                      do IAO=1,4
                         DIM2(iAO) = 0
                         batAOT(IAO) = AOTbatch(IAO)
                         DO Lbat = 1,nbat(IAO)
                            offset = (Lbat-1)*maxAng+(IAO-1)*maxAng*maxBat
                            batAOT(IAO) = batAOT(IAO)+1
                            nAng(IAO) = AOT(IAO)%p%BATCH(batAOT(IAO))%nAngmom
                            DO Langmom = 1,nAng(IAO)
                               nOrbD = lsao%nOrb(Langmom+offset)
                               DIM2(iAO) = DIM2(iAO) + nOrbD
                            ENDDO
                         ENDDO
                         IF(DIM2(iAO).EQ.0)DIM2(iAO)=1
                      ENDDO
                      DIM = DIM2(1)*DIM2(2)*DIM2(3)*DIM2(4)*ndim5
                      call mem_LSTpointer_alloc(lsao%elms,DIM)
                      lsao%nelms = DIM/ndim5
!                   ENDIF !screenatom
                   ENDIF !AC screen
                ENDIF !AB screen
                AOTbatch(1) = AOTbatch(1) + nbat(1)
             ENDDO
             AOTbatch(2) = AOTbatch(2) + nbat(2)
          ENDDO
       ENDIF
       AOTbatch(3) = AOTbatch(3) + nbat(3)
    ENDDO
    AOTbatch(4) = AOTbatch(4) + nbat(4)
 ENDDO
 !nullify the remaining terms
 DO J = I+1,TENSOR%nLSAO
    call lsaotensor_nullify(TENSOR%LSAO(J))
 ENDDO
 TENSOR%nLSAO = I
 call mem_dealloc(nonscreenAB)
 call mem_dealloc(nonscreenAC)
 call mem_dealloc(nonscreenCD)
 call FREE_EMPTY_AO(AOTE)
 
END SUBROUTINE INIT_LSTENSOR_5DIM

  SUBROUTINE Init_lstensor_1dim(TENSOR,ndim,lupri)
    implicit none
    TYPE(LSTENSOR)     :: TENSOR
    INTEGER            :: ndim,lupri
    !
    INTEGER :: imat,iAO,lstmem_index
    integer(kind=long) :: nmemsize,AllocInt,AllocRealk,AllocIntS,nsize
    TYPE(LSAOTENSOR),pointer    :: lsao
    call LSTENSOR_nullify(TENSOR)
    TENSOR%pChargetensor = .FALSE.
    TENSOR%Econtrib = .TRUE.
    TENSOR%Screentensor = .FALSE.
    TENSOR%Maggradienttensor = .FALSE.
    TENSOR%gradienttensor = .FALSE.
    tensor%screenoutput = .FALSE.
    CALL MEM_ALLOC(TENSOR%LSAO,1)
    AllocInt = 16
    AllocRealk = ndim
    AllocIntS = 0
    call init_lstensorMem(AllocInt,AllocRealk,AllocInts,lstmem_index)
    call zero_lstensorMem
    TENSOR%lstmem_index = lstmem_index
    call lsaotensor_nullify(TENSOR%LSAO(1))
    CALL MEM_ALLOC(TENSOR%INDEX,1,1,1,1)
    nsize = size(TENSOR%INDEX,KIND=long)*mem_intsize
    call mem_allocated_mem_lstensor(nsize)

    TENSOR%INDEX = 1
    DO IAO = 1,4
       TENSOR%natom(iAO) = 1
       TENSOR%nbast(iAO) = 1
       TENSOR%nbatches(iAO) = 1
    ENDDO
    TENSOR%ndim5 = ndim
    TENSOR%nLSAO = 1
    TENSOR%INDEX(1,1,1,1) = 1
    lsao => TENSOR%LSAO(1)
    lsao%ATOM(1) = 1; lsao%ATOM(2) = 1; lsao%ATOM(3) = 1; lsao%ATOM(4) = 1
    call mem_alloc(TENSOR%nAOBATCH,1,4)
    nsize = size(TENSOR%nAOBATCH,KIND=long)*mem_intsize
    call mem_allocated_mem_lstensor(nsize)

    TENSOR%nAOBATCH(1,1) = 1; TENSOR%nAOBATCH(1,2) = 1
    TENSOR%nAOBATCH(1,3) = 1; TENSOR%nAOBATCH(1,4) = 1
    lsao%AOBATCH(1) = 0; lsao%AOBATCH(2) = 0
    lsao%AOBATCH(3) = 0; lsao%AOBATCH(4) = 0
    lsao%nLocal=1
    CALL MEM_LSTPOINTER_ALLOC(lsao%nAngmom,4)
    lsao%nAngmom(1) = 1
    lsao%nAngmom(2) = 1
    lsao%nAngmom(3) = 1
    lsao%nAngmom(4) = 1
    lsao%maxBat = 1; lsao%maxAng = 1
    CALL MEM_LSTPOINTER_ALLOC(lsao%elms,ndim)
    lsao%nelms = 1
    DO IMAT = 1,ndim
       lsao%elms(IMAT)=0E0_realk
    ENDDO
    CALL MEM_LSTPOINTER_ALLOC(lsao%nOrb,4)
    CALL MEM_LSTPOINTER_ALLOC(lsao%startLocalOrb,4)
    CALL MEM_LSTPOINTER_ALLOC(lsao%startGlobalOrb,4)
    do IAO=1,4
       lsao%nOrb(iAO)=1
       lsao%startGlobalOrb(iAO)=1
       lsao%startLocalOrb(iAO)=1
    ENDDO

  end subroutine Init_lstensor_1dim

  !> \brief copy an lstensor to a new lstensor
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param TENSOR1 the original lstensor
  !> \param TENSOR2 the new lstensor
  SUBROUTINE copy_lstensor_to_lstensor(TENSOR1,TENSOR2)
    implicit none
    TYPE(LSTENSOR),intent(in)     :: TENSOR1
    TYPE(LSTENSOR),intent(inout)  :: TENSOR2
    
    call lstensor_copy(TENSOR1,TENSOR2,6)

  end SUBROUTINE copy_lstensor_to_lstensor

  !> \brief copy an lstensor to a new lstensor
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param TENSOR1 the original lstensor
  !> \param TENSOR2 the new sublstensor
  SUBROUTINE build_sublstensor_from_full_lstensor(TENSOR2,TENSOR1,&
       & natoms1,natoms2,natoms3,natoms4,atoms1,atoms2,atoms3,atoms4,&
       & nbast1,nbast2,nbast3,nbast4,sameFrag)
    implicit none
    TYPE(LSTENSOR)     :: TENSOR1,TENSOR2
    integer,intent(in) :: natoms1,natoms2,natoms3,natoms4,nbast1,nbast2,nbast3,nbast4
    integer,intent(in) :: atoms1(natoms1),atoms2(natoms2),atoms3(natoms3),atoms4(natoms4)  
    logical,intent(in) :: SameFrag
    !
    INTEGER   :: I1,Ibat,Jbat,Kbat,Lbat,ndim5,dim,ielms,I2
    INTEGER   :: IATOM1,IATOM2,JATOM1,JATOM2,natomG,iAO
    INTEGER   :: KATOM1,KATOM2,LATOM1,LATOM2,maxnAtom
    INTEGER   :: nAngmomA,nAngmomB,nAngmomC,nAngmomD,natom(4)
    Integer   :: ibatch2,ibatch1,jbatch1,jbatch2,nbatB,nbatA,IATOM,lstmem_index
    Integer   :: nbat1(natoms1),nbat2(natoms2),nstart1(natoms1),nstart2(natoms2)
    integer  :: nbatches1,nbatches2,nbatJ,nbatI,Jbatch,Ibatch,Lbatch,Kbatch
    integer(kind=long) :: nmemsize,AllocInt,AllocRealk,AllocIntS,nsize
    TYPE(LSAOTENSOR),pointer    :: lsao1
    TYPE(LSAOTENSOR),pointer    :: lsao2
!    print*,'TENSOR1 in build_sublstensor_from_full_lstensor'
!    call lstensor_print(TENSOR1,6)
    call LSTENSOR_nullify(TENSOR2)
    natom(1) = natoms1
    natom(2) = natoms2
    natom(3) = natoms3
    natom(4) = natoms4

    TENSOR2%pChargetensor = TENSOR1%pChargetensor
    TENSOR2%Econtrib = TENSOR1%Econtrib
    TENSOR2%Maggradienttensor = TENSOR1%Maggradienttensor
    TENSOR2%gradienttensor = TENSOR1%gradienttensor
    TENSOR2%ScreenTensor = TENSOR1%ScreenTensor
    TENSOR2%ndim5   = TENSOR1%ndim5
    ndim5 = TENSOR2%ndim5

    nullify(TENSOR2%nAOBATCH)
    IF(.NOT.TENSOR2%Econtrib)THEN
       IF(ASSOCIATED(TENSOR1%nAOBATCH))THEN
          maxnatom = MAX(natoms1,natoms2,natoms3,natoms4)
          call mem_alloc(TENSOR2%nAOBATCH,maxnatom,4)
          nsize = size(TENSOR2%nAOBATCH,KIND=long)*mem_intsize
          call mem_allocated_mem_lstensor(nsize)
          DO IATOM = 1,natoms1
             IATOM2 = atoms1(IATOM)
             TENSOR2%nAOBATCH(IATOM,1) = TENSOR1%nAOBATCH(IATOM2,1)
          ENDDO
          DO IATOM = 1,natoms2
             IATOM2 = atoms2(IATOM)
             TENSOR2%nAOBATCH(IATOM,2) = TENSOR1%nAOBATCH(IATOM2,2)
          ENDDO
          IF(.NOT.TENSOR1%ScreenTensor)THEN
             DO IATOM = 1,natoms3
                IATOM2 = atoms3(IATOM)
                TENSOR2%nAOBATCH(IATOM,3) = TENSOR1%nAOBATCH(IATOM2,3)
             ENDDO
             DO IATOM = 1,natoms4
                IATOM2 = atoms4(IATOM)
                TENSOR2%nAOBATCH(IATOM,4) = TENSOR1%nAOBATCH(IATOM2,4)
             ENDDO
          ENDIF
       ENDIF
    ENDIF
    IF(TENSOR2%maggradienttensor)THEN
       TENSOR2%natom(1) = 1
       TENSOR2%natom(2) = 1
       TENSOR2%natom(3) = 1
       TENSOR2%natom(4) = 1
       TENSOR2%nLSAO = 1
       CALL MEM_ALLOC(TENSOR2%LSAO,1)
       CALL MEM_ALLOC(TENSOR2%INDEX,1,1,1,1)
       nsize = size(TENSOR2%INDEX,KIND=long)*mem_intsize
       call mem_allocated_mem_lstensor(nsize)

       TENSOR2%INDEX(1,1,1,1) = 1
       TENSOR2%nbast(1) = 1
       TENSOR2%nbast(2) = TENSOR2%nLSAO
       TENSOR2%nbast(3) = 1
       TENSOR2%nbast(4) = 1
       AllocInt = 0
       AllocIntS = 0
       AllocRealk = 3
       call init_lstensorMem(AllocInt,AllocRealk,AllocInts,lstmem_index)
       call zero_lstensorMem
       TENSOR2%lstmem_index = lstmem_index
       call lsaotensor_nullify(TENSOR2%LSAO(1))
       lsao2 => TENSOR2%LSAO(1) 
       lsao2%ATOM = 1
       lsao2%AOBATCH = 0
       lsao2%nLocal = 0
       lsao2%maxBat = 0
       lsao2%maxAng = 0
       NULLIFY(lsao2%nOrb)
       NULLIFY(lsao2%startLocalOrb)
       NULLIFY(lsao2%startGlobalOrb)  
       NULLIFY(lsao2%nAngmom)
       CALL MEM_LSTPOINTER_ALLOC(lsao2%elms,3)
       lsao2%nelms = 3
    ELSEIF(TENSOR2%gradienttensor)THEN
       IF(SameFrag)THEN
          TENSOR2%natom(1) = natoms1
          TENSOR2%natom(2) = 0
          TENSOR2%natom(3) = 0
          TENSOR2%natom(4) = 0
          TENSOR2%nLSAO = TENSOR2%natom(1)
          CALL MEM_ALLOC(TENSOR2%LSAO,TENSOR2%natom(1))
          CALL MEM_ALLOC(TENSOR2%INDEX,TENSOR2%natom(1),1,1,4)
          nsize = size(TENSOR2%INDEX,KIND=long)*mem_intsize
          call mem_allocated_mem_lstensor(nsize)

          DO I2=1,4
             DO Iatom2 = 1,TENSOR2%natom(1)
                TENSOR2%INDEX(Iatom2,1,1,I2) = Iatom2
             ENDDO
          ENDDO
          TENSOR2%nbast(1) = 3
          TENSOR2%nbast(2) = TENSOR2%nLSAO
          TENSOR2%nbast(3) = 1
          TENSOR2%nbast(4) = 1
       ELSE
          maxnAtom=MAX(natoms1,natoms2,natoms3,natoms4)
          TENSOR2%nLSAO = natoms1+natoms2+natoms3+natoms4
          TENSOR2%nbast(1) = 3
          TENSOR2%nbast(2) = TENSOR2%nLSAO
          TENSOR2%nbast(3) = 1
          TENSOR2%nbast(4) = 1
          CALL MEM_ALLOC(TENSOR2%LSAO,natoms1+natoms2+natoms3+natoms4)
          CALL MEM_ALLOC(TENSOR2%INDEX,Maxnatom,1,1,4)
          nsize = size(TENSOR2%INDEX,KIND=long)*mem_intsize
          call mem_allocated_mem_lstensor(nsize)
          natomG = 0
          DO iAO=1,4
             DO Iatom2 = 1,natom(iAO)
                TENSOR2%INDEX(Iatom2,1,1,iAO) = natomG + Iatom2                
             ENDDO
             natomG = natomG + natom(iAO)
          ENDDO
          TENSOR2%natom(1) = nAtoms1
          TENSOR2%natom(2) = nAtoms2
          TENSOR2%natom(3) = nAtoms3
          TENSOR2%natom(4) = nAtoms4
       ENDIF
       AllocInt = 0
       AllocIntS = 0
       IF(SameFrag)THEN
          AllocRealk = 3*ndim5*natoms1
       ELSE
          AllocRealk = 3*ndim5*(natom(1)+natom(2)+natom(3)+natom(4))
       ENDIF
       call init_lstensorMem(AllocInt,AllocRealk,AllocInts,lstmem_index)
       call zero_lstensorMem
       TENSOR2%lstmem_index = lstmem_index
       DO Iatom2 = 1,natoms1
          I2=iatom2
          call lsaotensor_nullify(TENSOR2%LSAO(I2))
          lsao2 => TENSOR2%LSAO(I2) 
          lsao2%ATOM = 0
          lsao2%ATOM(1) = Iatom2
          lsao2%AOBATCH = 0
          NULLIFY(lsao2%nOrb)
          NULLIFY(lsao2%startLocalOrb)
          NULLIFY(lsao2%startGlobalOrb)    
          CALL MEM_LSTPOINTER_ALLOC(lsao2%elms,3*ndim5)
          lsao2%nelms = 3
       ENDDO
       natomG = natom(1)
       IF(.NOT.SameFrag)THEN
          DO IAO = 2,4
             DO Iatom2 = 1,natom(iAO)
                I2=iatom2
                call lsaotensor_nullify(TENSOR2%LSAO(natomG+I2))
                lsao2 => TENSOR2%LSAO(natomG+I2) 
                lsao2%ATOM = 0
                lsao2%ATOM(iAO) = Iatom2
                lsao2%AOBATCH = 0
                NULLIFY(lsao2%nAngmom)
                NULLIFY(lsao2%nOrb)
                NULLIFY(lsao2%startLocalOrb)
                NULLIFY(lsao2%startGlobalOrb)
                CALL MEM_LSTPOINTER_ALLOC(lsao2%elms,3*ndim5)
                lsao2%nelms = 3
             ENDDO
             natomG = natomG + natom(iAO)
          ENDDO
       ENDIF
    ELSEIF(TENSOR2%pChargetensor)THEN
       IF(nAtoms1.NE.TENSOR1%natom(1))THEN
          call lsquit('something is wrong thomas. WEEGEDbbfgf',-1)
       ENDIF
       CALL MEM_ALLOC(TENSOR2%LSAO,TENSOR1%natom(1))
       CALL MEM_ALLOC(TENSOR2%INDEX,TENSOR1%natom(1),1,1,1)
       nsize = size(TENSOR2%INDEX,KIND=long)*mem_intsize
       call mem_allocated_mem_lstensor(nsize)
       DO Iatom1 = 1,natoms1
          TENSOR2%INDEX(Iatom1,1,1,1) = Iatom1
       ENDDO
       AllocInt = 2*TENSOR1%natom(1)
       AllocIntS = 0
       AllocRealk = ndim5*TENSOR1%natom(1)
       call init_lstensorMem(AllocInt,AllocRealk,AllocInts,lstmem_index)
       call zero_lstensorMem
       TENSOR2%lstmem_index = lstmem_index
       TENSOR2%natom = TENSOR1%natom
       TENSOR2%nbast = TENSOR1%nbast
       TENSOR2%ndim5 = TENSOR1%ndim5
       DO Iatom1 = 1,TENSOR1%natom(1)
          I1=iatom1
          call lsaotensor_nullify(TENSOR2%LSAO(I1))
          lsao2 => TENSOR2%LSAO(I1) 
          lsao2%ATOM(1) = Iatom1
          lsao2%ATOM(2) = 0
          lsao2%ATOM(3) = 0
          lsao2%ATOM(4) = 0
          lsao2%AOBATCH = 0
          lsao2%AOBATCH(1) = 0
          lsao2%AOBATCH(2) = 0
          lsao2%AOBATCH(3) = 0
          lsao2%AOBATCH(4) = 0
          lsao2%nelms = ndim5
          NULLIFY(lsao2%elms)
          CALL MEM_LSTPOINTER_ALLOC(lsao2%elms,ndim5)
          lsao2%nLocal(1) = 1
          lsao2%nLocal(2) = 1
          lsao2%nLocal(3) = 1
          lsao2%nLocal(4) = 1
          lsao2%maxBat = 1
          lsao2%maxAng = 1
          CALL MEM_LSTPOINTER_ALLOC(lsao2%startLocalOrb,2)
          lsao2%startLocalOrb(1)=1
          lsao2%startLocalOrb(2)=1
       ENDDO
       TENSOR2%nLSAO = natoms1
    ELSEIF(TENSOR2%Screentensor)THEN
       TENSOR2%natom(1) = natoms1
       TENSOR2%natom(2) = natoms2
       TENSOR2%natom(3) = natoms3
       TENSOR2%natom(4) = natoms4
       TENSOR2%nbast(1) = nbast1
       TENSOR2%nbast(2) = nbast2 
       TENSOR2%nbast(3) = nbast3 
       TENSOR2%nbast(4) = nbast4 
       TENSOR2%ndim5   = TENSOR1%ndim5
       ndim5 = TENSOR2%ndim5
       IF(TENSOR1%nSLSAO.GT.0)THEN
          I2=0
          DO Jatom2 = 1,natoms2
             Jatom1 = atoms2(Jatom2)
             DO Iatom2 = 1,natoms1
                Iatom1 = atoms1(Iatom2)
                I1 = TENSOR1%INDEX(IATOM1,JATOM1,1,1)
                IF(I1.NE.0)I2=I2+1
             ENDDO
          ENDDO
          TENSOR2%nSLSAO  = I2
          CALL MEM_ALLOC(TENSOR2%SLSAO,TENSOR2%nSLSAO)
          CALL MEM_ALLOC(TENSOR2%INDEX,TENSOR2%natom(1),TENSOR2%natom(2),1,1)
          nsize = size(TENSOR2%INDEX,KIND=long)*mem_intsize
          call mem_allocated_mem_lstensor(nsize)
          call retrieve_lstMemVal(TENSOR1%lstmem_index,AllocInt,AllocRealk,AllocInts)
          call init_lstensorMem(AllocInt,AllocRealk,AllocInts,lstmem_index)
          call zero_lstensorMem
          TENSOR2%lstmem_index = lstmem_index
          I2=0
          Jbatch = 0
          DO Jatom2 = 1,natoms2
             Jatom1 = atoms2(Jatom2)             
             Ibatch = 0
             DO Iatom2 = 1,natoms1
                Iatom1 = atoms1(Iatom2)
                I1 = TENSOR1%INDEX(IATOM1,JATOM1,1,1)
                IF(I1.NE.0)THEN
                   I2=I2+1
                   TENSOR2%INDEX(IATOM2,JATOM2,1,1) = I2
                   call SLSAOTENSOR_copy(TENSOR1%SLSAO(I1),TENSOR2%SLSAO(I2),6)
                   TENSOR2%SLSAO(I2)%ATOM(1) = IATOM2
                   TENSOR2%SLSAO(I2)%ATOM(2) = JATOM2
                   TENSOR2%SLSAO(I2)%AOBATCH(1) = Ibatch
                   TENSOR2%SLSAO(I2)%AOBATCH(2) = Jbatch                   
                ELSE
                   TENSOR2%INDEX(IATOM2,JATOM2,1,1) = 0
                ENDIF
                Ibatch = Ibatch + TENSOR2%nAObatch(Iatom2,1)
             ENDDO
             Jbatch = Jbatch + TENSOR2%nAObatch(Jatom2,2)             
          ENDDO
       ELSE
          TENSOR2%nSLSAO = 0
       ENDIF
       IF(associated(TENSOR1%maxgab))THEN
          IF(natoms1.EQ.TENSOR1%natom(1).AND.natoms2.EQ.TENSOR1%natom(2))THEN
             call mem_alloc(TENSOR2%maxgab,TENSOR1%nbatches(1),TENSOR1%nbatches(2))
             nsize = size(TENSOR2%maxgab,KIND=long)*mem_shortintsize
             call mem_allocated_mem_lstensor(nsize)
             do jbat = 1,TENSOR1%nbatches(2) 
                do ibat = 1,TENSOR1%nbatches(1) 
                   TENSOR2%maxgab(ibat,jbat) = TENSOR1%maxgab(ibat,jbat)
                enddo
             enddo
             TENSOR2%nbatches(1) = TENSOR1%nbatches(1)
             TENSOR2%nbatches(2) = TENSOR1%nbatches(2)
             TENSOR2%nbatches(3) = TENSOR1%nbatches(3)
             TENSOR2%nbatches(4) = TENSOR1%nbatches(4)
          ELSE
             nbatches1 = 0
             DO Iatom = 1,natoms1
                nbatches1 = nbatches1 + TENSOR2%nAObatch(Iatom,1)
             ENDDO
             TENSOR2%nbatches(1) = nbatches1
             nbatches2 = 0
             DO Iatom = 1,natoms2
                nbatches2 = nbatches2 + TENSOR2%nAObatch(Iatom,2)
             ENDDO
             TENSOR2%nbatches(2) = nbatches2
             TENSOR2%nbatches(3) = 1
             TENSOR2%nbatches(4) = 1
             call mem_alloc(TENSOR2%maxgab,nbatches1,nbatches2)
             nsize = size(TENSOR2%maxgab,KIND=long)*mem_shortintsize
             call mem_allocated_mem_lstensor(nsize)
             call ls_sizero(TENSOR2%maxgab,nbatches1*nbatches2)
             DO Jatom2 = 1,natoms2 
                Jatom1 = atoms2(Jatom2)
                nbatJ = TENSOR2%nAOBATCH(Jatom2,2)
                DO Iatom2 = 1,natoms1
                   Iatom1 = atoms1(Iatom2)
                   nbatI = TENSOR2%nAOBATCH(Iatom2,1)
                   I2 = TENSOR2%INDEX(IATOM2,JATOM2,1,1)
                   I1 = TENSOR1%INDEX(IATOM1,JATOM1,1,1)
                   IF(I2.NE.0)THEN
                      IBATCH2 = TENSOR2%SLSAO(I2)%AOBATCH(1)
                      JBATCH2 = TENSOR2%SLSAO(I2)%AOBATCH(2)
                      IBATCH1 = TENSOR1%SLSAO(I1)%AOBATCH(1)
                      JBATCH1 = TENSOR1%SLSAO(I1)%AOBATCH(2)
                      DO Jbat = 1,nbatJ
                         DO Ibat = 1,nbatI
                            TENSOR2%maxgab(ibatch2+ibat,jbatch2+jbat) = TENSOR1%maxgab(ibatch1+ibat,jbatch1+jbat)
                         enddo
                      enddo
                   ENDIF
                enddo
             enddo
             !             call shortint_output(TENSOR2%maxgab,size(TENSOR2%maxgab,1),&
             !                  & size(TENSOR2%maxgab,2),6)
             !             call shortint_output(TENSOR1%maxgab,size(TENSOR1%maxgab,1),&
             !                  & size(TENSOR1%maxgab,2),6)
          ENDIF
          call set_lst_maxgabelms(TENSOR2)
       ELSE
          nullify(TENSOR2%maxgab)
          TENSOR2%maxgabelm = shortzero
       ENDIF
       IF(associated(TENSOR1%maxprimgab))THEN
          IF(natoms1.EQ.TENSOR1%natom(1).AND.natoms2.EQ.TENSOR1%natom(2))THEN
             call mem_alloc(TENSOR2%maxprimgab,TENSOR1%nbatches(1),TENSOR1%nbatches(2))
             nsize = size(TENSOR2%maxprimgab,KIND=long)*mem_shortintsize
             call mem_allocated_mem_lstensor(nsize)
             do jbat = 1,TENSOR1%nbatches(2) 
                do ibat = 1,TENSOR1%nbatches(1) 
                   TENSOR2%maxprimgab(ibat,jbat) = TENSOR1%maxprimgab(ibat,jbat)
                enddo
             enddo
             TENSOR2%nbatches(1) = TENSOR1%nbatches(1)
             TENSOR2%nbatches(2) = TENSOR1%nbatches(2)
             TENSOR2%nbatches(3) = TENSOR1%nbatches(3)
             TENSOR2%nbatches(4) = TENSOR1%nbatches(4)
          ELSE
             nbatches1 = 0
             DO Iatom = 1,natoms1
                nbatches1 = nbatches1 + TENSOR2%nAObatch(Iatom,1)
             ENDDO
             TENSOR2%nbatches(1) = nbatches1
             nbatches2 = 0
             DO Iatom = 1,natoms2
                nbatches2 = nbatches2 + TENSOR2%nAObatch(Iatom,2)
             ENDDO
             TENSOR2%nbatches(2) = nbatches2
             TENSOR2%nbatches(3) = 1
             TENSOR2%nbatches(4) = 1
             nullify(TENSOR2%maxprimgab)
             call mem_alloc(TENSOR2%maxprimgab,nbatches1,nbatches2)
             nsize = size(TENSOR2%maxprimgab,KIND=long)*mem_shortintsize
             call mem_allocated_mem_lstensor(nsize)
             call ls_sizero(TENSOR2%maxprimgab,nbatches1*nbatches2)
             DO Jatom2 = 1,natoms2 
                Jatom1 = atoms2(Jatom2)
                nbatJ = TENSOR2%nAOBATCH(Jatom2,2)
                DO Iatom2 = 1,natoms1
                   Iatom1 = atoms1(Iatom2)
                   nbatI = TENSOR2%nAOBATCH(Iatom2,1)
                   I2 = TENSOR2%INDEX(IATOM2,JATOM2,1,1)
                   I1 = TENSOR1%INDEX(IATOM1,JATOM1,1,1)
                   IF(I2.NE.0)THEN
                      IBATCH2 = TENSOR2%SLSAO(I2)%AOBATCH(1)
                      JBATCH2 = TENSOR2%SLSAO(I2)%AOBATCH(2)
                      IBATCH1 = TENSOR1%SLSAO(I1)%AOBATCH(1)
                      JBATCH1 = TENSOR1%SLSAO(I1)%AOBATCH(2)
                      DO Jbat = 1,nbatJ
                         DO Ibat = 1,nbatI
                            TENSOR2%maxprimgab(ibatch2+ibat,jbatch2+jbat) = TENSOR1%maxprimgab(ibatch1+ibat,jbatch1+jbat)
                         enddo
                      enddo
                   ENDIF
                enddo
             enddo
!             call shortint_output(TENSOR2%maxprimgab,size(TENSOR2%maxprimgab,1),&
!                  & size(TENSOR2%maxprimgab,2),6)
!             call shortint_output(TENSOR1%maxprimgab,size(TENSOR1%maxprimgab,1),&
!                  & size(TENSOR1%maxprimgab,2),6)
          ENDIF
          call set_lst_maxprimgabelms(TENSOR2)
       ELSE
          nullify(TENSOR2%maxprimgab)
          TENSOR2%maxprimgabelm = shortzero
       ENDIF
       TENSOR2%nMBIE = TENSOR1%nMBIE
       IF(associated(TENSOR1%MBIE))THEN
          IF(TENSOR2%nMBIE.NE.2)CALL LSQUIT('nMBIE not 2 in sublstensor',-1)
          IF(natoms1.EQ.TENSOR1%natom(1).AND.natoms2.EQ.TENSOR1%natom(2))THEN
             nullify(TENSOR2%MBIE)
             call mem_alloc(TENSOR2%MBIE,TENSOR1%nMBIE,TENSOR1%nbatches(1),TENSOR1%nbatches(2))
             nsize = size(TENSOR2%MBIE,KIND=long)*mem_realsize
             call mem_allocated_mem_lstensor(nsize)
             do jbat = 1,TENSOR1%nbatches(2) 
                do ibat = 1,TENSOR1%nbatches(1) 
                   TENSOR2%MBIE(1,ibat,jbat) = TENSOR1%MBIE(1,ibat,jbat)
                   TENSOR2%MBIE(2,ibat,jbat) = TENSOR1%MBIE(2,ibat,jbat)
                enddo
             enddo
          ELSE
             nbatches1 = 0
             DO Iatom = 1,natoms1
                nbatches1 = nbatches1 + TENSOR2%nAObatch(Iatom,1)
             ENDDO
             TENSOR2%nbatches(1) = nbatches1
             nbatches2 = 0
             DO Iatom = 1,natoms2
                nbatches2 = nbatches2 + TENSOR2%nAObatch(Iatom,2)
             ENDDO
             TENSOR2%nbatches(2) = nbatches2
             TENSOR2%nbatches(3) = 1
             TENSOR2%nbatches(4) = 1
             nullify(TENSOR2%MBIE)
             call mem_alloc(TENSOR2%MBIE,2,nbatches1,nbatches2)
             nsize = size(TENSOR2%MBIE,KIND=long)*mem_realsize
             call mem_allocated_mem_lstensor(nsize)
             call ls_dzero(TENSOR2%MBIE,2*nbatches1*nbatches2)
             DO Jatom2 = 1,natoms2 
                Jatom1 = atoms2(Jatom2)
                nbatJ = TENSOR2%nAOBATCH(Jatom2,2)
                DO Iatom2 = 1,natoms1
                   Iatom1 = atoms1(Iatom2)
                   nbatI = TENSOR2%nAOBATCH(Iatom2,1)
                   I2 = TENSOR2%INDEX(IATOM2,JATOM2,1,1)
                   I1 = TENSOR1%INDEX(IATOM1,JATOM1,1,1)
                   IF(I2.NE.0)THEN
                      IBATCH2 = TENSOR2%SLSAO(I2)%AOBATCH(1)
                      JBATCH2 = TENSOR2%SLSAO(I2)%AOBATCH(2)
                      IBATCH1 = TENSOR1%SLSAO(I1)%AOBATCH(1)
                      JBATCH1 = TENSOR1%SLSAO(I1)%AOBATCH(2)
                      DO Jbat = 1,nbatJ
                         DO Ibat = 1,nbatI
                            TENSOR2%MBIE(1,ibatch2+ibat,jbatch2+jbat) = TENSOR1%MBIE(1,ibatch1+ibat,jbatch1+jbat)
                            TENSOR2%MBIE(2,ibatch2+ibat,jbatch2+jbat) = TENSOR1%MBIE(2,ibatch1+ibat,jbatch1+jbat)
                         enddo
                      enddo
                   ENDIF
                enddo
             enddo
          ENDIF
       ELSE
          nullify(TENSOR2%MBIE)
       ENDIF
    ELSEIF(TENSOR2%Econtrib)THEN
       call Init_lstensor_1dim(TENSOR2,ndim5,6)
    ELSE

       TENSOR2%natom(1) = natoms1
       TENSOR2%natom(2) = natoms2
       TENSOR2%natom(3) = natoms3
       TENSOR2%natom(4) = natoms4
       TENSOR2%nbast(1) = nbast1
       TENSOR2%nbast(2) = nbast2 
       TENSOR2%nbast(3) = nbast3 
       TENSOR2%nbast(4) = nbast4 
       TENSOR2%ndim5   = TENSOR1%ndim5

       DO iAO=1,4
         TENSOR2%nbatches(iAO) = 0
         DO Iatom = 1,TENSOR2%natom(iAO)
            TENSOR2%nbatches(iAO) = TENSOR2%nbatches(iAO) + TENSOR2%nAObatch(Iatom,iAO)
         ENDDO
       ENDDO

       ndim5 = TENSOR2%ndim5
       TENSOR2%nLSAO  = natoms1*natoms2*natoms3*natoms4
       CALL MEM_ALLOC(TENSOR2%LSAO,TENSOR2%nLSAO)
       CALL MEM_ALLOC(TENSOR2%INDEX,natoms1,natoms2,natoms3,natoms4)
       nsize = size(TENSOR2%INDEX,KIND=long)*mem_intsize
       call mem_allocated_mem_lstensor(nsize)
       !FIXME OVERKILL TOO MUCH MEMORY FOR TENSOR2
       call retrieve_lstMemVal(TENSOR1%lstmem_index,AllocInt,AllocRealk,AllocInts)
       call init_lstensorMem(AllocInt,AllocRealk,AllocInts,lstmem_index)
       call zero_lstensorMem
       TENSOR2%lstmem_index = lstmem_index
       I2=0
       Lbatch = 0
       DO Latom2 = 1,natoms4
          Latom1 = atoms4(Latom2)
          Kbatch = 0
          DO Katom2 = 1,natoms3
             Katom1 = atoms3(Katom2)
             Jbatch = 0
             DO Jatom2 = 1,natoms2
                Jatom1 = atoms2(Jatom2)
                Ibatch = 0
                DO Iatom2 = 1,natoms1
                   Iatom1 = atoms1(Iatom2)
                   I1 = TENSOR1%INDEX(IATOM1,JATOM1,KATOM1,LATOM1)
                   IF(I1.NE.0)THEN
                      I2=I2+1
                      TENSOR2%INDEX(IATOM2,JATOM2,KATOM2,Latom2) = I2
                      call LSAOTENSOR_copy(TENSOR1%LSAO(I1),TENSOR2%LSAO(I2),6)
                      TENSOR2%LSAO(I2)%ATOM(1) = IATOM2
                      TENSOR2%LSAO(I2)%ATOM(2) = JATOM2
                      TENSOR2%LSAO(I2)%ATOM(3) = KATOM2
                      TENSOR2%LSAO(I2)%ATOM(4) = LATOM2
                      TENSOR2%LSAO(I2)%AOBATCH(1) = Ibatch
                      TENSOR2%LSAO(I2)%AOBATCH(2) = Jbatch                   
                      TENSOR2%LSAO(I2)%AOBATCH(3) = Kbatch
                      TENSOR2%LSAO(I2)%AOBATCH(4) = Lbatch                   
                   ELSE
                      TENSOR2%INDEX(IATOM2,JATOM2,KATOM2,Latom2) = 0
                   ENDIF
                   Ibatch = Ibatch + TENSOR2%nAObatch(Iatom2,1)             
                enddo
                Jbatch = Jbatch + TENSOR2%nAObatch(Jatom2,2)             
             enddo
             Kbatch = Kbatch + TENSOR2%nAObatch(Katom2,3)             
          enddo
          Lbatch = Lbatch + TENSOR2%nAObatch(Latom2,4)             
       enddo
       TENSOR2%nLSAO  = I2
    endif
  end SUBROUTINE build_sublstensor_from_full_lstensor

  !> \brief 
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param TENSOR1 dummy full lstensor
  !> \param TENSOR2 the new sublstensor
  SUBROUTINE memdist_lstensor_BuildFromScalapack(TENSOR,comm,mynum,numnodes,mat)
    implicit none
    TYPE(LSTENSOR)     :: TENSOR
    integer(kind=ls_mpik) :: comm,mynum,numnodes
    TYPE(MATRIX),optional,target  :: MAT
    !
    integer(kind=ls_mpik) :: sender,reciever
    TYPE(MATRIX),pointer  :: MAT2
    INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,Iangmom,Jangmom,Kangmom,Langmom,nOrbA
    INTEGER    :: sA,nOrbB,sB,nOrbC,sC,J,DIMI,DIMJ,SUM1,SUM2
    INTEGER    :: nOrbD,sD,IORB,JORB,KORB,LORB,IELM,IMAT
    INTEGER    :: nbast(4),MATINDEX(4),nrow2,ncol2,bufferproc
    INTEGER    :: IATOM,JATOM,KATOM,LATOM,nbat1,nbat2,nbat3,nbat4
    integer :: IATOMGL,JATOMGL,gi
    TYPE(LSAOTENSOR),pointer    :: lsao
    TYPE(GLOBALLSAOTENSOR),pointer    :: g_lsao
    REAL(REALK),pointer     :: elms(:)
#ifdef VAR_SCALAPACK
    logical :: CommunicationNodes(0:numnodes-1),debug
    integer :: localnrow(0:numnodes-1),localncol(0:numnodes-1)
    integer :: nsizeSend(0:numnodes-1),nsizeRecv(0:numnodes-1)
    integer :: bufferOffset(0:numnodes-1)
    type(lsmatrix) :: ASendbuffer(0:numnodes-1)
    type(lsmatrix) :: ARecvbuffer(0:numnodes-1)
    integer :: nbuffer,iproc,ndim1,ndim2,iproc2,nodes,maxAng,maxBat
    logical :: CommunicationNodesGlobalSend(0:numnodes-1)
    logical :: CommunicationNodesGlobalRecv(0:numnodes-1),LowerDiagZero
    real(realk),pointer :: dummy(:)
    integer :: IELMLSTENSOR_ORDERING,offset2,offset1,IELMMAT,sAtmp,sBtmp
    LowerDiagZero = TENSOR%LowerDiagZero
    nullify(dummy)
    call mem_alloc(dummy,TENSOR%G_maxnLocal*TENSOR%G_maxnLocal)
    nullify(MAT2)
    if(mynum.EQ.0)then
       IF(.NOT.present(mat))CALL LSQUIT('Master error in build_memdist_lstensor_fromScalapack',-1)
       call PDM_MATRIXSYNC(MAT)
       MAT2 => MAT
    else
       allocate(MAT2)
       call PDM_MATRIXSYNC(MAT2)
    endif
    bufferOffset = 0
    debug = .FALSE.
    nbuffer = numnodes-1
    nodes = numnodes-1
!    print*,'BUILD LSTENSOR the mat_scalapack_print_local(MAT2)',infpar%mynum
!    print*,'DIM:',MAT2%nrow,MAT2%ncol,'mynum',infpar%mynum
!    call sleep(infpar%mynum*10)
!    call mat_scalapack_print_local(MAT2)
    
    !====================================================================================
    !                                Step 1   
    !====================================================================================
    ! Determine who I (infpar%mynum) should send a block to. Loop through the full lstensor and if the 
    ! iproc that owns that lstensor block (lsao%mynum) needs something from me 
    ! (mynum=infpar%mynum) I collect and send this should maybe avoided - to avoid referencing
    !====================================================================================
    nsizeSend = 0
    bufferOffset = 0
    bufferproc = 0 
    CommunicationNodesGlobalSend=.FALSE.
    !THIS COULD BE DONE USING LOWER DIAGONAL MATRIX -> DMAT SPECIAL CASE TO REDUCE COMMUNICATION
    DO I = 1,TENSOR%G_nLSAO
       G_lsao => TENSOR%G_LSAO(I)
       elms => dummy
       ndim1 = G_lsao%G_nLocal(1)
       ndim2 = G_lsao%G_nLocal(2)
       sA = G_lsao%G_startGlobalOrb(1) 
       sB = G_lsao%G_startGlobalOrb(2) 
       !determine CommunicationNodes,localnrow,localncol
       call scalapack_transform_aux(MAT2,elms,ndim1,ndim2,ndim1,ndim2,sA,sB,&
            & CommunicationNodes,localnrow,localncol,ASendbuffer,nbuffer,bufferOffset,debug,bufferproc,nodes,&
            & STA_null)
       !if g_lsao%mynum needs a block from me (mynum) then I need to send it a block       
       IF(CommunicationNodes(mynum))THEN
          CommunicationNodesGlobalSend(g_lsao%mynum) = .TRUE.
          nsizeSend(g_lsao%mynum) = nsizeSend(g_lsao%mynum) + localnrow(mynum)*localncol(mynum)
       ENDIF
    ENDDO
    
    !call mem_alloc buffers for each processor I need to send a block to, including myself
    do iproc = 0,numnodes-1
       IF(CommunicationNodesGlobalSend(iproc))THEN
          CALL mem_alloc(ASendbuffer(Iproc)%elms,nsizeSend(iproc))
       ENDIF
    enddo

    !build buffers that I need to send to others (and myself)
    bufferOffset = 0
    DO I = 1,TENSOR%G_nLSAO
       G_lsao => TENSOR%G_LSAO(I)
       !g_lsao%mynum is the node that owns this lstensor block and I only need to 
       !build up a buffer for the nodes that need something from me.
       !if CommunicationNodesGlobalSend(g_lsao%mynum) is false I have previously determined
       !that I do not need to send anything to this node.
       IF(CommunicationNodesGlobalSend(g_lsao%mynum))THEN
          elms => dummy
          ndim1 = G_lsao%G_nLocal(1)
          ndim2 = G_lsao%G_nLocal(2)
          sA = G_lsao%G_startGlobalOrb(1) 
          sB = G_lsao%G_startGlobalOrb(2) 
          !build up the buffer Asendbuffer(G_lsao%mynum)%elms that I need to send to G_lsao%mynum  
          bufferproc = G_lsao%mynum  
          !determine localnrow,localncol
          call scalapack_transform_aux(MAT2,elms,ndim1,ndim2,ndim1,ndim2,sA,sB,&
               & CommunicationNodes,localnrow,localncol,ASendbuffer,&
               & nbuffer,bufferOffset,debug,bufferproc,nodes,STA_null)
          !localnrow is the first dimension of the block(s) residing on this processor
          !and therefor the rowsize of the buffer that should be sent to G_lsao%mynum 
          localnrow(bufferproc) = localnrow(mynum)
          !Build Buffer(g_lsao%mynum) from my local part of matrix MAT, that I will send to g_lsao%mynum
          call scalapack_transform_aux(MAT2,elms,ndim1,ndim2,ndim1,ndim2,sA,sB,&
               & CommunicationNodes,localnrow,localncol,ASendbuffer,&
               & nbuffer,bufferOffset,debug,bufferproc,nodes,STA_BuildBufferFromAlocal)
          !if g_lsao%mynum needed a block from me I update the bufferoffset
          IF(CommunicationNodes(mynum))THEN
             bufferOffset(bufferproc) = bufferOffset(bufferproc) + localnrow(mynum)*localncol(mynum)
          ENDIF
       ENDIF
    ENDDO
!    do iproc = 0,numnodes-1
!       Print*,'I am ',infpar%mynum,' should I send to ',Iproc,'?',CommunicationNodesGlobalSend(iproc)
!       IF(CommunicationNodesGlobalSend(iproc))THEN
!          Print*,'I am ',infpar%mynum,' and I send this block to ',Iproc
!          CALL ls_output(ASendbuffer(Iproc)%elms,1,nsizesend(iProc),1,1,nsizesend(iProc),1,1,6)
!       ENDIF
!    enddo

    !====================================================================================
    !                                Step 2   
    !====================================================================================
    ! Determine who I (infpar%mynum) should recieve from in order to get my part of the matrix
    ! Loop through the sub lstensor, determine which iprocs own parts of the matrix I need 
    ! determine sizes and call mem_alloc - ready to recieve 
    !====================================================================================

    !loop through the part of the lstensor I need
    nsizeRecv=0
    CommunicationNodesGlobalRecv = .FALSE.
    bufferOffset = 0
    bufferproc = 0
    DO I = 1,TENSOR%nLSAO
       lsao => TENSOR%LSAO(I)
       elms => lsao%elms
       ndim1 = lsao%nLocal(1)
       ndim2 = lsao%nLocal(2)
       IATOM = TENSOR%LSAO(I)%ATOM(1)
       JATOM = TENSOR%LSAO(I)%ATOM(2)
       IATOMGL = TENSOR%G_fullatoms1(IATOM)
       JATOMGL = TENSOR%G_fullatoms2(JATOM)
       gI = TENSOR%G_INDEX(IATOMGL,JATOMGL)
       sA = TENSOR%G_LSAO(gI)%G_startGlobalOrb(1) 
       sB = TENSOR%G_LSAO(gI)%G_startGlobalOrb(2) 
       !determine CommunicationNodes,localnrow,localncol
       call scalapack_transform_aux(MAT2,elms,ndim1,ndim2,ndim1,ndim2,sA,sB,&
            & CommunicationNodes,localnrow,localncol,ARecvbuffer,nbuffer,bufferOffset,debug,bufferproc,nodes,&
            & STA_null)
       !for all procs that have a matrix block I need to build elms I need to determine how much space I need to 
       !recieve this matrix block
       do iproc = 0,numnodes-1
          nsizeRecv(iproc) = nsizeRecv(iproc) + localnrow(iproc)*localncol(iproc)
          IF(CommunicationNodes(iproc)) CommunicationNodesGlobalRecv(iproc) = .TRUE.
       enddo
    ENDDO

    !call mem_alloc buffers one for each processor I need to recieve a block from 
    do iproc = 0,numnodes-1
       IF(CommunicationNodesGlobalRecv(iproc))THEN
          CALL mem_alloc(ARecvbuffer(Iproc)%elms,nsizeRecv(iproc))
       ENDIF
    enddo

    !send and recieve 
    IF(CommunicationNodesGlobalRecv(mynum))THEN
       !need to send myself a packages so I copy from send buffer to recieve buffer
       ARecvbuffer(mynum)%elms = ASendbuffer(mynum)%elms
       CALL DCOPY(size(ASendbuffer(mynum)%elms),ASendbuffer(mynum)%elms,1,ARecvbuffer(mynum)%elms,1)
       CALL mem_dealloc(ASendbuffer(mynum)%elms)
       CommunicationNodesGlobalSend(mynum) = .FALSE. 
    ENDIF
    DO Iproc=0,numnodes-1
       !send to or recieve from Iproc
       if(mynum.EQ.Iproc)THEN 
          !recieve from others Iproc2
          DO Iproc2=0,numnodes-1
             IF(Iproc2.EQ.mynum)CYCLE
             IF(CommunicationNodesGlobalRecv(Iproc2))THEN
!                call ls_mpisendrecv(ARecvbuffer(Iproc2)%elms,nsizeRecv(Iproc2),comm,Iproc2,mynum) !sender=Iproc2 reciever=mynum
                sender = Iproc2
                call ls_mpisendrecv(ARecvbuffer(Iproc2)%elms,nsizeRecv(Iproc2),comm,sender,mynum) !sender=Iproc2 reciever=mynum
             ENDIF
          ENDDO
       else
          !send to Iproc
          IF(CommunicationNodesGlobalSend(Iproc))THEN
!             call ls_mpisendrecv(ASendbuffer(Iproc)%elms,nsizeSend(Iproc),comm,mynum,Iproc) !sender=mynum reciever=Iproc
             reciever = Iproc
             call ls_mpisendrecv(ASendbuffer(Iproc)%elms,nsizeSend(Iproc),comm,mynum,reciever) !sender=mynum reciever=Iproc
          ENDIF
       endif
    ENDDO

    !call mem_dealloc buffer used for sending
    do iproc = 0,numnodes-1
       IF(CommunicationNodesGlobalSend(iproc))THEN
          CALL mem_dealloc(ASendbuffer(Iproc)%elms)
       ENDIF
    enddo

!    IF(mynum.EQ.0)THEN
!       do iproc = 0,numnodes-1
!          Print*,'I am ',infpar%mynum,'did I recieved from ',Iproc,'?',CommunicationNodesGlobalRecv(iproc)
!          IF(CommunicationNodesGlobalRecv(iproc))THEN
!             Print*,'I am ',infpar%mynum,'and I recieved this block from ',Iproc,'size',nsizeRecv(iProc)
!             CALL ls_output(ARecvbuffer(Iproc)%elms,1,nsizeRecv(iProc),1,1,nsizeRecv(iProc),1,1,6)
!          ENDIF
!       enddo
!    ENDIF

    !fill in lstensor from buffers
    bufferOffset = 0
    nbuffer = numnodes-1
    DO I = 1,TENSOR%nLSAO
       lsao => TENSOR%LSAO(I)
       elms => lsao%elms
       ndim1 = lsao%nLocal(1)
       ndim2 = lsao%nLocal(2)
       IATOM = TENSOR%LSAO(I)%ATOM(1)
       JATOM = TENSOR%LSAO(I)%ATOM(2)
       IATOMGL = TENSOR%G_fullatoms1(IATOM)
       JATOMGL = TENSOR%G_fullatoms2(JATOM)
       gI = TENSOR%G_INDEX(IATOMGL,JATOMGL)
       sA = TENSOR%G_LSAO(gI)%G_startGlobalOrb(1) 
       sB = TENSOR%G_LSAO(gI)%G_startGlobalOrb(2) 
       bufferproc = mynum    
       !determine localnrow,localncol
       call scalapack_transform_aux(MAT2,elms,ndim1,ndim2,ndim1,ndim2,sA,sB,&
            & CommunicationNodes,localnrow,localncol,ARecvbuffer,nbuffer,&
            & bufferOffset,debug,bufferproc,nodes,&
            & STA_Null)

       !In this case there may be many nodes that contribute to one elms 
       !and even 1 node that contributes with several matrix blocks so we update the 
       !buffer offset inside the subroutine

       call scalapack_transform_aux(MAT2,elms,ndim1,ndim2,ndim1,ndim2,sA,sB,&
            & CommunicationNodes,localnrow,localncol,ARecvbuffer,nbuffer,&
            & bufferOffset,debug,bufferproc,nodes,&
            & STA_BuildFullFromAllBuffer)

       do iproc = 0,numnodes-1
          bufferOffset(iproc) = bufferOffset(iproc) + localnrow(iproc)*localncol(iproc)
       enddo
       maxAng = lsao%maxAng
       IF(maxAng.GT.1)THEN
          !In the case of family type basis set we need to reorder the TENSOR%LSAO(I)%elms 
          !to correspond to the family basisset orientation. 
          IF(LowerDiagZero.AND.(sA.EQ.sB))THEN
             !the Matrix have had its lower tridiagonal part set to zero
             !when we reorder this could/will affect in a Matrix which nolonger have that 
             !property. This is only a problem on the diagonal so in this case we set
             !the lower diagnal part to the upper diagonal part and 
             !after the reordering we set it to zero again
             call FULL_SYMelms_FROM_TRIANGULARelms2(elms,ndim1,ndim2,1)             
          ENDIF
          maxBat = lsao%maxBat
          nbat1 = TENSOR%nAOBATCH(IATOM,1)
          nbat2 = TENSOR%nAOBATCH(JATOM,2)
          IELMLSTENSOR_ORDERING=0
          DO Jbat = 1,nbat2
           offset2 = (Jbat-1)*maxAng+maxBat*maxAng
           DO Jangmom = 1,lsao%nAngmom(Jbat+maxBat)
            nOrbB = lsao%nOrb(Jangmom+offset2)
            sBtmp = lsao%startGlobalOrb(Jangmom+offset2)-lsao%startGlobalOrb(1+maxBat*maxAng)
            DO JORB = 1,nOrbB
             DO Ibat = 1,nbat1
              offset1 = (Ibat-1)*maxAng
              DO Iangmom = 1,lsao%nAngmom(Ibat)
               nOrbA = lsao%nOrb(Iangmom+offset1)
               sAtmp = lsao%startGlobalOrb(Iangmom+offset1)-lsao%startGlobalOrb(1)
               DO IORB = 1,nOrbA
                IELMLSTENSOR_ORDERING = IELMLSTENSOR_ORDERING+1
                IELMMAT = sAtmp+IORB+(sBtmp+JORB-1)*ndim1
                DUMMY(IELMLSTENSOR_ORDERING) = elms(IELMMAT)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
          CALL DCOPY(nDim1*nDim2,dummy,1,elms,1)                    
          IF(LowerDiagZero.AND.(sA.EQ.sB))THEN
             call set_lowerelms_triangular_zero(elms,ndim1,ndim2,1)
          ENDIF
       ENDIF
    ENDDO
    
    !call mem_dealloc the buffers
    do iproc = 0,numnodes-1
       IF(CommunicationNodesGlobalRecv(iproc))THEN
          CALL mem_dealloc(ARecvbuffer(Iproc)%elms)
       ENDIF
    enddo    
    call mem_dealloc(dummy)
#else
      call lsquit('memdist_lstensor_BuildFromScalapack',-1)
#endif


    end SUBROUTINE memdist_lstensor_BuildFromScalapack

  !> \brief copy an lstensor to a new lstensor
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param TENSOR1 dummy full lstensor
  !> \param TENSOR2 the new sublstensor
  SUBROUTINE memdist_lstensor_BuildToScalapack(TENSOR,comm,mynum,numnodes,mat)
    implicit none
    TYPE(LSTENSOR)     :: TENSOR
    integer(kind=ls_mpik) :: comm,mynum,numnodes
    TYPE(MATRIX),optional,target  :: MAT
    !
    integer(kind=ls_mpik) :: sender,reciever
    TYPE(MATRIX),pointer  :: MAT2
    INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,Iangmom,Jangmom,Kangmom,Langmom,nOrbA
    INTEGER    :: sA,nOrbB,sB,nOrbC,sC,J,DIMI,DIMJ,SUM1,SUM2
    INTEGER    :: nOrbD,sD,IORB,JORB,KORB,LORB,IELM,IMAT
    INTEGER    :: nbast(4),MATINDEX(4),nrow2,ncol2,bufferproc
    INTEGER    :: IATOM,JATOM,KATOM,LATOM,nbat1,nbat2,nbat3,nbat4
    TYPE(LSAOTENSOR),pointer    :: lsao
    TYPE(GLOBALLSAOTENSOR),pointer    :: g_lsao
    REAL(REALK),pointer     :: elms(:)
    integer :: IATOMGL,JATOMGL,gi
#ifdef VAR_SCALAPACK
    logical :: CommunicationNodes(0:numnodes-1),debug
    integer :: localnrow(0:numnodes-1)
    integer :: localncol(0:numnodes-1)
    integer :: nsizeSend(0:numnodes-1)
    integer :: nsizeRecv(0:numnodes-1)
    integer :: bufferOffset(0:numnodes-1)
    type(lsmatrix) :: ASendbuffer(0:numnodes-1)
    type(lsmatrix) :: ARecvbuffer(0:numnodes-1)
    integer :: nbuffer,iproc,ndim1,ndim2,iproc2,nodes
    logical :: CommunicationNodesGlobalSend(0:numnodes-1),PermuteResultTensor
    logical :: CommunicationNodesGlobalRecv(0:numnodes-1)
    real(realk),pointer :: dummy(:)
    integer :: IELMLSTENSOR_ORDERING,offset2,offset1,IELMMAT,maxAng,maxBat,sAtmp,sBtmp
    PermuteResultTensor = TENSOR%PermuteResultTensor
!    call lsmpi_barrier(comm)
!    call sleep(infpar%mynum*20+1)
!    print*,'BUILDSCALAPACK: inside memdist_lstensor_BuildToScalapack',mynum
!    call lstensor_print(TENSOR,6)
!    call lsmpi_barrier(comm)

    nullify(dummy)
    call mem_alloc(dummy,TENSOR%G_maxnLocal*TENSOR%G_maxnLocal)
    nullify(MAT2)
    if(mynum.EQ.0)then
       IF(.NOT.present(mat))CALL LSQUIT('Master error in build_memdist_lstensor_fromScalapack',-1)
       call PDM_MATRIXSYNC(MAT)
       MAT2 => MAT
    else
       allocate(MAT2)
       call PDM_MATRIXSYNC(MAT2)
    endif
    bufferOffset = 0
    debug = .FALSE.
    nbuffer = numnodes-1
    nodes = numnodes-1

    !====================================================================================
    !                                Step 1   
    !====================================================================================
    ! Determine who I (infpar%mynum) should send to.
    ! Loop through the sub lstensor I have, determine which iprocs own parts of the 
    ! scalapack matrix I have determine sizes and call mem_alloc - ready to send
    !====================================================================================

    !loop through the part of the lstensor I have
    nsizeSend=0
    CommunicationNodesGlobalSend = .FALSE.
    bufferOffset = 0
    bufferproc = 0
    DO I = 1,TENSOR%nLSAO
       lsao => TENSOR%LSAO(I)
       elms => lsao%elms
       ndim1 = lsao%nLocal(1)
       ndim2 = lsao%nLocal(2)
       IATOM = TENSOR%LSAO(I)%ATOM(1)
       JATOM = TENSOR%LSAO(I)%ATOM(2)
       IATOMGL = TENSOR%G_fullatoms1(IATOM) !global atom index
       JATOMGL = TENSOR%G_fullatoms2(JATOM)
       gI = TENSOR%G_INDEX(IATOMGL,JATOMGL)  
       sA = TENSOR%G_LSAO(gI)%G_startGlobalOrb(1) 
       sB = TENSOR%G_LSAO(gI)%G_startGlobalOrb(2)
       !determine CommunicationNodes,localnrow,localncol
       call scalapack_transform_aux(MAT2,elms,ndim1,ndim2,ndim1,ndim2,sA,sB,&
            & CommunicationNodes,localnrow,localncol,ASendbuffer,nbuffer,bufferOffset,debug,bufferproc,nodes,&
            & STA_null)
       !the lstensor block I have, can have contributions to many matrix blocks belonging to other nodes
       !so if my lstensor block have contribution to a matrix block on node iproc I determine size and set 
       !CommunicationNodesGlobalSend(iproc) = .TRUE.
       do iproc = 0,numnodes-1
          nsizeSend(iproc) = nsizeSend(iproc) + localnrow(iproc)*localncol(iproc)
          IF(CommunicationNodes(iproc)) CommunicationNodesGlobalSend(iproc) = .TRUE.
       enddo
    ENDDO

    !call mem_alloc buffers one for each processor I need to send a block to
    do iproc = 0,numnodes-1
       IF(CommunicationNodesGlobalSend(iproc))THEN
 !         Print*,'BUILDSCALAPACK: I am ',infpar%mynum,' need to send  ',Iproc,' a full block of size',nsizeSend(iproc)
          CALL mem_alloc(ASendbuffer(Iproc)%elms,nsizeSend(iproc))
       ENDIF
    enddo

    !fill in buffers from lstensor
!    print*,'fill in the buffers from my lstensor'
    do iproc = 0,numnodes-1
       bufferOffset(iproc) = 0
    enddo
    nbuffer = numnodes-1
    DO I = 1,TENSOR%nLSAO
       lsao => TENSOR%LSAO(I)
       elms => lsao%elms
       ndim1 = lsao%nLocal(1)
       ndim2 = lsao%nLocal(2)
       IATOM = TENSOR%LSAO(I)%ATOM(1)
       JATOM = TENSOR%LSAO(I)%ATOM(2)
       IATOMGL = TENSOR%G_fullatoms1(IATOM)
       JATOMGL = TENSOR%G_fullatoms2(JATOM)
       gI = TENSOR%G_INDEX(IATOMGL,JATOMGL)
       bufferproc = TENSOR%G_LSAO(gI)%mynum
       sA = TENSOR%G_LSAO(gI)%G_startGlobalOrb(1) 
       sB = TENSOR%G_LSAO(gI)%G_startGlobalOrb(2)

       maxAng = lsao%maxAng
       IF(maxAng.GT.1)THEN
          !If this is a family basis set the lstensor and the matrix
          !is slightly different ordered.
          IF(PermuteResultTensor.AND.(sA.EQ.sB))THEN
             !the lower tridiagonal part of the Matrix have not been calculated 
             !and will be obtained from the upperdiagonal, but when we reorder we may 
             !put elements from the upper triangular part to the lower diagonal part 
             !and thereby overwrite some elements 
             call FULL_SYMelms_FROM_TRIANGULARelms2(elms,ndim1,ndim2,1)             
          ENDIF
          IELMLSTENSOR_ORDERING=0
          maxBat = lsao%maxBat
          nbat1 = TENSOR%nAOBATCH(IATOM,1)
          nbat2 = TENSOR%nAOBATCH(JATOM,2)
          DO Jbat = 1,nbat2
           offset2 = (Jbat-1)*maxAng+maxBat*maxAng
           DO Jangmom = 1,lsao%nAngmom(Jbat+maxBat)
            nOrbB = lsao%nOrb(Jangmom+offset2)
            sBtmp = lsao%startGlobalOrb(Jangmom+offset2)-lsao%startGlobalOrb(1+maxBat*maxAng)
            DO JORB = 1,nOrbB
             MATINDEX(2)=sB+JORB
             DO Ibat = 1,nbat1
              offset1 = (Ibat-1)*maxAng
              DO Iangmom = 1,lsao%nAngmom(Ibat)
               nOrbA = lsao%nOrb(Iangmom+offset1)
               sAtmp = lsao%startGlobalOrb(Iangmom+offset1)-lsao%startGlobalOrb(1)
               DO IORB = 1,nOrbA
                IELMLSTENSOR_ORDERING = IELMLSTENSOR_ORDERING+1
                IELMMAT = sAtmp+IORB+(sBtmp+JORB-1)*ndim1
                DUMMY(IELMMAT) = elms(IELMLSTENSOR_ORDERING)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
          elms => dummy
       ENDIF
       !determine localnrow,localncol
       !print*,'BUILDSCALAPACK I=',I,'mynum',mynum,'elms',elms
       call scalapack_transform_aux(MAT2,elms,ndim1,ndim2,ndim1,ndim2,sA,sB,&
            & CommunicationNodes,localnrow,localncol,ASendbuffer,nbuffer,bufferOffset,debug,bufferproc,nodes,&
            & STA_Null)
       !build up buffers from lstensor block. the lstensor block I have, can have contributions to 
       !many matrix blocks belonging to other nodes and even several matrix blocks belonging to another node
       !I therefor need to update the buffer (bufferOffset is updated inside STA_BuildAllBufferFromFull)
       call scalapack_transform_aux(MAT2,elms,ndim1,ndim2,ndim1,ndim2,sA,sB,&
            & CommunicationNodes,localnrow,localncol,ASendbuffer,nbuffer,bufferOffset,debug,bufferproc,nodes,&
            & STA_BuildAllBufferFromFull)
       
       do iproc = 0,numnodes-1
          bufferOffset(iproc) = bufferOffset(iproc) + localnrow(iproc)*localncol(iproc)
       enddo
    ENDDO

!    call sleep(infpar%mynum*10)
!    IF(mynum.EQ.0)THEN
!    do iproc = 0,numnodes-1
!       Print*,'BUILDSCALAPACK: I am ',infpar%mynum,' should I send to ',Iproc,'CommNodes',CommunicationNodesGlobalSend(iproc)
!       IF(CommunicationNodesGlobalSend(iproc))THEN
!          Print*,'BUILDSCALAPACK: I am ',infpar%mynum,'and I will send this block to ',Iproc
!          CALL ls_output(ASendbuffer(Iproc)%elms,1,nsizesend(iProc),1,1,nsizesend(iProc),1,1,6)
!       ENDIF
!    enddo
!    ENDIF

    !====================================================================================
    !                                Step 2   
    !====================================================================================
    ! Determine who I (infpar%mynum) should recv a block from. Loop through the full lstensor 
    ! and if the iproc owns a lstensor block that I need to build the scalapack matrix. 
    ! I call mem_alloc space for it and recieves
    !====================================================================================
    nsizeRecv = 0
    bufferOffset = 0
    bufferproc = 0 
    CommunicationNodesGlobalRecv=.FALSE.
    !THIS COULD BE DONE USING LOWER DIAGONAL MATRIX -> DMAT SPECIAL CASE TO REDUCE COMMUNICATION
    DO I = 1,TENSOR%G_nLSAO
       G_lsao => TENSOR%G_LSAO(I)
       elms => dummy
       ndim1 = G_lsao%G_nLocal(1)
       ndim2 = G_lsao%G_nLocal(2)
       sA = G_lsao%G_startGlobalOrb(1) 
       sB = G_lsao%G_startGlobalOrb(2) 
       !determine CommunicationNodes,localnrow,localncol
       call scalapack_transform_aux(MAT2,elms,ndim1,ndim2,ndim1,ndim2,sA,sB,&
            & CommunicationNodes,localnrow,localncol,ARecvbuffer,nbuffer,bufferOffset,debug,bufferproc,nodes,&
            & STA_null)
       !g_lsao%mynum owns a lstensor block I need, so I need to recv from   
       IF(CommunicationNodes(mynum))THEN
          CommunicationNodesGlobalRecv(g_lsao%mynum) = .TRUE.
          nsizeRecv(g_lsao%mynum) = nsizeRecv(g_lsao%mynum) + localnrow(mynum)*localncol(mynum)
       ENDIF
    ENDDO
    
    !call mem_alloc buffers for each processor I need to recv a block from, including myself
    do iproc = 0,numnodes-1
       IF(CommunicationNodesGlobalRecv(iproc))THEN
          CALL mem_alloc(Arecvbuffer(Iproc)%elms,nsizeRecv(iproc))
       ENDIF
    enddo

    !====================================================================================
    !                                Step 3   
    !====================================================================================
    !   Send and recieve 
    !====================================================================================
    IF(CommunicationNodesGlobalRecv(mynum))THEN
       !need to send myself a packages so
!       ARecvbuffer(mynum)%elms = ASendbuffer(mynum)%elms
       CALL DCOPY(size(ASendbuffer(mynum)%elms),ASendbuffer(mynum)%elms,1,ARecvbuffer(mynum)%elms,1)
       CALL mem_dealloc(ASendbuffer(mynum)%elms)
!       CommunicationNodesGlobalRecv(mynum) = .FALSE.
       CommunicationNodesGlobalSend(mynum) = .FALSE. 
    ENDIF
!    print*,'START SEND RECIEVE',infpar%mynum
    DO Iproc=0,numnodes-1
       !send to or recieve from Iproc
       if(mynum.EQ.Iproc)THEN
          !recieve from others Iproc2
          DO Iproc2=0,numnodes-1
             IF(Iproc2.EQ.mynum)CYCLE
             IF(CommunicationNodesGlobalRecv(Iproc2))THEN
!                print*,'Iproc',Iproc,'recieves from',Iproc2
!                call ls_mpisendrecv(ARecvbuffer(Iproc2)%elms,nsizeRecv(Iproc2),comm,Iproc2,mynum) !sender=Iproc2 reciever=mynum
                Sender = Iproc2
                call ls_mpisendrecv(ARecvbuffer(Iproc2)%elms,nsizeRecv(Iproc2),comm,sender,mynum) !sender=Iproc2 reciever=mynum

!                print*,'ARecvbuffer(Iproc2)%elms',ARecvbuffer(Iproc2)%elms
             ENDIF
          ENDDO
       else
          !send to Iproc
          IF(CommunicationNodesGlobalSend(Iproc))THEN
!             print*,'Iproc',mynum,'sends to',Iproc
!             call ls_mpisendrecv(ASendbuffer(Iproc)%elms,nsizeSend(Iproc),comm,mynum,Iproc) !sender=mynum reciever=Iproc
             reciever = Iproc
             call ls_mpisendrecv(ASendbuffer(Iproc)%elms,nsizeSend(Iproc),comm,mynum,Reciever) !sender=mynum reciever=Iproc
!             print*,'ASendbuffer(Iproc)%elms',ASendbuffer(Iproc)%elms
          ENDIF
       endif
    ENDDO
!    print*,'DONE SEND RECIEVE',infpar%mynum

    !call mem_dealloc buffer used for sending
    do iproc = 0,numnodes-1
       IF(CommunicationNodesGlobalSend(iproc))THEN
          CALL mem_dealloc(ASendbuffer(Iproc)%elms)
       ENDIF
    enddo

!    call sleep(infpar%mynum*10)
!    IF(mynum.EQ.0)THEN
!       do iproc = 0,numnodes-1
!          print*,'RECV:',CommunicationNodesGlobalRecv(iproc),'mynum',infpar%mynum
!          IF(CommunicationNodesGlobalRecv(iproc))THEN
!             Print*,'BUILDSCALAPACK: I am ',infpar%mynum,' and I recieved this block from ',Iproc
!             CALL ls_output(ARecvbuffer(Iproc)%elms,1,nsizeRecv(iProc),1,1,nsizeRecv(iProc),1,1,6)
!          ENDIF
!       enddo
!    ENDIF

    !====================================================================================
    !                                Step 4   
    !====================================================================================
    !   Fill in Matrix
    !====================================================================================

    !build Matrix from buffers
    bufferOffset = 0
    DO I = 1,TENSOR%G_nLSAO
       G_lsao => TENSOR%G_LSAO(I)
!       print*,'FILL IN MATRIX(',infpar%mynum,'): GLOBAL I=',I
!       print*,'FILL IN MATRIX(',infpar%mynum,')=',CommunicationNodesGlobalRecv(g_lsao%mynum)
       IF(CommunicationNodesGlobalRecv(g_lsao%mynum))THEN !Did I recieve something from g_lsao%mynum
          elms => dummy
          ndim1 = G_lsao%G_nLocal(1)
          ndim2 = G_lsao%G_nLocal(2)
          sA = G_lsao%G_startGlobalOrb(1) 
          sB = G_lsao%G_startGlobalOrb(2) 
          bufferproc = G_lsao%mynum  !build up Matrix from buffer Arecvbuffer(bufferproc)%elms
          !determine localnrow,localncol
          call scalapack_transform_aux(MAT2,elms,ndim1,ndim2,ndim1,ndim2,sA,sB,&
               & CommunicationNodes,localnrow,localncol,ARecvbuffer,&
               & nbuffer,bufferOffset,debug,bufferproc,nodes,STA_null)
          !localnrow is the first dimension of the block residing on this processor
          !and therefor the rowsize of the buffer that should be sent to lsao%mynum 
          localnrow(bufferproc) = localnrow(mynum)
          !Build my local part of matrix MAT from Buffer(lsao%mynum), that I recieved from G_lsao%mynum
          call scalapack_transform_aux(MAT2,elms,ndim1,ndim2,ndim1,ndim2,sA,sB,&
               & CommunicationNodes,localnrow,localncol,ARecvbuffer,&
               & nbuffer,bufferOffset,debug,bufferproc,nodes,STA_BuildAlocalFromBuffer)
          !up date offset 
          IF(CommunicationNodes(mynum))THEN
             bufferOffset(bufferproc) = bufferOffset(bufferproc) + localnrow(mynum)*localncol(mynum)
          ENDIF
       ENDIF
    ENDDO

    !call mem_dealloc the buffers
    do iproc = 0,numnodes-1
       IF(CommunicationNodesGlobalRecv(iproc))THEN
          CALL mem_dealloc(ARecvbuffer(Iproc)%elms)
       ENDIF
    enddo    
    call mem_dealloc(dummy)

!    print*,'BUILDSCA the mat_scalapack_print_local(MAT2)',infpar%mynum
!    call mat_scalapack_print_local(MAT2)
    !most likely not needed. 
    call lsmpi_barrier(comm)
#else
      call lsquit('memdist_lstensor_BuildToScalapack requires VAR_SCALAPACK',-1)
#endif

  end SUBROUTINE memdist_lstensor_BuildToScalapack

!> \brief set the lstensor to zero
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!FIXME WHAT ABOUT ANTISYMMETRY
SUBROUTINE memdist_lstensor_zero_lowerTriangularMat(TENSOR)
implicit none
TYPE(LSTENSOR)     :: TENSOR
!
#ifdef VAR_MPI
INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,Iangmom,Jangmom,Kangmom,Langmom,nOrbA
INTEGER    :: sA,nOrbB,sB,nOrbC,sC,dimA,dimB,jatom,iatom,I2
INTEGER    :: nOrbD,sD,IORB,JORB,KORB,LORB,IELM,IMAT,i1,dim,ielms,ndim5
TYPE(LSAOTENSOR),pointer    :: lsao,lsao2
INTEGER    :: G_Jatom,G_Iatom
!only for 2 dimensional matrices
IF(TENSOR%nbast(3).NE.1.OR.TENSOR%nbast(4).NE.1)THEN
   CALL LSQUIT('Error in lstensor_force_symMat_to_triangularMat1',-1)
ENDIF
!only for square matrices
IF(TENSOR%G_nbast(1).NE.TENSOR%G_nbast(2))THEN
   CALL LSQUIT('Error in lstensor_force_symMat_to_triangularMat2',-1)
ENDIF
IF(TENSOR%G_natom(1).NE.TENSOR%G_natom(2))THEN
   CALL LSQUIT('Error in lstensor_force_symMat_to_triangularMat2',-1)
ENDIF
ndim5 = TENSOR%ndim5
!set lower triangular to zero
DO Jatom = 1,TENSOR%natom(2)
   G_Jatom = TENSOR%G_fullatoms2(Jatom)
   DO Iatom = 1,TENSOR%natom(1)
      G_Iatom = TENSOR%G_fullatoms1(Iatom)
      IF(G_Iatom.GT.G_Jatom)THEN
         !set to zero
         I = TENSOR%INDEX(IATOM,JATOM,1,1)
         IF(I.NE.0)THEN
            lsao => TENSOR%LSAO(I)
            dim = ndim5*lsao%nelms
            call ls_dzero(lsao%elms,dim)
         ENDIF
      ELSEIF(G_Iatom.EQ.G_Jatom)THEN
         I = TENSOR%INDEX(IATOM,JATOM,1,1)
         IF(I.NE.0)THEN
            lsao => TENSOR%LSAO(I)
            dimA = lsao%nLocal(1)
            dimB = lsao%nLocal(2)
            IF(dimA.NE.dimB)call lsquit('errornjbhjbhjbkjbhj',-1)
            call set_lowerelms_triangular_zero(lsao%elms,dimA,dimB,ndim5)
         ENDIF
      ENDIF
   ENDDO
ENDDO
#else
    call lsquit('called memdist_lstensor_SetupFullinfo but no MPI',-1)
#endif
END SUBROUTINE memdist_lstensor_zero_lowerTriangularMat

!!$!> \brief set the lstensor to zero
!!$!> \author T. Kjaergaard
!!$!> \date 2010
!!$!> \param TENSOR the lstensor
!!$!FIXME WHAT ABOUT ANTISYMMETRY
!!$SUBROUTINE memdist_lstensor_dobbel_upperTriangularMat(TENSOR)
!!$implicit none
!!$TYPE(LSTENSOR)     :: TENSOR
!!$!
!!$INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,Iangmom,Jangmom,Kangmom,Langmom,nOrbA
!!$INTEGER    :: sA,nOrbB,sB,nOrbC,sC,dimA,dimB,jatom,iatom,I2
!!$INTEGER    :: nOrbD,sD,IORB,JORB,KORB,LORB,IELM,IMAT,i1,dim,ielms,ndim5
!!$TYPE(LSAOTENSOR),pointer    :: lsao,lsao2
!!$INTEGER    :: G_Jatom,G_Iatom
!!$!only for 2 dimensional matrices
!!$IF(TENSOR%nbast(3).NE.1.OR.TENSOR%nbast(4).NE.1)THEN
!!$   CALL LSQUIT('Error in lstensor_force_symMat_to_triangularMat1',-1)
!!$ENDIF
!!$!only for square matrices
!!$IF(TENSOR%G_nbast(1).NE.TENSOR%G_nbast(2))THEN
!!$   CALL LSQUIT('Error in lstensor_force_symMat_to_triangularMat2',-1)
!!$ENDIF
!!$IF(TENSOR%G_natom(1).NE.TENSOR%G_natom(2))THEN
!!$   CALL LSQUIT('Error in lstensor_force_symMat_to_triangularMat2',-1)
!!$ENDIF
!!$ndim5 = TENSOR%ndim5
!!$!set lower triangular to zero
!!$DO Jatom = 1,TENSOR%natom(2)
!!$   G_Jatom = TENSOR%G_fullatoms2(Jatom)
!!$   DO Iatom = 1,TENSOR%natom(1)
!!$      G_Iatom = TENSOR%G_fullatoms1(Iatom)
!!$      IF(G_Iatom.LT.G_Jatom)THEN
!!$         !set to zero
!!$         I = TENSOR%INDEX(IATOM,JATOM,1,1)
!!$         IF(I.NE.0)THEN
!!$            lsao => TENSOR%LSAO(I)
!!$            dim = ndim5*lsao%nelms
!!$            call dscal(dim,2.0E0_realk,lsao%elms,1)
!!$         ENDIF
!!$      ELSEIF(G_Iatom.EQ.G_Jatom)THEN
!!$         I = TENSOR%INDEX(IATOM,JATOM,1,1)
!!$         IF(I.NE.0)THEN
!!$            lsao => TENSOR%LSAO(I)
!!$            dimA = lsao%nLocal(1)
!!$            dimB = lsao%nLocal(2)
!!$            IF(dimA.NE.dimB)call lsquit('errornjbhjbhjbkjbhj',-1)
!!$            call set_upperelms_triangular_doobel(lsao%elms,dimA,dimB,ndim5)
!!$         ENDIF
!!$      ENDIF
!!$   ENDDO
!!$ENDDO
!!$END SUBROUTINE memdist_lstensor_dobbel_upperTriangularMat
!!$
!!$subroutine set_upperelms_triangular_dobbel(elms,dimenA,dimenB,nmat)
!!$  implicit none
!!$  integer,intent(in) :: dimenA,dimenB,nmat
!!$  real(realk) :: elms(dimenA,dimenB,nmat)
!!$  !
!!$  integer :: A,B,IMAT
!!$  Do IMAT = 1,NMAT
!!$     DO B=1,dimenB
!!$        DO A=1,B-1
!!$           elms(A,B,IMAT) = 2.0E0_realk*elms(A,B,IMAT) 
!!$        ENDDO
!!$     ENDDO
!!$  ENDDO
!!$end subroutine set_upperelms_triangular_dobbel

subroutine memdist_lstensor_SetupFullinfo(TENSOR,AO1,AO2,nbast1,nbast2,&
       & ndim5,useAO1,useAO2,ODscreen,lupri)
    implicit none
    INTEGER            :: nbast1,nbast2,lupri,ndim5
    TYPE(LSTENSOR)     :: TENSOR
    TYPE(AOITEM),target  :: AO1,AO2
    logical   :: useAO1,useAO2,ODscreen
    !
#ifdef VAR_MPI
    TYPE(AOITEMPOINTER) :: AOT(2)
    TYPE(AOITEM),target :: AOTE
    ! 
    INTEGER :: natom(2),nbat(2),AOTbatch(2),batAOT(2),nAng(2)
    !
    INTEGER :: nbastI,Ibat,Iangmom,Iatom,IORB,nOrbA,sA
    INTEGER :: nbastJ,Jbat,Jangmom,Jatom,JORB,nOrbB,sB
    ! COMMON
    INTEGER :: DIM,IELM,IMAT2,I,J,iAO,jAO,maxBat,maxAng,nLocal
    INTEGER :: nElms,startLocalOrb,nAngG,batAOTG,dim2(2),maxnatom
    INTEGER :: maxnlocal
    LOGICAL :: ERROR,Screen,ScreenAtom
    LOGICAL,pointer :: nonscreenAB(:,:)
    integer(kind=long) :: nmemsize
    TYPE(GLOBALLSAOTENSOR),pointer    :: g_lsao
    integer,pointer :: inttmp(:,:)
    integer,pointer :: inttmp2(:)
    !FIXME THIS SUBROUTINE SHOULD BE OMP parallel!
    call LSTENSOR_nullify(TENSOR)
    call SET_EMPTY_AO(AOTE)
    AOT(1)%p => AO1
    AOT(2)%p => AO2
    IF(.NOT.useAO1) AOT(1)%p => AOTE
    IF(.NOT.useAO2) AOT(2)%p => AOTE
    DO IAO = 1,2
       natom(iAO) = AOT(iAO)%p%natoms
    ENDDO
    maxnatom = MAX(natom(1),natom(2))
    nullify(TENSOR%G_nAOBATCH)
    call mem_alloc(TENSOR%G_nAOBATCH,maxnatom,2)
    inttmp => TENSOR%G_nAOBATCH
    DO IAO = 1,2
       inttmp2 => AOT(IAO)%p%ATOMICnBatch
       DO IATOM = 1,natom(iAO)
          inttmp(IATOM,IAO) = inttmp2(IATOM)
       ENDDO
    ENDDO
    DO IAO = 1,2
       TENSOR%G_natom(iAO) = natom(iAO)
    ENDDO
    TENSOR%G_nbast(1) = nbast1
    TENSOR%G_nbast(2) = nbast2
    DO IAO = 1,2
       TENSOR%nbatches(iAO) = AOT(iAO)%p%nbatches
    ENDDO
    TENSOR%ndim5=ndim5

    call mem_alloc(nonscreenAB,AOT(1)%p%natoms,AOT(2)%p%natoms)
    IF(ODscreen)THEN
       call build_atomicODScreen(AOT(1)%p,AOT(2)%p,nonscreenAB,AOT(1)%p%natoms,AOT(2)%p%natoms)
    else
       nonscreenAB = .true.
    ENDIF

    I = 0
    DO Jatom = 1,natom(2)
       DO Iatom = 1,natom(1)
          IF(nonscreenAB(Iatom,Jatom))THEN                          
             I = I +1
          ENDIF
       ENDDO
    ENDDO
    TENSOR%G_nLSAO = I
    NULLIFY(TENSOR%G_LSAO)
    NULLIFY(TENSOR%G_INDEX)
    CALL MEM_ALLOC(TENSOR%G_LSAO,I)
    CALL MEM_ALLOC(TENSOR%G_INDEX,natom(1),natom(2))

    DO Jatom = 1,natom(2)
       DO Iatom = 1,natom(1)
          TENSOR%G_INDEX(Iatom,Jatom) = 0 !if 0 lsaotensor not call mem_allocd 
       ENDDO
    ENDDO

    I = 0
    maxnlocal=0
    AOTbatch(2)=0
    DO Jatom = 1,natom(2)
       nbat(2) = AOT(2)%p%ATOMICnBatch(JATOM)
       AOTbatch(1)=0
       DO Iatom = 1,natom(1)
          nbat(1) = AOT(1)%p%ATOMICnBatch(IATOM)
          IF(nonscreenAB(Iatom,Jatom))THEN                          
             I=I+1
             TENSOR%G_INDEX(IATOM,JATOM) = I
             g_lsao => TENSOR%G_LSAO(I)
!                      g_lsao%ATOM(1) = Iatom
!                      g_lsao%ATOM(2) = Jatom
!                      g_lsao%AOBATCH(1) = AOTbatch(1)
!                      g_lsao%AOBATCH(2) = AOTbatch(2)
             do IAO=1,2
                g_lsao%G_startGlobalOrb(iAO)=0
             ENDDO
             do IAO=1,2
                nLocal = 0
                batAOTG = AOTbatch(IAO)
                g_lsao%G_startGlobalOrb(iAO)=AOT(iAO)%p%BATCH(batAOTG+1)%startOrbital(1)
                DO Ibat = 1,nbat(iAO)
                   batAOTG = batAOTG+1
                   nAngG = AOT(iAO)%p%BATCH(batAOTG)%nAngmom 
                   DO Iangmom = 1,nAngG
                      nLocal = nLocal + AOT(iAO)%p%BATCH(batAOTG)%nContracted(Iangmom)*&
                           & AOT(iAO)%p%BATCH(batAOTG)%nOrbComp(Iangmom)
                   ENDDO
                ENDDO
                g_lsao%G_nLocal(iAO) = nLocal
                maxnlocal = MAX(maxnlocal,nLocal)
             ENDDO
          ENDIF !AB screen
          AOTbatch(1) = AOTbatch(1) + nbat(1)
       ENDDO
       AOTbatch(2) = AOTbatch(2) + nbat(2)
    ENDDO
    TENSOR%G_maxnlocal = maxnlocal
    TENSOR%nLSAO = I
    call mem_dealloc(nonscreenAB)
    call FREE_EMPTY_AO(AOTE)

#else
    call lsquit('called memdist_lstensor_SetupFullinfo but no MPI',-1)
#endif
end subroutine memdist_lstensor_SetupFullinfo

  !> \brief 
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param TENSOR1 the original lstensor
  !> \param TENSOR2 the new sublstensor
  SUBROUTINE memdist_lstensor_SetupLocalinfo(TENSOR,AO1,AO2,nbast1,nbast2,&
       & useAO1,useAO2,natoms1,natoms2,atoms1,atoms2,lupri)
    implicit none
    TYPE(LSTENSOR)     :: TENSOR
    integer,intent(in) :: lupri,natoms1,natoms2,nbast1,nbast2
    integer,intent(in) :: atoms1(natoms1),atoms2(natoms2)
    TYPE(AOITEM),target  :: AO1,AO2
    logical   :: useAO1,useAO2
#ifdef VAR_MPI
    !
    TYPE(AOITEMPOINTER) :: AOT(2)
    TYPE(AOITEM),target :: AOTE
    INTEGER   :: I1,Ibat,Jbat,ndim5,dim,ielms,I2,iangmom,nLocal(2)
    INTEGER   :: IATOM1,IATOM2,JATOM1,JATOM2,natomG,iAO,nAngG,offset
    INTEGER   :: maxnAtom,nAngmomA,nAngmomB,batAOTG,startLocalOrb
    Integer   :: ibatch2,ibatch1,jbatch1,jbatch2,nbatB,nbatA,IATOM,lstmem_index
    Integer   :: nbat1(natoms1),nbat2(natoms2),nstart1(natoms1),nstart2(natoms2)
    integer  :: nbatches1,nbatches2,nbatJ,nbatI,Jbatch,Ibatch,Lbatch,Kbatch
    INTEGER :: natom(2),nbat(2),AOTbatch(2),batAOT(2),nAng(2),maxBat,maxAng
    integer(kind=long) :: nmemsize,AllocInt,AllocIntS,AllocRealk,nsize
    TYPE(LSAOTENSOR),pointer    :: lsao
!    REAL(REALK),pointer     :: elms(:)
!    TYPE(LSAOTENSOR),pointer    :: lsao2
!    REAL(REALK),pointer     :: elms2(:)
!    print*,'TENSOR in memdist_lstensor_SetupLocalinfo'
!    call lstensor_print(TENSOR1,6)

    call SET_EMPTY_AO(AOTE)
    AOT(1)%p => AO1
    AOT(2)%p => AO2
    IF(.NOT.useAO1) AOT(1)%p => AOTE
    IF(.NOT.useAO2) AOT(2)%p => AOTE

    natom(1) = natoms1
    natom(2) = natoms2

    call mem_alloc(TENSOR%G_fullatoms1,natoms1)
    TENSOR%G_fullatoms1 = atoms1
    call mem_alloc(TENSOR%G_fullatoms2,natoms2)
    TENSOR%G_fullatoms2 = atoms2

    nullify(TENSOR%nAOBATCH)
    IF(ASSOCIATED(TENSOR%G_nAOBATCH))THEN
       maxnatom = MAX(natoms1,natoms2)
       call mem_alloc(TENSOR%nAOBATCH,maxnatom,2)
       nsize = size(TENSOR%nAOBATCH,KIND=long)*mem_intsize
       call mem_allocated_mem_lstensor(nsize)
       DO IATOM = 1,natoms1
          IATOM2 = atoms1(IATOM)
          TENSOR%nAOBATCH(IATOM,1) = TENSOR%G_nAOBATCH(IATOM2,1)
       ENDDO
       DO IATOM = 1,natoms2
          IATOM2 = atoms2(IATOM)
          TENSOR%nAOBATCH(IATOM,2) = TENSOR%G_nAOBATCH(IATOM2,2)
       ENDDO
    ENDIF
    IF(TENSOR%maggradienttensor)THEN
       CALL LSQUIT(' TENSOR%maggradienttensor not implemtne memdist',-1)
    ELSEIF(TENSOR%gradienttensor)THEN
       CALL LSQUIT(' TENSOR%gradienttensor not implemtne memdist',-1)
    ELSEIF(TENSOR%pChargetensor)THEN
       CALL LSQUIT(' TENSOR%pChargetensor not implemtne memdist',-1)
    ELSEIF(TENSOR%Econtrib)THEN
       CALL LSQUIT(' TENSOR%Econtrib not implemtne memdist',-1)
    ELSEIF(TENSOR%Screentensor)THEN
       CALL LSQUIT(' TENSOR%Screentensor not implemtne memdist',-1)
    ELSE
       TENSOR%natom(1) = natoms1
       TENSOR%natom(2) = natoms2
       TENSOR%natom(3) = 1
       TENSOR%natom(4) = 1
       TENSOR%nbast(1) = nbast1
       TENSOR%nbast(2) = nbast2
       TENSOR%nbast(3) = 1
       TENSOR%nbast(4) = 1
       ndim5 = TENSOR%ndim5
       TENSOR%nLSAO  = natoms1*natoms2
       NULLIFY(TENSOR%INDEX)
       CALL MEM_ALLOC(TENSOR%INDEX,natoms1,natoms2,1,1)
       nsize = size(TENSOR%INDEX,KIND=long)*mem_intsize
       call mem_allocated_mem_lstensor(nsize)
       !memory estimation
       AllocInt = 0
       AllocIntS = 0
       AllocRealk = 0
       I2=0
       AOTbatch(2)=0
       DO Jatom2 = 1,natoms2
          Jatom1 = atoms2(Jatom2)
          Ibatch = 0
          nbat(2) = AOT(2)%p%ATOMICnBatch(JATOM2)
          AOTbatch(1)=0
          DO Iatom2 = 1,natoms1
             Iatom1 = atoms1(Iatom2)
             I1 = TENSOR%G_INDEX(IATOM1,JATOM1)             
             IF(I1.NE.0)THEN
                nbat(1) = AOT(1)%p%ATOMICnBatch(IATOM2)
                I2=I2+1
                dim = TENSOR%G_LSAO(I1)%G_nLocal(1)*TENSOR%G_LSAO(I1)%G_nLocal(2)*ndim5
                AllocRealk = AllocRealk + DIM
                maxBat = MAX(nBat(1),nBat(2))
                AllocInt = AllocInt + 2*maxBat
                maxAng = 0
                do IAO=1,2
                   batAOTG = AOTbatch(IAO)
                   DO Ibat = 1,nbat(iAO)
                      batAOTG = batAOTG+1
                      maxAng = MAX(maxAng,AOT(iAO)%p%BATCH(batAOTG)%nAngmom)
                   ENDDO
                enddo
                AllocInt = AllocInt + 6*maxBat*maxAng
             ENDIF
             AOTbatch(1) = AOTbatch(1) + nbat(1)
          enddo
          AOTbatch(2) = AOTbatch(2) + nbat(2)
       enddo
       TENSOR%nLSAO  = I2
       NULLIFY(TENSOR%LSAO)
       CALL MEM_ALLOC(TENSOR%LSAO,TENSOR%nLSAO)
       call init_lstensorMem(AllocInt,AllocRealk,AllocInts,lstmem_index)
       call zero_lstensorMem
       TENSOR%lstmem_index = lstmem_index

       I2=0
       AOTbatch(2)=0!Jbatch = 0
       DO Jatom2 = 1,natoms2
          Jatom1 = atoms2(Jatom2)
          Ibatch = 0
          nbat(2) = AOT(2)%p%ATOMICnBatch(JATOM2)
          AOTbatch(1)=0
          DO Iatom2 = 1,natoms1
             Iatom1 = atoms1(Iatom2)
             I1 = TENSOR%G_INDEX(IATOM1,JATOM1)             
             IF(I1.NE.0)THEN
                nbat(1) = AOT(1)%p%ATOMICnBatch(IATOM2)
                I2=I2+1
                TENSOR%INDEX(IATOM2,JATOM2,1,1) = I2
                lsao => TENSOR%LSAO(I2)
                !nLocal
                lsao%nLocal(1) = TENSOR%G_LSAO(I1)%G_nLocal(1)
                lsao%nLocal(2) = TENSOR%G_LSAO(I1)%G_nLocal(2)
                !nelms + elms
                dim = lsao%nLocal(1)*lsao%nLocal(2)
                TENSOR%LSAO(I2)%nelms = dim
                dim = dim*ndim5
                call mem_lstpointer_alloc(TENSOR%LSAO(I2)%elms,dim)
                !nAngmom
                maxBat = MAX(nBat(1),nBat(2))
                lsao%maxBat = maxBat
                call mem_lstpointer_alloc(lsao%nAngmom,maxBat*2)
                maxAng = 0
                do IAO=1,2
                   batAOTG = AOTbatch(IAO)
                   DO Ibat = 1,nbat(iAO)
                      batAOTG = batAOTG+1
                      lsao%nAngmom(Ibat+(IAO-1)*maxBat) = AOT(iAO)%p%BATCH(batAOTG)%nAngmom
                      maxAng = MAX(maxAng,AOT(iAO)%p%BATCH(batAOTG)%nAngmom)
                   ENDDO
                enddo
                !nOrb , startLocalOrb, startGlobalOrb
                NULLIFY(lsao%nOrb)
                NULLIFY(lsao%startLocalOrb)
                NULLIFY(lsao%startGlobalOrb)
                CALL MEM_LSTPOINTER_ALLOC(lsao%nOrb,maxAng*maxBat*2)
                CALL MEM_LSTPOINTER_ALLOC(lsao%startLocalOrb,maxAng*maxBat*2)
                CALL MEM_LSTPOINTER_ALLOC(lsao%startGlobalOrb,maxAng*maxBat*2)
                lsao%maxAng = maxAng
                do IAO=1,2
                   nLocal(iAO) = 0
                   batAOTG = AOTbatch(IAO)
                   startLocalOrb = 1
                   DO Ibat = 1,nbat(iAO)
                      batAOTG = batAOTG+1
                      nAngG = AOT(iAO)%p%BATCH(batAOTG)%nAngmom 
                      offset = (ibat-1)*maxAng+(iAO-1)*maxAng*maxBat
                      DO Iangmom = 1,nAngG
                         lsao%nOrb(Iangmom+offset)=AOT(iAO)%p%BATCH(batAOTG)%nContracted(Iangmom)*&
                              & AOT(iAO)%p%BATCH(batAOTG)%nOrbComp(Iangmom)
                         nLocal(iAO) = nLocal(iAO) + lsao%nOrb(Iangmom+offset)
                         lsao%startGlobalOrb(Iangmom+offset)=AOT(iAO)%p%BATCH(batAOTG)%startOrbital(Iangmom)
                         lsao%startLocalOrb(Iangmom+offset)=startLocalOrb
                         startLocalOrb = startLocalOrb + lsao%nOrb(Iangmom+offset)
                      ENDDO
                   ENDDO
                   lsao%nLocal(iAO) = nLocal(iAO)
                   IF(nLocal(iAO).NE.TENSOR%G_LSAO(I1)%G_nLocal(iAO))THEN
                      print*,'iAO',iAO,'infpar%mynum',infpar%mynum
                      print*,'nLocal(iAO)',nLocal(iAO),'I2',I2,'IATOM',IATOM2,JATOM2,&
                           & 'infpar%mynum',infpar%mynum
                      print*,'TENSOR%G_LSAO(I2)%G_nLocal(iAO)',TENSOR%G_LSAO(I2)%G_nLocal(iAO),'I1',&
     &                        I1,'IATOM',IATOM1,JATOM1,'infpar%mynum',infpar%mynum
                      CALL LSQUIT('nLocal mismatch localinfo',-1)
                   ENDIF
                ENDDO
                TENSOR%LSAO(I2)%ATOM(1) = IATOM2
                TENSOR%LSAO(I2)%ATOM(2) = JATOM2
                TENSOR%LSAO(I2)%ATOM(3) = 1
                TENSOR%LSAO(I2)%ATOM(4) = 1
                TENSOR%LSAO(I2)%AOBATCH(1) = AOTbatch(1)
                TENSOR%LSAO(I2)%AOBATCH(2) = AOTbatch(2)
                TENSOR%LSAO(I2)%AOBATCH(3) = 1
                TENSOR%LSAO(I2)%AOBATCH(4) = 1
             ELSE
                TENSOR%INDEX(IATOM2,JATOM2,1,1) = 0
             ENDIF
             AOTbatch(1) = AOTbatch(1) + nbat(1)
          enddo
          AOTbatch(2) = AOTbatch(2) + nbat(2)
       enddo
       IF(I2.NE.TENSOR%nLSAO)call lsquit('lstensor dim mismatch I2',-1)
    endif
    call FREE_EMPTY_AO(AOTE)
#else
    call lsquit('called memdist_lstensor_SetupFullinfo but no MPI',-1)
#endif
  end SUBROUTINE memdist_lstensor_SetupLocalinfo

  !> \brief 
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param TENSOR1 the original lstensor
  !> \param TENSOR2 the new sublstensor
  SUBROUTINE build_empty_sublstensor(TENSOR)
    implicit none
    TYPE(LSTENSOR)     :: TENSOR
#ifdef VAR_MPI
    !
    TYPE(AOITEMPOINTER) :: AOT(2)
    TYPE(AOITEM),target :: AOTE
    integer :: lstmem_index
    integer(kind=long) :: nmemsize,AllocInt,AllocIntS,AllocRealk,nsize

    call SET_EMPTY_AO(AOTE)
    AOT(1)%p => AOTE
    AOT(2)%p => AOTE

    call mem_alloc(TENSOR%G_fullatoms1,1)
    TENSOR%G_fullatoms1(1) = 1
    call mem_alloc(TENSOR%G_fullatoms2,1)
    TENSOR%G_fullatoms2(1) = 1

    nullify(TENSOR%nAOBATCH)
    IF(ASSOCIATED(TENSOR%G_nAOBATCH))THEN
       call mem_alloc(TENSOR%nAOBATCH,1,2)
       nsize = size(TENSOR%nAOBATCH,KIND=long)*mem_intsize
       call mem_allocated_mem_lstensor(nsize)
       TENSOR%nAOBATCH(1,1) = 1
       TENSOR%nAOBATCH(1,2) = 1
    ENDIF

    TENSOR%natom(1) = 0
    TENSOR%natom(2) = 0
    TENSOR%natom(3) = 0
    TENSOR%natom(4) = 0
    TENSOR%nbast(1) = 0
    TENSOR%nbast(2) = 0
    TENSOR%nbast(3) = 0
    TENSOR%nbast(4) = 0
    TENSOR%nLSAO  = 0
    NULLIFY(TENSOR%INDEX)
    CALL MEM_ALLOC(TENSOR%INDEX,1,1,1,1)
    nsize = size(TENSOR%INDEX,KIND=long)*mem_intsize
    call mem_allocated_mem_lstensor(nsize)
    !memory estimation
    CALL MEM_ALLOC(TENSOR%LSAO,1)
    AllocInt = 1
    AllocRealk = 1
    AllocInts = 1
    call init_lstensorMem(AllocInt,AllocRealk,AllocInts,lstmem_index)
    call zero_lstensorMem
    TENSOR%lstmem_index = lstmem_index
    call FREE_EMPTY_AO(AOTE)
#else
    call lsquit('called build_empty_sublstensor but no MPI',-1)
#endif
  end SUBROUTINE build_empty_sublstensor

  !> \brief copy an lstensor to a new lstensor
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param TENSOR1 the original lstensor
  !> \param TENSOR2 the new sublstensor
  SUBROUTINE alloc_build_empty_sublstensor(TENSOR2)
    implicit none
    TYPE(LSTENSOR)     :: TENSOR2
    integer(kind=long) :: nmemsize
    call LSTENSOR_nullify(TENSOR2)
  end SUBROUTINE alloc_build_empty_sublstensor

  !> \brief build full 5 dimensional array from an lstensor 
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param TENSOR the lstensor
  !> \param MAT134 the 5 dimensional array
  !> \param dim1 the size of the 1. dimension
  !> \param dim2 the size of the 2. dimension
  !> \param dim3 the size of the 3. dimension
  !> \param dim4 the size of the 4. dimension
  !> \param dim5 the size of the 5. dimension
  SUBROUTINE Build_full_5dim_from_lstensor(TENSOR,MAT5dim,dim1,dim2,dim3,dim4,dim5)
    implicit none
    integer          :: dim1,dim2,dim3,dim4,dim5
    TYPE(LSTENSOR)   :: TENSOR
    REAL(REALK),target :: Mat5dim(dim1,dim2,dim3,dim4,dim5)
    !
    INTEGER    :: I,IMAT
    LOGICAL    :: ERROR,MATELMZERO
    type(LSTAUXITEM) :: LSTAUX
    IF(TENSOR%pchargetensor)THEN
       IF(dim2+dim3+dim4.NE.3)then
          print*,'dim2,dim3,dim4',dim2,dim3,dim4,'dim1,TENSOR%nmat',dim1,TENSOR%ndim5
          call lsquit('error using pcharge in ',-1)
       endif
       DO I = 1,TENSOR%nAtom(1)
          DO IMAT = 1,dim5
             Mat5dim(I,1,1,1,IMAT)=TENSOR%LSAO(I)%elms(IMAT)
          ENDDO
       ENDDO
    ELSEIF(TENSOR%Econtrib)THEN
       DO IMAT = 1,dim5
          Mat5dim(1,1,1,1,IMAT)=TENSOR%LSAO(1)%elms(IMAT)
       ENDDO
    ELSE
       ERROR=.FALSE.
       IF(TENSOR%nbast(1) .NE. dim1) ERROR = .TRUE.
       IF(TENSOR%nbast(2) .NE. dim2) ERROR = .TRUE.
       IF(TENSOR%nbast(3) .NE. dim3) ERROR = .TRUE.
       IF(TENSOR%nbast(4) .NE. dim4) ERROR = .TRUE.
       IF(ERROR)THEN
          print*,'LSTENSOR%nbast1',TENSOR%nbast(1),'dim1',dim1
          print*,'LSTENSOR%nbast2',TENSOR%nbast(2),'dim2',dim2
          print*,'LSTENSOR%nbast3',TENSOR%nbast(3),'dim3',dim3
          print*,'LSTENSOR%nbast4',TENSOR%nbast(4),'dim4',dim4
          !   print*,'LSTENSOR%nbast5',TENSOR%nbast5,'dim5',dim5
          CALL lsQUIT('ERROR: Build_full_5dim_from_lstensor dim do not match',-1)
       ENDIF
       call ls_dzero(Mat5dim,dim1*dim2*dim3*dim4*dim5)
       call nullify_LSTAUXITEM(LSTAUX)
       LSTAUX%Mat5dim => Mat5dim
       MATELMZERO = .FALSE.
       call BuildFromlstensor_aux(TENSOR,0,0,0,MATELMZERO,LSTAUX,FiveDim_FuncLstSub)
    ENDIF
  END SUBROUTINE BUILD_FULL_5DIM_FROM_LSTENSOR

  !> \brief build full 4 dimensional array from an lstensor 
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param TENSOR the lstensor
  !> \param Mat4dim the 4 dimensional array
  !> \param dim1 the size of the 1. dimension
  !> \param dim2 the size of the 2. dimension
  !> \param dim3 the size of the 3. dimension
  !> \param dim4 the size of the 4. dimension
  SUBROUTINE Build_full_4dim_from_lstensor(TENSOR,Mat4dim,dim1,dim2,dim3,dim4)
    implicit none
    integer          :: dim1,dim2,dim3,dim4
    TYPE(LSTENSOR)   :: TENSOR
    REAL(REALK),target      :: Mat4dim(dim1,dim2,dim3,dim4)
    !
    type(LSTAUXITEM) :: LSTAUX
    INTEGER    :: I,Ibat,Jbat,Kbat,Lbat
    LOGICAL    :: ERROR,MATELMZERO
    IF(TENSOR%pchargetensor)CALL lsQUIT('pcharge Build_full_4dim_from_lstensor error',-1)
    IF(TENSOR%Econtrib)CALL lsQUIT('Econtrib Build_full_4dim_from_lstensor error',-1)
    IF(TENSOR%ndim5.NE.1)call lsquit('dim mismatch in Build_full_4dim_from_lstensor',-1)
    ERROR=.FALSE.
    IF(TENSOR%nbast(1) .NE. dim1) ERROR = .TRUE.
    IF(TENSOR%nbast(2) .NE. dim2) ERROR = .TRUE.
    IF(TENSOR%nbast(3) .NE. dim3) ERROR = .TRUE.
    IF(TENSOR%nbast(4) .NE. dim4) ERROR = .TRUE.
    IF(ERROR)THEN
       print*,'LSTENSOR%nbast1',TENSOR%nbast(1),'dim1',dim1
       print*,'LSTENSOR%nbast2',TENSOR%nbast(2),'dim2',dim2
       print*,'LSTENSOR%nbast3',TENSOR%nbast(3),'dim3',dim3
       print*,'LSTENSOR%nbast4',TENSOR%nbast(4),'dim4',dim4
       !   print*,'LSTENSOR%nbast5',TENSOR%nbast5,'dim5',dim5
       CALL lsQUIT('ERROR: Build_full_5dim_from_lstensor dim do not match',-1)
    ENDIF
    DO Lbat = 1,dim4
       DO Kbat = 1,dim3
          DO Jbat = 1,dim2
             DO Ibat = 1,dim1
                Mat4dim(Ibat,Jbat,Kbat,Lbat) = 0E0_realk
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    call nullify_LSTAUXITEM(LSTAUX)
    LSTAUX%Mat4dim => Mat4dim
    MATELMZERO=.FALSE.
    call BuildFromlstensor_aux(TENSOR,0,0,0,MATELMZERO,LSTAUX,FourDim_FuncLstSub)
  END SUBROUTINE Build_full_4dim_from_lstensor

  !> \brief build full 3 dimensional array from an lstensor 
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param TENSOR the lstensor
  !> \param Mat3dim the 3 dimensional array
  !> \param dim1 the size of the 1. dimension
  !> \param dim2 the size of the 2. dimension
  !> \param dim3 the size of the 3. dimension
  SUBROUTINE Build_full_3dim_from_lstensor(TENSOR,Mat3dim,dim1,dim2,dim3)
    implicit none
    integer          :: dim1,dim2,dim3
    TYPE(LSTENSOR)   :: TENSOR
    REAL(REALK),target      :: Mat3dim(dim1,dim2,dim3)
    !
    type(LSTAUXITEM) :: LSTAUX
    INTEGER    :: I,J,Ibat,Jbat,IMAT,nbast(4),DIMI,DIMJ,SUM1,SUM2
    LOGICAL    :: ERROR,MATELMZERO

    IF(TENSOR%pchargetensor)CALL lsQUIT('pcharge Build_full_3dim_from_lstensor error',-1)
    IF(TENSOR%Econtrib)CALL lsQUIT('Econtrib Build_full_3dim_from_lstensor error',-1)
    IF(TENSOR%ndim5.NE.1)call lsquit('dim mismatch in Build_full_3dim_from_lstensor',-1)
    IF(TENSOR%ndim5.NE.dim3)CALL lsQUIT('called Build_full_3dim_from_lstensor with TENSOR%ndim5.NE.dim3',-1)    
    nbast(1)=TENSOR%nbast(1)
    nbast(2)=TENSOR%nbast(2)
    nbast(3)=TENSOR%nbast(3)
    nbast(4)=TENSOR%nbast(4)
    DO I=1,4
       IF(dim1 .EQ.nbast(I))THEN
          DIMI = I
          EXIT
       ENDIF
    ENDDO
    DO J=1,4
       IF(J.NE.DIMI)THEN
          IF(dim2 .EQ.nbast(J))THEN
             DIMJ = J
             EXIT
          ENDIF
       ENDIF
    ENDDO
    SUM1=nbast(DIMI)+nbast(DIMJ)
    SUM2=dim1+dim2
    IF(SUM1.NE.SUM2)CALL lsQUIT('ERROR: Build_mat_from_lst dim do not match',-1)
    do IMAT=1,dim3
       DO Jbat = 1,dim2
          DO Ibat = 1,dim1
             Mat3dim(Ibat,Jbat,IMAT) = 0E0_realk
          ENDDO
       ENDDO
    ENDDO
    call nullify_LSTAUXITEM(LSTAUX)
    LSTAUX%Mat3dim => Mat3dim
    MATELMZERO = .FALSE.
    IF(DIMI.EQ.1.AND.DIMJ.EQ.2)THEN
       call BuildFromlstensor_aux(TENSOR,DIMI,DIMJ,0,MATELMZERO,LSTAUX,ThreeDim12_FuncLstSub)
    ELSEIF(DIMI.EQ.1.AND.DIMJ.EQ.3)THEN
       call BuildFromlstensor_aux(TENSOR,DIMI,DIMJ,0,MATELMZERO,LSTAUX,ThreeDim13_FuncLstSub)
    ELSE
       call BuildFromlstensor_aux(TENSOR,DIMI,DIMJ,0,MATELMZERO,LSTAUX,ThreeDim_FuncLstSub)
    ENDIF

  END SUBROUTINE Build_full_3dim_from_lstensor

  !> \brief build full 2 dimensional array from an lstensor 
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param TENSOR the lstensor
  !> \param MAT134 the 2 dimensional array
  !> \param dim1 the size of the 1. dimension
  !> \param dim2 the size of the 2. dimension
  SUBROUTINE Build_full_2dim_from_lstensor(TENSOR,MAT2dim,dim1,dim2)
    implicit none
    integer          :: dim1,dim2
    TYPE(LSTENSOR)   :: TENSOR
    REAL(REALK),target      :: MAT2dim(dim1,dim2)
    !
    type(LSTAUXITEM) :: LSTAUX
    INTEGER    :: I,J,Ibat,Jbat,DIMI,DIMJ,SUM1,SUM2
    LOGICAL    :: ERROR,MATELMZERO
    INTEGER    :: nbast(4)
    IF(TENSOR%pchargetensor)CALL lsQUIT('pcharge Build_full_2dim_from_lstensor error',-1)
    IF(TENSOR%Econtrib)CALL lsQUIT('Econtrib Build_full_2dim_from_lstensor error',-1)
    IF(Tensor%screenTensor)call lsquit('screen in Build_full_2dim_from_lstensor error',-1)
    ERROR=.FALSE.
    nbast(1)=TENSOR%nbast(1)
    nbast(2)=TENSOR%nbast(2)
    nbast(3)=TENSOR%nbast(3)
    nbast(4)=TENSOR%nbast(4)
    DO I=1,4
       IF(dim1 .EQ.nbast(I))THEN
          DIMI = I
          EXIT
       ENDIF
    ENDDO
    DO J=1,4
       IF(J.NE.DIMI)THEN
          IF(dim2 .EQ.nbast(J))THEN
             DIMJ = J
             EXIT
          ENDIF
       ENDIF
    ENDDO
    SUM1=nbast(DIMI)+nbast(DIMJ)
    SUM2=dim1+dim2
    IF(SUM1.NE.SUM2)CALL lsQUIT('ERROR: Build_mat_from_lst dim dont match',-1)
    DO Jbat = 1,dim2
       DO Ibat = 1,dim1
          MAT2dim(Ibat,Jbat) = 0E0_realk
       ENDDO
    ENDDO
    call nullify_LSTAUXITEM(LSTAUX)
    LSTAUX%Mat2dim => Mat2dim
    MATELMZERO = .FALSE.
    IF(DIMI.EQ.1.AND.DIMJ.EQ.2)THEN
       call BuildFromlstensor_aux(TENSOR,DIMI,DIMJ,0,MATELMZERO,LSTAUX,TwoDim12_FuncLstSub)
    ELSEIF(DIMI.EQ.1.AND.DIMJ.EQ.3)THEN
       call BuildFromlstensor_aux(TENSOR,DIMI,DIMJ,0,MATELMZERO,LSTAUX,TwoDim13_FuncLstSub)
    ELSE
       call BuildFromlstensor_aux(TENSOR,DIMI,DIMJ,0,MATELMZERO,LSTAUX,TwoDim_FuncLstSub)
    ENDIF
  END SUBROUTINE Build_full_2dim_from_lstensor

subroutine BuildFromlstensor_aux(TENSOR,DIMI,DIMJ,nROW,MATELMZERO,LSTAUX,FuncLstSub)
  implicit none
  TYPE(LSTAUXITEM) ::LSTAUX
  integer :: DIMI,DIMJ,nROW
  TYPE(LSTENSOR)   :: TENSOR
  logical :: MATELMZERO
  EXTERNAL FuncLstSub
  !
  integer :: I,IELM,IATOM,JATOM,KATOM,LATOM,maxBat,maxAng
  integer :: nbat1,nbat2,nbat3,nbat4,imat,MATINDEX(4),DIM
  integer :: IBAT,JBAT,KBAT,LBAT,offset1,offset2,offset3,offset4
  integer :: IANGMOM,JANGMOM,KANGMOM,LANGMOM,nOrbD,sD,nOrbC,sC
  integer :: LORB,KORB,JORB,IORB,sA,sB,nOrbA,nOrbB,nelms,nmat
  TYPE(LSAOTENSOR),pointer    :: lsao
  REAL(REALK),pointer     :: elms(:)
  !READY FOR OMP 
  MATINDEX=0
  IF(MATELMZERO)THEN
     nmat = TENSOR%ndim5/2
  ELSE
     nmat = TENSOR%ndim5
  ENDIF
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(I,lsao,elms,DIM,nelms,IELM,IATOM,&
!$OMP JATOM,KATOM,LATOM,maxBat,maxAng,nbat1,nbat2,nbat3,nbat4,IMAT,Lbat,&
!$OMP offset4,Langmom,nOrbD,sD,LORB,MATINDEX,Kbat,offset3,Kangmom,nOrbC,sC,&
!$OMP Korb,Jbat,offset2,Jangmom,nOrbB,sB,Jorb,Ibat,offset1,Iangmom,nOrbA,sA,&
!$OMP IORB) SHARED(MATELMZERO,LSTAUX,nROW,DIMI,DIMJ) DEFAULT(SHARED)
  DO I = 1,TENSOR%nLSAO
     lsao => TENSOR%LSAO(I) 
     elms => lsao%elms
     DIM = size(lsao%elms)
     nelms = lsao%nelms
     IELM = 0
     IATOM = lsao%ATOM(1)
     JATOM = lsao%ATOM(2)
     KATOM = lsao%ATOM(3)
     LATOM = lsao%ATOM(4)
     maxBat = lsao%maxBat
     maxAng = lsao%maxAng
     IF(maxAng.EQ.1)THEN 
        !nofamily basisset means:
        !no loop over angmom 
        !maxang hardcoded to 1
        nbat1 = TENSOR%nAOBATCH(IATOM,1)
        nbat2 = TENSOR%nAOBATCH(JATOM,2)
        nbat3 = TENSOR%nAOBATCH(KATOM,3)
        nbat4 = TENSOR%nAOBATCH(LATOM,4)
        DO IMAT = 1,nmat
         IF(MATELMZERO)IELM = 0
         DO Lbat = 1,nbat4
          offset4 = (Lbat-1)+3*maxBat
          nOrbD = lsao%nOrb(1+offset4)
          sD = lsao%startGlobalOrb(1+offset4)
          DO LORB = 0,nOrbD-1
           MATINDEX(4)=sD+LORB
           DO Kbat = 1,nbat3
            offset3 = (Kbat-1)+2*maxBat
            nOrbC = lsao%nOrb(1+offset3)
            sC = lsao%startGlobalOrb(1+offset3)    
            DO KORB = 0,nOrbC-1
             MATINDEX(3)=sC+KORB
             DO Jbat = 1,nbat2
              offset2 = (Jbat-1)+maxBat
              nOrbB = lsao%nOrb(1+offset2)
              sB = lsao%startGlobalOrb(1+offset2)
              DO JORB = 0,nOrbB-1
               MATINDEX(2)=sB+JORB
               DO Ibat = 1,nbat1
                nOrbA = lsao%nOrb(Ibat)
                sA = lsao%startGlobalOrb(Ibat)
                call FuncLstSub(nOrbA,sA,sB+JORB,sC+KORB,sD+LORB,IMAT,IELM,elms,dim,&
                     & DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
     ELSE
        nbat1 = TENSOR%nAOBATCH(IATOM,1)
        nbat2 = TENSOR%nAOBATCH(JATOM,2)
        nbat3 = TENSOR%nAOBATCH(KATOM,3)
        nbat4 = TENSOR%nAOBATCH(LATOM,4)
        DO IMAT = 1,nmat
         IF(MATELMZERO)IELM = 0
         DO Lbat = 1,nbat4
          offset4 = (Lbat-1)*maxAng+3*maxBat*maxAng
          DO Langmom = 1,lsao%nAngmom(Lbat+3*maxBat)
           nOrbD = lsao%nOrb(Langmom+offset4)
           sD = lsao%startGlobalOrb(Langmom+offset4)
           DO LORB = 0,nOrbD-1
            MATINDEX(4)=sD+LORB
            DO Kbat = 1,nbat3
             offset3 = (Kbat-1)*maxAng+2*maxBat*maxAng
             DO Kangmom = 1,lsao%nAngmom(Kbat+2*maxBat)
              nOrbC = lsao%nOrb(Kangmom+offset3)
              sC = lsao%startGlobalOrb(Kangmom+offset3)    
              DO KORB = 0,nOrbC-1
               MATINDEX(3)=sC+KORB
               DO Jbat = 1,nbat2
                offset2 = (Jbat-1)*maxAng+maxBat*maxAng
                DO Jangmom = 1,lsao%nAngmom(Jbat+maxBat)
                 nOrbB = lsao%nOrb(Jangmom+offset2)
                 sB = lsao%startGlobalOrb(Jangmom+offset2)
                 DO JORB = 0,nOrbB-1
                  MATINDEX(2)=sB+JORB
                  DO Ibat = 1,nbat1
                   offset1 = (Ibat-1)*maxAng
                   DO Iangmom = 1,lsao%nAngmom(Ibat)
                    nOrbA = lsao%nOrb(Iangmom+offset1)
                    sA = lsao%startGlobalOrb(Iangmom+offset1)
                    call FuncLstSub(nOrbA,sA,sB+JORB,sC+KORB,sD+LORB,IMAT,IELM,elms,dim,&
                        & DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
                   ENDDO
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
     ENDIF
  ENDDO
!$OMP END PARALLEL DO
END SUBROUTINE BuildFromlstensor_aux

subroutine nullify_LSTAUXITEM(LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) ::LSTAUX
  nullify(LSTAUX%sMAT)
  nullify(LSTAUX%aMAT)
  nullify(LSTAUX%COL1)
  nullify(LSTAUX%ROW1)
  nullify(LSTAUX%COL2)
  nullify(LSTAUX%ROW2)
  nullify(LSTAUX%NNZ1)
  nullify(LSTAUX%NNZ2)
  nullify(LSTAUX%VAL1)
  nullify(LSTAUX%VAL2)
  nullify(LSTAUX%MAT1dim)
  nullify(LSTAUX%MAT2dim)
  nullify(LSTAUX%MAT3dim)
  nullify(LSTAUX%MAT4dim)
  nullify(LSTAUX%MAT5dim)
end subroutine nullify_LSTAUXITEM

subroutine FiveDim_FuncLstSub(nOrbA,sA,IB,IC,ID,IMAT,IELM,elms,&
     & DIM,DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) ::LSTAUX
  integer,intent(in) :: sA,nOrbA,IB,IC,ID,IMAT,DIM,DIMI,DIMJ,nROW,nelms
  integer,intent(inout) :: IELM,MATINDEX(4)
  REAL(REALK),intent(in):: elms(DIM)
  !
  integer :: IORB
  DO IORB = 0,nOrbA-1
     IELM = IELM+1
     LSTAUX%MAT5dim(sA+IORB,IB,IC,ID,IMAT)=elms(IELM)
  ENDDO
end subroutine FiveDim_FuncLstSub

subroutine FourDim_FuncLstSub(nOrbA,sA,IB,IC,ID,IMAT,IELM,elms,&
     & DIM,DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) ::LSTAUX
  integer,intent(in) :: sA,nOrbA,IB,IC,ID,IMAT,DIM,DIMI,DIMJ,nROW,nelms
  integer,intent(inout) :: IELM,MATINDEX(4)
  REAL(REALK),intent(in):: elms(DIM)
  !
  integer :: IORB
  DO IORB = 0,nOrbA-1
     IELM = IELM+1
     LSTAUX%MAT4dim(sA+IORB,IB,IC,ID)=elms(IELM)
  ENDDO
end subroutine FourDim_FuncLstSub

subroutine ThreeDim12_FuncLstSub(nOrbA,sA,IB,IC,ID,IMAT,IELM,elms,&
     & DIM,DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) ::LSTAUX
  integer,intent(in) :: sA,nOrbA,IB,IC,ID,IMAT,DIM,DIMI,DIMJ,nROW,nelms
  integer,intent(inout) :: IELM,MATINDEX(4)
  REAL(REALK),intent(in):: elms(DIM)
  !
  integer :: IORB
  DO IORB = 0,nOrbA-1
     IELM = IELM+1
     LSTAUX%MAT3dim(sA+IORB,IB,IMAT)=elms(IELM)
  ENDDO
end subroutine ThreeDim12_FuncLstSub

subroutine ThreeDim13_FuncLstSub(nOrbA,sA,IB,IC,ID,IMAT,IELM,elms,&
     & DIM,DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) ::LSTAUX
  integer,intent(in) :: sA,nOrbA,IB,IC,ID,IMAT,DIM,DIMI,DIMJ,nROW,nelms
  integer,intent(inout) :: IELM,MATINDEX(4)
  REAL(REALK),intent(in):: elms(DIM)
  !
  integer :: IORB
  DO IORB = 0,nOrbA-1
     IELM = IELM+1
     LSTAUX%MAT3dim(sA+IORB,IC,IMAT)=elms(IELM)
  ENDDO
end subroutine ThreeDim13_FuncLstSub

subroutine ThreeDim_FuncLstSub(nOrbA,sA,IB,IC,ID,IMAT,IELM,elms,&
     & DIM,DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) ::LSTAUX
  integer,intent(in) :: sA,nOrbA,IB,IC,ID,IMAT,DIM,DIMI,DIMJ,nROW,nelms
  integer,intent(inout) :: IELM,MATINDEX(4)
  REAL(REALK),intent(in):: elms(DIM)
  !
  integer :: IORB
  DO IORB = 0,nOrbA-1
     IELM = IELM+1
     MATINDEX(1)=sA+IORB
     LSTAUX%MAT3dim(MATINDEX(DIMI),MATINDEX(DIMJ),IMAT)=elms(IELM)
  ENDDO
end subroutine ThreeDim_FuncLstSub

subroutine twoDim12_FuncLstSub(nOrbA,sA,IB,IC,ID,IMAT,IELM,elms,&
     & DIM,DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) ::LSTAUX
  integer,intent(in) :: sA,nOrbA,IB,IC,ID,IMAT,DIM,DIMI,DIMJ,nROW,nelms
  integer,intent(inout) :: IELM,MATINDEX(4)
  REAL(REALK),intent(in):: elms(DIM)
  !
  integer :: IORB
  DO IORB = 0,nOrbA-1
     IELM = IELM+1
     LSTAUX%MAT2dim(sA+IORB,IB)=elms(IELM)
  ENDDO
end subroutine twoDim12_FuncLstSub

subroutine twoDim13_FuncLstSub(nOrbA,sA,IB,IC,ID,IMAT,IELM,elms,&
     & DIM,DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) ::LSTAUX
  integer,intent(in) :: sA,nOrbA,IB,IC,ID,IMAT,DIM,DIMI,DIMJ,nROW,nelms
  integer,intent(inout) :: IELM,MATINDEX(4)
  REAL(REALK),intent(in):: elms(DIM)
  !
  integer :: IORB
  DO IORB = 0,nOrbA-1
     IELM = IELM+1
     LSTAUX%MAT2dim(sA+IORB,IC)=elms(IELM)
  ENDDO
end subroutine twoDim13_FuncLstSub

subroutine twoDim_FuncLstSub(nOrbA,sA,IB,IC,ID,IMAT,IELM,elms,&
     & DIM,DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) ::LSTAUX
  integer,intent(in) :: sA,nOrbA,IB,IC,ID,IMAT,DIM,DIMI,DIMJ,nROW,nelms
  integer,intent(inout) :: IELM,MATINDEX(4)
  REAL(REALK),intent(in):: elms(DIM)
  !
  integer :: IORB
  DO IORB = 0,nOrbA-1
     IELM = IELM+1
     MATINDEX(1)=sA+IORB
     LSTAUX%MAT2dim(MATINDEX(DIMI),MATINDEX(DIMJ))=elms(IELM)
  ENDDO
end subroutine twoDim_FuncLstSub

subroutine SMat_FuncLstSub(nOrbA,sA,IB,IC,ID,IMAT,IELM,elms,&
     & DIM,DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) ::LSTAUX
  integer,intent(in) :: sA,nOrbA,IB,IC,ID,IMAT,DIM,DIMI,DIMJ,nROW,nelms
  integer,intent(inout) :: IELM,MATINDEX(4)
  REAL(REALK),intent(in):: elms(DIM)
  !
  integer :: IORB
  DO IORB = 0,nOrbA-1
     IELM = IELM+1
     MATINDEX(1)=sA+IORB
     LSTAUX%sMAT%elms(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*nROW)=elms(IELM)
  ENDDO
end subroutine SMat_FuncLstSub

subroutine aMat_FuncLstSub(nOrbA,sA,IB,IC,ID,IMAT,IELM,elms,&
     & DIM,DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) ::LSTAUX
  integer,intent(in) :: sA,nOrbA,IB,IC,ID,IMAT,DIM,DIMI,DIMJ,nROW,nelms
  integer,intent(inout) :: IELM,MATINDEX(4)
  REAL(REALK),intent(in):: elms(DIM)
  !
  integer :: IORB
  DO IORB = 0,nOrbA-1
     IELM = IELM+1
     MATINDEX(1)=sA+IORB
     LSTAUX%aMAT(IMAT)%p%elms(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*nROW)=elms(IELM)
  ENDDO
end subroutine aMat_FuncLstSub

subroutine FromaMat_FuncLstSub(nOrbA,sA,IB,IC,ID,IMAT,IELM,elms,&
     & DIM,DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) ::LSTAUX
  integer,intent(in) :: sA,nOrbA,IB,IC,ID,IMAT,DIM,DIMI,DIMJ,nROW,nelms
  integer,intent(inout) :: IELM,MATINDEX(4)
  REAL(REALK),intent(inout):: elms(DIM)
  !
  integer :: IORB
  DO IORB = 0,nOrbA-1
     IELM = IELM+1
     MATINDEX(1)=sA+IORB
     elms(IELM) = LSTAUX%aMAT(IMAT)%p%elms(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*nROW)
  ENDDO
end subroutine FromaMat_FuncLstSub

subroutine nnz1_FuncLstSub(nOrbA,sA,IB,IC,ID,IMAT,IELM,elms,&
     & DIM,DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) ::LSTAUX
  integer,intent(in) :: sA,nOrbA,IB,IC,ID,IMAT,DIM,DIMI,DIMJ,nROW,nelms
  integer,intent(inout) :: IELM,MATINDEX(4)
  REAL(REALK),intent(in):: elms(DIM)
  !
  integer :: IORB
!$OMP CRITICAL
  DO IORB = 0,nOrbA-1
     IELM = IELM+1
     IF(ABS(elms(IELM)).GT.zeroCSR)THEN
        LSTAUX%NNZ1=LSTAUX%NNZ1+1
     ENDIF
  ENDDO
!$OMP END CRITICAL
end subroutine nnz1_FuncLstSub

subroutine CSR1_FuncLstSub(nOrbA,sA,IB,IC,ID,IMAT,IELM,elms,&
     & DIM,DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) :: LSTAUX
  integer,intent(in) :: sA,nOrbA,IB,IC,ID,IMAT,DIM,DIMI,DIMJ,nROW,nelms
  integer,intent(inout) :: IELM,MATINDEX(4)
  REAL(REALK),intent(in):: elms(DIM)
  !
  integer :: IORB
!$OMP CRITICAL
  DO IORB = 0,nOrbA-1
     IELM = IELM+1
     IF(ABS(elms(IELM)).GT.zeroCSR)THEN
        MATINDEX(1)=sA+IORB
        LSTAUX%NNZ1=LSTAUX%NNZ1+1
        LSTAUX%VAL1(LSTAUX%NNZ1) = elms(IELM) 
        LSTAUX%ROW1(LSTAUX%NNZ1) = MATINDEX(DIMI) 
        LSTAUX%COL1(LSTAUX%NNZ1) = MATINDEX(DIMJ) 
     ENDIF
  ENDDO
!$OMP END CRITICAL
end subroutine CSR1_FuncLstSub

subroutine nnz2_FuncLstSub(nOrbA,sA,IB,IC,ID,IMAT,IELM,elms,&
     & DIM,DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) ::LSTAUX
  integer,intent(in) :: sA,nOrbA,IB,IC,ID,IMAT,DIM,DIMI,DIMJ,nROW,nelms
  integer,intent(inout) :: IELM,MATINDEX(4)
  REAL(REALK),intent(in):: elms(DIM)
  !
  integer :: IORB
!$OMP CRITICAL
  DO IORB = 0,nOrbA-1
     IELM = IELM+1
     IF(ABS(elms(IELM)).GT.zeroCSR)THEN
        LSTAUX%NNZ2(IMAT)=LSTAUX%NNZ2(IMAT)+1
     ENDIF
  ENDDO
!$OMP END CRITICAL
end subroutine nnz2_FuncLstSub

subroutine CSR2_FuncLstSub(nOrbA,sA,IB,IC,ID,IMAT,IELM,elms,&
     & DIM,DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) :: LSTAUX
  integer,intent(in) :: sA,nOrbA,IB,IC,ID,IMAT,DIM,DIMI,DIMJ,nROW,nelms
  integer,intent(inout) :: IELM,MATINDEX(4)
  REAL(REALK),intent(in):: elms(DIM)
  !
  integer :: IORB
!$OMP CRITICAL
  DO IORB = 0,nOrbA-1
     IELM = IELM+1
     IF(ABS(elms(IELM)).GT.zeroCSR)THEN
        MATINDEX(1)=sA+IORB
        LSTAUX%NNZ2(IMAT)=LSTAUX%NNZ2(IMAT)+1
        LSTAUX%VAL2(LSTAUX%NNZ2(IMAT),IMAT) = elms(IELM) 
        LSTAUX%ROW2(LSTAUX%NNZ2(IMAT),IMAT) = MATINDEX(DIMI) 
        LSTAUX%COL2(LSTAUX%NNZ2(IMAT),IMAT) = MATINDEX(DIMJ) 
     ENDIF
  ENDDO
!$OMP END CRITICAL
end subroutine CSR2_FuncLstSub

subroutine opt1UnresMatS_FuncLstSub(nOrbA,sA,IB,IC,ID,IMAT,IELM,elms,&
     & DIM,DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) :: LSTAUX
  integer,intent(in) :: sA,nOrbA,IB,IC,ID,IMAT,DIM,DIMI,DIMJ,nROW,nelms
  integer,intent(inout) :: IELM,MATINDEX(4)
  REAL(REALK),intent(in):: elms(DIM)
  !
  integer :: IORB
  DO IORB = 0,nOrbA-1
     IELM = IELM+1
     MATINDEX(1)=sA+IORB
     LSTAUX%sMAT%elms(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*nrow)=elms(IELM)
     LSTAUX%sMAT%elmsb(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*nrow)=elms(IELM)
  ENDDO
end subroutine opt1UnresMatS_FuncLstSub

subroutine opt2UnresMatS_FuncLstSub(nOrbA,sA,IB,IC,ID,IMAT,IELM,elms,&
     & DIM,DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) :: LSTAUX
  integer,intent(in) :: sA,nOrbA,IB,IC,ID,IMAT,DIM,DIMI,DIMJ,nROW,nelms
  integer,intent(inout) :: IELM,MATINDEX(4)
  REAL(REALK),intent(in):: elms(DIM)
  !
  integer :: IORB
  DO IORB = 0,nOrbA-1
     IELM = IELM+1
     MATINDEX(1)=sA+IORB
     LSTAUX%sMAT%elms(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*nrow)=elms(IELM)
     LSTAUX%sMAT%elmsb(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*nrow)=elms(IELM+nElms)
  ENDDO
end subroutine opt2UnresMatS_FuncLstSub

! Coulomb type contributions
subroutine opt1UnresMatA_FuncLstSub(nOrbA,sA,IB,IC,ID,IMAT,IELM,elms,&
     & DIM,DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) :: LSTAUX
  integer,intent(in) :: sA,nOrbA,IB,IC,ID,IMAT,DIM,DIMI,DIMJ,nROW,nelms
  integer,intent(inout) :: IELM,MATINDEX(4)
  REAL(REALK),intent(in):: elms(DIM)
  !
  integer :: IORB
  DO IORB = 0,nOrbA-1
     IELM = IELM+1
     MATINDEX(1)=sA+IORB
     LSTAUX%aMAT(IMAT)%p%elms(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*nrow)=elms(IELM)
     LSTAUX%aMAT(IMAT)%p%elmsb(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*nrow)=elms(IELM)
  ENDDO
end subroutine opt1UnresMatA_FuncLstSub

! Exchange type contributions
subroutine opt2UnresMatA_FuncLstSub(nOrbA,sA,IB,IC,ID,IMAT,IELM,elms,&
     & DIM,DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) :: LSTAUX
  integer,intent(in) :: sA,nOrbA,IB,IC,ID,IMAT,DIM,DIMI,DIMJ,nROW,nelms
  integer,intent(inout) :: IELM,MATINDEX(4)
  REAL(REALK),intent(in):: elms(DIM)
  !
  integer :: IORB
  DO IORB = 0,nOrbA-1
     IELM = IELM+1
     MATINDEX(1)=sA+IORB
     LSTAUX%aMAT(IMAT)%p%elms(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*nrow)=elms(IELM+nElms*2*(IMAT-1))
     LSTAUX%aMAT(IMAT)%p%elmsb(MATINDEX(DIMI)+(MATINDEX(DIMJ)-1)*nrow)=elms(IELM+nElms*(2*IMAT-1))
  ENDDO
end subroutine opt2UnresMatA_FuncLstSub

subroutine add2dim_FuncLstSub(nOrbA,sA,IB,IC,ID,IMAT,IELM,elms,&
     & DIM,DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) :: LSTAUX
  integer,intent(in) :: sA,nOrbA,IB,IC,ID,IMAT,DIM,DIMI,DIMJ,nROW,nelms
  integer,intent(inout) :: IELM,MATINDEX(4)
  REAL(REALK),intent(inout):: elms(DIM)
  !
  integer :: IORB
  DO IORB = 0,nOrbA-1
     IELM = IELM+1
     MATINDEX(1)=sA+IORB
     elms(IELM) = elms(IELM) + LSTAUX%MAT2dim(MATINDEX(DIMI),MATINDEX(DIMJ))
  ENDDO
end subroutine add2dim_FuncLstSub

subroutine elmsFrom5dim_FuncLstSub(nOrbA,sA,IB,IC,ID,IMAT,IELM,elms,&
     & DIM,DIMI,DIMJ,nROW,MATINDEX,nelms,LSTAUX)
  implicit none
  TYPE(LSTAUXITEM) :: LSTAUX
  integer,intent(in) :: sA,nOrbA,IB,IC,ID,IMAT,DIM,DIMI,DIMJ,nROW,nelms
  integer,intent(inout) :: IELM,MATINDEX(4)
  REAL(REALK),intent(inout):: elms(DIM)
  !
  integer :: IORB
  DO IORB = 0,nOrbA-1
     IELM = IELM+1
     MATINDEX(1)=sA+IORB
     elms(IELM) = LSTAUX%MAT5dim(sA+IORB,IB,IC,ID,IMAT)
  ENDDO
end subroutine elmsFrom5dim_FuncLstSub

  !> \brief build full 2 dimensional array from an lstensor 
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param TENSOR the lstensor
  !> \param MAT134 the 2 dimensional array
  !> \param dim1 the size of the 1. dimension
  !> \param dim2 the size of the 2. dimension
  SUBROUTINE Build_full_shortint2dim_from_lstensor(TENSOR,MAT134,dim1,dim2,lupri)
    implicit none
    integer          :: dim1,dim2,lupri
    TYPE(LSTENSOR)   :: TENSOR
    integer(kind=short)      :: MAT134(dim1,dim2)
    !
    INTEGER    :: I,J,nbatches1,nbatches2,ndim5,jbatch,ibatch,nbatI,nbatJ
    INTEGER    :: ibat,jbat,s1,s2,n1,n2,nA,nB,nAngB,nAngA,iAng1,iAng2
    INTEGER    :: IATOM,JATOM,nbat1,nbat2,offset1,offset2,maxBat,maxAng
    LOGICAL    :: ERROR
    type(lsaotensor),pointer :: lsao
    integer(kind=long)  :: nmemsize,nsize
    integer(kind=short) :: maxgab,maxgab2
    ERROR=.FALSE.
    IF(dim1.NE.TENSOR%nbatches(1))ERROR=.TRUE.
    IF(dim2.NE.TENSOR%nbatches(2))ERROR=.TRUE.
    IF(ERROR)THEN
       print*,'ASSOCIATED(TENSOR%maxgab)',ASSOCIATED(TENSOR%maxgab)
       print*,'LSTENSOR%nbatches1',TENSOR%nbatches(1),'dim1',dim1
       print*,'LSTENSOR%nbatches2',TENSOR%nbatches(2),'dim2',dim2
       CALL lsQUIT('ERR:Build_full_shortint2dim_from_lstensor dimmismatch',-1)
    ENDIF
    IF(TENSOR%ScreenTensor)THEN
       DO J = 1,dim2
          DO I = 1,dim1
             MAT134(I,J) = TENSOR%maxgab(I,J)
          ENDDO
       ENDDO
    ELSE
       IF(.NOT.ASSOCIATED(TENSOR%maxgab))THEN
          nbatches1 = TENSOR%nbatches(1)
          nbatches2 = TENSOR%nbatches(2)
          NULLIFY(TENSOR%maxgab)
          CALL MEM_ALLOC(TENSOR%maxgab,nbatches1,nbatches2)
          nsize = size(TENSOR%maxgab,KIND=long)*mem_shortintsize
          call mem_allocated_mem_lstensor(nsize)
          call ls_sizero(TENSOR%maxgab,nbatches1*nbatches2)
          ndim5 = TENSOR%ndim5
          DO Jatom = 1,TENSOR%natom(2)
             nBatJ = TENSOR%nAOBATCH(JATOM,2)
             DO Iatom = 1,TENSOR%natom(1)       
                I = TENSOR%INDEX(IATOM,JATOM,1,1)
                IF(I.GT.0)THEN
                   lsao => TENSOR%LSAO(I)
                   Ibatch = lsao%AOBATCH(1)
                   Jbatch = lsao%AOBATCH(2)
                   maxBat = lsao%maxBat
                   maxAng = lsao%maxAng
                   nBatI = TENSOR%nAOBATCH(IATOM,1)

                   n1=TENSOR%LSAO(I)%nLocal(1)
                   n2=TENSOR%LSAO(I)%nLocal(2)
                   do jbat = 1,nBatJ
                      nAngB=TENSOR%LSAO(I)%nAngmom(jbat+maxBat)
                      offset2 = (jbat-1)*maxAng+maxAng*maxBat
                      do ibat = 1,nBatI
                         nAngA=TENSOR%LSAO(I)%nAngmom(ibat)
                         offset1 = (ibat-1)*maxAng
                         maxgab2 = shortzero  
                         do iAng2=1,nAngB
                            s2=TENSOR%LSAO(I)%startLocalOrb(iAng2+offset2)-1
                            nB=TENSOR%LSAO(I)%nOrb(iAng2+offset2)                      
                            do iAng1=1,nAngA
                               s1=TENSOR%LSAO(I)%startLocalOrb(iAng1+offset1)-1
                               nA=TENSOR%LSAO(I)%nOrb(iAng1+offset1)
                               call determine_maxabsval(maxgab,lsao%elms,n1,n2,s1,s2,nA,nB,ndim5,lupri)
                               IF(maxgab.GT.maxgab2)maxgab2 = maxgab  
                            enddo
                         enddo
                         TENSOR%maxgab(Ibatch + ibat,Jbatch + jbat) = maxgab2
                      enddo
                   enddo
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       DO J = 1,dim2
          DO I = 1,dim1
             MAT134(I,J) = TENSOR%maxgab(I,J)
          ENDDO
       ENDDO
    ENDIF
  END SUBROUTINE Build_full_shortint2dim_from_lstensor

  !> \brief build Atomic Gab matrix
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param TENSOR the lstensor
  !> \param MAT134 the 2 dimensional array
  !> \param dim1 the size of the 1. dimension
  !> \param dim2 the size of the 2. dimension
  SUBROUTINE lstensor_BuildAtomicGab(AtomGab,nAtom1,nAtom2,TENSOR,AO1,AO2,lupri)
    implicit none
    integer              :: nAtom1,nAtom2,lupri
    TYPE(LSTENSOR)       :: TENSOR
    integer(kind=short)  :: AtomGab(nAtom1,nAtom2)
    TYPE(AOITEM)         :: AO1,AO2
    !
    INTEGER    :: I,J,nbatches1,nbatches2,jbatch,ibatch,nbatI,nbatJ
    INTEGER    :: ibat,jbat
    INTEGER    :: IATOM,JATOM,nbat1,nbat2
    LOGICAL    :: ERROR
    integer(kind=short) :: maxgab,maxgab2

    ERROR=.FALSE.
    IF(nAtom1.NE.TENSOR%nAtom(1))ERROR=.TRUE.
    IF(nAtom2.NE.TENSOR%nAtom(2))ERROR=.TRUE.
    IF(ERROR)THEN
       print*,'ASSOCIATED(TENSOR%maxgab)',ASSOCIATED(TENSOR%maxgab)
       print*,'LSTENSOR%nAtom1',TENSOR%nAtom(1),'nAtom1',nAtom1
       print*,'LSTENSOR%nAtom2',TENSOR%nAtom(2),'nAtom2',nAtom2
       CALL lsQUIT('ERR: lstensor_BuildAtomicGab dim mismatch',-1)
    ENDIF
    IF(TENSOR%ScreenTensor)THEN
       DO Jatom = 1,TENSOR%natom(2)
          nBatJ = TENSOR%nAOBATCH(JATOM,2)
          DO Iatom = 1,TENSOR%natom(1)       
             maxgab2 = shortzero  
             I = TENSOR%INDEX(IATOM,JATOM,1,1)
             IF(I.GT.0)THEN
                Ibatch = TENSOR%SLSAO(I)%AOBATCH(1)
                Jbatch = TENSOR%SLSAO(I)%AOBATCH(2)                   
                nBatI = TENSOR%nAOBATCH(IATOM,1)
                do jbat = 1,nBatJ
                   do ibat = 1,nBatI
                      maxgab = TENSOR%maxgab(Ibatch + ibat,Jbatch + jbat)
                      IF(maxgab.GT.maxgab2) maxgab2 = maxgab  
                   enddo
                enddo
             ENDIF
             AtomGab(IAtom,Jatom) = maxgab2
          ENDDO
       ENDDO
    ENDIF
  END SUBROUTINE Lstensor_BuildAtomicGab

  subroutine determine_maxabsval(maxgab,elms,n1,n2,s1,s2,nA,nB,ndim5,lupri)
    implicit none
    integer(kind=short) :: maxgab
    integer     :: n1,n2,s1,s2,nA,nB,ndim5,lupri
    real(realk) :: elms(n1,n2,ndim5)
    !
    real(realk) :: maxElm
    integer :: IB,IA,iDmat
    maxElm = 0E0_realk
    DO IDMAT = 1,ndim5
       DO IB=1,nB
          DO IA=1,nA
             maxElm = MAX(maxElm,ABS(elms(s1+iA,s2+iB,IDMAT)))
          ENDDO
       ENDDO
    ENDDO
    !Beware when converting from double precision to short integer 
    !If double precision is less than 10^-33 then you can run into
    !problems with short integer overflow
    IF(maxElm.GT.shortintCRIT)THEN
       maxgab = CEILING(LOG10(maxElm))
    ELSE
       maxgab = shortzero
    ENDIF
  end subroutine determine_maxabsval

  !> \brief initiate the primitive screening lstensor structure
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param TENSOR the lstensor
  !> \param AO1 the Atomic orbital(AO) item for center 1 
  !> \param AO2 the AO item for center 2 
  !> \param AO3 the AO item for center 3 
  !> \param AO4 the AO item for center 4 
  !> \param nbast1 the size of the 1. dimension
  !> \param nbast2 the size of the 2. dimension
  !> \param nbast3 the size of the 3. dimension
  !> \param nbast4 the size of the 4. dimension
  !> \param nmat the number of density matrices or size of the 5. dimension
  !> \param useAO1 flag to describe if the AO1 item should be used or an empty AO 
  !> \param useAO2 flag to describe if the AO1 item should be used or an empty AO 
  !> \param useAO3 flag to describe if the AO1 item should be used or an empty AO 
  !> \param useAO4 flag to describe if the AO1 item should be used or an empty AO 
  !> \param ODscreen flag: use the overlap distribution screening 
  !> \param lupri the logical unit number for the output
  SUBROUTINE init_screening_lstensor(TENSOR,AO1,AO2,nbast1,nbast2,ODscreen,lupri,nonscreenAB)
    implicit none
    INTEGER            :: nbast1,nbast2,lupri
    TYPE(LSTENSOR)     :: TENSOR
    TYPE(AOITEM),target:: AO1,AO2
    LOGICAL            :: ODscreen
    LOGICAL,pointer,optional :: nonscreenAB(:,:)
    !
    TYPE(AOITEMPOINTER) :: AOT(2)
    ! 
    INTEGER :: natom(2),nbat(2),AOTbatch(2),batAOT(2),nAng(2)
    !
    LOGICAL,pointer :: nonscreenAB2(:,:)
    INTEGER :: nbastI,Ibat,Iangmom,Iatom,IORB,nOrbA,sA
    INTEGER :: nbastJ,Jbat,Jangmom,Jatom,JORB,nOrbB,sB
    INTEGER :: nbastK,Kbat,Kangmom,Katom,KORB,nOrbC,sC
    INTEGER :: nbastL,Lbat,Langmom,Latom,LORB,nOrbD,sD,lstmem_index
    INTEGER :: batAOT1,batAOT2,nbatJ,nbatI,maxnatom,offset
    ! COMMON
    INTEGER :: DIM,IELM,IMAT2,I,J,iAO,jAO,maxBat,maxAng,nLocal(2)
    INTEGER :: nElms,startLocalOrb,nAngG,batAOTG,nbatches1,nbatches2
    INTEGER :: batAOTGJ,batAOTGI
    LOGICAL :: ERROR,Screen,ScreenAtom,OVERALLscreen
    integer(kind=long) :: nmemsize,AllocInt,AllocRealk,AllocIntS,nsize
    TYPE(SLSAOTENSOR),pointer    :: slsao

    AOT(1)%p => AO1
    AOT(2)%p => AO2
    DO IAO = 1,2
       natom(iAO) = AOT(iAO)%p%natoms
    ENDDO

    IF(ODscreen.AND..NOT.present(nonscreenAB))THEN
       CALL MEM_ALLOC(nonscreenAB2,AO1%natoms,AO2%natoms)
       IF(ODscreen)THEN
          call build_atomicODScreen(AO1,AO2,nonscreenAB2,AO1%natoms,AO2%natoms)
       else
          nonscreenAB2 = .true.
       ENDIF
    ELSEIF(ODscreen)THEN !nonscreenAB present
       nonscreenAB2 => nonscreenAB
    ELSE
       CALL MEM_ALLOC(nonscreenAB2,AO1%natoms,AO2%natoms)
       nonscreenAB2 = .true.
    ENDIF


    I = 0
    DO Jatom = 1,natom(2)
       DO Iatom = 1,natom(1)
          IF(nonscreenAB2(Iatom,Jatom)) I=I+1
       ENDDO
    ENDDO
    TENSOR%nSLSAO = I
    TENSOR%pChargetensor = .FALSE.
    TENSOR%Econtrib = .FALSE.
    TENSOR%Screentensor = .TRUE.
    TENSOR%Screenoutput = .FALSE.
    TENSOR%MagGradienttensor = .FALSE.
    TENSOR%Gradienttensor = .FALSE.
    NULLIFY(TENSOR%SLSAO)
    CALL MEM_ALLOC(TENSOR%SLSAO,TENSOR%nSLSAO)
    NULLIFY(TENSOR%INDEX)
    CALL MEM_ALLOC(TENSOR%INDEX,natom(1),natom(2),1,1)
    nsize = size(TENSOR%INDEX,KIND=long)*mem_intsize
    call mem_allocated_mem_lstensor(nsize)
    DO Jatom = 1,natom(2)
       DO Iatom = 1,natom(1)
          TENSOR%INDEX(Iatom,Jatom,1,1) = 0 !if 0 lsaotensor not call mem_allocd 
       ENDDO
    ENDDO
    DO IAO = 1,2
       TENSOR%natom(iAO) = natom(iAO)
    ENDDO
    TENSOR%natom(3) = 1
    TENSOR%natom(4) = 1
    TENSOR%ndim5 = 1

    maxnatom = MAX(TENSOR%natom(1),TENSOR%natom(2))
    nullify(TENSOR%nAOBATCH)
    call mem_alloc(TENSOR%nAOBATCH,maxnatom,2)
    nsize = size(TENSOR%nAOBATCH,KIND=long)*mem_intsize
    call mem_allocated_mem_lstensor(nsize)

    DO IATOM = 1,TENSOR%natom(1)
       TENSOR%nAOBATCH(IATOM,1) = AO1%ATOMICnBatch(IATOM)
    ENDDO
    DO IATOM = 1,TENSOR%natom(2)
       TENSOR%nAOBATCH(IATOM,2) = AO2%ATOMICnBatch(IATOM)
    ENDDO

    nbatches1 = AO1%nbatches
    nbatches2 = AO2%nbatches
    TENSOR%nbatches(1) = nbatches1
    TENSOR%nbatches(2) = nbatches2
    TENSOR%nbatches(3) = 1
    TENSOR%nbatches(4) = 1
    TENSOR%nbast(1) = nbast1
    TENSOR%nbast(2) = nbast2
    TENSOR%nbast(3) = 1
    TENSOR%nbast(4) = 1
    !Mem Estimate
    AllocInt = 0
    AllocIntS = 0
    AllocRealk = 0 
    AOTbatch(2)=0
    DO Jatom = 1,natom(2)
       nbat(2) = AOT(2)%p%ATOMICnBatch(JATOM)
       AOTbatch(1)=0
       DO Iatom = 1,natom(1)
          nbat(1) = AOT(1)%p%ATOMICnBatch(IATOM)
          IF(nonscreenAB2(Iatom,Jatom))THEN
             maxBat = nBat(1)
             maxBat = MAX(maxBat,nBat(2))
             AllocInt = AllocInt + 4*maxBat
             DIM = 0
             batAOTGJ = AOTbatch(2)
             DO Jbat = 1,nbat(2)
                batAOTGJ = batAOTGJ+1
                nOrbB = AOT(2)%p%BATCH(batAOTGJ)%nPrimitives
                batAOTGI = AOTbatch(1)
                DO Ibat = 1,nbat(1)
                   batAOTGI = batAOTGI+1
                   nOrbA=AOT(1)%p%BATCH(batAOTGI)%nPrimitives
                   DIM = DIM + nOrbA*nOrbB
                ENDDO
             ENDDO
             AllocInts = AllocInts + DIM
          ENDIF
          AOTbatch(1) = AOTbatch(1) + nbat(1)
       ENDDO
       AOTbatch(2) = AOTbatch(2) + nbat(2)
    ENDDO
    call init_lstensorMem(AllocInt,AllocRealk,AllocInts,lstmem_index)
    call zero_lstensorMem
    TENSOR%lstmem_index = lstmem_index
    I = 0
    AOTbatch(2)=0
    DO Jatom = 1,natom(2)
       nbat(2) = AOT(2)%p%ATOMICnBatch(JATOM)
       AOTbatch(1)=0
       DO Iatom = 1,natom(1)
          nbat(1) = AOT(1)%p%ATOMICnBatch(IATOM)
          IF(nonscreenAB2(Iatom,Jatom))THEN
             I=I+1
             TENSOR%INDEX(IATOM,JATOM,1,1) = I
             slsao => TENSOR%SLSAO(I)
             call slsaotensor_nullify(slsao)
             slsao%ATOM(1) = Iatom
             slsao%ATOM(2) = Jatom
             slsao%AOBATCH(1) = AOTbatch(1)
             slsao%AOBATCH(2) = AOTbatch(2)
             maxBat = nBat(1)
             maxBat = MAX(maxBat,nBat(2))
             slsao%maxBat = maxBat
             NULLIFY(slsao%nOrb)
             NULLIFY(slsao%startLocalOrb)
             CALL MEM_LSTPOINTER_ALLOC(slsao%nOrb,maxBat*2)
             CALL MEM_LSTPOINTER_ALLOC(slsao%startLocalOrb,maxBat*2)
             do IAO=1,2
                nLocal(iAO) = 0
                batAOTG = AOTbatch(IAO)
                startLocalOrb = 1
                offset = (iAO-1)*maxBat
                DO Ibat = 1,nbat(iAO)
                   batAOTG = batAOTG+1
                   slsao%nOrb(IBat+offset)=AOT(iAO)%p%BATCH(batAOTG)%nPrimitives
                   nLocal(iAO) = nLocal(iAO) + slsao%nOrb(IBat+offset)
                   slsao%startLocalOrb(IBat+offset)=startLocalOrb
                   startLocalOrb = startLocalOrb + slsao%nOrb(IBat+offset)
                ENDDO
                slsao%nLocal(iAO) = nLocal(iAO)
             ENDDO
             NULLIFY(slsao%selms)
             slsao%nelms = 0
          ENDIF
          AOTbatch(1) = AOTbatch(1) + nbat(1)
       ENDDO
       AOTbatch(2) = AOTbatch(2) + nbat(2)
    ENDDO
    TENSOR%nSLSAO = I
    IF(ODscreen.AND..NOT.present(nonscreenAB))THEN
       CALL MEM_DEALLOC(nonscreenAB2)
    ELSEIF(ODscreen)THEN 
       !nonscreenAB present use nonscreenAB2 as pointer
    ELSE
       CALL MEM_DEALLOC(nonscreenAB2)
    ENDIF
  END SUBROUTINE INIT_SCREENING_LSTENSOR

  !> \brief initiate the primitive screening lstensor structure
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param TENSOR the lstensor
  !> \param AO1 the Atomic orbital(AO) item for center 1 
  !> \param AO2 the AO item for center 2 
  !> \param AO3 the AO item for center 3 
  !> \param AO4 the AO item for center 4 
  !> \param nbast1 the size of the 1. dimension
  !> \param nbast2 the size of the 2. dimension
  !> \param nbast3 the size of the 3. dimension
  !> \param nbast4 the size of the 4. dimension
  !> \param nmat the number of density matrices or size of the 5. dimension
  !> \param useAO1 flag to describe if the AO1 item should be used or an empty AO 
  !> \param useAO2 flag to describe if the AO1 item should be used or an empty AO 
  !> \param useAO3 flag to describe if the AO1 item should be used or an empty AO 
  !> \param useAO4 flag to describe if the AO1 item should be used or an empty AO 
  !> \param ODscreen flag: use the overlap distribution screening 
  !> \param lupri the logical unit number for the output
  SUBROUTINE init_ps_lstensor(TENSOR,AO1,AO2,nbast1,nbast2,ODscreen,lupri)
    implicit none
    INTEGER            :: nbast1,nbast2,lupri
    TYPE(LSTENSOR)     :: TENSOR
    TYPE(AOITEM),target  :: AO1,AO2
    LOGICAL            :: ODscreen
    !
    TYPE(AOITEMPOINTER) :: AOT(2)
    ! 
    INTEGER :: natom(2),nbat(2),AOTbatch(2),batAOT(2),nAng(2)
    !
    INTEGER :: nbastI,Ibat,Iangmom,Iatom,IORB,nOrbA,sA
    INTEGER :: nbastJ,Jbat,Jangmom,Jatom,JORB,nOrbB,sB
    INTEGER :: nbastK,Kbat,Kangmom,Katom,KORB,nOrbC,sC
    INTEGER :: nbastL,Lbat,Langmom,Latom,LORB,nOrbD,sD
    ! COMMON
    INTEGER :: DIM,IELM,IMAT2,I,J,iAO,jAO,maxBat,maxAng,nLocal(2)
    INTEGER :: nElms,startLocalOrb,nAngG,batAOTG,nbatches1,nbatches2
    INTEGER :: batAOT1,batAOT2,nbatI,nbatJ
    LOGICAL :: ERROR,Screen,ScreenAtom,OVERALLscreen
    LOGICAL,pointer :: nonscreenAB(:,:)
    integer(kind=long) :: nmemsize
    TYPE(SLSAOTENSOR),pointer    :: slsao

    AOT(1)%p => AO1
    AOT(2)%p => AO2
    tensor%screenoutput=.false.
    DO IAO = 1,2
       natom(iAO) = AOT(iAO)%p%natoms
    ENDDO

    CALL MEM_ALLOC(nonscreenAB,AO1%natoms,AO2%natoms)
    IF(ODscreen)THEN
       call build_atomicODScreen(AO1,AO2,nonscreenAB,AO1%natoms,AO2%natoms)
    else
       nonscreenAB = .true.
    ENDIF

    IF(.NOT.TENSOR%Screentensor)THEN
       !not initialised yet
       Call INIT_SCREENING_LSTENSOR(TENSOR,AO1,AO2,nbast1,nbast2,ODscreen,lupri,nonscreenAB)
    ELSE
       !already initialised 
    ENDIF
    nullify(TENSOR%maxprimgab)
    DO Jatom = 1,natom(2)
       nbat(2) = AOT(2)%p%ATOMICnBatch(JATOM)
       DO Iatom = 1,natom(1)
          nbat(1) = AOT(1)%p%ATOMICnBatch(IATOM)
          IF(nonscreenAB(Iatom,Jatom))THEN             
             I=TENSOR%INDEX(Iatom,Jatom,1,1)
             IF(I.NE.0)THEN
                slsao => TENSOR%SLSAO(I)
                maxbat = slsao%maxbat
                DIM = 0
                DO Jbat = 1,nbat(2)
                   nOrbB = slsao%nOrb(Jbat+maxBat)
                   DO Ibat = 1,nbat(1)
                      nOrbA = slsao%nOrb(Ibat)
                      DIM = DIM+nOrbA*nOrbB
                   ENDDO
                ENDDO
                NULLIFY(slsao%selms)
                CALL MEM_LSTPOINTER_ALLOC(slsao%selms,DIM)
                slsao%nelms = DIM
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    CALL MEM_DEALLOC(nonscreenAB)

  END SUBROUTINE INIT_PS_LSTENSOR

  !> \brief initiate the primitive screening lstensor structure
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param TENSOR the lstensor
  !> \param AO1 the Atomic orbital(AO) item for center 1 
  !> \param AO2 the AO item for center 2 
  !> \param AO3 the AO item for center 3 
  !> \param AO4 the AO item for center 4 
  !> \param nbast1 the size of the 1. dimension
  !> \param nbast2 the size of the 2. dimension
  !> \param nbast3 the size of the 3. dimension
  !> \param nbast4 the size of the 4. dimension
  !> \param nmat the number of density matrices or size of the 5. dimension
  !> \param useAO1 flag to describe if the AO1 item should be used or an empty AO 
  !> \param useAO2 flag to describe if the AO1 item should be used or an empty AO 
  !> \param useAO3 flag to describe if the AO1 item should be used or an empty AO 
  !> \param useAO4 flag to describe if the AO1 item should be used or an empty AO 
  !> \param ODscreen flag: use the overlap distribution screening 
  !> \param lupri the logical unit number for the output
  SUBROUTINE init_cs_lstensor(TENSOR,AO1,AO2,nbast1,nbast2,ODscreen,lupri)
    implicit none
    INTEGER            :: nbast1,nbast2,lupri
    LOGICAL            :: ODSCREEN
    TYPE(LSTENSOR)     :: TENSOR
    TYPE(AOITEM),target  :: AO1,AO2
    !
    INTEGER :: I,J,dim,nbatches1,nbatches2
    integer(kind=long) :: nmemsize,nsize
    IF(.NOT.TENSOR%Screentensor)THEN
       !not initialised yet
       Call INIT_SCREENING_LSTENSOR(TENSOR,AO1,AO2,nbast1,nbast2,ODscreen,lupri)
    ELSE
       !already initialised 
    ENDIF
    tensor%screenoutput=.false.
    nbatches1 = TENSOR%nbatches(1)
    nbatches2 = TENSOR%nbatches(2)
    nullify(TENSOR%maxgab)
    call mem_alloc(TENSOR%maxgab,nbatches1,nbatches2)
    nsize = size(TENSOR%maxgab,KIND=long)*mem_shortintsize
    call mem_allocated_mem_lstensor(nsize)

    dim = nbatches2*nbatches1
    call ls_sizero(TENSOR%maxgab,dim)

  END SUBROUTINE INIT_CS_LSTENSOR
  
  !> \brief set the maxgab values
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param TENSOR the lstensor
  !> \param lupri the logical unit number for the output
  SUBROUTINE set_lst_maxgabelms(TENSOR)
    implicit none
    TYPE(LSTENSOR)     :: TENSOR
    !
    integer :: I,J,nbatches1,nbatches2
    integer(kind=short) :: maxGab
    nbatches1 = TENSOR%nbatches(1)
    nbatches2 = TENSOR%nbatches(2)
    maxGab = shortzero
    IF(ASSOCIATED(TENSOR%maxgab))THEN
       DO J = 1,nbatches2
          DO I = 1,nbatches1
             IF(TENSOR%maxgab(I,J).GT.maxgab)THEN
                maxGab = TENSOR%maxgab(I,J)
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    TENSOR%maxGabelm = maxGab    
  END SUBROUTINE Set_lst_maxgabelms

  !> \brief set the maxprimgab values
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param TENSOR the lstensor
  !> \param lupri the logical unit number for the output
  SUBROUTINE set_lst_maxprimgabelms(TENSOR)
    implicit none
    TYPE(LSTENSOR)     :: TENSOR
    !
    integer :: I,J,nbatches1,nbatches2
    integer(kind=short) :: maxGab
    nbatches1 = TENSOR%nbatches(1)
    nbatches2 = TENSOR%nbatches(2)
    maxGab = shortzero
    IF(ASSOCIATED(TENSOR%maxprimgab))THEN
       DO J = 1,nbatches2
          DO I = 1,nbatches1
             IF(TENSOR%maxprimgab(I,J).GT.maxgab)THEN
                maxGab = TENSOR%maxprimgab(I,J)
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    TENSOR%maxprimGabelm = maxGab    
  END SUBROUTINE Set_lst_maxprimgabelms

  !> \brief set the maxprimgab values
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param TENSOR the lstensor
  !> \param lupri the logical unit number for the output
  SUBROUTINE set_lst_maxprimgab(TENSOR)
    implicit none
    TYPE(LSTENSOR)     :: TENSOR
    !
    integer :: I,J,nbatches1,nbatches2
    integer(kind=short) :: maxprimGab
    integer :: Jbatch,Ibatch,Iatom,Jatom,nBatI,nBatJ,jbat,ibat
    integer :: n1,n2,s1,s2,nA,nB,maxbat
    type(slsaotensor),pointer :: slsao
    integer(kind=long) :: nmemsize,nsize

    nbatches1 = TENSOR%nbatches(1)
    nbatches2 = TENSOR%nbatches(2)
    nullify(TENSOR%maxprimgab)
    call mem_alloc(TENSOR%maxprimgab,nbatches1,nbatches2)
    nsize = size(TENSOR%maxprimgab,KIND=long)*mem_shortintsize
    call mem_allocated_mem_lstensor(nsize)

    call ls_sizero(TENSOR%maxprimgab,nbatches1*nbatches2)

!    Jbatch = 0
    DO Jatom = 1,TENSOR%natom(2)
!     Ibatch = 0
     DO Iatom = 1,TENSOR%natom(1)       
      I = TENSOR%INDEX(IATOM,JATOM,1,1)
      IF(I.NE.0)THEN
         slsao => TENSOR%SLSAO(I)
         IF(associated(slsao%selms))THEN
            Ibatch = slsao%AOBATCH(1)
            Jbatch = slsao%AOBATCH(2)
            nBatI = TENSOR%nAOBATCH(IATOM,1)
            nBatJ = TENSOR%nAOBATCH(JATOM,2)
            n1=slsao%nLocal(1)
            n2=slsao%nLocal(2)
            maxbat=slsao%maxbat
            do jbat = 1,nBatJ
               s2=slsao%startLocalOrb(jbat+maxbat)-1
               nB=slsao%nOrb(jbat+maxbat)
               do ibat = 1,nBatI
                  s1=slsao%startLocalOrb(ibat)-1
                  nA=slsao%nOrb(ibat)
                  call determine_maxprimgab(maxprimgab,slsao%selms,n1,n2,s1,s2,nA,nB)
                  TENSOR%maxprimgab(Ibatch + ibat,Jbatch + jbat) = maxprimgab
               enddo
            enddo
         ENDIF
      ENDIF
!      Ibatch = Ibatch + nbatI
     ENDDO
!     Jbatch = Jbatch + nbatJ
    ENDDO
  END SUBROUTINE Set_lst_maxprimgab

  subroutine determine_maxprimgab(maxprimgab,elms,n1,n2,s1,s2,nA,nB)
    implicit none
    integer,intent(in) :: n1,n2,s1,s2,nA,nB
    integer(kind=short),intent(in)    :: elms(n1,n2)
    integer(kind=short),intent(inout) :: maxprimgab
    !
    integer :: I,J
    maxprimgab = shortzero
    DO J = 1,nB
       DO I = 1,nA
          IF(elms(s1+I,s2+J).GT.maxprimgab)THEN
             maxprimGab = elms(s1+I,s2+J)
          ENDIF
       ENDDO
    ENDDO
  END subroutine determine_maxprimgab

  !> \brief build nuclear lstensor
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param fullmat the full nuclear screening matrix 
  !> \param TENSOR the lstensor
  !> \param natoms the number of atoms
  SUBROUTINE build_Nuclearlstensor(fullmat,TENSOR,natoms)
    implicit none
    TYPE(LSTENSOR)     :: TENSOR
    INTEGER            :: nAtoms
    real(realk)        :: fullmat(natoms,1)
    !
    Integer :: I,J,lstmem_index
    integer(kind=long) :: nmemsize,AllocInt,AllocRealk,AllocIntS,nsize
    integer(kind=short) :: maxGab
    call lstensor_nullify(TENSOR)
    TENSOR%pChargetensor = .FALSE.
    TENSOR%Econtrib = .FALSE.
    TENSOR%Screentensor = .TRUE.
    TENSOR%MagGradienttensor = .FALSE.
    TENSOR%Gradienttensor = .FALSE.

    TENSOR%natom(1) = natoms
    TENSOR%natom(2) = 1
    TENSOR%natom(3) = 1
    TENSOR%natom(4) = 1
    TENSOR%nbast(1) = 1
    TENSOR%nbast(2) = 1
    TENSOR%nbast(3) = 1
    TENSOR%nbast(4) = 1
    TENSOR%nbatches(1) = natoms
    TENSOR%nbatches(2) = 1
    TENSOR%nbatches(3) = 1
    TENSOR%nbatches(4) = 1
    nullify(TENSOR%maxgab)
    call mem_alloc(TENSOR%maxgab,natoms,1)
    nsize = size(TENSOR%maxgab,KIND=long)*mem_shortintsize
    call mem_allocated_mem_lstensor(nsize)
    maxgab = shortzero
    DO I = 1,natoms
       IF(abs(fullmat(I,1)).GT.shortintCRIT)THEN
          TENSOR%maxgab(I,1) = CEILING(LOG10(ABS(fullmat(I,1))))  !ABS right???
          IF(TENSOR%maxgab(I,1).GT.maxgab)THEN
             maxgab = TENSOR%maxgab(I,1)
          ENDIF
       ELSE
          TENSOR%maxgab(I,1) = shortzero
       ENDIF
    ENDDO
    TENSOR%maxgabelm = maxgab
    TENSOR%maxprimgabelm = maxgab
    nullify(TENSOR%maxprimgab)
    call mem_alloc(TENSOR%maxprimgab,natoms,1)
    nsize = size(TENSOR%maxprimgab,KIND=long)*mem_shortintsize
    call mem_allocated_mem_lstensor(nsize)
    DO I = 1,natoms
       TENSOR%maxprimgab(I,1) = TENSOR%maxgab(I,1)
    ENDDO
    nullify(TENSOR%nAOBATCH)
    call mem_alloc(TENSOR%nAOBATCH,natoms,2)
    nsize = size(TENSOR%nAOBATCH,KIND=long)*mem_intsize
    call mem_allocated_mem_lstensor(nsize)
    DO J=1,2
       DO I = 1,natoms
          TENSOR%nAOBATCH(I,J) = 1
       ENDDO
    ENDDO
    nullify(TENSOR%INDEX)
    call mem_alloc(TENSOR%INDEX,natoms,1,1,1)
    nsize = size(TENSOR%INDEX,KIND=long)*mem_intsize
    call mem_allocated_mem_lstensor(nsize)
    nullify(TENSOR%SLSAO)
    call mem_alloc(TENSOR%SLSAO,natoms)
    AllocInt = natoms*2
    AllocRealk = 0
    AllocInts = natoms
    call init_lstensorMem(AllocInt,AllocRealk,AllocInts,lstmem_index)
    call zero_lstensorMem
    TENSOR%lstmem_index = lstmem_index
    DO I = 1,natoms
       TENSOR%INDEX(I,1,1,1) = I
       call slsaotensor_nullify(TENSOR%SLSAO(I))
       call mem_lstpointer_alloc(TENSOR%SLSAO(I)%selms,1)
       call mem_lstpointer_alloc(TENSOR%SLSAO(I)%startLocalOrb,2)
       TENSOR%SLSAO(I)%startLocalOrb(1) = 1
       TENSOR%SLSAO(I)%startLocalOrb(2) = 1
       TENSOR%SLSAO(I)%nLocal(1) = 1
       TENSOR%SLSAO(I)%nLocal(2) = 1
       TENSOR%SLSAO(I)%selms(1) = TENSOR%maxgab(I,1)
       TENSOR%SLSAO(I)%nelms = 1
       TENSOR%SLSAO(I)%atom = 1
       TENSOR%SLSAO(I)%maxbat = 1
    ENDDO
    TENSOR%nSLSAO = natoms 
    TENSOR%nLSAO = 0
    nullify(TENSOR%LSAO)

  END SUBROUTINE build_Nuclearlstensor

!> \brief wrapper routine for building a single type matrix from lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output
!> \param TENSOR the lstensor
!> \param MAT the type matrix
SUBROUTINE Build_singlemat_from_lst(lupri,TENSOR,MAT)
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(MATRIX)       :: MAT
integer :: lupri
!
real(realk),pointer  :: fullMAT(:,:,:,:,:)
integer              :: I,nbast(4)
select case(matrix_type)
case(mtype_dense)
   call Build_single_dense_mat_from_lst(TENSOR,MAT)
case(mtype_unres_dense)
   call Build_single_unres_mat_from_lst(TENSOR,MAT)
case(mtype_scalapack)
   call mem_alloc(fullMAT,TENSOR%nbast(1),TENSOR%nbast(2),TENSOR%nbast(3),TENSOR%nbast(4),1)
   call Build_full_5dim_from_lstensor(TENSOR,fullMAT,TENSOR%nbast(1),TENSOR%nbast(2),TENSOR%nbast(3),TENSOR%nbast(4),1)
   IF(TENSOR%nbast(2) .EQ. 1) call mat_set_from_full(fullmat,1E0_realk, MAT)
   IF(TENSOR%nbast(3) .EQ. 1) call mat_set_from_full(fullmat,1E0_realk, MAT)
   call mem_dealloc(fullMAT)
case(mtype_csr)
   call Build_single_csr_mat_from_lst(TENSOR,MAT)
case(mtype_pdmm)
   call mem_alloc(fullMAT,TENSOR%nbast(1),TENSOR%nbast(2),TENSOR%nbast(3),TENSOR%nbast(4),1)
   call Build_full_5dim_from_lstensor(TENSOR,fullMAT,TENSOR%nbast(1),TENSOR%nbast(2),TENSOR%nbast(3),TENSOR%nbast(4),1)
   IF(TENSOR%nbast(2) .EQ. 1) call mat_set_from_full(fullmat,1E0_realk, MAT)
   IF(TENSOR%nbast(3) .EQ. 1) call mat_set_from_full(fullmat,1E0_realk, MAT)
   call mem_dealloc(fullMAT)
case default
   stop "BUILD_SINGLEMAT_FROM_LST not implemented for this type of matrix"
end select
END SUBROUTINE BUILD_SINGLEMAT_FROM_LST

!> \brief build a single type dense matrix from lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT2 the type matrix
SUBROUTINE Build_single_dense_mat_from_lst(TENSOR,MAT,nrow,ncol)
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(MATRIX),target       :: MAT
integer, optional :: nrow,ncol
!
TYPE(LSTAUXITEM) ::LSTAUX
INTEGER    :: I,nrow2,ncol2,nbast(4),DIMI,DIMJ,SUM1,SUM2,J
INTEGER,pointer :: ROW(:,:),COL(:,:),NNZ(:)
REAL(REALK),pointer :: MAT2dim(:,:),MAT3dim(:,:,:),MAT4dim(:,:,:,:),MAT5dim(:,:,:,:,:)
LOGICAL :: MATELMZERO
nrow2 = MAT%nrow
ncol2 = MAT%ncol
IF(present(nrow))nrow2 = nrow
IF(present(nrow))ncol2 = ncol
nbast=TENSOR%nbast
!print*,'nbast',nbast
!print*,'MAT%nrow,MAT%ncol',MAT%nrow,MAT%ncol
IF(TENSOR%ndim5 .NE. 1)CALL lsQUIT('ERROR: Build_single_dense_mat_from_lst nmat > 1',-1)
DIMI = 0
DO I=1,4
   IF(nrow2 .EQ.nbast(I))THEN
      DIMI = I
      EXIT
   ENDIF
ENDDO
DIMJ=0
DO J=1,4
   IF(J.NE.DIMI)THEN
      IF(ncol2 .EQ.nbast(J))THEN
         DIMJ = J
         EXIT
      ENDIF
   ENDIF
ENDDO
!print*,'DIMI,DIMJ',DIMI,DIMJ
SUM1=nbast(DIMI)+nbast(DIMJ)
SUM2=nrow2+ncol2
IF(SUM1.NE.SUM2)CALL lsQUIT('ERROR: Build_single_dense_mat_from_lst dim do not match',-1)
MAT%elms = 0E0_realk
call nullify_LSTAUXITEM(LSTAUX)
LSTAUX%sMAT => MAT
MATELMZERO=.FALSE.
call BuildFromlstensor_aux(TENSOR,DIMI,DIMJ,MAT%nrow,MATELMZERO,LSTAUX,sMat_FuncLstSub)

END SUBROUTINE Build_single_dense_mat_from_lst

!> \brief build a single csr (compress sparse row) type matrix from lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT the type matrix
subroutine Build_single_csr_mat_from_lst(TENSOR,MAT)
!include 'mkl_spblas.fi'
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(MATRIX)       :: MAT
!
type(LSTAUXITEM) :: LSTAUX
INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,Iangmom,Jangmom,Kangmom,Langmom,nOrbA
INTEGER    :: sA,nOrbB,sB,nOrbC,sC,J,DIMI,DIMJ,SUM1,SUM2
INTEGER    :: nOrbD,sD,IORB,JORB,KORB,LORB,IMAT,IELMs
INTEGER    :: IATOM,JATOM,KATOM,LATOM,nbat1,nbat2,nbat3,nbat4
INTEGER    :: nbast(4),MATINDEX(4),job(8),n,info,l1,olda1,l2
INTEGER    :: offset1,offset2,offset3,offset4,maxBat,maxAng
INTEGER,target :: nnz
INTEGER,pointer    :: ROW(:),COL(:),IELM(:)
real(realk),pointer :: VAL(:),MAT1dim(:),MAT3dim(:,:,:),MAT4dim(:,:,:,:),Mat5dim(:,:,:,:,:)
TYPE(LSAOTENSOR),pointer    :: lsao
REAL(REALK),pointer     :: elms(:)
LOGICAL :: MATELMZERO
if(mat%ncol.NE.mat%nrow)then
   call lsquit('Build_single_csr_mat_from_lst requires a quadratic matrix',-1)
endif
nbast=TENSOR%nbast
IF(TENSOR%ndim5 .NE. 1)CALL lsQUIT('ERROR: Build_single_csr_mat_from_lst nmat > 1',-1)
DO I=1,4
   IF(MAT%nrow .EQ.nbast(I))THEN
      DIMI = I
      EXIT
   ENDIF
ENDDO
DO J=1,4
   IF(J.NE.DIMI)THEN
      IF(MAT%ncol .EQ.nbast(J))THEN
         DIMJ = J
         EXIT
      ENDIF
   ENDIF
ENDDO
SUM1=nbast(DIMI)+nbast(DIMJ)
SUM2=MAT%nrow+MAT%ncol
IF(SUM1.NE.SUM2)CALL lsQUIT('ERROR: Build_single_csr_mat_from_lst dim do not match',-1)

! PROPER WAY OF DOING THIS :
!
! ATTACH TO BASIS AND THEN TO AOBATCH A LOOPORDER FOR EACH TYPE
! ALSO ATTACH A TYPE FOR ALL ATOMS
! SO THAT I CAN DO
! do iatom = 1,tensor%natom1
!    I = tensor%index(Iatom,1,1,1)
!    ID = atomtype(iatom)
!    ib=1,length(order(ID)%batch(:))
!    ibat = order(ID)%batch(ib)
!    iang = order(ID)%ang(ib)
!    calc sA 
!    sA = 0
!    DO SAA=2,iang
!      sA = sA+nOrbComp(SAA-1)*lsao%BATCH(Ibat,1,1,1)%nOrbA(SAA-1)
!    ENDDO
!    nOrbCompA = lsao%BATCH(Ibat,1,1,1)%nOrbCompA(iang)
!    startA = lsao%BATCH(Ibat,1,1,1)%startOrbA(Iang)
!    do iContA = 1,lsao%BATCH(Ibat,1,1,1)%nOrbA(iang)
!     do iOrbCompA = 1,
!      ielm = sA+iOrbCompA+(icontA)*nOrbCompA
!      orbitalindex = startA + iOrbCompA+(icontA)*nOrbCompA
!     enddo
!    enddo
! enddo
! this should mean that orbitalindex is now consequtive

! FOR NOW WE JUST DO THIS
NNZ = 0
call nullify_LSTAUXITEM(LSTAUX)
LSTAUX%NNZ1 => NNZ
MATELMZERO=.FALSE.
call BuildFromlstensor_aux(TENSOR,DIMI,DIMJ,0,MATELMZERO,LSTAUX,NNZ1_FuncLstSub)

!call mat_csr_allocate(MAT,NNZ,MAT%nrow) Stinne, Rasmus 2/8-2010
call mat_csr_allocate(MAT,NNZ)

call mem_alloc(VAL,NNZ)
call mem_alloc(ROW,NNZ)
call mem_alloc(COL,NNZ)
NNZ = 0
call nullify_LSTAUXITEM(LSTAUX)
LSTAUX%NNZ1 => NNZ
LSTAUX%VAL1 => VAL
LSTAUX%COL1 => COL
LSTAUX%ROW1 => ROW
MATELMZERO = .FALSE.
call BuildFromlstensor_aux(TENSOR,DIMI,DIMJ,0,MATELMZERO,LSTAUX,CSR1_FuncLstSub)
!WRITE(6,*)'THE COORDINATE FORM   nnz=',nnz
!WRITE(6,'(2X,A4,5E13.3,/(6X,5E13.3))')'VAL:',(VAL(j),j=1,nnz)
!WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'COL:',(COL(j),j=1,nnz)
!WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'ROW:',(row(j),j=1,nnz)

call QSORT_csr(ROW,COL,VAL)

L1=1
olda1 = ROW(L1)
L2=1
do I=2,NNZ
   IF(olda1.EQ.ROW(I))THEN
      L2=L2+1
   ELSE
      call QSORT_csr(COL(L1:L2),ROW(L1:L2),VAL(L1:L2))
      L1=L2+1
      L2=L2+1
      olda1 = ROW(L1)
   ENDIF
   IF(I.EQ.NNZ)THEN
      call QSORT_csr(COL(L1:NNZ),ROW(L1:NNZ),VAL(L1:NNZ))
!      WRITE(6,*)'AFTER FINAL SORT'
!      WRITE(6,'(2X,A4,5E13.3,/(6X,5E13.3))')'VAL:',(VAL(j),j=1,nnz)
!      WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'COL:',(COL(j),j=1,nnz)
!      WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'ROW:',(ROW(j),j=1,nnz)
   ENDIF
enddo

job(1)=1
job(2)=1
job(3)=1
job(4)=1
job(5)=nnz
job(6)=0
job(7)=1
job(8)=1
n = mat%nrow
!print*,'MKL: n',n
!WRITE(6,*)'THE COORDINATE FORM   nnz=',nnz,'n',n
!WRITE(6,'(2X,A4,5E13.3,/(6X,5E13.3))')'VAL:',(VAL(j),j=1,nnz)
!WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'COL:',(COL(j),j=1,nnz)
!WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'ROW:',(row(j),j=1,nnz)
#ifdef VAR_CSR
call mkl_dcsrcoo(job,n,mat%val,mat%col,mat%row,nnz,VAL,ROW,COL,info)
#else
call lsquit('Build_single_csr_mat_from_lst requires VAR_CSR',-1)
#endif
!print*,'THE CSR from MKL  '
!call mat_print(mat,1,mat%nrow,1,mat%ncol,6)
call mem_dealloc(VAL)
call mem_dealloc(ROW)
call mem_dealloc(COL)

END subroutine Build_single_csr_mat_from_lst

!> \brief build an array of csr (compress sparse row) type matrix from lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT the array of type matrix
subroutine Build_array_csr_mat_from_lst(TENSOR,MAT)
use matrix_operations_csr
!include 'mkl_spblas.fi'
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(MATRIX)       :: MAT(:)
!
type(LSTAUXITEM) :: LSTAUX
INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,DIMI,DIMJ,SUM1,SUM2,J,nmat,nbast(4),IMAT,maxnnz
INTEGER    :: L1,job(8),n,info,olda1,l2
INTEGER,pointer    :: ROW(:,:),COL(:,:),IELM(:,:),nnz(:)
real(realk),pointer :: VAL(:,:)
REAL(REALK),pointer :: MAT1dim(:),MAT3dim(:,:,:),MAT4dim(:,:,:,:),MAT5dim(:,:,:,:,:)
LOGICAL :: MATELMZERO

if(mat(1)%ncol .NE. mat(1)%nrow)then
   call lsquit('Build_single_csr_mat_from_lst requires a quadratic matrix',-1)
endif
nbast=TENSOR%nbast
nmat = TENSOR%ndim5
IF(nmat.NE.size(MAT))call lsquit('dim mismatch in Build_array_csr_mat_from_lst',-1)
call mem_alloc(nnz,nmat)
DO I=1,4
   IF(MAT(1)%nrow .EQ.nbast(I))THEN
      DIMI = I
      EXIT
   ENDIF
ENDDO
DO J=1,4
   IF(J.NE.DIMI)THEN
      IF(MAT(1)%ncol .EQ.nbast(J))THEN
         DIMJ = J
         EXIT
      ENDIF
   ENDIF
ENDDO
SUM1=nbast(DIMI)+nbast(DIMJ)
SUM2=MAT(1)%nrow+MAT(1)%ncol
IF(SUM1.NE.SUM2)CALL lsQUIT('ERROR: Build_single_csr_mat_from_lst dim do not match',-1)

do IMAT=1,nmat
   NNZ(IMAT) = 0
enddo
call nullify_LSTAUXITEM(LSTAUX)
LSTAUX%NNZ2 => NNZ
MATELMZERO = .FALSE.
call BuildFromlstensor_aux(TENSOR,DIMI,DIMJ,0,MATELMZERO,LSTAUX,NNZ2_FuncLstSub)

do IMAT=1,nmat
   call mat_csr_allocate(MAT(IMAT),NNZ(IMAT))
enddo

maxnnz = 0
do IMAT=1,nmat
   maxnnz = max(maxnnz,nnz(IMAT))
   NNZ(IMAT) = 0
enddo

call mem_alloc(VAL,maxNNZ,NMAT)
call mem_alloc(ROW,maxNNZ,NMAT)
call mem_alloc(COL,maxNNZ,NMAT)
call nullify_LSTAUXITEM(LSTAUX)
LSTAUX%NNZ2 => NNZ
LSTAUX%VAL2 => VAL
LSTAUX%COL2 => COL
LSTAUX%ROW2 => ROW
MATELMZERO = .FALSE.
call BuildFromlstensor_aux(TENSOR,DIMI,DIMJ,0,MATELMZERO,LSTAUX,CSR2_FuncLstSub)

DO IMAT = 1,TENSOR%ndim5
   call QSORT_csr(ROW(1:NNZ(IMAT),IMAT),COL(1:NNZ(IMAT),IMAT),VAL(1:NNZ(IMAT),IMAT))

   L1=1
   olda1 = ROW(L1,IMAT)
   L2=1
   do I=2,NNZ(IMAT)
      IF(olda1.EQ.ROW(I,IMAT))THEN
         L2=L2+1
      ELSE
         call QSORT_csr(COL(L1:L2,IMAT),ROW(L1:L2,IMAT),VAL(L1:L2,IMAT))
         L1=L2+1
         L2=L2+1
         olda1 = ROW(L1,IMAT)
      ENDIF
      IF(I.EQ.NNZ(IMAT))THEN
         call QSORT_csr(COL(L1:NNZ(IMAT),IMAT),ROW(L1:NNZ(IMAT),IMAT),VAL(L1:NNZ(IMAT),IMAT))
         !      WRITE(6,*)'AFTER FINAL SORT'
         !      WRITE(6,'(2X,A4,5E13.3,/(6X,5E13.3))')'VAL:',(VAL(j),j=1,nnz)
         !      WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'COL:',(COL(j),j=1,nnz)
         !      WRITE(6,'(2X,A4,15I4,/(6X,15I4))')'ROW:',(ROW(j),j=1,nnz)
      ENDIF
   enddo

   job(1)=1
   job(2)=1
   job(3)=1
   job(4)=1
   job(5)=nnz(IMAT)
   job(6)=0
   job(7)=1
   job(8)=1
   n = mat(IMAT)%nrow
   !call mat_print(mat,1,mat%nrow,1,mat%nrow,6)
#ifdef VAR_CSR
   call mkl_dcsrcoo(job,n,mat(IMAT)%val,mat(IMAT)%col,mat(IMAT)%row,nnz(IMAT),&
        &VAL(1:NNZ(IMAT),IMAT),ROW(1:NNZ(IMAT),IMAT),COL(1:NNZ(IMAT),IMAT),info)
#else
   call lsquit('Build_array_csr_mat_from_lst requires VAR_CSR',-1)
#endif

enddo
!print*,'THE CSR from MKL  '
!call mat_print(mat(IMAT),1,mat(IMAT)%nrow,1,mat(IMAT)%nrow,6)
call mem_dealloc(nnz)
call mem_dealloc(VAL)
call mem_dealloc(ROW)
call mem_dealloc(COL)

END subroutine Build_array_csr_mat_from_lst

!> \brief recursive quick sort algorithm, smallest first
!> \author T. Kjaergaard
!> \date 2010
!> \param a1 the values to be sorted
!> \param a2 and array that should be sorted like a1
!> \param val and array that should be sorted like a1
RECURSIVE SUBROUTINE Qsort_csr(a1,a2,val)
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: a1(:),a2(:)
REAL(REALK) :: VAL(:)
INTEGER :: split
INTEGER :: I

  IF(size(a1) > 7) THEN 
     CALL Partition_csr(a1,a2,VAL,split)
     CALL Qsort_csr(a1(:split-1),a2(:split-1),VAL(:split-1))
     CALL Qsort_csr(a1(split:),a2(split:),VAL(split:))
  ELSEIF(size(a1) > 1)THEN
     CALL INSERTION_csr(a1,a2,VAL)
  END IF
 
END SUBROUTINE Qsort_csr

!> \brief insertion sort algorithm, smallest first
!> \author T. Kjaergaard
!> \date 2010
!> \param a the values to be sorted
!> \param a2 and array that should be sorted like a1
!> \param INDEXES and array that should be sorted like a1
SUBROUTINE INSERTION_csr(a,a2,INDEXES)
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: a(:),a2(:)
INTEGER                    :: I,J,temp,tempI2
real(realk) :: INDEXES(:),tempI

Do I = 2, Size(a)
   temp = a(I)
   tempI = INDEXES(I)
   tempI2 = a2(I)
   DO J = I-1, 1, -1
      IF(temp .LT. a(J)) Then
         a(J+1) = a(J)
         INDEXES(J+1) = INDEXES(J)
         a2(J+1) = a2(J)
      ELSE
         Exit
      ENDIF
   ENDDO
   a(J+1) = temp
   INDEXES(J+1) = tempI
   a2(J+1) = tempI2
End Do

END SUBROUTINE INSERTION_csr

!> \brief the worker routine used by the qucik sort (qsort_csr) routine, smallest first
!> \author T. Kjaergaard
!> \date 2010
!> \param a the values to be sorted
!> \param a2 and array that should be sorted like a1
!> \param INDEXES and array that should be sorted like a1
!> \param marker counter index
SUBROUTINE Partition_csr(a, a2,INDEXES, marker)
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: a(:),a2(:)
INTEGER, INTENT(OUT)       :: marker
INTEGER                    :: left, right, I, temp,pivot,tempI2
REAL(REALK)                :: INDEXES(:),tempI
 
pivot = (a(1) + a(size(a))) / 2 ! Average of first and last elements 
                                !to prevent quadratic behavior with
                                ! sorted or reverse sorted data
IF(a(1) .EQ. pivot)THEN
   DO I =size(a)-1,2,-1
      IF(a(I) .NE. pivot)THEN
         pivot = (a(1) + a(I)) / 2
         EXIT
      ENDIF
   ENDDO
ENDIF

left = 0                         ! behavior with sorted or reverse sorted data
right = size(a) + 1

DO WHILE (left < right)
   right = right - 1
   DO WHILE (a(right) > pivot)
      right = right-1
   END DO
   left = left + 1
   DO WHILE (a(left) < pivot)
      left = left + 1
   END DO
   IF (left < right) THEN 
      temp = a(left)
      a(left) = a(right)
      a(right) = temp
      tempI = INDEXES(left)
      tempI2 = a2(left)
      INDEXES(left) = INDEXES(right)
      a2(left) = a2(right)
      INDEXES(right) = tempI
      a2(right) = tempI2
   END IF
END DO

IF (left == right) THEN
   marker = left + 1
ELSE
   marker = left
END IF

END SUBROUTINE Partition_csr

!> \brief build a single unres_dense type matrix from lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT the array of type matrix
SUBROUTINE Build_single_unres_mat_from_lst(TENSOR,MAT)
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(MATRIX),target :: MAT
!
INTEGER    :: I,J,DIMI,DIMJ,nbast(4),SUM1,SUM2
TYPE(LSAOTENSOR),pointer    :: lsao
REAL(REALK),pointer     :: elms(:)
TYPE(LSTAUXITEM) :: LSTAUX
logical    :: option1
LOGICAL :: MATELMZERO

nbast=TENSOR%nbast
option1 =.true.
IF(TENSOR%ndim5 .EQ. 2)option1 =.false. !different alpha beta parts
IF(TENSOR%ndim5 .GT. 2)CALL lsQUIT('ERROR: Build_single_dense_mat_from_lst nmat > 1',-1)
DO I=1,4
   IF(MAT%nrow .EQ.nbast(I))THEN
      DIMI = I
      EXIT
   ENDIF
ENDDO
DO J=1,4
   IF(J.NE.DIMI)THEN
      IF(MAT%ncol .EQ.nbast(J))THEN
         DIMJ = J
         EXIT
      ENDIF
   ENDIF
ENDDO
SUM1=nbast(DIMI)+nbast(DIMJ)
SUM2=MAT%nrow+MAT%ncol
IF(SUM1.NE.SUM2)CALL lsQUIT('ERROR: Build_single_unres_mat_from_lst dim do not match',-1)
call mat_zero(MAT)

call nullify_LSTAUXITEM(LSTAUX)
LSTAUX%sMat => MAT
if(option1)then
   MATELMZERO=.FALSE.
   call BuildFromlstensor_aux(TENSOR,DIMI,DIMJ,MAT%nrow,MATELMZERO,LSTAUX,opt1UnresMatS_FuncLstSub)
else
   MATELMZERO=.FALSE.
   call BuildFromlstensor_aux(TENSOR,DIMI,DIMJ,MAT%nrow,MATELMZERO,LSTAUX,opt2UnresMatS_FuncLstSub)
endif
END SUBROUTINE Build_single_unres_mat_from_lst

!> \brief wrapper routine for building an array of type matrix from lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number of the output file
!> \param TENSOR the lstensor
!> \param MAT the array of type matrix
SUBROUTINE Build_matarray_from_lst(lupri,TENSOR,MAT)
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(MATRIX)       :: MAT(:)
integer            :: lupri
!
real(realk),pointer  :: fullMAT(:,:,:,:,:)
integer              :: I
!print*,'build_matrixarray matrix'

select case(matrix_type)
case(mtype_dense)
   call Build_dense_matarray_from_lst(TENSOR,MAT)
case(mtype_scalapack)
   call mem_alloc(fullMAT,TENSOR%nbast(1),TENSOR%nbast(2),TENSOR%nbast(3),TENSOR%nbast(4),TENSOR%ndim5)
   call Build_full_5dim_from_lstensor(TENSOR,fullMAT,TENSOR%nbast(1),TENSOR%nbast(2),TENSOR%nbast(3),TENSOR%nbast(4),TENSOR%ndim5)
   IF(TENSOR%nbast(2) .EQ. 1) THEN
      do I=1,TENSOR%ndim5
         call mat_set_from_full(fullmat(:,1,:,1,I),1E0_realk, MAT(I))
      enddo
!      call ls_output(fullmat,1,TENSOR%nbast(1),1,TENSOR%nbast(3),TENSOR%nbast(1),TENSOR%nbast(3),1,lupri)
   ENDIF
   IF(TENSOR%nbast(3) .EQ. 1) THEN
      do I=1,TENSOR%ndim5
         call mat_set_from_full(fullmat(:,:,1,1,I),1E0_realk, MAT(I))
      enddo
!      call ls_output(fullmat,1,TENSOR%nbast(1),1,TENSOR%nbast(2),TENSOR%nbast(1),TENSOR%nbast(2),1,lupri)
   ENDIF
   call mem_dealloc(fullMAT)
case(mtype_pdmm)
   call mem_alloc(fullMAT,TENSOR%nbast(1),TENSOR%nbast(2),TENSOR%nbast(3),TENSOR%nbast(4),TENSOR%ndim5)
   call Build_full_5dim_from_lstensor(TENSOR,fullMAT,TENSOR%nbast(1),TENSOR%nbast(2),TENSOR%nbast(3),TENSOR%nbast(4),TENSOR%ndim5)
   IF(TENSOR%nbast(2) .EQ. 1) THEN
      do I=1,TENSOR%ndim5
         call mat_set_from_full(fullmat(:,1,:,1,I),1E0_realk, MAT(I))
      enddo
!      call ls_output(fullmat,1,TENSOR%nbast(1),1,TENSOR%nbast(3),TENSOR%nbast(1),TENSOR%nbast(3),1,lupri)
   ENDIF
   IF(TENSOR%nbast(3) .EQ. 1) THEN
      do I=1,TENSOR%ndim5
         call mat_set_from_full(fullmat(:,:,1,1,I),1E0_realk, MAT(I))
      enddo
!      call ls_output(fullmat,1,TENSOR%nbast(1),1,TENSOR%nbast(2),TENSOR%nbast(1),TENSOR%nbast(2),1,lupri)
   ENDIF
   call mem_dealloc(fullMAT)
case(mtype_unres_dense)
   call Build_unres_matarray_from_lst(TENSOR,MAT)
case(mtype_csr)
   call Build_array_csr_mat_from_lst(TENSOR,MAT)
case default
   stop "BUILD_MATARRAY_FROM_LST not implemented for this type of matrix"
end select

END SUBROUTINE BUILD_MATARRAY_FROM_LST

!> \brief build an array of dense type matrix from lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT273 the array of type matrix
SUBROUTINE Build_dense_matarray_from_lst(TENSOR,MAT)
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(MATRIX),target :: MAT(:)
!
INTEGER    :: I,J,DIMI,DIMJ,SUM1,SUM2,IMAT,nbast(4)
TYPE(LSTAUXITEM) :: LSTAUX
LOGICAL :: MATELMZERO
nbast=TENSOR%nbast
DO I=1,4
   IF(MAT(1)%nrow .EQ.nbast(I))THEN
      DIMI = I
      EXIT
   ENDIF
ENDDO
DO J=1,4
   IF(J.NE.DIMI)THEN
      IF(MAT(1)%ncol .EQ.nbast(J))THEN
         DIMJ = J
         EXIT
      ENDIF
   ENDIF
ENDDO
SUM1=nbast(DIMI)+nbast(DIMJ)
SUM2=MAT(1)%nrow+MAT(1)%ncol
IF(SUM1.NE.SUM2)CALL lsQUIT('ERROR: Build_dense_matarray_from_lst dim do not match',-1)
call nullify_LSTAUXITEM(LSTAUX)
call mem_alloc(LSTAUX%aMAT,TENSOR%ndim5)
DO I = 1,TENSOR%ndim5
   MAT(I)%elms = 0E0_realk
   LSTAUX%aMAT(I)%p => MAT(I)
ENDDO
MATELMZERO=.FALSE.
call BuildFromlstensor_aux(TENSOR,DIMI,DIMJ,MAT(1)%nrow,MATELMZERO,LSTAUX,aMat_FuncLstSub)
call mem_dealloc(LSTAUX%aMAT)
END SUBROUTINE Build_dense_matarray_from_lst

!> \brief build an array of unres_dense type matrix from lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT273 the array of type matrix
SUBROUTINE Build_unres_matarray_from_lst(TENSOR,MAT)
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(MATRIX),target       :: MAT(:)
!
INTEGER    :: I,J,DIMI,DIMJ,SUM1,SUM2,IMAT,nbast(4),nmat
TYPE(LSTAUXITEM) :: LSTAUX
logical    :: option1
LOGICAL :: MATELMZERO

nbast=TENSOR%nbast
option1 = .true. !same alpha beta part
IF(TENSOR%ndim5 .EQ. (2*size(MAT))) option1 = .false. !different alpha beta part

IF (.NOT.option1) THEN
  nmat = TENSOR%ndim5/2
ELSE
  nmat = TENSOR%ndim5
ENDIF
IF(nmat.NE.size(MAT))call lsquit('nmat mismatch Build_unres_matarray_from_lst',-1)
DIMI=0
DIMJ=0
DO I=1,4
   IF(MAT(1)%nrow .EQ.nbast(I))THEN
      DIMI = I
      EXIT
   ENDIF
ENDDO
DO J=1,4
   IF(J.NE.DIMI)THEN
      IF(MAT(1)%ncol .EQ.nbast(J))THEN
         DIMJ = J
         EXIT
      ENDIF
   ENDIF
ENDDO
SUM1=nbast(DIMI)+nbast(DIMJ)
SUM2=MAT(1)%nrow+MAT(1)%ncol
IF(SUM1.NE.SUM2)CALL lsQUIT('ERROR: Build_dense_matarray_from_lst dim do not match',-1)
call nullify_LSTAUXITEM(LSTAUX)
call mem_alloc(LSTAUX%aMAT,nmat)
DO I = 1,nmat
   call mat_zero(MAT(I))
   LSTAUX%aMat(I)%p => MAT(I)
ENDDO
if(option1)then
   MATELMZERO = .FALSE.
   call BuildFromlstensor_aux(TENSOR,DIMI,DIMJ,MAT(1)%nrow,MATELMZERO,LSTAUX,opt1UnresMatA_FuncLstSub)
else
   MATELMZERO = .TRUE.
   call BuildFromlstensor_aux(TENSOR,DIMI,DIMJ,MAT(1)%nrow,MATELMZERO,LSTAUX,opt2UnresMatA_FuncLstSub)
endif
call mem_dealloc(LSTAUX%aMAT)

END SUBROUTINE Build_unres_matarray_from_lst

!> \brief wrapper routine for building an array of type matrix from lstensor
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number of the output file
!> \param TENSOR the lstensor
!> \param MAT the array of type matrix
SUBROUTINE Build_lst_from_matarray(TENSOR,MAT,AO1,AO2,nbast1,nbast2,nmat,useAO1,useAO2,ODscreen,lupri)
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(MATRIXp)       :: MAT(:)
TYPE(AOITEM),target :: AO1,AO2
INTEGER            :: nbast1,nbast2,nmat,lupri
logical            :: useAO1,useAO2,ODscreen
select case(matrix_type)
case(mtype_dense)
   call Build_lst_from_dense_matarray(TENSOR,MAT,AO1,AO2,nbast1,nbast2,nmat,useAO1,useAO2,ODscreen,lupri)
case(mtype_scalapack)
   call Build_lst_from_matarray_fallback(TENSOR,MAT,AO1,AO2,nbast1,nbast2,nmat,useAO1,useAO2,lupri)
case(mtype_pdmm)
   call Build_lst_from_matarray_fallback(TENSOR,MAT,AO1,AO2,nbast1,nbast2,nmat,useAO1,useAO2,lupri)
case(mtype_unres_dense)
   call lsquit('BUILD_LST_FROM_MATARRAY not implemented for unres',-1)
case default
   call Build_lst_from_matarray_fallback(TENSOR,MAT,AO1,AO2,nbast1,nbast2,nmat,useAO1,useAO2,lupri)
end select
CONTAINS 
  subroutine Build_lst_from_matarray_fallback(TENSOR,MAT,AO1,AO2,nbast1,nbast2,nmat,useAO1,useAO2,lupri)
    implicit none
    TYPE(LSTENSOR)     :: TENSOR
    TYPE(MATRIXP)       :: MAT(:)
    TYPE(AOITEM),target :: AO1,AO2
    INTEGER            :: nbast1,nbast2,nmat,lupri
    logical   :: useAO1,useAO2
    !
    integer :: I
    real(realk),pointer  :: fullMAT(:,:,:)
    call mem_alloc(fullMAT,nbast1,nbast2,nmat)
    do I=1,nmat
       call mat_to_full(MAT(I)%p, 1.0E0_realk,fullmat(:,:,I))
    enddo
    call Build_lstensor_from_full_3dim(TENSOR,fullMAT,AO1,AO2,nbast1,nbast2,nmat,useAO1,useAO2,lupri)
    call mem_dealloc(fullMAT)
  END subroutine BUILD_LST_FROM_MATARRAY_FALLBACK
END SUBROUTINE BUILD_LST_FROM_MATARRAY

!> \brief build a lstensor from full 3 dim array
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT2 the full 3dim array
!> \param AO1 the Atomic orbital(AO) item for center 1 
!> \param AO2 the AO item for center 2 
!> \param nbast1 the size of the 1. dimension
!> \param nbast2 the size of the 2. dimension
!> \param nmat the number of density matrices or size of the 5. dimension
!> \param useAO1 flag to describe if the AO1 item should be used or an empty AO 
!> \param useAO2 flag to describe if the AO1 item should be used or an empty AO 
!> \param lupri the logical unit number for the output
subroutine Build_lst_from_dense_matarray(TENSOR,MAT,AO1,AO2,nbast1,nbast2,nmat,useAO1,useAO2,ODscreen,lupri)
  implicit none
  TYPE(LSTENSOR)     :: TENSOR
  TYPE(MATRIXP)       :: MAT(:)
  TYPE(AOITEM),target :: AO1,AO2
  INTEGER            :: nbast1,nbast2,nmat,lupri
  logical   :: useAO1,useAO2,ODscreen
    !
  TYPE(AOITEM),target :: AOT
  TYPE(AOITEM),pointer :: AOT1,AOT2,AOT3,AOT4
  logical :: useAO3,useAO4
  LOGICAL :: MATELMZERO
  integer :: I
  TYPE(LSTAUXITEM) :: LSTAUX
  call SET_EMPTY_AO(AOT)
  AOT1 => AO1
  AOT2 => AO2
  AOT3 => AOT
  AOT4 => AOT
  IF(.NOT.useAO1) AOT1 => AOT
  IF(.NOT.useAO2) AOT2 => AOT
  useAO3=.TRUE.
  useAO4=.TRUE.
  call init_lstensor_5dim(TENSOR,AOT1,AOT2,AOT3,AOT4,nbast1,nbast2,&
       &1,1,nmat,useAO1,useAO2,useAO3,useAO4,ODscreen,.FALSE.,lupri)
  call FREE_EMPTY_AO(AOT)
  call nullify_LSTAUXITEM(LSTAUX)
  call mem_alloc(LSTAUX%aMAT,TENSOR%ndim5)
  DO I = 1,TENSOR%ndim5
     LSTAUX%aMAT(I)%p => MAT(I)%p
  ENDDO
  MATELMZERO = .FALSE.
  call BuildFromlstensor_aux(TENSOR,1,2,MAT(1)%p%nrow,MATELMZERO,LSTAUX,FromaMat_FuncLstSub)
  call mem_dealloc(LSTAUX%aMAT)
END subroutine BUILD_LST_FROM_DENSE_MATARRAY

!> \brief set the lstensor to zero
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!FIXME WHAT ABOUT ANTISYMMETRY
SUBROUTINE lstensor_force_symMat_to_triangularMat(TENSOR)
implicit none
TYPE(LSTENSOR)     :: TENSOR
!
INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,Iangmom,Jangmom,Kangmom,Langmom,nOrbA
INTEGER    :: sA,nOrbB,sB,nOrbC,sC,dimA,dimB,jatom,iatom,I2
INTEGER    :: nOrbD,sD,IORB,JORB,KORB,LORB,IELM,IMAT,i1,dim,ielms,ndim5
TYPE(LSAOTENSOR),pointer    :: lsao,lsao2

!only for 2 dimensional matrices
IF(TENSOR%nbast(3).NE.1.OR.TENSOR%nbast(4).NE.1)THEN
   CALL LSQUIT('Error in lstensor_force_symMat_to_triangularMat1',-1)
ENDIF
!only for square matrices
IF(TENSOR%nbast(1).NE.TENSOR%nbast(2))THEN
   CALL LSQUIT('Error in lstensor_force_symMat_to_triangularMat2',-1)
ENDIF
IF(TENSOR%natom(1).NE.TENSOR%natom(2))THEN
   CALL LSQUIT('Error in lstensor_force_symMat_to_triangularMat2',-1)
ENDIF
ndim5 = TENSOR%ndim5
!add lower triangular to upper triangular
DO Jatom = 1,TENSOR%natom(2)
 DO Iatom = 1,Jatom-1
  I = TENSOR%INDEX(IATOM,JATOM,1,1)
  I2 = TENSOR%INDEX(JATOM,IATOM,1,1)
  lsao => TENSOR%LSAO(I)
  lsao2 => TENSOR%LSAO(I2)
  IF(I.NE.0)THEN
     dimA = lsao%nLocal(1)
     dimB = lsao%nLocal(2)
     call force_symElms_to_triangularElms2(lsao%elms,lsao2%elms,dimA,dimB,ndim5)
  ENDIF
 ENDDO
 Iatom = Jatom
 I = TENSOR%INDEX(IATOM,JATOM,1,1)
 IF(I.NE.0)THEN
    lsao => TENSOR%LSAO(I)
    dimA = lsao%nLocal(1)
    dimB = lsao%nLocal(2)
    IF(dimA.NE.dimB)call lsquit('errornjbhjbhjbkjbhj',-1)
    call force_symElms_to_triangularElms(lsao%elms,dimA,dimB,ndim5)
 ENDIF
ENDDO
!set lower triangular to zero 
DO Jatom = 1,TENSOR%natom(2)
 Iatom = Jatom
 I = TENSOR%INDEX(IATOM,JATOM,1,1)
 IF(I.NE.0)THEN
    lsao => TENSOR%LSAO(I)
    dimA = lsao%nLocal(1)
    dimB = lsao%nLocal(2)
    IF(dimA.NE.dimB)call lsquit('errornjbhjbhjbkjbhj',-1)
    call set_lowerelms_triangular_zero(lsao%elms,dimA,dimB,ndim5)
 ENDIF
 DO Iatom = Jatom+1,TENSOR%natom(1)
    I = TENSOR%INDEX(IATOM,JATOM,1,1)
    IF(I.NE.0)THEN
       lsao => TENSOR%LSAO(I)
       dim = ndim5*lsao%nelms
       call ls_dzero(lsao%elms,dim)
    ENDIF
 ENDDO
ENDDO

END SUBROUTINE LSTENSOR_FORCE_SYMMAT_TO_TRIANGULARMAT

subroutine force_symElms_to_triangularElms2(elms,elms2,dimenA,dimenB,nmat)
  implicit none
  integer,intent(in) :: dimenA,dimenB,nmat
  real(realk) :: elms(dimenA,dimenB,nmat)
  real(realk),intent(in) :: elms2(dimenB,dimenA,nmat)
  !
  integer :: A,B,IMAT
  Do IMAT = 1,NMAT
     DO B=1,dimenB
        DO A=1,dimenA
           elms(A,B,IMAT) = elms(A,B,IMAT) + elms2(B,A,IMAT)
        ENDDO
     ENDDO
  ENDDO
end subroutine force_symElms_to_triangularElms2

subroutine force_symElms_to_triangularElms(elms,dimenA,dimenB,nmat)
  implicit none
  integer,intent(in) :: dimenA,dimenB,nmat
  real(realk) :: elms(dimenA,dimenB,nmat)
  !
  integer :: A,B,IMAT
  Do IMAT = 1,NMAT
     DO B=1,dimenB
        DO A=1,B-1
           elms(A,B,IMAT) = elms(A,B,IMAT) + elms(B,A,IMAT)
        ENDDO
     ENDDO
  ENDDO
end subroutine force_symElms_to_triangularElms

subroutine set_lowerelms_triangular_zero(elms,dimenA,dimenB,nmat)
  implicit none
  integer,intent(in) :: dimenA,dimenB,nmat
  real(realk) :: elms(dimenA,dimenB,nmat)
  !
  integer :: A,B,IMAT
  Do IMAT = 1,NMAT
     DO B=1,dimenB
        DO A=B+1,dimenA
           elms(A,B,IMAT) = 0.0E0_realk
        ENDDO
     ENDDO
  ENDDO
end subroutine set_lowerelms_triangular_zero

!> \brief set the lstensor to zero
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!FIXME WHAT ABOUT ANTISYMMETRY
SUBROUTINE lstensor_full_symMat_from_triangularMat(TENSOR)
implicit none
TYPE(LSTENSOR)     :: TENSOR
!
INTEGER    :: I,Ibat,Jbat,Kbat,Lbat,Iangmom,Jangmom,Kangmom,Langmom,nOrbA
INTEGER    :: sA,nOrbB,sB,nOrbC,sC,dimA,dimB,jatom,iatom,I2
INTEGER    :: nOrbD,sD,IORB,JORB,KORB,LORB,IELM,IMAT,i1,dim,ielms,nmat
TYPE(LSAOTENSOR),pointer    :: lsao,lsao2
!only for 2 dimensional matrices
IF(TENSOR%nbast(3).NE.1.OR.TENSOR%nbast(4).NE.1)THEN
   CALL LSQUIT('Error in lstensor_full_symMat_from_triangularMat1',-1)
ENDIF
IF(TENSOR%nbast(1).NE.TENSOR%nbast(2))THEN
   CALL LSQUIT('Error in lstensor_full_symMat_from_triangularMat2',-1)
ENDIF
IF(TENSOR%natom(1).NE.TENSOR%natom(2))THEN
   CALL LSQUIT('Error in lstensor_full_symMat_from_triangularMat3',-1)
ENDIF
nmat = TENSOR%ndim5
DO Jatom = 1,TENSOR%natom(2)
 DO Iatom = 1,Jatom-1
  I = TENSOR%INDEX(IATOM,JATOM,1,1)
  I2 = TENSOR%INDEX(JATOM,IATOM,1,1)
  IF(I.NE.0)THEN
     lsao => TENSOR%LSAO(I)
     IF(I2.EQ.0)call lsquit('erenkjbkjnlkj',-1)
     lsao2 => TENSOR%LSAO(I2)
     dimA = lsao%nLocal(1)
     dimB = lsao%nLocal(2)
     call FULL_SYMelms_FROM_TRIANGULARelms(lsao2%elms,lsao%elms,dimA,dimB,nmat)
  ENDIF
 ENDDO
 Iatom = Jatom
 I = TENSOR%INDEX(IATOM,JATOM,1,1)
 IF(I.NE.0)THEN
    lsao => TENSOR%LSAO(I)
    dimA = lsao%nLocal(1)
    dimB = lsao%nLocal(2)
    call FULL_SYMelms_FROM_TRIANGULARelms2(lsao%elms,&
         &dimA,dimB,nmat)
 ENDIF
ENDDO
END SUBROUTINE LSTENSOR_FULL_SYMMAT_FROM_TRIANGULARMAT

subroutine FULL_SYMelms_FROM_TRIANGULARelms(elms2,elms1,dimenA,dimenB,nmat)
  implicit none
  integer :: dimenA,dimenB,nmat
  real(realk),intent(inout) :: elms2(dimenB,dimenA,nmat)
  real(realk),intent(in) :: elms1(dimenA,dimenB,nmat)
  !
  integer :: A,B,IMAT
  Do IMAT = 1,NMAT
     DO B=1,dimenB
        DO A=1,dimenA
           elms2(B,A,IMAT) = elms1(A,B,IMAT)
        ENDDO
     ENDDO
  ENDDO
END subroutine FULL_SYMELMS_FROM_TRIANGULARELMS

subroutine FULL_SYMelms_FROM_TRIANGULARelms2(elms,dimenA,dimenB,nmat)
  implicit none
  integer :: dimenA,dimenB,nmat
  real(realk),intent(inout) :: elms(dimenA,dimenB,nmat)
  !
  integer :: A,B,IMAT
  Do IMAT = 1,NMAT
     DO B=1,dimenB
        DO A=1,B-1
           elms(B,A,IMAT) = elms(A,B,IMAT)
        ENDDO
     ENDDO
  ENDDO
END subroutine FULL_SYMELMS_FROM_TRIANGULARELMS2

!> \brief set the lstensor to zero
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
SUBROUTINE lstensor_zero_lowertriangular(TENSOR)
implicit none
TYPE(LSTENSOR)     :: TENSOR
!
INTEGER    :: I,Ibat,Jbat,dim,ielms,iatom,jatom,dimA,dimB,Jbatch,Ibatch
TYPE(SLSAOTENSOR),pointer    :: lsao

IF(.NOT.TENSOR%SCREENTENSOR)CALL LSQUIT('LSTENSOR_ZERO_LOWERTRIANGULAR on for screen',-1)
!IF(ASSOCIATED(TENSOR%INDEX))THEN
!   DO Jatom = 1,TENSOR%natom(2)
!      Iatom = Jatom
!      I = TENSOR%INDEX(IATOM,JATOM,1,1)
!      IF(I.NE.0)THEN
!         lsao => TENSOR%SLSAO(I)
!         dimA = lsao%nLocal(1)
!         dimB = lsao%nLocal(2)
!         IF(dimA.NE.dimB)call lsquit('errornjbhjbhjbkjbhj',-1)
!         call set_lowerSelms_triangular_zero(lsao%selms,dimA,dimB)
!      ENDIF
!      DO Iatom = Jatom+1,TENSOR%natom(1)
!         I = TENSOR%INDEX(IATOM,JATOM,1,1)
!         IF(I.NE.0)THEN
!            lsao => TENSOR%SLSAO(I)
!            dim = lsao%nelms
!            call ls_sizero(lsao%selms,dim)
!         ENDIF
!      ENDDO
!   ENDDO
!ENDIF
!we do not need to touch the PS gab matrix

IF(ASSOCIATED(TENSOR%maxgab))THEN
   DO Jbatch = 1,TENSOR%nbatches(2)
      DO Ibatch = Jbatch+1,TENSOR%nbatches(1)
         TENSOR%maxgab(Ibatch,Jbatch) = shortzero
      ENDDO
   ENDDO
ENDIF

IF(ASSOCIATED(TENSOR%maxprimgab))THEN
   DO Jbatch = 1,TENSOR%nbatches(2)
      DO Ibatch = Jbatch+1,TENSOR%nbatches(1)
         TENSOR%maxprimgab(Ibatch,Jbatch) = shortzero
      ENDDO
   ENDDO
ENDIF

END SUBROUTINE LSTENSOR_ZERO_LOWERTRIANGULAR

!subroutine set_lowerSelms_triangular_zero(Selms,dimenA,dimenB)
!  implicit none
!  integer,intent(in) :: dimenA,dimenB
!  integer(kind=short) :: Selms(dimenA,dimenB)
!  !
!  integer :: A,B
!  DO B=1,dimenB
!     DO A=B+1,dimenA
!        Selms(A,B) = 0.0E0_realk
!     ENDDO
!  ENDDO
!end subroutine set_lowerSelms_triangular_zero

SUBROUTINE add_sublstensor_to_full_lstensor(TENSOR1,TENSOR2,&
     & natoms1,natoms2,natoms3,natoms4,atoms1,atoms2,atoms3,atoms4,&
     & nbast1,nbast2,nbast3,nbast4,SameFrag)
  implicit none
  TYPE(LSTENSOR)     :: TENSOR2
  TYPE(LSTENSOR),intent(in)     :: TENSOR1
  integer,intent(in) :: natoms1,natoms2,natoms3,natoms4
  integer,intent(in) :: nbast1,nbast2,nbast3,nbast4
  integer,intent(in) :: atoms1(natoms1),atoms2(natoms2)
  integer,intent(in) :: atoms3(natoms3),atoms4(natoms4)  
  logical,intent(in) :: SameFrag
!
  INTEGER    :: I1,Ibat,Jbat,Kbat,Lbat,nmat,dim,ielms,I2
  INTEGER    :: IATOM1,IATOM2,JATOM1,JATOM2,nbatI,ibatch1
  INTEGER    :: KATOM1,KATOM2,LATOM1,LATOM2,nbatJ,ibatch2
  INTEGER    :: jbatch1,jbatch2,k
  INTEGER    :: nAngmomA,nAngmomB,nAngmomC,nAngmomD
  integer(kind=long) :: nmemsize
  TYPE(LSAOTENSOR),pointer     :: lsao1,lsao2
  REAL(REALK),pointer          :: elms1(:),elms2(:)
  TYPE(SLSAOTENSOR),pointer    :: slsao1,slsao2
  integer(kind=short),pointer  :: selms1(:),selms2(:)

  IF(TENSOR1%Maggradienttensor)THEN
     TENSOR2%LSAO(1)%elms(1:3) = &
          & TENSOR2%LSAO(1)%elms(1:3) &
          & + TENSOR1%LSAO(1)%elms(1:3)
  ELSEIF(TENSOR1%gradienttensor)THEN
     dim = TENSOR1%ndim5*3
     DO Iatom1 = 1,natoms1
        Iatom2 = atoms1(Iatom1)
        DO ielms=1,dim
           TENSOR2%LSAO(Iatom2)%elms(ielms) = &
                & TENSOR2%LSAO(Iatom2)%elms(ielms) &
                & + TENSOR1%LSAO(Iatom1)%elms(ielms) 
        ENDDO
     ENDDO
     IF(.NOT.SameFrag)THEN
        DO Iatom1 = 1,natoms2
           I1 = nAtoms1 + Iatom1
           Iatom2 = atoms2(Iatom1)
           DO ielms=1,dim
              TENSOR2%LSAO(Iatom2)%elms(ielms) = &
                   & TENSOR2%LSAO(Iatom2)%elms(ielms) &
                   & + TENSOR1%LSAO(I1)%elms(ielms) 
           ENDDO
        ENDDO
        DO Iatom1 = 1,natoms3
           I1 = nAtoms1 + nAtoms2 + Iatom1
           Iatom2 = atoms3(Iatom1)
           DO ielms=1,dim
              TENSOR2%LSAO(Iatom2)%elms(ielms) = &
                   & TENSOR2%LSAO(Iatom2)%elms(ielms) &
                   & + TENSOR1%LSAO(I1)%elms(ielms) 
           ENDDO
        ENDDO
        DO Iatom1 = 1,natoms4
           I1 = nAtoms1 + nAtoms2 + nAtoms3 + Iatom1
           Iatom2 = atoms4(Iatom1)
           DO ielms=1,dim
              TENSOR2%LSAO(Iatom2)%elms(ielms) = &
                   & TENSOR2%LSAO(Iatom2)%elms(ielms) &
                   & + TENSOR1%LSAO(I1)%elms(ielms) 
           ENDDO
        ENDDO
     ENDIF
  ELSEIF(TENSOR1%pChargetensor)THEN
     nmat = TENSOR1%ndim5
     DO Iatom1 = 1,natoms1
        elms1 => TENSOR1%LSAO(Iatom1)%elms
!       elms2 => TENSOR2%LSAO(Iatom1)%elms
       elms2 => TENSOR2%LSAO(atoms1(Iatom1))%elms
        DO ielms=1,nmat
           elms2(ielms) = elms2(ielms) + elms1(ielms)
        ENDDO
     ENDDO
  ELSEIF(TENSOR1%Econtrib)THEN
     nmat = TENSOR1%ndim5
     elms1 => TENSOR1%LSAO(1)%elms
     elms2 => TENSOR2%LSAO(1)%elms
     DO ielms=1,nmat
        elms2(ielms) = elms2(ielms) + elms1(ielms)
     ENDDO
  ELSEIF(TENSOR1%screenTensor)THEN
     IF(ASSOCIATED(TENSOR2%SLSAO))THEN
      DO Jatom1 = 1,natoms2
       Jatom2 = atoms2(Jatom1)
       DO Iatom1 = 1,natoms1     !ATOMINDEX IN SUB
        Iatom2 = atoms1(Iatom1)  !ATOMINDEX IN FULL
        I2 = TENSOR2%INDEX(IATOM2,JATOM2,1,1) !FULL
        I1 = TENSOR1%INDEX(IATOM1,JATOM1,1,1) !SUB
        IF(I1.NE.0)THEN    
         IF(I2.EQ.0)call lsquit('error ADD_SUBLSTENSOR_TO_FULL_LSTENSOR',-1)
         slsao2 => TENSOR2%SLSAO(I2) 
         slsao1 => TENSOR1%SLSAO(I1) 
         dim = slsao1%nelms
         selms1 => slsao1%selms
         selms2 => slsao2%selms
         do ielms = 1,dim
            selms2(ielms) = MAX(selms2(ielms),selms1(ielms))
         enddo
        endif
       enddo
      enddo
     endif
     IF(ASSOCIATED(TENSOR2%maxgab))THEN
      DO Jatom1 = 1,natoms2       !ATOMINDEX IN SUB
       Jatom2 = atoms2(Jatom1)    !ATOMINDEX IN FULL
       nbatJ = TENSOR2%nAOBATCH(Jatom2,2) !nbat from FULL 
       DO Iatom1 = 1,natoms1     !ATOMINDEX IN SUB
        Iatom2 = atoms1(Iatom1)  !ATOMINDEX IN FULL
        nbatI = TENSOR2%nAOBATCH(Iatom2,1) !nbat from FULL 
        I2 = TENSOR2%INDEX(IATOM2,JATOM2,1,1) !FULL
        I1 = TENSOR1%INDEX(IATOM1,JATOM1,1,1) !SUB
        IF(I2.NE.0)THEN
         IF(I1.NE.0)THEN            
!            CALL LSQUIT('maxgab screening error in add_sublstensor_to_full_lstensor',-1)
!         ENDIF
          IBATCH2 = TENSOR2%SLSAO(I2)%AOBATCH(1) !IN FULL
          JBATCH2 = TENSOR2%SLSAO(I2)%AOBATCH(2)
          IBATCH1 = TENSOR1%SLSAO(I1)%AOBATCH(1) !IN SUB
          JBATCH1 = TENSOR1%SLSAO(I1)%AOBATCH(2)
          DO Jbat = 1,nbatJ
            DO Ibat = 1,nbatI
               TENSOR2%maxgab(ibatch2+ibat,jbatch2+jbat) = &
                    & TENSOR1%maxgab(ibatch1+ibat,jbatch1+jbat)
            enddo
          enddo
         ENDIF
        ENDIF
       enddo
      enddo
     ENDIF
     IF(ASSOCIATED(TENSOR2%maxprimgab))THEN
      DO Jatom1 = 1,natoms2       !ATOMINDEX IN SUB
       Jatom2 = atoms2(Jatom1)    !ATOMINDEX IN FULL
       nbatJ = TENSOR2%nAOBATCH(Jatom2,2) !nbat from FULL 
       DO Iatom1 = 1,natoms1     !ATOMINDEX IN SUB
        Iatom2 = atoms1(Iatom1)  !ATOMINDEX IN FULL
        nbatI = TENSOR2%nAOBATCH(Iatom2,1) !nbat from FULL 
        I2 = TENSOR2%INDEX(IATOM2,JATOM2,1,1) !FULL
        I1 = TENSOR1%INDEX(IATOM1,JATOM1,1,1) !SUB
        IF(I2.NE.0)THEN
         IF(I1.NE.0)THEN
            !CALL LSQUIT('maxgab screening error in add_sublstensor_to_full_lstensor',-1)
          IBATCH2 = TENSOR2%SLSAO(I2)%AOBATCH(1) !IN FULL
          JBATCH2 = TENSOR2%SLSAO(I2)%AOBATCH(2)
          IBATCH1 = TENSOR1%SLSAO(I1)%AOBATCH(1) !IN SUB
          JBATCH1 = TENSOR1%SLSAO(I1)%AOBATCH(2)
          DO Jbat = 1,nbatJ
            DO Ibat = 1,nbatI
               TENSOR2%maxprimgab(ibatch2+ibat,jbatch2+jbat) = &
                    & TENSOR1%maxprimgab(ibatch1+ibat,jbatch1+jbat)
            enddo
          enddo
         ENDIF
        ENDIF
       enddo
      enddo
     ENDIF
     IF(ASSOCIATED(TENSOR2%MBIE))THEN
      DO Jatom1 = 1,natoms2       !ATOMINDEX IN SUB
       Jatom2 = atoms2(Jatom1)    !ATOMINDEX IN FULL
       nbatJ = TENSOR2%nAOBATCH(Jatom2,2) !nbat from FULL 
       DO Iatom1 = 1,natoms1     !ATOMINDEX IN SUB
        Iatom2 = atoms1(Iatom1)  !ATOMINDEX IN FULL
        nbatI = TENSOR2%nAOBATCH(Iatom2,1) !nbat from FULL 
        I2 = TENSOR2%INDEX(IATOM2,JATOM2,1,1) !FULL
        I1 = TENSOR1%INDEX(IATOM1,JATOM1,1,1) !SUB
        IF(I2.NE.0)THEN
!         IF(I1.EQ.0)CALL LSQUIT('maxgab screening error in add_sublstensor_to_full_lstensor',-1)
         IF(I1.NE.0)THEN 
          IBATCH2 = TENSOR2%SLSAO(I2)%AOBATCH(1) !IN FULL
          JBATCH2 = TENSOR2%SLSAO(I2)%AOBATCH(2)
          IBATCH1 = TENSOR1%SLSAO(I1)%AOBATCH(1) !IN SUB
          JBATCH1 = TENSOR1%SLSAO(I1)%AOBATCH(2)
          DO Jbat = 1,nbatJ
           DO Ibat = 1,nbatI
            DO K = 1,TENSOR1%nMBIE
              TENSOR2%MBIE(K,ibatch2+ibat,jbatch2+jbat) = &
                    & TENSOR1%MBIE(K,ibatch1+ibat,jbatch1+jbat)
            enddo
           enddo
          enddo
         ENDIF
        ENDIF
       enddo
      enddo
     ENDIF
  ELSE
   nmat = TENSOR2%ndim5
   DO Latom1 = 1,natoms4
    Latom2 = atoms4(Latom1)
    DO Katom1 = 1,natoms3
     Katom2 = atoms3(Katom1)
     DO Jatom1 = 1,natoms2
      Jatom2 = atoms2(Jatom1)
      DO Iatom1 = 1,natoms1     !ATOMINDEX IN SUB
       Iatom2 = atoms1(Iatom1)  !ATOMINDEX IN FULL
       I2 = TENSOR2%INDEX(IATOM2,JATOM2,KATOM2,Latom2) !FULL
       I1 = TENSOR1%INDEX(IATOM1,JATOM1,KATOM1,Latom1) !SUB
       IF(I1.NE.0)THEN    
          IF(I2.EQ.0)call lsquit('error ADD_SUBLSTENSOR_TO_FULL_LSTENSOR',-1)
          lsao2 => TENSOR2%LSAO(I2) 
          lsao1 => TENSOR1%LSAO(I1) 
          dim = lsao1%nelms*nmat
          elms1 => lsao1%elms
          elms2 => lsao2%elms
          do ielms = 1,dim
             elms2(ielms) = elms2(ielms) + elms1(ielms)
          enddo
       endif
      enddo
     enddo
    enddo
   enddo
  endif

END SUBROUTINE ADD_SUBLSTENSOR_TO_FULL_LSTENSOR

SUBROUTINE add_lstensor_to_lstensor(TENSOR1,TENSOR2)
  implicit none
  TYPE(LSTENSOR)     :: TENSOR1,TENSOR2
!
  Integer :: ndim5,I,dim,ielms
  TYPE(LSAOTENSOR),pointer    :: lsao1,lsao2
  real(realk),pointer         :: elms1(:),elms2(:)
  ndim5 = TENSOR2%ndim5
  DO I = 1,TENSOR1%nLSAO
     lsao1 => TENSOR1%LSAO(I) 
     lsao2 => TENSOR2%LSAO(I) 
     dim = ndim5*lsao1%nelms
     elms1 => lsao1%elms
     elms2 => lsao2%elms
     do ielms = 1,dim
        elms2(ielms) = elms2(ielms) + elms1(ielms)         
     enddo
  enddo
END SUBROUTINE ADD_LSTENSOR_TO_LSTENSOR

SUBROUTINE init_MBIE_lstensor_5dim(TENSOR,AO1,AO2,ODscreen,lupri)
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(AOITEM),target  :: AO1,AO2
INTEGER            :: lupri
logical   :: ODscreen
!
Integer :: I,J,K,nbatches1,nbatches2
real(realk),parameter :: D0 = 0.0E0_realk 
integer(kind=long) :: nmemsize,nsize

IF(TENSOR%Screentensor)THEN
   !already initialised 
ELSE
   TENSOR%Screentensor = .TRUE.
ENDIF
TENSOR%nMBIE = 2 ! FIXME CURRENTLY HARDCODED TO 2
nbatches1=AO1%nbatches
nbatches2=AO2%nbatches
nullify(TENSOR%MBIE)
call mem_alloc(TENSOR%MBIE,2,nbatches1,nbatches2)
nsize = size(TENSOR%MBIE,KIND=long)*mem_realsize
call mem_allocated_mem_lstensor(nsize)

DO K = 1,nbatches2
   DO J = 1,nbatches1
      TENSOR%MBIE(1,J,K) = D0
      TENSOR%MBIE(2,J,K) = D0
   ENDDO
ENDDO
end SUBROUTINE init_MBIE_lstensor_5dim

SUBROUTINE init_gradientlstensor(TENSOR,natom1,natom2,natom3,natom4,samefrag,nmat,dimtensor1,Maggrad,lupri)
implicit none
TYPE(LSTENSOR)     :: TENSOR
INTEGER            :: nmat,lupri,natom1,natom2,natom3,natom4,dimtensor1
LOGICAL            :: sameFrag,maggrad
!
integer :: maxnatom,I,iatom,iatom1,iatom2,iatom3,iatom4,lstmem_index
type(lsaotensor),pointer :: lsao
integer(kind=long) :: nmemsize,AllocInt,AllocRealk,AllocIntS,nsize
call lstensor_nullify(TENSOR)
IF(SameFrag)THEN
   IF (nAtom2.NE.nAtom1) CALL LSQUIT('Error in init_gradientlstensor nAtoms inconsist!',lupri)
   IF (nAtom3.NE.nAtom1) CALL LSQUIT('Error in init_gradientlstensor nAtoms inconsist!',lupri)
   IF (nAtom4.NE.nAtom1) CALL LSQUIT('Error in init_gradientlstensor nAtoms inconsist!',lupri)
   maxnAtom=natom1
   NULLIFY(TENSOR%LSAO)
   NULLIFY(TENSOR%INDEX)
   CALL MEM_ALLOC(TENSOR%LSAO,natom1)
   CALL MEM_ALLOC(TENSOR%INDEX,natom1,1,1,4)
   nsize = size(TENSOR%INDEX,KIND=long)*mem_intsize
   call mem_allocated_mem_lstensor(nsize)
   DO I=1,4
      DO Iatom = 1,maxnAtom
         TENSOR%INDEX(Iatom,1,1,I) = Iatom
      ENDDO
   ENDDO
   TENSOR%natom(1) = nAtom1
   TENSOR%natom(2) = 0
   TENSOR%natom(3) = 0
   TENSOR%natom(4) = 0
ELSE
   maxnAtom=MAX(natom1,natom2,natom3,natom4)
   NULLIFY(TENSOR%LSAO)
   NULLIFY(TENSOR%INDEX)
   CALL MEM_ALLOC(TENSOR%LSAO,natom1+natom2+natom3+natom4)
   CALL MEM_ALLOC(TENSOR%INDEX,Maxnatom,1,1,4)
   nsize = size(TENSOR%INDEX,KIND=long)*mem_intsize
   call mem_allocated_mem_lstensor(nsize)
   DO Iatom1 = 1,natom1
      TENSOR%INDEX(Iatom1,1,1,1) = Iatom1
   ENDDO
   DO Iatom2 = 1,natom2
      TENSOR%INDEX(Iatom2,1,1,2) = natom1+Iatom2
   ENDDO
   DO Iatom3 = 1,natom3
      TENSOR%INDEX(Iatom3,1,1,3) = natom1+natom2+Iatom3
   ENDDO
   DO Iatom4 = 1,natom4
      TENSOR%INDEX(Iatom4,1,1,4) = natom1+natom2+natom3+Iatom4
   ENDDO
   TENSOR%natom(1) = nAtom1
   TENSOR%natom(2) = nAtom2
   TENSOR%natom(3) = nAtom3
   TENSOR%natom(4) = nAtom4
ENDIF

!maxnatom = MAX(TENSOR%natom(1),TENSOR%natom(2),TENSOR%natom(3),TENSOR%natom(4))
nullify(TENSOR%nAOBATCH)
!call mem_alloc(TENSOR%nAOBATCH(maxnatom,4))
!DO IAO = 1,4
!   DO IATOM = 1,maxnatom
!      TENSOR%nAOBATCH(IATOM,IAO) = 0
!   ENDDO
!ENDDO

AllocInt = 0
AllocInts = 0
IF(SameFrag)THEN
   AllocRealk = 3*nmat*natom1
ELSE
   AllocRealk = 3*nmat*(natom1+natom2+natom3+natom4)
ENDIF
call init_lstensorMem(AllocInt,AllocRealk,AllocInts,lstmem_index)
call zero_lstensorMem
TENSOR%lstmem_index = lstmem_index

TENSOR%pChargetensor = .FALSE.
TENSOR%Econtrib = .FALSE.
TENSOR%ScreenTensor = .FALSE.
TENSOR%gradienttensor = .TRUE.
TENSOR%Maggradienttensor = maggrad
TENSOR%ndim5 = nmat
DO Iatom = 1,natom1
   I=iatom
   lsao => TENSOR%LSAO(I) 
   lsao%ATOM(1) = Iatom
   lsao%ATOM(2) = 0
   lsao%ATOM(3) = 0
   lsao%ATOM(4) = 0
   lsao%AOBATCH = 0
   lsao%nLocal = 0
   lsao%maxBat = 0
   lsao%maxAng = 0
   NULLIFY(lsao%nOrb)
   NULLIFY(lsao%startLocalOrb)
   NULLIFY(lsao%startGlobalOrb)
   NULLIFY(lsao%nAngmom)
   lsao%nelms = 3
   NULLIFY(lsao%elms)
   CALL MEM_LSTPOINTER_ALLOC(lsao%elms,3*nmat)
ENDDO
TENSOR%nLSAO = natom1
IF(.NOT.SameFrag)THEN
   DO Iatom = 1,natom2
      I=iatom
      lsao => TENSOR%LSAO(nAtom1+I) 
      lsao%ATOM(1) = 0
      lsao%ATOM(2) = Iatom
      lsao%ATOM(3) = 0
      lsao%ATOM(4) = 0
      lsao%AOBATCH = 0
      lsao%nLocal = 0
      lsao%maxBat = 0
      lsao%maxAng = 0
      NULLIFY(lsao%nOrb)
      NULLIFY(lsao%startLocalOrb)
      NULLIFY(lsao%startGlobalOrb)
      NULLIFY(lsao%nAngmom)
      lsao%nelms = 3
      NULLIFY(lsao%elms)
      CALL MEM_LSTPOINTER_ALLOC(lsao%elms,3*nmat)
   ENDDO
   DO Iatom = 1,natom3
      I=iatom
      lsao => TENSOR%LSAO(nAtom1+nAtom2+I) 
      lsao%ATOM(1) = 0
      lsao%ATOM(2) = 0
      lsao%ATOM(3) = Iatom
      lsao%ATOM(4) = 0
      lsao%AOBATCH = 0
      lsao%nLocal = 0
      lsao%maxBat = 0
      lsao%maxAng = 0
      NULLIFY(lsao%nOrb)
      NULLIFY(lsao%startLocalOrb)
      NULLIFY(lsao%startGlobalOrb)
      NULLIFY(lsao%nAngmom)
      lsao%nelms = 3
      NULLIFY(lsao%elms)
      CALL MEM_LSTPOINTER_ALLOC(lsao%elms,3*nmat)
   ENDDO
   DO Iatom = 1,natom4
      I=iatom
      lsao => TENSOR%LSAO(nAtom1+nAtom2+nAtom3+I) 
      lsao%ATOM(1) = 0
      lsao%ATOM(2) = 0
      lsao%ATOM(3) = 0
      lsao%ATOM(4) = Iatom
      lsao%AOBATCH = 0
      lsao%nLocal = 0
      lsao%maxBat = 0
      lsao%maxAng = 0
      NULLIFY(lsao%nOrb)
      NULLIFY(lsao%startLocalOrb)
      NULLIFY(lsao%startGlobalOrb)
      NULLIFY(lsao%nAngmom)
      lsao%nelms = 3
      NULLIFY(lsao%elms)
      CALL MEM_LSTPOINTER_ALLOC(lsao%elms,3*nmat)
   ENDDO
   TENSOR%nLSAO = natom1+natom2+natom3+natom4
ENDIF

TENSOR%nbast(1) = dimtensor1
TENSOR%nbast(2) = TENSOR%nLSAO
TENSOR%nbast(3) = 1
TENSOR%nbast(4) = 1

end SUBROUTINE init_gradientlstensor

SUBROUTINE init_pcharge_lstensor(TENSOR,natom1,nmat,lupri)
implicit none
TYPE(LSTENSOR)     :: TENSOR
INTEGER            :: lupri,natom1,nmat

INTEGER :: nElms,lstmem_index
INTEGER :: IELM,IMAT2,I,J,Iatom,iatom1,iatom2,iatom3,iatom4
integer(kind=long) :: nmemsize,AllocInt,AllocRealk,AllocIntS,nsize
TYPE(LSAOTENSOR),pointer    :: lsao
REAL(REALK),pointer     :: elms(:)
call lstensor_nullify(tensor)
NULLIFY(TENSOR%LSAO)
NULLIFY(TENSOR%INDEX)
CALL MEM_ALLOC(TENSOR%LSAO,natom1)
CALL MEM_ALLOC(TENSOR%INDEX,natom1,1,1,1)
nsize = size(TENSOR%INDEX,KIND=long)*mem_intsize
call mem_allocated_mem_lstensor(nsize)
DO Iatom = 1,natom1
   TENSOR%INDEX(Iatom,1,1,1) = Iatom
ENDDO
TENSOR%natom(1) = nAtom1
TENSOR%natom(2) = 0
TENSOR%natom(3) = 0
TENSOR%natom(4) = 0

!maxnatom = MAX(TENSOR%natom(1),TENSOR%natom(2),TENSOR%natom(3),TENSOR%natom(4))
nullify(TENSOR%nAOBATCH)
!call mem_alloc(TENSOR%nAOBATCH(maxnatom,4))
!DO IAO = 1,4
!   DO IATOM = 1,maxnatom
!      TENSOR%nAOBATCH(IATOM,IAO) = 0
!   ENDDO
!ENDDO

TENSOR%pchargetensor = .TRUE.
TENSOR%Econtrib = .FALSE.
TENSOR%ScreenTensor = .FALSE.
TENSOR%Maggradienttensor = .FALSE.
TENSOR%gradienttensor = .FALSE.
TENSOR%nbast(1) = natom1
TENSOR%nbast(2) = 1
TENSOR%nbast(3) = 1
TENSOR%nbast(4) = 1
TENSOR%ndim5 = 1

AllocInt = 2*natom1
AllocInts = 0
AllocRealk = nmat*natom1
call init_lstensorMem(AllocInt,AllocRealk,AllocInts,lstmem_index)
call zero_lstensorMem
TENSOR%lstmem_index = lstmem_index

DO Iatom = 1,natom1
   I=iatom
   call lsaotensor_nullify(TENSOR%LSAO(I))
   lsao => TENSOR%LSAO(I) 
   lsao%ATOM(1) = Iatom
   lsao%ATOM(2) = 0
   lsao%ATOM(3) = 0
   lsao%ATOM(4) = 0
   lsao%AOBATCH(1) = 0
   lsao%AOBATCH(2) = 0
   lsao%AOBATCH(3) = 0
   lsao%AOBATCH(4) = 0
   lsao%nelms = nmat
   NULLIFY(lsao%elms)
   CALL MEM_LSTPOINTER_ALLOC(lsao%elms,nmat)
   lsao%nLocal(1) = 1
   lsao%nLocal(2) = 1
   lsao%nLocal(3) = 1
   lsao%nLocal(4) = 1
   lsao%maxBat = 1
   lsao%maxAng = 1
   CALL MEM_LSTPOINTER_ALLOC(lsao%startLocalOrb,2)
   lsao%startLocalOrb(1)=1
   lsao%startLocalOrb(2)=1
ENDDO
TENSOR%nLSAO = natom1

end SUBROUTINE init_pcharge_lstensor

subroutine build_grad_from_gradlstensor(TENSOR,GRAD,natom,nmat,lupri)
implicit none
INTEGER,intent(in):: lupri,nmat,natom
TYPE(LSTENSOR),intent(in):: TENSOR
REAL(REALK)       :: GRAD(3,natom*nmat) 
!
INTEGER :: nElms,nmat2
INTEGER :: IELM,IMAT,I,J,Iatom,iatom1,iatom2,iatom3,iatom4
integer :: natom1,natom2,natom3!,natom4
REAL(REALK),pointer     :: elms(:)

nmat2 = TENSOR%ndim5
IF(nmat2.NE.nmat)call lsquit('dim mismatch in build_grad_from_gradlstensor',-1)
Iatom = 0
natom1 = TENSOR%natom(1)
natom2 = TENSOR%natom(2)
natom3 = TENSOR%natom(3)
!natom4 = TENSOR%natom(4)
DO Iatom1 = 1,natom1
   elms => TENSOR%LSAO(Iatom+Iatom1)%elms
   DO IMAT = 1,nmat2
      DO J=1,3
         GRAD(J,Iatom1+(IMAT-1)*natom) = elms(J+(IMAT-1)*3)
      ENDDO
   ENDDO
ENDDO
Iatom = Iatom + TENSOR%natom(1)
DO Iatom2 = 1,TENSOR%natom(2)
   elms => TENSOR%LSAO(Iatom+Iatom2)%elms
   DO IMAT = 1,nmat2
      DO J=1,3
         GRAD(J,natom1+Iatom2+(IMAT-1)*natom) = elms(J+(IMAT-1)*3)
      ENDDO
   ENDDO
ENDDO
Iatom = Iatom + TENSOR%natom(2)
DO Iatom3 = 1,TENSOR%natom(3)
   elms => TENSOR%LSAO(Iatom+Iatom3)%elms
   DO IMAT = 1,nmat2
      DO J=1,3
         GRAD(J,natom1+natom2+Iatom3+(IMAT-1)*natom) = elms(J+(IMAT-1)*3)
      ENDDO
   ENDDO
ENDDO
Iatom = Iatom + TENSOR%natom(3)
DO Iatom4 = 1,TENSOR%natom(4)
   elms => TENSOR%LSAO(Iatom+Iatom4)%elms
   DO IMAT = 1,nmat2
      DO J=1,3
         GRAD(J,natom1+natom2+natom3+Iatom4+(IMAT-1)*natom) = elms(J+(IMAT-1)*3)
      ENDDO
   ENDDO
ENDDO

end subroutine build_grad_from_gradlstensor

subroutine build_gradlstensor_from_grad(TENSOR,GRAD,natom,nmat,lupri)
implicit none
TYPE(LSTENSOR)    :: TENSOR
INTEGER,intent(in):: lupri,natom,nmat
REAL(REALK)       :: GRAD(3,natom*nmat) 
!
INTEGER :: IMAT,J,Iatom
REAL(REALK),pointer  :: lstensorelms(:)

IF(natom.NE.TENSOR%natom(1))CALL LSQUIT('build_gradlstensor_from_grad',-1)
IF(TENSOR%ndim5.EQ. 1)THEN
   DO Iatom = 1,natom
      lstensorelms => TENSOR%LSAO(Iatom)%elms
      DO J=1,3
          lstensorelms(J) = GRAD(J,Iatom)
      ENDDO
   ENDDO
ELSE
   DO Iatom = 1,natom
      lstensorelms => TENSOR%LSAO(Iatom)%elms
      DO IMAT = 1,nmat
         DO J=1,3
             lstensorelms(J+(IMAT-1)*3) = GRAD(J,Iatom+(IMAT-1)*natom)
         ENDDO
      ENDDO
   ENDDO
ENDIF

end subroutine build_gradlstensor_from_grad

!> \brief build a lstensor from full 5 dim array
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT2 the full 5dim array
!> \param AO1 the Atomic orbital(AO) item for center 1 
!> \param AO2 the AO item for center 2 
!> \param AO3 the AO item for center 3 
!> \param AO4 the AO item for center 4 
!> \param nbast1 the size of the 1. dimension
!> \param nbast2 the size of the 2. dimension
!> \param nbast3 the size of the 3. dimension
!> \param nbast4 the size of the 4. dimension
!> \param nmat the number of density matrices or size of the 5. dimension
!> \param useAO1 flag to describe if the AO1 item should be used or an empty AO 
!> \param useAO2 flag to describe if the AO1 item should be used or an empty AO 
!> \param useAO3 flag to describe if the AO1 item should be used or an empty AO 
!> \param useAO4 flag to describe if the AO1 item should be used or an empty AO 
!> \param lupri the logical unit number for the output
SUBROUTINE Build_lstensor_from_full_5dim(TENSOR,MAT2,AO1,AO2,AO3,AO4,nbast1,nbast2,nbast3,nbast4,&
     &nmat,useAO1,useAO2,useAO3,useAO4,lupri)
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(AOITEM),target  :: AO1,AO2,AO3,AO4
TYPE(AOITEM),pointer :: AOT1,AOT2,AOT3,AOT4
INTEGER            :: nbast1,nbast2,nbast3,nbast4,nmat,lupri
REAL(REALK),target  :: MAT2(nbast1,nbast2,nbast3,nbast4,nmat)
logical   :: useAO1,useAO2,useAO3,useAO4
!
LOGICAL :: MATELMZERO
TYPE(LSTAUXITEM) :: LSTAUX
call init_lstensor_5dim(TENSOR,AO1,AO2,AO3,AO4,nbast1,nbast2,&
     &nbast3,nbast4,nmat,useAO1,useAO2,useAO3,useAO4,.FALSE.,.FALSE.,lupri)
call nullify_LSTAUXITEM(LSTAUX)
LSTAUX%Mat5dim => MAT2
MATELMZERO = .FALSE.
call BuildFromlstensor_aux(TENSOR,0,0,0,MATELMZERO,LSTAUX,elmsFrom5dim_FuncLstSub)
END SUBROUTINE Build_lstensor_from_full_5dim

!> \brief build a lstensor from full 3 dim array
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT2 the full 3dim array
!> \param AO1 the Atomic orbital(AO) item for center 1 
!> \param AO2 the AO item for center 2 
!> \param nbast1 the size of the 1. dimension
!> \param nbast2 the size of the 2. dimension
!> \param nmat the number of density matrices or size of the 5. dimension
!> \param useAO1 flag to describe if the AO1 item should be used or an empty AO 
!> \param useAO2 flag to describe if the AO1 item should be used or an empty AO 
!> \param lupri the logical unit number for the output
SUBROUTINE Build_lstensor_from_full_3dim(TENSOR,MAT2,AO1,AO2,nbast1,nbast2,nmat,useAO1,useAO2,lupri)
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(AOITEM),target :: AO1,AO2
INTEGER            :: nbast1,nbast2,nmat,lupri
REAL(REALK)        :: MAT2(nbast1,nbast2,nmat)
logical   :: useAO1,useAO2
!
TYPE(AOITEM),target :: AOT
TYPE(AOITEM),pointer :: AOT1,AOT2,AOT3,AOT4
logical :: useAO3,useAO4
call SET_EMPTY_AO(AOT)
AOT1 => AO1
AOT2 => AO2
AOT3 => AOT
AOT4 => AOT
IF(.NOT.useAO1) AOT1 => AOT
IF(.NOT.useAO2) AOT2 => AOT
useAO3=.TRUE.
useAO4=.TRUE.
call Build_lstensor_from_full_5dim(TENSOR,MAT2,AOT1,AOT2,AOT3,AOT4,&
     & nbast1,nbast2,1,1,nmat,useAO1,useAO2,useAO3,useAO4,lupri)

call FREE_EMPTY_AO(AOT)
END SUBROUTINE Build_lstensor_from_full_3dim

!> \brief build a lstensor from full 2 dim array
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR the lstensor
!> \param MAT2 the full 3dim array
!> \param AO1 the Atomic orbital(AO) item for center 1 
!> \param AO2 the AO item for center 2 
!> \param nbast1 the size of the 1. dimension
!> \param nbast2 the size of the 2. dimension
!> \param useAO1 flag to describe if the AO1 item should be used or an empty AO 
!> \param useAO2 flag to describe if the AO1 item should be used or an empty AO 
!> \param lupri the logical unit number for the output
SUBROUTINE Build_lstensor_from_full_2dim(TENSOR,MAT2,AO1,AO2,nbast1,nbast2,useAO1,useAO2,lupri)
implicit none
TYPE(LSTENSOR)     :: TENSOR
TYPE(AOITEM),target :: AO1,AO2
INTEGER            :: nbast1,nbast2,lupri
REAL(REALK)        :: MAT2(nbast1,nbast2)
logical   :: useAO1,useAO2
!
logical   :: useAO3,useAO4
TYPE(AOITEM),pointer :: AOT1,AOT2,AOT3,AOT4
TYPE(AOITEM),target :: AOT
integer :: nmat
nmat=1
call SET_EMPTY_AO(AOT)
AOT1 => AO1
AOT2 => AO2
AOT3 => AOT
AOT4 => AOT
IF(.NOT.useAO1) AOT1 => AOT
IF(.NOT.useAO2) AOT2 => AOT
useAO3=.TRUE.
useAO4=.TRUE.
call Build_lstensor_from_full_5dim(TENSOR,MAT2,AOT1,AOT2,AOT3,AOT4,&
     & nbast1,nbast2,1,1,nmat,useAO1,useAO2,useAO3,useAO4,lupri)
call FREE_EMPTY_AO(AOT)
END SUBROUTINE Build_lstensor_from_full_2dim

SUBROUTINE add_full_2dim_to_lstensor(TENSOR,fullMAT,dim1,dim2,ndmat)
implicit none
TYPE(LSTENSOR)     :: TENSOR
real(realk),target :: fullMAT(:,:)
integer,intent(in) :: dim1,dim2,ndmat
!
TYPE(LSTAUXITEM) :: LSTAUX
INTEGER    :: I,J,DIMI,DIMJ,SUM1,SUM2,imat,nmat,nbast(4),natom,iatom
LOGICAL :: MATELMZERO
IF(ndmat.NE.TENSOR%ndim5)Call lsquit('dim mismatch add_full_2dim_to_lstensor',-1)
nbast=TENSOR%nbast
DO I=1,4
   IF(dim1 .EQ.nbast(I))THEN
      DIMI = I
      EXIT
   ENDIF
ENDDO
DO J=1,4
   IF(J.NE.DIMI)THEN
      IF(dim2 .EQ.nbast(J))THEN
         DIMJ = J
         EXIT
      ENDIF
   ENDIF
ENDDO
SUM1=nbast(DIMI)+nbast(DIMJ)
SUM2=dim1+dim2
IF(SUM1.NE.SUM2)THEN
   print*,'dim1',DIMI
   print*,'dim2',DIMJ
   print*,'dim1',dim1
   print*,'dim2',dim2
   print*,'nbast1',nbast(1)
   print*,'nbast2',nbast(2)
   print*,'nbast3',nbast(3)
   print*,'nbast4',nbast(4)
   print*,'nbast',nbast(DIMI)
   print*,'nbast',nbast(DIMJ)
   TENSOR%index(99,99,99,99) = -99
   CALL lsQUIT('ERROR: add_full_2dim_to_lstensor dim do not match',-1)
ENDIF
IF(TENSOR%gradienttensor)THEN
 natom = TENSOR%natom(1)
 DO Iatom = 1,natom
  DO IMAT = 1,ndmat
   DO J=1,3
    TENSOR%LSAO(Iatom)%elms(J+(IMAT-1)*3)=TENSOR%LSAO(Iatom)%elms(J+(IMAT-1)*3)&
         & + fullmat(J,Iatom+(IMAT-1)*natom)
   ENDDO
  ENDDO
 ENDDO
ELSE
   call nullify_LSTAUXITEM(LSTAUX)
   LSTAUX%Mat2dim => fullmat
   MATELMZERO=.FALSE.
   call BuildFromlstensor_aux(TENSOR,DIMI,DIMJ,0,MATELMZERO,LSTAUX,add2dim_FuncLstSub)
ENDIF

END SUBROUTINE ADD_FULL_2DIM_TO_LSTENSOR

!> \brief determine if the batch should be screened away  
!> \author T. Kjaergaard
!> \date 2010
!> \param nbatI number of batch for batch1
!> \param nbatJ number of batch for batch2
!> \param nbatK number of batch for batch3
!> \param nbatL number of batch for batch4
!> \param AOT1batch AO batch index 1
!> \param AOT2batch AO batch index 2
!> \param AOT3batch AO batch index 3
!> \param AOT4batch AO batch index 4
!> \param useAO1 use the AObatch 1 in contrast to use an empty 
!> \param useAO2 use the AObatch 2 in contrast to use an empty 
!> \param useAO3 use the AObatch 3 in contrast to use an empty 
!> \param useAO4 use the AObatch 4 in contrast to use an empty 
!> \param BATCHA3 the 1. AO batch 
!> \param BATCHB3 the 2. AO batch 
!> \param BATCHC3 the 3. AO batch 
!> \param BATCHD3 the 4. AO batch 
!> \param OVERALLscreen should the batch be screened
subroutine determineODscreening(nbatI,nbatJ,nbatK,nbatL,&
     & AOT1batch,AOT2batch,AOT3batch,AOT4batch,useAO1,useAO2,&
     & useAO3,useAO4,BATCHA3,BATCHB3,BATCHC3,BATCHD3,OVERALLscreen)
  use OD_Type, only: getODscreening
implicit none
integer,intent(in) :: nbatI,nbatJ,nbatK,nbatL,AOT1batch,AOT2batch,AOT3batch,AOT4batch
TYPE(AOBATCH),intent(in)  :: BATCHA3(:)
TYPE(AOBATCH),intent(in)  :: BATCHB3(:)
TYPE(AOBATCH),intent(in)  :: BATCHC3(:)
TYPE(AOBATCH),intent(in)  :: BATCHD3(:)
logical,intent(in) :: useAO1,useAO2,useAO3,useAO4
logical,intent(out) :: overallscreen
!
integer :: batAOT1,batAOT2,batAOT3,batAOT4
integer :: Ibat,Jbat,Kbat,Lbat
logical :: ODscreenLHS,ODscreenRHS

OVERALLscreen = .TRUE.

IF(useAO1 .OR. useAO2)THEN
 batAOT1 = AOT1batch
 DO Ibat = 1,nbatI
  batAOT1 = batAOT1+1 
  batAOT2 = AOT2batch
  DO Jbat = 1,nbatJ
   batAOT2 = batAOT2+1 
   ODscreenLHS = .FALSE.
   call getODscreening(BATCHA3(batAOT1),BATCHB3(batAOT2),ODscreenLHS) 
   IF(.NOT.ODscreenLHS)THEN
 !   print*,'ODscreenLHS',ODscreenLHS,'is false so we should calc'
    !if just one element is not screen - we do not screen any
    OVERALLscreen = .FALSE.
    EXIT
   ENDIF
  ENDDO
 ENDDO
ENDIF

IF(useAO3 .OR. useAO4)THEN
 IF(OVERALLscreen)THEN
  batAOT3 = AOT3batch
  DO Kbat = 1,nbatK
   batAOT3 = batAOT3 + 1 
   batAOT4 = AOT4batch
   DO Lbat = 1,nbatL
    batAOT4 = batAOT4 + 1  
    ODscreenRHS = .FALSE.
    call getODscreening(BATCHC3(batAOT3),BATCHD3(batAOT4),ODscreenRHS) 
    IF(.NOT.ODscreenRHS)THEN
     !if just one element is not screen - we do not screen any
     OVERALLscreen = .FALSE.
     EXIT
    ENDIF
   ENDDO
  ENDDO
 ENDIF
ENDIF

end subroutine determineODscreening

!> \brief 
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR1 the original full slstensor
!> \param TENSOR2 the new batch slstensor
SUBROUTINE build_BatchGab(AOfull1,AOfull2,AO1,AO2,iBatch1,jBatch1,&
     & ibatchStart,jbatchStart,dim1,dim2,TENSOR1,TENSOR2)
  implicit none
  integer,intent(in) :: dim1,dim2
  type(AOITEM),intent(in) :: AOfull1,AOfull2,AO1,AO2
  TYPE(LSTENSOR),intent(in)     :: TENSOR1
  TYPE(LSTENSOR),intent(inout)  :: TENSOR2
  integer,intent(in)   :: iBatchStart,jBatchStart
  integer,intent(in)   :: iBatch1,jBatch1
  !
  INTEGER    :: I1,I2,dim,ielms,nbatches1,nbatches2
  INTEGER    :: nAngmomA,nAngmomB,nAngmomC,nAngmomD
  integer(kind=long) :: nmemsize,nsize
  TYPE(SLSAOTENSOR),pointer    :: lsao1
  TYPE(SLSAOTENSOR),pointer    :: lsao2
  integer(kind=short),pointer     :: elms1(:),elms2(:)
  integer(kind=short) :: maxelm1
  integer :: natomJ,jatom_old,jBatch2,jatom,natomI,JbatchLoc,maxBat1
  integer :: iatom_old,iBatch2,iatom,jatom2,iatom2,IbatchLoc,maxBat2
  integer :: NIBAT,NJBAT,ibat,jbat,iBatchLoc2,jBatchLoc2,ibatch,jbatch
  logical :: CS_SCREEN,PS_SCREEN
  integer :: sA1,sA2,sB1,sB2,nOrbA1,nOrbA2,nOrbB1,nOrbB2,I,J,nA2,nB2,nA1,nB1
  CS_SCREEN = ASSOCIATED(TENSOR1%maxgab)
  PS_SCREEN = ASSOCIATED(TENSOR1%maxprimgab)
  call lstensor_nullify(TENSOR2)
  call init_ps_lstensor(TENSOR2,AO1,AO2,dim1,dim2,.FALSE.,6)

!  print*,'the FULL LSTENSOR'
!  call lstensor_print(TENSOR1,6)
!  print*,'the sub LSTENSOR'
!  call lstensor_print(TENSOR2,6)
!  print*,'the full AO'
!  call PRINT_AO(6,AOfull)
!  print*,'the AO1'
!  call PRINT_AO(6,AO1)
!  print*,'the AO2'
!  call PRINT_AO(6,AO2)
!  print*,'ibatchStart,jbatchStart',ibatchStart,jbatchStart

  nbatches1 = TENSOR2%nbatches(1)
  nbatches2 = TENSOR2%nbatches(2)
  dim = nbatches2*nbatches1
  IF(CS_SCREEN)THEN
     nullify(TENSOR2%maxgab)
     call mem_alloc(TENSOR2%maxgab,nbatches1,nbatches2)
     nsize = size(TENSOR2%maxgab,KIND=long)*mem_shortintsize
     call mem_allocated_mem_lstensor(nsize)

     call ls_sizero(TENSOR2%maxgab,dim)
  ENDIF
  IF(PS_SCREEN)THEN
     nullify(TENSOR2%maxprimgab)
     call mem_alloc(TENSOR2%maxprimgab,nbatches1,nbatches2)
     nsize = size(TENSOR2%maxprimgab,KIND=long)*mem_shortintsize
     call mem_allocated_mem_lstensor(nsize)
     call ls_sizero(TENSOR2%maxprimgab,dim)
  ENDIF

  Jatom2 = 0
  jatom_old = 0
  jBatch = jBatchStart
  jBatchLoc2 = 0
  DO jBatch2=1,nbatches2
     jBatch = jBatch + 1
     jBatchLoc2 = jBatchLoc2 + 1
     jatom = AOfull2%Batch(jBatch)%atom
     IF(jatom.NE.jatom_old)then
        jatom_old = jatom
        jatom2 = Jatom2 +1 
        jBatchLoc2 = 1
     ENDIF
     JbatchLoc = AOfull2%Batch(jBatch)%batch
!     IF(JbatchLoc2.NE.AO2%Batch(JBatch2)%batch)CALL LSQUIT('batchlocB',-1)
     
     Iatom2 = 0
     iatom_old = 0
     iBatch = iBatchStart
     iBatchLoc2 = 0
     DO iBatch2=1,nbatches1
        iBatch = iBatch + 1
        iBatchLoc2 = iBatchLoc2 + 1
        iatom = AOfull1%Batch(iBatch)%atom
        IF(iatom.NE.iatom_old)then
           iatom_old = iatom
           Iatom2 = Iatom2 +1 
           iBatchLoc2 = 1
        ENDIF
        IbatchLoc = AOfull1%Batch(iBatch)%batch
!        IF(IbatchLoc2.NE.AO1%Batch(iBatch2)%batch)CALL LSQUIT('batchlocA',-1)
        I1 = TENSOR1%INDEX(iatom,jatom,1,1)
        I2 = TENSOR2%INDEX(iatom2,jatom2,1,1)
        IF(I1.NE.0)THEN
           IF(I2.EQ.0)call lsquit('I1,I2 mismatch',-1)
           maxBat1 = TENSOR1%SLSAO(I1)%maxBat
           maxBat2 = TENSOR2%SLSAO(I2)%maxBat
           nOrbA1 = TENSOR1%SLSAO(I1)%nOrb(ibatchLoc)
           !        nOrbA2 = TENSOR2%SLSAO(I2)%nOrb(ibatchLoc2,1)
           !        IF(nOrbA1.NE.nOrbA2)call lsquit('errefbdbndb',-1)
           nOrbB1 = TENSOR1%SLSAO(I1)%nOrb(jbatchLoc+maxBat1)
           !        nOrbB2 = TENSOR2%SLSAO(I2)%nOrb(jbatchLoc2,2)
           !        IF(nOrbB1.NE.nOrbB2)call lsquit('errefbdbndb',-1)
           sA1 = TENSOR1%SLSAO(I1)%startLocalOrb(ibatchLoc)-1
           sA2 = TENSOR2%SLSAO(I2)%startLocalOrb(ibatchLoc2)-1
           sB1 = TENSOR1%SLSAO(I1)%startLocalOrb(jbatchLoc+maxBat1)-1
           sB2 = TENSOR2%SLSAO(I2)%startLocalOrb(jbatchLoc2+maxBat2)-1
           nA2 = TENSOR2%SLSAO(I2)%nLocal(1)
           nB2 = TENSOR2%SLSAO(I2)%nLocal(2)
           nA1 = TENSOR1%SLSAO(I1)%nLocal(1)
           nB1 = TENSOR1%SLSAO(I1)%nLocal(2)
           IF(Associated(TENSOR1%SLSAO(I1)%selms))THEN
              call copy_batch_selms(TENSOR2%SLSAO(I2)%selms,sA2,sB2,nA2,nB2,TENSOR1%SLSAO(I1)%selms,sA1,sB1,nA1,nB1,nOrbA1,nOrbB1)
           ELSE
              NULLIFY(TENSOR2%SLSAO(I2)%selms)
              TENSOR2%SLSAO(I2)%nelms=0
           ENDIF
           IF(CS_SCREEN)TENSOR2%maxgab(iBatch2,jBatch2) = TENSOR1%maxgab(iBatch,jBatch) 
           IF(PS_SCREEN)TENSOR2%maxprimgab(iBatch2,jBatch2) = TENSOR1%maxprimgab(iBatch,jBatch) 
        ENDIF
     ENDDO
  ENDDO
  IF(CS_SCREEN)call set_lst_maxgabelms(TENSOR2)
  IF(PS_SCREEN)call set_lst_maxprimgabelms(TENSOR2)
  
!  print*,'the FULL LSTENSOR AFTER build_BatchGab'
!  call lstensor_print(TENSOR1,6)
!  print*,'the sub LSTENSOR AFTER build_BatchGab'
!  call lstensor_print(TENSOR2,6)

end SUBROUTINE build_BatchGab

!!$SUBROUTINE build_BatchGab_from_full(TENSOR2,batchA,batchB,batchsizeA,batchSizeB,FullGab,nbatches)
!!$  implicit none
!!$  INTEGER,intent(in)            :: batchA,batchB
!!$  INTEGER,intent(in)            :: batchsizeA,batchSizeB  
!!$  integer(kind=short)           :: FullGab(nbatches,nbatches)
!!$  !
!!$  logical :: CS_SCREEN,PS_SCREEN
!!$  integer :: I,J
!!$  CS_SCREEN = ASSOCIATED(TENSOR1%maxgab)
!!$  PS_SCREEN = .FALSE.
!!$  call lstensor_nullify(TENSOR2)
!!$  TENSOR2%nSLSAO = 0
!!$  TENSOR2%pChargetensor = .FALSE.
!!$  TENSOR2%Econtrib = .FALSE.
!!$  TENSOR2%Screentensor = .TRUE.
!!$  TENSOR2%Screenoutput = .FALSE.
!!$  TENSOR2%MagGradienttensor = .FALSE.
!!$  TENSOR2%Gradienttensor = .FALSE.
!!$  NULLIFY(TENSOR%SLSAO)
!!$  NULLIFY(TENSOR%INDEX)
!!$  TENSOR2%natom(1) = -1 !not used
!!$  TENSOR2%natom(2) = -1
!!$  TENSOR2%natom(3) = 1
!!$  TENSOR2%natom(4) = 1
!!$  TENSOR2%ndim5 = 1
!!$
!!$  TENSOR2%Screentensor = .TRUE.
!!$  TENSOR2%nbatches(1) = batchsizeA
!!$  TENSOR2%nbatches(2) = batchSizeB  
!!$  TENSOR2%nbatches(3) = 0
!!$  TENSOR2%nbatches(4) = 0
!!$  dim = batchsizeA*batchsizeB
!!$  IF(CS_SCREEN)THEN
!!$     nullify(TENSOR2%maxgab)
!!$     call mem_alloc(TENSOR2%maxgab,batchsizeA,batchsizeB)
!!$     nsize = size(TENSOR2%maxgab,KIND=long)*mem_shortintsize
!!$     call mem_allocated_mem_lstensor(nsize)
!!$  ENDIF
!!$  nullify(TENSOR2%maxprimgab)
!!$  DO J=1,batchsizeB
!!$     DO I=1,batchsizeA
!!$        TENSOR2%maxgab(I,J) = FullGab(batchA+I-1,batchB+J-1)
!!$     ENDDO
!!$  ENDDO
!!$  IF(CS_SCREEN)call set_lst_maxgabelms(TENSOR2)  
!!$!  print*,'the sub LSTENSOR AFTER build_BatchGab_From_FULL'
!!$!  call lstensor_print(TENSOR2,6)
!!$end SUBROUTINE build_BatchGab_From_FULL

!> \brief 
!> \author T. Kjaergaard
!> \date 2010
!> \param TENSOR1 the original full slstensor
!> \param TENSOR2 the new batch slstensor
SUBROUTINE build_BatchGabK(AOfull,AO1,iBatch1,ibatchStart,dim1,TENSOR1,TENSOR2)
  implicit none
  integer,intent(in) :: dim1
  type(AOITEM),intent(in) :: AOfull,AO1
  TYPE(LSTENSOR),intent(in)     :: TENSOR1
  TYPE(LSTENSOR),intent(inout)  :: TENSOR2
  integer,intent(in)   :: iBatchStart
  integer,intent(in)   :: iBatch1
  !
  INTEGER    :: I1,I2,dim,ielms,nbatches1,nbatches2
  INTEGER    :: nAngmomA,nAngmomB,nAngmomC,nAngmomD
  integer(kind=long) :: nmemsize,nsize
  TYPE(SLSAOTENSOR),pointer    :: lsao1
  TYPE(SLSAOTENSOR),pointer    :: lsao2
  integer(kind=short),pointer     :: elms1(:),elms2(:)
  integer(kind=short) :: maxelm1
  integer :: natomJ,jatom_old,jBatch2,jatom,natomI,JbatchLoc,maxBat1
  integer :: iatom_old,iBatch2,iatom,jatom2,iatom2,IbatchLoc,maxBat2
  integer :: NIBAT,NJBAT,ibat,jbat,iBatchLoc2,jBatchLoc2,ibatch,jbatch
  logical :: CS_SCREEN,PS_SCREEN
  integer :: sA1,sA2,sB1,sB2,nOrbA1,nOrbA2,nOrbB1,nOrbB2,I,J,nA2,nB2,nA1,nB1
  CS_SCREEN = ASSOCIATED(TENSOR1%maxgab)
  PS_SCREEN = ASSOCIATED(TENSOR1%maxprimgab)
  call lstensor_nullify(TENSOR2)
  call init_ps_lstensor(TENSOR2,AO1,AOfull,dim1,AOfull%nbast,.FALSE.,6)

!  print*,'the FULL LSTENSOR'
!  call lstensor_print(TENSOR1,6)
!  print*,'the sub LSTENSOR'
!  call lstensor_print(TENSOR2,6)
!  print*,'the full AO'
!  call PRINT_AO(6,AOfull)
!  print*,'the AO1'
!  call PRINT_AO(6,AO1)
!  print*,'the AO2'
!  call PRINT_AO(6,AO2)
!  print*,'ibatchStart,jbatchStart',ibatchStart,jbatchStart

  nbatches1 = TENSOR2%nbatches(1)
  nbatches2 = TENSOR2%nbatches(2)
  dim = nbatches2*nbatches1
  IF(CS_SCREEN)THEN
     nullify(TENSOR2%maxgab)
     call mem_alloc(TENSOR2%maxgab,nbatches1,nbatches2)
     nsize = size(TENSOR2%maxgab,KIND=long)*mem_shortintsize
     call mem_allocated_mem_lstensor(nsize)
     call ls_sizero(TENSOR2%maxgab,dim)
  ENDIF
  IF(PS_SCREEN)THEN
     nullify(TENSOR2%maxprimgab)
     call mem_alloc(TENSOR2%maxprimgab,nbatches1,nbatches2)
     nsize = size(TENSOR2%maxprimgab,KIND=long)*mem_shortintsize
     call mem_allocated_mem_lstensor(nsize)
     call ls_sizero(TENSOR2%maxprimgab,dim)
  ENDIF

  jBatch = 0
  DO jBatch2=1,nbatches2
     jBatch = jBatch + 1
     jatom = AOfull%Batch(jBatch)%atom
     JbatchLoc = AOfull%Batch(jBatch)%batch
     
     Iatom2 = 0
     iatom_old = 0
     iBatch = iBatchStart
     iBatchLoc2 = 0
     DO iBatch2=1,nbatches1
        iBatch = iBatch + 1
        iBatchLoc2 = iBatchLoc2 + 1
        iatom = AOfull%Batch(iBatch)%atom
        IF(iatom.NE.iatom_old)then
           iatom_old = iatom
           Iatom2 = Iatom2 +1 
           iBatchLoc2 = 1
        ENDIF
        IbatchLoc = AOfull%Batch(iBatch)%batch
!        IF(IbatchLoc2.NE.AO1%Batch(iBatch2)%batch)CALL LSQUIT('batchlocA',-1)
        I1 = TENSOR1%INDEX(iatom,jatom,1,1)
        I2 = TENSOR2%INDEX(iatom2,jatom,1,1)
        IF(I1.NE.0)THEN
           IF(I2.EQ.0)call lsquit('I1,I2 mismatch',-1)
           maxBat1 = TENSOR1%SLSAO(I1)%maxBat
           maxBat2 = TENSOR2%SLSAO(I2)%maxBat
           nOrbA1 = TENSOR1%SLSAO(I1)%nOrb(ibatchLoc)
           nOrbB1 = TENSOR1%SLSAO(I1)%nOrb(jbatchLoc+maxBat1)
           sA1 = TENSOR1%SLSAO(I1)%startLocalOrb(ibatchLoc)-1
           sA2 = TENSOR2%SLSAO(I2)%startLocalOrb(ibatchLoc2)-1
           sB1 = TENSOR1%SLSAO(I1)%startLocalOrb(jbatchLoc+maxBat1)-1
           sB2 = TENSOR2%SLSAO(I2)%startLocalOrb(jbatchLoc+maxBat2)-1
           nA2 = TENSOR2%SLSAO(I2)%nLocal(1)
           nB2 = TENSOR2%SLSAO(I2)%nLocal(2)
           nA1 = TENSOR1%SLSAO(I1)%nLocal(1)
           nB1 = TENSOR1%SLSAO(I1)%nLocal(2)
           IF(Associated(TENSOR1%SLSAO(I1)%selms))THEN
              call copy_batch_selms(TENSOR2%SLSAO(I2)%selms,sA2,sB2,nA2,nB2,TENSOR1%SLSAO(I1)%selms,sA1,sB1,nA1,nB1,nOrbA1,nOrbB1)
           ELSE
              NULLIFY(TENSOR2%SLSAO(I2)%selms)
              TENSOR2%SLSAO(I2)%nelms=0
           ENDIF
           IF(CS_SCREEN)TENSOR2%maxgab(iBatch2,jBatch) = TENSOR1%maxgab(iBatch,jBatch) 
           IF(PS_SCREEN)TENSOR2%maxprimgab(iBatch2,jBatch) = TENSOR1%maxprimgab(iBatch,jBatch) 
        ENDIF
     ENDDO
  ENDDO
  IF(CS_SCREEN)call set_lst_maxgabelms(TENSOR2)
  IF(PS_SCREEN)call set_lst_maxprimgabelms(TENSOR2)
  
!  print*,'the FULL LSTENSOR AFTER build_BatchGabK'
!  call lstensor_print(TENSOR1,6)
!  print*,'the sub LSTENSOR AFTER build_BatchGabK'
!  call lstensor_print(TENSOR2,6)

end SUBROUTINE build_BatchGabK

subroutine copy_batch_selms(selms2,sA2,sB2,nA2,nB2,selms1,sA1,sB1,nA1,nB1,nOrbA1,nOrbB1)
  implicit none
  integer :: sA2,sB2,nA2,nB2,sA1,sB1,nA1,nB1,nOrbA1,nOrbB1 
  integer(kind=short),intent(inout) :: selms2(nA2,nB2)
  integer(kind=short),intent(in)    :: selms1(nA1,nB1)
!
  integer :: I,J
  DO J = 1,nOrbB1
     DO I = 1,nOrbA1
        selms2(sA2+I,sB2+J)= selms1(sA1+I,sB1+J)              
     ENDDO
  ENDDO
end subroutine copy_batch_selms

!!$!> \brief copy an lstensor to a new lstensor
!!$!> \author T. Kjaergaard
!!$!> \date 2010
!!$!> \param TENSOR1 the original full slstensor
!!$!> \param TENSOR2 the new batch slstensor
!!$SUBROUTINE build_singleBatchGab(AOfull,AO1,ibatch,nbatches,dim2,TENSOR1,TENSOR_LHS,TENSOR_RHS)
!!$  implicit none
!!$  integer,intent(in) :: dim2
!!$  type(AOITEM),intent(in) :: AOfull
!!$  type(AOITEM) :: AO1
!!$  TYPE(LSTENSOR),intent(in)     :: TENSOR1    ! full,full
!!$  TYPE(LSTENSOR),intent(inout)  :: TENSOR_LHS ! full,batch
!!$  TYPE(LSTENSOR),intent(inout)  :: TENSOR_RHS ! batch,full
!!$  integer,intent(in)   :: iBatch
!!$  !
!!$  INTEGER    :: I1,I2,dim,ielms,nbatches1,nbatches2,nbast
!!$  INTEGER    :: nAngmomA,nAngmomB,nAngmomC,nAngmomD
!!$  integer(kind=long) :: nmemsize
!!$  TYPE(SLSAOTENSOR),pointer    :: lsao1
!!$  TYPE(SLSAOTENSOR),pointer    :: lsao2
!!$  integer(kind=short),pointer     :: elms1(:),elms2(:)
!!$  integer(kind=short) :: maxelm1
!!$  integer :: natomJ,jatom_old,jBatch2,jatom,natomI,JbatchLoc,nbatches
!!$  integer :: iatom_old,iBatch2,iatom,jatom2,iatom2,IbatchLoc,I2_LHS,I2_RHS
!!$  integer :: NIBAT,NJBAT,ibat,jbat,iBatchLoc2,jBatchLoc2,jbatch,maxbat1,maxbat2
!!$  logical :: CS_SCREEN,PS_SCREEN
!!$  integer :: sA1,sA2,sB1,sB2,nOrbA1,nOrbA2,nOrbB1,nOrbB2,I,J,nA2,nB2,nA1,nB1
!!$  CS_SCREEN = ASSOCIATED(TENSOR1%maxgab)
!!$  PS_SCREEN = ASSOCIATED(TENSOR1%maxprimgab)
!!$  call lstensor_nullify(TENSOR_LHS)
!!$  call lstensor_nullify(TENSOR_RHS)
!!$  nbast = TENSOR1%nbast(1)
!!$  call init_ps_lstensor(TENSOR_LHS,AOfull,AO1,nbast,dim2,.FALSE.,6)
!!$  call init_ps_lstensor(TENSOR_RHS,AO1,AOfull,dim2,nbast,.FALSE.,6)
!!$  call LSTENSOR_mem_est(TENSOR_LHS,nmemsize)
!!$  call remove_mem_from_global(nmemsize)
!!$  call LSTENSOR_mem_est(TENSOR_RHS,nmemsize)
!!$  call remove_mem_from_global(nmemsize)
!!$
!!$  nbatches = TENSOR1%nbatches(1)
!!$  dim = nbatches
!!$  IF(CS_SCREEN)THEN
!!$     nullify(TENSOR_LHS%maxgab)
!!$     call mem_alloc(TENSOR_LHS%maxgab,nbatches,1)
!!$     call ls_sizero(TENSOR_LHS%maxgab,dim)
!!$     nullify(TENSOR_RHS%maxgab)
!!$     call mem_alloc(TENSOR_RHS%maxgab,1,nbatches)
!!$     call ls_sizero(TENSOR_RHS%maxgab,dim)
!!$  ENDIF
!!$  IF(PS_SCREEN)THEN
!!$     nullify(TENSOR_LHS%maxprimgab)
!!$     call mem_alloc(TENSOR_LHS%maxprimgab,nbatches,1)
!!$     call ls_sizero(TENSOR_LHS%maxprimgab,dim)
!!$     nullify(TENSOR_RHS%maxprimgab)
!!$     call mem_alloc(TENSOR_RHS%maxprimgab,1,nbatches)
!!$     call ls_sizero(TENSOR_RHS%maxprimgab,dim)
!!$  ENDIF
!!$  call LSTENSOR_mem_est(TENSOR_LHS,nmemsize)
!!$  call add_mem_to_global(nmemsize)
!!$  call LSTENSOR_mem_est(TENSOR_RHS,nmemsize)
!!$  call add_mem_to_global(nmemsize)
!!$
!!$
!!$  iatom = AOfull%Batch(iBatch)%atom
!!$  IbatchLoc = AOfull%Batch(iBatch)%batch
!!$
!!$  DO jBatch=1,nbatches
!!$     jatom = AOfull%Batch(jBatch)%atom
!!$     JbatchLoc = AOfull%Batch(jBatch)%batch
!!$
!!$     I1 = TENSOR1%INDEX(iatom,jatom,1,1)
!!$     IF(I1.NE.0)THEN
!!$        I1 = TENSOR1%INDEX(jatom,iatom,1,1)
!!$        IF(I1.EQ.0)call lsquit('GAB NOT SYM in buildbatchgab thksnjvd',-1)
!!$        maxBat1 = TENSOR1%SLSAO(I1)%maxBat
!!$        nOrbA1 = TENSOR1%SLSAO(I1)%nOrb(jbatchLoc)
!!$        nOrbB1 = TENSOR1%SLSAO(I1)%nOrb(ibatchLoc+maxBat1)
!!$        sA1 = TENSOR1%SLSAO(I1)%startLocalOrb(jbatchLoc)-1
!!$        sB1 = TENSOR1%SLSAO(I1)%startLocalOrb(ibatchLoc+maxBat1)-1
!!$        nA1 = TENSOR1%SLSAO(I1)%nLocal(1)
!!$        nB1 = TENSOR1%SLSAO(I1)%nLocal(2)
!!$
!!$        !LHS
!!$        I2_LHS = TENSOR_LHS%INDEX(jatom,1,1,1)
!!$        IF(I2_LHS.EQ.0)call lsquit('I1,I2_LHS mismatch',-1)
!!$        maxBat2 = TENSOR_LHS%SLSAO(I2_LHS)%maxBat
!!$        sA2 = TENSOR_LHS%SLSAO(I2_LHS)%startLocalOrb(jbatchLoc)-1
!!$        sB2 = TENSOR_LHS%SLSAO(I2_LHS)%startLocalOrb(1+maxBat2)-1
!!$        nA2 = TENSOR_LHS%SLSAO(I2_LHS)%nLocal(1)
!!$        nB2 = TENSOR_LHS%SLSAO(I2_LHS)%nLocal(2)
!!$        call copy_batch_selms(TENSOR_LHS%SLSAO(I2_LHS)%selms,sA2,sB2,nA2,&
!!$             & nB2,TENSOR1%SLSAO(I1)%selms,sA1,sB1,nA1,nB1,nOrbA1,nOrbB1)
!!$
!!$        IF(CS_SCREEN)TENSOR_LHS%maxgab(jBatch,1) = TENSOR1%maxgab(jBatch,iBatch)
!!$        IF(PS_SCREEN)TENSOR_LHS%maxprimgab(jBatch,1) = TENSOR1%maxprimgab(jBatch,iBatch) 
!!$
!!$        I1 = TENSOR1%INDEX(iatom,jatom,1,1)
!!$        maxBat1 = TENSOR1%SLSAO(I1)%maxBat
!!$        nOrbA1 = TENSOR1%SLSAO(I1)%nOrb(ibatchLoc)
!!$        nOrbB1 = TENSOR1%SLSAO(I1)%nOrb(jbatchLoc+maxbat1)
!!$        sA1 = TENSOR1%SLSAO(I1)%startLocalOrb(ibatchLoc)-1
!!$        sB1 = TENSOR1%SLSAO(I1)%startLocalOrb(jbatchLoc+maxbat1)-1
!!$        nA1 = TENSOR1%SLSAO(I1)%nLocal(1)
!!$        nB1 = TENSOR1%SLSAO(I1)%nLocal(2)
!!$
!!$
!!$        I2_RHS = TENSOR_RHS%INDEX(1,jatom,1,1)
!!$        IF(I2_RHS.EQ.0)call lsquit('I1,I2_RHS mismatch',-1)
!!$        maxBat2 = TENSOR_RHS%SLSAO(I2_RHS)%maxBat
!!$        sA2 = TENSOR_RHS%SLSAO(I2_RHS)%startLocalOrb(1)-1
!!$        sB2 = TENSOR_RHS%SLSAO(I2_RHS)%startLocalOrb(jbatchLoc+maxbat2)-1
!!$        nA2 = TENSOR_RHS%SLSAO(I2_RHS)%nLocal(1)
!!$        nB2 = TENSOR_RHS%SLSAO(I2_RHS)%nLocal(2)
!!$        IF(associated(TENSOR1%SLSAO(I1)%selms))THEN
!!$           call copy_batch_selms(TENSOR_RHS%SLSAO(I2_RHS)%selms,sA2,sB2,nA2,&
!!$                & nB2,TENSOR1%SLSAO(I1)%selms,sA1,sB1,nA1,nB1,nOrbA1,nOrbB1)
!!$        ELSE
!!$           nullify(TENSOR_RHS%SLSAO(I2_RHS)%selms)
!!$           TENSOR_RHS%SLSAO(I2_RHS)%nelms=0
!!$        ENDIF
!!$        IF(CS_SCREEN)TENSOR_RHS%maxgab(1,jBatch) = TENSOR1%maxgab(iBatch,jBatch)
!!$        IF(PS_SCREEN)TENSOR_RHS%maxprimgab(1,jBatch) = TENSOR1%maxprimgab(iBatch,jBatch) 
!!$     ENDIF
!!$  ENDDO
!!$  IF(CS_SCREEN)call set_lst_maxgabelms(TENSOR_LHS)
!!$  IF(PS_SCREEN)call set_lst_maxprimgabelms(TENSOR_LHS)
!!$  IF(CS_SCREEN)call set_lst_maxgabelms(TENSOR_RHS)
!!$  IF(PS_SCREEN)call set_lst_maxprimgabelms(TENSOR_RHS)
!!$end SUBROUTINE build_singleBatchGab
!!$
!!$!> \brief copy an lstensor to a new lstensor
!!$!> \author T. Kjaergaard
!!$!> \date 2010
!!$!> \param TENSOR1 the original full slstensor
!!$!> \param TENSOR2 the new batch slstensor
!!$SUBROUTINE build_SingleSingleBatchGab(AOfull,AO1,AO2,ibatch,jbatch,dim1,dim2,TENSOR1,TENSOR_RHS)
!!$  implicit none
!!$  integer,intent(in) :: dim1,dim2
!!$  type(AOITEM),intent(in) :: AOfull,AO1,AO2
!!$  TYPE(LSTENSOR),intent(in)     :: TENSOR1    ! full,full
!!$  TYPE(LSTENSOR),intent(inout)  :: TENSOR_RHS ! batch,full
!!$  integer,intent(in)   :: iBatch,jBatch
!!$  !
!!$  INTEGER    :: I1,I2,dim,maxbat1,maxbat2
!!$  integer(kind=long) :: nmemsize
!!$  TYPE(SLSAOTENSOR),pointer    :: lsao1
!!$  TYPE(SLSAOTENSOR),pointer    :: lsao2
!!$  integer(kind=short),pointer     :: elms1(:),elms2(:)
!!$  integer(kind=short) :: maxelm1
!!$  integer :: iatom,jatom,ibatchloc,jbatchloc,I2_RHS!,nbast
!!$  logical :: CS_SCREEN,PS_SCREEN
!!$  integer :: sA1,sA2,sB1,sB2,nOrbA1,nOrbA2,nOrbB1,nOrbB2,I,J,nA2,nB2,nA1,nB1
!!$  CS_SCREEN = ASSOCIATED(TENSOR1%maxgab)
!!$  PS_SCREEN = ASSOCIATED(TENSOR1%maxprimgab)
!!$!  nbast = TENSOR1%nbast(1)
!!$  call lstensor_nullify(TENSOR_RHS)
!!$  call init_ps_lstensor(TENSOR_RHS,AO1,AO2,dim1,dim2,.FALSE.,6)
!!$  call LSTENSOR_mem_est(TENSOR_RHS,nmemsize)
!!$  call remove_mem_from_global(nmemsize)
!!$  IF(CS_SCREEN)THEN
!!$     nullify(TENSOR_RHS%maxgab)
!!$     call mem_alloc(TENSOR_RHS%maxgab,1,1)
!!$     TENSOR_RHS%maxgab(1,1)=shortzero
!!$  ENDIF
!!$  IF(PS_SCREEN)THEN
!!$     nullify(TENSOR_RHS%maxprimgab)
!!$     call mem_alloc(TENSOR_RHS%maxprimgab,1,1)
!!$     TENSOR_RHS%maxprimgab(1,1)=shortzero
!!$  ENDIF
!!$  call LSTENSOR_mem_est(TENSOR_RHS,nmemsize)
!!$  call add_mem_to_global(nmemsize)
!!$  iatom = AOfull%Batch(iBatch)%atom
!!$  IbatchLoc = AOfull%Batch(iBatch)%batch
!!$  jatom = AOfull%Batch(jBatch)%atom
!!$  JbatchLoc = AOfull%Batch(jBatch)%batch
!!$  I1 = TENSOR1%INDEX(iatom,jatom,1,1)
!!$  IF(I1.NE.0)THEN
!!$     maxBat1 = TENSOR1%SLSAO(I1)%maxBat
!!$     nOrbA1 = TENSOR1%SLSAO(I1)%nOrb(ibatchLoc)
!!$     nOrbB1 = TENSOR1%SLSAO(I1)%nOrb(jbatchLoc+maxBat1)
!!$     sA1 = TENSOR1%SLSAO(I1)%startLocalOrb(ibatchLoc)-1
!!$     sB1 = TENSOR1%SLSAO(I1)%startLocalOrb(jbatchLoc+maxBat1)-1
!!$     nA1 = TENSOR1%SLSAO(I1)%nLocal(1)
!!$     nB1 = TENSOR1%SLSAO(I1)%nLocal(2)
!!$     I2_RHS = TENSOR_RHS%INDEX(1,1,1,1)
!!$     IF(I2_RHS.EQ.0)call lsquit('I1,I2_RHS mismatch',-1)
!!$     maxBat2 = TENSOR_RHS%SLSAO(I2_RHS)%maxBat
!!$     sA2 = TENSOR_RHS%SLSAO(I2_RHS)%startLocalOrb(1)-1
!!$     sB2 = TENSOR_RHS%SLSAO(I2_RHS)%startLocalOrb(1+maxBat2)-1
!!$     nA2 = TENSOR_RHS%SLSAO(I2_RHS)%nLocal(1)
!!$     nB2 = TENSOR_RHS%SLSAO(I2_RHS)%nLocal(2)
!!$     IF(associated(TENSOR1%SLSAO(I1)%selms))THEN
!!$        call copy_batch_selms(TENSOR_RHS%SLSAO(I2_RHS)%selms,sA2,sB2,nA2,nB2,TENSOR1%SLSAO(I1)%selms,sA1,sB1,nA1,nB1,nOrbA1,nOrbB1)
!!$     ELSE
!!$        NULLIFY(TENSOR_RHS%SLSAO(I2_RHS)%selms)
!!$        TENSOR_RHS%SLSAO(I2_RHS)%nelms=0
!!$     ENDIF
!!$     IF(CS_SCREEN)TENSOR_RHS%maxgab(1,1) = TENSOR1%maxgab(iBatch,jBatch) 
!!$     IF(PS_SCREEN)TENSOR_RHS%maxprimgab(1,1) = TENSOR1%maxprimgab(iBatch,jBatch) 
!!$  ENDIF
!!$  TENSOR_RHS%maxgabelm = TENSOR_RHS%maxgab(1,1)
!!$  TENSOR_RHS%maxprimgabelm = TENSOR_RHS%maxprimgab(1,1) 
!!$end SUBROUTINE build_SingleSingleBatchGab

!!$Subroutine cleanup_gabmatrix(GAB,CS_THRLOG,CS_SCREEN,PS_SCREEN,lupri)
!!$  implicit none
!!$  TYPE(LSTENSOR),intent(inout) :: GAB
!!$  integer :: lupri
!!$  integer(kind=short) :: CS_THRLOG
!!$  logical  :: cs_screen,ps_screen
!!$  !
!!$  TYPE(LSTENSOR) :: TENSOR
!!$  integer,pointer :: nbat1(:),nbat2(:)
!!$  integer(kind=long) :: nmemsize,AllocInt,AllocRealk,AllocIntS
!!$  logical, pointer :: nonscreenAB(:,:)
!!$  integer :: IATOM,JATOM,I,Ibatch,Jbatch,Ibat,Jbat,nbatJ,n1,n2,J,lstmem_index
!!$  integer :: nbatI,JbatchG,AOTbatch2,AOTbatch1,GABI2,maxbat,dim,I2
!!$  logical :: screen
!!$  TYPE(SLSAOTENSOR),pointer    :: slsao
!!$  call lstensor_nullify(TENSOR)
!!$!  print*,'GAB ',GAB%natom(1),GAB%natom(2),'nSLSAO',GAB%nSLSAO
!!$!  call lstensor_print(GAB,6)
!!$  IF(ps_screen.AND.cs_screen)THEN
!!$     call mem_alloc(nonscreenAB,GAB%natom(1),GAB%natom(2))
!!$     I = 0
!!$     DO JATOM = 1,GAB%nAtom(2)
!!$        nbatJ = GAB%nAObatch(JATOM,2)
!!$        DO IATOM = 1,GAB%nAtom(1)
!!$           screen = .TRUE. ! assume that we can screen away this contrib
!!$           I2 = GAB%INDEX(IATOM,JATOM,1,1)
!!$           IF(I2.NE.0)THEN
!!$              nbatI = GAB%nAObatch(IATOM,1)
!!$
!!$              Ibatch = GAB%SLSAO(I2)%AObatch(1)
!!$              Jbatch = GAB%SLSAO(I2)%AObatch(2)
!!$              screen = .TRUE.
!!$              loopJ: DO Jbat = 1,nbatJ
!!$                 JbatchG = Jbatch+Jbat
!!$                 DO Ibat = 1,nbatI
!!$                    IF(GAB%maxgab(Ibatch+Ibat,JbatchG).GE. CS_THRLOG)THEN
!!$                       screen = .FALSE.
!!$                       EXIT LoopJ
!!$                    ENDIF
!!$                 ENDDO
!!$              ENDDO loopJ
!!$              nonscreenAB(IATOM,JATOM) = .NOT. SCREEN
!!$              IF(nonscreenAB(IATOM,JATOM)) I = I + 1
!!$           ELSE
!!$              nonscreenAB(IATOM,JATOM) = .FALSE.
!!$           ENDIF
!!$        ENDDO
!!$     ENDDO
!!$
!!$
!!$     IF(I.NE.GAB%nSLSAO)THEN
!!$        TENSOR%nSLSAO = I
!!$        NULLIFY(TENSOR%SLSAO)
!!$        CALL MEM_ALLOC(TENSOR%SLSAO,I)
!!$        NULLIFY(TENSOR%INDEX)
!!$        CALL MEM_ALLOC(TENSOR%INDEX,GAB%nAtom(1),GAB%nAtom(2),1,1)
!!$        DO Jatom = 1,GAB%nAtom(2)
!!$           DO Iatom = 1,GAB%nAtom(1)
!!$              TENSOR%INDEX(Iatom,Jatom,1,1) = 0 !if 0 lsaotensor not call mem_allocd
!!$           ENDDO
!!$        ENDDO
!!$
!!$        NULLIFY(TENSOR%nAOBATCH)
!!$        CALL MEM_ALLOC(TENSOR%nAOBATCH,size(GAB%nAOBATCH,1),size(GAB%nAOBATCH,2))
!!$        TENSOR%nAOBATCH = GAB%nAOBATCH
!!$        !memory estimate
!!$        AllocInt = 0
!!$        AllocIntS = 0
!!$        AllocRealk = 0
!!$        DO JATOM = 1,GAB%nAtom(2)
!!$           nBatJ = GAB%nAOBATCH(Jatom,2)
!!$           DO IATOM = 1,GAB%nAtom(1)
!!$              IF(nonscreenAB(Iatom,Jatom))THEN
!!$                 GABI2 = GAB%INDEX(IATOM,JATOM,1,1)
!!$                 nBatI = GAB%nAOBATCH(Iatom,1)
!!$                 maxBat = MAX(nbatI,nbatJ)
!!$                 AllocInt = AllocInt + 4*maxBat
!!$                 DIM = GAB%SLSAO(GABI2)%nelms
!!$                 IF(associated(GAB%SLSAO(GABI2)%selms))THEN
!!$                    AllocIntS = AllocIntS + DIM
!!$                 ENDIF
!!$              ENDIF
!!$           ENDDO
!!$        ENDDO
!!$        call init_lstensorMem(AllocInt,AllocRealk,AllocInts,lstmem_index)
!!$        call zero_lstensorMem
!!$        TENSOR%lstmem_index = lstmem_index
!!$
!!$        I = 0
!!$        DO JATOM = 1,GAB%nAtom(2)
!!$           nBatJ = TENSOR%nAOBATCH(Jatom,2)
!!$           DO IATOM = 1,GAB%nAtom(1)
!!$              IF(nonscreenAB(Iatom,Jatom))THEN
!!$                 I=I+1
!!$                 TENSOR%INDEX(IATOM,JATOM,1,1) = I
!!$                 GABI2 = GAB%INDEX(IATOM,JATOM,1,1)
!!$                 slsao => TENSOR%SLSAO(I)
!!$                 slsao%ATOM(1) = Iatom
!!$                 slsao%ATOM(2) = Jatom
!!$                 nBatI = TENSOR%nAOBATCH(Iatom,1)
!!$                 slsao%AOBATCH = GAB%SLSAO(GABI2)%AOBATCH
!!$                 maxBat = MAX(nbatI,nbatJ)
!!$                 slsao%maxBat = maxBat
!!$                 NULLIFY(slsao%nOrb)
!!$                 NULLIFY(slsao%startLocalOrb)
!!$                 CALL MEM_LSTPOINTER_ALLOC(slsao%nOrb,maxBat*2)
!!$                 CALL MEM_LSTPOINTER_ALLOC(slsao%startLocalOrb,maxBat*2)
!!$                 slsao%nOrb = GAB%SLSAO(GABI2)%nOrb
!!$                 slsao%startLocalOrb = GAB%SLSAO(GABI2)%startLocalOrb
!!$                 slsao%nLocal = GAB%SLSAO(GABI2)%nLocal
!!$                 DIM = GAB%SLSAO(GABI2)%nelms
!!$                 slsao%nelms = DIM
!!$                 NULLIFY(slsao%selms)
!!$                 IF(associated(GAB%SLSAO(GABI2)%selms))THEN
!!$                    CALL MEM_LSTPOINTER_ALLOC(slsao%selms,DIM)
!!$                    DO J=1,DIM
!!$                       slsao%selms(J) = GAB%SLSAO(GABI2)%selms(J)
!!$                    ENDDO
!!$                 ENDIF
!!$              ENDIF
!!$           ENDDO
!!$        ENDDO
!!$        IF(I.NE.TENSOR%nSLSAO)call lsquit('cleanup_gabmatrix',-1)
!!$        NULLIFY(TENSOR%LSAO)
!!$        DO I=1,4
!!$           TENSOR%nAtom(I) = GAB%nAtom(I)
!!$        ENDDO
!!$        DO I=1,4
!!$           TENSOR%nbast(I) = GAB%nbast(I)
!!$        ENDDO
!!$        TENSOR%ndim5= GAB%ndim5
!!$        DO I=1,4
!!$           TENSOR%nbatches(I) = GAB%nbatches(I)
!!$        ENDDO
!!$        TENSOR%nLSAO= GAB%nLSAO
!!$        TENSOR%magGradienttensor= GAB%magGradienttensor
!!$        TENSOR%Gradienttensor= GAB%Gradienttensor
!!$        TENSOR%pChargetensor= GAB%pChargetensor
!!$        TENSOR%Econtrib= GAB%Econtrib
!!$        TENSOR%Screentensor= GAB%Screentensor
!!$        IF(ASSOCIATED(GAB%maxgab))THEN
!!$           n1 = size(GAB%maxgab,1)
!!$           n2 = size(GAB%maxgab,2)
!!$           call mem_alloc(TENSOR%maxgab,n1,n2)
!!$           DO J=1,n2
!!$              DO I=1,n1
!!$                 TENSOR%maxgab(I,J) = GAB%maxgab(I,J)
!!$              ENDDO
!!$           ENDDO
!!$        ELSE
!!$           NULLIFY(TENSOR%maxgab)
!!$        ENDIF
!!$        IF(ASSOCIATED(GAB%maxprimgab))THEN
!!$           n1 = size(GAB%maxprimgab,1)
!!$           n2 = size(GAB%maxprimgab,2)
!!$           call mem_alloc(TENSOR%maxprimgab,n1,n2)
!!$           DO J=1,n2
!!$              DO I=1,n1
!!$                 TENSOR%maxprimgab(I,J) = GAB%maxprimgab(I,J)
!!$              ENDDO
!!$           ENDDO
!!$        ELSE
!!$           NULLIFY(TENSOR%maxprimgab)
!!$        ENDIF
!!$        TENSOR%maxgabelm = GAB%maxgabelm
!!$        TENSOR%maxprimgabelm = GAB%maxprimgabelm
!!$        call lstensor_free(GAB)
!!$        call LSTENSOR_copy(TENSOR,GAB,lupri)
!!$        call LSTENSOR_mem_est(TENSOR,nmemsize)
!!$        call add_mem_to_global(nmemsize)
!!$        call lstensor_free(TENSOR)
!!$     ENDIF
!!$     call mem_dealloc(nonscreenAB)
!!$  ENDIF
!!$!  print*,'THE CLEANUP'
!!$!  call lstensor_print(GAB,6)
!!$
!!$End Subroutine Cleanup_gabmatrix
!!$
!!$Subroutine verify_gabmatrix(GAB,CS_THRLOG,lupri)
!!$  implicit none
!!$  TYPE(LSTENSOR),intent(inout) :: GAB
!!$  integer(kind=short) :: CS_THRLOG
!!$  integer :: lupri
!!$  !
!!$  TYPE(LSTENSOR) :: TENSOR
!!$  integer,pointer :: nbat1(:),nbat2(:)
!!$  integer(kind=long) :: nmemsize
!!$  integer :: IATOM,JATOM,I,Ibatch,Jbatch,Ibat,Jbat,nbatJ,n1,n2,J,Ibatch2
!!$  integer :: nbatI,JbatchG,AOTbatch2,AOTbatch1,GABI2,maxbat,dim,I2,Jbatch2
!!$  logical :: screen
!!$  TYPE(SLSAOTENSOR),pointer    :: slsao
!!$!  print*,'GAB ',GAB%natom(1),GAB%natom(2),'nSLSAO',GAB%nSLSAO
!!$!  call lstensor_print(GAB,6)
!!$  Jbatch = 0
!!$  DO JATOM = 1,GAB%nAtom(2)
!!$     nbatJ = GAB%nAObatch(JATOM,2)
!!$     Ibatch = 0
!!$     DO IATOM = 1,GAB%nAtom(1)
!!$        nbatI = GAB%nAObatch(IATOM,1)
!!$        I2 = GAB%INDEX(IATOM,JATOM,1,1)
!!$        IF(I2.EQ.0)THEN
!!$!           WRITE(lupri,*)'THIS BLOCK SHOULD BE ZERO  I=',I2,'IATOM,JATOM',IATOM
!!$           DO Jbat = 1,nbatJ
!!$              JbatchG = Jbatch+Jbat
!!$              DO Ibat = 1,nbatI
!!$                 IF(GAB%maxgab(Ibatch+Ibat,JbatchG).GE.CS_THRLOG)&
!!$                      &call lsquit('ERROR in verification',-1)
!!$              ENDDO
!!$           ENDDO
!!$           DO Jbat = 1,nbatJ
!!$              JbatchG = Jbatch+Jbat
!!$              DO Ibat = 1,nbatI
!!$                 IF(GAB%maxprimgab(Ibatch+Ibat,JbatchG).GE.CS_THRLOG)&
!!$                      &call lsquit('ERROR in verification',-1)
!!$              ENDDO
!!$           ENDDO
!!$        ENDIF
!!$        Ibatch = Ibatch + nbatI 
!!$     ENDDO
!!$     Jbatch = Jbatch + nbatJ 
!!$  ENDDO
!!$End Subroutine Verify_gabmatrix
!!$

subroutine build_atomicODScreen(AO1,AO2,nonscreenAB,natoms1,natoms2)
  implicit none
  TYPE(AOITEM)       :: AO1,AO2
  INTEGER            :: natoms1,natoms2
  logical            :: nonscreenAB(natoms1,natoms2)
  !
  INTEGER :: AOTbatch(2),Jatom,Iatom,nbat(2),batAOT1,batAOT2
  INTEGER :: Ibat,Jbat
  LOGICAL :: OVERALLscreen,screen
  
  IF(AO1%natoms.NE.natoms1)call lsquit('dim mismatch1 in build_atomicODScreen',-1)
  IF(AO2%natoms.NE.natoms2)call lsquit('dim mismatch2 in build_atomicODScreen',-1)
  AOTbatch(2)=0
  DO Jatom = 1,natoms2
     nbat(2) = AO2%ATOMICnBatch(JATOM)
     AOTbatch(1)=0
     DO Iatom = 1,natoms1
        nbat(1) = AO1%ATOMICnBatch(IATOM)
        
        OVERALLscreen = .TRUE.
        
        batAOT1 = AOTbatch(1)
        batchloop: DO Ibat = 1,nbat(1)
           batAOT1 = batAOT1+1 
           batAOT2 = AOTbatch(2)
           DO Jbat = 1,nbat(2)
              batAOT2 = batAOT2+1 
              screen = .FALSE.
              call getODscreening(AO1%BATCH(batAOT1),AO2%BATCH(batAOT2),screen)
              IF(.NOT.screen)THEN
                 !   print*,'ODscreenLHS',ODscreenLHS,'is false so we should calc'
                 !if just one element is not screen - we do not screen any
                 OVERALLscreen = .FALSE.
                 EXIT batchloop
              ENDIF
           ENDDO
        ENDDO batchloop
        nonscreenAB(Iatom,Jatom) = .NOT.OVERALLscreen
        AOTbatch(1) = AOTbatch(1) + nbat(1)
     ENDDO
     AOTbatch(2) = AOTbatch(2) + nbat(2)
  ENDDO
END subroutine BUILD_ATOMICODSCREEN

END MODULE LSTENSOR_OPERATIONSMOD
