MODULE ls_util
  use molecule_typetype
  use molecule_type
  use precision  
Integer,save :: LSIUNIT = 6

CONTAINS
SUBROUTINE LSHEADER(LUPRI,HEAD)
IMPLICIT NONE
  CHARACTER HEAD*(*)
  INTEGER :: LUPRI,LENGTH,INDENT,I
  
  LENGTH = LEN(HEAD)
  INDENT = (72 - LENGTH)/2 + 1
  
  WRITE (LUPRI, '(//,80A)') (' ',I=1,INDENT), HEAD
  WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), ('-',I=1,LENGTH)
  WRITE (LUPRI, '()')
  
END SUBROUTINE LSHEADER

SUBROUTINE LS_PRINT_GRADIENT(lupri,molecule,GRDMOL,natoms,TEXT)
  implicit none
  integer :: lupri,natoms
  type(moleculeinfo) :: molecule
  real(realk)       :: GRDMOL(3,natoms)
  CHARACTER*(*)     :: TEXT
  CHARACTER(len=15) :: PRINTTEXT
  CHARACTER(len=23) :: PRINTTEXT2
  CHARACTER(len=40) :: PRINTTEXT3
  CHARACTER(len=33) :: PRINTTEXT4
  INTEGER :: IUNIT,length,ioff,iatom,j
  IF (LUPRI.EQ.-1) THEN
     IUNIT = LSIUNIT
  ELSE
     IUNIT = LUPRI
  ENDIF

length = LEN(TEXT)
IF(length .GT. 15) CALL LSQUIT('TEXTLENGTH PROVIDED TO LS_PRINT_GRADIENT IS LIMITED TO 15',lupri)
IF (TEXT(1:5) .EQ. 'TOTAL') THEN
   PRINTTEXT2 = 'Molecular gradient (au)'
  CALL LSHEADER(iunit,PRINTTEXT2)
ELSEIF (TEXT .EQ. 'NUMGR') THEN
   PRINTTEXT4 = 'Numerical Molecular gradient (au)'
  CALL LSHEADER(iunit,PRINTTEXT4)
ELSE
   PRINTTEXT3 = TEXT(1:length)//' contribution to gradient'
   length = length+25
  CALL LSHEADER(iunit,PRINTTEXT3(1:length))
ENDIF

DO IATOM = 1, natoms
   WRITE (IUNIT,'(1X,A6,F17.10,2F24.10)') &
        &molecule%ATOM(IATOM)%Name, (GRDMOL(J,IATOM),J=1,3)
END DO

end SUBROUTINE LS_PRINT_GRADIENT

!> \brief Compute the 'almost' root-mean square of A-B, diff = RMS(A-B), where A,B are matrices
!> \author P. Merlot
!> \date 12-05-2010
!> \param diff The 'almost' RMS norm of the difference between A and B (almost, since divided by sqrt(nrow*ncol) and not by nrow*ncol) 
SUBROUTINE rms_Diff(A, B, nrow,ncol,diff)
  implicit none
  INTEGER, intent(IN)      :: nrow,ncol
  REAL(realk), intent(IN)  :: A(nrow,ncol), B(nrow,ncol)
  REAL(realk), intent(OUT) :: diff
!
  REAL(realk)              :: tempRms(nrow,ncol)    
  REAL(realk),pointer      :: WORK
  REAL(realk), external    :: dlange

  tempRms = 0.0E0_realk
  tempRms = abs(A)-abs(B)
  call lsquit('please fix the subroutine rms_DIFF: work not allocated',-1)
!  diff = dlange('F',nrow,ncol,tempRms,1,WORK) * sqrt(1E0_realk/(nrow*ncol))
END SUBROUTINE rms_Diff

!> Copy of dgemm in pdpack/gp_blas3.F modified for fortran 90,
!> this is to have a dgemm version which is always threadsafe,
!> independent of the linking.
!> Calculates: C = alpha*op(A)*op(B) + beta*C
!> See pdpack/gp_blas3.F for further details.
!> \author F90'ed by Kasper Kristensen
SUBROUTINE DGEMM_TS(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
  implicit none
  !> Scaling coefficients
  real(realk),intent(in) :: ALPHA,BETA
  !> Dimensions (see pdpack/gp_blas3.F)
  INTEGER,intent(in) :: K,LDA,LDB,LDC,M,N
  !> Transposition ('N' or 'T')
  CHARACTER(len=1),intent(in) :: TRANSA,TRANSB
  CHARACTER(len=1) :: TRANSAtmp,TRANSBtmp
  !> Input matrices
  real(realk),intent(in) :: A(lda,*), B(ldb,*)
  !> Output matrix
  real(realk),intent(inout) :: C(ldc,*)
  real(realk) :: TEMP,zero,one
  INTEGER :: I,INFO,J,L,NCOLA,NROWA,NROWB
  LOGICAL :: NOTA,NOTB

  ! Ensure capital letters for more transparent code
  TransAtmp = transA
  TransBtmp = transB
  if(TransAtmp=='n') TransAtmp='N'
  if(TransAtmp=='t') TransAtmp='T'
  if(TransBtmp=='n') TransBtmp='N'
  if(TransBtmp=='t') TransBtmp='T'
  zero = 0.0E0_realk
  one =  1.0E0_realk


  !     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
  !     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
  !     and  columns of  A  and the  number of  rows  of  B  respectively.
  NOTA = (TransAtmp=='N')
  NOTB = (TransBtmp=='N')
  IF (NOTA) THEN
     NROWA = M
     NCOLA = K
  ELSE
     NROWA = K
     NCOLA = M
  END IF
  IF (NOTB) THEN
     NROWB = K
  ELSE
     NROWB = N
  END IF

  !     Test the input parameters.
  INFO = 0
  IF ((.NOT.NOTA) .AND. (.NOT.TRANSATMP=='C') .AND. (.NOT.TRANSATMP=='T')) THEN
     INFO = 1
  ELSE IF ((.NOT.NOTB) .AND. (.NOT.TRANSBTMP=='C') .AND. (.NOT.TRANSBTMP=='T')) THEN
     INFO = 2
  ELSE IF (M.LT.0) THEN
     INFO = 3
  ELSE IF (N.LT.0) THEN
     INFO = 4
  ELSE IF (K.LT.0) THEN
     INFO = 5
  ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
     INFO = 8
  ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
     INFO = 10
  ELSE IF (LDC.LT.MAX(1,M)) THEN
     INFO = 13
  END IF
  IF (INFO.NE.0) THEN
     print *, 'INFO = ', INFO
     call lsquit('DGEMM_TS: Something wrong with input!',-1)
!     stop 'DGEMM_TS: Something wrong with input!'
  END IF

  !     Quick return if possible.
  IF ((M.EQ.0) .OR. (N.EQ.0) .OR.  &
       &  (((ALPHA.EQ.zero).OR. (K.EQ.0)).AND. (BETA.EQ.one))) RETURN

  IF (ALPHA.EQ.zero) THEN
     IF (BETA.EQ.zero) THEN
        DO J = 1,N
           DO I = 1,M
              C(I,J) = ZERO
           end DO
        end DO
     ELSE
        DO J = 1,N
           DO I = 1,M
              C(I,J) = BETA*C(I,J)
           end DO
        end DO
     END IF
     RETURN
  END IF


  ! *****************************
  !     Start the operations
  ! *****************************
  IF (NOTB) THEN
     IF (NOTA) THEN

        !           Form  C := alpha*A*B + beta*C.

        DO J = 1,N
           IF (BETA.EQ.ZERO) THEN
              DO I = 1,M
                 C(I,J) = ZERO
              end DO
           ELSE IF (BETA.NE.ONE) THEN
              DO I = 1,M
                 C(I,J) = BETA*C(I,J)
              end DO
           END IF
           DO L = 1,K
              IF (B(L,J).NE.ZERO) THEN
                 TEMP = ALPHA*B(L,J)
                 DO I = 1,M
                    C(I,J) = C(I,J) + TEMP*A(I,L)
                 end DO
              END IF
           end DO
        end DO
     ELSE

        !           Form  C := alpha*A**T*B + beta*C

        DO J = 1,N
           DO I = 1,M
              TEMP = ZERO
              DO L = 1,K
                 TEMP = TEMP + A(L,I)*B(L,J)
              end DO
              IF (BETA.EQ.ZERO) THEN
                 C(I,J) = ALPHA*TEMP
              ELSE
                 C(I,J) = ALPHA*TEMP + BETA*C(I,J)
              END IF
           end DO
        end DO
     END IF
  ELSE

     IF (NOTA) THEN

        !           Form  C := alpha*A*B**T + beta*C

        DO J = 1,N
           IF (BETA.EQ.ZERO) THEN
              DO I = 1,M
                 C(I,J) = ZERO
              end DO
           ELSE IF (BETA.NE.ONE) THEN
              DO I = 1,M
                 C(I,J) = BETA*C(I,J)
              end DO
           END IF
           DO L = 1,K
              IF (B(J,L).NE.ZERO) THEN
                 TEMP = ALPHA*B(J,L)
                 DO I = 1,M
                    C(I,J) = C(I,J) + TEMP*A(I,L)
                 end DO
              END IF
           end DO
        end DO
     ELSE

        !           Form  C := alpha*A**T*B**T + beta*C

        DO J = 1,N
           DO I = 1,M
              TEMP = ZERO
              DO L = 1,K
                 TEMP = TEMP + A(L,I)*B(J,L)
              end DO
              IF (BETA.EQ.ZERO) THEN
                 C(I,J) = ALPHA*TEMP
              ELSE
                 C(I,J) = ALPHA*TEMP + BETA*C(I,J)
              END IF
           end DO
        end DO
     END IF
  END IF


end SUBROUTINE DGEMM_TS

!> \brief Change all lower case letters in a string to capital letters
!> \author Kasper Kristensen
!> \date April 2013
subroutine capitalize_string(mystring)
  implicit none
  character(*) :: mystring
  integer :: length,gap,i

  ! Gap between lower case "a" and capital "A"
  gap=ICHAR('a')-ICHAR('A')

  length=len(mystring)

  if(length>0) then  ! only do something if string is not empty

     do i=1,length
        ! consider characters between "a" and "z"
        if(mystring(i:i) .le. 'z' .and. mystring(i:i) .ge. 'a') then
           mystring(i:i)=CHAR(ICHAR(mystring(i:i))-gap)
        end if
     end do

  end if

end subroutine capitalize_string

END MODULE LS_UTIL
