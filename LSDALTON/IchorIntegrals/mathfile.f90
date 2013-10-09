module math
  INTEGER, PARAMETER :: realk = 8
  contains
!> \brief calculate spherical to cartesian matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param L angular moment
!> \param SCMAT the transformation matrix
!> \param nLM the number of spherical components
!> \param nXYZ the number of cartesian components
!> \param lupri the logical unit number for the output file
!> \param iprint the printlevel, determining how much output should be generated
SUBROUTINE Sph_to_Cart_matrix(L,SCMAT,nLM,nXYZ,LUPRI,IPRINT)
implicit none
REAL(REALK) :: SCMAT(nLM,nXYZ)
INTEGER     :: L,nLM,nXYZ,LUPRI,IPRINT
!
Real(realk), parameter :: DM1 = -1.0E0_realk, D1 = 1.0E0_realk, D2 = 2.0E0_realk
INTEGER     :: M1,MADR,MABS,V0
REAL(realk) :: FACNRM,FAC3,FACTOR
INTEGER     :: T,U,V,A,B,C,P,Q,R,X,Y,Z,M
INTEGER     :: iLM,iXYZ

IF(L .EQ. 0)THEN ! S
 SCMAT(1,1)=1E0_realk
ELSEIF(L .EQ. 1)THEN ! P
 SCMAT=0E0_realk
 SCMAT(1,1)=1E0_realk
 SCMAT(2,2)=1E0_realk
 SCMAT(3,3)=1E0_realk
ELSEIF(L .GT. 1)THEN
 SCMAT=0E0_realk
 DO M1 = 0, 2*L 
  M = M1 - L
  IF (L.EQ. 1) THEN
    IF (M .EQ. -1) MADR =  0  
    IF (M .EQ.  0) MADR =  1 
    IF (M .EQ.  1) MADR = -1 
  ELSE
    MADR = M
  END IF
  MABS = ABS(M)
  V0 = 0
  IF (M .LT. 0) V0 = 1 
  FACNRM = D1
  IF (M .NE. 0) FACNRM = SQRT(D2*FACULT(LUPRI,L+MABS)*&
                         &FACULT(LUPRI,L-MABS))/(FACULT(LUPRI,L)*(D2**MABS))
  FACNRM = FACNRM*DM1**((0-MOD(0,2))/2)*D2**(-0)
  FACNRM = FACNRM/SQRT(FACUL2(LUPRI,2*L-1))
  DO T = 0, L - MABS, 2
   DO U = 0, T, 2
    DO V = V0, MABS, 2
          !        almost 6.4.48 in the book
     FAC3 = FACNRM*BINOM(LUPRI,L,T/2)*BINOM(LUPRI,L-T/2,MABS+T/2)&
          &                    *BINOM(LUPRI,T/2,U/2)*BINOM(LUPRI,MABS,V)
     DO A = 0, MIN(0,T+MABS-U-V) 
      DO B = 0, MIN(0,U+V)
       DO C = 0, MIN(0,L-T-MABS)
             !           6.4.47 in the book
        DO P = 0, - A, 2
         DO Q = 0, - B, 2
          DO R = 0, - C, 2
           FACTOR = DM1**(A+B+C+(T+V-V0-P-Q-R)/2)*&
                &   D2**(-A-B-C-P-Q-R-T)*FAC3
           X = T+MABS-U-V-2*A-P
           Y = U+V-2*B-Q
           Z = L-T-MABS-2*C-R
           iLM = 1 + L + MADR
           iXYZ = NCRT(X,Y,Z)
           SCMAT(iLM,iXYZ) = SCMAT(iLM,iXYZ) + FACTOR 
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

END SUBROUTINE Sph_to_Cart_matrix


SUBROUTINE Buildsphericaltransformation(Spherical,Spher1,Spher2,lm1,lm2,&
     &                                      ijk1,ijk2)
implicit none
Real(realk) :: Spherical(ijk1,ijk2,lm1,lm2)
Real(realk) :: Spher1(lm1,ijk1)
Real(realk) :: Spher2(lm2,ijk2)
Integer     :: ijk1,ijk2,lm1,lm2
!
Integer     :: indijk1,indijk2,indlm1,indlm2
DO indijk2=1,ijk2
   DO indijk1=1,ijk1
      DO indlm2=1,lm2
         DO indlm1=1,lm1
            Spherical(indijk1,indijk2,indlm1,indlm2) = &
                 & Spher1(indlm1,indijk1)*Spher2(indlm2,indijk2)
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE Buildsphericaltransformation


FUNCTION FACULT(LUPRI,N)
IMPLICIT NONE
real(realk), PARAMETER :: D1=1E0_realk
integer        :: N,I,LUPRI
real(realk)    :: FACULT
IF (N .LT. 0) THEN
   WRITE (LUPRI,'(/,A,I10,/A)')&
   &         ' Argument less than zero in FACULT:',N,&
   &         ' Program cannot continue.'
   STOP'Illegal argument in FACULT'
ELSE
   FACULT = D1
   DO I = 1, N
      FACULT = FACULT*I
   ENDDO
END IF
END FUNCTION FACULT

FUNCTION FACUL2(LUPRI,N)
IMPLICIT NONE
real(realk), PARAMETER :: D1=1E0_realk
real(realk)    :: FACUL2
integer :: N,I,LUPRI
IF (N .LT. 0) THEN
   FACUL2 = (N + 2)*1.0E0_realk
   DO I = N + 4, 1, 2
      FACUL2 = FACUL2*I
   END DO
   IF (FACUL2 .EQ. 0E0_realk) THEN
      WRITE (LUPRI,'(/,A,I10,/A)')&
      &            ' Double factorial undefined for ',N,&
      &            ' Program cannot continue.'
      STOP 'Illegal argument in FACUL2'
   ELSE
      FACUL2 = D1/FACUL2
   END IF
ELSE IF (N.EQ. 0) THEN
   FACUL2 = D1
ELSE ! N > 0
   FACUL2 = N*1.0E0_realk
   DO I = N - 2, 1, -2
      FACUL2 = FACUL2*I
   END DO
END IF
END FUNCTION FACUL2

FUNCTION BINOM(LUPRI,I,J)
real(realk), PARAMETER  :: D1=1E0_realk
INTEGER :: I,J,LUPRI
real(realk)    :: BINOM

IF (I .LT. J) THEN
   WRITE (LUPRI,'(/,A,2I5,/A)')&
   &         ' Second argument larger than first argument in BINOM:',&
   &         I,J,' Program cannot continue.'
   STOP 'Illegal arguments in BINOM'
ELSE
   BINOM = FACULT(LUPRI,I)/(FACULT(LUPRI,I-J)*FACULT(LUPRI,J))
END IF
END FUNCTION BINOM

FUNCTION NCRT(I,J,K)
IMPLICIT NONE
INTEGER  :: I,J,K,NCRT
NCRT = 1 + J + 2*K + (J + K)*(J + K - 1)/2
END FUNCTION NCRT

subroutine output(AMATRX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,NCTL,LUPRI)
implicit none 
      !> ROW NUMBER AT WHICH OUTPUT IS TO BEGIN
      integer,intent(in) :: ROWLOW
      !> ROW NUMBER AT WHICH OUTPUT IS TO END
      integer,intent(in) :: ROWHI
      !> COLUMN NUMBER AT WHICH OUTPUT IS TO BEGIN
      integer,intent(in) :: COLLOW
      !> COLUMN NUMBER AT WHICH OUTPUT IS TO END
      integer,intent(in) :: COLHI
      !> ROW DIMENSION OF AMATRX(',')
      integer,intent(in) :: ROWDIM
      !> COLUMN DIMENSION OF AMATRX(',')
      integer,intent(in) :: COLDIM
      !> CARRIAGE CONTROL FLAG; (1 FOR SINGLE, 2 FOR DOUBLE, 3 FOR TRIPLE SPACE)
      integer,intent(in) :: NCTL
      !> Logical unit number for output file
      integer,intent(in) :: lupri
      !> MATRIX TO BE OUTPUT
      real(realk), intent(in) :: AMATRX(ROWDIM,COLDIM)
      integer :: BEGIN, KCOL, j, i, k, mctl, last
      CHARACTER(len=1) :: ASA(3), BLANK, CTL
      CHARACTER :: PFMT*20, COLUMN*8
      real(realk) :: amax
      real(realk), PARAMETER :: ZERO=0.0E0_realk
      integer, parameter     :: KCOLP=4, KCOLN=6
      real(realk), PARAMETER :: FFMIN=1E-3_realk, FFMAX = 1E3_realk
      DATA COLUMN/'Column  '/, BLANK/' '/, ASA/' ', '0', '-'/
 
      IF (ROWHI.LT.ROWLOW) RETURN
      IF (COLHI.LT.COLLOW) RETURN
 
      AMAX = ZERO
      DO J = COLLOW,COLHI
        DO I = ROWLOW,ROWHI
            AMAX = MAX( AMAX, ABS(AMATRX(I,J)) )
        enddo
      enddo
      IF (AMAX == ZERO) THEN
         WRITE (LUPRI,'(/T6,A)') 'Zero matrix.'
         RETURN
      END IF
      IF (FFMIN <= AMAX .AND. AMAX <= FFMAX) THEN
!        use F output format
         PFMT = '(A1,I7,2X,8F15.8)'
      ELSE
!        use 1PD output format
         PFMT = '(A1,I7,2X,1P,8D15.6)'
      END IF
!
      IF (NCTL .LT. 0) THEN
         KCOL = KCOLN
      ELSE
         KCOL = KCOLP
      END IF
      MCTL = ABS(NCTL)
      IF ((MCTL.LE. 3).AND.(MCTL.GT. 0)) THEN
         CTL = ASA(MCTL)
      ELSE
         CTL = BLANK
      END IF
!
      LAST = MIN(COLHI,COLLOW+KCOL-1)
      DO 2 BEGIN = COLLOW,COLHI,KCOL
         WRITE (LUPRI,1000) (COLUMN,I,I = BEGIN,LAST)
         DO 1 K = ROWLOW,ROWHI
            DO 4 I = BEGIN,LAST
               IF (AMATRX(K,I).NE.ZERO) GO TO 5
    4       CONTINUE
         GO TO 1
    5       WRITE (LUPRI,PFMT) CTL,K,(AMATRX(K,I), I = BEGIN,LAST)
    1    CONTINUE
    2 LAST = MIN(LAST+KCOL,COLHI)
 1000 FORMAT (/12X,6(3X,A6,I4,2X),(3X,A6,I4))
      END subroutine output

end module math
