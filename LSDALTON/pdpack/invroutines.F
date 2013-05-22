* Routines for inversion
* From Jeppe - late november 03
      SUBROUTINE INVERT_BY_DIAG2(A,B,SCR,VEC,NDIM)
*
* Invert a symmetric  - hopefully nonsingular - matrix A 
* by diagonalization. 
*
* Jeppe Olsen, Oct 97 to check INVMAT
*              March 00 : Scale initial matrix to obtain unit diagonal
*
*. The input matrix is delivered in A in complete form, and 
*. the output matrix is returned in A, also in complete form
      IMPLICIT REAL*8(A-H,O-Z)
*. Input and output matrix
      DIMENSION A(NDIM*NDIM)       
*. Scratch matrices and vector
      DIMENSION B(NDIM*NDIM),SCR(NDIM*NDIM),VEC(NDIM)
*
      NTEST = 00
*. Reform a to symmetric packed form
      CALL LS_TRIPAK(A,SCR,1,NDIM,NDIM)
*. Extract diagonal
      CALL LS_COPDIA(SCR,VEC,NDIM,1)
*
*.scale 
*
      DO I = 1, NDIM
*. Scaling vector
        IF(VEC(I).EQ.0.0D0) THEN
          VEC(I) = 1.0D0
        ELSE
          VEC(I) = 1.0D0/SQRT(ABS(VEC(I)))
        END IF
      END DO
*. Scale matrix
      IJ = 0
      DO I = 1, NDIM
        DO J = 1, I
          IJ = IJ + 1
          SCR(IJ) = SCR(IJ)*VEC(I)*VEC(J)
        END DO
      END DO
C     DO I = 1, NDIM
C       VEC(I) = 1.0D0/VEC(I)
C     END DO
*. Diagonalize
      CALL LS_EIGEN(SCR,B,NDIM,0,1)
*. Scale eigenvectors
      DO IVEC = 1, NDIM 
        IOFF = 1 + (IVEC-1)*NDIM
        CALL LS_VVTOV(B(IOFF),VEC,B(IOFF),NDIM)
      END DO
*. 
      CALL LS_COPDIA(SCR,VEC,NDIM,1)
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Eigenvalues of scaled matrix : '
        CALL LS_WRTMAT(VEC,NDIM,1,NDIM,1)
      END IF
*. Invert diagonal elements 
      DO I = 1, NDIM
       IF(ABS(VEC(I)).GT.1.0D-15) THEN
         VEC(I) = 1.0D0/VEC(I)
       ELSE
         VEC(I) = 0.0D0
         WRITE(6,*) ' Singular mode activated '
       END IF
      END DO
*. and obtain inverse matrix by transformation
C     XDIAXT(XDX,X,DIA,NDIM,SCR)
      CALL XDIAXT(A,B,VEC,NDIM,SCR)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Inverse matrix from INVERSE_BY_DIAG'
        CALL LS_WRTMAT(A,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END 
      SUBROUTINE XDIAXT(XDX,X,DIA,NDIM,SCR)
*
* Obtain XDX = X * DIA * X(Transposed)
* where DIA is an diagonal matrix stored as a vector
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XDX(NDIM,NDIM)
      DIMENSION X(NDIM,NDIM),DIA(NDIM)
      DIMENSION SCR(NDIM,NDIM)
*
* DIA * X(transposed)
      DO 100 I=1,NDIM
        CALL LS_COPVEC(X(1,I),SCR(1,I),NDIM)
        CALL LS_SCALVE(SCR(1,I),DIA(I),NDIM)
  100 CONTINUE
* X * DIA * X(Transposed)
      CALL LS_MATML4(XDX,X,SCR,NDIM,NDIM,NDIM,NDIM,NDIM,NDIM,2)
*
      RETURN
      END

      SUBROUTINE LS_COPDIA(A,VEC,NDIM,IPACK)
*
* Copy diagonal of matrix A into vector VEC
*
*   IPACK = 0 : Full matrix
*   IPACK = 1 : Lower triangular matrix
*
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
      DIMENSION A(*),VEC(*)
*
      IF(IPACK .EQ. 0 ) THEN
        DO 100 I = 1,NDIM
          VEC(I) = A((I-1)*NDIM+I)
  100   CONTINUE
      ELSE
        DO 200 I = 1, NDIM
          VEC(I) = A(I*(I+1)/2)
  200   CONTINUE
      END IF
*
      RETURN
      END
C  /* Deck vvtov */
      SUBROUTINE LS_VVTOV(VECIN1,VECIN2,VECUT,NDIM)
C
C VECUT(I) = VECIN1(I) * VECIN2(I)
C
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
      DIMENSION VECIN1(NDIM),VECIN2(NDIM),VECUT(NDIM)
C
      DO 100 I = 1, NDIM
        VECUT(I) = VECIN1(I) * VECIN2(I)
  100 CONTINUE
C
      RETURN
      END
C  /* Deck scalve */
      SUBROUTINE LS_SCALVE(VECTOR,FACTOR,NDIM)
C
C CALCULATE SCALAR(FACTOR) TIMES VECTOR
C
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
      DIMENSION VECTOR(*)
C
      DO 100 I=1,NDIM
       VECTOR(I)=VECTOR(I)*FACTOR
  100 CONTINUE
C
      RETURN
      END
C  /* Deck matml4 */
      SUBROUTINE LS_MATML4(C,A,B,NCROW,NCCOL,NAROW,NACOL,
     &                  NBROW,NBCOL,ITRNSP )
C
C MULTIPLY A AND B TO GIVE C
C
C     C = A * B             FOR ITRNSP = 0
C
C     C = A(TRANSPOSED) * B FOR ITRNSP = 1
C
C     C = A * B(TRANSPOSED) FOR ITRNSP = 2
C
C... JEPPE OLSEN, LAST REVISION JULY 24 1987
C
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
      DIMENSION A(NAROW,NACOL),B(NBROW,NBCOL)
      DIMENSION C(NCROW,NCCOL)
C
      NTEST = 0
      IF ( NTEST .NE. 0 ) THEN
        WRITE(*,*)
        WRITE(*,*) ' A AND B MATRIX FROM LS_MATML4, ITRNSP =',ITRNSP
        WRITE(*,*)
        CALL LS_WRTMAT(A,NAROW,NACOL,NAROW,NACOL,0)
        CALL LS_WRTMAT(B,NBROW,NBCOL,NBROW,NBCOL,0)
      END IF
C
      CALL LS_SETVEC(C,0.0D0,NCROW*NCCOL)
C
      IF( ITRNSP .NE. 0 ) GOTO 001
        DO 50 J = 1,NCCOL
          DO 40 K = 1,NBROW
            BKJ = B(K,J)
            DO 30 I = 1, NCROW
              C(I,J) = C(I,J) + A(I,K)*BKJ
  30        CONTINUE
  40      CONTINUE
  50    CONTINUE
C
C
  001 CONTINUE
C
      IF ( ITRNSP .NE. 1 ) GOTO 101
C... C = A(T) * B
         DO 150 J = 1, NCCOL
           DO 140 K = 1, NBROW
             BKJ = B(K,J)
             DO 130 I = 1, NCROW
               C(I,J) = C(I,J) + A(K,I)*BKJ
  130        CONTINUE
  140      CONTINUE
  150    CONTINUE
C
  101 CONTINUE
C
      IF ( ITRNSP .NE. 2 ) GOTO 201
C... C = A*B(T)
        DO 250 J = 1,NCCOL
          DO 240 K = 1,NBCOL
            BJK = B(J,K)
            DO 230 I = 1, NCROW
              C(I,J) = C(I,J) + A(I,K)*BJK
 230        CONTINUE
 240      CONTINUE
 250    CONTINUE
C
C
  201 CONTINUE
C
      IF ( NTEST .NE. 0 ) THEN
        WRITE(*,*)
        WRITE(*,*) ' C MATRIX FROM LS_MATML4, ITRNSP =',ITRNSP
        WRITE(*,*)
        CALL LS_WRTMAT(C,NCROW,NCCOL,NCROW,NCCOL,0)
      END IF
C
      RETURN
      END
C  /* Deck copvec */
      SUBROUTINE LS_COPVEC(VECIN,VECOUT,NDIM)
C
C 880717 - HJAaJ - written based on a qualified guess
C                  about Jeppe's original
C
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
      DIMENSION VECIN(NDIM), VECOUT(NDIM)
      DO 100 I = 1,NDIM
         VECOUT(I) = VECIN(I)
  100 CONTINUE
      RETURN
      END
C  /* Deck setvec */
      SUBROUTINE LS_SETVEC(VECTOR,VALUE,NDIM)
C
C VECTOR(*) = VALUE
C
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
      DIMENSION VECTOR(NDIM)
C
      DO 100 I = 1, NDIM
         VECTOR(I) = VALUE
  100 CONTINUE
C
      RETURN
      END


