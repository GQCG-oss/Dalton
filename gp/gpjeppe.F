!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
#include <single.h>
C
C  /* Deck istvc2 */
      SUBROUTINE ISTVC2(IVEC,IBASE,IFACT,NDIM)
C
C IVEC(I) = IBASE + IFACT * I
C
      DIMENSION IVEC(NDIM)
C
      DO 100 I = 1,NDIM
         IVEC(I) = IBASE + IFACT*I
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ibion */
      FUNCTION IBION(M,N)
C
C BIONOMIAL COEFFICIENT (M / N ) = IFAC(M)/(IFAC(M-N)*IFAC(N))
C
C Revised 2-Jun-1989 Hans Joergen Aa. Jensen to use
C    floating point arithmetic, for larger range
C
#if defined (VAR_SECSEC)
C Alliant: experiments showed problems with -DAS (890502/hjaaj).
C The reason was found to be that integer arithmetic then was used
C in the following way:
C FNFAC = PRODUCT(<1:N1>); FBION = PRODUCT(<M-N1+1:M>) / FNFAC
C To get double precision accuracy one may either use the nonstandard
C DFLOAT(K) or FK (see below) to use default type conversion.
C (Note that nonstandard FLOAT is non-generic, so FLOAT(K)
C gives only single precision accuracy!).
C
#endif
#include <implicit.h>
      N1 = MIN(N,M-N)
      FNFAC = 1.0D0
      DO 50 K = 2,N1
         FK = K
C        ... to force use of floating point multiply
         FNFAC = FNFAC * FK
   50 CONTINUE
      FBION = 1.0D0 / FNFAC
      DO 100 K = (M-N1+1), M
         FK = K
C        ... to force use of floating point multiply
         FBION = FBION * FK
  100 CONTINUE
C
#if defined (VAR_ANOTHERSOLUTION)
      FK = M - N1
      FBION = 1
      DO 100 I = 1,N1
         FI  = I
         FBION = FBION + (FBION*FK) / FI
C              = FBION * (FK+FI) / FI
  100 CONTINUE
#endif
C
      IBION = NINT(FBION)
C
      RETURN
      END
C  /* Deck ifac */
      FUNCTION IFAC(N)
C
C N!
C
#include <priunit.h>
C
      IF ( N .GE. 0 ) THEN
         IFAC = 1
         DO 100 I = 2, N
            IFAC = IFAC * I
  100    CONTINUE
      ELSE
         WRITE (LUPRI,'(/A)')
     *   ' WARNING: FACULTY OF NEGATIVE NUMBER SET TO ZERO'
         IFAC = 0
      END IF
C
      RETURN
      END
C  /* Deck isetvc */
      SUBROUTINE ISETVC(IVEC,IVALUE,NDIM)
C
      DIMENSION IVEC(NDIM)
C
      DO 100 I = 1, NDIM
         IVEC(I) = IVALUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck iwrtma */
      SUBROUTINE IWRTMA(IMAT,NROW,NCOL,MAXROW,MAXCOL)
C  Jeppe Olsen
C  Rev. 910722/921202 hjaaj
      DIMENSION IMAT(MAXROW,MAXCOL)
#include <priunit.h>
C
      IMATMX = 0
      DO 120 J = 1,NCOL
         DO 110 I = 1, NROW
            IMATMX = MAX( IMATMX, ABS(IMAT(I,J)) )
  110    CONTINUE
  120 CONTINUE
C
      IF (IMATMX .EQ. 0) THEN
         WRITE(LUPRI,'(/A,2I8)')
     &      ' -- Zero integer matrix of dim. (row, col) :',NROW,NCOL
      ELSE
         NZROW = 0
         DO 300 I = 1, NROW
            DO 210 J = 1,NCOL
               IF (IMAT(I,J) .NE. 0) GO TO 220
  210       CONTINUE
            NZROW = NZROW + 1
            GO TO 300
  220       CONTINUE
            IF (IMATMX .LT. 1 000 000) THEN
               WRITE(LUPRI,1200) I,(IMAT(I,J),J= 1,NCOL)
            ELSE
               WRITE(LUPRI,1300) I,(IMAT(I,J),J= 1,NCOL)
            END IF
  300    CONTINUE
         IF (NZROW .GT. 0) WRITE(LUPRI,1400) NZROW
      END IF
 1200 FORMAT( I6,' :',10I7,/,(8X,10I7))
 1300 FORMAT( I6,' :',5I14,/,(8X,5I14))
 1400 FORMAT(/I6,' rows with only zeroes were not printed.')
C
      RETURN
      END
C  /* Deck icopve */
      SUBROUTINE ICOPVE(IFROM,ITO,NDIM)
C
C COPY INTEGER ARRAY
C
      DIMENSION IFROM(NDIM),ITO(NDIM)
C
      DO 100 I = 1,NDIM
        ITO(I) = IFROM(I)
  100 CONTINUE
C
      RETURN
      END
C  /* Deck iscave */
      SUBROUTINE ISCAVE(IVEC,ISCALE,NDIM)
C
      DIMENSION IVEC(NDIM)
C
      DO 100 I = 1,NDIM
        IVEC(I) = ISCALE*IVEC(I)
  100 CONTINUE
C
      RETURN
      END
C  /* Deck setvec */
      SUBROUTINE SETVEC(VECTOR,VALUE,NDIM)
C
C VECTOR(*) = VALUE
C
#include <implicit.h>
      DIMENSION VECTOR(NDIM)
C
      DO 100 I = 1, NDIM
         VECTOR(I) = VALUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck copvec */
      SUBROUTINE COPVEC(VECIN,VECOUT,NDIM)
C
C 880717 - HJAaJ - written based on a qualified guess
C                  about Jeppe's original
C
#include <implicit.h>
      DIMENSION VECIN(NDIM), VECOUT(NDIM)
      DO 100 I = 1,NDIM
         VECOUT(I) = VECIN(I)
  100 CONTINUE
      RETURN
      END
C  /* Deck wrtmat */
      SUBROUTINE WRTMAT(AMATRX,NRDIM,NCDIM,NRMAX,NCMAX,ITRANS)
C
#include <implicit.h>
      DIMENSION AMATRX(NRMAX,NCMAX)
#include <priunit.h>
C
      IF ( NRDIM .EQ. 1 ) THEN
         WRITE(LUPRI,1011) (AMATRX(1,J),J=1,NCDIM)
      ELSE IF ( ITRANS.EQ. 0 ) THEN
         DO 100 I = 1,NRDIM
            WRITE(LUPRI,1010) I,(AMATRX(I,J),J=1,NCDIM)
  100    CONTINUE
      ELSE
         DO 101 I = 1, NCDIM
            WRITE(LUPRI,1010) I,(AMATRX(J,I),J=1,NRDIM)
  101    CONTINUE
      END IF
C
 1010 FORMAT(/,I6,1P,4E16.8,/,(6X,1P,4E16.8) )
 1011 FORMAT(/,(6X,1P,4E16.8) )
      RETURN
      END
C  /* Deck vecsum */
      SUBROUTINE VECSUM(C,A,B,FACA,FACB,NDIM)
C
#include <implicit.h>
C
C C(*) = FACA*A(*) + FACB*B(*)
C
      DIMENSION C(NDIM),A(NDIM),B(NDIM)
C
      DO 100 I =1,NDIM
         C(I) = FACA*A(I) + FACB*B(I)
  100 CONTINUE
C
      RETURN
      END
C  /* Deck matvcd */
      SUBROUTINE MATVCD(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
#include <implicit.h>
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL) || defined (SYS_T90)
      REAL
#else
      DOUBLE PRECISION
#endif
     *        MATRIX(MATDIM,MATDIM),VECIN(*),VECOUT(*)
C
C     VECOUT=MATRIX*VECIN FOR ITRNSP=0
C     VECOUT=MATRIX(TRANSPOSED)*VECIN FOR ITRNSP .NE. 0
C
      IF(ITRNSP.EQ.0) THEN
       DO 10 I=1,NDIM
   10   VECOUT(I)=0.0D0
C
       DO 100 J=1,NDIM
        VECINJ=VECIN(J)
        DO 90 I=1,NDIM
         VECOUT(I)=VECOUT(I)+MATRIX(I,J)*VECINJ
   90   CONTINUE
  100  CONTINUE
C
      ELSE
       DO 200 I=1,NDIM
        X=0.0D0
        DO 190 J=1,NDIM
         X=X+MATRIX(J,I)*VECIN(J)
  190   CONTINUE
        VECOUT(I)=X
  200  CONTINUE
      END IF
      RETURN
      END
C  /* Deck memadd */
      SUBROUTINE MEMADD(KBASE,KADD,KFREE,IR)
C
C VERY SUBTLE ROUTINE FOR DYNAMIC ALLOCATION
C
#include <iratdef.h>
C
      IF (KADD .GE. 0) THEN
         KBASE = KFREE
         IF ( IR .EQ. 1 ) THEN
            KFREE = KFREE + (KADD-1)/IRAT + 1
         ELSE
            KFREE = KFREE + KADD
         END IF
      ELSE IF (KADD .LT. 0) THEN
         KBASE = 999 999 999
         IF ( IR .EQ. 1 ) THEN
            KFREE = KFREE - (-KADD-1)/IRAT - 1
         ELSE
            KFREE = KFREE + KADD
         END IF
      END IF
C
      RETURN
      END
C  /* Deck msaxty */
      SUBROUTINE MSAXTY(AX,A,X,TEST,NDIM,NVEC,INDEX,NVCEFF)
C
C AX(I) = SUM(L=1,NVEC) A(L)*X(I,INDEX(L))
C
#include <implicit.h>
      DIMENSION AX(NDIM),  X(NDIM,NVEC)
      DIMENSION A(*),      INDEX(*)
C     DIMENSION A(NVCEFF), INDEX(NVCEFF)
      IF (NVCEFF .LE. 0) THEN
         DO 10 I = 1,NDIM
            AX(I) = 0.0D0
   10    CONTINUE
         RETURN
      END IF
C
      DO 100 I = 1,NDIM
         AX(I) = A(1)*X(I,INDEX(1))
  100 CONTINUE
      DO 300 L = 2,NVCEFF
         DO 200 I = 1,NDIM
            AX(I) = AX(I) + A(L)*X(I,INDEX(L))
  200    CONTINUE
  300 CONTINUE
      RETURN
      END
C  /* Deck msaxpy */
      SUBROUTINE MSAXPY(AX,A,X,TEST,NDIM,NVEC,INDEX,NVCEFF)
C
C AX(I) = AX(I) + SUM(L=1,NVEC) A(L)*X(I,INDEX(L))
C
#include <implicit.h>
      DIMENSION AX(NDIM),  X(NDIM,NVEC)
      DIMENSION A(*),      INDEX(*)
C     DIMENSION A(NVCEFF), INDEX(NVCEFF)
C
      DO 300 L = 1,NVCEFF
         DO 200 I = 1,NDIM
            AX(I) = AX(I) + A(L)*X(I,INDEX(L))
  200    CONTINUE
  300 CONTINUE
      RETURN
      END
C  /* Deck reormt */
      SUBROUTINE REORMT(AIN,AOUT,NROW,NCOL,IROW,ICOL)
C
C REORDER MATRIX AIN TO GIVE AOUT
C
C  AOUT(I,J) = AIN(IROW(I),ICOL(J) )
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION AIN(NROW,NCOL),AOUT(NROW,NCOL)
      DIMENSION IROW(NROW),ICOL(NCOL)
C
      DO 200 J = 1, NCOL
        JJ = ICOL(J)
        DO 100 I = 1, NROW
         AOUT(I,J) = AIN(IROW(I),JJ)
  100   CONTINUE
  200 CONTINUE
C
#if defined (VAR_TSTREORMT)
      NTEST = 1
      IF( NTEST .NE. 0 ) THEN
        WRITE(LUPRI,*) ' input and output matrix from REORMT '
        CALL WRTMAT(AIN,NROW,NCOL,NROW,NCOL,0)
        WRITE(LUPRI,*)
        CALL WRTMAT(AOUT,NROW,NCOL,NROW,NCOL,0 )
      END IF
#endif
C
      RETURN
      END
C  /* Deck trpmat */
      SUBROUTINE TRPMAT(XIN,NROW,NCOL,XOUT)
C
C XOUT(I,J) = XIN(J,I)
C
#include <implicit.h>
      DIMENSION XIN(NROW,NCOL),XOUT(NCOL,NROW)
C
      DO 200 IROW =1, NROW
        DO 100 ICOL = 1, NCOL
          XOUT(ICOL,IROW) = XIN(IROW,ICOL)
  100   CONTINUE
  200 CONTINUE
C
      RETURN
      END
C  /* Deck scalve */
      SUBROUTINE SCALVE(VECTOR,FACTOR,NDIM)
C
C CALCULATE SCALAR(FACTOR) TIMES VECTOR
C
#include <implicit.h>
      DIMENSION VECTOR(*)
C
      DO 100 I=1,NDIM
       VECTOR(I)=VECTOR(I)*FACTOR
  100 CONTINUE
C
      RETURN
      END
C  /* Deck prsym */
      SUBROUTINE PRSYM(A,MATDIM)
C PRINT LOWER HALF OF A SYMMETRIC MATRIX OF DIMENSION MATDIM.
C THE LOWER HALF OF THE MATRIX IS SUPPOSED TO BE IN VECTOR A.
#include <implicit.h>
#include <priunit.h>
      DIMENSION A(*)
      JSTART=1
      JSTOP=0
      DO 100 I=1,MATDIM
        JSTART=JSTART+I-1
        JSTOP=JSTOP +I
        WRITE(LUPRI,1010) I,(A(J),J=JSTART,JSTOP)
  100 CONTINUE
      RETURN
 1010 FORMAT(/I5,':',5E14.7,/,(6X,5E14.7))
      END
C  /* Deck iswpve */
      SUBROUTINE ISWPVE(IVEC1,IVEC2,NDIM)
C
C SWAP INTEGER ARRAYS IVEC1 AND IVEC2
C
      DIMENSION IVEC1(*),IVEC2(*)
C
      DO 100 I = 1, NDIM
       IBUF = IVEC1(I)
       IVEC1(I) = IVEC2(I)
       IVEC2(I) = IBUF
  100 CONTINUE
C
      RETURN
      END
C  /* Deck prsm2 */
      SUBROUTINE PRSM2(A,NDIM)
C
C PRINT LOWER TRIANGULAR MATRIX PACKED IN COLUMN WISE FASHION
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION A(*)
C
      DO 100 I=1,NDIM
        WRITE(LUPRI,1010) I,
     &  (A((J-1)*NDIM-J*(J-1)/2+I),J=1,I)
  100 CONTINUE
      RETURN
 1010 FORMAT(/I5,':',5E14.7,/,(5X,5E14.7))
      END
C  /* Deck prsbl3 */
      SUBROUTINE PRSBL3(A,NRC,LRC,NCC,LCC,IBLK)
C
C PRINT BLOCKED MATRIX STORED AS LOWER HALF OF A MATRIX
C
C JANUARY 1989 , JEPPE OLSEN
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION A(*)
      DIMENSION LRC(NRC),LCC(NCC),IBLK(NRC,NCC)
C
      DO 200 IRC = 1, NRC
        DO 100 ICC = 1, NCC
C?        WRITE(LUPRI,*) ' IRC ICC IBLK ', IRC,ICC,IBLK(IRC,ICC)
          IF( IBLK(IRC,ICC).GT.0 .AND. LRC(IRC)*LCC(ICC).NE.0) THEN
            WRITE(LUPRI,'(A,2I3/A/)')
     &      '   BLOCK...',IRC,ICC,
     &      '  =================='
            IPNTR = IBLK(IRC,ICC)
            IF(IRC .NE. ICC ) THEN
              CALL WRTMAT(A(IPNTR),LRC(IRC),
     &                    LCC(ICC),LRC(IRC),LCC(ICC),0)
            ELSE
              CALL PRSM2(A(IPNTR),LRC(IRC) )
            END IF
            WRITE(LUPRI,'()')
          END IF
  100   CONTINUE
  200 CONTINUE
C
      RETURN
      END
C  /* Deck diavc2 */
      SUBROUTINE DIAVC2(VECOUT,VECIN,DIAG,SHIFT,NDIM)
C
C VECOUT(I)=VECIN(I)/(DIAG(I)+SHIFT)
C
#include <implicit.h>
      PARAMETER (THRES=1.0D-4)
      DIMENSION VECOUT(*),VECIN(*),DIAG(*)
C
      DO 100 I=1,NDIM
      DIVIDE=DIAG(I)+SHIFT
      IF(ABS(DIVIDE).LE.THRES) DIVIDE=THRES
      VECOUT(I)=VECIN(I)/DIVIDE
  100 CONTINUE
      RETURN
      END
C  /* Deck gatvec */
      SUBROUTINE GATVEC(VECO,VECI,INDEX,NDIM)
C
C GATHER VECTOR :
C VECO(I) = VECI(INDEX(I))
C
#include <implicit.h>
      DIMENSION VECI(*),VECO(*),INDEX(*)
C
      DO 100 I = 1, NDIM
  100 VECO(I) = VECI(INDEX(I))
C
      RETURN
      END
C  /* Deck scavec */
      SUBROUTINE SCAVEC(VECO,VECI,INDEX,NDIM)
C
C SCATTER VECTOR
C VECO(INDEX(I)) = VECI(I)
C
C NOTE: no dependence assumed for VECO
C
#include <implicit.h>
      DIMENSION VECI(*),VECO(*),INDEX(*)
C
#include <nosync.h>
      DO 100 I = 1, NDIM
  100 VECO(INDEX(I)) = VECI(I)
C
      RETURN
      END
C  /* Deck xtrcdi */
      SUBROUTINE XTRCDI(AMAT,DIAG,NDIM,ISYM)
C
C EXTRACT DIAGONAL OF A MATRIX
C
C IF ISYM .LE. 0 MATRIX IS ASSUMED STORED IN COMPLETE FORM
C IF ISYM .GT. 0 MATRIX IS ASSUMED PACKED ROWWISE IN
C                LOWER TRIANGULAR FORM
C
#include <implicit.h>
      DIMENSION AMAT(*),DIAG(*)
C
      DO 100 I = 1, NDIM
        IF ( ISYM .EQ. 0 ) THEN
          II = (I-1)*NDIM + I
        ELSE
          II = I*(I+1)/2
        END IF
        DIAG(I) = AMAT(II)
  100 CONTINUE
C
      RETURN
      END
C  /* Deck imnxvc */
      SUBROUTINE IMNXVC(IVEC,NDIM,MXMN,IVAL,IPLACE)
C
C MXMN = 1 : FIND LARGEST ELEMENT IN IVEC
C MXMN = 2 : FIND SMALLEST ELEMENT IN IVEC
C
C RESULTING VALUE : IVAL
C PLACE OF RESULTING VALUE : IPLACE
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION IVEC(*)
C
      IVAL = IVEC(1)
      IPLACE = 1
      IF( MXMN .EQ. 1 ) THEN
        DO 100 I = 2, NDIM
          IF(IVEC(I) .GE. IVAL ) THEN
            IVAL = IVEC(I)
            IPLACE = I
          END IF
  100   CONTINUE
      ELSE IF ( MXMN .EQ. 2 ) THEN
        DO 200 I = 2, NDIM
          IF(IVEC(I) .LE. IVAL ) THEN
            IVAL = IVEC(I)
            IPLACE = I
          END IF
  200   CONTINUE
      END IF
C
      NTEST = 1
      IF( NTEST .NE. 0 )
     &WRITE(LUPRI,*) ' MXMN IVAL IPLACE ' ,MXMN,IVAL,IPLACE
C
      RETURN
      END
C  /* Deck gatcsf */
      SUBROUTINE GATCSF(NDET,ADET,BCSF,IORD)
C
C 890216-hjaaj Gather determinant vector in CSF order
C              from determinant vector in string order,
C              with sign changes caused by switch from
C              string order to configuration order
C
#include <implicit.h>
      DIMENSION ADET(NDET), BCSF(NDET)
      INTEGER   IORD(NDET)
      DO 100 I = 1,NDET
         BCSF(I) = ADET( ABS(IORD(I)) )
  100 CONTINUE
      DO 200 I = 1,NDET
         IF (IORD(I) .LT. 0) BCSF(I) = -BCSF(I)
  200 CONTINUE
      RETURN
      END
C  /* Deck scacsf */
      SUBROUTINE SCACSF(NDET,ADET,BCSF,IORD)
C
C 890216-hjaaj Scatter determinant vector in CSF order
C              to determinant vector in string order,
C              with sign changes caused by switch from
C              string order to configuration order
C
#include <implicit.h>
      DIMENSION ADET(NDET), BCSF(NDET)
      INTEGER   IORD(NDET)
      DO 100 I = 1,NDET
         IF (IORD(I) .LT. 0) BCSF(I) = -BCSF(I)
  100 CONTINUE
#include <nosync.h>
      DO 200 I = 1,NDET
         ADET( ABS(IORD(I)) ) = BCSF(I)
  200 CONTINUE
      RETURN
      END
C  /* Deck matml4 */
      SUBROUTINE MATML4(C,A,B,NCROW,NCCOL,NAROW,NACOL,
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
#include <implicit.h>
#include <priunit.h>
      DIMENSION A(NAROW,NACOL),B(NBROW,NBCOL)
      DIMENSION C(NCROW,NCCOL)
C
      NTEST = 0
      IF ( NTEST .NE. 0 ) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' A AND B MATRIX FROM MATML4, ITRNSP =',ITRNSP
        WRITE(LUPRI,*)
        CALL WRTMAT(A,NAROW,NACOL,NAROW,NACOL,0)
        CALL WRTMAT(B,NBROW,NBCOL,NBROW,NBCOL,0)
      END IF
C
      CALL SETVEC(C,0.0D0,NCROW*NCCOL)
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
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' C MATRIX FROM MATML4, ITRNSP =',ITRNSP
        WRITE(LUPRI,*)
        CALL WRTMAT(C,NCROW,NCCOL,NCROW,NCCOL,0)
      END IF
C
      RETURN
      END
C  /* Deck vvtov */
      SUBROUTINE VVTOV(VECIN1,VECIN2,VECUT,NDIM)
C
C VECUT(I) = VECIN1(I) * VECIN2(I)
C
#include <implicit.h>
      DIMENSION VECIN1(NDIM),VECIN2(NDIM),VECUT(NDIM)
C
      DO 100 I = 1, NDIM
        VECUT(I) = VECIN1(I) * VECIN2(I)
  100 CONTINUE
C
      RETURN
      END
C  /* Deck dgmm2 */
      SUBROUTINE DGMM2 (AOUT,AIN,DIAG,IWAY,NRDIM,NCDIM)
C
C PRODUCT OF DIAGONAL MATRIX AND MATRIX :
C
C     IWAY = 1 : AOUT(I,J) = DIAG(I)*AIN(I,J)
C     IWAY = 2 : AOUT(I,J) = DIAG(J)*AIN(I,J)
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION AIN(NRDIM,NCDIM),DIAG(*)
      DIMENSION AOUT(NRDIM,NCDIM)
C
      IF ( IWAY .EQ. 1 ) THEN
         DO 100 J = 1, NCDIM
           CALL VVTOV(AIN(1,J),DIAG(1),AOUT(1,J),NRDIM)
  100    CONTINUE
C
      ELSE IF( IWAY .EQ. 2 ) THEN
        DO 200 J = 1, NCDIM
          FACTOR = DIAG(J)
          CALL VECSUM(AOUT(1,J),AOUT(1,J),AIN(1,J),0.0D0,
     &                FACTOR,NRDIM)
  200   CONTINUE
      END IF
C
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(LUPRI,*) ' *AIN* DIAG AOUT  FROM DGMTMT '
        CALL WRTMAT(AIN ,NRDIM,NCDIM,NRDIM,NCDIM,0)
        WRITE(LUPRI,*) ' AIN *DIAG* AOUT  FROM DGMTMT '
        IF(IWAY.EQ.1) THEN
          CALL WRTMAT(DIAG,1   ,NRDIM,1,NRDIM,0)
        ELSE
          CALL WRTMAT(DIAG,1   ,NCDIM,1,NCDIM,0)
        END IF
        WRITE(LUPRI,*) ' AIN DIAG *AOUT*  FROM DGMTMT '
        CALL WRTMAT(AOUT,NRDIM,NCDIM,NRDIM,NCDIM,0)
      END IF
C
      RETURN
      END
C  /* Deck imxvec */
      FUNCTION IMXVEC(IVEC,NELMNT)
c
c Largest element in Integer vector IVEC
c
#include <implicit.h>
#include <priunit.h>
      DIMENSION IVEC(*)
c
      IF(NELMNT.LE.0 ) THEN
        WRITE(LUPRI,*) ' >> IMXVEC << in problems '
        WRITE(LUPRI,*) ' Largest element of a vector with zero elements'
        WRITE(LUPRI,*) ' is hardly defined so STOP '
        CALL QUIT('IMXVEC ')
      END IF
c
      IMX = IVEC(1)
      DO 100 IELMNT = 2, NELMNT
        IMX = MAX(IMX,IVEC(IELMNT))
  100 CONTINUE
c
      IMXVEC = IMX
c
      RETURN
      END
C  /* Deck inprod */
      FUNCTION INPROD(A,B,NDIM)
C
C     PURPOSE: CALCULATE SCALAR PRODUCT BETWEEN TO VECTORS A,B
C
#include <implicit.h>
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL) || defined (SYS_T90)
      REAL
#else
      DOUBLE PRECISION
#endif
     *   INPROD
C
      DIMENSION A(*),B(*)
C
      INPROD=0.0D0
      DO 100 I=1,NDIM
       INPROD=INPROD+A(I)*B(I)
  100 CONTINUE
C
      RETURN
      END
C  /* Deck minprd */
      SUBROUTINE MINPRD(VU,A,VI,IP,NPROD,NROW)
C
C VU(I) = SUM(J) A(J,IP(I))*VI(J)
C
#include <implicit.h>
      DIMENSION VU(*),A(NROW,*),VI(*),IP(*)
      PARAMETER (D0 = 0.0D0)
C
      DO 100 I = 1, NPROD
         VU(I) = D0
#include <prefer_vector.h>
         DO 50 J = 1,NROW
            VU(I) = VU(I) + A(J,IP(I))*VI(J)
   50    CONTINUE
  100 CONTINUE
C
      RETURN
      END
