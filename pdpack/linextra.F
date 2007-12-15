C
C
C...   Dalton, Release 2.0
C...
C...   These routines are in the public domain and can be
C...   used freely in other programs.
C...
C
C FILE : pdpack/linextra.F
C
C Extra matrix and vector manipulation routines and driver routines,
C which not (always) are included in blas and linpack.
C
C Most of them are written by Hans Jørgen Aa. Jensen in the late 1980's.
C These routines may be freely used by anyone.
C
C 1) DNORM2
C 2) DAPTGE, etc. for changing matrix packing
C                 (AP: Antisymmetric Packed, SP, Symmetric Packed, etc.
C                 (ex: DAPTGE: Antisymmetric Packed To GEneral matrix)
C 3) DSUM, IDAMIN, IDMAX, IDMIN (supplement DASUM and IDAMAX)
C 4) "iblas": ICOPY, ISWAP, ...
C 5) DSPSOL, DGESOL, DGEINV : driver rouutines calling LINPACK
C
C
#include <single.h>
C  /* Deck dnorm2 */
#if !defined (VAR_ESSL) && !defined (VAR_DXML)
C DXML library on DEC-ALPHA
      FUNCTION DNORM2(N,DX,INCX)
C
C     Forms the two-norm of a vector.
C 19-Sep-1988 -- hjaaj -- based on DNRM2 from LINPACK
C     This version does not use extended precision for intermediates
C     as the LINPACK version does.
C     Equivalent to DNORM2 in IBM's ESSL library.
C
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     DNRM2: JACK DONGARRA, LINPACK, 3/11/78.
C
#include <implicit.h>
C
      DIMENSION DX(*)
      PARAMETER ( ZERO = 0.0D0 )
C
      DNORM2 = ZERO
      IF(N.LE.0)RETURN
      DTEMP  = ZERO
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DX(IX)
        IX = IX + INCX
   10 CONTINUE
      DNORM2 = SQRT(DTEMP)
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DX(I) + DX(I + 1)*DX(I + 1) +
     *   DX(I + 2)*DX(I + 2) + DX(I + 3)*DX(I + 3) + DX(I + 4)*DX(I + 4)
   50 CONTINUE
   60 DNORM2 = SQRT(DTEMP)
      RETURN
      END
#endif   /* not VAR_ESSL */
C  /* Deck daptge */
      SUBROUTINE DAPTGE(N,AAP,AGE)
C
C  8-Feb-1987 Hans Joergen Aa. Jensen
C  900720-hjaaj: NOSYNC directive
C
C Purpose: Transform from AP format to GE format, that is:
C          unpack antisymmetric,   packed (lower triangle) matrix AAP
C              to antisymmetric, unpacked matrix AGE.
C
#include <implicit.h>
      DIMENSION AAP(*), AGE(N,*)
C
#include <nosync.h>
      DO 200 J = 1,N
         JOFF = (J*J-J)/2
#include <nosync.h>
         DO 100 I = 1,J-1
            AGE(I,J) = - AAP(JOFF+I)
            AGE(J,I) =   AAP(JOFF+I)
  100    CONTINUE
         AGE(J,J) = AAP(JOFF+J)
C        ... is zero but included such that error may be detected.
  200 CONTINUE
C
      RETURN
      END
C  /* Deck dsptsi */
      SUBROUTINE DSPTSI(N,ASP,ASI)
C
C  8-Feb-1987 Hans Joergen Aa. Jensen
C  900720-hjaaj: inserted NOSYNC directive.
C
C Purpose: Transform from SP format to SI format, that is:
C          unpack symmetric,   packed matrix ASP
C              to symmetric, unpacked matrix ASI.
C
#include <implicit.h>
      DIMENSION ASP(*), ASI(N,*)
C
      ENTRY      DSPTGE(N,ASP,ASI)
C     ... equivalent subroutine name
C
#include <nosync.h>
      DO 200 J = 1,N
         JOFF = (J*J-J)/2
#include <nosync.h>
         DO 100 I = 1,J-1
            ASI(I,J) = ASP(JOFF+I)
            ASI(J,I) = ASP(JOFF+I)
  100    CONTINUE
         ASI(J,J) = ASP(JOFF+J)
  200 CONTINUE
C
      RETURN
      END
C  /* Deck dgetap */
      SUBROUTINE DGETAP(N,AGE,AAP)
C
C  8-Feb-1987 Hans Joergen Aa. Jensen
C
C Purpose: Transform from GE format to AP format, that is:
C          extract antisymmetric part of general matrix AGE
C          to antisymmetric, packed matrix AAP (lower
C          triangle saved).
C
#include <implicit.h>
      DIMENSION AGE(N,*), AAP(*)
      PARAMETER ( DP5 = 0.5D0 )
C
#include <nosync.h>
      DO 200 J = 1,N
         JOFF = (J*J-J)/2
         DO 100 I = 1,J
            AAP(JOFF+I) = DP5 * (AGE(J,I) - AGE(I,J))
  100    CONTINUE
  200 CONTINUE
C
      RETURN
      END
C  /* Deck dgetsp */
      SUBROUTINE DGETSP(N,AGE,ASP)
C
C  8-Feb-1987 Hans Joergen Aa. Jensen
C
C Purpose: Transform from GE format to SP format, that is:
C          extract symmetric part of general matrix AGE
C          to symmetric, packed matrix ASP.
C
#include <implicit.h>
      DIMENSION AGE(N,*), ASP(*)
      PARAMETER ( DP5 = 0.5D0 )
C
#include <nosync.h>
      DO 200 J = 1,N
         JOFF = (J*J-J)/2
         DO 100 I = 1,J
            ASP(JOFF+I) = DP5 * (AGE(I,J) + AGE(J,I))
  100    CONTINUE
  200 CONTINUE
C
      RETURN
      END
C  /* Deck dgefsp */
      SUBROUTINE DGEFSP(N,AGE,ASP)
C
C  3-Nov-1989 Hans Joergen Aa. Jensen
C
C Purpose: Fold from GE format to SP format, that is:
C          ASP(ij) = AGE(I,J) + (1 - DELTAij) AGE(J,I)
C
#include <implicit.h>
      DIMENSION AGE(N,*), ASP(*)
C
#include <nosync.h>
      DO 200 J = 1,N
         JOFF = (J*J-J)/2
         DO 100 I = 1,J-1
            ASP(JOFF+I) = AGE(I,J) + AGE(J,I)
  100    CONTINUE
         ASP(JOFF+J) = AGE(J,J)
  200 CONTINUE
C
      RETURN
      END
C  /* Deck dgeasp */
      SUBROUTINE DGEASP(N,AGE,ASP)
C
C 4-Dec-1991 : = DGETSP but adds to SP matrix
C
C Purpose: Transform from GE format to SP format, that is:
C          extract symmetric part of general matrix AGE
C          to symmetric, packed matrix ASP.
C
#include <implicit.h>
      DIMENSION AGE(N,*), ASP(*)
      PARAMETER ( DP5 = 0.5D0 )
C
#include <nosync.h>
      DO 200 J = 1,N
         JOFF = (J*J-J)/2
         DO 100 I = 1,J
            ASP(JOFF+I) = ASP(JOFF+I) + DP5 * (AGE(I,J) + AGE(J,I))
  100    CONTINUE
  200 CONTINUE
C
      RETURN
      END
C  /* Deck dunfld */
      SUBROUTINE DUNFLD(N,ASP,AGE)
C
C  2-Dec-1991 Hans Agren
C
C Purpose: Unfold from SP format to GE format, that is:
C          AGE(I,J) = AGE(J,I) = ASP(ij)/(2.D0 - Delta(I,J))
C
#include <implicit.h>
      DIMENSION AGE(N,*), ASP(*)
      PARAMETER ( DP5 = 0.5D0)
C
#include <nosync.h>
      DO 200 J = 1,N
         JOFF = (J*J-J)/2
         DO 100 I = 1,J-1
            X = ASP(JOFF+I)*DP5
            AGE(I,J) = X
            AGE(J,I) = X
  100    CONTINUE
         AGE(J,J) = ASP(JOFF+J)
  200 CONTINUE
C
      RETURN
      END
C  /* Deck dsitsp */
      SUBROUTINE DSITSP(N,ASI,ASP)
C
C  8-Feb-1987 Hans Joergen Aa. Jensen
C
C Purpose: Transform from SI format to SP format, that is:
C          copy upper triangle of symmetric matrix ASI
C          to symmetric, packed matrix ASP.
C
#include <implicit.h>
      DIMENSION ASI(N,*), ASP(*)
C
#include <nosync.h>
      DO 200 J = 1,N
         JOFF = (J*J-J)/2
         DO 100 I = 1,J
            ASP(JOFF+I) = ASI(I,J)
  100    CONTINUE
  200 CONTINUE
C
      RETURN
      END
C  /* Deck dsifsp */
      SUBROUTINE DSIFSP(N,ASI,ASP)
C
C  3-Nov-1989 Hans Joergen Aa. Jensen
C
C Purpose: Fold from SI format to SP format, that is:
C          ASP(ij) = ASI(I,J) + (1 - DELTAij) ASI(J,I)
C                  = (2 - DELTAij) * ASI(I,J)
C
#include <implicit.h>
      DIMENSION ASI(N,*), ASP(*)
      PARAMETER (D2 = 2.0D0)
C
#include <nosync.h>
      DO 200 J = 1,N
         JOFF = (J*J-J)/2
         DO 100 I = 1,J-1
            ASP(JOFF+I) = D2*ASI(I,J)
  100    CONTINUE
         ASP(JOFF+J) = ASI(J,J)
  200 CONTINUE
C
      RETURN
      END
C  /* Deck dgetsi */
      SUBROUTINE DGETSI(N,AGE,ASI)
C
C  3-Nov-1997 Hans Joergen Aa. Jensen
C
C Purpose: extract symmetric part of general matrix AGE
C          to symmetric matrix ASI.
C
#include <implicit.h>
      DIMENSION AGE(N,*), ASI(N,*)
      PARAMETER ( DP5 = 0.5D0 )
C
#include <nosync.h>
      DO 200 J = 1,N
         ASI(J,J) = AGE(J,J)
         DO 100 I = 1,J-1
            ASI(J,I) = DP5 * (AGE(I,J) + AGE(J,I))
            ASI(I,J) = ASI(J,I)
  100    CONTINUE
  200 CONTINUE
C
      RETURN
      END
C  /* Deck dgetrn */
      SUBROUTINE DGETRN(AGE,NROWA,NRDIMA)
C
C 3-Apr-1987 HJAaJ (DGETRS)
C 900108-hjaaj: block with NBLK for reduced paging
C   when virtual memory
C new name 971103-hjaaj (old name DGETRS was same as 
C   a routine in LAPACK for solving linear equations;
C   when linking with complib on SGI/IRIX the LAPACK routine
C   was loaded instead of this one).
C
C Replace AGE by AGE(transposed)
C
#include <implicit.h>
      DIMENSION AGE(NRDIMA,*)
      PARAMETER (NBLK = 128)
      DO 400 JBLK = 1,NROWA,NBLK
         JEND = MIN(NROWA,JBLK-1+NBLK)
         DO 300 IBLK = 1,JBLK,NBLK
            IEND = MIN(NROWA,IBLK-1+NBLK)
            DO 200 J = JBLK,JEND
               DO 100 I = IBLK,MIN(J-1,IEND)
                  SWAP     = AGE(I,J)
                  AGE(I,J) = AGE(J,I)
                  AGE(J,I) = SWAP
  100          CONTINUE
  200       CONTINUE
  300    CONTINUE
  400 CONTINUE
      RETURN
      END
#if !defined (VAR_DXML)
C  DXML library on DEC-ALPHA
C  /* Deck dsum */
      FUNCTION DSUM(N,DA,INCA)
C
C     Sums elements of a vector.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C 30-Apr-1984 -- hjaaj -- based on DDOT from LINPACK
C     DDOT: JACK DONGARRA, LINPACK, 3/11/78.
C
#include <implicit.h>
C
      DIMENSION DA(*)
      PARAMETER ( D0 = 0.0D0 )
C
      DSUM  = D0
      IF(N.LE.0)RETURN
      DTEMP = D0
      IF(INCA.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IA = 1
      IF(INCA.LT.0)IA = (-N+1)*INCA + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DA(IA)
        IA = IA + INCA
   10 CONTINUE
      DSUM = DTEMP
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DA(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
         DTEMP = DTEMP     + DA(I)     + DA(I + 1)
     *         + DA(I + 2) + DA(I + 3) + DA(I + 4)
   50 CONTINUE
   60 DSUM = DTEMP
      RETURN
      END
C  /* Deck idamin */
      INTEGER FUNCTION IDAMIN(N,DX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MIN. ABSOLUTE VALUE.
C     890927-Hans Joergen Aa. Jensen
C     Based on IDAMAX by
C     JACK DONGARRA, LINPACK, 3/11/78.
C
#include <implicit.h>
C
      DIMENSION DX(*)
C
      IDAMIN = 0
      IF( N .LT. 1 ) RETURN
      IDAMIN = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      DMIN = DABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
         IF(DABS(DX(IX)).GE.DMIN) GO TO 5
         IDAMIN = I
         DMIN = DABS(DX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DMIN = DABS(DX(1))
      DO 30 I = 2,N
         IF(DABS(DX(I)).GE.DMIN) GO TO 30
         IDAMIN = I
         DMIN = DABS(DX(I))
   30 CONTINUE
      RETURN
      END
C  /* Deck idmax */
      INTEGER FUNCTION IDMAX(N,DX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. VALUE.
C     890105 hjaaj, based on IDAMAX by JACK DONGARRA, LINPACK, 3/11/78.
C
#include <implicit.h>
C
      DIMENSION DX(*)
C
      IDMAX = 0
      IF( N .LT. 1 ) RETURN
      IDMAX = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      DMAX = DX(1)
      IX = IX + INCX
      DO 10 I = 2,N
         IF(DX(IX).LE.DMAX) GO TO 5
         IDMAX = I
         DMAX = DX(IX)
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DMAX = DX(1)
      DO 30 I = 2,N
         IF(DX(I).LE.DMAX) GO TO 30
         IDMAX = I
         DMAX = DX(I)
   30 CONTINUE
      RETURN
      END
C  /* Deck idmin */
      INTEGER FUNCTION IDMIN(N,DX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MIN. VALUE.
C     890105 hjaaj, based on IDAMAX by JACK DONGARRA, LINPACK, 3/11/78.
C
#include <implicit.h>
C
      DIMENSION DX(*)
C
      IDMIN = 0
      IF( N .LT. 1 ) RETURN
      IDMIN = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      DMIN = DX(1)
      IX = IX + INCX
      DO 10 I = 2,N
         IF(DX(IX).GE.DMIN) GO TO 5
         IDMIN = I
         DMIN = DX(IX)
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DMIN = DX(1)
      DO 30 I = 2,N
         IF(DX(I).GE.DMIN) GO TO 30
         IDMIN = I
         DMIN = DX(I)
   30 CONTINUE
      RETURN
      END
#endif
C  /* Deck iblas */
      SUBROUTINE ICOPY(N,IX,INCX,IY,INCY)
C
C     COPY integer IX TO integer IY.
C     FOR I = 0 TO N-1, COPY IX(LX+I*INCX) TO IY(LY+I*INCY),
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
C     (860516 - hjaaj - based on BLAS DCOPY)
C
      INTEGER IX(*),IY(*)
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) THEN
         IF(INCX.EQ.1) GOTO 20
         IF(INCX.GT.1) GOTO 60
      END IF
    5 CONTINUE
C
C        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      JX = 1
      JY = 1
      IF(INCX.LT.0)JX = (-N+1)*INCX + 1
      IF(INCY.LT.0)JY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        IY(JY) = IX(JX)
        JX = JX + INCX
        JY = JY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7.
C
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        IY(I) = IX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        IY(I) = IX(I)
        IY(I + 1) = IX(I + 1)
        IY(I + 2) = IX(I + 2)
        IY(I + 3) = IX(I + 3)
        IY(I + 4) = IX(I + 4)
        IY(I + 5) = IX(I + 5)
        IY(I + 6) = IX(I + 6)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS=N*INCX
          DO 70 I=1,NS,INCX
          IY(I) = IX(I)
   70     CONTINUE
      RETURN
      END
      SUBROUTINE ISCAL(N,IA,IX,INCX)
C
C     Scale integer vector IX with IA
C     FOR I = 0 TO N-1, SCALE IX(LX+I*INCX) WITH IA
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N
C
C     (901219 - hjaaj - based on ICOPY)
C
      INTEGER IX(*)
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1) GOTO 20
      IF(INCX.GT.1) GOTO 60
    5 CONTINUE
C
C        CODE FOR NONPOSITIVE INCREMENT.
C
      JX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        IX(JX) = IA*IX(JX)
        JX = JX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7.
C
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        IX(I) = IA*IX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        IX(I) = IA*IX(I)
        IX(I + 1) = IA*IX(I + 1)
        IX(I + 2) = IA*IX(I + 2)
        IX(I + 3) = IA*IX(I + 3)
        IX(I + 4) = IA*IX(I + 4)
        IX(I + 5) = IA*IX(I + 5)
        IX(I + 6) = IA*IX(I + 6)
   50 CONTINUE
      RETURN
C
C        CODE FOR  POSITIVE, NONUNIT INCREMENT.
C
   60 CONTINUE
      NS=N*INCX
          DO 70 I=1,NS,INCX
          IX(I) = IA*IX(I)
   70     CONTINUE
      RETURN
      END
      SUBROUTINE ISWAP(N,IX,INCX,IY,INCY)
C
C     Swap integer arrays IX and IY.
C     FOR I = 0 TO N-1, SWAP IX(LX+I*INCX) WITH IY(LY+I*INCY),
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
C     (901219 - hjaaj - based on ICOPY)
C
      INTEGER IX(*),IY(*)
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY .AND. INCX .GT. 0) GO TO 60
C
C        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      JX = 1
      JY = 1
      IF(INCX.LT.0)JX = (-N+1)*INCX + 1
      IF(INCY.LT.0)JY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        IHOLD  = IY(JY)
        IY(JY) = IX(JX)
        IX(JY) = IHOLD
        JX = JX + INCX
        JY = JY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE INCREMENTS.
C
   60 CONTINUE
      NS=N*INCX
      DO 70 I=1,NS,INCX
         IHOLD = IY(I)
         IY(I) = IX(I)
         IX(I) = IHOLD
   70 CONTINUE
      RETURN
      END
C  /* Deck dtmdv */
      SUBROUTINE  DTMDV (N,DA,INCA,DX,INCX,DB,INCB)
C
C     (TRIANGULAR MATRIX DOTTED INTO VECTOR)
C     (06-OCT-80) (1-May-1984/HJAaJ, taken from MOLPRO by NHFB)
C
C     THIS ROUTINE IS AN ADDITION TO THE BASIC LINEAR ALGEBRA
C     SUBROUTINES (THE BLAS) DESIGNED TO EFFICIENTLY COMPUTE
C     THE PRODUCT OF A LOWER-TRIANGULAR PACKED SYMMETRIC MATRIX
C     WITH A VECTOR.  THE ALGORITHM MAKES A PURELY SEQUENTIAL
C     TRAVERSAL OF THE SYMMETRIC MATRIX AND OF BOTH VECTORS.
C     SPECIFICALLY, IT COMPUTES B = A*X FROM
C
C     B(K) = SUM(L=1,N) A(KL) * X(L)      (K = 1,N)
C
C     WHERE THE (K,L) ELEMENT OF A(*,*) IS FOUND AT POSITION
C     KL = (K*(K-1))/2 + L (K .GE. L).  THE VECTORS B(*) AND X(*)
C     ARE STORED IN DB(*) AND DX(*) ACCORDING TO
C
C     B(K) IS AT DB(1 + (K-1)*INCB) IF INCB .GE. 0
C          OR AT DB(1 + (K-N)*INCB) IF INCB .LT. 0
C
C     X(L) IS AT DX(1 + (L-1)*INCX) IF INCX .GE. 0
C          OR AT DX(1 + (L-N)*INCX) IF INCX .LT. 0
C
C     A(KL) IS AT DA(1 + (KL-1)*INCA) AND INCA MUST BE NON-NEGATIVE.
C
      INTEGER          IB,          IBINIT,      INCA
      INTEGER          INCB,        INCX,        IX,          IXINIT
      INTEGER          JX,          K,           KL,          KM1
      INTEGER          L,           M,           MOD,         MP1
      INTEGER          N
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
      REAL DA(*),       DB(*),       DX(*),       DXK
      REAL SUM1
#else
      DOUBLE PRECISION DA(*),       DB(*),       DX(*),       DXK
      DOUBLE PRECISION SUM1
#endif
      IF (N .LE. 0) GO TO 90
      IF ((INCX .NE. 1) .OR. (INCA .NE. 1) .OR. (INCB .NE. 1)) GO TO 60
C
C     CASE FOR ALL INCREMENTS = 1
C
      KL = 1
      DB(1) = DA(1)*DX(1)
      IF (N .LE. 1) GO TO 90
      DO 50 K = 2,N
           SUM1 = 0.0
           DXK = DX(K)
           KM1 = K - 1
           M = MOD(KM1,3)
           IF (M .EQ. 0) GO TO 20
           DO 10 L = 1,M
                KL = KL + 1
                DB(L) = DB(L) + DA(KL)*DXK
                SUM1 = SUM1 + DA(KL)*DX(L)
   10      CONTINUE
           IF (KM1 .LT. 3) GO TO 40
   20      MP1 = M + 1
           DO 30 L = MP1,KM1,3
                DB(L  ) = DB(L  ) + DA(KL+1)*DXK
                DB(L+1) = DB(L+1) + DA(KL+2)*DXK
                DB(L+2) = DB(L+2) + DA(KL+3)*DXK
                SUM1 = SUM1 + DA(KL+1)*DX(L) + DA(KL+2)*DX(L+1) +
     X                        DA(KL+3)*DX(L+2)
                KL = KL + 3
   30      CONTINUE
   40      KL = KL + 1
           DB(K) = SUM1 + DA(KL)*DXK
   50 CONTINUE
      GO TO 90
C
C     CASE FOR INCREMENTS .NE. 1
C
   60 IB = 1
      IX = 1
      IF (INCX .LT. 0) IX = (1-N)*INCX + 1
      IF (INCB .LT. 0) IB = (1-N)*INCB + 1
      IBINIT = IB
      IXINIT = IX
      KL = 1
      DB(IB) = DA(KL)*DX(IX)
      IF (N .LE. 1) GO TO 90
      DO 80 K = 2,N
           SUM1 = 0.0
           IX = IX + INCX
           DXK = DX(IX)
           KM1 = K - 1
           JX = IXINIT
           IB = IBINIT
           DO 70 L = 1,KM1
                KL = KL + INCA
                DB(IB) = DB(IB) + DA(KL)*DXK
                SUM1 = SUM1 + DA(KL)*DX(JX)
                JX = JX + INCX
                IB = IB + INCB
   70      CONTINUE
           KL = KL + INCA
           DB(IB) = SUM1 + DA(KL)*DXK
   80 CONTINUE
   90 RETURN
      END
C  /* Deck ndxgta */
      INTEGER FUNCTION NDXGTA(N,A,DX,INCX)
C
C 900319-hjaaj
C
C Return number of elements with absolute value .gt. A
C
#include <implicit.h>
      DIMENSION DX(N)
      IF (A .LT. 0.0D0) THEN
         NUM = N
      ELSE IF (INCX .EQ. 1) THEN
         NUM = 0
         DO 200 I = 1,N
            IF (ABS(DX(I)) .GT. A) NUM = NUM + 1
  200    CONTINUE
      ELSE
         NUM = 0
         IF (INCX.GT.0) THEN
            IX = 1 - INCX
         ELSE
            IX = 1 - N*INCX
         END IF
         DO 300 I = 1,N
            IF (ABS(DX(IX+I*INCX)) .GT. A) NUM = NUM + 1
  300    CONTINUE
      END IF
      NDXGTA = NUM
      RETURN
      END
C  /* Deck dgeinv */
      SUBROUTINE DGEINV(N,A,AINV,IPVT,WRK,INFO)
C
C 850618-HJAAJ
C
C Call Linpack routines to calculate the inverse of a
C general matrix A.
C
      INTEGER N, IPVT(*), N2
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL) || defined (SYS_T90)
      REAL
#else
      DOUBLE PRECISION
#endif
     &   A(*), AINV(*), WRK(*), DET(2)
C
      N2 = N*N
      CALL DCOPY(N2,A,1,AINV,1)
      CALL DGEFA(AINV,N,N,IPVT,INFO)
      IF (INFO .EQ. 0) CALL DGEDI(AINV,N,N,IPVT,DET,WRK,01)
      RETURN
      END
C  /* Deck dgesol */
      SUBROUTINE DGESOL (NSIM,N,LDA,LDB,A,B,KPVT,INFO)
C
C Written 22-Mar-1985 Hans Joergen Aa. Jensen
C No revisions.
C
C Purpose:
C  Solve the NSIM simultaneous eqautions:
C
C     B(n,nsim) := A(n,n) inverse * B(n,nsim)
C
C  using LINPACK routines DGEFA and DGESL.
C
C  A    is the matrix (general, non-singular)
C  KPVT is a scratch array of length N.
C
#include <implicit.h>
      DIMENSION A(LDA,*),B(LDB,*),KPVT(*)
C
      CALL DGEFA (A,LDA,N,KPVT,INFO)
      IF (INFO.NE.0) RETURN
C
      DO 100 J = 1,NSIM
         CALL DGESL (A,LDA,N,KPVT,B(1,J),0)
  100 CONTINUE
C
      RETURN
      END
C  /* Deck dspsol */
      SUBROUTINE DSPSOL (N,NSIM,AP,B,KPVT,INFO)
C
C Written 8-Feb-1985 Hans Joergen Aa. Jensen
C No revisions.
C
C Purpose:
C  Solve the NSIM simultaneous eqautions:
C
C     B(n,nsim) := A(n,n) inverse * B(n,nsim)
C
C  AP is A in packed format.
C  KPVT is a scratch array of length N.
C
#include <implicit.h>
      DIMENSION AP(*),B(N,*),KPVT(*)
C
      CALL DSPFA (AP,N,KPVT,INFO)
      IF (INFO.NE.0) RETURN
C
      DO 100 J = 1,NSIM
        CALL DSPSL (AP,N,KPVT,B(1,J))
  100 CONTINUE
C
      RETURN
      END
C  /* Deck dspsli */
      SUBROUTINE DSPSLI (N,NSIM,AP,B,KPVT,INFO,DET,INERT)
C
C Written 24-Feb-1989 Hans Joergen Aa. Jensen, based on DSPSOL.
C No revisions.
C
C Purpose:
C  Solve the NSIM simultaneous eqautions:
C
C     B(n,nsim) := A(n,n) inverse * B(n,nsim)
C
C  AP is A in packed format.
C  KPVT is a scratch array of length N.
C
#include <implicit.h>
      DIMENSION AP(*),B(N,*),KPVT(*),DET(2),INERT(3)
C
      CALL DSPFA (AP,N,KPVT,INFO)
      IF (INFO.NE.0) RETURN
C
      CALL DSPDI(AP,N,KPVT,DET,INERT,DUMMY,110)
C     CALL DSPDI(AP,N,KPVT,DET,INERT,WORK,JOB)
C
      DO 100 J = 1,NSIM
        CALL DSPSL (AP,N,KPVT,B(1,J))
  100 CONTINUE
C
      RETURN
      END
C --- end of linextra.F ---
