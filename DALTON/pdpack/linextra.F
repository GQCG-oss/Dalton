C
C
C...   Dalton, Release DALTON2016
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
C 1) DNORM2 (emulate ESSL DNORM2: do not use extended precision for intermediates)
C 2) DAPTGE, DGETAP etc. for triangular packing and unpacking of matrices
C            (AP: Antisymmetric Packed, SP, Symmetric Packed, etc.
C            (ex: DAPTGE: Antisymmetric Packed To GEneral matrix)
C            NOTE: packed part is always LOWER triangle (important for AP)
C 3) DSUM, IDAMIN, IDMAX, IDMIN (supplement DASUM and IDAMAX)
C 4) "iblas": ICOPY, ISWAP, ...
C 5) other extra routines: DGEZERO, DSPZERO, NDXGTA
C            (zero block of matrix; find next element .gt. a)
C 6) DSPSOL, DSPSLI, DGESOL, DGEINV : driver routines calling LINPACK
C 7) BLAS and LAPACK routines missing for ESSL:
C            LSAME; ILAENV, DSPTRI, DSPCON
C
C
C
C Block 1:
C DNORM2 (emulate ESSL DNORM2: do not use extended precision for intermediates)
C
#if !defined (VAR_ESSL) && !defined (VAR_DXML)
C DNORM2 is in ESSL library from IBM and DXML library on DEC-ALPHA
C  /* Deck dnorm2 */
      FUNCTION DNORM2(N,DX,INCX)
C
C     Forms the two-norm of a vector.
C 19-Sep-1988 -- hjaaj -- based on DNRM2 from LINPACK
C     This version does not use extended precision for intermediates
C     as the LINPACK version does.C 1) DNORM2 (emulate ESSL DNORM2: do not use extended precision for intermediates)

C     Equivalent to DNORM2 in IBM's ESSL library.
C
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     DNRM2: JACK DONGARRA, LINPACK, 3/11/78.
C
#include "implicit.h"
C
      DIMENSION DX(*)
      PARAMETER ( ZERO = 0.0D0 )
C
      IF (N.LE.0) THEN
         DNORM2 = ZERO
         RETURN
      END IF
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
C Block 2:
C   DAPTGE, DGETAP etc. for triangular packing and unpacking of matrices
C           (AP: Antisymmetric Packed, SP, Symmetric Packed, etc.
C           (ex: DAPTGE: Antisymmetric Packed To GEneral matrix)
C           NOTE: packed part is always LOWER triangle (important for AP)
C  /* Deck daptge */
      SUBROUTINE DAPTGE(N,AAP,AGE)
C
C  8-Feb-1987 Hans Joergen Aa. Jensen
C
C Purpose: Transform from AP format to GE format, that is:
C          unpack antisymmetric,   packed (lower triangle) matrix AAP
C              to antisymmetric, unpacked matrix AGE.
C
#include "implicit.h"
      DIMENSION AAP(*), AGE(N,*)
C
      DO 200 J = 1,N
         JOFF = (J*J-J)/2
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
C
C Purpose: Transform from SP format to SI format, that is:
C          unpack symmetric,   packed matrix ASP
C              to symmetric, unpacked matrix ASI.
C
#include "implicit.h"
      DIMENSION ASP(*), ASI(N,*)
C
      ENTRY      DSPTGE(N,ASP,ASI)
C     ... equivalent subroutine name
C
      DO 200 J = 1,N
         JOFF = (J*J-J)/2
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
Cdgetap
C Purpose: Transform from GE format to AP format, that is:
C          extract antisymmetric part of general matrix AGE
C          to antisymmetric, packed matrix AAP (lower
C          triangle saved).
C
#include "implicit.h"
      DIMENSION AGE(N,*), AAP(*)
      PARAMETER ( DP5 = 0.5D0 )
C
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
#include "implicit.h"
      DIMENSION AGE(N,*), ASP(*)
      PARAMETER ( DP5 = 0.5D0 )
#ifdef LUCI_DEBUG
#include <priunit.h>
#endif
C
      DO 200 J = 1,N
         JOFF = (J*J-J)/2
         DO 100 I = 1,J
            ASP(JOFF+I) = DP5 * (AGE(I,J) + AGE(J,I))
#ifdef LUCI_DEBUG
            write(lupri,*) 'joff,i,j ==> age(i,j), age(j,i,)',joff,
     &                      i,j,AGE(I,J), AGE(J,I),ASP(JOFF+I)
#endif
  100    CONTINUE
  200 CONTINUE
#ifdef LUCI_DEBUG
      write(lupri,*) 'end of DGETSP'
#endif
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
#include "implicit.h"
      DIMENSION AGE(N,*), ASP(*)
C
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
#include "implicit.h"
      DIMENSION AGE(N,*), ASP(*)
      PARAMETER ( DP5 = 0.5D0 )
C
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
#include "implicit.h"
      DIMENSION AGE(N,*), ASP(*)
      PARAMETER ( DP5 = 0.5D0)
C
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
#include "implicit.h"
      DIMENSION ASI(N,*), ASP(*)
C
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
#include "implicit.h"
      DIMENSION ASI(N,*), ASP(*)
      PARAMETER (D2 = 2.0D0)
C
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
C          Can be used in-place ( CALL DGETSI(N,A,A) )
C
#include "implicit.h"
      DIMENSION AGE(N,*), ASI(N,*)
      PARAMETER ( DP5 = 0.5D0 )
C
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
C Replace AGE by AGE(transposed)
C
C 3-Apr-1987 HJAaJ
C 900108-hjaaj: block with NBLK for reduced paging
C   when virtual memory
C new name 971103-hjaaj (old name DGETRS was same as
C   a routine in LAPACK for solving linear equations;
C   when linking with complib on SGI/IRIX the LAPACK routine
C   was loaded instead of this one).
C
#include "implicit.h"
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

C Block 3:
C DSUM, IDAMIN, IDMAX, IDMIN (supplement DASUM and IDAMAX)
C  /* Deck dsum */
      FUNCTION DSUM(N,DA,INCA)
C
C     Sums elements of a vector.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C 30-Apr-1984 -- hjaaj -- based on DDOT from LINPACK
C     DDOT: JACK DONGARRA, LINPACK, 3/11/78.
C
#include "implicit.h"
C
      DIMENSION DA(*)
      PARAMETER ( D0 = 0.0D0 )
C
      IF (N.LE.0) THEN
        DSUM  = D0
        RETURN
      END IF
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
#include "implicit.h"
C
      DIMENSION DX(*)
C
      IF( N .LT. 1 ) THEN
        IDAMIN = 0
        RETURN
      END IF
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
#include "implicit.h"
C
      DIMENSION DX(*)
C
      IF( N .LT. 1 ) THEN
        IDMAX = 0
        RETURN
      END IF
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
#include "implicit.h"
C
      DIMENSION DX(*)
C
      IF( N .LT. 1 ) THEN
        IDMIN = 0
        RETURN
      END IF

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

C Block 4: "iblas": ICOPY, ISCAL, ISWAP

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
CC 5) DSPSOL, DGESOL, DGEINV : driver routines calling LINPACK
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

C Block 5: other routines
C   DGEZERO, DSPZERO, NDXGTA, DUNIT, DZERO, ISUM, IZERO
C

      SUBROUTINE DGEZERO(N,AGE,NRSTA,NREND,NCSTA,NCEND)
C
C     Oct. 09, Hans Joergen Aa. Jensen
C
C     Zero block AGE(NRSTA:NREND,NCSTA:NCEND)
C
      IMPLICIT  NONE
      INTEGER N, NRSTA, NREND, NCSTA, NCEND
      REAL*8  AGE(N,N)
      INTEGER IC, IR

      DO IC = NCSTA, NCEND
         DO IR = NRSTA, NREND
            AGE(IR,IC) = 0.0D0
         END DO
      END DO
      RETURN
      END

      SUBROUTINE DSPZERO(N,ASP,NRSTA,NREND,NCSTA,NCEND)
C
C     Oct. 09, Hans Joergen Aa. Jensen
C
C     Zero block ASP(NRSTA:NREND,NCSTA:NCEND),
C     however, ASP is Symmetric Packed (triangular packed)
C
      IMPLICIT  NONE
      REAL*8  ASP(*)
      INTEGER N, NRSTA, NREND, NCSTA, NCEND
      INTEGER IC, IR, IROFF, MCEND, ICOFF, MREND

C     First the part placed in lower triangle of the full matrix

     
      IF (NCSTA .LE. NREND) THEN
      DO IR = NRSTA, NREND
         IROFF = (IR*IR-IR)/2
         MCEND = MIN(NCEND, IR)
         DO IC = NCSTA, MCEND
            ASP(IROFF+IC) = 0.0D0
         END DO
      END DO
      END IF

C     and then the part placed in uppper triangle

      IF (NRSTA .LE. NCEND) THEN
      DO IC = NCSTA, NCEND
         ICOFF = (IC*IC-IC)/2
         MREND = MIN(NREND, IC)
         DO IR = NRSTA, MREND
            ASP(ICOFF+IR) = 0.0D0
         END DO
      END DO
      END IF

      RETURN
      END

C  /* Deck ndxgta */
      INTEGER FUNCTION NDXGTA(N,A,DX,INCX)
C
C 900319-hjaaj
C
C Return number of elements with absolute value .gt. A
C
#include "implicit.h"
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
C  /* Deck dunit */
      SUBROUTINE DUNIT(A,N)
C
C  SUBROUTINE DUNIT SETS THE REAL SQUARE MATRIX A EQUAL
C  TO A UNIT MATRIX.
C  /VER 2/ 14-Sep-1983 hjaaj
C
#include "implicit.h"
      DIMENSION A(*)
      PARAMETER (D1=1.0D00, D0=0.0D00)
C
      NN = N*N
      DO 100 I = 1,NN
         A(I) = D0
  100 CONTINUE
      N1 = N + 1
      DO 200 I = 1,NN,N1
         A(I) = D1
  200 CONTINUE
      RETURN
      END
C  /* Deck dzero */
      SUBROUTINE DZERO(DX,LENGTH)
#include "implicit.h"
C
C Last revision 5-May-1984 by Hans Jorgen Aa. Jensen
C
C   Subroutine DZERO sets a real array of length *LENGTH*
C   to zero.
C...................................................................
      DIMENSION DX(*)
C
      IF (LENGTH.LE.0) RETURN
C
      DO I = 1,LENGTH
         DX(I) = 0.0D00
      END DO
C
      RETURN
      END
C  /* Deck isum */
      FUNCTION ISUM(N,IA,INCA)
C
C     8-Feb-1987 hjaaj
C     Sums elements of a integer vector.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     -- based on DDOT from LINPACK
C     DDOT: JACK DONGARRA, LINPACK, 3/11/78.
C
      INTEGER ISUM,  IA(*), ITEMP
      INTEGER I,INCA,JA,M,MP1,N
C
      IF(N.LE.0) THEN
        ISUM  = 0
        RETURN
      END IF

      ITEMP = 0
      IF(INCA.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      JA = 1
      IF(INCA.LT.0)JA = (-N+1)*INCA + 1
      DO 10 I = 1,N
        ITEMP = ITEMP + IA(JA)
        JA = JA + INCA
   10 CONTINUE
      ISUM = ITEMP
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
        ITEMP = ITEMP + IA(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
         ITEMP = ITEMP     + IA(I)     + IA(I + 1)
     *         + IA(I + 2) + IA(I + 3) + IA(I + 4)
   50 CONTINUE
   60 ISUM = ITEMP
      RETURN
      END
C  /* Deck izero */
      SUBROUTINE IZERO(INTVEC,LENGTH)
C...................................................................
C Written 5-May-1984 by Hans Jorgen Aa. Jensen
C
C   Subroutine IZERO sets an integer array of length *LENGTH*
C   to zero.
C...................................................................
      INTEGER LENGTH, INTVEC(*)
C
      IF (LENGTH.LE.0) RETURN
C
      DO I=1,LENGTH
         INTVEC(I) = 0
      END DO
C
      END


C Block 6:
C   DGEINV, DGESOL, DSPSOL, DSPSLI : driver routines calling LINPACK

C  /* Deck dgeinv */
      SUBROUTINE DGEINV(N,A,AINV,IPVT,WRK,INFO)
C
C 850618-HJAAJ
C
C Call Linpack routines to calculate the inverse of a
C general matrix A.
C
      INTEGER N, IPVT(*), N2
      REAL*8  A(*), AINV(*), WRK(*), DET(2)
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
#include "implicit.h"
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
#include "implicit.h"
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
#include "implicit.h"
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
C Block 7: BLAS and LAPACK routines missing for ESSL
C            LSAME; ILAENV, DSPTRI, DSPCON
#if defined (VAR_ESSL)
      LOGICAL FUNCTION LSAME(CA,CB)
*
*  -- Reference BLAS level1 routine (version 3.1) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER CA,CB
*     ..
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER INTA,INTB,ZCODE
*     ..
*
*     Test if the characters are equal
*
      LSAME = CA .EQ. CB
      IF (LSAME) RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR('Z')
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR(CA)
      INTB = ICHAR(CB)
*
      IF (ZCODE.EQ.90 .OR. ZCODE.EQ.122) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
          IF (INTA.GE.97 .AND. INTA.LE.122) INTA = INTA - 32
          IF (INTB.GE.97 .AND. INTB.LE.122) INTB = INTB - 32
*
      ELSE IF (ZCODE.EQ.233 .OR. ZCODE.EQ.169) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
          IF (INTA.GE.129 .AND. INTA.LE.137 .OR.
     +        INTA.GE.145 .AND. INTA.LE.153 .OR.
     +        INTA.GE.162 .AND. INTA.LE.169) INTA = INTA + 64
          IF (INTB.GE.129 .AND. INTB.LE.137 .OR.
     +        INTB.GE.145 .AND. INTB.LE.153 .OR.
     +        INTB.GE.162 .AND. INTB.LE.169) INTB = INTB + 64
*
      ELSE IF (ZCODE.EQ.218 .OR. ZCODE.EQ.250) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
          IF (INTA.GE.225 .AND. INTA.LE.250) INTA = INTA - 32
          IF (INTB.GE.225 .AND. INTB.LE.250) INTB = INTB - 32
      END IF
      LSAME = INTA .EQ. INTB
*
*     RETURN
*
*     End of LSAME
*
      END
      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
*
*  -- LAPACK auxiliary routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IC, IZ, NB, NBMIN, NX
      LOGICAL            CNAME, SNAME
      CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*6
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. External Functions ..
      INTEGER            IEEECK, IPARMQ
      EXTERNAL           IEEECK, IPARMQ
*     ..
*     .. Executable Statements ..
*
      GO TO ( 10, 10, 10, 80, 90, 100, 110, 120,
     $        130, 140, 150, 160, 160, 160, 160, 160 )ISPEC
*
*     Invalid value for ISPEC
*
      ILAENV = -1
      RETURN
*
   10 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1: 1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC+64 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I:
     $             I ) = CHAR( IC+64 )
   30       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 40 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   40       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1: 1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      C4 = C3( 2: 3 )
*
      GO TO ( 50, 60, 70 )ISPEC
*
   50 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
*
   60 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.
     $       'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
   70 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.
     $       'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
*
   80 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ILAENV = 6
      RETURN
*
   90 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ILAENV = 2
      RETURN
*
  100 CONTINUE
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  110 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ILAENV = 1
      RETURN
*
  120 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ILAENV = 50
      RETURN
*
  130 CONTINUE
*
*     ISPEC = 9:  maximum size of the subproblems at the bottom of the
*                 computation tree in the divide-and-conquer algorithm
*                 (used by xGELSD and xGESDD)
*
      ILAENV = 25
      RETURN
*
  140 CONTINUE
*
*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
*
*     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 )
      END IF
      RETURN
*
  150 CONTINUE
*
*     ISPEC = 11: infinity arithmetic can be trusted not to trap
*
*     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 )
      END IF
      RETURN
*
  160 CONTINUE
*
*     12 <= ISPEC <= 16: xHSEQR or one of its subroutines. 
*
      ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      RETURN
*
*     End of ILAENV
*
      END
      SUBROUTINE DSPTRI( UPLO, N, AP, IPIV, WORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   AP( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, K, KC, KCNEXT, KP, KPC, KSTEP, KX, NPP
      DOUBLE PRECISION   AK, AKKP1, AKP1, D, T, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DSPMV, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSPTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Check that the diagonal matrix D is nonsingular.
*
      IF( UPPER ) THEN
*
*        Upper triangular storage: examine D from bottom to top
*
         KP = N*( N+1 ) / 2
         DO 10 INFO = N, 1, -1
            IF( IPIV( INFO ).GT.0 .AND. AP( KP ).EQ.ZERO )
     $         RETURN
            KP = KP - INFO
   10    CONTINUE
      ELSE
*
*        Lower triangular storage: examine D from top to bottom.
*
         KP = 1
         DO 20 INFO = 1, N
            IF( IPIV( INFO ).GT.0 .AND. AP( KP ).EQ.ZERO )
     $         RETURN
            KP = KP + N - INFO + 1
   20    CONTINUE
      END IF
      INFO = 0
*
      IF( UPPER ) THEN
*
*        Compute inv(A) from the factorization A = U*D*U**T.
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         K = 1
         KC = 1
   30    CONTINUE
*
*        If K > N, exit from loop.
*
         IF( K.GT.N )
     $      GO TO 50
*
         KCNEXT = KC + K
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Invert the diagonal block.
*
            AP( KC+K-1 ) = ONE / AP( KC+K-1 )
*
*           Compute column K of the inverse.
*
            IF( K.GT.1 ) THEN
               CALL DCOPY( K-1, AP( KC ), 1, WORK, 1 )
               CALL DSPMV( UPLO, K-1, -ONE, AP, WORK, 1, ZERO, AP( KC ),
     $                     1 )
               AP( KC+K-1 ) = AP( KC+K-1 ) -
     $                        DDOT( K-1, WORK, 1, AP( KC ), 1 )
            END IF
            KSTEP = 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Invert the diagonal block.
*
            T = ABS( AP( KCNEXT+K-1 ) )
            AK = AP( KC+K-1 ) / T
            AKP1 = AP( KCNEXT+K ) / T
            AKKP1 = AP( KCNEXT+K-1 ) / T
            D = T*( AK*AKP1-ONE )
            AP( KC+K-1 ) = AKP1 / D
            AP( KCNEXT+K ) = AK / D
            AP( KCNEXT+K-1 ) = -AKKP1 / D
*
*           Compute columns K and K+1 of the inverse.
*
            IF( K.GT.1 ) THEN
               CALL DCOPY( K-1, AP( KC ), 1, WORK, 1 )
               CALL DSPMV( UPLO, K-1, -ONE, AP, WORK, 1, ZERO, AP( KC ),
     $                     1 )
               AP( KC+K-1 ) = AP( KC+K-1 ) -
     $                        DDOT( K-1, WORK, 1, AP( KC ), 1 )
               AP( KCNEXT+K-1 ) = AP( KCNEXT+K-1 ) -
     $                            DDOT( K-1, AP( KC ), 1, AP( KCNEXT ),
     $                            1 )
               CALL DCOPY( K-1, AP( KCNEXT ), 1, WORK, 1 )
               CALL DSPMV( UPLO, K-1, -ONE, AP, WORK, 1, ZERO,
     $                     AP( KCNEXT ), 1 )
               AP( KCNEXT+K ) = AP( KCNEXT+K ) -
     $                          DDOT( K-1, WORK, 1, AP( KCNEXT ), 1 )
            END IF
            KSTEP = 2
            KCNEXT = KCNEXT + K + 1
         END IF
*
         KP = ABS( IPIV( K ) )
         IF( KP.NE.K ) THEN
*
*           Interchange rows and columns K and KP in the leading
*           submatrix A(1:k+1,1:k+1)
*
            KPC = ( KP-1 )*KP / 2 + 1
            CALL DSWAP( KP-1, AP( KC ), 1, AP( KPC ), 1 )
            KX = KPC + KP - 1
            DO 40 J = KP + 1, K - 1
               KX = KX + J - 1
               TEMP = AP( KC+J-1 )
               AP( KC+J-1 ) = AP( KX )
               AP( KX ) = TEMP
   40       CONTINUE
            TEMP = AP( KC+K-1 )
            AP( KC+K-1 ) = AP( KPC+KP-1 )
            AP( KPC+KP-1 ) = TEMP
            IF( KSTEP.EQ.2 ) THEN
               TEMP = AP( KC+K+K-1 )
               AP( KC+K+K-1 ) = AP( KC+K+KP-1 )
               AP( KC+K+KP-1 ) = TEMP
            END IF
         END IF
*
         K = K + KSTEP
         KC = KCNEXT
         GO TO 30
   50    CONTINUE
*
      ELSE
*
*        Compute inv(A) from the factorization A = L*D*L**T.
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         NPP = N*( N+1 ) / 2
         K = N
         KC = NPP
   60    CONTINUE
*
*        If K < 1, exit from loop.
*
         IF( K.LT.1 )
     $      GO TO 80
*
         KCNEXT = KC - ( N-K+2 )
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Invert the diagonal block.
*
            AP( KC ) = ONE / AP( KC )
*
*           Compute column K of the inverse.
*
            IF( K.LT.N ) THEN
               CALL DCOPY( N-K, AP( KC+1 ), 1, WORK, 1 )
               CALL DSPMV( UPLO, N-K, -ONE, AP( KC+N-K+1 ), WORK, 1,
     $                     ZERO, AP( KC+1 ), 1 )
               AP( KC ) = AP( KC ) - DDOT( N-K, WORK, 1, AP( KC+1 ), 1 )
            END IF
            KSTEP = 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Invert the diagonal block.
*
            T = ABS( AP( KCNEXT+1 ) )
            AK = AP( KCNEXT ) / T
            AKP1 = AP( KC ) / T
            AKKP1 = AP( KCNEXT+1 ) / T
            D = T*( AK*AKP1-ONE )
            AP( KCNEXT ) = AKP1 / D
            AP( KC ) = AK / D
            AP( KCNEXT+1 ) = -AKKP1 / D
*
*           Compute columns K-1 and K of the inverse.
*
            IF( K.LT.N ) THEN
               CALL DCOPY( N-K, AP( KC+1 ), 1, WORK, 1 )
               CALL DSPMV( UPLO, N-K, -ONE, AP( KC+( N-K+1 ) ), WORK, 1,
     $                     ZERO, AP( KC+1 ), 1 )
               AP( KC ) = AP( KC ) - DDOT( N-K, WORK, 1, AP( KC+1 ), 1 )
               AP( KCNEXT+1 ) = AP( KCNEXT+1 ) -
     $                          DDOT( N-K, AP( KC+1 ), 1,
     $                          AP( KCNEXT+2 ), 1 )
               CALL DCOPY( N-K, AP( KCNEXT+2 ), 1, WORK, 1 )
               CALL DSPMV( UPLO, N-K, -ONE, AP( KC+( N-K+1 ) ), WORK, 1,
     $                     ZERO, AP( KCNEXT+2 ), 1 )
               AP( KCNEXT ) = AP( KCNEXT ) -
     $                        DDOT( N-K, WORK, 1, AP( KCNEXT+2 ), 1 )
            END IF
            KSTEP = 2
            KCNEXT = KCNEXT - ( N-K+3 )
         END IF
*
         KP = ABS( IPIV( K ) )
         IF( KP.NE.K ) THEN
*
*           Interchange rows and columns K and KP in the trailing
*           submatrix A(k-1:n,k-1:n)
*
            KPC = NPP - ( N-KP+1 )*( N-KP+2 ) / 2 + 1
            IF( KP.LT.N )
     $         CALL DSWAP( N-KP, AP( KC+KP-K+1 ), 1, AP( KPC+1 ), 1 )
            KX = KC + KP - K
            DO 70 J = K + 1, KP - 1
               KX = KX + N - J + 1
               TEMP = AP( KC+J-K )
               AP( KC+J-K ) = AP( KX )
               AP( KX ) = TEMP
   70       CONTINUE
            TEMP = AP( KC )
            AP( KC ) = AP( KPC )
            AP( KPC ) = TEMP
            IF( KSTEP.EQ.2 ) THEN
               TEMP = AP( KC-N+K-1 )
               AP( KC-N+K-1 ) = AP( KC-N+KP-1 )
               AP( KC-N+KP-1 ) = TEMP
            END IF
         END IF
*
         K = K - KSTEP
         KC = KCNEXT
         GO TO 60
   80    CONTINUE
      END IF
*
      RETURN
*
*     End of DSPTRI
*
      END
      SUBROUTINE DSPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, IWORK,
     $                   INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, N
      DOUBLE PRECISION   ANORM, RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), IWORK( * )
      DOUBLE PRECISION   AP( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, IP, KASE
      DOUBLE PRECISION   AINVNM
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLACN2, DSPTRS, XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSPCON', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.LE.ZERO ) THEN
         RETURN
      END IF
*
*     Check that the diagonal matrix D is nonsingular.
*
      IF( UPPER ) THEN
*
*        Upper triangular storage: examine D from bottom to top
*
         IP = N*( N+1 ) / 2
         DO 10 I = N, 1, -1
            IF( IPIV( I ).GT.0 .AND. AP( IP ).EQ.ZERO )
     $         RETURN
            IP = IP - I
   10    CONTINUE
      ELSE
*
*        Lower triangular storage: examine D from top to bottom.
*
         IP = 1
         DO 20 I = 1, N
            IF( IPIV( I ).GT.0 .AND. AP( IP ).EQ.ZERO )
     $         RETURN
            IP = IP + N - I + 1
   20    CONTINUE
      END IF
*
*     Estimate the 1-norm of the inverse.
*
      KASE = 0
   30 CONTINUE
      CALL DLACN2( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
*
*        Multiply by inv(L*D*L**T) or inv(U*D*U**T).
*
         CALL DSPTRS( UPLO, N, 1, AP, IPIV, WORK, N, INFO )
         GO TO 30
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM.NE.ZERO )
     $   RCOND = ( ONE / AINVNM ) / ANORM
*
      RETURN
*
*     End of DSPCON
*
      END
#endif
      real*8 function ddot_128(n,da,inca,db,incb)
!
!     Purpose: accumulate in quadruple precision
!              to avoid round-off errors for big vectors
!     hjaaj July 2012
!
      implicit none
      real*8 da(*), db(*), ddot
      integer n, inca, incb, i
      ! only real*16 for debugging (e.g. PGI compilers do not allow real*16)
      ! real*16 ddot_value 
      real*8 ddot_value

      if (inca .ne. 1 .or. incb .ne. 1) then
         ddot_128 =  ddot(n,da,inca,db,incb)
         return
      end if
      ddot_value = 0
      do i = 1,n
         ddot_value = ddot_value + da(i)*db(i)
      end do
      ddot_128 = ddot_value
      return
      end
C --- end of linextra.F ---
