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
C  /* Deck xerox */
      SUBROUTINE XEROX(N,DX,DY)
C
C     dy(i) = dx(i)
C
#include <implicit.h>
      DIMENSION DX(N), DY(N)
C
      DO I = 1, N
         DY(I) = DX(I)
      END DO
C
      RETURN
      END
C  /* Deck dxerox */
      SUBROUTINE DXEROX(N,DX,DY)
C
C     dy(i) = dy(i) + dx(i)
C
#include <implicit.h>
      DIMENSION DX(N), DY(N)
C
      DO I = 1, N
         DY(I) = DY(I) + DX(I)
      END DO
C
      RETURN
      END
C  /* Deck ixerox */
      SUBROUTINE IXEROX(N,IX,IY)
C
C     iy(j) = iy(j) + ix(j)
C
#include <implicit.h>
      DIMENSION IX(N), IY(N)
C
      DO J = 1, N
         IY(J) = IY(J) + IX(J)
      END DO
C
      RETURN
      END
C  /* Deck draise */
      SUBROUTINE DRAISE(N,DX,K,DY)
C
C     dy(i) = dx(i)^k
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D00)
      DIMENSION DX(N), DY(N)
C
      IF (K .EQ. 0) THEN
         DO I = 1, N
            DY(I) = D1
         END DO
      ELSE IF (K .EQ. 1) THEN
         DO I = 1, N
            DY(I) = DX(I)
         END DO
      ELSE IF (K .GT. 1) THEN
         DO I = 1, N
            DY(I) = DX(I)**K
         END DO
      ELSE IF (K .LT. 0) THEN
         WRITE(LUPRI,'(A,I5)')
     &      'Error in DRAISE: Exponent is negative. K=',K
         CALL QUIT('Error in DRAISE: Exponent is negative')
      END IF
C
      RETURN
      END
C  /* Deck caxpi */
      SUBROUTINE CAXPI(N,DA,DX,DY)
C
C     dy(i) = da*dx(i)
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION DX(N), DY(N)
C
      DO I = 1, N
         DY(I) = DA*DX(I)
      END DO
C
      RETURN
      END
C  /* Deck cady */
      SUBROUTINE CADY(N,DA,DX1,DX2,DY)
C
C     dy(i) = da*dx1(i)*dx2(i)
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION DX1(N), DX2(N), DY(N)
C
      DO I = 1, N
         DY(I) = DA*DX1(I)*DX2(I)
      END DO
C
      RETURN
      END
C  /* Deck dady */
      SUBROUTINE DADY(N,DA,DX1,DX2,DY)
C
C     dy(i) = dy(i) + da*dx1(i)*dx2(i)
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION DX1(N), DX2(N), DY(N)
C
      DO I = 1, N
         DY(I) = DY(I) + DA*DX1(I)*DX2(I)
      END DO
C
      RETURN
      END
C  /* Deck caddy */
      SUBROUTINE CADDY(N,DA,DX1,K,DX2,DY)
C
C     dy(i) = da*dx1(i)^k*dx2(i)
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION DX1(N), DX2(N), DY(N)
C
      IF (K .EQ. 0) THEN
         DO I = 1, N
            DY(I) = DA*DX2(I)
         END DO
      ELSE IF (K .EQ. 1) THEN
         DO I = 1, N
            DY(I) = DA*DX2(I)*DX1(I)
         END DO
      ELSE IF (K .GT. 1) THEN
         DO I = 1, N
            DY(I) = DA*DX2(I)*DX1(I)**K
         END DO
      ELSE IF (K .LT. 0) THEN
         WRITE(LUPRI,'(A,I5)')
     &      'Error in CADDY: Exponent is negative. K=',K
         CALL QUIT('Error in CADDY: Exponent is negative')
      END IF
C
      RETURN
      END
C  /* Deck daddy */
      SUBROUTINE DADDY(N,DA,DX1,K,DX2,DY)
C
C     dy(i) = dy(i) + da*dx1(i)^k*dx2(i)
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION DX1(N), DX2(N), DY(N)
C
      IF (K .EQ. 0) THEN
         DO I = 1, N
            DY(I) = DY(I) + DA*DX2(I)
         END DO
      ELSE IF (K .EQ. 1) THEN
         DO I = 1, N
            DY(I) = DY(I) + DA*DX2(I)*DX1(I)
         END DO
      ELSE IF (K .GT. 1) THEN
         DO I = 1, N
            DY(I) = DY(I) + DA*DX2(I)*DX1(I)**K
         END DO
      ELSE IF (K .LT. 0) THEN
         WRITE(LUPRI,'(A,I5)')
     &      'Error in DADDY: Exponent is negative. K=',K
         CALL QUIT('Error in DADDY: Exponent is negative')
      END IF
C
      RETURN
      END
C  /* Deck outpk1 */
      SUBROUTINE OUTPK1(AMATRX,IVW,NVW,NTOT,NROW,KCOL,LUPRI)
C.......................................................................
C Revised 04-Jul-1998 by Paal Dahle
C
C OUTPK1 PRINTS A REAL SYMMETRIC MATRIX STORED IN ROW-PACKED UPPER
C TRIANGULAR FORM (SEE DIAGRAM BELOW) IN FORMATTED FORM WITH NUMBERED
C ROWS AND COLUMNS.  THE INPUT IS AS FOLLOWS:
C
C        AMATRX(NVW,NTOT)....PACKED MATRIX
C
C        IVW.................WHICH NVW TO OUTPUT
C
C        NROW................NUMBER OF ROWS TO BE OUTPUT
C
C        KCOL................NUMBER OF COLUMNS IN OUTPUT (4,5 OR 6)
C
C
C THE MATRIX ELEMENTS ARE ARRANGED IN STORAGE AS FOLLOWS:
C
C        1    2    4    7   11   16   22   29
C             3    5    8   12   17   23   30
C                  6    9   13   18   24   31
C                      10   14   19   25   32
C                           15   20   26    . 
C                                21   27    .
C                                     28    
C
C
C AUTHOR:  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
C          FLORIDA, GAINESVILLE
C..........VERSION = 09/05/73/03
C.......................................................................
C
#include <implicit.h>
      PARAMETER (ZERO=0.D00)
      INTEGER   BEGIN
      LOGICAL   TEST
      CHARACTER COLUMN*6
      CHARACTER*20 KFMT, PFMT, QFMT(7)
      DIMENSION AMATRX(NVW,NTOT)
      DATA COLUMN/'Column'/
C
      J = NROW*(NROW+1)/2
      AMAX = ZERO
      DO I=1,J
         AMAX = MAX( AMAX, ABS(AMATRX(IVW,I)) )
      END DO
      IF (AMAX .EQ. ZERO) THEN
         WRITE (LUPRI,'(/T6,A)') 'Zero matrix.'
         RETURN
      END IF
C
      IF (KCOL .EQ. 3) THEN
         FFMIN = 1.D-10
         FFMAX = 1.D3
         TEST  = FFMIN.LE.AMAX .AND. AMAX.LE.FFMAX
c         test  = .false.
         KFMT  = '(/17X,3(A6,I4,13X))'
         IF (TEST) THEN
C           use F output format
            QFMT(1) = '(I5, 2X,8F23.15)'
            QFMT(2) = '(I5,25X,8F23.15)'
            QFMT(3) = '(I5,48X,8F23.15)'
         ELSE
C           use D output format
            QFMT(1) = '(I5, 2X,8D23.15)'
            QFMT(2) = '(I5,25X,8D23.15)'
            QFMT(3) = '(I5,48X,8D23.15)'
         END IF
      ELSE IF (KCOL .EQ. 4) THEN
         FFMIN = 1.D-8
         FFMAX = 1.D6
         TEST  = FFMIN.LE.AMAX .AND. AMAX.LE.FFMAX
         KFMT  = '(/14X,4(A6,I4,7X))'
         IF (TEST) THEN
C           use F output format
            QFMT(1) = '(I5, 2X,8F17.10)'
            QFMT(2) = '(I5,19X,8F17.10)'
            QFMT(3) = '(I5,36X,8F17.10)'
            QFMT(4) = '(I5,53X,8F17.10)'
         ELSE
C           use D output format
            QFMT(1) = '(I5, 2X,8D17.10)'
            QFMT(2) = '(I5,19X,8D17.10)'
            QFMT(3) = '(I5,36X,8D17.10)'
            QFMT(4) = '(I5,53X,8D17.10)'
         END IF
      ELSE IF (KCOL .EQ. 5) THEN
         FFMIN = 1.D-7
         FFMAX = 1.D6
         TEST  = FFMIN.LE.AMAX .AND. AMAX.LE.FFMAX
         KFMT  = '(/9X,5(A6,I4,4X))'
         IF (TEST) THEN
C           use F output format
            QFMT(1) = '(I4, 1X,8F14.7)'
            QFMT(2) = '(I4,15X,8F14.7)'
            QFMT(3) = '(I4,29X,8F14.7)'
            QFMT(4) = '(I4,43X,8F14.7)'
            QFMT(5) = '(I4,57X,8F14.7)'
         ELSE
C           use D output format
            QFMT(1) = '(I4, 1X,8D14.6)'
            QFMT(2) = '(I4,15X,8D14.6)'
            QFMT(3) = '(I4,29X,8D14.6)'
            QFMT(4) = '(I4,43X,8D14.6)'
            QFMT(5) = '(I4,57X,8D14.6)'
         END IF
      ELSE IF (KCOL .EQ. 6) THEN
         FFMIN = 1.D-6
         FFMAX = 1.D4
         TEST  = FFMIN.LE.AMAX .AND. AMAX.LE.FFMAX
         KFMT  = '(/8X,6(A6,I3,3X))'
         IF (TEST) THEN
C           use F output format
            QFMT(1) = '(I4, 1X,8F12.6)'
            QFMT(2) = '(I4,13X,8F12.6)'
            QFMT(3) = '(I4,25X,8F12.6)'
            QFMT(4) = '(I4,37X,8F12.6)'
            QFMT(5) = '(I4,49X,8F12.6)'
            QFMT(6) = '(I4,61X,8F12.6)'
         ELSE
C           use D output format
            QFMT(1) = '(I4, 1X,8D12.5)'
            QFMT(2) = '(I4,13X,8D12.5)'
            QFMT(3) = '(I4,25X,8D12.5)'
            QFMT(4) = '(I4,37X,8D12.5)'
            QFMT(5) = '(I4,49X,8D12.5)'
            QFMT(6) = '(I4,61X,8D12.5)'
         END IF
      ELSE IF (KCOL .EQ. 7) THEN
         FFMIN = 1.D-5
         FFMAX = 1.D2
         TEST  = FFMIN .LE. AMAX .AND. AMAX .LE. FFMAX
       test = .true.
         KFMT  = '(/7X,7(A6,I2,2X))'
         IF (TEST) THEN
C           use F output format
            QFMT(1) = '(I4, 1X,8F10.5)'
            QFMT(2) = '(I4,11X,8F10.5)'
            QFMT(3) = '(I4,21X,8F10.5)'
            QFMT(4) = '(I4,31X,8F10.5)'
            QFMT(5) = '(I4,41X,8F10.5)'
            QFMT(6) = '(I4,51X,8F10.5)'
            QFMT(7) = '(I4,61X,8F10.5)'
         ELSE
C           use D output format
            QFMT(1) = '(I4, 1X,8D10.3)'
            QFMT(2) = '(I4,11X,8D10.3)'
            QFMT(3) = '(I4,21X,8D10.3)'
            QFMT(4) = '(I4,31X,8D10.3)'
            QFMT(5) = '(I4,41X,8D10.3)'
            QFMT(6) = '(I4,51X,8D10.3)'
            QFMT(7) = '(I4,61X,8D10.3)'
         END IF
      ELSE
         WRITE(LUPRI,'(/1X,A)')
     &        ' Error in OUTPK1: KCOL is different from 4,5,6 or 7'
      END IF
C
C     LAST  is the last  column number in the row currently being printed
C     BEGIN is the first column number in the row currently being printed.
C
      KMAX  = MIN(NROW,KCOL)
      LAST  = MIN(NROW,KCOL)
      BEGIN = 1 
C
      DO WHILE (BEGIN .LE. NROW)
         NCOL  = 1
         WRITE (LUPRI,KFMT) (COLUMN,I,I = BEGIN,LAST)
         DO K = 1, LAST
            PFMT = QFMT(1)        
            IF (K.GT.BEGIN) PFMT = QFMT(K-BEGIN+1)        
            WRITE(LUPRI,PFMT) K,
     &           (AMATRX(IVW,((BEGIN+J-1)*(BEGIN+J-2))/2+K),J=NCOL,KMAX)
            IF (K .GE. BEGIN) NCOL = NCOL + 1
         END DO
         LAST = LAST + KCOL
         IF (LAST.GT.NROW) THEN
            KMAX = KCOL - (LAST-NROW)
            LAST = NROW
         END IF
         BEGIN = BEGIN + KCOL
      END DO
C
      RETURN
      END
C  /* Deck outpu1 */
      SUBROUTINE OUTPU1 (AMATRX,IVW,NVW,ROWDIM,COLDIM,KCOL,LUPRI)
C.......................................................................
C Revised 07-Jul-1998 by Paal Dahle
C
C OUTPUT PRINTS A REAL MATRIX IN FORMATTED FORM WITH NUMBERED ROWS
C AND COLUMNS.  THE INPUT IS AS FOLLOWS;
C
C        AMATRX(',').........MATRIX TO BE OUTPUT
C 
C        IVW.................WHICH NVW TO PRINT
C
C        ROWDIM..............ROW DIMENSION OF AMATRX(',')
C
C        COLDIM..............COLUMN DIMENSION OF AMATRX(',')
C
C        KCOL................NUMBER OF COLUMNS TO USE IN OUTPUT
C
C
C THE PARAMETERS THAT FOLLOW MATRIX ARE ALL OF TYPE INTEGER.  THE
C PROGRAM IS SET UP TO HANDLE 5 COLUMNS/PAGE WITH A 1P,5D24.15 FORMAT
C FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED,
C CHANGE FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER
C OF COLUMNS.
C
C AUTHOR;  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
C          FLORIDA, GAINESVILLE
C REVISED; FEBRUARY 26, 1971
C
C.......................................................................
C
#include <implicit.h>
      PARAMETER (ZERO=0.D00)
      INTEGER   ROWDIM,COLDIM,KCOL,BEGIN
      LOGICAL   TEST
      CHARACTER COLUMN*6
      CHARACTER*20 KFMT, PFMT
      DIMENSION AMATRX(NVW,ROWDIM,COLDIM)
      DATA COLUMN/'Column'/
C
      AMAX = ZERO
      DO J = 1, COLDIM
         DO I = 1, ROWDIM
            AMAX = MAX( AMAX, ABS(AMATRX(IVW,I,J)) )
         END DO
      END DO
      IF (AMAX .EQ. ZERO) THEN
         WRITE (LUPRI,'(/T6,A)') 'Zero matrix.'
         RETURN
      END IF
C
      IF (KCOL .EQ. 3) THEN
         FFMIN = 1.D-10
         FFMAX = 1.D6
         TEST  = FFMIN .LE. AMAX .AND. AMAX .LE. FFMAX
         KFMT  = '(/15X,3(A6,I4,13X))'
         IF (TEST) THEN
C           use F output format
            PFMT = '(I5,8F23.15)'
         ELSE
C           use D output format
            PFMT = '(I5, 2X,8D23.15)'
         END IF
      ELSE IF (KCOL .EQ. 4) THEN
         FFMIN = 1.D-8
         FFMAX = 1.D6
         TEST  = FFMIN .LE. AMAX .AND. AMAX .LE. FFMAX
         KFMT  = '(/15X,4(A6,I4,5X))'
         IF (TEST) THEN
C           use F output format
            PFMT = '(I8, 2X,8F15.8)'
         ELSE
C           use D output format
            PFMT = '(I8, 2X,8D15.8)'
         END IF
      ELSE IF (KCOL .EQ. 5) THEN
         FFMIN = 1.D-7
         FFMAX = 1.D6
         TEST  = FFMIN .LE. AMAX .AND. AMAX .LE. FFMAX
         KFMT  = '(/9X,5(A6,I4,4X))'
         IF (TEST) THEN
C           use F output format
            PFMT = '(I4, 1X,8F14.7)'
         ELSE
C           use D output format
            PFMT = '(I4, 1X,8D14.6)'
         END IF
      ELSE IF (KCOL .EQ. 6) THEN
         FFMIN = 1.D-6
         FFMAX = 1.D4
         TEST  = FFMIN .LE. AMAX .AND. AMAX .LE. FFMAX
         KFMT  = '(/8X,6(A6,I3,3X))'
         IF (TEST) THEN
C           use F output format
            PFMT = '(I4, 1X,8F12.6)'
         ELSE
C           use D output format
            PFMT = '(I4, 1X,8D12.5)'
         END IF
      ELSE IF (KCOL .EQ. 7) THEN
         FFMIN = 1.D-5
         FFMAX = 1.D2
         TEST  = FFMIN .LE. AMAX .AND. AMAX .LE. FFMAX
         KFMT  = '(/7X,7(A6,I2,2X))'
         IF (TEST) THEN
C           use F output format
            PFMT = '(I4, 1X,8F10.5)'
         ELSE
C           use D output format
            PFMT = '(I4, 1X,8D10.3)'
         END IF
      ELSE
         WRITE(LUPRI,'(/1X,A)')
     &        ' Error in OUTPU1: KCOL is different from 4,5,6 or 7'
      END IF
C
      LAST = MIN(COLDIM,KCOL)
      DO BEGIN = 1, COLDIM,KCOL
         WRITE (LUPRI,KFMT) (COLUMN,I,I = BEGIN,LAST)
         DO K = 1, ROWDIM
            WRITE (LUPRI,PFMT) K,(AMATRX(IVW,K,I), I = BEGIN,LAST)
         END DO
         LAST = MIN(LAST+KCOL,COLDIM)
      END DO
C
      RETURN
      END
C  /* Deck getname */
      SUBROUTINE GETNAME(HSTNAM)
#include <implicit.h>
      CHARACTER*(20) HSTNAM
C
#if defined (SYS_HAL) || defined (VAR_G95)
      IRET = HOSTNM(HSTNAM)
#else
      CALL HOSTNM(HSTNAM)
#endif
C
      RETURN
      END
C  /* Deck getdate */
      SUBROUTINE GETDATE
#include <implicit.h>
#include <priunit.h>
#if defined (SYS_UNIX) || defined (SYS_DEC) || defined (SYS_IRIX) || defined (SYS_NEC) || defined (SYS_SUN) || defined (SYS_LINUX) || defined (SYS_AIX)
      CHARACTER*(24) FDATE
#endif
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
      CDATE = DATE()
      CTIME = CLOCK()
      WRITE (LUPRI,'(T6,A,A8,1X,A8)')
     *   'Date and time (CRAY) : ',CDATE,CTIME
#endif
#if defined (SYS_HAL)
      CHARACTER*(24) CFDATE
      CALL FDATE(CFDATE)
      WRITE (LUPRI,'(T6,2A)') 'Date and time (HAL)    : ',CFDATE
#endif
#if 0 && defined (SYS_HPUX)
      WRITE (LUPRI,'(T6,2A)') 'Date and time (HP-UX)  : ',FDATE()
#endif
#if defined (SYS_DEC)
      WRITE (LUPRI,'(T6,2A)') 'Date and time (DEC)    : ',FDATE()
#endif
#if defined (SYS_IRIX)
      WRITE (LUPRI,'(T6,2A)') 'Date and time (IRIX)   : ',FDATE()
#endif
#if defined (SYS_NEC)
      WRITE (LUPRI,'(T6,2A)') 'Date and time (NEC-SX4): ',FDATE()
#endif
#if defined (SYS_SUN)
      WRITE (LUPRI,'(T6,2A)') 'Date and time (SUN)    : ',FDATE()
#endif
#if defined (SYS_LINUX)
      WRITE (LUPRI,'(T6,2A)') 'Date and time (Linux)  : ',FDATE()
#endif
#if defined (SYS_AIX)
      WRITE (LUPRI,'(T6,2A)') 'Date and time (IBM-AIX): ',FDATE()
#endif
C
      RETURN
      END
