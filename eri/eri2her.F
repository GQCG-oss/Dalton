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
C  /* Deck erigam */
      SUBROUTINE ERIGAM(RJ000,FACINT,ALPHA,COORPQ,WORK,LWORK,IPRINT)
C     Two-electron Dirac delta added (WK/UniKA/22-11-2002).
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
      DIMENSION RJ000(NPPX,0:JMAX), FACINT(NPPX),
     &          ALPHA(NPPX), COORPQ(NPPX,3),
     &          WORK(LWORK)
#include <ericom.h>
#include <drw2el.h>
C-----------------------------------------------------------------------
C
C     Allocate ERIGAM
C
      KWVALU = 1
      KNDADR = KWVALU +  NPPX
      KLAST  = KNDADR + (NPPX - 1)/IRAT + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('ERIGAM',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
C
      IF (DO2DAR .OR. S4CENT) THEN
         CALL ERIGA3(RJ000,FACINT,ALPHA,COORPQ,WORK(KWVALU),
     &               WORK(KNDADR),WORK(KLAST),LWRK,IPRINT)
      ELSE IF (AD2DAR) THEN
         CALL ERIGA2(RJ000,FACINT,ALPHA,COORPQ,WORK(KWVALU),
     &               WORK(KNDADR),WORK(KLAST),LWRK,IPRINT)
      ELSE
         CALL ERIGA1(RJ000,FACINT,ALPHA,COORPQ,WORK(KWVALU),
     &               WORK(KNDADR),WORK(KLAST),LWRK,IPRINT)
      END IF
C
      RETURN
      END
C  /* Deck eriga1 */
      SUBROUTINE ERIGA1(RJ000,FACINT,ALPHA,COORPQ,
     &                  WVALU,INDADR,WORK,LWORK,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
      PARAMETER (D0 = 0.D0, D1 = 1.D0, DP5 = -0.5D0)
C
      DIMENSION RJ000(NPPX,0:JMAX),FACINT(NPPX),ALPHA(NPPX),
     &          COORPQ(NPPX,3),WVALU(NPPX),INDADR(NPPX),WORK(LWORK)
C
#include <ericom.h>
#include <erithr.h>
C-----------------------------------------------------------------------
C
C     One-center Integrals
C     ====================
C
C     Note: There should be no testing for small integrals since
C     this may in the case of one-center integrals introduce
C     numerical instabilities for large exponents.
C
      IF (IPQXYZ.EQ.7) THEN
         CALL DCOPY(NPPX,FACINT,1,RJ000(1,0),1)
         IF (JMAX .GT. 0) THEN
            CALL DCOPY(NPPX,ALPHA,1,FACINT,1)
            DO 500 J = 1, JMAX
               FAC = D1/(2.0D0*J + 1.0D0)
               DO 510 I = 1, NPPX
                  RJ000(I,J) = FAC*FACINT(I)*RJ000(I,0)
                  FACINT(I)  = ALPHA(I)*FACINT(I)
  510          CONTINUE
  500       CONTINUE
         END IF
C
C     Multicenter Integrals
C     =====================
C
      ELSE
         NODS  = 0
         DO 600 I = 1, NPPX
            IF (ABS(FACINT(I)) .GT. THRSH) THEN
               NODS = NODS + 1
               WVALU(NODS) = DP5*ALPHA(I)
     &            *(COORPQ(I,1)**2 + COORPQ(I,2)**2 + COORPQ(I,3)**2)
               INDADR(NODS) = I
            ELSE
               FACINT(I) = D0
               DO J = 0, JMAX
                  RJ000(I,J) = D0
               END DO
            END IF
  600    CONTINUE
C
C        Calculate gamma function
C        ========================
C
C        Allocate GETGAM
C
         KINDAD = 1
         KWVALS = KINDAD + (3*NPPX - 1)/IRAT + 1
         KFJW   = KWVALS +  3*NPPX
         KREXPW = KFJW   +    NPPX*(JMAX + 1)
         KLAST  = KREXPW +    NPPX
         IF (KLAST .GT. LWORK) CALL STOPIT('ERIGA1',' ',KLAST,LWORK)
         CALL GETGAM(NODS,INDADR,WVALU,RJ000,JMAX,NPPX,WORK(KFJW),
     &               WORK(KINDAD),WORK(KWVALS),WORK(KREXPW),IPRINT)
C
C        Scale gamma function
C        ====================
C
         DO 800 J = 0, JMAX
         DO 800 I = 1, NPPX
            RJ000(I,J) = FACINT(I)*RJ000(I,J)
            FACINT(I)  = ALPHA(I)*FACINT(I)
  800    CONTINUE
      END IF
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Output from ERIGAM',-1)
         WRITE (LUPRI,'(2X,A,I10)') 'JMAX:  ', JMAX
         WRITE (LUPRI,'(2X,A,I10)') 'NPPX:  ', NPPX
         WRITE (LUPRI,'(2X,A,I10)') 'IPQXYZ:', IPQXYZ
         IF (IPRINT .GT. 20) THEN
            CALL HEADER('COORPQ in ERIGAM',-1)
            CALL OUTPUT(COORPQ,1,NPPX,1,3,NPPX,3,1,LUPRI)
            CALL HEADER('Scaled incomplete gamma function in ERIGAM',-1)
            CALL OUTPUT(RJ000(1,0),1,NPPX,1,JMAX+1,NPPX,JMAX+1,1,LUPRI)
         END IF
      END IF
      RETURN
      END
C  /* Deck eriga2 */
      SUBROUTINE ERIGA2(RJ000,FACINT,ALPHA,COORPQ,
     &                  WVALU,INDADR,WORK,LWORK,IPRINT)
C
C     This subroutine adds two-electron Darwin integrals,
C     ---weighted with the perturbation parameter DARFAC---
C     onto the standard repulsion integrals.
C     See also: J. Comput. Chem. 18, 20 (1997), Eq. (31).
C
C     Wim Klopper, University of Karlsruhe, 22 November 2002.
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
      PARAMETER (D0 = 0.D0, D1 = 1.D0, DP5 = -0.5D0, DP25 = 0.25D0)
C
      DIMENSION RJ000(NPPX,0:JMAX),FACINT(NPPX),ALPHA(NPPX),
     &          COORPQ(NPPX,3),WVALU(NPPX),INDADR(NPPX),WORK(LWORK)
C
#include <codata.h>
#include <drw2el.h>
#include <ericom.h>
#include <erithr.h>
C
      DRWFAC=DP25*ALPHA2*DARFAC
C
C     One-center Integrals
C     ====================
C
C     Note: There should be no testing for small integrals since
C     this may in the case of one-center integrals introduce
C     numerical instabilities for large exponents.
C
      IF (IPQXYZ.EQ.7) THEN
         DO 500 J = 0, JMAX
            FAC = D1/(2.0D0*J + 1.0D0)
            DO 510 I = 1, NPPX
               RJ000(I,J) = (FAC + DRWFAC*ALPHA(I))*FACINT(I)
               FACINT(I)  = ALPHA(I)*FACINT(I)
  510       CONTINUE
  500    CONTINUE
C
C     Multicenter Integrals
C     =====================
C
      ELSE
         NODS  = 0
         DO 600 I = 1, NPPX
            NODS = NODS + 1
            WVALU(NODS) = DP5*ALPHA(I)
     &         *(COORPQ(I,1)**2 + COORPQ(I,2)**2 + COORPQ(I,3)**2)
            INDADR(NODS) = I
  600    CONTINUE
C
C        Calculate gamma function
C        ========================
C
C        Allocate GETGAM
C
         KINDAD = 1
         KWVALS = KINDAD + (3*NPPX - 1)/IRAT + 1
         KFJW   = KWVALS +  3*NPPX
         KREXPW = KFJW   +    NPPX*(JMAX + 1)
         KLAST  = KREXPW +    NPPX
         IF (KLAST .GT. LWORK) CALL STOPIT('ERIGA1',' ',KLAST,LWORK)
         CALL GETGAM(NODS,INDADR,WVALU,RJ000,JMAX,NPPX,WORK(KFJW),
     &               WORK(KINDAD),WORK(KWVALS),WORK(KREXPW),IPRINT)
C
C        Scale gamma function
C        ====================
C
         DO 900 I = 1, NPPX
            WVALU(I) = DRWFAC*ALPHA(I)*EXP(-WVALU(I))
  900    CONTINUE
         DO 800 J = 0, JMAX
         DO 800 I = 1, NPPX
            RJ000(I,J) = FACINT(I)*(RJ000(I,J) + WVALU(I))
            FACINT(I)  = ALPHA(I)*FACINT(I)
  800    CONTINUE
      END IF
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Output from ERIGAM',-1)
         WRITE (LUPRI,'(2X,A,I10)') 'JMAX:  ', JMAX
         WRITE (LUPRI,'(2X,A,I10)') 'NPPX:  ', NPPX
         WRITE (LUPRI,'(2X,A,I10)') 'IPQXYZ:', IPQXYZ
         IF (IPRINT .GT. 20) THEN
            CALL HEADER('COORPQ in ERIGAM',-1)
            CALL OUTPUT(COORPQ,1,NPPX,1,3,NPPX,3,1,LUPRI)
            CALL HEADER('Scaled incomplete gamma function in ERIGAM',-1)
            CALL OUTPUT(RJ000(1,0),1,NPPX,1,JMAX+1,NPPX,JMAX+1,1,LUPRI)
         END IF
      END IF
      RETURN
      END
C  /* Deck eriga3 */
      SUBROUTINE ERIGA3(RJ000,FACINT,ALPHA,COORPQ,
     &                  WVALU,INDADR,WORK,LWORK,IPRINT)
C
C     Computation of the two-electron Dirac delta function
C
C     The subroutine either computes the two-electron
C     Darwin integral (DO2DAR = .true.) or the one-electron
C     four-center overlap integral (S4CENT = .true.)
C     See also: J. Comput. Chem. 18, 20 (1997), Eq. (31).
C
C     Wim Klopper, University of Karlsruhe, 22 November 2002.
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
      PARAMETER (D0 = 0.D0, D1 = 1.D0, DP5 = -0.5D0, DP25 = 0.25D0)
C
      DIMENSION RJ000(NPPX,0:JMAX),FACINT(NPPX),ALPHA(NPPX),
     &          COORPQ(NPPX,3),WVALU(NPPX),INDADR(NPPX),WORK(LWORK)
C
#include <drw2el.h>
#include <codata.h>
#include <ericom.h>
#include <erithr.h>
C
      IF (DO2DAR) DRWFAC=DP25*ALPHA2
      IF (S4CENT) DRWFAC=-DP25/PI
C
C     One-center Integrals
C     ====================
C
C     Note: There should be no testing for small integrals since
C     this may in the case of one-center integrals introduce
C     numerical instabilities for large exponents.
C
      IF (IPQXYZ.EQ.7) THEN
         DO 500 J = 1, JMAX
            DO 510 I = 1, NPPX
               RJ000(I,J) = DRWFAC*ALPHA(I)*FACINT(I)
               FACINT(I)  = ALPHA(I)*FACINT(I)
  510       CONTINUE
  500    CONTINUE
C
C     Multicenter Integrals
C     =====================
C
      ELSE
         NODS  = 0
         DO 600 I = 1, NPPX
            NODS = NODS + 1
            WVALU(NODS) = DP5*ALPHA(I)
     &         *(COORPQ(I,1)**2 + COORPQ(I,2)**2 + COORPQ(I,3)**2)
            INDADR(NODS) = I
  600    CONTINUE
         DO 810 I = 1, NPPX
            WVALU(I) = DRWFAC*ALPHA(I)*EXP(-WVALU(I))
  810    CONTINUE
         DO 800 J = 0, JMAX
         DO 800 I = 1, NPPX
            RJ000(I,J) = FACINT(I)*WVALU(I)
            FACINT(I)  = ALPHA(I)*FACINT(I)
  800    CONTINUE
      END IF
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Output from ERIGAM',-1)
         WRITE (LUPRI,'(2X,A,I10)') 'JMAX:  ', JMAX
         WRITE (LUPRI,'(2X,A,I10)') 'NPPX:  ', NPPX
         WRITE (LUPRI,'(2X,A,I10)') 'IPQXYZ:', IPQXYZ
         IF (IPRINT .GT. 20) THEN
            CALL HEADER('COORPQ in ERIGAM',-1)
            CALL OUTPUT(COORPQ,1,NPPX,1,3,NPPX,3,1,LUPRI)
            CALL HEADER('Scaled incomplete gamma function in ERIGAM',-1)
            CALL OUTPUT(RJ000(1,0),1,NPPX,1,JMAX+1,NPPX,JMAX+1,1,LUPRI)
         END IF
      END IF
      RETURN
      END
C  /* Deck erituv */
      SUBROUTINE ERITUV(HERINT,RJ000,COORPQ,INDHER,IODDHR,
     &                  WORK,LWORK,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
      INTEGER T, U, V, TUV
      DIMENSION HERINT(NPPX,NTUV), RJ000(NPPX,0:JMAX), COORPQ(NPPX,3),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP), IODDHR(NRTOP),
     &          WORK(LWORK)
#include <ericom.h>
#include <hertop.h>
C
C     The final R(T,U,V) integrals are arranged as follows:
C
C     R(000)
C     R(100) R(010) R(001)
C     R(200) R(110) R(101) R(020) R(011) R(002)
C     R(300) R(210) R(201) R(120) R(111) R(102) R(030) R(021) R(012)
C                                                             R(003)
C     Special case JMAX = 0
C     =====================
C
      IF (JMAX .EQ. 0) THEN
         CALL DCOPY(NPPX,RJ000,1,HERINT,1)
      ELSE
C
C        Allocate work space
C        ===================
C
         KHRWRK = 1
         KLAST  = KHRWRK + NTUV*NPPX
         IF (KLAST .GT. LWORK) CALL STOPIT('HERI',' ',KLAST,LWORK)
C
C        Recursion loop for Hermite integrals
C        ====================================
C
         IPQX = IPQ0(1) + 1
         IPQY = IPQ0(2) + 1
         IPQZ = IPQ0(3) + 1
C        CALL DZERO(HERINT,NPPX*NTUV)
         DO 200 JVAL = 1, JMAX
            IF (MOD(JMAX-JVAL,2).EQ.0) THEN
               CALL ERIHRC(HERINT,WORK(KHRWRK),JVAL,RJ000,
     &                     COORPQ(1,1),COORPQ(1,2),COORPQ(1,3),
     &                     INDHER,JMAX,MAXDER,NPPX,NTUV,
     &                     IPQX,IPQY,IPQZ)
            ELSE
               CALL ERIHRC(WORK(KHRWRK),HERINT,JVAL,RJ000,
     &                     COORPQ(1,1),COORPQ(1,2),COORPQ(1,3),
     &                     INDHER,JMAX,MAXDER,NPPX,NTUV,
     &                     IPQX,IPQY,IPQZ)
            END IF
  200    CONTINUE
      END IF
C
C     Print section
C     =============
C
      IF (IPRINT .GE. 10) THEN
         CALL TITLER('Output from ERITUV','*',103)
         WRITE (LUPRI,'(2X,A,I10)') 'JMAX  ', JMAX
         WRITE (LUPRI,'(2X,A,I10)') 'NPPX  ', NPPX
         WRITE (LUPRI,'(2X,A,I10)') 'NTUV  ', NTUV
         IF (IPRINT .GE. 20) THEN
            CALL HEADER('Hermite integrals R(t,u,v)',1)
            DO 300 J = 0, JMAX
              DO 320 T = J, 0, -1
                DO 330 U = J - T, 0, -1
                  V = J - T - U
                  TUV = INDHER(T,U,V)
                  IF (IODDHR(TUV) .EQ. 0) THEN
                    WRITE (LUPRI,'(2X,3(A,I1),A,2X,5F12.8/,
     &                                                 (12X,5F12.8))')
     &              'R(',T,',',U,',',V,')', (HERINT(I,TUV),I=1,NPPX)
                  WRITE (LUPRI,'()')
                  END IF
  330           CONTINUE
  320         CONTINUE
  300      CONTINUE
         END IF
      END IF
      RETURN
      END
C  /* Deck erihrc */
      SUBROUTINE ERIHRC(CUR,OLD,JVAL,RJ000,PQX,PQY,PQZ,INDHER,JMAX,
     &                  MAXDER,NPPX,NTUV,IPQX,IPQY,IPQZ)
C
#include <implicit.h>
      INTEGER T, U, V, TUV
      LOGICAL PQXGT0, PQYGT0, PQZGT0
      DIMENSION CUR(NPPX,NTUV), OLD(NPPX,NTUV),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          PQX(NPPX), PQY(NPPX), PQZ(NPPX),
     &          RJ000(NPPX,0:JMAX)
#include <doxyz.h>
#include <hertop.h>
C
#include <ibtfun.h>
C
      PQXGT0 = IPQX .EQ. 1
      PQYGT0 = IPQY .EQ. 1
      PQZGT0 = IPQZ .EQ. 1
C
C     JVAL = 1
C     ========
C
      IF (JVAL .EQ. 1) THEN
         CALL DCOPY(NPPX,RJ000(1,JMAX-1),1,CUR(1,1),1)
         IF (PQXGT0) THEN
            DO 110 I = 1, NPPX
               CUR(I,2) = PQX(I)*RJ000(I,JMAX)
  110       CONTINUE
         END IF
         IF (PQYGT0) THEN
            DO 120 I = 1, NPPX
               CUR(I,3) = PQY(I)*RJ000(I,JMAX)
  120       CONTINUE
         END IF
         IF (PQZGT0) THEN
            DO 130 I = 1, NPPX
               CUR(I,4) = PQZ(I)*RJ000(I,JMAX)
  130       CONTINUE
         END IF
C
C     JVAL > 1
C     ========
C
      ELSE
         MAXT   = JMAX
         MAXU   = JMAX
         MAXV   = JMAX
         IF (.NOT.DOX) MAXT = JMAX - MAXDER
         IF (.NOT.DOY) MAXU = JMAX - MAXDER
         IF (.NOT.DOZ) MAXV = JMAX - MAXDER
C
C        R(0,0,0)
C
         CALL DCOPY(NPPX,RJ000(1,JMAX-JVAL),1,CUR,1)
C
C        R(T,0,0)
C
         IF (PQXGT0) THEN
            DO 200 I = 1, NPPX
               CUR(I,2) = PQX(I)*OLD(I,1)
  200       CONTINUE
            DO 300 T = 2, MIN(MAXT,JVAL)
               TMIN1 = T - 1.0D0
               TUV   = INDHER(T  ,0,0)
               M1T   = INDHER(T-1,0,0)
               M2T   = INDHER(T-2,0,0)
               DO 310 I = 1, NPPX
                  CUR(I,TUV) = PQX(I)*OLD(I,M1T) + TMIN1*OLD(I,M2T)
  310          CONTINUE
  300       CONTINUE
         ELSE
            DO 400 T = 2, MIN(MAXT,JVAL), 2
               TMIN1 = T - 1.0D0
               TUV   = INDHER(T  ,0,0)
               M2T   = INDHER(T-2,0,0)
               DO 410 I = 1, NPPX
                  CUR(I,TUV) = TMIN1*OLD(I,M2T)
  410          CONTINUE
  400       CONTINUE
         END IF
C
C        R(T,U,0)
C
         IF (PQYGT0) THEN
            DO 500 T = 0, MIN(MAXT,JVAL - 1), IPQX
               TUV = INDHER(T,1,0)
               M1U = INDHER(T,0,0)
               DO 510 I = 1, NPPX
                  CUR(I,TUV) = PQY(I)*OLD(I,M1U)
  510          CONTINUE
  500       CONTINUE
            DO 600 U = 2, MIN(MAXU,JVAL)
               UMIN1  = U - 1.0D0
               DO 610 T = 0, MIN(MAXT,JVAL - U), IPQX
                  TUV = INDHER(T,U  ,0)
                  M1U = INDHER(T,U-1,0)
                  M2U = INDHER(T,U-2,0)
                  DO 620 I = 1, NPPX
                     CUR(I,TUV) = PQY(I)*OLD(I,M1U) + UMIN1*OLD(I,M2U)
  620             CONTINUE
  610          CONTINUE
  600       CONTINUE
         ELSE
            DO 700 U = 2, MIN(MAXU,JVAL), 2
               UMIN1  = U - 1.0D0
               DO 710 T = 0, MIN(MAXT,JVAL - U), IPQX
                  TUV = INDHER(T,U  ,0)
                  M2U = INDHER(T,U-2,0)
                  DO 720 I = 1, NPPX
                     CUR(I,TUV) = UMIN1*OLD(I,M2U)
  720             CONTINUE
  710          CONTINUE
  700       CONTINUE
         END IF
C
C        R(T,U,V)
C
         IF (PQZGT0) THEN
            IUMAX  = JVAL - 1
            DO 800 U = 0, MIN(MAXU,IUMAX), IPQY
               DO 810 T = 0, MIN(MAXT,IUMAX - U), IPQX
                  TUV = INDHER(T,U,1)
                  M1V = INDHER(T,U,0)
                  DO 820 I = 1, NPPX
                     CUR(I,TUV) = PQZ(I)*OLD(I,M1V)
  820             CONTINUE
  810          CONTINUE
  800       CONTINUE
            DO 900 V = 2, MIN(MAXV,JVAL)
               VMIN1  = V - 1.0D0
               IUMAX  = JVAL - V
               DO 910 U = 0, MIN(MAXU,IUMAX), IPQY
                  DO 920 T = 0, MIN(MAXT,IUMAX - U), IPQX
                     TUV = INDHER(T,U,V  )
                     M1V = INDHER(T,U,V-1)
                     M2V = INDHER(T,U,V-2)
                     DO 930 I = 1, NPPX
                        CUR(I,TUV) = PQZ(I)*OLD(I,M1V)+VMIN1*OLD(I,M2V)
  930                CONTINUE
  920             CONTINUE
  910          CONTINUE
  900       CONTINUE
         ELSE
            DO 1000 V = 2, MIN(MAXV,JVAL), 2
               VMIN1  = V - 1.0D0
               IUMAX  = JVAL - V
               DO 1010 U = 0, MIN(MAXU,IUMAX), IPQY
                  DO 1020 T = 0, MIN(MAXT,IUMAX - U), IPQX
                     TUV = INDHER(T,U,V  )
                     M2V = INDHER(T,U,V-2)
                     DO 1030 I = 1, NPPX
                        CUR(I,TUV) = VMIN1*OLD(I,M2V)
 1030                CONTINUE
 1020             CONTINUE
 1010          CONTINUE
 1000       CONTINUE
         END IF
      END IF
      RETURN
      END
C  /* Deck eriswp */
      SUBROUTINE ERISWP(HEROLD,HERNEW,INDHER,IODDHR,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
      INTEGER T, U, V, TUV
      DIMENSION INDHER(0:JTOP,0:JTOP,0:JTOP), IODDHR(NRTOP)
      DIMENSION HEROLD(NPQBCX,NPRFAB,NPRFCD,NTUV),
     &          HERNEW(NPQBCX,NPRFCD,NPRFAB,NTUV)
#include <ericom.h>
#include <hertop.h>
C
      DO 100 J = 0, JMAX
         DO 110 T = J, 0, -1
            DO 120 U = J - T, 0, -1
               V = J - T - U
               TUV = INDHER(T,U,V)
               IF (IODDHR(TUV) .EQ. 0) THEN
                  DO 200 K = 1, NPRFAB
                  DO 200 L = 1, NPRFCD
                     DO 210 I = 1, NPQBCX
                        HERNEW(I,L,K,TUV) = HEROLD(I,K,L,TUV)
  210                CONTINUE
  200             CONTINUE
               END IF
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C     Print section
C     =============
C
      IF (IPRINT .GE. 10) THEN
         CALL TITLER('Output from ERISWP','*',103)
         WRITE (LUPRI,'(2X,A,I10)') 'JMAX  ', JMAX
         WRITE (LUPRI,'(2X,A,I10)') 'NPPX  ', NPPX
         WRITE (LUPRI,'(2X,A,I10)') 'NTUV  ', NTUV
         IF (IPRINT .GE. 20) THEN
            CALL HEADER('Transposed Hermite integrals R(t,u,v)',1)
            DO 300 J = 0, JMAX
              DO 320 T = J, 0, -1
                DO 330 U = J - T, 0, -1
                  V = J - T - U
                  TUV = INDHER(T,U,V)
                  IF (IODDHR(TUV) .EQ. 0) THEN
                    WRITE (LUPRI,'(2X,3(A,I1),A,2X,5F12.8/,
     &                                                 (12X,5F12.8))')
     &              'R(',T,',',U,',',V,')', (HERNEW(I,1,1,TUV),I=1,NPPX)
                  WRITE (LUPRI,'()')
                  END IF
  330           CONTINUE
  320         CONTINUE
  300      CONTINUE
         END IF
      END IF
      RETURN
      END
