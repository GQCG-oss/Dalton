C
C...   Copyright (c) 1997 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dalton, an ab initio electronic structure program, Release 1.0
C...   (1997), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, T. Andersen, K. L. Bak, V. Bakken,
C...   O. Christiansen, P. Dahle, E. K. Dalskov, T. Enevoldsen,
C...   H. Heiberg, D. Jonsson, S. Kirpekar, R. Kobayashi, H. Koch,
C...   K. V. Mikkelsen, P. Norman, M. J. Packer, T.Saue,
C...   P. R. Taylor, and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be 
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may 
C...   be distributed outside the research group of the licence holder. 
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence. 
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html 
C
C
C$Id: eri2her.F,v 1.1.1.1 2001-02-08 13:33:13 hjj Exp $
C
#include <single.h>
C  /* Deck erigam */
      SUBROUTINE ERIGAM(RJ000,FACINT,EXPPQ,PQXYZ,WORK,LWORK,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
      DIMENSION RJ000(NPPX,0:JMAX), FACINT(NPPX),
     &          EXPPQ(NPPX,2), PQXYZ(NPPX,3),
     &          WORK(LWORK)
#include <ericom.h>
C-----------------------------------------------------------------------
C
C     Allocate ERIGAM
C
      KALPHA = 1
      KWVALU = KALPHA +  NPPX
      KNDADR = KWVALU +  NPPX
      KLAST  = KNDADR + (NPPX - 1)/IRAT + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('ERIGAM',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
C
      CALL ERIGA1(RJ000,FACINT,EXPPQ(1,1),EXPPQ(1,2),PQXYZ,
     &            WORK(KALPHA),WORK(KWVALU),WORK(KNDADR),
     &            WORK(KLAST),LWRK,IPRINT)
C
      RETURN
      END
C  /* Deck eriga1 */
      SUBROUTINE ERIGA1(RJ000,FACINT,EXPP,EXPQ,PQXYZ,TALPHA,
     &                  WVALU,INDADR,WORK,LWORK,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
      PARAMETER (D0 = 0.D0, D1 = 1.D0, D2 = 2.D0, DP5M = -0.5D0)
C
      DIMENSION RJ000(NPPX,0:JMAX),FACINT(NPPX),
     &          EXPP(NPPX),EXPQ(NPPX),
     &          PQXYZ(NPPX,3),TALPHA(NPPX),
     &          WVALU(NPPX),INDADR(NPPX),WORK(LWORK)
C
#include <ericom.h>
#include <erithr.h>
#include <sdpre.h>
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
         DO 300 I = 1, NPPX
            RJ000(I,0) = FACINT(I)/SQRT(EXPP(I) + EXPQ(I))
  300    CONTINUE
         IF (JMAX .GT. 0) THEN
            DO 400 I = 1, NPPX
               TALPHA(I) = -D2*EXPP(I)*EXPQ(I)/(EXPP(I) + EXPQ(I))
  400       CONTINUE
            CALL DCOPY(NPPX,TALPHA,1,FACINT,1)
            DO 500 J = 1, JMAX
               FAC = D1/SDPRE(2*J + 1)
               DO 510 I = 1, NPPX
                  RJ000(I,J) = FAC*FACINT(I)*RJ000(I,0)
                  FACINT(I)  = TALPHA(I)*FACINT(I)
  510          CONTINUE
  500       CONTINUE
         END IF
C
C     Multicenter Integrals
C     =====================
C
      ELSE
         DO 600 I = 1, NPPX
            TALPHA(I) = -D2*EXPP(I)*EXPQ(I)/(EXPP(I) + EXPQ(I))
            FACINT(I) =       FACINT(I)/SQRT(EXPP(I) + EXPQ(I))
  600    CONTINUE
         NODS  = 0
C
C        Make sure this test is OK!
C
         DO 700 I = 1, NPPX
            IF (ABS(FACINT(I)) .GT. THRSH) THEN
               NODS = NODS + 1
               WVALU(NODS) = DP5M*TALPHA(I)
     &            *(PQXYZ(I,1)**2 + PQXYZ(I,2)**2 + PQXYZ(I,3)**2)
               INDADR(NODS) = I
            ELSE
               TALPHA(I) = D0
            END IF
  700    CONTINUE
C
         CALL DZERO(RJ000(1,0),(JMAX+1)*NPPX)
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
            FACINT(I)  = TALPHA(I)*FACINT(I)
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
            CALL HEADER('TALPHA in ERIGAM',-1)
            CALL OUTPUT(TALPHA,1,1,1,NPPX,1,NPPX,1,LUPRI)
            CALL HEADER('PQXYZ in ERIGAM',-1)
            CALL OUTPUT(PQXYZ,1,NPPX,1,3,NPPX,3,1,LUPRI)
            CALL HEADER('Scaled incomplete gamma function in ERIGAM',-1)
            CALL OUTPUT(RJ000(1,0),1,NPPX,1,JMAX+1,NPPX,JMAX+1,1,LUPRI)
         END IF
      END IF
      RETURN
      END
C  /* Deck erituv */
      SUBROUTINE ERITUV(HERINT,RJ000,PQXYZ,INDHER,IODDHR,
     &                  WORK,LWORK,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
      INTEGER T, U, V, TUV
      DIMENSION HERINT(NPPX,NTUV), RJ000(NPPX,0:JMAX), PQXYZ(NPPX,3),
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
         CALL DZERO(HERINT,NPPX*NTUV)
         DO 200 JVAL = 1, JMAX
            IF (MOD(JMAX-JVAL,2).EQ.0) THEN
               CALL ERIHRC(HERINT,WORK(KHRWRK),JVAL,RJ000,
     &                     PQXYZ(1,1),PQXYZ(1,2),PQXYZ(1,3),
     &                     INDHER,JMAX,MAXDER,NPPX,NTUV,
     &                     IPQX,IPQY,IPQZ)
            ELSE
               CALL ERIHRC(WORK(KHRWRK),HERINT,JVAL,RJ000,
     &                     PQXYZ(1,1),PQXYZ(1,2),PQXYZ(1,3),
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
#include <sdpre.h>
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
               TMIN1 = SDPRE(T - 1)
               TUV   = INDHER(T  ,0,0)
               M1T   = INDHER(T-1,0,0)
               M2T   = INDHER(T-2,0,0)
               DO 310 I = 1, NPPX
                  CUR(I,TUV) = PQX(I)*OLD(I,M1T) + TMIN1*OLD(I,M2T)
  310          CONTINUE
  300       CONTINUE
         ELSE
            DO 400 T = 2, MIN(MAXT,JVAL), 2
               TMIN1 = SDPRE(T - 1)
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
               UMIN1  = SDPRE(U - 1)
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
               UMIN1  = SDPRE(U - 1)
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
               VMIN1  = SDPRE(V - 1)
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
               VMIN1  = SDPRE(V - 1)
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
