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
C$Id: eri2car1.F,v 1.1.1.1 2001-02-08 13:33:13 hjj Exp $
C
#include <single.h>
C  /* Deck cr1drv */
      SUBROUTINE CR1DRV(HERINT,HCINT,INDHER,IODDHH,IODDHC,INDHSQ,
     &                  LMNPWR,IPNTUV,COOR12,COORP,PQINV,
     &                  WORK,LWORK,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
C
      LOGICAL DOREC, FRSTAO
      DIMENSION HERINT(*), HCINT(*),
     &          INDHER(*), INDHSQ(*), IODDHH(*), IODDHC(*),
     &          LMNPWR(*), IPNTUV(*),
     &          COOR12(*), COORP(*), PQINV(*),
     &          WORK(LWORK)
#include <cbieri.h>
#include <ericom.h>
#include <eriao.h>
      COMMON /DHCINF/ IHCADR(10), IHCSYM(10)
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from CR1DRV','*',103)
C
      IF (OLDCR1) THEN
C
C        Allocate CR1TWO
C
C        HCCONT  | HCSINT | LFCSNT | LFCINT
C                |
C                | ECOEF | HPI | P1 | P2
C                        |
C                        | EUV | ETUV | HCPRIM
C
         LCCONT = 0
         LCSINT = 0
         LFCSNT = 0
         LFCINT = 0
         LECOEF = 0
         LHPI   = 0
         LP1    = 0
         LP2    = 0
         LEUV   = 0
         LETUV  = 0
         LHCPRM = 0
         LFIRST = 0
C
         IF (KHKT12 .GT. 1) THEN
            IF (SPHR12) THEN
               LCCONT = NCCPP*NTUV34*KCKT12
               LFCINT = NTUV34*KHKT12
            END IF
            IF (SPHR1 .AND. SPHR2) THEN
               LCSINT = NCCPP*NTUV34*KHKT2
               LFCSNT = NTUV34*KHKT2
            END IF
C
            LECOEF = 3*NPPPP*(JMAX12 + 1)*(JMAX1 + 1)*(JMAX2+1)
            LHPI   =   NPPPP
            LP1    = 3*NPPPP
            LP2    = 3*NPPPP
C
            LEUV   = NPPPP
            LETUV  = NPPPP
C
            LHCPRM = NPPPP*NTUV34
            LFIRST = NTUV34
         END IF
C
         LTOTAL = LCCONT
     &          + MAX(
     &                 LCSINT + LFCSNT + LFCINT,
     &                 LECOEF + MAX(
     &                              LHPI+LP1+LP2,
     &                              LEUV+LETUV+LHCPRM+LFIRST
     &                             )
     &               )
         IF (LTOTAL.GT.LWORK) CALL STOPIT('CR1DRV','CR1TWO',KLAST,LWORK)
C
         KCCONT = 1
         KCSINT = KCCONT + LCCONT
         KFCSNT = KCSINT + LCSINT
         KFCINT = KFCSNT + LFCSNT
         KECOEF = KCCONT + LCCONT
         KHPI   = KECOEF + LECOEF
         KP1    = KHPI   + LHPI
         KP2    = KP1    + LP1
         KEUV   = KECOEF + LECOEF
         KETUV  = KEUV   + LEUV
         KHCPRM = KETUV  + LETUV
         KFIRST = KHCPRM + LHCPRM
C
         CALL CR1TWO(HERINT,INDHER,INDHER,IODDHH,IODDHC,INDHSQ,
     &               WORK(KECOEF),WORK(KEUV),WORK(KETUV),
     &               LMNPWR,IPNTUV,COOR12,COORP,PQINV,
     &               WORK(KHPI),WORK(KP1),WORK(KP2),
     &               HCINT,WORK(KFCINT),WORK(KHCPRM),WORK(KFIRST),
     &               WORK(KCCONT),WORK(KCSINT),WORK(KFCSNT),
     &               IPRINT)
      ELSE
C
C        Allocate CR1ONE
C
C        HCCONT  | XDIF
C                | HCSINT
C                |
C                | ECOEF | HPI | P1 | P2
C                        |
C                        | EUV | ETUV | HCPRIM
C
ckr         LCCONT = 0
ckr         LXDIF  = 0
ckr         LCSINT = 0
ckr         LECOEF = 0
ckr         LHPI   = 0
ckr         LP1    = 0
ckr         LP2    = 0
ckr         LEUV   = 0
ckr         LETUV  = 0
ckr         LHCPRM = 0
ckr         LFIRST = 0
C
ckr         DOREC  = MIN(NHKT1,NHKT2) .NE. 1
ckr         FRSTAO = NHKT1 .LE. NHKT2
C
ckr         IF (KHKT12 .GT. 1) THEN
ckr            IF (DOREC) THEN
ckr               LCCONT = NCCPP*NTUV34*KCREC1
ckr               LXDIF  = 3*NCCPP
ckr            END IF
C           IF (SPHR1 .AND. SPHR2) LCSINT = NCCPP*NTUV34*KHKT2
C
ckr            LECOEF = 3*NPPPP*(JMAX12 + 1)*(JMAX12 + 1)
ckr            LHPI   =   NPPPP
ckr            IF (FRSTAO) THEN
ckr               LP2 = 3*NPPPP
ckr            ELSE
ckr               LP1 = 3*NPPPP
ckr            END IF
C
ckr            LEUV   = NPPPP
ckr            LETUV  = NPPPP
C
ckr            LHCPRM = NPPPP*NTUV34
ckr            LFIRST = NTUV34
ckr         END IF
C
ckr         LTOTAL = LCCONT
ckr     &          + MAX( LXDIF,
ckr     &                 LCSINT,
ckr     &                 LECOEF + MAX(
ckr     &                              LHPI+LP1+LP2,
ckr     &                              LEUV+LETUV+LHCPRM+LFIRST
ckr     &                             )
ckr     &               )
ckr         IF (LTOTAL.GT.LWORK) CALL STOPIT('CR1DRV','CR1ONE',KLAST,LWORK)
C
ckr         KCCONT = 1
ckr         KCSINT = KCCONT + LCCONT
ckr         KXDIF  = KCCONT + LCCONT
ckr         KECOEF = KCCONT + LCCONT
ckr         KHPI   = KECOEF + LECOEF
ckr         KP1    = KHPI   + LHPI
ckr         KP2    = KP1    + LP1
ckr         KEUV   = KECOEF + LECOEF
ckr         KETUV  = KEUV   + LEUV
ckr         KHCPRM = KETUV  + LETUV
ckr         KFIRST = KHCPRM + LHCPRM
ckr         CALL CR1ONE(HERINT,INDHER,INDHER,IODDHH,IODDHC,INDHSQ,
ckr     &               WORK(KECOEF),WORK(KEUV),WORK(KETUV),
ckr     &               LMNPWR,IPNTUV,COOR12,COORP,PQINV,
ckr     &               WORK(KHPI),WORK(KP1),WORK(KP2),
ckr     &               HCINT,WORK(KFCINT),WORK(KHCPRM),WORK(KFIRST),
ckr     &               WORK(KCCONT),WORK(KXDIF),
ckr     &               WORK(KCSINT),WORK(KFCSNT),DOREC,FRSTAO,IPRINT)
      END IF
      RETURN
      END
C  /* Deck cr1two */
      SUBROUTINE CR1TWO(HERINT,INDHER,INDHVC,IODDHH,IODDHC,INDHSQ,
     &                  ECOEF,EUV,ETUV,LMNPWR,IPNTUV,
     &                  COOR12,COORP,PQINV,HPI,P1,P2,
     &                  HCINT,FCINT,HCPRIM,FRSTUV,HCCONT,
     &                  HCSINT,FCSINT,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
C
      INTEGER T, U, V, TUV
      LOGICAL SAMMRY, FRSTUV(NTUV34),
     &        FCINT(NTUV34,KHKT12), FCSINT(NTUV34,KHKT2)
C
      COMMON /DHCINF/ IHCADR(10), IHCSYM(10)
C
      DIMENSION HERINT(NPPPP,NTUV),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP), INDHVC(0:*),
     &          IODDHH(NRTOP), IODDHC(NRTOP), INDHSQ(NRTOP),
C
     &          HCINT(NCCPP,NTUV34,KHKT12),
C
     &          ECOEF(NPPPP,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2,3),
     &          ETUV(NPPPP), EUV(NPPPP),
     &          LMNPWR(KCKMAX,NHKMAX,3), IPNTUV(KC2MAX,2,2),
     &          COOR12(NPPPP,3,2), COORP(NPPPP,3), PQINV(NPPPP),
     &          HPI(NPPPP), P1(NPPPP,3), P2(NPPPP,3),
C
     &          HCPRIM(NPPPP,NTUV34),
     &          HCCONT(NCCPP,NTUV34,KCKT12),
     &          HCSINT(NCCPP,NTUV34,KHKT2)
C
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
#include <sphtrm.h>
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from CR1TWO','*',103)
C
C     COMMON /DHCINF/
C     ---------------
C
      IHCADR(1) = 1
      IHCSYM(1) = 0
C
C     *********************************
C     ***** Special Case: (ss|xy) *****
C     *********************************
C
      IF (KHKT12 .EQ. 1) THEN
         SAMMRY = .FALSE.
         IODD12 = 0
         DO 10 TUV = 1, NTUV34
            FRSTUV(TUV) = IODDHH(TUV) .NE. 0
   10    CONTINUE
         CALL CR1SEG(HERINT,HCINT,IODDHC,IODD12,SAMMRY,FRSTUV,IPRINT)
C
C     *********************************
C     ***** General Case: (xy|zw) *****
C     *********************************
C
      ELSE
C
C        Hermite-to-Cartesian expansion coefficients
C
         IELCTR = IPATH
         CALL EXPCFT(ECOEF,NPPPP,NPPPP,JMAX1,JMAX2,
     &               COOR12(1,1,1),COOR12(1,1,2),COORP,I120,PQINV,
     &               HPI,P1,P2,IELCTR,JMAX1,JMAX2,NCNT12,IPRINT)
C
C        Run over components
C        ===================
C
         ICMP12 = 0
C
         DO 200 ICOMP1 = 1, KCKT1
            L1 = LMNPWR(ICOMP1,NHKT1,1)
            M1 = LMNPWR(ICOMP1,NHKT1,2)
            N1 = LMNPWR(ICOMP1,NHKT1,3)
            MAX2 = KCKT2
            IF (TCMP12) MAX2 = ICOMP1
            DO 300 ICOMP2 = 1, MAX2
               L2 = LMNPWR(ICOMP2,NHKT2,1)
               M2 = LMNPWR(ICOMP2,NHKT2,2)
               N2 = LMNPWR(ICOMP2,NHKT2,3)
C
               ICMP12 = ICMP12 + 1
               L = L1 + L2
               M = M1 + M2
               N = N1 + N2
C
C              Primitive integrals
C              ===================
C
               INCT = I120(1) + 1
               INCU = I120(2) + 1
               INCV = I120(3) + 1
               MINT = IBTAND(L,I120(1))
               MINU = IBTAND(M,I120(2))
               MINV = IBTAND(N,I120(3))
               IF (IPRINT .GT. 25) THEN
                  WRITE(LUPRI,'(/,1X,A,3I5/)')' ICOMP1, ICOMP2, ICMP12',
     &                                          ICOMP1, ICOMP2, ICMP12
                  WRITE(LUPRI,'(1X,A,15X,3I5)')' T loop:',MINT,L,INCT
                  WRITE(LUPRI,'(1X,A,15X,3I5)')' U loop:',MINU,M,INCU
                  WRITE(LUPRI,'(1X,A,15X,3I5)')' V loop:',MINV,N,INCV
               END IF
C
               DO 350 TUV = 1, NTUV34
                  FRSTUV(TUV) = .TRUE.
  350          CONTINUE
C
               DO 400 V = MINV, N, INCV
               DO 400 U = MINU, M, INCU
                  DO 450 I = 1, NPPPP
                     EUV(I) = ECOEF(I,V,N1,N2,3)
     &                      * ECOEF(I,U,M1,M2,2)
  450             CONTINUE
               DO 400 T = MINT, L, INCT
C
                  DO 500 I = 1, NPPPP
                     ETUV(I) = ECOEF(I,T,L1,L2,1)*EUV(I)
  500             CONTINUE
C
                  ITUV = INDHER(T,U,V)
                  INDS = INDHSQ(ITUV)
                  DO 600 TUV = 1, NTUV34
                  IF (IODDHH(ITUV) .EQ. IODDHH(TUV)) THEN
#if defined (SYS_AIX)
C                    code due to AIX xlf version 2.2 bug
                     INDT = INDS + INDHSQ(TUV)
                     INDT = INDHVC(INDT)
#else
                     INDT = INDHVC(INDS + INDHSQ(TUV))
#endif
                     IF (FRSTUV(TUV)) THEN
                        FRSTUV(TUV) = .FALSE.
                        DO 700 I = 1, NPPPP
                           HCPRIM(I,TUV) = ETUV(I)*HERINT(I,INDT)
  700                   CONTINUE
                     ELSE
                        DO 750 I = 1, NPPPP
                           HCPRIM(I,TUV) = HCPRIM(I,TUV)
     &                          + ETUV(I)*HERINT(I,INDT)
  750                   CONTINUE
                     END IF
                  END IF
C
  600             CONTINUE
  400          CONTINUE
C
C              Contracted integrals
C              ====================
C
               SAMMRY = .FALSE.
               IODD12 = IODDHC(INDHER(L,M,N))
               IF (SPHR12) THEN
                  CALL CR1SEG(HCPRIM,HCCONT(1,1,ICMP12),IODDHC,IODD12,
     &                        SAMMRY,FRSTUV,IPRINT)
               ELSE
                  CALL CR1SEG(HCPRIM,HCINT(1,1,ICMP12),IODDHC,IODD12,
     &                        SAMMRY,FRSTUV,IPRINT)
               END IF
C
  300       CONTINUE
  200    CONTINUE
C
C        Spherical integrals
C        ===================
C
         IF (SPHR12) THEN
            CALL CR1SPH(HCCONT,HCINT,FCINT,HCSINT,FCSINT,
     &                  CSP(ISPADR(NHKT1)),CSP(ISPADR(NHKT2)),
     &                  IODDHC,IPNTUV,IPRINT)
         END IF
      END IF
      IF (IPRINT .GT. 20) THEN
         CALL HCPRINT(HCINT,IODDHC,IPNTUV,IPRINT)
      END IF
      RETURN
      END
C  /* Deck cr1one */
      SUBROUTINE CR1ONE(HERINT,INDHER,INDHVC,IODDHH,IODDHC,INDHSQ,
     &                  ECOEF,EUV,ETUV,LMNPWR,IPNTUV,
     &                  COOR12,COORP,PQINV,HPI,P1,P2,
     &                  HCINT,FCINT,HCPRIM,FRSTUV,HCCONT,XDIF,
     &                  HCSINT,FCSINT,DOREC,FRSTAO,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
C
      INTEGER T, U, V, TUV
      LOGICAL FRSTAO, DOREC, SAMMRY, FRSTUV(NTUV34),
     &        FCINT(NTUV34,KHKT12), FCSINT(NTUV34,KHKT2)
C
C
      COMMON /DHCINF/ IHCADR(10), IHCSYM(10)
C
      DIMENSION HERINT(NPPPP,NTUV),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP), INDHVC(0:*),
     &          IODDHH(NRTOP), IODDHC(NRTOP), INDHSQ(NRTOP),
C
     &          HCINT(NCCPP,NTUV34,KHKT12),
C
     &          ECOEF(NPPPP,0:JMAX1+JMAX2,0:JMAX1+JMAX2,3),
     &          ETUV(NPPPP), EUV(NPPPP),
     &          LMNPWR(KCKMAX,NHKMAX,3), IPNTUV(KC2MAX,2,2),
     &          COOR12(NPPPP,3,2), COORP(NPPPP,3), PQINV(NPPPP),
     &          HPI(NPPPP), P1(NPPPP,3), P2(NPPPP,3),
C
     &          HCPRIM(NPPPP,NTUV34),
     &          HCCONT(NCCPP,NTUV34,KCREC1),
     &          XDIF(NCCPP,3),
     &          HCSINT(*)
C
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
#include <sphtrm.h>
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from CR1ONE','*',103)
C
C
C     COMMON /DHCINF/
C     ---------------
C
      IHCADR(1) = 1
      IHCSYM(1) = 0
C
C     *********************************
C     ***** Special Case: (ss|xy) *****
C     *********************************
C
      IF (KHKT12 .EQ. 1) THEN
         SAMMRY = .FALSE.
         IODD12 = 0
         DO 100 TUV = 1, NTUV34
            FRSTUV(TUV) = .FALSE.
  100    CONTINUE
         CALL CR1SEG(HERINT,HCINT,IODDHC,IODD12,SAMMRY,FRSTUV,IPRINT)
C
C     *********************************
C     ***** General Case: (xy|zw) *****
C     *********************************
C
      ELSE
C
C        Hermite-to-Cartesian expansion coefficients
C        ===========================================
C
         IELCTR = IPATH
C
         IF (FRSTAO) THEN
            JMX1 = 0
            JMX2 = JMAX1 + JMAX2
         ELSE
            JMX1 = JMAX1 + JMAX2
            JMX2 = 0
         END IF
         CALL EXPCFT(ECOEF,NPPPP,NPPPP,JMX1,JMX2,
     &               COOR12(1,1,1),COOR12(1,1,2),COORP,I120,PQINV,
     &               HPI,P1,P2,IELCTR,JMX1,JMX2,NCNT12,IPRINT)
C
C        Run over components
C        ===================
C
         ICMP12 = 0
         DO 200 J = MAX(JMAX1,JMAX2), JMAX1 + JMAX2
C
            DO 300 L = J,     0, -1
            DO 300 M = J - L, 0, -1
               N = J - L - M
C
               ICMP12 = ICMP12 + 1
C
C              Primitive integrals
C              ===================
C
               INCT = I120(1) + 1
               INCU = I120(2) + 1
               INCV = I120(3) + 1
               MINT = IBTAND(L,I120(1))
               MINU = IBTAND(M,I120(2))
               MINV = IBTAND(N,I120(3))
               DO 350 TUV = 1, NTUV34
                  FRSTUV(TUV) = .TRUE.
  350          CONTINUE
C
               DO 400 V = MINV, N, INCV
               DO 400 U = MINU, M, INCU
                  DO 450 I = 1, NPPPP
                     EUV(I) = ECOEF(I,V,N,3)
     &                      * ECOEF(I,U,M,2)
  450             CONTINUE
               DO 400 T = MINT, L, INCT
C
                  DO 500 I = 1, NPPPP
                     ETUV(I) = ECOEF(I,T,L,1)*EUV(I)
  500             CONTINUE
C
                  ITUV = INDHER(T,U,V)
                  INDS = INDHSQ(ITUV)
                  DO 600 TUV = 1, NTUV34
                  IF (IODDHH(ITUV) .EQ. IODDHH(TUV)) THEN
#if defined (SYS_AIX)
C                    code due to AIX xlf version 2.2 bug
                     INDT = INDS + INDHSQ(TUV)
                     INDT = INDHVC(INDT)
#else
                     INDT = INDHVC(INDS + INDHSQ(TUV))
#endif
                     IF (FRSTUV(TUV)) THEN
                        FRSTUV(TUV) = .FALSE.
                        DO 700 I = 1, NPPPP
                           HCPRIM(I,TUV) = ETUV(I)*HERINT(I,INDT)
  700                   CONTINUE
                     ELSE
                        DO 750 I = 1, NPPPP
                           HCPRIM(I,TUV) = HCPRIM(I,TUV)
     &                          + ETUV(I)*HERINT(I,INDT)
  750                   CONTINUE
                     END IF
                  END IF
C
  600             CONTINUE
  400          CONTINUE
C
C              Contracted integrals
C              ====================
C
               SAMMRY = .FALSE.
               IODD12 = IODDHC(INDHER(L,M,N))
               IF (DOREC) THEN
                  CALL CR1SEG(HCPRIM,HCCONT(1,1,ICMP12),IODDHC,IODD12,
     &                        SAMMRY,FRSTUV,IPRINT)
               ELSE
                  CALL CR1SEG(HCPRIM,HCINT(1,1,ICMP12),IODDHC,IODD12,
     &                        SAMMRY,FRSTUV,IPRINT)
               END IF
C
  300       CONTINUE
  200    CONTINUE
C
C        Orbital transfer recursion
C        ==========================
C
         IF (DOREC) THEN
            CALL CR1OTR(HCCONT,HCINT,INDHER,IODDHC,COOR12,XDIF,
     &                  FRSTAO,IPRINT)
         END IF
C
C        Spherical integrals
C        ===================
C
         IF (SPHR12) THEN
            CALL CR1SPH(HCCONT,HCINT,FCINT,HCSINT,FCSINT,
     &                  CSP(ISPADR(NHKT1)),CSP(ISPADR(NHKT2)),
     &                  IODDHC,IPNTUV,NTUV34,IPRINT)
         END IF
      END IF
      IF (IPRINT .GT. 20) THEN
         CALL HCPRINT(HCINT,IODDHC,IPNTUV,IPRINT)
      END IF
      RETURN
      END
C  /* Deck cr1otr */
      SUBROUTINE CR1OTR(OINT,TINT,INDHER,IODDHC,COOR12,XDIF,FRSTAO,
     &                  IPRINT)
#include <implicit.h>
#include <priunit.h>
      LOGICAL FRSTAO
      INTEGER X
      DIMENSION OINT(NCCPP,NTUV34,*),
     &          TINT(NCCPP,NTUV34,*),
     &          IODDHC(NRTOP), INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          COOR12(NPPPP,3,2), XDIF(NCCPP,3)
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
C
      ITRI1(I) = (I + 1)*(I + 2)/2
      ITRI2(I) = (I + 1)*(I + 2)*(I + 3)/6
C
      IF (FRSTAO) THEN
         IA = 2
         IB = 1
      ELSE
         IA = 1
         IB = 2
      END IF
C
C     XDIF
C     ====
C
      DO 100 X = 1, 3
      IF (I120(X) .EQ. 0) THEN
         DO 200 I = 1, NCCPP
            XDIF(I,X) = COOR12(I,X,IA) - COOR12(I,X,IB)
  200    CONTINUE
      END IF
  100 CONTINUE
C
C     Recursion
C     =========
C
      NRCURS = MIN(JMAX1,JMAX2)
      N2OFF = ITRI2(MAX(JMAX1,JMAX2) - 1)
      NC1 = ITRI1(0)
      NC2 = ITRI2(JMAX1 + JMAX2) - N2OFF
      IF (FRSTAO) THEN
         NO = NC2
      ELSE
         NO = NC1
      END IF
      NDIM1 = NCCPP*NTUV34*NC1*NC2
      IF (MOD(NRCURS,2).EQ.0) CALL DCOPY(NDIM,OINT,1,TINT,1)
C
      DO 300 JVAL = 1, NRCURS
         NC1 = ITRI1(JVAL)
         NC2 = ITRI2(JMAX1 + JMAX2 - JVAL) - N2OFF
         IF (FRSTAO) THEN
            NC = NC2
         ELSE
            NC = NC1
         END IF
         IF (MOD(NRCURS - JVAL,2).EQ.0) THEN
            CALL CR1RCR(JVAL,OINT,TINT,NC,NO,XDIF,
     &                  INDHER,IODDHC,FRSTAO,IPRINT)
         ELSE
            CALL CR1RCR(JVAL,TINT,OINT,NC,NO,XDIF,
     &                  INDHER,IODDHC,FRSTAO,IPRINT)
         END IF
         NO = NC
  300 CONTINUE
C
      IF (TCMP12) THEN
         IJSQR = 0
         IJTRI = 0
         DO 400 I = 1, KCKT1
         DO 400 J = 1, KCKT1
            IJSQR = IJSQR + 1
            IF (I .GE. J) THEN
               IJTRI = IJTRI + 1
               IF (IJSQR .NE. IJTRI) THEN
                  CALL DCOPY(NCCPP*NTUV34,TINT(1,1,IJSQR),1,
     &                                    TINT(1,1,IJTRI),1)
               END IF
            END IF
  400    CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck cr1rcr */
      SUBROUTINE CR1RCR(JVAL,OLD,CUR,NC,NO,XDIF,
     &                  INDHER,IODDHC,FRSTAO,IPRINT)
#include <implicit.h>
#include <priunit.h>
C
      INTEGER X, Y, Z, TUV, P
      LOGICAL NOFAC, FRSTAO
      DIMENSION CUR(NCCPP,NTUV34,*), OLD(NCCPP,NTUV34,*),
     &          XDIF(NCCPP,3), IODDHC(NRTOP),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP)
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
#include <ibtfun.h>
C
C     Run over components for first electron
C     ======================================
C
      DO 100 I = JVAL,     0, -1
      DO 100 J = JVAL - I, 0, -1
         K = JVAL - I - J
C
         X = 0
         Y = 0
         Z = 0
         NOFAC = .FALSE.
         IF      (I120(1).EQ.1 .AND. I.GT.0) THEN
            ICOOR = 1
            X     = 1
            NOFAC = .TRUE.
         ELSE IF (I120(2).EQ.1 .AND. J.GT.0) THEN
            ICOOR = 2
            Y     = 1
            NOFAC = .TRUE.
         ELSE IF (I120(3).EQ.1 .AND. K.GT.0) THEN
            ICOOR = 3
            Z     = 1
            NOFAC = .TRUE.
         ELSE IF (I .GT. 0) THEN
            ICOOR = 1
            X     = 1
         ELSE IF (J .GT. 0) THEN
            ICOOR = 2
            Y = 1
         ELSE
            ICOOR = 3
            Z = 1
         END IF
C
         IJK    = INDHER(I,J,K)       - INDHER(JVAL,  0,0) + 1
         IJKD   = INDHER(I-X,J-Y,K-Z) - INDHER(JVAL-1,0,0) + 1
C
C        Run over components for second orbital
C        ======================================
C
         JSTART = MAX(JMAX1,JMAX2)
         DO 200 JVAL2 = JSTART, JMAX1 + JMAX2 - JVAL
         DO 200 L = JVAL2,     0, -1
         DO 200 M = JVAL2 - L, 0, -1
            N = JVAL2 - L - M
            LMN  = INDHER(L,M,N)       - INDHER(JSTART,0,0) + 1
            LMNI = INDHER(L+X,M+Y,N+Z) - INDHER(JSTART,0,0) + 1
C
            IF (FRSTAO) THEN
               ICUR  = (IJK  - 1)*NC + LMN
               IOLD1 = (IJKD - 1)*NO + LMN
               IOLD2 = (IJKD - 1)*NO + LMNI
            ELSE
               ICUR  = (LMN  - 1)*NC + IJK
               IOLD1 = (LMN  - 1)*NO + IJKD
               IOLD2 = (LMNI - 1)*NO + IJKD
            END IF
C
C           PQ symmetry
C           ===========
C
            IF (IHCXYZ.GT.0) THEN
               IODD = IODDHC(INDHER(I+L,J+M,K+N))
               IF (NOFAC) THEN
                  DO 300 TUV = 1, NTUV34
                  IF (IODDHC(TUV) .EQ. IODD) THEN
                     DO 310 P = 1, NCCPP
                       CUR(P,TUV,ICUR) = OLD(P,TUV,IOLD2)
  310                CONTINUE
                  END IF
  300             CONTINUE
               ELSE
                  DO 400 TUV = 1, NTUV34
                  IF (IODDHC(TUV) .EQ. IODD) THEN
                     DO 410 P = 1, NCCPP
                        CUR(P,TUV,ICUR) = XDIF(P,ICOOR)*OLD(P,TUV,IOLD1)
     &                                                + OLD(P,TUV,IOLD2)
  410                CONTINUE
                  END IF
  400             CONTINUE
               END IF
C
C           No PQ symmetry
C           ==============
C
            ELSE
               IF (NOFAC) THEN
                  CALL DCOPY(NCCPP*NTUV34,OLD(1,1,ICUR),1,
     &                                    CUR(1,1,IOLD2),1)
               ELSE
                  DO 500 TUV = 1, NTUV34
                  DO 500 P = 1, NCCPP
                     CUR(P,TUV,ICUR) = XDIF(P,ICOOR)*OLD(P,TUV,IOLD1)
     &                                             + OLD(P,TUV,IOLD2)
  500             CONTINUE
               END IF
            END IF
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck hcprint */
      SUBROUTINE HCPRINT(HCINT,IODDHC,IPNTUV,IPRINT)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (NHCTYP = 1)
      INTEGER TUV
      DIMENSION HCINT(NPQBCX,NPRF34,NCTF12,NTUV34,KHKT12,NHCTYP),
     &          IODDHC(NRTOP), IPNTUV(KC2MAX,2,2)
C
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
      COMMON /DHCINF/ IHCADR(10), IHCSYM(10)
#include <ibtfun.h>
C
C     Print Hermite-Spherical integrals
C     =================================
C
      CALL HEADER('Hermite-spherical integrals in HCPRINT',-1)
C
      WRITE (LUPRI,'(7X,A, I5)')'NPQBCX            ',NPQBCX
      WRITE (LUPRI,'(7X,A,2I5)')'NPRF34,NCTF12:    ',NPRF34,NCTF12
      WRITE (LUPRI,'(7X,A, I5)')'NTUV34            ',NTUV34
      WRITE (LUPRI,'(7X,A,3I5)')'KHKT1,KHKT2,KHKT12',KHKT1,KHKT2,KHKT12
      WRITE (LUPRI,'(7X,A, I5)')'Num. of integrals:',
     &                       NPQBCX*NPRF34*NCTF12*NTUV34*KHKT12*NHCTYP
C
      DO 100 ITYPE = 1, NHCTYP
         ICMP12 = 0
         DO 200 ICOMP1 = 1, KHKT1
            MAX2 = KHKT2
            IF (TKMP12) MAX2 = ICOMP1
            DO 300 ICOMP2 = 1, MAX2
               ICMP12 = ICMP12 + 1
               DO 400 TUV = 1, NTUV34
               IF (IODDHC(TUV) .EQ. IODDHC(IPNTUV(ICMP12,2,1))) THEN
                  WRITE (LUPRI,'(/,1X,A,I3,1X,A,2I3,1X,A,I3)' )
     &               'Integral type:', ITYPE,
     &               'Components:   ', ICOMP1,ICOMP2,
     &               'TUV:          ', TUV
                  CALL OUTPUT(HCINT(1,1,1,TUV,ICMP12,ITYPE),
     &                        1,NPQBCX,1,NPRF34*NCTF12,
     &                        NPQBCX,NPRF34*NCTF12,
     &                        1,LUPRI)
               END IF
  400          CONTINUE
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck cr1seg */
      SUBROUTINE CR1SEG(PRMINT,CNTINT,IODDHC,IODD12,SAMMRY,FRSTUV,
     &                  IPRINT)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0)
      INTEGER TUV
      LOGICAL SAMMRY, FRSTUV(NTUV34)
      DIMENSION PRMINT(NCCPP,NPRF12,NTUV34),
     &          CNTINT(NCCPP,NTUV34), IODDHC(NRTOP)
#include <hertop.h>
#include <ericom.h>
#include <eriao.h>
C
      IF (NPRF12.EQ.1 .AND. .NOT.SAMMRY) THEN
         DO 100 TUV = 1, NTUV34
            IF (IODDHC(TUV) .EQ. IODD12) THEN
               IF (FRSTUV(TUV)) THEN
                  DO 200 J = 1, NCCPP
                     CNTINT(J,TUV) = D0
  200             CONTINUE
               ELSE
                  DO 300 J = 1, NCCPP
                     CNTINT(J,TUV) = PRMINT(J,1,TUV)
  300             CONTINUE
               END IF
            END IF
  100    CONTINUE
      ELSE
         DO 400 TUV = 1, NTUV34
         IF (IODDHC(TUV) .EQ. IODD12) THEN
            IF (FRSTUV(TUV)) THEN
               DO 500 J = 1, NCCPP
                  CNTINT(J,TUV) = D0
  500          CONTINUE
            ELSE
               DO 600 J = 1, NCCPP
                  CNTINT(J,TUV) = PRMINT(J,1,TUV) + PRMINT(J,2,TUV)
  600          CONTINUE
               DO 700 I = 3, NPRF12
               DO 700 J = 1, NCCPP
                  CNTINT(J,TUV) = CNTINT(J,TUV) + PRMINT(J,I,TUV)
  700          CONTINUE
            END IF
         END IF
  400    CONTINUE
      END IF
C
      IF (IPRINT .GT. 25) THEN
         CALL HEADER('Output from CR1SEG',-1)
         WRITE (LUPRI,'(2X,A,2I5)') 'NCCPP, NPRF12', NCCPP, NPRF12
         DO 800 TUV = 1, NTUV34
            IF (IODDHC(TUV) .EQ. IODD12) THEN
               CALL HEADER('PRMINT in CR1SEG',-1)
               WRITE (LUPRI,'(2X,A,I5)') 'TUV: ',TUV
               CALL OUTPUT(PRMINT(1,1,TUV),1,NCCPP,1,NPRF12,
     &                     NCCPP,NPRF12,1,LUPRI)
               CALL HEADER('CNTINT in CR1SEG',-1)
               WRITE (LUPRI,'(2X,A,I5)') 'TUV: ',TUV
               CALL OUTPUT(CNTINT(1,TUV),1,1,1,NCCPP,1,NCCPP,1,LUPRI)
            END IF
  800    CONTINUE
      END IF
C
      RETURN
      END
