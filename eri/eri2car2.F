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
C$Id: eri2car2.F,v 1.1.1.1 2001-02-08 13:33:13 hjj Exp $
C
#include <single.h>
C  /* Deck cr2drv */
      SUBROUTINE CR2DRV(AOINT,HCINT,INDHER,IODDHC,IODDCC,
     &                  LMNPWR,IPNTUV,COOR34,COORP,PQINV,
     &                  WORK,LWORK,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      DIMENSION AOINT(NCCCC,KHKT12,KHKT34),
     &          HCINT(NCCPP,NTUV34,KHKT12),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          IODDHC(NRTOP), IODDCC(NRTOP),
     &          LMNPWR(KCKMAX,NHKMAX,3), IPNTUV(KC2MAX,2,2),
     &          COOR34(NPPPP,3,2), COORP(NPPPP,3), PQINV(NPPPP),
     &          WORK(LWORK)
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from CR2DRV','*',103)
C
C        Work space allocation in CR2TWO
C
C        CCONT  | CSINT  | FCSINT | FAOINT
C               |
C               | ECOEF  | HPI    | P3     | P4
C                        | EUV    | ETUV   | CCPRIM
C
         LCCONT = 0
         LCSINT = 0
         LFCSNT = 0
         LFAONT = 0
         LECOEF = 0
         LHPI   = 0
         LP3    = 0
         LP4    = 0
         LEUV   = 0
         LETUV  = 0
         LCCPRM = 0
         LFIRST = 0
C
         IF (KHKT34 .GT. 1) THEN
            IF (SPHR34) THEN
               LCCONT = NCCCC*KHKT12*KCKT34
               LFAONT = KHKT12*KHKT34
            END IF
            IF (SPHR3 .AND. SPHR4) THEN
               LCSINT = NCCCC*KHKT12*KHKT4
               LFCSNT = KHKT12*KHKT4
            END IF
C
            LECOEF = 3*NCCPP*(JMAX34 + 1)*(JMAX3 + 1)*(JMAX4 + 1)
            LHPI   =   NCCPP
            LP3    = 3*NCCPP
            LP4    = 3*NCCPP
C
            LEUV   = NCCPP
            LETUV  = NCCPP
C
            LCCPRM = NCCPP*KHKT12
            LFIRST = KHKT12
         END IF
C
         LTOTAL = LCCONT + MAX(
     &                         LCSINT + LFCSNT + LFAONT,
     &                         LECOEF + MAX(
     &                                      LHPI + LP3 + LP4,
     &                                      LEUV+LETUV+LCCPRM+LFIRST
     &                                     )
     &                        )
         IF (LTOTAL.GT.LWORK) CALL STOPIT('CR2DRV',' ',LTOTAL,LWORK)
C
         KCCONT = 1
         KCSINT = KCCONT + LCCONT
         KFCSNT = KCSINT + LCSINT
         KFAONT = KFCSNT + LFCSNT
         KECOEF = KCCONT + LCCONT
         KHPI   = KECOEF + LECOEF
         KP3    = KHPI   + LHPI
         KP4    = KP3    + LP3
         KEUV   = KECOEF + LECOEF
         KETUV  = KEUV   + LEUV
         KCCPRM = KETUV  + LETUV
         KFIRST = KCCPRM + LCCPRM
C
         CALL CR2TWO(AOINT,WORK(KFAONT),
     &               HCINT,INDHER,IODDHC,IODDCC,WORK(KCCPRM),
     &               WORK(KFIRST),WORK(KCCONT),
     &               WORK(KCSINT),WORK(KFCSNT),
     &               WORK(KECOEF),WORK(KEUV),WORK(KETUV),LMNPWR,IPNTUV,
     &               COOR34,COORP,PQINV,WORK(KHPI),WORK(KP3),WORK(KP4),
     &               IPRINT)
C
      RETURN
      END
C  /* Deck cr2two */
      SUBROUTINE CR2TWO(AOINT,FAOINT,
     &                  HCINT,INDHER,IODDHC,IODDCC,
     &                  CCPRIM,FRSTUV,CCONT,CSINT,FCSINT,
     &                  ECOEF,EUV,ETUV,LMNPWR,IPNTUV,
     &                  COOR34,COORP,PQINV,HPI,P3,P4,
     &                  IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      INTEGER T, U, V
      LOGICAL SAMMRY, FRSTUV(KHKT12),
     &        FAOINT(KHKT12,KHKT34), FCSINT(KHKT12,KHKT4)
      DIMENSION AOINT(NCCCC,KHKT12,KHKT34),
     &          HCINT(NCCPP,NTUV34,KHKT12),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          IODDHC(NRTOP), IODDCC(NRTOP),
     &          CCPRIM(NCCPP,KHKT12),
     &          CCONT(NCCCC,KHKT12,KCKT34),
     &          CSINT(NCCCC,KHKT12,KHKT4),
     &          ECOEF(NCCPP,0:JMAX3+JMAX4,0:JMAX3,0:JMAX4,3),
     &          EUV(NCCPP), ETUV(NCCPP),
     &          LMNPWR(KCKMAX,NHKMAX,3), IPNTUV(KC2MAX,2,2),
     &          COOR34(NPPPP,3,2), COORP(NPPPP,3), PQINV(NPPPP),
     &          HPI(NCCPP), P3(NCCPP,3), P4(NCCPP,3)
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
#include <sphtrm.h>
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from CR2TWO','*',103)
C
      IF (KHKT34 .EQ. 1) THEN
         SAMMRY = .FALSE.
         IODD34 = 0
         DO 10 ICMP12 = 1, KHKT12
            FRSTUV(ICMP12) = IODDHC(IPNTUV(ICMP12,2,1)) .NE. 0
   10    CONTINUE
         CALL CR2SEG(HCINT,AOINT,IODDCC,IPNTUV,IODD34,SAMMRY,FRSTUV,
     &               IPRINT)
      ELSE
C
C        Hermite-to-Cartesian expansion coefficients
C
         IELCTR = 2
         IF (IPATH .EQ. 2) IELCTR = 1
         CALL EXPCFT(ECOEF,NCCPP,NPPPP,JMAX3,JMAX4,
     &               COOR34(1,1,1),COOR34(1,1,2),COORP,I340,PQINV,
     &               HPI,P3,P4,IELCTR,JMAX3,JMAX4,NCNT34,IPRINT)
C
C        Run over components
C
         ICMP34 = 0
         DO 100 ICOMP3 = 1, KCKT3
            L3 = LMNPWR(ICOMP3,NHKT3,1)
            M3 = LMNPWR(ICOMP3,NHKT3,2)
            N3 = LMNPWR(ICOMP3,NHKT3,3)
C
            MAX4 = KCKT4
            IF (TCMP34) MAX4 = ICOMP3
            DO 200 ICOMP4 = 1, MAX4
               ICMP34 = ICMP34 + 1
               L4 = LMNPWR(ICOMP4,NHKT4,1)
               M4 = LMNPWR(ICOMP4,NHKT4,2)
               N4 = LMNPWR(ICOMP4,NHKT4,3)
C
               INCT = I340(1) + 1
               INCU = I340(2) + 1
               INCV = I340(3) + 1
               MAXT = L3 + L4
               MAXU = M3 + M4
               MAXV = N3 + N4
               MINT = IBTAND(MAXT,I340(1))
               MINU = IBTAND(MAXU,I340(2))
               MINV = IBTAND(MAXV,I340(3))
C
               IF (IPRINT .GT. 25) THEN
                  WRITE(LUPRI,'(/,1X,A,3I5/)')' ICOMP3, ICOMP4, ICMP34',
     &                                          ICOMP3, ICOMP4, ICMP34
                  WRITE(LUPRI,'(1X,A,15X,3I5)')' T loop:',MINT,MAXT,INCT
                  WRITE(LUPRI,'(1X,A,15X,3I5)')' U loop:',MINU,MAXU,INCU
                  WRITE(LUPRI,'(1X,A,15X,3I5)')' V loop:',MINV,MAXV,INCV
               END IF
C
               DO 250 ICMP12 = 1, KHKT12
                  FRSTUV(ICMP12) = .TRUE.
  250          CONTINUE
C
               DO 300 V = MINV, MAXV, INCV
               DO 300 U = MINU, MAXU, INCU
                  DO 350 I = 1, NCCPP
                     EUV(I) = ECOEF(I,V,N3,N4,3)
     &                      * ECOEF(I,U,M3,M4,2)
  350             CONTINUE
               DO 300 T = MINT, MAXT, INCT
C
C                 Expansion coefficients
C                 ======================
C
                     DO 310 I = 1, NCCPP
                        ETUV(I) = ECOEF(I,T,L3,L4,1)*EUV(I)
  310                CONTINUE
C
C                    Cartesian integrals
C                    ===================
C
                     ITUV = INDHER(T,U,V)
                     IODD = IODDHC(ITUV)
                     DO 400 ICMP12 = 1, KHKT12
                     IF (IODDHC(IPNTUV(ICMP12,2,1)) .EQ. IODD) THEN
                        IF (FRSTUV(ICMP12)) THEN
                           FRSTUV(ICMP12) = .FALSE.
                           DO 410 I = 1, NCCPP
                              CCPRIM(I,ICMP12) =
     &                               ETUV(I)*HCINT(I,ITUV,ICMP12)
  410                      CONTINUE
                        ELSE
                           DO 420 I = 1, NCCPP
                              CCPRIM(I,ICMP12) = CCPRIM(I,ICMP12)
     &                             + ETUV(I)*HCINT(I,ITUV,ICMP12)
  420                      CONTINUE
                        END IF
                     END IF
  400                CONTINUE
  300          CONTINUE
C
C              ********************************
C              ***** Contracted Integrals *****
C              ********************************
C
               SAMMRY = .FALSE.
               IODD34 = IODDCC(IPNTUV(ICMP34,1,2))
               IF (SPHR34) THEN
                  CALL CR2SEG(CCPRIM,CCONT(1,1,ICMP34),IODDCC,IPNTUV,
     &                        IODD34,SAMMRY,FRSTUV,IPRINT)
               ELSE
                  CALL CR2SEG(CCPRIM,AOINT(1,1,ICMP34),IODDCC,IPNTUV,
     &                        IODD34,SAMMRY,FRSTUV,IPRINT)
               END IF
  200       CONTINUE
  100    CONTINUE
C
C        Spherical integrals
C        ===================
C
         IF (SPHR34) THEN
            CALL CR2SPH(CCONT,CSINT,FCSINT,AOINT,FAOINT,
     &                  CSP(ISPADR(NHKT3)),CSP(ISPADR(NHKT4)),
     &                  IODDCC,IPNTUV,IPRINT)
         END IF
      END IF
C
C     *************************
C     ***** Print Section *****
C     *************************
C
      IF (IPRINT .GE. 15) THEN
         CALL HEADER('Final spherical integrals - CR2TWO',-1)
         DO 900 I = 1, KHKT12
         DO 900 J = 1, KHKT34
         IF (IODDCC(IPNTUV(I,2,1)) .EQ. IODDCC(IPNTUV(J,2,2))) THEN
            WRITE (LUPRI, '(/1X,A,2I3)') ' ICMP12, ICMP34 ', I, J
            CALL OUTPUT(AOINT(1,I,J),1,1,1,NCCCC,1,NCCCC,1,LUPRI)
         END IF
  900    CONTINUE
      END IF
      RETURN
      END
C  /* Deck cr2seg */
      SUBROUTINE CR2SEG(PRMINT,CNTINT,IODDCC,IPNTUV,IODD34,SAMMRY,
     &                  FRSTUV,IPRINT)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0)
      LOGICAL SAMMRY, FRSTUV(KHKT12)
      DIMENSION PRMINT(NCCCC,NPRF34,KHKT12),
     &          CNTINT(NCCCC,KHKT12),
     &          IPNTUV(KC2MAX,2,2), IODDCC(NRTOP)
#include <hertop.h>
#include <ericom.h>
#include <eriao.h>
C
      IF (NPRF34.EQ.1 .AND. .NOT.SAMMRY) THEN
         DO 100 ICMP12 = 1, KHKT12
         IF (IODDCC(IPNTUV(ICMP12,2,1)) .EQ. IODD34) THEN
            IF (FRSTUV(ICMP12)) THEN
               DO 200 J = 1, NCCCC
                  CNTINT(J,ICMP12) = D0
  200          CONTINUE
            ELSE
               DO 300 J = 1, NCCCC
                  CNTINT(J,ICMP12) = PRMINT(J,1,ICMP12)
  300          CONTINUE
            END IF
         END IF
  100    CONTINUE
      ELSE
         DO 400 ICMP12 = 1, KHKT12
         IF (IODDCC(IPNTUV(ICMP12,2,1)) .EQ. IODD34) THEN
            IF (FRSTUV(ICMP12)) THEN
               DO 500 J = 1, NCCCC
                  CNTINT(J,ICMP12) = D0
  500          CONTINUE
            ELSE
               DO 600 J = 1, NCCCC
                  CNTINT(J,ICMP12) = PRMINT(J,1,ICMP12)
     &                             + PRMINT(J,2,ICMP12)
  600          CONTINUE
               DO 700 I = 3, NPRF34
               DO 700 J = 1, NCCCC
                  CNTINT(J,ICMP12) = CNTINT(J,ICMP12)
     &                             + PRMINT(J,I,ICMP12)
  700          CONTINUE
            END IF
         END IF
  400    CONTINUE
      END IF
C
      IF (IPRINT .GT. 25) THEN
         CALL HEADER('Output from CR2SEG',-1)
         WRITE (LUPRI,'(2X,A,2I5)') 'NCCCC, NPRF34', NCCCC, NPRF34
         DO 800 ICMP12 = 1, KHKT12
            IF (IODDCC(IPNTUV(ICMP12,2,1)) .EQ. IODD34) THEN
               CALL HEADER('PRMINT in CR2SEG',-1)
               WRITE (LUPRI,'(2X,A,I5)') 'ICMP12: ',ICMP12
               CALL OUTPUT(PRMINT(1,1,ICMP12),1,NCCCC,1,NPRF34,
     &                     NCCCC,NPRF34,1,LUPRI)
               CALL HEADER('CNTINT in CR2SEG',-1)
               WRITE (LUPRI,'(2X,A,I5)') 'ICMP12: ',ICMP12
               CALL OUTPUT(CNTINT(1,ICMP12),1,1,1,NCCCC,1,NCCCC,1,LUPRI)
            END IF
  800    CONTINUE
      END IF
C
      RETURN
      END
