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
C$Id: eri2car1.F,v 1.2 2001-02-12 18:17:53 vebjornb Exp $
C
#include <single.h>
C  /* Deck cr1drv */
      SUBROUTINE CR1DRV(HERINT,HCINT,INDHER,IODDHH,IODDHC,INDHSQ,
     &                  LMNPWR,IPNTUV,COOR12,EXP12,CSQ,
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
     &          COOR12(NPP12,3), EXP12(NPP12,3), CSQ(*),
     &          WORK(LWORK)
#include <cbieri.h>
#include <ericom.h>
#include <eriao.h>
      COMMON /DHCINF/ IHCADR(10), IHCSYM(10)
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from CR1DRV','*',103)
C
      ICORAB = 1
      IEXPAB = 1
      IF (IELCT1 .EQ. 2) THEN
         ICORAB = 2
         IEXPAB = 3
      END IF
C
      IF (.TRUE.) THEN
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
            LECOEF = 3*NPP12*(JMAX12 + 1)*(JMAX1 + 1)*(JMAX2+1)
            LHPI   =   NPP12
            LP1    = 3*NPP12
            LP2    = 3*NPP12
C
            LEUV   = NPP12
            LETUV  = NPP12
C
            LHCPRM = NPPPP*NTUV34
         END IF
C
         LFIRST = NTUV34
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
         KODDKC = 1
         KCCONT = KODDKC + KC2MAX
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
     &               LMNPWR,IPNTUV,WORK(KODDKC),
     &               COOR12,EXP12,CSQ,
     &               WORK(KHPI),WORK(KP1),WORK(KP2),
     &               HCINT,WORK(KFCINT),WORK(KHCPRM),WORK(KFIRST),
     &               WORK(KCCONT),WORK(KCSINT),WORK(KFCSNT),
     &               IPRINT)
      END IF
      RETURN
      END
C  /* Deck cr1two */
      SUBROUTINE CR1TWO(HERINT,INDHER,INDHVC,IODDHH,IODDHC,INDHSQ,
     &                  ECOEF,EUV,ETUV,LMNPWR,IPNTUV,IODDKC,
     &                  COOR12,EXP12,CSQ,HPI,P1,P2,
     &                  HCINT,FCINT,HCPRIM,FRSTUV,HCCONT,
     &                  HCSINT,FCSINT,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
C
      INTEGER TUV
      LOGICAL SAMMRY, FRSTUV(NTUV34),
     &        FCINT(NTUV34,KHKT12), FCSINT(NTUV34,KHKT2)
C
      COMMON /DHCINF/ IHCADR(10), IHCSYM(10)
C
      DIMENSION HERINT(NPP12,NPRF34,NTUV),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP), INDHVC(0:*),
     &          IODDHH(NRTOP), IODDHC(NRTOP), INDHSQ(NRTOP),
     &          HCINT(NCCPP,NTUV34,KHKT12),
     &          ECOEF(NPP12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2,3),
     &          ETUV(NPP12), EUV(NPP12),
     &          LMNPWR(KCKMAX,NHKMAX,3), IPNTUV(KC2MAX,0:NRDER,2),
     &          COOR12(NPP12,3), EXP12(NPP12,3), 
     &          HPI(NPP12), P1(NPP12,3), P2(NPP12,3),
     &          IODDKC(KC2MAX),
     &          CSQ(NCSQ1,NCSQ2),
     &          HCPRIM(NPP12,NPRF34,NTUV34),
     &          HCCONT(NCCPP,NTUV34,KCKT12),
     &          HCSINT(NCCPP,NTUV34,KHKT2)
C
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
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
         IODKC = 0
         DO 10 TUV = 1, NTUV34
            FRSTUV(TUV) = IODDHH(TUV) .NE. 0
   10    CONTINUE
         CALL CR1SEG(HERINT,HCINT,NPQBCX,NPRF12,NPRF34,
     &               IODDHC,IODKC,SAMMRY,FRSTUV,IPRINT)
C
C     *********************************
C     ***** General Case: (xy|zw) *****
C     *********************************
C
      ELSE
         J1D = 0
         J1U = 0
         J1I = 1
         J2D = 0
         J2U = 0
         J2I = 1
         DO 100 J1 = MAX(1,NHKT1 - J1D), NHKT1 + J1U, J1I
         DO 100 J2 = MAX(1,NHKT2 - J2D), NHKT2 + J1U, J2I 
            KCKT1 = J1*(J1 + 1)/2 
            KCKT2 = J2*(J2 + 1)/2 
C
C           Hermite-to-Cartesian expansion coefficients
C           ===========================================
C
            IRUTIN = 1
            IDL    = 0
            ISCAL1 = 0
            ISCAL2 = 0
            JMAX1 = J1 - 1
            JMAX2 = J2 - 1
            CALL EXPCFT(ECOEF,NPP12,JMAX1,JMAX2,COOR12,EXP12,I120,HPI,
     &                  P1,P2,IRUTIN,IELCT1,JMAX1,JMAX2,NCNT12,IDL,
     &                  ISCAL1,ISCAL2,IPRINT)
C
C           Run over components
C           ===================
C
            ICMP12 = 0
            DO 300 ICOMP1 = 1, KCKT1
               MAX2 = KCKT2
               IF (TCMP12) MAX2 = ICOMP1
               DO 400 ICOMP2 = 1, MAX2
                  ICMP12 = ICMP12 + 1
C
                  L1 = LMNPWR(ICOMP1,NHKT1,1)
                  M1 = LMNPWR(ICOMP1,NHKT1,2)
                  N1 = LMNPWR(ICOMP1,NHKT1,3)
                  L2 = LMNPWR(ICOMP2,NHKT2,1)
                  M2 = LMNPWR(ICOMP2,NHKT2,2)
                  N2 = LMNPWR(ICOMP2,NHKT2,3)
                  IODDKC(ICMP12) = IODDHC(INDHER(L1+L2,M1+M2,N1+N2))
C
C                 Primitive integrals
C                 ===================
C
                  CALL CR1UND(L1,M1,N1,L2,M2,N2,
     &                        HERINT,INDHER,INDHVC,IODDHH,INDHSQ,
     &                        ECOEF,EUV,ETUV,HCPRIM,FRSTUV,
     &                        ICOMP1,ICOMP2,IPRINT)
C
C                 Contracted integrals
C                 ====================
C
                  SAMMRY = .FALSE.
                  IF (SPHR12) THEN
                     CALL CR1SEG(HCPRIM,HCCONT(1,1,ICMP12),NPQBCX,
     &                           NPRF12,NPRF34,IODDHC,IODDKC(ICMP12),
     &                           SAMMRY,FRSTUV,IPRINT)
                  ELSE
                     CALL CR1SEG(HCPRIM,HCINT(1,1,ICMP12),NPQBCX,NPRF12,
     &                           NPRF34,IODDHC,IODDKC(ICMP12),SAMMRY,
     &                           FRSTUV,IPRINT)
                  END IF
C
  400          CONTINUE
  300       CONTINUE
C
C           Spherical integrals
C           ===================
C
            NDMIN1 = 0
            NDMAX1 = 0
            NDMIN2 = 0
            NDMAX2 = 0
            IF (SPHR12) THEN
               DO 500 NDER1 = NDMIN1, NDMAX1 
               DO 500 IX1 = NDER1, 0, -1   
               DO 500 IY1 = NDER1 - IX1, 0, -1
                  DO 600 NDER2 = NDMIN2, NDMAX2 
                  DO 600 IX2 = NDER2, 0, -1
                  DO 600 IY2 = NDER2 - IX2, 0, -1 
                     IZ1 = NDER1 - IX1 - IY1
                     IZ2 = NDER2 - IX2 - IY2
                     IC1 = J1 - NHKT1 
                     IC2 = J2 - NHKT2 
                     CALL CR1SPH(HCCONT,HCINT,FCINT,HCSINT,FCSINT,
     &                           CSQ(KSQADR(NHKT1-1,IX1,IY1,IZ1,IC1),1),
     &                           CSQ(KSQADR(NHKT2-1,IX2,IY2,IZ2,IC2),1),
     &                           IODDHC,IPNTUV,IODDKC,IPRINT)
  600             CONTINUE
  500          CONTINUE
            END IF
  100    CONTINUE
      END IF
      IF (IPRINT .GT. 20) THEN
         CALL HCPRINT(HCINT,IODDHC,IPNTUV,IPRINT)
      END IF
      RETURN
      END
C  /* Deck hcprint */
      SUBROUTINE HCPRINT(HCINT,IODDHC,IPNTUV,IPRINT)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (NHCTYP = 1)
      INTEGER TUV
      DIMENSION HCINT(NPQBCX,NPRF34,NCTF12,NTUV34,KHKT12,NHCTYP),
     &          IODDHC(NRTOP), IPNTUV(KC2MAX,0:NRDER,2)
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
               IF (IODDHC(TUV).EQ.IODDHC(IPNTUV(ICMP12,0,IELCT1)))THEN
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
C  /* Deck cr1und */
      SUBROUTINE CR1UND(L1,M1,N1,L2,M2,N2,
     &                  HERINT,INDHER,INDHVC,IODDHH,INDHSQ,
     &                  ECOEF,EUV,ETUV,HCPRIM,FRSTUV,
     &                  ICOMP1,ICOMP2,IPRINT)
C
C     T. Helgaker
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
C
      INTEGER T, U, V, TUV
      LOGICAL FRSTUV(NTUV34)
C
      DIMENSION HERINT(NPP12,NPRF34,NTUV),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP), INDHVC(0:*),
     &          IODDHH(NRTOP), INDHSQ(NRTOP),
     &          ECOEF(NPP12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2,3),
     &          ETUV(NPP12), EUV(NPP12),
     &          HCPRIM(NPP12,NPRF34,NTUV34)
C
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from CR1UND','*',103)
C
      INCT = I120(1) + 1
      INCU = I120(2) + 1
      INCV = I120(3) + 1
      MAXT = L1 + L2
      MAXU = M1 + M2
      MAXV = N1 + N2
      MINT = IBTAND(MAXT,I120(1))
      MINU = IBTAND(MAXU,I120(2))
      MINV = IBTAND(MAXV,I120(3))
C
      IF (IPRINT .GT. 25) THEN
         WRITE(LUPRI,'(/,1X,A,3I5/)')' ICOMP1, ICOMP2, ICMP12',
     &                                 ICOMP1, ICOMP2, ICMP12
         WRITE(LUPRI,'(1X,A,15X,3I5)')' T loop:',MINT,MAXT,INCT
         WRITE(LUPRI,'(1X,A,15X,3I5)')' U loop:',MINU,MAXU,INCU
         WRITE(LUPRI,'(1X,A,15X,3I5)')' V loop:',MINV,MAXV,INCV
      END IF
C
      DO 100 TUV = 1, NTUV34
         FRSTUV(TUV) = .TRUE.
  100 CONTINUE
C
      DO 200 V = MINV, MAXV, INCV
      DO 200 U = MINU, MAXU, INCU
         DO 210 I = 1, NPP12 
            EUV(I) = ECOEF(I,V,N1,N2,3) 
     &             * ECOEF(I,U,M1,M2,2)
  210    CONTINUE
         DO 300 T = MINT, MAXT, INCT
            DO 310 I = 1, NPP12
               ETUV(I) = ECOEF(I,T,L1,L2,1)*EUV(I)
  310       CONTINUE
C
            ITUV = INDHER(T,U,V)
            INDS = INDHSQ(ITUV)
            DO 400 TUV = 1, NTUV34
            IF (IODDHH(ITUV) .EQ. IODDHH(TUV)) THEN
#if defined (SYS_AIX)
C              code due to AIX xlf version 2.2 bug
               INDT = INDS + INDHSQ(TUV)
               INDT = INDHVC(INDT)
#else
               INDT = INDHVC(INDS + INDHSQ(TUV))
#endif
               IF (FRSTUV(TUV)) THEN
                  FRSTUV(TUV) = .FALSE.
                  DO 500 J = 1, NPRF34
                  DO 500 I = 1, NPP12
                     HCPRIM(I,J,TUV) = ETUV(I)*HERINT(I,J,INDT)
  500             CONTINUE
               ELSE
                  DO 600 J = 1, NPRF34
                  DO 600 I = 1, NPP12
                     HCPRIM(I,J,TUV) = HCPRIM(I,J,TUV)
     &                       + ETUV(I)*HERINT(I,J,INDT)
  600             CONTINUE
               END IF
            END IF
  400       CONTINUE
C
  300    CONTINUE
  200 CONTINUE
C
      RETURN
      END
C  /* Deck cr1seg */
      SUBROUTINE CR1SEG(PRMINT,CNTINT,NINNER,NPRIMS,NOUTER,
     &                  IODDHC,IODKC,SAMMRY,FRSTUV,IPRINT)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0)
      INTEGER TUV
      LOGICAL SAMMRY, FRSTUV(NTUV34)
      DIMENSION PRMINT(NINNER,NPRIMS,NOUTER,NTUV34),
     &          CNTINT(NINNER,NOUTER,NTUV34),
     &          IODDHC(NRTOP)
#include <hertop.h>
#include <ericom.h>
#include <eriao.h>
C
      IF (NPRIMS.EQ.1 .AND. .NOT.SAMMRY) THEN
         DO 100 TUV = 1, NTUV34
            IF (IODDHC(TUV) .EQ. IODKC) THEN
               IF (FRSTUV(TUV)) THEN
                  DO 200 K = 1, NOUTER
                  DO 200 J = 1, NINNER
                     CNTINT(J,K,TUV) = D0
  200             CONTINUE
               ELSE
                  DO 300 K = 1, NOUTER
                  DO 300 J = 1, NINNER
                     CNTINT(J,K,TUV) = PRMINT(J,1,K,TUV)
  300             CONTINUE
               END IF
            END IF
  100    CONTINUE
      ELSE
         DO 400 TUV = 1, NTUV34
         IF (IODDHC(TUV) .EQ. IODKC) THEN
            IF (FRSTUV(TUV)) THEN
               DO 500 K = 1, NOUTER
               DO 500 J = 1, NINNER
                  CNTINT(J,K,TUV) = D0
  500          CONTINUE
            ELSE
               DO 600 K = 1, NOUTER
                  DO 700 J = 1, NINNER
                     CNTINT(J,K,TUV)=PRMINT(J,1,K,TUV)+PRMINT(J,2,K,TUV)
  700             CONTINUE
                  DO 800 I = 3, NPRIMS
                  DO 800 J = 1, NINNER
                     CNTINT(J,K,TUV)=CNTINT(J,K,TUV) + PRMINT(J,I,K,TUV)
  800             CONTINUE
  600          CONTINUE
            END IF
         END IF
  400    CONTINUE
      END IF
C
      IF (IPRINT .GT. 25) THEN
         CALL HEADER('Output from CR1SEG',-1)
         WRITE (LUPRI,'(2X,A,3I5)') 
     &         'NINNER, NPRIMS, NOUTER', NINNER, NPRIMS, NOUTER
         DO 900 TUV = 1, NTUV34
            IF (IODDHC(TUV) .EQ. IODKC) THEN
               IF (.NOT.FRSTUV(TUV)) THEN
                  CALL HEADER('PRMINT in CR1SEG',-1)
                  WRITE (LUPRI,'(2X,A,I5)') 'TUV: ',TUV
                  CALL OUTPUT(PRMINT(1,1,1,TUV),
     &                        1,NINNER*NPRIMS,1,NOUTER,
     &                        NINNER*NPRIMS,NOUTER,1,LUPRI)
               END IF
               CALL HEADER('CNTINT in CR1SEG',-1)
               WRITE (LUPRI,'(2X,A,I5)') 'TUV: ',TUV
               CALL OUTPUT(CNTINT(1,1,TUV),
     &                     1,NINNER,1,NOUTER,
     &                     NINNER,NOUTER,1,LUPRI)
            END IF
  900    CONTINUE
      END IF
C
      RETURN
      END
