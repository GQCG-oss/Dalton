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
C$Id: eri2car2.F,v 1.2 2001-02-12 18:17:53 vebjornb Exp $
C
#include <single.h>
C
C  /* Deck cr2drv */
      SUBROUTINE CR2DRV(AOINT,HCINT,INDHER,IODDHC,IODD34,IODDCC,
     &                  LMNPWR,IPNTUV,COOR34,EXP34,CSQ,NCENTR,
     &                  NDER3,NDER4,NDIMD,WORK,LWORK,IPRINT)
C
C     A. Halkier & T. Helgaker 4/2 1999. Modified to do
C     derivative integrals as well, as determined by MAXDER.
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      DIMENSION AOINT(NCCT,KHKTAB,KHKTCD,NDIMD),
     &          HCINT(NCCPP,NTUV34,KHKT12),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          IODDHC(NRTOP), IODD34(NRTOP), IODDCC(NRTOP),
     &          LMNPWR(KCKMAX,NHKMAX,3), IPNTUV(KC2MAX,0:NRDER,2),
     &          COOR34(NPP34,3), EXP34(NPP34,3), CSQ(*),
     &          NCENTR(NCCT,4),
     &          WORK(LWORK)
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from CR2DRV','*',103)
C
C        Work space allocation in CR2TWX
C
C        CCONT  | CSINT  | FCSINT | FAOINT
C               |
C               | ECOEF  | HPI    | P3     | P4
C                        | EUV    | ETUV   | CCPRIM
C
         IF (NEWDER) THEN
            JMXD3  = JMAX3  + NDER3
            JMXD4  = JMAX4  + NDER4
            JMXD34 = JMAX34 + NDER3 + NDER4 
            KCKT34 = MAX(NPAIRS(JMAX3+1)*NPAIRS(JMXD4+1),
     &                   NPAIRS(JMXD3+1)*NPAIRS(JMAX4+1))
            KHC3 = NPAIRS(JMAX3+1)*NPAIRS(JMXD3+1)
            KHC4 = NPAIRS(JMAX4+1)*NPAIRS(JMXD4+1)
         ELSE
            JMXD3  = JMAX3
            JMXD4  = JMAX4
            JMXD34 = JMAX34
         END IF
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
         LCCPRD = 0
         LCCPRM = 0
         LCA3   = 0
         LCV3   = 0
         LCA4   = 0
         LCV4   = 0
C
         LFIRST = 6*KHKT12
C
         IF (KHKT34 .GT. 1 .OR. MAXDER .GT.0) THEN
            IF (SPHR34) THEN
               LCCONT = NDIMD*NCCCC*KHKT12*KCKT34
               LFAONT = KHKT12*KHKT34
            END IF
            IF (SPHR3 .AND. SPHR4) THEN
               LCSINT = NCCCC*KHKT12*KHKT4
               LFCSNT = KHKT12*KHKT4
            END IF
C
            LECOEF = (MAXDER+1)*3*NPP34*(JMXD34+1)*(JMXD3+1)*(JMXD4+1)
            LHPI   =   NPP34
            LP3    = 3*NPP34
            LP4    = 3*NPP34
C
            LEUV   = (2*MAXDER+1)*NPP34
            LETUV  = (3*MAXDER+1)*NPP34
C
            LCCPRD = 3*MAXDER*NCCPP*KHKT12
            LCCPRM = NDIMD*NCCPP*KHKT12
C
            IF (BDER) THEN
               LCA3 = KHC3
               LCA4 = KHC4
               LCV3 = NCCCC*LCA3
               LCV4 = NCCCC*LCA4
               LC34 = LCA3 + LCV3 + LCA4 + LCV4 
C              IF (MAXDER .GT. 1) THEN
C                 LC34 = LCA3 + LCV3 + LCA4 + LCV4 
C              ELSE
C                 LC34 = MAX(LCA3 + LCV3, LCA4 + LCV4)
C              END IF
            END IF 
         END IF
C
         LTOTAL = LCCONT + MAX(
     &                         LCSINT + LFCSNT + LFAONT + LC34,
     &                         LECOEF + MAX(
     &                                      LHPI + LP3 + LP4,
     &                                      LEUV+LETUV+LCCPRM+LFIRST
     &                                     )
     &                        )
         IF (LTOTAL.GT.LWORK) CALL STOPIT('CR2DRV',' ',LTOTAL,LWORK)
C
         KODDKC = 1
         KCCONT = KODDKC + 4*KC2MAX
         KCSINT = KCCONT + LCCONT
         KFCSNT = KCSINT + LCSINT
         KFAONT = KFCSNT + LFCSNT
         KCA3   = KFAONT + LFAONT
         KCV3   = KCA3   + LCA3
         KCA4   = KCV3   + LCV3 
         KCV4   = KCA4   + LCA4
         KECOEF = KCCONT + LCCONT
         KHPI   = KECOEF + LECOEF
         KP3    = KHPI   + LHPI
         KP4    = KP3    + LP3
         KEUV   = KECOEF + LECOEF
         KETUV  = KEUV   + LEUV
         KCCPRM = KETUV  + LETUV
         KCCPDE = KCCPRM + LCCPRM
         KCCRDE = KCCPDE + LCCPRD
         KFIRST = KCCRDE + LCCPRD
C
         IF (IELCT2 .EQ. 1) THEN
            ICORCD = 1
            IEXPCD = 1
         ELSE
            ICORCD = 2
            IEXPCD = 3
         END IF
C
         IF (NEWDER) THEN
            CALL CR2TWO(AOINT,WORK(KFAONT),
     &                  HCINT,INDHER,IODDHC,IODD34,IODDCC,WORK(KODDKC),
     &                  WORK(KCCPDE),WORK(KCCRDE),WORK(KCCPRM),
     &                  WORK(KFIRST),WORK(KCCONT),
     &                  WORK(KCSINT),WORK(KFCSNT),
     &                  WORK(KECOEF),WORK(KEUV),WORK(KETUV),
     &                  LMNPWR,IPNTUV,COOR34,EXP34,
     &                  CSQ,WORK(KCA3),WORK(KCV3),WORK(KCA4),WORK(KCV4),
     &                  NCENTR,WORK(KHPI),WORK(KP3),WORK(KP4),
     &                  NDER3,NDER4,NDIMD,IPRINT)
         ELSE
            CALL CR2TWX(AOINT,WORK(KFAONT),
     &                  HCINT,INDHER,IODDHC,IODD34,IODDCC,WORK(KODDKC),
     &                  WORK(KCCPDE),WORK(KCCRDE),WORK(KCCPRM),
     &                  WORK(KFIRST),WORK(KCCONT),
     &                  WORK(KCSINT),WORK(KFCSNT),
     &                  WORK(KECOEF),WORK(KEUV),WORK(KETUV),
     &                  LMNPWR,IPNTUV,COOR34,EXP34,
     &                  CSQ,WORK(KCA3),WORK(KCV3),WORK(KCA4),WORK(KCV4),
     &                  NCENTR,WORK(KHPI),WORK(KP3),WORK(KP4),
     &                  NDIMD,IPRINT)
         END IF
C
      RETURN
      END
C  /* Deck cr2two */
      SUBROUTINE CR2TWO(AOINT,FAOINT,
     &                  HCINT,INDHER,IODDHC,IODD34,IODDCC,IODDKC,
     &                  CCPDER,CCRDER,
     &                  CCPRIM,FRSTUV,CCONT,CSINT,FCSINT,
     &                  ECOEF,EUV,ETUV,LMNPWR,IPNTUV,
     &                  COOR34,EXP34,CSQ,
     &                  CSA3,CSV3,CSA4,CSV4,NCENTR,
     &                  HPI,P3,P4,
     &                  NDER3,NDER4,NDIMD,IPRINT)
C
C     tuh
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      LOGICAL SAMMRY, FLEVEL, FRSTUV(KHKT12),
     &        FAOINT(KHKT12,KHKT34), FCSINT(KHKT12,KHKT4),
     &        BDER3, BDER4
      DIMENSION AOINT(NCCT,KHKTAB,KHKTCD,NDIMD),
     &          HCINT(NPP34,NCTF12,NTUV34,KHKT12),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP), IODDKC(KC2MAX),
     &          IODDHC(NRTOP), IODD34(NRTOP), IODDCC(NRTOP),
     &          CCPRIM(NPP34,NCTF12,KHKT12),
     &          CCONT(NCCCC,KHKT12,KCKT34),
     &          CSINT(NCCCC,KHKT12,KHKT4),
     &          ECOEF(NPP34,0:JMXD3+JMXD4,0:JMXD3,0:JMXD4,3),
     &          EUV(NPP34,2*MAXDER+1), ETUV(NPP34,0:3*MAXDER),
     &          LMNPWR(KCKMAX,NHKMAX,3), IPNTUV(KC2MAX,0:NRDER,2),
     &          COOR34(NPP34,3), EXP34(NPP34,3), 
     &          HPI(NPP34), P3(NPP34,3), P4(NPP34,3),
     &          CSQ(NCSQ1,NCSQ2),  NCENTR(NCCT,4),
     &          CSA3(KHKT3,KCKT3), CSV3(NCCCC,KHKT3,KCKT3), 
     &          CSA4(KHKT4,KCKT4), CSV4(NCCCC,KHKT4,KCKT4)
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from CR2TWO','*',103)
C
      IF (MAXDER.EQ.0 .AND. KHKT34 .EQ. 1) THEN
         SAMMRY = .FALSE.
         DO 20 ICMP12 = 1, KHKT12
            FRSTUV(ICMP12) = IODDHC(IPNTUV(ICMP12,0,IELCT1)).NE.0
   20    CONTINUE
         IODKC = 0
         NDIM2 = 1
         CALL CR2SEG(HCINT,AOINT,IODDCC,IPNTUV,IODKC,
     &               SAMMRY,FRSTUV,NCCT,NDIM2,IPRINT)
      ELSE
         NDER = MAX(NDER3,NDER4)
         IF (NDER.EQ.0) THEN
            J3D = 0 
            J3U = 0 
            J3I = 1
            J4D = 0 
            J4U = 0 
            J4I = 1
         ELSE
            IF (GDER) THEN
               J3D = NDER3 
               J3U = NDER3 
               J3I = 2
               J4D = NDER4 
               J4U = NDER4 
               J4I = 2
            END IF
            IF (BDER) THEN
               J3U = NDER3 
               J4U = NDER4 
               J3I = 2
               J4I = 2
               J3D = 0 
               J4D = 0 
            END IF 
         END IF
         NDMIN3 = NDER3
         NDMAX3 = NDER3
         NDMIN4 = NDER4
         NDMAX4 = NDER4
         FLEVEL = .TRUE.
         DO 100 J3 = NHKT3 + J3U, MAX(NHKT3 - J3D,1), -J3I
         DO 100 J4 = NHKT4 + J4U, MAX(NHKT4 - J4D,1), -J4I
C
            IF (IPRINT.GT. 10) THEN
               WRITE (LUPRI,'(//,2X,A,2I5/,)') ' Angular levels:',
     &               J3 - NHKT3, J4 - NHKT4
            END IF
            KCKT3 = J3*(J3 + 1)/2
            KCKT4 = J4*(J4 + 1)/2
C
C           Hermite-to-Cartesian expansion coefficients
C
            IDL    = 0
            IRUTIN = 2
            ISCAL1 = (J3 - NHKT3 + J3D)/2
            ISCAL2 = (J4 - NHKT4 + J4D)/2
            CALL EXPCFT(ECOEF,NPP34,JMXD3,JMXD4,COOR34,EXP34,I340,HPI,
     &                  P3,P4,IRUTIN,IELCT2,J3-1,J4-1,NCNT34,IDL,
     &                  ISCAL1,ISCAL2,IPRINT)
C
C           Run over components
C
            ICMP34 = 0
            DO 200 ICOMP3 = 1, KCKT3
               MAX4 = KCKT4
               IF (TCMP34) MAX4 = ICOMP3
               DO 300 ICOMP4 = 1, MAX4
                  ICMP34 = ICMP34 + 1
C
                  L3 = LMNPWR(ICOMP3,J3,1)
                  M3 = LMNPWR(ICOMP3,J3,2)
                  N3 = LMNPWR(ICOMP3,J3,3)
                  L4 = LMNPWR(ICOMP4,J4,1)
                  M4 = LMNPWR(ICOMP4,J4,2)
                  N4 = LMNPWR(ICOMP4,J4,3)
                  IODDKC(ICMP34) = IODDCC(INDHER(L3+L4,M3+M4,N3+N4))
C
C                 *****************************************
C                 ***** Cartesian Primitive Integrals *****
C                 *****************************************
C
                  CALL CR2UND(L3,M3,N3,L4,M4,N4,
     &                        HCINT,INDHER,IODDHC,CCPRIM,
     &                        FRSTUV,ECOEF,EUV,ETUV,IPNTUV,
     &                        ICOMP3,ICOMP4,IPRINT)
C
C                 ********************************
C                 ***** Contracted Integrals *****
C                 ********************************
C
                  SAMMRY = .FALSE.
                  IF (SPHR34) THEN
                     NDIM2 = 1
                     CALL CR2SEG(CCPRIM,CCONT(1,1,ICMP34),IODDCC,IPNTUV,
     &                           IODDKC(ICMP34),SAMMRY,FRSTUV,NCCCC,
     &                           NDIM2,IPRINT)
                  ELSE
                     IF (IELCT2.EQ.2) THEN
                        NDIM2 = 1
                        CALL CR2SEG(CCPRIM,AOINT(1,1,ICMP34,1),IODDCC,
     &                              IPNTUV,IODDKC(ICMP34),SAMMRY,FRSTUV,
     &                              NCCT,NDIM2,IPRINT)
                     ELSE
                        NDIM2 = KHKT34
                        CALL CR2SEG(CCPRIM,AOINT(1,ICMP34,1,1),IODDCC,
     &                              IPNTUV,IODDKC(ICMP34),SAMMRY,FRSTUV,
     &                              NCCT,NDIM2,IPRINT)
                     END IF
                  END IF
  300          CONTINUE
  200       CONTINUE
C
C           Spherical integrals
C           ===================
C
            INDX = 0
            IDER = 0
            IF (SPHR34) THEN
               DO 400 NDR3 = NDMIN3, NDMAX3
               DO 400 IX3 = NDR3, 0, -1
               DO 400 IY3 = NDR3 - IX3, 0, -1
                  DO 500 NDR4 = NDMIN4, NDMAX4
                  DO 500 IX4 = NDR4, 0, -1
                  DO 500 IY4 = NDR4 - IX4, 0, -1
                     IZ3 = NDR3 - IX3 - IY3
                     IZ4 = NDR4 - IX4 - IY4
                     IC3 = J3 - NHKT3
                     IC4 = J4 - NHKT4
                     IDER = IDER + 1
                     IF (.NOT.BDER) THEN
                        IODX = IODDCC(INDHER(IX3+IX4,IY3+IY4,IZ3+IZ4))
                     ELSE
                        IODX = IODDCC(INDHER(IY3+IY4+IZ3+IZ4,
     &                                       IX3+IX4+IZ3+IZ4,
     &                                       IX3+IX4+IY3+IY4))
                     END IF
                     BDER3 = BDER .AND. NDR3 .GT. 0 
                     BDER4 = BDER .AND. NDR4 .GT. 0 
                     CALL CR2SPH(CCONT,CSINT,FCSINT,
     &                           AOINT(1,1,1,IDER),FAOINT,FLEVEL,
     &                           CSQ(KSQADR(NHKT3-1,IX3,IY3,IZ3,IC3),1),
     &                           CSQ(KSQADR(NHKT4-1,IX4,IY4,IZ4,IC4),1),
     &                           CSA3,CSA4,CSV3,CSV4,BDER3,BDER4,NCENTR,
     &                           IODDCC,IPNTUV,IODDKC,INDX,IODX,IPRINT)
  500             CONTINUE
  400          CONTINUE
            END IF
C
            FLEVEL = .FALSE.
  100    CONTINUE
      END IF
C
C     *************************
C     ***** Print Section *****
C     *************************
C
      IF (IPRINT .GE. 15) THEN
         CALL HEADER('Final spherical integrals - CR2TWO',-1)
         DO 900 K = 1, NDIMD
            IF (NDER.GT.0) THEN
               INDX = MOD(K-1,3)+1
            ELSE
               INDX = 0 
            END IF 
            DO 910 I = 1, KHKTAB
            DO 910 J = 1, KHKTCD
            IF (IODDCC(IPNTUV(I,0,1)) .EQ. 
     &          IODDCC(IPNTUV(J,INDX,2))) THEN
               WRITE (LUPRI, '(/1X,A,3I3)') ' K,ICMPAB, ICMPCD ', K,I,J
               CALL OUTPUT(AOINT(1,I,J,K),1,1,1,NCCCC,1,NCCCC,1,LUPRI)
            END IF
  910       CONTINUE
  900    CONTINUE
      END IF
      RETURN
      END
C  /* Deck cr2twx */
      SUBROUTINE CR2TWX(AOINT,FAOINT,
     &                  HCINT,INDHER,IODDHC,IODD34,IODDCC,IODDKC,
     &                  CCPDER,CCRDER,
     &                  CCPRIM,FRSTUV,CCONT,CSINT,FCSINT,
     &                  ECOEF,EUV,ETUV,LMNPWR,IPNTUV,
     &                  COOR34,EXP34,CSQ,
     &                  CSA3,CSV3,CSA4,CSV4,NCENTR,
     &                  HPI,P3,P4,
     &                  NDIMD,IPRINT)
C
C     A. Halkier & T. Helgaker 4/2 1999. Modified to do
C     derivative integrals as well, as determined by MAXDER.
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      PARAMETER (D1 = 1.0D0)
      INTEGER T, U, V
      LOGICAL SAMMRY,  FLEVEL, FRSTUV(KHKT12,6),
     &        FAOINT(KHKT12,KHKT34), FCSINT(KHKT12,KHKT4),BDER3,BDER4
      DIMENSION AOINT(NCCT,KHKTAB,KHKTCD,NDIMD),
     &          HCINT(NPP34,NCTF12,NTUV34,KHKT12),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP), IODDKC(KC2MAX,0:3),
     &          IODDHC(NRTOP), IODD34(NRTOP), IODDCC(NRTOP),
     &          CCPRIM(NPP34,NCTF12,KHKT12,NDIMD),
     &          CCPDER(NPP34,NCTF12,KHKT12,3),
     &          CCRDER(NPP34,NCTF12,KHKT12,3),
     &          CCONT(NCCCC,KHKT12,KCKT34,NDIMD),
     &          CSINT(NCCCC,KHKT12,KHKT4),
     &          ECOEF(NPP34,0:JMXD3+JMXD4,0:JMXD3,0:JMXD4,3,MAXDER+1),
     &          EUV(NPP34,2*MAXDER+1), ETUV(NPP34,0:3*MAXDER),
     &          LMNPWR(KCKMAX,NHKMAX,3), IPNTUV(KC2MAX,0:NRDER,2),
     &          COOR34(NPP34,3), EXP34(NPP34,3), 
     &          HPI(NPP34), P3(NPP34,3), P4(NPP34,3),
     &          CSQ(NCSQ1,NCSQ2), NCENTR(NCCT,4),
     &          CSA3(KHKT3,KCKT3), CSV3(NCCCC,KHKT3,KCKT3), 
     &          CSA4(KHKT4,KCKT4), CSV4(NCCCC,KHKT4,KCKT4)
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from CR2TWX','*',103)
C
      IF (MAXDER.EQ.0 .AND. KHKT34 .EQ. 1) THEN
         SAMMRY = .FALSE.
         IODKC = 0
         DO 20 ICMP12 = 1, KHKT12
            FRSTUV(ICMP12,1) = IODDHC(IPNTUV(ICMP12,0,IELCT1)).NE.0
   20    CONTINUE
         NDIM2 = 1
         CALL CR2SEG(HCINT,AOINT,IODDCC,IPNTUV,IODKC,
     &               SAMMRY,FRSTUV,NCCT,NDIM2,IPRINT)
      ELSE
C
C        Hermite-to-Cartesian expansion coefficients
C
         IRUTIN = 2
         ISCAL1 = 0
         ISCAL2 = 0
         CALL EXPCFT(ECOEF,NPP34,JMXD3,JMXD4,COOR34,EXP34,I340,HPI,
     &               P3,P4,IRUTIN,IELCT2,JMXD3,JMXD4,NCNT34,MAXDER,
     &               ISCAL1,ISCAL2,IPRINT)
C
C        Run over components
C
         FLEVEL = .TRUE.
         ICMP34 = 0
         DO 100 ICOMP3 = 1, KCKT3
            MAX4 = KCKT4
            IF (TCMP34) MAX4 = ICOMP3
            DO 200 ICOMP4 = 1, MAX4
               ICMP34 = ICMP34 + 1
C
               L3 = LMNPWR(ICOMP3,NHKT3,1)
               M3 = LMNPWR(ICOMP3,NHKT3,2)
               N3 = LMNPWR(ICOMP3,NHKT3,3)
               L4 = LMNPWR(ICOMP4,NHKT4,1)
               M4 = LMNPWR(ICOMP4,NHKT4,2)
               N4 = LMNPWR(ICOMP4,NHKT4,3)
C
C              *****************************************
C              ***** Cartesian Primitive Integrals *****
C              *****************************************
C
               IF (MAXDER .EQ. 0) THEN
                  IODDKC(ICMP34,0) = IODD34(INDHER(L3+L4,M3+M4,N3+N4))
                  CALL CR2UND(L3,M3,N3,L4,M4,N4,
     &                        HCINT,INDHER,IODDHC,CCPRIM,
     &                        FRSTUV,ECOEF,EUV,ETUV,IPNTUV,
     &                        ICOMP3,ICOMP4,IPRINT)
               ELSE IF (MAXDER .EQ. 1) THEN
                  IODDKC(ICMP34,0) = IODD34(INDHER(L3+L4,M3+M4,N3+N4))
                  IODDKC(ICMP34,1) = IODDCC(INDHER(L3+L4+1,M3+M4,N3+N4))
                  IODDKC(ICMP34,2) = IODDCC(INDHER(L3+L4,M3+M4+1,N3+N4))
                  IODDKC(ICMP34,3) = IODDCC(INDHER(L3+L4,M3+M4,N3+N4+1))
                  CALL CR2DER(L3,M3,N3,L4,M4,N4,
     &                        HCINT,INDHER,IODDHC,IODD34,CCPDER,CCRDER,
     &                        FRSTUV(1,1),FRSTUV(1,4),ECOEF,EUV,ETUV,
     &                        IPNTUV,ICOMP3,ICOMP4,IODDKC(ICMP34,0),
     &                        IPRINT)
                  CALL GABDER(CCPDER,CCRDER,CCPRIM(1,1,1,1),
     &                        CCPRIM(1,1,1,4),EXP34,IODDCC,IPNTUV,
     &                        ICMP34,IODDKC,FRSTUV(1,1),FRSTUV(1,4),
     &                        IPRINT)
               ENDIF
C
C              ********************************
C              ***** Contracted Integrals *****
C              ********************************
C
               DO 275 IDIMD = 1,NDIMD
                  IF (MAXDER .EQ. 0) THEN
                     IODKC = IODDKC(ICMP34,0)
                  ELSE
                     IODKC = IODDKC(ICMP34,MOD(IDIMD-1,3)+1) 
                  END IF 
                  SAMMRY = .FALSE.
                  IF (SPHR34) THEN
                     NDIM2 = 1
                     CALL CR2SEG(CCPRIM(1,1,1,IDIMD),
     &                           CCONT(1,1,ICMP34,IDIMD),IODDCC,IPNTUV,
     &                           IODKC,SAMMRY,FRSTUV(1,IDIMD),
     &                           NCCCC,NDIM2,IPRINT)
                  ELSE
                     IF (IELCT2.EQ.2) THEN
                        NDIM2 = 1
                        CALL CR2SEG(CCPRIM(1,1,1,IDIMD),
     &                              AOINT(1,1,ICMP34,IDIMD),IODDCC,
     &                              IPNTUV,IODKC,SAMMRY,
     &                              FRSTUV(1,IDIMD),NCCT,NDIM2,IPRINT)
                     ELSE
                        NDIM2 = KHKT34
                        CALL CR2SEG(CCPRIM(1,1,1,IDIMD),
     &                              AOINT(1,ICMP34,1,IDIMD),IODDCC,
     &                              IPNTUV,IODKC,SAMMRY,
     &                              FRSTUV(1,IDIMD),NCCT,NDIM2,IPRINT)
                     END IF
                  END IF
  275          CONTINUE
  200       CONTINUE
  100    CONTINUE
C
C        Spherical integrals
C        ===================
C
         DO 110 IDIMD = 1,NDIMD
            INDX = MOD(IDIMD-1,3)+1
            IF (SPHR34) THEN
               BDER3 = BDER .AND. NDR3 .GT. 0 
               BDER4 = BDER .AND. NDR4 .GT. 0 
               CALL CR2SPH(CCONT(1,1,1,IDIMD),CSINT,FCSINT,
     &                     AOINT(1,1,1,IDIMD),FAOINT,FLEVEL,
     &                     CSQ(KSQADR(NHKT3-1,0,0,0,0),1),
     &                     CSQ(KSQADR(NHKT4-1,0,0,0,0),1),
     &                     CSA3,CSA4,CSV3,CSV4,BDER3,BDER4,NCENTR,
     &                     IODDCC,IPNTUV,IODDKC(1,INDX),INDX,0,IPRINT)
            END IF
  110    CONTINUE
      END IF
C
C     *************************
C     ***** Print Section *****
C     *************************
C
      IF (IPRINT .GE. 15) THEN
         CALL HEADER('Final spherical integrals - CR2TWX',-1)
         DO 900 K = 1, NDIMD
            INDX = MOD(K-1,3)+1
            DO 910 I = 1, KHKTAB
            DO 910 J = 1, KHKTCD
            IF (IODDCC(IPNTUV(I,0,1)) .EQ. 
     &          IODDCC(IPNTUV(J,INDX,2))) THEN
               WRITE (LUPRI, '(/1X,A,3I3)') ' K,ICMPAB, ICMPCD ', K,I,J
               CALL OUTPUT(AOINT(1,I,J,K),1,1,1,NCCCC,1,NCCCC,1,LUPRI)
            END IF
  910       CONTINUE
  900    CONTINUE
      END IF
      RETURN
      END
C  /* Deck cr2und */
      SUBROUTINE CR2UND(L3,M3,N3,L4,M4,N4,
     &                  HCINT,INDHER,IODDHC,CCPRIM,FRSTUV,
     &                  ECOEF,EUV,ETUV,IPNTUV,
     &                  ICOMP3,ICOMP4,IPRINT)
C
C     A. Halkier and T. Helgaker 4/2 - 1999. Calculate contraction
C     of Hermite-integrals with hermite/cartesian expansion
C     coefficients for undifferentiated integrals.
C
C     Based on old ERI/HERMIT routines
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      INTEGER T, U, V
      LOGICAL FRSTUV(KHKT12)
      DIMENSION HCINT(NPP34,NCTF12,NTUV34,KHKT12),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          IODDHC(NRTOP),
     &          CCPRIM(NPP34,NCTF12,KHKT12),
     &          ECOEF(NPP34,0:JMXD3+JMXD4,0:JMXD3,0:JMXD4,3),
     &          EUV(NPP34), ETUV(NPP34),
     &          IPNTUV(KC2MAX,0:NRDER,2)
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from CR2UND','*',103)
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
         WRITE(LUPRI,'(/,1X,A,2I5/)')' ICOMP3, ICOMP4',
     &                                 ICOMP3, ICOMP4
         WRITE(LUPRI,'(1X,A,15X,3I5)')' T loop:',MINT,MAXT,INCT
         WRITE(LUPRI,'(1X,A,15X,3I5)')' U loop:',MINU,MAXU,INCU
         WRITE(LUPRI,'(1X,A,15X,3I5)')' V loop:',MINV,MAXV,INCV
      END IF
C
      DO 100 ICMP12 = 1, KHKT12
         FRSTUV(ICMP12) = .TRUE.
  100 CONTINUE
C
      DO 200 V = MINV, MAXV, INCV
      DO 200 U = MINU, MAXU, INCU
         DO 210 I = 1, NPP34
            EUV(I) = ECOEF(I,V,N3,N4,3)*ECOEF(I,U,M3,M4,2)
  210    CONTINUE
         DO 300 T = MINT, MAXT, INCT
            DO 310 I = 1, NPP34
               ETUV(I) = ECOEF(I,T,L3,L4,1)*EUV(I)
  310       CONTINUE
C
            ITUV = INDHER(T,U,V)
            IODD = IODDHC(ITUV)
            DO 400 ICMP12 = 1, KHKT12
               IF (IODDHC(IPNTUV(ICMP12,0,IELCT1)) .EQ. IODD) THEN
                  IF (FRSTUV(ICMP12)) THEN
                     FRSTUV(ICMP12) = .FALSE.
                     DO 410 J = 1, NCTF12
                     DO 410 I = 1, NPP34
                        CCPRIM(I,J,ICMP12) =
     &                    ETUV(I)*HCINT(I,J,ITUV,ICMP12)
  410                CONTINUE
                  ELSE
                     DO 440 J = 1, NCTF12
                     DO 440 I = 1, NPP34
                        CCPRIM(I,J,ICMP12) = CCPRIM(I,J,ICMP12)
     &                         + ETUV(I)*HCINT(I,J,ITUV,ICMP12)
  440                CONTINUE
                  END IF
               END IF
  400       CONTINUE
C
  300    CONTINUE
  200 CONTINUE
C
      RETURN
      END
C  /* Deck cr2der */
      SUBROUTINE CR2DER(L3,M3,N3,L4,M4,N4,
     &                  HCINT,INDHER,IODDHC,IODD34,CCP,CCR,
     &                  FIRSTR,FIRSTP,ECOEF,EUV,ETUV,IPNTUV,
     &                  ICOMP3,ICOMP4,IOD340,IPRINT)
C
C     A. Halkier and T. Helgaker 4/2 - 1999. 
C
C     Calculate contraction of Hermite-integrals with 
C     hermite/cartesian expansion coefficients for derivative 
C     integrals (order given by MAXDER).
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      INTEGER T, U, V, X, TUV(3), IOD12X(3), IOD34(0:3)
      LOGICAL FIRSTR(KHKT12,3), FIRSTP(KHKT12,3)
      DIMENSION HCINT(NPP34,NCTF12,NTUV34,KHKT12),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          IODDHC(NRTOP),IODD34(NRTOP),
     &          CCP(NPP34,NCTF12,KHKT12,3),
     &          CCR(NPP34,NCTF12,KHKT12,3),
     &          ECOEF(NPP34,0:JMXD3+JMXD4,0:JMXD3,0:JMXD4,3,2),
     &          EUV(NPP34,3), ETUV(NPP34,0:3),
     &          IPNTUV(KC2MAX,0:NRDER,2)
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from CR2DER','*',103)
C
      IOD12X(1) = IODDHC(INDHER(1,0,0))
      IOD12X(2) = IODDHC(INDHER(0,1,0))
      IOD12X(3) = IODDHC(INDHER(0,0,1))
C
      L34 = L3 + L4
      M34 = M3 + M4
      N34 = N3 + N4
C
      DO 100 ICMP12 = 1, KHKT12
      DO 100 X = 1, 3
         FIRSTP(ICMP12,X) = .TRUE.
         FIRSTR(ICMP12,X) = .TRUE.
  100 CONTINUE
C
      DO 200 V = 0, N34
      DO 200 U = 0, M34
         DO 220 I = 1, NPP34
            EUV(I,1) = ECOEF(I,V,N3,N4,3,1)*ECOEF(I,U,M3,M4,2,1)
            EUV(I,2) = ECOEF(I,V,N3,N4,3,1)*ECOEF(I,U,M3,M4,2,2)
            EUV(I,3) = ECOEF(I,V,N3,N4,3,2)*ECOEF(I,U,M3,M4,2,1)
  220    CONTINUE
         DO 300 T = 0, L34
C
            ITUV   = INDHER(T,U,V)
            TUV(1) = INDHER(T+1,U,V)
            TUV(2) = INDHER(T,U+1,V)
            TUV(3) = INDHER(T,U,V+1)
C
            IOD34(0) = IBTXOR(IOD340,IODD34(ITUV))
            IOD34(1) = IBTXOR(IOD340,IODD34(TUV(1)))
            IOD34(2) = IBTXOR(IOD340,IODD34(TUV(2)))
            IOD34(3) = IBTXOR(IOD340,IODD34(TUV(3)))
C
C           Expansion coefficients
C           ======================
C
            IF (IOD34(0).EQ.0) THEN
               DO 400 I = 1, NPP34
                  ETUV(I,0) = ECOEF(I,T,L3,L4,1,1)*EUV(I,1)
  400          CONTINUE
            END IF
            IF (IOD34(1).EQ.0) THEN
               DO 410 I = 1, NPP34
                  ETUV(I,1) = ECOEF(I,T,L3,L4,1,2)*EUV(I,1)
  410          CONTINUE
            END IF
            IF (IOD34(2).EQ.0) THEN
               DO 420 I = 1, NPP34
                  ETUV(I,2) = ECOEF(I,T,L3,L4,1,1)*EUV(I,2)
  420          CONTINUE
            END IF
            IF (IOD34(3).EQ.0) THEN
               DO 430 I = 1, NPP34
                  ETUV(I,3) = ECOEF(I,T,L3,L4,1,1)*EUV(I,3)
  430          CONTINUE
            END IF
C
C           Cartesian integrals
C           ===================
C
            IODHC  = IODDHC(ITUV)
            DO 500 K = 1, KHKT12
               IOD12 = IBTXOR(IODDHC(IPNTUV(K,0,IELCT1)),IODHC)
C
C              ***** P derivatives *****
C
               IF (IOD34(0).EQ.0) THEN
                  DO 600 X = 1, 3
                  IF (IOD12.EQ.IOD12X(X)) THEN
                     ITUVX = TUV(X)
                     IF (FIRSTP(K,X)) THEN
                        FIRSTP(K,X) = .FALSE.
                        DO 610 J = 1, NCTF12
                        DO 610 I = 1, NPP34
                           CCP(I,J,K,X) = ETUV(I,0)*HCINT(I,J,ITUVX,K)
  610                   CONTINUE
                     ELSE
                        DO 620 J = 1, NCTF12
                        DO 620 I = 1, NPP34
                           CCP(I,J,K,X) = CCP(I,J,K,X)
     &                                  + ETUV(I,0)*HCINT(I,J,ITUVX,K)
  620                   CONTINUE
                     END IF
                  END IF
  600             CONTINUE
               END IF 
C
C              ***** R derivatives *****
C
               IF (IOD12.EQ.0) THEN
                  DO 700 X = 1, 3
                  IF (IOD34(X).EQ.0) THEN
                     IF (FIRSTR(K,X)) THEN
                        FIRSTR(K,X) = .FALSE.
                        DO 710 J = 1, NCTF12
                        DO 710 I = 1, NPP34
                           CCR(I,J,K,X) = ETUV(I,X)*HCINT(I,J,ITUV,K)
  710                   CONTINUE
                     ELSE
                        DO 720 J = 1, NCTF12
                        DO 720 I = 1, NPP34
                           CCR(I,J,K,X) = CCR(I,J,K,X)
     &                                  + ETUV(I,X)*HCINT(I,J,ITUV,K)
  720                   CONTINUE
                     END IF
                  END IF
  700             CONTINUE
               END IF
C
  500       CONTINUE
  300    CONTINUE
  200 CONTINUE
C
      RETURN
      END
C  /* Deck gabder */
      SUBROUTINE GABDER(CCPDER,CCRDER,CCADER,CCBDER,EXP34,IODDCC,IPNTUV,
     &                  ICMP34,IODDKC,FIRSTR,FIRSTP,IPRINT)
C
C     Asger Halkier and T. Helgaker 3/2 1999. Transform derivative
C     integrals from P & R basis to A & B basis.
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      PARAMETER (D1 = 1.0D0)
      INTEGER X
      LOGICAL FRST, FIRSTR(KHKT12,3), FIRSTP(KHKT12,3)
      DIMENSION CCPDER(NPP34,NCTF12,KHKT12,3),
     &          CCRDER(NPP34,NCTF12,KHKT12,3),
     &          CCADER(NPP34,NCTF12,KHKT12,3),
     &          CCBDER(NPP34,NCTF12,KHKT12,3),
     &          EXP34(NPP34,3), IODDCC(NRTOP),
     &          IPNTUV(KC2MAX,0:NRDER,2),
     &          IODDKC(KC2MAX,0:3)
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from GABDER','*',103)
C
      IF (IPRINT .GT. 25) THEN
         CALL HEADER('CCPDER and CCRDER in GABDER',-1)
         DO 800 X = 1,3
            IODKC = IODDKC(ICMP34,X) 
            DO 810 ICMP12 = 1, KHKT12
            IF (IODDCC(IPNTUV(ICMP12,0,IELCT1)).EQ.IODKC) THEN
               WRITE (LUPRI,'(2X,A,2I5)') ' X, ICMP12 ', X, ICMP12
               CALL OUTPUT(CCPDER(1,1,ICMP12,X),
     &                1,NPP34,1,NCTF12,NPP34,NCTF12,1,LUPRI)
               CALL OUTPUT(CCRDER(1,1,ICMP12,X),
     &                1,NPP34,1,NCTF12,NPP34,NCTF12,1,LUPRI)
            END IF
  810       CONTINUE
  800    CONTINUE
      END IF
C
      SGN = D1
      IF (IELCT2.EQ.2) SGN = -D1
C
      IF (.TRUE.) THEN
      DO 100 X = 1,3
         IODKC = IODDKC(ICMP34,X) 
         DO 110 ICMP12 = 1, KHKT12
            FRST = .TRUE.
            IF (IODDCC(IPNTUV(ICMP12,0,IELCT1)).EQ.IODKC) THEN
               IF (.NOT.(FIRSTP(ICMP12,X).OR.FIRSTR(ICMP12,X))) THEN
                  DO 120 J = 1, NCTF12
                  DO 120 I = 1, NPP34
                     CCRINT = CCRDER(I,J,ICMP12,X)
                     CCPINT = SGN*CCPDER(I,J,ICMP12,X)
                     CCADER(I,J,ICMP12,X) = EXP34(I,1)*CCPINT + CCRINT 
                     CCBDER(I,J,ICMP12,X) = EXP34(I,2)*CCPINT - CCRINT 
  120             CONTINUE
                  FRST = .FALSE.
               ELSE IF (.NOT.FIRSTR(ICMP12,X)) THEN
                  DO 121 J = 1, NCTF12
                  DO 121 I = 1, NPP34
                     CCRINT = CCRDER(I,J,ICMP12,X)
                     CCADER(I,J,ICMP12,X) =   CCRINT 
                     CCBDER(I,J,ICMP12,X) = - CCRINT 
  121             CONTINUE
                  FRST = .FALSE.
               ELSE IF (.NOT.FIRSTP(ICMP12,X)) THEN
                  DO 122 J = 1, NCTF12
                  DO 122 I = 1, NPP34
                     CCPINT = SGN*CCPDER(I,J,ICMP12,X)
                     CCADER(I,J,ICMP12,X) = EXP34(I,1)*CCPINT
                     CCBDER(I,J,ICMP12,X) = EXP34(I,2)*CCPINT
  122             CONTINUE
                  FRST = .FALSE.
               END IF
            END IF
            FIRSTP(ICMP12,X) = FRST
            FIRSTR(ICMP12,X) = FRST
  110    CONTINUE
  100 CONTINUE
C
C     Another way of obtaining the B derivative...
C
      ELSE
      DO 130 X = 1,3
         IODKC = IODDKC(ICMP34,X) 
         DO 140 ICMP12 = 1, KHKT12
         IF (IODDCC(IPNTUV(ICMP12,0,IELCT1)).EQ.IODKC) THEN
            DO 150 J = 1, NCTF12
            DO 150 I = 1, NPP34
               CCADER(I,J,ICMP12,X) = EXP34(I,1)*CCPDER(I,J,ICMP12,X)
     *                                         + CCRDER(I,J,ICMP12,X)
               CCBDER(I,J,ICMP12,X) = CCPDER(I,J,ICMP12,X)
     *                              - CCADER(I,J,ICMP12,X)
  150       CONTINUE
         END IF
  140    CONTINUE
  130 CONTINUE
      END IF
C
      IF (IPRINT .GT. 25) THEN
         CALL HEADER('CCADER and CCBDER in GABDER',-1)
         DO 900 X = 1,3
            IODKC = IODDKC(ICMP34,X) 
            DO 910 ICMP12 = 1, KHKT12
            IF (IODDCC(IPNTUV(ICMP12,0,IELCT1)).EQ.IODKC) THEN
               WRITE (LUPRI,'(2X,A,2I5,//)') ' X, ICMP12 ', X, ICMP12
               CALL OUTPUT(CCADER(1,1,ICMP12,X),
     &                1,NPP34,1,NCTF12,NPP34,NCTF12,1,LUPRI)
               CALL OUTPUT(CCBDER(1,1,ICMP12,X),
     &                1,NPP34,1,NCTF12,NPP34,NCTF12,1,LUPRI)
            END IF
  910       CONTINUE
  900    CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck cr2seg */
      SUBROUTINE CR2SEG(PRMINT,CNTINT,IODDCC,IPNTUV,IODKC,SAMMRY,
     &                  FRSTUV,NDIM1,NDIM2,IPRINT)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0)
      LOGICAL SAMMRY, FRSTUV(KHKT12)
      DIMENSION PRMINT(NCCCC,NPRF34,KHKT12),
     &          CNTINT(NDIM1,NDIM2,KHKT12),
     &          IPNTUV(KC2MAX,0:NRDER,2), IODDCC(NRTOP)
#include <hertop.h>
#include <ericom.h>
#include <eriao.h>
C
      IF (NPRF34.EQ.1 .AND. .NOT.SAMMRY) THEN
         DO 100 ICMP12 = 1, KHKT12
         IF (IODDCC(IPNTUV(ICMP12,0,IELCT1)) .EQ. IODKC) THEN
            IF (FRSTUV(ICMP12)) THEN
               DO 200 J = 1, NCCCC
                  CNTINT(J,1,ICMP12) = D0
  200          CONTINUE
            ELSE
               DO 300 J = 1, NCCCC
                  CNTINT(J,1,ICMP12) = PRMINT(J,1,ICMP12)
  300          CONTINUE
            END IF
         END IF
  100    CONTINUE
      ELSE
         DO 400 ICMP12 = 1, KHKT12
         IF (IODDCC(IPNTUV(ICMP12,0,IELCT1)) .EQ. IODKC) THEN
            IF (FRSTUV(ICMP12)) THEN
               DO 500 J = 1, NCCCC
                  CNTINT(J,1,ICMP12) = D0
  500          CONTINUE
            ELSE
               DO 600 J = 1, NCCCC
                  CNTINT(J,1,ICMP12) = PRMINT(J,1,ICMP12)
     &                               + PRMINT(J,2,ICMP12)
  600          CONTINUE
               DO 700 I = 3, NPRF34
               DO 700 J = 1, NCCCC
                  CNTINT(J,1,ICMP12) = CNTINT(J,1,ICMP12)
     &                               + PRMINT(J,I,ICMP12)
  700          CONTINUE
            END IF
         END IF
  400    CONTINUE
      END IF
C
      IF (IPRINT .GT. 25) THEN
         CALL HEADER('Output from CR2SEG',-1)
         WRITE (LUPRI,'(2X,A,3I5)') 'NCCCC, NPRF34, NDIM2', 
     &                               NCCCC, NPRF34, NDIM2
         DO 800 ICMP12 = 1, KHKT12
            IF (IODDCC(IPNTUV(ICMP12,0,IELCT1)) .EQ. IODKC) THEN
               CALL HEADER('PRMINT in CR2SEG',-1)
               WRITE (LUPRI,'(2X,A,I5)') 'ICMP12: ',ICMP12
               CALL OUTPUT(PRMINT(1,1,ICMP12),1,NCCCC,1,NPRF34,
     &                     NCCCC,NPRF34,1,LUPRI)
               CALL HEADER('CNTINT in CR2SEG',-1)
               WRITE (LUPRI,'(2X,A,I5)') 'ICMP12: ',ICMP12
               CALL OUTPUT(CNTINT(1,1,ICMP12),1,1,1,NCCCC,
     &                     1,NCCCC,1,LUPRI)
            END IF
  800    CONTINUE
      END IF
C
      RETURN
      END
