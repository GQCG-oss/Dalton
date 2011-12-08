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
C FILE : eri2car2.F
C
#include <single.h>
C
C  /* Deck cr2drv */
      SUBROUTINE CR2DRV(AOINT,HCINT,INDHER,IODDHC,IODD34,IODDCC,
     &                  LMNPWR,IPNTUV,COOR34,EXP34,CSQ,CCFBT,NCENTR,
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
     &          COOR34(NPP34,3), EXP34(NPP34,3), CSQ(*), CCFBT(*),
     &          NCENTR(NPQBCX,2),
     &          WORK(LWORK)
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
#include <r12int.h>
C     LOCU12 is a local copy of U12INT, which is set .false.
C     for the transformation from Hermite to Cartesian functions
C     for electron 2 (WK/UniKA/04-11-2002).
      LOGICAL LOCU12
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
         JMXD3  = JMAX3  + NDER3
         JMXD4  = JMAX4  + NDER4
         JMXD34 = JMAX34 + NDER3 + NDER4 
         KCKT34 = MAX(NPAIRS(JMAX3+1)*NPAIRS(JMXD4+1),
     &                NPAIRS(JMXD3+1)*NPAIRS(JMAX4+1))
         KHC3 = NPAIRS(JMAX3+1)*NPAIRS(JMXD3+1)
         KHC4 = NPAIRS(JMAX4+1)*NPAIRS(JMXD4+1)
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
         LC34   = 0
         LCPINT = 0
C        Make local copy of U12INT (WK/UniKA/04-11-2002).
         LOCU12 = U12INT
         U12INT = .FALSE.
C
         LFIRST = KHKT12
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
         IF (GCON34) LCPINT = NPRF4*NCTF3*NCTF12*NPQBCX
C
         LTOTAL = LFIRST + LCCONT 
     &                   + MAX(
     &                         LCSINT + LFCSNT + LFAONT + LC34,
     &                         LECOEF + MAX(LHPI + LP3 + LP4,
     &                                      LEUV+LETUV+LCCPRM+LCPINT
     &                                     )
     &                        )
         IF (LTOTAL.GT.LWORK) CALL STOPIT('CR2DRV',' ',LTOTAL,LWORK)
C
         KFIRST = 1
         KODDKC = KFIRST + LFIRST 
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
         KCPINT = KCCPRM + LCCPRM
         KCCPDE = KCCPRM + LCCPRM
         KCCRDE = KCCPDE + LCCPRD
C 
         IF (IELCT2 .EQ. 1) THEN
            ICORCD = 1
            IEXPCD = 1
         ELSE
            ICORCD = 2
            IEXPCD = 3
         END IF
C
         CALL CR2TWO(AOINT,WORK(KFAONT),
     &               HCINT,INDHER,IODDHC,IODD34,IODDCC,WORK(KODDKC),
     &               WORK(KCCPDE),WORK(KCCRDE),WORK(KCCPRM),
     &               WORK(KFIRST),WORK(KCCONT),
     &               WORK(KCSINT),WORK(KFCSNT),
     &               WORK(KECOEF),WORK(KEUV),WORK(KETUV),
     &               LMNPWR,IPNTUV,COOR34,EXP34,CSQ,
     &               WORK(KCA3),WORK(KCV3),WORK(KCA4),WORK(KCV4),
     &               NCENTR,CCFBT,WORK(KCPINT),
     &               WORK(KHPI),WORK(KP3),WORK(KP4),
     &               NDER3,NDER4,NDIMD,IPRINT,
     &               WORK,WORK,WORK,WORK)
         U12INT = LOCU12
C
      RETURN
      END
C  /* Deck cr2two */
      SUBROUTINE CR2TWO(AOINT,FAOINT,
     &                  HCINT,INDHER,IODDHC,IODD34,IODDCC,IODDKC,
     &                  CCPDER,CCRDER,
     &                  CCPRIM,FRSTUV,CCONT,CSINT,FCSINT,
     &                  ECOEF,EUV,ETUV,LMNPWR,IPNTUV,
     &                  COOR34,EXP34,
     &                  CSQ,CSA3,CSV3,CSA4,CSV4,NCENTR,
     &                  CCFBT,CPINT,
     &                  HPI,P3,P4,
     &                  NDER3,NDER4,NDIMD,IPRINT,
     &                  TEXPA,TEXPB,AOVERP,BOVERP)
C
C     tuh
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      LOGICAL FLEVEL, FRSTUV(KHKT12),
     &        FAOINT(KHKT12,KHKT34), FCSINT(KHKT12,KHKT4),
     &        BDER3, BDER4
      DIMENSION AOINT(NCCT,KHKTAB,KHKTCD,NDIMD),
     &          HCINT(NCTF12,NPP34,NTUV34,KHKT12),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP), IODDKC(KC2MAX),
     &          IODDHC(NRTOP), IODD34(NRTOP), IODDCC(NRTOP),
     &          CCPRIM(NCTF12,NPP34,KHKT12),
     &          CCONT(NCCCC,KHKT12,KCKT34),
     &          CSINT(NCCCC,KHKT12,KHKT4),
     &          ECOEF(NPP34,0:JMXD3+JMXD4,0:JMXD3,0:JMXD4,3),
     &          EUV(NPP34,2*MAXDER+1), ETUV(NPP34,0:3*MAXDER),
     &          LMNPWR(KCKMAX,NHKMAX,3), IPNTUV(KC2MAX,0:NRDER,2),
     &          COOR34(NPP34,3), EXP34(NPP34,3), 
     &          HPI(NPP34), P3(NPP34,3), P4(NPP34,3),
     &          CSQ(NCSQ1,NCSQ2),  NCENTR(NPQBCX,2),
     &          CSA3(KHKT3,KCKT3), CSV3(NCCCC,KHKT3,KCKT3), 
     &          CSA4(KHKT4,KCKT4), CSV4(NCCCC,KHKT4,KCKT4),
     &          CCFBT(*), CPINT(*),
     &          TEXPA(*), TEXPB(*), AOVERP(*), BOVERP(*)
C               TEXPA, TEXPB, AOVERP, and BOVERP are dummies (WK/UniKA/04-11-2002).
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from CR2TWO','*',103)
C
      IF (MAXDER.EQ.0 .AND. KHKT34 .EQ. 1) THEN
         DO ICMP12 = 1, KHKT12
            FRSTUV(ICMP12) = IODDHC(IPNTUV(ICMP12,0,IELCT1)).NE.0
         END DO
         CALL ERICT2(HCINT,AOINT,CCFBT(ICMAT3),CCFBT(ICMAT4),CPINT,
     &               NCCT,1,IODDCC,IPNTUV,0,FRSTUV,IPRINT)
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
               WRITE (LUPRI,'(//,2X,A,2I5/)') ' Angular levels:',
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
     &                  ISCAL1,ISCAL2,IPRINT,
     &                  TEXPA,TEXPB,AOVERP,BOVERP)
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
                  IF (SPHR34) THEN
                     CALL ERICT2(CCPRIM,CCONT(1,1,ICMP34),CCFBT(ICMAT3),
     &                           CCFBT(ICMAT4),CPINT,NCCCC,1,IODDCC,
     &                           IPNTUV,IODDKC(ICMP34),FRSTUV,IPRINT)
                  ELSE IF (IELCT2.EQ.2) THEN
                     CALL ERICT2(CCPRIM,AOINT(1,1,ICMP34,1),
     &                           CCFBT(ICMAT3),CCFBT(ICMAT4),CPINT,
     &                           NCCT,1,IODDCC,IPNTUV,IODDKC(ICMP34),
     &                           FRSTUV,IPRINT)
                  ELSE
                     CALL ERICT2(CCPRIM,AOINT(1,1,ICMP34,1),
     &                           CCFBT(ICMAT3),CCFBT(ICMAT4),CPINT,
     &                           NCCT,KHKT34,IODDCC,IPNTUV,
     &                           IODDKC(ICMP34),FRSTUV,IPRINT)
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
                     IF (BDER) THEN
                        IODX = IODDCC(INDHER(IY3+IY4+IZ3+IZ4,
     &                                       IX3+IX4+IZ3+IZ4,
     &                                       IX3+IX4+IY3+IY4))
                     ELSE
                        IODX = IODDCC(INDHER(IX3+IX4,IY3+IY4,IZ3+IZ4))
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
      DIMENSION HCINT(NCTF12,NPP34,NTUV34,KHKT12),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          IODDHC(NRTOP),
     &          CCPRIM(NCTF12,NPP34,KHKT12),
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
      DO ICMP12 = 1, KHKT12
         FRSTUV(ICMP12) = .TRUE.
      END DO
C
      DO 200 V = MINV, MAXV, INCV
      DO 200 U = MINU, MAXU, INCU
         DO I = 1, NPP34
            EUV(I) = ECOEF(I,V,N3,N4,3)*ECOEF(I,U,M3,M4,2)
         END DO
         DO 300 T = MINT, MAXT, INCT
            DO I = 1, NPP34
               ETUV(I) = ECOEF(I,T,L3,L4,1)*EUV(I)
            END DO
C
            ITUV = INDHER(T,U,V)
            IODD = IODDHC(ITUV)
C
            IF (NCTF12.GT.1) THEN
               DO ICMP12 = 1, KHKT12
               IF (IODDHC(IPNTUV(ICMP12,0,IELCT1)) .EQ. IODD) THEN
                  IF (FRSTUV(ICMP12)) THEN
                     FRSTUV(ICMP12) = .FALSE.
                     DO J = 1, NPP34
                     DO I = 1, NCTF12
                        CCPRIM(I,J,ICMP12) =
     &                    ETUV(J)*HCINT(I,J,ITUV,ICMP12)
                     END DO
                     END DO
                  ELSE
                     DO J = 1, NPP34
                     DO I = 1, NCTF12
                        CCPRIM(I,J,ICMP12) = CCPRIM(I,J,ICMP12)
     &                         + ETUV(J)*HCINT(I,J,ITUV,ICMP12)
                     END DO 
                     END DO
                  END IF
               END IF
               END DO
            ELSE
               DO ICMP12 = 1, KHKT12
               IF (IODDHC(IPNTUV(ICMP12,0,IELCT1)) .EQ. IODD) THEN
                  IF (FRSTUV(ICMP12)) THEN
                     FRSTUV(ICMP12) = .FALSE.
                     DO J = 1, NPP34
                        CCPRIM(1,J,ICMP12) =
     &                    ETUV(J)*HCINT(1,J,ITUV,ICMP12)
                     END DO
                  ELSE
                     DO J = 1, NPP34
                        CCPRIM(1,J,ICMP12) = CCPRIM(1,J,ICMP12)
     &                         + ETUV(J)*HCINT(1,J,ITUV,ICMP12)
                     END DO
                  END IF
               END IF
               END DO
            END IF
C
  300    CONTINUE
  200 CONTINUE
C
      IF (IPRINT .GT. 50) THEN
         CALL HEADER('CCPRIM - CR2UND',-1)
         DO ICMP12 = 1, KHKT12
            WRITE (LUPRI,'(1X,A,I5)') ' CCPRIM for ICMP12 =',ICMP12
            CALL OUTPUT(CCPRIM(1,1,ICMP12),1,NCTF12,1,NPP34,
     &                  NCTF12,NPP34,1,LUPRI)
         END DO
      END IF
C
      RETURN
      END
C  /* Deck erict2 */
      SUBROUTINE ERICT2(PPINT,CCINT,CMAT3,CMAT4,CPINT,NDIM1,NDIM2,
     &                  IODDCC,IPNTUV,IODKC,FRSTUV,IPRINT)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      LOGICAL FRSTUV(KHKT12)
      DIMENSION CMAT3(NPRF3,NCTF3)
      DIMENSION CMAT4(NPRF4,NCTF4)
      DIMENSION PPINT(NCTF12*NPQBCX,NPRF34,KHKT12)
      DIMENSION CCINT(NDIM1,NDIM2,KHKT12)
      DIMENSION CPINT(NCTF3,NCTF12*NPQBCX,NPRF4)
      DIMENSION IPNTUV(KC2MAX,0:NRDER,2), IODDCC(NRTOP)
#include <hertop.h>
#include <ericom.h>
#include <eriao.h>
C
      NCC12 = NCTF12*NPQBCX
C
      IF (NPRF34.EQ.1 .AND. .NOT.GCON34) THEN
         DO ICMP12 = 1, KHKT12
         IF (IODDCC(IPNTUV(ICMP12,0,IELCT1)) .EQ. IODKC) THEN
            IF (FRSTUV(ICMP12)) THEN
               DO J = 1, NCCCC
                  CCINT(J,1,ICMP12) = D0
               END DO
            ELSE
               DO J = 1, NCCCC
                  CCINT(J,1,ICMP12) = PPINT(J,1,ICMP12)
               END DO
            END IF
         END IF
         END DO
      ELSE
         DO 400 ICMP12 = 1, KHKT12
         IF (IODDCC(IPNTUV(ICMP12,0,IELCT1)) .EQ. IODKC) THEN
            IF (FRSTUV(ICMP12)) THEN
               DO J = 1, NCCCC
                  CCINT(J,1,ICMP12) = D0
               END DO
            ELSE
               IF (GCON34) THEN
                  CALL CONT34(NCC12,PPINT(1,1,ICMP12),
     &                 CCINT(1,1,ICMP12),CMAT3,CMAT4,CPINT)
               ELSE
                  DO J = 1, NCCCC
                     CCINT(J,1,ICMP12) = PPINT(J,1,ICMP12)
     &                                 + PPINT(J,2,ICMP12)
                  END DO
                  DO I = 3, NPRF34 
                  DO J = 1, NCCCC
                     CCINT(J,1,ICMP12) = CCINT(J,1,ICMP12)
     &                                 + PPINT(J,I,ICMP12)
                  END DO
                  END DO
               END IF
            END IF
         END IF
  400    CONTINUE
      END IF
C
      IF (IPRINT .GT. 25) THEN
         CALL HEADER('Output from ERICT2',-1)
         WRITE (LUPRI,'(2X,A,3I5)') 'NCCCC, NPRF34, NDIM2', 
     &                               NCCCC, NPRF34, NDIM2
         WRITE (LUPRI,'(2X,A,2I5)') 'NPRF3,NCTF3',NPRF3,NCTF3
         WRITE (LUPRI,'(2X,A,2I5)') 'NPRF4,NCTF4',NPRF4,NCTF4
         WRITE (LUPRI,'(2X,A,2I5)') 'NCC12,KHKT12',NCC12,KHKT12
         IF (GCON3) THEN
            CALL HEADER('CMAT3 in ERICT2',-1)
            CALL OUTPUT(CMAT3,1,NPRF3,1,NCTF3,NPRF3,NCTF3,1,LUPRI)
         END IF
         IF (GCON4) THEN
            CALL HEADER('CMAT4 in ERICT2',-1)
            CALL OUTPUT(CMAT4,1,NPRF4,1,NCTF4,NPRF4,NCTF4,1,LUPRI)
         END IF
         DO 800 ICMP12 = 1, KHKT12
            IF (IODDCC(IPNTUV(ICMP12,0,IELCT1)) .EQ. IODKC) THEN
               CALL HEADER('PPINT in ERICT2',-1)
               WRITE (LUPRI,'(2X,A,I5)') 'ICMP12: ',ICMP12
               CALL OUTPUT(PPINT(1,1,ICMP12),1,NCCCC,1,NPRF34,
     &                     NCCCC,NPRF34,1,LUPRI)
               CALL HEADER('CCINT in ERICT2',-1)
               WRITE (LUPRI,'(2X,A,I5)') 'ICMP12: ',ICMP12
               CALL OUTPUT(CCINT(1,1,ICMP12),1,1,1,NCCCC,
     &                     1,NCCCC,1,LUPRI)
            END IF
  800    CONTINUE
      END IF
C
      RETURN
      END

C  /* Deck cont34 */
      SUBROUTINE CONT34(NCC12,PPINT,CCINT,CMAT3,CMAT4,CPINT)
#include <implicit.h>
      DIMENSION CMAT3(NPRF3,NCTF3)
      DIMENSION CMAT4(NPRF4,NCTF4)
      DIMENSION PPINT(NCC12*NPRF4,NPRF3)
      DIMENSION CCINT(NCTF4,NCC12*NCTF3)
      DIMENSION CPINT(NCTF3,NCC12,NPRF4)
#include <ericom.h>
#include <eriao.h>
      CALL DZERO(CPINT,NCTF3*NCC12*NPRF4)
      DO 200 K = 1, NPRF3
         DO 200 M = 1, NCTF3
            DO 200 N = 1, NCC12*NPRF4
C           ... "dirty trick" for CPINT to get long loop length
C               on last two indices (will go out of bounds on 2. index)
               CPINT(M,N,1) = CPINT(M,N,1) +
     &              CMAT3(K,M)*PPINT(N,K)
 200  CONTINUE
      CALL DZERO(CCINT(1,1),NCTF3*NCC12*NCTF4)
      DO 300 K = 1, NPRF4
         DO 300 M = 1, NCTF4
            DO 300 N = 1, NCTF3*NCC12
C           ... "dirty trick" for CPINT to get long loop length
C               on first two indices (will go out of bounds on 1. index)
               CCINT(M,N) = CCINT(M,N) + 
     &              CMAT4(K,M)*CPINT(N,1,K)
 300  CONTINUE
      RETURN
      END
C --- end of eri2car2.F ---
