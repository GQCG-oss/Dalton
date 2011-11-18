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
C FILE: eri2car1.F
C
#include <single.h>
C
C  /* Deck cr1drv */
      SUBROUTINE CR1DRV(HERINT,HERR12,HCINT,INDHER,IODDHH,IODDHC,INDHSQ,
     &                  LMNPWR,IPNTUV,COOR12,EXP12,CSQ,CCFBT,
     &                  WORK,LWORK,IPRINT)
C
C     HERR12 integrals have been added (WK/UniKA/04-11-2002).
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
C
      LOGICAL DOREC, FRSTAO, CONTRA
      DIMENSION HERINT(*), HCINT(*), HERR12(*),
     &          INDHER(*), INDHSQ(*), IODDHH(*), IODDHC(*),
     &          LMNPWR(*), IPNTUV(*), 
     &          COOR12(NPP12,3), EXP12(NPP12,3), CSQ(*), CCFBT(*),
     &          WORK(LWORK)
C     Call to subroutine CR1UND or CR1U12 in CR1TWO (WK/UniKA/04-11-2002).
      EXTERNAL CR1UND, CR1U12
#include <cbieri.h>
#include <ericom.h>
#include <eriao.h>
#include <r12int.h>
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
C     HCCONT  | HCSINT | LFCSNT | LFCINT
C             |
C             | ECOEF | HPI | P1 | P2
C                     |
C                     | EUV | ETUV | HCPRIM
C
      LCCONT = 0
      LCSINT = 0
      LFCSNT = 0
      LFCINT = 0
      LECOEF = 0
      LHPI   = 0
      LP1    = 0
      LP2    = 0
C     Space for R12 method (WK/UniKA/04-11-2002).
      LXA    = 0
      LEUV   = 0
      LETUV  = 0
      LHCPRM = 0
      LCPINT = 0
C
c     IF (KHKT12 .GT. 1 .OR. U12INT. OR. MIAU12) THEN
      IF (KHKT12 .GT. 1 .OR. U12INT) THEN
C        No symmetry for [T1,r12] integrals WK/UniKA/04-11-2002).
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
         IF (U12INT) THEN
C           Space for [T1,r12] integrals (WK/UniKA/04-11-2002).
            LECOEF = 2*LECOEF
            LXA    =   NPP12
         END IF
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
      IF (GCON12) LCPINT = NPP12
C
      LFIRST = NTUV34
C
C     KC2MAX and LXA have been added (WK/UniKA/04-11-2002).
      LTOTAL = LCCONT + KC2MAX
     &       + MAX(LCSINT + LFCSNT + LFCINT,
     &             LECOEF + 2*LXA + MAX(LHPI + LP1 + LP2 + 2*LXA,
     &             LEUV + LETUV + LHCPRM + LFIRST + LCPINT ) )
      IF (LTOTAL.GT.LWORK) CALL STOPIT('CR1DRV','CR1TWO',KLAST,LWORK)
C
      KODDKC = 1
      KCCONT = KODDKC + KC2MAX
      KCSINT = KCCONT + LCCONT
      KFCSNT = KCSINT + LCSINT
      KFCINT = KFCSNT + LFCSNT
C
C     KAOP and KBOP are addresses for R12 method (WK/UniKA/04-11-2002).
      KAOP   = KCCONT + LCCONT
      KBOP   = KAOP   + LXA
      KECOEF = KBOP   + LXA
C
      KHPI   = KECOEF + LECOEF
      KP1    = KHPI   + LHPI
      KP2    = KP1    + LP1
C
C     KXA and KXB are addresses for R12 method (WK/UniKA/04-11-2002).
      KXA    = KP2    + LP2
      KXB    = KXA    + LXA
C
      KEUV   = KECOEF + LECOEF
      KETUV  = KEUV   + LEUV
      KHCPRM = KETUV  + LETUV
      KFIRST = KHCPRM + LHCPRM
      KCPINT = KFIRST + LFIRST 
C
      IF (U12INT) THEN
C        CALL CR1TWO with call to CR1U12 (WK/UniKA/04-11-2002).
         CALL CR1TWO(HERINT,HERR12,INDHER,INDHER,IODDHH,IODDHC,INDHSQ,
     &               WORK(KECOEF),WORK(KEUV),WORK(KETUV),
     &               LMNPWR,IPNTUV,WORK(KODDKC),
     &               COOR12,EXP12,CSQ,CCFBT,
     &               WORK(KHPI),WORK(KP1),WORK(KP2),
     &               HCINT,WORK(KFCINT),WORK(KHCPRM),WORK(KFIRST),
     &               WORK(KCCONT),WORK(KCSINT),WORK(KFCSNT),
     &               WORK(KCPINT),IPRINT,
     &               WORK(KXA),WORK(KXB),WORK(KAOP),WORK(KBOP),CR1U12)
      ELSE
C        CALL CR1TWO with call to CR1UND (WK/UniKA/04-11-2002).
         CALL CR1TWO(HERINT,HERR12,INDHER,INDHER,IODDHH,IODDHC,INDHSQ,
     &               WORK(KECOEF),WORK(KEUV),WORK(KETUV),
     &               LMNPWR,IPNTUV,WORK(KODDKC),
     &               COOR12,EXP12,CSQ,CCFBT,
     &               WORK(KHPI),WORK(KP1),WORK(KP2),
     &               HCINT,WORK(KFCINT),WORK(KHCPRM),WORK(KFIRST),
     &               WORK(KCCONT),WORK(KCSINT),WORK(KFCSNT),
     &               WORK(KCPINT),IPRINT,
     &               WORK(KXA),WORK(KXB),WORK(KAOP),WORK(KBOP),CR1UND)
      END IF
      RETURN
      END
C  /* Deck cr1two */
      SUBROUTINE CR1TWO(HERINT,HERR12,INDHER,INDHVC,IODDHH,IODDHC,
     &                  INDHSQ,ECOEF,EUV,ETUV,LMNPWR,IPNTUV,IODDKC,
     &                  COOR12,EXP12,CSQ,CCFBT,HPI,P1,P2,
     &                  HCINT,FCINT,HCPRIM,FRSTUV,HCCONT,
     &                  HCSINT,FCSINT,CPINT,IPRINT,
     &                  TEXPA,TEXPB,AOVERP,BOVERP,CR1SUB)
C     TEXPA,TEXPB,AOVERP,BOVERP and CR1SUB have been added (WK/UniKA/04-11-2002).
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
C
      EXTERNAL CR1SUB
      INTEGER TUV
      LOGICAL CONTRA, FRSTUV(NTUV34),
     &        FCINT(NTUV34,KHKT12), FCSINT(NTUV34,KHKT2)
C
      DIMENSION HERINT(NPP12,NPRF34,NTUV),HERR12(NPP12,NPRF34,NTUV),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP), INDHVC(0:*),
     &          IODDHH(NRTOP), IODDHC(NRTOP), INDHSQ(NRTOP),
C
     &          HCPRIM(NPP12,NPRF34,NTUV34),
     &          HCCONT(NCCPP*NTUV34,KCKT12),
     &          HCINT(NCCPP*NTUV34,KHKT12),
C
     &          ECOEF(NPP12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2,3),
     &          ETUV(NPP12), EUV(NPP12),
     &          LMNPWR(KCKMAX,NHKMAX,3), IPNTUV(KC2MAX,0:NRDER,2),
     &          COOR12(NPP12,3), EXP12(NPP12,3), 
     &          HPI(NPP12), P1(NPP12,3), P2(NPP12,3),
     &          IODDKC(KC2MAX), CSQ(NCSQ1,NCSQ2), CCFBT(*), CPINT(*),
     &          HCSINT(NCCPP,NTUV34,KHKT2),
     &          TEXPA(*),TEXPB(*),AOVERP(*),BOVERP(*)
C
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
#include <r12int.h>
      LOGICAL LOCMIA
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from CR1TWO','*',103)
C
C     *********************************
C     ***** Special Case: (ss|xy) *****
C     *********************************
C
      CONTRA = GCON12 .OR. NPRF12.NE.1
C
      IF (KHKT12 .EQ. 1 .AND. .NOT. (U12INT)) THEN
C        No symmetry for [T1,r12] integrals (WK/UniKA/04-11-2002).
         DO TUV = 1, NTUV34
            FRSTUV(TUV) = IODDHH(TUV) .NE. 0
         END DO
C        if (.not.contra) call quit('fix .not.contra...')
         CALL ERICT1(HERINT,HCINT,CCFBT(ICMAT1),CCFBT(ICMAT2),CPINT,
     &               IODDHC,0,FRSTUV,IPRINT)
C
C     *********************************
C     ***** General Case: (xy|zw) *****
C     *********************************
C
      ELSE
        LOCMIA=MIAU12
        MIAU12=.false.
c       the variable miau12, needed for the computation of the
c       (ab|[T1,r12]|cd) integral derivatives
c           
c           ((a-b)/(a+b) d/dA (ab|1/r12|cd) and
c            (a-b)/(a+b) d/dB (ab|1/r12|cd) integrals)
c
c       must be set to .false. when computing the first electron
c       coefficients transformation, i.e. before the first call
c       to subroutine expcft.
c       The expansion coefficients are scaled more conveniently
c       at the second call of subroutine expcft for the (a-b)/(a+b)
c       factor (call from subroutine cr2drv, computing second
c       electron transformation)
c
c       (C. Villani,   Feb 2004)
c
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
     &                  ISCAL1,ISCAL2,IPRINT,
     &                  TEXPA,TEXPB,AOVERP,BOVERP)
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
C
                  IODDKC(ICMP12) = IODDHC(INDHER(L1+L2,M1+M2,N1+N2))
C
C                 contracted spherical harmonics
C                 ==============================
C
                  IF (CONTRA .AND. SPHR12) THEN
                     CALL CR1SUB(HERINT,HERR12,HCPRIM,L1,M1,N1,L2,M2,N2,
     &                           INDHER,INDHVC,IODDHH,INDHSQ,
     &                           ECOEF,EUV,ETUV,FRSTUV,ICOMP1,ICOMP2,
     &                           AOVERP,IPRINT)
                     CALL ERICT1(HCPRIM,HCCONT(1,ICMP12),
     &                           CCFBT(ICMAT1),CCFBT(ICMAT2),CPINT,
     &                           IODDHC,IODDKC(ICMP12),FRSTUV,IPRINT)
C
C                 contracted Cartesians
C                 =====================
C
                  ELSE IF (CONTRA) THEN
                     CALL CR1SUB(HERINT,HERR12,HCPRIM,L1,M1,N1,L2,M2,N2,
     &                           INDHER,INDHVC,IODDHH,INDHSQ,
     &                           ECOEF,EUV,ETUV,FRSTUV,ICOMP1,ICOMP2,
     &                           AOVERP,IPRINT)
                     CALL ERICT1(HCPRIM,HCINT(1,ICMP12),
     &                           CCFBT(ICMAT1),CCFBT(ICMAT2),CPINT,
     &                           IODDHC,IODDKC(ICMP12),FRSTUV,IPRINT)
C
C                 primitive spherical harmonics
C                 =============================
C
                  ELSE IF (SPHR12) THEN
                     CALL CR1SUB(HERINT,HERR12,HCCONT(1,ICMP12),
     &                           L1,M1,N1,L2,M2,N2,
     &                           INDHER,INDHVC,IODDHH,INDHSQ,
     &                           ECOEF,EUV,ETUV,FRSTUV,ICOMP1,ICOMP2,
     &                           AOVERP,IPRINT)
                     CALL ERICT1(HCCONT(1,ICMP12),HCCONT(1,ICMP12),
     &                           CCFBT(ICMAT1),CCFBT(ICMAT2),CPINT,
     &                           IODDHC,IODDKC(ICMP12),FRSTUV,IPRINT)
C
C                 primitive Cartesians
C                 ====================
C
                  ELSE
                     CALL CR1SUB(HERINT,HERR12,HCINT(1,ICMP12),
     &                           L1,M1,N1,L2,M2,N2,
     &                           INDHER,INDHVC,IODDHH,INDHSQ,
     &                           ECOEF,EUV,ETUV,FRSTUV,ICOMP1,ICOMP2,
     &                           AOVERP,IPRINT)
                     CALL ERICT1(HCINT(1,ICMP12),HCINT(1,ICMP12),
     &                           CCFBT(ICMAT1),CCFBT(ICMAT2),CPINT,
     &                           IODDHC,IODDKC(ICMP12),FRSTUV,IPRINT)
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
c       giving miau12 its old value back (C. Villani, Uni-Ka, Feb 2005)
        MIAU12=LOCMIA
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
      DIMENSION HCINT(NCTF12,NPQBCX,NPRF34,NTUV34,KHKT12,NHCTYP),
     &          IODDHC(NRTOP), IPNTUV(KC2MAX,0:NRDER,2)
C
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
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
     &                        1,NCTF12*NPQBCX,1,NPRF34,
     &                        NCTF12*NPQBCX,NPRF34,
     &                        1,LUPRI)
               END IF
  400          CONTINUE
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck cr1und */
      SUBROUTINE CR1UND(HERINT,HERR12,HCPRIM,L1,M1,N1,L2,M2,N2,
     &                  INDHER,INDHVC,IODDHH,INDHSQ,
     &                  ECOEF,EUV,ETUV,FRSTUV,
     &                  ICOMP1,ICOMP2,AOVERP,IPRINT)
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
      DIMENSION HERINT(NPP12,NPRF34,NTUV),HERR12(NPP12,NPRF34,NTUV),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP), INDHVC(0:*),
     &          IODDHH(NRTOP), INDHSQ(NRTOP),
     &          ECOEF(NPP12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2,3),
     &          ETUV(NPP12), EUV(NPP12),
     &          HCPRIM(NPP12,NPRF34,NTUV34),AOVERP(NPP12)
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
         WRITE(LUPRI,'(/,1X,A,2I5/)')' ICOMP1, ICOMP2', ICOMP1,ICOMP2
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
C  /* Deck erict1 */
      SUBROUTINE ERICT1(PPINT,CCINT,CMAT1,CMAT2,CPINT,IODDHC,IODKC,
     &                  FRSTUV,IPRINT)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      INTEGER TUV
      LOGICAL FRSTUV(NTUV34)
      DIMENSION CMAT1(NPRF1,NCTF1)
      DIMENSION CMAT2(NPRF2,NCTF2)
      DIMENSION PPINT(NPQBCX,NPRF12,NPRF34,NTUV34)
      DIMENSION CPINT(NCTF1,NPQBCX,NPRF2)
      DIMENSION CCINT(NCTF12*NPQBCX,NPRF34,NTUV34)
      DIMENSION IODDHC(NRTOP)
#include <hertop.h>
#include <ericom.h>
#include <eriao.h>
C
      IF (IPRINT .GT. 25) THEN
         CALL HEADER('Output from ERICT1',-1)
         CALL HEADER('CMAT1 in ERICT1',-1)
         WRITE (LUPRI,'(2X,A,2I5)') 'NPRF1,NCTF1 ',NPRF1,NCTF1
         CALL OUTPUT(CMAT1,1,NPRF1,1,NCTF1,NPRF1,NCTF1,1,LUPRI)
         CALL HEADER('CMAT2 in ERICT1',-1)
         WRITE (LUPRI,'(2X,A,2I5)') 'NPRF2,NCTF2 ',NPRF2,NCTF2
         CALL OUTPUT(CMAT2,1,NPRF2,1,NCTF2,NPRF2,NCTF2,1,LUPRI)
         WRITE (LUPRI,'(2X,A,2I5)') 'NPRF34,NTUV34',NPRF34,NTUV34
      END IF
C
      IF (NPRF12.EQ.1 .AND. .NOT.GCON12) THEN
         DO 100 TUV = 1, NTUV34
         IF (IODDHC(TUV) .EQ. IODKC) THEN
            IF (FRSTUV(TUV)) THEN
               DO K = 1, NPRF34 
                 DO J = 1, NPQBCX
                  CCINT(J,K,TUV) = D0
                 END DO 
               END DO
            ELSE
               DO K = 1, NPRF34
                 DO J = 1, NPQBCX
                  CCINT(J,K,TUV) = PPINT(J,1,K,TUV)
                 END DO
               END DO
            END IF
         END IF
  100    CONTINUE
      ELSE
         DO 200 TUV = 1, NTUV34
         IF (IODDHC(TUV) .EQ. IODKC) THEN
            IF (FRSTUV(TUV)) THEN
               DO K = 1, NPRF34 
                 DO I = 1, NCTF12*NPQBCX
                  CCINT(I,K,TUV) = D0
                 END DO 
               END DO
            ELSE
               IF (GCON12) THEN
                  CALL CONT12(PPINT(1,1,1,TUV),CCINT(1,1,TUV),
     &                 CMAT1,CMAT2,CPINT)
               ELSE
                  DO K = 1, NPRF34
                     DO J = 1, NPQBCX
                        CCINT(J,K,TUV) = PPINT(J,1,K,TUV)
     &                                 + PPINT(J,2,K,TUV)
                     END DO
                     DO I = 3, NPRF12 
                       DO J = 1, NPQBCX
                        CCINT(J,K,TUV) = CCINT(J,K,TUV) 
     &                                 + PPINT(J,I,K,TUV)
                       END DO 
                     END DO
                  END DO
               END IF
            END IF
         END IF
  200    CONTINUE
      END IF
C
      IF (IPRINT .GT. 25) THEN
         CALL HEADER('Output from CR1SEG',-1)
         WRITE (LUPRI,'(2X,A,2I5)') 'NPRF1, NCTF1 ',NPRF1, NCTF1
         WRITE (LUPRI,'(2X,A,2I5)') 'NPRF2, NCTF2 ',NPRF2, NCTF2
         WRITE (LUPRI,'(2X,A,2I5)') 'NPRF34,NTUV34',NPRF34,NTUV34
         WRITE (LUPRI,'(2X,A,2I5)') 'NPQBCX,NPRF12',NPQBCX,NPRF12
         WRITE (LUPRI,'(2X,A,3L5)') 'GCON1, GCON2,',GCON1, GCON2, GCON12
         IF (GCON1) THEN
            CALL HEADER('CMAT1 in ERICT1',-1)
            CALL OUTPUT(CMAT1,1,NPRF1,1,NCTF1,NPRF1,NCTF1,1,LUPRI)
         END IF
         IF (GCON2) THEN
            CALL HEADER('CMAT2 in ERICT1',-1)
            CALL OUTPUT(CMAT2,1,NPRF2,1,NCTF2,NPRF2,NCTF2,1,LUPRI)
         END IF
         DO TUV = 1, NTUV34 
           IF (IODDHC(TUV) .EQ. IODKC) THEN
            IF (.NOT.FRSTUV(TUV)) THEN
               CALL HEADER('PPINT in CR1SEG',-1)
               WRITE (LUPRI,'(2X,A,I5)') 'TUV: ',TUV
               CALL OUTPUT(PPINT(1,1,1,TUV),1,NPQBCX*NPRF12,1,NPRF34,
     &                     NPQBCX*NPRF12,NPRF34,1,LUPRI)
            END IF
            CALL HEADER('CCINT in CR1SEG',-1)
            WRITE (LUPRI,'(2X,A,I5)') 'TUV: ',TUV
            CALL OUTPUT(CCINT(1,1,TUV),1,NPQBCX,1,NPRF34,
     &                  NPQBCX,NPRF34,1,LUPRI)
           END IF 
         END DO
      END IF
C
      RETURN
      END

C  /* Deck cont12 */
      SUBROUTINE CONT12(PPINT,CCINT,CMAT1,CMAT2,CPINT)
#include <implicit.h>
      DIMENSION CMAT1(NPRF1,NCTF1)
      DIMENSION CMAT2(NPRF2,NCTF2)
      DIMENSION PPINT(NPQBCX*NPRF2,NPRF1,NPRF34)
      DIMENSION CPINT(NCTF1,NPQBCX,NPRF2)
      DIMENSION CCINT(NCTF2,NCTF1*NPQBCX,NPRF34)
#include <ericom.h>
#include <eriao.h>
      DO L = 1, NPRF34
         CALL DZERO(CPINT,NCTF1*NPQBCX*NPRF2)
         DO 200 K = 1, NPRF1
            DO 200 M = 1, NCTF1
               DO 200 N = 1, NPQBCX*NPRF2
C              ... "dirty trick" for CPINT to get long loop length
C                  on last two indices (will go out of bounds on 2. index)
                  CPINT(M,N,1) = CPINT(M,N,1) +
     &                 CMAT1(K,M)*PPINT(N,K,L)
 200     CONTINUE
         CALL DZERO(CCINT(1,1,L),NCTF2*NCTF1*NPQBCX)
         DO 300 K = 1, NPRF2
            DO 300 M = 1, NCTF2
               DO 300 N = 1, NCTF1*NPQBCX
C              ... "dirty trick" for CPINT to get long loop length
C                  on first two indices (will go out of bounds on 1. index)
                  CCINT(M,N,L) = CCINT(M,N,L) + 
     &                 CMAT2(K,M)*CPINT(N,1,K)
 300     CONTINUE
      END DO
      RETURN
      END
c     Part for [T1,r12] integral derivatives, based on subroutine CR1DRV
      SUBROUTINE CR1MIA(HERINT,HERR12,HCINT,INDHER,IODDHH,IODDHC,INDHSQ,
     &                  LMNPWR,IPNTUV,COOR12,EXP12,CSQ,CCFBT,
     &                  WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
      LOGICAL DOREC, FRSTAO, CONTRA
      DIMENSION HERINT(*), HCINT(*), HERR12(*),
     &          INDHER(*), INDHSQ(*), IODDHH(*), IODDHC(*),
     &          LMNPWR(*), IPNTUV(*), 
     &          COOR12(NPP12,3), EXP12(NPP12,3), CSQ(*), CCFBT(*),
     &          WORK(LWORK)
      EXTERNAL  CMIASB
#include <cbieri.h>
#include <ericom.h>
#include <eriao.h>
#include <r12int.h>
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
C     HCCONT  | HCSINT | LFCSNT | LFCINT
C             |
C             | ECOEF | HPI | P1 | P2
C                     |
C                     | EUV | ETUV | HCPRIM
C
      LCCONT = 0
      LCSINT = 0
      LFCSNT = 0
      LFCINT = 0
      LECOEF = 0
      LHPI   = 0
      LP1    = 0
      LP2    = 0
C     Space for R12 method (WK/UniKA/04-11-2002).
      LXA    = 0
      LEUV   = 0
      LETUV  = 0
      LHCPRM = 0
      LCPINT = 0
C
c     IF (KHKT12 .GT. 1 .OR. U12INT) THEN
C        No symmetry for [T1,r12] integrals WK/UniKA/04-11-2002).
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
c        IF (U12INT) THEN
C           Space for [T1,r12] integrals (WK/UniKA/04-11-2002).
c           if(MAXDER.eq.0) then
c             LECOEF = 2*LECOEF
c           else
c             Additional space for E(tij);2 coefficients,
c             needed for [T1,r12] derivatives
c             (C. Villani,   Winter 2004)
              LECOEF = 3*LECOEF
c           end if
            LXA    =   NPP12
c        END IF
         LHPI   =   NPP12
         LP1    = 3*NPP12
         LP2    = 3*NPP12
C
         LEUV   = NPP12
         LETUV  = NPP12
C
         LHCPRM = NPPPP*NTUV34
c     END IF
C
      IF (GCON12) LCPINT = NPP12
C
      LFIRST = NTUV34
C
C     KC2MAX and LXA have been added (WK/UniKA/04-11-2002).
c     if(.not.u12int.or.maxder.eq.0) then
c       LTOTAL = LCCONT + KC2MAX
c    &         + MAX(LCSINT + LFCSNT + LFCINT,
c    &               LECOEF + 2*LXA + MAX(LHPI + LP1 + LP2 + 2*LXA,
c    &               LEUV + LETUV + LHCPRM + LFIRST + LCPINT ) )
c     else
c       Additional space needed for E(tij);2 coefficients,
c       needed for [T1,r12] derivatives   [To be checked!]
c       (C. Villani,   Winter 2004)
        LTOTAL = LCCONT + KC2MAX
     &         + MAX(LCSINT + LFCSNT + LFCINT,
     &               LECOEF + 3*LXA + MAX(LHPI + LP1 + LP2 + 3*LXA,
     &               LEUV + LETUV + LHCPRM + LFIRST + LCPINT ) )
c     end if
      IF (LTOTAL.GT.LWORK) CALL STOPIT('CR1DRV','CR1TWO',KLAST,LWORK)
C
      KODDKC = 1
      KCCONT = KODDKC + KC2MAX
      KCSINT = KCCONT + LCCONT
      KFCSNT = KCSINT + LCSINT
      KFCINT = KFCSNT + LFCSNT
C
C     KAOP and KBOP are addresses for R12 method (WK/UniKA/04-11-2002).
      KAOP   = KCCONT + LCCONT
      KBOP   = KAOP   + LXA
      KECOEF = KBOP   + LXA
C
      KHPI   = KECOEF + LECOEF
      KP1    = KHPI   + LHPI
      KP2    = KP1    + LP1
C
C     KXA and KXB are addresses for R12 method (WK/UniKA/04-11-2002).
      KXA    = KP2    + LP2
      KXB    = KXA    + LXA
C
      KEUV   = KECOEF + LECOEF
      KETUV  = KEUV   + LEUV
      KHCPRM = KETUV  + LETUV
      KFIRST = KHCPRM + LHCPRM
      KCPINT = KFIRST + LFIRST 
C
C     CALL CR1TWO with call to CMIASB (C. Villani Uni-Ka Feb 2005).
      CALL CR1TWO(HERINT,HERR12,INDHER,INDHER,IODDHH,IODDHC,INDHSQ,
     &            WORK(KECOEF),WORK(KEUV),WORK(KETUV),
     &            LMNPWR,IPNTUV,WORK(KODDKC),
     &            COOR12,EXP12,CSQ,CCFBT,
     &            WORK(KHPI),WORK(KP1),WORK(KP2),
     &            HCINT,WORK(KFCINT),WORK(KHCPRM),WORK(KFIRST),
     &            WORK(KCCONT),WORK(KCSINT),WORK(KFCSNT),
     &            WORK(KCPINT),IPRINT,
     &            WORK(KXA),WORK(KXB),WORK(KAOP),WORK(KBOP),CMIASB)
      RETURN
      END
      SUBROUTINE CMIASB(HERINT,HERR12,HCPRIM,L1,M1,N1,L2,M2,N2,
     &                  INDHER,INDHVC,IODDHH,INDHSQ,
     &                  ECOEF,EUV,ETUV,FRSTUV,
     &                  ICOMP1,ICOMP2,ABVERP,IPRINT)
C
C     This subroutine is a variant of CR1U12 for derivative integrals
c     over [T1,r12].
C     Written by Cristian Villani (University of Karlsruhe, winter 2004).
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
C
      PARAMETER (DP5 = 0.5D0)
      INTEGER T, U, V, TUV
      LOGICAL FRSTUV(NTUV34)
      INTEGER MIAXYZ(3)
C
      DIMENSION HERINT(NPP12,NPRF34,NTUV),HERR12(NPP12,NPRF34,NTUV),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP), INDHVC(0:*),
     &          IODDHH(NRTOP), INDHSQ(NRTOP),
     &          ECOEF(NPP12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2,3,3),
     &          ETUV(NPP12), EUV(NPP12),
     &          HCPRIM(NPP12,NPRF34,NTUV34),ABVERP(NPP12)
C
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
#include <r12int.h>
#include <ibtfun.h>

      DO I=1,3
        MIAXYZ(I)=0
      END DO
      MIAXYZ(XMIADR)=1
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from CMIASB','*',103)
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
      MINTP = IBTAND(MAXT+1,I120(1))
      MINUP = IBTAND(MAXU+1,I120(2))
      MINVP = IBTAND(MAXV+1,I120(3))
      MINTP1 = IBTAND(MAXT+MIAXYZ(1)+1,I120(1))
      MINTP2 = IBTAND(MAXT+MIAXYZ(1),I120(1))
      MINUP1 = IBTAND(MAXU+MIAXYZ(2)+1,I120(2))
      MINUP2 = IBTAND(MAXU+MIAXYZ(2),I120(2))
      MINVP1 = IBTAND(MAXV+MIAXYZ(3)+1,I120(3))
      MINVP2 = IBTAND(MAXV+MIAXYZ(3),I120(3))
C
      IF (IPRINT .GT. 25) THEN
         WRITE(LUPRI,'(/,1X,A,2I5/)')' ICOMP1, ICOMP2', ICOMP1,ICOMP2
         WRITE(LUPRI,'(1X,A,15X,3I5)')' T loop:',MINT,MAXT,INCT
         WRITE(LUPRI,'(1X,A,15X,3I5)')' U loop:',MINU,MAXU,INCU
         WRITE(LUPRI,'(1X,A,15X,3I5)')' V loop:',MINV,MAXV,INCV
      END IF
C
      DO 100 TUV = 1, NTUV34
         FRSTUV(TUV) = .TRUE.
  100 CONTINUE
C
C     Expansion coefficients := R(Rx*Ey*Ez)
C     ======================================
C
      DO 201 V = MINVP2, MAXV, INCV
      DO 201 U = MINUP2, MAXU, INCU
         DO 211 I = 1, NPP12 
            EUV(I) = ECOEF(I,V,N1,N2,3,MIAXYZ(3)+1) 
     &             * ECOEF(I,U,M1,M2,2,MIAXYZ(2)+1)
  211    CONTINUE
         DO 301 T = MINTP1, MAXT, INCT
            DO 311 I = 1, NPP12
               ETUV(I) = ECOEF(I,T,L1,L2,1,MIAXYZ(1)+2) * EUV(I)
  311       CONTINUE
            ITUV = INDHER(T+1,U,V)
            INDS = INDHSQ(ITUV)
            DO 401 TUV = 1, NTUV34
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
                  DO 501 J = 1, NPRF34
                  DO 501 I = 1, NPP12
                     HCPRIM(I,J,TUV) = ETUV(I)*HERR12(I,J,INDT)
  501             CONTINUE
               ELSE
                  DO 601 J = 1, NPRF34
                  DO 601 I = 1, NPP12
                     HCPRIM(I,J,TUV) = HCPRIM(I,J,TUV)
     &                       + ETUV(I)*HERR12(I,J,INDT)
  601             CONTINUE
               END IF
            END IF
  401       CONTINUE
  301    CONTINUE
  201 CONTINUE
C
C     Expansion coefficients := (A/B)/p * Rx*Ey*Ez
C     =============================================
C
      DO V = MINV, MAXV, INCV
        DO U = MINU, MAXU, INCU
          DO I = 1, NPP12
            EUV(I) = ECOEF(I,V,N1,N2,3,1) * ECOEF(I,U,M1,M2,2,1)
          END DO
          DO T = MINTP, MAXT, INCT
            DO I = 1, NPP12
               ETUV(I) = ABVERP(I)*ECOEF(I,T,L1,L2,1,2) * EUV(I)
            END DO
            ITUV = INDHER(MIAXYZ(1)+T+1,MIAXYZ(2)+U,MIAXYZ(3)+V)
            INDS = INDHSQ(ITUV)
            DO TUV = 1, NTUV34
              IF(IODDHH(ITUV) .EQ. IODDHH(TUV)) THEN
#if defined (SYS_AIX)
C               code due to AIX xlf version 2.2 bug
                INDT = INDS + INDHSQ(TUV)
                INDT = INDHVC(INDT)
#else
                INDT = INDHVC(INDS + INDHSQ(TUV))
#endif
                IF(FRSTUV(TUV)) THEN
                  FRSTUV(TUV) = .FALSE.
                  DO J = 1, NPRF34
                    DO I = 1, NPP12
                      HCPRIM(I,J,TUV) = ETUV(I)*HERR12(I,J,INDT)
                    END DO
                  END DO
                ELSE
                  DO J = 1, NPRF34
                    DO I = 1, NPP12
                      HCPRIM(I,J,TUV) = HCPRIM(I,J,TUV)
     &                                + ETUV(I)*HERR12(I,J,INDT)
                    END DO
                  END DO
                END IF
              END IF
            END DO
          END DO
        END DO
      END DO
C
C     Expansion coefficients := R(Ex*Ry*Ez)
C     ======================================
C
      DO 202 V = MINVP2, MAXV, INCV
      DO 202 U = MINUP1, MAXU, INCU
         DO 212 I = 1, NPP12 
            EUV(I) = ECOEF(I,V,N1,N2,3,MIAXYZ(3)+1)
     &             * ECOEF(I,U,M1,M2,2,MIAXYZ(2)+2)
  212    CONTINUE
         DO 302 T = MINTP2, MAXT, INCT
            DO 312 I = 1, NPP12
               ETUV(I) = ECOEF(I,T,L1,L2,1,MIAXYZ(1)+1) * EUV(I)
  312       CONTINUE
            ITUV = INDHER(T,U+1,V)
            INDS = INDHSQ(ITUV)
            DO 402 TUV = 1, NTUV34
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
                  DO 502 J = 1, NPRF34
                  DO 502 I = 1, NPP12
                     HCPRIM(I,J,TUV) = ETUV(I)*HERR12(I,J,INDT)
  502             CONTINUE
               ELSE
                  DO 602 J = 1, NPRF34
                  DO 602 I = 1, NPP12
                     HCPRIM(I,J,TUV) = HCPRIM(I,J,TUV)
     &                       + ETUV(I)*HERR12(I,J,INDT)
  602             CONTINUE
               END IF
            END IF
  402       CONTINUE
  302    CONTINUE
  202 CONTINUE
C
C     Expansion coefficients := (A/B)/p * Ex*Ry*Ez
C     =============================================
C
      DO V = MINV, MAXV, INCV
      DO U = MINUP, MAXU, INCU
         DO I = 1, NPP12
            EUV(I) = ECOEF(I,V,N1,N2,3,1) * ECOEF(I,U,M1,M2,2,2)
         END DO
         DO T = MINT, MAXT, INCT
            DO I = 1, NPP12
               ETUV(I) = ABVERP(I) * ECOEF(I,T,L1,L2,1,1) * EUV(I)
            END DO
            ITUV = INDHER(MIAXYZ(1)+T,MIAXYZ(2)+U+1,MIAXYZ(3)+V)
            INDS = INDHSQ(ITUV)
            DO TUV = 1, NTUV34
              IF (IODDHH(ITUV) .EQ. IODDHH(TUV)) THEN
#if defined (SYS_AIX)
C               code due to AIX xlf version 2.2 bug
                INDT = INDS + INDHSQ(TUV)
                INDT = INDHVC(INDT)
#else
                INDT = INDHVC(INDS + INDHSQ(TUV))
#endif
                IF(FRSTUV(TUV)) THEN
                  FRSTUV(TUV) = .FALSE.
                  DO J = 1, NPRF34
                    DO I = 1, NPP12
                      HCPRIM(I,J,TUV) = ETUV(I)*HERR12(I,J,INDT)
                    END DO
                  END DO
                ELSE
                  DO J = 1, NPRF34
                    DO I = 1, NPP12
                      HCPRIM(I,J,TUV) = HCPRIM(I,J,TUV)
     &                                + ETUV(I)*HERR12(I,J,INDT)
                    END DO
                  END DO
                END IF
              END IF
            END DO
          END DO
        END DO
      END DO
C
C     Expansion coefficients := R(Ex*Ey*Rz)
C     ==================================
C
      DO 203 V = MINVP1, MAXV, INCV
      DO 203 U = MINUP2, MAXU, INCU
         DO 213 I = 1, NPP12 
            EUV(I) = ECOEF(I,V,N1,N2,3,MIAXYZ(3)+2) *
     .               ECOEF(I,U,M1,M2,2,MIAXYZ(2)+1)
  213    CONTINUE
         DO 303 T = MINTP2, MAXT, INCT
            DO 313 I = 1, NPP12
               ETUV(I) = ECOEF(I,T,L1,L2,1,MIAXYZ(1)+1) * EUV(I)
  313       CONTINUE
            ITUV = INDHER(T,U,V+1)
            INDS = INDHSQ(ITUV)
            DO 403 TUV = 1, NTUV34
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
                  DO 503 J = 1, NPRF34
                  DO 503 I = 1, NPP12
                     HCPRIM(I,J,TUV) = ETUV(I)*HERR12(I,J,INDT)
  503             CONTINUE
               ELSE
                  DO 603 J = 1, NPRF34
                  DO 603 I = 1, NPP12
                     HCPRIM(I,J,TUV) = HCPRIM(I,J,TUV)
     &                       + ETUV(I)*HERR12(I,J,INDT)
  603             CONTINUE
               END IF
            END IF
  403       CONTINUE
  303    CONTINUE
  203 CONTINUE
C
C     Expansion coefficients := (A/B)/p * (Ex*Ey*Rz)
C     ==================================
C
      DO V = MINVP, MAXV, INCV
        DO U = MINU, MAXU, INCU
          DO I = 1, NPP12 
            EUV(I) = ECOEF(I,V,N1,N2,3,2) * ECOEF(I,U,M1,M2,2,1)
          END DO
          DO T = MINT, MAXT, INCT
            DO I = 1, NPP12
               ETUV(I) = ABVERP(I) * ECOEF(I,T,L1,L2,1,1) * EUV(I)
            END DO
C
            ITUV = INDHER(MIAXYZ(1)+T,MIAXYZ(2)+U,MIAXYZ(3)+V+1)
            INDS = INDHSQ(ITUV)
            DO TUV = 1, NTUV34
              IF (IODDHH(ITUV) .EQ. IODDHH(TUV)) THEN
#if defined (SYS_AIX)
C               code due to AIX xlf version 2.2 bug
                INDT = INDS + INDHSQ(TUV)
                INDT = INDHVC(INDT)
#else
                INDT = INDHVC(INDS + INDHSQ(TUV))
#endif
                IF (FRSTUV(TUV)) THEN
                  FRSTUV(TUV) = .FALSE.
                  DO J = 1, NPRF34
                    DO I = 1, NPP12
                      HCPRIM(I,J,TUV) = ETUV(I)*HERR12(I,J,INDT)
                    END DO
                  END DO
                ELSE
                  DO J = 1, NPRF34
                    DO I = 1, NPP12
                      HCPRIM(I,J,TUV) = HCPRIM(I,J,TUV)
     &                                + ETUV(I)*HERR12(I,J,INDT)
                    END DO
                  END DO
                END IF
              END IF
            END DO
          END DO
        END DO
      END DO

      RETURN
      END
C --- end of eri2car1.F ---
