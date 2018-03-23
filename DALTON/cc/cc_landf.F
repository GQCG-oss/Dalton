!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2016 (2015), see http://daltonprogram.org"
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
C
C  /* Deck cc_FHF   */
      SUBROUTINE CC_FHF(C1AM,RHO1,ISYMV,WORK,LWORK)
C
C     Written by Ove Christiansen june 1996 
C
C     Purpose: Calculate <HF| [[H,tau,mu],tau,ny]|HF> contribution
C              to contraction of F-matrix with one vector.
C
C     RHO1(ai)  = sum(bj)[c1am(bj)*L(ia,jb)]
C                                  (is packed ai,bj)
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccinftap.h"
#include "ccsdsym.h"
#include "cclr.h"
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION RHO1(*),C1AM(*),WORK(LWORK)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC_FHF')
C
      KI = 1
C
      IF (LWORK .LT. NT2AM(ISYMOP)) 
     * CALL QUIT('Too little workspace in CC_FHF ')
C
      REWIND(LUIAJB)
      CALL READI(LUIAJB,IRAT*NT2AM(ISYMOP),WORK(KI))
C
      IF (IPRINT .GT. 40 ) THEN
         RHO2N = DDOT(NT2AM(ISYMOP),WORK(KI),1,WORK(KI),1)
         WRITE(LUPRI,1) 'Norm of Integrals (ia|jb) in CC_FHF',RHO2N
         CALL AROUND( 'In CC_FHF:  Integrals (ia|jb) ' )
         CALL CC_PRP(DUM,WORK(KI),1,0,1)
      ENDIF
C
      IOPTTCME = 1
      CALL CCSD_TCMEPK(WORK(KI),1.0D0,ISYMOP,IOPTTCME)
C
      ISYMAI = ISYMV
      ISYMBJ = ISYMV
C
      DO 100 NAI = 1, NT1AM(ISYMAI)
C
         DO 200 NBJ = 1, NT1AM(ISYMBJ)
C
            NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *            + INDEX(NAI,NBJ) 
            RHO1(NAI) = RHO1(NAI) + TWO*WORK(NAIBJ)*C1AM(NBJ)
C
  200    CONTINUE
C
  100 CONTINUE
C
      CALL QEXIT('CC_FHF')
C
   1  FORMAT(1x,A35,1X,E20.10)
      RETURN
      END
C  /* Deck cc_c1t2c */
      SUBROUTINE CC_C1T2C(C1T2C,VLTR1,RTR2,XLAMDP,XLAMDH,XLMPC,XLMHC,
     *                    WORK,LWORK,ISYL,ISYR,IOPT)
C
C     Written by Asger Halkier 14/10 - 1995
C
C     Purpose: Calculate the contraction of the trial
C              vector VLTR1 with packed RTR2-amplitudes and
C              transformation of the result to AO-basis.
C
C     If IOPT equals 1 the whole routine is used, and the result is
C     returned back-transformed to AO-basis, whereas, if IOPT
C     equals 2, then the result is returned directly in MO-basis.
C     ISYL is the symmetry of the vector that is to be contracted
C     with the T2-amplitude-array!
C
C     Ove Christiansen 6-6-1996:
C     Generalised for F-matrix calculation; IOPT=3
C     Debugged for general sym. of L and R.
C
C     Xdl = Sum(bk)L(bk){2R(dl,bk)-R(dk,bl)} (iopt2)
C
C     C1T2C = sum(dl)[Xdl*LamH(beta,d)*LamP(alfa,l)]
C            +sum(bk)[VLTR1(bk)*LampC(alfa,b)*LamH(beta,k)] 
C            +sum(bk)[VLTR1(bk)*Lamp(alfa,b)*LamHC(beta,k)] 
C
C     For IOPT=1 is the third term omitted. 
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XLAMDH(*),XLAMDP(*),XLMPC(*),XLMHC(*)
      DIMENSION C1T2C(*),VLTR1(*),RTR2(*),WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "cclr.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      CALL QENTER('CC_C1T2C')
C
      ISYRES = MULD2H(ISYR,ISYL)
C
C-----------------------------
C     Initialize result array.
C-----------------------------
C
      IF ((IOPT .EQ. 1).OR.(IOPT.EQ.3)) THEN
         CALL DZERO(C1T2C,N2BST(ISYRES))
      ELSE IF (IOPT .EQ. 2) THEN
         CALL DZERO(C1T2C,NT1AM(ISYRES))
      ENDIF
C
      ISYMBK = ISYL  
      ISYMDL = MULD2H(ISYMBK,ISYR)
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      KLDINT = 1
      KEND1  = KLDINT + NT1AM(ISYMDL)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient work space in CC_C1T2C-1')
      ENDIF
C
      CALL DZERO(WORK(KLDINT),NT1AM(ISYMDL))
C
C-------------------------------------------
C     Calculate coulomb part of contraction.
C-------------------------------------------
C
      IF (.NOT. CCS) THEN
       IF (ISYMDL .LT. ISYMBK) THEN
C
         KOFF1  = IT2AM(ISYMDL,ISYMBK) + 1
C
         NTOTDL = MAX(NT1AM(ISYMDL),1)
C
         CALL DGEMV('N',NT1AM(ISYMDL),NT1AM(ISYMBK),TWO,RTR2(KOFF1),
     *              NTOTDL,VLTR1,1,ZERO,WORK(KLDINT),1)
C
       ELSE IF (ISYMDL .GT. ISYMBK) THEN
C
         KOFF1  = IT2AM(ISYMBK,ISYMDL) + 1
C
         NTOTBK = MAX(NT1AM(ISYMBK),1)
C
         CALL DGEMV('T',NT1AM(ISYMBK),NT1AM(ISYMDL),TWO,RTR2(KOFF1),
     *              NTOTBK,VLTR1,1,ZERO,WORK(KLDINT),1)
C
       ELSE IF (ISYMDL .EQ. ISYMBK) THEN
C
         DO 100 ISYML = 1,NSYM
C
            ISYMD = MULD2H(ISYML,ISYMDL)
C
C----------------------------------------------------
C           First intermediate work space allocation.
C----------------------------------------------------
C
            KCOUL  = KEND1
            KENDI1 = KCOUL + NT1AM(ISYMBK)*NVIR(ISYMD)
            LWRKI1 = LWORK - KENDI1
C
            IF (LWRKI1 .LT. 0) THEN
               CALL QUIT('Insufficient work space in CC_C1T2C-2')
            ENDIF
C
            CALL DZERO(WORK(KCOUL),NT1AM(ISYMBK)*NVIR(ISYMD))
C
            DO 110 L = 1,NRHF(ISYML)
C
               DO 120 D = 1,NVIR(ISYMD)
C
                  NDL = IT1AM(ISYMD,ISYML) + NVIR(ISYMD)*(L - 1) + D
C
                  DO 130 NBK = 1,NT1AM(ISYMBK)
C
                     NDBK  = KCOUL + NVIR(ISYMD)*(NBK - 1) + D - 1
                     NBKDL = IT2AM(ISYMBK,ISYMDL) + INDEX(NBK,NDL)
C
                     WORK(NDBK) = RTR2(NBKDL)
C
  130             CONTINUE
  120          CONTINUE
C
               KOFF1 = KLDINT + IT1AM(ISYMD,ISYML)
     *               + NVIR(ISYMD)*(L - 1)
C
               NTOTD = MAX(NVIR(ISYMD),1)
C
               CALL DGEMV('N',NVIR(ISYMD),NT1AM(ISYMBK),TWO,
     *                    WORK(KCOUL),NTOTD,VLTR1,1,ZERO,WORK(KOFF1),1)
C
  110       CONTINUE
  100    CONTINUE
C
       ENDIF
C
C---------------------------------------------
C      Calculate exchange part of contraction.
C---------------------------------------------
C
       DO 140 ISYMDK = 1,NSYM
C
         ISYMBL = MULD2H(ISYMDK,ISYR  )
C
         DO 150 ISYML = 1,NSYM
C
            ISYMB = MULD2H(ISYML,ISYMBL)
            ISYMK = MULD2H(ISYMB,ISYL  )
            ISYMD = MULD2H(ISYMK,ISYMDK)
C
C-----------------------------------------------------
C           Second intermediate work space allocation.
C-----------------------------------------------------
C
            KEXCH  = KEND1
            KENDI2 = KEXCH + NVIR(ISYMD)*NT1AM(ISYL  )
            LWRKI2 = LWORK - KENDI2
C
            IF (LWRKI2 .LT. 0) THEN
               CALL QUIT('Insufficient work space in CC_C1T2C-3')
            ENDIF
C
            CALL DZERO(WORK(KEXCH),NVIR(ISYMD)*NT1AM(ISYL  ))
C
            DO 160 L = 1,NRHF(ISYML)
C
               DO 170 B = 1,NVIR(ISYMB)
C
                  NBL = IT1AM(ISYMB,ISYML) + NVIR(ISYMB)*(L - 1) + B
C
                  IF (ISYMDK .EQ. ISYMBL) THEN
C
                     DO 180 K = 1,NRHF(ISYMK)
C
                        NBK = IT1AM(ISYMB,ISYMK) + NVIR(ISYMB)*(K - 1)
     *                      + B
C
                        DO 190 D = 1,NVIR(ISYMD)
C
                           NDBK  = KEXCH + NVIR(ISYMD)*(NBK - 1) + D - 1
                           NDK   = IT1AM(ISYMD,ISYMK)
     *                           + NVIR(ISYMD)*(K - 1) + D
                           NDKBL = IT2AM(ISYMDK,ISYMBL) + INDEX(NDK,NBL)
C
                           WORK(NDBK) = RTR2(NDKBL)
C
  190                   CONTINUE
  180                CONTINUE
C
                  ELSE IF (ISYMDK .LT. ISYMBL) THEN
C
                     DO 200 K = 1,NRHF(ISYMK)
C
                        NBK = IT1AM(ISYMB,ISYMK) + NVIR(ISYMB)*(K - 1)
     *                      + B
C
                        DO 210 D = 1,NVIR(ISYMD)
C
                           NDBK  = KEXCH + NVIR(ISYMD)*(NBK - 1) + D - 1
                           NDK   = IT1AM(ISYMD,ISYMK)
     *                           + NVIR(ISYMD)*(K - 1) + D
                           NDKBL = IT2AM(ISYMDK,ISYMBL)
     *                           + NT1AM(ISYMDK)*(NBL - 1) + NDK
C
                           WORK(NDBK) = RTR2(NDKBL)
C
  210                   CONTINUE
  200                CONTINUE
C
                  ELSE IF (ISYMDK .GT. ISYMBL) THEN
C
                     DO 220 K = 1,NRHF(ISYMK)
C
                        NBK = IT1AM(ISYMB,ISYMK) + NVIR(ISYMB)*(K - 1)
     *                      + B
C
                        DO 230 D = 1,NVIR(ISYMD)
C
                           NDBK  = KEXCH + NVIR(ISYMD)*(NBK - 1) + D - 1
                           NDK   = IT1AM(ISYMD,ISYMK)
     *                           + NVIR(ISYMD)*(K - 1) + D
                           NDKBL = IT2AM(ISYMBL,ISYMDK)
     *                           + NT1AM(ISYMBL)*(NDK - 1) + NBL
C
                           WORK(NDBK) = RTR2(NDKBL)
C
  230                   CONTINUE
  220                CONTINUE
C
                  ENDIF
C
  170          CONTINUE
C
               KOFF2 = KLDINT + IT1AM(ISYMD,ISYML)
     *               + NVIR(ISYMD)*(L - 1)
C
               NTOTD = MAX(NVIR(ISYMD),1)
C
               CALL DGEMV('N',NVIR(ISYMD),NT1AM(ISYL  ),-ONE,
     *                    WORK(KEXCH),NTOTD,VLTR1,1,ONE,WORK(KOFF2),1)
C
  160       CONTINUE
  150    CONTINUE
  140  CONTINUE
C
      ENDIF
C
      IF ((IOPT .EQ. 1) .OR. (IOPT .EQ. 3)) THEN
C
C===================================
C        Transformation to AO-basis.
C===================================
C
         ISYMDL = MULD2H(ISYL,ISYR)
         ISALBE = ISYMDL
C
C-------------------------------------------------
C        Calculate first transformed contribution.
C-------------------------------------------------
C
         DO 240 ISYMAL = 1,NSYM
C
            ISYMBE = MULD2H(ISYMAL,ISALBE)
            ISYML  = ISYMAL
            ISYVID = ISYMBE
C
C-----------------------------------
C           Work-space allocation 2.
C-----------------------------------
C
            KSCR2 = KEND1
            KEND2 = KSCR2 + NRHF(ISYML)*NBAS(ISYMBE)
            LWRK2 = LWORK - KEND2
C
            IF (LWRK2 .LT. 0) THEN
               CALL QUIT('Insufficient work-space in CC_C1T2C-4')
            ENDIF
C
C-------------------------------------
C           Do the actual calculation.
C-------------------------------------
C
            KOFF1 = ILMVIR(ISYVID) + 1
            KOFF2 = KLDINT + IT1AM(ISYVID,ISYML)
            KOFF3 = KSCR2
C
            NTOTBE = MAX(NBAS(ISYMBE),1)
            NTOVID = MAX(NVIR(ISYVID),1)
C
            CALL DGEMM('N','N',NBAS(ISYMBE),NRHF(ISYML),NVIR(ISYVID),
     *                 ONE,XLAMDH(KOFF1),NTOTBE,WORK(KOFF2),NTOVID,
     *                 ZERO,WORK(KOFF3),NTOTBE)
C
            KOFF1 = ILMRHF(ISYML) + 1
            KOFF2 = KSCR2
            KOFF3 = IAODIS(ISYMAL,ISYMBE) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTBE = MAX(NBAS(ISYMBE),1)
C
            CALL DGEMM('N','T',NBAS(ISYMAL),NBAS(ISYMBE),NRHF(ISYML),
     *                 ONE,XLAMDP(KOFF1),NTOTAL,WORK(KOFF2),NTOTBE,
     *                 ONE,C1T2C(KOFF3),NTOTAL)
C
  240    CONTINUE
C
C--------------------------------------------
C        Add second transformed contribution.
C        First for lampc*lamh
C--------------------------------------------
C
         ISALBE = ISYRES
C
         DO 250 ISYMAL = 1,NSYM
C
            ISYMBE = MULD2H(ISYMAL,ISALBE)
            ISYMB  = MULD2H(ISYMAL,ISYR)
            ISYMK  = ISYMBE
C
C-----------------------------------
C           Work-space allocation 3.
C-----------------------------------
C
            KSCR3 = KEND1
            KEND3 = KSCR3 + NBAS(ISYMAL)*NRHF(ISYMK)
            LWRK3 = LWORK - KEND3
C
            IF (LWRK3 .LT. 0) THEN
               CALL QUIT('Insufficient work-space in CC_C1T2C-5')
            ENDIF
C
C-------------------------------------
C           Do the actual calculation.
C-------------------------------------
C
            KOFF1 = IGLMVI(ISYMAL,ISYMB) + 1
            KOFF2 = IT1AM(ISYMB,ISYMK) + 1
            KOFF3 = KSCR3
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTB = MAX(NVIR(ISYMB),1)
C
            CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMK),NVIR(ISYMB),
     *                 ONE,XLMPC(KOFF1),NTOTAL,VLTR1(KOFF2),NTOTB,
     *                 ZERO,WORK(KOFF3),NTOTAL)
C
            KOFF1 = KSCR3
            KOFF2 = ILMRHF(ISYMK) + 1
            KOFF3 = IAODIS(ISYMAL,ISYMBE) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTBE = MAX(NBAS(ISYMBE),1)
C
            CALL DGEMM('N','T',NBAS(ISYMAL),NBAS(ISYMBE),NRHF(ISYMK),
     *                 ONE,WORK(KOFF1),NTOTAL,XLAMDH(KOFF2),NTOTBE,
     *                 ONE,C1T2C(KOFF3),NTOTAL)
C
  250    CONTINUE
C
      ENDIF
C
      IF (IOPT .EQ. 3) THEN
C
C--------------------------------------------
C        lamp*lamhc
C--------------------------------------------
C
         ISALBE = ISYRES
C
         DO 260 ISYMAL = 1,NSYM
C
            ISYMBE = MULD2H(ISYMAL,ISALBE)
            ISYMB  = ISYMAL
            ISYMK  = MULD2H(ISYMBE,ISYR)
C
C-----------------------------------
C           Work-space allocation 3.
C-----------------------------------
C
            KSCR3 = KEND1
            KEND3 = KSCR3 + NBAS(ISYMAL)*NRHF(ISYMK)
            LWRK3 = LWORK - KEND3
C
            IF (LWRK3 .LT. 0) THEN
               CALL QUIT('Insufficient work-space in CC_C1T2C-6')
            ENDIF
C
C-------------------------------------
C           Do the actual calculation.
C-------------------------------------
C
            KOFF1 = ILMVIR(ISYMB) + 1
            KOFF2 = IT1AM(ISYMB,ISYMK) + 1
            KOFF3 = KSCR3
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTB = MAX(NVIR(ISYMB),1)
C
            CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMK),NVIR(ISYMB),
     *                 ONE,XLAMDP(KOFF1),NTOTAL,VLTR1(KOFF2),NTOTB,
     *                 ZERO,WORK(KOFF3),NTOTAL)
C
            KOFF1 = KSCR3
            KOFF2 = IGLMRH(ISYMBE,ISYMK) + 1
            KOFF3 = IAODIS(ISYMAL,ISYMBE) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTBE = MAX(NBAS(ISYMBE),1)
C
            CALL DGEMM('N','T',NBAS(ISYMAL),NBAS(ISYMBE),NRHF(ISYMK),
     *                 ONE,WORK(KOFF1),NTOTAL,XLMHC(KOFF2),NTOTBE,
     *                 ONE,C1T2C(KOFF3),NTOTAL)
C
  260    CONTINUE
C
      ENDIF
C
      IF (IOPT. EQ. 2) THEN
C
         CALL DCOPY(NT1AM(ISYMDL),WORK(KLDINT),1,C1T2C,1)
C
      ENDIF
C
      CALL QEXIT('CC_C1T2C')
C
      RETURN
      END
C  /* Deck cc_L1FCK */
      SUBROUTINE CC_L1FCK(RHO2,C1AM,FCKMO,ISYMV,ISYMF,WORK,LWORK)
C
C     Written by Asger Halkier 23/7 - 1995.
C
C     Purpose: To calculate the 12A contribution to rho2.
C
C     Generalised for general symmetries of C1AM and FCKMO,
C     to calculate F-terms
C     rho2(ai,bj) = Paibj(2L1am(ai)FCKMO(jb)-L1am(aj)*FCKMO(ib)
C
C     Ove Christiansen 13-6-1966
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION RHO2(*),C1AM(*),FCKMO(*),WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC_L1FCK')
C
      ISYRES = MULD2H(ISYMV,ISYMF)
C
      IF (LWORK .LT. NT1AM(ISYMF)) THEN
         CALL QUIT('Insufficient work-space-area in CC_L1FCK')
      ENDIF
C
C-----------------------------------
C     Extract the fock matrix F(jb).
C-----------------------------------
C
      DO 50 ISYMC = 1,NSYM
C
         ISYMK = MULD2H(ISYMC,ISYMF)
C
         DO 60 K = 1,NRHF(ISYMK)
C
            DO 70 C = 1,NVIR(ISYMC)
C
               KOFF1 = IFCVIR(ISYMK,ISYMC) + NORB(ISYMK)*(C - 1) + K
               KOFF2 = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
C
               WORK(KOFF2) = FCKMO(KOFF1)
C
   70       CONTINUE
   60    CONTINUE
   50 CONTINUE
C
C--------------------------------------------------------------
C     Calculate coulomb part of contraction of C1AM with FCKMO.
C--------------------------------------------------------------
C
      ISYMAI = ISYMV 
      ISYMBJ = ISYMF 
C
      IF (ISYMAI .EQ. ISYMBJ) THEN
C
         DO 100 NBJ = 1,NT1AM(ISYMBJ)
C
            DO 110 NAI = 1,NBJ
C
               NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
C
               RHO2(NAIBJ) = RHO2(NAIBJ) + TWO*C1AM(NAI)*WORK(NBJ)
     *                                   + TWO*C1AM(NBJ)*WORK(NAI)
C
  110       CONTINUE
C
  100    CONTINUE
C
      ENDIF
C
      IF (ISYMAI .LT. ISYMBJ) THEN
C
         DO 120 NBJ = 1,NT1AM(ISYMBJ)
C
            DO 130 NAI = 1,NT1AM(ISYMAI)
C
               NAIBJ = IT2AM(ISYMAI,ISYMBJ) + NT1AM(ISYMAI)*(NBJ - 1)
     *               + NAI
C
               RHO2(NAIBJ) = RHO2(NAIBJ) + TWO*C1AM(NAI)*WORK(NBJ)
C
  130       CONTINUE
C
  120    CONTINUE
C
      ENDIF
C
      IF (ISYMAI .GT. ISYMBJ) THEN
C
         DO 140 NBJ = 1,NT1AM(ISYMBJ)
C
            DO 150 NAI = 1,NT1AM(ISYMAI)
C
               NAIBJ = IT2AM(ISYMBJ,ISYMAI) + NT1AM(ISYMBJ)*(NAI - 1)
     *               + NBJ
C
               RHO2(NAIBJ) = RHO2(NAIBJ) + TWO*C1AM(NAI)*WORK(NBJ)
C
  150       CONTINUE
C
  140    CONTINUE
C
      ENDIF
C
C---------------------------------------------------------------
C     Calculate exchange part of contraction of C1AM with FCKMO.
C---------------------------------------------------------------
C
      DO 160 ISYMI = 1,NSYM
C
         ISYMB = MULD2H(ISYMF,ISYMI)
C
         DO 170 ISYMA = 1,NSYM
C
            ISYMAI = MULD2H(ISYMA,ISYMI)
            ISYMJ  = MULD2H(ISYMA,ISYMV)
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
            DO 180 J = 1,NRHF(ISYMJ)
C
               DO 190 B = 1,NVIR(ISYMB)
C
                  NBJ = IT1AM(ISYMB,ISYMJ)
     *                + NVIR(ISYMB)*(J - 1) + B
C
                  DO 200 I = 1,NRHF(ISYMI)
C
                     NBI = IT1AM(ISYMB,ISYMI)
     *                   + NVIR(ISYMB)*(I - 1) + B
C
                        IF (ISYMAI .EQ. ISYMBJ) THEN
C
                           DO 210 A = 1,NVIR(ISYMA)
C
                              NAI = IT1AM(ISYMA,ISYMI)
     *                            + NVIR(ISYMA)*(I - 1) + A
C
                              IF (NAI .GT. NBJ) GOTO 210
C
                              NAJ = IT1AM(ISYMA,ISYMJ)
     *                            + NVIR(ISYMA)*(J - 1) + A
C
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + INDEX(NAI,NBJ)
C
                              RHO2(NAIBJ) = RHO2(NAIBJ) -
     *                                    C1AM(NAJ)*WORK(NBI)
                              RHO2(NAIBJ) = RHO2(NAIBJ) -
     *                                    C1AM(NBI)*WORK(NAJ)
  210                      CONTINUE
C
                        ELSE IF (ISYMAI .LT. ISYMBJ) THEN
C
                           DO 220 A = 1,NVIR(ISYMA)
C
                              NAJ = IT1AM(ISYMA,ISYMJ)
     *                            + NVIR(ISYMA)*(J - 1) + A
                              NAI = IT1AM(ISYMA,ISYMI)
     *                            + NVIR(ISYMA)*(I - 1) + A
C
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + NT1AM(ISYMAI)*(NBJ - 1) + NAI
C
                              RHO2(NAIBJ) = RHO2(NAIBJ) -
     *                                    C1AM(NAJ)*WORK(NBI)
C
  220                      CONTINUE
C
                        ELSE IF (ISYMAI .GT. ISYMBJ) THEN
C
                           DO 230 A = 1,NVIR(ISYMA)
C
                              NAJ = IT1AM(ISYMA,ISYMJ)
     *                            + NVIR(ISYMA)*(J - 1) + A
                              NAI = IT1AM(ISYMA,ISYMI)
     *                            + NVIR(ISYMA)*(I - 1) + A
C
                              NAIBJ = IT2AM(ISYMBJ,ISYMAI)
     *                              + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
C
                              RHO2(NAIBJ) = RHO2(NAIBJ) -
     *                                    C1AM(NAJ)*WORK(NBI)
C
  230                      CONTINUE
C
                        ENDIF
C
  200             CONTINUE
  190          CONTINUE
  180       CONTINUE
  170    CONTINUE
  160 CONTINUE
C
      CALL QEXIT('CC_L1FCK')
C
      RETURN
      END
C  /* Deck cc_12c */
      SUBROUTINE CC_21B12C(RHO1,RHO2,CTR1,ISYMV,XLAMDH,ISYMH,
     *                     XMAT,ISYMX,ISINT,WORK,LWORK,
     *                     LUO3,O3FIL,IOPT)
C
C     Written by Asger Halkier & Henrik Koch 13/9 - 1995.
C
C     Purpose: To calculate the 12C contribution to rho2 and the
C              21B contribution to rho1.
C
C     Ove Christiansen 13-6-1996:
C
C     Generalised to calculate contributions to F matrix * vector.
C     ISYMV is symmetry of CTR1.
C     ISINT is symmetry of intermediate on disk.
C     ISYMX is symmetry of intermediate on disk and XMAT.
C     ISYMH is symmetry of lambda matrix.
C
C     rho2(ai,bj) = P(ai,bj)(-sum(k)CTR1(a,k)*L(jbik)
C                   where k may be a barred index.
C
C     IOPT controls if both 21B and 12 C are calculated:
C     if IOPT = 2 both are calculated if IOPT = 1 only
C     21B is calculated.
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION RHO1(*),RHO2(*),CTR1(*),XLAMDH(*),WORK(LWORK),XMAT(*)
      CHARACTER O3FIL*(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
#include "cclr.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC_21B12C')
C
      ISYMLI = MULD2H(ISINT,ISYMH)
      ISYRES = MULD2H(ISYMV,ISYMLI)
      IF (IOPT .EQ. 1) ISYRES = MULD2H(ISYMX,ISYMLI)
      IF ((IOPT .EQ. 2) .AND. (ISYMV .NE. ISYMX)) THEN
         WRITE(LUPRI,*) ' Symmetry mismatch in CC_21B12C '
         CALL QUIT(' Symmetry mismatch in CC_21B12C ')
      ENDIF
C
      DO 100 ISYMD = 1,NSYM
C
         IF (NBAS(ISYMD) .EQ. 0) GOTO 100
C
         ISYIKJ = MULD2H(ISINT,ISYMD)
         ISYMB  = MULD2H(ISYMH,ISYMD)
C
C---------------------------------
C        Allocation of work space.
C---------------------------------
C
         KSCR1 = 1
         KSCR2 = KSCR1 + NMAIJK(ISYIKJ)*NVIR(ISYMB)
         KSCR3 = KSCR2 + NMAIJK(ISYIKJ)
         KEND1 = KSCR3 + NMAIJK(ISYIKJ)*NBAS(ISYMD)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient work space area in CC_21B12C')
         ENDIF
C
C------------------------------------
C        Read integrals from disc.
C------------------------------------
C
         NTOT = NMAIJK(ISYIKJ)*NBAS(ISYMD)
         IOFF = I3ODEL(ISYIKJ,ISYMD) + 1
C
         CALL GETWA2(LUO3,O3FIL,WORK(KSCR3),IOFF,NTOT)
C
C---------------------------
C        Transform AO index.
C---------------------------
C
         KOFF1  = ILMVIR(ISYMB) + 1
C
         NTOIKJ = MAX(NMAIJK(ISYIKJ),1)
         NTOTD  = MAX(NBAS(ISYMD),1)
C
         CALL DGEMM('N','N',NMAIJK(ISYIKJ),NVIR(ISYMB),NBAS(ISYMD),
     *              ONE,WORK(KSCR3),NTOIKJ,XLAMDH(KOFF1),NTOTD,
     *              ZERO,WORK(KSCR1),NTOIKJ)
C
C-------------------------------------------------------------
C        Calculate 2 coulomb - exchange for (ik|jb) integrals.
C-------------------------------------------------------------
C
         DO 110 B = 1,NVIR(ISYMB)
C
            DO 120 ISYMJ = 1,NSYM
C
               ISYMIK = MULD2H(ISYIKJ,ISYMJ)
C
               DO 130 J = 1,NRHF(ISYMJ)
C
                  DO 140 ISYMK = 1,NSYM
C
                     ISYMI  = MULD2H(ISYMIK,ISYMK)
                     ISYMJK = MULD2H(ISYMJ,ISYMK)
C
                     DO 150 K = 1,NRHF(ISYMK)
C
                        DO 160 I = 1,NRHF(ISYMI)
C
                           NIKJ  = IMAIJK(ISYMIK,ISYMJ)
     *                           + NMATIJ(ISYMIK)*(J - 1)
     *                           + IMATIJ(ISYMI,ISYMK)
     *                           + NRHF(ISYMI)*(K - 1) + I
                           NJKI  = IMAIJK(ISYMJK,ISYMI)
     *                           + NMATIJ(ISYMJK)*(I - 1)
     *                           + IMATIJ(ISYMJ,ISYMK)
     *                           + NRHF(ISYMJ)*(K - 1) + J
                           NBIKJ = NMAIJK(ISYIKJ)*(B - 1) + NIKJ
                           NBJKI = NMAIJK(ISYIKJ)*(B - 1) + NJKI
C
                           WORK(KSCR2 + NIKJ - 1) =
     *                              TWO*WORK(KSCR1 + NBIKJ - 1) -
     *                              WORK(KSCR1 + NBJKI - 1)
C
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
C
            DO 170 ISYMJ = 1,NSYM
C
               ISYMBJ = MULD2H(ISYMJ,ISYMB)
               ISYMAI = MULD2H(ISYMBJ,ISYRES)
               ISYMIK = MULD2H(ISYMJ,ISYIKJ)
C
               LAVL   = LWORK - KSCR3 - NT1AM(ISYMAI)
C
               IF (LAVL .LT. 0) THEN
                  CALL QUIT('Insufficient space in CC_21B12C')
               ENDIF
C
               DO 180 J = 1,NRHF(ISYMJ)
C
C------------------------------------------------
C                 Calculate the 21B contribution.
C------------------------------------------------
C
                  NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
                  IF (.NOT.(CC2)) THEN
C
                     IF (ISYMBJ .EQ. ISYRES) THEN
C
                        KOFF1 = KSCR2 + IMAIJK(ISYMIK,ISYMJ)
     *                                + NMATIJ(ISYMIK)*(J - 1)
C
                        RHO1(NBJ) = RHO1(NBJ) - DDOT(NMATIJ(ISYMIK),
     *                              XMAT,1,WORK(KOFF1),1)
C
                     ENDIF
C
                  ENDIF
C
C-----------------------------------------------------------------
C                 IF only 21B contributions then goto end of loop.
C-----------------------------------------------------------------
C
                  IF ( IOPT .EQ. 1 ) GOTO 180
C
C-----------------------------------------------------------
C                 Contract integrals with trial vector CTR1.
C-----------------------------------------------------------
C
                  DO 190 ISYMK = 1,NSYM
C
                     ISYMI  = MULD2H(ISYMK,ISYMIK)
                     ISYMA  = MULD2H(ISYMK,ISYMV)
                     ISYMAJ = MULD2H(ISYMA,ISYMJ)
                     ISYMBI = MULD2H(ISYMB,ISYMI)
C
                     KOFF1 = IT1AM(ISYMA,ISYMK) + 1
                     KOFF2 = KSCR2 + IMAIJK(ISYMIK,ISYMJ)
     *                             + NMATIJ(ISYMIK)*(J - 1)
     *                             + IMATIJ(ISYMI,ISYMK)
                     KOFF3 = KSCR3 + IT1AM(ISYMA,ISYMI)
C
                     NTOTA = MAX(NVIR(ISYMA),1)
                     NTOTI = MAX(NRHF(ISYMI),1)
C
                     CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),
     *                          NRHF(ISYMK),-ONE,CTR1(KOFF1),NTOTA,
     *                          WORK(KOFF2),NTOTI,ZERO,WORK(KOFF3),
     *                          NTOTA)
C
  190             CONTINUE
C
C------------------------------------------------
C                 Calculate the 12C contribution.
C------------------------------------------------
C
                  IF (ISYMAI .EQ. ISYMBJ)
     *                      WORK(KSCR3+NBJ-1) = TWO*WORK(KSCR3+NBJ-1)
C
                  DO 200 NAI = 1,NT1AM(ISYMAI)
C
                     IF (ISYMAI .EQ. ISYMBJ) THEN
                        NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                     ELSE IF (ISYMAI .LT. ISYMBJ) THEN
                        NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                        + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                     ELSE IF (ISYMAI .GT. ISYMBJ) THEN
                        NAIBJ = IT2AM(ISYMBJ,ISYMAI)
     *                        + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
                     ENDIF
C
                     RHO2(NAIBJ) = RHO2(NAIBJ) + WORK(KSCR3+NAI-1)
C
  200             CONTINUE
  180          CONTINUE
  170       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CC_21B12C')
C
      RETURN
      END
      SUBROUTINE CC_INT3O(X3OINT,DSRHF,XLAMDH,ISYMH,XLAMDP,
     *                    WORK,LWORK,IDEL,ISYMD,LUO3,O3FIL)
C
C     Written by Henrik Koch & Asger Halkier 27/7 - 1995
C
C     Purpose: To calculate and write to disc an integral batch with
C              three occupied indices for a given delta.
C
C     Modified by Asger Halkier to return integrals in X3OINT as
C     well as writing them to disc 28/10 - 1995.
C
C     Generalized to also calculated integrals (j del i k) where k is barred
C     thus lamda is not total symmetric.
C
C     Ove Christiansen 13-6-1996
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION X3OINT(*),DSRHF(*),WORK(LWORK)
      DIMENSION XLAMDH(*),XLAMDP(*)
      CHARACTER O3FIL*(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "cclr.h"
      LOGICAL LWTD
C
      CALL QENTER('CC_INT3O')
C
C-------------------------------------------
C     Work space allocation 1 * outer loops.
C-------------------------------------------
C
      ISYABJ = MULD2H(ISYMOP,ISYMD)
      ISYIKJ = MULD2H(ISYABJ,ISYMH)
C
      DO 100 ISYMJ = 1,NSYM
C
         ISALBE = MULD2H(ISYABJ,ISYMJ)
C
         DO 110 J = 1,NRHF(ISYMJ)
C
C------------------------------------------------------------
C           Work space allocation 1 * unpacking of integrals.
C------------------------------------------------------------
C
            KUNINT = 1
            KEND1  = KUNINT + N2BST(ISALBE)
            LWRK1  = LWORK  - KEND1
C
            IF (LWRK1 .LT. 0) THEN
               CALL QUIT('1-Insufficient work space area in CC_INT3O')
            ENDIF
C
            KOFF1 = IDSRHF(ISALBE,ISYMJ) + NNBST(ISALBE)*(J - 1) + 1
C
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KUNINT))
C
            DO 120 ISYMK = 1,NSYM
C
C-----------------------------------------------------------------------
C              Transform remaining AO-indices of integrals to occ. space
C-----------------------------------------------------------------------
C
               ISYMBE = MULD2H(ISYMK,ISYMH)
               ISYMAL = MULD2H(ISALBE,ISYMBE)
               ISYMI  = ISYMAL
               ISYMIK = MULD2H(ISYMI,ISYMK)
C
               KINMA1 = KEND1
               KEND2  = KINMA1 + NBAS(ISYMAL)*NRHF(ISYMK)
               LWRK2  = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  CALL QUIT('2-Insufficient work space area in CCINT3O')
               ENDIF
C
               KOFF1 = KUNINT + IAODIS(ISYMAL,ISYMBE)
               KOFF2 = IGLMRH(ISYMBE,ISYMK) + 1
               KOFF3 = KINMA1
C
               NTOTA = MAX(NBAS(ISYMAL),1)
               NTOTB = MAX(NBAS(ISYMBE),1)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMK),NBAS(ISYMBE),
     *                    ONE,WORK(KOFF1),NTOTA,XLAMDH(KOFF2),NTOTB,
     *                    ZERO,WORK(KOFF3),NTOTA)
C
               KOFF1 = ILMRHF(ISYMI) + 1
               KOFF2 = KINMA1
               KOFF3 = IMAIJK(ISYMIK,ISYMJ) + NMATIJ(ISYMIK)*(J - 1)
     *               + IMATIJ(ISYMI,ISYMK) + 1
C
               NTOTA = MAX(NBAS(ISYMAL),1)
               NTOTI = MAX(NRHF(ISYMI),1)
C
               CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMK),NBAS(ISYMAL),
     *                    ONE,XLAMDP(KOFF1),NTOTA,WORK(KOFF2),NTOTA,
     *                    ZERO,X3OINT(KOFF3),NTOTI)
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
      LWTD = .TRUE.
      IF (LUO3 .LT. 0) LWTD = .FALSE.
C
      IF ((LWTD) .AND. (.NOT. MP2)) THEN
C
C--------------------------------
C        Write integrals on disc.
C--------------------------------
C
         D    = IDEL - IBAS(ISYMD)
C
         NTOT = NMAIJK(ISYIKJ)
         IOFF = I3ODEL(ISYIKJ,ISYMD) + NTOT*(D - 1) + 1
C
         IF (NTOT .GT. 0) THEN
            CALL PUTWA2(LUO3,O3FIL,X3OINT,IOFF,NTOT)
         ENDIF
      ENDIF
C
      CALL QEXIT('CC_INT3O')
C
      RETURN
      END
C  /* Deck ccll_lamtra */
      SUBROUTINE CCLL_LAMTRA(CTR1,ISYMV,XLAMDP,ISYMLP,CT1AO)
C
C     Written by Asger Halkier 24/7 - 1995
C
C     Purpose: To calculate the contraction of the trialvector CTR1
C              with the lampda-matrix XLAMDP to give zeta(Alfa,J)
C
C     Slight generalisation for input of symmetries
C     ISYMV = symmetry of singles vector
C     ISYMLP = symmetry of general Lambdamatrix.
C     Ove Christiansen 14-6-1996
C
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION CTR1(*),XLAMDP(*),CT1AO(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
C
      CALL QENTER('CCLL_LAMTRA')
C
C---------------------------------------------
C     Half-transformation to AO-basis of CTR1.
C---------------------------------------------
C
      ISYMAO = MULD2H(ISYMV,ISYMLP)
      ISYMCJ = ISYMV 
C
      CALL DZERO(CT1AO,NGLMDT(ISYMAO))
C
      DO 100 ISYMAL = 1,NSYM
C
         ISYMC = MULD2H(ISYMAL,ISYMLP)
         ISYMJ = MULD2H(ISYMC,ISYMCJ)
C
         KOFF1 = IGLMVI(ISYMAL,ISYMC) + 1
         KOFF2 = IT1AM(ISYMC,ISYMJ) + 1
         KOFF3 = IGLMRH(ISYMAL,ISYMJ) + 1
C
         NTOBAL = MAX(NBAS(ISYMAL),1)
         NTOVIC = MAX(NVIR(ISYMC),1)
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMJ),NVIR(ISYMC),ONE,
     *              XLAMDP(KOFF1),NTOBAL,CTR1(KOFF2),NTOVIC,ZERO,
     *              CT1AO(KOFF3),NTOBAL)
C
  100 CONTINUE
C
      CALL QEXIT('CCLL_LAMTRA')
C
      RETURN
      END
C  /* Deck cc_21c */
      SUBROUTINE CC_21DC(RHO1,CTR2,ISYMV,SCRM,ISYMM,SCRCM,ISYMCM,
     *                   XINT,XLMH1,ISYMH1,
     *                   XLMP1,ISYMP1,XLMH2,ISYMH2,XLMP2,ISYMP2,
     *                   WORK,LWORK,IDEL,ISYMD,IOPT,ICON)
C
C     Written by Asger Halkier 10/10 - 1995.
C
C     Purpose: To calculate the 21C contribution to rho1
C
C     Generalised to calculate the 21D contribution to rho1 as well
C     by Asger Halkier 21/12 - 1995, in order to speed up CC2-
C     calculations!
C
C     IMPORTANT NOTE: This routine assumes AO symmetric
C     integrals and hence cannot be used directly with
C     complex basis functions!!!!!!!
C
C     Ove Christiansen 14-6-1996:
C         Significant modifications to included 
C         non-tot sym CTR2 - ISYMV
C         and introduction of new sections 
C         (IOPT = 1) section with 
C         XLMHC and XLMPC with sym ISYMLC
C
C     Ordinary Contributions:
C
C        if (icon .eq. 2) rho1(a,i) = -sum(dkl)[t(dk,al)*(dk|il)]
C
C        rho1(b,j) = -sum(del)[t(dl,ej)*(eb|dl)]
C
C
C     IOPT .EQ. 2 
C
C        if (icon .eq. 2) rho1 = -sum(dkl)[t(dk,al)*((d-bar k|il)+(dk-bar|il)]
C        if (icon .eq. 2) rho1 = -sum(dkl)[t(dk,al)*(dk|il-bar)]
C
C        rho1(b,j) = -sum(del)[t(dl,ej)*(e-bar b|dl)]
C    
C        rho1(b,j) = -sum(del)[t(dl,ej)*(eb|dl-bar)]
C                             (calculate by making (alf, beta|l-bar delta))
C
#include "implicit.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "cclr.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION RHO1(*),CTR2(*),SCRM(*),SCRCM(*),XINT(*), 
     *          XLMH1(*),XLMP1(*),XLMH2(*),XLMP2(*),
     *          WORK(LWORK)
C
      CALL QENTER('CC_21DC')
C
      ISYML1 = MULD2H(ISYMH1,ISYMP1)
      ISYML2 = MULD2H(ISYMH2,ISYMP2)
      IF (ISYML1 .NE. ISYML2 ) THEN
         CALL QUIT('Something is wrong in cc_21dc')
      ENDIF
      ISYRES = MULD2H(ISYMV,ISYML1)
      IF (ISYRES .NE. MULD2H(ISYMCM,ISYMD)) THEN
         CALL QUIT('Symmetry mismatch in cc_21dc')
      ENDIF
C
      ISYINT = MULD2H(ISYMD,ISYMOP)
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      KDSRHF = 1
      KEND1  = KDSRHF + NDSRHF(ISYMD)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient work space in CC_21DC-1')
      ENDIF
C
C-------------------------------------------------------------
C     Transform one integral distribution index to occ. space.
C-------------------------------------------------------------
C
      CALL CCTRBT(XINT,WORK(KDSRHF),XLMH1,ISYMH1,WORK(KEND1),LWRK1,
     *            ISYINT)
C
      ISYMI  = MULD2H(ISYMD,ISYMP2)
      ISYMA  = MULD2H(ISYMI,ISYRES)
C
      DO 100 ISYML = 1,NSYM
C
         ISALBE = MULD2H(ISYML,ISYINT)
         ISYMDK = MULD2H(ISALBE,ISYML1)
         ISYMAL = MULD2H(ISYML,ISYMA)
         ISYMEB = MULD2H(ISALBE,ISYML1)
C
C----------------------------------
C        Work space allocation two.
C----------------------------------
C
         IF (ICON .EQ. 2) THEN
            KMOINT = KEND1
            KUPINT = KMOINT + NT1AM(ISYMDK)*NRHF(ISYML)
         ELSE
            KUPINT = KEND1
         ENDIF
         KEND2  = KUPINT + N2BST(ISALBE)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_21DC-2')
         ENDIF
C
         IF ( ICON .EQ. 2) CALL DZERO(WORK(KMOINT),
     *                  NT1AM(ISYMDK)*NRHF(ISYML))
C
C----------------------------------------------------
C        Unpack integral distribution (al be | L de).
C----------------------------------------------------
C
         DO 110 L = 1,NRHF(ISYML)
C
            KOFF1 = KDSRHF + IDSRHF(ISALBE,ISYML)
     *            + NNBST(ISALBE)*(L - 1)
C
            CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KUPINT))
C
C-------------------------------------------------------------------------
C           Transform remaining unpacked distribution indices to MO-basis.
C-------------------------------------------------------------------------
C
            DO 120 ISYALF = 1,NSYM
C
               ISYMBE = MULD2H(ISYALF,ISALBE)
C
               ISYME  = MULD2H(ISYALF,ISYMP1)
               ISYMB  = MULD2H(ISYMBE,ISYMH1)
               ISYMJ  = MULD2H(ISYMB,ISYRES)
C
               ISYEJL = ISYMM
               ISYMEJ = MULD2H(ISYML,ISYEJL)
C
C------------------------------------------
C              Work space allocation three.
C------------------------------------------
C
               KSCRD = KEND2
               KSCRI = KSCRD + NBAS(ISYALF)*NVIR(ISYMB)
               KENDD = KSCRI + NVIR(ISYMB)*NVIR(ISYME)
               LWRKD = LWORK - KENDD
C
               IF (LWRKD .LT. 0) THEN
                  CALL QUIT('Insufficient work space in CC_21DC-3')
               ENDIF
C
               KOFF10 = KUPINT + IAODIS(ISYALF,ISYMBE)
               KOFF11 = IGLMVI(ISYMBE,ISYMB) + 1
C
               NTOTAL = MAX(NBAS(ISYALF),1)
               NTOTBE = MAX(NBAS(ISYMBE),1)
C
               CALL DGEMM('N','N',NBAS(ISYALF),NVIR(ISYMB),
     *                    NBAS(ISYMBE),ONE,WORK(KOFF10),NTOTAL,
     *                    XLMH1(KOFF11),NTOTBE,ZERO,WORK(KSCRD),
     *                    NTOTAL)
C
               KOFF12 = IGLMVI(ISYALF,ISYME) + 1
               NTOTAL = MAX(NBAS(ISYALF),1)
               NTOTB  = MAX(NVIR(ISYMB),1)
C
               CALL DGEMM('T','N',NVIR(ISYMB),NVIR(ISYME),
     *                    NBAS(ISYALF),ONE,WORK(KSCRD),NTOTAL,
     *                    XLMP1(KOFF12),NTOTAL,
     *                    ZERO,WORK(KSCRI),NTOTB)
C
               KOFF13 = IT2BCD(ISYMEJ,ISYML) 
     *                + NT1AM(ISYMEJ)*(L - 1)
     *                + IT1AM(ISYME,ISYMJ) + 1
               KOFF14 = IT1AM(ISYMB,ISYMJ) + 1
C
               NTOTB  = MAX(NVIR(ISYMB),1)
               NTOTE  = MAX(NVIR(ISYME),1)
C
C------------------------------------------------------------------
C              Contraction with SCRM and add D-term to RHO1-vector.
C------------------------------------------------------------------
C
               CALL DGEMM('N','N',NVIR(ISYMB),NRHF(ISYMJ),NVIR(ISYME),
     *                    ONE,WORK(KSCRI),NTOTB,SCRM(KOFF13),NTOTE,
     *                    ONE,RHO1(KOFF14),NTOTB)
C
C----------------------------------------------------------
C              Extra section: T2AO with C1. times integral.
C----------------------------------------------------------
C
               IF (IOPT . EQ. 2 ) THEN
C
                 ISYME  = MULD2H(ISYALF,ISYMP2)
                 ISYMB  = MULD2H(ISYMBE,ISYMH1)
                 ISYMJ  = MULD2H(ISYMB,ISYRES)
C
                 ISYEJL = ISYMCM
                 ISYMEJ = MULD2H(ISYML,ISYEJL)
C
C----------------------------------------
C                Work space allocation 4.
C----------------------------------------
C
                 KSCRD = KEND2
                 KSCRI = KSCRD + NBAS(ISYALF)*NVIR(ISYMB)
                 KENDD = KSCRI + NVIR(ISYMB)*NVIR(ISYME)
                 LWRKD = LWORK - KENDD
C
                 IF (LWRKD .LT. 0) THEN
                  CALL QUIT('Insufficient work space in CC_21DC-4')
                 ENDIF
C
                 KOFF10 = KUPINT + IAODIS(ISYALF,ISYMBE)
                 KOFF11 = IGLMVI(ISYMBE,ISYMB) + 1
C
                 NTOTAL = MAX(NBAS(ISYALF),1)
                 NTOTBE = MAX(NBAS(ISYMBE),1)
C
                 CALL DGEMM('N','N',NBAS(ISYALF),NVIR(ISYMB),
     *                    NBAS(ISYMBE),ONE,WORK(KOFF10),NTOTAL,
     *                    XLMH1(KOFF11),NTOTBE,ZERO,WORK(KSCRD),
     *                    NTOTAL)
C
                 KOFF12 = IGLMVI(ISYALF,ISYME) + 1
                 NTOTAL = MAX(NBAS(ISYALF),1)
                 NTOTB  = MAX(NVIR(ISYMB),1)
C
                 CALL DGEMM('T','N',NVIR(ISYMB),NVIR(ISYME),
     *                      NBAS(ISYALF),ONE,WORK(KSCRD),NTOTAL,
     *                      XLMP2(KOFF12),NTOTAL,
     *                      ZERO,WORK(KSCRI),NTOTB)
C
                 KOFF13 = IT2BCD(ISYMEJ,ISYML) 
     *                  + NT1AM(ISYMEJ)*(L - 1)
     *                  + IT1AM(ISYME,ISYMJ) + 1
                 KOFF14 = IT1AM(ISYMB,ISYMJ) + 1
C
                 NTOTB  = MAX(NVIR(ISYMB),1)
                 NTOTE  = MAX(NVIR(ISYME),1)
C
C---------------------------------------------------------------------
C                Contraction with SCRCM and add D-term to RHO1-vector.
C---------------------------------------------------------------------
C
                 CALL DGEMM('N','N',NVIR(ISYMB),NRHF(ISYMJ),
     *                      NVIR(ISYME),ONE,WORK(KSCRI),NTOTB,
     *                      SCRCM(KOFF13),NTOTE,
     *                      ONE,RHO1(KOFF14),NTOTB)
C 
               ENDIF
C
               IF (ICON .EQ. 2) THEN
C
                  ISYVID = MULD2H(ISYALF,ISYMP1)
                  ISYMK  = MULD2H(ISYMBE,ISYMH1)
C
C-----------------------------------------
C                 Work space allocation 5.
C-----------------------------------------
C
                  KSCR1 = KEND2
                  KEND3 = KSCR1 + NBAS(ISYALF)*NRHF(ISYMK)
                  LWRK3 = LWORK - KEND3
C
                  IF (LWRK3 .LT. 0) THEN
                     CALL QUIT('Insufficient work space in CC_21DC-5')
                  ENDIF
C
                  KOFF2  = KUPINT + IAODIS(ISYALF,ISYMBE)
                  KOFF3  = IGLMRH(ISYMBE,ISYMK) + 1
C
                  NTOTAL = MAX(NBAS(ISYALF),1)
                  NTOTBE = MAX(NBAS(ISYMBE),1)
C
                  CALL DGEMM('N','N',NBAS(ISYALF),NRHF(ISYMK),
     *                       NBAS(ISYMBE),ONE,WORK(KOFF2),NTOTAL,
     *                       XLMH1(KOFF3),NTOTBE,ZERO,WORK(KSCR1),
     *                       NTOTAL)
C
                  KOFF4  = IGLMVI(ISYALF,ISYVID) + 1
                  KOFF5  = KMOINT + NT1AM(ISYMDK)*(L - 1)
     *                            + IT1AM(ISYVID,ISYMK)
C
                  NTOTAL = MAX(NBAS(ISYALF),1)
                  NTOTVD = MAX(NVIR(ISYVID),1)
C
                  CALL DGEMM('T','N',NVIR(ISYVID),NRHF(ISYMK),
     *                       NBAS(ISYALF),ONE,XLMP1(KOFF4),NTOTAL,
     *                       WORK(KSCR1),NTOTAL,ONE ,WORK(KOFF5),
     *                       NTOTVD)
C
C---------------------------------------
C                 Extra section (dk-bar| 
C---------------------------------------
C
                  IF (IOPT .EQ. 2) THEN
C
                    ISYVID = MULD2H(ISYALF,ISYMP2)
                    ISYMK  = MULD2H(ISYMBE,ISYMH2)
C
                    KSCR1 = KEND2
                    KEND3 = KSCR1 + NBAS(ISYALF)*NRHF(ISYMK)
                    LWRK3 = LWORK - KEND3
C
                    IF (LWRK3 .LT. 0) THEN
                       CALL QUIT('Insufficient work space in CC_21DC-6')
                    ENDIF
C
                    KOFF2  = KUPINT + IAODIS(ISYALF,ISYMBE)
                    KOFF3  = IGLMRH(ISYMBE,ISYMK) + 1
C
                    NTOTAL = MAX(NBAS(ISYALF),1)
                    NTOTBE = MAX(NBAS(ISYMBE),1)
C
                    CALL DGEMM('N','N',NBAS(ISYALF),NRHF(ISYMK),
     *                       NBAS(ISYMBE),ONE,WORK(KOFF2),NTOTAL,
     *                       XLMH2(KOFF3),NTOTBE,ZERO,WORK(KSCR1),
     *                       NTOTAL)
C
                    KOFF4  = IGLMVI(ISYALF,ISYVID) + 1
                    KOFF5  = KMOINT + NT1AM(ISYMDK)*(L - 1)
     *                            + IT1AM(ISYVID,ISYMK)
C
                    NTOTAL = MAX(NBAS(ISYALF),1)
                    NTOTVD = MAX(NVIR(ISYVID),1)
C
                    CALL DGEMM('T','N',NVIR(ISYVID),NRHF(ISYMK),
     *                       NBAS(ISYALF),ONE,XLMP2(KOFF4),NTOTAL,
     *                       WORK(KSCR1),NTOTAL,ONE ,WORK(KOFF5),
     *                       NTOTVD)
C
                  ENDIF
C
               ENDIF
C
  120       CONTINUE
C
  110    CONTINUE
C
         IF (ICON .EQ. 2) THEN
C
C---------------------------------------------
C           Contraction with trialvector CTR2.
C---------------------------------------------
C
            LNEED = LWRK2 - NVIR(ISYMA)
C
            IF (LNEED .LT. ZERO) THEN
               CALL QUIT('Insufficient work space in CC_21DC-7')
            ENDIF
C
            CALL DZERO(WORK(KEND2),NVIR(ISYMA))
C
            DO 130 NDK = 1,NT1AM(ISYMDK)
C
               KOFF6 = IT2SQ(ISYMAL,ISYMDK) + NT1AM(ISYMAL)*(NDK - 1)
     *               + IT1AM(ISYMA,ISYML) + 1
               KOFF7 = KMOINT + NDK - 1
C
               NTOTA = MAX(NVIR(ISYMA),1)
C
               CALL DGEMV('N',NVIR(ISYMA),NRHF(ISYML),ONE,CTR2(KOFF6),
     *                    NTOTA,WORK(KOFF7),NT1AM(ISYMDK),ONE,
     *                    WORK(KEND2),1)
C
  130       CONTINUE
C
C-----------------------------------------------------------------
C           Final scaling of C-term with XLAMDP * storage in RHO1.
C-----------------------------------------------------------------
C
            DO 140 I = 1,NRHF(ISYMI)
C
               KOFF8 = IGLMRH(ISYMD,ISYMI) + NBAS(ISYMD)*(I - 1)
     *               + IDEL - IBAS(ISYMD)
               KOFF9 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
C
               CALL DAXPY(NVIR(ISYMA),-XLMP2(KOFF8),WORK(KEND2),1,
     *                    RHO1(KOFF9),1)
C
  140       CONTINUE
C
         ENDIF
C
  100 CONTINUE
C
C---------------------------------------------------------------------------
C     Section for calc of 
C        rho1(b,j) = -sum(del)[t(dl,ej)*(eb|dl-bar)]
C        rho1(a,i) = -sum(del)[t(dk,al)*(dk|il-bar)]
C                             (calculated by making (alf, beta|l-bar delta))
C---------------------------------------------------------------------------
C   
      IF ( IOPT .EQ. 2 ) THEN
C
         ISYINT = MULD2H(ISYMD,ISYMOP)
         ISTINT = MULD2H(ISYMD,ISYMH2)
C
C--------------------------------
C        Work space allocation 8.
C--------------------------------
C
         KDSRHF = 1
         KEND1  = KDSRHF + NDSRHF(ISTINT) 
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_21DC-8 ')
         ENDIF
C
C----------------------------------------------------------------
C        Transform one integral distribution index to occ. space.
C----------------------------------------------------------------
C
         CALL CCTRBT(XINT,WORK(KDSRHF),XLMH2,ISYMH2,WORK(KEND1),
     *               LWRK1,ISYINT)
C
         DO 102 ISYML = 1,NSYM
C
            ISALBE = MULD2H(ISYML,ISTINT)
            ISYMEB = ISALBE
            ISYMDK = MULD2H(ISALBE,ISYMOP)
C
C-----------------------------------
C           Work space allocation 9.
C-----------------------------------
C
            IF (ICON .EQ. 2) THEN
               KMOINT = KEND1
               KUPINT = KMOINT + NT1AM(ISYMDK)*NRHF(ISYML)
            ELSE
               KUPINT = KEND1
            ENDIF
            KEND2  = KUPINT + N2BST(ISALBE)
            LWRK2  = LWORK  - KEND2
C
            IF (LWRK2 .LT. 0) THEN
               CALL QUIT('Insufficient work space in CC_21DC-9')
            ENDIF
C
            IF ( ICON .EQ. 2) CALL DZERO(WORK(KMOINT),
     *                  NT1AM(ISYMDK)*NRHF(ISYML))
C
C-----------------------------------------------------------
C           Unpack integral distribution (al be | L-bar de).
C-----------------------------------------------------------
C
            DO 112 L = 1,NRHF(ISYML)
C
               KOFF1 = KDSRHF + IDSRHF(ISALBE,ISYML)
     *               + NNBST(ISALBE)*(L - 1)
C
               CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KUPINT))
C
C----------------------------------------------------------------------------
C              Transform remaining unpacked distribution indices to MO-basis.
C----------------------------------------------------------------------------
C
               DO 122 ISYALF = 1,NSYM
C
                  ISYMBE = MULD2H(ISYALF,ISALBE)
C
                  ISYME  = MULD2H(ISYALF,ISYMP2)
                  ISYMB  = MULD2H(ISYMEB,ISYME)
                  ISYMJ  = MULD2H(ISYMB,ISYRES)
C
                  ISYEJL = ISYMM
                  ISYMEJ = MULD2H(ISYML,ISYEJL)
C
C------------------------------------------
C                 Work space allocation 10.
C------------------------------------------
C
                  KSCRD = KEND2
                  KSCRI = KSCRD + NBAS(ISYALF)*NVIR(ISYMB)
                  KENDD = KSCRI + NVIR(ISYMB)*NVIR(ISYME)
                  LWRKD = LWORK - KENDD
C
                  IF (LWRKD .LT. 0) THEN
                     CALL QUIT('Insufficient work space in CC_21DC-10')
                  ENDIF
C
                  KOFF10 = KUPINT + IAODIS(ISYALF,ISYMBE)
                  KOFF11 = IGLMVI(ISYMBE,ISYMB) + 1
C
                  NTOTAL = MAX(NBAS(ISYALF),1)
                  NTOTBE = MAX(NBAS(ISYMBE),1)
C
                  CALL DGEMM('N','N',NBAS(ISYALF),NVIR(ISYMB),
     *                    NBAS(ISYMBE),ONE,WORK(KOFF10),NTOTAL,
     *                    XLMH1(KOFF11),NTOTBE,ZERO,WORK(KSCRD),
     *                    NTOTAL)
C
                  KOFF12 = IGLMVI(ISYALF,ISYME) + 1
                  NTOTAL = MAX(NBAS(ISYALF),1)
                  NTOTB  = MAX(NVIR(ISYMB),1)
C
                  CALL DGEMM('T','N',NVIR(ISYMB),NVIR(ISYME),
     *                    NBAS(ISYALF),ONE,WORK(KSCRD),
     *                    NTOTAL,XLMP2(KOFF12),NTOTAL,
     *                    ZERO,WORK(KSCRI),NTOTB)
C
                  KOFF13 = IT2BCD(ISYMEJ,ISYML)
     *                   + NT1AM(ISYMEJ)*(L - 1)
     *                   + IT1AM(ISYME,ISYMJ) + 1
                  KOFF14 = IT1AM(ISYMB,ISYMJ) + 1
C
                  NTOTB  = MAX(NVIR(ISYMB),1)
                  NTOTE  = MAX(NVIR(ISYME),1)
C
C---------------------------------------------------------------------
C                 Contraction with SCRM and add D-term to RHO1-vector.
C---------------------------------------------------------------------
C
                  CALL DGEMM('N','N',NVIR(ISYMB),NRHF(ISYMJ),
     *                       NVIR(ISYME),ONE,WORK(KSCRI),NTOTB,
     *                       SCRM(KOFF13),NTOTE,
     *                       ONE,RHO1(KOFF14),NTOTB)
                  IF (ICON .EQ. 2) THEN
C
C-------------------------------------------------------------
C                    Transform remaining unpacked distribution 
C                    indices to MO-basis.
C                    Now with usual Lambda matrices: LH1, LP2
C-------------------------------------------------------------
C
                     ISYVID = MULD2H(ISYALF,ISYMP2)
                     ISYMBE = MULD2H(ISYALF,ISALBE)
C
                     ISYMK  = MULD2H(ISYMBE,ISYMH1)
C
C---------------------------------------------
C                    Work space allocation 11.
C---------------------------------------------
C
                     KSCR1 = KEND2
                     KEND3 = KSCR1 + NBAS(ISYALF)*NRHF(ISYMK)
                     LWRK3 = LWORK - KEND3
C
                     IF (LWRK3 .LT. 0) THEN
                        CALL QUIT('Insufficient work '//
     &                            'space in CC_21DC-11')
                     ENDIF
C
                     KOFF2  = KUPINT + IAODIS(ISYALF,ISYMBE)
                     KOFF3  = IGLMRH(ISYMBE,ISYMK) + 1
C
                     NTOTAL = MAX(NBAS(ISYALF),1)
                     NTOTBE = MAX(NBAS(ISYMBE),1)
C
                     CALL DGEMM('N','N',NBAS(ISYALF),NRHF(ISYMK),
     *                       NBAS(ISYMBE),ONE,WORK(KOFF2),NTOTAL,
     *                       XLMH1(KOFF3),NTOTBE,ZERO,WORK(KSCR1),
     *                       NTOTAL)
C
                     KOFF4  = IGLMVI(ISYALF,ISYVID) + 1
                     KOFF5  = KMOINT + NT1AM(ISYMDK)*(L - 1)
     *                         + IT1AM(ISYVID,ISYMK)
C
                     NTOTAL = MAX(NBAS(ISYALF),1)
                     NTOTVD = MAX(NVIR(ISYVID),1)
C
                     CALL DGEMM('T','N',NVIR(ISYVID),NRHF(ISYMK),
     *                       NBAS(ISYALF),ONE,XLMP2(KOFF4),NTOTAL,
     *                       WORK(KSCR1),NTOTAL,ONE ,WORK(KOFF5),
     *                       NTOTVD)
C
                  ENDIF
C
  122          CONTINUE
C
  112       CONTINUE
C
            IF (ICON .EQ. 2) THEN
C
C------------------------------------------------
C              Contraction with trialvector CTR2.
C------------------------------------------------
C
               ISYMAL = MULD2H(ISYML,ISYMA)
               ISYMDK = ISALBE
C
               LNEED = LWRK2 - NVIR(ISYMA)
C
               IF (LNEED .LT. ZERO) THEN
                  CALL QUIT('Insufficient work space in CC_21DC-12')
               ENDIF
C
               CALL DZERO(WORK(KEND2),NVIR(ISYMA))
C
               DO 132 NDK = 1,NT1AM(ISYMDK)
C
                  KOFF6 = IT2SQ(ISYMAL,ISYMDK) + NT1AM(ISYMAL)*(NDK-1)
     *               + IT1AM(ISYMA,ISYML) + 1
                  KOFF7 = KMOINT + NDK - 1
C
                  NTOTA = MAX(NVIR(ISYMA),1)
C
               CALL DGEMV('N',NVIR(ISYMA),NRHF(ISYML),ONE,CTR2(KOFF6),
     &                    NTOTA,WORK(KOFF7),NT1AM(ISYMDK),ONE,
     &                    WORK(KEND2),1)
C
  132          CONTINUE
C
C-----------------------------------------------------------------
C              Final scaling of C-term with XLAMDP * storage in RHO1.
C-----------------------------------------------------------------
C
               DO 142 I = 1,NRHF(ISYMI)
C
                  KOFF8 = IGLMRH(ISYMD,ISYMI) + NBAS(ISYMD)*(I - 1)
     *               + IDEL - IBAS(ISYMD)
                  KOFF9 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
C
                  CALL DAXPY(NVIR(ISYMA),-XLMP2(KOFF8),WORK(KEND2),1,
     *                    RHO1(KOFF9),1)
C
  142          CONTINUE
C
             ENDIF
C
  102    CONTINUE
C
      ENDIF
C
      CALL QEXIT('CC_21DC')
C
      RETURN
      END
C  /* Deck cc_11a */
      SUBROUTINE CC_11A(RHO1,FOCKG,ISYMFG,XLAMDH,XLAMDP,WORK,LWORK)
C
C     Written by Henrik Koch & Asger Halkier 27/6 - 1995
C
C     Version: 1.0
C
C     Purpose: To calculate the Fock-matrix-alike term (11A)
C              entering the 1.1-block of the left hand side
C              transformation of a trial vector!
C
C     Ove Christiansen 6-10-1996
C              Input symmtery of FOCK matrix: ISYMFG
C
C
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION RHO1(*),FOCKG(*)
      DIMENSION XLAMDH(*),XLAMDP(*),WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
C
      CALL QENTER('CC_11A')
C
      ISYRHO = MULD2H(ISYMFG,ISYMOP)
C
      IF (LWORK .LT. NT1AO(ISYRHO)) THEN
         CALL QUIT('Insufficient space in CC_11A')
      ENDIF
C
C-----------------------------------
C     Transformation of gamma-index.
C-----------------------------------
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMG = ISYMI
         ISYMD = MULD2H(ISYMI,ISYRHO)
C
         NTOTG = MAX(NBAS(ISYMG),1)
         NTOTD = MAX(NBAS(ISYMD),1)
C
         KOFF1 = IAODIS(ISYMG,ISYMD) + 1
         KOFF2 = ILMRHF(ISYMI) + 1
         KOFF3 = IT1AO(ISYMD,ISYMI) + 1
C
         CALL DGEMM('T','N',NBAS(ISYMD),NRHF(ISYMI),NBAS(ISYMG),
     *              ONE,FOCKG(KOFF1),NTOTG,XLAMDP(KOFF2),NTOTG,
     *              ZERO,WORK(KOFF3),NTOTD)
C
  100 CONTINUE
C
C-----------------------------------
C     Transformation of delta-index.
C-----------------------------------
C
      DO 110 ISYMI = 1,NSYM
C
         ISYMA = MULD2H(ISYMI,ISYRHO)
         ISYMD = ISYMA
C
         NTOTD = MAX(NBAS(ISYMD),1)
         NTOTA = MAX(NVIR(ISYMA),1)
C
         KOFF1 = ILMVIR(ISYMA) + 1
         KOFF2 = IT1AO(ISYMD,ISYMI) + 1
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMD),
     *              ONE,XLAMDH(KOFF1),NTOTD,WORK(KOFF2),NTOTD,
     *              ONE,RHO1(KOFF3),NTOTA)
C
  110 CONTINUE
C
      CALL QEXIT('CC_11A')
C
      RETURN
      END
C  /* Deck cc_12b */
      SUBROUTINE CC_12B(RHO2,DSRHF,C1AO,ISC1AO,
     *                  XLAMDH,WORK,LWORK,IDEL,ISYMD)
C
C     Written by Asger Halkier 22/11 - 1995
C
C     Version: 1.0
C
C     Purpose: To calculate the 12B contribution to RHO2 in a
C     CC2-calculation, where it's not possible to use the BF-routine!
C
C     Ove Christiansen 7-10-1996:
C
C     Very small generalization to input symmetry of C1AO.
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION RHO2(*), DSRHF(*), C1AO(*), XLAMDH(*), WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC_12B')
C
      ISYINT = MULD2H(ISYMD,ISYMOP)
      ISYBCO = ISYMD
      ISYAEX = ISYMD
      ISYRES = MULD2H(ISYMOP,ISC1AO)
C
      DO 100 ISYMJ = 1,NSYM
C
         ISALBE = MULD2H(ISYMJ,ISYINT)
         ISBJCO = MULD2H(ISYMJ,ISYBCO)
         ISAICO = MULD2H(ISBJCO,ISYRES)
         ISAJEX = MULD2H(ISYMJ,ISYAEX)
         ISBIEX = MULD2H(ISAJEX,ISYRES)
C
C----------------------------------
C        Work space allocation one.
C----------------------------------
C
         KAOINT = 1
         KSCR1C = KAOINT + N2BST(ISALBE)
         KSCR1E = KSCR1C + NT1AM(ISAICO)
         KEND1  = KSCR1E + NT1AM(ISBIEX)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_12B')
         ENDIF
C
         DO 110 J = 1,NRHF(ISYMJ)
C
C-------------------------------------------
C           Square up integral distribution.
C-------------------------------------------
C
            KOFF1 = IDSRHF(ISALBE,ISYMJ) + NNBST(ISALBE)*(J - 1) + 1
C
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KAOINT))
C
C-----------------------------------------------------------------
C           Contraction of integrals with lampda matrices (1 gen).
C-----------------------------------------------------------------
C
            DO 120 ISYMAL = 1,NSYM
C
               ISYMI  = MULD2H(ISYMAL,ISC1AO)
               ISYMBE = MULD2H(ISYMAL,ISALBE)
               ISYACO = ISYMBE
               ISYBEX = ISYMBE
C
C----------------------------------------
C              Work space allocation two.
C----------------------------------------
C
               KSCR2 = KEND1
               KEND2 = KSCR2 + NBAS(ISYMBE)*NRHF(ISYMI)
               LWRK2 = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  CALL QUIT('Insufficient work space in CC_12B')
               ENDIF
C
               KOFF2  = KAOINT + IAODIS(ISYMAL,ISYMBE)
               KOFF3  = IGLMRH(ISYMAL,ISYMI) + 1
C
               NTOTAL = MAX(NBAS(ISYMAL),1)
               NTOTBE = MAX(NBAS(ISYMBE),1)
C
               CALL DGEMM('T','N',NBAS(ISYMBE),NRHF(ISYMI),
     &                    NBAS(ISYMAL),ONE,WORK(KOFF2),NTOTAL,
     &                    C1AO(KOFF3),NTOTAL,ZERO,WORK(KSCR2),NTOTBE)
C
               KOFF4  = ILMVIR(ISYACO) + 1
               KOFF5  = KSCR1C + IT1AM(ISYACO,ISYMI)
C
               NTOTBE = MAX(NBAS(ISYMBE),1)
               NTOTAC = MAX(NVIR(ISYACO),1)
C
               CALL DGEMM('T','N',NVIR(ISYACO),NRHF(ISYMI),
     &                    NBAS(ISYMBE),ONE,XLAMDH(KOFF4),NTOTBE,
     &                    WORK(KSCR2),NTOTBE,ZERO,WORK(KOFF5),NTOTAC)
C
               KOFF6  = ILMVIR(ISYBEX) + 1
               KOFF7  = KSCR1E + IT1AM(ISYBEX,ISYMI)
C
               NTOTBE = MAX(NBAS(ISYMBE),1)
               NTOTBX = MAX(NVIR(ISYBEX),1)
C
               CALL DGEMM('T','N',NVIR(ISYBEX),NRHF(ISYMI),
     &                    NBAS(ISYMBE),ONE,XLAMDH(KOFF6),NTOTBE,
     &                    WORK(KSCR2),NTOTBE,ZERO,WORK(KOFF7),NTOTBX)
C
  120       CONTINUE
C
C---------------------------------------------------------------------
C           Scaling with XLAMDH-element and storage of coulomb-result.
C---------------------------------------------------------------------
C
            DO 130 B = 1,NVIR(ISYBCO)
C
               NBJ = IT1AM(ISYBCO,ISYMJ) + NVIR(ISYBCO)*(J - 1) + B
               NDB = ILMVIR(ISYBCO) + NBAS(ISYMD)*(B - 1)
     &             + IDEL - IBAS(ISYMD)
C
               FACT = TWO*XLAMDH(NDB)
C
               IF (ISAICO .EQ. ISBJCO) THEN
C
                  DO 140 NAI = 1,NT1AM(ISAICO)
C
                     NAIBJ = IT2AM(ISAICO,ISBJCO) + INDEX(NAI,NBJ)
C
                     RHO2(NAIBJ) = RHO2(NAIBJ)
     &                           + FACT*WORK(KSCR1C + NAI - 1)
C
                     IF (NAI .EQ. NBJ) THEN
C
                        RHO2(NAIBJ) = RHO2(NAIBJ)
     &                              + FACT*WORK(KSCR1C + NAI - 1)
C
                     ENDIF
C
  140             CONTINUE
C
               ELSE IF (ISAICO .LT. ISBJCO) THEN
C
                  NAIBJ = IT2AM(ISAICO,ISBJCO)
     &                  + NT1AM(ISAICO)*(NBJ - 1) + 1
C
                  CALL DAXPY(NT1AM(ISAICO),FACT,WORK(KSCR1C),1,
     &                       RHO2(NAIBJ),1)
C
               ELSE IF (ISAICO .GT. ISBJCO) THEN
C
                  NAIBJ = IT2AM(ISBJCO,ISAICO) + NBJ
C
                  CALL DAXPY(NT1AM(ISAICO),FACT,WORK(KSCR1C),1,
     &                       RHO2(NAIBJ),NT1AM(ISBJCO))
C
               ENDIF
C
  130       CONTINUE
C
C----------------------------------------------------------------------
C           Scaling with XLAMDH-element and storage of exchange-result.
C----------------------------------------------------------------------
C
            DO 150 A = 1,NVIR(ISYAEX)
C
               NDA = ILMVIR(ISYAEX) + NBAS(ISYMD)*(A - 1)
     &             + IDEL - IBAS(ISYMD)
C
               FACT = -XLAMDH(NDA)
C
               DO 160 ISYBEX = 1,NSYM
C
                  ISYMI  = MULD2H(ISYBEX,ISBIEX)
                  ISYMAI = MULD2H(ISYAEX,ISYMI)
                  ISYMBJ = MULD2H(ISYMAI,ISYRES)
C
                  IF (ISYMAI .EQ. ISYMBJ) THEN
C
                     DO 170 I = 1,NRHF(ISYMI)
C
                        NAI = IT1AM(ISYAEX,ISYMI)
     &                      + NVIR(ISYAEX)*(I - 1) + A
C
                        DO 180 B = 1,NVIR(ISYBEX)
C
                           NBI   = IT1AM(ISYBEX,ISYMI)
     &                           + NVIR(ISYBEX)*(I - 1) + B
                           NBJ   = IT1AM(ISYBEX,ISYMJ)
     &                           + NVIR(ISYBEX)*(J - 1) + B
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     &                           + INDEX(NAI,NBJ)
C
                           RHO2(NAIBJ) = RHO2(NAIBJ)
     &                                 + FACT*WORK(KSCR1E + NBI - 1)
C
                           IF (NAI .EQ. NBJ) THEN
C
                              RHO2(NAIBJ) = RHO2(NAIBJ)
     &                                    + FACT*WORK(KSCR1E + NBI - 1)
C
                           ENDIF
C
  180                   CONTINUE
  170                CONTINUE
C
                  ELSE IF (ISYMAI .LT. ISYMBJ) THEN
C
                     DO 190 I = 1,NRHF(ISYMI)
C
                        NAI = IT1AM(ISYAEX,ISYMI)
     &                      + NVIR(ISYAEX)*(I - 1) + A
C
                        DO 200 B = 1,NVIR(ISYBEX)
C
                           NBI   = IT1AM(ISYBEX,ISYMI)
     &                           + NVIR(ISYBEX)*(I - 1) + B
                           NBJ   = IT1AM(ISYBEX,ISYMJ)
     &                           + NVIR(ISYBEX)*(J - 1) + B
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     &                           + NT1AM(ISYMAI)*(NBJ - 1) + NAI
C
                           RHO2(NAIBJ) = RHO2(NAIBJ)
     &                                 + FACT*WORK(KSCR1E + NBI - 1)
C
  200                   CONTINUE
  190                CONTINUE
C
                  ELSE IF (ISYMAI .GT. ISYMBJ) THEN
C
                     DO 210 I = 1,NRHF(ISYMI)
C
                        NAI   = IT1AM(ISYAEX,ISYMI)
     &                        + NVIR(ISYAEX)*(I - 1) + A
                        NAIBJ = IT2AM(ISYMBJ,ISYMAI)
     &                        + NT1AM(ISYMBJ)*(NAI - 1)
     &                        + IT1AM(ISYBEX,ISYMJ)
     &                        + NVIR(ISYBEX)*(J - 1) + 1
                        NBI   = KSCR1E + IT1AM(ISYBEX,ISYMI)
     &                        + NVIR(ISYBEX)*(I - 1)
C
                        CALL DAXPY(NVIR(ISYBEX),FACT,WORK(NBI),1,
     &                             RHO2(NAIBJ),1)
C
  210                CONTINUE
C
                  ENDIF
C
  160          CONTINUE
  150       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CC_12B')
C
      RETURN
      END
C  /* Deck cclt_p21i */
      SUBROUTINE CCLT_P21I(X21INT,ISYINT,WORK,LWORK,
     &                     IT2BCD,NT2BCD,IT1AM,NT1AM,MVIR)
C
C     Written by Asger Halkier 26/10 - 1995
C
C     Version: 1.0
C
C     Purpose: To transpose the occupied indices of the intermediates
C              entering some of the terms in the 2.1-block.
C
C
C#include "implicit.h"
      IMPLICIT NONE
#include "priunit.h"
#include "ccorb.h"
C#include "ccsdsym.h"

      INTEGER LWORK,ISYINT
      INTEGER IT2BCD(8,8),NT2BCD(8),IT1AM(8,8),NT1AM(8),MVIR(8)
      INTEGER ISYMJ, J, ISYMI, I, ISYMA, ISYMAI, ISYMAJ, NAIJ, NAJI
#if defined (SYS_CRAY)
      REAL X21INT(*), WORK(LWORK)
#else
      DOUBLE PRECISION X21INT(*), WORK(LWORK)
#endif
C
      CALL QENTER('CCLT_P21I')
C
      IF (LWORK .LT. NT2BCD(ISYINT)) THEN
         CALL QUIT('Insufficient work space in CCLT_I21P')
      ENDIF
C
      CALL DCOPY(NT2BCD(ISYINT),X21INT,1,WORK,1)
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMAI = MULD2H(ISYMJ,ISYINT)
C
         DO 110 ISYMI = 1,NSYM
C
            ISYMA  = MULD2H(ISYMI,ISYMAI)
            ISYMAJ = MULD2H(ISYMA,ISYMJ)
C
C----------------------------------------------------
C           Do the transposition of occupied indices.
C----------------------------------------------------
C
            DO 120 J = 1,NRHF(ISYMJ)
C
               DO 130 I = 1,NRHF(ISYMI)
C
                  NAIJ = IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J - 1)
     &                 + IT1AM(ISYMA,ISYMI) + MVIR(ISYMA)*(I - 1) + 1
                  NAJI = IT2BCD(ISYMAJ,ISYMI) + NT1AM(ISYMAJ)*(I - 1)
     &                 + IT1AM(ISYMA,ISYMJ) + MVIR(ISYMA)*(J - 1) + 1
C
                  CALL DCOPY(MVIR(ISYMA),WORK(NAIJ),1,X21INT(NAJI),1)
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CCLT_P21I')
C
      RETURN
      END
C  /* Deck cclt_ct1ao */
      SUBROUTINE CCLT_CT1AO(CTR1,XLAMDP,CT1AO)
C
C     Written by Asger Halkier 24/7 - 1995
C
C     Version 1.0
C
C     Purpose: To calculate the contraction of the trialvector CTR1
C              with the lampda-matrix XLAMDP to give zeta(Alfa,J)!
C
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION CTR1(*),XLAMDP(*),CT1AO(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
C
      CALL QENTER('CCLT_CT1AO')
C
C---------------------------------------------
C     Half-transformation to AO-basis of CTR1.
C---------------------------------------------
C
      ISYMCJ = ISYMTR
C
      DO 100 ISYMAL = 1,NSYM
C
         ISYMC = ISYMAL
         ISYMJ = MULD2H(ISYMC,ISYMCJ)
C
         KOFF1 = ILMVIR(ISYMC) + 1
         KOFF2 = IT1AM(ISYMC,ISYMJ) + 1
         KOFF3 = IGLMRH(ISYMAL,ISYMJ) + 1
C
         NTOBAL = MAX(NBAS(ISYMAL),1)
         NTOVIC = MAX(NVIR(ISYMC),1)
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMJ),NVIR(ISYMC),ONE,
     &              XLAMDP(KOFF1),NTOBAL,CTR1(KOFF2),NTOVIC,ZERO,
     &              CT1AO(KOFF3),NTOBAL)
C
  100 CONTINUE
C
      CALL QEXIT('CCLT_CT1AO')
C
      RETURN
      END
C  /* Deck cclt_s21i */
      SUBROUTINE CCLT_S21I(XINTS,XINT1,XINT2,ISYINT,FACT)
C
C     Written by Asger Halkier 30/10 - 1995
C
C     Version: 1.0
C
C     Purpose: To resort and add together 2 21-block-
C              intermediates!
C
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XINTS(*), XINT1(*), XINT2(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
C
      CALL QENTER('CCLT_S21I')
C
C-----------------------------
C     Initialize result array.
C-----------------------------
C
      CALL DZERO(XINTS,NT2BCD(ISYINT))
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMAL = MULD2H(ISYMJ,ISYINT)
C
         DO 110 ISYML = 1,NSYM
C
            ISYMA  = MULD2H(ISYML,ISYMAL)
            ISYMJL = MULD2H(ISYMJ,ISYML)
C
C------------------------------
C           Do the calculation.
C------------------------------
C
            DO 120 J = 1,NRHF(ISYMJ)
C
               DO 130 L = 1,NRHF(ISYML)
C
                  NJL  = IMATIJ(ISYMJ,ISYML) + NRHF(ISYMJ)*(L - 1) + J
                  NALJ = IT2BCD(ISYMAL,ISYMJ) + NT1AM(ISYMAL)*(J - 1)
     &                 + IT1AM(ISYMA,ISYML) + NVIR(ISYMA)*(L - 1) + 1
                  NAJL = IT2AIJ(ISYMA,ISYMJL)
     &                 + NVIR(ISYMA)*(NJL - 1) + 1
C
                  CALL DCOPY(NVIR(ISYMA),XINT1(NALJ),1,XINTS(NAJL),1)
C
                  CALL DAXPY(NVIR(ISYMA),FACT,XINT2(NALJ),1,
     &                       XINTS(NAJL),1)
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CCLT_S21I')
C
      RETURN
      END
C  /* Deck cclt_pi3o */
      SUBROUTINE CCLT_PI3O(X3OTR,X3OINT,ISYINT)
C
C     Written by Asger Halkier 30/10 - 1995
C
C     Version: 1.0
C
C     Purpose: To resort integral batch (il|jdel) to
C              (jl|idel)!
C
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION X3OTR(*), X3OINT(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
C
      CALL QENTER('CCLT_PI3O')
C
C-----------------------------
C     Initialize result array.
C-----------------------------
C
      CALL DZERO(X3OTR,NMAIJK(ISYINT))
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMIL = MULD2H(ISYMJ,ISYINT)
C
         DO 110 ISYMI = 1,NSYM
C
            ISYML  = MULD2H(ISYMI,ISYMIL)
            ISYMJL = MULD2H(ISYMJ,ISYML)
C
C----------------------------
C           Do the resorting.
C----------------------------
C
            DO 120 J = 1,NRHF(ISYMJ)
C
               DO 130 L = 1,NRHF(ISYML)
C
                  NILJ = IMAIJK(ISYMIL,ISYMJ) + NMATIJ(ISYMIL)*(J - 1)
     &                 + IMATIJ(ISYMI,ISYML) + NRHF(ISYMI)*(L - 1) + 1
                  NJL  = IMATIJ(ISYMJ,ISYML) + NRHF(ISYMJ)*(L - 1) + J
                  NJLI = IMAIJK(ISYMJL,ISYMI) + NJL
C
                  CALL DCOPY(NRHF(ISYMI),X3OINT(NILJ),1,
     &                       X3OTR(NJLI),NMATIJ(ISYMJL))
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CCLT_PI3O')
C
      RETURN
      END
C  /* Deck ccpinaje */
      SUBROUTINE CCPINAJE(RES,SUB,ISYMA,ISYMJ,ISYME)
C
C     Written by Asger Halkier 31/10 - 1995
C
C     Version: 1.0
C
C     Purpose: To resort integrals (aj,e del) to (a,ej del)!
C
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION RES(*), SUB(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
C
      CALL QENTER('CCPINAJE')
C
      DO 100 E = 1,NVIR(ISYME)
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            NEJ  = IT1AM(ISYME,ISYMJ) + NVIR(ISYME)*(J - 1) + E
            NAEJ = NVIR(ISYMA)*(NEJ - 1) + 1
            NAJE = NVIR(ISYMA)*NRHF(ISYMJ)*(E - 1)
     &           + NVIR(ISYMA)*(J - 1) + 1
C
            CALL DCOPY(NVIR(ISYMA),SUB(NAJE),1,RES(NAEJ),1)
C
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CCPINAJE')
C
      RETURN
      END
C  /* Deck ccrinaej */
      SUBROUTINE CCRINAEJ(RES,SUB,ISYMA,ISYME,ISYMJ)
C
C     Written by Asger Halkier 31/10 - 1995
C
C     Version: 1.0
C
C     Purpose: To resort integrals (ae,j del) to (a,ej del)!
C
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION RES(*), SUB(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
C
      CALL QENTER('CCRINAEJ')
C
      DO 100 J = 1,NRHF(ISYMJ)
C
         DO 110 E = 1,NVIR(ISYME)
C
            NEJ  = IT1AM(ISYME,ISYMJ) + NVIR(ISYME)*(J - 1) + E
            NAEJ = NVIR(ISYMA)*(NEJ - 1) + 1
            NAJE = NVIR(ISYMA)*NVIR(ISYME)*(J - 1)
     &           + NVIR(ISYMA)*(E - 1) + 1
C
            CALL DCOPY(NVIR(ISYMA),SUB(NAJE),1,RES(NAEJ),1)
C
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CCRINAEJ')
C
      RETURN
      END
C  /* Deck cclt_11bc */
      SUBROUTINE CCLT_11BC(RHO1,CTR1,FCKMO,WORK,LWORK)
C
C     Written by Asger Halkier & Henrik Koch 19/6-1995
C
C     Purpose: To calculate the 11B and 11C contributions to rho1.
C
#include "implicit.h"
#include "dummy.h"
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RHO1(*),CTR1(*),FCKMO(*)
      DIMENSION WORK(LWORK)
C
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
C
      CALL QENTER('CCLT_11BC')
C
C---------------------------
C     Work space allocation.
C---------------------------
C
      KRMAT = 1
      KSMAT = KRMAT + NMATAB(ISYMOP)
      KEND1 = KSMAT + NMATIJ(ISYMOP)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCLT_11BC')
      ENDIF
C
C-----------------------------------------------
C     Read RMAT and SMAT written by energy code.
C-----------------------------------------------
C
      LUE1 = -1
      CALL GPOPEN(LUE1,'CC_E1IM','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUE1)
      READ(LUE1) (WORK(KRMAT+I-1),I = 1,NMATAB(ISYMOP))
      CALL GPCLOSE(LUE1,'KEEP')
C
      LUE2 = -1
      CALL GPOPEN(LUE2,'CC_E2IM','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUE2)
      READ(LUE2) (WORK(KSMAT+I-1),I = 1,NMATIJ(ISYMOP))
      CALL GPCLOSE(LUE2,'KEEP')
C
C-----------------------------------------------------------------------
C     Contraction with trial vector RHO(ai) = c1(bi)*r(ba)-c1(aj)*s(ij).
C-----------------------------------------------------------------------
C
      DO 160 ISYMJ = 1,NSYM
C
         ISYMA = MULD2H(ISYMTR,ISYMJ)
         ISYMI = MULD2H(ISYMJ,ISYMOP)
C
         KOFF1 = IT1AM(ISYMA,ISYMJ) + 1
         KOFF2 = KSMAT + IMATIJ(ISYMI,ISYMJ)
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTI = MAX(NRHF(ISYMI),1)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &              -ONE,CTR1(KOFF1),NTOTA,WORK(KOFF2),NTOTI,ONE,
     &              RHO1(KOFF3),NTOTA)
C
  160 CONTINUE
C
      DO 170 ISYMA = 1,NSYM
C
         ISYMB = MULD2H(ISYMA,ISYMOP)
         ISYMI = MULD2H(ISYMTR,ISYMB)
C
         KOFF1 = KRMAT + IMATAB(ISYMB,ISYMA)
         KOFF2 = IT1AM(ISYMB,ISYMI) + 1
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTB = MAX(NVIR(ISYMB),1)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &              ONE,WORK(KOFF1),NTOTB,CTR1(KOFF2),NTOTB,ONE,
     &              RHO1(KOFF3),NTOTA)
C
  170 CONTINUE
C
      CALL QEXIT('CCLT_11BC')
C
      RETURN
      END
C  /* Deck cc_eitr */
      SUBROUTINE CC_EITR(E1INT,E2INT,WORK,LWORK,ISYINT)
C
C     Written by Asger Halkier 21/11 - 1995
C
C     Version: 1.0
C
C     Purpose: To transpose the indices of the two E-intermediate-
C              matrices!
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION E1INT(*), E2INT(*), WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
C
      CALL QENTER('CC_EITR')
C
      IF ((LWORK .LT. NMATAB(ISYINT)) .OR. (LWORK .LT. NMATIJ(ISYINT)))
     &   CALL QUIT('Insufficient work space in CC_EITR')
C
C------------------------------------
C     Transpose virtual intermediate.
C------------------------------------
C
      CALL DCOPY(NMATAB(ISYINT),E1INT,1,WORK,1)
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMB = MULD2H(ISYMA,ISYINT)
C
         DO 110 B = 1,NVIR(ISYMB)
C
            NORIG = IMATAB(ISYMA,ISYMB) + NVIR(ISYMA)*(B - 1) + 1
            NTRAN = IMATAB(ISYMB,ISYMA) + B
C
            CALL DCOPY(NVIR(ISYMA),WORK(NORIG),1,E1INT(NTRAN),
     &                 NVIR(ISYMB))
C
  110    CONTINUE
  100 CONTINUE
C
C-------------------------------------
C     Transpose occupied intermediate.
C-------------------------------------
C
      CALL DCOPY(NMATIJ(ISYINT),E2INT,1,WORK,1)
C
      DO 120 ISYMI = 1,NSYM
C
         ISYMJ = MULD2H(ISYMI,ISYINT)
C
         DO 130 J = 1,NRHF(ISYMJ)
C
            NORIG = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + 1
            NTRAN = IMATIJ(ISYMJ,ISYMI) + J
C
            CALL DCOPY(NRHF(ISYMI),WORK(NORIG),1,E2INT(NTRAN),
     &                 NRHF(ISYMJ))
C
  130    CONTINUE
  120 CONTINUE
C
      CALL QEXIT('CC_EITR')
C
      RETURN
      END
C  /* Deck cc2_cc2ffc */
      SUBROUTINE CC_RCC2FF(OMEGA2,T2AM,XLAMDP,XLAMDH,XLPC1,XLHC1,
     *                     WORK,LWORK,ISYMTR) 
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     CC2 finite diff. contribution 
C
C     Ove Christiansen - march 1997 
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      USE PELIB_INTERFACE, ONLY: USE_PELIB
C
#include "implicit.h"
C
      DIMENSION OMEGA2(*),T2AM(*),WORK(LWORK)
      DIMENSION XLAMDP(*),XLAMDH(*),XLPC1(*),XLHC1(*)
C
      PARAMETER (ONE=1.0D00)
C
      LOGICAL ETRAN,FCKCON
C
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
#include "ccsections.h"
#include "qm3.h"
!#include "qmmm.h"
#include "ccfield.h"
#include "leinf.h"
      REAL*8, ALLOCATABLE :: FOCKMAT(:), FOCKTEMP(:)
C
C      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC_RCC2FF')
C
      IF (((NFIELD.GT.0).OR.CCSLV.OR.USE_PELIB()).AND.NONHF) THEN
C
         KFOCK1 = 1
         KFOCK2 = KFOCK1 + MAX(N2BST(ISYMTR),N2BST(ISYMOP))
         KEMAT1 = KFOCK2 + MAX(N2BST(ISYMTR),N2BST(ISYMOP))
         KEMAT2 = KEMAT1 + NEMAT1(ISYMTR)
         KEND1  = KEMAT2 + NMATIJ(ISYMTR)
         LWRK1  = LWORK  - KEND1
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Needed:', KEND1, 'Available:', LWORK
            CALL QUIT('Insufficient memory for allocation in cc_rcc2ff')
         ENDIF
C
         CALL DZERO(WORK(KFOCK2),MAX(N2BST(ISYMTR),N2BST(ISYMOP)))
         CALL DZERO(WORK(KEMAT1),NEMAT1(ISYMTR))
         CALL DZERO(WORK(KEMAT2),NMATIJ(ISYMTR))
C
         DO 13 IF = 1, NFIELD
            FF =  EFIELD(IF)
            CALL CC_ONEP(WORK(KFOCK2),WORK(KEND1),LWRK1,FF,1,LFIELD(IF))
 13      CONTINUE
C
C-------------------------------------
C        Solvent contribution.
C        Put into one-electron integrals.
C       SLV02,OC
C-------------------------------------
C
         IF (CCSLV .AND. (.NOT. CCMM )) THEN
            DTIME = SECOND()
            CALL CCSL_RHSTG(WORK(KFOCK2),WORK(KEND1),LWRK1)
            IF (IPRINT .GT. 5) THEN
              WRITE(LUPRI,*)'Time used in CCSL_RHSTG (CC_RCC2FF) 1',
     *                       SECOND() - DTIME
            END IF
         ENDIF
C
C-------------------------------------
C        Solvent contribution.
C        Put into one-electron integrals.
C       CCMM02,JA+AO
C-------------------------------------
C
         IF (CCMM) THEN
            DTIME = SECOND()
            IF (.NOT. NYQMMM) THEN
               CALL CCMM_RHSTG(WORK(KFOCK2),WORK(KEND1),LWRK1)
            ELSE IF (NYQMMM) THEN
              IF (HFFLD) THEN
                CALL CCMM_ADDGHF(WORK(KFOCK2),WORK(KEND1),LWRK1)
              ELSE 
                CALL CCMM_ADDG(WORK(KFOCK2),WORK(KEND1),LWRK1)
              END IF
            END IF
            IF (IPRINT .GT. 5) THEN
              WRITE(LUPRI,*)'Time used in CCMM_RHSTG (CC_RCC2FF) 1',
     *                       SECOND() - DTIME
            END IF
         ENDIF
         IF (USE_PELIB()) THEN
             ALLOCATE(FOCKMAT(NNBASX),
     &                FOCKTEMP(MAX(N2BST(ISYMTR),N2BST(ISYMOP))))
             IF (HFFLD) THEN
                 CALL GET_FROM_FILE('FOCKMHF',NNBASX,FOCKMAT)
             ELSE
             CALL GET_FROM_FILE('FOCKMAT',NNBASX,FOCKMAT)
             END IF
             CALL DSPTSI(NBAS,FOCKMAT,FOCKTEMP)
             CALL DAXPY(MAX(N2BST(ISYMTR),N2BST(ISYMOP)),1.0d0,
     &                  FOCKTEMP,1,WORK(KFOCK2),1)
             DEALLOCATE(FOCKMAT,FOCKTEMP)
         END IF
C
C
C--------------------------------------
C
         CALL DCOPY(N2BST(ISYMOP),WORK(KFOCK2),1,WORK(KFOCK1),1)
         ISYFAO = 1
         ISYMPA = ISYMTR
         ISYMHO = 1
         CALL CC_FCKMO(WORK(KFOCK1),XLPC1,XLAMDH,
     *                 WORK(KEND1),LWRK1,ISYFAO,ISYMPA,ISYMHO)
C
         ISYFAO = 1
         ISYMPA = 1
         ISYMHO = ISYMTR
         CALL CC_FCKMO(WORK(KFOCK2),XLAMDP,XLHC1,
     *                 WORK(KEND1),LWRK1,ISYFAO,ISYMPA,ISYMHO)
C
         CALL DAXPY(N2BST(ISYMTR),ONE,WORK(KFOCK1),1,WORK(KFOCK2),1)
C
         ETRAN  = .FALSE.
         FCKCON = .TRUE.
         ISYMV  = ISYMOP
         ISYMEI = ISYMTR
         CALL CCRHS_EFCK(WORK(KEMAT1),WORK(KEMAT2),XLAMDH,
     *                   WORK(KFOCK2),WORK(KEND1),LWRK1,
     *                   FCKCON,ETRAN,ISYMEI)
         CALL CCRHS_E(OMEGA2,T2AM,WORK(KEMAT1),WORK(KEMAT2),
     *                WORK(KEND1),LWRK1,ISYMV,ISYMEI)
C
      ENDIF
C
      CALL QEXIT('CC_RCC2FF')
C
      RETURN
      END
      SUBROUTINE CC_LCC2FF(RHO1,XLAMDP,XLAMDH,XLPC1,XLHC1,ISYMLA, 
     *                     XI,YI,ISYMIM,WORK,LWORK,IOPT) 
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     CC2 finite diff. contribution 
C
C     Ove Christiansen - march 1997 
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      USE PELIB_INTERFACE, ONLY: USE_PELIB
#include "implicit.h"
C
      DIMENSION RHO1(*),WORK(LWORK)
      DIMENSION XLAMDP(*),XLAMDH(*),XI(*),YI(*),XLPC1(*),XLHC1(*)
C
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
#include "ccsections.h"
#include "ccfield.h"
#include "leinf.h"
#include "qm3.h"
      REAL*8, ALLOCATABLE :: FOCKMAT(:), FOCKTEMP(:)
!#include "qmmm.h"
C      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC_LCC2FF')
C
      ONE = 1.0D0
C
      KFOCK1 = 1
      KFOCK2 = KFOCK1 + MAX(N2BST(ISYMOP),N2BST(ISYMLA))
      KEND1  = KFOCK2 + MAX(N2BST(ISYMOP),N2BST(ISYMLA))
      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Needed:', KEND1, 'Available:', LWORK
           CALL QUIT('Insufficient memory for allocation in cc_lcc2ff')
      ENDIF
C
      CALL DZERO(WORK(KFOCK1),MAX(N2BST(ISYMOP),N2BST(ISYMLA)))
C
      DO 13 IF = 1, NFIELD
         FF =  EFIELD(IF)
         CALL CC_ONEP(WORK(KFOCK1),WORK(KEND1),LWRK1,FF,1,LFIELD(IF))
 13   CONTINUE
C
C-------------------------------------
C        Solvent contribution.
C        Put into one-electron integrals.
C       SLV02,OC
C-------------------------------------
C
         IF (CCSLV .AND. (.NOT. CCMM )) THEN
            DTIME = SECOND()
            CALL CCSL_RHSTG(WORK(KFOCK1),WORK(KEND1),LWRK1)
            IF (IPRINT .GT. 5) THEN
              WRITE(LUPRI,*)'Time used in CCSL_RHSTG (CC_RCC2FF) 1',
     *                       SECOND() - DTIME
            END IF
         ENDIF
C
C-------------------------------------
C        Solvent contribution.
C        Put into one-electron integrals.
C       CCMM02,JA+AO
C-------------------------------------
C
         IF (CCMM) THEN
            DTIME = SECOND()
            IF (.NOT. NYQMMM) THEN
               CALL CCMM_RHSTG(WORK(KFOCK1),WORK(KEND1),LWRK1)
            ELSE IF (NYQMMM) THEN
              IF (HFFLD) THEN
                CALL CCMM_ADDGHF(WORK(KFOCK1),WORK(KEND1),LWRK1)
              ELSE
               CALL CCMM_ADDG(WORK(KFOCK1),WORK(KEND1),LWRK1)
              END IF
            END IF
            IF (IPRINT .GT. 5) THEN
              WRITE(LUPRI,*)'Time used in CCMM_RHSTG (CC_RCC2FF) 1',
     *                       SECOND() - DTIME
            END IF
         ENDIF
         IF (USE_PELIB()) THEN
             ALLOCATE(FOCKMAT(NNBASX),
     &                FOCKTEMP(MAX(N2BST(ISYMOP),N2BST(ISYMLA))))
             IF (HFFLD) THEN
                 CALL GET_FROM_FILE('FOCKMHF',NNBASX,FOCKMAT)
             ELSE
             CALL GET_FROM_FILE('FOCKMAT',NNBASX,FOCKMAT)
             END IF
             CALL DSPTSI(NBAS,FOCKMAT,FOCKTEMP)
             CALL DAXPY(MAX(N2BST(ISYMOP),N2BST(ISYMLA)),1.0d0,FOCKTEMP,
     &                  1,WORK(KFOCK1),1)
             DEALLOCATE(FOCKMAT,FOCKTEMP)
         END IF
C
C--------------------------------------
C
      IF (IOPT .EQ.1) THEN
         ISYFAO = 1
         ISYMPA = 1
         ISYMHO = 1
         CALL CC_FCKMO(WORK(KFOCK1),XLAMDP,XLAMDH,
     *                 WORK(KEND1),LWRK1,ISYFAO,ISYMPA,ISYMHO)
      ENDIF 
C
      CALL CC_21EFM(RHO1,WORK(KFOCK1),ISYMLA,XI,YI,
     *              ISYMIM,WORK(KEND1),LWRK1)
C
      CALL QEXIT('CC_LCC2FF')
C
      RETURN
      END
      SUBROUTINE CC_FCC2FF(OMEGA2,C2AM,ISYMTR,XLAMDP,XLAMDH,XLPC1,
     *                     XLHC1,ISYMLA,WORK,LWORK,ISIDE) 
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     CC2 finite diff. contribution to F-matrix transformation
C
C     Ove Christiansen - march 1997 
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      USE PELIB_INTERFACE, ONLY: USE_PELIB
#include "implicit.h"
C
      DIMENSION OMEGA2(*),C2AM(*),WORK(LWORK)
      DIMENSION XLAMDP(*),XLAMDH(*),XLPC1(*),XLHC1(*)
C
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
#include "ccsections.h"
#include "ccfield.h"
#include "leinf.h"
#include "qm3.h"
      REAL*8, ALLOCATABLE :: FOCKMAT(:), FOCKTEMP(:)
!#include "qmmm.h"
C
      LOGICAL ETRAN,FCKCON
C
C      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC_FCC2FF')
C
      ONE = 1.0D0
C
      KFOCK1 = 1
      KFOCK2 = KFOCK1 + MAX(N2BST(ISYMOP),N2BST(ISYMLA))
      KEND1  = KFOCK2 + MAX(N2BST(ISYMOP),N2BST(ISYMLA))
      KEMAT1 = KEND1  
      KEMAT2 = KEMAT1 + NEMAT1(ISYMLA)
      KEND1  = KEMAT2 + NMATIJ(ISYMLA)
      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Needed:', KEND1, 'Available:', LWORK
         CALL QUIT('Insufficient memory for allocation in cc_fcc2ff')
      ENDIF
C
      CALL DZERO(WORK(KFOCK1),MAX(N2BST(ISYMOP),N2BST(ISYMLA)))
      CALL DZERO(WORK(KEMAT1),NEMAT1(ISYMLA))
      CALL DZERO(WORK(KEMAT2),NMATIJ(ISYMLA))
C
      DO 13 IF = 1, NFIELD
         FF =  EFIELD(IF)
         CALL CC_ONEP(WORK(KFOCK1),WORK(KEND1),LWRK1,FF,1,LFIELD(IF))
 13   CONTINUE
C
C
C-------------------------------------
C     Solvent contribution.
C     Put into one-electron integrals.
C SLV02,OC
C-------------------------------------
C
      IF (CCSLV .AND. (.NOT. CCMM )) THEN
         DTIME = SECOND()
         CALL CCSL_RHSTG(WORK(KFOCK1),WORK(KEND1),LWRK1)
         IF (IPRINT .GT. 5) THEN
           WRITE(LUPRI,*)'Time used in CCSL_RHSTG (CC_FCC2FF) 1',
     *                    SECOND() - DTIME
         END IF
      ENDIF
C
C-------------------------------------
C     Solvent contribution.
C     Put into one-electron integrals.
C CCMM02,JA+AO
C-------------------------------------
C
      IF (CCMM) THEN
         DTIME = SECOND()
         IF (.NOT. NYQMMM ) THEN
            CALL CCMM_RHSTG(WORK(KFOCK1),WORK(KEND1),LWRK1)
         ELSE IF (NYQMMM) THEN
            IF (HFFLD) THEN
              CALL CCMM_ADDGHF(WORK(KFOCK1),WORK(KEND1),LWRK1)
            ELSE
              CALL CCMM_ADDG(WORK(KFOCK1),WORK(KEND1),LWRK1)
            END IF
         END IF
         IF (IPRINT .GT. 5) THEN
           WRITE(LUPRI,*)'Time used in CCMM_RHSTG (CC_FCC2FF) 1',
     *                    SECOND() - DTIME
         END IF

      ENDIF
      IF (USE_PELIB()) THEN
          ALLOCATE(FOCKMAT(NNBASX),
     &             FOCKTEMP(MAX(N2BST(ISYMOP),N2BST(ISYMLA))))
          IF (HFFLD) THEN
              CALL GET_FROM_FILE('FOCKMHF',NNBASX,FOCKMAT)
          ELSE
          CALL GET_FROM_FILE('FOCKMAT',NNBASX,FOCKMAT)
          END IF
          CALL DSPTSI(NBAS,FOCKMAT,FOCKTEMP)
          CALL DAXPY(MAX(N2BST(ISYMOP),N2BST(ISYMLA)),1.0d0,FOCKTEMP,
     &               1,WORK(KFOCK1),1)
          DEALLOCATE(FOCKMAT,FOCKTEMP)
      END IF
C
C-------------------------------------
C
      CALL DCOPY(N2BST(ISYMOP),WORK(KFOCK1),1,WORK(KFOCK2),1)
      ISYFAO = 1
      ISYMPA = ISYMLA
      ISYMHO = 1
      CALL CC_FCKMO(WORK(KFOCK1),XLPC1,XLAMDH,
     *              WORK(KEND1),LWRK1,ISYFAO,ISYMPA,ISYMHO)
      ISYFAO = 1
      ISYMPA = 1
      ISYMHO = ISYMLA
      CALL CC_FCKMO(WORK(KFOCK2),XLAMDP,XLHC1,
     *              WORK(KEND1),LWRK1,ISYFAO,ISYMPA,ISYMHO)
      CALL DAXPY(N2BST(ISYMLA),ONE,WORK(KFOCK2),1,WORK(KFOCK1),1)
C
      ETRAN  = .FALSE.
      FCKCON = .TRUE.
      ISYMEI = ISYMLA

      CALL CCRHS_EFCK(WORK(KEMAT1),WORK(KEMAT2),XLAMDH,
     *                WORK(KFOCK1),WORK(KEND1),LWRK1,
     *                FCKCON,ETRAN,ISYMEI)
C
      IF (ISIDE .EQ. -1 ) THEN
        CALL CC_EITR(WORK(KEMAT1),WORK(KEMAT2),WORK(KEND1),LWRK1,
     *               ISYMEI)
      ENDIF
C
      CALL CCRHS_E(OMEGA2,C2AM,WORK(KEMAT1),WORK(KEMAT2),
     *             WORK(KEND1),LWRK1,ISYMTR,ISYMEI)
C
      CALL QEXIT('CC_FCC2FF')
C
      RETURN
      END
