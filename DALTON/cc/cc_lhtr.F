!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
C
C  /* Deck cc_lhtr */
      SUBROUTINE CC_LHTR(ECURR,
     *                   FRHO1,LUFR1,FRHO2,LUFR2,FRHO12,LUFR12,
     *                   FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                   RHO1,RHO2,CTR1,CTR2,WORK,LWORK,
     *                   NSIMTR,IVEC,ITR,LRHO1,APROXR12)
C
C     Written by Asger Halkier & Henrik Koch summer/fall 1995.
C
C     Version 1.1
C
C     Purpose:
C
C     Calculate left hand side transformation of a trial vector
C     in an AO-integral direct fashion. The biorthonormal basis
C     is used.
C
C     Generalised for transforming more than one vector per
C     integral calculation by Asger Halkier & Ove Christiansen
C     medio October 1996. New input list and new name 4 nov. 1996.
C     (IVEC is first number for first vector on files, FRHO1,FC1AM,FC2AM,
C           FRHO12,FC12AM;
C      ITR is first vector on FRHO2)
C     OC,CC2 FF bugfix, March 1997.
C     Some changes for CC2 and NONHF=.true. to allow for finite diff.
C     with respect to the orbtial coeff. vector. Ch. Haettig, march 2000
C     Some changes for CC2-R12 and NONHF=.true. Chr. Neiss, April 2005
C
C     Current models are: CC2, CCD, CCSD
C
C
C
      USE PELIB_INTERFACE, ONLY: USE_PELIB, PELIB_IFC_TRANSFORMER
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxorb.h"
#include "maxash.h"
#include "mxcent.h"
#include "aovec.h"
#include "iratdef.h"
#include "ccorb.h"
C
C   #include "infind.h" replaced by: #include <ccisao.h>
C   C. Haettig 13.08.2004
C
#include "ccisao.h"
#include "blocks.h"
#include "ccsdinp.h"
#include "ccsections.h"
#include "ccfield.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccinftap.h"
#include "distcl.h"
#include "cbieri.h"
#include "cclr.h"
#include "eritap.h"
#include "ccslvinf.h"
#include "qm3.h"
!#include "qmmm.h"
#include "ccnoddy.h"
#include "r12int.h"
#include "ccr12int.h"
C
      CHARACTER*10 MODEL
      CHARACTER*8  FRHO2, FRHO1, FC2AM, FC1AM, FN3SRT, FNTOC, FRHO12,
     & FC12AM
      CHARACTER*7  FN3FOP2X
      CHARACTER*6  CFIL, DFIL, PQFIL, FN3V, FNCKJD, FN3VI, FN3FOPX
      CHARACTER*5  FNDKBC3
      CHARACTER*4  FN4V, FNDKBC
      CHARACTER*3  APROXR12
      CHARACTER*1  LR
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (FOURTH = 0.25D0, THREE = 3.0D0)
      DIMENSION INDEXA(MXCORB_CC)
      DIMENSION RHO1(LRHO1,NSIMTR), CTR1(LRHO1,NSIMTR)
      DIMENSION RHO2(*), CTR2(*)
      DIMENSION WORK(LWORK)
      DIMENSION IADRPQ(MXCORB_CC,MAXSIM)
      REAL*8, ALLOCATABLE :: FOCKMAT(:), FOCKTEMP(:)
C
      LOGICAL ETRAN,FCKCON,CCSDR12,SKIPMI
C
      CALL QENTER('CC_LHTR')
C
      THIRD = ONE/THREE
C
C
      !set CCSDR12
      CCSDR12 = CCR12 .AND. .NOT.(CC2.OR.CCS)
C
      SKIPMI = .FALSE.
C
C---------------------------------
C     Initialze timing parameters.
C---------------------------------
C
      TIMTOT = ZERO
      TIMTOT = SECOND()
      TIMHE1 = ZERO
      TIMHE2 = ZERO
      TIRDAO = ZERO
      TIMIO  = ZERO
      TIMEYI = ZERO
      TIMEXI = ZERO
      TIMEMI = ZERO
      TIMETI = ZERO
      TIMEZ1 = ZERO
      TIMEZ2 = ZERO
      TIMFCK = ZERO
      TIMEC  = ZERO
      TIMEA  = ZERO
      TIMEH  = ZERO
      TIMEI  = ZERO
      TIME3O = ZERO
      TIME1O = ZERO
      TIMETP = ZERO
      TIMEBF = ZERO
      TIM2BF = ZERO
      TIMEEM = ZERO
      TIMEB  = ZERO
      TIMEG  = ZERO
      TIMEC2 = ZERO
      TIM12B = ZERO
      TIM12A = ZERO
      TIM11B = ZERO
      TIMEAM = ZERO
      TIM2EM = ZERO
      TIM22E = ZERO
      TIM22A = ZERO
      TIM22C = ZERO
      TIM22D = ZERO
      TIM2AO = ZERO
      TIM2MO = ZERO
C
C-----------------------------
C     Work-space allocation 1.
C-----------------------------
C
      ISYRES = MULD2H(ISYMTR,ISYMOP)
C
      KLAMDP = 1
      KLAMDH = KLAMDP + NLAMDT
      KC1T2  = KLAMDH + NLAMDT
      KCT1AO = KC1T2  + N2BST(ISYRES)*NSIMTR
      KDENSI = KCT1AO + NGLMRH(ISYMTR)*NSIMTR
      KFOCK  = KDENSI + N2BST(1)
      KFOCKG = KFOCK  + N2BST(ISYMOP)
      IF (CC2.AND.(.NOT.NONHF)) THEN
         KEND0 = KFOCKG + N2BST(ISYMTR)*NSIMTR
      ELSE
         KYMAT = KFOCKG + N2BST(ISYMTR)*NSIMTR
         KYDEN = KYMAT  + NMATAB(ISYRES)*NSIMTR
         KXMAT = KYDEN  + N2BST(ISYRES)*NSIMTR
         KEND0 = KXMAT  + NMATIJ(ISYRES)*NSIMTR
      ENDIF
C
      IF (CCR12) THEN
        KRHOR12 = KEND0
        KEND0 = KRHOR12 + NTR12SQ(ISYMTR)
      END IF
C
      IF (CC3) THEN
         KOVVO = KEND0
         KOOVV = KOVVO + NT2SQ(1)
         KEND0 = KOOVV + NT2SQ(1)
      ENDIF
C
      IF (CC2) THEN
         IF (NONHF) THEN
            KFCKHF = KEND0
            KEND0  = KFCKHF + N2BST(1)
         END IF
         KT2AM = KEND0
         KT1AM = KT2AM  + MAX(NT2AMX,NT2AM(ISYMOP))
         KEND1 = KT1AM  + MAX(NT1AMX,NT1AM(ISYMOP))
      ELSE
         KT2AM = KEND0
         KT1AM = KT2AM  + MAX(NT2AMX,NT2AM(ISYMOP),NT2AM(ISYRES))
         KEND1 = KT1AM  + MAX(NT1AMX,NT1AM(ISYMOP))
      ENDIF
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CC_LHTR')
      ENDIF
C
C-------------------------------
C     Open scratch file CCINT3O.
C-------------------------------
C
      LUIN30 = -1
      CALL WOPEN2(LUIN30,'CCINT3O',64,0)
C
C-------------------------------------------------------------
C     Open files with C- & D-intermediates needed in 22-block.
C-------------------------------------------------------------
C
      LUCIM = -1
      LUDIM = -1
C
      CFIL = 'PMAT_C'
      DFIL = 'PMAT_D'
C
      CALL WOPEN2(LUCIM,CFIL,64,0)
      CALL WOPEN2(LUDIM,DFIL,64,0)
C
C-------------------------------------------------------------
C     Open files with P- & Q-intermediates needed in 21-block.
C-------------------------------------------------------------
C
      IF (.TRUE.) THEN

         LUPQ = -1
         PQFIL = 'CCPQIM'
C
         CALL WOPEN2(LUPQ,PQFIL,64,0)

      END IF
C
C--------------------------------------------------
C     Open files needed for CC3 and initialise
C     intermediates allocated earlier
C--------------------------------------------------
C
      IF (CC3) THEN
C
         LU3SRT   = -1
         LUCKJD   = -1
         LUDKBC   = -1
         LUTOC    = -1
         LU3VI    = -1
         LUDKBC3  = -1
         LU3FOPX  = -1
         LU3FOP2X = -1
         LU3V     = -1
         LU4V     = -1
C
         FN3SRT   = 'CC3_SORT'
         FNCKJD   = 'CKJDEL'
         FNDKBC   = 'DKBC'
         FNTOC    = 'CCSDT_OC'
         FN3VI    = 'CC3_VI'
         FNDKBC3  = 'DKBC3'
         FN3FOPX  = 'PTFOPX'
         FN3FOP2X = 'PTFOP2X'
         FN3V     = 'ABCDEL'
         FN4V     = 'ABCD'
C
         CALL WOPEN2(LU3SRT,FN3SRT,64,0)
         CALL WOPEN2(LUCKJD,FNCKJD,64,0)
         CALL WOPEN2(LUDKBC,FNDKBC,64,0)
         CALL WOPEN2(LUTOC,FNTOC,64,0)
         CALL WOPEN2(LU3VI,FN3VI,64,0)
         CALL WOPEN2(LUDKBC3,FNDKBC3,64,0)
         CALL WOPEN2(LU3FOPX,FN3FOPX,64,0)
         CALL WOPEN2(LU3FOP2X,FN3FOP2X,64,0)
         IF (LVVVV) THEN
            CALL WOPEN2(LU4V,FN4V,64,0)
            CALL WOPEN2(LU3V,FN3V,64,0)
         END IF
C
         CALL DZERO(WORK(KOVVO),NT2SQ(1))
         CALL DZERO(WORK(KOOVV),NT2SQ(1))
C
      ENDIF
C-------------------------------------------
C     Read zero'th order cluster amplitudes.
C-------------------------------------------
C
      DTIME = SECOND()

      IOPT = 3
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),WORK(KT2AM))

      DTIME = SECOND() - DTIME
      TIMIO = TIMIO + DTIME
C
C----------------------------------------------------------------------
C     Initialize result-arrays and zero out single-vectors in CCD-calc.
C----------------------------------------------------------------------
C
      NRHO2 = MAX(NT2AM(ISYRES),2*NT2ORT(ISYRES))
      IF (CC2) NRHO2 = NT2AM(ISYRES)
C
      CALL DZERO(RHO1(1,1),NT1AM(ISYRES)*NSIMTR)
      CALL DZERO(RHO2,NRHO2)
C
      IF (CCD) THEN
         CALL DZERO(WORK(KT1AM),NT1AMX)
         CALL DZERO(CTR1(1,1),NT1AM(ISYMTR)*NSIMTR)
      ENDIF
C
      CALL DZERO(WORK(KFOCKG),N2BST(ISYMTR)*NSIMTR)
C
C----------------------------------
C     Calculate the lamda matrices.
C----------------------------------
C
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND1),
     &            LWRK1)
C
C------------------------------------------
C     Regain work space from T1-amplitudes.
C------------------------------------------
C
      KEND1 = KT1AM
      LWRK1 = LWORK - KEND1
C
C----------------------------------
C     Calculate the density matrix.
C----------------------------------
C
      ISYMH = 1
      IC    = 1
      CALL CC_AODENS(WORK(KLAMDP),WORK(KLAMDH),WORK(KDENSI),ISYMH,
     &                 IC,WORK(KEND1),LWRK1)
C
C--------------------------------------------------------
C        initialize start adress for P & Q intermediates:
C--------------------------------------------------------
C
      IADRSTRT = 1
C
C------------------------------------------------------------
C     Loop over trial vectors for constructing intermediates.
C------------------------------------------------------------
C
      DO 95 IV = 1,NSIMTR
C
C--------------------------------------------------------
C        Calculate contraction of CTR1 with lamda-matrix.
C--------------------------------------------------------
C
         KOFFAO = KCT1AO + NGLMRH(ISYMTR)*(IV - 1)
C
         CALL CCLT_CT1AO(CTR1(1,IV),WORK(KLAMDP),WORK(KOFFAO))
C
C-----------------------------------------------
C        Calculate contraction of CTR1 and T2AM.
C-----------------------------------------------
C
         KOFFCT = KC1T2 + N2BST(ISYRES)*(IV - 1)
C
         IOPT = 1
         CALL CC_C1T2C(WORK(KOFFCT),CTR1(1,IV),WORK(KT2AM),
     *                 WORK(KLAMDP),WORK(KLAMDH),WORK(KLAMDP),
     *                 WORK(KLAMDH),WORK(KEND1),LWRK1,ISYMTR,
     *                 ISYMOP,IOPT)
C
C--------------------------------------
C        Read CTR2 from disc into RHO2.
C--------------------------------------
C
         DTIME = SECOND()
         CALL CC_RVEC(LUFC2,FC2AM,NT2AM(ISYMTR),NT2AM(ISYMTR),
     *                IVEC+IV-1,RHO2)
         DTIME = SECOND() - DTIME
         TIMIO = TIMIO + DTIME
C
C-----------------------
C        Square up CTR2.
C-----------------------
C
         CALL CC_T2SQ(RHO2,CTR2,ISYMTR)
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYMTR),CTR1(1,IV),1,CTR1(1,IV),1)
            RHO2N = DDOT(NT2SQ(ISYMTR),CTR2,1,CTR2,1)
            WRITE(LUPRI,1) 'Norm of CTR1 -Read in before loop:  ',RHO1N
            WRITE(LUPRI,1) 'Norm of CTR2 -Read in before loop:  ',RHO2N
         ENDIF
C
         IF (.NOT. (CC2 .AND.(.NOT.NONHF))) THEN
C
C--------------------------------------------------
C        Calculate Y-intermediate of CTR2 and T2AM.
C--------------------------------------------------
C
            KOFFYI = KYMAT + NMATAB(ISYRES)*(IV - 1)
            KOFFYD = KYDEN + N2BST(ISYRES)*(IV - 1)
C
            TIMEYI = SECOND()
            CALL CC_YI(WORK(KOFFYI),CTR2,ISYMTR,WORK(KT2AM),ISYMOP,
     *                 WORK(KEND1),LWRK1)
            CALL CC_YD(WORK(KOFFYD),WORK(KOFFYI),ISYRES,WORK(KLAMDH),
     *                 WORK(KLAMDP),ISYMOP,WORK(KEND1),LWRK1)
            TIMEYI = SECOND() - TIMEYI
C
C-----------------------------------------------------
C           Calculate X-intermediate of CTR2 and T2AM.
C-----------------------------------------------------
C
            KOFFXI = KXMAT + NMATIJ(ISYRES)*(IV - 1)
C
            TIMEXI = SECOND()
            CALL CC_XI(WORK(KOFFXI),CTR2,ISYMTR,WORK(KT2AM),ISYMOP,
     *                 WORK(KEND1),LWRK1)
            TIMEXI = SECOND() - TIMEXI
         END IF
C
C------------------------------------------------------------
C           Precalculate P and Q intermediates
C           and restore CTR2 which is overwritten in CC_PQI:
C------------------------------------------------------------
C
         IF (.NOT. CC2) THEN
C
            IF ( DEBUG ) THEN
               RHO1N = DDOT(NT1AM(ISYMTR),CTR1(1,IV),1,CTR1(1,IV),1)
               RHO2N = DDOT(NT2SQ(ISYMTR),CTR2,1,CTR2,1)
               WRITE(LUPRI,1) 'Norm of CTR1 - before CC_PQI:  ',RHO1N
               WRITE(LUPRI,1) 'Norm of CTR2 - before CC_PQI:  ',RHO2N
               RHO2N = DDOT(NT2AM(ISYMOP),WORK(KT2AM),1,WORK(KT2AM),1)
               WRITE(LUPRI,1) 'Norm of T2AM - before CC_PQI:  ',RHO2N
            ENDIF

            CALL CC_PQI(CTR2,ISYMTR,WORK(KT2AM),ISYMOP,PQFIL,LUPQ,
     *                  IADRPQ,IADRSTRT,IV,WORK(KEND1),LWRK1)
C
            DTIME = SECOND()
            CALL CC_RVEC(LUFC2,FC2AM,NT2AM(ISYMTR),NT2AM(ISYMTR),
     *                   IVEC+IV-1,RHO2)
            DTIME = SECOND() - DTIME
            TIMIO = TIMIO + DTIME
C
            CALL CC_T2SQ(RHO2,CTR2,ISYMTR)
C
            IF ( DEBUG ) THEN
               RHO1N = DDOT(NT1AM(ISYMTR),CTR1(1,IV),1,CTR1(1,IV),1)
               RHO2N = DDOT(NT2SQ(ISYMTR),CTR2,1,CTR2,1)
               WRITE(LUPRI,1) 'Norm of CTR1 - before CC_PQI:  ',RHO1N
               WRITE(LUPRI,1) 'Norm of CTR2 - before CC_PQI:  ',RHO2N
               RHO2N = DDOT(NT2AM(ISYMOP),WORK(KT2AM),1,WORK(KT2AM),1)
               WRITE(LUPRI,1) 'Norm of T2AM - before CC_PQI:  ',RHO2N
            ENDIF

         ENDIF
C
  95  CONTINUE
C
C
C---------------------------------------------------------------
C     Calculate transposed CTR2 if t2tcor and reinitialize RHO2.
C---------------------------------------------------------------
C
      IF (MINSCR) THEN
C
         CALL DZERO(RHO2,NRHO2)
C
         IF (((.NOT. DIRECT) .AND. T2TCOR) .AND. (.NOT. CC2)) THEN
C
            KCTR2T = KEND1
            KEND1  = KCTR2T + NT2SQ(ISYMTR)
            LWRK1  = LWORK  - KEND1
            IF (LWRK1 .LT. 0) THEN
               CALL QUIT('Insufficient core in CCRHSN')
            END IF
C
            JSYM = ISYMTR
            CALL DCOPY(NT2SQ(ISYMTR),CTR2,1,WORK(KCTR2T),1)
            AUTIME = SECOND()
            CALL CCSD_T2TP(WORK(KCTR2T),WORK(KEND1),LWRK1,JSYM)
            TIMETP = SECOND() - AUTIME
C
         ENDIF
      ENDIF
C
C-------------------------------------------------------------
C     Regain work space from T2-amplitudes in CC2-calculation.
C-------------------------------------------------------------
C
      IF (CC2) THEN
C
         KEND1 = KT2AM
         LWRK1 = LWORK - KEND1
C
      ENDIF
C
C-----------------------------------
C     Start the loop over integrals.
C-----------------------------------
C
      KENDS2 = KEND1
      LWRKS2 = LWRK1
C
      IF (DIRECT) THEN
         DTIME  = SECOND()
         IF (HERDIR) THEN
            CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
         ELSE
            KCCFB1 = KEND1
            KINDXB = KCCFB1 + MXPRIM*MXCONT
            KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
            LWRK1  = LWORK  - KEND1
            CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     *                  KODPP1,KODPP2,KRDPP1,KRDPP2,
     *                  KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     *                  WORK(KEND1),LWRK1,IPRERI)
            KEND1 = KFREE
            LWRK1 = LFREE
         ENDIF
         DTIME  = SECOND() - DTIME
         TIMHE1 = TIMHE1 + DTIME
         NTOSYM = 1
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      KENDSV = KEND1
      LWRKSV = LWRK1
C
      ICDEL1 = 0
      DO 100 ISYMD1 = 1,NTOSYM
C
         IF (DIRECT) THEN
            IF (HERDIR) THEN
               NTOT = MAXSHL
            ELSE
               NTOT = MXCALL
            ENDIF
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
C
         DO 110 ILLL = 1,NTOT
C
C-----------------------------------------------------------------
C           If direct calculate the integrals and transposed CTR2.
C-----------------------------------------------------------------
C
            IF (DIRECT) THEN
C
               KEND1 = KENDSV
               LWRK1 = LWRKSV
C
               DTIME  = SECOND()
               IF (HERDIR) THEN
                  CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                        IPRINT)
               ELSE
                  CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     *                        WORK(KODCL1),WORK(KODCL2),
     *                        WORK(KODBC1),WORK(KODBC2),
     *                        WORK(KRDBC1),WORK(KRDBC2),
     *                        WORK(KODPP1),WORK(KODPP2),
     *                        WORK(KRDPP1),WORK(KRDPP2),
     *                        WORK(KCCFB1),WORK(KINDXB),
     *                        WORK(KEND1),LWRK1,IPRERI)
               ENDIF
               DTIME  = SECOND() - DTIME
               TIMHE2 = TIMHE2 + DTIME
C
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient core in CC_LHTR')
               END IF
C
               IF ((MINSCR .AND. T2TCOR) .AND. (.NOT. CC2)) THEN
C
                  KCTR2T = KEND1
                  KEND1  = KCTR2T + NT2SQ(ISYMTR)
                  LWRK1  = LWORK  - KEND1
                  IF (LWRK1 .LT. 0) THEN
                     CALL QUIT('Insufficient core in CC_LHTR')
                  END IF
C
                  JSYM = ISYMTR
                  CALL DCOPY(NT2SQ(ISYMTR),CTR2,1,WORK(KCTR2T),1)
                  AUTIME = SECOND()
                  CALL CCSD_T2TP(WORK(KCTR2T),WORK(KEND1),LWRK1,JSYM)
                  AUTIME = SECOND() - AUTIME
                  TIMETP = TIMETP + AUTIME
C
               END IF
C
            ELSE
               NUMDIS = 1
            ENDIF
C
C----------------------------------------------------
C           Loop over number of trial vectors nsimtr.
C----------------------------------------------------
C
            DO 115 IV = 1, NSIMTR
C
               IF (.NOT. MINSCR) THEN
C
C-------------------------------------------------------------
C                 Read CTR2 from disc into RHO2 and square up.
C-------------------------------------------------------------
C
                  DTIME = SECOND()
                  CALL CC_RVEC(LUFC2,FC2AM,NT2AM(ISYMTR),NT2AM(ISYMTR),
     *                         IVEC+IV-1,RHO2)
                  DTIME = SECOND() - DTIME
                  TIMIO = TIMIO + DTIME
C
                  CALL CC_T2SQ(RHO2,CTR2,ISYMTR)
C
C----------------------------------------------
C                 Read result vector from disc.
C----------------------------------------------
C
                  DTIME = SECOND()
                  CALL CC_RVEC(LUFR2,FRHO2,NRHO2,NRHO2,
     *                         IV+ITR-1,RHO2)
                  DTIME = SECOND() - DTIME
                  TIMIO = TIMIO + DTIME
C
C-----------------------------------------------------
C                 Calculate transposed CTR2 if t2tcor.
C-----------------------------------------------------
C
                  IF (T2TCOR) THEN
C
                     KCTR2T = KEND1
                     KEND1  = KCTR2T + NT2SQ(ISYMTR)
                     LWRK1  = LWORK  - KEND1
                     IF (LWRK1 .LT. 0) THEN
                        CALL QUIT('Insufficient core in CC_LHTR')
                     END IF
C
                     JSYM = ISYMTR
                     CALL DCOPY(NT2SQ(ISYMTR),CTR2,1,WORK(KCTR2T),1)
                     AUTIME = SECOND()
                     CALL CCSD_T2TP(WORK(KCTR2T),WORK(KEND1),
     *                              LWRK1,JSYM)
                     AUTIME = SECOND() - AUTIME
                     TIMETP = TIMETP + AUTIME
C
                  ENDIF
               ENDIF
C
C--------------------------------------------------------------
C           Calculate adresses for CTR-dependent intermediates.
C--------------------------------------------------------------
C
               KOFFAO = KCT1AO + NGLMRH(ISYMTR)*(IV - 1)
               KOFFCT = KC1T2  + N2BST(ISYRES)*(IV - 1)
               KOFFYI = KYMAT  + NMATAB(ISYRES)*(IV - 1)
               KOFFYD = KYDEN  + N2BST(ISYRES)*(IV - 1)
               KOFFFG = KFOCKG + N2BST(ISYMTR)*(IV - 1)
C
C-----------------------------------------------------
C           Loop over number of distributions in disk.
C-----------------------------------------------------
C
            DO 120 IDEL2 = 1,NUMDIS
C
               IF (DIRECT) THEN
                  IDEL  = INDEXA(IDEL2)
                  IF (NOAUXB) THEN
                     IDUM = 1
                     CALL IJKAUX(IDEL,IDUM,IDUM,IDUM)
                  END IF
                  ISYMD = ISAO(IDEL)
               ELSE
                  IDEL  = IBAS(ISYMD1) + ILLL
                  ISYMD = ISYMD1
               ENDIF
C
               ISYDIS = MULD2H(ISYMD,ISYMOP)
C
C------------------------------------------
C              Work space allocation no. 2.
C------------------------------------------
C
               KXINT  = KEND1
               KEND2  = KXINT + NDISAO(ISYDIS)
               LWRK2  = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CC_LHTR')
               ENDIF
C
C--------------------------------------------
C              Read AO integral distribution.
C--------------------------------------------
C
               AUTIME = SECOND()
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
               AUTIME = SECOND() - AUTIME
               TIRDAO = TIRDAO + AUTIME
C
C-------------------------------------------
C              Calculate the AO-fock-matrix.
C-------------------------------------------
C
               AUTIME = SECOND()
               ISYDEN = ISYMTR
               CALL CC_AOFOCK(WORK(KXINT),WORK(KOFFCT),WORK(KOFFFG),
     *                        WORK(KEND2),LWRK2,IDEL,ISYMD,.FALSE.,
     *                        DUMMY,ISYDEN)
               AUTIME = SECOND() - AUTIME
               TIMFCK = TIMFCK + AUTIME
C
C------------------------------------------
C              Work space allocation no. 3.
C------------------------------------------
C
               ISYMTI = MULD2H(ISYMD,ISYMOP)
               ISY21I = MULD2H(ISYMD,ISYRES)
C
               KDSRHF = KEND2
               K3OINT = KDSRHF + NDSRHF(ISYMD)
               IF (CC2) THEN
                  KEND3  = K3OINT + NMAIJK(ISYMTI)
               ELSE
                  KSCRTI = K3OINT + NMAIJK(ISYMTI)
                  KSCRZI = KSCRTI + NT2BCD(ISYMTI)
                  KSCRWI = KSCRZI + NT2BCD(ISY21I)
                  KSCRVI = KSCRWI + NT2BCD(ISY21I)
                  KEND3  = KSCRVI + NT2BCD(ISY21I)
               ENDIF
               LWRK3  = LWORK  - KEND3
C
               IF (LWRK3 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND3,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CC_LHTR')
               ENDIF
C
C--------------------------------------------------------
C              Transform one index in the integral batch.
C--------------------------------------------------------
C
               AUTIME = SECOND()
               ISYMLP = 1
               CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KLAMDP),
     *                     ISYMLP,WORK(KEND3),LWRK3,ISYDIS)
               AUTIME = SECOND() - AUTIME
               TIME1O = TIME1O + AUTIME
C
C-------------------------------------------------------------------
C              Calculate integral batch with three occupied indices.
C-------------------------------------------------------------------
C
               AUTIME = SECOND()
               CALL CC_INT3O(WORK(K3OINT),WORK(KDSRHF),WORK(KLAMDH),
     *                       ISYMOP,WORK(KLAMDP),WORK(KEND3),LWRK3,
     *                       IDEL,ISYMD,LUIN30,'CCINT3O')
               AUTIME = SECOND() - AUTIME
               TIME3O = TIME3O + AUTIME
C
C--------------------------------------------------------------------
C              Calculate integrals for CC3.
C              Integrals with 3 virtual indices are stored on disc
C--------------------------------------------------------------------
C
               IF (CC3) THEN
                  IF (LVVVV) THEN
                   CALL CC3_INTDEL(WORK(KXINT),1,LU3V,FN3V,WORK(KLAMDP),
     *                             1,WORK(KLAMDH),1,1,
     *                             WORK(KEND3),LWRK3,IDEL,ISYMD)
                  END IF
C
                  CALL CC3_2O2V(WORK(KXINT),1,WORK(KDSRHF),1,
     *                          WORK(KOVVO),WORK(KOOVV),
     *                          WORK(KLAMDP),1,WORK(KLAMDH),1,
     *                          WORK(KLAMDP),1,WORK(KLAMDH),1,1,
     *                          WORK(KEND3),LWRK3,IDEL,ISYMD)
               ENDIF
C
C----------------------------------------------------------
C              Calculate the LT12B term in CC2-calculation.
C----------------------------------------------------------
C
               IF (CC2) THEN
C
                  AUTIME = SECOND()
                  CALL CC_12B(RHO2,WORK(KDSRHF),WORK(KOFFAO),
     *                        ISYRES,WORK(KLAMDH),WORK(KEND3),
     *                        LWRK3,IDEL,ISYMD)
                  AUTIME = SECOND() - AUTIME
                  TIM12B = TIM12B + AUTIME
C
               ENDIF
C
C--------------------------------------------------------------
C              Calculate intermediates needed for the 21-block.
C--------------------------------------------------------------
C
               IF (.TRUE.) THEN

               IF (.NOT. CC2) THEN

                  IADR = IADRPQ(IDEL,IV)
                  CALL GETWA2(LUPQ,PQFIL,WORK(KSCRZI),
     *                        IADR,NT2BCD(ISY21I))

                  IADR = IADRPQ(IDEL,IV) + NT2BCD(ISY21I)
                  CALL GETWA2(LUPQ,PQFIL,WORK(KSCRVI),
     *                        IADR,NT2BCD(ISY21I))

C                 CALL DSCAL(NT2BCD(ISY21I),TWO,WORK(KSCRVI),1)
C                 CALL DSCAL(NT2BCD(ISY21I),TWO,WORK(KSCRZI),1)
                  CALL DZERO(WORK(KSCRWI),NT2BCD(ISY21I))

                  IF ( DEBUG ) THEN
                   PN=DDOT(NT2BCD(ISY21I),WORK(KSCRZI),1,WORK(KSCRZI),1)
                   QN=DDOT(NT2BCD(ISY21I),WORK(KSCRVI),1,WORK(KSCRVI),1)
                   WRITE(LUPRI,*) 'IV,IDEL,IADR,P,Q',IV,IDEL,IADR,PN,QN
                  ENDIF

               END IF

               ELSE
               IF (.NOT. CC2) THEN
C
                  AUTIME = SECOND()
                  CALL CC_TI(WORK(KSCRTI),ISYMTI,WORK(KT2AM),ISYMOP,
     *                         WORK(KLAMDH),1,WORK(KEND3),LWRK3,
     *                         IDEL,ISYMD)
                  AUTIME = SECOND() - AUTIME
                  TIMETI = TIMETI + AUTIME
C
                  IOPT = 1
C
                  AUTIME = SECOND()
                  CALL CC_ZWVI(WORK(KSCRZI),CTR2,ISYMTR,WORK(KSCRTI),
     *                         ISYMTI,WORK(KEND3),LWRK3,IOPT)
                  AUTIME = SECOND() - AUTIME
                  TIMEZ1 = TIMEZ1 + AUTIME
C
                  IOPT = 2
C
                  AUTIME = SECOND()
                  CALL CC_ZWVI(WORK(KSCRWI),CTR2,ISYMTR,WORK(KSCRTI),
     *                         ISYMTI,WORK(KEND3),LWRK3,IOPT)
                  AUTIME = SECOND() - AUTIME
                  TIMEZ2 = TIMEZ2 + AUTIME
C
                  IF (.NOT. T2TCOR) THEN
C
                     AUTIME = SECOND()
                     CALL CCSD_T2TP(CTR2,WORK(KEND3),LWRK3,ISYMTR)
                     AUTIME = SECOND() - AUTIME
                     TIMETP = TIMETP + AUTIME
C
                     IOPT = 2
C
                     AUTIME = SECOND()
                     CALL CC_ZWVI(WORK(KSCRVI),CTR2,ISYMTR,WORK(KSCRTI),
     &                            ISYMTI,WORK(KEND3),LWRK3,IOPT)
                     AUTIME = SECOND() - AUTIME
                     TIMEZ2 = TIMEZ2 + AUTIME
C
                     AUTIME = SECOND()
                     CALL CCSD_T2TP(CTR2,WORK(KEND3),LWRK3,ISYMTR)
                     AUTIME = SECOND() - AUTIME
                     TIMETP = TIMETP + AUTIME
C
                  ELSE IF (T2TCOR) THEN
C
                     IOPT = 2
C
                     AUTIME = SECOND()
                     CALL CC_ZWVI(WORK(KSCRVI),WORK(KCTR2T),ISYMTR,
     &                            WORK(KSCRTI),ISYMTI,WORK(KEND3),
     &                            LWRK3,IOPT)
                     AUTIME = SECOND() - AUTIME
                     TIMEZ2 = TIMEZ2 + AUTIME
C
                  ENDIF
C
               ENDIF
               ENDIF
C
C--------------------------------------------------
C                 Calculate the LT21I contribution.
C--------------------------------------------------
C
               IF (.NOT. CC2) THEN

                  IOPT   = 1
                  ISYLRD = MULD2H(ISYMD,ISYRES)
                  AUTIME = SECOND()
                  CALL CC_21I2(RHO1(1,IV),WORK(KXINT),ISYDIS,DUMMY,0,
     *                         WORK(KSCRZI),WORK(KSCRVI),ISYLRD,
     *                         DUMMY,       DUMMY,       0,
     *                         WORK(KLAMDP),WORK(KLAMDH),ISYMOP,
     *                         WORK(KLAMDP),ISYMOP,
     *                         WORK(KEND3),LWRK3,IOPT,.FALSE.,.FALSE.,
     *                         .FALSE.)
                  AUTIME = SECOND() - AUTIME
                  TIMEI = TIMEI + AUTIME

                  IF ( DEBUG ) THEN
                    RHO1N=DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
                    WRITE(LUPRI,*) 'Norm of RHO1 after CC_21I2:',IV,
     &                   RHO1N,IDEL
                  ENDIF

               ENDIF
C
C--------------------------------------------------
C                 Calculate the LT21A contribution.
C--------------------------------------------------
C
               IF (.NOT. CC2) THEN
                  AUTIME = SECOND()
                  CALL CC_21A(RHO1(1,IV),WORK(KDSRHF),WORK(KOFFYI),
     *                        ISYRES,WORK(KOFFYD),ISYRES,WORK(KLAMDH),
     *                        WORK(KLAMDP),ISYMOP,
     *                        WORK(KEND3),LWRK3,IDEL,ISYMD)
                  AUTIME = SECOND() - AUTIME
                  TIMEA = TIMEA + AUTIME

                  IF ( DEBUG ) THEN
                    RHO1N=DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
                    WRITE(LUPRI,*) 'Norm of RHO1 after CC_21A:',IV,RHO1N
                  ENDIF

               ENDIF
C
C--------------------------------------------------
C                 Calculate the LT21H contribution.
C--------------------------------------------------
C
               IF (.NOT. CC2) THEN
                  AUTIME = SECOND()
                  CALL CC_21H(RHO1(1,IV),ISYRES,WORK(KSCRVI),
     *                        WORK(KSCRWI),WORK(KSCRZI),ISYLRD,
     *                        WORK(K3OINT),ISYMOP,WORK(KEND3),
     *                        LWRK3,ISYMD)
                  AUTIME = SECOND() - AUTIME
                  TIMEH = TIMEH + AUTIME

                  IF ( DEBUG ) THEN
                    RHO1N=DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
                    WRITE(LUPRI,*) 'Norm of RHO1 after CC_21H:',IV,RHO1N
                  ENDIF

               ENDIF
C
C------------------------------------------
C              Work space allocation no. 4.
C------------------------------------------
C
               ISSCRM = MULD2H(ISYMD,ISYMTR)
C
               KSCRM = KEND3
               KEND4 = KSCRM  + NT2BCD(ISSCRM)
               LWRK4 = LWORK  - KEND4
C
               IF (LWRK4 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND4,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CC_LHTR')
               ENDIF
C
C--------------------------------------------------------------
C              Construct the partially transformed CTR2-vector.
C--------------------------------------------------------------
C
               AUTIME = SECOND()
               ICON = 1
               ISYMLP = 1
               CALL CC_T2AO(CTR2,WORK(KLAMDP),ISYMLP,WORK(KSCRM),
     *                        WORK(KEND4),LWRK4,IDEL,ISYMD,
     *                        ISYMTR,ICON)
               AUTIME = SECOND() - AUTIME
               TIM2AO = TIM2AO + AUTIME
C
C-------------------------------------------------------
C              Calculate the LT21C- and D contributions.
C-------------------------------------------------------
C
               AUTIME = SECOND()
               IOPT  = 1
               IF ( CC2 ) THEN
                 ICON  = 2
               ELSE
                 ICON  = 1
               ENDIF
C
               ISYMM = MULD2H(ISYMD,ISYMTR)
               CALL CC_21DC(RHO1(1,IV),CTR2,ISYMTR,WORK(KSCRM),ISYMM,
     *                      WORK(KSCRM),ISYMM,WORK(KXINT),
     *                      WORK(KLAMDH),ISYMOP,WORK(KLAMDP),ISYMOP,
     *                      WORK(KLAMDH),ISYMOP,WORK(KLAMDP),ISYMOP,
     *                      WORK(KEND4),LWRK4,IDEL,ISYMD,IOPT,ICON)
               AUTIME = SECOND() - AUTIME
               TIMEC  = TIMEC + AUTIME

                  IF ( DEBUG ) THEN
                    RHO1N=DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
                    WRITE(LUPRI,*) 'Norm of RHO1 after CC_21DC:',
     &                   IV,RHO1N
                  ENDIF

C
C--------------------------------------------------------
C              Calculate the LT12B & LT22B contributions.
C--------------------------------------------------------
C
               IF (.NOT. CC2) THEN
C
                  AUTIME = SECOND()
                  IOPT   = 2
                  ISYMM  = MULD2H(ISYMD,ISYMTR)
                  CALL CC_BF(WORK(KXINT),RHO2,WORK(KLAMDP),1,
     *                         WORK(KOFFAO),ISYMTR,WORK(KLAMDP),1,
     *                         WORK(KSCRM),ISYMM,DUMMY,1,WORK(KEND4),
     *                         LWRK4,IDEL,ISYMD,IOPT)
                  AUTIME = SECOND() - AUTIME
                  TIM2BF = TIM2BF + AUTIME
C
               ENDIF
C
  120       CONTINUE
C
               IF (.NOT. MINSCR) THEN
C
C-----------------------------------------
C                 Write out result vector.
C-----------------------------------------
C
                  DTIME   = SECOND()
                  CALL CC_WVEC(LUFR2,FRHO2,NRHO2,NRHO2,IV + ITR -1,
     *                         RHO2)
                  DTIME   = SECOND() - DTIME
                  TIMIO   = TIMIO    + DTIME
C
               ENDIF
C
  115    CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C------------------------
C     Recover work space.
C------------------------
C
      KEND1 = KENDS2
      LWRK1 = LWRKS2
C
C--------------------------------------------------------------
C     calculate R12-contribution to rho_BZ (<--> LT12B & LT22B)
C--------------------------------------------------------------
C
      IF (CCSDR12) THEN
        DO IV = 1, NSIMTR
          IF (.NOT. MINSCR) THEN
            DTIME = SECOND()
            CALL CC_RVEC(LUFR2,FRHO2,NRHO2,NRHO2,IV + ITR -1,
     *                   RHO2)
            DTIME = SECOND() - DTIME
            TIMIO = TIMIO + DTIME
          END IF
C  
          !calculate B''-contribution:
          IOPT = 0
          IAMP = 1
          IFILE = IVEC+IV-1
          CALL CCRHS_BP(RHO2,ISYMTR,IOPT,IAMP,FC12AM,LUFC12,
     &                  IFILE,DUMMY,IDUMMY,2.0D0*BRASCL,
     &                  WORK(KEND1),LWRK1)
C 
          IF (.NOT. MINSCR) THEN
            DTIME = SECOND()
            CALL CC_WVEC(LUFR2,FRHO2,NRHO2,NRHO2,IV + ITR -1,
     *                   RHO2)
            DTIME = SECOND() - DTIME
            TIMIO = TIMIO + DTIME
          END IF
        END DO
      END IF
C
C-----------------------------
C     Loop over trial vectors.
C-----------------------------
C
      DO 125 IV = 1,NSIMTR
C
        IF (CCR12) CALL DZERO(WORK(KRHOR12),NTR12SQ(ISYMTR))
C
C--------------------------------------------------------------------
C     for CC2 and NONHF=.true. calculate Fock matrix entering E-term:
C       the SCF Fock matrix is in principle given by the SCF orbital
C       energies, but in recomputing it here from the SCF AO-Fock
C       matrix computed in CCSD_IAJB allows to do finite difference
C       on the vector function with respect to the CMO vector
C       (see CC_FDXI & CC_FDETA). Note the SCF AO-Fock matrix read
C       from file includes the `relaxed' external fields, so we
C       only have to add the unrelaxed fields.
C--------------------------------------------------------------------
C
         IF ((CC2.OR.CCR12) .AND. NONHF) THEN
           KCMO   = KEND1
           KFIELD = KCMO   + MAX(NLAMDT,NLAMDS)
           KEND2  = KFIELD + N2BAST
           IF (CCR12) THEN
             if (isymop.ne.1) call quit('Symmetry problem in CC_LHTR')
             kvxintsq = kend2
             kxint    = kvxintsq + nr12r12sq(isymop)
             kxintsq  = kxint + nr12r12p(1)
             kctr12sq = kxintsq + nr12r12sq(1)
             ktr12sq  = kctr12sq + ntr12sq(isymtr)
             kscr     = ktr12sq + ntr12sq(1)
             kend2    = kscr + ntr12sq(isymtr)
           END IF
           LWRK2  = LWORK  - KEND2
           IF (LWRK2 .LT. 0) THEN
             CALL QUIT('Insufficient memory in CCRHSN.')
           END IF

           CALL DZERO(WORK(KFIELD),N2BST(1))
           IF (CCR12) THEN
             CALL DZERO(WORK(KVXINTSQ),NR12R12SQ(1))
           END IF
           DO IF = 1, NFIELD
C            WRITE(LUPRI,*) IF, NHFFIELD(IF), EFIELD(IF), LFIELD(IF)
             IF ( NHFFIELD(IF) ) THEN
               CALL CC_ONEP(WORK(KFIELD),WORK(KEND2),LWRK2,EFIELD(IF),1,
     *                      LFIELD(IF))
               IF (CCR12) THEN
                 CALL CC_R12RDVXINT(WORK(KVXINTSQ),WORK(KEND2),LWRK2,
     &                            EFIELD(IF),1, LFIELD(IF))
               END IF
             ELSE IF (.NOT. NHFFIELD(IF) .AND. CCR12) THEN
               CALL QUIT('CCR12 response can only handle unrelaxed '//
     &                   'orbitals (w.r.t. the perturbation)')
             END IF
           END DO

           LUSIFC = -1
           CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     *                 IDUMMY,.FALSE.)
           REWIND(LUSIFC)
           CALL MOLLAB('TRCCINT ',LUSIFC,LUERR)
           READ(LUSIFC)
           READ(LUSIFC)
           READ(LUSIFC) (WORK(KCMO+I-1),I=1,NLAMDS)
           CALL GPCLOSE(LUSIFC,'KEEP')

           CALL CMO_REORDER(WORK(KCMO),WORK(KEND2),LWRK2)

           IF (CCR12) THEN
              ! read R12 trial vector and reorder to full square
              ifile = IV+IVEC-1
              call cc_r12getct(work(kctr12sq),isymtr,1,2.0D0*brascl,
     &                         .false.,'N',lufc12,fc12am,ifile,
     &                         cdummy,idummy,work(kend2),lwrk2)
C
              ! read R12 overlap matrix and reorder to full square
              lunit = -1
              call gpopen(lunit,fccr12x,'old',' ','unformatted',idummy,
     &                    .false.)
              rewind(lunit)
 8888         read(lunit) ian
              read(lunit) (work(kxint-1+i), i=1, nr12r12p(1))
              if (ian.ne.ianr12) goto 8888 
              call gpclose(lunit,'KEEP')
              iopt = 2
              call ccr12unpck2(work(kxint),1,work(kxintsq),'N',iopt)
C
              ! calculate R12 response contributions to Omega_{kilj}:
              ! first contribution: 
              CALL CC_R12XI(work(kscr),isymtr,'T',work(kctr12sq),isymtr,
     &                      work(kxintsq),work(kvxintsq),isymop,
     &                      work(kfield),work(klamdh),work(klamdp),'N',
     &                      work(kend2),lwrk2)
C
C             ! transpose: in CC_R12XI the r12-pair index (kl) is leading,
C             ! in RHOR12 the occ. index pair (ij) is leading!!! 
C             call cclr_trsqr12(work(kscr),isymtr)
C
C             ! add it to RHOR12 
              call daxpy(ntr12sq(isymtr),one,work(kscr),1,work(krhor12),
     &                   1)
C 
              ! second contribution:
              ! need R12 amplitudes
              CALL cc_r12getct(work(ktr12sq),1,0,dummy,.false.,'N',
     &                         idummy,cdummy,idummy,cdummy,
     &                         idummy,work(kend2),lwrk2)
C 
              call cc_r12etaa(RHO1(1,IV),ISYMTR,work(kctr12sq),isymtr,
     &                       work(ktr12sq),1,work(kxintsq),work(kfield),
     &                       1,work(kcmo),work(kcmo),.false.,
     &                       work(kend2),lwrk2)
 
           END IF
C
           IF (CC2) THEN
             LUFCK = -1
             CALL GPOPEN(LUFCK,'CC_FCKREF','UNKNOWN',' ','UNFORMATTED',
     *                   IDUMMY,.FALSE.)
             REWIND(LUFCK)
             READ(LUFCK)(WORK(KFCKHF + I-1),I = 1,N2BST(ISYMOP))
             CALL GPCLOSE(LUFCK,'KEEP' )

             ! SCF  Fock matrix in transformed using CMO vector
             CALL CC_FCKMO(WORK(KFCKHF),WORK(KLAMDP),WORK(KLAMDH),
     &                  WORK(KEND2),LWRK2,1,1,1)

C-------------------------------------
C     Solvent contribution.
C     Put into one-electron integrals.
C SLV98,OC
C-------------------------------------
C
C
        IF (CCSLV .AND. (.NOT. CCMM )) THEN
           DTIME = SECOND()
           CALL CCSL_RHSTG(WORK(KFIELD),WORK(KEND2),LWRK2)
c          TIMMLT = TIMMLT + SECOND() - DTIME
           IF (IPRINT .GT. 5) THEN
             WRITE(LUPRI,*)'Time used in CCSL_RHSTG (CC_LHTR) 1',
     *                      SECOND() - DTIME
           END IF
        ENDIF
C
C
C-------------------------------------
C     Solvent contribution.
C     Put into one-electron integrals.
C CCMM02,JA+AO
C-------------------------------------
C
        IF (CCMM) THEN
           DTIME = SECOND()
           IF (.NOT. NYQMMM) THEN
              CALL CCMM_RHSTG(WORK(KFIELD),WORK(KEND2),LWRK2)
           ELSE IF (NYQMMM) THEN
              IF (HFFLD) THEN  !Note we are inside cc2 
                 CALL CCMM_ADDGHF(WORK(KFIELD),WORK(KEND2),LWRK2)
              ELSE
                 CALL CCMM_ADDG(WORK(KFIELD),WORK(KEND2),LWRK2)
              ENDIF
           END IF
c          TIMMLT = TIMMLT + SECOND() - DTIME
           IF (IPRINT .GT. 5) THEN
             WRITE(LUPRI,*)'Time used in CCMM_RHSTG (CC_LHTR) 1',
     *                      SECOND() - DTIME
           END IF
        ENDIF
C
        IF (USE_PELIB()) THEN
            ALLOCATE(FOCKMAT(NNBASX),FOCKTEMP(N2BAST))
            IF (HFFLD) THEN
                CALL GET_FROM_FILE('FOCKMHF',NNBASX,FOCKMAT)
            ELSE
                CALL GET_FROM_FILE('FOCKMAT',NNBASX,FOCKMAT)
            END IF
            CALL DSPTSI(NBAS,FOCKMAT,FOCKTEMP)
            CALL DAXPY(N2BAST,1.0d0,FOCKTEMP,1,WORK(KFIELD),1)
            DEALLOCATE(FOCKMAT,FOCKTEMP)
        END IF
C
C====================================================
C-------------------------------------

             ! unrelaxed fields are transformed using the Lambda matrices
             CALL CC_FCKMO(WORK(KFIELD),WORK(KLAMDP),WORK(KLAMDH),
     *                     WORK(KEND2),LWRK2,1,1,1)

             CALL DAXPY(N2BST(1),ONE,WORK(KFIELD),1,WORK(KFCKHF),1)
           END IF 
         END IF
C
         IF (.NOT. MINSCR) THEN
C
C------------------------------
C           Recover work space.
C------------------------------
C
            KEND1 = KENDS2
            LWRK1 = LWRKS2
C
C-------------------------------------------------------
C           Read CTR2 from disc into RHO2 and square up.
C-------------------------------------------------------
C
            DTIME = SECOND()
            CALL CC_RVEC(LUFC2,FC2AM,NT2AM(ISYMTR),NT2AM(ISYMTR),
     *                   IVEC+IV-1,RHO2)
            DTIME = SECOND() - DTIME
            TIMIO = TIMIO + DTIME
C
            CALL CC_T2SQ(RHO2,CTR2,ISYMTR)
C
            IF ( DEBUG ) THEN
               RHO1N = DDOT(NT1AM(ISYMTR),CTR1(1,IV),1,CTR1(1,IV),1)
               RHO2N = DDOT(NT2SQ(ISYMTR),CTR2,1,CTR2,1)
               WRITE(LUPRI,1) 'Norm of CTR1 -Read in after loop:  ',
     &              RHO1N
               WRITE(LUPRI,1) 'Norm of CTR2 -Read in after loop:  ',
     &              RHO2N
            ENDIF
C
C
C----------------------------------------
C           Read result vector from disc.
C----------------------------------------
C
            DTIME = SECOND()
            CALL CC_RVEC(LUFR2,FRHO2,NRHO2,NRHO2,IV+ITR-1,RHO2)
            DTIME = SECOND() - DTIME
            TIMIO = TIMIO + DTIME
C
            IF (.NOT. CC2) THEN
C
C------------------------------------------------------
C              Reread zero'th order cluster amplitudes.
C------------------------------------------------------
C
               DTIME = SECOND()

               IOPT = 3
               CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KEND1),
     *                       WORK(KT2AM))

               DTIME = SECOND() - DTIME
               TIMIO = TIMIO + DTIME
            ENDIF
C
         ENDIF
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYRES),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NRHO2,RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of RHO1 loop over vect. 1.   ', RHO1N
            WRITE(LUPRI,1) 'Norm of RHO2 loop over vect. 1.   ', RHO2N
         ENDIF
C
C-----------------------------------------------------------------------
C        Compute R12 contributions
C-----------------------------------------------------------------------
C
         IF (CCR12) THEN
           IF (IANR12.NE.1) THEN
              CALL QUIT('NOT YET UPDATED FOR ANSAETZE 2/3')
           END IF
           IF ( DEBUG ) THEN
               RHO1N = DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
               RHO2N = DDOT(NRHO2,RHO2,1,RHO2,1)
               RHO12N = DDOT(NTR12SQ(ISYMTR),WORK(KRHOR12),1,
     &                       WORK(KRHOR12),1)
               WRITE(LUPRI,*) 'Trial vector:',IV
               WRITE(LUPRI,1) 'Norm of Rho1 -before R12:  ',RHO1N
               WRITE(LUPRI,1) 'Norm of Rho2 -before R12:  ',RHO2N
               WRITE(LUPRI,1) "Norm of Rho2' -before R12:  ",RHO12N
               WRITE(LUPRI,1) 'Norm Lambda^h',
     &          DDOT(NLAMDT,WORK(KLAMDH),1,WORK(KLAMDH),1)
               WRITE(LUPRI,1) 'Norm Lambda^p',
     &          DDOT(NLAMDT,WORK(KLAMDP),1,WORK(KLAMDP),1)
           ENDIF

C          ------------------------------
C          Calculate F' term:
C          ------------------------------
           IFILE = IV+IVEC-1
           CALL CCLHTR_FP(RHO1(1,IV),WORK(KLAMDH),WORK(KEND1),LWRK1,
     &                   ISYMTR,FC12AM,LUFC12,IFILE)

C          ------------------------------
C          Calculate G' term:
C          ------------------------------
           CALL CCLHTR_GP(CTR1(1,IV),ISYMTR,WORK(KLAMDP),1,
     &                    WORK(KRHOR12),ISYMTR,WORK(KEND1),LWRK1)

C          ------------------------------
C          Calculate CCSD contributions:
C          ------------------------------
           IF (CCSDR12) THEN
C            ------------
C            add B' term:
C            ------------
             CALL CCRHS_BPP(WORK(KRHOR12),CTR2,ISYMTR,.FALSE.,
     *                      'R12VCTDTKL',1,WORK(KEND1),LWRK1)
C            ------------
C            add A' term:
C            ------------
             !calculate M-intermediate:
             KMINT = KEND1
             KEND1 = KMINT + N3ORHF(ISYRES)
             LWRK1 = LWORK - KEND1
             IF (LWRK1 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
               CALL QUIT('Insufficient memory for M-intermediate '//
     &                   'in CC_LHTR')
             ENDIF
C
             SKIPMI = .TRUE.
             TIMEMI = SECOND()
             CALL CC_MI(WORK(KMINT),CTR2,ISYMTR,WORK(KT2AM),ISYMOP,
     *                  WORK(KEND1),LWRK1)
             TIMEMI = SECOND() - TIMEMI
C
             IF ( DEBUG ) THEN
               WRITE(LUPRI,1) 'Norm of M-INT before CCLHTR_AP: ',
     &           DDOT(N3ORHF(ISYRES),WORK(KMINT),1,WORK(KMINT),1)
             ENDIF
C
             !read V-intermediate from disk:
             KVINT = KEND1
             KEND2 = KVINT + NTR12AM(1)
             LWRK2 = LWORK - KEND2
             IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
               CALL QUIT('Insufficient memory for V-intermediate '//
     &                   'in CC_LHTR')
             ENDIF
             lunit = -1
             call gpopen(lunit,fccr12v,'old',' ','unformatted',
     &                       idum,ldum)
 9999        read(lunit) ian
             read(lunit) (work(kvint+i), i=0, ntr12am(1)-1)
             if (ian.ne.ianr12) goto 9999 
             call gpclose(lunit,'KEEP')
C
             CALL CCLHTR_AP(WORK(KRHOR12),WORK(KVINT),WORK(KMINT),
     &                      ISYMTR,WORK(KEND2),LWRK2)
C            -----------------------
C            add additional E' term:
C            -----------------------
             KOFFXI = KXMAT + NMATIJ(ISYRES)*(IV - 1)
             CALL CCLHTR_EP(WORK(KRHOR12),WORK(KVINT),WORK(KOFFXI),
     &                      ISYMTR,WORK(KEND2),LWRK2)
           END IF

C          -----------------------------------
C          Calculate E' term and add to rho12:
C          -----------------------------------
           IFILE = IV+IVEC-1
           CALL CCRHS_EP(WORK(KRHOR12),ISYMTR,.FALSE.,dummy,
     &                   WORK(KEND1),LWRK1,1,FC12AM,LUFC12,
     &                   FRHO12,LUFR12,IFILE,APROXR12,
     &                   0.5D0*KETSCL,2.0D0*BRASCL)
C
           IF ( DEBUG ) THEN
               RHO1N = DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
               RHO2N = DDOT(NRHO2,RHO2,1,RHO2,1)
               RHO12N = DDOT(NTR12SQ(ISYMTR),WORK(KRHOR12),1,
     &                       WORK(KRHOR12),1)
               WRITE(LUPRI,*) 'after CCRHS_EP:'
               WRITE(LUPRI,*) 'Trial vector:',IV
               WRITE(LUPRI,1) 'Norm of Rho1 -after R12:  ',RHO1N
               WRITE(LUPRI,1) 'Norm of Rho2 -after R12:  ',RHO2N
               WRITE(LUPRI,1) "Norm of Rho2' -after R12:  ",RHO12N
               WRITE(LUPRI,1) 'Norm Lambda^h',
     &          DDOT(NLAMDT,WORK(KLAMDH),1,WORK(KLAMDH),1)
               WRITE(LUPRI,1) 'Norm Lambda^p',
     &          DDOT(NLAMDT,WORK(KLAMDP),1,WORK(KLAMDP),1)
           ENDIF
         END IF
C
C-----------------------------------------------------------
C        Calculate adresses for CTR-dependent intermediates.
C        Skip large section for CC2.
C-----------------------------------------------------------
C
       KOFFFG = KFOCKG + N2BST(ISYMTR)*(IV - 1)
       KOFFYI = KYMAT  + NMATAB(ISYRES)*(IV - 1)
       KOFFXI = KXMAT  + NMATIJ(ISYRES)*(IV - 1)
C
       IF (.NOT. CC2) THEN
C
C------------------------------
C        Work space allocation.
C------------------------------
C
         IF (.NOT.SKIPMI) THEN
           KMINT = KENDS2 
           KEND1 = KMINT + N3ORHF(ISYRES)
           LWRK1 = LWORK - KEND1
C
           IF (LWRK1 .LT. 0) THEN
              WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
              CALL QUIT('Insufficient memory for M-intermediate '//
     &                  'in CC_LHTR')
           ENDIF
C
C--------------------------------------------------
C        Calculate M-intermediate of CTR2 and T2AM.
C--------------------------------------------------
C
           TIMEMI = SECOND()
           CALL CC_MI(WORK(KMINT),CTR2,ISYMTR,WORK(KT2AM),ISYMOP,
     *                  WORK(KEND1),LWRK1)
           TIMEMI = SECOND() - TIMEMI
C
         END IF
C
         IF ( DEBUG ) THEN
            WRITE(LUPRI,1) 'Norm of CTR2           : ',
     &        DDOT(NT2SQ(ISYMTR),CTR2,1,CTR2,1)
            WRITE(LUPRI,1) 'Norm of T2AM           : ',
     &        DDOT(NT2AMX,WORK(KT2AM),1,WORK(KT2AM),1)
            WRITE(LUPRI,1) 'Norm of M-INT          : ',
     &        DDOT(N3ORHF(ISYRES),WORK(KMINT),1,WORK(KMINT),1)
         ENDIF
C
C-----------------------------------------
C        Calculate the LT21G contribution.
C-----------------------------------------
C
         TIMEG = SECOND()
         CALL CC_21G(RHO1(1,IV),WORK(KMINT),ISYRES,WORK(KLAMDH),
     *               WORK(KEND1),LWRK1,ISYMOP,LUIN30,'CCINT3O')
         TIMEG = SECOND() - TIMEG
C
C----------------------------------------------
C        Transform the RHO2 vector to MO basis.
C----------------------------------------------
C
         AUTIME = SECOND()
         NEWGAM = .FALSE.
         CALL CC_T2MO(PHONEY,FAKE,ISYMOP,
     *                RHO2,WORK(KT2AM),DUMMY,WORK(KLAMDH),
     *                WORK(KLAMDH),1,WORK(KEND1),LWRK1,ISYRES,1)
         NEWGAM = .TRUE.
         TIM2MO = SECOND() - AUTIME
C
C-----------------------------------------------------
C        Copy the MO vector back to the result vector.
C-----------------------------------------------------
C
         CALL DCOPY(NT2AM(ISYRES),WORK(KT2AM),1,RHO2,1)
C
         IF (IPRINT .GT. 120) THEN
            CALL AROUND('Transformed vector after B-TERM')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYRES,0,1)
         ENDIF
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYRES),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of RHO1 after B-TERM:        ', RHO1N
            WRITE(LUPRI,1) 'Norm of RHO2 after B-TERM:        ', RHO2N
         ENDIF
C
C-------------------------------------------------------------------
C        Read integrals ( ma | nb ) stored as ( am | bn ) from disc.
C        Ove: CCSD_IAJB is assumed open through the complete coupled
C             cluster calculation.
C-------------------------------------------------------------------
C
         REWIND(LUIAJB)
         READ(LUIAJB) (WORK(KT2AM + J - 1), J = 1,NT2AM(ISYMOP))
C
C------------------------------------------
C        Calculate the LT22AM contribution.
C------------------------------------------
C
         AUTIME = SECOND()
         CALL CC_22AM(RHO2,WORK(KT2AM),WORK(KMINT),ISYRES,
     *                WORK(KEND1),LWRK1)
         TIMEAM = SECOND() - AUTIME
C
C------------------------------------------
C        Calculate the LT22EM contribution.
C------------------------------------------
C
         AUTIME = SECOND()
         CALL CC_22EC(RHO2,WORK(KT2AM),WORK(KOFFXI),WORK(KOFFYI),
     *                ISYRES,WORK(KEND1),LWRK1)
         CALL CC_22EE(RHO2,WORK(KT2AM),WORK(KOFFXI),WORK(KOFFYI),
     *                ISYRES,WORK(KEND1),LWRK1)
         TIM2EM = SECOND() - AUTIME
C
         IF (IPRINT .GT. 50) THEN
            CALL AROUND('Transformed vector after EM-TERM')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYRES,0,1)
         ENDIF
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYRES),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of RHO1 after EM-TERM:       ', RHO1N
            WRITE(LUPRI,1) 'Norm of RHO2 after EM-TERM:       ', RHO2N
         ENDIF
C
C---------------------------------------------
C        Regain work space from T2-amplitudes.
C---------------------------------------------
C
         KEND1 = KT2AM
         LWRK1 = LWORK - KEND1
C
C---------------------------------------------
C        Save RHO2 on disc to gain work space.
C---------------------------------------------
C
         LURHO2 = -1
         CALL GPOPEN(LURHO2,'LSRHO2','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LURHO2)
         WRITE(LURHO2) (RHO2(I), I = 1,NT2AM(ISYRES))
C
C-----------------------------------------------------------
C        Read Omega(albe,ij) written to disc by energy code.
C-----------------------------------------------------------
C
         TIMEBF = SECOND()
         LUOM = -1
         CALL GPOPEN(LUOM,'CC_BFIM','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LUOM)
         READ(LUOM) (RHO2(I),I = 1,2*NT2ORT(1) )
         CALL GPCLOSE(LUOM,'KEEP')
C
C------------------------------------------
C        Calculate the LT21BF contribution.
C------------------------------------------
C
         ISYM = 1
         ICON = 3
C
         NEWGAM = .FALSE.
         CALL CC_T2MO(RHO1(1,IV),CTR2,ISYMTR,RHO2,PHONEY,DUMMY,
     *                WORK(KLAMDP),WORK(KLAMDP),ISYM,WORK(KEND1),
     *                LWRK1,ISYMOP,ICON)
         NEWGAM = .TRUE.
         TIMEBF = SECOND() - TIMEBF
C
C-----------------------------------
C        Restore RHO2 result vector.
C-----------------------------------
C
         REWIND(LURHO2)
         READ(LURHO2) (RHO2(I), I = 1,NT2AM(ISYRES))
         CALL GPCLOSE(LURHO2,'DELETE')
C
         IF (IPRINT .GT. 50) THEN
            CALL AROUND('Transformed vectors after 21BF-TERM')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYRES,1,1)
         ENDIF
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYRES),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of RHO1 after 21BF-TERM:     ', RHO1N
            WRITE(LUPRI,1) 'Norm of RHO2 after 21BF-TERM:     ', RHO2N
         ENDIF
C
C---------------------------------------
C        Read E-intermediates from disc.
C---------------------------------------
C
         KE1INT = KEND1
         KE2INT = KE1INT + NMATAB(ISYMOP)
         KENDTE = KE2INT + NMATIJ(ISYMOP)
         LWRKTE = LWORK  - KENDTE
C
         IF (LWRKTE .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_LHTR')
         ENDIF
C
         AUTIME = SECOND()
         LUE1 = -1
         CALL GPOPEN(LUE1,'CC_E1IM','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LUE1)
         READ(LUE1) (WORK(KE1INT + J - 1), J = 1,NMATAB(ISYMOP))
         CALL GPCLOSE(LUE1,'KEEP')
C
         LUE2 = -1
         CALL GPOPEN(LUE2,'CC_E2IM','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LUE2)
         READ(LUE2) (WORK(KE2INT + J - 1), J = 1,NMATIJ(ISYMOP))
         CALL GPCLOSE(LUE2,'KEEP')
C
C--------------------------------------------------------------
C        Prepare the E-intermediates for contraction with CTR2.
C--------------------------------------------------------------
C
         CALL CC_EITR(WORK(KE1INT),WORK(KE2INT),WORK(KENDTE),LWRKTE,
     &                ISYMOP)
C
C-----------------------------------------
C        Calculate the LT22E contribution.
C-----------------------------------------
C
         CALL CCRHS_E(RHO2,CTR2,WORK(KE1INT),WORK(KE2INT),WORK(KENDTE),
     &                LWRKTE,ISYMTR,ISYMOP)
         TIM22E = SECOND() - AUTIME
C
         IF (IPRINT .GT. 50) THEN
            CALL AROUND('Transformed vector after E-TERM')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYRES,0,1)
         ENDIF
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYRES),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of RHO1 after E-TERM:        ', RHO1N
            WRITE(LUPRI,1) 'Norm of RHO2 after E-TERM:        ', RHO2N
         ENDIF
C
C------------------------------------------
C        Read Gamma-intermediate from disc.
C------------------------------------------
C
         KGAMMI = KEND1
         KENDGI = KGAMMI + NGAMMA(ISYMOP)
         LWRKGI = LWORK  - KENDGI
C
         IF (LWRKGI .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_LHTR')
         ENDIF
C
         AUTIME = SECOND()
         LUGI = -1
         CALL GPOPEN(LUGI,'CC_GAMIM','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LUGI)
         READ(LUGI) (WORK(KGAMMI + J - 1), J = 1,NGAMMA(ISYMOP))
         CALL GPCLOSE(LUGI,'KEEP')
C
C-----------------------------------------
C        Calculate the LT22A contribution.
C-----------------------------------------
C
         ISYG = ISYMOP
         ISYV = ISYMTR
         IOPT = 2
C
         CALL CCRHS_A(RHO2,CTR2,WORK(KGAMMI),WORK(KENDGI),LWRKGI,
     &                ISYG,ISYV,IOPT)
         TIM22A = SECOND() - AUTIME
C
C-----------------------------------------
C        Calculate the LT22D contribution.
C-----------------------------------------
C
         AUTIME = SECOND()
         IOPT = 1
         ISYDIM = 1
         CALL CC_22D(RHO2,CTR2,ISYMTR,WORK(KLAMDH),WORK(KEND1),
     *               LWRK1,ISYDIM,LUDIM,DFIL,IV,IOPT)
         TIM22D = SECOND() - AUTIME
C
C-----------------------------------------
C        Calculate the LT22C contribution.
C-----------------------------------------
C
         AUTIME = SECOND()
C
         CALL CCRHS_T2BT(CTR2,WORK(KEND1),LWRK1,ISYMTR)
         CALL DSCAL(NT2SQ(ISYMTR),THREE,CTR2,1)
         CALL CCSD_T2TP(CTR2,WORK(KEND1),LWRK1,ISYMTR)
C
         IOPT = 1
         ISYCIM = 1
         CALL CC_22C(RHO2,CTR2,ISYMTR,WORK(KLAMDH),WORK(KEND1),
     *               LWRK1,ISYCIM,LUCIM,CFIL,IV,IOPT)
C
         TIM22C = SECOND() - AUTIME
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYRES),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of RHO1 after CC_22C:        ', RHO1N
            WRITE(LUPRI,1) 'Norm of RHO2 after CC_22C:        ', RHO2N
         ENDIF
C
       ENDIF
C
C--------------------------------------
C     Calculate the LT11A contribution.
C--------------------------------------
C
      AUTIME = SECOND()
      CALL CC_11A(RHO1(1,IV),WORK(KOFFFG),ISYRES,WORK(KLAMDH),
     *            WORK(KLAMDP),WORK(KEND1),LWRK1)
C
      IF ( DEBUG ) THEN
         RHO1N = DDOT(NT1AM(ISYRES),RHO1(1,IV),1,RHO1(1,IV),1)
         RHO2N = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
         WRITE(LUPRI,1) 'Norm of RHO1 after CC_11A:        ', RHO1N
         WRITE(LUPRI,1) 'Norm of RHO2 after CC_11A:        ', RHO2N
      ENDIF
C
C------------------------------------------
C     Ove: Read in AO Fock.
C     Transform AO Fock matrix to MO basis.
C------------------------------------------
C
      LFOCK = -1
      CALL GPOPEN(LFOCK,'CC_FCKH','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LFOCK)
      READ(LFOCK) (WORK(KFOCK + I - 1) , I = 1,N2BST(ISYMOP))
      CALL GPCLOSE(LFOCK,'KEEP')
C
      IHELP = 1
C
      CALL CC_FCKMO(WORK(KFOCK),WORK(KLAMDP),WORK(KLAMDH),
     *                 WORK(KEND1),LWRK1,IHELP,IHELP,IHELP)
C
C------------------------------------------------
C     Calculate the LT11B and LT11C contribution.
C------------------------------------------------
C
      CALL CCLT_11BC(RHO1(1,IV),CTR1(1,IV),WORK(KFOCK),WORK(KEND1),
     *               LWRK1)
      TIM11B = SECOND() - AUTIME
C
C--------------------------------------
C     Calculate the LT12A contribution.
C--------------------------------------
C
      AUTIME = SECOND()
      CALL CC_L1FCK(RHO2,CTR1(1,IV),WORK(KFOCK),ISYMTR,ISYMOP,
     *              WORK(KEND1),LWRK1)
      TIM12A = SECOND() - AUTIME
C
      IF ( DEBUG ) THEN
         RHO1N = DDOT(NT1AM(ISYRES),RHO1(1,IV),1,RHO1(1,IV),1)
         RHO2N = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
         WRITE(LUPRI,1) 'Norm of RHO1 after _11BC and _12A:', RHO1N
         WRITE(LUPRI,1) 'Norm of RHO2 after _11BC and _12A:', RHO2N
      ENDIF
C
C
C-----------------------------------------------
C     Calculate the LT12C & LT21B contributions.
C-----------------------------------------------
C
      TIMEB = SECOND()
      IOPT  = 2
      CALL CC_21B12C(RHO1(1,IV),RHO2,CTR1(1,IV),ISYMTR,WORK(KLAMDH),
     *               ISYMOP,WORK(KOFFXI),ISYRES,ISYMOP,
     *               WORK(KEND1),LWRK1,LUIN30,'CCINT3O',IOPT)
      TIMEB = SECOND() - TIMEB
C
      IF ( DEBUG ) THEN
         RHO1N = DDOT(NT1AM(ISYRES),RHO1(1,IV),1,RHO1(1,IV),1)
         RHO2N = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
         WRITE(LUPRI,1) 'Norm of RHO1 after CC_12C         ', RHO1N
         WRITE(LUPRI,1) 'Norm of RHO2 after CC_12C         ', RHO2N
      ENDIF
C
C----------------------------------------------------------------
C     Calculate the diagonal 2.2-contribution in CC2-calculation.
C----------------------------------------------------------------
C
      IF (CC2) THEN
C
         TIMEC2 = SECOND()
         IF (.NOT. NONHF) THEN
           ISIDE =  -1
           CALL CC2_FCK(RHO2,CTR2,WORK(KEND1),LWRK1,ISYMTR,
     *                  WORK(KLAMDP),WORK(KLAMDH),ISIDE)
         ELSE
           KEMAT1 = KEND1
           KEMAT2 = KEMAT1 + NEMAT1(1)
           KEND2  = KEMAT2 + NMATIJ(1)
           LWRK2  = LWORK  - KEND2
           IF (LWRK2 .LT. 0) THEN
              CALL QUIT('Insufficient work space in CC_LHTR')
           ENDIF

           ETRAN  = .FALSE.
           FCKCON = .TRUE.
           CALL CCRHS_EFCK(WORK(KEMAT1),WORK(KEMAT2),WORK(KLAMDH),
     *                   WORK(KFCKHF),WORK(KEND2),LWRK2,FCKCON,ETRAN,1)
           CALL CC_EITR(WORK(KEMAT1),WORK(KEMAT2),WORK(KEND2),LWRK2,1)
           CALL CCRHS_E(RHO2,CTR2,WORK(KEMAT1),WORK(KEMAT2),WORK(KEND2),
     *                  LWRK2,ISYMTR,1)
         END IF
         TIMEC2 = SECOND() -TIMEC2
C
      ENDIF
C
C----------------------------------------
C     Calculate the LT21EFM contribution.
C----------------------------------------
C
      IF (.NOT. CC2) THEN
C
         KOFFYI = KYMAT  + NMATAB(ISYRES)*(IV - 1)
         KOFFXI = KXMAT  + NMATIJ(ISYRES)*(IV - 1)
         TIMEEM = SECOND()
         CALL CC_21EFM(RHO1(1,IV),WORK(KFOCK),ISYMOP,WORK(KOFFXI),
     *                 WORK(KOFFYI),ISYRES,WORK(KEND1),LWRK1)
         TIMEEM = SECOND() - TIMEEM
C
      ENDIF
C
C-------------------------------------------
C     SLV98,OC
C     Calculate T^{gB} solvent contribution.
C     Unsquared T2 is required!
C-------------------------------------------
C
      IF (CCSLV.AND.LSTBTR.AND. (.NOT. CCMM)) THEN
         CALL CC_RVEC(LUFC2,FC2AM,NT2AM(ISYMTR),NT2AM(ISYMTR),
     *                IVEC+IV-1,CTR2)
         LR = 'L'
         DTIME = SECOND()
         CALL CCSL_LTRB(RHO1(1,IV),RHO2,CTR1(1,IV),CTR2,ISYMTR,LR,
     *                  WORK(KEND1),LWRK1)
c        TIMMLT = TIMMLT + SECOND() - DTIME
         IF (IPRINT .GT. 5) THEN
            WRITE(LUPRI,*)'Time used in CCSL_LTRB (CC_LHTR) 1',
     *                     SECOND() - DTIME
         END IF
      ENDIF
C
C
C-------------------------------------------
C     CCMM02,JA+AO
C     Calculate T^{gB} solvent contribution.
C     Unsquared T2 is required!
C-------------------------------------------
C
      IF (CCMM.AND.LSTBTR) THEN
         CALL CC_RVEC(LUFC2,FC2AM,NT2AM(ISYMTR),NT2AM(ISYMTR),
     *                IVEC+IV-1,CTR2)
         LR = 'L'
         DTIME = SECOND()
         IF (.NOT. NYQMMM) THEN
            CALL CCMM_LTRB(RHO1(1,IV),RHO2,CTR1(1,IV),CTR2,
     *                  ISYMTR,LR,WORK(KEND1),LWRK1)
         ELSE IF (NYQMMM) THEN
             CALL CCMM_TRANSFORMER(RHO1(1,IV),RHO2,CTR1(1,IV),CTR2,
     *                            MODEL,ISYMTR,LR,WORK(KEND1),LWRK1)
         END IF
c        TIMMLT = TIMMLT + SECOND() - DTIME
         IF (IPRINT .GT. 5) THEN
            WRITE(LUPRI,*)'Time used in CCMM_LTRB (CC_LHTR) 1',
     *                     SECOND() - DTIME
         END IF
      ENDIF
C
      IF (USE_PELIB().AND.LSTBTR.AND.(.NOT.HFFLD)) THEN
         CALL CC_RVEC(LUFC2,FC2AM,NT2AM(ISYMTR),NT2AM(ISYMTR),
     &                IVEC+IV-1,CTR2)
         LR = 'L'
         DTIME = SECOND()
         CALL PELIB_IFC_TRANSFORMER(RHO1(1,IV),RHO2,CTR1(1,IV),CTR2,
     &                         MODEL,ISYMTR,LR,WORK(KEND1),LWRK1)
         IF (IPRINT .GT. 5) THEN
            WRITE(LUPRI,*)'Time used in PECC TRANSFORMER (CC_LHTR) 1',
     &                     SECOND() - DTIME
         END IF
      ENDIF
C
C------------------------------------------------------------
C     Calculate the triples part of the left hand side
C     eigenvalue equation.
C     Start by calculating different integrals.
C     CC3_T3_LHTR calculate the terms from T3
C     whereas CC3_L3_LHTR calculate the terms from L3
C------------------------------------------------------------
C
      IF (CC3) THEN
C
         KOOOO = KEND1
         KEND1 = KOOOO + N3ORHF(1)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Out of memory in CC_LHTR (CC3)')
         ENDIF
C
         CALL DZERO(WORK(KOOOO),N3ORHF(1))
C
C        Since g^0 use lamH^0 twice in call
C
         IF (LVVVV) THEN
            IOPT = 3
         ELSE
            IOPT = 1
         END IF
         CALL CC3_INTSTORE(LUIN30,'CCINT3O',WORK(KOOOO),1,
     *                     WORK(KLAMDH),ISYMOP,WORK(KLAMDH),ISYMOP,
     *                     LU3V,FN3V,LU4V,FN4V,ISYMOP,
     *                     WORK(KEND1),LWRK1,IOPT)
C
         CALL CC3_SORT4O(WORK(KOOOO),1,WORK(KEND1),LWRK1)
C
         KT1AM = KEND1
         KT2AM = KT1AM + NT1AM(ISYMOP)
         KEND2 = KT2AM + NT2SQ(ISYMOP)
         LWRK2 = LWORK - KEND2
C
         IF (LWRK2 .LT. NT2AM(ISYMOP)) THEN
            CALL QUIT('Out of workspace in CC_LHTR (CC3 part)')
         ENDIF
C
         IOPT = 3
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),WORK(KEND2))
C
         CALL CC_T2SQ(WORK(KEND2),WORK(KT2AM),1)
C
C
         CALL CC_RVEC(LUFC2,FC2AM,NT2AM(ISYMTR),NT2AM(ISYMTR),
     *                IVEC+IV-1,WORK(KEND2))
         DTIME = SECOND() - DTIME
         TIMIO = TIMIO + DTIME
C
         CALL CC_T2SQ(WORK(KEND2),CTR2,ISYMTR)
C
C-----------------------
C        Calculate
C-----------------------
C
         IF (.NOT.NODDY_LHTR) THEN

          CALL CC3_T3_LHTR(0.0D0,RHO1(1,IV),WORK(KT1AM),ISYMOP,
     *                    WORK(KT2AM),ISYMOP,CTR2,ISYMTR,
     *                    WORK(KLAMDP),WORK(KLAMDH),
     *                    WORK(KEND2),LWRK2,LU3SRT,FN3SRT,
     *                    LUCKJD,FNCKJD)

          IF (IPRINT .GT. 55) THEN
            RHO1N = DDOT(NT1AM(ISYRES),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of RHO1 after CC3_T3_LHTR :', RHO1N
            WRITE(LUPRI,1) 'Norm of RHO2 after CC3_T3_LHTR :', RHO2N
          ENDIF

         END IF
C
C----------------------------------------------------------------
C        L2 and T2 are destroyed by CC3_T3_LHTR.
C        Readin again.
C----------------------------------------------------------------
C
         CALL CC_RVEC(LUFC2,FC2AM,NT2AM(ISYMTR),NT2AM(ISYMTR),
     *                IVEC+IV-1,WORK(KEND2))
         DTIME = SECOND() - DTIME
         TIMIO = TIMIO + DTIME
C
         CALL CC_T2SQ(WORK(KEND2),CTR2,ISYMTR)
C
         IOPT = 2
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KEND2))
C
         CALL CC_T2SQ(WORK(KEND2),WORK(KT2AM),1)
C
         IF (.NOT.NODDY_LHTR) THEN
          CALL CC3_L3_LHTR(ECURR,CTR1(1,IV),ISYMTR,CTR2,ISYMTR,
     *                    WORK(KT2AM),ISYMOP,
     *                    RHO1(1,IV),RHO2,ISYRES,
     *                    WORK(KOOOO),WORK(KOVVO),WORK(KOOVV),
     *                    WORK(KLAMDP),WORK(KLAMDH),
     *                    WORK(KEND2),LWRK2,
     *                    LUCKJD,FNCKJD,LUDKBC,FNDKBC,LUTOC,FNTOC,
     *                    LU3VI,FN3VI,LU4V,FN4V,LUDKBC3,FNDKBC3,
     *                    LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X)
         END IF
C
         IF (NODDY_LHTR) THEN
             IF (DEBUG) CALL AROUND('Before CC_LHTR_NODDY')
             CALL CC_LHTR_NODDY(ECURR,RHO1(1,IV),RHO2,WORK(KT1AM),
     *                          WORK(KT2AM),CTR1(1,IV),CTR2,
     *                          WORK(KLAMDP),WORK(KLAMDH),
     *                          WORK(KEND2),LWRK2,0)
         ENDIF
C
         IF (IPRINT .GT. 55) THEN
            RHO1N = DDOT(NT1AM(ISYRES),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of RHO1 after CC3_L3_LHTR :', RHO1N
            WRITE(LUPRI,1) 'Norm of RHO2 after CC3_L3_LHTR :', RHO2N
         ENDIF

      ENDIF
C--------------------------------------------------------------------
C     Print out result vectors - zero out single vectors in CCD-calc.
C--------------------------------------------------------------------
C
      IF (CCD) CALL DZERO(RHO1(1,IV),NT1AM(ISYRES))
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND('Transformed vectors coming out of CC_LHTR')
         CALL CC_PRP(RHO1(1,IV),RHO2,ISYRES,1,1)
      ENDIF
C
      IF ( DEBUG ) THEN
         RHO1N = DDOT(NT1AM(ISYRES),RHO1(1,IV),1,RHO1(1,IV),1)
         RHO2N = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
         WRITE(LUPRI,1) 'Norm of RHO1 coming out of CC_LHTR:', RHO1N
         WRITE(LUPRI,1) 'Norm of RHO2 coming out of CC_LHTR:', RHO2N
      ENDIF
C
C----------------------------------
C     Write out transformed vector.
C     Write out in all cases.
C----------------------------------

         CALL CC_WVEC(LUFR1,FRHO1,NT1AM(ISYMTR),NT1AM(ISYMTR),
     *                IV + IVEC -1,RHO1(1,IV))
C
         DTIME   = SECOND()
         CALL CC_WVEC(LUFR2,FRHO2,NRHO2,NT2AM(ISYMTR),
     *                   IV + ITR -1,RHO2)
         DTIME   = SECOND() - DTIME
         TIMIO   = TIMIO    + DTIME
C
  125 CONTINUE
C
C------------------------------
C     Close intermediate files.
C------------------------------
C
      CALL WCLOSE2(LUCIM,CFIL,'KEEP')
      CALL WCLOSE2(LUDIM,DFIL,'KEEP')
C
C----------------------------------------
C     Close intermediate files for P & Q.
C----------------------------------------
C
      IF (.TRUE.) THEN
         CALL WCLOSE2(LUPQ,PQFIL,'KEEP')
      END IF
C
C----------------------------------------
C     Close files for CC3.
C----------------------------------------
C
      IF (CC3) THEN
C
         IF (LVVVV) THEN
            CALL WCLOSE2(LU4V,FN4V,'DELETE')
            CALL WCLOSE2(LU3V,FN3V,'KEEP')
         END IF
C
         CALL WCLOSE2(LU3SRT,FN3SRT,'KEEP')
         CALL WCLOSE2(LUCKJD,FNCKJD,'KEEP')
         CALL WCLOSE2(LUDKBC,FNDKBC,'KEEP')
         CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
         CALL WCLOSE2(LU3VI,FN3VI,'KEEP')
         CALL WCLOSE2(LUDKBC3,FNDKBC3,'KEEP')
         CALL WCLOSE2(LU3FOPX,FN3FOPX,'KEEP')
         CALL WCLOSE2(LU3FOP2X,FN3FOP2X,'KEEP')
C
      ENDIF
C
C-------------------------------
C     Delete intermediate files.
C-------------------------------
C
      CALL WCLOSE2(LUIN30,'CCINT3O','DELETE')
C
      TIMTOT = SECOND() - TIMTOT
C
C-------------------------------
C     Write out program timings.
C-------------------------------
C
      IF (IPRINT .GT. 3) THEN
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Time used in CC_LHTR :',TIMTOT
         WRITE(LUPRI,*) ' '
      ENDIF
      IF (IPRINT .GT. 9) THEN
         WRITE(LUPRI,*) 'Time used as follows:'
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Time used in HERMIT1:',TIMHE1
         WRITE(LUPRI,*) 'Time used in HERMIT2:',TIMHE2
         WRITE(LUPRI,*) 'Time used in CCRDAO :',TIRDAO
         WRITE(LUPRI,*) 'Time used in Vect.IO:',TIMIO
         WRITE(LUPRI,*) 'Time used in CCLT_YI:',TIMEYI
         WRITE(LUPRI,*) 'Time used in CCLT_XI:',TIMEXI
         WRITE(LUPRI,*) 'Time used in CCLT_MI:',TIMEMI
         WRITE(LUPRI,*) 'Time used in CCLT_TI:',TIMETI
         WRITE(LUPRI,*) 'Time used in CCLT_Z1:',TIMEZ1
         WRITE(LUPRI,*) 'Time used in CCLT_Z2:',TIMEZ2
         WRITE(LUPRI,*) 'Time used in CCLT_C :',TIMEC
         WRITE(LUPRI,*) 'Time used in CCLT_A :',TIMEA
         WRITE(LUPRI,*) 'Time used in CCLT_I :',TIMEI
         WRITE(LUPRI,*) 'Time used in CCINT3O:',TIME3O
         WRITE(LUPRI,*) 'Time used in CCTRBT :',TIME1O
         WRITE(LUPRI,*) 'Time used in AOFOCK :',TIMFCK
         WRITE(LUPRI,*) 'Time used in CCT2TP :',TIMETP
         WRITE(LUPRI,*) 'Time used in CCT2AO :',TIM2AO
         WRITE(LUPRI,*) 'Time used in CCT2MO :',TIM2MO
         WRITE(LUPRI,*) 'Time used in CCLT_H :',TIMEH
         WRITE(LUPRI,*) 'Time used in LT_21BF:',TIMEBF
         WRITE(LUPRI,*) 'Time used in LT_22BF:',TIM2BF
         WRITE(LUPRI,*) 'Time used in CCLT_EM:',TIMEEM
         WRITE(LUPRI,*) 'Time used in CCLT_G :',TIMEG
         WRITE(LUPRI,*) 'Time used in 21B&12C:',TIMEB
         WRITE(LUPRI,*) 'Time used in CC2_FCK:',TIMEC2
         WRITE(LUPRI,*) 'Time used in LT_12B :',TIM12B
         WRITE(LUPRI,*) 'Time used in LT_12A :',TIM12A
         WRITE(LUPRI,*) 'Time used in 11BLOCK:',TIM11B
         WRITE(LUPRI,*) 'Time used in LT_22AM:',TIMEAM
         WRITE(LUPRI,*) 'Time used in LT_22EM:',TIM2EM
         WRITE(LUPRI,*) 'Time used in LT_22E :',TIM22E
         WRITE(LUPRI,*) 'Time used in LT_22A :',TIM22A
         WRITE(LUPRI,*) 'Time used in LT_22C :',TIM22C
         WRITE(LUPRI,*) 'Time used in LT_22D :',TIM22D
      ENDIF
C
   1  FORMAT(1x,A35,1X,E20.10)
C
C
      CALL QEXIT('CC_LHTR')
      RETURN
      END
