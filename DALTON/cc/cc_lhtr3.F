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
C  /* Deck cc_lhtr */
      SUBROUTINE CC_LHTR3(ECURR,
     *                   FRHO1,LUFR1,FRHO2,LUFR2,
     *                   FC1AM,LUFC1,FC2AM,LUFC2,
     *                   RHO1,RHO2,CTR1,CTR2,WORK,LWORK,
     *                   NSIMTR,IVEC,ITR,LRHO1)
C
C     Written by Rasmus Faber, summer 2017
C
C     Based on CC_LHTR by Asger Halkier & Henrik Koch.
C
C     Purpose:
C
C     Calculate left hand side transformation of a triplet trial vector
C     in an AO-integral direct fashion.
C
C     (IVEC is first number for first vector on files, FRHO1,FC1AM,FC2AM,
C           FRHO12,FC12AM;
C      ITR is first vector on FRHO2)
C
C     Current models are: CCSD..?
C
C
C
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
      CHARACTER(LEN=*), PARAMETER :: myname = 'CC_LHTR3'
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
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0,
     &           FOUR = 4.0D0)
      PARAMETER (FOURTH = 0.25D0, THREE = 3.0D0)
      PARAMETER (ONEM = -1.0D0, TWOM = -TWO)
      DIMENSION INDEXA(MXCORB_CC)
      DIMENSION RHO1(LRHO1,NSIMTR), CTR1(LRHO1,NSIMTR)
      DIMENSION RHO2(*), CTR2(*)
      DIMENSION WORK(LWORK)
      DIMENSION IADRPQ(MXCORB_CC,MAXSIM)
C
      LOGICAL ETRAN,FCKCON,CCSDR12,SKIPMI
C
      CALL QENTER(myname)
C
      THIRD = ONE/THREE
C
C
      !set CCSDR12
      CCSDR12 = .FALSE.
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
         CALL QUIT('Insufficient space in '//myname)
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
      NRHO2 = MAX(NT2SQ(ISYRES),NT2ORT(ISYRES)+NT2ORT3(ISYRES))
      IF (CC2) NRHO2 = NT2SQ(ISYRES)
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
         IOPT = -1
         CALL CC_C1T2C(WORK(KOFFCT),CTR1(1,IV),WORK(KT2AM),
     *                 WORK(KLAMDP),WORK(KLAMDH),WORK(KLAMDP),
     *                 WORK(KLAMDH),WORK(KEND1),LWRK1,ISYMTR,
     *                 ISYMOP,IOPT)
C
C--------------------------------------
C        Read CTR2+ from disc into RHO2.
C--------------------------------------
C
         DTIME = SECOND()
         CALL CC_RVEC3(LUFC2,FC2AM,2*NT2AM(ISYMTR),NT2AM(ISYMTR),
     *                 IVEC+IV-1,0,RHO2)
         DTIME = SECOND() - DTIME
         TIMIO = TIMIO + DTIME
C
C-----------------------
C        Square up CTR2.
C-----------------------
C
         CALL CCRHS3_R2IJ(RHO2,WORK(KEND1),LWRK1,ISYMTR)
         CALL CC_T2SQ(RHO2,CTR2,ISYMTR)
C
C--------------------------------------
C        Read CTR2- from disc into RHO2.
C--------------------------------------
C
         DTIME = SECOND()
         CALL CC_RVEC3(LUFC2,FC2AM,2*NT2AM(ISYMTR),NT2AM(ISYMTR),
     *                 IVEC+IV-1,NT2AM(ISYMTR),RHO2)
         DTIME = SECOND() - DTIME
         TIMIO = TIMIO + DTIME
C
         CALL CC_T2SQ3A(RHO2,CTR2,ISYMTR,ONEM)
C
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
            KOFFYD = KYDEN + N2BST(ISYRES)*(IV - 1) ! Not needed
C
            TIMEYI = SECOND()
            CALL CC_YI(WORK(KOFFYI),CTR2,ISYMTR,WORK(KT2AM),ISYMOP,
     *                 WORK(KEND1),LWRK1)
            CALL CC_YD(WORK(KOFFYD),WORK(KOFFYI),ISYRES,WORK(KLAMDH),
     *                 WORK(KLAMDP),ISYMOP,WORK(KEND1),LWRK1)
            TIMEYI = SECOND() - TIMEYI

            CALL DAXPY(N2BST(ISYRES),ONE,
     &                 WORK(KOFFYD),1,WORK(KOFFCT),1)
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
            CALL CC_XD3(WORK(KOFFCT),WORK(KOFFXI),ISYRES,WORK(KLAMDH),
     &                  WORK(KLAMDP),1,WORK(KEND1),LWRK1)
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
C
C   Calculate the triplet P-matrix
C
            CALL CC_PQI3(CTR2,ISYMTR,WORK(KT2AM),ISYMOP,PQFIL,LUPQ,
     *                  IADRPQ,IADRSTRT,IV,WORK(KLAMDH),0,
     &                  WORK(KEND1),LWRK1)
C
C    We need to rearrange CTR2 in order to calculate triplet Q matrix:
C    First get -(+)L - (-)L from current (+)L - (-)L by transposing,
C    then scaling

            CALL CC_T2SQTRANSP(CTR2,ISYMTR)
C
            CALL DSCAL(NT2SQ(ISYMTR),ONEM,CTR2,1)
C
C    RHO2 still contains (-)L, reorder it so (-)L(ai<bj) to -(-)L(aj<bi)
C
            CALL CCSD_TCMEPK3(RHO2,ISYMTR)
C
C    Create L'(ai,bj) = -(+)L(ai,bj)-(-)L(ai,bj) + 2(-)L(aj,bi)
            CALL CC_T2SQ3A(RHO2,CTR2,ISYMTR,TWO)
C
            CALL CC_PQI3(CTR2,ISYMTR,WORK(KT2AM),ISYMOP,PQFIL,LUPQ,
     *                  IADRPQ,IADRSTRT,IV,WORK(KLAMDH),1,
     &                  WORK(KEND1),LWRK1)
C
C    Create L'(ai,bj) = -(+)L(ai,bj)-(-)L(ai,bj) - 2(-)L(aj,bi)
C
            CALL CC_T2SQ3A(RHO2,CTR2,ISYMTR,-FOUR)
C
            CALL CC_PQI3(CTR2,ISYMTR,WORK(KT2AM),ISYMOP,PQFIL,LUPQ,
     *                  IADRPQ,IADRSTRT,IV,WORK(KLAMDH),2,
     &                  WORK(KEND1),LWRK1)

C
            IF (MINSCR) THEN
C
C    We better clean up the mess....
C    Restore (+)L + (-)L on CTR2
               CALL CC_T2SQ3A(RHO2,CTR2,ISYMTR,TWO)
               CALL DSCAL(NT2SQ(ISYMTR),ONEM,CTR2,1)

            END IF
C            DTIME = SECOND()
C            CALL CC_RVEC(LUFC2,FC2AM,NT2AM(ISYMTR),NT2AM(ISYMTR),
C     *                   IVEC+IV-1,RHO2)
C            DTIME = SECOND() - DTIME
C            TIMIO = TIMIO + DTIME
C
C            CALL CC_T2SQ(RHO2,CTR2,ISYMTR)
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
C       In the integral loop, we need the + combinations of (+) and (-)
C       We currently have the - combination, the (-) part is in RHO2
C
C         CALL CC_T2SQ3A(RHO2,CTR2,ISYMTR,TWO)
         CALL DZERO(RHO2,NRHO2)
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

                  CALL QUIT(
     &               'Triplet LHTR without MINSCR not implemented'
     &               )
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
C                  IF (NOAUXB) THEN
c                     IDUM = 1
C                     CALL IJKAUX(IDEL,IDUM,IDUM,IDUM)
C                  END IF
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
C              Keep Exchange part only
C-------------------------------------------
C
               AUTIME = SECOND()
               ISYDEN = ISYMTR
               CALL CC_AOFOCK3(WORK(KXINT),WORK(KOFFCT),WORK(KOFFFG),
     *                        WORK(KEND2),LWRK2,IDEL,ISYMD,ISYDEN)
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
                  KSCRTI = K3OINT + NMAIJK(ISYMTI) ! is KSCRTI ever used?
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
C--------------------------------------------------------------
C              Calculate intermediates needed for the 21-block.
C--------------------------------------------------------------
C
               IF (.NOT. CC2) THEN
C
                  IADR = IADRPQ(IDEL,IV)
                  CALL GETWA2(LUPQ,PQFIL,WORK(KSCRZI),
     *                        IADR,NT2BCD(ISY21I))
C
                  IADR = IADRPQ(IDEL,IV) + NT2BCD(ISY21I)
                  CALL GETWA2(LUPQ,PQFIL,WORK(KSCRVI),
     *                        IADR,NT2BCD(ISY21I))
C
                  IADR = IADRPQ(IDEL,IV) + 2*NT2BCD(ISY21I)
                  CALL GETWA2(LUPQ,PQFIL,WORK(KSCRWI),
     *                        IADR,NT2BCD(ISY21I))
C
                  IF ( DEBUG ) THEN
                   PN=DDOT(NT2BCD(ISY21I),WORK(KSCRZI),1,WORK(KSCRZI),1)
                   QN=DDOT(NT2BCD(ISY21I),WORK(KSCRVI),1,WORK(KSCRVI),1)
                   WRITE(LUPRI,*) 'IV,IDEL,IADR,P,Q',IV,IDEL,IADR,PN,QN
                  ENDIF

C
C--------------------------------------------------
C                 Calculate the LT21I contribution.
C--------------------------------------------------
C

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

C
C--------------------------------------------------
C                 Calculate the LT21H contribution.
C-------------------------------------------------
C
                  CALL DZERO(WORK(KSCRVI),NT2BCD(ISY21I))
                  AUTIME = SECOND()
                  CALL CC_21H(RHO1(1,IV),ISYRES,WORK(KSCRWI),
     *                        WORK(KSCRVI),WORK(KSCRZI),ISYLRD,
     *                        WORK(K3OINT),ISYMOP,WORK(KEND3),
     *                        LWRK3,ISYMD)
                  AUTIME = SECOND() - AUTIME
                  TIMEH = TIMEH + AUTIME

                  IF ( DEBUG ) THEN
                    RHO1N=DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
                    WRITE(LUPRI,*) 'Norm of RHO1 after CC_21H:',IV,RHO1N
                  ENDIF

               ENDIF! (.NOT. CC)
C
C------------------------------------------
C              Work space allocation no. 4.
C------------------------------------------
C
               ISSCRM = MULD2H(ISYMD,ISYMTR)
C
               KSCRM = KEND3 ! KSCR{T,V,W,Z} is not uused beyond this!
                             ! Reset to KSCRTI?
               KSCRM2 = KSCRM  + NT2BCD(ISSCRM)
               KEND4 = KSCRM2  + NT2BCD(ISSCRM)
               LWRK4 = LWORK  - KEND4
C
               IF (LWRK4 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND4,'Available : ',LWORK
                  CALL QUIT('Insufficient space in'//myname)
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

               CALL CC_T2AO3(CTR2,WORK(KLAMDP),ISYMLP,WORK(KSCRM2),
     *                        WORK(KEND4),LWRK4,IDEL,ISYMD,
     *                        ISYMTR)
C
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
                  CALL CC_BF3(WORK(KXINT),RHO2,WORK(KLAMDP),1,
     *                        WORK(KOFFAO),ISYMTR,WORK(KLAMDP),1,
     *                        WORK(KSCRM),ISYMM,WORK(KSCRM2),ISYMM,
     *                        WORK(KEND4),LWRK4,IDEL,ISYMD,IOPT)
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
                  DTIME = SECOND()
                  CALL CC_WVEC(LUFR2,FRHO2,NRHO2,NRHO2,IV + ITR -1,
     *                         RHO2)
                  DTIME = SECOND() - DTIME
                  TIMIO = TIMIO    + DTIME
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
C-----------------------------
C     Loop over trial vectors.
C-----------------------------
C
      DO 125 IV = 1,NSIMTR
C
C
         IF (.NOT. MINSCR) THEN
C
C           Read in CTR2+ and CTR. and
C           construct CTR+ - CTR- intermediate
C           T2 vectors are overwritten further down
C           and needs to be read in once more if
C           not first iteration of the loop.
C           Alternatively we could calculate and
C           write out the M intermediate before the
C           loop over integrals.
         ELSE
C
C---------------------------------------
C        We have in memory CTR2+ + CTR2-
C        in the following we need
C        CTR2+ - CTR-
C---------------------------------------
C
            CALL CC_T2SQTRANSP(CTR2,ISYMTR)

         ENDIF
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYRES),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NRHO2,RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of RHO1 loop over vect. 1.   ', RHO1N
            WRITE(LUPRI,1) 'Norm of RHO2 loop over vect. 1.   ', RHO2N
         ENDIF
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
           KMINT = KT2AM + MAX(NT2SQ(ISYRES),NT2AM(1))
           KEND1 = KMINT + N3ORHF(ISYRES)
           LWRK1 = LWORK - KEND1
C
           IF (LWRK1 .LT. 0) THEN
              WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
              CALL QUIT('Insufficient memory for M-intermediate '//
     &                  'in CC_LHTR')
           ENDIF
C
C-----------------------------------------------------------
C        Calculate transposed M-intermediate of CTR2 and T2AM.
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
         CALL CC_T2SQTRANSP(CTR2,ISYMTR)
C
C----------------------------------------------
C        Transform the RHO2 vector to MO basis.
C----------------------------------------------
C
         AUTIME = SECOND()
         CALL CC_T2MOT(PHONEY,FAKE,ISYMOP,
     *                 RHO2,WORK(KT2AM),DUMMY,WORK(KLAMDH),
     *                 WORK(KLAMDH),1,WORK(KEND1),LWRK1,ISYRES,1)
         TIM2MO = SECOND() - AUTIME
C
C-----------------------------------------------------
C        Copy the MO vector back to the result vector.
C-----------------------------------------------------
C
         CALL DCOPY(NT2SQ(ISYRES),WORK(KT2AM),1,RHO2,1)
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
C        This seems to be called the E contribution
C        in the article..?
C------------------------------------------
C
         AUTIME = SECOND()
         CALL CC_22AM3(RHO2,WORK(KT2AM),WORK(KMINT),ISYRES,
     *                WORK(KEND1),LWRK1)
         TIMEAM = SECOND() - AUTIME
C
C------------------------------------------
C        Read Gamma-intermediate from disc.
C------------------------------------------
C
         KGAMMI = KENDS2
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
C
         CALL CCRHS_ASQ(RHO2,CTR2,WORK(KGAMMI),WORK(KENDGI),LWRKGI,
     &                  ISYG,ISYV,IOPT)
         TIM22A = SECOND() - AUTIME
C
C-----------------------------------------------------------------------
C        The above terms carries a factor of 1/2 on the (+) contribution
C        relative to the (-) contribution.
C        Take care of that by scaling the symmetric part py 1/2
C-----------------------------------------------------------------------
C
         CALL CC_T2SQSYMSCAL(RHO2,ISYMTR,HALF)
C
C------------------------------------------
C        Calculate the LT22EM contribution.
C        This seems to be called the F term
C        the article?
C------------------------------------------
C
         AUTIME = SECOND()
         CALL CC_22E3(RHO2,WORK(KT2AM),WORK(KOFFXI),WORK(KOFFYI),
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
C        Calculate the LT11B and LT11C contribution.
C------------------------------------------------
C
         CALL CCLT_11BC(RHO1(1,IV),CTR1(1,IV),WORK(KFOCK),WORK(KEND1),
     *                  LWRK1)
         TIM11B = SECOND() - AUTIME
C
C--------------------------------------
C        Calculate the LT12A (H1) contribution.
C--------------------------------------
C
         AUTIME = SECOND()
         CALL CC_L1FCK3(RHO2,CTR1(1,IV),WORK(KFOCK),ISYMTR,ISYMOP,
     *                  WORK(KEND1),LWRK1)
         TIM12A = SECOND() - AUTIME
C
C-----------------------------------------
C        LT12C / H3 contribution
C-----------------------------------------
C
         CALL CC_12C3(RHO2,CTR1(1,IV),ISYMTR,WORK(KLAMDH),1,
     &                1,WORK(KEND1),LWRK1,LUIN30,'CCINT3O')
C
C---------------------------------------------
C        Save RHO2 on disc to gain work space.
C---------------------------------------------
C
         LURHO2 = -1
         CALL GPOPEN(LURHO2,'LSRHO2','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LURHO2)
         WRITE(LURHO2) (RHO2(I), I = 1,NT2SQ(ISYRES))
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
C        Need again C2(+) - C2(-)
         CALL CC_T2SQTRANSP(CTR2,ISYMTR)
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
         READ(LURHO2) (RHO2(I), I = 1,NT2SQ(ISYRES))
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
C-----------------------------------------
C        Calculate the LT22D contribution.
C-----------------------------------------
C
         AUTIME = SECOND()
         IOPT = 1
         ISYDIM = 1
         FACT = HALF
         CALL CC_22CD3(RHO2,CTR2,ISYMTR,WORK(KLAMDH),WORK(KEND1),
     *                 LWRK1,ISYDIM,LUDIM,DFIL,1,IOPT,FACT)
         TIM22D = SECOND() - AUTIME
         CALL CC_T2SQTRANSP(CTR2,ISYMTR)
C
C-----------------------------------------
C        Calculate the LT22C contribution.
C-----------------------------------------
C
         AUTIME = SECOND()
         IOPT = 1
         ISYCIM = 1
         FACT = -HALF

         CALL CC_T2SQTRANSP(RHO2,ISYMTR)
         CALL CC_22CD3(RHO2,CTR2,ISYMTR,WORK(KLAMDH),WORK(KEND1),
     *                 LWRK1,ISYCIM,LUCIM,CFIL,1,IOPT,FACT)
         CALL CC_T2SQTRANSP(RHO2,ISYMTR)
C
         TIM22C = SECOND() - AUTIME
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYRES),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of RHO1 after CC_22C:        ', RHO1N
            WRITE(LUPRI,1) 'Norm of RHO2 after CC_22C:        ', RHO2N
         ENDIF
C----------------------------------------------------------------
C     Extract the plus combination from the mixed doubles vector,
C----------------------------------------------------------------
C     use WORK(KEND1) as scratch
         CALL CC_TRIP_EXTRACT(RHO2,WORK(KEND1),ISYRES,.FALSE.)
C        Write out (+) vector
         CALL CC_WVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     &                 NT2AM(ISYMTR),IV+ITR-1,0,WORK(KEND1))
C
C   Remove (+) parts from RHO2 and CTR, now they contain (-) parts only
C
         CALL CC_T2SQSYMSCAL(RHO2,ISYMTR,ZERO)
         CALL CC_T2SQSYMSCAL(CTR2,ISYMTR,ZERO)
C
C------------------------------------------
C        Calculate the LT22C- contribution.
C------------------------------------------
C
         CALL CCSD_T2TP(CTR2,WORK(KEND1),LWRK1,ISYMTR)
         CALL CCSD_T2TP(RHO2,WORK(KEND1),LWRK1,ISYMTR)
C
         IOPT = 1
         ISYCIM = 1
         FACT = ONE
C
         CALL CC_22CD3(RHO2,CTR2,ISYMTR,WORK(KLAMDH),WORK(KEND1),
     *                 LWRK1,ISYCIM,LUCIM,CFIL,1,IOPT,FACT)
C
         CALL CCSD_T2TP(CTR2,WORK(KEND1),LWRK1,ISYMTR)
         CALL CCSD_T2TP(RHO2,WORK(KEND1),LWRK1,ISYMTR)
C
C------------------------------------------
C        Pack the RHO2(-) vector as (ai<bj)
C------------------------------------------
C
         CALL CC_TRIP_EXTRACT(RHO2,WORK(KEND1),ISYRES,.TRUE.)
         CALL DCOPY(NT2AM(ISYMTR),WORK(KEND1),1,RHO2,1)
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
C        Calculate the LT22E (-) contribution.
C-----------------------------------------
C
         CALL CCRHS_E3(DUMMY,.FALSE.,CTR2,WORK(KE1INT),WORK(KE2INT),
     &                 WORK(KENDTE),LWRKTE,ISYMTR,ISYMOP,
     &                 RHO2,.TRUE.)
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
C        For now, explicitly zero (-) diagonal!
         IF (ISYMTR.EQ.1) CALL CCLR_DIASCL(RHO2,ZERO,ISYMTR)
C------------------------------------------
C        RHO2(-) done: write out (-) vector
C------------------------------------------
C
         CALL CC_WVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     &                 NT2AM(ISYMTR),IV+ITR-1,NT2AM(ISYMTR),RHO2)
C
C--------------------------------------
C        Read CTR2+ from disc into RHO2.
C--------------------------------------
C
         DTIME = SECOND()
         CALL CC_RVEC3(LUFC2,FC2AM,2*NT2AM(ISYMTR),NT2AM(ISYMTR),
     *                 IVEC+IV-1,0,RHO2)
         DTIME = SECOND() - DTIME
         TIMIO = TIMIO + DTIME
C
C--------------------------
C        Square up CTR2(+).
C--------------------------
C
         CALL CCRHS3_R2IJ(RHO2,WORK(KENDTE),LWRKTE,ISYMTR)
         CALL CC_T2SQ(RHO2,CTR2,ISYMTR)
C-----------------------------
C        Read RHO2 (+) back in
C-----------------------------
         CALL CC_RVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     &                 NT2AM(ISYMTR),IV+ITR-1,0,RHO2)
C
C---------------------------------------------
C        Calculate the LT22E (+) contribution.
C---------------------------------------------
C        This term carries a factor of 1/2: Scale E1, E2
         CALL DSCAL(NMATAB(1)+NMATIJ(1),HALF,WORK(KE1INT),1)
C
         AUTIME = SECOND()
         CALL CCRHS_E(RHO2,CTR2,WORK(KE1INT),WORK(KE2INT),
     &                WORK(KENDTE),LWRKTE,ISYMTR,ISYMOP)
         TIM22E = TIM22E + SECOND() - AUTIME

C        Permute indices
         CALL CCRHS3_IJ(RHO2,WORK(KEND1),
     &                  LWRK1,ISYRES)
C        Write out (+) vector
         CALL CC_WVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     &                 NT2AM(ISYMTR),IV+ITR-1,0,RHO2)

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
C         WRITE(LUPRI,1) 'Norm of RHO2 after CC_11A:        ', RHO2N
      ENDIF
C
      IF ( DEBUG ) THEN
         RHO1N = DDOT(NT1AM(ISYRES),RHO1(1,IV),1,RHO1(1,IV),1)
         RHO2N = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
         WRITE(LUPRI,1) 'Norm of RHO1 after CC_12C         ', RHO1N
         WRITE(LUPRI,1) 'Norm of RHO2 after CC_12C         ', RHO2N
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
C--------------------------------------------------------------------
C     Print out result vectors - zero out single vectors in CCD-calc.
C--------------------------------------------------------------------
C
      IF (CCD) CALL DZERO(RHO1(1,IV),NT1AM(ISYRES))

         DTIME   = SECOND()
         CALL CC_WVEC(LUFR1,FRHO1,NT1AM(ISYMTR),NT1AM(ISYMTR),
     *                IV + IVEC -1,RHO1(1,IV))
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


C
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
         WRITE(LUPRI,*) 'Time used in 12C:',TIMEB
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
      CALL QEXIT(myname)
      RETURN
      END
