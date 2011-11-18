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
#include <single.h>
C
c*DECK CC_RHTR3 
      SUBROUTINE CC_RHTR3(ECURR,FRHO1,LUFR1,FRHO2,LUFR2,
     &                      FC1AM,LUFC1,FC2AM,LUFC2,
     &                      WORK, LWORK,NSIMTR,IVEC,ITR)
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!     CCS part:
!     Written by Kasper Hald and Christof Haettig.
!     February 1. 1999
!
!     CCSD part:
!     Written by Kasper Hald.
!     Spring 1999
!     Optimized by Kasper Hald and Christof Haettig
!     Summer 1999.
!
!     Biorthonormal basis for the triplet case:
!
!     1-electron: R: |(ai)> = T(ai) |HF>
!
!                 L: <(ai|  =<HF| T(ia) / 2
!
!     2-electron: R: |+(ai,bj)> = (T(ai)E(bj) + T(bj)E(ai)) |HF>
!
!                 R: |-(ai,bj)> = (T(ai)E(bj) - T(bj)E(ai)) |HF>
!
!                 L: <(ai,bj)+| = (1/8)<HF|(E(jb)T(ia)+E(ia)T(jb))
!
!                 L: <(ai,bj)+| = (1/8)<HF|(E(jb)T(ia)-E(ia)T(jb))
!
!
!     Purpose :   Calculate the right transformation of a trial
!                 vector for the triplet case.
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
      IMPLICIT NONE
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
#include <ccorb.h>
#include <ccisao.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccfield.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
#include <cclr.h>
#include <blocks.h>
#include <ccsdio.h>
#include <second.h>
#include <r12int.h>
!
      INTEGER LWORK, LUC, LUD, LUCS, LUP, LUM, IERRC, IERRD, IERRCS
      INTEGER IERRP, IERRM, LUDTIM, LUCSIM, IERRDTIM, LUB3
      INTEGER IERRCSIM, NC2, NRHO2, NRHO2P, NRHO2M
      INTEGER NRHO2TOT, NC2AM, NRHO1, KRHO1, KRHO2, KC1AM
      INTEGER KC2AM, KEND1, LWRK1, IV, KOFF1, KOFF2, NSIMTR, IVEC
      INTEGER NR1, NR2, NE1TOT, NE2TOT, NGAMTO, KT1AM, IERRC2
      INTEGER KLAMDP, KLAMDH, KDENS0, KEMAT1, KEMAT2, KFOCK0
      INTEGER KLAMPC, KLAMHC, KDENSC, KFOCKC, KFOCKT, KFOCKC2
      INTEGER ISYMH, IC, KAMC1, KLPC1, KLHC1, KDNSC1, IF
      INTEGER KENDS2, LWRKS2, KRHO22, KC22AM, NTOSYM,IERRB3
      INTEGER KCCFB1, KINDXB, KODCL1, KODCL2, KODBC1, KODBC2
      INTEGER KRDBC1, KRDBC2, KODPP1, KODPP2, KRDPP1, KRDPP2
      INTEGER KFREE, LFREE, ICDEL1, ICDEL2, KENDSV, LWRKSV
      INTEGER ISYMD1, NTOT, ILLL, NUMDIS, KRECNR
      INTEGER KFCKC1, KFCKC2, KEI1, KEI2, IDEL, IDEL2, ISYMD
      INTEGER ISYMT12, ISYMC2, ISYMT1, ISYMT2
      INTEGER ISYDIS, ISYAIK, ITR, KXINT, KEND2, LWRK2, ISYDEN
      INTEGER KDSRHF, KEND3, LWRK3, ISYMLP, KSCRCM, KSCRCM2
      INTEGER KEND4, LWRK4, ICON, ISYMLH, IOPT, ISYMM, ISYMPC
      INTEGER ISYMHC, ISYMP1, ISYMH1, KGAMC, ISYFAO
      INTEGER ISYMPA, ISYMHO, LFOCK, KGAMC1, ISYMBF, ISYFCK, LUE1
      INTEGER LUE2, ISYMEI, ISYMIM, ISYVEC, LUGAM, ISYGAM, ISYCIM
      INTEGER IVECNR, ISYDIM, KGAMMC, KFCKAO, LEND2
      INTEGER KGAMI, LUBF, KHSCOUNTER, ISYM0
      INTEGER INDEXA(MXCORB)
      INTEGER LUFC1, LUFC2, LUFR1, LUFR2, ISIDE
      INTEGER KRHO3, KC2AMP, KC2AMM, KENDCC3, LWRKCC3
      INTEGER LUSRTT1, LUSRTT2, LUSRTT3, LUCKJDT1, LUCKJDT2, LUCKJDT3
      INTEGER LU3VI, LU3VI2, LUTOC
      INTEGER IDUM
!
#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL ZERO, ONE, TWO, HALF, FACTD, RHO1N, RHO2N, FACTE1, FACTE2
      REAL FACTC, FACG, FACH, FAKE, PHONEY, NFACC, PFACC
      REAL XNGAM, XMTWO, XMONE, XMHALF, FACTD2, ECURR
      REAL TIMALL, TIMER1, TIMER2, TIMIO, TIMRDAO, FF, DTIME
      REAL RHO1N, DDOT, TIMFCK, TIMAOD, TIMT2SQ, TIMLAM, TIMTRBT
      REAL TIMT2AO, TIMD, TIMC, TIMEI, TIMBF, TIMG, TIMH, TIMFCKMO
      REAL TIMT2MO, TIMA, TIMCIO, TIMI, TIMDIO, TIMCCSD, TIME1C1
      REAL TIMT2TR, ECURR2, TIMOM32, TIMT3I
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ZERO, ONE, TWO, HALF, FACTD, RHO1N, RHO2N
      DOUBLE PRECISION FACTE1, FACTE2, XMHALF, FACTD2, ECURR
      DOUBLE PRECISION FACTC, FACG, FACH, FAKE, PHONEY, NFACC, PFACC
      DOUBLE PRECISION XNGAM, XMTWO, XMONE
      DOUBLE PRECISION TIMALL, TIMER1, TIMER2, TIMIO, TIMRDAO, FF
      DOUBLE PRECISION DDOT, DTIME, TIMFCK, TIMAOD
      DOUBLE PRECISION TIMT2SQ, TIMLAM, TIMTRBT, TIMT2AO, TIMD
      DOUBLE PRECISION TIMC, TIMEI, TIMBF, TIMG, TIMH, TIMFCKMO
      DOUBLE PRECISION TIMT2MO, TIMA, TIMCIO, TIMI, TIMDIO, TIMCCSD
      DOUBLE PRECISION TIME1C1, TIMT2TR, ECURR2, TIMOM32, TIMT3I
#endif
!
      PARAMETER (ZERO = 0.0D00, ONE = 1.0D00, TWO = 2.0D00)
      PARAMETER (HALF= 0.5D00, XMTWO= -2.0D00, XMONE = -1.0D00)
      PARAMETER (XMHALF = -0.5D00)
!
      CHARACTER*8 FRHO1, FRHO2, FC1AM, FC2AM
      CHARACTER*8 FNCKJDT1, FNCKJDT2, FNCKJDT3, FNTOC, FN3VI2
      CHARACTER*7 PFIL, CTFIL, DTFIL, CSFIL
      CHARACTER*7 DTIMFIL, BFFIL 
      CHARACTER*6 CSIMFIL, FN3VI
      CHARACTER*10 MODEL, FNSRTT1, FNSRTT2, FNSRTT3
!
      LOGICAL FCKCON, ETRAN, ANTISYM, NORMALCONT, PIJCONT
      LOGICAL NEWGS, KHLOG, LDEBUG
!
      CALL QENTER('CC_RHTR3')
!
      LDEBUG = .FALSE.
!
!---------------------------------------------------------
!
      IF (IPRINT .GT. 15) THEN
        CALL AROUND(' START OF CC_RHTR3 ')
        IF (DIRECT) WRITE(LUPRI,*) ' Atomic direct calculation '
      ENDIF
!
      IF ( .NOT. (CCS .OR. CC2 .OR. CCSD .OR. CC3) ) THEN
!
        WRITE(LUPRI,*) 'Unknown Coupled Cluster calculation (CC_RHTR3)'
        WRITE(LUPRI,*) 'Triplet calc. can only (for the moment) '
        WRITE(LUPRI,*) 'be made in CCS, CC2, CCSD or CC3! '
	CALL QUIT(' ')
      ENDIF
!
!----------------
!     Open files.
!----------------
!
      IF (DUMPCD .AND. (.NOT. CCS )) THEN
         LUB3 = -1
         LUC  = -1
         LUD  = -1
         LUCS = -1
         LUP  = -1
!
         BFFIL = 'CCLR_B3'
         CTFIL = 'CCLR_CT'
         DTFIL = 'CCLR_DT'
         CSFIL = 'CCLR_CS'
         PFIL  = 'CCLR_PT'
!
         CALL WOPEN2(LUB3,BFFIL,64,0)
         CALL WOPEN2(LUC,CTFIL,64,0)
         CALL WOPEN2(LUD,DTFIL,64,0)
         CALL WOPEN2(LUCS,CSFIL,64,0)
         CALL WOPEN2(LUP,PFIL,64,0)
!
         IF (.NOT. CC2) THEN
!
            LUDTIM  = -1
            LUCSIM  = -1
!
            DTIMFIL  = 'PMAT_DT'
            CSIMFIL  = 'PMAT_C'
!
            CALL WOPEN2(LUDTIM,DTIMFIL,64,0)
            CALL WOPEN2(LUCSIM,CSIMFIL,64,0)
!
         ENDIF
!
      END IF
!
      IF (CCSDT) THEN
!
         LUSRTT1  = -1
         LU3VI    = -1
         LU3VI2   = -1
         LUCKJDT1 = -1
         LUCKJDT2 = -1
         LUCKJDT3 = -1
         LUSRTT2  = -1
         LUSRTT3  = -1
         LUTOC    = -1
!
         FNSRTT1  = 'CC3_SORTT1'
         FN3VI    = 'CC3_VI'
         FN3VI2   = 'CC3_VI12'
         FNCKJDT1 = 'CKJDELT1'
         FNCKJDT2 = 'CKJDELT2'
         FNCKJDT3 = 'CKJDELT3'
         FNSRTT2  = 'CC3_SORTT2'
         FNSRTT3  = 'CC3_SORTT3'
         FNTOC    = 'CCSDT_OC'
!
         CALL WOPEN2(LUSRTT1,FNSRTT1,64,0)
         CALL WOPEN2(LU3VI,FN3VI,64,0)
         CALL WOPEN2(LU3VI2,FN3VI2,64,0)
         CALL WOPEN2(LUCKJDT1,FNCKJDT1,64,0)
         CALL WOPEN2(LUCKJDT2,FNCKJDT2,64,0)
         CALL WOPEN2(LUCKJDT3,FNCKJDT3,64,0)
         CALL WOPEN2(LUSRTT2,FNSRTT2,64,0)
         CALL WOPEN2(LUSRTT3,FNSRTT3,64,0)
         CALL WOPEN2(LUTOC,FNTOC,64,0)
!
      ENDIF
!
!-----------------------------------------------
!     Time Initialiation.
!     UPDATE these times to triplet case.
!-----------------------------------------------
!
      TIMALL    = SECOND()
!
      TIMA      = 0.0D00
      TIMAOD    = 0.0D00
      TIMBF     = 0.0D00
      TIMC      = 0.0D00
      TIMCIO    = 0.0D00
      TIMD      = 0.0D00
      TIMDIO    = 0.0D00
      TIME1C1   = 0.0D00
      TIMEI     = 0.0D00
      TIMER1    = 0.0D00
      TIMER2    = 0.0D00
      TIMFCK    = 0.0D00
      TIMFCKMO  = 0.0D00
      TIMG      = 0.0D00
      TIMH      = 0.0D00
      TIMI      = 0.0D00
      TIMIO     = 0.0D00
      TIMLAM    = 0.0D00
      TIMOM32   = 0.0D00
      TIMRDAO   = 0.0D00
      TIMT2AO   = 0.0D00
      TIMT2MO   = 0.0D00
      TIMT2SQ   = 0.0D00
      TIMT2TR   = 0.0D00
      TIMT3I    = 0.0D00
      TIMTRBT   = 0.0D00
!
!--------------------------------
! For CCS there is no C2 vectors.
!--------------------------------
!
      NC2 = 1
      IF ( CCS ) NC2 = 0
!
!--------------------------------------------------------------
!     Work space allocation for rho and C.
!--------------------------------------------------------------
!
      NRHO2 = MAX(NT2AM(ISYMTR),NT2AMA(ISYMTR),
     *            2*NT2ORT(ISYMTR),NT2AM(1))
      IF (.NOT. CCS) THEN
         NRHO2P = MAX(NT2AM(ISYMTR),NT2ORT(ISYMTR))
         NRHO2M = MAX(NT2AMA(ISYMTR),2*NT2ORT(ISYMTR))
         NRHO2TOT = NRHO2P + NRHO2M
      ENDIF
      IF ( CCS ) NRHO2 = 0
!
      NC2AM = MAX(NT2SQ(ISYMTR),NT2SQ(1),
     *        NT2AM(ISYMTR)+2*NT2ORT(1))
      IF ( CCS ) NC2AM = 0
!
      NRHO1 = NT1AM(ISYMTR)*NSIMTR
!
      KC1AM = 1
      KRHO2 = KC1AM + NRHO1
      KRHO1 = KRHO2 + NRHO2
      KC2AM = KRHO1 + NRHO1
      KEND1 = KC2AM + NC2AM
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LE. 0 )
     *    CALL QUIT('Too little workspace in CC_RHTR3 (1) ')
!
!------------------------------------------------------
!        Read the CC trial vectors from disk.
!        First C2+ and then C2-
!------------------------------------------------------
!
      DO 25 IV = 1, NSIMTR
         KAMC1  = KC1AM + NT1AM(ISYMTR)*(IV - 1)
         CALL CC_RVEC(LUFC1,FC1AM,NT1AM(ISYMTR),
     *                NT1AM(ISYMTR),
     *                IV+IVEC-1,WORK(KAMC1))
   25 CONTINUE
!
      IF ((.NOT.CCS).AND.( MINSCR).AND.(IPRINT.GT.45)) THEN
         CALL CC_RVEC3(LUFC2,FC2AM,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                 NT2AM(ISYMTR),IVEC,0,WORK(KRHO2))
      ENDIF
!
!-----------------------------------------------
!        Print (C2+) amplitudes (packed) for
!        (ai)<(bj) AND i<j
!-----------------------------------------------
!
      IF (IPRINT .GT. 45) THEN
         DO 30 IV = 1, NSIMTR
            KAMC1  = KC1AM + NT1AM(ISYMTR)*(IV - 1)
            CALL AROUND('CC_RHTR3: (C2+) packed. Both (ai<bj AND i<j)')
            WRITE(LUPRI,*) 'Vector nr is = ',IV+IVEC-1
            CALL CC_PRP(WORK(KAMC1),
     *                  WORK(KRHO2),ISYMTR,0,NC2)
   30    CONTINUE
      ENDIF
!
!-----------------------------------------------
!        "Square up" C2+ to only (ai)<bj)
!-----------------------------------------------
!
      IF ((.NOT. CCS).AND.(IPRINT.GT.45)) THEN
         CALL CCRHS3_R2IJ(WORK(KRHO2),WORK(KEND1),
     *                    LWRK1,ISYMTR)
      ENDIF
!
!-----------------------------------------------
!        Print (C2+) amplitudes (packed) for
!        (ai)<(bj).
!-----------------------------------------------
!
      IF (IPRINT .GT. 45) THEN
         DO IV = 1, NSIMTR
            KAMC1  = KC1AM + NT1AM(ISYMTR)*(IV - 1)
            CALL AROUND('CC_RHTR3: (C2+) packed. ONLY (ai<bj)')
            WRITE(LUPRI,*) 'Vector nr is = ',IV+IVEC-1
            CALL CC_PRP(WORK(KAMC1),
     *                  WORK(KRHO2),ISYMTR,0,NC2)
         ENDDO
      ENDIF
!
!-----------------------------------------------
!           Prepare the (C1,C2+)-amplitudes.
!-----------------------------------------------
!
      IF ((.NOT.CCS).AND.(IPRINT.GT.45)) THEN
         CALL CCLR_DIASCL(WORK(KRHO2),ZERO,ISYMTR)
         CALL CC_T2SQ(WORK(KRHO2),WORK(KC2AM),ISYMTR)
      ENDIF
!
!----------------------------------------------------
!           Print the squared vectors (C2+).
!----------------------------------------------------
!
      IF (IPRINT.GT.45) THEN
         CALL AROUND('CC_RHTR3: (C2+) squared ')
         CALL CC_PRSQ(WORK(KC1AM),WORK(KC2AM),ISYMTR,0,NC2)
      ENDIF
!
      IF (LDEBUG) THEN
         DO IV = 1, NSIMTR
           KAMC1 = KC1AM + NT1AM(ISYMTR)*(IV - 1)
           RHO1N = DDOT(NT1AM(ISYMTR),WORK(KAMC1),1,WORK(KAMC1),1)
           WRITE(LUPRI,1) 'Norm of C1AM - first in CC_RHTR3 :',RHO1N
         ENDDO
!
         RHO2N = DDOT(NT2SQ(ISYMTR),WORK(KC2AM),1,WORK(KC2AM),1)
         WRITE(LUPRI,1) 'Norm of C2+  - first in CC_RHTR3 :',RHO2N
!     
      ENDIF
!
!------------------------------------------------
!           Read only (C2-), since we have C1
!------------------------------------------------
!
      IF ((.NOT.CCS).AND.( MINSCR)) THEN
         CALL CC_RVEC3(LUFC2,FC2AM,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                 NT2AMA(ISYMTR),IVEC,NT2AM(ISYMTR),WORK(KRHO2))
      ENDIF
!
!-----------------------------------------------
!        Print (C1,C2-) amplitudes (Packed)
!-----------------------------------------------
!
      IF (.NOT. CCS) THEN
         IF (IPRINT .GT. 45) THEN
            DO 40 IV = 1, NSIMTR
               KAMC1  = KC1AM + NT1AM(ISYMTR)*(IV - 1)
               CALL AROUND('CC_RHTR3: (C2-) packed ')
               WRITE(LUPRI,*) 'Vector nr is = ',IV+IVEC-1
               CALL CC_PRP(WORK(KAMC1),
     *                     WORK(KRHO2),ISYMTR,1,NC2)
   40       CONTINUE
         ENDIF
!
!-----------------------------------------------
!           Prepare the (C-)-amplitudes.
!-----------------------------------------------
!
         CALL CCLR_DIASCL(WORK(KRHO2),ZERO,ISYMTR)
         CALL CC_T2SQ3(WORK(KRHO2),WORK(KC2AM),ISYMTR)
!
!---------------------------------------------------
!           Print the squared (C2-) vectors.
!---------------------------------------------------
!
         IF (IPRINT.GT.45) THEN
            CALL AROUND('CC_RHTR3: (C1,C2-) squared ')
            CALL CC_PRSQ(WORK(KC1AM),WORK(KC2AM),ISYMTR,1,NC2)
         ENDIF
!
         IF (LDEBUG) THEN
           RHO2N = DDOT(NT2SQ(ISYMTR),WORK(KC2AM),1,
     *                   WORK(KC2AM),1)
           WRITE(LUPRI,1) 'Norm of C2-  first in CC_RHTR3   :',RHO2N
         ENDIF
!
      ENDIF
!
!---------------------------
!           Zero rho vector.
!---------------------------
!
      CALL DZERO(WORK(KRHO1),NT1AM(ISYMTR)*NSIMTR)
      CALL DZERO(WORK(KRHO2),NRHO2)
!
!------------------------------
!           File read and save.
!------------------------------
!
            DO 47 IV = 1, NSIMTR
               NR1 = IV + IVEC - 1
!
               IF (.NOT. (CCS .OR. CC2)) THEN
                  NR2 = IV
               ELSE
                  NR2 = NR1
               ENDIF
               KOFF1 = KRHO1 + NT1AM(ISYMTR)*(IV - 1)
               CALL CC_WVEC(LUFR1,FRHO1,NT1AM(ISYMTR),
     *                   NT1AM(ISYMTR),NR1,WORK(KOFF1))
               IF (.NOT.CCS) THEN
                 CALL CC_WVEC3(LUB3,BFFIL,NRHO2TOT,
     *                         NRHO2P,NR2,0,WORK(KRHO2))
                 CALL CC_WVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                         NT2AM(ISYMTR),NR2,0,WORK(KRHO2))
                 CALL CC_WVEC3(LUB3,BFFIL,NRHO2TOT,
     *                         NRHO2M,NR2,NRHO2P,WORK(KRHO2))
                 CALL CC_WVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                         NT2AMA(ISYMTR),NR2,NT2AM(ISYMTR),
     *                         WORK(KRHO2))
               ENDIF
   47       CONTINUE
!
!----------------------------------------------------------------
!     Work space allocation for C1 independent arrays:
!----------------------------------------------------------------
!
      NE1TOT = MAX(NEMAT1(ISYMTR)*NSIMTR,NEMAT1(1))
      NE2TOT = MAX(NMATIJ(ISYMTR)*NSIMTR,NMATIJ(1))
!
      IF (CCS) THEN
             NGAMTO = 2
             NE1TOT = 2
             NE2TOT = 2
      ENDIF
!
      IF ( IPRINT .GT. 25)
     &   WRITE(LUPRI,*) ' In CC_RHTR3 after omg/rho work is:',LWORK

      KT1AM   = KEND1
      KLAMDP  = KT1AM  + NT1AM(ISYMOP)
      KLAMDH  = KLAMDP + NLAMDT
      KDENS0  = KLAMDH + NLAMDT
      KFOCK0  = KDENS0 + N2BAST
      KEMAT1  = KFOCK0 + N2BST(ISYMOP)
      KEMAT2  = KEMAT1 + NE1TOT
      KEND1   = KEMAT2 + NE2TOT
      LWRK1   = LWORK  - KEND1
!
!------------------------------------------------------------
!     Work space allocation for C1 dependent intermediates.
!------------------------------------------------------------
!
      KLAMPC = KEND1
      KLAMHC = KLAMPC + NGLMDT(ISYMTR)*NSIMTR
      KDENSC = KLAMHC + NGLMDT(ISYMTR)*NSIMTR
      KFOCKC = KDENSC + N2BST(ISYMTR)*NSIMTR
      KFOCKT = KFOCKC + N2BST(ISYMTR)*NSIMTR
      KFOCKC2= KFOCKT + MAX(N2BST(ISYMTR),N2BST(ISYMOP))
      KEND1  = KFOCKC2 + N2BST(ISYMTR)*NSIMTR
      LWRK1  = LWORK  - KEND1
!
      IF ( IPRINT .GT. 25)
     &   WRITE(LUPRI,*)
     &        ' In CC_RHTR3 after the 1. alloc. work. left:',LWRK1
!
      IF (LWRK1 .LE. 0)
     &        CALL QUIT('Too little workspace in CC_RHTR3 (2) ')
!
!---------------------------------------------------
!     Initialize some C1 independent intermediates.
!---------------------------------------------------
!
      CALL DZERO(WORK(KT1AM),NT1AM(1))
      CALL DZERO(WORK(KEMAT1),NE1TOT)
      CALL DZERO(WORK(KEMAT2),NE2TOT)
      CALL DZERO(WORK(KDENS0),N2BST(ISYMOP))
      CALL DZERO(WORK(KFOCK0),N2BST(ISYMOP))
!
!-------------------------------------------
!     Initialize C1 dependent intermediates.
!-------------------------------------------
!
      CALL DZERO(WORK(KFOCKC),N2BST(ISYMTR)*NSIMTR)
      CALL DZERO(WORK(KFOCKC2),N2BST(ISYMTR)*NSIMTR)
      CALL DZERO(WORK(KFOCKT),MAX(N2BST(ISYMTR),N2BST(ISYMOP)))
      CALL DZERO(WORK(KDENSC),N2BST(ISYMTR)*NSIMTR)
!
!----------------------------------------------
!    The T1 amp. are zero for the CCS model.
!    For CCSD read the T1-amplitudes from file.
!----------------------------------------------
!
!
      IF (.NOT. CCS) THEN
!
         DTIME = SECOND()
         IOPT = 1
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),DUMMY)
         DTIME = SECOND() - DTIME
         TIMIO = TIMIO + DTIME
!
      ENDIF
!
!----------------------------------
!     Calculate the lamda matrices.
!----------------------------------
!
      DTIME = SECOND()
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     &            WORK(KEND1),LWRK1)
      DTIME = SECOND() - DTIME
      TIMLAM = DTIME + TIMLAM
!
      IF (IPRINT .GT.100) THEN
         CALL AROUND('Usual triplet Lambda matrices ')
         CALL CC_PRLAM(WORK(KLAMDP),WORK(KLAMDH),1)
      ENDIF
!
!----------------------------------
!     Calculate the density matrix.
!----------------------------------
!
      DTIME = SECOND()
      ISYMH  = 1
      IC     = 1
      CALL CC_AODENS(WORK(KLAMDP),WORK(KLAMDH),WORK(KDENS0),ISYMH,
     &               IC,WORK(KEND1),LWRK1)
      DTIME = SECOND() - DTIME
      TIMAOD = TIMAOD + DTIME
!
      IF (IPRINT .GT. 45) THEN
         CALL AROUND('CC_RHTR3: Usual triplet Lamda density matrix')
         CALL CC_PRAODEN(WORK(KDENS0),1)
      ENDIF
!
      DO 50 IV = 1, NSIMTR
!
!----------------------------------------------------
!        Allocation specifications.
!----------------------------------------------------
!
         KLPC1  = KLAMPC + NGLMDT(ISYMTR)*(IV-1)
         KLHC1  = KLAMHC + NGLMDT(ISYMTR)*(IV-1)
         KDNSC1 = KDENSC + N2BST(ISYMTR)*(IV-1)
         KAMC1  = KC1AM  + NT1AM(ISYMTR)*(IV-1)
         KOFF1  = KRHO1  + NT1AM(ISYMTR)*(IV - 1)
!
!-----------------------------------------
!        Transform lamda matrices with C1.
!-----------------------------------------
!
         CALL CCLR_LAMTRA(WORK(KLAMDP),WORK(KLPC1),
     &                    WORK(KLAMDH),WORK(KLHC1),
     &                    WORK(KAMC1),ISYMTR)
!
         IF (IPRINT .GT. 45) THEN
            CALL AROUND('C1 transformed triplet Lambda matrices ')
            CALL CC_PRLAM(WORK(KLPC1),WORK(KLHC1),ISYMTR)
         ENDIF
!
!-----------------------------------------------
!        Print C1 and rho1 amplitudes.
!-----------------------------------------------
!
         IF (IPRINT .GT. 50) THEN
            CALL AROUND('START OF CC_RHTR3: (C1) ')
            CALL CC_PRSQ(WORK(KAMC1),WORK,ISYMTR,1,0)
            CALL AROUND('START OF CC_RHTR3: (RHO1) ')
            CALL CC_PRSQ(WORK(KOFF1),WORK,ISYMTR,1,0)
         ENDIF
!
!----------------------------------------------------
!        Calculate the C1 transformed density matrix.
!        no core contribution ic = 0
!----------------------------------------------------
!
         DTIME = SECOND()
         ISYMH  = ISYMTR
         IC     = 0
         CALL CC_AODENS(WORK(KLAMDP),WORK(KLHC1),WORK(KDNSC1),ISYMH,
     &               IC,WORK(KEND1),LWRK1)
         DTIME = SECOND() - DTIME
         TIMAOD = TIMAOD + DTIME
!
         IF (IPRINT .GT. 50) THEN
            CALL AROUND('CC_RHTR3: C1 trans. Lamda density matrix')
            CALL CC_PRAODEN(WORK(KDNSC1),ISYMTR)
         ENDIF
!
  50  CONTINUE
!
!
!------------------------------------------------
!     Read one-electron integrals in Fock-matrix.
!------------------------------------------------
!
      TIMFCK = SECOND()
      CALL CCRHS_ONEAO(WORK(KFOCK0),WORK(KEND1),LWRK1)
      TIMFCK = SECOND() - TIMFCK
!
!-------------------------------------------------------
!     Read one-electron integrals into Fock-matrix for
!     finite field.
!-------------------------------------------------------
!
      DO IF = 1, NFIELD
         DTIME = SECOND()
         FF    = EFIELD(IF)
         CALL CC_ONEP(WORK(KFOCK0),WORK(KEND1),LWRK1,FF,1,LFIELD(IF))
         DTIME = DTIME - SECOND()
         TIMFCK= TIMFCK + DTIME
      ENDDO
!
!
!====================================================
!     Initialize integral program ERI:
!====================================================
!
      IF (.NOT. (CC2 .OR. CCS)) THEN
         KENDS2 = KEND1 + NT2AM(ISYMTR)
         LWRKS2 = LWRK1 - NT2AM(ISYMTR)
      ELSE
         KENDS2 = KEND1
         LWRKS2 = LWRK1
      ENDIF
!
!--------------------------------------------------------
!      Allocate space for a second RHO2 vector.
!      Allocate space for a second C2 vector.
!      Initialize second RHO2 and second C2.
!      Only to be used in the loop over delta.
!--------------------------------------------------------
!
      IF (.NOT. (CCS)) THEN
         KRHO22 = KEND1
         KC22AM = KRHO22 + NRHO2P
         KEND1  = KC22AM + NC2AM
         LWRK1  = LWORK - KEND1
!
         IF (LWRK1 .LT. 0) THEN
             CALL QUIT('Insufficient core in CC_RHTR3 (RHO22)')
         ENDIF
!
         CALL DZERO(WORK(KRHO22),NRHO2P)
!
!--------------------------------------------------
!                  Readin C2+ in RHO2
!--------------------------------------------------
!
         IF (MINSCR) THEN
            DTIME = SECOND()
!
            CALL CC_RVEC3(LUFC2,FC2AM,
     *                    NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                    NT2AM(ISYMTR),IVEC,0,WORK(KRHO2))
            DTIME = SECOND() - DTIME
            TIMIO = TIMIO + DTIME
!
!---------------------------------------------------
!                 Square up C2+
!---------------------------------------------------
!
            DTIME = SECOND()
            CALL CCRHS3_R2IJ(WORK(KRHO2),WORK(KEND1),
     *                       LWRK1,ISYMTR)
            CALL CCLR_DIASCL(WORK(KRHO2),ZERO,ISYMTR)
            CALL CC_T2SQ(WORK(KRHO2),WORK(KC22AM),ISYMTR)
!
!-------------------------------------
!        Zero the rho2 vector.
!-------------------------------------
!
            CALL DZERO(WORK(KRHO2),NRHO2)
!
            DTIME = SECOND() - DTIME
            TIMT2SQ = TIMT2SQ + DTIME
!
            IF (IPRINT.GT.50) THEN
               CALL AROUND('CC_RHTR3: C2+ vector readin/sq.')
               CALL CC_PRSQ(WORK(KAMC1),WORK(KC22AM),
     *                     ISYMTR,0,1)
            ENDIF
!
         ENDIF
!
      ENDIF
!
!------------------------------------------------------
!      Continue integral program ERI
!------------------------------------------------------
!
      IF (DIRECT) THEN
         DTIME  = SECOND()
!
         IF (HERDIR) THEN
           CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
         ELSE
           KCCFB1 = KEND1
           KINDXB = KCCFB1 + MXPRIM*MXCONT
           KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
           LWRK1  = LWORK  - KEND1
!
           CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     &                 KODPP1,KODPP2,KRDPP1,KRDPP2,
     &                 KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     &                 WORK(KEND1),LWRK1,IPRERI)
           KEND1 = KFREE
           LWRK1 = LFREE
         END IF

         DTIME  = SECOND() - DTIME
         TIMER1 = TIMER1 + DTIME

         NTOSYM = 1
      ELSE
         NTOSYM = NSYM
      ENDIF
!
      ICDEL1 = 0
      ICDEL2 = 0
      KENDSV = KEND1
      LWRKSV = LWRK1
!
!====================================================
!     Start the loop over distributions of integrals.
!====================================================
!
      DO 100 ISYMD1 = 1,NTOSYM
!
         IF (DIRECT) THEN
            IF (HERDIR) THEN
               NTOT = MAXSHL
            ELSE
               NTOT = MXCALL
            END IF
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
!
         DO 110 ILLL = 1,NTOT
!
!------------------------------------------
!        If direct calculate the integrals.
!------------------------------------------
!
            IF (DIRECT) THEN
!
               KEND1 = KENDSV
               LWRK1 = LWRKSV
!
               DTIME  = SECOND()
               IF (HERDIR) THEN
                 CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                       IPRERI)
               ELSE
                 CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     &                       WORK(KODCL1),WORK(KODCL2),
     &                       WORK(KODBC1),WORK(KODBC2),
     &                       WORK(KRDBC1),WORK(KRDBC2),
     &                       WORK(KODPP1),WORK(KODPP2),
     &                       WORK(KRDPP1),WORK(KRDPP2),
     &                       WORK(KCCFB1),WORK(KINDXB),
     &                       WORK(KEND1), LWRK1,IPRERI)
               END IF
               DTIME  = SECOND() - DTIME
               TIMER2 = TIMER2 + DTIME
!
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient core in CC_RHTR3')
               END IF
!
            ELSE
               NUMDIS = 1
            ENDIF
!
!----------------------------------------------------
!           Loop over number of trial vectors NSIMTR.
!           Calculate C2+/C2- tilde terms AND
!           C2-/C2+ independent tilde terms
!----------------------------------------------------
!
            DO 120 IV = 1, NSIMTR
!
               KAMC1  = KC1AM   + NT1AM(ISYMTR) *(IV-1)
               KOFF1  = KRHO1   + NT1AM(ISYMTR) *(IV-1)
               KLPC1  = KLAMPC  + NGLMDT(ISYMTR)*(IV-1)
               KLHC1  = KLAMHC  + NGLMDT(ISYMTR)*(IV-1)
               KDNSC1 = KDENSC  + N2BST(ISYMTR) *(IV-1)
               KFCKC1 = KFOCKC  + N2BST(ISYMTR) *(IV-1)
               KFCKC2 = KFOCKC2 + N2BST(ISYMTR) *(IV-1)
               KEI1   = KEMAT1  + NEMAT1(ISYMTR)*(IV-1)
               KEI2   = KEMAT2  + NMATIJ(ISYMTR)*(IV-1)
!
               IF ((.NOT.MINSCR).AND.(.NOT.CCS)) THEN
!
!-----------------------------------------------
!                 Readin C2- amplitude in RHO2.
!-----------------------------------------------
!
                  DTIME = SECOND()
                  CALL CC_RVEC3(LUFC2,FC2AM,
     *                          NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                          NT2AMA(ISYMTR),IVEC+IV-1,
     *                          NT2AM(ISYMTR),WORK(KRHO2))
                  DTIME = SECOND() - DTIME
                  TIMIO = TIMIO + DTIME
!
!-----------------------------------------
!                 Square up C2- amplitudes.
!-----------------------------------------
!
                  DTIME = SECOND()
                  CALL CCLR_DIASCL(WORK(KRHO2),ZERO,ISYMTR)
                  CALL CC_T2SQ3(WORK(KRHO2),WORK(KC2AM),ISYMTR)
                  DTIME = SECOND() - DTIME
                  TIMT2SQ = TIMT2SQ + DTIME
!
                  IF (IPRINT.GT.50) THEN
                     CALL AROUND('CC_RHTR3: C2- vector readin')
                     CALL CC_PRSQ(WORK(KAMC1),WORK(KC2AM),
     &                           ISYMTR,0,1)
                  ENDIF
!
!--------------------------------------------------
!                  Readin C2+ in RHO2
!--------------------------------------------------
!
                  DTIME = SECOND()
!
                  CALL CC_RVEC3(LUFC2,FC2AM,
     *                          NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                          NT2AM(ISYMTR),IVEC+IV-1,0,WORK(KRHO2))
                  DTIME = SECOND() - DTIME
                  TIMIO = TIMIO + DTIME
!
!---------------------------------------------------
!                 Square up C2+
!---------------------------------------------------
!
                  DTIME = SECOND()
                  CALL CCRHS3_R2IJ(WORK(KRHO2),WORK(KEND1),
     *                             LWRK1,ISYMTR)
                  CALL CCLR_DIASCL(WORK(KRHO2),ZERO,ISYMTR)
                  CALL CC_T2SQ(WORK(KRHO2),WORK(KC22AM),ISYMTR)
                  DTIME = SECOND() - DTIME
                  TIMT2SQ = TIMT2SQ + DTIME
!
                  IF (IPRINT.GT.50) THEN
                     CALL AROUND('CC_RHTR3: C2+ vector readin/sq.')
                     CALL CC_PRSQ(WORK(KAMC1),WORK(KC22AM),
     &                           ISYMTR,0,1)
                  ENDIF
!
!---------------------------------------------
!                 Read in result vectorS.
!---------------------------------------------
!
                  DTIME = SECOND()
!
                  CALL CC_RVEC3(LUB3,BFFIL,NRHO2TOT,NRHO2M,
     &                          IV+ITR-1,NRHO2P,WORK(KRHO2))
                  CALL CC_RVEC3(LUB3,BFFIL,NRHO2TOT,NRHO2P,
     *                          IV+ITR-1,0,WORK(KRHO22))
                  DTIME = SECOND() - DTIME
                  TIMIO = TIMIO + DTIME
!
                  IF (IPRINT .GT. 50) THEN
                     CALL AROUND('CC_RHTR3 : BF- vector readin')
                     CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,1,1)
                     CALL AROUND('CC_RHTR3 : BF+ vector readin')
                     CALL CC_PRP(WORK(KOFF1),WORK(KRHO22),ISYMTR,0,1)
                  ENDIF
!
               ENDIF
!
!------------------------------------------------------
!           Loop over number of distributions in disk.
!------------------------------------------------------
!
               DO 125 IDEL2 = 1,NUMDIS
!
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
!
                  ISYDIS = MULD2H(ISYMD,ISYMOP)
                  ISYAIK = MULD2H(ISYDIS,ISYMTR)
!
!----------------------------------------------------------
!              Calculate adresses for c,cio,d,dio routines.
!----------------------------------------------------------
!
                  IF ( IV .EQ. 1) THEN
                     IT2DEL(IDEL) = ICDEL1
                     ICDEL1 = ICDEL1 + NT2BCD(ISYDIS)
                  ENDIF
!
                  IT2DLR(IDEL,IV) = ICDEL2
                  ICDEL2 = ICDEL2 + NT2BCD(ISYAIK)
!
!------------------------------------------
!              Work space allocation 
!------------------------------------------
!
                  KXINT  = KEND1
                  KEND2  = KXINT + NDISAO(ISYDIS)
                  LWRK2  = LWORK - KEND2
!
                  IF ( IPRINT .GT. 55)
     &               WRITE(LUPRI,*) 
     &               ' In CC_RHTR3 LWRK2 after 2. alloc.:',LWRK2
!
                  IF (LWRK2 .LT. 0) THEN
                     WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',
     &                              LWORK
                     CALL QUIT('Insufficient space in CC_RHTR3 (1)')
                  ENDIF
!
!-----------------------------------------
!              Read in batch of integrals.
!-----------------------------------------
!
                  DTIME   = SECOND()
                  CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     &                        WORK(KRECNR),DIRECT)
                  DTIME   = SECOND() - DTIME
                  TIMRDAO = TIMRDAO  + DTIME
!
!-----------------------------------------------------------
!              Calculate transformed integrals used in t3am.
!-----------------------------------------------------------
!
                  IF (CCSDT) THEN
!
                     DTIME = SECOND()
                     CALL CC3_T3INT3(WORK(KXINT),WORK(KLAMDP),
     *                              WORK(KLAMDH),WORK(KAMC1),ISYMTR,
     *                              WORK(KEND2),LWRK2,IDEL,ISYMD,
     *                              LUSRTT1,FNSRTT1,LUSRTT2,FNSRTT2,
     *                              LUSRTT3,FNSRTT3,LUCKJDT1,FNCKJDT1,
     *                              LUCKJDT2,FNCKJDT2,LUCKJDT3,FNCKJDT3)
                     DTIME   = SECOND() - DTIME
                     TIMT3I  = TIMT3I   + DTIME
!
               ENDIF
!
!------------------------------------------------------------------
!           Calculate the AO-Fock matrix.
!           Not neccesary if RSPIM is used. So only for CCS
!------------------------------------------------------------------
!
                  IF (CCS) THEN
                     ISYDEN  = 1
                     IF (IV .EQ. 1) THEN
                        CALL CC_AOFOCK(WORK(KXINT),WORK(KDENS0),
     &                                 WORK(KFOCK0),WORK(KEND2),
     &                                 LWRK2,IDEL,ISYMD,.FALSE.,DUMMY,
     &                                 ISYDEN)
                     ENDIF
                  ENDIF
!
!------------------------------------------------------------------
!              Calculate an AO-Fock matrix with C1. trans. density.
!------------------------------------------------------------------
!
                  ISYDEN  = ISYMTR
                  DTIME   = SECOND()
                  CALL CC_AOFOCK3(WORK(KXINT),WORK(KDNSC1),WORK(KFCKC1),
     &                         WORK(KEND2),LWRK2,IDEL,ISYMD,ISYDEN)
                  DTIME  = SECOND() - DTIME
                  TIMFCK = TIMFCK + DTIME
!
!-----------------------------------------
!              IF CCS jump to end of loop.
!-----------------------------------------
!
                  IF (.NOT. CCS) THEN
!
!------------------------------------------
!              Work space allocation no. 3.
!------------------------------------------
!
                     KDSRHF = KEND2
                     KEND3  = KDSRHF + NDSRHF(ISYMD)
                     LWRK3  = LWORK  - KEND3
!
                     IF (LWRK3 .LT. 0) THEN
                        WRITE(LUPRI,*) 'Need : ',KEND3,
     &                                 'Available : ',LWORK
                        CALL QUIT('Insufficient space in CC_RHTR3 (2)')
                     ENDIF
!
!--------------------------------------------------------
!              Transform one index in the integral batch.
!--------------------------------------------------------
!
                     DTIME   = SECOND()
                     ISYMLP  = 1
                     CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KLAMDP),
     &                                ISYMLP,WORK(KEND3),LWRK3,ISYDIS)
                     DTIME   = SECOND() - DTIME
                     TIMTRBT = TIMTRBT + DTIME
!
                     IF (IPRINT .GT. 200) THEN
                        CALL CC_PRINT(WORK(KXINT),WORK(KDSRHF),ISYDIS)
                     ENDIF
!
!--------------------------------------------------------------------
!              Construct the partially transformed C2-
!              amplitudes, and the partially transformed C2+
!              amplitudes.
!--------------------------------------------------------------------
!
                     KSCRCM  = KEND3
                     KSCRCM2 = KSCRCM + NT2MMO(ISYMD,ISYMTR)
                     KEND4   = KSCRCM2 + NT2MMO(ISYMD,ISYMTR)
                     LWRK4   = LWORK  - KEND4
!
                     DTIME   = SECOND()
                     ICON = 1
                     ISYMLH = 1
                     CALL CC_T2AO(WORK(KC2AM),WORK(KLAMDH),ISYMLH,
     *                            WORK(KSCRCM),WORK(KEND4),LWRK4,IDEL,
     *                            ISYMD,ISYMTR,ICON)
                     CALL CC_T2AO(WORK(KC22AM),WORK(KLAMDH),ISYMLH,
     *                            WORK(KSCRCM2),WORK(KEND4),LWRK4,IDEL,
     *                            ISYMD,ISYMTR,ICON)
                     DTIME   = SECOND() - DTIME
                     TIMT2AO = TIMT2AO + DTIME
!
                     IF (IPRINT .GT. 100) THEN
                        CALL AROUND('The partially transformed C2-')
                        CALL CC_PRTM(WORK(KSCRCM),ISYMD,ISYMTR)
                        CALL AROUND('The partially transformed C2+')
                        CALL CC_PRTM(WORK(KSCRCM2),ISYMD,ISYMTR)
                     ENDIF
!
!-------------------------------------
!              Calculate the B- term.
!-------------------------------------
!
                     IF (.NOT. CC2) THEN
!
                     IOPT = 5
                     DTIME   = SECOND()
                     ISYMM  = MULD2H(ISYMD,ISYMTR)
                     CALL CC_BF3(WORK(KXINT),WORK(KRHO2),WORK(KLAMDH),
     *                          1,WORK(KLHC1),ISYMTR,WORK(KLAMDH),1,
     *                          WORK(KSCRCM),ISYMM,DUMMY,1,
     *                          WORK(KEND4),LWRK4,IDEL,
     *                          ISYMD,IOPT)
                     DTIME   = SECOND() - DTIME
                     TIMBF   = TIMBF    + DTIME
!
                     ELSE
!
                     IOPT  = 2
                     ANTISYM = .TRUE.
                     DTIME = SECOND()
                     CALL CC_MOFCON3(WORK(KXINT),WORK(KRHO2),
     *                              WORK(KLAMDP),
     *                              WORK(KLAMDH),WORK(KLPC1),
     *                              WORK(KLHC1),WORK(KEND4),
     *                              LWRK4,IDEL,ISYMD,ISYMTR,IOPT,
     *                              ANTISYM)
                     DTIME = SECOND() - DTIME
                     TIMBF = TIMBF + DTIME
! 
                     ENDIF
!
!----------------------------------------
!              Calculate B+
!----------------------------------------
!
                     IF (.NOT. CC2) THEN
                     IOPT=6
                     DTIME   = SECOND()
                     ISYMM  = MULD2H(ISYMD,ISYMTR)
                     CALL CC_BF3(WORK(KXINT),WORK(KRHO22),WORK(KLAMDH),
     *                           1,WORK(KLHC1),ISYMTR,WORK(KLAMDH),1,
     *                           WORK(KSCRCM2),ISYMM,DUMMY,1,
     *                           WORK(KEND4),LWRK4,IDEL,
     *                           ISYMD,IOPT)
                     DTIME   = SECOND() - DTIME
                     TIMBF   = TIMBF    + DTIME
!
                     ELSE
                     IOPT = 2
                     ANTISYM = .FALSE.
                     DTIME = SECOND()
                     CALL CC_MOFCON3(WORK(KXINT),WORK(KRHO22),
     *                               WORK(KLAMDP),
     *                               WORK(KLAMDH),WORK(KLPC1),
     *                               WORK(KLHC1),WORK(KEND4),
     *                               LWRK4,IDEL,ISYMD,ISYMTR,IOPT,
     *                               ANTISYM)
                     DTIME = SECOND() - DTIME
                     TIMBF = TIMBF + DTIME
!
                     ENDIF
!
!---------------------------------------------------------------
!              Calculate the (1)D-tilde local intermediate.
!              Note: No c2am contribution for D-tilde.
!---------------------------------------------------------------
!
                     IF (.NOT. CC2) THEN
                     ISYMPC = ISYMTR
                     ISYMHC = ISYMTR
                     ICON   = 1  
                     FACTD  = 0.0D00
!
                     DTIME   = SECOND()
                     CALL CCRHS_D3(WORK(KXINT),WORK(KDSRHF),WORK(KRHO2),
     *                             WORK(KC2AM),ISYMTR,
     *                             WORK(KLAMDP),DUMMY,WORK(KLAMDH),
     *                             WORK(KLPC1),ISYMPC,WORK(KLHC1),
     *                             ISYMHC,WORK(KSCRCM),WORK(KEND4),
     *                             LWRK4,IDEL,ISYMD,FACTD,ICON,LUD,
     *                             DTFIL,IV)
!
                     DTIME   = SECOND() - DTIME
                     TIMD    = TIMD     + DTIME
!
                     ENDIF
!
!-----------------------------------------------------------
!              Calculate the (1)C-tilde local intermediate.
!              Note: No c2am contribution for C-tilde.
!-----------------------------------------------------------
!
                     IF (.NOT. CC2) THEN
                     FACTC = 0.0D00
                     ICON = 1
                     ISYMPC = ISYMTR
                     ISYMHC = ISYMTR
!
                     DTIME   = SECOND()
                     CALL CCRHS_C3(WORK(KXINT),WORK(KDSRHF),WORK(KRHO2),
     *                            WORK(KC2AM),ISYMTR,WORK(KLAMDP),DUMMY,
     *                            WORK(KLAMDH),WORK(KLPC1),ISYMPC,
     *                            WORK(KLHC1),ISYMHC,
     *                            WORK(KSCRCM),WORK(KEND4),
     *                            LWRK4,IDEL,ISYMD,FACTC,ICON,LUCS,
     *                            CSFIL,IV)
!
                     DTIME   = SECOND() - DTIME
                     TIMC    = TIMC     + DTIME
!
                     ENDIF
!
!--------------------------------------------------
!               Calculate (3)D-tilde.
!               Note: No C2 cont.
!--------------------------------------------------
!
                     IF (.NOT. CC2) THEN
                     ISYMPC = ISYMTR
                     ISYMHC = ISYMTR
                     ICON   = 4
                     FACTD  = 0.0D00
!
                     DTIME   = SECOND()
                     CALL CCRHS_D3(WORK(KXINT),WORK(KDSRHF),
     *                             WORK(KRHO22),WORK(KC22AM),ISYMTR,
     *                             WORK(KLAMDP),DUMMY,WORK(KLAMDH),
     *                             WORK(KLPC1),ISYMPC,WORK(KLHC1),
     *                             ISYMHC,WORK(KSCRCM2),WORK(KEND4),
     *                             LWRK4,IDEL,ISYMD,FACTD,ICON,LUP,
     *                             PFIL,IV)
!
                     DTIME   = SECOND() - DTIME
                     TIMD    = TIMD     + DTIME
!
                     ENDIF
!
!-----------------------------------------------------------
!              Calculate the (3)C-tilde local intermediate.
!              Note: No c2am contribution for C-tilde.
!-----------------------------------------------------------
!
                     IF (.NOT. CC2) THEN
                     FACTC = 0.0D00
                     ICON = 4 
!
                     DTIME   = SECOND()
                     CALL CCRHS_C3(WORK(KXINT),WORK(KDSRHF),WORK(KRHO2),
     *                            WORK(KC2AM),ISYMTR,WORK(KLAMDP),DUMMY,
     *                            WORK(KLAMDH),WORK(KLPC1),ISYMTR,
     *                            WORK(KLHC1),ISYMTR,
     *                            WORK(KSCRCM),WORK(KEND4),
     *                            LWRK4,IDEL,ISYMD,FACTC,ICON,LUC,
     *                            CTFIL,IV)
!
                     DTIME   = SECOND() - DTIME
                     TIMC    = TIMC     + DTIME
!
                     ENDIF
!
!-----------------------------------------------
!            Calculate E(tilde)-intermediates.
!            First the cont. from C2-
!            No trans to 2(C2-) -(C2-) so have 
!            g instead of L
!-----------------------------------------------
!
                     IF (.NOT. CC2) THEN
                     FACTE1 = XMTWO
                     FACTE2 = XMTWO
                     DTIME   = SECOND()
                     CALL CCRHS3_EI(WORK(KDSRHF),WORK(KEI1),
     *                              WORK(KEI2),WORK(KC2AM),
     *                              WORK(KSCRCM),WORK(KLAMDP),
     *                              WORK(KLAMDH),WORK(KEND4),LWRK4,
     *                              IDEL,ISYMD,ISYDIS,ISYMTR,
     *                              FACTE1,FACTE2)
!
                     DTIME   = SECOND() - DTIME
                     TIMEI   = TIMEI    + DTIME
!
                     ENDIF
!
!----------------------------------------------
!         Calculate the term from (+)C2
!         to E(tilde)-intermediates.
!----------------------------------------------
!
                     IF (.NOT. CC2) THEN
                     DTIME   = SECOND()
!
                     FACTE1  = TWO
                     FACTE2  = TWO
!
                     CALL CCRHS3_EI(WORK(KDSRHF),WORK(KEI1),
     *                             WORK(KEI2),WORK(KC22AM),
     *                             WORK(KSCRCM2),WORK(KLAMDP),
     *                             WORK(KLAMDH),WORK(KEND4),LWRK4,
     *                             IDEL,ISYMD,ISYDIS,ISYMTR,
     *                             FACTE1,FACTE2)
!
                     DTIME   = SECOND() - DTIME
                     TIMEI   = TIMEI    + DTIME
!
                     ENDIF
!
!---------------------------------------------
!         Calculate the (G-)-term.
!         No 2(C2-)-(C2-) so g instead of L
!---------------------------------------------
!   
                        ISYMP1 = 1
                        ISYMH1 = 1
                        FACG = TWO
                        DTIME   = SECOND()
                        CALL CCRHS_G3(WORK(KDSRHF),WORK(KOFF1),
     *                                WORK(KLAMDP),ISYMP1,WORK(KLAMDH),
     *                                ISYMH1,WORK(KSCRCM),WORK(KEND4),
     *                                LWRK4,ISYDIS,ISYMD,ISYMTR,FACG)
                        DTIME   = SECOND() - DTIME
                        TIMG    = TIMG     + DTIME
!
!----------------------------------------------
!          Calculate the (G+)-term
!----------------------------------------------
!
                        ISYMP1 = 1
                        ISYMH1 = 1
                        FACG = TWO
                        DTIME   = SECOND()
                        CALL CCRHS_G3(WORK(KDSRHF),WORK(KOFF1),
     *                                WORK(KLAMDP),ISYMP1,WORK(KLAMDH),
     *                                ISYMH1,WORK(KSCRCM2),WORK(KEND4),
     *                                LWRK4,ISYDIS,ISYMD,ISYMTR,FACG)
                        DTIME   = SECOND() - DTIME
                        TIMG    = TIMG     + DTIME
!
!---------------------------------------------
!         Calculate the (H-)-term.
!         No 2(C2-) - (C2-) so g instead of L
!---------------------------------------------
!
                        FACH = TWO
                        DTIME   = SECOND()
                        CALL CCRHS_H3(WORK(KDSRHF),WORK(KOFF1),
     *                                WORK(KLAMDP),WORK(KLAMDH),
     *                                WORK(KSCRCM),WORK(KEND4),LWRK4,
     *                                ISYDIS,ISYMD,ISYMTR,FACH)
                        DTIME   = SECOND() - DTIME
                        TIMH    = TIMH     + DTIME
!
!-------------------------------------------------
!          Calculate the (H+)-term.
!-------------------------------------------------
!
                        FACH = TWO
                        DTIME   = SECOND()
                        CALL CCRHS_H3(WORK(KDSRHF),WORK(KOFF1),
     *                                WORK(KLAMDP),WORK(KLAMDH),
     *                                WORK(KSCRCM2),WORK(KEND4),LWRK4,
     *                                ISYDIS,ISYMD,ISYMTR,FACH)
                        DTIME   = SECOND() - DTIME
                        TIMH    = TIMH     + DTIME
!
!-----------------------------------------------
!            (.NOT. CCS) loop ends here.
!-----------------------------------------------
!
                  ENDIF
!
!---------------------------------------------
!           Norms at the end of the loop.
!---------------------------------------------
!
                  IF (IPRINT .GT. 150) THEN
                     RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,
     *                            WORK(KOFF1),1)
                     WRITE(LUPRI,*)
     &                      'Norm of Rho1 end of loop:       ',RHO1N
                     IF (.NOT. CCS) THEN
                         RHO2N = DDOT(NRHO2M,WORK(KRHO2),1,
     *                            WORK(KRHO2),1)
                         WRITE(LUPRI,*)
     &                          'Norm BF- (end of loop)   :       ',
     *                               RHO2N
                         RHO2N = DDOT(NRHO2P,WORK(KRHO22),1,
     *                                WORK(KRHO22),1)
                         WRITE(LUPRI,*)
     &                          'Norm of BF+ (end of loop) :      ',
     *                               RHO2N
                     ENDIF
!
                  ENDIF
!
  125 	       CONTINUE
!
!--------------------------------------
!              write out result vector.
!--------------------------------------
!
               IF (.NOT. CCS) THEN
!
                  DTIME   = SECOND()
                  CALL CC_WVEC3(LUB3,BFFIL,NRHO2TOT,NRHO2M,
     *                          IV + ITR - 1,NRHO2P,WORK(KRHO2))
                  CALL CC_WVEC3(LUB3,BFFIL,NRHO2TOT,NRHO2P,
     *                          IV + ITR - 1,0,WORK(KRHO22))
!
                  DTIME   = SECOND() - DTIME
                  TIMIO   = TIMIO    + DTIME
!
               ENDIF
!
  120       CONTINUE
!
  110    CONTINUE
!
  100 CONTINUE
!
!==========================================
!
!     Calculation of RHO2-
!
!==========================================
!
!------------------------
!     Recover work space.
!------------------------
!
      KEND1 = KENDS2
      LWRK1 = LWRKS2
!
!--------------------------------------------
!     Allocate space for the gamma matrix.
!--------------------------------------------
!
      IF ((NEWGAM) .AND. (.NOT. (CCS .OR. CC2))) THEN
!
         KGAMMC = KEND1
         KEND1  = KGAMMC + NGAMMA(ISYMTR)*NSIMTR
         LWRK1  = LWORK  - KEND1
!
         IF (LWRK1 .LT. 0)
     &           CALL QUIT('Insufficient memory in cclr_rhtr3')
!
!--------------------------------------------
!     Initialize the gamma vector
!--------------------------------------------
!
         CALL DZERO(WORK(KGAMMC),NGAMMA(ISYMTR)*NSIMTR)
!
      END IF
!
!----------------------------------
!     Usual Fock Matrix.
!     Save AO Fock matrix in FOCKT
!----------------------------------
!
      ISYFAO = 1
      ISYMPA = 1
      ISYMHO = 1
!
      IF ( .NOT. CCS) THEN
         LFOCK = -1
         CALL GPOPEN(LFOCK,'CC_FCKH','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LFOCK)
         READ(LFOCK) (WORK(KFOCK0 + I - 1) , I = 1,N2BST(ISYMOP))
         CALL GPCLOSE(LFOCK,'KEEP')
      ENDIF
!
      IF (IPRINT .GT.140) THEN
         CALL AROUND( 'Usual Fock AO matrix' )
         CALL CC_PRFCKAO(WORK(KFOCK0),ISYFAO)
      ENDIF
!
      CALL DCOPY(N2BST(ISYMOP),WORK(KFOCK0),1,WORK(KFOCKT),1)
!
      DTIME = SECOND()
      CALL CC_FCKMO(WORK(KFOCK0),WORK(KLAMDP),WORK(KLAMDH),
     &              WORK(KEND1),LWRK1,ISYFAO,ISYMPA,ISYMHO)
      DTIME = SECOND() - DTIME
      TIMFCKMO = DTIME + TIMFCKMO
!
      IF (IPRINT .GT. 50) THEN
         CALL AROUND( 'Usual Fock MO matrix' )
         CALL CC_PRFCKMO(WORK(KFOCK0),ISYFAO)
      ENDIF
!
!
!----------------------------------------------------
!     Loop over number of trial vectors nsimtr.
!----------------------------------------------------
!
      DO 150 IV = 1, NSIMTR
!
         KOFF1  = KRHO1   + NT1AM(ISYMTR) *(IV-1)
         KAMC1  = KC1AM   + NT1AM(ISYMTR) *(IV -1)
         KLPC1  = KLAMPC  + NGLMDT(ISYMTR)*(IV-1)
         KLHC1  = KLAMHC  + NGLMDT(ISYMTR)*(IV-1)
         KFCKC1 = KFOCKC  + N2BST(ISYMTR) *(IV-1)
         KFCKC2 = KFOCKC2 + N2BST(ISYMTR) *(IV-1)
      IF (.NOT. (CCS)) THEN
         KGAMC1 = KGAMMC  + NGAMMA(ISYMTR)*(IV-1)
!
!===============================================================
!     After Loop CCSD Contributions with C2SQ in memory.
!===============================================================
!
!---------------------------------
!           Read in result vector.
!---------------------------------
!
         IF (.NOT. MINSCR) THEN
            DTIME = SECOND()
            CALL CC_RVEC3(LUB3,BFFIL,NRHO2TOT,NRHO2M,
     *                    IV + ITR -1,NRHO2P,WORK(KRHO2))
            DTIME = SECOND() - DTIME
            TIMIO = TIMIO + DTIME
         ENDIF
!
!--------------------------------------------------------------
!        Transform the Omega2- vector to the MO basis.
!--------------------------------------------------------------
!
         IF (.NOT. CC2) THEN
!
         ICON    = 1
         ISYMPC  = 1
         ISYMBF  = ISYMTR
         ANTISYM = .TRUE.
         KHLOG = .FALSE.
!
!--------------------------------------------------------
!        Calculate the (-)RHO-cont to Omega2-
!        Scale by two since calculate half the int.
!        that is needed.
!--------------------------------------------------------
!
         CALL DSCAL(NRHO2M,TWO,WORK(KRHO2),1)
!
         DTIME = SECOND()
         CALL CC_T2MOTRIP(FAKE,PHONEY,ISYMOP,WORK(KRHO2),WORK(KC2AM),
     *                 WORK(KGAMC1),WORK(KLAMDP),WORK(KLAMDP),ISYMPC,
     *                 WORK(KEND1),LWRK1,ISYMBF,ICON,DUMMY,
     *                 KHLOG,ANTISYM)
         DTIME = SECOND() -DTIME
         TIMT2MO = TIMT2MO + DTIME
!
         CALL DCOPY(NT2AMA(ISYMTR),WORK(KC2AM),1,WORK(KRHO2),1)
!
         IF ( LDEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
            RHO2N = DDOT(NT2AMA(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
            WRITE(LUPRI,*) 'Norm of Rho1  after loop in MO:    ',RHO1N
            WRITE(LUPRI,*) 'Norm of Rho2- after loop in MO:    ',RHO2N
         ENDIF
!
         IF ( LDEBUG ) THEN
            XNGAM = DDOT(NGAMMA(ISYMTR),WORK(KGAMC1),1,WORK(KGAMC1),1)
            WRITE(LUPRI,*) 'Norm of GamC :                  ',XNGAM
         ENDIF
!
         ENDIF
!
         IF (IPRINT.GT.120 .OR. LDEBUG) THEN
            CALL AROUND('After  T2MO-1:BF3(C1,C2) ')
            CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,1,1)
         ENDIF
!
!--------------------------------
!        write out result vector.
!--------------------------------
!
         DTIME   = SECOND()
         CALL CC_WVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                 NT2AMA(ISYMTR),IV + ITR -1,NT2AM(ISYMTR),
     *                 WORK(KRHO2))
         DTIME   = SECOND() - DTIME
         TIMIO   = TIMIO    + DTIME
!
      ENDIF
!
!---------------------------------------------
!        Continue the calculation of RHO2-
!---------------------------------------------
!
      IF (.NOT. CCS) THEN
!
!-------------------------------------------
!        Readin C2- amplitude in RHO2.
!-------------------------------------------
!
         DTIME = SECOND()
         CALL CC_RVEC3(LUFC2,FC2AM,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                 NT2AMA(ISYMTR),IVEC+IV-1,NT2AM(ISYMTR),
     *                 WORK(KRHO2))
         DTIME = SECOND() - DTIME
         TIMIO = TIMIO + DTIME
!
!--------------------------------
!        Square up C2- amplitudes.
!--------------------------------
!
         CALL CCLR_DIASCL(WORK(KRHO2),ZERO,ISYMTR)
         CALL CC_T2SQ3(WORK(KRHO2),WORK(KC2AM),ISYMTR)
!
         IF (IPRINT.GT.50) THEN
            CALL AROUND('CC_RHTR3 : (C1,C2-) vector readin')
            CALL CC_PRSQ(WORK(KAMC1),WORK(KC2AM),ISYMTR,1,1)
         ENDIF
!
!---------------------------------
!           Read in result vector.
!---------------------------------
!
         DTIME = SECOND()
         CALL CC_RVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                 NT2AMA(ISYMTR),IV + ITR -1,NT2AM(ISYMTR),
     *                 WORK(KRHO2))
         CALL CC_RVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                 NT2AM(ISYMTR),IV + ITR -1,0,
     *                 WORK(KRHO22))
         DTIME = SECOND() - DTIME
         TIMIO = TIMIO + DTIME
!
      ENDIF
!
      IF ( LDEBUG ) THEN
         RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
         WRITE(LUPRI,*) 'Norm of Rho1 after loop in mo:     ',RHO1N
         IF ( .NOT. CCS) THEN
            RHO2N = DDOT(NT2AMA(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
            WRITE(LUPRI,*) 'Norm of Rho2- after loop in mo:     ',RHO2N
         ENDIF
      ENDIF
!
!------------------------------------------
!     Transform AO Fock matrix to MO basis.
!------------------------------------------
!
!-------------------------------------
!     C1 transformed Fock Matrix.
!     F-star: F*ai = -SUM-k g(ak-bar,ki).
!-------------------------------------
!
      ISYFCK = ISYMTR
      ISYMPA = 1
      ISYMHO = 1
!
      IF (IPRINT .GT.140) THEN
         CALL AROUND( 'Fock AO matrix calc. from C1 transf. dens.' )
         CALL CC_PRFCKAO(WORK(KFCKC1),ISYFCK)
      ENDIF
!
      DTIME = SECOND()
      CALL CC_FCKMO(WORK(KFCKC1),WORK(KLAMDP),WORK(KLAMDH),
     *           WORK(KEND1),LWRK1,ISYFCK,ISYMPA,ISYMHO)
      DTIME = SECOND() - DTIME
      TIMFCKMO = DTIME + TIMFCKMO
!
      IF (IPRINT .GT. 50) THEN
         CALL AROUND( 'Fock MO matrix calc. from C1 transf. dens.' )
         CALL CC_PRFCKMO(WORK(KFCKC1),ISYFCK)
      ENDIF
!
!----------------------------------------------------------
!     Add -(SUM-k(g(ak-bar,ki)))
!----------------------------------------------------------
!
      CALL CCLR_1C1F(WORK(KOFF1),WORK(KFCKC1),ISYMTR)
!
      IF (IPRINT .GT. 120) THEN
         CALL AROUND('After  1C1F:  CCLR:RHO ')
         CALL CC_PRP(WORK(KOFF1),WORK,ISYMTR,1,0)
      ENDIF
!
      IF (LDEBUG) THEN
         RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
         WRITE(LUPRI,*) 'Norm of Rho1 after -g(ak-bar,ki) :',RHO1N
      ENDIF
!
!---------------------------------------------
!     Save -sum(k)g(ak-bar,ki)
!     Need this in calculation of rho1
!     KFCKC1 is used for F-tilde.
!---------------------------------------------
!
      IF (.NOT. CCS) THEN
        CALL DCOPY(N2BST(ISYMTR),WORK(KFCKC1),1,WORK(KFCKC2),1)
        CALL DSCAL(N2BST(ISYMTR),XMONE,WORK(KFCKC2),1)
      ENDIF
!
!---------------------------------------------------------------------
!     Make Fock tilde in the triplet case.
!     F(tilde)(ai) = -sum(k)g(ak-bar k i) + F(a-bar i) + F(a i-bar)
!---------------------------------------------------------------------
!
      IF ( .NOT. CCS ) THEN
C
         KFCKAO = KEND1
         KEND2  = KFCKAO + MAX(N2BST(ISYMOP),N2BST(ISYMTR))
         LWRK2  = LWORK  - KEND2
C
         IF (IPRINT .GT.140) THEN
            CALL AROUND( 'Fock tilde -1.' )
            CALL CC_PRFCKMO(WORK(KFCKC1),ISYMTR)
         ENDIF
C
         ISYFAO = 1
         ISYMPA = ISYMTR
         ISYMHO = 1
C
         CALL DCOPY(N2BST(ISYMOP),WORK(KFOCKT),1,WORK(KFCKAO),1)
         DTIME  = SECOND()
         CALL CC_FCKMO(WORK(KFCKAO),WORK(KLPC1),WORK(KLAMDH),
     *                 WORK(KEND2),LWRK2,ISYFAO,ISYMPA,ISYMHO)
         DTIME  = SECOND() - DTIME
         TIMFCKMO = DTIME  +  TIMFCKMO
C
         CALL DAXPY(N2BST(ISYMTR),ONE,WORK(KFCKAO),1,WORK(KFCKC1),1)
C
         IF (IPRINT .GT.140) THEN
            CALL AROUND( 'Fock tilde - 2.' )
            CALL CC_PRFCKMO(WORK(KFCKC1),ISYMTR)
         ENDIF
C
         ISYFAO = 1
         ISYMPA = 1
         ISYMHO = ISYMTR
C
         CALL DCOPY(N2BST(ISYMOP),WORK(KFOCKT),1,WORK(KFCKAO),1)
         DTIME  = SECOND()
         CALL CC_FCKMO(WORK(KFCKAO),WORK(KLAMDP),WORK(KLHC1),
     *                 WORK(KEND2),LWRK2,ISYFAO,ISYMPA,ISYMHO)
         DTIME  = SECOND() - DTIME
         TIMFCKMO = DTIME  +  TIMFCKMO
C
         CALL DAXPY(N2BST(ISYMTR),ONE,WORK(KFCKAO),1,WORK(KFCKC1),1)
C
         IF (IPRINT .GT. 50) THEN
            CALL AROUND( 'Fock-Tilde MO matrix ' )
            CALL CC_PRFCKMO(WORK(KFCKC1),ISYMTR)
         ENDIF
C
      ENDIF
!
!---------------------------------------------------
!     Readin E-intermediates: EI1(t1,t2),EI2(t1,t2).
!---------------------------------------------------
!
      KEI1   = KEND1
      KEI2   = KEI1 + NMATAB(ISYMOP)
      KEND2  = KEI2 + NMATIJ(ISYMOP)
      LWRK2  = LWORK - KEND2
!
      IF ( RSPIM.AND.((.NOT.(ONLY21)).AND.(.NOT.CCS))) THEN
!
         DTIME = SECOND() - DTIME
!
         LUE1 = -1
         CALL GPOPEN(LUE1,'CC_E1IM','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LUE1)
         READ (LUE1)(WORK(KEI1 + J-1),J = 1,NMATAB(ISYMOP))
         CALL GPCLOSE(LUE1,'KEEP' )
!
         LUE2 = -1
         CALL GPOPEN(LUE2,'CC_E2IM','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LUE2)
         READ (LUE2)(WORK(KEI2 + J-1),J = 1,NMATIJ(ISYMOP))
         CALL GPCLOSE(LUE2,'KEEP' )
!
         DTIME = SECOND() - DTIME
         TIMIO = TIMIO + DTIME
!
         IF (IPRINT.GT.150) THEN
            CALL AROUND( 'E-intermediates read from disk ')
            CALL CC_PREI(WORK(KEI1),WORK(KEI2),ISYMOP,1)
         ENDIF
!
      ENDIF
!
!---------------------------------------------------------
!        If ccs then put fock matrix into the EI matrices.
!---------------------------------------------------------
!
      IF (CCS) THEN
!
         FCKCON = .TRUE.
         ETRAN  = .FALSE.
         ISYMEI = 1
         DTIME = SECOND()
         CALL CCRHS_EFCK(WORK(KEI1),WORK(KEI2),WORK(KLAMDH),
     &                   WORK(KFOCK0),WORK(KEND2),LWRK2,FCKCON,
     &                   ETRAN,ISYMEI)
         DTIME  = SECOND() - DTIME
         TIMFCK = TIMFCK + DTIME
!
         IF (IPRINT.GT.150) THEN
            CALL AROUND( 'E-intermediates from EFCK ')
            CALL CC_PREI(WORK(KEI1),WORK(KEI2),ISYMOP,1)
         ENDIF
!
      ENDIF
!
!-----------------------------------------------------
!     Calculate E contributions to rho1: C1*EI(t1,t2).
!-----------------------------------------------------
!
      DTIME = SECOND()
      ISYMIM = 1
      CALL CCLR_E1C1(WORK(KOFF1),WORK(KAMC1),WORK(KEI1),WORK(KEI2),
     &               WORK(KEND2),LWRK2,ISYMTR,ISYMIM,'N')
!
      DTIME = SECOND() - DTIME
      TIME1C1 = TIME1C1 + DTIME
 
      IF (IPRINT .GT. 120) THEN
         CALL AROUND(' CC_RHTR3:  AFTER E1C1 EI*(C1) RHO1 ')
         CALL CC_PRP(WORK(KOFF1),WORK,ISYMTR,1,0)
      ENDIF
!
      IF ( LDEBUG ) THEN
         RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
         WRITE(LUPRI,*) 'Norm of Rho1 -after Ei*(c1,c2):    ',RHO1N
      ENDIF
!
!---------------------------------------------------------------
!     Calculate E contributions to rho2-: (C2-)*EI(t1,t2).
!     Scale by 2 to get the right sign/size.
!---------------------------------------------------------------
!
      IF (.NOT.(CCS)) THEN
!
        IF (.NOT. CC2) THEN
         CALL DSCAL(NMATAB(ISYMOP),TWO,WORK(KEI1),1)
         CALL DSCAL(NMATIJ(ISYMOP),TWO,WORK(KEI2),1)
         ANTISYM = .TRUE.
         KHLOG   = .FALSE.
         ISYVEC = ISYMTR
         ISYMIM = ISYMOP
         DTIME = SECOND()
         CALL CCRHS_E3(DUMMY,KHLOG,WORK(KC2AM),WORK(KEI1),WORK(KEI2),
     *                 WORK(KEND2),LWRK2,ISYVEC,ISYMIM,
     *                 WORK(KRHO2),ANTISYM)
         DTIME = SECOND() - DTIME
         TIMEI  = TIMEI + DTIME
!
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' CC_RHTR3:  AFTER EI*(C2) RHO2 ')
            CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,0,NC2)
         ENDIF
!
         IF ( LDEBUG ) THEN
           RHO2N = DDOT(NT2AMA(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
           WRITE(LUPRI,*) 'Norm of Rho2- after Ei*(c1,c2-):    ',RHO2N
         ENDIF
!
        ELSE
        DTIME = SECOND()
        ISIDE = 1
        IF (ISYMTR .EQ. 1) THEN
        CALL DSCAL(NC2AM,XMTWO,WORK(KC2AM),1)
        ELSE
        CALL DSCAL(NC2AM,TWO,WORK(KC2AM),1)
        ENDIF
        CALL CC2_FCK(WORK(KRHO2),WORK(KC2AM),WORK(KEND2),LWRK2,
     *               ISYMTR,WORK(KLAMDP),WORK(KLAMDH),ISIDE)
        IF (ISYMTR .EQ. 1) THEN
        CALL DSCAL(NC2AM,XMHALF,WORK(KC2AM),1)
        ELSE
        CALL DSCAL(NC2AM,HALF,WORK(KC2AM),1)
        ENDIF
        DTIME = SECOND() - DTIME
        TIMBF = TIMBF + DTIME 
!
        ENDIF
      ENDIF
!
!----------------------------------------------------
!     Readin  Gamma(t1,t2) vector: CCLR intermediate.
!----------------------------------------------------
!
      IF ((RSPIM) .AND. (.NOT.(CCS .OR. CC2))) THEN
!
        DTIME = SECOND()
!
        KGAMI  = KEND1
        KEND2  = KGAMI + NGAMMA(ISYMOP)
        LWRK2  = LWORK - KEND2
        LUGAM  = -1
        CALL GPOPEN(LUGAM,'CC_GAMIM','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &              .FALSE.)
        REWIND(LUGAM)
        READ (LUGAM)(WORK(KGAMI+J-1),J = 1,NGAMMA(ISYMOP))
        CALL GPCLOSE(LUGAM,'KEEP')
!
        DTIME = SECOND() - DTIME
        TIMIO = TIMIO + DTIME
!
!-------------------------------------------------
!        Calculate A-term: (C2-)*GAMMA(t1,t2).
!        Multiply with 2 to get the right sign/size.
!-------------------------------------------------
!
        ISYGAM = 1
        ISYVEC = ISYMTR
        IOPT    = 1
!
        CALL DSCAL(NGAMMA(ISYMOP),TWO,WORK(KGAMI),1)
 
        DTIME  = SECOND()
        CALL CCRHS_A(WORK(KRHO2),WORK(KC2AM),WORK(KGAMI),WORK(KEND2),
     *               LWRK2,ISYGAM,ISYVEC,IOPT)
        DTIME  = SECOND() - DTIME
        TIMA   = DTIME  + TIMA
!
        IF (IPRINT .GT. 120) THEN
           CALL AROUND(' AFTER C2*GAMMA(t1,t2) A-term RHO is ')
           KOFF1 = KRHO1 + NT1AM(ISYMTR)*(IV - 1)
           CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,0,1)
        ENDIF
!
        IF ( LDEBUG ) THEN
           RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
           RHO2N = DDOT(NT2AMA(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
           WRITE(LUPRI,*) 'Norm of Rho1  after C2*GAM:        ',RHO1N
           WRITE(LUPRI,*) 'Norm of Rho2- after C2*GAM:        ',RHO2N
        ENDIF
!
      ENDIF
!
!---------------------------------------------------
!     Calculate the C(triplet)-term: (C2-)*(1)C
!---------------------------------------------------
!
      IF  (.NOT. (CCS .OR. CC2)) THEN
!
         ISYCIM = 1
!
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' CC_RHTR:  AFTER EI*(C1,C2) RHO= ')
            CALL CC_PRSQ(WORK(KOFF1),WORK(KRHO2),ISYMTR,1,NC2)
         ENDIF
!
         ANTISYM= .TRUE.
         NORMALCONT = .FALSE.
         NFACC = ZERO    ! Doesn't matter.
         PIJCONT = .TRUE.
         PFACC = XMTWO
         IOPT  = 1
         IVECNR = 1
         DTIME  = SECOND()
         CALL CCSD_T2TP(WORK(KC2AM),WORK(KEND2),LWRK2,ISYMTR)
         CALL CCRHS3_CIO(WORK(KRHO2),WORK(KC2AM),WORK(KLAMDH),
     *                   WORK(KEND2),LWRK2,ISYMTR,ISYCIM,LUCSIM,
     *                   CSIMFIL,IVECNR,IOPT,NFACC,NORMALCONT,
     *                   PFACC,PIJCONT,ANTISYM)
         CALL CCSD_T2TP(WORK(KC2AM),WORK(KEND2),LWRK2,ISYMTR)
         DTIME = SECOND() - DTIME
         TIMCIO  = DTIME + TIMCIO
!
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' AFTER C2*CIM(t1,t2) C-term RHO is ')
            CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,0,1)
         ENDIF
!
         IF ( LDEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
            RHO2N = DDOT(NT2AMA(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
            WRITE(LUPRI,*) 'Norm of Rho1  after C2*CIM:        ',RHO1N
            WRITE(LUPRI,*) 'Norm of Rho2- after C2*CIM:        ',RHO2N
         ENDIF
!
      ENDIF
!
!----------------------------------------------
!     Calculate the D-term: C2*(-(1)C)
!     Don't transform to 2C2-C2
!----------------------------------------------
!
      IF  (.NOT. (CCS .OR. CC2)) THEN
!
         ISYDIM = 1
!
         IOPT    = 1
         IVECNR  = 1
         FACTD   = ZERO
         FACTD2  = XMTWO
         ANTISYM = .TRUE.
         KHLOG   = .FALSE.
         DTIME   = SECOND()
         CALL CCRHS3_DIO(DUMMY,KHLOG,WORK(KC2AM),WORK(KLAMDH),
     *                   WORK(KEND2),LWRK2,ISYMTR,ISYDIM,LUCSIM,
     *                   CSIMFIL,IVECNR,IOPT,FACTD,ANTISYM,
     *                   WORK(KRHO2),FACTD2)
         DTIME = SECOND() - DTIME
         TIMDIO  = DTIME + TIMDIO
!
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' AFTER C2*DIM(t1,t2) D-term RHO is ')
            CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,0,1)
         ENDIF
!
         IF ( LDEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
            RHO2N = DDOT(NT2AMA(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
            WRITE(LUPRI,*) 'Norm of Rho1  after C2*DIM:        ',RHO1N
            WRITE(LUPRI,*) 'Norm of Rho2- after C2*DIM:        ',RHO2N
         ENDIF
!
      ENDIF
!
!----------------------------------------------------------------------
!     Readin  BF(t1,t2) and contract with lampc and lamdp  and add.
!     ANTISYM is false since the global rho(BF) is symmetric in (al i be j)
!     Calculate both the (+) and (-) terms 
!     NB This destroys C2AM.
!----------------------------------------------------------------------
!
      IF (.NOT. (CCS .OR. CC2)) THEN
!
         IF ( RSPIM ) THEN
!
            DTIME = SECOND()
!
            LUBF = -1
            CALL GPOPEN(LUBF,'CC_BFIM','OLD',' ','UNFORMATTED',IDUMMY,
     &                  .FALSE.)
            READ(LUBF) (WORK(KC2AM +J - 1),J = 1,2*NT2ORT(1))
            CALL GPCLOSE(LUBF,'KEEP')
!
            DTIME = SECOND() - DTIME
            TIMIO = TIMIO + DTIME
!
         ENDIF
!
         ISYMBF = 1
         NEWGS  = NEWGAM
         NEWGAM = .FALSE.
         ICON = 5
         ANTISYM= .FALSE.
         KHLOG  = .TRUE.
!
         DTIME = SECOND()
!
         KOFF2 = KC2AM + 2*NT2ORT(1)
         CALL CC_T2MOTRIP(FAKE,PHONEY,ISYMOP,WORK(KC2AM),WORK(KOFF2),
     *                    DUMMY,WORK(KLAMDP),WORK(KLPC1),ISYMTR,
     *                    WORK(KEND1),LWRK1,ISYMBF,ICON,WORK(KRHO22),
     *                    KHLOG,ANTISYM)
         CALL DAXPY(NT2AM(ISYMTR),ONE,WORK(KOFF2),1,WORK(KRHO2),1)
         DTIME = SECOND() -DTIME
         TIMT2MO = TIMT2MO + DTIME
!
         NEWGAM = NEWGS
!
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' AFTER C1*BFIM(t1,t2) BF-tilde term RHO is ')
            CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,0,1)
         ENDIF
!
         IF ( LDEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
            RHO2N = DDOT(NT2AMA(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
            WRITE(LUPRI,*) 'Norm of Rho1  after C1*BFIM:       ',RHO1N
            WRITE(LUPRI,*) 'Norm of Rho2- after C1*BFIM:       ',RHO2N
         ENDIF
!
      ENDIF
!
!
!-----------------------------
!     write out result vector.
!-----------------------------
!
      CALL CC_WVEC(LUFR1,FRHO1,NT1AM(ISYMTR),NT1AM(ISYMTR),
     &             IV + IVEC -1,WORK(KOFF1))
!
      IF (.NOT. CCS) THEN
         DTIME   = SECOND()
         CALL CC_WVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                 NT2AMA(ISYMTR),IV + ITR -1,NT2AM(ISYMTR),
     *                 WORK(KRHO2))
         CALL CC_WVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                 NT2AM(ISYMTR),IV + ITR -1,0,
     *                 WORK(KRHO22))
         DTIME   = SECOND() - DTIME
         TIMIO   = TIMIO    + DTIME
      ENDIF
!
150   CONTINUE
!
!=========================================================
!     MO basis (+)R + (-)R section.
!     Contract global intermediates with (+)R + (-)R
!=========================================================
!
      IF (.NOT. CCS) THEN
        DO 155 IV = 1, NSIMTR
!
          KOFF1 = KRHO1 + NT1AM(ISYMTR)*(IV-1)
!
!----------------------------------------------
!     Readin and square up the C2- amplitudes.
!----------------------------------------------
!
          DTIME = SECOND()
          CALL CC_RVEC3(LUFC2,FC2AM,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                  NT2AMA(ISYMTR),IVEC+IV-1,NT2AM(ISYMTR),
     *                  WORK(KRHO2))
          DTIME = SECOND() - DTIME
          TIMIO = TIMIO + DTIME
!
!
          CALL CCLR_DIASCL(WORK(KRHO2),ZERO,ISYMTR)
          CALL CC_T2SQ3(WORK(KRHO2),WORK(KC2AM),ISYMTR)
!
!-----------------------------------------------
!     Readin and square up the C2+ amplitudes.
!-----------------------------------------------
!
          DTIME = SECOND()
          CALL CC_RVEC3(LUFC2,FC2AM,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                  NT2AM(ISYMTR),IVEC+IV-1,0,
     *                  WORK(KRHO2))
          DTIME = SECOND() - DTIME
          TIMIO = TIMIO + DTIME
!
!
          CALL CCRHS3_R2IJ(WORK(KRHO2),WORK(KEND1),LWRK1,ISYMTR)
          CALL CC_T2SQ(WORK(KRHO2),WORK(KEND1),ISYMTR)
!
!-----------------------------------------------
!     Calculate the (C2-) + (C2+)
!-----------------------------------------------
!
          CALL DAXPY(NC2AM,ONE,WORK(KEND1),1,WORK(KC2AM),1)
!
!-----------------------------------------------
!     Readin the result vectors.
!     The (-) result vector in KRHO2
!     The (+) result vector in KRHO22
!-----------------------------------------------
!
          IF (.NOT. CC2) THEN
            DTIME = SECOND()
            CALL CC_RVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                    NT2AMA(ISYMTR),IV + ITR -1,NT2AM(ISYMTR),
     *                    WORK(KRHO2))
            CALL CC_RVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                    NT2AM(ISYMTR),IV + ITR -1,0,
     *                    WORK(KRHO22))
            DTIME = SECOND() - DTIME
            TIMIO = TIMIO + DTIME
          ENDIF
!
!---------------------------------------------
!      Calculate I-term to Rho1.
!---------------------------------------------
!
      IF (.NOT. CCS) THEN
         CALL DSCAL(N2BST(ISYMOP),TWO,WORK(KFOCK0),1)
!
         DTIME    = SECOND()
         CALL CCRHS_I(WORK(KOFF1),WORK(KC2AM),WORK(KFOCK0),WORK(KEND2),
     *                LWRK2,ISYMTR,1)
         CALL DSCAL(N2BST(ISYMOP),HALF,WORK(KFOCK0),1)
         DTIME    = SECOND()-DTIME
         TIMI     = DTIME + TIMI
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' AFTER I-term :RHO ')
            CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,1,0)
         ENDIF
         IF (LDEBUG) THEN
            RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
            RHO2N = DDOT(NT2AMA(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
            WRITE(LUPRI,*) 'Norm of Rho1  after I-term:    ',RHO1N
            WRITE(LUPRI,*) 'Norm of Rho2- after I-term:    ',RHO2N
         ENDIF
      ENDIF
!
!--------------------------------------------------
!      Calculate the (C2(+) + C2(-)) * (3)D term
!--------------------------------------------------
!
          IF (.NOT. CC2) THEN
!
            ISYDIM = 1
            IOPT    = 1
            IVECNR  = 1
            FACTD   = TWO
            FACTD2  = TWO
            ANTISYM = .TRUE.
            KHLOG   = .TRUE.
            DTIME   = SECOND()
            CALL CCRHS3_DIO(WORK(KRHO22),KHLOG,WORK(KC2AM),
     *                      WORK(KLAMDH),WORK(KEND2),LWRK2,
     *                      ISYMTR,ISYDIM,LUDTIM,DTIMFIL,IVECNR,
     *                      IOPT,FACTD,ANTISYM,WORK(KRHO2),FACTD2)
            DTIME = SECOND() - DTIME
            TIMDIO  = DTIME + TIMDIO
!
!
!-----------------------------------------------------
!      Write out the result vectors.
!-----------------------------------------------------
!
            DTIME   = SECOND()
            CALL CC_WVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                    NT2AMA(ISYMTR),IV + ITR -1,NT2AM(ISYMTR),
     *                    WORK(KRHO2))
            CALL CC_WVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                    NT2AM(ISYMTR),IV + ITR -1,0,WORK(KRHO22))
!
          ENDIF
!
          DTIME   = SECOND() - DTIME
          TIMIO   = TIMIO    + DTIME 
!
  155   CONTINUE
      ENDIF
!
!==========================================================
!     MO basis T2SQ section.
!     Contract intermediates constructed in loop with T2SQ.
!==========================================================
!
      IF (.NOT. CCS) THEN 
!
!
!----------------------------------------------
!     Loop over number of trial vectors nsimtr.
!----------------------------------------------
!
      DO 160 IV = 1, NSIMTR
!
            KOFF1  = KRHO1   + NT1AM(ISYMTR) *(IV-1)
            KAMC1  = KC1AM   + NT1AM(ISYMTR) *(IV-1)
            KLPC1  = KLAMPC  + NGLMDT(ISYMTR)*(IV-1)
            KLHC1  = KLAMHC  + NGLMDT(ISYMTR)*(IV-1)
            KFCKC1 = KFOCKC  + N2BST(ISYMTR) *(IV-1)
            KFCKC2 = KFOCKC2 + N2BST(ISYMTR) *(IV-1)
            KEI1   = KEMAT1  + NEMAT1(ISYMTR)*(IV-1)
            KEI2   = KEMAT2  + NMATIJ(ISYMTR)*(IV-1)
            KGAMC1 = KGAMMC  + NGAMMA(ISYMTR)*(IV-1)
!
!------------------------------------
!        Readin T2 amplitude in RHO2.
!------------------------------------
!
            DTIME = SECOND()
!
            IOPT = 2
            CALL CC_RDRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KRHO2))
!
            DTIME = SECOND() - DTIME
            TIMIO = TIMIO + DTIME
!
!--------------------------------
!        Square up T2 amplitudes.
!--------------------------------
!
            DTIME = SECOND()
            CALL CC_T2SQ(WORK(KRHO2),WORK(KC2AM),1)
            DTIME = SECOND() - DTIME
            TIMT2SQ = TIMT2SQ + DTIME
!
            IF (IPRINT.GT.50) THEN
               CALL AROUND('CC_RHTR: (T1,T2) vector readin')
               CALL CC_PRSQ(WORK(KT1AM),WORK(KC2AM),1,1,1)
            ENDIF
!
!--------------------------------------
!        Read in result vectors.
!--------------------------------------
!
            DTIME = SECOND()
            CALL CC_RVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                    NT2AMA(ISYMTR),IV + ITR  -1,NT2AM(ISYMTR),
     *                    WORK(KRHO2))
            CALL CC_RVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                    NT2AM(ISYMTR),IV + ITR  -1,0,
     *                    WORK(KRHO22))
            DTIME = SECOND() - DTIME
            TIMIO = TIMIO + DTIME
!
!----------------------------------------------------
!     Calculate E-tilde cont. to RHO2-
!----------------------------------------------------
!
            IF (.NOT. CC2) THEN
!
            IF (IPRINT.GT.150) THEN
               CALL AROUND('E-intermdiates calc. EI(C1,C2+,C2-)-ao1')
               CALL CC_PREI(WORK(KEI1),WORK(KEI2),ISYMTR,0)
            ENDIF
!
            ETRAN  = .TRUE.
            FCKCON = .TRUE.
            ISYMEI = ISYMTR
            DTIME = SECOND()
!
            CALL CCRHS_EFCK(WORK(KEI1),WORK(KEI2),WORK(KLAMDH),
     *                      WORK(KFCKC1),WORK(KEND1),LWRK1,FCKCON,
     *                      ETRAN,ISYMTR)
!
            DTIME = SECOND() - DTIME
            TIMFCK = TIMFCK + DTIME
!
            IF (IPRINT.GT.150) THEN
               CALL AROUND( 'E-intermdiates calc EI(C1,C2+,C2-) -mo' )
               CALL CC_PREI(WORK(KEI1),WORK(KEI2),ISYMTR,1)
            ENDIF
!
            KHLOG   = .TRUE.
            ANTISYM = .TRUE.
            ISYVEC = 1
            ISYMIM = ISYMTR
            DTIME = SECOND()
            CALL CCRHS_E3(WORK(KRHO22),KHLOG,WORK(KC2AM),WORK(KEI1),
     *                    WORK(KEI2),WORK(KEND1),LWRK1,ISYVEC,
     *                    ISYMTR,WORK(KRHO2),ANTISYM)
            DTIME = SECOND() - DTIME
            TIMEI  = TIMEI + DTIME
!
            IF (IPRINT .GT. 120) THEN
               CALL AROUND(' AFTER CCRHS_E3 cont.  CCLR:RHO(-) ')
               CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,0,1)
            ENDIF
!
            IF ( LDEBUG ) THEN
               RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
               RHO2N = DDOT(NT2AMA(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
               WRITE(LUPRI,*) 'Norm of Rho1  after T2*(-)EI:         ',
     &                        RHO1N
               WRITE(LUPRI,*) 'Norm of Rho2- after T2*(-)EI:         ',
     &                        RHO2N
            ENDIF
            ENDIF
!
!-----------------------------------
!     Calculate T2*GAMMA(C2-) A-term.
!-----------------------------------
!
            IF (.NOT. CC2) THEN
            ISYGAM  = ISYMTR
            ISYVEC  = 1
            IOPT    = 1
!
            CALL DSCAL(NGAMMA(ISYMTR),XMONE,WORK(KGAMC1),1)
!
            DTIME   = SECOND()
            ANTISYM = .TRUE.
            CALL CCRHS_A3(WORK(KRHO2),WORK(KC2AM),WORK(KGAMC1),
     *                    WORK(KEND1),LWRK1,ISYGAM,ISYVEC,IOPT,
     *                    ANTISYM)
            DTIME   = SECOND() - DTIME
            TIMA    = DTIME + TIMA
!
            IF (IPRINT .GT. 120) THEN
               CALL AROUND(' AFTER T2*GAMMA(C1,C2) A-term RHO ' )
               CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,0,1)
            ENDIF
!
            IF ( LDEBUG ) THEN
               RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
               RHO2N = DDOT(NT2AMA(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
               WRITE(LUPRI,*) 'Norm of Rho1  after T2*GAM:        ',
     &                        RHO1N
               WRITE(LUPRI,*) 'Norm of Rho2- after T2*GAM:        ',
     &                        RHO2N
            ENDIF
            ENDIF
!
!------------------------------------
!     Recover work from gamma matrix.
!     IF not minscr!
!------------------------------------
!
            IF (MINSCR) THEN
               KEND1 = KENDS2
               LWRK1 = LWRKS2
            ENDIF 
!
!--------------------------------------------
!     Calculate the C-term: T2*(T)C(C1).
!--------------------------------------------
!
            IF (.NOT. CC2) THEN
            ISYCIM = ISYMTR
            ISYVEC = 1
!
            IOPT = 2
            NORMALCONT = .FALSE.
            NFACC = ZERO      ! Doesn't matter.
            PIJCONT = .TRUE.
            PFACC = ONE
            ANTISYM = .TRUE.
            DTIME  = SECOND()
            CALL CCSD_T2TP(WORK(KC2AM),WORK(KEND1),LWRK1,ISYVEC)
            CALL CCRHS3_CIO(WORK(KRHO2),WORK(KC2AM),WORK(KLAMDH),
     *                      WORK(KEND1),LWRK1,ISYVEC,ISYCIM,LUC,CTFIL,
     *                      IV,IOPT,NFACC,NORMALCONT,PFACC,PIJCONT,
     *                      ANTISYM)
            CALL CCSD_T2TP(WORK(KC2AM),WORK(KEND1),LWRK1,ISYVEC)
            DTIME   = SECOND() - DTIME
            TIMCIO  = DTIME + TIMCIO
!
            IF (IPRINT .GT. 120) THEN
               CALL AROUND(' AFTER T2*CIM(C1) C-term RHO is ')
               CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,0,1)
            ENDIF
!
            IF ( LDEBUG ) THEN
               RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
               RHO2N = DDOT(NT2AMA(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
               WRITE(LUPRI,*) 'Norm of Rho1  after T2*CI:         ',
     &                        RHO1N
               WRITE(LUPRI,*) 'Norm of Rho2- after T2*CI:         ',
     &                        RHO2N
            ENDIF
!
            ENDIF
!
!-------------------------------------------------------
!     Calculate (T2)*(1)D(C1).
!-------------------------------------------------------
!
            IF (.NOT. CC2) THEN
            ISYDIM  = ISYMTR
            ISYVEC  = 1
            IOPT    = 2
            ANTISYM = .TRUE.
            KHLOG   = .TRUE.
            FACTD   = ONE
            FACTD2  = XMONE
!
            DTIME   = SECOND()
            CALL CCRHS3_DIO(WORK(KRHO22),KHLOG,WORK(KC2AM),
     *                      WORK(KLAMDH),WORK(KEND1),LWRK1,
     *                      ISYVEC,ISYDIM,LUD,DTFIL,
     *                      IV,IOPT,FACTD,ANTISYM,
     *                      WORK(KRHO2),FACTD2)
            DTIME   = SECOND() - DTIME
            TIMDIO  = DTIME + TIMDIO
!
            IF (IPRINT .GT. 120) THEN
               CALL AROUND(' AFTER (2T2-T2)*DI(C1) D-term RHO is ')
               CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,0,1)
            ENDIF
!
            IF ( LDEBUG ) THEN
               RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
               RHO2N = DDOT(NT2AMA(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
               WRITE(LUPRI,*) 'Norm of Rho1  after T2*DI:         ',
     &                        RHO1N
               WRITE(LUPRI,*) 'Norm of Rho2- after T2*DI:         ',
     &                        RHO2N
            ENDIF
!
            ENDIF
!
!------------------------------------------
!     Transform the T2amp. to - exchange.
!------------------------------------------
!
            CALL CCSD_T2TP(WORK(KC2AM),WORK(KEND1),LWRK1,1)
!
!-------------------------------------------------------
!     Calculate (T2)*(3)D(C1).
!-------------------------------------------------------
!
            IF (.NOT. CC2) THEN
            ISYDIM = ISYMTR
            ISYVEC = 1
            IOPT = 2
            ANTISYM = .TRUE.
            KHLOG   = .TRUE.
            FACTD  = XMONE
            FACTD2 = ONE
!
            DTIME   = SECOND()
            CALL CCRHS3_DIO(WORK(KRHO22),KHLOG,WORK(KC2AM),
     *                      WORK(KLAMDH),WORK(KEND1),LWRK1,
     *                      ISYVEC,ISYDIM,LUP,PFIL,IV,
     *                      IOPT,FACTD,ANTISYM,WORK(KRHO2),FACTD2)
            DTIME   = SECOND() - DTIME
            TIMDIO  = DTIME + TIMDIO
!
            IF (IPRINT .GT. 120) THEN
               CALL AROUND(' AFTER (2T2-T2)*DI(C1) D-term RHO is ')
               CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,0,1)
            ENDIF
!
            IF ( LDEBUG ) THEN
               RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
               RHO2N = DDOT(NT2AMA(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
               WRITE(LUPRI,*) 'Norm of Rho1  after T2*DI:         ',
     &                        RHO1N
               WRITE(LUPRI,*) 'Norm of Rho2- after T2*DI:         ',
     &                        RHO2N
            ENDIF
!
            ENDIF
!
!----------------------------------------------------------
!     Calculate I-tilde-term. Have only the exchange part
!     in F-tilde. No Coulomb, F(a-bar,i) or F(a,i-bar)
!     Use the saved -sum(k)g(pk-bar,kq).
!     For T2 we also only have the exchange part.
!----------------------------------------------------------
!
            DTIME    = SECOND()
            CALL CCRHS_I(WORK(KOFF1),WORK(KC2AM),WORK(KFCKC2),
     *                   WORK(KEND1),LWRK1,1,ISYMTR)
            DTIME   = SECOND() - DTIME
            TIMI    = DTIME + TIMI
            IF (IPRINT .GT. 120) THEN
               CALL AROUND(' AFTER CCLR_I tilde contribution: RHO ')
               CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,1,0)
            ENDIF
!
            IF ( LDEBUG ) THEN
               RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
               RHO2N = DDOT(NT2AMA(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
               WRITE(LUPRI,*) 'Norm of Rho1  after I-tilde:   ',RHO1N
               WRITE(LUPRI,*) 'Norm of Rho2- after I-tilde:   ',RHO2N
            ENDIF
!
!---------------------------------------
!        Diagonal of RHO2- is zero.
!---------------------------------------
!
            CALL CCLR_DIASCL(WORK(KRHO2),ZERO,ISYMTR)
!
!====================================
!        Scale RHO2-
!====================================
!
            CALL DSCAL(NT2AMA(ISYMTR),HALF,WORK(KRHO2),1)
!
!------------------------------------
!        write out result vectors.
!------------------------------------
!
         DTIME   = SECOND()
         CALL CC_WVEC(LUFR1,FRHO1,NT1AM(ISYMTR),NT1AM(ISYMTR),
     *                IV+IVEC-1,WORK(KOFF1))
         IF (.NOT. CCS) THEN
            CALL CC_WVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                    NT2AMA(ISYMTR),IV + ITR -1,NT2AM(ISYMTR),
     *                    WORK(KRHO2))
            CALL CC_WVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                    NT2AM(ISYMTR),IV + ITR -1,0,
     *                    WORK(KRHO22))
            DTIME   = SECOND() - DTIME
            TIMIO   = TIMIO    + DTIME
         ENDIF
!
!=============================================================
!     End section for RHO2-.
!=============================================================
!
         IF (MINSCR) THEN
           KEND1 = KENDS2
           LWRK1 = LWRKS2
         END IF
!
         IF (IPRINT .GT. 50 ) THEN
            CALL AROUND('END OF FIRST PART OF CC_RHTR3 :RHO1, RHO2- ')
            WRITE (LUPRI,*) 'number of doubles vector on file:',
     &                      IV+ITR-1
            NC2 = 1
            IF ( CCS ) NC2 = 0
            CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYMTR,1,NC2)
!
         ENDIF
!
         IF ( LDEBUG ) THEN
             RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
             WRITE(LUPRI,*) 'Norm of Rho1 for IV = ',IV -1 + IVEC,
     *                  ' -end of cclr_rhtr3 (-part1):    ',RHO1N
             IF (.NOT. CCS) THEN
                 RHO2N = DDOT(NT2AMA(ISYMTR),WORK(KRHO2),1,
     *                        WORK(KRHO2),1)
                 WRITE(LUPRI,*) 
     *           'Norm of Rho2- (end of cclr_rhtr3 (part1)):    ',RHO2N
             ENDIF
         ENDIF
!
  160 CONTINUE
!
!++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     Calculation of RHO2+
!
!++++++++++++++++++++++++++++++++++++++++++++++++++
!
!--------------------------------------------
!     Allocate space for the gamma matrix.
!--------------------------------------------
!
         IF (NEWGAM .AND. (.NOT. CC2)) THEN
!
            KGAMMC = KEND1
            KEND1  = KGAMMC + NGAMMA(ISYMTR)*NSIMTR
            LWRK1  = LWORK  - KEND1
!
            IF (LWRK1 .LT. 0)
     &           CALL QUIT('Insufficient memory in cclr_rhtr3_1')
!
!--------------------------------------------
!     Initialize the gamma vector
!--------------------------------------------
!
            CALL DZERO(WORK(KGAMMC),NGAMMA(ISYMTR)*NSIMTR)
!
         ENDIF
!
!----------------------------------------------------
!     Loop over number of trial vectors nsimtr.
!----------------------------------------------------
!
         DO 170 IV = 1, NSIMTR
!
            KOFF1  = KRHO1  + NT1AM(ISYMTR)*(IV-1)
            KAMC1  = KC1AM  + NT1AM(ISYMTR)*(IV-1)
            KLPC1  = KLAMPC + NGLMDT(ISYMTR)*(IV-1)
            KLHC1  = KLAMHC + NGLMDT(ISYMTR)*(IV-1)
            KFCKC1 = KFOCKC + N2BST(ISYMTR)*(IV-1)
            KGAMC1 = KGAMMC + NGAMMA(ISYMTR) * (IV-1)
!
!===============================================================
!     After Loop CCSD Contributions with C2SQ in memory.
!===============================================================
!
!-------------------------------------------
!           Read in local BF-int. vector.
!-------------------------------------------
!
               DTIME = SECOND()
               CALL DZERO(WORK(KRHO2),NRHO2)
               CALL CC_RVEC3(LUB3,BFFIL,NRHO2TOT,NRHO2P,
     *                       IV + ITR -1,0,WORK(KRHO2))
!
               IF (.NOT. CC2) THEN
!
               DTIME = SECOND() - DTIME
               TIMIO = TIMIO + DTIME
!
!--------------------------------------------------------------
!        Transform the Omega2+ vector to the MO basis.
!--------------------------------------------------------------
!
            ICON    = 6
            ISYMPC  = 1
            ISYMBF  = ISYMTR
            ANTISYM = .FALSE.
            KHLOG   = .FALSE.
!
            DTIME = SECOND()
            CALL CC_T2MOTRIP(FAKE,PHONEY,ISYMOP,WORK(KRHO2),WORK(KC2AM),
     *                       WORK(KGAMC1),WORK(KLAMDP),WORK(KLAMDP),
     *                       ISYMPC,WORK(KEND1),LWRK1,ISYMBF,
     *                       ICON,DUMMY,KHLOG,ANTISYM)
!
            TIMT2MO = TIMT2MO + DTIME
!
!---------------------------------------------------------------
!         Readin the result vector and add the BF-tilde cont.
!---------------------------------------------------------------
!
            CALL CC_RVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                    NT2AM(ISYMTR),IV+ITR-1,0,
     *                    WORK(KRHO2))
!
            CALL DAXPY(NT2AM(ISYMTR),ONE,WORK(KC2AM),1,WORK(KRHO2),1)
!
            IF (LDEBUG) THEN
               RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
               RHO2N = DDOT(NT2AM(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
               WRITE(LUPRI,*) 'Norm of Rho1  after loop in MO:    ',
     &                        RHO1N
               WRITE(LUPRI,*) 'Norm of Rho2+ after loop in MO:    ',
     &                        RHO2N
            ENDIF
!
            IF (IPRINT .GT. 120) THEN
               CALL AROUND('After  T2MO-1:BF3(C1,C2+) ')
               CALL CC_PRP(WORK(KAMC1),WORK(KRHO2),ISYMTR,0,1)
            ENDIF
!
            IF ( LDEBUG ) THEN
               XNGAM = DDOT(NGAMMA(ISYMTR),WORK(KGAMC1),1,
     *                      WORK(KGAMC1),1)
               WRITE(LUPRI,*) 
     *               'Norm of GamC -after calc.:               ',XNGAM
            ENDIF
!
            ENDIF
!
!--------------------------------
!        write out result vector.
!--------------------------------
!
            DTIME   = SECOND()
            CALL CC_WVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                    NT2AM(ISYMTR),IV + ITR -1,0,WORK(KRHO2))
            DTIME   = SECOND() - DTIME
            TIMIO   = TIMIO    + DTIME
!
!---------------------------------------------
!        Continue the calculation of RHO2+
!---------------------------------------------
!
!-----------------------------------
!        Readin C2+ amplitude in RHO.
!-----------------------------------
!
            DTIME = SECOND()
            CALL CC_RVEC3(LUFC2,FC2AM,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                    NT2AM(ISYMTR),IVEC+IV-1,0,WORK(KRHO2))
            DTIME = SECOND() - DTIME
            TIMIO = TIMIO + DTIME
!
!--------------------------------
!        Square up C2+ amplitudes.
!--------------------------------
!
            CALL CCRHS3_R2IJ(WORK(KRHO2),WORK(KEND1),LWRK1,ISYMTR)
            CALL CCLR_DIASCL(WORK(KRHO2),ZERO,ISYMTR)
            CALL CC_T2SQ(WORK(KRHO2),WORK(KC2AM),ISYMTR)
!
            IF (IPRINT.GT.50) THEN
               CALL AROUND('CC_RHTR3: C2+ vector readin/sq')
               CALL CC_PRSQ(WORK(KAMC1),WORK(KC2AM),ISYMTR,0,1)
            ENDIF
!
!---------------------------------
!           Read in result vector.
!---------------------------------
!
               DTIME = SECOND()
               CALL CC_RVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                       NT2AM(ISYMTR),IV + ITR -1,0,WORK(KRHO2))
               DTIME = SECOND() - DTIME
               TIMIO = TIMIO + DTIME
!
            IF ( LDEBUG ) THEN
               RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
               WRITE(LUPRI,*) 'Norm of Rho1 after loop in mo:     ',
     &                        RHO1N
               IF ( .NOT. CCS) THEN
                  RHO2N = DDOT(NT2AM(ISYMTR),WORK(KRHO2),1,
     *                         WORK(KRHO2),1)
                  WRITE(LUPRI,*) 'Norm of Rho2+ after loop in mo:     ',
     *                        RHO2N
               ENDIF
            ENDIF
            IF (IPRINT .GT. 120) THEN
               CALL AROUND('After reading of file:')
               CALL CC_PRP(WORK(KAMC1),WORK(KRHO2),ISYMTR,0,1)
            ENDIF
!
!---------------------------------------------------
!     Readin E-intermediates: EI1(t1,t2),EI2(t1,t2).
!---------------------------------------------------
!
            IF (.NOT. CC2) THEN
!
            KEI1   = KEND1
            KEI2   = KEI1 + NMATAB(ISYMOP)
            KEND2  = KEI2 + NMATIJ(ISYMOP)
            LWRK2  = LWORK - KEND2
!
            IF ( RSPIM.AND.((.NOT.(ONLY21)))) THEN
!
               DTIME = SECOND() - DTIME
!
               CALL GPOPEN(LUE1,'CC_E1IM','UNKNOWN',' ','UNFORMATTED',
     *                     IDUMMY,.FALSE.)
               REWIND(LUE1)
               READ (LUE1)(WORK(KEI1 + J-1),J = 1,NMATAB(ISYMOP))
               CALL GPCLOSE(LUE1,'KEEP' )
!
               CALL GPOPEN(LUE2,'CC_E2IM','UNKNOWN',' ','UNFORMATTED',
     *                     IDUMMY,.FALSE.)
               REWIND(LUE2)
               READ (LUE2)(WORK(KEI2 + J-1),J = 1,NMATIJ(ISYMOP))
               CALL GPCLOSE(LUE2,'KEEP' )
!
               DTIME = SECOND() - DTIME
               TIMIO = TIMIO + DTIME
!
               IF (IPRINT.GT.150) THEN
                  CALL AROUND( 'E-intermediates read from disk ')
                  CALL CC_PREI(WORK(KEI1),WORK(KEI2),ISYMOP,1)
               ENDIF
!
            ENDIF
!
!---------------------------------------------------------------
!     Calculate E contributions to rho2+: (C2+)*EI(t1,t2).
!---------------------------------------------------------------
!
            KHLOG   = .TRUE.
            ANTISYM = .FALSE.
            ISYVEC = ISYMTR
            ISYMIM = ISYMOP
            DTIME = SECOND()
            CALL CCRHS_E3(WORK(KRHO2),KHLOG,WORK(KC2AM),WORK(KEI1),
     *                    WORK(KEI2),WORK(KEND2),LWRK2,ISYVEC,
     *                    ISYMIM,DUMMY,ANTISYM)
            DTIME = SECOND() - DTIME
            TIMEI  = TIMEI + DTIME
!
            IF (IPRINT .GT. 120) THEN
               CALL AROUND(' CC_RHTR:  AFTER E EI*(C2) RHO2 ')
               KOFF1 = KRHO1 + NT1AM(ISYMTR)*(IV - 1)
               CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,0,NC2)
            ENDIF
!
            IF ( LDEBUG ) THEN
               RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
               WRITE(LUPRI,*) 'Norm of Rho1 -after Ei*(c1,c2):    ',
     &                        RHO1N
               IF ( .NOT. CCS) THEN
                 RHO2N = DDOT(NT2AM(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
                 WRITE(LUPRI,*) 'Norm of Rho2+ -after Ei*(c1,c2):    ',
     &                          RHO2N
               ENDIF
            ENDIF
!
            ENDIF
!
!------------------------------------------
!    Calculate the C2+ dep. term for CC2.
!------------------------------------------
!
            IF (CC2) THEN
!
            DTIME = SECOND()
            ISIDE = 1
            CALL CC2_FCK(WORK(KRHO2),WORK(KC2AM),WORK(KEND2),LWRK2,
     *                   ISYMTR,WORK(KLAMDP),WORK(KLAMDH),ISIDE)
!
            ENDIF
!
!
!----------------------------------------------------
!     Readin  Gamma(t1,t2) vector: CCLR intermediate.
!----------------------------------------------------
!
            IF (RSPIM .AND. (.NOT. CC2)) THEN
!
               DTIME = SECOND()
!
               KGAMI  = KEND1
               KEND2  = KGAMI + NGAMMA(ISYMOP)
               LWRK2  = LWORK - KEND2
               CALL GPOPEN(LUGAM,'CC_GAMIM','UNKNOWN',' ','UNFORMATTED',
     *                     IDUMMY,.FALSE.)
               REWIND(LUGAM)
               READ (LUGAM)(WORK(KGAMI+J-1),J = 1,NGAMMA(ISYMOP))
               CALL GPCLOSE(LUGAM,'KEEP')
!
               DTIME = SECOND() - DTIME
               TIMIO = TIMIO + DTIME
!
!-------------------------------------------------
!        Calculate A-term: (C2+)*GAMMA(t1,t2).
!-------------------------------------------------
!
               ISYGAM = 1
               ISYVEC = ISYMTR
               IOPT    = 1
!
               DTIME  = SECOND()
               CALL CCRHS_A(WORK(KRHO2),WORK(KC2AM),WORK(KGAMI),
     *                      WORK(KEND2),LWRK2,ISYGAM,ISYVEC,IOPT)
               DTIME  = SECOND() - DTIME
               TIMA   = DTIME  + TIMA
!
               IF (IPRINT .GT. 120) THEN
                  CALL AROUND(' AFTER C2*GAMMA(t1,t2) A-term RHO is ')
                  KOFF1 = KRHO1 + NT1AM(ISYMTR)*(IV - 1)
                  CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,0,1)
               ENDIF
!
               IF ( LDEBUG ) THEN
                  RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,
     *                         WORK(KOFF1),1)
                  RHO2N = DDOT(NT2AM(ISYMTR),WORK(KRHO2),1,
     *                         WORK(KRHO2),1)
                  WRITE(LUPRI,*) 'Norm of Rho1 -after C2*GAM:        ',
     &                           RHO1N
                  WRITE(LUPRI,*) 'Norm of Rho2+ -after C2*GAM:       ',
     &                           RHO2N
               ENDIF
!
            ENDIF
!
!----------------------------------------------
!     Calculate the C-term: (C2+)*(-(1)C)
!----------------------------------------------
!
            IF (.NOT. CC2) THEN
!
            ISYDIM = 1
!
            FACTD = XMTWO
            FACTD2 = ZERO
            IOPT = 1
            IVECNR = 1
            ANTISYM = .FALSE.
            KHLOG = .TRUE.
            DTIME = SECOND()
!
            CALL CCRHS3_DIO(WORK(KRHO2),KHLOG,WORK(KC2AM),WORK(KLAMDH),
     *                      WORK(KEND1),LWRK1,ISYMTR,ISYDIM,LUCSIM,
     *                      CSIMFIL,IVECNR,IOPT,FACTD,ANTISYM,
     *                      DUMMY,FACTD2)
            DTIME = SECOND() - DTIME
            TIMDIO  = DTIME + TIMDIO
!
            IF (IPRINT .GT. 120) THEN
               CALL AROUND(' AFTER C2*DIM(t1,t2) D-term RHO is ')
               CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,0,1)
            ENDIF
!
            IF ( LDEBUG ) THEN
               RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
               RHO2N = DDOT(NT2AM(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
               WRITE(LUPRI,*) 'Norm of Rho1 -after C2*DIM:    ',RHO1N
               WRITE(LUPRI,*) 'Norm of Rho2+ -after C2*DIM:   ',RHO2N
            ENDIF
!
            ENDIF
!
!-------------------------------------
!     write out result vector.
!-------------------------------------
!
            DTIME   = SECOND()
            CALL CC_WVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                    NT2AM(ISYMTR),IV + ITR -1,0,WORK(KRHO2))
            DTIME   = SECOND() - DTIME
            TIMIO   = TIMIO    + DTIME
!
  170    CONTINUE
!
!==========================================================
!     MO basis T2SQ section.
!     Contract intermediates constructed in loop with T2SQ.
!==========================================================
!
!
!----------------------------------------------
!     Loop over number of trial vectors nsimtr.
!----------------------------------------------
!
         DO 180 IV = 1, NSIMTR
!
            KLPC1  = KLAMPC  + NGLMDT(ISYMTR)*(IV-1)
            KLHC1  = KLAMHC  + NGLMDT(ISYMTR)*(IV-1)
            KFCKC1 = KFOCKC  + N2BST(ISYMTR)*(IV-1)
            KEI1   = KEMAT1  + NEMAT1(ISYMTR)*(IV-1)
            KEI2   = KEMAT2  + NMATIJ(ISYMTR)*(IV-1)
            KGAMC1 = KGAMMC  + NGAMMA(ISYMTR)*(IV-1)
            KAMC1  = KC1AM   + NT1AM(ISYMTR)*(IV-1)
            KOFF1  = KRHO1   + NT1AM(ISYMTR)*(IV-1)
!
!------------------------------------
!        Readin T2 amplitude in RHO2.
!------------------------------------
!
            IF (.NOT. CC2) THEN
!
            DTIME = SECOND()
!
            IOPT = 2
            CALL CC_RDRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KRHO2))
!
            DTIME = SECOND() - DTIME
            TIMIO = TIMIO + DTIME
!
!--------------------------------
!        Square up T2 amplitudes.
!--------------------------------
!
            DTIME = SECOND()
            CALL CC_T2SQ(WORK(KRHO2),WORK(KC2AM),1)
            DTIME = SECOND() - DTIME
            TIMT2SQ = TIMT2SQ + DTIME
!
            IF (IPRINT.GT.50) THEN
               CALL AROUND('CC_RHTR: (T1,T2) vector readin')
               CALL CC_PRSQ(WORK(KT1AM),WORK(KC2AM),1,1,1)
            ENDIF
!
!--------------------------------------
!        Read in result vector.
!--------------------------------------
!
            DTIME = SECOND()
            CALL CC_RVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                    NT2AM(ISYMTR),IV + ITR  -1,0,WORK(KRHO2))
            DTIME = SECOND() - DTIME
            TIMIO = TIMIO + DTIME
!
!-----------------------------------
!     Calculate T2*GAMMA(C2+) A-term.
!-----------------------------------
!
            ISYGAM  = ISYMTR
            ISYVEC  = 1
            IOPT    = 1
!
            DTIME   = SECOND()
            CALL CCRHS_A(WORK(KRHO2),WORK(KC2AM),WORK(KGAMC1),
     *                   WORK(KEND1),LWRK1,ISYGAM,ISYVEC,IOPT)
            DTIME   = SECOND() - DTIME
            TIMA    = DTIME + TIMA
!
            IF (IPRINT .GT. 120) THEN
               CALL AROUND(' AFTER T2*GAMMA(C1,C2) A-term RHO ' )
               KOFF1 = KRHO1 + NT1AM(ISYMTR)*(IV - 1)
               CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,0,1)
            ENDIF
C
            IF ( LDEBUG ) THEN
               RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
               RHO2N = DDOT(NT2AM(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
               WRITE(LUPRI,*) 'Norm of Rho1 -after T2*GAM:        ',
     &                        RHO1N
               WRITE(LUPRI,*) 'Norm of Rho2+ -after T2*GAM:        ',
     &                        RHO2N
            ENDIF
!
!------------------------------------
!     Recover work from gamma matrix.
!------------------------------------
!
            IF (MINSCR) THEN
               KEND1 = KENDS2
               LWRK1 = LWRKS2
            ENDIF
!
!--------------------------------------------------
!     Calculate the C-term: T2*(1)C(tilde)(C1).
!--------------------------------------------------
!
            ISYCIM = ISYMTR
            ISYVEC = 1
!
            IOPT = 2
            NORMALCONT = .FALSE.
            NFACC = ZERO      ! Doesn't matter.
            PIJCONT = .TRUE.
            PFACC =  XMONE
            ANTISYM = .FALSE.
            DTIME  = SECOND()
            CALL CCSD_T2TP(WORK(KC2AM),WORK(KEND1),LWRK1,ISYVEC)
            CALL CCRHS3_CIO(WORK(KRHO2),WORK(KC2AM),WORK(KLAMDH),
     *                      WORK(KEND1),LWRK1,ISYVEC,ISYCIM,LUCS,CSFIL,
     *                      IV,IOPT,NFACC,NORMALCONT,PFACC,PIJCONT,
     *                      ANTISYM)
            CALL CCSD_T2TP(WORK(KC2AM),WORK(KEND1),LWRK1,ISYVEC)
            DTIME   = SECOND() - DTIME
            TIMCIO  = DTIME + TIMCIO
!
            IF (IPRINT .GT. 120) THEN
               CALL AROUND(' AFTER T2*CIM(C1) C-term RHO is ')
               CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,0,1)
            ENDIF
!
            IF ( LDEBUG ) THEN
               RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
               RHO2N = DDOT(NT2AM(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
               WRITE(LUPRI,*) 'Norm of Rho1 -after T2*CI:     ',RHO1N
               WRITE(LUPRI,*) 'Norm of Rho2+ -after T2*CI:    ',RHO2N
            ENDIF
!
!-----------------------------
!   (.NOT. CC2) ends here.
!-----------------------------
!
            ENDIF
!
!---------------------------------
!           For CC2 read in the
!           result vector.
!---------------------------------
!
            IF (CC2) THEN
!
              DTIME = SECOND()
!
              CALL CC_RVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                      NT2AM(ISYMTR),IV + ITR -1,0,WORK(KRHO2))
              DTIME = SECOND() - DTIME
              TIMIO = TIMIO + DTIME
!
              IF ( LDEBUG ) THEN
                 RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,
     *                                      WORK(KOFF1),1)
                 WRITE(LUPRI,*) 
     *              'Norm of Rho1 after loop in mo:     ',RHO1N
                 IF ( .NOT. CCS) THEN
                 RHO2N = DDOT(NT2AM(ISYMTR),WORK(KRHO2),1,
     *                                      WORK(KRHO2),1)
                 WRITE(LUPRI,*) 
     *              'Norm of Rho2+ after loop in mo:    ',RHO2N
                 ENDIF
              ENDIF
!
            ENDIF
!
!--------------------------------------------
!        For CC3 calc. the correction.
!        The Omega2- cont. is saved to file
!        inside the subroutine.
!--------------------------------------------
!
            IF (CC3) THEN
!
               IF (MINSCR) THEN
                  KEND1 = KENDS2
                  LWRK1 = LWRKS2
               ENDIF
!
               KENDCC3  = KEND1
               LWRKCC3  = LWRK1
!
!---------------------------------------------
!     Allocate two R2
!---------------------------------------------
!
               KC2AMP = KEND1
               KC2AMM = KC2AMP + NC2AM
               KRHO22 = KC2AMM + NC2AM
               KRHO3  = KRHO22 + NT2AMA(ISYMTR)
               KEND1  = KRHO3  + NRHO2
               LWRK1  = LWORK  - KEND1
!
               IF (LWRK1 .LE. 0) THEN
                 WRITE(LUPRI,*)'LWORK :',LWORK
                 WRITE(LUPRI,*)'KEND1 :',KEND1
                 CALL QUIT('Insufficient space in CC3 (cc_rhtr3.F)')
               ENDIF
!
!------------------------------------
!        Readin T2 amplitude in KRHO3.
!------------------------------------
!
               DTIME = SECOND()
!
               IOPT = 2
               CALL CC_RDRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KRHO3))
!
               DTIME = SECOND() - DTIME
               TIMIO = TIMIO + DTIME
!
!--------------------------------
!        Square up T2 amplitudes.
!--------------------------------
!
               DTIME = SECOND()
               CALL CC_T2SQ(WORK(KRHO3),WORK(KC2AM),1)
               DTIME = SECOND() - DTIME
               TIMT2SQ = TIMT2SQ + DTIME
!
               IF (IPRINT.GT.50) THEN
                  CALL AROUND('CC_RHTR: (T1,T2) vector readin')
                  CALL CC_PRSQ(WORK(KT1AM),WORK(KC2AM),1,1,1)
               ENDIF
!
!---------------------------------------------
!                 Read in (-) result vector.
!---------------------------------------------
!
               DTIME = SECOND()
               CALL CC_RVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                       NT2AMA(ISYMTR),IV + ITR -1,NT2AM(ISYMTR),
     *                       WORK(KRHO22))
               DTIME = SECOND() - DTIME
               TIMIO = TIMIO + DTIME
!
!------------------------------------------
!     Approximative triples correction:
!     <mu2|[[H,R1],T3(T2)]|HF>
!     T2AM is in C2AM.
!------------------------------------------
!
               CC3LR   = .TRUE.
               ECURR2  = ECURR
               ECURR   = 0.0D0
               ISYMT1  = ISYMTR
               ISYMT2  = 1
               DTIME   = SECOND()
!
               CALL CC3_OMEG3(0.0D0,WORK(KOFF1),WORK(KRHO2),
     *                        WORK(KRHO22),
     *                        WORK(KAMC1),ISYMT1,WORK(KC2AM),ISYMT2,
     *                        WORK(KFOCK0),WORK(KLAMDP),WORK(KLAMDH),
     *                        WORK(KEND1),LWRK1,LUSRTT1,FNSRTT1,
     *                        LUTOC,FNTOC,LUCKJDT1,FNCKJDT1)
!
               DTIME   = SECOND() - DTIME
               TIMOM32 = TIMOM32 + DTIME
               ECURR   = ECURR2
!
               IF (IPRINT .GT. 55) THEN
                 RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
                 WRITE(LUPRI,1) 'Norm of Rho1  after cc3_omeg3: ',RHO1N
                 RHO2N = DDOT(NT2AM(ISYMTR),WORK(KRHO2),1,
     *                                      WORK(KRHO2),1)
                 WRITE(LUPRI,1) 'Norm of Rho2+ after cc3_omeg3: ',RHO2N
                 RHO2N = DDOT(NT2AM(ISYMTR),WORK(KRHO22),1,
     *                                      WORK(KRHO22),1)
                 WRITE(LUPRI,1) 'Norm of Rho2- after cc3_omeg3: ',RHO2N
               ENDIF
!
!------------------------------------
!        Readin T2 amplitude in KRHO3.
!------------------------------------
!
               DTIME = SECOND()
!
               IOPT = 2
               CALL CC_RDRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KRHO3))
!
               DTIME = SECOND() - DTIME
               TIMIO = TIMIO + DTIME
!
!--------------------------------
!        Square up T2 amplitudes.
!--------------------------------
!
               DTIME = SECOND()
               CALL CC_T2SQ(WORK(KRHO3),WORK(KC2AM),1)
               DTIME = SECOND() - DTIME
               TIMT2SQ = TIMT2SQ + DTIME
!
               IF (IPRINT.GT.50) THEN
                  CALL AROUND('CC_RHTR: (T1,T2) vector readin')
                  CALL CC_PRSQ(WORK(KT1AM),WORK(KC2AM),1,1,1)
               ENDIF
!
!-----------------------------------------------
!     Readin and square up the C2- amplitudes.
!-----------------------------------------------
!
               DTIME = SECOND()
               CALL CC_RVEC3(LUFC2,FC2AM,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                       NT2AMA(ISYMTR),IVEC+IV-1,NT2AM(ISYMTR),
     *                       WORK(KRHO3))
               DTIME = SECOND() - DTIME
               TIMIO = TIMIO + DTIME
!
!
               CALL CCLR_DIASCL(WORK(KRHO3),ZERO,ISYMTR)
               CALL CC_T2SQ3(WORK(KRHO3),WORK(KC2AMM),ISYMTR) 
!
!----------------------------------------------
!     Readin and square up the C2+ amplitudes.
!----------------------------------------------
!
               DTIME = SECOND()
               CALL CC_RVEC3(LUFC2,FC2AM,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                       NT2AM(ISYMTR),IVEC+IV-1,0,
     *                       WORK(KRHO3))
               DTIME = SECOND() - DTIME
               TIMIO = TIMIO + DTIME
!
!
               CALL CCRHS3_R2IJ(WORK(KRHO3),WORK(KEND1),LWRK1,ISYMTR)
               CALL CC_T2SQ(WORK(KRHO3),WORK(KC2AMP),ISYMTR)
!
!-----------------------------------------------
!     Calc. the cont.
!-----------------------------------------------
!
               CC3LR   = .FALSE.
               ISYMT1  = ISYMTR
               ISYMT12 = 1
               ISYMT2  = 1
               ISYMC2  = ISYMTR
!
               CALL CC3_OMEG33(ECURR,WORK(KOFF1),
     *                        WORK(KRHO2),WORK(KRHO22),
     *                        WORK(KAMC1),ISYMT1,WORK(KC2AM),ISYMT2,
     *                        WORK(KC2AMP),WORK(KC2AMM),ISYMC2,ISYMT12,
     *                        WORK(KFOCK0),WORK(KLAMDP),WORK(KLAMDH),
     *                        WORK(KEND1),LWRK1,LUSRTT1,FNSRTT1,
     *                        LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                        LUCKJDT1,FNCKJDT1,LUCKJDT2,FNCKJDT2,
     *                        LUCKJDT3,FNCKJDT3,LUSRTT2,FNSRTT2,
     *                        LUSRTT3,FNSRTT3,LUTOC,FNTOC)
!
               IF (IPRINT .GT. 55) THEN
                 RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,WORK(KOFF1),1)
                 WRITE(LUPRI,1) 'Norm of Rho1  after cc3_omeg33: ',RHO1N
                 RHO2N = DDOT(NT2AM(ISYMTR),WORK(KRHO2),1,
     *                                      WORK(KRHO2),1)
                 WRITE(LUPRI,1) 'Norm of Rho2+ after cc3_omeg33: ',RHO2N
                 RHO2N = DDOT(NT2AM(ISYMTR),WORK(KRHO22),1,
     *                                      WORK(KRHO22),1)
                 WRITE(LUPRI,1) 'Norm of Rho2- after cc3_omeg33: ',RHO2N
               ENDIF
!
               IF (IPRINT .GT. 110) THEN
                  CALL AROUND('Omega1 & Omega2+ after CC3')
                  KOFF1 = KRHO1 + NT1AM(ISYMTR)*(IV - 1)
                  CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,1,1)
                  CALL AROUND('Omega2- after CC3')
                  KOFF1 = KRHO1 + NT1AM(ISYMTR)*(IV - 1)
                  CALL CC_PRP(WORK(KOFF1),WORK(KRHO22),ISYMTR,0,1)
               ENDIF
!
!----------------------------------------
!     Write omega1 and omega2- to disk.
!----------------------------------------
!
               CALL CC_WVEC(LUFR1,FRHO1,NT1AM(ISYMTR),NT1AM(ISYMTR),
     *                      IV + IVEC -1,WORK(KOFF1)) 
               CALL CC_WVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                       NT2AMA(ISYMTR),IV + ITR -1,NT2AM(ISYMTR),
     *                       WORK(KRHO22))
!
!----------------------------------
!     CC3 and friends ends here
!----------------------------------
!
               KEND1 = KENDCC3
               LWRK1 = LWRKCC3
!
            ENDIF
!
!----------------------------------------
!        Finally scale RHO2+ by 1/2 
!        and use P(tilde)(ij)
!----------------------------------------
!
            CALL DSCAL(NT2AM(ISYMTR),HALF,WORK(KRHO2),1)
!
            IF ( LDEBUG ) THEN
               RHO2N = DDOT(NT2AM(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
               IF (.NOT. CC3) THEN
                  WRITE(LUPRI,*) 
     *                 'Norm of Rho2+ -after scaling:   ',RHO2N
               ELSE
                  WRITE(LUPRI,*) 
     *                 'Norm of Rho2+ -after CC3:       ',RHO2N
               ENDIF
            ENDIF
!
            IF (.NOT. CCP2) THEN
               CALL CCRHS3_IJ(WORK(KRHO2),WORK(KEND1),LWRK1,ISYMTR)
            ENDIF
!
            IF (IPRINT .GT. 150) THEN
               CALL AROUND(' After CCRHS3_IJ RHO2+ is :')
               KOFF1 = KRHO1 + NT1AM(ISYMTR)*(IV - 1)
               CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,0,1)
            ENDIF
!
!
            IF ( LDEBUG ) THEN
               RHO2N = DDOT(NT2AM(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
               WRITE(LUPRI,*) 'Norm of Rho2+ -after CCRHS3_IJ:       ',
     &                        RHO2N
            ENDIF
!
!
!--------------------------------
!        write out result vector.
!--------------------------------
!
            DTIME   = SECOND()
            CALL CC_WVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                    NT2AM(ISYMTR),IV + ITR -1,0,WORK(KRHO2))
            DTIME   = SECOND() - DTIME
            TIMIO   = TIMIO    + DTIME
!
!=============================================================
!     End section for RHO2+.
!=============================================================
!
!
            IF (IPRINT .GT. 50 ) THEN
               CALL AROUND('END OF CC_RHTR3 :RHO ')
               NC2 = 1
               IF ( CCS ) NC2 = 0
               CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYMTR,NSIMTR,NC2)
            ENDIF
!
            IF ( LDEBUG ) THEN
                RHO1N = DDOT(NT1AM(ISYMTR),WORK(KOFF1),1,
     *                       WORK(KOFF1),1)
                WRITE(LUPRI,*) 'Norm of Rho1 for IV = ',IV,
     *                  ' -end of cclr_rhtr3 (+part2):    ',RHO1N
                IF (.NOT. CCS) THEN
                    RHO2N = DDOT(NT2AM(ISYMTR),WORK(KRHO2),1,
     *                           WORK(KRHO2),1)
                    WRITE(LUPRI,*) 'Norm of Rho2+ for IV = ',IV,
     *                         'end of cclr_rhtr3 (+part2):    ',RHO2N
                ENDIF
            ENDIF
!
  180    CONTINUE
!
      ENDIF
!
!-----------------------
!     Timings.
!-----------------------
!
      TIMALL  = SECOND() - TIMALL
!
      IF (IPRINT .GT. 4) THEN
         TIMCCSD = TIMA + TIMAOD + TIMBF + TIMC + TIMCIO + TIMD
     *           + TIMDIO + TIMEI + TIMER1 + TIMER2 + TIMLAM
     *           + TIMFCK + TIMFCKMO + TIMG + TIMH + TIMI
     *           + TIMIO + TIMRDAO + TIMT2AO + TIMT2MO + TIMT2SQ
     *           + TIMTRBT
       ENDIF
!
      IF (IPRINT .GT. 4) THEN
         WRITE(LUPRI,'(/)')
         WRITE(LUPRI,9998) 'CC_RHTR3 in total ', TIMALL
         WRITE(LUPRI,9998) 'CC part           ', TIMCCSD
         WRITE(LUPRI,9998) 'Int. calc. & read ', TIMER1 + TIMER2
     &                  + TIMRDAO
         WRITE(LUPRI,9998) 'IO of vectors/I.M.', TIMIO
         WRITE(LUPRI,*) ' '
      ENDIF
!
      IF (IPRINT .GT. 9) THEN
         WRITE(LUPRI,'(A)')
         WRITE(LUPRI,9999) 'A               ', TIMA
         WRITE(LUPRI,9999) 'AOD             ', TIMAOD
         WRITE(LUPRI,9999) 'BF              ', TIMBF
         WRITE(LUPRI,9999) 'C               ', TIMC
         WRITE(LUPRI,9999) 'CIO             ', TIMCIO
         WRITE(LUPRI,9999) 'D               ', TIMD
         WRITE(LUPRI,9999) 'DIO             ', TIMDIO
         WRITE(LUPRI,9999) 'E               ', TIMEI
         WRITE(LUPRI,9999) 'FCK             ', TIMFCK
         WRITE(LUPRI,9999) 'FCKMO           ', TIMFCKMO
         WRITE(LUPRI,9999) 'G               ', TIMG
         WRITE(LUPRI,9999) 'H               ', TIMH
         WRITE(LUPRI,9999) 'I               ', TIMI
         WRITE(LUPRI,9999) 'LAMBDA          ', TIMLAM
         WRITE(LUPRI,9999) 'T2AO            ', TIMT2AO
         WRITE(LUPRI,9999) 'T2MO            ', TIMT2MO
         WRITE(LUPRI,9999) 'T2SQ            ', TIMT2SQ
         IF (CCSDT) THEN
            WRITE(LUPRI,9999) 'Triples corr.   ', TIMOM32
            WRITE(LUPRI,9999) 'Triples integr. ', TIMOM32
         ENDIF
!
         WRITE(LUPRI,'(A)')
         WRITE(LUPRI,9999) 'ERIDI1          ', TIMER1
         WRITE(LUPRI,9999) 'ERIDI2          ', TIMER2
         WRITE(LUPRI,9999) 'RDAO            ', TIMRDAO
         WRITE(LUPRI,9999) 'TR. BT.         ', TIMTRBT
!
         WRITE(LUPRI,*) ' '
      ENDIF
      IF (IPRINT .GT. 15) THEN
         IF (DIRECT) WRITE(LUPRI,*) ' Atomic direct calculation'
         CALL AROUND(' END OF CC_RHTR3 ')
      ENDIF
!
!--------------------------------------------
!     Close files 
!--------------------------------------------
!
      IF (DUMPCD .AND. (.NOT. CCS)) THEN
         CALL WCLOSE2(LUB3,BFFIL,'DELETE')
         CALL WCLOSE2(LUC,CTFIL,'DELETE')
         CALL WCLOSE2(LUD,DTFIL,'DELETE')
         CALL WCLOSE2(LUCS,CSFIL,'DELETE')
         CALL WCLOSE2(LUP,PFIL,'DELETE')
!
         IF (.NOT. CC2) THEN
           CALL WCLOSE2(LUDTIM,DTIMFIL,'KEEP')
           CALL WCLOSE2(LUCSIM,CSIMFIL,'KEEP')
         ENDIF
!
      ENDIF
!
      IF (CCSDT) THEN
         CALL WCLOSE2(LUSRTT1,FNSRTT1,'DELETE')
         CALL WCLOSE2(LU3VI,FN3VI,'KEEP')
         CALL WCLOSE2(LU3VI2,FN3VI2,'KEEP')
         CALL WCLOSE2(LUCKJDT1,FNCKJDT1,'DELETE')
         CALL WCLOSE2(LUCKJDT2,FNCKJDT2,'DELETE')
         CALL WCLOSE2(LUCKJDT3,FNCKJDT3,'DELETE')
         CALL WCLOSE2(LUSRTT2,FNSRTT2,'DELETE')
         CALL WCLOSE2(LUSRTT3,FNSRTT3,'DELETE')
         CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
      ENDIF
!
!-------------------------------------
!     Formats.
!-------------------------------------
!
      CALL QEXIT('CC_RHTR3')
!
   1  FORMAT(1x,A35,1X,E20.10)
9998  FORMAT(1x,'Time used in',2x,A20,2x,': ',f10.2,' seconds')
9999  FORMAT(1x,'Time used in',2x,A18,2x,': ',f10.2,' seconds')
!
	RETURN
	END
