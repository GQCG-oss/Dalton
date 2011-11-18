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
*=====================================================================*
      SUBROUTINE CC_FMATTST(LLIST,ILLNR,RLIST,IRLNR,WORK,LWORK)
*--------------------------------------------------------------------*
*
*     Purpose: transformation of trial vectors with the CC F-matrix
*
*                 Gamma(mu) = <L|[[H,R],Tau,mu]|CC>
*
*     Left hand vector is read in according to LLIST,ILLNR.
*     Right vector according to RLIST,IRLNR.
*    
*     For LLIST .NE. 'L0' skip the HF term.
*
*     Result vector is returned in the beginning of the work space
* 
*     Written by Ove Christiansen April-october 1996
*     Heavily revised to reduce operation count by C. Haettig July 1998
*     `Left' B intermediate revised in November 1998, C. Haettig
*---------------------------------------------------------------------*
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxash.h>
#include <maxorb.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
#include <ccorb.h>
#include <ccisao.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccfield.h>
#include <cclr.h>
#include <blocks.h>
#include <ccsdio.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
#include <ccroper.h>
#include <ccinftap.h>
#include <r12int.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, XMTWO= -2.0D0)
      PARAMETER (XMHALF=-0.5D0,XHALF=0.5D0,THREE=3.0D0,XMONE=-1.0D0)
C
      LOGICAL NEWZWV, ITST, ATST, E2TST, HTST, E22TST, C12TST, A11TST
      LOGICAL D2TST, C2TST, DOUBLETST, SINGLETST, G21TST, B21TST
      PARAMETER (NEWZWV = .TRUE.)
C
      INTEGER LUTR, LUPQR0, LUPQRR,  LUBFI, LUBFD
C
      CHARACTER*(6) TRFIL, FILPQIM, FILXYM, FNBFI
      CHARACTER*(8) FNBFD
      PARAMETER (FILPQIM = 'CCPQ00', TRFIL = 'CC_TRA', FILXYM='CC_XYM')
      PARAMETER (FNBFI   = 'CCFBFI', FNBFD = 'CCBFDENS')
C
      INTEGER ISYM0
      PARAMETER ( ISYM0 = 1 )
C
      INTEGER LWORK
      DIMENSION INDEXA(MXCORB)
      DIMENSION IADRPQ0(MXCORB+IRAT), IADRPQR(MXCORB+IRAT)
      DIMENSION IADRBFI(MXCORB), IADRBFD(MXCORB)
      DIMENSION WORK(*)
      CHARACTER*6 CFIL,DFIL,CTFIL,DTFIL
      CHARACTER*7 O3FIL,O3TFIL
      CHARACTER LLIST*(*),RLIST*(*),MODEL*10,MODELX*10
      LOGICAL FCKCON,ETRAN,L1TST,L2TST
      LOGICAL T2TSAV,ORSAVE
      LOGICAL LGAMMA, LGIM, LO3BF, LBFZETA
      INTEGER ISYMR, ISYML, ISYRES, ISYMH, ISTARTBFI, ISTARTBFD
      CHARACTER*6 FILPQR0, FILPQRR
      INTEGER NBFRHF(8), IBFRHF(8,8), NULL
C
      CALL QENTER('CC_FMATTST')
C
C-----------------------------------
C     initialize overall timing
C-----------------------------------
C
      TIMALL  = SECOND()
      TIMIO   = ZERO
      TIMT2SQ = ZERO
      NULL    = 0
C
C-----------------------------------
C     echo input
C-----------------------------------
C
      IF (DEBUG) THEN
        WRITE (LUPRI,*) 'entered CC_FMAT:'
        WRITE (LUPRI,*) 'RLIST, IRLNR :',RLIST,IRLNR
        WRITE (LUPRI,*) 'LLIST, ILLNR :',LLIST,ILLNR
        WRITE (LUPRI,*) 'WORK(1):',WORK(1)
        WRITE (LUPRI,*) 'LWORK:',LWORK
      END IF
C
      DEBUG = .TRUE.
C
C---------------------------------------------
C     Set symmetries, test flags and CC model:
C---------------------------------------------
C
      ISYML = ILSTSYM(LLIST,ILLNR) 
      ISYMR = ILSTSYM(RLIST,IRLNR) 
      ISYRES = MULD2H(ISYML,ISYMR) 
C
C     ----------------------------------------------------------
C     initialize file addresses for BF density and intermediate:
C     ----------------------------------------------------------
C
      ISTARTBFI = 1
      ISTARTBFD = 1
C
C     --------------------------------------
C     precalculate symmetry array for BFRHF:
C     --------------------------------------
C
      DO ISYM = 1, NSYM
        ICOUNT = 0
        DO ISYMAK = 1, NSYM
           ISYBET = MULD2H(ISYMAK,ISYM)
           IBFRHF(ISYMAK,ISYBET) = ICOUNT
           ICOUNT = ICOUNT + NT1AO(ISYMAK)*NBAS(ISYBET)
        END DO
        NBFRHF(ISYM) = ICOUNT
      END DO
C
      L1TST = .FALSE.
      L2TST = .FALSE.
C
      MODEL = 'CCSD      '
      IF (CCS) MODEL = 'CCS       '
      IF (CC2) MODEL = 'CC2       '
C
      IF (ISYRES .NE. MULD2H(ISYML,ISYMR)) THEN
         WRITE(LUPRI,*) 'Symmetry mismatch in CC_FMAT'
         CALL QUIT('Symmetry mismatch in CC_FMAT')
      ENDIF
C
      IF (IPRINT .GT. 15) THEN
         CALL AROUND(' START OF CC_FMATTST ')
         IF (DIRECT) WRITE(LUPRI,*) ' Atomic direct transformation'
      ENDIF
C
      IF ( CCSDT ) THEN
         WRITE(LUPRI,*) ' F-Matrix transformations not implemented '
     *              //'for triples yet: Go and do it.'
         CALL QUIT(' no triples F-matrix')
      ENDIF 
C
C---------------------------------------
C     Set common block control variables
C---------------------------------------
C
      T2TSAV = T2TCOR
      IF (CCS .OR. CC2) T2TCOR = .FALSE.
C     COMMENTED OUT! OBSOLETE SUBROUTINE!
C     IF (CC2 .AND.(NSIDE .EQ. 2)) T2TCOR = T2TSAV
C
      ORSAVE = OMEGOR
      IF (CCS .OR. CC2) THEN
         OMEGOR = .FALSE.
      ELSE
         OMEGOR = .TRUE.
      ENDIF
C
      IF  (.NOT.DUMPCD) THEN
         CALL QUIT('CC_FMAT requires DUMPCD=.TRUE.')
      END IF
C
C--------------------------------------------------------------
C     open files for result vector and for BFZeta intermediate:
C--------------------------------------------------------------
C
      IERRTR = 0
      IERRBF = 0
C
      LUTR  = -1
      IF ( .NOT. CCS ) THEN
         CALL WOPEN2(LUTR, TRFIL,64,0)
      ENDIF
C
      LUBFI = -1
      LUBFD = -1
      IF (CCSD) THEN
         CALL WOPEN2(LUBFI,FNBFI,64,0)
         CALL WOPEN2(LUBFD,FNBFD,64,0)
      END IF
C
*--------------------------------------------------------------------*
*     calculate P and Q intermediates, write them to file and
*     open the file for read within the loop over AO integrals
*--------------------------------------------------------------------*
      TIMEZ = SECOND()

      IF (.NOT. (CCS.OR.CC2)) THEN

        IF (LLIST(1:2).EQ.'L0') THEN
          FILPQR0 = FILPQIM
        ELSE
          FILPQR0 = 'CCPQR0'
          CALL CC_PQIM(LLIST,ILLNR,'R0',0,FILPQR0,WORK,LWORK)
        END IF

        FILPQRR = 'CCPQRR'
        CALL CC_PQIM(LLIST,ILLNR,RLIST,IRLNR,FILPQRR,WORK,LWORK)

C       -----------------------------------------
C       open file for P and Q intermediates:
C       -----------------------------------------
        LUPQR0 = -1
        LUPQRR = -1
        CALL WOPEN2(LUPQR0,FILPQR0,64,0)
        CALL WOPEN2(LUPQRR,FILPQRR,64,0)

        CALL GETWA2(LUPQR0,FILPQR0,IADRPQ0,1,(MXCORB+IRAT)/IRAT)
        CALL GETWA2(LUPQRR,FILPQRR,IADRPQR,1,(MXCORB+IRAT)/IRAT)

      END IF

      TIMEZ = SECOND() - TIMEZ

C
C--------------------------------------------------------------------
C     Allocate space for solution vector and right hand trial vector.
C--------------------------------------------------------------------
C
      IF (CCS) THEN
         NRHO2 = 2
      ELSE IF (CC2) THEN
         NRHO2 = MAX(NT2AM(ISYRES),NT2AM(ISYMR),NT2AM(ISYML),NT2AM(1))
      ELSE
         NRHO2 = MAX(NT2AM(ISYRES),NT2AM(ISYMR),NT2AM(ISYML),NT2AM(1),
     *               2*NT2ORT(ISYRES),2*NT2ORT(ISYMR),NT2AO(ISYRES)) 
      ENDIF 
C
      NC2AM = MAX(NT2SQ(ISYMR),NT2SQ(1),NT2AM(ISYMR)+2*NT2ORT(ISYMR),
     *            NT2AM(1)+2*NT2ORT(1),NT2AM(ISYRES))
      IF ( CC2 ) NC2AM = MAX(NT2SQ(ISYMR),NT2SQ(1))
      IF ( CCS ) NC2AM = 2
C
      KRHO1   = 1
      KRHO2   = KRHO1   + NT1AM(ISYRES)
      KC1AM   = KRHO2   + NRHO2
      KC2AM   = KC1AM   + NT1AM(ISYMR)
      KL1AM   = KC2AM   + NC2AM
      KL2AM   = KL1AM   + NT1AM(ISYML)
      KL1R2   = KL2AM   + NT2SQ(ISYML)
      KL1R2PK = KL1R2   + N2BST(ISYRES)
      KT1AM   = KL1R2PK + NNBST(ISYRES)
      KEND0   = KT1AM   + NT1AM(ISYM0)
      LWRK0   = LWORK   - KEND0
      IF (LWRK0 .LE. 0 ) THEN
         CALL QUIT('Too little workspace in CC_FMAT.')
      END IF
C
C-------------------------------------------------------------------
C     Read the lagrangian multiplier vector from disk and square up:
C-------------------------------------------------------------------
C
      IF (.NOT. ( CCS .AND. LLIST(1:2).EQ.'L0') ) THEN
        DTIME = SECOND()
        IOPT = 3
        CALL CC_RDRSP(LLIST,ILLNR,ISYML,IOPT,MODEL,
     *                WORK(KL1AM),WORK(KRHO2))
        DTIME = SECOND() - DTIME
        TIMIO = TIMIO + DTIME
      END IF

      IF (.NOT. CCS)  THEN
         DTIME = SECOND()
         CALL CC_T2SQ(WORK(KRHO2),WORK(KL2AM),ISYML)
         DTIME = SECOND() - DTIME
         TIMT2SQ = TIMT2SQ + DTIME
      END IF
C
C     ------------------
C     Test options.
C     ------------------
C
      IF (L1TST .AND. (.NOT. CCS)) CALL DZERO(WORK(KL2AM),NT2SQ(ISYML))
      IF (L2TST)                   CALL DZERO(WORK(KL1AM),NT1AM(ISYML)) 
 
      IF ( DEBUG ) THEN
         RHO1N = DDOT(NT1AM(ISYML),WORK(KL1AM),1,WORK(KL1AM ),1)
         WRITE(LUPRI,1) 'Norm of L1AM as read from file:     ',RHO1N
      ENDIF
      IF ( DEBUG .AND. (.NOT. CCS)) THEN
         RHO2N = DDOT(NT2SQ(ISYML),WORK(KL2AM),1,WORK(KL2AM),1)
         WRITE(LUPRI,1) 'Norm of L2AM -after  square:       ',RHO2N
      ENDIF
C
C-----------------------------------------------------------------
C     Read the CC reference amplitudes from disk, T2AMP0 in KRHO2
C     (double excitation part only needed if LLIST .NE. 'L0',
C      i.e. if X, Y and M intermediates have to be calculated)
C-----------------------------------------------------------------
C
      IOPT = 1
      IF (LLIST(1:2).NE.'L0') IOPT = 3

      DTIME = SECOND()
      CALL CC_RDRSP('R0',0,1,IOPT,MODELX,WORK(KT1AM),WORK(KRHO2))
      DTIME = SECOND() - DTIME
      TIMIO = TIMIO + DTIME
C
C-------------------------------------------------------------------
C     initialize C1 and C2 vectors:
C     we start the calculations with C2 vector packed in WORK(KC2AM)
C-------------------------------------------------------------------
C
      DTIME = SECOND()
      IOPT = 3
      CALL CC_RDRSP(RLIST,IRLNR,ISYMR,IOPT,MODELX,
     *              WORK(KC1AM),WORK(KC2AM))
      DTIME = SECOND() - DTIME
      TIMIO = TIMIO + DTIME

      IF (.NOT. CCS) CALL CCLR_DIASCL(WORK(KC2AM),TWO,ISYMR)

      IF ( DEBUG ) THEN
         RHO1N = DDOT(NT1AM(ISYMR),WORK(KC1AM),1,WORK(KC1AM),1)
         WRITE(LUPRI,1) 'CC_FMAT: Norm of T1-response vector:',RHO1N
         IF (.NOT. CCS) THEN
            RHO2N = DDOT(NT2AM(ISYMR),WORK(KC2AM),1,WORK(KC2AM),1)
            WRITE(LUPRI,1) 'CC_FMAT: Norm of T2-response vector:',RHO2N
         ENDIF
      ENDIF
C
C-------------------------------------------------------------------
C     initialize result vector:
C     (note that double-excitation result vector is kept on file...)
C-------------------------------------------------------------------
C
      CALL DZERO(WORK(KRHO1),NT1AM(ISYRES))
C
C------------------------------------------------------------
C     Calculate <HF| [[H,tau,mu],R1]|HF> contribution.
C     Note inside FHF a Integral(ai,bj) is allocated. 
C     Of course released afterwards.
C------------------------------------------------------------
C
      IF ((LLIST(1:2).EQ.'L0')) THEN
 
         CALL CC_FHF(WORK(KC1AM),WORK(KRHO1),ISYMR,WORK(KEND0),LWRK0)
 
         IF (DEBUG) THEN
            WRITE(LUPRI,*) 
            RHO1N = DDOT(NT1AM(ISYRES),WORK(KRHO1),1,WORK(KRHO1),1)
            WRITE(LUPRI,1) 'Norm of Rho1 -after FHF:           ',RHO1N
         ENDIF
 
      ENDIF
C
C------------------------------------------
C     If CCS all is done and we can return.
C------------------------------------------
C
      IF ( CCS .AND. (LLIST(1:2).EQ.'L0')) THEN
         RETURN
      ENDIF
C
C----------------
C     Open files.
C----------------
C
      IF (.NOT.(CCS .OR. CC2)) THEN
         LUC = -1
         LUD = -1
C
         CTFIL = 'CCLR_C'
         DTFIL = 'CCLR_D'
C
         CALL WOPEN2(LUC,CTFIL,64,0)
         CALL WOPEN2(LUD,DTFIL,64,0)
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
      END IF
C
C-----------------------------------------
C     Open scratch file CCINT3O and O3TFIL
C-----------------------------------------
C
      LUO3 = -1
      O3FIL = 'CCINT3O'
C
      CALL WOPEN2(LUO3,O3FIL,64,0)
C
      LUO3T =  -1
      O3TFIL = 'CCO3INT'
C
      CALL WOPEN2(LUO3T,O3TFIL,64,0)
C
C------------------------
C     Time Initialiation.
C------------------------
C
      TIMER1    = 0.0D00
      TIMER2    = 0.0D00
      TIMRDAO   = 0.0D00
C
      TIMA      = 0.0D00
      TIMBF     = 0.0D00
      TIMC      = 0.0D00
      TIMD      = 0.0D00
      TIME      = 0.0D00
      TIMF      = 0.0D00
      TIMG      = 0.0D00
      TIMH      = 0.0D00
      TIMI      = 0.0D00
      TIMEI     = 0.0D00
      TIMEXI    = 0.0D00
      TIMEYI    = 0.0D00
      TIMEMI    = 0.0D00
      TIMEH     = 0.0D00
      TIMEFCK   = 0.0D00
      TIME1C1   = 0.0D00
      TIM1C1F   = 0.0D00
      TIMCIO    = 0.0D00
      TIMDIO    = 0.0D00
      TIM22D    = 0.0D00
      TIMEC     = 0.0D00
      TIMEG     = 0.0D00
      TIM2EM    = 0.0D00
      TIME3O    = 0.0D00
C
      TIMLAM    = 0.0D00
      TIMAOD    = 0.0D00
      TIMFCK    = 0.0D00
      TIMTRBT   = 0.0D00
      TIMT2AO   = 0.0D00
      TIMTCME   = 0.0D00
      TIMT2TR   = 0.0D00
      TIMT2SQ   = 0.0D00
      TIMT2BT   = 0.0D00
      TIMT2MO   = 0.0D00
      TIMFCKMO  = 0.0D00
      TIMC1T2   = 0.0D00
      TIM12B    = 0.0D00
C
      TIMT3I    = 0.0D00
      TIMO31    = 0.0D00
      TIMO32    = 0.0D00
      TIMO33    = 0.0D00
C
      TIMIO     = 0.0D00
C
C-----------------------------------------------------
C     Work space allocation for general intermediates.
C-----------------------------------------------------
C
      NE1TOT = MAX(NEMAT1(ISYMR),NEMAT1(1))
      NE2TOT = MAX(NMATIJ(ISYMR),NMATIJ(1))
      N2C2C2 = 0
      IF ( T2TCOR ) N2C2C2 = NT2SQ(ISYMR)
      IF ( CCS ) THEN
         N2C2C2 = 2
      ENDIF
C
      IF ( IPRINT .GT. 25)
     *   WRITE(LUPRI,*) ' In CC_FMAT work in start is:',LWORK
C
      K2C2C2  = KEND0
      KLAMDP  = K2C2C2  + N2C2C2
      KLAMDH  = KLAMDP  + NLAMDT
      KDENSI  = KLAMDH  + NLAMDT
      KDNSPKI = KDENSI  + N2BAST
      KFOCK   = KDNSPKI + NNBST(ISYMOP)
      KEMAT1  = KFOCK   + N2BST(ISYMOP)
      KEMAT2  = KEMAT1  + NE1TOT
      KEND1   = KEMAT2  + NE2TOT
      LWRK1   = LWORK   - KEND1
C
C----------------------------------------------------------------
C     Extra Work space allocation for C1 dependent intermediates.
C----------------------------------------------------------------
C
      KLAMPC  = KEND1
      KLAMHC  = KLAMPC  + NGLMDT(ISYMR)
      KDENSC  = KLAMHC  + NGLMDT(ISYMR)
      KDNSPKC = KDENSC  + N2BST(ISYMR)
      KFOCKC  = KDNSPKC + NNBST(ISYMR)
      KFOCKT  = KFOCKC  + N2BST(ISYMR)
      KFCKLR  = KFOCKT  + MAX(N2BST(ISYMR),N2BST(ISYMOP))
      KL1AOC  = KFCKLR  + N2BST(ISYRES)
      KRHO1AO = KL1AOC  + NGLMDT(ISYRES)
      KEND1   = KRHO1AO + NT1AO(ISYRES)
      LWRK1   = LWORK   - KEND1
C
      NRHOR  = NT2AOIJ(ISYMR)
      NMGD   = NT2AOIJ(ISYRES)
C
      IF (.NOT.((CC2.AND.(.NOT.NONHF)).OR.CCS)) THEN
         KXMAT  = KEND1
         KXMATX = KXMAT  + NMATIJ(ISYML)
         KYMAT  = KXMATX + NMATIJ(ISYRES)
         KYMATX = KYMAT  + NMATAB(ISYML)
         KYDENX = KYMATX + NMATAB(ISYRES)
         KYDENB = KYDENX + N2BST(ISYRES)
         KMINT  = KYDENB + N2BST(ISYRES)
         KMINTX = KMINT  + N3ORHF(ISYML)
         KCHIM  = KMINTX + N3ORHF(ISYRES)
         KRHOR  = KCHIM  + NMATIJ(ISYRES)
         KMGDL  = KRHOR  + NRHOR
         KEND1  = KMGDL  + NMGD
         LWRK1  = LWORK  - KEND1
      ENDIF
C
      IF ( IPRINT .GT. 25)
     *   WRITE(LUPRI,*) ' In CC_FMAT 1. alloc. work. left:',LWRK1
C
      IF (LWRK1 .LE. 0) CALL QUIT('Too little workspace in CC_FMAT')
C
C---------------------------------------
C     Initialize t1am and intermediates.
C---------------------------------------
C
C     CALL DZERO(WORK(KT1AM),NT1AM(1))
      CALL DZERO(WORK(KEMAT1),NE1TOT)
      CALL DZERO(WORK(KEMAT2),NE2TOT)
      CALL DZERO(WORK(KDENSI),N2BST(ISYMOP))
      CALL DZERO(WORK(KFOCK),N2BST(ISYMOP))
C
C-------------------------------------------
C     Initialize C1 dependent intermediates.
C-------------------------------------------
C
      CALL DZERO(WORK(KL1R2),N2BST(ISYRES))
      CALL DZERO(WORK(KFCKLR),N2BST(ISYRES))
      CALL DZERO(WORK(KFOCKC),N2BST(ISYMR))
      CALL DZERO(WORK(KDENSC),N2BST(ISYMR))
      CALL DZERO(WORK(KRHO1AO),NT1AO(ISYRES))
C
C----------------------------------
C     Calculate the lamda matrices:
C----------------------------------
C
      DTIME = SECOND()
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     *            WORK(KEND1),LWRK1)
 
      CALL CCLR_LAMTRA(WORK(KLAMDP),WORK(KLAMPC),
     *                 WORK(KLAMDH),WORK(KLAMHC),WORK(KC1AM),ISYMR)
 
      DTIME = SECOND() - DTIME
      TIMLAM  = DTIME + TIMLAM
 
      IF (IPRINT .GT.45) THEN
         CALL AROUND('Usual Lambda matrices ')
         CALL CC_PRLAM(WORK(KLAMDP),WORK(KLAMDH),ISYM0)
         CALL AROUND('C1 transformed Lambda matrices ')
         CALL CC_PRLAM(WORK(KLAMPC),WORK(KLAMHC),ISYMR)
      ENDIF
C
      IF (.NOT.(CCS.OR.(CC2.AND.(.NOT.NONHF)))) THEN
C
         IF ( LLIST(1:2) .EQ. 'L0' ) THEN
C
C           ---------------------------------------------------
C           Read X, Y and M intermediates from
C           zeroth-order Zeta and T2AMP0 from file
C           ---------------------------------------------------
C
            DTIME = SECOND()
            LUXYM = -1
            CALL GPOPEN(LUXYM,FILXYM,'OLD',' ','UNFORMATTED',IDUMMY,
     &                  .FALSE.)

            REWIND(LUXYM)

            IF (.NOT.(CCS.OR.(CC2.AND.(.NOT.NONHF)))) THEN
               READ(LUXYM,ERR=999) (WORK(KYMAT-1+I),I=1,NMATAB(ISYML))
               READ(LUXYM,ERR=999) (WORK(KXMAT-1+I),I=1,NMATIJ(ISYML))
            END IF

            IF (.NOT.(CCS.OR.(CC2.AND.(.NOT.NONHF)))) THEN
               READ(LUXYM,ERR=999) (WORK(KMINT-1+I),I=1,N3ORHF(ISYML))
            END IF
            
            CALL GPCLOSE(LUXYM,'KEEP')
            DTIME = SECOND() - DTIME
            TIMIO = TIMIO + DTIME

         ELSE
C
C           ---------------------------------------------------
C           Calculate X, Y and M intermediates, from
C           response Zeta and T2AMP0 
C           ---------------------------------------------------
C
            TIMEYI = SECOND()
            CALL CC_YI(WORK(KYMAT),WORK(KL2AM),ISYML,WORK(KRHO2),ISYM0,
     &                 WORK(KEND1),LWRK1)
            TIMEYI = SECOND() - TIMEYI
 
            TIMEXI = SECOND()
            CALL CC_XI(WORK(KXMAT),WORK(KL2AM),ISYML,WORK(KRHO2),ISYM0,
     &                 WORK(KEND1),LWRK1)
            TIMEXI = SECOND() - TIMEXI
 
            IF (.NOT.(CCS.OR.(CC2.AND.(.NOT.NONHF)))) THEN
               TIMEMI = SECOND()
               CALL CC_MI(WORK(KMINT),WORK(KL2AM),ISYML,WORK(KRHO2),
     &                    ISYM0,WORK(KEND1),LWRK1)
               TIMEMI = SECOND() - TIMEMI
            ENDIF

         END IF
C
C        -----------------------------------
C        Calculate X, Y and M intermediates:
C        from Zeta and response T2AMP
C        -----------------------------------
C
         DTIME = SECOND()
         CALL CC_YI(WORK(KYMATX),WORK(KL2AM),ISYML,
     *              WORK(KC2AM),ISYMR,WORK(KEND1),LWRK1)
         DTIME = SECOND() - DTIME
         TIMEYI = TIMEYI + DTIME
C
         DTIME = SECOND()
         CALL CC_XI(WORK(KXMATX),WORK(KL2AM),ISYML,
     *              WORK(KC2AM),ISYMR,WORK(KEND1),LWRK1)
         DTIME = SECOND() - DTIME
         TIMEXI = TIMEXI + DTIME
C
         IF (.NOT.(CCS.OR.(CC2.AND.(.NOT.NONHF)))) THEN
            DTIME = SECOND()
            CALL CC_MI(WORK(KMINTX),WORK(KL2AM),ISYML,
     *                 WORK(KC2AM),ISYMR,WORK(KEND1),LWRK1)
            DTIME = SECOND() - DTIME
            TIMEMI = TIMEMI + DTIME
         ENDIF
C 
C        ----------------------
C        calculate Y densities:
C        ----------------------
C
         DTIME = SECOND()
         CALL CC_YD(WORK(KYDENB),WORK(KYMAT),ISYML,WORK(KLAMDH),
     *              WORK(KLAMPC),ISYMR,WORK(KEND1),LWRK1)
C
         CALL CC_YD(WORK(KYDENX),WORK(KYMATX),ISYRES,WORK(KLAMDH),
     *              WORK(KLAMDP),ISYMOP,WORK(KEND1),LWRK1)
         DTIME = SECOND() - DTIME
         TIMEYI = TIMEYI + DTIME
C 
C        ------------------------------
C        calculate response Chi matrix:
C        (not backtransformed to AO)
C        ------------------------------
C
         DO ISYMK = 1, NSYM
           ISYMC = MULD2H(ISYMK,ISYMR)
           ISYMI = MULD2H(ISYMC,ISYML)

           KOFF1 = KC1AM + IT1AM(ISYMC,ISYMK) 
           KOFF2 = KL1AM + IT1AM(ISYMC,ISYMI) 
           KOFF3 = KCHIM + IMATIJ(ISYMK,ISYMI)

           NRHFK = MAX(NRHF(ISYMK),1)
           NVIRC = MAX(NVIR(ISYMC),1)

           CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMI),NVIR(ISYMC),
     *                -ONE,WORK(KOFF1),NVIRC,WORK(KOFF2),NVIRC,
     *                ZERO,WORK(KOFF3),NRHFK)

         END DO
         CALL DAXPY(NMATIJ(ISYRES),-ONE,WORK(KXMATX),1,WORK(KCHIM),1)
C
C        -------------------------------------------------------
C        precalculate effective density for left B intermediate:
C        -------------------------------------------------------
C
         IF (CCSD) THEN

            DTIME = SECOND()
            CALL CC_BFDENF(WORK(KL2AM),ISYML,WORK(KMINTX),ISYRES,
     *                     WORK(KLAMDP),ISYMOP,WORK(KCHIM),ISYRES,
     *                     WORK(KC1AM),ISYMR,WORK(KMGDL),
     *                     WORK(KEND1),LWRK1)
            TIMBF = TIMBF + SECOND() - DTIME

         END IF

      ENDIF
C
C--------------------------------
C     Calculate L1R2 contraction:
C--------------------------------
C
      IF (DEBUG) THEN
         RHO1N = DDOT(NT1AM(ISYML),WORK(KL1AM),1,WORK(KL1AM ),1)
         WRITE(LUPRI,1) 'Norm of L1AM before  CC_C1T2C:      ',RHO1N
         RHO1N = DDOT(NT2AM(ISYMR),RHO2,1,RHO2,1)
         WRITE(LUPRI,1) 'Norm of c2am before C1T2C:         ',RHO1N
         RHO1N = DDOT(NGLMDT(ISYMOP),WORK(KLAMDH),1,WORK(KLAMDH),1)
         WRITE(LUPRI,1) 'Norm of lamdaH before C1T2C:       ',RHO1N
         RHO1N = DDOT(NGLMDT(ISYMR ),WORK(KLAMPC),1,WORK(KLAMPC),1)
         WRITE(LUPRI,1) 'Norm of lamdPC before C1T2C:       ',RHO1N
         RHO1N = DDOT(NGLMDT(ISYMR ),WORK(KLAMHC),1,WORK(KLAMHC),1)
         WRITE(LUPRI,1) 'Norm of lamdHC before C1T2C:       ',RHO1N
         RHO1N = DDOT(NGLMDT(ISYMOP),WORK(KLAMDP),1,WORK(KLAMDP),1)
         WRITE(LUPRI,1) 'Norm of lamdaP before C1T2C:       ',RHO1N
      ENDIF
      IOPT = 3
      TIMC1T2 = SECOND()
CTST
C     SINGLETST = .TRUE.
C     IF (SINGLETST) THEN
C        CALL DZERO(WORK(KC2AM),NT2AM(ISYMR))
C     END IF
CTST
      CALL CC_C1T2C(WORK(KL1R2),WORK(KL1AM),WORK(KC2AM),
     *              WORK(KLAMDP),WORK(KLAMDH),
     *              WORK(KLAMPC),WORK(KLAMHC),
     *              WORK(KEND1),LWRK1,ISYML,ISYMR,IOPT)
      TIMC1T2 = SECOND() - TIMC1T2
      IF (DEBUG) THEN
         RHO2N = DDOT(NT2AM(ISYMR),WORK(KC2AM),1,WORK(KC2AM),1)
         RHO1N = DDOT(N2BST(ISYRES),WORK(KL1R2),1,WORK(KL1R2),1)
         WRITE(LUPRI,1) 'Norm of RHO2  after  CC_C1T2C:      ',RHO2N
         WRITE(LUPRI,1) 'Norm of FCKR2 after  CC_C1T2C:      ',RHO1N
      ENDIF
C
C     ------------------------------------------------------------
C     Include the two LT21A contributions here by 2N^2 operations:
C     ------------------------------------------------------------
C
      IF (.NOT.(CC2.OR.CCS)) THEN

C       SINGLETST = .FALSE.
C       IF (SINGLETST) THEN
          CALL DAXPY(N2BST(ISYRES),ONE,WORK(KYDENX),1,WORK(KL1R2),1)
C       END IF

          CALL DAXPY(N2BST(ISYRES),ONE,WORK(KYDENB),1,WORK(KL1R2),1)

      END IF

      CALL CC_DNSPK(WORK(KL1R2),WORK(KL1R2PK),ISYRES)
C
C------------------------------------------------------------------
C     Calculate the density matrices:
C     IC=1 include core contribution for zeroth-order matrix DENSI
C     IC=0 no core contribution for C1 transformed matrix DENSC
C------------------------------------------------------------------
C
      DTIME  = SECOND()
C
      IC     = 1
      CALL CC_AODENS(WORK(KLAMDP),WORK(KLAMDH),WORK(KDENSI),ISYM0,
     *               IC,WORK(KEND1),LWRK1)
      CALL CC_DNSPK(WORK(KDENSI),WORK(KDNSPKI),ISYM0)
C
      IC     = 0
      CALL CC_AODENS(WORK(KLAMDP),WORK(KLAMHC),WORK(KDENSC),ISYMR,
     *               IC,WORK(KEND1),LWRK1)
      CALL CC_DNSPK(WORK(KDENSC),WORK(KDNSPKC),ISYMR)
C
      DTIME  = SECOND() - DTIME
      TIMAOD = TIMAOD + DTIME
C
      IF (IPRINT .GT. 45) THEN
         CALL AROUND('CC_FMAT : Usual Lamda density matrix')
         CALL CC_PRAODEN(WORK(KDENSI),ISYM0)
         CALL AROUND('CC_FMAT : C1 trans. Lamda density matrix')
         CALL CC_PRAODEN(WORK(KDENSC),ISYMR)
      ENDIF
C
C--------------------------------------------
C     Calculate L1LamdahaC1 Lambda vector.
C--------------------------------------------
C
      CALL CCLL_LAMTRA(WORK(KL1AM),ISYML,WORK(KLAMPC),
     *                 ISYMR,WORK(KL1AOC))
C
      IF (IPRINT .GT. 45) THEN
         CALL AROUND('L1 transformed Lambda C1 P matrix ')
         CALL CC_PRLAM(WORK(KL1AOC),WORK(KLAMHC),ISYRES)
      ENDIF
C
C--------------------------------------------------------------
C     Prepare C2 vector and its transposed counterpart in core.
C--------------------------------------------------------------
C
      IF ( .NOT. CCS ) THEN

         IOPT = 2
         CALL CC_RDRSP(RLIST,IRLNR,ISYMR,IOPT,MODELX,
     *                 DUMMY,WORK(KRHO2))
         CALL CCLR_DIASCL(WORK(KRHO2),TWO,ISYMR)
         CALL CC_T2SQ(WORK(KRHO2),WORK(KC2AM),ISYMR)

         IF (T2TCOR) THEN
            DTIME = SECOND()
            CALL DCOPY(NT2SQ(ISYMR),WORK(KC2AM),1,WORK(K2C2C2),1)
            CALL CCSD_T2TP(WORK(K2C2C2),WORK(KEND1),LWRK1,ISYMR)
            DTIME = SECOND() - DTIME
            TIMT2TR = TIMT2TR + DTIME
         ENDIF

      ENDIF
C
C     --------------------------------------------------------
C     precalculate effective density for right B intermediate:
C     --------------------------------------------------------
C
      IF (CCSD) THEN

          IVEC  = 1
          IOPT  = 3
          DTIME = SECOND()
          CALL CC_BFDEN(WORK(KC2AM),ISYMR,DUMMY,IDUMMY,
     *                  WORK(KLAMDH),ISYMOP,WORK(KLAMDH),ISYMOP,
     *                  WORK(KLAMHC),ISYMR,DUMMY,IDUMMY,
     *                  FNBFD,LUBFD,IADRBFD,ISTARTBFD,
     *                  IVEC, IOPT, .FALSE., WORK(KEND1), LWRK1)
          TIMBF = TIMBF + SECOND() - DTIME

      END IF
C
C------------------------------------------------
C     Read one-electron integrals in Fock-matrix.
C------------------------------------------------
C
      TIMFCK = SECOND()
      CALL CCRHS_ONEAO(WORK(KFOCK),WORK(KEND1),LWRK1)
      DO IF = 1, NFIELD
         FF = EFIELD(IF)
         CALL CC_ONEP(WORK(KFOCK),WORK(KEND1),LWRK1,FF,1,LFIELD(IF))
      END DO
      TIMFCK = SECOND() - TIMFCK
C
C-------------------------------------------
C     initialize B intermediates with zeros:
C-------------------------------------------
C
      IF (.NOT. CCS) THEN
        CALL DZERO(WORK(KRHO2),NRHO2)
        CALL CC_WVEC(LUTR,TRFIL,NRHO2,NRHO2,1,WORK(KRHO2))
      ENDIF
      IF (.NOT. (CC2.OR.CCS)) THEN
         CALL DZERO(WORK(KRHOR),NRHOR)
      ENDIF
C
C====================================================
C     Start the loop over distributions of integrals.
C====================================================
C
      KENDS2 = KEND1
      LWRKS2 = LWRK1
C
      IF (DIRECT) THEN
         NTOSYM = 1
         DTIME  = SECOND()
C
         IF (HERDIR) THEN
           CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
         ELSE
           KCCFB1 = KEND1
           KINDXB = KCCFB1 + MXPRIM*MXCONT
           KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
           LWRK1  = LWORK  - KEND1
C
           CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     &                 KODPP1,KODPP2,KRDPP1,KRDPP2,
     &                 KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     &                 WORK(KEND1),LWRK1,IPRERI)
         END IF

         KEND1 = KFREE
         LWRK1 = LFREE
         DTIME  = SECOND() - DTIME
         TIMER1 = TIMER1 + DTIME
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      KENDSV = KEND1
      LWRKSV = LWRK1
C
      ICDEL1 = 0
      ICDEL2 = 0
C
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
C------------------------------------------
C        If direct calculate the integrals.
C------------------------------------------
C
            IF (DIRECT) THEN
C
               KEND1 = KENDSV
               LWRK1 = LWRKSV
C
               DTIME  = SECOND()
               IF (HERDIR) THEN
                 CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                       IPRINT)                             
               ELSE
                 CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     &                       WORK(KODCL1),WORK(KODCL2),
     &                       WORK(KODBC1),WORK(KODBC2),
     &                       WORK(KRDBC1),WORK(KRDBC2),
     &                       WORK(KODPP1),WORK(KODPP2),
     &                       WORK(KRDPP1),WORK(KRDPP2),
     &                       WORK(KCCFB1),WORK(KINDXB),
     &                       WORK(KEND1), LWRK1,IPRERI)
               ENDIF
               DTIME  = SECOND() - DTIME
               TIMER2 = TIMER2 + DTIME
C
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient core in CC_FMAT')
               END IF
C
            ELSE
               NUMDIS = 1
            ENDIF
C
C--------------------------------------------------
C           Loop over number of distributions in disk.
C--------------------------------------------------
C
            DO 130 IDEL2 = 1,NUMDIS
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
               ISYAIK = MULD2H(ISYDIS,ISYMR)
C
C----------------------------------------------------------
C              Calculate adresses for c,cio,d,dio routines.
C----------------------------------------------------------
C
               IT2DEL(IDEL) = ICDEL1
               ICDEL1 = ICDEL1 + NT2BCD(ISYDIS)
C
               IT2DLR(IDEL,1) = ICDEL2
               ICDEL2 = ICDEL2 + NT2BCD(ISYAIK)
C
C------------------------------------------
C              Work space allocation no. 2.
C------------------------------------------
C
               KXINT  = KEND1
               KEND2  = KXINT + NDISAO(ISYDIS)
               LWRK2  = LWORK - KEND2
C
               IF ( IPRINT .GT. 55)
     *            WRITE(LUPRI,*) ' In CC_FMAT 2. alloc. work left:',
     *            LWRK2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CC_FMAT-2')
               ENDIF
C
C-----------------------------------------
C              Read in batch of integrals.
C-----------------------------------------
C
               DTIME   = SECOND()
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
               DTIME   = SECOND() - DTIME
               TIMRDAO = TIMRDAO  + DTIME
C
C------------------------------------------------------------------
C              Calculate an AO-Fock matrix with C1. trans. density.
C------------------------------------------------------------------
C
               ISYDEN  = ISYMR
               DTIME   = SECOND()
               IOPT    = 1
               CALL CC_AOFOCK2(WORK(KXINT),WORK(KDENSC),WORK(KDNSPKC),
     *                         WORK(KFOCKC),WORK(KEND2),LWRK2,
     *                         IDEL,ISYDIS,ISYMD,ISYDEN,IOPT)
               DTIME  = SECOND() - DTIME
               TIMFCK = TIMFCK + DTIME
C
C---------------------------------------------------------
C              Calculate AO-Fock matrix with L1R2 density.
C---------------------------------------------------------
C
               ISYDEN  = ISYRES 
               DTIME   = SECOND()
               IOPT    = 1
               CALL CC_AOFOCK2(WORK(KXINT),WORK(KL1R2),WORK(KL1R2PK),
     *                         WORK(KFCKLR),WORK(KEND2),LWRK2,
     *                         IDEL,ISYDIS,ISYMD,ISYDEN,IOPT)
               DTIME  = SECOND() - DTIME
               TIMFCK = TIMFCK + DTIME
C
C-------------------------------------------------------------------
C              IF CCS calculate fock matrix since it is not on file.
C              IF CCS jump to end of loop.
C-------------------------------------------------------------------
C
               IF ( CCS ) THEN
                  ISYDEN  = 1
                  DTIME   = SECOND()
                  IOPT    = 1
                  CALL CC_AOFOCK2(WORK(KXINT),WORK(KDENSI),
     *                            WORK(KDNSPKI),WORK(KFOCK),
     *                            WORK(KEND2),LWRK2,
     *                            IDEL,ISYDIS,ISYMD,ISYDEN,IOPT)
                  DTIME  = SECOND() - DTIME
                  TIMFCK = TIMFCK + DTIME
               ENDIF
C
C-------------------------------
C              IF CCS then JUMP.
C-------------------------------
C
               IF (CCS) GOTO 130
C
C-------------------------------------------------------
C              Work space allocation no. 3.
C-------------------------------------------------------
C
               ISAIJK = MULD2H(ISYMD,ISYMR)
               IAIJK2 = MULD2H(ISYMD,ISYMOP)
C
               KDSRHF = KEND2
               K3OINT = KDSRHF + NDSRHF(ISYMD)
               KO3T   = K3OINT + NMAIJK(ISAIJK)
               KEND3  = KO3T   + NMAIJK(IAIJK2)
C
               LWRK3  = LWORK  - KEND3
               IF ( IPRINT .GT. 55)
     *            WRITE(LUPRI,*) ' In CC_FMAT 3. alloc. work left:',
     *            LWRK3
C
               IF (LWRK3 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND3,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CC_FMAT-3 ')
               ENDIF
C
C--------------------------------------------------------
C              Transform one index in the integral batch.
C--------------------------------------------------------
C
               DTIME  = SECOND()
               CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KLAMDP),
     *                     ISYM0,WORK(KEND3),LWRK3,ISYDIS)
               TIMTRBT = TIMTRBT + SECOND() - DTIME
C
C-------------------------------------------------------------------
C              calculate the BF intermediates:
C-------------------------------------------------------------------
C
               IF ( CCSD ) THEN
                  ISYMGDR = MULD2H(ISYMD,ISYMR)
                  NMGDR   = NT2BGD(ISYMGDR)
C
                  KMGDR  = KEND3
                  KBFRHF = KMGDR  + NMGDR
                  KEND4  = KBFRHF + NBFRHF(ISYMD)
                  LWRK4  = LWORK  - KEND4
C
                  IF (LWRK4 .LT. 0) THEN
                     CALL QUIT('Insufficient space in CC_FMAT-4 ')
                  ENDIF
C
                  IADR  = IADRBFD(IDEL)
                  CALL GETWA2(LUBFD,FNBFD,WORK(KMGDR),IADR,NMGDR)
C
                  CALL CC_BFBSORT(WORK(KDSRHF),WORK(KBFRHF),ISYDIS)
C
                  CALL CC_BFIB(WORK(KRHOR),WORK(KBFRHF),ISYDIS,
     *                         WORK(KMGDR),ISYMGDR,WORK(KEND4),LWRK4)
C
                  CALL CC_BFIF(WORK(KBFRHF),ISYDIS,WORK(KMGDL),ISYRES,
     *                         LUBFI,FNBFI,IADRBFI,ISTARTBFI,IDEL,
     *                         WORK(KEND4),LWRK4)

                 XNORM = DDOT(NMGDR,WORK(KMGDR),1,WORK(KMGDR),1)
                 WRITE (LUPRI,*) 'IDEL,Norm(MGD) = ',IDEL,XNORM
                 XNORM = DDOT(NT2AOIJ(ISYMR),WORK(KRHOR),1,
     *                                       WORK(KRHOR),1)
                 WRITE (LUPRI,*) 'IDEL,Norm(BF) = ',IDEL,XNORM
               END IF
C
C-------------------------------------------------------------------
C              Calculate integral batch with three occupied indices.
C-------------------------------------------------------------------
C
               DTIME = SECOND()
               CALL CC_INT3O(WORK(K3OINT),WORK(KDSRHF),WORK(KLAMHC),
     *                       ISYMR,WORK(KLAMDP),WORK(KEND3),
     *                       LWRK3,IDEL,ISYMD,LUO3,O3FIL)
               TIME3O = TIME3O + SECOND() - DTIME
C
C-------------------------------------------------------------------
C              Calculate integral batch with three occupied indices.
C              but C1 independent.
C-------------------------------------------------------------------
C
               IF (.NOT.(CC2.OR.CCS)) THEN
                  DTIME = SECOND()
                  CALL CC_INT3O(WORK(KO3T),WORK(KDSRHF),WORK(KLAMDH),
     *                          ISYMOP,WORK(KLAMDP),WORK(KEND3),
     *                          LWRK3,IDEL,ISYMD,LUO3T,O3TFIL)
                  TIME3O = TIME3O + SECOND() - DTIME
               ENDIF
C
C--------------------------------------------------------------
C              Calculate intermediates needed for the 21-block.
C              This section is for both perturbed an unperturbed.
C--------------------------------------------------------------
C
               IF (.NOT. (CC2.OR.CCS)) THEN
C
C                 ---------------------------------------------------
C                 calculate first 21I and 21H terms with zeroth-order
C                 T2 amplitudes and first-order integrals (or LamdaP)
C                 ---------------------------------------------------
C
                  ISYMTD = MULD2H(ISYMOP,ISYMD)
                  ISYLTD = MULD2H(ISYMTD,ISYML)
C
                  LEN0   = NT2BCD(ISYLTD)
C
                  ISYMRD = MULD2H(ISYMR,ISYMD)
                  ISYLRD = MULD2H(ISYMRD,ISYML)
C
                  LEN1   = NT2BCD(ISYLRD)
C    
                  KSCTZI = KEND3
                  KSCTVI = KSCTZI + NT2BCD(ISYLTD)
                  KSCRZI = KSCTVI + NT2BCD(ISYLTD)
                  KSCRVI = KSCRZI + NT2BCD(ISYLRD)
                  KSCRWI = KSCRVI + NT2BCD(ISYLRD)
                  KSCTWI = KSCRWI + NT2BCD(ISYLRD)
                  KEND4  = KSCTWI + NT2BCD(ISYLTD)
                  LWRK4  = LWORK  - KEND4
                  IF (LWRK4 .LT. 0) THEN
                     CALL QUIT('Insufficient memory in CC_FMAT (4).')
                  ENDIF
C
                  CALL DZERO(WORK(KSCRWI),NT2BCD(ISYLRD))
                  CALL DZERO(WORK(KSCTWI),NT2BCD(ISYLTD))
C
C                 ---------------------------------------
C                 read P0 and Q0 intermediates from file:
C                 ---------------------------------------
C
                  DTIME = SECOND()
                  CALL GETWA2(LUPQR0,FILPQR0,WORK(KSCTZI),
     &                        IADRPQ0(IDEL),NT2BCD(ISYLTD))

                  CALL GETWA2(LUPQR0,FILPQR0,WORK(KSCTVI),
     &                   IADRPQ0(IDEL)+NT2BCD(ISYLTD),NT2BCD(ISYLTD))
                  TIMIO = TIMIO + SECOND() - DTIME
C
C                 ---------------------------------------
C                 read PR and QR intermediates from file:
C                 ---------------------------------------
C
                  DTIME  = SECOND()
                  CALL GETWA2(LUPQRR,FILPQRR,WORK(KSCRZI),
     &                        IADRPQR(IDEL),NT2BCD(ISYLRD))

                  CALL GETWA2(LUPQRR,FILPQRR,WORK(KSCRVI),
     &                   IADRPQR(IDEL)+NT2BCD(ISYLRD),NT2BCD(ISYLRD))
                  TIMIO = TIMIO + SECOND() - DTIME
C
C                 ---------------------------------
C                 Calculate the LT21I contribution: 
C                 ---------------------------------
C
                  IOPT = 2
                  DTIME = SECOND()
                  CALL CC_21I2(WORK(KRHO1AO), 
     *                         WORK(KXINT),ISYDIS,DUMMY,0,
     *                         WORK(KSCTZI),WORK(KSCTVI),ISYLTD,
     *                         WORK(KSCRZI),WORK(KSCRVI),ISYLRD,
     *                         WORK(KLAMDP),WORK(KLAMDH),ISYMOP,
     *                         WORK(KLAMPC),ISYMR,
     *                         WORK(KEND4),LWRK4,IOPT,
     *                         .TRUE.,.FALSE.,.FALSE.)
                  TIMEI = TIMEI + SECOND() - DTIME 

         WRITE (LUPRI,*)'IDEL = ',IDEL,IADRPQ0(IDEL),IADRPQR(IDEL),
     &                 LEN0,LEN1
         WRITE (LUPRI,*) 'Norm(PINT0) = ',
     *             DDOT(LEN0,WORK(KSCTZI),1,WORK(KSCTZI),1)
         WRITE (LUPRI,*) 'Norm(QINT0) = ',
     *             DDOT(LEN0,WORK(KSCTVI),1,WORK(KSCTVI),1)
         WRITE (LUPRI,*) 'Norm(PINT1) = ',
     *             DDOT(LEN1,WORK(KSCRZI),1,WORK(KSCRZI),1)
         WRITE (LUPRI,*) 'Norm(QINT1) = ',
     *             DDOT(LEN1,WORK(KSCRVI),1,WORK(KSCRVI),1)
         WRITE (LUPRI,*) 'Norm(WINT1) = ',
     *             DDOT(LEN1,WORK(KSCRWI),1,WORK(KSCRWI),1)
         WRITE (LUPRI,*) 'Norm(RHO1AO) = ',
     *              DDOT(NT1AO(ISYRES),WORK(KRHO1AO),1,WORK(KRHO1AO),1)
C
C                 ---------------------------------
C                 Calculate the LT21H contribution:
C                 ---------------------------------
C
                  ISYVWZ = ISYLTD 
                  ISYINT = ISYMR
                  DTIME  = SECOND()
                  CALL CC_21H(WORK(KRHO1),ISYRES,WORK(KSCTVI),
     *                        WORK(KSCTWI),WORK(KSCTZI),ISYVWZ,
     *                        WORK(K3OINT),ISYINT,
     *                        WORK(KEND4),LWRK4,ISYMD,NEWZWV)
                  TIMEH = TIMEH + SECOND() - DTIME
         WRITE (LUPRI,*) 'Norm(RHO1) = ',
     *              DDOT(NT1AM(ISYRES),WORK(KRHO1),1,WORK(KRHO1),1)
C
C                 ----------------------------------------------------
C                 Calculate the LT21H contribution with first-order
C                 T2 amplitudes and zeroth-order integrals 
C                 ----------------------------------------------------
C
                  ISYVWZ = ISYLRD 
                  ISYINT = ISYM0
                  DTIME  = SECOND()
                  CALL CC_21H(WORK(KRHO1),ISYRES,WORK(KSCRVI),
     *                        WORK(KSCRWI),WORK(KSCRZI),ISYVWZ,
     *                        WORK(KO3T),ISYINT,
     *                        WORK(KEND4),LWRK4,ISYMD,NEWZWV)
                  TIMEH = TIMEH + SECOND() - DTIME
         WRITE (LUPRI,*) 'Norm(RHO1) = ',
     *              DDOT(NT1AM(ISYRES),WORK(KRHO1),1,WORK(KRHO1),1)

               ENDIF
C
C----------------------------------------------------------
C              Calculate the LT12B term in CC2-calculation.
C              rho(ai,bj)=P(ai,bj)(sum(c)(l(c,j)*~cbia))
C----------------------------------------------------------
C
               IF (CC2) THEN
 
                  DTIME = SECOND()
                  CALL CC_12B(WORK(KRHO2),WORK(KDSRHF),WORK(KL1AOC),
     *                        ISYRES,WORK(KLAMDH),WORK(KEND3),LWRK3,
     *                        IDEL,ISYMD)
                  TIM12B = TIM12B + SECOND() - DTIME
 
               ENDIF
C
C--------------------------------------------------------------------
C              Construct the partially transformed L2-amplitudes
C--------------------------------------------------------------------
C
               ISYDVI  = MULD2H(ISYMD,ISYMR)
               KSCRLM  = KEND3
               KSCRLC  = KSCRLM + NT2MMO(ISYMD,ISYML)
               KEND41  = KSCRLC + NT2MMO(ISYDVI,ISYML)
               LWRK41  = LWORK  - KEND41
C
               DTIME  = SECOND()
               ICON   = 1
               ISYMLP = 1
               CALL CC_T2AO(WORK(KL2AM),WORK(KLAMDP),ISYMLP,
     *                      WORK(KSCRLM),WORK(KEND41),LWRK41,
     *                      IDEL,ISYMD,ISYML,ICON)
               TIMT2AO = TIMT2AO + SECOND() - DTIME
C
               DTIME  = SECOND()
               ICON   = 1
               ISYMLP = ISYMR
               CALL CC_T2AO(WORK(KL2AM),WORK(KLAMPC),ISYMLP,
     *                      WORK(KSCRLC),WORK(KEND41),LWRK41,
     *                      IDEL,ISYMD,ISYML,ICON)
               TIMT2AO = TIMT2AO + SECOND() - DTIME
C
C----------------------------------------------------------------
C              Calculate the 21C- and D contributions.
C              If not CC2 skip 21C calculation and calculate this
C              as part of the 21B terms.
C----------------------------------------------------------------
C
               DTIME = SECOND()
               IOPT  = 2
               IF ( CC2 ) THEN
                 ICON  = 2
               ELSE
                 ICON  = 1
               ENDIF

               ISYMM  = MULD2H(ISYMD,ISYML)
               ISYMMC = MULD2H(ISYDVI,ISYML)
               CALL CC_21DC(WORK(KRHO1),WORK(KL2AM),ISYML,
     *                      WORK(KSCRLM),ISYMM,
     *                      WORK(KSCRLC),ISYMMC,WORK(KXINT),
     *                      WORK(KLAMDH),1,WORK(KLAMPC),ISYMR,
     *                      WORK(KLAMHC),ISYMR,WORK(KLAMDP),1,
     *                      WORK(KEND41),LWRK41,IDEL,ISYMD,IOPT,ICON)
               TIMEC = TIMEC + SECOND() - DTIME
C
C--------------------------------------------------------------------
C              Construct the partially transformed C2-amplitudes, CM.
C--------------------------------------------------------------------
C
               KSCRCM = KEND41
               KEND4  = KSCRCM + NT2MMO(ISYMD,ISYMR)
               LWRK4  = LWORK  - KEND4
C
               DTIME  = SECOND()
               ICON   = 1
               ISYMLH = 1
               CALL CC_T2AO(WORK(KC2AM),WORK(KLAMDH),ISYMLH,
     *                      WORK(KSCRCM),WORK(KEND4),LWRK4,
     *                      IDEL,ISYMD,ISYMR,ICON)
               TIMT2AO = TIMT2AO +  SECOND() - DTIME
C
C---------------------------------------
C              Transform CM to 2CM - CM.
C---------------------------------------
C
               DTIME = SECOND()
               CALL CC_MTCME(WORK(KSCRCM),WORK(KEND4),LWRK4,
     *                       ISYMD,ISYMR)
               TIMTCME = TIMTCME  + SECOND() - DTIME
C
C-------------------------------------------------------
C              Calculate the C-tilde local intermediate.
C-------------------------------------------------------
C
               IF ( .NOT. (CC2.OR.CCS)) THEN
C
                  FACTC = XMONE 
                  ICON  = 3
                  IOPTR12 = 0
                  IOPTE = 0
C
                  DTIME   = SECOND()
                  IF (.NOT. T2TCOR) THEN
                     CALL CCRHS_C(WORK(KXINT),WORK(KDSRHF),
     *                            DUMMY, WORK(KC2AM),ISYMR,
     *                            WORK(KLAMDP),WORK(KLAMDP),
     *                            WORK(KLAMDH),WORK(KLAMPC),ISYMR,
     *                            WORK(KLAMHC),ISYMR,
     *                            DUMMY,DUMMY,WORK(KEND4),
     *                            LWRK4,IDEL,ISYMD,FACTC,ICON,IOPTR12,
     *                            IOPTE,LUC,CTFIL,IDUMMY,DUMMY,1)
                  ELSE
                     CALL CCRHS_C(WORK(KXINT),WORK(KDSRHF),
     *                            DUMMY, WORK(K2C2C2),ISYMR,
     *                            WORK(KLAMDP),WORK(KLAMDP),
     *                            WORK(KLAMDH),WORK(KLAMPC),ISYMR,
     *                            WORK(KLAMHC),ISYMR,
     *                            DUMMY,DUMMY,WORK(KEND4),
     *                            LWRK4,IDEL,ISYMD,FACTC,ICON,IOPTR12,
     *                            IOPTE,LUC,CTFIL,DUMMY,DUMMY,1)
                  ENDIF
                  TIMC    = TIMC + SECOND() - DTIME

               ENDIF
C
C---------------------------------------
C              Transform C2 to 2C2 - C2.
C---------------------------------------
C
               DTIME   = SECOND()
               IF (T2TCOR) THEN
                  CALL DSCAL(NT2SQ(ISYMR),TWO,WORK(KC2AM),1)
                  CALL DAXPY(NT2SQ(ISYMR),-ONE,WORK(K2C2C2),1,
     *                                         WORK(KC2AM),1)
               ELSE 
                  CALL CCRHS_T2TR(WORK(KC2AM),WORK(KEND4),LWRK4,ISYMR)
               ENDIF
               TIMT2TR = TIMT2TR  + SECOND() - DTIME
C
C-------------------------------------------------------
C              Calculate the D-tilde local intermediate.
C-------------------------------------------------------
C
               IF ( .NOT. (CC2.OR.CCS)) THEN
C
                  ICON   = 3
                  FACTD  = 1.0D00
                  IOPTR12 = 0
                  IOPTE  = 0
C
                  DTIME   = SECOND()
                  CALL CCRHS_D(WORK(KXINT),WORK(KDSRHF),DUMMY,
     *                         WORK(KC2AM),ISYMR,
     *                         WORK(KLAMDP),DUMMY,WORK(KLAMDH),
     *                         WORK(KLAMPC),ISYMR,WORK(KLAMHC),ISYMR,
     *                         DUMMY,DUMMY,WORK(KEND4),LWRK4,IDEL,
     *                         ISYMD,FACTD,ICON,IOPTR12,IOPTE,
     *                         LUD,DTFIL,IDUMMY,DUMMY,1)
                  TIMD    = TIMD     + SECOND() - DTIME
C
               ENDIF
C
C----------------------------------------
C              Calculate E-intermediates.
C----------------------------------------
C
               IF (.NOT.CCS) THEN
                  DTIME   = SECOND()
                  CALL CCRHS_EI(WORK(KDSRHF),WORK(KEMAT1),
     *                          WORK(KEMAT2),WORK(KC2AM),
     *                          WORK(KSCRCM),WORK(KLAMDP),
     *                          WORK(KLAMDH),WORK(KEND4),LWRK4,
     *                          IDEL,ISYMD,ISYDIS,ISYMR)
               ENDIF
               TIMEI   = TIMEI    + SECOND() - DTIME
C
C---------------------------------------------
C              BackTransform C2 from 2C2 - C2.
C---------------------------------------------
C
               DTIME   = SECOND()
               IF (T2TCOR) THEN
                  CALL DAXPY(NT2SQ(ISYMR),ONE,WORK(K2C2C2),1,
     *                       WORK(KC2AM),1)
                  CALL DSCAL(NT2SQ(ISYMR),XHALF,WORK(KC2AM),1)
               ELSE
                  CALL CCRHS_T2BT(WORK(KC2AM),WORK(KEND4),LWRK4,ISYMR)
               END IF
               TIMT2BT = TIMT2BT  + SECOND() - DTIME
C
c              IF (DEBUG) THEN
c                 CALL AROUND('END of LOOP  CCLR:RHO ')
c                 CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,1,1)
c              ENDIF
C
               IF (IPRINT .GT. 20 .OR. DEBUG) THEN
                  RHO1N=DDOT(NT1AM(ISYRES),WORK(KRHO1),1,WORK(KRHO1),1)
                  IF (OMEGOR) THEN
                    RHO2N = DDOT(2*NT2ORT(ISYRES),WORK(KRHO2),1,
     *                           WORK(KRHO2),1)
                  ELSE
                    RHO2N = DDOT(NT2AM(ISYRES),WORK(KRHO2),1,
     *                           WORK(KRHO2),1)
                  ENDIF
                  WRITE (LUPRI,1)'Norm of Rho1 - end of loop:        ',
     &                 RHO1N
                  WRITE (LUPRI,1)'Norm of Rho2 - end of loop:        ',
     &                 RHO2N
               ENDIF
C
  130       CONTINUE
C
C--------------------------------------
C              write out result vector.
C--------------------------------------
C
               IF (.NOT. CCS) THEN
                  DTIME   = SECOND()
                  CALL CC_WVEC(LUTR,TRFIL,NRHO2,NRHO2,1,WORK(KRHO2))
                  DTIME   = SECOND() - DTIME
                  TIMIO   = TIMIO    + DTIME
               ENDIF
C
  110    CONTINUE
  100 CONTINUE
C
      WRITE (LUPRI,*) 'after loop over AO WORK(0) = ',WORK(NULL)
      CALL AROUND('After LOOP over AO integrals:RHO ')
      CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,1,0)
C
C------------------------
C     Recover work space.
C------------------------
C
      KT2AMP = KENDS2
      KENDS2 = KT2AMP + MAX(NT2AM(ISYMOP),2*NT2ORT(ISYMOP))
      LWRKS2 = LWORK  - KENDS2

      KEND1 = KENDS2
      LWRK1 = LWRKS2
C
C--------------------------------------------
C     Allocate space for the gamma matrix.
C--------------------------------------------
C
      IF (NEWGAM) THEN
C
         KGAMMC = KEND1
         KEND1  = KGAMMC + NGAMMA(ISYMR)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) CALL QUIT('Insufficient memory in CC_FMAT')
C
      END IF
C
C----------------------------------
C     Usual Fock Matrix.
C     Save AO fock matric in fockt.
C----------------------------------
C
      ISYFAO = 1
      ISYMPA = 1
      ISYMHO = 1
C
      IF ( .NOT. CCS) THEN
         LFOCK = -1
         CALL GPOPEN(LFOCK,'CC_FCKH','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     *               .FALSE.)
         REWIND(LFOCK)
         READ(LFOCK) (WORK(KFOCK + I - 1) , I = 1,N2BST(ISYMOP))
         CALL GPCLOSE(LFOCK,'KEEP')
      ENDIF
C
      IF (IPRINT .GT.140) THEN
         CALL AROUND( 'Usual Fock AO matrix' )
         CALL CC_PRFCKAO(WORK(KFOCK),ISYFAO)
      ENDIF
C
      CALL DCOPY(N2BST(ISYMOP),WORK(KFOCK),1,WORK(KFOCKT),1)
C
      DTIME  = SECOND()
      CALL CC_FCKMO(WORK(KFOCK),WORK(KLAMDP),WORK(KLAMDH),
     *              WORK(KEND1),LWRK1,ISYFAO,ISYMPA,ISYMHO)
      DTIME  = SECOND() - DTIME
      TIMFCKMO = DTIME  +  TIMFCKMO
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND( 'Usual Fock MO matrix' )
         CALL CC_PRFCKMO(WORK(KFOCK),ISYFAO)
      ENDIF
C
C-----------------------------------------
C        Calculate the LT11A contribution.
C-----------------------------------------
C
      DTIME = SECOND()
      IF (DEBUG) THEN
         RHO1N = DDOT(N2BST(ISYRES),WORK(KFCKLR),1,WORK(KFCKLR),1)
         WRITE(LUPRI,1) 'Norm of FLR1 before CC_11A:        ',RHO1N
      ENDIF

      SINGLETST = .TRUE.
      IF (.NOT.SINGLETST) THEN

      A11TST = .FALSE.
      IF (A11TST) CALL DZERO(WORK(KRHO1),NT1AM(ISYRES))

      CALL CC_11A(WORK(KRHO1),WORK(KFCKLR),ISYRES,WORK(KLAMDH),
     *            WORK(KLAMDP),WORK(KEND1),LWRK1)

      END IF

      IF (DEBUG) THEN
         RHO1N = DDOT(NT1AM(ISYRES),WORK(KRHO1),1,WORK(KRHO1),1)
         WRITE(LUPRI,1) 'Norm of Rho1 -after CC_11A:        ',RHO1N
         WRITE (LUPRI,*) 'Rho1 after CC_11A contribution:'
         CALL CC_PRP(WORK(KRHO1),WORK,ISYRES,1,0)
      ENDIF
C
C--------------------------------------------------------------
C     transform CC_21I contribution to MO and add to RHO1:
C--------------------------------------------------------------
C
      ITST = .FALSE.
      IF (ITST)   CALL DZERO(WORK(KRHO1),NT1AM(ISYRES))
      CALL CC_T1AM(WORK(KRHO1),ISYRES,WORK(KRHO1AO),ISYRES,
     *             WORK(KLAMDH),ISYMOP,ONE)

      IF (DEBUG.OR.ITST) THEN
         WRITE (LUPRI,*) 'CC_21I contribution in AO:'
         WRITE (LUPRI,'(5F12.8)') (WORK(KRHO1AO+I),I=0,NT1AO(ISYRES)-1)
         WRITE (LUPRI,*) 'Rho1 after CC_21I contribution:'
         CALL CC_PRP(WORK(KRHO1),WORK,ISYRES,1,0)
         WRITE (LUPRI,*) 'WORK(0) = ',WORK(NULL)
      ENDIF
C
C===============================================================
C     After Loop CCSD and CC2 Contributions with C2SQ in memory.
C===============================================================
C
      IF ( .NOT. (CC2.OR.CCS)) THEN
C
C--------------------------------------------------------------
C        Transform the Omega2 vector to the MO basis.
C        We thus have the B(L2) term.
C--------------------------------------------------------------
C
         ICON   = 1
         IOPTG  = 0
         ISYMHC = 1
         ISYMBF = ISYRES
C
         LGAMMA  = .FALSE.
         LGIM    = .FALSE.
         LO3BF   = .FALSE.
         LBFZETA = .TRUE.
         DTIME   = SECOND()

         CALL DZERO(WORK(KC2AM),NT2AM(ISYRES))

         CALL CC_BFIFSORT(WORK(KRHO2),ISYRES,LUBFI,FNBFI,IADRBFI,
     *                    WORK(KEND1),LWRK1)

         CALL CC_T2MO3(DUMMY,DUMMY,ISYMOP,WORK(KRHO2),
     *                 WORK(KC2AM),DUMMY,DUMMY,DUMMY,
     *                 WORK(KLAMDH),ISYM0,WORK(KLAMDH),ISYM0,
     *                 WORK(KEND1),LWRK1,ISYMBF,
     *                 ICON,LGAMMA,IOPTG,LO3BF,LBFZETA)
         DTIME   = SECOND() -DTIME
         TIMT2MO = TIMT2MO + DTIME
         NEWGAM = .TRUE.
C
         CALL DCOPY(NT2AM(ISYRES),WORK(KC2AM),1,WORK(KRHO2),1)
C
         IF (DEBUG) THEN
            RHO1N = DDOT(NT1AM(ISYRES),WORK(KRHO1),1,WORK(KRHO1),1)
            RHO2N = DDOT(NT2AM(ISYRES),WORK(KRHO2),1,WORK(KRHO2),1)
            WRITE(LUPRI,1) 'Norm of Rho1 -after loop in MO:    ',RHO1N
            WRITE(LUPRI,1) 'Norm of Rho2 -after loop in MO:    ',RHO2N
            CALL AROUND('RHO2 after B(L2) term:')
            CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,0,1)
         ENDIF
C
         IF (IPRINT .GT. 120) THEN
            CALL AROUND('After  T2MO-1:BF(C1,C2) ')
            CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,1,1)
         ENDIF
C
C--------------------------------------
C     Calculate the LT21G contribution.
C--------------------------------------
C
         IF (DEBUG) THEN
           xn  = DDOT(N3ORHF(ISYRES),WORK(KMINTX),1,WORK(KMINTX),1)
           WRITE(LUPRI,1) 'Norm of MX   -before CC_21G:        ',xn   
           xn  = DDOT(N3ORHF(ISYML),WORK(KMINT),1,WORK(KMINT),1)
           WRITE(LUPRI,1) 'Norm of M    -before CC_21G:        ',xn   
         ENDIF
C
         TIMEG = SECOND()
C
         G21TST = .FALSE.
         IF (G21TST) CALL DZERO(WORK(KRHO1),NT1AM(ISYRES))
         CALL CC_21G(WORK(KRHO1),WORK(KMINTX),ISYRES,WORK(KLAMDH),
     *               WORK(KEND1),LWRK1,ISYMOP,LUO3T,O3TFIL)
C      
         IF (DEBUG.OR.G21TST) THEN
            CALL AROUND('RHO1 after first 21G contribution:')
            CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,1,0)
         ENDIF
C
         G21TST = .FALSE.
         IF (G21TST) CALL DZERO(WORK(KRHO1),NT1AM(ISYRES))
         CALL CC_21G(WORK(KRHO1),WORK(KMINT),ISYML,WORK(KLAMDH),
     *               WORK(KEND1),LWRK1,ISYMR,LUO3,O3FIL)
         TIMEG = SECOND() - TIMEG
C
         IF (DEBUG.OR.G21TST) THEN
            RHO1N = DDOT(NT1AM(ISYRES),WORK(KRHO1),1,WORK(KRHO1),1)
            RHO2N = DDOT(NT2AM(ISYRES),WORK(KRHO2),1,WORK(KRHO2),1)
            WRITE(LUPRI,1) 'Norm of Rho1 -after CC_21G:        ',RHO1N
            WRITE(LUPRI,1) 'Norm of Rho2 -after CC_21G:        ',RHO2N
            CALL AROUND('RHO1 after second 21G contribution:')
            CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,1,0)
         ENDIF
C
C--------------------------------------------------------------
C        Transform the Omega2 vector to the MO basis.
C        We thus have the B(C2) and (a,i-bar|bj) contributions.
C--------------------------------------------------------------
C
         IF (DEBUG) THEN
            RHO2N = DDOT(NRHOR,WORK(KRHOR),1,WORK(KRHOR),1)
            WRITE(LUPRI,1) 'Norm of BF(C1,C2) intermediate:    ',RHO2N
            RHO2N = DDOT(NT2SQ(ISYML),WORK(KL2AM),1,WORK(KL2AM),1)
            WRITE(LUPRI,1) 'Norm of L2:                        ',RHO2N
         ENDIF
C
         CALL DZERO(WORK(KGAMMC),NGAMMA(ISYMR))
C
         ICON   = 3
         IOPTG  = 0
         LGAMMA = .TRUE.
         LGIM   = .FALSE.
         LO3BF  = .TRUE.
         LBFZETA= .FALSE.
         ISYMPC = 1
         ISYMBF = ISYMR
C
         XNORM = DDOT(NT2AOIJ(ISYMR),WORK(KRHOR),1,WORK(KRHOR),1)
         WRITE (LUPRI,*) 'Norm(BF) = ',XNORM
         XNORM = DDOT(NT2SQ(ISYML),WORK(KL2AM),1,WORK(KL2AM),1)
         WRITE (LUPRI,*) 'Norm(ZETA2) = ',XNORM
C
         E2TST = .FALSE.
         IF (E2TST) CALL DZERO(WORK(KRHO1),NT1AM(ISYRES))
C
         DTIME = SECOND()
         CALL CC_T2MO3(WORK(KRHO1),WORK(KL2AM),ISYML,
     *                 WORK(KRHOR),DUMMY,WORK(KGAMMC),DUMMY,DUMMY,
     *                 WORK(KLAMDP),ISYM0,WORK(KLAMDP),ISYMPC,
     *                 WORK(KEND1),LWRK1,ISYMBF,
     *                 ICON,LGAMMA,IOPTG,LO3BF,LBFZETA)
         DTIME = SECOND() -DTIME
         TIMT2MO = TIMT2MO + DTIME
C
         IF ( DEBUG ) THEN
            XNGAM = DDOT(NGAMMA(ISYMR),WORK(KGAMMC),1,WORK(KGAMMC),1)
            RHO1N = DDOT(NT1AM(ISYRES),WORK(KRHO1),1,WORK(KRHO1),1)
            RHO2N = DDOT(NT2AM(ISYRES),WORK(KRHO2),1,WORK(KRHO2),1)
            WRITE(LUPRI,1) 'Norm of GamC -after T2MO:          ',XNGAM
            WRITE(LUPRI,1) 'Norm of Rho1 -after L2*BF(C1,C2):  ',RHO1N
            WRITE(LUPRI,1) 'Norm of Rho2 -after L2*BF(C1,C2):  ',RHO2N
            CALL AROUND('RHO after L2*BF(C1,C2):')
            CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,1,1)
            WRITE (LUPRI,*) 'WORK(0) = ',WORK(NULL)
         ENDIF
C
C----------------------------------
C        Read in integrals (ia|jb).
C----------------------------------
C
         REWIND(LUIAJB)
         READ(LUIAJB) (WORK(KC2AM-1+J), J = 1,NT2AM(ISYMOP))
C
C-----------------------------------------------------------------
C        Calculate 22EM contributions:
C        calculate 2(ia|jb) - (ja|ib) combination of the integrals
C        and calculate coulomb and exchange contributions together
C-----------------------------------------------------------------
C
         DTIME = SECOND()

         KXIM  = KEND1
         KYIM  = KXIM  + NMATIJ(ISYRES)
         KEND1 = KYIM  + NMATAB(ISYRES)
         LWRK1 = LWORK - KEND1

         E22TST = .FALSE.
         IF (E22TST) CALL DZERO(WORK(KRHO2),NT2AM(ISYRES))
     
         CALL DZERO(WORK(KXIM),NMATIJ(ISYRES))
         CALL DCOPY(NMATAB(ISYRES),WORK(KYMATX),1,WORK(KYIM),1)
         CALL DSCAL(NMATAB(ISYRES),XHALF,WORK(KYIM),1)
         IOPTTCME = 1
         CALL CCSD_TCMEPK(WORK(KC2AM),ONE,ISYMOP,IOPTTCME)

         CALL CC_22EC(WORK(KRHO2),WORK(KC2AM),WORK(KXIM),
     *                WORK(KYIM),ISYRES,WORK(KEND1),LWRK1)

         KEND1 = KXIM
         LWRK1 = LWORK - KEND1
C
         TIM2EM = SECOND() - DTIME

         IF (DEBUG) THEN
            RHO1N = DDOT(NT1AM(ISYRES),WORK(KRHO1),1,WORK(KRHO1),1)
            RHO2N = DDOT(NT2AM(ISYRES),WORK(KRHO2),1,WORK(KRHO2),1)
            WRITE(LUPRI,1) 'Norm of Rho1 -after Int*X + Int*Y: ',RHO1N
            WRITE(LUPRI,1) 'Norm of Rho2 -after Int*X + Int*Y: ',RHO2N
         ENDIF
C
         IF (IPRINT .GT. 120 .OR. E22TST .OR. DEBUG) THEN
            CALL AROUND('After Int*X and Int*Y E-terms) ')
            CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,1,1)
            CALL AROUND('KYIM:')
            CALL CC_PREI(WORK(KYIM),WORK(KXIM),ISYRES,1)
            CALL AROUND('Integrals:')
            CALL CC_PRP(WORK,WORK(KC2AM),ISYRES,0,1)
            WRITE (LUPRI,*) 'WORK(0) = ',WORK(NULL)
         ENDIF
C
C--------------------------------
C        write out result vector.
C--------------------------------
C
         DTIME   = SECOND()
         CALL CC_WVEC(LUTR,TRFIL,NRHO2,NT2AM(ISYRES),1,WORK(KRHO2))
         DTIME   = SECOND() - DTIME
         TIMIO   = TIMIO    + DTIME
C
      ENDIF
C
      IF (.NOT.(CCS.OR.CC2)) THEN
C
C--------------------------------
C        Read Omega intermediate.
C--------------------------------
C
         TIMEBF = SECOND()
         LUOM = -1
         CALL GPOPEN(LUOM,'CC_BFIM','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LUOM)
         READ(LUOM) (WORK(KT2AMP+I-1),I = 1,2*NT2ORT(1))
         CALL GPCLOSE(LUOM,'KEEP')
C
C------------------------------------------
C        Calculate the LT21BF contribution.
C------------------------------------------
C
         ICON = 3
C
         IF (DEBUG) THEN
            RHO2N=DDOT(2*NT2ORT(ISYMOP),WORK(KT2AMP),1,WORK(KT2AMP),1)
            WRITE(LUPRI,1) 'Norm of BF intermediate from disk: ',RHO2N
            RHO2N = DDOT(NGLMDT(ISYMR),WORK(KLAMPC),1,WORK(KLAMPC),1)
            WRITE(LUPRI,1) 'Norm of LMPC1: ',RHO2N
         ENDIF
C
         E2TST = .FALSE.
         IF (E2TST) CALL DZERO(WORK(KRHO1),NT1AM(ISYRES))
C
         NEWGAM = .FALSE.
         CALL CC_T2MO(WORK(KRHO1),WORK(KL2AM),ISYML,
     *                WORK(KT2AMP),PHONEY,DUMMY,
     *                WORK(KLAMDP),WORK(KLAMPC),ISYMR,WORK(KEND1),
     *                LWRK1,ISYMOP,ICON)
         NEWGAM = .TRUE.
         TIMEBF = SECOND() - TIMEBF
C
         IF (DEBUG) THEN
            RHO1N = DDOT(NT1AM(ISYRES),WORK(KRHO1),1,WORK(KRHO1),1)
            RHO2N = DDOT(NT2AM(ISYRES),WORK(KRHO2),1,WORK(KRHO2),1)
            WRITE(LUPRI,1) 'Norm of Rho1 -after CC_T2MO 21BF:  ',RHO1N
            WRITE(LUPRI,1) 'Norm of Rho2 -after CC_T2MO 21BF:  ',RHO2N
            CALL AROUND('RHO after L2*BF(T1,T2):')
            CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,1,1)
            WRITE (LUPRI,*) 'WORK(0) = ',WORK(NULL)
         ENDIF
C
C-----------------------------------
C        Readin C2 amplitude in RHO.
C-----------------------------------
C
         DTIME = SECOND()
         IOPT = 2
         CALL CC_RDRSP(RLIST,IRLNR,ISYMR,IOPT,MODELX,
     *                 DUMMY,WORK(KRHO2))
         DTIME = SECOND() - DTIME
         TIMIO = TIMIO + DTIME
C
C--------------------------------
C        Square up C2 amplitudes.
C--------------------------------
C
         DTIME = SECOND()
         CALL CCLR_DIASCL(WORK(KRHO2),TWO,ISYMR)
         CALL CC_T2SQ(WORK(KRHO2),WORK(KC2AM),ISYMR)
         DTIME = SECOND() - DTIME
         TIMT2SQ = TIMT2SQ + DTIME
C
         IF (IPRINT.GT.50) THEN
            CALL AROUND('CC_FMAT: (C1,C2) vector readin')
            CALL CC_PRSQ(WORK(KC1AM),WORK(KC2AM),ISYMR,1,1)
         ENDIF
C
C------------------------------
C        Read in result vector.
C------------------------------
C
         DTIME = SECOND()
         CALL CC_RVEC(LUTR,TRFIL,NRHO2,NT2AM(ISYRES),1,WORK(KRHO2))
         DTIME = SECOND() - DTIME
         TIMIO = TIMIO + DTIME
C
      ENDIF
C
      IF (IPRINT .GT. 15) THEN
         RHO1N = DDOT(NT1AM(ISYRES),WORK(KRHO1),1,WORK(KRHO1),1)
         WRITE(LUPRI,1) 'Norm of Rho1 -after loop in mo:     ',RHO1N
         IF (.NOT. CCS) THEN
            RHO2N = DDOT(NT2AM(ISYRES),WORK(KRHO2),1,WORK(KRHO2),1)
            WRITE(LUPRI,1) 'Norm of Rho2 -after loop in mo:     ',RHO2N
         ENDIF
      ENDIF
C
C------------------------------------------
C     Transform AO Fock matrix to MO basis.
C------------------------------------------
C
C
C-------------------------------------
C     C1 transformed Fock Matrix.
C     first make
C     F-dot: F`ai = SUM-k-(Lai,kk-bar).
C-------------------------------------
C
      ISYFCK = ISYMR
      ISYMPA = 1
      ISYMHO = 1
C
      IF (IPRINT .GT.140) THEN
         CALL AROUND( 'Fock AO matrix calc. from C1 transf. dens.' )
         CALL CC_PRFCKAO(WORK(KFOCKC),ISYFCK)
      ENDIF
C
      DTIME  = SECOND()
      CALL CC_FCKMO(WORK(KFOCKC),WORK(KLAMDP),WORK(KLAMDH),
     *              WORK(KEND1),LWRK1,ISYFCK,ISYMPA,ISYMHO)
      DTIME  = SECOND() - DTIME
      TIMFCKMO = DTIME  +  TIMFCKMO
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND( 'Fock MO matrix calc. from C1 transf. dens.' )
         CALL CC_PRFCKMO(WORK(KFOCKC),ISYFCK)
      ENDIF
C
C-----------------------------------------------------
C     Make Fock tilde:
C     F~ai = SUM-k-(Lai,kk-bar) + Fa-bar,i + Fa,i-bar.
C-----------------------------------------------------
C
      KFCKAO = KEND1
      KEND2  = KFCKAO + MAX(N2BST(ISYMOP),N2BST(ISYMR))
      LEND2  = LWORK  - KEND2
C
      IF (IPRINT .GT.140) THEN
         CALL AROUND( 'Fock tilde -1.' )
         CALL CC_PRFCKMO(WORK(KFOCKC),ISYMR)
      ENDIF
C
      ISYFAO = 1
      ISYMPA = ISYMR
      ISYMHO = 1
C
      CALL DCOPY(N2BST(ISYMOP),WORK(KFOCKT),1,WORK(KFCKAO),1)
      DTIME  = SECOND()
      CALL CC_FCKMO(WORK(KFCKAO),WORK(KLAMPC),WORK(KLAMDH),
     *              WORK(KEND2),LWRK2,ISYFAO,ISYMPA,ISYMHO)
      DTIME  = SECOND() - DTIME
      TIMFCKMO = DTIME  +  TIMFCKMO
C
      CALL DAXPY(N2BST(ISYMR),ONE,WORK(KFCKAO),1,WORK(KFOCKC),1)
C
      IF (IPRINT .GT.140) THEN
         CALL AROUND( 'Fock tilde - 2.' )
         CALL CC_PRFCKMO(WORK(KFOCKC),ISYMR)
      ENDIF
C
      ISYFAO = 1
      ISYMPA = 1
      ISYMHO = ISYMR
C
      CALL DCOPY(N2BST(ISYMOP),WORK(KFOCKT),1,WORK(KFCKAO),1)
      DTIME  = SECOND()
      CALL CC_FCKMO(WORK(KFCKAO),WORK(KLAMDP),WORK(KLAMHC),
     *              WORK(KEND2),LWRK2,ISYFAO,ISYMPA,ISYMHO)
      DTIME  = SECOND() - DTIME
      TIMFCKMO = DTIME  +  TIMFCKMO
C
      CALL DAXPY(N2BST(ISYMR),ONE,WORK(KFCKAO),1,WORK(KFOCKC),1)
C
      IF (DEBUG) THEN
         XE1 = DDOT(N2BST(ISYMR),WORK(KFOCKC),1,WORK(KFOCKC),1)
         WRITE(LUPRI,1) 'Norm of FCKC1 after construction:  ',XE1
      ENDIF
      IF (IPRINT .GT. 50) THEN
         CALL AROUND( 'Fock-Tilde MO matrix ' )
         CALL CC_PRFCKMO(WORK(KFOCKC),ISYMR)
      ENDIF
C
C-------------------------------------------------------------
C     Calculate simple fock contributions:
C     rho2(ai,bj) = Paibj(2L1am(ai)FCKMO(jb)-L1am(aj)*FCKMO(ib)
C-------------------------------------------------------------
C
      IF (.NOT. CCS) THEN
         IF (IPRINT .GT. 15) THEN
            RHO1N = DDOT(NT1AM(ISYRES),WORK(KRHO1),1,WORK(KRHO1),1)
            WRITE(LUPRI,1) 'Norm of Rho1 before L1FCK:         ',RHO1N
            RHO2N = DDOT(NT2AM(ISYRES),WORK(KRHO2),1,WORK(KRHO2),1)
            WRITE(LUPRI,1) 'Norm of Rho2 before L1FCK:         ',RHO2N
         ENDIF
C 
         HTST = .FALSE.
         IF (HTST) CALL DZERO(WORK(KRHO2),NT2AM(ISYRES))
C
         CALL CC_L1FCK(WORK(KRHO2),WORK(KL1AM),WORK(KFOCKC),ISYML,ISYMR,
     *                 WORK(KEND2),LWRK2)
C
         IF (HTST .OR. DEBUG) THEN
            CALL AROUND("Fock type contribution to H term:")
            CALL CC_PRP(WORK,WORK(KRHO2),ISYRES,0,1)
         END IF
C
         IF (DEBUG) THEN
            RHO1N = DDOT(NT1AM(ISYRES),WORK(KRHO1),1,WORK(KRHO1),1)
            WRITE(LUPRI,1) 'Norm of Rho1 -after L1FCK:         ',RHO1N
            RHO2N = DDOT(NT2AM(ISYRES),WORK(KRHO2),1,WORK(KRHO2),1)
            WRITE(LUPRI,1) 'Norm of Rho2 -after L1FCK:         ',RHO2N
            WRITE(LUPRI,*) '        '
            WRITE (LUPRI,*) 'WORK(0) = ',WORK(NULL)
         ENDIF
      ENDIF
C
C----------------------------------------------------
C     Calculate the 12C contribution:
C     rho2(ai,bj) = P(ai.bj)(-sum(k)CTR1(a,k)*L(jbik)
C     If .not.cc2 LT21B contributions to rho1.
C----------------------------------------------------
C
      CALL AROUND('before LT21B contribution:')
      CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,1,0)
C
C        SINGLETST = .TRUE.
C        IF (SINGLETST) THEN
C           CALL DZERO(WORK(KXMAT),NMATIJ(ISYML))
C           CALL DZERO(WORK(KXMATX),NMATIJ(ISYRES))
C           WRITE (LUPRI,*) 'Warning: zeroed 1.th order X intermediate!!!'
C        END IF
C
C
      TIME12C = SECOND()
      IF (.NOT. CCS) THEN

        C12TST = .FALSE.
        IF (C12TST) CALL DZERO(WORK(KRHO2),NT2AM(ISYRES))
        B21TST = .FALSE.
        IF (B21TST) CALL DZERO(WORK(KRHO1),NT1AM(ISYRES))

        IOPT = 2
        CALL CC_21B12C(WORK(KRHO1),WORK(KRHO2),WORK(KL1AM),ISYML,
     *                 WORK(KLAMDH),ISYMOP,WORK(KXMAT),ISYML,ISYMR,
     *                 WORK(KEND2),LWRK2,LUO3,O3FIL,IOPT)
        IF (C12TST .OR. DEBUG) THEN
           CALL AROUND('after 12C contribution:')
           CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,0,1)
           WRITE (LUPRI,*) 'WORK(0) = ',WORK(NULL)
        END IF

      ENDIF
      IF (.NOT.(CC2.OR.CCS)) THEN

        IOPT = 1
        CALL CC_21B12C(WORK(KRHO1),WORK(KRHO2),WORK(KL1AM),ISYML,
     *                 WORK(KLAMDH),ISYMOP,WORK(KXMATX),ISYRES,ISYMOP,
     *                 WORK(KEND2),LWRK2,LUO3T,O3TFIL,IOPT)
        
      ENDIF
      TIME12C = SECOND() - TIME12C
C
      IF (DEBUG .AND. (.NOT. CCS)) THEN
         RHO1N = DDOT(NT1AM(ISYRES),WORK(KRHO1),1,WORK(KRHO1),1)
         RHO2N = DDOT(NT2AM(ISYRES),WORK(KRHO2),1,WORK(KRHO2),1)
         WRITE(LUPRI,1) 'Norm of Rho1 -after CC_12C:        ',RHO1N
         WRITE(LUPRI,1) 'Norm of Rho2 -after CC_12C:        ',RHO2N
         CALL AROUND('after 21B contribution:')
         CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,1,0)
      ENDIF
C
C----------------------------------------
C     Calculate the LT21EFM contribution.
C----------------------------------------
C
      IF ((NFIELD.GT.0).AND.(CC2.AND.NONHF)) THEN
         ISYMXY = MULD2H(ISYML,ISYMR)
         IOPT = 1
         CALL CC_LCC2FF(WORK(KRHO1),WORK(KLAMDP),WORK(KLAMDH), 
     *                  WORK(KLAMPC),WORK(KLAMHC),ISYMOP,
     *                  WORK(KXMATX),WORK(KYMATX),ISYMXY,
     *                  WORK(KEND1),LWRK1,IOPT)
C
         CALL AROUND('after CC2 finite field contribution:')
         CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,1,0)
C
      ENDIF
      IF (.NOT.(CC2.OR.CCS)) THEN
C
         SINGLETST = .TRUE.
         IF (SINGLETST) THEN
C           CALL DZERO(WORK(KYMAT),NMATAB(ISYML))
C           CALL DZERO(WORK(KYMATX),NMATAB(ISYRES))
C           WRITE (LUPRI,*) 'Warning: zeroed 1.th order Y intermediate!!!'
C           CALL DZERO(WORK(KXMAT),NMATIJ(ISYML))
C           CALL DZERO(WORK(KXMATX),NMATIJ(ISYRES))
C           WRITE (LUPRI,*) 'Warning: zeroed 1.th order X intermediate!!!'
         END IF
C
         DTIME  = SECOND()
         ISYFCK = ISYMR
         ISYMXY = ISYML
         CALL CC_21EFM(WORK(KRHO1),WORK(KFOCKC),ISYFCK,WORK(KXMAT),
     *                 WORK(KYMAT),ISYMXY,
     *                 WORK(KEND1),LWRK1)
         DTIME  = SECOND() - DTIME 
         TIMEEM = DTIME    
         DTIME  = SECOND()
C
         ISYFCK = ISYMOP
         ISYMXY = MULD2H(ISYML,ISYMR)
         CALL CC_21EFM(WORK(KRHO1),WORK(KFOCK),ISYFCK,WORK(KXMATX),
     *                 WORK(KYMATX),ISYMXY,
     *                 WORK(KEND1),LWRK1)
         DTIME  = SECOND() - DTIME 
         TIMEEM = TIMEEM   + DTIME    
C
         CALL AROUND('after 21EFM contribution:')
         CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,1,0)
C
      ENDIF
C
C-----------------------------
C     write out result vector.
C-----------------------------
C
      IF (.NOT. CCS) THEN
         DTIME   = SECOND()
         CALL CC_WVEC(LUTR,TRFIL,NRHO2,NT2AM(ISYRES),1,WORK(KRHO2))
         DTIME   = SECOND() - DTIME
         TIMIO   = TIMIO    + DTIME
      ENDIF
C
C
C==========================================================
C     MO basis T2SQ section.
C     Contract intermediates constructed in loop with T2SQ.
C==========================================================
C
C-------------------------------------
C        Read in T2 amplitude in RHO2.
C-------------------------------------
C
         DTIME = SECOND()
C
         IF (.NOT. CCS) THEN
            IOPT = 2
            CALL CC_RDRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KRHO2))
C
            DTIME = SECOND() - DTIME
            TIMIO = TIMIO + DTIME
C
C-----------------------------------
C           Square up T2 amplitudes.
C-----------------------------------
C
            DTIME = SECOND()
            CALL CC_T2SQ(WORK(KRHO2),WORK(KC2AM),1)
            DTIME = SECOND() - DTIME
            TIMT2SQ = TIMT2SQ + DTIME
C
            IF (IPRINT.GT.50) THEN
               CALL AROUND('CC_FMAT: (T1,T2) vector readin')
               CALL CC_PRSQ(WORK(KL1AM),WORK(KC2AM),1,0,1)
            ENDIF
C
C---------------------------------
C           Read in result vector.
C---------------------------------
C
            DTIME = SECOND()
            CALL CC_RVEC(LUTR,TRFIL,NRHO2,NT2AM(ISYRES),1,WORK(KRHO2))
            DTIME = SECOND() - DTIME
            TIMIO = TIMIO + DTIME
         ENDIF
C
C----------------------------------------------------
C     Calculate E-term: T2*E(C1,C2)
C----------------------------------------------------
C
         IF ( IPRINT .GT. 50) THEN
            CALL AROUND( 'E-intermdiates calc. EI(C1,C2) - ao.')
            CALL CC_PREI(WORK(KEMAT1),WORK(KEMAT2),ISYMR,0)
         ENDIF
         IF ( DEBUG ) THEN
            XE1 = DDOT(NEMAT1(ISYMR),WORK(KEMAT1),1,WORK(KEMAT1),1)
            XE2 = DDOT(NMATIJ(ISYMR),WORK(KEMAT2),1,WORK(KEMAT2),1)
            WRITE(LUPRI,1) 'Norm of EI1  -ao:                  ',XE1
            WRITE(LUPRI,1) 'Norm of EI2  -ao:                  ',XE2  
            XE1 = DDOT(N2BST(ISYMR),WORK(KFOCKC),1,WORK(KFOCKC),1)
            WRITE(LUPRI,1) 'Norm of FCKCC1:                    ',XE1
         ENDIF
C
         FCKCON = .TRUE.
         ETRAN  = .TRUE.
         ISYMEI = ISYMR
         DTIME  = SECOND()
         CALL CCRHS_EFCK(WORK(KEMAT1),WORK(KEMAT2),WORK(KLAMDH),
     *                   WORK(KFOCKC),WORK(KEND1),LWRK1,FCKCON,
     *                   ETRAN,ISYMR)
         DTIME = SECOND() - DTIME
         TIMEFCK = TIMEFCK + DTIME
C
         IF ( DEBUG ) THEN
            XE1 = DDOT(NMATAB(ISYMR),WORK(KEMAT1),1,WORK(KEMAT1),1)
            XE2 = DDOT(NMATIJ(ISYMR),WORK(KEMAT2),1,WORK(KEMAT2),1)
            XL1 = DDOT(NT1AM(ISYML),WORK(KL1AM),1,WORK(KL1AM),1)
            WRITE(LUPRI,1) 'Norm of EI1  -mo:                  ',XE1
            WRITE(LUPRI,1) 'Norm of EI2  -mo:                  ',XE2  
         ENDIF
C
         IF (IPRINT.GT.54 .OR. DEBUG) THEN
            CALL AROUND( 'E-intermdiates calc EI(C1,C2) -mo' )
            CALL CC_PREI(WORK(KEMAT1),WORK(KEMAT2),ISYMR,1)
         ENDIF
C
         CALL CCLR_E1C1(WORK(KRHO1),WORK(KL1AM),
     *                  WORK(KEMAT1),WORK(KEMAT2),
     *                  WORK(KEND1),LWRK1,ISYML,ISYMR,'T')
C
         CALL AROUND("after CCLR_E1C1 (B2) contribution:")
         CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,1,0)
C
      IF (.NOT. (CCS.OR.CC2))  THEN
C
         ISYVEC = ISYML
         ISYMIM = ISYMR
         DTIME = SECOND()
C
C------------------------------------------------------------
C        Prepare the E-intermediates for contraction with L2. 
C------------------------------------------------------------
C
         CALL CC_EITR(WORK(KEMAT1),WORK(KEMAT2),WORK(KEND1),LWRK1,
     *                ISYMIM)
         IF ( DEBUG ) THEN
            XE1 = DDOT(NMATAB(ISYMR),WORK(KEMAT1),1,WORK(KEMAT1),1)
            XE2 = DDOT(NMATIJ(ISYMR),WORK(KEMAT2),1,WORK(KEMAT2),1)
            WRITE(LUPRI,1) 'Norm of EI1  -after EITR:          ',XE1
            WRITE(LUPRI,1) 'Norm of EI2  -after EITR:          ',XE2  
         ENDIF
C
C
C---------------------------------------
C        L2*EI(T1,T2,R1,R2) contraction.
C---------------------------------------
C
         CALL CCRHS_E(WORK(KRHO2),WORK(KL2AM),
     *                WORK(KEMAT1),WORK(KEMAT2),
     *                WORK(KEND1),LWRK1,ISYVEC,ISYMIM)
C
         DTIME = SECOND() - DTIME
         TIME  = TIME + DTIME
C
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' AFTER CCRHS_E 1.cont.  CCLR:RHO ')
            CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,0,1)
         ENDIF
C
      ENDIF
C
      IF ((NFIELD .GT.0).AND.(CC2.AND.NONHF)) THEN
         IOPT = -1
         CALL CC_FCC2FF(WORK(KRHO2),WORK(KL2AM),ISYML,
     *                  WORK(KLAMDP),WORK(KLAMDH), 
     *                  WORK(KLAMPC),WORK(KLAMHC),ISYMR,
     *                  WORK(KEND1),LWRK1,IOPT)
      ENDIF
C
      IF (DEBUG) THEN
         RHO1N = DDOT(NT1AM(ISYRES),WORK(KRHO1),1,WORK(KRHO1),1)
         WRITE(LUPRI,1) 'Norm of Rho1 -after EI(C) terms:   ',RHO1N
         IF (.NOT. CCS) THEN
            RHO2N = DDOT(NT2AM(ISYRES),WORK(KRHO2),1,WORK(KRHO2),1)
            WRITE(LUPRI,1) 'Norm of Rho2 -after EI(C) terms:   ',RHO2N
         ENDIF
         WRITE(LUPRI,*) '  '
      ENDIF
C
C----------------------------------------
C     Calculate A-term from gamma matrix.
C----------------------------------------
C
      IF  (.NOT. (CCS .OR. CC2)) THEN
C
         IF ( DEBUG ) THEN
            XGAM  = DDOT(NGAMMA(ISYMR),WORK(KGAMMC),1,WORK(KGAMMC),1)
            WRITE(LUPRI,1) 'Norm of Gamma matrix from T2MO:    ',XGAM 
         ENDIF
C
         ATST = .FALSE.
         IF (ATST) CALL DZERO(WORK(KRHO2),NT2AM(ISYRES))
C
         DTIME = SECOND()
         IOPT  = 2
         CALL CCRHS_A(WORK(KRHO2),WORK(KL2AM),WORK(KGAMMC),
     *                WORK(KEND1),LWRK1,
     *                ISYMR,ISYML,IOPT)
         DTIME  = SECOND() - DTIME
         TIMA   = TIMA     + DTIME
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYRES),WORK(KRHO1),1,WORK(KRHO1),1)
            RHO2N = DDOT(NT2AM(ISYRES),WORK(KRHO2),1,WORK(KRHO2),1)
            WRITE(LUPRI,1) 'Norm of Rho1 -after L*Gam(R,T):    ',RHO1N
            WRITE(LUPRI,1) 'Norm of Rho2 -after L*Gam(R,T):    ',RHO2N
            CALL AROUND('RHO2 after double excit. A term:')
            CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,0,1)
            WRITE (LUPRI,*) 'WORK(0) = ',WORK(NULL)
         ENDIF
C
      ENDIF
C
C------------------------------------
C     Recover work from gamma matrix.
C------------------------------------
C
      KEND1 = KENDS2
      LWRK1 = LWRKS2
C
C-------------------------------------------
C     Calculate the D-term: L2*DI(C1,C2,T1).
C-------------------------------------------
C
      IF  (.NOT. (CCS .OR. CC2)) THEN
C
         ISYDIM = ISYMR
         ISYVEC = ISYML
         IOPT = 2
C
         IF ( DEBUG ) THEN
            RHO2N = DDOT(NT2SQ(ISYML),WORK(KL2AM),1,WORK(KL2AM),1)
            WRITE(LUPRI,1) 'Norm of L2AM -before D-term(R,T):   ',RHO2N
         ENDIF
C
         D2TST = .FALSE.
         IF (D2TST) CALL DZERO(WORK(KRHO2),NT2AM(ISYRES))
C
         DTIME = SECOND()
         CALL CC_22D(WORK(KRHO2),WORK(KL2AM),ISYVEC,WORK(KLAMDH),
     *               WORK(KEND1),LWRK1,
     *               ISYDIM,LUD,DTFIL,1,IOPT)
         DTIME  = SECOND() - DTIME
         TIMDIO = TIM22D   + DTIME
C
         IF (IPRINT .GT. 120 .OR. D2TST .OR. DEBUG) THEN
            CALL AROUND(' AFTER (2T2-T2)*DI(C1) D-term RHO is ')
            CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,0,1)
            WRITE (LUPRI,*) 'WORK(0) = ',WORK(NULL)
         ENDIF
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYRES),WORK(KRHO1),1,WORK(KRHO1),1)
            RHO2N = DDOT(NT2AM(ISYRES),WORK(KRHO2),1,WORK(KRHO2),1)
            WRITE(LUPRI,1) 'Norm of Rho1 -after T2*DI(R,T):    ',RHO1N
            WRITE(LUPRI,1) 'Norm of Rho2 -after T2*DI(R,T):    ',RHO2N
         ENDIF
C
      ENDIF
C
C-------------------------------------------
C     Calculate the C-term: L2*CI(C1,C2,T1).
C-------------------------------------------
C
      IF  (.NOT. (CCS .OR. CC2)) THEN
C
         ISYCIM = ISYMR
         ISYVEC = ISYML
C
         IOPT = 2
C
C-------------------
C        Prepare L2.
C-------------------
C
         CALL CCRHS_T2BT(WORK(KL2AM),WORK(KEND1),LWRK1,ISYML)
         CALL DSCAL(NT2SQ(ISYML),THREE,WORK(KL2AM),1)
         CALL CCSD_T2TP(WORK(KL2AM),WORK(KEND1),LWRK1,ISYML)
C
C------------------
C        Calculate.
C------------------
C
         C2TST = .FALSE.
         IF (C2TST) CALL DZERO(WORK(KRHO2),NT2AM(ISYRES))
C
         DTIME   = SECOND()
         CALL CC_22C(WORK(KRHO2),WORK(KL2AM),ISYVEC,WORK(KLAMDH),
     *               WORK(KEND1),LWRK1,ISYCIM,
     *               LUC,CTFIL,1,IOPT)
         DTIME   = SECOND() - DTIME
         TIMCIO  = TIMCIO   + DTIME
C
         IF (IPRINT .GT. 120 .OR. C2TST .OR. DEBUG) THEN
            CALL AROUND(' AFTER T2*CIM(C1,C2) C-term RHO is ')
            CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,0,1)
            WRITE (LUPRI,*) 'WORK(0) = ',WORK(NULL)
         ENDIF
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYRES),WORK(KRHO1),1,WORK(KRHO1),1)
            RHO2N = DDOT(NT2AM(ISYRES),WORK(KRHO2),1,WORK(KRHO2),1)
            WRITE(LUPRI,1) 'Norm of Rho1 -after L2*CI(R,T):    ',RHO1N
            WRITE(LUPRI,1) 'Norm of Rho2 -after L2*CI(R,T):    ',RHO2N
            WRITE(LUPRI,*) '  '
         ENDIF
C
      ENDIF
C
C--------------------------------------
C     End loop over vectors in scratch.
C--------------------------------------
C
      KEND1 = KENDS2
      LWRK1 = LWRKS2
C
C=============================================================
C
C     End section.
C
C=============================================================
C
      IF (IPRINT .GT. 50 .OR. DEBUG) THEN
         CALL AROUND('END OF CC_FMAT :RHO ')
         NC2 = 1
         IF ( CCS ) NC2 = 0
         CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,1,NC2)
      ENDIF
C
      IF ((IPRINT .GT.15).OR.DEBUG) THEN
         RHO1N = DDOT(NT1AM(ISYRES),WORK(KRHO1),1,WORK(KRHO1),1)
         RHO2N = DDOT(NT2AM(ISYRES),WORK(KRHO2),1,WORK(KRHO2),1)
         WRITE(LUPRI,1) 'Norm of Rho1 -end of CC_FMAT:    ',RHO1N
         WRITE(LUPRI,1) 'Norm of Rho2 -end of CC_FMAT:    ',RHO2N
      ENDIF
C
C-------------
C     Timings.
C-------------
C
      TIMALL  = SECOND() - TIMALL
C
      IF (IPRINT .GT. 4) THEN
         TIMCCSD = TIMA    + TIMBF    + TIMC    + TIMD    + TIME  +
     *             TIMF    + TIMG     + TIMH    + TIMI    + TIMEI +
     *             TIMEFCK + TIME1C1  + TIM1C1F + TIMCIO  + TIMEZ +
     *             TIMDIO  + TIMLAM   + TIMAOD  + TIMTRBT + TIMEXI +
     *             TIMT2AO + TIMTCME  + TIMT2TR + TIMT2BT + TIMEYI +
     *             TIMT2MO + TIMFCKMO + TIMFCK  + TIMT2SQ + TIMEMI +
     *             TIMEH   + TIMC1T2  + TIM12B  + TIMEG   + TIM2EM +
     *             TIMEEM
         WRITE(LUPRI,'(/)')
         WRITE(LUPRI,9998) 'CC_FMAT   in total ', TIMALL
         WRITE(LUPRI,9998) 'CC part            ', TIMCCSD
         WRITE(LUPRI,9998) 'Int. calc. & read  ', TIMER1 + TIMER2
     *                  + TIMRDAO
         WRITE(LUPRI,9998) 'IO of vectors/I.M. ', TIMIO
         WRITE(LUPRI,*) ' '
      ENDIF
C
      IF (IPRINT .GT. 9) THEN
         WRITE(LUPRI,'(A)')
         IF ( .NOT. (CC2.OR.CCS)) THEN
            WRITE(LUPRI,9999) 'A               ', TIMA
            WRITE(LUPRI,9999) 'BF              ', TIMBF
            WRITE(LUPRI,9999) 'C               ', TIMC
            WRITE(LUPRI,9999) 'CIO             ', TIMCIO
            WRITE(LUPRI,9999) 'D               ', TIMD
            WRITE(LUPRI,9999) 'DIO             ', TIMDIO
            WRITE(LUPRI,9999) 'ZWV             ', TIMEZ
         ENDIF
         IF ( .NOT. CCS ) THEN
            WRITE(LUPRI,9999) 'E               ', TIME
            IF (CC2) WRITE(LUPRI,9999) 'F               ', TIMF
            WRITE(LUPRI,9999) 'G               ', TIMG
            WRITE(LUPRI,9999) '21G             ', TIMEG
            WRITE(LUPRI,9999) '22EM            ', TIM2EM
            WRITE(LUPRI,9999) 'LT21EFM         ', TIMEEM
            WRITE(LUPRI,9999) '12C             ', TIME12C
            WRITE(LUPRI,9999) 'H+21H           ', TIMH+TIMEH
            WRITE(LUPRI,9999) 'I               ', TIMI
            WRITE(LUPRI,9999) 'EI+21I          ', TIMEI
            WRITE(LUPRI,9999) 'X,Y,M           ', TIMEXI+TIMEYI+TIMEMI
            WRITE(LUPRI,9999) 'INT3O           ', TIME3O
            WRITE(LUPRI,9999) 'CC_C1T2C        ', TIMC1T2
         ENDIF
         WRITE(LUPRI,9999) 'EFCK            ', TIMEFCK
         WRITE(LUPRI,9999) 'E1C1            ', TIME1C1
         WRITE(LUPRI,9999) '1C1F            ', TIM1C1F
         WRITE(LUPRI,9999) 'LAMTRA          ', TIMLAM
         WRITE(LUPRI,9999) 'AODENS          ', TIMAOD
         WRITE(LUPRI,9999) 'FCK             ', TIMFCK
         WRITE(LUPRI,9999) 'TRBT            ', TIMTRBT
         WRITE(LUPRI,9999) '21DC            ', TIMEC
         WRITE(LUPRI,9999) '12B             ', TIM12B
         IF ( .NOT. CCS) THEN
            WRITE(LUPRI,9999) 'T2AO            ', TIMT2AO
            WRITE(LUPRI,9999) 'TCME            ', TIMTCME
            WRITE(LUPRI,9999) 'T2TR            ', TIMT2TR
            WRITE(LUPRI,9999) 'T2BT            ', TIMT2BT
            WRITE(LUPRI,9999) 'T2SQ            ', TIMT2SQ
            WRITE(LUPRI,9999) 'T2MO            ', TIMT2MO
         ENDIF
         WRITE(LUPRI,9999) 'FCKMO           ', TIMFCKMO
C
         WRITE(LUPRI,'(A)')
         WRITE(LUPRI,9999) 'RDAO            ', TIMRDAO
         WRITE(LUPRI,9999) 'ERIDI1          ', TIMER1
         WRITE(LUPRI,9999) 'ERIDI2          ', TIMER2
C
         IF (CCSDT) THEN
            WRITE(LUPRI,'(A)')
            WRITE(LUPRI,9999) 'T3INT           ', TIMT3I
            WRITE(LUPRI,9999) 'OMEG-1          ', TIMO31
            WRITE(LUPRI,9999) 'OMEG-2          ', TIMO32
            WRITE(LUPRI,9999) 'OMEG-3          ', TIMO33
         ENDIF
         WRITE(LUPRI,*) ' '
      ENDIF
      IF (IPRINT .GT. 15) THEN
         IF (DIRECT) WRITE(LUPRI,*) ' Atomic direct calculation'
         CALL AROUND(' END OF CC_FMAT ')
      ENDIF
C
   1  FORMAT(1x,A35,1X,E20.10)
9998  FORMAT(1x,'Time used in',2x,A20,2x,': ',f10.2,' seconds')
9999  FORMAT(1x,'Time used in',2x,A18,2x,': ',f10.2,' seconds')
C
C
C-----------------
C     Close files.
C-----------------
C
      IF (.NOT.(CCS.OR.CC2)) THEN
         CALL WCLOSE2(LUC,CTFIL,'DELETE')
         CALL WCLOSE2(LUD,DTFIL,'DELETE')
C
         CALL WCLOSE2(LUCIM,CFIL,'KEEP')
         CALL WCLOSE2(LUDIM,DFIL,'KEEP')
C
      ENDIF
C
C------------------------------------------------
C     Close (and delete) files for intermediates:
C------------------------------------------------
C
      IF (.NOT.(CCS.OR.CC2)) THEN
         IF (LLIST(1:2).EQ.'L0') THEN
            CALL WCLOSE2(LUPQR0,FILPQIM,'KEEP')
         ELSE
            CALL WCLOSE2(LUPQR0,FILPQR0,'DELETE')
         END IF
         CALL WCLOSE2(LUPQRR,FILPQRR,'DELETE')
      END IF
C
      CALL WCLOSE2(LUO3, O3FIL, 'DELETE')
      CALL WCLOSE2(LUO3T,O3TFIL,'DELETE')
C
      IF ( .NOT. CCS) THEN
         CALL WCLOSE2(LUTR, TRFIL,'DELETE')
      ENDIF
C
      IF (CCSD) THEN
         CALL WCLOSE2(LUBFI,FNBFI,'DELETE')
         CALL WCLOSE2(LUBFD,FNBFD,'DELETE')
      END IF
C
C-------------------------------
C     Restore T2TCOR and OMEGOR.
C-------------------------------
C
      T2TCOR = T2TSAV
      OMEGOR = ORSAVE
C
      RETURN

*---------------------------------------------------------------------*
* handle I/O error:
*---------------------------------------------------------------------*
999   CONTINUE
      WRITE (LUPRI,*) 'I/O error in CC_FMAT.'
      CALL QUIT('I/O error in CC_FMAT.')

      CALL QEXIT('CC_FMATTST')
C
      END
*=====================================================================*
*                      END OF SUBROUTINE CC_FMAT                   
*=====================================================================*
