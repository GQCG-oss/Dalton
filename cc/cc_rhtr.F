C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
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
#include <single.h>
C
c*DECK CC_RHTR
      SUBROUTINE CC_RHTR(ECURR,FRHO1,LUFR1,FRHO2,LUFR2,
     *                   FC1AM,LUFC1,FC2AM,LUFC2,
     *                   RHO1,RHO2,C1AM,C2AM,WORK,
     *                   LWORK,NSIMTR,IVEC,ITR,LRHO1)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Written by Ove Christiansen May 1994.
C
C     AO-integrals are used directly from disk.
C
C     Intermediate files are used:
C
C     BF,GAM,E1,E2
C
C     The biorthonormal basis is used.
C
C     (L: 1/2*Eia,(2*EiaEjb + EjaEib)/[6*(1 + delta(ab)delta(ij))] )
C     (R: Eai, EaiEbj)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxash.h>
#include <maxorb.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, XMTWO= -2.0D0)
      PARAMETER (XMHALF= -0.5D0, XHALF= 0.5D0)
C
#include <ccorb.h>
#include <ccisao.h>
#include <ccsdsym.h>
#include <ccsections.h>
#include <ccsdinp.h>
#include <ccfield.h>
#include <cclr.h>
#include <blocks.h>
#include <ccsdio.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
#include <ccslvinf.h>
#include <ccnoddy.h>
C
      DIMENSION INDEXA(MXCORB)
      DIMENSION RHO1(LRHO1,NSIMTR),C1AM(LRHO1,NSIMTR)
      DIMENSION RHO2(*),C2AM(*),WORK(LWORK)
      CHARACTER*10 MODEL, FN3SRTR
      CHARACTER*8 FRHO1, FRHO2, FC1AM, FC2AM, FNDELDR
      CHARACTER*8 FNCKJDR, FN3VI2, FNTOC
      CHARACTER*6 CFIL, DFIL, CTFIL, DTFIL, FNDKBCR, FN3VI
      CHARACTER*1 LR
C
      LOGICAL FCKCON,NEWGS,CC1BSA,ETRAN
C
      IF (IPRINT .GT. 15) THEN
         CALL AROUND(' START OF CC_RHTR ')
         IF (DIRECT) WRITE(LUPRI,*) ' Atomic direct calculation'
      ENDIF
C
      CC1BSA  = CC1B
      IF ( CC1A ) THEN
         CC1B    = .TRUE.
      ENDIF
      CC3LR = .FALSE.
C
      IF ( IPRINT .GT. 10 ) THEN
C
         WRITE(LUPRI,*) ' In CC_RHTR Right transformation routine:'
         WRITE(LUPRI,*) ' CCSDT, CC1A, CC1B, CC3: ',CCSDT, CC1A,
     &        CC1B, CC3
C
      ENDIF
C
C---------------------------
C     For ccs no c2 vectors.
C---------------------------
C
      NC2 = 1
      IF ( CCS ) NC2 = 0
C
      NRHO2 = MAX(NT2AM(ISYMTR),2*NT2ORT(ISYMTR))
      IF (CC2) NRHO2 = NT2AM(ISYMTR)
      IF (CCS) NRHO2 = 0
C
C----------------
C     Open files.
C----------------
C
      IF (DUMPCD .AND. (.NOT.(CCS .OR. CC2))) THEN
         LUC = 0
         LUD = 0
C
         CTFIL = 'CCLR_C'
         DTFIL = 'CCLR_D'
C
         CALL WOPEN2(LUC,CTFIL,64,0)
         CALL WOPEN2(LUD,DTFIL,64,0)
C
         LUCIM = 0
         LUDIM = 0
C
         CFIL = 'PMAT_C'
         DFIL = 'PMAT_D'
C
         CALL WOPEN2(LUCIM,CFIL,64,0)
         CALL WOPEN2(LUDIM,DFIL,64,0)
C
      END IF
C
      IF (CCSDT) THEN
         LU3SRTR = 0
         LUDELDR = 0
         LUCKJDR = 0
         LUDKBCR = 0
         LUTOC   = 0
         LU3VI   = 0
         LU3VI2  = 0
         FN3SRTR = 'CC3_SORT_R'
         FNDELDR = 'CKDELD_R'
         FNCKJDR = 'CKJDEL_R'
         FNDKBCR = 'DKBC_R'
         FNTOC   = 'CCSDT_OC'
         FN3VI   = 'CC3_VI'
         FN3VI2  = 'CC3_VI12'
         CALL WOPEN2(LU3SRTR,FN3SRTR,64,0)
         CALL WOPEN2(LUDELDR,FNDELDR,64,0)
         CALL WOPEN2(LUCKJDR,FNCKJDR,64,0)
         CALL WOPEN2(LUDKBCR,FNDKBCR,64,0)
         CALL WOPEN2(LUTOC,FNTOC,64,0)
         CALL WOPEN2(LU3VI,FN3VI,64,0)
         CALL WOPEN2(LU3VI2,FN3VI2,64,0)
      ENDIF
C
C------------------------
C     Time Initialiation.
C------------------------
C
      TIMALL    = SECOND()
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
      TIMEFCK   = 0.0D00
      TIME1C1   = 0.0D00
      TIM1C1F   = 0.0D00
      TIMCIO    = 0.0D00
      TIMDIO    = 0.0D00
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
C
      TIMT3I    = 0.0D00
      TIMTRIP   = 0.0D00
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
      NE1TOT = MAX(NEMAT1(ISYMTR)*NSIMTR,NEMAT1(1))
      NE2TOT = MAX(NMATIJ(ISYMTR)*NSIMTR,NMATIJ(1))
      N2C2C2 = 0
C
C-----------------------------------------------------------
C     Transposed T2 not needed in right transformation.
C     For CCSD it may be advantageous to keep 2C2-C in core.
C     For CC2 this is not the case.
C-----------------------------------------------------------
C
      IF ( T2TCOR ) N2C2C2 = NT2SQ(ISYMTR)
      IF (CCS.OR.CC2 ) NGAMTO = 2
      IF (CCS.OR.CC2) THEN
         NE1TOT = 2
         NE2TOT = 2
      ENDIF
C
      IF ( IPRINT .GT. 25)
     *   WRITE(LUPRI,*) ' In CC_RHTR work in start is:',LWORK
      K2C2C2 = 1
      KT1AM  = K2C2C2 + N2C2C2
      KLAMDP = KT1AM  + NT1AM(ISYMOP)
      KLAMDH = KLAMDP + NLAMDT
      KDENSI = KLAMDH + NLAMDT
      KFOCK  = KDENSI + N2BAST
      KEMAT1 = KFOCK  + N2BST(ISYMOP)
      KEMAT2 = KEMAT1 + NE1TOT
      KEND1  = KEMAT2 + NE2TOT
      LWRK1  = LWORK  - KEND1
C
C----------------------------------------------------------------
C     Extra Work space allocation for C1 dependent intermediates.
C----------------------------------------------------------------
C
      KLAMPC = KEND1
      KLAMHC = KLAMPC + NGLMDT(ISYMTR)*NSIMTR
      KDENSC = KLAMHC + NGLMDT(ISYMTR)*NSIMTR
      KFOCKC = KDENSC + N2BST(ISYMTR)*NSIMTR
      KFOCKT = KFOCKC + N2BST(ISYMTR)*NSIMTR
      KEND1  = KFOCKT + MAX(N2BST(ISYMTR),N2BST(ISYMOP))
      LWRK1  = LWORK  - KEND1
C
      IF ( IPRINT .GT. 25)
     *   WRITE(LUPRI,*) ' In CC_RHTR 1. alloc. work. left:',LWRK1
C
      IF (LWRK1 .LE. 0) CALL QUIT('Too little workspace in cclr_rhtr')
C
C---------------------------------------
C     Initialize t1am and intermediates.
C---------------------------------------
C
      CALL DZERO(WORK(KT1AM),NT1AM(1))
      CALL DZERO(WORK(KEMAT1),NE1TOT)
      CALL DZERO(WORK(KEMAT2),NE2TOT)
      CALL DZERO(WORK(KDENSI),N2BST(ISYMOP))
      CALL DZERO(WORK(KFOCK),N2BST(ISYMOP))
C
C-------------------------------------------
C     Initialize C1 dependent intermediates.
C-------------------------------------------
C
      CALL DZERO(WORK(KFOCKC),N2BST(ISYMTR)*NSIMTR)
      CALL DZERO(WORK(KDENSC),N2BST(ISYMTR)*NSIMTR)
C
C------------------------------------------------
C     Read the CC reference amplitudes from disk.
C------------------------------------------------
C
      DTIME = SECOND()
      IF (.NOT.(CCS.OR.CCP2)) THEN
         IOPT = 1
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),DUMMY)
      ENDIF
      DTIME = SECOND() - DTIME
      TIMIO = TIMIO + DTIME
C
C----------------------------------
C     Calculate the lamda matrices.
C----------------------------------
C
      DTIME = SECOND()
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     *            WORK(KEND1),LWRK1)
      DTIME = SECOND() - DTIME
      TIMLAM  = DTIME + TIMLAM
C
      IF (IPRINT .GT.100) THEN
         CALL AROUND('Usual Lambda matrices ')
         CALL CC_PRLAM(WORK(KLAMDP),WORK(KLAMDH),1)
      ENDIF
C
C----------------------------------
C     Calculate the density matrix.
C----------------------------------
C
      DTIME  = SECOND()
      ISYMH  = 1
      IC     = 1
      CALL CC_AODENS(WORK(KLAMDP),WORK(KLAMDH),WORK(KDENSI),ISYMH,
     *               IC,WORK(KEND1),LWRK1)
      DTIME  = SECOND() - DTIME
      TIMAOD = TIMAOD + DTIME
C
      IF (IPRINT .GT. 45) THEN
         CALL AROUND('CC_RHTR: Usual Lamda density matrix')
         CALL CC_PRAODEN(WORK(KDENSI),1)
      ENDIF
C
      DO 50 IV = 1, NSIMTR
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYMTR),C1AM(1,IV),1,C1AM(1,IV),1)
            WRITE(LUPRI,1) 'Norm of C1AM -first in CC_RHTR:  ',RHO1N
            IF (.NOT. CCS) THEN
               RHO2N = DDOT(NT2SQ(ISYMTR),C2AM,1,C2AM,1)
               WRITE(LUPRI,1) 'Norm of C2AM -first in CC_RHTR:  ',RHO2N
            ENDIF
         ENDIF
C
         KLPC1  = KLAMPC + NGLMDT(ISYMTR)*(IV-1)
         KLHC1  = KLAMHC + NGLMDT(ISYMTR)*(IV-1)
         KDNSC1 = KDENSC + N2BST(ISYMTR)*(IV-1)
C
C-----------------------------------------
C        Transform lamda matrices with C1.
C-----------------------------------------
C
         DTIME = SECOND()
         CALL CCLR_LAMTRA(WORK(KLAMDP),WORK(KLPC1),
     *                    WORK(KLAMDH),WORK(KLHC1),C1AM(1,IV),ISYMTR)
         DTIME = SECOND() - DTIME
         TIMLAM  = DTIME + TIMLAM
C
         IF (IPRINT .GT. 45) THEN
            CALL AROUND('C1 transformed Lambda matrices ')
            CALL CC_PRLAM(WORK(KLPC1),WORK(KLHC1),ISYMTR)
         ENDIF
C
C---------------------------
C        Print C amplitudes.
C---------------------------
C
         IF (IPRINT .GT. 50) THEN
            CALL AROUND('START OF CC_RHTR: (C1) ')
            CALL CC_PRSQ(C1AM(1,IV),C2AM,ISYMTR,1,0)
         ENDIF
         IF (IPRINT .GT. 50) THEN
            CALL AROUND('START OF CC_RHTR: RHO1 ')
            CALL CC_PRSQ(RHO1(1,IV),C2AM,ISYMTR,1,0)
         ENDIF
C
C----------------------------------------------------
C        Calculate the C1 transformed density matrix.
C        no core contribution ic = 0
C----------------------------------------------------
C
         ISYMH  = ISYMTR
         DTIME  = SECOND()
         IC     = 0
         CALL CC_AODENS(WORK(KLAMDP),WORK(KLHC1),WORK(KDNSC1),ISYMH,
     *               IC,WORK(KEND1),LWRK1)
         DTIME  = SECOND() - DTIME
         TIMAOD = TIMAOD + DTIME
C
         IF (IPRINT .GT. 50) THEN
            CALL AROUND('CC_RHTR: C1 trans. Lamda density matrix')
            CALL CC_PRAODEN(WORK(KDNSC1),ISYMTR)
         ENDIF
C
  50  CONTINUE
C
C--------------------------------------------------
C     If t2tcor then save 2C2 - C2 in core.
C     note that is it NOT the transposed vector.
C--------------------------------------------------
C
      IF ( T2TCOR .AND. MINSCR) THEN
         DTIME = SECOND()
         CALL DCOPY(NT2SQ(ISYMTR),C2AM,1,WORK(K2C2C2),1)
         CALL CCRHS_T2TR(WORK(K2C2C2),WORK(KEND1),LWRK1,ISYMTR)
         DTIME = SECOND() - DTIME
         TIMT2TR = TIMT2TR + DTIME
      ENDIF
C
C------------------------------------------------
C     Read one-electron integrals in Fock-matrix.
C------------------------------------------------
C
      TIMFCK = SECOND()
      CALL CCRHS_ONEAO(WORK(KFOCK),WORK(KEND1),LWRK1)
      TIMFCK = SECOND() - TIMFCK
C
C------------------------------------------------
C     Read one-electron integrals into Fock-matrix for
C     finite field.
C------------------------------------------------
C
      DO 13 IF = 1, NFIELD
         DTIME  = SECOND()
         FF = EFIELD(IF)
         CALL CC_ONEP(WORK(KFOCK),WORK(KEND1),LWRK1,FF,1,LFIELD(IF))
         DTIME  = DTIME - SECOND()
         TIMFCK = TIMFCK + DTIME
  13  CONTINUE
C
C----------------------------------------------------
C     Solvent contribution to one-electron
C     integrals. T^g contribution to transformation.
C SLV98,OC
C----------------------------------------------------
C
      IF (CCSLV) THEN
         CALL CCSL_RHSTG(WORK(KFOCK),WORK(KEND1),LWRK1)
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
         DTIME  = SECOND()
C
         IF (HERDIR) THEN
            CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
         ELSE
            KCCFB1 = KEND1
            KINDXB = KCCFB1 + MXPRIM*MXCONT
            KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
            LWRK1  = LWORK  - KEND1
            CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     &                  KODPP1,KODPP2,KRDPP1,KRDPP2,
     &                  KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     &                  WORK(KEND1),LWRK1,IPRERI)
            KEND1 = KFREE
            LWRK1 = LFREE
         ENDIF
         DTIME  = SECOND() - DTIME
         TIMER1 = TIMER1 + DTIME
         NTOSYM = 1
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      KENDSV = KEND1
      LWRKSV = LWRK1
C
      ICDEL1 = 0
      ICDEL2 = 0
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
     &                        IPRERI)
               ELSE
                  CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     &                        WORK(KODCL1),WORK(KODCL2),
     &                        WORK(KODBC1),WORK(KODBC2),
     &                        WORK(KRDBC1),WORK(KRDBC2),
     &                        WORK(KODPP1),WORK(KODPP2),
     &                        WORK(KRDPP1),WORK(KRDPP2),
     &                        WORK(KCCFB1),WORK(KINDXB),
     &                        WORK(KEND1), LWRK1,IPRERI)
               ENDIF
               DTIME  = SECOND() - DTIME
               TIMER2 = TIMER2 + DTIME
C
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient core in CC_RHTR')
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
            DO 120 IV = 1, NSIMTR
C
               KLPC1  = KLAMPC + NGLMDT(ISYMTR)*(IV-1)
               KLHC1  = KLAMHC + NGLMDT(ISYMTR)*(IV-1)
               KDNSC1 = KDENSC + N2BST(ISYMTR)*(IV-1)
               KFCKC1 = KFOCKC + N2BST(ISYMTR)*(IV-1)
               KEI1   = KEMAT1 + NEMAT1(ISYMTR)*(IV-1)
               KEI2   = KEMAT2 + NMATIJ(ISYMTR)*(IV-1)
C
               IF ((.NOT.MINSCR).AND.(.NOT.CCS)) THEN
C
C--------------------------------------------
C                 Readin C2 amplitude in RHO.
C--------------------------------------------
C
                  DTIME = SECOND()
                  CALL CC_RVEC(LUFC2,FC2AM,NT2AM(ISYMTR),NT2AM(ISYMTR),
     *                         IVEC+IV-1,RHO2)
                  DTIME = SECOND() - DTIME
                  TIMIO = TIMIO + DTIME
C
C-----------------------------------------
C                 Square up C2 amplitudes.
C-----------------------------------------
C
                  DTIME = SECOND()
                  CALL CCLR_DIASCL(RHO2,TWO,ISYMTR)
                  CALL CC_T2SQ(RHO2,C2AM,ISYMTR)
                  DTIME = SECOND() - DTIME
                  TIMT2SQ = TIMT2SQ + DTIME
C
                  IF (IPRINT.GT.50) THEN
                     CALL AROUND('CC_RHTR: (C1,C2) vector readin')
                     CALL CC_PRSQ(C1AM(1,IV),C2AM,ISYMTR,1,1)
                  ENDIF
C
C---------------------------------------
C                 Read in result vector.
C---------------------------------------
C
                  DTIME = SECOND()
                  CALL CC_RVEC(LUFR2,FRHO2,NRHO2,NRHO2,
     *                         IV+ITR -1,RHO2)
                  DTIME = SECOND() - DTIME
                  TIMIO = TIMIO + DTIME
C
                  IF (IPRINT.GT.50) THEN
                     CALL AROUND('CC_RHTR: RHO2  ) vector readin')
                     CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,1,1)
                  ENDIF
C
               ENDIF
C
C--------------------------------------------------------
C              If t2tcor then save 2C2 - C2 in core.
C              note that is it NOT the transposed vector.
C--------------------------------------------------------
C
               IF ( T2TCOR .AND.(.NOT. MINSCR)) THEN
                  DTIME   = SECOND()
                  CALL DCOPY(NT2SQ(ISYMTR),C2AM,1,WORK(K2C2C2),1)
                  CALL CCRHS_T2TR(WORK(K2C2C2),WORK(KEND1),
     *                            LWRK1,ISYMTR)
                  DTIME   = SECOND() - DTIME
                  TIMT2TR = TIMT2TR  + DTIME
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
                  ISYMD = ISAO(IDEL)
               ELSE
                  IDEL  = IBAS(ISYMD1) + ILLL
                  ISYMD = ISYMD1
               ENDIF
C
               ISYDIS = MULD2H(ISYMD,ISYMOP)
               ISYAIK = MULD2H(ISYDIS,ISYMTR)
C
C----------------------------------------------------------
C              Calculate adresses for c,cio,d,dio routines.
C----------------------------------------------------------
C
               IF ( IV .EQ. 1) THEN
                  IT2DEL(IDEL) = ICDEL1
                  ICDEL1 = ICDEL1 + NT2BCD(ISYDIS)
               ENDIF
C
               IT2DLR(IDEL,IV) = ICDEL2
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
     *            WRITE(LUPRI,*) ' In CC_RHTR 2. alloc. work left:',
     *            LWRK2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CC_RHTR ')
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
C-----------------------------------------------------------
C              Calculate transformed integrals used in t3am.
C-----------------------------------------------------------
C
               IF (CCSDT) THEN
C
                  DTIME = SECOND()
                  CALL CC3_T3INT(WORK(KXINT),WORK(KLAMDP),WORK(KLAMDH),
     *                           DUMM,1,WORK(KEND2),LWRK2,IDEL,ISYMD,1,
     *                           LU3SRTR,FN3SRTR,LUCKJDR,FNCKJDR)
                  DTIME   = SECOND() - DTIME
                  TIMT3I  = TIMT3I   + DTIME
C
               ENDIF
C
C---------------------------------------------
C              Calculate the AO-Fock matrix.
C              Not neccesary if RSPIM is used.
C---------------------------------------------
C
               IF ( CCS ) THEN
                  ISYDEN  = 1
                  IF (IV .EQ. 1) THEN
                     CALL CC_AOFOCK(WORK(KXINT),WORK(KDENSI),
     *                              WORK(KFOCK),WORK(KEND2),
     *                              LWRK2,IDEL,ISYMD,ISYDEN)
                  ENDIF
               ENDIF
C
C------------------------------------------------------------------
C              Calculate an AO-Fock matrix with C1. trans. density.
C------------------------------------------------------------------
C
               ISYDEN  = ISYMTR
               DTIME   = SECOND()
               CALL CC_AOFOCK(WORK(KXINT),WORK(KDNSC1),WORK(KFCKC1),
     *                        WORK(KEND2),LWRK2,IDEL,ISYMD,ISYDEN)
               DTIME  = SECOND() - DTIME
               TIMFCK = TIMFCK + DTIME
C
C-----------------------------------------
C              IF CCS jump to end of loop.
C-----------------------------------------
C
               IF (CCS) GOTO 130
C
C------------------------------------------
C              Work space allocation no. 3.
C------------------------------------------
C
               KDSRHF = KEND2
               KEND3  = KDSRHF + NDSRHF(ISYMD)
               LWRK3  = LWORK  - KEND3
C
               IF (LWRK3 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND3,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CC_RHTR ')
               ENDIF
C
C--------------------------------------------------------
C              Transform one index in the integral batch.
C--------------------------------------------------------
C
               DTIME   = SECOND()
               ISYMLP  = 1
               CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KLAMDP),
     *                          ISYMLP,WORK(KEND3),LWRK3,ISYDIS)
               DTIME   = SECOND() - DTIME
               TIMTRBT = TIMTRBT + DTIME
C
               IF (IPRINT .GT. 200) THEN
                  CALL CC_PRINT(WORK(KXINT),WORK(KDSRHF),ISYDIS)
               ENDIF
C
C--------------------------------------------------------------------
C              Construct the partially transformed C2-amplitudes, CM.
C--------------------------------------------------------------------
C
               KSCRCM = KEND3
               KEND4  = KSCRCM + NT2MMO(ISYMD,ISYMTR)
               LWRK4  = LWORK  - KEND4
C
               DTIME   = SECOND()
               ICON = 1
               ISYMLH = 1
               CALL CC_T2AO(C2AM,WORK(KLAMDH),ISYMLH,WORK(KSCRCM),
     *                         WORK(KEND4),LWRK4,IDEL,ISYMD,
     *                         ISYMTR,ICON)
               DTIME   = SECOND() - DTIME
               TIMT2AO = TIMT2AO + DTIME
               IF (IPRINT .GT. 100) THEN
                  CALL CC_PRTM(WORK(KSCRCM),ISYMD,ISYMTR)
               ENDIF
C
C-----------------------------------
C              Calculate the B-term.
C-----------------------------------
C
               IF (.NOT.(CCS .OR. CC2)) THEN
                  IOPT = 3
                  DTIME   = SECOND()
                  ISYMM  = MULD2H(ISYMD,ISYMTR)
                  CALL CC_BF(WORK(KXINT),RHO2,WORK(KLAMDH),
     *                       1,WORK(KLHC1),ISYMTR,WORK(KLAMDH),1,
     *                       WORK(KSCRCM),ISYMM,DUMMY,1,
     *                       WORK(KEND4),LWRK4,IDEL,
     *                       ISYMD,IOPT)
                  DTIME   = SECOND() - DTIME
                  TIMBF   = TIMBF    + DTIME
               ENDIF
C
C------------------------------------------------------------------
C              Calculate (ai|bj) tilde: MO 2C1 contribution in CC2.
C------------------------------------------------------------------
C
               IF ( CC2 ) THEN
C
                  IOPT = 2
                  DTIME   = SECOND()
                  CALL CC_MOFCON(WORK(KXINT),RHO2,WORK(KLAMDP),
     *                           WORK(KLAMDH),WORK(KLPC1),
     *                           WORK(KLHC1),WORK(KEND4),LWRK4,
     *                           IDEL,ISYMD,ISYMTR,IOPT)
                  DTIME   = SECOND() - DTIME
                  TIMF    = TIMF     + DTIME
C
                  IF (IPRINT .GT. 120) THEN
                     CALL AROUND(' AFTER CC_MOFCON in loop RHTR: RHO')
                     CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,0,1)
                  ENDIF
C
               ENDIF
C
C---------------------------------------
C              Transform CM to 2CM - CM.
C---------------------------------------
C
               DTIME = SECOND()
               CALL CC_MTCME(WORK(KSCRCM),WORK(KEND4),LWRK4,
     *                         ISYMD,ISYMTR)
               DTIME   = SECOND() - DTIME
               TIMTCME = TIMTCME  + DTIME
C
C-------------------------------------------------------
C              Calculate the D-tilde local intermediate.
C              Note: No c2am contribution for D-tilde.
C-------------------------------------------------------
C
               IF (.NOT.(CCS .OR. CC2)) THEN
C
                  ISYMPC = ISYMTR
                  ISYMHC = ISYMTR
                  ICON   = 1
                  FACTD  = 0.0D00
C
                  DTIME   = SECOND()
                  CALL CCRHS_D(WORK(KXINT),WORK(KDSRHF),RHO2,
     *                         C2AM,ISYMTR,
     *                         WORK(KLAMDP),DUMMY,WORK(KLAMDH),
     *                         WORK(KLPC1),ISYMPC,WORK(KLHC1),ISYMHC,
     *                         WORK(KSCRCM),WORK(KEND4),LWRK4,IDEL,
     *                         ISYMD,FACTD,ICON,LUD,DTFIL,IV)
                  DTIME   = SECOND() - DTIME
                  TIMD    = TIMD     + DTIME
C
               ENDIF
C
C-------------------------------------------------------
C              Calculate the C-tilde local intermediate.
C              Note: No c2am contribution for C-tilde.
C-------------------------------------------------------
C
               IF (.NOT.(CCS .OR. CC2)) THEN
C
                  FACTC = 0.0D00
                  ICON = 1
C
                  DTIME   = SECOND()
                  CALL CCRHS_C(WORK(KXINT),WORK(KDSRHF),RHO2,
     *                         C2AM,ISYMTR,WORK(KLAMDP),DUMMY,
     *                         WORK(KLAMDH),WORK(KLPC1),ISYMTR,
     *                         WORK(KLHC1),ISYMTR,
     *                         WORK(KSCRCM),WORK(KEND4),
     *                         LWRK4,IDEL,ISYMD,FACTC,ICON,LUC,
     *                         CTFIL,IV)
                  DTIME   = SECOND() - DTIME
                  TIMC    = TIMC     + DTIME
               ENDIF
C
C---------------------------------------
C              Transform C2 to 2C2 - C2.
C---------------------------------------
C
               IF (( .NOT. CC2 ).AND.(.NOT. T2TCOR)) THEN
                  DTIME   = SECOND()
                  CALL CCRHS_T2TR(C2AM,WORK(KEND4),LWRK4,ISYMTR)
                  DTIME   = SECOND() - DTIME
                  TIMT2TR = TIMT2TR  + DTIME
               ENDIF
C
C----------------------------------------
C              Calculate E-intermediates.
C----------------------------------------
C
               DTIME   = SECOND()
               IF (.NOT.(CCS .OR. CC2)) THEN
                  IF (.NOT. T2TCOR) THEN
                     CALL CCRHS_EI(WORK(KDSRHF),WORK(KEI1),
     *                             WORK(KEI2),C2AM,
     *                             WORK(KSCRCM),WORK(KLAMDP),
     *                             WORK(KLAMDH),WORK(KEND4),LWRK4,
     *                             IDEL,ISYMD,ISYDIS,ISYMTR)
                  ELSE
                     CALL CCRHS_EI(WORK(KDSRHF),WORK(KEI1),
     *                             WORK(KEI2),WORK(K2C2C2),
     *                             WORK(KSCRCM),WORK(KLAMDP),
     *                             WORK(KLAMDH),WORK(KEND4),LWRK4,
     *                             IDEL,ISYMD,ISYDIS,ISYMTR)
                  ENDIF
               ENDIF
               DTIME   = SECOND() - DTIME
               TIMEI   = TIMEI    + DTIME
C
C-----------------------------------
C              Calculate the G-term.
C-----------------------------------
C
               IF (.NOT.ONLY21) THEN
                  ISYMP1 = 1
                  ISYMH1 = 1
                  DTIME   = SECOND()
                  CALL CCRHS_G(WORK(KDSRHF),RHO1(1,IV),WORK(KLAMDP),
     *                         ISYMP1,WORK(KLAMDH),ISYMH1,
     *                         WORK(KSCRCM),WORK(KEND4),LWRK4,
     *                         ISYDIS,ISYMD,ISYMTR)
                  DTIME   = SECOND() - DTIME
                  TIMG    = TIMG     + DTIME
                  IF (IPRINT .GT. 120) THEN
                     CALL AROUND('After G-term  IN LOOP CCLR:RHO ')
                     CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,1,0)
                  ENDIF
               ENDIF
C
C-----------------------------------
C              Calculate the H-term.
C-----------------------------------
C
               IF (.NOT.ONLY21) THEN
                  DTIME   = SECOND()
                  CALL CCRHS_H(WORK(KDSRHF),RHO1(1,IV),WORK(KLAMDP),
     *                         WORK(KLAMDH),WORK(KSCRCM),WORK(KEND4),
     *                         LWRK4,ISYDIS,ISYMD,ISYMTR)
                  DTIME   = SECOND() - DTIME
                  TIMH    = TIMH     + DTIME
                  IF (IPRINT .GT. 120) THEN
                     CALL AROUND('After H-term IN LOOP  CCLR:RHO ')
                     CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,1,0)
                  ENDIF
               ENDIF
C
C---------------------------------------------
C              BackTransform C2 from 2C2 - C2.
C---------------------------------------------
C
               IF (( .NOT. CC2 ).AND. (.NOT. T2TCOR)) THEN
                  DTIME   = SECOND()
                  CALL CCRHS_T2BT(C2AM,WORK(KEND4),LWRK4,ISYMTR)
                  DTIME   = SECOND() - DTIME
                  TIMT2BT = TIMT2BT  + DTIME
               ENDIF
C
               IF (IPRINT .GT. 120) THEN
                  CALL AROUND('END of LOOP  CCLR:RHO ')
                  CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,1,1)
               ENDIF
C
               IF (IPRINT .GT. 50) THEN
                  RHO1N = DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
                  IF (OMEGOR .AND. (.NOT. CCS)) THEN
                    RHO2N = DDOT(2*NT2ORT(ISYMTR),RHO2,1,RHO2,1)
                  ELSE IF (.NOT. CCS) THEN
                    RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
                  ENDIF
                  WRITE(LUPRI,1) 'Norm of Rho1 -in loop:             ',
     &                 RHO1N
                  WRITE(LUPRI,1) 'Norm of Rho2 -in loop:             ',
     &                 RHO2N
               ENDIF
C
  130       CONTINUE
C
C--------------------------------------
C              write out result vector.
C--------------------------------------
C
               DTIME   = SECOND()
               CALL CC_WVEC(LUFR2,FRHO2,NRHO2,NRHO2,IV + ITR -1,RHO2)
               DTIME   = SECOND() - DTIME
               TIMIO   = TIMIO    + DTIME
C
  120       CONTINUE
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
C--------------------------------------------
C     Allocate space for the gamma matrix.
C--------------------------------------------
C
      IF (NEWGAM) THEN
C
         KGAMMC = KEND1
         KEND1  = KGAMMC + NGAMMA(ISYMTR)*NSIMTR
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) CALL QUIT('Insufficient memory in cclr_rhtr')
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
         LFOCK = 0
         CALL GPOPEN(LFOCK,'CC_FCKH','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
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
C----------------------------------------------------
C     Loop over number of trial vectors nsimtr.
C----------------------------------------------------
C
      DO 150 IV = 1, NSIMTR
C
         KLPC1  = KLAMPC + NGLMDT(ISYMTR)*(IV-1)
         KLHC1  = KLAMHC + NGLMDT(ISYMTR)*(IV-1)
         KFCKC1 = KFOCKC + N2BST(ISYMTR)*(IV-1)
         KGAMC1 = KGAMMC + NGAMMA(ISYMTR)*(IV-1)
C
C
C===============================================================
C     After Loop CCSD and CC2 Contributions with C2SQ in memory.
C===============================================================
C
      IF ( .NOT. ( CCS .OR. CC2 ) ) THEN
C
C---------------------------------
C           Read in result vector.
C---------------------------------
C
         IF ( .NOT. MINSCR ) THEN
C
            DTIME = SECOND()
            CALL CC_RVEC(LUFR2,FRHO2,NRHO2,NRHO2,IV + ITR -1,RHO2)
            DTIME = SECOND() - DTIME
            TIMIO = TIMIO + DTIME
C
         ENDIF
C
C--------------------------------------------------------------
C        Transform the Omega2 vector to the MO basis.
C        We thus have the B(C2) and (a,i-bar|bj) contributions.
C--------------------------------------------------------------
C
         ICON = 1
         ISYMPC = 1
         ISYMBF = ISYMTR
C
         DTIME = SECOND()
         CALL CC_T2MO(FAKE,PHONEY,ISYMOP,RHO2,C2AM,WORK(KGAMC1),
     *                WORK(KLAMDP),WORK(KLAMDP),ISYMPC,
     *                WORK(KEND1),LWRK1,ISYMBF,ICON)
         DTIME = SECOND() -DTIME
         TIMT2MO = TIMT2MO + DTIME
C
         CALL DCOPY(NT2AM(ISYMTR),C2AM,1,RHO2,1)
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of Rho1 -after loop in MO:    ',RHO1N
            WRITE(LUPRI,1) 'Norm of Rho2 -after loop in MO:    ',RHO2N
         ENDIF
C
         IF (IPRINT .GT. 120) THEN
            CALL AROUND('After  T2MO-1:BF(C1,C2) ')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,1,1)
         ENDIF
C
         IF ( DEBUG ) THEN
            XNGAM = DDOT(NGAMMA(ISYMTR),WORK(KGAMC1),1,WORK(KGAMC1),1)
            WRITE(LUPRI,1) 'Norm of GamC -after T2MO:          ',XNGAM
         ENDIF
C
C--------------------------------
C        write out result vector.
C--------------------------------
C
         DTIME   = SECOND()
         CALL CC_WVEC(LUFR2,FRHO2,NRHO2,NT2AM(ISYMTR),
     *                IV + ITR -1,RHO2)
         DTIME   = SECOND() - DTIME
         TIMIO   = TIMIO    + DTIME
C
      ENDIF
C
      IF (((.NOT. MINSCR).AND. (.NOT. CCS))
     *        .OR.(.NOT.(CCS.OR.CC2))) THEN
C
C-----------------------------------
C        Readin C2 amplitude in RHO.
C-----------------------------------
C
         DTIME = SECOND()
         CALL CC_RVEC(LUFC2,FC2AM,NT2AM(ISYMTR),NT2AM(ISYMTR),
     *                IVEC+IV-1,RHO2)
         DTIME = SECOND() - DTIME
         TIMIO = TIMIO + DTIME
C
C--------------------------------
C        Square up C2 amplitudes.
C--------------------------------
C
         CALL CCLR_DIASCL(RHO2,TWO,ISYMTR)
         CALL CC_T2SQ(RHO2,C2AM,ISYMTR)
C
         IF (IPRINT.GT.40) THEN
            CALL AROUND('CC_RHTR: (C1,C2) vector readin')
            CALL CC_PRSQ(C1AM(1,IV),C2AM,ISYMTR,1,1)
         ENDIF
C
C------------------------------
C        Read in result vector.
C------------------------------
C
         DTIME = SECOND()
         CALL CC_RVEC(LUFR2,FRHO2,NRHO2,NT2AM(ISYMTR),
     *                IV + ITR -1,RHO2)
         DTIME = SECOND() - DTIME
         TIMIO = TIMIO + DTIME
C
      ENDIF
C
      IF ( DEBUG ) THEN
         RHO1N = DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
         WRITE(LUPRI,1) 'Norm of Rho1 -after loop in mo:     ',RHO1N
         IF ( .NOT. CCS) THEN
            RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of Rho2 -after loop in mo:     ',RHO2N
         ENDIF
         IF (IPRINT.GT.40) THEN
            CALL AROUND(' RHO: ')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,0,1)
         END IF
      ENDIF
C
C-----------------------------------
C     For CC2 add Fock contribution.
C-----------------------------------
C
      IF (CC2.AND..NOT.(ONLY21)) THEN
C
         DTIME = SECOND()
         ISIDE = 1
         CALL CC2_FCK(RHO2,C2AM,WORK(KEND1),LWRK1,ISYMTR,
     *                WORK(KLAMDP),WORK(KLAMDH),ISIDE)
         DTIME = SECOND() - DTIME
         TIMFCK = TIMFCK + DTIME
C
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' AFTER CC2_FCK cont. CCLR:RHO ')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,0,1)
         ENDIF
C
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
      ISYFCK = ISYMTR
      ISYMPA = 1
      ISYMHO = 1
C
      IF (IPRINT .GT.140) THEN
         CALL AROUND( 'Fock AO matrix calc. from C1 transf. dens.' )
         CALL CC_PRFCKAO(WORK(KFCKC1),ISYFCK)
      ENDIF
C
      DTIME  = SECOND()
      CALL CC_FCKMO(WORK(KFCKC1),WORK(KLAMDP),WORK(KLAMDH),
     *              WORK(KEND1),LWRK1,ISYFCK,ISYMPA,ISYMHO)
      DTIME  = SECOND() - DTIME
      TIMFCKMO = DTIME  +  TIMFCKMO
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND( 'Fock MO matrix calc. from C1 transf. dens.' )
         CALL CC_PRFCKMO(WORK(KFCKC1),ISYFCK)
      ENDIF
C
C----------------------------------------------------------
C     Add ´F`ai = SUM-k-(Lai,kk-bar) To 1C1 transformation.
C----------------------------------------------------------
C
      IF (CCS) CALL DZERO(RHO1(1,IV),NT1AM(ISYMTR))
C
      DTIME = SECOND()
      CALL CCLR_1C1F(RHO1(1,IV),WORK(KFCKC1),ISYMTR)
      DTIME = SECOND() - DTIME
      TIM1C1F = TIM1C1F + DTIME
C
      IF (IPRINT .GT. 120) THEN
         CALL AROUND('After  1C1F:  CCLR:RHO ')
         CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,1,0)
      ENDIF
C
C-----------------------------------------------------
C     Make Fock tilde:
C     F~ai = SUM-k-(Lai,kk-bar) + Fa-bar,i + Fa,i-bar.
C-----------------------------------------------------
C
      IF ( .NOT. CCS ) THEN
C
         KFCKAO = KEND1
         KEND2  = KFCKAO + MAX(N2BST(ISYMOP),N2BST(ISYMTR))
         LEND2  = LWORK  - KEND2
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
C
C---------------------------------------------------
C     Readin E-intermediates: EI1(t1,t2),EI2(t1,t2).
C---------------------------------------------------
C
      KEI1   = KEND1
      KEI2   = KEI1 + NMATAB(ISYMOP)
      KEND2  = KEI2 + NMATIJ(ISYMOP)
      LWRK2  = LWORK - KEND2
C
      IF ( RSPIM.AND.((.NOT.(ONLY21)).AND.(.NOT.CCS))) THEN
C
         DTIME = SECOND() - DTIME
C
         LUE1 = 0
         CALL GPOPEN(LUE1,'CC_E1IM','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LUE1 )
         READ (LUE1 )(WORK(KEI1+ J-1),J = 1,NMATAB(ISYMOP))
         CALL GPCLOSE(LUE1,'KEEP')
C
         LUE2 = 0
         CALL GPOPEN(LUE2,'CC_E2IM','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LUE2 )
         READ (LUE2 )(WORK(KEI2+ J-1),J = 1,NMATIJ(ISYMOP))
         CALL GPCLOSE(LUE2,'KEEP')
C
         DTIME = SECOND() - DTIME
         TIMIO = TIMIO + DTIME
C
         IF (IPRINT.GT.150) THEN
            CALL AROUND( 'E-intermediates read from disk ')
            CALL CC_PREI(WORK(KEI1),WORK(KEI2),ISYMOP,1)
         ENDIF
C
      ENDIF
C
      IF (CCS) THEN
C
C---------------------------------------------------------
C        If ccs then put fock matrix into the EI matrices.
C---------------------------------------------------------
C
         FCKCON = .TRUE.
         ETRAN  = .FALSE.
         ISYMEI = 1
         DTIME = SECOND()
         CALL CCRHS_EFCK(WORK(KEI1),WORK(KEI2),WORK(KLAMDH),
     *                   WORK(KFOCK),WORK(KEND2),LWRK2,FCKCON,
     *                   ETRAN,ISYMEI)
         DTIME = SECOND() - DTIME
         TIMEFCK = TIMEFCK + DTIME
C
         IF (IPRINT.GT.150) THEN
            CALL AROUND( 'E-intermediates from from EFCK ')
            CALL CC_PREI(WORK(KEI1),WORK(KEI2),ISYMOP,1)
         ENDIF
C
      ENDIF
C
C-----------------------------------------------------
C     Calculate E contributions to rho1: C1*EI(t1,t2).
C-----------------------------------------------------
C
      IF (.NOT.(ONLY21)) THEN
C
         DTIME = SECOND()
         ISYMIM = 1
         CALL CCLR_E1C1(RHO1(1,IV),C1AM(1,IV),WORK(KEI1),WORK(KEI2),
     *                  WORK(KEND2),LWRK2,ISYMTR,ISYMIM,'N')
         DTIME = SECOND() - DTIME
         TIME1C1 = TIME1C1 + DTIME
C
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' CC_RHTR:  AFTER E1C1 EI*(C1) RHO1 ')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,1,0)
         ENDIF
C
      ENDIF
C
C-----------------------------------------------------------
C     Calculate first E contributions to rho2: C2*EI(t1,t2).
C-----------------------------------------------------------
C
      IF (.NOT.(CCS.OR.CC2)) THEN
C
         ISYVEC = ISYMTR
         ISYMIM = ISYMOP
         DTIME = SECOND()
         CALL CCRHS_E(RHO2,C2AM,WORK(KEI1),WORK(KEI2),
     *                WORK(KEND2),LWRK2,ISYVEC,ISYMIM)
         DTIME = SECOND() - DTIME
         TIME  = TIME + DTIME
C
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' CC_RHTR:  AFTER E EI*(C2) RHO2 ')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,0,NC2)
         ENDIF
C
      ENDIF
C
      IF ( DEBUG ) THEN
         RHO1N = DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
         WRITE(LUPRI,1) 'Norm of Rho1 -after Ei*(c1,c2):    ',RHO1N
         IF ( .NOT. CCS) THEN
           RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
           WRITE(LUPRI,1) 'Norm of Rho2 -after Ei*(c1,c2):    ',RHO2N
         ENDIF
      ENDIF
C
C----------------------------------------------------
C     Readin  Gamma(t1,t2) vector: CCLR intermediate.
C----------------------------------------------------
C
      IF ( RSPIM.AND..NOT.(CC2.OR.CCS)) THEN
C
         DTIME = SECOND()
C
         KGAMI  = KEND1
         KEND2  = KGAMI + NGAMMA(ISYMOP)
         LWRK2  = LWORK - KEND2
         LUGAM = 0
         CALL GPOPEN(LUGAM,'CC_GAMIM','UNKNOWN',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND(LUGAM)
         READ (LUGAM)(WORK(KGAMI+J-1),J = 1,NGAMMA(ISYMOP))
         CALL GPCLOSE(LUGAM,'KEEP')
C
         DTIME = SECOND() - DTIME
         TIMIO = TIMIO + DTIME
C
C------------------------------------------
C        Calculate A-term: C2*GAMMA(t1,t2).
C------------------------------------------
C
         ISYGAM = 1
         ISYVEC = ISYMTR
         IOPT    = 1
C
         DTIME  = SECOND()
         CALL CCRHS_A(RHO2,C2AM,WORK(KGAMI),WORK(KEND2),LWRK2,
     *                ISYGAM,ISYVEC,IOPT)
         DTIME  = SECOND() - DTIME
         TIMA   = DTIME  + TIMA
C
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' AFTER C2*GAMMA(t1,t2) A-term RHO is ')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,0,1)
         ENDIF
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of Rho1 -after C2*GAM:        ',RHO1N
            WRITE(LUPRI,1) 'Norm of Rho2 -after C2*GAM:        ',RHO2N
         ENDIF
C
      ENDIF
C
C----------------------------------------
C     Calculate the C-term: C2*CI(t1,t2).
C----------------------------------------
C
      IF  (.NOT. (CCS .OR. CC2)) THEN
C
         ISYCIM = 1
C
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' CC_RHTR:  AFTER EI*(C1,C2) RHO= ')
            CALL CC_PRSQ(RHO1(1,IV),RHO2,ISYMTR,1,NC2)
         ENDIF
C
         IOPT  = 1
         IVECNR = 1
         DTIME  = SECOND()
         CALL CCSD_T2TP(C2AM,WORK(KEND2),LWRK2,ISYMTR)
         CALL CCRHS_CIO(RHO2,C2AM,WORK(KLAMDH),WORK(KEND2),LWRK2,
     *                  ISYMTR,ISYCIM,LUCIM,CFIL,IVECNR,IOPT)
         CALL CCSD_T2TP(C2AM,WORK(KEND2),LWRK2,ISYMTR)
         DTIME = SECOND() - DTIME
         TIMCIO  = DTIME + TIMCIO
C
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' AFTER C2*CIM(t1,t2) C-term RHO is ')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,0,1)
         ENDIF
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of Rho1 -after C2*CIM:        ',RHO1N
            WRITE(LUPRI,1) 'Norm of Rho2 -after C2*CIM:        ',RHO2N
         ENDIF
C
      ENDIF
C
C------------------------------
C     Transform C2 to 2C2 - C2.
C------------------------------
C
      IF ( .NOT. CCS ) THEN
         DTIME    = SECOND()
         CALL CCRHS_T2TR(C2AM,WORK(KEND1),LWRK1,ISYMTR)
         DTIME    = SECOND() - DTIME
         TIMT2TR  = TIMT2TR + DTIME
      ENDIF
C
C----------------------------------------
C     Calculate the D-term: C2*DI(t1,t2).
C----------------------------------------
C
      IF  (.NOT. (CCS .OR. CC2)) THEN
C
         ISYDIM = 1
C
         IOPT = 1
         IVECNR = 1
         DTIME = SECOND()
         CALL CCRHS_DIO(RHO2,C2AM,WORK(KLAMDH),WORK(KEND2),LWRK2,
     *                  ISYMTR,ISYDIM,LUDIM,DFIL,IVECNR,IOPT)
         DTIME = SECOND() - DTIME
         TIMDIO  = DTIME + TIMDIO
C
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' AFTER C2*DIM(t1,t2) D-term RHO is ')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,0,1)
         ENDIF
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of Rho1 -after C2*DIM:        ',RHO1N
            WRITE(LUPRI,1) 'Norm of Rho2 -after C2*DIM:        ',RHO2N
         ENDIF
C
      ENDIF
C
C----------------------
C     Calculate I-term.
C----------------------
C
      IF ( .NOT. CCS ) THEN
         DTIME    = SECOND()
         CALL CCRHS_I(RHO1(1,IV),C2AM,WORK(KFOCK),WORK(KEND1),LWRK1,
     *                ISYMTR,1)
         DTIME    = SECOND()-DTIME
         TIMI     = DTIME + TIMI
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' AFTER I-term :RHO ')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,1,0)
         ENDIF
      ENDIF
C
C---------------------------------------------------------------------
C        Readin  BF(t1,t2) and contract with lampc and lamdp  and add.
C        NB This destroys C2AM.
C---------------------------------------------------------------------
C
      IF ( .NOT. (CC2 .OR. CCS)) THEN
C
         IF ( RSPIM ) THEN
C
            DTIME = SECOND()
C
            LUBF = 0
            CALL GPOPEN(LUBF,'CC_BFIM','OLD',' ','UNFORMATTED',IDUMMY,
     &                  .FALSE.)
            READ(LUBF) (C2AM(J),J = 1,2*NT2ORT(1))
            CALL GPCLOSE(LUBF,'KEEP')
C
            DTIME = SECOND() - DTIME
            TIMIO = TIMIO + DTIME
C
         ENDIF
C
         ISYMBF = 1
         NEWGS  = NEWGAM
         NEWGAM = .FALSE.
         ICON = 2
C
         DTIME = SECOND()
         CALL CC_T2MO(FAKE,PHONEY,ISYMOP,C2AM,C2AM(1+2*NT2ORT(1)),DUM,
     *                WORK(KLAMDP),WORK(KLPC1),ISYMTR,
     *                WORK(KEND1),LWRK1,ISYMBF,ICON)
         CALL DAXPY(NT2AM(ISYMTR),ONE,C2AM(1+2*NT2ORT(1)),1,RHO2,1)
         DTIME = SECOND() -DTIME
         TIMT2MO = TIMT2MO + DTIME
C
         NEWGAM = NEWGS
C
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' AFTER C1*BFIM(t1,t2) BF-tilde term RHO is ')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,0,1)
         ENDIF
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of Rho1 -after C1*BFIM:       ',RHO1N
            WRITE(LUPRI,1) 'Norm of Rho2 -after C1*BFIM:       ',RHO2N
         ENDIF
C
      ENDIF
C
C-----------------------------
C     write out result vector.
C-----------------------------
C
      CALL CC_WVEC(LUFR1,FRHO1,NT1AM(ISYMTR),NT1AM(ISYMTR),
     *             IV + IVEC -1,RHO1(1,IV))
      IF (.NOT. CCS) THEN
         DTIME   = SECOND()
         CALL CC_WVEC(LUFR2,FRHO2,NRHO2,NT2AM(ISYMTR),
     *                IV + ITR -1,RHO2)
         DTIME   = SECOND() - DTIME
         TIMIO   = TIMIO    + DTIME
      ENDIF
C
  150 CONTINUE
C
C==========================================================
C     MO basis T2SQ section.
C     Contract intermediates constructed in loop with T2SQ.
C==========================================================
C
C
C----------------------------------------------
C     Loop over number of trial vectors nsimtr.
C----------------------------------------------
C
      DO 160 IV = 1, NSIMTR
C
         IF (CCS) GOTO 160
C
         KLPC1  = KLAMPC + NGLMDT(ISYMTR)*(IV-1)
         KLHC1  = KLAMHC + NGLMDT(ISYMTR)*(IV-1)
         KFCKC1 = KFOCKC + N2BST(ISYMTR)*(IV-1)
         KEI1   = KEMAT1 + NEMAT1(ISYMTR)*(IV-1)
         KEI2   = KEMAT2 + NMATIJ(ISYMTR)*(IV-1)
         KGAMC1 = KGAMMC + NGAMMA(ISYMTR)*(IV-1)
C
C------------------------------------
C        Readin T2 amplitude in RHO2.
C------------------------------------
C
         DTIME = SECOND()
C
         IOPT = 2
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,DUMMY,RHO2)
C
         DTIME = SECOND() - DTIME
         TIMIO = TIMIO + DTIME
C
C--------------------------------
C        Square up T2 amplitudes.
C--------------------------------
C
         DTIME = SECOND()
         CALL CC_T2SQ(RHO2,C2AM,1)
         DTIME = SECOND() - DTIME
         TIMT2SQ = TIMT2SQ + DTIME
C
         IF (IPRINT.GT.50) THEN
            CALL AROUND('CC_RHTR: (T1,T2) vector readin')
            CALL CC_PRSQ(WORK(KT1AM),C2AM,1,1,1)
         ENDIF
C
C------------------------------
C        Read in result vector.
C------------------------------
C
         DTIME = SECOND()
         CALL CC_RVEC(LUFR2,FRHO2,NRHO2,NT2AM(ISYMTR),
     *                IV + ITR  -1,RHO2)
         DTIME = SECOND() - DTIME
         TIMIO = TIMIO + DTIME
C
C----------------------------------------------------
C     Calculate E-term: T2*E(C1,C2)
C----------------------------------------------------
C
      IF ( .NOT.(CC2.OR.CCS)) THEN
C
         IF (IPRINT.GT.150) THEN
            CALL AROUND( 'E-intermdiates calc. EI(C1,C2) - ao.')
            CALL CC_PREI(WORK(KEI1),WORK(KEI2),ISYMTR,0)
         ENDIF
C
         ETRAN  = .TRUE.
         FCKCON = .TRUE.
         ISYMEI = ISYMTR
         DTIME = SECOND()
         CALL CCRHS_EFCK(WORK(KEI1),WORK(KEI2),WORK(KLAMDH),
     *                   WORK(KFCKC1),WORK(KEND1),LWRK1,FCKCON,
     *                   ETRAN,ISYMTR)
         DTIME = SECOND() - DTIME
         TIMEFCK = TIMEFCK + DTIME
C
         IF (IPRINT.GT.150) THEN
            CALL AROUND( 'E-intermdiates calc EI(C1,C2) -mo' )
            CALL CC_PREI(WORK(KEI1),WORK(KEI2),ISYMTR,1)
         ENDIF
C
         ISYVEC = 1
         ISYMIM = ISYMTR
         DTIME = SECOND()
         CALL CCRHS_E(RHO2,C2AM,WORK(KEI1),WORK(KEI2),
     *                WORK(KEND1),LWRK1,ISYVEC,ISYMTR)
         DTIME = SECOND() - DTIME
         TIME  = TIME + DTIME
C
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' AFTER CCRHS_E 1.cont.  CCLR:RHO ')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,0,1)
         ENDIF
C
      ENDIF
C
      IF (CC2 .AND. (NFIELD .GT. 0)) THEN
         CALL CC_RCC2FF(RHO2,C2AM,WORK(KLAMDP),WORK(KLAMDH),
     *                  WORK(KLPC1),WORK(KLHC1),
     *                  WORK(KEND1),LWRK1,ISYMTR)
      ENDIF
C
      IF ( DEBUG ) THEN
         RHO1N = DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
         RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
         WRITE(LUPRI,1) 'Norm of Rho1 -after T2*EI:         ',RHO1N
         WRITE(LUPRI,1) 'Norm of Rho2 -after T2*EI:         ',RHO2N
      ENDIF
C
C-----------------------------------
C     Calculate T2*GAMMA(C2) A-term.
C-----------------------------------
C
      IF (.NOT.(CCS .OR. CC2)) THEN
C
         ISYGAM  = ISYMTR
         ISYVEC  = 1
         IOPT    = 1
C
         DTIME   = SECOND()
         CALL CCRHS_A(RHO2,C2AM,WORK(KGAMC1),WORK(KEND1),LWRK1,
     *                ISYGAM,ISYVEC,IOPT)
         DTIME   = SECOND() - DTIME
         TIMA    = DTIME + TIMA
C
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' AFTER T2*GAMMA(C1,C2) A-term RHO ' )
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,0,1)
         ENDIF
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of Rho1 -after T2*GAM:        ',RHO1N
            WRITE(LUPRI,1) 'Norm of Rho2 -after T2*GAM:        ',RHO2N
         ENDIF
C
      ENDIF
C
C------------------------------------
C     Recover work from gamma matrix.
C     IF not minscr!
C------------------------------------
C
      IF ( MINSCR) THEN
         KEND1 = KENDS2
         LWRK1 = LWRKS2
      ENDIF
C
C-------------------------------------
C     Calculate the C-term: T2*CI(C1).
C-------------------------------------
C
      IF  (.NOT. (CCS .OR. CC2)) THEN
C
         ISYCIM = ISYMTR
         ISYVEC = 1
C
         IOPT = 2
         DTIME  = SECOND()
         CALL CCSD_T2TP(C2AM,WORK(KEND1),LWRK1,ISYVEC)
         CALL CCRHS_CIO(RHO2,C2AM,WORK(KLAMDH),WORK(KEND1),LWRK1,
     *                  ISYVEC,ISYCIM,LUC,CTFIL,IV,IOPT)
         CALL CCSD_T2TP(C2AM,WORK(KEND1),LWRK1,ISYVEC)
         DTIME   = SECOND() - DTIME
         TIMCIO  = DTIME + TIMCIO
C
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' AFTER T2*CIM(C1) C-term RHO is ')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,0,1)
         ENDIF
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of Rho1 -after T2*CI:         ',RHO1N
            WRITE(LUPRI,1) 'Norm of Rho2 -after T2*CI:         ',RHO2N
         ENDIF
C
      ENDIF
C
C------------------------------
C     Transform T2 to 2T2 - T2.
C------------------------------
C
      IF ( .NOT. CCS) THEN
         DTIME   = SECOND()
         CALL CCRHS_T2TR(C2AM,WORK(KEND1),LWRK1,1)
         DTIME   = SECOND() - DTIME
         TIMT2TR = TIMT2TR  + DTIME
      ENDIF
C
C-------------------------------
C     Calculate (2T2-T2)*DI(C1).
C-------------------------------
C
      IF  (.NOT. (CCS .OR. CC2)) THEN
C
         ISYDIM = ISYMTR
         ISYVEC = 1
         IOPT = 2
C
         DTIME   = SECOND()
         CALL CCRHS_DIO(RHO2,C2AM,WORK(KLAMDH),WORK(KEND1),LWRK1,
     *                  ISYVEC,ISYDIM,LUD,DTFIL,IV,IOPT)
         DTIME   = SECOND() - DTIME
         TIMDIO  = DTIME + TIMDIO
C
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' AFTER (2T2-T2)*DI(C1) D-term RHO is ')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,0,1)
         ENDIF
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of Rho1 -after T2*DI:         ',RHO1N
            WRITE(LUPRI,1) 'Norm of Rho2 -after T2*DI:         ',RHO2N
         ENDIF
C
      ENDIF
C
C----------------------------
C     Calculate I-tilde-term.
C----------------------------
C
      IF ( .NOT. CCS ) THEN
         DTIME    = SECOND()
         CALL CCRHS_I(RHO1(1,IV),C2AM,WORK(KFCKC1),WORK(KEND1),LWRK1,
     *                1,ISYMTR)
         DTIME   = SECOND() - DTIME
         TIMI    = DTIME + TIMI
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' AFTER CCLR_I tilde contribution: RHO ')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,1,0)
         ENDIF
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of Rho1 -after I-tilde:       ',RHO1N
            WRITE(LUPRI,1) 'Norm of Rho2 -after I-tilde:       ',RHO2N
         ENDIF
C
      ENDIF
C
C-------------------------------------------
C     SLV98,OC
C     Calculate T^{gC} solvent contribution.
C     Unsquared T2 is required!
C-------------------------------------------
C
      IF (CCSLV.AND.LSTBTR.AND.(.NOT.JACTST)) THEN
         CALL CC_RVEC(LUFC2,FC2AM,NT2AM(ISYMTR),NT2AM(ISYMTR),
     *                IVEC+IV-1,C2AM)
         LR = 'R'
         CALL CCSL_LTRB(RHO1(1,IV),RHO2,C1AM(1,IV),
     *                  C2AM,ISYMTR,LR,WORK(KEND1),LWRK1)
      ENDIF
C
C---------------------------------------------------------
C        And finally scale diagonal (if not triples cont.)
C---------------------------------------------------------
C
      IF (.NOT. CCSDT) CALL CCLR_DIASCL(RHO2,XHALF,ISYMTR)
C
      CALL CC_WVEC(LUFR1,FRHO1,NT1AM(ISYMTR),NT1AM(ISYMTR),
     *             IV + IVEC -1,RHO1(1,IV))
      IF (.NOT. CCS   ) THEN
C
C--------------------------------
C        write out result vector.
C--------------------------------
C
         DTIME   = SECOND()
         CALL CC_WVEC(LUFR2,FRHO2,NRHO2,NT2AM(ISYMTR),
     *                IV + ITR -1,RHO2)
         DTIME   = SECOND() - DTIME
         TIMIO   = TIMIO    + DTIME
C
      ENDIF
C
C--------------------------------------
C     End loop over vectors in scratch.
C--------------------------------------
C
  160 CONTINUE
C
      KEND1 = KENDS2
      LWRK1 = LWRKS2
C
C=============================================================
C
C     If CCSDT = true
C     Partitioned triples contributions to the transformation.
C
C=============================================================
C
c     only for one vector at a time in triples code.
c
      IV = 1
      IF (CCSDT) THEN
 
       TIMTRIP = SECOND()

       IF (JACTST) ECURR = 0.0D00

       IF (NODDY_RHTR) THEN

         KFCKC1 = KFOCKC + N2BST(ISYMTR)*(IV-1)
         CALL CC_RHTR_NODDY(IVEC+IV-1,LUFC2,FC2AM,ECURR,
     &                      C1AM(1,IV),WORK(KFCKC1),ISYMTR,
     &                      WORK(KFOCK),WORK(KLAMDP),WORK(KLAMDH),
     &                      RHO1(1,IV),RHO2,WORK(KEND1),LWRK1)

       ELSE
C
C---------------------------------------
C        Readin C2 amplitudes from disk.
C---------------------------------------
C
         KT2AM = KEND1
         KEND2 = KT2AM + MAX(NT2SQ(1),NT2SQ(ISYMTR),NT1AM(ISYMTR))
         LWRK2 = LWORK - KEND2
C
         DTIME = SECOND()
         CALL CC_RVEC(LUFC2,FC2AM,NT2AM(ISYMTR),NT2AM(ISYMTR),
     *                IVEC+IV-1,WORK(KT2AM))
         DTIME = SECOND() - DTIME
         TIMIO = TIMIO + DTIME
C
C---------------------------------
C        Prepare the C-amplitudes.
C---------------------------------
C
         CALL CCLR_DIASCL(WORK(KT2AM),TWO,ISYMTR)
         CALL CC_T2SQ(WORK(KT2AM),C2AM,ISYMTR)
C
C-----------------------------------------------------------
C        Calculate terms depending on C2 occuring in energy.
C        <mui|[H,T3(C2)]|hf>
C        T3(C2) = ome-1*< |[H,C2]| >
C-----------------------------------------------------------
C
         CC3LR = .FALSE.
         ISYMT1 = 1
         ISYMT2 = ISYMTR
         DTIME   = SECOND()
         CALL CC3_OMEG(ECURR,RHO1(1,IV),RHO2,C1AM(1,IV),ISYMT1,C2AM,
     *                 ISYMT2,WORK(KFOCK),WORK(KLAMDP),WORK(KLAMDH),
     *                 WORK(KEND1),LWRK1,LU3SRTR,FN3SRTR,LUDELDR,
     *                 FNDELDR,LUCKJDR,FNCKJDR,LUDKBCR,FNDKBCR,
     *                 LUTOC,FNTOC,LU3VI,FN3VI,LU3VI2,FN3VI2)
C
         DTIME   = SECOND() - DTIME
         TIMOM31 = TIMO31 + DTIME
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' RHO after <mui|[H,T3(C2)]|hf>  ')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,1,1)
         ENDIF
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of Rho1 -after cc3_omeg-1:    ',RHO1N
            WRITE(LUPRI,1) 'Norm of Rho2 -after cc3_omeg-1:    ',RHO2N
         ENDIF
C
C-------------------------------------------------
C        If CC1A no t1 in triples part, else
C        avoid also for ccsdr(t)
C-------------------------------------------------
C
       IF (.NOT. CC1A ) THEN
C
        IF (.NOT. CCRT) THEN
C
C------------------------------
C        Readin T2SQ from disk.
C------------------------------
C
         DTIME = SECOND()
C
         IOPT = 2
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KT2AM))
C
         DTIME = SECOND() - DTIME
         TIMIO = TIMIO + DTIME
C
         DTIME = SECOND()
         CALL CC_T2SQ(WORK(KT2AM),C2AM,1)
         DTIME = SECOND() - DTIME
         TIMT2SQ = TIMT2SQ + DTIME
C
C-------------------------------------------------
C        Calculate c1 term not occuring in energy.
C        <mu2|[[H,C1],T3(t2)]|hf>
C         T2am is in C2AM.
C-------------------------------------------------
C
         CC3LR   = .TRUE.
         ECURR2  = ECURR
         ECURR   = 0.0D0
         ISYMT1  = ISYMTR
         ISYMT2  = 1
         DTIME   = SECOND()
         CALL CC3_OMEG(0.0D0,RHO1(1,IV),RHO2,C1AM(1,IV),ISYMT1,C2AM,
     *                 ISYMT2,WORK(KFOCK),WORK(KLAMDP),WORK(KLAMDH),
     *                 WORK(KEND1),LWRK1,LU3SRTR,FN3SRTR,LUDELDR,
     *                 FNDELDR,LUCKJDR,FNCKJDR,LUDKBCR,FNDKBCR,
     *                 LUTOC,FNTOC,LU3VI,FN3VI,LU3VI2,FN3VI2)
         DTIME   = SECOND() - DTIME
         TIMOM32 = TIMO32 + DTIME
         ECURR   = ECURR2
         CC3LR   = .FALSE.
C
         IF (IPRINT .GT. 120) THEN
            CALL AROUND(' RHO after <mu2|[[H,C1],T3(t2)]|hf> ')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,1,1)
         ENDIF
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of Rho1 -after cc3_omeg-2:    ',RHO1N
            WRITE(LUPRI,1) 'Norm of Rho2 -after cc3_omeg-2:    ',RHO2N
         ENDIF
C
        ENDIF
C
C--------------------------------------------------
C        Calculate Transformed integrals with C1.
C        Not needed if CC1A or CC1B - only for CC3.
C--------------------------------------------------
C
        IF (.NOT.CC1B) THEN
C
C==============================================================
C        Start the second loop over distributions of integrals.
C==============================================================
C
         KENDS2 = KEND1
         LWRKS2 = LWRK1
C
         IF (DIRECT) THEN
            NTOSYM = 1
            DTIME  = SECOND()
C
            KCCFB1 = KEND1
            KINDXB = KCCFB1 + MXPRIM*MXCONT
            KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
            LWRK1  = LWORK  - KEND1
C
            CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     &                  KODPP1,KODPP2,KRDPP1,KRDPP2,
     &                  KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     &                  WORK(KEND1),LWRK1,IPRERI)
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
         IF (IPRINT.GT.50) THEN
            WRITE(LUPRI,*) '  '
            WRITE(LUPRI,*) 'IN CC_RHTR -2.loop: NTOSYM = ',NTOSYM
         ENDIF
         ICDEL1 = 0
         ICDEL2 = 0
         DO 200 ISYMD1 = 1,NTOSYM
C
            IF (DIRECT) THEN
               NTOT = MXCALL
            ELSE
               NTOT = NBAS(ISYMD1)
            ENDIF
C
            DO 210 ILLL = 1,NTOT
C
C------------------------------------------------
C              If direct calculate the integrals.
C------------------------------------------------
C
               IF (DIRECT) THEN
C
                  KEND1 = KENDSV
                  LWRK1 = LWRKSV
C
                  DTIME  = SECOND()
                  CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     &                        WORK(KODCL1),WORK(KODCL2),
     &                        WORK(KODBC1),WORK(KODBC2),
     &                        WORK(KRDBC1),WORK(KRDBC2),
     &                        WORK(KODPP1),WORK(KODPP2),
     &                        WORK(KRDPP1),WORK(KRDPP2),
     &                        WORK(KCCFB1),WORK(KINDXB),
     &                        WORK(KEND1), LWRK1,IPRERI)
                  DTIME  = SECOND() - DTIME
                  TIMER2 = TIMER2 + DTIME
C
                  KRECNR = KEND1
                  KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
                  LWRK1  = LWORK  - KEND1
                  IF (LWRK1 .LT. 0) THEN
                     CALL QUIT('Insufficient core in CC_RHTR')
                  END IF
C
               ELSE
                  NUMDIS = 1
               ENDIF
C
C--------------------------------------------------------
C              Loop over number of distributions in disk.
C--------------------------------------------------------
C
               DO 220 IDEL2 = 1,NUMDIS
C
                  IF (DIRECT) THEN
                     IDEL  = INDEXA(IDEL2)
                     ISYMD = ISAO(IDEL)
                  ELSE
                     IDEL  = IBAS(ISYMD1) + ILLL
                     ISYMD = ISYMD1
                  ENDIF
C
                  ISYDIS = MULD2H(ISYMD,ISYMOP)
                  ISYAIK = MULD2H(ISYDIS,ISYMTR)
C
C---------------------------------------------
C                 Work space allocation no. 2.
C---------------------------------------------
C
                  KXINT  = KEND1
                  KEND2  = KXINT + NDISAO(ISYDIS)
                  LWRK2  = LWORK - KEND2
C
                  IF ( IPRINT .GT. 55)
     *            WRITE(LUPRI,*) ' In CC_RHTR 2. alloc. work left:',
     *            LWRK2
C
                  IF (LWRK2 .LT. 0) THEN
                     WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                     CALL QUIT('Insufficient space in CCRHSN')
                  ENDIF
C
C--------------------------------------------
C                 Read in batch of integrals.
C--------------------------------------------
C
                  DTIME   = SECOND()
                  CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                        WORK(KRECNR),DIRECT)
                  DTIME   = SECOND() - DTIME
                  TIMRDAO = TIMRDAO  + DTIME
C
C--------------------------------------------------------------
C                 Calculate transformed integrals used in t3am.
C--------------------------------------------------------------
C
                  DTIME = SECOND()
                  CALL CC3_T3INT(WORK(KXINT),WORK(KLAMDP),
     *                           WORK(KLAMDH),C1AM(1,IV),ISYMTR,
     *                           WORK(KEND2),LWRK2,IDEL,ISYMD,2,
     *                           LU3SRTR,FN3SRTR,LUCKJDR,FNCKJDR)
                  DTIME = SECOND() - DTIME
                  TIMT3I = TIMT3I + DTIME
C
  220          CONTINUE
  210       CONTINUE
  200    CONTINUE
C
C---------------------------
C        Recover work space.
C---------------------------
C
         KEND1 = KENDS2
         LWRK1 = LWRKS2
C
C------------------------------
C        Readin T2SQ from disk.
C------------------------------
C
         KT2AM = KEND1
         KEND2 = KT2AM + MAX(NT2SQ(1),NT2SQ(ISYMTR),NT1AM(ISYMTR))
         LWRK2 = LWORK - KEND2
C
         DTIME = SECOND()
C
         IOPT = 2
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KT2AM))
C
         DTIME = SECOND() - DTIME
         TIMIO = TIMIO + DTIME
C
         DTIME = SECOND()
         CALL CC_T2SQ(WORK(KT2AM),C2AM,1)
         DTIME = SECOND() - DTIME
         TIMT2SQ = TIMT2SQ + DTIME
C
C-------------------------------------------
C        Calculate C1 in t3 terms.
C        <mui|[H,T3(C1,T2)]|hf>
C        T3(C1,T2) = ome-1*< |[[H,C1],T2]| >
C        T2 is in C2AM
C-------------------------------------------
C
         CC3LR = .FALSE.
         ISYMT1 = ISYMTR
         ISYMT2 = 1
         DTIME = SECOND()
         CALL CC3_OMEG(ECURR,RHO1(1,IV),RHO2,WORK(KT1AM),ISYMT1,C2AM,
     *                 ISYMT2,WORK(KFOCK),WORK(KLAMDP),WORK(KLAMDH),
     *                 WORK(KEND1),LWRK1,LU3SRTR,FN3SRTR,LUDELDR,
     *                 FNDELDR,LUCKJDR,FNCKJDR,LUDKBCR,FNDKBCR,
     *                 LUTOC,FNTOC,LU3VI,FN3VI,LU3VI2,FN3VI2)
         DTIME = SECOND() - DTIME
         TIMOM33 = TIMO33 + DTIME
C
         IF (IPRINT .GT. 120) THEN
            CALL AROUND('RHO after <mui|[H,T3(C1,T2)]|hf> ')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,1,1)
         ENDIF
C
         IF ( DEBUG ) THEN
            RHO1N = DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of Rho1 -after cc3_omeg-3:    ',RHO1N
            WRITE(LUPRI,1) 'Norm of Rho2 -after cc3_omeg-3:    ',RHO2N
         ENDIF
C
        ENDIF! (.NOT.CC1B)
       ENDIF ! (.NOT. CC1A)
       ENDIF ! (NODDY_RHTR)
C
C---------------------
C      Scale Diagonal.
C---------------------
C
       CALL CCLR_DIASCL(RHO2,XHALF,ISYMTR)
C
C------------------------------
C      write out result vector.
C------------------------------
C
       DTIME = SECOND()
       CALL CC_WVEC(LUFR1,FRHO1,NT1AM(ISYMTR),NT1AM(ISYMTR),
     *              IV + IVEC -1,RHO1(1,IV))
       CALL CC_WVEC(LUFR2,FRHO2,NRHO2,NT2AM(ISYMTR),
     *              IV + ITR  -1,RHO2)
       DTIME = SECOND() - DTIME
       TIMIO = TIMIO  + DTIME
C
C-----------------------------
C      END OF TRIPLES SECTION.
C-----------------------------
C
       TIMTRIP = SECOND() - TIMTRIP
C
      ENDIF
C
C=============================================================
C
C     End section.
C
C=============================================================
C
      IF (IPRINT .GT. 50 ) THEN
         CALL AROUND('END OF CC_RHTR :RHO ')
         NC2 = 1
         IF ( CCS ) NC2 = 0
         CALL CC_PRP(RHO1(1,IV),RHO2,ISYMTR,1,NC2)
      ENDIF
C
      IF ( DEBUG ) THEN
         RHO1N = DDOT(NT1AM(ISYMTR),RHO1(1,IV),1,RHO1(1,IV),1)
         WRITE(LUPRI,1) 'Norm of Rho1 -end of cclr_rhtr:    ',RHO1N
         IF ( .NOT. CCS) THEN
            RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of Rho2 -end of cclr_rhtr:    ',RHO2N
         ENDIF
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
     *             TIMEFCK + TIME1C1  + TIM1C1F + TIMCIO  +
     *             TIMDIO  + TIMLAM   + TIMAOD  + TIMTRBT +
     *             TIMT2AO + TIMTCME  + TIMT2TR + TIMT2BT +
     *             TIMT2MO + TIMFCKMO + TIMFCK  + TIMT2SQ
         WRITE(LUPRI,'(/)')
         WRITE(LUPRI,9998) 'CC_RHTR in total ', TIMALL
         WRITE(LUPRI,9998) 'CC part            ', TIMCCSD
         WRITE(LUPRI,9998) 'Int. calc. & read  ', TIMER1 + TIMER2
     *                  + TIMRDAO
         WRITE(LUPRI,9998) 'IO of vectors/I.M. ', TIMIO
         IF (CCSDT) THEN
            WRITE(LUPRI,9998) 'Triples corr.      ',TIMTRIP
         ENDIF
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
         ENDIF
         IF ( .NOT. CCS ) THEN
            WRITE(LUPRI,9999) 'E               ', TIME
            IF (CC2) WRITE(LUPRI,9999) 'F               ', TIMF
            WRITE(LUPRI,9999) 'G               ', TIMG
            WRITE(LUPRI,9999) 'H               ', TIMH
            WRITE(LUPRI,9999) 'I               ', TIMI
            WRITE(LUPRI,9999) 'EI              ', TIMEI
         ENDIF
         WRITE(LUPRI,9999) 'EFCK            ', TIMEFCK
         WRITE(LUPRI,9999) 'E1C1            ', TIME1C1
         WRITE(LUPRI,9999) '1C1F            ', TIM1C1F
         WRITE(LUPRI,9999) 'LAMTRA          ', TIMLAM
         WRITE(LUPRI,9999) 'AODENS          ', TIMAOD
         WRITE(LUPRI,9999) 'FCK             ', TIMFCK
         WRITE(LUPRI,9999) 'TRBT            ', TIMTRBT
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
         CALL AROUND(' END OF CC_RHTR ')
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
      IF (DUMPCD.AND.(.NOT.(CCS.OR.CC2))) THEN
         CALL WCLOSE2(LUC,CTFIL,'DELETE')
         CALL WCLOSE2(LUD,DTFIL,'DELETE')
C
         CALL WCLOSE2(LUCIM,CFIL,'KEEP')
         CALL WCLOSE2(LUDIM,DFIL,'KEEP')
C
      ENDIF
C
      IF (CCSDT) THEN
         CALL WCLOSE2(LU3SRTR,FN3SRTR,'DELETE')
         CALL WCLOSE2(LUDELDR,FNDELDR,'DELETE')
         CALL WCLOSE2(LUCKJDR,FNCKJDR,'DELETE')
         CALL WCLOSE2(LUDKBCR,FNDKBCR,'DELETE')
         CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
         CALL WCLOSE2(LU3VI,FN3VI,'KEEP')
         CALL WCLOSE2(LU3VI2,FN3VI2,'KEEP')
      ENDIF
C
C-------------------------------
C     Restore logical variables.
C-------------------------------
C
      CC1B  = CC1BSA
C
      RETURN
      END
