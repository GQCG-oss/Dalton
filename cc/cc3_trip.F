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
C  /* Deck cc3_triplet */
      SUBROUTINE CC3_TRIPLET(ECURR,ISYMTR,OMEGA1,OMEGA2P,T2AM,
     *                       C2AMP,C2AMM,FOCK,
     *                       XLAMDP,XLAMDH,XLAMPC,XLAMHC,WORK,LWORK,
     *                       LUFR2,FRHO2,IV,ITR,C1AM)
!
!     Written by Kasper Hald
!
!     Purpose:  Calculate the CC3 corrections to the CCSD
!               triplet excitation energies.
!
!     NB: The T2 amplitudes are assumed to be a full square,
!         C2AMP is (C2+) and C2AMM is (C2-).
!         The C2AMP and C2AMM are also assumed to be full squares
!
!     NB: It is assumed that the vectors are allocated in the following
!     order:
!         OMEGA1(*), OMEGA2(*), T2AM(*), SCR(*), WRK(*).
!
!     NB: This is ONLY a noddy code ... no intermediates etc.
!
      IMPLICIT NONE
!
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <maxash.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
#include <eribuf.h>
#include <ccorb.h>
#include <inftap.h>
#include <ccisao.h>
#include <r12int.h>
!
      INTEGER INDEX, LWORK, KCMO, KSCR1, KFOCKD, KINT1T
      INTEGER KINT2T, KINT3T, KINT4T, KINT1S, KINT2S, KINT3S, KINT4S
      INTEGER KINT5S, KINT6S, KT3AM, KOME1, KOME2P, KOME2M, KIAJB
      INTEGER KYIAJB, KEND1, LWRK1, NTOSYM
      INTEGER ISYMD1, NTOT, ILLL, NUMDIS, IDEL2, IDEL, ISYMD
      INTEGER ISYDIS, KXINT, KEND2, LWRK2, LUFR2, KOME2, IV, ITR
      INTEGER IJ, NIJ, KRECNR, ISYMTR
      INTEGER INDEXA(MXCORB)
!
      INTEGER KOME22, IERRR1, IERRR2
      INTEGER LUCC3R1, LUCC3R2, TEMP
      INTEGER KINTTOT, KINTTOT2
      INTEGER KOFF1, KOFF2
      INTEGER KINTTMP1, KINTTMP2, KTRTMP1
      INTEGER KTRTMP2, KTRTMP3, KTRTMP4
!
      INTEGER IDUM
!
#if defined (SYS_CRAY)
      REAL ZERO, HALF, ONE, TWO, DTIME ! , TIMHER1, TIMHER2
      REAL TIMRDAO, OMEGA1(*), OMEGA2P(*), T2AM(*)
      REAL FOCK(*), XLAMDP(*), XLAMDH(*), XLAMPC(*)
      REAL XLAMHC(*), WORK(LWORK), C2AMP(*), C2AMM(*)
      REAL DDOT, RHO2N
      REAL XONE, XHALF, ECURR
      REAL C1AM(*)
#else
      DOUBLE PRECISION ZERO, HALF, ONE, TWO, DTIME ! , TIMHER1, TIMHER2
      DOUBLE PRECISION TIMRDAO, OMEGA1(*), OMEGA2P(*), T2AM(*)
      DOUBLE PRECISION FOCK(*), XLAMDP(*), XLAMDH(*), XLAMPC(*)
      DOUBLE PRECISION XLAMHC(*), WORK(LWORK), C2AMP(*), C2AMM(*)
      DOUBLE PRECISION DDOT, RHO2N
      DOUBLE PRECISION XONE, XHALF, ECURR
      DOUBLE PRECISION C1AM(*)
#endif
!
      LOGICAL LOCDBG
!
      CHARACTER*8 FRHO2
      CHARACTER*8 RHO1FIL, RHO2FIL
!
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (XONE = -1.0D00, XHALF = -0.5D0)
      PARAMETER (LOCDBG = .FALSE.)
!
!
#include <infind.h>
#include <blocks.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <second.h>
!
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
!
      IF (NSYM .NE. 1) THEN
         WRITE(LUPRI,*)'CC3 Triplet excitation energies using the'
         WRITE(LUPRI,*)'noddy code can only be performed for NSYM=1'
         CALL QUIT('No symmetry in this part yet')
      ENDIF
!
      IF (LOCDBG) THEN
        WRITE(LUPRI,*)'ENTERING THE CC3 FILE.'
        WRITE(LUPRI,*)'The energy (ECURR) in this iteration is ',ECURR
        RHO2N = DDOT(NT2SQ(ISYMTR),C2AMP(1),1,C2AMP(1),1)
        WRITE(LUPRI,*)
     *     'The norm of the C2(+) vector in CC3 :',RHO2N
        RHO2N = DDOT(NT2SQ(ISYMTR),C2AMM(1),1,C2AMM(1),1)
        WRITE(LUPRI,*)
     *     'The norm of the C2(-) vector in CC3 :',RHO2N
        RHO2N = DDOT(NT2SQ(1),T2AM(1),1,T2AM(1),1)
        WRITE(LUPRI,*)
     *     'The norm of the T2 vector in CC3 :',RHO2N
        WRITE(LUPRI,*)'  '
        CALL FLSHFO(LUPRI)
      ENDIF
!
!------------------------
!     Dynamic Allocation.
!------------------------
!
      KCMO   = 1
      KSCR1  = KCMO   + NLAMDS
      KFOCKD = KSCR1  + NT1AMX
      KINT1T = KFOCKD + NORBTS
      KINT2T = KINT1T + NT1AMX*NVIRT*NVIRT
      KINT3T = KINT2T + NRHFT*NRHFT*NT1AMX
      KINT4T = KINT3T + NT1AMX*NVIRT*NVIRT
      KINT1S = KINT4T + NRHFT*NRHFT*NT1AMX
      KINT2S = KINT1S + NT1AMX*NVIRT*NVIRT
      KINT3S = KINT2S + NRHFT*NRHFT*NT1AMX
      KINT4S = KINT3S + NT1AMX*NVIRT*NVIRT
      KINT5S = KINT4S + NRHFT*NRHFT*NT1AMX
      KINT6S = KINT5S + NT1AMX*NVIRT*NVIRT
      KT3AM  = KINT6S + NRHFT*NRHFT*NT1AMX
      KOME1  = KT3AM  + NT1AMX*NT1AMX*NT1AMX
      KOME2P = KOME1  + NT1AMX
      KOME2M = KOME2P + NT1AMX*NT1AMX
      KIAJB  = KOME2M + NT1AMX*NT1AMX
      KYIAJB = KIAJB  + NT1AMX*NT1AMX
      KINTTOT= KYIAJB + NT1AMX*NT1AMX
      KINTTOT2= KINTTOT+ NORBT*NORBT*NORBT*NORBT
!
!     This is intermediates needed for the integrals.
!     Not needed at the same time
!
      KINTTMP1= KINTTOT2 + NORBT*NORBT*NORBT*NORBT
      KINTTMP2= KINTTMP1 + NBAST*NBAST*NBAST*NBAST
!
      KTRTMP1= KINTTOT2 + NORBT*NORBT*NORBT*NORBT
      KTRTMP2= KTRTMP1 + NRHFT*NORBT*NORBT*NORBT
      KTRTMP3= KTRTMP2 + NVIRT*NORBT*NORBT*NORBT
      KTRTMP4= KTRTMP3 + NRHFT*NORBT*NORBT*NORBT
!
      IF (NBAST*NBAST*NBAST*NBAST .GT.
     *    (NRHFT+NVIRT)*NORBT*NORBT*NORBT) THEN
!
          KEND1 = KINTTMP2 + NBAST*NBAST*NBAST*NBAST
      ELSE
          KEND1 = KTRTMP4 + NVIRT*NORBT*NORBT*NORBT
      ENDIF
!
      LWRK1  = LWORK  - KEND1
!
!------------------------------------------------------
!     We sum the minus up in KOME2,
!     which is put at the same place as KIAJB.
!------------------------------------------------------
!
      KOME2  = KIAJB
!
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CC3_TRIPLET')
      ENDIF
!
!--------------------------------
!     Initialize integral arrays.
!--------------------------------
!
      CALL DZERO(WORK(KINT1T),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2T),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KINT3T),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT4T),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KINT1S),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2S),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KINT3S),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT4S),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KINT5S),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT6S),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KT3AM),NT1AMX*NT1AMX*NT1AMX)
      CALL DZERO(WORK(KOME1),NT1AMX)
      CALL DZERO(WORK(KOME2P),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KOME2M),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KIAJB),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KYIAJB),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KINTTOT),NORBT*NORBT*NORBT*NORBT)
      CALL DZERO(WORK(KINTTOT2),NORBT*NORBT*NORBT*NORBT)
!
      IF (NBAST*NBAST*NBAST*NBAST .GT.
     *    (NRHFT+NVIRT)*NORBT*NORBT*NORBT) THEN
!
         CALL DZERO(WORK(KINTTMP1),NBAST*NBAST*NBAST*NBAST)
         CALL DZERO(WORK(KINTTMP2),NBAST*NBAST*NBAST*NBAST)
      ELSE
         CALL DZERO(WORK(KTRTMP1),NRHFT*NORBT*NORBT*NORBT)
         CALL DZERO(WORK(KTRTMP2),NVIRT*NORBT*NORBT*NORBT)
         CALL DZERO(WORK(KTRTMP3),NRHFT*NORBT*NORBT*NORBT)
         CALL DZERO(WORK(KTRTMP4),NVIRT*NORBT*NORBT*NORBT)
      ENDIF
!
!----------------------------------------------
!     Read caconical orbital energies &
!     MO-coefficients from interface file.
!----------------------------------------------
!
      IF (LUSIFC .LE. 0) THEN
        LUSIFC = -1
        CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     *              IDUMMY,.FALSE.)
      ENDIF
!
      REWIND (LUSIFC)
!
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
!
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBTS)
      READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDS)
!
      CALL GPCLOSE(LUSIFC,'KEEP')
!
!===============================================
!     Frozen orbital stuff.
!===============================================
!
       CALL CMO_REORDER(WORK(KCMO),WORK(KEND1),LWRK1)
!
      IF (FROIMP .OR. FROEXP)
     *   CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND1),LWRK1)
!
!====================================================
!     Start the loop over distributions of integrals.
!====================================================
!
      IF (DIRECT) THEN
         NTOSYM = 1
!        DTIME  = SECOND()
!!!!        CALL HERDI1(WORK(KEND1),LWRK1,IPRINT)
         CALL QUIT('Direct not implemented')
!        DTIME  = SECOND() - DTIME
!        TIMHER1 = TIMHER1 + DTIME
      ELSE
         NTOSYM = NSYM
      ENDIF
!
      DO 100 ISYMD1 = 1,NTOSYM
!
         IF (DIRECT) THEN
            NTOT = MAXSHL
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
!
         DO 110 ILLL = 1,NTOT
!
!-----------------------------------------------------------------
!           If direct calculate the integrals and transposed t2am.
!-----------------------------------------------------------------
!
            IF (DIRECT) THEN
!              DTIME  = SECOND()
!!!!              CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,IPRINT)
               CALL QUIT('Direct not implemented')
!              DTIME  = SECOND() - DTIME
!              TIMHER2 = TIMHER2 + DTIME
!
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK - KEND1
               IF (LWRK1 .LT. 0) THEN
                 CALL QUIT('Insufficient space in cc3_trip.F (KRECNR)')
               ENDIF
            ELSE
               NUMDIS = 1
            ENDIF
!
!-----------------------------------------------------
!           Loop over number of distributions in disk.
!-----------------------------------------------------
!
            DO 120 IDEL2 = 1,NUMDIS
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
!
!------------------------------------------
!              Work space allocation no. 2.
!------------------------------------------
!
               KXINT  = KEND1
               KEND2  = KXINT + NDISAO(ISYDIS)
               LWRK2  = LWORK - KEND2
!
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,
     *                           'Available : ',LWORK
                  CALL QUIT('Insufficient space in CC3_TRIPLET')
               ENDIF
!
!-----------------------------------------
!              Read in batch of integrals.
!-----------------------------------------
!
               DTIME   = SECOND()
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
               DTIME   = SECOND() - DTIME
               TIMRDAO = TIMRDAO  + DTIME
!
!----------------------------------------------------
!              Calculate the integrals that are
!              needed in CC3. INTTOT contains T_{1}
!              while INTTOT2 doesnt.
!----------------------------------------------------
!
               CALL CC3EXCI_INT1(WORK(KCMO),WORK(KXINT),
     *                           IDEL,WORK(KINTTOT2),
     *                           WORK(KINTTMP1),WORK(KINTTMP2))
!
               CALL CC3EXCI_INT2(WORK(KINTTOT),XLAMDP,XLAMDH,
     *                           WORK(KXINT),IDEL,
     *                           WORK(KINTTMP1),WORK(KINTTMP2))
!
!
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
!
       IF (LOCDBG) THEN
!
          DO I = 1, NORBT*NORBT*NORBT*NORBT
             IF (ABS(WORK(KINTTOT - 1 + I)) .GT. 1.0D-8) THEN
               WRITE(LUPRI,*)'Hat{g}(I) ',WORK(KINTTOT - 1 + I)
             ENDIF
          ENDDO
!
             WRITE(LUPRI,*)'       '
!
          DO I = 1, NORBT*NORBT*NORBT*NORBT
             IF (ABS(WORK(KINTTOT2 - 1 + I)) .GT. 1.0D-8) THEN
               WRITE(LUPRI,*)'g(I) ',WORK(KINTTOT2 - 1 + I)
             ENDIF
          ENDDO
       ENDIF
!
!
      IF (LOCDBG) THEN
         WRITE(LUPRI,*)'After the integral section'
      ENDIF
!
!-------------------------------------------------------
!     Transform the integrals from g(p,q,r,s) to
!     whatever is needed
!-------------------------------------------------------
!
      CALL CC33_TRANSTEMP(WORK(KINTTOT),C1AM,WORK(KINT1S),
     *                    WORK(KINT2S),WORK(KINT3S),
     *                    WORK(KINT4S),WORK(KINT5S),
     *                    WORK(KINT6S),WORK(KINTTOT2),
     *                    WORK(KIAJB),WORK(KYIAJB),
     *                    WORK(KINT1T),WORK(KINT2T),
     *                    WORK(KINT3T),WORK(KINT4T),
     *                    WORK(KTRTMP1),WORK(KTRTMP2),
     *                    WORK(KTRTMP3),WORK(KTRTMP4))
!
!
!---------------------------------------
!     Calculate the triple amplitudes.
!---------------------------------------
!
      IF (LOCDBG) THEN
         WRITE(LUPRI,*)'Before _T3AM'
         CALL FLSHFO(LUPRI)
      ENDIF
!
      CALL CCSDT_T3AM(WORK(KT3AM),WORK(KINT1S),WORK(KINT2S),T2AM,
     *                WORK(KSCR1),WORK(KFOCKD))
!
!---------------------------------------------------------
!     Calculate the corrections from the triple ampl..
!---------------------------------------------------------
!
      IF (LOCDBG) THEN
         WRITE(LUPRI,*)'Before OMEGA2T3AM'
         CALL FLSHFO(LUPRI)
      ENDIF
!
      CALL CC33_OMEGA2T3(WORK(KOME2P),WORK(KOME2M),WORK(KINT3T),
     *                   WORK(KINT4T),WORK(KT3AM))
!
!
      IF (LOCDBG) THEN
         WRITE(LUPRI,*)'After corr. with t3; norm of 2 omegas :'
         TEMP = NT1AMX*NT1AMX
         RHO2N = DDOT(TEMP,WORK(KOME2P),1,WORK(KOME2P),1)
         WRITE(LUPRI,*)'The (+) :',RHO2N
         RHO2N = DDOT(TEMP,WORK(KOME2M),1,WORK(KOME2M),1)
         WRITE(LUPRI,*)'The (-) :',RHO2N
      ENDIF
!
!--------------------------------------------------------------
!     Reset the array and then
!     calculate the triplet triple trialvector amplitudes.
!--------------------------------------------------------------
!
      IF (LOCDBG) THEN
         WRITE(LUPRI,*)'Before _C3AM'
         CALL FLSHFO(6)
      ENDIF
!
      CALL DZERO(WORK(KT3AM),NT1AMX*NT1AMX*NT1AMX)
!
      CALL CC33_C3AM(WORK(KT3AM),WORK(KINT1S),WORK(KINT2S),
     *               WORK(KINT3S),WORK(KINT4S),WORK(KINT5S),
     *               WORK(KINT6S),T2AM,C2AMP,C2AMM,WORK(KSCR1),
     *               WORK(KFOCKD))
!
      IF (LOCDBG) THEN
         TEMP = NT1AMX*NT1AMX*NT1AMX
         RHO2N = DDOT(TEMP,WORK(KT3AM),1,WORK(KT3AM),1)
         WRITE(LUPRI,*)'The norm of the C3 vector is:',RHO2N
      ENDIF
!
!-------------------------------------------------------------
!     Calculate the corrections from the triple trialvector
!     to the 1 excited result vector.
!-------------------------------------------------------------
!
      IF (LOCDBG) THEN
         RHO2N = DDOT(NT1AMX,WORK(KOME1),1,WORK(KOME1),1)
         WRITE(LUPRI,*)'The norm of omega1 at the beginning.:',RHO2N
         WRITE(LUPRI,*)'Before OMEGA1C3AM'
         CALL FLSHFO(LUPRI)
      ENDIF
!
      CALL CC33_OMEGA1C3(WORK(KOME1),WORK(KIAJB),
     *                   WORK(KYIAJB),WORK(KT3AM))
!
      IF (LOCDBG) THEN
         RHO2N = DDOT(NT1AMX,WORK(KOME1),1,WORK(KOME1),1)
         WRITE(LUPRI,*)'The norm of omega1 after the corr. :',RHO2N
      ENDIF
!
!------------------------------------------------------------
!     Sum the T3-results into the CCSD results and go
!     on to calculate the C3 results.
!------------------------------------------------------------
!
!
      CALL DZERO(WORK(KOME2),NT1AMX*NT1AMX)
!
      CALL CC_RVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *              NT2AMA(ISYMTR),IV+ITR-1,NT2AM(ISYMTR),
     *              WORK(KOME2))
!
      IF (LOCDBG) THEN
         TEMP = NT2AM(ISYMTR)
         RHO2N = DDOT(TEMP,WORK(KOME2),1,WORK(KOME2),1)
         WRITE(LUPRI,*)'Norm of (-)-result from CCSD',RHO2N
         RHO2N = DDOT(TEMP,OMEGA2P(1),1,OMEGA2P(1),1)
         WRITE(LUPRI,*)'Norm of (+)-result from CCSD',RHO2N
         WRITE(LUPRI,*)'The Omega2 after the readin of Rho2-.'
         TEMP = NT1AMX*NT1AMX
         RHO2N = DDOT(TEMP,WORK(KOME2P),1,WORK(KOME2P),1)
         WRITE(LUPRI,*)'The norm of the (+)-vector:',RHO2N
         RHO2N = DDOT(TEMP,WORK(KOME2M),1,WORK(KOME2M),1)
         WRITE(LUPRI,*)'The norm of the (-)-vector:',RHO2N
      ENDIF
!
!------------------------------------
!     Sum Omega2 (+ & -) for T3.
!------------------------------------
!
      DO 310 I = 1,NT1AMX
         DO 320 J = 1,I
            IJ = NT1AMX*(I-1) + J
            NIJ = INDEX(I,J)
!
            OMEGA2P(NIJ) = OMEGA2P(NIJ) + WORK(KOME2P+IJ-1)
!
            WORK(KOME2-1+NIJ) = WORK(KOME2-1+NIJ) + WORK(KOME2M+IJ-1)
!
  320    CONTINUE
  310 CONTINUE 
!
      IF (LOCDBG) THEN
         TEMP = NT2AM(ISYMTR)
         RHO2N = DDOT(TEMP,OMEGA2P,1,OMEGA2P,1)
         WRITE(LUPRI,*)
     *        'The norm of the (+) vector after CCSD & T3',RHO2N
         TEMP = NT2AMA(ISYMTR)
         RHO2N = DDOT(TEMP,WORK(KOME2),1,WORK(KOME2),1)
         WRITE(LUPRI,*)
     *        'The norm of the (-) vector after CCSD & T3',RHO2N
      ENDIF
!       
!------------------------------------------------------------
!     Proceed with the calculation omega2. Calculate the
!     C3 dependent terms, but first zero the results vectors.
!------------------------------------------------------------
!
      IF (LOCDBG) THEN
         WRITE(LUPRI,*)'Before OMEGA2C3AM'
         CALL FLSHFO(6)
      ENDIF
!
      CALL DZERO(WORK(KOME2P),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KOME2M),NT1AMX*NT1AMX)
!
      CALL CC33_OMEGA2C3(WORK(KOME2P),WORK(KOME2M),WORK(KINT1T),
     *                   WORK(KINT2T),WORK(KT3AM),FOCK)
!
       IF (LOCDBG) THEN
         WRITE(LUPRI,*)'The Omega2 after the C3 corr.'
         TEMP = NT1AMX*NT1AMX
         RHO2N = DDOT(TEMP,WORK(KOME2P),1,WORK(KOME2P),1)
         WRITE(LUPRI,*)'The norm of the (+)-vector:',RHO2N
         RHO2N = DDOT(TEMP,WORK(KOME2M),1,WORK(KOME2M),1)
         WRITE(LUPRI,*)'The norm of the (-)-vector:',RHO2N
       ENDIF
!
!--------------------------------
!     Sum Omega1
!--------------------------------
!
      DO 400 I = 1,NT1AMX
         OMEGA1(I) = OMEGA1(I) + WORK(KOME1+I-1)
  400 CONTINUE
!
      IF (LOCDBG) THEN
         RHO2N = DDOT(NT1AMX,WORK(KOME1),1,WORK(KOME1),1)
         WRITE(LUPRI,*)
     *        'The norm of omega1 after the summation. :',RHO2N
      ENDIF
!
!------------------------------------------------------
!     Sum Omega2 (+ & -) for the C3 dependent terms
!------------------------------------------------------
!
      IF (LOCDBG) THEN
         WRITE(LUPRI,*)'Before the final summation'
         WRITE(LUPRI,*)' '
         CALL FLSHFO(LUPRI)
      ENDIF
!
      DO 410 I = 1,NT1AMX
         DO 420 J = 1,I
            IJ = NT1AMX*(I-1) + J
            NIJ = INDEX(I,J)
!
            OMEGA2P(NIJ) = OMEGA2P(NIJ) + WORK(KOME2P+IJ-1)
!
            WORK(KOME2-1+NIJ) = WORK(KOME2-1+NIJ) + WORK(KOME2M+IJ-1)
!
  420    CONTINUE
  410 CONTINUE
!
      IF (LOCDBG) THEN
         TEMP = NT2AM(ISYMTR)
         RHO2N = DDOT(TEMP,OMEGA2P,1,OMEGA2P,1)
         WRITE(LUPRI,*)
     *        'The norm of the (+) vector at the end',RHO2N
         TEMP = NT2AMA(ISYMTR)
         RHO2N = DDOT(TEMP,WORK(KOME2),1,WORK(KOME2),1)
         WRITE(LUPRI,*)
     *        'The norm of the (-) vector at the end',RHO2N
         WRITE(LUPRI,*)'End of triple part of the CC3 calc.'
      ENDIF
!
!-------------------------------
!     Write Omega2- to file.
!-------------------------------
!
      CALL CC_WVEC3(LUFR2,FRHO2,NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *              NT2AMA(ISYMTR),IV + ITR -1,NT2AM(ISYMTR),
     *              WORK(KOME2))
!
!--------------------------
!     End the calc.
!--------------------------
!
      RETURN
      END
C  /* Deck cc3_c3am */
      SUBROUTINE CC33_C3AM(ECURR,C3AM,
     *                     XINT1S,XINT2S,XINT3S,XINT4S,XINT5S,
     *                     XINT6S,T2AM,C2AMP,C2AMM,SCR1,FOCKD)
!
!     Written by Kasper Hald Jan 2000.
!
      IMPLICIT NONE
!
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
!
      INTEGER NAI, NCK, NCJ, NBK, NBJ, NDI, NDJ, NDK
      INTEGER NBL, NAL, NCL, NAK, NAJ, NCI, NBI
!
#if defined (SYS_CRAY)
      REAL XINT1S(NT1AMX,NVIRT,NVIRT), XINT2S(NT1AMX,NRHFT,NRHFT)
      REAL XINT3S(NT1AMX,NVIRT,NVIRT), XINT4S(NT1AMX,NRHFT,NRHFT)
      REAL XINT5S(NT1AMX,NVIRT,NVIRT), XINT6S(NT1AMX,NRHFT,NRHFT)
      REAL C3AM(NT1AMX,NT1AMX,NT1AMX),SCR1(NT1AMX),FOCKD(*)
      REAL T2AM(NT1AMX,NT1AMX), HALF, QUART, ONE, ECURR
      REAL AIBJCK, C2AMP(NT1AMX,NT1AMX), C2AMM(NT1AMX,NT1AMX)
#else
      DOUBLE PRECISION XINT1S(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT2S(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION XINT3S(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT4S(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION XINT5S(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT6S(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION C3AM(NT1AMX,NT1AMX,NT1AMX),SCR1(NT1AMX),FOCKD(*)
      DOUBLE PRECISION T2AM(NT1AMX,NT1AMX), HALF, QUART, AIBJCK
      DOUBLE PRECISION C2AMP(NT1AMX,NT1AMX), C2AMM(NT1AMX,NT1AMX)
      DOUBLE PRECISION ONE, ECURR
#endif
!
      PARAMETER (QUART = 0.25D00, HALF=0.5D00, ONE=1.0D00)
!
      DO I = 1,NRHFT
         DO A = 1,NVIRT
            NAI = NVIRT*(I-1) + A
            SCR1(NAI) = FOCKD(NRHFT+A) - FOCKD(I)
         ENDDO
      ENDDO
!
!
      DO 100 K = 1, NRHFT
      DO 110 C = 1, NVIRT
         NCK = NVIRT*(K-1) + C
!
         DO 120 J = 1,NRHFT
!
         NCJ = NVIRT*(J-1) + C
         DO 130 B = 1,NVIRT
!
            NBJ = NVIRT*(J-1) + B
            NBK = NVIRT*(K-1) + B
!
            DO 140 I = 1,NRHFT
            NBI = NVIRT*(I-1) + B
            NCI = NVIRT*(I-1) + C
            DO 150 A = 1,NVIRT
               NAI = NVIRT*(I-1) + A
               NAJ = NVIRT*(J-1) + A
               NAK = NVIRT*(K-1) + A
!
               AIBJCK = 0.0D0
               DO 160 D = 1,NVIRT
!
                  NDI = NVIRT*(I-1) + D
                  NDJ = NVIRT*(J-1) + D
                  NDK = NVIRT*(K-1) + D
!
                  AIBJCK = AIBJCK + HALF*XINT1S(NAI,C,D)*
     *                        C2AMP(NDK,NBJ)
!
                  AIBJCK = AIBJCK + ONE*XINT1S(NBJ,A,D)*
     *                        C2AMM(NCK,NDI)
!
                  AIBJCK = AIBJCK + ONE*XINT1S(NCJ,B,D)*
     *                        C2AMM(NAI,NDK)
!
                  AIBJCK = AIBJCK + HALF*XINT3S(NCJ,B,D)*T2AM(NAI,NDK)
!
                  AIBJCK = AIBJCK + HALF*XINT5S(NAI,C,D)*T2AM(NBJ,NDK)
!
                  AIBJCK = AIBJCK + HALF*XINT5S(NCJ,A,D)*T2AM(NBK,NDI)
!
  160          CONTINUE
!
               DO 170 L = 1,NRHFT
!
                  NAL = NVIRT*(L-1) + A
                  NBL = NVIRT*(L-1) + B
                  NCL = NVIRT*(L-1) + C
!
                 AIBJCK = AIBJCK - HALF*XINT2S(NAI,L,K)*
     *                        C2AMP(NCL,NBJ)
!
                  AIBJCK = AIBJCK + ONE*XINT2S(NBJ,L,I)*
     *                        C2AMM(NAL,NCK)
!
                  AIBJCK = AIBJCK + ONE*XINT2S(NCJ,L,K)*
     *                        C2AMM(NBL,NAI)
!
                  AIBJCK = AIBJCK - HALF*XINT4S(NBK,L,J)*T2AM(NAI,NCL)
!
                  AIBJCK = AIBJCK - HALF*XINT6S(NAI,L,K)*T2AM(NBJ,NCL)
!
                  AIBJCK = AIBJCK - HALF*XINT6S(NCJ,L,I)*T2AM(NBK,NAL)
!
  170          CONTINUE
!
               AIBJCK = AIBJCK/(ECURR-SCR1(NAI)-SCR1(NBJ)-SCR1(NCK))
!
!
               C3AM(NAI,NBJ,NCK) = C3AM(NAI,NBJ,NCK) + AIBJCK
               C3AM(NAI,NBK,NCJ) = C3AM(NAI,NBK,NCJ) - AIBJCK
               C3AM(NAI,NCK,NBJ) = C3AM(NAI,NCK,NBJ) + AIBJCK
               C3AM(NAI,NCJ,NBK) = C3AM(NAI,NCJ,NBK) - AIBJCK
!
  150       CONTINUE
  140       CONTINUE
  130    CONTINUE
  120    CONTINUE
  110 CONTINUE
  100 CONTINUE
!
!
      RETURN
      END
C  /* Deck cc3_omega1c3 */
      SUBROUTINE CC33_OMEGA1C3(OMEG1,XIAJB,YIAJB,C3AM)
!
!     Written by Kasper Hald Jan. 2000
!
      IMPLICIT NONE
!
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
!
      INTEGER NAI, NAJ, NBJ, NBK, NCK, NCI
      INTEGER VIR1, VIR2, OCC1, OCC2
!
#if defined (SYS_CRAY)
      REAL TWO
      REAL OMEG1(NT1AMX),XIAJB(NT1AMX,NT1AMX)
      REAL YIAJB(NT1AMX,NT1AMX)
      REAL C3AM(NT1AMX,NT1AMX,NT1AMX)
#else
      DOUBLE PRECISION TWO
      DOUBLE PRECISION OMEG1(NT1AMX),XIAJB(NT1AMX,NT1AMX)
      DOUBLE PRECISION YIAJB(NT1AMX,NT1AMX)
      DOUBLE PRECISION C3AM(NT1AMX,NT1AMX,NT1AMX)
#endif
!
      PARAMETER (TWO = 2.0D0)
!
!
      DO 100 I = 1,NRHFT
      DO 110 A = 1,NVIRT
         NAI = NVIRT*(I-1) + A
!
         DO 120 J = 1,NRHFT
         NAJ = NVIRT*(J-1) + A
         DO 130 B = 1,NVIRT
!
            NBJ = NVIRT*(J-1) + B
!
            DO 140 K = 1,NRHFT
            NBK = NVIRT*(K-1) + B
            DO 150 C = 1,NVIRT
               NCK = NVIRT*(K-1) + C
               NCI = NVIRT*(I-1) + C
!
!
               OMEG1(NAI) = OMEG1(NAI) + TWO*XIAJB(NCK,NBJ)*
     *                       C3AM(NBJ,NAI,NCK)
     *                       + TWO*YIAJB(NCK,NBJ)*
     *                       (C3AM(NAJ,NBK,NCI)+
     *                        C3AM(NCI,NBK,NAJ))
!
!
  150       CONTINUE
  140       CONTINUE
  130    CONTINUE
  120    CONTINUE
  110 CONTINUE
  100 CONTINUE
!
      RETURN
      END
C  /* Deck cc33_omega2t3 */
      SUBROUTINE CC33_OMEGA2T3(OMEGA2P,OMEGA2M,XINT3T,XINT4T,T3AM)
!
!     Written by Kasper Hald Jan. 2000.
!
!     Calculate the T3,R1 contribution to omega2+:
!     .5*P(ab,ij) (tilde)P(ij)
!        [sum(dlm) t(ai,dj,bl)*g(lm-bar,md) + 
!         sum(dlm)(t(am,bl,di)+t(ai,dl,bm)-2t(ai,dm,bl))*g(lj-bar,md) +
!         sum(del)(2t(ai,dl,ej)-t(al,di,ej)-t(ai,el,dj))*g(b-bar e,ld)
!        ]
!     N.B. Do not calculate with (Tilde)P(ij).
!
!     Calculate the T3,R1 contribution to omega2-:
!     .5*(tilde)P(ab,ij) 
!        [sum(dlm) t(ai,dj,bl)*g(lm-bar,md) + 
!         sum(dlm)(t(am,bl,di)+t(ai,dl,bm)-2t(ai,dm,bl))*g(lj-bar,md) +
!         sum(del)(2t(ai,dl,ej)-t(al,di,ej)-t(ai,el,dj))*g(b-bar e,ld)
!        ]
!
      IMPLICIT NONE
!
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
!
      INTEGER NBK, NAK, NAI, NBJ, NCK, NCJ, NCI, NDJ
      INTEGER NDK, NDI, NBL, NCL, NAL
      INTEGER NAJ, NBI
!
#if defined (SYS_CRAY)
      REAL XINT3T(NT1AMX,NVIRT,NVIRT)
      REAL XINT4T(NT1AMX,NRHFT,NRHFT)
      REAL OMEGA2P(NT1AMX,NT1AMX)
      REAL OMEGA2M(NT1AMX,NT1AMX)
      REAL T3AM(NT1AMX,NT1AMX,NT1AMX)
      REAL HALF,TWO, XAIBJ, ZERO
#else
      DOUBLE PRECISION XINT3T(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT4T(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION OMEGA2P(NT1AMX,NT1AMX)
      DOUBLE PRECISION OMEGA2M(NT1AMX,NT1AMX)
      DOUBLE PRECISION T3AM(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION HALF,TWO, XAIBJ, ZERO, XHALF
#endif
!
      PARAMETER (ZERO = 0.0D00, HALF = 0.5D00, TWO = 2.0D0)
      PARAMETER (XHALF=-0.5D0)
!
!
      DO A = 1, NVIRT
      DO I = 1, NRHFT
         NAI = NVIRT*(I-1) + A
!
         DO J = 1, NRHFT
         NAJ = NVIRT*(J-1) + A
         DO K = 1, NRHFT
            NAK = NVIRT*(K-1) + A
            T3AM(NAI,NAJ,NAK) = 0.0D0
            T3AM(NAI,NAK,NAJ) = 0.0D0
            T3AM(NAJ,NAI,NAK) = 0.0D0
            T3AM(NAJ,NAK,NAI) = 0.0D0
            T3AM(NAK,NAJ,NAI) = 0.0D0
            T3AM(NAK,NAI,NAJ) = 0.0D0
         ENDDO
         ENDDO
!
         DO B = 1, NVIRT
         NBI = NVIRT*(I-1)+B
         DO C = 1, NVIRT
            NCI = NVIRT*(I-1)+C
!
            T3AM(NAI,NBI,NCI) = 0.0D0
            T3AM(NAI,NCI,NBI) = 0.0D0
            T3AM(NBI,NAI,NCI) = 0.0D0
            T3AM(NBI,NCI,NAI) = 0.0D0
            T3AM(NCI,NAI,NBI) = 0.0D0
            T3AM(NCI,NBI,NAI) = 0.0D0
         ENDDO
         ENDDO
      ENDDO
      ENDDO
!
!
      DO 100 I = 1,NRHFT
         DO 110 A = 1,NVIRT
            NAI = NVIRT*(I-1) + A
!
            DO 120 J = 1,NRHFT
               DO 130 B = 1,NVIRT
                  NBJ = NVIRT*(J-1) + B
!
                  IF (NAI .NE. NBJ) THEN
!
                  DO 140 K = 1,NRHFT
                     NBK = NVIRT*(K-1) + B
                     NAK = NVIRT*(K-1) + A
                     DO 150 C = 1,NVIRT
!
                        NCK = NVIRT*(K-1) + C
                        NCJ = NVIRT*(J-1) + C
                        NCI = NVIRT*(I-1) + C
!
                        DO 160 D = 1,NVIRT
!
                           NDJ = NVIRT*(J-1) + D
                           NDK = NVIRT*(K-1) + D
                           NDI = NVIRT*(I-1) + D
!
                           XAIBJ = (TWO*T3AM(NAI,NDK,NCJ) - 
     *                              T3AM(NAK,NDI,NCJ) - 
     *                              T3AM(NAI,NCK,NDJ))*XINT3T(NDK,B,C)
!
                           OMEGA2M(NAI,NBJ) = OMEGA2M(NAI,NBJ) + XAIBJ
!
                           IF ((A .NE. B) .AND. (I .NE. J)) THEN
                           OMEGA2P(NAI,NBJ) = OMEGA2P(NAI,NBJ) + XAIBJ
                           ENDIF
!
  160                   CONTINUE
!
                        DO 170 L = 1,NRHFT
!
                           NBL = NVIRT*(L-1) + B
                           NCL = NVIRT*(L-1) + C
                           NAL = NVIRT*(L-1) + A
!
                           XAIBJ = T3AM(NAI,NCJ,NBL)*XINT4T(NCK,L,K)
     *                           +(T3AM(NAK,NBL,NCI) +
     *                             T3AM(NAI,NCL,NBK) -
     *                           TWO*T3AM(NAI,NCK,NBL))*XINT4T(NCK,L,J)
!
                           OMEGA2M(NAI,NBJ) = OMEGA2M(NAI,NBJ) + XAIBJ
!
                           IF ((A .NE. B) .AND. (I .NE. J)) THEN
                           OMEGA2P(NAI,NBJ) = OMEGA2P(NAI,NBJ) + XAIBJ
                           ENDIF
!
  170                   CONTINUE
!
  150                CONTINUE
  140             CONTINUE
!
                  ENDIF
!
  130          CONTINUE
  120       CONTINUE
!
  110    CONTINUE
  100 CONTINUE
!
      DO 200 NAI = 1,NT1AMX
         DO 210 NBJ = 1,NAI
!
!   The scaling of (+) with 1/2 is done in the end
!   of the CCSD routine.
!
            IF (NAI .NE. NBJ) THEN
               XAIBJ = OMEGA2P(NAI,NBJ) + OMEGA2P(NBJ,NAI)
               OMEGA2P(NAI,NBJ) = XAIBJ
               OMEGA2P(NBJ,NAI) = XAIBJ
               XAIBJ = HALF*(OMEGA2M(NAI,NBJ) - OMEGA2M(NBJ,NAI))
               OMEGA2M(NBJ,NAI) = XAIBJ
               OMEGA2M(NAI,NBJ) = -XAIBJ
            ELSE
               OMEGA2M(NAI,NBJ) = ZERO
            ENDIF
!
  210    CONTINUE
  200 CONTINUE
!
      RETURN
      END
C  /* Deck cc33_omega2c3 */
      SUBROUTINE CC33_OMEGA2C3(OMEGA2P,OMEGA2M,XINT1T,XINT2T,C3AM,FOCK)
!
!     Written by Kasper Hald Jan. 2000.
!
!     Calc. the cont. : 
!
!     (+): .5*P(ab,ij)[sum(dl)( c3(dl,bj,ai)+2*c3(ai,bj,dl)+
!                               c3(di,bl,aj)+c3(bl,di,aj) )*F(ld)
!                     +sum(del)( 2*c3(dl,ei,bj)-c3(el,di,bj)+
!                                2*c3(bj,ei,dl)+2*c3(ei,bj,dl)+
!                                2*c3(dj,bi,el)-c3(bl,dj,ei) )*g(ld,ae)
!                     +sum(dlm)( 2*c3(dl,am,bi)-c3(dm,al,bi)+
!                                2*c3(bi,am,dl)+2*c3(am,bi,dl) +
!                                2*c3(al,di,bm)-c3(di,bl,am) )*g(ld,mj)
!                     ]
!
!     (-): (Tilde)P(ab,ij)[sum(ld) c3((ai,bj,dl)*F(ld)+
!                          sum(del) (c3(bj,di,el)+c3(ei,bj,dl))*g(ae,ld)
!                          sum(dlm) (c3(ai,dm,bl)+c3(bm,ai,dl))*g(ld,mj)
!                         ]
!
      IMPLICIT NONE
!
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
!
      INTEGER NAI,NAJ,NBI,NBJ,NBK,NCK,NCI,NDJ,NDK,NDI
      INTEGER NBL, NAL, NCL, NAK
!
#if defined (SYS_CRAY)
      REAL XINT1T(NT1AMX,NVIRT,NVIRT)
      REAL XINT2T(NT1AMX,NRHFT,NRHFT)
      REAL OMEGA2P(NT1AMX,NT1AMX), HALF, ZERO
      REAL OMEGA2M(NT1AMX,NT1AMX), TWO, XAIBJ
      REAL C3AM(NT1AMX,NT1AMX,NT1AMX),FOCK(NORBT,NORBT)
#else
      DOUBLE PRECISION XINT1T(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT2T(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION OMEGA2P(NT1AMX,NT1AMX), HALF, ZERO
      DOUBLE PRECISION OMEGA2M(NT1AMX,NT1AMX), TWO, XAIBJ
      DOUBLE PRECISION C3AM(NT1AMX,NT1AMX,NT1AMX),FOCK(NORBT,NORBT)
#endif
!
      PARAMETER (ZERO = 0.0D00, HALF = 0.5D00, TWO = 2.0D0)
!
!
      DO 100 I = 1,NRHFT
         DO 110 A = 1,NVIRT
            NAI = NVIRT*(I-1) + A
!
            DO 120 J = 1,NRHFT
               NAJ = NVIRT*(J-1) + A
               DO 130 B = 1,NVIRT
                  NBI = NVIRT*(I-1) + B
                  NBJ = NVIRT*(J-1) + B
!
                  IF (NAI .NE. NBJ) THEN
!
                  DO 140 K = 1,NRHFT
                     NAK = NVIRT*(K-1) + A
                     NBK = NVIRT*(K-1) + B
                     DO 150 C = 1,NVIRT
!
                        NCK = NVIRT*(K-1) + C
                        NCI = NVIRT*(I-1) + C
!
!===================================
!     The plus part of omega2
!===================================
!
                     IF ((A .NE. B) .AND. (I .NE. J)) THEN
!
                        OMEGA2P(NAI,NBJ) = OMEGA2P(NAI,NBJ) +
     *      (C3AM(NCK,NBJ,NAI)+TWO*C3AM(NAI,NBJ,NCK)+
     *       C3AM(NCI,NBK,NAJ)+C3AM(NBK,NCI,NAJ))*FOCK(K,NRHFT+C)
!
                     ENDIF
!
!===================================
!     The minus part of omega2
!===================================
!
                        OMEGA2M(NAI,NBJ) = OMEGA2M(NAI,NBJ) +
     *       C3AM(NAI,NBJ,NCK)*FOCK(K,NRHFT+C)
!
!-----------------------------------
!-----------------------------------
!
                        DO 160 D = 1,NVIRT
!
                           NDJ = NVIRT*(J-1) + D
                           NDK = NVIRT*(K-1) + D
                           NDI = NVIRT*(I-1) + D
!
!===================================
!     The plus part of omega2
!===================================
!
                     IF ((A .NE. B) .AND. (I .NE. J)) THEN
!
                           OMEGA2P(NAI,NBJ) = OMEGA2P(NAI,NBJ) +
     *  (TWO*C3AM(NDK,NCI,NBJ)-C3AM(NCK,NDI,NBJ)+TWO*C3AM(NBJ,NCI,NDK)
     *  +TWO*C3AM(NCI,NBJ,NDK)+TWO*C3AM(NDJ,NBI,NCK)-C3AM(NBK,NDJ,NCI))
     *  *XINT1T(NDK,A,C)
!
                     ENDIF
!
!===================================
!     The minus part of omega2
!===================================
!
                           OMEGA2M(NAI,NBJ) = OMEGA2M(NAI,NBJ) +
     *  (C3AM(NBJ,NDI,NCK)+C3AM(NCI,NBJ,NDK))*XINT1T(NDK,A,C)
!
  160                   CONTINUE
!
!-----------------------------------
!-----------------------------------
!
                        DO 170 L = 1,NRHFT
!
                           NBL = NVIRT*(L-1) + B
                           NCL = NVIRT*(L-1) + C
                           NAL = NVIRT*(L-1) + A
!
!===================================
!     The plus part of omega2
!===================================
!
                     IF ((A .NE. B) .AND. (I .NE. J)) THEN
!
                           OMEGA2P(NAI,NBJ) = OMEGA2P(NAI,NBJ) +
     *  (TWO*C3AM(NCL,NAK,NBI)-C3AM(NCK,NAL,NBI)+TWO*C3AM(NBI,NAK,NCL)
     *  +TWO*C3AM(NAK,NBI,NCL)+TWO*C3AM(NAL,NCI,NBK)-C3AM(NCI,NBL,NAK))
     *  *XINT2T(NCL,K,J)
!
                     ENDIF
!
!===================================
!     The minus part of omega2
!===================================
!
                           OMEGA2M(NAI,NBJ) = OMEGA2M(NAI,NBJ) +
     *  (C3AM(NAI,NCK,NBL)+C3AM(NBK,NAI,NCL))*XINT2T(NCL,K,J)
!
  170                   CONTINUE
!
  150                CONTINUE
  140             CONTINUE
!
                  ENDIF
!
  130          CONTINUE
  120       CONTINUE
!
  110    CONTINUE
  100 CONTINUE
!
      DO 200 NAI = 1,NT1AMX
         DO 210 NBJ = 1,NAI
!
!     (Anti)symmetrisation
!
!     The Half for + is done in the CCSD part.
!
            IF (NAI .NE. NBJ) THEN
               XAIBJ = OMEGA2P(NAI,NBJ)+OMEGA2P(NBJ,NAI)
               OMEGA2P(NAI,NBJ) = XAIBJ
               OMEGA2P(NBJ,NAI) = XAIBJ
               XAIBJ = OMEGA2M(NAI,NBJ) - OMEGA2M(NBJ,NAI)
               OMEGA2M(NAI,NBJ) = -XAIBJ
               OMEGA2M(NBJ,NAI) = XAIBJ
            ELSE
               OMEGA2P(NAI,NBJ) = ZERO
               OMEGA2M(NAI,NBJ) = ZERO
            ENDIF
!
  210    CONTINUE
  200 CONTINUE
!
      RETURN
      END
C  /* Deck ccsdt_t3am */
      SUBROUTINE CCSDT_T3AM(T3AM,XINT1,XINT2,T2AM,SCR1,FOCKD)
!
!     Written by Henrik Koch.
!
!     Calculates : T3 
!
!     N.B.
!     This routine is also present in the file ccsd_triple.F
!     which however is not a normal part of the cc-program.
!
      IMPLICIT NONE
!
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
!
      INTEGER NAI, NCK, NBJ, NDJ, NBL
!
#if defined (SYS_CRAY)
      REAL XINT1(NT1AMX,NVIRT,NVIRT), XINT2(NT1AMX,NRHFT,NRHFT)
      REAL T3AM(NT1AMX,NT1AMX,NT1AMX),SCR1(NT1AMX),FOCKD(*)
      REAL T2AM(NT1AMX,NT1AMX), AIBJCK
#else
      DOUBLE PRECISION XINT1(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT2(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION T3AM(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION SCR1(NT1AMX),FOCKD(*)
      DOUBLE PRECISION T2AM(NT1AMX,NT1AMX), AIBJCK
#endif
!
!
      DO 50 I = 1,NRHFT
         DO 60 A = 1,NVIRT
            NAI = NVIRT*(I-1) + A
            SCR1(NAI) = FOCKD(NRHFT+A) - FOCKD(I)
   60    CONTINUE
   50 CONTINUE
!
      DO 100 NCK = 1,NT1AMX
!
         DO 110 J = 1,NRHFT
            DO 120 B = 1,NVIRT
!
               NBJ = NVIRT*(J-1) + B
!
               DO 130 NAI = 1,NT1AMX
!
                  AIBJCK = 0.0D0
                  DO 140 D = 1,NVIRT
!
                     NDJ = NVIRT*(J-1) + D
!
                     AIBJCK = AIBJCK + XINT1(NCK,B,D)*T2AM(NDJ,NAI)
!
  140             CONTINUE
!
                  DO 150 L = 1,NRHFT
!
                     NBL = NVIRT*(L-1) + B
!
                     AIBJCK = AIBJCK - XINT2(NCK,L,J)*T2AM(NBL,NAI)
!
  150             CONTINUE
!
                  AIBJCK = AIBJCK/(SCR1(NAI) + SCR1(NBJ) + SCR1(NCK))
!
                  T3AM(NAI,NBJ,NCK) = T3AM(NAI,NBJ,NCK) - AIBJCK
                  T3AM(NAI,NCK,NBJ) = T3AM(NAI,NCK,NBJ) - AIBJCK
                  T3AM(NBJ,NAI,NCK) = T3AM(NBJ,NAI,NCK) - AIBJCK
                  T3AM(NCK,NAI,NBJ) = T3AM(NCK,NAI,NBJ) - AIBJCK
                  T3AM(NBJ,NCK,NAI) = T3AM(NBJ,NCK,NAI) - AIBJCK
                  T3AM(NCK,NBJ,NAI) = T3AM(NCK,NBJ,NAI) - AIBJCK
!
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
!
      RETURN
      END
C  /* Deck cc3exci_int2 */
      SUBROUTINE CC3EXCI_INT2(XINT,XLAMDP,XLAMDH,AOINT,IDEL,
     *                        XINTTEMP,XINTTEMP2)
!
!     Written by Kasper Hald, April 2000
!
      IMPLICIT NONE
!
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
!
      INTEGER IB, INDEX, IDEL, NAB
!
#if defined (SYS_CRAY)
      REAL XINT(NORBT,NORBT,NORBT,NORBT)
      REAL XLAMDP(NBAST,NORBT), XLAMDH(NBAST,NORBT)
      REAL AOINT((NBAST*(NBAST+1))/2,NBAST), ZERO
      REAL XINTTEMP(NBAST,NBAST,NBAST,NBAST)
      REAL XINTTEMP2(NBAST,NBAST,NBAST,NBAST)
#else
      DOUBLE PRECISION XINT(NORBT,NORBT,NORBT,NORBT)
      DOUBLE PRECISION XLAMDP(NBAST,NORBT), XLAMDH(NBAST,NORBT)
      DOUBLE PRECISION AOINT((NBAST*(NBAST+1))/2,NBAST), ZERO
      DOUBLE PRECISION XINTTEMP(NBAST,NBAST,NBAST,NBAST)
      DOUBLE PRECISION XINTTEMP2(NBAST,NBAST,NBAST,NBAST)
#endif
!
      PARAMETER (ZERO = 0.0D0)
!
!
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
!
      CALL DZERO(XINTTEMP,NBAST*NBAST*NBAST*NBAST)
      CALL DZERO(XINTTEMP2,NBAST*NBAST*NBAST*NBAST)
!
!============================================
!     Calculate the needed integrals.
!============================================
!
      DO G = 1,NBAST
         DO IB = 1,NBAST
            DO A = 1,NBAST
               NAB = INDEX(A,IB)
!
               IF (.NOT. (AOINT(NAB,G) .EQ. ZERO)) THEN
!
               DO P = 1, NORBT
                  XINTTEMP(P,IB,G,IDEL) = XINTTEMP(P,IB,G,IDEL)
     *                + AOINT(NAB,G)*XLAMDP(A,P)
               ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!
      DO P = 1, NORBT
         DO G = 1, NBAST
            DO IB = 1, NBAST
               IF (.NOT. (XINTTEMP(P,IB,G,IDEL) .EQ. ZERO)) THEN
               DO Q = 1, NORBT
                  XINTTEMP2(P,Q,G,IDEL) = XINTTEMP2(P,Q,G,IDEL)
     *                + XINTTEMP(P,IB,G,IDEL)*XLAMDH(IB,Q)
               ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!
      CALL DZERO(XINTTEMP,NBAST*NBAST*NBAST*NBAST)
! 
      DO P = 1, NORBT
         DO Q = 1, NORBT
            DO G = 1, NBAST
               IF (.NOT. (XINTTEMP2(P,Q,G,IDEL) .EQ. ZERO)) THEN
               DO R = 1, NORBT
                  XINTTEMP(P,Q,R,IDEL) = XINTTEMP(P,Q,R,IDEL)
     *                + XINTTEMP2(P,Q,G,IDEL)*XLAMDP(G,R)
               ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO
! 
      DO P = 1, NORBT
         DO Q = 1, NORBT
            DO R = 1, NORBT
               IF (.NOT. (XINTTEMP(P,Q,R,IDEL) .EQ. ZERO)) THEN
               DO S = 1, NORBT
                  XINT(P,Q,R,S) = XINT(P,Q,R,S) +
     *                  XINTTEMP(P,Q,R,IDEL)*XLAMDH(IDEL,S) 
               ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!
      RETURN
      END
C  /* Deck cc33_transtemp */
      SUBROUTINE CC33_TRANSTEMP(XINT,C1AM,XINT1S,XINT2S,XINT3S,
     *                          XINT4S,XINT5S,XINT6S,XINT2,
     *                          XIAJB,YIAJB,XINT1T,XINT2T,
     *                          XINT3T,XINT4T,XITRAN1,XITRAN2,
     *                          XITRAN3,XITRAN4)
!
!     Written by Kasper Hald January 2000.
!
!
      IMPLICIT NONE
!
#include <priunit.h>
#include <maxorb.h>
#include <maxash.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
#include <eribuf.h>
#include <ccsdsym.h>
#include <ccorb.h>
!
      INTEGER NCK, NCL, NDK, NAI, NBJ, NBK
!
#if defined (SYS_CRAY)
      REAL XINT(NORBT,NORBT,NORBT,NORBT)
      REAL XINT2(NORBT,NORBT,NORBT,NORBT)
      REAL XITRAN1(NRHFT,NORBT,NORBT,NORBT)
      REAL XITRAN2(NVIRT,NORBT,NORBT,NORBT)
      REAL XITRAN3(NORBT,NRHFT,NORBT,NORBT)
      REAL XITRAN4(NORBT,NVIRT,NORBT,NORBT)
      REAL XINT1S(NT1AMX,NVIRT,NVIRT)
      REAL XINT2S(NT1AMX,NRHFT,NRHFT)
      REAL XINT3S(NT1AMX,NVIRT,NVIRT)
      REAL XINT4S(NT1AMX,NRHFT,NRHFT)
      REAL XINT5S(NT1AMX,NVIRT,NVIRT)
      REAL XINT6S(NT1AMX,NRHFT,NRHFT)
      REAL XINT1T(NT1AMX,NVIRT,NVIRT)
      REAL XINT2T(NT1AMX,NRHFT,NRHFT)
      REAL XINT3T(NT1AMX,NVIRT,NVIRT)
      REAL XINT4T(NT1AMX,NRHFT,NRHFT)
      REAL TEMP
      REAL C1AM(NT1AMX), TWO
      REAL XIAJB(NT1AMX,NT1AMX), YIAJB(NT1AMX,NT1AMX)
#else
      DOUBLE PRECISION XINT(NORBT,NORBT,NORBT,NORBT)
      DOUBLE PRECISION XINT2(NORBT,NORBT,NORBT,NORBT)
      DOUBLE PRECISION XITRAN1(NRHFT,NORBT,NORBT,NORBT)
      DOUBLE PRECISION XITRAN2(NVIRT,NORBT,NORBT,NORBT)
      DOUBLE PRECISION XITRAN3(NORBT,NRHFT,NORBT,NORBT)
      DOUBLE PRECISION XITRAN4(NORBT,NVIRT,NORBT,NORBT)
      DOUBLE PRECISION XINT1S(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT2S(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION XINT3S(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT4S(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION XINT5S(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT6S(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION XINT1T(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT2T(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION XINT3T(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT4T(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION TEMP
      DOUBLE PRECISION C1AM(NT1AMX), TWO
      DOUBLE PRECISION XIAJB(NT1AMX,NT1AMX), YIAJB(NT1AMX,NT1AMX)
#endif
!
      LOGICAL LOCDBG
!
      PARAMETER (TWO = 2.0D0)
      PARAMETER (LOCDBG = .FALSE.)
!
!
      DO A = 1, NVIRT
         DO I = 1, NRHFT
            DO B = 1, NVIRT
               DO C = 1, NVIRT
                  NAI = NVIRT*(I-1)+A
                  XINT1S(NAI,B,C) = XINT1S(NAI,B,C) +
     *                              XINT(NRHFT+A,I,NRHFT+B,NRHFT+C)
                  XINT1T(NAI,B,C) = XINT1T(NAI,B,C) +
     *                              XINT(I,NRHFT+A,NRHFT+B,NRHFT+C)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
      DO I = 1, NRHFT
         DO A = 1, NVIRT
            DO J = 1, NRHFT
               DO K = 1, NRHFT
                  NAI = NVIRT*(I-1)+A
                  XINT2S(NAI,J,K) = XINT2S(NAI,J,K) +
     *                              XINT(NRHFT+A,I,J,K)
                  XINT2T(NAI,J,K) = XINT2T(NAI,J,K) +
     *                              XINT(I,NRHFT+A,J,K)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
      DO P = 1, NORBT
!
         DO Q = 1, NORBT
!
            DO R = 1, NORBT
!
               DO K = 1, NRHFT
!
                  TEMP = 0.0D0
                  DO D = 1, NVIRT
!
                     NDK = NVIRT*(K-1) + D
                     TEMP = TEMP + C1AM(NDK)*XINT(P,D+NRHFT,Q,R)
!
                  ENDDO ! The D-loop stops here.
!
                     XITRAN3(P,K,Q,R) = TEMP
!
               ENDDO ! The K-loop stops here.
!
               DO C = 1, NVIRT
!
                  TEMP = 0.0D0
                  DO L = 1, NRHFT
!
                     NCL = NVIRT*(L-1) + C
                     TEMP = TEMP + C1AM(NCL)*XINT(L,P,Q,R)
!
                  ENDDO ! The L-loop stops here
!
                  XITRAN2(C,P,Q,R) = -TEMP
!
               ENDDO ! The C-loop ends here.
!
            ENDDO ! The Q-loop stops here.
         ENDDO ! The R-loop stops here.
      ENDDO ! The S-loop stops here.
!
!====================================
!     Debug information.
!====================================
!
      IF (LOCDBG) THEN
        DO A = 1, NVIRT
           DO K = 1, NORBT
              DO B = 1, NORBT
                 DO C = 1, NORBT
                   IF (ABS(XITRAN2(A,K,B,C)) .GT. 1.0D-8) THEN
                     WRITE(LUPRI,*)'C, P, R, S :',A,K,B,C
                     WRITE(LUPRI,*)'XITRAN2(C,P,R,S) = ',
     *                              XITRAN2(A,K,B,C)
                   ENDIF
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        WRITE(LUPRI,*)'                          '
!
        DO A = 1, NORBT
           DO K = 1, NRHFT
              DO B = 1, NORBT
                 DO C = 1, NORBT
                   IF (ABS(XITRAN3(A,K,B,C)) .GT. 1.0D-8) THEN
                     WRITE(LUPRI,*)'P, K, R, S :',A,K,B,C
                     WRITE(LUPRI,*)'XITRAN3(P,K,R,S) = ',
     *                              XITRAN3(A,K,B,C)
                   ENDIF
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        WRITE(LUPRI,*)'                          '
      ENDIF
!
!=======================
!
      DO A = 1, NVIRT
         DO B = 1, NVIRT
            DO C = 1, NVIRT
               DO K = 1, NRHFT
!
                  NCK = NVIRT*(K-1) + C
                  XINT3S(NCK,A,B) = XINT3S(NCK,A,B) 
     *                            + XITRAN2(C,B+NRHFT,A+NRHFT,K)
     *                            - XITRAN2(C,K,A+NRHFT,B+NRHFT) 
     *                            - XITRAN3(C+NRHFT,K,A+NRHFT,B+NRHFT)
                  XINT5S(NCK,A,B) = XINT5S(NCK,A,B) 
     *                            + XITRAN2(A,B+NRHFT,C+NRHFT,K)
     *                            - XITRAN2(C,K,A+NRHFT,B+NRHFT) 
     *                            - XITRAN3(C+NRHFT,K,A+NRHFT,B+NRHFT)
!
               ENDDO ! The K-loop ends here
            ENDDO ! The C-loop ends here.
         ENDDO ! The B-loop ends here.
      ENDDO ! The A-loop ends here.
!
      DO I = 1, NRHFT
         DO J = 1, NRHFT
            DO C = 1, NVIRT
               DO K = 1, NRHFT
!
                  NCK = NVIRT*(K-1) + C
                  XINT4S(NCK,I,J) = XINT4S(NCK,I,J) 
     *                            + XITRAN3(I,K,C+NRHFT,J)
     *                            - XITRAN3(C+NRHFT,K,I,J) 
     *                            - XITRAN2(C,K,I,J)
!
                  XINT6S(NCK,I,J) = XINT6S(NCK,I,J) 
     *                            + XITRAN3(I,J,C+NRHFT,K)
     *                            - XITRAN3(C+NRHFT,K,I,J) 
     *                            - XITRAN2(C,K,I,J)
!
               ENDDO ! The K-loop ends here
            ENDDO ! The C-loop ends here
         ENDDO ! The J-loop ends here
      ENDDO ! The I-loop ends here
!
!==============================
!     More debug information.
!==============================
!
      IF (LOCDBG) THEN
        DO A = 1, NT1AMX
           DO I = 1, NRHFT
              DO J = 1, NRHFT
                 IF (ABS(XINT3S(A,I,J)) .GT. 1.0D-8) THEN
                   WRITE(LUPRI,*)'XINT3S(NAI,J,K) ',XINT3S(A,I,J)
                   WRITE(LUPRI,*)'NAI, J, K ',A,I,J
                 ENDIF
                 IF (ABS(XINT5S(A,I,J)) .GT. 1.0D-8) THEN
                   WRITE(LUPRI,*)'XINT5S(NAI,J,K) ',XINT5S(A,I,J)
                   WRITE(LUPRI,*)'NAI, J, K ',A,I,J
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
!
        DO A = 1, NT1AMX
           DO C = 1, NVIRT
              DO D = 1, NVIRT
                 IF (ABS(XINT4S(A,C,D)) .GT. 1.0D-8) THEN
                   WRITE(LUPRI,*)'XINT4S(NAI,B,C) ',XINT4S(A,C,D)
                   WRITE(LUPRI,*)'NAI, B, C ',A,C,D
                 ENDIF
                 IF (ABS(XINT6S(A,C,D)) .GT. 1.0D-8) THEN
                   WRITE(LUPRI,*)'XINT6S(NAI,B,C) ',XINT6S(A,C,D)
                   WRITE(LUPRI,*)'NAI, B, C ',A,C,D
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
      ENDIF
!
!===========================
!
      DO A = 1, NVIRT
      DO I = 1, NRHFT
         NAI = NVIRT*(I-1) + A
!
         DO B = 1, NVIRT
         DO J = 1, NRHFT
            NBJ = NVIRT*(J-1) + B
!
            XIAJB(NAI,NBJ) = TWO*XINT2(A+NRHFT,I,B+NRHFT,J)
     *                      - XINT2(A+NRHFT,J,B+NRHFT,I)
            YIAJB(NAI,NBJ) = XINT2(A+NRHFT,I,B+NRHFT,J)
!
         ENDDO
         ENDDO
      ENDDO
      ENDDO
!
!
      DO A = 1, NVIRT
          DO I = 1, NRHFT
             NAI = NVIRT*(I-1) + A
             DO B = 1, NVIRT
                DO C = 1, NVIRT
                   TEMP = 0.0D0
                   DO K = 1, NRHFT
                      NCK  = NVIRT*(K-1) + C
                      TEMP = TEMP + C1AM(NCK)*XINT(I,A+NRHFT,K,B+NRHFT)
                   ENDDO
                   XINT3T(NAI,C,B) = -TEMP
                ENDDO
             ENDDO
          ENDDO
      ENDDO
!
!
      DO A = 1, NVIRT
         DO I = 1, NRHFT
            NAI = NVIRT*(I-1) + A
            DO J = 1, NRHFT
               DO K = 1, NRHFT
                  TEMP = 0.0D0
                  DO B = 1, NVIRT
!
                     NBK = NVIRT*(K-1) + B
                     TEMP = TEMP + C1AM(NBK)*XINT(I,A+NRHFT,J,B+NRHFT)
!
                  ENDDO
!
                  XINT4T(NAI,J,K) =  TEMP
!
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
!
      RETURN
      END
C  /* Deck intcompare */
      SUBROUTINE INTCOMPARE(XIAJB,XIAJB2,XIAJB3,YIAJB,YIAJB2,YIAJB3)
!
!     Written by K. Hald, Jan. 2000.
!     
!
      IMPLICIT NONE
!
#include <priunit.h>
#include <ccsdsym.h>
!
#if defined (SYS_CRAY)
      REAL XIAJB(NT1AMX,NT1AMX), XIAJB2(NT1AMX,NT1AMX)
      REAL YIAJB(NT1AMX,NT1AMX), YIAJB2(NT1AMX,NT1AMX)
      REAL XIAJB3(NT1AMX,NT1AMX), YIAJB3(NT1AMX,NT1AMX)
#else
      DOUBLE PRECISION XIAJB(NT1AMX,NT1AMX), XIAJB2(NT1AMX,NT1AMX)
      DOUBLE PRECISION YIAJB(NT1AMX,NT1AMX), YIAJB2(NT1AMX,NT1AMX)
      DOUBLE PRECISION XIAJB3(NT1AMX,NT1AMX), YIAJB3(NT1AMX,NT1AMX)
#endif
!
      DO A = 1, NT1AMX
         DO I = 1, NT1AMX
            IF (ABS(XIAJB(A,I)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'IAJB         (',A,',',I,') = ',
     *                         XIAJB(A,I)
            ENDIF
            IF (ABS(XIAJB2(A,I)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'IAJBTEMP     (',A,',',I,') = ',
     *                         XIAJB2(A,I)
            ENDIF
            IF (ABS(XIAJB3(A,I)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'IAJBTEMP(CMO)(',A,',',I,') = ',
     *                         XIAJB3(A,I)
            ENDIF
            IF ((ABS(XIAJB(A,I)-XIAJB2(A,I)) .GT. 1.0D-8)) THEN
                WRITE(LUPRI,*)'WARNING THE TWO FIRST VALUES DIFFER!!'
            ENDIF
            IF ((ABS(XIAJB(A,I)-XIAJB3(A,I)) .GT. 1.0D-8)) THEN
                WRITE(LUPRI,*)'WARNING THE 1.ST AND 3.RD VALUE DIFFER!'
            ENDIF
            IF ((ABS(XIAJB2(A,I)-XIAJB3(A,I)) .GT. 1.0D-8)) THEN
                WRITE(LUPRI,*)'WARNING THE 2.ND AND 3.RD VALUE DIFFER!'
            ENDIF
         ENDDO
      ENDDO
!
      DO A = 1, NT1AMX
         DO I = 1, NT1AMX
            IF (ABS(YIAJB(A,I)) .GT. 1.0D-8) THEN
               WRITE(LUPRI,*)'YIAJB         (',A,',',I,') = ',
     *                        YIAJB(A,I)
            ENDIF
            IF (ABS(YIAJB2(A,I)) .GT. 1.0D-8) THEN
               WRITE(LUPRI,*)'YIAJBTEMP     (',A,',',I,') = ',
     *                        YIAJB2(A,I)
            ENDIF
            IF (ABS(YIAJB3(A,I)) .GT. 1.0D-8) THEN
               WRITE(LUPRI,*)'YIAJBTEMP(CMO)(',A,',',I,') = ',
     *                        YIAJB3(A,I)
            ENDIF
!
            CALL FLSHFO(LUPRI)
!
            IF (ABS(YIAJB(A,I)-YIAJB2(A,I)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'WARNING 1. AND 2. VALUES DIFFER!!'
            ENDIF
            IF (ABS(YIAJB(A,I)-YIAJB3(A,I)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'WARNING 1. AND 3. VALUES DIFFER!!'
            ENDIF
            IF (ABS(YIAJB2(A,I)-YIAJB3(A,I)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'WARNING 2. AND 3. VALUES DIFFER!!'
            ENDIF
         ENDDO
      ENDDO
!
      RETURN
      END
C  /* Deck gcompare */
      SUBROUTINE GCOMPARE(INT1,INT2)
!
!     Written by K. Hald, April 2000.
!
!
      IMPLICIT NONE
!
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
!
#if defined (SYS_CRAY)
      REAL INT1(NORBT,NORBT,NORBT,NORBT)
      REAL INT2(NORBT,NORBT,NORBT,NORBT)
#else
      DOUBLE PRECISION INT1(NORBT,NORBT,NORBT,NORBT)
      DOUBLE PRECISION INT2(NORBT,NORBT,NORBT,NORBT)
#endif
!
      DO A = 1, NORBT
         DO I = 1, NORBT
            DO B = 1, NORBT
               DO J = 1, NORBT
                  IF (ABS(INT1(A,I,B,J)) .GT. 1.0D-8) THEN
                     WRITE(LUPRI,*)'G(lambda) = ',INT1(A,I,B,J)
                  ENDIF
                  IF (ABS(INT2(A,I,B,J)) .GT. 1.0D-8) THEN
                     WRITE(LUPRI,*)'G(CMO)    = ',INT2(A,I,B,J)
                  ENDIF
                  IF (ABS(INT1(A,I,B,J)-INT2(A,I,B,J)) .GT. 1.0D-8) THEN
                     WRITE(LUPRI,*)'THE TWO VALUES DIFFER'
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
      RETURN
      END
C  /* Deck intcompare2 */
      SUBROUTINE INTCOMPARE2(XINT1S,XINT1STEMP,XINT2S,XINT2STEMP,
     *                       XINT3S,XINT3STEMP,XINT4S,XINT4STEMP,
     *                       XINT5S,XINT5STEMP,XINT6S,XINT6STEMP)
!
!     Written by K. Hald, April 2000.
!
!
      IMPLICIT NONE
!
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
!
#if defined (SYS_CRAY)
      REAL XINT1S(NT1AMX,NVIRT,NVIRT)
      REAL XINT1STEMP(NT1AMX,NVIRT,NVIRT)
      REAL XINT2S(NT1AMX,NRHFT,NRHFT)
      REAL XINT2STEMP(NT1AMX,NRHFT,NRHFT)
      REAL XINT3S(NT1AMX,NVIRT,NVIRT)
      REAL XINT3STEMP(NT1AMX,NVIRT,NVIRT)
      REAL XINT4S(NT1AMX,NRHFT,NRHFT)
      REAL XINT4STEMP(NT1AMX,NRHFT,NRHFT)
      REAL XINT5S(NT1AMX,NVIRT,NVIRT)
      REAL XINT5STEMP(NT1AMX,NVIRT,NVIRT)
      REAL XINT6S(NT1AMX,NRHFT,NRHFT)
      REAL XINT6STEMP(NT1AMX,NRHFT,NRHFT)
#else
      DOUBLE PRECISION XINT1S(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT1STEMP(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT2S(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION XINT2STEMP(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION XINT3S(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT3STEMP(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT4S(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION XINT4STEMP(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION XINT5S(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT5STEMP(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT6S(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION XINT6STEMP(NT1AMX,NRHFT,NRHFT)
#endif
!
      DO A = 1, NT1AMX
        DO B = 1, NVIRT
          DO I = 1, NVIRT
             IF (ABS(XINT1S(A,B,I)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'XINT1S     = ',XINT1S(A,B,I)
             ENDIF
             IF (ABS(XINT1STEMP(A,B,I)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'XINT1STEMP = ',XINT1STEMP(A,B,I)
             ENDIF
             IF (ABS(XINT1S(A,B,I)-XINT1STEMP(A,B,I)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'For AI, B, C',A,B,I
                WRITE(LUPRI,*)'XINT1S & XINT1STEMP',XINT1S(A,B,I),
     *                                      XINT1STEMP(A,B,I)
                WRITE(LUPRI,*)'XINT1S AND XINT1STEMP DIFFER'
             ENDIF
          ENDDO
        ENDDO
      ENDDO
!
      DO A = 1, NT1AMX
        DO B = 1, NRHFT
          DO I = 1, NRHFT
             IF (ABS(XINT2S(A,B,I)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'XINT2S     = ',XINT2S(A,B,I)
             ENDIF
             IF (ABS(XINT2STEMP(A,B,I)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'XINT2STEMP = ',XINT2STEMP(A,B,I)
             ENDIF
             IF (ABS(XINT2S(A,B,I)-XINT2STEMP(A,B,I)) .GT. 1.0D-8)
     *                                                THEN
                WRITE(LUPRI,*)'For AI, B, C',A,B,I
                WRITE(LUPRI,*)'XINT2S & XINT2STEMP',XINT2S(A,B,I),
     *                                      XINT2STEMP(A,B,I)
                WRITE(LUPRI,*)'XINT2S AND XINT2STEMP DIFFER'
             ENDIF
          ENDDO
        ENDDO
      ENDDO
!
      DO A = 1, NT1AMX
        DO B = 1, NVIRT
          DO I = 1, NVIRT
             IF (ABS(XINT3S(A,B,I)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'XINT3S     = ',XINT3S(A,B,I)
             ENDIF
             IF (ABS(XINT3STEMP(A,B,I)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'XINT3STEMP = ',XINT3STEMP(A,B,I)
             ENDIF
             IF (ABS(XINT3S(A,B,I)-XINT3STEMP(A,B,I)) .GT. 1.0D-8)
     *                                                THEN
                WRITE(LUPRI,*)'For AI, B, C',A,B,I
                WRITE(LUPRI,*)'XINT3S & XINT3STEMP',XINT3S(A,B,I),
     *                                      XINT3STEMP(A,B,I)
                WRITE(LUPRI,*)'XINT3S AND XINT3STEMP DIFFER'
             ENDIF
          ENDDO
        ENDDO
      ENDDO
!
      DO A = 1, NT1AMX
        DO B = 1, NRHFT
          DO I = 1, NRHFT
             IF (ABS(XINT4S(A,B,I)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'XINT4S     = ',XINT4S(A,B,I)
             ENDIF
             IF (ABS(XINT4STEMP(A,B,I)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'XINT4STEMP = ',XINT4STEMP(A,B,I)
             ENDIF
             IF (ABS(XINT4S(A,B,I)-XINT4STEMP(A,B,I)) .GT. 1.0D-8)
     *                                                THEN
                WRITE(LUPRI,*)'For AI, B, C',A,B,I
                WRITE(LUPRI,*)'XINT4S & XINT4STEMP',XINT4S(A,B,I),
     *                                      XINT4STEMP(A,B,I)
                WRITE(LUPRI,*)'XINT4S AND XINT4STEMP DIFFER'
             ENDIF
          ENDDO
        ENDDO
      ENDDO
!
      DO A = 1, NT1AMX
        DO B = 1, NVIRT
          DO I = 1, NVIRT
             IF (ABS(XINT5S(A,B,I)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'XINT5S     = ',XINT5S(A,B,I)
             ENDIF
             IF (ABS(XINT5STEMP(A,B,I)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'XINT5STEMP = ',XINT5STEMP(A,B,I)
             ENDIF
             IF (ABS(XINT5S(A,B,I)-XINT5STEMP(A,B,I)) .GT. 1.0D-8)
     *                                                THEN
                WRITE(LUPRI,*)'For AI, B, C',A,B,I
                WRITE(LUPRI,*)'XINT5S & XINT5STEMP',XINT5S(A,B,I),
     *                                      XINT5STEMP(A,B,I)
                WRITE(LUPRI,*)'XINT5S AND XINT5STEMP DIFFER'
             ENDIF
          ENDDO
        ENDDO
      ENDDO
!
      DO A = 1, NT1AMX
        DO B = 1, NRHFT
          DO I = 1, NRHFT
             IF (ABS(XINT6S(A,B,I)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'XINT6S     = ',XINT6S(A,B,I)
             ENDIF
             IF (ABS(XINT6STEMP(A,B,I)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'XINT6STEMP = ',XINT6STEMP(A,B,I)
             ENDIF
             IF (ABS(XINT6S(A,B,I)-XINT6STEMP(A,B,I)) .GT. 1.0D-8)
     *                                                THEN
                WRITE(LUPRI,*)'For AI, B, C',A,B,I
                WRITE(LUPRI,*)'XINT6S & XINT6STEMP',XINT6S(A,B,I),
     *                                      XINT6STEMP(A,B,I)
                WRITE(LUPRI,*)'XINT6S AND XINT6STEMP DIFFER'
             ENDIF
          ENDDO
        ENDDO
      ENDDO
!
      RETURN
      END
C  /* Deck tcompare */
      SUBROUTINE TCOMPARE(XINT1T,XINT1TTEMP,XINT2T,XINT2TTEMP,
     *                    XINT3T,XINT3TTEMP,XINT4T,XINT4TTEMP)
!
!     Written by K. Hald, April 2000.
!
!
      IMPLICIT NONE
!
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
!
#if defined (SYS_CRAY)
      REAL XINT1T(NT1AMX,NVIRT,NVIRT)
      REAL XINT1TTEMP(NT1AMX,NVIRT,NVIRT)
      REAL XINT2T(NT1AMX,NRHFT,NRHFT)
      REAL XINT2TTEMP(NT1AMX,NRHFT,NRHFT)
      REAL XINT3T(NT1AMX,NVIRT,NVIRT)
      REAL XINT3TTEMP(NT1AMX,NVIRT,NVIRT)
      REAL XINT4T(NT1AMX,NRHFT,NRHFT)
      REAL XINT4TTEMP(NT1AMX,NRHFT,NRHFT)
#else
      DOUBLE PRECISION XINT1T(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT1TTEMP(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT2T(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION XINT2TTEMP(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION XINT3T(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT3TTEMP(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT4T(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION XINT4TTEMP(NT1AMX,NRHFT,NRHFT)
#endif
!
      DO A = 1, NT1AMX
         DO B = 1, NVIRT
            DO D = 1, NVIRT
            IF (ABS(XINT1T(A,B,D)) .GT. 1.0D-8) THEN
               WRITE(LUPRI,*)'XINT1T     = ',XINT1T(A,B,D)
            ENDIF
            IF (ABS(XINT1TTEMP(A,B,D)) .GT. 1.0D-8) THEN
               WRITE(LUPRI,*)'XINT1TTEMP = ',XINT1TTEMP(A,B,D)
            ENDIF
            IF (ABS(XINT1T(A,B,D)-XINT1TTEMP(A,B,D)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'For NAI, B, D',A,B,D
                WRITE(LUPRI,*)'XINT1T & XINT1TTEMP DIFFER!!!'
            ENDIF
            IF (ABS(XINT3T(A,B,D)) .GT. 1.0D-8) THEN
               WRITE(LUPRI,*)'XINT3T     = ',XINT3T(A,B,D)
            ENDIF
            IF (ABS(XINT3TTEMP(A,B,D)) .GT. 1.0D-8) THEN
               WRITE(LUPRI,*)'XINT3TTEMP = ',XINT3TTEMP(A,B,D)
            ENDIF
            IF (ABS(XINT3T(A,B,D)-XINT3TTEMP(A,B,D)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'For NAI, B, D',A,B,D
                WRITE(LUPRI,*)'XINT3T     = ',XINT3T(A,B,D)
                WRITE(LUPRI,*)'XINT3TTEMP = ',XINT3TTEMP(A,B,D)
                WRITE(LUPRI,*)'XINT3T & XINT3TTEMP DIFFER!!!'
            ENDIF
            ENDDO
         ENDDO
      ENDDO
!
      DO A = 1, NT1AMX
         DO B = 1, NRHFT
            DO D = 1, NRHFT
            IF (ABS(XINT2T(A,B,D)) .GT. 1.0D-8) THEN
               WRITE(LUPRI,*)'XINT2T     = ',XINT2T(A,B,D)
            ENDIF
            IF (ABS(XINT2TTEMP(A,B,D)) .GT. 1.0D-8) THEN
               WRITE(LUPRI,*)'XINT2TTEMP = ',XINT2TTEMP(A,B,D)
            ENDIF
            IF (ABS(XINT2T(A,B,D)-XINT2TTEMP(A,B,D)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'For NAI, B, D',A,B,D
                WRITE(LUPRI,*)'XINT2T     = ',XINT2T(A,B,D)
                WRITE(LUPRI,*)'XINT2TTEMP = ',XINT2TTEMP(A,B,D)
                WRITE(LUPRI,*)'XINT2T & XINT2TTEMP DIFFER!!!'
            ENDIF
            IF (ABS(XINT4T(A,B,D)) .GT. 1.0D-8) THEN
               WRITE(LUPRI,*)'XINT4T     = ',XINT4T(A,B,D)
            ENDIF
            IF (ABS(XINT4TTEMP(A,B,D)) .GT. 1.0D-8) THEN
               WRITE(LUPRI,*)'XINT4TTEMP = ',XINT4TTEMP(A,B,D)
            ENDIF
            IF (ABS(XINT4T(A,B,D)-XINT4TTEMP(A,B,D)) .GT. 1.0D-8) THEN
                WRITE(LUPRI,*)'For NAI, B, D',A,B,D
                WRITE(LUPRI,*)'XINT4T     = ',XINT4T(A,B,D)
                WRITE(LUPRI,*)'XINT4TTEMP = ',XINT4TTEMP(A,B,D)
                WRITE(LUPRI,*)'XINT4T & XINT4TTEMP DIFFER!!!'
            ENDIF
            ENDDO
         ENDDO
      ENDDO
!
      RETURN
      END
C  /* Deck cc3exci_int1 */
      SUBROUTINE CC3EXCI_INT1(CMO,AOINT,IDEL,INTTOT,
     *                        INTTOTTEMP,INTTOTTEMP2)
!
!     Written by K. Hald, April 2000
!
!     Calculate g(p,q,r,s) 
!
      IMPLICIT NONE
!
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h> 
      INTEGER INDEX, NAB, NDL, NCK, IDEL
!
#if defined (SYS_CRAY)
      REAL TWO, AOINT((NBAST*(NBAST+1))/2,NBAST)
      REAL CMO(NBAST,NORBT), ZERO
      REAL INTTOT(NORBT,NORBT,NORBT,NORBT)
      REAL INTTOTTEMP(NBAST,NBAST,NBAST,NBAST)
      REAL INTTOTTEMP2(NBAST,NBAST,NBAST,NBAST)
#else
      DOUBLE PRECISION TWO, AOINT((NBAST*(NBAST+1))/2,NBAST)
      DOUBLE PRECISION CMO(NBAST,NORBT), ZERO
      DOUBLE PRECISION INTTOT(NORBT,NORBT,NORBT,NORBT)
      DOUBLE PRECISION INTTOTTEMP(NBAST,NBAST,NBAST,NBAST)
      DOUBLE PRECISION INTTOTTEMP2(NBAST,NBAST,NBAST,NBAST)
#endif 
!
      PARAMETER (TWO = 2.0D0, ZERO = 0.0D0)
!
!
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
!
!
      CALL DZERO(INTTOTTEMP,NBAST*NBAST*NBAST*NBAST)
      CALL DZERO(INTTOTTEMP2,NBAST*NBAST*NBAST*NBAST)
!
      DO G = 1,NBAST
         DO B = 1,NBAST
            DO A = 1,NBAST
               NAB = INDEX(A,B)
!
               IF (.NOT.(AOINT(NAB,G) .EQ. ZERO)) THEN
               DO C = 1, NORBT
                 INTTOTTEMP(C,B,G,IDEL) = INTTOTTEMP(C,B,G,IDEL)
     *                     + AOINT(NAB,G)*CMO(A,C)
               ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!
      DO C = 1, NORBT
         DO G = 1,NBAST
            DO B = 1,NBAST
               IF (.NOT.(INTTOTTEMP(C,B,G,IDEL) .EQ. ZERO)) THEN
                 DO K = 1, NORBT
                 INTTOTTEMP2(C,K,G,IDEL) = INTTOTTEMP2(C,K,G,IDEL)
     *                     + INTTOTTEMP(C,B,G,IDEL)*CMO(B,K)
                 ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!
      CALL DZERO(INTTOTTEMP,NBAST*NBAST*NBAST*NBAST)
!
      DO C = 1, NORBT
         DO K = 1, NORBT
            DO G = 1, NBAST
               IF (.NOT.(INTTOTTEMP2(C,K,G,IDEL) .EQ. ZERO)) THEN
               DO D = 1, NORBT
                  INTTOTTEMP(C,K,D,IDEL) = INTTOTTEMP(C,K,D,IDEL)
     *                     +  INTTOTTEMP2(C,K,G,IDEL)*CMO(G,D)
               ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!
      DO C = 1, NORBT
         DO K = 1, NORBT
            DO D = 1, NORBT
               IF (.NOT. (INTTOTTEMP(C,K,D,IDEL) .EQ. ZERO)) THEN
               DO L = 1, NORBT
                  INTTOT(C,K,D,L) = INTTOT(C,K,D,L) +
     *                   INTTOTTEMP(C,K,D,IDEL)*CMO(IDEL,L)
               ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!
!
!      DO C = 1, NORBT
!        DO K = 1, NORBT
!          WRITE(LUPRI,*)'Output for C,K = ',C,K
!          CALL OUTPUT(INTTOT(C,K,1,1),1,NORBT,1,NORBT,NORBT,NORBT,1,6)
!        ENDDO
!      ENDDO
!
!
      RETURN
      END
