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
      SUBROUTINE CCSDT_ETA_NODDY(LISTL,IDLSTL,IOPTRES,
     &                           FOCKA,FOCKA_AO,FREQE,FOCK0,
     &                           OMEGA1,OMEGA2,
     &                           OMEGA1EFF,OMEGA2EFF,
     &                           IDOTS,DOTPROD,LISTDP,ITRAN,
     &                           NEATRAN,MXVEC,WORK,LWORK )
*---------------------------------------------------------------------*
*
*    Purpose: compute triples contribution to Eta vector
*
*             Eta^eff_1,2 = Eta_1,2(CCSD) + Eta_1,2(T^0_3)
*                               - Eta_3 A_3;1,2 (w_3 - w)^1 
*
*        
*     Written by Christof Haettig, April 2002 
*     based on different other noddy codes
*
*=====================================================================*
      IMPLICIT NONE  
#include <dummy.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccfield.h>
#include <ccorb.h>
#include <ccnoddy.h>

      LOGICAL LOCDBG, FD_TEST, XI_ONLY
      PARAMETER (LOCDBG=.FALSE., FD_TEST=.false., 
     &           XI_ONLY = .FALSE.)

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      CHARACTER*3 LISTL, LISTDP
      INTEGER LWORK, IDLSTL, IOPTRES, ITRAN, MXVEC, NEATRAN
      INTEGER IDOTS(MXVEC,NEATRAN)

#if defined (SYS_CRAY)
      REAL DOTPROD(MXVEC,NEATRAN), DDOT
      REAL WORK(LWORK), FREQL, FREQE, FREQC, SIGN
      REAL FOCKA_AO(NORBT,NORBT)
      REAL OMEGA1(*), OMEGA2(*), FOCKA(NORBT,NORBT)
      REAL OMEGA1EFF(*), OMEGA2EFF(*), FOCK0(NORBT,NORBT)
      REAL SIXTH, TWO, TCON, DCON, SCON, FF
#else
      DOUBLE PRECISION DOTPROD(MXVEC,NEATRAN), DDOT
      DOUBLE PRECISION FOCKA_AO(NORBT,NORBT)
      DOUBLE PRECISION WORK(LWORK), FREQL, FREQE, FREQC, SIGN
      DOUBLE PRECISION OMEGA1(*), OMEGA2(*), FOCKA(NORBT,NORBT)
      DOUBLE PRECISION OMEGA1EFF(*), OMEGA2EFF(*), FOCK0(NORBT,NORBT)
      DOUBLE PRECISION SIXTH, TWO, TCON, DCON, SCON, FF
#endif
      PARAMETER(SIXTH=1.0D0/6.0D0, TWO=2.0D0)

      CHARACTER*10 MODEL
      LOGICAL L2INCL
      INTEGER INDEX, LUSIFC, IOPT, ISYMD, ILLL, IDEL, ISYDIS, NIJ, IJ,
     &        IVEC, IDLSTC, ISYMC, LUFOCK, ILSTSYM, ISYML, LUTEMP, IDX
      INTEGER KSCR1, KFOCKD, KEND1, KT1AMP0, KLAMP0, KLAMH0, KFIELD,
     &        KINT1T0, KINT2T0, KINT1S0, KINT2S0, KXIAJB, KYIAJB,
     &        K0IOVVO, K0IOOVV, K0IOOOO, K0IVVVV, KOME1, KOME2, KDUM,
     &        KXINT, KEND2, LWRK2, KL1AM, KL2AM, KL3AM, KT3AM, KT2AM,
     &        K3AM, KEND3, LWRK3, KINT1SC, KINT2SC, KLAMPC, KLAMHC,
     &        KFOCKC, LWRK1, KE3AM, KTC3AM, KTC1AM, KTC2AM,
     &        KFDETA1, KFDETA2, KFDETA3, KEND4, LWRK4, KMMAT, KDIAM,
     &        KT3SCR, KFIELDAO, KFOCK0

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J 

      CALL QENTER('CCSDT_ETA_NODDY')
      KDUM = 1

      IF (DIRECT) CALL QUIT('CCSDT_ETA_NODDY: DIRECT NOT IMPLEMENTED')

      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'CCSDT_ETA_NODDY> Eta vector on entry:'
        CALL CC_PRP(OMEGA1,OMEGA2,1,1,1)
      END IF

*---------------------------------------------------------------------*
*     Memory allocation:
*---------------------------------------------------------------------*
      KSCR1   = 1
      KFOCKD  = KSCR1  + NT1AMX
      KFOCK0  = KFOCKD + NORBT
      KEND1   = KFOCK0 + NORBT*NORBT

      IF (NONHF) THEN
        KFIELD   = KEND1
        KFIELDAO = KFIELD   + NBAST*NBAST
        KEND1    = KFIELDAO + NBAST*NBAST
      END IF

      KT1AMP0 = KEND1
      KLAMP0  = KT1AMP0 + NT1AMX
      KLAMH0  = KLAMP0  + NLAMDT
      KEND1   = KLAMH0  + NLAMDT

      KINT1T0 = KEND1
      KINT2T0 = KINT1T0 + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2T0 + NRHFT*NRHFT*NT1AMX

      KINT1S0 = KEND1
      KINT2S0 = KINT1S0 + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2S0 + NRHFT*NRHFT*NT1AMX

      KXIAJB  = KEND1
      KYIAJB  = KXIAJB  + NT1AMX*NT1AMX
      KEND1   = KYIAJB  + NT1AMX*NT1AMX

      K0IOVVO = KEND1
      K0IOOVV = K0IOVVO + NRHFT*NVIRT*NVIRT*NRHFT
      K0IOOOO = K0IOOVV + NRHFT*NVIRT*NVIRT*NRHFT
      K0IVVVV = K0IOOOO + NRHFT*NRHFT*NRHFT*NRHFT
      KEND1   = K0IVVVV + NVIRT*NVIRT*NVIRT*NVIRT 

      KOME1   = KEND1
      KOME2   = KOME1  + NT1AMX
      KEND1   = KOME2  + NT1AMX*NT1AMX

      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_ETA_NODDY')
      ENDIF

*---------------------------------------------------------------------*
*     Get zeroth-order Lambda matrices:
*---------------------------------------------------------------------*
      IOPT   = 1
      Call CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AMP0),WORK(KDUM))

      Call LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AMP0),
     &            WORK(KEND1),LWRK1)

*---------------------------------------------------------------------*
*     Read precalculated integrals from file:
*---------------------------------------------------------------------*
      CALL CCSDT_READ_NODDY(.TRUE.,WORK(KFOCKD),WORK(KFOCK0),
     &                             WORK(KFIELD),WORK(KFIELDAO),
     &                      .TRUE.,WORK(KXIAJB),WORK(KYIAJB),
     &                      .TRUE.,WORK(KINT1S0),WORK(KINT2S0),
     &                      .TRUE.,WORK(KINT1T0),WORK(KINT2T0),
     &                      .TRUE.,WORK(K0IOVVO),WORK(K0IOOVV),
     &                             WORK(K0IOOOO),WORK(K0IVVVV),
     &                      NORBT,NLAMDT,NRHFT,NVIRT,NT1AMX)

*---------------------------------------------------------------------*
*     Some more memory allocations:
*---------------------------------------------------------------------*
      KL1AM  = KEND1
      KL2AM  = KL1AM + NT1AMX
      KL3AM  = KL2AM + NT1AMX*NT1AMX
      KEND2  = KL3AM + NT1AMX*NT1AMX*NT1AMX
      LWRK2  = LWORK - KEND2

      KT3AM  = KEND2
      KT2AM  = KT3AM + NT1AMX*NT1AMX*NT1AMX
      KEND3  = KT2AM + NT1AMX*NT1AMX
      LWRK3  = LWORK - KEND3

      KEND4  = KEND3
      IF (NONHF) THEN
        ! allocate scratch array needed in CCSDT_T03AM and 
        ! also in CCSDT_L03AM in case of finite difference runs
        KT3SCR = KEND4
        KEND4  = KT3SCR + NT1AMX*NT1AMX*NT1AMX
      END IF
      LWRK4  = LWORK - KEND4

      IF (LWRK4 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_ETA_NODDY')
      ENDIF
 
      ! read T^0 doubles amplitudes from file and square up
      IOPT  = 2
      Call CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KDUM),WORK(KT3AM))
      CALL CC_T2SQ(WORK(KT3AM),WORK(KT2AM),ISYM0)   

      ! read L^0 multipliers from file and square up doubles part
      ISYML = ILSTSYM(LISTL,IDLSTL)
      IOPT  = 3
      Call CC_RDRSP(LISTL,IDLSTL,ISYML,IOPT,MODEL,
     &              WORK(KL1AM),WORK(KT3AM))
      CALL CC_T2SQ(WORK(KT3AM),WORK(KL2AM),ISYM0)   

*---------------------------------------------------------------------*
*     Compute zeroth-order triples amplitudes:
*---------------------------------------------------------------------*
      LUTEMP = -1
      CALL GPOPEN(LUTEMP,FILNODT30,'UNKNOWN',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      READ(LUTEMP) (WORK(KT3AM+I-1), I=1,NT1AMX*NT1AMX*NT1AMX)
      CALL GPCLOSE(LUTEMP,'KEEP')

      IF (NONHF) THEN
        LUTEMP = -1
        CALL GPOPEN(LUTEMP,'T3AMPL','UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        REWIND LUTEMP
        WRITE (LUTEMP) (WORK(KT3AM-1+IDX),IDX=1,NT1AMX*NT1AMX*NT1AMX)
        CALL GPCLOSE(LUTEMP,'KEEP')
      END IF
*---------------------------------------------------------------------*
*     Compute L^0_3 multipliers:
*---------------------------------------------------------------------*
      IF (LISTL(1:3).EQ.'L0 ') THEN

        ! remember that CCSDT_L03AM returns -L3 !!
        CALL CCSDT_L03AM(WORK(KL3AM),WORK(KINT1T0),WORK(KINT2T0),
     *                  WORK(KXIAJB),FOCK0,WORK(KL1AM),
     *                  WORK(KL2AM),WORK(KSCR1),WORK(KFOCKD),
     *                  WORK(KFIELD),WORK(KT3SCR))

        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KL3AM),1)

      ELSE IF (LISTL(1:3).EQ.'L1 ' .OR. LISTL(1:3).EQ.'LE ' .OR.
     &         LISTL(1:3).EQ.'M1 ' .OR. LISTL(1:3).EQ.'N2 ' .OR.
     &         LISTL(1:3).EQ.'E0 '                              ) THEN

        CALL CCSDT_TBAR31_NODDY(WORK(KL3AM),FREQL,LISTL,IDLSTL,
     &                        WORK(KLAMP0),WORK(KLAMH0),
     &                        FOCK0,WORK(KFOCKD),WORK(KSCR1),
     &                        WORK(KXIAJB),WORK(KINT1T0),WORK(KINT2T0),
     &                        WORK(KEND3),LWRK3)

        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KL3AM),1)

      ELSE

        CALL QUIT('CCSDT_ETA_NODDY> LISTL NOT AVAILABLE:'//LISTL)
      
      END IF

*---------------------------------------------------------------------*
*     Compute contribution from <L_3|[[A,T^0_3],\tau_nu_1|HF>:
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KOME1),NT1AMX)

      CALL CCSDT_E1AM(WORK(KOME1),WORK(KL3AM),WORK(KT3AM),FOCKA)

      DO I = 1,NT1AMX
         OMEGA1(I) = OMEGA1(I) + WORK(KOME1+I-1)
      END DO

      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'CCSDT_ETA_NODDY> Contribution to Eta1:'
        CALL CC_PRP(WORK(KOME1),WORK,1,1,0)
      END IF

*---------------------------------------------------------------------*
*     Compute contribution from <L_3|[[A,T^0_2],\tau_nu_2]|HF>
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KOME2),NT1AMX*NT1AMX)

      CALL CCSDT_E2AM(WORK(KOME2),WORK(KL3AM),WORK(KT2AM),FOCKA)

      DO I = 1,NT1AMX
         DO J = 1,I
            IJ = NT1AMX*(I-1) + J
            NIJ = INDEX(I,J)
            OMEGA2(NIJ) = OMEGA2(NIJ) + WORK(KOME2+IJ-1)
         END DO
      END DO

*---------------------------------------------------------------------*
*     finite difference test:
*---------------------------------------------------------------------*
      IF (FD_TEST) THEN
         KFDETA1 = KEND3
         KFDETA2 = KFDETA1 + NT1AMX
         KFDETA3 = KFDETA2 + NT1AMX*NT1AMX
         KEND4   = KFDETA3 + NT1AMX*NT1AMX*NT1AMX
         LWRK4   = LWORK - KEND4
         IF (LWRK4 .LT. 0) THEN
           CALL QUIT('Insufficient space in CCSDT_ETA_NODDY')
         ENDIF

         CALL CCSDT_ETA_FD(WORK(KT1AMP0),WORK(KT2AM),WORK(KT3AM),
     &                     WORK(KL3AM),WORK(KL2AM),
     &                     WORK(KFDETA1),WORK(KFDETA2),WORK(KFDETA3),
     &                     FOCKA,FOCKA_AO,
     &                     WORK(KSCR1),WORK(KEND4),LWRK4)
 
        WRITE(LUPRI,*) 'CCSDT_ETA_NODDY> my ETA1:'
        CALL OUTPUT(WORK(KOME1),1,NVIRT,1,NRHFT,
     &              NVIRT,NRHFT,1,LUPRI)
        WRITE(LUPRI,*) 'CCSDT_ETA_NODDY> finite difference Eta1:'
        CALL OUTPUT(WORK(KFDETA1),1,NVIRT,1,NRHFT,
     &              NVIRT,NRHFT,1,LUPRI)
        CALL DAXPY(NT1AMX,-1.0D0,WORK(KOME1),1,WORK(KFDETA1),1)
        WRITE(LUPRI,*) 'CCSDT_ETA_NODDY> norm of difference:',
     &   DDOT(NT1AMX,WORK(KFDETA1),1,WORK(KFDETA1),1)

        WRITE(LUPRI,*) 'CCSDT_ETA_NODDY> my ETA2:'
        CALL OUTPUT(WORK(KOME2),1,NT1AMX,1,NT1AMX,
     &              NT1AMX,NT1AMX,1,LUPRI)
        WRITE(LUPRI,*) 'CCSDT_ETA_NODDY> finite difference Eta2:'
        CALL OUTPUT(WORK(KFDETA2),1,NT1AMX,1,NT1AMX,
     &              NT1AMX,NT1AMX,1,LUPRI)
        CALL DAXPY(NT1AMX*NT1AMX,-1.0D0,WORK(KOME2),1,WORK(KFDETA2),1)
        WRITE(LUPRI,*) 'CCSDT_ETA_NODDY> norm of difference:',
     &   DDOT(NT1AMX*NT1AMX,WORK(KFDETA2),1,WORK(KFDETA2),1)

        WRITE(LUPRI,*) 'CCSDT_ETA_NODDY> finite difference Eta3:'
        CALL OUTPUT(WORK(KFDETA3),1,NT1AMX*NT1AMX,1,NT1AMX,
     &              NT1AMX*NT1AMX,NT1AMX,1,LUPRI)
 
      END IF

*---------------------------------------------------------------------*
*     Compute triples result vector 
*       <L_2|[A,\tau_nu_3]|HF> + <L_3|[A,\tau_nu_3]|HF> ,
*     (CCSDT_E3AM accounts for the wrong sign of L_3)
*---------------------------------------------------------------------*
      ! overwrite T3 vector
      KE3AM  = KT3AM
  
      CALL DZERO(WORK(KE3AM),NT1AMX*NT1AMX*NT1AMX)

      L2INCL = .TRUE.
      CALL CCSDT_E3AM(WORK(KE3AM),WORK(KL2AM),WORK(KL3AM),FOCKA,L2INCL)

      IF (FD_TEST) THEN
        WRITE(LUPRI,*) 'CCSDT_ETA_NODDY> my Eta3:'
        CALL OUTPUT(WORK(KE3AM),1,NT1AMX*NT1AMX,1,NT1AMX,
     &              NT1AMX*NT1AMX,NT1AMX,1,LUPRI)
        CALL DAXPY(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KE3AM),1,
     &                 WORK(KFDETA3),1)
        WRITE(LUPRI,*) 'CCSDT_ETA_NODDY> norm of difference:',
     &   DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KFDETA3),1,WORK(KFDETA3),1)
      END IF
*---------------------------------------------------------------------*
*     Now we split:
*       for IOPTRES < 5 we compute the effective Eta vector
*       for IOPTRES = 5 we compute the contractions Eta^A T^B
*---------------------------------------------------------------------*
      IF (IOPTRES.GE.1 .AND. IOPTRES.LE.4) THEN

        CALL DCOPY(NT1AMX,OMEGA1,1,OMEGA1EFF,1)
        CALL DCOPY(NT2AMX,OMEGA2,1,OMEGA2EFF,1)

        CALL CC_LHPART_NODDY(OMEGA1EFF,OMEGA2EFF,WORK(KE3AM),-FREQE,
     &                       WORK(KFOCKD),WORK(KFIELD),
     &                       WORK(K0IOOOO),WORK(K0IOVVO),
     &                       WORK(K0IOOVV),WORK(K0IVVVV),
     &                       WORK(KT2AM),WORK(KINT1S0),WORK(KINT2S0),
     &                       WORK(KEND3),LWRK3)


      ELSE IF (IOPTRES.EQ.5) THEN

        SIGN = +1.0D0
        CALL CCDOTRSP_NODDY(WORK(KOME1),WORK(KOME2),WORK(KE3AM),SIGN,
     &                      ITRAN,LISTDP,IDOTS,DOTPROD,MXVEC,
     &                      WORK(KLAMP0),WORK(KLAMH0),
     &                      FOCK0,WORK(KFOCKD),
     &                      WORK(KXIAJB), WORK(KYIAJB),
     &                      WORK(KINT1T0),WORK(KINT2T0),
     &                      WORK(KINT1S0),WORK(KINT2S0),
     &                      'CCSDT_ETA_NODDY',LOCDBG,LOCDBG,.FALSE.,
     &                      WORK(KEND3),LWRK3)

      ELSE
        CALL QUIT('Illegal value for IOPTRES IN CCSDT_ETA_NODDY')
      END IF

      CALL QEXIT('CCSDT_ETA_NODDY')

      RETURN
      END 
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_ETA_NODDY                      *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_E3AM(E3AM,L2AM,L3AM,FOCKA,L2INCL)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples exictation amplitudes of eta vector
*
*             Eta_nu_3 = <L_2|[A,tau_nu_3]|HF>
*                          + <L_3|[A,tau_nu_3]|HF>
*
*     if L2INCL=.false. the first contribution (from L2) is skipped
*        
*     Written by Christof Haettig, April 2002 
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccfield.h>
#include <ccorb.h>

      LOGICAL L2INCL

#if defined (SYS_CRAY)
      REAL AIBJCK
      REAL E3AM(NT1AMX,NT1AMX,NT1AMX)
      REAL L2AM(NT1AMX,NT1AMX)
      REAL L3AM(NT1AMX,NT1AMX,NT1AMX)
      REAL FOCKA(NORBT,NORBT)
      REAL HALF
#else
      DOUBLE PRECISION AIBJCK
      DOUBLE PRECISION E3AM(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION L2AM(NT1AMX,NT1AMX)
      DOUBLE PRECISION L3AM(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION FOCKA(NORBT,NORBT)
      DOUBLE PRECISION HALF, THIRD
#endif
      PARAMETER ( HALF = 0.5D0, THIRD = 1.0D0/3.0D0 )

      INTEGER NK, NC, NCK, NJ, NB, NBJ, NBK, NI, NA, NAI, NL, ND, 
     &        NDJ, NBL, NCI, NAK, NBI, NAJ
 
      DO NK = 1,NRHFT
       DO NC = 1,NVIRT
        NCK = NVIRT*(NK-1) + NC
        
        DO NJ = 1,NRHFT
         DO NB = 1,NVIRT
          NBJ = NVIRT*(NJ-1) + NB
          NBK = NVIRT*(NK-1) + NB
        
          DO NI = 1,NRHFT
           DO NA = 1,NVIRT
             NAI = NVIRT*(NI-1) + NA
        
             AIBJCK = 0.0D0
        
             IF (L2INCL) THEN
                AIBJCK = AIBJCK + 
     &            ( + FOCKA(NK,NRHFT+NC) * L2AM(NAI,NBJ)
     &              - FOCKA(NJ,NRHFT+NC) * L2AM(NAI,NBK) )
             END IF
        
             DO ND = 1, NVIRT
                NDJ = NVIRT*(NJ-1) + ND
                AIBJCK = AIBJCK
     &            + HALF*L3AM(NAI,NDJ,NCK) * FOCKA(NRHFT+ND,NRHFT+NB)
             END DO
 
             DO NL = 1, NRHFT
                NBL = NVIRT*(NL-1) + NB
                AIBJCK = AIBJCK - HALF*L3AM(NAI,NBL,NCK) * FOCKA(NJ,NL)
             END DO

             E3AM(NAI,NBJ,NCK) = E3AM(NAI,NBJ,NCK) + AIBJCK
             E3AM(NAI,NCK,NBJ) = E3AM(NAI,NCK,NBJ) + AIBJCK
             E3AM(NBJ,NAI,NCK) = E3AM(NBJ,NAI,NCK) + AIBJCK
             E3AM(NCK,NAI,NBJ) = E3AM(NCK,NAI,NBJ) + AIBJCK
             E3AM(NBJ,NCK,NAI) = E3AM(NBJ,NCK,NAI) + AIBJCK
             E3AM(NCK,NBJ,NAI) = E3AM(NCK,NBJ,NAI) + AIBJCK
 
           END DO
          END DO
         END DO
        END DO
       END DO
      END DO
C
C------------------------------------------------------
C     Get rid of amplitudes that are not allowed.
C------------------------------------------------------

      CALL CCSDT_CLEAN_T3(E3AM,NT1AMX,NVIRT,NRHFT)

      RETURN
      END 
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_E3AM                           *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_E2AM(E2AM,L3AM,T2AM,FOCKA)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples correction to doubles Eta vector
*
*             Eta_nu_2(L3) = <L_3|[[A,T^0_2],tau_nu_2]|HF>
*
*     Written by Christof Haettig, April 2002 
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccorb.h>

#if defined (SYS_CRAY)
      REAL E2AM(NT1AMX,NT1AMX)
      REAL L3AM(NT1AMX,NT1AMX,NT1AMX)
      REAL T2AM(NT1AMX,NT1AMX)
      REAL FOCKA(NORBT,NORBT)
      REAL CONTRIB, HALF
#else
      DOUBLE PRECISION E2AM(NT1AMX,NT1AMX)
      DOUBLE PRECISION L3AM(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION T2AM(NT1AMX,NT1AMX)
      DOUBLE PRECISION FOCKA(NORBT,NORBT)
      DOUBLE PRECISION CONTRIB, HALF
#endif
      PARAMETER (HALF = 0.5D0)

      INTEGER NAI, NCK, NJ, NL, NB, ND, NBJ, NBL, NDL, NDJ

      DO NAI = 1, NT1AMX
       DO NCK = 1, NT1AMX
        DO NJ = 1, NRHFT
         DO NL = 1, NRHFT
          DO NB = 1, NVIRT
           DO ND = 1, NVIRT
             NBJ = NVIRT*(NJ-1) + NB
             NBL = NVIRT*(NL-1) + NB
             NDL = NVIRT*(NL-1) + ND
             NDJ = NVIRT*(NJ-1) + ND

             CONTRIB = 
     &        ( L3AM(NAI,NCK,NBL)*T2AM(NCK,NDL)*FOCKA(NJ,NRHFT+ND) +
     &          L3AM(NAI,NCK,NDJ)*T2AM(NCK,NDL)*FOCKA(NL,NRHFT+NB)   )

             E2AM(NAI,NBJ) = E2AM(NAI,NBJ) - CONTRIB
             E2AM(NBJ,NAI) = E2AM(NBJ,NAI) - CONTRIB

           END DO
          END DO
         END DO
        END DO
       END DO
      END DO


      RETURN
      END 
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_E2AM                           *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_E2AXC(XMMAT,T2CAM,DIA,FOCKA)
*---------------------------------------------------------------------*
*
*    Purpose: noddy code for M matrix intermediate used in the `real`
*             code to compute the contraction Eta x R
*
*     Written by Christof Haettig, May 2002 
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccorb.h>

#if defined (SYS_CRAY)
      REAL XMMAT(NRHFT,NRHFT,NT1AMX)
      REAL T2CAM(NT1AMX,NT1AMX)
      REAL FOCKA(NORBT,NORBT)
      REAL DIA(NT1AMX), CONTRIB
#else
      DOUBLE PRECISION XMMAT(NRHFT,NRHFT,NT1AMX)
      DOUBLE PRECISION T2CAM(NT1AMX,NT1AMX)
      DOUBLE PRECISION FOCKA(NORBT,NORBT)
      DOUBLE PRECISION DIA(NT1AMX), CONTRIB
#endif

      INTEGER NFN, NM, NI, NA, NAI, NAM

      CALL DZERO(DIA,NT1AMX)
      
      WRITE(LUPRI,*) 'M complete matrix:'
      CALL OUTPUT(XMMAT,1,NRHFT*NRHFT,1,NT1AMX,NRHFT*NRHFT,NT1AMX,
     &            1,LUPRI)

      WRITE(LUPRI,*) 'response amplitudes:'
      CALL OUTPUT(T2CAM,1,NT1AMX,1,NT1AMX,NT1AMX,NT1AMX,
     &            1,LUPRI)

      DO NFN = 1, NT1AMX

        DO NM = 1, NRHFT
          DO NI = 1, NRHFT
            DO NA = 1, NVIRT
              NAI = NVIRT*(NI-1) + NA
              NAM = NVIRT*(NM-1) + NA
              DIA(NAI) = DIA(NAI) + 
     &          XMMAT(NI,NM,NFN) * T2CAM(NAM,NFN)
            END DO
          END DO
        END DO

        write(lupri,*) 'NFN:',NFN
        write(lupri,*) 'amplitudes:'
        call output(t2cam(1,nfn),1,nvirt,1,nrhft,nvirt,nrhft,1,lupri)
        write(lupri,*) 'm matrix:'
        call output(xmmat(1,1,nfn),1,nrhft,1,nrhft,nrhft,nrhft,1,lupri)
        write(lupri,*) 'density:'
        call output(dia,1,nvirt,1,nrhft,nvirt,nrhft,1,lupri)

      END DO

      WRITE(LUPRI,*) 'CCSDT_E2AXC> Density:'
      CALL OUTPUT(DIA,1,NVIRT,1,NRHFT,NVIRT,NRHFT,1,LUPRI)

      CONTRIB = 0.0D0
      DO NI = 1, NRHFT
        DO NA = 1, NVIRT
          NAI = NVIRT*(NI-1) + NA
          CONTRIB = CONTRIB + DIA(NAI) * FOCKA(NI,NRHFT+NA)
        END DO
      END DO

      WRITE(LUPRI,*) 'CCSDT_E2AXC> Contribution:',CONTRIB

      RETURN
      END 
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_E2AXC                          *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_E1AM(E1AM,L3AM,T3AM,FOCKA)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples correction to singles Eta vector
*
*             Eta_nu_1(L3,T3) = <L_3|[[A,T^0_3],tau_nu_1]|HF>
*
*     Written by Christof Haettig, April 2002 
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccorb.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

#if defined (SYS_CRAY)
      REAL E1AM(NT1AMX), DDOT
      REAL L3AM(NT1AMX,NT1AMX,NT1AMX)
      REAL T3AM(NT1AMX,NT1AMX,NT1AMX)
      REAL FOCKA(NORBT,NORBT)
#else
      DOUBLE PRECISION E1AM(NT1AMX), DDOT
      DOUBLE PRECISION L3AM(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION T3AM(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION FOCKA(NORBT,NORBT)
#endif
   
      INTEGER NAI, NBJ, NK, NL, NCK, NC, ND, NCL, NDL, NDK
     
      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'CCSDT_E1AM> Norm^(E1AM) on input:',
     &    DDOT(NT1AMX,E1AM,1,E1AM,1)
        WRITE(LUPRI,*) 'CCSDT_E1AM> Norm^(L3AM) on input:',
     &    DDOT(NT1AMX*NT1AMX*NT1AMX,L3AM,1,L3AM,1)
        WRITE(LUPRI,*) 'CCSDT_E1AM> Norm^(T3AM) on input:',
     &    DDOT(NT1AMX*NT1AMX*NT1AMX,T3AM,1,T3AM,1)
        WRITE(LUPRI,*) 'CCSDT_E1AM> Norm^(FOCKA) on input:',
     &    DDOT(NORBT*NORBT,FOCKA,1,FOCKA,1)
      END IF

      DO NAI = 1, NT1AMX
       DO NBJ = 1, NT1AMX
        DO NK = 1, NRHFT
         DO NL = 1, NRHFT
          DO NC = 1, NVIRT
           DO ND = 1, NVIRT
             NCL = NVIRT*(NL-1) + NC
             NCK = NVIRT*(NK-1) + NC
             NDL = NVIRT*(NL-1) + ND
             NDK = NVIRT*(NK-1) + ND

             E1AM(NDL) = E1AM(NDL) - 0.5D0 * 
     &        ( L3AM(NAI,NBJ,NDK) * FOCKA(NL,NRHFT+NC) +
     &          L3AM(NAI,NBJ,NCL) * FOCKA(NK,NRHFT+ND) ) *
     &          T3AM(NAI,NBJ,NCK)

           END DO
          END DO
         END DO
        END DO
       END DO
      END DO

      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'CCSDT_E1AM> E1AM on output:'
        CALL CC_PRP(E1AM,E1AM,1,1,0)
      END IF

      RETURN
      END 
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_E1AM                           *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_ETA_FD(T1AM,T2AM,T3AM,L3AM,L2AM,ETA1,ETA2,ETA3,
     &                        FOCKA,FOCKA_AO,SCR1,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: Construct Eta vector by finite difference on Xi vector
*
*    Written by Christof Haettig, April 2002 
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccorb.h>

#if defined (SYS_CRAY)
      REAL DELTA, ZERO, HALF, SIXTH
#else
      DOUBLE PRECISION DELTA, ZERO, HALF, SIXTH
#endif
      PARAMETER(DELTA=1.0D-6, ZERO=0.0D0, HALF=0.5D0, SIXTH=1.0D0/6.0D0)

      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL WORK(LWORK), SCR1(NT1AMX),
     &   T1AM(NT1AMX), T2AM(NT1AMX,NT1AMX), T3AM(NT1AMX,NT1AMX,NT1AMX),
     &   L2AM(NT1AMX,NT1AMX), L3AM(NT1AMX,NT1AMX,NT1AMX),
     &   ETA1(NT1AMX), ETA2(NT1AMX,NT1AMX), ETA3(NT1AMX,NT1AMX,NT1AMX),
     &   FOCKA(NORBT,NORBT), FOCKA_AO(NORBT,NORBT),
     &   TAIBJ, TAIBJCK, DDOT, EAIBJ
#else
      DOUBLE PRECISION WORK(LWORK), SCR1(NT1AMX),
     &   T1AM(NT1AMX), T2AM(NT1AMX,NT1AMX), T3AM(NT1AMX,NT1AMX,NT1AMX),
     &   L2AM(NT1AMX,NT1AMX), L3AM(NT1AMX,NT1AMX,NT1AMX),
     &   ETA1(NT1AMX), ETA2(NT1AMX,NT1AMX), ETA3(NT1AMX,NT1AMX,NT1AMX),
     &   FOCKA(NORBT,NORBT), FOCKA_AO(NORBT,NORBT),
     &   TAIBJ, TAIBJCK, DDOT, EAIBJCK, EAIBJ
#endif
  
      INTEGER KT2AM, KXI2, KXI3, NAI, NBJ, NCK, KEND1, LWRK1, KT1AM,
     &        KFOCKA, KLAMP0, KLAMH0, KFOCKD, NA, NI, NB, NBI, NC, 
     &        NCI, NJ, NAJ, NK, NAK, KL3AM


      CALL QUIT('CCSDT_ETA_FD needs to be addapted for intermediates.')
*---------------------------------------------------------------------*
*     Allocations:
*---------------------------------------------------------------------*
      KL3AM = 1
      KFOCKD= KL3AM  + NT1AMX*NT1AMX*NT1AMX
      KFOCKA= KFOCKD + NORBT
      KLAMP0= KFOCKA + NORBT*NORBT
      KLAMH0= KLAMP0 + NLAMDT
      KT1AM = KLAMH0 + NLAMDT
      KT2AM = KT1AM + NT1AMX
      KXI2  = KT2AM + NT1AMX*NT1AMX
      KXI3  = KXI2  + NT1AMX*NT1AMX
      KEND1 = KXI3  + NT1AMX*NT1AMX*NT1AMX
      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_ETA_FD')
      ENDIF

      DO I = 1, NRHFT
         WORK(KFOCKD-1+I) = 0.0d0
      ENDDO
      DO I = 1, NVIRT
         WORK(KFOCKD-1+NRHFT+I) = 1.0d0/3.0d0
      ENDDO

      ! turn sign of T3, because Xi vector routines work with -T3
      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,T3AM,1)

C     WRITE(LUPRI,*) 'in finite difference norm^2(L3AM):',
C    &  DDOT(NT1AMX*NT1AMX*NT1AMX,L3AM,1,L3AM,1)
C     WRITE(LUPRI,*) 'in finite difference norm^2(T3AM):',
C    &  DDOT(NT1AMX*NT1AMX*NT1AMX,T3AM,1,T3AM,1)
C     WRITE(LUPRI,*) 'in finite difference norm^2(FOCKA):',
C    &  DDOT(NORBT*NORBT,FOCKA,1,FOCKA,1)

C     write(lupri,*) 'in finite difference focka:'
C     call output(focka,1,norbt,1,norbt,norbt,norbt,1,lupri)
C     write(lupri,*) 'in finite difference T3AM:'
C     call output(t3am,1,NT1AMX*NT1AMX,1,NT1AMX,
C    &             NT1AMX*NT1AMX,NT1AMX,1,lupri)
*---------------------------------------------------------------------*
*     take the derivative of Xksi_3 w.r.t. T1
*     use, that Xksi_3 is linear in T1
*---------------------------------------------------------------------*
      CALL DCOPY(NT1AMX,T1AM,1,WORK(KT1AM),1) 
      CALL DZERO(WORK(KT2AM),NT1AMX*NT1AMX) 
      CALL DZERO(ETA1,NT1AMX)

      DO NAI = 1, NT1AMX
      
        CALL DZERO(WORK(KT1AM),NT1AMX) 
        WORK(KT1AM-1+NAI) = T1AM(NAI) - DELTA*HALF

        CALL LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AM),
     &              WORK(KEND1),LWRK1)
        CALL DCOPY(NORBT*NORBT,FOCKA_AO,1,WORK(KFOCKA),1)
        CALL CC_FCKMO(WORK(KFOCKA),WORK(KLAMP0),WORK(KLAMH0),
     &                WORK(KEND1),LWRK1,1,1,1)

        CALL DZERO(WORK(KXI3),NT1AMX*NT1AMX*NT1AMX)
        CALL CCSDT_XKSI3(WORK(KXI3),WORK(KFOCKA),T3AM,T2AM)
        
         ETA1(NAI) = ETA1(NAI) -
     &      SIXTH*DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KXI3),1,L3AM,1)

        WORK(KT1AM-1+NAI) = T1AM(NAI) + DELTA*HALF

        CALL LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AM),
     &              WORK(KEND1),LWRK1)
        CALL DCOPY(NORBT*NORBT,FOCKA_AO,1,WORK(KFOCKA),1)
        CALL CC_FCKMO(WORK(KFOCKA),WORK(KLAMP0),WORK(KLAMH0),
     &                 WORK(KEND1),LWRK1,1,1,1)

        CALL DZERO(WORK(KXI3),NT1AMX*NT1AMX*NT1AMX)
        CALL CCSDT_XKSI3(WORK(KXI3),WORK(KFOCKA),T3AM,T2AM)

         ETA1(NAI) = ETA1(NAI) +
     &      SIXTH*DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KXI3),1,L3AM,1)

        WORK(KT1AM-1+NAI) = T1AM(NAI)
      END DO
 

*---------------------------------------------------------------------*
*     take the derivative of Xksi2 & Xksi_3 w.r.t. T3
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KT2AM),NT1AMX*NT1AMX)
      CALL DZERO(ETA3,NT1AMX*NT1AMX*NT1AMX)

C     CALL DZERO(WORK(KL3AM),NT1AMX*NT1AMX*NT1AMX)
C     WORK(KL3AM-1+13) = L3AM(13,1,1)
      CALL DCOPY(NT1AMX*NT1AMX*NT1AMX,L3AM,1,WORK(KL3AM),1)

      if (.true.) then

      DO NAI = 1, NT1AMX
        DO NBJ = 1, NT1AMX
          DO NCK = 1, NT1AMX
            TAIBJCK = T3AM(NAI,NBJ,NCK)
            T3AM(NAI,NBJ,NCK) = TAIBJCK - DELTA*HALF

            CALL DZERO(WORK(KXI2),NT1AMX*NT1AMX)
            CALL CCSDT_XKSI2(WORK(KXI2),FOCKA,T3AM)

            CALL DZERO(WORK(KXI3),NT1AMX*NT1AMX*NT1AMX)
            CALL CCSDT_XKSI3(WORK(KXI3),FOCKA,T3AM,WORK(KT2AM))

            EAIBJCK =  +
     &         HALF*DDOT(NT1AMX*NT1AMX,WORK(KXI2),1,L2AM,1) +
     &         SIXTH*DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KXI3),1,L3AM,1)

            T3AM(NAI,NBJ,NCK) = TAIBJCK + DELTA*HALF

            CALL DZERO(WORK(KXI2),NT1AMX*NT1AMX)
            CALL CCSDT_XKSI2(WORK(KXI2),FOCKA,T3AM)

            CALL DZERO(WORK(KXI3),NT1AMX*NT1AMX*NT1AMX)
            CALL CCSDT_XKSI3(WORK(KXI3),FOCKA,T3AM,WORK(KT2AM))

            EAIBJCK = EAIBJCK -
     &         HALF*DDOT(NT1AMX*NT1AMX,WORK(KXI2),1,L2AM,1) -
     &         SIXTH*DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KXI3),1,L3AM,1) 

            ETA3(NAI,NBJ,NCK) = ETA3(NAI,NBJ,NCK) + EAIBJCK
            ETA3(NAI,NCK,NBJ) = ETA3(NAI,NCK,NBJ) + EAIBJCK
            ETA3(NBJ,NAI,NCK) = ETA3(NBJ,NAI,NCK) + EAIBJCK
            ETA3(NCK,NAI,NBJ) = ETA3(NCK,NAI,NBJ) + EAIBJCK
            ETA3(NBJ,NCK,NAI) = ETA3(NBJ,NCK,NAI) + EAIBJCK
            ETA3(NCK,NBJ,NAI) = ETA3(NCK,NBJ,NAI) + EAIBJCK

            T3AM(NAI,NBJ,NCK) = TAIBJCK
          END DO
        END DO
      END DO
      end if

*---------------------------------------------------------------------*
*     take the derivative of Xksi_3 w.r.t. T2
*---------------------------------------------------------------------*
      CALL DCOPY(NT1AMX*NT1AMX,T2AM,1,WORK(KT2AM),1)
      CALL DZERO(ETA2,NT1AMX*NT1AMX)

      if (.true.) then

      DO NAI = 1, NT1AMX
        DO NBJ = 1, NT1AMX
            TAIBJ = T2AM(NAI,NBJ)
            WORK(KT2AM-1+(NBJ-1)*NT1AMX+NAI) = TAIBJ - DELTA*HALF

            CALL DZERO(WORK(KXI3),NT1AMX*NT1AMX*NT1AMX)
            CALL CCSDT_XKSI3(WORK(KXI3),FOCKA,T3AM,WORK(KT2AM))

            EAIBJ = -
     &         SIXTH*DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KXI3),1,L3AM,1)


            WORK(KT2AM-1+(NBJ-1)*NT1AMX+NAI) = TAIBJ + DELTA*HALF

            CALL DZERO(WORK(KXI3),NT1AMX*NT1AMX*NT1AMX)
            CALL CCSDT_XKSI3(WORK(KXI3),FOCKA,T3AM,WORK(KT2AM))

            EAIBJ = EAIBJ +
     &         SIXTH*DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KXI3),1,L3AM,1)

            ETA2(NAI,NBJ) = ETA2(NAI,NBJ) + EAIBJ
            ETA2(NBJ,NAI) = ETA2(NBJ,NAI) + EAIBJ

            WORK(KT2AM-1+(NBJ-1)*NT1AMX+NAI) = TAIBJ
        END DO
      END DO

      end if

*---------------------------------------------------------------------*
*     scale result with 1/DELTA and remove forbidden triples elements
*---------------------------------------------------------------------*
      CALL DSCAL(NT1AMX,              1.0D0/DELTA,ETA1,1)
      CALL DSCAL(NT1AMX*NT1AMX,       1.0D0/DELTA,ETA2,1)
      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,1.0D0/DELTA,ETA3,1)

      CALL CCSDT_CLEAN_T3(ETA3,NT1AMX,NVIRT,NRHFT)

      ! restor correct sign of T3
      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,T3AM,1)

      RETURN
      END 
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_ETA_FD                         *
*---------------------------------------------------------------------*
