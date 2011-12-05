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
      SUBROUTINE CCSDT_GMAT_NODDY(LISTL,IDLSTL,LISTB,IDLSTB,
     &                            LISTC,IDLSTC,
     &                            OMEGA1,OMEGA2,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples contribution to G matrix transformation
*
*  (G T^B T^C)^eff_1,2 = (G T^B T^C)_1,2(CCSD) 
*                            + (G T^B T^C)_1,2(L_3)
*        
*     Written by Christof Haettig, April 2002 
*     based on CCSDT_FMAT_NODDY
*
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

      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.FALSE.)

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      CHARACTER*3 LISTL, LISTB, LISTC
      INTEGER LWORK, IDLSTL, IDLSTB, IDLSTC
#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL OMEGA1(NT1AMX), OMEGA2(NT2AMX)
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION OMEGA1(NT1AMX), OMEGA2(NT2AMX)
#endif

      CHARACTER*10 MODEL
      INTEGER KXINT, KXIAJB, KYIAJB, KT1AMP0, KLAMP0, KLAMH0, KEND1,
     &        KFOCK0, LWRK1, KLAMPB, KLAMHB, KLAMPC, KLAMHC, KT1AMPB,
     &        KT1AMPC, KINT1T0, KINT2T0, KINT1SBC, KINT2SBC,
     &        KBIOVVO, KBIOOVV, KBIOOOO, KBIVVVV,
     &        KCIOVVO, KCIOOVV, KCIOOOO, KCIVVVV,
     &        KBCIOVVO, KBCIOOVV, KBCIOOOO, KBCIVVVV,
     &        KOME1, KOME2, KL1AM, KL2AM, KL3AM, KT2AM, KFOCKD,
     &        KDUM, IRECNR, KSCR1
      INTEGER ISYMD, ILLL, IDEL, ISYDIS, IDUMMY, IJ, NIJ, LUSIFC, INDEX,
     &        ISYMC, ISYMB, LUFOCK, KEND2, LWRK2, IOPT, ILSTSYM

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J 

      CALL QENTER('CCSDT_GMAT_NODDY')

      IF (DIRECT) CALL QUIT('CCSDT_GMAT_NODDY: DIRECT NOT IMPLEMENTED')

*---------------------------------------------------------------------*
*     Memory allocation:
*---------------------------------------------------------------------*
      KSCR1   = 1
      KFOCKD  = KSCR1  + NT1AMX
      KEND1   = KFOCKD + NORBT

      KFOCK0  = KEND1
      KEND1   = KFOCK0  + NORBT*NORBT

      KT1AMP0 = KEND1
      KLAMP0  = KT1AMP0 + NT1AMX
      KLAMH0  = KLAMP0  + NLAMDT
      KEND1   = KLAMH0  + NLAMDT

      KT1AMPB = KEND1
      KLAMPB  = KT1AMPB + NT1AMX
      KLAMHB  = KLAMPB  + NLAMDT
      KEND1   = KLAMHB  + NLAMDT

      KT1AMPC = KEND1
      KLAMPC  = KT1AMPC + NT1AMX
      KLAMHC  = KLAMPC  + NLAMDT
      KEND1   = KLAMHC  + NLAMDT

      KXIAJB  = KEND1
      KYIAJB  = KXIAJB  + NT1AMX*NT1AMX
      KEND1   = KYIAJB  + NT1AMX*NT1AMX

      KINT1T0 = KEND1
      KINT2T0 = KINT1T0 + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2T0 + NRHFT*NRHFT*NT1AMX

      KBIOVVO = KEND1
      KBIOOVV = KBIOVVO + NRHFT*NVIRT*NVIRT*NRHFT
      KBIOOOO = KBIOOVV + NRHFT*NVIRT*NVIRT*NRHFT
      KBIVVVV = KBIOOOO + NRHFT*NRHFT*NRHFT*NRHFT
      KEND1   = KBIVVVV + NVIRT*NVIRT*NVIRT*NVIRT

      KCIOVVO = KEND1
      KCIOOVV = KCIOVVO + NRHFT*NVIRT*NVIRT*NRHFT
      KCIOOOO = KCIOOVV + NRHFT*NVIRT*NVIRT*NRHFT
      KCIVVVV = KCIOOOO + NRHFT*NRHFT*NRHFT*NRHFT
      KEND1   = KCIVVVV + NVIRT*NVIRT*NVIRT*NVIRT

      KBCIOVVO = KEND1
      KBCIOOVV = KBCIOVVO + NRHFT*NVIRT*NVIRT*NRHFT
      KBCIOOOO = KBCIOOVV + NRHFT*NVIRT*NVIRT*NRHFT
      KBCIVVVV = KBCIOOOO + NRHFT*NRHFT*NRHFT*NRHFT
      KEND1    = KBCIVVVV + NVIRT*NVIRT*NVIRT*NVIRT

      KINT1SBC = KEND1
      KINT2SBC = KINT1SBC + NT1AMX*NVIRT*NVIRT
      KEND1    = KINT2SBC + NRHFT*NRHFT*NT1AMX

      KOME1   = KEND1
      KOME2   = KOME1  + NT1AMX
      KEND1   = KOME2  + NT1AMX*NT1AMX

      KL1AM   = KEND1
      KL2AM   = KL1AM  + NT1AMX
      KL3AM   = KL2AM  + NT1AMX*NT1AMX
      KEND1   = KL3AM + NT1AMX*NT1AMX*NT1AMX

      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_GMAT_NODDY')
      ENDIF

*---------------------------------------------------------------------*
*     Read SCF orbital energies from file:
*---------------------------------------------------------------------*
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBT)
      CALL GPCLOSE(LUSIFC,'KEEP')

*---------------------------------------------------------------------*
*     Get zeroth-order Lambda matrices:
*---------------------------------------------------------------------*
      IOPT   = 1
      KDUM = KEND1
      Call CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AMP0),WORK(KDUM))

      Call LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AMP0),
     &            WORK(KEND1),LWRK1)

*---------------------------------------------------------------------*
*     Get response  Lambda matrices:
*---------------------------------------------------------------------*
      ISYMB = ILSTSYM(LISTB,IDLSTB)
      IOPT  = 1
      Call CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &              WORK(KT1AMPB),WORK(KDUM))

      CALL CCLR_LAMTRA(WORK(KLAMP0),WORK(KLAMPB),
     &                 WORK(KLAMH0),WORK(KLAMHB),WORK(KT1AMPB),ISYMB)

      ISYMC = ILSTSYM(LISTC,IDLSTC)
      IOPT  = 1
      Call CC_RDRSP(LISTC,IDLSTC,ISYMC,IOPT,MODEL,
     &              WORK(KT1AMPC),WORK(KDUM))

      CALL CCLR_LAMTRA(WORK(KLAMP0),WORK(KLAMPC),
     &                 WORK(KLAMH0),WORK(KLAMHC),WORK(KT1AMPC),ISYMC)

*---------------------------------------------------------------------*
*     read zeroth-order AO Fock matrix from file: 
*---------------------------------------------------------------------*
      LUFOCK = -1
      CALL GPOPEN(LUFOCK,'CC_FCKH','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUFOCK)
      READ(LUFOCK) (WORK(KFOCK0-1+I),I=1,N2BST(ISYM0))
      CALL GPCLOSE(LUFOCK,'KEEP')

*---------------------------------------------------------------------*
*     Compute integrals needed for the following contributions:
*---------------------------------------------------------------------*

      CALL DZERO(WORK(KBIOVVO),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KBIOOVV),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KBIOOOO),NRHFT*NRHFT*NRHFT*NRHFT)
      CALL DZERO(WORK(KBIVVVV),NVIRT*NVIRT*NVIRT*NVIRT)

      CALL DZERO(WORK(KCIOVVO),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KCIOOVV),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KCIOOOO),NRHFT*NRHFT*NRHFT*NRHFT)
      CALL DZERO(WORK(KCIVVVV),NVIRT*NVIRT*NVIRT*NVIRT)

      CALL DZERO(WORK(KBCIOVVO),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KBCIOOVV),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KBCIOOOO),NRHFT*NRHFT*NRHFT*NRHFT)
      CALL DZERO(WORK(KBCIVVVV),NVIRT*NVIRT*NVIRT*NVIRT)

      CALL DZERO(WORK(KXIAJB),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KYIAJB),NT1AMX*NT1AMX)

      CALL DZERO(WORK(KINT1T0),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2T0),NRHFT*NRHFT*NT1AMX)

      CALL DZERO(WORK(KINT1SBC),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2SBC),NRHFT*NRHFT*NT1AMX)

      DO ISYMD = 1, NSYM
         DO ILLL = 1,NBAS(ISYMD)
            IDEL   = IBAS(ISYMD) + ILLL
            ISYDIS = MULD2H(ISYMD,ISYMOP)
 
C           ----------------------------
C           Work space allocation no. 2.
C           ----------------------------
            KXINT  = KEND1
            KEND2  = KXINT + NDISAO(ISYDIS)
            LWRK2  = LWORK - KEND2
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
               CALL QUIT('Insufficient space in CCSDT_GMAT_NODDY')
            ENDIF
 
C           ---------------------------
C           Read in batch of integrals.
C           ---------------------------
            CALL CCRDAO(WORK(KXINT),IDEL,1,WORK(KEND2),LWRK2,
     *                  IRECNR,DIRECT)
 
C           ----------------------------------
C           Calculate integrals needed in CC3:
C           ----------------------------------
            CALL CCSDT_TRAN1(WORK(KINT1T0),WORK(KINT2T0),
     &                       WORK(KLAMP0),WORK(KLAMH0),
     &                       WORK(KXINT),IDEL)

            CALL CC3_TRAN2(WORK(KXIAJB),WORK(KYIAJB),
     &                     WORK(KLAMP0),WORK(KLAMH0),
     &                     WORK(KXINT),IDEL)

            ! XINT1S = XINT1S + (C-barB K-barC|B D)
            ! XINT2S = XINT2S + (C-barB K-barC|L J)
            CALL CCSDT_TRAN3_R(WORK(KINT1SBC),WORK(KINT2SBC),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPB),WORK(KLAMHC),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KXINT),IDEL)

            ! XINT1S = XINT1S + (C-barB K|B-barC D)
            ! XINT2S = XINT2S + (C-barB K|L J-barC)
            CALL CCSDT_TRAN3_R(WORK(KINT1SBC),WORK(KINT2SBC),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPB),WORK(KLAMH0),
     &                         WORK(KLAMPC),WORK(KLAMHC),
     &                         WORK(KXINT),IDEL)

            ! XINT1S = XINT1S + (C-barC K-barB|B D)
            ! XINT2S = XINT2S + (C-barC K-barB|L J)
            CALL CCSDT_TRAN3_R(WORK(KINT1SBC),WORK(KINT2SBC),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPC),WORK(KLAMHB),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KXINT),IDEL)

            ! XINT1S = XINT1S + (C K-barB|B-barC D)
            ! XINT2S = XINT2S + (C K-barB|L J-barC)
            CALL CCSDT_TRAN3_R(WORK(KINT1SBC),WORK(KINT2SBC),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMP0),WORK(KLAMHB),
     &                         WORK(KLAMPC),WORK(KLAMHC),
     &                         WORK(KXINT),IDEL)

            ! XINT1S = XINT1S + (C-barC K|B-barB D)
            ! XINT2S = XINT2S + (C-barC K|L J-barB)
            CALL CCSDT_TRAN3_R(WORK(KINT1SBC),WORK(KINT2SBC),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPC),WORK(KLAMH0),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KXINT),IDEL)

            ! XINT1S = XINT1S + (C K-barC|B-barB D)
            ! XINT2S = XINT2S + (C K-barC|L J-barB)
            CALL CCSDT_TRAN3_R(WORK(KINT1SBC),WORK(KINT2SBC),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMP0),WORK(KLAMHC),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KXINT),IDEL)

C           ----------------------------------------------
C             (OV|VO)-B  (OO|VV)-B  (OO|OO)-B  (VV|VV)-B
C           ----------------------------------------------
            CALL CCFOP_TRAN1_R(WORK(KBIOVVO),WORK(KBIOOVV),
     &                         WORK(KBIOOOO),WORK(KBIVVVV),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KXINT),IDEL)

            CALL CCFOP_TRAN1_R(WORK(KBIOVVO),WORK(KBIOOVV),
     &                         WORK(KBIOOOO),WORK(KBIVVVV),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KXINT),IDEL)

C           ----------------------------------------------
C             (OV|VO)-C  (OO|VV)-C  (OO|OO)-C  (VV|VV)-C
C           ----------------------------------------------
            CALL CCFOP_TRAN1_R(WORK(KCIOVVO),WORK(KCIOOVV),
     &                         WORK(KCIOOOO),WORK(KCIVVVV),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPC),WORK(KLAMHC),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KXINT),IDEL)

            CALL CCFOP_TRAN1_R(WORK(KCIOVVO),WORK(KCIOOVV),
     &                         WORK(KCIOOOO),WORK(KCIVVVV),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPC),WORK(KLAMHC),
     &                         WORK(KXINT),IDEL)

C           ----------------------------------------------
C           (OV|VO)-BC  (OO|VV)-BC  (OO|OO)-BC  (VV|VV)-BC
C           ----------------------------------------------
            CALL CCFOP_TRAN1_R(WORK(KBCIOVVO),WORK(KBCIOOVV),
     &                         WORK(KBCIOOOO),WORK(KBCIVVVV),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KLAMPC),WORK(KLAMHC),
     &                         WORK(KXINT),IDEL)

            CALL CCFOP_TRAN1_R(WORK(KBCIOVVO),WORK(KBCIOOVV),
     &                         WORK(KBCIOOOO),WORK(KBCIVVVV),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPC),WORK(KLAMHC),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KXINT),IDEL)

         END DO   
      END DO  

*---------------------------------------------------------------------*
*     Compute L^0_3 multipliers:
*---------------------------------------------------------------------*

      IF (LISTL(1:3).EQ.'L0 ') THEN

        CALL CCSDT_L3AM(WORK(KL3AM),WORK(KINT1T0),WORK(KINT2T0),
     *                  WORK(KXIAJB),WORK(KFOCK0),WORK(KL1AM),
     *                  WORK(KL2AM),WORK(KSCR1),WORK(KFOCKD))

      ELSE

        CALL QUIT('CCSDT_GMAT_NODDY> LISTL NOT AVAILABLE:'//LISTL)
      
      END IF

*---------------------------------------------------------------------*
*     Compute contribution from  <L_3|[[H^BC,T^0_2],\tau_nu_1]|HF>
*                          and   <L_3|[[H^B, T^C_2],\tau_nu_1]|HF>
*                          and   <L_3|[[H^C, T^B_2],\tau_nu_1]|HF>
*---------------------------------------------------------------------*
      KT2AM  = KEND1
      KEND1  = KT2AM + NT1AMX*NT1AMX
      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. NT2AMX) THEN
         CALL QUIT('Insufficient space in CCSDT_GMAT_NODDY')
      ENDIF

      CALL DZERO(WORK(KOME1),NT1AMX)

      ! read T^0 doubles amplitudes from file and square up 
      IOPT   = 2
      KDUM = KEND1
      Call CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KDUM),WORK(KEND1))
      CALL CC_T2SQ(WORK(KEND1),WORK(KT2AM),ISYM0)

      CALL CC3_L3_OMEGA1_NODDY(WORK(KOME1),WORK(KL3AM),
     &                         WORK(KBCIOOOO),WORK(KBCIOVVO),
     &                         WORK(KBCIOOVV),WORK(KBCIVVVV),
     &                         WORK(KT2AM))


      ! read T^C doubles amplitudes from file and square up 
      IOPT   = 2
      Call CC_RDRSP(LISTC,IDLSTC,ISYMC,IOPT,MODEL,
     &              WORK(KDUM),WORK(KEND1))
      CALL CC_T2SQ(WORK(KEND1),WORK(KT2AM),ISYMC)

      CALL CC3_L3_OMEGA1_NODDY(WORK(KOME1),WORK(KL3AM),
     &                         WORK(KBIOOOO),WORK(KBIOVVO),
     &                         WORK(KBIOOVV),WORK(KBIVVVV),
     &                         WORK(KT2AM))


      ! read T^B doubles amplitudes from file and square up 
      IOPT   = 2
      Call CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &              WORK(KDUM),WORK(KEND1))
      CALL CC_T2SQ(WORK(KEND1),WORK(KT2AM),ISYMB)

      CALL CC3_L3_OMEGA1_NODDY(WORK(KOME1),WORK(KL3AM),
     &                         WORK(KCIOOOO),WORK(KCIOVVO),
     &                         WORK(KCIOOVV),WORK(KCIVVVV),
     &                         WORK(KT2AM))

      DO I = 1,NT1AMX
         OMEGA1(I) = OMEGA1(I) + WORK(KOME1+I-1)
      END DO

*---------------------------------------------------------------------*
*     Compute contribution from  <L_3|[H^BC,\tau_nu_2]|HF>
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KOME2),NT1AMX*NT1AMX)

      CALL CC3_L3_OMEGA2_NODDY(WORK(KOME2),WORK(KL3AM),
     *                         WORK(KINT1SBC),WORK(KINT2SBC))

      DO I = 1,NT1AMX
         DO J = 1,I
            IJ = NT1AMX*(I-1) + J
            NIJ = INDEX(I,J)
            OMEGA2(NIJ) = OMEGA2(NIJ) + WORK(KOME2+IJ-1)
         END DO
      END DO

      CALL QEXIT('CCSDT_GMAT_NODDY')

      RETURN
      END 
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_GMAT_NODDY                     *
*---------------------------------------------------------------------*
