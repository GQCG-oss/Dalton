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
      SUBROUTINE CCSDT_FMAT_NODDY(LISTL,IDLSTL,LISTB,IDLSTB,
     &                            IOPTRES,NEW_NODDY_F_CONT,
     &                            OMEGA1,OMEGA2,
     &                            OMEGA1EFF,OMEGA2EFF,
     &                            IDOTS,DOTPROD,LISTDP,ITRAN,
     &                            NFTRAN,MXVEC,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples contribution to F matrix transformation
*
*             (F T^B)^eff_1,2 = (F T^B)_1,2(CCSD) 
*                               + (F T^B)_1,2(L3,T^B3)
*                               - (F T^B)_3 A_3;1,2 (w_3 - w)^1 
*
*        
*     Written by Christof Haettig, April 2002 
*     based on different other noddy codes
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccfield.h>
#include <ccorb.h>
#include <inftap.h>
#include <dummy.h>
#include <ccnoddy.h>

      LOGICAL LOCDBG, FD_TEST, XI_ONLY
      PARAMETER (LOCDBG=.FALSE., FD_TEST=.FALSE., XI_ONLY=.FALSE.)

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      LOGICAL NEW_NODDY_F_CONT
      CHARACTER*3 LISTL, LISTB, LISTDP
      INTEGER LWORK, IDLSTL, IDLSTB, ITRAN, NFTRAN, MXVEC
      INTEGER IDOTS(MXVEC,NFTRAN), IOPTRES
#if defined (SYS_CRAY)
      REAL DOTPROD(MXVEC,NFTRAN), DDOT
      REAL WORK(LWORK), FREQB, FREQL, FREQF, FREQC
      REAL OMEGA1(NT1AMX), OMEGA2(NT2AMX)
      REAL OMEGA1EFF(NT1AMX), OMEGA2EFF(NT2AMX)
      REAL SIXTH, ONE, TWO, TCON, SCON, DCON, FF, SIGN
#else
      DOUBLE PRECISION DOTPROD(MXVEC,NFTRAN), DDOT
      DOUBLE PRECISION WORK(LWORK), FREQB, FREQL, FREQF, FREQC
      DOUBLE PRECISION OMEGA1(NT1AMX), OMEGA2(NT2AMX)
      DOUBLE PRECISION OMEGA1EFF(NT1AMX), OMEGA2EFF(NT2AMX)
      DOUBLE PRECISION SIXTH, ONE, TWO, TCON, SCON, DCON, FF, SIGN
#endif
      PARAMETER(SIXTH=1.0D0/6.0D0, ONE=1.0D0, TWO=2.0D0)

      CHARACTER*10 MODEL
      LOGICAL SKIP_T3_CONT, L2INCL, SKIP_F3
      INTEGER KXINT, KXIAJB, KYIAJB, KT1AMP0, KLAMP0, KLAMH0, KEND1,
     &        LWRK1, KINT1SB, KINT2SB, KEND2, LWRK2, KTB3AM, KL1AM,
     &        KL2AM, KINT1T0, KINT2T0, KINT1TB, KINT2TB, KEND3, LWRK3,
     &        K0IOVVO, K0IOOVV, K0IOOOO, K0IVVVV,
     &        KBIOVVO, KBIOOVV, KBIOOOO, KBIVVVV,
     &        KL3AM, KT02AM, KF3AM, IRECNR, KINT1S0, KINT2S0, KSCR1,
     &        KFOCK0, KFOCKB, KLAMPB, KLAMHB, KOME1, KOME2, KDUM,
     &        LUFOCK, ISYML, KFOCKD, ISYMB, KEND4, LWRK4, 
     &        KTC3AM, KINT1SC, KINT2SC, KLAMPC, KLAMHC, KFOCKC, IDLSTC,
     &        IVEC, ISYMC, KFFD1AM, KFFD2AM, KFFD3AM, KTB1AM, KTB2AM,
     &        ISYMD, ILLL, IDEL, ISYDIS, IJ, NIJ, INDEX, LUTEMP, IDX,
     &        IOPT, ILSTSYM, KTC1AM, KTC2AM, KFIELD, KFLDB1, KFIELDAO,
     &        KT03AM, KT3SCR

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J 

      CALL QENTER('CCSDT_FMAT_NODDY')

      IF (DIRECT) CALL QUIT('CCSDT_FMAT_NODDY: DIRECT NOT IMPLEMENTED')

      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'CCSDT_FMAT_NODDY> F matrix result on entry:'
        CALL CC_PRP(OMEGA1,OMEGA2,1,1,1)
      END IF

      SKIP_T3_CONT = .FALSE.
      IF (NEW_NODDY_F_CONT) THEN
       SKIP_T3_CONT = (IOPTRES.EQ.5)
       IF (LOCDBG .OR. DEBUG) THEN
        WRITE(LUPRI,*) 'Use new noddy F matrix implementation...'
        IF (SKIP_T3_CONT)
     &   WRITE(LUPRI,*) 'Contributions from T3 will be skipped here...'
       END IF
      ENDIF
*---------------------------------------------------------------------*
*     Memory allocation:
*---------------------------------------------------------------------*
      KSCR1   = 1
      KFOCKD  = KSCR1  + NT1AMX
      KEND1   = KFOCKD + NORBT

      KFOCK0  = KEND1
      KEND1   = KFOCK0  + NORBT*NORBT
     
      IF (NONHF) THEN
        KFIELD   = KEND1
        KFLDB1   = KFIELD   + NORBT*NORBT
        KFIELDAO = KFLDB1   + NORBT*NORBT
        KEND1    = KFIELDAO + NORBT*NORBT
      END IF

      KT1AMP0 = KEND1
      KLAMP0  = KT1AMP0 + NT1AMX
      KLAMH0  = KLAMP0  + NLAMDT
      KEND1   = KLAMH0  + NLAMDT

      KTB1AM  = KEND1
      KFOCKB  = KTB1AM + NT1AMX
      KLAMPB  = KFOCKB + NORBT*NORBT
      KLAMHB  = KLAMPB + NLAMDT
      KEND1   = KLAMHB + NLAMDT

      KOME1   = KEND1
      KOME2   = KOME1  + NT1AMX
      KEND1   = KOME2  + NT1AMX*NT1AMX

      KL1AM   = KEND1
      KL2AM   = KL1AM  + NT1AMX
      KEND1   = KL2AM  + NT1AMX*NT1AMX

      KINT1T0 = KEND1
      KINT2T0 = KINT1T0 + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2T0 + NRHFT*NRHFT*NT1AMX

      KINT1S0 = KEND1
      KINT2S0 = KINT1S0 + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2S0 + NRHFT*NRHFT*NT1AMX

      KXIAJB  = KEND1
      KYIAJB  = KXIAJB  + NT1AMX*NT1AMX
      KEND1   = KYIAJB  + NT1AMX*NT1AMX

      KINT1SB = KEND1
      KINT2SB = KINT1SB + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2SB + NRHFT*NRHFT*NT1AMX

      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         write(lupri,*) 'Need:',kend1,' have:',lwork
         CALL QUIT('Insufficient space in CCSDT_FMAT_NODDY (1)')
      ENDIF

      CALL DZERO(WORK(KOME1),NT1AMX)
      CALL DZERO(WORK(KOME2),NT1AMX*NT1AMX)

*---------------------------------------------------------------------*
*     Get zeroth-order Lambda matrices:
*---------------------------------------------------------------------*
      IOPT   = 1
      KDUM = KEND1
      Call CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AMP0),WORK(KDUM))

      Call LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AMP0),
     &            WORK(KEND1),LWRK1)

*---------------------------------------------------------------------*
*     Read precalculated integrals from file:
*           XINT1S0 =  (CK|BD)
*           XINT2S0 =  (CK|LJ)
*           XINT1T0 =  (KC|BD)
*           XINT2T0 =  (KC|LJ)
*           XIAJB   = 2(IA|JB) - (IB|JA)
*           YIAJB   =  (IA|JB)
*---------------------------------------------------------------------*
      CALL CCSDT_READ_NODDY(.TRUE.,WORK(KFOCKD),WORK(KFOCK0),
     &                             WORK(KFIELD),WORK(KFIELDAO),
     &                      .TRUE.,WORK(KXIAJB),WORK(KYIAJB),
     &                      .TRUE.,WORK(KINT1S0),WORK(KINT2S0),
     &                      .TRUE.,WORK(KINT1T0),WORK(KINT2T0),
     &                      .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &                      NORBT,NLAMDT,NRHFT,NVIRT,NT1AMX)

*---------------------------------------------------------------------*
*     If needed, compute zeroth-order triples amplitudes:
*---------------------------------------------------------------------*
      IF (NONHF) THEN
        KT03AM = KEND1
        KEND2  = KT03AM + NT1AMX*NT1AMX*NT1AMX

        LWRK2  = LWORK  - KEND2
        IF (LWRK2 .LT. 0) THEN
          CALL QUIT('Insufficient space in CCSDT_FMAT_NODDY (T03)')
        END IF

        LUTEMP = -1
        CALL GPOPEN(LUTEMP,FILNODT30,'UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        READ(LUTEMP) (WORK(KT03AM+I-1), I=1,NT1AMX*NT1AMX*NT1AMX)
        CALL GPCLOSE(LUTEMP,'KEEP')


        ! store triples amplitudes on file
        LUTEMP = -1
        CALL GPOPEN(LUTEMP,'T3AMPL','UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        REWIND LUTEMP
        WRITE (LUTEMP) (WORK(KT03AM-1+IDX),IDX=1,NT1AMX*NT1AMX*NT1AMX)
        CALL GPCLOSE(LUTEMP,'KEEP')
      END IF

*---------------------------------------------------------------------*
*     Compute T^B_3 amplitudes and some more integrals:
*
*           LISTB = "R1" --> first-order response amplitudes
*                 = "RE" --> right eigenvectors
*
*           XINT1SB =  (CK|BD)-bar
*           XINT2SB =  (CK|LJ)-bar
*---------------------------------------------------------------------*
      KTB3AM = KEND1
      KEND2  = KTB3AM + NT1AMX*NT1AMX*NT1AMX

      ! for finite difference we need T3 until the end
      IF (FD_TEST) KEND1  = KEND2

      LWRK2  = LWORK  - KEND2
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_FMAT_NODDY (2)')
      ENDIF

      IF      (LISTB(1:3).EQ.'R1 ' .OR. LISTB(1:3).EQ.'RE ' .OR.
     &         LISTB(1:3).EQ.'RC '                              ) THEN
        KDUM = KEND2
        CALL CCSDT_T31_NODDY(WORK(KTB3AM),LISTB,IDLSTB,FREQB,XI_ONLY,
     &                       .FALSE.,WORK(KINT1S0),WORK(KINT2S0),
     &                       .FALSE.,WORK(KINT1T0),WORK(KINT2T0),
     &                       .FALSE.,WORK(KXIAJB),WORK(KYIAJB),
     &                               WORK(KINT1SB),WORK(KINT2SB),
     &                       WORK(KLAMPB),WORK(KLAMHB),WORK(KFOCKB),
     &                       WORK(KLAMP0),WORK(KLAMH0),WORK(KFOCK0),
     &                       WORK(KDUM),WORK(KFOCKD),
     &                       WORK(KEND2),LWRK2)
      ELSE IF (LISTB(1:3).EQ.'R2 ' .OR. LISTB(1:3).EQ.'ER1') THEN
        CALL CCSDT_T32_NODDY(WORK(KTB3AM),LISTB,IDLSTB,FREQB,
     &                       WORK(KINT1S0),WORK(KINT2S0),
     &                       WORK(KLAMP0),WORK(KLAMH0),WORK(KFOCK0),
     &                       WORK(KFOCKD),WORK(KFIELDAO),WORK(KFIELD),
     &                       WORK(KSCR1),WORK(KEND2),LWRK2)
        IOPT = 1
        ISYMB = ILSTSYM(LISTB,IDLSTB)
        CALL CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,WORK(KTB1AM),DUMMY)
        CALL CCLR_LAMTRA(WORK(KLAMP0),WORK(KLAMPB),WORK(KLAMH0),
     &                   WORK(KLAMHB),WORK(KTB1AM),ISYMB)
        CALL CCSDT_INTS1_NODDY(.TRUE.,WORK(KINT1SB),WORK(KINT2SB),
     &                         .FALSE.,DUMMY,DUMMY,
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KEND2),LWRK2)
      ELSE
        CALL QUIT('Unknown or illegal list in CCSDT_FMAT_NODDY.')
      END IF

      IF (NONHF .OR. FD_TEST) THEN
        LUTEMP = -1
        CALL GPOPEN(LUTEMP,'T3BARAM','UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        REWIND LUTEMP
        WRITE (LUTEMP) (WORK(KTB3AM-1+IDX),IDX=1,NT1AMX*NT1AMX*NT1AMX)
        CALL GPCLOSE(LUTEMP,'KEEP')
        IF (LOCDBG) WRITE(LUPRI,*) 'Dumped R1_3... norm^2=',
     &    DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KTB3AM),1,WORK(KTB3AM),1)
      END IF

*---------------------------------------------------------------------*
*     Compute contribution from <L_2|[[H,T^B_3],\tau_nu_1|HF>:
*---------------------------------------------------------------------*
      IF (LWRK2 .LT. NT2AMX) THEN
         CALL QUIT('Insufficient space in CCSDT_FMAT_NODDY (3)')
      ENDIF

      ISYML = ILSTSYM(LISTL,IDLSTL)

      ! read left amplitudes from file and square up the doubles part
      IOPT   = 3
      Call CC_RDRSP(LISTL,IDLSTL,ISYML,IOPT,MODEL,
     &              WORK(KL1AM),WORK(KEND2))
      CALL CC_T2SQ(WORK(KEND2),WORK(KL2AM),ISYML)

      IF (.NOT.SKIP_T3_CONT) THEN ! might be done via CCSDT_FBC_NODDY
        CALL T3_LEFT1(WORK(KOME1),WORK(KYIAJB),WORK(KXIAJB),
     &                WORK(KL2AM),WORK(KTB3AM))
 
        CALL T3_LEFT2(WORK(KOME1),WORK(KYIAJB),WORK(KXIAJB),
     &                WORK(KL2AM),WORK(KTB3AM))
 
        CALL T3_LEFT3(WORK(KOME1),WORK(KXIAJB),WORK(KL2AM),WORK(KTB3AM))
      END IF

*---------------------------------------------------------------------*
*     Compute integrals needed for the following contributions:
*---------------------------------------------------------------------*
      KEND2 = KEND1
     
      K0IOVVO = KEND2
      K0IOOVV = K0IOVVO + NRHFT*NVIRT*NVIRT*NRHFT
      K0IOOOO = K0IOOVV + NRHFT*NVIRT*NVIRT*NRHFT
      K0IVVVV = K0IOOOO + NRHFT*NRHFT*NRHFT*NRHFT
      KEND2   = K0IVVVV + NVIRT*NVIRT*NVIRT*NVIRT

      KBIOVVO = KEND2
      KBIOOVV = KBIOVVO + NRHFT*NVIRT*NVIRT*NRHFT
      KBIOOOO = KBIOOVV + NRHFT*NVIRT*NVIRT*NRHFT
      KBIVVVV = KBIOOOO + NRHFT*NRHFT*NRHFT*NRHFT
      KEND2   = KBIVVVV + NVIRT*NVIRT*NVIRT*NVIRT

      KINT1TB = KEND2
      KINT2TB = KINT1TB + NT1AMX*NVIRT*NVIRT
      KEND2   = KINT2TB + NRHFT*NRHFT*NT1AMX

      LWRK2  = LWORK  - KEND2
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
         CALL QUIT('Insufficient space in CCSDT_FMAT_NODDY (4)')
      ENDIF

      CALL CCSDT_READ_NODDY(.FALSE.,WORK(KFOCKD),WORK(KFOCK0),
     &                              WORK(KFIELD),WORK(KFIELDAO),
     &                      .FALSE.,WORK(KXIAJB),WORK(KYIAJB),
     &                      .FALSE.,WORK(KINT1S0),WORK(KINT2S0),
     &                      .FALSE.,WORK(KINT1T0),WORK(KINT2T0),
     &                      .TRUE.,WORK(K0IOVVO),WORK(K0IOOVV),
     &                             WORK(K0IOOOO),WORK(K0IVVVV),
     &                      NORBT,NLAMDT,NRHFT,NVIRT,NT1AMX)

      CALL DZERO(WORK(KBIOVVO),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KBIOOVV),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KBIOOOO),NRHFT*NRHFT*NRHFT*NRHFT)
      CALL DZERO(WORK(KBIVVVV),NVIRT*NVIRT*NVIRT*NVIRT)

      CALL DZERO(WORK(KINT1TB),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2TB),NRHFT*NRHFT*NT1AMX)

      DO ISYMD = 1, NSYM
         DO ILLL = 1,NBAS(ISYMD)
            IDEL   = IBAS(ISYMD) + ILLL
            ISYDIS = MULD2H(ISYMD,ISYMOP)
 
C           ----------------------------
C           Work space allocation no. 2.
C           ----------------------------
            KXINT  = KEND2
            KEND3  = KXINT + NDISAO(ISYDIS)
            LWRK3  = LWORK - KEND3
            IF (LWRK3 .LT. 0) THEN
               WRITE(LUPRI,*) 'Need : ',KEND3,'Available : ',LWORK
               CALL QUIT('Insufficient space in CCSDT_FMAT_NODDY (5)')
            ENDIF
 
C           ---------------------------
C           Read in batch of integrals.
C           ---------------------------
            CALL CCRDAO(WORK(KXINT),IDEL,1,WORK(KEND3),LWRK3,
     *                  IRECNR,DIRECT)
 
C           ----------------------------------
C           Calculate integrals needed in CC3:
C           ----------------------------------

            CALL CCSDT_TRAN1_R(WORK(KINT1TB),WORK(KINT2TB),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KXINT),IDEL)

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

         END DO   
      END DO  

*---------------------------------------------------------------------*
*     Compute L_3 multipliers:
*---------------------------------------------------------------------*
      KEND3  = KEND2

      KL3AM  = KEND3
      KEND3  = KL3AM + NT1AMX*NT1AMX*NT1AMX

      LWRK3  = LWORK  - KEND3
      IF (LWRK3 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_FMAT_NODDY (6)')
      ENDIF


      IF      (LISTL(1:3).EQ.'L0 ') THEN

        FREQL = 0.0D0

        IF (NONHF .AND. LWRK3.LT.NT1AMX*NT1AMX*NT1AMX)
     *    CALL QUIT('Out of memory in CCSDT_F_NODDY.')

        CALL CCSDT_L03AM(WORK(KL3AM),WORK(KINT1T0),WORK(KINT2T0),
     *                   WORK(KXIAJB),WORK(KFOCK0),WORK(KL1AM),
     *                   WORK(KL2AM),WORK(KSCR1),WORK(KFOCKD),
     *                   WORK(KFIELD),WORK(KEND3))

        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KL3AM),1)

      ELSE IF (LISTL(1:3).EQ.'L1 ' .OR. LISTL(1:3).EQ.'LE ' .OR.
     &         LISTL(1:3).EQ.'M1 ' .OR. LISTL(1:3).EQ.'N2 ' .OR.
     &         LISTL(1:3).EQ.'E0 '                              ) THEN

        CALL CCSDT_TBAR31_NODDY(WORK(KL3AM),FREQL,LISTL,IDLSTL,
     &                        WORK(KLAMP0),WORK(KLAMH0),
     &                        WORK(KFOCK0),WORK(KFOCKD),WORK(KSCR1),
     &                        WORK(KXIAJB),WORK(KINT1T0),WORK(KINT2T0),
     &                        WORK(KEND3),LWRK3)

        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KL3AM),1)

      ELSE

        ! FREQL = ??

        CALL QUIT('CCSDT_FMAT_NODDY> LISTL NOT AVAILABLE:'//LISTL)
      
      END IF

*---------------------------------------------------------------------*
*     Compute contribution from  <L_3|[[H^B,T^0_2],\tau_nu_1]|HF>
*                          and   <L_3|[[H^0,T^B_2],\tau_nu_1]|HF>
*     for finite difference include <L_3|[V,T^B_3],\tau_nu_1]|HF>
*---------------------------------------------------------------------*
      KT02AM = KEND3
      KTB2AM = KT02AM + NT1AMX*NT1AMX
      KEND3  = KTB2AM + NT1AMX*NT1AMX
      LWRK3  = LWORK  - KEND3
      IF (LWRK3 .LT. NT2AMX) THEN
         CALL QUIT('Insufficient space in CCSDT_FMAT_NODDY (7)')
      ENDIF

      ! read T^0 doubles amplitudes from file and square up 
      IOPT   = 2
      KDUM = KEND3
      Call CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KDUM),WORK(KEND3))
      CALL CC_T2SQ(WORK(KEND3),WORK(KT02AM),ISYM0)

      CALL CC3_L3_OMEGA1_NODDY(WORK(KOME1),WORK(KL3AM),
     &                         WORK(KBIOOOO),WORK(KBIOVVO),
     &                         WORK(KBIOOVV),WORK(KBIVVVV),
     &                         WORK(KT02AM))



      ! read T^B amplitudes from file and square doubles up 
      ISYMB = ILSTSYM(LISTB,IDLSTB)
      IOPT  = 3
      Call CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &              WORK(KTB1AM),WORK(KEND3))
      Call CCLR_DIASCL(WORK(KEND3),TWO,ISYMB)
      CALL CC_T2SQ(WORK(KEND3),WORK(KTB2AM),ISYMB)

      CALL CC3_L3_OMEGA1_NODDY(WORK(KOME1),WORK(KL3AM),
     &                         WORK(K0IOOOO),WORK(K0IOVVO),
     &                         WORK(K0IOOVV),WORK(K0IVVVV),
     &                         WORK(KTB2AM))


      IF (NONHF) THEN
        KTB3AM  = KEND3
        KEND4   = KTB3AM + NT1AMX*NT1AMX*NT1AMX
        LWRK4   = LWORK  - KEND4
        IF (LWRK4 .LT. 0) THEN
          CALL QUIT('Insufficient space in CCSDT_FMAT_NODDY (nhf1)')
        ENDIF

        LUTEMP = -1
        CALL GPOPEN(LUTEMP,'T3BARAM','UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        REWIND LUTEMP
        READ (LUTEMP) (WORK(KTB3AM-1+IDX),IDX=1,NT1AMX*NT1AMX*NT1AMX)
        CALL GPCLOSE(LUTEMP,'KEEP')

        CALL CCSDT_E1AM(WORK(KOME1),WORK(KL3AM),
     &                  WORK(KTB3AM),WORK(KFIELD))
      END IF


      DO I = 1,NT1AMX
         OMEGA1(I) = OMEGA1(I) + WORK(KOME1+I-1)
      END DO

*---------------------------------------------------------------------*
*     Compute contribution from  <L_3|[H^B,\tau_nu_2]|HF>
*     for finite difference include <L_3|[[V,T^B_2],\tau_nu_2]|HF>
*---------------------------------------------------------------------*

      CALL CC3_L3_OMEGA2_NODDY(WORK(KOME2),WORK(KL3AM),
     *                         WORK(KINT1SB),WORK(KINT2SB))

      IF (NONHF) THEN
         CALL CCSDT_E2AM(WORK(KOME2),WORK(KL3AM),
     *                   WORK(KTB2AM),WORK(KFIELD))
      END IF

      DO I = 1,NT1AMX
         DO J = 1,I
            IJ = NT1AMX*(I-1) + J
            NIJ = INDEX(I,J)
            OMEGA2(NIJ) = OMEGA2(NIJ) + WORK(KOME2+IJ-1)
         END DO
      END DO

      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'CCSDT_FMAT_NODDY> F matrix result on output:'
        CALL CC_PRP(OMEGA1,OMEGA2,1,1,1)
      END IF

*---------------------------------------------------------------------*
*     Compute response Fock matrix F^B(kc) = sum_ia L_kcia T^B_ia
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KFOCKB),NORBT*NORBT)
      
      CALL CCSDT_FCK_R(WORK(KFOCKB),WORK(KXIAJB),WORK(KTB1AM))

*---------------------------------------------------------------------*
*     Compute triples result vector <L_2|[H^B,\tau_nu_3]|HF>:
*---------------------------------------------------------------------*
      KEND3  = KEND2
     
      IF (NONHF) THEN
        ! in case of non-HF fields we need to keep L3
        KEND3  = KL3AM + NT1AMX*NT1AMX*NT1AMX
      END IF

      KF3AM  = KEND3
      KEND3  = KF3AM + NT1AMX*NT1AMX*NT1AMX

      LWRK3  = LWORK  - KEND3
      IF (LWRK3 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_FMAT_NODDY (8)')
      ENDIF
  
      CALL DZERO(WORK(KF3AM),NT1AMX*NT1AMX*NT1AMX)

      IF (.NOT.SKIP_T3_CONT) THEN ! might be done via CCSDT_FBC_NODDY
        CALL CCSDT_F3AM(WORK(KF3AM),WORK(KFOCKB),WORK(KINT1TB),
     &                  WORK(KINT2TB),WORK(KL2AM))
      END IF

      IF (NONHF) THEN
c       -------------------------------------------
c       compute one-index transformed field [V,R1]:
c       -------------------------------------------
        CALL DCOPY(NORBT*NORBT,WORK(KFIELDAO),1,WORK(KFLDB1),1)
        CALL CC_FCKMO(WORK(KFLDB1),WORK(KLAMP0),WORK(KLAMHB),
     *                WORK(KEND3),LWRK3,1,1,1)
        CALL CC_FCKMO(WORK(KFIELDAO),WORK(KLAMPB),WORK(KLAMH0),
     *                WORK(KEND3),LWRK3,1,1,1)
        CALL DAXPY(NORBT*NORBT,ONE,WORK(KFIELDAO),1,WORK(KFLDB1),1)
        IF (LOCDBG) WRITE(LUPRI,*) 'Norm^2(FLDB1):',
     *    DDOT(NORBT*NORBT,WORK(KFLDB1),1,WORK(KFLDB1),1)

        L2INCL = .FALSE.
        KDUM = KEND3
        CALL CCSDT_E3AM(WORK(KF3AM),WORK(KDUM),WORK(KL3AM),
     &                  WORK(KFLDB1),L2INCL)
      END IF

*---------------------------------------------------------------------*
*     Compute finite difference result and compare:
*---------------------------------------------------------------------*
      IF (FD_TEST) THEN
        KTB2AM  = KEND3
        KTB3AM  = KTB2AM  + NT1AMX*NT1AMX
        KFFD1AM = KTB3AM  + NT1AMX*NT1AMX*NT1AMX
        KFFD2AM = KFFD1AM + NT1AMX
        KFFD3AM = KFFD2AM + NT1AMX*NT1AMX
        KEND4   = KFFD3AM + NT1AMX*NT1AMX*NT1AMX
        LWRK4   = LWORK  - KEND4
        IF (LWRK4 .LT. NT2AMX) THEN
           CALL QUIT('Insufficient space in CCSDT_FMAT_NODDY (9)')
        ENDIF

        LUTEMP = -1
        CALL GPOPEN(LUTEMP,'T3BARAM','UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        REWIND LUTEMP
        READ (LUTEMP) (WORK(KTB3AM-1+IDX),IDX=1,NT1AMX*NT1AMX*NT1AMX)
        CALL GPCLOSE(LUTEMP,'DELETE')

        ! read T^B doubles amplitudes from file and square up 
        ISYMB = ILSTSYM(LISTB,IDLSTB)
        IOPT  = 3
        Call CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &                WORK(KTB1AM),WORK(KEND4))
        Call CCLR_DIASCL(WORK(KEND4),TWO,ISYMB)
        CALL CC_T2SQ(WORK(KEND4),WORK(KTB2AM),ISYMB)
 
        CALL CCSDT_FMAT_FD(WORK(KTB1AM),WORK(KTB2AM),WORK(KTB3AM),
     &                     WORK(KFFD1AM),WORK(KFFD2AM),WORK(KFFD3AM),
     &                     WORK(KFOCKB),WORK(KFOCK0),WORK(KEND4),LWRK4)

 
        WRITE(LUPRI,*) 'CCSDT_F_NODDY> my F1AM:'
        CALL OUTPUT(WORK(KOME1),1,NVIRT,1,NRHFT,
     &              NVIRT,NRHFT,1,LUPRI)
        WRITE(LUPRI,*) 'CCSDT_F_NODDY> finite difference F1AM:'
        CALL OUTPUT(WORK(KFFD1AM),1,NVIRT,1,NRHFT,
     &              NVIRT,NRHFT,1,LUPRI)
        CALL DAXPY(NT1AMX,-1.0D0,WORK(KOME1),1,WORK(KFFD1AM),1)
        WRITE(LUPRI,*) 'CCSDT_F_NODDY> norm of difference:',
     &   DDOT(NT1AMX,WORK(KFFD1AM),1,WORK(KFFD1AM),1)

        WRITE(LUPRI,*) 'CCSDT_F_NODDY> my F2AM:'
        CALL OUTPUT(WORK(KOME2),1,NT1AMX,1,NT1AMX,
     &              NT1AMX,NT1AMX,1,LUPRI)
        WRITE(LUPRI,*) 'CCSDT_F_NODDY> finite difference F2AM:'
        CALL OUTPUT(WORK(KFFD2AM),1,NT1AMX,1,NT1AMX,
     &              NT1AMX,NT1AMX,1,LUPRI)
        CALL DAXPY(NT1AMX*NT1AMX,-1.0D0,WORK(KOME2),1,WORK(KFFD2AM),1)
        WRITE(LUPRI,*) 'CCSDT_F_NODDY> norm of difference:',
     &   DDOT(NT1AMX*NT1AMX,WORK(KFFD2AM),1,WORK(KFFD2AM),1)

        WRITE(LUPRI,*) 'CCSDT_F_NODDY> my F3AM:'
        CALL OUTPUT(WORK(KF3AM),1,NT1AMX*NT1AMX,1,NT1AMX,
     &              NT1AMX*NT1AMX,NT1AMX,1,LUPRI)
        WRITE(LUPRI,*) 'CCSDT_F_NODDY> finite difference F3AM:'
        CALL OUTPUT(WORK(KFFD3AM),1,NT1AMX*NT1AMX,1,NT1AMX,
     &              NT1AMX*NT1AMX,NT1AMX,1,LUPRI)
        CALL DAXPY(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KF3AM),1,
     &                 WORK(KFFD3AM),1)
        WRITE(LUPRI,*) 'CCSDT_F_NODDY> norm of difference:',
     &   DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KFFD3AM),1,WORK(KFFD3AM),1)

        WRITE(LUPRI,*) 'CCSDT_F_NODDY> difference vector:'
        CALL OUTPUT(WORK(KFFD3AM),1,NT1AMX*NT1AMX,1,NT1AMX,
     &              NT1AMX*NT1AMX,NT1AMX,1,LUPRI)
 
      END IF

*---------------------------------------------------------------------*
*     Now we split:
*       for IOPTRES < 5 we compute the effective F transformed vector
*       for IOPTRES = 5 we compute the contractions F T^B T^C
*---------------------------------------------------------------------*
      IF (IOPTRES.GE.1 .AND. IOPTRES.LE.4) THEN

        KT02AM = KEND3
        KEND3  = KT02AM + NT1AMX*NT1AMX
        LWRK3  = LWORK  - KEND3
        IF (LWRK3 .LT. NT2AMX) THEN
           CALL QUIT('Insufficient space in CCSDT_FMAT_NODDY (10)')
        ENDIF

        ! read T^0 doubles amplitudes from file and square up 
        IOPT   = 2
        KDUM = KEND3
        Call CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KDUM),WORK(KEND3))
        CALL CC_T2SQ(WORK(KEND3),WORK(KT02AM),ISYM0)

        CALL DCOPY(NT1AMX,OMEGA1,1,OMEGA1EFF,1)
        CALL DCOPY(NT2AMX,OMEGA2,1,OMEGA2EFF,1)

        FREQF = FREQL + FREQB 

        CALL CC_LHPART_NODDY(OMEGA1EFF,OMEGA2EFF,WORK(KF3AM),-FREQF,
     &                       WORK(KFOCKD),WORK(KFIELD),
     &                       WORK(K0IOOOO),WORK(K0IOVVO),
     &                       WORK(K0IOOVV),WORK(K0IVVVV),
     &                       WORK(KT02AM),WORK(KINT1S0),WORK(KINT2S0),
     &                       WORK(KEND3),LWRK3)


      ELSE IF (IOPTRES.EQ.5) THEN

        SIGN = +1.0D0

        SKIP_F3 = SKIP_T3_CONT .AND. (.NOT.NONHF)

        CALL CCDOTRSP_NODDY(WORK(KOME1),WORK(KOME2),WORK(KF3AM),SIGN,
     &                      ITRAN,LISTDP,IDOTS,DOTPROD,MXVEC,
     &                      WORK(KLAMP0),WORK(KLAMH0),
     &                      WORK(KFOCK0),WORK(KFOCKD),
     &                      WORK(KXIAJB), WORK(KYIAJB),
     &                      WORK(KINT1T0),WORK(KINT2T0),
     &                      WORK(KINT1S0),WORK(KINT2S0),
     &                      'CCSDT_F_NODDY',LOCDBG,LOCDBG,SKIP_F3,
     &                      WORK(KEND3),LWRK3)

      ELSE
        CALL QUIT('Illegal value for IOPTRES IN CCSDT_FMAT_NODDY')
      END IF

      CALL QEXIT('CCSDT_FMAT_NODDY')

      RETURN
      END 
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_FMAT_NODDY                     *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_F3AM(F3AM,FOCK,XINT1,XINT2,L2AM)  
*---------------------------------------------------------------------*
*
*    Purpose: compute triples exictation amplitudes of
*             F matrix transformed vector
*
*             (F T^B)_nu_3 = <L_2|[H^B,tau_nu_3]|HF>
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

#if defined (SYS_CRAY)
      REAL AIBJCK
      REAL F3AM(NT1AMX,NT1AMX,NT1AMX)
      REAL L2AM(NT1AMX,NT1AMX)
      REAL XINT1(NT1AMX,NVIRT,NVIRT)
      REAL XINT2(NT1AMX,NRHFT,NRHFT)
      REAL FOCK(NORBT,NORBT)
      REAL TWO
#else
      DOUBLE PRECISION AIBJCK
      DOUBLE PRECISION F3AM(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION L2AM(NT1AMX,NT1AMX)
      DOUBLE PRECISION XINT1(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT2(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION FOCK(NORBT,NORBT)
      DOUBLE PRECISION TWO
#endif
      PARAMETER ( TWO = 2.0D0 )

      INTEGER NI, NA, NAI, NK, NC, NCK, NJ, NCJ, NB, NBJ, NBK, NBI,
     &        NAJ, ND, NDJ, NDK, NL, NAL, NBL, NCI, NAK
C
      DO NK = 1,NRHFT
       DO NC = 1,NVIRT
        NCK = NVIRT*(NK-1) + NC
C       
        DO NJ = 1,NRHFT
         NCJ = NVIRT*(NJ-1) + NC

         DO NB = 1,NVIRT
          NBJ = NVIRT*(NJ-1) + NB
          NBK = NVIRT*(NK-1) + NB
C       
          DO NI = 1,NRHFT
           NBI = NVIRT*(NI-1) + NB
        
           DO NA = 1,NVIRT
             NAI = NVIRT*(NI-1) + NA
             NAJ = NVIRT*(NJ-1) + NA
C       
             AIBJCK = 0.0D0

             AIBJCK = AIBJCK + FOCK(NK,NRHFT+NC) * L2AM(NAI,NBJ)
C
             AIBJCK = AIBJCK - FOCK(NJ,NRHFT+NC) * L2AM(NAI,NBK)
C       
             DO ND = 1,NVIRT
                NDJ = NVIRT*(NJ-1) + ND
                NDK = NVIRT*(NK-1) + ND
C       
                AIBJCK = AIBJCK + 
     &                   (TWO*XINT1(NCK,ND,NB)-XINT1(NBK,ND,NC))*
     &                   L2AM(NDJ,NAI)
C       
                AIBJCK = AIBJCK - XINT1(NBI,ND,NC) * L2AM(NDK,NAJ)
             END DO
C       
             DO NL = 1,NRHFT
                NAL = NVIRT*(NL-1) + NA
                NBL = NVIRT*(NL-1) + NB
C       
                AIBJCK = AIBJCK - 
     &                   (TWO*XINT2(NCK,NJ,NL)-XINT2(NCJ,NK,NL))*
     &                   L2AM(NBL,NAI)
C       
                AIBJCK = AIBJCK + XINT2(NCJ,NI,NL) * L2AM(NBK,NAL)
             END DO
C       
             F3AM(NAI,NBJ,NCK) = F3AM(NAI,NBJ,NCK) + AIBJCK
             F3AM(NAI,NCK,NBJ) = F3AM(NAI,NCK,NBJ) + AIBJCK
             F3AM(NBJ,NAI,NCK) = F3AM(NBJ,NAI,NCK) + AIBJCK
             F3AM(NCK,NAI,NBJ) = F3AM(NCK,NAI,NBJ) + AIBJCK
             F3AM(NBJ,NCK,NAI) = F3AM(NBJ,NCK,NAI) + AIBJCK
             F3AM(NCK,NBJ,NAI) = F3AM(NCK,NBJ,NAI) + AIBJCK
C
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

      CALL CCSDT_CLEAN_T3(F3AM,NT1AMX,NVIRT,NRHFT)

      RETURN
      END 
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_F3AM                           *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_FMAT_FD(TB1AM,TB2AM,TB3AM,F1AM,F2AM,F3AM,
     &                         FOCKA,FOCK0,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: Construct F matrix transformed vector 
*             by finite difference on left hand transformation
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
      REAL DELTA, ZERO
#else
      DOUBLE PRECISION DELTA, ZERO, ONE
#endif
      PARAMETER( DELTA = 1.0D-6, ONE = 1.0D0, ZERO = 0.0D0)

      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL WORK(LWORK), 
     &  FOCK0(NORBT,NORBT), FOCKA(NORBT,NORBT), DDOT,
     &  TB1AM(NT1AMX),TB2AM(NT1AMX,NT1AMX),TB3AM(NT1AMX,NT1AMX,NT1AMX),
     &  F1AM(NT1AMX), F2AM(NT1AMX,NT1AMX), F3AM(NT1AMX,NT1AMX,NT1AMX)
#else
      DOUBLE PRECISION WORK(LWORK), 
     &  FOCK0(NORBT,NORBT), FOCKA(NORBT,NORBT), DDOT,
     &  TB1AM(NT1AMX),TB2AM(NT1AMX,NT1AMX),TB3AM(NT1AMX,NT1AMX,NT1AMX),
     &  F1AM(NT1AMX), F2AM(NT1AMX,NT1AMX), F3AM(NT1AMX,NT1AMX,NT1AMX)
#endif
  
      CHARACTER*10 MODEL
      INTEGER KT1AM, KT2AM, KT3AM, KL1AM, KL2AM, KL3AM, KFOCK,
     &        KLAMDP, KLAMDH, KEND1, LWRK1, LUSIFC, IOPT, IDUMMY

      CALL QUIT('CCSDT_FMAT_FD not adapted for intermediates!')
*---------------------------------------------------------------------*
*     memory allocations:
*---------------------------------------------------------------------*
      KT1AM = 1
      KT2AM = KT1AM + NT1AMX
      KT3AM = KT2AM + NT1AMX*NT1AMX
      KL1AM = KT3AM + NT1AMX*NT1AMX*NT1AMX
      KL2AM = KL1AM + NT1AMX
      KL3AM = KL2AM + NT1AMX*NT1AMX
      KFOCK = KL3AM + NT1AMX*NT1AMX*NT1AMX
      KLAMDP= KFOCK + NORBT*NORBT
      KLAMDH= KLAMDP+ NLAMDT
      KEND1 = KLAMDH+ NLAMDT
      LWRK1 = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_FD')
      ENDIF

      WRITE(LUPRI,*) 'CCSDT_FMAT_FD> norm^2(TB3AM):',
     &   DDOT(NT1AMX**3,TB3AM,1,TB3AM,1)
 
*---------------------------------------------------------------------*
*     Get T0 singles & doubles  and L0 singles & doubles and
*     zeroth-order lambda matrices:
*---------------------------------------------------------------------*
      IOPT   = 3
      Call CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),WORK(KT3AM))
      CALL CC_T2SQ(WORK(KT3AM),WORK(KT2AM),1)

      IOPT   = 3
      Call CC_RDRSP('L0',0,1,IOPT,MODEL,WORK(KL1AM),WORK(KL3AM))
      CALL CC_T2SQ(WORK(KL3AM),WORK(KL2AM),1)

      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     &            WORK(KEND1),LWRK1)

      IOPT = 1
      CALL CC_LHTR_NODDY(0.0D0,F1AM,F2AM,WORK(KT1AM),WORK(KT2AM),
     &                   WORK(KL1AM),WORK(KL2AM),
     &                   WORK(KLAMDP),WORK(KLAMDH),
     &                   WORK(KEND1),LWRK1,IOPT)

*---------------------------------------------------------------------*
*     Get T0 singles, doubles & triples
*---------------------------------------------------------------------*
      IOPT   = 3
      Call CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),WORK(KT3AM))
      CALL CC_T2SQ(WORK(KT3AM),WORK(KT2AM),1)

      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'T3AMPL','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
      READ (LUSIFC) (WORK(KT3AM-1+I), I=1,NT1AMX*NT1AMX*NT1AMX)
      CALL GPCLOSE(LUSIFC,'KEEP')

      WRITE(LUPRI,*) 'CCSDT_FMAT_FD> norm^2(T3AM):',
     &   DDOT(NT1AMX**3,WORK(KT3AM),1,WORK(KT3AM),1)
 
*---------------------------------------------------------------------*
*     Get L0 singles, doubles & triples
*---------------------------------------------------------------------*
      IOPT   = 3
      Call CC_RDRSP('L0',0,1,IOPT,MODEL,WORK(KL1AM),WORK(KL3AM))
      CALL CC_T2SQ(WORK(KL3AM),WORK(KL2AM),1)

      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'L3AMPL','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
      READ (LUSIFC) (WORK(KL3AM-1+I), I=1,NT1AMX*NT1AMX*NT1AMX)
      CALL GPCLOSE(LUSIFC,'KEEP')
 
      WRITE(LUPRI,*) 'CCSDT_FMAT_FD> norm^2(L3AM):',
     &   DDOT(NT1AMX**3,WORK(KL3AM),1,WORK(KL3AM),1)

*---------------------------------------------------------------------*
*     Compute finite difference vector:
*---------------------------------------------------------------------*
      CALL DZERO(F1AM,NT1AMX              )
      CALL DZERO(F2AM,NT1AMX*NT1AMX       )
      CALL DZERO(F3AM,NT1AMX*NT1AMX*NT1AMX)

      CALL DCOPY(NORBT*NORBT,FOCK0,1,WORK(KFOCK),1)

      CALL DAXPY(NT1AMX,              -DELTA,TB1AM,1,WORK(KT1AM),1)
      CALL DAXPY(NT1AMX*NT1AMX,       -DELTA,TB2AM,1,WORK(KT2AM),1)
      CALL DAXPY(NT1AMX*NT1AMX*NT1AMX,-DELTA,TB3AM,1,WORK(KT3AM),1)

      CALL DAXPY(NORBT*NORBT,-DELTA,FOCKA,1,WORK(KFOCK),1)

      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     &            WORK(KEND1),LWRK1)


      CALL CC_LHTR_NODDY2(F1AM,F2AM,F3AM,
     &                    WORK(KT1AM),WORK(KT2AM),WORK(KT3AM),
     &                    WORK(KL1AM),WORK(KL2AM),WORK(KL3AM),
     &                    WORK(KLAMDP),WORK(KLAMDH),
     &                    WORK(KFOCK),WORK(KEND1),LWRK1)

      CALL DSCAL(NT1AMX,              -ONE,F1AM,1)
      CALL DSCAL(NT1AMX*NT1AMX,       -ONE,F2AM,1)
      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-ONE,F3AM,1)

      CALL DAXPY(NT1AMX,              2.D0*DELTA,TB1AM,1,WORK(KT1AM),1)
      CALL DAXPY(NT1AMX*NT1AMX,       2.D0*DELTA,TB2AM,1,WORK(KT2AM),1)
      CALL DAXPY(NT1AMX*NT1AMX*NT1AMX,2.D0*DELTA,TB3AM,1,WORK(KT3AM),1)

      CALL DAXPY(NORBT*NORBT,2.D0*DELTA,FOCKA,1,WORK(KFOCK),1)

      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     &            WORK(KEND1),LWRK1)

      CALL CC_LHTR_NODDY2(F1AM,F2AM,F3AM,
     &                    WORK(KT1AM),WORK(KT2AM),WORK(KT3AM),
     &                    WORK(KL1AM),WORK(KL2AM),WORK(KL3AM),
     &                    WORK(KLAMDP),WORK(KLAMDH),
     &                    WORK(KFOCK),WORK(KEND1),LWRK1)


      CALL DSCAL(NT1AMX              ,0.5D0/DELTA,F1AM,1)
      CALL DSCAL(NT1AMX*NT1AMX       ,0.5D0/DELTA,F2AM,1)
      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,0.5D0/DELTA,F3AM,1)

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_F_FD                           *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_TBAR31_NODDY(TBAR3AM,FREQB,LISTB,IDLSTB,
     &                              XLAMDP,XLAMDH,FOCK0,FOCKD,SCR1,
     &                              XIAJB,XINT1T0,XINT2T0,
     &                              WORK,LWORK)
*---------------------------------------------------------------------*
*
*  Purpose: compute triples part of first-order response multipliers
*
*   Implemented for L1, X1B, E0, N2, M1 and LE vectors.
*
*   Input:   LISTB,   IDLSTB
*            XLAMDP,  XLAMDH
*            FOCK0,   FOCKD
*            XIAJB
*            XINT1T0, XINT2T0
*            
*   Output:  TBAR3AM
*            FREQB
*
*   Scratch: SCR1
*
*  Written by Christof Haettig, April 2002
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccl1rsp.h>
#include <ccn2rsp.h>
#include <cclrmrsp.h>
#include <ccexgr.h>
#include <ccexci.h>
#include <ccfield.h>

      LOGICAL LOCDBG, PRINT_T3
      PARAMETER (LOCDBG=.false.,PRINT_T3=.false.)


      CHARACTER LISTB*3
      LOGICAL RHS_ONLY
      INTEGER LWORK, IDLSTB
#if defined (SYS_CRAY)
      REAL WORK(LWORK), FOCKD(*), FOCK0(*), SCR1(*)
      REAL FREQB, DUMMY, DDOT, ZERO, FF, ONE
      REAL XLAMDP(*), XLAMDH(*)
      REAL TBAR3AM(NT1AMX,NT1AMX,NT1AMX)
      REAL XINT1T0(*), XINT2T0(*), XIAJB(*), YIAJB(*)
#else
      DOUBLE PRECISION WORK(LWORK), FOCKD(*), FOCK0(*), SCR1(*)
      DOUBLE PRECISION FREQB, DUMMY, DDOT, ZERO, FF, ONE
      DOUBLE PRECISION XLAMDP(*), XLAMDH(*)
      DOUBLE PRECISION TBAR3AM(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION XINT1T0(*), XINT2T0(*), XIAJB(*)
#endif
      LOGICAL L2INCL
      INTEGER ISYM0
      PARAMETER( L2INCL=.TRUE. , ISYM0=1, ZERO=0.0D0, ONE=1.0D0)

      CHARACTER MODEL*10, LABELB*8, LISTR*3, LISTL*3
      LOGICAL LORXB, DO_INTR
      INTEGER ISYMB, KFOCKB, KL1AM, KL2AM, KL3AM, KEND1, LWRK1,
     &        IOPT, IRREP, IERR, NAI, NBJ, NCK, NAIBJCK, ILSTNR,
     &        KTB1AM, KLAMPB, KLAMHB, KINT1TB, KINT2TB, ISYMD, ILLL,
     &        IDEL, ISYDIS, KXINT, KEND2, LWRK2, IRECNR, IR1TAMP,
     &        ILSTSYM, ILSTNL, ISYML, ISYMR, IMAT, KFIELD, KFIELDAO,
     &        KFLDB1, KT3SCR
      

      CALL QENTER('CCTBAR31')
*---------------------------------------------------------------------*
*     Do some initial tests and memory allocations:
*---------------------------------------------------------------------*
      if (locdbg) then
        write(lupri,*) 'entered ccsdt_tbar31_noddy for "',LISTB(1:3),
     &    '" ',idlstb
      end if

      if (NSYM.NE.1) CALL QUIT('NO SYMMETRY FOR CCSDT_TBAR31_NODDY')

      ISYMB = ILSTSYM(LISTB,IDLSTB)

c     ------------------------------------------------------------
c     first-order response of ground state lagrangian multipliers:
c     ------------------------------------------------------------
      IF (LISTB(1:3).EQ.'L1 '.OR.LISTB(1:3).EQ.'X1B') THEN
        RHS_ONLY = LISTB.EQ.'X1B'

        ! get symmetry, frequency and integral label from common blocks
        ! defined in ccl1rsp.h
        FREQB  = FRQLRZ(IDLSTB)
        LABELB = LRZLBL(IDLSTB)
        LORXB  = LORXLRZ(IDLSTB)

        IF (LORXB) CALL QUIT('NO ORBITAL RELAX. IN CCSDT_TBAR31_NODDY')

        ! find corresponding first-order right amplitudes
        LISTR  = 'R1 '
        ILSTNR = IR1TAMP(LABELB,LORXB,FREQB,ISYMB)

        ! set zero-order left vector
        LISTL  = 'L0 '
        ILSTNL = 0

c     -------------------------------------------------------
c     zero-order lagrangian multipliers for ground to excited 
c     state transitions:
c     -------------------------------------------------------
      ELSE IF (LISTB(1:3).EQ.'M1 ') THEN
        RHS_ONLY = .FALSE.

        FREQB  = FRQLRM(IDLSTB)
        LABELB = '- none -'

        ! find corresponding right eigenvector
        LISTR  = 'RE '
        ILSTNR = ILRM(IDLSTB)

        ! set zero-order left vector
        LISTL  = 'L0 '
        ILSTNL = 0

c     --------------------------------------------------------
c     zero-order lagrangian multipliers for excited states:
c     (Diagonal case for N2, but with Tbar^(0) added. The 
c      latter does here only contribute indirectly via the
c      singles and doubles part of the solution vector "E0")
c     --------------------------------------------------------
      ELSE IF (LISTB(1:3).EQ.'E0 ') THEN
        RHS_ONLY = .FALSE.

        FREQB  = ZERO
        LABELB = '- none -'

        ! find corresponding right eigenvector
        LISTR  = 'RE '
        ILSTNR = IXGRST(IDLSTB)

        ! set zero-order left vector
        LISTL  = 'LE '
        ILSTNL = IXGRST(IDLSTB)

c     --------------------------------------------------------
c     zero-order lagrangian multipliers for excited to excited
c     state transitions:
c     --------------------------------------------------------
      ELSE IF (LISTB(1:3).EQ.'N2 ') THEN
        RHS_ONLY = .FALSE.

        FREQB  = FRQIN2(IDLSTB) + FRQFN2(IDLSTB)
        LABELB = '- none -'

        ! find corresponding right eigenvector
        LISTR  = 'RE '
        ILSTNR = IFN2(IDLSTB)

        ! set zero-order left vector
        LISTL  = 'LE '
        ILSTNL = IIN2(IDLSTB)

c     -----------------------------
c     zero-order left eigenvectors:
c     -----------------------------
      ELSE IF (LISTB(1:3).EQ.'LE ') THEN
        RHS_ONLY = .FALSE.

        FREQB = -EIGVAL(IDLSTB)
        LABELB = '- none -'

        ! we don't have any "right" vector entering a right hand side
        LISTR  = '---'
        ILSTNR = -99

        ! we don't have any "zero-order left" vector entering a rh side
        LISTL  = '---'
        ILSTNL = -99

      ELSE
        WRITE(LUPRI,*) 'LISTB,IDLSTB:',LISTB(1:3),IDLSTB
        CALL QUIT('Unknown/Illegal list in CCSDT_TBAR31_NODDY.')
      END IF

      IF (LISTB(1:3).NE.'LE ') THEN
        DO_INTR = .TRUE.
        ISYML   = ILSTSYM(LISTL,ILSTNL)
        ISYMR   = ILSTSYM(LISTR,ILSTNR)
        IF (MULD2H(ISYML,ISYMR).NE.ISYMB) CALL QUIT(
     &       'Symmetry mismatch in CCSDT_TBAR31_NODDY.')
      ELSE
        DO_INTR = .FALSE.
      END IF


      KL1AM  = 1
      KL2AM  = KL1AM  + NT1AMX
      KL3AM  = KL2AM  + NT1AMX*NT1AMX
      KEND1  = KL3AM  + NT1AMX*NT1AMX*NT1AMX

      KTB1AM  = KEND1
      KFOCKB  = KTB1AM + NT1AMX
      KLAMPB  = KFOCKB + NORBT*NORBT
      KLAMHB  = KLAMPB + NLAMDT
      KEND1   = KLAMHB + NLAMDT

      IF (NONHF) THEN
        KFIELD   = KEND1
        KFLDB1   = KFIELD   + NORBT*NORBT
        KFIELDAO = KFLDB1   + NORBT*NORBT
        KEND1    = KFIELDAO + NORBT*NORBT
      END IF

      KINT1TB = KEND1
      KINT2TB = KINT1TB + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2TB + NRHFT*NRHFT*NT1AMX

      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. NT2AMX) THEN
         CALL QUIT('Insufficient space in CCSDT_TBAR31_NODDY')
      ENDIF
      
*---------------------------------------------------------------------*
*     Calculate response Lambda matrices from the right vectors:
*---------------------------------------------------------------------*
      IF (DO_INTR) THEN
        IOPT = 1 
        CALL CC_RDRSP(LISTR,ILSTNR,ISYMB,IOPT,MODEL,WORK(KTB1AM),WORK)

        CALL CCLR_LAMTRA(XLAMDP,WORK(KLAMPB),
     &                   XLAMDH,WORK(KLAMHB),WORK(KTB1AM),ISYMB)
      END IF

*---------------------------------------------------------------------*
*     Get the one electron integrals from the fields
*---------------------------------------------------------------------*
      IF (NONHF) THEN
         CALL DZERO(WORK(KFIELDAO),NORBT*NORBT)
         DO I = 1, NFIELD
            FF = EFIELD(I)
            CALL CC_ONEP(WORK(KFIELDAO),WORK(KEND1),LWRK1,
     *                   FF,1,LFIELD(I))
         ENDDO
         CALL DCOPY(NORBT*NORBT,WORK(KFIELDAO),1,WORK(KFIELD),1)

         ! calculate external field in zero-order lambda basis
         CALL CC_FCKMO(WORK(KFIELD),XLAMDP,XLAMDH,
     *                 WORK(KEND1),LWRK1,1,1,1)

         ! one-index transformed field integrals for [V,R1]
         CALL DCOPY(NORBT*NORBT,WORK(KFIELDAO),1,WORK(KFLDB1),1)
         CALL CC_FCKMO(WORK(KFLDB1),XLAMDP,WORK(KLAMHB),
     *                 WORK(KEND1),LWRK1,1,1,1)

         CALL CC_FCKMO(WORK(KFIELDAO),WORK(KLAMPB),XLAMDH,
     *                 WORK(KEND1),LWRK1,1,1,1)
         CALL DAXPY(NORBT*NORBT,ONE,WORK(KFIELDAO),1,WORK(KFLDB1),1)
      ENDIF

*---------------------------------------------------------------------*
*     Compute integrals needed for the following contributions:
*---------------------------------------------------------------------*
      IF (DO_INTR) THEN

        CALL DZERO(WORK(KINT1TB),NT1AMX*NVIRT*NVIRT)
        CALL DZERO(WORK(KINT2TB),NRHFT*NRHFT*NT1AMX)

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
               CALL QUIT('Insufficient space in CCSDT_TBAR31_NODDY')
            ENDIF
 
C           ---------------------------
C           Read in batch of integrals.
C           ---------------------------
            CALL CCRDAO(WORK(KXINT),IDEL,1,WORK(KEND2),LWRK2,
     &                  IRECNR,DIRECT)
 
C           ----------------------------------
C           Calculate integrals needed in CC3:
C           ----------------------------------
            CALL CCSDT_TRAN1_R(WORK(KINT1TB),WORK(KINT2TB),
     &                         XLAMDP,XLAMDH,
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KXINT),IDEL)
          END DO   
        END DO  
      END IF

*---------------------------------------------------------------------*
*     triples part of the zero-order lagrangian multipliers:
*     (note: in case of non-HF external fields TBAR3AM is used
*            here as scratch array)
*---------------------------------------------------------------------*
      IF ( LISTB(1:3).EQ.'L1 ' .OR. LISTB(1:3).EQ.'X1B' .OR.
     &    (LISTB(1:3).NE.'LE '.AND.NONHF) ) THEN

C       IF (.NOT.( LISTB(1:3).EQ.'L1 ' .OR. LISTB(1:3).EQ.'X1B' ))
C    &    CALL QUIT('CCSDT_TBAR31_NODDY not tested for this case!')

        IF (LWRK1 .LT. NT2AMX) 
     &    CALL QUIT('Insufficient space in CCSDT_TBAR31_NODDY')

        ! read L^0 multipliers from file and square up doubles part
        IOPT  = 3
        Call CC_RDRSP('L0 ',0,ISYM0,IOPT,MODEL,
     &                WORK(KL1AM),WORK(KEND1))
        CALL CC_T2SQ(WORK(KEND1),WORK(KL2AM),ISYM0)   

        ! remember that CCSDT_L03AM returns -L3 !!
        CALL CCSDT_L03AM(WORK(KL3AM),XINT1T0,XINT2T0,XIAJB,FOCK0,
     &                  WORK(KL1AM),WORK(KL2AM),SCR1,FOCKD,
     &                  WORK(KFIELD),TBAR3AM)

        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KL3AM),1)
      END IF

*---------------------------------------------------------------------*
*     Initialize result vector TBAR3AM:
*---------------------------------------------------------------------*
      CALL DZERO(TBAR3AM,NT1AMX*NT1AMX*NT1AMX)

*=====================================================================*
*     Compute Eta_3 contribution to Tbar^B_3 (only for "L1" vectors):
*     for this contribution we need to compute first the triples part
*     of the ground state zero-order lagrangian multipliers L^0_3
*     and one-electron property integral matrix:
*=====================================================================*
      IF (LISTB(1:3).EQ.'L1 '.OR.LISTB(1:3).EQ.'X1B') THEN

c       -------------------------------------------
c       Matrix with property integrals in MO basis:
c       -------------------------------------------
        ! read property integrals from file:
        CALL CCPRPAO(LABELB,.TRUE.,WORK(KFOCKB),IRREP,IMAT,IERR,
     &               WORK(KEND1),LWRK1)
        IF ((IERR.GT.0) .OR. (IERR.EQ.0 .AND. IRREP.NE.ISYMB)) THEN
          WRITE(LUPRI,*) 'IERR,IRREP,ISYMB:',IERR,IRREP,ISYMB
          CALL QUIT('CCSDT_TBAR31_NODDY: error reading oper. '//LABELB)
        ELSE IF (IERR.LT.0) THEN
          CALL DZERO(WORK(KFOCKB),N2BST(ISYMB))
        END IF
 
        ! transform property integrals to Lambda-MO basis
        CALL CC_FCKMO(WORK(KFOCKB),XLAMDP,XLAMDH,
     &                WORK(KEND1),LWRK1,ISYMB,1,1)

c       -------------------------------
c       Eta_3 contribution to Tbar^B_3:
c       -------------------------------
        CALL CCSDT_E3AM(TBAR3AM,WORK(KL2AM),WORK(KL3AM),
     &                  WORK(KFOCKB),L2INCL)

        IF (LOCDBG) THEN
          WRITE(LUPRI,*) 'CCSDT_TBAR31_NODDY> norm^2(E3AM):',
     &        DDOT(NT1AMX*NT1AMX*NT1AMX,TBAR3AM,1,TBAR3AM,1)
        END IF

      END IF

*=====================================================================*
*     Compute F matrix (FT^B)_3 contribution to Tbar^B_3 (for all 
*     vectors, apart from the eigenvectors "LE"):
*     (Note, that here we overwrite FOCKB with a new matrix!!!)
*=====================================================================*
      IF (LISTB(1:3).NE.'LE ') THEN
        ! read zero-order left vector from file and square up 
        IOPT  = 3
        Call CC_RDRSP(LISTL,ILSTNL,ISYML,IOPT,MODEL,
     &                WORK(KL1AM),WORK(KEND1))
        CALL CC_T2SQ(WORK(KEND1),WORK(KL2AM),ISYM0)   

        CALL DZERO(WORK(KFOCKB),NORBT*NORBT)
      
        CALL CCSDT_FCK_R(WORK(KFOCKB),XIAJB,WORK(KTB1AM))

        CALL CCSDT_F3AM(TBAR3AM,WORK(KFOCKB),WORK(KINT1TB),
     &                  WORK(KINT2TB),WORK(KL2AM))

        IF (NONHF) THEN
          CALL CCSDT_E3AM(TBAR3AM,DUMMY,WORK(KL3AM),
     &                    WORK(KFLDB1),.FALSE.)
        END IF

        IF (LOCDBG) THEN
          WRITE(LUPRI,*) 'CCSDT_TBAR31_NODDY> norm^2(F3AM+E3AM):',
     &        DDOT(NT1AMX*NT1AMX*NT1AMX,TBAR3AM,1,TBAR3AM,1)
        END IF

      END IF

*=====================================================================*
*     Compute contribution from jacobian transformation:
*                 <Tbar^B_1|[H,tau_nu_3]|HF>/eps_3
*                + <Tbar^B_2|[H,tau_nu_3]|HF>/eps_3
*     Note: here we overwrite L1AM and L2AM with response multipliers
*=====================================================================*
      IF (.NOT. RHS_ONLY) THEN
        IOPT  = 3
        Call CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &                WORK(KL1AM),WORK(KEND1))
        if (locdbg) write(lupri,*)'ccsdt_tbar31_noddy> norm^2 of L2:',
     &     ddot(nt2amx,work(kend1),1,work(kend1),1) 
        CALL CC_T2SQ(WORK(KEND1),WORK(KL2AM),ISYMB)   

        IF (LOCDBG) THEN
          WRITE(LUPRI,*) 'CCSDT_TBAR31_NODDY> norm^2(F3AM+E3AM)-2:',
     &          DDOT(NT1AMX*NT1AMX*NT1AMX,TBAR3AM,1,TBAR3AM,1)
        END IF

        CALL CCSDT_L3AM_R(TBAR3AM,-FREQB,XINT1T0,XINT2T0,XIAJB,FOCK0,
     &                    WORK(KL1AM),WORK(KL2AM),SCR1,FOCKD,.FALSE.)

        IF (LOCDBG) THEN
          WRITE(LUPRI,*) 'CCSDT_TBAR31_NODDY> FREQB:',FREQB
          WRITE(LUPRI,*) 'CCSDT_TBAR31_NODDY> norm^2(XIAJB):',
     &        DDOT(NT1AMX*NT1AMX,XIAJB,1,XIAJB,1)
          WRITE(LUPRI,*) 'CCSDT_TBAR31_NODDY> norm^2(XINT1T0):',
     &        DDOT(NT1AMX*NVIRT*NVIRT,XINT1T0,1,XINT1T0,1)
          WRITE(LUPRI,*) 'CCSDT_TBAR31_NODDY> norm^2(XINT2T0):',
     &        DDOT(NT1AMX*NRHFT*NRHFT,XINT1T0,1,XINT1T0,1)
          WRITE(LUPRI,*) 'CCSDT_TBAR31_NODDY> norm^2(FOCK0):',
     &        DDOT(NORBT*NORBT,FOCK0,1,FOCK0,1)
          WRITE(LUPRI,*) 
     &     'CCSDT_TBAR31_NODDY> norm^2(F3AM+E3AM+L3AM(Tbar)):',
     &        DDOT(NT1AMX*NT1AMX*NT1AMX,TBAR3AM,1,TBAR3AM,1)
        END IF
      END IF

*=====================================================================*
*     Divide by orbital energy differences and frequency:
*     (for finite difference we solve the equations iteratively)
*=====================================================================*
      IF (RHS_ONLY) THEN
        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,TBAR3AM,1)
      ELSE
        KT3SCR = KL3AM ! scratch array for NONHF case

        CALL CCSDT_3AM(TBAR3AM,-FREQB,SCR1,FOCKD,
     &                 NONHF,WORK(KFIELD),.TRUE.,WORK(KT3SCR))

      END IF
*---------------------------------------------------------------------*
*     clean forbidden elements, print, if required and return:
*---------------------------------------------------------------------*
      CALL CCSDT_CLEAN_T3(TBAR3AM,NT1AMX,NVIRT,NRHFT)

      IF (PRINT_T3) THEN
        WRITE(LUPRI,*)'CCSDT_TBAR31_AM> first-order T3-bar vector:'
        WRITE(LUPRI,*)'CCSDT_TBAR31_AM> list,idlst:',LISTB,IDLSTB
        WRITE(LUPRI,*)'CCSDT_TBAR31_AM> freq,label:',FREQB,LABELB
        CALL PRINT_PT3_NODDY(TBAR3AM)
      END IF
      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'CCSDT_TBAR31_AM> LISTB,IDLSTB:',LISTB,IDLSTB
        WRITE(LUPRI,*) 'CCSDT_TBAR31_AM> FREQB,RHS_ONLY:',FREQB,RHS_ONLY
        WRITE(LUPRI,*) 'CCSDT_TBAR31_AM> LISTL,LISTR:',LISTL,LISTR
        WRITE(LUPRI,*) 'CCSDT_TBAR31_AM> norm^2(TBAR3AM):',
     &   DDOT(NT1AMX**3,TBAR3AM,1,TBAR3AM,1)
      END IF

      CALL QEXIT('CCTBAR31')

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_TBAR31_NODDY                   *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_FBC_NODDY(NEW_NODDY_F_CONT,
     &                           IBTRAN,IDOTS,DOTPROD,NBTRAN,MXVEC,
     &                           LISTL,LISTB,LISTC,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: compute a triples contribution to F/B matrix contraction
*
*             BCON = <L2|[[H,R_1],R_3]|HF>
*        
*     Written by Christof Haettig, May 2002 
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
#include <ccr1rsp.h>
#include <ccexci.h>
#include <dummy.h>
#include <inftap.h>

      LOGICAL LOCDBG, XI_ONLY
      PARAMETER (LOCDBG=.FALSE., XI_ONLY=.FALSE.)

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      LOGICAL NEW_NODDY_F_CONT
      CHARACTER*3 LISTL, LISTB, LISTC
      INTEGER LWORK, NBTRAN, MXVEC
      INTEGER IDOTS(MXVEC,NBTRAN), IBTRAN(3,NBTRAN)
#if defined (SYS_CRAY)
      REAL DOTPROD(MXVEC,NBTRAN), DDOT
      REAL WORK(LWORK)
      REAL FREQB, TCON
      REAL FREQLST
#else
      DOUBLE PRECISION DOTPROD(MXVEC,NBTRAN), DDOT
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION FREQB, TCON
      DOUBLE PRECISION FREQLST
#endif

      CHARACTER MODEL*(10)
      INTEGER KEND1, ITRAN, IDLSTB, IDLSTC, ISYMRB, KSCR1,
     &        KFOCKD, KFOCK0, KT1AMP0, KLAMP0, KLAMH0, KFOCKB, 
     &        KLAMPB, KLAMHB, KOME1, KLAM2, KINT1T0, KINT2T0, KDUM,
     &        KXIAJB, KYIAJB, KINT1SB, KINT2SB, LWRK1, IOPT, KTB3AM,
     &        KEND2, LWRK2, IVEC, IDLSTL, ISYML, ISYMRC, KL2AM,
     &        LUFOCK, KRC1AM, KINT1S0, KINT2S0, KFIELD, KFIELDAO

      INTEGER ILSTSYM

      IF (.NOT. NEW_NODDY_F_CONT) THEN
        RETURN
      ENDIF

      CALL QENTER('CCSDT_FBC_NODDY')

      IF (DIRECT) CALL QUIT('CCSDT_FBC_NODDY: DIRECT NOT IMPLEMENTED')

      IF (LOCDBG) WRITE(LUPRI,*) 'entered CCSDT_FBC_NODDY...'

*---------------------------------------------------------------------*
*     Memory allocation:
*---------------------------------------------------------------------*
      KSCR1   = 1
      KFOCKD  = KSCR1  + NT1AMX
      KEND1   = KFOCKD + NORBT

      KFOCK0  = KEND1
      KEND1   = KFOCK0  + NORBT*NORBT

      IF (NONHF) THEN
        KFIELD   = KEND1
        KFIELDAO = KFIELD   + NBAST*NBAST
        KEND1    = KFIELDAO + NBAST*NBAST
      END IF

      KT1AMP0 = KEND1
      KLAMP0  = KT1AMP0 + NT1AMX
      KLAMH0  = KLAMP0  + NLAMDT
      KEND1   = KLAMH0  + NLAMDT

      KFOCKB  = KEND1
      KLAMPB  = KFOCKB + NORBT*NORBT
      KLAMHB  = KLAMPB + NLAMDT
      KEND1   = KLAMHB + NLAMDT

      KOME1   = KEND1
      KRC1AM  = KOME1  + NT1AMX
      KEND1   = KRC1AM + NT1AMX

      KL2AM   = KEND1 
      KEND1   = KL2AM  + NT1AMX*NT1AMX

      KINT1T0 = KEND1
      KINT2T0 = KINT1T0 + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2T0 + NRHFT*NRHFT*NT1AMX

      KXIAJB  = KEND1
      KYIAJB  = KXIAJB  + NT1AMX*NT1AMX
      KEND1   = KYIAJB  + NT1AMX*NT1AMX

      KINT1S0 = KEND1
      KINT2S0 = KINT1S0 + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2S0 + NRHFT*NRHFT*NT1AMX

      KINT1SB = KEND1
      KINT2SB = KINT1SB + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2SB + NRHFT*NRHFT*NT1AMX

      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_FBC_NODDY')
      ENDIF

*---------------------------------------------------------------------*
*     Get zeroth-order Lambda matrices:
*---------------------------------------------------------------------*
      IOPT   = 1
      KDUM = KEND1
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
     &                      .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &                      NORBT,NLAMDT,NRHFT,NVIRT,NT1AMX)

*---------------------------------------------------------------------*
*     Loop right amplitude vectors:
*---------------------------------------------------------------------*
      DO ITRAN = 1, NBTRAN
        IDLSTB = IBTRAN(1,ITRAN)
        IDLSTC = IBTRAN(2,ITRAN)

        ISYMRB = ILSTSYM(LISTB,IDLSTB)
        ISYMRC = ILSTSYM(LISTC,IDLSTC)

        FREQB  = FREQLST(LISTB,IDLSTB)

        IF (LOCDBG) THEN
          WRITE(LUPRI,*) 'IDLSTB:',IDLSTB
          WRITE(LUPRI,*) 'FREQB:',FREQB
          WRITE(LUPRI,*) 'ISYMRB:',ISYMRB
          WRITE(LUPRI,*) 'IDLSTC:',IDLSTC
          WRITE(LUPRI,*) 'ISYMRC:',ISYMRC
        END IF

*---------------------------------------------------------------------*
*     Compute response amplitudes:
*---------------------------------------------------------------------*
      KTB3AM = KEND1
      KEND2  = KTB3AM + NT1AMX*NT1AMX*NT1AMX

      LWRK2  = LWORK  - KEND2
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_FBC_NODDY')
      ENDIF

      IF (LISTB(1:3).EQ.'R1 ' .OR. LISTB(1:3).EQ.'RE ' .OR.
     &    LISTB(1:3).EQ.'RC '                               ) THEN
        KDUM = KEND2
        CALL CCSDT_T31_NODDY(WORK(KTB3AM),LISTB,IDLSTB,FREQB,XI_ONLY,
     &                       .FALSE.,WORK(KINT1S0),WORK(KINT2S0),
     &                       .FALSE.,WORK(KINT1T0),WORK(KINT2T0),
     &                       .FALSE.,WORK(KXIAJB),WORK(KYIAJB),
     &                               WORK(KINT1SB),WORK(KINT2SB),
     &                       WORK(KLAMPB),WORK(KLAMHB),WORK(KFOCKB),
     &                       WORK(KLAMP0),WORK(KLAMH0),WORK(KFOCK0),
     &                       WORK(KDUM),WORK(KFOCKD),
     &                       WORK(KEND2),LWRK2)
      ELSE IF (LISTB(1:3).EQ.'R2 ' .OR. LISTB(1:3).EQ.'ER1') THEN
        CALL CCSDT_T32_NODDY(WORK(KTB3AM),LISTB,IDLSTB,FREQB,
     &                       WORK(KINT1S0),WORK(KINT2S0),
     &                       WORK(KLAMP0),WORK(KLAMH0),WORK(KFOCK0),
     &                       WORK(KFOCKD),DUMMY,DUMMY,
     &                       WORK(KSCR1),WORK(KEND2),LWRK2)
      ELSE
        CALL QUIT('Unknown or illegal list in CCSDT_FBC_NODDY.')
      END IF

*---------------------------------------------------------------------*
*     Compute contribution from <L_2|[[H,T^B_3],\tau_nu_1|HF>:
*---------------------------------------------------------------------*
      IVEC = 1
      DO WHILE(IDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
        IDLSTL = IDOTS(IVEC,ITRAN)

        IF (LWRK2 .LT. NT2AMX) THEN
           CALL QUIT('Insufficient space in CCSDT_FBC_NODDY')
        ENDIF

        ISYML = ILSTSYM(LISTL,IDLSTL)

        ! read left amplitudes from file and square up the doubles part
        IOPT   = 2
        Call CC_RDRSP(LISTL,IDLSTL,ISYML,IOPT,MODEL,DUMMY,WORK(KEND2))
        CALL CC_T2SQ(WORK(KEND2),WORK(KL2AM),ISYML)

        CALL DZERO(WORK(KOME1),NT1AMX)

        CALL T3_LEFT1(WORK(KOME1),WORK(KYIAJB),WORK(KXIAJB),
     &                WORK(KL2AM),WORK(KTB3AM))
 
        CALL T3_LEFT2(WORK(KOME1),WORK(KYIAJB),WORK(KXIAJB),
     &                WORK(KL2AM),WORK(KTB3AM))
 
        CALL T3_LEFT3(WORK(KOME1),WORK(KXIAJB),WORK(KL2AM),WORK(KTB3AM))

        ! read RC_1 amplitudes from file 
        ISYMRC = ILSTSYM(LISTC,IDLSTC)
        IOPT = 1
        Call CC_RDRSP(LISTC,IDLSTC,ISYMRC,IOPT,MODEL,WORK(KRC1AM),DUMMY)

        IF (ISYML.NE.MULD2H(ISYMRB,ISYMRC)) 
     &    CALL QUIT('Symmetry mismatch in CCSDT_FBC_NODDY.')

        TCON = DDOT(NT1AM(ISYMRC),WORK(KRC1AM),1,WORK(KOME1),1)
        DOTPROD(IVEC,ITRAN) = DOTPROD(IVEC,ITRAN) + TCON

        IF (LOCDBG) THEN
          WRITE(LUPRI,*) 'ITRAN,IVEC,IDLSTL:',ITRAN,IVEC,IDLSTL
          WRITE(LUPRI,*) 'TCON:',TCON
        END IF
      
        IVEC = IVEC + 1
      END DO ! WHILE

      END DO ! ITRAN

      CALL QEXIT('CCSDT_FBC_NODDY')

      RETURN
      END 
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_FBC_NODDY                      *
*---------------------------------------------------------------------*
