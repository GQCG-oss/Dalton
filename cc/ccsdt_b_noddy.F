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
      SUBROUTINE CCSDT_BMAT_NODDY(LISTA,IDLSTA,LISTB,IDLSTB,IOPTRES,
     &                            XLAMDP0,XLAMDH0,
     &                            OMEGA1,OMEGA2,OMEGA1EFF,OMEGA2EFF,
     &                            IDOTS,DOTPROD,LISTDP,ITRAN,
     &                            NBTRAN,MXVEC,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples contribution to B matrix transformation
*
*     (B T^A T^B)^eff_1,2 = (B T^A T^B)_1,2(CCSD) + (B T^A T^B)_1,2(T3)
*                            - A_1,2;3 (w_3 - w)^1 (B T^A T^B)_3
*
*        
*     Written by Christof Haettig, April 2002, based on CCSDT_TRIPLE.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <dummy.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccfield.h>
#include <ccorb.h>

      LOGICAL LOCDBG, PRINT_T3
      PARAMETER (LOCDBG=.FALSE., PRINT_T3=.FALSE.)
      INTEGER ISYM0
      PARAMETER (ISYM0=1)

      CHARACTER*3 LISTDP, LISTA, LISTB
      INTEGER LWORK, IOPTRES, ITRAN, MXVEC, NBTRAN, IDLSTA, IDLSTB
      INTEGER IDOTS(MXVEC,NBTRAN)

#if defined (SYS_CRAY)
      REAL FREQA, FREQB, FREQ, DDOT, FF
      REAL DOTPROD(MXVEC,NBTRAN), WORK(LWORK), TWO, ONE
      REAL XLAMDP0(*), XLAMDH0(*), SIGN
      REAL OMEGA1(*), OMEGA2(*)
      REAL OMEGA1EFF(*), OMEGA2EFF(*)
#else
      DOUBLE PRECISION FREQA, FREQB, FREQ, DDOT, FF
      DOUBLE PRECISION DOTPROD(MXVEC,NBTRAN), WORK(LWORK), TWO, ONE
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*), SIGN
      DOUBLE PRECISION OMEGA1(*), OMEGA2(*)
      DOUBLE PRECISION OMEGA1EFF(*), OMEGA2EFF(*)
#endif
      PARAMETER(ONE = 1.0D0,  TWO = 2.0D0 )

      CHARACTER*10 MODEL
      INTEGER KSCR1, KFOCKD, KEND1, KINT1T0, KINT2T0, KINT1TA, KINT2TA,
     &        KINT1TB, KINT2TB, KINT1SA, KINT2SA, KINT1SB, KINT2SB, 
     &        KINT1SAB, KINT2SAB, KTA3AM, KTB3AM, KFOCKA, KFOCKB,
     &        KLAMHB, KLAMPB, KLAMHA, KLAMPA, KT2AM0, KDUM, KXIAJB,
     &        KYIAJB, KFOCK0, KOMEGA2, LWRK1, KXINT, KEND2, IDX,
     &        LWRK2, KB3AM, IRECNR, KFIELDAO, KFLDA1, KFLDB1, KFLDBUF
      INTEGER IJ, NIJ, LUSIFC, INDEX, ILSTSYM, IOPT, LUFOCK, LUTEMP,
     &        ILLL, IDEL, ISYDIS, ISYMA, ISYMB, IVEC, ISYMC, IDLSTC,
     &        ISYMD, KFIELD, KINT1S0, KINT2S0, KFCKAR, KFCKBR
      INTEGER KTA1AM, KTB1AM, KEND1A, LWRK1A, KT2AMB, KT2AMA, KT3AM

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J 

      CALL QENTER('CCSDT_BMAT_NODDY')

      IF(DIRECT)CALL QUIT('DIRECT NOT IMPLEMENTED IN CCSDT_BMAT_NODDY')

*---------------------------------------------------------------------*
*     Memory allocation:
*---------------------------------------------------------------------*
      KSCR1   = 1
      KFOCKD  = KSCR1 + NT1AMX
      KFOCK0  = KFOCKD + NORBT
      KEND1   = KFOCK0 + NORBT*NORBT

      IF (NONHF) THEN
        KFIELD   = KEND1
        KFIELDAO = KFIELD   + NORBT*NORBT
        KFLDB1   = KFIELDAO + NORBT*NORBT
        KFLDA1   = KFLDB1   + NORBT*NORBT
        KFLDBUF  = KFLDA1   + NORBT*NORBT
        KEND1    = KFLDBUF  + NORBT*NORBT
      END IF

      KTA1AM  = KEND1
      KFOCKA  = KTA1AM + NT1AMX
      KFCKAR  = KFOCKA + NORBT*NORBT
      KLAMPA  = KFCKAR + NORBT*NORBT
      KLAMHA  = KLAMPA + NLAMDT
      KEND1   = KLAMHA + NLAMDT

      KTB1AM  = KEND1
      KFOCKB  = KTB1AM + NT1AMX
      KFCKBR  = KFOCKB + NORBT*NORBT
      KLAMPB  = KFCKBR + NORBT*NORBT
      KLAMHB  = KLAMPB + NLAMDT
      KEND1   = KLAMHB + NLAMDT

      KXIAJB   = KEND1
      KYIAJB   = KXIAJB + NT1AMX*NT1AMX
      KEND1    = KYIAJB + NT1AMX*NT1AMX

      KINT1T0 = KEND1
      KINT2T0 = KINT1T0 + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2T0 + NRHFT*NRHFT*NT1AMX

      KINT1S0 = KEND1
      KINT2S0 = KINT1S0 + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2S0 + NRHFT*NRHFT*NT1AMX

      KB3AM   = KEND1
      KEND1   = KB3AM + NT1AMX*NT1AMX*NT1AMX

      ! what is above has to be kept until the end...
      ! everything below might be overwritten 
      ! in CC_RHPART_NODDY / CCDOTRSP_NODDY
      KEND1A  = KEND1
      LWRK1A  = LWORK  - KEND1A
 
      KINT1TA = KEND1
      KINT2TA = KINT1TA + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2TA + NRHFT*NRHFT*NT1AMX

      KINT1TB = KEND1
      KINT2TB = KINT1TB + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2TB + NRHFT*NRHFT*NT1AMX

      KINT1SA = KEND1
      KINT2SA = KINT1SA + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2SA + NRHFT*NRHFT*NT1AMX

      KINT1SB = KEND1
      KINT2SB = KINT1SB + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2SB + NRHFT*NRHFT*NT1AMX

      KINT1SAB = KEND1
      KINT2SAB = KINT1SAB + NT1AMX*NVIRT*NVIRT
      KEND1    = KINT2SAB + NRHFT*NRHFT*NT1AMX

      KT3AM    = KEND1
      KEND1    = KT3AM + NT1AMX*NT1AMX*NT1AMX

      KOMEGA2  = KEND1
      KEND1    = KOMEGA2 + NT1AMX*NT1AMX

      KT2AMB   = KEND1
      KT2AMA   = KT2AMB + NT1AMX*NT1AMX
      KEND1    = KT2AMA + NT1AMX*NT1AMX

      KT2AM0 = KOMEGA2
      KTA3AM = KB3AM
      KTB3AM = KT3AM

      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_BMAT_NODDY')
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
*     read zeroth-order AO Fock matrix from file: 
*---------------------------------------------------------------------*
      LUFOCK = -1
      CALL GPOPEN(LUFOCK,'CC_FCKH','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUFOCK)
      READ(LUFOCK) (WORK(KFOCK0-1+I),I=1,N2BST(ISYM0))
      CALL GPCLOSE(LUFOCK,'KEEP')

      CALL CC_FCKMO(WORK(KFOCK0),XLAMDP0,XLAMDH0,
     &              WORK(KEND1),LWRK1,ISYM0,ISYM0,ISYM0)

*---------------------------------------------------------------------*
*     If needed get external field:
*---------------------------------------------------------------------*
      IF (NONHF .AND. NFIELD.GT.0) THEN
         CALL DZERO(WORK(KFIELDAO),NORBT*NORBT)
         DO I = 1, NFIELD
          FF = EFIELD(I)
          CALL CC_ONEP(WORK(KFIELDAO),WORK(KEND1),LWRK1,FF,1,LFIELD(I))
         ENDDO

         CALL DCOPY(NORBT*NORBT,WORK(KFIELDAO),1,WORK(KFIELD),1)
         CALL CC_FCKMO(WORK(KFIELD),XLAMDP0,XLAMDH0,
     *                 WORK(KEND1),LWRK1,1,1,1)
      ENDIF

*---------------------------------------------------------------------*
*     Compute some integrals:
*           XINT1S0 =  (CK|BD)
*           XINT2S0 =  (CK|LJ)
*           XINT1T0 =  (KC|BD)
*           XINT2T0 =  (KC|LJ)
*           XIAJB   = 2(IA|JB) - (IB|JA)
*           YIAJB   =  (IA|JB)
*---------------------------------------------------------------------*
      CALL CCSDT_INTS0_NODDY(.TRUE.,WORK(KXIAJB), WORK(KYIAJB),
     &                       .TRUE.,WORK(KINT1S0),WORK(KINT2S0),
     &                       .TRUE.,WORK(KINT1T0),WORK(KINT2T0),
     &                       XLAMDP0,XLAMDH0,
     &                       WORK(KEND1),LWRK1)

*---------------------------------------------------------------------*
*     Compute corrections to triples vector T^A_3, and corresponding 
*     lambda matrices and the XINT1SA,XINT2SA integrals and set FREQA:
*---------------------------------------------------------------------*
      IF (LISTA(1:3).EQ.'R1 ' .or. LISTA(1:3).EQ.'RE ' .or.
     &    LISTA(1:3).EQ.'RC '                              ) THEN

        KDUM = KEND1
        CALL CCSDT_T31_NODDY(WORK(KTA3AM),LISTA,IDLSTA,FREQA,.FALSE.,
     &                       .FALSE.,WORK(KINT1S0),WORK(KINT2S0),
     &                       .FALSE.,WORK(KINT1T0),WORK(KINT2T0),
     &                       .FALSE.,WORK(KXIAJB), WORK(KYIAJB),
     &                               WORK(KINT1SA),WORK(KINT2SA),
     &                       WORK(KLAMPA),WORK(KLAMHA),WORK(KFOCKA),
     &                       XLAMDP0,XLAMDH0,WORK(KFOCK0),
     &                       WORK(KDUM),WORK(KFOCKD),
     &                       WORK(KEND1),LWRK1)

      ELSE IF (LISTA(1:3).EQ.'R2 ' .or. LISTA(1:3).EQ.'ER1') THEN

        CALL CCSDT_T32_NODDY(WORK(KTA3AM),LISTA,IDLSTA,FREQA,
     &                       WORK(KINT1S0),WORK(KINT2S0),
     &                       XLAMDP0,XLAMDH0,WORK(KFOCK0),
     &                       WORK(KFOCKD),WORK(KFIELDAO),WORK(KFIELD),
     &                       WORK(KSCR1),WORK(KEND1),LWRK1)

        CALL CCLR_LAMTRA(XLAMDP0,WORK(KLAMPA),XLAMDH0,WORK(KLAMHA),
     &                   WORK(KTA1AM),ISYMA)

        CALL CCSDT_INTS1_NODDY(.TRUE.,WORK(KINT1SA),WORK(KINT2SA),
     &                         .FALSE.,DUMMY,DUMMY,
     &                         XLAMDP0,XLAMDH0,
     &                         WORK(KLAMPA),WORK(KLAMHA),
     &                         WORK(KEND1),LWRK1)

      ELSE
        CALL QUIT('Unknown LISTB in CCSDT_BMAT_NODDY.')
      END IF

      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KTA3AM),1)

      IF (NONHF) THEN
        LUTEMP = -1
        CALL GPOPEN(LUTEMP,'T3AMPA','UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        REWIND LUTEMP
        WRITE (LUTEMP) (WORK(KTA3AM-1+IDX),IDX=1,NT1AMX*NT1AMX*NT1AMX)
        CALL GPCLOSE(LUTEMP,'KEEP')
      END IF

*---------------------------------------------------------------------*
*     Compute corrections to triples vector T^B_3, and corresponding 
*     lambda matrices and the XINT1SB,XINT2SB integrals and set FREQB:
*---------------------------------------------------------------------*
      IF (LISTB(1:3).EQ.'R1 ' .or. LISTB(1:3).EQ.'RE ' .or.
     &    LISTB(1:3).EQ.'RC '                              ) THEN

        KDUM = KEND1
        CALL CCSDT_T31_NODDY(WORK(KTB3AM),LISTB,IDLSTB,FREQB,.FALSE.,
     &                       .FALSE.,WORK(KINT1S0),WORK(KINT2S0),
     &                       .FALSE.,WORK(KINT1T0),WORK(KINT1T0),
     &                       .FALSE.,WORK(KXIAJB), WORK(KYIAJB),
     &                               WORK(KINT1SB),WORK(KINT2SB),
     &                       WORK(KLAMPB),WORK(KLAMHB),WORK(KFOCKB),
     &                       XLAMDP0,XLAMDH0,WORK(KFOCK0),
     &                       WORK(KDUM),WORK(KFOCKD),
     &                       WORK(KEND1),LWRK1)

      ELSE IF (LISTB(1:3).EQ.'R2 ' .or. LISTB(1:3).EQ.'ER1') THEN

        CALL CCSDT_T32_NODDY(WORK(KTB3AM),LISTB,IDLSTB,FREQB,
     &                       WORK(KINT1S0),WORK(KINT2S0),
     &                       XLAMDP0,XLAMDH0,WORK(KFOCK0),
     &                       WORK(KFOCKD),WORK(KFIELDAO),WORK(KFIELD),
     &                       WORK(KSCR1),WORK(KEND1),LWRK1)

        CALL CCLR_LAMTRA(XLAMDP0,WORK(KLAMPB),XLAMDH0,WORK(KLAMHB),
     &                   WORK(KTB1AM),ISYMB)

        CALL CCSDT_INTS1_NODDY(.TRUE.,WORK(KINT1SB),WORK(KINT2SB),
     &                         .FALSE.,DUMMY,DUMMY,
     &                         XLAMDP0,XLAMDH0,
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KEND1),LWRK1)

      ELSE
        CALL QUIT('Unknown LISTB in CCSDT_BMAT_NODDY.')
      END IF

      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KTB3AM),1)

      IF (NONHF) THEN
        LUTEMP = -1
        CALL GPOPEN(LUTEMP,'T3AMPB','UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        REWIND LUTEMP
        WRITE (LUTEMP) (WORK(KTB3AM-1+IDX),IDX=1,NT1AMX*NT1AMX*NT1AMX)
        CALL GPCLOSE(LUTEMP,'KEEP')
      END IF

*---------------------------------------------------------------------*
*     Compute required integrals:
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KINT1SAB),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2SAB),NT1AMX*NRHFT*NRHFT)

      CALL DZERO(WORK(KINT1TA),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2TA),NT1AMX*NRHFT*NRHFT)

      CALL DZERO(WORK(KINT1TB),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2TB),NT1AMX*NRHFT*NRHFT)

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
               CALL QUIT('Insufficient space in CCSDT_BMAT_NODDY')
            ENDIF
 
C           ---------------------------
C           Read in batch of integrals.
C           ---------------------------
            CALL CCRDAO(WORK(KXINT),IDEL,1,WORK(KEND2),LWRK2,
     *                  IRECNR,DIRECT)
 
C           ----------------------------------
C           Calculate integrals needed in CC3:
C           ----------------------------------
            CALL CCSDT_TRAN1_R(WORK(KINT1TA),WORK(KINT2TA),
     &                         XLAMDP0,XLAMDH0,
     &                         WORK(KLAMPA),WORK(KLAMHA),
     &                         WORK(KXINT),IDEL)

            CALL CCSDT_TRAN1_R(WORK(KINT1TB),WORK(KINT2TB),
     &                         XLAMDP0,XLAMDH0,
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KXINT),IDEL)


            ! XINT1SAB = XINT1SAB + (C-barB K-barA|B D)
            ! XINT2SAB = XINT2SAB + (C-barB K-barA|L J)
            CALL CCSDT_TRAN3_R(WORK(KINT1SAB),WORK(KINT2SAB),
     &                         XLAMDP0,XLAMDH0,
     &                         WORK(KLAMPB),WORK(KLAMHA),
     &                         XLAMDP0,XLAMDH0,
     &                         WORK(KXINT),IDEL)

            ! XINT1SAB = XINT1SAB + (C-barA K-barB|B D)
            ! XINT2SAB = XINT2SAB + (C-barA K-barB|L J)
            CALL CCSDT_TRAN3_R(WORK(KINT1SAB),WORK(KINT2SAB),
     &                         XLAMDP0,XLAMDH0,
     &                         WORK(KLAMPA),WORK(KLAMHB),
     &                         XLAMDP0,XLAMDH0,
     &                         WORK(KXINT),IDEL)

            ! XINT1SAB = XINT1SAB + (C-barB K|B-barA D)
            ! XINT2SAB = XINT2SAB + (C-barB K|L J-barA)
            CALL CCSDT_TRAN3_R(WORK(KINT1SAB),WORK(KINT2SAB),
     &                         XLAMDP0,XLAMDH0,
     &                         WORK(KLAMPB),XLAMDH0,
     &                         WORK(KLAMPA),WORK(KLAMHA),
     &                         WORK(KXINT),IDEL)

            ! XINT1SAB = XINT1SAB + (C-barA K|B-barB D)
            ! XINT2SAB = XINT2SAB + (C-barA K|L J-barB)
            CALL CCSDT_TRAN3_R(WORK(KINT1SAB),WORK(KINT2SAB),
     &                         XLAMDP0,XLAMDH0,
     &                         WORK(KLAMPA),XLAMDH0,
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KXINT),IDEL)

            ! XINT1SAB = XINT1SAB + (C K-barB|B-barA D)
            ! XINT2SAB = XINT2SAB + (C K-barB|L J-barA)
            CALL CCSDT_TRAN3_R(WORK(KINT1SAB),WORK(KINT2SAB),
     &                         XLAMDP0,XLAMDH0,
     &                         XLAMDP0,WORK(KLAMHB),
     &                         WORK(KLAMPA),WORK(KLAMHA),
     &                         WORK(KXINT),IDEL)

            ! XINT1SAB = XINT1SAB + (C K-barA|B-barB D)
            ! XINT2SAB = XINT2SAB + (C K-barA|L J-barB)
            CALL CCSDT_TRAN3_R(WORK(KINT1SAB),WORK(KINT2SAB),
     &                         XLAMDP0,XLAMDH0,
     &                         XLAMDP0,WORK(KLAMHA),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KXINT),IDEL)

         END DO   
      END DO  

C     CALL CCSDT_INTS2_NODDY(WORK(KXINT1SAB),WORK(KXINT2SAB),
C    &                       XLAMDP0,XLAMDH0,
C    &                       WORK(KLAMPB),WORK(KLAMHB),
C    &                       WORK(KLAMPA),WORK(KLAMHA),
C    &                       WORK(KEND1),LWRK1)

C     CALL CCSDT_INTS1_NODDY(.FALSE.,WORK(KINT1SA),WORK(KINT2SA),
C    &                       .TRUE.,WORK(KINT1TA),WORK(KINT2TA),
C    &                       XLAMDP0,XLAMDH0,
C    &                       WORK(KLAMPA),WORK(KLAMHA),
C    &                       WORK(KEND1),LWRK1)

C     CALL CCSDT_INTS1_NODDY(.FALSE.,WORK(KINT1SB),WORK(KINT2SB),
C    &                       .TRUE.,WORK(KINT1TB),WORK(KINT2TB),
C    &                       XLAMDP0,XLAMDH0,
C    &                       WORK(KLAMPB),WORK(KLAMHB),
C    &                       WORK(KEND1),LWRK1)

*---------------------------------------------------------------------*
*     Compute corrections to doubles result vector from T^A_3, T^B_3:
*          omega2 +=  P^AB <mu_2|[[H,T^A_1],T^B_3]|HF>
*---------------------------------------------------------------------*
      ! read RA_1 amplitudes from file 
      ISYMA = ILSTSYM(LISTA,IDLSTA)
      IOPT = 1
      Call CC_RDRSP(LISTA,IDLSTA,ISYMA,IOPT,MODEL,WORK(KTA1AM),DUMMY)

      ! read RB_1 amplitudes from file 
      ISYMB = ILSTSYM(LISTB,IDLSTB)
      IOPT = 1
      Call CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,WORK(KTB1AM),DUMMY)

      ! compute fock matrix like contribution to [H,T^A_1]
      CALL DZERO(WORK(KFCKAR),NORBT*NORBT)
      CALL CCSDT_FCK_R(WORK(KFCKAR),WORK(KXIAJB),WORK(KTA1AM))
 
      ! compute fock matrix like contribution to [H,T^B_1]
      CALL DZERO(WORK(KFCKBR),NORBT*NORBT)
      CALL CCSDT_FCK_R(WORK(KFCKBR),WORK(KXIAJB),WORK(KTB1AM))


      ! <mu_2|[[H,T^A_1],T^B_3]|HF>
      CALL DZERO(WORK(KOMEGA2),NT1AMX*NT1AMX)
      CALL CCSDT_OMEGA2(WORK(KOMEGA2),WORK(KINT1TA),WORK(KINT2TA),
     &                  WORK(KTB3AM),WORK(KFCKAR))

      DO I = 1,NT1AMX
         DO J = 1,I
            IJ = NT1AMX*(I-1) + J
            NIJ = INDEX(I,J)
            OMEGA2(NIJ) = OMEGA2(NIJ) + WORK(KOMEGA2+IJ-1)
         END DO
      END DO


      ! <mu_2|[[H,T^B_1],T^A_3]|HF>
      CALL DZERO(WORK(KOMEGA2),NT1AMX*NT1AMX)

      CALL CCSDT_OMEGA2(WORK(KOMEGA2),WORK(KINT1TB),WORK(KINT2TB),
     &                  WORK(KTA3AM),WORK(KFCKBR))

      DO I = 1,NT1AMX
         DO J = 1,I
            IJ = NT1AMX*(I-1) + J
            NIJ = INDEX(I,J)
            OMEGA2(NIJ) = OMEGA2(NIJ) + WORK(KOMEGA2+IJ-1)
         END DO
      END DO

*---------------------------------------------------------------------*
*     Compute triples vector (B T^A T^B)_3:
*---------------------------------------------------------------------*

      CALL DZERO(WORK(KB3AM),NT1AMX*NT1AMX*NT1AMX)

      if (.true. ) then 

      IF (LWRK1 .LT. NT2AMX) THEN
         CALL QUIT('Insufficient space in CCSDT_BMAT_NODDY')
      ENDIF

      ! add <nu_3|[H^AB,T^0_2]|HF>
      IOPT = 2
      KDUM = KEND1
      CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KDUM),WORK(KEND1))
      CALL CC_T2SQ(WORK(KEND1),WORK(KT2AM0),ISYM0)

      CALL CCSDT_T3AM_R(WORK(KB3AM),0.0D0,
     &                  WORK(KINT1SAB),WORK(KINT2SAB),WORK(KT2AM0),
     &                  WORK(KSCR1),WORK(KFOCKD),
     &                  .FALSE.,DUMMY,.FALSE.)


      ! add <nu_3|[H^A,T^B_2]|HF>
      ISYMB = ILSTSYM(LISTB,IDLSTB)
      IOPT  = 2
      KDUM = KEND1
      CALL CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &              WORK(KDUM),WORK(KEND1))
      CALL CCLR_DIASCL(WORK(KEND1),TWO,ISYMB) 
      CALL CC_T2SQ(WORK(KEND1),WORK(KT2AMB),ISYMB)

      CALL CCSDT_T3AM_R(WORK(KB3AM),0.0D0,
     &                  WORK(KINT1SA),WORK(KINT2SA),WORK(KT2AMB),
     &                  WORK(KSCR1),WORK(KFOCKD),
     &                  .FALSE.,DUMMY,.FALSE.)

      ! add <nu_3|[H^B,T^A_2]|HF>
      ISYMA = ILSTSYM(LISTA,IDLSTA)
      IOPT  = 2
      KDUM = KEND1
      CALL CC_RDRSP(LISTA,IDLSTA,ISYMA,IOPT,MODEL,
     &              WORK(KDUM),WORK(KEND1))
      CALL CCLR_DIASCL(WORK(KEND1),TWO,ISYMA) 
      CALL CC_T2SQ(WORK(KEND1),WORK(KT2AMA),ISYMA)

      CALL CCSDT_T3AM_R(WORK(KB3AM),0.0D0,
     &                  WORK(KINT1SB),WORK(KINT2SB),WORK(KT2AMA),
     &                  WORK(KSCR1),WORK(KFOCKD),
     &                  .FALSE.,DUMMY,.FALSE.)

      IF (.TRUE. .AND. NONHF .AND. NFIELD.GT.0) THEN
cch
        write(lupri,*) 'norm^2(field) before FF:',
     &    ddot(norbt*norbt,work(kfield),1,work(kfield),1)
        write(lupri,*) 'norm^2(b3am) before FF:',
     &    ddot(nt1amx**3,work(kb3am),1,work(kb3am),1)
cch
        CALL CCSDT_XKSI3_1(WORK(KB3AM),WORK(KFIELD),
     &                     WORK(KT2AMB),WORK(KT2AMA),ONE)
        CALL CCSDT_XKSI3_1(WORK(KB3AM),WORK(KFIELD),
     &                     WORK(KT2AMA),WORK(KT2AMB),ONE)
cch
        write(lupri,*) 'norm^2(b3am) after FF-1:',
     &    ddot(nt1amx**3,work(kb3am),1,work(kb3am),1)
cch

        ! calculate one-index transf. field integrals: V^A = [V,T1^A]
        CALL DCOPY(NORBT*NORBT,WORK(KFIELDAO),1,WORK(KFLDBUF),1)
        CALL DCOPY(NORBT*NORBT,WORK(KFIELDAO),1,WORK(KFLDA1),1)
        CALL CC_FCKMO(WORK(KFLDBUF),WORK(KLAMPA),XLAMDH0,
     &                WORK(KEND1),LWRK1,ISYMA,ISYM0,ISYMA)
        CALL CC_FCKMO(WORK(KFLDA1),XLAMDP0,WORK(KLAMHA),
     &                WORK(KEND1),LWRK1,ISYM0,ISYMA,ISYMA)
        CALL DAXPY(NORBT*NORBT,ONE,WORK(KFLDBUF),1,WORK(KFLDA1),1)
ctest
        write(lupri,*) 'norm^2(flda1) before FF:',
     &    ddot(norbt*norbt,work(kflda1),1,work(kflda1),1)

C       CALL DSCAL(NORBT*NORBT,-1.0D0,WORK(KFLDA1),1)
ctest

        LUTEMP = -1
        CALL GPOPEN(LUTEMP,'T3AMPB','UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        REWIND LUTEMP
        READ (LUTEMP) (WORK(KT3AM-1+IDX),IDX=1,NT1AMX*NT1AMX*NT1AMX)
        CALL GPCLOSE(LUTEMP,'DELETE')

        CALL CCSDT_XKSI3_2(WORK(KB3AM),WORK(KFLDA1),WORK(KT3AM))

        ! calculate one-index transf. field integrals: V^B = [V,T1^B]
        CALL DCOPY(NORBT*NORBT,WORK(KFIELDAO),1,WORK(KFLDBUF),1)
        CALL DCOPY(NORBT*NORBT,WORK(KFIELDAO),1,WORK(KFLDB1),1)
        CALL CC_FCKMO(WORK(KFLDBUF),WORK(KLAMPB),XLAMDH0,
     &                WORK(KEND1),LWRK1,ISYMB,ISYM0,ISYMB)
        CALL CC_FCKMO(WORK(KFLDB1),XLAMDP0,WORK(KLAMHB),
     &                WORK(KEND1),LWRK1,ISYM0,ISYMB,ISYMB)
        CALL DAXPY(NORBT*NORBT,ONE,WORK(KFLDBUF),1,WORK(KFLDB1),1)
ctest
        write(lupri,*) 'norm^2(fldb1) before FF:',
     &    ddot(norbt*norbt,work(kfldb1),1,work(kfldb1),1)

C       CALL DSCAL(NORBT*NORBT,-1.0D0,WORK(KFLDB1),1)
ctest

        LUTEMP = -1
        CALL GPOPEN(LUTEMP,'T3AMPA','UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        REWIND LUTEMP
        READ (LUTEMP) (WORK(KT3AM-1+IDX),IDX=1,NT1AMX*NT1AMX*NT1AMX)
        CALL GPCLOSE(LUTEMP,'DELETE')

        CALL CCSDT_XKSI3_2(WORK(KB3AM),WORK(KFLDB1),WORK(KT3AM))
cch
        write(lupri,*) 'norm^2(b3am) after FF:',
     &    ddot(nt1amx**3,work(kb3am),1,work(kb3am),1)
cch
      END IF

      else 

      CALL CCSDT_B3AM(WORK(KB3AM),
     &                WORK(KINT1SAB),WORK(KINT2SAB),WORK(KFOCKD),
     &                LISTA,IDLSTA,WORK(KINT1SA),WORK(KINT2SA),
     &                LISTB,IDLSTB,WORK(KINT1SB),WORK(KINT2SB),
     &                WORK(KSCR1),WORK(KT2AM0),WORK(KEND1),LWRK1)

      end if

C     if (ioptres.eq.5) then
C       write(lupri,*) 'IOPTRES=5... set B3AM to zero!!!'
C       call dzero(work(kb3am),nt1amx*nt1amx*nt1amx)
C     end if

      if (print_t3) then
        write(lupri,*)'CCSDT_B_NODDY> vector types:',lista,listb
        write(lupri,*)'CCSDT_B_NODDY> idlsts:',idlsta,idlstb
        write(lupri,*)'CCSDT_B_NODDY> freq:',freqa+freqb
        call print_pt3_noddy(WORK(KB3AM))
      end if

*---------------------------------------------------------------------*
*     Now we split:
*       for IOPTRES < 5 we compute an effective rhs vector
*       for IOPTRES = 5 we compute contractions L B T^A T^B
*---------------------------------------------------------------------*
      IF (IOPTRES.GE.1 .AND. IOPTRES.LE.4) THEN

         CALL DCOPY(NT1AMX,OMEGA1,1,OMEGA1EFF,1)
         CALL DCOPY(NT2AMX,OMEGA2,1,OMEGA2EFF,1)

         FREQ = FREQA + FREQB

         CALL CC_RHPART_NODDY(OMEGA1EFF,OMEGA2EFF,WORK(KB3AM),FREQ,
     &                        WORK(KFOCKD),WORK(KFOCK0),WORK(KFIELD),
     &                        WORK(KXIAJB),WORK(KINT1T0),WORK(KINT2T0),
     &                        WORK(KEND1A),LWRK1A)
      
      ELSE IF (IOPTRES.EQ.5) THEN

        SIGN = -1.0D0
        CALL CCDOTRSP_NODDY(DUMMY,DUMMY,WORK(KB3AM),SIGN,
     &                      ITRAN,LISTDP,IDOTS,DOTPROD,MXVEC,
     &                      XLAMDP0,XLAMDH0,WORK(KFOCK0),WORK(KFOCKD),
     &                      WORK(KXIAJB), WORK(KYIAJB),
     &                      WORK(KINT1T0),WORK(KINT2T0),
     &                      WORK(KINT1S0),WORK(KINT2S0),
     &                      'CCSDT_B_NODDY',LOCDBG,.FALSE.,.FALSE.,
     &                      WORK(KEND1A),LWRK1A)

      ELSE
        CALL QUIT('Illegal value for IOPTRES IN CCSDT_BMAT_NODDY')
      END IF

      CALL QEXIT('CCSDT_BMAT_NODDY')
      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_BMAT_NODDY                     *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_B3AM(B3AM,XINT1SAB,XINT2SAB,FOCKD,
     &                      LISTA,IDLSTA,XINT1SA,XINT2SA,
     &                      LISTB,IDLSTB,XINT1SB,XINT2SB,
     &                      SCR1,T2AM,WORK,LWORK)
*---------------------------------------------------------------------*
*    Purpose: compute triples result of B matrix transformation
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <dummy.h>
#include <ccsdsym.h>

      INTEGER ISYM0
      PARAMETER (ISYM0=1)

      CHARACTER*3 LISTA, LISTB
      INTEGER LWORK, IDLSTB, IDLSTA

#if defined (SYS_CRAY)
      REAL B3AM(*), SCR1(*), T2AM(*), WORK(*)
      REAL XINT1SA(*), XINT2SA(*), XINT1SB(*), XINT2SB(*)
      REAL XINT1SAB(*), XINT2SAB(*), FOCKD(*)
      REAL TWO, ZERO
#else
      DOUBLE PRECISION B3AM(*), SCR1(*), T2AM(*), WORK(*)
      DOUBLE PRECISION XINT1SA(*), XINT2SA(*), XINT1SB(*), XINT2SB(*)
      DOUBLE PRECISION XINT1SAB(*), XINT2SAB(*), FOCKD(*)
      DOUBLE PRECISION TWO, ZERO
#endif
      PARAMETER( TWO = 2.0D0, ZERO = 0.0D0 )

      CHARACTER*10 MODEL
      INTEGER IOPT, ISYMA, ISYMB, ILSTSYM

      IF (LWORK .LT. NT2AMX) THEN
         CALL QUIT('Insufficient space in CCSDT_B3AM')
      ENDIF

      ! add <nu_3|[H^AB,T^0_2]|HF>
      IOPT = 2
      CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,DUMMY,WORK)
      CALL CC_T2SQ(WORK,T2AM,ISYM0)

      CALL CCSDT_T3AM_R(B3AM,ZERO,XINT1SAB,XINT2SAB,T2AM,
     &                  SCR1,FOCKD,.FALSE.,DUMMY,.FALSE.)


      ! add <nu_3|[H^A,T^B_2]|HF>
      ISYMB = ILSTSYM(LISTB,IDLSTB)
      IOPT  = 2
      CALL CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,DUMMY,WORK)
      CALL CCLR_DIASCL(WORK,TWO,ISYMB) 
      CALL CC_T2SQ(WORK,T2AM,ISYMB)

      CALL CCSDT_T3AM_R(B3AM,ZERO,XINT1SA,XINT2SA,T2AM,
     &                  SCR1,FOCKD,.FALSE.,DUMMY,.FALSE.)

      ! add <nu_3|[H^B,T^A_2]|HF>
      ISYMA = ILSTSYM(LISTA,IDLSTA)
      IOPT  = 2
      CALL CC_RDRSP(LISTA,IDLSTA,ISYMA,IOPT,MODEL,DUMMY,WORK)
      CALL CCLR_DIASCL(WORK,TWO,ISYMA) 
      CALL CC_T2SQ(WORK,T2AM,ISYMA)

      CALL CCSDT_T3AM_R(B3AM,ZERO,XINT1SB,XINT2SB,T2AM,
     &                  SCR1,FOCKD,.FALSE.,DUMMY,.FALSE.)

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_B3AM
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_INTS2_NODDY(XINT1SAB,XINT2SAB,XLAMDP0,XLAMDH0,
     &                             XLAMDPB,XLAMDHB,XLAMDPA,XLAMDHA,
     &                             WORK,LWORK)
*---------------------------------------------------------------------*
*  Purpose:
*     Loop over distributions of integrals and compute second-order
*     response versions XINT1SAB=(ck|db)-AB and XINT2SAB=(ck|lj)-AB
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccorb.h>

      INTEGER LWORK
#if defined (SYS_CRAY)
      REAL WORK(*), XINT1SAB(*), XINT2SAB(*)
      REAL XLAMDP0(*), XLAMDH0(*)
      REAL XLAMDPA(*), XLAMDHA(*)
      REAL XLAMDPB(*), XLAMDHB(*)
#else
      DOUBLE PRECISION WORK(*), XINT1SAB(*), XINT2SAB(*)
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*)
      DOUBLE PRECISION XLAMDPA(*), XLAMDHA(*)
      DOUBLE PRECISION XLAMDPB(*), XLAMDHB(*)
#endif

      INTEGER ISYMD, ILLL, IDEL, ISYDIS, KXINT, KEND2, LWRK2, IRECNR

      CALL DZERO(XINT1SAB,NT1AMX*NVIRT*NVIRT)
      CALL DZERO(XINT2SAB,NT1AMX*NRHFT*NRHFT)


      DO ISYMD = 1, NSYM
         DO ILLL = 1,NBAS(ISYMD)
            IDEL   = IBAS(ISYMD) + ILLL
            ISYDIS = ISYMD
 
C           ----------------------------
C           Work space allocation no. 2.
C           ----------------------------
            KXINT  = 1
            KEND2  = KXINT + NDISAO(ISYDIS)
            LWRK2  = LWORK - KEND2
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
               CALL QUIT('Insufficient space in CCSDT_INTS2_NODDY')
            ENDIF
 
C           ---------------------------
C           Read in batch of integrals.
C           ---------------------------
            CALL CCRDAO(WORK(KXINT),IDEL,1,WORK(KEND2),LWRK2,
     *                  IRECNR,DIRECT)
 
C           ----------------------------------
C           Calculate integrals needed in CC3:
C           ----------------------------------
            ! XINT1SAB = XINT1SAB + (C-barB K-barA|B D)
            ! XINT2SAB = XINT2SAB + (C-barB K-barA|L J)
            CALL CCSDT_TRAN3_R(XINT1SAB,XINT2SAB,XLAMDP0,XLAMDH0,
     &                         XLAMDPB,XLAMDHA,  XLAMDP0,XLAMDH0,
     &                         WORK(KXINT),IDEL)

            ! XINT1SAB = XINT1SAB + (C-barA K-barB|B D)
            ! XINT2SAB = XINT2SAB + (C-barA K-barB|L J)
            CALL CCSDT_TRAN3_R(XINT1SAB,XINT2SAB,XLAMDP0,XLAMDH0,
     &                         XLAMDPA,XLAMDHB,  XLAMDP0,XLAMDH0,
     &                         WORK(KXINT),IDEL)

            ! XINT1SAB = XINT1SAB + (C-barB K|B-barA D)
            ! XINT2SAB = XINT2SAB + (C-barB K|L J-barA)
            CALL CCSDT_TRAN3_R(XINT1SAB,XINT2SAB,XLAMDP0,XLAMDH0,
     &                         XLAMDPB,XLAMDH0,  XLAMDPA,XLAMDHA,
     &                         WORK(KXINT),IDEL)

            ! XINT1SAB = XINT1SAB + (C-barA K|B-barB D)
            ! XINT2SAB = XINT2SAB + (C-barA K|L J-barB)
            CALL CCSDT_TRAN3_R(XINT1SAB,XINT2SAB,XLAMDP0,XLAMDH0,
     &                         XLAMDPA,XLAMDH0,  XLAMDPB,XLAMDHB,
     &                         WORK(KXINT),IDEL)

            ! XINT1SAB = XINT1SAB + (C K-barB|B-barA D)
            ! XINT2SAB = XINT2SAB + (C K-barB|L J-barA)
            CALL CCSDT_TRAN3_R(XINT1SAB,XINT2SAB,XLAMDP0,XLAMDH0,
     &                         XLAMDP0,XLAMDHB,  XLAMDPA,XLAMDHA,
     &                         WORK(KXINT),IDEL)

            ! XINT1SAB = XINT1SAB + (C K-barA|B-barB D)
            ! XINT2SAB = XINT2SAB + (C K-barA|L J-barB)
            CALL CCSDT_TRAN3_R(XINT1SAB,XINT2SAB,XLAMDP0,XLAMDH0,
     &                         XLAMDP0,XLAMDHA,  XLAMDPB,XLAMDHB,
     &                         WORK(KXINT),IDEL)

         END DO   
      END DO  

      RETURN
      END
*---------------------------------------------------------------------*
*             END OF SUBROUTINE CCSDT_INTS2_NODDY                     *
*---------------------------------------------------------------------*
