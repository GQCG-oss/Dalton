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

      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.FALSE.)
      INTEGER ISYM0
      PARAMETER (ISYM0=1)

      CHARACTER*3 LISTDP, LISTA, LISTB
      INTEGER LWORK, IOPTRES, ITRAN, MXVEC, NBTRAN, IDLSTA, IDLSTB
      INTEGER IDOTS(MXVEC,NBTRAN)

#if defined (SYS_CRAY)
      REAL FREQA, FREQB, FREQ, DDOT, FF
      REAL DOTPROD(MXVEC,NBTRAN), WORK(LWORK), TWO
      REAL XLAMDP0(*), XLAMDH0(*), SIGN
      REAL OMEGA1(*), OMEGA2(*)
      REAL OMEGA1EFF(*), OMEGA2EFF(*)
#else
      DOUBLE PRECISION FREQA, FREQB, FREQ, DDOT, FF
      DOUBLE PRECISION DOTPROD(MXVEC,NBTRAN), WORK(LWORK), TWO
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*), SIGN
      DOUBLE PRECISION OMEGA1(*), OMEGA2(*)
      DOUBLE PRECISION OMEGA1EFF(*), OMEGA2EFF(*)
#endif
      PARAMETER( TWO = 2.0D0 )

      CHARACTER*10 MODEL
      INTEGER KSCR1, KFOCKD, KEND1, KINT1T0, KINT2T0, KINT1TA, KINT2TA,
     &        KINT1TB, KINT2TB, KINT1SA, KINT2SA, KINT1SB, KINT2SB, 
     &        KINT1SAB, KINT2SAB, KTA3AM, KTB3AM, KFOCKA, KFOCKB,
     &        KLAMHB, KLAMPB, KLAMHA, KLAMPA, KT2AM, KDUM, KXIAJB,
     &        KYIAJB, KFOCK0, KOMEGA2, LWRK1, KXINT, KEND2,
     &        LWRK2, KB3AM, IRECNR
      INTEGER IJ, NIJ, LUSIFC, INDEX, ILSTSYM, IOPT, LUFOCK,
     &        ILLL, IDEL, ISYDIS, ISYMA, ISYMB, IVEC, ISYMC, IDLSTC,
     &        ISYMD, KFIELD, KINT1S0, KINT2S0, KFCKAR, KFCKBR
      INTEGER KTA1AM, KTB1AM, KEND1A, LWRK1A

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
        KFIELD = KEND1
        KEND1  = KFIELD + NORBT*NORBT
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

      KTB3AM   = KEND1
      KEND1    = KTB3AM + NT1AMX*NT1AMX*NT1AMX

      KOMEGA2  = KEND1
      KEND1    = KOMEGA2 + NT1AMX*NT1AMX

      KT2AM  = KOMEGA2
      KTA3AM = KB3AM

      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_BMAT_NODDY')
      ENDIF

*---------------------------------------------------------------------*
*     Read SCF orbital energies from file:
*---------------------------------------------------------------------*
      LUSIFC = 0
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
      LUFOCK = 0
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
      IF ((NONHF) .AND. NFIELD .GT. 0) THEN
         CALL DZERO(WORK(KFIELD),NORBT*NORBT)
         DO I = 1, NFIELD
            FF = EFIELD(I)
            CALL CC_ONEP(WORK(KFIELD),WORK(KEND1),LWRK1,FF,1,LFIELD(I))
         ENDDO
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
      IF (LISTA(1:3).EQ.'R1 ' .or. LISTA(1:3).EQ.'RE ') THEN

        CALL CCSDT_T31_NODDY(WORK(KTA3AM),LISTA,IDLSTA,FREQA,.FALSE.,
     &                       .FALSE.,WORK(KINT1S0),WORK(KINT2S0),
     &                       .FALSE.,WORK(KINT1T0),WORK(KINT2T0),
     &                       .FALSE.,WORK(KXIAJB), WORK(KYIAJB),
     &                               WORK(KINT1SA),WORK(KINT2SA),
     &                       WORK(KLAMPA),WORK(KLAMHA),WORK(KFOCKA),
     &                       XLAMDP0,XLAMDH0,WORK(KFOCK0),
     &                       WORK(KDUM),WORK(KFOCKD),
     &                       WORK(KEND1),LWRK1)

        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KTA3AM),1)

      ELSE
        CALL QUIT('Unkown LISTB in CCSDT_BMAT_NODDY.')
      END IF

*---------------------------------------------------------------------*
*     Compute corrections to triples vector T^B_3, and corresponding 
*     lambda matrices and the XINT1SB,XINT2SB integrals and set FREQB:
*---------------------------------------------------------------------*
      IF (LISTB(1:3).EQ.'R1 ' .or. LISTB(1:3).EQ.'RE ') THEN

        CALL CCSDT_T31_NODDY(WORK(KTB3AM),LISTB,IDLSTB,FREQB,.FALSE.,
     &                       .FALSE.,WORK(KINT1S0),WORK(KINT2S0),
     &                       .FALSE.,WORK(KINT1T0),WORK(KINT1T0),
     &                       .FALSE.,WORK(KXIAJB), WORK(KYIAJB),
     &                               WORK(KINT1SB),WORK(KINT2SB),
     &                       WORK(KLAMPB),WORK(KLAMHB),WORK(KFOCKB),
     &                       XLAMDP0,XLAMDH0,WORK(KFOCK0),
     &                       WORK(KDUM),WORK(KFOCKD),
     &                       WORK(KEND1),LWRK1)

        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KTB3AM),1)

      ELSE
        CALL QUIT('Unkown LISTB in CCSDT_BMAT_NODDY.')
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

C     CALL CCSDT_INTS1_NODDY(.TRUE.,WORK(KINT1SA),WORK(KINT2SA),
C    &                       .TRUE.,WORK(KINT1TA),WORK(KINT2TA),
C    &                       XLAMDP0,XLAMDH0,
C    &                       WORK(KLAMPA),WORK(KLAMHA),
C    &                       WORK(KEND1),LWRK1)

C     CALL CCSDT_INTS1_NODDY(.TRUE.,WORK(KINT1SB),WORK(KINT2SB),
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
      IF (LWRK1 .LT. NT2AMX) THEN
         CALL QUIT('Insufficient space in CCSDT_BMAT_NODDY')
      ENDIF

      CALL DZERO(WORK(KB3AM),NT1AMX*NT1AMX*NT1AMX)

      ! add <nu_3|[H^AB,T^0_2]|HF>
      IOPT = 2
      CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KDUM),WORK(KEND1))
      CALL CC_T2SQ(WORK(KEND1),WORK(KT2AM),ISYM0)

      CALL CCSDT_T3AM_R(WORK(KB3AM),0.0D0,
     &                  WORK(KINT1SAB),WORK(KINT2SAB),WORK(KT2AM),
     &                  WORK(KSCR1),WORK(KFOCKD),
     &                  .FALSE.,DUMMY,.FALSE.)


      ! add <nu_3|[H^A,T^B_2]|HF>
      ISYMB = ILSTSYM(LISTB,IDLSTB)
      IOPT  = 2
      CALL CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &              WORK(KDUM),WORK(KEND1))
      CALL CCLR_DIASCL(WORK(KEND1),TWO,ISYMB) 
      CALL CC_T2SQ(WORK(KEND1),WORK(KT2AM),ISYMB)

      CALL CCSDT_T3AM_R(WORK(KB3AM),0.0D0,
     &                  WORK(KINT1SA),WORK(KINT2SA),WORK(KT2AM),
     &                  WORK(KSCR1),WORK(KFOCKD),
     &                  .FALSE.,DUMMY,.FALSE.)

      ! add <nu_3|[H^B,T^A_2]|HF>
      ISYMA = ILSTSYM(LISTA,IDLSTA)
      IOPT  = 2
      CALL CC_RDRSP(LISTA,IDLSTA,ISYMA,IOPT,MODEL,
     &              WORK(KDUM),WORK(KEND1))
      CALL CCLR_DIASCL(WORK(KEND1),TWO,ISYMA) 
      CALL CC_T2SQ(WORK(KEND1),WORK(KT2AM),ISYMA)

      CALL CCSDT_T3AM_R(WORK(KB3AM),0.0D0,
     &                  WORK(KINT1SB),WORK(KINT2SB),WORK(KT2AM),
     &                  WORK(KSCR1),WORK(KFOCKD),
     &                  .FALSE.,DUMMY,.FALSE.)

*---------------------------------------------------------------------*
*     Now we split:
*       for IOPTRES < 5 we compute an effective rhs vector
*       for IOPTRES = 5 we compute contractions L B T^A T^B
*---------------------------------------------------------------------*
      IF (IOPTRES.GE.1 .AND. IOPTRES.LE.4) THEN

         CALL DCOPY(NT1AMX,OMEGA1,1,OMEGA1EFF,1)
         CALL DCOPY(NT2AMX,OMEGA1,1,OMEGA2EFF,1)

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
