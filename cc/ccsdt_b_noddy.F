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
      REAL FREQA, FREQB, FREQ, DDOT, TCON, FREQC
      REAL DOTPROD(MXVEC,NBTRAN), WORK(LWORK), SIXTH, TWO
      REAL XLAMDP0(*), XLAMDH0(*)
      REAL OMEGA1(*), OMEGA2(*)
      REAL OMEGA1EFF(*), OMEGA2EFF(*)
#else
      DOUBLE PRECISION FREQA, FREQB, FREQ, DDOT, TCON, FREQC
      DOUBLE PRECISION DOTPROD(MXVEC,NBTRAN), WORK(LWORK), SIXTH, TWO
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*)
      DOUBLE PRECISION OMEGA1(*), OMEGA2(*)
      DOUBLE PRECISION OMEGA1EFF(*), OMEGA2EFF(*)
#endif
      PARAMETER( SIXTH = 1.0D0/6.0D0, TWO = 2.0D0 )

      CHARACTER*10 MODEL
      INTEGER KSCR1, KFOCKD, KEND1, KINT1T0, KINT2T0, KINT1TA, KINT2TA,
     &        KINT1TB, KINT2TB, KINT1SA, KINT2SA, KINT1SB, KINT2SB, 
     &        KINT1SAB, KINT2SAB, KTA3AM, KTB3AM, KFOCKA, KFOCKB,
     &        KLAMHB, KLAMPB, KLAMHA, KLAMPA, KT2AM, KDUM, KXIAJB,
     &        KYIAJB, KFOCK0, KOMEGA1, KOMEGA2, LWRK1, KXINT, KEND2,
     &        LWRK2, KB3AM, KLC3AM, KRECNR, KL1AM, KL2AM
      INTEGER IJ, NIJ, LUSIFC, INDEX, ILSTSYM, IOPT, LUFOCK,
     &        ILLL, IDEL, ISYDIS, ISYMA, ISYMB, IVEC, ISYMC, IDLSTC,
     &        ISYMD

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

      KFOCKA  = KEND1
      KLAMPA  = KFOCKA + NORBT*NORBT
      KLAMHA  = KLAMPA + NLAMDT
      KEND1   = KLAMHA + NLAMDT

      KFOCKB  = KEND1
      KLAMPB  = KFOCKB + NORBT*NORBT
      KLAMHB  = KLAMPB + NLAMDT
      KEND1   = KLAMHB + NLAMDT

      KINT1T0 = KEND1
      KINT2T0 = KINT1T0 + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2T0 + NRHFT*NRHFT*NT1AMX

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

      KXIAJB   = KEND1
      KYIAJB   = KXIAJB + NT1AMX*NT1AMX
      KEND1    = KYIAJB + NT1AMX*NT1AMX

      KTA3AM   = KEND1
      KEND1    = KTA3AM + NT1AMX*NT1AMX*NT1AMX

      KTB3AM   = KEND1
      KEND1    = KTB3AM + NT1AMX*NT1AMX*NT1AMX

      KOMEGA1  = KEND1
      KOMEGA2  = KOMEGA1 + NT1AMX
      KEND1    = KOMEGA2 + NT1AMX*NT1AMX

      KT2AM  = KYIAJB
      KB3AM  = KTA3AM
      KLC3AM = KTB3AM
      KL2AM  = KT2AM
  
      IF (LISTDP(1:3).EQ.'L0 ') THEN
        KL1AM = KEND1
        KEND1 = KL1AM + NT1AMX
      END IF

      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_BMAT_NODDY')
      ENDIF

*---------------------------------------------------------------------*
*     Read SCF orbital energies from file:
*---------------------------------------------------------------------*
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
*     Compute corrections to triples vector T^A_3, and corresponding 
*     lambda matrices and the XINT1SA,XINT2SA integrals and set FREQA:
*---------------------------------------------------------------------*
      IF (LISTA(1:3).EQ.'R1 ' .or. LISTA(1:3).EQ.'RE ') THEN
        CALL CCSDT_T31_NODDY(WORK(KTA3AM),LISTA,IDLSTA,FREQA,.FALSE.,
     &                       .TRUE. ,WORK(KINT1SA),WORK(KINT2SA),
     &                       .FALSE.,WORK(KDUM),WORK(KDUM),
     &                       .FALSE.,WORK(KDUM),WORK(KDUM),
     &                               WORK(KINT1SA),WORK(KINT2SA),
     &                       WORK(KLAMPA),WORK(KLAMHA),WORK(KFOCKA),
     &                       XLAMDP0,XLAMDH0,WORK(KFOCK0),
     &                       WORK(KDUM),WORK(KFOCKD),
     &                       WORK(KEND1),LWRK1)
      ELSE
        CALL QUIT('Unkown LISTB in CCSDT_BMAT_NODDY.')
      END IF

*---------------------------------------------------------------------*
*     Compute corrections to triples vector T^B_3, and corresponding 
*     lambda matrices and the XINT1SB,XINT2SB integrals and set FREQB:
*---------------------------------------------------------------------*
      IF (LISTB(1:3).EQ.'R1 ' .or. LISTB(1:3).EQ.'RE ') THEN
        CALL CCSDT_T31_NODDY(WORK(KTB3AM),LISTB,IDLSTB,FREQB,.FALSE.,
     &                       .TRUE. ,WORK(KINT1SB),WORK(KINT2SB),
     &                       .FALSE.,WORK(KDUM),WORK(KDUM),
     &                       .FALSE.,WORK(KDUM),WORK(KDUM),
     &                               WORK(KINT1SB),WORK(KINT2SB),
     &                       WORK(KLAMPB),WORK(KLAMHB),WORK(KFOCKB),
     &                       XLAMDP0,XLAMDH0,WORK(KFOCK0),
     &                       WORK(KDUM),WORK(KFOCKD),
     &                       WORK(KEND1),LWRK1)
      ELSE
        CALL QUIT('Unkown LISTB in CCSDT_BMAT_NODDY.')
      END IF

*---------------------------------------------------------------------*
*     Compute required integrals:
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KINT1SAB),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2SAB),NT1AMX*NRHFT*NRHFT)

      CALL DZERO(WORK(KINT1T0),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2T0),NT1AMX*NRHFT*NRHFT)

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
     *                  WORK(KRECNR),DIRECT)
 
C           ----------------------------------
C           Calculate integrals needed in CC3:
C           ----------------------------------
            CALL CC3_TRAN2(WORK(KXIAJB),WORK(KYIAJB),XLAMDP0,XLAMDH0,
     &                     WORK(KXINT),IDEL)

            CALL CCSDT_TRAN1_R(WORK(KINT1T0),WORK(KINT2T0),
     &                         XLAMDP0,XLAMDH0,
     &                         XLAMDP0,XLAMDH0,
     &                         WORK(KXINT),IDEL)

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

            ! XINT1SAB = XINT1SAB + (C-barB K|B-barC D)
            ! XINT2SAB = XINT2SAB + (C-barB K|L J-barC)
            CALL CCSDT_TRAN3_R(WORK(KINT1SAB),WORK(KINT2SAB),
     &                         XLAMDP0,XLAMDH0,
     &                         WORK(KLAMPB),XLAMDH0,
     &                         WORK(KLAMPA),WORK(KLAMHA),
     &                         WORK(KXINT),IDEL)

            ! XINT1SAB = XINT1SAB + (C-barC K-barB|B D)
            ! XINT2SAB = XINT2SAB + (C-barC K-barB|L J)
            CALL CCSDT_TRAN3_R(WORK(KINT1SAB),WORK(KINT2SAB),
     &                         XLAMDP0,XLAMDH0,
     &                         WORK(KLAMPA),WORK(KLAMHB),
     &                         XLAMDP0,XLAMDH0,
     &                         WORK(KXINT),IDEL)

            ! XINT1SAB = XINT1SAB + (C K-barB|B-barC D)
            ! XINT2SAB = XINT2SAB + (C K-barB|L J-barC)
            CALL CCSDT_TRAN3_R(WORK(KINT1SAB),WORK(KINT2SAB),
     &                         XLAMDP0,XLAMDH0,
     &                         XLAMDP0,WORK(KLAMHB),
     &                         WORK(KLAMPA),WORK(KLAMHA),
     &                         WORK(KXINT),IDEL)

            ! XINT1SAB = XINT1SAB + (C-barC K|B-barB D)
            ! XINT2SAB = XINT2SAB + (C-barC K|L J-barB)
            CALL CCSDT_TRAN3_R(WORK(KINT1SAB),WORK(KINT2SAB),
     &                         XLAMDP0,XLAMDH0,
     &                         WORK(KLAMPA),XLAMDH0,
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KXINT),IDEL)

            ! XINT1SAB = XINT1SAB + (C K-barC|B-barB D)
            ! XINT2SAB = XINT2SAB + (C K-barC|L J-barB)
            CALL CCSDT_TRAN3_R(WORK(KINT1SAB),WORK(KINT2SAB),
     &                         XLAMDP0,XLAMDH0,
     &                         XLAMDP0,WORK(KLAMHA),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KXINT),IDEL)

         END DO   
      END DO  

*---------------------------------------------------------------------*
*     Compute corrections to doubles result vector from T^A_3, T^B_3:
*---------------------------------------------------------------------*
      CALL CCSDT_OMEGA2(WORK(KOMEGA2),WORK(KINT1TA),WORK(KINT2TA),
     &                  WORK(KTB3AM),WORK(KFOCK0))

      CALL CCSDT_OMEGA2(WORK(KOMEGA2),WORK(KINT1TB),WORK(KINT2TB),
     &                  WORK(KTA3AM),WORK(KFOCK0))

*---------------------------------------------------------------------*
*     Compute triples vector (B T^A T^B)_3:
*---------------------------------------------------------------------*

      CALL DZERO(WORK(KB3AM),NT1AMX*NT1AMX*NT1AMX)

      ! add <nu_3|[H^AB,T^0_2]|HF>
      IOPT = 2
      Call CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KDUM),WORK(KT2AM))
      CALL CCSDT_B3AM(WORK(KB3AM),WORK(KT2AM),
     &                WORK(KINT1SAB),WORK(KINT2SAB))

      ! add <nu_3|[H^A,T^B_2]|HF>
      ISYMB = ILSTSYM(LISTB,IDLSTB)
      IOPT  = 2
      Call CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &              WORK(KDUM),WORK(KT2AM))
      Call CCLR_DIASCL(WORK(KT2AM),TWO,ISYMB) 
      CALL CCSDT_B3AM(WORK(KB3AM),WORK(KT2AM),
     &                WORK(KINT1SA),WORK(KINT2SA))

      ! add <nu_3|[H^B,T^A_2]|HF>
      ISYMA = ILSTSYM(LISTA,IDLSTA)
      IOPT  = 2
      Call CC_RDRSP(LISTA,IDLSTA,ISYMA,IOPT,MODEL,
     &              WORK(KDUM),WORK(KT2AM))
      Call CCLR_DIASCL(WORK(KT2AM),TWO,ISYMA) 
      CALL CCSDT_B3AM(WORK(KB3AM),WORK(KT2AM),
     &                WORK(KINT1SB),WORK(KINT2SB))

*---------------------------------------------------------------------*
*     Now we split:
*       for IOPTRES < 5 we compute an effective rhs vector
*       for IOPTRES = 5 we compute contractions L B T^A T^B
*---------------------------------------------------------------------*
      IF (IOPTRES.GE.1 .AND. IOPTRES.LE.4) THEN

C        -----------------------------------------------------------
C        Calculate -B_nu_3 = Bi_nu_3 / (eps_nu_3 - w)
C        -----------------------------------------------------------
         IF (LWRK1.LT.NT1AMX) CALL QUIT('OUT OF MEMORY - CCSDT_XI3')

         FREQ = FREQA + FREQB

         CALL CCSDT_3AM(WORK(KB3AM),FREQ,WORK(KEND1),WORK(KFOCKD))
         CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KB3AM),1)

C        --------------------------------------------------------
C        Add contribution of xksi_3 to the effective xksi vector:
C        --------------------------------------------------------
         CALL DZERO(WORK(KOMEGA2),NT1AMX*NT1AMX)
         CALL DZERO(WORK(KOMEGA1),NT1AMX)

         CALL CCSDT_OMEGA1(WORK(KOMEGA1),WORK(KXIAJB),WORK(KB3AM))
 
         CALL CCSDT_OMEGA2(WORK(KOMEGA2),WORK(KINT1T0),WORK(KINT2T0),
     &                     WORK(KB3AM),WORK(KFOCK0))
 
         DO I = 1,NT1AMX
            OMEGA1EFF(I) = OMEGA1(I) + WORK(KOMEGA1+I-1)
         END DO
 
         DO I = 1,NT1AMX
            DO J = 1,I
               IJ = NT1AMX*(I-1) + J
               NIJ = INDEX(I,J)
               OMEGA2EFF(NIJ) = OMEGA2(NIJ) + WORK(KOMEGA2+IJ-1)
            END DO
         END DO

         IF (LOCDBG) THEN
           WRITE(LUPRI,*), 'CCSDT_BMAT_NODDY> RESULT VECTOR ON OUTPUT:'
           WRITE(LUPRI,*), 'NORM^2 = ',
     &            DDOT(NT1AMX+NT2AMX,OMEGA1,1,OMEGA1,1)
           ! CALL CC_PRP(OMEGA1,OMEGA2,1,1,1)
           WRITE(LUPRI,*), 'CCSDT_BMAT_NODDY> EFFECTIVE RESULT VECTOR:'
           WRITE(LUPRI,*), 'NORM^2 = ',
     &            DDOT(NT1AMX+NT2AMX,OMEGA1EFF,1,OMEGA1EFF,1)
           ! CALL CC_PRP(OMEGA1EFF,OMEGA2EFF,1,1,1)
         END IF

      ELSE IF (IOPTRES.EQ.5) THEN
*---------------------------------------------------------------------*
*       for IOPTRES = 5 we compute the contractions : L B T^A T^B
*       here in the triples code we compute only the contribution 
*       from the triples result vector: 
*---------------------------------------------------------------------*
        IVEC = 1
        DO WHILE (IDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)

          IDLSTC = IDOTS(IVEC,ITRAN)
          ISYMC  = ILSTSYM(LISTDP,IDLSTC)


          IF (LISTDP(1:3).EQ.'L0 ') THEN

             CALL DZERO(WORK(KLC3AM),NT1AMX*NT1AMX*NT1AMX)

             CALL CCSDT_L3AM(WORK(KLC3AM),WORK(KINT1T0),WORK(KINT2T0),
     *                       WORK(KXIAJB),WORK(KFOCK0),WORK(KL1AM),
     *                       WORK(KL2AM),WORK(KSCR1),WORK(KFOCKD),0.0D0)

             CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KLC3AM),1)


      
          ELSE IF (LISTDP(1:3).EQ.'L1 ' .OR. LISTDP(1:3).EQ.'LE ') THEN

             CALL CCSDT_TBAR31_NODDY(WORK(KLC3AM),FREQC,LISTDP,IDLSTC,
     &                             XLAMDP0,XLAMDH0,
     &                             WORK(KFOCK0),WORK(KFOCKD),
     &                             WORK(KSCR1),
     &                             .FALSE.,WORK(KXIAJB),DUMMY,
     &                             .FALSE.,WORK(KINT1T0),WORK(KINT2T0),
     &                             WORK(KEND1),LWRK1)
          ELSE
            CALL QUIT('Unknown or illegal list in CCSDT_BMAT_NODDY.')
          END IF

          IF (LOCDBG) THEN
            WRITE(LUPRI,*) 'CCSDT_BMAT_NODDY> IVEC,ITRAN:',IVEC,ITRAN
            WRITE(LUPRI,*) 'NORM^2(B3AM):',
     &        DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KB3AM),1,WORK(KB3AM),1)
            WRITE(LUPRI,*) 'NORM^2(TC3AM):',
     &        DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KLC3AM),1,WORK(KLC3AM),1)
          END IF

          TCON = SIXTH *
     &      DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KLC3AM),1,WORK(KB3AM),1)

          DOTPROD(IVEC,ITRAN) = DOTPROD(IVEC,ITRAN) + TCON

          IF (LOCDBG) THEN
            WRITE(LUPRI,*) 'CCSDT_BMAT_NODDY:'
            WRITE(LUPRI,*) 'L B T^A T^B:',TCON
          END IF


          IVEC = IVEC + 1

        END DO


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
*=====================================================================*
      SUBROUTINE CCSDT_B3AM(B3AM,T2AM,XINT1,XINT2)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples result to B matrix transformation
*
*      <nu_3| [H^AB,T^0_2] + [H^A,T^B_2] + [H^B,T^A_2] |HF>
*
*     to compute the full contributions, the routine has to be called
*     three times with the 3 different combinations of H and T_2
*
*     Written by Christof Haettig, April 2002
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
#include <ccorb.h>

#if defined (SYS_CRAY)
      REAL B3AM(NT1AMX,NT1AMX,NT1AMX), T2AM(NT1AMX,NT1AMX)
      REAL XINT1(NT1AMX,NVIRT,NVIRT)
      REAL XINT2(NT1AMX,NRHFT,NRHFT)
      REAL AIBJCK
#else
      DOUBLE PRECISION B3AM(NT1AMX,NT1AMX,NT1AMX), T2AM(NT1AMX,NT1AMX)
      DOUBLE PRECISION XINT1(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT2(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION AIBJCK
#endif
      INTEGER NCK, NBJ, NAI, NDJ, NBL

      DO NCK = 1,NT1AMX
         DO J = 1,NRHFT
            DO B = 1,NVIRT
               NBJ = NVIRT*(J-1) + B
               DO NAI = 1,NT1AMX
 
                  AIBJCK = 0.0D0
                  DO D = 1,NVIRT
                     NDJ = NVIRT*(J-1) + D
                     AIBJCK = AIBJCK + XINT1(NCK,B,D)*T2AM(NDJ,NAI)
                  END DO
 
                  DO L = 1,NRHFT
                     NBL = NVIRT*(L-1) + B
                     AIBJCK = AIBJCK - XINT2(NCK,L,J)*T2AM(NBL,NAI)
                  END DO

                  B3AM(NAI,NBJ,NCK) = B3AM(NAI,NBJ,NCK) + AIBJCK
                  B3AM(NAI,NCK,NBJ) = B3AM(NAI,NCK,NBJ) + AIBJCK
                  B3AM(NBJ,NAI,NCK) = B3AM(NBJ,NAI,NCK) + AIBJCK
                  B3AM(NCK,NAI,NBJ) = B3AM(NCK,NAI,NBJ) + AIBJCK
                  B3AM(NBJ,NCK,NAI) = B3AM(NBJ,NCK,NAI) + AIBJCK
                  B3AM(NCK,NBJ,NAI) = B3AM(NCK,NBJ,NAI) + AIBJCK
 
               END DO
            END DO
         END DO
      END DO

      CALL CCSDT_CLEAN_T3(B3AM,NT1AMX,NVIRT,NRHFT)

      RETURN
      END
*=====================================================================*
