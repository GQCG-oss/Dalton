
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
*=====================================================================*
      SUBROUTINE CCSDT_BAMAT_NODDY(IOPTRES,LABELB,FREQB,ISYMB,
     &                             LISTC,IDLSTC,
     &                             LISTD,IDLSTD,
     &                             OMEGA1,OMEGA2,
     &                             OMEGA1EFF,OMEGA2EFF,
     &                             IDOTS,DOTPROD,LISTDP,ITRAN,
     &                             NXTRAN,MXVEC,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples contribution to A{B} transformed vector
*
*    (A{B} T^C)^eff_1,2 = (A{B} T^C)_1,2(CCSD) + (A{B} T^C)_1,2(T3)
*                            - A_1,2;3 (w_3 - w)^1 (A{B} T^C)_3
*
*        
*   Written by Christof Haettig, April 2002, based on CCSDT_XI_NODDY
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

      CHARACTER LISTDP*3, LABELB*8, LISTC*3, LISTD*3
      INTEGER LWORK, IOPTRES, ITRAN, MXVEC, NXTRAN,ISYMB,IDLSTC,IDLSTD
      INTEGER IDOTS(MXVEC,NXTRAN)

#if defined (SYS_CRAY)
      REAL DOTPROD(MXVEC,NXTRAN), FREQC, FREQE, FREQ, FREQB
      REAL WORK(LWORK), SIXTH, TWO, DUMMY, TCON, DDOT
      REAL OMEGA1(*), OMEGA2(*), FREQD
      REAL OMEGA1EFF(*), OMEGA2EFF(*)
#else
      DOUBLE PRECISION DOTPROD(MXVEC,NXTRAN), FREQC, FREQE, FREQ, FREQB
      DOUBLE PRECISION WORK(LWORK), SIXTH, TWO, DUMMY, TCON, DDOT
      DOUBLE PRECISION OMEGA1(*), OMEGA2(*), FREQD
      DOUBLE PRECISION OMEGA1EFF(*), OMEGA2EFF(*)
#endif
      PARAMETER( SIXTH = 1.0D0/6.0D0, TWO = 2.0D0 )

      CHARACTER*10 MODEL
      INTEGER KFOCKB, KT3AM, KT2AMD, KT2AMC, KEND2, KEND1, KOMEGA1,
     &        KOMEGA2, KOMEGA3, KSCR1, KFOCKD, KLAMP0, KLAMH0, KFOCK0,
     &        KFOCKB_AO, KFOCKBC, KLAMPC, KLAMHC, KPERTC, KINT1SC, 
     &        KINT2SC, KINT1SD, KINT2SD, KXIAJB, KYIAJB, KINT1T, KINT2T,
     &        LWRK1, KDUM, KXINT, LWRK2, KRECNR, KL3AM, KT1AMP0,
     &        KFOCKBD, KTC1AM, KLAMPD, KLAMHD, KPERTD, KTD1AM
      INTEGER IJ, NIJ, LUSIFC, INDEX, IDUMMY, ILSTSYM, ISYMC, LUFOCK, 
     &        IRREP, IERR, ILLL, IDEL, ISYDIS, IOPT, ISYMD, IVEC,
     &        IDLSTE

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J 

      CALL QENTER('CCSDT_BAMAT_NODDY')

      IF(DIRECT)CALL QUIT('DIRECT NOT IMPLEMENTED IN CCSDT_BAMAT_NODDY')

*---------------------------------------------------------------------*
*     Memory allocation:
*---------------------------------------------------------------------*

      KOMEGA1 = KEND1
      KOMEGA2 = KOMEGA1 + NT1AMX
      KOMEGA3 = KOMEGA2 + NT1AMX*NT1AMX
      KEND1   = KOMEGA3 + NT1AMX*NT1AMX*NT1AMX

      KSCR1   = KEND1 
      KFOCKD  = KSCR1  + NT1AMX
      KLAMP0  = KSCR1  + NORBT
      KLAMH0  = KLAMP0 + NLAMDT
      KFOCK0  = KLAMH0 + NLAMDT
      KT1AMP0 = KFOCK0 + NORBT*NORBT
      KEND1   = KT1AMP0+ NT1AMX

      KFOCKB    = KEND1 
      KFOCKB_AO = KFOCKB    + NORBT*NORBT
      KFOCKBC   = KFOCKB_AO + NORBT*NORBT
      KFOCKBD   = KFOCKBC   + NORBT*NORBT
      KEND1     = KFOCKBD   + NORBT*NORBT

      KLAMPC  = KEND1
      KLAMHC  = KLAMPC + NLAMDT
      KPERTC  = KLAMHC + NLAMDT
      KTC1AM  = KPERTC + NORBT*NORBT
      KEND1   = KTC1AM + NT1AMX

      KLAMPD  = KEND1
      KLAMHD  = KLAMPD + NLAMDT
      KPERTD  = KLAMHD + NLAMDT
      KTD1AM  = KPERTD + NORBT*NORBT
      KEND1   = KTD1AM + NT1AMX

      KINT1SC = KEND1
      KINT2SC = KINT1SC + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2SC + NRHFT*NRHFT*NT1AMX 

      KINT1SD = KINT1SC
      KINT2SD = KINT2SC

      KXIAJB  = KEND1
      KYIAJB  = KXIAJB  + NT1AMX*NT1AMX
      KEND1   = KYIAJB  + NT1AMX*NT1AMX

      KINT1T = KEND1
      KINT2T = KINT1T + NT1AMX*NVIRT*NVIRT
      KEND1  = KINT2T + NRHFT*NRHFT*NT1AMX 

      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_BAMAT_NODDY')
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
*     Get zeroth-order Lambda matrices:
*---------------------------------------------------------------------*
      IOPT   = 1
      Call CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AMP0),WORK(KDUM))

      Call LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AMP0),
     &            WORK(KEND1),LWRK1)

*---------------------------------------------------------------------*
*     Calculate response Lambda matrices:
*---------------------------------------------------------------------*
      ISYMC = ILSTSYM(LISTC,IDLSTC)
      IOPT = 1
      CALL CC_RDRSP(LISTC,IDLSTC,ISYMC,IOPT,MODEL,WORK(KTC1AM),DUMMY)

      CALL CCLR_LAMTRA(WORK(KLAMP0),WORK(KLAMPC),
     &                 WORK(KLAMH0),WORK(KLAMHC),WORK(KTC1AM),ISYMC)

      ISYMD = ILSTSYM(LISTD,IDLSTD)
      IOPT = 2
      CALL CC_RDRSP(LISTD,IDLSTD,ISYMD,IOPT,MODEL,WORK(KTD1AM),DUMMY)

      CALL CCLR_LAMTRA(WORK(KLAMP0),WORK(KLAMPD),
     &                 WORK(KLAMH0),WORK(KLAMHD),WORK(KTD1AM),ISYMD)

*---------------------------------------------------------------------*
*     read zeroth-order AO Fock matrix from file: 
*---------------------------------------------------------------------*
      LUFOCK = 0
      CALL GPOPEN(LUFOCK,'CC_FCKH','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUFOCK)
      READ(LUFOCK) (WORK(KFOCK0-1+I),I=1,N2BST(ISYM0))
      CALL GPCLOSE(LUFOCK,'KEEP')

      CALL CC_FCKMO(WORK(KFOCK0),WORK(KLAMP0),WORK(KLAMH0),
     &              WORK(KEND1),LWRK1,ISYM0,ISYM0,ISYM0)

*---------------------------------------------------------------------*
*     Matrix with property integrals in MO basis:
*---------------------------------------------------------------------*
      ! read property integrals from file:
      CALL CCPRPAO(LABELB,WORK(KFOCKB_AO),IRREP,ISYMB,IERR,
     &             WORK(KEND1),LWRK1)
      CALL DCOPY(NORBT*NORBT,WORK(KFOCKB_AO),1,WORK(KFOCKB),1)
      IF ((IERR.GT.0) .OR. (IERR.EQ.0 .AND. IRREP.NE.ISYMB)) THEN
        CALL QUIT('CCSDT_BAMAT_NODDY: error reading operator '//LABELB)
      ELSE IF (IERR.LT.0) THEN
        CALL DZERO(WORK(KFOCKB),N2BST(ISYMB))
      END IF
 
      ! transform property integrals to Lambda-MO basis
      CALL CC_FCKMO(WORK(KFOCKB),WORK(KLAMP0),WORK(KLAMH0),
     &              WORK(KEND1),LWRK1,ISYMB,1,1)

*---------------------------------------------------------------------*
*     Loop over distributions of integrals:
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KINT1T),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2T),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KXIAJB),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KYIAJB),NT1AMX*NT1AMX)

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
               CALL QUIT('Insufficient space in CCSDT_XI3_NODDY')
            ENDIF
 
C           ---------------------------
C           Read in batch of integrals.
C           ---------------------------
            CALL CCRDAO(WORK(KXINT),IDEL,1,WORK(KEND2),LWRK2,
     *                  WORK(KRECNR),DIRECT)
 
C           ----------------------------------
C           Calculate integrals needed in CC3:
C           ----------------------------------
            CALL CCSDT_TRAN1(WORK(KINT1T),WORK(KINT2T),
     &                       WORK(KLAMP0),WORK(KLAMH0),WORK(KXINT),IDEL)

            CALL CC3_TRAN2(WORK(KXIAJB),WORK(KYIAJB),WORK(KLAMP0),
     &                     WORK(KLAMH0),WORK(KXINT),IDEL)
 
         END DO   
      END DO  

*---------------------------------------------------------------------*
*     End Loop over distributions of integrals.
*---------------------------------------------------------------------*
      KT3AM  = KEND1
      KT2AMD = KT3AM  + NT1AMX*NT1AMX*NT1AMX
      KT2AMC = KT2AMD + NT1AMX*NT1AMX
      KEND2  = KT2AMC + NT1AMX*NT1AMX
      LWRK2  = LWORK  - KEND2
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_BAMAT_NODDY')
      ENDIF

      ISYMC = ILSTSYM(LISTC,IDLSTC)
      IOPT = 2
      CALL CC_RDRSP(LISTC,IDLSTC,ISYMC,IOPT,MODEL,DUMMY,WORK(KT3AM))
      CALL CC_T2SQ(WORK(KT3AM),WORK(KT2AMC),ISYMC)

      ISYMD = ILSTSYM(LISTD,IDLSTD)
      IOPT = 2
      CALL CC_RDRSP(LISTD,IDLSTD,ISYMD,IOPT,MODEL,DUMMY,WORK(KT3AM))
      CALL CC_T2SQ(WORK(KT3AM),WORK(KT2AMD),ISYMD)

C     ---------------------------------
C     Initialize triples result vector:
C     ---------------------------------
      CALL DZERO(WORK(KOMEGA3),NT1AMX*NT1AMX*NT1AMX)

C     ---------------------------------------------------
C     add triples contribution: <mu_3|[[B,T2^D],T2^C]|HF>
C     ---------------------------------------------------
      CALL CCSDT_XKSI3_1(WORK(KOMEGA3),WORK(KFOCKB),
     &                   WORK(KT2AMD),WORK(KT2AMC),TWO)

C     -------------------------------------------------
C     calculate the first-order triples amplitudes T^C:
C     -------------------------------------------------
      CALL CCSDT_T31_NODDY(WORK(KT3AM),LISTC,IDLSTC,FREQC,.FALSE.,
     &                     WORK(KINT1SC),WORK(KINT2SC),
     &                     .FALSE.,WORK(KDUM),WORK(KDUM),
     &                     .FALSE.,WORK(KDUM),WORK(KDUM),
     &                     WORK(KLAMPC),WORK(KLAMHC),WORK(KPERTC),
     &                     WORK(KLAMP0),WORK(KLAMH0),WORK(KFOCK0),
     &                     WORK(KDUM),WORK(KFOCKD),
     &                     WORK(KEND2),LWRK2)

C     --------------------------------------------------------------
C     calculate one-index transf. property integrals: B^D = [B,T1^D]
C     --------------------------------------------------------------
      CALL DCOPY(NORBT*NORBT,WORK(KFOCKB_AO),1,WORK(KFOCKBD),1)
      CALL CC_FCKMO(WORK(KFOCKBD),WORK(KLAMPD),WORK(KLAMHD),
     &              WORK(KEND2),LWRK2,ISYMB,ISYMD,ISYMD)

C     ----------------------------------------------
C     add triples contribution: <mu_3|[B^D,T3^C]|HF>
C     ----------------------------------------------
      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KT3AM))
      CALL CCSDT_XKSI3_2(WORK(KOMEGA3),WORK(KFOCKBD),WORK(KT3AM))

C     -------------------------------------------------
C     calculate the first-order triples amplitudes T^D:
C     -------------------------------------------------
      CALL CCSDT_T31_NODDY(WORK(KT3AM),LISTD,IDLSTD,FREQD,.FALSE.,
     &                     WORK(KINT1SD),WORK(KINT2SD),
     &                     .FALSE.,WORK(KDUM),WORK(KDUM),
     &                     .FALSE.,WORK(KDUM),WORK(KDUM),
     &                     WORK(KLAMPD),WORK(KLAMHD),WORK(KPERTD),
     &                     WORK(KLAMP0),WORK(KLAMH0),WORK(KFOCK0),
     &                     WORK(KDUM),WORK(KFOCKD),
     &                     WORK(KEND2),LWRK2)

C     --------------------------------------------------------------
C     calculate one-index transf. property integrals: B^C = [B,T1^C]
C     --------------------------------------------------------------
      CALL DCOPY(NORBT*NORBT,WORK(KFOCKB_AO),1,WORK(KFOCKBC),1)
      CALL CC_FCKMO(WORK(KFOCKBC),WORK(KLAMPC),WORK(KLAMHC),
     &              WORK(KEND2),LWRK2,ISYMB,ISYMC,ISYMC)

C     ----------------------------------------------
C     add triples contribution: <mu_3|[B^C,T3^D]|HF>
C     ----------------------------------------------
      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KT3AM))
      CALL CCSDT_XKSI3_2(WORK(KOMEGA3),WORK(KFOCKBC),WORK(KT3AM))

*---------------------------------------------------------------------*
*     Now we split:
*       for IOPTRES < 5 we compute the result vector
*       for IOPTRES = 5 we compute the contractions Tbar^E B{B} T^C T^D
*---------------------------------------------------------------------*
      IF (IOPTRES.GE.1 .AND. IOPTRES.LE.4) THEN

C        -----------------------------------------------------------
C        Calculate -O_nu_3 = O_nu_3 / (eps_nu_3 - w)
C        -----------------------------------------------------------
         FREQ = FREQB + FREQC + FREQD

         IF (LWRK1.LT.NT1AMX) CALL QUIT('OUT OF MEMORY - CCSDT_BAMAT')
         CALL CCSDT_3AM(WORK(KOMEGA3),FREQ,WORK(KEND1),WORK(KFOCKD))

C        --------------------------------------------------------
C        Add contributions to the effective result vector:
C        --------------------------------------------------------
         CALL DZERO(WORK(KOMEGA2),NT1AMX*NT1AMX)
         CALL DZERO(WORK(KOMEGA1),NT1AMX)

         CALL CCSDT_OMEGA1(WORK(KOMEGA1),WORK(KXIAJB),WORK(KOMEGA3))
 
         CALL CCSDT_OMEGA2(WORK(KOMEGA2),WORK(KINT1T),WORK(KINT2T),
     *                     WORK(KOMEGA3),WORK(KFOCK0))
 
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

      ELSE IF (IOPTRES.EQ.5) THEN
*---------------------------------------------------------------------*
*       for IOPTRES = 5 we compute the contractions L^E B{B} T^C T^D:
*       here in the triples code we compute only the contribution 
*       from the triples result vector: 
*---------------------------------------------------------------------*
        KT3AM = KEND1
        KEND2 = KT3AM + NT1AMX*NT1AMX*NT1AMX
        LWRK2 = LWORK - KEND2
        IF (LWRK2 .LT. 0) THEN
           CALL QUIT('Insufficient space in CCSDT_BAMAT_NODDY')
        ENDIF

        IVEC = 1
        DO WHILE (IDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)

          IDLSTE = IDOTS(IVEC,ITRAN)
          ISYMD  = ILSTSYM(LISTDP,IDLSTE)

          IF (LISTDP(1:3).EQ.'L1 ' .OR. LISTDP(1:3).EQ.'LE ') THEN

             CALL CCSDT_TBAR31_NODDY(WORK(KL3AM),FREQE,LISTDP,IDLSTE,
     &                               WORK(KLAMP0),WORK(KLAMH0),
     &                               WORK(KFOCK0),WORK(KFOCKD),
     &                               WORK(KSCR1),
     &                               WORK(KXIAJB),WORK(KINT1T),
     &                               WORK(KINT2T),WORK(KEND2),LWRK2)
          ELSE
            CALL QUIT('Unknown or illegal list in CCSDT_BAMAT_NODDY.')
          END IF

          IF (LOCDBG) THEN
            WRITE(LUPRI,*) 'CCSDT_BAMAT_NODDY> IVEC,ITRAN:',IVEC,ITRAN
            WRITE(LUPRI,*) 'NORM^2(OMEGA3):',
     &        DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KOMEGA3),1,WORK(KOMEGA3),1)
            WRITE(LUPRI,*) 'NORM^2(L3AM):',
     &        DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KL3AM),1,WORK(KL3AM),1)
          END IF

          TCON = - SIXTH *
     &      DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KL3AM),1,WORK(KOMEGA3),1)

          DOTPROD(IVEC,ITRAN) = DOTPROD(IVEC,ITRAN) + TCON

          IF (LOCDBG) THEN
            WRITE(LUPRI,*) 'CCSDT_BAMAT_NODDY:'
            WRITE(LUPRI,*) 'L_3  x (B{B} T^C T^D)_3:',TCON
          END IF

          IVEC = IVEC + 1

        END DO


      ELSE
        CALL QUIT('Illegal value for IOPTRES IN CCSDT_BAMAT_NODDY')
      END IF

      CALL QEXIT('CCSDT_BAMAT_NODDY')
      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_BAMAT_NODDY                    *
*---------------------------------------------------------------------*
