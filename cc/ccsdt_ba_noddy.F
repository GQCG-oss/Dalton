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
      SUBROUTINE CCSDT_BAMAT_NODDY(IOPTRES,FREQRES,LABELB,ISYMB,
     &                             LISTC,IDLSTC,
     &                             LISTD,IDLSTD,
     &                             OMEGA1,OMEGA2,
     &                             OMEGA1EFF,OMEGA2EFF,
     &                             IDOTS,DOTPROD,LISTDP,ITRAN,
     &                             NXTRAN,MXVEC,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples contribution to B{A} transformed vector
*
*    (B{A} T^C T^D)^eff_1,2 = (B{A} T^C T^D)_1,2(CCSD) 
*                            - A_1,2;3 (w_3 - w_res)^1 (B{A} T^C T^D)_3
*
*        
*   Written by Christof Haettig, April 2002, based on CCSDT_AAMAT_NODDY
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
      REAL DOTPROD(MXVEC,NXTRAN), FREQC, FREQE, FREQ, FREQRES
      REAL WORK(LWORK), SIXTH, ONE, TWO, DUMMY, TCON, DDOT
      REAL OMEGA1(*), OMEGA2(*), FREQD, FF, SIGN
      REAL OMEGA1EFF(*), OMEGA2EFF(*)
#else
      DOUBLE PRECISION DOTPROD(MXVEC,NXTRAN), FREQC, FREQE, FREQRES
      DOUBLE PRECISION WORK(LWORK), SIXTH, ONE, TWO, DUMMY, TCON, DDOT
      DOUBLE PRECISION OMEGA1(*), OMEGA2(*), FREQD, FF, SIGN
      DOUBLE PRECISION OMEGA1EFF(*), OMEGA2EFF(*)
#endif
      PARAMETER( SIXTH = 1.0D0/6.0D0, ONE = 1.0D0, TWO = 2.0D0 )

      CHARACTER*10 MODEL
      INTEGER KFOCKB, KT3AM, KT2AMD, KT2AMC, KEND2, KEND1, KOMEGA1,
     &        KOMEGA2, KOMEGA3, KSCR1, KFOCKD, KLAMP0, KLAMH0, KFOCK0,
     &        KFOCKB_AO, KFOCKBC, KLAMPC, KLAMHC, KPERTC, KINT1SC, 
     &        KINT2SC, KINT1SD, KINT2SD, KXIAJB, KYIAJB, KINT1T, KINT2T,
     &        LWRK1, KDUM, KXINT, LWRK2, KL3AM, KT1AMP0,
     &        KFOCKBD, KTC1AM, KLAMPD, KLAMHD, KPERTD, KTD1AM
      INTEGER IJ, NIJ, LUSIFC, INDEX, IDUMMY, ILSTSYM, ISYMC, LUFOCK, 
     &        IRREP, IERR, ILLL, IDEL, ISYDIS, IOPT, ISYMD, IVEC,
     &        IDLSTE, KFCKBUF, KINT1S, KINT2S, KFIELD, KFIELDAO

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J 

      CALL QENTER('CCSDT_BAMAT_NODDY')

      IF(DIRECT)CALL QUIT('DIRECT NOT IMPLEMENTED IN CCSDT_BAMAT_NODDY')

*---------------------------------------------------------------------*
*     Memory allocation:
*---------------------------------------------------------------------*
      KEND1   = 1

      KOMEGA1 = KEND1
      KOMEGA2 = KOMEGA1 + NT1AMX
      KOMEGA3 = KOMEGA2 + NT1AMX*NT1AMX
      KEND1   = KOMEGA3 + NT1AMX*NT1AMX*NT1AMX

      KSCR1   = KEND1 
      KFOCKD  = KSCR1  + NT1AMX
      KLAMP0  = KFOCKD + NORBT
      KLAMH0  = KLAMP0 + NLAMDT
      KFOCK0  = KLAMH0 + NLAMDT
      KT1AMP0 = KFOCK0 + NORBT*NORBT
      KEND1   = KT1AMP0+ NT1AMX
 
      IF (NONHF) THEN
        KFIELD   = KEND1
        KFIELDAO = KFIELD   + NORBT*NORBT
        KEND1    = KFIELDAO + NORBT*NORBT
      END IF 

      KFOCKB    = KEND1 
      KFOCKB_AO = KFOCKB    + NORBT*NORBT
      KFOCKBC   = KFOCKB_AO + NORBT*NORBT
      KFOCKBD   = KFOCKBC   + NORBT*NORBT
      KFCKBUF   = KFOCKBD   + NORBT*NORBT
      KEND1     = KFCKBUF   + NORBT*NORBT

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

      KINT1S = KEND1
      KINT2S = KINT1S + NT1AMX*NVIRT*NVIRT
      KEND1  = KINT2S + NRHFT*NRHFT*NT1AMX 

      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_BAMAT_NODDY')
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
*     Calculate response Lambda matrices:
*---------------------------------------------------------------------*
      ISYMC = ILSTSYM(LISTC,IDLSTC)
      IOPT = 1
      CALL CC_RDRSP(LISTC,IDLSTC,ISYMC,IOPT,MODEL,WORK(KTC1AM),DUMMY)

      CALL CCLR_LAMTRA(WORK(KLAMP0),WORK(KLAMPC),
     &                 WORK(KLAMH0),WORK(KLAMHC),WORK(KTC1AM),ISYMC)

      ISYMD = ILSTSYM(LISTD,IDLSTD)
      IOPT = 1
      CALL CC_RDRSP(LISTD,IDLSTD,ISYMD,IOPT,MODEL,WORK(KTD1AM),DUMMY)

      CALL CCLR_LAMTRA(WORK(KLAMP0),WORK(KLAMPD),
     &                 WORK(KLAMH0),WORK(KLAMHD),WORK(KTD1AM),ISYMD)

*---------------------------------------------------------------------*
*     read zeroth-order AO Fock matrix from file: 
*---------------------------------------------------------------------*
      LUFOCK = -1
      CALL GPOPEN(LUFOCK,'CC_FCKH','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUFOCK)
      READ(LUFOCK) (WORK(KFOCK0-1+I),I=1,N2BST(ISYM0))
      CALL GPCLOSE(LUFOCK,'KEEP')

      CALL CC_FCKMO(WORK(KFOCK0),WORK(KLAMP0),WORK(KLAMH0),
     &              WORK(KEND1),LWRK1,ISYM0,ISYM0,ISYM0)

*---------------------------------------------------------------------*
*     If needed get external field:
*---------------------------------------------------------------------*
      IF ((NONHF) .AND. NFIELD .GT. 0) THEN
         CALL DZERO(WORK(KFIELDAO),NORBT*NORBT)
         DO I = 1, NFIELD
          FF = EFIELD(I)
          CALL CC_ONEP(WORK(KFIELDAO),WORK(KEND1),LWRK1,FF,1,LFIELD(I))
         ENDDO
 
         CALL DCOPY(NORBT*NORBT,WORK(KFIELDAO),1,WORK(KFIELD),1)
         CALL CC_FCKMO(WORK(KFIELD),WORK(KLAMP0),WORK(KLAMH0),
     *                 WORK(KEND1),LWRK1,1,1,1)
      ENDIF

*---------------------------------------------------------------------*
*     Matrix with property integrals in MO basis:
*---------------------------------------------------------------------*
      ! read property integrals from file:
      CALL CCPRPAO(LABELB,.TRUE.,WORK(KFOCKB_AO),IRREP,ISYMB,IERR,
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
*     Compute some integrals:
*           XINT1T  =  (KC|BD)
*           XINT2T  =  (KC|LJ)
*           XIAJB   = 2(IA|JB) - (IB|JA)
*           YIAJB   =  (IA|JB)
*---------------------------------------------------------------------*
      CALL CCSDT_INTS0_NODDY(.TRUE.,WORK(KXIAJB),WORK(KYIAJB),
     &                       .TRUE.,WORK(KINT1S),WORK(KINT2S),
     &                       .TRUE.,WORK(KINT1T),WORK(KINT2T),
     &                       WORK(KLAMP0),WORK(KLAMH0),
     &                       WORK(KEND1),LWRK1)

*---------------------------------------------------------------------*
*     allocate work space for a triples and two doubles vectors and
*     read zeroth-order and response doubles amplitudes from file,
*     compute response Lambda matrices:
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
      IOPT  = 2
      CALL CC_RDRSP(LISTC,IDLSTC,ISYMC,IOPT,MODEL,DUMMY,WORK(KT3AM))
      Call CCLR_DIASCL(WORK(KT3AM),TWO,ISYMC)
      CALL CC_T2SQ(WORK(KT3AM),WORK(KT2AMC),ISYMC)

      ISYMD = ILSTSYM(LISTD,IDLSTD)
      IOPT  = 2
      CALL CC_RDRSP(LISTD,IDLSTD,ISYMD,IOPT,MODEL,DUMMY,WORK(KT3AM))
      Call CCLR_DIASCL(WORK(KT3AM),TWO,ISYMD)
      CALL CC_T2SQ(WORK(KT3AM),WORK(KT2AMD),ISYMD)

*---------------------------------------------------------------------*
*     compute the three contributions to the triples result vector:
*---------------------------------------------------------------------*
      ! initialize triples result vector:
      CALL DZERO(WORK(KOMEGA3),NT1AMX*NT1AMX*NT1AMX)

C     ---------------------------------------------------
C     add triples contribution: <mu_3|[[B,T2^D],T2^C]|HF>
C     ---------------------------------------------------
      CALL CCSDT_XKSI3_1(WORK(KOMEGA3),WORK(KFOCKB),
     &                   WORK(KT2AMD),WORK(KT2AMC),ONE)
      CALL CCSDT_XKSI3_1(WORK(KOMEGA3),WORK(KFOCKB),
     &                   WORK(KT2AMC),WORK(KT2AMD),ONE)

C     -------------------------------------------------
C     calculate the first-order triples amplitudes T^C:
C     -------------------------------------------------
      KDUM = KEND2
      CALL CCSDT_T31_NODDY(WORK(KT3AM),LISTC,IDLSTC,FREQC,.FALSE.,
     &                     .FALSE.,WORK(KINT1S),WORK(KINT2S),
     &                     .FALSE.,WORK(KDUM),WORK(KDUM),
     &                     .FALSE.,WORK(KDUM),WORK(KDUM),
     &                             WORK(KINT1SC),WORK(KINT2SC),
     &                     WORK(KLAMPC),WORK(KLAMHC),WORK(KPERTC),
     &                     WORK(KLAMP0),WORK(KLAMH0),WORK(KFOCK0),
     &                     WORK(KDUM),WORK(KFOCKD),
     &                     WORK(KEND2),LWRK2)

      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KT3AM),1)

C     --------------------------------------------------------------
C     calculate one-index transf. property integrals: B^D = [B,T1^D]
C     --------------------------------------------------------------
      CALL DCOPY(NORBT*NORBT,WORK(KFOCKB_AO),1,WORK(KFCKBUF),1)
      CALL DCOPY(NORBT*NORBT,WORK(KFOCKB_AO),1,WORK(KFOCKBD),1)
      CALL CC_FCKMO(WORK(KFCKBUF),WORK(KLAMPD),WORK(KLAMH0),
     &              WORK(KEND2),LWRK2,ISYMB,ISYMD,ISYM0)
      CALL CC_FCKMO(WORK(KFOCKBD),WORK(KLAMP0),WORK(KLAMHD),
     &              WORK(KEND2),LWRK2,ISYMB,ISYM0,ISYMD)
      CALL DAXPY(NORBT*NORBT,ONE,WORK(KFCKBUF),1,WORK(KFOCKBD),1)

C     ----------------------------------------------
C     add triples contribution: <mu_3|[B^D,T3^C]|HF>
C     ----------------------------------------------
      CALL CCSDT_XKSI3_2(WORK(KOMEGA3),WORK(KFOCKBD),WORK(KT3AM))

C     -------------------------------------------------
C     calculate the first-order triples amplitudes T^D:
C     -------------------------------------------------
      KDUM = KEND2
      CALL CCSDT_T31_NODDY(WORK(KT3AM),LISTD,IDLSTD,FREQD,.FALSE.,
     &                     .FALSE.,WORK(KINT1S),WORK(KINT2S),
     &                     .FALSE.,WORK(KDUM),WORK(KDUM),
     &                     .FALSE.,WORK(KDUM),WORK(KDUM),
     &                             WORK(KINT1SD),WORK(KINT2SD),
     &                     WORK(KLAMPD),WORK(KLAMHD),WORK(KPERTD),
     &                     WORK(KLAMP0),WORK(KLAMH0),WORK(KFOCK0),
     &                     WORK(KDUM),WORK(KFOCKD),
     &                     WORK(KEND2),LWRK2)

      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KT3AM),1)

C     --------------------------------------------------------------
C     calculate one-index transf. property integrals: B^C = [B,T1^C]
C     --------------------------------------------------------------
      CALL DCOPY(NORBT*NORBT,WORK(KFOCKB_AO),1,WORK(KFCKBUF),1)
      CALL DCOPY(NORBT*NORBT,WORK(KFOCKB_AO),1,WORK(KFOCKBC),1)
      CALL CC_FCKMO(WORK(KFCKBUF),WORK(KLAMPC),WORK(KLAMH0),
     &              WORK(KEND2),LWRK2,ISYMB,ISYMC,ISYMC)
      CALL CC_FCKMO(WORK(KFOCKBC),WORK(KLAMP0),WORK(KLAMHC),
     &              WORK(KEND2),LWRK2,ISYMB,ISYMC,ISYMC)
      CALL DAXPY(NORBT*NORBT,ONE,WORK(KFCKBUF),1,WORK(KFOCKBC),1)

C     ----------------------------------------------
C     add triples contribution: <mu_3|[B^C,T3^D]|HF>
C     ----------------------------------------------
      CALL CCSDT_XKSI3_2(WORK(KOMEGA3),WORK(KFOCKBC),WORK(KT3AM))

C     ---------------------------------------------------------
C     we have no contributions to the singles and doubles part:
C     ---------------------------------------------------------
      CALL DZERO(WORK(KOMEGA1),NT1AMX)
      CALL DZERO(WORK(KOMEGA2),NT1AMX*NT1AMX)

*---------------------------------------------------------------------*
*     Now we split:
*       for IOPTRES < 5 we compute the result vector
*       for IOPTRES = 5 we compute the contractions Tbar^E B{B} T^C T^D
*---------------------------------------------------------------------*
      IF (IOPTRES.GE.1 .AND. IOPTRES.LE.4) THEN

         CALL DCOPY(NT1AMX,OMEGA1,1,OMEGA1EFF,1)
         CALL DCOPY(NT2AMX,OMEGA2,1,OMEGA2EFF,1)
 
         CALL CC_RHPART_NODDY(OMEGA1EFF,OMEGA2EFF,WORK(KOMEGA3),FREQRES,
     &                        WORK(KFOCKD),WORK(KFOCK0),WORK(KFIELD),
     &                        WORK(KXIAJB),WORK(KINT1T),WORK(KINT2T),
     &                        WORK(KEND1),LWRK1)

      ELSE IF (IOPTRES.EQ.5) THEN

 
        SIGN = -1.0D0
        CALL CCDOTRSP_NODDY(WORK(KOMEGA1),WORK(KOMEGA2),
     &                      WORK(KOMEGA3),SIGN,
     &                      ITRAN,LISTDP,IDOTS,DOTPROD,MXVEC,
     &                      WORK(KLAMP0),WORK(KLAMH0),
     &                      WORK(KFOCK0),WORK(KFOCKD),
     &                      WORK(KXIAJB),WORK(KYIAJB),
     &                      WORK(KINT1T),WORK(KINT2T),
     &                      WORK(KINT1S),WORK(KINT2S),
     &                      'CCSDT_BAMAT_NODDY',LOCDBG,LOCDBG,.FALSE.,
     &                      WORK(KEND1),LWRK1)

      END IF

      CALL QEXIT('CCSDT_BAMAT_NODDY')
      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_BAMAT_NODDY                    *
*---------------------------------------------------------------------*
