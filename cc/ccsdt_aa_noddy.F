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
      SUBROUTINE CCSDT_AAMAT_NODDY(IOPTRES,FREQRES,LABELB,ISYMB,
     &                             LISTC,IDLSTC,OEFF_INIT,
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
      IMPLICIT NONE  
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccfield.h>
#include <ccorb.h>

      LOGICAL LOCDBG, PRINT_T3
      PARAMETER (LOCDBG=.FALSE., PRINT_T3=.FALSE.)
      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      CHARACTER LISTDP*3, LABELB*8, LISTC*3
      LOGICAL OEFF_INIT
      INTEGER LWORK, IOPTRES, ITRAN, MXVEC, NXTRAN, ISYMB, IDLSTC
      INTEGER IDOTS(MXVEC,NXTRAN)

#if defined (SYS_CRAY)
      REAL DOTPROD(MXVEC,NXTRAN), FREQC, FREQRES
      REAL WORK(LWORK), ONE, TWO, DUMMY, FF, SIGN, DDOT
      REAL OMEGA1(*),OMEGA2(*)
      REAL OMEGA1EFF(*),OMEGA2EFF(*)
#else
      DOUBLE PRECISION DOTPROD(MXVEC,NXTRAN), FREQC, FREQRES
      DOUBLE PRECISION WORK(LWORK), ONE, TWO, DUMMY, FF, SIGN, DDOT
      DOUBLE PRECISION OMEGA1(*),OMEGA2(*)
      DOUBLE PRECISION OMEGA1EFF(*),OMEGA2EFF(*)
#endif
      PARAMETER( ONE = 1.0D0, TWO = 2.0D0 )

      CHARACTER*10 MODEL
      INTEGER KFOCKB, KT3AM, KT2AM0, KT2AMC, KEND2, KEND1, KOMEGA1,
     &        KOMEGA2, KOMEGA3, KSCR1, KFOCKD, KFOCK0, KT1AMC,
     &        KFOCKB_AO, KFOCKBC, KLAMPC, KLAMHC, KFOCKC, KINT1SC, 
     &        KINT2SC, KINT1S, KINT2S, KXIAJB, KYIAJB, KINT1T, KINT2T,
     &        LWRK1, KDUM, LWRK2, KT1AMP0, KLAMP0, KLAMH0
      INTEGER IJ, NIJ, LUSIFC, INDEX, IDUMMY, ILSTSYM, ISYMC, LUFOCK, 
     &        IRREP, IERR, ILLL, IDEL, ISYDIS, IOPT, ISYMD, IVEC,
     &        IDLSTD, KFCKBUF, KFIELD, KFIELDAO

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J 

*---------------------------------------------------------------------*
*     Begin:
*---------------------------------------------------------------------*
      CALL QENTER('CCSDT_AAMAT_NODDY')

      IF (LOCDBG) WRITE(LUPRI,*) 'Entered CCSDT_AAMAT_NODDY...'

      IF(DIRECT)CALL QUIT('DIRECT NOT IMPLEMENTED IN CCSDT_AAMAT_NODDY')

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


      KFOCKB    = KEND1 
      KFOCKB_AO = KFOCKB    + NORBT*NORBT
      KFOCKBC   = KFOCKB_AO + NORBT*NORBT
      KFCKBUF   = KFOCKBC   + NORBT*NORBT
      KEND1     = KFCKBUF   + NORBT*NORBT

      IF (NONHF) THEN
        KFIELD   = KEND1
        KFIELDAO = KFIELD   + NORBT*NORBT
        KEND1    = KFIELDAO + NORBT*NORBT
      END IF

      KLAMPC  = KEND1
      KLAMHC  = KLAMPC + NLAMDT
      KFOCKC  = KLAMHC + NLAMDT
      KEND1   = KFOCKC + NORBT*NORBT

      KINT1S  = KEND1
      KINT2S  = KINT1S  + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2S  + NRHFT*NRHFT*NT1AMX 

      KINT1SC = KEND1
      KINT2SC = KINT1SC + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2SC + NRHFT*NRHFT*NT1AMX 

      KXIAJB  = KEND1
      KYIAJB  = KXIAJB  + NT1AMX*NT1AMX
      KEND1   = KYIAJB  + NT1AMX*NT1AMX

      KINT1T = KEND1
      KINT2T = KINT1T + NT1AMX*NVIRT*NVIRT
      KEND1  = KINT2T + NRHFT*NRHFT*NT1AMX 

      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_AAMAT_NODDY')
      ENDIF

      CALL DZERO(WORK(KOMEGA1),NT1AMX)

*---------------------------------------------------------------------*
*     Get zeroth-order Lambda matrices:
*---------------------------------------------------------------------*
      IOPT   = 1
      KDUM = KEND1
      Call CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AMP0),WORK(KDUM))

      Call LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AMP0),
     &            WORK(KEND1),LWRK1)

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
        CALL QUIT('CCSDT_AAMAT_NODDY: error reading operator '//LABELB)
      ELSE IF (IERR.LT.0) THEN
        CALL DZERO(WORK(KFOCKB),N2BST(ISYMB))
      END IF
 
      ! transform property integrals to Lambda-MO basis
      CALL CC_FCKMO(WORK(KFOCKB),WORK(KLAMP0),WORK(KLAMH0),
     &              WORK(KEND1),LWRK1,ISYMB,1,1)

*---------------------------------------------------------------------*
*     Compute some integrals:
*           XINT1S  =  (CK|BD)
*           XINT2S  =  (CK|LJ)
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
      KT2AM0 = KT3AM  + NT1AMX*NT1AMX*NT1AMX
      KT2AMC = KT2AM0 + NT1AMX*NT1AMX
      KT1AMC = KT2AMC + NT1AMX*NT1AMX
      KEND2  = KT1AMC + NT1AMX
      LWRK2  = LWORK  - KEND2
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_AAMAT_NODDY')
      ENDIF

      IOPT = 2
      CALL CC_RDRSP('R0',0,ISYMOP,IOPT,MODEL,DUMMY,WORK(KT3AM))
      CALL CC_T2SQ(WORK(KT3AM),WORK(KT2AM0),ISYMOP)

      ISYMC = ILSTSYM(LISTC,IDLSTC)
      IOPT = 3
      CALL CC_RDRSP(LISTC,IDLSTC,ISYMC,IOPT,MODEL,
     &              WORK(KT1AMC),WORK(KT3AM))
      CALL CCLR_DIASCL(WORK(KT3AM),TWO,ISYMC)
      CALL CC_T2SQ(WORK(KT3AM),WORK(KT2AMC),ISYMC)

      CALL CCLR_LAMTRA(WORK(KLAMP0),WORK(KLAMPC),
     &                 WORK(KLAMH0),WORK(KLAMHC),WORK(KT1AMC),ISYMC)

*---------------------------------------------------------------------*
*     compute zeroth-order triples response amplitudes and its 
*     contributions to the result vector
*---------------------------------------------------------------------*

C     --------------------------------------------------------------
C     compute zeroth-order triples; note: CCSDT_T3AM uses in non-hf
C     case WORK(KOMEGA3) as scratch vector for the iterative 
C     solution of the triples equations!
C     --------------------------------------------------------------
      CALL DZERO(WORK(KT3AM),NT1AMX*NT1AMX*NT1AMX) 
      CALL CCSDT_T03AM(WORK(KT3AM),WORK(KINT1S),WORK(KINT2S),
     &                 WORK(KT2AM0),WORK(KSCR1),WORK(KFOCKD),
     &                 WORK(KFIELD),WORK(KOMEGA3))

C     --------------------------------------------------------------
C     calculate one-index transf. property integrals: B^C = [B,T1^C]
C     --------------------------------------------------------------
      CALL DCOPY(NORBT*NORBT,WORK(KFOCKB_AO),1,WORK(KFCKBUF),1)
      CALL DCOPY(NORBT*NORBT,WORK(KFOCKB_AO),1,WORK(KFOCKBC),1)
      CALL CC_FCKMO(WORK(KFCKBUF),WORK(KLAMPC),WORK(KLAMH0),
     &              WORK(KEND2),LWRK2,ISYMB,ISYMC,ISYM0)
      CALL CC_FCKMO(WORK(KFOCKBC),WORK(KLAMP0),WORK(KLAMHC),
     &              WORK(KEND2),LWRK2,ISYMB,ISYM0,ISYMC)
      CALL DAXPY(NORBT*NORBT,ONE,WORK(KFCKBUF),1,WORK(KFOCKBC),1)

C     ---------------------------------
C     Initialize triples result vector:
C     ---------------------------------
      CALL DZERO(WORK(KOMEGA3),NT1AMX*NT1AMX*NT1AMX)

C     ----------------------------------------------
C     add triples contribution: <mu_3|[B^C,T3^0]|HF>
C     ----------------------------------------------
      CALL CCSDT_XKSI3_2(WORK(KOMEGA3),WORK(KFOCKBC),WORK(KT3AM))

*---------------------------------------------------------------------*
*     compute contribution from doubles amplitudes:
*---------------------------------------------------------------------*

C     ---------------------------------------------------
C     add triples contribution: <mu_3|[[B,T2^C],T2^0]|HF>
C     ---------------------------------------------------
      CALL CCSDT_XKSI3_1(WORK(KOMEGA3),WORK(KFOCKB),
     &                   WORK(KT2AM0),WORK(KT2AMC),ONE)
      CALL CCSDT_XKSI3_1(WORK(KOMEGA3),WORK(KFOCKB),
     &                   WORK(KT2AMC),WORK(KT2AM0),ONE)

*---------------------------------------------------------------------*
*     compute triples response amplitudes and its contributions
*     to the result vector
*---------------------------------------------------------------------*
      IF      (LISTC(1:3).EQ.'R1 ' .OR. LISTC(1:3).EQ.'RE ' .OR.
     &         LISTC(1:3).EQ.'RC '                              ) THEN

C       -----------------------------------------
C       calculate first-order triples amplitudes:
C       -----------------------------------------
        KDUM = KEND2
        CALL CCSDT_T31_NODDY(WORK(KT3AM),LISTC,IDLSTC,FREQC,.FALSE.,
     &                       .FALSE.,WORK(KINT1S),WORK(KINT2S),
     &                       .FALSE.,WORK(KDUM),WORK(KDUM),
     &                       .FALSE.,WORK(KDUM),WORK(KDUM),
     &                               WORK(KINT1SC),WORK(KINT2SC),
     &                       WORK(KLAMPC),WORK(KLAMHC),WORK(KFOCKC),
     &                       WORK(KLAMP0),WORK(KLAMH0),WORK(KFOCK0),
     &                       WORK(KDUM),WORK(KFOCKD),
     &                       WORK(KEND2),LWRK2)

      ELSE IF (LISTC(1:3).EQ.'R2 ' .OR. LISTC(1:3).EQ.'ER1') THEN

C       ------------------------------------------
C       calculate second-order triples amplitudes:
C       ------------------------------------------
        CALL CCSDT_T32_NODDY(WORK(KT3AM),LISTC,IDLSTC,FREQC,
     &                       WORK(KINT1S),WORK(KINT2S),
     &                       WORK(KLAMP0),WORK(KLAMH0),WORK(KFOCK0),
     &                       WORK(KFOCKD),WORK(KFIELDAO),WORK(KFIELD),
     &                       WORK(KSCR1),WORK(KEND2),LWRK2)

      ELSE

        CALL QUIT('Unknown list '//LISTC//' in CCSDT_AA_NODDY.')

      END IF

      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KT3AM),1)

      IF (LOCDBG) THEN
         WRITE(LUPRI,*) 'CCSDT_AAMAT_NODDY> T^C:'
         WRITE(LUPRI,*) 'CCSDT_AAMAT_NODDY> LISTC,IDLSTC:',LISTC,IDLSTC
         CALL PRINT_PT3_NODDY(WORK(KT3AM))
       END IF

C     -----------------------------------------------
C     add contribution to doubles: <mu_2|[B,T3^C]|HF>
C     -----------------------------------------------
      CALL DZERO(WORK(KOMEGA2),NT1AMX*NT1AMX)
      CALL CCSDT_XKSI2_2(WORK(KOMEGA2),WORK(KFOCKB),WORK(KT3AM))

      DO I = 1,NT1AMX
         DO J = 1,I
            IJ = NT1AMX*(I-1) + J
            NIJ = INDEX(I,J)
            OMEGA2(NIJ) = OMEGA2(NIJ) + WORK(KOMEGA2+IJ-1)
            IF (.NOT. OEFF_INIT) 
     &        OMEGA2EFF(NIJ) = OMEGA2EFF(NIJ) + WORK(KOMEGA2+IJ-1)
         END DO
      END DO

C     -------------------------------------------
C     contribution to triples: <mu_3|[B,T3^C]|HF>
C     -------------------------------------------
      CALL CCSDT_XKSI3_2(WORK(KOMEGA3),WORK(KFOCKB),WORK(KT3AM))

      if (print_t3) then
        write(lupri,*)'CCSDT_AA_NODDY> vector type:',listc
        write(lupri,*)'CCSDT_AA_NODDY> idlst:',idlstc
        write(lupri,*)'CCSDT_AA_NODDY> freq:',freqres
        call ccsdt_clean_t3(work(komega3),nt1amx,nvirt,nrhft)
        call print_pt3_noddy(work(komega3))
      end if

*---------------------------------------------------------------------*
*     Now we split:
*       for IOPTRES < 5 we compute the effective result vector
*       for IOPTRES = 5 we compute the contractions Tbar^D A{B} T^C
*---------------------------------------------------------------------*
      IF (IOPTRES.GE.1 .AND. IOPTRES.LE.4) THEN

         IF (OEFF_INIT) THEN
           CALL DCOPY(NT1AMX,OMEGA1,1,OMEGA1EFF,1)
           CALL DCOPY(NT2AMX,OMEGA2,1,OMEGA2EFF,1)
         END IF

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
     &                      'CCSDT_AAMAT_NODDY',LOCDBG,LOCDBG,.FALSE.,
     &                      WORK(KEND1),LWRK1)

      ELSE
        CALL QUIT('Illegal value for IOPTRES IN CCSDT_AAMAT_NODDY')
      END IF

      CALL QEXIT('CCSDT_AAMAT_NODDY')
      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_AAMAT_NODDY                    *
*---------------------------------------------------------------------*
