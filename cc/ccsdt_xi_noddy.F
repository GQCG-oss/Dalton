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
      SUBROUTINE CCSDT_XI_NODDY(IOPTRES,XI1,XI2,XI1EFF,XI2EFF,
     &                          LABELB,FOCKB,FREQ,ICAU,LCAUCHY,
     &                          XLAMDP,XLAMDH,
     &                          IDOTS,DOTPROD,LISTDP,ITRAN,
     &                          NXTRAN,MXVEC,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples contribution to XKSI vector
*
*             XKSI^eff_1,2 = XSKI_1,2(CCSD) + XKSI_1,2(T3)
*                            - A_1,2;3 (w_3 - w)^1 XKSI_3
*
*        
*     Written by Christof Haettig, April 2002, based on CCSDT_TRIPLE.
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccfield.h>
#include <ccorb.h>
#include <dummy.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.FALSE.)

      CHARACTER LISTDP*3, LABELB*8
      LOGICAL LCAUCHY
      INTEGER LWORK, IOPTRES, ITRAN, MXVEC, NXTRAN, ICAU, ISYMB
      INTEGER IDOTS(MXVEC,NXTRAN)

#if defined (SYS_CRAY)
      REAL DOTPROD(MXVEC,NXTRAN)
      REAL WORK(LWORK)
      REAL XI1(NT1AMX), XI2(NT2AMX), FREQ
      REAL XI1EFF(NT1AMX), XI2EFF(NT2AMX), DDOT, FF, SIGN, FREQB
      REAL FOCKB(*), CMO, XLAMDP(*), XLAMDH(*)
#else
      DOUBLE PRECISION DOTPROD(MXVEC,NXTRAN)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION XI1(NT1AMX), XI2(NT2AMX), FREQ, FF, SIGN
      DOUBLE PRECISION XI1EFF(NT1AMX), XI2EFF(NT2AMX), DDOT, FREQB
      DOUBLE PRECISION FOCKB(*), CMO, XLAMDP(*), XLAMDH(*)
#endif

      CHARACTER*10 MODEL
      INTEGER KEND1, LWRK1, KFOCKD, KIAJB, KYIAJB, KINT1T, KINT2T, 
     &        KXKSI3, KINT1S, KINT2S, IVEC, KEND2, LWRK2, NDOTS, KFIELD
      INTEGER IJ, NIJ, LUSIFC, INDEX, ILSTSYM, KFIELDAO, KFOCK0,
     &        KINT1SB, KINT2SB, KLAMPB, KLAMHB, KFCKB, IDLSTB

* external functions:
      INTEGER ILRCAMP

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J 

*---------------------------------------------------------------------*
*     Do some tests:
*---------------------------------------------------------------------*
      CALL QENTER('CCSDT_XI_NODDY')

      IF (DIRECT) CALL QUIT('DIRECT NOT IMPLEMENTED IN CCSDT_XI_NODDY')

      IF (IOPTRES.EQ.5) THEN
        ! check if anything to do at all
        NDOTS = 0
        IVEC  = 1
        DO WHILE (IDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
          NDOTS = NDOTS + 1
          IVEC  = IVEC  + 1
        END DO
        IF (NDOTS.EQ.0) THEN
          CALL QEXIT('CCSDT_XI_NODDY')
          RETURN
        END IF
      END IF

      IF (LOCDBG .AND. LCAUCHY) THEN
        WRITE(LUPRI,*) 'ICAU in noddy:',ICAU
      END IF

*---------------------------------------------------------------------*
*     Memory allocation:
*---------------------------------------------------------------------*
      KFOCKD = 1
      KFOCK0 = KFOCKD + NORBT
      KINT1T = KFOCK0 + NORBT*NORBT
      KINT2T = KINT1T + NT1AMX*NVIRT*NVIRT
      KIAJB  = KINT2T + NRHFT*NRHFT*NT1AMX
      KYIAJB = KIAJB  + NT1AMX*NT1AMX
      KINT1S = KYIAJB + NT1AMX*NT1AMX
      KINT2S = KINT1S + NT1AMX*NVIRT*NVIRT
      KEND1  = KINT2S + NRHFT*NRHFT*NT1AMX

      IF (NONHF) THEN
        KFIELD   = KEND1
        KFIELDAO = KFIELD   + NORBT*NORBT
        KEND1    = KFIELDAO + NORBT*NORBT
      END IF

      KXKSI3 = KEND1 
      KEND1  = KXKSI3 + NT1AMX*NT1AMX*NT1AMX

      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         write(lupri,*)'LWORK,KEND1 ', LWORK,KEND1
         CALL QUIT('Insufficient space in CCSDT_XI_NODDY')
      ENDIF

*---------------------------------------------------------------------*
*     Read precalculated integrals from file:
*---------------------------------------------------------------------*
      CALL CCSDT_READ_NODDY(.TRUE.,WORK(KFOCKD),WORK(KFOCK0),
     &                             WORK(KFIELD),WORK(KFIELDAO),
     &                      .TRUE.,WORK(KIAJB),WORK(KYIAJB),
     &                      .TRUE.,WORK(KINT1S),WORK(KINT2S),
     &                      .TRUE.,WORK(KINT1T),WORK(KINT2T),
     &                      .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &                      NORBT,NLAMDT,NRHFT,NVIRT,NT1AMX)

*---------------------------------------------------------------------*
*     Compute xksi_3 and add xksi_1,2(T3) contribution:
*---------------------------------------------------------------------*
      IF ((.NOT.LCAUCHY) .OR. ICAU.EQ.-1) THEN

        CALL CCSDT_XI3_NODDY(WORK(KXKSI3),FOCKB,.TRUE.,XI2,
     &                       WORK(KEND1),LWRK1)

      ELSE IF (LCAUCHY .AND. ICAU.GE.0) THEN

        KINT1SB = KEND1
        KINT2SB = KINT1SB + NT1AMX*NVIRT*NVIRT
        KLAMPB  = KINT2SB + NT1AMX*NRHFT*NRHFT
        KLAMHB  = KLAMPB  + NLAMDT
        KFCKB   = KLAMHB  + NLAMDT
        KEND2   = KFCKB   + NBAST*NBAST
        LWRK2   = LWORK   - KEND2
        IF (LWRK1 .LT. 0) THEN
           CALL QUIT('Insufficient space in CCSDT_XI_NODDY')
        ENDIF

        ISYMB  = 1
        IDLSTB = ILRCAMP(LABELB,ICAU,ISYMB)
        CALL CCSDT_T31_NODDY(WORK(KXKSI3),'RC ',IDLSTB,FREQB,.FALSE.,
     &                       .FALSE.,WORK(KINT1S),WORK(KINT2S),
     &                       .FALSE.,WORK(KINT1T),WORK(KINT2T),
     &                       .FALSE.,WORK(KIAJB), WORK(KYIAJB),
     &                       WORK(KINT1SB),WORK(KINT2SB),
     &                       WORK(KLAMPB),WORK(KLAMHB),WORK(KFCKB),
     &                       XLAMDP,XLAMDH,WORK(KFOCK0),DUMMY,
     &                       WORK(KFOCKD),WORK(KEND2),LWRK2)

      ELSE
         CALL QUIT('Error in CCSDT_XI_NODDY.')
      END IF

*---------------------------------------------------------------------*
*     Now we split:
*       for IOPTRES < 5 we compute the effective Xi vector
*       for IOPTRES = 5 we compute the contractions Tbar^C Xi
*---------------------------------------------------------------------*
      IF (IOPTRES.GE.1 .AND. IOPTRES.LE.4) THEN

         CALL DCOPY(NT1AMX,XI1,1,XI1EFF,1)
         CALL DCOPY(NT2AMX,XI2,1,XI2EFF,1)

         CALL CC_RHPART_NODDY(XI1EFF,XI2EFF,WORK(KXKSI3),FREQ,
     &                        WORK(KFOCKD),WORK(KFOCK0),WORK(KFIELD),
     &                        WORK(KIAJB),WORK(KINT1T),WORK(KINT2T),
     &                        WORK(KEND1),LWRK1)

      ELSE IF (IOPTRES.EQ.5) THEN
        IF (LCAUCHY) CALL QUIT('Inconsistency in CCSDT_XI_NODDY')

        SIGN = -1.0D0
        CALL CCDOTRSP_NODDY(DUMMY,DUMMY,WORK(KXKSI3),SIGN,
     &                      ITRAN,LISTDP,IDOTS,DOTPROD,MXVEC,
     &                      XLAMDP,XLAMDH,WORK(KFOCK0),WORK(KFOCKD),
     &                      WORK(KIAJB), WORK(KYIAJB),
     &                      WORK(KINT1T),WORK(KINT2T),
     &                      WORK(KINT1S),WORK(KINT2S),
     &                      'CCSDT_XI_NODDY',LOCDBG,.FALSE.,.FALSE.,
     &                      WORK(KEND1),LWRK1)

      ELSE
        CALL QUIT('Illegal value for IOPTRES IN CCSDT_XI_NODDY')
      END IF

      CALL QEXIT('CCSDT_XI_NODDY')
      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_XI_NODDY                       *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_XI3_NODDY(XI3,FOCKB,ADD_DOUBLE,XI2,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples part of XKSI vector 
*
*             if (add_double) compute also
*                XKSI_1,2(CC3) = XSKI_1,2(CCSD) + XKSI_1,2(T3)
*        
*     Written by Christof Haettig, April 2002, based on CCSDT_TRIPLE.
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccnoddy.h>
#include <dummy.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.FALSE.)

      LOGICAL ADD_DOUBLE
      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL XI2(NT2AMX), DDOT
      REAL XI3(NT1AMX,NT1AMX,NT1AMX)
      REAL FOCKB(*)
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION XI2(NT2AMX), DDOT
      DOUBLE PRECISION XI3(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION FOCKB(*)
#endif

      CHARACTER*10 MODEL
      INTEGER KEND1, LWRK1, KEND2, LWRK2, KT3AM,
     &        KSCR1, KT2AM, KXKSI2, IOPT, LUTEMP
      INTEGER IJ, NIJ, INDEX

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J 

      CALL QENTER('CCXI3NOD')

*---------------------------------------------------------------------*
*     Memory allocation:
*---------------------------------------------------------------------*
      KSCR1  = 1
      KT3AM  = KSCR1  + NT1AMX
      KT2AM  = KT3AM  + NT1AMX*NT1AMX*NT1AMX
      KEND1  = KT2AM  + NT1AMX*NT1AMX

      IF (ADD_DOUBLE) THEN
        KXKSI2 = KEND1
        KEND1  = KXKSI2 + NT1AMX*NT1AMX
      END IF

      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_XI3_NODDY')
      ENDIF

*---------------------------------------------------------------------*
*     End Loop over distributions of integrals.
*---------------------------------------------------------------------*
      IOPT = 2
      CALL CC_RDRSP('R0',0,ISYMOP,IOPT,MODEL,DUMMY,WORK(KT3AM))
      CALL CC_T2SQ(WORK(KT3AM),WORK(KT2AM),ISYMOP)

C     ---------------------------------------
C     read the zero-order triples amplitudes:
C     ---------------------------------------
      LUTEMP = -1
      CALL GPOPEN(LUTEMP,FILNODT30,'UNKNOWN',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      READ(LUTEMP) (WORK(KT3AM+I-1), I=1,NT1AMX*NT1AMX*NT1AMX)
      CALL GPCLOSE(LUTEMP,'KEEP')
      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KT3AM),1)


C     ------------------------------------------------------
C     Calculate triples contribution to doubles xksi vector:
C      xksi_mu_2 = "xksi_mu_2(CCSD)" + <mu_2|[V,T^0_3]|HF>
C     ------------------------------------------------------
      IF (ADD_DOUBLE) THEN
        CALL DZERO(WORK(KXKSI2),NT1AMX*NT1AMX)
        CALL CCSDT_XKSI2_2(WORK(KXKSI2),FOCKB,WORK(KT3AM))

        DO I = 1,NT1AMX
          DO J = 1,I
            IJ = NT1AMX*(I-1) + J
            NIJ = INDEX(I,J)
            XI2(NIJ) = XI2(NIJ) + WORK(KXKSI2+IJ-1)
          END DO
        END DO 
      END IF

C     ----------------------------------------------
C     Calculate the triples part of the xksi vector:
C     ----------------------------------------------
      CALL DZERO(XI3,NT1AMX*NT1AMX*NT1AMX)

      CALL CCSDT_XKSI3(XI3,FOCKB,WORK(KT3AM),WORK(KT2AM))
 
      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'CCSDT_XI3_NODDY> norm^2(T3AM):',
     *    DDOT(NT1AMX**3,WORK(KT3AM),1,WORK(KT3AM),1)
        WRITE(LUPRI,*) 'CCSDT_XI3_NODDY> norm^2(T2AM):',
     *    DDOT(NT1AMX**2,WORK(KT2AM),1,WORK(KT2AM),1)
        WRITE(LUPRI,*) 'CCSDT_XI3_NODDY> norm^2(XI3):',
     *    DDOT(NT1AMX**3,XI3,1,XI3,1)
      END IF

      CALL QEXIT('CCXI3NOD')
      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_XI3_NODDY                      *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_T31_NODDY(TB3AM,LISTB,IDLSTB,FREQB,XI_ONLY,
     &                           DO_INT0,XINT1S,XINT2S,
     &                           DO_INT1,XINT1T,XINT2T,
     &                           DO_INT2,XIAJB,YIAJB,
     &                           XINT1SB,XINT2SB,
     &                           XLAMPB,XLAMHB,FOCKB,
     &                           XLAMDP,XLAMDH,FOCK0,
     &                           CMO,FOCKD,
     &                           WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples part of first-order response amplitudes
*             and the property integrals as well as the integrals 
*              XINT1S = (ck|bd)-bar and XINT2S = (ck|lj)-bar
*
*             if (do_int0) compute also
*                 XINT1S = (ck|bd) and 
*                 XINT2S = (ck|lj)
*             else these integrals must be input
*             if (listb.eq.'RE ') do_int0 must be .false.
*
*             if (do_int1) compute also
*                 XINT1T = (kc|bd) and 
*                 XINT2T = (kc|lj)
*
*             if (do_int2) compute also
*                 XIAJB  = 2(kc|ld) - (kd|lc)
*                 YIAJB  =  (kc|ld)
*
*             if (xi_only) compute only right hand side (Xi) vector, 
*             not the solution (T1) vector; for RE return zero 
* 
*
*   Input:    LISTB, IDLSTB
*             XI_ONLY, DO_INT1, DO_INT2
*             XLAMDP,  XLAMDH
*             FOCKD
*             LWORK
*             
*   Output:   FREQB
*             XINT1T, XINT2T
*             XIAJB,  YIAJB
*             XLAMPB, XLAMHB
*             FOCKB
*
*   Scratch:  WORK
*
*   not needed (?): FOCK0, CMO
*
*   Implemented for R1 and RE vectors
*
*   Written by Christof Haettig, April 2002, based on CCSDT_XI3_NODDY.
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccr1rsp.h>
#include <ccrc1rsp.h>
#include <ccexci.h>
#include <ccfield.h>
#include <ccnoddy.h>
#include <dummy.h>

      LOGICAL LOCDBG, PRINT_T3
      PARAMETER (LOCDBG=.FALSE.,PRINT_T3=.FALSE.)

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      LOGICAL ADD_DOUBLE, DO_INT0, DO_INT1, DO_INT2, XI_ONLY, LCAUCHY

      CHARACTER LISTB*3
      INTEGER LWORK, IDLSTB
#if defined (SYS_CRAY)
      REAL WORK(LWORK), FOCKD(*)
      REAL FREQB, DDOT
      REAL TB3AM(NT1AMX,NT1AMX,NT1AMX)
      REAL XINT1T(NT1AMX*NVIRT*NVIRT)
      REAL XINT2T(NRHFT*NRHFT*NT1AMX)
      REAL XINT1S(NT1AMX*NVIRT*NVIRT)
      REAL XINT2S(NRHFT*NRHFT*NT1AMX)
      REAL XINT1SB(NT1AMX*NVIRT*NVIRT)
      REAL XINT2SB(NRHFT*NRHFT*NT1AMX)
      REAL XIAJB(NT1AMX*NT1AMX)
      REAL YIAJB(NT1AMX*NT1AMX)
      REAL FOCKB(*), CMO(*), XLAMDP(*), XLAMDH(*), FOCK0(*)
      REAL XLAMPB(*), XLAMHB(*), ZERO, ONE, TWO, FF
#else
      DOUBLE PRECISION WORK(LWORK), FOCKD(*)
      DOUBLE PRECISION FREQB, DDOT
      DOUBLE PRECISION TB3AM(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION XINT1T(NT1AMX*NVIRT*NVIRT)
      DOUBLE PRECISION XINT2T(NRHFT*NRHFT*NT1AMX)
      DOUBLE PRECISION XINT1S(NT1AMX*NVIRT*NVIRT)
      DOUBLE PRECISION XINT2S(NRHFT*NRHFT*NT1AMX)
      DOUBLE PRECISION XINT1SB(NT1AMX*NVIRT*NVIRT)
      DOUBLE PRECISION XINT2SB(NRHFT*NRHFT*NT1AMX)
      DOUBLE PRECISION XIAJB(NT1AMX*NT1AMX)
      DOUBLE PRECISION YIAJB(NT1AMX*NT1AMX)
      DOUBLE PRECISION FOCKB(*), CMO(*), XLAMDP(*), XLAMDH(*), FOCK0(*)
      DOUBLE PRECISION XLAMPB(*), XLAMHB(*), ZERO, ONE, TWO, FF
#endif
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)

      CHARACTER*10 MODEL, LABELB*8
      INTEGER KTB1AM, KEND0, KSCR1, KTB2AM, KEND1, NAI, NBJ, NCK,
     &        LWRK1, IOPT, KXINT, ISYMD, ILLL, IDEL, KEND2, LWRK2,
     &        KT02AM, ISYDIS, ISYMB, IERR, IRREP, NAIBJCK,
     &        ISYMM, ILSTSYM, LWRK0, KFIELD, KFLDB1, NCAU, MCAU,
     &        KFIELDAO, KT03AM, KT3SCR, LUTEMP, IDLSTBM

* external functions:
      INTEGER ILRCAMP


      CALL QENTER('CCT31NOD')
*---------------------------------------------------------------------*
*     Do some initializations:
*---------------------------------------------------------------------*
      IF (NSYM.NE.1) CALL QUIT('NO SYMMETRY FOR CCSDT_T31_NODDY')

      LCAUCHY = .FALSE.
      NCAU    = 0

      KEND0 = 1

      IF (NONHF) THEN
        KFIELD   = KEND0
        KFLDB1   = KFIELD   + NORBT*NORBT
        KFIELDAO = KFLDB1   + NORBT*NORBT
        KEND0    = KFIELDAO + NORBT*NORBT
      END IF

      LWRK0 = LWORK  - KEND0
      IF (LWRK0 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_T31_NODDY')
      ENDIF

c     ----------------------------------------------
c     Get the one electron integrals from the fields
c     ----------------------------------------------
      IF (NONHF .AND. NFIELD .GT. 0) THEN
        CALL DZERO(WORK(KFIELDAO),N2BST(ISYM0))
        DO I = 1, NFIELD
          FF = EFIELD(I)
          CALL CC_ONEP(WORK(KFIELDAO),WORK(KEND0),LWRK0,FF,1,LFIELD(I))
        ENDDO
        CALL DCOPY(NORBT*NORBT,WORK(KFIELDAO),1,WORK(KFIELD),1)

        ! calculate external field in zero-order lambda basis
        CALL CC_FCKMO(WORK(KFIELD),XLAMDP,XLAMDH,
     *                WORK(KEND0),LWRK0,1,1,1)
      ENDIF

*---------------------------------------------------------------------*
* R1: Compute xksi_3/eps_3 contribution to T^B_3 and some integrals:
*---------------------------------------------------------------------*
      IF (LISTB(1:3).EQ.'R1 '.OR.LISTB(1:3).EQ.'RC ') THEN

        ! get symmetry, frequency and integral label from common blocks
        ! defined in ccr1rsp.h
        IF (LISTB(1:3).EQ.'R1 ') THEN
          ISYMB  = ISYLRT(IDLSTB)
          FREQB  = FRQLRT(IDLSTB)
          LABELB = LRTLBL(IDLSTB)
          NCAU   = 0
          IF (LORXLRT(IDLSTB)) CALL QUIT(
     &      'NO ORBITAL RELAXED PERTURBATION IN CCSDT_T31_NODDY')
        ELSE IF (LISTB(1:3).EQ.'RC ') THEN
          ISYMB  = ISYLRC(IDLSTB)
          FREQB  = ZERO
          LABELB = LRCLBL(IDLSTB)
          NCAU   = ILRCAU(IDLSTB)
          LCAUCHY= .TRUE.
        ELSE
          CALL QUIT('Error in CCSDT_T31_NODDY.')
        END IF

        ! read property integrals from file:
        CALL CCPRPAO(LABELB,.TRUE.,FOCKB,IRREP,ISYMM,IERR,
     &               WORK(KEND0),LWRK0)
        IF ((IERR.GT.0) .OR. (IERR.EQ.0 .AND. IRREP.NE.ISYMB)) THEN
          CALL QUIT('CCSDT_T31_NODDY: error reading operator '//LABELB)
        ELSE IF (IERR.LT.0) THEN
          CALL DZERO(FOCKB,N2BST(ISYMB))
        END IF
 
        ! transform property integrals to Lambda-MO basis
        CALL CC_FCKMO(FOCKB,XLAMDP,XLAMDH,WORK(KEND0),LWRK0,ISYMB,1,1)

        IF (DO_INT0 .OR. DO_INT1 .OR. DO_INT2)
     &    CALL QUIT('CCSDT_XI3_NODDY does no longer compute integrals')

        CALL CCSDT_XI3_NODDY(TB3AM,FOCKB,.FALSE.,DUMMY,
     &                       WORK(KEND0),LWRK0)

        IF (LOCDBG) THEN
          WRITE(LUPRI,*) 'CCSDT_T31_NODDY: LISTB,IDLSTB =',LISTB,IDLSTB
          WRITE(LUPRI,*) 'CCSDT_T31_NODDY: LABELB =',LABELB
          WRITE(LUPRI,*) 'CCSDT_T31_NODDY:  FREQB =',FREQB
          WRITE(LUPRI,*) 'CCSDT_T31_NODDY:LCAUCHY =',LCAUCHY
          WRITE(LUPRI,*) 'CCSDT_T31_NODDY:   NCAU =',NCAU
          WRITE(LUPRI,*) 'NORM^2(TB3AM/XKI3)',
     &     DDOT(NT1AMX*NT1AMX*NT1AMX,TB3AM,1,TB3AM,1)
          WRITE(LUPRI,*) 'NORM^2(FOCKB):',
     &     DDOT(NORBT*NORBT,FOCKB,1,FOCKB,1)
C         WRITE(LUPRI,*) 'TB3AM after CCSDT_XI3_NODDY:'
C         CALL OUTPUT(TB3AM,1,NT1AMX*NT1AMX,1,NT1AMX,
C    &                        NT1AMX*NT1AMX,NT1AMX,1,LUPRI)
        END IF

*---------------------------------------------------------------------*
* RE: for excited states the rhs vector is zero:
*---------------------------------------------------------------------*
      ELSE IF (LISTB(1:3).EQ.'RE ') THEN

        ISYMB  = ILSTSYM(LISTB,IDLSTB)
        FREQB  = EIGVAL(IDLSTB)
        NCAU   = 0
        CALL DZERO(TB3AM,NT1AMX*NT1AMX*NT1AMX)

        IF (DO_INT0) THEN
          CALL QUIT('Need XINT1S/XINT2S as input in CCSDT_T31_NODDY.')
        END IF

      ELSE
        CALL QUIT('Unknown/Illegal LIST in CCSDT_T31_NODDY.')
      END IF

      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'NORM^2(TB3AM) at 1:',
     &   DDOT(NT1AMX*NT1AMX*NT1AMX,TB3AM,1,TB3AM,1)
      END IF

*---------------------------------------------------------------------*
*     For Cauchy vectors start now loop over A x C type contributions:
*---------------------------------------------------------------------*

      DO MCAU = 0, NCAU

        IDLSTBM = IDLSTB
        IF (LCAUCHY) THEN
          IDLSTBM = ILRCAMP(LABELB,MCAU,ISYMB)
          IF (MCAU.EQ.NCAU .AND. IDLSTBM.NE.IDLSTB) 
     &     CALL QUIT('Inconsistency in Cauchy loop in CCSDT_T31_NODDY')
        END IF

*---------------------------------------------------------------------*
*     Some allocations:
*---------------------------------------------------------------------*
      KSCR1  = KEND0
      KTB1AM = KSCR1  + NT1AMX
      KTB2AM = KTB1AM + NT1AMX
      KT02AM = KTB2AM + NT1AMX*NT1AMX
      KEND1  = KT02AM + NT1AMX*NT1AMX
      LWRK1  = LWORK  - KEND1

      LWRK1  = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_T31_NODDY (1) ')
      ENDIF

*---------------------------------------------------------------------*
*     Read zero-order amplitudes:
*---------------------------------------------------------------------*
      IF (LWRK1.LT.NT2AMX) 
     &   CALL QUIT('Insufficient space in CCSDT_T31_NODDY (2) ')

      IOPT = 2
      CALL CC_RDRSP('R0',0,ISYMOP,IOPT,MODEL,DUMMY,WORK(KEND1))
      CALL CC_T2SQ(WORK(KEND1),WORK(KT02AM),ISYMOP)

*---------------------------------------------------------------------*
*     Read response amplitudes and compute response lambda matrices:
*---------------------------------------------------------------------*
      IOPT = 3 ! read singles and doubles amplitudes
      CALL CC_RDRSP(LISTB,IDLSTBM,ISYMB,IOPT,MODEL,
     &              WORK(KTB1AM),WORK(KEND1))
      Call CCLR_DIASCL(WORK(KEND1),TWO,ISYMB)
      CALL CC_T2SQ(WORK(KEND1),WORK(KTB2AM),ISYMB)

      CALL CCLR_LAMTRA(XLAMDP,XLAMPB,XLAMDH,XLAMHB,WORK(KTB1AM),ISYMB)

*---------------------------------------------------------------------*
*     Compute contribution from jacobian transformation:
*       <mu_3|[H^B,T^0_2]|HF> + <mu_3|[H,T^B_2]|HF>   
*---------------------------------------------------------------------*
      IF (.NOT. XI_ONLY) THEN

        IF ((NONHF) .AND. (NFIELD .GT. 0)) THEN
c         ----------------------------------
c         get zero-order triples amplitudes:
c         ----------------------------------
          KT03AM = KEND1
          KEND2  = KT03AM + NT1AMX*NT1AMX*NT1AMX
          LWRK2  = LWORK  - KEND2
          IF (LWRK2 .LT. 0) THEN
             CALL QUIT('Insufficient space in CCSDT_T31_NODDY (nonhf)')
          ENDIF

          LUTEMP = -1
          CALL GPOPEN(LUTEMP,FILNODT30,'UNKNOWN',' ','UNFORMATTED',
     &                IDUMMY,.FALSE.)
          READ(LUTEMP) (WORK(KT03AM+I-1), I=1,NT1AMX*NT1AMX*NT1AMX)
          CALL GPCLOSE(LUTEMP,'KEEP')
          CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KT03AM),1)

        ELSE
          KEND2 = KEND1
          LWRK2 = LWRK1
        END IF

        CALL CCSDT_A3AM(TB3AM,WORK(KTB1AM),WORK(KTB2AM),FREQB,ISYMB,
     &                  WORK(KT02AM),WORK(KT03AM),
     &                  XINT1S,XINT2S,XINT1SB,XINT2SB,
     &                  FOCKD,XLAMDP,XLAMDH,
     &                  WORK(KFIELDAO),WORK(KFIELD),WORK(KSCR1),
     &                  WORK(KEND2),LWRK2)

      END IF

*---------------------------------------------------------------------*
*     Solve triples equations:
*      T_3 = { Xksi_3 + <mu_3|[H^B,T^0_2]|HF>
*                     + <mu_3|[H,T^B_2]|HF>   } / eps_3
*---------------------------------------------------------------------*
      KSCR1 = KEND0
      KEND1 = KSCR1 + NT1AMX
      IF (NONHF) THEN
        KT3SCR = KEND1
        KEND1  = KT3SCR + NT1AMX*NT1AMX*NT1AMX
      END IF
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_T31_NODDY (solv)')
      ENDIF

      CALL CCSDT_3AM(TB3AM,FREQB,WORK(KSCR1),FOCKD,
     &               NONHF,WORK(KFIELD),.FALSE.,WORK(KT3SCR))

*---------------------------------------------------------------------*
*     End loop over A x C contributions to Cauchy vectors
*---------------------------------------------------------------------*
      END DO ! MCAU

*---------------------------------------------------------------------*
*     turn the sign, since
*      T_3 = - { Xksi_3 + <mu_3|[H^B,T^0_2]|HF>
*                       + <mu_3|[H,T^B_2]|HF>   } / eps_3
*---------------------------------------------------------------------*
      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,TB3AM,1)

      CALL CCSDT_CLEAN_T3(TB3AM,NT1AMX,NVIRT,NRHFT)

*---------------------------------------------------------------------*
*     print debug output and return
*---------------------------------------------------------------------*
      IF (LOCDBG.OR.PRINT_T3) THEN
        WRITE(LUPRI,*)'CCSDT_T31_NODDY> first-order T3 vector:'
        WRITE(LUPRI,*)'CCSDT_T31_NODDY> xi_only:',xi_only
        WRITE(LUPRI,*)'CCSDT_T31_NODDY> list,idlst:',listb,idlstb
        WRITE(LUPRI,*)'CCSDT_T31_NODDY> freq,label:',freqb,labelb
        WRITE(LUPRI,*) 'NORM^2(TB3AM)',
     &   DDOT(NT1AMX*NT1AMX*NT1AMX,TB3AM,1,TB3AM,1)
      END IF
      IF (PRINT_T3) THEN
        CALL PRINT_PT3_NODDY(TB3AM)
      END IF

      CALL QEXIT('CCT31NOD')
      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_T31_NODDY                      *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_INTS1_NODDY(DO_INTS1,XINT1SB,XINT2SB,
     &                             DO_INTT1,XINT1TB,XINT2TB,
     &                             XLAMDP,XLAMDH,XLAMPB,XLAMHB,
     &                             WORK,LWORK)
*---------------------------------------------------------------------*
*  Purpose:
*     Loop over distributions of integrals and compute first-order
*     response versions XINT1SB=(ck|db)-bar and XINT2SB=(ck|lj)-bar
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccorb.h>

      LOGICAL DO_INTS1, DO_INTT1
      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL XINT1SB(NT1AMX*NVIRT*NVIRT)
      REAL XINT2SB(NRHFT*NRHFT*NT1AMX)
      REAL PRECISION XINT1TB(NT1AMX*NVIRT*NVIRT)
      REAL PRECISION XINT2TB(NRHFT*NRHFT*NT1AMX)  
      REAL XLAMDP(*), XLAMDH(*), XLAMPB(*), XLAMHB(*)
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION XINT1SB(NT1AMX*NVIRT*NVIRT)
      DOUBLE PRECISION XINT2SB(NRHFT*NRHFT*NT1AMX)
      DOUBLE PRECISION XINT1TB(NT1AMX*NVIRT*NVIRT)
      DOUBLE PRECISION XINT2TB(NRHFT*NRHFT*NT1AMX)  
      DOUBLE PRECISION XLAMDP(*), XLAMDH(*), XLAMPB(*), XLAMHB(*)
#endif  

      INTEGER ISYMD, ILLL, IDEL, ISYDIS, KEND2, LWRK2, KXINT, IRECNR

      IF (DO_INTS1) THEN
        CALL DZERO(XINT1SB,NT1AMX*NVIRT*NVIRT)
        CALL DZERO(XINT2SB,NT1AMX*NRHFT*NRHFT)
      END IF
      IF (DO_INTT1) THEN
        CALL DZERO(XINT1TB,NT1AMX*NVIRT*NVIRT)
        CALL DZERO(XINT2TB,NT1AMX*NRHFT*NRHFT)
      END IF

      DO ISYMD = 1, NSYM
         DO ILLL = 1,NBAS(ISYMD)
            IDEL   = IBAS(ISYMD) + ILLL
            ISYDIS = MULD2H(ISYMD,ISYMOP)
 
C           ----------------------------
C           Work space allocation no. 2.
C           ----------------------------
            KXINT  = 1
            KEND2  = KXINT + NDISAO(ISYDIS)
            LWRK2  = LWORK - KEND2
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
               CALL QUIT('Insufficient space in CCSDT_INTS1_NODDY')
            ENDIF
 
C           ---------------------------
C           Read in batch of integrals.
C           ---------------------------
            CALL CCRDAO(WORK(KXINT),IDEL,1,WORK(KEND2),LWRK2,
     *                  IRECNR,DIRECT)
 
C           ---------------------------------------
C           Calculate transformed integrals needed: 
C           ---------------------------------------
            IF (DO_INTS1) THEN

              ! XINT1SB = XINT1SB + (C-bar K|B D)
              ! XINT2SB = XINT2SB + (C-bar K|L J)
              CALL CCSDT_TRAN3_R(XINT1SB,XINT2SB,XLAMDP,XLAMDH,
     &                           XLAMPB,XLAMDH,XLAMDP,XLAMDH,
     &                           WORK(KXINT),IDEL)
             
              ! XINT1SB = XINT1SB + (C K-bar|B D)
              ! XINT2SB = XINT2SB + (C K-bar|L J)
              CALL CCSDT_TRAN3_R(XINT1SB,XINT2SB,XLAMDP,XLAMDH,
     &                           XLAMDP,XLAMHB,XLAMDP,XLAMDH,
     &                           WORK(KXINT),IDEL)
             
              ! XINT1SB = XINT1SB + (C K|B-bar D)
              ! XINT2SB = XINT2SB + (C K|L J-bar)
              CALL CCSDT_TRAN3_R(XINT1SB,XINT2SB,XLAMDP,XLAMDH,
     &                           XLAMDP,XLAMDH,XLAMPB,XLAMHB,
     &                           WORK(KXINT),IDEL)
 
            END IF

            IF (DO_INTT1) THEN
              CALL CCSDT_TRAN1_R(XINT1TB,XINT2TB,
     &                           XLAMDP,XLAMDH,XLAMPB,XLAMHB,
     &                           WORK(KXINT),IDEL)
            END IF

         END DO   
      END DO  

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_INTS1_NODDY                    *
*---------------------------------------------------------------------* 
*=====================================================================*
      SUBROUTINE CCSDT_XKSI3(XKSI3,FOCKB,T3AM,T2AM)
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
 
#if defined (SYS_CRAY)
      REAL ONE
#else
      DOUBLE PRECISION ONE
#endif
      PARAMETER (ONE = 1.0D0)
 
#if defined (SYS_CRAY)
      REAL XKSI3(NT1AMX,NT1AMX,NT1AMX)
      REAL T2AM(NT1AMX,NT1AMX), FOCKB(NORBT,NORBT)
      REAL T3AM(NT1AMX,NT1AMX,NT1AMX)
#else
      DOUBLE PRECISION XKSI3(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION T2AM(NT1AMX,NT1AMX), FOCKB(NORBT,NORBT)
      DOUBLE PRECISION T3AM(NT1AMX,NT1AMX,NT1AMX)
#endif

      CALL CCSDT_XKSI3_1(XKSI3,FOCKB,T2AM,T2AM,ONE)
    
      CALL CCSDT_XKSI3_2(XKSI3,FOCKB,T3AM)

      CALL CCSDT_CLEAN_T3(T3AM,NT1AMX,NVIRT,NRHFT)

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_XKSI3                          *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_XKSI3_1(XKSI3,FOCKA,T2AMB,T2AMC,FAC)
*---------------------------------------------------------------------*
* Purpose: add fac/2 *<mu_3|[[F^A,T2^B],T2^C]|HF> 
*          to XKSI3 result vector
*          Note: assumes T2AMB and T2AMC to be identical, it this
*                is not the case, the routine must be called twice
*                with T2AMB and T2AMC interchanged in second call
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
 
#if defined (SYS_CRAY)
      REAL HALF
#else
      DOUBLE PRECISION HALF
#endif
      PARAMETER (HALF = 0.5D0)
 
#if defined (SYS_CRAY)
      REAL AIBJCK, FAC
      REAL XKSI3(NT1AMX,NT1AMX,NT1AMX)
      REAL FOCKA(NORBT,NORBT)
      REAL T2AMB(NT1AMX,NT1AMX)
      REAL T2AMC(NT1AMX,NT1AMX)
#else
      DOUBLE PRECISION AIBJCK, FAC
      DOUBLE PRECISION XKSI3(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION FOCKA(NORBT,NORBT)
      DOUBLE PRECISION T2AMB(NT1AMX,NT1AMX)
      DOUBLE PRECISION T2AMC(NT1AMX,NT1AMX)
#endif
    
      INTEGER NCK, NBJ, NAI, NBI, NCI, NAJ, NAK, NDJ, NBL
 
      DO NCK = 1,NT1AMX
         DO J = 1,NRHFT
            DO B = 1,NVIRT
               NBJ = NVIRT*(J-1) + B
               DO NAI = 1,NT1AMX
                  AIBJCK = 0.0D0
 
                  DO D = 1, NVIRT
                     NDJ = NVIRT*(J-1) + D
                     DO L = 1, NRHFT
                        NBL = NVIRT*(L-1) + B
                        AIBJCK = AIBJCK
     *                         - T2AMB(NAI,NBL)
     *                          *T2AMC(NCK,NDJ)
     *                          *FOCKA(L,NRHFT+D)
                     END DO
                  END DO

                  AIBJCK = AIBJCK * FAC

                  XKSI3(NAI,NBJ,NCK) = XKSI3(NAI,NBJ,NCK) + AIBJCK
                  XKSI3(NAI,NCK,NBJ) = XKSI3(NAI,NCK,NBJ) + AIBJCK
                  XKSI3(NBJ,NAI,NCK) = XKSI3(NBJ,NAI,NCK) + AIBJCK
                  XKSI3(NCK,NAI,NBJ) = XKSI3(NCK,NAI,NBJ) + AIBJCK
                  XKSI3(NBJ,NCK,NAI) = XKSI3(NBJ,NCK,NAI) + AIBJCK
                  XKSI3(NCK,NBJ,NAI) = XKSI3(NCK,NBJ,NAI) + AIBJCK
 
               END DO
            END DO
         END DO
      END DO

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_XKSI3_1                        *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_XKSI3_2(XKSI3,FOCKB,T3AM)
*---------------------------------------------------------------------*
* Purpose: add <mu_3|[F^B,T3]|HF> to XKSI3 vector
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
 
#if defined (SYS_CRAY)
      REAL HALF
#else
      DOUBLE PRECISION HALF
#endif
      PARAMETER (HALF = 0.5D0)
 
#if defined (SYS_CRAY)
      REAL AIBJCK
      REAL XKSI3(NT1AMX,NT1AMX,NT1AMX)
      REAL T3AM(NT1AMX,NT1AMX,NT1AMX), FOCKB(NORBT,NORBT)
#else
      DOUBLE PRECISION AIBJCK
      DOUBLE PRECISION XKSI3(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION T3AM(NT1AMX,NT1AMX,NT1AMX), FOCKB(NORBT,NORBT)
#endif
    
      INTEGER NCK, NBJ, NAI, NBI, NCI, NAJ, NAK, NDJ, NBL
 
      DO NCK = 1,NT1AMX
         DO J = 1,NRHFT
            DO B = 1,NVIRT
               NBJ = NVIRT*(J-1) + B
               DO NAI = 1,NT1AMX
                  AIBJCK = 0.0D0
 
                  ! remember: T3 has turned sign
                  DO D = 1, NVIRT
                     NDJ = NVIRT*(J-1) + D
                     AIBJCK = AIBJCK
     *                      - HALF*T3AM(NAI,NDJ,NCK)
     *                            *FOCKB(NRHFT+B,NRHFT+D)
                  END DO
 
                  ! remember: T3 has turned sign
                  DO L = 1, NRHFT
                     NBL = NVIRT*(L-1) + B
                     AIBJCK = AIBJCK
     *                      + HALF*T3AM(NAI,NBL,NCK)
     *                         *FOCKB(L,J)
                  END DO
 
                  XKSI3(NAI,NBJ,NCK) = XKSI3(NAI,NBJ,NCK) + AIBJCK
                  XKSI3(NAI,NCK,NBJ) = XKSI3(NAI,NCK,NBJ) + AIBJCK
                  XKSI3(NBJ,NAI,NCK) = XKSI3(NBJ,NAI,NCK) + AIBJCK
                  XKSI3(NCK,NAI,NBJ) = XKSI3(NCK,NAI,NBJ) + AIBJCK
                  XKSI3(NBJ,NCK,NAI) = XKSI3(NBJ,NCK,NAI) + AIBJCK
                  XKSI3(NCK,NBJ,NAI) = XKSI3(NCK,NBJ,NAI) + AIBJCK
 
               END DO
            END DO
         END DO
      END DO

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_XKSI3_2                        *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_XKSI2(XKSI2,FOCKB,T3AM)
*---------------------------------------------------------------------*
*     Purpose: compute triples contribution to double excitation 
*              part of the Xski vector
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>

#if defined (SYS_CRAY)
      REAL XKSI2(NT1AMX,NT1AMX)
      REAL FOCKB(NORBT,NORBT), CONTRIB
      REAL T3AM(NT1AMX,NT1AMX,NT1AMX)
#else
      DOUBLE PRECISION XKSI2(NT1AMX,NT1AMX)
      DOUBLE PRECISION FOCKB(NORBT,NORBT), CONTRIB
      DOUBLE PRECISION T3AM(NT1AMX,NT1AMX,NT1AMX)
#endif
    
      INTEGER NAI, NBJ, NBK, NCK, NCJ, INDEX

      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J 

      DO NAI = 1, NT1AMX
        DO J = 1,NRHFT
          DO B = 1,NVIRT
            NBJ = NVIRT*(J-1) + B
            DO K = 1,NRHFT
              NBK = NVIRT*(K-1) + B
              DO C = 1,NVIRT
                NCK = NVIRT*(K-1) + C
                NCJ = NVIRT*(J-1) + C
                CONTRIB = FOCKB(K,NRHFT+C) *
     &                        (T3AM(NAI,NBJ,NCK) - T3AM(NAI,NBK,NCJ))
                XKSI2(NBJ,NAI) = XKSI2(NBJ,NAI) - CONTRIB
                XKSI2(NAI,NBJ) = XKSI2(NAI,NBJ) - CONTRIB
 
              END DO
            END DO
          END DO
        END DO
      END DO
 
      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_XKSI2                          *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_XKSI2_2(OMEGA2,FOCK,T3AM)
*---------------------------------------------------------------------*
*     Purpose: compute triples contribution to double excitation 
*              part of the Xski vector
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>

#if defined (SYS_CRAY)
      REAL OMEGA2(NT1AMX,NT1AMX)
      REAL FOCK(NORBT,NORBT), XAIBJ
      REAL T3AM(NT1AMX,NT1AMX,NT1AMX)
#else
      DOUBLE PRECISION OMEGA2(NT1AMX,NT1AMX)
      DOUBLE PRECISION FOCK(NORBT,NORBT), XAIBJ
      DOUBLE PRECISION T3AM(NT1AMX,NT1AMX,NT1AMX)
#endif
      
      INTEGER NAI, NBJ, NBK, NAK, NCK, NCJ, NCI
    
      DO 100 I = 1,NRHFT
         DO 110 A = 1,NVIRT
            NAI = NVIRT*(I-1) + A
C
            DO 120 J = 1,NRHFT
               DO 130 B = 1,NVIRT
                  NBJ = NVIRT*(J-1) + B
C
                  DO 140 K = 1,NRHFT
                     NBK = NVIRT*(K-1) + B
                     NAK = NVIRT*(K-1) + A
                     DO 150 C = 1,NVIRT
C
                        NCK = NVIRT*(K-1) + C
                        NCJ = NVIRT*(J-1) + C
                        NCI = NVIRT*(I-1) + C
C
                        OMEGA2(NBJ,NAI) = OMEGA2(NBJ,NAI) -
     *      (T3AM(NAI,NBJ,NCK) - T3AM(NAI,NBK,NCJ))*FOCK(K,NRHFT+C)
C
C
  150                CONTINUE
  140             CONTINUE
C
  130          CONTINUE
  120       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      DO 200 NAI = 1,NT1AMX
         DO 210 NBJ = 1,NAI
C
            XAIBJ = OMEGA2(NAI,NBJ) + OMEGA2(NBJ,NAI)
            OMEGA2(NAI,NBJ) = XAIBJ
            OMEGA2(NBJ,NAI) = XAIBJ
C
  210    CONTINUE
  200 CONTINUE
 
      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_XKSI2                          *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_FCK_R(FOCK,XLIAJB,T1AM)
*---------------------------------------------------------------------*
C
C     Compute a response Fock matrix:
C        FOCK(kc) = FOCK(kc) + sum_ia L_kcia t_ia
C
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
      DIMENSION XLIAJB(NT1AMX,NT1AMX)
      DIMENSION T1AM(NT1AMX)
      DIMENSION FOCK(NORBT,NORBT)

      DO NK = 1, NRHFT
        DO NC = 1, NVIRT
          NCK = (NK-1)*NVIRT + NC
          DO NAI = 1, NT1AMX
            FOCK(NK,NRHFT+NC) = FOCK(NK,NRHFT+NC) + 
     &         XLIAJB(NCK,NAI) * T1AM(NAI)
          END DO
        END DO
      END DO

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_FCK_R                          *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_XI_CONT_NODDY(LISTL,DOTPROD,
     &                               IDXL_XIDEN,NXIDEN,
     &                               IXDOTS,IXETRAN,MXVEC,NXETRAN,
     &                               EMULATE_CCSDPT_DENS2,
     &                               FNDPTIA,FNDPTIA2,FNDPTAB,FNDPTIJ,
     &                               WORK,LWORK)
*---------------------------------------------------------------------*
*     Purpose: compute density corresponding to a contraction of      *
*              the triples parts of a Xi vector and (first-order)     *
*              lagrangian multipliers                                 *
*                                                                     *
*     Written by Christof Haettig, Nov 2002, Friedrichstal            *
*---------------------------------------------------------------------*
      IMPLICIT NONE  

#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccfield.h>
#include <ccorb.h>
#include <cclists.h>
#include <ccroper.h>

      CHARACTER*12 FNDEFF
      PARAMETER( FNDEFF = 'TMP-XIETADEN' )

      CHARACTER*(*) FNDPTIA, FNDPTAB, FNDPTIJ, FNDPTIA2

C     CHARACTER*5 FNDPTIA, FNDPTAB, FNDPTIJ
C     CHARACTER*6 FNDPTIA2 
C     PARAMETER(FNDPTIA='DPTIA', FNDPTIA2 = 'DPTIA2',
C    *          FNDPTAB='DPTAB' ,FNDPTIJ  = 'DPTIJ'  )

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYM0, ISYDEN
      PARAMETER( ISYM0=1, ISYDEN=1 )

      CHARACTER*3 LISTL
      LOGICAL EMULATE_CCSDPT_DENS2
      INTEGER LWORK

      INTEGER NXETRAN, MXVEC, NXIDEN, IRECORD
      INTEGER IXETRAN(MXDIM_XEVEC,NXETRAN)
      INTEGER IXDOTS(MXVEC,NXETRAN)
      INTEGER IDXL_XIDEN(NXIDEN)

#if defined (SYS_CRAY)
      REAL WORK(LWORK), DUMMY, FREQL, TRIPCON, FF
      REAL DOTPROD(MXVEC,NXETRAN)
#else
      DOUBLE PRECISION WORK(LWORK), DUMMY, FREQL, TRIPCON, FF
      DOUBLE PRECISION DOTPROD(MXVEC,NXETRAN)
#endif

      CHARACTER MODEL*10, LABELH*8
      LOGICAL LTWOEL, LRELAX, SKIPXI
      INTEGER KINT1T, KINT2T, KXIAJB, KYIAJB, KINT1S, KINT2S, KEND1,
     &        ISYMD, ILLL, IDEL, ISYDIS, KXINT, KEND2, LWRK2, LWRK1,
     &        IRECNR, KT2AM, KT3AM, IOPT, LUSIFC, IDLSTL, KFOCK0,
     &        LUDEFF, M2BST, KFOCKD, ISYM, IDUMMY, LUFOCK, KSCR1, IF,
     &        KLC3AM, KDIJ, KDAB, KDIA, KMINT, IDENS, IADRIJ, IADRAB,
     &        IADRIA, ITRAN, IOPER, IRELAX, ISYHOP, ISYML, IDX, IVEC,
     &        KLAMP0, KLAMH0, KT1AMP, KLC2AM, KLC1AM, KDIA2, KFIELD,
     &        KT3AM2, KEND1A, LWRK1A, LUT3AM

      INTEGER ILSTSYM
#if defined (SYS_CRAY)
      REAL CC_XIETA_DENCON
#else
      DOUBLE PRECISION CC_XIETA_DENCON
#endif

      
*---------------------------------------------------------------------*
*     some initializations:
*---------------------------------------------------------------------*
      CALL QENTER('CCSDT_XI_CONT_NODDY')

      IF (.NOT. EMULATE_CCSDPT_DENS2) THEN
        LUDEFF = -1
        CALL WOPEN2(LUDEFF,FNDEFF,64,0)
      END IF

      M2BST = 0
      DO ISYM = 1, NSYM
        M2BST = MAX(M2BST,N2BST(ISYM))
      END DO 

*---------------------------------------------------------------------*
*     Memory allocation:
*---------------------------------------------------------------------*
      KLAMP0 = 1
      KLAMH0 = KLAMP0 + NLAMDT
      KT1AMP = KLAMH0 + NLAMDT
      KFOCK0 = KT1AMP + NT1AMX
      KFOCKD = KFOCK0 + NORBT*NORBT
      KINT1T = KFOCKD + NORBT
      KINT2T = KINT1T + NT1AMX*NVIRT*NVIRT
      KXIAJB = KINT2T + NRHFT*NRHFT*NT1AMX
      KYIAJB = KXIAJB + NT1AMX*NT1AMX
      KINT1S = KYIAJB + NT1AMX*NT1AMX
      KINT2S = KINT1S + NT1AMX*NVIRT*NVIRT
      KEND1  = KINT2S + NRHFT*NRHFT*NT1AMX

      IF ((NONHF) .AND. (NFIELD .GT. 0)) THEN
         KFIELD = KEND1
         KEND1  = KFIELD + N2BST(ISYM0)
      ENDIF

      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_XI_CONT_NODDY')
      ENDIF

      IF (DIRECT) CALL QUIT(
     &         'DIRECT NOT IMPLEMENTED IN CCSDT_XI_CONT_NODDY')

*---------------------------------------------------------------------*
*     initialize 0.th-order Lambda:
*---------------------------------------------------------------------*
      IOPT = 1
      CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AMP),DUMMY)
 
      CALL LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AMP),
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
*     Get the one electron integrals for external fields
*---------------------------------------------------------------------*
      IF ((NONHF) .AND. NFIELD .GT. 0) THEN
         CALL DZERO(WORK(KFIELD),N2BST(ISYM0))
         DO IF = 1, NFIELD
           FF = EFIELD(IF)
           CALL CC_ONEP(WORK(KFIELD),WORK(KEND1),LWRK1,FF,1,LFIELD(IF))
         ENDDO
         CALL CC_FCKMO(WORK(KFIELD),WORK(KLAMP0),WORK(KLAMH0),
     *                 WORK(KEND1),LWRK1,1,1,1)
      ENDIF

*---------------------------------------------------------------------*
*     Loop over distributions of integrals:
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KINT1S),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2S),NT1AMX*NRHFT*NRHFT)

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
               CALL QUIT('Insufficient space in CCSDT_XI_CONT_NODDY')
            ENDIF
 
C           ---------------------------
C           Read in batch of integrals.
C           ---------------------------
            CALL CCRDAO(WORK(KXINT),IDEL,1,WORK(KEND2),LWRK2,
     *                  IRECNR,DIRECT)
 
C           ----------------------------------
C           Calculate integrals needed in CC3:
C           ----------------------------------
            CALL CCSDT_TRAN1(WORK(KINT1T),WORK(KINT2T),
     &                       WORK(KLAMP0),WORK(KLAMH0),
     *                       WORK(KXINT),IDEL)

            CALL CC3_TRAN2(WORK(KXIAJB),WORK(KYIAJB),
     &                     WORK(KLAMP0),WORK(KLAMH0),
     *                     WORK(KXINT),IDEL)
 
            CALL CCSDT_TRAN3(WORK(KINT1S),WORK(KINT2S),
     &                       WORK(KLAMP0),WORK(KLAMH0),
     *                       WORK(KXINT),IDEL)
 
         END DO   
      END DO  

*---------------------------------------------------------------------*
*     get zero-order cluster amplitudes:
*---------------------------------------------------------------------*
      KSCR1 = KEND1
      KT2AM = KSCR1 + NT1AMX
      KT3AM = KT2AM + NT1AMX*NT1AMX
      KEND1 = KT3AM + NT1AMX*NT1AMX*NT1AMX
      LWRK1 = LWORK - KEND1

      IF ((NONHF) .AND. (NFIELD .GT. 0)) THEN
        KT3AM2 = KEND1
        KEND1A = KT3AM2 + NT1AMX*NT1AMX*NT1AMX
      ELSE
        KEND1A = KEND1
      END IF

      LWRK1A = LWORK - KEND1A
      IF (LWRK1A .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_XI_CONT_NODDY')
      ENDIF

      IOPT = 2
      CALL CC_RDRSP('R0',0,ISYMOP,IOPT,MODEL,DUMMY,WORK(KT3AM))
      CALL CC_T2SQ(WORK(KT3AM),WORK(KT2AM),ISYMOP)

C     --------------------------------------------
C     Calculate the zero-order triples amplitudes:
C     --------------------------------------------
      CALL CCSDT_T03AM(WORK(KT3AM),WORK(KINT1S),WORK(KINT2S),
     *                 WORK(KT2AM),WORK(KSCR1),WORK(KFOCKD),
     *                 WORK(KFIELD),WORK(KT3AM2))

      IF (NONHF .AND. NFIELD.GT.0) THEN
         LUT3AM = -1
         CALL GPOPEN(LUT3AM,'T3AMPL','UNKNOWN',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND LUT3AM
         WRITE (LUT3AM) (WORK(KT3AM-1+IDX),IDX=1,NT1AMX*NT1AMX*NT1AMX)
         CALL GPCLOSE(LUT3AM,'KEEP')
      ENDIF

*---------------------------------------------------------------------*
*     precalculate density matrices:
*---------------------------------------------------------------------*
      KLC1AM = KEND1
      KLC2AM = KLC1AM + NT1AMX
      KLC3AM = KLC2AM + NT1AMX*NT1AMX
      KDIJ   = KLC3AM + NT1AMX*NT1AMX*NT1AMX
      KDAB   = KDIJ   + NRHFT*NRHFT
      KDIA   = KDAB   + NVIRT*NVIRT
      KDIA2  = KDIA   + NVIRT*NRHFT
      KMINT  = KDIA2  + NVIRT*NRHFT
      KEND2  = KMINT  + NT1AMX*NRHFT*NRHFT

      LWRK2  = LWORK - KEND2
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_XI_CONT_NODDY')
      ENDIF

   
      DO IDENS = 1, NXIDEN
     
        IDLSTL = IDXL_XIDEN(IDENS)
        ISYML  = ILSTSYM(LISTL,IDLSTL)

        IOPT = 3
        CALL CC_RDRSP(LISTL,IDLSTL,ISYML,IOPT,MODEL,
     *                WORK(KLC1AM),WORK(KLC3AM))
        CALL CC_T2SQ(WORK(KLC3AM),WORK(KLC2AM),ISYML)

        IF (LISTL(1:3).EQ.'L0 ') THEN

          CALL CCSDT_L03AM(WORK(KLC3AM),WORK(KINT1T),WORK(KINT2T),
     *                     WORK(KXIAJB),WORK(KFOCK0),WORK(KLC1AM),
     *                     WORK(KLC2AM),WORK(KSCR1),WORK(KFOCKD),
     *                     WORK(KFIELD),WORK(KT3AM))

          IF (NONHF .AND. NFIELD.GT.0) THEN
             ! restor T3 amplitudes which in case of finite difference 
             ! calculations have been overwritten in CCSDT_L03AM 
             LUT3AM = -1
             CALL GPOPEN(LUT3AM,'T3AMPL','UNKNOWN',' ','UNFORMATTED',
     &                   IDUMMY,.FALSE.)
             REWIND LUT3AM
             READ(LUT3AM)(WORK(KT3AM-1+IDX),IDX=1,NT1AMX*NT1AMX*NT1AMX)
             CALL GPCLOSE(LUT3AM,'KEEP')
          ENDIF

        ELSE IF (LISTL(1:3).EQ.'L1 ' .OR. LISTL(1:3).EQ.'LE ' .OR.
     &           LISTL(1:3).EQ.'M1 ' .OR. LISTL(1:3).EQ.'E0 ' .OR.
     &           LISTL(1:3).EQ.'N2 '                              )THEN

          CALL CCSDT_TBAR31_NODDY(WORK(KLC3AM),FREQL,LISTL,IDLSTL,
     &                          WORK(KLAMP0),WORK(KLAMH0),
     &                          WORK(KFOCK0),WORK(KFOCKD),WORK(KSCR1),
     &                          WORK(KXIAJB),WORK(KINT1T),WORK(KINT2T),
     &                          WORK(KEND2),LWRK2)

        ELSE
          WRITE(LUPRI,*) 'LISTL,IDLSTL:',LISTL(1:3),IDLSTL
          CALL QUIT('Unknown/Illegal list in CCSDT_XI_CONT_NODDY.')
        END IF


        CALL CCSDT_XI_DEN_NODDY(WORK(KLC2AM),WORK(KLC3AM),
     &                          WORK(KT2AM), WORK(KT3AM),
     &                          WORK(KDIJ),WORK(KDAB),
     &                          WORK(KDIA),WORK(KDIA2),
     &                          WORK(KMINT),NT1AMX,NVIRT,NRHFT)

        IF (EMULATE_CCSDPT_DENS2) THEN

          IRECORD = IDLSTL

c         -----------------------------------------------
c         write different blocks of the density matrix to 
c         special files as is also done by CCSDPT_DENS2:
c         (they will be read in by CC_D1AO)
c         -----------------------------------------------
          LUDEFF = -1
          IADRIA = 1 + NT1AM(ISYDEN)*(IRECORD-1)
          CALL WOPEN2(LUDEFF,FNDPTIA,64,0)  
          CALL PUTWA2(LUDEFF,FNDPTIA,WORK(KDIA),IADRIA,NT1AM(ISYDEN))
          CALL WCLOSE2(LUDEFF,FNDPTIA,'KEEP') 

          LUDEFF = -1
          IADRAB = 1 + NMATAB(ISYDEN)*(IRECORD-1)
          CALL WOPEN2(LUDEFF,FNDPTAB,64,0)  
          CALL PUTWA2(LUDEFF,FNDPTAB,WORK(KDAB),IADRAB,NMATAB(ISYDEN))
          CALL WCLOSE2(LUDEFF,FNDPTAB,'KEEP') 
 
          LUDEFF = -1
          IADRIJ = 1 + NMATIJ(ISYDEN)*(IRECORD-1)
          CALL WOPEN2(LUDEFF,FNDPTIJ,64,0)  
          CALL PUTWA2(LUDEFF,FNDPTIJ,WORK(KDIJ),IADRIJ,NMATIJ(ISYDEN))
          CALL WCLOSE2(LUDEFF,FNDPTIJ,'KEEP') 
 
          LUDEFF = -1
          IADRIA = 1 + NT1AM(ISYDEN)*(IRECORD-1)
          CALL WOPEN2(LUDEFF,FNDPTIA2,64,0)  
          CALL PUTWA2(LUDEFF,FNDPTIA2,WORK(KDIA2),IADRIA,NT1AM(ISYDEN))
          CALL WCLOSE2(LUDEFF,FNDPTIA2,'KEEP') 

        ELSE

c         ---------------------------------------------------------
c         write the blocks of the density matrix to a scratch file:
c         ---------------------------------------------------------

          IADRIJ = 1 + M2BST*(IDENS-1)
          CALL PUTWA2(LUDEFF,FNDEFF,WORK(KDIJ),IADRIJ,NMATIJ(ISYDEN))

          IADRAB = IADRIJ + NMATIJ(ISYDEN)
          CALL PUTWA2(LUDEFF,FNDEFF,WORK(KDAB),IADRAB,NMATAB(ISYDEN))

          CALL DAXPY(NT1AM(ISYDEN),1.0D0,WORK(KDIA2),1,WORK(KDIA),1)
          IADRIA = IADRAB + NMATAB(ISYDEN)
          CALL PUTWA2(LUDEFF,FNDEFF,WORK(KDIA),IADRIA,NT1AM(ISYDEN))

        END IF

      END DO

*---------------------------------------------------------------------*
*     calculate triples contributions to L x Xi{O} contractions:
*---------------------------------------------------------------------*
      IF (.NOT. EMULATE_CCSDPT_DENS2) THEN
       DO ITRAN = 1, NXETRAN
        IOPER  = IXETRAN(1,ITRAN)  ! operator index
        IRELAX = IXETRAN(5,ITRAN)  ! flag for relax. contrib.

        ISYHOP = ISYOPR(IOPER)     ! symmetry of O operator
        LABELH = LBLOPR(IOPER)     ! operator label
        LTWOEL = LPDBSOP(IOPER)    ! two-electron contribution to O?

        SKIPXI = ( IXETRAN(3,ITRAN) .EQ. -1 )

        LRELAX = LTWOEL .OR. (IRELAX.GE.1)
        IF (LRELAX) CALL QUIT(
     &    'Relaxed contributions not yet implemented in CCSDT_XI_CONT')

        IF (.NOT.SKIPXI) THEN
          IVEC = 1
          DO WHILE(IXDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
            IDLSTL = IXDOTS(IVEC,ITRAN)
            ISYML  = ILSTSYM(LISTL,IDLSTL)

            IF (ISYML.EQ.ISYHOP) THEN
                
              ! find index of density
              IDENS = 0
              DO IDX = 1, NXIDEN
                IF (IDLSTL.EQ.IDXL_XIDEN(IDX)) THEN
                  IDENS = IDX
                END IF
              END DO

              IF (IDENS.EQ.0) 
     &          CALL QUIT('Density not found in ccsdt_xi_cont_noddy')

              TRIPCON = CC_XIETA_DENCON(IDENS,LABELH,ISYHOP,
     &                     WORK(KLAMP0),WORK(KLAMH0),
     &                     LUDEFF,FNDEFF,M2BST,WORK(KEND1),LWRK1)

              IF (LOCDBG) THEN
                WRITE(LUPRI,*) 'CCSDT_XI_CONT_NODDY>'
                WRITE(LUPRI,*) 'IVEC,ITRAN,LABELH:',IVEC,ITRAN,LABELH
                WRITE(LUPRI,*) 'tbar3 x xi3 :',TRIPCON 
              END IF

              DOTPROD(IVEC,ITRAN) = DOTPROD(IVEC,ITRAN) + TRIPCON

            END IF 

            IVEC = IVEC + 1
          END DO

        END IF
       END DO ! ITRAN

       CALL WCLOSE2(LUDEFF,FNDEFF,'DELETE')
      END IF ! (.NOT. EMULATE_CCSDPT_DENS2)

      CALL QEXIT('CCSDT_XI_CONT_NODDY')
      RETURN
      END
*---------------------------------------------------------------------*
*             END OF SUBROUTINE CCSDT_XI_CONT_NODDY                   *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_XI_DEN_NODDY(T2BAR,T3BAR,T2AM,T3AM,
     &                              DIJ,DAB,DIA,DIA2,XMINT,
     &                              NT1AMX,NVIRT,NRHFT)
*---------------------------------------------------------------------*
*     Purpose: compute density corresponding to a contraction of      *
*              the triples parts of a Xi vector and (first-order)     *
*              lagrangian multipliers for a given multiplier vector   *
*                                                                     *
*        D^\xi_ab = 1/2 sum_delmn  tbar^dea_lmn t^deb_lmn             *
*        D^\xi_ij =-1/2 sum_delmf  tbar^def_lmj t^def_lmi             *
*        D^\xi_ia =-1/2 sum_delmfn tbar^def_lmn t^df_li t^ea_mn +     *
*                        sum_delm  tbar^de_lm (t^dea_lmi - t^dea_lim) *
*                                                                     *
*     Written by Christof Haettig, Nov 2002, Friedrichstal            *
*---------------------------------------------------------------------*
      IMPLICIT NONE  
#include <priunit.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER NT1AMX, NRHFT, NVIRT

#if defined (SYS_CRAY)
      REAL DAB(NVIRT,NVIRT), DIJ(NRHFT,NRHFT)
      REAL DIA(NVIRT,NRHFT), DIA2(NVIRT,NRHFT)
      REAL T3BAR(NT1AMX,NT1AMX,NT1AMX) T2BAR(NT1AMX,NT1AMX)
      REAL T3AM(NT1AMX,NT1AMX,NT1AMX), T2AM(NT1AMX,NT1AMX)
      REAL XMINT(NT1AMX,NRHFT,NRHFT), HALF, DDOT
#else
      DOUBLE PRECISION DAB(NVIRT,NVIRT), DIJ(NRHFT,NRHFT)
      DOUBLE PRECISION DIA(NVIRT,NRHFT), DIA2(NVIRT,NRHFT)
      DOUBLE PRECISION T3BAR(NT1AMX,NT1AMX,NT1AMX),T2BAR(NT1AMX,NT1AMX)
      DOUBLE PRECISION T3AM(NT1AMX,NT1AMX,NT1AMX), T2AM(NT1AMX,NT1AMX)
      DOUBLE PRECISION XMINT(NT1AMX,NRHFT,NRHFT), HALF, DDOT
#endif
      PARAMETER( HALF = 0.5D0 )

      INTEGER NDL, NEM, N, NXN, A, B, I, NXI, J, NXJ, C, F, M, E, NXM

      CALL DZERO(DIJ,NRHFT*NRHFT)
      CALL DZERO(DAB,NVIRT*NVIRT)
      CALL DZERO(DIA,NVIRT*NRHFT)
      CALL DZERO(DIA2,NVIRT*NRHFT)
      CALL DZERO(XMINT,NT1AMX*NRHFT*NRHFT)

      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'CCSDT_XI_DEN_NODDY>'
        WRITE(LUPRI,*) 'NORM^2(T3AM):',
     &     DDOT(NT1AMX*NT1AMX*NT1AMX,T3AM,1,T3AM,1)
        WRITE(LUPRI,*) 'NORM^2(T2AM):',
     &     DDOT(NT1AMX*NT1AMX,T2AM,1,T2AM,1)
        WRITE(LUPRI,*) 'NORM^2(T3BAR):',
     &     DDOT(NT1AMX*NT1AMX*NT1AMX,T3BAR,1,T3BAR,1)
      END IF

      DO NDL = 1, NT1AMX
        DO NEM = 1, NT1AMX

           DO N = 1, NRHFT
             NXN = NVIRT*(N-1)
             DO A = 1, NVIRT
               DO B = 1, NVIRT

                 DAB(A,B) = DAB(A,B) + 
     &             HALF * T3BAR(NDL,NEM,NXN+A) * T3AM(NDL,NEM,NXN+B)

               END DO
             END DO
           END DO


           DO I = 1, NRHFT
             NXI = NVIRT*(I-1)
             DO J = 1, NRHFT
               NXJ = NVIRT*(J-1)
               DO C = 1, NVIRT

                 DIJ(I,J) = DIJ(I,J) - 
     &             HALF * T3BAR(NDL,NEM,NXJ+C) * T3AM(NDL,NEM,NXI+C)
      
               END DO
             END DO
           END DO

           
           DO N = 1, NRHFT
             NXN = NVIRT*(N-1)
             DO I = 1, NRHFT
               NXI = NVIRT*(I-1)
               DO F = 1, NVIRT

                 XMINT(NEM,I,N) = XMINT(NEM,I,N) +
     &              T3BAR(NDL,NEM,NXN+F) * T2AM(NDL,NXI+F)
   
               END DO
             END DO
           END DO
   
          
        END DO

        DO M = 1, NRHFT
          NXM = NVIRT*(M-1)
          DO I = 1, NRHFT
            NXI = NVIRT*(I-1)
            DO E = 1, NVIRT
              DO A = 1, NVIRT
                DIA(A,I) = DIA(A,I) - T2BAR(NDL,NXM+E) * 
     &            (T3AM(NDL,NXM+E,NXI+A)-T3AM(NDL,NXI+E,NXM+A))
              END DO
            END DO
          END DO
        END DO

      END DO

   
      DO NEM = 1, NT1AMX
        DO N = 1, NRHFT
          NXN = NVIRT*(N-1)
          DO I = 1, NRHFT
            DO A = 1, NVIRT

              DIA2(A,I) = DIA2(A,I) + XMINT(NEM,I,N) * T2AM(NEM,NXN+A)

            END DO
          END DO
        END DO
      END DO

      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'NORM^2(DIJ):',DDOT(NRHFT*NRHFT,DIJ,1,DIJ,1)
        WRITE(LUPRI,*) 'NORM^2(DAB):',DDOT(NVIRT*NVIRT,DAB,1,DAB,1)
        WRITE(LUPRI,*) 'NORM^2(DIA):',DDOT(NVIRT*NRHFT,DIA,1,DIA,1)
        WRITE(LUPRI,*) 'NORM^2(DIA2):',DDOT(NVIRT*NRHFT,DIA2,1,DIA2,1)
        WRITE(LUPRI,*) 'NORM^2(XMINT):',
     &                   DDOT(NT1AMX*NVIRT*NRHFT,XMINT,1,XMINT,1)
      END IF

      RETURN
      END
*---------------------------------------------------------------------*
*             END OF SUBROUTINE CCSDT_XI_DEN_NODDY                    *
*---------------------------------------------------------------------*
