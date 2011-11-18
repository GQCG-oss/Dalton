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
*---------------------------------------------------------------------*
c/* Deck CC_XETST */
*=====================================================================*
       SUBROUTINE CC_XETST(WORK,LWORK)
*---------------------------------------------------------------------*
*
* Purpose : provide some tests for the CC_XIETA routine
*           for more detailed information of the tests see below
*
*           noddy implementation for programmers use only...
*           ... to switch between the different test or for changing
*           the parameters the program must be recompiled...
*
* Christof Haettig, 1998/99
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccorb.h"
#include "maxorb.h"
#include "ccroper.h"
#include "ccr1rsp.h"
#include "cco1rsp.h"
#include "ccx1rsp.h"
#include "ccfro.h"
#include "cclists.h"
#include "dummy.h"

* local parameters:
      CHARACTER MSGDBG*(18)
      PARAMETER (MSGDBG='[debug] CC_XETST> ')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER MXXETRAN, MXVEC
      PARAMETER (MXXETRAN = 20, MXVEC = 10)
      INTEGER LWORK
#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL FREQ, WRKDLM
      REAL XCONS(MXVEC,MXXETRAN), ECONS(MXVEC,MXXETRAN)
      REAL DDOT, PROPRSP, PROPAVE, FF
      REAL ZERO, ONE, TWO, FOUR, FIVE, FACUP
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION FREQ, FREQ1, FREQ2, WRKDLM
      DOUBLE PRECISION XCONS(MXVEC,MXXETRAN), ECONS(MXVEC,MXXETRAN)
      DOUBLE PRECISION DDOT, PROPRSP, PROPAVE, FF
      DOUBLE PRECISION ZERO, ONE, TWO, FOUR, FIVE, FACUP
#endif 
      PARAMETER (ZERO = 0.0D0, ONE  = 1.0D0, TWO = 2.0D0)
      PARAMETER (FOUR = 4.0D0, FIVE = 5.0D0)

      CHARACTER*(3) FILXI, FILETA, LISTL, APROXR12
      CHARACTER*(8) LABEL, LABEL1, LABEL2
      CHARACTER*(10) MODEL, CROUT
      LOGICAL LORX, LORX1, LORX2, LTWO, LTWO1, LTWO2
      INTEGER IOPTRES, N2VEC
      INTEGER IXETRAN(MXDIM_XEVEC,MXXETRAN), NXETRAN, IDXR1HF
      INTEGER IXDOTS(MXVEC,MXXETRAN), IEDOTS(MXVEC,MXXETRAN)
      INTEGER IADRLQ(MXXETRAN), IADRDQ(MXXETRAN), ISYAMP, IOPTCC2
      INTEGER ISYM0, ISYHOP, IOPT, IOPER, ITRAN, IVEC, IDUM, ISYETA
      INTEGER KEND0, ITEST, KEND1A, LWRK1A, IORDER, ILSTETA, ISIGN
      INTEGER KIT2DEL0, KIT2DELQ, KEND1, KDUM, KEND2, LWRK1, LWRK2
      INTEGER KXKSI1, KXKSI2, KRHS1, KRHS2, KT1AMP0, KT2AMP0,KT0AMP0
      INTEGER NDENSQ, NLAMDQ, IPRS, KOMEGA1, KOMEGA2, KSCR2, KT0
      INTEGER KCMO,KCMO0,KCMOPQ,KCMOHQ,IDLSTL,ISYCTR,IDXR1,INUM,IREAL
      INTEGER KKAPPA, KRMAT, KLEFT1, KLEFT2, KETA1, KETA2, KRHO1, KRHO2
      INTEGER IRELAX, IRELAX1, IRELAX2, ISYM1, ISYM2, KCTR1, KCTR2
      INTEGER KETA12


* external function:
      INTEGER ILSTSYM
      INTEGER IROPER
      INTEGER IR1TAMP
      INTEGER IR1KAPPA
      INTEGER IL1ZETA
      INTEGER IRHSR1
      INTEGER IETA1

*---------------------------------------------------------------------*
* set up information for rhs vectors for the different tests:
*---------------------------------------------------------------------*
*  number of simultaneous transformations:
      NXETRAN = 1
*---------------------------------------------------------------------*
*    some initialization
      LORX1 = .FALSE.
      LORX2 = .FALSE.
*---------------------------------------------------------------------*
*  test 1: use zeroth-order Hamiltonian as perturbation. this gives
*          a xi vector which is 5 x vector function (omega) and a
*          a eta vector which is 5 x (eta(0) + zeta(0) A). if the
*          cluster amplitude and lagrangian multiplier equations are
*          converged, both vectors should be zero.
*          (does not require that anything is done before CC_XETST)
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      ITEST = 1
      IORDER = 1
      LABEL  = 'HAM0    '
      LORX   = .TRUE.
      LTWO   = .TRUE.
      FREQ   = ZERO
      ISYHOP = 1
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
*
*  test 2: test the xi and eta vector for a non-relaxed one-elctron 
*          perturbation against the old implemenation.
*          (does only require that the integrals for the operator
*           are available on file)
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c     ITEST = 2
c     IORDER = 1
c     LABEL  = 'ZDIPLEN '
c     LORX   = .FALSE.
c     LTWO   = .FALSE.
c     FREQ   = ZERO
c     ISYHOP = 1
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
*
*  test 3: set differentiated integrals to zero and test only the
*          orbital relaxation & reorthogonalization part of the xi
*          and eta vectors.
*          (requires that the CPHF equations for the `reference' 
*           operator, specified here, have been solved and the Kappa
*           vector is available on disc)
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c     ITEST = 3
c     IORDER = 1
c     LABEL  = 'ZDIPLEN '
c     LORX   = .TRUE.
c     LTWO   = .FALSE.
c     FREQ   = ZERO
c     ISYHOP = 1
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
*
*  test 4: test the xi and eta vectors for a `orbital-relaxed'
*          one-electron perturbation.
*          (requires that the CPHF equations for the operator specified 
*           have been solved and the Kappa vector is available on disc)
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c     ITEST = 4
c     IORDER = 1
c     LABEL  = 'ZDIPLEN '
c     LORX   = .TRUE.
c     LTWO   = .FALSE.
c     FREQ   = ZERO
c     ISYHOP = 1
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
*
*  test 5: test the xi and eta vectors for a second-order perturbation
*          operator made up from a relaxed field with pert.-dep. basis
*          and an unrelaxed one-electron perturbation
*          (requires that the CPHF equations for the operator specified 
*           have been solved and the Kappa vector is available on disc)
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c     ITEST  = 5
c     IORDER = 2
c     LABEL1 = '1DHAM003'
c     LORX1  = .TRUE.
c     LTWO1  = .TRUE.
c     FREQ1  = ZERO
c     LABEL2 = 'ZDIPLEN '
c     LORX2  = .FALSE.
c     LTWO2  = .FALSE.
c     FREQ2  = ZERO
c     CALL CC_FIND_SO_OP(LABEL1,LABEL2,LABEL,ISYHOP,ISIGN,
c    &                   INUM,WORK,LWORK)
c     LORX   = ( LORX1 .AND. LORX2 )
c     LTWO   = ( LTWO1 .AND. LTWO2 )
c     FREQ   = FREQ1 + FREQ2
*---------------------------------------------------------------------*

      ! allow extensions in the vector lists for the next few lines
      LOPROPN = .TRUE.
      LO1OPN  = .TRUE.
      LX1OPN  = .TRUE.

      IRELAX = 0
      IF (LORX) THEN
        IRELAX = IR1TAMP(LABEL,LORX,FREQ,ISYHOP)
      END IF


      LPDBSOP(IROPER(LABEL,ISYHOP)) = LTWO

      ! set the IXETRAN array for one (XI,ETA) pair
      IXETRAN(1,1) = IROPER(LABEL,ISYHOP)
      IXETRAN(2,1) = 0
      IXETRAN(3,1) = IRHSR1(LABEL,LORX,FREQ,ISYHOP)
      IXETRAN(4,1) = IETA1(LABEL,LORX,FREQ,ISYHOP)
      IXETRAN(5,1) = IRELAX
 
      IF (IORDER.EQ.2) THEN
        IRELAX1 = 0
        IRELAX2 = 0
        IF (LORX1) IRELAX1 = IR1TAMP(LABEL1,LORX1,FREQ1,ISYM1)
        IF (LORX2) IRELAX2 = IR1TAMP(LABEL2,LORX2,FREQ2,ISYM2)
        IXETRAN(5,1) = IRELAX1
        IXETRAN(6,1) = IRELAX2
      END IF
 
      ! disallow again extension in the vector lists
      LOPROPN = .FALSE.
      LO1OPN  = .FALSE.
      LX1OPN  = .FALSE.
 
      ! for test 3, replace now the operator labels on the common blocks
      ! 'DUMMYOP ', interpreted inside of CC_XIETA as a zero one-electr.
      ! operator, but associated with a orb.-relax. (kappa) vector.
      IF (ITEST.EQ.3) THEN
C
         IDXR1   = IR1TAMP(LABEL,LORX,FREQ,ISYHOP)
         IDXR1HF = IR1KAPPA(LABEL,FREQ,ISYHOP)
         KT0   = 1
         KEND1 = KT0   + 2*NALLAI(ISYHOP)
         LWRK1 = LWORK - KEND1
         CALL CC_RDHFRSP('R1 ',IDXR1HF,ISYHOP,WORK(KT0))
C
         LABEL  = 'DUMMYOP '
         LBLOPR(IXETRAN(1,1)) = LABEL
         LBLO1(IXETRAN(3,1))  = LABEL
         LBLX1(IXETRAN(4,1))  = LABEL
         LRTLBL(IDXR1)        = LABEL
C
         IOPT = 4
         CALL CC_WRRSP('R1 ',IDXR1,ISYHOP,IOPT,MODEL,WORK(KT0),
     &                 DUMMY,DUMMY,WORK(KEND1),LWRK1)
C
      END IF
C
      CALL AROUND('CC_XETST: test of CC_XIETA subroutine')
      WRITE (LUPRI,*) 'ITEST  =',ITEST
      WRITE (LUPRI,*) 'LABEL  =',LABEL
      WRITE (LUPRI,*) 'LTWO   =',LTWO
      WRITE (LUPRI,*) 'LORX   =',LORX
      WRITE (LUPRI,*) 'FREQ   =',FREQ
      WRITE (LUPRI,*) 'IROPER =',IXETRAN(1,1)
      WRITE (LUPRI,*) 'ILEFT  =',IXETRAN(2,1)
      WRITE (LUPRI,*) 'IRHSR1 =',IXETRAN(3,1)
      WRITE (LUPRI,*) 'IETA1  =',IXETRAN(4,1)
      WRITE (LUPRI,*) 'IRELAX =',IXETRAN(5,1)
C
      DO I  = 2, NXETRAN
         IXETRAN(1,I) = IXETRAN(1,1)
         IXETRAN(2,I) = IXETRAN(2,1)
         IXETRAN(3,I) = IXETRAN(3,1) + I-1
         IXETRAN(4,I) = IXETRAN(4,1) + I-1
         IXETRAN(5,I) = IXETRAN(5,1)
         IXETRAN(6,I) = IXETRAN(6,1)
         LBLO1(IXETRAN(3,1)+I-1)  = LABEL
         ISYO1(IXETRAN(3,1)+I-1)  = ISYHOP
         LBLX1(IXETRAN(4,1)+I-1)  = LABEL
         ISYX1(IXETRAN(4,1)+I-1)  = ISYHOP
      END DO
      NO1LBL = MAX(NO1LBL,IXETRAN(3,1)+NXETRAN-1)
      NX1LBL = MAX(NX1LBL,IXETRAN(4,1)+NXETRAN-1)
C 
      N2VEC = 1
      IF (CCS) N2VEC = 0

*---------------------------------------------------------------------*
* call CC_XIETA to calculate the Xi and Eta vectors:
*---------------------------------------------------------------------*
      LISTL  = 'L0 '
      FILXI  = 'O1 '
      FILETA = 'X1 '
      IOPTRES = 3

      KEND0 = 1
      KEND1 = KEND0 + 1
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient work space in CC_XETST. (1)')
      END IF

      IDUM = 0
      WORK(KEND1-1) = WORK(IDUM)
      WRITE (LUPRI,*) 'WORK(0) before entering CC_XIETA:',WORK(KEND1-1)

      CALL CC_XIETA(IXETRAN, NXETRAN, IOPTRES, IORDER, LISTL, 
     &              FILXI,  IXDOTS, XCONS,
     &              FILETA, IEDOTS, ECONS,
     &              .FALSE.,0,  WORK(KEND1), LWRK1 )

      WRITE (LUPRI,*) 'returned from CC_XIETA... WORK(0) :',WORK(KEND1)
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
*     calculate the reference values xi vector and compare 
*     with the results from the CC_XIETA routine:
*---------------------------------------------------------------------*
      KKAPPA  = 1
      KRMAT   = KKAPPA  + 2*NALLAI(ISYHOP)
      KCMOPQ  = KRMAT   + N2BST(ISYHOP)
      KCMOHQ  = KCMOPQ  + NGLMDT(ISYHOP)
      KT1AMP0 = KCMOHQ  + NGLMDT(ISYHOP)
      KOMEGA1 = KT1AMP0 + NT1AMX
      KOMEGA2 = KOMEGA1 + NT1AM(ISYHOP)
      KRHS1   = KOMEGA2 + NT2AM(ISYHOP)
      KRHS2   = KRHS1   + NT1AM(ISYHOP)
      KEND1   = KRHS2   + NT2AM(ISYHOP)
      LWRK1   = LWORK   - KEND1

      KT2AMP0 = KOMEGA2 + MAX(NT2AMX,2*NT2ORT(1),NT2AO(1))
      KSCR2   = KT2AMP0 + NT2AMX
      KEND1A  = KSCR2   + NT2AMX + NT1AMX
      LWRK1A  = LWORK   - KEND1A

      KCTR1   = KEND1
      KCTR2   = KCTR1 + NT1AMX
      KEND2   = KCTR2 + NT2AMX
      LWRK2   = LWORK - KEND2

      IF (LWRK1.LT.0 .OR. LWRK1A.LT.0 .OR. LWRK2.LT.0) THEN
         CALL QUIT('Insufficient memory in CC_XETST.')
      END IF
 
C     ------------------------------
C     get orbital relaxation vector:
C     ------------------------------
      IF (LORX) THEN
        IF (LABEL.EQ.'HAM0    ') THEN
           CALL DZERO(WORK(KKAPPA),2*NALLAI(ISYHOP))
        ELSE
           IDXR1HF = IR1KAPPA(LABEL,FREQ,ISYHOP)
           CALL CC_RDHFRSP('R1 ',IDXR1HF,ISYHOP,WORK(KKAPPA))
        END IF
      ELSE
        CALL DZERO(WORK(KKAPPA),2*NALLAI(ISYHOP))
      END IF

C     ------------------------------
C     get orbital connection matrix:
C     ------------------------------
      IF (LORX.OR.LTWO) THEN
        IOPER = IROPER(LABEL,ISYHOP)
        CALL CC_GET_RMAT(WORK(KRMAT),IOPER,1,ISYHOP,WORK(KEND1),LWRK1)
      END IF 

C     ------------------------------------
C     construct the derivative CMO vector:
C     ------------------------------------
      IF (LORX.OR.LTWO) THEN
         IOPT  = 0
         IREAL = ISYMAT(IOPER)
         CALL CC_LAMBDAQ(DUMMY,DUMMY,WORK(KCMOPQ),WORK(KCMOHQ),ISYHOP,
     &                   DUMMY,WORK(KKAPPA),WORK(KRMAT),IREAL,IOPT,
     &                   WORK(KEND1),LWRK1)
         WRITE (LUPRI,*) 'CC_XIETA> CMOPQ vector:'
         CALL CC_PRONELAO(WORK(KCMOPQ),ISYHOP)
         WRITE (LUPRI,*) 'CC_XIETA> CMOHQ vector:'
         CALL CC_PRONELAO(WORK(KCMOHQ),ISYHOP)
      END IF

      IF (LORX.OR.LTWO) THEN
        IF (ITEST.EQ.1) THEN
C         ------------------------------------------------------------
C         call the vector function routine to calculate the reference 
C         result for the xi vector:
C         ------------------------------------------------------------
          IOPT = 3
          CALL CC_RDRSP('R0',0,1,IOPT,MODEL,
     &                  WORK(KT1AMP0),WORK(KT2AMP0))
          RSPIM = .FALSE.
          CALL CCRHSN(WORK(KOMEGA1),WORK(KOMEGA2),WORK(KT1AMP0),
     &                WORK(KT2AMP0),WORK(KEND1A),LWRK1A,'XXX')
          CALL DSCAL(NT1AMX,5.0D0,WORK(KOMEGA1),1)
          IF (.NOT.CCS) CALL DSCAL(NT2AMX,5.0D0,WORK(KOMEGA2),1)
          CALL AROUND('result from CCRHSN (x 5):')
        ELSE
C         ------------------------------------------------------------
C         evaluate orbital relaxation and reorthogonalization contrib.
C         to the xi vector by finite difference on the vector function 
C         (does only work for totally symmetric perturbations... 
C          actually only tested without symmetry...)
C         ------------------------------------------------------------
          IF (ISYHOP.NE.1) THEN
             WRITE (LUPRI,*) 'finite difference test of Xi vector not '
             WRITE (LUPRI,*) 'available for non-totally sym. '//
     &                       'perturbations.'
             CALL QUIT(
     &            'CC_FDXI only implemented for tot. sym. perturb.')
          END IF
          IF (IREAL.NE.1) THEN
             WRITE (LUPRI,*) 'finite difference test of Xi vector not '
             WRITE (LUPRI,*) 'available for imaginary perturbations.'
             CALL QUIT('CC_FDXI only implemented for real perturb.')
          END IF
          CALL CC_FDXI(WORK(KCMOHQ),WORK(KOMEGA1),WORK(KEND1),LWRK1)
          CALL AROUND('fin. diff. orbital relaxation Xi vector:')
        END IF
        CALL CC_PRP(WORK(KOMEGA1),WORK(KOMEGA2),ISYHOP,1,N2VEC)
      ELSE
         WRITE (LUPRI,*) 'no orb. relax. or reorth. contribution to Xi.'
         CALL DZERO(WORK(KOMEGA1),NT1AM(ISYHOP))
         IF (.NOT.CCS) CALL DZERO(WORK(KOMEGA2),NT2AM(ISYHOP))
      END IF

C     ----------------------------------------------------
C     calculate the contribution from one-electron h^1:
C     ----------------------------------------------------
      IF (LABEL(1:8).NE.'DUMMYOP '.AND.LABEL(1:8).NE.'HAM0    ') THEN
         IOPTCC2 = 0
         IF (LORX) IOPTCC2 = 1
         CALL CC_XKSI(WORK(KRHS1),LABEL,ISYHOP,IOPTCC2,DUMMY,
     &                WORK(KEND1),LWRK1)
         CALL AROUND('one-electron h^1 contribution to Xi vector:')
         CALL CC_PRP(WORK(KRHS1),WORK(KRHS2),ISYHOP,1,N2VEC)
         CALL DAXPY(NT1AM(ISYHOP),ONE,WORK(KRHS1),1,WORK(KOMEGA1),1)
         IF (.NOT.CCS) 
     &      CALL DAXPY(NT2AM(ISYHOP),ONE,WORK(KRHS2),1,WORK(KOMEGA2),1)
      END IF

      CALL AROUND('num. orb.-relax. + analyt. one-electr. Xi vector:')
      CALL CC_PRP(WORK(KOMEGA1),WORK(KOMEGA2),ISYHOP,1,N2VEC)

C     -------------------------------------------------
C     read the result vector from CC_XIETA and compare:
C     -------------------------------------------------
      DO ITRAN = 1, NXETRAN
        IOPT = 3
        IVEC = IXETRAN(3,ITRAN)
        WRITE (LUPRI,*) 'CC_XETST: IVEC = ',IVEC
        CALL CC_RDRSP('O1 ',IVEC,ISYHOP,IOPT,MODEL,
     &                WORK(KRHS1),WORK(KRHS2))

        CALL AROUND('analytical Xi vector:')
        CALL CC_PRP(WORK(KRHS1),WORK(KRHS2),ISYHOP,1,N2VEC)

        CALL DAXPY(NT1AM(ISYHOP),-1.0d0,WORK(KOMEGA1),1,WORK(KRHS1),1)
        IF (.NOT.CCS)
     &   CALL DAXPY(NT2AM(ISYHOP),-1.0d0,WORK(KOMEGA2),1,WORK(KRHS2),1)

        WRITE (LUPRI,*) 'Norm of difference between analytical Xi '
     &             // 'vector and the numerical result:'
        WRITE (LUPRI,*) 'singles excitation part:',
     &     DSQRT(DDOT(NT1AM(ISYHOP),WORK(KRHS1),1,WORK(KRHS1),1))
        IF (.NOT.CCS) WRITE (LUPRI,*) 'double excitation part: ',
     &     DSQRT(DDOT(NT2AM(ISYHOP),WORK(KRHS2),1,WORK(KRHS2),1))

        WRITE (LUPRI,*) 'difference vector:'
        CALL CC_PRP(WORK(KRHS1),WORK(KRHS2),ISYHOP,1,N2VEC)

C       -------------------------------------------------
C       test first-order property result:
C       -------------------------------------------------
        IF (ISYHOP.EQ.1) THEN
          IOPT = 3
          IVEC = IXETRAN(3,ITRAN)
          CALL CC_RDRSP('O1 ',IVEC,ISYHOP,IOPT,MODEL,
     &                  WORK(KRHS1),WORK(KRHS2))
          IOPT  = 3
          IVEC  = 0
          ISYM0 = 1
          CALL CC_RDRSP('L0',IVEC,ISYM0,IOPT,MODEL,
     &                  WORK(KCTR1),WORK(KCTR2))
          PROPRSP = DDOT(NT1AMX+NT2AMX,WORK(KRHS1),1,WORK(KCTR1),1)
          IOPER   = IROPER(LABEL,ISYHOP)
          IF (LORX.OR.LTWO .OR. IORDER.GT.1) THEN
            ILSTETA = IETA1(LABEL,LORX,FREQ,ISYHOP)
            PROPAVE = AVEX1(ILSTETA)
            WRITE (LUPRI,*) '...average contribution taken from '//
     &           'AVEX1...'
          ELSE
            FF = 1.0D0
            CALL CC_AVE(PROPAVE,LABEL,WORK(KEND1),LWRK1,FF)
            WRITE (LUPRI,*) '...average contribution calculated '//
     &                      'by CC_AVE...'
          END IF
          WRITE (LUPRI,*) ' OPERATOR : ',LABEL
          WRITE (LUPRI,*) ' AVERAGE  CONTRIBUTION :',PROPAVE
          WRITE (LUPRI,*) ' RESPONSE CONTRIBUTION :',PROPRSP
          WRITE (LUPRI,*) ' TOTAL ELECTRONIC CONTRIBUTION :',
     &         PROPRSP + PROPAVE
        END IF
    
      END DO


      IF (ITEST.EQ.1) THEN
          ! recalculate the response intermediates
          IOPT = 3
          CALL CC_RDRSP('R0',0,1,IOPT,MODEL,
     &                  WORK(KT1AMP0),WORK(KT2AMP0))
          RSPIM = .TRUE.
          CALL CCRHSN(WORK(KOMEGA1),WORK(KOMEGA2),WORK(KT1AMP0),
     &                WORK(KT2AMP0),WORK(KEND1A),LWRK1A,'XXX')
      END IF

*---------------------------------------------------------------------*
*     calculate the reference value for the eta vector and compare 
*     with the results from the CC_XIETA routine:
*---------------------------------------------------------------------*
      ISYCTR = ILSTSYM(LISTL,IDLSTL)
      ISYETA = MULD2H(ISYCTR,ISYHOP)

      KKAPPA  = 1
      KRMAT   = KKAPPA  + 2*NALLAI(ISYHOP)
      KCMOPQ  = KRMAT   + N2BST(ISYHOP)
      KCMOHQ  = KCMOPQ  + NGLMDT(ISYHOP)
      KLEFT1  = KCMOHQ  + NGLMDT(ISYHOP)
      KLEFT2  = KLEFT1  + NT1AM(ISYCTR)
      KRHO1   = KLEFT2  + NT2AM(ISYCTR)
      KRHO2   = KRHO1   + NT1AM(ISYETA)
      KETA1   = KRHO2   + NT2AM(ISYETA)
      KETA2   = KETA1   + NT1AM(ISYETA)
      KEND1   = KETA2   + NT2AM(ISYETA)
      if (CCR12) then
        keta12 = kend1
        kend1 = keta12 + ntr12am(isyeta)
      end if
      LWRK1   = LWORK   - KEND1

      IF (LWRK1.LT.0 .OR. LWRK1A.LT.0) THEN
         CALL QUIT('Insufficient memory in CC_XETST.')
      END IF
 
      IOPT   = 3
      Call CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,
     &              WORK(KLEFT1),WORK(KLEFT2))

      IF (LORX.OR.LTWO) THEN
         IF (ITEST.EQ.1) THEN
C          ------------------------------------------------------------
C          call the left transformation to calculate the reference 
C          result for the eta vector:
C          ------------------------------------------------------------
           CALL DCOPY(NT1AM(ISYCTR),WORK(KLEFT1),1,WORK(KRHO1),1)
           CALL DCOPY(NT2AM(ISYCTR),WORK(KLEFT2),1,WORK(KRHO2),1)
           CALL CC_ATRR(0.0D0,ISYCTR,-1,WORK(KRHO1),LWRK1,.FALSE.,DUMMY,
     &                  APROXR12)
           IF (LISTL(1:2).EQ.'L0') THEN
             CALL CC_ETA(WORK(KETA1),WORK(KEND1),LWRK1)
             CALL DAXPY(NT1AMX,1.0D0,WORK(KETA1),1,WORK(KRHO1),1)
             IF (.NOT.CCS) 
     &         CALL DAXPY(NT2AMX,1.0D0,WORK(KETA2),1,WORK(KRHO2),1)
           END IF
           CALL DSCAL(NT1AMX,5.0D0,WORK(KRHO1),1)
           IF (.NOT.CCS) CALL DSCAL(NT2AMX,5.0D0,WORK(KRHO2),1)
           CALL AROUND('result of CC_ETA & CC_LHTR (x 5):')
         ELSE
C          ---------------------------------------------------------
C          evaluate orbital relaxation and reorthogonalization contrib.
C          to the eta vector by finite difference on the jacobian
C          left transformation + zeroth-order eta vector.
C          (does only work for totally symmetric perturbations... 
C           actually only tested without symmetry...)
C          ---------------------------------------------------------
           IF (ISYETA.NE.1) THEN
             CALL QUIT
     *          ('CC_FDETA only implemented for total symmetric VECL.')
           END IF
           IF (IREAL.NE.1) THEN
             WRITE (LUPRI,*) 'finite difference test of Xi vector '//
     *                       'not available for imaginary'
             WRITE (LUPRI,*) 'perturbations.'
             CALL QUIT('CC_FDETA only implemented for real perturb.')
           END IF
           IF (ISYHOP.NE.1) THEN
             WRITE (LUPRI,*) 'finite difference test of Xi '//
     &             'vector not available'
             WRITE (LUPRI,*) 'for non-totally symmetric perturbations.'
             CALL QUIT('CC_FDETA only implemented for tot. '//
     &                 'sym. perturb.')
           END IF
           CALL CC_FDETA(WORK(KCMOHQ),LISTL,WORK(KLEFT1),WORK(KLEFT2),
     &                 WORK(KRHO1),WORK(KEND1),LWRK1)
           CALL AROUND('fin. diff. orbital relaxation Eta vector:')
         END IF
         CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYETA,1,N2VEC)
      ELSE
         WRITE (LUPRI,*) 'no orb. relax. or reorth. '//
     &        'contribution to Eta.'
         CALL DZERO(WORK(KRHO1),NT1AM(ISYETA))
         IF (.NOT.CCS) CALL DZERO(WORK(KRHO2),NT2AM(ISYETA))
      END IF

C     ----------------------------------------------------
C     calculate the contribution from one-electron h^1:
C     ----------------------------------------------------
      IF (LABEL(1:8).NE.'DUMMYOP '.AND.LABEL(1:8).NE.'HAM0    ') THEN
         IOPTCC2 = 0
         IF (LORX) IOPTCC2 = 1
         CALL CC_ETAC(ISYHOP,LABEL,WORK(KETA1),LISTL,IDLSTL,IOPTCC2,
     &                DUMMY,WORK(KEND1),LWRK1)
         CALL AROUND('one-electron h^1 contribution to Eta vector:')
         CALL CC_PRP(WORK(KETA1),WORK(KETA2),ISYETA,1,N2VEC)
         CALL DAXPY(NT1AM(ISYETA),ONE,WORK(KETA1),1,WORK(KRHO1),1)
         IF (.NOT.CCS) 
     &      CALL DAXPY(NT2AM(ISYETA),ONE,WORK(KETA2),1,WORK(KRHO2),1)
      END IF

      CALL AROUND('num. orb.-relax. + analyt. one-electr. Eta vector:')
      CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYETA,1,N2VEC)

C     -------------------------------------------------
C     read the result vector from CC_XIETA and compare:
C     -------------------------------------------------
      DO ITRAN = 1, NXETRAN
        IOPT = 3
        IVEC = IXETRAN(4,ITRAN) 
        WRITE (LUPRI,*) 'CC_XETST: IVEC = ',IVEC
        CALL CC_RDRSP('X1 ',IVEC,ISYETA,IOPT,MODEL,
     &                WORK(KETA1),WORK(KETA2))

        CALL AROUND('analytical Eta vector:')
        CALL CC_PRP(WORK(KETA1),WORK(KETA2),ISYETA,1,N2VEC)

        CALL DAXPY(NT1AM(ISYETA),-1.0d0,WORK(KRHO1),1,WORK(KETA1),1)
        IF (.NOT.CCS) 
     &    CALL DAXPY(NT2AM(ISYETA),-1.0d0,WORK(KRHO2),1,WORK(KETA2),1)

        WRITE (LUPRI,*) 'Norm of difference between analytical '
     &             // 'Eta vector and the numerical result:'
        WRITE (LUPRI,*) 'singles excitation part:',
     &     DSQRT(DDOT(NT1AM(ISYETA),WORK(KETA1),1,WORK(KETA1),1))
        IF (.NOT.CCS) WRITE (LUPRI,*) 'double excitation part: ',
     &     DSQRT(DDOT(NT2AM(ISYETA),WORK(KETA2),1,WORK(KETA2),1))

        WRITE (LUPRI,*) 'difference vector:'
        CALL CC_PRP(WORK(KETA1),WORK(KETA2),ISYETA,1,N2VEC)
      END DO

      RETURN
      END 
*=====================================================================*
*=====================================================================*
      SUBROUTINE CC_FDXI(CMO1,OMEGADIF,WORK,LWORK)
C
C---------------------------------------------------------------------
C Test routine for calculating the orbital relaxation contribution
C to the CC Xi vector by finite difference on the vector function.
C Ch. Haettig, januar 1999
C---------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <iratdef.h>
#include <ccorb.h>
#include <aovec.h>
#include <ccsdinp.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <inftap.h>
#include <ccinftap.h>
#include <leinf.h>
#include <ccfield.h>
C
      DIMENSION WORK(LWORK), OMEGADIF(*), CMO1(*)
      PARAMETER (DELTA = 1.0D-07, DELINV = 1.0D+7)
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0, TWO = 2.0d0, XHALF = 0.5D0)
      LOGICAL NONHFSAVE,LHTF
      CHARACTER MODEL*10
C
      IF (IPRINT.GT.10) THEN
         CALL AROUND( 'IN CC_FDXI  : MAKING FINITE DIFF. CC Xi Vector')
      END IF
C
      IF (CCR12) CALL QUIT('Finite-difference Xi-vector for CCR12 '//
     &                     'not adapted')
C
      IF (CC2) THEN
        IF (NFIELD.GT.0 .AND. (.NOT.NONHF)) THEN
          CALL QUIT('CC_FDXI incompatible with relaxed finite fields.')
        END IF
        NONHFSAVE = NONHF
        NONHF     = .TRUE.
      END IF
C
C----------------------------
C     Work space allocations.
C----------------------------
C
      ISYMTR     = 1
      ISYMOP     = 1
      ISYM0      = 1
C
      NTAMP      = NT1AMX + NT2AMX
      IF (CCS) NTAMP = NT1AMX
C
      KCMOREF    = 1
      KCMO       = KCMOREF  + NLAMDS
      KOMREF1    = KCMO     + NLAMDS
      KOMREF2    = KOMREF1  + NT1AMX
      KXIMAT     = KOMREF2  + NT2AMX
      KFCKDIF    = KXIMAT   + NTAMP  * NLAMDS
      KFCKREF    = KFCKDIF  + N2BST(ISYM0) 
      KFCKMAT    = KFCKREF  + N2BST(ISYM0)
      KEND0      = KFCKMAT  + N2BST(ISYM0) * NLAMDS
C
      KT1AMP0    = KEND0
      KOMEGA1    = KT1AMP0  + NT1AMX
      KOMEGA2    = KOMEGA1  + NT1AMX
      KT2AMP0    = KOMEGA2  + MAX(NT2AMX,2*NT2ORT(1),NT2AO(1))
      KSCR2      = KT2AMP0  + NT2AMX
      KEND1      = KSCR2    + NT2AMX + NT1AMX
      LWRK1      = LWORK    - KEND1
C
      IF (LWRK1.LT.0 ) THEN
         WRITE(LUPRI,*) 'Too little work space in CC_FDXI '
         WRITE(LUPRI,*) 'AVAILABLE: LWORK   =  ',LWORK
         WRITE(LUPRI,*) 'NEEDED (AT LEAST)  =  ',KEND1
         CALL QUIT('TOO LITTLE WORKSPACE IN CC_FDXI ')
      ENDIF
C
C---------------------
C     Initializations.
C---------------------
C
      X1 = 0.0D0
      X2 = 0.0D0
      CALL DZERO(OMEGADIF,NTAMP)
C
C---------------------------------------------
C     Read the reference CMO vector from disk:
C---------------------------------------------
C
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ(LUSIFC)
      READ(LUSIFC)
      READ(LUSIFC) (WORK(KCMOREF-1+I),I=1,NLAMDS)
      CALL GPCLOSE(LUSIFC,'KEEP') 
C
      CALL DZERO(WORK(KT1AMP0),NT1AMX)
      LHTF = .FALSE.
      CALL CCSD_IAJB(WORK(KT2AMP0),WORK(KT1AMP0),LHTF,
     &               .FALSE.,.FALSE.,WORK(KEND1),LWRK1)
      REWIND(LUIAJB)
      CALL WRITI(LUIAJB,IRAT*NT2AMX,WORK(KT2AMP0))
C
C----------------------------
C     Read the CC amplitudes:
C----------------------------
C
      IOPT = 3
      CALL CC_RDRSP('R0 ',0,1,IOPT,MODEL,WORK(KT1AMP0),WORK(KT2AMP0))
C
C---------------------------------------------------------
C     Calculate reference vector function and fock matrix:
C---------------------------------------------------------
C
      RSPIM = .FALSE.
      CALL CCRHSN(WORK(KOMEGA1),WORK(KOMEGA2),
     &            WORK(KT1AMP0),WORK(KT2AMP0),
     &            WORK(KEND1),LWRK1,'XXX')  
C
      IF (CCS) CALL DZERO(WORK(KOMEGA2),NT2AMX)
C
      OMEGA1N = DDOT(NT1AMX,WORK(KOMEGA1),1,WORK(KOMEGA1),1)
      OMEGA2N = DDOT(NT2AMX,WORK(KOMEGA2),1,WORK(KOMEGA2),1)
C
      IF (IPRINT.GT.10) THEN
        WRITE (LUPRI,*) 'CC_FDXI: norm of reference OMEGA1:',OMEGA1N
        WRITE (LUPRI,*) 'CC_FDXI: norm of reference OMEGA2:',OMEGA2N
      END IF
C
      CALL DCOPY(NT1AMX,WORK(KOMEGA1),1,WORK(KOMREF1),1)
      CALL DCOPY(NT2AMX,WORK(KOMEGA2),1,WORK(KOMREF2),1)
C
      LUFOCK = -1
      CALL GPOPEN(LUFOCK,'CC_FCKH','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUFOCK)
      READ(LUFOCK) (WORK(KFCKREF-1+I),I=1,N2BST(ISYM0))
      CALL GPCLOSE(LUFOCK,'KEEP') 
 
*----------------------------------------------------------------------*
*     Calculate the derivative of the vector function with respect
*     to the CMO vector by finite difference:
*----------------------------------------------------------------------*
 
      DO IDXDIF = 1, NLAMDS
C
C        -------------------------------
C        add finite displacement to CMO:
C        -------------------------------
C
         CALL DCOPY(NLAMDS,WORK(KCMOREF),1,WORK(KCMO),1)
         WORK(KCMO-1 + IDXDIF) = WORK(KCMO-1 + IDXDIF) + DELTA
C
         CMON = DDOT(NLAMDS,WORK(KCMO),1,WORK(KCMO),1)
C
         CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
         READ(LUSIFC)
         READ(LUSIFC)
         WRITE(LUSIFC) (WORK(KCMO-1+I),I=1,NLAMDS)
         CALL GPCLOSE(LUSIFC,'KEEP') 
C
         CALL DZERO(WORK(KT1AMP0),NT1AMX)
         LHTF = .FALSE.
         CALL CCSD_IAJB(WORK(KT2AMP0),WORK(KT1AMP0),
     &                   LHTF,.FALSE.,.FALSE.,WORK(KEND1),LWRK1)
         REWIND(LUIAJB)
         CALL WRITI(LUIAJB,IRAT*NT2AMX,WORK(KT2AMP0))
C
C        ------------------------------
C        calculate the vector function:
C        ------------------------------
C
         IOPT = 3
         CALL CC_RDRSP('R0 ',0,1,IOPT,MODEL,WORK(KT1AMP0),WORK(KT2AMP0))

         RSPIM = .FALSE.
         CALL CCRHSN(WORK(KOMEGA1),WORK(KOMEGA2),
     &               WORK(KT1AMP0),WORK(KT2AMP0),
     &               WORK(KEND1),LWRK1,'XXX')  
C
         IF (CCS) CALL DZERO(WORK(KOMEGA2),NT2AMX)
C
         CALL GPOPEN(LUFOCK,'CC_FCKH','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LUFOCK)
         READ(LUFOCK) (WORK(KFCKDIF-1+I),I=1,N2BST(ISYM0))
         CALL GPCLOSE(LUFOCK,'KEEP') 
C
C        --------------------------------------------------
C        construct the row nb. IDXDIF of the result matrix:
C        --------------------------------------------------
C
         OMEGA1N = DDOT(NT1AMX,WORK(KOMEGA1),1,WORK(KOMEGA1),1)
         OMEGA2N = DDOT(NT2AMX,WORK(KOMEGA2),1,WORK(KOMEGA2),1)
C
         IF (IPRINT.GT.10) THEN
            WRITE (LUPRI,*) 'CMO index:',IDXDIF
            WRITE (LUPRI,*) 'Norm of OMEGA1: ',OMEGA1N
            WRITE (LUPRI,*) 'Norm of OMEGA2: ',OMEGA2N
         END IF
C
         CALL DAXPY(NT1AMX,-1.0D0,WORK(KOMREF1),1,WORK(KOMEGA1),1)
         CALL DAXPY(NT2AMX,-1.0D0,WORK(KOMREF2),1,WORK(KOMEGA2),1)
         CALL DSCAL(NT1AMX,DELINV,WORK(KOMEGA1),1)
         CALL DSCAL(NT2AMX,DELINV,WORK(KOMEGA2),1)
C
         IF (IPRINT.GT.100) THEN
            CALL CC_PRP(WORK(KOMEGA1),WORK(KOMEGA2),1,1,1)
            WRITE (LUPRI,*) 'CMO1(index) = ', CMO1(IDXDIF)
            CALL DAXPY(NT1AMX,CMO1(IDXDIF),WORK(KOMEGA1),1,
     &                                     OMEGADIF(1),1)
            CALL DAXPY(NT2AMX,CMO1(IDXDIF),WORK(KOMEGA2),1,
     &                                     OMEGADIF(1+NT1AMX),1)
            WRITE (LUPRI,*) 'accumulated result vector:'
            CALL CC_PRP(OMEGADIF,OMEGADIF(NT1AMX+1),1,1,1)
         ENDIF
C
         KOFF1 = KXIMAT + NTAMP*(IDXDIF-1)
         KOFF2 = KOFF1  + NT1AMX
         CALL DCOPY(NT1AMX,WORK(KOMEGA1),1,WORK(KOFF1),1)
         CALL DCOPY(NT2AMX,WORK(KOMEGA2),1,WORK(KOFF2),1)
C
         X1 = X1 + DDOT(NT1AMX,WORK(KOMEGA1),1,WORK(KOMEGA1),1)
         X2 = X2 + DDOT(NT2AMX,WORK(KOMEGA2),1,WORK(KOMEGA2),1)
C
         CALL DAXPY(N2BST(ISYM0),-1.0D0,WORK(KFCKREF),1,WORK(KFCKDIF),1)
         CALL DSCAL(N2BST(ISYM0),1.0D0/DELTA,WORK(KFCKDIF),1)
C
         KOFF1 = KFCKMAT + N2BST(ISYM0)*(IDXDIF-1)
         CALL DCOPY(N2BST(ISYM0),WORK(KFCKDIF),1,WORK(KOFF1),1)
C
      END DO
C
      IF (IPRINT .GT. 10) THEN
         WRITE(LUPRI,*) '    '
         WRITE(LUPRI,*) '**  FINITE DIFF WITH DELTA ',DELTA, '**'
         WRITE(LUPRI,*) '    '
         IF (IPRINT .GT. 20) THEN
            CALL AROUND( 'FINITE DIFF. CC Xi-Matrix - 10 & 20 PART' )
            CALL OUTPUT(WORK(KXIMAT),1,NTAMP,1,NLAMDS,NTAMP,NLAMDS,1,
     &                  LUPRI)
         END IF
         XN = X1 + X2 
         WRITE(LUPRI,*) '  '
         WRITE(LUPRI,*) ' NORM OF 10 PART OF FD. Xi-mat.: ', SQRT(X1)
         WRITE(LUPRI,*) ' NORM OF 20 PART OF FD. Xi-mat.: ', SQRT(X2)
         WRITE(LUPRI,*) ' NORM OF  COMPLETE  FD. Xi-mat.: ', SQRT(XN)
      ENDIF
C
C-------------------------------------------
C     calculate Xi-matrix times CMO1 vector:
C-------------------------------------------
C
      CALL DGEMV('N',NTAMP,NLAMDS,ONE,WORK(KXIMAT),NTAMP,CMO1,1,
     *           ZERO,OMEGADIF,1)
C
C-----------------------------
C     scale diagonal with 1/2:
C-----------------------------
      CALL CCLR_DIASCL(OMEGADIF(NT1AM(1)+1),XHALF,1)

C
C----------------------------------------------
C     restore the reference CMO vector on disk:
C----------------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ(LUSIFC)
      READ(LUSIFC)
      WRITE(LUSIFC) (WORK(KCMOREF-1+I),I=1,NLAMDS)
      CALL GPCLOSE(LUSIFC,'KEEP') 
      IF (CC2) THEN
        NONHF = NONHFSAVE
      END IF
C
C------------------------------------------------------------------
C     make sure that all intermediates on file are calculated with
C     the reference CMO vector:
C------------------------------------------------------------------
C
C
      CALL DZERO(WORK(KT1AMP0),NT1AMX)
      LHTF = .FALSE.
      CALL CCSD_IAJB(WORK(KT2AMP0),WORK(KT1AMP0),LHTF,
     &               .FALSE.,.FALSE.,WORK(KEND1),LWRK1)
      REWIND(LUIAJB)
      CALL WRITI(LUIAJB,IRAT*NT2AMX,WORK(KT2AMP0))
C
      RSPIM = .TRUE.
      CALL CCRHSN(WORK(KOMEGA1),WORK(KOMEGA2),
     &            WORK(KT1AMP0),WORK(KT2AMP0),
     &            WORK(KEND1),LWRK1,'XXX')  
C
      IF (IPRINT.GT.10) THEN
         CALL AROUND(' END OF CC_FDXI ')
      END IF
C
      RETURN
      END
*=====================================================================*
*======================================================================*
      SUBROUTINE CC_FDETA(CMO1,LISTL,VECL1,VECL2,RHODIF,WORK,LWORK)
C
C----------------------------------------------------------------------
C Test routine for calculating the orbital relaxation contribution
C to the CC ETA{O} vector by finite difference on the jacob. left tran.
C Ch. Haettig, januar 1999
C----------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <iratdef.h>
#include <ccorb.h>
#include <aovec.h>
#include <ccsdinp.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccinftap.h>
#include <inftap.h>
#include <leinf.h>
#include <ccfield.h>
C
C------------------------------------------------------------
C     the displacement for the finite difference calculation:
C------------------------------------------------------------
      PARAMETER (DELTA = 1.0D-07, DELINV = 1.0D+7)
C------------------------------------------------------------
C
      DIMENSION WORK(LWORK), RHODIF(*), CMO1(*), VECL1(*), VECL2(*)
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0, TWO = 2.0d0, XHALF = 0.5D0)
      CHARACTER MODEL*10, LISTL*(*)
      LOGICAL NONHFSAVE,LHTF
C
      IF (IPRINT.GT.10) THEN
        CALL AROUND('IN CC_FDETA : MAKING FIN. DIFF. CC Eta{O} Vector')
      END IF
C
      IF (CC2) THEN
        IF (NFIELD.GT.0 .AND. (.NOT.NONHF)) THEN
          CALL QUIT('CC_FDETA incompatible with relaxed finite fields.')
        END IF
        NONHFSAVE = NONHF
        NONHF     = .TRUE.
      END IF
C
C----------------------------
C     Work space allocations.
C----------------------------
C
      ISYMTR     = 1
      ISYMOP     = 1
      ISYM0      = 1
C
      NTAMP      = NT1AMX + NT2AMX
      IF (CCS) NTAMP = NT1AMX
C
      KCMOREF    = 1
      KCMO       = KCMOREF  + NLAMDS
      KRHREF1    = KCMO     + NLAMDS
      KRHREF2    = KRHREF1  + NT1AMX
      KXIMAT     = KRHREF2  + NT2AMX
      KEND0      = KXIMAT   + NTAMP  * NLAMDS
      LWRK0      = LWORK    - KEND0
C
      KT1AMP0    = KEND0
      KOMEGA1    = KT1AMP0  + NT1AMX
      KOMEGA2    = KOMEGA1  + NT1AMX
      KT2AMP0    = KOMEGA2  + MAX(NT2AMX,2*NT2ORT(1),NT2AO(1))
      KSCR2      = KT2AMP0  + NT2AMX
      KEND1      = KSCR2    + NT2AMX + NT1AMX
      LWRK1      = LWORK    - KEND1
C
      KRHO1      = KEND0
      KRHO2      = KRHO1    + NT1AMX
      KEND1A     = KRHO2    + NT2AMX
      LWRK1A     = LWORK    - KEND1A
C
      KETA1      = KEND1A   
      KETA2      = KETA1    + NT1AMX
      KEND1B     = KETA2    + NT2AMX
      if (CCR12) then
        keta12   = kend1b
        kend1b   = keta12   + ntr12am(1)
      end if
      LWRK1B     = LWORK    - KEND1B
C
      IF ( LWRK1.LT.0 .OR. LWRK1B.LT.0 .OR. LWRK1A.LT.0) THEN
         WRITE(LUPRI,*) 'Too little work space in CC_FDETA '
         WRITE(LUPRI,*) 'AVAILABLE: LWORK   =  ',LWORK
         WRITE(LUPRI,*) 'NEEDED (AT LEAST)  =  ',MAX(KEND1,KEND1B)
         CALL QUIT('TOO LITTLE WORKSPACE IN CC_FDETA ')
      ENDIF
C
C---------------------
C     Initializations.
C---------------------
C
      X1 = 0.0D0
      X2 = 0.0D0
      CALL DZERO(RHODIF,NTAMP)
C
C---------------------------------------------
C     Read the reference CMO vector from disk:
C---------------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ(LUSIFC)
      READ(LUSIFC)
      READ(LUSIFC) (WORK(KCMOREF-1+I),I=1,NLAMDS)
      CALL GPCLOSE(LUSIFC,'KEEP') 
C
C-------------------------------------------------------------------
C     make sure that the correct response intermediates are on disc:
C-------------------------------------------------------------------
C
      CALL DZERO(WORK(KT1AMP0),NT1AMX)
      LHTF = .FALSE.
      CALL CCSD_IAJB(WORK(KT2AMP0),WORK(KT1AMP0),LHTF,
     &               .FALSE.,.FALSE.,WORK(KEND1),LWRK1)
      REWIND(LUIAJB)
      CALL WRITI(LUIAJB,IRAT*NT2AMX,WORK(KT2AMP0))
C
      IOPT = 3
      CALL CC_RDRSP('R0 ',0,1,IOPT,MODEL,WORK(KT1AMP0),WORK(KT2AMP0))
C
      RSPIM = .TRUE.
      CALL CCRHSN(WORK(KRHO1),WORK(KRHO2),
     &            WORK(KT1AMP0),WORK(KT2AMP0),
     &            WORK(KEND1),LWRK1,'XXX')  
C
C-------------------------------------------------------------------
C     calculate the reference vector:
C-------------------------------------------------------------------
C
      CALL DCOPY(NT1AMX,VECL1,1,WORK(KRHO1),1)
      IF (.NOT.CCS) CALL DCOPY(NT2AMX,VECL2,1,WORK(KRHO2),1)
C
      CALL CC_ATRR(0.0D0,ISYM0,-1,WORK(KRHO1),LWRK0,.FALSE.,DUMMY,
     &                  APROXR12)
C
      IF (LISTL(1:2).EQ.'L0') THEN
         CALL CC_ETA(WORK(KETA1),WORK(KEND1B),LWRK1B)
         CALL DAXPY(NT1AMX,ONE,WORK(KETA1),1,WORK(KRHO1),1)
         CALL DAXPY(NT2AMX,ONE,WORK(KETA2),1,WORK(KRHO2),1)
      END IF
C
      IF (CCS) CALL DZERO(WORK(KRHO2),NT2AMX)
C
      RHO1N = DDOT(NT1AMX,WORK(KRHO1),1,WORK(KRHO1),1)
      RHO2N = DDOT(NT2AMX,WORK(KRHO2),1,WORK(KRHO2),1)
C
      IF (IPRINT.GT.10) THEN
         WRITE (LUPRI,*) 'CC_FDETA: norm of reference RHO1:',RHO1N
         WRITE (LUPRI,*) 'CC_FDETA: norm of reference RHO2:',RHO2N
      END IF
C
      CALL DCOPY(NT1AMX,WORK(KRHO1),1,WORK(KRHREF1),1)
      CALL DCOPY(NT2AMX,WORK(KRHO2),1,WORK(KRHREF2),1)
 
*----------------------------------------------------------------------*
*     Calculate the derivative of the vector function with respect
*     to the CMO vector by finite difference:
*----------------------------------------------------------------------*
 
      DO IDXDIF = 1, NLAMDS
C
C        -------------------------------
C        add finite displacement to CMO:
C        -------------------------------
C
         CALL DCOPY(NLAMDS,WORK(KCMOREF),1,WORK(KCMO),1)
         WORK(KCMO-1 + IDXDIF) = WORK(KCMO-1 + IDXDIF) + DELTA
C
         CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
         READ(LUSIFC)
         READ(LUSIFC)
         WRITE(LUSIFC) (WORK(KCMO-1+I),I=1,NLAMDS)
         CALL GPCLOSE(LUSIFC,'KEEP') 
C
C        -------------------------------------
C        calculate new response intermediates:
C        -------------------------------------
C
         CALL DZERO(WORK(KT1AMP0),NT1AMX)
         LHTF = .FALSE.
         CALL CCSD_IAJB(WORK(KT2AMP0),WORK(KT1AMP0),
     &                  LHTF,.FALSE.,.FALSE.,WORK(KEND1),LWRK1)
         REWIND(LUIAJB)
         CALL WRITI(LUIAJB,IRAT*NT2AMX,WORK(KT2AMP0))
C
         IOPT = 3
         CALL CC_RDRSP('R0 ',0,1,IOPT,MODEL,WORK(KT1AMP0),WORK(KT2AMP0))
C
         RSPIM = .TRUE.
         CALL CCRHSN(WORK(KRHO1),WORK(KRHO2),
     &               WORK(KT1AMP0),WORK(KT2AMP0),
     &               WORK(KEND1),LWRK1,'XXX')  
C
C        ---------------------------------
C        calculate the transformed vector:
C        ---------------------------------
C
         CALL DCOPY(NT1AMX,VECL1,1,WORK(KRHO1),1)
         IF (.NOT.CCS) CALL DCOPY(NT2AMX,VECL2,1,WORK(KRHO2),1)
C
         CALL CC_ATRR(0.0D0,ISYM0,-1,WORK(KRHO1),LWRK0,.FALSE.,DUMMY,
     &                  APROXR12)
C
         IF (LISTL(1:2).EQ.'L0') THEN
            CALL CC_ETA(WORK(KETA1),WORK(KEND1B),LWRK1B)
            CALL DAXPY(NT1AMX,ONE,WORK(KETA1),1,WORK(KRHO1),1)
            CALL DAXPY(NT2AMX,ONE,WORK(KETA2),1,WORK(KRHO2),1)
         END IF
C
         IF (CCS) CALL DZERO(WORK(KRHO2),NT2AMX)
C
C        --------------------------------------------------
C        construct the row nb. IDXDIF of the result matrix:
C        --------------------------------------------------
C
         RHO1N = DDOT(NT1AMX,WORK(KRHO1),1,WORK(KRHO1),1)
         RHO2N = DDOT(NT2AMX,WORK(KRHO2),1,WORK(KRHO2),1)
C
         IF (IPRINT.GT.10) THEN
            WRITE (LUPRI,*) 'CMO index:',IDXDIF
            WRITE (LUPRI,*) 'Norm of RHO1: ',RHO1N
            WRITE (LUPRI,*) 'Norm of RHO2: ',RHO2N
         END IF
C
         CALL DAXPY(NT1AMX,-1.0D0,WORK(KRHREF1),1,WORK(KRHO1),1)
         CALL DAXPY(NT2AMX,-1.0D0,WORK(KRHREF2),1,WORK(KRHO2),1)
         CALL DSCAL(NT1AMX,DELINV,WORK(KRHO1),1)
         CALL DSCAL(NT2AMX,DELINV,WORK(KRHO2),1)
C
         IF (IPRINT.GT.100) THEN
            CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),1,1,1)
            WRITE (LUPRI,*) 'CMO1(index) = ', CMO1(IDXDIF)
            CALL DAXPY(NT1AMX,CMO1(IDXDIF),WORK(KRHO1),1,
     &                                     RHODIF(1),1)
            CALL DAXPY(NT2AMX,CMO1(IDXDIF),WORK(KRHO2),1,
     &                                     RHODIF(1+NT1AMX),1)
            WRITE (LUPRI,*) 'accumulated result vector:'
            CALL CC_PRP(RHODIF,RHODIF(NT1AMX+1),1,1,1)
         ENDIF
C
         KOFF1 = KXIMAT + NTAMP*(IDXDIF-1)
         KOFF2 = KOFF1  + NT1AMX
         CALL DCOPY(NT1AMX,WORK(KRHO1),1,WORK(KOFF1),1)
         CALL DCOPY(NT2AMX,WORK(KRHO2),1,WORK(KOFF2),1)
C
         X1 = X1 + DDOT(NT1AMX,WORK(KRHO1),1,WORK(KRHO1),1)
         X2 = X2 + DDOT(NT2AMX,WORK(KRHO2),1,WORK(KRHO2),1)
C
      END DO
C
      IF (IPRINT.GT.10) THEN
         WRITE(LUPRI,*) '    '
         WRITE(LUPRI,*) '**  FINITE DIFF WITH DELTA ',DELTA, '**'
         WRITE(LUPRI,*) '    '
         IF (IPRINT .GT. 20) THEN
            CALL AROUND( 'FINITE DIFF. CC Eta-Matrix - 10 & 20 PART' )
            CALL OUTPUT(WORK(KXIMAT),1,NTAMP,1,NLAMDS,NTAMP,NLAMDS,1,
     &                  LUPRI)
         ENDIF
         XN = X1 + X2 
         WRITE(LUPRI,*) '  '
         WRITE(LUPRI,*) ' NORM OF 10 PART OF FD. Eta-mat.: ', SQRT(X1)
         WRITE(LUPRI,*) ' NORM OF 20 PART OF FD. Eta-mat.: ', SQRT(X2)
         WRITE(LUPRI,*) ' NORM OF  COMPLETE  FD. Eta-mat.: ', SQRT(XN)
      END IF
C
C-------------------------------------------
C     calculate Eta-matrix times CMO1 vector:
C-------------------------------------------
C
      CALL DGEMV('N',NTAMP,NLAMDS,ONE,WORK(KXIMAT),NTAMP,CMO1,1,
     *           ZERO,RHODIF,1)
C
C----------------------------------------------
C     restore the reference CMO vector on disk:
C----------------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ(LUSIFC)
      READ(LUSIFC)
      WRITE(LUSIFC) (WORK(KCMOREF-1+I),I=1,NLAMDS)
      CALL GPCLOSE(LUSIFC,'KEEP') 
C
C------------------------------------------------------------------
C     make sure that all intermediates on file are calculated with
C     the reference CMO vector:
C------------------------------------------------------------------
C
      CALL DZERO(WORK(KT1AMP0),NT1AMX)
      LHTF = .FALSE.
      CALL CCSD_IAJB(WORK(KT2AMP0),WORK(KT1AMP0),LHTF,
     &               .FALSE.,.FALSE.,WORK(KEND1),LWRK1)
      REWIND(LUIAJB)
      CALL WRITI(LUIAJB,IRAT*NT2AMX,WORK(KT2AMP0))
C
      IOPT = 3
      CALL CC_RDRSP('R0 ',0,1,IOPT,MODEL,WORK(KT1AMP0),WORK(KT2AMP0))
C
      IF (CC2) THEN
        NONHF = NONHFSAVE
      END IF
C
      RSPIM = .TRUE.
      CALL CCRHSN(WORK(KRHO1),WORK(KRHO2),
     &            WORK(KT1AMP0),WORK(KT2AMP0),
     &            WORK(KEND1),LWRK1,'XXX')  
C
      IF (IPRINT.GT.10) THEN
        CALL AROUND(' END OF CC_FDETA ')
      END IF
C
      RETURN
      END
*=====================================================================*
