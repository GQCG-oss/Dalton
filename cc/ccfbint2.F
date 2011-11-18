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
      SUBROUTINE CCFBINT2(ITRAN,   LABELH, TA2AMP, 
     &                    DENSPKQ, ONEHQ,   FOCKQ, DENSQ,
     &                    DENSA,   DENSPKA, FOCKA,
     &                    DENSQA,  DENPKQA, FOCKQA,
     &                    XLAMDH,  ISYM0,
     &                    XLAMDHQ, ISYHOP,
     &                    XLAMDHA, ISYMTA,
     &                    XLAMHQA, ISYHTA,
     &                    FNBFDA,  LUBFDA, IADRBFA,  IADBFA,
     &                    FNBFDQA, LUBFDQA, IADRBFQA, IADBFQA,
     &                    LRELAX,  LTWOEL,  LNEWTA,   LNEWOP,
     &                    WORK,    LWORK)
*---------------------------------------------------------------------*
* Purpose:
*
*     Precalculate some intermediates for F^BT^A vector depending
*     on T^A and/or IOPER (No Zeta vectors required):
*     -- one electron part of B operator AO integrals (ONEHQ)
*     -- AO-FOCKQ (initialized with ONEHQ)
*     -- 
*     -- The packed densities for FOCK(Q,A,QA) intermediates
*     -- The effective density of the rho^BFA intermediate (FNBFDA)
*     -- The effective density of the rho^BFQA intermediate (FNBFDA)
*
*     BFA density only computed for LNEWTA and if
*     LTWOEL or LRELAX are set
*
*     BFQA density computed for LNEWTA or LNEWOP and LRELAX
*     (to be checked)
*
*     Fock, OneHam & density intermediates computed always
*
*     Sonia Coriani, February 1999. Based on CCXIINT1
*     
* The actual calculation of the Fock densities could be moved inside here!
* (OBS: the routine is not called for CCS!!!)
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccfield.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL LRELAX, LTWOEL, LZERO, LNEWTA, LNEWOP
      CHARACTER*(*) FNBFDA, FNBFDQA
      CHARACTER*(8) LABELH,LABTEST
      INTEGER ITRAN, ISYM0, ISYMTA,ISYHOP, ISYHTA,LWORK
      INTEGER LUBFDA,  IADBFA,  IADRBFA(MXCORB,*)
      INTEGER LUBFDQA, IADBFQA, IADRBFQA(MXCORB,*)
      
#if defined (SYS_CRAY)
      REAL DENSPKQ(*), ONEHQ(*), FOCKQ(*)
      REAL FOCKA(*), FOCKQA(*)
      REAL XLAMDH(*), XLAMDHQ(*)
      REAL XLAMDHA(*), XLAMHQA(*)
      REAL DENSQ(*), TA2AMP(*), WORK(*)
      REAL DENSA(*), DENSPKA(*)
      REAL DENSQA(*), DENPKQA(*)
      REAL ZERO, THREE, DUMMY
#else
      DOUBLE PRECISION DENSPKQ(*), ONEHQ(*), FOCKQ(*)
      DOUBLE PRECISION FOCKA(*), FOCKQA(*)
      DOUBLE PRECISION XLAMDH(*), XLAMDHQ(*)
      DOUBLE PRECISION XLAMDHA(*), XLAMHQA(*)
      DOUBLE PRECISION DENSQ(*), TA2AMP(*), WORK(*)
      DOUBLE PRECISION DENSA(*), DENSPKA(*)
      DOUBLE PRECISION DENSQA(*), DENPKQA(*)
      DOUBLE PRECISION ZERO, THREE, DUMMY
#endif
      PARAMETER (ZERO = 0.0D0, THREE = 3.0D0)

      CHARACTER MODEL*(10)
      INTEGER IOPT, IDEL, IDUMMY, IFIELD, IRREP, ISYM, IERR
      INTEGER LFOCKQMO

*---------------------------------------------------------------------*
* generate lower triangular packed density matrices for Fock densities:
*---------------------------------------------------------------------*
      CALL CC_DNSPK(DENSQ,DENSPKQ,ISYHOP)
c
      IF (LNEWTA) THEN
         CALL CC_DNSPK(DENSA,DENSPKA,ISYMTA)
      END IF
c
      CALL CC_DNSPK(DENSQA,DENPKQA,ISYHTA)

*---------------------------------------------------------------------*
* get AO one-electron integrals h^X (in ONEHQ)
*---------------------------------------------------------------------*
      IF ( LABELH(1:8) .EQ. 'HAM0    ' ) THEN

        CALL CCRHS_ONEAO(ONEHQ,WORK,LWORK)
*       for zeroth-order Hamiltonian add finite fields:
        DO IFIELD = 1, NFIELD
          CALL CC_ONEP(ONEHQ,WORK,LWORK,
     &                 EFIELD(IFIELD),ISYHOP,LFIELD(IFIELD)  ) 
        END DO

C       --------------------------------------------
C       scale the one-electron integrals with three:
C       --------------------------------------------
        IF (LRELAX) THEN
           CALL DSCAL(N2BST(ISYHOP),THREE,ONEHQ,1)
           WRITE (LUPRI,*) 'Warning: multiply ONEHQ with 3 ...'
        END IF

      ELSE IF ( LABELH(1:8) .EQ. 'DUMMYOP ' ) THEN
        CALL DZERO(ONEHQ,N2BST(ISYHOP))
      ELSE 
* check what ISYM is
        CALL CCPRPAO(LABELH,.TRUE.,ONEHQ,IRREP,ISYM,IERR,WORK,LWORK)
        IF (IERR.NE.0 .OR. IRREP.NE.ISYHOP) THEN
          CALL QUIT('CCFBINT2: error while reading operator '//LABELH)
        END IF

      END IF

*---------------------------------------------------------------------*
* initialize derivative AO Fock matrix with h^x integrals (FOCKQ)
* and the others FOCKA and FOCKQA with zero's
*---------------------------------------------------------------------*
c FOCKB reused in ccfbtaf, clean up possible exceeding space!!!
c
      LFOCKQMO = MAX(N2BST(ISYHOP),N2BST(ISYHTA))
      CALL DZERO(FOCKQ,LFOCKQMO)
      CALL DCOPY(N2BST(ISYHOP),ONEHQ,1,FOCKQ,1)
c
      CALL DZERO(FOCKA,N2BST(ISYMTA))
      CALL DZERO(FOCKQA,N2BST(ISYHTA))

*---------------------------------------------------------------------*
* calculate effective density matrices for the rho^BFA, rho^BFA inter-
* mediates:
*---------------------------------------------------------------------*
      IF (CCSD) THEN

* a) BFA-density: for every NEW T^A, written on file inside called routine

         IF (LNEWTA .AND. (LRELAX.OR.LTWOEL) ) THEN
            IOPT = 3
            CALL CC_BFDEN(TA2AMP, ISYMTA, DUMMY,  IDUMMY,
     *                    XLAMDH, ISYM0,  XLAMDH, ISYM0,
     *                    XLAMDHA, ISYMTA, DUMMY,  IDUMMY,
     *                    FNBFDA,LUBFDA,IADRBFA, IADBFA,
     *                    ITRAN, IOPT, .FALSE., WORK, LWORK)
         ELSE IF (LRELAX) THEN
            DO IDEL = 1, NBAST
              IADRBFA(IDEL,ITRAN) = IADRBFA(IDEL,ITRAN-1)
            END DO
c         ELSE
c            DO IDEL = 1, NBAST
c              !rho^BFA non calculated if NOT relaxed/twoel case
c              IADRBFA(IDEL,ITRAN) = -999999    
c            END DO
         END IF


* b) BFQA-density: for every new T^A or IOPER

         IF ((LNEWTA .OR. LNEWOP).AND.LRELAX) THEN
            IOPT = 7
            CALL CC_BFDEN(TA2AMP, ISYMTA, DUMMY,   IDUMMY,
     *                    XLAMDHQ,ISYHOP, XLAMDH,  ISYM0,
     *                    XLAMDHA,ISYMTA, XLAMHQA, ISYHTA,
     *                    FNBFDQA,LUBFDQA,IADRBFQA, IADBFQA,
     *                    ITRAN, IOPT, .FALSE., WORK, LWORK)
         ELSE IF (LRELAX) THEN
            DO IDEL = 1, NBAST
              IADRBFQA(IDEL,ITRAN) = IADRBFQA(IDEL,ITRAN-1)
            END DO
c         ELSE
c            DO IDEL = 1, NBAST
c              IADRBFQA(IDEL,ITRAN) = -999999
c            END DO
          END IF
c
      END IF

*---------------------------------------------------------------------*
* that's it; return:
*---------------------------------------------------------------------*
      RETURN

      END 
*=====================================================================*
