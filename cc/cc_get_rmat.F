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
      SUBROUTINE CC_GET_RMAT(RMAT, IOPER, ORDER, ISYHOP, WORK, LWORK )
*---------------------------------------------------------------------*
*
*     Purpose: retrieve the orbital connection matrix for the
*              perturbation operator specified in IOPER. 
*              
*              IOPER  -- operator index on IROPER/IROPER2 list
*              ORDER  -- derivative order of the operator
*              ISYHOP -- symmetry of operator and RMAT
*              
*              dimension of RMAT should be N2BST(ISYHOP)
*
*
*     Christof Haettig, March 1999
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccroper.h>
#include <ccropr2.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      INTEGER IOPER, ISYHOP, LWORK, ORDER

#if defined (SYS_CRAY)
      REAL RMAT(*), WORK(LWORK) 
      REAL HALF, ONE, ZERO
#else
      DOUBLE PRECISION RMAT(*), WORK(LWORK) 
      DOUBLE PRECISION HALF, ONE, ZERO
#endif
      PARAMETER(HALF=0.5D0, ONE=1.0D0, ZERO=0.0D0)

      CHARACTER*8 LABELH, LABELS
      LOGICAL LPDBS
      INTEGER ISYM, IRREP, IERR, ISYOPH

*---------------------------------------------------------------------*
* get operator label and check symmetry:
*---------------------------------------------------------------------*
      IF (ORDER.EQ.0 .OR. ORDER.EQ.1) THEN
         LABELH = LBLOPR(IOPER)
         ISYOPH = ISYOPR(IOPER)
         LPDBS  = LPDBSOP(IOPER)
         IF (LOCDBG) THEN
           WRITE(LUPRI,*) 'CC_GET_RMAT> ORDER:',ORDER
           WRITE(LUPRI,*) 'CC_GET_RMAT> IOPER,LABELH,ISYOPH,LPDBS:',
     &                                  IOPER,LABELH,ISYOPH,LPDBS
         END IF
      ELSE IF (ORDER.EQ.2) THEN
         LABELH = LBLOP2(IOPER,3)
         ISYOPH = ISYOP2(IOPER)
         LPDBS  = LPDBSOP2(IOPER)
         IF (LOCDBG) THEN
           WRITE(LUPRI,*) 'CC_GET_RMAT> ORDER:',ORDER
           WRITE(LUPRI,*) 'CC_GET_RMAT> IOPER,LABELH,ISYOPH,LPDBS:',
     &                                  IOPER,LABELH,ISYOPH,LPDBS
         END IF
      ELSE
         WRITE (LUPRI,*) 'CC_GET_RMAT> illegal operator order:',ORDER
         CALL QUIT('CC_GET_RMAT> illegal value for operator order.')
      END IF

      IF ( ISYHOP .NE. ISYOPH ) THEN
         WRITE (LUPRI,*) 'Symmetry mismatch in CC_GET_RMAT:'
         WRITE (LUPRI,*) 'Operator label:',LABELH
         WRITE (LUPRI,*) 'input symmetry:',ISYHOP
         WRITE (LUPRI,*) 'symmetry found:',ISYOPH
         CALL QUIT('Symmetry mismatch in CC_GET_RMAT.')
      END IF

*---------------------------------------------------------------------*
* case 1: basis set does not depend on this perturbation -->
*         connection matrix is zero
*---------------------------------------------------------------------*
      IF ( .NOT. LPDBS ) THEN

         CALL DZERO(RMAT,N2BST(ISYHOP))

*---------------------------------------------------------------------*
* case 2: test case 'HAM0    ' 
*         connection matrix is set to undifferentiated overlap matrix 
*---------------------------------------------------------------------*
      ELSE IF (LABELH.EQ.'HAM0    ') THEN

         IF (LWORK .LT. NBAST) THEN
            CALL QUIT('Insufficient work space in CC_GET_RMAT.')
         END IF
         CALL RDONEL('OVERLAP ',.TRUE.,WORK,NBAST)
         CALL CCSD_SYMSQ(WORK,ISYM0,RMAT) 

*---------------------------------------------------------------------*
* case 3: '1DHAM' first derivatives w.r.p to nuclear coordinates
*         connection matrix is set to 
*            a) symmetric connection : differentiated overlap matrix 
*            b) natural   connection : not yet available
*
*         (for some strange reason I (CH) can only get correct results
*          for the dipole gradient with symmetric connection if RMAT
*          is set to +1/2 S^(1) instead of -1/2 S^(1)...)
*
*---------------------------------------------------------------------*
      ELSE IF (LABELH(1:5).EQ.'1DHAM') THEN

         IF ( CONNECTION .EQ. 'SYMMETR' ) THEN
            WRITE(LABELS,'(A5,A3)') '1DOVL', LABELH(6:8)
            CALL CCPRPAO(LABELS,.TRUE.,RMAT,IRREP,ISYM,IERR,WORK,LWORK)
            IF (IERR.GT.0) THEN
              WRITE (LUPRI,'(A,A8,1X,A)') 
     &              'Warning:',LABELS,'Integrals missing!'
              WRITE (LUPRI,'(A,A8,1X,A)') 
     &              'Connection matrix for operator ',LABELH,' ignored.'
              CALL DZERO(RMAT,N2BST(ISYHOP))
            ELSE IF (IERR.LT.0) THEN
              CALL DZERO(RMAT,N2BST(ISYHOP))
            END IF
            CALL DSCAL(N2BST(ISYHOP),HALF,RMAT,1) 
         ELSE IF ( CONNECTION .EQ. 'NATURAL' ) THEN
            WRITE(LABELS,'(A6,A2)') 'SQHDOR',LABELH(7:8)
            CALL CCPRPAO(LABELS,.TRUE.,RMAT,IRREP,ISYM,IERR,WORK,LWORK)
            IF (IERR.GT.0) THEN
              WRITE (LUPRI,'(A,A8,1X,A)') 
     &            'Warning:',LABELS,'Integrals missing!'
              WRITE (LUPRI,'(A,A8,1X,A)') 
     &            'Connection matrix for operator ',LABELH,' ignored.'
              CALL DZERO(RMAT,N2BST(ISYHOP))
            ELSE IF (IERR.LT.0) THEN
              CALL DZERO(RMAT,N2BST(ISYHOP))
            END IF
            CALL DSCAL(N2BST(ISYHOP),ONE,RMAT,1)
            WRITE (LUPRI,*) 'Natural connection used for 1DHAM.'
         ELSE
            WRITE (LUPRI,*) 
     &         'Required connection not available for 1DHAM.'
            CALL QUIT('Required connection not available for 1DHAM.')
         END IF

*---------------------------------------------------------------------*
* case 4: 'dh/dB' first derivatives w.r.p to magnetic field
*         connection matrix is set to 
*            a) symmetric connection : differentiated overlap matrix 
*            b) natural   connection : not yet available
*
*---------------------------------------------------------------------*
      ELSE IF (LABELH(1:5).EQ.'dh/dB') THEN

         IF ( CONNECTION .EQ. 'SYMMETR' ) THEN
            WRITE(LABELS,'(A5,A3)') 'dS/dB',LABELH(6:8)
            CALL CCPRPAO(LABELS,.TRUE.,RMAT,IRREP,ISYM,IERR,WORK,LWORK)
            IF (IERR.GT.0) THEN
              WRITE (LUPRI,'(A,A8,1X,A)')
     &            'Warning:',LABELS,'Integrals missing!'
              WRITE (LUPRI,'(A,A8,1X,A)')
     &            'Connection matrix for operator ',LABELH,' ignored.'
              CALL DZERO(RMAT,N2BST(ISYHOP))
            ELSE IF (IERR.LT.0) THEN
              CALL DZERO(RMAT,N2BST(ISYHOP))
            END IF
            CALL DSCAL(N2BST(ISYHOP),-HALF,RMAT,1) 
C           WRITE (LUPRI,*) 'Symmetric connection for dh/dB used.'
         ELSE IF ( CONNECTION .EQ. 'NATURAL' ) THEN
            WRITE(LABELS,'(A7,A1)') 'd|S>/dB',LABELH(6:6)
            CALL CCPRPAO(LABELS,.TRUE.,RMAT,IRREP,ISYM,IERR,WORK,LWORK)
            IF (IERR.GT.0) THEN
              WRITE (LUPRI,'(A,A8,1X,A)')
     &            'Warning:',LABELS,'Integrals missing!'
              WRITE (LUPRI,'(A,A8,1X,A)')
     &            'Connection matrix for operator ',LABELH,' ignored.'
              CALL DZERO(RMAT,N2BST(ISYHOP))
            ELSE IF (IERR.LT.0) THEN
              CALL DZERO(RMAT,N2BST(ISYHOP))
            END IF
            CALL DSCAL(N2BST(ISYHOP),-ONE,RMAT,1) 
C           WRITE (LUPRI,*) 'Natural connection for dh/dB used.'
         ELSE
            WRITE (LUPRI,*) 
     &            'Required connection not available for dh/dB.'
            CALL QUIT('Required connection not available for dh/dB.')
         END IF

*---------------------------------------------------------------------*
* unknown operator: print error message and stop
*---------------------------------------------------------------------*
      ELSE
       WRITE (LUPRI,*) 'Error in CC_GET_RMAT:'
       WRITE (LUPRI,*) 'No connection matrix available for ',LABELH,
     &                 'operator.'
       CALL QUIT('Unknown operator/connection matrix in CC_GET_RMAT.')
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'CC_GET_RMAT> connection matrix R for ',LABELH
        WRITE (LUPRI,*) 'CC_GET_RMAT> connection used : ',CONNECTION
        CALL CC_PRONELAO(RMAT,ISYHOP)
      END IF

      RETURN
      END
*=====================================================================*
