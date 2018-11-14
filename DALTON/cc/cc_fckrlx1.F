!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
C
*======================================================================
      SUBROUTINE CC_FCKRLX1(FOCK1,FOCK0,ISYFCK1,ISYFCK0,
     &                     XLAMDP0,XLAMDH0,ISYMP0,ISYMH0,
     &                     XLAMDP1,XLAMDH1,ISYMP1,ISYMH1,
     &                     LRELAX,WORK,LWORK)
*----------------------------------------------------------------------
*
*     Purpose: transform derivative AO fock matrix to MO basis and
*              add relaxation contributions coming from the 
*              derivatives of the transformation matrices
* 
*              if LRELAX, add relaxation contributions 
*
*              FOCK1 : derivative fock matrix, replaced on output
*              FOCK0 : zeroth-order fock matrix, unchanged on output
*
*   WARNING: Symmetry of result is NOT necessarily ISYFCK1!!!
*
*     Christof Haettig, July 1998
*
*     Generalized for pairs of XLAMDH and XLAMDP of different symmetry
*     Symmetry of FOCK1, FOCK0 in input must be specified outside
*     Sonia Coriani, September 1999
*     (symmetry tests to be added) 
*======================================================================
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL LRELAX
      INTEGER ISYMP0, ISYMP1, ISYMH0, ISYMH1, ISYPH0, ISYPH1
      INTEGER ISYFCK1,ISYFCK0
      INTEGER LWORK, ISYRES, KEND1, KSCR, LWRK1

#if defined (SYS_CRAY)
      REAL FOCK1(*), FOCK0(*), WORK(*)
      REAL XLAMDP0(*), XLAMDH0(*), XLAMDP1(*), XLAMDH1(*)
      REAL ONE, XNORM, DDOT, DNRM2
#else
      DOUBLE PRECISION FOCK1(*), FOCK0(*), WORK(*)
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*), XLAMDP1(*), XLAMDH1(*)
      DOUBLE PRECISION ONE, XNORM, DDOT, DNRM2
#endif
      PARAMETER(ONE=1.0D0)

 
*---------------------------------------------------------------------*
*       if debug flag set, print input matrices in AO:
*---------------------------------------------------------------------*
        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'CC_FCKRLX1> FOCK1 in AO:'
          CALL CC_PRFCKAO(FOCK1,ISYFCK1)
          XNORM = DNRM2(N2BST(ISYFCK1),FOCK1,1)
          WRITE (LUPRI,*) 'Norm of AO FOCK1 matrix:', XNORM
          WRITE (LUPRI,*) 'CC_FCKRLX1> FOCK0 in AO:'
          CALL CC_PRFCKAO(FOCK0,ISYFCK0)
          XNORM = DNRM2(N2BST(ISYFCK0),FOCK0,1)
          WRITE (LUPRI,*) 'Norm of AO FOCK0 matrix:', XNORM
          CALL FLSHFO(LUPRI)
        END IF

*---------------------------------------------------------------------*
*       transform derivative AO Fock matrix to MO using XLAMDP0/XLAMDH0
*---------------------------------------------------------------------*
        ISYPH0 = MULD2H(ISYMP0,ISYMH0)
        ISYPH1 = MULD2H(ISYMP1,ISYMH1)
        ISYRES = MULD2H(ISYFCK1,ISYPH0)
       
        IF (ISYRES.NE.ISYFCK1) THEN
          WRITE (LUPRI,*) 
     *        'Warning:ISYRES.NE.ISYFCK1. Replace FOCK1 in output?'
        END IF

*  FOCK1 must be allocated outside as MAX(N2BST(ISYFCK1),N2BST(ISYRES))!!

        CALL CC_FCKMO(FOCK1,XLAMDP0,XLAMDH0,
     *                WORK,LWORK,ISYFCK1,ISYMP0,ISYMH0)

*---------------------------------------------------------------------*
*       transform zero AO Fock matrix to MO using XLAMDP0/XLAMDH1
*                                                 XLAMDP1/XLAMDH0
*---------------------------------------------------------------------*
        IF (LRELAX) THEN
          KSCR  = 1                                         !contains Fock_MO
          KEND1 = KSCR + MAX(N2BST(ISYFCK0),N2BST(ISYRES))
          LWRK1 = LWORK - KEND1
         
          IF ( LWRK1 .LT. 0 ) THEN
            CALL QUIT('Insufficient work space in CC_FCKRLX.')
          END IF

*         duplicate zeroth-order AO Fock matrix in WORK:
          CALL DCOPY(N2BST(ISYFCK0),FOCK0,1,WORK(KSCR),1)

*         transform zeroth-order AO FOCK with XLAMDP1 and XLAMDH0,
*         and add to transformed derivative Fock matrix:
          CALL CC_FCKMO(WORK(KSCR),XLAMDP1,XLAMDH0,
     &                  WORK(KEND1),LWRK1,ISYFCK0,ISYMP1,ISYMH0)

          CALL DAXPY(N2BST(ISYRES),ONE,WORK(KSCR),1,FOCK1,1)

*         duplicate zeroth-order AO Fock matrix in WORK:
          CALL DCOPY(N2BST(ISYFCK0),FOCK0,1,WORK(KSCR),1)

*         transform zeroth-order AO FOCK with XLAMDP0 and XLAMDH1,
*         and add to transformed derivative Fock matrix:
          CALL CC_FCKMO(WORK(KSCR),XLAMDP0,XLAMDH1,
     &                  WORK(KEND1),LWRK1,ISYFCK0,ISYMP0,ISYMH1)

          CALL DAXPY(N2BST(ISYRES),ONE,WORK(KSCR),1,FOCK1,1)

        END IF

*---------------------------------------------------------------------*
*       print debug output and return:
*---------------------------------------------------------------------*

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'CC_FCKRLX1> FOCK1 in MO:'
c          CALL CC_PRFCKMO(FOCK1,ISYRES)
          XNORM = DDOT(N2BST(ISYRES),FOCK1,1,FOCK1,1)
          WRITE (LUPRI,*) 'Norm of MO FOCK1 matrix:', XNORM
          CALL FLSHFO(LUPRI)
        END IF

      RETURN

      END

*======================================================================
