C
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren.
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
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>

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
