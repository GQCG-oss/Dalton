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
*=====================================================================*
      SUBROUTINE CC_GET_LAMBDA1(IRELAX, IOPER, IOPT, XLAMDPQ, XLAMDHQ, 
     &                          CMOPQ,CMOHQ,ISYMQ,T1AMP0,WORK,LWORK )
*---------------------------------------------------------------------*
*
*     Purpose: Calculate the derivative LambdaQ matrices 
*              (for more information see CC_LAMBDAQ routine)
*
*     Christof Haettig, summer 1999
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccfro.h>
#include <ccsdsym.h>
#include <ccroper.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      INTEGER ISYMQ, LWORK, IOPER, IRELAX, IOPT

#if defined (SYS_CRAY)
      REAL XLAMDPQ(*), XLAMDHQ(*), CMOPQ(*), CMOHQ(*), T1AMP0(*)
      REAL WORK(LWORK) 
#else
      DOUBLE PRECISION XLAMDPQ(*),XLAMDHQ(*),CMOPQ(*),CMOHQ(*),T1AMP0(*)
      DOUBLE PRECISION WORK(LWORK) 
#endif

      CHARACTER*(10) MODEL
      INTEGER KAPPA, KRMAT, KEND, LWRK, IREAL

*---------------------------------------------------------------------*
*     memory allocation:
*---------------------------------------------------------------------*
      KAPPA = 1
      KRMAT = KAPPA + 2*NALLAI(ISYMQ)
      KEND  = KRMAT + N2BST(ISYMQ)
      LWRK  = LWORK - KEND

      IF (LWRK .LT. 0) THEN
         CALL QUIT('Insufficient work space in CC_GET_LAMBDA1.')
      END IF

      IF (LBLOPR(IOPER).EQ.'HAM0    ' .AND. IRELAX.GE.1) THEN
         WRITE (LUPRI,*) 'Test case "HAM0"... no relaxation '//
     &        'vector used.'
         CALL DZERO(WORK(KAPPA),2*NALLAI(ISYMQ))
      ELSE IF (IRELAX.GE.1) THEN
         CALL CC_RDHFRSP('R1 ',IRELAX,ISYMQ,WORK(KAPPA))
      ELSE
         CALL DZERO(WORK(KAPPA),2*NALLAI(ISYMQ))
      END IF

      CALL CC_GET_RMAT(WORK(KRMAT),IOPER,1,ISYMQ,WORK(KEND),LWRK)

      IOPT  = 1
      IREAL = ISYMAT(IOPER)
      CALL CC_LAMBDAQ(XLAMDPQ,XLAMDHQ,CMOPQ,CMOHQ,ISYMQ,
     &                T1AMP0,WORK(KAPPA),WORK(KRMAT),
     &                IREAL,IOPT,WORK(KEND),LWRK)
      

      RETURN
      END
*=====================================================================*
