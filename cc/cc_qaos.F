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
      SUBROUTINE CC_QAOS(QAOS,RMAT,XKAPPA,IREAL,ISYMQ,SAO,WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: calculate the product of Q^{p,ao} matrix with the 
*              AO overlap matrix:
*              
*         QAOS   -- result matrix: CMO Q^p CMO^T S^AO
*         RMAT   -- orbital connection matrix in AO basis
*         XKAPPA -- orbital relaxation vector in MO basis
*         IREAL  -- flag for real/imaginary R and kappa
*         ISYMQ  -- symmetry of XKAPPA, RMAT, and QAOS
*         SAO    -- overlap matrix
*
*     Christof Haettig, March 1999
*
*     N.B.: not yet fully adapted non-antisymmetric kappa
*           and/or non-symmetric R !!!
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <dummy.h>
#include <ccorb.h>
#include <ccsdsym.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      INTEGER IREAL, ISYMQ, LWORK

#if defined (SYS_CRAY)
      REAL QAOS(*), RMAT(*), XKAPPA(*), SAO(*), WORK(LWORK) 
      REAL ONE, ZERO
#else
      DOUBLE PRECISION QAOS(*), RMAT(*), XKAPPA(*), SAO(*), WORK(LWORK) 
      DOUBLE PRECISION ONE, ZERO
#endif
      PARAMETER(ONE=1.0D0, ZERO=0.0D0)

      LOGICAL NOKAPPA
      INTEGER ISYALP, ISYBET, ISYGAM, ISYMP, NBASA, NBASB
      INTEGER KQMOP, KQMOH, KCMOQ
      INTEGER KCMO, KQAO, KEND1, LWRK1, KOFF1, KOFF2, KOFF3, NORBSA
      INTEGER NCMO(8), ICMO(8,8), ISYM, ICOUNT, ISYM2, ISYM1

*---------------------------------------------------------------------*
*     set ICMO & NCMO arrays:
*---------------------------------------------------------------------*
      DO ISYM = 1, NSYM
         ICOUNT = 0
         DO ISYM2 = 1, NSYM
            ISYM1 = MULD2H(ISYM,ISYM2)
            ICMO(ISYM1,ISYM2) = ICOUNT
            ICOUNT = ICOUNT + NBAS(ISYM1)*NORBS(ISYM2)
         END DO
         NCMO(ISYM) = ICOUNT
      END DO 

*---------------------------------------------------------------------*
*     memory allocation:
*---------------------------------------------------------------------*
      KCMO  = 1
      KCMOQ = KCMO  + NLAMDS
      KQMOP = KCMOQ + NCMO(ISYMQ)
      KQMOH = KQMOP + N2BST(ISYMQ)
      KQAO  = KQMOH + N2BST(ISYMQ)
      KEND1 = KQAO  + N2BST(ISYMQ)
      LWRK1 = LWORK - KEND1

      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient work space in CC_QAOS.')
      END IF

*---------------------------------------------------------------------*
*     read (undifferentiated) MO coefficients from file:
*---------------------------------------------------------------------*
      CALL CC_GET_CMO(WORK(KCMO))

*---------------------------------------------------------------------*
*     build Q matrix in MO representation:
*---------------------------------------------------------------------*
      NOKAPPA = .FALSE.
      CALL CC_QMAT(WORK(KQMOP),WORK(KQMOH),RMAT,XKAPPA,
     &             IREAL,ISYMQ,NOKAPPA,WORK(KCMO),WORK(KEND1),LWRK1)

*---------------------------------------------------------------------*
*     transform to leading index to contravariant AO basis:
*             CMOQ = CMO x Q
*---------------------------------------------------------------------*
      DO ISYALP = 1, NSYM
         ISYBET = MULD2H(ISYALP,ISYMQ)

         NBASA  = MAX(NBAS(ISYALP),1)
         NORBSA = MAX(NORBS(ISYALP),1)

         KOFF1 = KCMO  + ICMO(ISYALP,ISYALP)
         KOFF2 = KQMOP + IAODIS(ISYALP,ISYBET)
         KOFF3 = KCMOQ + ICMO(ISYALP,ISYBET)

         CALL DGEMM('N','N',NBAS(ISYALP),NORBS(ISYBET),NORBS(ISYALP),
     &              ONE,WORK(KOFF1),NBASA,WORK(KOFF2),NORBSA,
     &              ZERO,WORK(KOFF3),NBASA)

      END DO  

*---------------------------------------------------------------------*
*     transform to second index to contravariant AO basis:
*             Q^ao = CMOQ x CMO^T
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KQAO),N2BST(ISYMQ))

      DO ISYALP = 1, NSYM

         ISYBET = MULD2H(ISYALP,ISYMQ)
         ISYMP  = ISYBET

         NBASA = MAX(NBAS(ISYALP),1)
         NBASB = MAX(NBAS(ISYBET),1)

         KOFF1 = KCMOQ + ICMO(ISYALP,ISYMP)
         KOFF2 = KCMO  + ICMO(ISYBET,ISYMP)
         KOFF3 = KQAO  + IAODIS(ISYALP,ISYBET)

         CALL DGEMM('N','T',NBAS(ISYALP),NBAS(ISYBET),NORBS(ISYMP),
     &              ONE,WORK(KOFF1),NBASA,WORK(KOFF2),NBASB,
     &              ZERO,WORK(KOFF3),NBASA)

      END DO

*---------------------------------------------------------------------*
*     multiply with the overlap matrix:
*---------------------------------------------------------------------*
      CALL CC_MAOMAO('N','N',ONE,WORK(KQAO),ISYMQ,SAO,ISYM0,
     &               ZERO,QAOS,ISYMQ)

      RETURN
      END
*=====================================================================*
