C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
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
*=====================================================================*
      SUBROUTINE CC_EFFCKMO(FOCK,ISYFCK,CMO,OVERLP,WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: transform an effective Fock matrix from the mixed
*              covariant AO / contravariant AO basis to the MO basis.
*              the transformation is done in place.
*
*              FOCK    --  effective Fock matrix
*              ISYFCK  --  symmetry of FOCK
*              CMO     --  orbital coefficient matrix
*              OVERLP  --  AO overlap matrix
*
*     Christof Haettig 6-2-1999
*
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit_dec.h>
#include <priunit.h>
#include <ccorb_dec.h>
#include <ccorb.h>
#include <ccsdsym_dec.h>
#include <ccsdsym.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
 
      INTEGER ISYM0
      PARAMETER( ISYM0 = 1) 

      INTEGER LWORK, ISYFCK

#if defined (SYS_CRAY)
      REAL FOCK(*), CMO(*), OVERLP(*), WORK(LWORK)
      REAL ONE, ZERO
#else
      DOUBLE PRECISION FOCK(*), CMO(*), OVERLP(*), WORK(LWORK)
      DOUBLE PRECISION ONE, ZERO
#endif
      PARAMETER( ZERO = 0.0D0, ONE = 1.0D0 )

      INTEGER ICMO(8,8), KSCR1, KSCR2, KEND1, LWRK1, ICOUNT, ISYM1,
     &        ISYM2, KOFF1, KOFF2, KOFF3, KOFF4, NBASP, NORBP, NBASR,
     &        ISYMP, ISYMR, ISYM

*---------------------------------------------------------------------*
*     memory allocation & initializations:
*---------------------------------------------------------------------*
      KSCR1 = 1
      KSCR2 = KSCR1 + N2BST(ISYFCK)
      KEND1 = KSCR2 + N2BST(ISYFCK)
      LWRK1 = LWORK - KEND1

      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient work space in CC_EFFCKMO.')
      END IF

      ! set symmetry offsets for CMO matrix:
      DO ISYM = 1, NSYM
         ICOUNT = 0
         DO ISYM2 = 1, NSYM
            ISYM1 = MULD2H(ISYM,ISYM2)
            ICMO(ISYM1,ISYM2) = ICOUNT
            ICOUNT = ICOUNT + NBAS(ISYM1)*NORBS(ISYM2)
         END DO
      END DO
      
      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'CC_EFFCKMO> effective Fock matrix '//
     &        'in ao/AO basis:'
         CALL CC_PRONELAO(FOCK,ISYFCK)
      END IF

*---------------------------------------------------------------------*
*     transform the leading index from the contravariant ao to the
*     covariant AO basis using the overlap matrix:
*---------------------------------------------------------------------*
      DO ISYMP = 1, NSYM

         ISYMR = MULD2H(ISYMP,ISYFCK)

         KOFF1 = IAODIS(ISYMP,ISYMP) + 1
         KOFF2 = IAODIS(ISYMP,ISYMR) + 1
         KOFF3 = KSCR1 + IAODIS(ISYMP,ISYMR)

         NBASP = MAX(1,NBAS(ISYMP))

         CALL DGEMM('N','N',NBAS(ISYMP),NBAS(ISYMR),NBAS(ISYMP),
     *              ONE,OVERLP(KOFF1),NBASP,FOCK(KOFF2),NBASP,
     *              ZERO,WORK(KOFF3),NBASP)
      END DO

      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'CC_EFFCKMO> effective Fock matrix '//
     &        'in AO/AO basis:'
         CALL CC_PRONELAO(WORK(KSCR1),ISYFCK)
      END IF

*---------------------------------------------------------------------*
*     transform the leading index from covariant AO to MO:
*---------------------------------------------------------------------*
      DO ISYMP = 1, NSYM

         ISYMR = MULD2H(ISYMP,ISYFCK)

         KOFF1 = ICMO(ISYMP,ISYMP) + 1
         KOFF3 = KSCR1 + IAODIS(ISYMP,ISYMR)
         KOFF4 = KSCR2 + IAODIS(ISYMP,ISYMR)

         NBASP = MAX(1,NBAS(ISYMP))
         NORBP = MAX(1,NORB(ISYMP))

         CALL DGEMM('T','N',NORB(ISYMP),NBAS(ISYMR),NBAS(ISYMP),
     *              ONE,CMO(KOFF1),NBASP,WORK(KOFF3),NBASP,
     *              ZERO,WORK(KOFF4),NORBP)
      END DO

*---------------------------------------------------------------------*
*     transform the second index from covariant AO to MO:
*---------------------------------------------------------------------*
      ! initialize output vector 
      CALL DZERO(FOCK,N2BST(ISYFCK))
      
      DO ISYMP = 1, NSYM

         ISYMR = MULD2H(ISYMP,ISYFCK)

         KOFF1 = ICMO(ISYMR,ISYMR) + 1
         KOFF2 = IAODIS(ISYMP,ISYMR) + 1
         KOFF4 = KSCR2 + IAODIS(ISYMP,ISYMR)

         NBASR = MAX(1,NBAS(ISYMR))
         NORBP = MAX(1,NORB(ISYMP))

         CALL DGEMM('N','N',NORB(ISYMP),NORB(ISYMR),NBAS(ISYMR),
     *              ONE,WORK(KOFF4),NORBP,CMO(KOFF1),NBASR,
     *              ZERO,FOCK(KOFF2),NORBP)
      END DO
 
*---------------------------------------------------------------------*
*     print to output & return:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'CC_EFFCKMO> effective Fock matrix '//
     &        'in the MO basis:'
         CALL CC_PRONELAO(FOCK,ISYFCK)
      END IF

      RETURN
      END
*======================================================================*
