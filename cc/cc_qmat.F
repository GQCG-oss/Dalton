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
      SUBROUTINE CC_QMAT(QMATP, QMATH, RMAT, XKAPPA,
     &                   IREAL, ISYMQ, NOKAPPA, CMO, WORK, LWORK )
*---------------------------------------------------------------------*
*
*     Purpose: build the Q matrices, which are basically the connection
*              matrix transformed to the MO representation minus 
*              the orbital relaxation matrix:
*
*              Q^h :  R   - kappa
*              Q^p :  R^* + kappa^T
*
*              QMATP  : Q^p matrices
*              QMATH  : Q^h matrices
*              RMAT   : orbital connection matrix in AO basis
*              XKAPPA : orbital relaxation vector (packed)
*              CMO    : orbital coefficient matrix, Sirius ordering
*
*              NOKAPPA = .TRUE.  --> skip contribution of kappa
*                                  (used f.x. for SCF kappa rhs vectors)
*
*     Christof Haettig, March 1999
*     generalized for non-symmetric R and Kappa in November 1999
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#  include <priunit_dec.h>
#  include <ccorb_dec.h>
#  include <ccsdsym_dec.h>
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccfro.h>
#include <ccsdsym.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      LOGICAL NOKAPPA
      INTEGER IREAL, ISYMQ, LWORK

#if defined (SYS_CRAY)
      REAL QMATP(*), QMATH(*), RMAT(*), XKAPPA(*), CMO(*)
      REAL WORK(LWORK), ONE, ZERO, SIGN
#else
      DOUBLE PRECISION QMATP(*), QMATH(*), RMAT(*), XKAPPA(*), CMO(*)
      DOUBLE PRECISION WORK(LWORK), ONE, ZERO, SIGN
#endif
      PARAMETER(ONE=1.0D0, ZERO=0.0D0)

      INTEGER ISYALP, ISYBET, NBASA, NBASB, NORBSA, IORBI, IORBA
      INTEGER KOFF1,KOFF2,KOFF3,KOFF4,KKAI,KRAI,KRIA,KSCR1,KSCR2,KEND
      INTEGER NCMO(8), ICMO(8,8), ISYM, ICOUNT, ISYM2, ISYM1, LEN

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
*     put the orbital relaxation (kappa) vector into the Q matrices:
*---------------------------------------------------------------------*
      IF (.NOT. NOKAPPA) THEN
        CALL CCKAPPASQ(QMATH,XKAPPA,ISYMQ,'N')
        CALL CCKAPPASQ(QMATP,XKAPPA,ISYMQ,'T')
      ELSE
        CALL DZERO(QMATH,N2BST(ISYMQ))
        CALL DZERO(QMATP,N2BST(ISYMQ))
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'debug CC_QMAT> NOKAPPA :',NOKAPPA
        WRITE (LUPRI,*) 'debug CC_QMAT> IREAL   :',IREAL
        WRITE (LUPRI,*) 'debug CC_QMAT> kappa matrix in MO:'
        CALL CC_PRONELAO(QMATH,ISYMQ)
        WRITE (LUPRI,*) 'debug CC_QMAT> kappa^T matrix in MO:'
        CALL CC_PRONELAO(QMATP,ISYMQ)
        WRITE (LUPRI,*) 'debug CC_QMAT> R matrix in AO:'
        CALL CC_PRONELAO(RMAT,ISYMQ)
      END IF

*---------------------------------------------------------------------*
*     transform the connection matrix R to MO using the CMO vector,
*     which is not resorted to the CC standard ordering. 
*     put the result into the Q matrices.
*
*     THE FOLLOWING HAS TO BE COUNTERCHECKED!
*     WHAT HAPPENS FOR FROZEN AND DELETED(!) ORBITALS!
*---------------------------------------------------------------------*
      DO ISYALP = 1, NSYM
         ISYBET = MULD2H(ISYALP,ISYMQ)

         KSCR1  = 1
         KSCR2  = KSCR1 +  NBAS(ISYALP)*NORBS(ISYBET)
         KEND   = KSCR2 + NORBS(ISYALP)*NORBS(ISYBET)

         IF ( KEND .GT. LWORK ) THEN
            CALL QUIT('Insufficient memory in CC_QMAT.')
         END IF

         NBASA  = MAX(NBAS(ISYALP),1)
         NBASB  = MAX(NBAS(ISYBET),1)
         NORBSA = MAX(NORBS(ISYALP),1)

         ! transform second index of R(alpha,beta) to MO
         KOFF1 = IAODIS(ISYALP,ISYBET) + 1
         KOFF2 = ICMO(ISYBET,ISYBET)   + 1
         CALL DGEMM('N','N',NBAS(ISYALP),NORBS(ISYBET),NBAS(ISYBET),
     &              ONE,RMAT(KOFF1),NBASA,CMO(KOFF2),NBASB,
     &              ZERO,WORK(KSCR1),NBASA)


         ! transform leading index of R(alpha,beta) to MO
         KOFF3 = ICMO(ISYALP,ISYALP)   + 1
         CALL DGEMM('T','N',NORBS(ISYALP),NORBS(ISYBET),NBAS(ISYALP),
     &              ONE,CMO(KOFF3),NBASA,WORK(KSCR1),NBASA,
     &              ZERO,WORK(KSCR2),NORBSA)

         LEN   = NORBS(ISYALP) * NORBS(ISYBET)
         KOFF4 = IAODIS(ISYALP,ISYBET) + 1

         ! Q^p(a,b) = R(a,b) - kappa(a,b)
         CALL DSCAL(LEN,-ONE,QMATH(KOFF4),1)
         CALL DAXPY(LEN,+ONE,WORK(KSCR2),1,QMATH(KOFF4),1)


         ! Q^p(a,b) = R^*(a,b) + kappa^T(a,b)
         SIGN = DBLE(IREAL)
         CALL DAXPY(LEN,SIGN,WORK(KSCR2),1,QMATP(KOFF4),1)

      END DO


      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'debug CC_QMAT> Q^h matrix in MO:'
        CALL CC_PRONELAO(QMATH,ISYMQ)
        WRITE (LUPRI,*) 'debug CC_QMAT> Q^p matrix in MO:'
        CALL CC_PRONELAO(QMATP,ISYMQ)
      END IF

      RETURN
      END
*=====================================================================*
