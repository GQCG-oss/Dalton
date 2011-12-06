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
      SUBROUTINE CC_HFR1DEN(R1DEN,IOPER,IORDER,ISYMOPR,WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: calculate HF density one-index transformed with the
*              connection matrix for the operator specified by IOPER
*
*     Christof Haettig, October 1999
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
#include <dummy.h>
#include <inftap.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER LUSIRG

      INTEGER IOPER, ISYMOPR, IORDER, LWORK

#if defined (SYS_CRAY)
      REAL R1DEN(*), WORK(LWORK) 
      REAL TEMP, SIGN, TWO, ONE, HALF, ZERO
#else
      DOUBLE PRECISION R1DEN(*), WORK(LWORK) 
      DOUBLE PRECISION TEMP, SIGN, TWO, ONE, HALF, ZERO
#endif
      PARAMETER(HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0, ZERO = 0.0D0)

      CHARACTER*8 LABEL
      LOGICAL NOKAPPA
      INTEGER KCMO, KCMOHQ, KRMAT, KQMATP, KQMATH, KEND1, LWRK1
      INTEGER ICMO(8,8), NCMO(8), ISYM, ISYM1, ISYM2, ICOUNT
      INTEGER ISYALP, ISYBET, ISYMI, NBASA, NBASB, NORBSA
      INTEGER KOFF1, KOFF2, KOFF3, IDXAB, IDXBA, IREAL

*---------------------------------------------------------------------*
*     short cut for the test case 'HAM0    ':
*---------------------------------------------------------------------*
      CALL QENTER('CC_HFR1DEN')
      LABEL = LBLOPR(IOPER)
      IF (LABEL.EQ.'HAM0    ') THEN
         CALL DZERO(R1DEN,N2BST(ISYMOPR))
         CALL QEXIT('CC_HFR1DEN')
         RETURN
      END IF                       

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'CC_HFR1DEN> entered... LABEL:',LABEL
        WRITE (LUPRI,*) 'CC_HFR1DEN> LWORK',LWORK
      END IF

*---------------------------------------------------------------------*
*     allocate work space and initialize index array ICMO:
*---------------------------------------------------------------------*
C VB mar 2001
C     NCMO has to be initialized before it can be used to allocate
C     memory. The allocation of work space has been moved down below
C     the initialization.
C
      DO ISYM = 1, NSYM
         ICOUNT = 0
         DO ISYM2 = 1, NSYM
            ISYM1 = MULD2H(ISYM,ISYM2)
            ICMO(ISYM1,ISYM2) = ICOUNT
            ICOUNT = ICOUNT + NBAS(ISYM1)*NORBS(ISYM2)
         END DO
         NCMO(ISYM) = ICOUNT
      END DO                             

      KCMO   = 1
      KCMOHQ = KCMO   + NLAMDS
      KRMAT  = KCMOHQ + NCMO(ISYMOPR)
      KQMATP = KRMAT  + N2BST(ISYMOPR)
      KQMATH = KQMATP + N2BST(ISYMOPR)
      KEND1  = KQMATH + N2BST(ISYMOPR)
      LWRK1  = LWORK - KEND1

      IF (LWRK1 .LT. 0) THEN
         CALL STOPIT('CC_HFR1DEN',' ',KEND1,LWORK)
      END IF

*---------------------------------------------------------------------*
*     read SCF orbital coefficient matrix from file: 
*---------------------------------------------------------------------*
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     *            IDUMMY,.FALSE.)
      REWIND(LUSIFC)
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ(LUSIFC)
      READ(LUSIFC)
      READ(LUSIFC) (WORK(KCMO+I-1),I=1,NLAMDS)
      CALL GPCLOSE(LUSIFC,'KEEP')

*---------------------------------------------------------------------*
*     get MO connection matrix: R -> Q^h; (R^* -> Q^p)
*---------------------------------------------------------------------*
      ! get the connection matrix
      CALL CC_GET_RMAT(WORK(KRMAT),IOPER,IORDER,ISYMOPR,
     &                 WORK(KEND1),LWRK1)

      ! transform connection matrix to MO representation
      NOKAPPA = .TRUE.
      IREAL   = ISYMAT(IOPER)
      CALL CC_QMAT(WORK(KQMATP),WORK(KQMATH),WORK(KRMAT),DUMMY,
     &             IREAL,ISYMOPR,NOKAPPA,WORK(KCMO),WORK(KEND1),LWRK1)

*---------------------------------------------------------------------*
*     transform leading index to contravariant AO: CMOQ^h = CMO x Q^h
*---------------------------------------------------------------------*
      DO ISYALP = 1, NSYM
         ISYBET = MULD2H(ISYALP,ISYMOPR)

         NBASA  = MAX(NBAS(ISYALP),1)
         NORBSA = MAX(NORBS(ISYALP),1)

         KOFF1 = KCMO   + ICMO(ISYALP,ISYALP)
         KOFF2 = KQMATH + IAODIS(ISYALP,ISYBET)
         KOFF3 = KCMOHQ + ICMO(ISYALP,ISYBET)

         CALL DGEMM('N','N',NBAS(ISYALP),NORBS(ISYBET),NORBS(ISYALP),
     &              ONE,WORK(KOFF1),NBASA,WORK(KOFF2),NORBSA,
     &              ZERO,WORK(KOFF3),NBASA)

      END DO
           
*---------------------------------------------------------------------*
*     calculate R1DEN(alp,bet) = sum_i CMOQ(alp,i) x CMO(bet,i)
*---------------------------------------------------------------------*
      DO ISYALP = 1, NSYM

         ISYMI  = MULD2H(ISYALP,ISYMOPR)
         ISYBET = ISYMI

         NBASA  = MAX(NBAS(ISYALP),1)
         NBASB  = MAX(NBAS(ISYBET),1)

         KOFF1 = KCMOHQ + ICMO(ISYALP,ISYMI)
         KOFF2 = KCMO  + ICMO(ISYBET,ISYMI)
         KOFF3 = IAODIS(ISYALP,ISYBET) + 1

         CALL DGEMM('N','T',NBAS(ISYALP),NBAS(ISYBET),NRHFS(ISYMI),
     &              ONE,WORK(KOFF1),NBASA,WORK(KOFF2),NBASB,
     &              ZERO,R1DEN(KOFF3),NBASA)
      END DO

*---------------------------------------------------------------------*
*     Add second contribution sum_i CMO(bet,i) x CMOQ(alp,i)^*
*     by symmetrization / antisymmetrization
*---------------------------------------------------------------------*
      SIGN = DBLE(ISYMAT(IOPER))

      DO ISYALP = 1, NSYM
        ISYBET = MULD2H(ISYALP,ISYMOPR)
        IF      (ISYBET .GT. ISYALP) THEN
          DO A = 1, NBAS(ISYALP)
          DO B = 1, NBAS(ISYBET)
            IDXAB = IAODIS(ISYALP,ISYBET) + NBAS(ISYALP)*(B-1) + A
            IDXBA = IAODIS(ISYBET,ISYALP) + NBAS(ISYBET)*(A-1) + B
            TEMP  = SIGN * R1DEN(IDXAB) + R1DEN(IDXBA)
            R1DEN(IDXAB) = TEMP
            R1DEN(IDXBA) = TEMP * SIGN
          END DO
          END DO
        ELSE IF (ISYBET .EQ. ISYALP) THEN
          DO A = 1, NBAS(ISYALP)
          DO B = 1, A
            IDXAB = IAODIS(ISYALP,ISYBET)+NBAS(ISYALP)*(B-1)+A
            IDXBA = IAODIS(ISYBET,ISYALP)+NBAS(ISYBET)*(A-1)+B
            TEMP  = SIGN * R1DEN(IDXAB) + R1DEN(IDXBA)
            R1DEN(IDXAB) = TEMP
            R1DEN(IDXBA) = TEMP * SIGN
          END DO
          END DO
        END IF
      END DO                 
      
*---------------------------------------------------------------------*
*     that's it... print some debug output and return:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'R1DEN for operator ',LABEL,IOPER
        CALL CC_PRONELAO(R1DEN,ISYMOPR)
        CALL FLSHFO(LUPRI)
      END IF

      CALL QEXIT('CC_HFR1DEN')
      RETURN
      END
*=====================================================================*
