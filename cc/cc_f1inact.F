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
      SUBROUTINE CC_F1INACT(FOCK1,ISYFCK,IREAL,LABEL,CMO,WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: calculate the SCF F^[1] matrix which is equal to the 
*              inactive Fock matrix calculated from the derivative 
*              integrals multiplied with the density matrix (i.e.
*              projected to the occupied space)
*
*              FOCK1   --  the F^[1] matrix (output)
*              ISYFCK  --  symmetry of the F^[1] matrix (input)
*              LABEL   --  operator label for the perturbation (input)
*              CMO     --  HF orbital coefficients (input)
*
*     Christof Haettig 5-2-1999, restructured 24-5-99
*
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccexpfck.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
 
      INTEGER ISYM0, LUFCK
      PARAMETER( ISYM0 = 1 ) 

      CHARACTER*(8) LABEL
      INTEGER ISYFCK, LWORK

#if defined (SYS_CRAY)
      REAL FOCK1(*), CMO(*), WORK(LWORK)
      REAL ONE, TWO, ZERO, FTRACE, TEMP, SIGN
#else
      DOUBLE PRECISION FOCK1(*), CMO(*), WORK(LWORK)
      DOUBLE PRECISION HALF, ONE, TWO, ZERO, FTRACE, TEMP, SIGN
#endif
      PARAMETER( HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0, ZERO = 0.0D0 )

      INTEGER IFOCK, IADRF, IOPT, ISYM, IEXPV, IREAL
      INTEGER KFOCK1,KOVERLP,KEND1,LWRK1
      INTEGER INDF, INDE1, INDE2, KOFF1, IORBB, KSAB
      INTEGER ISYMA, ISYMI, KSAI, KSIA, IOCCI, IORBA, IVIRA, IVIRB

* external functions:
      INTEGER IEFFFOCK
      INTEGER IEXPECT

*---------------------------------------------------------------------*
*     set indeces for SCF eff. Fock matrices and expect. values
*---------------------------------------------------------------------*
      INDF  = 2
      INDE1 = 3
      INDE2 = 4

*---------------------------------------------------------------------*
*     allocate memory for effective Fock and a dummy overlap matrix:
*---------------------------------------------------------------------*
      KFOCK1  = 1
      KOVERLP = KFOCK1  + N2BST(ISYFCK)
      KEND1   = KOVERLP + N2BST(ISYM0)
      LWRK1   = LWORK   - KEND1

      IF (LWRK1 .LT. NBAST) THEN
         CALL QUIT('Insufficient work space in CC_F1INACT.')
      END IF
      
*---------------------------------------------------------------------*
*     read first-order effective Fock matrix from file
*---------------------------------------------------------------------*
      IFOCK = IEFFFOCK(LABEL,ISYM,1)
      IADRF = IADRFCK(INDF,IFOCK)

      LUFCK = -1
      CALL WOPEN2(LUFCK,FILFCKEFF,64,0)
      CALL GETWA2(LUFCK,FILFCKEFF,WORK(KFOCK1),IADRF,N2BST(ISYFCK))
      CALL WCLOSE2(LUFCK,FILFCKEFF,'KEEP')

      IF (LOCDBG) THEN
         FTRACE = ZERO
         IF (ISYFCK.EQ.1) THEN
            DO ISYM = 1, NSYM
               KOFF1 = KFOCK1 + IAODIS(ISYM,ISYM)
               DO I = 1, NBAS(ISYM)
                 FTRACE = FTRACE + WORK(KOFF1+NBAS(ISYM)*(I-1)+I-1)
               END DO
            END DO
         END IF
         IEXPV = IEXPECT(LABEL,ISYM)
         WRITE (LUPRI,*) 'LABEL:',LABEL
         WRITE (LUPRI,*) 'ISYFCK,IFOCK,IEXPV:',ISYFCK,IFOCK,IEXPV
         WRITE (LUPRI,*) 'FTRACE of read matrix:',FTRACE
         WRITE (LUPRI,*) 'one-electron expect:',EXPVALUE(INDE1,IEXPV)
         WRITE (LUPRI,*) 'two-electron expect:',EXPVALUE(INDE2,IEXPV)
         WRITE (LUPRI,*) 'CC_F1INACT> F^[1] matrix in SO basis:'
         CALL CC_PRONELAO(FOCK1,ISYFCK)
      END IF
      
*---------------------------------------------------------------------*
*     transform to MO using unit overlap matrix:
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KOVERLP),N2BST(ISYM0))
      DO ISYM = 1, NSYM
         KOFF1 = KOVERLP + IAODIS(ISYM,ISYM)
         DO I = 1, NBAS(ISYM)
           WORK(KOFF1+NBAS(ISYM)*(I-1)+I-1) = ONE
         END DO
      END DO

      ! transform effective Fock matrix to MO basis
      CALL CC_EFFCKMO(WORK(KFOCK1),ISYFCK,CMO,WORK(KOVERLP),
     &                WORK(KEND1),LWRK1)

      IF (LOCDBG) THEN
         FTRACE = ZERO
         IF (ISYFCK.EQ.1) THEN
            DO ISYM = 1, NSYM
               KOFF1 = KFOCK1 + IAODIS(ISYM,ISYM)
               DO I = 1, NBAS(ISYM)
                 FTRACE = FTRACE + WORK(KOFF1+NBAS(ISYM)*(I-1)+I-1)
               END DO
            END DO
         END IF
         WRITE (LUPRI,*) 'LABEL:',LABEL
         WRITE (LUPRI,*) 'ISYFCK:',ISYFCK
         WRITE (LUPRI,*) 'FTRACE of matrix generated in CC_EFCKMO:',
     &                    FTRACE
         WRITE (LUPRI,*) 'CC_F1INACT> F^[1] matrix in MO basis:'
         CALL CC_PRONELAO(FOCK1,ISYFCK)
      END IF

*---------------------------------------------------------------------*
*     put into result matrix and project onto occupied space:
*---------------------------------------------------------------------*
      ! initialize output vector
      CALL DZERO(FOCK1,N2BST(ISYFCK))

      ! add 2 times to F^[1] matrix
      CALL DAXPY(N2BST(ISYFCK),TWO,WORK(KFOCK1),1,FOCK1,1)

      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'CC_F1INACT> direct contribution '//
     &        'to F^[1] matrix:'
         CALL CC_PRONELAO(FOCK1,ISYFCK)
      END IF

      ! delete the virtual block of X
      IF (.FALSE.) THEN
      DO ISYMA = 1, NSYM
        ISYMI  = MULD2H(ISYFCK,ISYMA)
        DO IVIRA = 1, NVIRS(ISYMA)
        DO IOCCI = 1, NBAS(ISYMI)
          IORBA = NRHFS(ISYMA) + IVIRA
          KSAI = IAODIS(ISYMA,ISYMI) + (IOCCI-1)*NORBS(ISYMA) + IORBA
          KSIA = IAODIS(ISYMI,ISYMA) + (IORBA-1)*NORBS(ISYMI) + IOCCI
          FOCK1(KSAI) = ZERO
C
C         quick hack for imaginary operators and/or natural connection
C
          SIGN = DBLE(IREAL)
          TEMP = HALF * ( FOCK1(KSAI) + SIGN * FOCK1(KSIA) )
          FOCK1(KSAI) = TEMP
          FOCK1(KSIA) = TEMP * SIGN
C
        END DO
        DO IVIRB = 1, NBAS(ISYMI)
          IORBA = NRHFS(ISYMA) + IVIRA
          IORBB = NRHFS(ISYMI) + IVIRB
          KSAB = IAODIS(ISYMA,ISYMI) + (IORBB-1)*NORBS(ISYMA) + IORBA
          FOCK1(KSAB) = ZERO
        END DO
        END DO
      END DO
      END IF


      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'CC_F1INACT> final result for F^[1] matrix:',
     *                    LABEL
         WRITE (LUPRI,*) 'trace of F^[1] matrix:',FTRACE
         WRITE (LUPRI,*) IREAL,SIGN
         CALL CC_PRONELAO(FOCK1,ISYFCK)
         FTRACE = ZERO
         DO ISYM = 1, NSYM
            KOFF1 = IAODIS(ISYM,ISYM)
            DO I = 1, NBAS(ISYM)
              FTRACE = FTRACE + FOCK1(KOFF1+NBAS(ISYM)*(I-1)+I)
            END DO
         END DO
      END IF

      RETURN
      END
*======================================================================*
