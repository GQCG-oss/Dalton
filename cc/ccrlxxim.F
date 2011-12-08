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
      SUBROUTINE CCRLXXIM(XIM,ISYXIM,LABEL,LORX,LPDBS,FREQ,
     &                    CMO,WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: Calculate the X intermediate needed to calculate
*              the orbital relaxation contributions to second-
*              and higher-order derivatives
*
*              XIM     --  the X intermediate (output)
*              ISYXIM  --  symmetry of the X intermediate (input)
*              LABEL   --  operator label for the perturbation (input)
*              LORX    --  orbital relaxation flag (input)
*              LPDBS   --  flag for reorthogonalization contrib. (input)
*              FREQ    --  frequency of the perturbation (input)
*              CMO     --  HF orbital coefficients
*
*     Christof Haettig 5-2-1999
*
*     N.B.: this routine is not yet adapted for QMATP diff. from QMATH
*
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccexpfck.h>
#include <cc1dxfck.h>
#include <ccfro.h>
#include <ccroper.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
 
      INTEGER ISYM0, LUFCK
      PARAMETER( ISYM0 = 1 ) 
      CHARACTER LABEL0*(8)
      PARAMETER( LABEL0 = 'HAM0    ' )

      CHARACTER*(8) LABEL
      LOGICAL LORX, LPDBS
      INTEGER ISYXIM, LWORK

#if defined (SYS_CRAY)
      REAL FREQ, XIM(*), CMO(*), WORK(LWORK)
      REAL ONE, TWO, ZERO, FTRACE
#else
      DOUBLE PRECISION FREQ, XIM(*), CMO(*), WORK(LWORK)
      DOUBLE PRECISION ONE, TWO, ZERO, FTRACE
#endif
      PARAMETER( ONE = 1.0D0, TWO = 2.0D0, ZERO = 0.0D0 )

      CHARACTER*(10) MODEL
      LOGICAL NOKAPPA, new
      INTEGER IFOCK, IADRF, IOPT, IKAPPA, ISYM, IEXPV, IREAL
      INTEGER KFOCK1,KFOCK0,KOVERLP,KAPPA,KEND1,KSCR1,LWRK1
      INTEGER IFCK1, KRMAT, KQMATP, KQMATH, IOPER, INDF, INDE1, INDE2
      INTEGER ISYMA, ISYMI, ISYALP, ISYBET, ISYAL0, ISYBT0
      INTEGER KCMOQ, KFCKMO, KLAMDPQ, KLAMDHQ, KT1AM, KSCR0
      INTEGER NBASA, NBASB, NBSA0, NBSB0, NMOA, ISYM1, ISYM2
      INTEGER KOFF1, KOFF2, KOFF3, KOFF4, KOFF5, KOFF6, KOFF7
      INTEGER ICMO(8,8), NCMO(8), ICOUNT
      INTEGER KSAI, KSIA, IORBI, IORBA, IVIRA

* external functions:
      INTEGER IEFFFOCK
      INTEGER IEXPECT
      INTEGER IR1KAPPA
      INTEGER IROPER
      INTEGER I1DXFCK

*---------------------------------------------------------------------*
*     some settings:
*---------------------------------------------------------------------*
      INDF  = 1
      INDE1 = 1
      INDE2 = 2

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
*     allocate memory for some intermediates used localy:
*---------------------------------------------------------------------*
      KFOCK1  = 1
      KOVERLP = KFOCK1  + N2BST(ISYXIM)
      KEND1   = KOVERLP + N2BST(ISYM0)
      LWRK1   = LWORK   - KEND1

      IF (LWRK1 .LT. NBAST) THEN
         CALL QUIT('Insufficient work space in CCRLXXIM.')
      END IF
      
      ! initialize output vector
      CALL DZERO(XIM,N2BST(ISYXIM))

      ! read first-order effective Fock matrix from file
      IFOCK = IEFFFOCK(LABEL,ISYM,1)
      IEXPV = IEXPECT(LABEL,ISYM)
      IADRF = IADRFCK(INDF,IFOCK)

      LUFCK = -1
      CALL WOPEN2(LUFCK,FILFCKEFF,64,0)
      CALL GETWA2(LUFCK,FILFCKEFF,WORK(KFOCK1),IADRF,N2BST(ISYXIM))
      CALL WCLOSE2(LUFCK,FILFCKEFF,'KEEP')

      IF (LOCDBG) THEN
         FTRACE = ZERO
         IF (ISYXIM.EQ.1) THEN
            DO ISYM = 1, NSYM
               KOFF1 = KFOCK1 + IAODIS(ISYM,ISYM)
               DO I = 1, NBAS(ISYM)
                 FTRACE = FTRACE + WORK(KOFF1+NBAS(ISYM)*(I-1)+I-1)
               END DO
            END DO
         END IF
         WRITE (LUPRI,*) 'LABEL:',LABEL
         WRITE (LUPRI,*) 'ISYXIM,IFOCK,IEXPV:',ISYXIM,IFOCK,IEXPV
         WRITE (LUPRI,*) 'FTRACE of read matrix:',FTRACE
         WRITE (LUPRI,*) 'one-electron expect:',EXPVALUE(INDE1,IEXPV)
         WRITE (LUPRI,*) 'two-electron expect:',EXPVALUE(INDE2,IEXPV)
      END IF
      
      ! read zero-order overlap matrix from file and square up
      CALL RDONEL('OVERLAP ',.TRUE.,WORK(KEND1),NBAST)
      CALL CCSD_SYMSQ(WORK(KEND1),ISYM0,WORK(KOVERLP))


      ! transform effective Fock matrix to MO basis
      CALL CC_EFFCKMO(WORK(KFOCK1),ISYXIM,CMO,WORK(KOVERLP),
     &                WORK(KEND1),LWRK1)

      IF (LOCDBG) THEN
         FTRACE = ZERO
         IF (ISYXIM.EQ.1) THEN
            DO ISYM = 1, NSYM
               KOFF1 = KFOCK1 + IAODIS(ISYM,ISYM)
               DO I = 1, NBAS(ISYM)
                 FTRACE = FTRACE + WORK(KOFF1+NBAS(ISYM)*(I-1)+I-1)
               END DO
            END DO
         END IF
         WRITE (LUPRI,*) 'LABEL:',LABEL
         WRITE (LUPRI,*) 'ISYXIM:',ISYXIM
         WRITE (LUPRI,*) 'FTRACE of matrix generated in CC_EFCKMO:',
     &        FTRACE
      END IF

      ! add 2 times to X intermediate
      CALL DAXPY(N2BST(ISYXIM),TWO,WORK(KFOCK1),1,XIM,1)

      IF (LOCDBG) THEN
         WRITE (LUPRI,*)
     &        'CCRLXXIM> direct contribution to X intermediate:'
         CALL CC_PRONELAO(XIM,ISYXIM)
      END IF


      IF (LORX .OR. LPDBS) THEN

         IOPER   = IROPER(LABEL,ISYM)
         IREAL   = ISYMAT(IOPER)
         IFCK1   = I1DXFCK('HAM0    ','R1 ',LABEL,FREQ,ISYM)

         KFOCK0  = KEND1
         KFOCK1  = KFOCK0  + N2BST(ISYM0)
         KAPPA   = KFOCK1  + N2BST(ISYXIM)
         KRMAT   = KAPPA   + 2*NALLAI(ISYXIM)
         KQMATP  = KRMAT   + N2BST(ISYXIM)
         KQMATH  = KQMATP  + N2BST(ISYXIM)
         KSCR1   = KQMATH  + N2BST(ISYXIM)
         KEND1   = KSCR1   + N2BST(ISYXIM)
         LWRK1   = LWORK   - KEND1

         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CCRLXXIM.')
         END IF
      
         IF (LORX) THEN
            IKAPPA = IR1KAPPA(LABEL,FREQ,ISYM)
            CALL CC_RDHFRSP('R1 ',IKAPPA,ISYM,WORK(KAPPA))
         ELSE
            CALL DZERO(WORK(KAPPA),2*NALLAI(ISYXIM))
         END IF

         CALL CC_GET_RMAT(WORK(KRMAT),IOPER,1,ISYXIM,WORK(KEND1),LWRK1)
     
       new = .true.
       if (.not. new) then
         CALL DSCAL(N2BST(ISYXIM),-1.0d0,WORK(KRMAT),1)
       end if

         NOKAPPA = .FALSE.
         CALL CC_QMAT(WORK(KQMATP),WORK(KQMATH),WORK(KRMAT),WORK(KAPPA),
     &                IREAL,ISYXIM,NOKAPPA,CMO,WORK(KEND1),LWRK1)

       if (new) then
         DO ISYM1 = 1, NSYM
            ISYM2 = MULD2H(ISYM1,ISYXIM)
            KOFF1 = KQMATP + IAODIS(ISYM1,ISYM2)
            KOFF2 = KQMATH + IAODIS(ISYM2,ISYM1)
            CALL TRSREC(NBAS(ISYM1),NBAS(ISYM2),
     &                  WORK(KOFF1),WORK(KOFF2))
         END DO          
         CALL DSCAL(N2BST(ISYXIM),-ONE,WORK(KQMATH),1)
       end if

         IFOCK = IEFFFOCK(LABEL0,ISYM,1)
         IEXPV = IEXPECT(LABEL0,ISYM)
         IADRF = IADRFCK(INDF,IFOCK)
         CALL WOPEN2(LUFCK,FILFCKEFF,64,0)
         CALL GETWA2(LUFCK,FILFCKEFF,WORK(KFOCK0),IADRF,N2BST(ISYM0))
         CALL WCLOSE2(LUFCK,FILFCKEFF,'KEEP')

         IF (LOCDBG) THEN
            FTRACE = ZERO
            DO ISYM = 1, NSYM
               KOFF1 = KFOCK0 + IAODIS(ISYM,ISYM)
               DO I = 1, NBAS(ISYM)
                 FTRACE = FTRACE + WORK(KOFF1+NBAS(ISYM)*(I-1)+I-1)
               END DO
            END DO
            WRITE (LUPRI,*) 'LABEL:',LABEL0
            WRITE (LUPRI,*) 'ISYXIM,IFOCK,IEXPV:',ISYXIM,IFOCK,IEXPV
            WRITE (LUPRI,*) 'FTRACE of read matrix :',FTRACE
            WRITE (LUPRI,*) 'one-electron expect:',EXPVALUE(INDE1,IEXPV)
            WRITE (LUPRI,*) 'two-electron expect:',EXPVALUE(INDE2,IEXPV)
         END IF

         CALL CC_EFFCKMO(WORK(KFOCK0),ISYM0,CMO,WORK(KOVERLP),
     &                   WORK(KEND1),LWRK1)
  
         IF (LOCDBG) THEN
            FTRACE = ZERO
            DO ISYM = 1, NSYM
               KOFF1 = KFOCK0 + IAODIS(ISYM,ISYM)
               DO I = 1, NBAS(ISYM)
                 FTRACE = FTRACE + WORK(KOFF1+NBAS(ISYM)*(I-1)+I-1)
               END DO
            END DO
            WRITE (LUPRI,*) 'LABEL:',LABEL0
            WRITE (LUPRI,*) 'ISYXIM:',ISYM0
            WRITE (LUPRI,*)
     &           'FTRACE of matrix generated by CC_EFFCKMO:',FTRACE
         END IF

         ! calculate F^0 x Q, result is put into scr
         CALL CC_MMOMMO('N','N',ONE,WORK(KFOCK0),ISYM0,
     &                  WORK(KQMATH),ISYXIM,ZERO,WORK(KSCR1),ISYXIM)

         IF (LOCDBG) THEN
            WRITE (LUPRI,*)
     &           'CCRLXXIM> relax. contrib. 1 to X intermediate:'
            CALL CC_PRONELAO(WORK(KSCR1),ISYXIM)
            FTRACE = ZERO
            DO ISYM = 1, NSYM
               KOFF1 = KSCR1 + IAODIS(ISYM,ISYM)
               DO I = 1, NBAS(ISYM)
                 FTRACE = FTRACE + WORK(KOFF1+NBAS(ISYM)*(I-1)+I-1)
               END DO
            END DO
            WRITE (LUPRI,*) 'trace:',FTRACE
         END IF

         ! add to result matrix:
         CALL DAXPY(N2BST(ISYXIM),ONE,WORK(KSCR1),1,XIM,1)

         ! read contribution from 'one-index' transformed density
         ! from file
         IADRF = IADR1DXF(1,IFCK1)
         CALL WOPEN2(LUFCK,FIL1DXFCK,64,0)
         CALL GETWA2(LUFCK,FIL1DXFCK,WORK(KFOCK1),IADRF,N2BST(ISYXIM))
         CALL WCLOSE2(LUFCK,FIL1DXFCK,'KEEP')

         CALL CC_EFFCKMO(WORK(KFOCK1),ISYXIM,CMO,WORK(KOVERLP),
     &                   WORK(KEND1),LWRK1)

         IF (LOCDBG) THEN
            WRITE (LUPRI,*) 
     &           'CCRLXXIM> relax. contrib. 2 to X intermediate:'
            CALL CC_PRONELAO(WORK(KFOCK1),ISYXIM)
            FTRACE = ZERO
            DO ISYM = 1, NSYM
               KOFF1 = KFOCK1 + IAODIS(ISYM,ISYM)
               DO I = 1, NBAS(ISYM)
                 FTRACE = FTRACE + WORK(KOFF1+NBAS(ISYM)*(I-1)+I-1)
               END DO
            END DO
            WRITE (LUPRI,*) 'trace:',FTRACE
         END IF
  
         ! add to result matrix:
         CALL DAXPY(N2BST(ISYXIM),ONE,WORK(KFOCK1),1,XIM,1)
      END IF
 
      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'CCRLXXIM> final result for X intermediate:'
         CALL CC_PRONELAO(XIM,ISYXIM)
         FTRACE = ZERO
         DO ISYM = 1, NSYM
            KOFF1 = IAODIS(ISYM,ISYM)
            DO I = 1, NBAS(ISYM)
              FTRACE = FTRACE + XIM(KOFF1+NBAS(ISYM)*(I-1)+I)
            END DO
         END DO
         WRITE (LUPRI,*) 'trace of X intermediate:',FTRACE
      END IF

      RETURN
      END
*======================================================================*
