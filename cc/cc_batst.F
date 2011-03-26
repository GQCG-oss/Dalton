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
*---------------------------------------------------------------------*
c/* Deck CC_BATST */
*=====================================================================*
       SUBROUTINE CC_BATST(WORK,LWORK)
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccorb.h>

* local parameters:
      CHARACTER MSGDBG*(18)
      PARAMETER (MSGDBG='[debug] CC_BATST> ')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      LOGICAL NEW_CODE
      PARAMETER (NEW_CODE = .TRUE.)
      INTEGER MXCTRAN
      PARAMETER (MXCTRAN = 2)

      INTEGER LWORK
#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL DDOT, RDUM
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION DDOT, RDUM
#endif 

      LOGICAL LORXA, LORXB, LORXC
      CHARACTER*(3) LISTA, LISTB, LISTC
      CHARACTER*(8) FILCMA, LABELA, LABELB, LABELC
      CHARACTER*(10) MODEL
      INTEGER IOPTRES, IDUM
      INTEGER ICTRAN(4,MXCTRAN), NCTRAN
      INTEGER IDLSTA, IDLSTB, IDLSTC, ISYMA, ISYMB, ISYMC, ISYMABC
      INTEGER KTHETA1, KTHETA2, KT1AMPC, KT2AMPC, KRESLT1, KRESLT2
      INTEGER KEND1, LEND1, IOPT, IRELAXA

* external function:
      INTEGER IROPER
      INTEGER IR1TAMP
C     INTEGER IL1ZETA
      INTEGER ILSTSYM



*---------------------------------------------------------------------*
* call B{O} matrix transformation:
*---------------------------------------------------------------------*
      LISTA   = 'o1'
      LISTB   = 'R1'
      LISTC   = 'R1'
      LABELA  = 'ZDIPLEN '
      LABELB  = 'ZDIPLEN '
      LABELC  = 'ZDIPLEN '
      LORXA   = .FALSE.
      LORXB   = .FALSE.
      LORXC   = .FALSE.
      IDLSTA  = IROPER(LABELA,ISYMA)
      IDLSTB  = IR1TAMP(LABELB,LORXB,0.0D0,ISYMB)
      IDLSTC  = IR1TAMP(LABELC,LORXC,0.0D0,ISYMC)
      IRELAXA = 0
      IF (LORXA) IRELAXA = IR1TAMP(LABELA,LORXA,0.0D0,ISYMA)

      ICTRAN(1,1) = IDLSTA
      ICTRAN(2,1) = IDLSTB
      ICTRAN(3,1) = IDLSTC
      NCTRAN = 1

      IOPTRES = 1
      FILCMA  = 'CCCMAT'

      CALL CC_BAMAT(ICTRAN,  NCTRAN,  LISTA,  LISTB, LISTC,
     &              IOPTRES, FILCMA, IDUM, RDUM, 0, WORK, LWORK )


      ISYMA  = ILSTSYM(LISTA,IDLSTA)
      ISYMB  = ILSTSYM(LISTB,IDLSTB)
      ISYMC  = ILSTSYM(LISTC,IDLSTC)
      ISYMABC = MULD2H(MULD2H(ISYMA,ISYMB),ISYMC)

      KTHETA1 = ICTRAN(4,1)
      KTHETA2 = KTHETA1 + NT1AM(ISYMABC)

      IF (NSYM.EQ.1 .AND. LOCDBG) THEN
        KT1AMPC = KTHETA2 + NT2AM(ISYMABC)
        KT2AMPC = KT1AMPC + NT1AM(ISYMC)
        KRESLT1 = KT2AMPC + NT2AM(ISYMC)
        KRESLT2 = KRESLT1 + NT1AM(ISYMABC)
        KEND1   = KRESLT2 + NT2AM(ISYMABC)
        LEND1   = LWORK - KEND1

        IF (LEND1 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_BATST.')
        END IF

        IOPT = 3
        Call CC_RDRSP(LISTC,IDLSTC,ISYMC,IOPT,MODEL,
     &                WORK(KT1AMPC),WORK(KT2AMPC))

        WRITE (LUPRI,*) 'CC_BATST: C vector:'
        Call CC_PRP(WORK(KT1AMPC),WORK(KT2AMPC),ISYMC,1,1)

        ! zero singles or doubles C vector:
C       CALL DZERO(WORK(KT1AMPC),NT1AM(ISYMC))
C       CALL DZERO(WORK(KT2AMPC),NT2AM(ISYMC))
        CALL DZERO(WORK(KRESLT1),NT1AM(ISYMABC)+NT2AM(ISYMABC))

        IF (NEW_CODE) THEN
           CALL CC_FDBAMAT(WORK(KRESLT1),WORK(KRESLT2),
     >                     LISTB,IDLSTB,LISTC,IDLSTC,
     >                     LABELA,IRELAXA,WORK(KEND1),LEND1)
        ELSE
           IPRINT  = 6
           CALL CC_FDBA(NT1AM(ISYMABC),NT2AM(ISYMABC),
     >                  LISTA,IDLSTA,LISTB,IDLSTB,  
     >                  WORK(KT1AMPC), WORK(KRESLT1),
     >                  WORK(KEND1), LEND1)
        ENDIF

        IPRINT  = 0

        IF (.TRUE.) THEN
          WRITE (LUPRI,*)
          WRITE (LUPRI,*)
          WRITE (LUPRI,*) 'FINITE DIFFERENCE TEST FOR B{O} MATRIX:'
          WRITE (LUPRI,*) '---------------------------------------'
          WRITE (LUPRI,*) 
     &          'LISTA, IDLSTA, ISYMA:',LISTA(1:2),IDLSTA,ISYMA
          WRITE (LUPRI,*) 
     &          'LISTB, IDLSTB, ISYMB:',LISTB(1:2),IDLSTB,ISYMB
          WRITE (LUPRI,*) 
     &          'LISTC, IDLSTC, ISYMC:',LISTC(1:2),IDLSTC,ISYMC
          WRITE (LUPRI,*) 'ISYMABC:',ISYMABC
          WRITE (LUPRI,*)
          WRITE (LUPRI,*) 'finite difference Theta vector:'
          Call CC_PRP(WORK(KRESLT1),WORK(KRESLT2),ISYMABC,1,1)
          WRITE (LUPRI,*) 'analytical Theta vector:'
          Call CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYMABC,1,1)
        END IF

        Call DAXPY(NT1AM(ISYMABC),-1.0d0,WORK(KTHETA1),1,
     &                                  WORK(KRESLT1),1)
        IF (.NOT.CCS) THEN
          Call DAXPY(NT2AM(ISYMABC),-1.0d0,WORK(KTHETA2),1,
     &                                    WORK(KRESLT2),1)
        ELSE
          Call DZERO(WORK(KRESLT2),NT2AM(ISYMABC))
        END IF

        WRITE (LUPRI,*)
        WRITE (LUPRI,*) 'FINITE DIFFERENCE TEST FOR B{O} MATRIX:'
        WRITE (LUPRI,*) '---------------------------------------'
        WRITE (LUPRI,*) 'Norm of difference between analytical THETA '
     >           // 'vector and the numerical result:'
        WRITE (LUPRI,*) 'singles excitation part:',
     >   DSQRT(DDOT(NT1AM(ISYMA),WORK(KRESLT1),1,WORK(KRESLT1),1))
        WRITE (LUPRI,*) 'double excitation part: ',
     >   DSQRT(DDOT(NT2AM(ISYMA),WORK(KRESLT2),1,WORK(KRESLT2),1))

        WRITE (LUPRI,*) 'difference vector:'
        Call CC_PRP(WORK(KRESLT1),WORK(KRESLT2),ISYMABC,1,1)

        CALL FLSHFO(LUPRI)


      ELSE IF (NSYM.NE.1 .AND. LOCDBG) THEN
       WRITE (LUPRI,*) 'CC_BATST> can not calculate finite diff. '//
     &                 'B{O} matrix'
       WRITE (LUPRI,*) 'CC_BATST> with symmetry.'
      END IF

      RETURN
      END 
*=====================================================================*
