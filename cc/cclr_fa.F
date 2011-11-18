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
c /* deck cclr_fa */
*=====================================================================*
       SUBROUTINE CCLR_FA ( LABELA, ISYMTA,  ! inp: label/symmetry A
     &                      LISTB,  ITAMPB,  ! inp: B resp. amplit.
     &                      LISTC,  IZETVC,  ! inp: C resp. zeta vec.
     &                      XINT,            ! inp: integrals of operator
     &                      WORK,   LWORK   )! work space
*---------------------------------------------------------------------*
*
*    Purpose: transformation of a response vector with a F matrix
*             where the hamiltonian has been substituted by a 
*             perturbation operator
*
*             F^C{A} * t^B = <lambda^C|[[A,t^B],tau_nu]|CC>
*
*             JK+OC. CCMM: Allow for input of integrals if
*             LABELA .eq. 'GIVE INT' 
*
*    symmetries/variables:
*              
*           ISYRES : result vector GAMMA1, GAMMA2
*           ISYCTR : lagrangian multipliers (zeta vector) CTR1, CTR2
*           ISYMTA : A perturbation
*           ISYMTB : B response amplitudes 
*
*    Note: the single and double excitation parts of the result GAMMA2 
*          are returned at the beginning of the work space in
*          WORK(1)... WORK(NT1AM(ISYRES))
*          WORK(NT1AM(ISYRES)+1)... WORK(NT1AM(ISYRES)+NT2AM(ISYRES))
*          (double excitation part will be stored in packed form)
*
*     Written by Christof Haettig, October 1996.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <iratdef.h>
#include <cbieri.h>
#include <mxcent.h>
#include <eribuf.h>
#include <maxorb.h>
#include <distcl.h>
#include <ccorb.h>
#include <ccisao.h>
#include <ccsdsym.h>
#include <ccsdinp.h>

* local parameters:
      CHARACTER MSGDBG*(17)
      PARAMETER (MSGDBG='[debug] CCLR_FA> ')
#if defined (SYS_CRAY)
      REAL ONE, TWO, THREE 
#else
      DOUBLE PRECISION ONE, TWO, THREE 
#endif
      PARAMETER (ONE = 1.0d0, TWO = 2.0d0, THREE = 3.0d0)
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER KDUM
      PARAMETER (KDUM = +99 999 999) ! dummy address

      CHARACTER*8 LABELA
      CHARACTER*10 MODEL
      CHARACTER LISTB*(*), LISTC*(*)
      INTEGER ISYRES, ISYCTR, ISYMTA, ISYMTB, LWORK
      INTEGER ITAMPB, IZETVC

      INTEGER ISYTATB, ISYMJ, ISYMB, ISYMXY, ISYMI, ISYMA, IRREP, ISYM
      INTEGER KBTAOO, KBTAVV, KPERTA, KT1AMPB, KT2AMPB, KCTR1, KCTR2
      INTEGER KCTMO, KT1AMP0, KLAMDP0, KLAMDH0, KOFF1, KOFF2, KSCR
      INTEGER KEND1, KEND2, KEND3, KEND1A, KXBMAT, KYBMAT, IERR
      INTEGER LEND1, LEND2, LEND3, LEND1A, KGAMMA1, KGAMMA2, KEND0
      INTEGER IOPT, MAXJ, NIJ, NJI, NAB, NBA, KEMAT1, KEMAT2

#if defined (SYS_CRAY)
      REAL FREQB
      REAL SWAP
      REAL WORK(LWORK), XINT(*)
#else
      DOUBLE PRECISION FREQB
      DOUBLE PRECISION SWAP
      DOUBLE PRECISION WORK(LWORK), XINT(*)
#endif

      INTEGER ILSTSYM
  
*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (CCSDT) THEN
        WRITE (LUPRI,'(/1x,a)') 'F{A} matrix transformations not '
     &          //'implemented for triples yet...'
        CALL QUIT('Triples not implemented for G '//
     &       'matrix transformations')
      END IF

      IF ( .not. (CCS .or. CC2 .or. CCSD) ) THEN
        WRITE (LUPRI,'(/1x,a)') 'CCLR_FA called for a Coupled Cluster '
     &          //'method not implemented in CCLR_FA...'
        CALL QUIT('Unknown CC method in CCLR_FA.')
      END IF

*---------------------------------------------------------------------*
* set & check symmetries:
*---------------------------------------------------------------------*
      ISYMTB  = ILSTSYM(LISTB,ITAMPB)   ! B
      ISYCTR  = ILSTSYM(LISTC,IZETVC)   ! C
      ISYTATB = MULD2H(ISYMTA,ISYMTB)   ! A x B
      ISYRES  = MULD2H(ISYCTR,ISYTATB)  ! A x B x C

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'LISTB,ITAMPB,ISYMTB:',LISTB,ITAMPB,ISYMTB
        WRITE (LUPRI,*) 'LISTC,IZETVC,ISYCTR:',LISTC,IZETVC,ISYCTR
      END IF


      IF (ISYMOP .NE. 1) THEN
        WRITE (LUPRI,'(/1x,a)') 'non-total-symmetric MO integrals?!... '
     &          //'CCLR_G has never been debugged for that!...'
      END IF

      IF (MULD2H(ISYCTR,ISYTATB) .NE. ISYRES) THEN
        CALL QUIT('Symmetry mismatch in CCLR_FA.')
      END IF

*---------------------------------------------------------------------*
* flush print unit
*---------------------------------------------------------------------*
      Call FLSHFO(LUPRI)

      IF (LOCDBG) THEN
        WRITE (LUPRI,'(/1x,a,i15)') 'work space in CCLR_FA:',LWORK
      END IF
*---------------------------------------------------------------------*
* initialize pointer for work space and allocate memory for
*  1) single excitation part of the result vector 
*  2) one-index transformed perturbation integrals A^B (occ/occ block)
*  3) one-index transformed perturbation integrals A^B (vir/vir block)
*  4) perturbation integrals A
*  5) singles part of response amplitudes T1^B
*  6) singles part of zeroth order lagrangian multipliers
*---------------------------------------------------------------------*
      KGAMMA1 = 1
      KBTAOO  = KGAMMA1 + NT1AM(ISYRES)
      KBTAVV  = KBTAOO  + NMATIJ(ISYTATB)
      KPERTA  = KBTAVV  + NMATAB(ISYTATB)
      KT1AMPB = KPERTA  + NT1AM(ISYMTA)
      KCTR1   = KT1AMPB + NT1AM(ISYMTB)
      KEND1   = KCTR1   + NT1AM(ISYCTR)
      LEND1   = LWORK - KEND1

      IF (LEND1 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CCLR_FA.')
      END IF

*---------------------------------------------------------------------*
* initialize single excitation part of result vector GAMMA1:
*---------------------------------------------------------------------*
      Call DZERO (WORK(KGAMMA1), NT1AM(ISYRES))

*---------------------------------------------------------------------*
* for CCS and zeroth-order zeta vector all contributions vanish:
*---------------------------------------------------------------------*
      IF (CCS .AND. LISTC(1:2).EQ.'L0') RETURN
      
*---------------------------------------------------------------------*
* read singles parts for B response amplitudes and zeta vector:
*---------------------------------------------------------------------*
      IOPT = 1
      CALL CC_RDRSP(LISTB,ITAMPB,ISYMTB,IOPT,MODEL,
     &                  WORK(KT1AMPB),WORK(KDUM)  )


      IOPT = 1
      Call CC_RDRSP(LISTC,IZETVC,ISYCTR,IOPT,MODEL,
     &                  WORK(KCTR1),WORK(KDUM))

      IF (LOCDBG) THEN
        CAll AROUND('response T amplitudes B:')
        WRITE (LUPRI,*) 'LIST/INDEX:',LISTB,ITAMPB
        WRITE (LUPRI,*) 'Symmetry:      ',ISYMTB
        CAll CC_PRP(WORK(KT1AMPB),WORK(KDUM),ISYMTB,1,0)
        CALL AROUND('CC lagrange multipliers')
        CALL CC_PRP(WORK(KCTR1), WORK(KDUM),  ISYCTR, 1, 0)
      END IF

*---------------------------------------------------------------------*
* read & resort one-electron integrals for operator A:
*---------------------------------------------------------------------*
      KCTMO   = KEND1   
      KT1AMP0 = KCTMO   + N2BST(ISYMTA)
      KLAMDP0 = KT1AMP0 + NT1AM(ISYMOP)
      KLAMDH0 = KLAMDP0 + NLAMDT
      KEND1A  = KLAMDH0 + NLAMDT
      LEND1A  = LWORK - KEND1A

      IF (LEND1A .LT. 0) THEN
        CALL QUIT('Insufficient work space in CCLR_FA.')
      END IF
C
C     JK+OC, CCSLV
      IF (LABELA.EQ.'GIVE INT') THEN 
        CALL DCOPY(N2BST(ISYMTA),XINT(1),1,WORK(KCTMO),1)
      ELSE
* read the AO integrals:
       CALL CCPRPAO(LABELA,.TRUE.,WORK(KCTMO),IRREP,ISYM,IERR,
     &             WORK(KEND1A),LEND1A)
       IF (IERR.NE.0 .OR. IRREP.NE.ISYMTA) THEN
        CALL QUIT('CCLR_FA: error while reading operator '//LABELA)
       END IF
      END IF
C
* get MO coefficients:
      CALL DZERO(WORK(KT1AMP0),NT1AMX)
      CALL LAMMAT(WORK(KLAMDP0),WORK(KLAMDH0),WORK(KT1AMP0),
     &            WORK(KEND1A),LEND1A)

* transform one-electron integrals in place:
      CALL CC_FCKMO(WORK(KCTMO),WORK(KLAMDP0),WORK(KLAMDH0),
     &              WORK(KEND1A),LEND1A,ISYMTA,1,1)

* resort occupied/virtual block to T1 like storage:
      CALL DZERO(WORK(KPERTA),NT1AM(ISYMTA))
      DO ISYMJ = 1, NSYM
        ISYMB = MULD2H(ISYMJ,ISYMTA)

        DO J = 1, NRHF(ISYMJ)
        DO B = 1, NVIR(ISYMB)
          KOFF1 = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1) + B
          KOFF2 = IFCVIR(ISYMJ,ISYMB) + NORB(ISYMJ)*(B-1) + J

          WORK(KPERTA-1+KOFF1) = WORK(KCTMO-1+KOFF2)
        END DO
        END DO
      END DO

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG,' A integrals in MO basis:'
        WRITE (LUPRI,*) MSGDBG,' label, symmetry:',LABELA,ISYMTA
        Call CC_PRP(WORK(KPERTA),WORK(KDUM),ISYMTA,1,0)
      END IF

*---------------------------------------------------------------------*
* calculate A perturbation integrals one-index transformed with
* the B response amplitudes T1^B:
*---------------------------------------------------------------------*
* occ/occ block:
      Call CCG_1ITROO(WORK(KBTAOO), ISYTATB,
     &                WORK(KPERTA), ISYMTA,
     &                WORK(KT1AMPB),ISYMTB  )

* vir/vir block:
      Call CCG_1ITRVV(WORK(KBTAVV), ISYTATB,
     &                WORK(KPERTA), ISYMTA,
     &                WORK(KT1AMPB),ISYMTB  )

*=====================================================================*
*   CCS part:  < Zeta_1 | [tA^B, tau_1] | HF>
*=====================================================================*
* do one-index transformation with Zeta vector:
      IOPT  = 2
      Call CCG_1ITRVO(WORK(KGAMMA1),ISYRES,WORK(KBTAOO),WORK(KBTAVV),
     &                ISYTATB,WORK(KCTR1),ISYCTR,IOPT          )

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG,'one-index trans. A (occ/occ block):'
        WRITE (LUPRI,'(5f12.6)') (WORK(KBTAOO-1+I),I=1,NMATIJ(ISYTATB))
        WRITE (LUPRI,*) MSGDBG,'one-index trans. A (vir/vir block):'
        WRITE (LUPRI,'(2f12.6)') (WORK(KBTAVV-1+I),I=1,NMATAB(ISYTATB))
        WRITE (LUPRI,*) MSGDBG, 
     *            'contrib. of one-index trans. A to GAMMA:'
        Call CC_PRP(WORK(KGAMMA1),WORK(KDUM),ISYRES,1,0)
      END IF

*---------------------------------------------------------------------*
* end of CCS part
*---------------------------------------------------------------------*
      
      If (CCS) Return

*=====================================================================*
* CC2/CCSD part for the singles: <Zeta_2| [[A,T2B], tau_1] |HF>
*=====================================================================*

*---------------------------------------------------------------------*
* memory allocation:
* 1) double excitation part of response amplitudes T2B (packed)
* 2) double excitation part of zeta vector (squared)
* 3) double excitation part of zeta vector (packed)
* N.B. we account here for the fact, that the packed double excitation 
* part of the result vector will be returned at the beginning of the
* work space, so we make sure, that there is enough space before
* the zeta vector to store there later on GAMMA2
*---------------------------------------------------------------------*
       KT2AMPB = KEND1
       KCTR2   = KT2AMPB + MAX( NT2AM(ISYMTB), NT2AM(ISYRES) )
       KEND2   = KCTR2 + NT2SQ(ISYCTR)
       LEND2   = LWORK - KEND2

       IF (LEND2 .LT. NT2AM(ISYCTR) ) THEN
         CALL QUIT('Insufficient work space in CCLR_FA.')
       END IF

*---------------------------------------------------------------------*
* read response amplitudes T2B and scale the diagonal:
*---------------------------------------------------------------------*
      IOPT = 2
      CALL CC_RDRSP(LISTB,ITAMPB,ISYMTB,IOPT,MODEL,
     &                  WORK(KDUM),WORK(KT2AMPB)  )

      CAll CCLR_DIASCL(WORK(KT2AMPB),TWO,ISYMTB)

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG, 'B response amplitudes:'
        Call CC_PRP(WORK(KT1AMPB),WORK(KT2AMPB),ISYMTB,1,1)
      END IF

*---------------------------------------------------------------------*
* read packed lagrangian multipliers and square them up:
*---------------------------------------------------------------------*
      IOPT = 2  
      Call CC_RDRSP(LISTC,IZETVC,ISYCTR,IOPT,MODEL,
     &                  WORK(KDUM),WORK(KEND2))

      CALL CC_T2SQ (WORK(KEND2), WORK(KCTR2), ISYCTR)

*---------------------------------------------------------------------*
* calculate X^B and Y^B intermediates:
*---------------------------------------------------------------------*
      ISYMXY  = MULD2H(ISYCTR,ISYMTB)

      KXBMAT  = KEND2
      KYBMAT  = KXBMAT  + NMATIJ(ISYMXY)
      KSCR    = KYBMAT  + NMATAB(ISYMXY)
      KEND3   = KSCR    + NT1AM(ISYRES)
      LEND3   = LWORK - KEND3

      If (LEND3 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CCLR_FA.')
      END IF

* calculate X^C & Y^C intermediate:
      Call CC_XI(WORK(KXBMAT),WORK(KCTR2), ISYCTR,
     &           WORK(KT2AMPB),ISYMTB,WORK(KEND3),LEND3)

      Call CC_YI(WORK(KYBMAT),WORK(KCTR2), ISYCTR,
     &           WORK(KT2AMPB),ISYMTB,WORK(KEND3),LEND3)

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG,'response X intermediate:'
        WRITE (LUPRI,'(5f12.6)') (WORK(KXBMAT-1+I),I=1,NMATIJ(ISYMXY))
        WRITE (LUPRI,*) MSGDBG,'response Y intermediate:'
        WRITE (LUPRI,'(2f12.6)') (WORK(KYBMAT-1+I),I=1,NMATAB(ISYMXY))
      END IF

* calculate XY^B:  XY^B_ij = X^B_ji,  XY^B_bd = -Y^B_bd
* 1.) transpose X^B intermediate
      DO ISYMI = 1, NSYM
        ISYMJ = MULD2H(ISYMI,ISYMXY)
        IF (ISYMJ .LE. ISYMI) THEN
          DO I = 1, NRHF(ISYMI)
            MAXJ =  NRHF(ISYMJ)
            IF (ISYMJ .EQ. ISYMI) MAXJ = I-1
          DO J = 1, MAXJ
            NIJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J-1) + I
            NJI = IMATIJ(ISYMJ,ISYMI) + NRHF(ISYMJ)*(I-1) + J
            SWAP = WORK(KXBMAT-1+NIJ)
            WORK(KXBMAT-1+NIJ) = WORK(KXBMAT-1+NJI)
            WORK(KXBMAT-1+NJI) = SWAP
          END DO
          END DO
        END IF
      END DO

* 2.) multiply Y^B intermediate with -1:
      Call DSCAL(NMATAB(ISYMXY), -ONE, WORK(KYBMAT), 1)


* do one-index transformation of XY^B with A integrals:
      IOPT  = 2
      Call CCG_1ITRVO(WORK(KSCR),ISYRES,
     &                WORK(KXBMAT),WORK(KYBMAT),ISYMXY,
     &                WORK(KPERTA),ISYMTA,      IOPT    )

* add contribution to GAMMA1:
      Call DAXPY (NT1AM(ISYRES), ONE, WORK(KSCR),1, WORK(KGAMMA1), 1)

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG,'A integrals:'
        WRITE (LUPRI,'(5f12.6)') (WORK(KPERTA-1+I),I=1,NT1AM(ISYMTA))
        WRITE (LUPRI,*) MSGDBG, 'CC2/CCSD contribution to singles part:'
        Call CC_PRP(WORK(KSCR),WORK(KDUM),ISYRES,1,0)
        WRITE (LUPRI,*) MSGDBG, 'GAMMA1 now:'
        Call CC_PRP(WORK(KGAMMA1),WORK(KDUM),ISYRES,1,0)
      END IF

*---------------------------------------------------------------------*

*=====================================================================*
* CC2/CCSD part for the doubles: <Zeta_2| [[A,T1B], tau_2] |HF>
*=====================================================================*

*---------------------------------------------------------------------*
* reorganize work space, so that the result vector GAMMA2 can be
* stored at the early beginning of the work space
*---------------------------------------------------------------------*
      KGAMMA1 = 1
      KGAMMA2 = KGAMMA1 + NT1AM(ISYRES)
      KEND0   = KGAMMA2 + NT2AM(ISYRES)

      IF (KEND0 .GT. KCTR2) THEN
        CALL QUIT('memory organization mixed up in CCLR_FA.')
      END IF

      KEMAT1 = KCTR2  + NT2SQ(ISYCTR) 
      KEMAT2 = KEMAT1 + NMATAB(ISYTATB)
      KEND3  = KEMAT2 + NMATIJ(ISYTATB)
      LEND3  = LWORK - KEND3
 
      IF ( LEND3 .LT. 0 ) THEN
        CALL QUIT('Insufficient work space in CCLR_FA.')
      END IF

*---------------------------------------------------------------------*
* transpose tA^B(a b) --> EMAT1(b a)
*---------------------------------------------------------------------*
      DO ISYMA = 1, NSYM
        ISYMB = MULD2H(ISYMA,ISYTATB)
        DO A = 1, NVIR(ISYMA)
        DO B = 1, NVIR(ISYMB)
          NAB = IMATAB(ISYMA,ISYMB) + NVIR(ISYMA)*(B-1) + A
          NBA = IMATAB(ISYMB,ISYMA) + NVIR(ISYMB)*(A-1) + B
         
          WORK(KEMAT1 - 1 + NBA) = WORK(KBTAVV - 1 + NAB)
        END DO
        END DO
      END DO


*---------------------------------------------------------------------*
* transpose tA^B(i j) --> EMAT2(j i)
*---------------------------------------------------------------------*
      DO ISYMI = 1, NSYM
        ISYMJ = MULD2H(ISYMI,ISYTATB)
        DO I = 1, NRHF(ISYMI)
        DO J = 1, NRHF(ISYMJ)
          NIJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J-1) + I
          NJI = IMATIJ(ISYMJ,ISYMI) + NRHF(ISYMJ)*(I-1) + J
         
          WORK(KEMAT2 - 1 + NJI) = WORK(KBTAOO - 1 + NIJ)
        END DO
        END DO
      END DO

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG,'one-index trans. A (occ/occ block):'
        WRITE (LUPRI,'(5f12.6)') (WORK(KBTAOO-1+I),I=1,NMATIJ(ISYTATB))
        WRITE (LUPRI,*) MSGDBG,'one-index trans. A (vir/vir block):'
        WRITE (LUPRI,'(2f12.6)') (WORK(KBTAVV-1+I),I=1,NMATAB(ISYTATB))
        WRITE (LUPRI,*) MSGDBG,'EMAT2:'
        WRITE (LUPRI,'(5f12.6)') (WORK(KEMAT2-1+I),I=1,NMATIJ(ISYTATB))
        WRITE (LUPRI,*) MSGDBG,'EMAT1:'
        WRITE (LUPRI,'(2f12.6)') (WORK(KEMAT1-1+I),I=1,NMATAB(ISYTATB))
      END IF

*---------------------------------------------------------------------*
* combine EMAT1/EMAT2 with lagrangian multipliers:
* (note: this overwrites the intermedites stored at the beginning
*        of the work space...)
*---------------------------------------------------------------------*
* initialize GAMMA2:
      CALL DZERO(WORK(KGAMMA2),NT2AM(ISYRES))

* do the caculation:
      Call CCRHS_E(WORK(KGAMMA2),WORK(KCTR2),WORK(KEMAT1), 
     &             WORK(KEMAT2), WORK(KEND3), LEND3, ISYCTR, ISYTATB)

*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG, 'GAMMA:'
        Call CC_PRP(WORK(KGAMMA1),WORK(KGAMMA2),ISYRES,1,1)
      END IF
C
      RETURN
      END
*=====================================================================*
*                  END OF SUBROUTINE CCLR_FA                          *
*=====================================================================*
