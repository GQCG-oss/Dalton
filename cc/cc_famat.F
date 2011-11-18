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
c /* deck cc_famat */
*=====================================================================*
       SUBROUTINE CC_FAMAT( LABELA, ISYMTA,  ! inp: label/symmetry A
     &                      LISTB,  ITAMPB,  ! inp: B resp. amplit.
     &                      LISTC,  IZETVC,  ! inp: C resp. zeta vec.
     &                      IOPTRES,IFATRAN, ! output option
     &                      IFADOTS,FACON,   ! indeces/dotproducts
     &                      FILFA,  ITRAN,   ! list,index for dotprods
     &                      NFATRAN,MXVEC,   ! dimensions for dotprods
     &                      NODDY_CCSDT,     ! flag for noddy triples 
     &                      WORK,   LWORK   )! work space
*---------------------------------------------------------------------*
*
*    Purpose: transformation of a response vector with a F matrix
*             where the hamiltonian has been substituted by a 
*             perturbation operator
*
*             F^C{A} * t^B = <lambda^C|[[A,t^B],tau_nu]|CC>
*
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
*     CC3 noddy version, April 2002, Christof Haettig
*     Added CCR12, June 2005, Christian Neiss
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
#include <cclists.h>
#include <r12int.h>
#include <ccr12int.h>

* local parameters:
      CHARACTER MSGDBG*(18)
      PARAMETER (MSGDBG='[debug] CC_FAMAT> ')
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
      CHARACTER*10 MODEL, MODELW
      CHARACTER LISTB*(*), LISTC*(*), FILFA*(*)
      LOGICAL NODDY_CCSDT 
      INTEGER ISYRES, ISYCTR, ISYMTA, ISYMTB, LWORK
      INTEGER ITAMPB, IZETVC

      INTEGER ISYTATB, ISYMJ, ISYMB, ISYMXY, ISYMI, ISYMA, IRREP, ISYM
      INTEGER KBTAOO, KBTAVV, KPERTA, KT1AMPB, KT2AMPB, KCTR1, KCTR2
      INTEGER KCTMO, KT1AMP0, KLAMDP0, KLAMDH0, KOFF1, KOFF2, KSCR
      INTEGER KEND1, KEND2, KEND3, KEND1A, KXBMAT, KYBMAT, IERR
      INTEGER LEND1, LEND2, LEND3, LEND1A, KGAMMA1, KGAMMA2, KEND0
      INTEGER IOPT, MAXJ, NIJ, NJI, NAB, NBA, KEMAT1, KEMAT2, IVEC,
     &        LEND0, KGAMMA1EFF, KGAMMA2EFF, IOPTW, IOPTWE, IFILE

      INTEGER NFATRAN, MXVEC, IOPTRES, ITRAN
      INTEGER IFATRAN(MXDIM_FATRAN,NFATRAN), IFADOTS(MXVEC,NFATRAN)
      INTEGER KGAMMAR12, KCTRR12, KT12AMB, KXINTTRI, KXINTSQ, IDUM
      INTEGER LUNIT, IAN, KLAMDPB, KLAMDHB, LENMOD, IOPTWR12

#if defined (SYS_CRAY)
      REAL FREQB, FREQA
      REAL SWAP, DUMMY
      REAL WORK(LWORK), FACON(MXVEC,NFATRAN)
#else
      DOUBLE PRECISION FREQB, FREQA
      DOUBLE PRECISION SWAP, DUMMY
      DOUBLE PRECISION WORK(LWORK), FACON(MXVEC,NFATRAN)
#endif

      INTEGER ILSTSYM
      CHARACTER*(3) APROXR12

  
      CALL QENTER('CC_FAMAT')
*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF ( .not. (CCS .or. CC2 .or. CCSD .or. CC3) ) THEN
        WRITE (LUPRI,'(/1x,a)') 'CC_FAMAT called for a Coupled Cluster '
     &          //'method not implemented in CC_FAMAT...'
        CALL QUIT('Unknown CC method in CC_FAMAT.')
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
        CALL QUIT('Symmetry mismatch in CC_FAMAT.')
      END IF

*---------------------------------------------------------------------*
* flush print unit
*---------------------------------------------------------------------*
      Call FLSHFO(LUPRI)

      IF (LOCDBG) THEN
        WRITE (LUPRI,'(/1x,a,i15)') 'work space in CC_FAMAT:',LWORK
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
      KEND0   = KGAMMA1 + NT1AM(ISYRES)
      LEND0   = LWORK   - KEND0

      KBTAOO  = KEND0
      KBTAVV  = KBTAOO  + NMATIJ(ISYTATB)
      KPERTA  = KBTAVV  + NMATAB(ISYTATB)
      KT1AMPB = KPERTA  + NT1AM(ISYMTA)
      KCTR1   = KT1AMPB + NT1AM(ISYMTB)
      KEND1   = KCTR1   + NT1AM(ISYCTR)
      LEND1   = LWORK - KEND1

      IF (LEND1 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CC_FAMAT.')
      END IF

*---------------------------------------------------------------------*
* initialize single excitation part of result vector GAMMA1:
*---------------------------------------------------------------------*
      Call DZERO (WORK(KGAMMA1), NT1AM(ISYRES))

*---------------------------------------------------------------------*
* for CCS and zeroth-order zeta vector all contributions vanish:
*---------------------------------------------------------------------*
      IF (CCS .AND. LISTC(1:2).EQ.'L0') GOTO 9999
      
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
        CALL QUIT('Insufficient work space in CC_FAMAT.')
      END IF

* read the AO integrals:
      CALL CCPRPAO(LABELA,.TRUE.,WORK(KCTMO),IRREP,ISYM,IERR,
     &             WORK(KEND1A),LEND1A)
      IF (IERR.NE.0 .OR. IRREP.NE.ISYMTA) THEN
         CALL QUIT('CC_FAMAT: error while reading operator '//LABELA)
      END IF


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
      
      If (CCS) GOTO 9999

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
         CALL QUIT('Insufficient work space in CC_FAMAT.')
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
        CALL QUIT('Insufficient work space in CC_FAMAT.')
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
      LEND0   = LWORK   - KEND0

      IF (KEND0 .GT. KCTR2) THEN
        CALL QUIT('memory organization mixed up in CC_FAMAT.')
      END IF

      KEMAT1 = KCTR2  + NT2SQ(ISYCTR) 
      KEMAT2 = KEMAT1 + NMATAB(ISYTATB)
      KEND3  = KEMAT2 + NMATIJ(ISYTATB)
      LEND3  = LWORK - KEND3
 
      IF ( LEND3 .LT. 0 ) THEN
        CALL QUIT('Insufficient work space in CC_FAMAT.')
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

*=====================================================================*
* CC-R12 part:
*=====================================================================*
      IF (CCR12) THEN 
        IF (.NOT.IANR12.EQ.1) THEN
          CALL QUIT('Sorry, only CC-R12 Ansatz 1 implemented yet!')
        ENDIF
* allocate memory:
        KGAMMAR12 = KEND0
        KEND0     = KGAMMAR12 + NTR12AM(ISYRES)
        LEND0     = LWORK - KEND0
C
        KCTRR12 = KEND0
        KT12AMB = KCTRR12 + NTR12SQ(ISYCTR)
        KXINTTRI= KT12AMB + MAX(NTR12SQ(ISYMTB),NTR12SQ(ISYRES))
        KXINTSQ = KXINTTRI+ NR12R12P(1)
        KSCR    = KXINTSQ + NR12R12SQ(1)
        KCTMO   = KSCR    + NTR12SQ(ISYRES)
        KT1AMP0 = KCTMO   + N2BST(ISYMTA)
        KLAMDP0 = KT1AMP0 + NT1AMX
        KLAMDH0 = KLAMDP0 + NLAMDT
        KT1AMPB = KLAMDH0 + NLAMDT 
        KLAMDPB = KT1AMPB + NT1AM(ISYMTB)
        KLAMDHB = KLAMDPB + NGLMDT(ISYMTB)
        KEND1   = KLAMDHB + NGLMDT(ISYMTB)
        LEND1   = LWORK - KEND1
        IF (LEND1.LT.0) THEN
          CALL QUIT('Insufficient work space for R12 in CC_FAMAT')
        END IF
        
*---------------------------------------------------------------------*
* read in everything needed: 
*---------------------------------------------------------------------*
* read the AO integrals:
        CALL CCPRPAO(LABELA,.TRUE.,WORK(KCTMO),IRREP,ISYM,IERR,
     &               WORK(KEND1),LEND1)
        IF (IERR.NE.0 .OR. IRREP.NE.ISYMTA) THEN
          CALL QUIT('CC_FAMAT: error while reading operator '//LABELA)
        END IF

* get MO coefficients:
        CALL DZERO(WORK(KT1AMP0),NT1AMX)
        CALL LAMMAT(WORK(KLAMDP0),WORK(KLAMDH0),WORK(KT1AMP0),
     &              WORK(KEND1),LEND1)        
* transform MO matrices with C1:
        IOPT = 1
        Call CC_RDRSP(LISTB,ITAMPB,ISYMTB,IOPT,MODEL,WORK(KT1AMPB),
     &                DUMMY)
        CALL CCLR_LAMTRA(WORK(KLAMDP0),WORK(KLAMDPB),WORK(KLAMDH0),
     &                   WORK(KLAMDHB),WORK(KT1AMPB),ISYMTB)

* read R12 Lagrangian multipliers:
        CALL CC_R12GETCT(WORK(KCTRR12),ISYCTR,2,2.0D0*BRASCL,.FALSE.,
     &                   'N',IDUM,IDUM,IDUM,LISTC,IZETVC,
     &                   WORK(KEND1),LEND1) 

* read R12 response amplitudes:
        CALL CC_R12GETCT(WORK(KT12AMB),ISYMTB,2,KETSCL,.FALSE.,'N',
     &                   IDUM,IDUM,IDUM,LISTB,ITAMPB,WORK(KEND1),LEND1)

* read R12 overlap matrix:
        LUNIT = -1
        CALL GPOPEN(LUNIT,FCCR12X,'OLD',' ','UNFORMATTED',
     &              IDUM,.FALSE.)
        REWIND(LUNIT)
 8888   READ(LUNIT) IAN
        READ(LUNIT) (WORK(KXINTTRI-1+I), I=1, NR12R12P(1))
        IF (IAN.NE.IANR12) GOTO 8888
        CALL GPCLOSE(LUNIT,'KEEP')
        IOPT = 2
        CALL CCR12UNPCK2(WORK(KXINTTRI),1,WORK(KXINTSQ),'N',IOPT) 
        
*---------------------------------------------------------------------*
* Calculate singles contribution:
*---------------------------------------------------------------------*
       CALL CC_R12ETAA(WORK(KGAMMA1),ISYRES,WORK(KCTRR12),ISYCTR,
     &                  WORK(KT12AMB),ISYMTB,WORK(KXINTSQ),WORK(KCTMO),
     &                  ISYMTA,WORK(KLAMDP0),WORK(KLAMDH0),.FALSE.,
     &                  WORK(KEND1),LEND1)

*---------------------------------------------------------------------*
* Calculate R12 doubles contribution:
*---------------------------------------------------------------------*
        CALL DZERO(WORK(KT12AMB),NTR12SQ(ISYRES))
        CALL CC_R12XI2A(WORK(KSCR),ISYRES,WORK(KCTRR12),ISYCTR,
     &                  WORK(KCTMO),ISYMTA,WORK(KLAMDP0),WORK(KLAMDHB),
     &                  ISYMTB,'T',WORK(KEND1),LEND1)
        CALL CC_R12XI2B(WORK(KT12AMB),'N',WORK(KXINTSQ),1,'N',
     &                  WORK(KSCR),ISYRES,'N',-1.0D0)
        IOPT = 1
        CALL CCR12PCK2(WORK(KGAMMAR12),ISYRES,.FALSE.,WORK(KT12AMB),
     &                 'N',IOPT)
        CALL CCLR_DIASCLR12(WORK(KGAMMAR12),0.5D0*KETSCL,ISYRES)

      END IF 

*=====================================================================*
* CC3 part:
*=====================================================================*
      IF (CCSDT) THEN
        IF (IOPTRES.LT.5) THEN
          KGAMMA1EFF = KEND0
          KGAMMA2EFF = KGAMMA1EFF + NT1AM(ISYRES)
          KEND0      = KGAMMA2EFF + NT2AM(ISYRES)
        END IF
        LEND0 = LWORK - KEND0

        IF ( LEND0 .LT. 0 ) THEN
          CALL QUIT('Insufficient work space in CC_FAMAT.')
        END IF

        IF (IOPTRES.GE.1 .AND. IOPTRES.LE.4) THEN
          CALL QUIT('CC_FAMAT needs to be fixed for this case')
          ! in this case we need to find the frequency associated
          ! with the A perturbation, such that we can construct
          ! the correct effective rhs vector...
          ! FREQA = 
        ELSE IF (IOPTRES.EQ.5) THEN
          FREQA = 0.0D0
        ELSE
          CALL QUIT('Illegal value for IOPTRES in CC_FAMAT.')
        END IF

        IF (IOPTRES.NE.5 .OR. NODDY_CCSDT) THEN

          CALL CCSDT_FAMAT_NODDY(LISTC,IZETVC,LISTB,ITAMPB,IOPTRES,
     &                           LABELA,FREQA,
     &                           WORK(KGAMMA1),WORK(KGAMMA2),
     &                           WORK(KGAMMA1EFF),WORK(KGAMMA2EFF),
     &                           IFADOTS,FACON,FILFA,ITRAN,
     &                           NFATRAN,MXVEC,WORK(KEND0),LEND0 )

        END IF

      END IF

*=====================================================================*
* final section: depending on IOPTRES store vector in memory or on
*                file or contract it with some other vectors:
*
*      the memory has to be organized as follows:
*         kgamma1     singles result vector
*         kgamma2     doubles result vector
*         kgammar12   r12 doubles result vector
*         kgamma1eff  effective singles result vector for CC3
*         kgamma2eff  effective doubles result vector for CC3
*         kend0       start of unused work space
*=====================================================================*
9999  CONTINUE 

      LEND0 = LWORK - KEND0

      IF (CCS) THEN
         MODELW = 'CCS       '
         IOPTW  = 1
      ELSE IF (CC2) THEN
         MODELW = 'CC2       '
         IOPTW  = 3
      ELSE IF (CCSD) THEN
         MODELW = 'CCSD      '
         IOPTW  = 3
      ELSE IF (CC3) THEN
         MODELW = 'CC3       '
         IOPTW  = 3
         IOPTWE = 24
      ELSE
         CALL QUIT('Unknown coupled cluster model in CC_FAMAT.')
      END IF
      IF (CCR12) THEN
        APROXR12 = '   '
        CALL CCSD_MODEL(MODELW,LENMOD,10,MODELW,10,APROXR12)
        IOPTWR12 = 32
      END IF


      IF (IOPTRES.EQ.0 .OR. IOPTRES.EQ.1) THEN
         CALL QUIT('IOPTRES=0,1 not implemented in CC_FAMAT.')

      ELSE IF (IOPTRES.EQ.3) THEN
       IFILE  = IFATRAN(4,ITRAN)
       IF (ILSTSYM(FILFA,IFILE).NE.ISYRES) THEN
         CALL QUIT('Symmetry mismatch for result vector in CC_FAMAT.')
       END IF
       CALL CC_WRRSP(FILFA,IFILE,ISYRES,IOPTW,MODELW,DUMMY,
     &               WORK(KGAMMA1),WORK(KGAMMA2),WORK(KEND0),LEND0)
       IF (CCR12) THEN
         CALL CC_WRRSP(FILFA,IFILE,ISYRES,IOPTWR12,
     &                 MODELW,DUMMY,DUMMY,WORK(KGAMMAR12),
     &                 WORK(KEND0),LEND0)
       END IF
       IF (CCSDT) THEN
         CALL CC_WRRSP(FILFA,IFILE,ISYRES,IOPTWE,MODELW,DUMMY,
     &             WORK(KGAMMA1EFF),WORK(KGAMMA2EFF),WORK(KEND0),LEND0)
       END IF
      ELSE IF (IOPTRES.EQ.4) THEN
       IFILE  = IFATRAN(4,ITRAN)
       IF (ILSTSYM(FILFA,IFILE).NE.ISYRES) THEN
         CALL QUIT('Symmetry mismatch for result vector in CC_FAMAT.')
       END IF
       CALL CC_WARSP(FILFA,IFILE,ISYRES,IOPTW,MODELW,DUMMY,
     &               WORK(KGAMMA1),WORK(KGAMMA2),WORK(KEND0),LEND0)
       IF (CCR12) THEN
         CALL CC_WARSP(FILFA,IFILE,ISYRES,IOPTWR12,
     &                 MODELW,DUMMY,DUMMY,WORK(KGAMMAR12),
     &                 WORK(KEND0),LEND0)
       END IF
       IF (CCSDT) THEN
        CALL CC_WARSP(FILFA,IFILE,ISYRES,IOPTWE,MODELW,DUMMY,
     &             WORK(KGAMMA1EFF),WORK(KGAMMA2EFF),WORK(KEND0),LEND0)
       END IF
      ELSE IF (IOPTRES.EQ.5) THEN
       IF (LOCDBG) THEN
         IVEC = 1
         WRITE(LUPRI,*) 'FACON TRIPLES CONTRIBUTION:'
         DO WHILE (IFADOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
            WRITE (LUPRI,*) 
     &            'FACON:',IVEC,ITRAN,FACON(IVEC,ITRAN),IOPTW
            IVEC = IVEC + 1
         END DO
       END IF                

       IF (.NOT.CCS) CALL CCLR_DIASCL(WORK(KGAMMA2),TWO,ISYRES)
       CALL CCDOTRSP(IFADOTS,FACON,IOPTW,FILFA,ITRAN,NFATRAN,MXVEC,
     &               WORK(KGAMMA1),WORK(KGAMMA2),ISYRES,
     &               WORK(KEND0),LEND0)

       IF (CCR12) THEN
         CALL CCDOTRSP(IFADOTS,FACON,IOPTWR12,FILFA,ITRAN,NFATRAN,MXVEC,
     &                 DUMMY,WORK(KGAMMAR12),ISYRES,
     &                 WORK(KEND0),LEND0)
       END IF

       IF (LOCDBG) THEN
         IVEC = 1
         DO WHILE (IFADOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
            WRITE (LUPRI,*) 
     &            'FACON:',IVEC,ITRAN,FACON(IVEC,ITRAN),IOPTW
            IVEC = IVEC + 1
         END DO
       END IF                
      ELSE
       CALL QUIT('Illegal value for IOPTRES in CC_FAMAT.')
      END IF

      CALL QEXIT('CC_FAMAT')
      RETURN
      END
*=====================================================================*
*                  END OF SUBROUTINE CC_FAMAT                         *
*=====================================================================*
*=====================================================================*
      SUBROUTINE CCSDT_FAMAT_NODDY(LISTL,IDLSTL,LISTB,IDLSTB,IOPTRES,
     &                             LABELA,FREQA,
     &                             OMEGA1,OMEGA2,
     &                             OMEGA1EFF,OMEGA2EFF,
     &                             IDOTS,DOTPROD,LISTDP,ITRAN,
     &                             NFATRAN,MXVEC,WORK,LWORK )
*---------------------------------------------------------------------*
*
*    Purpose: compute triples contribution to F{A} transformed vector
*
*    (F{A} T^B)^eff_1,2 = (F{A} T^B)_1,2(CCSD) + (F{A} T^B)_1,2(T^B_3)
*                               - (F{A} T^B)_3 A_3;1,2 (w_3 - w)^1 
*
*        
*     Written by Christof Haettig, April 2002 
*     based on CCSDT_ETA_NODDY
*     Extensions for cubic response, CCH, May 2003
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccfield.h>
#include <ccorb.h>
#include <dummy.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.FALSE.)

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      CHARACTER LABELA*8
      CHARACTER*3 LISTL, LISTDP, LISTB
      INTEGER LWORK, IDLSTL, IDLSTB, IOPTRES, ITRAN, MXVEC, NFATRAN
      INTEGER IDOTS(MXVEC,NFATRAN)

#if defined (SYS_CRAY)
      REAL DOTPROD(MXVEC,NFATRAN), DDOT
      REAL WORK(LWORK), FREQL, FREQA, FREQB, FREQF, FREQC
      REAL OMEGA1(*), OMEGA2(*)
      REAL OMEGA1EFF(*), OMEGA2EFF(*)
      REAL SIXTH, ONE, TWO, TCON, DCON, SCON, FF, SIGN
#else
      DOUBLE PRECISION DOTPROD(MXVEC,NFATRAN), DDOT
      DOUBLE PRECISION WORK(LWORK), FREQL, FREQA, FREQB, FREQE, FREQC
      DOUBLE PRECISION OMEGA1(*), OMEGA2(*)
      DOUBLE PRECISION OMEGA1EFF(*), OMEGA2EFF(*)
      DOUBLE PRECISION SIXTH, ONE, TWO, TCON, DCON, SCON, FF, SIGN
#endif
      PARAMETER(SIXTH=1.0D0/6.0D0, ONE=1.0D0, TWO=2.0D0)

      CHARACTER*10 MODEL
      LOGICAL L2INCL
      INTEGER INDEX, LUSIFC, IOPT, ISYMD, ILLL, IDEL, ISYDIS, NIJ, IJ,
     &        IVEC, IDLSTC, ISYMC, LUFOCK, ILSTSYM, ISYML
      INTEGER KSCR1, KFOCKD, KEND1, KT1AMP0, KLAMP0, KLAMH0,
     &        KINT1T0, KINT2T0, KINT1S0, KINT2S0, KXIAJB, KYIAJB,
     &        K0IOVVO, K0IOOVV, K0IOOOO, K0IVVVV, KOME1, KOME2, KDUM,
     &        KXINT, KEND2, LWRK2, KL1AM, KL2AM, KL3AM, KT3AM, KT2AM,
     &        KEND3, LWRK3, KINT1SC, KINT2SC, KLAMPC, KLAMHC,KFOCKAB,
     &        KFOCKC, LWRK1, KE3AM, KTC3AM, KTC1AM, KTC2AM,
     &        ISYMA, KINT1SB, KINT2SB, KLAMPB, KLAMHB, KFOCKB, ISYMB,
     &        KFOCKA, KFOCKA_AO, KFOCK0, IRREP, ISYM, IERR, KFOCKAB1,
     &        KFIELD, KFIELDAO, KT1AMB

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J 

      CALL QENTER('CCSDT_FAMAT_NODDY')

      KDUM = 1
      IF (DIRECT) CALL QUIT('CCSDT_FAMAT_NODDY: DIRECT NOT IMPLEMENTED')

      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'CCSDT_FAMAT_NODDY> left  vector:',LISTL,IDLSTL
        WRITE(LUPRI,*) 'CCSDT_FAMAT_NODDY> right vector:',LISTB,IDLSTB
        WRITE(LUPRI,*) 'CCSDT_FAMAT_NODDY> operator    :',LABELA
        WRITE(LUPRI,*) 'CCSDT_FAMAT_NODDY> frequency   :',FREQA
        WRITE(LUPRI,*) 'CCSDT_FAMAT_NODDY> result vector on entry:'
        CALL CC_PRP(OMEGA1,OMEGA2,1,1,1)
      END IF

*---------------------------------------------------------------------*
*     Memory allocation:
*---------------------------------------------------------------------*
      KSCR1   = 1
      KFOCKD  = KSCR1  + NT1AMX
      KFOCK0  = KFOCKD + NORBT
      KEND1   = KFOCK0 + NORBT*NORBT

      KFOCKA    = KEND1  
      KFOCKA_AO = KFOCKA    + NORBT*NORBT
      KEND1     = KFOCKA_AO + NORBT*NORBT

      IF (NONHF) THEN
        KFIELD   = KEND1
        KFIELDAO = KFIELD   + NORBT*NORBT
        KEND1    = KFIELDAO + NORBT*NORBT
      END IF

      KT1AMP0 = KEND1
      KLAMP0  = KT1AMP0 + NT1AMX
      KLAMH0  = KLAMP0  + NLAMDT
      KEND1   = KLAMH0  + NLAMDT

      KINT1T0 = KEND1
      KINT2T0 = KINT1T0 + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2T0 + NRHFT*NRHFT*NT1AMX

      KINT1S0 = KEND1
      KINT2S0 = KINT1S0 + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2S0 + NRHFT*NRHFT*NT1AMX

      KXIAJB  = KEND1
      KYIAJB  = KXIAJB  + NT1AMX*NT1AMX
      KEND1   = KYIAJB  + NT1AMX*NT1AMX

      K0IOVVO = KEND1
      K0IOOVV = K0IOVVO + NRHFT*NVIRT*NVIRT*NRHFT
      K0IOOOO = K0IOOVV + NRHFT*NVIRT*NVIRT*NRHFT
      K0IVVVV = K0IOOOO + NRHFT*NRHFT*NRHFT*NRHFT
      KEND1   = K0IVVVV + NVIRT*NVIRT*NVIRT*NVIRT 

      KOME1   = KEND1
      KOME2   = KOME1  + NT1AMX
      KEND1   = KOME2  + NT1AMX*NT1AMX

      KFOCKAB1= KEND1
      KFOCKAB = KFOCKAB1+ NORBT*NORBT
      KEND1   = KFOCKAB + NORBT*NORBT

      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_FAMAT_NODDY')
      ENDIF

*---------------------------------------------------------------------*
*     Read SCF orbital energies from file:
*---------------------------------------------------------------------*
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBT)
      CALL GPCLOSE(LUSIFC,'KEEP')

*---------------------------------------------------------------------*
*     Get zeroth-order Lambda matrices:
*---------------------------------------------------------------------*
      IOPT   = 1
      Call CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AMP0),WORK(KDUM))

      Call LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AMP0),
     &            WORK(KEND1),LWRK1)

*---------------------------------------------------------------------*
*     read zeroth-order AO Fock matrix from file: 
*---------------------------------------------------------------------*
      LUFOCK = -1
      CALL GPOPEN(LUFOCK,'CC_FCKH','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUFOCK)
      READ(LUFOCK) (WORK(KFOCK0-1+I),I=1,N2BST(ISYM0))
      CALL GPCLOSE(LUFOCK,'KEEP')

      CALL CC_FCKMO(WORK(KFOCK0),WORK(KLAMP0),WORK(KLAMH0),
     &              WORK(KEND1),LWRK1,ISYM0,ISYM0,ISYM0)

*---------------------------------------------------------------------*
*     If needed get external field:
*---------------------------------------------------------------------*
      IF (NONHF) THEN
        CALL DZERO(WORK(KFIELDAO),NORBT*NORBT)
        DO I = 1, NFIELD
          FF = EFIELD(I)
          CALL CC_ONEP(WORK(KFIELDAO),WORK(KEND1),LWRK1,FF,1,LFIELD(I))
        ENDDO
        CALL DCOPY(NORBT*NORBT,WORK(KFIELDAO),1,WORK(KFIELD),1)

        ! calculate external field in zero-order lambda basis
        CALL CC_FCKMO(WORK(KFIELD),WORK(KLAMP0),WORK(KLAMH0),
     *                WORK(KEND1),LWRK1,1,1,1)

        IF (LOCDBG) WRITE(LUPRI,*) 'NORM^2(FIELD):',
     &     DDOT(NORBT*NORBT,WORK(KFIELD),1,WORK(KFIELD),1)
      ENDIF

*---------------------------------------------------------------------*
*     Get property integrals and transform them to the MO basis:
*---------------------------------------------------------------------*
      ISYMA = 1 ! since this code is limited to C1 symmetry...

      CALL CCPRPAO(LABELA,.TRUE.,WORK(KFOCKA_AO),IRREP,ISYM,IERR,
     &             WORK(KEND1),LWRK1)
      IF (IERR.NE.0 .OR. IRREP.NE.ISYMA) THEN
       CALL QUIT('CCSDT_FA_NODDY: error reading operator '//LABELA)
      END IF

      CALL DCOPY(NORBT*NORBT,WORK(KFOCKA_AO),1,WORK(KFOCKA),1)

      CALL CC_FCKMO(WORK(KFOCKA),WORK(KLAMP0),WORK(KLAMH0),
     &              WORK(KEND1),LWRK1,ISYMA,ISYM0,ISYM0)

*---------------------------------------------------------------------*
*     Compute some integral intermediates:
*---------------------------------------------------------------------*

      CALL DZERO(WORK(KINT1T0),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2T0),NRHFT*NRHFT*NT1AMX)

      CALL DZERO(WORK(KINT1S0),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2S0),NRHFT*NRHFT*NT1AMX)

      CALL DZERO(WORK(KXIAJB), NT1AMX*NT1AMX)
      CALL DZERO(WORK(KYIAJB), NT1AMX*NT1AMX)

      CALL DZERO(WORK(K0IOVVO),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(K0IOOVV),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(K0IOOOO),NRHFT*NRHFT*NRHFT*NRHFT)
      CALL DZERO(WORK(K0IVVVV),NVIRT*NVIRT*NVIRT*NVIRT )


      DO ISYMD = 1, NSYM
         DO ILLL = 1,NBAS(ISYMD)
            IDEL   = IBAS(ISYMD) + ILLL
            ISYDIS = MULD2H(ISYMD,ISYMOP)
 
C           ----------------------------
C           Work space allocation no. 2.
C           ----------------------------
            KXINT  = KEND1
            KEND2  = KXINT + NDISAO(ISYDIS)
            LWRK2  = LWORK - KEND2
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
               CALL QUIT('Insufficient space in CCSDT_FAMAT_NODDY')
            ENDIF
 
C           ---------------------------
C           Read in batch of integrals.
C           ---------------------------
            CALL CCRDAO(WORK(KXINT),IDEL,1,WORK(KEND2),LWRK2,
     *                  WORK(KDUM),DIRECT)
 
C           ----------------------------------
C           Calculate integrals needed in CC3:
C           ----------------------------------
            CALL CCSDT_TRAN1(WORK(KINT1T0),WORK(KINT2T0),
     &                       WORK(KLAMP0),WORK(KLAMH0),
     &                       WORK(KXINT),IDEL)

            CALL CC3_TRAN2(WORK(KXIAJB),WORK(KYIAJB),
     &                     WORK(KLAMP0),WORK(KLAMH0),
     &                     WORK(KXINT),IDEL)

            CALL CCSDT_TRAN3(WORK(KINT1S0),WORK(KINT2S0),WORK(KLAMP0),
     &                       WORK(KLAMH0),WORK(KXINT),IDEL)

            CALL CCFOP_TRAN1_R(WORK(K0IOVVO),WORK(K0IOOVV),
     &                         WORK(K0IOOOO),WORK(K0IVVVV),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KXINT),IDEL)

         END DO   
      END DO  

*---------------------------------------------------------------------*
*     Some more memory allocations:
*---------------------------------------------------------------------*
      KL1AM  = KEND1
      KL2AM  = KL1AM + NT1AMX
      KL3AM  = KL2AM + NT1AMX*NT1AMX
      KEND2  = KL3AM + NT1AMX*NT1AMX*NT1AMX
      LWRK2  = LWORK - KEND2

      KT3AM  = KEND2
      KT2AM  = KT3AM + NT1AMX*NT1AMX*NT1AMX
      KEND3  = KT2AM + NT1AMX*NT1AMX
      LWRK3  = LWORK - KEND3

      IF (LWRK3 .LT. NT2AMX) THEN
         CALL QUIT('Insufficient space in CCSDT_FAMAT_NODDY')
      ENDIF
 
      ! read T^B doubles amplitudes from file and square up
      ISYMB = ILSTSYM(LISTB,IDLSTB)
      IOPT  = 2
      Call CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &              WORK(KDUM),WORK(KEND3))
      Call CCLR_DIASCL(WORK(KEND3),TWO,ISYMB) 
      CALL CC_T2SQ(WORK(KEND3),WORK(KT2AM),ISYMB)

      ! read L^0 multipliers from file and square up doubles part
      ISYML = ILSTSYM(LISTL,IDLSTL)
      IOPT  = 3
      Call CC_RDRSP(LISTL,IDLSTL,ISYML,IOPT,MODEL,
     &              WORK(KL1AM),WORK(KT3AM))
      CALL CC_T2SQ(WORK(KT3AM),WORK(KL2AM),ISYM0)   

*---------------------------------------------------------------------*
*     Compute triples amplitude response:
*---------------------------------------------------------------------*
      KINT1SB = KEND3
      KINT2SB = KINT1SB + NT1AMX*NVIRT*NVIRT
      KEND3   = KINT2SB + NRHFT*NRHFT*NT1AMX

      KT1AMB  = KEND3
      KLAMPB  = KT1AMB + NT1AMX
      KLAMHB  = KLAMPB + NLAMDT
      KFOCKB  = KLAMHB + NLAMDT
      KEND3   = KFOCKB + NORBT*NORBT

      LWRK3  = LWORK  - KEND3
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_FAMAT_NODDY')
      ENDIF

      IF      (LISTB(1:3).EQ.'R1 ' .OR. LISTB(1:3).EQ.'RE ' .OR.
     &         LISTB(1:3).EQ.'RC '                              ) THEN

        CALL CCSDT_T31_NODDY(WORK(KT3AM),LISTB,IDLSTB,FREQB,.FALSE.,
     &                       .FALSE.,WORK(KINT1S0),WORK(KINT2S0),
     &                       .FALSE.,WORK(KDUM),WORK(KDUM),
     &                       .FALSE.,WORK(KDUM),WORK(KDUM),
     &                               WORK(KINT1SB),WORK(KINT2SB),
     &                       WORK(KLAMPB),WORK(KLAMHB),WORK(KFOCKB),
     &                       WORK(KLAMP0),WORK(KLAMH0),WORK(KFOCK0),
     &                       WORK(KDUM),WORK(KFOCKD),
     &                       WORK(KEND3),LWRK3)

      ELSE IF (LISTB(1:3).EQ.'R2 ' .OR. LISTB(1:3).EQ.'ER1') THEN

        CALL CCSDT_T32_NODDY(WORK(KT3AM),LISTB,IDLSTB,FREQB,
     &                       WORK(KINT1S0),WORK(KINT2S0),
     &                       WORK(KLAMP0),WORK(KLAMH0),WORK(KFOCK0),
     &                       WORK(KFOCKD),WORK(KFIELDAO),WORK(KFIELD),
     &                       WORK(KSCR1),WORK(KEND3),LWRK3)

        IOPT = 1
        Call CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,WORK(KT1AMB),DUMMY)

        CALL CCLR_LAMTRA(WORK(KLAMP0),WORK(KLAMPB),WORK(KLAMH0),
     &                   WORK(KLAMHB),WORK(KT1AMB),ISYMB)
      ELSE
        
        CALL QUIT('Unknown list '//LISTB//' in CCSDT_FAMAT_NODDY.')
    
      END IF

*---------------------------------------------------------------------*
*     Compute triples multipliers L_3:
*---------------------------------------------------------------------*
      IF (LISTL(1:3).EQ.'L0 ') THEN

        FREQL = 0.0D0

        CALL DZERO(WORK(KL3AM),NT1AMX*NT1AMX*NT1AMX)

        IF (NONHF .AND. LWRK3.LT.NT1AMX*NT1AMX*NT1AMX)
     *    CALL QUIT('Out of memory in CCSDT_FAMAT_NODDY.')

        ! remember that CCSDT_L03AM returns -L3 !!
        CALL CCSDT_L03AM(WORK(KL3AM),WORK(KINT1T0),WORK(KINT2T0),
     *                   WORK(KXIAJB),WORK(KFOCK0),WORK(KL1AM),
     *                   WORK(KL2AM),WORK(KSCR1),WORK(KFOCKD),
     *                   WORK(KFIELD),WORK(KEND3))

        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KL3AM),1)

      ELSE IF (LISTL(1:3).EQ.'L1 ' .OR. LISTL(1:3).EQ.'LE ' .OR.
     &         LISTL(1:3).EQ.'M1 ' .OR. LISTL(1:3).EQ.'N2 ' .OR.
     &         LISTL(1:3).EQ.'E0 '                              ) THEN

        CALL DZERO(WORK(KL3AM),NT1AMX*NT1AMX*NT1AMX)

        CALL CCSDT_TBAR31_NODDY(WORK(KL3AM),FREQL,LISTL,IDLSTL,
     &                        WORK(KLAMP0),WORK(KLAMH0),
     &                        WORK(KFOCK0),WORK(KFOCKD),WORK(KSCR1),
     &                        WORK(KXIAJB),WORK(KINT1T0),WORK(KINT2T0),
     &                        WORK(KEND3),LWRK3)

        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KL3AM),1)

      ELSE

        ! FREQL = ??

        CALL QUIT('CCSDT_FAMAT_NODDY> LISTL NOT AVAILABLE:'//LISTL)
      
      END IF

*---------------------------------------------------------------------*
*     Compute contribution from <L_3|[[A,T^B_3],\tau_nu_1|HF>:
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KOME1),NT1AMX)

      CALL CCSDT_E1AM(WORK(KOME1),WORK(KL3AM),WORK(KT3AM),WORK(KFOCKA))

      DO I = 1,NT1AMX
         OMEGA1(I) = OMEGA1(I) + WORK(KOME1+I-1)
      END DO

ctest
C     LUFOCK = -1
C     CALL GPOPEN(LUFOCK,'CCTEST'//LISTL(1:2)//LISTB(1:2),
C    &            'UNKNOWN',' ','FORMATTED',IDUMMY,.FALSE.)
C     REWIND(LUFOCK)
C     WRITE(LUFOCK,'(A)') 'OMEGA1:'
C     WRITE(LUFOCK,'(F20.16)') (WORK(KOME1-1+I),I=1,NT1AMX)
C     WRITE(LUFOCK,'(A)') 'FOCKA:'
C     WRITE(LUFOCK,'(F20.16)') (WORK(KFOCKA-1+I),I=1,NORBT*NORBT)
C     WRITE(LUFOCK,'(A)') 'L3AM:'
C     WRITE(LUFOCK,'(F20.16)')(WORK(KL3AM-1+I),I=1,NT1AMX*NT1AMX*NT1AMX)
C     WRITE(LUFOCK,'(A)') 'T3AM:'
C     WRITE(LUFOCK,'(F20.16)')(WORK(KT3AM-1+I),I=1,NT1AMX*NT1AMX*NT1AMX)
C     CALL GPCLOSE(LUFOCK,'KEEP')
ctest
      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'CCSDT_FAMAT_NODDY> Contribution to F{A} T^B:'
        CALL CC_PRP(WORK(KOME1),WORK,1,1,0)
      END IF

*---------------------------------------------------------------------*
*     Compute contribution from <L_3|[[A,T^B_2],\tau_nu_2]|HF>
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KOME2),NT1AMX*NT1AMX)

      CALL CCSDT_E2AM(WORK(KOME2),WORK(KL3AM),WORK(KT2AM),WORK(KFOCKA))

      DO I = 1,NT1AMX
         DO J = 1,I
            IJ = NT1AMX*(I-1) + J
            NIJ = INDEX(I,J)
            OMEGA2(NIJ) = OMEGA2(NIJ) + WORK(KOME2+IJ-1)
         END DO
      END DO

*---------------------------------------------------------------------*
*     Compute [A,T^B_1] by AO-to-MO transformation of A with
*     the response Lambda matrices:
*---------------------------------------------------------------------*
      CALL DCOPY(NORBT*NORBT,WORK(KFOCKA_AO),1,WORK(KFOCKAB1),1)
      CALL DCOPY(NORBT*NORBT,WORK(KFOCKA_AO),1,WORK(KFOCKAB),1)

      CALL CC_FCKMO(WORK(KFOCKAB1),WORK(KLAMPB),WORK(KLAMH0),
     &              WORK(KEND3),LWRK3,ISYMA,ISYMB,ISYM0)

      CALL CC_FCKMO(WORK(KFOCKAB),WORK(KLAMP0),WORK(KLAMHB),
     &              WORK(KEND3),LWRK3,ISYMA,ISYM0,ISYMB)
     
      CALL DAXPY(NORBT*NORBT,ONE,WORK(KFOCKAB1),1,WORK(KFOCKAB),1)

*---------------------------------------------------------------------*
*     Compute triples result vector 
*       <L_3|[[A,T^B_1],\tau_nu_3]|HF> ,
*---------------------------------------------------------------------*
      ! overwrite T3 vector
      KE3AM  = KT3AM
  
      CALL DZERO(WORK(KE3AM),NT1AMX*NT1AMX*NT1AMX)

      L2INCL = .FALSE.
      CALL CCSDT_E3AM(WORK(KE3AM),WORK(KDUM),WORK(KL3AM),
     &                WORK(KFOCKAB),L2INCL)

*---------------------------------------------------------------------*
*     Now we split:
*       for IOPTRES < 5 we compute the effective result vector
*       for IOPTRES = 5 we compute the contractions F{A} T^B T^C
*---------------------------------------------------------------------*
      IF (IOPTRES.GE.1 .AND. IOPTRES.LE.4) THEN

        IOPT  = 2
        Call CC_RDRSP('R0 ',0,ISYM0,IOPT,MODEL,WORK(KDUM),WORK(KEND3))
        CALL CC_T2SQ(WORK(KEND3),WORK(KT2AM),ISYM0)

        CALL DCOPY(NT1AMX,OMEGA1,1,OMEGA1EFF,1)
        CALL DCOPY(NT2AMX,OMEGA2,1,OMEGA2EFF,1)

        FREQE = FREQL + FREQA + FREQB

        CALL CC_LHPART_NODDY(OMEGA1EFF,OMEGA2EFF,WORK(KE3AM),-FREQE,
     &                       WORK(KFOCKD),WORK(KFIELD),
     &                       WORK(K0IOOOO),WORK(K0IOVVO),
     &                       WORK(K0IOOVV),WORK(K0IVVVV),
     &                       WORK(KT2AM),WORK(KINT1S0),WORK(KINT2S0),
     &                       WORK(KEND3),LWRK3)

      ELSE IF (IOPTRES.EQ.5) THEN

        SIGN = +1.0D0

        CALL CCDOTRSP_NODDY(WORK(KOME1),WORK(KOME2),WORK(KE3AM),SIGN,
     &                      ITRAN,LISTDP,IDOTS,DOTPROD,MXVEC,
     &                      WORK(KLAMP0),WORK(KLAMH0),
     &                      WORK(KFOCK0),WORK(KFOCKD),
     &                      WORK(KXIAJB), WORK(KYIAJB),
     &                      WORK(KINT1T0),WORK(KINT2T0),
     &                      WORK(KINT1S0),WORK(KINT2S0),
     &                      'CCSDT_FAMAT_NODDY',LOCDBG,LOCDBG,
     &                      .FALSE.,WORK(KEND3),LWRK3)

      ELSE
        CALL QUIT('Illegal value for IOPTRES IN CCSDT_FAMAT_NODDY')
      END IF

      CALL QEXIT('CCSDT_FAMAT_NODDY')

      RETURN
      END 
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_FAMAT_NODDY                    *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_FA_SETUP(IFATRAN,IFADOTS,NFATRAN,MXVEC,
     &                          IDXL_FADEN,IDXB_FADEN,IDXC_FADEN,
     &                          NFADEN,MXFADEN,
     &                          IDXL_INTER,IDXR_INTER,LSTR_INTER,
     &                          NINTER,MXINTER,
     &                          LISTL,LISTO,LISTB,LISTC)
*---------------------------------------------------------------------*
*
* Purpose: setup loop structures to compute intermediates and 
*          densities needed to calculate contractions of the
*          form F{A} t^B t^C
*
* Christof Haettig, May 2003
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclists.h>
#include <ccroper.h>

      CHARACTER*(25) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCSDT_FA_SETUP> ')  
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
 
* input:
      CHARACTER*3 LISTL, LISTB, LISTC, LISTO
      INTEGER MXFADEN, NFATRAN, MXVEC, MXINTER
      INTEGER IFATRAN(MXDIM_FATRAN,NFATRAN), IFADOTS(MXVEC,NFATRAN)

* output:
      CHARACTER*3 LSTR_INTER(*)
      INTEGER IDXL_FADEN(*), IDXB_FADEN(*), IDXC_FADEN(*),
     &        IDXL_INTER(*), IDXR_INTER(*)
      INTEGER NFADEN, NINTER

* local:
      CHARACTER*8 LABELA
      LOGICAL LPDBSA
      INTEGER ITRAN, IDLSTL, IOPERA, IDLSTB, IFILE, IRELAX, ISYRES,
     &        IVEC, IDLSTC, ISYMTC, ISYCTR, ISYMTA, ISYMTB, 
     &        IFADEN, IDXLB, IDXLC, IDX

* external functions:
      INTEGER ILSTSYM

C     ------
C     begin:
C     ------
C     IF (LISTL.NE.'L0') THEN
C       CALL QUIT('CCSDT_FA_DEN can not yet treat '//LISTL//
C    &            ' type vectors.')
C     END IF

C     ------------------------------------------------
C     set up list of effective F{A} t^B t^C densities:
C     ------------------------------------------------
      NFADEN = 0
      DO ITRAN = 1, NFATRAN

        IDLSTL = IFATRAN(1,ITRAN)
        IOPERA = IFATRAN(2,ITRAN)
        IDLSTB = IFATRAN(3,ITRAN)
        IFILE  = IFATRAN(4,ITRAN)
        IRELAX = IFATRAN(5,ITRAN)

        LABELA = LBLOPR(IOPERA)
        LPDBSA = LPDBSOP(IOPERA)

        ISYCTR = ILSTSYM(LISTL,IDLSTL)
        ISYMTA = ILSTSYM(LISTO,IOPERA)
        ISYMTB = ILSTSYM(LISTB,IDLSTB)

        ISYRES = MULD2H(MULD2H(ISYMTA,ISYMTB),ISYCTR)

        IF ( IRELAX.GT.0 .OR. LPDBSA ) THEN
          CALL QUIT('Relaxed perturbations not yet implemented in '//
     &              'CCSDT_FA_DEN.')
        END IF

        IVEC = 1
        DO WHILE(IFADOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC) 
          IDLSTC = IFADOTS(IVEC,ITRAN)
          ISYMTC = ILSTSYM(LISTC,IDLSTC)
 
          IF (ISYMTC.NE.ISYRES) THEN
            CALL QUIT('Symmetry mismatch in CCSDT_FA_DEN')
          END IF

          ! check if density is already on to-do list
          IFADEN = 0
          DO IDX = 1, NFADEN
            IF (IDLSTL .EQ.IDXL_FADEN(IDX) .AND.
     &          IDLSTB .EQ.IDXB_FADEN(IDX) .AND.
     &          IDLSTC .EQ.IDXC_FADEN(IDX)       ) THEN
              IFADEN = IDX
            END IF
            IF (IFADEN.EQ.0 .AND. LISTB.EQ.LISTC .AND.
     &          IDLSTL .EQ.IDXL_FADEN(IDX) .AND.
     &          IDLSTC .EQ.IDXB_FADEN(IDX) .AND.
     &          IDLSTB .EQ.IDXC_FADEN(IDX)       ) THEN
              IFADEN = IDX
            END IF
          END DO

          ! if not found, put on to-do list for densities
          IF (IFADEN.EQ.0) THEN
            NFADEN = NFADEN + 1
            IF (NFADEN.GT.MXFADEN) CALL QUIT('NFADEN out of range')
            IDXL_FADEN(NFADEN) = IDLSTL
            IDXB_FADEN(NFADEN) = IDLSTB
            IDXC_FADEN(NFADEN) = IDLSTC
          END IF

          IVEC = IVEC + 1
        END DO

      END DO

C     ------------------------------------------------
C     set up list of intermediates that depend on 
C     pairs (L,B) and (L,C):
C     ------------------------------------------------
      NINTER = 0
      DO IFADEN = 1, NFADEN
        
        ! check if pair (L,B) is already on the list
        IDXLB = 0
        DO IDX = 1, NINTER
          IF (IDXL_FADEN(IFADEN).EQ.IDXL_INTER(IDX) .AND.
     &        IDXB_FADEN(IFADEN).EQ.IDXR_INTER(IDX) .AND.
     &                    LISTB .EQ.LSTR_INTER(IDX)      ) THEN
            IDXLB = IDX
          END IF
        END DO

        ! if not found, put on to-do list for intermediates
        IF (IDXLB.EQ.0) THEN
          NINTER = NINTER + 1
          IF (NINTER.GT.MXINTER) CALL QUIT('NINTER out of range')
          IDXL_INTER(NINTER) = IDXL_FADEN(IFADEN)
          IDXR_INTER(NINTER) = IDXB_FADEN(IFADEN)
          LSTR_INTER(NINTER) = LISTB
        END IF

        ! check if pair (L,C) is already on the list
        IDXLC = 0
        DO IDX = 1, NINTER
          IF (IDXL_FADEN(IFADEN).EQ.IDXL_INTER(IDX) .AND.
     &        IDXC_FADEN(IFADEN).EQ.IDXR_INTER(IDX) .AND.
     &                    LISTC .EQ.LSTR_INTER(IDX)      ) THEN
            IDXLC = IDX
          END IF
        END DO

        ! if not found, put on to-do list for intermediates
        IF (IDXLC.EQ.0) THEN
          NINTER = NINTER + 1
          IF (NINTER.GT.MXINTER) CALL QUIT('NINTER out of range')
          IDXL_INTER(NINTER) = IDXL_FADEN(IFADEN)
          IDXR_INTER(NINTER) = IDXC_FADEN(IFADEN)
          LSTR_INTER(NINTER) = LISTC
        END IF

      END DO

C     -------------------------
C     if requested print lists:
C     -------------------------
      IF (LOCDBG) THEN
        WRITE(LUPRI,'(//,A)') 'list of F{A} T^B T^C densities:'
        WRITE(LUPRI,'(A)') '-------------------------------'
        WRITE(LUPRI,'(5X,A)') ' IDX   L    B    C   '
        DO IFADEN = 1, NFADEN
          WRITE(LUPRI,'(5X,4I5)') IFADEN, IDXL_FADEN(IFADEN),
     &             IDXB_FADEN(IFADEN), IDXC_FADEN(IFADEN)
        END DO
        WRITE(LUPRI,'(//,A)') 'list intermediates:'
        WRITE(LUPRI,'(A)')    '-------------------'
        WRITE(LUPRI,'(5X,A)') ' IDX   L       R     '
        DO IFADEN = 1, NINTER
          WRITE(LUPRI,'(5X,2I5,3X,A3,I3)') IFADEN, IDXL_INTER(IFADEN),
     &             LSTR_INTER(IFADEN), IDXR_INTER(IFADEN)
        END DO
        WRITE(LUPRI,'(//)')
      END IF

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_FA_SETUP                       *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_FA_DEN(LISTL,LISTO,LISTB,LISTC,FADEN_NODDY,
     &                        NFATRAN,MXVEC,IFATRAN,IFADOTS,FACON,
     &                        FNDELD,FNCKJD,FNDKBC,FNTOC,
     &                        FN3VI,FNDKBC3,FN3FOPX,FN3FOP2X,
     &                        IDXL_INTER,IDXR_INTER,
     &                        LSTR_INTER,NINTETA,
     &                        IDXL_FADEN,IDXB_FADEN,
     &                        IDXC_FADEN,NFADEN,
     &                        IADR_INTER,IADR_DEN,
     &                        WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: compute triples contributions to to F{A} matrix 
*              contractions via effective densities
*
*     Written by Christof Haettig, Mai 2003, Friedrichstal
*
*=====================================================================*
      IMPLICIT NONE
#include <ccsdsym.h>
#include <priunit.h>
#include <ccorb.h>
#include <dummy.h>
#include <ccr1rsp.h>
#include <ccroper.h>
#include <ccexci.h>
#include <cclists.h>

      LOGICAL LOCDBG
      PARAMETER ( LOCDBG = .FALSE. )

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      CHARACTER FNINTER*12, FNR2TP*12, FNDEFF*12
      PARAMETER (FNINTER='CCFADENINTER', FNR2TP = 'CCFADEN_R2TP',
     *           FNDEFF ='CCFADEN_DEFF')

* input/output:
      CHARACTER*3 LISTL, LISTO, LISTB, LISTC, LSTR_INTER(*)
      CHARACTER*(*) FNDELD, FNCKJD, FNDKBC, FNTOC, FN3VI, FNDKBC3
      CHARACTER*(*) FN3FOPX, FN3FOP2X 
      LOGICAL FADEN_NODDY
      INTEGER NINTETA, NFADEN, MXVEC, NFATRAN, LWORK
      INTEGER IDXL_INTER(*), IDXR_INTER(*), IDXL_FADEN(*),
     *        IDXB_FADEN(*), IDXC_FADEN(*),
     *        IADR_INTER(3,*), IADR_DEN(*), 
     *        IFADOTS(MXVEC,*), IFATRAN(MXDIM_FATRAN,*)

#if defined (SYS_CRAY)
      REAL FACON(MXVEC,NFATRAN)
      REAL WORK(*)
#else
      DOUBLE PRECISION FACON(MXVEC,NFATRAN)
      DOUBLE PRECISION WORK(*)
#endif

* local:
      CHARACTER MODEL*(10), LISTR*3, LABELR*8, LABELA*8
      LOGICAL LPDBSA, CALL_ETA_FA_DEN, NEED_L1AM, NEED_L2TP,
     &        NEED_T2TP, NEED_R2TP
      INTEGER ISYM, ILEN, ISYMFN, ISYMIM, NIMFN(8)
      INTEGER LUFOCK, LUINTER, LUR2TP,
     &        LUDELD, LUDKBC, LUTOC, LU3VI,
     &        LUDKBC3, LU3FOPX, LU3FOP2X, LUDEFF, LUCKJD
      INTEGER KEND1, KT1AM, KLAMP0, KLAMH0, KFOCK0, KT2TP, KL1AM,
     &        KL2TP, IOPT, IDLSTL, ISYML, KEND2, LWRK2, IADRINT,
     &        IDLSTR, ISYMR, INTETA, KR2TP, KEND3, LWRK3, ISYETA,
     &        KDAB, KDIJ, KMMAT, ISYDEN, KDIA1, KDIA2, KDIA3, KDIA4,
     &        KDIA5, KDIA6, IFADEN, IADRIA, ITRAN, IOPERA, IDLSTB,
     &        IFILE, IRELAX, ISYMA, ISYMB, ISYRES, IVEC, IDLSTC,
     &        ISYMC, KFOCK, KFOCKIA, IRREP, IERR, IMAT, LWRK1, IDX,
     &        KEND0

#if defined (SYS_CRAY)
      REAL FREQR, TRIPCON
      REAL TWO
#else
      DOUBLE PRECISION FREQR, TRIPCON
      DOUBLE PRECISION TWO
#endif
      PARAMETER ( TWO = 2.0D0 )

* external functions:
      INTEGER ILSTSYM
#if defined (SYS_CRAY)
      REAL DDOT
#else
      DOUBLE PRECISION DDOT
#endif


*---------------------------------------------------------------------*
*     some initializations:
*---------------------------------------------------------------------*
      DO ISYM = 1, NSYM
        ILEN = 0
        DO ISYMFN = 1, NSYM
          ISYMIM = MULD2H(ISYM,ISYMFN)
          ILEN   = ILEN + NMATIJ(ISYMIM)*NT1AM(ISYMFN)
        END DO
        NIMFN(ISYM) = ILEN
      END DO 

      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'Entered CCSDT_FA_DEN...'
        WRITE(LUPRI,*) 'LISTL,LISTB,LISTC:',LISTL,LISTB,LISTC
      END IF
*---------------------------------------------------------------------*
*     allocate some zeroth-order intermediates:
*---------------------------------------------------------------------*
      KEND0 = 1

      KT1AM  = KEND0
      KLAMP0 = KT1AM  + NT1AMX
      KLAMH0 = KLAMP0 + NLAMDT
      KEND0  = KLAMH0 + NLAMDT

      KFOCK0 = KEND0
      KEND1  = KFOCK0 + N2BAST
 
      LWRK1 = LWORK - KEND1

      IF (LWRK1 .LT. 0) THEN
        CALL QUIT('Insufficient work space CCSDT_FA_DEN (1)')
      END IF

*---------------------------------------------------------------------*
*     read zeroth-order amplitudes and compute lambda matrices:
*---------------------------------------------------------------------*
      IOPT = 1
      CALL CC_RDRSP('R0 ',0,ISYM0,IOPT,MODEL,WORK(KT1AM),DUMMY) 
      
      CALL LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AM),
     &            WORK(KEND1),LWRK1)

*---------------------------------------------------------------------*
*     get zeroth-order AO Fock in Lambda-MO basis:
*---------------------------------------------------------------------*
      LUFOCK = -1
      CALL GPOPEN(LUFOCK,'CC_FCKH','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUFOCK)
      READ(LUFOCK) (WORK(KFOCK0-1+I),I=1,N2BST(ISYM0))
      CALL GPCLOSE(LUFOCK,'KEEP')
 
      CALL CC_FCKMO(WORK(KFOCK0),WORK(KLAMP0),WORK(KLAMH0),
     &              WORK(KEND1),LWRK1,ISYM0,ISYM0,ISYM0) 

*---------------------------------------------------------------------*
*     open some files:
*---------------------------------------------------------------------*
      LUINTER = -1
      CALL WOPEN2(LUINTER,FNINTER,64,0)

      LUDELD  = -1
      LUCKJD  = -1
      LUDKBC  = -1
      LUTOC   = -1
      LU3VI   = -1
      LUDKBC3 = -1
      LU3FOPX = -1
      LU3FOP2X= -1
 
      CALL WOPEN2(LUDELD,FNDELD,64,0)
      CALL WOPEN2(LUCKJD,FNCKJD,64,0)
      CALL WOPEN2(LUDKBC,FNDKBC,64,0)
      CALL WOPEN2(LUTOC,FNTOC,64,0)
      CALL WOPEN2(LU3VI,FN3VI,64,0)
      CALL WOPEN2(LUDKBC3,FNDKBC3,64,0)
      CALL WOPEN2(LU3FOPX,FN3FOPX,64,0)
      CALL WOPEN2(LU3FOP2X,FN3FOP2X,64,0)  

*---------------------------------------------------------------------*
*     loop over the intermediates that need to be computed:
*---------------------------------------------------------------------*
      IADRINT = 1

      DO INTETA = 1, NINTETA

        IF (LISTL.EQ.'L0 ') THEN
          IDLSTL = 0
          ISYML  = 1
        ELSE
          IDLSTL = IDXL_INTER(INTETA)
          ISYML  = ILSTSYM(LISTL,IDLSTL)
        END IF
 
        LISTR  = LSTR_INTER(INTETA)
        IDLSTR = IDXR_INTER(INTETA)
        ISYMR  = ILSTSYM(LISTR,IDLSTR)

        ! check for modules which needed extra arrays / preparations
        CALL_ETA_FA_DEN = ( (LISTL.EQ.'L0 ' .OR. LISTL.EQ.'LE ') .AND. 
     *                     (LISTR.EQ.'R1 ' .OR. LISTR.EQ.'RE '
     *                     .OR. LISTR.EQ.'R2 ') )

        NEED_L1AM = CALL_ETA_FA_DEN
        NEED_L2TP = CALL_ETA_FA_DEN
        NEED_T2TP = CALL_ETA_FA_DEN
        NEED_R2TP = CALL_ETA_FA_DEN


c       -------------------------------
c       allocate space for multipliers:
c       -------------------------------
        KEND2 = KEND1
    
        IF (NEED_T2TP) THEN
          KT2TP = KEND2
          KEND2 = KT2TP + NT2SQ(ISYM0)
        END IF

        IF (NEED_L1AM) THEN
          KL1AM = KEND2 
          KEND2 = KL1AM + NT1AM(ISYML)
        END IF
          
        IF (NEED_L2TP) THEN
          KL2TP = KEND2
          KEND2 = KL2TP + NT2SQ(ISYML)
        END IF

        LWRK2 = LWORK - KEND2
        IF (LWRK2 .LT. NT2SQ(ISYML)) THEN
          CALL QUIT('Insufficient work space CCSDT_FA_DEN (3)')
        END IF

c       ------------------------------------------------
c       read and resort zeroth-order doubles amplitudes:
c       ------------------------------------------------
        IF (NEED_T2TP) THEN
          IF (LWRK2 .LT. NT2SQ(ISYM0)) 
     &      CALL QUIT('Insufficient work space CCSDT_FA_DEN (2a)')
          
          IOPT = 2
          CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,DUMMY,WORK(KT2TP))
          CALL CC_T2SQ(WORK(KT2TP),WORK(KEND2),ISYM0)
          CALL CC3_T2TP(WORK(KT2TP),WORK(KEND2),ISYM0)
          IF (LOCDBG) WRITE(LUPRI,*) 'Norm of T2TP ',
     *        DDOT(NT2SQ(ISYM0),WORK(KT2TP),1,WORK(KT2TP),1)
        END IF

c       -------------------------
c       read singles multipliers:
c       -------------------------
        IF (NEED_L1AM) THEN
          IOPT = 1
          CALL CC_RDRSP(LISTL,IDLSTL,ISYML,IOPT,MODEL,WORK(KL1AM),DUMMY)
        END IF
       
c       ------------------------------------
c       read and resort doubles multipliers:
c       ------------------------------------
        IF (NEED_L2TP) THEN
          IF (LWRK2 .LT. NT2SQ(ISYM0)) 
     &      CALL QUIT('Insufficient work space CCSDT_FA_DEN (2b)')

          IOPT = 3
          CALL CC_RDRSP(LISTL,IDLSTL,ISYML,IOPT,MODEL,
     *                  WORK(KL1AM),WORK(KL2TP))
          CALL CC_T2SQ(WORK(KL2TP),WORK(KEND2),ISYML)
          CALL CC3_T2TP(WORK(KL2TP),WORK(KEND2),ISYML)
          IF (LOCDBG) WRITE(LUPRI,*) 'Norm of L2TP ',
     *        DDOT(NT2SQ(ISYML),WORK(KL2TP),1,WORK(KL2TP),1)  
        END IF

c       ---------------------------------------------
c       read and resort doubles response multipliers:
c       ---------------------------------------------
        IF (NEED_R2TP) THEN
          LUR2TP = -1
          CALL WOPEN2(LUR2TP,FNR2TP,64,0)

          KR2TP = KEND2
          KEND3 = KR2TP + NT2SQ(ISYMR)
          LWRK3 = LWORK - KEND3
          IF (LWRK3.LT.NT2AM(ISYMR))
     &      CALL QUIT('Out of memory in CCSDT_FA_DEN. (4a)')
         
          IOPT = 2
          CALL CC_RDRSP(LISTR,IDLSTR,ISYMR,IOPT,MODEL,
     &                  DUMMY,WORK(KEND3))
          CALL CCLR_DIASCL(WORK(KEND3),TWO,ISYMR)
          CALL CC_T2SQ(WORK(KEND3),WORK(KR2TP),ISYMR)
         
          CALL PUTWA2(LUR2TP,FNR2TP,WORK(KR2TP),1,NT2SQ(ISYMR))
        END IF

c       ------------------------------------------------------
c       allocate memory for intermediates and initialize them:
c       ------------------------------------------------------
        ISYETA = MULD2H(ISYML,ISYMR)

        KDAB  = KEND2
        KDIJ  = KDAB  + NMATAB(ISYETA)
        KMMAT = KDIJ  + NMATIJ(ISYETA)
        KEND3 = KMMAT + NIMFN(ISYETA)

        LWRK3 = LWORK - KEND3
        IF (LWRK3.LT.0) CALL QUIT('Out of memory in CCSDT_FA_DEN. (4)')

        CALL DZERO(WORK(KDAB), NMATAB(ISYETA))
        CALL DZERO(WORK(KDIJ), NMATIJ(ISYETA))
        CALL DZERO(WORK(KMMAT),NIMFN(ISYETA))

c       ------------------------------------------------------
c       compute triples intermediates: D(ab), D(ij), M(imfn)
c       ------------------------------------------------------
        IF ( ( (LISTL(1:3).EQ.'L0 ') .OR. (LISTL(1:3).EQ.'LE ') ) .AND. 
     *      ( (LISTR(1:3).EQ.'R1 ') .OR. (LISTR(1:3).EQ.'RE ') ) ) THEN

          IF (LISTR(1:3).EQ.'R1 ') THEN
            FREQR  = FRQLRT(IDLSTR)
            LABELR = LRTLBL(IDLSTR)
            IF (LORXLRT(IDLSTR)) CALL QUIT('NO ORBITAL RELAXED '//
     &         'PERTURBATION IMPLEMENTED IN CCSDT_FA_DEN.')
          ELSE IF (LISTR(1:3).EQ.'RE ') THEN
            FREQR  = EIGVAL(IDLSTR)
            LABELR = '- none -'
          ELSE
            CALL QUIT('Illegal right vector in CCSDT_FA_DEN.')
          END IF


          IF (FADEN_NODDY) THEN
           WRITE(LUPRI,*) 'No noddy code for the case in CCSDT_FA_DEN:'
           WRITE(LUPRI,*) 'LISTL,LISTR:',LISTL,LISTR
           WRITE(LUPRI,*) 'the real code will be used instead...'
          END IF

          !frequency of LISTL is handled inside
          CALL CCSDT_ETA_FA_DEN(LISTL,IDLSTL,ISYML,
     *                          LISTR,IDLSTR,ISYMR,FREQR,LABELR,
     *                          .FALSE.,0,
     *                          WORK(KLAMP0),WORK(KLAMH0),WORK(KFOCK0),
     *                          WORK(KT2TP),WORK(KL2TP),WORK(KL1AM),
     *                          ISYETA,WORK(KDAB),WORK(KDIJ),
     *                          .FALSE.,DUMMY,
     *                          .FALSE.,DUMMY,.TRUE.,WORK(KMMAT),
     *                          WORK(KEND3),LWRK3,
     *                          LUDELD,  FNDELD,  LUCKJD, FNCKJD,
     *                          LUDKBC,  FNDKBC,  LUTOC,  FNTOC,
     *                          LU3VI,   FN3VI,   LUDKBC3,FNDKBC3,
     *                          LU3FOPX, FN3FOPX, LU3FOP2X,FN3FOP2X,
     *                          LUR2TP, FNR2TP  )
C

        ELSE IF (  LISTL.EQ.'L0 '.AND.
     *            (LISTR.EQ.'R2 '.OR. LISTR.EQ.'ER1')       ) THEN

C

          IF (FADEN_NODDY) THEN

            CALL CCSDT_ADEN_NODDY(LISTL,IDLSTL,LISTR,IDLSTR,
     *                            WORK(KLAMP0),WORK(KLAMH0),
     *                            WORK(KFOCK0),
     *                            WORK(KDIJ),WORK(KDAB),
     *                            .FALSE.,DUMMY,.TRUE.,WORK(KMMAT),
     *                            WORK(KEND3),LWRK3,
     *                            LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,
     *                            FNDKBC,LUTOC,FNTOC,LU3VI,FN3VI,
     *                            LUDKBC3,FNDKBC3,LU3FOPX,FN3FOPX,
     *                            LU3FOP2X,FN3FOP2X)
C
          ELSE

       !FIX LISTL WITH IF STATEMENTS INSIDE !!!!
             CALL CC3_ADEN_CUB_T0('L1 ',1,LISTR,IDLSTR,
     *                     WORK(KLAMP0),WORK(KLAMH0),WORK(KFOCK0),
     *                     WORK(KDIJ),WORK(KDAB),ISYETA,
     *                     .TRUE.,WORK(KMMAT),
     *                     WORK(KEND3),LWRK3,
     *                     LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,
     *                     FNDKBC,LUTOC,FNTOC,LU3VI,FN3VI,
     *                     LUDKBC3,FNDKBC3,LU3FOPX,FN3FOPX,
     *                     LU3FOP2X,FN3FOP2X)
C
          END IF

        ELSE IF ( (LISTL.EQ.'L1 '.OR. LISTL.EQ.'LE ') .AND. 
     *            (LISTR.EQ.'R1 '.OR. LISTR.EQ.'RE ')       ) THEN

          IF (FADEN_NODDY) THEN
            CALL CCSDT_ADEN_NODDY(LISTL,IDLSTL,LISTR,IDLSTR,
     *                            WORK(KLAMP0),WORK(KLAMH0),
     *                            WORK(KFOCK0),
     *                            WORK(KDIJ),WORK(KDAB),
     *                            .FALSE.,DUMMY,.TRUE.,WORK(KMMAT),
     *                            WORK(KEND3),LWRK3,
     *                            LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,
     *                            FNDKBC,LUTOC,FNTOC,LU3VI,FN3VI,
     *                            LUDKBC3,FNDKBC3,LU3FOPX,FN3FOPX,
     *                            LU3FOP2X,FN3FOP2X)
C

          ELSE
             CALL CC3_ADEN(LISTL,IDLSTL,LISTR,IDLSTR,
     *                     WORK(KLAMP0),WORK(KLAMH0),WORK(KFOCK0),
     *                     WORK(KDIJ),WORK(KDAB),
     *                     .FALSE.,DUMMY,ISYETA,
     *                     .TRUE.,WORK(KMMAT),
     *                     WORK(KEND3),LWRK3,
     *                     LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,
     *                     FNDKBC,LUTOC,FNTOC,LU3VI,FN3VI,
     *                     LUDKBC3,FNDKBC3,LU3FOPX,FN3FOPX,
     *                     LU3FOP2X,FN3FOP2X)
C
          END IF

        ELSE
          CALL QUIT('Encountered non-implemented case in CCSDT_FA_DEN.')
        END IF

c       -------------------------------------------------
c       put intermediates on file and remember addresses:
c       -------------------------------------------------
        IADR_INTER(1,INTETA) = IADRINT
        CALL PUTWA2(LUINTER,FNINTER,WORK(KDAB),IADRINT,NMATAB(ISYETA))
        IADRINT = IADRINT + NMATAB(ISYETA)

        IADR_INTER(2,INTETA) = IADRINT
        CALL PUTWA2(LUINTER,FNINTER,WORK(KDIJ),IADRINT,NMATIJ(ISYETA))
        IADRINT = IADRINT + NMATIJ(ISYETA)

        IADR_INTER(3,INTETA) = IADRINT
        CALL PUTWA2(LUINTER,FNINTER,WORK(KMMAT),IADRINT,NIMFN(ISYETA))
        IADRINT = IADRINT + NIMFN(ISYETA)


        IF (NEED_R2TP) CALL WCLOSE2(LUR2TP,FNR2TP,'DELETE')
 
      END DO

*---------------------------------------------------------------------*
*     close/delete some file:
*---------------------------------------------------------------------*
      CALL WCLOSE2(LUDELD,FNDELD,'KEEP')
      CALL WCLOSE2(LUCKJD,FNCKJD,'KEEP')
      CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
      CALL WCLOSE2(LUDKBC,FNDKBC,'KEEP')
      CALL WCLOSE2(LU3VI,FN3VI,'KEEP')
      CALL WCLOSE2(LUDKBC3,FNDKBC3,'KEEP')
      CALL WCLOSE2(LU3FOPX,FN3FOPX,'KEEP')
      CALL WCLOSE2(LU3FOP2X,FN3FOP2X,'KEEP') 

*---------------------------------------------------------------------*
*     now compute the F{A} densities from the intermediates:
*---------------------------------------------------------------------*
      LUDEFF = -1
      CALL WOPEN2(LUDEFF,FNDEFF,64,0)

      IADRIA = 1

      DO IFADEN = 1, NFADEN

        IDLSTL = IDXL_FADEN(IFADEN)
        IDLSTB = IDXB_FADEN(IFADEN)
        IDLSTC = IDXC_FADEN(IFADEN)

        ISYML  = ILSTSYM(LISTL,IDLSTL)
        ISYMB  = ILSTSYM(LISTB,IDLSTB)
        ISYMC  = ILSTSYM(LISTC,IDLSTC)

        ISYDEN = MULD2H(MULD2H(ISYML,ISYMB),ISYMC)

        KDIA1  = KEND0
        KEND1  = KDIA1 + NT1AM(ISYDEN)

        IF (LOCDBG) THEN
          KDIA2 = KEND1
          KDIA3 = KDIA2 + NT1AM(ISYDEN)
          KDIA4 = KDIA3 + NT1AM(ISYDEN)
          KDIA5 = KDIA4 + NT1AM(ISYDEN)
          KDIA6 = KDIA5 + NT1AM(ISYDEN)
          KEND1 = KDIA6 + NT1AM(ISYDEN)
        ELSE
          KDIA2 = KDIA1
          KDIA3 = KDIA1
          KDIA4 = KDIA1
          KDIA5 = KDIA1
          KDIA6 = KDIA1
        END IF

        LWRK1 = LWORK - KEND1
        IF (LWRK1 .LT. 0) THEN
          CALL QUIT('Insufficient work space CCSDT_FA_DEN (5)')
        END IF

        CALL DZERO(WORK(KDIA1),NT1AM(ISYDEN))
        IF (LOCDBG) THEN
          CALL DZERO(WORK(KDIA2),NT1AM(ISYDEN))
          CALL DZERO(WORK(KDIA3),NT1AM(ISYDEN))
          CALL DZERO(WORK(KDIA4),NT1AM(ISYDEN))
          CALL DZERO(WORK(KDIA5),NT1AM(ISYDEN))
          CALL DZERO(WORK(KDIA6),NT1AM(ISYDEN))
        END IF

        CALL CCSDT_FA_DEN1(WORK(KDIA1),WORK(KDIA2),WORK(KDIA3),
     &                     IDLSTL,ISYML,
     &                     LISTB,IDLSTB,ISYMB,
     &                     LISTC,IDLSTC,ISYMC,
     &                     IDXL_INTER,IDXR_INTER,LSTR_INTER,
     &                     IADR_INTER,LUINTER,FNINTER,
     &                     NIMFN,NINTETA,WORK(KEND1),LWRK1)

        CALL CCSDT_FA_DEN1(WORK(KDIA4),WORK(KDIA5),WORK(KDIA6),
     &                     IDLSTL,ISYML,
     &                     LISTC,IDLSTC,ISYMC,
     &                     LISTB,IDLSTB,ISYMB,
     &                     IDXL_INTER,IDXR_INTER,LSTR_INTER,
     &                     IADR_INTER,LUINTER,FNINTER,
     &                     NIMFN,NINTETA,WORK(KEND1),LWRK1)

        IF (LOCDBG) THEN
          WRITE(LUPRI,*) 'triples contribution 1 to F{A} density:',
     &      DDOT(NT1AM(ISYDEN),WORK(KDIA1),1,WORK(KDIA1),1)
          CALL CC_PRP(WORK(KDIA1),DUMMY,1,1,0)
          WRITE(LUPRI,*) 'triples contribution 2 to F{A} density:',
     &      DDOT(NT1AM(ISYDEN),WORK(KDIA2),1,WORK(KDIA2),1)
          CALL CC_PRP(WORK(KDIA2),DUMMY,1,1,0)
          WRITE(LUPRI,*) 'triples contribution 3 to F{A} density:',
     &      DDOT(NT1AM(ISYDEN),WORK(KDIA3),1,WORK(KDIA3),1)
          CALL CC_PRP(WORK(KDIA3),DUMMY,1,1,0)
          WRITE(LUPRI,*) 'triples contribution 4 to F{A} density:',
     &      DDOT(NT1AM(ISYDEN),WORK(KDIA4),1,WORK(KDIA4),1)
          CALL CC_PRP(WORK(KDIA4),DUMMY,1,1,0)
          WRITE(LUPRI,*) 'triples contribution 5 to F{A} density:',
     &      DDOT(NT1AM(ISYDEN),WORK(KDIA5),1,WORK(KDIA5),1)
          CALL CC_PRP(WORK(KDIA5),DUMMY,1,1,0)
          WRITE(LUPRI,*) 'triples contribution 6 to F{A} density:',
     &      DDOT(NT1AM(ISYDEN),WORK(KDIA6),1,WORK(KDIA6),1)
          CALL CC_PRP(WORK(KDIA6),DUMMY,1,1,0)

          ! sum up all 6 contributions
          CALL DAXPY(NT1AM(ISYDEN),1.0D0,WORK(KDIA2),1,WORK(KDIA1),1)
          CALL DAXPY(NT1AM(ISYDEN),1.0D0,WORK(KDIA3),1,WORK(KDIA1),1)
          CALL DAXPY(NT1AM(ISYDEN),1.0D0,WORK(KDIA4),1,WORK(KDIA1),1)
          CALL DAXPY(NT1AM(ISYDEN),1.0D0,WORK(KDIA5),1,WORK(KDIA1),1)
          CALL DAXPY(NT1AM(ISYDEN),1.0D0,WORK(KDIA6),1,WORK(KDIA1),1)

          WRITE(LUPRI,*) 'complete F{A} density:',
     &      DDOT(NT1AM(ISYDEN),WORK(KDIA1),1,WORK(KDIA1),1)
          CALL CC_PRP(WORK(KDIA1),DUMMY,1,1,0)
        END IF
          
        ! store density on file
        IADR_DEN(IFADEN) = IADRIA
        CALL PUTWA2(LUDEFF,FNDEFF,WORK(KDIA1),IADRIA,NT1AM(ISYDEN))
        IADRIA = IADRIA + NT1AM(ISYDEN)

      END DO

      CALL WCLOSE2(LUINTER,FNINTER,'DELETE')

*---------------------------------------------------------------------*
*     now compute from the densities the contractions F{A} t^B t^C:
*---------------------------------------------------------------------*
      DO ITRAN = 1, NFATRAN

        IDLSTL = IFATRAN(1,ITRAN)
        IOPERA = IFATRAN(2,ITRAN)
        IDLSTB = IFATRAN(3,ITRAN)
        IFILE  = IFATRAN(4,ITRAN)
        IRELAX = IFATRAN(5,ITRAN)

        LABELA = LBLOPR(IOPERA)
        LPDBSA = LPDBSOP(IOPERA)

        ISYML = ILSTSYM(LISTL,IDLSTL)
        ISYMA = ILSTSYM(LISTO,IOPERA)
        ISYMB = ILSTSYM(LISTB,IDLSTB)

        ISYRES = MULD2H(MULD2H(ISYMA,ISYMB),ISYML)

        IF ( IRELAX.GT.0 .OR. LPDBSA ) THEN
          CALL QUIT('Relaxed perturbations not yet implemented in '//
     &              'CCSDT_FA_DEN.')
        END IF

        IVEC = 1
        DO WHILE(IFADOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC) 
          IDLSTC = IFADOTS(IVEC,ITRAN)
          ISYMC  = ILSTSYM(LISTC,IDLSTC)
 
          IF (ISYMC.NE.ISYRES) THEN
            CALL QUIT('Symmetry mismatch in CCSDT_FA_DEN')
          END IF

          ! find index of the density needed
          IFADEN = 0
          DO IDX = 1, NFADEN
            IF (IDLSTL .EQ.IDXL_FADEN(IDX) .AND.
     &          IDLSTB .EQ.IDXB_FADEN(IDX) .AND.
     &          IDLSTC .EQ.IDXC_FADEN(IDX)       ) THEN
              IFADEN = IDX
            END IF
            IF (IFADEN.EQ.0 .AND. LISTB.EQ.LISTC .AND.
     &          IDLSTL .EQ.IDXL_FADEN(IDX) .AND.
     &          IDLSTC .EQ.IDXB_FADEN(IDX) .AND.
     &          IDLSTB .EQ.IDXC_FADEN(IDX)       ) THEN
              IFADEN = IDX
            END IF
          END DO
          IF (IFADEN.LE.0) CALL QUIT('Fatal error in CCSDT_FA_DEN.')

          ISYDEN = MULD2H(MULD2H(ISYML,ISYMB),ISYMC)

          IF (ISYDEN.EQ.ISYMA) THEN
            KDIA1   = KEND0
            KFOCK   = KDIA1   + NT1AM(ISYDEN)
            KFOCKIA = KFOCK   + N2BST(ISYDEN)
            KEND1   = KFOCKIA + NT1AM(ISYDEN)
            LWRK1   = LWORK - KEND1

            IADRIA = IADR_DEN(IFADEN)
            CALL GETWA2(LUDEFF,FNDEFF,WORK(KDIA1),IADRIA,NT1AM(ISYDEN))

            CALL CCPRPAO(LABELA,.TRUE.,WORK(KFOCK),IRREP,IMAT,IERR,
     *                   WORK(KEND1),LWRK1)
            IF ((IERR.GT.0) .OR. (IERR.EQ.0 .AND. IRREP.NE.ISYMA)) THEN
              WRITE(LUPRI,*) 'ISYMA :',ISYMA
              WRITE(LUPRI,*) 'IRREP :',IRREP
              WRITE(LUPRI,*) 'IERR  :',IERR
              WRITE(LUPRI,*) 'LABEL :',LABELA
              CALL QUIT('CCSDT_FA_DEN: error reading operator ')
            ELSE IF (IERR.LT.0) THEN
              CALL DZERO(WORK(KFOCK),N2BST(ISYMA))
            END IF

            ! transform property integrals to Lambda-MO basis and resort
            CALL CC_FCKMO(WORK(KFOCK),WORK(KLAMP0),WORK(KLAMH0),
     &                    WORK(KEND1),LWRK1,ISYMA,1,1) 
            CALL CC_FOCK_RESORT(DUMMY,.FALSE.,WORK(KFOCKIA),.TRUE.,
     &                          DUMMY,.FALSE.,DUMMY,.FALSE.,
     &                          WORK(KFOCK),ISYMA)

            TRIPCON = DDOT(NT1AM(ISYMA),WORK(KDIA1),1,WORK(KFOCKIA),1)

            FACON(IVEC,ITRAN) = FACON(IVEC,ITRAN) + TRIPCON
          END IF

          IVEC = IVEC + 1
        END DO

      END DO

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_FA_DEN                         *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_FA_DEN1(DIA1,DIA2,DIA3,IDLSTL,ISYML,
     &                         LISTB,IDLSTB,ISYMB,
     &                         LISTC,IDLSTC,ISYMC,
     &                         IDXL_INTER,IDXR_INTER,LSTR_INTER,
     &                         IADR_INTER,LUINTER,FNINTER,
     &                         NIMFN,NINTER,WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: compute contribution to F{A} density from
*              precomputed intermediates on file
*
*     Written by Christof Haettig, Mai 2003, Friedrichstal
*
*=====================================================================*
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dummy.h>

      CHARACTER FNINTER*(*)
      CHARACTER*3 LSTR_INTER(*), LISTB, LISTC

* input/output:
      INTEGER LWORK, NINTER, NIMFN(8)
      INTEGER IDLSTL, ISYML, IDLSTB, ISYMB, IDLSTC, ISYMC, LUINTER
      INTEGER IDXL_INTER(*), IDXR_INTER(*), IADR_INTER(3,*)

#if defined (SYS_CRAY)
      REAL DIA1(*), DIA2(*), DIA3(*)
      REAL WORK(*)
      REAL ONE, TWO
#else
      DOUBLE PRECISION DIA1(*), DIA2(*), DIA3(*)
      DOUBLE PRECISION WORK(*)
      DOUBLE PRECISION ONE, TWO
#endif
      PARAMETER ( ONE=1.0D0, TWO=2.0D0 )

* local:
      CHARACTER CDUMMY*1, MODEL*10
      INTEGER ISYETA, ISYDEN, KDAE, KDIJ, KMMAT, KC1AM, KC2TP, KEND1,
     &        LWRK1, IDXLB, IDX, IADRINT, IOPT, IOPTT2, ISYMA, ISYMI,
     &        ISYME, ISYMJ, KOFFDEA, KOFFDIJ, KOFFDIA, KOFFTEI,
     &        KOFFTAJ, NVIRA, NVIRE, NRHFI

   

c     ------------------
c     memory allocation:
c     ------------------
      ISYETA = MULD2H(ISYML, ISYMB) 
      ISYDEN = MULD2H(ISYETA,ISYMC) 

      KDAE  = 1
      KDIJ  = KDAE  + NMATAB(ISYETA)
      KMMAT = KDIJ  + NMATIJ(ISYETA)
      KC1AM = KMMAT + NIMFN(ISYETA)
      KC2TP = KC1AM + NT1AM(ISYMC)
      KEND1 = KC2TP + NT2SQ(ISYMC)

      LWRK1 = LWORK - KEND1
      IF (LWRK1.LT.NT2SQ(ISYMC)) THEN
        CALL QUIT('Out of memory in CCSDT_FA_DEN1.')
      END IF

c     -----------------------------------------------------
c     restore intermediates depending on L and B from file
c     -----------------------------------------------------
      ! get index for the intermediates depending on L and B:
      IDXLB = 0
      DO IDX = 1, NINTER
        IF (IDLSTL.EQ.IDXL_INTER(IDX) .AND.
     &      IDLSTB.EQ.IDXR_INTER(IDX) .AND.
     &      LISTB .EQ.LSTR_INTER(IDX)      ) THEN
          IDXLB = IDX
        END IF
      END DO
      IF (IDXLB.LE.0) CALL QUIT('Fatal error in CCSDT_FA_DEN1.')

      ! read intermediates
      IADRINT = IADR_INTER(1,IDXLB)
      CALL GETWA2(LUINTER,FNINTER,WORK(KDAE),IADRINT,NMATAB(ISYETA))

      IADRINT = IADR_INTER(2,IDXLB) 
      CALL GETWA2(LUINTER,FNINTER,WORK(KDIJ),IADRINT,NMATIJ(ISYETA))

      IADRINT = IADR_INTER(3,IDXLB) 
      CALL GETWA2(LUINTER,FNINTER,WORK(KMMAT),IADRINT,NIMFN(ISYETA))

c     --------------------------------------------
c     read C vector and resort doubles amplitudes:
c     --------------------------------------------
      IOPT = 3
      CALL CC_RDRSP(LISTC,IDLSTC,ISYMC,IOPT,MODEL,
     &              WORK(KC1AM),WORK(KC2TP))
      CALL CCLR_DIASCL(WORK(KC2TP),TWO,ISYMC)
      CALL CC_T2SQ(WORK(KC2TP),WORK(KEND1),ISYMC)
      CALL CC3_T2TP(WORK(KC2TP),WORK(KEND1),ISYMC)


C     ---------------------------------------------
C     D(ia) <-- D(ia) + sum_fnm t(am,fn) y^M(imfn):
C     ---------------------------------------------
C     ! turn sign for this contribution to D(ia)
C     CALL DSCAL(NT2SQ(ISYMC),-1.0D0,WORK(KC2TP),1)

      IOPTT2 = 0
      CALL CCSDT_ETA_TM2(DIA3,ISYDEN,WORK(KMMAT),ISYETA,
     &                   WORK(KC2TP),IDUMMY,CDUMMY,ISYMC,IOPTT2,
     &                   WORK(KEND1),LWRK1)

 
C     ----------------------------------------------------------------
C     D(ia) <-- D(ia) - sum_e y^t(ei) D^0(ea) + sum_j y^t(aj) D^0(ij):
C     ----------------------------------------------------------------
      DO ISYMA = 1, NSYM
        ISYMI = MULD2H(ISYDEN,ISYMA)
        ISYME = MULD2H(ISYETA,ISYMA)
        ISYMJ = MULD2H(ISYETA,ISYMI)
        
        KOFFDEA = KDAE  + IMATAB(ISYME,ISYMA)
        KOFFDIJ = KDIJ  + IMATIJ(ISYMI,ISYMJ)
        KOFFDIA = 1     + IT1AM(ISYMA,ISYMI)
        KOFFTEI = KC1AM + IT1AM(ISYME,ISYMI)
        KOFFTAJ = KC1AM + IT1AM(ISYMA,ISYMJ)
        
        NVIRA   = MAX(NVIR(ISYMA),1)
        NVIRE   = MAX(NVIR(ISYME),1)
        NRHFI   = MAX(NRHF(ISYMI),1)

        CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYME),
     &             -ONE,WORK(KOFFDEA),NVIRE,WORK(KOFFTEI),NVIRE,
     &              ONE,DIA1(KOFFDIA),NVIRA)

        CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &              ONE,WORK(KOFFTAJ),NVIRA,WORK(KOFFDIJ),NRHFI,
     &              ONE,DIA2(KOFFDIA),NVIRA)

      END DO
                    
      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_FA_DEN1                        *
*---------------------------------------------------------------------*
! -- end of cc_famat.F --
