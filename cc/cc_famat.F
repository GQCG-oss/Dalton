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
#include <single.h>
C
c /* deck cc_famata */
*=====================================================================*
       SUBROUTINE CC_FAMAT( LABELA, ISYMTA,  ! inp: label/symmetry A
     &                      LISTB,  ITAMPB,  ! inp: B resp. amplit.
     &                      LISTC,  IZETVC,  ! inp: C resp. zeta vec.
     &                      IOPTRES,IFATRAN, ! output option
     &                      IFADOTS,FACON,   ! indeces/dotproducts
     &                      FILFA,  ITRAN,   ! list,index for dotprods
     &                      NFATRAN,MXVEC,   ! dimensions for dotprods
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
C#include <ccr1rsp.h>

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
      INTEGER IFATRAN(5,NFATRAN), IFADOTS(MXVEC,NFATRAN)

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

  
      CALL QENTER('CC_FAMAT')
*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF ( .not. (CCS .or. CC2 .or. CCSD .or. CC3) ) THEN
        WRITE (LUPRI,'(/1x,a)') 'CC_FAMAT called for a Coupled Cluster '
     &          //'method not implemented in CC_FAMAT...'
        CALL QUIT('Unkown CC method in CC_FAMAT.')
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
      CALL CCPRPAO(LABELA,WORK(KCTMO),IRREP,ISYM,IERR,
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

        CALL CCSDT_FAMAT_NODDY(LISTC,IZETVC,LISTB,ITAMPB,IOPTRES,
     &                         LABELA,FREQA,
     &                         WORK(KGAMMA1),WORK(KGAMMA2),
     &                         WORK(KGAMMA1EFF),WORK(KGAMMA2EFF),
     &                         IFADOTS,FACON,FILFA,ITRAN,
     &                         NFATRAN,MXVEC,WORK(KEND0),LEND0 )
      END IF

*=====================================================================*
* final section: depending on IOPTRES store vector in memory or on
*                file or contract it with some other vectors:
*
*      the memory has to be organized as follows:
*         kgamma1     singles result vector
*         kgamma2     doubles result vector
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
         CALL QUIT('Unkown coupled cluster model in CC_FAMAT.')
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
     &        KFOCKC, LWRK1, KRECNR, KE3AM, KTC3AM, KTC1AM, KTC2AM,
     &        ISYMA, KINT1SB, KINT2SB, KLAMPB, KLAMHB, KFOCKB, ISYMB,
     &        KFOCKA, KFOCKA_AO, KFOCK0, IRREP, ISYM, IERR, KFOCKAB1,
     &        KFIELD, KFIELDAO

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J 

      CALL QENTER('CCSDT_FAMAT_NODDY')

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
      LUFOCK = 0
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

      CALL CCPRPAO(LABELA,WORK(KFOCKA_AO),IRREP,ISYM,IERR,
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
     *                  WORK(KRECNR),DIRECT)
 
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
*     Compute first-order response amplitudes:
*---------------------------------------------------------------------*
      KINT1SB = KEND3
      KINT2SB = KINT1SB + NT1AMX*NVIRT*NVIRT
      KEND3   = KINT2SB + NRHFT*NRHFT*NT1AMX

      KLAMPB  = KEND3
      KLAMHB  = KLAMPB + NLAMDT
      KFOCKB  = KLAMHB + NLAMDT
      KEND3   = KFOCKB + NORBT*NORBT

      LWRK3  = LWORK  - KEND3
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_FAMAT_NODDY')
      ENDIF

      CALL CCSDT_T31_NODDY(WORK(KT3AM),LISTB,IDLSTB,FREQB,.FALSE.,
     &                     .FALSE.,WORK(KINT1S0),WORK(KINT2S0),
     &                     .FALSE.,WORK(KDUM),WORK(KDUM),
     &                     .FALSE.,WORK(KDUM),WORK(KDUM),
     &                             WORK(KINT1SB),WORK(KINT2SB),
     &                     WORK(KLAMPB),WORK(KLAMHB),WORK(KFOCKB),
     &                     WORK(KLAMP0),WORK(KLAMH0),WORK(KFOCK0),
     &                     WORK(KDUM),WORK(KFOCKD),
     &                     WORK(KEND3),LWRK3)

*---------------------------------------------------------------------*
*     Compute L^0_3 multipliers:
*---------------------------------------------------------------------*
      IF (LISTL(1:3).EQ.'L0 ') THEN

        FREQL = 0.0D0

        CALL DZERO(WORK(KL3AM),NT1AMX*NT1AMX*NT1AMX)

        ! remember that CCSDT_L03AM returns -L3 !!
        CALL CCSDT_L03AM(WORK(KL3AM),WORK(KINT1T0),WORK(KINT2T0),
     *                   WORK(KXIAJB),WORK(KFOCK0),WORK(KL1AM),
     *                   WORK(KL2AM),WORK(KSCR1),WORK(KFOCKD),
     *                   DUMMY,DUMMY)

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
