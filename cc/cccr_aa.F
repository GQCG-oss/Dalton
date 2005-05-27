C
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren.
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
c /* deck cccr_aa */
*=====================================================================*
       SUBROUTINE CCCR_AA ( LABELA, ISYMTA,  ! inp: label/symmetry A
     &                      LISTB,  ITAMPB,  ! inp: B resp. amplit.
     &                      XINT,
     &                      WORK,   LWORK   )! work space
*---------------------------------------------------------------------*
*
*    Purpose: transformation of a response vector with a Jacobian
*             where the hamiltonian has been substituted by a 
*             perturbation operator
*
*             A{A} * T^B = <mu|[A,T^B]|CC>
*
*    Note: the single and double excitation parts of the result RHO 
*          are returned at the beginning of the work space in
*          WORK(1)... WORK(NT1AM(ISYRES))
*          WORK(NT1AM(ISYRES)+1)... WORK(NT1AM(ISYRES)+NT2AM(ISYRES))
*          (double excitation part will be stored in packed form)
*
*    symmetries:
*           ISYMTA : A perturbation
*           ISYMTB : B response amplitudes 
*           ISYRES : result vector RHO1, RHO2, i.e. ISYMTA x ISYMTB
*
*     Written by Christof Haettig, Februar 1997.
*
*             JK+OC. CCSLV: Allow for input of integrals if
*             LABELA .eq. 'GIVE INT'
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
C#include <ccr1rsp.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER ISYM0
      PARAMETER (ISYM0 = 1) ! symmetry of the reference state
      INTEGER KDUM
      PARAMETER (KDUM =  99 999 999) ! dummy address in work space

#if defined (SYS_CRAY)
      REAL ONE, TWO, THREE, HALF
#else
      DOUBLE PRECISION ONE, TWO, THREE, HALF
#endif
      PARAMETER (ONE = 1.0d0, TWO = 2.0d0, THREE = 3.0d0, HALF = 0.5d0)

      CHARACTER*8 LABELA
      CHARACTER*10 MODEL
      CHARACTER LISTB*(*)
      INTEGER ISYRES, ISYMTA, ISYMTB, ITAMPB, LWORK, IRREP,ISYM,IERR
      INTEGER KRHO1, KRHO2, KEND0, KAOO, KAOV, KAVV, KBTAOO, KBTAVV
      INTEGER KT1AMPB, KT2AMPB, KEND1, LEND1, KEND2, LEND2, IOPT
      INTEGER KPERTA, KT1AMP0, KT2AMP0, KLAMDP0, KLAMDH0
      


#if defined (SYS_CRAY)
      REAL FREQB
      REAL WORK(LWORK), XINT(*)
#else
      DOUBLE PRECISION FREQB
      DOUBLE PRECISION WORK(LWORK), XINT(*)
#endif

* external functions:
      INTEGER ILSTSYM

  
*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF ( .not. (CCS .or. CC2 .or. CCSD .or. CC3) ) THEN
        WRITE(LUPRI,'(/1x,a)') 'CCCR_AA called for a Coupled Cluster '
     &          //'method not implemented in CCCR_AA...'
        CALL QUIT('Unkown CC method in CCCR_AA.')
      END IF

*---------------------------------------------------------------------*
* set & check symmetries:
*---------------------------------------------------------------------*
      ISYMTB  = ILSTSYM(LISTB,ITAMPB)  ! B
      ISYRES  = MULD2H(ISYMTA,ISYMTB)  ! A x B

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'LISTB,ITAMPB,ISYMTB:',LISTB,ITAMPB,ISYMTB
        WRITE (LUPRI,*) 'LABELA,ISYMTA:',LABELA,ISYMTA
      END IF

*---------------------------------------------------------------------*
* flush print unit
*---------------------------------------------------------------------*
      Call FLSHFO(LUPRI)

      IF (LOCDBG) THEN
        WRITE(LUPRI,'(/1x,a,i15)') 'work space in CCCR_AA:',LWORK
      END IF
*---------------------------------------------------------------------*
* initialize pointer for work space and allocate memory for
*  1) single excitation part of the result vector 
*  2) double excitation part of the result vector 
*  3) AO/MO perturbation integrals A (complete matrix)
*  4) MO transformed perturbation integrals A (occ/occ block)
*  5) MO transformed perturbation integrals A (vir/vir block)
*  6) MO transformed perturbation integrals A (occ/vir block)
*  7) one-index transformed perturbation integrals A^B (occ/occ block)
*  8) one-index transformed perturbation integrals A^B (vir/vir block)
*  9) singles part of response amplitudes T1^B
*---------------------------------------------------------------------*
      KRHO1   = 1
      KRHO2   = KRHO1   + NT1AM(ISYRES)
      KEND0   = KRHO2   + NT2AM(ISYRES)

      IF (CCS) THEN ! no double excitation result vector for CCS
        KEND0 = KRHO2
        KRHO2 = KDUM
      END IF

      KPERTA  = KEND0
      KAOO    = KPERTA  + N2BST(ISYMTA)
      KAOV    = KAOO    + NMATIJ(ISYMTA)
      KAVV    = KAOV    + NT1AM(ISYMTA)
      KBTAOO  = KAVV    + NMATAB(ISYMTA)
      KBTAVV  = KBTAOO  + NMATIJ(ISYRES)
      KT1AMPB = KBTAVV  + NMATAB(ISYRES)
      KEND1   = KT1AMPB + NT1AM(ISYMTB)
      LEND1   = LWORK - KEND1

      IF (LEND1 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CCCR_AA.')
      END IF

*---------------------------------------------------------------------*
* initialize single excitation part of result vector RHO1:
*---------------------------------------------------------------------*
      Call DZERO (WORK(KRHO1), NT1AM(ISYRES))

*---------------------------------------------------------------------*
* read singles parts for B response amplitudes:
*---------------------------------------------------------------------*
      IOPT = 1
      CALL CC_RDRSP(LISTB,ITAMPB,ISYMTB,IOPT,MODEL,
     &                  WORK(KT1AMPB),WORK(KDUM)  )

      IF (LOCDBG) THEN
        CAll AROUND('response T amplitudes B:')
        WRITE (LUPRI,*) 'LIST/INDEX:',LISTB,ITAMPB
        WRITE (LUPRI,*) 'Symmetry:      ',ISYMTB
        CAll CC_PRP(WORK(KT1AMPB),WORK(KDUM),ISYMTB,1,0)
      END IF

*---------------------------------------------------------------------*
* read & resort one-electron integrals for operator A:
*---------------------------------------------------------------------*
      KT1AMP0 = KEND1
      KLAMDP0 = KT1AMP0 + NT1AM(ISYM0)
      KLAMDH0 = KLAMDP0 + NLAMDT
      KEND2   = KLAMDH0 + NLAMDT
      LEND2   = LWORK - KEND2

      IF (LEND2 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CCCR_AA.')
      END IF
C
C     JK+OC, CCSLV 
      IF (LABELA.EQ.'GIVE INT') THEN
        CALL DCOPY(N2BST(ISYMTA),XINT(1),1,WORK(KPERTA),1)
      ELSE
* read the AO integrals:
        CALL CCPRPAO(LABELA,WORK(KPERTA),IRREP,ISYM,IERR,
     &               WORK(KEND2),LEND2)
        IF (IERR.NE.0 .OR. IRREP.NE.ISYMTA) THEN
          CALL QUIT('CCCR_AA: error while reading operator '//LABELA)
        END IF
      END IF
* get single excitation part of the zeroth order amplitudes:
      IOPT = 1
      CALL CC_RDRSP('R0',1,ISYM0,IOPT,MODEL,WORK(KT1AMP0),WORK(KDUM))

* get zeroth order Lambda matrices:
      CALL LAMMAT(WORK(KLAMDP0),WORK(KLAMDH0),WORK(KT1AMP0),
     &            WORK(KEND2),LEND2)

* transform one-electron integrals in place:
      CALL CC_FCKMO(WORK(KPERTA),WORK(KLAMDP0),WORK(KLAMDH0),
     &              WORK(KEND2),LEND2,ISYMTA,1,1)

* gather occ/occ, occ/vir, and vir/vir blocks:
      CALL CC_GATHEROO(WORK(KPERTA),WORK(KAOO),ISYMTA)
      CALL CC_GATHEROV(WORK(KPERTA),WORK(KAOV),ISYMTA)
      CALL CC_GATHERVV(WORK(KPERTA),WORK(KAVV),ISYMTA)

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CCCR_AA> A integrals in '//
     &        'MO basis (occ/vir):'
        WRITE (LUPRI,*) 'DEBUG_CCCR_AA> label, symmetry:',LABELA,ISYMTA
        Call CC_PRP(WORK(KAOV),WORK(KDUM),ISYMTA,1,0)
        WRITE (LUPRI,*)'DEBUG_CCCR_AA> A integrals (occ/occ) '//
     &       'and (vir/vir):'
        CALL CC_PREI(WORK(KAVV),WORK(KAOO),ISYMTA,1)
      END IF

*---------------------------------------------------------------------*
* calculate A perturbation integrals one-index transformed with
* the B response amplitudes T1^B:
*---------------------------------------------------------------------*
      IF (.NOT.(CCS .OR. CCSTST)) THEN

*       occ/occ block:
        Call CCG_1ITROO(WORK(KBTAOO), ISYRES,
     &                  WORK(KAOV),   ISYMTA,
     &                  WORK(KT1AMPB),ISYMTB  )

*       vir/vir block:
        Call CCG_1ITRVV(WORK(KBTAVV), ISYRES,
     &                  WORK(KAOV),   ISYMTA,
     &                  WORK(KT1AMPB),ISYMTB  )

      END IF

*=====================================================================*
*   CCS part:  < mu_1 | [A, T1^B] | HF>
*=====================================================================*
      IOPT  = 1
      Call CCG_1ITRVO(WORK(KRHO1),ISYRES,WORK(KAOO),WORK(KAVV),
     &                ISYMTA,WORK(KT1AMPB),ISYMTB,IOPT          )

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CCCR_AA> CCS contribution to RHO1:'
        Call CC_PRP(WORK(KRHO1),WORK(KDUM),ISYRES,1,0)
      END IF

*---------------------------------------------------------------------*
* end of CCS part
*---------------------------------------------------------------------*
      If (CCS) Return

*---------------------------------------------------------------------*
* CC2/CCSD part
*---------------------------------------------------------------------*

* initialize RHO2:
      CALL DZERO(WORK(KRHO2),NT2AM(ISYRES))

      If (CCSTST) Return

*---------------------------------------------------------------------*
* E term like contribution to doubles part: <mu_2|[A,T2^B]|HF>
*---------------------------------------------------------------------*
       KT2AMPB = KEND1
       KEND2   = KT2AMPB + NT2SQ(ISYMTB)
       LEND2   = LWORK - KEND2

       IF (LEND2 .LT. NT2AM(ISYMTB) ) THEN
         CALL QUIT('Insufficient work space in CCCR_AA.')
       END IF

* read amplitudes in scratch space, scale diagonal and square up:
      IOPT = 2
      CALL CC_RDRSP(LISTB,ITAMPB,ISYMTB,IOPT,MODEL,
     &                  WORK(KDUM),WORK(KEND2)  )

      CAll CCLR_DIASCL(WORK(KEND2),TWO,ISYMTB)

      CALL CC_T2SQ(WORK(KEND2), WORK(KT2AMPB), ISYMTB)

* do the caculation:
      Call CCRHS_E(WORK(KRHO2),WORK(KT2AMPB),WORK(KAVV), 
     &             WORK(KAOO), WORK(KEND2), LEND2, ISYMTB, ISYMTA)

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CCCR_AA> RHO after first E contribution:'
        Call CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,1,1)
      END IF

*---------------------------------------------------------------------*
* I term like contribution to singles part: <mu_1| [A,T2^B] |HF>
*---------------------------------------------------------------------*

* calculate 2*t^B(iajb) - t^B(ibja) in place
      CALL CCRHS_T2TR(WORK(KT2AMPB),WORK(KEND2),LEND2,ISYMTB)

* calculate the I term like contribution:
      CALL CCRHS_I(WORK(KRHO1),WORK(KT2AMPB),WORK(KPERTA),
     &             WORK(KEND2), LEND2, ISYMTB, ISYMTA)

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CCCR_AA> RHO1 after '//
     &        'CC2/CCSD contribution:'
        Call CC_PRP(WORK(KRHO1),WORK(KDUM),ISYRES,1,0)
      END IF

*---------------------------------------------------------------------*
* calculate second E term like contribution <mu_2|[A^B,T2^0]|HF>
*---------------------------------------------------------------------*
       KT2AMP0 = KEND1
       KEND2   = KT2AMP0 + NT2SQ(ISYM0)
       LEND2   = LWORK - KEND2

       IF (LEND2 .LT. NT2AM(ISYM0) ) THEN
         CALL QUIT('Insufficient work space in CCCR_AA.')
       END IF

* read amplitudes in scratch space and square up:
      IOPT = 2
      CALL CC_RDRSP('R0',1,ISYM0,IOPT,MODEL,WORK(KDUM),WORK(KEND2))

      CALL CC_T2SQ(WORK(KEND2), WORK(KT2AMP0), ISYM0)

* do the caculation:
      Call CCRHS_E(WORK(KRHO2),WORK(KT2AMP0),WORK(KBTAVV), 
     &             WORK(KBTAOO), WORK(KEND2), LEND2, ISYM0, ISYRES)

*---------------------------------------------------------------------*
      CAll CCLR_DIASCL(WORK(KRHO2),HALF,ISYRES)

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'Final result of CCCR_AA: RHO:'
        Call CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,1,1)
      END IF

      RETURN
      END
*=====================================================================*
*                  END OF SUBROUTINE CCCR_AA                          *
*=====================================================================*
