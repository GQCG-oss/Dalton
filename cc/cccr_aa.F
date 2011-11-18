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
*    Note: the single, double and r12-double excitation parts of the 
*          result RHO are returned at the beginning of the work space in
*          WORK(1)... WORK(NT1AM(ISYRES))
*          WORK(NT1AM(ISYRES)+1)... WORK(NT1AM(ISYRES)+NT2AM(ISYRES))
*          WORK(WORK(NT1AM(ISYRES)+NT2AM(ISYRES)+1)... WORK(NT1AM(ISYRES)+NT2AM(ISYRES)+NTR12AM(ISYRES))
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
#include <r12int.h>
#include <ccr12int.h>

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
      INTEGER LUNIT,IDUM,IAN,KLAMDPB,KLAMDHB,KRHO12SQ,KRHOR12
      INTEGER KXINTTRI,KXINTSQ,KVXINTSQ,KT12AM0,KT12AMB,KSCR
      

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
        WRITE(LUPRI,'(/a)') ' CCCR_AA called for a Coupled Cluster '
     &          //'method not implemented in CCCR_AA...'
        CALL QUIT('Unknown CC method in CCCR_AA.')
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

      IF (CCR12) THEN
        KRHOR12 = KEND0 
        KEND0   = KRHOR12 + NTR12AM(ISYRES)
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
      KEND1   = KLAMDH0 + NLAMDT
      LEND1   = LWORK - KEND1

      IF (LEND1 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CCCR_AA.')
      END IF
C
C     JK+OC, CCSLV 
      IF (LABELA.EQ.'GIVE INT') THEN
        CALL DCOPY(N2BST(ISYMTA),XINT(1),1,WORK(KPERTA),1)
      ELSE
* read the AO integrals:
        CALL CCPRPAO(LABELA,.TRUE.,WORK(KPERTA),IRREP,ISYM,IERR,
     &               WORK(KEND1),LEND1)
        IF (IERR.NE.0 .OR. IRREP.NE.ISYMTA) THEN
          CALL QUIT('CCCR_AA: error while reading operator '//LABELA)
        END IF
      END IF
* get single excitation part of the zeroth order amplitudes:
      IOPT = 1
      CALL CC_RDRSP('R0',1,ISYM0,IOPT,MODEL,WORK(KT1AMP0),WORK(KDUM))

* get zeroth order Lambda matrices:
      CALL LAMMAT(WORK(KLAMDP0),WORK(KLAMDH0),WORK(KT1AMP0),
     &            WORK(KEND1),LEND1)

* transform one-electron integrals in place:
      CALL CC_FCKMO(WORK(KPERTA),WORK(KLAMDP0),WORK(KLAMDH0),
     &              WORK(KEND1),LEND1,ISYMTA,1,1)

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

*---------------------------------------------------------------------*
* calculate R12 contributions:
*
* C. Neiss,  june 2005
*---------------------------------------------------------------------*
      IF (CCR12) THEN
        CALL DZERO(WORK(KRHOR12),NTR12AM(ISYRES))
C
        KT12AM0 = KEND1
        KT12AMB = KT12AM0 + NTR12SQ(1)
        KXINTTRI= KT12AMB + NTR12SQ(ISYMTB)
        KXINTSQ = KXINTTRI+ NR12R12P(1)
        KVXINTSQ= KXINTSQ + NR12R12SQ(1)
        KSCR    = KVXINTSQ+ NR12R12SQ(ISYMTA)
        KRHO12SQ= KSCR    + NTR12SQ(ISYRES)
        KLAMDPB = KRHO12SQ+ NTR12SQ(ISYRES)
        KLAMDHB = KLAMDPB + NGLMDT(ISYMTB)
        KEND2   = KLAMDHB + NGLMDT(ISYMTB)
        LEND2   = LWORK - KEND2
        IF (LEND2.LT.0) THEN
          CALL QUIT('Insufficient work space for R12 in CCCR_AA')
        END IF

* transform MO matrices with singles response amplitudes:
        CALL CCLR_LAMTRA(WORK(KLAMDP0),WORK(KLAMDPB),WORK(KLAMDH0),
     &                   WORK(KLAMDHB),WORK(KT1AMPB),ISYMTB)

* read R12 ground state amplitudes:
        CALL CC_R12GETCT(WORK(KT12AM0),1,0,IDUM,.FALSE.,'N',
     &                   IDUM,IDUM,IDUM,IDUM,IDUM,WORK(KEND2),LEND2)

* read R12 response amplitudes:
        CALL CC_R12GETCT(WORK(KT12AMB),ISYMTB,2,KETSCL,.FALSE.,'N',
     &                   IDUM,IDUM,IDUM,LISTB,ITAMPB,WORK(KEND2),LEND2)

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

* read R12 VXINT:
        CALL DZERO(WORK(KVXINTSQ),NR12R12SQ(ISYMTA))
        CALL CC_R12RDVXINT(WORK(KVXINTSQ),WORK(KEND2),LEND2,ONE,
     &                     ISYMTA,LABELA)

* read the AO integrals:
      CALL CCPRPAO(LABELA,.TRUE.,WORK(KPERTA),IRREP,ISYM,IERR,
     &             WORK(KEND2),LEND2)
      IF (IERR.NE.0 .OR. IRREP.NE.ISYMTA) THEN
        CALL QUIT('CCCR_AA: error while reading operator '//LABELA)
      END IF

* calculate first contribution:
        CALL CC_R12XI(WORK(KRHO12SQ),ISYRES,'N',WORK(KT12AMB),ISYMTB,
     &                WORK(KXINTSQ),WORK(KVXINTSQ),ISYMTA,
     &                WORK(KPERTA),WORK(KLAMDP0),WORK(KLAMDH0),
     &                'N',WORK(KEND2),LEND2)

* calculate second contribution:
        CALL CC_R12XI2A(WORK(KSCR),ISYRES,WORK(KT12AM0),1,
     &                  WORK(KPERTA),ISYMTA,WORK(KLAMDP0),
     &                  WORK(KLAMDHB),ISYMTB,'N',WORK(KEND2),LEND2)
        CALL CC_R12XI2B(WORK(KRHO12SQ),'N',WORK(KXINTSQ),1,'N',
     &                  WORK(KSCR),ISYRES,'N',-ONE)

* pack to triangular format:
        IOPT = 1
        CALL CCR12PCK2(WORK(KRHOR12),ISYRES,.FALSE.,WORK(KRHO12SQ),
     &                 'N',IOPT)
        CALL CCLR_DIASCLR12(WORK(KRHOR12),BRASCL,ISYRES)
 
      END IF
*---------------------------------------------------------------------*

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'Final result of CCCR_AA: RHO:'
        Call CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYRES,1,1)
        IF (CCR12) CALL CC_PRPR12(WORK(KRHOR12),ISYRES,1,.TRUE.)
      END IF

      RETURN
      END
*=====================================================================*
*                  END OF SUBROUTINE CCCR_AA                          *
*=====================================================================*
