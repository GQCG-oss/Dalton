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
c/* Deck CC_BAMAT */
*=====================================================================*
      SUBROUTINE CC_BAMAT(IBATRAN,NBATRAN,LISTO,LISTA,LISTB,IOPTRES,
     &                    FILBAM, IBDOTS, BCONS,MXVEC,WORK,LWORK   )
*---------------------------------------------------------------------*
*
*    Purpose: linear transformation of two CC amplitude vectors, 
*             T^A and T^B, with the CC B{O} matrix
*          
*             The linear transformations are calculated for a list
*             of operators and T^A and T^B vectors: 
*
*               LISTO        -- 'o1'
*               LISTA        -- type of T^A vectors
*               LISTB        -- type of T^B vectors
*               IBATRAN(1,*) -- indeces of the operators
*               IBATRAN(2,*) -- indeces of T^A vectors
*               IBATRAN(3,*) -- indeces of T^B vectors
*               IBATRAN(4,*) -- indeces or addresses of result vectors
*               NBATRAN      -- number of requested transformations
*               FILBAM       -- file name / list type of result vectors
*                               or list type of vectors to be dotted on
*                IBDOTS      -- indeces of vectors to be dotted on
*                BCONS       -- contains the dot products on return
*
*    return of the result vectors:
*
*          IOPTRES = 0 :  all result vectors are written to a direct
*                         access file, FILBAM is used as file name
*                         the start addresses of the vectors are
*                         returned in IBATRAN(4,*)
*
*          IOPTRES = 1 :  the vectors are kept and returned in WORK
*                         if possible, start addresses returned in
*                         IBATRAN(4,*). N.B.: if WORK is not large
*                         enough iopt is automatically reset to 0!!!
*
*          IOPTRES = 3 :  each result vector is written to its own
*                         file by a call to CC_WRRSP, FILBAM is used
*                         as list type and IBATRAN(4,*) as list index
*                         NOTE that IBATRAN(4,*) is in this case input!
*
*          IOPTRES = 4 :  each result vector is written/added to a
*                         file by a call to CC_WARSP, FILBAM is used
*                         as list type and IBATRAN(4,*) as list index
*                         NOTE that IBATRAN(4,*) is in this case input!
*
*          IOPTRES = 5 :  the result vectors are dotted on a array
*                         of vectors, the type of the arrays given
*                         by FILBAM and the indeces from IBDOTS
*                         the result of the dot products is returned
*                         in the BCONS array
*
*     Written by Christof Haettig, Maj 1997.
*
*     Adapted for CC-R12: Chrsitian Neiss, July 2005
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <ccsdinp.h>
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccroper.h>
#include <cclists.h>
#include <r12int.h>
#include <ccr12int.h>

* local parameters:
      CHARACTER MSGDBG*(18)
      PARAMETER (MSGDBG='[debug] CC_BAMAT> ')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER KDUM
      PARAMETER ( KDUM = +99 999 999 ) ! dummy address for work space

      INTEGER ISYMT0
      PARAMETER ( ISYMT0 = 1 ) ! symmetry of the reference state

      INTEGER LUBMAT
      
      CHARACTER*(*) LISTO, LISTA, LISTB, FILBAM
      INTEGER IOPTRES
      INTEGER NBATRAN, MXVEC, LWORK
      INTEGER IBATRAN(MXDIM_BATRAN,NBATRAN)
      INTEGER IBDOTS(MXVEC,NBATRAN)

#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL BCONS(MXVEC,NBATRAN)
      REAL ZERO, ONE, TWO, HALF
      REAL XNORM
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION BCONS(MXVEC,NBATRAN)
      DOUBLE PRECISION ZERO, ONE, TWO, HALF
      DOUBLE PRECISION XNORM
#endif 
      PARAMETER (ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0, HALF = 0.5d0)

      CHARACTER*(10) MODEL, MODELW
      CHARACTER*(8)  LABELO

      INTEGER KEND1, KEND2, LEND2, LENALL, KEND3, LEND3
      INTEGER ISYMTA, ISYMTB, ISYMO, ISYRES, ISYMAO, ISYMBO, IOPTW
      INTEGER IDLSTO, IDLSTA, IDLSTB, IOPTWE
      INTEGER KPERT, KOO, KOV, KVV, KAOO, KBOO, KAVV, KBVV, KSCR
      INTEGER KT1AMPA, KT1AMPB, KT1AMP0, KLAMDP0, KLAMDH0
      INTEGER KTHETA0, KTHETA1, KTHETA2, KT2AMPA, KT2AMPB
      INTEGER IRREP, IERR, IERRB, IOPT, ISYM, LEN, ITRAN, IADRTH,
     &        KTHETA1EFF, KTHETA2EFF
      INTEGER IOPTWR12,LENMOD,KTHETAR12,KATRANR12
      INTEGER KLAMDPB,KLAMDHB,KLAMDPA,KLAMDHA,KT12AMP,KXINTTRI,KXINTSQ
      INTEGER LUNIT,IAN,DUMMY,ISYM1,ISYM2,IDLST1,IDLST2
      CHARACTER APROXR12*3,LIST1*3,LIST2*3

* external functions:
      INTEGER ILSTSYM
#if defined (SYS_CRAY) 
      REAL DDOT, FREQLST, FREQB
#else 
      DOUBLE PRECISION DDOT, FREQLST, FREQB
#endif 

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        Call AROUND('ENTERED CC_BAMAT')
        WRITE (LUPRI,*) 'LISTO : ',LISTO
        WRITE (LUPRI,*) 'LISTA : ',LISTA
        WRITE (LUPRI,*) 'LISTB : ',LISTB
        WRITE (LUPRI,*) 'FILBAM: ',FILBAM
        WRITE (LUPRI,*) 'NBATRAN: ',NBATRAN
        WRITE (LUPRI,*) 'IOPTRES:',IOPTRES
        CALL FLSHFO(LUPRI)
      END IF
      
      ! well, this is no longer true, since CC3 is implemented
      ! but it is not yet debugged... though, be aware!!
      IF (CCSDT) THEN
        WRITE(LUPRI,'(/1x,a)') 'B{O} matrix transformations not '
     &          //'implemented for triples yet...'
        CALL QUIT('Triples not implemented for B '//
     &            'matrix transformations')
      END IF

      IF ( .NOT. (CCS .OR. CC2 .OR. CCSD .OR. CC3) ) THEN
        WRITE(LUPRI,'(/1x,a)') 'CC_BAMAT called for a Coupled Cluster '
     &          //'method not implemented in CC_BAMAT...'
        CALL QUIT('Unknown CC method in CC_BAMAT.')
      END IF

      IF (LISTO(1:2).NE.'o1') THEN
        WRITE (LUPRI,*) 'LISTO must refer to operator list '//
     &        'o1 in CC_BAMAT.'
        CALL QUIT('Illegal LISTO in CC_BAMAT.')
      END IF

      IF (LISTA(1:1).NE.'R' .OR. LISTB(1:1).NE.'R') THEN
        WRITE(LUPRI,*) 'LISTA and LISTB must refer to t-amplitude',
     &                    ' vectors in CC_BAMAT.'
        CALL QUIT('Illegal LISTA or LISTB in CC_BAMAT.')
      END IF

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
         CALL QUIT('Unknown coupled cluster model in CC_BAMAT.')
      END IF
      IF (CCR12) THEN
        APROXR12 = '   '
        CALL CCSD_MODEL(MODELW,LENMOD,10,MODELW,10,APROXR12)
        IOPTWR12 = 32
      END IF

* check return option for the result vectors:
      IF (IOPTRES .EQ. 0 .OR. IOPTRES .EQ. 1) THEN

        LUBMAT = -1
        CALL WOPEN2(LUBMAT,FILBAM,64,0)

      ELSE IF (IOPTRES.EQ.3 .OR. IOPTRES.EQ.4) THEN
        CONTINUE
      ELSE IF (IOPTRES.EQ.5) THEN
        IF (MXVEC*NBATRAN.NE.0) CALL DZERO(BCONS,MXVEC*NBATRAN)
      ELSE
        CALL QUIT('Illegal value of IOPTRES in CC_BAMAT.')
      END IF

*=====================================================================*
* calculate B matrix transformations:
*=====================================================================*

      KEND1  = 1

      IADRTH = 1

      DO ITRAN = 1, NBATRAN

        IDLSTO = IBATRAN(1,ITRAN)
        IDLSTA = IBATRAN(2,ITRAN)
        IDLSTB = IBATRAN(3,ITRAN)

        LABELO = LBLOPR(IDLSTO)

        ISYMO  = ILSTSYM(LISTO,IDLSTO)
        ISYMTA = ILSTSYM(LISTA,IDLSTA)
        ISYMTB = ILSTSYM(LISTB,IDLSTB)
        ISYMAO = MULD2H(ISYMTA,ISYMO)
        ISYMBO = MULD2H(ISYMTB,ISYMO)
        ISYRES = MULD2H(MULD2H(ISYMTA,ISYMTB),MULD2H(ISYMT0,ISYMO))

*---------------------------------------------------------------------*
* allocate work space for the result vector(s):
*---------------------------------------------------------------------*
        IF (CCS) THEN
          KTHETA1 = KEND1
          KTHETA2 = KDUM
          KEND2   = KTHETA1 + NT1AM(ISYRES)
        ELSE 
          KTHETA1 = KEND1
          KTHETA2 = KTHETA1 + NT1AM(ISYRES)
          KEND2   = KTHETA2 + NT2AM(ISYRES)
          IF (CCR12) THEN
            KTHETAR12 = KTHETA2 + NT2AM(ISYRES)
            KEND2     = KTHETAR12 + NTR12AM(ISYRES)
          END IF
        END IF
        IF (CCSDT .AND. IOPTRES.NE.5) THEN
          KTHETA1EFF = KEND2
          KTHETA2EFF = KTHETA1EFF + NT1AM(ISYRES)
          KEND2      = KTHETA2EFF + NT2AM(ISYRES)
        END IF

        IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'B{O} matrix transformation for ITRAN,',ITRAN
         WRITE (LUPRI,*) 'IADRTH:',IADRTH
         WRITE (LUPRI,*) 'LISTO,IDLSTO:',LISTO,IDLSTO
         WRITE (LUPRI,*) 'LISTA,IDLSTA:',LISTA,IDLSTA
         WRITE (LUPRI,*) 'LISTB,IDLSTB:',LISTB,IDLSTB
         WRITE (LUPRI,*) 'ISYMO,ISYMTA,ISYMTB:',ISYMO,ISYMTA,ISYMTB
         CALL FLSHFO(LUPRI)
        END IF

*---------------------------------------------------------------------*
* allocate memory for property integrals and response vectors: 
*  1) AO/MO perturbation integrals (complete matrix)
*  2) MO transformed perturbation integrals (occ/occ block)
*  3) MO transformed perturbation integrals (vir/vir block)
*  4) MO transformed perturbation integrals (occ/vir block)
*  5) singles excitation part of T^A
*  6) singles excitation part of T^B
*  7) singles excitation part of T^0
*  8) zeroth-order lambda particle matrix
*  9) zeroth-order lambda hole matrix
* 10) a scratch vector with the size of the singles result vector
*---------------------------------------------------------------------*
        KPERT   = KEND2
        KOO     = KPERT   + N2BST(ISYMO)
        KOV     = KOO     + NMATIJ(ISYMO)
        KVV     = KOV     + NT1AM(ISYMO)
        KAOO    = KVV     + NMATAB(ISYMO)
        KBOO    = KAOO    + NMATIJ(ISYMAO)
        KAVV    = KBOO    + NMATIJ(ISYMBO)
        KBVV    = KAVV    + NMATAB(ISYMAO)
        KT1AMPA = KBVV    + NMATAB(ISYMBO)
        KT1AMPB = KT1AMPA + NT1AM(ISYMTA)
        KT1AMP0 = KT1AMPB + NT1AM(ISYMTB)
        KLAMDP0 = KT1AMP0 + NT1AM(ISYMT0)
        KLAMDH0 = KLAMDP0 + NGLMDT(ISYMT0)
        KSCR    = KLAMDH0 + NGLMDT(ISYMT0)
        KEND2   = KSCR    + NT1AM(ISYRES)
        LEND2   = LWORK - KEND2

        IF (LEND2 .LE. 0) THEN
          CALL QUIT('Insufficient work space in CC_BAMAT. (8)')
        END IF

* read single excitation part of T^A:
        IOPT = 1
        CALL CC_RDRSP(LISTA,IDLSTA,ISYMTA,IOPT,MODEL,
     &                WORK(KT1AMPA),WORK(KDUM))

* read single excitation part of T^B:
        IOPT = 1
        CALL CC_RDRSP(LISTB,IDLSTB,ISYMTB,IOPT,MODEL,
     &                WORK(KT1AMPB),WORK(KDUM))

* read single excitation part of zeroth-order coupled cluster vector:
        IOPT = 1
        CALL CC_RDRSP('R0',0,ISYMT0,IOPT,MODEL,
     &                WORK(KT1AMP0),WORK(KDUM))

* get zeroth-order Lambda matrices:
        CALL LAMMAT(WORK(KLAMDP0),WORK(KLAMDH0),WORK(KT1AMP0),
     &              WORK(KEND2),LEND2)

* read the AO integrals:
      CALL CCPRPAO(LABELO,.TRUE.,WORK(KPERT),IRREP,ISYM,IERR,
     &             WORK(KEND2),LEND2)
      IF (IERR.NE.0 .OR. IRREP.NE.ISYMO) THEN
        CALL QUIT('CC_BAMAT: error while reading operator '//LABELO)
      END IF

* transform one-electron integrals in place:
      CALL CC_FCKMO(WORK(KPERT),WORK(KLAMDP0),WORK(KLAMDH0),
     &              WORK(KEND2),LEND2,ISYMO,1,1)

* gather occ/occ, occ/vir, and vir/vir blocks:
      CALL CC_GATHEROO(WORK(KPERT),WORK(KOO),ISYMO)
      CALL CC_GATHEROV(WORK(KPERT),WORK(KOV),ISYMO)
      CALL CC_GATHERVV(WORK(KPERT),WORK(KVV),ISYMO)

*---------------------------------------------------------------------*
* calculate O^{A} = [O, T1^A]
*---------------------------------------------------------------------*
* Ftilde^{A}, occupied/occupied blocks:
        CALL CCG_1ITROO(WORK(KAOO),ISYMAO,
     &                  WORK(KOV), ISYMO, WORK(KT1AMPA),ISYMTA )

* Ftilde^{B}, occupied/occupied blocks:
        CALL CCG_1ITROO(WORK(KBOO),ISYMBO,
     &                  WORK(KOV), ISYMO, WORK(KT1AMPB),ISYMTB )

* Ftilde^{A}, virtual/virtual blocks:
        CALL CCG_1ITRVV(WORK(KAVV),ISYMAO, 
     &                  WORK(KOV), ISYMO, WORK(KT1AMPA),ISYMTA  )

* Ftilde^{B}, virtual/virtual blocks:
        CALL CCG_1ITRVV(WORK(KBVV),ISYMBO, 
     &                  WORK(KOV), ISYMO, WORK(KT1AMPB),ISYMTB  )


        IF (LOCDBG) THEN
          XNORM=DDOT(NMATIJ(ISYMTA),WORK(KAOO),1,WORK(KAOO),1)
          WRITE (LUPRI,*) 'Norm of O^AOO:',XNORM
          XNORM=DDOT(NMATAB(ISYMTA),WORK(KAVV),1,WORK(KAVV),1)
          WRITE (LUPRI,*) 'Norm of O^AVV:',XNORM
          XNORM=DDOT(NMATIJ(ISYMTB),WORK(KBOO),1,WORK(KBOO),1)
          WRITE (LUPRI,*) 'Norm of O^BOO:',XNORM
          XNORM=DDOT(NMATAB(ISYMTB),WORK(KBVV),1,WORK(KBVV),1)
          WRITE (LUPRI,*) 'Norm of O^BVV:',XNORM
          WRITE (LUPRI,*) 'T^A (singles only):'
          Call CC_PRP(WORK(KT1AMPA),WORK,ISYMTA,1,0)
          WRITE (LUPRI,*) 'T^B (singles only):'
          Call CC_PRP(WORK(KT1AMPB),WORK,ISYMTB,1,0)
          CALL FLSHFO(LUPRI)
        END IF

*---------------------------------------------------------------------*
* initialize the singles part of the result vector THETA:
*---------------------------------------------------------------------*
        CALL DZERO(WORK(KTHETA1),NT1AM(ISYRES))

*---------------------------------------------------------------------*
* J contribution: vir/occ block of double transformed integrals:
*---------------------------------------------------------------------*

* 1. contribution [O^A,T^B]:
        IOPT = 1
        CALL CCG_1ITRVO(WORK(KSCR),ISYRES, WORK(KBOO),
     &                  WORK(KBVV),ISYMBO,
     &                  WORK(KT1AMPA),ISYMTA,IOPT )

        CALL DAXPY(NT1AM(ISYRES),HALF,WORK(KSCR),1,WORK(KTHETA1),1)

* 2. contribution [O^B,T^A]:
        IOPT = 1
        CALL CCG_1ITRVO(WORK(KSCR),ISYRES, WORK(KAOO),
     &                  WORK(KAVV),ISYMAO,
     &                  WORK(KT1AMPB),ISYMTB,IOPT )

        CALL DAXPY(NT1AM(ISYRES),HALF,WORK(KSCR),1,WORK(KTHETA1),1)


        IF (LOCDBG) THEN
          XNORM=DDOT(NT1AM(ISYRES),WORK(KTHETA1),1,WORK(KTHETA1),1)
          WRITE (LUPRI,*) 'Norm of O^ABOV (J contribution):',XNORM
C         WRITE (LUPRI,'(/5X,A)') 'AVV / AOO:'
C         CALL CC_PREI(WORK(KAVV),WORK(KAOO),ISYMAO,1)
C         WRITE (LUPRI,'(/5X,A)') 'BVV / BOO:'
C         CALL CC_PREI(WORK(KBVV),WORK(KBOO),ISYMBO,1)
          CALL FLSHFO(LUPRI)
        END IF

     
*---------------------------------------------------------------------*
* initialize the doubles part of the result vector THETA:
*---------------------------------------------------------------------*
      IF (.NOT. CCS ) THEN
        CALL DZERO( WORK(KTHETA2), NT2AM(ISYRES) )
      END IF

*----------------------------------------
* first E1 & E2 contributions, T^B x F^A:
*----------------------------------------
C     IF (.NOT. (CCS .OR. CC2 .OR. CCSTST) ) THEN
      IF (.NOT. (CCS .OR. CCSTST) ) THEN
        KT2AMPB = KEND2
        KEND3   = KT2AMPB + NT2SQ(ISYMTB)
        LEND3   = LWORK - KEND3 

        IF (LEND3 .LT. NT2AM(ISYMTB)) THEN
          CALL QUIT('Insufficient work space in CC_BAMAT. (15)')
        END IF

        IOPT = 2
        CALL CC_RDRSP(LISTB,IDLSTB,ISYMTB,IOPT,
     &                MODEL,WORK(KDUM),WORK(KEND3))

        CAll CCLR_DIASCL(WORK(KEND3),TWO,ISYMTB)

        CALL CC_T2SQ(WORK(KEND3),WORK(KT2AMPB),ISYMTB)

* calculate the contribution to THETA2:
        CALL CCRHS_E(WORK(KTHETA2),WORK(KT2AMPB),WORK(KAVV),
     &               WORK(KAOO),WORK(KEND3),LEND3,ISYMTB,ISYMAO)

        IF (LOCDBG) THEN
          XNORM=DDOT(NMATIJ(ISYMAO),WORK(KAOO),1,WORK(KAOO),1)
          WRITE (LUPRI,*) 'Norm of KAOO intermediate:',XNORM
          XNORM=DDOT(NMATAB(ISYMAO),WORK(KAVV),1,WORK(KAVV),1)
          WRITE (LUPRI,*) 'Norm of KAVV intermediate:',XNORM
          XNORM=DDOT(NT2AM(ISYRES),WORK(KTHETA2),1,WORK(KTHETA2),1)
          WRITE (LUPRI,*) 'Norm of THETA2 after first E contribution:',
     &         XNORM
          Call AROUND('final result for B{O} matrix transformation:')
          Call CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYRES,1,1)
          CALL FLSHFO(LUPRI)
        END IF

      END IF ! (.NOT. (CCS .OR. CC2 .OR. CCSTST))

*-----------------------------------------
* second E1 & E2 contributions, T^A x F^B:
*-----------------------------------------
      IF (.NOT. (CCS .OR. CCSTST) ) THEN
        KT2AMPA = KEND2
        KEND3   = KT2AMPA + NT2SQ(ISYMTA)
        LEND3   = LWORK - KEND3 

        IF (LEND3 .LT. NT2AM(ISYMTA)) THEN
          CALL QUIT('Insufficient work space in CC_BAMAT. (16)')
        END IF

        IOPT = 2
        CALL CC_RDRSP(LISTA,IDLSTA,ISYMTA,IOPT,
     &                MODEL,WORK(KDUM),WORK(KEND3))

        CAll CCLR_DIASCL(WORK(KEND3),TWO,ISYMTA)

        CALL CC_T2SQ(WORK(KEND3),WORK(KT2AMPA),ISYMTA)

* calculate the contribution to THETA2:
        CALL CCRHS_E(WORK(KTHETA2),WORK(KT2AMPA),WORK(KBVV),
     &               WORK(KBOO),WORK(KEND3),LEND3,ISYMTA,ISYMBO)

        IF (LOCDBG) THEN
          XNORM=DDOT(NMATIJ(ISYMBO),WORK(KBOO),1,WORK(KBOO),1)
          WRITE (LUPRI,*) 'Norm of KBOO intermediate:',XNORM
          XNORM=DDOT(NMATAB(ISYMBO),WORK(KBVV),1,WORK(KBVV),1)
          WRITE (LUPRI,*) 'Norm of KBVV intermediate:',XNORM
          XNORM=DDOT(NT2AM(ISYRES),WORK(KTHETA2),1,WORK(KTHETA2),1)
          WRITE (LUPRI,*) 'Norm of THETA2 after second E contribution:',
     &         XNORM
          Call AROUND('final result for B{O} matrix transformation:')
          Call CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYRES,1,1)
          CALL FLSHFO(LUPRI)
        END IF

      END IF ! (.NOT. (CCS .OR. CCSTST))

*---------------------------------------------------------------------*
* initialize the R12 doubles part of the result vector THETA:
*---------------------------------------------------------------------*
      IF ( CCR12 ) THEN
        CALL DZERO(WORK(KTHETAR12),NTR12AM(ISYRES))

        KLAMDPB = KEND2
        KLAMDHB = KLAMDPB + NGLMDT(ISYMTB)
        KLAMDPA = KLAMDHB + NGLMDT(ISYMTB)
        KLAMDHA = KLAMDPA + NGLMDT(ISYMTA)
        KT12AMP = KLAMDHA + NGLMDT(ISYMTA)
        KXINTTRI= KT12AMP + MAX(NTR12SQ(ISYMTA),NTR12SQ(ISYMTB))
        KXINTSQ = KXINTTRI + MAX(NR12R12P(1),NTR12SQ(ISYRES))
        KSCR    = KXINTSQ + NR12R12SQ(1)
        KEND3   = KSCR + NTR12SQ(ISYRES)
        LEND3   = LWORK - KEND3
        IF (LEND3 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_BAMAT')
        END IF

        CALL CCPRPAO(LABELO,.TRUE.,WORK(KPERT),IRREP,ISYM,IERR,
     &             WORK(KEND2),LEND2)
        IF (IERR.NE.0 .OR. IRREP.NE.ISYMO) THEN
          CALL QUIT('CC_BAMAT: error while reading operator '//LABELO)
        END IF

        CALL CCLR_LAMTRA(WORK(KLAMDP0),WORK(KLAMDPA), WORK(KLAMDH0),
     &                   WORK(KLAMDHA),WORK(KT1AMPA),ISYMTA)

        CALL CCLR_LAMTRA(WORK(KLAMDP0),WORK(KLAMDPB), WORK(KLAMDH0),
     &                   WORK(KLAMDHB),WORK(KT1AMPB),ISYMTB)

        LUNIT = -1
        CALL GPOPEN(LUNIT,FCCR12X,'OLD',' ','UNFORMATTED',
     &              DUMMY,.FALSE.)
        REWIND(LUNIT)
 8888   READ(LUNIT) IAN
        READ(LUNIT) (WORK(KXINTTRI-1+I), I=1, NR12R12P(1))
        IF (IAN.NE.IANR12) GOTO 8888
        CALL GPCLOSE(LUNIT,'KEEP')
        IOPT = 2
        CALL CCR12UNPCK2(WORK(KXINTTRI),1,WORK(KXINTSQ),'N',IOPT)
C
        DO I = 1, 2
          IF (I.EQ.1) THEN
            ISYM1 = ISYMTA
            ISYM2 = ISYMTB
            LIST1 = LISTA
            IDLST1 = IDLSTA
          ELSE IF (I.EQ.2) THEN
            ISYM1 = ISYMTB
            ISYM2 = ISYMTA
            LIST1 = LISTB
            IDLST1 = IDLSTB
          END IF
          !read R12 response amplitudes
          CALL CC_R12GETCT(WORK(KT12AMP),ISYM1,2,KETSCL,.FALSE.,'N',
     &                 DUMMY,DUMMY,DUMMY,LIST1,IDLST1,WORK(KEND3),LEND3)
          !calculate....
          IF (I.EQ.1) THEN
            CALL CC_R12XI2A(WORK(KSCR),ISYRES,WORK(KT12AMP),ISYM1,
     &                      WORK(KPERT),ISYMO,WORK(KLAMDP0),
     &                      WORK(KLAMDHB),ISYM2,'N',
     &                      WORK(KEND3),LEND3)
            CALL DCOPY(NTR12SQ(ISYRES),WORK(KSCR),1,WORK(KXINTTRI),1)
          ELSE IF (I.EQ.2) THEN
            CALL CC_R12XI2A(WORK(KSCR),ISYRES,WORK(KT12AMP),ISYM1,
     &                      WORK(KPERT),ISYMO,WORK(KLAMDP0),
     &                      WORK(KLAMDHA),ISYM2,'N',
     &                      WORK(KEND3),LEND3)
          END IF
        END DO       
C
        CALL DAXPY(NTR12SQ(ISYRES),ONE,WORK(KXINTTRI),1,WORK(KSCR),1)
        CALL DZERO(WORK(KXINTTRI),NTR12SQ(ISYRES))
        CALL CC_R12XI2B(WORK(KXINTTRI),'N',WORK(KXINTSQ),1,'N',
     &                  WORK(KSCR),ISYRES,'N',-ONE)
C
        IOPT = 1
        CALL CCR12PCK2(WORK(KTHETAR12),ISYRES,.FALSE.,WORK(KXINTTRI),
     &                 'N',IOPT)
        CALL CCLR_DIASCLR12(WORK(KTHETAR12),BRASCL,ISYRES)
C
      END IF

*---------------------------------------------------------------------*
* compute CC3 contributions:
*---------------------------------------------------------------------*
      IF (CCSDT) THEN

         IF (IOPTRES.EQ.5) THEN
           FREQB = 0.0D0
         ELSE
           FREQB = FREQLST(FILBAM,IBATRAN(4,ITRAN))
         END IF

         CALL CCSDT_BAMAT_NODDY(IOPTRES,FREQB,LABELO,ISYMO,
     &                          LISTA,IDLSTA,
     &                          LISTB,IDLSTB,
     &                          WORK(KTHETA1),WORK(KTHETA2),
     &                          WORK(KTHETA1EFF),WORK(KTHETA2EFF),
     &                          IBDOTS,BCONS,FILBAM,ITRAN,
     &                          NBATRAN,MXVEC,WORK(KEND2),LEND2)

      END IF

*---------------------------------------------------------------------*
* write result vector to output:
*---------------------------------------------------------------------*
      IF (IOPTRES .EQ. 0  .OR. IOPTRES .EQ. 1) THEN

*       write to a common direct access file, 
*       store start address in IBATRAN(4,ITRAN)

        IBATRAN(4,ITRAN) = IADRTH

        CALL PUTWA2(LUBMAT,FILBAM,WORK(KTHETA1),IADRTH,NT1AM(ISYRES))
        IADRTH = IADRTH + NT1AM(ISYRES)

        IF (.NOT.CCS) THEN
          CALL PUTWA2(LUBMAT,FILBAM,WORK(KTHETA2),IADRTH,NT2AM(ISYRES))
          IADRTH = IADRTH + NT2AM(ISYRES)
        END IF
        
        IF (CCR12) THEN
          CALL PUTWA2(LUBMAT,FILBAM,WORK(KTHETAR12),IADRTH,
     &                NTR12AM(ISYRES))
          IADRTH = IADRTH + NTR12AM(ISYRES)
        END IF

        IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'B{O} matrix transform. nb. ',ITRAN,
     &          ' saved on file.'
         WRITE (LUPRI,*) 'ADRESS, LENGTH:',
     &        IBATRAN(4,ITRAN),IADRTH-IBATRAN(4,ITRAN)
         XNORM = DDOT(NT1AM(ISYRES),WORK(KTHETA1),1,WORK(KTHETA1),1)
         IF (.NOT.CCS) XNORM = XNORM +
     &           DDOT(NT2AM(ISYRES),WORK(KTHETA2),1,WORK(KTHETA2),1)
         IF (CCR12) XNORM = XNORM +
     &        DDOT(NTR12AM(ISYRES),WORK(KTHETAR12),1,WORK(KTHETAR12),1)
         WRITE (LUPRI,*) 'Norm:', XNORM

         Call AROUND('B{O} matrix transformation written to file:')
         Call CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYRES,1,1)
         IF (CCR12) CALL CC_PRPR12(WORK(KTHETAR12),ISYRES,1,.TRUE.)
        END IF

      ELSE IF ( IOPTRES .EQ. 3 .OR. IOPTRES .EQ. 4 ) THEN

*        write to a sequential file by a call to CC_WRRSP/CC_WARSP,
*        use FILBAM as LIST type and IBATRAN(4,ITRAN) as index

         KTHETA0 = -999999

         IF (IOPTRES.EQ.3) THEN
           CALL CC_WRRSP(FILBAM,IBATRAN(4,ITRAN),ISYRES,IOPTW,MODELW,
     &                   WORK(KTHETA0),WORK(KTHETA1),WORK(KTHETA2),
     &                   WORK(KEND2),LEND2)
           IF (CCR12) THEN
             CALL CC_WRRSP(FILBAM,IBATRAN(4,ITRAN),ISYRES,IOPTWR12,
     &                     MODELW,DUMMY,DUMMY,WORK(KTHETAR12),
     &                     WORK(KEND2),LEND2)
           END IF
           IF (CCSDT) THEN
            CALL CC_WRRSP(FILBAM,IBATRAN(4,ITRAN),ISYRES,IOPTWE,MODEL,
     &                    WORK(KTHETA0),WORK(KTHETA1EFF),
     &                    WORK(KTHETA2EFF),WORK(KEND2),LEND2)
           END IF
         ELSE IF (IOPTRES.EQ.4) THEN
           CALL CC_WARSP(FILBAM,IBATRAN(4,ITRAN),ISYRES,IOPTW,MODELW,
     &                   WORK(KTHETA0),WORK(KTHETA1),WORK(KTHETA2),
     &                   WORK(KEND2),LEND2)
           IF (CCR12) THEN
             CALL CC_WARSP(FILBAM,IBATRAN(4,ITRAN),ISYRES,IOPTWR12,
     &                     MODELW,DUMMY,DUMMY,WORK(KTHETAR12),
     &                     WORK(KEND2),LEND2)
           END IF
           IF (CCSDT) THEN
            CALL CC_WARSP(FILBAM,IBATRAN(4,ITRAN),ISYRES,IOPTWE,MODELW,
     &                    WORK(KTHETA0),WORK(KTHETA1EFF),
     &                    WORK(KTHETA2EFF),WORK(KEND2),LEND2)
           END IF
         END IF

         IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'Write B{O} * ',LISTA,' * ',LISTB,
     &           ' transformation',
     &              ' as ',FILBAM,' type vector to file.'
           WRITE (LUPRI,*) 'index of inp. O integr:',IBATRAN(1,ITRAN)
           WRITE (LUPRI,*) 'index of inp. A vector:',IBATRAN(2,ITRAN)
           WRITE (LUPRI,*) 'index of inp. B vector:',IBATRAN(3,ITRAN)
           WRITE (LUPRI,*) 'index of result vector:',IBATRAN(4,ITRAN)
           XNORM = DDOT(NT1AM(ISYRES),WORK(KTHETA1),1,WORK(KTHETA1),1)
           IF (.NOT.CCS) XNORM = XNORM +
     &         DDOT(NT2AM(ISYRES),WORK(KTHETA2),1,WORK(KTHETA2),1)
           IF (CCR12) XNORM = XNORM +
     &         DDOT(NTR12AM(ISYRES),WORK(KTHETAR12),1,WORK(KTHETAR12),1)
           WRITE (LUPRI,*) 'norm^2 of result vector:',XNORM
         END IF

      ELSE IF ( IOPTRES .EQ. 5 ) THEN

         CALL CCDOTRSP(IBDOTS,BCONS,IOPTW,FILBAM,ITRAN,NBATRAN,MXVEC,
     &                 WORK(KTHETA1),WORK(KTHETA2),ISYRES,
     &                 WORK(KEND2),LEND2)
         IF (CCR12) THEN
           CALL CCDOTRSP(IBDOTS,BCONS,IOPTWR12,FILBAM,ITRAN,NBATRAN,
     &                   MXVEC,DUMMY,WORK(KTHETAR12),ISYRES,
     &                   WORK(KEND2),LEND2)
         END IF

      ELSE
        CALL QUIT('Illegal value for IOPTRES in CC_BAMAT.')
      END IF
*---------------------------------------------------------------------*
* End of loop over B matrix transformations
*---------------------------------------------------------------------*
      END DO

*---------------------------------------------------------------------*
* if IOPTRES=1 and enough work space available, read result
* vectors back into memory:
*---------------------------------------------------------------------*

* check size of work space:
      IF (IOPTRES .EQ. 1) THEN
        LENALL = IADRTH-1
        IF (LENALL .GT. LWORK) IOPTRES = 0
      END IF

* read the result vectors back into memory:
      IF (IOPTRES .EQ. 1) THEN

        CALL GETWA2(LUBMAT,FILBAM,WORK(1),1,LENALL)

        IF (LOCDBG) THEN
         DO ITRAN = 1, NBATRAN
           IF (ITRAN.LT.NBATRAN) THEN
             LEN     = IBATRAN(4,ITRAN+1)-IBATRAN(4,ITRAN)
           ELSE
             LEN     = IADRTH-IBATRAN(4,NBATRAN)
           END IF
           KTHETA1 = IBATRAN(4,ITRAN)
           XNORM   = DDOT(LEN, WORK(KTHETA1),1, WORK(KTHETA1),1)
           WRITE (LUPRI,*) 'Read B matrix transformation nb. ',NBATRAN
           WRITE (LUPRI,*) 'Adress, length, NORM:',IBATRAN(4,NBATRAN),
     &          LEN,XNORM
          END DO
          CALL FLSHFO(LUPRI)
        END IF
      END IF 

*---------------------------------------------------------------------*
* close B matrix file & return
*---------------------------------------------------------------------*
* check return option for the result vectors:
      IF (IOPTRES .EQ. 0 .OR. IOPTRES .EQ. 1) THEN
        CALL WCLOSE2(LUBMAT,FILBAM,'DELETE')

      ELSE IF (IOPTRES.EQ.3 .OR. IOPTRES.EQ.4 .OR. IOPTRES.EQ.5) THEN
        CONTINUE
      ELSE
        CALL QUIT('Illegal value of IOPTRES in CC_BAMAT.')
      END IF


*=====================================================================*

      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE CC_BAMAT
*=====================================================================*
