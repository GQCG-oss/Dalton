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
c /* deck ccfbtaf */
*=====================================================================*
      SUBROUTINE CCFBTAF(THETA1,THETA2,
     &                   RAIM,RQAIM,FAIM,FQAIM,
     &                   XAIM,YAIM,ZETA1,
     &                   T1AMP0,T2AMP0,
     &                   XIJK0,
     &                   ONEH0,ONEHB,
     &                   FOCK0,FOCKB, FOCKA, FOCKBA,
     &                   LUBFA, FILBFA,IBFA,                         
     &                   LUBFQA,FILBFQA,IBFQA,
     &                   LUBFZA, FILBFZA, IBFZA,
     &                   LUBFZQA, FILBFZQA, IBFZQA,
     &                   LU0IAJB,FN0IAJB,IT2DEL0,
     &                   LU1IAJB,FN1IAJB,IT2DELB,
     &                   LU0IJCB, FN0IJCB,
     &                   LU0CJIB, FN0CJIB, IT2DELA0,
     &                   LU1IJCB, FN1IJCB,
     &                   LU1CJIB, FN1CJIB, IT2DELAB,
     &                   LU1IJKL, FN1IJKL,
     &                   LUPQMO, FILPQMO,IADRPQMO,
     &                   LUPQAMO, FILPQAMO,IADPQAMO,
     &                   IADRMO,
     &                   XLAMDP0,XLAMDH0,XLAMDQP,XLAMDQH,
     &                   XLAMDAP,XLAMDAH, XLAMQAP,XLAMQAH,
     &                   LISTL,IDLSTL,LABELH,ISYHOP,
     &                   LISTR,IDLSTR,
     &                   LTWOEL,LRELAX,ITRAN, WORK, LWORK)
*---------------------------------------------------------------------*
*
*   Purpose: calculate final F^B T^A vector (THETA) from intermediates
*   OBS: LRELAX and LTWOEL are SEPARATE logicals in here
*
*   Written by Sonia Coriani, February 1999
*   Version: 02.03.2000
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccisao.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      LOGICAL LFMATST
      PARAMETER (LFMATST = .FALSE.)   !scale 1/5 if test against FMAT

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)
      INTEGER KDUM
      PARAMETER (KDUM = +99 999 999) ! dummy address in work space
      INTEGER IDUM
      PARAMETER (IDUM = 0)         ! dummy symmetry index

* variables:
      CHARACTER*(*) LISTL, LISTR, LABELH 
      CHARACTER*(*) FN1IAJB,FN0IAJB
      CHARACTER*(*) FN0IJCB,FN0CJIB,FN1IJCB,FN1CJIB
      CHARACTER*(*) FILPQMO, FILPQAMO
      CHARACTER*(*) FILBFA,FILBFQA,FILBFZA,FILBFZQA
      CHARACTER*(*) FN1IJKL
      LOGICAL LRELAX, LTWOEL
      INTEGER LWORK, IDLSTL, IDLSTR, LU1IAJB, LU0IAJB
      INTEGER LU0IJCB, LU0CJIB, LU1IJCB, LU1CJIB
      INTEGER LU1IJKL
      INTEGER LUPQMO, LUPQAMO
      INTEGER ISYHOP, ISYMTA, ISYZTA, ISYHZE, ISYHTA
      INTEGER IT2DEL0(*), IT2DELB(*), IADRPQMO(*)
      INTEGER IT2DELA0(*), IT2DELAB(*), IADPQAMO(*)
c
      INTEGER LUBFA,IBFA,LUBFQA,IBFQA,IADRMO
      INTEGER LUBFZA,IBFZA(*),LUBFZQA,IBFZQA(*)    !not sure, maybe IBFZA(MAXORB)

#if defined (SYS_CRAY)
      REAL XLAMDP0(*), XLAMDH0(*), XLAMDQP(*), XLAMDQH(*)
      REAL XLAMDAP(*), XLAMDAH(*), XLAMQAP(*), XLAMQAH(*)
      REAL THETA1(*), THETA2(*), XAIM(*), YAIM(*)
      REAL RAIM(*), RQAIM(*), FAIM(*), FQAIM(*)
      REAL ONEH0(*),ONEHB(*)
      REAL FOCK0(*),FOCKB(*), FOCKA(*),FOCKBA(*)
      REAL WORK(LWORK), ZETA1(*)
      REAL T1AMP0(*), T2AMP0(*), XIJK0(*)
      REAL HALF, ONE, ZERO, FIFTH, DUM, XNORM, FAC
      REAL DDOT, DNRM2
#else
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*), XLAMDQP(*), XLAMDQH(*)
      DOUBLE PRECISION XLAMDAP(*), XLAMDAH(*), XLAMQAP(*), XLAMQAH(*)
      DOUBLE PRECISION THETA1(*), THETA2(*), XAIM(*), YAIM(*) 
      DOUBLE PRECISION RAIM(*), RQAIM(*), FAIM(*), FQAIM(*)
      DOUBLE PRECISION ONEH0(*),ONEHB(*)
      DOUBLE PRECISION FOCK0(*),FOCKB(*), FOCKA(*),FOCKBA(*)
      DOUBLE PRECISION WORK(LWORK), ZETA1(*)
      DOUBLE PRECISION T1AMP0(*), T2AMP0(*), XIJK0(*)
      DOUBLE PRECISION HALF, ONE, ZERO, TWO, DUM,FIFTH, XNORM, FAC
      DOUBLE PRECISION DDOT, DNRM2
#endif
      PARAMETER (HALF = 0.5d0, ONE = 1.0d0, TWO = 2.0d0, ZERO = 0.0d0)
      PARAMETER (FIFTH = 0.2d0)

      CHARACTER CDUM, MODEL*(10)
      LOGICAL LGAMMA, LO3BF, LBFZETA, LRCON, LGCON, FCKCON, LOCC
      LOGICAL LE1, LE2, LDUMMY, LDERRLX, LTRSPF
      INTEGER ISYCTR, ISYRES, KZETA2, KZETPK, KBF
      INTEGER ISYMI, ISYMA, KOFF1, KOFF2, LENBF0, LENBF1, KEND1, LWRK1
      INTEGER IOPT, ICON, KA2IM, KA2RS, KXPCK, KXIAJB, KT2AMP, IOPTG
      INTEGER KB2CON, KEMAT1, KEMAT2, KCTR2, KEND2, LWRK2, MT2BCD
      INTEGER ISYIFJL, KXIFJL, KYJLIF, KPINT0, KQINT0, KWINT0
      INTEGER ISYMF, ISYALJ, ISYJLI, IVIRF, KYPS0, KYPS1, LEN
      INTEGER NTAIKJ(9), ITAIKJ(8,8), ICOUNT, ISYM, ISYMAIK, ISYMJ
      INTEGER ISYMBF, KCDBAR, IADR, IORDER
      INTEGER LENBFA, LENBFQA, KEMAT1A
      INTEGER KGAIM, KGQAI, KGAMQA, KBFA, KEND0, LWRK0, KEND3, LWRK3
      INTEGER KPAINT, KQAINT, KWAINT,KEND4,LWRK4,KTA1AM,KTA2AM,KHFPRJ
      INTEGER KGZAIM, KGZQAI, KBFZA, KBFZQA, KCTILD, KDTILD
      INTEGER KXIJCB, KXCJIB, KGAMMAQ, LENGAM, KYPSA,KYPSQA,KXBICJ
      INTEGER ISYMCJ, ISYMBI, NDIM1, NDIM2, ISYIAJB, ITRAN, KTEST
      INTEGER IOPTTCME, IOPTR12
* external 
      INTEGER ILSTSYM

      CALL QENTER('CCFBTAF')

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
* set symmetries of various intermediates and/or contributions
*
      ISYCTR = ILSTSYM(LISTL,IDLSTL)
      ISYMTA = ILSTSYM(LISTR,IDLSTR)
      ISYZTA = MULD2H(ISYMTA,ISYCTR)
      ISYHZE = MULD2H(ISYHOP,ISYCTR)
      ISYHTA = MULD2H(ISYMTA,ISYHOP)
      ISYRES = MULD2H(ISYHOP,ISYZTA)
*
* set NTAIKJ and ITAIKJ arrays:
*
      DO ISYM = 1, NSYM
        ICOUNT = 0
        DO ISYMAIK = 1, NSYM
           ISYMJ  = MULD2H(ISYMAIK,ISYM)
           ITAIKJ(ISYMAIK,ISYMJ) = ICOUNT
           ICOUNT = ICOUNT + NT2BCD(ISYMAIK)*NRHF(ISYMJ)
        END DO
        NTAIKJ(ISYM) = ICOUNT
      END DO
*    
* initialize single excitation part of result vector (THETA1):
* for the given ITRAN
*
      CALL DZERO(THETA1,NT1AM(ISYRES))  
*
* initialize double excitation part of the result vector (THETA2):
*
      IF (.NOT.CCS)
     &    CALL DZERO(THETA2,NT2AM(ISYRES))         

*---------------------------------------------------------------------*
*     transform one-electron hamiltonian and ``barred'' Fock matrix
*     (FOCKB) to MO basis and add extra relaxation contributions, as well
*     as FOCKA and FOCKBA contributions:
*     Save on FOCKB. Keep FOCKBA for later use
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
c         WRITE (LUPRI,*) 'CCFBTAF> ONEHB matrix in AO:'
c         CALL CC_PRFCKAO(ONEHB,ISYHOP)
c         WRITE (LUPRI,*) 'CCFBTAF> ONEH0 matrix in AO:'
c         CALL CC_PRFCKAO(ONEH0,ISYM0)
         XNORM = DNRM2(N2BST(ISYM0),FOCK0,1)
         WRITE (LUPRI,*) 'CCFBTAF> FOCK0 matrix in AO:', XNORM
c         CALL CC_PRFCKAO(FOCK0,ISYM0)

         XNORM = DNRM2(N2BST(ISYHOP),FOCKB,1)
         WRITE (LUPRI,*) 'CCFBTAF> FOCKB matrix in AO:', XNORM
c         CALL CC_PRFCKAO(FOCKB,ISYHOP)

         XNORM = DNRM2(N2BST(ISYMTA),FOCKA,1)
         WRITE (LUPRI,*) 'CCFBTAF> FOCKA matrix in AO:', XNORM
c         CALL CC_PRFCKAO(FOCKA,ISYMTA)

         XNORM = DNRM2(N2BST(ISYHTA),FOCKBA,1)
         WRITE (LUPRI,*) 'CCFBTAF> FOCKBA matrix in AO:', XNORM
c         CALL CC_PRFCKAO(FOCKBA,ISYHTA)
      END IF
*
      IF (.NOT.CC2) THEN

         IORDER = 1
         CALL CC_FCKRLX(ONEH0,ONEHB,DUM,DUM,
     &                  XLAMDP0,XLAMDH0,ISYM0,
     &                  XLAMDQP,XLAMDQH,ISYHOP,LRELAX,
     &                  DUM,DUM,IDUM,LDUMMY,
     &                  DUM,DUM,IDUM,LDUMMY,
     &                  IORDER,WORK,LWORK)

         CALL CC_FCKRLX1(FOCKBA,FOCKA,ISYHTA,ISYMTA,
     &                   XLAMDP0,XLAMDH0,ISYM0,ISYM0,
     &                   XLAMDQP,XLAMDQH,ISYHOP,ISYHOP,
     &                   LRELAX,WORK,LWORK)

*
*  add extra terms 
*
*         IF (CCS) THEN
*
* two steps: first get F^(1)_pq, then transform with T^A
*
*            IORDER = 1
*            CALL CC_FCKRLX(FOCK0,FOCKB,DUM,DUM,
*     &                  XLAMDP0,XLAMDH0,ISYM0,
*     &                  XLAMDQP,XLAMDQH,ISYHOP,LRELAX,
*     &                  DUM,DUM,IDUM,LDUMMY,
*     &                  DUM,DUM,IDUM,LDUMMY,
*     &                  IORDER,WORK,LWORK)
*
*            LE1   = .TRUE.                       
*            LE2   = .FALSE.                       
*            CALL CC_HPQTA(FOCKB,ISYHOP,FOCKB,DUM,ISYHTA,
*     &                      LISTR,IDLSTR,ISYMTA,LE1,LE2,
*     &                      WORK,LWORK)
*         ELSE
*
           CALL CC_FCKRLX1(FOCKB,FOCK0,ISYHOP,ISYM0,
     &                  XLAMDAP,XLAMDH0,ISYMTA,ISYM0,
     &                  XLAMQAP,XLAMDQH,ISYHTA,ISYHOP,
     &                  LRELAX,WORK,LWORK)
*         END IF

         CALL DAXPY(N2BST(ISYHTA),ONE,FOCKBA,1,FOCKB,1)

      ELSE
          CALL QUIT('CC2 not yet available')
      END IF
*
* Note that I keep FOCKBA for later use (L0). The global Fock
* transformed matrix is in FOCKB
*
      IF (LOCDBG) THEN
c         WRITE (LUPRI,*) 'LRELAX: ', LRELAX
c         WRITE (LUPRI,*) 'Lambda0'
c         CALL CC_PRLAM(XLAMDP0,XLAMDH0,ISYM0)
c         WRITE (LUPRI,*) 'LambdaQ'
c         CALL CC_PRLAM(XLAMDQP,XLAMDQH,ISYHOP)
c         WRITE (LUPRI,*) 'LambdaA'
c         CALL CC_PRLAM(XLAMDAP,XLAMDAH,ISYMTA)
c         WRITE (LUPRI,*) 'LambdaQA'
c         CALL CC_PRLAM(XLAMQAP,XLAMQAH,ISYHTA)
      END IF  
      IF (LOCDBG) THEN
c         WRITE (LUPRI,*) 'CCFBTAF> ONEHB matrix in MO:'
c         CALL CC_PRFCKMO(ONEHB,ISYHOP)
c         WRITE (LUPRI,*) 'CCFBTAF> FOCKB(+TA contribution) matrix in MO:'
c         CALL CC_PRFCKMO(FOCKB,ISYHTA)
c
c         WRITE (LUPRI,*) 'CCFBTAF> FOCKBA matrix in MO:'
c         CALL CC_PRFCKMO(FOCKBA,ISYHTA)
      END IF
*
* ---------------------------------------------------------------------*
* (move here <HF| contribution if calculated from FOCKBA)
* Use FOCKBA instead (store ai)
* ---------------------------------------------------------------------*
*
      IF ( (.NOT.CCS) .AND. (LTWOEL.OR.LRELAX) ) THEN
        IF ( LISTL(1:2).EQ.'L0') THEN
          DO ISYMI = 1, NSYM
             ISYMA = MULD2H(ISYMI,ISYHTA)    !ISYHTA = ISYRES here
             DO I = 1, NRHF(ISYMI)
                DO A = 1, NVIR(ISYMA)
                   KOFF1 = IT1AM(ISYMA,ISYMI)  + NVIR(ISYMA)*(I-1)+A
                   KOFF2 = IFCVIR(ISYMI,ISYMA) + NORB(ISYMI)*(A-1)+I
                   THETA1(KOFF1) = TWO * FOCKBA(KOFF2)
                END DO
             END DO
          END DO
          IF (LOCDBG) THEN
             CALL AROUND('FBTAF: THETA1 after <HF| term: from FOCKBA')
             CALL CC_PRP(THETA1,THETA2,ISYRES,1,0)
             XNORM = DNRM2(NT1AM(ISYRES),THETA1,1)
             WRITE (LUPRI,*) 'Norm of Theta(<HF|)_ai:', XNORM
          END IF
        END IF                
      END IF
*
*---------------------------------------------------------------------*
*     calculate D1 contributions:
*     both relaxed/2-el and unrelaxed/one-el cases
*---------------------------------------------------------------------*
*
      IF (.NOT. CCS) THEN

        CALL CC_21EFM(THETA1,ONEHB,ISYHOP,XAIM,YAIM,ISYZTA,WORK,LWORK)

        IF (LOCDBG) THEN
          CALL AROUND('CCFBTAF: THETA1 after D1 contribution')
          CALL CC_PRP(THETA1,THETA2,ISYRES,1,0)
          XNORM = DNRM2(NT1AM(ISYRES),THETA1,1)
          WRITE (LUPRI,*) 'Norm of Theta(+D1)_ai:', XNORM
        END IF

      END IF
 
*---------------------------------------------------------------------*
*     From the BF(A) and BF(QA) calculate:
*   -  the E2"_ai contribution (in 2 shots) 
*   -  the GA and GQA intermediates (2BFA-BFA)
*   -  the tilde-Gamma(QA) intermediate for AG2_aibj (in 2 shots)
*      (it is assumed that the 2 contributions to E2" and Gamma are
*       automatically cumulated in THETA1 and KGAMQA.
*---------------------------------------------------------------------*
* Allocate for Zeta2(squared), calG's and Gamma(tilde,QA)

      IF ((LRELAX.OR.LTWOEL).AND.(.NOT.(CCS.OR.CC2))) THEN
        KZETA2 = 1
        KEND0  = KZETA2 + NT2SQ(ISYCTR)
        LWRK0  = LWORK  - KEND0
c
        KGAMQA = KEND0
        KGAIM  = KGAMQA + NGAMMA(ISYHTA)
        KGQAI  = KGAIM  + NT1AO(ISYMTA)
        KEND1  = KGQAI  + NT1AO(ISYHTA)

        LWRK1  = LWORK - KEND1
        IF (LWRK1 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CCFBTAF. (CTR2)')
        END IF
* initialize
        CALL DZERO(WORK(KGAIM),NT1AO(ISYMTA))
        CALL DZERO(WORK(KGQAI),NT1AO(ISYHTA))
        CALL DZERO(WORK(KGAMQA),NGAMMA(ISYHTA))
* relaxed case
c      IF ( (.NOT. (CCS .OR. CC2)) .AND. (LTWOEL.OR.LRELAX)) THEN

        LENBFA   = NT2AOIJ(ISYMTA)
        LENBFQA  = NT2AOIJ(ISYHTA)

        KBFA    = KEND1
        KEND2   = KBFA   + MAX(LENBFA,LENBFQA)
        LWRK2   = LWORK  - KEND2

        IF (LWRK2 .LT. NT2AM(ISYCTR)) THEN
          CALL QUIT('Insufficient work space in CCFBTAF. (CTR2)')
        END IF

C       --------------------------------------------------
C       read lagrange multipliers Zeta_2 from file and square up:
C       --------------------------------------------------
        IOPT = 2
        CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,WORK(KDUM),
     &                WORK(KEND2))

        CALL CC_T2SQ(WORK(KEND2),WORK(KZETA2),ISYCTR)

C       ---------------------------------
C       read BF(A) intermediate from file:
C       ---------------------------------
        CALL GETWA2(LUBFA,FILBFA,WORK(KBFA),IBFA,LENBFA)
        IF (LOCDBG) THEN
          XNORM = DNRM2(LENBFA,WORK(KBFA),1)
          WRITE (LUPRI,*) 'CCFBTAF> Norm of portion BFA', XNORM
        END IF

C       ------------------------------------------------------
C       transform to MO representation using XLAMDPQ matrices:
C       ------------------------------------------------------
        ICON    = 3
        IOPTG   = 2
        LGAMMA  = .TRUE.
        LO3BF   = .TRUE.                    !use BF with 3MO
        LBFZETA = .FALSE.

        CALL CC_T2MO3(THETA1,WORK(KZETA2),ISYCTR,WORK(KBFA),
     &                DUM,WORK(KGAMQA),WORK(KGAIM),DUM,
     &                XLAMDQP,ISYHOP,XLAMDQP,ISYHOP,
     &                WORK(KEND2),LWRK2,ISYMTA,
     &                ICON,LGAMMA,IOPTG,LO3BF,LBFZETA)

        IF (LOCDBG) THEN
           CALL AROUND('THETA1 after E2 contributions (part 1):')
           CALL CC_PRP(THETA1,WORK,ISYRES,1,0)
           XNORM = DNRM2(NT1AM(ISYRES),THETA1,1)
           WRITE (LUPRI,*) 'Norm of Theta(+E2,1)_ai:', XNORM
c           WRITE (LUPRI,*) 'the GAIM intermediate:'
c           WRITE (LUPRI,'(5F12.6)'),(WORK(KGAIM-1+I),I=1,NT1AO(ISYMTA))
c           WRITE (LUPRI,*) 'Norm(GAMQA) = ',
c     &       DDOT(NGAMMA(ISYHTA),WORK(KGAMQA),1,WORK(KGAMQA),1)
        END IF


**        IF ( LRELAX ) THEN            !I don't think we need this distinction
*
C       ---------------------------------
C       read BF(QA) intermediate from file:
C       ---------------------------------
           CALL GETWA2(LUBFQA,FILBFQA,WORK(KBFA),IBFQA,LENBFQA)   
           IF (LOCDBG) THEN
             XNORM = DNRM2(LENBFQA,WORK(KBFA),1)
             WRITE (LUPRI,*) 'Norm of portion BFQA', XNORM
           END IF


C          ------------------------------------
C          transform to MO representation using 
C          XLAMDP0 and XLAMDQP matrices:
C          ------------------------------------
           ICON    = 3                             
           IOPTG   = 2
           LGAMMA  = .TRUE.
           LO3BF   = .TRUE.
           LBFZETA = .FALSE.
           CALL CC_T2MO3(THETA1,WORK(KZETA2),ISYCTR,WORK(KBFA),
     &                     DUM,WORK(KGAMQA),WORK(KGQAI),DUM,
     &                     XLAMDP0,ISYM0,XLAMDP0,ISYM0,
     &                     WORK(KEND2),LWRK2,ISYHTA,
     &                     ICON,LGAMMA,IOPTG,LO3BF,LBFZETA)

**        END IF
        IF (LOCDBG) THEN
          CALL AROUND('FBTAF: THETA1 after E2 (sec) term contrib.s:')
          CALL CC_PRP(THETA1,THETA2,ISYRES,1,0)
          XNORM = DNRM2(NT1AM(ISYRES),THETA1,1)
          WRITE (LUPRI,*) 'Norm of Theta(+E2,2)_ai:', XNORM
        END IF

      END IF
*
*---------------------------------------------------------------------*
*     Calculate calE1* intermediate from h(1)_C^,a RA and RQA
*               calE2  intermediate from h(1)_i,K^ GA and GQA
*    ------->   (both relaxed/2e and unrelaxed/oneel
*               B2_ai from Zeta_1, calE1* and calE2 
*               (TGammaQA-delta*calE2) for AG2_aibj term
*---------------------------------------------------------------------*
*
      KEMAT2  = KEND1
      KEMAT1  = KEMAT2  + NMATIJ(ISYHTA)
      KEND2   = KEMAT1  + NMATAB(ISYHTA)
      LWRK2   = LWORK   - KEND2
*
      IF (LWRK2 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_FBTA. (B2)')
      END IF
*
* relaxed 2-electron operator
*
      IF ( (LTWOEL.OR.LRELAX) ) THEN
         IF (.NOT.CCS) THEN
            LRCON  = .TRUE.                      ! include contrib. from R 
            LGCON  = .TRUE.                      ! include contrib. from G
            FCKCON = .FALSE.                     ! esclude contrib. from Oneham
            IOPT   =   2
            CALL CC_EIM(WORK(KEMAT1),WORK(KEMAT2),
     &                 RAIM,RQAIM,WORK(KGAIM),WORK(KGQAI),DUM,DUM,
     &                 XLAMDH0,XLAMDP0,ISYM0,
     &                 XLAMDQH,XLAMDQP,ISYHOP,
     &                 FCKCON,LRCON,LGCON,LRELAX,IOPT,ISYHTA)

            LE1   = .TRUE.             !calculate and add Oneh^ta to E1* and E2
            LE2   = .TRUE.                       
            CALL CC_HPQTA(ONEHB,ISYHOP,WORK(KEMAT1),WORK(KEMAT2),ISYHTA,
     &                    LISTR,IDLSTR,ISYMTA,LE1,LE2,WORK(KEND2),LWRK2)
*
* CCS relaxed case: the full E1 and E2 intermediates (FOCKB. NB: FOCKB_ij is WRONG)
*
         ELSE
            LRCON  = .FALSE.                      ! esclude contrib. from R
            LGCON  = .FALSE.                      ! esclude contrib. from G
            FCKCON = .TRUE.                       ! include contrib. from FOCKB
            IOPT   =   1
            CALL CC_EIM(WORK(KEMAT1),WORK(KEMAT2),
     &                 DUM,DUM,DUM,DUM,FOCKB,FOCKB,
     &                 DUM,DUM,IDUM,
     &                 DUM,DUM,IDUM,
     &                 FCKCON,LRCON,LGCON,LRELAX,IOPT,ISYHTA)
         END IF

      ELSE

          CALL DZERO(WORK(KEMAT2),NMATIJ(ISYHTA))
          CALL DZERO(WORK(KEMAT1),NMATAB(ISYHTA))
          LE1   = .TRUE.        !calculate and add Oneh^ta to both E1* and E2
          LE2   = .TRUE.
          CALL CC_HPQTA(ONEHB,ISYHOP,WORK(KEMAT1),WORK(KEMAT2),ISYHTA,
     &                  LISTR,IDLSTR,ISYMTA,LE1,LE2,WORK(KEND2),LWRK2)
*
      END IF

      IF (LOCDBG) THEN
         CALL AROUND('E1* and E2 intermediates:')
         CALL CC_PREI(WORK(KEMAT1),WORK(KEMAT2),ISYRES,1)
         CALL FLSHFO(LUPRI)
      END IF
*
*-----------------------------------------------------------------*
* Calculated the B2 contribution to RES from E1* and E2 and Zeta1
* in both relaxed/2-e and unrel/1-e
* skip if unrel/1-e and CCS and L0 ?
*-----------------------------------------------------------------*
*

      CALL CCLR_E1C1(THETA1,ZETA1,WORK(KEMAT1),WORK(KEMAT2),
     &               WORK(KEND2),LWRK2,ISYCTR,ISYHTA,'T')


      IF (LOCDBG) THEN
          CALL AROUND('FBTAF: THETA1 after B2 contribution:')
          CALL CC_PRP(THETA1,THETA2,ISYRES,1,0)
          XNORM = DNRM2(NT1AM(ISYRES),THETA1,1)
          WRITE (LUPRI,*) 'Norm of Theta1(+B2) :',XNORM
          CALL FLSHFO(LUPRI)
      END IF

*---------------------------------------------------------------------*
*     for non-orbital-relaxed one-electron perturbations calculate
*     here the contribution of the E intermediates to the THETA2,bj
*     (dvs A and G contributions)
*---------------------------------------------------------------------*
      IF ( (.NOT.(LTWOEL.OR.LRELAX)) .AND. (.NOT.CCS) ) THEN

         IF (LOCDBG)
     &     WRITE (LUPRI,*) 'Unrelaxed 1-e oper B: A and G contributions'
         KZETA2 = KEND2
         KCTR2  = KZETA2 + NT2SQ(ISYCTR)
         KEND3  = KCTR2  + NT2AM(ISYCTR)
         LWRK3  = LWORK  - KEND3
   
         IF (LWRK3 .LT. 0) THEN
             CALL QUIT('Insufficient work space in CCFBTAF (A+G unrel)')
         END IF

         IOPT = 2
         CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,WORK(KDUM),
     &                 WORK(KCTR2))
         CALL CC_T2SQ(WORK(KCTR2),WORK(KZETA2),ISYCTR)

         CALL CC_EITR(WORK(KEMAT1),WORK(KEMAT2),
     &                                  WORK(KEND3),LWRK3,ISYHOP)

         CALL CCRHS_E(THETA2,WORK(KZETA2),WORK(KEMAT1),WORK(KEMAT2),
     &                WORK(KEND3),LWRK3,ISYCTR,ISYHOP)

         IF (LOCDBG) THEN
            WRITE (LUPRI,*) 
     &         'E term contribution from CCRHS_E: (A+G unrel)'
            CALL CC_PRP(THETA1,THETA2,ISYRES,0,1)
         END IF

      END IF
*
*-------------------------------------------------------------------*
* scale diagonal of Tilde-Gamma-QA with calE2 for later AG2_aibj
*                   (TGammaQA-delta*calE2) 
*-------------------------------------------------------------------*
*
      IF ( (.NOT. (CCS .OR. CC2)) .AND. (LTWOEL.OR.LRELAX)) THEN
         CALL CC_GAMMA2(WORK(KGAMQA),WORK(KEMAT2),ISYHTA)
      ENDIF
*                                             Finished with E1* and E2
*---------------------------------------------------------------------*
* calculate AG2-term contribution:
*       requires:  -- Gamma intermediate in core (KGAMQA)
*                  -- ZETA2 amplitudes squared in core (read them in)
*---------------------------------------------------------------------*
*
      IF ( (.NOT. (CCS .OR. CC2)) .AND. (LTWOEL.OR.LRELAX)) THEN

         KZETA2 = KGAIM                                    !reuse GAIM space
         KEND3  = KZETA2 + NT2SQ(ISYCTR)
         LWRK3  = LWORK - KEND3
c
         IF (LWRK3 .LT. NT2AM(ISYCTR)) THEN
            CALL QUIT('Insufficient work space in CCFBTAF. (AG2)')
         END IF

* read in Zeta_2 and square up

         IOPT = 2
         CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,WORK(KDUM),
     &                                              WORK(KEND3))
         CALL CC_T2SQ(WORK(KEND3),WORK(KZETA2),ISYCTR)
*
* calculate the A-term:
         IOPT = 2
         CALL CCRHS_A(THETA2,WORK(KZETA2),WORK(KGAMQA),
     *                WORK(KEND3),LWRK3,ISYHTA,ISYCTR,IOPT)
c
         IF (LOCDBG) THEN
            CALL AROUND('CC_FBTAF> THETA after AG2-term section:')
            CALL CC_PRP(THETA1,THETA2,ISYRES,1,1)
            XNORM = DNRM2(NT1AM(ISYRES),THETA1,1)
            WRITE (LUPRI,*) 'Norm of Theta_ai:', XNORM
            XNORM = DNRM2(NT2AM(ISYRES),THETA2,1)
            WRITE (LUPRI,*) 'Norm of Theta(AG2)_aibj:', XNORM

c            WRITE (LUPRI,*) 'ZETA2:'
c            CALL CC_PRSQ(WORK,WORK(KZETA2),ISYRES,0,1)
c            CALL FLSHFO(LUPRI)
c            WRITE (LUPRI,*) 'GAMMAQA:',(WORK(KGAMQA-1+I),I=1,NGAMMA(ISYHTA))
         END IF

      END IF
*
*---------------------------------------------------------------------*
*       calculate BFZeta-term (BEHF1)  and GZeta intermediate:
*---------------------------------------------------------------------*
*
      IF ( (.NOT. (CCS .OR. CC2)) .AND. (LTWOEL.OR.LRELAX)) THEN

         KGZAIM = KEND0                          !reuse/redefine G vectors
         KGZQAI = KGZAIM + NT1AO(ISYZTA)
         KEND1  = KGZQAI + NT1AO(ISYRES)
         KBFZA  = KEND1
         KBFZQA = KBFZA  + NT2AO(ISYZTA)
         KEND1  = KBFZQA + NT2AO(ISYRES)
         LWRK1  = LWORK  - KEND1

         CALL DZERO(WORK(KGZQAI),NT1AO(ISYRES))               
         CALL DZERO(WORK(KGZAIM),NT1AO(ISYZTA))
*
* Read BFZA and BFZQA from file and resort from varrho_beta,ij(delta)
* to varrho_i del, j bet (store del i, bet j)
*
         CALL CC_BFIFSORT(WORK(KBFZA),ISYZTA,LUBFZA,FILBFZA,IBFZA,
     &                                          WORK(KEND1),LWRK1)

         CALL CC_BFIFSORT(WORK(KBFZQA),ISYRES,LUBFZQA,FILBFZQA,IBFZQA,
     &                                              WORK(KEND1),LWRK1)
*
* transform BFZQA intermediate to MO representation using 
* XLAMDH0 matrices (first contribution to THETA2 and to GZQA):
*
         ICON   = 1
         IOPTG  = 1
         LGAMMA  = .FALSE.
         LO3BF   = .FALSE.
         LBFZETA = .TRUE.

         ISYMBF = ISYRES
         CALL CC_T2MO3(DUM,DUM,1,WORK(KBFZQA),THETA2,                   
     &                DUM,WORK(KGZQAI),DUM,
     &                XLAMDH0,ISYM0,XLAMDH0,ISYM0,WORK(KEND1),LWRK1,
     &                ISYMBF,ICON,LGAMMA,IOPTG,LO3BF,LBFZETA)
   
         IF ( LRELAX ) THEN

            ICON   = 4
            ISYMBF = ISYZTA
            CALL CC_T2MO3(DUM,DUM,1,WORK(KBFZA),THETA2,                
     &                  DUM,WORK(KGZAIM),WORK(KGZQAI),
     &                  XLAMDH0,ISYM0,XLAMDQH,ISYHOP,
     &                  WORK(KEND1),LWRK1,ISYMBF,
     &                  ICON,LGAMMA,IOPTG,LO3BF,LBFZETA)

        END IF
c           
         IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'Norm of THETA1 after B+E term:',
     &               DNRM2(NT1AM(ISYRES),THETA1,1)
            IF (.NOT.CCS)
     &        WRITE (LUPRI,*) 'Norm of THETA2 after BEHF1 term:',
     &                 DNRM2(NT2AM(ISYRES),THETA2,1)
           CALL AROUND('THETA2 after BEHF1 term:(B+E)')
C           CALL CC_PRP(THETA1,THETA2,ISYRES,0,1)
           CALL FLSHFO(LUPRI)
         END IF


      END IF
*
*---------------------------------------------------------------------*
*     calculate (E1')'' term from GZeta intermediate:
*     can also be done with 2 calls to CC_T1AM
*---------------------------------------------------------------------*
*
      IF ( (LTWOEL.OR.LRELAX) .AND. .NOT.(CCS.OR.CC2) ) THEN

         LGCON  = .TRUE.
         LRCON  = .FALSE.
         FCKCON = .FALSE.
         LTRSPF = .FALSE.
         IOPT   =   2
         CALL CC_GHJ(THETA1,DUM,DUM,WORK(KGZAIM),WORK(KGZQAI),
     *               DUM,
     *               XLAMDP0,XLAMDH0,ISYM0,
     *               XLAMDQP,XLAMDQH,ISYHOP,
     *               FCKCON,LRCON,LGCON,LRELAX,LTRSPF,IOPT,ISYRES)
*
*          note that use of hole & particle matrices is exchanged
*          compared to routine description (XI vector)

         IF (LOCDBG) THEN
           CALL AROUND("THETA1 after (E1')'' term:")
           CALL CC_PRP(THETA1,THETA2,ISYRES,1,0)
           XNORM = DDOT(NT1AM(ISYRES),THETA1,1,THETA1,1)
           WRITE (LUPRI,*) 'Norm of THETA1 after E1 term:',XNORM
           CALL FLSHFO(LUPRI)
         END IF

      END IF
*====================================================================*
* Calculate E1A intermediate from RAIM,RQAIM,FOCKB
*--------------------------------------------------------------------*
*
      KEMAT1A = KEND1                           
      KEND2   = KEMAT1A + NMATAB(ISYHTA)
      LWRK2   = LWORK   - KEND2
*
      IF (LWRK2 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_FBTA. (E1A interm)')
      END IF


      IF ( (.NOT.CCS) .AND. (LTWOEL.OR.LRELAX) ) THEN

         LRCON   = .TRUE.
         FCKCON  = .TRUE.
         IOPT    = 2
c
         CALL CC_E1AIM(WORK(KEMAT1A),RAIM,RQAIM,FOCKB,
     *                 XLAMDH0,ISYM0,XLAMDQH,ISYHOP,
     *                 FCKCON,LRCON,LRELAX,IOPT,ISYHTA)
*
      ELSE
         CALL QUIT('E1A intermediate: Alternative not yet available')
      END IF
*
      IF (LOCDBG) THEN
c         CALL AROUND('E1A-intermediate out from CC_E1AIM:')
c         CALL CC_PREI(WORK(KEMAT1A),WORK(KEMAT2),ISYHTA,1)
      END IF
*
*--------------------------------------------------------------------*
* Calculate the IAJB(1) integrals (4MO) (and save on file)
*--------------------------------------------------------------------*
*
      IF ( (.NOT. (CCS)) .AND. (LRELAX.OR.LTWOEL)) THEN  !do always? 

         ISYIAJB = ISYHOP
         KXIAJB  = KEND2
         KEND3   = KXIAJB + NT2SQ(ISYIAJB)
         LWRK3   = LWORK  - KEND1
c
         IF (LWRK3 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CCFBTAAO. (IAJB)')
         END IF
c
         CALL DZERO(WORK(KXIAJB),NT2SQ(ISYIAJB))
c
         IOPT = 2
         LOCC = .FALSE.
         CALL CC_IAJB2(WORK(KXIAJB),ISYIAJB,IOPT,LRELAX,LOCC,.FALSE.,
     &                LU0IAJB,FN0IAJB,IT2DEL0,XLAMDH0,ISYM0,
     &                LU1IAJB,FN1IAJB,IT2DELB,XLAMDQH,ISYHOP,
     &                WORK(KEND3),LWRK3)

*
* write (ia|jb)-bar integrals back to file for later use:
* (append after (ia|j delta) integrals)
*
         CALL PUTWA2(LU0IAJB,FN0IAJB,WORK(KXIAJB),
     &                               IADRMO,NT2SQ(ISYIAJB)) 
         IF (LOCDBG) THEN
c            WRITE (LUPRI,*) 'CCFBTAF: (ia|jb)(1) integrals:(PUTWA2)'
c            CALL CC_PRSQ(DUM,WORK(KXIAJB),ISYRES,0,1)
         END IF

*      END IF
*---------------------------------------------------------------------*
*     calculate A2_ai contrib. and project. vs <HF| for <HF| contrib.
*     requires:  -- TA2 amplitudes packed, (ia|jb)(1) squared in core
*     (reread integrals from file avoided by now)
*---------------------------------------------------------------------* 
*      IF ( (.NOT.CCS) .AND. (LTWOEL.OR.LRELAX) ) THEN  

          KA2IM  = KEND3                        !A2 interm
          KA2RS  = KA2IM  + NT1AM(ISYZTA)       !A2 contrib
          KTA2AM = KA2RS  + NT2AM(ISYMTA)       !TA2 amplitudes (packed)
          KEND4  = KTA2AM + NT2AM(ISYMTA)
          LWRK4  = LWORK  - KEND4
          IF (LWRK4 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CCFBTA. (CTR2)')
          END IF                     
*         --------------------------------------------------
*         read double excitation response amplitudes,
*         construct 2TA2 - TA2 in packed form
*         --------------------------------------------------
          IOPT = 2
          CALL CC_RDRSP(LISTR,IDLSTR,ISYMTA,IOPT,MODEL,
     &                               WORK(KDUM),WORK(KTA2AM))
          CALL CCLR_DIASCL(WORK(KTA2AM),TWO,ISYMTA)
          IOPTTCME = 1
          CALL CCSD_TCMEPK(WORK(KTA2AM),ONE,ISYMTA,IOPTTCME)

*         -----------------------------------------
*         contract: A2IM_ai = sum_bj Zeta1_bj (2TA-TA)_aibj
*         -----------------------------------------
          IOPT = 0                              !packed form
          CALL CCG_LXD(WORK(KA2IM),ISYZTA,ZETA1,ISYCTR,
     &                             WORK(KTA2AM),ISYMTA,IOPT)  
          IF (LOCDBG) THEN
             WRITE (LUPRI,*) 'CCFBTAF: The A2IM_ai intermediate for A2'
             CALL CC_PRP(WORK(KA2IM),WORK(KDUM),ISYZTA,1,0)
          END IF
*         ----------------------------------------------------------
*         calculate L(ia,jb)(1) from the (ia|jb)(1) integrals
*         done with squared integrals, (ia|jb)(1) overwritten
*         Checked input-output indices of L(ia,jb)!!!!!!!!!!!!!!!!!
*         ----------------------------------------------------------
c
c          CALL GETWA2(LU0IAJB,FN0IAJB,WORK(KXIAJB),
c     &                                IADRMO,NT2SQ(ISYHOP)) 
c
          KEND4 = KTA2AM
          LWRK4 = LWORK  - KEND4

          CALL CCRHS_T2TR(WORK(KXIAJB),WORK(KEND4),LWRK4,ISYHOP)
*
*         -------------------------------------------
*         contract: theta_ai = sum_bj L(ia,jb) A2IM(jb) --> RES
*         -------------------------------------------
          IOPT = 1                                    !squared
          CALL CCG_LXD(WORK(KA2RS),ISYRES,WORK(KA2IM),ISYZTA,
     &                               WORK(KXIAJB),ISYHOP,IOPT)
          
          CALL DAXPY(NT1AM(ISYRES),ONE,WORK(KA2RS),1,THETA1,1)  

          IF (LOCDBG) THEN
              CALL AROUND('CCFBTAF: THETA1 after A2 contribution:')
              CALL CC_PRP(THETA1,THETA2,ISYRES,1,0)
              XNORM = DNRM2(NT1AM(ISYRES),THETA1,1)
              WRITE (LUPRI,*) 
     *           'Norm of THETA1 after A2 contribution:',XNORM
          END IF
*
      END IF 
*
*---------------------------------------------------------------------*
*       calculate C-term contribution:
*         requires:  -- (ia|jb)(1) integrals squared in core (X_ai,bj)
*                    -- Special (ij|cb)(1) integrals
*                    -- TA2 amplitudes packed in core
*                    -- E1A intermediate
*                    -- C-tilde intermediate in core (squared matrix)
*                    -- result vector in core (packed)
*---------------------------------------------------------------------*
      IF ( (.NOT. (CCS .OR. CC2)) .AND. (LTWOEL.OR.LRELAX)) THEN

         KCTILD = KEND2
         KXBICJ = KCTILD + NT2SQ(ISYHTA)      !must be initialized outside call
         KXIJCB = KXBICJ + NT2SQ(ISYHTA)
         KEND2  = KXIJCB + NT2SQ(ISYHTA)
         LWRK2  = LWORK  - KEND2
*  -----------------------------------------------------------
*  read special (ij|c del) integrals into memory (stored cj,i,delta)
*  transform to (ij|cb) (stored cj,bi) using XLAMDH0/XLAMDQH
*  resort to Y_bi,cj, add to CTILDE interm.
*  -----------------------------------------------------------
       
         CALL DZERO(WORK(KXIJCB),NT2SQ(ISYHTA))

         IOPT = 2
         LOCC = .FALSE.
         CALL CC_IAJB2(WORK(KXIJCB),ISYHTA,IOPT,LRELAX,LOCC,.FALSE.,
     &                  LU0IJCB,FN0IJCB,IT2DELA0,XLAMDH0,ISYM0,
     &                  LU1IJCB,FN1IJCB,IT2DELAB,XLAMDQH,ISYHOP,
     &                  WORK(KEND2),LWRK2)
         IF (LOCDBG) THEN
c            WRITE (LUPRI,*) 'CCFBTAF:(ij|cb)(1) special integrals:(cj,bi)'
c            CALL CC_PRSQ(DUM,WORK(KXIJCB),ISYHTA,0,1)
         END IF

* transpose to (bi,cj)

         DO ISYMCJ = 1, NSYM
            ISYMBI = MULD2H(ISYMCJ,ISYHTA)
            NDIM1  = NT1AM(ISYMCJ)
            NDIM2  = NT1AM(ISYMBI)
            KOFF1  = KXIJCB + IT2SQ(ISYMCJ,ISYMBI)
            KOFF2  = KXBICJ + IT2SQ(ISYMBI,ISYMCJ)

            CALL TRSREC(NDIM1,NDIM2,WORK(KOFF1),WORK(KOFF2))
         
         END DO
         IF (LOCDBG) THEN
c            WRITE (LUPRI,*) 'CCFBTAF:(ij|cb)(1) special integrals:(bi,cj)'
c            CALL CC_PRSQ(DUM,WORK(KXBICJ),ISYHTA,0,1)
         END IF

* second contribution (reuse KXIJCB)

         KXIAJB = KXIJCB
         KTA2AM = KXIAJB + NT2SQ(ISYIAJB)
         KEND2  = KTA2AM + NT2AM(ISYMTA)
         LWRK2  = LWORK  - KEND2
*
         IF (LWRK2 .LT. 0) THEN
           CALL QUIT('Insufficient work space in CCFBTAF. (C tilde)')
         END IF
* read amplitudes from file and scale with 2 (pck):
         IOPT = 2
         CALL CC_RDRSP(LISTR,IDLSTR,ISYMTA,IOPT,MODEL,
     &                              WORK(KDUM),WORK(KTA2AM)) 
         CALL CCLR_DIASCL(WORK(KTA2AM),TWO,ISYMTA)
* read 4MO (ia|jb)(1) integrals from file (sqr):
         CALL GETWA2(LU0IAJB,FN0IAJB,WORK(KXIAJB),
     &                               IADRMO,NT2SQ(ISYHOP))
         IF (LOCDBG) THEN
c            WRITE (LUPRI,*) 'CCFBTAF: (ia|jb)(1) integrals: (GETWA2)'
c            CALL CC_PRSQ(DUM,WORK(KXIAJB),ISYRES,0,1)
         END IF
*
* contract 4MO (ia|jb)(1) integrals with TA2 amplitudes and add to KCTILD
*
         ISYIAJB = ISYHOP
         IOPT = 1                                        !result in memory
         CALL CCB_CDBAR('C',WORK(KXIAJB),ISYIAJB, WORK(KTA2AM),ISYMTA,
     &                      WORK(KCTILD),ISYHTA, WORK(KEND2),LWRK2,
     &                      CDUM, IDUM, IDUM, IOPT)
*
*         ------------------------------------------------------------
*         add E1A intermediate to the diagonal of the C intermediate:
*         ------------------------------------------------------------
*
         FAC = -0.5D0
         CALL CC_CDBAR2(WORK(KCTILD),WORK(KEMAT1A),FAC,.FALSE.,ISYHTA)

*        ------------------------------------------------------------
*        add special integrals
*        ------------------------------------------------------------

         CALL DAXPY(NT2SQ(ISYHTA),ONE,WORK(KXBICJ),1,WORK(KCTILD),1)
*
* Read in Zeta_2 (packed)
*
         KZETA2 = KXBICJ
         KEND2  = KZETA2 + NT2AM(ISYCTR)
         LWRK2  = LWORK  - KEND2
*
         IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_FBTAF. (CTR2)')
         END IF                                                
* read lagrange multipliers from file:
         IOPT = 2
         CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,WORK(KDUM),
     &                                             WORK(KZETA2))
*
          IF (LOCDBG) THEN
             WRITE (LUPRI,*) 'KZETA2:',KZETA2
             WRITE (LUPRI,*) 'ISYCTR:',ISYCTR
             WRITE (LUPRI,*) 'ZETA2:'
c             CALL CC_PRP(WORK,WORK(KZETA2),ISYCTR,0,1)
             CALL FLSHFO(LUPRI)
          END IF  
*
*         ------------------------------------------
*         calculate C term and add to result vector:
*         ------------------------------------------
          ioptr12 = 0
          CALL CC_CD('C',-1,ioptr12,THETA2,ISYRES,WORK(KZETA2),
     &               ISYCTR,WORK(KCTILD),ISYHTA,WORK(KEND2),LWRK2)

          IF (LOCDBG) THEN
             CALL AROUND ('FBTAF: THETA2 after the C contribution:')
c             CALL CC_PRP(THETA1,THETA2,ISYRES,1,1)
             WRITE (LUPRI,*) 'Norm of Theta_ai: ', 
     &                     DNRM2(NT1AM(ISYRES),THETA1,1) 
             WRITE (LUPRI,*) 'Norm of Theta_aibj: ', 
     &                     DNRM2(NT2AM(ISYRES),THETA2,1) 
             CALL FLSHFO(LUPRI)
          END IF            

      END IF                            
*-----------------------------------------------------------*
* D term
*-----------------------------------------------------------*
      IF ( (.NOT. (CCS .OR. CC2)) .AND. (LTWOEL.OR.LRELAX)) THEN

         KDTILD = KCTILD
         KXBICJ = KDTILD + NT2SQ(ISYHTA)
         KXCJIB = KXBICJ + NT2SQ(ISYHTA)
         KEND2  = KXIJCB + NT2SQ(ISYHTA)
         LWRK2  = LWORK  - KEND2

         CALL DZERO(WORK(KXCJIB),NT2SQ(ISYHTA))
         CALL DZERO(WORK(KXBICJ),NT2SQ(ISYHTA))
*  -----------------------------------------------------------
*  read special L(cj|i del) integrals into memory, transform to
*  (cj|bi) using XLAMDH0/QH, resort to Y_bi,cj, add to DBAR interm.
*  -----------------------------------------------------------
         IOPT = 2
         LOCC = .FALSE.
         CALL CC_IAJB2(WORK(KXCJIB),ISYHTA,IOPT,LRELAX,LOCC,.FALSE.,
     &                  LU0CJIB,FN0CJIB,IT2DELA0,XLAMDH0,ISYM0,
     &                  LU1CJIB,FN1CJIB,IT2DELAB,XLAMDQH,ISYHOP,
     &                  WORK(KEND2),LWRK2)
         IF (LOCDBG) THEN
c            WRITE (LUPRI,*) 'CCFBTAF:L(cj|ib)(1) special integrals:(cj,bi)'
c            CALL CC_PRSQ(DUM,WORK(KXCJIB),ISYRES,0,1)
         END IF
* transpose integrals
         DO ISYMCJ = 1, NSYM
            ISYMBI = MULD2H(ISYMCJ,ISYHTA)
            NDIM1  = NT1AM(ISYMCJ)
            NDIM2  = NT1AM(ISYMBI)
            KOFF1  = KXCJIB + IT2SQ(ISYMCJ,ISYMBI)
            KOFF2  = KXBICJ + IT2SQ(ISYMBI,ISYMCJ)

            CALL TRSREC(NDIM1,NDIM2,WORK(KOFF1),WORK(KOFF2))

         END DO
*
         IF (LOCDBG) THEN
c            WRITE (LUPRI,*) 'CCFBTAF:L(cj|ib)(1) special integrals:(bi,cj)'
c            CALL CC_PRSQ(DUM,WORK(KXBICJ),ISYRES,0,1)
         END IF
*
         ISYIAJB = ISYHOP
         KXIAJB  = KXCJIB
         KTA2AM  = KXIAJB + NT2SQ(ISYIAJB)
         KEND2   = KTA2AM + NT2AM(ISYMTA)
         LWRK2   = LWORK - KEND2
 
* read amplitudes from file and scale with 2 (pck):
         IOPT = 2
         CALL CC_RDRSP(LISTR,IDLSTR,ISYMTA,IOPT,MODEL,
     &                              WORK(KDUM),WORK(KTA2AM))
         CALL CCLR_DIASCL(WORK(KTA2AM),TWO,ISYMTA)

* read 4MO (ia|jb)(1) integrals from file:
         CALL GETWA2(LU0IAJB,FN0IAJB,WORK(KXIAJB),
     &                              IADRMO,NT2SQ(ISYIAJB))
*
* contract (ia|jb)-bar integrals with T2 amplitudes:
* (this overwrites T2 with 2T2(ia;jb) - T2(ja;ib)...)
         IOPT = 1
         CALL CCB_CDBAR('D',WORK(KXIAJB),ISYIAJB, WORK(KTA2AM),ISYMTA,
     &                      WORK(KDTILD),ISYHTA, WORK(KEND2),LWRK2,
     &                      CDUM, IDUM, IDUM, IOPT)
*                                                                    
*        ------------------------------------------------------------
*        add E1A intermediate to the diagonal of the D intermediate:
*        ------------------------------------------------------------
*
         FAC = 0.5D0
         CALL CC_CDBAR2(WORK(KDTILD),WORK(KEMAT1A),FAC,.FALSE.,ISYHTA)

*        ------------------------------------------------------------
*        add special integrals
*        ------------------------------------------------------------

         CALL DAXPY(NT2SQ(ISYHTA),ONE,WORK(KXBICJ),1,WORK(KDTILD),1)

*         ------------------------------------------------------------
*         calculate D term and add to result vector:
*         ------------------------------------------------------------
* reread in Zeta_2
*
         KZETA2 = KXBICJ
         KEND2  = KZETA2 + NT2AM(ISYCTR)
         LWRK2  = LWORK  - KEND2
*
         IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_FBTAF. (CTR2)')
         END IF
* read lagrange multipliers from file:
         IOPT = 2
         CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,WORK(KDUM),
     &                                               WORK(KZETA2))
c
          IF (LOCDBG) THEN
c             WRITE (LUPRI,*) 'KZETA2:',KZETA2
c             WRITE (LUPRI,*) 'ISYCTR:',ISYCTR
c             WRITE (LUPRI,*) 'ZETA2:'
c             CALL CC_PRP(WORK,WORK(KZETA2),ISYCTR,0,1)
c             CALL FLSHFO(LUPRI)
          END IF
*
          ioptr12 = 0
          CALL CC_CD('D',-1,ioptr12,THETA2,ISYRES,WORK(KZETA2),ISYCTR,
     &                      WORK(KDTILD),ISYHTA,WORK(KEND2),LWRK2)
*
          IF (LOCDBG) THEN
             CALL AROUND('FBTAF: Finished D-term section...')
             WRITE (LUPRI,*) 'THETA2:'
c             CALL CC_PRP(WORK,THETA2,ISYRES,0,1)
             WRITE (LUPRI,*) 'Norm of Theta_ai: ', 
     &                     DNRM2(NT1AM(ISYRES),THETA1,1) 
             WRITE (LUPRI,*) 'Norm of Theta_aibj: ', 
     &                     DNRM2(NT2AM(ISYRES),THETA2,1) 
             CALL FLSHFO(LUPRI)
          END IF
c
      END IF  
*                                                              
*---------------------------------------------------------------------*
*       (note that ZETA2 vector has been overwritten in C/D section)
*---------------------------------------------------------------------*

*---------------------------------------------------------------------*
*       calculate 21F and 21G contributions (the F and G terms):
*---------------------------------------------------------------------*
 
      IF ( (LTWOEL.OR.LRELAX) .AND. .NOT.(CCS .OR. CC2) ) THEN

        CALL CC_T1AM(THETA1,ISYRES,FAIM,ISYZTA,XLAMDQH,ISYHOP,ONE)
        CALL CC_T1AM(THETA1,ISYRES,FQAIM,ISYRES,XLAMDH0,ISYM0,ONE)

        IF (LOCDBG) THEN
          CALL AROUND ('The total F contrib. (21F) + E to THETA1')
          CALL CC_PRP(THETA1,WORK,ISYRES,1,0)
             WRITE (LUPRI,*) 'Norm of Theta_ai: ', 
     &                     DNRM2(NT1AM(ISYRES),THETA1,1) 
             WRITE (LUPRI,*) 'Norm of Theta_aibj: ', 
     &                     DNRM2(NT2AM(ISYRES),THETA2,1) 
          CALL FLSHFO(LUPRI)
        END IF
*
* The two G contributions 
*
        MT2BCD = 0
        DO ISYM = 1, NSYM
          MT2BCD = MAX(MT2BCD,NT2BCD(ISYM))   
        END DO
*
* 1) the G''
*
        ISYIFJL = ISYHOP
        KXIFJL  = KEND2
        KYJLIF  = KXIFJL + NTAIKJ(ISYIFJL)
        KPAINT  = KYJLIF + NTAIKJ(ISYIFJL)
        KQAINT  = KPAINT + MT2BCD
        KWAINT  = KQAINT + MT2BCD
        KEND3   = KWAINT + MT2BCD
        LWRK3   = LWORK  - KEND3
*
        IF (LWRK3 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_FBTAF. (21G-term)')
        END IF
*
*       -------------------------------------------------
*       calculate (if|jl) integrals from the 2 (if|jdelta), 
*       stored as X(fi,l,j), and resort to Y(jli,f):
*       integrals read from file inside CC_IAJB2
*       -------------------------------------------------
        CALL DZERO(WORK(KXIFJL),NTAIKJ(ISYIFJL))

        IOPT = 2              !also derivative/relaxed integrals
        LOCC = .TRUE.         !delta is transformed to occupied
        CALL CC_IAJB2(WORK(KXIFJL),ISYIFJL,IOPT,LRELAX,LOCC,.FALSE.,
     &                LU0IAJB,FN0IAJB,IT2DEL0,XLAMDH0,ISYM0,
     &                LU1IAJB,FN1IAJB,IT2DELB,XLAMDQH,ISYHOP,
     &                WORK(KEND3),LWRK3)

        CALL CCG_SORT1(WORK(KXIFJL),WORK(KYJLIF),ISYIFJL,5)

*       --------------------------------------------------------
*       contract with calP and calQ intermediates to G''
*       --------------------------------------------------------
        CALL DZERO(WORK(KWAINT),MT2BCD)

        DO ISYMF = 1, NSYM                           !sym f

           ISYALJ = MULD2H(ISYZTA,ISYMF)             !sym alj in PA,QA
           ISYJLI = MULD2H(ISYIFJL,ISYMF)            !sym jli in Y
*
           DO F = 1, NVIR(ISYMF)
              IVIRF  = IVIR(ISYMF) + F
              IADR   = IADPQAMO(IVIRF)       !get address for PA^f,QA^f on file
              LEN    = NT2BCD(ISYALJ)
              CALL GETWA2(LUPQAMO,FILPQAMO,WORK(KPAINT),IADR    ,LEN)   !get PA
              CALL GETWA2(LUPQAMO,FILPQAMO,WORK(KQAINT),IADR+LEN,LEN)   !get QA
*
              KOFF1 = ISJIKA(ISYJLI,ISYMF) + NMAIJK(ISYJLI)*(F-1)
     &                + KYJLIF

              CALL CC_21H(THETA1,ISYRES,WORK(KQAINT),
     &                    WORK(KWAINT),WORK(KPAINT),ISYALJ,
     &                    WORK(KOFF1),ISYIFJL,WORK(KEND3),LWRK3,ISYMF)
           END DO

        END DO 
*
        IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'Finished first 21G term section... THETA1:'
           CALL CC_PRP(THETA1,WORK,ISYRES,1,0)
             WRITE (LUPRI,*) 'Norm of Theta_ai: ', 
     &                     DNRM2(NT1AM(ISYRES),THETA1,1) 
             WRITE (LUPRI,*) 'Norm of Theta_aibj: ', 
     &                     DNRM2(NT2AM(ISYRES),THETA2,1) 
           CALL FLSHFO(LUPRI)
        END IF
*
* 2) the G'
*
        ISYIFJL = ISYHTA
        KXIFJL  = KEND2
        KYJLIF  = KXIFJL + NTAIKJ(ISYIFJL)
        KPINT0  = KYJLIF + NTAIKJ(ISYIFJL)
        KQINT0  = KPINT0 + MT2BCD
        KWINT0  = KQINT0 + MT2BCD
        KEND3   = KWINT0 + MT2BCD
        LWRK3   = LWORK  - KEND3

        IF (LWRK3 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_FBTAF. (21G-term)')
        END IF                                                    
*       -------------------------------------------------
*       calculate (if|j l^) integrals from the 2 (if|jdelta),
*       stored as X(fi,l^,j), and resort to Y(j(l^)i,f):
*       -------------------------------------------------
        CALL DZERO(WORK(KXIFJL),NTAIKJ(ISYIFJL))

        IOPT = 2
        LOCC = .TRUE.
        CALL CC_IAJB2(WORK(KXIFJL),ISYIFJL,IOPT,LRELAX,LOCC,.FALSE.,
     &                LU0IAJB,FN0IAJB,IT2DEL0,XLAMDAH,ISYMTA,
     &                LU1IAJB,FN1IAJB,IT2DELB,XLAMQAH,ISYHTA,
     &                WORK(KEND3),LWRK3)
*
        CALL CCG_SORT1(WORK(KXIFJL),WORK(KYJLIF),ISYIFJL,5)            
*       --------------------------------------------------------
*       contract with P and Q intermediates to G':
*       --------------------------------------------------------
        CALL DZERO(WORK(KWINT0),MT2BCD)

        DO ISYMF = 1, NSYM                           !sym f

           ISYALJ = MULD2H(ISYCTR,ISYMF)             !sym alj in P,Q
           ISYJLI = MULD2H(ISYIFJL,ISYMF)            !sym jli in Y

           DO F = 1, NVIR(ISYMF)
              IVIRF  = IVIR(ISYMF) + F
              IADR   = IADRPQMO(IVIRF)            !get address for P,Q on file
              LEN    = NT2BCD(ISYALJ)
              CALL GETWA2(LUPQMO,FILPQMO,WORK(KPINT0),IADR    ,LEN)
              CALL GETWA2(LUPQMO,FILPQMO,WORK(KQINT0),IADR+LEN,LEN)

              KOFF1 = KYJLIF + ISJIKA(ISYJLI,ISYMF)
     &                       + NMAIJK(ISYJLI)*(F-1)

              CALL CC_21H(THETA1,ISYRES,WORK(KQINT0),
     &                    WORK(KWINT0),WORK(KPINT0),ISYALJ,
     &                    WORK(KOFF1),ISYIFJL,WORK(KEND3),LWRK3,ISYMF)
           END DO

        END DO

        IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'Finished second 21G term section... THETA1:'
           CALL CC_PRP(THETA1,WORK,ISYRES,1,0)
           WRITE (LUPRI,*) 'Norm of Theta_ai: ', 
     &                     DNRM2(NT1AM(ISYRES),THETA1,1) 
           WRITE (LUPRI,*) 'Norm of Theta_aibj: ', 
     &                     DNRM2(NT2AM(ISYRES),THETA2,1) 
           CALL FLSHFO(LUPRI)
        END IF                       

      END IF
*---------------------------------------------------------------------*
*       calculate the E2' contribution to THETA1
*       requires:  the 3OCC integrals
*                  the IAJB integrals (squared)
*                  the t20 amplitudes (packed) --> passed outside ITRAN loop
*                  the t1a amplitudes
*                  the Zeta_2 vectors (packed o squared?)
*---------------------------------------------------------------------*
      IF ( (.NOT.CCS).AND.(LTWOEL.OR.LRELAX) ) THEN
         KGAMMAQ = 1                              !initialized in CC_MI
         KXIAJB  = KGAMMAQ + N3ORHF(ISYHOP)
         KEND1   = KXIAJB  + NT2SQ(ISYHOP)      
         LWRK1   = LWORK - KEND1
         IF (LWRK1 .LT. 0) THEN
           CALL QUIT(
     *       'Insufficient work space in CCFBTAF. (E2prim term)')
         END IF

* read 4MO (ia|jb)(1) integrals from file:
         CALL GETWA2(LU0IAJB,FN0IAJB,WORK(KXIAJB),
     &                              IADRMO,NT2SQ(ISYHOP))
* contract with T(0)_aibj to M-type intermediate:       
         CALL CC_MI(WORK(KGAMMAQ),WORK(KXIAJB),ISYHOP,T2AMP0,
     &                               ISYM0,WORK(KEND1),LWRK1)
         IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'Norm(GAMMABFQ) after CC_MI',
     &      DNRM2(N3ORHF(ISYHOP),WORK(KGAMMAQ),1)
         END IF
* calculate and add the contribution from the 4occ integrals
         LDERRLX = .FALSE.
         IF (LTWOEL.OR.LRELAX) LDERRLX = .TRUE.
         CALL CC_E21CON(WORK(KGAMMAQ),ISYHOP,
     &                XIJK0,ISYM0,LU1IJKL,FN1IJKL, 
     &                ISYHOP,XLAMDH0,ISYM0,XLAMDQH,ISYHOP,LDERRLX,
     &                ITRAN,WORK(KEND1),LWRK1)
         IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'Norm(GAMMABFQ) after CC_E21CON',
     &      DNRM2(N3ORHF(ISYHOP),WORK(KGAMMAQ),1)
         END IF
*
* read in ZA2 and square up,  read in TA1
*
         KZETA2 = KXIAJB
         KTA1AM = KZETA2 + NT2SQ(ISYCTR)
         KEND1  = KTA1AM + NT1AM(ISYMTA)
         LWRK1  = LWORK - KEND1
         IOPT = 2
         CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,WORK(KDUM),
     &                WORK(KEND1))
         CALL CC_T2SQ(WORK(KEND1),WORK(KZETA2),ISYCTR)
c
         IOPT = 1
         CALL CC_RDRSP(LISTR,IDLSTR,ISYMTA,IOPT,MODEL,
     &                             WORK(KTA1AM),WORK(KDUM))
*
* transform b to k of Z2, resort to jkl,a order (as in gamma), 
* matrix multiply and add  to result THETA1
*
         CALL CC_E22CON(WORK(KZETA2),ISYCTR,WORK(KTA1AM),ISYMTA,
     &         WORK(KGAMMAQ),ISYHOP,THETA1,ISYRES,WORK(KEND1),LWRK1)

         IF (LOCDBG) THEN
            CALL AROUND('THETA after E2primo-term:')
c            CALL CC_PRP(THETA1,THETA2,ISYRES,1,1)
             WRITE (LUPRI,*) 'Norm of Theta_ai: ', 
     &                     DNRM2(NT1AM(ISYRES),THETA1,1) 
             WRITE (LUPRI,*) 'Norm of Theta_aibj: ', 
     &                     DNRM2(NT2AM(ISYRES),THETA2,1) 
         END IF

       END IF 
*
*---------------------------------------------------------------------*
*       calculate the H contribution (HF2):
*---------------------------------------------------------------------*
*
      IF ((.NOT. CCS).AND.(LTWOEL.OR.LRELAX)) THEN

         KYPSA   = 1
         KYPSQA  = KYPSA  + NGLMDT(ISYZTA)
         KXPCK   = KYPSQA + NGLMDT(ISYRES)
         KXIAJB  = KXPCK  + NT2AM(ISYRES)          !the sum of the last 2 terms in HF2
         KEND3   = KXIAJB + NT2SQ(ISYRES)          
         LWRK3   = LWORK  - KEND3
*    
         IF (LWRK3 .LT. 0) THEN
           CALL QUIT('Insufficient work space in CCFBTAF. (H-term)')
         END IF
*
*        the (P^ab_ij) (2Coul-Exc) CTR1 x L(1)_jbkk-hat (in FOCKBA) contribution:
*
         CALL CC_L1FCK(THETA2,ZETA1,FOCKBA,ISYCTR,ISYHTA,
     &                                       WORK(KEND3),LWRK3)

         IF (LOCDBG) THEN
           CALL AROUND('CC_L1FCK (F2 of HF2) contribution to THETA2')
C           CALL CC_PRP(WORK,THETA2,ISYRES,0,1)
             WRITE (LUPRI,*) 'Norm of Theta_aibj: ', 
     &                     DNRM2(NT2AM(ISYRES),THETA2,1) 
         END IF

         IF (LTWOEL.OR.LRELAX) THEN
*
*           calculate the Ypsilon matrices:
*
            CALL CCLT_YPS1(ZETA1,ISYCTR,YAIM,ISYZTA,XLAMDAH,ISYMTA, 
     &                                   XLAMDH0,ISYM0,WORK(KYPSA))
            CALL CCLT_YPS1(ZETA1,ISYCTR,YAIM,ISYZTA,XLAMQAH,ISYHTA,
     &                                 XLAMDQH,ISYHOP,WORK(KYPSQA))
*
*          transform the (ia|j del) integrals with Ypsilon matrices:
*          the Ypsilon matrices are similar to Lambda matrices, so
*          we can use CC_IAJB2. Output array must be initialized here

            CALL DZERO(WORK(KXIAJB),NT2SQ(ISYRES))
            CALL DZERO(WORK(KXPCK),NT2AM(ISYRES))

            IOPT = 2
            LOCC = .FALSE.
            CALL CC_IAJB2(WORK(KXIAJB),ISYRES,IOPT,LRELAX,LOCC,.FALSE.,
     &                    LU0IAJB,FN0IAJB,IT2DEL0,WORK(KYPSA),ISYZTA,
     &                    LU1IAJB,FN1IAJB,IT2DELB,WORK(KYPSQA),ISYRES,
     &                    WORK(KEND3),LWRK3)
*
*           apply P^ab_ij, pack the integrals and add to result vector:
*           sort of 2 coulomb - exchange.....
*
            IOPT = 1
            CALL CC_T2PK(WORK(KXPCK),WORK(KXIAJB),ISYRES,1) 
            IOPTTCME = 1
            CALL CCSD_TCMEPK(WORK(KXPCK),ONE,ISYRES,IOPTTCME)
            !copy on RES vect
            CALL DAXPY(NT2AM(ISYRES),-ONE,WORK(KXPCK),1,THETA2,1)
c
            IF (LOCDBG) THEN
c               WRITE (LUPRI,*) 'L(jb|ia^Y) integrals in H-term:'
c               CALL CC_PRP(WORK,WORK(KXPCK),ISYRES,0,1)
               CALL AROUND('THETA after H(HF2)-term:')
c                 CALL CC_PRP(THETA1,THETA2,ISYRES,1,1)
               WRITE (LUPRI,*) 'Norm of Theta_aibj: ', 
     &                     DNRM2(NT2AM(ISYRES),THETA2,1) 
            END IF
c
         END IF
c
      END IF
c
*---------------------------------------------------------------------*
*       That's it; return
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
         CALL AROUND('final result vector in CCFBTAF ')
          IF (LFMATST) THEN
             CALL AROUND('Test vs FMAT: scale with 1/5')
             CALL DSCAL(NT1AM(ISYHOP),FIFTH,THETA1,1)
             IF (.NOT.CCS)
     &         CALL DSCAL(NT2AM(ISYHOP),FIFTH,THETA2,1)
          END IF
         CALL CC_PRP(THETA1,THETA2,ISYRES,1,1) 
      END IF
     
      CALL QEXIT('CCFBTAF')

      RETURN
      END 
*=====================================================================*
*                      END OF SUBROUTINE CCFBTAF                      *
*=====================================================================*
