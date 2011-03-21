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
c /* deck ccfbtaao1 */
*=====================================================================*
      SUBROUTINE CCFBTAAO1(X0INT,ISY0DIS,X1INT,ISY1DIS, 
     &                    BF0RHF,BF1RHF,               
     &                    DENS0,DPCK0,FOCK0,
     &                    DENSQ,DPCKQ,FOCKQ, 
     &                    DENSA,DPCKA,FOCKA,
     &                    DENSQA,DPCKQA,FOCKQA,
     &                    XLAMD0P,XLAMD0H,XLAMDQP,XLAMDQH,
     &                    XLAMDAP,XLAMDAH,XLAMDQAP,XLAMDQAH, 
     &                    RHOBFA,RHOBFQA,                   
     &                    LUBFDA,FNBFDA,IADRBFA,           
     &                    LUBFDQA,FNBFDQA,IADRBFQA,       
     &                    LU0IAJB,FN0IAJB,
     &                    LU1IAJB,FN1IAJB,
     &                    IT2DEL0,IADR0,IT2DELB,IADRB,
     &                    LU0IJCB,LU0CJIB,
     &                    FN0IJCB,FN0CJIB,
     &                    LU1IJCB,LU1CJIB,
     &                    FN1IJCB,FN1CJIB,
     &                    IT2DEL0A,IADR0A,
     &                    IT2DELBA,IADRBA,
     &                    IDEL,LZERO,LNEWTA,
     &                    LRELAX,LTWOEL,LX1ISQ,IREAL,
     &                    ISYHOP,ISYMTA,WORK, LWORK)
*
*---------------------------------------------------------------------*
*
*   Purpose: calculate intermediates for FbTa vector which depend on 
*            the AO integrals and at most TA and IOPER (No Zeta)
*
*     contrib. depending on X1INT/D1SRHF are only computed for LTWOEL ?
*     contrib. depending on X0INT/D0SRHF are only computed for LRELAX ?
*
*     (only exception: the (ia|jb), (ij|cb), (cj|ib) integrals)
*     zeroth-order MO integrals are only computed for LZERO
*
*   Written by Sonia Coriani, February 1999
*   Version: 08/10-1999
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

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

* variables:
      LOGICAL LZERO, LNEWTA, LRELAX, LTWOEL, LX1ISQ
      CHARACTER*(*) FNBFDA, FNBFDQA
      CHARACTER*(*) FN0IAJB,FN1IAJB
      CHARACTER*(*) FN0IJCB,FN0CJIB,FN1IJCB,FN1CJIB
      INTEGER ISY0DIS, ISY1DIS, IDEL, ISYHOP, IREAL
      INTEGER LU0IAJB, LU1IAJB
      INTEGER LU0IJCB, LU0CJIB, LU1IJCB, LU1CJIB
      INTEGER LUBFDA, LUBFDQA, LWORK
      INTEGER IADRBFA(*), IADRBFQA(*)
      INTEGER IT2DEL0(*), IT2DELB(*)
      INTEGER IT2DEL0A(*), IT2DELBA(*)
      INTEGER KDUM, IDUMMY
      PARAMETER (KDUM = +99 999 999) ! dummy address


#if defined (SYS_CRAY)
      REAL X0INT(*), X1INT(*), BF0RHF(*), BF1RHF(*)
      REAL XLAMD0P(*), XLAMD0H(*), XLAMDQP(*), XLAMDQH(*)
      REAL XLAMDAP(*), XLAMDAH(*), XLAMDQAP(*), XLAMDQAH(*)
      REAL DENS0(*), DPCK0(*), FOCK0(*)
      REAL DENSQ(*), DPCKQ(*), FOCKQ(*)
      REAL DENSA(*), DPCKA(*), FOCKA(*)
      REAL DENSQA(*), DPCKQA(*), FOCKQA(*)
      REAL RHOBFA(*),RHOBFQA(*)
      REAL WORK(LWORK)
      REAL ONE, ZERO, TWO, XNORM, DDOT, DNRM2
#else
      DOUBLE PRECISION X0INT(*), X1INT(*), BF0RHF(*), BF1RHF(*)
      DOUBLE PRECISION XLAMD0P(*), XLAMD0H(*), XLAMDQP(*), XLAMDQH(*)
      DOUBLE PRECISION XLAMDAP(*), XLAMDAH(*), XLAMDQAP(*), XLAMDQAH(*)
      DOUBLE PRECISION DENS0(*), DPCK0(*),FOCK0(*) 
      DOUBLE PRECISION DENSQ(*), DPCKQ(*), FOCKQ(*)
      DOUBLE PRECISION DENSA(*), DPCKA(*), FOCKA(*)
      DOUBLE PRECISION DENSQA(*), DPCKQA(*), FOCKQA(*)
      DOUBLE PRECISION RHOBFA(*),RHOBFQA(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ONE, ZERO, TWO, XNORM, DDOT, DNRM2
#endif
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0, TWO = 2.0d0)

      INTEGER ISYDEL, ISYMM1, ISYMM2, NMGD, KEND4, LWRK4, IADR, KMGD
      INTEGER KSCRCM2, KSCRCM1, KDSRHF
      INTEGER KX1IAJB, KX0IAJB, KXA1IJCB, KXA1CJIB, KXA0IJCB, KXA0CJIB 
      INTEGER LEN0, LEN1, LEN0A, LEN1A, ISYGAM, ISY0ALBE, ISY1ALBE
      INTEGER JGAM, KOFF0, KOFF1, ISYSRH1, KEND5, LWRK5, IOPT, ISYM
      INTEGER ISYMM0, ISYBF0, ISYBF1, ISYMTA
      INTEGER IADR0, IADRB, IADR0A, IADRBA
      INTEGER ISY0IAJ, ISY1IAJ, ISYA0IJC, ISYA1IJC,ISYHTA

*---------------------------------------------------------------------*
*     begin:
*---------------------------------------------------------------------*
      ISYDEL = ISAO(IDEL)
      D      = IDEL - IBAS(ISYDEL)
      
      ISYHTA = MULD2H(ISYHOP,ISYMTA)

*---------------------------------------------------------------------*
*     add 2-electr. contribution to AO Fock matrix Fbar (FOCKQ):
*     For CCS add the 2 electron part to FOCK0 too!!!!!!!!
*---------------------------------------------------------------------*

      IF (LRELAX) THEN
         CALL CC_AOFOCK2(X0INT,DENSQ,DPCKQ,FOCKQ,WORK,LWORK,
     &                   IDEL,ISY0DIS,ISYDEL,ISYHOP,.FALSE.)
         CALL CC_AOFOCK2(X0INT,DENSA,DPCKA,FOCKA,WORK,LWORK,
     &                   IDEL,ISY0DIS,ISYDEL,ISYMTA,.FALSE.)
         CALL CC_AOFOCK2(X0INT,DENSQA,DPCKQA,FOCKQA,WORK,LWORK,
     &                   IDEL,ISY0DIS,ISYDEL,ISYHTA,.FALSE.)
         IF (CCS.AND.LZERO) THEN
            CALL CC_AOFOCK2(X0INT,DENS0,DPCK0,FOCK0,WORK,LWORK,
     &                      IDEL,ISY0DIS,ISYDEL,ISYM0,.FALSE.)
         END IF
      END IF

      IF (LTWOEL) THEN
         CALL CC_AOFOCK2(X1INT,DENS0,DPCK0,FOCKQ,WORK,LWORK,
     &                   IDEL,ISY1DIS,ISYDEL,ISYM0,LX1ISQ)
         CALL CC_AOFOCK2(X1INT,DENSA,DPCKA,FOCKQA,WORK,LWORK,
     &                   IDEL,ISY1DIS,ISYDEL,ISYMTA,LX1ISQ)
      END IF

*---------------------------------------------------------------------*
*     for CCS we are done ...
*---------------------------------------------------------------------*
      IF (CCS) RETURN

*---------------------------------------------------------------------*
*     for CCSD calculate the first-order BF intermediates
*     the BF(A) intermediate only depends on TA
*     the BF(QA) intermediate depends on TA and IOPER 
*
*     for CC2 the F term and the G intermediate (skip)
*---------------------------------------------------------------------*
      ISYMM0 = MULD2H(ISYDEL,ISYM0)
      ISYMM1 = MULD2H(ISYDEL,ISYMTA)
      ISYMM2 = MULD2H(ISYDEL,ISYHTA)
      ISYBF0 = ISYMM0
      ISYBF1 = MULD2H(ISYDEL,ISYHOP)

*     -------------------------------------------
*     CCSD contributions: the BF intermediates...
*     -------------------------------------------
      IF (.NOT. CC2) THEN

*       --------------------------------------------------------
*       allocate an array for the different effective densities:
*       --------------------------------------------------------
        NMGD = 0
        DO ISYM = 1, NSYM
           NMGD = MAX(NMGD,NT2BGD(ISYM))            !max length
        END DO
       
        KMGD  = 1
        KEND4 = KMGD  + NMGD
        LWRK4 = LWORK - KEND4

        IF (LWRK4 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CCFBTAAO. (4)')
        END IF
*
* -------------------------------------------------------
* read in the BF(A) effective density and contract:
* with the PRESORTED g(1)(al-m,gam;del) --> result in RHOBFQA
*
* with the PRESORTED g(0)(al-m,gam;del) --> result in RHOBFA
* (only for a new T^A)
* LTWOEL/LRELAX not carried thru
* The BF intermediates are written on file OUTSIDE
* -------------------------------------------------------

*   read delta batch of the effective density DeltaA for BF(A) and BF(QA):

         IADR = IADRBFA(IDEL)
         NMGD = NT2BGD(ISYMM1)
         CALL GETWA2(LUBFDA,FNBFDA,WORK(KMGD),IADR,NMGD)

*   update BF(A) intermediate (RHOBFA_al i,kj, sym ISYMD*ISYMM):

         CALL CC_BFIB(RHOBFA,BF0RHF,ISYBF0,WORK(KMGD),ISYMM1,
     *                                  WORK(KEND4),LWRK4)

*   update BF(QA) intermediate:

         CALL CC_BFIB(RHOBFQA,BF1RHF,ISYBF1,WORK(KMGD),ISYMM1,
     *                                      WORK(KEND4),LWRK4)

*   read idelta batch of the effective density DeltaQA for BF(QA):

         IADR = IADRBFQA(IDEL)
         NMGD = NT2BGD(ISYMM2)
         CALL GETWA2(LUBFDQA,FNBFDQA,WORK(KMGD),IADR,NMGD)

*   update BF(QA) intermediate (add to previous contribution):
*   (added inside)

         CALL CC_BFIB(RHOBFQA,BF0RHF,ISYBF0,WORK(KMGD),ISYMM2,
     *                                      WORK(KEND4),LWRK4)

       END IF
*      ELSE
*       ---------------------------------------------------------
*       CC2 contributions: the F term and the G intermediate...
*       (the G term is here calculated in a certainly very clumsy
*       way using the one-index backtransformed amplitudes ...)
*       CC2 NOT YET IMPLEMENTED
*       ---------------------------------------------------------
C
C        KSCRCM1 = 1
C        KEND4   = KSCRCM1 + NT2BCD(ISYMM1)
C        LWRK4   = LWORK - KEND4
C
C        IF (LWRK4 .LT. 0) THEN
C           CALL QUIT('Insufficient work space in CCXIINTAO. (4)')
C        END IF
C
*       calculate one-index backtransformed amplitudes:
*       scrm1 - backtransformed with XLAMDH0 matrix
C        IOPT = 1
C        CALL CC_T2AO(T2AMP0,XLAMDH0,ISYM0,
C     &               WORK(KSCRCM1), WORK(KEND4),LWRK4,
C     &               IDEL, ISYDEL, ISYM0, IOPT )
C
C        IF (LTWOEL) THEN
C
*          ------------------------------------------
*          for CC2 the F term and the G intermediate:
*          ------------------------------------------
C           IOPT = 0
C           CALL CC_MOFCON2(X1INT,RHO2,XLAMDP0,XLAMDH0,
C     &                     XLAMDP0,XLAMDH0,XLAMDP0,XLAMDH0,
C     &                     ISYM0,ISYM0,ISYM0,ISYM0,
C     &                     WORK(KEND4),LWRK4,IDEL,
C     &                     ISYDEL,ISYHOP,ISYHOP,IOPT)
C
C           CALL CC_GIM(D1SRHF,ISY1DIS,WORK(KSCRCM1),ISYMM1,
C     &             XLAMDH0,ISYM0,GBIM,WORK(KEND4),LWRK4)
C
C        END IF
C
C        IF (LRELAX) THEN
C
*          -----------------------------------------------------
*          add the contributions from the relax. Lambda matrices
*          -----------------------------------------------------
C           KSCRCM2 = KEND4
C           KEND4   = KSCRCM2 + NT2BCD(ISYMM2)
C           LWRK4   = LWORK - KEND4
C
C           IF (LWRK4 .LT. 0) THEN
C             CALL QUIT('Insufficient work space in CCXIINTAO. (4b)')
C           END IF
C
*          calculate one-index backtransformed amplitudes:
*          scrm2 - backtransformed with XLAMDQH matrix
C           IOPT = 1
C           CALL CC_T2AO(T2AMP0,XLAMDQH,ISYHOP,
C     &                  WORK(KSCRCM2), WORK(KEND4),LWRK4,
C     &                  IDEL, ISYDEL, ISYM0, IOPT )
C
*          ------------------------------------------
*          for CC2 the F term and the G intermediate:
*          ------------------------------------------
C           IOPT = 0
C           CALL CC_MOFCON2(X0INT,RHO2,XLAMDQP,XLAMDQH,
C     &                     XLAMDP0,XLAMDH0,XLAMDP0,XLAMDH0,
C     &                     ISYHOP,ISYM0,ISYM0,ISYM0,
C     &                     WORK(KEND4),LWRK4,IDEL,
C     &                     ISYDEL,ISYHOP,ISYM0,IOPT)
C
C           IF (LZERO) THEN
*             ...without relaxation...
C              CALL CC_GIM(D0SRHF,ISY0DIS,WORK(KSCRCM1),ISYMM1,
C     &                    XLAMDH0,ISYM0,G0IM,WORK(KEND4),LWRK4)
C           END IF
C
*          ...relaxation of the XLAMDH used inside of CC_GIM...
C           CALL CC_GIM(D0SRHF,ISY0DIS,WORK(KSCRCM1),ISYMM1,
C     &                 XLAMDQH,ISYHOP,GBIM,WORK(KEND4),LWRK4)
C
*          ...relaxation of the XLAMDH used for T2 backtransf....
C           CALL CC_GIM(D0SRHF,ISY0DIS,WORK(KSCRCM2),ISYMM2,
C     &                 XLAMDH0,ISYM0,GBIM,WORK(KEND4),LWRK4)
C
C
C           ISYSRH1 = MULD2H(ISY0DIS,ISYHOP)
C           KDSRHF  = KEND4
C           KEND5   = KDSRHF + NDSRHF(ISYSRH1)
C           LWRK5   = LWORK - KEND5
C
C           IF (LWRK5 .LT. 0) THEN
C             CALL QUIT('Insufficient work space in CCXIINTAO. (5)')
C           END IF
C
*          ...relaxation of the XLAMDP used in CCTRBT....
C           CALL CCTRBT(X0INT,WORK(KDSRHF),XLAMDQP,
C     &                 ISYHOP,WORK(KEND5),LWRK5,ISY0DIS)
C
C           CALL CC_GIM(WORK(KDSRHF),ISYSRH1,WORK(KSCRCM1),ISYMM1,
C     &                 XLAMDH0,ISYM0,GBIM,WORK(KEND5),LWRK5)
C
C        END IF
C
C      END IF       
*---------------------------------------------------------------------*
*    calculate 3-index transformed integrals: 
*          (ia|j del), (ia|j del)-bar, 
*---------------------------------------------------------------------*

      ISY0IAJ = MULD2H(ISY0DIS,ISYM0)   !ISY0DIS * 3Lambda0
      ISY1IAJ = MULD2H(ISY0DIS,ISYHOP)  !ISY0DIS*2Lamda0*1LamdaQ=ISY1DIS*3Lamda0

C     -------------------------------------
C     allocate memory for integral batches:
C     -------------------------------------
      KX1IAJB = 1
      KEND4   = KX1IAJB + NT2BCD(ISY1IAJ)
*
      IF (LZERO) THEN
         KX0IAJB = KEND4
         KEND4   = KX0IAJB + NT2BCD(ISY0IAJ)
      END IF
*
      LWRK4   = LWORK - KEND4
*
      IF (LWRK4 .LT. 0) THEN
         CALL QUIT('Insufficient work space in CCFBTAAO1. (4b)')
      END IF
*
      IF (LZERO) THEN
         CALL DZERO(WORK(KX0IAJB),NT2BCD(ISY0IAJ))
      END IF
*
      CALL DZERO(WORK(KX1IAJB),NT2BCD(ISY1IAJ))
*
C     ---------------------------------------------------
C     do the 3-index transformation in a loop over gamma:
C     ---------------------------------------------------
      DO ISYGAM = 1, NSYM

        ISY0ALBE = MULD2H(ISY0DIS,ISYGAM)
        ISY1ALBE = MULD2H(ISY1DIS,ISYGAM)

        DO G = 1, NBAS(ISYGAM)
          JGAM = G + IBAS(ISYGAM)  !absolute index for gamma as IDEL
                   
          KOFF0 = IDSAOG(ISYGAM,ISY0DIS)+NNBST(ISY0ALBE)*(G-1)+1
          IF (LX1ISQ) THEN
            KOFF1 = IDSAOGSQ(ISYGAM,ISY1DIS)+N2BST(ISY1ALBE)*(G-1)+1
          ELSE
            KOFF1 = IDSAOG(ISYGAM,ISY1DIS)+NNBST(ISY1ALBE)*(G-1)+1
          END IF

          IOPT = 0
          CALL CC_IAJB(X0INT(KOFF0),ISY0ALBE,
     &                 X1INT(KOFF1),ISY1ALBE,
     &                 IDEL,JGAM,.FALSE.,IDUMMY,
     &                 WORK(KX0IAJB),WORK(KDUM),WORK(KDUM),
     &                 WORK(KX1IAJB),WORK(KDUM),WORK(KDUM),
     &                 XLAMD0P,XLAMD0H,ISYM0,XLAMDQP,XLAMDQH,ISYHOP,
     &                 XLAMD0P,XLAMD0H,ISYM0,XLAMDQP,XLAMDQH,ISYHOP,
     &                 WORK(KEND4),LWRK4,IOPT,LTWOEL,LRELAX,LZERO,
     &                 .TRUE.,LX1ISQ,IREAL)

        END DO

      END DO

*     --------------------------------------------
*     write 3-index transformed integrals to disk:
*     --------------------------------------------
      IF (LZERO) THEN
         LEN0 = NT2BCD(ISY0IAJ)

         CALL PUTWA2(LU0IAJB, FN0IAJB, WORK(KX0IAJB), IADR0, LEN0)

         IT2DEL0(IDEL) = IADR0
         IADR0 = IADR0 + LEN0
      END IF
*
      LEN1 = NT2BCD(ISY1IAJ)

      CALL PUTWA2(LU1IAJB, FN1IAJB, WORK(KX1IAJB), IADRB, LEN1)

      IT2DELB(IDEL) = IADRB
      IADRB = IADRB + LEN1

*---------------------------------------------------------------------*
*    calculate 3-index transformed integrals: 
*          (ij^|cb) + (ij|c^b), (ij^|cb) + (ij|c^b)-bar
*          (cj^|ib) + (c^j|ib), (cj^|ib) + (c^j|ib)-bar
*   for the C and D intermediates
*---------------------------------------------------------------------*
*
      ISYA0IJC = MULD2H(ISY0DIS,ISYMTA)  
      ISYA1IJC = MULD2H(ISY0DIS,MULD2H(ISYM0,ISYHTA))

C     -------------------------------------
C     allocate memory for integral batches:
C     -------------------------------------
      KXA1IJCB = 1
      KXA1CJIB = KXA1IJCB + NT2BCD(ISYA1IJC)
      KEND5    = KXA1CJIB + NT2BCD(ISYA1IJC)
*
      IF (LNEWTA) THEN              
         KXA0IJCB = KEND5
         KXA0CJIB = KXA0IJCB + NT2BCD(ISYA0IJC)
         KEND5    = KXA0CJIB + NT2BCD(ISYA0IJC)
      END IF
*
      LWRK5   = LWORK - KEND5
*
      IF (LWRK5 .LT. 0) THEN
         CALL QUIT('Insufficient work space in CCFBTAAO. (4b)')
      END IF

*  Initialize memory areas

      IF (LZERO.OR.LNEWTA) THEN
         CALL DZERO(WORK(KXA0IJCB),NT2BCD(ISYA0IJC))
         CALL DZERO(WORK(KXA0CJIB),NT2BCD(ISYA0IJC))
      END IF
      CALL DZERO(WORK(KXA1IJCB),NT2BCD(ISYA1IJC))
      CALL DZERO(WORK(KXA1CJIB),NT2BCD(ISYA1IJC))

C     ---------------------------------------------------
C     do the 3-index transformation in a loop over gamma:
C     ---------------------------------------------------
      DO ISYGAM = 1, NSYM

        ISY0ALBE = MULD2H(ISY0DIS,ISYGAM)
        ISY1ALBE = MULD2H(ISY1DIS,ISYGAM)

        DO G = 1, NBAS(ISYGAM)
          JGAM = G + IBAS(ISYGAM)
                   
          KOFF0 = IDSAOG(ISYGAM,ISY0DIS)+NNBST(ISY0ALBE)*(G-1)+1

          IF (LX1ISQ) THEN
            KOFF1 = IDSAOGSQ(ISYGAM,ISY1DIS)+
     &                                     N2BST(ISY1ALBE)*(G-1)+1
          ELSE
            KOFF1 = IDSAOG(ISYGAM,ISY1DIS)+NNBST(ISY1ALBE)*(G-1)+1
          END IF

          IOPT = 1                                                
          CALL CC_IJCB(X0INT(KOFF0),ISY0ALBE,X1INT(KOFF1),ISY1ALBE,
     &              IDEL,JGAM,
     &              WORK(KXA0IJCB),
     &              WORK(KXA0CJIB),
     &              WORK(KXA1IJCB),
     &              WORK(KXA1CJIB),
     &              XLAMD0P,XLAMD0H,ISYM0,XLAMDQP,XLAMDQH,ISYHOP,
     &              XLAMDAP,XLAMDAH,ISYMTA,XLAMDQAP,XLAMDQAH,ISYHTA,
     &              WORK(KEND5),LWRK5,
     &              IOPT,LTWOEL,LRELAX,LZERO,LNEWTA,LX1ISQ)

        END DO

      END DO

C     ------------------------------------
C     transform (cj|i del) to L(cj|i del):
C     ------------------------------------
      IF (LNEWTA) THEN
         CALL DSCAL(NT2BCD(ISYA0IJC),TWO,WORK(KXA0CJIB),1)
         CALL DAXPY(NT2BCD(ISYA0IJC),-ONE,WORK(KXA0IJCB),1,
     &                                   WORK(KXA0CJIB),1)
      END IF

      CALL DSCAL(NT2BCD(ISYA1IJC),TWO,WORK(KXA1CJIB),1)
      CALL DAXPY(NT2BCD(ISYA1IJC),-ONE,WORK(KXA1IJCB),1,
     &                                WORK(KXA1CJIB),1)


C     --------------------------------------------
C     write 3-index transformed integrals to disk:
C     --------------------------------------------
      IF (LOCDBG) THEN
        XNORM = DNRM2(NT2BCD(ISYA0IJC),WORK(KXA0IJCB),1)
        WRITE(LUPRI,*)'CCFBTAAO1> IDEL: ', idel
        WRITE(LUPRI,*)'Norm special integrals (0A):', XNORM
        XNORM = DNRM2(NT2BCD(ISYA1IJC),WORK(KXA1IJCB),1)
        WRITE(LUPRI,*)'Norm special integrals (BA):', XNORM
      END IF                                                  

      IF (LNEWTA) THEN
         LEN0A = NT2BCD(ISYA0IJC)

         CALL PUTWA2(LU0IJCB, FN0IJCB, WORK(KXA0IJCB), IADR0A, LEN0A)
         CALL PUTWA2(LU0CJIB, FN0CJIB, WORK(KXA0CJIB), IADR0A, LEN0A)

         IT2DEL0A(IDEL) = IADR0A
         IADR0A = IADR0A + LEN0A
* else?
      END IF


      LEN1A = NT2BCD(ISYA1IJC)

      CALL PUTWA2(LU1IJCB, FN1IJCB, WORK(KXA1IJCB), IADRBA, LEN1A)
      CALL PUTWA2(LU1CJIB, FN1CJIB, WORK(KXA1CJIB), IADRBA, LEN1A)

      IT2DELBA(IDEL) = IADRBA
      IADRBA = IADRBA + LEN1A

*---------------------------------------------------------------------*
* that's it;  return:
*---------------------------------------------------------------------*
      RETURN
      END 
*=====================================================================*
*                  END OF SUBROUTINE CCFBTAAO1                        *
*=====================================================================*
