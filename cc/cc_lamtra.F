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
C=================================================================*
C  /* Deck cc_lamtr */
      SUBROUTINE CC_LAMTRA(XLAMDP,ISYMLP,XLAMDPR,ISYMLPR,
     *                     XLAMDH,ISYMLH,XLAMDHR,ISYMLHR,
     *                     R1AM,ISYMR1)
C=================================================================*
C
C     PURPOSE:
C             transform general symmetry lambda matrices 
C             with a general symmetry R(c j) vector (R1AM)
C             occupied in XLAMDP transformed to virtual
C             virtual in XLAMDH transformed to occupied
C
C     Sonia Coriani 25-11-1998
C     Based on Ove's CCLR_LAMTRA
C     Debugged 9.8.99
C
C     XLAMDP XLAMDH   = Lambda^p(alp,k) and Lambda^h(alp,c)
C     XLAMDPR XLAMDHR = Lambda^{R,p}(alp,c) and Lambda^{R,h}(alp,k) 
C                       transformed
C==================================================================*
#include <implicit.h>
#include <iratdef.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, XMONE= -1.0D00)
      DIMENSION R1AM(*),XLAMDP(*),XLAMDH(*),XLAMDPR(*),XLAMDHR(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
C
      IF (IPRINT .GT.25) THEN
         CALL AROUND('IN CC_LAMTR  ')
      ENDIF
      IF (MULD2H(ISYMLH,ISYMR1).NE.ISYMLHR) THEN
         CALL QUIT('Symmetry mismatch for ISYMLHR in CC_LAMTRA.')
      END IF
      IF (MULD2H(ISYMLP,ISYMR1).NE.ISYMLPR) THEN
         CALL QUIT('Symmetry mismatch for ISYMLPR in CC_LAMTRA.')
      END IF
C
C-----------------------------------------
C     Transform general lambda particle matrix.
C     LaP~(al,a) = -sum(k)[ LaH(al,k)*R(a,k)]
C NB!! note the minus sign.
C Note that transformed LambdaP is always (alpha,VIRTUAL) 
C     the (alpha,OCCUPIED) block is = ZERO
C-----------------------------------------
C
      CALL DZERO(XLAMDPR,NGLMDT(ISYMLPR))  
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMK  = MULD2H(ISYMR1,ISYMA)
         ISYALF = MULD2H(ISYMK,ISYMLP)
C
         NBASALF = MAX(NBAS(ISYALF),1)
         NBASA   = MAX(NVIR(ISYMA),1)
C
         KOFF1  = IGLMRH(ISYALF,ISYMK) + 1 
         KOFF2  = IT1AM(ISYMA,ISYMK)   + 1
         KOFF3  = IGLMVI(ISYALF,ISYMA) + 1
C
         CALL DGEMM('N','T',NBAS(ISYALF),NVIR(ISYMA),NRHF(ISYMK),
     *              XMONE,XLAMDP(KOFF1),NBASALF,R1AM(KOFF2),NBASA,
     *              ZERO,XLAMDPR(KOFF3),NBASALF)
C
  100 CONTINUE
C
C-----------------------------------------
C     Transform Lambda hole matrix.
C     LaH~(al,i) = + sum(c)[ LaH(al,c)*C(c,i)]
C     Note that transformed LambdaH is always (alpha,OCC.) 
C     the (alpha,VIRTUAL) block is = ZERO
C-----------------------------------------
C
      CALL DZERO(XLAMDHR,NGLMDT(ISYMLHR))
C
      DO 200 ISYMI = 1,NSYM
C
         ISYMC  = MULD2H(ISYMR1,ISYMI)
         ISYALF = MULD2H(ISYMC,ISYMLH)
C
         NBASALF = MAX(NBAS(ISYALF),1)      
         NBASC   = MAX(NVIR(ISYMC),1)
C
         KOFF1  = IGLMVI(ISYALF,ISYMC) + 1 
         KOFF2  = IT1AM(ISYMC,ISYMI) + 1  
         KOFF3  = IGLMRH(ISYALF,ISYMI) + 1
C
         CALL DGEMM('N','N',NBAS(ISYALF),NRHF(ISYMI),NVIR(ISYMC),
     *              ONE,XLAMDH(KOFF1),NBASALF,R1AM(KOFF2),NBASC,
     *              ZERO,XLAMDHR(KOFF3),NBASALF)
C
  200 CONTINUE
C
      RETURN
      END
*=================================================================*
C  /* Deck cclt_z1a */
      SUBROUTINE CCLT_Z1A(CTR1,ISYCTR,TA1,ISYMTA1,ISYMZA,ZAKJ)
*-----------------------------------------------------------------*
C
C     Purpose: To calculate the Zeta1^A intermediate:
C
C     ZetaA(k j)  = - sum_c CTR1(c j) t^A(c k) 
C     
C     or in general 
C
C     Result(k j) = - sum_c Left(c j) Right(c k)
C
C     Sonia Coriani, November 1998
C     Debug 16.08.1999 OK
*-----------------------------------------------------------------*
#include <implicit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION CTR1(*),TA1(*),ZAKJ(*)


      ISYMCJ = ISYCTR           ! symmetry of Left vector
      ISYMCK = ISYMTA1          ! symmetry of Right vector
      ISYRES = MULD2H(ISYCTR,ISYMTA1)
c      IF (ISYRES.NE.ISYMZA) 

      DO ISYMK = 1, NSYM
         ISYMC = MULD2H(ISYMK,ISYMCK)
         ISYMJ = MULD2H(ISYMC,ISYMCJ)

         KOFF1 = 1 + IT1AM(ISYMC,ISYMK)    ! for right vector
         KOFF2 = 1 + IT1AM(ISYMC,ISYMJ)    ! for left vector
         KOFF3 = 1 + IMATIJ(ISYMK,ISYMJ)   ! for result (occ,occ)

         NRHFK = MAX(NRHF(ISYMK),1)        ! total # occupied K
         NVIRC = MAX(NVIR(ISYMC),1)        ! total # virtual C

         CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMJ),NVIR(ISYMC),
     *                -ONE,TA1(KOFF1),NVIRC,CTR1(KOFF2),NVIRC,
     *                ZERO,ZAKJ(KOFF3),NRHFK)

      END DO
      RETURN
      END
C-------------------------------------------------------------------C
