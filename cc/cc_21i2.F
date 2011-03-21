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
C  /* Deck cc_21i2 */
      SUBROUTINE CC_21I2(RHO1, XINT, ISYDIS, XINT2, ISYDIS2,
     *                   PINT0,QINT0,ISYVWZ0,PINT2,QINT2,ISYVWZ2,
     *                   XLAMP0,XLAMH0,ISYXL0,XLAMP2,ISYXL2,
     *                   WORK,LWORK,IOPT,LRHOAO,LDERIV,LX2ISQ)
C
C     Written by Asger Halkier 31/10 - 1995
C     restructured by Christof Haettig July 1998
C     X2INT already squared in input, Sonia Coriani November 1999
C
C     Version: 2.0 
C
C     XINT,  ISYDIS  -- (* *|* delta) batch of integrals, its symmetry
C     XINT2, ISYDIS2 -- the same for derivative integrals (for IOPT=3)
C     PINT0, QINT0   -- P, Q linear combination of ZWV intermediates
C     PINT2, QINT2   -- P, Q calculated from response amplitudes 
C
C     IOPT =  1 : calculate only usual LHTR contributions
C             2 : include response contrib. from PINT2, QINT2 (F mat.)
C             3 : include additional orbital relaxation 
C                 and derivative contributions (for ETA/RHS vectors)
C
C     LRHOAO = .TRUE.    return result with a index in AO
C     LRHOAO = .FALSE.   return result with a index in MO
C     LDERIV = .FALSE.   omit contribution from derivative integrals
C     LX2ISQ = .TRUE.    (a b|* *) of XINT2 is a full matrix (fx LAO)
C
C     Purpose: To calculate the 21I contribution to RHO1!
C
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#   include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      LOGICAL LRHOAO, LDERIV, LX2ISQ
      INTEGER LWORK,ISYDIS,ISYDIS2,ISYVWZ0,ISYVWZ2,ISYXL0,ISYXL2,IOPT
C
#if defined (SYS_CRAY)
      REAL RHO1(*), PINT0(*), PINT2(*), QINT0(*), QINT2(*), 
     *     XINT(*), XINT2(*), WORK(LWORK),
     *     XLAMP0(*), XLAMH0(*),XLAMP2(*)
      REAL DUMMY, ZERO, ONE, TWO
#else
      DOUBLE PRECISION RHO1(*), PINT0(*), PINT2(*), QINT0(*), QINT2(*),
     *                 XINT(*), XINT2(*), WORK(LWORK),
     *                 XLAMP0(*), XLAMH0(*),XLAMP2(*)
      DOUBLE PRECISION DUMMY, ZERO, ONE, TWO
#endif
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      INTEGER NT3BAO(8), IT3BAO(8,8)
      INTEGER ISYGIJ, ISYRES, ISY3AO, ISYMG
      INTEGER ISYALI, ISYDEN
      INTEGER KRHOAO, KPQAO, KPQMO, KPQDEN, NPQMO
      INTEGER KEND1,LWRK1
      INTEGER ICOUNT
      INTEGER KPQAO0, KPQDEN0, ISYGIJ0, ISYDEN0, IOPTDEN, IOPTCON
C
C     set some symmetries use globally in this routine:
C
      ISYGIJ  = MULD2H(ISYVWZ0,ISYXL2)
      ISYDEN  = MULD2H(ISYGIJ,ISYXL0)
      ISYRES  = MULD2H(ISYDEN,ISYDIS)
      ISYGIJ0 = MULD2H(ISYVWZ0,ISYXL0)
      ISYDEN0 = MULD2H(ISYGIJ0,ISYXL0)
C
      IF ( IOPT.GE.2 .AND. MULD2H(ISYVWZ2,ISYXL0).NE.ISYGIJ ) THEN
            CALL QUIT('Symmetry mismatch in CC_21I2.')
      END IF
C
C     calculate index arrays for intermediates with 3 indeces in AO:
C
      DO ISY3AO = 1, NSYM
        ICOUNT = 0
        DO ISYMG = 1, NSYM
           ISYALI = MULD2H(ISY3AO,ISYMG)
           IT3BAO(ISYALI,ISYMG) = ICOUNT
           ICOUNT = ICOUNT + NBAS(ISYMG) * NT1AO(ISYALI)
        END DO
        NT3BAO(ISY3AO) = ICOUNT
      END DO
C
C----------------------------------------------------------------
C     allocate work space for 
C          --  intermediate with one index backtransformed
C          --  intermediate with one more index backtransformed
C          --  intermediate with two more indeces backtransformed
C----------------------------------------------------------------
C
      NPQMO  = NT2BCD(ISYVWZ0)
      IF (IOPT.GE.2) NPQMO = MAX(NPQMO,NT2BCD(ISYVWZ2))
C
      KRHOAO = 1
      KPQMO  = KRHOAO + NT1AO(ISYRES)
      KPQAO  = KPQMO  + NPQMO
      KPQDEN = KPQAO  + NT2BGD(ISYGIJ)
      KEND1  = KPQDEN + NT3BAO(ISYDEN)
      LWRK1  = LWORK  - KEND1
C
      IF (IOPT.EQ.3) THEN
        KPQAO0  = KEND1
        KEND1   = KPQAO0  + NT2BGD(ISYGIJ0)
        IF (LDERIV) THEN
           KPQDEN0 = KEND1
           KEND1   = KPQDEN0 + NT3BAO(ISYDEN0)
        END IF
        LWRK1   = LWORK   - KEND1
      END IF
C
      IF ( LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient work space in CC_21I.')
      ENDIF
C
C--------------------------------------------------------------------
C     1. Exchange like contribution: add P and Q intermediates
C     together and transform virtual index back to AO
C--------------------------------------------------------------------
C
      CALL DCOPY(NT2BCD(ISYVWZ0),PINT0,1,WORK(KPQMO),1)
      CALL DAXPY(NT2BCD(ISYVWZ0),ONE,QINT0,1,WORK(KPQMO),1)
C
      CALL CC_BFAO(WORK(KPQAO),WORK(KPQMO),ISYVWZ0,XLAMP2,ISYXL2,
     &             1,1,DUMMY,1,DUMMY,1,ZERO,0)
C
c     XNORM = DDOT(NT2BGD(ISYGIJ),WORK(KPQAO),1,WORK(KPQAO),1)
c     WRITE (LUPRI,*) 'Norm(PQAO)1=',XNORM
C
C--------------------------------------------------------------------
C     if IOPT = 3 calculate PQAO0, backtransformed with XLAMP0
C     for if IOPT = 2/3 include response contribution:
C--------------------------------------------------------------------
C
      IF (IOPT.EQ.3) THEN
        CALL CC_BFAO(WORK(KPQAO0),WORK(KPQMO),ISYVWZ0,XLAMP0,ISYXL0,
     &               1,1,DUMMY,1,DUMMY,1,ZERO,0)
      END IF
C
      IF (IOPT.EQ.2 .OR. IOPT.EQ.3) THEN
 
         CALL DCOPY(NT2BCD(ISYVWZ2),PINT2,1,WORK(KPQMO),1)
         CALL DAXPY(NT2BCD(ISYVWZ2),ONE,QINT2,1,WORK(KPQMO),1)
C
         CALL CC_BFAO(WORK(KPQAO),WORK(KPQMO),ISYVWZ2,XLAMP0,ISYXL0,
     &                1,1,DUMMY,1,DUMMY,1,ONE,0)

      END IF
C
c     XNORM = DDOT(NT2BGD(ISYGIJ),WORK(KPQAO),1,WORK(KPQAO),1)
c     WRITE (LUPRI,*) 'Norm(PQAO)2=',XNORM
C
C--------------------------------------------------------------------
C     transform occupied index j back to AO (alpha)
C          for IOPT=3 add extra relaxation and transform also
C          the j index of the KPQAO0 intermediate
C--------------------------------------------------------------------
C
      IOPTDEN = 0
      CALL CC_21IDEN(WORK(KPQAO),ISYGIJ,XLAMP0,ISYXL0,WORK(KPQDEN),
     *               ISYDEN,IT3BAO,WORK(KEND1),LWRK1,'EXC',IOPTDEN)

      IF (IOPT.EQ.3) THEN
         IOPTDEN = 1
         CALL CC_21IDEN(WORK(KPQAO0),ISYGIJ0,XLAMP2,ISYXL2,
     *                  WORK(KPQDEN),ISYDEN,IT3BAO,
     *                  WORK(KEND1),LWRK1,'EXC',IOPTDEN)

         IF (LDERIV) THEN
            IOPTDEN = 0
            CALL CC_21IDEN(WORK(KPQAO0),ISYGIJ0,XLAMP0,ISYXL0,
     *                     WORK(KPQDEN0),ISYDEN,IT3BAO,
     *                     WORK(KEND1),LWRK1,'EXC',IOPTDEN)
         END IF
      END IF
C
c     XNORM = DDOT(NT3BAO(ISYGIJ),WORK(KPQDEN),1,WORK(KPQDEN),1)
c     WRITE (LUPRI,*) 'Norm(PQDEN)1=',XNORM
C
C--------------------------------------------------------------------
C     2. Coulomb like contribution: scale P intermediate with -2
C     and transform virutal e back to AO
C--------------------------------------------------------------------
C
      CALL DCOPY(NT2BCD(ISYVWZ0),PINT0,1,WORK(KPQMO),1)
      CALL DSCAL(NT2BCD(ISYVWZ0),-TWO,WORK(KPQMO),1)

      CALL CC_BFAO(WORK(KPQAO),WORK(KPQMO),ISYVWZ0,XLAMP2,ISYXL2,
     &             1,1,DUMMY,1,DUMMY,1,ZERO,0)
C
c     XNORM = DDOT(NT2BGD(ISYGIJ),WORK(KPQAO),1,WORK(KPQAO),1)
c     WRITE (LUPRI,*) 'Norm(PQAO)3=',XNORM
C
C--------------------------------------------------------------------
C     if IOPT = 3 calculate PQAO0, backtransformed with XLAMP0
C     for IOPT = 2/3 include response contribution:
C--------------------------------------------------------------------
C
      IF (IOPT.EQ.3) THEN
        CALL CC_BFAO(WORK(KPQAO0),WORK(KPQMO),ISYVWZ0,XLAMP0,ISYXL0,
     &               1,1,DUMMY,1,DUMMY,1,ZERO,0)
      END IF
C
      IF (IOPT.EQ.2 .OR. IOPT.EQ.3) THEN

        CALL DCOPY(NT2BCD(ISYVWZ2),PINT2,1,WORK(KPQMO),1)
        CALL DSCAL(NT2BCD(ISYVWZ2),-TWO,WORK(KPQMO),1)

        CALL CC_BFAO(WORK(KPQAO),WORK(KPQMO),ISYVWZ2,XLAMP0,ISYXL0,
     &               1,1,DUMMY,1,DUMMY,1,ONE,0)

      END IF
C
c     XNORM = DDOT(NT2BGD(ISYGIJ),WORK(KPQAO),1,WORK(KPQAO),1)
c     WRITE (LUPRI,*) 'Norm(PQAO)4=',XNORM
C
C--------------------------------------------------------------------
C     transform occupied index j back to AO (gamma), add result to
C     the effective density stored in KPQDEN
C--------------------------------------------------------------------
C
      IOPTDEN = 1
      CALL CC_21IDEN(WORK(KPQAO),ISYGIJ,XLAMP0,ISYXL0,WORK(KPQDEN),
     *             ISYDEN,IT3BAO,WORK(KEND1),LWRK1,'COU',IOPTDEN)

      IF (IOPT.EQ.3) THEN
         IOPTDEN = 1
         CALL CC_21IDEN(WORK(KPQAO0),ISYGIJ0,XLAMP2,ISYXL2,
     *                  WORK(KPQDEN),ISYDEN,IT3BAO,
     *                  WORK(KEND1),LWRK1,'COU',IOPTDEN)

         IF (LDERIV) THEN
            IOPTDEN = 1
            CALL CC_21IDEN(WORK(KPQAO0),ISYGIJ0,XLAMP0,ISYXL0,
     *                     WORK(KPQDEN0),ISYDEN,IT3BAO,
     *                     WORK(KEND1),LWRK1,'COU',IOPTDEN)
         END IF
      END IF
C
C--------------------------------------------------------------------
C     contract the effective density with the 2-el integrals:
C--------------------------------------------------------------------
C
c     XNORM = DDOT(NT3BAO(ISYGIJ),WORK(KPQDEN),1,WORK(KPQDEN),1)
c     WRITE (LUPRI,*) 'Norm(PQDEN)2=',XNORM
C
      CALL DZERO(WORK(KRHOAO),NT1AO(ISYRES))

      CALL CC_21ICON(WORK(KRHOAO),ISYRES,WORK(KPQDEN),ISYDEN,
     &                XINT,ISYDIS,IT3BAO,WORK(KEND1),LWRK1,1)

      IF (IOPT.EQ.3 .AND. LDERIV) THEN

         IF (LX2ISQ) THEN
            IOPTCON = 2             !X2INT has full (a b| already
         ELSE 
            IOPTCON = 1             !X2INT is squared inside CC_21ICON
         END IF
        
         CALL CC_21ICON(WORK(KRHOAO),ISYRES,WORK(KPQDEN0),ISYDEN0,
     &              XINT2,ISYDIS2,IT3BAO,WORK(KEND1),LWRK1,IOPTCON)
      END IF
C
c     XNORM = DDOT(NT1AO(ISYRES),RHO1,1,RHO1,1)
c     WRITE (LUPRI,*) 'XNORM=',XNORM
c
C---------------------------------------------
C     transformation and/or storage in result:
C---------------------------------------------
C
      IF (LRHOAO) THEN
        CALL DAXPY(NT1AO(ISYRES),ONE,WORK(KRHOAO),1,RHO1,1)
      ELSE
        CALL CC_T1AM(RHO1,ISYRES,WORK(KRHOAO),ISYRES,XLAMH0,ISYXL0,ONE)
      END IF
C
      RETURN
      END 
C
*=====================================================================*
C  /* Deck cc_t1am */
      SUBROUTINE CC_T1AM(RHO1MO,ISYMO,RHO1AO,ISYAO,XLAMD,ISYLAM,FAC)
C
C  Purpose: transform AO index of a two-index array with 
C           IT1AO structure to MO and store in IT1AM structure
C
C  RHO1MO(ai) = FAC*RHO1MO(ai) + sum_alp XLAMD(alp a) RHO1AO(alp i)
C
C  Written by Christof Haettig 16 July 1998
C
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      DIMENSION RHO1MO(*), RHO1AO(*), XLAMD(*)
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
C
      IF (MULD2H(ISYAO,ISYLAM).NE.ISYMO ) THEN
         CALL QUIT('Symmetry mismatch in CC_T1AM.')
      END IF
C
      DO ISYMA = 1, NSYM
C
         ISYMI  = MULD2H(ISYMO,ISYMA)
         ISYMBE = MULD2H(ISYLAM,ISYMA)
C
         KOFF1 = IGLMVI(ISYMBE,ISYMA) + 1
         KOFF2 = IT1AO(ISYMBE,ISYMI)  + 1
         KOFF3 = IT1AM(ISYMA,ISYMI)   + 1
C
         NTOTA  = MAX(NVIR(ISYMA),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMBE),
     *              ONE,XLAMD(KOFF1),NTOTBE,RHO1AO(KOFF2),NTOTBE,
     *              FAC,RHO1MO(KOFF3),NTOTA)
C
      END DO
C
      RETURN
      END
*=====================================================================*
C  /* Deck cc_21iden */
      SUBROUTINE CC_21IDEN(PQAO,ISYGIJ,XLAMDP,ISYMXL,PQDEN,ISYDEN,
     *                     IT3BAO,WORK,LWORK,FLAG,IOPT)
C
C     Purpose: transform and resort the two-index backtransformed
C              P+Q intermediate used in CC_21I to an density like 
C              intermediate with three indeces in AO and one in MO
C
C     FLAG = 'EXC' do transformation/resort for exchange like contrib.
C            'COU' do transformation for coulomb like contrib.
C
C     IOPT = 0 intialize PQDEN with the actual contribution
C            1 add the actual contribution to what is already in PQDEN
C
C     Written by Christof Haettig July 1998
C
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#   include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CHARACTER*(*) FLAG
      INTEGER LWORK, ISYGIJ, ISYDEN, ISYMXL, IOPT, IT3BAO(8,8)
C
#if defined (SYS_CRAY)
      REAL PQAO(*), PQDEN(*), XLAMDP(*), WORK(LWORK)
      REAL ZERO, ONE, FACT
#else
      DOUBLE PRECISION PQAO(*), PQDEN(*), XLAMDP(*), WORK(LWORK)
      DOUBLE PRECISION ZERO, ONE, FACT
#endif
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
C
      INTEGER ISYMJ,ISYMGI,ISYMAL,ISYMI,ISYALI,ISYGAM,KOFF1,KOFF2,KOFF3
      INTEGER NTOTGA, NTOTAL, NTOTGI, NGI, NAGI, NAIG, KOFF4, ISYMG
C
C--------------------------------------------------------------------
C     transform occupied index j back to AO (alpha)
C--------------------------------------------------------------------
C
      IF (FLAG(1:3).EQ.'COU') THEN

         FACT = ZERO
         IF (IOPT.EQ.1) FACT = ONE

         DO ISYMJ = 1, NSYM

            ISYMGI = MULD2H(ISYGIJ,ISYMJ)
            ISYGAM = MULD2H(ISYMXL,ISYMJ)

            KOFF1  = IGLMRH(ISYGAM,ISYMJ)  + 1
            KOFF2  = IT2BGD(ISYMGI,ISYMJ)  + 1
            KOFF3  = IT3BAO(ISYMGI,ISYGAM) + 1

            NTOTGA = MAX(NBAS(ISYGAM),1)
            NTOTGI = MAX(NT1AO(ISYMGI),1)

            CALL DGEMM('N','T',NT1AO(ISYMGI),NBAS(ISYGAM),NRHF(ISYMJ),
     *                 ONE,PQAO(KOFF2),NTOTGI,XLAMDP(KOFF1),NTOTGA,
     *                 FACT,PQDEN(KOFF3),NTOTGI)

         END DO 
C
C--------------------------------------------------------------------
C     transform occupied index j back to AO (alpha)
C--------------------------------------------------------------------
C
      ELSE IF (FLAG(1:3).EQ.'EXC') THEN

         DO ISYMJ = 1, NSYM

            ISYMGI = MULD2H(ISYGIJ,ISYMJ)
            ISYMAL = MULD2H(ISYMXL,ISYMJ)

            KOFF1  = IGLMRH(ISYMAL,ISYMJ) + 1
            KOFF2  = IT2BGD(ISYMGI,ISYMJ) + 1

            IF ( LWORK .LT. NBAS(ISYMAL)*NT1AO(ISYMGI) ) THEN
               CALL QUIT('Insufficient work space in CC_21IDEN.')
            ENDIF

            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTGI = MAX(NT1AO(ISYMGI),1)

            CALL DGEMM('N','T',NBAS(ISYMAL),NT1AO(ISYMGI),NRHF(ISYMJ),
     *                 ONE,XLAMDP(KOFF1),NTOTAL,PQAO(KOFF2),NTOTGI,
     *                 ZERO,WORK,NTOTAL)
C
C           ------------------------------------------------------------
C           resort from PQDEN(alpha, gamma i) to PQDEN(alpha i, gamma)
C           ------------------------------------------------------------
C
            DO ISYMG = 1, NSYM

               ISYMI  = MULD2H(ISYMGI,ISYMG)
               ISYALI = MULD2H(ISYMAL,ISYMI)

               DO G = 1, NBAS(ISYMG)

                  KOFF4 = IT3BAO(ISYALI,ISYMG) 
     *                  + NT1AO(ISYALI)*(G-1) + IT1AO(ISYMAL,ISYMI)

                  IF (IOPT.EQ.1) THEN
                     DO I = 1, NRHF(ISYMI)
                        NGI  = IT1AO(ISYMG,ISYMI) + NBAS(ISYMG)*(I-1)+G
                        NAGI = NBAS(ISYMAL)*(NGI-1) + 1
                        NAIG = KOFF4 + NBAS(ISYMAL)*(I-1) + 1

                        CALL DAXPY(NBAS(ISYMAL),ONE,WORK(NAGI),1,
     *                                              PQDEN(NAIG),1)
                     END DO
                  ELSE
                     DO I = 1, NRHF(ISYMI)
                        NGI  = IT1AO(ISYMG,ISYMI) + NBAS(ISYMG)*(I-1)+G
                        NAGI = NBAS(ISYMAL)*(NGI-1) + 1
                        NAIG = KOFF4 + NBAS(ISYMAL)*(I-1) + 1

                        CALL DCOPY(NBAS(ISYMAL),WORK(NAGI),1,
     *                                          PQDEN(NAIG),1)
                     END DO
                  END IF

               END DO

            END DO

         END DO

      ELSE
         CALL QUIT('Illegal FLAG in CC_21IDEN.')
      END IF

      RETURN

      END
*=====================================================================*
C  /* Deck cc_21icon */
      SUBROUTINE CC_21ICON(RHOAO,ISYRES,PQDEN,ISYDEN,XINT,ISYDIS,
     *                     IT3BAO,WORK,LWORK,IOPT)
*---------------------------------------------------------------------*
C
C     Purpose: contract the effective density build in CC_21I
C              with the 2-el integrals to RHO1 contribution
C
C     Written by Christof Haettig 16 July 1998
C     Modified by Sonia Coriani 07 November 1999 to handle 
C     full (a b| g d) in input (IOPT = 2)
C     IOPT = 1, usual packed a>= b
C
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER ISYRES, ISYDEN, ISYDIS, IT3BAO(8,8), LWORK
      INTEGER IOPT 
C
#if defined (SYS_CRAY)
      REAL RHOAO(*), PQDEN(*), XINT(*), WORK(LWORK)
      REAL ZERO, ONE
#else
      DOUBLE PRECISION RHOAO(*), PQDEN(*), XINT(*), WORK(LWORK)
      DOUBLE PRECISION ZERO, ONE
#endif
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
C
      INTEGER ISYMG, ISALBE, ISYMAL, ISYMBE, ISYMI, ISYALI
      INTEGER KOFF1, KOFF3, KOFF4, KOFF5, NTOTAL, NTOTBE
C
C Start: check option
C      
      IF ((IOPT.NE.1) .AND. (IOPT.NE.2))
     &      CALL QUIT('Unknown IOPT in CC_21ICON (start)')
C
      DO ISYMG = 1,NSYM
C
         ISALBE = MULD2H(ISYDIS,ISYMG)
C
         IF ((IOPT.EQ.1).AND.(LWORK .LT. N2BST(ISALBE))) THEN
            CALL QUIT('Insufficient work space in CC_21ICON.')
         ENDIF
C 
         DO G = 1, NBAS(ISYMG)

C           ----------------------------------------------------------
C           If IOPT = 1 Square up integral distribution (al be| ga de).
C           ----------------------------------------------------------

            IF (IOPT.EQ.1) THEN

              KOFF1 = IDSAOG(ISYMG,ISYDIS) + 
     &                                 NNBST(ISALBE)*(G - 1) + 1
              CALL CCSD_SYMSQ(XINT(KOFF1),ISALBE,WORK)

            END IF
C
C           ------------------------------------------------
C           Final contraction and storage in result.
C           ------------------------------------------------
C
            DO ISYMAL = 1, NSYM

              ISYMBE = MULD2H(ISALBE,ISYMAL)
              ISYMI  = MULD2H(ISYRES,ISYMBE)
              ISYALI = MULD2H(ISYMAL,ISYMI)

              IF (MULD2H(ISYALI,ISYMG).NE.ISYDEN) THEN
                 WRITE (LUPRI,*) 'Symmetry mismatch in CC_21ICON:'
                 WRITE (LUPRI,*) ISYALI,ISYMG,ISYDEN
                 CALL QUIT('Symmetry mismatch in CC_21ICON.')
              END IF

              KOFF4 = IT3BAO(ISYALI,ISYMG) + NT1AO(ISYALI)*(G-1) 
     &                                     + IT1AO(ISYMAL,ISYMI) + 1
              KOFF5 = IT1AO(ISYMBE,ISYMI) + 1
 
              NTOTAL = MAX(NBAS(ISYMAL),1)
              NTOTBE = MAX(NBAS(ISYMBE),1)
 
              IF (IOPT.EQ.1) THEN

                 KOFF3 = IAODIS(ISYMAL,ISYMBE) + 1
 
                 CALL DGEMM('T','N',NBAS(ISYMBE),NRHF(ISYMI),
     &                      NBAS(ISYMAL),
     &                     -ONE,WORK(KOFF3),NTOTAL,PQDEN(KOFF4),NTOTAL,
     &                      ONE,RHOAO(KOFF5),NTOTBE)

              ELSE IF (IOPT.EQ.2) THEN

                 KOFF3 = IDSAOGSQ(ISYMG,ISYDIS)+ N2BST(ISALBE)*(G - 1)
     &                   + IAODIS(ISYMAL,ISYMBE) + 1

                 CALL DGEMM('T','N',NBAS(ISYMBE),NRHF(ISYMI),
     &                      NBAS(ISYMAL),
     &                     -ONE,XINT(KOFF3),NTOTAL,PQDEN(KOFF4),NTOTAL,
     &                      ONE,RHOAO(KOFF5),NTOTBE)
              ELSE
                
                 CALL QUIT('Unknown option in CC_21ICON')

              END IF

            END DO

         END DO 
 
      END DO

      RETURN
      END 
*=====================================================================*
