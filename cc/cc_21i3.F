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
C  /* Deck cc_21i3 */
*--------------------------------------------------------------
      SUBROUTINE CC_21I3(RHO1, XINT0, ISYDIS0, XINT1, ISYDIS1,
     *                   PINT0,QINT0,ISYPQ0,PINT1,QINT1,ISYPQ1,
     *                   PAINT0,QAINT0,ISYPQA0,
     *                   PAINT1,QAINT1,ISYPQA1,
     *                   XLAMP0,XLAMH0, ISYML0,XLAMPQ,ISYMLQ,
     *                   XLAMPA,ISYMLA,XLAMPQA,ISYMLQA,
     *                   WORK,LWORK,IOPT,LRHOAO,LDERIV,LXI1SQ)
*--------------------------------------------------------------
*
*     Generalization of the CC_21I and CC_21I2 routines 
*     for the F contribution (21I contr) to the F^bT^A result vector
*     (Calculation of the FQA intermediate)
*
*     Sonia Coriani February 1999
*
*     Version: 2.0 
*
*     XINT0, ISYDIS0  -- (* *|* delta) batch of integrals, its symmetry
*     XINT1, ISYDIS1 -- the same for derivative integrals
*     PINT0, QINT0   -- P_ci,j (del), Q_ci,j (del)      (ISYPQ0)
*     PINT1, QINT1   -- Pbar_ci,j (del), Qbar_ci,j(del) (ISYPQ1)
*     PAINT0,QAINT0  -- PA_ci,j (del) QA_ci,j(del)      (ISYPQA0)
*     PAINT1,QAINT1  -- PAbar_ci,j QAbar_ci,j (del)     (ISYPQA1)
*
*     LRHOAO = .TRUE.    return result with a index in AO
*     LRHOAO = .FALSE.   return result with a index in MO
*     LDERIV = .FALSE.   omit contribution from derivative integrals
*     LXI1SQ = .TRUE.    (a b|* *) of XINT1 is a full matrix (fx LAO)
*
*     IOPT = 1   The FA intermediate (test case)
*     IOPT = 2   The FQA intermediate
*------------------------------------------------------------------
*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#   include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      LOGICAL LRHOAO, LDERIV, LXI1SQ
      INTEGER LWORK,ISYDIS0,ISYDIS1,IOPT
      INTEGER ISYPQ0,ISYPQ1,ISYPQA0,ISYPQA1, ip
      INTEGER ISYML0, ISYMLQ, ISYMLA, ISYMLQA, ISYM
C
#if defined (SYS_CRAY)
      REAL RHO1(*), PINT0(*), PINT1(*), QINT0(*), QINT1(*) 
      REAL PAINT0(*), PAINT1(*), QAINT0(*), QAINT1(*)
      REAL XINT0(*), XINT1(*), WORK(LWORK)
      REAL XLAMP0(*), XLAMPA(*), XLAMPQ(*), XLAMPQA(*)
      REAL XLAMH0(*)
      REAL DUMMY, ZERO, ONE, TWO, DDOT, XNORM
#else
      DOUBLE PRECISION RHO1(*), PINT0(*), PINT1(*), QINT0(*), QINT1(*)
      DOUBLE PRECISION PAINT0(*), PAINT1(*), QAINT0(*), QAINT1(*)
      DOUBLE PRECISION XINT0(*), XINT1(*), WORK(LWORK)
      DOUBLE PRECISION XLAMP0(*), XLAMPA(*), XLAMPQ(*), XLAMPQA(*)
      DOUBLE PRECISION XLAMH0(*)
      DOUBLE PRECISION DUMMY, ZERO, ONE, TWO, DDOT, XNORM
#endif
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      INTEGER NT3BAO(8), IT3BAO(8,8)
      INTEGER ISYGIJ1, ISYGIJ2, ISYDENBA, ISYDENA
      INTEGER ISYRES, ISYMJ, ISYMGI, ISYMAL, ISY3AO, ISYMG
      INTEGER ISYALI, ISYMI, ISALBE, ISYMBE, ISYMA, ISYGAM
      INTEGER KRHOAO, KPQMO, KPQONE, KPQTWO, KDENBA, KDENA
      INTEGER KAOINT,NTOTA,NPQMO, NPQMO01, NPQMOA1, NTOTGA
      INTEGER KEND1,LWRK1,KOFF1,KOFF2,KOFF3,KOFF4,KOFF5
      INTEGER ICOUNT, NTOTAL, NTOTGI, NGI, NAGI, NAIG, NTOTBE, NTOTAJ
      INTEGER IOPTDEN, IOPTCON
C
C     set some symmetries used globally in this routine:
C
      ISYGIJ1   = MULD2H(ISYPQ0,ISYMLA)
*      IF (IOPT.EQ.2) THEN
*      ISYGIJ2   = MULD2H(ISYPQ0,ISYMLQA)
*      END IF
      ISYGIJ2   = MULD2H(ISYPQA0,ISYMLQ)
      ISYDENBA  = MULD2H(ISYGIJ1,ISYMLQ)
      ISYDENA   = MULD2H(ISYGIJ1,ISYML0)
      ISYRES    = MULD2H(ISYDENBA,ISYDIS0)
C
*      IF ( IOPT.GE.2 .AND. MULD2H(ISYVWZ2,ISYXL0).NE.ISYGIJ ) THEN
*            CALL QUIT('Symmetry mismatch in CC_21I3NEW.')
*      END IF
C
C     calculate index arrays for intermediates with 3 indices in AO:
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
C          --  intermediate with two more indices backtransformed
C----------------------------------------------------------------
C
      NPQMO = 1
      DO ISYM = 1, NSYM
        NPQMO = MAX(NPQMO,NT2BCD(ISYM))
      END DO
        
C
C      IF (IOPT.GE.2) NPQMO = MAX(NPQMO,NT2BCD(ISYVWZ2))
C
      KRHOAO = 1
      KPQMO  = KRHOAO + NT1AO(ISYRES)
      KPQONE = KPQMO  + NPQMO
      IF (IOPT.EQ.2) THEN
        KPQTWO = KPQONE + NT2BGD(ISYGIJ1)
      ELSE 
        KPQTWO = KPQONE
*        ISYGIJ2 = ISYGIJ1 
      END IF

      KDENBA = KPQTWO + NT2BGD(ISYGIJ2)
      KEND1  = KDENBA + NT3BAO(ISYDENBA)

      LWRK1  = LWORK  - KEND1
 
      IF ((IOPT.EQ.2).AND.LDERIV) THEN
         KDENA = KEND1
         KEND1 = KDENA + NT3BAO(ISYDENA)
      END IF
      LWRK1   = LWORK - KEND1
C
      IF ( LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient work space in CC_21I3.')
      ENDIF
C
C--------------------------------------------------------------------
C     1. Exchange like contribution: add P and Q intermediates
C     together and transform virtual index "c" back to AO
C     with the LambdaA^p matrix
C     add PA and QA intermediates together and backtransform  
C     "c" with the Lambda^p matrix. Add to previous
C     Repeat for Pbar (P1) and Qbar (Q1)
C     Repeat for PAbar (PA1) and QAbar (QA1)
C--------------------------------------------------------------------
C
      CALL DCOPY(NT2BCD(ISYPQ0),PINT0,1,WORK(KPQMO),1)
      CALL DAXPY(NT2BCD(ISYPQ0),ONE,QINT0,1,WORK(KPQMO),1)
C
      !LambdaA*PQ -> KPQONE
      CALL CC_BFAO(WORK(KPQONE),WORK(KPQMO),ISYPQ0,XLAMPA,ISYMLA,
     &             1,1,DUMMY,1,DUMMY,1,ZERO,0)
c
c      XNORM = DDOT(NT2BGD(ISYGIJ1),WORK(KPQONE),1,WORK(KPQONE),1)
c      WRITE (LUPRI,*) 'CC_21I3> Norm(PQ-AO)1=',XNORM
c
      IF (IOPT.EQ.2) THEN
        !LambdaQA*PQ -> KPQTWO
        CALL CC_BFAO(WORK(KPQTWO),WORK(KPQMO),ISYPQ0,XLAMPQA,ISYMLQA,
     &                    1,1,DUMMY,1,DUMMY,1,ZERO,0)
      END IF

      CALL DCOPY(NT2BCD(ISYPQA0),PAINT0,1,WORK(KPQMO),1)                 
c      WRITE(LUPRI,*) 'CC_21I3> test PQ'
c      do ip = 1, 5
c       WRITE(LUPRI,*) PAINT0(ip)
c      end do
      CALL DAXPY(NT2BCD(ISYPQA0),ONE,QAINT0,1,WORK(KPQMO),1)

      !Lambda0*calPQ -> + KPQONE
      CALL CC_BFAO(WORK(KPQONE),WORK(KPQMO),ISYPQA0,XLAMP0,ISYML0,
     &             1,1,DUMMY,1,DUMMY,1,ONE,0)
c
c      XNORM = DDOT(NT2BGD(ISYGIJ1),WORK(KPQONE),1,WORK(KPQONE),1)
c      WRITE (LUPRI,*) 'CC_21I3> Norm(PQA-AO)1=',XNORM
c

      IF (IOPT.EQ.2) THEN
        !LambdaQ*calPQ -> + KPQTWO
        CALL CC_BFAO(WORK(KPQTWO),WORK(KPQMO),ISYPQA0,XLAMPQ,ISYMLQ,
     &               1,1,DUMMY,1,DUMMY,1,ONE,0)
C--------------------------------------------------------------------
C     Backtransform also barred PQ  and calPQ and add to KPQTWO
C--------------------------------------------------------------------

        CALL DCOPY(NT2BCD(ISYPQ1),PINT1,1,WORK(KPQMO),1)           
        CALL DAXPY(NT2BCD(ISYPQ1),ONE,QINT1,1,WORK(KPQMO),1)
C
        !LambdaA*barPQ -> + KPQTWO
        CALL CC_BFAO(WORK(KPQTWO),WORK(KPQMO),ISYPQ1,XLAMPA,ISYMLA,
     &                1,1,DUMMY,1,DUMMY,1,ONE,0)

        CALL DCOPY(NT2BCD(ISYPQA1),PAINT1,1,WORK(KPQMO),1)           
        !LambdaA*bar-calPQ -> + KPQTWO
        CALL DAXPY(NT2BCD(ISYPQA1),ONE,QAINT1,1,WORK(KPQMO),1)

        CALL CC_BFAO(WORK(KPQTWO),WORK(KPQMO),ISYPQA1,XLAMP0,ISYML0,
     &                1,1,DUMMY,1,DUMMY,1,ONE,0)
      END IF
C
C--------------------------------------------------------------------
C     transform occupied index j back to AO (alpha)
C--------------------------------------------------------------------
C
      IOPTDEN = 0   !initialize density with actual contribution
      !backt j with XLAMPQ
      CALL CC_21IDEN(WORK(KPQONE),ISYGIJ1,XLAMPQ,ISYMLQ,WORK(KDENBA),
     *               ISYDENBA,IT3BAO,WORK(KEND1),LWRK1,'EXC',IOPTDEN)

      IF (IOPT.EQ.2) THEN
        IOPTDEN = 1   !add to density
        !backt j with XLAMP0
        CALL CC_21IDEN(WORK(KPQTWO),ISYGIJ2,XLAMP0,ISYML0,
     *                  WORK(KDENBA),ISYDENBA,IT3BAO,
     *                  WORK(KEND1),LWRK1,'EXC',IOPTDEN)

        IF (LDERIV) THEN
          IOPTDEN = 0
          CALL CC_21IDEN(WORK(KPQONE),ISYGIJ1,XLAMP0,ISYML0,
     *                     WORK(KDENA),ISYDENA,IT3BAO,
     *                     WORK(KEND1),LWRK1,'EXC',IOPTDEN)
        END IF
      END IF
C
c     XNORM = DDOT(NT3BAO(ISYGIJ1),WORK(KDENBA),1,WORK(KPQDEN),1)
c     WRITE (LUPRI,*) 'Norm(PQDEN)1=',XNORM
C
C--------------------------------------------------------------------
C     2. Coulomb like contribution: scale P intermediate with -2
C     and transform virtual "c" back to AO
C--------------------------------------------------------------------
C
      CALL DCOPY(NT2BCD(ISYPQ0),PINT0,1,WORK(KPQMO),1)
      CALL DSCAL(NT2BCD(ISYPQ0),-TWO,WORK(KPQMO),1)

      CALL CC_BFAO(WORK(KPQONE),WORK(KPQMO),ISYPQ0,XLAMPA,ISYMLA,
     &             1,1,DUMMY,1,DUMMY,1,ZERO,0)
C
c     XNORM = DDOT(NT2BGD(ISYGIJ),WORK(KPQAO),1,WORK(KPQAO),1)
c     WRITE (LUPRI,*) 'Norm(PQAO)3=',XNORM
C
C
      IF (IOPT.EQ.2) THEN
        CALL CC_BFAO(WORK(KPQTWO),WORK(KPQMO),ISYPQ0,XLAMPQA,
     &               ISYMLQA,1,1,DUMMY,1,DUMMY,1,ZERO,0)
      END IF
C

      CALL DCOPY(NT2BCD(ISYPQA0),PAINT0,1,WORK(KPQMO),1)
      CALL DSCAL(NT2BCD(ISYPQA0),-TWO,WORK(KPQMO),1)

      CALL CC_BFAO(WORK(KPQONE),WORK(KPQMO),ISYPQA0,XLAMP0,ISYML0,
     &                                  1,1,DUMMY,1,DUMMY,1,ONE,0)

      IF (IOPT.EQ.2) THEN
         CALL CC_BFAO(WORK(KPQTWO),WORK(KPQMO),ISYPQA0,XLAMPQ,ISYMLQ,
     &                                     1,1,DUMMY,1,DUMMY,1,ONE,0)

C
c     XNORM = DDOT(NT2BGD(ISYGIJ),WORK(KPQAO),1,WORK(KPQAO),1)
c     WRITE (LUPRI,*) 'Norm(PQAO)4=',XNORM
C

         CALL DCOPY(NT2BCD(ISYPQ1),PINT1,1,WORK(KPQMO),1)
         CALL DSCAL(NT2BCD(ISYPQ1),-TWO,WORK(KPQMO),1)

         CALL CC_BFAO(WORK(KPQTWO),WORK(KPQMO),ISYPQ1,XLAMPA,ISYMLA,
     &                                    1,1,DUMMY,1,DUMMY,1,ONE,0)

         CALL DCOPY(NT2BCD(ISYPQA1),PAINT1,1,WORK(KPQMO),1)
         CALL DSCAL(NT2BCD(ISYPQA1),-TWO,WORK(KPQMO),1)

         CALL CC_BFAO(WORK(KPQTWO),WORK(KPQMO),ISYPQA1,XLAMP0,ISYML0,
     &                                     1,1,DUMMY,1,DUMMY,1,ONE,0)

      END IF
C--------------------------------------------------------------------
C     transform occupied index j back to AO (gamma), add result to
C     the effective density stored in KPQDEN
C--------------------------------------------------------------------
C

      IOPTDEN = 1
      CALL CC_21IDEN(WORK(KPQONE),ISYGIJ1,XLAMPQ,ISYMLQ,WORK(KDENBA),
     *             ISYDENBA,IT3BAO,WORK(KEND1),LWRK1,'COU',IOPTDEN)

      IF (IOPT.EQ.2) THEN
      IOPTDEN = 1
      CALL CC_21IDEN(WORK(KPQTWO),ISYGIJ2,XLAMP0,ISYML0,
     *                  WORK(KDENBA),ISYDENBA,IT3BAO,
     *                  WORK(KEND1),LWRK1,'COU',IOPTDEN)

      IF (LDERIV) THEN
         IOPTDEN = 1
         CALL CC_21IDEN(WORK(KPQONE),ISYGIJ1,XLAMP0,ISYML0,
     *                     WORK(KDENA),ISYDENA,IT3BAO,
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

      CALL CC_21ICON(WORK(KRHOAO),ISYRES,WORK(KDENBA),ISYDENBA,
     *               XINT0,ISYDIS0,IT3BAO,WORK(KEND1),LWRK1,1)

      IF ((IOPT.EQ.2).AND.LDERIV) THEN

         IF (LXI1SQ) THEN
            IOPTCON = 2             !X2INT has full (a b| already
         ELSE
            IOPTCON = 1             !X2INT is squared inside CC_21ICON1
         END IF

         CALL CC_21ICON(WORK(KRHOAO),ISYRES,WORK(KDENA),ISYDENA,
     *            XINT1,ISYDIS1,IT3BAO,WORK(KEND1),LWRK1,IOPTCON)
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
        CALL CC_T1AM(RHO1,ISYRES,WORK(KRHOAO),ISYRES,XLAMH0,ISYML0,ONE)
      END IF
C
      RETURN
      END 
C
*=====================================================================*
