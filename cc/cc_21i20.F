C  /* Deck cc_21i2 */
      SUBROUTINE CC_21I20(RHO1, XINT, ISYDIS, XINT2, ISYDIS2,
     *                   PINT0,QINT0,ISYVWZ0,PINT2,QINT2,ISYVWZ2,
     *                   XLAMP0,XLAMH0,ISYXL0,XLAMP2,ISYXL2,
     *                   WORK,LWORK,IOPT,LRHOAO,LDERIV)
C
C     Written by Asger Halkier 31/10 - 1995
C     restructured by Christof Haettig July 1998
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
C
C     Purpose: To calculate the 21I contribution to RHO1!
C
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#   include <priunit_dec.h>
#   include <ccorb_dec.h>
#   include <ccsdsym_dec.h>
#else
#   include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      LOGICAL LRHOAO, LDERIV
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
      INTEGER KPQAO0, KPQDEN0, ISYGIJ0, ISYDEN0, IOPTDEN
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
C     together and transform virtual index e back to AO
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

      CALL CC_21ICON0(WORK(KRHOAO),ISYRES,WORK(KPQDEN),ISYDEN,
     *               XINT,ISYDIS,IT3BAO,WORK(KEND1),LWRK1)

      IF (IOPT.EQ.3 .AND. LDERIV) THEN

         CALL CC_21ICON0(WORK(KRHOAO),ISYRES,WORK(KPQDEN0),ISYDEN0,
     *              XINT2,ISYDIS2,IT3BAO,WORK(KEND1),LWRK1)
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
*=====================================================================*
C  /* Deck cc_21icon */
      SUBROUTINE CC_21ICON0(RHOAO,ISYRES,PQDEN,ISYDEN,XINT,ISYDIS,
     *                     IT3BAO,WORK,LWORK)
C
C     Purpose: contract the effective density build in CC_21I
C              with the 2-el integrals to RHO1 contribution
C
C     Written by Christof Haettig 16 July 1998
C
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#   include <priunit_dec.h>
#   include <ccorb_dec.h>
#   include <ccsdsym_dec.h>
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER ISYRES, ISYDEN, ISYDIS, IT3BAO(8,8), LWORK
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
      DO ISYMG = 1,NSYM
C
         ISALBE = MULD2H(ISYDIS,ISYMG)
C
         IF ( LWORK .LT. N2BST(ISALBE) ) THEN
            CALL QUIT('Insufficient work space in CC_21I.')
         ENDIF
C
         DO G = 1, NBAS(ISYMG)
C
C           Square up integral distribution (al be| ga de).
C
            KOFF1 = IDSAOG(ISYMG,ISYDIS) + NNBST(ISALBE)*(G - 1) + 1
C
            CALL CCSD_SYMSQ(XINT(KOFF1),ISALBE,WORK)
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
                 WRITE (LUPRI,*) 'Symmetry mismatch in CC_21ICON0:'
                 WRITE (LUPRI,*) ISYALI,ISYMG,ISYDEN
                 CALL QUIT('Symmetry mismatch in CC_21ICON0.')
              END IF

              KOFF3 = IAODIS(ISYMAL,ISYMBE) + 1
              KOFF4 = IT3BAO(ISYALI,ISYMG) + NT1AO(ISYALI)*(G-1) 
     *                                     + IT1AO(ISYMAL,ISYMI) + 1
              KOFF5 = IT1AO(ISYMBE,ISYMI) + 1
 
              NTOTAL = MAX(NBAS(ISYMAL),1)
              NTOTBE = MAX(NBAS(ISYMBE),1)
 
              CALL DGEMM('T','N',NBAS(ISYMBE),NRHF(ISYMI),NBAS(ISYMAL),
     *                   -ONE,WORK(KOFF3),NTOTAL,PQDEN(KOFF4),NTOTAL,
     *                    ONE,RHOAO(KOFF5),NTOTBE)

            END DO

         END DO 
 
      END DO

      RETURN
      END 
*=====================================================================*
