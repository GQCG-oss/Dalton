c /* deck cc_iabj */
*=====================================================================*
       SUBROUTINE CC_IAJBNEW( X1INT,   ISY1ALBE,
     &                     X2INT,   ISY2ALBE, 
     &                     IDEL,    IGAM, 
     &                     XIAJB,   XIABJ,   XIJBA, 
     &                     X2IAJB,  X2IABJ,  X2IJBA,
     &                     XLAMDP1, XLAMDH1, ISYM1,
     &                     XLAMDP2, XLAMDH2, ISYM2,
     &                     WORK,    LWORK,   
     &                     IOPT,    LDERIV,  LRELAX, LZERO, 
     &                     LX2ISQ,  IREAL )
*---------------------------------------------------------------------*
*
*   Purpose: generalized transformation to (ia|jb) & (ia|bj) integrals
*            for the two-index (**|gam del) approach
*            assumes three-index arrays XIAJB & XIABJ in core
*
*            this routine drives the transformation of the indices 
*            ia and j, the transformation of the delta index to b has 
*            to be done from the outside.  
*
*            X1INT, XIAJB,  XIABJ,  XIJBA  : usual integrals
*            X2INT, X2IAJB, X2IABJ, X2IJBA : derivative integrals
*            
*            IOPT=0: (ia|j del) only
*
*            IOPT=1: (ia|j del) and (ia|del j) integrals
*
*            IOPT=2: (ia|j del), (ia|del j), and (ij|del a) integrals
*
*            IF LDERIV=.TRUE. transform also the derivative integrals:
*                IOPT=0: (ia|j del)-bar
*                IOPT=1: (ia|j del)-bar, (ia|del j)-bar
*                IOPT=2: (ia|j del)-bar, (ia|del j)-bar, (ij|del a)-bar 
*                If LX2ISQ=.TRUE. The derivative integrals are already 
*                                 squared in input
*                if IREAL = -1 The derivative integrals are LAO 
*                               i.e. antisymmetric (and already squared)
*                if IREAL = +1 derivative integrals assumed symmetric
*
*            IF LRELAX=.TRUE. include relaxation contribution to the
*            derivative integrals using XLAMDP2, XLAMDH2 matrices:
*
*            IF LZERO=.FALSE. skip calculation of zero-order integrals
*                              (ia|j del), (ia|del j), and (ij|del a)
*
*                IOPT=0: (ia|j del)-bar += (i-bar a|j del) 
*                                    + (i a-bar|j del) + (ia|j-bar del)
*
*                IOPT=1: as for IOPT=0, but add also
*                        (ia|del j)-bar += (i-bar a|del j) 
*                                    + (i a-bar|del j) + (ia|del j-bar)
*
*                IOPT=2: as for IOPT=1, but add also
*                        (ij|del a)-bar= (i-bar j|del j) 
*                                   + (i j|del a-bar) + (i j-bar|del a)
*
*            (ia|jb) integrals:
*                  i,j          transform. with XLAMDP1 with sym. ISYM1
*                  a            transform. with XLAMDH1 with sym. ISYM1
*                  i-bar, j-bar transform. with XLAMDP2 with sym. ISYM2
*                  a-bar        transform. with XLAMDH2 with sym. ISYM2
*             
*            (ia|bj) integrals:
*                  i            transform. with XLAMDP1 with sym. ISYM1
*                  a,j          transform. with XLAMDH1 with sym. ISYM1
*                  i-bar        transform. with XLAMDP2 with sym. ISYM2
*                  a-bar, j-bar transform. with XLAMDH2 with sym. ISYM2
*             
*            (ij|ba) integrals:
*                  i            transform. with XLAMDP1 with sym. ISYM1
*                  a,j          transform. with XLAMDH1 with sym. ISYM1
*                  i-bar        transform. with XLAMDP2 with sym. ISYM2
*                  a-bar, j-bar transform. with XLAMDH2 with sym. ISYM2
*             
*            for IOPT < 10 are XLAMDP2, XLAMDH2, ISYM2,
*                              X2IABJ, X2IJBA dummy input
*            for IOPT=0, 10 or 20 XIABJ, XIJBA are dummy input
*            for IOPT=1, 11 or 21 XIJBA is dummy input
*            for IOPT < 20 X2INT, ISY2ALBE are dummy input
*
*
*    Written by Christof Haettig, May 1998.
*    Restructured by Sonia Coriani, October 1999 in order to handle
*    squared X2INT in input and LAO's
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#  include <priunit_dec.h>
#  include <ccorb_dec.h>
#  include <ccsdsym_dec.h>
#  include <mxorb_dec.h>
#  include <ccisao_dec.h>
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <mxorb.h>
#include <ccisao.h>

#if defined (SYS_CRAY)
      REAL ONE, ZERO
#else
      DOUBLE PRECISION ONE, ZERO
#endif
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0)

      LOGICAL LDERIV, LRELAX, LZERO, LX2ISQ, LLAO
      INTEGER IDEL,IGAM,ISY1ALBE,ISY2ALBE,ISYM1,ISYM2,LWORK,IOPT,IREAL
      
#if defined (SYS_CRAY)
      REAL XLAMDP1(*), XLAMDH1(*)
      REAL XLAMDP2(*), XLAMDH2(*)
      REAL X1INT(*), X2INT(*)
      REAL XIAJB(*),  XIABJ(*),  XIJBA(*)
      REAL X2IAJB(*), X2IABJ(*), X2IJBA(*)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION XLAMDP1(*), XLAMDH1(*)
      DOUBLE PRECISION XLAMDP2(*), XLAMDH2(*)
      DOUBLE PRECISION X1INT(*), X2INT(*)
      DOUBLE PRECISION XIAJB(*),  XIABJ(*),  XIJBA(*)
      DOUBLE PRECISION X2IAJB(*), X2IABJ(*), X2IJBA(*)
      DOUBLE PRECISION WORK(LWORK)
#endif

      INTEGER ISYM11, ISYM12, ISYMAI, ISYMJ, ISYMB, ISYMBJ, ISYMIJ
      INTEGER ISYMI, ISYMA, ISYALP, ISYBET, ISYDEL, ISYGAM, KOFF0
      INTEGER KSCR1, KEND1, LWRK1, KLAMJ, KLAMB
      INTEGER KINT, KLAMD, KIAJB, KOFF1, KOFF2, KOFF3, KOFF4, KSCR6
      INTEGER NNBSAB, NBASG, NBASA, NBASB, NVIRA, NAI, KOFF6
      INTEGER NRHFI, KSCR3, KOFF5

*---------------------------------------------------------------------*
*     set some symmetries
*---------------------------------------------------------------------*
      ISYM11  = MULD2H(ISYM1,ISYM1)
      ISYM12  = MULD2H(ISYM1,ISYM2)

      IF      ( (IOPT.EQ. 0) .OR. (IOPT.EQ. 1) .OR. (IOPT.EQ. 2) ) THEN
        ISYMAI  = MULD2H(ISY1ALBE,ISYM11)
      ELSE IF ( (IOPT.EQ.10) .OR. (IOPT.EQ.11) .OR. (IOPT.EQ.12) ) THEN
        ISYMAI  = MULD2H(ISY1ALBE,ISYM12)
      ELSE
        CALL QUIT('Unknown IOPT in CC_IAJB.')
      END IF

      ISYDEL  = ISAO(IDEL)
      ISYGAM  = ISAO(IGAM)

      ISYMJ   = MULD2H(ISYM1,ISYGAM)
      ISYMB   = MULD2H(ISYM1,ISYDEL)
      ISYMBJ  = MULD2H(ISYMB,ISYMJ)

      IF ( (IOPT.EQ.0) .AND. (ISYMAI.GT.ISYMBJ) ) RETURN


      IF      (IREAL.EQ.+1 .OR. (.NOT.LDERIV)) THEN
        LLAO = .FALSE.
      ELSE IF (IREAL.EQ.-1 .AND. LDERIV) THEN
        LLAO = .TRUE.
      ELSE
        CALL QUIT('ILLEGAL VALUE OF IREAL IN CC_IAJBNEW.')
      END IF                                               

*---------------------------------------------------------------------*
*     work space allocation:
*
*     KSCR1 --  I^{del,gam}(alp bet)
*     KSCR3 --  I^{del,gam}(alp bet)[1] (only if LDERIV and .NOT.LX2ISQ)
*
*---------------------------------------------------------------------*
      KSCR1   = 1
      KEND1   = KSCR1 + N2BST(ISY1ALBE)

      IF (LDERIV.AND.(.NOT.LX2ISQ)) THEN
        KSCR3 = KEND1
        KEND1 = KSCR3 + N2BST(ISY2ALBE)
      ELSE
C
C VB mar 2001:
C     KDUM should not be used as an address (it was also undefined)
C        KSCR3 = KDUM
        KSCR3 = KEND1
      END IF

      LWRK1   = LWORK - KEND1

      IF ( LWRK1 .LT. 0) THEN
        CALL QUIT('Insufficient memory in CC_IAJB.')
      END IF

*---------------------------------------------------------------------*
*     square integral matrix/matrices up 
*---------------------------------------------------------------------*

      CALL CCSD_SYMSQ(X1INT,ISY1ALBE,WORK(KSCR1))
      
      IF ((LDERIV).AND.(.NOT.LX2ISQ)) THEN
        CALL CCSD_SYMSQ(X2INT,ISY2ALBE,WORK(KSCR3))
      END IF
*---------------------------------------------------------------------*
*     call routine for actual calculation of transformed integrals
*---------------------------------------------------------------------*
      IF (.NOT.LX2ISQ) THEN
         CALL CCIAJB(WORK(KSCR1),ISY1ALBE,WORK(KSCR3),ISY2ALBE,
     &               IDEL, IGAM, XIAJB,   XIABJ,   XIJBA, 
     &               X2IAJB,  X2IABJ,  X2IJBA, XLAMDP1, XLAMDH1, ISYM1,
     &               XLAMDP2, XLAMDH2, ISYM2, WORK(KEND1),  LWRK1,
     &               IOPT,    LDERIV,  LRELAX, LZERO, LLAO )
      ELSE
         CALL CCIAJB(WORK(KSCR1),ISY1ALBE,X2INT,ISY2ALBE,
     &               IDEL, IGAM, XIAJB,   XIABJ,   XIJBA, 
     &               X2IAJB,  X2IABJ,  X2IJBA, XLAMDP1, XLAMDH1, ISYM1,
     &               XLAMDP2, XLAMDH2, ISYM2, WORK(KEND1),  LWRK1,
     &               IOPT,    LDERIV,  LRELAX, LZERO, LLAO )
      END IF

*---------------------------------------------------------------------*
*     return
*---------------------------------------------------------------------*

      RETURN
      END
*=====================================================================*
*                 END OF SUBROUTINE CC_IAJB                          *
*=====================================================================*
c /* deck cciabj */
*=====================================================================*
       SUBROUTINE CCIAJB( XINT1,   ISY1ALBE,
     &                    XINT2,   ISY2ALBE, 
     &                    IDEL,    IGAM, 
     &                    XIAJB,   XIABJ,   XIJBA, 
     &                    X2IAJB,  X2IABJ,  X2IJBA,
     &                    XLAMDP1, XLAMDH1, ISYM1,
     &                    XLAMDP2, XLAMDH2, ISYM2,
     &                    WORK,    LWORK,   
     &                    IOPT,    LDERIV,  LRELAX, LZERO, LLAO )
*---------------------------------------------------------------------*
*
*   Purpose: generalized transformation to (ia|jb) & (ia|bj) integrals
*            for the two-index (**|gam del) approach
*            assumes three-index arrays XIAJB & XIABJ in core
*
*
*   Note:    The (ia|bj) and (ij|ba) integrals are calculated by 
*            transforming the gamma index with XLAMDH matrices, 
*            assuming symmetric AO integrals 
*            For antisymmetric integrals (LAO), beta and alpha
*            have to be switched, and sign changed
*
*            See CC_IAJB for documentation
*            LLAO = .TRUE. Derivative integrals are London AO (antisym)
*            If is assumed that LLAO=TRUE implies LDERIV=TRUE
*
*    Written by Christof Haettig, May 1998.
*    Restructured by Sonia Coriani to handle LAO's, October 1999
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#  include <priunit_dec.h>
#  include <ccorb_dec.h>
#  include <ccsdsym_dec.h>
#  include <mxorb_dec.h>
#  include <ccisao_dec.h>
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <mxorb.h>
#include <ccisao.h>

#if defined (SYS_CRAY)
      REAL ONE, ZERO
#else
      DOUBLE PRECISION ONE, ZERO
#endif
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0)

      LOGICAL LDERIV, LRELAX, LZERO, LLAO
      INTEGER IDEL, IGAM, ISY1ALBE, ISY2ALBE, ISYM1, ISYM2, LWORK, IOPT
      
#if defined (SYS_CRAY)
      REAL XLAMDP1(*), XLAMDH1(*)
      REAL XLAMDP2(*), XLAMDH2(*)
      REAL XINT1(*), XINT2(*)
      REAL XIAJB(*),  XIABJ(*),  XIJBA(*)
      REAL X2IAJB(*), X2IABJ(*), X2IJBA(*)
      REAL WORK(LWORK)
      REAL FAC
#else
      DOUBLE PRECISION XLAMDP1(*), XLAMDH1(*)
      DOUBLE PRECISION XLAMDP2(*), XLAMDH2(*)
      DOUBLE PRECISION XINT1(*), XINT2(*)
      DOUBLE PRECISION XIAJB(*),  XIABJ(*),  XIJBA(*)
      DOUBLE PRECISION X2IAJB(*), X2IABJ(*), X2IJBA(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION FAC
#endif

      INTEGER ISYM11, ISYM12, ISYMAI, ISYMJ, ISYMB, ISYMBJ, ISYMIJ
      INTEGER ISYMI, ISYMA, ISYALP, ISYBET, ISYDEL, ISYGAM, KOFF0
      INTEGER KSCR2, KSCR4, KSCR5, KEND1, LWRK1, KLAMJ, KLAMB
      INTEGER KINT, KLAMD, KIAJB, KOFF1, KOFF2, KOFF3, KOFF4, KSCR6
      INTEGER NNBSAB, NBASG, NBASA, NBASB, NVIRA, NAI, KOFF6
      INTEGER NRHFI, KSCR7, KOFF7, KOFF5, KSCR8, KOFF8, KSCR68

*---------------------------------------------------------------------*
* set some symmetries
*---------------------------------------------------------------------*
      ISYM11  = MULD2H(ISYM1,ISYM1)
      ISYM12  = MULD2H(ISYM1,ISYM2)

      ISYDEL  = ISAO(IDEL)
      ISYGAM  = ISAO(IGAM)

*---------------------------------------------------------------------*
*     work space allocation:
*
* KSCR2 - I^{del,gam}(i bet) 
*         (If LRELAX reuse for      I^{del,gam}(i-bar bet)) 
*         (If LDERIV add/reuse for  I^[1]{del,gam}(i bet))
* KSCR4 - I^{del,gam}(i a)
* KSCR5 - I^{del,gam}(i j) (if IOPT=2)
* KSCR6 - I^{del,gam}(i a-bar) + I^{del,gam}(i-bar a) + I^[1]{del,gam}(i a)
* KSCR7 - I^{del,gam}(i j-bar) + I^{del,gam}(i-bar j) + I^[1]{del,gam}(i j)
*
* If IOPT >= 1 .and. LLAO:
* KSCR8 - The IA derivative integrals  contribution to X2IAJB
*---------------------------------------------------------------------*
      KSCR2   = 1
      KSCR4   = KSCR2 + NBAST*NRHFT
      KEND1   = KSCR4 + NVIRT*NRHFT

      IF (IOPT.EQ.2) THEN
         KSCR5 = KEND1
         KEND1 = KSCR5 + NRHFT*NRHFT
      END IF

      IF (LDERIV .OR. LRELAX) THEN
         KSCR6 = KEND1
         KEND1 = KSCR6 + NVIRT*NRHFT

         IF (IOPT.EQ.2) THEN
            KSCR7 = KEND1
            KEND1 = KSCR7 + NRHFT*NRHFT
         END IF
*
         IF ((IOPT.GE.1).AND.LLAO) THEN
            KSCR8 = KEND1
            KEND1 = KSCR8 + NVIRT*NRHFT
         END IF


      END IF

      LWRK1   = LWORK - KEND1

      IF ( LWRK1 .LT. 0) THEN
        CALL QUIT('Insufficient memory in CCIAJB.')
      END IF

*---------------------------------------------------------------------*
* transform alpha index to I using XLAMDP1
*     -- store (i bet|gam del) in SCR2 
*---------------------------------------------------------------------*
      KOFF2 = KSCR2

      DO ISYMI = 1, NSYM
              
        ISYALP = MULD2H(ISYM1,ISYMI)
        ISYBET = MULD2H(ISYALP,ISY1ALBE)

        KOFF1 = 1 + IAODIS(ISYALP,ISYBET) 
        KLAMD = IGLMRH(ISYALP,ISYMI) + 1

        NBASA = MAX(NBAS(ISYALP),1)
        NBASB = MAX(NBAS(ISYBET),1)

        CALL DGEMM('T','N',NBAS(ISYBET),NRHF(ISYMI),NBAS(ISYALP),
     *             ONE,XINT1(KOFF1),NBASA,XLAMDP1(KLAMD),
     *             NBASA,ZERO,WORK(KOFF2),NBASB)

        KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)
          
      END DO

*---------------------------------------------------------------------*
* transform beta index to A using XLAMDH1 
*      -- store (i a|gam del) in SCR4 
*---------------------------------------------------------------------*
      KOFF2 = KSCR2

      DO ISYMI = 1, NSYM
              
        ISYALP = MULD2H(ISYM1,ISYMI)
        ISYBET = MULD2H(ISYALP,ISY1ALBE)
        ISYMA  = MULD2H(ISYM1,ISYBET)

        KLAMD = IGLMVI(ISYBET,ISYMA) + 1
        KOFF4 = KSCR4 + IT1AM(ISYMA,ISYMI)

        NBASB = MAX(NBAS(ISYBET),1)
        NVIRA = MAX(NVIR(ISYMA),1)

        CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYBET),
     *             ONE,XLAMDH1(KLAMD),NBASB,WORK(KOFF2),
     *             NBASB,ZERO,WORK(KOFF4),NVIRA)

        KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)

      END DO

*---------------------------------------------------------------------*
* if LRELAX transform beta index to a-bar using XLAMDH2 
*     -- store (i a-bar|gam del) in SCR6 
*---------------------------------------------------------------------*
      IF ( LRELAX ) THEN

         KOFF2 = KSCR2

         DO ISYMI = 1, NSYM
              
           ISYALP = MULD2H(ISYM1,ISYMI)
           ISYBET = MULD2H(ISYALP,ISY1ALBE)
           ISYMA  = MULD2H(ISYM2,ISYBET)
 
           KLAMD = IGLMVI(ISYBET,ISYMA) + 1
           KOFF6 = KSCR6 + IT1AM(ISYMA,ISYMI)

           NBASB = MAX(NBAS(ISYBET),1)
           NVIRA = MAX(NVIR(ISYMA),1)

           CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYBET),
     *                ONE,XLAMDH2(KLAMD),NBASB,WORK(KOFF2),
     *                NBASB,ZERO,WORK(KOFF6),NVIRA)

           KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)

         END DO

      END IF

*---------------------------------------------------------------------*
* For IOPT=2 transform beta index to j using XLAMDH1
*     -- store (i j|gam del) in SCR5
*---------------------------------------------------------------------*
      IF ( IOPT.EQ.2 ) THEN

        KOFF2 = KSCR2

        DO ISYMI = 1, NSYM
    
          ISYALP = MULD2H(ISYM1,ISYMI)
          ISYBET = MULD2H(ISYALP,ISY1ALBE)
          ISYMJ  = MULD2H(ISYM1,ISYBET)

          KLAMD = IGLMRH(ISYBET,ISYMJ) + 1
          KOFF5 = KSCR5 + IMATIJ(ISYMI,ISYMJ) 

          NBASB = MAX(NBAS(ISYBET),1)
          NRHFI = MAX(NRHF(ISYMI),1)

          CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NBAS(ISYBET),
     *               ONE,WORK(KOFF2),NBASB,XLAMDH1(KLAMD),NBASB,
     *               ZERO,WORK(KOFF5),NRHFI)

          KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)

        END DO

      END IF

*---------------------------------------------------------------------*
* If LRELAX transform beta index to j-bar using XLAMDH2
*     -- store (i j-bar| gam del) in SCR7
*---------------------------------------------------------------------*
      IF ( IOPT.EQ.2 .AND. LRELAX ) THEN

        KOFF2 = KSCR2

        DO ISYMI = 1, NSYM
    
          ISYALP = MULD2H(ISYM1,ISYMI)
          ISYBET = MULD2H(ISYALP,ISY1ALBE)
          ISYMJ  = MULD2H(ISYM2,ISYBET)

          KLAMD = IGLMRH(ISYBET,ISYMJ) + 1
          KOFF7 = KSCR7 + IMATIJ(ISYMI,ISYMJ) 

          NBASB = MAX(NBAS(ISYBET),1)
          NRHFI = MAX(NRHF(ISYMI),1)

          CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NBAS(ISYBET),
     *               ONE,WORK(KOFF2),NBASB,XLAMDH2(KLAMD),NBASB,
     *               ZERO,WORK(KOFF7),NRHFI)

          KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)

        END DO

      END IF

*---------------------------------------------------------------------*
* For LDERIV/LRELAX add extra contributions from
*         (i beta|gam del)[1] and/or (i-bar beta|gam del) :
*---------------------------------------------------------------------*
      IF ( LDERIV .OR. LRELAX ) THEN

*---------------------------------------------------------------------*
* If LDERIV add contribution from the derivative integrals 
*                                          (alp beta|gam del)[1]:
*
* transform alpha index to i using XLAMDP1 
*     -- put (i bet|gam del)[1] to in SCR2 
*---------------------------------------------------------------------*
         IF ( LDERIV ) THEN

            IF (.NOT. LRELAX) THEN
              IF (IOPT.EQ.2) CALL DZERO(WORK(KSCR7),NRHFT*NRHFT)
            END IF

            KOFF2 = KSCR2

            DO ISYMI = 1, NSYM
              
               ISYALP = MULD2H(ISYM1,ISYMI)
               ISYBET = MULD2H(ISYALP,ISY2ALBE)

               KOFF3 = 1 + IAODIS(ISYALP,ISYBET) 
               KLAMD = IGLMRH(ISYALP,ISYMI) + 1

               NBASA = MAX(NBAS(ISYALP),1)
               NBASB = MAX(NBAS(ISYBET),1)

               CALL DGEMM('T','N',NBAS(ISYBET),NRHF(ISYMI),NBAS(ISYALP),
     *                  ONE,XINT2(KOFF3),NBASA,XLAMDP1(KLAMD),
     *                  NBASA,ZERO,WORK(KOFF2),NBASB)

               KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)
            END DO
*-----------------------------------------------------------------------*
* For IOPT>=1.AND.LLAO we keep (i beta|gam del)[1] separate from 
* and only add the 3 (i a|-bar contribution together after 
* transformation of beta index 
* Now transform beta of (i beta|[1] to (i a|[1] and put in SCR8
*-----------------------------------------------------------------------*
            IF ((IOPT.GE.1).AND.(LLAO)) THEN

               KOFF2 = KSCR2

               DO ISYMI = 1, NSYM
                  ISYALP = MULD2H(ISYM1,ISYMI)
                  ISYBET = MULD2H(ISYALP,ISY2ALBE)
                  ISYMA  = MULD2H(ISYM1,ISYBET)

                  KLAMD = IGLMVI(ISYBET,ISYMA) + 1
                  KOFF8 = KSCR8 + IT1AM(ISYMA,ISYMI)
 
                  NBASB = MAX(NBAS(ISYBET),1)
                  NVIRA = MAX(NVIR(ISYMA),1)

                  CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     *                        NBAS(ISYBET),
     *                        ONE,XLAMDH1(KLAMD),NBASB,WORK(KOFF2),
     *                        NBASB,ZERO,WORK(KOFF8),NVIRA)
 
                  KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)
               END DO

            END IF

         END IF                                             !LDERIV
*----------------------------------------------------------------------*
* transform alpha index to i-bar using XLAMDP2
*     -- store (i-bar bet|gam del) in SCR2  (reuse or add)
*  NB:   For (IOPT.GE.1).AND.(LLAO) SCR2 will then only contain
*        (i-bar bet|gam del) (no derivative contribution)
*----------------------------------------------------------------------*
         IF ( LRELAX ) THEN

            IF (.NOT.LDERIV) CALL DZERO(WORK(KSCR2),NBAST*NRHFT)
            IF ((IOPT.GE.1).AND.(LLAO)) THEN
c               CALL DZERO(WORK(KSCR2),NBAST*NRHFT)
               FAC = ZERO
            ELSE
               FAC = ONE
            END IF

            KOFF2 = KSCR2
            DO ISYMI = 1, NSYM
              
               ISYALP = MULD2H(ISYM2,ISYMI)
               ISYBET = MULD2H(ISYALP,ISY1ALBE)

               KOFF1 = 1 + IAODIS(ISYALP,ISYBET) 
               KLAMD = IGLMRH(ISYALP,ISYMI) + 1

               NBASA = MAX(NBAS(ISYALP),1)
               NBASB = MAX(NBAS(ISYBET),1)

               CALL DGEMM('T','N',NBAS(ISYBET),NRHF(ISYMI),NBAS(ISYALP),
     *                     ONE,XINT1(KOFF1),NBASA,XLAMDP2(KLAMD),
     *                     NBASA,FAC,WORK(KOFF2),NBASB)

               KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)
            END DO

         END IF                                             !LRELAX
*---------------------------------------------------------------------*
* transform beta index to A using XLAMDH1
*  -- add (i-bar{+i[1]}) a|gam del) to what already SCR6,i.e. (i a-bar|
*---------------------------------------------------------------------*

         KOFF2 = KSCR2
        
         DO ISYMI = 1, NSYM

            ISYALP = MULD2H(ISYM2,ISYMI)
            ISYBET = MULD2H(ISYALP,ISY1ALBE)
            ISYMA  = MULD2H(ISYM1,ISYBET)
 
            KLAMD = IGLMVI(ISYBET,ISYMA) + 1
            KOFF6 = KSCR6 + IT1AM(ISYMA,ISYMI)

            NBASB = MAX(NBAS(ISYBET),1)
            NVIRA = MAX(NVIR(ISYMA),1)

            CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYBET),
     *                 ONE,XLAMDH1(KLAMD),NBASB,WORK(KOFF2),
     *                 NBASB,ONE,WORK(KOFF6),NVIRA)

            KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)

         END DO

*---------------------------------------------------------------------*
* IF (LLAO).AND.(IOPT.GE.1) 
*---------------------------------------------------------------------*        
*
* Add SCR6 to SCR8, (i-bar a|+(i a-bar|+(i a|[1]
*
       IF ((IOPT.GE.1).AND.(LLAO)) THEN
          DO ISYMI = 1, NSYM
             ISYALP = MULD2H(ISYM2,ISYMI)
             ISYBET = MULD2H(ISYALP,ISY1ALBE)
             ISYMA  = MULD2H(ISYM1,ISYBET)

             DO I = 1, NRHF(ISYMI)
               KOFF6 = KSCR6 + IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1)
               KOFF8 = KSCR8 + IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1)
               CALL DAXPY(NVIR(ISYMA),ONE,WORK(KOFF6),1,WORK(KOFF8),1)
             END DO
          END DO
*
* For IOPT=2, transform (i-bar beta| (in SCR2) to (i-bar j|, add to SCR7
*
          IF (IOPT.EQ.2) THEN
             KOFF2 = KSCR2
             DO ISYMI = 1, NSYM
               ISYALP = MULD2H(ISYM2,ISYMI)
               ISYBET = MULD2H(ISYALP,ISY1ALBE)
               ISYMJ  = MULD2H(ISYM1,ISYBET)

               KLAMD = IGLMRH(ISYBET,ISYMJ) + 1
               KOFF7 = KSCR7 + IMATIJ(ISYMI,ISYMJ)

               NBASB = MAX(NBAS(ISYBET),1)
               NRHFI = MAX(NRHF(ISYMI),1)

               CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NBAS(ISYBET),
     *                    ONE,WORK(KOFF2),NBASB,XLAMDH1(KLAMD),NBASB,
     *                    ONE,WORK(KOFF7),NRHFI)
               KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)
             END DO
          END IF
*
* Consider the (bet alp|del gam)[1] = -(alp bet|gam del)[1] integrals
c Transform beta of X2INT to I with XLAMDP1 (overwrite SCR2) and
* change sign ====>  -(alp i_p | gam del)[1]
*
          KOFF2 = KSCR2

          DO ISYMI = 1, NSYM
             ISYBET = MULD2H(ISYM1,ISYMI)
             ISYALP = MULD2H(ISYBET,ISY2ALBE)
             KOFF3  = IAODIS(ISYALP,ISYBET) + 1
             KLAMD  = IGLMRH(ISYBET,ISYMI)  + 1

             NBASA = MAX(NBAS(ISYALP),1)
             NBASB = MAX(NBAS(ISYBET),1)

             CALL DGEMM('N','N',NBAS(ISYALP),NRHF(ISYMI),NBAS(ISYBET),
     *                   -ONE,XINT2(KOFF3),NBASA,XLAMDP1(KLAMD),
     *                   NBASB,ZERO,WORK(KOFF2),NBASA)
             KOFF2 = KOFF2 + NBAS(ISYALP)*NRHF(ISYMI)

          END DO

*
* For IOPT >= 1, transform -alpha to A with XLAMDH1 ==> add to SCR6
* -(a_h i_p|[1]
*
          KOFF2 = KSCR2
          DO ISYMI = 1, NSYM
             ISYBET = MULD2H(ISYM1,ISYMI)
             ISYALP = MULD2H(ISYBET,ISY2ALBE)
             ISYMA  = MULD2H(ISYM1,ISYALP)

             KLAMD = IGLMVI(ISYALP,ISYMA) + 1
             KOFF6 = KSCR6 + IT1AM(ISYMA,ISYMI)

             NBASA = MAX(NBAS(ISYALP),1)
             NVIRA = MAX(NVIR(ISYMA),1)

             CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     *            NBAS(ISYALP),ONE,XLAMDH1(KLAMD),NBASA,
     *            WORK(KOFF2),NBASA,ONE,WORK(KOFF6),NVIRA)

             KOFF2 = KOFF2 + NBAS(ISYALP)*NRHF(ISYMI)
          END DO
*
* If IOPT=2 transform -alpha to J with XLAMDH1 ===> add to SCR7
* -(j_h i_p|[1], stored -(i_p j_h|[1]
*
             IF (IOPT.EQ.2) THEN
                KOFF2 = KSCR2
                DO ISYMI = 1, NSYM
                   ISYBET = MULD2H(ISYMI,ISYM1)
                   ISYALP = MULD2H(ISYBET,ISY2ALBE)
                   ISYMJ  = MULD2H(ISYM1,ISYALP)
   
                   KLAMD = IGLMRH(ISYALP,ISYMJ) + 1
                   KOFF7 = KSCR7 + IMATIJ(ISYMI,ISYMJ)
   
                   NBASA = MAX(NBAS(ISYALP),1)
                   NRHFI = MAX(NRHF(ISYMI),1)
   
                   CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),
     *                  NBAS(ISYALP),ONE,WORK(KOFF2),NBASA,
     *                  XLAMDH1(KLAMD),NBASA,ONE,WORK(KOFF7),NRHFI)

                   KOFF2 = KOFF2 + NBAS(ISYALP)*NRHF(ISYMI)
                END DO

             END IF


       END IF     
* 
* End of special module for LLAO.and.(IOPT.GE.1)
*
*---------------------------------------------------------------------*
* For IOPT=2 (NOT.LLAO) transform beta index to j using XLAMDH1
*     -- add ((i-bar+i[1]) j|gam del) add to (i j-bar| in SCR7
*---------------------------------------------------------------------*
       IF (( IOPT.EQ.2 ).AND.(.NOT.LLAO)) THEN

           KOFF2 = KSCR2
           DO ISYMI = 1, NSYM
    
              ISYALP = MULD2H(ISYM2,ISYMI)
              ISYBET = MULD2H(ISYALP,ISY1ALBE)
              ISYMJ  = MULD2H(ISYM1,ISYBET)

              KLAMD = IGLMRH(ISYBET,ISYMJ) + 1
              KOFF7 = KSCR7 + IMATIJ(ISYMI,ISYMJ) 

              NBASB = MAX(NBAS(ISYBET),1)
              NRHFI = MAX(NRHF(ISYMI),1)

              CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NBAS(ISYBET),
     *                    ONE,WORK(KOFF2),NBASB,XLAMDH1(KLAMD),NBASB,
     *                    ONE,WORK(KOFF7),NRHFI)

              KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)
           END DO
        END IF
*---------------------------------------------------------------------*
      END IF                                         !LRELAX.OR.LDERIV
*---------------------------------------------------------------------*
*     Add the contribution to the result XIAJB vector
*     transform thereby gamma to j using XLAMDP1 and 
*---------------------------------------------------------------------*

      IF ( LZERO ) THEN
C        -------------------------
C        add (i a|j del) to XIAJB:
C        -------------------------
         ISYMAI = MULD2H(ISY1ALBE,ISYM11)

         CALL CC_IAJB1(IGAM, WORK(KSCR4), ISYMAI, ISYGAM,
     &                 XLAMDP1, ISYM1, XIAJB)
      END IF


      IF ( LRELAX ) THEN
C        ------------------------------
C        add (i a|j-bar del) to X2IAJB:
C        ------------------------------
         ISYMAI = MULD2H(ISY1ALBE,ISYM11)

         CALL CC_IAJB1(IGAM, WORK(KSCR4), ISYMAI, ISYGAM,
     &                 XLAMDP2, ISYM2, X2IAJB)
      END IF

      IF ( LDERIV .OR. LRELAX ) THEN
C        ------------------------------------------------
C        add (i-bar a|j del) + (i a-bar|j del) to X2IAJB:
C        ------------------------------------------------
         ISYMAI = MULD2H(ISY1ALBE,ISYM12)

         IF ((IOPT.GE.1).AND.(LLAO)) THEN
            KSCR68 = KSCR8
         ELSE
            KSCR68 = KSCR6
         END IF

         CALL CC_IAJB1(IGAM, WORK(KSCR68), ISYMAI, ISYGAM,
     &                 XLAMDP1, ISYM1, X2IAJB)
      END IF

*---------------------------------------------------------------------*
*     Add the contributions to the (ia|del j) integrals in XIABJ vector
*     and to the (ia|del j)-bar integrals in X2IABJ vector
*     transform thereby gamma to j/j-bar using XLAMDH1/XLAMDH2 
*---------------------------------------------------------------------*
      IF ( IOPT.EQ.1 .OR. IOPT.EQ.2 ) THEN

         IF ( LZERO ) THEN
C           -------------------------
C           add (i a|del j) to XIABJ:
C           -------------------------
            ISYMAI = MULD2H(ISY1ALBE,ISYM11)

            CALL CC_IAJB1(IGAM, WORK(KSCR4), ISYMAI, ISYGAM,
     &                    XLAMDH1, ISYM1, XIABJ)
         END IF


         IF ( LRELAX ) THEN
C           ------------------------------
C           add (i a|del j-bar) to X2IAJB:
C           ------------------------------
            ISYMAI = MULD2H(ISY1ALBE,ISYM11)

            CALL CC_IAJB1(IGAM, WORK(KSCR4), ISYMAI, ISYGAM,
     &                    XLAMDH2, ISYM2, X2IABJ)
         END IF

         IF ( LDERIV .OR. LRELAX ) THEN
C           ------------------------------------------------
C           add (i-bar a|del j) + (i a-bar|del j) to X2IABJ:
C           ------------------------------------------------
            ISYMAI = MULD2H(ISY1ALBE,ISYM12)

            CALL CC_IAJB1(IGAM, WORK(KSCR6), ISYMAI, ISYGAM,
     &                                XLAMDH1, ISYM1, X2IABJ)
         END IF

      END IF
*---------------------------------------------------------------------*
*     Add the contributions to the result XIJBA/X2IJBA vector
*     transform thereby gamma to a/a-bar using XLAMDH1/XLAMDH2 
*---------------------------------------------------------------------*
      IF ( IOPT.EQ.2 ) THEN

         IF ( LZERO ) THEN
C           -------------------------
C           add (i j|del a) to XIJBA:
C           -------------------------
            ISYMIJ = MULD2H(ISY1ALBE,ISYM11)

            CALL CC_IJBA1(IGAM, WORK(KSCR5), ISYMIJ, ISYGAM,
     &                    XLAMDH1, ISYM1, XIJBA)
         END IF


         IF ( LRELAX ) THEN
C           ------------------------------
C           add (i j|del a-bar) to X2IJBA:
C           ------------------------------
            ISYMIJ = MULD2H(ISY1ALBE,ISYM11)

            CALL CC_IJBA1(IGAM, WORK(KSCR5), ISYMIJ, ISYGAM,
     &                    XLAMDH2, ISYM2, X2IJBA)
         END IF

         IF ( LDERIV .OR. LRELAX ) THEN
C           ------------------------------------------------
C           add (i-bar j|del a) + (i j-bar|del a) to X2IJBA:
C           ------------------------------------------------
            ISYMIJ = MULD2H(ISY1ALBE,ISYM12)

            CALL CC_IJBA1(IGAM, WORK(KSCR7), ISYMIJ, ISYGAM,
     &                    XLAMDH1, ISYM1, X2IJBA)
         END IF

      END IF 

      RETURN
      END
*=====================================================================*
*                 END OF SUBROUTINE CCIAJB                            *
*=====================================================================*
