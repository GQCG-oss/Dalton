c /* deck cciabj0 */
*======================================================================*
       SUBROUTINE CCIAJB0( XINT1,   ISY1ALBE,
     &                     XINT2,   ISY2ALBE, 
     &                     IDEL,    IGAM, 
     &                     X1IA, X2GDIA, X2DGIA, X1IJ, X2IJ,
     &                     XLAMDP1, XLAMDH1, ISYM1,
     &                     XLAMDP2, XLAMDH2, ISYM2,
     &                     WORK,    LWORK,   IOPT,
     &                     LDERIV,  LRELAX,  LX2ISQ, IREAL )
*----------------------------------------------------------------------*
*
*   Purpose: perform two first quarter transformations for 
*            for the two-index (**|gam del) approach:
*
*            X1IA   = (i a|gam del)^(0)
*            X2GDIA = (i a|gam del)^(1)
*            X2DGIA = (i a|del gam)^(1) 
*            X1IJ   = (i j|gam del)^(0) 
*            X2IJ   = (i j|gam del)^(1)
*
*            iopt = 0 : only X1IA and X2GDIA
*            iopt = 1 : in addiation X2DGIA if needed 
*            iopt = 2 : in addiation X1IJ and X2IJ
*
*            X2DGIA is only needed if LDERIV is set and IREAL=-1
*
*            LRELAX : include relaxation/reorthogonalization contrib.
*                     from XLAMDP2/XLAMDH2 matrices
*
*            LDERIV : include contributions from derivative integrals
*                     passed on the XINT2 array
*
*            LX2ISQ : derivative integrals are stored in squared form
*                     (mandatory if IREAL=-1 )
*
*            IREAL = +1 : derivative integrals are real
*            IREAL = -1 : derivative integrals are pure imaginary
*
*    Written by Christof Haettig, May 1998.
*    Restructured by Sonia Coriani to handle LAO's, October 1999
*    First two quarter transformation put into own subroutine for
*    CC2 option in CC_XIETA and better readabiliy, Ch. H., spring 2000
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

      LOGICAL LDERIV, LRELAX, LAO, LX2ISQ
      INTEGER IDEL, IGAM, ISY1ALBE, ISY2ALBE, ISYM1, ISYM2, LWORK, IOPT
      INTEGER IREAL
      
#if defined (SYS_CRAY)
      REAL XLAMDP1(*), XLAMDH1(*)
      REAL XLAMDP2(*), XLAMDH2(*)
      iEAL XINT1(*), XINT2(*)
      REAL X1IA(*), X2GDIA(*), X2DGIA(*), X1IJ(*), X2IJ(*)
      REAL WORK(LWORK)
      REAL FAC
#else
      DOUBLE PRECISION XLAMDP1(*), XLAMDH1(*)
      DOUBLE PRECISION XLAMDP2(*), XLAMDH2(*)
      DOUBLE PRECISION XINT1(*), XINT2(*)
      DOUBLE PRECISION X1IA(*), X2GDIA(*), X2DGIA(*), X1IJ(*), X2IJ(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION FAC
#endif

      INTEGER ISYM11, ISYM12, ISYMJ, ISYMB, ISYRES1, ISYRES2
      INTEGER ISYMI, ISYMA, ISYALP, ISYBET, ISYDEL, ISYGAM, LENHLF
      INTEGER KLAMD, KOFF1, KOFF2, KOFF3, KOFF4, KHALF1, KHALF2, KXAO
      INTEGER NBASA, NBASB, NVIRA, KOFF6, KEND1, LWRK1
      INTEGER NRHFI, KOFF7, KOFF5, KOFF8, ISYHLF

*----------------------------------------------------------------------*
* set some symmetries and flags, check work space:
*----------------------------------------------------------------------*
      ISYM11  = MULD2H(ISYM1,ISYM1)
      ISYM12  = MULD2H(ISYM1,ISYM2)

      ISYDEL  = ISAO(IDEL)
      ISYGAM  = ISAO(IGAM)

      ISYRES1 = MULD2H(ISYM11,ISY1ALBE)

      IF ( LRELAX .AND. LDERIV ) THEN
        IF ( MULD2H(ISYM12,ISY1ALBE) .NE. MULD2H(ISYM11,ISY2ALBE) ) THEN
          CALL QUIT('Symmetry mismatch in CC_IAJB0.')
        END IF
      END IF

      IF (LRELAX) ISYRES2 = MULD2H(ISYM12,ISY1ALBE)
      IF (LDERIV) ISYRES2 = MULD2H(ISYM11,ISY2ALBE)

      IF      (IREAL.EQ.+1 .OR. (.NOT.LDERIV)) THEN
        LAO = .FALSE.
      ELSE IF (IREAL.EQ.-1 .AND. LDERIV) THEN
        LAO = .TRUE.
      ELSE
        CALL QUIT('ILLEGAL VALUE OF IREAL IN CC_IAJB0.')
      END IF                                               

      ! maximum dimension of half-transformed integrals
      LENHLF  =            NT1AO(MULD2H(ISYM1,ISY1ALBE))
      LENHLF  = MAX(LENHLF,NT1AO(MULD2H(ISYM2,ISY1ALBE)))
      LENHLF  = MAX(LENHLF,NT1AO(MULD2H(ISYM1,ISY2ALBE)))

      ! scratch area for squared AO integrals
      KXAO  = 1
      KEND1 = KXAO + N2BST(ISY1ALBE)
      IF (LDERIV .AND. (.NOT.LX2ISQ)) THEN
        KEND1 = KXAO + MAX(N2BST(ISY1ALBE),N2BST(ISY2ALBE))
      END IF

      ! scratch area for 1-index transformed integrals
      KHALF1 = KEND1
      KEND1  = KHALF1 + LENHLF

      ! special scratch area only needed for imaginary integrals
      IF (LAO) THEN
        KHALF2 = KEND1
        KEND1  = KHALF2 + LENHLF
      END IF

      LWRK1  = LWORK  - KEND1
      IF ( LWRK1 .LT. 0) THEN
        CALL QUIT('Insufficient memory in CCIAJB0.')
      END IF

*----------------------------------------------------------------------*
* square up usual integrals and store in WORK(KXAO)
*----------------------------------------------------------------------*

      CALL CCSD_SYMSQ(XINT1,ISY1ALBE,WORK(KXAO))

*----------------------------------------------------------------------*
* transform alpha index to i using XLAMDP1
*     -- store (i bet|gam del) in WORK(KHALF1) 
*----------------------------------------------------------------------*

      KOFF2 = KHALF1

      DO ISYMI = 1, NSYM
              
        ISYALP = MULD2H(ISYM1,ISYMI)
        ISYBET = MULD2H(ISYALP,ISY1ALBE)

        KOFF1 = KXAO + IAODIS(ISYALP,ISYBET) 
        KLAMD = IGLMRH(ISYALP,ISYMI) + 1

        NBASA = MAX(NBAS(ISYALP),1)
        NBASB = MAX(NBAS(ISYBET),1)

        CALL DGEMM('T','N',NBAS(ISYBET),NRHF(ISYMI),NBAS(ISYALP),
     *             ONE,WORK(KOFF1),NBASA,XLAMDP1(KLAMD),
     *             NBASA,ZERO,WORK(KOFF2),NBASB)

        KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)
          
      END DO

*----------------------------------------------------------------------*
* transform beta index to A using XLAMDH1 
*      -- store (i a|gam del) in X1IA 
*----------------------------------------------------------------------*
      KOFF2 = KHALF1

      DO ISYMI = 1, NSYM
              
        ISYALP = MULD2H(ISYM1,ISYMI)
        ISYBET = MULD2H(ISYALP,ISY1ALBE)
        ISYMA  = MULD2H(ISYM1,ISYBET)

        KLAMD = IGLMVI(ISYBET,ISYMA) + 1
        KOFF4 = IT1AM(ISYMA,ISYMI) + 1

        NBASB = MAX(NBAS(ISYBET),1)
        NVIRA = MAX(NVIR(ISYMA),1)

        CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYBET),
     *             ONE,XLAMDH1(KLAMD),NBASB,WORK(KOFF2),
     *             NBASB,ZERO,X1IA(KOFF4),NVIRA)

        KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)

      END DO

*----------------------------------------------------------------------*
* if LRELAX transform beta index to a-bar using XLAMDH2 
*     -- store (i a-bar|gam del) in X2GDIA 
*----------------------------------------------------------------------*
      IF ( LRELAX ) THEN

         KOFF2 = KHALF1

         DO ISYMI = 1, NSYM
              
           ISYALP = MULD2H(ISYM1,ISYMI)
           ISYBET = MULD2H(ISYALP,ISY1ALBE)
           ISYMA  = MULD2H(ISYM2,ISYBET)
 
           KLAMD = IGLMVI(ISYBET,ISYMA) + 1
           KOFF6 = IT1AM(ISYMA,ISYMI) + 1

           NBASB = MAX(NBAS(ISYBET),1)
           NVIRA = MAX(NVIR(ISYMA),1)

           CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYBET),
     *                ONE,XLAMDH2(KLAMD),NBASB,WORK(KOFF2),
     *                NBASB,ZERO,X2GDIA(KOFF6),NVIRA)

           KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)

         END DO

      END IF

*----------------------------------------------------------------------*
* For IOPT=2 transform beta index to j using XLAMDH1
*     -- store (i j|gam del) in X1IJ
*----------------------------------------------------------------------*
      IF ( IOPT.EQ.2 ) THEN

        KOFF2 = KHALF1

        DO ISYMI = 1, NSYM
    
          ISYALP = MULD2H(ISYM1,ISYMI)
          ISYBET = MULD2H(ISYALP,ISY1ALBE)
          ISYMJ  = MULD2H(ISYM1,ISYBET)

          KLAMD = IGLMRH(ISYBET,ISYMJ) + 1
          KOFF5 = IMATIJ(ISYMI,ISYMJ) + 1

          NBASB = MAX(NBAS(ISYBET),1)
          NRHFI = MAX(NRHF(ISYMI),1)

          CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NBAS(ISYBET),
     *               ONE,WORK(KOFF2),NBASB,XLAMDH1(KLAMD),NBASB,
     *               ZERO,X1IJ(KOFF5),NRHFI)

          KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)

        END DO

      END IF

*----------------------------------------------------------------------*
* If LRELAX transform beta index to j-bar using XLAMDH2
*     -- store (i j-bar| gam del) in X2IJ
*----------------------------------------------------------------------*
      IF ( IOPT.EQ.2 .AND. LRELAX ) THEN

        KOFF2 = KHALF1

        DO ISYMI = 1, NSYM
    
          ISYALP = MULD2H(ISYM1,ISYMI)
          ISYBET = MULD2H(ISYALP,ISY1ALBE)
          ISYMJ  = MULD2H(ISYM2,ISYBET)

          KLAMD = IGLMRH(ISYBET,ISYMJ) + 1
          KOFF7 = IMATIJ(ISYMI,ISYMJ) + 1

          NBASB = MAX(NBAS(ISYBET),1)
          NRHFI = MAX(NRHF(ISYMI),1)

          CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NBAS(ISYBET),
     *               ONE,WORK(KOFF2),NBASB,XLAMDH2(KLAMD),NBASB,
     *               ZERO,X2IJ(KOFF7),NRHFI)

          KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)

        END DO

      END IF

*----------------------------------------------------------------------*
* For LDERIV/LRELAX add extra contributions from
*         (i beta|gam del)[1] and/or (i-bar beta|gam del) :
*----------------------------------------------------------------------*
      IF ( LDERIV .OR. LRELAX ) THEN

         IF (.NOT. LRELAX) THEN
           CALL DZERO(X2GDIA,NT1AM(ISYRES2))
           IF (IOPT.EQ.2) CALL DZERO(X2IJ,NMATIJ(ISYRES2))
         END IF

         IF (LAO .AND. IOPT.GE.1) THEN
           CALL DCOPY(NT1AM(ISYRES2),X2GDIA,1,X2DGIA,1)
         END IF
*----------------------------------------------------------------------*
* transform alpha index to i-bar using XLAMDP2
*     -- store (i-bar bet|gam del) in WORK(KHALF1)
*        if LAO flag set save an extra copy in WORK(KHALF2)
* (We initialize here WORK(KHALF1) and WORK(KHALF2) if LAO flag is set)
*----------------------------------------------------------------------*
         IF ( LRELAX ) THEN

            KOFF2 = KHALF1

            DO ISYMI = 1, NSYM
           
               ISYALP = MULD2H(ISYM2,ISYMI)
               ISYBET = MULD2H(ISYALP,ISY1ALBE)

               KOFF1 = KXAO + IAODIS(ISYALP,ISYBET) 
               KLAMD = IGLMRH(ISYALP,ISYMI) + 1

               NBASA = MAX(NBAS(ISYALP),1)
               NBASB = MAX(NBAS(ISYBET),1)

               CALL DGEMM('T','N',NBAS(ISYBET),NRHF(ISYMI),NBAS(ISYALP),
     *                     ONE,WORK(KOFF1),NBASA,XLAMDP2(KLAMD),
     *                     NBASA,ZERO,WORK(KOFF2),NBASB)

               KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)
            END DO

         ELSE
            ISYHLF = MULD2H(ISY1ALBE,ISYM2)
            CALL DZERO(WORK(KHALF1),NT1AO(ISYHLF))
         END IF  

         IF (LAO) THEN
            ISYHLF = MULD2H(ISY1ALBE,ISYM2)
            CALL DCOPY(NT1AO(ISYHLF),WORK(KHALF1),1,WORK(KHALF2),1)
         END IF

*----------------------------------------------------------------------*
* If LDERIV add contribution from the derivative integrals:
* transform alpha index to i using XLAMDP1 
*     -- put (i bet|gam del)[1] + (i-bar bet|gam del) in WORK(KHALF1) 
*----------------------------------------------------------------------*
         IF ( LDERIV ) THEN

            IF (.NOT. LX2ISQ) THEN
              CALL CCSD_SYMSQ(XINT2,ISY2ALBE,WORK(KXAO))
            END IF

            KOFF2 = KHALF1

            DO ISYMI = 1, NSYM
              
              ISYALP = MULD2H(ISYM1,ISYMI)
              ISYBET = MULD2H(ISYALP,ISY2ALBE)

              KLAMD = IGLMRH(ISYALP,ISYMI) + 1

              NBASA = MAX(NBAS(ISYALP),1)
              NBASB = MAX(NBAS(ISYBET),1)

              IF (LX2ISQ) THEN
               KOFF3 = 1 + IAODIS(ISYALP,ISYBET) 
               CALL DGEMM('T','N',NBAS(ISYBET),NRHF(ISYMI),NBAS(ISYALP),
     *                    ONE,XINT2(KOFF3),NBASA,XLAMDP1(KLAMD),
     *                    NBASA,ONE,WORK(KOFF2),NBASB)
              ELSE
               KOFF3 = KXAO + IAODIS(ISYALP,ISYBET) 
               CALL DGEMM('T','N',NBAS(ISYBET),NRHF(ISYMI),NBAS(ISYALP),
     *                    ONE,WORK(KOFF3),NBASA,XLAMDP1(KLAMD),
     *                    NBASA,ONE,WORK(KOFF2),NBASB)
              END IF

              KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)
            END DO

         END IF

*----------------------------------------------------------------------*
* transform beta index to A using XLAMDH1
*  -- add (i-bar{+i[1]}) a|gam del) to what is already in X2GDIA
*----------------------------------------------------------------------*

         KOFF2 = KHALF1
        
         DO ISYMI = 1, NSYM

            ISYALP = MULD2H(ISYM2,ISYMI)
            ISYBET = MULD2H(ISYALP,ISY1ALBE)
            ISYMA  = MULD2H(ISYM1,ISYBET)
 
            KLAMD = IGLMVI(ISYBET,ISYMA) + 1
            KOFF6 = IT1AM(ISYMA,ISYMI) + 1

            NBASB = MAX(NBAS(ISYBET),1)
            NVIRA = MAX(NVIR(ISYMA),1)

            CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYBET),
     *                 ONE,XLAMDH1(KLAMD),NBASB,WORK(KOFF2),
     *                 NBASB,ONE,X2GDIA(KOFF6),NVIRA)

            KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)

         END DO

*----------------------------------------------------------------------*
* If LAO flag set and IOPT.GE.1 calculate and put int WORK(KHALF1) now
* (i alp|del gam)^(1) = -(alp i|gam del)[1] + (i-bar alp|gam del)
* i.e. transform beta index to i using XLAMDP1
*     -- add to what is in KHALF2 and put result in KHALF1
*----------------------------------------------------------------------*
         IF (LAO .AND. IOPT.GE.1) THEN

           IF (LAO .AND. (.NOT. LX2ISQ) ) THEN
             CALL QUIT('Illegal option combination in CC_IAJB0.')
           END IF

           ISYHLF = MULD2H(ISY1ALBE,ISYM2)
           CALL DCOPY(NT1AO(ISYHLF),WORK(KHALF2),1,WORK(KHALF1),1)

           KOFF2 = KHALF1

           DO ISYMI = 1, NSYM
              ISYBET = MULD2H(ISYM1,ISYMI)
              ISYALP = MULD2H(ISYBET,ISY2ALBE)
              KOFF3  = 1 + IAODIS(ISYALP,ISYBET) 
              KLAMD  = IGLMRH(ISYBET,ISYMI)  + 1
 
              NBASA = MAX(NBAS(ISYALP),1)
              NBASB = MAX(NBAS(ISYBET),1)
              CALL DGEMM('N','N',NBAS(ISYALP),NRHF(ISYMI),NBAS(ISYBET),
     *                   -ONE,XINT2(KOFF3),NBASA,XLAMDP1(KLAMD),NBASB,
     *                    ONE,WORK(KOFF2),NBASA)
              KOFF2 = KOFF2 + NBAS(ISYALP)*NRHF(ISYMI)

           END DO

         END IF

*----------------------------------------------------------------------*
* If LAO flag set and IOPT.GE.1 complete (i a|del gam)[1] integrals:
* transform alpha index to a using XLAMDH1 
*     -- add (i a|del gam)[1] + (i-bar a|del gam) to X2DGIA
* (if LAO flag not set X2DGIA = X2GDIA, and we don't need to compute it)
*----------------------------------------------------------------------*
         IF (LAO .AND. IOPT.GE.1) THEN

           KOFF2 = KHALF1
           DO ISYMI = 1, NSYM
              ISYBET = MULD2H(ISYM1,ISYMI)
              ISYALP = MULD2H(ISYBET,ISY2ALBE)
              ISYMA  = MULD2H(ISYM1,ISYALP)
 
              KLAMD = IGLMVI(ISYALP,ISYMA) + 1
              KOFF6 = IT1AM(ISYMA,ISYMI) + 1

              NBASA = MAX(NBAS(ISYALP),1)
              NVIRA = MAX(NVIR(ISYMA),1)

              CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     *             NBAS(ISYALP),ONE,XLAMDH1(KLAMD),NBASA,
     *             WORK(KOFF2),NBASA,ONE,X2DGIA(KOFF6),NVIRA)
 
              KOFF2 = KOFF2 + NBAS(ISYALP)*NRHF(ISYMI)
            END DO

          END IF

*----------------------------------------------------------------------*
* For IOPT=2 transform beta index to j using XLAMDH1
*     -- add ((i-bar+i[1]) j|del gam) add to (i j-bar| in X2IJ
*----------------------------------------------------------------------*
       IF ( IOPT.EQ.2 ) THEN

           KOFF2 = KHALF1
           DO ISYMI = 1, NSYM
    
              ISYALP = MULD2H(ISYM2,ISYMI)
              ISYBET = MULD2H(ISYALP,ISY1ALBE)
              ISYMJ  = MULD2H(ISYM1,ISYBET)

              KLAMD = IGLMRH(ISYBET,ISYMJ) + 1
              KOFF7 = IMATIJ(ISYMI,ISYMJ) + 1

              NBASB = MAX(NBAS(ISYBET),1)
              NRHFI = MAX(NRHF(ISYMI),1)

              CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NBAS(ISYBET),
     *                    ONE,WORK(KOFF2),NBASB,XLAMDH1(KLAMD),NBASB,
     *                    ONE,X2IJ(KOFF7),NRHFI)

              KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)
           END DO

        END IF

*----------------------------------------------------------------------*
      END IF                                         !LRELAX.OR.LDERIV

      RETURN
      END
*======================================================================*
*                 END OF SUBROUTINE CCIAJB                             *
*======================================================================*
