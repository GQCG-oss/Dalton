*=====================================================================*
C  /* Deck cc_sort4o */
      SUBROUTINE CC_SORT4O(XKJIL,ISYM4O,XJKLI,IOPT)
*---------------------------------------------------------------------*
*     Purpose: resort I(kj,i;l) to I(jk,l;i)
*
*     IOPT = 1 : XJKLI area is initialized here
*     IOPT = 2 : XJKLI is added to gamma intermediate
*                       already stored as I(jk,l;i)
*
*     Sonia Coriani, 14/09-1999
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#  include <ccorb_dec.h>
#  include <maxorb_dec.h>
#  include <ccsdsym_dec.h>
#else
#  include <implicit.h>
#endif
#include <ccorb.h>
#include <maxorb.h>
#include <ccsdsym.h>
      INTEGER ISYM4O, IOPT
#if defined (SYS_CRAY)
      REAL XKJIL(*), XJKLI(*)
      REAL ZERO, ONE, HALF, DDOT, XNORM, FAC
#else
      DOUBLE PRECISION XKJIL(*), XJKLI(*)
      DOUBLE PRECISION ZERO, ONE, HALF, DDOT, XNORM, FAC
#endif
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0)
*
      INTEGER ISYML, ISYKJI, ISYMI, ISYMKJ, ISYJKL, ISYMJK
      INTEGER ISYMK, ISYMJ,  KJIL,  JKLI
*
*     --------------------------------------
*     Initialize result area with zero's or add to previous
*     --------------------------------------
      IF (IOPT.EQ.1) THEN
        FAC = ZERO
      ELSE
        FAC = ONE 
      END IF
*     --------------------------------------
*     Reorder thru loops on all 4 indices
*     --------------------------------------
      DO ISYML = 1,NSYM
         ISYKJI = MULD2H(ISYM4O,ISYML)
         DO L = 1, NRHF(ISYML)
            DO ISYMI = 1, NSYM
               ISYMKJ = MULD2H(ISYKJI,ISYMI)
               ISYJKL = MULD2H(ISYM4O,ISYMI)
               ISYMJK = MULD2H(ISYJKL,ISYML)
               DO I = 1, NRHF(ISYMI)
                  DO ISYMK = 1, NSYM
                     ISYMJ = MULD2H(ISYMJK,ISYMK)
                     DO K = 1, NRHF(ISYMK)
                     DO J = 1, NRHF(ISYMJ)

                KJIL = I3ORHF(ISYKJI,ISYML) + NMAIJK(ISYKJI)*(L-1)+
     &                 IMAIJK(ISYMKJ,ISYMI) + NMATIJ(ISYMKJ)*(I-1)+
     &                 IMATIJ(ISYMK,ISYMJ)  + NRHF(ISYMK)*(J-1) + K

                JKLI = I3ORHF(ISYJKL,ISYMI) + NMAIJK(ISYJKL)*(I-1)+
     &                 IMAIJK(ISYMJK,ISYML) + NMATIJ(ISYMJK)*(L-1)+
     &                 IMATIJ(ISYMJ,ISYMK)  + NRHF(ISYMJ)*(K-1) + J

               XJKLI(JKLI) = FAC*XJKLI(JKLI) + XKJIL(KJIL)

                     END DO         !J
                     END DO         !K
                  END DO            !ISYMK
               END DO               !I
            END DO                  !ISYMI
         END DO                     !L
      END DO                        !ISYML
*     ---------------------------------------
*     Finished, return
*     ---------------------------------------
      RETURN
      END
*=====================================================================*
