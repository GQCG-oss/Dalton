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
c /* deck cc_ijcb */
*=====================================================================*
       SUBROUTINE CC_IJCB( XINT1,   ISY1ALBE,
     &                     XINT2,   ISY2ALBE, 
     &                     IDEL,    IGAM, 
     &                     X1IJCB,  X1CJIB,
     &                     X2IJCB,  X2CJIB,
     &                     XLAMDP1, XLAMDH1, ISYML1,
     &                     XLAMDP2, XLAMDH2, ISYML2,
     &                     XLAMDP3, XLAMDH3, ISYML3,
     &                     XLAMDP4, XLAMDH4, ISYML4,
     &                     WORK,    LWORK,   
     &                     IOPT,    LDERIV,  LRELAX, LZERO,
     &                     LNEWTA,  LX2ISQ )
*---------------------------------------------------------------------*
*
* Purpose: drive the generalized transformation to 
*          (ij^|cb) + (ij|c^b) integrals and
*          (cj^|ib) + (c^j|ib) integrals
*          for the two-index (**|gam del) approach
*          assumes three-index arrays XIJCB & XCJIB in core
*
*          this routine transforms the indices ij and c, the 
*          transformation of the delta index to b has to be done
*          from the outside.  
*
*        XINT1 zero order AO integrals
*        XINT2 derivative AO integrals
*        XA1IJCB,  XA1CJIB : TA dependent (zero ord) integrals
*        XA2IJCB,  XA2CJIB : TA dependent deriv & relax integrals
*
*        IOPT=0: (ij^|c del) + (ij|c^ del)  only
*
*        IOPT=1: (ij^|c del) + (ij|c^ del) and 
*                (cj^|i del) + (c^j|i del) integrals
*
*        IF LDERIV=.TRUE. transform also the derivative integrals g[1]:
*                           with the XLAMD_1 and XLAMD_3
*           IOPT=0: (ij^|c del) + (ij|c^ del)
*           IOPT=1: (ij^|c del) + (ij|c^ del) and (cj^|i del) + (c^j|i del)
*
*        IF LRELAX=.TRUE. include relaxation contribution to the
*                         derivative integrals from the transformation 
*                         of g[0] with XLAMD_1 * XLAMD_2 * XLAMD_3 * XLAMD_4
*                         (or just reorthonormalization if IRELAX = 0)
*
*        IF LZERO=.FALSE. skip calculation of zero-order integrals
*        IF LX2ISQ = .TRUE. the (al bet| part of X2INT is already full matr.
*
*    Written by Sonia Coriani, February 1999
*    based  on Christof's CC_IAJB
*    Note that no special case is needed for LAO apart LX2ISQ = TRUE
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <maxorb.h>
#include <ccisao.h>

#if defined (SYS_CRAY)
      REAL ONE, ZERO
#else
      DOUBLE PRECISION ONE, ZERO
#endif
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0)

      LOGICAL LDERIV, LRELAX, LZERO, LNEWTA, LX2ISQ
      INTEGER IDEL, IGAM, ISY1ALBE, ISY2ALBE
      INTEGER ISYML1, ISYML2, ISYML3, ISYML4, LWORK, IOPT
      
#if defined (SYS_CRAY)
      REAL XLAMDP1(*), XLAMDH1(*)
      REAL XLAMDP2(*), XLAMDH2(*)
      REAL XLAMDP3(*), XLAMDH3(*)
      REAL XLAMDP4(*), XLAMDH4(*)
      REAL XINT1(*),  XINT2(*)
      REAL X1IJCB(*), X1CJIB(*)
      REAL X2IJCB(*), X2CJIB(*)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION XLAMDP1(*), XLAMDH1(*)
      DOUBLE PRECISION XLAMDP2(*), XLAMDH2(*)
      DOUBLE PRECISION XLAMDP3(*), XLAMDH3(*)
      DOUBLE PRECISION XLAMDP4(*), XLAMDH4(*)
      DOUBLE PRECISION XINT1(*), XINT2(*)
      DOUBLE PRECISION X1IJCB(*), X1CJIB(*)
      DOUBLE PRECISION X2IJCB(*), X2CJIB(*)
      DOUBLE PRECISION WORK(LWORK)
#endif

* local 

      INTEGER ISYL11, ISYL12, ISYL13, ISYL14, ISYL23
      INTEGER ISYDEL, ISYGAM
      INTEGER KSCR1, KSCR3, KEND1, LWRK1, KDUM
      PARAMETER (KDUM = + 999 999 999)

*---------------------------------------------------------------------*
*     set some symmetries
*---------------------------------------------------------------------*

      ISYDEL  = ISAO(IDEL)
      ISYGAM  = ISAO(IGAM)

*---------------------------------------------------------------------*
*   work space allocation:
*
*   KSCR1  --  I^{del,gam}(alp bet)
*
*   KSCR3  --  I^{del,gam}(alp bet)[1]
*
*---------------------------------------------------------------------*
      KSCR1   = 1
      KEND1   = KSCR1  + N2BST(ISY1ALBE)           

c      IF ((LDERIV .OR. LRELAX).AND.(.NOT.LX2ISQ)) THEN

      IF ((LDERIV).AND.(.NOT.LX2ISQ)) THEN
        KSCR3  = KEND1
        KEND1  = KSCR3  + N2BST(ISY2ALBE)
      ELSE
        KSCR3  = KDUM
      END IF

      LWRK1   = LWORK - KEND1

      IF ( LWRK1 .LT. 0) THEN
        CALL QUIT('Insufficient memory in CC_IJCB.')
      END IF

*---------------------------------------------------------------------*
*     square zero-order integral matrix up (alp and bet)
*     and derivative also if (LDERIV).AND.(.NOT.LX2ISQ)
*---------------------------------------------------------------------*

      CALL CCSD_SYMSQ(XINT1,ISY1ALBE,WORK(KSCR1))

      IF ((LDERIV).AND.(.NOT.LX2ISQ)) THEN
        CALL CCSD_SYMSQ(XINT2,ISY2ALBE,WORK(KSCR3))
      END IF

*---------------------------------------------------------------------*
*     call routine for actual calculation of transformed integrals
*---------------------------------------------------------------------*
      IF (.NOT.LX2ISQ) THEN
         CALL CCIJCB(WORK(KSCR1),ISY1ALBE,WORK(KSCR3),ISY2ALBE,
     &              IDEL, IGAM, X1IJCB,  X1CJIB, X2IJCB,  X2CJIB,
     &              XLAMDP1, XLAMDH1, ISYML1, XLAMDP2, XLAMDH2, ISYML2,
     &              XLAMDP3, XLAMDH3, ISYML3, XLAMDP4, XLAMDH4, ISYML4,
     &              WORK(KEND1),  LWRK1, 
     &              IOPT,    LDERIV,  LRELAX, LZERO, LNEWTA )
      ELSE
         CALL CCIJCB(WORK(KSCR1),ISY1ALBE,XINT2,ISY2ALBE,
     &              IDEL, IGAM, X1IJCB,  X1CJIB, X2IJCB,  X2CJIB,
     &              XLAMDP1, XLAMDH1, ISYML1, XLAMDP2, XLAMDH2, ISYML2,
     &              XLAMDP3, XLAMDH3, ISYML3, XLAMDP4, XLAMDH4, ISYML4,
     &              WORK(KEND1),  LWRK1, 
     &              IOPT,    LDERIV,  LRELAX, LZERO, LNEWTA )
      END IF

*---------------------------------------------------------------------*
*     return
*---------------------------------------------------------------------*

      RETURN
      END
*=====================================================================*
*                 END OF SUBROUTINE CC_IJCB                           *
*=====================================================================*
c /* deck ccijcb */
*=====================================================================*
       SUBROUTINE CCIJCB( X1INT,   ISY1ALBE,
     &                    X2INT,   ISY2ALBE, 
     &                    IDEL,    IGAM, 
     &                    X1IJCB,  X1CJIB,
     &                    X2IJCB,  X2CJIB,
     &                    XLAMDP1, XLAMDH1, ISYML1,
     &                    XLAMDP2, XLAMDH2, ISYML2,
     &                    XLAMDP3, XLAMDH3, ISYML3,
     &                    XLAMDP4, XLAMDH4, ISYML4,
     &                    WORK,    LWORK,   
     &                    IOPT,    LDERIV,  LRELAX, LZERO,
     &                    LNEWTA )
*---------------------------------------------------------------------*
*
* Purpose: realise the generalized transformation to 
*          (ij^|cb) + (ij|c^b) integrals and
*          (cj^|ib) + (c^j|ib) integrals
*          for the two-index (**|gam del) approach
*          assumes three-index arrays XIJCB & XCJIB in core
*
* See CC_IJCB for details
*
*    Written by Sonia Coriani, February 1999
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <maxorb.h>
#include <ccisao.h>

#if defined (SYS_CRAY)
      REAL ONE, ZERO
#else
      DOUBLE PRECISION ONE, ZERO
#endif
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0)

      LOGICAL LDERIV, LRELAX, LZERO, LNEWTA
      INTEGER IDEL, IGAM, ISY1ALBE, ISY2ALBE, ISYML1, ISYML2
      INTEGER ISYML3, ISYML4, LWORK, IOPT, IDUMMY
      
#if defined (SYS_CRAY)
      REAL XLAMDP1(*), XLAMDH1(*)
      REAL XLAMDP2(*), XLAMDH2(*)
      REAL XLAMDP3(*), XLAMDH3(*)
      REAL XLAMDP4(*), XLAMDH4(*)
      REAL X1INT(*), X2INT(*)
      REAL X1IJCB(*), X1CJIB(*)
      REAL X2IJCB(*), X2CJIB(*)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION XLAMDP1(*), XLAMDH1(*)
      DOUBLE PRECISION XLAMDP2(*), XLAMDH2(*)
      DOUBLE PRECISION XLAMDP3(*), XLAMDH3(*)
      DOUBLE PRECISION XLAMDP4(*), XLAMDH4(*)
      DOUBLE PRECISION X1INT(*), X2INT(*)
      DOUBLE PRECISION X1IJCB(*), X1CJIB(*)
      DOUBLE PRECISION X2IJCB(*), X2CJIB(*)
      DOUBLE PRECISION WORK(LWORK)
#endif

* local 

      INTEGER ISYL11, ISYL12, ISYL13, ISYL14, ISYL23
      INTEGER ISYALP, ISYBET, ISYDEL, ISYGAM
      INTEGER ISYALP1, ISYBET1
      INTEGER ISYMI, ISYMJ, ISYMIJ, ISYMC
      INTEGER KSCR2, KSCR4, KSCR41, KSCR5
      INTEGER KSCR6, KSCR61, KSCR7, KEND1, LWRK1
      INTEGER KOFF1, KOFF2, KOFF3, KOFF4, KOFF41, KLAMD
      INTEGER KOFF5, KOFF6, KOFF61, KOFF7
      INTEGER NBASA, NBASB, NVIRC, NRHFI
      INTEGER ISYIJ1,ISYIJ2,ISYIJ3,ISYIJ4
      INTEGER ISYCJ2,ISYCJ4,ISYCJ3

*---------------------------------------------------------------------*
*     set some symmetries
*---------------------------------------------------------------------*

      ISYDEL  = ISAO(IDEL)
      ISYGAM  = ISAO(IGAM)

*---------------------------------------------------------------------*
*   work space allocation:
*
*   KSCR2  --  I^{del,gam}(alp j); (alp j^); (alp j-bar); (alp j^-bar)
*   KSCR4  --  I^{del,gam}(i j)
*   KSCR41 --  I^{del,gam}(i j^)
*
*   KSCR6  --  I^{del,gam}([i-bar j + i j-bar]) + I^[1]^{del,gam}(i j)
*   KSCR61 --  I^{del,gam}([i-bar j^ + i j^-bar]) + I^[1]^{del,gam}(i j^)
*
*   KSCR5 --  I^{del,gam}([c^ j  + c j^])
*   KSCR7 --  I^{del,gam}([c^-bar j + c-bar j^ + c j^-bar + c^ j-bar])
*             + I^[1]^{del,gam}(c^ j + c j^)
*
*---------------------------------------------------------------------*
      KSCR2   = 1
      KSCR4   = KSCR2  + NBAST*NRHFT
      KSCR41  = KSCR4  + NRHFT*NRHFT
      KEND1   = KSCR41 + NRHFT*NRHFT

      IF (IOPT.EQ.1) THEN
        KSCR5 = KEND1
        KEND1 = KSCR5 + NVIRT*NRHFT
      END IF

      IF (LDERIV .OR. LRELAX) THEN
        KSCR6  = KEND1
        KSCR61 = KSCR6  + NRHFT*NRHFT
        KEND1  = KSCR61 + NRHFT*NRHFT

        IF (IOPT.EQ.1) THEN
          KSCR7 = KEND1
          KEND1 = KSCR7 + NVIRT*NRHFT
        END IF
      END IF

      LWRK1   = LWORK - KEND1

      IF ( LWRK1 .LT. 0) THEN
        CALL QUIT('Insufficient memory in CC_IJCB.')
      END IF

*---------------------------------------------------------------------*
*     transform beta index to J using XLAMDH1
*      -- store (alf J|gam del)  in SCR2 
*---------------------------------------------------------------------*
c**      
c      ISYSCR2 = MULD2H(ISY1ALBE,ISYML1)
c      CALL DZERO(NT1AO(ISYSCR2),WORK(KSCR2))
c**
      KOFF2 = KSCR2
      DO ISYMJ = 1, NSYM
              
        ISYBET = MULD2H(ISYML1,ISYMJ)
        ISYALP = MULD2H(ISYBET,ISY1ALBE)

        KOFF1 = 1 + IAODIS(ISYALP,ISYBET) 
        KLAMD = IGLMRH(ISYBET,ISYMJ) + 1

        NBASA = MAX(NBAS(ISYALP),1)
        NBASB = MAX(NBAS(ISYBET),1)

        CALL DGEMM('N','N',NBAS(ISYALP),NRHF(ISYMJ),NBAS(ISYBET),
     *             ONE,X1INT(KOFF1),NBASA,XLAMDH1(KLAMD),
     *             NBASB,ZERO,WORK(KOFF2),NBASA)

        KOFF2 = KOFF2 + NBAS(ISYALP)*NRHF(ISYMJ)
          
      END DO
*---------------------------------------------------------------------*
*     transform alpha index to I using XLAMDP1 
*      -- store (i j|gam del) in SCR4 
*---------------------------------------------------------------------*
      KOFF2 = KSCR2
      DO ISYMJ = 1, NSYM
              
        ISYBET = MULD2H(ISYML1,ISYMJ)
        ISYALP = MULD2H(ISYBET,ISY1ALBE)
        ISYMI  = MULD2H(ISYML1,ISYALP)

        KLAMD = IGLMRH(ISYALP,ISYMI) + 1
        KOFF4 = KSCR4 + IMATIJ(ISYMI,ISYMJ)

        NBASA = MAX(NBAS(ISYALP),1)
        NRHFI = MAX(NRHF(ISYMI),1)

        CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NBAS(ISYALP),
     *             ONE,XLAMDP1(KLAMD),NBASA,WORK(KOFF2),
     *             NBASA,ZERO,WORK(KOFF4),NRHFI)

        KOFF2 = KOFF2 + NBAS(ISYALP)*NRHF(ISYMJ)

      END DO

*---------------------------------------------------------------------*
* LRELAX   transform alpha index to i-bar using XLAMDP2
*          -- store (i-bar j|gam del) in SCR6
*---------------------------------------------------------------------*
      IF ( LRELAX ) THEN

         KOFF2 = KSCR2

         DO ISYMJ = 1, NSYM

           ISYBET = MULD2H(ISYML1,ISYMJ)
           ISYALP = MULD2H(ISYBET,ISY1ALBE)
           ISYMI  = MULD2H(ISYML2,ISYALP)

           KLAMD = IGLMRH(ISYALP,ISYMI) + 1
           KOFF6 = KSCR6 + IMATIJ(ISYMI,ISYMJ)

           NBASA = MAX(NBAS(ISYALP),1)
           NRHFI = MAX(NRHF(ISYMI),1)

           CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),
     *                NBAS(ISYALP),ONE,XLAMDP2(KLAMD),NBASA,
     *                WORK(KOFF2),NBASA,ZERO,WORK(KOFF6),NRHFI)

           KOFF2 = KOFF2 + NBAS(ISYALP)*NRHF(ISYMJ)

         END DO

      END IF

*---------------------------------------------------------------------*
*     for IOPT=1 transform alpha index to C^ using XLAMDP3
*        -- store (c^ j|gam del) in SCR5
*---------------------------------------------------------------------*
      IF ( IOPT.EQ.1 ) THEN

        KOFF2 = KSCR2

        DO ISYMJ = 1, NSYM

          ISYBET = MULD2H(ISYML1,ISYMJ)
          ISYALP = MULD2H(ISYBET,ISY1ALBE)
          ISYMC  = MULD2H(ISYML3,ISYALP)

          KLAMD = IGLMVI(ISYALP,ISYMC) + 1
          KOFF5 = KSCR5 + IT1AM(ISYMC,ISYMJ)

          NBASA = MAX(NBAS(ISYALP),1)
          NVIRC = MAX(NVIR(ISYMC),1)

          CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMJ),NBAS(ISYALP),
     *               ONE,XLAMDP3(KLAMD),NBASA,WORK(KOFF2),NBASA,
     *               ZERO,WORK(KOFF5),NVIRC)

          KOFF2 = KOFF2 + NBAS(ISYALP)*NRHF(ISYMJ)
c
        END DO

      END IF
*---------------------------------------------------------------------*
*     for IOPT=1 and LRELAX transform alpha index to C^(bar) 
*     using XLAMDP4 -- store (c^-bar j|gam del) in SCR7
*---------------------------------------------------------------------*
      IF (( IOPT.EQ.1 ).AND.LRELAX) THEN

        KOFF2 = KSCR2

        DO ISYMJ = 1, NSYM

          ISYBET = MULD2H(ISYML1,ISYMJ)
          ISYALP = MULD2H(ISYBET,ISY1ALBE)
          ISYMC = MULD2H(ISYML4,ISYALP)

          KLAMD = IGLMVI(ISYALP,ISYMC) + 1
          KOFF7 = KSCR7 + IT1AM(ISYMC,ISYMJ)

          NBASA = MAX(NBAS(ISYALP),1)
          NVIRC = MAX(NVIR(ISYMC),1)

          CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMJ),NBAS(ISYALP),
     *               ONE,XLAMDP4(KLAMD),NBASA,WORK(KOFF2),NBASA,
     *               ZERO,WORK(KOFF7),NVIRC)

          KOFF2 = KOFF2 + NBAS(ISYALP)*NRHF(ISYMJ)

        END DO

      END IF

*---------------------------------------------------------------------*
*     transform beta index to J^ using XLAMDH3
*      -- store (alf J^|gam del)  in SCR2 
*---------------------------------------------------------------------*
      KOFF2 = KSCR2

      DO ISYMJ = 1, NSYM
              
        ISYBET = MULD2H(ISYML3,ISYMJ)
        ISYALP = MULD2H(ISYBET,ISY1ALBE)

        KOFF1 = 1 + IAODIS(ISYALP,ISYBET) 
        KLAMD = IGLMRH(ISYBET,ISYMJ) + 1

        NBASA = MAX(NBAS(ISYALP),1)
        NBASB = MAX(NBAS(ISYBET),1)

        CALL DGEMM('N','N',NBAS(ISYALP),NRHF(ISYMJ),NBAS(ISYBET),
     *             ONE,X1INT(KOFF1),NBASA,XLAMDH3(KLAMD),
     *             NBASB,ZERO,WORK(KOFF2),NBASA)

        KOFF2 = KOFF2 + NBAS(ISYALP)*NRHF(ISYMJ)
          
      END DO

*---------------------------------------------------------------------*
*     transform alpha index to I using XLAMDP1 
*      -- store (i j^|gam del) in SCR41 
*---------------------------------------------------------------------*
      KOFF2 = KSCR2

      DO ISYMJ = 1, NSYM
              
        ISYBET = MULD2H(ISYML3,ISYMJ)
        ISYALP = MULD2H(ISYBET,ISY1ALBE)
        ISYMI  = MULD2H(ISYML1,ISYALP)

        KLAMD  = IGLMRH(ISYALP,ISYMI) + 1
        KOFF41 = KSCR41 + IMATIJ(ISYMI,ISYMJ)

        NBASA = MAX(NBAS(ISYALP),1)
        NRHFI = MAX(NRHF(ISYMI),1)

        CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NBAS(ISYALP),
     *             ONE,XLAMDP1(KLAMD),NBASA,WORK(KOFF2),
     *             NBASA,ZERO,WORK(KOFF41),NRHFI)

        KOFF2 = KOFF2 + NBAS(ISYALP)*NRHF(ISYMJ)

      END DO

*---------------------------------------------------------------------*
*     if LRELAX transform alpha index to i-bar using XLAMDP2
*        -- store (i-bar j^| gam del) in SCR61
*---------------------------------------------------------------------*
      IF ( LRELAX ) THEN

        KOFF2 = KSCR2

        DO ISYMJ = 1, NSYM
    
          ISYBET = MULD2H(ISYML3,ISYMJ)
          ISYALP = MULD2H(ISYBET,ISY1ALBE)
          ISYMI  = MULD2H(ISYML2,ISYALP)

          KLAMD  = IGLMRH(ISYALP,ISYMI) + 1
          KOFF61 = KSCR61 + IMATIJ(ISYMI,ISYMJ) 

          NBASA = MAX(NBAS(ISYALP),1)
          NRHFI = MAX(NRHF(ISYMI),1)

          CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NBAS(ISYALP),
     *               ONE,XLAMDP2(KLAMD),NBASA,WORK(KOFF2),NBASA,
     *               ZERO,WORK(KOFF61),NRHFI)

          KOFF2 = KOFF2 + NBAS(ISYALP)*NRHF(ISYMJ)

        END DO

      END IF

*---------------------------------------------------------------------*
*     for IOPT=1 transform alpha index to C using XLAMDP1
*     -- add (c j^|gam del) to (c^ j|gam del) in SCR5
*---------------------------------------------------------------------*
      IF ( IOPT.EQ.1 ) THEN

        KOFF2 = KSCR2

        DO ISYMJ = 1, NSYM
    
          ISYBET = MULD2H(ISYML3,ISYMJ)
          ISYALP = MULD2H(ISYBET,ISY1ALBE)
          ISYMC  = MULD2H(ISYML1,ISYALP)

          KLAMD = IGLMVI(ISYALP,ISYMC) + 1
          KOFF5 = KSCR5 + IT1AM(ISYMC,ISYMJ)

          NBASA = MAX(NBAS(ISYALP),1)
          NVIRC = MAX(NVIR(ISYMC),1)

          CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMJ),NBAS(ISYALP),
     *               ONE,XLAMDP1(KLAMD),NBASA,WORK(KOFF2),NBASA,
     *               ONE,WORK(KOFF5),NVIRC)

          KOFF2 = KOFF2 + NBAS(ISYALP)*NRHF(ISYMJ)

        END DO
      END IF
*---------------------------------------------------------------------*
*     for IOPT=1.and.LRELAX  transform alpha to C-bar using XLAMDP2
*     -- add (c-bar j^|gam del) to (c^-bar j|gam del) in SCR7
*---------------------------------------------------------------------*
      IF (( IOPT.EQ.1 ).AND.(LRELAX)) THEN

        KOFF2 = KSCR2

        DO ISYMJ = 1, NSYM
    
          ISYBET = MULD2H(ISYML3,ISYMJ)
          ISYALP = MULD2H(ISYBET,ISY1ALBE)
          ISYMC  = MULD2H(ISYML2,ISYALP)

          KLAMD = IGLMVI(ISYALP,ISYMC) + 1
          KOFF7 = KSCR7 + IT1AM(ISYMC,ISYMJ)

          NBASA = MAX(NBAS(ISYALP),1)
          NVIRC = MAX(NVIR(ISYMC),1)

          CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMJ),NBAS(ISYALP),
     *               ONE,XLAMDP2(KLAMD),NBASA,WORK(KOFF2),NBASA,
     *               ONE,WORK(KOFF7),NVIRC)

          KOFF2 = KOFF2 + NBAS(ISYALP)*NRHF(ISYMJ)

        END DO

      END IF
C
C Finished with SCR2 again
C
*---------------------------------------------------------------------*
*     for LRELAX add extra contributions from (alp j-bar|gam del),
*     (alp ^j-bar|gam del) etc
*---------------------------------------------------------------------*
      IF ( LRELAX ) THEN
*---------------------------------------------------------------------*
*     transform beta to J-bar using XLAMDH2
*     -- store (alp j-bar|gam del) in SCR2 
*     If (LDERIV) add also (alp j|gam del)[1]
*---------------------------------------------------------------------*

         KOFF2 = KSCR2

         DO ISYMJ = 1, NSYM
              
            ISYBET = MULD2H(ISYML2,ISYMJ)
            ISYALP = MULD2H(ISYBET,ISY1ALBE)

            KOFF1 = 1 + IAODIS(ISYALP,ISYBET) 
            KLAMD = IGLMRH(ISYBET,ISYMJ) + 1

            NBASA = MAX(NBAS(ISYALP),1)
            NBASB = MAX(NBAS(ISYBET),1)

            CALL DGEMM('N','N',NBAS(ISYALP),NRHF(ISYMJ),NBAS(ISYBET),
     *                 ONE,X1INT(KOFF1),NBASA,XLAMDH2(KLAMD),
     *                 NBASB,ZERO,WORK(KOFF2),NBASA)

            IF (LDERIV) THEN

              ISYBET1 = MULD2H(ISYML1,ISYMJ)
              ISYALP1 = MULD2H(ISYBET1,ISY2ALBE)
              IF (ISYALP1.NE.ISYALP) 
     *          CALL QUIT('Symmetry mismatch in CC_IJCB')

              KOFF3  = 1 + IAODIS(ISYALP1,ISYBET1) 
              KLAMD  = IGLMRH(ISYBET1,ISYMJ) + 1
              NBASA = MAX(NBAS(ISYALP1),1)
              NBASB = MAX(NBAS(ISYBET1),1)
              CALL DGEMM('N','N',NBAS(ISYALP1),NRHF(ISYMJ),
     *                  NBAS(ISYBET1),ONE,X2INT(KOFF3),NBASA,
     *                  XLAMDH1(KLAMD),NBASB,ONE,WORK(KOFF2),NBASA)
            END IF


            KOFF2 = KOFF2 + NBAS(ISYALP)*NRHF(ISYMJ)
          
          END DO

*---------------------------------------------------------------------*
*     transform alpha to I using XLAMDP1
*     -- add (i j-bar|gam del) to (i-bar j|gam del) in SCR6
*---------------------------------------------------------------------*

          KOFF2 = KSCR2

          DO ISYMJ = 1, NSYM

             ISYBET = MULD2H(ISYML2,ISYMJ)
             ISYALP = MULD2H(ISYBET,ISY1ALBE)
             ISYMI  = MULD2H(ISYML1,ISYALP)

             KLAMD = IGLMRH(ISYALP,ISYMI) + 1
             KOFF6 = KSCR6 + IMATIJ(ISYMI,ISYMJ)

             NBASA = MAX(NBAS(ISYALP),1)
             NRHFI = MAX(NRHF(ISYMI),1)

             CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),
     *                  NBAS(ISYALP),ONE,XLAMDP1(KLAMD),NBASA,
     *                  WORK(KOFF2), NBASA,ONE,WORK(KOFF6),NRHFI)

             KOFF2 = KOFF2 + NBAS(ISYALP)*NRHF(ISYMJ)

           END DO
*---------------------------------------------------------------------*
*          for IOPT=1  transform alpha to C^ using XLAMDP3
*          -- add (c^ j-bar | gam del) in SCR7
*---------------------------------------------------------------------*
           IF ( IOPT.EQ.1 ) THEN

           KOFF2 = KSCR2

           DO ISYMJ = 1, NSYM

             ISYBET = MULD2H(ISYML2,ISYMJ)
             ISYALP = MULD2H(ISYBET,ISY1ALBE)
             ISYMC  = MULD2H(ISYML3,ISYALP)

             KLAMD = IGLMVI(ISYALP,ISYMC) + 1
             KOFF7 = KSCR7 + IT1AM(ISYMC,ISYMJ)

             NBASA = MAX(NBAS(ISYALP),1)
             NVIRC = MAX(NVIR(ISYMC),1)

             CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMJ),
     *               NBAS(ISYALP),ONE,XLAMDP3(KLAMD),NBASA,
     *               WORK(KOFF2),NBASA,ONE,WORK(KOFF7),NVIRC)

             KOFF2 = KOFF2 + NBAS(ISYALP)*NRHF(ISYMJ)

           END DO
           END IF

*---------------------------------------------------------------------*
*     transform beta to J^-bar using XLAMDH4
*     -- store (alp j^-bar|gam del) in SCR2
* if (LDERIV) add derivative contribution (al j^|gam del)(1)
*---------------------------------------------------------------------*

         KOFF2 = KSCR2

         DO ISYMJ = 1, NSYM

            ISYBET = MULD2H(ISYML4,ISYMJ)
            ISYALP = MULD2H(ISYBET,ISY1ALBE)

            KOFF1 = 1 + IAODIS(ISYALP,ISYBET)
            KLAMD = IGLMRH(ISYBET,ISYMJ) + 1

            NBASA = MAX(NBAS(ISYALP),1)
            NBASB = MAX(NBAS(ISYBET),1)

            CALL DGEMM('N','N',NBAS(ISYALP),NRHF(ISYMJ),NBAS(ISYBET),
     *                 ONE,X1INT(KOFF1),NBASA,XLAMDH4(KLAMD),
     *                 NBASB,ZERO,WORK(KOFF2),NBASA)

            IF (LDERIV) THEN

              ISYBET1 = MULD2H(ISYML3,ISYMJ)
              ISYALP1 = MULD2H(ISYBET1,ISY2ALBE)
              IF (ISYALP1.NE.ISYALP) 
     *           CALL QUIT('Symmetry mismatch in CC_IJCB')

              KOFF3  = 1 + IAODIS(ISYALP1,ISYBET1)
              KLAMD  = IGLMRH(ISYBET1,ISYMJ) + 1
              NBASA = MAX(NBAS(ISYALP1),1)
              NBASB = MAX(NBAS(ISYBET1),1)
              CALL DGEMM('N','N',NBAS(ISYALP1),NRHF(ISYMJ),
     *                  NBAS(ISYBET1),ONE,X2INT(KOFF3),NBASA,
     *                  XLAMDH3(KLAMD),NBASB,ONE,WORK(KOFF2),NBASA)
            END IF

            KOFF2 = KOFF2 + NBAS(ISYALP)*NRHF(ISYMJ)

          END DO
*---------------------------------------------------------------------*
*     transform alpha to I using XLAMDP1
*     -- add (i j^-bar|gam del) to (i-bar j^|gam del) in SCR61
*---------------------------------------------------------------------*

          KOFF2 = KSCR2

          DO ISYMJ = 1, NSYM

             ISYBET = MULD2H(ISYML4,ISYMJ)
             ISYALP = MULD2H(ISYBET,ISY1ALBE)
             ISYMI  = MULD2H(ISYML1,ISYALP)

             KLAMD  = IGLMRH(ISYALP,ISYMI) + 1
             KOFF61 = KSCR61 + IMATIJ(ISYMI,ISYMJ)

             NBASA = MAX(NBAS(ISYALP),1)
             NRHFI = MAX(NRHF(ISYMI),1)

             CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),
     *                 NBAS(ISYALP),ONE,XLAMDP1(KLAMD),NBASA,
     *                 WORK(KOFF2),NBASA,ONE,WORK(KOFF61),NRHFI)

             KOFF2 = KOFF2 + NBAS(ISYALP)*NRHF(ISYMJ)

           END DO
*---------------------------------------------------------------------*
*          for IOPT=1  transform alpha to C using XLAMDP1
*          -- add (c j^-bar | gam del) in SCR7
*---------------------------------------------------------------------*
           IF ( IOPT.EQ.1 ) THEN

           KOFF2 = KSCR2

           DO ISYMJ = 1, NSYM

             ISYBET = MULD2H(ISYML4,ISYMJ)
             ISYALP = MULD2H(ISYBET,ISY1ALBE)
             ISYMC  = MULD2H(ISYML1,ISYALP)

             KLAMD = IGLMVI(ISYALP,ISYMC) + 1
             KOFF7 = KSCR7 + IT1AM(ISYMC,ISYMJ)

             NBASA = MAX(NBAS(ISYALP),1)
             NVIRC = MAX(NVIR(ISYMC),1)

             CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMJ),
     *               NBAS(ISYALP),ONE,XLAMDP1(KLAMD),NBASA,
     *               WORK(KOFF2),NBASA,ONE,WORK(KOFF7),NVIRC)

             KOFF2 = KOFF2 + NBAS(ISYALP)*NRHF(ISYMJ)

           END DO
           END IF
      END IF
*---------------------------------------------------------------------*
*     Add the contribution to the result X1IJCB and X2IJCB vector:
*     (transform the gamma index to VIRTUAL)
*     Note that IJCB integrals are sorted as I_cj,i,del
*---------------------------------------------------------------------*
      ISYL11 = MULD2H(ISYML1,ISYML1)
      ISYL12 = MULD2H(ISYML1,ISYML2)
      ISYL13 = MULD2H(ISYML1,ISYML3)
      ISYL14 = MULD2H(ISYML1,ISYML4)
      ISYL23 = MULD2H(ISYML2,ISYML3)

      IF ( LNEWTA ) THEN
C        -------------------------
C        add (ij|c^ del) to X1IJCB:
C        -------------------------
         ISYIJ1 = MULD2H(ISY1ALBE,ISYL11)

         CALL CC_IJCB2(IGAM, WORK(KSCR4), ISYIJ1, ISYGAM,
     &                 XLAMDP3, ISYML3, X1IJCB)
C        -------------------------
C        add (ij^|c del) to X1IJCB:
C        -------------------------
         ISYIJ3 = MULD2H(ISY1ALBE,ISYL13)
         CALL CC_IJCB2(IGAM, WORK(KSCR41), ISYIJ3, ISYGAM,
     &                 XLAMDP1, ISYML1, X1IJCB)
      END IF


      IF ( LRELAX ) THEN
C        ------------------------------
C        add (ij|c^-bar del) to X2IJCB:
C        ------------------------------
         ISYIJ1 = MULD2H(ISY1ALBE,ISYL11)

         CALL CC_IJCB2(IGAM, WORK(KSCR4), ISYIJ1, ISYGAM,
     &                 XLAMDP4, ISYML4, X2IJCB)
C        ------------------------------
C        add (ij^|c-bar del) to X2IJCB:
C        ------------------------------                 
         ISYIJ3 = MULD2H(ISY1ALBE,ISYL13)

         CALL CC_IJCB2(IGAM, WORK(KSCR41), ISYIJ3, ISYGAM,
     &                 XLAMDP2, ISYML2, X2IJCB)
      END IF
      IF ( LDERIV .OR. LRELAX ) THEN
C        ------------------------------------------------
C        add ([i-bar j^ + i j^-bar]|c del) to X2IJCB:
C        ------------------------------------------------
         ISYIJ4 = MULD2H(ISY1ALBE,ISYL14)

         CALL CC_IJCB2(IGAM, WORK(KSCR61), ISYIJ4, ISYGAM,
     &                 XLAMDP1, ISYML1, X2IJCB)
C        ------------------------------------------------
C        add ([i-bar j + i j-bar]|c^ del) to X2IJCB:
C        ------------------------------------------------
         ISYIJ2 = MULD2H(ISY1ALBE,ISYL12)

         CALL CC_IJCB2(IGAM, WORK(KSCR6), ISYIJ2, ISYGAM,
     &                 XLAMDP3, ISYML3, X2IJCB)
      END IF

*---------------------------------------------------------------------*
*     Add the contribution to the result X1CJIB and X2CJIB vector:
*     transform gamma to OCCUPIED
*     Integrals sorted as I_cj,i,del
*---------------------------------------------------------------------*
      IF ( IOPT.EQ.1 ) THEN

         IF ( LNEWTA ) THEN
C           -------------------------
C           add ([c^ j + c j^]|i del) to X1CJIB:
C           -------------------------
c***
            ISYCJ3 = MULD2H(ISY1ALBE,ISYL13)

            CALL CC_IAJB1(IGAM, WORK(KSCR5), ISYCJ3, ISYGAM,
     &                    XLAMDP1, ISYML1, X1CJIB, .FALSE., IDUMMY)
         END IF


         IF ( LRELAX ) THEN
C           ------------------------------
C           add  ([c^ j + c j^]|i-bar del) to X2CJIB:
C           ------------------------------
            ISYCJ3 = MULD2H(ISY1ALBE,ISYL13)

            CALL CC_IAJB1(IGAM, WORK(KSCR5), ISYCJ3, ISYGAM,
     &                    XLAMDP2, ISYML2, X2CJIB, .FALSE., IDUMMY)
         END IF

         IF ( LDERIV .OR. LRELAX ) THEN
C           ------------------------------------------------
C           add ([c-bar j^ + c j^-bar + c^-bar j + c^ j-bar] | i del) to X2IABJ:
C           ------------------------------------------------
            ISYCJ4 = MULD2H(ISY1ALBE,ISYL14)

            CALL CC_IAJB1(IGAM, WORK(KSCR7), ISYCJ4, ISYGAM,
     &                    XLAMDP1, ISYML1, X2CJIB, .FALSE., IDUMMY)
         END IF

      END IF

      RETURN
      END
*=====================================================================*
*                 END OF SUBROUTINE CCIJCB                            *
*=====================================================================*
c /* deck cc_ijcb2 */
*=====================================================================*
      SUBROUTINE CC_IJCB2(IGAM, XIJG, ISYMIJ, ISYGAM, 
     &                    XLAMDA, ISYLAM, XIJCB  )
*---------------------------------------------------------------------*
*
*   Purpose: transform (ij|gam del) to (ij|c del), with sorting
*            I_cj,i,del
*   Sonia, March 1999
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <ccsdsym.h>
#include <ccorb.h>
#include <maxorb.h>
#include <ccisao.h>

#if defined (SYS_CRAY)
      REAL XIJG(*), XIJCB(*), XLAMDA(*)
#else
      DOUBLE PRECISION XIJG(*), XIJCB(*), XLAMDA(*)
#endif

      INTEGER IGAM, ISYMIJ, ISYGAM, ISYLAM, ISYMC, KLAMD, KRES
      INTEGER ISYMCJ, ISYMI, ISYMJ, KOFF5, NBASG
      INTEGER KXIJ

* transform integral batch:
      ISYMC  = MULD2H(ISYGAM,ISYLAM)
      G      = IGAM - IBAS(ISYGAM)
      NBASG  = MAX(NBAS(ISYGAM),1)

      DO ISYMJ = 1, NSYM
         ISYMI  = MULD2H(ISYMIJ,ISYMJ)
         ISYMCJ = MULD2H(ISYMC,ISYMJ)

         DO J = 1, NRHF(ISYMJ)
           DO I = 1, NRHF(ISYMI)

              KLAMD = IGLMVI(ISYGAM,ISYMC) + G               !offs Lambda_gam,c
              !offset I_ij,gam
              KXIJ  = IMATIJ(ISYMI,ISYMJ)  + NRHF(ISYMI)*(J-1) + I
              !offset I_cj,i
              KRES  = IT2BCD(ISYMCJ,ISYMI) + NT1AM(ISYMCJ)*(I-1) 
     &                + IT1AM(ISYMC,ISYMJ) + NVIR(ISYMC)*(J-1) + 1
              
              CALL DAXPY(NVIR(ISYMC),XIJG(KXIJ),XLAMDA(KLAMD),
     &                                        NBASG,XIJCB(KRES),1)
           END DO
         END DO

      END DO

      RETURN
      END
*=====================================================================*
*                 END OF SUBROUTINE CC_IJCB2                          *
*=====================================================================*
