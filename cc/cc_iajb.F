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
c /* deck cc_iabj */
*======================================================================*
       SUBROUTINE CC_IAJB( X1INT,   ISY1ALBE,
     &                     X2INT,   ISY2ALBE, 
     &                     IDEL, IGAM, LAUXG, IBASX,
     &                     XIAJB,   XIABJ,   XIJBA,  
     &                     X2IAJB,  X2IABJ,  X2IJBA, 
     &                     XLAMDP1, XLAMDH1, ISYM1,
     &                     XLAMDP2, XLAMDH2, ISYM2,
     &                     XLMDPB1, XLMDHB1, ISYB1,
     &                     XLMDPB2, XLMDHB2, ISYB2,
     &                     WORK,    LWORK,   IOPT,    
     &                     LDERIV,  LRELAX,  LZERO,   
     &                     LTRIANG, LX2ISQ,  IREAL )
*----------------------------------------------------------------------*
*
*   Purpose: generalized transformation to (ia|jb), (ia|bj), and 
*            (ij|ab) integrals for the two-index (**|gam del) approach
*            assumes three-index arrays XIAJB, XIABJ, XIJBA in core
*
*            this routine drives the transformation of the indices 
*            ia and j, the transformation of the delta index to b has 
*            to be done from the outside.  
*
*            X1INT, XIAJB,  XIABJ,  XIJBA  : usual integrals
*            X2INT, X2IAJB, X2IABJ, X2IJBA : derivative integrals
*
*            XLAMDP1, XLAMDH1, XLAMDP2, XLAMDH2 are used for the first
*            two quarter transformations, the Lambda matrices in
*            XLMDPB1, XLMDHB1, XLMDPB2, XLMDHB2 for the 3. quartertrans.
*            
*            
*            IOPT=0: (ia|j del) only
*
*            IOPT=1: (ia|j del) and (ia|del j) integrals
*
*            IOPT=2: (ia|j del), (ia|del j), and (ij|del a) integrals
*
*            IOPT=3: (ia|del j), and (ij|del a) integrals
*
*            IF LDERIV=.TRUE. transform also the derivative integrals:
*                IOPT=0: (ia|j del)-bar
*                IOPT=1: (ia|j del)-bar, (ia|del j)-bar
*                IOPT=2: (ia|j del)-bar, (ia|del j)-bar, (ij|del a)-bar 
*
*                If LX2ISQ=.TRUE. The derivative integrals are already 
*                                 squared in input
*                if IREAL = -1 The derivative integrals are LAO 
*                               i.e. antisymmetric (and already squared)
*                if IREAL = +1 derivative integrals assumed symmetric
*
*            IF LZERO=.FALSE. skip calculation of zero-order integrals
*                              (ia|j del), (ia|del j), and (ij|del a)
*
*            IF LRELAX=.TRUE. include relaxation contribution to the
*                derivative integrals using XLAMDP2, XLAMDH2 matrices:
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
*            IF LTRIANG=.TRUE. and iopt=0 compute only those (ia|j del) 
*                integrals needed for lower triangle of (ia|jb)
*                     
*
*            IF LX2ISQ=.TRUE. derivative integrals are assumed to be
*              stored in squared form (default is packed)
*
*            IF LAUXG=.TRUE. substract IBASX(ISYM) when calculating
*              from IGAM the index within the symmetry class (irrep)
*
*            (ia|jb) integrals:
*                  i            transform. with XLAMDP1 with sym. ISYM1
*                  a            transform. with XLAMDH1 with sym. ISYM1
*                  j            transform. with XLMDPB1 with sym. ISYB1
*                  i-bar        transform. with XLAMDP2 with sym. ISYM2
*                  a-bar        transform. with XLAMDH2 with sym. ISYM2
*                  j-bar        transform. with XLMDPB2 with sym. ISYB2
*             
*            (ia|bj) integrals:
*                  i            transform. with XLAMDP1 with sym. ISYM1
*                  a            transform. with XLAMDH1 with sym. ISYM1
*                  j            transform. with XLMDHB1 with sym. ISYB1
*                  i-bar        transform. with XLAMDP2 with sym. ISYM2
*                  a-bar        transform. with XLAMDH2 with sym. ISYM2
*                  j-bar        transform. with XLMDHB2 with sym. ISYB2
*             
*            (ij|ba) integrals:
*                  i            transform. with XLAMDP1 with sym. ISYM1
*                  j            transform. with XLAMDH1 with sym. ISYM1
*                  a            transform. with XLMDBH1 with sym. ISYB1
*                  i-bar        transform. with XLAMDP2 with sym. ISYM2
*                  j-bar        transform. with XLAMDH2 with sym. ISYM2
*                  a-bar        transform. with XLMDBH2 with sym. ISYB2
*             
*
*    Written by Christof Haettig, May 1998.
*
*    Restructured by Sonia Coriani, October 1999 in order to handle
*    squared X2INT in input and LAO's
*    First two quarter transformations put into extra subroutine and
*    separate Lambda matrices for third quarter transformation 
*    introduced in spring 2000, Ch. H.
*======================================================================*
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

      LOGICAL LDERIV, LRELAX, LZERO, LTRIANG, LX2ISQ, LLAO, LAUXG
      INTEGER IDEL,IGAM,ISY1ALBE,ISY2ALBE,ISYM1,ISYM2,LWORK,IOPT,IREAL
      INTEGER ISYB1, ISYB2, IBASX(*)
      
#if defined (SYS_CRAY)
      REAL XLAMDP1(*), XLAMDH1(*), XLMDPB1(*), XLMDHB1(*)
      REAL XLAMDP2(*), XLAMDH2(*), XLMDPB2(*), XLMDHB2(*)
      REAL X1INT(*), X2INT(*)
      REAL XIAJB(*),  XIABJ(*),  XIJBA(*)
      REAL X2IAJB(*), X2IABJ(*), X2IJBA(*)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION XLAMDP1(*), XLAMDH1(*), XLMDPB1(*), XLMDHB1(*)
      DOUBLE PRECISION XLAMDP2(*), XLAMDH2(*), XLMDPB2(*), XLMDHB2(*)
      DOUBLE PRECISION X1INT(*), X2INT(*)
      DOUBLE PRECISION XIAJB(*),  XIABJ(*),  XIJBA(*)
      DOUBLE PRECISION X2IAJB(*), X2IABJ(*), X2IJBA(*)
      DOUBLE PRECISION WORK(LWORK)
#endif

      INTEGER ISYM11, ISYM12, ISYRES1, ISYRES2, ISYGAM, ISYDEL
      INTEGER KX1IA, KX2GDIA, KX2DGIA, KX1IJ, KX2IJ, KEND1, LWRK1
      INTEGER ISYMAI, ISYMBJ

*---------------------------------------------------------------------*
*     set some symmetries
*---------------------------------------------------------------------*
      ISYM11  = MULD2H(ISYM1,ISYM1)
      ISYM12  = MULD2H(ISYM1,ISYM2)

      ISYRES1 = MULD2H(ISYM11,ISY1ALBE)

      IF ( LRELAX .AND. LDERIV ) THEN
        IF ( MULD2H(ISYM12,ISY1ALBE) .NE. MULD2H(ISYM11,ISY2ALBE) ) THEN
          CALL QUIT('Symmetry mismatch in CC_IAJB. (1)')
        END IF
      END IF

      IF (LRELAX) ISYRES2 = MULD2H(ISYM12,ISY1ALBE)
      IF (LDERIV) ISYRES2 = MULD2H(ISYM11,ISY2ALBE)

      IF ( LRELAX .OR. LDERIV ) THEN
        IF ( MULD2H(ISYRES1,ISYB2) .NE. MULD2H(ISYRES2,ISYB1) ) THEN
          CALL QUIT('Symmetry mismatch in CC_IAJB. (2)')
        END IF
      END IF

      LLAO = ( LDERIV .AND. (IREAL.EQ.-1) )

      ISYDEL = ISAO(IDEL)
      ISYGAM = ISAO(IGAM)

      ISYMAI = MULD2H(ISY1ALBE,ISYM11)
      ISYMBJ = MULD2H(MULD2H(ISYDEL,ISYGAM),ISYM11)

      IF ( LTRIANG .AND. (IOPT.EQ.0) .AND. (ISYMAI.GT.ISYMBJ) ) RETURN

*---------------------------------------------------------------------*
*     work space allocation:
*---------------------------------------------------------------------*
      KX1IA   = 1
      KEND1   = KX1IA + NT1AM(ISYRES1)

      IF (LDERIV.OR.LRELAX) THEN
        KX2GDIA = KEND1
        KEND1   = KX2GDIA + NT1AM(ISYRES2)
      END IF

      IF (IOPT.GE.1 .AND. LLAO) THEN
        KX2DGIA = KEND1
        KEND1   = KX2DGIA + NT1AM(ISYRES2)
      END IF

      IF (IOPT.GE.2) THEN
         KX1IJ = KEND1
         KEND1 = KX1IJ + NMATIJ(ISYRES1)
         IF (LDERIV.OR.LRELAX) THEN
           KX2IJ = KEND1
           KEND1 = KX2IJ + NMATIJ(ISYRES2)
         END IF
      END IF

      LWRK1   = LWORK - KEND1

      IF ( LWRK1 .LT. 0) THEN
        CALL QUIT('Insufficient memory in CC_IAJB.')
      END IF

*---------------------------------------------------------------------*
*     do the first two quarter transformations using the matrices
*     XLAMDP1, XLAMDH1, XLAMDP2, and XLAMDH2
*---------------------------------------------------------------------*
      CALL CC_IAJB0(X1INT,ISY1ALBE,X2INT,ISY2ALBE,IDEL,IGAM, 
     &              WORK(KX1IA),WORK(KX2GDIA),WORK(KX2DGIA),
     &              WORK(KX1IJ),WORK(KX2IJ),
     &              XLAMDP1, XLAMDH1, ISYM1, XLAMDP2, XLAMDH2, ISYM2, 
     &              WORK(KEND1),LWRK1,IOPT,LDERIV,LRELAX,LX2ISQ,IREAL)

*---------------------------------------------------------------------*
*     Add the contribution to the result XIAJB and X2IAJB vectors
*     transform thereby gamma to j using XLMDPB1 and XLMDPB2
*---------------------------------------------------------------------*
      IF (IOPT.NE.3) THEN

         IF ( LZERO ) THEN
C           -------------------------
C           add (i a|j del) to XIAJB:
C           -------------------------
            CALL CC_IAJB1(IGAM, WORK(KX1IA), ISYRES1, ISYGAM,
     &                    XLMDPB1, ISYB1, XIAJB, LAUXG, IBASX)
         END IF
       
       
         IF ( LRELAX ) THEN
C           ------------------------------
C           add (i a|j-bar del) to X2IAJB:
C           ------------------------------
            CALL CC_IAJB1(IGAM, WORK(KX1IA), ISYRES1, ISYGAM,
     &                    XLMDPB2, ISYB2, X2IAJB, LAUXG, IBASX)
         END IF
       
         IF ( LDERIV .OR. LRELAX ) THEN
C           ------------------------------------------------
C           add (i-bar a|j del) + (i a-bar|j del) to X2IAJB:
C           ------------------------------------------------
            CALL CC_IAJB1(IGAM, WORK(KX2GDIA), ISYRES2, ISYGAM,
     &                    XLMDPB1, ISYB1, X2IAJB, LAUXG, IBASX)
         END IF

      END IF
*---------------------------------------------------------------------*
*     Add the contributions to the (ia|del j) integrals in XIABJ vector
*     and to the (ia|del j)-bar integrals in X2IABJ vector
*     transform thereby gamma to j/j-bar using XLMDHB1/XLMDHB2 
*---------------------------------------------------------------------*
      IF ( IOPT.EQ.1 .OR. IOPT.EQ.2 .OR. IOPT.EQ.3) THEN

         IF ( LZERO ) THEN
C           -------------------------
C           add (i a|del j) to XIABJ:
C           -------------------------
            CALL CC_IAJB1(IGAM, WORK(KX1IA), ISYRES1, ISYGAM,
     &                    XLMDHB1, ISYB1, XIABJ, LAUXG, IBASX)
         END IF


         IF ( LRELAX ) THEN
C           ------------------------------
C           add (i a|del j-bar) to X2IAJB:
C           ------------------------------
            CALL CC_IAJB1(IGAM, WORK(KX1IA), ISYRES1, ISYGAM,
     &                    XLMDHB2, ISYB2, X2IABJ, LAUXG, IBASX)
         END IF

         IF ( LDERIV .OR. LRELAX ) THEN
C           ------------------------------------------------
C           add (i-bar a|del j) + (i a-bar|del j) to X2IABJ:
C           ------------------------------------------------
            IF (LLAO) THEN
              CALL CC_IAJB1(IGAM, WORK(KX2DGIA), ISYRES2, ISYGAM,
     &                      XLMDHB1, ISYB1, X2IABJ, LAUXG, IBASX)
            ELSE
              CALL CC_IAJB1(IGAM, WORK(KX2GDIA), ISYRES2, ISYGAM,
     &                      XLMDHB1, ISYB1, X2IABJ, LAUXG, IBASX)
            END IF
         END IF

      END IF
*---------------------------------------------------------------------*
*     Add the contributions to the result XIJBA/X2IJBA vectors
*     transform thereby gamma to a/a-bar using XLMDHB1/XLMDHB2 
*---------------------------------------------------------------------*
      IF ( IOPT.EQ.2 .OR. IOPT.EQ.3) THEN

         IF ( LZERO ) THEN
C           -------------------------
C           add (i j|del a) to XIJBA:
C           -------------------------
            CALL CC_IJBA1(IGAM, WORK(KX1IJ), ISYRES1, ISYGAM,
     &                    XLMDHB1, ISYB1, XIJBA, LAUXG, IBASX)
         END IF


         IF ( LRELAX ) THEN
C           ------------------------------
C           add (i j|del a-bar) to X2IJBA:
C           ------------------------------
            CALL CC_IJBA1(IGAM, WORK(KX1IJ), ISYRES1, ISYGAM,
     &                    XLMDHB2, ISYB2, X2IJBA, LAUXG, IBASX)
         END IF

         IF ( LDERIV .OR. LRELAX ) THEN
C           ------------------------------------------------
C           add (i-bar j|del a) + (i j-bar|del a) to X2IJBA:
C           ------------------------------------------------
            CALL CC_IJBA1(IGAM, WORK(KX2IJ), ISYRES2, ISYGAM,
     &                    XLMDHB1, ISYB1, X2IJBA, LAUXG, IBASX)
         END IF

      END IF 

*---------------------------------------------------------------------*
*     return
*---------------------------------------------------------------------*

      RETURN
      END
*=====================================================================*
*                 END OF SUBROUTINE CC_IAJB                          *
*=====================================================================*
c /* deck cc_iabj0 */
*======================================================================*
       SUBROUTINE CC_IAJB0(XINT1,   ISY1ALBE,
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
*            iopt = 1 : in addition X2DGIA if needed 
*            iopt = 2 : in addition X1IJ and X2IJ
*            iopt = 3 : in addition X1IJ (added on 30.05.2006)
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

      LOGICAL LDERIV, LRELAX, LAO, LX2ISQ
      INTEGER IDEL, IGAM, ISY1ALBE, ISY2ALBE, ISYM1, ISYM2, LWORK, IOPT
      INTEGER IREAL
      
#if defined (SYS_CRAY)
      REAL XLAMDP1(*), XLAMDH1(*)
      REAL XLAMDP2(*), XLAMDH2(*)
      REAL XINT1(*), XINT2(*)
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
      INTEGER NRHFI, KOFF7, KOFF5, KOFF8, ISYHLF, KSCR4

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
          CALL QUIT('Symmetry mismatch in CC_IAJB0.' )
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
        CALL QUIT('Insufficient memory in CC_IAJB0.')
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
CCN   added iopt.eq.3 (30.05.2006)
      IF ( IOPT.EQ.2 .OR. IOPT.EQ.3 ) THEN

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

         ENDIF

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
c /* deck cc_iajb1 */
*=====================================================================*
      SUBROUTINE CC_IAJB1(IGAM, XIAG, ISYMAI, ISYGAM, 
     &                    XLAMDA, ISYLAM, XIABJ, LAUXG, IBASX  )
*---------------------------------------------------------------------*
*
*   Purpose: transform (ia|gam del) to (ia|j del)
*
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
      REAL XIAG(*), XIABJ(*), XLAMDA(*)
#else
      DOUBLE PRECISION XIAG(*), XIABJ(*), XLAMDA(*)
#endif

      LOGICAL LAUXG
      INTEGER IGAM, ISYMAI, ISYGAM, ISYLAM, ISYMJ, KLAMD, KOFF0,
     &        IBASX(*)

* transform integral batch:
      ISYMJ  = MULD2H(ISYGAM,ISYLAM)
      G      = IGAM - IBAS(ISYGAM)
      IF (LAUXG) G = G - IBASX(ISYGAM)

      DO J = 1, NRHF(ISYMJ)

        KLAMD = IGLMRH(ISYGAM,ISYMJ) + NBAS(ISYGAM)*(J-1) + G

        KOFF0 = IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J-1) + 1
              
        CALL DAXPY(NT1AM(ISYMAI),XLAMDA(KLAMD),XIAG,1,XIABJ(KOFF0),1)

      END DO

      RETURN
      END
*=====================================================================*
*                 END OF SUBROUTINE CC_IAJB1                          *
*=====================================================================*
c /* deck cc_ijba1 */
*=====================================================================*
      SUBROUTINE CC_IJBA1(IGAM, XIJG, ISYMIJ, ISYGAM, 
     &                    XLAMDA, ISYLAM, XIJBA, LAUXG, IBASX  )
*---------------------------------------------------------------------*
*
*   Purpose: transform (ij|gam del) to (ij|a del)
*
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
      REAL XIJG(*), XIJBA(*), XLAMDA(*)
#else
      DOUBLE PRECISION XIJG(*), XIJBA(*), XLAMDA(*)
#endif

      LOGICAL LAUXG
      INTEGER IGAM, ISYMIJ, ISYGAM, ISYLAM, ISYMA, KLAMD, KOFF0
      INTEGER ISYMAI, ISYMI, ISYMJ, KOFF5, NVIRA, IBASX(*)

* transform integral batch:
      ISYMA  = MULD2H(ISYGAM,ISYLAM)
      G      = IGAM - IBAS(ISYGAM)
      IF (LAUXG) G = G - IBASX(ISYGAM)
      NVIRA  = MAX(NVIR(ISYMA),1)

      DO A = 1, NVIR(ISYMA)

        KLAMD = IGLMVI(ISYGAM,ISYMA) + NBAS(ISYGAM)*(A-1) + G

        DO ISYMJ = 1, NSYM
           ISYMI  = MULD2H(ISYMIJ,ISYMJ)
           ISYMAI = MULD2H(ISYMA,ISYMI)

           DO J = 1, NRHF(ISYMJ)

              KOFF5 = IMATIJ(ISYMI,ISYMJ)  + NRHF(ISYMI)*(J-1)   + 1
              KOFF0 = IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J-1) 
     &              + IT1AM(ISYMA,ISYMI) + A
              
              CALL DAXPY(NRHF(ISYMI),XLAMDA(KLAMD),XIJG(KOFF5),1,
     &                                             XIJBA(KOFF0),NVIRA)
           END DO
        END DO

      END DO

      RETURN
      END
*=====================================================================*
*                 END OF SUBROUTINE CC_IJBA1                          *
*=====================================================================*
