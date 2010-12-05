C
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren.
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
#include <single.h>
C
*=====================================================================*
C  /* Deck cc_int4o */
      SUBROUTINE CC_INT4O(XIJK0,ISYIJK0,XIJK1,ISYIJK1,
     &                    XLAMDA0,ISYLAM0,
     &                    XLAMDA1,ISYLAM1, ISYMD,
     &                    XIJKL,LRELAX,WORK,LWORK,IOPT)
*---------------------------------------------------------------------*
*
*     Purpose: transform the del index of (jk|ldel) to occupied L.
*
*     IOPT = 1 --> XIJKdel0 * XLAMDA0                (LRELAX = .FALSE.)
*     IOPT = 2 --> XIJKdel0 * XLAMDA1 + XIJKdel1 * XLAMDA0
*     XIJKL assumed initialized OUTSIDE
*
*     Sonia Coriani, 10/09-1999
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <ccorb.h>
#include <maxorb.h>
#include <ccsdsym.h>

      INTEGER ISYIJK0,ISYIJK1,ISYLAM0,ISYLAM1,IOPT, ISYMD,LWORK
      LOGICAL LRELAX

#if defined (SYS_CRAY)
      REAL XIJK0(*), XIJK1(*), XIJKL(*)
      REAL XLAMDA0(*), XLAMDA1(*), WORK(LWORK)
      REAL ZERO, ONE, HALF, DDOT, XNORM
#else
      DOUBLE PRECISION XIJK0(*), XIJK1(*), XIJKL(*)
      DOUBLE PRECISION XLAMDA0(*), XLAMDA1(*), WORK(LWORK)
      DOUBLE PRECISION ZERO, ONE, HALF, DDOT, XNORM
#endif
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0)
      INTEGER KOFF1, KOFF2, KOFF3
      INTEGER ISYIJKL, ISYML, NBASD, NTOIJK
      INTEGER ISYIJKDE
*
*     --------------------------------------
*     Begin: symmetry of result IJKL integrals
*     --------------------------------------
*
      IF (IOPT.EQ.1) THEN
         ISYIJKDE = MULD2H(ISYIJK0,ISYMD)
         ISYIJKL  = MULD2H(ISYIJKDE,ISYLAM0)
      ELSE
         ISYIJKDE = MULD2H(ISYIJK1,ISYMD)
         ISYIJKL  = MULD2H(ISYIJKDE,ISYLAM0)
         IF (ISYIJKL .NE. MULD2H(ISYLAM1,MULD2H(ISYIJK0,ISYMD)))
     &       CALL QUIT('Symmetry mismatch in CC_INT4O' )
      END IF
*
*     -------------------------------------------------------*
*        Transform AO integral index to occupied space. (L)
*     -------------------------------------------------------*
*
      IF (IOPT.EQ.1) THEN               
         ISYML = MULD2H(ISYLAM0,ISYMD)
         KOFF1 = 1 
         KOFF2 = IGLMRH(ISYMD,ISYML) + 1
         KOFF3 = I3ORHF(ISYIJK0,ISYMD) + 1 
         
         NBASD   = MAX(NBAS(ISYMD),1)
         NTOIJK  = MAX(NMAIJK(ISYIJK0),1)

         CALL DGEMM('N','N',NMAIJK(ISYIJK0),NRHF(ISYML),NBAS(ISYMD),
     &              ONE,XIJK0(KOFF1),NTOIJK,XLAMDA0(KOFF2),NBASD,
     &              ONE,XIJKL(KOFF3),NTOIJK)
       ELSE 

         ISYML  = MULD2H(ISYLAM0,ISYMD)
         KOFF1  = 1
         KOFF2  = IGLMRH(ISYMD,ISYML) + 1
         KOFF3  = I3ORHF(ISYIJK1,ISYML) + 1
         NTOIJK = MAX(NMAIJK(ISYIJK1),1)
         NBASD  = MAX(NBAS(ISYMD),1)
         
         CALL DGEMM('N','N',NMAIJK(ISYIJK1),NRHF(ISYML),NBAS(ISYMD),
     &              ONE,XIJK1(KOFF1),NTOIJK,XLAMDA0(KOFF2),NBASD,
     &              ONE,XIJKL(KOFF3),NTOIJK)
c
         ISYML  = MULD2H(ISYLAM1,ISYMD)
         KOFF1  = 1
         KOFF2  = IGLMRH(ISYMD,ISYML) + 1
         KOFF3  = I3ORHF(ISYIJK0,ISYML) + 1
         NTOIJK = MAX(NMAIJK(ISYIJK0),1)
         NBASD  = MAX(NBAS(ISYMD),1)

         CALL DGEMM('N','N',NMAIJK(ISYIJK0),NRHF(ISYML),NBAS(ISYMD),
     &              ONE,XIJK0(KOFF1),NTOIJK,XLAMDA1(KOFF2),NBASD,
     &              ONE,XIJKL(KOFF3),NTOIJK)

      END IF
      RETURN
      END
*=====================================================================*
