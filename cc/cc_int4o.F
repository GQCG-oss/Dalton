C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
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
#  include <ccorb_dec.h>
#  include <maxorb_dec.h>
#  include <ccsdsym_dec.h>
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
