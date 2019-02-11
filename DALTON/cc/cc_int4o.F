!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
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
#  include "implicit.h"
#endif
#include "ccorb.h"
#include "maxorb.h"
#include "ccsdsym.h"

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
