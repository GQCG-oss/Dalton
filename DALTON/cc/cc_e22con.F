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
C  /* Deck cc_e22con */
      SUBROUTINE CC_E22CON(CTR2,ISYCTR,TAMP1,ISYTAM,QGAMMA,ISYGAM,
     &                     RHO1,ISYRHO,WORK,LWORK)
                           
*---------------------------------------------------------------------*
*
*     Purpose:  lead calculation of E2' contribution to FBTA
*               transformed vector (second part)
*
*     Sonia Coriani, 14/09-1999
*
*  Transform ZA2_bj,ai to ZA2_kj,ai with TA1_bk
*  Resort ZA2_kj,ai to ZA2_jki,a
*  Compute final result sum_jki ZA2_jki,a * Gamma_jki,i
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "ccorb.h"
#include "maxorb.h"
#include "ccsdsym.h"

      INTEGER ISYCTR,ISYTAM,ISYGAM,ISYRHO,LWORK
      LOGICAL LRELAX
  
#if defined (SYS_CRAY)
      REAL CTR2(*), TAMP1(*), QGAMMA(*), RHO1(*), WORK(LWORK)
      REAL ZERO, ONE, HALF, DDOT, XNORM
#else
      DOUBLE PRECISION CTR2(*),TAMP1(*),QGAMMA(*),RHO1(*),WORK(LWORK)
      DOUBLE PRECISION ZERO, ONE, HALF, DDOT, XNORM
#endif
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0)
*
      INTEGER ISYMI, ISYJKM, ISYMA,ISYRES,ISYMAI,ISYZTA
      INTEGER KOFFZ,KOFFG,KOFFR,KEND1,KZKJAM,KZJKMA
      INTEGER NVIRA,NTOJKM,LWRK1,IOPT
*
* Symmetry checks
*
      ISYZTA = MULD2H(ISYCTR,ISYTAM)
      ISYRES = MULD2H(ISYZTA,ISYGAM)
      IF (ISYRES.NE.ISYRHO) CALL QUIT('Symmetry mismatch in E2 2nd')
*
* allocate room for transformed Zeta's
*
      KZKJAM = 1
      KZJKMA = KZKJAM + N3OVIR(ISYZTA)
      KEND1  = KZJKMA + N3OVIR(ISYZTA)
      LWRK1  = LWORK - KEND1
*
* transform ZA2_bj,am to ZA2_kj,am with TA1_bk
*
      CALL CC_ZKJAM(CTR2,ISYCTR,TAMP1,ISYTAM,WORK(KZKJAM))
*
* resort to ZA2_jkm,a
*
      IOPT = 1
      CALL CC_SORTZ2(WORK(KZKJAM),WORK(KZJKMA),ISYZTA,IOPT)
*
* contract  sum_{jkm} ZA2_jkm,a GammaQ_jkm,i = rho_ai
*
      DO ISYMI = 1, NSYM
         ISYJKM = MULD2H(ISYGAM,ISYMI)
         ISYMA  = MULD2H(ISYJKM,ISYZTA)
* check
         ISYMAI = MULD2H(ISYMA,ISYMI)
         IF (ISYMAI.NE.ISYRHO) 
     *        CALL QUIT('Symmetry mismatch 2 in E2 2nd')

         KOFFZ  = I3OVIR(ISYJKM,ISYMA) + KZJKMA
         KOFFG  = I3ORHF(ISYJKM,ISYMI) + 1
         KOFFR  = IT1AM(ISYMA,ISYMI)   + 1

         NVIRA   = MAX(NVIR(ISYMA),1)
         NTOJKM  = MAX(NMAIJK(ISYJKM),1)

         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMAIJK(ISYJKM),
     &              ONE,WORK(KOFFZ),NTOJKM,QGAMMA(KOFFG),NTOJKM,
     &              ONE,RHO1(KOFFR),NVIRA) 

      END DO

      RETURN
      END
*=====================================================================*
