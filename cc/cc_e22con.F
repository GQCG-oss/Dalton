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
#  include <implicit.h>
#endif
#include <ccorb.h>
#include <maxorb.h>
#include <ccsdsym.h>

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
