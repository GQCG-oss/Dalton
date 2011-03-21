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
C  /* Deck cclt_yps1 */
      SUBROUTINE CCLT_YPS1(CTR1,ISYCTR,YI,ISYMYI,XLAMDH1,
     &                     ISYLM1,XLAMDH2,ISYLM2,YPS)
C
C     Purpose: To calculate the Ypsilon-type intermediates:
C
C     Yps(alpha a)  =   sum_k XLAMDH(alpha k) CTR1(a k) 
C                     + sum_f XLAMDH(alpha f)   YI(f a) 
C
C     ISYCTR : symmetry of CTR1, YI               (Zeta_1)
C     ISYLAM : symmetry of XLAMDH
C
C     Christof Haettig, October 1998
C    
C     Generalized for FbTa transformation:
C
C     YpsA(alpha a)  =   sum_k XLAMDH1(alpha k) CTR1(a k)
C                      + sum_f XLAMDH2(alpha f)   YI(f a)
C     All vectors and matrices can have general symmetry
C     but the two contributions to YpsA must in total have
C     the same symmetry. 
C     If ISYYPS is given in input, ISYLM1 and ISYLM2 not needed
C
C     Sonia Coriani, February 1999
C
#include <implicit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION CTR1(*),XLAMDH1(*),XLAMDH2(*),YPS(*),YI(*)
C
C---------------------------------------------------
C     Half-transformation to AO-basis of CTR1 and YI
C---------------------------------------------------
C
      ISYMAO1 = MULD2H(ISYCTR,ISYLM1)
      ISYMAO2 = MULD2H(ISYMYI,ISYLM2)
      IF (ISYMAO1.NE.ISYMAO2) CALL QUIT('Symmetry mismatch in CCLT_YPS')
      ISYYPS  = ISYMAO1
C
      CALL DZERO(YPS,NGLMDT(ISYYPS))
C
      DO ISYMAL = 1,NSYM                   !alpha
C
         ISYMA = MULD2H(ISYMAL,ISYYPS)      
         ISYMK = MULD2H(ISYMA,ISYCTR)
         ISYMF = MULD2H(ISYMA,ISYMYI)
C
         KOFF1 = IGLMRH(ISYMAL,ISYMK) + 1    !offset LambdaH1_al,k
         KOFF2 = IT1AM(ISYMA,ISYMK)   + 1    !offset Zeta1_ak
         KOFF3 = IGLMVI(ISYMAL,ISYMA) + 1    !offset Yps_al,a
C
         NTOTBA = MAX(NBAS(ISYMAL),1)
         NTOTVI = MAX(NVIR(ISYMA),1)
C
C backtransformation of Zeta1 --> Zeta_al,a
C
         CALL DGEMM('N','T',NBAS(ISYMAL),NVIR(ISYMA),NRHF(ISYMK),
     *               ONE,XLAMDH1(KOFF1),NTOTBA,CTR1(KOFF2),NTOTVI,
     *               ONE,YPS(KOFF3),NTOTBA)
C
         KOFF4 = IMATAB(ISYMF,ISYMA)  + 1    !offset YI_fa
         KOFF5 = IGLMVI(ISYMAL,ISYMF) + 1    !offset LambdaH2_al,f 
C
         NTOTVI = MAX(NVIR(ISYMF),1)
C
C backtranformation of YI  --> YI_al,a
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NVIR(ISYMA),NVIR(ISYMF),
     *               ONE,XLAMDH2(KOFF5),NTOTBA,YI(KOFF4),NTOTVI,
     *               ONE,YPS(KOFF3),NTOTBA)
C
      END DO
C
      RETURN
      END
