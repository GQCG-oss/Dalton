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
      SUBROUTINE CCRHS_IPM1(XINT,XINTP,XINTM,SCRAB,INDV1,INDV2,
     *                     ISYMAB,ISYMG,NUMG,IG1,IG2,IOPT)
C
C     Written by Henrik Koch 17-aug-1994.
C
C     Purpose: Making plus and minus combination of integrals.
C              (a>=g|bd) -> K+ and K- where
C                           K+- = (ag|bd) +- (bg|ad) a<=b,g<=d
C
C     Modified by Sonia Coriani 26-oct-1999 to 
C              handle XINT with squared (ag|bd) (ag part): 
C              If IOPT = 0, XINT in input is packed  (a>=g|bd)
C              If IOPT = 1, XINT in input is squared (ag|bd)
C        
C
#include <implicit.h>
#include <maxorb.h>
      PARAMETER(ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XINT(*),XINTP(*),XINTM(*),SCRAB(*)
      DIMENSION INDV1(*), INDV2(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <symsq.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ISYDIS = MULD2H(ISYMAB,ISYMG)
C
C
      DO 100 G = IG1,IG2
C
         IG = G - IG1 + 1
C
         DO 110 ISYMB = 1,NSYM
C
            ISYMA  = MULD2H(ISYMB,ISYMAB)
            ISYMAG = MULD2H(ISYMA,ISYMG)
C
            NTOTA  = MAX(NBAS(ISYMA),1)
cs
            IF (IOPT.EQ.0) THEN
               NTOTAG = MAX(NNBST(ISYMAG),1)
            ELSE IF (IOPT.EQ.1) THEN
               NTOTAG = MAX(N2BST(ISYMAG),1)
            ELSE
               CALL QUIT('Unknown option in CCRHS_IPM1')
            END IF
C
            DO 120 A = 1,NBAS(ISYMA)
C
               IF (IOPT.EQ.0) THEN
                  IF (ISYMA .EQ. ISYMG) THEN
                     KOFF1 = IDSAOG(ISYMB,ISYDIS) + IAODPK(ISYMA,ISYMG)
     *                     + INDEX(G,A)
                  ELSE IF (ISYMA .LT. ISYMG) THEN
                     KOFF1 = IDSAOG(ISYMB,ISYDIS) + IAODPK(ISYMA,ISYMG)
     *                     + NBAS(ISYMA)*(G - 1) + A
                  ELSE
                     KOFF1 = IDSAOG(ISYMB,ISYDIS) + IAODPK(ISYMA,ISYMG)
     *                     + NBAS(ISYMG)*(A - 1) + G
                  ENDIF
               ELSE IF (IOPT.EQ.1) THEN
                  KOFF1 = IDSAOGSQ(ISYMB,ISYDIS) + IAODIS(ISYMA,ISYMG)
     &                  + NBAS(ISYMA)*(G - 1) + A
               END IF
C
               KOFF2 = IAODIS(ISYMA,ISYMB) + A
C
               CALL DCOPY(NBAS(ISYMB),XINT(KOFF1),NTOTAG,
     *                    SCRAB(KOFF2),NTOTA)
C
  120       CONTINUE
C
  110    CONTINUE
C
         KOFF = NNBST(ISYMAB)*(IG - 1)
C
#if !defined (SYS_CRAY)
         DO 130 I = 1,NNBST(ISYMAB)
C
            XINTP(KOFF + I) = SCRAB(INDV1(I))
            XINTM(KOFF + I) = SCRAB(INDV2(I))
C
  130    CONTINUE
#endif
#if defined (SYS_CRAY)
         CALL GATHER(NNBST(ISYMAB),XINTP(KOFF + 1),SCRAB,INDV1)
         CALL GATHER(NNBST(ISYMAB),XINTM(KOFF + 1),SCRAB,INDV2)
#endif
C
  100 CONTINUE
C
C
      NTOT = NNBST(ISYMAB)*NUMG
C
      CALL DAXPY(NTOT,ONE,XINTM,1,XINTP,1)
      CALL DSCAL(NTOT,-TWO,XINTM,1)
      CALL DAXPY(NTOT,ONE,XINTP,1,XINTM,1)
C
      RETURN
      END
