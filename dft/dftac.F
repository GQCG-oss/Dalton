C
C...   Copyright (c) 1997 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dalton, an ab initio electronic structure program, Release 1.0
C...   (1997), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, T. Andersen, K. L. Bak, V. Bakken,
C...   O. Christiansen, P. Dahle, E. K. Dalskov, T. Enevoldsen,
C...   H. Heiberg, D. Jonsson, S. Kirpekar, R. Kobayashi, H. Koch,
C...   K. V. Mikkelsen, P. Norman, M. J. Packer, T.Saue,
C...   P. R. Taylor, and O. Vahtras"
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
#include <single.h>
C /* Deck dftac */ 
      SUBROUTINE DFTAC(VA,DST,CRX,CRY,CRZ)
C
C     T. Helgaker
C
#include <implicit.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0) 
      LOGICAL INSIDE 
      DIMENSION DST(NATOMS)
#include <dfterg.h>
#include <dftcom.h>
#include <dftacb.h>
#include <nuclei.h>
#include <codata.h>
C
      DIMENSION RADIUS(0:18)
      DATA RADIUS /0.0D0,0.35D0,0.35D0,
     &      1.45D0,1.05D0,0.85D0,0.70D0,0.65D0,0.60D0,0.50D0,0.45D0,
     &      1.80D0,1.50D0,1.25D0,1.10D0,1.00D0,1.00D0,1.00D0,1.00D0/
C
      SHIFT = DFTIPT + DFTHMO
C
      INSIDE = .FALSE.
      DO 100 I = 1, NATOMS
         BRAG   = DFTBR1*RADIUS(NINT(CHARGE(I)))/XTANG
         DST(I) = DSQRT((CRX - CORD(1,I))**2
     &                + (CRY - CORD(2,I))**2
     &                + (CRZ - CORD(3,I))**2)
         IF (DST(I).LT.BRAG) INSIDE = .TRUE.
  100 CONTINUE
C
      IF (.NOT.INSIDE) THEN
         RHOT = D0
         RINV = D0 
         DO 200 I = 1, NATOMS
            RHOT = RHOT + DFTCHR(I)
            RINV = RINV + DFTCHR(I)/DST(I)
  200    CONTINUE
         ASYMP = -WDFTX*RINV/RHOT + SHIFT
C
         NFRC = 0
         FRAC = D0 
         DO 300 I = 1, NATOMS
            BRAG  = RADIUS(NINT(CHARGE(I)))/XTANG
            BRAG1 = DFTBR1*BRAG
            BRAG2 = DFTBR2*BRAG
            IF (DST(I) .LT. BRAG2) THEN
               NFRC = NFRC + 1
               FRAC = FRAC + (DST(I) - BRAG1)/(BRAG2 - BRAG1)
            END IF
  300    CONTINUE
C
         IF (NFRC.EQ.0) THEN
            VA = ASYMP
         ELSE
            VA = VA + FRAC*(ASYMP - VA)/NFRC
         END IF
      END IF
C
      VA = VA - SHIFT
      RETURN
      END
