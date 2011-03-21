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
C /* Deck dftac */ 
      SUBROUTINE DFTAC(VA,DST,CRX,CRY,CRZ)
C
C     T. Helgaker
C
#include <implicit.h>
#include <mxcent.h>
#include <nuclei.h>
      PARAMETER (D0 = 0.0D0) 
      LOGICAL INSIDE 
      DIMENSION DST(NATOMS)
#include <dfterg.h>
#include <dftcom.h>
#include <dftacb.h>
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
