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
#include <single.h>
C
C /* Deck dftac */ 
      SUBROUTINE DFTAC(VA,CRX,CRY,CRZ,WGHT)
C
C     T. Helgaker
C
#include <implicit.h>
#include <mxcent.h>
#include <nuclei.h>
      PARAMETER (D0 = 0.0D0) 
      LOGICAL INSIDE 
      DIMENSION DST(MXCENT)
#include <dfterg.h>
#include <dftcom.h>
#include <dftacb.h>
#include <codata.h>
C
      DIMENSION RADIUS(0:18)
      DATA RADIUS/ 0.75D0,
     &  0.35D0,0.35D0,
     &  1.45D0,1.05D0,0.85D0,0.70D0,0.65D0,0.60D0,0.50D0,0.45D0,
     &  1.80D0,1.50D0,1.25D0,1.10D0,1.00D0,1.00D0,1.00D0,1.00D0/
C
C     Skip at first iteration
C
      IF (EHOMO.EQ.D0) RETURN
C
C     Test if any correction should be done
C
      INSIDE = .FALSE.
      DO I = 1, NATOMS
         BRAG   = DFTBR1*RADIUS(NINT(CHARGE(I)))/XTANG
         DST(I) = DSQRT((CRX - CORD(1,I))**2
     &                + (CRY - CORD(2,I))**2
     &                + (CRZ - CORD(3,I))**2)
         IF (DST(I).LT.BRAG) INSIDE = .TRUE.
      END DO
C
C     Apply correction if necessary
C
      IF (.NOT.INSIDE) THEN
C
C        Calculate correction
C
         RHOT = D0
         RINV = D0 
         DO I = 1, NATOMS
            RHOT = RHOT + DFTCH0(I)
            RINV = RINV + DFTCH0(I)/DST(I)
         END DO
         ASYMP = WGHT*(-RINV/RHOT + DFTIPT + EHOMO)
C
         NFRC = 0
         FRAC = D0 
         DO I = 1, NATOMS
            BRAG  = RADIUS(NINT(CHARGE(I)))/XTANG
            BRAG1 = DFTBR1*BRAG
            BRAG2 = DFTBR2*BRAG
            IF (DST(I) .LT. BRAG2) THEN
               NFRC = NFRC + 1
               FRAC = FRAC + (DST(I) - BRAG1)/(BRAG2 - BRAG1)
            END IF
         END DO
C
C        Add correction
C
         IF (NFRC.EQ.0) THEN
            VA = ASYMP
         ELSE
            VA = VA + FRAC*(ASYMP - VA)/NFRC
         END IF
C
      END IF
      RETURN
      END
