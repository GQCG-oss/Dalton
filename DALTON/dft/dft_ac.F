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
C /* Deck dftac */ 
      SUBROUTINE DFTAC(VA,VFA,DST,CRX,CRY,CRZ,GRD,RHO43,RHO13)
C
C     T. Helgaker
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
      PARAMETER (D0 = 0.0D0) 
      LOGICAL INSIDE 
      DIMENSION DST(NATOMS)
#include "inforb.h"
#include "dfterg.h"
#include "dftcom.h"
#include "dftacb.h"
#include "codata.h"
#include "gnrinf.h"
C
      DIMENSION RADIUS(0:54)

      DATA RADIUS/0.00D0,0.35D0,0.35D0,
     &1.45D0,1.05D0,0.85D0,0.70D0,0.65D0,0.60D0,0.50D0,0.45D0,
     &1.80D0,1.50D0,1.25D0,1.10D0,1.00D0,1.00D0,1.00D0,1.00D0,
     &2.20D0,1.80D0,
     &1.60D0,1.40D0,1.35D0,1.40D0,1.40D0,
     &1.40D0,1.35D0,1.35D0,1.35D0,1.35D0,
     &1.30D0,1.25D0,1.15D0,1.15D0,1.15D0,1.15D0,
     &18*1.3D0/


      IF (NASHT .GT. 0 ) THEN
         WRITE (LUPRI,*) 
     &    'Fatal error: AC not implemented for 
     &     restricted open-shell DFT'
         CALL QUIT(
     &    'Fatal error: is not implemented for 
     &     restricted open-shell DFT')
      ENDIF

      SHIFT = DFTIPTA + EHOMO  

      IF (DOLB94) THEN
        BETA_LB = 0.05d0
        XGRD = GRD / RHO43
        XASINH = dlog(XGRD+DSQRT(1.0d0+XGRD*XGRD))
        VLDA = 2.0d0*((2.0d0**(1.0d0/3.0d0))
     &         *((3.0d0/PI)**(1.0d0/3.0d0))
     &         *RHO13) 
        VFA = (BETA_LB*RHO13*(XGRD*XGRD)
     &         /(1.0d0+(3.0d0*BETA_LB*XGRD*XASINH)))+VLDA
      ENDIF

      IF (LGRAC) THEN
         XGRD = GRD / RHO43
         BETA = DFTBR2 
         ALPHA = DFTBR1 
         ASYMP = (HFXFAC-1.0d0)*VFA + SHIFT
         FRAC = 1.0d0 / (1.0d0 + exp(-ALPHA*(XGRD-BETA))) 
         VA = (1.0d0-FRAC)*VA + FRAC*ASYMP
      ELSE
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
!If the point is defined as part of the core potential -- shift down
      IF (.NOT.INSIDE) THEN

        IF (DOMPOLE) THEN
          RINV=0.0d0
          XNELEC=0.0d0
          DO I = 1, NATOMS
            RINV = RINV + (CHARGE(I)+(KCHARG/NATOMS))/DST(I) 
            XNELEC = XNELEC + (CHARGE(I)+(KCHARG/NATOMS))
          ENDDO
          VFA = RINV/XNELEC
        ENDIF

         ASYMP = (HFXFAC-1.0d0)*VFA + SHIFT

         NFRC = 0
         FRAC = D0 
         IF (LTAN) FRAC = 1.0d0
         DO 300 I = 1, NATOMS
            BRAG  = RADIUS(NINT(CHARGE(I)))/XTANG
            BRAG1 = DFTBR1*BRAG
            BRAG2 = DFTBR2*BRAG
            IF (DST(I) .LT. BRAG2) THEN
               IF(DST(I).LT.DISTMIN) DISTMIN = DST(I)
               IF(DST(I).GT.DISTMAX) DISTMAX = DST(I)
               IF (LLIN) THEN
                NFRC = NFRC + 1
                FRAC = FRAC + (DST(I) - BRAG1)/(BRAG2 - BRAG1)
               ELSE IF (LTAN) THEN
                NFRC = NFRC + 1
                WA = (1.0d0 / ((DFTBR2-DFTBR1)*BRAG))
     &               *LOG((1.0d0+0.998d0)/(1.0d0-0.998d0))
                CA = ((DFTBR1+DFTBR2)*BRAG)/2.0d0
                TANHA = 0.5d0*(tanh(WA*(DST(I)-CA))+1.0d0)  
                FRAC = FRAC * TANHA
               ENDIF
            END IF
  300    CONTINUE
         IF (NFRC.EQ.0) THEN
            VA = ASYMP
         ELSE IF (LLIN) THEN
             VA = VA + FRAC*(ASYMP - VA)/NFRC
         ELSE IF (LTAN) THEN
             VA = VA + FRAC*(ASYMP - VA)
         END IF
      END IF
      END IF !LGRAC ENDIF
      RETURN
      END

