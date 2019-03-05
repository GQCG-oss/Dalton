!
!...   Copyright (c) 2014 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2014 (2014), see http://daltonprogram.org"
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


C     **************************
C     ***** Dirac exchange *****
C     **************************

C*****************************************************************************
      PURE SUBROUTINE EDRC(ED,RHO,RHO13)
C*****************************************************************************
      implicit none
#include "pi.h"
      real*8, intent(in) :: RHO, RHO13
      real*8, intent(out) :: ED
      real*8, parameter :: D1 = 1.0d0 
      real*8, parameter :: D3 = 3.0d0 
      real*8, parameter :: D4 = 4.0d0 
      ED = -(D3/D4)*((D3/PI)**(D1/D3))*RHO*RHO13
      RETURN
      END


C*****************************************************************************
      PURE SUBROUTINE VDRC(VD,RHO13)
C*****************************************************************************
      implicit none
#include "pi.h"
      real*8, intent(in) :: RHO13
      real*8, intent(out) :: VD
      real*8, parameter :: D1 = 1.0d0 
      real*8, parameter :: D3 = 3.0d0 
      VD = -((D3/PI)**(D1/D3))*RHO13
      RETURN
      END


C*****************************************************************************
      PURE SUBROUTINE V1DRC(VD0,VD1,RHO,RHO13)
C*****************************************************************************
      implicit none
#include "pi.h"
      real*8, intent(in) :: RHO, RHO13
      real*8, intent(out) :: VD0, VD1
      real*8, parameter :: D1 = 1.0d0 
      real*8, parameter :: D2 = 2.0d0 
      real*8, parameter :: D3 = 3.0d0 
      real*8, parameter :: D13 = D1/D3
      VD0 = (-(D3/PI)**D13)*RHO13
      VD1 = (-(D3/PI)**D13)*D13*RHO13/RHO
      RETURN
      END


C     ********************************
C     ***** Vosko, Wilk & Nusair ***** 
C     ********************************

C*****************************************************************************
      SUBROUTINE EVWN(ECV,RHO,RHO13)
C*****************************************************************************
C
C     Correlation energy (G. J. Laming)
C
C*****************************************************************************
#include "implicit.h"
      CALL WVWN(VA1,RHO,RHO13,E1,.FALSE.,.TRUE.)
      ECV = E1*RHO
      RETURN
      END


C*****************************************************************************
      SUBROUTINE VVWN(VA,RHO,RHO13)
C*****************************************************************************
#include "implicit.h"
      CALL WVWN(VA,RHO,RHO13,E1,.TRUE.,.TRUE.)
      RETURN
      END


C*****************************************************************************
      SUBROUTINE WVWN(VA1,RHO,RHO13,E1,POTENT,SPNCMP)
C*****************************************************************************
C
C     C. W. MURRAY & G. J. LAMING
C
C*****************************************************************************
#include "implicit.h" 
#include "pi.h"
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, D4 = 4.0D0,
     &           D6 = 6.0D0)
      PARAMETER (CV = 0.5D0) 
c      PARAMETER (DCRS = (D3/(D4*PI))**(D1/D6))
C    &           AI   = 0.0621814D0,
C    &           BI   = 13.0720D0,
C    &           CI   = 42.7198D0,
C    &           X0I  = -0.409286D0,
C    &           AI=0.0621814D0,
C    &           BI=3.72744D0,
C    &           CI=12.9352D0,
C    &           X0I=-0.10498D0,
C    &           Q    = sqrt(D4*CI - BI**2),
C    &           XF0I = X0I*X0I + BI*X0I + CI,
C    &           YF0I = Q/(BI + 2*X0I))
      LOGICAL POTENT,SPNCMP
      DCRS = (D3/(D4*PI))**(D1/D6)
C
C     See Table 5 in VWN paper
C
      IF (SPNCMP) THEN
         AI  =  0.0621814D0
         X0I = -0.10498D0
         BI  =  3.72744D0
         CI  = 12.9352D0
      ELSE
         AI  =  0.0310907D0
         X0I = -0.32500D0
         BI  =  7.06042D0 
         CI  = 18.0578D0
      END IF
      Q    = sqrt(D4*CI - BI**2)
      XF0I = X0I*X0I + BI*X0I + CI
      YF0I = Q/(BI + 2*X0I)
C
      X   = DCRS/sqrt(RHO13)
      B   = X0I/XF0I
      C   = XF0I*YF0I
      XF  = X*X + BI*X + CI
      XFX = D2*X + BI
      YF  = Q/XFX
      YFX = -D2*YF/XFX
      E1  = AI*CV*(D2*log(X) + log(XF)*(BI*B - D1)
     &            +D2*BI*(D1/Q - X0I/C)*DATAN(YF)
     &            -D2*(BI*B)*log(X - X0I))
C
      IF (POTENT) THEN
         XRHO = -X/(RHO*D6)
         EX1 = AI*CV*(D2/X + XFX*(BI*B - D1)/XF
     &               +D2*BI*(D1/Q-X0I/C)*YFX/(D1+YF*YF)
     &               -D2*(BI*B)/(X-X0I))
         VA1 = E1 + RHO*EX1*XRHO
      ENDIF
      RETURN
      END


C*****************************************************************************
      SUBROUTINE V1VWN(VD0,VD1,RHO,RHO13)
C*****************************************************************************
C
C     WRITTEN   By    G. J. Laming
C     Some adjustments by A. M. Lee
C
C*****************************************************************************
#include "implicit.h"
#include "pi.h"
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, D4 = 4.0D0,
     &           D6 = 6.0D0, D7 = 7.0D0, D8 = 8.0D0)
      PARAMETER (CV = D1/D2) 
c      PARAMETER (DCRS = (D3/(D4*PI))**(D1/D6),
C    &           AI   = 0.0621814D0,
C    &           BI   = 13.0720D0,
C    &           CI   = 42.7198D0,
C    &           X0I  = -0.409286D0,
      PARAMETER (AI=0.0621814D0,
     &           BI=3.72744D0,
     &           CI=12.9352D0,
     &           X0I=-0.10498D0)

      Q    = sqrt(D4*CI - BI**2)
      XF0I = X0I*X0I + BI*X0I + CI
      YF0I = Q/(BI + 2*X0I)
      DCRS = (D3/(D4*PI))**(D1/D6)
C
      X    = DCRS/sqrt(RHO13)
      B    = X0I/XF0I
      C    = XF0I*YF0I
      XF   = X*X + BI*X+CI
      XFX  = D2*X + BI
      YF   = Q/XFX
      ACON = B*BI - D1
      BCON = D2*ACON + D2
      CCON = D2*BI*(D1/Q - X0I/C)
C
      E1   = D2*log(X) 
     &     + ACON*log(XF)
     &     - BCON*log(X - X0I)
     &     + CCON*DATAN(YF)
      EX1  = D2/X 
     &     + ACON*XFX/XF
     &     - BCON/(X - X0I)
     &     - CCON*(D2*YF/XFX)/(D1 + YF*YF)
      EXX1 = D2/(X*X) 
     &     - ACON*(D2/XF - XFX*XFX/(XF*XF))
     &     - BCON/((X - X0I)**2)
     &     - CCON*D8*Q*XFX/((Q*Q + XFX*XFX)**2)
C
      XRHO = -X/(D6*RHO)
      VD0  = AI*CV*(E1 + RHO*EX1*XRHO)
      VD1  = AI*CV*XRHO*(D2*EX1 - RHO*(EXX1*XRHO + D7*EX1/(D6*RHO)))
C
      RETURN
      END


C*****************************************************************************
      SUBROUTINE VTVWN(VD1,RHO,RHO13)
C*****************************************************************************
#include "implicit.h"
      PARAMETER(D1 = 1.0D0, D2 = 2.0D0, POW13 = 1.0D0/3.0D0, 
     &          D8D9 = 8.0D0/9.0D0)
      FAC = -D8D9/(D2**POW13 - D1)
      CALL WVWN(VA,RHO,RHO13,E0,.FALSE.,.TRUE.)
      CALL WVWN(VA,RHO,RHO13,E1,.FALSE.,.FALSE.)
      VD1 = FAC*(E0 - E1)/RHO
      RETURN
      END


C*****************************************************************************
      SUBROUTINE EBCK(BECKI,RHO,RHO13,RHOGRD)
C*****************************************************************************
C     Evaluates only the energy
C     *************WRITTEN BY G. J. LAMING********************
C*****************************************************************************
#include "implicit.h"
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0)
      PARAMETER (BETAP = 0.0042D0)
      DF = (D1/D2)**(D1/D3)
C
      XALPHA = RHOGRD/(DF*RHO*RHO13)
      SINHIA = log(XALPHA+sqrt((XALPHA*XALPHA) + D1))
      DENOMA = (D1 + (6.0D0*BETAP*XALPHA*SINHIA))
      BECKI  = -BETAP*((DF*RHO*RHO13)*XALPHA*XALPHA/DENOMA)
      RETURN
      END


C*****************************************************************************
      SUBROUTINE GBCK(VB1,VB2,RHO,RHO13,RHOGRD)
C*****************************************************************************
C     Evaluates the first derivatives VB1, VB2 for KS scheme
C     VB1 = dE/dRHO
C     VB2 = dE/dRHOGRDa2 = 4*dE/dRHOGRD2
C     where RHOGRD   : gradient of density
C           RHOGRD2  : square of gradient of density
C           RHOGRDa  : gradient of spin-alpha density
C           RHOGRDa2 : square of gradient of spin-alpha density
C     *************WRITTEN BY G. J. LAMING********************
C*****************************************************************************
#include "implicit.h"
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, DP5 = 0.5D0, D6 = 6.0D0)
      PARAMETER (FT = 4.0D0/3.0D0, BETAP = 0.0042D0)
      DF = D2**FT
C
C     AMMENDED FEBRUARY 1994 G. J. Laming
C
      XL    = DF/(RHO*RHO13)
      AAA   = FT*RHOGRD/RHO
      ALPHA = DP5*RHOGRD*XL
      EXV   = sqrt(D1 + ALPHA*ALPHA)
      ASH   = log(ALPHA + EXV)
      B     = BETAP/(D1 + D6*BETAP*ALPHA*ASH)
      BBB   = D6*ALPHA*B*(ASH+ALPHA/EXV)
      VB2   = (BBB - D2)*B*XL
      VB1   = (D1 - BBB)*ALPHA*B*AAA
      RETURN
      END


C*****************************************************************************
      SUBROUTINE VBCK(VA,RHO,RHO13,RHOGRD,RHOLAP,RHOGHG)
C*****************************************************************************
C     Evaluates directly the potential 
C*****************************************************************************
#include "implicit.h"
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, D4 = 4.0D0, 
     &           D6 = 6.0D0, D12 = 12.0D0, DP5 = 0.5D0)
      PARAMETER (BETAP=0.0042D0)
      DF = D2**(D1/D3)
C
      XL    =  DF*D2/(RHO*RHO13)
      XK    =  DP5*RHOGRD
      G     =  XL*XL*XK*XK
      EXV   =  sqrt(D1+G)
      EXV2  =  D6*BETAP*XK*XL
      ASH   =  log(XK*XL+EXV)
      A     = -BETAP*XK*XK*XL
      B     =  D1+EXV2*ASH
      BSQ   =  B*B
      BCD   =  BSQ*B
      REC   =  D1/EXV
      RECD  =  REC*REC*REC
      AK    = -EXV2/D3
      BK    =  EXV2/XK*ASH+EXV2*XL*REC
      AL    =  A/XL
      BL    =  EXV2/XL*ASH+EXV2*XK*REC
      BKK   =  D12*BETAP*XL*XL*REC-EXV2*XK*XL**3*RECD
      BKL   =  D6*BETAP*ASH+D3*EXV2*REC - EXV2*G*RECD
      AKL   = -D2*BETAP*XK
      AKK   = -D2*BETAP*XL
      U     =  (AK/B-A*BK/BSQ)/XK
      UL    =  (AKL/B-BL*AK/BSQ-AL*BK/BSQ-A*BKL/BSQ+D2*A*BK*BL/BCD)/XK
      UK    = -U/XK+(AKK/B-D2*AK*BK/BSQ-A*BKK/BSQ+D2*A*BK*BK/BCD)/XK
      VAG   =  D2*D4/(D3*RHO)
      TERM1 = -(AL/B-A*BL/BSQ)*XL*VAG
      TERM2 = -DP5*U*RHOLAP
      TERM3 = -UK*(DP5**3)*RHOGHG/XK
      TERM4 =  UL*VAG*G/XL
      VA    =  TERM1+TERM2+TERM3+TERM4
      RETURN
      END


C*****************************************************************************
      SUBROUTINE V1BCK(DF1000,DF0010,DF2000,DF1010,DF0020,
     &                 RHOINP,RHOGRD)
C*****************************************************************************
C     (C) A. M. Lee, Cambridge 1993.
C
C     Calculates the second and lower derivatives of the Becke '88
C     exchange functional, wrt the density and related functions.
C     work with alpha-spin density, RHO; and length of the
C     alpha-spin gradient vector, GRHO. 
C
C*****************************************************************************
#include "implicit.h"
      PARAMETER (D0 = 0.0D0, DP5 = 0.50D0, D1 = 1.0D0, D2 = 2.0D0, 
     &           D3 = 3.0D0, D4 = 4.0D0, D43 = D4/D3, D6 = 6.0D0, 
     &           D7 = 7.0D0, BETAP = 0.0042D0)
      GRHO = DP5*RHOGRD
      RHO  = DP5*RHOINP
      IF (RHO .LT. 1D-10) THEN
         DF1000 = D0 
         DF0010 = D0
         DF1010 = D0
         DF2000 = D0
         DF0020 = D0
      ELSE 
         X = GRHO*RHO**(-D43)
         DENX2 = D1/(D1 + X*X)
         DENX1 = sqrt(DENX2)
         DAS   = log(X + sqrt(D1 + X*X))
         FX = -(DAS + X*DENX1)
         HX = -(D2 + X*X)*DENX1*DENX2
         GX = (D4 + X*X)*X*DENX1*DENX2*DENX2
C
         X10 = -D43/RHO
         X01 = D1/GRHO
C
         F2   = D1/(D1 + D6*BETAP*X*DAS)
         F    = -BETAP*RHO**(D43)*X*X*F2
         F3   = D6*BETAP*X*F2*FX
         TMP  = F*(F3 + D1)
         DF10 = X10*TMP
         DF01 = X01*(TMP + F)
         F4   = F3*F3 + F3 + D6*BETAP*F2*X*X*HX
         DF310=X10*F4
C
         DF1000 = DF10
         DF0010 = DF01
         DF1010 = DP5*((D2 + F3)*DF10 + F*DF310)*X01
         DF2000 = DP5*(DF10*(D7/D4 + F3) + F*DF310)*X10
         DF0020 = DP5*(DF01*(D1 + F3) + F*X01*F4)*X01
      END IF
      RETURN
      END


C     ****************************
C     ***** Lee, Yang & Parr *****
C     ****************************

C*****************************************************************************
      SUBROUTINE ELYP(CLYPCI,RHO,RHO13,RHOGRD)
C*****************************************************************************
C     *************WRITTEN BY G. J. LAMING********************
C*****************************************************************************
#include "implicit.h"
#include "pi.h"
      PARAMETER (DP5 = 0.50D0,
     &           ALL=0.04918D0,
     &           BLL=0.132D0,
     &           CLL=0.2533D0,
     &           DLL=0.349D0,
     &           ELL=(47.0D0/18.0D0),
     &           FLL=(7.0D0/18.0D0),
     &           GLL=(5.0D0/2.0D0),
     &           HLL=(1.0D0/18.0D0))
      DF  = DP5**(1.0D0/3.0D0)
      C333=2.0D0**(11.0D0/3.0D0)
      CFFF=0.3D0*((3.0D0*PI*PI)**(2.0D0/3.0D0))*C333
C     LEE,PARR AND YANG'S 'CORRELATION FUNCTIONAL:
      AMOSQ  = RHOGRD**2
      RHOA   = DP5*RHO
      AMASQ  = DP5*DP5*AMOSQ
      RHOA13 = DF*RHO13
      TM    = 1.0D0/RHO13
      AU    = 1.0D0+DLL*TM
      RRR   = RHO*RHO
      W     = (RHO13/RRR)*exp(-CLL*TM)/(RRR*AU)
      DT    = CLL*TM+DLL*TM/AU
      C1    = ALL*BLL*W
      C2    = C1*RHOA*RHOA
      C5    = (AMASQ+AMASQ)
      C6    = RHOA*AMASQ/RHO
      C7    = C6
      CCC   = (2.0D0/3.0D0)*RRR
      C8    = (CCC-RHOA*RHOA)*AMASQ
      C9    = (CCC-RHOA*RHOA)*AMASQ
      TERM1 = -ALL*(4.0D0/AU)*RHOA*RHOA/RHO
      TERM2 = -C2*CFFF*((RHOA*RHOA*RHOA/RHOA13)+(RHOA*RHOA*RHOA/RHOA13))
      TERM3 = -C2*(ELL-FLL*DT)*AMOSQ
      TERM4 =  C2*(GLL-HLL*DT)*C5
      TERM5 =  C2*((DT-11.0D0)/9.0D0)*(C6+C7)
      TERM6 =  C1*(2.0D0/3.0D0)*RRR*AMOSQ
      TERM7 = -C1*C8
      TERM8 = -C1*C9
      CLYPCI=TERM1+TERM2+TERM3
     &      +TERM4+TERM5+TERM6
     &      +TERM7+TERM8
      RETURN
      END


C*****************************************************************************
      SUBROUTINE GLYP(TERM1,U,RHO,RHO13,RHOGRD)
C*****************************************************************************
C
C***********WRITTEN   BY   G. J. LAMING******************
C
C*****************************************************************************
#include "implicit.h"
#include "pi.h"
      PARAMETER (ALL=0.04918D0,
     &           BLL=0.132D0,
     &           CLL=0.2533D0,
     &           DLL=0.349D0,
     &           ELL=(47.0D0/18.0D0),
     &           FLL=(7.0D0/18.0D0),
     &           GLL=(5.0D0/2.0D0),
     &           HLL=(1.0D0/18.0D0),
     &           PLL=(5.0D0/18.0D0))
      CFFFF=0.3D0*((3.0D0*PI*PI)**(2.0D0/3.0D0))
      
      RRR   = RHO*RHO
      BG    = ALL*BLL
      RRRR  = RRR*RHO
      XXX   = BG*CFFFF
      RHOT  = 1.0D0/RHO13
      XB    = 1.0D0+DLL*RHOT
      AG    = RHOT/(3.0D0*RHO)
      XB2   = XB*XB
      XBRHO = -DLL*AG
      XD    = exp(-CLL*RHOT)
      XDRHO = CLL*XD*AG
      XA    = -ALL*RHO-XXX*RHO*XD
      XARHO = -ALL-XXX*XD-XXX*RHO*XDRHO
      TERM1 = XARHO/XB-XBRHO*XA/XB2
      XE    = FLL*(CLL/(RRR)+DLL/(RRR*XB))
     &      + RHO13/(6.0D0*RRR)
      XERHO = FLL*(-2.0D0*CLL/RRRR-2.0D0*DLL/(RRRR*XB)
     &      - XBRHO*DLL/(RRR*XB2)) - PLL*RHO13/RRRR
      XF    = (BG/4.0D0)*XD*XE/XB
      XFRHO = (BG/4.0D0)*(XDRHO*XE/XB+XD*XERHO/XB-XD*XE*XBRHO/XB2)
      TERM1 = TERM1+RHOGRD*RHOGRD*XFRHO
      U     = 4.0D0*XF
      RETURN
      END


C*****************************************************************************
      SUBROUTINE VLYP(VA,RHO,RHO13,RHOGRD,RHOLAP)
C*****************************************************************************
C
C     Lee
C  
C*****************************************************************************
#include "implicit.h"
#include "pi.h"
      PARAMETER (ALL=0.04918D0,BLL=0.132D0,CLL=0.2533D0,DLL=0.349D0,
     &           FLL=(7.0D0/18.0D0),PLL=(5.0D0/18.0D0))
      CFFFF=0.3D0*((3.0D0*PI*PI)**(2.0D0/3.0D0))
      
      RRR    = RHO*RHO
      RRRR   = RHO*RRR
      RHOT   = 1.0D0/RHO13
      XB     = 1.0D0+DLL*RHOT
      AG     = RHOT/(3.0D0*RHO)
      BG     = ALL*BLL
      BG2    = BG*CFFFF
      XB2    = XB*XB
      XBRHO  = -DLL*AG
      XD     = exp(-CLL*RHOT)
      XDRHO  = CLL*XD*AG
      XA     = -ALL*RHO-BG2*RHO*XD
      XARHO  = -ALL-BG2*XD-BG2*RHO*XDRHO
      TERM1  = XARHO/XB-XBRHO*XA/XB2
      XE     = FLL*(CLL/RRR+DLL/(RRR*XB)) + RHO13/(6.0D0*RRR)
      XERHO  = FLL*(-2.0D0*CLL/RRRR-2.0D0*DLL/(RRRR*XB)
     &       - XBRHO*DLL/(RRR*XB2)) - PLL*RHO13/RRRR
      XF     = (BG/4.0D0)*XD*XE/XB
      XFRHO  = (BG/4.0D0)*(XDRHO*XE/XB+XD*XERHO/XB-XD*XE*XBRHO/XB2)
      TERM24 = -RHOGRD*RHOGRD*XFRHO
      TERM3  = -2.0D0*RHOLAP*XF
      VA     = TERM1+TERM24+TERM3
      RETURN
      END


C*****************************************************************************
      SUBROUTINE V1LYP(DF1000,DF0010,DF2000,DF1010,DF0020,RHOINP,RHOGRD)
C*****************************************************************************
C     (C) A. M. Lee Cambridge 1993/94
C
C     The second and lower derivatives of the LYP functional,wrt
C     the density and related quantities.
C 
C     The following is closed shell code and cannot be generalised.
C     Care has been taken to avoid numerical instabilities.
C
C*****************************************************************************
#include "implicit.h"
#include "pi.h"
      PARAMETER (DP5=0.050,
     &           ALL=0.04918D0,BLL=0.132D0,CLL=0.2533D0,DLL=0.349D0)
      PARAMETER(C13=1.0D0/3.0D0,C19=1.0D0/9.0D0,C49=4.0D0/9.0D0,
     &          C23=2.0D0/3.0D0,C43=4.0D0/3.0D0,C512=5.0D0/12.0D0,
     &          C524=5.0D0/24.0D0,C548=5.0D0/48.0D0,C506=5.0D0/6.0D0,
     &          C172=1.0D0/72.0D0,C136=1.0D0/36.0D0,C118=1.0D0/18.0D0)
      C333=2.0D0**(11.0D0/3.0D0)
      A=0.04918D0
      B=0.132D0
      C=0.2533D0
      D=0.349D0
      A=ALL
      B=BLL
      C=CLL
      D=DLL
      AB=A*B
      CF=0.3D0*(3.0D0*PI*PI)**C23
      DS2=sqrt(2.0D0)
C
      rhoa=0.50D0*RHOINP
      IF (abs(RHOA).LT.1D-10 .OR. abs(RHOGRD).LT.1D-10) THEN
      DF1000=0.0D0
      DF0010=0.0D0
      DF1010=0.0D0
      DF2000=0.0D0
      DF0020=0.0D0
      ELSE 
      RHO=RHOA+RHOA
      GRHOA=0.50D0*RHOGRD
      grho=grhoa*2.0D0
C
      rhoa2=rhoa*rhoa
      rhon1=1.0D0/rho
      rhon2=rhon1*rhon1
      rho2=rho*rho
C
      rho13=rho**(1.0D0/3.0D0)
      rho13=1.0D0/rho13
      rho23=rho13*rho13
      rho43=rho23*rho23
      rho113=rho23*rhon2*rhon1
      crho13=c*rho13
      drho13=d*rho13
      crho43=c*rho43
      drho43=d*rho43
C
      grho2=grho*grho
C
      X10=-4.0D0/3.0D0/rho
      X01=1.0D0/grho
C
      eps=1.0D0/(1.0D0+drho13)
      eps2=eps*eps
C
C     F=-a*f1-a*b*(omega*(rhoa*rhoa*G2+G3))
C     Functions of rho(rhoa) only: eps,delta,omega,f1,f2,f3
C     Functions of rho,grho(rho,grhoa): f4-f9
C
      Deps10=eps2*c13*drho43
      Deps20=c13*drho43*eps*(2.0D0*deps10+X10*eps)
C
C     Deps20/eps-Deps10/eps/eps Appears unstable. Use:
C
      omega1=exp(-crho13)
      omega2=rho23*rhon2*rhon1
      omega=eps*omega1*omega2

C     Examine carefully:

      Dw10=(c13*crho43+Deps10/eps
     &-11.0D0/3.0D0*rhon1)

      Dw10x2=omega*Dw10*Dw10
      Dw10=omega*Dw10

C    Dw10x2=Dw10*dw10/omega is also unstable:
C    lim rho->big, dw10/omega->0/0

      alt=c13*drho43*(deps10+X10*eps)
      Dw20=Dw10x2+omega*(
     &X10*c13*crho43+11.0D0/3.0D0*rhon2+alt)



      delta=drho13*eps+crho13
      Ddel10=-c13*drho43*eps-c13*crho43+
     &drho13*Deps10
      Ddel20=c49*rhon1*rho43*(d*eps+c)
     &-2.0D0*c13*drho43*Deps10+drho13*Deps20

      f1=rho*eps
      Df110=(eps+rho*Deps10)
      Df120=2.0D0*Deps10+rho*Deps20

      f2=4.0D0*Cf*rho**(8.0D0/3.0D0)
      Df210=8.0D0/3.0D0*rhon1*f2
      Df220=40.0D0/9.0D0*rhon2*f2

      parent=(47.0D0-7.0D0*delta)/18.0D0
      f3=parent*grho2
      Df301=parent*2.0D0*grho
      Df302=Df301/grho
      Df310=-7.0D0/18.0D0*grho2*Ddel10
      Df320=-7.0D0/18.0D0*grho2*Ddel20
      Df311=-7.0D0/9.0D0*grho*Ddel10

      parent=(delta*c136-1.25D0)
      f4=parent*grho2
      Df410=grho2*Ddel10*c136
      Df420=grho2*Ddel20*c136
      Df401=2.0D0*parent*grho
      Df402=2.0D0*parent
      Df411=c118*grho*Ddel10

C     f5=f5+f6

      f5=-C136*(delta-11.0D0)*grho2
      Df501=-C118*(delta-11.0D0)*grho
      Df502=-C118*(delta-11.0D0)
      Df510=-C136*grho2*Ddel10
      Df520=-C136*grho2*Ddel20
      Df511=-C118*grho*Ddel10
      
      f7=-c23*rho2*grho2
      Df710=2.0D0*f7/rho
      Df720=2.0D0*f7/rho2
      Df701=2.0D0*f7/grho
      Df702=2.0D0*f7/grho2
      Df711=-4.0D0*c23*rho*grho

C     f8=f8+f9
      f8=c524*rho2*grho2
      Df810=c512*rho*grho2
      Df820=c512*grho2
      Df801=c512*rho2*grho
      Df802=c512*rho2
      Df811=c506*rho*grho

      G2=f2+f3+f4+f5
      G3=f7+f8

      DG210=Df210+Df310+Df410+Df510
      DG220=Df220+Df320+Df420+Df520
      DG201=Df301+Df401+Df501
      DG202=Df302+Df402+Df502
      DG211=Df311+Df411+Df511

      DG310=Df710+Df810
      DG320=Df720+Df820
      DG301=Df701+Df801
      DG302=Df702+Df802
      DG311=Df711+Df811

      DF10=omega*(DG310+rhoa*(G2+rhoa*DG210))+
     &Dw10*(G3+G2*rhoa2)
      DF10=-(a*DF110+ab*DF10)

      DF11=
     &    omega*(DG311+rhoa*(DG201+rhoa*DG211))+
     &Dw10*(DG301+DG201*rhoa2)
      DF11=-(ab*DF11)

      DF20=DG320*omega+DG310*Dw10+DG310*Dw10+G3*Dw20+
     &DG210*rhoa2*dw10+G2*rhoa*dw10+dw20*rhoa2*G2+
     &(G2+rhoa*DG210)*(0.50D0*omega+rhoa*Dw10)+
     &rhoa*omega*(DG210+0.50D0*DG210+rhoa*DG220)
      DF20=-(a*DF120+ab*DF20)

      DF01=(DG301+rhoa2*DG201)*omega
      DF01=-(ab*DF01)

      DF02=(DG302+rhoa2*DG202)*omega
      DF02=-(ab*DF02)

      F=-a*f1-a*b*(omega*(rhoa2*G2+G3))

C     DF1000 scaled for tozer

      DF1000=DF10
      DF0010=DF01
      DF1010=DF11
      DF2000=DF20
      DF0020=DF02
      END IF
      RETURN
      END


C*****************************************************************************
      subroutine vtlyp(F20000,F02000,F11000,F10100,F01010,
     &   F10010,F01100,F10001,F01001,pa,pb,gaa,gbb,gab)
C*****************************************************************************

      implicit none

      double precision a,c,d,cod,ab,st,x,y,xop,yop,gt
      double precision gaapbb,gaambb,xgaapygbb,q,r,qd
      double precision rd,qdd,rdd,hm1,k,kd,kdd,del,deld
      double precision deldd,w,part,imp,wd,wdd,s1,s2,s3
      double precision s1d,s2d,s3d,s1dd,s2dd,s3dd,F01001
      double precision s12,s12d,i1,i2,i3,i4,i5,i6,rp11
      double precision F10100,F01100,F01010,F10010,pa23
      double precision i9,i10,i11,i12,bs,j2,j3,bsda,bsdb
      double precision bsdadb,bsdada,j1,i13,i14,i15,bt,F02000
      double precision i8,i7,btda,btdb,btdadb,btdada,qwe
      double precision phm1,phm1sq,phm1cu,z,z2,trans,same
      double precision jab,jaa,abw,abwd,F20000,F11000,same2
      double precision  F10001,Co,p2,gaa,gbb,p,ta,tb,ert,yui 
      double precision pa,pb,gab,pa53,pb53,pa83,pb83,pa2,pb2
      double precision btdbdb

c     ****data for Lyp
      Co=  36.46239897876477 d0
      a=    0.04918          d0
c     b=    0.132            d0
      c=    0.2533           d0
      d=    0.349            d0
      cod=  0.72578796561604 d0
      ab=   0.00649176       d0

c     ***some definitions using the basic quantities

      p=pa+pb
      st=pa*pb
      x=pa/p
      y=pb/p
      xop=x/p
      yop=y/p

      gt=gaa+gbb+2.0d0*gab
      gaapbb=gaa+gbb
      gaambb=gaa-gbb
      xgaapygbb=x*gaa+y*gbb

c     **** powers of pa,pb,p that are needed

      p2=p*p
      
      pa53=pa**(5.0d0/3.0d0)
      pb53=pb**(5.0d0/3.0d0)
      pa23=pa53/pa
c     pb23=pb53/p
      pa83=pa53*pa
      pb83=pb53*pb
      pa2=pa*pa
      pb2=pb*pb
      rp11=p**(-11.0d0/3.0d0)
      
c     ****delta and its first and second derivatives

      q=p**(-1.0d0/3.0d0)
      r=d*q
      qd=-1.0d0/3.0d0*q/p
      rd=d*qd
      qdd=-4.0d0/3.0d0*qd/p
      rdd=d*qdd
      hm1=1/(1+r) 

      k=r*hm1
      kd=rd*hm1
      kdd=rdd*hm1

      del=cod*r+k
      deld=cod*rd+kd*hm1
      deldd=cod*rdd+hm1*(kdd-2.0d0*kd*kd)

c     ****omega and its first and second derivatives

      w=exp(-c*q)*rp11*hm1
      part=11/q
      imp=(part-c-d*hm1)
      wd=w*qd*imp
c     s implify 2nd deriv with intermediates
      wdd=wd*imp*qd+w*qdd*imp+w*qd*qd*(-part/q+d*d*hm1*hm1)

c     ****sigma 1,2,3 and their first and second derivatives

      s1=(47.0d0/18.0d0-7.0d0/18.0d0*del)
      s2=(1.0d0/18.0d0*del-5.0d0/2.0d0)
      s3=(11.0d0/9.0d0-del/9.0d0)

      s1d=-7.0d0/18.0d0*deld
      s2d=1.0d0/18.0d0*deld
      s3d=-deld/9.0d0
      s1dd=-7.0d0/18.0d0*deldd
      s2dd=1.0d0/18.0d0*deldd
      s3dd=-deldd/9.0d0

c     ****derivatives of density followed by a gradient
      
      abw=ab*w

      s12=s1+s2
      s12d=s1d+s2d

      i1=s12+s3*x
      i2=s12d+s3d*x
      i3=s12+s3*y
      i4=s12d+s3d*y 
      i5=s3*yop
      i6=s3*xop
      
      F10100=-ab*(st*(wd*i1+w*(i2+i5))+pb*w*i1-wd*pb2)
      F01010=-ab*(st*(wd*i3+w*(i4+i6))+pa*w*i3-wd*pa2)
      F01100=-ab*(st*(wd*i1+w*(i2-i6))+pa*w*i1-wd*pb2-2.0d0*w*pb)  
      F10010=-ab*(st*(wd*i3+w*(i4-i5))+pb*w*i3-wd*pa2-2.0d0*w*pa)
      
c     F10100=-ab*(st*(wd*(i1)+w*(i2+i5))+pb*w*(i1)+wd*pb*pb)
c     F01010=-ab*(st*(wd*(i3)+w*(i4+i6))+pa*w*(i3)+wd*pa*pa)
c     F01100=-ab*(st*(wd*(i1)+w*(i2-i6))+pa*w*(i1)+wd*pb*pb+2*w*pb)  
c     F10010=-ab*(st*(wd*(i3)+w*(i4-i5))+pb*w*i3+wd*pa*pa+2*w*pa)      

c     ****definition of S and its first and second derivatives

      i9=Co*(pa83+pb83)
      i10=s1*gt
      i11=s2*gaapbb
      i12=s3*xgaapygbb

      bs=-ab*(i9+i10+i11+i12)
      ert=Co*(8.0d0/3.0d0)
      yui=s1d*gt+s2d*(gaapbb)+s3d*(xgaapygbb)
      j2=ert*pa53+yui
      j3=s1dd*gt+s2dd*gaapbb+s3dd*(xgaapygbb)
      bsda=-ab*(ert*pa53+yui+s3*y/p*gaambb)
      bsdb=-ab*(ert*pb53+yui+s3*x/p*(-gaambb))
      bsdadb=-ab*(j3+s3d*(y-x)*gaambb/p+(x-y)/p/p*s3*gaambb)
      bsdada=-ab*(Co*40/9*pa23+j3+s3d*2*y/p*(gaambb)
     &       +s3*(gaambb)*(-2.0d0*y/p/p))

c Timings only
      bsdada=-ab*(Co*40/9*pa23+j3+s3d*2*y/p*(gaambb)
     &       +s3*(gaambb)*(-2.0d0*y/p/p))

c     ****definition of T and its first and second derivatives

      j1=(2.0d0/3.0d0)*p2
      i13=-j1*gt
      i14=(j1-pa2)*gbb
      i15=(j1-pb2)*gaa
      qwe=(i13+i14+i15)
      bt=-ab*qwe

      i8=-8.0d0/3.0d0*gab
      i7=i8*p

      btda=-ab*(i7-2.0d0*pa*gbb)
      btdb=-ab*(i7-2.0d0*pb*gaa)
      btdadb=-ab*i8
      btdada=-ab*(i8-2.0d0*gbb)
c     btda=-ab*((i7-2.0d0*pa*gbb)*w+qwe*wd)
c     btdb=-ab*((i7-2.0d0*pb*gaa)*w+qwe*wd)
c     btdada=-ab*(wdd*qwe+2.0d0*(i7-2.0d0*pa*gbb)*wd+w*(i8-2.0d0*gbb))
c     btdadb=-ab*(wdd*qwe+wd*((i7-2.0d0*pb*gaa)+(i7-2.0d0*pa*gbb))
c    1            +w*(i8))
c Timings only
      btdbdb=i8-2.0d0*gaa

c     ****definition of F and its second derivatives

      phm1=hm1/p
      phm1sq=phm1*phm1
      phm1cu=phm1sq*phm1
      z=1.0d0+2.0d0/3.0d0*r
      z2=z*z
c      trans constant takes other constant into account -dubious
      trans=2.0d0/3.0d0*rd
      same=-(p*z+st*trans)*phm1sq+2*st*z2*phm1cu
      same2=-(2.0d0*pb*z+st*trans)*phm1sq+2*st*z2*phm1cu
      jab=-4.0d0*a*(phm1+same)
      jaa=-4.0d0*a*same2
      ta=-pb*z*phm1sq
      tb=phm1-pb*z*phm1sq
c     ****Assembly of second derivatives involving pa,pb

      abwd=ab*wd
      F20000=jaa+wdd*(st*bs+bt)+2.0d0*wd*(pb*bs+st*bsda+btda)
     &       +w*(2.0D0*pb*bsda+st*bsdada+btdada)

c For timings only

      F02000=jaa+wdd*(st*bs+bt)+2.0d0*wd*(pb*bs+st*bsda+btda)
     &       +w*(2.0D0*pb*bsda+st*bsdada+btdada)

c End
 
      F11000=jab+wdd*(st*bs+bt)
     &       +wd*(pa*bs+st*bsdb+btdb+pb*bs+st*bsda+btda)
     &       +w*(bs+pb*bsdb+pa*bsda+st*bsdadb+btdadb)

c    ****Second derivatives involving gab and a density

      F10001=-abw*2.0d0*(pb*s1+st*s1d-4.0d0/3.0d0*p)
     &       -ab*wd*2.0d0*(st*s1-2.0d0/3.0d0*p2)
      F01001=-abw*2.0d0*(pa*s1+st*s1d-4.0d0/3.0d0*p)
     &       -ab*wd*2.0d0*(st*s1-2.0d0/3.0d0*p2)


c    ****All other derivatives are aero for Lyp
      return
      end


C*****************************************************************************
      SUBROUTINE GLYPCO(EASYGA,EASYGB,HARDGA,HARDGB,HARDER,
     &                  RHO,RHO13,RHOGRD,UHF)
C*****************************************************************************
C
C
C     OPENSHELL LYP CORRELATION  GRADIENT CODE
C
C***********WRITTEN   BY   G. J. LAMING******************
C
C*****************************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (DP5 = 0.50D0)
      LOGICAL UHF
C
      IF (.NOT.UHF) THEN
         CALL GCSLYP(TERM1,U,RHO,RHO13,RHOGRD)
         EASYGA=TERM1
         HARDGA=2.0D0*U
      ELSE
         RHOA   = DP5*RHO
         RHOB   = RHOA 
         RHOA13 = RHOA**(1.0D0/3.0D0)
         RHOB13 = RHOA13
         RHOA53 = RHOA*RHOA/RHOA13
         RHOB53 = RHOB*RHOB/RHOB13
         RHOA83 = RHOA53*RHOA
         RHOB83 = RHOB53*RHOB
         RHO43  = RHO*RHO13
         RHO23  = RHO/RHO13
         RHO113 = RHO*RHO*RHO*RHO/RHO13
         RHO143 = RHO113*RHO
         AMA    = DP5*RHOGRD
         AMB    = AMA 
         AMASQ  = AMA*AMA 
         AMBSQ  = AMASQ 
         AMAAMB = AMASQ 
         CALL GOSLYP(A1,A2,A3,A4,A5,RHOA,RHOB,RHOA53,RHOB53,RHOA83,
     &               RHOB83,RHO13,RHO23,RHO43,RHO113,RHO143,
     &               AMASQ,AMBSQ,AMA,AMB,AMAAMB)
C
C
         EASYGA=A1
         EASYGB=A2
         HARDGA=A3
         HARDGB=A4
         HARDER=A5
      ENDIF
C
      RETURN
      END


C*****************************************************************************
      SUBROUTINE GCSLYP(TERM1,U,RHO,RHO13,GMOD)
C*****************************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
#include "pi.h"
      COMMON/DFTCON/BCX,CVX,AL,BL,CL,DL,CF,C3,CRS,PRE,THIRD,CINF,FF,AA
     &,BB,GG,DD,AK1,AK2,AK3,PP,FT,QQ,BETAP,CV,AI,BI,CI,X0I,AII,BII,CII
     &,X0II,AIII,BIII,CIII,X0III,FPPZ,FTT,FTTT,QQQI,QQQII,QQQIII,XF0I
     &,YF0I,XF0II,YF0II,XF0III,YF0III,DCRS,PD,Q,P,ALL,BLL,CLL,DLL,CFFF
     &,ELL,FLL,GLL,HLL,PLL,CFFFF,ONE,TWO,THREE,FOUR,FIVE,EIGHT,TEN
     &,ELEVEN,VAR6,POW1,VAR14,VAR2,CFG,C2G,PRODAB,C3G,C33G,C4G,C44G
     &,C5G,C55G,C6G,C66G,C7G,C77G,C8G,C88G,C9G,C99G,C11G,C1111G
      RRR=RHO*RHO
      BG=ALL*BLL
      RRRR=RRR*RHO
      XXX=BG*CFFFF
      RHOT=1.0D0/RHO13
      XB=1.0D0+DLL*RHOT
      AG=RHOT/(3.0D0*RHO)
      XB2=XB*XB
      XBRHO=-DLL*AG
      XD=exp(-CLL*RHOT)
      XDRHO=CLL*XD*AG
      XA=-ALL*RHO-XXX*RHO*XD
      XARHO=-ALL-XXX*XD-XXX*RHO*XDRHO
      TERM1=XARHO/XB-XBRHO*XA/XB2
      XE=FLL*(CLL/(RRR)+DLL/(RRR*XB))
     &  +RHO13/(6.0D0*RRR)
      XERHO=FLL*(-2.0D0*CLL/RRRR-2.0D0*DLL/(RRRR*XB)
     &                    -XBRHO*DLL/(RRR*XB2))
     &     -PLL*RHO13/RRRR
      XF=(BG/4.0D0)*XD*XE/XB
      XFRHO=(BG/4.0D0)*(XDRHO*XE/XB+XD*XERHO/XB-XD*XE*XBRHO/XB2)
C
      TERM1=TERM1+GMOD*GMOD*XFRHO
      U=2.0D0*XF
      RETURN
      END


C*****************************************************************************
      SUBROUTINE GOSLYP(A1,A2,A3,A4,A5,RHOA,RHOB,RHOA53,RHOB53
C*****************************************************************************
C     OPENSHELL POTENTIAL FOR THE GRADIENT OF LYP FUNCTIONAL ENERGY
C     CODED UP BY G.J.LAMING, JUNE 1991
C*****************************************************************************
     &,RHOA83,RHOB83,RHO13,RHO23,RHO43
     &,RHO113,RHO143,AMASQ,AMBSQ,AMA,AMB,AMAAMB)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
#include "pi.h"
      COMMON/DFTCON/BCX,CVX,AL,BL,CL,DL,CF,C3,CRS,PRE,THIRD,CINF,FF,AA
     &,BB,GG,DD,AK1,AK2,AK3,PP,FT,QQ,BETAP,CV,AI,BI,CI,X0I,AII,BII,CII
     &,X0II,AIII,BIII,CIII,X0III,FPPZ,FTT,FTTT,QQQI,QQQII,QQQIII,XF0I
     &,YF0I,XF0II,YF0II,XF0III,YF0III,DCRS,PD,Q,P,ALL,BLL,CLL,DLL,CFFF
     &,ELL,FLL,GLL,HLL,PLL,CFFFF,ONE,TWO,THREE,FOUR,FIVE,EIGHT,TEN
     &,ELEVEN,VAR6,POW1,VAR14,VAR2,CFG,C2G,PRODAB,C3G,C33G,C4G,C44G
     &,C5G,C55G,C6G,C66G,C7G,C77G,C8G,C88G,C9G,C99G,C11G,C1111G
C
      RHO=RHOA+RHOB
      RR=RHO*RHO
      PAB=RHOA*RHOB
      VAR=1.0D0/RHO13
      VAR3=1.0D0/RHO43
      VAR4=1.0D0/RHO23
      VAR5=RHO13/RR
      VAR9=RHOA83+RHOB83
      VAR10=VAR6*RHOA53
C     CALCULATION OF W ABD IT'S DERIVATIVE WRT DENSITY
      GGGG=exp(-CLL*VAR)
      AW=GGGG/RHO113
      BW=ONE+DLL*VAR
      W=AW/BW
      AWDA=(CLL/THREE)*VAR3*AW-GGGG*(ELEVEN/THREE)/RHO143
      BWDA=-DLL*VAR3/THREE
      WDA=(AWDA/BW)-(AW*BWDA/(BW*BW))
      WDB=WDA
C     CALCULATION OF DELTA AND IT'S DERIVATIVE WRT DENSITY
      DTA=CLL*VAR+(DLL*VAR)/BW
      ADTA=(CLL+DLL)*VAR+DLL*CLL*VAR4
      BDTA=BW
      ADTADA=-((CLL+DLL)*VAR3/THREE)-(POW1*DLL*CLL*VAR5)
      BDTADA=BWDA
      DTADA=(ADTADA/BDTA)-(ADTA*BDTADA/(BDTA*BDTA))
      DTADB=DTADA
      ZA=(POW1*RR)-RHOA*RHOA
      ZB=(POW1*RR)-RHOB*RHOB
      Q1=RHOA/RHO
      Q2=RHOB/RHO
      Q1DA=(RHOB)/RR
      Q1DB=(-RHOA)/RR
      Q2DA=-Q1DA
      Q2DB=-Q1DB
      XDA=AMASQ*Q1DA-AMBSQ*Q1DA
      XDB=AMASQ*Q1DB-AMBSQ*Q2DA
      ZADA=POW1*(TWO*RHOB-RHOA)
      ZADB=FT*RHO
      ZBDA=ZADB
      ZBDB=POW1*(TWO*RHOA-RHOB)
      F3=C3G*PAB*W
      F4=C4G*PAB*W*DTA
      F5=C5G*PAB*W
      F6=C6G*PAB*W*DTA
      F7A=C7G*PAB*W*DTA*Q1
      F7B=C7G*PAB*W*DTA*Q2
      F8A=C8G*PAB*W*Q1
      F8B=C8G*PAB*W*Q2
      F9=C9G*W*RHO*RHO
      F10=PRODAB*W*ZA
      F11=C11G*W*ZB
      FB=-F3+F4+F5-F6+F7A-F8A+F9-F11
      FC=-F3+F4+F5-F6+F7B-F8B+F9-F10
      FD=TWO*(-F3+F4+F9)
      CALL NGLYP(A1,RHOA,RHOB,AMA,AMB,AMAAMB,W,WDA,WDB,DTA,DTADA,
     & RHOA53,RHOB53,RHOA83,RHOB83,
     & RHO13,RHO23,RHO43,
     & RHO113,RHO143,
     & F3,F4,F5,F6,
     &  VAR,VAR3,VAR4,VAR5,VAR9,RR,PAB,RHO,ZA,ZB,
     & Q1,Q2,Q1DA,Q2DA,XDA,XDB,ZADA,ZBDA,
     & AMASQ,AMBSQ)
      CALL NGLYP(A2,RHOB,RHOA,AMB,AMA,AMAAMB,W,WDA,WDB,DTA,DTADA,
     & RHOB53,RHOA53,RHOB83,RHOA83,
     & RHO13,RHO23,RHO43,
     & RHO113,RHO143,
     & F3,F4,F5,F6,
     & VAR,VAR3,VAR4,VAR5,VAR9,RR,PAB,RHO,ZB,ZA,
     & Q2,Q1,Q2DB,Q1DB,XDB,XDA,ZBDB,ZADB,
     & AMBSQ,AMASQ)
      A3=TWO*FB
      A4=TWO*FC
      A5=FD
      RETURN
      END


C*****************************************************************************
      SUBROUTINE NGLYP(T,RHOA,RHOB,AMA,AMB,AMAAMB,W,WDA,WDB,DTA,DTADA,
C*****************************************************************************
C     OPENSHELL POTENTIAL FOR THE GRADIENT OF LYP FUNCTIONAL ENERGY
C
C**********CODED   UP   BY   G.  J.  LAMING,   JUNE 1991*************
C
C*****************************************************************************
     & RHOA53,RHOB53,RHOA83,RHOB83,
     & RHO13,RHO23,RHO43,
     & RHO113,RHO143,
     & F3,F4,F5,F6,
     &  VAR,VAR3,VAR4,VAR5,VAR9,RR,PAB,RHO,ZA,ZB,
     & Q1,Q2,Q1DA,Q2DA,XDA,XDB,ZADA,ZBDA,
     & AMASQ,AMBSQ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
#include "pi.h"
      COMMON/DFTCON/BCX,CVX,AL,BL,CL,DL,CF,C3,CRS,PRE,THIRD,CINF,FF,AA
     &,BB,GG,DD,AK1,AK2,AK3,PP,FT,QQ,BETAP,CV,AI,BI,CI,X0I,AII,BII,CII
     &,X0II,AIII,BIII,CIII,X0III,FPPZ,FTT,FTTT,QQQI,QQQII,QQQIII,XF0I
     &,YF0I,XF0II,YF0II,XF0III,YF0III,DCRS,PD,Q,P,ALL,BLL,CLL,DLL,CFFF
     &,ELL,FLL,GLL,HLL,PLL,CFFFF,ONE,TWO,THREE,FOUR,FIVE,EIGHT,TEN
     &,ELEVEN,VAR6,POW1,VAR14,VAR2,CFG,C2G,PRODAB,C3G,C33G,C4G,C44G
     &,C5G,C55G,C6G,C66G,C7G,C77G,C8G,C88G,C9G,C99G,C11G,C1111G
      VAR10=VAR6*RHOA53
      A1=FOUR*ALL*PAB
      A1D=FOUR*ALL*RHOB
      B1=RHO+DLL*RHO23
      B1D=ONE+POW1*DLL*VAR
      F1D=(A1D/B1) - (A1*B1D/(B1*B1))
      E1D=F1D
      F2D=C2G*(RHOB*VAR9*W+PAB*VAR10*W+PAB*VAR9*WDA)
      E2D=F2D
      VV=(RHOA/RR)*PAB
      PAD=(TWO*PAB/RHO)-VV
      PBD=(RHOA*RHOA/RHO)-VV
      F7A=C7G*PAB*W*DTA*Q1
      F7B=C7G*PAB*W*DTA*Q2
      F8A=C8G*PAB*W*Q1
      F8B=C8G*PAB*W*Q2
      F9=C9G*W*RHO*RHO
      F10=PRODAB*W*ZA
      F11=C11G*W*ZB
      T3A=C3G*(RHOB*W+PAB*WDA)
      T4A=C4G*(WDA*PAB*DTA+W*RHOB*DTA+W*PAB*DTADA)
      T5A=C5G*(WDA*PAB+W*RHOB)
      T6A=C6G*(WDA*PAB*DTA+W*RHOB*DTA+W*PAB*DTADA)
      T7A=C7G*
     & (WDA*PAB*DTA*Q1+W*RHOB*DTA*Q1+W*PAB*DTADA*Q1+W*PAB*DTA*Q1DA)
      T7AA=C7G*
     &  (WDA*PAB*DTA*Q2+W*RHOB*DTA*Q2+W*PAB*DTADA*Q2+W*PAB*DTA*Q2DA)
      T8A=C8G*(WDA*PAB*Q1+W*RHOB*Q1+W*PAB*Q1DA)
      T8AA=C8G*(WDA*PAB*Q2+W*RHOB*Q2+W*PAB*Q2DA)
      T9A=C9G*(WDA*RR+W*TWO*RHO)
      T10A=PRODAB*(WDA*ZA+W*ZADA)
      T11A=C11G*(WDA*ZB+W*ZBDA)
      E1DA=E1D
      E2DA=E2D
      P1=(-E1DA-E2DA)
      SUMA1=-T3A+T4A+T5A-T6A+T7A-T8A+T9A-T11A
      P2=AMASQ*SUMA1
      SUMA2=-T3A+T4A+T5A-T6A+T7AA-T8AA+T9A-T10A
      P3=AMBSQ*SUMA2
      SUMA3=TWO*(-T3A+T4A+T9A)
      P4=AMAAMB*SUMA3
      T=P1+P2+P3+P4
      RETURN
      END


C*****************************************************************************
      SUBROUTINE ECSRLDA(ECV,RHO,RHO13,CHI,MULOCAL,ERFEXP)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Apr. 2003
C
C     Description :
C     Short-range LDA correlation-functional suitable for the erf 
C     and erf+exp type two-electron operator.
C
C     Based on implementation in Molpro by A. Savin
C
C*****************************************************************************
#include "implicit.h"
      LOGICAL MULOCAL(0:2),ERFEXP(*)
      REAL*8  ZKSR(3)
      CALL VCSRLDA(ZKSR,RHO,RHO13,CHI,MULOCAL,0,ERFEXP)
      ECV = ZKSR(1)
      RETURN
      END


C*****************************************************************************
      SUBROUTINE VCLALDA(ZKSR,RHO,RHO13,VLAMBDA)
C*****************************************************************************
C Kamal SHARKAS 17-05-11
C linear complement non-scaled LDA correlation-functional 
C*****************************************************************************
#include "implicit.h"
      REAL*8  ZKSR(3), RHO, RHO13, VLAMBDA
      LOGICAL MULOCAL(0:2),ERFEXP(0:2)

      MULOCAL(0:2) = .FALSE.
      ERFEXP(0:2)  = .FALSE.

      CALL VCSRLDA(ZKSR,RHO,RHO13,0.D0,MULOCAL,1,ERFEXP)
      ZKSR(1) = (1.D0 - VLAMBDA**2) * ZKSR(1)
      ZKSR(2) = (1.D0 - VLAMBDA**2) * ZKSR(2)
      ZKSR(3) = (1.D0 - VLAMBDA**2) * ZKSR(3)
      RETURN
      END


C*****************************************************************************
      SUBROUTINE VCLASLDA(ZKSR,RHO,RHO13,VLAMBDA)
C*****************************************************************************
C linear complement scaled LDA correlation-functional
C*****************************************************************************
#include "implicit.h"
      REAL*8  ZKSR(3),ZKSRCOUL(3),ZKSRSCALED(3)
      REAL*8  RHOSCALED,RHO13SCALED
      LOGICAL MULOCAL(0:2),ERFEXP(0:2)

      MULOCAL(0:2)  = .FALSE.
      ERFEXP(0:2)   = .FALSE.
      ZKSRCOUL(1:3) = 0.d0
      ZKSRSCALED(1:3) = 0.d0
      RHOSCALED     = RHO/VLAMBDA**3
      RHO13SCALED   = (RHOSCALED)**(1.D0/3.D0)

      CALL VCSRLDA(ZKSRCOUL,RHO,RHO13,0.D0,MULOCAL,1,ERFEXP)

      CALL VCSRLDA(ZKSRSCALED,RHOSCALED,RHO13SCALED,
     &   0.0D0,MULOCAL,1,ERFEXP)

      ZKSR(1) = ZKSRCOUL(1) - (VLAMBDA**5) * ZKSRSCALED(1)
      ZKSR(2) = ZKSR(2) + (ZKSRCOUL(2) - (VLAMBDA**2) * ZKSRSCALED(2))
      ZKSR(3) = ZKSR(3) + (ZKSRCOUL(3) - (1/VLAMBDA)* ZKSRSCALED(3))
      RETURN
      END


C*****************************************************************************
      SUBROUTINE VCSRLDA(ZKSR,RHO,RHO13,CHI,MULOCAL,NORDER,ERFEXP)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Aug. 2003
C
C     Description :
C     Short-range LDA correlation-functional suitable for the 
C     erf and erf+exp type two-electron operators.
C     Exact mu -> infinity limit imposed
C     and fitted from RPA values with rs=0.2,0.5,1,2,3,4,5,6,7,8,9,10
C
C*****************************************************************************
#include "implicit.h"
      LOGICAL MULOCAL(0:2),ERFEXP(0:2),POTENT
      DIMENSION DENOM(3),ZKSR(3),ZVWN(3)
C
      POTENT = (NORDER.GT.0)
      CALL WVWN(ZVWN(2),RHO,RHO13,E1,POTENT,.TRUE.)
      ZVWN(1) = E1*RHO
      IF (NORDER.GT.1) CALL V1VWN(VVWN0,ZVWN(3),RHO,RHO13)
      CALL SRLDAFAC(DENOM,RHO,RHO13,ZVWN,CHI,ERFEXP,NORDER,MULOCAL)
C
      ZKSR(1)= ZVWN(1)/DENOM(1)
      IF (NORDER.GT.0) THEN
         ZKSR(2) = (ZVWN(2)-ZKSR(1)*DENOM(2))/DENOM(1)
      END IF
      IF (NORDER.GT.1) THEN
         ZKSR(3) = (ZVWN(3)-2.0D0*ZKSR(2)*DENOM(2) - 
     &              ZKSR(1)*DENOM(3))/DENOM(1)
      ENDIF
      RETURN
      END


C*****************************************************************************
      SUBROUTINE SRLDAFAC(DENOM,RHO,RHO13,ZVWN,CHI,
     &                 ERFEXP,NORDER,MULOCAL)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Sep. 2003
C
C     Purpose   : Calculate the factor used for fitting the LDA
C                 functionals to the short-range two-electron operator.
C
C*****************************************************************************
#include "implicit.h"
#include "pi.h"
#include "dftcom.h"
      PARAMETER (D0 = 0.0D0, DP25=0.25D0, DP5= 0.5D0, D1 = 1.0D0,
     &           D2 = 2.0D0, D3 = 3.0D0 , D4 = 4.0D0, D6 = 6.0D0,
     &           D9 = 9.0D0,
     &           D32= 32.0D0, THIRD = D1/D3, DP75=0.75D0,D1P5=1.5D0, 
C                RSFAC = (D3/(D4*PI))**THIRD 
     &           RSFAC = 0.62035 04908 99400 08660 D0, 
     &           D1P62 = 1.62D0)
C     Constants for g(0)
      PARAMETER (D = 3.39530545262710070631D0, A = 3.2581D0,
     &           BET = 163.44D0, GAM = 4.7125D0)
C     For ERFEXP we need :
C                CERFGAU = D1 + D6*SQRT(D3)
      PARAMETER (CERFGAU = 11.39230484541326376112)
      LOGICAL ERFEXP(0:2),MULOCAL(0:2), LMULOCAL, LMULOCAL2
      DIMENSION DENOM(3),ZVWN(3)
C
C     Operator specific parameters :
C
      IF (ERFEXP(1)) THEN
         AMU  = CHI
         U1=0.4795415558493149D0
         U2=1.0094069779513646D0
         V1=10.124678249931538D0
C        Factor for SRCMULO_A functional
         XMULFAC=4.3D0  
C        Factor for the SRCMULO_B functional fitted to the ERF operator
C        must be changed if used with ERFGAU operator
         XMULFAC2=0.814516076948D0
      ELSE IF (ERFEXP(2)) THEN
         ! we use the erf correlation functional
         AMU = CHI/D1P62 ! mu(erf) = mu(erfexp) / 1.62
         U1=1.0270741452992294D0
         U2=-0.230160617208092D0
         V1=0.6196884832404359D0
C        Factor for SRCMULO_A functional
         IF (SRCMULOFAC) THEN
            XMULFAC = XMULFAC_READIN
         ELSE
            XMULFAC=1.3D0
         END IF
C        Factor for the SRCMULO_B functional fitted to the ERF operator
C        must be changed if used with ERFEXP operator
         XMULFAC2=0.814516076948D0
      ELSE
         AMU = CHI
         U1=1.0270741452992294D0
         U2=-0.230160617208092D0
         V1=0.6196884832404359D0
C        Factor for SRCMULO_A functional
         IF (SRCMULOFAC) THEN
            XMULFAC = XMULFAC_READIN
         ELSE
            XMULFAC=1.3D0
         END IF
C        Factor for the SRCMULO_B functional
         XMULFAC2=0.814516076948D0
         IF (SRCMULOFAC) THEN
            XFAC = XMULFAC_READIN * XMULFAC2
         ELSE
            XFAC= XMULFAC2
         END IF
      ENDIF
C
            RHO2 = RHO*RHO
            RS   = RSFAC/RHO13
            RS2  = RS*RS
            DRS  = -RS/(D3*RHO)
            DRS2 = DRS*DRS
            D2RS = D4*RS/(D9*RHO2)
            LMULOCAL = MULOCAL(1)
            LMULOCAL2 = MULOCAL(2)

      IF(MULOCAL(0)) THEN
         IF(LMULOCAL) THEN
               ! SRCMULO_A
C            WRITE(LUPRI,*) 'XMULFAC ==>', XMULFAC
            AMUL = XMULFAC/RS
             print *, 'AMUL, AMU', AMUL, AMU
C            write(lupri,*) 'amul,mu==>',AMUL,AMU
            IF (AMUL .GT. AMU) THEN
               AMU = AMUL
               DAMU = -DRS*AMUL/RS
               D2AMU = -D2RS*AMUL/RS + DRS2*D2*AMUL/RS2
            ELSE
C               ... global \mu used
               LMULOCAL = .FALSE.
            END IF
         ELSE IF(LMULOCAL2) THEN
             AMUL = XFAC/SQRT(RS)
             IF (AMUL .GT. AMU) THEN
                AMU = AMUL
                DAMU = -AMUL*DRS*(D1/(D2*RS))
                D2AMU = -(D2RS*AMUL)/(D2*RS) + (D3*DRS2*AMUL)/(D4/RS2)
             ELSE
C               ... global \mu used
               LMULOCAL2 = .FALSE.
             END IF
         END IF
      END IF
C
      AMU2 = AMU*AMU
      IF (ERFEXP(1)) AMU2 = AMU2 / CERFGAU
C     ... in this case
C         C2 * AMU2 -> (C2/CERFGAU) * AMU2 = C2 * (AMU2 / CERFGAU)
C         we can thus get the right numbers by correcting AMU2
C
C     Factor for the energy
C
      GRS32 = (GAM+RS)**D1P5
      GRS12 = SQRT(GAM+RS)
      EXPGRS12 = EXP(-A*GRS12)
      G0  = D*(GRS32+BET)*EXPGRS12
      C1N = U1*RS + U2*RS2
      C1D = D1 + V1*RS
      C1  = C1N/C1D
      C2D = DP5*PI*RHO2*(G0-DP5)
      C2  = ZVWN(1)/C2D
      DENOM(1) = D1 + C1*AMU + C2*AMU2
C
C     1'st derivative of factor (used for potential)
C
      IF(NORDER.GT.0) THEN
c
      DG0  = D*D1P5*GRS12*DRS*EXPGRS12
     &     - G0*A*DP5*DRS/GRS12
      C1D2 = C1D*C1D
      DC1  = (U1 + D2*U2*RS + U2*V1*RS2)*DRS/C1D2
      DC2D = (D2*C2D)/RHO + DP5*PI*RHO2*DG0
      DC2 = (ZVWN(2) - C2*DC2D)/C2D
      DENOM(2) = DC1*AMU +DC2*AMU2
      IF (LMULOCAL) THEN
         DENOM(2) = DENOM(2) + DAMU*(C1 +C2*D2*AMU2/AMU)
      END IF
      IF (LMULOCAL2) THEN
         DENOM(2) = DENOM(2) + DAMU*(C1 +C2*D2*AMU2/AMU)
      END IF
c
      ENDIF
C
C     2'nd derivative of factor (used for hessian)
C
      IF (NORDER.GT.1) THEN
c
C   % DENOM(3)  = D2C1*AMU + D2C2*AMU2
C   % Get D2C1 : Rewrite DC1 to keep track of terms !
C   % DC1 = C1FAC*DRS/C1D2, C1FAC=(U1 + D2*U2*RS + U2*V1*RS2)
C   % So we need :
      C1FAC=(U1 + D2*U2*RS + U2*V1*RS2)
      DC1FAC = D2*U2*DRS + D2*U2*V1*RS*DRS
      C1D4   = C1D2*C1D2
      DC1D2  = D2*C1D*V1*DRS
C   % All in all
      D2C1 = ((DC1FAC*DRS + C1FAC*D2RS)*C1D2 - C1FAC*DRS*DC1D2)/C1D4
C   % Get D2C2, first derive D2G0
c     DG0  = D*D1P5*GRS12*DRS*EXPGRS12
c    &     - G0*A*DP5*DRS/GRS12
      D2G0 = D*EXPGRS12*(DRS2*DP75/GRS12
     &     +   D2RS*D1P5*GRS12-DP75*A*DRS2)
     &     + DG0*(-DP5*A*DRS/GRS12)
     &     + G0*(DP25*A*DRS2/GRS32-DP5*A*D2RS/GRS12)
      D2C2D = -D2*C2D/RHO2 + D2*DC2D/RHO + PI*RHO*DG0 + DP5*PI*RHO2*D2G0
      D2C2  = (ZVWN(3) - D2*DC2*DC2D - C2*D2C2D)/C2D
C   % And finally !
      DENOM(3)  = D2C1*AMU + D2C2*AMU2
      IF (LMULOCAL) THEN
         DENOM(3) = DENOM(3) + D2AMU*(C1 +C2*D2*AMU2/AMU)
     &                       + DAMU*DAMU*C2*D2*AMU2/(AMU*AMU)
      END IF
      IF (LMULOCAL2) THEN
         DENOM(3) = DENOM(3) + D2AMU*(C1 +C2*D2*AMU2/AMU)
     &                       + DAMU*DAMU*C2*D2*AMU2/(AMU*AMU)
      END IF
c
      ENDIF
C
      RETURN
      END


C*****************************************************************************
      SUBROUTINE EXSRLDA(ECV,RHO,RHO13,CHI,ERFEXP)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Mar. 2003
C
C     Description :
C     Short-range LDA (Dirac) exchange-functional suitable for the 
C     erf and erf+exp type short-range two-electron
C     - MNP has also implemented srLDA exchange-functional for the modified
C     erf+exp type short-range two electron
C
C     Based on implementation in Molpro by J. Toulouse and A. Savin
C
C*****************************************************************************
#include "implicit.h"
      LOGICAL ERFEXP(*)
      REAL*8  ECV, RHO, RHO13, CHI, XKSR(3)
      XKSR(:) = 0.0D0
      CALL VXSRLDA (XKSR,RHO,RHO13,CHI,0,ERFEXP)
      ECV = XKSR(1)
      RETURN
      END


C*****************************************************************************
      SUBROUTINE VXLALDA(XKSR,RHO,RHO13,VLAMBDA)
C*****************************************************************************
C Kamal SHARKAS 17-05-11
C linear complement LDA exchange-functional
C*****************************************************************************
#include "implicit.h"
      LOGICAL ERFEXP(0:2)
      REAL*8  XKSR(3)

      CALL VXSRLDA(XKSR,RHO,RHO13,0.D0,1,ERFEXP)

      XKSR(1) = (1.D0 - VLAMBDA) * XKSR(1)
      XKSR(2) = (1.D0 - VLAMBDA) * XKSR(2)
      XKSR(3) = (1.D0 - VLAMBDA) * XKSR(3)
      RETURN
      END


C*****************************************************************************
      SUBROUTINE VXSRLDA(XKSR,RHO,RHO13,CHI,NORDER,ERFEXP)
C*****************************************************************************
C
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
      PARAMETER (D0 = 0.0D0, DP25=0.25D0, DP5 = 0.5D0, D1 = 1.0D0,
     &           D2 = 2.0D0, D3 = 3.0D0, D4 = 4.0D0 , D5 = 5.0D0,
     &           D6 = 6.0D0, D7 = 7.0D0, D8 = 8.0D0 , D9 = 9.0D0,
     &           D10 = 10.0, D15=15.0D0, D16=16.0D0 , D20 = 20.0D0,
     &           D24=24.0D0, D25=25.0D0, D27 = 27.0D0, D30 = 30.D0,
     &           D40 = 40.0D0, D50 = 50.0D0, D60 = 60.0D0, D81 = 81.0D0,
     &           D96 =96.0D0, THIRD = D1/D3, SIXTH = D1/D6,
     &           FIVSIX = D5/D6, D1P5 = 1.5D0)
      LOGICAL ERFEXP(0:2)
      DIMENSION XKSR(3)
C
      CHISQ = CHI*CHI
      AMU = CHI
      COEF = - D2*CHI/SQRTPI
      COEFM = - (D10*CHI)/(D9*SQRTPI)
      EXPO = CHISQ/D3
      FACGAU  = 0.07106785043055514894D0
      DFACGAU = 16.3241942781079606553D0
      FACEXP  = 0.029428248996D0
      DFACEXP = 3.7553506259D0
      AKF = RHO13*((D3*PI2)**THIRD)
      A = AMU/(D2*AKF)
      A2 = A*A
      A3 = A2*A
      B = SQRT(EXPO)/(D2*AKF)
      B2 = B*B
      RHO2 = RHO*RHO
      EXPQ = (D3*CHISQ)/D5
      C = SQRT(EXPQ)/(D2*AKF)
      C2 = C*C
      B3 = B2*B
      C3 = C2*C
      D = D30*SQRT(D5/D3)
C
C     ... <<< Energy >>>
C        ================
C        ERF CONTRIBUTION
C        ================
         IF (A .LT. 1.D-9) THEN
            XKSR(1) = -D3/D8 * RHO*RHO13 *(D24/PI)**THIRD
         ELSE IF (A.LE.1.D+2) THEN
            XKSR(1) = - (RHO*RHO13*(D24/PI)**THIRD)
     &           *(D3/D8-A*(SQRTPI*erf(DP5/A)+
     &           (D2*A-D4*A3)*EXP(-DP25/A2)-D3*A+D4*A3))
         ELSE IF (A.LT.1.D+9) THEN
            XKSR(1) = - (RHO*RHO13*(D24/PI)**THIRD)/(D96*A2)
         ELSE
            XKSR(1) = D0
         ENDIF
         IF (ERFEXP(1)) THEN
C           ================
C           EXP CONTRIBUTION to erfgau operator
C           ================
            IF (B .LT. 1.D-9) THEN
C           ... nothing done at this point ???
            ELSE IF (B.LE.1.D+2) THEN
                XKSR(1) = XKSR(1) + RHO*
     &                D6*COEF*((PI/EXPO)**D1P5)
     &                *RHO*(D2*B3*(SQRTPI*erf(DP5/B)
     &                +(D2*B-D16*B3)*EXP(-DP25/B2)-D6*B
     &                +D16*B3))
            ELSE IF (B.LT.1.D+9) THEN
               XKSR(1) = XKSR(1) - FACGAU*RHO13*RHO/B2 ! The old form was missing RHO -MNP
            ENDIF
         ELSE IF (ERFEXP(2)) THEN

C           The modified exp contribution
            IF (C .LT. 1.D-9) THEN
C           ... nothing done at this point ???
            ELSE IF (C.LE.1.D+2) THEN
                  XKSR(1) = XKSR(1) - RHO*
     &                RHO*(D10*D10*SQRT(D15)*PI)/(D27*AMU*AMU)
     &                *(D2*C3*(SQRTPI*erf(DP5/C)
     &                  +(D2*C-D16*C3)*EXP(-DP25/C2)-D6*C
     &                  +D16*C3))
            ELSEIF (C.LT.1.D+9) THEN
               XKSR(1) = XKSR(1) - FACEXP*RHO13*RHO/C2
            ENDIF
         ENDIF
      IF (NORDER.GT.0) THEN
C     ...<<< Potential >>>
C        =================
C        ERF CONTRIBUTION
C        =================
         IF (A.LT.1.D-9) THEN
            XKSR(2) = - RHO13*(D3/PI)**THIRD
         ELSE IF (A.LE.1.D+2) THEN
            XKSR(2) = - RHO13*(D3/PI)**THIRD
     &      + D2*A*AMU/PI*(EXP(-DP25/A2)-D1)+AMU/SQRTPI
     &      * erf(DP5/A)
         ELSEIF (A.LT.1.D+9) THEN
            XKSR(2) = - PI*RHO/(D2*AMU*AMU)
         ELSE
            XKSR(2)=D0
         ENDIF
         IF (ERFEXP(1)) THEN
C           ================
C           EXP CONTRIBUTION
C           ================
c           ------- Old erf+exp
            IF (B.LT.1.D-9) THEN
C           ... nothing done at this point ???
            ELSE IF (B.LE.1.D+2) THEN
               XKSR(2) = XKSR(2) + COEF*(D2*B/SQRTPI*
     &              (EXP(-DP25/B2)-D1)+DP5*erf(DP5/B))
            ELSE IF (B.LT.1.D+9) THEN
               XKSR(2) = XKSR(2) - DFACGAU*RHO/(AMU*AMU)
            ENDIF
         ELSE IF (ERFEXP(2)) THEN
c           ------- Modified erf+exp
            IF (C.LT.1.D-9) THEN
C           ... nothing done at this point ???
            ELSE IF (C.LE.1.D+2) THEN
            XKSR(2) = XKSR(2) + COEFM*(D2*C/SQRTPI*
     &              (EXP(-DP25/C2)-D1)+erf(DP5/C)*DP5)
            ELSE IF (C.LT.1.D+9) THEN
               XKSR(2) = XKSR(2) - DFACEXP*RHO/(AMU*AMU)
            ENDIF
         ENDIF
      ENDIF
      IF (NORDER.GT.1) THEN
C     ...<<< 2'nd derivative (for Hessian) >>>
C        =================
C        ERF CONTRIBUTION
C        =================
         IF (A.LT.1.D-9) THEN
            XKSR(3) = - (THIRD/RHO)*XKSR(2)
         ELSE IF (A.LE.1.D+2) THEN
             XKSR(3) =
     &        D4*A2*(PI**(-THIRD))*((D3*RHO)**(-D2/D3))
     &       *(D1-EXP(-DP25/A2))-(PI**(-THIRD))*((D3*RHO)**(-D2/D3))
         ELSEIF (A.LT.1.D+9) THEN
            XKSR(3) = - PI/(D2*AMU*AMU)
         ELSE
            XKSR(3)=D0
         ENDIF
         IF (ERFEXP(1)) THEN
C           ================
C           EXP CONTRIBUTION
C           ================
c            ------- Old erf+exp
            IF (B.LT.1.D-9) THEN
C           ... nothing done at this point ???
            ELSE IF (B.LE.1.D+2) THEN
               XKSR(3) = XKSR(3) + D4*B2*D2*
     &                 (EXP(-DP25/B2)*(D1+(DP25/B2)) - D1)/
     &                 (((RHO2*PI)**THIRD)*(D3**SIXTH))
            ELSE IF (B.LT.1.D+9) THEN
               XKSR(3) = XKSR(3) - DFACGAU/(AMU*AMU)
            ENDIF
         ELSE IF (ERFEXP(2)) THEN
c           --------- Modified erf+exp
            IF (C.LT.1.D-9) THEN
C           ... nothing done at this point ???
            ELSE IF (C.LE.1.D+2) THEN
               XKSR(3) = XKSR(3) -
     &                ((D40*SQRT(D5)*D3**FIVSIX)/(D81*(PI*RHO2)**THIRD))
     &                *C2*(EXP(-DP25/C2)*(D1+(DP25/C2)) - D1)
            ELSE IF (C.LT.1.D+9) THEN
               XKSR(3) = XKSR(3) - DFACEXP/(AMU*AMU)
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END


C*****************************************************************************
      SUBROUTINE VWN_EPS(EPS,XARRAY,A,X0,B,C)
C*****************************************************************************
#include "implicit.h"
      XX(X)=X**2+B*X+C
      Q=SQRT(4*C-B**2)
      X=XARRAY
      EPS=A*( LOG(X**2/XX(X)) - B*(X0/XX(X0))*LOG((X-X0)**2/XX(X))
     &       +(2*B/Q)*(1-(X0*(2*X0+B)/XX(X0))) * ATAN(Q/(2*X+B)) )
      RETURN
      END


C*****************************************************************************
      SUBROUTINE VWN_DEPS(DEPS,XARRAY,A,X0,B,C)
C*****************************************************************************
#include "implicit.h"
      XX(X)=X**2+B*X+C
      Q=SQRT(4*C-B**2)
      X=XARRAY
      DEPS=A*( 2/X - (2*X+B)/XX(X) - 4*B/((2*X+B)**2+Q**2)
     & - (B*X0/XX(X0)) * (
     &  2/(X-X0) - (2*X+B)/XX(X) - 4*(2*X0+B) / ((2*X+B)**2+Q**2)
     & ) )
      RETURN
      END


C*****************************************************************************
      SUBROUTINE VXSRGGA(XKSR,RHO,RHO13,RHOGRD,RHOLAP,
     &                   RHOGHG,CHI,ERFEXP)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Sep. 2003
C
C
C     Purpose   :'Estimate' a gradient correction to the 
C                 short-range LDA exchange functional. (Dirac-exhange)
C
C     Energy    : E_x,GGA ~ A * (E_x,LDA + E_x,B)
C                 where :
C                           A = (E_x,SRLDA/E_x,LDA)
C
C     Potential : V_x.GGA ~ B * (E_x,LDA + E_x,B) 
C                         + A * (V_x,LDA + V_x,B)
C                 where :
C                           B = V_x,SRLDA/E_x,LDA - (E_X,SRLDA*V_x,LDA)/(E_x,LDA)^2
C
C*****************************************************************************
#include "implicit.h"
C
      REAL*8  XKSR(3)
      LOGICAL ERFEXP(0:2)
C     ... Get regular LDA exchange Energy and Potential (Dirac exchange)
          CALL EDRC(EXLDA,RHO,RHO13)
          CALL VDRC(VXLDA,RHO13)
C
C     ... Get Short-range LDA exchange Energy and Potential
          CALL VXSRLDA(XKSR,RHO,RHO13,CHI,1,ERFEXP)
C
C     ... Get Becke gradient correction to LDA exchange
          CALL EBCK(EXB,RHO,RHO13,RHOGRD)
          CALL VBCK(VXB,RHO,RHO13,RHOGRD,RHOLAP,RHOGHG)
C     next two lines are for debugging (should give LDA values)
c         EXB = 0.0D0
c         VXB = 0.0D0
C
C     ... Combine to get gradient corrected SR-LDA E and V
C
          A       = XKSR(1)/EXLDA
          B       = XKSR(2)/EXLDA - VXLDA*A/EXLDA
c         XKSR(1) = A * (EXLDA + EXB)
          XKSR(1) = XKSR(1) + A * EXB
c         XKSR(2) = B * (EXLDA + EXB) + A * (VXLDA + VXB)
          XKSR(2) = XKSR(2) + B * EXB + A * VXB
      RETURN
      END


C*****************************************************************************
      SUBROUTINE VLAXBEK(XKSR,RHO,RHO13,RHOGRD,RHOLAP,RHOGHG,
     &                 VLAMBDA)
C*****************************************************************************
C Kamal SHARKAS 23-05-11
C linear complement Becke exchange-functional
C*****************************************************************************
#include "implicit.h"
      REAL*8  XKSR(3)
      LOGICAL ERFEXP(0:2)

      ERFEXP(:) = .FALSE.
      CALL VXSRGGA(XKSR,RHO,RHO13,RHOGRD,RHOLAP,
     &             RHOGHG,0.D0,ERFEXP)

      XKSR(1) = (1.D0 - VLAMBDA) * XKSR(1)
      XKSR(2) = (1.D0 - VLAMBDA) * XKSR(2)
      RETURN
      END


C*****************************************************************************
      SUBROUTINE VXSRGGA2(XKSR,RHO,RHO13,RHOGRD,RHOLAP,
     &                    RHOGHG,CHI,ERFEXP)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Sep. 2003
C
C
C     Purpose   :'Estimate' a gradient correction to the 
C                 short-range LDA exchange functional. (Dirac-exhange)
C
C     Energy    : E_x,GGA ~  A * (E_x,LDA + E_x,B)
C     Potential : V_x.GGA ~ dA * (E_x,LDA + E_x,B) 
C                         +  A * (V_x,LDA + V_x,B)
C                 Where we make the assumption that A is the same factor
C                 as used for the SR-LDA cor. func. (E_sr,lda ~ A * E_lda)
C
C*****************************************************************************
#include "implicit.h"
C
      LOGICAL ERFEXP(0:2), MULOCAL(0:2)
      REAL*8  DENOM(3),ZVWN(3),XKSR(3)
C     ... Get regular LDA exchange Energy and Potential (Dirac exchange)
          CALL EDRC(EXLDA,RHO,RHO13)
          CALL VDRC(VXLDA,RHO13)
C
C     ... Get Becke gradient correction to LDA exchange
          CALL EBCK(EXB,RHO,RHO13,RHOGRD)
          CALL VBCK(VXB,RHO,RHO13,RHOGRD,RHOLAP,RHOGHG)
C
C     ... Get Short-range factors
          CALL WVWN(ZVWN(2),RHO,RHO13,ECFAC,.TRUE.,.TRUE.)
          ZVWN(1) = ECFAC*RHO
          MULOCAL(:) = .FALSE.
          CALL SRLDAFAC(DENOM,RHO,RHO13,ZVWN,CHI,ERFEXP,2,MULOCAL)
C
C     ... Combine to get gradient corrected SR-LDA E and V
C
          EX = EXLDA + EXB
          XKSR(1) = EX/DENOM(1)
c
          VX = VXLDA + VXB
          XKSR(2) = (VX-XKSR(1)*DENOM(2))/DENOM(1)
      RETURN
      END


C*****************************************************************************
      SUBROUTINE VCSRGGA(ZKSR,RHO,RHO13,RHOGRD,RHOLAP,CHI,ERFEXP)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Sep. 2003
C
C     Purpose   :'Estimate' a gradient correction to the 
C                 short-range LDA correlation functional. (LYP)
C
C     Energy    : E_c,GGA ~ A * E_lyp
C     Potential : V_c.GGA ~ A * V_lyp + dA * E_lyp
C                 Where we make the assumption that A is the same factor
C                 as used for the SR-LDA func. (E_sr,lda ~ A * E_lda)
C
C*****************************************************************************
#include "implicit.h"
      LOGICAL ERFEXP(0:2),  MULOCAL(0:2)
      REAL*8  ZKSR(3),DENOM(3),ZVWN(3)
C
      CALL ELYP(ECLYP,RHO,RHO13,RHOGRD)
      CALL VLYP(VCLYP,RHO,RHO13,RHOGRD,RHOLAP)
      CALL WVWN(ZVWN(2),RHO,RHO13,ECFAC,.TRUE.,.TRUE.)
      ZVWN(1) = ECFAC*RHO
      MULOCAL(0:2) = .FALSE.
      CALL SRLDAFAC(DENOM,RHO,RHO13,ZVWN,CHI,ERFEXP,2,MULOCAL)
C     next two lines are for debugging (should give LDA values)
c     ECLYP = ECFAC
c     VCLYP = VCFAC
C
C     ===============
C     ENERGY
C     ===============
C
      ZKSR(1) = ECLYP/DENOM(1)
C
C     ===============
C     POTENTIAL
C     ===============
C
      ZKSR(2) = (VCLYP-ZKSR(1)*DENOM(2))/DENOM(1)
C
      RETURN
      END


C*****************************************************************************
      SUBROUTINE VNCLALYP(ZKSR,RHO,RHO13,RHOGRD,RHOLAP,VLAMBDA)
C*****************************************************************************
C Kamal SHARKAS 23-05-11
C linear Non-scaled comlement LYP functional  
C*****************************************************************************
#include "implicit.h"
      REAL*8  ZKSR(3),DENOM(3),ZVWN(3)
      LOGICAL ERFEXP(0:2)

      ERFEXP(:) = .FALSE.
      CALL VCSRGGA(ZKSR,RHO,RHO13,RHOGRD,RHOLAP,
     &             0.0D0,ERFEXP)
      ZKSR(1) = (1.D0 - VLAMBDA**2) * ZKSR(1)
      ZKSR(2) = (1.D0 - VLAMBDA**2) * ZKSR(2)
      RETURN
      END


C*****************************************************************************
      SUBROUTINE VSCLALYP(ZKSR,RHO,RHO13,RHOGRD,RHOLAP,VLAMBDA)
C*****************************************************************************
C Scaled complement LYP functional
C*****************************************************************************
#include "implicit.h"
      REAL*8  RHOSCA,RHO13SCA,RHOGRDSCA,RHOLAPSCA
      REAL*8  ZKSR(3),ZKSRCOUL(3),ZKSRSCA(3)
      LOGICAL ERFEXP(0:2)

      ZKSRCOUL(1:3) = 0.D0
      ERFEXP(0:2)   = .FALSE.

      CALL VCSRGGA(ZKSRCOUL,RHO,RHO13,RHOGRD,RHOLAP,
     &             0.0D0,ERFEXP)

      RHOSCA   = RHO/VLAMBDA**3
      RHO13SCA = (RHOSCA)**(1.D0/3.D0)
      RHOGRDSCA= RHOGRD/VLAMBDA**4
      RHOLAPSCA= RHOLAP/VLAMBDA**5

      ZKSRSCA(1) = 0.D0
      ZKSRSCA(2) = 0.D0
      CALL VCSRGGA(ZKSRSCA,RHOSCA,RHO13SCA,RHOGRDSCA,RHOLAPSCA,
     &             0.0D0,ERFEXP)

      ZKSR(1) = ZKSRCOUL(1) - (VLAMBDA**5) * ZKSRSCA(1)
      ZKSR(2) = ZKSR(2) + (ZKSRCOUL(2) - (VLAMBDA**2) * ZKSRSCA(2))
      RETURN
      END


C*****************************************************************************
      SUBROUTINE VSRWINT(ZKSR,ZMU,RHO,RHO13,CHI,NORDER,ERFEXP)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Sep. 2004
C
C     Purpose   : Generate short-range functional as a weighted interpolation
C                 between a regular (Kohn-Sham) functional in the \mu=0
C                 limit and the srlda functional.
C                 (As suggested by A. Savin and J. Toulouse)
C
C     Energy    : E ~ (E_KS - E_LDA)erfc(Rs*\mu) + E_srlda,mu
C     Potential : V ~ (V_KS - V_LDA)erfc(Rs*\mu) + 
C                     (E_KS - E_LDA)(\mu*(dRs/drho)(erfc(Rs*\mu)/drho)) + 
C                      V_srlda,mu
C
C     NOTE : On entry ZKSR contains the energy and potential
C            from the (E_KS - E_LDA) and (V_KS - V_LDA) terms.
C
C*****************************************************************************
#include "implicit.h"
#include "pi.h"
#include "dftcom.h"
      LOGICAL ERFEXP(0:2)
      DIMENSION ZKSR(3),ZMU(3)
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, D4 = 4.0D0,
     &           D9 = 9.0D0,
     &          RSFAC = 0.62035 04908 99400 08660 D0)
C               RSFAC = (D3/(D4*PI))**THIRD
c
      EDIFF  = ZKSR(1)
      VDIFF  = ZKSR(2)
      DVDIFF = ZKSR(3)
c
      RS   = RSFAC/RHO13
      DRS  = -RS/(D3*RHO)
      CHIM = CHI*COPFAC
      X  = CHIM*RS
      DERFCX = ERFC(X)
      DDERFCX = -(D2/SQRTPI)*EXP(-X**2)
C
c
      ZKSR(1) = EDIFF*DERFCX + ZMU(1)
      ZKSR(2) = VDIFF*DERFCX + EDIFF*DDERFCX*CHIM*DRS + ZMU(2)
      IF (NORDER.GT.1) THEN
         D2RS = D4*RS/(D9*RHO*RHO)
         D2DERFCX = -D2*X*DDERFCX
         ZKSR(3) = DVDIFF*DERFCX + D2*VDIFF*DDERFCX*CHIM*DRS +
     &          EDIFF*(D2DERFCX*(CHIM*DRS)**2 + DDERFCX*CHIM*D2RS) +
     &          ZMU(3)
      END IF
      RETURN
      END


C*****************************************************************************
      SUBROUTINE VSRWINT2(ZKSR,ZMU,RHO,RHO13,CHI,NORDER,ERFEXP)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Sep. 2004
C
C     Purpose   : Generate short-range functional as a combination of a 
C                 (gradient corrected) functional for small \mu
C                 and the usual short-range LDA for large \mu. The
C                 combination is controller by the error-function.
C
C     Energy    : E ~ E_KS*erfc(Rs*mu/c) + E_srlda,mu*erf(Rs*mu/c)
C     Potential : V ~ V_KS*erfc(Rs*mu/c) + E_KS*erfc'(Rs*mu/c) +
C                     V_srlda,mu*erf(Rs*mu/c) + E_srlda,mu*erf'(Rs*mu/c)
C     Hessian   : V'~ V'_KS*erfc(Rs*mu/c) + 2*V_KS*erfc'(Rs*mu/c) +
C                     E_KS,mu*erfc''(Rs*mu/c) + 
C                     V'_srlda,mu*erf(Rs*mu/c) + 2*V_srlda,mu*erf'(Rs*mu/c) + 
C                     E_srlda,mu*erf''(Rs*mu/c)
C
C*****************************************************************************
#include "implicit.h"
#include "pi.h"
#include "dftcom.h"
      LOGICAL ERFEXP(0:2)
      DIMENSION ZKSR(3),ZMU(3)
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, D4 = 4.0D0,
     &           D9 = 9.0D0, D13= 1.0D0/3.0D0,D23=2.0D0/3.0D0,
     &          RSFAC = 0.62035 04908 99400 08660 D0)
C               RSFAC = (D3/(D4*PI))**THIRD
c
      EKS  = ZKSR(1)
      VKS  = ZKSR(2)
      DVKS = ZKSR(3)
      ESR  = ZMU(1)
      VSR  = ZMU(2)
      DVSR = ZMU(3)
      RS   = RSFAC/RHO13
      RS2  = RS*RS
      DRS  = -RS/(D3*RHO)
      DRS2 = DRS*DRS
      CHIM = CHI*COPFAC
      X  = CHIM*RS
      DERFX = erf(X)
      DDERFX = (D2/SQRTPI)*EXP(-X**2)
c 
      IF (IWINT.EQ.3) THEN
         GAUX=(D2/SQRTPI)*X*EXP(-D13*X**2)
         GFAC=(D1/RS)*DRS - D23*CHIM*X*DRS
         DGAUX=GFAC*GAUX
      END IF
c
      DERFCX = ERFC(X)
      DDERFCX = -(D2/SQRTPI)*EXP(-X**2)
c
c    --- Energy
      ZKSR(1) = EKS*DERFCX + ESR*DERFX
c     ...Exponential term for the erfgau operator
      IF (IWINT.EQ.3) ZKSR(1) = ZKSR(1) + EKS*GAUX - ESR*GAUX
c
c    --- First derivative
      ZKSR(2) = VKS*DERFCX + EKS*DDERFCX*CHIM*DRS +
     &          VSR*DERFX  + ESR*DDERFX*CHIM*DRS
      IF (IWINT.EQ.3) THEN
c     ...Exponential term for the erfgau operator
         ZKSR(2) = ZKSR(2)  + VKS*GAUX + EKS*DGAUX -
     &             VSR*GAUX - ESR*DGAUX
      ENDIF
c    --- Second derivative
      IF (NORDER.GT.1) THEN
         D2RS = D4*RS/(D9*RHO*RHO)
         D2DERFCX = -D2*X*DDERFCX
         D2DERFX  = -D2*X*DDERFX
         ZKSR(3)  = DVKS*DERFCX + D2*VKS*DDERFCX*CHIM*DRS +
     &       EKS*(D2DERFCX*(CHIM*DRS)**2 + DDERFCX*CHIM*D2RS) +
     &       DVSR*DERFX + D2*VSR*DDERFX*CHIM*DRS +
     &       ESR*(D2DERFX*(CHIM*DRS)**2 + DDERFX*CHIM*D2RS) 
c        ...Exponential term for the erfgau operator
         IF (IWINT.EQ.3) THEN
            DGFAC = ((-D1/RS2) - D23*CHIM**2)*DRS2 
     &            + (D1/RS - D23*CHIM*X)*D2RS
            D2GAUX=DGFAC*GAUX + GFAC*DGAUX
            ZKSR(3) = ZKSR(3) + DVKS*GAUX + D2*VKS*DGAUX + EKS*D2GAUX -
     &                          DVSR*GAUX - D2*VSR*DGAUX - ESR*D2GAUX
         END IF
      END IF
      RETURN
      END


C*****************************************************************************
      SUBROUTINE CONDFT
C*****************************************************************************
C     DFT CONSTANTS
C*****************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
#include "pi.h"
      logical vwn5
      COMMON/DFTCON/BCX,CVX,AL,BL,CL,DL,CF,C3,CRS,PRE,THIRD,CINF,FF,AA
     &,BB,GG,DD,AK1,AK2,AK3,PP,FT,QQ,BETAP,CV,AI,BI,CI,X0I,AII,BII,CII
     &,X0II,AIII,BIII,CIII,X0III,FPPZ,FTT,FTTT,QQQI,QQQII,QQQIII,XF0I
     &,YF0I,XF0II,YF0II,XF0III,YF0III,DCRS,PD,Q,P,ALL,BLL,CLL,DLL,CFFF
     &,ELL,FLL,GLL,HLL,PLL,CFFFF,ONE,TWO,THREE,FOUR,FIVE,EIGHT,TEN
     &,ELEVEN,VAR6,POW1,VAR14,VAR2,CFG,C2G,PRODAB,C3G,C33G,C4G,C44G
     &,C5G,C55G,C6G,C66G,C7G,C77G,C8G,C88G,C9G,C99G,C11G,C1111G
C
      THIRD=1.0D0/3.0D0
      PRE=(2.0D0**THIRD)
      BCX=-PRE*0.75D0*((3.0D0/PI)**THIRD)
      CVX=-PRE*((3.0D0/PI)**THIRD)
C
C
      ONE=1.0D0
      TWO=2.0D0
      THREE=3.0D0
      FOUR=4.0D0
      FIVE=5.0D0
      EIGHT=8.0D0
      TEN=10.0D0
      ELEVEN=11.0D0
C
      POW1=TWO/THREE
      ALL=0.04918D0
      BLL=0.132D0
      CLL=0.2533D0
      DLL=0.349D0
      ELL=(47.0D0/18.0D0)
      FLL=(7.0D0/18.0D0)
      GLL=(5.0D0/2.0D0)
      HLL=(1.0D0/18.0D0)
      PLL=(5.0D0/18.0D0)
      C333=2.0D0**(11.0D0/3.0D0)
      CFFF=0.3D0*((3.0D0*PI*PI)**(2.0D0/3.0D0))*C333
      CFFFF=0.3D0*((3.0D0*PI*PI)**(2.0D0/3.0D0))
      CRS=(3.0D0/(4.0D0*PI))**THIRD
      DCRS=sqrt(CRS)
      PRODAB=ALL*BLL
      C3G=PRODAB*ELL
      C33G=TWO*C3G
      C4G=PRODAB*FLL
      C44G=TWO*C4G
      C5G=GLL*PRODAB
      C55G=FIVE*PRODAB
      C6G=PRODAB*HLL
      C66G=PRODAB/9.0D0
      C7G=PRODAB/9.0D0
      C77G=TWO*PRODAB/9.0D0
      C8G=ELEVEN*PRODAB/9.0D0
      C88G=TWO*ELEVEN*PRODAB/9.0D0
      C9G=PRODAB*POW1
      C99G=C9G*TWO
      C11G=PRODAB
      C1111G=TWO*PRODAB
      VAR14=(ELEVEN/THREE)
      VAR2=TWO**VAR14
      CFG=(THREE/TEN)*(THREE*PI*PI)**POW1
      C2G=PRODAB*CFG*VAR2
      VAR6=EIGHT/THREE
C
C     CONSTANTS USED IN SP86
      CINF=0.004235D0
      FF=1.745D0*0.11D0
      AA=0.023266D0
      BB=7.389D-6
      GG=8.723D0
      DD=0.472D0
      AK1=0.001667D0
      AK2=0.002568D0
      AK3=10000D0
      PP=5.0D0/3.0D0
      FT=4.0D0/3.0D0
      QQ=2.0D0**(-5.0D0/6.0D0)
C
C     BECKE FUNCTIONAL CONSTANTS
      BETAP=0.0042D0
C
C     VWN CONSTANTS
c     CV=0.5D0
      CV=0.5D0

C
      vwn5=.true.
      if(vwn5)then
      AI=0.0621814D0
      BI=3.72744D0
      CI=12.9352D0
      X0I=-0.10498D0
      else
      ai=0.0621814D0
      bi=13.0720d0
      ci=42.7198d0
      x0i=-0.409286d0
      endif

      QQQI=sqrt(4.0D0*CI-BI*BI)
      XF0I=X0I*X0I+BI*X0I+CI
      YF0I=QQQI/(2.0D0*X0I+BI)
C
      AII=0.0310907D0
      BII=7.06042D0
      CII=18.0578D0
      X0II=-0.32500D0
      QQQII=sqrt(4.0D0*CII-BII*BII)
      XF0II=X0II*X0II+BII*X0II+CII
      YF0II=QQQII/(2.0D0*X0II+BII)
C
      AIII=-0.033774D0
      BIII=1.131071D0
      CIII=13.0045D0
      X0III=-0.0047584D0
      QQQIII=sqrt(4.0D0*CIII-BIII*BIII)
      XF0III=X0III*X0III+BIII*X0III+CIII
      YF0III=QQQIII/(2.0D0*X0III+BIII)
C
      FPPZ=8.0D0/(9.0D0*((2.0D0**FT)-2.0D0))
      FTT=1.0D0/3.0D0
      FTTT=(2.0D0**FT)-2.0D0
C     CONSTANTS FOR WOP86
      PD=2.0D0/3.0D0
      Q=-1.5D0
      P=5.0D0/3.0D0
      RETURN
      END
