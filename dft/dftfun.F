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
C     **************************
C     ***** Dirac exchange *****
C     **************************
C
C  /* Deck edrc */
      SUBROUTINE EDRC(ED,RHO,RHO13)
#include <implicit.h>
#include <pi.h>
      PARAMETER (D1 = 1.0D0, D3 = 3.0D0, D4 = 4.0D0)
c     &           CX = -(D3/D4)*((D3/PI)**(D1/D3)))
      ED = -(D3/D4)*((D3/PI)**(D1/D3))*RHO*RHO13
      RETURN
      END
C
C  /* Deck vdrc */
      SUBROUTINE VDRC(VD,RHO13)
#include <implicit.h>
#include <pi.h>
      PARAMETER (D1 = 1.0D0, D3 = 3.0D0)
c      , CXV = -(D3/PI)**(D1/D3))
      VD = -((D3/PI)**(D1/D3))*RHO13
      RETURN
      END
C
C  /* Deck V1DRC */ 
      SUBROUTINE V1DRC(VD0,VD1,RHO,RHO13)
#include <implicit.h>
#include <pi.h>
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, 
     &           D13 = D1/D3)
c, CXV = -(D3/PI)**D13)
C
      VD0 = (-(D3/PI)**D13)*RHO13
      VD1 = (-(D3/PI)**D13)*D13*RHO13/RHO
      RETURN
      END
C
C     ********************************
C     ***** Vosko, Wilk & Nusair ***** 
C     ********************************
C
C  /* Deck evwn */
      SUBROUTINE EVWN(ECV,RHO,RHO13)
C
C     Correlation energy (G. J. Laming)
C
#include <implicit.h>
C
      CALL WVWN(VA1,RHO,RHO13,E1,.FALSE.,.TRUE.,.FALSE.)
      ECV = E1*RHO
      RETURN
      END
C
C  /* Deck vvwn */
      SUBROUTINE VVWN(VA,RHO,RHO13)
#include <implicit.h>
      CALL WVWN(VA,RHO,RHO13,E1,.TRUE.,.TRUE.,.FALSE.)
      RETURN
      END
C
C
C  /* Deck wvwn */
      SUBROUTINE WVWN(VA1,RHO,RHO13,E1,POTENT,SPNCMP,SPNINT)
C
C     C. W. MURRAY & G. J. LAMING
C
#include <implicit.h> 
#include <pi.h>
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, D4 = 4.0D0,
     &           D6 = 6.0D0)
      PARAMETER (CV = D1/D2) 
c      PARAMETER (DCRS = (D3/(D4*PI))**(D1/D6))
C    &           AI   = 0.0621814D0,
C    &           BI   = 13.0720D0,
C    &           CI   = 42.7198D0,
C    &           X0I  = -0.409286D0,
C    &           AI=0.0621814D0,
C    &           BI=3.72744D0,
C    &           CI=12.9352D0,
C    &           X0I=-0.10498D0,
C    &           Q    = DSQRT(D4*CI - BI**2),
C    &           XF0I = X0I*X0I + BI*X0I + CI,
C    &           YF0I = Q/(BI + 2*X0I))
      LOGICAL POTENT,SPNCMP, SPNINT
#include <dftcom.h>
      DCRS = (D3/(D4*PI))**(D1/D6)
C
      IF (SPNINT) THEN
         AI  = -0.0337737D0
         BI  =  1.13107D0
         CI  = 13.0045D0
         X0I = -0.0047584D0
      ELSE 
C
C        Gaussian version
C
         IF (DOVWN3) THEN
            IF (SPNCMP) THEN
               AI   =  0.0621814D0
               BI   = 13.0720D0
               CI   = 42.7198D0
               X0I  = -0.409286D0
            ELSE
               AI  =   0.0310907D0
               BI  =  20.1231D0 
               CI  = 101.578D0
               X0I =  -0.743294D0
            END IF
C
C        standard version
C
         ELSE
            IF (SPNCMP) THEN
               AI  =  0.0621814D0
               BI  =  3.72744D0
               CI  = 12.9352D0
               X0I = -0.10498D0
            ELSE
               AI  =  0.0310907D0
               BI  =  7.06042D0 
               CI  = 18.0578D0
               X0I = -0.32500D0
            END IF
         END IF
      END IF
      Q    = DSQRT(D4*CI - BI**2)
      XF0I = X0I*X0I + BI*X0I + CI
      YF0I = Q/(BI + 2*X0I)
C
      X   = DCRS/DSQRT(RHO13)
      B   = X0I/XF0I
      C   = XF0I*YF0I
      XF  = X*X + BI*X + CI
      XFX = D2*X + BI
      YF  = Q/XFX
      YFX = -D2*YF/XFX
      E1  = AI*CV*(D2*DLOG(X) + DLOG(XF)*(BI*B - D1)
     &            +D2*BI*(D1/Q - X0I/C)*DATAN(YF)
     &            -D2*(BI*B)*DLOG(X - X0I))
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
C
C  /* Deck V1VWN */ 
      SUBROUTINE V1VWN(VD0,VD1,RHO,RHO13)
C
C     WRITTEN   By    G. J. Laming
C     Some adjustments by A. M. Lee
C
#include <implicit.h>
#include <pi.h>
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

      Q    = DSQRT(D4*CI - BI**2)
      XF0I = X0I*X0I + BI*X0I + CI
      YF0I = Q/(BI + 2*X0I)
      DCRS = (D3/(D4*PI))**(D1/D6)
C
      X    = DCRS/DSQRT(RHO13)
      B    = X0I/XF0I
      C    = XF0I*YF0I
      XF   = X*X + BI*X+CI
      XFX  = D2*X + BI
      YF   = Q/XFX
      ACON = B*BI - D1
      BCON = D2*ACON + D2
      CCON = D2*BI*(D1/Q - X0I/C)
C
      E1   = D2*DLOG(X) 
     &     + ACON*DLOG(XF)
     &     - BCON*DLOG(X - X0I)
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
C
C  /* Deck vtvwn */
      SUBROUTINE VTVWN(VD1,RHO,RHO13)
#include <implicit.h>
      PARAMETER(D1 = 1.0D0, D2 = 2.0D0, DP3 = 1.0D0/3.0D0, 
     &          D8D9 = 8.0D0/9.0D0)
#include <dftcom.h>
      IF (DOVWNI) THEN
         FAC = -D8D9/(D2**DP3 - D1)
         CALL WVWN(VA,RHO,RHO13,E0,.FALSE.,.TRUE.,.FALSE.)
         CALL WVWN(VA,RHO,RHO13,E1,.FALSE.,.FALSE.,.FALSE.)
         VD1 = FAC*(E0 - E1)/RHO
      ELSE
         CALL WVWN(VA,RHO,RHO13,EI,.FALSE.,.FALSE.,.TRUE.)
         VD1 = D2*EI/RHO
      END IF
      RETURN
      END
C     *****************
C
C /* Deck ebck */
      SUBROUTINE EBCK(BECKI,RHO,RHO13,RHOGRD)
C     *************WRITTEN BY G. J. LAMING********************
#include <implicit.h>
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0)
      PARAMETER (BETAP = 0.0042D0)
      DF = (D1/D2)**(D1/D3)
C
      XALPHA = RHOGRD/(DF*RHO*RHO13)
      SINHIA = DLOG(XALPHA+DSQRT((XALPHA*XALPHA) + D1))
      DENOMA = (D1 + (6.0D0*BETAP*XALPHA*SINHIA))
      BECKI  = -BETAP*((DF*RHO*RHO13)*XALPHA*XALPHA/DENOMA)
      RETURN
      END
C /* Deck gbck */ 
      SUBROUTINE GBCK(VB1,VB2,RHO,RHO13,RHOGRD)
#include <implicit.h>
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, DP5 = 0.5D0, D6 = 6.0D0)
      PARAMETER (FT = 4.0D0/3.0D0, BETAP = 0.0042D0)
      DF = D2**FT
C
C     AMMENDED FEBRUARY 1994 G. J. Laming
C
      XL    = DF/(RHO*RHO13)
      AAA   = FT*RHOGRD/RHO
      ALPHA = DP5*RHOGRD*XL
      EXV   = DSQRT(D1 + ALPHA*ALPHA)
      ASH   = DLOG(ALPHA + EXV)
      B     = BETAP/(D1 + D6*BETAP*ALPHA*ASH)
      BBB   = D6*ALPHA*B*(ASH+ALPHA/EXV)
      VB2   = (BBB - D2)*B*XL
      VB1   = (D1 - BBB)*ALPHA*B*AAA
      RETURN
      END
c /* Deck vbck */
      SUBROUTINE VBCK(VA,RHO,RHO13,RHOGRD,RHOLAP,RHOGHG)
#include <implicit.h>
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, D4 = 4.0D0, 
     &           D6 = 6.0D0, D12 = 12.0D0, DP5 = 0.5D0)
      PARAMETER (BETAP=0.0042D0)
      DF = D2**(D1/D3)
C
      XL    =  DF*D2/(RHO*RHO13)
      XK    =  DP5*RHOGRD
      G     =  XL*XL*XK*XK
      EXV   =  DSQRT(D1+G)
      EXV2  =  D6*BETAP*XK*XL
      ASH   =  DLOG(XK*XL+EXV)
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
C  /* Deck B88xd2 */
      SUBROUTINE V1BCK(DF1000,DF0010,DF2000,DF1010,DF0020,
     &                 RHOINP,RHOGRD)
C
#include <implicit.h>
C
      PARAMETER (D0 = 0.0D0, DP5 = 0.50D0, D1 = 1.0D0, D2 = 2.0D0, 
     &           D3 = 3.0D0, D4 = 4.0D0, D43 = D4/D3, D6 = 6.0D0, 
     &           D7 = 7.0D0, BETAP = 0.0042D0)
C
C--------------------------------------------------------------------
C     (C) A. M. Lee, Cambridge 1993.
C
C     Calculates the second and lower derivatives of the Becke '88
C     exchange functional, wrt the density and related functions.
C     work with alpha-spin density, RHO; and length of the
C     alpha-spin gradient vector, GRHO. 
C
C--------------------------------------------------------------------
C
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
         DENX1 = DSQRT(DENX2)
         DAS   = DLOG(X + DSQRT(D1 + X*X))
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
C
C     ****************************
C     ***** Lee, Yang & Parr *****
C     ****************************
C
C /* Deck elyp */
      SUBROUTINE ELYP(CLYPCI,RHO,RHO13,RHOGRD)
C     *************WRITTEN BY G. J. LAMING********************
#include <implicit.h>
#include <pi.h>
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
      W     = (RHO13/RRR)*DEXP(-CLL*TM)/(RRR*AU)
      DT    = CLL*TM+DLL*TM/AU
      C2    = ALL*BLL*W*RHOA*RHOA
      C1    = ALL*BLL*W
      C5    = (AMASQ+AMASQ)
      C6    = RHOA*AMASQ/RHO
      C7    = RHOA*AMASQ/RHO
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
C /* Deck glyp */
      SUBROUTINE GLYP(TERM1,U,RHO,RHO13,RHOGRD)
C
C***********WRITTEN   BY   G. J. LAMING******************
C
#include <implicit.h>
#include <pi.h>
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
C
      RRR   = RHO*RHO
      BG    = ALL*BLL
      RRRR  = RRR*RHO
      XXX   = BG*CFFFF
      RHOT  = 1.0D0/RHO13
      XB    = 1.0D0+DLL*RHOT
      AG    = RHOT/(3.0D0*RHO)
      XB2   = XB*XB
      XBRHO = -DLL*AG
      XD    = DEXP(-CLL*RHOT)
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
c /* Deck vlyp */
      SUBROUTINE VLYP(VA,RHO,RHO13,RHOGRD,RHOLAP)
C
C     Lee
C  
#include <implicit.h>
#include <pi.h>
      PARAMETER (ALL=0.04918D0,BLL=0.132D0,CLL=0.2533D0,DLL=0.349D0,
     &           FLL=(7.0D0/18.0D0),PLL=(5.0D0/18.0D0))
      CFFFF=0.3D0*((3.0D0*PI*PI)**(2.0D0/3.0D0))
C
      RRR    = RHO*RHO
      RRRR   = RHO*RRR
      RHOT   = 1.0D0/RHO13
      XB     = 1.0D0+DLL*RHOT
      AG     = RHOT/(3.0D0*RHO)
      BG     = ALL*BLL
      BG2    = BG*CFFFF
      XB2    = XB*XB
      XBRHO  = -DLL*AG
      XD     = DEXP(-CLL*RHOT)
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
C  /* Deck v1lyp */
      SUBROUTINE V1LYP(DF1000,DF0010,DF2000,DF1010,DF0020,RHOINP,RHOGRD)
#include <implicit.h>
#include <pi.h>
      PARAMETER (DP5=0.050,
     &           ALL=0.04918D0,BLL=0.132D0,CLL=0.2533D0,DLL=0.349D0)
C
C-------------------------------------------------------------------
C     (C) A. M. Lee Cambridge 1993/94
C
C     The second and lower derivatives of the LYP functional,wrt
C     the density and related quantities.
C 
C     The following is closed shell code and cannot be generalised.
C     Care has been taken to avoid numerical instabilities.
C
C-------------------------------------------------------------------
C
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
      DS2=DSQRT(2.0D0)
C
      rhoa=0.50D0*RHOINP
      IF (DABS(RHOA).LT.1D-10 .OR. DABS(RHOGRD).LT.1D-10) THEN
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
      epsilon=eps
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
      omega1=DEXP(-crho13)
      omega2=rho23*rhon2*rhon1
      omega=eps*omega1*omega2

C     Examine carefully:

      Dw10=(c13*crho43+Deps10/eps
     1-11.0D0/3.0D0*rhon1)

      Dw10x2=omega*Dw10*Dw10
      Dw10=omega*Dw10

C    Dw10x2=Dw10*dw10/omega is also unstable:
C    lim rho->big, dw10/omega->0/0

      alt=c13*drho43*(deps10+X10*eps)
      Dw20=Dw10x2+omega*(
     1X10*c13*crho43+11.0D0/3.0D0*rhon2+alt)



      delta=drho13*eps+crho13
      Ddel10=-c13*drho43*eps-c13*crho43+
     1drho13*Deps10
      Ddel20=c49*rhon1*rho43*(d*eps+c)
     1-2.0D0*c13*drho43*Deps10+drho13*Deps20

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
     1Dw10*(G3+G2*rhoa2)
      DF10=-(a*DF110+ab*DF10)

      DF11=
     1    omega*(DG311+rhoa*(DG201+rhoa*DG211))+
     1Dw10*(DG301+DG201*rhoa2)
      DF11=-(ab*DF11)

      DF20=DG320*omega+DG310*Dw10+DG310*Dw10+G3*Dw20+
     1DG210*rhoa2*dw10+G2*rhoa*dw10+dw20*rhoa2*G2+
     1(G2+rhoa*DG210)*(0.50D0*omega+rhoa*Dw10)+
     1rhoa*omega*(DG210+0.50D0*DG210+rhoa*DG220)
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
      subroutine vtlyp(F20000,F02000,F11000,F10100,F01010,
     1   F10010,F01100,F10001,F01001,pa,pb,gaa,gbb,gab)

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

      w=dexp(-c*q)*rp11*hm1
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
     1       +s3*(gaambb)*(-2.0d0*y/p/p))

c Timings only
      bsdada=-ab*(Co*40/9*pa23+j3+s3d*2*y/p*(gaambb)
     1       +s3*(gaambb)*(-2.0d0*y/p/p))

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
     1       +w*(2.0D0*pb*bsda+st*bsdada+btdada)

c For timings only

      F02000=jaa+wdd*(st*bs+bt)+2.0d0*wd*(pb*bs+st*bsda+btda)
     1       +w*(2.0D0*pb*bsda+st*bsdada+btdada)

c End
 
      F11000=jab+wdd*(st*bs+bt)
     1       +wd*(pa*bs+st*bsdb+btdb+pb*bs+st*bsda+btda)
     2       +w*(bs+pb*bsdb+pa*bsda+st*bsdadb+btdadb)

c    ****Second derivatives involving gab and a density

      F10001=-abw*2.0d0*(pb*s1+st*s1d-4.0d0/3.0d0*p)
     1       -ab*wd*2.0d0*(st*s1-2.0d0/3.0d0*p2)
      F01001=-abw*2.0d0*(pa*s1+st*s1d-4.0d0/3.0d0*p)
     1       -ab*wd*2.0d0*(st*s1-2.0d0/3.0d0*p2)


c    ****All other derivatives are aero for Lyp

c         F00200=F00020=F00110=F00101=F00011=F00002=0
c     write(6,*)'w,wd,wdd',w,wd,wdd
c     write(6,*)'del,deld,deldd',del,deld,deldd
c     write(6,*)'bsda,bsdb',bsda,bsdb
c     write(6,*)'jaa,jab',jaa,jab
c     write(6,*)'pagbb,pbgbb',F10010,F01010
c     write(6,*)'papa',F20000
c     write(6,*)'papb',F11000
c     write(6,*)'pagab',F10001
c     write(6,*)'pagaa,bs',F10100,bs
c     write(6,*)'btda,btdb,btdadb,btdada',btda,btdb,btdadb,btdada
c     write(6,*)'pbgaa,ta,tb,abwb2',F01100,ta,tb,ab*wd*pb*pb
c     write(6,*)'F10001,bt,qwewdab',F10001,bt,wd*qwe*ab
c     write(6,*)'bsdada,bsdadb,F01001',bsdada,bsdadb,F01001

c     write(6,*)'DF1010,DF2000',
c    1           F10100,F20000
      return
      end
      SUBROUTINE GLYPCO(EASYGA,EASYGB,HARDGA,HARDGB,HARDER,
     &                  RHO,RHO13,RHOGRD,UHF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (DP5 = 0.50D0)
C
C
C     OPENSHELL LYP CORRELATION  GRADIENT CODE
C
C***********WRITTEN   BY   G. J. LAMING******************
C
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
C
      SUBROUTINE GCSLYP(TERM1,U,RHO,RHO13,GMOD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DFTCON/BCX,CVX,AL,BL,CL,DL,CF,C3,CRS,PRE,THIRD,CINF,FF,AA
     1,BB,GG,DD,AK1,AK2,AK3,PP,FT,QQ,BETAP,CV,AI,BI,CI,X0I,AII,BII,CII
     1,X0II,AIII,BIII,CIII,X0III,FPPZ,FTT,FTTT,QQQI,QQQII,QQQIII,XF0I
     1,YF0I,XF0II,YF0II,XF0III,YF0III,DCRS,PD,Q,P,ALL,BLL,CLL,DLL,CFFF
     1,ELL,FLL,GLL,HLL,PLL,CFFFF,ONE,TWO,THREE,FOUR,FIVE,EIGHT,TEN
     1,ELEVEN,VAR6,POW1,VAR14,VAR2,CFG,C2G,PRODAB,C3G,C33G,C4G,C44G
     1,C5G,C55G,C6G,C66G,C7G,C77G,C8G,C88G,C9G,C99G,C11G,C1111G
#include <pi.h>
      RRR=RHO*RHO
      BG=ALL*BLL
      RRRR=RRR*RHO
      XXX=BG*CFFFF
      RHOT=1.0D0/RHO13
      XB=1.0D0+DLL*RHOT
      AG=RHOT/(3.0D0*RHO)
      XB2=XB*XB
      XBRHO=-DLL*AG
      XD=DEXP(-CLL*RHOT)
      XDRHO=CLL*XD*AG
      XA=-ALL*RHO-XXX*RHO*XD
      XARHO=-ALL-XXX*XD-XXX*RHO*XDRHO
      TERM1=XARHO/XB-XBRHO*XA/XB2
      XE=FLL*(CLL/(RRR)+DLL/(RRR*XB))
     1  +RHO13/(6.0D0*RRR)
      XERHO=FLL*(-2.0D0*CLL/RRRR-2.0D0*DLL/(RRRR*XB)
     1                    -XBRHO*DLL/(RRR*XB2))
     1     -PLL*RHO13/RRRR
      XF=(BG/4.0D0)*XD*XE/XB
      XFRHO=(BG/4.0D0)*(XDRHO*XE/XB+XD*XERHO/XB-XD*XE*XBRHO/XB2)
C
      TERM1=TERM1+GMOD*GMOD*XFRHO
      U=2.0D0*XF
      RETURN
      END

C     OPENSHELL POTENTIAL FOR THE GRADIENT OF LYP FUNCTIONAL ENERGY
C     CODED UP BY G.J.LAMING, JUNE 1991
      SUBROUTINE GOSLYP(A1,A2,A3,A4,A5,RHOA,RHOB,RHOA53,RHOB53
     1,RHOA83,RHOB83,RHO13,RHO23,RHO43
     1,RHO113,RHO143,AMASQ,AMBSQ,AMA,AMB,AMAAMB)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
#include <pi.h>
      COMMON/DFTCON/BCX,CVX,AL,BL,CL,DL,CF,C3,CRS,PRE,THIRD,CINF,FF,AA
     1,BB,GG,DD,AK1,AK2,AK3,PP,FT,QQ,BETAP,CV,AI,BI,CI,X0I,AII,BII,CII
     1,X0II,AIII,BIII,CIII,X0III,FPPZ,FTT,FTTT,QQQI,QQQII,QQQIII,XF0I
     1,YF0I,XF0II,YF0II,XF0III,YF0III,DCRS,PD,Q,P,ALL,BLL,CLL,DLL,CFFF
     1,ELL,FLL,GLL,HLL,PLL,CFFFF,ONE,TWO,THREE,FOUR,FIVE,EIGHT,TEN
     1,ELEVEN,VAR6,POW1,VAR14,VAR2,CFG,C2G,PRODAB,C3G,C33G,C4G,C44G
     1,C5G,C55G,C6G,C66G,C7G,C77G,C8G,C88G,C9G,C99G,C11G,C1111G
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
      GGGG=DEXP(-CLL*VAR)
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
     1 RHOA53,RHOB53,RHOA83,RHOB83,
     1 RHO13,RHO23,RHO43,
     1 RHO113,RHO143,
     1 F3,F4,F5,F6,
     1  VAR,VAR3,VAR4,VAR5,VAR9,RR,PAB,RHO,ZA,ZB,
     1 Q1,Q2,Q1DA,Q2DA,XDA,XDB,ZADA,ZBDA,
     1 AMASQ,AMBSQ)
      CALL NGLYP(A2,RHOB,RHOA,AMB,AMA,AMAAMB,W,WDA,WDB,DTA,DTADA,
     1 RHOB53,RHOA53,RHOB83,RHOA83,
     1 RHO13,RHO23,RHO43,
     1 RHO113,RHO143,
     1 F3,F4,F5,F6,
     1 VAR,VAR3,VAR4,VAR5,VAR9,RR,PAB,RHO,ZB,ZA,
     1 Q2,Q1,Q2DB,Q1DB,XDB,XDA,ZBDB,ZADB,
     1 AMBSQ,AMASQ)
      A3=TWO*FB
      A4=TWO*FC
      A5=FD
      RETURN
      END
C     OPENSHELL POTENTIAL FOR THE GRADIENT OF LYP FUNCTIONAL ENERGY
C
C**********CODED   UP   BY   G.  J.  LAMING,   JUNE 1991*************
C
      SUBROUTINE NGLYP(T,RHOA,RHOB,AMA,AMB,AMAAMB,W,WDA,WDB,DTA,DTADA,
     1 RHOA53,RHOB53,RHOA83,RHOB83,
     1 RHO13,RHO23,RHO43,
     1 RHO113,RHO143,
     1 F3,F4,F5,F6,
     1  VAR,VAR3,VAR4,VAR5,VAR9,RR,PAB,RHO,ZA,ZB,
     1 Q1,Q2,Q1DA,Q2DA,XDA,XDB,ZADA,ZBDA,
     1 AMASQ,AMBSQ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DFTCON/BCX,CVX,AL,BL,CL,DL,CF,C3,CRS,PRE,THIRD,CINF,FF,AA
     1,BB,GG,DD,AK1,AK2,AK3,PP,FT,QQ,BETAP,CV,AI,BI,CI,X0I,AII,BII,CII
     1,X0II,AIII,BIII,CIII,X0III,FPPZ,FTT,FTTT,QQQI,QQQII,QQQIII,XF0I
     1,YF0I,XF0II,YF0II,XF0III,YF0III,DCRS,PD,Q,P,ALL,BLL,CLL,DLL,CFFF
     1,ELL,FLL,GLL,HLL,PLL,CFFFF,ONE,TWO,THREE,FOUR,FIVE,EIGHT,TEN
     1,ELEVEN,VAR6,POW1,VAR14,VAR2,CFG,C2G,PRODAB,C3G,C33G,C4G,C44G
     1,C5G,C55G,C6G,C66G,C7G,C77G,C8G,C88G,C9G,C99G,C11G,C1111G
#include <pi.h>
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
     1 (WDA*PAB*DTA*Q1+W*RHOB*DTA*Q1+W*PAB*DTADA*Q1+W*PAB*DTA*Q1DA)
      T7AA=C7G*
     1  (WDA*PAB*DTA*Q2+W*RHOB*DTA*Q2+W*PAB*DTADA*Q2+W*PAB*DTA*Q2DA)
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
C
      SUBROUTINE CONDFT
C     DFT CONSTANTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      logical vwn5
#include <pi.h>
      COMMON/DFTCON/BCX,CVX,AL,BL,CL,DL,CF,C3,CRS,PRE,THIRD,CINF,FF,AA
     1,BB,GG,DD,AK1,AK2,AK3,PP,FT,QQ,BETAP,CV,AI,BI,CI,X0I,AII,BII,CII
     1,X0II,AIII,BIII,CIII,X0III,FPPZ,FTT,FTTT,QQQI,QQQII,QQQIII,XF0I
     1,YF0I,XF0II,YF0II,XF0III,YF0III,DCRS,PD,Q,P,ALL,BLL,CLL,DLL,CFFF
     1,ELL,FLL,GLL,HLL,PLL,CFFFF,ONE,TWO,THREE,FOUR,FIVE,EIGHT,TEN
     1,ELEVEN,VAR6,POW1,VAR14,VAR2,CFG,C2G,PRODAB,C3G,C33G,C4G,C44G
     1,C5G,C55G,C6G,C66G,C7G,C77G,C8G,C88G,C9G,C99G,C11G,C1111G
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
      DCRS=DSQRT(CRS)
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

      QQQI=DSQRT(4.0D0*CI-BI*BI)
      XF0I=X0I*X0I+BI*X0I+CI
      YF0I=QQQI/(2.0D0*X0I+BI)
C
      AII=0.0310907D0
      BII=7.06042D0
      CII=18.0578D0
      X0II=-0.32500D0
      QQQII=DSQRT(4.0D0*CII-BII*BII)
      XF0II=X0II*X0II+BII*X0II+CII
      YF0II=QQQII/(2.0D0*X0II+BII)
C
      AIII=-0.033774D0
      BIII=1.131071D0
      CIII=13.0045D0
      X0III=-0.0047584D0
      QQQIII=DSQRT(4.0D0*CIII-BIII*BIII)
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
C
C     *************************
C     ***** KT functional *****
C     *************************
C
C  /* Deck ekt */
      SUBROUTINE EKT(EK,RHO,RHO13,RHOGRD)
#include <implicit.h>
#include <pi.h>
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, D4 = 4.0D0,
     &           GAMMA = -0.006D0, DELTA = 0.1D0)
      GAMBAR = D2**(D1/D3)*GAMMA
      DLTBAR = D2**(D4/D3)*DELTA
      EK = GAMBAR*(RHOGRD**2)/(RHO*RHO13 + DLTBAR) 
      RETURN
      END
C  /* Deck gkt */
      SUBROUTINE GKT(VK1,VK2,RHO,RHO13,RHOGRD)
#include <implicit.h>
#include <pi.h>
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, D4 = 4.0D0,
     &           GAMMA = -0.006D0, DELTA = 0.1D0)
      GAMBAR = D2**(D1/D3)*GAMMA
      DLTBAR = D2**(D4/D3)*DELTA
      VK1 = -(D4/D3)*GAMBAR*(RHOGRD**2)*RHO13/(RHO*RHO13 + DLTBAR)**2
      VK2 = D2*D2*GAMBAR/(RHO*RHO13 + DLTBAR)
      RETURN
      END
C  /* Deck vkt */
      SUBROUTINE VKT(VK1,RHO,RHO13,RHOGRD,RHOLAP)
#include <implicit.h>
#include <pi.h>
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, D4 = 4.0D0,
     &           GAMMA = -0.006D0, DELTA = 0.1D0)
      GAMBAR = D2**(D1/D3)*GAMMA
      DLTBAR = D2**(D4/D3)*DELTA
      VK1 = (D4/D3)*GAMBAR*(RHOGRD**2)*RHO13/(RHO*RHO13 + DLTBAR)**2
     &    - D2*GAMBAR*RHOLAP/(RHO*RHO13 + DLTBAR)
      RETURN
      END
C  /* Deck v1kt */
      SUBROUTINE V1KT(FR0,FZ0,FRR,FRZ,FZZ,RHO,RHO13,RHOGRD)
#include <implicit.h>
#include <pi.h>
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, D4 = 4.0D0,
     &           D8 = 8.0D0, D9 = 9.0D0,
     &           GAMMA = -0.006D0, DELTA = 0.1D0)
      GAMBAR = D2**(D1/D3)*GAMMA
      DLTBAR = D2**(D4/D3)*DELTA
      RHODLT = D1/(RHO*RHO13 + DLTBAR)
      FR0 = -(D4/D3)*GAMBAR*(RHOGRD**2)*RHO13*(RHODLT**2)
      FZ0 = D2*GAMBAR*RHODLT*RHOGRD
      FRR = -(D4/D9)*GAMBAR*(RHOGRD**2)
     &      *(RHO13/RHO - D8*RHO13*RHODLT)*(RHODLT**2)
      FRZ = -(D8/D3)*GAMBAR*RHOGRD*RHO13*(RHODLT**2)
      FZZ = D2*GAMBAR*RHODLT
      RETURN
      END
