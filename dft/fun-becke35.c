/*


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

*/
/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* fun-becke.c:

   REFERENCE: R.D. Adamson, P.M.W. Gill and J.A. Pople, Chem. Phys. Lett., 284, 6 (1998).
 
   Implemented by David Wilson (david.wilson@latrobe.edu.au), Jun 2005.
   Includes modified beta value of 0.0035 (compared to 0.0042 in Becke functional).

   Based on implementation of Becke(88) functional and its derivatives.
   or exactly: Becke GGA correction to the functional. (total Becke(88)
   energy is E_LDA+E_BCK).
   (c) Pawel Salek, pawsa@theochem.kth.se, aug 2001
   Z. Rinkevicius adapted for open shell systems: energy, first derivatives.

   NOTE:
   this file may seem unnecessarily complex but the structure does pay off
   when implementing multiple functionals depending on different parameters.
*/

/* strictly conform to XOPEN ANSI C standard */
#if !defined(SYS_DEC)
/* XOPEN compliance is missing on old Tru64 4.0E Alphas */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1
#endif

#include <math.h>
#include <stddef.h>

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int  becke35_isgga(void) { return 1; }
static int  becke35_read(const char* conf_line);
static real becke35_energy(const FunDensProp* dens_prop);
static void becke35_first(FunFirstFuncDrv *ds, real factor, 
                        const FunDensProp* dens_prop);
static void becke35_second(FunSecondFuncDrv *ds, real factor,
                         const FunDensProp* dens_prop);
static void becke35_third(FunThirdFuncDrv *ds, real factor,
                        const FunDensProp* dens_prop);
static void becke35_fourth(FunFourthFuncDrv *ds, real factor,
                         const FunDensProp* dens_prop);

Functional mBeckeFunctional = {
    "mBecke",      /* name */
    becke35_isgga,  /* gga-corrected */
    becke35_read,   /* set bloody common blocks */
    NULL,         /* reporter */
    becke35_energy, 
    becke35_first,
    becke35_second,
    becke35_third,
    becke35_fourth
};

/* IMPLEMENTATION PART */

static int
becke35_read(const char* conf_line)
{
    fun_set_hf_weight(0.0);
    return 1;
}

/* becke35_energy:
   note that in reality E_BCK = E_BCK,alpha + E_BCK,beta
   i.e the is linear in alpha and beta densities.

   Becke threshold is needed to avoid numerical problems on 0/0
   divisions.  The problems are small but it is better to be on the
   safe side.
*/
static const real BECKE_THRESHOLD = 1e-14;
static const real BETA = 0.0035;
static real
becke35_energy(const FunDensProp* dp)
{
   real ea,eb;
   if (dp->rhob<BECKE_THRESHOLD)
     eb = 0.0;
   else {
     real xb = dp->gradb*pow(dp->rhob,-4.0/3.0);
     real rb = pow(dp->rhob,4.0/3.0);
     real denomb = 1.0 +6.0*xb*BETA*asinh(xb);
     eb = rb*xb*xb/denomb; 
   } 
   if (dp->rhoa<BECKE_THRESHOLD) 
     ea=0;
   else {
       real xa = dp->grada*pow(dp->rhoa,-4.0/3.0);
       real ra = pow(dp->rhoa,4.0/3.0);
       real denoma = 1.0 +6.0*BETA*xa*asinh(xa);
       ea = ra*xa*xa/denoma;
   }
   return -BETA*(ea+eb);
}



static void
becke35_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real alpha, aa2, asha, sq1aa2;
    real alphb, ab2, ashb, sq1ab2;
    real denoma, denoma2;
    real denomb, denomb2; 
    real alphaa10, alphaa01;
    real alphab10, alphab01;
    real ffa, ff1a;
    real ffb, ff1b;

    if (dp->rhoa >BECKE_THRESHOLD) {
      alpha = dp->grada*pow(dp->rhoa,-4.0/3.0);
      aa2 = alpha*alpha;
      asha = asinh(alpha);
      alphaa10 = -4.0/3.0*alpha/dp->rhoa;   /* dalpha/drhoa   */
      alphaa01 = pow(dp->rhoa,-4.0/3.0);    /* dalpha/dgrho  */
      sq1aa2 = sqrt(1 + aa2);
      denoma= 1 + 6*alpha*BETA*asha;
      denoma2 = denoma*denoma;
      ffa  = -alpha*BETA/(1+6*alpha*BETA*asha);
      ff1a = BETA*(6*aa2*BETA - sq1aa2)/(sq1aa2*denoma2);
      ds->df1000 += factor*dp->grada* ff1a *alphaa10;
      ds->df0010 += factor*(ffa + dp->grada*ff1a*alphaa01);
    }
    if (dp->rhob >BECKE_THRESHOLD) {
        alphb = dp->gradb*pow(dp->rhob,-4.0/3.0);   
        ab2 = alphb*alphb;
        ashb = asinh(alphb);
        alphab10 = -4.0/3.0*alphb/dp->rhob; 
        alphab01 = pow(dp->rhob,-4.0/3.0);
        sq1ab2 = sqrt(1 + ab2);
        denomb= 1 + 6*alphb*BETA*ashb;      
        denomb2 = denomb*denomb;  
        ffb  = -alphb*BETA/(1+6*alphb*BETA*ashb);   
        ff1b = BETA*(6*ab2*BETA - sq1ab2)/(sq1ab2*denomb2); 
        ds->df0100 += factor*dp->gradb* ff1b *alphab10;
        ds->df0001 += factor*(ffb + dp->gradb*ff1b*alphab01); 
    } 
}

static void
becke35_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    const real BETA2 = BETA*BETA;
    real alpha, a2, a3, a4, asha, asha2, sq1a2;
    real denom, denom2, denom3;
    real alpha10, alpha20, alpha01, alpha11, alpha10_2;
    real ff, ff1, ff2;
    
    if (dp->rhoa >BECKE_THRESHOLD) {    
      alpha = dp->grada*pow(dp->rhoa,-4.0/3.0);
      a2 = alpha*alpha;
      a3 = alpha*a2;
      a4 = alpha*a3;
      asha = asinh(alpha);
      asha2= asha*asha;
      alpha10 = -4.0/3.0*alpha/dp->rhoa;   /* dalpha/drhoa   */
      alpha20 = -7.0/3.0*alpha10/dp->rhoa; /* d²alpha/drhoa² */
      alpha01 = pow(dp->rhoa,-4.0/3.0);    /* dalpha/dgrho  */
      alpha11 = -4.0/3.0*alpha01/dp->rhoa; /* d²alpha/(drho dgrho) */
      alpha10_2 = alpha10*alpha10;
      sq1a2 = sqrt(1 + a2);
      denom= 1 + 6*alpha*BETA*asha;
      denom2 = denom*denom;
      denom3 = denom2*denom;

      ff  = -alpha*BETA/(1+6*alpha*BETA*asha);
      ff1 = BETA*(6*a2*BETA - sq1a2)/(sq1a2*denom2);

      ff2 = (6*BETA2*(4*alpha + a3*(3 - 12*sq1a2*BETA) + 
	     2*(pow(1 + a2,1.5) - 3*a4*BETA)*asha))/
	    (pow(1 + a2,1.5)*denom3);
    
      ds->df1000 += factor*dp->grada* ff1 *alpha10;
      ds->df0010 += factor*(ff + dp->grada*ff1*alpha01);
      ds->df1010 += factor*(ff1*alpha10 + 
                    dp->grada*(ff2*alpha10*alpha01+ff1*alpha11));
      ds->df2000 += factor*dp->grada*(ff2* alpha10_2 + ff1 * alpha20);
      ds->df0020 += factor*(2*ff1*alpha01 + 
		    dp->grada*(ff2*alpha01*alpha01));
  }
   /* note: reuse of variables for beta part, maybe a bit .. */
  if (dp->rhob >BECKE_THRESHOLD) {    
      alpha = dp->gradb*pow(dp->rhob,-4.0/3.0);
      a2 = alpha*alpha;
      a3 = alpha*a2;
      a4 = alpha*a3;
      asha = asinh(alpha);
      asha2= asha*asha;
      alpha10 = -4.0/3.0*alpha/dp->rhob;   
      alpha20 = -7.0/3.0*alpha10/dp->rhob; 
      alpha01 = pow(dp->rhob,-4.0/3.0);    
      alpha11 = -4.0/3.0*alpha01/dp->rhob; 
      alpha10_2 = alpha10*alpha10;
      sq1a2 = sqrt(1 + a2);
      denom= 1 + 6*alpha*BETA*asha;
      denom2 = denom*denom;
      denom3 = denom2*denom;

      ff  = -alpha*BETA/(1+6*alpha*BETA*asha);
      ff1 = BETA*(6*a2*BETA - sq1a2)/(sq1a2*denom2);

      ff2 = (6*BETA2*(4*alpha + a3*(3 - 12*sq1a2*BETA) + 
	     2*(pow(1 + a2,1.5) - 3*a4*BETA)*asha))/
	    (pow(1 + a2,1.5)*denom3);
    
      ds->df0100 += factor*dp->gradb* ff1 *alpha10;
      ds->df0001 += factor*(ff + dp->gradb*ff1*alpha01);
      ds->df0101 += factor*(ff1*alpha10 + 
                    dp->gradb*(ff2*alpha10*alpha01+ff1*alpha11));
      ds->df0200 += factor*dp->gradb*(ff2* alpha10_2 + ff1 * alpha20);
      ds->df0002 += factor*(2*ff1*alpha01 + 
		    dp->gradb*(ff2*alpha01*alpha01));
    }  
}

/* becke35_third:
   Becke functional derivatives.
   Input: rho, rhogrd.
   Output:
   df1000 - d/drho F
   df1010 - d/dgrd F
   df2000 - d^2/drho^2 F
   df1010 - d^2/drho d/dgrd F
   df0020 - d^2/dgrd^2 F
   df3000 - d^3/drho^3 F

   NOTES: Instead of passing 6+ numbers, a pointer to a structure
   ggaSecDrv could be passed instead. The numbers are closely related
   to each other.

   The Becke functional is expressed as a function of grho and alpha:
   F_B = grho*f
 
   The functional derivatives are expressed through the partial derivatives.
   For example (g==grho):
   df/drho   = g df/dalpha * dalpha/drho
   d²f/drho² = g(d²f/dalpha² * (dalpha/drho)² + df/dalpha * d²alpha/drho²)
   d³f/drho³ = g(d³f/dalpha³ * (dalpha/drho)³ 
             + 3d²f/dalpha²*(dalpha/drho)*d²alpha/drho²
             + df/dalpha * d³alpha/drho³)

   The derivatives with respect to grho==g are:
   df/dgrho  = f + g df/dalpha * dalpha/dgrho
   d²f/dgrho²= 2*df/dalpha*dalpha/dg + g(d²f/dalpha²(dalpha/dg)²)

   NOTE that d²alpha/dg² = 0, and some terms are missing above.
   d³f/dgrho³= 3*d²f/dalpha²*(dalpha/dg)² (+0)
             + (+0)+ g(d³f/dalpha³(dalpha/dg)³ (+0))

   NOTE: it adds to ds, not sets it.
   This routine has unrestricted interface but works only for
   restricted case.
*/
static void
becke35_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    const real BETA2 = BETA*BETA;
    real alpha, a2, a3, a4, a5, asha, asha2, sq1a2;
    real denom, denom2, denom3, denom4;
    real alpha10, alpha20, alpha30, alpha01, alpha11, alpha21, alpha10_2;
    real ff, ff1, ff2, ff3;

    if (dp->rhoa > BECKE_THRESHOLD) {
      alpha = dp->grada*pow(dp->rhoa,-4.0/3.0);
      a2 = alpha*alpha;
      a3 = alpha*a2;
      a4 = alpha*a3;
      a5 = alpha*a4;
      asha = asinh(alpha);
      asha2= asha*asha;
      alpha10 = -4.0/3.0*alpha/dp->rhoa;   /* dalpha/drhoa   */
      alpha20 = -7.0/3.0*alpha10/dp->rhoa; /* d²alpha/drhoa² */
      alpha30 = -10./3.0*alpha20/dp->rhoa; /* d³alpha/drho³ */
      alpha01 = pow(dp->rhoa,-4.0/3.0);    /* dalpha/dgrho  */
      alpha11 = -4.0/3.0*alpha01/dp->rhoa; /* d²alpha/(drho dgrho) */
      alpha21 = alpha20/dp->grada;         /* d³alpha/(drho² dgrho) */
      alpha10_2 = alpha10*alpha10;
      sq1a2 = sqrt(1 + a2);
      denom= 1 + 6*alpha*BETA*asha;
      denom2 = denom*denom;
      denom3 = denom2*denom;
      denom4 = denom3*denom;

      ff  = -alpha*BETA/(1+6*alpha*BETA*asha);
      ff1 = BETA*(6*a2*BETA - sq1a2)/(sq1a2*denom2);

      ff2 = (6*BETA2*(4*alpha + a3*(3 - 12*sq1a2*BETA) + 
            2*(pow(1 + a2,1.5) - 3*a4*BETA)*asha))/
	    (pow(1 + a2,1.5)*denom3);
    
      ff3 = (6*BETA2*(6 + 5*a2 + 2*a4 
	    - 12*alpha*(6 + 16*a2 + 7*a4)*BETA*asha + 
	    36*sq1a2*BETA*(-a2*(3 + 2*a2) + 6*a5*BETA*asha - 
	    pow(1 + a2,2.0)*asha2) + 
	    36*a4*BETA2*(6*(1 + a2) + (-1 + 2*a2)*asha2)))/
            (pow(1 + a2,2.5)*denom4);

      ds->df1000 += factor*dp->grada* ff1 *alpha10;
      ds->df0010 += factor*(ff + dp->grada*ff1*alpha01);
      ds->df1010 += 
	factor*(ff1*alpha10 + dp->grada*(ff2*alpha10*alpha01+ff1*alpha11));
      ds->df2000 += factor*dp->grada*(ff2* alpha10_2 + ff1 * alpha20);
      ds->df0020 += factor*(2*ff1*alpha01 + dp->grada*(ff2*alpha01*alpha01));
    

      ds->df2010 += factor*(ff2*alpha10_2 + ff1 * alpha20 +
                          dp->grada*(ff3*alpha10_2*alpha01 + 
                                 2*ff2*alpha11*alpha10 +
                                 ff2*alpha20*alpha01 + ff1*alpha21));

      ds->df1020 += factor*(2*ff2*alpha10*alpha01 + 2*ff1*alpha11 +
                          dp->grada*(ff3*alpha01*alpha01*alpha10 + 
                                 2*ff2*alpha11*alpha01));
    
      ds->df3000 += factor*dp->grada*(ff3 * alpha10_2*alpha10 +
                                3*ff2*alpha10*alpha20 +
                                ff1*alpha30);

      ds->df0030 += factor*(3*ff2*alpha01*alpha01 +
                          dp->grada*ff3*pow(alpha01,3.0));
    }
 if (dp->rhob > BECKE_THRESHOLD) {
      alpha = dp->gradb*pow(dp->rhob,-4.0/3.0);
      a2 = alpha*alpha;
      a3 = alpha*a2;
      a4 = alpha*a3;
      a5 = alpha*a4;
      asha = asinh(alpha);
      asha2= asha*asha;
      alpha10 = -4.0/3.0*alpha/dp->rhob;   /* dalpha/drhoa   */
      alpha20 = -7.0/3.0*alpha10/dp->rhob; /* d²alpha/drhoa² */
      alpha30 = -10./3.0*alpha20/dp->rhob; /* d³alpha/drho³ */
      alpha01 = pow(dp->rhob,-4.0/3.0);    /* dalpha/dgrho  */
      alpha11 = -4.0/3.0*alpha01/dp->rhob; /* d²alpha/(drho dgrho) */
      alpha21 = alpha20/dp->gradb;         /* d³alpha/(drho² dgrho) */
      alpha10_2 = alpha10*alpha10;
      sq1a2 = sqrt(1 + a2);
      denom= 1 + 6*alpha*BETA*asha;
      denom2 = denom*denom;
      denom3 = denom2*denom;
      denom4 = denom3*denom;

      ff  = -alpha*BETA/(1+6*alpha*BETA*asha);
      ff1 = BETA*(6*a2*BETA - sq1a2)/(sq1a2*denom2);

      ff2 = (6*BETA2*(4*alpha + a3*(3 - 12*sq1a2*BETA) + 
            2*(pow(1 + a2,1.5) - 3*a4*BETA)*asha))/
	    (pow(1 + a2,1.5)*denom3);
    
      ff3 = (6*BETA2*(6 + 5*a2 + 2*a4 
	    - 12*alpha*(6 + 16*a2 + 7*a4)*BETA*asha + 
	    36*sq1a2*BETA*(-a2*(3 + 2*a2) + 6*a5*BETA*asha - 
	    pow(1 + a2,2.0)*asha2) + 
	    36*a4*BETA2*(6*(1 + a2) + (-1 + 2*a2)*asha2)))/
            (pow(1 + a2,2.5)*denom4);

      ds->df0100 += factor*dp->gradb* ff1 *alpha10;
      ds->df0001 += factor*(ff + dp->gradb*ff1*alpha01);
      ds->df0101 += 
	factor*(ff1*alpha10 + dp->gradb*(ff2*alpha10*alpha01+ff1*alpha11));
      ds->df0200 += factor*dp->gradb*(ff2* alpha10_2 + ff1 * alpha20);
      ds->df0002 += factor*(2*ff1*alpha01 + dp->gradb*(ff2*alpha01*alpha01));
    

      ds->df0201 += factor*(ff2*alpha10_2 + ff1 * alpha20 +
                          dp->gradb*(ff3*alpha10_2*alpha01 + 
                                 2*ff2*alpha11*alpha10 +
                                 ff2*alpha20*alpha01 + ff1*alpha21));

      ds->df0102 += factor*(2*ff2*alpha10*alpha01 + 2*ff1*alpha11 +
                          dp->gradb*(ff3*alpha01*alpha01*alpha10 + 
                                 2*ff2*alpha11*alpha01));
    
      ds->df0300 += factor*dp->gradb*(ff3 * alpha10_2*alpha10 +
                                3*ff2*alpha10*alpha20 +
                                ff1*alpha30);

      ds->df0003 += factor*(3*ff2*alpha01*alpha01 +
                          dp->gradb*ff3*pow(alpha01,3.0));
    }
}

/* Different paritioning has been used for Fourth Derivatives

 Becke functional of the form Becke=-Beta*F
 where F(rho,grad)=rho^(4/3) * X^2 / ( 1 + 6*Beta*X*ArcSinh(X) )
 
 F(rho, grad) has been partitoned in following manner:

 F(rho, grad)	  = Nom(rho, grad)*Denom(rho, grad)

 Nom(rho, grad)   = grad^2 / rho^(4/3)
 Denom(rho, grad) = 1/(1 + 6*Beta*XARX(rho, grad))
 XARX(rho, grad)  = X(rho, grad)*ArcSinh(X(rho, grad))
 X(rho, grad)     = grad/rho^(4/3)

 *************************************************************

 NAME       :	Variable name
 NAMEXY     : 	X-th partial derivative with respect to rho
 		Y-th partial derivative with respect to grad
		of Variable (function) NAME
 NAME_pX    :	Xth power of NAME  (NAME^X)
 NAME_pmX   :	NAME^(-X)
 NAME_pmXfY :	NAME^(-X/Y)
 INAME	    :	1/NAME
 
 *************************************************************
 (c)2002 by B. Jansik
 */	
	
static void
becke35_fourth(FunFourthFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real ASINH_X, B, BXARX;
    real BXARX_pm2, BXARX_pm3, BXARX_pm4, BXARX_pm5;
    real B_p2, B_p3, B_p4;
    real DENX, DENX_pm1f2, DENX_pm3f2, DENX_pm5f2, DENX_pm7f2;
    real Denom, Denom01, Denom02, Denom03, Denom04, Denom10;
    real Denom11, Denom12, Denom13, Denom20, Denom21, Denom22;
    real Denom30, Denom31, Denom40, IDENX;
    real Nom, Nom01, Nom02, Nom10, Nom11, Nom12, Nom20, Nom21, Nom22;
    real Nom30, Nom31, Nom40;
    real X, X01, X01_p2, X01_p3, X01_p4;
    real X10, X10_p2, X10_p3, X10_p4;
    real X11, X11_p2, X20,  X20_p2;
    real  X21, X30, X31, X40, XARX01, XARX01_p2, XARX01_p3, XARX01_p4;
    real  XARX02, XARX02_p2, XARX03, XARX04, XARX10;
    real  XARX10_p2, XARX10_p3, XARX10_p4, XARX11, XARX11_p2, XARX12;
    real  XARX13, XARX20, XARX20_p2, XARX21, XARX22, XARX30, XARX31, XARX40;
    real  X_p2, X_p3, X_p4, gro, gro_p2;
    real  iro, ro, ro_pm10f3, ro_pm13f3, ro_pm16f3, ro_pm4f3,  ro_pm7f3;
    real  out;
    FunThirdFuncDrv ds_third;

/* Setting up lower order derivatives
 * BeckeFunctional.third calculate third and also lower order derivatives  */

    drv3_clear(&ds_third);
    BeckeFunctional.third(&ds_third, factor, dp);

    ds->df1000 += ds_third.df1000;
    ds->df0010 += ds_third.df0010;
    ds->df2000 += ds_third.df2000;
    ds->df0020 += ds_third.df0020;
    ds->df2010 += ds_third.df2010;
    ds->df1020 += ds_third.df1020;
    ds->df3000 += ds_third.df3000;
    ds->df0030 += ds_third.df0030;

    ds->df0100 += ds_third.df0100;
    ds->df0001 += ds_third.df0001;
    ds->df0200 += ds_third.df0200;
    ds->df0002 += ds_third.df0002;
    ds->df0201 += ds_third.df0201;
    ds->df0102 += ds_third.df0102;
    ds->df0300 += ds_third.df0300;
    ds->df0003 += ds_third.df0003;
    ds->df1010 += ds_third.df1010;
    ds->df0101 += ds_third.df0101;
   
/* Powers of Beta */
    B=BETA;
    B_p2=B*B;
    B_p3=B_p2*B;
    B_p4=B_p3*B;
	
    if (dp->rhoa > BECKE_THRESHOLD) {

        ro=dp->rhoa;
        gro=dp->grada;


/* Powers of ro, gro */
        iro=1.0/ro;
        ro_pm4f3=pow(ro,(-4.0/3.0));
/* ro_pm4f3=1.0/cbrt(ro*ro*ro*ro); */
        ro_pm7f3=iro*ro_pm4f3;
        ro_pm10f3=iro*ro_pm7f3;
        ro_pm13f3=iro*ro_pm10f3;
        ro_pm16f3=iro*ro_pm13f3;

        gro_p2=gro*gro;
/* End of Powers of ro, gro */


/*  X and derivatives */
        X=gro*(ro_pm4f3);

        X01=(ro_pm4f3);

        X10=-(4.0/3.0)*gro*ro_pm7f3;
        X20=(28.0/9.0)*gro*ro_pm10f3;
        X30=-(280.0/27.0)*gro*ro_pm13f3;
        X40=(3640.0/81.0)*gro*ro_pm16f3;

        X11=-(4.0/3.0)*ro_pm7f3;
        X21=(28.0/9.0)*ro_pm10f3;
        X31=-(280.0/27.0)*ro_pm13f3;

/*  X deriv. powers */
        X_p2=X*X;
        X_p3=X_p2*X;
        X_p4=X_p3*X;

        X01_p2=X01*X01;
        X01_p3=X01_p2*X01;
        X01_p4=X01_p3*X01;

        X10_p2=X10*X10;
        X10_p3=X10_p2*X10;
        X10_p4=X10_p3*X10;

        X20_p2=X20*X20;

        X11_p2=X11*X11;

/*  End of X deriv powers */
/*  End of X derivatives */


        DENX=1.0+X_p2;
        IDENX=1.0/DENX;

        DENX_pm1f2=1.0/sqrt(DENX);
        DENX_pm3f2=IDENX*(DENX_pm1f2);
        DENX_pm5f2=IDENX*(DENX_pm3f2);
        DENX_pm7f2=IDENX*(DENX_pm5f2);

        ASINH_X=asinh(X);
        BXARX=1.0+(6.0*B*X*ASINH_X);
        Denom=1.0/BXARX;

        BXARX_pm2=Denom*Denom;
        BXARX_pm3=BXARX_pm2*Denom;
        BXARX_pm4=BXARX_pm3*Denom;
        BXARX_pm5=BXARX_pm4*Denom;

/*  XARX derivatives */

        XARX01=ASINH_X*X01  + X*X01*DENX_pm1f2;

        XARX02=-X_p2*X01_p2*DENX_pm3f2 + 
            2.0*(X01_p2)* DENX_pm1f2;

        XARX03=3.0*(X_p3)*(X01_p3)* DENX_pm5f2 - 
            4.0*X*(X01_p3)* DENX_pm3f2;

        XARX04=-15.0*(X_p4)*(X01_p4) * DENX_pm7f2 + 
            21.0*(X_p2)*(X01_p4)* DENX_pm5f2 - 
            4.0*(X01_p4) * DENX_pm3f2;

        XARX10=ASINH_X*X10  + X*X10*DENX_pm1f2;

        XARX20=-X_p2*X10_p2*DENX_pm3f2 + 
            2.0*(X10_p2)* DENX_pm1f2 + 
            ASINH_X*X20 + 
            X*X20*DENX_pm1f2;

        XARX30=3.0*(X_p3)*(X10_p3)* DENX_pm5f2 - 
            4.0*X*(X10_p3)* DENX_pm3f2 - 
            3.0*X_p2*X10*X20 * DENX_pm3f2 + 
            6.0*X10*X20*DENX_pm1f2 + 
            ASINH_X*X30 + 
            X*X30* DENX_pm1f2;

        XARX40=-15.0*(X_p4)*(X10_p4) * DENX_pm7f2 + 
            21.0*(X_p2)*(X10_p4)* DENX_pm5f2 - 
            4.0*(X10_p4) * DENX_pm3f2 + 
            18.0*X_p3*X10_p2*X20*DENX_pm5f2 - 
            24.0*X*X10_p2*X20*DENX_pm3f2 - 
            3.0*X_p2*X20_p2*DENX_pm3f2 + 
            6.0*X20_p2*DENX_pm1f2 - 
            4.0*X_p2*X10*X30*DENX_pm3f2 + 
            8.0*X10*X30*DENX_pm1f2 + 
            ASINH_X*X40 + 
            X*X40*DENX_pm1f2;

        XARX11=-X_p2*X01*X10*DENX_pm3f2 + 
            2.0*X01*X10*DENX_pm1f2 + 
            ASINH_X*X11 + 
            X*X11*DENX_pm1f2;

        XARX21=3.0*X_p3*X01*X10_p2*DENX_pm5f2 - 
            4.0*X*X01*X10_p2*DENX_pm3f2 - 
            2.0*X_p2*X10*X11*DENX_pm3f2 + 
            4.0*X10*X11*DENX_pm1f2 - 
            X_p2*X01*X20*DENX_pm3f2 + 
            2.0*X01*X20*DENX_pm1f2 + 
            ASINH_X*X21 + 
            X*X21*DENX_pm1f2;

        XARX31=-15.0*X_p4*X01*X10_p3*DENX_pm7f2 + 
            21.0*X_p2*X01*X10_p3*DENX_pm5f2 - 
            4.0*X01*X10_p3*DENX_pm3f2 + 
            9.0*X_p3*X10_p2*X11*DENX_pm5f2 - 
            12.0*X*X10_p2*X11*DENX_pm3f2 + 
            9.0*X_p3*X01*X10*X20*DENX_pm5f2 - 
            12.0*X*X01*X10*X20*DENX_pm3f2 - 
            3.0*X_p2*X11*X20*DENX_pm3f2 + 
            6.0*X11*X20*DENX_pm1f2 - 
            3.0*X_p2*X10*X21*DENX_pm3f2 + 
            6.0*X10*X21*DENX_pm1f2 - 
            X_p2*X01*X30*DENX_pm3f2 + 
            2.0*X01*X30*DENX_pm1f2 + 
            ASINH_X*X31 + 
            X*X31*DENX_pm1f2;
	
        XARX12=3.0*X_p3*X01_p2*X10*DENX_pm5f2 - 
            4.0*X*X01_p2*X10*DENX_pm3f2 - 
            2.0*X_p2*X01*X11*DENX_pm3f2 + 
            4.0*X01*X11*DENX_pm1f2;

        XARX22=-15.0*X_p4*X01_p2*X10_p2*DENX_pm7f2 + 
            21.0*X_p2*X01_p2*X10_p2*DENX_pm5f2 - 
            4.0*X01_p2*X10_p2*DENX_pm3f2 + 
            12.0*X_p3*X01*X10*X11*DENX_pm5f2 - 
            16.0*X*X01*X10*X11*DENX_pm3f2 - 
            2.0*X_p2*X11_p2*DENX_pm3f2 + 
            4.0*X11_p2*DENX_pm1f2 + 
            3.0*X_p3*X01_p2*X20*DENX_pm5f2 - 
            4.0*X*X01_p2*X20*DENX_pm3f2 - 
            2.0*X_p2*X01*X21*DENX_pm3f2 + 
            4.0*X01*X21*DENX_pm1f2;

        XARX13=-15.0*X_p4*X01_p3*X10*DENX_pm7f2 + 
            21.0*X_p2*X01_p3*X10*DENX_pm5f2 - 
            4.0*X01_p3*X10*DENX_pm3f2 + 
            9.0*X_p3*X01_p2*X11*DENX_pm5f2 - 
            12.0*X*X01_p2*X11*DENX_pm3f2;

/* Xarx der powers */
        XARX10_p2=XARX10*XARX10;
        XARX10_p3=XARX10_p2*XARX10;
        XARX10_p4=XARX10_p3*XARX10;

        XARX20_p2=XARX20*XARX20;

        XARX01_p2=XARX01*XARX01;
        XARX01_p3=XARX01_p2*XARX01;
        XARX01_p4=XARX01_p3*XARX01;

        XARX02_p2=XARX02*XARX02;

        XARX11_p2=XARX11*XARX11;
/* End of Xarx der powers */

/* End of XARX derivatives */

/* Nom derivatives */

        Nom=(gro_p2)*(ro_pm4f3);

        Nom01=2.0*gro*(ro_pm4f3);
        Nom02=2.0*(ro_pm4f3);

        Nom10=-(4.0/3.0)*gro_p2*ro_pm7f3;
        Nom20=(28.0/9.0)*gro_p2*ro_pm10f3;
        Nom30=-(280.0/27.0)*gro_p2*ro_pm13f3;
        Nom40=(3640.0/81.0)*gro_p2*ro_pm16f3;

        Nom11=-(8.0/3.0)*gro*ro_pm7f3;
        Nom21=(56.0/9.0)*gro*ro_pm10f3;
        Nom31=-(560.0/27.0)*gro*ro_pm13f3;

        Nom12=-(8.0/3.0)*ro_pm7f3;
        Nom22=(56.0/9.0)*ro_pm10f3;

/* End of Nom derivatives */

/* Denom derivatives */


        Denom01=-6.0*B*XARX01*BXARX_pm2;

        Denom02=72.0*B_p2*XARX01_p2*BXARX_pm3 - 
            6.0*B*XARX02*BXARX_pm2;

        Denom03=-1296.0*B_p3*XARX01_p3*BXARX_pm4 + 
            216.0*B_p2*XARX01*XARX02*BXARX_pm3 - 
            6.0*B*XARX03*BXARX_pm2;


        Denom04=31104.0*B_p4*XARX01_p4*BXARX_pm5 - 
            7776.0*B_p3*XARX01_p2*XARX02*BXARX_pm4 + 
            216.0*B_p2*XARX02_p2*BXARX_pm3 + 
            288.0*B_p2*XARX01*XARX03*BXARX_pm3 - 
            6.0*B*XARX04*BXARX_pm2;


        Denom10=-6.0*B*XARX10*BXARX_pm2;

        Denom20=72.0*B_p2*XARX10_p2*BXARX_pm3 - 
            6.0*B*XARX20*BXARX_pm2;

        Denom30=-1296.0*B_p3*XARX10_p3*BXARX_pm4 + 
            216.0*B_p2*XARX10*XARX20*BXARX_pm3 - 
            6.0*B*XARX30*BXARX_pm2;


        Denom40=31104.0*B_p4*XARX10_p4*BXARX_pm5 - 
            7776.0*B_p3*XARX10_p2*XARX20*BXARX_pm4 + 
            216.0*B_p2*XARX20_p2*BXARX_pm3 + 
            288.0*B_p2*XARX10*XARX30*BXARX_pm3 - 
            6.0*B*XARX40*BXARX_pm2;

        Denom11=72.0*B_p2*XARX01*XARX10*BXARX_pm3 - 
            6.0*B*XARX11*BXARX_pm2;

        Denom21=-1296.0*B_p3*XARX01*XARX10_p2*BXARX_pm4 + 
            144.0*B_p2*XARX10*XARX11*BXARX_pm3 + 
            72.0*B_p2*XARX01*XARX20*BXARX_pm3 - 
            6.0*B*XARX21*BXARX_pm2;

        Denom31=31104.0*B_p4*XARX01*XARX10_p3*BXARX_pm5 - 
            3888.0*B_p3*XARX10_p2*XARX11*BXARX_pm4 - 
            3888.0*B_p3*XARX01*XARX10*XARX20*BXARX_pm4 + 
            216.0*B_p2*XARX11*XARX20*BXARX_pm3 + 
            216.0*B_p2*XARX10*XARX21*BXARX_pm3 + 
            72.0*B_p2*XARX01*XARX30*BXARX_pm3 - 
            6.0*B*XARX31*BXARX_pm2;

        Denom12=-1296.0*B_p3*XARX01_p2*XARX10*BXARX_pm4 + 
            144.0*B_p2*XARX01*XARX11*BXARX_pm3 + 
            72.0*B_p2*XARX02*XARX10*BXARX_pm3 - 
            6.0*B*XARX12*BXARX_pm2;

        Denom22=31104.0*B_p4*XARX01_p2*XARX10_p2*BXARX_pm5 - 
            1296.0*B_p3*XARX02*XARX10_p2*BXARX_pm4 - 
            5184.0*B_p3*XARX01*XARX10*XARX11*BXARX_pm4 + 
            144.0*B_p2*XARX11_p2*BXARX_pm3 + 
            144.0*B_p2*XARX10*XARX12*BXARX_pm3 - 
            1296.0*B_p3*XARX01_p2*XARX20*BXARX_pm4 + 
            72.0*B_p2*XARX02*XARX20*BXARX_pm3 + 
            144.0*B_p2*XARX01*XARX21*BXARX_pm3 - 
            6.0*B*XARX22*BXARX_pm2;

        Denom13=31104.0*B_p4*XARX01_p3*XARX10*BXARX_pm5 - 
            3888.0*B_p3*XARX01*XARX02*XARX10*BXARX_pm4  + 
            72.0*B_p2*XARX03*XARX10*BXARX_pm3 - 
            3888.0*B_p3*XARX01_p2*XARX11*BXARX_pm4 + 
            216.0*B_p2*XARX02*XARX11*BXARX_pm3 + 
            216.0*B_p2*XARX01*XARX12*BXARX_pm3 - 
            6.0*B*XARX13*BXARX_pm2;
/*  End of Denom derivatives */

/* d4/ dro4 */
        out = (6.0*Denom20*Nom20)+(4.0*Nom10*Denom30)+(4.0*Denom10*Nom30) + 
            (Nom*Denom40) + (Denom*Nom40);
        ds->df4000 += -B*out*factor;

/*  d4/ dro3 dgro */
        out = 3.0*Nom11*Denom20 + 3.0*Denom11*Nom20 + 3.0*Nom10*Denom21 + 
            3.0*Denom10*Nom21 + Nom01*Denom30 + 
            Denom01*Nom30 + Nom*Denom31 + Denom*Nom31;
        ds->df3010 += -B*out*factor;

/* d4 / dro2 dgro2 */
        out =  4.0*Denom11*Nom11 + 2.0*Nom10*Denom12 + 2.0*Denom10*Nom12 + 
            Nom02*Denom20 + Denom02*Nom20 + 2.0*Nom01*Denom21 + 
            2.0*Denom01*Nom21 + Nom*Denom22 + Denom*Nom22;
        ds->df2020 += -B*out*factor;

/*  d4 / dro dgro3 */
        out = Denom03*Nom10 + 3.0*Nom02*Denom11 + 3.0*Denom02*Nom11 + 
            3*Nom01*Denom12 + 3.0*Denom01*Nom12 + Nom*Denom13;
        ds->df1030 += -B*out*factor;

/*  d4 / dgro4 */
        out = (6.0*Denom02*Nom02)+(4.0*Nom01*Denom03)+(Nom*Denom04);
        ds->df0040 += -B*out*factor;

    } /* endif */

	
    if (dp->rhob > BECKE_THRESHOLD) {

        ro=dp->rhob;
        gro=dp->gradb;


/* Powers of ro, gro */
        iro=1.0/ro;
        ro_pm4f3=pow(ro,(-4.0/3.0));
/* ro_pm4f3=1.0/cbrt(ro*ro*ro*ro); */
        ro_pm7f3=iro*ro_pm4f3;
        ro_pm10f3=iro*ro_pm7f3;
        ro_pm13f3=iro*ro_pm10f3;
        ro_pm16f3=iro*ro_pm13f3;

        gro_p2=gro*gro;
/* End of Powers of ro, gro */


/*  X and derivatives */
        X=gro*(ro_pm4f3);

        X01=(ro_pm4f3);

        X10=-(4.0/3.0)*gro*ro_pm7f3;
        X20=(28.0/9.0)*gro*ro_pm10f3;
        X30=-(280.0/27.0)*gro*ro_pm13f3;
        X40=(3640.0/81.0)*gro*ro_pm16f3;

        X11=-(4.0/3.0)*ro_pm7f3;
        X21=(28.0/9.0)*ro_pm10f3;
        X31=-(280.0/27.0)*ro_pm13f3;

/*  X deriv. powers */
        X_p2=X*X;
        X_p3=X_p2*X;
        X_p4=X_p3*X;

        X01_p2=X01*X01;
        X01_p3=X01_p2*X01;
        X01_p4=X01_p3*X01;

        X10_p2=X10*X10;
        X10_p3=X10_p2*X10;
        X10_p4=X10_p3*X10;

        X20_p2=X20*X20;

        X11_p2=X11*X11;

/*  End of X deriv powers */
/*  End of X derivatives */


        DENX=1.0+X_p2;
        IDENX=1.0/DENX;

        DENX_pm1f2=1.0/sqrt(DENX);
        DENX_pm3f2=IDENX*(DENX_pm1f2);
        DENX_pm5f2=IDENX*(DENX_pm3f2);
        DENX_pm7f2=IDENX*(DENX_pm5f2);

        ASINH_X=asinh(X);
        BXARX=1.0+(6.0*B*X*ASINH_X);
        Denom=1.0/BXARX;

        BXARX_pm2=Denom*Denom;
        BXARX_pm3=BXARX_pm2*Denom;
        BXARX_pm4=BXARX_pm3*Denom;
        BXARX_pm5=BXARX_pm4*Denom;

/*  XARX derivatives */

        XARX01=ASINH_X*X01  + X*X01*DENX_pm1f2;

        XARX02=-X_p2*X01_p2*DENX_pm3f2 + 
            2.0*(X01_p2)* DENX_pm1f2;

        XARX03=3.0*(X_p3)*(X01_p3)* DENX_pm5f2 - 
            4.0*X*(X01_p3)* DENX_pm3f2;

        XARX04=-15.0*(X_p4)*(X01_p4) * DENX_pm7f2 + 
            21.0*(X_p2)*(X01_p4)* DENX_pm5f2 - 
            4.0*(X01_p4) * DENX_pm3f2;

        XARX10=ASINH_X*X10  + X*X10*DENX_pm1f2;

        XARX20=-X_p2*X10_p2*DENX_pm3f2 + 
            2.0*(X10_p2)* DENX_pm1f2 + 
            ASINH_X*X20 + 
            X*X20*DENX_pm1f2;

        XARX30=3.0*(X_p3)*(X10_p3)* DENX_pm5f2 - 
            4.0*X*(X10_p3)* DENX_pm3f2 - 
            3.0*X_p2*X10*X20 * DENX_pm3f2 + 
            6.0*X10*X20*DENX_pm1f2 + 
            ASINH_X*X30 + 
            X*X30* DENX_pm1f2;

        XARX40=-15.0*(X_p4)*(X10_p4) * DENX_pm7f2 + 
            21.0*(X_p2)*(X10_p4)* DENX_pm5f2 - 
            4.0*(X10_p4) * DENX_pm3f2 + 
            18.0*X_p3*X10_p2*X20*DENX_pm5f2 - 
            24.0*X*X10_p2*X20*DENX_pm3f2 - 
            3.0*X_p2*X20_p2*DENX_pm3f2 + 
            6.0*X20_p2*DENX_pm1f2 - 
            4.0*X_p2*X10*X30*DENX_pm3f2 + 
            8.0*X10*X30*DENX_pm1f2 + 
            ASINH_X*X40 + 
            X*X40*DENX_pm1f2;

        XARX11=-X_p2*X01*X10*DENX_pm3f2 + 
            2.0*X01*X10*DENX_pm1f2 + 
            ASINH_X*X11 + 
            X*X11*DENX_pm1f2;

        XARX21=3.0*X_p3*X01*X10_p2*DENX_pm5f2 - 
            4.0*X*X01*X10_p2*DENX_pm3f2 - 
            2.0*X_p2*X10*X11*DENX_pm3f2 + 
            4.0*X10*X11*DENX_pm1f2 - 
            X_p2*X01*X20*DENX_pm3f2 + 
            2.0*X01*X20*DENX_pm1f2 + 
            ASINH_X*X21 + 
            X*X21*DENX_pm1f2;

        XARX31=-15.0*X_p4*X01*X10_p3*DENX_pm7f2 + 
            21.0*X_p2*X01*X10_p3*DENX_pm5f2 - 
            4.0*X01*X10_p3*DENX_pm3f2 + 
            9.0*X_p3*X10_p2*X11*DENX_pm5f2 - 
            12.0*X*X10_p2*X11*DENX_pm3f2 + 
            9.0*X_p3*X01*X10*X20*DENX_pm5f2 - 
            12.0*X*X01*X10*X20*DENX_pm3f2 - 
            3.0*X_p2*X11*X20*DENX_pm3f2 + 
            6.0*X11*X20*DENX_pm1f2 - 
            3.0*X_p2*X10*X21*DENX_pm3f2 + 
            6.0*X10*X21*DENX_pm1f2 - 
            X_p2*X01*X30*DENX_pm3f2 + 
            2.0*X01*X30*DENX_pm1f2 + 
            ASINH_X*X31 + 
            X*X31*DENX_pm1f2;
	
        XARX12=3.0*X_p3*X01_p2*X10*DENX_pm5f2 - 
            4.0*X*X01_p2*X10*DENX_pm3f2 - 
            2.0*X_p2*X01*X11*DENX_pm3f2 + 
            4.0*X01*X11*DENX_pm1f2;

        XARX22=-15.0*X_p4*X01_p2*X10_p2*DENX_pm7f2 + 
            21.0*X_p2*X01_p2*X10_p2*DENX_pm5f2 - 
            4.0*X01_p2*X10_p2*DENX_pm3f2 + 
            12.0*X_p3*X01*X10*X11*DENX_pm5f2 - 
            16.0*X*X01*X10*X11*DENX_pm3f2 - 
            2.0*X_p2*X11_p2*DENX_pm3f2 + 
            4.0*X11_p2*DENX_pm1f2 + 
            3.0*X_p3*X01_p2*X20*DENX_pm5f2 - 
            4.0*X*X01_p2*X20*DENX_pm3f2 - 
            2.0*X_p2*X01*X21*DENX_pm3f2 + 
            4.0*X01*X21*DENX_pm1f2;

        XARX13=-15.0*X_p4*X01_p3*X10*DENX_pm7f2 + 
            21.0*X_p2*X01_p3*X10*DENX_pm5f2 - 
            4.0*X01_p3*X10*DENX_pm3f2 + 
            9.0*X_p3*X01_p2*X11*DENX_pm5f2 - 
            12.0*X*X01_p2*X11*DENX_pm3f2;

/* Xarx der powers */
        XARX10_p2=XARX10*XARX10;
        XARX10_p3=XARX10_p2*XARX10;
        XARX10_p4=XARX10_p3*XARX10;

        XARX20_p2=XARX20*XARX20;

        XARX01_p2=XARX01*XARX01;
        XARX01_p3=XARX01_p2*XARX01;
        XARX01_p4=XARX01_p3*XARX01;

        XARX02_p2=XARX02*XARX02;

        XARX11_p2=XARX11*XARX11;
/* End of Xarx der powers */

/* End of XARX derivatives */

/* Nom derivatives */

        Nom=(gro_p2)*(ro_pm4f3);

        Nom01=2.0*gro*(ro_pm4f3);
        Nom02=2.0*(ro_pm4f3);

        Nom10=-(4.0/3.0)*gro_p2*ro_pm7f3;
        Nom20=(28.0/9.0)*gro_p2*ro_pm10f3;
        Nom30=-(280.0/27.0)*gro_p2*ro_pm13f3;
        Nom40=(3640.0/81.0)*gro_p2*ro_pm16f3;

        Nom11=-(8.0/3.0)*gro*ro_pm7f3;
        Nom21=(56.0/9.0)*gro*ro_pm10f3;
        Nom31=-(560.0/27.0)*gro*ro_pm13f3;

        Nom12=-(8.0/3.0)*ro_pm7f3;
        Nom22=(56.0/9.0)*ro_pm10f3;

/* End of Nom derivatives */

/* Denom derivatives */


        Denom01=-6.0*B*XARX01*BXARX_pm2;

        Denom02=72.0*B_p2*XARX01_p2*BXARX_pm3 - 
            6.0*B*XARX02*BXARX_pm2;

        Denom03=-1296.0*B_p3*XARX01_p3*BXARX_pm4 + 
            216.0*B_p2*XARX01*XARX02*BXARX_pm3 - 
            6.0*B*XARX03*BXARX_pm2;


        Denom04=31104.0*B_p4*XARX01_p4*BXARX_pm5 - 
            7776.0*B_p3*XARX01_p2*XARX02*BXARX_pm4 + 
            216.0*B_p2*XARX02_p2*BXARX_pm3 + 
            288.0*B_p2*XARX01*XARX03*BXARX_pm3 - 
            6.0*B*XARX04*BXARX_pm2;


        Denom10=-6.0*B*XARX10*BXARX_pm2;

        Denom20=72.0*B_p2*XARX10_p2*BXARX_pm3 - 
            6.0*B*XARX20*BXARX_pm2;

        Denom30=-1296.0*B_p3*XARX10_p3*BXARX_pm4 + 
            216.0*B_p2*XARX10*XARX20*BXARX_pm3 - 
            6.0*B*XARX30*BXARX_pm2;


        Denom40=31104.0*B_p4*XARX10_p4*BXARX_pm5 - 
            7776.0*B_p3*XARX10_p2*XARX20*BXARX_pm4 + 
            216.0*B_p2*XARX20_p2*BXARX_pm3 + 
            288.0*B_p2*XARX10*XARX30*BXARX_pm3 - 
            6.0*B*XARX40*BXARX_pm2;

        Denom11=72.0*B_p2*XARX01*XARX10*BXARX_pm3 - 
            6.0*B*XARX11*BXARX_pm2;

        Denom21=-1296.0*B_p3*XARX01*XARX10_p2*BXARX_pm4 + 
            144.0*B_p2*XARX10*XARX11*BXARX_pm3 + 
            72.0*B_p2*XARX01*XARX20*BXARX_pm3 - 
            6.0*B*XARX21*BXARX_pm2;

        Denom31=31104.0*B_p4*XARX01*XARX10_p3*BXARX_pm5 - 
            3888.0*B_p3*XARX10_p2*XARX11*BXARX_pm4 - 
            3888.0*B_p3*XARX01*XARX10*XARX20*BXARX_pm4 + 
            216.0*B_p2*XARX11*XARX20*BXARX_pm3 + 
            216.0*B_p2*XARX10*XARX21*BXARX_pm3 + 
            72.0*B_p2*XARX01*XARX30*BXARX_pm3 - 
            6.0*B*XARX31*BXARX_pm2;

        Denom12=-1296.0*B_p3*XARX01_p2*XARX10*BXARX_pm4 + 
            144.0*B_p2*XARX01*XARX11*BXARX_pm3 + 
            72.0*B_p2*XARX02*XARX10*BXARX_pm3 - 
            6.0*B*XARX12*BXARX_pm2;

        Denom22=31104.0*B_p4*XARX01_p2*XARX10_p2*BXARX_pm5 - 
            1296.0*B_p3*XARX02*XARX10_p2*BXARX_pm4 - 
            5184.0*B_p3*XARX01*XARX10*XARX11*BXARX_pm4 + 
            144.0*B_p2*XARX11_p2*BXARX_pm3 + 
            144.0*B_p2*XARX10*XARX12*BXARX_pm3 - 
            1296.0*B_p3*XARX01_p2*XARX20*BXARX_pm4 + 
            72.0*B_p2*XARX02*XARX20*BXARX_pm3 + 
            144.0*B_p2*XARX01*XARX21*BXARX_pm3 - 
            6.0*B*XARX22*BXARX_pm2;

        Denom13=31104.0*B_p4*XARX01_p3*XARX10*BXARX_pm5 - 
            3888.0*B_p3*XARX01*XARX02*XARX10*BXARX_pm4  + 
            72.0*B_p2*XARX03*XARX10*BXARX_pm3 - 
            3888.0*B_p3*XARX01_p2*XARX11*BXARX_pm4 + 
            216.0*B_p2*XARX02*XARX11*BXARX_pm3 + 
            216.0*B_p2*XARX01*XARX12*BXARX_pm3 - 
            6.0*B*XARX13*BXARX_pm2;

/*  End of Denom derivatives */

/* d4/ dro4 */
        out=(6.0*Denom20*Nom20)+(4.0*Nom10*Denom30)+(4.0*Denom10*Nom30) + 
            (Nom*Denom40) + (Denom*Nom40);
        ds->df0400 += -B*out*factor;

/*  d4/ dro3 dgro */
        out=3.0*Nom11*Denom20 + 3.0*Denom11*Nom20 + 3.0*Nom10*Denom21 + 
            3.0*Denom10*Nom21 + Nom01*Denom30 + 
            Denom01*Nom30 + Nom*Denom31 + Denom*Nom31;
        ds->df0301 += -B*out*factor;

/* d4 / dro2 dgro2 */
        out=  4.0*Denom11*Nom11 + 2.0*Nom10*Denom12 + 2.0*Denom10*Nom12 + 
            Nom02*Denom20 + Denom02*Nom20 + 2.0*Nom01*Denom21 + 
            2.0*Denom01*Nom21 + Nom*Denom22 + Denom*Nom22;
        ds->df0202 += -B*out*factor;

/*  d4 / dro dgro3 */
        out=Denom03*Nom10 + 3.0*Nom02*Denom11 + 3.0*Denom02*Nom11 + 
            3.0*Nom01*Denom12 + 3.0*Denom01*Nom12 + Nom*Denom13;
        ds->df0103 += -B*out*factor;

/*  d4 / dgro4 */
        out=(6.0*Denom02*Nom02)+(4.0*Nom01*Denom03)+(Nom*Denom04);
        ds->df0004 += -B*out*factor;

    } /* endif */

}
