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
/* fun-lypr.c:

   Modified lyp functional for EDF1 functional, with: 
   constants: A = 0.055, B = 0.158, C = 0.25, D = 0.3505;
   Implemented by David Wilson (david.wilson@latrobe.edu.au), Jun 2005.

   Modified from implementation of LYP functional and its derivatives 
   (c) Pawel Salek, pawsa@theochem.kth.se, oct 2001
    Z. Rinkevicius modification for open-shell, general 5 variables formalism.
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
static int lypr_isgga(void) { return 1; }
static int lypr_read(const char* conf_line);
static real lypr_energy(const FunDensProp* dp);
static void lypr_first(FunFirstFuncDrv *ds,   real fac, const FunDensProp* dp);
static void lypr_second(FunSecondFuncDrv *ds, real fac, const FunDensProp* dp);
static void lypr_third(FunThirdFuncDrv *ds,   real fac, const FunDensProp* dp);
static void lypr_fourth(FunFourthFuncDrv *ds, real fac, const FunDensProp* dp);

Functional LYPrFunctional = {
    "LYPr",      /* name */
    lypr_isgga,  /* gga-corrected */
    lypr_read,   /* no extra input expected, just set the common block */
    NULL, /* reporter */
    lypr_energy, 
    lypr_first,
    lypr_second,
    lypr_third,
    lypr_fourth
};

/* IMPLEMENTATION PART */
static int
lypr_read(const char* conf_line)
{
    fun_set_hf_weight(0.0);
    return 1;
}

/**
   The LYP formulas are based on Miehlich et al. article (CPL 157,
   p. 200, 1989).
   The implementation works also for unrestricted case (which is more
   important than you think).
*/
static real
lypr_energy(const FunDensProp* dp)
{
    const real A = 0.055, B = 0.158, C = 0.25, D = 0.3505;
    const real CF = 0.3*pow(3*M_PI*M_PI,2.0/3.0);
    
    real rho, rhom13, denom, omega, delta, ret, ngrad2, ngrada2, ngradb2;
    real rho2, t1, t2, t3, t4, t5, t6;
    real rhoa = dp->rhoa;
    real rhob = dp->rhob;
    real grada = dp->grada;
    real gradb = dp->gradb;

    rho = rhoa+rhob;
    rho2 = rho*rho;
    ngrada2 = grada*grada;
    ngradb2 = gradb*gradb;
    ngrad2  = ngrada2+ngradb2+2*dp->gradab; 
    rhom13 = pow(rho,-1.0/3.0);
    denom = 1+D*rhom13;
    omega = exp(-C*rhom13)/denom*pow(rho,-11.0/3.0);
    delta = rhom13*(C + D/denom);
 
    t1 =  pow(2.0,11.0/3.0)*CF*(pow(rhoa,8.0/3.0) +pow(rhob,8.0/3.0));
    t2 =  (47.0 - 7.0*delta)*ngrad2/18.0;
    t3 = -(2.5 -delta/18.0)*(ngrada2+ngradb2);
    t4 =  (11.0-delta)/9.0*(rhoa*ngrada2 + rhob*ngradb2)/rho;
    t5 = -2.0/3.0*rho2*ngrad2;
    t6 = ((2.0/3.0*rho2-rhoa*rhoa)*ngradb2 +
          (2.0/3.0*rho2-rhob*rhob)*ngrada2);
    ret = -A*(4*rhoa*rhob/(denom*rho)
	      +B*omega*(rhoa*rhob*(t1+t2+t3+t4)+t5+t6)); 

    return ret;
}


/* lypr_first:
   the derivatives of LYP functional wrt alpha density and gradient, suitable
   for unrestricted calculations.
   lyprho : dF/drhoa
   lypgrad: dF/dngrada
   See Doc/dft/functionals.tex for details.
*/
static void
lypr_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
    const real A = 0.055, B = 0.158, C = 0.25, D = 0.3505;
    const real CF = 0.3*pow(3*M_PI*M_PI,2.0/3.0);
    
    real rho, expcr, grad, om, dl, f1, f2, rho2, rhoa2, rhob2;
    real rho13, drho3_2, drho13, grad2, grada2, gradb2, om_1, dl_1;
    real f0_1000, f1_1000, f1_0010, f2_1000, f2_0010, sA, sA10;
    real f0_0100, f1_0100, f1_0001, f2_0100, f2_0001, sA01;
    real f1_00001, f2_00001;
    real rhoa = dp->rhoa;
    real rhob = dp->rhob;
    real grada = dp->grada;
    real gradb = dp->gradb;

    rho = rhoa + rhob;
    rho2 = rho*rho;
    rhoa2 = rhoa*rhoa;
    rhob2 = rhob*rhob;
    grad = grada+gradb;
    rho13  = pow(rho,1.0/3.0);
    drho13 = D+rho13;
    drho3_2= drho13*drho13;
    expcr  = exp(-C/rho13);
    grada2 = grada*grada;
    gradb2 = gradb*gradb;
    grad2  = grada2+gradb2+2*dp->gradab; 
    sA   = (rhoa*grada2+rhob*gradb2)/rho;
    sA10 = (grada-gradb)*grad*rhob/rho2;
    sA01 = (gradb-grada)*grad*rhoa/rho2;
     
    om = expcr*pow(rho,-11.0/3.0)/(1.0+D/rho13); 
    om_1 = (-11*rho*rho13+rho*(C-10*D)+rho/rho13*C*D)/
	(3*pow(rho,16.0/3.0)*drho3_2)*expcr;

    dl = C/rho13 + D/(rho13*(1+D/rho13));
    dl_1 = (-C*drho3_2-D*rho/rho13)/(3*rho*rho13*drho3_2);

    f0_1000 = 4*rhob*(D*rhoa+3*rhob*drho13)*rho13
        /(3*rho*rho*drho3_2);
    f0_0100 = 4*rhoa*(D*rhob+3*rhoa*drho13)*rho13
        /(3*rho*rho*drho3_2);

    f1 = pow(2.0, 11.0/3.0)*CF*(pow(rhoa,8.0/3.0)+pow(rhob,8.0/3.0))
	+ (47-7*dl)*grad2/18.0 +(dl-45)*(grada2+gradb2)/18.0
	  + (11-dl)*sA/9.0;

    f1_1000 = pow(2.0,11.0/3.0)*CF*8.0/3.0*pow(rhoa,5.0/3.0)
    	+ (grada2+gradb2-7*grad2-2*sA)*dl_1/18.0
    	+ (11-dl)*sA10/9.0;
    f1_0100 = pow(2.0,11.0/3.0)*CF*8.0/3.0*pow(rhob,5.0/3.0)
    	+ (grada2+gradb2-7*grad2-2*sA)*dl_1/18.0
    	+ (11-dl)*sA01/9.0;

    f1_0010 = (47-7*dl)*grada/9.0 +
	(dl-45+(22-2*dl)*rhoa/rho)*grada/9.0;
    f1_0001 = (47-7*dl)*gradb/9.0 +
	(dl-45+(22-2*dl)*rhob/rho)*gradb/9.0;
        
    f2 = -2.0/3.0*rho2*grad2
	+(2.0/3.0*rho2-rhoa2)*gradb2+(2.0/3.0*rho2-rhob2)*grada2;
    f2_1000 = -8.0/3.0*rho*dp->gradab - 2*rhoa*gradb2;
    f2_0100 = -8.0/3.0*rho*dp->gradab - 2*rhob*grada2;
    f2_0010 = -2.0*rhob2*grada;
    f2_0001 = -2.0*rhoa2*gradb;

    f1_00001 = (47-7*dl)/9.0;
    f2_00001 = -4.0/3.0*rho2;
   
 
    /* the final section: begin */
    ds->df1000 += factor*(-A*f0_1000-A*B*om_1*(rhoa*rhob*f1+f2)
                          -A*B*om*(rhob*f1+rhoa*rhob*f1_1000+f2_1000)); 
    ds->df0100 += factor*(-A*f0_0100-A*B*om_1*(rhoa*rhob*f1+f2)
                          -A*B*om*(rhoa*f1+rhoa*rhob*f1_0100+f2_0100)); 
    ds->df0010 += factor*(-A*B*om*(rhoa*rhob*f1_0010+f2_0010));
    ds->df0001 += factor*(-A*B*om*(rhoa*rhob*f1_0001+f2_0001));
    ds->df00001+= factor*(-A*B*om*(rhoa*rhob*f1_00001+f2_00001));
}

/* lypr_second:
   second order derivatives of LYP functional with respect to 
   alpha density and gradient.
*/
static void
lypr_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    const real A = 0.055, B = 0.158, C = 0.25, D = 0.3505;
    const real CF = 0.3*pow(3*M_PI*M_PI,2.0/3.0);
    
    real rho, expcr, grad, om, dl, f1, f2, rho2, rhoa2, rhob2, rhoab, de, de1;
    real rho13, drho3_2, drho13, grad2, grada2, gradb2, om_1, dl_1, dl_2;
    real omx, omx1, omx2, om_2, drho3_3;
    real sA, sA1000, sA0100, sA2000, sA0200, sA1100;   
    real f0_1000, f1_1000, f1_0010, f2_1000, f2_0010; 
    real f0_0100, f1_0100, f1_0001, f2_0100, f2_0001;
    real f1_2000, f1_0200, f1_0020, f1_0002;
    real f2_2000, f2_0200, f2_0020, f2_0002;    
    real f1_1100, f1_1010, f1_1001, f1_0101, f1_0110, f1_10001,  f1_01001;
    real f2_1100, f2_1010, f2_1001, f2_0101, f2_0110, f2_10001,  f2_01001;  
    real f1_00001, f2_00001;
    real f0_2000, f0_0200, f0_1100;  
    real rff, rff_1000, rff_0100;
    real rff_2000, rff_0200, rff_0010, rff_0001, rff_00001;
    real rff_1100, rff_1010, rff_0101, rff_10001, rff_01001;

    real rhoa = dp->rhoa;
    real rhob = dp->rhob;
    real grada = dp->grada;
    real gradb = dp->gradb;
     
     
    rho = rhoa + rhob;
    rho2 = rho*rho;
    rhoa2 = rhoa*rhoa;
    rhob2 = rhob*rhob;
    rhoab = rhoa*rhob; 
    grad = grada+gradb;
    rho13  = pow(rho,1.0/3.0);
    drho13 = D+rho13;
    drho3_2= drho13*drho13;
    drho3_3= drho13*drho3_2;
    expcr  = exp(-C/rho13);
    de     = C/(3*rho*rho13); 
    de1    = -4.0*C/(9*rho2*rho13);
    grada2 = grada*grada;
    gradb2 = gradb*gradb;
    grad2  = grada2+gradb2+2*dp->gradab; 
    /* s derivatives */
    sA   = (rhoa*grada2+rhob*gradb2)/rho;
    sA1000 = (grada-gradb)*grad*rhob/rho2;
    sA0100 = (gradb-grada)*grad*rhoa/rho2;
    sA2000 =  -2*sA1000/rho;
    sA0200 =  -2*sA0100/rho;
    sA1100 = (rhoa-rhob)*(grada-gradb)*grad/(rho2*rho);    
    /* omega derivatives */ 
    omx  = pow(rho,-11.0/3.0)/(1.0+D/rho13);
    omx1 = (-11*rho13-10*D)/(3*rho2*rho2*rho13*drho3_2);
    omx2 = 2*(77*rho/rho13+141*D*rho13+65*D*D)/
        (9*pow(rho, 16.0/3.0)*drho3_3);
    om   = expcr*omx;
    om_1 = expcr*(omx1 + omx*de);
    om_2 = expcr*(omx2 + 2*omx1*de + omx*(de1+de*de));    
    /* delta derivatives */
    dl = C/rho13 + D/(rho13*(1+D/rho13));
    dl_1 = (-C*drho3_2-D*rho/rho13)/(3*rho*rho13*drho3_2);
    dl_2 = (4*rho*(C+D) + 2*D*(6*rho13*C*D +2*C*D*D + rho/rho13*(6*C+D)))/
        (9*rho2*rho13*drho3_3);
    /* f0 derivatives */
    f0_1000 = 4*rhob*(D*rhoa+3*rhob*drho13)*rho13
        /(3*rho*rho*drho3_2);
    f0_0100 = 4*rhoa*(D*rhob+3*rhoa*drho13)*rho13
        /(3*rho*rho*drho3_2);
    f0_2000 = -8*rhob*(rhoa*D*(D+2*rho13)+3*rhob*drho13*(2*D+3*rho13))/
        (9*rho2*rho/rho13*drho3_2*drho13);
    f0_0200 = -8*rhoa*(rhob*D*(D+2*rho13)+3*rhoa*drho13*(2*D+3*rho13))/
        (9*rho2*rho/rho13*drho3_2*drho13);
    f0_1100 = 4*(18*rhoa*rhob*rho/rho13 + rho13*(3*rhoa2+32*rhoab+3*rhob2)*D +
              (3*rhoa2+16*rhoab+3*rhob2)*D*D)/(9*rho2*rho/rho13*drho3_3);
    /* f1 derivatives */
    f1 = pow(2.0, 11.0/3.0)*CF*(pow(rhoa,8.0/3.0)+pow(rhob,8.0/3.0))
	+ (47-7*dl)*grad2/18.0 +(dl-45)*(grada2+gradb2)/18.0
	  + (11-dl)*sA/9.0;
    f1_1000 = pow(2.0,11.0/3.0)*CF*8.0/3.0*pow(rhoa,5.0/3.0)
    	+ (grada2+gradb2-7*grad2-2*sA)*dl_1/18.0
    	+ (11-dl)*sA1000/9.0;
    f1_0100 = pow(2.0,11.0/3.0)*CF*8.0/3.0*pow(rhob,5.0/3.0)
    	+ (grada2+gradb2-7*grad2-2*sA)*dl_1/18.0
    	+ (11-dl)*sA0100/9.0;
    f1_0010 = (47-7*dl)*grada/9.0 +
	(dl-45+(22-2*dl)*rhoa/rho)*grada/9.0;
    f1_0001 = (47-7*dl)*gradb/9.0 +
	(dl-45+(22-2*dl)*rhob/rho)*gradb/9.0;
    f1_00001 = (47-7*dl)/9.0;    
    f1_2000 = pow(2.0,11.0/3.0)*CF*40.0/9.0*pow(rhoa,2.0/3.0)
        - 2*sA1000*dl_1/9.0
        + (grada2+gradb2-7*grad2-2*sA)*dl_2/18.0
        + (11-dl)*sA2000/9.0;
    f1_0200 = pow(2.0,11.0/3.0)*CF*40.0/9.0*pow(rhob,2.0/3.0)
        - 2*sA0100*dl_1/9.0
        + (grada2+gradb2-7*grad2-2*sA)*dl_2/18.0
        + (11-dl)*sA0200/9.0;
    f1_0020 = (47-7*dl)/9.0+(dl-45+(22-2*dl)*rhoa/rho)/9.0; 
    f1_0002 = (47-7*dl)/9.0+(dl-45+(22-2*dl)*rhob/rho)/9.0;
    f1_1100 = -2*sA0100*dl_1/18.0+(grada2+gradb2-7*grad2-2*sA)*dl_2/18.0
              -dl_1*sA1000/9.0+(11-dl)*sA1100/9.0;
    f1_1010 = (grada*(1-2*rhoa/rho)-7.0*grada)*dl_1/9.0 +
              (11-dl)*rhob/rho2*grada/4.5;
    f1_0101 = (gradb*(1-2*rhob/rho)-7.0*gradb)*dl_1/9.0 +
              (11-dl)*rhoa/rho2*gradb/4.5;
    f1_1001 = (gradb*(1-2*rhob/rho)-7.0*gradb)*dl_1/9.0 -
              (11-dl)*rhob/rho2*gradb/4.5;    
    f1_0110 = (grada*(1-2*rhoa/rho)-7.0*grada)*dl_1/9.0 -
              (11-dl)*rhoa/rho2*grada/4.5;  
    f1_10001 = -7.0*dl_1/9.0;
    f1_01001 = -7.0*dl_1/9.0; 
    /* f2 derivatives */
    f2 = -2.0/3.0*rho2*grad2
	+(2.0/3.0*rho2-rhoa2)*gradb2+(2.0/3.0*rho2-rhob2)*grada2;
    f2_1000  = -8.0/3.0*rho*dp->gradab - 2.0*rhoa*gradb2;
    f2_0100  = -8.0/3.0*rho*dp->gradab - 2.0*rhob*grada2;
    f2_0010  = -2.0*rhob2*grada;
    f2_0001  = -2.0*rhoa2*gradb;
    f2_00001 = -4.0/3.0*rho2;
    f2_2000  = -8.0/3.0*dp->gradab - 2.0*gradb2;
    f2_0200  = -8.0/3.0*dp->gradab - 2.0*grada2;
    f2_0020  = -2.0*rhob2;
    f2_0002  = -2.0*rhoa2;  
    f2_1100  = -8.0/3.0*dp->gradab; 
    f2_1010  =  0.0;
    f2_0101  =  0.0;
    f2_1001  = -4.0*rhoa*gradb; 
    f2_0110  = -4.0*rhob*grada; 
    f2_10001  = -8.0/3.0*rho;
    f2_01001  = -8.0/3.0*rho;
    /* derivatives sums */
    rff      = rhoa*rhob*f1+f2;
    rff_1000 = rhob*f1+rhoa*rhob*f1_1000+f2_1000;
    rff_0100 = rhoa*f1+rhoa*rhob*f1_0100+f2_0100;
    rff_2000 = 2*rhob*f1_1000+rhoa*rhob*f1_2000+f2_2000;
    rff_0200 = 2*rhoa*f1_0100+rhoa*rhob*f1_0200+f2_0200;
    rff_1100 = f1+rhob*f1_0100+rhoa*f1_1000+rhoa*rhob*f1_1100+f2_1100;
    rff_0010 = rhoa*rhob*f1_0010+f2_0010;
    rff_0001 = rhoa*rhob*f1_0001+f2_0001;
    rff_1010 = rhob*f1_0010+rhoa*rhob*f1_1010+f2_1010;
    rff_0101 = rhoa*f1_0001+rhoa*rhob*f1_0101+f2_0101;
    /* derivatives sum with respect grada*gradb */
    rff_00001 = rhoa*rhob*f1_00001+f2_00001;
    rff_10001 = rhob*f1_00001+rhoa*rhob*f1_10001+f2_10001;
    rff_01001 = rhoa*f1_00001+rhoa*rhob*f1_01001+f2_01001;

 
    /* the final section: first derivatives */
    ds->df1000  += factor*(-A*f0_1000-A*B*om_1*(rhoa*rhob*f1+f2)
                           -A*B*om*(rhob*f1+rhoa*rhob*f1_1000+f2_1000)); 
    ds->df0100  += factor*(-A*f0_0100-A*B*om_1*(rhoa*rhob*f1+f2)
                           -A*B*om*(rhoa*f1+rhoa*rhob*f1_0100+f2_0100)); 
    ds->df0010  += factor*(-A*B*om*(rhoa*rhob*f1_0010+f2_0010));
    ds->df0001  += factor*(-A*B*om*(rhoa*rhob*f1_0001+f2_0001));
    ds->df00001 += factor*(-A*B*om*(rhoa*rhob*f1_00001+f2_00001));
   /* the final section: second derivatives */
    ds->df2000 += factor*(-A*f0_2000
                          -A*B*(om_2*rff+2*om_1*rff_1000+om*rff_2000));
    ds->df0200 += factor*(-A*f0_0200
                          -A*B*(om_2*rff+2*om_1*rff_0100+om*rff_0200));
    ds->df0020 += factor*(-A*B*om*(rhoa*rhob*f1_0020+f2_0020));
    ds->df0002 += factor*(-A*B*om*(rhoa*rhob*f1_0002+f2_0002));
    /* the mixed derivatives */ 
    ds->df1100 += factor*(-A*f0_1100-A*B*(om_2*rff+om_1*rff_0100+
			   om_1*rff_1000 + om*rff_1100));
    ds->df1010 += factor*(-A*B*(om_1*rff_0010+om*rff_1010));
    ds->df1001 += factor*(-A*B*(om_1*(rhoa*rhob*f1_0001+f2_0001) +
                           om*(rhob*f1_0001+rhoa*rhob*f1_1001+f2_1001)));
    ds->df0101 += factor*(-A*B*(om_1*rff_0001+om*rff_0101));
    ds->df0110 += factor*(-A*B*(om_1*(rhoa*rhob*f1_0010+f2_0010) +
                           om*(rhoa*f1_0010+rhoa*rhob*f1_0110+f2_0110)));
    ds->df10001 += factor*(-A*B*(om_1*rff_00001+om*rff_10001));
    ds->df01001 += factor*(-A*B*(om_1*rff_00001+om*rff_01001));
}

static void
lypr_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    const real A = 0.055, B = 0.158, C = 0.25, D = 0.3505;
    const real CF = 0.3*pow(3*M_PI*M_PI,2.0/3.0);
    
    real rho, expcr, grad, om, dl, f1, f2, rho2, rhoa2, rhob2, rhoab;
    real de, de1, de2;
    real rho13, drho3_2, drho13, grad2, grada2, gradb2, om_1, dl_1, dl_2, dl_3;
    real omx, omx1, omx2,omx3, om_2,om_3, drho3_3;
    real sA, sA1000, sA0100, sA2000, sA0200, sA1100, sA3000, sA0300, sA2100;
    real sA1200;   
    real f0_1000, f1_1000, f1_0010, f2_1000, f2_0010; 
    real f0_0100, f1_0100, f1_0001, f2_0100, f2_0001;
    real f1_2000, f1_0200, f1_0020, f1_0002;
    real f2_2000, f2_0200, f2_0020, f2_0002;    
    real f1_1100, f1_1010, f1_1001, f1_0101, f1_0110, f1_10001,  f1_01001;
    real f2_1100, f2_1010, f2_1001, f2_0101, f2_0110, f2_10001,  f2_01001;  
    real f1_00001, f2_00001;
    real f0_2000, f0_0200, f0_1100, f0_3000, f0_0300, f0_2100, f0_1200;
    real f1_3000, f1_0300, f1_2100, f1_1200;  
    real f1_2010, f1_2001, f1_0201, f1_0210, f1_1020, f1_0102,f1_1002, f1_0120;
    real f1_1110, f1_1101, f1_20001, f1_02001, f1_11001; 
    real f2_2001, f2_0210, f2_1002, f2_0120, f2_20001, f2_02001, f2_11001;
    real rff, rff_1000, rff_0100;
    real rff_2000, rff_0200, rff_0010, rff_0001, rff_00001;
    real rff_1100, rff_1010, rff_0101, rff_10001,rff_01001, rff_20001, rff_02001;
    real rff_3000, rff_0300, rff_2100, rff_1200, rff_1110, rff_1101, rff_11001;
    real rff_2010, rff_2001, rff_1001, rff_0201, rff_0210, rff_0110; 

    real rhoa = dp->rhoa;
    real rhob = dp->rhob;
    real grada = dp->grada;
    real gradb = dp->gradb;
     
     
    rho = rhoa + rhob;
    rho2 = rho*rho;
    rhoa2 = rhoa*rhoa;
    rhob2 = rhob*rhob;
    rhoab = rhoa*rhob; 
    grad = grada+gradb;
    rho13  = pow(rho,1.0/3.0);
    drho13 = D+rho13;
    drho3_2= drho13*drho13;
    drho3_3= drho13*drho3_2;
    expcr  = exp(-C/rho13);
    de     = C/(3*rho*rho13); 
    de1    = -4.0*C/(9*rho2*rho13);
    de2    = -de1*7.0/(3*rho);
    grada2 = grada*grada;
    gradb2 = gradb*gradb;
    grad2  = grada2+gradb2+2*dp->gradab; 
    /* s derivatives */
    sA   = (rhoa*grada2+rhob*gradb2)/rho;
    sA1000 = (grada-gradb)*grad*rhob/rho2;
    sA0100 = (gradb-grada)*grad*rhoa/rho2;
    sA2000 =  -2*sA1000/rho;
    sA0200 =  -2*sA0100/rho;
    sA1100 = (rhoa-rhob)*(grada-gradb)*grad/(rho2*rho);    
    sA2100 = -2*sA1100/rho+2*sA1000/rho2;
    sA1200 = -2*sA1100/rho+2*sA0100/rho2;
    sA3000 = -3*sA2000/rho;
    sA0300 = -3*sA0200/rho; 
    /* omega derivatives */ 
    omx  = pow(rho,-11.0/3.0)/(1.0+D/rho13);
    omx1 = (-11*rho13-10*D)/(3*rho2*rho2*rho13*drho3_2);
    omx2 = 2*(77*rho/rho13+141*D*rho13+65*D*D)/
        (9*pow(rho, 16.0/3.0)*drho3_3);
    omx3 = -2*(1309*rho+3616*rho/rho13*D+3350*rho13*D*D+1040*D*D*D)/
        (27*pow(rho,19.0/3.0)*(drho3_2*drho3_2));
    om   = expcr*omx;
    om_1 = expcr*(omx1 + omx*de);
    om_2 = expcr*(omx2 + 2*omx1*de + omx*(de1+de*de));    
    om_3 = expcr*(omx3+3*omx2*de+3*omx1*(de1+de*de)
                     +omx*(de2 + 3*de1*de+ de*de*de));
    /* delta derivatives */
    dl = C/rho13 + D/(rho13*(1+D/rho13));
    dl_1 = (-C*drho3_2-D*rho/rho13)/(3*rho*rho13*drho3_2);
    dl_2 = (4*rho*(C+D) + 2*D*(6*rho13*C*D +2*C*D*D + rho/rho13*(6*C+D)))/
        (9*rho2*rho13*drho3_3);
    dl_3 = (-2*rho/rho13*D*D*(84*C+5*D)-4*D*(rho*(28*C+8*D)+7*C*D*D*D)
               -28*rho13*(4*C*D*D*D+rho*(C+D)))/
        (27*rho2*rho*rho13*drho3_2*drho3_2);
    /* f0 derivatives */
    f0_1000 = 4*rhob*(D*rhoa+3*rhob*drho13)*rho13
        /(3*rho*rho*drho3_2);
    f0_0100 = 4*rhoa*(D*rhob+3*rhoa*drho13)*rho13
        /(3*rho*rho*drho3_2);
    f0_2000 = -8*rhob*(rhoa*D*(D+2*rho13)+3*rhob*drho13*(2*D+3*rho13))/
        (9*rho2*rho/rho13*drho3_2*drho13);
    f0_0200 = -8*rhoa*(rhob*D*(D+2*rho13)+3*rhoa*drho13*(2*D+3*rho13))/
        (9*rho2*rho/rho13*drho3_2*drho13);
    f0_1100 = 4*(18*rhoa*rhob*rho/rho13 + rho13*(3*rhoa2+32*rhoab+3*rhob2)*D +
              (3*rhoa2+16*rhoab+3*rhob2)*D*D)/(9*rho2*rho/rho13*drho3_3);
    f0_3000 = 8*rhob*(81*rhoa*rhob+81*rhob2+2*rho/rho13*(7*rhoa+99*rhob)*D
                      +2*rho13*(8*rhoa+81*rhob)*D*D+(5*rhoa+45*rhob)*D*D*D)/
        (27*pow(rho,11.0/3.0)*drho3_2*drho3_2);
    f0_0300 = 8*rhoa*(81*rhob*rhoa+81*rhoa2+2*rho/rho13*(7*rhob+99*rhoa)*D
                      +2*rho13*(8*rhob+81*rhoa)*D*D+(5*rhob+45*rhoa)*D*D*D)/
        (27*pow(rho,11.0/3.0)*drho3_2*drho3_2);
    f0_2100 = 8*(-54*rhoa2*rhob-27*rhoa*rhob2+27*rhob2*rhob
                -2*rho/rho13*(3*rhoa2+65*rhoab-30*rhob2)*D
                -rho13*(9*rhoa2+110*rhoab-45*rhob2)*D*D
                +(12*rhob2-3*rhoa2-31*rhoab)*D*D*D)/
        (27*rho2*rho2/rho13*drho3_3*drho13);
    f0_1200 = 8*(-54*rhob2*rhoa-27*rhob*rhoa2+27*rhoa2*rhoa
                -2*rho/rho13*(3*rhob2+65*rhoab-30*rhoa2)*D
                -rho13*(9*rhob2+110*rhoab-45*rhoa2)*D*D
                +(12*rhoa2-3*rhob2-31*rhoab)*D*D*D)/
        (27*rho2*rho2/rho13*drho3_3*drho13);
    /* f1 derivatives */
    f1 = pow(2.0, 11.0/3.0)*CF*(pow(rhoa,8.0/3.0)+pow(rhob,8.0/3.0))
	+ (47-7*dl)*grad2/18.0 +(dl-45)*(grada2+gradb2)/18.0
	  + (11-dl)*sA/9.0;
    f1_1000 = pow(2.0,11.0/3.0)*CF*8.0/3.0*pow(rhoa,5.0/3.0)
    	+ (grada2+gradb2-7*grad2-2*sA)*dl_1/18.0
    	+ (11-dl)*sA1000/9.0;
    f1_0100 = pow(2.0,11.0/3.0)*CF*8.0/3.0*pow(rhob,5.0/3.0)
    	+ (grada2+gradb2-7*grad2-2*sA)*dl_1/18.0
    	+ (11-dl)*sA0100/9.0;
    f1_0010 = (47-7*dl)*grada/9.0 +
	(dl-45+(22-2*dl)*rhoa/rho)*grada/9.0;
    f1_0001 = (47-7*dl)*gradb/9.0 +
	(dl-45+(22-2*dl)*rhob/rho)*gradb/9.0;
    f1_00001 = (47-7*dl)/9.0;    
    f1_2000 = pow(2.0,11.0/3.0)*CF*40.0/9.0*pow(rhoa,2.0/3.0)
        - 2*sA1000*dl_1/9.0
        + (grada2+gradb2-7*grad2-2*sA)*dl_2/18.0
        + (11-dl)*sA2000/9.0;
    f1_0200 = pow(2.0,11.0/3.0)*CF*40.0/9.0*pow(rhob,2.0/3.0)
        - 2*sA0100*dl_1/9.0
        + (grada2+gradb2-7*grad2-2*sA)*dl_2/18.0
        + (11-dl)*sA0200/9.0;
    f1_0020 = (47-7*dl)/9.0+(dl-45+(22-2*dl)*rhoa/rho)/9.0; 
    f1_0002 = (47-7*dl)/9.0+(dl-45+(22-2*dl)*rhob/rho)/9.0;
    f1_1100 = -sA0100*dl_1/9.0+(grada2+gradb2-7*grad2-2*sA)*dl_2/18.0
              -dl_1*sA1000/9.0+(11-dl)*sA1100/9.0;
    f1_1010 = (grada*(1-2*rhoa/rho)-7.0*grada)*dl_1/9.0 +
              (11-dl)*rhob/rho2*grada/4.5;
    f1_0101 = (gradb*(1-2*rhob/rho)-7.0*gradb)*dl_1/9.0 +
              (11-dl)*rhoa/rho2*gradb/4.5;
    f1_1001 = (gradb*(1-2*rhob/rho)-7.0*gradb)*dl_1/9.0 -
              (11-dl)*rhob/rho2*gradb/4.5;    
    f1_0110 = (grada*(1-2*rhoa/rho)-7.0*grada)*dl_1/9.0 -
              (11-dl)*rhoa/rho2*grada/4.5;  
    f1_10001 = -7.0*dl_1/9.0;
    f1_01001 = -7.0*dl_1/9.0; 
    f1_3000 = pow(2.0,11.0/3.0)*CF*80.0/27.0*pow(rhoa,-1.0/3.0)
        - 2*sA2000*dl_1/9.0 - 2*sA1000*dl_2/9.0
        + (grada2+gradb2-7*grad2-2*sA)*dl_3/18.0
        - sA1000*dl_2/9.0
        -dl_1*sA2000/9.0 + (11-dl)*sA3000/9.0 ;
    f1_0300 = pow(2.0,11.0/3.0)*CF*80.0/27.0*pow(rhob,-1.0/3.0)
        - 2*sA0200*dl_1/9.0 - 2*sA0100*dl_2/9.0
        + (grada2+gradb2-7*grad2-2*sA)*dl_3/18.0
        - sA0100*dl_2/9.0
        -dl_1*sA0200/9.0 + (11-dl)*sA0300/9.0 ;
    f1_2100 = - 2*(sA1100*dl_1+sA1000*dl_2)/9.0
        + (grada2+gradb2-7*grad2-2*sA)*dl_3/18.0-sA0100*dl_2/9.0
        + (11-dl)*sA2100/9.0-dl_1*sA2000/9.0;
    f1_1200 = - 2*(sA1100*dl_1+sA0100*dl_2)/9.0
        + (grada2+gradb2-7*grad2-2*sA)*dl_3/18.0-sA1000*dl_2/9.0
     + (11-dl)*sA1200/9.0-dl_1*sA0200/9.0;
    f1_2010 = (grada-7.0*grada)*dl_2/9.0
        +(-dl_1*rhob/rho2-2*(11-dl)*rhob/(rho2*rho)
          -dl_2*rhoa/rho-dl_1*rhob/rho2)*grada/4.5;
    f1_0201 = (gradb-7.0*gradb)*dl_2/9.0
        +(-dl_1*rhoa/rho2-2*(11-dl)*rhoa/(rho2*rho)
          -dl_2*rhob/rho-dl_1*rhoa/rho2)*gradb/4.5;
    f1_2001 =  2*gradb*rhob/rho2*dl_1/4.5
        + (gradb-7.0*gradb-2.0*gradb*rhob/rho)*dl_2/9.0
        + (11-dl)*2.0*gradb*rhob/(4.5*rho2*rho);
    f1_0210 =  2*grada*rhoa/rho2*dl_1/4.5
        + (grada-7.0*grada-2.0*grada*rhoa/rho)*dl_2/9.0
        + (11-dl)*2.0*grada*rhoa/(4.5*rho2*rho);
    f1_1020 = ((1-2*rhoa/rho)-7.0)*dl_1/9.0 +
              (11-dl)*rhob/rho2/4.5;
    f1_0102 = ((1-2*rhob/rho)-7.0)*dl_1/9.0 +
              (11-dl)*rhoa/rho2/4.5;
    f1_0120 = ((1-2*rhoa/rho)-7.0)*dl_1/9.0 -
              (11-dl)*rhoa/rho2/4.5;  
    f1_1002 = ((1-2*rhob/rho)-7.0)*dl_1/9.0 -
              (11-dl)*rhob/rho2/4.5;  
    f1_1110= (2*grada*(rhoa-rhob)*(11-dl+rho*dl_1)-
              rho2*(2*(4*rhoa+3*rhob)*grada)*dl_2)/(9*rho2*rho);
    f1_1101= (2*gradb*(rhob-rhoa)*(11-dl+rho*dl_1)-
              rho2*(2*(4*rhob+3*rhoa)*gradb)*dl_2)/(9*rho2*rho);
    f1_20001 = -7.0*dl_2/9.0;
    f1_02001 = -7.0*dl_2/9.0;    
    f1_11001 = -7.0*dl_2/9.0; 
    /* f2 derivatives */
    f2 = -2.0/3.0*rho2*grad2
	+(2.0/3.0*rho2-rhoa2)*gradb2+(2.0/3.0*rho2-rhob2)*grada2;
    f2_1000  = -8.0/3.0*rho*dp->gradab - 2.0*rhoa*gradb2;
    f2_0100  = -8.0/3.0*rho*dp->gradab - 2.0*rhob*grada2;
    f2_0010  = -2.0*rhob2*grada;
    f2_0001  = -2.0*rhoa2*gradb;
    f2_00001 = -4.0/3.0*rho2;
    f2_2000  = -8.0/3.0*dp->gradab - 2.0*gradb2;
    f2_0200  = -8.0/3.0*dp->gradab - 2.0*grada2;
    f2_0020  = -2.0*rhob2;
    f2_0002  = -2.0*rhoa2;  
    f2_1100  = -8.0/3.0*dp->gradab; 
    f2_1010  =  0.0;
    f2_0101  =  0.0;
    f2_1001  = -4.0*rhoa*gradb; 
    f2_0110  = -4.0*rhob*grada; 
    f2_10001 = -8.0/3.0*rho;
    f2_01001 = -8.0/3.0*rho;
    f2_2001  = -4.0*gradb;
    f2_0210  = -4.0*grada;
    f2_0120  = -4.0*rhob;
    f2_1002  = -4.0*rhoa;
    f2_20001 = -8.0/3.0;
    f2_02001 = -8.0/3.0;
    f2_11001 = -8.0/3.0;
    /* derivatives sums */
    rff      = rhoa*rhob*f1+f2;
    rff_1000 = rhob*f1+rhoa*rhob*f1_1000+f2_1000;
    rff_0100 = rhoa*f1+rhoa*rhob*f1_0100+f2_0100;
    rff_2000 = 2*rhob*f1_1000+rhoa*rhob*f1_2000+f2_2000;
    rff_0200 = 2*rhoa*f1_0100+rhoa*rhob*f1_0200+f2_0200;
    rff_1100 = f1+rhob*f1_0100+rhoa*f1_1000+rhoa*rhob*f1_1100+f2_1100;
    rff_0010 = rhoa*rhob*f1_0010+f2_0010;
    rff_0001 = rhoa*rhob*f1_0001+f2_0001;
    rff_1010 = rhob*f1_0010+rhoa*rhob*f1_1010+f2_1010;
    rff_0101 = rhoa*f1_0001+rhoa*rhob*f1_0101+f2_0101;
    /* derivatives sum with respect grada*gradb */
    rff_00001 = rhoa*rhob*f1_00001+f2_00001;
    rff_10001 = rhob*f1_00001+rhoa*rhob*f1_10001+f2_10001;
    rff_01001 = rhoa*f1_00001+rhoa*rhob*f1_01001+f2_01001;
    rff_20001 = 2*rhob*f1_10001+rhoab*f1_20001+f2_20001;
    rff_02001 = 2*rhoa*f1_01001+rhoab*f1_02001+f2_02001;
    rff_11001 = f1_00001+rhob*f1_01001+rhoa*f1_10001+rhoab*f1_11001+f2_11001;
    /* THIRD DERIVATIVES */
    rff_3000 = 2*rhob*f1_2000+rhob*(f1_2000+rhoa*f1_3000); 
    rff_0300 = 2*rhoa*f1_0200+rhoa*(f1_0200+rhob*f1_0300);
    rff_2100 = 2*(f1_1000+rhob*f1_1100)+rhoa*f1_2000+rhoa*rhob*f1_2100;
    rff_1200 = 2*(f1_0100+rhoa*f1_1100)+rhob*f1_0200+rhoa*rhob*f1_1200;
    rff_2010 = 2*rhob*f1_1010+rhoab*f1_2010;
    rff_1001 = rhob*f1_0001+rhoab*f1_1001+f2_1001;
    rff_2001 = 2*rhob*f1_1001+rhoab*f1_2001+f2_2001;
    rff_0201 = 2*rhoa*f1_0101+rhoab*f1_0201;
    rff_0110 = rhoa*f1_0010+rhoab*f1_0110+f2_0110;
    rff_0210 = 2*rhoa*f1_0110+rhoab*f1_0210+f2_0210;
    rff_1110 = f1_0010+rhob*f1_0110+rhoa*f1_1010+rhoab*f1_1110;
    rff_1101 = f1_0001+rhob*f1_0101+rhoa*f1_1001+rhoab*f1_1101;


 
    /* the final section: first derivatives */
    ds->df1000  += factor*(-A*f0_1000-A*B*om_1*(rhoa*rhob*f1+f2)
                           -A*B*om*(rhob*f1+rhoa*rhob*f1_1000+f2_1000)); 
    ds->df0100  += factor*(-A*f0_0100-A*B*om_1*(rhoa*rhob*f1+f2)
                           -A*B*om*(rhoa*f1+rhoa*rhob*f1_0100+f2_0100)); 
    ds->df0010  += factor*(-A*B*om*(rhoa*rhob*f1_0010+f2_0010));
    ds->df0001  += factor*(-A*B*om*(rhoa*rhob*f1_0001+f2_0001));
    ds->df00001 += factor*(-A*B*om*(rhoa*rhob*f1_00001+f2_00001));
   /* the final section: second derivatives */
    ds->df2000 += factor*(-A*f0_2000
                          -A*B*(om_2*rff+2*om_1*rff_1000+om*rff_2000));
    ds->df0200 += factor*(-A*f0_0200
                          -A*B*(om_2*rff+2*om_1*rff_0100+om*rff_0200));
    ds->df0020 += factor*(-A*B*om*(rhoa*rhob*f1_0020+f2_0020));
    ds->df0002 += factor*(-A*B*om*(rhoa*rhob*f1_0002+f2_0002));
    /* the mixed derivatives */ 
    ds->df1100 += factor*(-A*f0_1100-A*B*(om_2*rff+om_1*rff_0100+
			   om_1*rff_1000 + om*rff_1100));
    ds->df1010 += factor*(-A*B*(om_1*rff_0010+om*rff_1010));
    ds->df1001 += factor*(-A*B*(om_1*(rhoa*rhob*f1_0001+f2_0001) +
                           om*(rhob*f1_0001+rhoa*rhob*f1_1001+f2_1001)));
    ds->df0101 += factor*(-A*B*(om_1*rff_0001+om*rff_0101));
    ds->df0110 += factor*(-A*B*(om_1*(rhoa*rhob*f1_0010+f2_0010) +
                           om*(rhoa*f1_0010+rhoa*rhob*f1_0110+f2_0110)));
    ds->df10001 += factor*(-A*B*(om_1*rff_00001+om*rff_10001));
    ds->df01001 += factor*(-A*B*(om_1*rff_00001+om*rff_01001));
    /* third order derivative: density only dependent */
    ds->df3000 += factor*(-A*f0_3000
                          -A*B*(om_3*rff+3*om_2*rff_1000 +
                                3*om_1*rff_2000+om*rff_3000));
    ds->df0300 += factor*(-A*f0_0300
                          -A*B*(om_3*rff+3*om_2*rff_0100 +
                                3*om_1*rff_0200+om*rff_0300));
    ds->df2100 += factor*(-A*f0_2100
			  -A*B*(om_3*rff+om_2*(2*rff_1000+rff_0100)+
				om_1*(rff_2000+2*rff_1100)+om*rff_2100));
    ds->df1200 += factor*(-A*f0_1200
			  -A*B*(om_3*rff+om_2*(2*rff_0100+rff_1000)+
				om_1*(rff_0200+2*rff_1100)+om*rff_1200));
    /* third order derivative: mixed  */    
    ds->df2010 += factor*(-A*B*(om_2*rff_0010+2*om_1*rff_1010+om*rff_2010));     
    ds->df2001 += factor*(-A*B*(om_2*rff_0001+2*om_1*rff_1001+om*rff_2001));
    ds->df0201 += factor*(-A*B*(om_2*rff_0001+2*om_1*rff_0101+om*rff_0201));
    ds->df0210 += factor*(-A*B*(om_2*rff_0010+2*om_1*rff_0110+om*rff_0210));
    ds->df1020 += factor*(-A*B*(om_1*(rhoa*rhob*f1_0020+f2_0020)+
                               om*rhob*(f1_0020+rhoa*f1_1020)));
    ds->df0102 += factor*(-A*B*(om_1*(rhoa*rhob*f1_0002+f2_0002)+
                               om*rhoa*(f1_0002+rhob*f1_0102)));
    ds->df1002 += factor*(-A*B*(om_1*(rhoab*f1_0002+f2_0002)+
                               om*(rhob*f1_0002+rhoab*f1_1002+f2_1002)));
    ds->df0120 += factor*(-A*B*(om_1*(rhoab*f1_0020+f2_0020)+
                               om*(rhoa*f1_0020+rhoab*f1_0120+f2_0120)));
    ds->df1110 += factor*(-A*B*(om_2*rff_0010+om_1*rff_0110+om_1*rff_1010 +
                               om*rff_1110));
    ds->df1101 += factor*(-A*B*(om_2*rff_0001+om_1*rff_0101+om_1*rff_1001 +
                               om*rff_1101));
    ds->df20001 += factor*(-A*B*(om_2*rff_00001+2*om_1*rff_10001+om*rff_20001)); 
    ds->df02001 += factor*(-A*B*(om_2*rff_00001+2*om_1*rff_01001+om*rff_02001)); 
    ds->df11001 += factor*(-A*B*(om_2*rff_00001+om_1*rff_01001+om_1*rff_10001 +
                               om*rff_11001));

}

static 
void lypr_fourth(FunFourthFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real  ComDenom, ComDenom_p2, ComDenom_p3, ComDenom_p4, ComDenom_p5;
    real  Delta, Delta_01, Delta_02, Delta_03, Delta_04, Delta_10;
    real  Delta_11, Delta_12, Delta_13, Delta_20, Delta_21, Delta_22;
    real  Delta_30, Delta_31, Delta_40, EXPC;
    real  F0_04, F0_11, F0_12, F0_13, F0_21, F0_22, F0_31, F0_40;
    real  F1, F1_00001, F1_00010, F1_00020, F1_00100, F1_00200;
    real  F1_01000, F1_01001, F1_01010, F1_01020, F1_01100, F1_01200;
    real  F1_02000, F1_02001, F1_02010, F1_02020, F1_02100, F1_02200;
    real  F1_03000, F1_03001, F1_03010, F1_03100, F1_04000, F1_10000;
    real  F1_10001, F1_10010, F1_10020, F1_10100, F1_10200, F1_11000;
    real  F1_11001, F1_11010, F1_11020, F1_11100, F1_11200, F1_12000;
    real  F1_12001, F1_12010, F1_12100, F1_13000, F1_20000, F1_20001;
    real  F1_20010, F1_20020, F1_20100, F1_20200, F1_21000, F1_21001;
    real  F1_21010, F1_21100, F1_22000, F1_30000, F1_30001, F1_30010;
    real  F1_30100, F1_31000, F1_40000, F1term;
    real  F2, F2_00001, F2_00010, F2_00020, F2_00100, F2_00200;
    real  F2_01000, F2_01001, F2_01100, F2_01200, F2_02000, F2_02001;
    real  F2_02100, F2_02200, F2_10000, F2_10001, F2_10010, F2_10020;
    real  F2_11000, F2_11001, F2_20000, F2_20001, F2_20010, F2_20020;
    real  S, S_0001, S_0002, S_0010, S_0020, S_0100;
    real  S_0101, S_0102, S_0110, S_0120, S_0200, S_0201;
    real  S_0202, S_0210, S_0220, S_0300, S_0301, S_0310;
    real  S_0400, S_1000, S_1001, S_1002, S_1010, S_1020;
    real  S_1100, S_1101, S_1102, S_1110, S_1120, S_1200;
    real  S_1201, S_1210, S_1300, S_2000, S_2001, S_2002;
    real  S_2010, S_2020, S_2100, S_2101, S_2110, S_2200;
    real  S_3000, S_3001, S_3010, S_3100, S_4000, W;
    real  W_01, W_02, W_03, W_04, W_10, W_11;
    real  W_12, W_13, W_20, W_21, W_22, W_30;
    real  W_31, W_40, gmix, groa, groa_grob, groa_grob_p2;
    real groa_p2, grob, grob_p2, roa, roa_p1f3, roa_p2;
    real roa_p2f3, roa_p3, roa_p4, roa_p4f3, roa_p5f3, roa_p8f3;
    real roa_pm1f3, roa_pm4f3, roa_rob, roa_rob_p10f3, roa_rob_p11f3;
    real roa_rob_p13f3, roa_rob_p14f3, roa_rob_p1f3, roa_rob_p2;
    real roa_rob_p20f3, roa_rob_p22f3, roa_rob_p26f3, roa_rob_p2f3;
    real roa_rob_p3, roa_rob_p4, roa_rob_p4f3, roa_rob_p5;
    real roa_rob_p5f3, roa_rob_p6, roa_rob_p7, roa_rob_p7f3, roa_rob_p8;
    real roa_rob_p8f3, roa_rob_pm1, roa_rob_pm1f3, roa_rob_pm2, rob, rob_p1f3;
    real rob_p2, rob_p2f3, rob_p3, rob_p4, rob_p4f3, rob_p5f3;
    real rob_p8f3, rob_pm1f3, rob_pm4f3;

/* constants */
	
    const real a=0.04918;
    const real b=0.132;
    const real c=0.2533;
    const real d=0.3490;
    const real Cf=(3.0/10.0)*pow((3.0*M_PI*M_PI),(2.0/3.0));
    const real Cx=pow(2.0,(11.0/3.0))*Cf;

/* Powers of constants */
    const real d_p2=d*d;
    const real d_p3=d_p2*d;
    const real d_p4=d_p3*d;
    const real d_p5=d_p4*d;

    const real c_p2=c*c;
    const real c_p3=c_p2*c;
    const real c_p4=c_p3*c;


/* Setting up lower order derivatives */
    FunThirdFuncDrv ds_third;

    drv3_clear(&ds_third);
    LYPFunctional.third(&ds_third, factor, dp);

    /* first derivatives */
    ds->df1000  += ds_third.df1000;
    ds->df0100  += ds_third.df0100;
    ds->df0010  += ds_third.df0010;
    ds->df0001  += ds_third.df0001;
    ds->df00001 += ds_third.df00001;
    /*  second derivatives */
    ds->df2000  += ds_third.df2000;
    ds->df0200  += ds_third.df0200;
    ds->df0020  += ds_third.df0020;
    ds->df0002  += ds_third.df0002;
    /* the mixed derivatives */ 
    ds->df1100  += ds_third.df1100;
    ds->df1010  += ds_third.df1010;
    ds->df1001  += ds_third.df1001;
    ds->df0101  += ds_third.df0101;
    ds->df0110  += ds_third.df0110;
    ds->df10001 += ds_third.df10001;
    ds->df01001 += ds_third.df01001;
    /* third order derivative: density only dependent */
    ds->df3000  += ds_third.df3000;
    ds->df0300  += ds_third.df0300;
    ds->df2100  += ds_third.df2100;
    ds->df1200  += ds_third.df1200;
    /* third order derivative: mixed  */    
    ds->df2010  += ds_third.df2010;
    ds->df2001  += ds_third.df2001;
    ds->df0201  += ds_third.df0201;
    ds->df0210  += ds_third.df0210;
    ds->df1020  += ds_third.df1020;
    ds->df0102  += ds_third.df0102;
    ds->df1002  += ds_third.df1002;
    ds->df0120  += ds_third.df0120;
    ds->df1110  += ds_third.df1110;
    ds->df1101  += ds_third.df1101;
    ds->df20001 += ds_third.df20001;
    ds->df02001 += ds_third.df02001;
    ds->df11001 += ds_third.df11001;
    
/* Setting up roa, rob, groa, grob */

    roa=dp->rhoa;
    rob=dp->rhob;
    groa=dp->grada;
    grob=dp->gradb;
    gmix=dp->gradab;

/* Powers of roa, rob */
    roa_p2=roa*roa;
    roa_p3=roa_p2*roa;
    roa_p4=roa_p2*roa;

    roa_p1f3=pow(roa,(1.0/3.0));
    roa_p2f3=roa_p1f3*roa_p1f3;
    roa_p4f3=roa*roa_p1f3;
    roa_p5f3=roa*roa_p2f3;
    roa_p8f3=roa_p2*roa_p2f3;
    roa_pm1f3=1.0/roa_p1f3;
    roa_pm4f3=(1.0/roa)*roa_pm1f3;


    rob_p2=rob*rob;
    rob_p3=rob_p2*rob;
    rob_p4=rob_p2*rob;

    rob_p1f3=pow(rob,(1.0/3.0));
    rob_p2f3=rob_p1f3*rob_p1f3;
    rob_p4f3=rob*rob_p1f3;
    rob_p5f3=rob*rob_p2f3;
    rob_p8f3=rob_p2*rob_p2f3;
    rob_pm1f3=1.0/rob_p1f3;
    rob_pm4f3=(1.0/rob)*rob_pm1f3;

/* Powers of groa, grob, groa+grob */
    groa_p2=groa*groa;
    grob_p2=grob*grob;
    groa_grob=groa+grob;
    groa_grob_p2=groa_grob*groa_grob;


/* Powers of roa_rob  */
    roa_rob=(roa+rob);

    roa_rob_pm1=1.0/roa_rob;
    roa_rob_p2=roa_rob*roa_rob;
    roa_rob_pm2=1.0/roa_rob_p2;
    roa_rob_p3=roa_rob_p2*roa_rob;
    roa_rob_p4=roa_rob_p2*roa_rob_p2;
    roa_rob_p5=roa_rob_p4*roa_rob;
    roa_rob_p6=roa_rob_p4*roa_rob_p2;
    roa_rob_p7=roa_rob_p6*roa_rob;
    roa_rob_p8=roa_rob_p4*roa_rob_p4;

    roa_rob_p1f3=pow(roa_rob,(1.0/3.0));
    roa_rob_p2f3=roa_rob_p1f3*roa_rob_p1f3;
    roa_rob_p4f3=roa_rob*roa_rob_p1f3;
    roa_rob_p5f3=roa_rob_p4f3*roa_rob_p1f3;
    roa_rob_p7f3=roa_rob_p2*roa_rob_p1f3;
    roa_rob_p8f3=roa_rob_p7f3*roa_rob_p1f3;

    roa_rob_p10f3=roa_rob_p3*roa_rob_p1f3;
    roa_rob_p11f3=roa_rob_p10f3*roa_rob_p1f3;
    roa_rob_p13f3=roa_rob_p4*roa_rob_p1f3;
    roa_rob_p14f3=roa_rob_p4*roa_rob_p2f3;
    roa_rob_p20f3=roa_rob_p6*roa_rob_p2f3;
    roa_rob_p22f3=roa_rob_p7*roa_rob_p1f3;
    roa_rob_p26f3=roa_rob_p8*roa_rob_p2f3;

    roa_rob_pm1f3=1.0/roa_rob_p1f3;

/* Powers of Common Denominator ComDenom */
    ComDenom=(d + roa_rob_p1f3);
    ComDenom_p2=ComDenom*ComDenom;
    ComDenom_p3=ComDenom_p2*ComDenom;
    ComDenom_p4=ComDenom_p2*ComDenom_p2;
    ComDenom_p5=ComDenom_p4*ComDenom;


    EXPC=exp(-c*roa_rob_pm1f3);



/*
  %F0_4=-(320.0/81.0)*d_p4*roa*rob*roa_rob_pm14f3*ComDenom_pm5 - 
  %	(1120.0/81.0)*d*roa_p2*rob*roa_rob_pm14f3*ComDenom_pm5 - 
  %	(1280.0/27.0)*d_p4*rob_p2*roa_rob_pm14f3*ComDenom_pm5 - 
  %	(26560.0/81.0)*d*roa*rob_p2*roa_rob_pm14f3*ComDenom_pm5 - 
  %	(8480.0/27.0)*d*rob_p3*roa_rob_pm14f3*ComDenom_pm5 - 
  %	(1376.0/81.0)*d_p3*roa*rob*roa_rob_pm13f3*ComDenom_pm5 - 
  %	(1984.0/9.0)*d_p3*rob_p2*roa_rob_pm13f3*ComDenom_pm5 - 
  %	(96.0)*roa*rob_p2*roa_rob_pm13f3*ComDenom_pm5 - 
  %	(96.0)*rob_p3*roa_rob_pm13f3*ComDenom_pm5 - 
  %	(2080.0/81.0)*d_p2*roa*rob*roa_rob_pm4*ComDenom_pm5 - 
  %	(3520.0/9.0)*d_p2*rob_p2*roa_rob_pm4*ComDenom_pm5

  %F0_4=-(32.0*rob*(243.0*rob*roa_rob_p4f3 + 10.0*d_p4*(roa+12.0*rob) + 
  %	5.0*d_p2*roa_rob_p2f3*(13.0*roa + 198.0*rob) + 
  %	d_p3*roa_rob_p1f3*(43.0*roa + 558.0*rob) + 
  %	5.0*d*(7.0*roa_p2 + 166.0*roa*rob + 159.0*rob_p2))) * 
  %	(1.0/81.0)*roa_rob_pm14f3*ComDenom_pm5

*/

    F0_40=((-(((32.0*rob*((243.0*rob*roa_rob_p4f3+10.0*d_p4*((roa+12.0*rob))+5.0*d_p2*roa_rob_p2f3*((13.0*roa+198.0*rob))+d_p3*roa_rob_p1f3*((43.0*roa+558.0*rob))+5.0*d*((7.0*roa_p2+166.0*roa*rob+159.0*rob_p2))))))/((81.0*roa_rob_p14f3*ComDenom_p5)))));


    F0_31=(((8.0*((d_p3*roa_rob_p1f3*((63.0*roa_p2+1070.0*roa*rob-1053.0*rob_p2))+10.0*d_p2*roa_rob_p2f3*((9.0*roa_p2+190.0*roa*rob-189.0*rob_p2))+5.0*d_p4*((3.0*roa_p2+46.0*roa*rob-45.0*rob_p2))+486.0*rob*roa_rob_p1f3*((roa_p2-rob_p2))+2.0*d*((21.0*roa_p3+788.0*roa_p2*rob-7.0*roa*rob_p2-774.0*rob_p3))))))/((81.0*roa_rob_p14f3*ComDenom_p5)));


    F0_21=((-(((8.0*((d_p2*roa_rob_p1f3*((9.0*roa_p2+110.0*roa*rob-45.0*rob_p2))+2.0*d*roa_rob_p2f3*((3.0*roa_p2+65.0*roa*rob-30.0*rob_p2))+d_p3*((3.0*roa_p2+31.0*roa*rob-12.0*rob_p2))+27.0*rob*((2.0*roa_p2+roa*rob-rob_p2))))))/((27.0*roa_rob_p11f3*ComDenom_p4)))));


    F0_11=((4.0*((18.0*roa*rob*roa_rob_p2f3+d_p2*((3.0*roa_p2+16.0*roa*rob+3.0*rob_p2))+d*roa_rob_p1f3*((3.0*roa_p2+32.0*roa*rob+3.0*rob_p2)))))/(9.0*roa_rob_p8f3*ComDenom_p3));


    F0_22=((-(((16.0*((10.0*d_p4*((3.0*roa_p2-16.0*roa*rob+3.0*rob_p2))+10.0*d_p2*roa_rob_p2f3*((27.0*roa_p2-131.0*roa*rob+27.0*rob_p2))+2.0*d_p3*roa_rob_p1f3*((72.0*roa_p2-371.0*roa*rob+72.0*rob_p2))+81.0*roa_rob_p1f3*((roa_p3-3.0*roa_p2*rob-3.0*roa*rob_p2+rob_p3))+d*((237.0*roa_p3-809.0*roa_p2*rob-809.0*roa*rob_p2+237.0*rob_p3))))))/((81.0*roa_rob_p14f3*ComDenom_p5)))));


    F0_12=(((8.0*((d_p2*roa_rob_p1f3*((45.0*roa_p2-110.0*roa*rob-9.0*rob_p2))+2.0*d*roa_rob_p2f3*((30.0*roa_p2-65.0*roa*rob-3.0*rob_p2))+d_p3*((12.0*roa_p2-31.0*roa*rob-3.0*rob_p2))+27.0*roa*((roa_p2-roa*rob-2.0*rob_p2))))))/((27.0*roa_rob_p11f3*ComDenom_p4)));


    F0_13=((-(((8.0*((d_p3*roa_rob_p1f3*((1053.0*roa_p2-1070.0*roa*rob-63.0*rob_p2))+10.0*d_p2*roa_rob_p2f3*((189.0*roa_p2-190.0*roa*rob-9.0*rob_p2))+5.0*d_p4*((45.0*roa_p2-46.0*roa*rob-3.0*rob_p2))+486.0*roa*roa_rob_p1f3*((roa_p2-rob_p2))+2.0*d*((774.0*roa_p3+7.0*roa_p2*rob-788.0*roa*rob_p2-21.0*rob_p3))))))/((81.0*roa_rob_p14f3*ComDenom_p5)))));


    F0_04=((-(((32.0*roa*((243.0*roa*roa_rob_p4f3+10.0*d_p4*((12.0*roa+rob))+5.0*d_p2*roa_rob_p2f3*((198.0*roa+13.0*rob))+d_p3*roa_rob_p1f3*((558.0*roa+43.0*rob))+5.0*d*((159.0*roa_p2+166.0*roa*rob+7.0*rob_p2))))))/((81.0*roa_rob_p14f3*ComDenom_p5)))));


/*
  %W_4=(((EXPC*((c_p4*((d_p4+4.0*d_p3*roa_rob_p1f3+6.0*d_p2*roa_rob_p2f3+4.0*d*roa_rob+roa_rob_p4f3))*-8.0*c*roa_rob*((1605.0*d_p4+6730.0*d_p3*roa_rob_p1f3+10602.0*d_p2*roa_rob_p2f3+7437.0*d*roa_rob+1960.0*roa_rob_p4f3))*-4.0*c_p3*((16.0*d_p4*roa_rob_p1f3+65.0*d_p3*roa_rob_p2f3+99.0*d_p2*roa_rob+67.0*d*roa_rob_p4f3+17.0*roa_rob_p5f3))+8.0*roa_rob*((4940.0*d_p4*roa_rob_p1f3+21055.0*d_p3*roa_rob_p2f3+33793.0*d_p2*roa_rob+24220.0*d*roa_rob_p4f3+6545.0*roa_rob_p5f3))+4.0*c_p2*((355.0*d_p4*roa_rob_p2f3+1465.0*d_p3*roa_rob+2268.0*d_p2*roa_rob_p4f3+1561.0*d*roa_rob_p5f3+403.0*roa_rob_p2))))))/((81.0*roa_rob_p26f3*ComDenom_p5)))
*/

    W_40=(((EXPC*((c_p4*((d_p4+4.0*d_p3*roa_rob_p1f3+6.0*d_p2*roa_rob_p2f3+ 
                          4.0*d*roa_rob+roa_rob_p4f3))-8.0*c*roa_rob*((1605.0*d_p4+6730.0*d_p3*roa_rob_p1f3+ 
                                                                       10602.0*d_p2*roa_rob_p2f3+7437.0*d*roa_rob+1960.0*roa_rob_p4f3))- 
                   4.0*c_p3*((16.0*d_p4*roa_rob_p1f3+65.0*d_p3*roa_rob_p2f3+99.0*d_p2*roa_rob+ 
                              67.0*d*roa_rob_p4f3+17.0*roa_rob_p5f3))+8.0*roa_rob*((4940.0*d_p4*roa_rob_p1f3+ 
                                                                                    21055.0*d_p3*roa_rob_p2f3+33793.0*d_p2*roa_rob+ 
                                                                                    24220.0*d*roa_rob_p4f3+6545.0*roa_rob_p5f3))+ 
                   4.0*c_p2*((355.0*d_p4*roa_rob_p2f3+1465.0*d_p3*roa_rob+ 
                              2268.0*d_p2*roa_rob_p4f3+1561.0*d*roa_rob_p5f3+403.0*roa_rob_p2))))))/((81.0*roa_rob_p26f3*ComDenom_p5)));


    W_30=(((EXPC*((c_p3*((d_p3+roa+rob+3.0*d_p2*roa_rob_p1f3+3.0*d*roa_rob_p2f3))-2.0*roa_rob*((1040.0*d_p3+3350.0*d_p2*roa_rob_p1f3+3616.0*d*roa_rob_p2f3+1309.0*roa_rob))-3.0*c_p2*((14.0*d_p3*roa_rob_p1f3+43.0*d_p2*roa_rob_p2f3+44.0*d*roa_rob+15.0*roa_rob_p4f3))+2.0*c*((269.0*d_p3*roa_rob_p2f3+846.0*d_p2*roa_rob+888.0*d*roa_rob_p4f3+311.0*roa_rob_p5f3))))))/((27.0*roa_rob_p22f3*ComDenom_p4)));



    W_20=(((EXPC*(((-2.0)*c*roa_rob*((12.0*d_p2+25.0*d*roa_rob_p1f3+13.0*roa_rob_p2f3))+c_p2*((d_p2*roa_rob_p2f3+2.0*d*roa_rob+roa_rob_p4f3))+2.0*((77.0*roa_p2+rob*((77.0*rob+65.0*d_p2*roa_rob_p1f3+141.0*d*roa_rob_p2f3))+roa*((154.0*rob+65.0*d_p2*roa_rob_p1f3+141.0*d*roa_rob_p2f3))))))))/((9.0*roa_rob_p20f3*ComDenom_p3)));



    W_10=((EXPC*(((-10.0)*d*roa_rob_p1f3-11.0*roa_rob_p2f3+c*ComDenom)))/(3.0*roa_rob_p14f3*ComDenom_p2));


/* It looks like W mixed derivatives are same as not mixed ones

%W_31=(((EXPC*((c_p4*((d_p4+4.0*d_p3*roa_rob_p1f3+6.0*d_p2*roa_rob_p2f3+4.0*d*roa_rob+roa_rob_p4f3))-8.0*c*roa_rob*((1605.0*d_p4+6730.0*d_p3*roa_rob_p1f3+10602.0*d_p2*roa_rob_p2f3+7437.0*d*roa_rob+1960.0*roa_rob_p4f3))-4.0*c_p3*((16.0*d_p4*roa_rob_p1f3+65.0*d_p3*roa_rob_p2f3+99.0*d_p2*roa_rob+67.0*d*roa_rob_p4f3+17.0*roa_rob_p5f3))+8.0*roa_rob*((4940.0*d_p4*roa_rob_p1f3+21055.0*d_p3*roa_rob_p2f3+33793.0*d_p2*roa_rob+24220.0*d*roa_rob_p4f3+6545.0*roa_rob_p5f3))+4.0*c_p2*((355.0*d_p4*roa_rob_p2f3+1465.0*d_p3*roa_rob+2268.0*d_p2*roa_rob_p4f3+1561.0*d*roa_rob_p5f3+403.0*roa_rob_p2))))))/((81.0*roa_rob_p26f3*ComDenom_p5)))


%W_21=(((EXPC*((c_p3*((d_p3+roa+rob+3.0*d_p2*roa_rob_p1f3+3.0*d*roa_rob_p2f3))-2.0*roa_rob*((1040.0*d_p3+3350.0*d_p2*roa_rob_p1f3+3616.0*d*roa_rob_p2f3+1309.0*roa_rob))-3.0*c_p2*((14.0*d_p3*roa_rob_p1f3+43.0*d_p2*roa_rob_p2f3+44.0*d*roa_rob+15.0*roa_rob_p4f3))+2.0*c*((269.0*d_p3*roa_rob_p2f3+846.0*d_p2*roa_rob+888.0*d*roa_rob_p4f3+311.0*roa_rob_p5f3))))))/((27.0*roa_rob_p22f3*ComDenom_p4)))

*/

    W_31=W_40;
    W_21=W_30;
    W_11=W_20;

    W_22=W_40;
    W_04=W_40;

    W_13=W_40;
    W_12=W_30;

    W_03=W_30;
    W_02=W_20;
    W_01=W_10;

    W=EXPC/(roa_rob_p10f3*ComDenom);


    Delta_40=(((8.0*((35.0*c*((d_p5+5.0*d_p4*roa_rob_p1f3+10.0*d_p3*roa_rob_p2f3+10.0*d_p2*roa_rob+5.0*d*roa_rob_p4f3+roa_rob_p5f3))+d*((10.0*d_p3*roa_rob_p2f3+43.0*d_p2*roa_rob+65.0*d*roa_rob_p4f3+35.0*roa_rob_p5f3))))))/((81.0*roa_rob_p13f3*ComDenom_p5)));


    Delta_30=((-(((2.0*((14.0*c*((d_p4+4.0*d_p3*roa_rob_p1f3+6.0*d_p2*roa_rob_p2f3+4.0*d*roa_rob+roa_rob_p4f3))+d*((5.0*d_p2*roa_rob_p2f3+16.0*d*roa_rob+14.0*roa_rob_p4f3))))))/((27.0*roa_rob_p10f3*ComDenom_p4)))));


    Delta_20=((2.0*((d*((2.0*roa+2.0*rob+d*roa_rob_p2f3))+2.0*c*((d_p3+roa+rob+3.0*d_p2*roa_rob_p1f3+3.0*d*roa_rob_p2f3)))))/(9.0*roa_rob_p7f3*ComDenom_p3));


    Delta_10=(((-d)*roa_rob_p2f3-c*ComDenom_p2)/(3.0*roa_rob_p4f3*ComDenom_p2));

/* Same hold for Delta, mixed derivatives are equal to total derivatie order */

    Delta_31=Delta_40;
    Delta_21=Delta_30;
    Delta_11=Delta_20;

    Delta_22=Delta_40;
    Delta_12=Delta_30;
    Delta_13=Delta_40;
    Delta_04=Delta_40;
    Delta_03=Delta_30;
    Delta_02=Delta_20;
    Delta_01=Delta_10;

    Delta=c*roa_rob_pm1f3 + d*roa_rob_pm1f3*(1.0/(1.0+d*roa_rob_pm1f3));

/* S auxaliary function and derivatives */

    S=(roa*groa_p2*roa_rob_pm1) + (rob*grob_p2*roa_rob_pm1);

    S_4000=24.0*(grob_p2-groa_p2)*rob/roa_rob_p5;

    S_3000=6.0*(groa_p2-grob_p2)*rob/roa_rob_p4;

    S_2000=2.0*(grob_p2-groa_p2)*rob/roa_rob_p3;

    S_1000=(groa_p2-grob_p2)*rob/roa_rob_p2;

    S_0400=24.0*(groa_p2-grob_p2)*roa/roa_rob_p5;

    S_0300=6.0*(grob_p2-groa_p2)*roa/roa_rob_p4;

    S_0200=2.0*(groa_p2-grob_p2)*roa/roa_rob_p3;

    S_0100=(grob_p2-groa_p2)*roa/roa_rob_p2;

    S_3100=6.0*(groa_p2-grob_p2)*(roa-3.0*rob)/roa_rob_p5;

    S_2100=-2.0*(groa_p2-grob_p2)*(roa-2.0*rob)/roa_rob_p4;

    S_1100=(groa_p2-grob_p2)*(roa-rob)/roa_rob_p3;

    S_1200=-2.0*(groa_p2-grob_p2)*(2.0*roa-rob)/roa_rob_p4; 

    S_2200=12.0*(groa_p2-grob_p2)*(roa-rob)/roa_rob_p5;

    S_1300=6.0*(groa_p2-grob_p2)*(3.0*roa-rob)/roa_rob_p5;

    S_0010=2.0*groa*roa/roa_rob;
    S_0001=2.0*grob*rob/roa_rob;

    S_1010=2.0*groa*rob/roa_rob_p2;
    S_0110=-2.0*groa*roa/roa_rob_p2;
    S_0101=2.0*grob*roa/roa_rob_p2;
    S_1001=-2.0*grob*rob/roa_rob_p2;

    S_2010=-4.0*groa*rob/roa_rob_p3;
    S_2001=4.0*grob*rob/roa_rob_p3;

    S_0210=4.0*groa*roa/roa_rob_p3;


    S_3010=12.0*groa*rob/roa_rob_p4;
    S_3001=-12.0*grob*rob/roa_rob_p4;

    S_1110=2.0*groa*(roa-rob)/roa_rob_p3;
    S_1101=2.0*grob*(rob-roa)/roa_rob_p3;

    S_2110=-4.0*groa*(roa-2.0*rob)/roa_rob_p4;
    S_2101=4.0*grob*(roa-2.0*rob)/roa_rob_p4;

    S_0020=2.0*roa/roa_rob;
    S_0002=2.0*rob/roa_rob;

    S_1020=2.0*rob/roa_rob_p2;
    S_1002=-2.0*rob/roa_rob_p2;

    S_2020=-4.0*rob/roa_rob_p3;
    S_2002=4.0*rob/roa_rob_p3;

    S_1210=4.0*groa*(rob-2.0*roa)/roa_rob_p4;

    S_0201=-4.0*grob*roa/roa_rob_p3;

    S_1201=grob*(8.0*roa-4.0*rob)/roa_rob_p4;

    S_1120=2.0*(roa-rob)/roa_rob_p3;

    S_0120=-2.0*roa/roa_rob_p2;

    S_0102=-S_0120;

    S_1102=-S_1120;

    S_0310=-12.0*groa*roa/roa_rob_p4;

    S_0301=12.0*grob*roa/roa_rob_p4;

    S_0220=4.0*roa/roa_rob_p3;

    S_0202=-S_0220;

/* Often used term in F1 derivatives */ 
    F1term=(7.0*gmix + 3.0*(groa_p2 + grob_p2) + S);

/* F1 auxaliary function and its derivatives */
    F1=Cx*(roa_p8f3 + rob_p8f3) + (1.0/18.0)*(groa_p2 + 2.0*gmix + grob_p2)*(47.0-7.0*Delta) - 
	(1.0/18.0)*(groa_p2+grob_p2)*(45.0-Delta) - 
	(1.0/9.0)*(Delta-11.0)*S;

    F1_40000=(1.0/81.0)*(-80.0*Cx*roa_pm4f3 - 9.0*(7.0*gmix + 3.0*(groa_p2 + grob_p2) + S)*Delta_40 - 
                         36.0*Delta_30*S_1000 - 54.0*Delta_20*S_2000 -36.0*Delta_10*S_3000 + 99.0*S_4000 - 9.0*Delta*S_4000);

    F1_30000=(1.0/27.0)*(80.0*Cx*roa_pm1f3 -3.0*Delta_30*(7.0*gmix + 3.0*(groa_p2 + grob_p2) + S) - 
                         9.0*Delta_20*S_1000 -9.0*Delta_10*S_2000 + 33.0*S_3000 -3.0*Delta*S_3000);

    F1_20000=(1.0/9.0)*(40.0*Cx*roa_p2f3 - Delta_20*(7.0*gmix + 3.0*(groa_p2 + grob_p2) + S) - 
                        2.0*Delta_10*S_1000 + 11.0*S_2000 - Delta*S_2000);

    F1_10000=(1.0/9.0)*(24.0*Cx*roa_p5f3 - Delta_10*(7.0*gmix + 3.0*(groa_p2 + grob_p2) + S) - (-11.0 + Delta)*S_1000);


    F1_31000=(1.0/9.0)*(-Delta_31*F1term - Delta_30*S_0100 - 3.0*Delta_21*S_1000 - 
                        3.0*Delta_20*S_1100 -3.0*Delta_11*S_2000 - 3.0*Delta_10*S_2100 - Delta_01*S_3000 + 
                        11.0*S_3100 - Delta*S_3100);

    F1_21000=(1.0/9.0)*(-Delta_21*F1term - Delta_20*S_0100 - 2.0*Delta_11*S_1000 - 2.0*Delta_10*S_1100 - 
                        Delta_01*S_2000 + 11.0*S_2100 - Delta*S_2100);

    F1_11000=(1.0/9.0)*(-Delta_11*F1term - Delta_10*S_0100 - Delta_01*S_1000 + 11.0*S_1100 - Delta*S_1100);

    F1_01000=(1.0/9.0)*(24.0*Cx*rob_p5f3 - Delta_01*F1term - (-11.0 + Delta)*S_0100);


    F1_02000=(1.0/9.0)*(40.0*Cx*rob_p2f3 - Delta_02*F1term - 2.0*Delta_01*S_0100 + 11.0*S_0200 - Delta*S_0200);

    F1_12000=(1.0/9.0)*(-Delta_12*F1term -2.0*Delta_11*S_0100 - Delta_10*S_0200 - Delta_02*S_1000 - 
                        2.0*Delta_01*S_1100 + 11.0*S_1200 - Delta*S_1200);

    F1_22000=(1.0/9.0)*(-Delta_22*F1term - 2.0*Delta_21*S_0100 - 
                        Delta_20*S_0200 - 2.0*Delta_12*S_1000 - 4.0*Delta_11*S_1100 - 
                        2.0*Delta_10*S_1200 - Delta_02*S_2000 - 2.0*Delta_01*S_2100 + 
                        11.0*S_2200 - Delta*S_2200);

    F1_03000=(1.0/27.0)*(80.0*Cx*rob_pm1f3 -3.0*Delta_03*(7.0*gmix + 3.0*(groa_p2 + grob_p2) + S) - 
                         9.0*Delta_02*S_0100 -9.0*Delta_01*S_0200 + 33.0*S_0300 -3.0*Delta*S_0300);

    F1_13000=(1.0/9.0)*(-Delta_13*F1term - Delta_03*S_1000 - 3.0*Delta_12*S_0100 - 
                        3.0*Delta_02*S_1100 -3.0*Delta_11*S_0200 - 3.0*Delta_01*S_1200 - Delta_10*S_0300 + 
                        11.0*S_1300 - Delta*S_1300);

    F1_04000=(1.0/81.0)*(-80.0*Cx*rob_pm4f3 - 9.0*F1term*Delta_04 - 
                         36.0*Delta_03*S_0100 - 54.0*Delta_02*S_0200 -36.0*Delta_01*S_0300 + 99.0*S_0400 - 9.0*Delta*S_0400);


    F1_00100=(1.0/9.0)*(2.0*groa + 11.0*S_0010 - Delta*(6.0*groa+S_0010));
    F1_00010=(1.0/9.0)*(2.0*grob + 11.0*S_0001 - Delta*(6.0*grob+S_0001));

    F1_10100=(1.0/9.0)*(-Delta_10*(6.0*groa + S_0010) - (-11.0 + Delta)*S_1010 );
    F1_01100=(1.0/9.0)*(-Delta_01*(6.0*groa + S_0010) - (-11.0 + Delta)*S_0110 );
    F1_01010=(1.0/9.0)*(-Delta_01*(6.0*grob + S_0001) - (-11.0 + Delta)*S_0101 ); 
    F1_10010=(1.0/9.0)*(-Delta_10*(6.0*grob + S_0001) - (-11.0 + Delta)*S_1001 );

    F1_01001=-(7.0/9.0)*Delta_01;

    F1_20100=(1.0/9.0)*(-Delta_20*(6.0*groa + S_0010) - 2.0*Delta_10*S_1010 - (-11.0 + Delta)*S_2010 );
    F1_20010=(1.0/9.0)*(-Delta_20*(6.0*grob + S_0001) - 2.0*Delta_10*S_1001 - (-11.0 + Delta)*S_2001 );

    F1_02100=(1.0/9.0)*(-Delta_02*(6.0*groa + S_0010) - 2.0*Delta_01*S_0110 - (-11.0 + Delta)*S_0210 );
    F1_02010=(1.0/9.0)*(-Delta_02*(6.0*grob + S_0001) - 2.0*Delta_01*S_0101 - (-11.0 + Delta)*S_0201 );

    F1_30100=(1.0/9.0)*(-Delta_30*(6.0*groa + S_0010) - 3.0*Delta_20*S_1010 - 
                        3.0*Delta_10*S_2010 - (-11.0 + Delta)*S_3010);
    F1_30010=(1.0/9.0)*(-Delta_30*(6.0*grob + S_0001) - 3.0*Delta_20*S_1001 - 
                        3.0*Delta_10*S_2001 - (-11.0 + Delta)*S_3001);

    F1_03100=(1.0/9.0)*(-Delta_03*(6.0*groa + S_0010) - 3.0*Delta_02*S_0110 - 
                        3.0*Delta_01*S_0210 - (-11.0 + Delta)*S_0310);

    F1_03010=(1.0/9.0)*(-Delta_03*(6.0*grob + S_0001) - 3.0*Delta_02*S_0101 - 
                        3.0*Delta_01*S_0201 - (-11.0 + Delta)*S_0301);

    F1_30001=-(7.0/9.0)*Delta_30;

    F1_20001=-(7.0/9.0)*Delta_20;

    F1_10001=-(7.0/9.0)*Delta_10;

    F1_00001=(1.0/9.0)*(47.0-7.0*Delta);

    F1_11100=(1.0/9.0)*(-Delta_11*(6.0*groa + S_0010) - Delta_10*S_0110 - 
                        Delta_01*S_1010 - (-11.0 + Delta)*S_1110);

    F1_11010=(1.0/9.0)*(-Delta_11*(6.0*grob + S_0001) - Delta_10*S_0101 - 
                        Delta_01*S_1001 - (-11.0 + Delta)*S_1101);

    F1_21100=(1.0/9.0)*(-Delta_21*(6.0*groa + S_0010) - 
                        Delta_20*S_0110 - 2.0*Delta_11*S_1010  - 
                        2.0*Delta_10*S_1110 - Delta_01*S_2010 + 
                        + 11.0*S_2110 - Delta*S_2110);

    F1_21010=(1.0/9.0)*(-Delta_21*(6.0*grob + S_0001) - 
                        Delta_20*S_0101 - 2.0*Delta_11*S_1001  - 
                        2.0*Delta_10*S_1101 - Delta_01*S_2001 + 
                        + 11.0*S_2101 - Delta*S_2101);

    F1_11001=-(7.0/9.0)*Delta_11;

    F1_21001=-(7.0/9.0)*Delta_21;

    F1_00200=(1.0/9.0)*(2.0 + 11.0*S_0020 -Delta*(6.0 + S_0020));
    F1_00020=(1.0/9.0)*(2.0 + 11.0*S_0002 -Delta*(6.0 + S_0002));

    F1_10200=(1.0/9.0)*(-Delta_10*(6.0+S_0020) - S_1020*(-11.0 + Delta));
    F1_10020=(1.0/9.0)*(-Delta_10*(6.0+S_0002) - S_1002*(-11.0 + Delta));

    F1_20200=(1.0/9.0)*(-Delta_20*(6.0+S_0020) - 2.0*Delta_10*S_1020 - S_2020*(-11.0 + Delta));
    F1_20020=(1.0/9.0)*(-Delta_20*(6.0+S_0002) - 2.0*Delta_10*S_1002 - S_2002*(-11.0 + Delta));

    F1_02200=(1.0/9.0)*(-Delta_02*(6.0+S_0020) - 2.0*Delta_01*S_0120 - S_0220*(-11.0 + Delta));
    F1_02020=(1.0/9.0)*(-Delta_02*(6.0+S_0002) - 2.0*Delta_01*S_0102 - S_0202*(-11.0 + Delta));

    F1_12100=(1.0/9.0)*(-Delta_12*(6.0*groa + S_0010) - 2.0*Delta_11*S_0110 - Delta_10*S_0210 - 
                        Delta_02*S_1010 -2.0*Delta_01*S_1110 + 11.0*S_1210 - Delta*S_1210);

    F1_12010=(1.0/9.0)*(-Delta_12*(6.0*grob + S_0001) - 2.0*Delta_11*S_0101 - Delta_10*S_0201 - 
                        Delta_02*S_1001 -2.0*Delta_01*S_1101 + 11.0*S_1201 - Delta*S_1201);


    F1_02001=-(7.0/9.0)*Delta_02;

    F1_12001=-(7.0/9.0)*Delta_12;

    F1_01200=(1.0/9.0)*(-Delta_01*(6.0+S_0020) - (-11.0 + Delta)*S_0120);

    F1_01020=(1.0/9.0)*(-Delta_01*(6.0+S_0002) - (-11.0 + Delta)*S_0102);

    F1_11200=(1.0/9.0)*(-Delta_11*(6.0+ S_0020) - Delta_10*S_0120 - Delta_01*S_1020 - 
                        (-11.0 + Delta)*S_1120);

    F1_11020=(1.0/9.0)*(-Delta_11*(6.0+ S_0002) - Delta_10*S_0102 - Delta_01*S_1002 - 
                        (-11.0 + Delta)*S_1102);

    F1_03001=-(7.0/9.0)*Delta_03;

/* F2 Function and derivatives */

    F2=-grob_p2*roa_p2 - groa_p2*rob_p2 - (4.0/3.0)*roa_rob_p2*gmix;

    F2_20000=-(8.0/3.0)*gmix - 2.0*grob_p2;

    F2_10000=-2.0*grob_p2*roa - (8.0/3.0)*gmix*roa_rob;

    F2_01000=-2.0*groa_p2*rob - (8.0/3.0)*gmix*roa_rob;

    F2_11000=-(8.0/3.0)*gmix;

    F2_02000=-(8.0/3.0)*gmix - 2.0*groa_p2;

    F2_00100=-2.0*groa*rob_p2;

    F2_01100=-4.0*groa*rob;

    F2_00001=-(4.0/3.0)*roa_rob_p2;

    F2_01200=-4.0*rob;

    F2_10020=-4.0*roa;

    F2_20020=-4.0;

    F2_00010=-2.0*grob*roa_p2;

    F2_00001=-(4.0/3.0)*roa_rob_p2;

    F2_10001=-(8.0/3.0)*roa_rob;

    F2_20001=-(8.0/3.0);

    F2_10010=-4.0*grob*roa;

    F2_00020=-2.0*roa_p2;

    F2_00200=-2.0*rob_p2;

    F2_02100=-4.0*groa;

    F2_01001=-(8.0/3.0)*roa_rob;

    F2_02001=-(8.0/3.0);

    F2_11001=-(8.0/3.0);

    F2_10010=-4.0*grob*roa;

    F2_20010=-4.0*grob;

    F2_11001=-(8.0/3.0);

    F2_01200=-4.0*rob;
    F2_02200=-4.0;

/*Derivatives */

    ds->df4000 += factor*
        ( -a*F0_40 - a*b*W_40*(roa*rob*F1 + F2) - 
          4.0*a*b*W_30*(rob*F1 + roa*rob*F1_10000 + F2_10000) - 
          6.0*a*b*W_20*(2.0*rob*F1_10000 + roa*rob*F1_20000 + F2_20000) - 
          4.0*a*b*W_10*(3.0*rob*F1_20000 + roa*rob*F1_30000) - 
          a*b*W*(4.0*rob*F1_30000 + roa*rob*F1_40000) 
        );


    ds->df3100 += factor*
        ( -a*F0_31 - 
          a*b*W_31*(roa*rob*F1 + F2) - 
          a*b*W_30*(roa*F1 + roa*rob*F1_01000 + F2_01000) - 
          3.0*a*b*W_21*(rob*F1 + roa*rob*F1_10000 + F2_10000) - 
          3.0*a*b*W_20*(F1 + rob*F1_01000 + roa*F1_10000 + roa*rob*F1_11000 + F2_11000) - 
          3.0*a*b*W_11*(2.0*rob*F1_10000 + roa*rob*F1_20000 + F2_20000) - 
          3.0*a*b*W_10*(2.0*F1_10000 + 2.0*rob*F1_11000 + roa*F1_20000 + roa*rob*F1_21000) -
          a*b*W_01*(3.0*rob*F1_20000 + roa*rob*F1_30000  ) - 
          a*b*W*(3.0*F1_20000 + 3.0*rob*F1_21000 + roa*F1_30000 + roa*rob*F1_31000 ) 
        );


    ds->df2200 += factor*
        ( -a*F0_22 - 
          a*b*W_22*(roa*rob*F1 + F2) - 
          2.0*a*b*W_21*(roa*F1 + roa*rob*F1_01000 + F2_01000) - 
          a*b*W_20*(2.0*roa*F1_01000 + roa*rob*F1_02000 + F2_02000) - 
          2.0*a*b*W_12*(rob*F1 + roa*rob*F1_10000 + F2_10000) - 
          4.0*a*b*W_11*(F1 + rob*F1_01000 + roa*F1_10000 + roa*rob*F1_11000 + F2_11000) - 
          2.0*a*b*W_10*(2.0*F1_01000 + rob*F1_02000 + 2.0*roa*F1_11000 + roa*rob*F1_12000 ) - 
          a*b*W_02*(2.0*rob*F1_10000 + roa*rob*F1_20000 + F2_20000) - 
          2.0*a*b*W_01*(2.0*F1_10000 + 2.0*rob*F1_11000 + roa*F1_20000 + roa*rob*F1_21000 ) - 
          a*b*W*(4.0*F1_11000 + 2.0*rob*F1_12000 + 2.0*roa*F1_21000 + roa*rob*F1_22000)
        );


    ds->df1300 += factor*
        ( -a*F0_13 - 
          a*b*W_13*(roa*rob*F1 + F2) - 
          3.0*a*b*W_12*(roa*F1 + roa*rob*F1_01000 + F2_01000) - 
          3.0*a*b*W_11*(2.0*roa*F1_01000 + roa*rob*F1_02000 + F2_02000) - 
          a*b*W_10*(3.0*roa*F1_02000 + roa*rob*F1_03000 ) - 
          a*b*W_03*(rob*F1 + roa*rob*F1_10000 + F2_10000) - 
          3.0*a*b*W_02*(F1 + rob*F1_01000 + roa*F1_10000 + roa*rob*F1_11000 + F2_11000) - 
          3.0*a*b*W_01*(2.0*F1_01000 + rob*F1_02000 + 2.0*roa*F1_11000 + roa*rob*F1_12000 ) - 
          a*b*W*(3.0*F1_02000 + rob*F1_03000 + 3.0*roa*F1_12000 + roa*rob*F1_13000 ) 
        );
	


    ds->df0400 += factor*
        ( -a*F0_04 - 
          a*b*W_04*(roa*rob*F1 + F2) - 
          4.0*a*b*W_03*(roa*F1 + roa*rob*F1_01000 + F2_01000) - 
          6.0*a*b*W_02*(2.0*roa*F1_01000 + roa*rob*F1_02000 + F2_02000) - 
          4.0*a*b*W_01*(3.0*roa*F1_02000 + roa*rob*F1_03000 ) - 
          a*b*W*(4.0*roa*F1_03000 + roa*rob*F1_04000) 
        );


    ds->df3010 += factor*
        ( -a*b*W_30*(roa*rob*F1_00100 + F2_00100) - 
          3.0*a*b*W_20*(rob*F1_00100 + roa*rob*F1_10100 ) - 
          3.0*a*b*W_10*(2.0*rob*F1_10100 + roa*rob*F1_20100) - 
          a*b*W*(3.0*rob*F1_20100 + roa*rob*F1_30100) 
        );


    ds->df3001 += factor*
        ( -a*b*W_30*(roa*rob*F1_00010 + F2_00010) - 
          3.0*a*b*W_20*(rob*F1_00010 + roa*rob*F1_10010 + F2_10010) - 
          3.0*a*b*W_10*(2.0*rob*F1_10010 + roa*rob*F1_20010 + F2_20010) - 
          a*b*W*(3.0*rob*F1_20010 + roa*rob*F1_30010) 
        );


    ds->df30001 += factor*
        ( -a*b*W_30*(roa*rob*F1_00001 + F2_00001) - 
          3.0*a*b*W_20*(rob*F1_00001 + roa*rob*F1_10001 + F2_10001) - 
          3.0*a*b*W_10*(2.0*rob*F1_10001 + roa*rob*F1_20001 + F2_20001) - 
          a*b*W*(3.0*rob*F1_20001 + roa*rob*F1_30001) 
        );


    ds->df2110 += factor*
        ( -a*b*W_21*(roa*rob*F1_00100 + F2_00100) - 
          a*b*W_20*(roa*F1_00100 + roa*rob*F1_01100 + F2_01100) - 
          2.0*a*b*W_11*(rob*F1_00100 + roa*rob*F1_10100 ) - 
          2.0*a*b*W_10*(F1_00100 + rob*F1_01100 + roa*F1_10100 + roa*rob*F1_11100) - 
          a*b*W_01*(2.0*rob*F1_10100 + roa*rob*F1_20100) - 
          a*b*W*(2.0*F1_10100 + 2.0*rob*F1_11100 + roa*F1_20100 + roa*rob*F1_21100) 
        );


    ds->df2101 += factor*
        ( -a*b*W_21*(roa*rob*F1_00010 + F2_00010) - 
          a*b*W_20*(roa*F1_00010 + roa*rob*F1_01010 ) - 
          2.0*a*b*W_11*(rob*F1_00010 + roa*rob*F1_10010 + F2_10010) - 
          2.0*a*b*W_10*(F1_00010 + rob*F1_01010 + roa*F1_10010 + roa*rob*F1_11010 ) - 
          a*b*W_01*(2.0*rob*F1_10010 + roa*rob*F1_20010 + F2_20010) - 
          a*b*W*(2.0*F1_10010 + 2.0*rob*F1_11010 + roa*F1_20010 + roa*rob*F1_21010) 
        );

    ds->df21001 += factor*
        ( -a*b*W_21*(roa*rob*F1_00001 + F2_00001) - 
          a*b*W_20*(roa*F1_00001 + roa*rob*F1_01001 + F2_01001) - 
          2.0*a*b*W_11*(rob*F1_00001 + roa*rob*F1_10001 + F2_10001) - 
           2.0*a*b*W_10*(F1_00001 + rob*F1_01001 + roa*F1_10001 + roa*rob*F1_11001 + F2_11001) - 
          a*b*W_01*(2.0*rob*F1_10001 + roa*rob*F1_20001 + F2_20001) - 
          a*b*W*(2.0*F1_10001 + 2.0*rob*F1_11001 + roa*F1_20001 + roa*rob*F1_21001 ) 
        );

    ds->df2020 += factor*
        ( -a*b*W_20*(roa*rob*F1_00200 + F2_00200) - 
          2.0*a*b*W_10*(rob*F1_00200 + roa*rob*F1_10200 ) - 
          a*b*W*(2.0*rob*F1_10200 + roa*rob*F1_20200 ) 
        );

    ds->df2002 += factor*
        (  -a*b*W_20*(roa*rob*F1_00020 + F2_00020) - 
           2.0*a*b*W_10*(rob*F1_00020 + roa*rob*F1_10020 + F2_10020 ) - 
           a*b*W*(2.0*rob*F1_10020 + roa*rob*F1_20020 + F2_20020) 
        );

    ds->df1210 += factor*
        ( -a*b*W_12*(roa*rob*F1_00100 + F2_00100) - 
          2.0*a*b*W_11*(roa*F1_00100 + roa*rob*F1_01100 + F2_01100) - 
          a*b*W_10*(2.0*roa*F1_01100 + roa*rob*F1_02100 + F2_02100) - 
          a*b*W_02*(rob*F1_00100 + roa*rob*F1_10100 ) - 
          2.0*a*b*W_01*(F1_00100 + rob*F1_01100 + roa*F1_10100 + roa*rob*F1_11100 ) - 
          a*b*W*(2.0*F1_01100 + rob*F1_02100 + 2.0*roa*F1_11100 + roa*rob*F1_12100 ) 
        );

    ds->df1201 += factor*
        ( -a*b*W_12*(roa*rob*F1_00010 + F2_00010) - 
          2.0*a*b*W_11*(roa*F1_00010 + roa*rob*F1_01010 ) - 
          a*b*W_10*(2.0*roa*F1_01010 + roa*rob*F1_02010 ) - 
          a*b*W_02*(rob*F1_00010 + roa*rob*F1_10010 + F2_10010) - 
          2.0*a*b*W_01*(F1_00010 + rob*F1_01010 + roa*F1_10010 + roa*rob*F1_11010 ) - 
          a*b*W*(2.0*F1_01010 + rob*F1_02010 + 2.0*roa*F1_11010 + roa*rob*F1_12010 ) 
        );


    ds->df12001 += factor*
        ( -a*b*W_12*(roa*rob*F1_00001 + F2_00001) - 
          2.0*a*b*W_11*(roa*F1_00001 + roa*rob*F1_01001 + F2_01001) - 
          a*b*W_10*(2.0*roa*F1_01001 + roa*rob*F1_02001 + F2_02001) - 
          a*b*W_02*(rob*F1_00001 + roa*rob*F1_10001 + F2_01001) - 
          2.0*a*b*W_01*(F1_00001 + rob*F1_01001 + roa*F1_10001 + roa*rob*F1_11001 + F2_11001) - 
          a*b*W*(2.0*F1_01001 + rob*F1_02001 + 2.0*roa*F1_11001 + roa*rob*F1_12001 ) 
        );


    ds->df1120 += factor*
        (  -a*b*W_11*(roa*rob*F1_00200 + F2_00200) - 
           a*b*W_10*(roa*F1_00200 + roa*rob*F1_01200 + F2_01200 ) - 
           a*b*W_01*(rob*F1_00200 + roa*rob*F1_10200 ) - 
           a*b*W*(F1_00200 + rob*F1_01200 + roa*F1_10200 + roa*rob*F1_11200) 
        );


    ds->df1102 += factor*
        (  -a*b*W_11*(roa*rob*F1_00020 + F2_00020) - 
           a*b*W_10*(roa*F1_00020 + roa*rob*F1_01020 ) - 
           a*b*W_01*(rob*F1_00020 + roa*rob*F1_10020 + F2_10020 ) - 
           a*b*W*(F1_00020 + rob*F1_01020 + roa*F1_10020 + roa*rob*F1_11020) 
        );


    ds->df0310 += factor*
        ( -a*b*W_03*(roa*rob*F1_00100 + F2_00100) - 
          3.0*a*b*W_02*(roa*F1_00100 + roa*rob*F1_01100 + F2_01100) - 
          3.0*a*b*W_01*(2.0*roa*F1_01100 + roa*rob*F1_02100 + F2_02100 ) - 
          a*b*W*(3.0*roa*F1_02100 + roa*rob*F1_03100  ) );


    ds->df0301 += factor*
        (   -a*b*W_03*(roa*rob*F1_00010 + F2_00010) - 
            3.0*a*b*W_02*(roa*F1_00010 + roa*rob*F1_01010 ) - 
            3.0*a*b*W_01*(2.0*roa*F1_01010 + roa*rob*F1_02010) - 
            a*b*W*(3.0*roa*F1_02010 + roa*rob*F1_03010) );


    ds->df03001 += factor*
        ( -a*b*W_03*(roa*rob*F1_00001 + F2_00001) - 
          3.0*a*b*W_02*(roa*F1_00001 + roa*rob*F1_01001 + F2_01001) - 
          3.0*a*b*W_01*(2.0*roa*F1_01001 + roa*rob*F1_02001 + F2_02001) - 
          a*b*W*(3.0*roa*F1_02001 + roa*rob*F1_03001) );


    ds->df0220 += factor*
        ( -a*b*W_02*(roa*rob*F1_00200 + F2_00200) - 
          2.0*a*b*W_01*(roa*F1_00200 + roa*rob*F1_01200 + F2_01200 ) - 
          a*b*W*(2.0*roa*F1_01200 + roa*rob*F1_02200 + F2_02200 ) );

    ds->df0202 += factor*
        ( -a*b*W_02*(roa*rob*F1_00020 + F2_00020) - 
          2.0*a*b*W_01*(roa*F1_00020 + roa*rob*F1_01020 ) - 
          a*b*W*(2.0*roa*F1_01020 + roa*rob*F1_02020  ) );

}
