/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* fun-lyp.c:
   implementation of LYP functional and its derivatives 
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
static int lyp_isgga(void) { return 1; }
static int lyp_read(const char* conf_line);
static real lyp_energy(const FunDensProp* dp);
static void lyp_first(FunFirstFuncDrv *ds,   real factor, const FunDensProp* dp);
static void lyp_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp);
static void lyp_third(FunThirdFuncDrv *ds,   real factor, const FunDensProp* dp);

Functional LYPFunctional = {
    "LYP",      /* name */
    lyp_isgga,  /* gga-corrected */
    lyp_read,   /* no extra input expected, just set the common block */
    NULL, /* reporter */
    lyp_energy, 
    lyp_first,
    lyp_second,
    lyp_third
};

/* IMPLEMENTATION PART */
static int
lyp_read(const char* conf_line)
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
lyp_energy(const FunDensProp* dp)
{
    const real A = 0.04918, B = 0.132, C = 0.2533, D = 0.349;
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


/* lyp_first:
   the derivatives of LYP functional wrt alpha density and gradient, suitable
   for unrestricted calculations.
   lyprho : dF/drhoa
   lypgrad: dF/dngrada
   See Doc/dft/functionals.tex for details.
*/
static void
lyp_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
    const real A = 0.04918, B = 0.132, C = 0.2533, D = 0.349;
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

/* lyp_second:
   second order derivatives of LYP functional with respect to 
   alpha density and gradient.
*/
static void
lyp_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    const real A = 0.04918, B = 0.132, C = 0.2533, D = 0.349;
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
lyp_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    const real A = 0.04918, B = 0.132, C = 0.2533, D = 0.349;
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
