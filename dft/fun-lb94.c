/* fun-lb94.c:
   implementation of Exchange-correlation potential with correct 
   asymptotic behavior by R. van Leeuwen and E. J. Baerends;
   (The LB94 functional) and its derivatives.
   (c) B. Jansik, brano@theochem.kth.se, dec 2002
*/

#include <math.h>
#include <stddef.h>

#define __CVERSION__
/* define FOURTH_ORDER_DERIVATIVES */

#include "functionals.h"
#include "dftcom.h"

/* INTERFACE PART */
static int  lb94pot_isgga(void) { return 1; }
static int  lb94pot_read(const char* conf_line);
static real lb94pot_energy(const DftDensProp* dens_prop);
static void lb94pot_first(FirstFuncDrv *ds, real factor, 
                       const DftDensProp* dens_prop);
static void lb94pot_second(SecondFuncDrv *ds, real factor,
                        const DftDensProp* dens_prop);

static void lb94pot_third(ThirdFuncDrv *ds, real factor,
                       const DftDensProp* dens_prop);

#ifdef FOURTH_ORDER_DERIVATIVES
static void lb94pot_fourth(FourthFuncDrv *ds, real factor,
                        const DftDensProp* dens_prop);
#endif

static int  lb94high_isgga(void) { return 1; }
static int  lb94high_read(const char* conf_line);
static real lb94high_energy(const DftDensProp* dens_prop);
static void lb94high_first(FirstFuncDrv *ds, real factor,
                       const DftDensProp* dens_prop);
static void lb94high_second(SecondFuncDrv *ds, real factor,
                        const DftDensProp* dens_prop);
                                                                                                 
static void lb94high_third(ThirdFuncDrv *ds, real factor,
                       const DftDensProp* dens_prop);

#ifdef FOURTH_ORDER_DERIVATIVES
static void lb94high_fourth(FourthFuncDrv *ds, real factor,
                        const DftDensProp* dens_prop);
#endif

static int  lb94_isgga(void) { return 1; }
static int  lb94_read(const char* conf_line);
static real lb94_energy(const DftDensProp* dens_prop);
static void lb94_first(FirstFuncDrv *ds, real factor,
                       const DftDensProp* dens_prop);
static void lb94_second(SecondFuncDrv *ds, real factor,
                        const DftDensProp* dens_prop);
                                                                                                                    
static void lb94_third(ThirdFuncDrv *ds, real factor,
                       const DftDensProp* dens_prop);
                                                                                                                    
#ifdef FOURTH_ORDER_DERIVATIVES
static void lb94_fourth(FourthFuncDrv *ds, real factor,
                        const DftDensProp* dens_prop);
#endif


Functional LB94PotFunctional = {"LB94pot",  /* name */
                             lb94pot_isgga,  /* gga-corrected */
                             lb94pot_read,   /* set bloody common blocks */
                             NULL,         /* reporter */
                             lb94pot_energy, 
                             lb94pot_first,
                             lb94pot_second,
                             lb94pot_third,
#ifdef FOURTH_ORDER_DERIVATIVES
                             lb94pot_fourth
#endif
};

Functional LB94HighFunctional = {"LB94high",  /* name */
                             lb94high_isgga,  /* gga-corrected */
                             lb94high_read,   /* set bloody common blocks */
                             NULL,         /* reporter */
                             lb94high_energy,
                             lb94high_first,
                             lb94high_second,
                             lb94high_third,
#ifdef FOURTH_ORDER_DERIVATIVES
                             lb94high_fourth
#endif
};

Functional LB94Functional = {"LB94",  /* name */
                             lb94_isgga,  /* gga-corrected */
                             lb94_read,   /* set bloody common blocks */
                             NULL,        /* reporter */
                             lb94_energy,
                             lb94_first,
                             lb94_second,
                             lb94_third,
#ifdef FOURTH_ORDER_DERIVATIVES
                             lb94_fourth
#endif
};

/* IMPLEMENTATION PART */

static int
lb94pot_read(const char* conf_line)
{
    dft_set_hf_weight(0.0);
    return 1;
}

/* lb94_energy:
   note that in reality E_LB94 = E_LB94,alpha + E_LB94,beta
   i.e the is linear in alpha and beta densities.

   lb94 threshold is needed to avoid numerical problems on 0/0
   divisions.  The problems are small but it is better to be on the
   safe side.
*/
static const real LB94_THRESHOLD = 1e-14;
static const real BETA = 0.050;
/* BETA:
 * note that BETA parameter in LB94 functional is
 * different to BETA in Becke(88) functional.
 */

static real
lb94pot_energy(const DftDensProp* dp)
{
 /* No energy functional for LB94 */
 return (0.0);
}

/* LB94 functional derivatives
 * 
 * Following partitioning of LB94 has been used:
 * lb94= -BETA*rho^(1/3)*F( X(rho, grad) )
 *
 * where auxaliary functions X, F 
 * are defined as:
 *  
 * X(rho,grad) = grad/rho^(4/3)
 * F(X) = (X^2) / 1 + 3*BETA*X*ArcSinh(X)
 *
 * Derivatives were then expressed in terms of derivatives of
 * auxaliary functions.
 *
 * **********************************************************
 * Variable naming convention:
 *
 * NAME       :   Variable name
 * NAME_XY    :   X-th partial derivative with respect to rho
 * 	 	  Y-th partial derivative with respect to grad
 * 	 	  of Variable (function) NAME
 * 	 	  
 * NAME_pX    :   Xth power of NAME  (NAME^X)
 * NAME_pXfY  :   NAME^(X/Y)
 *
 * NAME_X'X'  :   X-th partial derivative with respect to X function
 * 	 	  
 * Exception  :   ASINH_X does not follow these rules
 * 		  ASINH_X = asinh(X)
 * **********************************************************
 */

inline double
max(double a, double b) { return a>b ? a : b; }

static void
lb94pot_first(FirstFuncDrv *ds, real factor, const DftDensProp* dp)
{
    real rho    = dp->rhoa + dp->rhob;
    real rho13 = pow(rho, 1.0/3.0);
    real grad = dp->grada + dp->gradb;
    real X = grad/max(rho*rho13,1e-13);
    real X_p2   = X*X;
                                                                                
    real vx = -BETA*rho13*X_p2/
        (1+3*BETA*X*asinh(X));
                                                                                
    ds->df1000 += vx*factor;
    ds->df0100 += vx*factor;

}

static void
lb94pot_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp)
{
}

 
static void
lb94pot_third(ThirdFuncDrv *ds, real factor, const DftDensProp* dp) 
{
}

#ifdef FOURTH_ORDER_DERIVATIVES  
static void
lb94pot_fourth(FourthFuncDrv *ds, real factor, const DftDensProp* dp) 
{
}
#endif

static int
lb94high_read(const char* conf_line)
{
    dft_set_hf_weight(0.0);
    return 1;
}

static real
lb94high_energy(const DftDensProp* dp)
{
 /* No energy functional for LB94 */
 return (0.0);
}

static void
lb94high_first(FirstFuncDrv *ds, real factor, const DftDensProp* dp)
{
/* No potential here, only higher order derivatives */
}

static void
lb94high_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp)
{
real ro = dp->rhoa + dp->rhob;
real gro = dp->grada + dp->gradb;

real vx, vx_10, vx_01;
real ro_p1f3, ro_p2f3, ro_p4f3 , ro_p2, ro_p7f3;
real X, X_p2;
real X_10, X_01;
real SU1X2, SQ1X2;
real ASINH_X;
real ComDenom, ComDenom_p2;
real F, F_1X;

const real B = BETA;
/* Powers of ro */
        ro_p1f3=pow(ro,1.0/3.0);
        ro_p2f3=ro_p1f3*ro_p1f3;
        ro_p4f3=ro*ro_p1f3;
        ro_p2=ro*ro;
        ro_p7f3=ro_p2*ro_p1f3;
/* X auxaliary function and its powers */
        X=gro/ro_p4f3;
                                                                                                 
        X_p2=X*X;
/* Derivatives of  X */
        X_10=-(4.0/3.0)*gro/ro_p7f3;
        X_01=1.0/ro_p4f3;
/* F aux function */
/* (1+ X*X) and its powers (often used term) */
        SU1X2=1.0 + X_p2;
        SQ1X2=sqrt(SU1X2);
                                                                                                 
/* asinh(X) and its powers */
        ASINH_X=asinh(X);
                                                                                                 
/* ComDenom - frequently used Common Denominator from LB94 energy expression
 * and its powers */
        ComDenom=1.0+ 3.0*B*X*ASINH_X;
        ComDenom_p2=ComDenom*ComDenom;
/* F(X) auxaliary function */
        F=X_p2/ComDenom;
                                                                                                 
/* F(X) derivatives with respect to X */
        F_1X=X*(2.0 - (3.0*X_p2*B/SQ1X2) + (ComDenom - 1.0) )/ ComDenom_p2;

/* derivatives of potential */
        vx    = -BETA*(ro_p1f3*F);

        vx_10 = -BETA*(F/(3*ro_p2f3) + F_1X*X_10*ro_p1f3);
        vx_01 = -BETA*F_1X*X_01*ro_p1f3;

/* output derivatives */
       ds->df1000 += factor*vx;
       ds->df0100 += factor*vx;

       ds->df2000 += factor*vx_10;
       ds->df0200 += factor*vx_10;
       ds->df1100 += factor*vx_10;
       ds->df1010 += factor*vx_01;
       ds->df1001 += factor*vx_01;
       ds->df0110 += factor*vx_01;
       ds->df0101 += factor*vx_01;
       

}

 
static void
lb94high_third(ThirdFuncDrv *ds, real factor, const DftDensProp* dp) {

/* Declarations */
    real  ASINH_X;
    real  B, B_p2;
    real  ComDenom;
    real  ComDenom_p2;
    real  ComDenom_p3;
    real  F;
    real  F_1X;
    real  F_2X;
    real  SQ1X2;
    real  SU1X2;
    real  SU1X2_p2, SU1X2_p3;
    real  SU1X2_p3f2;
    real  X;
    real  X_01;
    real  X_10;
    real  X_10_p2, X_01_p2;
    real  X_11;
    real  X_20;
    real  X_p2;
    real  X_p3;
    real  X_p4;
    real  gro;
    real  ro;
    real  ro_p10f3;
    real  ro_p1f3;
    real  ro_p2;
    real  ro_p2f3;
    real  ro_p5f3;
    real  ro_p4f3;
    real  ro_p7f3;
    real  ro_p8f3;
    real  vx, vx_10, vx_01, vx_20, vx_11, vx_02;

/* Setting up Beta (In the following, BETA is called just 'B' ) and its powers */
    B=BETA;
    B_p2=B*B;

/* Setting up ro(density) and gro(gradient of density) 
 * This is done since the code reuses variables for alpha and beta 
 * densities  */
        ro = dp->rhoa + dp->rhob;
        gro = dp->grada + dp->gradb;

/* Powers of ro */
        ro_p1f3=pow(ro,1.0/3.0);
        ro_p2f3=ro_p1f3*ro_p1f3;
        ro_p4f3=ro*ro_p1f3;
        ro_p5f3=ro*ro_p2f3;
        ro_p2=ro*ro;
        ro_p7f3=ro_p2*ro_p1f3;
        ro_p8f3=ro_p7f3*ro_p1f3;
        ro_p10f3=ro_p7f3*ro;

/* X auxaliary function and its powers */
        X=gro/ro_p4f3;

        X_p2=X*X;
        X_p3=X_p2*X;
        X_p4=X_p3*X;

/* (1+ X*X) and its powers (often used term) */
        SU1X2=1.0 + X_p2;
        SQ1X2=sqrt(SU1X2);
        SU1X2_p3f2=SU1X2*SQ1X2;
        SU1X2_p2 = SU1X2*SU1X2;
        SU1X2_p3 = SU1X2_p2*SU1X2;

/* asinh(X) and its powers */
        ASINH_X=asinh(X);

/* ComDenom - frequently used Common Denominator from LB94 energy expression
 * and its powers */
        ComDenom=1.0+ 3.0*B*X*ASINH_X;
        ComDenom_p2=ComDenom*ComDenom;
        ComDenom_p3=ComDenom_p2*ComDenom;

/* Derivatives of X and powers of derivatives */
        X_10=-(4.0/3.0)*gro/ro_p7f3;
        X_20=(28.0/9.0)*gro/ro_p10f3;
        X_01=1.0/ro_p4f3;
        X_11=-(4.0/3.0)/ro_p7f3;

        X_10_p2=X_10*X_10;
        X_01_p2=X_01*X_01;

/* F(X) auxaliary function */
        F=X_p2/ComDenom;

/* F(X) derivatives with respect to X */
        F_1X=X*(2.0 - (3.0*X_p2*B/SQ1X2) + (ComDenom - 1.0) )/ ComDenom_p2;

        F_2X=(2.0*SQ1X2 + 2*X_p2*(SQ1X2 -9.0*B) + 3*X_p4*B*(-5.0 +6.0*B*SQ1X2) - 
              9.0*X_p3*(2+X_p2)*B_p2*ASINH_X ) / (SU1X2_p3f2*ComDenom_p3);

/* derivatives of potential */
        vx    = -BETA*ro_p1f3*F;

        vx_10 = -BETA*(F/(3*ro_p2f3) + F_1X*X_10*ro_p1f3);
        vx_01 = -BETA*F_1X*X_01*ro_p1f3;

        vx_20 = -BETA*((-2.0*F)/(9.0*ro_p5f3) + (2.0*F_1X*X_10)/(3.0*ro_p2f3) +
                F_2X*X_10_p2*ro_p1f3 + F_1X*X_20*ro_p1f3);

        vx_11 = -BETA*((F_1X*X_01)/(3.0*ro_p2f3) + F_2X*X_01*X_10*ro_p1f3 +
                F_1X*X_11*ro_p1f3);

        vx_02 = -BETA*(F_2X*X_01_p2*ro_p1f3);

/* output derivatives expressions */
        ds->df1000 += factor*vx;
        ds->df0100 += factor*vx;

        ds->df2000 += factor*vx_10;
        ds->df0200 += factor*vx_10;
        ds->df1100 += factor*vx_10;
        ds->df1010 += factor*vx_01;
        ds->df1001 += factor*vx_01;
        ds->df0110 += factor*vx_01;
        ds->df0101 += factor*vx_01;

        
        ds->df3000 += factor*vx_20;
        ds->df2100 += factor*vx_20;
        ds->df1200 += factor*vx_20;
        ds->df0300 += factor*vx_20;

        ds->df2010 += factor*vx_11;
        ds->df0210 += factor*vx_11;
        ds->df1110 += factor*vx_11;
        ds->df2001 += factor*vx_11;
        ds->df0201 += factor*vx_11;
        ds->df1101 += factor*vx_11;

        ds->df1020 += factor*vx_02;
        ds->df1002 += factor*vx_02;
        ds->df1011 += factor*vx_02;
        ds->df0120 += factor*vx_02;
        ds->df0102 += factor*vx_02;
        ds->df0111 += factor*vx_02;

}

#ifdef FOURTH_ORDER_DERIVATIVES  
static void
lb94high_fourth(FourthFuncDrv *ds, real factor, const DftDensProp* dp) {

/* Declarations */
    real  ASINH_X, ASINH_X_p2;
    real  B, B_p2;
    real  ComDenom;
    real  ComDenom_p2;
    real  ComDenom_p3;
    real  ComDenom_p4;
    real  F;
    real  F_1X;
    real  F_2X;
    real  F_3X;
    real  SQ1X2;
    real  SU1X2;
    real  SU1X2_p2;
    real  SU1X2_p3;
    real  SU1X2_p3f2;
    real  SU1X2_p5f2;
    real  X;
    real  X_01;
    real  X_01_p2, X_01_p3;
    real  X_10;
    real  X_10_p2, X_10_p3;
    real  X_11;
    real  X_20, X_30, X_21;
    real  X_p2;
    real  X_p3;
    real  X_p4, X_p5, X_p6;
    real  gro;
    real  ro;
    real  ro_p10f3;
    real  ro_p13f3;
    real  ro_p1f3;
    real  ro_p2;
    real  ro_p2f3;
    real  ro_p5f3;
    real  ro_p4f3;
    real  ro_p7f3;
    real  ro_p8f3;
    real  vx, vx_10, vx_01, vx_20, vx_11, vx_02, vx_30, vx_21, vx_12, vx_03;

/* Setting up Beta (In the following, BETA is called just 'B' ) and its powers */
    B=BETA;
    B_p2=B*B;

/* Setting up ro(density) and gro(gradient of density) 
 * This is done since the code reuses variables for alpha and beta 
 * densities  */
        ro = dp->rhoa + dp->rhob;
        gro = dp->grada + dp->gradb;

/* Powers of ro */
        ro_p1f3=pow(ro,1.0/3.0);
        ro_p2f3=ro_p1f3*ro_p1f3;
        ro_p4f3=ro*ro_p1f3;
        ro_p5f3=ro*ro_p2f3;
        ro_p2=ro*ro;
        ro_p7f3=ro_p2*ro_p1f3;
        ro_p8f3=ro_p7f3*ro_p1f3;
        ro_p10f3=ro_p7f3*ro;
        ro_p13f3=ro_p10f3*ro;

/* X auxaliary function and its powers */
        X=gro/ro_p4f3;

        X_p2=X*X;
        X_p3=X_p2*X;
        X_p4=X_p3*X;
        X_p5=X_p4*X;
        X_p6=X_p5*X;

/* (1+ X*X) and its powers (often used term) */
        SU1X2=1.0 + X_p2;
        SQ1X2=sqrt(SU1X2);
        SU1X2_p3f2=SU1X2*SQ1X2;
        SU1X2_p5f2=SU1X2_p3f2*SU1X2;
        SU1X2_p2 = SU1X2*SU1X2;
        SU1X2_p3 = SU1X2_p2*SU1X2;

/* asinh(X) and its powers */
        ASINH_X=asinh(X);
        ASINH_X_p2=ASINH_X*ASINH_X;

/* ComDenom - frequently used Common Denominator from LB94 energy expression
 * and its powers */
        ComDenom=1.0+ 3.0*B*X*ASINH_X;
        ComDenom_p2=ComDenom*ComDenom;
        ComDenom_p3=ComDenom_p2*ComDenom;
        ComDenom_p4=ComDenom_p3*ComDenom;

/* Derivatives of X and powers of derivatives */
        X_10=-(4.0/3.0)*gro/ro_p7f3;
        X_20=(28.0/9.0)*gro/ro_p10f3;
        X_30=-(280.0/27.0)*gro/ro_p13f3;
        X_01=1.0/ro_p4f3;
        X_11=-(4.0/3.0)/ro_p7f3;
        X_21=(28.0/9.0)/ro_p10f3;

        X_10_p2=X_10*X_10;
        X_10_p3=X_10_p2*X_10;
        
	X_01_p2=X_01*X_01;
        X_01_p3=X_01_p2*X_01;

/* F(X) auxaliary function */
        F=X_p2/ComDenom;

/* F(X) derivatives with respect to X */
        F_1X=X*(2.0 - (3.0*X_p2*B/SQ1X2) + (ComDenom - 1.0) )/ ComDenom_p2;

        F_2X=(2.0*SQ1X2 + 2*X_p2*(SQ1X2 -9.0*B) + 3*X_p4*B*(-5.0 +6.0*B*SQ1X2) - 
              9.0*X_p3*(2+X_p2)*B_p2*ASINH_X ) / (SU1X2_p3f2*ComDenom_p3);

        F_3X=(3.0*B*(-X*(18.0 + 54.0*X_p6*B_p2 + X_p2*(26.0 - 72.0*B*SQ1X2) + 
              X_p4*(11.0 -54.0*B*SQ1X2 + 54.0*B_p2)) + 6.0*(-SQ1X2 - 
              2.0*X_p2*SQ1X2 + 4.0*X_p6*B + X_p4*(-SQ1X2 + 7.0*B + 9.0*SQ1X2*B_p2))*ASINH_X +
              9.0*X_p5*(4.0 + X_p2)*B_p2*ASINH_X_p2)) / (SU1X2_p5f2*ComDenom_p4);

/* derivatives of potential */
        vx    = -BETA*ro_p1f3*F;

        vx_10 = -BETA*(F/(3*ro_p2f3) + F_1X*X_10*ro_p1f3);
        vx_01 = -BETA*F_1X*X_01*ro_p1f3;

        vx_20 = -BETA*((-2.0*F)/(9.0*ro_p5f3) + (2.0*F_1X*X_10)/(3.0*ro_p2f3) +
                F_2X*X_10_p2*ro_p1f3 + F_1X*X_20*ro_p1f3);

        vx_11 = -BETA*((F_1X*X_01)/(3.0*ro_p2f3) + F_2X*X_01*X_10*ro_p1f3 +
                F_1X*X_11*ro_p1f3);

        vx_02 = -BETA*(F_2X*X_01_p2*ro_p1f3);

        vx_30 = -BETA*((10.0*F)/(27.0*ro_p8f3) - (2.0*F_1X*X_10)/(3.0*ro_p5f3) +
                (F_2X*X_10_p2)/ro_p2f3 + (F_1X*X_20)/ro_p2f3 +
                F_3X*X_10_p3*ro_p1f3 + 3.0*F_2X*X_10*X_20*ro_p1f3 +
                F_1X*X_30*ro_p1f3);

        vx_21 = -BETA*((-2.0*F_1X*X_01)/(9.0*ro_p5f3) + (2.0*F_2X*X_01*X_10)/
                (3.0*ro_p2f3) + (2.0*F_1X*X_11)/(3.0*ro_p2f3) +
                F_3X*X_01*X_10_p2*ro_p1f3 + 2.0*F_2X*X_10*X_11*ro_p1f3 +
                F_2X*X_01*X_20*ro_p1f3 + F_1X*X_21*ro_p1f3);

        vx_12 = -BETA*((F_2X*X_01_p2)/(3.0*ro_p2f3) + F_3X*X_01_p2*X_10*ro_p1f3 +
                2.0*F_2X*X_01*X_11*ro_p1f3);

        vx_03 = -BETA*F_3X*X_01_p3*ro_p1f3;

/* output derivatives expressions */
        ds->df10000 += factor*vx;
        ds->df01000 += factor*vx;

        ds->df20000 += factor*vx_10;
        ds->df02000 += factor*vx_10;
        ds->df11000 += factor*vx_10;
        ds->df10100 += factor*vx_01;
        ds->df10010 += factor*vx_01;
        ds->df01100 += factor*vx_01;
        ds->df01010 += factor*vx_01;

        
        ds->df30000 += factor*vx_20;
        ds->df21000 += factor*vx_20;
        ds->df12000 += factor*vx_20;
        ds->df03000 += factor*vx_20;

        ds->df20100 += factor*vx_11;
        ds->df02100 += factor*vx_11;
        ds->df11100 += factor*vx_11;
        ds->df20010 += factor*vx_11;
        ds->df02010 += factor*vx_11;
        ds->df11010 += factor*vx_11;

        ds->df10200 += factor*vx_02;
        ds->df10020 += factor*vx_02;
        ds->df10110 += factor*vx_02;
        ds->df01200 += factor*vx_02;
        ds->df01020 += factor*vx_02;
        ds->df01110 += factor*vx_02;

        ds->df40000 += factor*vx_30;
        ds->df04000 += factor*vx_30;
        ds->df31000 += factor*vx_30;
        ds->df22000 += factor*vx_30;
        ds->df13000 += factor*vx_30;

        ds->df30100 += factor*vx_21;
        ds->df21100 += factor*vx_21;
        ds->df12100 += factor*vx_21;
        ds->df03100 += factor*vx_21;
        ds->df30010 += factor*vx_21;
        ds->df21010 += factor*vx_21;
        ds->df12010 += factor*vx_21;
        ds->df03010 += factor*vx_21;

        ds->df20200 += factor*vx_12;
        ds->df11200 += factor*vx_12;
        ds->df02200 += factor*vx_12;
        ds->df20110 += factor*vx_12;
        ds->df11110 += factor*vx_12;
        ds->df02020 += factor*vx_12;
        ds->df02110 += factor*vx_12;
        ds->df20020 += factor*vx_12;
        ds->df11020 += factor*vx_12;

        ds->df10300 += factor*vx_03;
        ds->df10210 += factor*vx_03;
        ds->df10120 += factor*vx_03;
        ds->df10030 += factor*vx_03;
        ds->df01300 += factor*vx_03;
        ds->df01210 += factor*vx_03;
        ds->df01120 += factor*vx_03;
        ds->df01030 += factor*vx_03;


  

}
#endif

#if 0
Functional *lb94energy_func = NULL;
void fort_print(const char* format, ...);
void dftreset_(void);
int dftinput_(const char* line, int * res, int len);
real dft_get_hf_weight(void);

static int
lb94mix_read(const char* conf_line)
{
    /* Save current selected functional */
    Functional *save = selected_func;
    int res;
    
    /* Read and set lb94energy functional (via global selected_func) */
    dft_set_hf_weight(0.0);
    dftreset_();
    dftinput_(&conf_line[1], &res, 80);    
    lb94energy_func = selected_func;

    /* quit if no functional selected */
    if (! res)
        fort_print
        ("--- ERROR:  LB94mix: Additional functional for energy must be selected! ---");

    /* Check weather lb94energy functional is hybrid functional */
    if (dft_get_hf_weight() > 0.0 ) {
        fort_print
        ("--- ERROR:  LB94mix: It is not possible to select hybrid functional! ---");
        res = 0;
    }

    /* Restore functional */
    selected_func = save;
    dft_set_hf_weight(0.0);
    return res;
}

static void
lb94mix_report(void)
{
     fort_print
     ("    LB94mix: Functional used for energy: %s", lb94energy_func->name);
     if (lb94energy_func->report) lb94energy_func->report();
}
#endif

static int
lb94_read(const char* conf_line) 
{ 
     dft_set_hf_weight(0.0);
     return 1;
}


static real
lb94_energy(const DftDensProp* dp)
{
 /* return (lb94energy_func->func(dp)); */
    return (
             DiracFunctional.func(dp) +
             VWNFunctional.func  (dp)
    );
}

static void                                                                                                      
lb94_first(FirstFuncDrv *ds, real factor, const DftDensProp* dp)
{
    DiracFunctional.first   (ds, factor, dp);
    VWNFunctional.first     (ds, factor, dp);
    LB94PotFunctional.first (ds, factor, dp);
}

static void
lb94_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp)
{
    DiracFunctional.second   (ds, factor, dp);
    VWNFunctional.second     (ds, factor, dp);
    /* LB94HighFunctional.second (ds, factor, dp); */
}
                                                                                                                    
                                                                                                                    
static void
lb94_third(ThirdFuncDrv *ds, real factor, const DftDensProp* dp)
{
    DiracFunctional.third   (ds, factor, dp);
    VWNFunctional.third     (ds, factor, dp);
    /* LB94HighFunctional.third (ds, factor, dp); */
}
                                                                                                                    
#ifdef FOURTH_ORDER_DERIVATIVES
static void
lb94_fourth(FourthFuncDrv *ds, real factor, const DftDensProp* dp)
{
    DiracFunctional.fourth   (ds, factor, dp);
    VWNFunctional.fourth     (ds, factor, dp);
    /* LB94HighFunctional.fourth (ds, factor, dp); */
}
#endif
