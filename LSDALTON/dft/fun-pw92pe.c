/* fun-pw92.c:
   implementation of the electron-gas correlation energy
   by J. P. Perdew and Y. Wang;
   J.P.Perdew,Y. Wang; Phys. Rew. B; 40, 13244, (1992) [Ref 1]
   
   The PW92 functional and its derivatives up to fourth order.
   Implementation (c)by B. Jansik, brano@theochem.kth.se, jan 2003

   NOTE: Improvement over VWN:
   While we confirm the practical accuracy of the VWN and PZ reprezentations,
   we eliminate some minor problems with these forms. [1]
   See ref [1] for details.

   NOTE: This functional returns energy and derivatives of
         'per electron' PW92 functional as given in ref. [1]

*/

#if !defined(SYS_DEC)
/* XOPEN compliance is missing on old Tru64 4.0E Alphas and pow() prototype
 * is not specified. */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1
#endif

#include <math.h>
#include <stddef.h>

#define __CVERSION__

#include "lsdalton_functionals.h"

#define PW92THR 1e-14 /* PW92 Threshold, denities below this are ignored */

/* INTERFACE PART */
static integer  pw92_isgga(void) { return 0; }
static integer  pw92_read(const char* conf_line, real *hfweight);
static void pw92_report(integer lupri);
static real pw92_energy(const DftDensProp* dens_prop);
static void pw92_first(FirstFuncDrv *ds, real factor,
                       const DftDensProp* dens_prop);
static void pw92_second(SecondFuncDrv *ds, real factor,
                        const DftDensProp* dens_prop);

static void pw92_third(ThirdFuncDrv *ds, real factor,
                       const DftDensProp* dens_prop);
/* #define FOURTH_ORDER_DERIVATIVES */
#ifdef FOURTH_ORDER_DERIVATIVES
static void pw92_fourth(FourthFuncDrv *ds, real factor,
                        const DftDensProp* dens_prop);
#endif
Functional PW92peFunctional = {
    "PW92pe",      /* name */
    pw92_isgga,  /* gga-corrected */
    pw92_read,   /* set bloody common blocks */
    pw92_report,         /* reporter */
    pw92_energy, 
    pw92_first,
    pw92_second,
    pw92_third
#ifdef FOURTH_ORDER_DERIVATIVES
    ,pw92_fourth
#endif
};

/* Parameters for PW92 functional.
 * For the parameter values see J. P. Perdew and Y. Wang, 
 * Phys. Rev. B. 45, 13244 (1992). Table I 
 */
struct GfitParams_ {
    real  p;
    real  A;
    real  alpha1;
    real  beta1;
    real  beta2;
    real  beta3;
    real  beta4;
};

typedef struct GfitParams_ GfitParams;
/* Declarations of aux functions for pw92 functional */
real Gfit    (real *pRS, const GfitParams *COEF);
real Gfit_1RS(real *pRS, const GfitParams *COEF);
real Gfit_2RS(real *pRS, const GfitParams *COEF);
real Gfit_3RS(real *pRS, const GfitParams *COEF);
real Gfit_4RS(real *pRS, const GfitParams *COEF);



/* in RPA approx */         /*   p,        A,   alpha1,  beta1,  beta2,   beta3,   beta4 */
const GfitParams E0_rpa   = { 0.75, 0.031091, 0.082477, 5.1486, 1.6483, 0.23647, 0.20614 };
const GfitParams E1_rpa   = { 0.75, 0.015545, 0.035374, 6.4869, 1.3083, 0.15180, 0.082349 };
const GfitParams minA_rpa = { 1.0,  0.016887, 0.028829, 10.357, 3.6231, 0.47990, 0.12279 };

/* beyond RPA approx */     /*   p,          A,   alpha1,   beta1,  beta2,   beta3,    beta4 */
const GfitParams E0_brpa  = { 1.00, 0.03109070, 0.213700, 7.59570, 3.5876, 1.63820, 0.492940 };
const GfitParams E1_brpa  = { 1.00, 0.01554535, 0.205480, 14.1189, 6.1977, 3.36620, 0.625170 };
const GfitParams minA_brpa= { 1.00, 0.01688690, 0.111250, 10.3570, 3.6231, 0.88026, 0.496710 };

/* IMPLEMENTATION PART */

static integer
pw92_read(const char* conf_line, real *hfweight)
{
/*  dft_set_hf_weight(0.0); */
    return 1;
}


static void
pw92_report(integer lupri)
{
#ifndef TESTER
  lsfort_print(lupri,"\n     WARNING: PW92pe functional does not give total energy,");
   lsfort_print(lupri,"             it provides energy per electron instead!");
   lsfort_print(lupri,"             Functional is only for internal purposes!!");
   lsfort_print(lupri,"             Use PW92 instead! (PW92=num_of_electrons*PW92pe)\n");
   exit(1);
#endif
}

static real
pw92_energy(const DftDensProp* dp)
{
/* Declarations */

    real  E0;
    real  E1;
    real  F;
    real  RS;
    real  ZETA;
    real  ZETA_p2;
    real  ZETA_p4;
    real  minA;
    real  roa;
    real  rob;


/* Numbers and constants */

    const real Fdenom=(pow(2.0,(4.0/3.0)) - 2.0);
    const real F0=1.709921;

/* Setting up alpha and beta densities roa and rob */
    roa=dp->rhoa;
    rob=dp->rhob;

/* Ignore small densities */
   if((roa+rob)<PW92THR) return(0.0);

/* ZETA and RS */
    ZETA=(roa-rob)/(roa+rob);
    RS=pow((3.0/(4.0*M_PI*(roa+rob))),(1.0/3.0));

/* Powers of ZETA */
    ZETA_p2=ZETA*ZETA;
    ZETA_p4=ZETA_p2*ZETA_p2;

/* Spin function F */

    F=(  pow((1.0+ZETA),(4.0/3.0)) + pow((1.0-ZETA),(4.0/3.0)) - 2.0)/Fdenom;

/* E0, F1, minA according to G fit */
    E0=Gfit(&RS, &E0_brpa);
    E1=Gfit(&RS, &E1_brpa);
    minA=Gfit(&RS, &minA_brpa);

/* Electron gas Correlation energy
 * epsilon(RS,ZETA) */
    return (E0 - minA*(F/F0)*(1.0-ZETA_p4) + (E1 - E0)*F*ZETA_p4);

}

static void
pw92_first(FirstFuncDrv *ds, real factor, const DftDensProp* dp)
{
/* Declarations */
    real  E0;
    real  E0_1RS;
    real  E1;
    real  E1_1RS;
    real  F;
    real  F_1ZETA;
    real  M1ZETA;
    real  M1ZETA_p1f3;
    real  M1ZETA_p4f3;
    real  P1ZETA;
    real  P1ZETA_p1f3;
    real  P1ZETA_p4f3;
    real  RS;
    real  RS_01;
    real  RS_10;
    real  ZETA;
    real  ZETA_01;
    real  ZETA_10;
    real  ZETA_p2;
    real  ZETA_p3;
    real  ZETA_p4;
    real  minA;
    real  minA_1RS;
    real  roa;
    real  roa_rob;
    real  roa_rob_p2;
    real  roa_rob_p1f3;
    real  roa_rob_p4f3;
    real  rob;

/* Numbers and constants */

    const real  Fdenom=(pow(2.0,(4.0/3.0)) - 2.0);
    const real  iF0=1.0/1.709921;

/* Numbers in expressions */
    const real  nRS_10=1.0/(pow(6.0,(2.0/3.0))*pow(M_PI,(1.0/3.0)));

/* Setting up alpha and beta densities roa and rob */
    roa=dp->rhoa;
    rob=dp->rhob;

/* Ignore small densities */
   if((roa+rob)<PW92THR) return;

/* (roa + rob) and its powers */
    roa_rob=roa+rob;
    roa_rob_p2=roa_rob*roa_rob;
    roa_rob_p1f3=pow(roa_rob,(1.0/3.0));
    roa_rob_p4f3=roa_rob*roa_rob_p1f3;

/* Auxaliary functions RS and ZETA */

    ZETA=(roa-rob)/(roa+rob);
    RS=pow((3.0/(4.0*M_PI*(roa+rob))),(1.0/3.0));


/* Powers of ZETA */
    ZETA_p2=ZETA*ZETA;
    ZETA_p3=ZETA_p2*ZETA;
    ZETA_p4=ZETA_p2*ZETA_p2;

    P1ZETA=1.0 + ZETA;
    M1ZETA=1.0 - ZETA;

    P1ZETA_p1f3=pow(P1ZETA,(1.0/3.0));
    P1ZETA_p4f3=P1ZETA_p1f3*P1ZETA;
    M1ZETA_p1f3=pow(M1ZETA,(1.0/3.0));
    M1ZETA_p4f3=M1ZETA_p1f3*M1ZETA;

/* Derivatives od RS */

    RS_10=-(1.0/roa_rob_p4f3)*nRS_10;
    RS_01=RS_10;

/* Derivatives of ZETA */

    ZETA_10=-((roa - rob)/roa_rob_p2) + (1.0/roa_rob);
    ZETA_01=-((roa - rob)/roa_rob_p2) - (1.0/roa_rob);

/* Derivatives of F (with respect to ZETA) */
    F=(P1ZETA_p4f3 + M1ZETA_p4f3 - 2.0)/Fdenom;

    F_1ZETA=((4.0/3.0)*P1ZETA_p1f3 - (4.0/3.0)*M1ZETA_p1f3)/Fdenom;

/* Functions E0, E1, minA */

    E0=Gfit(&RS, &E0_brpa);
    E1=Gfit(&RS, &E1_brpa);
    minA=Gfit(&RS, &minA_brpa);

/* Derivatives of  E0, E1, minA (with respect to RS) */

    E0_1RS=Gfit_1RS(&RS, &E0_brpa);

    E1_1RS=Gfit_1RS(&RS, &E1_brpa);

    minA_1RS=Gfit_1RS(&RS, &minA_brpa);


/* Derivatives */

    ds->df1000 += factor*( 
                           E0_1RS*RS_10 - 
                           iF0*F*(1.0-ZETA_p4)*minA_1RS*RS_10 + 
                           F*ZETA_p4*(-E0_1RS*RS_10 + E1_1RS*RS_10) + 
                           4.0*(-E0 + E1)*F*ZETA_p3*ZETA_10 + 
                           4.0*iF0*F*ZETA_p3*minA*ZETA_10 + 
                           (-E0+E1)*ZETA_p4*F_1ZETA*ZETA_10 - 
                           iF0*(1.0-ZETA_p4)*minA*F_1ZETA*ZETA_10
	);

    ds->df0100 += factor*( 
                           E0_1RS*RS_01 - 
                           iF0*F*(1.0-ZETA_p4)*minA_1RS*RS_01 + 
                           F*ZETA_p4*(-E0_1RS*RS_01 + E1_1RS*RS_01) + 
                           4.0*(-E0 + E1)*F*ZETA_p3*ZETA_01 + 
                           4.0*iF0*F*ZETA_p3*minA*ZETA_01 + 
                           (-E0+E1)*ZETA_p4*F_1ZETA*ZETA_01 - 
                           iF0*(1.0-ZETA_p4)*minA*F_1ZETA*ZETA_01
	);

}

static void
pw92_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp)
{
/* Declarations */
    real  E0, E0_1RS, E0_2RS;
    real  E1, E1_1RS, E1_2RS;
    real  F,  F_1ZETA, F_2ZETA;
    real  M1ZETA, M1ZETA_p1f3, M1ZETA_p2f3, M1ZETA_p4f3;
    real  P1ZETA, P1ZETA_p1f3, P1ZETA_p2f3, P1ZETA_p4f3;
    real  RS, RS_01, RS_01_p2, RS_02, RS_10, RS_10_p2, RS_11, RS_20;
    real  ZETA;
    real  ZETA_01;
    real  ZETA_01_p2;
    real  ZETA_02;
    real  ZETA_10;
    real  ZETA_10_p2;
    real  ZETA_11;
    real  ZETA_20;
    real  ZETA_p2;
    real  ZETA_p3;
    real  ZETA_p4;
    real  minA, minA_1RS, minA_2RS;
    real  roa;
    real  roa_rob;
    real  roa_rob_p1f3;
    real  roa_rob_p2;
    real  roa_rob_p3;
    real  roa_rob_p4f3;
    real  roa_rob_p7f3;
    real  rob;

/* Numbers and constants */

    const real  Fdenom=(pow(2.0,(4.0/3.0)) - 2.0);
    const real  iF0=1.0/1.709921;

/* Numbers in expressions */
    const real  nRS_10=1.0/(pow(6.0,(2.0/3.0))*pow(M_PI,(1.0/3.0)));
    const real  nRS_20=2.0*(pow((2.0/M_PI),(1.0/3.0)))/(3.0*pow(3.0,(2.0/3.0)));

/* Setting up alpha and beta densities roa and rob */
    roa=dp->rhoa;
    rob=dp->rhob;

/* Ignore small densities */
   if((roa+rob)<PW92THR) return;

/* (roa + rob) and its powers */
    roa_rob=roa+rob;
    roa_rob_p2=roa_rob*roa_rob;
    roa_rob_p3=roa_rob_p2*roa_rob;
    roa_rob_p1f3=pow(roa_rob,(1.0/3.0));
    roa_rob_p4f3=roa_rob*roa_rob_p1f3;
    roa_rob_p7f3=roa_rob_p4f3*roa_rob;

/* Auxaliary functions RS and ZETA */

    ZETA=(roa-rob)/(roa+rob);
    RS=pow((3.0/(4.0*M_PI*(roa+rob))),(1.0/3.0));


/* Powers of ZETA */
    ZETA_p2=ZETA*ZETA;
    ZETA_p3=ZETA_p2*ZETA;
    ZETA_p4=ZETA_p2*ZETA_p2;

    P1ZETA=1.0 + ZETA;
    M1ZETA=1.0 - ZETA;

    P1ZETA_p1f3=pow(P1ZETA,(1.0/3.0));
    P1ZETA_p2f3=P1ZETA_p1f3*P1ZETA_p1f3;
    P1ZETA_p4f3=P1ZETA_p1f3*P1ZETA;
    M1ZETA_p1f3=pow(M1ZETA,(1.0/3.0));
    M1ZETA_p2f3=M1ZETA_p1f3*M1ZETA_p1f3;
    M1ZETA_p4f3=M1ZETA_p1f3*M1ZETA;

/* Derivatives od RS */

    RS_10=-(1.0/roa_rob_p4f3)*nRS_10;
    RS_01=RS_10;

    RS_20=(1.0/roa_rob_p7f3)*nRS_20;
    RS_02=RS_20;
    RS_11=RS_20;

/* Powers of RS */
    RS_10_p2=RS_10*RS_10;
    RS_01_p2=RS_01*RS_01;


/* Derivatives of ZETA */

    ZETA_10=-((roa - rob)/roa_rob_p2) + (1.0/roa_rob);
    ZETA_01=-((roa - rob)/roa_rob_p2) - (1.0/roa_rob);

    ZETA_20=2.0*((roa-rob)/roa_rob_p3) - (2.0/roa_rob_p2);
    ZETA_02=2.0*((roa-rob)/roa_rob_p3) + (2.0/roa_rob_p2);
    ZETA_11=2.0*((roa-rob)/roa_rob_p3);

/* Powers of ZETA */

    ZETA_10_p2=ZETA_10*ZETA_10;
    ZETA_01_p2=ZETA_01*ZETA_01;

/* Derivatives of F (with respect to ZETA) */
    F=(P1ZETA_p4f3 + M1ZETA_p4f3 - 2.0)/Fdenom;

    F_1ZETA=((4.0/3.0)*P1ZETA_p1f3 - (4.0/3.0)*M1ZETA_p1f3)/Fdenom;

    F_2ZETA=(4.0/(9.0*Fdenom))*( (1.0/M1ZETA_p2f3) + (1.0/P1ZETA_p2f3)   );

/* Functions E0, E1, minA */

    E0=Gfit(&RS, &E0_brpa);
    E1=Gfit(&RS, &E1_brpa);
    minA=Gfit(&RS, &minA_brpa);

/* Derivatives of  E0, E1, minA (with respect to RS) */

    E0_1RS=Gfit_1RS(&RS, &E0_brpa);
    E0_2RS=Gfit_2RS(&RS, &E0_brpa);

    E1_1RS=Gfit_1RS(&RS, &E1_brpa);
    E1_2RS=Gfit_2RS(&RS, &E1_brpa);

    minA_1RS=Gfit_1RS(&RS, &minA_brpa);
    minA_2RS=Gfit_2RS(&RS, &minA_brpa);

/* Derivatives */

    ds->df1000 += factor*( 
                           E0_1RS*RS_10 - 
                           iF0*F*(1.0-ZETA_p4)*minA_1RS*RS_10 + 
                           F*ZETA_p4*(-E0_1RS*RS_10 + E1_1RS*RS_10) + 
                           4.0*(-E0 + E1)*F*ZETA_p3*ZETA_10 + 
                           4.0*iF0*F*ZETA_p3*minA*ZETA_10 + 
                           (-E0+E1)*ZETA_p4*F_1ZETA*ZETA_10 - 
                           iF0*(1.0-ZETA_p4)*minA*F_1ZETA*ZETA_10
	);

    ds->df0100 += factor*( 
                           E0_1RS*RS_01 - 
                           iF0*F*(1.0-ZETA_p4)*minA_1RS*RS_01 + 
                           F*ZETA_p4*(-E0_1RS*RS_01 + E1_1RS*RS_01) + 
                           4.0*(-E0 + E1)*F*ZETA_p3*ZETA_01 + 
                           4.0*iF0*F*ZETA_p3*minA*ZETA_01 + 
                           (-E0+E1)*ZETA_p4*F_1ZETA*ZETA_01 - 
                           iF0*(1.0-ZETA_p4)*minA*F_1ZETA*ZETA_01
	);

    ds->df2000 += factor*( 
                           E0_2RS*RS_10_p2 + 
                           iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_10_p2 + 
                           8.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*RS_10*ZETA_10 + 
                           2.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*RS_10*ZETA_10 + 
                           8.0*iF0*F*ZETA_p3*minA_1RS*RS_10*ZETA_10 + 
                           2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_10*ZETA_10 + 
                           12.0*(-E0 + E1)*F*ZETA_p2*ZETA_10_p2 + 
                           12.0*iF0*F*ZETA_p2*minA*ZETA_10_p2 + 
                           8.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_10_p2 + 
                           8.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_10_p2 + 
                           (-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_10_p2 + 
                           iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_10_p2 + 
                           E0_1RS*RS_20 + 
                           iF0*F*(-1.0 + ZETA_p4)*minA_1RS*RS_20 + 
                           F*ZETA_p4*((-E0_2RS + E1_2RS)*RS_10_p2 + 
                                      (-E0_1RS + E1_1RS)*RS_20) + 
                           4.0*(-E0 + E1)*F*ZETA_p3*ZETA_20 + 
                           4.0*iF0*F*ZETA_p3*minA*ZETA_20 + 
                           (-E0 + E1)*ZETA_p4*F_1ZETA*ZETA_20 + 
                           iF0*(-1.0 + ZETA_p4)*minA*F_1ZETA*ZETA_20
	);

    ds->df0200 += factor*( 
                           E0_2RS*RS_01_p2 + 
                           iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_01_p2 + 
                           8.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*RS_01*ZETA_01 + 
                           2.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_01 + 
                           8.0*iF0*F*ZETA_p3*minA_1RS*RS_01*ZETA_01 + 
                           2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_01*ZETA_01 + 
                           12.0*(-E0 + E1)*F*ZETA_p2*ZETA_01_p2 + 
                           12.0*iF0*F*ZETA_p2*minA*ZETA_01_p2 + 
                           8.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_01_p2 + 
                           8.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_01_p2 + 
                           (-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_01_p2 + 
                           iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_01_p2 + 
                           E0_1RS*RS_02 + 
                           iF0*F*(-1.0 + ZETA_p4)*minA_1RS*RS_02 + 
                           F*ZETA_p4*((-E0_2RS + E1_2RS)*RS_01_p2 + 
                                      (-E0_1RS + E1_1RS)*RS_02) + 
                           4.0*(-E0 + E1)*F*ZETA_p3*ZETA_02 + 
                           4.0*iF0*F*ZETA_p3*minA*ZETA_02 + 
                           (-E0 + E1)*ZETA_p4*F_1ZETA*ZETA_02 + 
                           iF0*(-1.0 + ZETA_p4)*minA*F_1ZETA*ZETA_02
	);

    ds->df1100 += factor*( 
                           E0_2RS*RS_01*RS_10+iF0*F*(-1.0+ZETA_p4)*minA_2RS*RS_01*RS_10+
                           4.0*F*ZETA_p3*(-E0_1RS+E1_1RS)*ZETA_01*RS_10+
                           ZETA_p4*(-E0_1RS+E1_1RS)*F_1ZETA*ZETA_01*RS_10+
                           4.0*iF0*F*ZETA_p3*minA_1RS*ZETA_01*RS_10+
                           iF0*(-1.0+ZETA_p4)*F_1ZETA*minA_1RS*ZETA_01*RS_10+
                           4.0*F*ZETA_p3*(-E0_1RS+E1_1RS)*RS_01*ZETA_10+
                           ZETA_p4*(-E0_1RS+E1_1RS)*F_1ZETA*RS_01*ZETA_10+
                           4.0*iF0*F*ZETA_p3*minA_1RS*RS_01*ZETA_10+
                           iF0*(-1.0+ZETA_p4)*F_1ZETA*minA_1RS*RS_01*ZETA_10+12.0*(-E0+E1)
                           *F*ZETA_p2*ZETA_01*ZETA_10+12.0*iF0*F*ZETA_p2*minA*ZETA_01*ZETA_10+
                           8.0*(-E0+E1)*ZETA_p3*F_1ZETA*ZETA_01*ZETA_10+
                           8.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_01*ZETA_10+
                           (-E0+E1)*ZETA_p4*F_2ZETA*ZETA_01*ZETA_10+
                           iF0*(-1.0+ZETA_p4)*minA*F_2ZETA*ZETA_01*ZETA_10+
                           E0_1RS*RS_11+iF0*F*(-1.0+ZETA_p4)*minA_1RS*RS_11+
                           F*ZETA_p4*((-E0_2RS+E1_2RS)*RS_01*RS_10+(-E0_1RS+E1_1RS)*RS_11)+
                           4.0*(-E0+E1)*F*ZETA_p3*ZETA_11+4.0*iF0*F*ZETA_p3*minA*ZETA_11+
                           (-E0+E1)*ZETA_p4*F_1ZETA*ZETA_11+
                           iF0*(-1.0+ZETA_p4)*minA*F_1ZETA*ZETA_11
	);

}
 
static void
pw92_third(ThirdFuncDrv *ds, real factor, const DftDensProp* dp) {
/* Declarations */
    real  E0, E0_1RS,  E0_2RS,  E0_3RS;
    real  E1, E1_1RS,  E1_2RS,  E1_3RS;
    real  F,  F_1ZETA, F_2ZETA, F_3ZETA;
    real  M1ZETA, M1ZETA_p1f3, M1ZETA_p2f3, M1ZETA_p4f3, M1ZETA_p5f3;
    real  P1ZETA, P1ZETA_p1f3, P1ZETA_p2f3, P1ZETA_p4f3, P1ZETA_p5f3;
    real  RS;
    real  RS_01;
    real  RS_01_p2;
    real  RS_01_p3;
    real  RS_02;
    real  RS_03;
    real  RS_10;
    real  RS_10_p2;
    real  RS_10_p3;
    real  RS_11;
    real  RS_12;
    real  RS_20;
    real  RS_21;
    real  RS_30;
    real  ZETA;
    real  ZETA_01;
    real  ZETA_01_p2;
    real  ZETA_01_p3;
    real  ZETA_02;
    real  ZETA_03;
    real  ZETA_10;
    real  ZETA_10_p2;
    real  ZETA_10_p3;
    real  ZETA_10_p4;
    real  ZETA_11;
    real  ZETA_12;
    real  ZETA_20;
    real  ZETA_21;
    real  ZETA_30;
    real  ZETA_p2;
    real  ZETA_p3;
    real  ZETA_p4;
    real  minA, minA_1RS, minA_2RS, minA_3RS;
    real  roa;
    real  roa_rob;
    real  roa_rob_p10f3;
    real  roa_rob_p1f3;
    real  roa_rob_p2;
    real  roa_rob_p3;
    real  roa_rob_p4;
    real  roa_rob_p4f3;
    real  roa_rob_p7f3;
    real  rob;

/* Numbers and constants */

    const real  Fdenom=(pow(2.0,(4.0/3.0)) - 2.0);
    const real  iF0=1.0/1.709921;

/* Numbers in expressions */
    const real  nRS_10=1.0/(pow(6.0,(2.0/3.0))*pow(M_PI,(1.0/3.0)));
    const real  nRS_20=2.0*(pow((2.0/M_PI),(1.0/3.0)))/(3.0*pow(3.0,(2.0/3.0)));
    const real  nRS_30=14.0*(pow((2.0/M_PI),(1.0/3.0)))/(9.0*pow(3.0,(2.0/3.0)));

/* Setting up alpha and beta densities roa and rob */
    roa=dp->rhoa;
    rob=dp->rhob;

/* Ignore small densities */
   if((roa+rob)<PW92THR) return;

/* (roa + rob) and its powers */
    roa_rob=roa+rob;
    roa_rob_p2=roa_rob*roa_rob;
    roa_rob_p3=roa_rob_p2*roa_rob;
    roa_rob_p4=roa_rob_p2*roa_rob_p2;
    roa_rob_p1f3=pow(roa_rob,(1.0/3.0));
    roa_rob_p4f3=roa_rob*roa_rob_p1f3;
    roa_rob_p7f3=roa_rob_p4f3*roa_rob;
    roa_rob_p10f3=roa_rob_p7f3*roa_rob;

/* Auxaliary functions RS and ZETA */

    ZETA=(roa-rob)/(roa+rob);
    RS=pow((3.0/(4.0*M_PI*(roa+rob))),(1.0/3.0));


/* Powers of ZETA */
    ZETA_p2=ZETA*ZETA;
    ZETA_p3=ZETA_p2*ZETA;
    ZETA_p4=ZETA_p2*ZETA_p2;

    P1ZETA=1.0 + ZETA;
    M1ZETA=1.0 - ZETA;

    P1ZETA_p1f3=pow(P1ZETA,(1.0/3.0));
    P1ZETA_p2f3=P1ZETA_p1f3*P1ZETA_p1f3;
    P1ZETA_p4f3=P1ZETA_p1f3*P1ZETA;
    P1ZETA_p5f3=P1ZETA_p2f3*P1ZETA;
    M1ZETA_p1f3=pow(M1ZETA,(1.0/3.0));
    M1ZETA_p2f3=M1ZETA_p1f3*M1ZETA_p1f3;
    M1ZETA_p4f3=M1ZETA_p1f3*M1ZETA;
    M1ZETA_p5f3=M1ZETA_p2f3*M1ZETA;

/* Derivatives od RS */

    RS_10=-(1.0/roa_rob_p4f3)*nRS_10;
    RS_01=RS_10;

    RS_20=(1.0/roa_rob_p7f3)*nRS_20;
    RS_02=RS_20;
    RS_11=RS_20;

    RS_30=-(1.0/roa_rob_p10f3)*nRS_30;
    RS_21=RS_30;
    RS_12=RS_30;
    RS_03=RS_30;

/* Powers of RS */
    RS_10_p2=RS_10*RS_10;
    RS_10_p3=RS_10_p2*RS_10;
    RS_01_p2=RS_01*RS_01;
    RS_01_p3=RS_01_p2*RS_01;

/* Derivatives of ZETA */

    ZETA_10=-((roa - rob)/roa_rob_p2) + (1.0/roa_rob);
    ZETA_01=-((roa - rob)/roa_rob_p2) - (1.0/roa_rob);

    ZETA_20=2.0*((roa-rob)/roa_rob_p3) - (2.0/roa_rob_p2);
    ZETA_02=2.0*((roa-rob)/roa_rob_p3) + (2.0/roa_rob_p2);
    ZETA_11=2.0*((roa-rob)/roa_rob_p3);

    ZETA_30=-6.0*((roa-rob)/roa_rob_p4) + (6.0/roa_rob_p3);
    ZETA_03=-6.0*((roa-rob)/roa_rob_p4) - (6.0/roa_rob_p3);
    ZETA_21=-6.0*((roa-rob)/roa_rob_p4) + (2.0/roa_rob_p3);
    ZETA_12=-6.0*((roa-rob)/roa_rob_p4) - (2.0/roa_rob_p3);

/* Powers of ZETA */

    ZETA_10_p2=ZETA_10*ZETA_10;
    ZETA_10_p3=ZETA_10_p2*ZETA_10;
    ZETA_10_p4=ZETA_10_p3*ZETA_10;
    ZETA_01_p2=ZETA_01*ZETA_01;
    ZETA_01_p3=ZETA_01_p2*ZETA_01;


/* Derivatives of F (with respect to ZETA) */
    F=(P1ZETA_p4f3 + M1ZETA_p4f3 - 2.0)/Fdenom;

    F_1ZETA=((4.0/3.0)*P1ZETA_p1f3 - (4.0/3.0)*M1ZETA_p1f3)/Fdenom;

    F_2ZETA=(4.0/(9.0*Fdenom))*( (1.0/M1ZETA_p2f3) + (1.0/P1ZETA_p2f3)   );

    F_3ZETA=(8.0/(27.0*Fdenom))*( (1.0/M1ZETA_p5f3) - (1.0/P1ZETA_p5f3)  );


/* Functions E0, E1, minA */

    E0=Gfit(&RS, &E0_brpa);
    E1=Gfit(&RS, &E1_brpa);
    minA=Gfit(&RS, &minA_brpa);

/* Derivatives of  E0, E1, minA (with respect to RS) */

    E0_1RS=Gfit_1RS(&RS, &E0_brpa);
    E0_2RS=Gfit_2RS(&RS, &E0_brpa);
    E0_3RS=Gfit_3RS(&RS, &E0_brpa);

    E1_1RS=Gfit_1RS(&RS, &E1_brpa);
    E1_2RS=Gfit_2RS(&RS, &E1_brpa);
    E1_3RS=Gfit_3RS(&RS, &E1_brpa);

    minA_1RS=Gfit_1RS(&RS, &minA_brpa);
    minA_2RS=Gfit_2RS(&RS, &minA_brpa);
    minA_3RS=Gfit_3RS(&RS, &minA_brpa);


/* Derivatives */

    ds->df1000 += factor*( 
                           E0_1RS*RS_10 - 
                           iF0*F*(1.0-ZETA_p4)*minA_1RS*RS_10 + 
                           F*ZETA_p4*(-E0_1RS*RS_10 + E1_1RS*RS_10) + 
                           4.0*(-E0 + E1)*F*ZETA_p3*ZETA_10 + 
                           4.0*iF0*F*ZETA_p3*minA*ZETA_10 + 
                           (-E0+E1)*ZETA_p4*F_1ZETA*ZETA_10 - 
                           iF0*(1.0-ZETA_p4)*minA*F_1ZETA*ZETA_10
	);

    ds->df0100 += factor*( 
                           E0_1RS*RS_01 - 
                           iF0*F*(1.0-ZETA_p4)*minA_1RS*RS_01 + 
                           F*ZETA_p4*(-E0_1RS*RS_01 + E1_1RS*RS_01) + 
                           4.0*(-E0 + E1)*F*ZETA_p3*ZETA_01 + 
                           4.0*iF0*F*ZETA_p3*minA*ZETA_01 + 
                           (-E0+E1)*ZETA_p4*F_1ZETA*ZETA_01 - 
                           iF0*(1.0-ZETA_p4)*minA*F_1ZETA*ZETA_01
	);

    ds->df2000 += factor*( 
                           E0_2RS*RS_10_p2 + 
                           iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_10_p2 + 
                           8.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*RS_10*ZETA_10 + 
                           2.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*RS_10*ZETA_10 + 
                           8.0*iF0*F*ZETA_p3*minA_1RS*RS_10*ZETA_10 + 
                           2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_10*ZETA_10 + 
                           12.0*(-E0 + E1)*F*ZETA_p2*ZETA_10_p2 + 
                           12.0*iF0*F*ZETA_p2*minA*ZETA_10_p2 + 
                           8.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_10_p2 + 
                           8.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_10_p2 + 
                           (-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_10_p2 + 
                           iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_10_p2 + 
                           E0_1RS*RS_20 + 
                           iF0*F*(-1.0 + ZETA_p4)*minA_1RS*RS_20 + 
                           F*ZETA_p4*((-E0_2RS + E1_2RS)*RS_10_p2 + 
                                      (-E0_1RS + E1_1RS)*RS_20) + 
                           4.0*(-E0 + E1)*F*ZETA_p3*ZETA_20 + 
                           4.0*iF0*F*ZETA_p3*minA*ZETA_20 + 
                           (-E0 + E1)*ZETA_p4*F_1ZETA*ZETA_20 + 
                           iF0*(-1.0 + ZETA_p4)*minA*F_1ZETA*ZETA_20
	);

    ds->df0200 += factor*( 
                           E0_2RS*RS_01_p2 + 
                           iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_01_p2 + 
                           8.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*RS_01*ZETA_01 + 
                           2.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_01 + 
                           8.0*iF0*F*ZETA_p3*minA_1RS*RS_01*ZETA_01 + 
                           2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_01*ZETA_01 + 
                           12.0*(-E0 + E1)*F*ZETA_p2*ZETA_01_p2 + 
                           12.0*iF0*F*ZETA_p2*minA*ZETA_01_p2 + 
                           8.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_01_p2 + 
                           8.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_01_p2 + 
                           (-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_01_p2 + 
                           iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_01_p2 + 
                           E0_1RS*RS_02 + 
                           iF0*F*(-1.0 + ZETA_p4)*minA_1RS*RS_02 + 
                           F*ZETA_p4*((-E0_2RS + E1_2RS)*RS_01_p2 + 
                                      (-E0_1RS + E1_1RS)*RS_02) + 
                           4.0*(-E0 + E1)*F*ZETA_p3*ZETA_02 + 
                           4.0*iF0*F*ZETA_p3*minA*ZETA_02 + 
                           (-E0 + E1)*ZETA_p4*F_1ZETA*ZETA_02 + 
                           iF0*(-1.0 + ZETA_p4)*minA*F_1ZETA*ZETA_02
	);

    ds->df1100 += factor*( 
                           E0_2RS*RS_01*RS_10+iF0*F*(-1.0+ZETA_p4)*minA_2RS*RS_01*RS_10+
                           4.0*F*ZETA_p3*(-E0_1RS+E1_1RS)*ZETA_01*RS_10+
                           ZETA_p4*(-E0_1RS+E1_1RS)*F_1ZETA*ZETA_01*RS_10+
                           4.0*iF0*F*ZETA_p3*minA_1RS*ZETA_01*RS_10+
                           iF0*(-1.0+ZETA_p4)*F_1ZETA*minA_1RS*ZETA_01*RS_10+
                           4.0*F*ZETA_p3*(-E0_1RS+E1_1RS)*RS_01*ZETA_10+
                           ZETA_p4*(-E0_1RS+E1_1RS)*F_1ZETA*RS_01*ZETA_10+
                           4.0*iF0*F*ZETA_p3*minA_1RS*RS_01*ZETA_10+
                           iF0*(-1.0+ZETA_p4)*F_1ZETA*minA_1RS*RS_01*ZETA_10+
                           12.0*(-E0+E1)*F*ZETA_p2*ZETA_01*ZETA_10+
                           12.0*iF0*F*ZETA_p2*minA*ZETA_01*ZETA_10+
                           8.0*(-E0+E1)*ZETA_p3*F_1ZETA*ZETA_01*ZETA_10+
                           8.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_01*ZETA_10+
                           (-E0+E1)*ZETA_p4*F_2ZETA*ZETA_01*ZETA_10+
                           iF0*(-1.0+ZETA_p4)*minA*F_2ZETA*ZETA_01*ZETA_10+
                           E0_1RS*RS_11+iF0*F*(-1.0+ZETA_p4)*minA_1RS*RS_11+
                           F*ZETA_p4*((-E0_2RS+E1_2RS)*RS_01*RS_10+
                                      (-E0_1RS+E1_1RS)*RS_11)+
                           4.0*(-E0+E1)*F*ZETA_p3*ZETA_11+4.0*iF0*F*ZETA_p3*minA*ZETA_11+
                           (-E0+E1)*ZETA_p4*F_1ZETA*ZETA_11+iF0*(-1.0+ZETA_p4)*minA*F_1ZETA*ZETA_11
	);



    ds->df3000 += factor* ( 
                            E0_3RS*RS_10_p3 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_3RS*RS_10_p3 + 
                            12.0*iF0*F*ZETA_p3*minA_2RS*RS_10_p2*ZETA_10 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_10_p2*ZETA_10 + 
                            36.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*RS_10*ZETA_10_p2 + 
                            24.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*RS_10*ZETA_10_p2 + 
                            36.0*iF0*F*ZETA_p2*minA_1RS*RS_10*ZETA_10_p2  + 
                            24.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*RS_10*ZETA_10_p2 + 
                            3.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*RS_10*ZETA_10_p2 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*RS_10*ZETA_10_p2 + 
                            24.0*(-E0 + E1)*F*ZETA*ZETA_10_p3 + 
                            24.0*iF0*F*ZETA*minA*ZETA_10_p3 + 
                            36.0*(-E0 + E1)*ZETA_p2*F_1ZETA*ZETA_10_p3 + 
                            36.0*iF0*ZETA_p2*minA*F_1ZETA*ZETA_10_p3 + 
                            12.0*(-E0 + E1)*ZETA_p3*F_2ZETA*ZETA_10_p3 + 
                            12.0*iF0*ZETA_p3*minA*F_2ZETA*ZETA_10_p3 + 
                            (-E0 + E1)*ZETA_p4*F_3ZETA*ZETA_10_p3 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_3ZETA*ZETA_10_p3 + 
                            3.0*E0_2RS*RS_10*RS_20 + 
                            3.0*iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_10*RS_20 + 
                            12.0*iF0*F*ZETA_p3*minA_1RS*ZETA_10*RS_20 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*ZETA_10*RS_20 + 
                            12.0*F*ZETA_p3*ZETA_10*((-E0_2RS + E1_2RS)*RS_10_p2 + 
                                                    (-E0_1RS + E1_1RS)*RS_20) + 
                            3.0*ZETA_p4*F_1ZETA*ZETA_10*((-E0_2RS + E1_2RS)*RS_10_p2 + 
                                                         (-E0_1RS + E1_1RS)*RS_20) + 
                            12.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*RS_10*ZETA_20 + 
                            3.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*RS_10*ZETA_20 + 
                            12.0*iF0*F*ZETA_p3*minA_1RS*RS_10*ZETA_20 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_10*ZETA_20 + 
                            36.0*(-E0 + E1)*F*ZETA_p2*ZETA_10*ZETA_20 + 
                            36.0*iF0*F*ZETA_p2*minA*ZETA_10*ZETA_20 + 
                            24.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_10*ZETA_20 + 
                            24.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_10*ZETA_20 + 
                            3.0*(-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_10*ZETA_20 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_10*ZETA_20 + 
                            E0_1RS*RS_30 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_1RS*RS_30 + 
                            F*ZETA_p4*(-((E0_3RS-E1_3RS)*RS_10_p3)-3.0*(E0_2RS-E1_2RS)*RS_10*RS_20 + 
                                       (-E0_1RS + E1_1RS)*RS_30) + 
                            4.0*(-E0 + E1)*F*ZETA_p3*ZETA_30 + 
                            4.0*iF0*F*ZETA_p3*minA*ZETA_30 + 
                            (-E0 + E1)*ZETA_p4*F_1ZETA*ZETA_30 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_1ZETA*ZETA_30
	);


    ds->df0300 += factor* ( 
                            E0_3RS*RS_01_p3 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_3RS*RS_01_p3 + 
                            12.0*iF0*F*ZETA_p3*minA_2RS*RS_01_p2*ZETA_01 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_01_p2*ZETA_01 + 
                            36.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*RS_01*ZETA_01_p2 + 
                            24.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_01_p2 + 
                            36.0*iF0*F*ZETA_p2*minA_1RS*RS_01*ZETA_01_p2  + 
                            24.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*RS_01*ZETA_01_p2 + 
                            3.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*RS_01*ZETA_01_p2 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*RS_01*ZETA_01_p2 + 
                            24.0*(-E0 + E1)*F*ZETA*ZETA_01_p3 + 
                            24.0*iF0*F*ZETA*minA*ZETA_01_p3 + 
                            36.0*(-E0 + E1)*ZETA_p2*F_1ZETA*ZETA_01_p3 + 
                            36.0*iF0*ZETA_p2*minA*F_1ZETA*ZETA_01_p3 + 
                            12.0*(-E0 + E1)*ZETA_p3*F_2ZETA*ZETA_01_p3 + 
                            12.0*iF0*ZETA_p3*minA*F_2ZETA*ZETA_01_p3 + 
                            (-E0 + E1)*ZETA_p4*F_3ZETA*ZETA_01_p3 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_3ZETA*ZETA_01_p3 + 
                            3.0*E0_2RS*RS_01*RS_02 + 
                            3.0*iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_01*RS_02 + 
                            12.0*iF0*F*ZETA_p3*minA_1RS*ZETA_01*RS_02 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*ZETA_01*RS_02 + 
                            12.0*F*ZETA_p3*ZETA_01*((-E0_2RS + E1_2RS)*RS_01_p2 + 
                                                    (-E0_1RS + E1_1RS)*RS_02) + 
                            3.0*ZETA_p4*F_1ZETA*ZETA_01*((-E0_2RS + E1_2RS)*RS_01_p2 + 
                                                         (-E0_1RS + E1_1RS)*RS_02) + 
                            12.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*RS_01*ZETA_02 + 
                            3.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_02 + 
                            12.0*iF0*F*ZETA_p3*minA_1RS*RS_01*ZETA_02 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_01*ZETA_02 + 
                            36.0*(-E0 + E1)*F*ZETA_p2*ZETA_01*ZETA_02 + 
                            36.0*iF0*F*ZETA_p2*minA*ZETA_01*ZETA_02 + 
                            24.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_01*ZETA_02 + 
                            24.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_01*ZETA_02 + 
                            3.0*(-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_01*ZETA_02 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_01*ZETA_02 + 
                            E0_1RS*RS_03 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_1RS*RS_03 + 
                            F*ZETA_p4*(-((E0_3RS-E1_3RS)*RS_01_p3)-3.0*(E0_2RS-E1_2RS)*RS_01*RS_02 + 
                                       (-E0_1RS + E1_1RS)*RS_03) + 
                            4.0*(-E0 + E1)*F*ZETA_p3*ZETA_03 + 
                            4.0*iF0*F*ZETA_p3*minA*ZETA_03 + 
                            (-E0 + E1)*ZETA_p4*F_1ZETA*ZETA_03 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_1ZETA*ZETA_03
	);


    ds->df2100 += factor* (
                           E0_3RS*RS_01*RS_10_p2 + 
                           iF0*F*(-1.0 + ZETA_p4)*minA_3RS*RS_01*RS_10_p2 + 
                           4.0*iF0*F*ZETA_p3*minA_2RS*ZETA_01*RS_10_p2 + 
                           iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*ZETA_01*RS_10_p2 + 
                           8.0*iF0*F*ZETA_p3*minA_2RS*RS_01*RS_10*ZETA_10 + 
                           2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_01*RS_10*ZETA_10 + 
                           24.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*ZETA_01*RS_10*ZETA_10 + 
                           16.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*ZETA_01*RS_10*ZETA_10 + 
                           24.0*iF0*F*ZETA_p2*minA_1RS*ZETA_01*RS_10*ZETA_10 + 
                           16.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*ZETA_01*RS_10*ZETA_10 + 
                           2.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*ZETA_01*RS_10*ZETA_10 + 
                           2.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*ZETA_01*RS_10*ZETA_10 + 
                           12.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*RS_01*ZETA_10_p2 + 
                           8.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_10_p2 + 
                           12.0*iF0*F*ZETA_p2*minA_1RS*RS_01*ZETA_10_p2 + 
                           8.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*RS_01*ZETA_10_p2 + 
                           ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*RS_01*ZETA_10_p2 + 
                           iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*RS_01*ZETA_10_p2 + 
                           24.0*(-E0 + E1)*F*ZETA*ZETA_01*ZETA_10_p2 + 
                           24.0*iF0*F*ZETA*minA*ZETA_01*ZETA_10_p2 + 
                           36.0*(-E0 + E1)*ZETA_p2*F_1ZETA*ZETA_01*ZETA_10_p2 + 
                           36.0*iF0*ZETA_p2*minA*F_1ZETA*ZETA_01*ZETA_10_p2 + 
                           12.0*(-E0 + E1)*ZETA_p3*F_2ZETA*ZETA_01*ZETA_10_p2 + 
                           12.0*iF0*ZETA_p3*minA*F_2ZETA*ZETA_01*ZETA_10_p2 + 
                           (-E0 + E1)*ZETA_p4*F_3ZETA*ZETA_01*ZETA_10_p2 + 
                           iF0*(-1.0 + ZETA_p4)*minA*F_3ZETA*ZETA_01*ZETA_10_p2 + 
                           2.0*E0_2RS*RS_10*RS_11 + 
                           2.0*iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_10*RS_11 + 
                           8.0*iF0*F*ZETA_p3*minA_1RS*ZETA_10*RS_11 + 
                           2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*ZETA_10*RS_11 + 
                           8.0*F*ZETA_p3*ZETA_10*((-E0_2RS + 
                                                   E1_2RS)*RS_01*RS_10 + (-E0_1RS + E1_1RS)*RS_11) + 
                           2.0*ZETA_p4*F_1ZETA*ZETA_10*((-E0_2RS + E1_2RS)*RS_01*RS_10 + 
                                                        (-E0_1RS + E1_1RS)*RS_11) + 
                           8.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*RS_10*ZETA_11 + 
                           2.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*RS_10*ZETA_11 + 
                           8.0*iF0*F*ZETA_p3*minA_1RS*RS_10*ZETA_11 + 
                           2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_10*ZETA_11 + 
                           24.0*(-E0 + E1)*F*ZETA_p2*ZETA_10*ZETA_11 + 
                           24.0*iF0*F*ZETA_p2*minA*ZETA_10*ZETA_11 + 
                           16.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_10*ZETA_11 + 
                           16.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_10*ZETA_11 + 
                           2.0*(-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_10*ZETA_11 + 
                           2.0*iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_10*ZETA_11 + 
                           E0_2RS*RS_01*RS_20 + 
                           iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_01*RS_20 + 
                           4.0*iF0*F*ZETA_p3*minA_1RS*ZETA_01*RS_20 + 
                           iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*ZETA_01*RS_20 + 
                           4.0*F*ZETA_p3*ZETA_01*((-E0_2RS + E1_2RS)*RS_10_p2 + 
                                                  (-E0_1RS + E1_1RS)*RS_20) + 
                           ZETA_p4*F_1ZETA*ZETA_01*((-E0_2RS + E1_2RS)*RS_10_p2 + 
                                                    (-E0_1RS + E1_1RS)*RS_20) + 
                           4.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*RS_01*ZETA_20 + 
                           ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_20 + 
                           4.0*iF0*F*ZETA_p3*minA_1RS*RS_01*ZETA_20 + 
                           iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_01*ZETA_20 + 
                           12.0*(-E0 + E1)*F*ZETA_p2*ZETA_01*ZETA_20 + 
                           12.0*iF0*F*ZETA_p2*minA*ZETA_01*ZETA_20 + 
                           8.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_01*ZETA_20 + 
                           8.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_01*ZETA_20 + 
                           (-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_01*ZETA_20 + 
                           iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_01*ZETA_20 + 
                           E0_1RS*RS_21 + 
                           iF0*F*(-1.0 + ZETA_p4)*minA_1RS*RS_21 + 
                           F*ZETA_p4*(-(E0_3RS*RS_01*RS_10_p2) + 
                                      E1_3RS*RS_01*RS_10_p2-(E0_2RS-E1_2RS)*(2.0*RS_10*RS_11 + 
                                                                             RS_01*RS_20) + 
                                      (-E0_1RS + E1_1RS)*RS_21) + 
                           4.0*(-E0 + E1)*F*ZETA_p3*ZETA_21 + 
                           4.0*iF0*F*ZETA_p3*minA*ZETA_21 + 
                           (-E0 + E1)*ZETA_p4*F_1ZETA*ZETA_21 + 
                           iF0*(-1.0 + ZETA_p4)*minA*F_1ZETA*ZETA_21
	);



    ds->df1200 += factor* ( 
                            E0_3RS*RS_01_p2*RS_10 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_3RS*RS_01_p2*RS_10 + 
                            8.0*iF0*F*ZETA_p3*minA_2RS*RS_01*ZETA_01*RS_10 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_01*ZETA_01*RS_10 + 
                            12.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*ZETA_01_p2*RS_10 + 
                            8.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*ZETA_01_p2*RS_10 + 
                            12.0*iF0*F*ZETA_p2*minA_1RS*ZETA_01_p2*RS_10 + 
                            8.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*ZETA_01_p2*RS_10 + 
                            ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*ZETA_01_p2*RS_10 + 
                            iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*ZETA_01_p2*RS_10 + 
                            E0_2RS*RS_02*RS_10 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_02*RS_10 + 
                            4.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*ZETA_02*RS_10 + 
                            ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*ZETA_02*RS_10 + 
                            4.0*iF0*F*ZETA_p3*minA_1RS*ZETA_02*RS_10 + 
                            iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*ZETA_02*RS_10 + 
                            4.0*iF0*F*ZETA_p3*minA_2RS*RS_01_p2*ZETA_10 + 
                            iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_01_p2*ZETA_10 + 
                            24.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*RS_01*ZETA_01*ZETA_10 + 
                            16.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_01*ZETA_10 + 
                            24.0*iF0*F*ZETA_p2*minA_1RS*RS_01*ZETA_01*ZETA_10 + 
                            16.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*RS_01*ZETA_01*ZETA_10 + 
                            2.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*RS_01*ZETA_01*ZETA_10 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*RS_01*ZETA_01*ZETA_10 + 
                            24.0*(-E0 + E1)*F*ZETA*ZETA_01_p2*ZETA_10 + 
                            24.0*iF0*F*ZETA*minA*ZETA_01_p2*ZETA_10 + 
                            36.0*(-E0 + E1)*ZETA_p2*F_1ZETA*ZETA_01_p2*ZETA_10 + 
                            36.0*iF0*ZETA_p2*minA*F_1ZETA*ZETA_01_p2*ZETA_10 + 
                            12.0*(-E0 + E1)*ZETA_p3*F_2ZETA*ZETA_01_p2*ZETA_10 + 
                            12.0*iF0*ZETA_p3*minA*F_2ZETA*ZETA_01_p2*ZETA_10 + 
                            (-E0 + E1)*ZETA_p4*F_3ZETA*ZETA_01_p2*ZETA_10 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_3ZETA*ZETA_01_p2*ZETA_10 + 
                            4.0*iF0*F*ZETA_p3*minA_1RS*RS_02*ZETA_10 + 
                            iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_02*ZETA_10 + 
                            4.0*F*ZETA_p3*((-E0_2RS + E1_2RS)*RS_01_p2 + 
                                           (-E0_1RS + E1_1RS)*RS_02)*ZETA_10 + 
                            ZETA_p4*F_1ZETA*((-E0_2RS + E1_2RS)*RS_01_p2 + 
                                             (-E0_1RS + E1_1RS)*RS_02)*ZETA_10 + 
                            12.0*(-E0 + E1)*F*ZETA_p2*ZETA_02*ZETA_10 + 
                            12.0*iF0*F*ZETA_p2*minA*ZETA_02*ZETA_10 + 
                            8.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_02*ZETA_10 + 
                            8.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_02*ZETA_10 + 
                            (-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_02*ZETA_10 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_02*ZETA_10 + 
                            2.0*E0_2RS*RS_01*RS_11 + 
                            2.0*iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_01*RS_11 + 
                            8.0*iF0*F*ZETA_p3*minA_1RS*ZETA_01*RS_11 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*ZETA_01*RS_11 + 
                            8.0*F*ZETA_p3*ZETA_01*((-E0_2RS + E1_2RS)*RS_01*RS_10 + 
                                                   (-E0_1RS + E1_1RS)*RS_11) + 
                            2.0*ZETA_p4*F_1ZETA*ZETA_01*((-E0_2RS + E1_2RS)*RS_01*RS_10 + 
                                                         (-E0_1RS + E1_1RS)*RS_11) + 
                            8.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*RS_01*ZETA_11 + 
                            2.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_11 + 
                            8.0*iF0*F*ZETA_p3*minA_1RS*RS_01*ZETA_11 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_01*ZETA_11 + 
                            24.0*(-E0 + E1)*F*ZETA_p2*ZETA_01*ZETA_11 + 
                            24.0*iF0*F*ZETA_p2*minA*ZETA_01*ZETA_11 + 16.0*(-E0 + 
                                                                            E1)*ZETA_p3*F_1ZETA*ZETA_01*ZETA_11 + 
                            16.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_01*ZETA_11 + 
                            2.0*(-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_01*ZETA_11 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_01*ZETA_11 + 
                            E0_1RS*RS_12 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_1RS*RS_12 + 
                            F*ZETA_p4*(-(((E0_3RS-E1_3RS)*RS_01_p2 + 
                                          (E0_2RS-E1_2RS)*RS_02)*RS_10)-2.0*(E0_2RS-E1_2RS)*RS_01*RS_11 + 
                                       (-E0_1RS + E1_1RS)*RS_12) + 
                            4.0*(-E0 + E1)*F*ZETA_p3*ZETA_12 + 
                            4.0*iF0*F*ZETA_p3*minA*ZETA_12 + 
                            (-E0 + E1)*ZETA_p4*F_1ZETA*ZETA_12 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_1ZETA*ZETA_12
	);

}
  
#ifdef FOURTH_ORDER_DERIVATIVES
static void
pw92_fourth(FourthFuncDrv *ds, real factor, const DftDensProp* dp) {
/* Declarations */
    real  E0;
    real  E0_1RS;
    real  E0_2RS;
    real  E0_3RS;
    real  E0_4RS;
    real  E1;
    real  E1_1RS;
    real  E1_2RS;
    real  E1_3RS;
    real  E1_4RS;
    real  F;
    real  F_1ZETA;
    real  F_2ZETA;
    real  F_3ZETA;
    real  F_4ZETA;
    real  M1ZETA;
    real  M1ZETA_p1f3;
    real  M1ZETA_p2f3;
    real  M1ZETA_p4f3;
    real  M1ZETA_p5f3;
    real  M1ZETA_p8f3;
    real  P1ZETA;
    real  P1ZETA_p1f3;
    real  P1ZETA_p2f3;
    real  P1ZETA_p4f3;
    real  P1ZETA_p5f3;
    real  P1ZETA_p8f3;
    real  RS;
    real  RS_01;
    real  RS_01_p2;
    real  RS_01_p3;
    real  RS_01_p4;
    real  RS_02;
    real  RS_02_p2;
    real  RS_03;
    real  RS_04;
    real  RS_10;
    real  RS_10_p2;
    real  RS_10_p3;
    real  RS_10_p4;
    real  RS_11;
    real  RS_11_p2;
    real  RS_12;
    real  RS_13;
    real  RS_20;
    real  RS_20_p2;
    real  RS_21;
    real  RS_22;
    real  RS_30;
    real  RS_31;
    real  RS_40;
    real  ZETA;
    real  ZETA_01;
    real  ZETA_01_p2;
    real  ZETA_01_p3;
    real  ZETA_01_p4;
    real  ZETA_02;
    real  ZETA_02_p2;
    real  ZETA_03;
    real  ZETA_04;
    real  ZETA_10;
    real  ZETA_10_p2;
    real  ZETA_10_p3;
    real  ZETA_10_p4;
    real  ZETA_11;
    real  ZETA_11_p2;
    real  ZETA_12;
    real  ZETA_13;
    real  ZETA_20;
    real  ZETA_20_p2;
    real  ZETA_21;
    real  ZETA_22;
    real  ZETA_30;
    real  ZETA_31;
    real  ZETA_40;
    real  ZETA_p2;
    real  ZETA_p3;
    real  ZETA_p4;
    real  minA;
    real  minA_1RS;
    real  minA_2RS;
    real  minA_3RS;
    real  minA_4RS;
    real  roa;
    real  roa_rob;
    real  roa_rob_p10f3;
    real  roa_rob_p13f3;
    real  roa_rob_p1f3;
    real  roa_rob_p2;
    real  roa_rob_p3;
    real  roa_rob_p4;
    real  roa_rob_p4f3;
    real  roa_rob_p5;
    real  roa_rob_p7f3;
    real  rob;

/* Numbers and constants */

    const real  Fdenom=(pow(2.0,(4.0/3.0)) - 2.0);
    const real  iF0=1.0/1.709921;

/* Numbers in expressions */
    const real  nRS_10=1.0/(pow(6.0,(2.0/3.0))*pow(M_PI,(1.0/3.0)));
    const real  nRS_20=2.0*(pow((2.0/M_PI),(1.0/3.0)))/(3.0*pow(3.0,(2.0/3.0)));
    const real  nRS_30=14.0*(pow((2.0/M_PI),(1.0/3.0)))/(9.0*pow(3.0,(2.0/3.0)));
    const real  nRS_40=(140.0/27.0)*(1.0/(pow(3.0,(2.0/3.0))))*(pow((2.0/M_PI),(1.0/3.0)));

/* Setting up alpha and beta densities roa and rob */
    roa=dp->rhoa;
    rob=dp->rhob;

/* Ignore small densities */
   if((roa+rob)<PW92THR) return;

/* (roa + rob) and its powers */
    roa_rob=roa+rob;
    roa_rob_p2=roa_rob*roa_rob;
    roa_rob_p3=roa_rob_p2*roa_rob;
    roa_rob_p4=roa_rob_p2*roa_rob_p2;
    roa_rob_p5=roa_rob_p4*roa_rob;
    roa_rob_p1f3=pow(roa_rob,(1.0/3.0));
    roa_rob_p4f3=roa_rob*roa_rob_p1f3;
    roa_rob_p7f3=roa_rob_p4f3*roa_rob;
    roa_rob_p10f3=roa_rob_p7f3*roa_rob;
    roa_rob_p13f3=roa_rob_p10f3*roa_rob;

/* Auxaliary functions RS and ZETA */

    ZETA=(roa-rob)/(roa+rob);
    RS=pow((3.0/(4.0*M_PI*(roa+rob))),(1.0/3.0));


/* Powers of ZETA */
    ZETA_p2=ZETA*ZETA;
    ZETA_p3=ZETA_p2*ZETA;
    ZETA_p4=ZETA_p2*ZETA_p2;

    P1ZETA=1.0 + ZETA;
    M1ZETA=1.0 - ZETA;

    P1ZETA_p1f3=pow(P1ZETA,(1.0/3.0));
    P1ZETA_p2f3=P1ZETA_p1f3*P1ZETA_p1f3;
    P1ZETA_p4f3=P1ZETA_p1f3*P1ZETA;
    P1ZETA_p5f3=P1ZETA_p2f3*P1ZETA;
    P1ZETA_p8f3=P1ZETA_p5f3*P1ZETA;
    M1ZETA_p1f3=pow(M1ZETA,(1.0/3.0));
    M1ZETA_p2f3=M1ZETA_p1f3*M1ZETA_p1f3;
    M1ZETA_p4f3=M1ZETA_p1f3*M1ZETA;
    M1ZETA_p5f3=M1ZETA_p2f3*M1ZETA;
    M1ZETA_p8f3=M1ZETA_p5f3*M1ZETA;

/* Derivatives od RS */

    RS_10=-(1.0/roa_rob_p4f3)*nRS_10;
    RS_01=RS_10;

    RS_20=(1.0/roa_rob_p7f3)*nRS_20;
    RS_02=RS_20;
    RS_11=RS_20;

    RS_30=-(1.0/roa_rob_p10f3)*nRS_30;
    RS_21=RS_30;
    RS_12=RS_30;
    RS_03=RS_30;

    RS_40=(1.0/roa_rob_p13f3)*nRS_40;
    RS_31=RS_40;
    RS_22=RS_40;
    RS_13=RS_40;
    RS_04=RS_40;

/* Powers of RS */
    RS_10_p2=RS_10*RS_10;
    RS_10_p3=RS_10_p2*RS_10;
    RS_10_p4=RS_10_p3*RS_10;
    RS_01_p2=RS_01*RS_01;
    RS_01_p3=RS_01_p2*RS_01;
    RS_01_p4=RS_01_p3*RS_01;

    RS_20_p2=RS_20*RS_20;
    RS_02_p2=RS_02*RS_02;
    RS_11_p2=RS_11*RS_11;

/* Derivatives of ZETA */

    ZETA_10=-((roa - rob)/roa_rob_p2) + (1.0/roa_rob);
    ZETA_01=-((roa - rob)/roa_rob_p2) - (1.0/roa_rob);

    ZETA_20=2.0*((roa-rob)/roa_rob_p3) - (2.0/roa_rob_p2);
    ZETA_02=2.0*((roa-rob)/roa_rob_p3) + (2.0/roa_rob_p2);
    ZETA_11=2.0*((roa-rob)/roa_rob_p3);

    ZETA_30=-6.0*((roa-rob)/roa_rob_p4) + (6.0/roa_rob_p3);
    ZETA_03=-6.0*((roa-rob)/roa_rob_p4) - (6.0/roa_rob_p3);
    ZETA_21=-6.0*((roa-rob)/roa_rob_p4) + (2.0/roa_rob_p3);
    ZETA_12=-6.0*((roa-rob)/roa_rob_p4) - (2.0/roa_rob_p3);

    ZETA_40=24.0*((roa-rob)/roa_rob_p5) - (24.0/roa_rob_p4);
    ZETA_31=24.0*((roa-rob)/roa_rob_p5) - (12.0/roa_rob_p4);
    ZETA_22=24.0*((roa-rob)/roa_rob_p5);
    ZETA_13=24.0*((roa-rob)/roa_rob_p5) + (12.0/roa_rob_p4);
    ZETA_04=24.0*((roa-rob)/roa_rob_p5) + (24.0/roa_rob_p4);


/* Powers of ZETA */

    ZETA_10_p2=ZETA_10*ZETA_10;
    ZETA_10_p3=ZETA_10_p2*ZETA_10;
    ZETA_10_p4=ZETA_10_p3*ZETA_10;
    ZETA_01_p2=ZETA_01*ZETA_01;
    ZETA_01_p3=ZETA_01_p2*ZETA_01;
    ZETA_01_p4=ZETA_01_p3*ZETA_01;

    ZETA_20_p2=ZETA_20*ZETA_20;
    ZETA_02_p2=ZETA_02*ZETA_02;
    ZETA_11_p2=ZETA_11*ZETA_11;

/* Derivatives of F (with respect to ZETA) */
    F=(P1ZETA_p4f3 + M1ZETA_p4f3 - 2.0)/Fdenom;

    F_1ZETA=((4.0/3.0)*P1ZETA_p1f3 - (4.0/3.0)*M1ZETA_p1f3)/Fdenom;

    F_2ZETA=(4.0/(9.0*Fdenom))*( (1.0/M1ZETA_p2f3) + (1.0/P1ZETA_p2f3)   );

    F_3ZETA=(8.0/(27.0*Fdenom))*( (1.0/M1ZETA_p5f3) - (1.0/P1ZETA_p5f3)  );

    F_4ZETA=(40.0/(81.0*Fdenom))*( (1.0/M1ZETA_p8f3) + (1.0/P1ZETA_p8f3)  );


/* Functions E0, E1, minA */

    E0=Gfit(&RS, &E0_brpa);
    E1=Gfit(&RS, &E1_brpa);
    minA=Gfit(&RS, &minA_brpa);

/* Derivatives of  E0, E1, minA (with respect to RS) */

    E0_1RS=Gfit_1RS(&RS, &E0_brpa);
    E0_2RS=Gfit_2RS(&RS, &E0_brpa);
    E0_3RS=Gfit_3RS(&RS, &E0_brpa);
    E0_4RS=Gfit_4RS(&RS, &E0_brpa);

    E1_1RS=Gfit_1RS(&RS, &E1_brpa);
    E1_2RS=Gfit_2RS(&RS, &E1_brpa);
    E1_3RS=Gfit_3RS(&RS, &E1_brpa);
    E1_4RS=Gfit_4RS(&RS, &E1_brpa);

    minA_1RS=Gfit_1RS(&RS, &minA_brpa);
    minA_2RS=Gfit_2RS(&RS, &minA_brpa);
    minA_3RS=Gfit_3RS(&RS, &minA_brpa);
    minA_4RS=Gfit_4RS(&RS, &minA_brpa);


/* Derivatives */

    ds->df10000 += factor*( 
                            E0_1RS*RS_10 - 
                            iF0*F*(1.0-ZETA_p4)*minA_1RS*RS_10 + 
                            F*ZETA_p4*(-E0_1RS*RS_10 + E1_1RS*RS_10) + 
                            4.0*(-E0 + E1)*F*ZETA_p3*ZETA_10 + 
                            4.0*iF0*F*ZETA_p3*minA*ZETA_10 + 
                            (-E0+E1)*ZETA_p4*F_1ZETA*ZETA_10 - 
                            iF0*(1.0-ZETA_p4)*minA*F_1ZETA*ZETA_10
	);

    ds->df01000 += factor*( 
                            E0_1RS*RS_01 - 
                            iF0*F*(1.0-ZETA_p4)*minA_1RS*RS_01 + 
                            F*ZETA_p4*(-E0_1RS*RS_01 + E1_1RS*RS_01) + 
                            4.0*(-E0 + E1)*F*ZETA_p3*ZETA_01 + 
                            4.0*iF0*F*ZETA_p3*minA*ZETA_01 + 
                            (-E0+E1)*ZETA_p4*F_1ZETA*ZETA_01 - 
                            iF0*(1.0-ZETA_p4)*minA*F_1ZETA*ZETA_01
	);

    ds->df20000 += factor*( 
                            E0_2RS*RS_10_p2 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_10_p2 + 
                            8.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*RS_10*ZETA_10 + 
                            2.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*RS_10*ZETA_10 + 
                            8.0*iF0*F*ZETA_p3*minA_1RS*RS_10*ZETA_10 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_10*ZETA_10 + 
                            12.0*(-E0 + E1)*F*ZETA_p2*ZETA_10_p2 + 
                            12.0*iF0*F*ZETA_p2*minA*ZETA_10_p2 + 
                            8.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_10_p2 + 
                            8.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_10_p2 + 
                            (-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_10_p2 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_10_p2 + 
                            E0_1RS*RS_20 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_1RS*RS_20 + 
                            F*ZETA_p4*((-E0_2RS + E1_2RS)*RS_10_p2 + 
                                       (-E0_1RS + E1_1RS)*RS_20) + 
                            4.0*(-E0 + E1)*F*ZETA_p3*ZETA_20 + 
                            4.0*iF0*F*ZETA_p3*minA*ZETA_20 + 
                            (-E0 + E1)*ZETA_p4*F_1ZETA*ZETA_20 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_1ZETA*ZETA_20
	);

    ds->df02000 += factor*( 
                            E0_2RS*RS_01_p2 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_01_p2 + 
                            8.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*RS_01*ZETA_01 + 
                            2.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_01 + 
                            8.0*iF0*F*ZETA_p3*minA_1RS*RS_01*ZETA_01 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_01*ZETA_01 + 
                            12.0*(-E0 + E1)*F*ZETA_p2*ZETA_01_p2 + 
                            12.0*iF0*F*ZETA_p2*minA*ZETA_01_p2 + 
                            8.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_01_p2 + 
                            8.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_01_p2 + 
                            (-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_01_p2 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_01_p2 + 
                            E0_1RS*RS_02 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_1RS*RS_02 + 
                            F*ZETA_p4*((-E0_2RS + E1_2RS)*RS_01_p2 + 
                                       (-E0_1RS + E1_1RS)*RS_02) + 
                            4.0*(-E0 + E1)*F*ZETA_p3*ZETA_02 + 
                            4.0*iF0*F*ZETA_p3*minA*ZETA_02 + 
                            (-E0 + E1)*ZETA_p4*F_1ZETA*ZETA_02 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_1ZETA*ZETA_02
	);

    ds->df11000 += factor*( 
                            E0_2RS*RS_01*RS_10+iF0*F*(-1.0+ZETA_p4)*minA_2RS*RS_01*RS_10+4.0*F*ZETA_p3*(-E0_1RS+E1_1RS)*ZETA_01*RS_10+ZETA_p4*(-E0_1RS+E1_1RS)*F_1ZETA*ZETA_01*RS_10+4.0*iF0*F*ZETA_p3*minA_1RS*ZETA_01*RS_10+iF0*(-1.0+ZETA_p4)*F_1ZETA*minA_1RS*ZETA_01*RS_10+4.0*F*ZETA_p3*(-E0_1RS+E1_1RS)*RS_01*ZETA_10+ZETA_p4*(-E0_1RS+E1_1RS)*F_1ZETA*RS_01*ZETA_10+4.0*iF0*F*ZETA_p3*minA_1RS*RS_01*ZETA_10+iF0*(-1.0+ZETA_p4)*F_1ZETA*minA_1RS*RS_01*ZETA_10+12.0*(-E0+E1)*F*ZETA_p2*ZETA_01*ZETA_10+12.0*iF0*F*ZETA_p2*minA*ZETA_01*ZETA_10+8.0*(-E0+E1)*ZETA_p3*F_1ZETA*ZETA_01*ZETA_10+8.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_01*ZETA_10+(-E0+E1)*ZETA_p4*F_2ZETA*ZETA_01*ZETA_10+iF0*(-1.0+ZETA_p4)*minA*F_2ZETA*ZETA_01*ZETA_10+E0_1RS*RS_11+iF0*F*(-1.0+ZETA_p4)*minA_1RS*RS_11+F*ZETA_p4*((-E0_2RS+E1_2RS)*RS_01*RS_10+(-E0_1RS+E1_1RS)*RS_11)+4.0*(-E0+E1)*F*ZETA_p3*ZETA_11+4.0*iF0*F*ZETA_p3*minA*ZETA_11+(-E0+E1)*ZETA_p4*F_1ZETA*ZETA_11+iF0*(-1.0+ZETA_p4)*minA*F_1ZETA*ZETA_11
	);



    ds->df30000 += factor* ( 
                             E0_3RS*RS_10_p3 + 
                             iF0*F*(-1.0 + ZETA_p4)*minA_3RS*RS_10_p3 + 
                             12.0*iF0*F*ZETA_p3*minA_2RS*RS_10_p2*ZETA_10 + 
                             3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_10_p2*ZETA_10 + 
                             36.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*RS_10*ZETA_10_p2 + 
                             24.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*RS_10*ZETA_10_p2 + 
                             36.0*iF0*F*ZETA_p2*minA_1RS*RS_10*ZETA_10_p2  + 
                             24.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*RS_10*ZETA_10_p2 + 
                             3.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*RS_10*ZETA_10_p2 + 
                             3.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*RS_10*ZETA_10_p2 + 
                             24.0*(-E0 + E1)*F*ZETA*ZETA_10_p3 + 
                             24.0*iF0*F*ZETA*minA*ZETA_10_p3 + 
                             36.0*(-E0 + E1)*ZETA_p2*F_1ZETA*ZETA_10_p3 + 
                             36.0*iF0*ZETA_p2*minA*F_1ZETA*ZETA_10_p3 + 
                             12.0*(-E0 + E1)*ZETA_p3*F_2ZETA*ZETA_10_p3 + 
                             12.0*iF0*ZETA_p3*minA*F_2ZETA*ZETA_10_p3 + 
                             (-E0 + E1)*ZETA_p4*F_3ZETA*ZETA_10_p3 + 
                             iF0*(-1.0 + ZETA_p4)*minA*F_3ZETA*ZETA_10_p3 + 
                             3.0*E0_2RS*RS_10*RS_20 + 
                             3.0*iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_10*RS_20 + 
                             12.0*iF0*F*ZETA_p3*minA_1RS*ZETA_10*RS_20 + 
                             3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*ZETA_10*RS_20 + 
                             12.0*F*ZETA_p3*ZETA_10*((-E0_2RS + E1_2RS)*RS_10_p2 + 
                                                     (-E0_1RS + E1_1RS)*RS_20) + 
                             3.0*ZETA_p4*F_1ZETA*ZETA_10*((-E0_2RS + E1_2RS)*RS_10_p2 + 
                                                          (-E0_1RS + E1_1RS)*RS_20) + 
                             12.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*RS_10*ZETA_20 + 
                             3.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*RS_10*ZETA_20 + 
                             12.0*iF0*F*ZETA_p3*minA_1RS*RS_10*ZETA_20 + 
                             3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_10*ZETA_20 + 
                             36.0*(-E0 + E1)*F*ZETA_p2*ZETA_10*ZETA_20 + 
                             36.0*iF0*F*ZETA_p2*minA*ZETA_10*ZETA_20 + 
                             24.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_10*ZETA_20 + 
                             24.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_10*ZETA_20 + 
                             3.0*(-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_10*ZETA_20 + 
                             3.0*iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_10*ZETA_20 + 
                             E0_1RS*RS_30 + 
                             iF0*F*(-1.0 + ZETA_p4)*minA_1RS*RS_30 + 
                             F*ZETA_p4*(-((E0_3RS-E1_3RS)*RS_10_p3)-3.0*(E0_2RS-E1_2RS)*RS_10*RS_20 + 
                                        (-E0_1RS + E1_1RS)*RS_30) + 
                             4.0*(-E0 + E1)*F*ZETA_p3*ZETA_30 + 
                             4.0*iF0*F*ZETA_p3*minA*ZETA_30 + 
                             (-E0 + E1)*ZETA_p4*F_1ZETA*ZETA_30 + 
                             iF0*(-1.0 + ZETA_p4)*minA*F_1ZETA*ZETA_30
	);


    ds->df03000 += factor* ( 
                             E0_3RS*RS_01_p3 + 
                             iF0*F*(-1.0 + ZETA_p4)*minA_3RS*RS_01_p3 + 
                             12.0*iF0*F*ZETA_p3*minA_2RS*RS_01_p2*ZETA_01 + 
                             3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_01_p2*ZETA_01 + 
                             36.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*RS_01*ZETA_01_p2 + 
                             24.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_01_p2 + 
                             36.0*iF0*F*ZETA_p2*minA_1RS*RS_01*ZETA_01_p2  + 
                             24.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*RS_01*ZETA_01_p2 + 
                             3.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*RS_01*ZETA_01_p2 + 
                             3.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*RS_01*ZETA_01_p2 + 
                             24.0*(-E0 + E1)*F*ZETA*ZETA_01_p3 + 
                             24.0*iF0*F*ZETA*minA*ZETA_01_p3 + 
                             36.0*(-E0 + E1)*ZETA_p2*F_1ZETA*ZETA_01_p3 + 
                             36.0*iF0*ZETA_p2*minA*F_1ZETA*ZETA_01_p3 + 
                             12.0*(-E0 + E1)*ZETA_p3*F_2ZETA*ZETA_01_p3 + 
                             12.0*iF0*ZETA_p3*minA*F_2ZETA*ZETA_01_p3 + 
                             (-E0 + E1)*ZETA_p4*F_3ZETA*ZETA_01_p3 + 
                             iF0*(-1.0 + ZETA_p4)*minA*F_3ZETA*ZETA_01_p3 + 
                             3.0*E0_2RS*RS_01*RS_02 + 
                             3.0*iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_01*RS_02 + 
                             12.0*iF0*F*ZETA_p3*minA_1RS*ZETA_01*RS_02 + 
                             3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*ZETA_01*RS_02 + 
                             12.0*F*ZETA_p3*ZETA_01*((-E0_2RS + E1_2RS)*RS_01_p2 + 
                                                     (-E0_1RS + E1_1RS)*RS_02) + 
                             3.0*ZETA_p4*F_1ZETA*ZETA_01*((-E0_2RS + E1_2RS)*RS_01_p2 + 
                                                          (-E0_1RS + E1_1RS)*RS_02) + 
                             12.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*RS_01*ZETA_02 + 
                             3.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_02 + 
                             12.0*iF0*F*ZETA_p3*minA_1RS*RS_01*ZETA_02 + 
                             3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_01*ZETA_02 + 
                             36.0*(-E0 + E1)*F*ZETA_p2*ZETA_01*ZETA_02 + 
                             36.0*iF0*F*ZETA_p2*minA*ZETA_01*ZETA_02 + 
                             24.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_01*ZETA_02 + 
                             24.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_01*ZETA_02 + 
                             3.0*(-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_01*ZETA_02 + 
                             3.0*iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_01*ZETA_02 + 
                             E0_1RS*RS_03 + 
                             iF0*F*(-1.0 + ZETA_p4)*minA_1RS*RS_03 + 
                             F*ZETA_p4*(-((E0_3RS-E1_3RS)*RS_01_p3)-3.0*(E0_2RS-E1_2RS)*RS_01*RS_02 + 
                                        (-E0_1RS + E1_1RS)*RS_03) + 
                             4.0*(-E0 + E1)*F*ZETA_p3*ZETA_03 + 
                             4.0*iF0*F*ZETA_p3*minA*ZETA_03 + 
                             (-E0 + E1)*ZETA_p4*F_1ZETA*ZETA_03 + 
                             iF0*(-1.0 + ZETA_p4)*minA*F_1ZETA*ZETA_03
	);


    ds->df21000 += factor* (
                            E0_3RS*RS_01*RS_10_p2 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_3RS*RS_01*RS_10_p2 + 
                            4.0*iF0*F*ZETA_p3*minA_2RS*ZETA_01*RS_10_p2 + 
                            iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*ZETA_01*RS_10_p2 + 
                            8.0*iF0*F*ZETA_p3*minA_2RS*RS_01*RS_10*ZETA_10 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_01*RS_10*ZETA_10 + 
                            24.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*ZETA_01*RS_10*ZETA_10 + 
                            16.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*ZETA_01*RS_10*ZETA_10 + 
                            24.0*iF0*F*ZETA_p2*minA_1RS*ZETA_01*RS_10*ZETA_10 + 
                            16.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*ZETA_01*RS_10*ZETA_10 + 
                            2.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*ZETA_01*RS_10*ZETA_10 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*ZETA_01*RS_10*ZETA_10 + 
                            12.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*RS_01*ZETA_10_p2 + 
                            8.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_10_p2 + 
                            12.0*iF0*F*ZETA_p2*minA_1RS*RS_01*ZETA_10_p2 + 
                            8.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*RS_01*ZETA_10_p2 + 
                            ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*RS_01*ZETA_10_p2 + 
                            iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*RS_01*ZETA_10_p2 + 
                            24.0*(-E0 + E1)*F*ZETA*ZETA_01*ZETA_10_p2 + 
                            24.0*iF0*F*ZETA*minA*ZETA_01*ZETA_10_p2 + 
                            36.0*(-E0 + E1)*ZETA_p2*F_1ZETA*ZETA_01*ZETA_10_p2 + 
                            36.0*iF0*ZETA_p2*minA*F_1ZETA*ZETA_01*ZETA_10_p2 + 
                            12.0*(-E0 + E1)*ZETA_p3*F_2ZETA*ZETA_01*ZETA_10_p2 + 
                            12.0*iF0*ZETA_p3*minA*F_2ZETA*ZETA_01*ZETA_10_p2 + 
                            (-E0 + E1)*ZETA_p4*F_3ZETA*ZETA_01*ZETA_10_p2 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_3ZETA*ZETA_01*ZETA_10_p2 + 
                            2.0*E0_2RS*RS_10*RS_11 + 
                            2.0*iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_10*RS_11 + 
                            8.0*iF0*F*ZETA_p3*minA_1RS*ZETA_10*RS_11 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*ZETA_10*RS_11 + 
                            8.0*F*ZETA_p3*ZETA_10*((-E0_2RS + 
                                                    E1_2RS)*RS_01*RS_10 + (-E0_1RS + E1_1RS)*RS_11) + 
                            2.0*ZETA_p4*F_1ZETA*ZETA_10*((-E0_2RS + E1_2RS)*RS_01*RS_10 + 
                                                         (-E0_1RS + E1_1RS)*RS_11) + 
                            8.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*RS_10*ZETA_11 + 
                            2.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*RS_10*ZETA_11 + 
                            8.0*iF0*F*ZETA_p3*minA_1RS*RS_10*ZETA_11 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_10*ZETA_11 + 
                            24.0*(-E0 + E1)*F*ZETA_p2*ZETA_10*ZETA_11 + 
                            24.0*iF0*F*ZETA_p2*minA*ZETA_10*ZETA_11 + 
                            16.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_10*ZETA_11 + 
                            16.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_10*ZETA_11 + 
                            2.0*(-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_10*ZETA_11 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_10*ZETA_11 + 
                            E0_2RS*RS_01*RS_20 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_01*RS_20 + 
                            4.0*iF0*F*ZETA_p3*minA_1RS*ZETA_01*RS_20 + 
                            iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*ZETA_01*RS_20 + 
                            4.0*F*ZETA_p3*ZETA_01*((-E0_2RS + E1_2RS)*RS_10_p2 + 
                                                   (-E0_1RS + E1_1RS)*RS_20) + 
                            ZETA_p4*F_1ZETA*ZETA_01*((-E0_2RS + E1_2RS)*RS_10_p2 + 
                                                     (-E0_1RS + E1_1RS)*RS_20) + 
                            4.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*RS_01*ZETA_20 + 
                            ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_20 + 
                            4.0*iF0*F*ZETA_p3*minA_1RS*RS_01*ZETA_20 + 
                            iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_01*ZETA_20 + 
                            12.0*(-E0 + E1)*F*ZETA_p2*ZETA_01*ZETA_20 + 
                            12.0*iF0*F*ZETA_p2*minA*ZETA_01*ZETA_20 + 
                            8.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_01*ZETA_20 + 
                            8.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_01*ZETA_20 + 
                            (-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_01*ZETA_20 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_01*ZETA_20 + 
                            E0_1RS*RS_21 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_1RS*RS_21 + 
                            F*ZETA_p4*(-(E0_3RS*RS_01*RS_10_p2) + 
                                       E1_3RS*RS_01*RS_10_p2-(E0_2RS-E1_2RS)*(2.0*RS_10*RS_11 + 
                                                                              RS_01*RS_20) + 
                                       (-E0_1RS + E1_1RS)*RS_21) + 
                            4.0*(-E0 + E1)*F*ZETA_p3*ZETA_21 + 
                            4.0*iF0*F*ZETA_p3*minA*ZETA_21 + 
                            (-E0 + E1)*ZETA_p4*F_1ZETA*ZETA_21 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_1ZETA*ZETA_21
	);



    ds->df12000 += factor* ( 
                             E0_3RS*RS_01_p2*RS_10 + 
                             iF0*F*(-1.0 + ZETA_p4)*minA_3RS*RS_01_p2*RS_10 + 
                             8.0*iF0*F*ZETA_p3*minA_2RS*RS_01*ZETA_01*RS_10 + 
                             2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_01*ZETA_01*RS_10 + 
                             12.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*ZETA_01_p2*RS_10 + 
                             8.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*ZETA_01_p2*RS_10 + 
                             12.0*iF0*F*ZETA_p2*minA_1RS*ZETA_01_p2*RS_10 + 
                             8.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*ZETA_01_p2*RS_10 + 
                             ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*ZETA_01_p2*RS_10 + 
                             iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*ZETA_01_p2*RS_10 + 
                             E0_2RS*RS_02*RS_10 + 
                             iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_02*RS_10 + 
                             4.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*ZETA_02*RS_10 + 
                             ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*ZETA_02*RS_10 + 
                             4.0*iF0*F*ZETA_p3*minA_1RS*ZETA_02*RS_10 + 
                             iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*ZETA_02*RS_10 + 
                             4.0*iF0*F*ZETA_p3*minA_2RS*RS_01_p2*ZETA_10 + 
                             iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_01_p2*ZETA_10 + 
                             24.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*RS_01*ZETA_01*ZETA_10 + 
                             16.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_01*ZETA_10 + 
                             24.0*iF0*F*ZETA_p2*minA_1RS*RS_01*ZETA_01*ZETA_10 + 
                             16.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*RS_01*ZETA_01*ZETA_10 + 
                             2.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*RS_01*ZETA_01*ZETA_10 + 
                             2.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*RS_01*ZETA_01*ZETA_10 + 
                             24.0*(-E0 + E1)*F*ZETA*ZETA_01_p2*ZETA_10 + 
                             24.0*iF0*F*ZETA*minA*ZETA_01_p2*ZETA_10 + 
                             36.0*(-E0 + E1)*ZETA_p2*F_1ZETA*ZETA_01_p2*ZETA_10 + 
                             36.0*iF0*ZETA_p2*minA*F_1ZETA*ZETA_01_p2*ZETA_10 + 
                             12.0*(-E0 + E1)*ZETA_p3*F_2ZETA*ZETA_01_p2*ZETA_10 + 
                             12.0*iF0*ZETA_p3*minA*F_2ZETA*ZETA_01_p2*ZETA_10 + 
                             (-E0 + E1)*ZETA_p4*F_3ZETA*ZETA_01_p2*ZETA_10 + 
                             iF0*(-1.0 + ZETA_p4)*minA*F_3ZETA*ZETA_01_p2*ZETA_10 + 
                             4.0*iF0*F*ZETA_p3*minA_1RS*RS_02*ZETA_10 + 
                             iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_02*ZETA_10 + 
                             4.0*F*ZETA_p3*((-E0_2RS + E1_2RS)*RS_01_p2 + 
                                            (-E0_1RS + E1_1RS)*RS_02)*ZETA_10 + 
                             ZETA_p4*F_1ZETA*((-E0_2RS + E1_2RS)*RS_01_p2 + 
                                              (-E0_1RS + E1_1RS)*RS_02)*ZETA_10 + 
                             12.0*(-E0 + E1)*F*ZETA_p2*ZETA_02*ZETA_10 + 
                             12.0*iF0*F*ZETA_p2*minA*ZETA_02*ZETA_10 + 
                             8.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_02*ZETA_10 + 
                             8.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_02*ZETA_10 + 
                             (-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_02*ZETA_10 + 
                             iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_02*ZETA_10 + 
                             2.0*E0_2RS*RS_01*RS_11 + 
                             2.0*iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_01*RS_11 + 
                             8.0*iF0*F*ZETA_p3*minA_1RS*ZETA_01*RS_11 + 
                             2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*ZETA_01*RS_11 + 
                             8.0*F*ZETA_p3*ZETA_01*((-E0_2RS + E1_2RS)*RS_01*RS_10 + 
                                                    (-E0_1RS + E1_1RS)*RS_11) + 
                             2.0*ZETA_p4*F_1ZETA*ZETA_01*((-E0_2RS + E1_2RS)*RS_01*RS_10 + 
                                                          (-E0_1RS + E1_1RS)*RS_11) + 
                             8.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*RS_01*ZETA_11 + 
                             2.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_11 + 
                             8.0*iF0*F*ZETA_p3*minA_1RS*RS_01*ZETA_11 + 
                             2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_01*ZETA_11 + 
                             24.0*(-E0 + E1)*F*ZETA_p2*ZETA_01*ZETA_11 + 
                             24.0*iF0*F*ZETA_p2*minA*ZETA_01*ZETA_11 + 16.0*(-E0 + 
                                                                             E1)*ZETA_p3*F_1ZETA*ZETA_01*ZETA_11 + 
                             16.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_01*ZETA_11 + 
                             2.0*(-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_01*ZETA_11 + 
                             2.0*iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_01*ZETA_11 + 
                             E0_1RS*RS_12 + 
                             iF0*F*(-1.0 + ZETA_p4)*minA_1RS*RS_12 + 
                             F*ZETA_p4*(-(((E0_3RS-E1_3RS)*RS_01_p2 + 
                                           (E0_2RS-E1_2RS)*RS_02)*RS_10)-2.0*(E0_2RS-E1_2RS)*RS_01*RS_11 + 
                                        (-E0_1RS + E1_1RS)*RS_12) + 
                             4.0*(-E0 + E1)*F*ZETA_p3*ZETA_12 + 
                             4.0*iF0*F*ZETA_p3*minA*ZETA_12 + 
                             (-E0 + E1)*ZETA_p4*F_1ZETA*ZETA_12 + 
                             iF0*(-1.0 + ZETA_p4)*minA*F_1ZETA*ZETA_12
	);



    ds->df40000 += factor* (
                            E0_4RS*RS_10_p4 + 
                            iF0*F*(-1.0 +ZETA_p4)*minA_4RS*RS_10_p4 + 
                            16.0*iF0*F*ZETA_p3*minA_3RS*RS_10_p3*ZETA_10 + 
                            4.0*iF0*(-1.0 +ZETA_p4)*F_1ZETA*minA_3RS*RS_10_p3*ZETA_10 + 
                            72.0*iF0*F*ZETA_p2*minA_2RS*RS_10_p2*ZETA_10_p2 + 
                            48.0*iF0*ZETA_p3*F_1ZETA*minA_2RS*RS_10_p2*ZETA_10_p2 + 
                            6.0*iF0*(-1.0 +ZETA_p4)*F_2ZETA*minA_2RS*RS_10_p2*ZETA_10_p2 + 
                            96.0*F*ZETA*(-E0_1RS +E1_1RS)*RS_10*ZETA_10_p3 + 
                            144.0*ZETA_p2*(-E0_1RS +E1_1RS)*F_1ZETA*RS_10*ZETA_10_p3 + 
                            96.0*iF0*F*ZETA*minA_1RS*RS_10*ZETA_10_p3 + 
                            144.0*iF0*ZETA_p2*F_1ZETA*minA_1RS*RS_10*ZETA_10_p3 + 
                            48.0*ZETA_p3*(-E0_1RS +E1_1RS)*F_2ZETA*RS_10*ZETA_10_p3 + 
                            48.0*iF0*ZETA_p3*minA_1RS*F_2ZETA*RS_10*ZETA_10_p3 + 
                            4.0*ZETA_p4*(-E0_1RS +E1_1RS)*F_3ZETA*RS_10*ZETA_10_p3 + 
                            4.0*iF0*(-1.0 +ZETA_p4)*minA_1RS*F_3ZETA*RS_10*ZETA_10_p3 + 
                            24.0*(-E0 +E1)*F*ZETA_10_p4 + 
                            24.0*iF0*F*minA*ZETA_10_p4 + 
                            96.0*(-E0 +E1)*ZETA*F_1ZETA*ZETA_10_p4 + 
                            96.0*iF0*ZETA*minA*F_1ZETA*ZETA_10_p4 + 
                            72.0*(-E0 +E1)*ZETA_p2*F_2ZETA*ZETA_10_p4 + 
                            72.0*iF0*ZETA_p2*minA*F_2ZETA*ZETA_10_p4 + 
                            16.0*(-E0 +E1)*ZETA_p3*F_3ZETA*ZETA_10_p4 + 
                            16.0*iF0*ZETA_p3*minA*F_3ZETA*ZETA_10_p4 + 
                            (-E0 +E1)*ZETA_p4*F_4ZETA*ZETA_10_p4 + 
                            iF0*(-1.0 +ZETA_p4)*minA*F_4ZETA*ZETA_10_p4 + 
                            6.0*E0_3RS*RS_10_p2*RS_20 + 
                            6.0*iF0*F*(-1.0 +ZETA_p4)*minA_3RS*RS_10_p2*RS_20 + 
                            48.0*iF0*F*ZETA_p3*minA_2RS*RS_10*ZETA_10*RS_20 + 
                            12.0*iF0*(-1.0 +ZETA_p4)*F_1ZETA*minA_2RS*RS_10*ZETA_10*RS_20 + 
                            72.0*iF0*F*ZETA_p2*minA_1RS*ZETA_10_p2*RS_20 + 
                            48.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*ZETA_10_p2*RS_20 + 
                            6.0*iF0*(-1.0 +ZETA_p4)*minA_1RS*F_2ZETA*ZETA_10_p2*RS_20 + 
                            3.0*E0_2RS*RS_20_p2 + 
                            3.0*iF0*F*(-1.0 +ZETA_p4)*minA_2RS*RS_20_p2 + 
                            72.0*F*ZETA_p2*ZETA_10_p2*((-E0_2RS + 
                                                        E1_2RS)*RS_10_p2 + 
                                                       (-E0_1RS +E1_1RS)*RS_20) + 
                            48.0*ZETA_p3*F_1ZETA*ZETA_10_p2*((-E0_2RS + 
                                                              E1_2RS)*RS_10_p2 + 
                                                             (-E0_1RS +E1_1RS)*RS_20) + 
                            6.0*ZETA_p4*F_2ZETA*ZETA_10_p2*((-E0_2RS + 
                                                             E1_2RS)*RS_10_p2 + 
                                                            (-E0_1RS +E1_1RS)*RS_20) + 
                            24.0*iF0*F*ZETA_p3*minA_2RS*RS_10_p2*ZETA_20 + 
                            6.0*iF0*(-1.0 +ZETA_p4)*F_1ZETA*minA_2RS*RS_10_p2*ZETA_20 + 
                            144.0*F*ZETA_p2*(-E0_1RS +E1_1RS)*RS_10*ZETA_10*ZETA_20 + 
                            96.0*ZETA_p3*(-E0_1RS +E1_1RS)*F_1ZETA*RS_10*ZETA_10*ZETA_20 + 
                            144.0*iF0*F*ZETA_p2*minA_1RS*RS_10*ZETA_10*ZETA_20 + 
                            96.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*RS_10*ZETA_10*ZETA_20 + 
                            12.0*ZETA_p4*(-E0_1RS +E1_1RS)*F_2ZETA*RS_10*ZETA_10*ZETA_20 + 
                            12.0*iF0*(-1.0 +ZETA_p4)*minA_1RS*F_2ZETA*RS_10*ZETA_10*ZETA_20 + 
                            144.0*(-E0 +E1)*F*ZETA*ZETA_10_p2*ZETA_20 + 
                            144.0*iF0*F*ZETA*minA*ZETA_10_p2*ZETA_20 + 
                            216.0*(-E0 +E1)*ZETA_p2*F_1ZETA*ZETA_10_p2*ZETA_20 + 
                            216.0*iF0*ZETA_p2*minA*F_1ZETA*ZETA_10_p2*ZETA_20 + 
                            72.0*(-E0 +E1)*ZETA_p3*F_2ZETA*ZETA_10_p2*ZETA_20 + 
                            72.0*iF0*ZETA_p3*minA*F_2ZETA*ZETA_10_p2*ZETA_20 + 
                            6.0*(-E0 +E1)*ZETA_p4*F_3ZETA*ZETA_10_p2*ZETA_20 + 
                            6.0*iF0*(-1.0 +ZETA_p4)*minA*F_3ZETA*ZETA_10_p2*ZETA_20 + 
                            24.0*iF0*F*ZETA_p3*minA_1RS*RS_20*ZETA_20 + 
                            6.0*iF0*(-1.0 +ZETA_p4)*F_1ZETA*minA_1RS*RS_20*ZETA_20 + 
                            24.0*F*ZETA_p3*((-E0_2RS + 
                                             E1_2RS)*RS_10_p2 + 
                                            (-E0_1RS +E1_1RS)*RS_20)*ZETA_20 + 
                            6.0*ZETA_p4*F_1ZETA*((-E0_2RS + 
                                                  E1_2RS)*RS_10_p2 + 
                                                 (-E0_1RS +E1_1RS)*RS_20)*ZETA_20 + 
                            36.0*(-E0 +E1)*F*ZETA_p2*ZETA_20_p2 + 
                            36.0*iF0*F*ZETA_p2*minA*ZETA_20_p2 + 
                            24.0*(-E0 +E1)*ZETA_p3*F_1ZETA*ZETA_20_p2 + 
                            24.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_20_p2 + 
                            3.0*(-E0 +E1)*ZETA_p4*F_2ZETA*ZETA_20_p2 + 
                            3.0*iF0*(-1.0 +ZETA_p4)*minA*F_2ZETA*ZETA_20_p2 + 
                            4.0*E0_2RS*RS_10*RS_30 + 
                            4.0*iF0*F*(-1.0 +ZETA_p4)*minA_2RS*RS_10*RS_30 + 
                            16.0*iF0*F*ZETA_p3*minA_1RS*ZETA_10*RS_30 + 
                            4.0*iF0*(-1.0 +ZETA_p4)*F_1ZETA*minA_1RS*ZETA_10*RS_30 + 
                            16.0*F*ZETA_p3*ZETA_10*(-((E0_3RS-E1_3RS)*RS_10_p3)-3.0*(E0_2RS-E1_2RS)*RS_10*RS_20 + 
                                                    (-E0_1RS +E1_1RS)*RS_30) + 
                            4.0*ZETA_p4*F_1ZETA*ZETA_10*(-((E0_3RS-E1_3RS)*RS_10_p3)-3.0*(E0_2RS-E1_2RS)*RS_10*RS_20 + 
                                                         (-E0_1RS +E1_1RS)*RS_30) + 
                            16.0*F*ZETA_p3*(-E0_1RS +E1_1RS)*RS_10*ZETA_30 + 
                            4.0*ZETA_p4*(-E0_1RS +E1_1RS)*F_1ZETA*RS_10*ZETA_30 + 
                            16.0*iF0*F*ZETA_p3*minA_1RS*RS_10*ZETA_30 + 
                            4.0*iF0*(-1.0 +ZETA_p4)*F_1ZETA*minA_1RS*RS_10*ZETA_30 + 
                            48.0*(-E0 +E1)*F*ZETA_p2*ZETA_10*ZETA_30 + 
                            48.0*iF0*F*ZETA_p2*minA*ZETA_10*ZETA_30 + 
                            32.0*(-E0 +E1)*ZETA_p3*F_1ZETA*ZETA_10*ZETA_30 + 
                            32.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_10*ZETA_30 + 
                            4.0*(-E0 +E1)*ZETA_p4*F_2ZETA*ZETA_10*ZETA_30 + 
                            4.0*iF0*(-1.0 +ZETA_p4)*minA*F_2ZETA*ZETA_10*ZETA_30 + 
                            E0_1RS*RS_40 + 
                            iF0*F*(-1.0 +ZETA_p4)*minA_1RS*RS_40 + 
                            F*ZETA_p4*(-(E0_4RS*RS_10_p4) + 
                                       E1_4RS*RS_10_p4-6.0*(E0_3RS-E1_3RS)*RS_10_p2*RS_20-3.0*(E0_2RS-E1_2RS)*RS_20_p2-4.0*(E0_2RS-E1_2RS)*RS_10*RS_30 + 
                                       (-E0_1RS +E1_1RS)*RS_40) + 
                            4.0*(-E0 +E1)*F*ZETA_p3*ZETA_40 + 
                            4.0*iF0*F*ZETA_p3*minA*ZETA_40 + 
                            (-E0 +E1)*ZETA_p4*F_1ZETA*ZETA_40 + 
                            iF0*(-1.0 +ZETA_p4)*minA*F_1ZETA*ZETA_40
	);

    ds->df04000 += factor* (
                            E0_4RS*RS_01_p4 + 
                            iF0*F*(-1.0 +ZETA_p4)*minA_4RS*RS_01_p4 + 
                            16.0*iF0*F*ZETA_p3*minA_3RS*RS_01_p3*ZETA_01 + 
                            4.0*iF0*(-1.0 +ZETA_p4)*F_1ZETA*minA_3RS*RS_01_p3*ZETA_01 + 
                            72.0*iF0*F*ZETA_p2*minA_2RS*RS_01_p2*ZETA_01_p2 + 
                            48.0*iF0*ZETA_p3*F_1ZETA*minA_2RS*RS_01_p2*ZETA_01_p2 + 
                            6.0*iF0*(-1.0 +ZETA_p4)*F_2ZETA*minA_2RS*RS_01_p2*ZETA_01_p2 + 
                            96.0*F*ZETA*(-E0_1RS +E1_1RS)*RS_01*ZETA_01_p3 + 
                            144.0*ZETA_p2*(-E0_1RS +E1_1RS)*F_1ZETA*RS_01*ZETA_01_p3 + 
                            96.0*iF0*F*ZETA*minA_1RS*RS_01*ZETA_01_p3 + 
                            144.0*iF0*ZETA_p2*F_1ZETA*minA_1RS*RS_01*ZETA_01_p3 + 
                            48.0*ZETA_p3*(-E0_1RS +E1_1RS)*F_2ZETA*RS_01*ZETA_01_p3 + 
                            48.0*iF0*ZETA_p3*minA_1RS*F_2ZETA*RS_01*ZETA_01_p3 + 
                            4.0*ZETA_p4*(-E0_1RS +E1_1RS)*F_3ZETA*RS_01*ZETA_01_p3 + 
                            4.0*iF0*(-1.0 +ZETA_p4)*minA_1RS*F_3ZETA*RS_01*ZETA_01_p3 + 
                            24.0*(-E0 +E1)*F*ZETA_01_p4 + 
                            24.0*iF0*F*minA*ZETA_01_p4 + 
                            96.0*(-E0 +E1)*ZETA*F_1ZETA*ZETA_01_p4 + 
                            96.0*iF0*ZETA*minA*F_1ZETA*ZETA_01_p4 + 
                            72.0*(-E0 +E1)*ZETA_p2*F_2ZETA*ZETA_01_p4 + 
                            72.0*iF0*ZETA_p2*minA*F_2ZETA*ZETA_01_p4 + 
                            16.0*(-E0 +E1)*ZETA_p3*F_3ZETA*ZETA_01_p4 + 
                            16.0*iF0*ZETA_p3*minA*F_3ZETA*ZETA_01_p4 + 
                            (-E0 +E1)*ZETA_p4*F_4ZETA*ZETA_01_p4 + 
                            iF0*(-1.0 +ZETA_p4)*minA*F_4ZETA*ZETA_01_p4 + 
                            6.0*E0_3RS*RS_01_p2*RS_02 + 
                            6.0*iF0*F*(-1.0 +ZETA_p4)*minA_3RS*RS_01_p2*RS_02 + 
                            48.0*iF0*F*ZETA_p3*minA_2RS*RS_01*ZETA_01*RS_02 + 
                            12.0*iF0*(-1.0 +ZETA_p4)*F_1ZETA*minA_2RS*RS_01*ZETA_01*RS_02 + 
                            72.0*iF0*F*ZETA_p2*minA_1RS*ZETA_01_p2*RS_02 + 
                            48.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*ZETA_01_p2*RS_02 + 
                            6.0*iF0*(-1.0 +ZETA_p4)*minA_1RS*F_2ZETA*ZETA_01_p2*RS_02 + 
                            3.0*E0_2RS*RS_02_p2 + 
                            3.0*iF0*F*(-1.0 +ZETA_p4)*minA_2RS*RS_02_p2 + 
                            72.0*F*ZETA_p2*ZETA_01_p2*((-E0_2RS + 
                                                        E1_2RS)*RS_01_p2 + 
                                                       (-E0_1RS +E1_1RS)*RS_02) + 
                            48.0*ZETA_p3*F_1ZETA*ZETA_01_p2*((-E0_2RS + 
                                                              E1_2RS)*RS_01_p2 + 
                                                             (-E0_1RS +E1_1RS)*RS_02) + 
                            6.0*ZETA_p4*F_2ZETA*ZETA_01_p2*((-E0_2RS + 
                                                             E1_2RS)*RS_01_p2 + 
                                                            (-E0_1RS +E1_1RS)*RS_02) + 
                            24.0*iF0*F*ZETA_p3*minA_2RS*RS_01_p2*ZETA_02 + 
                            6.0*iF0*(-1.0 +ZETA_p4)*F_1ZETA*minA_2RS*RS_01_p2*ZETA_02 + 
                            144.0*F*ZETA_p2*(-E0_1RS +E1_1RS)*RS_01*ZETA_01*ZETA_02 + 
                            96.0*ZETA_p3*(-E0_1RS +E1_1RS)*F_1ZETA*RS_01*ZETA_01*ZETA_02 + 
                            144.0*iF0*F*ZETA_p2*minA_1RS*RS_01*ZETA_01*ZETA_02 + 
                            96.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*RS_01*ZETA_01*ZETA_02 + 
                            12.0*ZETA_p4*(-E0_1RS +E1_1RS)*F_2ZETA*RS_01*ZETA_01*ZETA_02 + 
                            12.0*iF0*(-1.0 +ZETA_p4)*minA_1RS*F_2ZETA*RS_01*ZETA_01*ZETA_02 + 
                            144.0*(-E0 +E1)*F*ZETA*ZETA_01_p2*ZETA_02 + 
                            144.0*iF0*F*ZETA*minA*ZETA_01_p2*ZETA_02 + 
                            216.0*(-E0 +E1)*ZETA_p2*F_1ZETA*ZETA_01_p2*ZETA_02 + 
                            216.0*iF0*ZETA_p2*minA*F_1ZETA*ZETA_01_p2*ZETA_02 + 
                            72.0*(-E0 +E1)*ZETA_p3*F_2ZETA*ZETA_01_p2*ZETA_02 + 
                            72.0*iF0*ZETA_p3*minA*F_2ZETA*ZETA_01_p2*ZETA_02 + 
                            6.0*(-E0 +E1)*ZETA_p4*F_3ZETA*ZETA_01_p2*ZETA_02 + 
                            6.0*iF0*(-1.0 +ZETA_p4)*minA*F_3ZETA*ZETA_01_p2*ZETA_02 + 
                            24.0*iF0*F*ZETA_p3*minA_1RS*RS_02*ZETA_02 + 
                            6.0*iF0*(-1.0 +ZETA_p4)*F_1ZETA*minA_1RS*RS_02*ZETA_02 + 
                            24.0*F*ZETA_p3*((-E0_2RS + 
                                             E1_2RS)*RS_01_p2 + 
                                            (-E0_1RS +E1_1RS)*RS_02)*ZETA_02 + 
                            6.0*ZETA_p4*F_1ZETA*((-E0_2RS + 
                                                  E1_2RS)*RS_01_p2 + 
                                                 (-E0_1RS +E1_1RS)*RS_02)*ZETA_02 + 
                            36.0*(-E0 +E1)*F*ZETA_p2*ZETA_02_p2 + 
                            36.0*iF0*F*ZETA_p2*minA*ZETA_02_p2 + 
                            24.0*(-E0 +E1)*ZETA_p3*F_1ZETA*ZETA_02_p2 + 
                            24.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_02_p2 + 
                            3.0*(-E0 +E1)*ZETA_p4*F_2ZETA*ZETA_02_p2 + 
                            3.0*iF0*(-1.0 +ZETA_p4)*minA*F_2ZETA*ZETA_02_p2 + 
                            4.0*E0_2RS*RS_01*RS_03 + 
                            4.0*iF0*F*(-1.0 +ZETA_p4)*minA_2RS*RS_01*RS_03 + 
                            16.0*iF0*F*ZETA_p3*minA_1RS*ZETA_01*RS_03 + 
                            4.0*iF0*(-1.0 +ZETA_p4)*F_1ZETA*minA_1RS*ZETA_01*RS_03 + 
                            16.0*F*ZETA_p3*ZETA_01*(-((E0_3RS-E1_3RS)*RS_01_p3)-3.0*(E0_2RS-E1_2RS)*RS_01*RS_02 + 
                                                    (-E0_1RS +E1_1RS)*RS_03) + 
                            4.0*ZETA_p4*F_1ZETA*ZETA_01*(-((E0_3RS-E1_3RS)*RS_01_p3)-3.0*(E0_2RS-E1_2RS)*RS_01*RS_02 + 
                                                         (-E0_1RS +E1_1RS)*RS_03) + 
                            16.0*F*ZETA_p3*(-E0_1RS +E1_1RS)*RS_01*ZETA_03 + 
                            4.0*ZETA_p4*(-E0_1RS +E1_1RS)*F_1ZETA*RS_01*ZETA_03 + 
                            16.0*iF0*F*ZETA_p3*minA_1RS*RS_01*ZETA_03 + 
                            4.0*iF0*(-1.0 +ZETA_p4)*F_1ZETA*minA_1RS*RS_01*ZETA_03 + 
                            48.0*(-E0 +E1)*F*ZETA_p2*ZETA_01*ZETA_03 + 
                            48.0*iF0*F*ZETA_p2*minA*ZETA_01*ZETA_03 + 
                            32.0*(-E0 +E1)*ZETA_p3*F_1ZETA*ZETA_01*ZETA_03 + 
                            32.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_01*ZETA_03 + 
                            4.0*(-E0 +E1)*ZETA_p4*F_2ZETA*ZETA_01*ZETA_03 + 
                            4.0*iF0*(-1.0 +ZETA_p4)*minA*F_2ZETA*ZETA_01*ZETA_03 + 
                            E0_1RS*RS_04 + 
                            iF0*F*(-1.0 +ZETA_p4)*minA_1RS*RS_04 + 
                            F*ZETA_p4*(-(E0_4RS*RS_01_p4) + 
                                       E1_4RS*RS_01_p4-6.0*(E0_3RS-E1_3RS)*RS_01_p2*RS_02-3.0*(E0_2RS-E1_2RS)*RS_02_p2-4.0*(E0_2RS-E1_2RS)*RS_01*RS_03 + 
                                       (-E0_1RS +E1_1RS)*RS_04) + 
                            4.0*(-E0 +E1)*F*ZETA_p3*ZETA_04 + 
                            4.0*iF0*F*ZETA_p3*minA*ZETA_04 + 
                            (-E0 +E1)*ZETA_p4*F_1ZETA*ZETA_04 + 
                            iF0*(-1.0 +ZETA_p4)*minA*F_1ZETA*ZETA_04
	);

    ds->df31000 += factor* (
                            E0_4RS*RS_01*RS_10_p3 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_4RS*RS_01*RS_10_p3 + 
                            4.0*iF0*F*ZETA_p3*minA_3RS*ZETA_01*RS_10_p3 + 
                            iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_3RS*ZETA_01*RS_10_p3 + 
                            12.0*iF0*F*ZETA_p3*minA_3RS*RS_01*RS_10_p2*ZETA_10 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_3RS*RS_01*RS_10_p2*ZETA_10 + 
                            36.0*iF0*F*ZETA_p2*minA_2RS*ZETA_01*RS_10_p2*ZETA_10 + 
                            24.0*iF0*ZETA_p3*F_1ZETA*minA_2RS*ZETA_01*RS_10_p2*ZETA_10 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_2ZETA*minA_2RS*ZETA_01*RS_10_p2*ZETA_10 + 
                            36.0*iF0*F*ZETA_p2*minA_2RS*RS_01*RS_10*ZETA_10_p2 + 
                            24.0*iF0*ZETA_p3*F_1ZETA*minA_2RS*RS_01*RS_10*ZETA_10_p2 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_2ZETA*minA_2RS*RS_01*RS_10*ZETA_10_p2 + 
                            72.0*F*ZETA*(-E0_1RS + E1_1RS)*ZETA_01*RS_10*ZETA_10_p2 + 
                            108.0*ZETA_p2*(-E0_1RS + E1_1RS)*F_1ZETA*ZETA_01*RS_10*ZETA_10_p2 + 
                            72.0*iF0*F*ZETA*minA_1RS*ZETA_01*RS_10*ZETA_10_p2 + 
                            108.0*iF0*ZETA_p2*F_1ZETA*minA_1RS*ZETA_01*RS_10*ZETA_10_p2 + 
                            36.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_2ZETA*ZETA_01*RS_10*ZETA_10_p2 + 
                            36.0*iF0*ZETA_p3*minA_1RS*F_2ZETA*ZETA_01*RS_10*ZETA_10_p2 + 
                            3.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_3ZETA*ZETA_01*RS_10*ZETA_10_p2 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_3ZETA*ZETA_01*RS_10*ZETA_10_p2 + 
                            24.0*F*ZETA*(-E0_1RS + E1_1RS)*RS_01*ZETA_10_p3 + 
                            36.0*ZETA_p2*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_10_p3 + 
                            24.0*iF0*F*ZETA*minA_1RS*RS_01*ZETA_10_p3 + 
                            36.0*iF0*ZETA_p2*F_1ZETA*minA_1RS*RS_01*ZETA_10_p3 + 
                            12.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_2ZETA*RS_01*ZETA_10_p3 + 
                            12.0*iF0*ZETA_p3*minA_1RS*F_2ZETA*RS_01*ZETA_10_p3 + 
                            ZETA_p4*(-E0_1RS + E1_1RS)*F_3ZETA*RS_01*ZETA_10_p3 + 
                            iF0*(-1.0 + ZETA_p4)*minA_1RS*F_3ZETA*RS_01*ZETA_10_p3 + 
                            24.0*(-E0 + E1)*F*ZETA_01*ZETA_10_p3 + 
                            24.0*iF0*F*minA*ZETA_01*ZETA_10_p3 + 
                            96.0*(-E0 + E1)*ZETA*F_1ZETA*ZETA_01*ZETA_10_p3 + 
                            96.0*iF0*ZETA*minA*F_1ZETA*ZETA_01*ZETA_10_p3 + 
                            72.0*(-E0 + E1)*ZETA_p2*F_2ZETA*ZETA_01*ZETA_10_p3 + 
                            72.0*iF0*ZETA_p2*minA*F_2ZETA*ZETA_01*ZETA_10_p3 + 
                            16.0*(-E0 + E1)*ZETA_p3*F_3ZETA*ZETA_01*ZETA_10_p3 + 
                            16.0*iF0*ZETA_p3*minA*F_3ZETA*ZETA_01*ZETA_10_p3 + 
                            (-E0 + E1)*ZETA_p4*F_4ZETA*ZETA_01*ZETA_10_p3 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_4ZETA*ZETA_01*ZETA_10_p3 + 
                            3.0*E0_3RS*RS_10_p2*RS_11 + 
                            3.0*iF0*F*(-1.0 + ZETA_p4)*minA_3RS*RS_10_p2*RS_11 + 
                            24.0*iF0*F*ZETA_p3*minA_2RS*RS_10*ZETA_10*RS_11 + 
                            6.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_10*ZETA_10*RS_11 + 
                            36.0*iF0*F*ZETA_p2*minA_1RS*ZETA_10_p2*RS_11 + 
                            24.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*ZETA_10_p2*RS_11 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*ZETA_10_p2*RS_11 + 
                            36.0*F*ZETA_p2*ZETA_10_p2*((-E0_2RS + E1_2RS)*RS_01*RS_10 + 
                                                       (-E0_1RS + E1_1RS)*RS_11) + 
                            24.0*ZETA_p3*F_1ZETA*ZETA_10_p2*((-E0_2RS + E1_2RS)*RS_01*RS_10 + 
                                                             (-E0_1RS + E1_1RS)*RS_11) + 
                            3.0*ZETA_p4*F_2ZETA*ZETA_10_p2*((-E0_2RS + E1_2RS)*RS_01*RS_10 + 
                                                            (-E0_1RS + E1_1RS)*RS_11) + 
                            12.0*iF0*F*ZETA_p3*minA_2RS*RS_10_p2*ZETA_11 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_10_p2*ZETA_11 + 
                            72.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*RS_10*ZETA_10*ZETA_11 + 
                            48.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*RS_10*ZETA_10*ZETA_11 + 
                            72.0*iF0*F*ZETA_p2*minA_1RS*RS_10*ZETA_10*ZETA_11 + 
                            48.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*RS_10*ZETA_10*ZETA_11 + 
                            6.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*RS_10*ZETA_10*ZETA_11 + 
                            6.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*RS_10*ZETA_10*ZETA_11 + 
                            72.0*(-E0 + E1)*F*ZETA*ZETA_10_p2*ZETA_11 + 
                            72.0*iF0*F*ZETA*minA*ZETA_10_p2*ZETA_11 + 
                            108.0*(-E0 + E1)*ZETA_p2*F_1ZETA*ZETA_10_p2*ZETA_11 + 
                            108.0*iF0*ZETA_p2*minA*F_1ZETA*ZETA_10_p2*ZETA_11 + 
                            36.0*(-E0 + E1)*ZETA_p3*F_2ZETA*ZETA_10_p2*ZETA_11 + 
                            36.0*iF0*ZETA_p3*minA*F_2ZETA*ZETA_10_p2*ZETA_11 + 
                            3.0*(-E0 + E1)*ZETA_p4*F_3ZETA*ZETA_10_p2*ZETA_11 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA*F_3ZETA*ZETA_10_p2*ZETA_11 + 
                            3.0*E0_3RS*RS_01*RS_10*RS_20 + 
                            3.0*iF0*F*(-1.0 + ZETA_p4)*minA_3RS*RS_01*RS_10*RS_20 + 
                            12.0*iF0*F*ZETA_p3*minA_2RS*ZETA_01*RS_10*RS_20 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*ZETA_01*RS_10*RS_20 + 
                            12.0*iF0*F*ZETA_p3*minA_2RS*RS_01*ZETA_10*RS_20 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_01*ZETA_10*RS_20 + 
                            36.0*iF0*F*ZETA_p2*minA_1RS*ZETA_01*ZETA_10*RS_20 + 
                            24.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*ZETA_01*ZETA_10*RS_20 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*ZETA_01*ZETA_10*RS_20 + 
                            3.0*E0_2RS*RS_11*RS_20 + 
                            3.0*iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_11*RS_20 + 
                            12.0*iF0*F*ZETA_p3*minA_1RS*ZETA_11*RS_20 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*ZETA_11*RS_20 + 
                            36.0*F*ZETA_p2*ZETA_01*ZETA_10*((-E0_2RS + E1_2RS)*RS_10_p2 + 
                                                            (-E0_1RS + E1_1RS)*RS_20) + 
                            24.0*ZETA_p3*F_1ZETA*ZETA_01*ZETA_10*((-E0_2RS + E1_2RS)*RS_10_p2 + 
                                                                  (-E0_1RS + E1_1RS)*RS_20) + 
                            3.0*ZETA_p4*F_2ZETA*ZETA_01*ZETA_10*((-E0_2RS + E1_2RS)*RS_10_p2 + 
                                                                 (-E0_1RS + E1_1RS)*RS_20) + 
                            12.0*F*ZETA_p3*ZETA_11*((-E0_2RS + E1_2RS)*RS_10_p2 + 
                                                    (-E0_1RS + E1_1RS)*RS_20) + 
                            3.0*ZETA_p4*F_1ZETA*ZETA_11*((-E0_2RS + E1_2RS)*RS_10_p2 + 
                                                         (-E0_1RS + E1_1RS)*RS_20) + 
                            12.0*iF0*F*ZETA_p3*minA_2RS*RS_01*RS_10*ZETA_20 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_01*RS_10*ZETA_20 + 
                            36.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*ZETA_01*RS_10*ZETA_20 + 
                            24.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*ZETA_01*RS_10*ZETA_20 + 
                            36.0*iF0*F*ZETA_p2*minA_1RS*ZETA_01*RS_10*ZETA_20 + 
                            24.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*ZETA_01*RS_10*ZETA_20 + 
                            3.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*ZETA_01*RS_10*ZETA_20 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*ZETA_01*RS_10*ZETA_20 + 
                            36.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*RS_01*ZETA_10*ZETA_20 + 
                            24.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_10*ZETA_20 + 
                            36.0*iF0*F*ZETA_p2*minA_1RS*RS_01*ZETA_10*ZETA_20 + 
                            24.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*RS_01*ZETA_10*ZETA_20 + 
                            3.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*RS_01*ZETA_10*ZETA_20 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*RS_01*ZETA_10*ZETA_20 + 
                            72.0*(-E0 + E1)*F*ZETA*ZETA_01*ZETA_10*ZETA_20 + 
                            72.0*iF0*F*ZETA*minA*ZETA_01*ZETA_10*ZETA_20 + 
                            108.0*(-E0 + E1)*ZETA_p2*F_1ZETA*ZETA_01*ZETA_10*ZETA_20 + 
                            108.0*iF0*ZETA_p2*minA*F_1ZETA*ZETA_01*ZETA_10*ZETA_20 + 
                            36.0*(-E0 + E1)*ZETA_p3*F_2ZETA*ZETA_01*ZETA_10*ZETA_20 + 
                            36.0*iF0*ZETA_p3*minA*F_2ZETA*ZETA_01*ZETA_10*ZETA_20 + 
                            3.0*(-E0 + E1)*ZETA_p4*F_3ZETA*ZETA_01*ZETA_10*ZETA_20 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA*F_3ZETA*ZETA_01*ZETA_10*ZETA_20 + 
                            12.0*iF0*F*ZETA_p3*minA_1RS*RS_11*ZETA_20 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_11*ZETA_20 + 
                            12.0*F*ZETA_p3*((-E0_2RS + E1_2RS)*RS_01*RS_10 + 
                                            (-E0_1RS + E1_1RS)*RS_11)*ZETA_20 + 
                            3.0*ZETA_p4*F_1ZETA*((-E0_2RS + E1_2RS)*RS_01*RS_10 + 
                                                 (-E0_1RS + E1_1RS)*RS_11)*ZETA_20 + 
                            36.0*(-E0 + E1)*F*ZETA_p2*ZETA_11*ZETA_20 + 
                            36.0*iF0*F*ZETA_p2*minA*ZETA_11*ZETA_20 + 
                            24.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_11*ZETA_20 + 
                            24.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_11*ZETA_20 + 
                            3.0*(-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_11*ZETA_20 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_11*ZETA_20 + 
                            3.0*E0_2RS*RS_10*RS_21 + 
                            3.0*iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_10*RS_21 + 
                            12.0*iF0*F*ZETA_p3*minA_1RS*ZETA_10*RS_21 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*ZETA_10*RS_21 + 
                            12.0*F*ZETA_p3*ZETA_10*(-(E0_3RS*RS_01*RS_10_p2) + 
                                                    E1_3RS*RS_01*RS_10_p2-(E0_2RS-E1_2RS)*(2.0*RS_10*RS_11 + 
                                                                                           RS_01*RS_20) + 
                                                    (-E0_1RS + E1_1RS)*RS_21) + 
                            3.0*ZETA_p4*F_1ZETA*ZETA_10*(-(E0_3RS*RS_01*RS_10_p2) + 
                                                         E1_3RS*RS_01*RS_10_p2-(E0_2RS-E1_2RS)*(2.0*RS_10*RS_11 + 
                                                                                                RS_01*RS_20) + 
                                                         (-E0_1RS + E1_1RS)*RS_21) + 
                            12.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*RS_10*ZETA_21 + 
                            3.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*RS_10*ZETA_21 + 
                            12.0*iF0*F*ZETA_p3*minA_1RS*RS_10*ZETA_21 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_10*ZETA_21 + 
                            36.0*(-E0 + E1)*F*ZETA_p2*ZETA_10*ZETA_21 + 
                            36.0*iF0*F*ZETA_p2*minA*ZETA_10*ZETA_21 + 
                            24.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_10*ZETA_21 + 
                            24.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_10*ZETA_21 + 
                            3.0*(-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_10*ZETA_21 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_10*ZETA_21 + 
                            E0_2RS*RS_01*RS_30 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_01*RS_30 + 
                            4.0*iF0*F*ZETA_p3*minA_1RS*ZETA_01*RS_30 + 
                            iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*ZETA_01*RS_30 + 
                            4.0*F*ZETA_p3*ZETA_01*(-((E0_3RS-E1_3RS)*RS_10_p3)-3.0*(E0_2RS-E1_2RS)*RS_10*RS_20 + 
                                                   (-E0_1RS + E1_1RS)*RS_30) + 
                            ZETA_p4*F_1ZETA*ZETA_01*(-((E0_3RS-E1_3RS)*RS_10_p3)-3.0*(E0_2RS-E1_2RS)*RS_10*RS_20 + 
                                                     (-E0_1RS + E1_1RS)*RS_30) + 
                            4.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*RS_01*ZETA_30 + 
                            ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_30 + 
                            4.0*iF0*F*ZETA_p3*minA_1RS*RS_01*ZETA_30 + 
                            iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_01*ZETA_30 + 
                            12.0*(-E0 + E1)*F*ZETA_p2*ZETA_01*ZETA_30 + 
                            12.0*iF0*F*ZETA_p2*minA*ZETA_01*ZETA_30 + 
                            8.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_01*ZETA_30 + 
                            8.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_01*ZETA_30 + 
                            (-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_01*ZETA_30 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_01*ZETA_30 + 
                            E0_1RS*RS_31 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_1RS*RS_31 + 
                            F*ZETA_p4*(-(E0_4RS*RS_01*RS_10_p3) + 
                                       E1_4RS*RS_01*RS_10_p3 + 
                                       3.0*(-(E0_3RS*RS_10*(RS_10*RS_11 + 
                                                            RS_01*RS_20)) + 
                                            E1_3RS*RS_10*(RS_10*RS_11 + 
                                                          RS_01*RS_20)-(E0_2RS-E1_2RS)*(RS_11*RS_20 + 
                                                                                        RS_10*RS_21))-(E0_2RS-E1_2RS)*RS_01*RS_30 + 
                                       (-E0_1RS + E1_1RS)*RS_31) + 
                            4.0*(-E0 + E1)*F*ZETA_p3*ZETA_31 + 
                            4.0*iF0*F*ZETA_p3*minA*ZETA_31 + 
                            (-E0 + E1)*ZETA_p4*F_1ZETA*ZETA_31 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_1ZETA*ZETA_31
	);

    ds->df22000 += factor* (
                            E0_4RS*RS_01_p2*RS_10_p2 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_4RS*RS_01_p2*RS_10_p2 + 
                            8.0*iF0*F*ZETA_p3*minA_3RS*RS_01*ZETA_01*RS_10_p2 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_3RS*RS_01*ZETA_01*RS_10_p2 + 
                            12.0*iF0*F*ZETA_p2*minA_2RS*ZETA_01_p2*RS_10_p2 + 
                            8.0*iF0*ZETA_p3*F_1ZETA*minA_2RS*ZETA_01_p2*RS_10_p2 + 
                            iF0*(-1.0 + ZETA_p4)*F_2ZETA*minA_2RS*ZETA_01_p2*RS_10_p2 + 
                            E0_3RS*RS_02*RS_10_p2 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_3RS*RS_02*RS_10_p2 + 
                            4.0*iF0*F*ZETA_p3*minA_2RS*ZETA_02*RS_10_p2 + 
                            iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*ZETA_02*RS_10_p2 + 
                            8.0*iF0*F*ZETA_p3*minA_3RS*RS_01_p2*RS_10*ZETA_10 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_3RS*RS_01_p2*RS_10*ZETA_10 + 
                            48.0*iF0*F*ZETA_p2*minA_2RS*RS_01*ZETA_01*RS_10*ZETA_10 + 
                            32.0*iF0*ZETA_p3*F_1ZETA*minA_2RS*RS_01*ZETA_01*RS_10*ZETA_10 + 
                            4.0*iF0*(-1.0 + ZETA_p4)*F_2ZETA*minA_2RS*RS_01*ZETA_01*RS_10*ZETA_10 + 
                            48.0*F*ZETA*(-E0_1RS + E1_1RS)*ZETA_01_p2*RS_10*ZETA_10 + 
                            72.0*ZETA_p2*(-E0_1RS + E1_1RS)*F_1ZETA*ZETA_01_p2*RS_10*ZETA_10 + 
                            48.0*iF0*F*ZETA*minA_1RS*ZETA_01_p2*RS_10*ZETA_10 + 
                            72.0*iF0*ZETA_p2*F_1ZETA*minA_1RS*ZETA_01_p2*RS_10*ZETA_10 + 
                            24.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_2ZETA*ZETA_01_p2*RS_10*ZETA_10 + 
                            24.0*iF0*ZETA_p3*minA_1RS*F_2ZETA*ZETA_01_p2*RS_10*ZETA_10 + 
                            2.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_3ZETA*ZETA_01_p2*RS_10*ZETA_10 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_3ZETA*ZETA_01_p2*RS_10*ZETA_10 + 
                            8.0*iF0*F*ZETA_p3*minA_2RS*RS_02*RS_10*ZETA_10 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_02*RS_10*ZETA_10 + 
                            24.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*ZETA_02*RS_10*ZETA_10 + 
                            16.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*ZETA_02*RS_10*ZETA_10 + 
                            24.0*iF0*F*ZETA_p2*minA_1RS*ZETA_02*RS_10*ZETA_10 + 
                            16.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*ZETA_02*RS_10*ZETA_10 + 
                            2.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*ZETA_02*RS_10*ZETA_10 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*ZETA_02*RS_10*ZETA_10 + 
                            12.0*iF0*F*ZETA_p2*minA_2RS*RS_01_p2*ZETA_10_p2 + 
                            8.0*iF0*ZETA_p3*F_1ZETA*minA_2RS*RS_01_p2*ZETA_10_p2 + 
                            iF0*(-1.0 + ZETA_p4)*F_2ZETA*minA_2RS*RS_01_p2*ZETA_10_p2 + 
                            48.0*F*ZETA*(-E0_1RS + E1_1RS)*RS_01*ZETA_01*ZETA_10_p2 + 
                            72.0*ZETA_p2*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_01*ZETA_10_p2 + 
                            48.0*iF0*F*ZETA*minA_1RS*RS_01*ZETA_01*ZETA_10_p2 + 
                            72.0*iF0*ZETA_p2*F_1ZETA*minA_1RS*RS_01*ZETA_01*ZETA_10_p2 + 
                            24.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_2ZETA*RS_01*ZETA_01*ZETA_10_p2 + 
                            24.0*iF0*ZETA_p3*minA_1RS*F_2ZETA*RS_01*ZETA_01*ZETA_10_p2 + 
                            2.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_3ZETA*RS_01*ZETA_01*ZETA_10_p2 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_3ZETA*RS_01*ZETA_01*ZETA_10_p2 + 
                            24.0*(-E0 + E1)*F*ZETA_01_p2*ZETA_10_p2 + 
                            24.0*iF0*F*minA*ZETA_01_p2*ZETA_10_p2 + 
                            96.0*(-E0 + E1)*ZETA*F_1ZETA*ZETA_01_p2*ZETA_10_p2 + 
                            96.0*iF0*ZETA*minA*F_1ZETA*ZETA_01_p2*ZETA_10_p2 + 
                            72.0*(-E0 + E1)*ZETA_p2*F_2ZETA*ZETA_01_p2*ZETA_10_p2 + 
                            72.0*iF0*ZETA_p2*minA*F_2ZETA*ZETA_01_p2*ZETA_10_p2 + 
                            16.0*(-E0 + E1)*ZETA_p3*F_3ZETA*ZETA_01_p2*ZETA_10_p2 + 
                            16.0*iF0*ZETA_p3*minA*F_3ZETA*ZETA_01_p2*ZETA_10_p2 + 
                            (-E0 + E1)*ZETA_p4*F_4ZETA*ZETA_01_p2*ZETA_10_p2 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_4ZETA*ZETA_01_p2*ZETA_10_p2 + 
                            12.0*iF0*F*ZETA_p2*minA_1RS*RS_02*ZETA_10_p2 + 
                            8.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*RS_02*ZETA_10_p2 + 
                            iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*RS_02*ZETA_10_p2 + 
                            12.0*F*ZETA_p2*((-E0_2RS + E1_2RS)*RS_01_p2 + 
                                            (-E0_1RS + E1_1RS)*RS_02)*ZETA_10_p2 + 
                            8.0*ZETA_p3*F_1ZETA*((-E0_2RS + E1_2RS)*RS_01_p2 + 
                                                 (-E0_1RS + E1_1RS)*RS_02)*ZETA_10_p2 + 
                            ZETA_p4*F_2ZETA*((-E0_2RS + E1_2RS)*RS_01_p2 + 
                                             (-E0_1RS + E1_1RS)*RS_02)*ZETA_10_p2 + 
                            24.0*(-E0 + E1)*F*ZETA*ZETA_02*ZETA_10_p2 + 
                            24.0*iF0*F*ZETA*minA*ZETA_02*ZETA_10_p2 + 
                            36.0*(-E0 + E1)*ZETA_p2*F_1ZETA*ZETA_02*ZETA_10_p2 + 
                            36.0*iF0*ZETA_p2*minA*F_1ZETA*ZETA_02*ZETA_10_p2 + 
                            12.0*(-E0 + E1)*ZETA_p3*F_2ZETA*ZETA_02*ZETA_10_p2 + 
                            12.0*iF0*ZETA_p3*minA*F_2ZETA*ZETA_02*ZETA_10_p2 + 
                            (-E0 + E1)*ZETA_p4*F_3ZETA*ZETA_02*ZETA_10_p2 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_3ZETA*ZETA_02*ZETA_10_p2 + 
                            4.0*E0_3RS*RS_01*RS_10*RS_11 + 
                            4.0*iF0*F*(-1.0 + ZETA_p4)*minA_3RS*RS_01*RS_10*RS_11 + 
                            16.0*iF0*F*ZETA_p3*minA_2RS*ZETA_01*RS_10*RS_11 + 
                            4.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*ZETA_01*RS_10*RS_11 + 
                            16.0*iF0*F*ZETA_p3*minA_2RS*RS_01*ZETA_10*RS_11 + 
                            4.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_01*ZETA_10*RS_11 + 
                            48.0*iF0*F*ZETA_p2*minA_1RS*ZETA_01*ZETA_10*RS_11 + 
                            32.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*ZETA_01*ZETA_10*RS_11 + 
                            4.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*ZETA_01*ZETA_10*RS_11 + 
                            2.0*E0_2RS*RS_11_p2 + 
                            2.0*iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_11_p2 + 
                            48.0*F*ZETA_p2*ZETA_01*ZETA_10*((-E0_2RS + E1_2RS)*RS_01*RS_10 + 
                                                            (-E0_1RS + E1_1RS)*RS_11) + 
                            32.0*ZETA_p3*F_1ZETA*ZETA_01*ZETA_10*((-E0_2RS + E1_2RS)*RS_01*RS_10 + 
                                                                  (-E0_1RS + E1_1RS)*RS_11) + 
                            4.0*ZETA_p4*F_2ZETA*ZETA_01*ZETA_10*((-E0_2RS + E1_2RS)*RS_01*RS_10 + 
                                                                 (-E0_1RS + E1_1RS)*RS_11) + 
                            16.0*iF0*F*ZETA_p3*minA_2RS*RS_01*RS_10*ZETA_11 + 
                            4.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_01*RS_10*ZETA_11 + 
                            48.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*ZETA_01*RS_10*ZETA_11 + 
                            32.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*ZETA_01*RS_10*ZETA_11 + 
                            48.0*iF0*F*ZETA_p2*minA_1RS*ZETA_01*RS_10*ZETA_11 + 
                            32.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*ZETA_01*RS_10*ZETA_11 + 
                            4.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*ZETA_01*RS_10*ZETA_11 + 
                            4.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*ZETA_01*RS_10*ZETA_11 + 
                            48.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*RS_01*ZETA_10*ZETA_11 + 
                            32.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_10*ZETA_11 + 
                            48.0*iF0*F*ZETA_p2*minA_1RS*RS_01*ZETA_10*ZETA_11 + 
                            32.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*RS_01*ZETA_10*ZETA_11 + 
                            4.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*RS_01*ZETA_10*ZETA_11 + 
                            4.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*RS_01*ZETA_10*ZETA_11 + 
                            96.0*(-E0 + E1)*F*ZETA*ZETA_01*ZETA_10*ZETA_11 + 
                            96.0*iF0*F*ZETA*minA*ZETA_01*ZETA_10*ZETA_11 + 
                            144.0*(-E0 + E1)*ZETA_p2*F_1ZETA*ZETA_01*ZETA_10*ZETA_11 + 
                            144.0*iF0*ZETA_p2*minA*F_1ZETA*ZETA_01*ZETA_10*ZETA_11 + 
                            48.0*(-E0 + E1)*ZETA_p3*F_2ZETA*ZETA_01*ZETA_10*ZETA_11 + 
                            48.0*iF0*ZETA_p3*minA*F_2ZETA*ZETA_01*ZETA_10*ZETA_11 + 
                            4.0*(-E0 + E1)*ZETA_p4*F_3ZETA*ZETA_01*ZETA_10*ZETA_11 + 
                            4.0*iF0*(-1.0 + ZETA_p4)*minA*F_3ZETA*ZETA_01*ZETA_10*ZETA_11 + 
                            16.0*iF0*F*ZETA_p3*minA_1RS*RS_11*ZETA_11 + 
                            4.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_11*ZETA_11 + 
                            16.0*F*ZETA_p3*((-E0_2RS + E1_2RS)*RS_01*RS_10 + 
                                            (-E0_1RS + E1_1RS)*RS_11)*ZETA_11 + 
                            4.0*ZETA_p4*F_1ZETA*((-E0_2RS + E1_2RS)*RS_01*RS_10 + 
                                                 (-E0_1RS + E1_1RS)*RS_11)*ZETA_11 + 
                            24.0*(-E0 + E1)*F*ZETA_p2*ZETA_11_p2 + 
                            24.0*iF0*F*ZETA_p2*minA*ZETA_11_p2 + 
                            16.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_11_p2 + 
                            16.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_11_p2 + 
                            2.0*(-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_11_p2 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_11_p2 + 
                            2.0*E0_2RS*RS_10*RS_12 + 
                            2.0*iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_10*RS_12 + 
                            8.0*iF0*F*ZETA_p3*minA_1RS*ZETA_10*RS_12 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*ZETA_10*RS_12 + 
                            8.0*F*ZETA_p3*ZETA_10*(-(((E0_3RS-E1_3RS)*RS_01_p2 + 
                                                      (E0_2RS-E1_2RS)*RS_02)*RS_10)-2.0*(E0_2RS-E1_2RS)*RS_01*RS_11 + 
                                                   (-E0_1RS + E1_1RS)*RS_12) + 
                            2.0*ZETA_p4*F_1ZETA*ZETA_10*(-(((E0_3RS-E1_3RS)*RS_01_p2 + 
                                                            (E0_2RS-E1_2RS)*RS_02)*RS_10)-2.0*(E0_2RS-E1_2RS)*RS_01*RS_11 + 
                                                         (-E0_1RS + E1_1RS)*RS_12) + 
                            8.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*RS_10*ZETA_12 + 
                            2.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*RS_10*ZETA_12 + 
                            8.0*iF0*F*ZETA_p3*minA_1RS*RS_10*ZETA_12 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_10*ZETA_12 + 
                            24.0*(-E0 + E1)*F*ZETA_p2*ZETA_10*ZETA_12 + 
                            24.0*iF0*F*ZETA_p2*minA*ZETA_10*ZETA_12 + 
                            16.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_10*ZETA_12 + 
                            16.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_10*ZETA_12 + 
                            2.0*(-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_10*ZETA_12 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_10*ZETA_12 + 
                            E0_3RS*RS_01_p2*RS_20 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_3RS*RS_01_p2*RS_20 + 
                            8.0*iF0*F*ZETA_p3*minA_2RS*RS_01*ZETA_01*RS_20 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_01*ZETA_01*RS_20 + 
                            12.0*iF0*F*ZETA_p2*minA_1RS*ZETA_01_p2*RS_20 + 
                            8.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*ZETA_01_p2*RS_20 + 
                            iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*ZETA_01_p2*RS_20 + 
                            E0_2RS*RS_02*RS_20 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_02*RS_20 + 
                            4.0*iF0*F*ZETA_p3*minA_1RS*ZETA_02*RS_20 + 
                            iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*ZETA_02*RS_20 + 
                            12.0*F*ZETA_p2*ZETA_01_p2*((-E0_2RS + E1_2RS)*RS_10_p2 + 
                                                       (-E0_1RS + E1_1RS)*RS_20) + 
                            8.0*ZETA_p3*F_1ZETA*ZETA_01_p2*((-E0_2RS + E1_2RS)*RS_10_p2 + 
                                                            (-E0_1RS + E1_1RS)*RS_20) + 
                            ZETA_p4*F_2ZETA*ZETA_01_p2*((-E0_2RS + E1_2RS)*RS_10_p2 + 
                                                        (-E0_1RS + E1_1RS)*RS_20) + 
                            4.0*F*ZETA_p3*ZETA_02*((-E0_2RS + E1_2RS)*RS_10_p2 + 
                                                   (-E0_1RS + E1_1RS)*RS_20) + 
                            ZETA_p4*F_1ZETA*ZETA_02*((-E0_2RS + E1_2RS)*RS_10_p2 + 
                                                     (-E0_1RS + E1_1RS)*RS_20) + 
                            4.0*iF0*F*ZETA_p3*minA_2RS*RS_01_p2*ZETA_20 + 
                            iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_01_p2*ZETA_20 + 
                            24.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*RS_01*ZETA_01*ZETA_20 + 
                            16.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_01*ZETA_20 + 
                            24.0*iF0*F*ZETA_p2*minA_1RS*RS_01*ZETA_01*ZETA_20 + 
                            16.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*RS_01*ZETA_01*ZETA_20 + 
                            2.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*RS_01*ZETA_01*ZETA_20 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*RS_01*ZETA_01*ZETA_20 + 
                            24.0*(-E0 + E1)*F*ZETA*ZETA_01_p2*ZETA_20 + 
                            24.0*iF0*F*ZETA*minA*ZETA_01_p2*ZETA_20 + 
                            36.0*(-E0 + E1)*ZETA_p2*F_1ZETA*ZETA_01_p2*ZETA_20 + 
                            36.0*iF0*ZETA_p2*minA*F_1ZETA*ZETA_01_p2*ZETA_20 + 
                            12.0*(-E0 + E1)*ZETA_p3*F_2ZETA*ZETA_01_p2*ZETA_20 + 
                            12.0*iF0*ZETA_p3*minA*F_2ZETA*ZETA_01_p2*ZETA_20 + 
                            (-E0 + E1)*ZETA_p4*F_3ZETA*ZETA_01_p2*ZETA_20 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_3ZETA*ZETA_01_p2*ZETA_20 + 
                            4.0*iF0*F*ZETA_p3*minA_1RS*RS_02*ZETA_20 + 
                            iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_02*ZETA_20 + 
                            4.0*F*ZETA_p3*((-E0_2RS + E1_2RS)*RS_01_p2 + 
                                           (-E0_1RS + E1_1RS)*RS_02)*ZETA_20 + 
                            ZETA_p4*F_1ZETA*((-E0_2RS + E1_2RS)*RS_01_p2 + 
                                             (-E0_1RS + E1_1RS)*RS_02)*ZETA_20 + 
                            12.0*(-E0 + E1)*F*ZETA_p2*ZETA_02*ZETA_20 + 
                            12.0*iF0*F*ZETA_p2*minA*ZETA_02*ZETA_20 + 
                            8.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_02*ZETA_20 + 
                            8.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_02*ZETA_20 + 
                            (-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_02*ZETA_20 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_02*ZETA_20 + 
                            2.0*E0_2RS*RS_01*RS_21 + 
                            2.0*iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_01*RS_21 + 
                            8.0*iF0*F*ZETA_p3*minA_1RS*ZETA_01*RS_21 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*ZETA_01*RS_21 + 
                            8.0*F*ZETA_p3*ZETA_01*(-(E0_3RS*RS_01*RS_10_p2) + 
                                                   E1_3RS*RS_01*RS_10_p2-(E0_2RS-E1_2RS)*(2.0*RS_10*RS_11 + 
                                                                                          RS_01*RS_20) + 
                                                   (-E0_1RS + E1_1RS)*RS_21) + 
                            2.0*ZETA_p4*F_1ZETA*ZETA_01*(-(E0_3RS*RS_01*RS_10_p2) + 
                                                         E1_3RS*RS_01*RS_10_p2-(E0_2RS-E1_2RS)*(2.0*RS_10*RS_11 + 
                                                                                                RS_01*RS_20) + 
                                                         (-E0_1RS + E1_1RS)*RS_21) + 
                            8.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*RS_01*ZETA_21 + 
                            2.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_21 + 
                            8.0*iF0*F*ZETA_p3*minA_1RS*RS_01*ZETA_21 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_01*ZETA_21 + 
                            24.0*(-E0 + E1)*F*ZETA_p2*ZETA_01*ZETA_21 + 
                            24.0*iF0*F*ZETA_p2*minA*ZETA_01*ZETA_21 + 
                            16.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_01*ZETA_21 + 
                            16.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_01*ZETA_21 + 
                            2.0*(-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_01*ZETA_21 + 
                            2.0*iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_01*ZETA_21 + 
                            E0_1RS*RS_22 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_1RS*RS_22 + 
                            F*ZETA_p4*(-(E0_4RS*RS_01_p2*RS_10_p2) + 
                                       E1_4RS*RS_01_p2*RS_10_p2-E0_3RS*RS_02*RS_10_p2 + 
                                       E1_3RS*RS_02*RS_10_p2-4.0*E0_3RS*RS_01*RS_10*RS_11 + 
                                       4.0*E1_3RS*RS_01*RS_10*RS_11-2.0*E0_2RS*RS_11_p2 + 
                                       2.0*E1_2RS*RS_11_p2-2.0*E0_2RS*RS_10*RS_12 + 
                                       2.0*E1_2RS*RS_10*RS_12-E0_3RS*RS_01_p2*RS_20 + 
                                       E1_3RS*RS_01_p2*RS_20-E0_2RS*RS_02*RS_20 + 
                                       E1_2RS*RS_02*RS_20-2.0*(E0_2RS-E1_2RS)*RS_01*RS_21 + 
                                       (-E0_1RS + E1_1RS)*RS_22) + 
                            4.0*(-E0 + E1)*F*ZETA_p3*ZETA_22 + 
                            4.0*iF0*F*ZETA_p3*minA*ZETA_22 + 
                            (-E0 + E1)*ZETA_p4*F_1ZETA*ZETA_22 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_1ZETA*ZETA_22
	);

    ds->df13000 += factor* (
                            E0_4RS*RS_01_p3*RS_10 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_4RS*RS_01_p3*RS_10 + 
                            12.0*iF0*F*ZETA_p3*minA_3RS*RS_01_p2*ZETA_01*RS_10 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_3RS*RS_01_p2*ZETA_01*RS_10 + 
                            36.0*iF0*F*ZETA_p2*minA_2RS*RS_01*ZETA_01_p2*RS_10 + 
                            24.0*iF0*ZETA_p3*F_1ZETA*minA_2RS*RS_01*ZETA_01_p2*RS_10 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_2ZETA*minA_2RS*RS_01*ZETA_01_p2*RS_10 + 
                            24.0*F*ZETA*(-E0_1RS + E1_1RS)*ZETA_01_p3*RS_10 + 
                            36.0*ZETA_p2*(-E0_1RS + E1_1RS)*F_1ZETA*ZETA_01_p3*RS_10 + 
                            24.0*iF0*F*ZETA*minA_1RS*ZETA_01_p3*RS_10 + 
                            36.0*iF0*ZETA_p2*F_1ZETA*minA_1RS*ZETA_01_p3*RS_10 + 
                            12.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_2ZETA*ZETA_01_p3*RS_10 + 
                            12.0*iF0*ZETA_p3*minA_1RS*F_2ZETA*ZETA_01_p3*RS_10 + 
                            ZETA_p4*(-E0_1RS + E1_1RS)*F_3ZETA*ZETA_01_p3*RS_10 + 
                            iF0*(-1.0 + ZETA_p4)*minA_1RS*F_3ZETA*ZETA_01_p3*RS_10 + 
                            3.0*E0_3RS*RS_01*RS_02*RS_10 + 
                            3.0*iF0*F*(-1.0 + ZETA_p4)*minA_3RS*RS_01*RS_02*RS_10 + 
                            12.0*iF0*F*ZETA_p3*minA_2RS*ZETA_01*RS_02*RS_10 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*ZETA_01*RS_02*RS_10 + 
                            12.0*iF0*F*ZETA_p3*minA_2RS*RS_01*ZETA_02*RS_10 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_01*ZETA_02*RS_10 + 
                            36.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*ZETA_01*ZETA_02*RS_10 + 
                            24.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*ZETA_01*ZETA_02*RS_10 + 
                            36.0*iF0*F*ZETA_p2*minA_1RS*ZETA_01*ZETA_02*RS_10 + 
                            24.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*ZETA_01*ZETA_02*RS_10 + 
                            3.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*ZETA_01*ZETA_02*RS_10 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*ZETA_01*ZETA_02*RS_10 + 
                            E0_2RS*RS_03*RS_10 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_03*RS_10 + 
                            4.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*ZETA_03*RS_10 + 
                            ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*ZETA_03*RS_10 + 
                            4.0*iF0*F*ZETA_p3*minA_1RS*ZETA_03*RS_10 + 
                            iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*ZETA_03*RS_10 + 
                            4.0*iF0*F*ZETA_p3*minA_3RS*RS_01_p3*ZETA_10 + 
                            iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_3RS*RS_01_p3*ZETA_10 + 
                            36.0*iF0*F*ZETA_p2*minA_2RS*RS_01_p2*ZETA_01*ZETA_10 + 
                            24.0*iF0*ZETA_p3*F_1ZETA*minA_2RS*RS_01_p2*ZETA_01*ZETA_10 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_2ZETA*minA_2RS*RS_01_p2*ZETA_01*ZETA_10 + 
                            72.0*F*ZETA*(-E0_1RS + E1_1RS)*RS_01*ZETA_01_p2*ZETA_10 + 
                            108.0*ZETA_p2*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_01_p2*ZETA_10 + 
                            72.0*iF0*F*ZETA*minA_1RS*RS_01*ZETA_01_p2*ZETA_10 + 
                            108.0*iF0*ZETA_p2*F_1ZETA*minA_1RS*RS_01*ZETA_01_p2*ZETA_10 + 
                            36.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_2ZETA*RS_01*ZETA_01_p2*ZETA_10 + 
                            36.0*iF0*ZETA_p3*minA_1RS*F_2ZETA*RS_01*ZETA_01_p2*ZETA_10 + 
                            3.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_3ZETA*RS_01*ZETA_01_p2*ZETA_10 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_3ZETA*RS_01*ZETA_01_p2*ZETA_10 + 
                            24.0*(-E0 + E1)*F*ZETA_01_p3*ZETA_10 + 
                            24.0*iF0*F*minA*ZETA_01_p3*ZETA_10 + 
                            96.0*(-E0 + E1)*ZETA*F_1ZETA*ZETA_01_p3*ZETA_10 + 
                            96.0*iF0*ZETA*minA*F_1ZETA*ZETA_01_p3*ZETA_10 + 
                            72.0*(-E0 + E1)*ZETA_p2*F_2ZETA*ZETA_01_p3*ZETA_10 + 
                            72.0*iF0*ZETA_p2*minA*F_2ZETA*ZETA_01_p3*ZETA_10 + 
                            16.0*(-E0 + E1)*ZETA_p3*F_3ZETA*ZETA_01_p3*ZETA_10 + 
                            16.0*iF0*ZETA_p3*minA*F_3ZETA*ZETA_01_p3*ZETA_10 + 
                            (-E0 + E1)*ZETA_p4*F_4ZETA*ZETA_01_p3*ZETA_10 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_4ZETA*ZETA_01_p3*ZETA_10 + 
                            12.0*iF0*F*ZETA_p3*minA_2RS*RS_01*RS_02*ZETA_10 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_01*RS_02*ZETA_10 + 
                            36.0*iF0*F*ZETA_p2*minA_1RS*ZETA_01*RS_02*ZETA_10 + 
                            24.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*ZETA_01*RS_02*ZETA_10 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*ZETA_01*RS_02*ZETA_10 + 
                            36.0*F*ZETA_p2*ZETA_01*((-E0_2RS + E1_2RS)*RS_01_p2 + 
                                                    (-E0_1RS + E1_1RS)*RS_02)*ZETA_10 + 
                            24.0*ZETA_p3*F_1ZETA*ZETA_01*((-E0_2RS + E1_2RS)*RS_01_p2 + 
                                                          (-E0_1RS + E1_1RS)*RS_02)*ZETA_10 + 
                            3.0*ZETA_p4*F_2ZETA*ZETA_01*((-E0_2RS + E1_2RS)*RS_01_p2 + 
                                                         (-E0_1RS + E1_1RS)*RS_02)*ZETA_10 + 
                            36.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*RS_01*ZETA_02*ZETA_10 + 
                            24.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_02*ZETA_10 + 
                            36.0*iF0*F*ZETA_p2*minA_1RS*RS_01*ZETA_02*ZETA_10 + 
                            24.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*RS_01*ZETA_02*ZETA_10 + 
                            3.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*RS_01*ZETA_02*ZETA_10 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*RS_01*ZETA_02*ZETA_10 + 
                            72.0*(-E0 + E1)*F*ZETA*ZETA_01*ZETA_02*ZETA_10 + 
                            72.0*iF0*F*ZETA*minA*ZETA_01*ZETA_02*ZETA_10 + 
                            108.0*(-E0 + E1)*ZETA_p2*F_1ZETA*ZETA_01*ZETA_02*ZETA_10 + 
                            108.0*iF0*ZETA_p2*minA*F_1ZETA*ZETA_01*ZETA_02*ZETA_10 + 
                            36.0*(-E0 + E1)*ZETA_p3*F_2ZETA*ZETA_01*ZETA_02*ZETA_10 + 
                            36.0*iF0*ZETA_p3*minA*F_2ZETA*ZETA_01*ZETA_02*ZETA_10 + 
                            3.0*(-E0 + E1)*ZETA_p4*F_3ZETA*ZETA_01*ZETA_02*ZETA_10 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA*F_3ZETA*ZETA_01*ZETA_02*ZETA_10 + 
                            4.0*iF0*F*ZETA_p3*minA_1RS*RS_03*ZETA_10 + 
                            iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_03*ZETA_10 + 
                            4.0*F*ZETA_p3*(-((E0_3RS-E1_3RS)*RS_01_p3)-3.0*(E0_2RS-E1_2RS)*RS_01*RS_02 + 
                                           (-E0_1RS + E1_1RS)*RS_03)*ZETA_10 + 
                            ZETA_p4*F_1ZETA*(-((E0_3RS-E1_3RS)*RS_01_p3)-3.0*(E0_2RS-E1_2RS)*RS_01*RS_02 + 
                                             (-E0_1RS + E1_1RS)*RS_03)*ZETA_10 + 
                            12.0*(-E0 + E1)*F*ZETA_p2*ZETA_03*ZETA_10 + 
                            12.0*iF0*F*ZETA_p2*minA*ZETA_03*ZETA_10 + 
                            8.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_03*ZETA_10 + 
                            8.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_03*ZETA_10 + 
                            (-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_03*ZETA_10 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_03*ZETA_10 + 
                            3.0*E0_3RS*RS_01_p2*RS_11 + 
                            3.0*iF0*F*(-1.0 + ZETA_p4)*minA_3RS*RS_01_p2*RS_11 + 
                            24.0*iF0*F*ZETA_p3*minA_2RS*RS_01*ZETA_01*RS_11 + 
                            6.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_01*ZETA_01*RS_11 + 
                            36.0*iF0*F*ZETA_p2*minA_1RS*ZETA_01_p2*RS_11 + 
                            24.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*ZETA_01_p2*RS_11 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*ZETA_01_p2*RS_11 + 
                            3.0*E0_2RS*RS_02*RS_11 + 
                            3.0*iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_02*RS_11 + 
                            12.0*iF0*F*ZETA_p3*minA_1RS*ZETA_02*RS_11 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*ZETA_02*RS_11 + 
                            36.0*F*ZETA_p2*ZETA_01_p2*((-E0_2RS + E1_2RS)*RS_01*RS_10 + 
                                                       (-E0_1RS + E1_1RS)*RS_11) + 
                            24.0*ZETA_p3*F_1ZETA*ZETA_01_p2*((-E0_2RS + E1_2RS)*RS_01*RS_10 + 
                                                             (-E0_1RS + E1_1RS)*RS_11) + 
                            3.0*ZETA_p4*F_2ZETA*ZETA_01_p2*((-E0_2RS + E1_2RS)*RS_01*RS_10 + 
                                                            (-E0_1RS + E1_1RS)*RS_11) + 
                            12.0*F*ZETA_p3*ZETA_02*((-E0_2RS + E1_2RS)*RS_01*RS_10 + 
                                                    (-E0_1RS + E1_1RS)*RS_11) + 
                            3.0*ZETA_p4*F_1ZETA*ZETA_02*((-E0_2RS + E1_2RS)*RS_01*RS_10 + 
                                                         (-E0_1RS + E1_1RS)*RS_11) + 
                            12.0*iF0*F*ZETA_p3*minA_2RS*RS_01_p2*ZETA_11 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_2RS*RS_01_p2*ZETA_11 + 
                            72.0*F*ZETA_p2*(-E0_1RS + E1_1RS)*RS_01*ZETA_01*ZETA_11 + 
                            48.0*ZETA_p3*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_01*ZETA_11 + 
                            72.0*iF0*F*ZETA_p2*minA_1RS*RS_01*ZETA_01*ZETA_11 + 
                            48.0*iF0*ZETA_p3*F_1ZETA*minA_1RS*RS_01*ZETA_01*ZETA_11 + 
                            6.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_2ZETA*RS_01*ZETA_01*ZETA_11 + 
                            6.0*iF0*(-1.0 + ZETA_p4)*minA_1RS*F_2ZETA*RS_01*ZETA_01*ZETA_11 + 
                            72.0*(-E0 + E1)*F*ZETA*ZETA_01_p2*ZETA_11 + 
                            72.0*iF0*F*ZETA*minA*ZETA_01_p2*ZETA_11 + 
                            108.0*(-E0 + E1)*ZETA_p2*F_1ZETA*ZETA_01_p2*ZETA_11 + 
                            108.0*iF0*ZETA_p2*minA*F_1ZETA*ZETA_01_p2*ZETA_11 + 
                            36.0*(-E0 + E1)*ZETA_p3*F_2ZETA*ZETA_01_p2*ZETA_11 + 
                            36.0*iF0*ZETA_p3*minA*F_2ZETA*ZETA_01_p2*ZETA_11 + 
                            3.0*(-E0 + E1)*ZETA_p4*F_3ZETA*ZETA_01_p2*ZETA_11 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA*F_3ZETA*ZETA_01_p2*ZETA_11 + 
                            12.0*iF0*F*ZETA_p3*minA_1RS*RS_02*ZETA_11 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_02*ZETA_11 + 
                            12.0*F*ZETA_p3*((-E0_2RS + E1_2RS)*RS_01_p2 + 
                                            (-E0_1RS + E1_1RS)*RS_02)*ZETA_11 + 
                            3.0*ZETA_p4*F_1ZETA*((-E0_2RS + E1_2RS)*RS_01_p2 + 
                                                 (-E0_1RS + E1_1RS)*RS_02)*ZETA_11 + 
                            36.0*(-E0 + E1)*F*ZETA_p2*ZETA_02*ZETA_11 + 
                            36.0*iF0*F*ZETA_p2*minA*ZETA_02*ZETA_11 + 
                            24.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_02*ZETA_11 + 
                            24.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_02*ZETA_11 + 
                            3.0*(-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_02*ZETA_11 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_02*ZETA_11 + 
                            3.0*E0_2RS*RS_01*RS_12 + 
                            3.0*iF0*F*(-1.0 + ZETA_p4)*minA_2RS*RS_01*RS_12 + 
                            12.0*iF0*F*ZETA_p3*minA_1RS*ZETA_01*RS_12 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*ZETA_01*RS_12 + 
                            12.0*F*ZETA_p3*ZETA_01*(-(((E0_3RS-E1_3RS)*RS_01_p2 + 
                                                       (E0_2RS-E1_2RS)*RS_02)*RS_10)-2.0*(E0_2RS-E1_2RS)*RS_01*RS_11 + 
                                                    (-E0_1RS + E1_1RS)*RS_12) + 
                            3.0*ZETA_p4*F_1ZETA*ZETA_01*(-(((E0_3RS-E1_3RS)*RS_01_p2 + 
                                                            (E0_2RS-E1_2RS)*RS_02)*RS_10)-2.0*(E0_2RS-E1_2RS)*RS_01*RS_11 + 
                                                         (-E0_1RS + E1_1RS)*RS_12) + 
                            12.0*F*ZETA_p3*(-E0_1RS + E1_1RS)*RS_01*ZETA_12 + 
                            3.0*ZETA_p4*(-E0_1RS + E1_1RS)*F_1ZETA*RS_01*ZETA_12 + 
                            12.0*iF0*F*ZETA_p3*minA_1RS*RS_01*ZETA_12 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*F_1ZETA*minA_1RS*RS_01*ZETA_12 + 
                            36.0*(-E0 + E1)*F*ZETA_p2*ZETA_01*ZETA_12 + 
                            36.0*iF0*F*ZETA_p2*minA*ZETA_01*ZETA_12 + 
                            24.0*(-E0 + E1)*ZETA_p3*F_1ZETA*ZETA_01*ZETA_12 + 
                            24.0*iF0*ZETA_p3*minA*F_1ZETA*ZETA_01*ZETA_12 + 
                            3.0*(-E0 + E1)*ZETA_p4*F_2ZETA*ZETA_01*ZETA_12 + 
                            3.0*iF0*(-1.0 + ZETA_p4)*minA*F_2ZETA*ZETA_01*ZETA_12 + 
                            E0_1RS*RS_40 + 
                            iF0*F*(-1.0 + ZETA_p4)*minA_1RS*RS_40 + 
                            F*ZETA_p4*(-(E0_4RS*RS_01_p3*RS_10) + 
                                       E1_4RS*RS_01_p3*RS_10-3.0*E0_3RS*RS_01*RS_02*RS_10 + 
                                       3.0*E1_3RS*RS_01*RS_02*RS_10-E0_2RS*RS_03*RS_10 + 
                                       E1_2RS*RS_03*RS_10-3.0*E0_3RS*RS_01_p2*RS_11 + 
                                       3.0*E1_3RS*RS_01_p2*RS_11-3.0*E0_2RS*RS_02*RS_11 + 
                                       3.0*E1_2RS*RS_02*RS_11-3.0*(E0_2RS-E1_2RS)*RS_01*RS_12 + 
                                       (-E0_1RS + E1_1RS)*RS_40) + 
                            4.0*(-E0 + E1)*F*ZETA_p3*ZETA_13 + 
                            4.0*iF0*F*ZETA_p3*minA*ZETA_13 + 
                            (-E0 + E1)*ZETA_p4*F_1ZETA*ZETA_13 + 
                            iF0*(-1.0 + ZETA_p4)*minA*F_1ZETA*ZETA_13
	);

}
#endif

/* Auxaliary functions specific for pw92 */


real Gfit(real *pRS, const GfitParams *COEF)
{
    /* Declarations */
    real  A;
    real  G1;
    real  G2;
    real  RS_p1f2;
    real  RS_p3f2;
    real  alpha1;
    real  beta1;
    real  beta2;
    real  beta3;
    real  beta4;
    real  p;

    real  RS;

    /* Setting up G parameters */	
    A=COEF->A;
    alpha1=COEF->alpha1;
    beta1=COEF->beta1;
    beta2=COEF->beta2;
    beta3=COEF->beta3;
    beta4=COEF->beta4;
    p=COEF->p;
	
    /* Powers of RS */
    RS=*pRS;
    RS_p1f2=sqrt(RS);
    RS_p3f2=RS_p1f2*RS;

    /* Aux functions G */
    G1=-2.0*A*(1.0+alpha1*RS);
    G2=2.0*A*(beta1*RS_p1f2 + beta2*RS + beta3*RS_p3f2 + beta4*pow(RS,(p+1.0)));

    /*G fit 
     *  G=G1*log(1.0+(1.0/G2)) */
    return ( G1*log(1.0+(1.0/G2)) );
}


real Gfit_1RS(real *pRS, const GfitParams *COEF) {
	
    /* Declarations */
    real  A;
    real  G1;
    real  G2;
    real  G1_1RS, G2_1RS;
    real  RS_p1f2;
    real  RS_p3f2;
    real  RS_p1f2p;
    real  RS_pp;
    real  alpha1;
    real  beta1;
    real  beta2;
    real  beta3;
    real  beta4;
    real  p;

    real  RS;

    /* Setting up G parameters */	
    A=COEF->A;
    alpha1=COEF->alpha1;
    beta1=COEF->beta1;
    beta2=COEF->beta2;
    beta3=COEF->beta3;
    beta4=COEF->beta4;
    p=COEF->p;
	
    /* Powers of RS */
    RS=*pRS;
    RS_p1f2=sqrt(RS);
    RS_pp=pow(RS,p);
    RS_p3f2=RS_p1f2*RS;
    RS_p1f2p=RS_pp*RS_p1f2; /* RS_p1f2p=pow(RS,(0.5+p));  */

    /* Aux functions G */
    G1=-2.0*A*(1.0+alpha1*RS);
    G2=2.0*A*(beta1*RS_p1f2 + beta2*RS + beta3*RS_p3f2 + beta4*RS_pp*RS);

    /* Derivatives of G */
    G1_1RS=-2.0*A*alpha1;
    G2_1RS=2.0*A*(beta2 + (beta1/(2.0*RS_p1f2)) + (3.0/2.0)*beta3*RS_p1f2 + beta4*(1.0+p)*RS_pp);
	
    /* G fit derivative vith respect to RS */
    return( log(1.0+(1.0/G2))*G1_1RS - (G1*G2_1RS/(G2*(1.0 + G2))) );
}


real
Gfit_2RS(real *pRS, const GfitParams *COEF)
{

    /* Declarations */
    real  A;
    real  G1;
    real  G2;
    real  G1_1RS, G2_1RS, G2_2RS;
    real  RS_p1f2;
    real  RS_p3f2;
    real  RS_pp, RS_p1f2p;
    real  G2_p2, G2_p3, G2_1RS_p2;
    real  alpha1;
    real  beta1;
    real  beta2;
    real  beta3;
    real  beta4;
    real  p;

    real  RS;

    /* Setting up G parameters */	
    A=COEF->A;
    alpha1=COEF->alpha1;
    beta1=COEF->beta1;
    beta2=COEF->beta2;
    beta3=COEF->beta3;
    beta4=COEF->beta4;
    p=COEF->p;
	
    /* Powers of RS */
    RS=*pRS;
    RS_p1f2=sqrt(RS);
    RS_pp=pow(RS,p);
    RS_p3f2=RS_p1f2*RS;
    RS_p1f2p=RS_pp*RS_p1f2; /* RS_p1f2p=RS^(0.5+p); */

    /* Aux functions G */
    G1=-2.0*A*(1.0+alpha1*RS);
    G2=2.0*A*(beta1*RS_p1f2 + beta2*RS + beta3*RS_p3f2 + beta4*RS*RS_pp);

    /* Derivatives of G */
    G1_1RS=-2.0*A*alpha1;
	
    G2_1RS=2.0*A*(beta2 + (beta1/(2.0*RS_p1f2)) + (3.0/2.0)*beta3*RS_p1f2 + beta4*(1.0+p)*RS_pp);
    G2_2RS=2.0*A*(-(beta1/(4.0*RS_p3f2)) + ((3.0/4.0)*beta3/RS_p1f2) + beta4*p*(1.0+p)*(RS_pp/RS));
	
    /* Powers of G1, G2, G3 functions */
    G2_p2=G2*G2;
    G2_p3=G2_p2*G2;
    G2_1RS_p2=G2_1RS*G2_1RS;
	

    /* G fit second derivative vith respect to RS */
    return ( 
             -(G1*G2_1RS_p2/(G2_p2*(1.0+G2)*(1.0+G2))) + 
             (2.0*G1*G2_1RS_p2/(G2_p2 + G2_p3)) - 
             (2.0*G1_1RS*G2_1RS + G1*G2_2RS)/(G2*(1.0+G2)) 
	);
}


real
Gfit_3RS(real *pRS, const GfitParams *COEF)
{

    /* Declarations */
    real  A;
    real  G1;
    real  G1_1RS;
    real  G2;
    real  G2_1RS;
    real  G2_1RS_p2;
    real  G2_1RS_p3;
    real  G2_2RS;
    real  G2_3RS;
    real  G2_p2;
    real  G2_p3;
    real  RS_p1f2;
    real  RS_p1f2p;
    real  RS_p2;
    real  RS_p3f2;
    real  RS_p5f2;
    real  RS_pp;
    real  alpha1;
    real  beta1;
    real  beta2;
    real  beta3;
    real  beta4;
    real  p;
    real  SU1G2, SU1G2_p2, SU1G2_p3;

    real  RS;

    /* Setting up G parameters */	
    A=COEF->A;
    alpha1=COEF->alpha1;
    beta1=COEF->beta1;
    beta2=COEF->beta2;
    beta3=COEF->beta3;
    beta4=COEF->beta4;
    p=COEF->p;
	
    /* Powers of RS */
    RS=*pRS;
    RS_p1f2=sqrt(RS);
    RS_pp=pow(RS,p);
    RS_p3f2=RS_p1f2*RS;
    RS_p5f2=RS_p3f2*RS;
    RS_p1f2p=RS_p1f2*RS_pp; /* RS_p1f2p=RS^(0.5+p); */
    RS_p2=RS*RS;

    G1=-2.0*A*(1.0+alpha1*RS);
    G2=2.0*A*(beta1*RS_p1f2 + beta2*RS + beta3*RS_p3f2 + beta4*RS*RS_pp);

    G1_1RS=-2.0*A*alpha1;
	
    G2_1RS=2.0*A*(beta2 + (beta1/(2.0*RS_p1f2)) + (3.0/2.0)*beta3*RS_p1f2 +
                  beta4*(1.0+p)*RS_pp);
    G2_2RS=2.0*A*(-(beta1/(4.0*RS_p3f2)) + ((3.0/4.0)*beta3/RS_p1f2) + 
                  beta4*p*(1.0+p)*(RS_pp/RS));
    G2_3RS=2.0*A*( ((3.0/8.0)*beta1/RS_p5f2) - ((3.0/8.0)*beta3/RS_p3f2) +
                   beta4*(-1.0+p)*p*(1.0+p)*(RS_pp/RS_p2));

    /* Powers of G1, G2, G3 functions */
    G2_p2=G2*G2;
    G2_p3=G2_p2*G2;
    G2_1RS_p2=G2_1RS*G2_1RS;
    G2_1RS_p3=G2_1RS_p2*G2_1RS;

    /* Often used term (and powers)*/
    SU1G2 = 1.0 + G2;
    SU1G2_p2=SU1G2*SU1G2;
    SU1G2_p3=SU1G2_p2*SU1G2;

    /* G fit third derivative vith respect to RS */
    return( 
            (6.0*G1_1RS*G2_1RS_p2/(G2_p2 + G2_p3)) - 
            (2.0*G1*(1.0 + 3.0*G2 +3.0*G2_p2)*G2_1RS_p3/(G2_p3* SU1G2_p3  )) + 
            (6.0*G1*G2_1RS*G2_2RS/(G2_p2 + G2_p3)) - 
            (3.0*G2_1RS*(G1_1RS*G2_1RS + G1*G2_2RS)/(  G2_p2*SU1G2_p2   ) ) - 
            ((3.0*G1_1RS*G2_2RS + G1*G2_3RS)/( G2*SU1G2 )) 
	);
}


real
Gfit_4RS(real *pRS, const GfitParams *COEF)
{
	
    /* Declarations */
    real  A;
    real  G1;
    real  G1_1RS;
    real  G2;
    real  G2_1RS;
    real  G2_1RS_p2;
    real  G2_1RS_p3;
    real  G2_1RS_p4;
    real  G2_2RS;
    real  G2_2RS_p2;
    real  G2_3RS;
    real  G2_4RS;
    real  G2_p2;
    real  G2_p3;
    real  G2_p4;
    real  G2_p5;
    real  G2_p6;
    real  RS_p1f2;
    real  RS_p1f2p;
    real  RS_p2;
    real  RS_p3;
    real  RS_p3f2;
    real  RS_p5f2;
    real  RS_p7f2;
    real  RS_pp;
    real  alpha1;
    real  beta1;
    real  beta2;
    real  beta3;
    real  beta4;
    real  p;
    real  SU1G2, SU1G2_p2, SU1G2_p4;

    real  RS;

    /* Setting up G parameters */
    A=COEF->A;
    alpha1=COEF->alpha1;
    beta1=COEF->beta1;
    beta2=COEF->beta2;
    beta3=COEF->beta3;
    beta4=COEF->beta4;
    p=COEF->p;

    /* Powers of RS */
    RS=*pRS;
    RS_p1f2=sqrt(RS);
    RS_pp=pow(RS,p);
    RS_p3f2=RS_p1f2*RS;
    RS_p5f2=RS_p3f2*RS;
    RS_p7f2=RS_p5f2*RS;
    RS_p1f2p=RS_p1f2*RS_pp; /* RS_p1f2p=RS^(0.5.0+p); */
    RS_p2=RS*RS;
    RS_p3=RS_p2*RS;

    /* Aux functions G */
    G1=-2.0*A*(1.0+alpha1*RS);
    G2=2.0*A*(beta1*RS_p1f2 + beta2*RS + beta3*RS_p3f2 + beta4*RS*RS_pp);

    /* Derivatives of G */
    G1_1RS=-2.0*A*alpha1;
	
    G2_1RS=2.0*A*(beta2 + (beta1/(2.0*RS_p1f2)) + (3.0/2.0)*beta3*RS_p1f2 + beta4*(1.0+p)*RS_pp);
    G2_2RS=2.0*A*(-(beta1/(4.0*RS_p3f2)) + ((3.0/4.0)*beta3/RS_p1f2) + beta4*p*(1.0+p)*(RS_pp/RS));
    G2_3RS=2.0*A*( ((3.0/8.0)*beta1/RS_p5f2) - ((3.0/8.0)*beta3/RS_p3f2) + beta4*(-1.0+p)*p*(1.0+p)*(RS_pp/RS_p2));
    G2_4RS=2.0*A*(-((15.0/16.0)*beta1/RS_p7f2) + ((9.0/16.0)*beta3/RS_p5f2) + beta4*(-2.0+p)*(-1.0+p)*p*(1.0+p)*(RS_pp/RS_p3));

    /* Powers of G1, G2 functions */
    G2_p2=G2*G2;
    G2_p3=G2_p2*G2;
    G2_p4=G2_p2*G2_p2;
    G2_p5=G2_p3*G2_p2;
    G2_p6=G2_p3*G2_p3;

    G2_1RS_p2=G2_1RS*G2_1RS;
    G2_1RS_p3=G2_1RS_p2*G2_1RS;
    G2_1RS_p4=G2_1RS_p3*G2_1RS;

    G2_2RS_p2=G2_2RS*G2_2RS;

    /* Often used term (and powers)*/
    SU1G2 = 1.0 + G2;
    SU1G2_p2=SU1G2*SU1G2;
    SU1G2_p4=SU1G2_p2*SU1G2_p2;

    /* G fit fourth derivative vith respect to RS */

    return ( 
             (6.0*G1*G2_1RS_p4-4.0*G2*G2_1RS_p2*(2.0*G1_1RS*G2_1RS+3.0*G1*(-2.0*G2_1RS_p2+G2_2RS))+G2_p2*(G1_1RS*(-32.0*G2_1RS_p3+12.0*G2_1RS*G2_2RS)+G1*(36.0*G2_1RS_p4-48.0*G2_1RS_p2*G2_2RS+3.0*G2_2RS_p2+4.0*G2_1RS*G2_3RS))+G2_p5*(12.0*G1_1RS*(2.0*G2_1RS*G2_2RS-G2_3RS)+G1*(6.0*G2_2RS_p2+8.0*G2_1RS*G2_3RS-3.0*G2_4RS))+G2_p4*(-12.0*G1_1RS*(2.0*G2_1RS_p3-5.0*G2_1RS*G2_2RS+G2_3RS)+G1*(-36.0*G2_1RS_p2*G2_2RS+15.0*G2_2RS_p2+20.0*G2_1RS*G2_3RS-3.0*G2_4RS))+G2_p3*(-4.0*G1_1RS*(12.0*G2_1RS_p3-12.0*G2_1RS*G2_2RS+G2_3RS)+G1*(24.0*G2_1RS_p4-72.0*G2_1RS_p2*G2_2RS+12.0*G2_2RS_p2+16.0*G2_1RS*G2_3RS-G2_4RS))-G2_p6*(4.0*G1_1RS*G2_3RS+G1*G2_4RS))/(G2_p4*SU1G2_p4) 
        );
	
}
