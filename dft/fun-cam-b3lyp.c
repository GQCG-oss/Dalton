/* -*- mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* fun-camb3lyp.c:

Pawel Salek, 2004.06, Himmelbjerg.
*/

 
/* strictly conform to XOPEN ANSI C standard */
#if !defined(SYS_DEC)
/* XOPEN compliance is missing on old Tru64 4.0E Alphas and pow() prototype
 * is not specified. */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1
#endif
/* Use BSD's strncasecmp() */
#define _BSD_SOURCE 1

#include <math.h>

#include <ctype.h>
#include <string.h>
#include <stdio.h>
 
#define __CVERSION__
 
#include "functionals.h"

static real CamMuFactor = 0.33, CamAlpha = 0.19, CamBeta = 0.46;
#define MU    CamMuFactor
#define ALPHA CamAlpha
#define BETA  CamBeta

#define LYP_WEIGHT 0.81
#define VWN_WEIGHT 0.19
#define ADD_CORRELATION 1

#if 1
/* yanai-consistent */
#define BECKE88_CORR_WEIGHT 1
#else
/* B3-lyp consistent */
#define BECKE88_CORR_WEIGHT 0.9
#endif

/* RGFirstDrv, RGSecondDrv and RGThirdDrv hold derivatives of a function
 * with respect to density R (R) and density gradient (G).
 * Field dfAB contains d^A/dR^a d^B/dG^B f.
 */
typedef struct {
    real df10, df01;
} RGFirstDrv;
typedef struct {
    real df10, df01, df20, df11, df02;
} RGSecondDrv;
typedef struct {
    real df10, df01, df20, df11, df02;
    real df30, df21, df12, df03;
} RGThirdDrv;
 
/* INTERFACE PART */
static int camb3lyp_isgga(void) { return 1; }
static int camb3lyp_read(const char *conf_line);
static void camb3lyp_report(void);
static real camb3lyp_energy(const DftDensProp* dp);
static void camb3lyp_first(FirstFuncDrv *ds,   real factor,
                           const DftDensProp* dp);
static void camb3lyp_second(SecondFuncDrv *ds, real factor,
                            const DftDensProp* dp);
static void camb3lyp_third(ThirdFuncDrv *ds, real factor,
                           const DftDensProp* dp);

Functional Camb3lypFunctional = {
  "Camb3lyp",       /* name */
  camb3lyp_isgga,   /* gga-corrected */
  camb3lyp_read,
  camb3lyp_report,
  camb3lyp_energy,
  camb3lyp_first,
  camb3lyp_second,
  camb3lyp_third
};
 
/* IMPLEMENTATION PART */
static int
parse_table(const char *func, const char *str,
            int cnt, const char *keywords[], float *weights)
{
  int res=1;
  while(*str) {
    int i;
    while(*str && isspace((int)*str)) str++; /* skip whitespace */
    if(*str =='\0') break; /* line ended by whitespace */
    for(i=0; i<cnt; i++) {
      int len = strlen(keywords[i]);
      if(strncasecmp(keywords[i], str, len)==0 &&
         str[len] == '=') {
        if(sscanf(str+len+1,"%g", &weights[i]) != 1) {
          fort_print("%s: %s not followed by the weight: ",
                     func, keywords[i]);
          res = 0;
        }
        break;
      }
    }
    if(i==cnt) {
      fort_print("%s: unknown string: '%s'", func, str);
      res = 0;
    }
    while(*str && !isspace((int)*str)) str++; /* skip nonws */
  }
  return res;
}

static const char *cam_keywords[] = { "alpha", "beta", "mu" };
static int
camb3lyp_read(const char *conf_line)
{
    float weights[ELEMENTS(cam_keywords)];

    weights[0] = CamAlpha;
    weights[1] = CamBeta;
    weights[2] = CamMuFactor;
    if(!parse_table("CAM-B3LYP", conf_line,
                    ELEMENTS(cam_keywords), cam_keywords, weights))
        return 0;
    dft_set_hf_weight(weights[0]);
    dft_set_cam_param(weights[2], weights[1]);
    CamMuFactor = weights[2]; CamBeta = weights[1];

    CamAlpha = weights[0];
    return 1;
}


static void
camb3lyp_report(void)
{
    fort_print("CAM-B3LYP functional with alpha=%5.3f beta=%5.3f mu=%5.3f",
               ALPHA, BETA, MU);
}

/* ===================================================================
 * compute a and its derivatives.
 * =================================================================== */
static real
fun_a(real rho, real grad)
{
  DftDensProp dp; real a, k, kd, kb;
  dp.rhoa = dp.rhob = rho;
  dp.grada = dp.gradb = grad;
  kd = DiracFunctional.func(&dp);
  kb = BeckeFunctional.func(&dp)*BECKE88_CORR_WEIGHT;
  k = -0.5*( kd + kb) * pow(rho,-4.0/3.0);
  a = MU*sqrt(k)/(6.0*sqrt(M_PI)*pow(rho,1.0/3.0));
  return a;
}

static void
fun_a_first(real rho, real grad, RGFirstDrv *res)
{
    FirstFuncDrv ds;
    DftDensProp dp;
    real ex, a;

    dp.rhoa  = dp.rhob  = rho;
    dp.grada = dp.gradb = grad;
    ex = 0.5*(DiracFunctional.func(&dp) +
              BeckeFunctional.func(&dp)*BECKE88_CORR_WEIGHT);
    a = fun_a(rho, grad);
    memset(&ds, 0, sizeof(ds));
    if(fabs(a)<1e-40 || fabs(ex)<1e-40) return;
    DiracFunctional.first(&ds, 1, &dp);
    BeckeFunctional.first(&ds, BECKE88_CORR_WEIGHT, &dp);
    res->df10 = (ds.df1000/(2*ex)-1.0/rho)*a;
    res->df01 = ds.df0010/(2*ex)*a;
}

static void
fun_a_second(real rho, real grad, RGSecondDrv *res)
{
    SecondFuncDrv ds;
    DftDensProp dp;
    real ex, a, f10, f20, f01;

    dp.rhoa  = dp.rhob  = rho;
    dp.grada = dp.gradb = grad;
    ex = 0.5*(DiracFunctional.func(&dp) +
              BeckeFunctional.func(&dp)*BECKE88_CORR_WEIGHT);
    a = fun_a(rho, grad);
    memset(&ds, 0, sizeof(ds));
    if(fabs(a)<1e-40 || fabs(ex)<1e-40) return;
    DiracFunctional.second(&ds, 1, &dp);
    BeckeFunctional.second(&ds, BECKE88_CORR_WEIGHT, &dp);
    f10 = ds.df1000/(2*ex)-1.0/rho;
    f20 = ds.df2000/(2*ex) - ds.df1000*ds.df1000/(2*ex*ex) +1/(rho*rho);
    f01 = ds.df0010/(2*ex);
    res->df10 = a*f10;
    res->df01 = a*f01;
    res->df20 = a*(f10*f10 + f20);
    res->df11 = a*(f10*f01 + ds.df1010/(2*ex) -
                   ds.df1000*ds.df0010/(2*ex*ex));
    res->df02 = a*(ds.df0020/(2*ex)-ds.df0010*ds.df0010/(4*ex*ex)); 
}

static void
fun_a_third(real rho, real grad, RGThirdDrv *res)
{
    ThirdFuncDrv ds;
    DftDensProp dp;
    real ex, a, f10, f01, f20, f11, f02, f30, f21, f12, f03;

    dp.rhoa  = dp.rhob  = rho;
    dp.grada = dp.gradb = grad;
    ex = 0.5*(DiracFunctional.func(&dp) +
              BeckeFunctional.func(&dp)*BECKE88_CORR_WEIGHT);
    a = fun_a(rho, grad);
    memset(&ds, 0, sizeof(ds));
    if(fabs(a)<1e-15 || fabs(ex)<1e-15) return;
    DiracFunctional.third(&ds, 1, &dp);
    BeckeFunctional.third(&ds, BECKE88_CORR_WEIGHT, &dp);
    f10 = ds.df1000/(2*ex)-1.0/rho;
    f01 = ds.df0010/(2*ex);
    f20 = ds.df2000/(2*ex) - ds.df1000*ds.df1000/(2*ex*ex) +1/(rho*rho);
    f11 = ds.df1010/(2*ex) - ds.df1000*ds.df0010/(2*ex*ex);
    f02 = ds.df0020/(2*ex) - ds.df0010*ds.df0010/(2*ex*ex);
    f30 = ds.df3000/(2*ex) - 1.5*ds.df1000*ds.df2000/(ex*ex)
        + ds.df1000*ds.df1000*ds.df1000/(ex*ex*ex) - 2/(rho*rho*rho);
    f21 = ds.df2010/(2*ex)
        - ds.df2000*ds.df0010/(2*ex*ex) - 2*ds.df1000*ds.df1010/(2*ex*ex)
        + ds.df1000*ds.df1000*ds.df0010/(ex*ex*ex);
    f12 = ds.df1020/(2*ex)
        - 2*ds.df1010*ds.df0010/(2*ex*ex) - ds.df1000*ds.df0020/(2*ex*ex) 
        + ds.df1000*ds.df0010*ds.df0010/(ex*ex*ex);
    f03 = ds.df0030/(2*ex) - 1.5*ds.df0020*ds.df0010/(ex*ex)
        + ds.df0010*ds.df0010*ds.df0010/(ex*ex*ex);
    res->df10 = a*f10;
    res->df01 = a*f01;
    res->df20 = a*(f10*f10 + f20);
    res->df11 = a*(f10*f01 + f11);
    res->df02 = a*(f01*f01 + f02);

    /* terms will partially cancel out - think how to simplify them */
    res->df30 = a*(f10*f10*f10 + 3*f10*f20 + f30);
    res->df21 = a*(f10*f10*f01 + f20*f01 + 2*f10*f11 + f21);
    res->df12 = a*(f10*f01*f01 + 2*f11*f01 + f10*f02 + f12);
    res->df03 = a*(f01*f01*f01 + 3*f01*f02 + f03);
}

/* ===================================================================
 * The expansion of the B-factor and its first derivative for small
 * values of a.
 * =================================================================== */
#define SQRT_PI 1.77245385090552
static __inline__ real
camb3lyp_b_energy_small(real a)
{
    real res;
    a = 2*a; /* the expension derived for different a; correct for this. */
    res = 1-4.0/3.0*SQRT_PI*a + 2 * a*a - 2.0/3.0*a*a*a*a;
    return 1-ALPHA -BETA*(1-res);
}


static __inline__ real
camb3lyp_b_first_small(real a)
{
    real res;
    a = 2*a; /* the expension derived for different a; correct for this. */
    res = 4.0/3.0*(-SQRT_PI + 3 * a +(2*exp(-1/(a*a)) - 2.0)*a*a*a);
    return 2*BETA*res;
}

/* ===================================================================
 * The expansion of the B factor and its first derivative for large
 * values of a.
 * =================================================================== */
#define MAX_LARGE_COEFS 5
static const real large_coefs[] = { 9, 60, 420, 3240, 27720 };
static __inline__ real
camb3lyp_b_energy_large(real a)
{
    real res, ac, a2;
    int i;
    
    a = 2*a; /* the expension derived for different a; correct for this. */
    a2 = a*a;
    res = 0;
    for(i=0, ac = a2; i<MAX_LARGE_COEFS; i++, ac *= -a2)
        res += 1.0/(large_coefs[i]*ac);
    return 1-ALPHA - BETA*(1-res);
}


static const real large_coefs1[] = { 4.5, 15, 70, 405 };
static __inline__ real
camb3lyp_b_first_large(real a)
{
    real tmp;
    real ac, a2;
    int i;

    a = 2*a; /* the expension derived for different a; correct for this. */
    a2  = a*a;
    tmp = 0;
    for(i=0, ac = -a2*a; i<MAX_LARGE_COEFS-1; i++) {
        tmp += 1.0/(large_coefs1[i]*ac);
        ac *= -a2;
    }
    return 2*tmp*BETA;
}

/* ===================================================================
 * The expansion of the B factor and its first, second and third
 * derivatives for medium values of a. This part uses the full
 * expression.
 * 8/3*a*(sqrt(%PI)*erf(1/(2*a))+2*a*(b-c))
 * =================================================================== */
static __inline__ real
camb3lyp_b_energy_medium(real a)
{
    real b = exp(-1/(4*a*a))-1;
    real c = 2*a*a*b + 0.5;
    real res= 1-ALPHA - BETA*8.0/3.0*a*(SQRT_PI*erf(1/(2*a))+2*a*(b-c));
    return res;
}

/* tested version */
static real
camb3lyp_b_first_medium(real a)
{
    real t1 = 1/a;
    real t2 = a*a;
    real t3 = 1/t2;
    real t4 = exp(-0.25*t3);
    real t5 = t4-1;
    real t6 = t4-2*t2*t5-1.5;
    real res = -2.666666666666667*a
        *(2*a*(t4/(2*pow(a,3.0))-4*a*t5-t1*t4)+2*t6-t3*t4)
        -2.666666666666667*(2*a*t6+SQRT_PI*erf(0.5*t1));
    return BETA*res;
}

static real
evaluate_series(int n, const real*coefs, real lambda)
{
    real res = 0, ac =1.0;
    int i;
    for(i=0; i<n; i++) {
        res += 1.0/(coefs[i]*ac);
        ac *= lambda;
    }
    return res;
}
static __inline__ real
camb3lyp_b_second_medium(real a)
{
    real t1, t2, t3, t4, t5, t6, t7;

    if(fabs(a)<1e-40) return 16.0;
    t1 = a*a;
    if(a>=5)  {
        static const double large_coefs[] = {  6, -48, 640,  -11520 };
        return evaluate_series(ELEMENTS(large_coefs), large_coefs, t1)/
            (t1*a);
    }
    t2 = 1/t1;
    t3 = exp(-0.25*t2);
    t4 = pow(a,-3.0);
    t5 = t3-1.0;
    t6 = -t2*t3;
    t7 = -t3/a+0.5*t4*t3-4*a*t5;
    return -(8*a*(2*a*(t3/(4*pow(a,6.0))-2*t3*pow(a,-4.0)+t6-4*t5)
                  -t3/(2*pow(a,5.0))+4*t7+2*t4*t3)/3.0
             +16*(2.0*a*t7+2.0*(t3-2*t1*t5-1.5)+t6)/3.0)*BETA;
}

static __inline__ real
camb3lyp_b_third_medium(real a)
{
    real t1, t2, t3, t4, t5, t6, t7, t8;

    if(fabs(a)<1e-40) return 0;
    if(a>=5)  {
        static const double large_coefs[] = {  -1.5, 8, -80,  1152 };
        real a2 = a*a;
        return evaluate_series(ELEMENTS(large_coefs), large_coefs, a2)/
            (a2*a2*a);
    }
    t1 = pow(a,-2.0);
    t2 = exp(-0.25*t1);
    t3 = pow(a,-6.0);
    t4 = pow(a,-4.0);
    t5 = pow(a,-5.0);
    t6 = t2-1;
    t7 = -t1*t2-2*t4*t2+t3*t2/4-4*t6;
    t8 = pow(a,-3.0);
    return -8*BETA*
        (a*(2*a*(t2/(8*pow(a,9))-5*t2/(2*pow(a,7))+15*t5*t2/2)
                 -t2/(4*pow(a,8))+6*t7-6*t4*t2+7*t3*t2/2)/3
         +(4*(-t2/a+t8*t2/2-4*a*t6)+2*a*t7+2*t8*t2-t5*t2/2));
}

static __inline__ real
camb3lyp_b_fourth_medium(real a)
{
    real t1, t2, t3, t4, t5, t6, t7, t8, t9;

    if(a < 0.05) return -256;
    if(a>=5)  {
        static const double large_coefs[] =
            { 0.3, -8.0/7.0, 80.0/9.0, -1152.0/11.0 };
        real a2 = a*a;
        return evaluate_series(ELEMENTS(large_coefs), large_coefs, a2)/
            (a2*a2*a2);
    }
    t1 = pow(a,-2.0);
    t2 = exp(-0.25*t1);
    t3 = pow(a,-9.0);
    t4 = pow(a,-7.0);
    t5 = pow(a,-5.0);
    t6 = pow(a,-8.0);
    t7 = pow(a,-6.0);
    t8 = 15*t5*t2/2-5*t4*t2/2+t3*t2/8;
    t9 = pow(a,-4.0);

    return -8*a*(2.0*a*(t2/(16*pow(a,12))-19.0*t2/(8.0*pow(a,10.0))
                      -75*t7*t2/2+85*t6*t2*0.25)-t2/(8*pow(a,11.0))
                 +8*t8+24*t5*t2-24*t4*t2+15*0.25*t3*t2)/3.0
        -32*(6*(-t1*t2-2*t9*t2+0.25*t7*t2-4.0*(t2-1))+2*a*t8-6*t9*t2
             +7*t7*t2/2-t6*t2/4)/3;
}

#define FAC M_SQRT2
#define EVALUATOR(a,type) \
  ((a<0.1) ? camb3lyp_b_ ## type ## _small(a*FAC)  : \
  ((a<3)   ? camb3lyp_b_ ## type ## _medium(a*FAC) : \
   camb3lyp_b_ ## type ## _large(a*FAC)))

static real
camb3lyp_energy_sigma(real rho, real grad, real a)
{
    real bfactor = EVALUATOR(a, energy);
    return -2*18*M_PI/(MU*MU)*rho*rho*a*a*bfactor;
}

static real
camb3lyp_energy(const DftDensProp *dp)
{
    real res, ea, eb, a = fun_a(dp->rhoa, dp->grada);
    
    ea = camb3lyp_energy_sigma(dp->rhoa, dp->grada, a);
    if(fabs(dp->rhoa-dp->rhob)>1e-40 || fabs(dp->grada-dp->gradb)>1e-40) {
        a = fun_a(dp->rhob, dp->gradb);
        eb = camb3lyp_energy_sigma(dp->rhob, dp->gradb, a);
    } else eb = ea;
    res = ea + eb;
#if ADD_CORRELATION
    res +=
        LYPFunctional.func(dp)*LYP_WEIGHT +
        VWNFunctional.func(dp)*VWN_WEIGHT;
#endif
    return res;
}

static void
camb3lyp_first_sigma(real rho, real grad, real a, RGFirstDrv *res)
{
    real bfactor       = EVALUATOR(a, energy);
    real bfactor_first = EVALUATOR(a, first);
    RGFirstDrv ader;
    FirstFuncDrv ds;
    DftDensProp dp;
    real ex;

    fun_a_first(rho, grad, &ader);
    ader.df10 *= FAC; ader.df01 *= FAC;
    dp.rhoa = dp.rhob = rho;
    dp.grada = dp.gradb = grad;
    ex = 0.5*(DiracFunctional.func(&dp) +
              BeckeFunctional.func(&dp)*BECKE88_CORR_WEIGHT);
    memset(&ds, 0, sizeof(ds));
    DiracFunctional.first(&ds, 1, &dp);
    BeckeFunctional.first(&ds, BECKE88_CORR_WEIGHT, &dp);

    res->df10 = ds.df1000*bfactor + ex*bfactor_first*ader.df10;
    res->df01 = ds.df0010*bfactor + ex*bfactor_first*ader.df01;
}

static void
camb3lyp_first(FirstFuncDrv *ds, real factor, const DftDensProp *dp)
{
    RGFirstDrv res;
    real a = fun_a(dp->rhoa, dp->grada);

    camb3lyp_first_sigma(dp->rhoa, dp->grada, a, &res);
    ds->df1000 += factor*res.df10;
    ds->df0010 += factor*res.df01;

    if(dp->rhob>10e-13) {
        if(fabs(dp->rhoa-dp->rhob)>1e-40 || fabs(dp->grada-dp->gradb)>1e-40) {
            a = fun_a(dp->rhob, dp->gradb);
            camb3lyp_first_sigma(dp->rhob, dp->gradb, a, &res);
        }
        ds->df0100 += factor*res.df10;
        ds->df0001 += factor*res.df01;
    }
#if ADD_CORRELATION
    LYPFunctional.first(ds, LYP_WEIGHT*factor, dp);
    VWNFunctional.first(ds, VWN_WEIGHT*factor, dp);
#endif
}

static void
camb3lyp_second_sigma(real rho, real grad, real a, RGSecondDrv *res)
{
    real bfactor  = EVALUATOR(a, energy);
    real b_first  = EVALUATOR(a, first);
    real b_second = camb3lyp_b_second_medium(a*FAC);
    RGSecondDrv ader;
    SecondFuncDrv ds;
    DftDensProp dp;
    real ex;

    fun_a_second(rho, grad, &ader);
    ader.df10 *= FAC; ader.df01 *= FAC;
    ader.df20 *= FAC; ader.df11 *= FAC; ader.df02 *= FAC;
    
    dp.rhoa  = dp.rhob  = rho;
    dp.grada = dp.gradb = grad;
    ex = 0.5*(DiracFunctional.func(&dp) +
              BeckeFunctional.func(&dp)*BECKE88_CORR_WEIGHT);
    memset(&ds, 0, sizeof(ds));
    DiracFunctional.second(&ds, 1, &dp);
    BeckeFunctional.second(&ds, BECKE88_CORR_WEIGHT, &dp);

    res->df10 = ds.df1000*bfactor + ex*b_first*ader.df10;
    res->df01 = ds.df0010*bfactor + ex*b_first*ader.df01;

    res->df20 = ds.df2000*bfactor + 2*ds.df1000*b_first*ader.df10 +
        ex*b_second*ader.df10*ader.df10 + ex*b_first*ader.df20;
    res->df11 = ds.df1010*bfactor + ds.df1000*b_first*ader.df01 + 
        ds.df0010*b_first*ader.df10 + 
        ex*b_second*ader.df10*ader.df01 +
        ex*b_first*ader.df11;
    res->df02 = ds.df0020*bfactor + 2*ds.df0010*b_first*ader.df01 +
        ex*b_second*ader.df01*ader.df01 + ex*b_first*ader.df02;
}

static void
camb3lyp_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp)
{
    RGSecondDrv res;
    real a = fun_a(dp->rhoa, dp->grada);

    camb3lyp_second_sigma(dp->rhoa, dp->grada, a, &res);

    ds->df1000 += factor*res.df10;
    ds->df0010 += factor*res.df01;
    ds->df2000 += factor*res.df20;
    ds->df1010 += factor*res.df11;
    ds->df0020 += factor*res.df02;

    if(fabs(dp->rhoa-dp->rhob)>1e-40 || fabs(dp->grada-dp->gradb)>1e-40) {
        a = fun_a(dp->rhob, dp->gradb);
        camb3lyp_second_sigma(dp->rhob, dp->gradb, a, &res);
    }

    ds->df0100 += factor*res.df10;
    ds->df0001 += factor*res.df01;
    ds->df0200 += factor*res.df20;
    ds->df0101 += factor*res.df11;
    ds->df0002 += factor*res.df02;

#if ADD_CORRELATION
    LYPFunctional.second(ds, LYP_WEIGHT*factor, dp);
    VWNFunctional.second(ds, VWN_WEIGHT*factor, dp);
#endif
}

/* ===================================================================
   Third order derivatives specific code.
   =================================================================== */

static void
camb3lyp_third_sigma(real rho, real grad, real a, RGThirdDrv *res)
{
    real bfactor  = EVALUATOR(a, energy);
    real b_first  = EVALUATOR(a, first);
    real b_second = camb3lyp_b_second_medium(a*FAC);
    real b_third  = camb3lyp_b_third_medium(a*FAC);
    RGThirdDrv ader;
    ThirdFuncDrv ds;
    DftDensProp dp;
    real ex;

    fun_a_third(rho, grad, &ader);
    ader.df10 *= FAC; ader.df01 *= FAC;
    ader.df20 *= FAC; ader.df11 *= FAC; ader.df02 *= FAC;
    ader.df30 *= FAC; ader.df21 *= FAC; ader.df12 *= FAC; ader.df03 *= FAC;

    dp.rhoa  = dp.rhob  = rho;
    dp.grada = dp.gradb = grad;
    ex = 0.5*(DiracFunctional.func(&dp) +
              BeckeFunctional.func(&dp)*BECKE88_CORR_WEIGHT);
    memset(&ds, 0, sizeof(ds));
    DiracFunctional.third(&ds, 1, &dp);
    BeckeFunctional.third(&ds, BECKE88_CORR_WEIGHT, &dp);

    res->df10 = ds.df1000*bfactor + ex*b_first*ader.df10;
    res->df01 = ds.df0010*bfactor + ex*b_first*ader.df01;

    res->df20 = ds.df2000*bfactor + 2*ds.df1000*b_first*ader.df10 +
        ex*b_second*ader.df10*ader.df10 + ex*b_first*ader.df20;
    res->df11 = ds.df1010*bfactor + ds.df1000*b_first*ader.df01 + 
        ds.df0010*b_first*ader.df10 + 
        ex*b_second*ader.df10*ader.df01 +
        ex*b_first*ader.df11;
    res->df02 = ds.df0020*bfactor + 2*ds.df0010*b_first*ader.df01 +
        ex*b_second*ader.df01*ader.df01 + ex*b_first*ader.df02;

    res->df30 = ds.df3000*bfactor + 3*ds.df2000*b_first*ader.df10
        + 3*ds.df1000*b_second*ader.df10*ader.df10
        + 3*ds.df1000*b_first *ader.df20 + 3*ex*b_second*ader.df10*ader.df20
        + ex*b_third*ader.df10*ader.df10*ader.df10
        + ex*b_first*ader.df30;

    res->df21 = ds.df2010*bfactor + ex*b_first*ader.df21 +
        ex*b_third*ader.df10*ader.df10*ader.df01 +
        b_first*(ds.df2000*ader.df01 + ds.df0010*ader.df20 + 
                 2*ds.df1010*ader.df10 + 2*ds.df1000*ader.df11) + 
        b_second*(ex*ader.df20*ader.df01 + ds.df0010*ader.df10*ader.df10 +
                  2*ex*ader.df10*ader.df11 + 2*ds.df1000*ader.df10*ader.df01);

    res->df12 = ds.df1020*bfactor + ex*b_first*ader.df12 +
        ex*b_third*ader.df10*ader.df01*ader.df01 +
        b_first*(ds.df0020*ader.df10 + ds.df1000*ader.df02 + 
                 2*ds.df1010*ader.df01 + 2*ds.df0010*ader.df11) +
        b_second*(ex*ader.df10*ader.df02 + ds.df1000*ader.df01*ader.df01 +
                  2*ex*ader.df01*ader.df11 + 2*ds.df0010*ader.df10*ader.df01);

    res->df03 = ds.df0030*bfactor + 3*ds.df0020*b_first*ader.df01
        + 3*ds.df0010*b_second*ader.df01*ader.df01
        + 3*ds.df0010*b_first *ader.df02 + 3*ex*b_second*ader.df01*ader.df02
        + ex*b_third*ader.df01*ader.df01*ader.df01
        + ex*b_first*ader.df03;
}


static void
camb3lyp_third(ThirdFuncDrv *ds, real factor, const DftDensProp* dp)
{
    RGThirdDrv res;
    real a = fun_a(dp->rhoa, dp->grada);

    camb3lyp_third_sigma(dp->rhoa, dp->grada, a, &res);

    ds->df1000 += factor*res.df10;
    ds->df0010 += factor*res.df01;
    ds->df2000 += factor*res.df20;
    ds->df1010 += factor*res.df11;
    ds->df0020 += factor*res.df02;

    ds->df3000 += factor*res.df30;
    ds->df2010 += factor*res.df21;
    ds->df1020 += factor*res.df12;
    ds->df0030 += factor*res.df03;

    if(fabs(dp->rhoa-dp->rhob)>1e-40 || fabs(dp->grada-dp->gradb)>1e-40) {
        a = fun_a(dp->rhob, dp->gradb);
        camb3lyp_third_sigma(dp->rhob, dp->gradb, a, &res);
    }

    ds->df0100 += factor*res.df10;
    ds->df0001 += factor*res.df01;
    ds->df0200 += factor*res.df20;
    ds->df0101 += factor*res.df11;
    ds->df0002 += factor*res.df02;

    ds->df0300 += factor*res.df30;
    ds->df0201 += factor*res.df21;
    ds->df0102 += factor*res.df12;
    ds->df0003 += factor*res.df03;

#if ADD_CORRELATION
    LYPFunctional.third(ds, LYP_WEIGHT*factor, dp);
    VWNFunctional.third(ds, VWN_WEIGHT*factor, dp);
#endif
}

/* ===================================================================
   TEST CODE
   =================================================================== */
#ifdef TEST
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char *argv[])
{
    DftDensProp dp;
    RGFirstDrv fa;
    real a, bfactor, bfactor_first, bfactor_second;
    if(argc<4) { printf("rhoa grada mu\n"); return 1;}
    dp.rhoa  = dp.rhob  = atof(argv[1]);
    dp.grada = dp.gradb = atof(argv[2]);
    CamMuFactor = atof(argv[3]);
    CamAlpha = 0.19;
    CamBeta  = 0.46;
    printf("mu=%f: energy(%g,%g) = %g\n", CamMuFactor, dp.rhoa*2, dp.grada*2,
           camb3lyp_energy(&dp));
    fun_a_first(1,1,&fa);
    printf("fun a derivatives: %g %g\n", fa.df10, fa.df01);
    for(a=0.1; a<10; a += 0.1) {
        bfactor_first  = camb3lyp_b_second_medium(a*FAC);//EVALUATOR(a, second);
        real b1 = EVALUATOR((a+1e-6), first);
        real b2 = EVALUATOR((a-1e-6), first);
        printf("%20.15g %20.15g %20.15g %20.15g\n", a,
               bfactor, bfactor_first*FAC, (b1-b2)/2e-6);
    }
    return 1;
}

void dftsethf_(real* s) {}
void dftsetcam_(real* s, real *b) {}
void fort_print(const char *fmt, ...){}
#endif


