/* Automatically generated functional code: pbex
   Maxima input:
    >> pi:3.14159265358979312;
    >> 
    >> // parameters for pbex 
    >> R:0.804;
    >> d:0.066725;
    >> mu:d*pi^2/3;
    >> Sa:xa/(2*(6*pi^2)^(1/3)); 
    >> Sb:xb/(2*(6*pi^2)^(1/3)); 
    >> 
    >> // functions for pbex 
    >> F(S):=1+R-R/(1+mu*S^2/R);
    >> Ea(n):=-3/(4*pi)*(3*pi^2)^(1/3)*n^(4/3)*F(Sa);
    >> Eb(n):=-3/(4*pi)*(3*pi^2)^(1/3)*n^(4/3)*F(Sb);
    >> 
    >> // kernel 
    >> K(rhoa,grada,rhob,gradb,gradab):=0.5*(Ea(2*rhoa)+Eb(2*rhob));
    >> 
    >> 
*/

// add "extern Functional pbexFunctional;" to 'functionals.h'
// add "&pbexFunctional," to 'functionals.c'
// add "fun-pbex.c" to 'Makefile.in'

#include <math.h>
#include <stddef.h>

#define __CVERSION__

#include "functionals.h"
#define LOG log
#define ABS fabs
#define ASINH asinh
#define SQRT sqrt

/* INTERFACE PART */
static int pbex_read(const char* conf_line);
static real pbex_energy(const DftDensProp* dp);
static void pbex_first(FirstFuncDrv *ds, real factor, 
                       const DftDensProp* dp);
static void pbex_second(SecondFuncDrv *ds, real factor,
                        const DftDensProp* dp);
static void pbex_third(ThirdFuncDrv *ds, real factor,
                       const DftDensProp* dp);

Functional PbexFunctional = {
  "pbex",
  fun_true,
  pbex_read,
  NULL,
  pbex_energy,
  pbex_first,
  pbex_second,
  pbex_third
};

/* IMPLEMENTATION PART */
static int
pbex_read(const char* conf_line)
{
    dft_set_hf_weight(0);
    return 1;
}


static real
pbex_energy(const DftDensProp* dp)
{
    real t[2],zk;
    real rhoa = dp->rhoa;
    real rhob = dp->rhob;
    real grada = dp->grada;
    real gradb = dp->gradb;

    t[1] = pow(2.0,0.33333333333333);
    zk = 0.5*(-1.477117532764045*t[1]*pow(rhob,1.333333333333333)*(1.804-0.804/(0.00449279948023*pow(gradb,2.0)/pow(rhob,2.666666666666667)+1.0))-1.477117532764045*t[1]*pow(rhoa,1.333333333333333)*(1.804-0.804/(0.00449279948023*pow(grada,2.0)/pow(rhoa,2.666666666666667)+1.0)));
    return zk;
}

static void
pbex_first(FirstFuncDrv *ds, real factor, const DftDensProp* dp)
{
    real t[8];
    real dfdra, dfdrb, dfdga, dfdgb, dfdab;
    real rhoa = dp->rhoa;
    real rhob = dp->rhob;
    real grada = dp->grada;
    real gradb = dp->gradb;

    t[1] = pow(2.0,0.33333333333333);
    t[2] = pow(grada,2.0);
    t[3] = 0.00449279948023*t[2]/pow(rhoa,2.666666666666667)+1.0;
    t[4] = 1/pow(t[3],2.0);
    t[5] = pow(gradb,2.0);
    t[6] = 0.00449279948023*t[5]/pow(rhob,2.666666666666667)+1.0;
    t[7] = 1/pow(t[6],2.0);
    dfdra = 0.5*(0.0142284263421*t[1]*t[2]*t[4]/pow(rhoa,2.333333333333334)-1.969490043685393*t[1]*(1.804-0.804/t[3])*pow(rhoa,0.33333333333333));
    dfdrb = 0.5*(0.0142284263421*t[1]*t[5]*t[7]/pow(rhob,2.333333333333334)-1.969490043685393*t[1]*(1.804-0.804/t[6])*pow(rhob,0.33333333333333));
    dfdga = -0.00533565987829*t[1]*t[4]*grada/pow(rhoa,1.333333333333333);
    dfdgb = -0.00533565987829*t[1]*t[7]*gradb/pow(rhob,1.333333333333333);
    dfdab = 0.0;
    ds->df1000 += factor*dfdra;
    ds->df0100 += factor*dfdrb;
    ds->df0010 += factor*dfdga;
    ds->df0001 += factor*dfdgb;
    ds->df00001 += factor*dfdab;
}

static void
pbex_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp)
{
    real t[16];
    real dfdra, dfdrb, dfdga, dfdgb, dfdab;
    real d2fdraga, d2fdrara, d2fdrarb, d2fdragb, d2fdrbrb;
    real d2fdrbgb, d2fdgaga, d2fdgbgb, d2fdrbga;
    real d2fdraab, d2fdrbab;
    real rhoa = dp->rhoa;
    real rhob = dp->rhob;
    real grada = dp->grada;
    real gradb = dp->gradb;

    t[1] = pow(2.0,0.33333333333333);
    t[2] = pow(grada,2.0);
    t[3] = 0.00449279948023*t[2]/pow(rhoa,2.666666666666667)+1.0;
    t[4] = 1/pow(t[3],2.0);
    t[5] = 1/pow(rhoa,2.333333333333334);
    t[6] = 1.804-0.804/t[3];
    t[7] = pow(gradb,2.0);
    t[8] = 0.00449279948023*t[7]/pow(rhob,2.666666666666667)+1.0;
    t[9] = 1/pow(t[8],2.0);
    t[10] = 1/pow(rhob,2.333333333333334);
    t[11] = 1.804-0.804/t[8];
    t[12] = 1/pow(rhoa,1.333333333333333);
    t[13] = 1/pow(rhob,1.333333333333333);
    t[14] = 1/pow(t[3],3.0);
    t[15] = 1/pow(t[8],3.0);
    dfdra = 0.5*(0.0142284263421*t[1]*t[2]*t[4]*t[5]-1.969490043685393*t[1]*t[6]*pow(rhoa,0.33333333333333));
    dfdrb = 0.5*(0.0142284263421*t[1]*t[7]*t[9]*t[10]-1.969490043685393*t[1]*t[11]*pow(rhob,0.33333333333333));
    dfdga = -0.00533565987829*t[1]*grada*t[4]*t[12];
    dfdgb = -0.00533565987829*t[1]*gradb*t[9]*t[13];
    dfdab = 0.0;
    d2fdrara = 0.5*(3.409358211962494E-4*t[1]*t[14]*pow(grada,4.0)/pow(rhoa,6.0)-0.0142284263421*t[1]*t[2]*t[4]/pow(rhoa,3.333333333333334)-0.65649668122846*t[1]*t[6]/pow(rhoa,0.66666666666667));
    d2fdrarb = 0.0;
    d2fdraga = 0.5*(0.0142284263421*t[1]*grada*t[4]*t[5]-2.557018658971871E-4*t[1]*t[14]*pow(grada,3.0)/pow(rhoa,5.0));
    d2fdragb = 0.0;
    d2fdrbrb = 0.5*(3.409358211962494E-4*t[1]*t[15]*pow(gradb,4.0)/pow(rhob,6.0)-0.0142284263421*t[1]*t[7]*t[9]/pow(rhob,3.333333333333334)-0.65649668122846*t[1]*t[11]/pow(rhob,0.66666666666667));
    d2fdraab = 0.0;
    d2fdrbab = 0.0;
    d2fdgaga = 9.588819971144516E-5*t[1]*t[14]*t[2]/pow(rhoa,4.0)-0.00533565987829*t[1]*t[4]*t[12];
    d2fdgbgb = 9.588819971144516E-5*t[1]*t[15]*t[7]/pow(rhob,4.0)-0.00533565987829*t[1]*t[9]*t[13];
    d2fdrbga = 0.0;
    d2fdrbgb = 0.5*(0.0142284263421*t[1]*gradb*t[9]*t[10]-2.557018658971871E-4*t[1]*t[15]*pow(gradb,3.0)/pow(rhob,5.0));
    ds->df1000 += factor*dfdra;
    ds->df0100 += factor*dfdrb;
    ds->df0010 += factor*dfdga;
    ds->df0001 += factor*dfdgb;
    ds->df00001 += factor*dfdab;
    ds->df2000 += factor*d2fdrara;
    ds->df1100 += factor*d2fdrarb;
    ds->df1010 += factor*d2fdraga;
    ds->df1001 += factor*d2fdragb;
    ds->df10001 += factor*d2fdraab;
    ds->df0200 += factor*d2fdrbrb;
    ds->df0110 += factor*d2fdrbga;
    ds->df0101 += factor*d2fdrbgb;
    ds->df01001 += factor*d2fdrbab;
    ds->df0020 += factor*d2fdgaga;
    ds->df0002 += factor*d2fdgbgb;
}

static void
pbex_third(ThirdFuncDrv *ds, real factor, const DftDensProp* dp)
{
    real t[29];
    real dfdra, dfdrb, dfdga, dfdgb, dfdab;
    real d2fdraga, d2fdrara, d2fdrarb, d2fdragb, d2fdrbrb;
    real d2fdrbgb, d2fdgaga, d2fdgbgb, d2fdrbga;
    real d2fdraab, d2fdrbab;
    real d3fdraraga, d3fdraragb, d3fdraraab, d3fdrbrbab;
    real d3fdrarara, d3fdrararb, d3fdragaga, d3fdrarbrb;
    real d3fdragbgb, d3fdrarbgb, d3fdrarbab, d3fdgagaga;
    real d3fdrbrbrb, d3fdrbrbga, d3fdrbrbgb, d3fdrbgbgb;
    real d3fdrbgbga, d3fdrarbga, d3fdrbgaga;
    real rhoa = dp->rhoa;
    real rhob = dp->rhob;
    real grada = dp->grada;
    real gradb = dp->gradb;

    t[1] = pow(2.0,0.33333333333333);
    t[2] = pow(grada,2.0);
    t[3] = 0.00449279948023*t[2]/pow(rhoa,2.666666666666667)+1.0;
    t[4] = 1/pow(t[3],2.0);
    t[5] = 1/pow(rhoa,2.333333333333334);
    t[6] = 1.804-0.804/t[3];
    t[7] = pow(gradb,2.0);
    t[8] = 0.00449279948023*t[7]/pow(rhob,2.666666666666667)+1.0;
    t[9] = 1/pow(t[8],2.0);
    t[10] = 1/pow(rhob,2.333333333333334);
    t[11] = 1.804-0.804/t[8];
    t[12] = 1/pow(rhoa,1.333333333333333);
    t[13] = 1/pow(rhob,1.333333333333333);
    t[14] = pow(grada,4.0);
    t[15] = 1/pow(t[3],3.0);
    t[16] = 1/pow(rhoa,6.0);
    t[17] = 1/pow(rhoa,3.333333333333334);
    t[18] = pow(grada,3.0);
    t[19] = 1/pow(rhoa,5.0);
    t[20] = pow(gradb,4.0);
    t[21] = 1/pow(t[8],3.0);
    t[22] = 1/pow(rhob,6.0);
    t[23] = 1/pow(rhob,3.333333333333334);
    t[24] = 1/pow(rhoa,4.0);
    t[25] = pow(gradb,3.0);
    t[26] = 1/pow(rhob,5.0);
    t[27] = 1/pow(t[3],4.0);
    t[28] = 1/pow(t[8],4.0);
    dfdra = 0.5*(0.0142284263421*t[1]*t[2]*t[4]*t[5]-1.969490043685393*t[1]*t[6]*pow(rhoa,0.33333333333333));
    dfdrb = 0.5*(0.0142284263421*t[1]*t[7]*t[9]*t[10]-1.969490043685393*t[1]*t[11]*pow(rhob,0.33333333333333));
    dfdga = -0.00533565987829*t[1]*grada*t[4]*t[12];
    dfdgb = -0.00533565987829*t[1]*gradb*t[9]*t[13];
    dfdab= 0.0;
    d2fdrara = 0.5*(-0.65649668122846*t[1]*t[6]/pow(rhoa,0.66666666666667)-0.0142284263421*t[1]*t[2]*t[4]*t[17]+3.409358211962494E-4*t[1]*t[14]*t[15]*t[16]);
    d2fdrarb = 0.0;
    d2fdraga = 0.5*(0.0142284263421*t[1]*grada*t[4]*t[5]-2.557018658971871E-4*t[1]*t[18]*t[15]*t[19]);
    d2fdragb = 0.0;
    d2fdrbrb = 0.5*(-0.65649668122846*t[1]*t[11]/pow(rhob,0.66666666666667)-0.0142284263421*t[1]*t[7]*t[9]*t[23]+3.409358211962494E-4*t[1]*t[20]*t[21]*t[22]);
    d2fdraab = 0.0;
    d2fdrbab = 0.0;
    d2fdgaga = 9.588819971144516E-5*t[1]*t[2]*t[15]*t[24]-0.00533565987829*t[1]*t[4]*t[12];
    d2fdgbgb = 9.588819971144516E-5*t[1]*t[21]*t[7]/pow(rhob,4.0)-0.00533565987829*t[1]*t[9]*t[13];
    d2fdrbga = 0.0;
    d2fdrbgb = 0.5*(0.0142284263421*t[1]*gradb*t[9]*t[10]-2.557018658971871E-4*t[1]*t[25]*t[21]*t[26]);
    d3fdrararb = 0.0;
    d3fdraraga = 0.5*(-9.190537681576017E-6*t[1]*t[27]*pow(grada,5.0)/pow(rhoa,8.666666666666666)-0.0331996614649*t[1]*grada*t[4]*t[17]+0.00161944515068*t[1]*t[18]*t[15]*t[16]);
    d3fdraragb =0.0;
    d3fdrbrbab = 0.0;
    d3fdraraab = 0.0;
    d3fdrarbrb = 0.0;
    d3fdrarbga = 0.0;
    d3fdrarbgb = 0.0;
    d3fdrarbab = 0.0;
    d3fdragaga = 0.5*(6.892903261182013E-6*t[1]*t[14]*t[27]/pow(rhoa,7.666666666666667)+0.0142284263421*t[1]*t[4]*t[5]-0.00102280746359*t[1]*t[2]*t[15]*t[19]);
    d3fdragbgb = 0.0;
    d3fdrarara = 0.5*(1.225405024210135E-5*t[1]*t[27]*pow(grada,6.0)/pow(rhoa,9.666666666666666)-0.00238655074837*t[1]*t[14]*t[15]/pow(rhoa,7.0)+0.05375183284794*t[1]*t[2]*t[4]/pow(rhoa,4.333333333333333)+0.43766445415231*t[1]*t[6]/pow(rhoa,1.666666666666667));
    d3fdrbrbrb = 0.5*(1.225405024210135E-5*t[1]*t[28]*pow(gradb,6.0)/pow(rhob,9.666666666666666)-0.00238655074837*t[1]*t[20]*t[21]/pow(rhob,7.0)+0.05375183284794*t[1]*t[7]*t[9]/pow(rhob,4.333333333333333)+0.43766445415231*t[1]*t[11]/pow(rhob,1.666666666666667));
    d3fdrbrbga = 0.0;
    d3fdrbrbgb = 0.5*(-9.190537681576017E-6*t[1]*t[28]*pow(gradb,5.0)/pow(rhob,8.666666666666666)-0.0331996614649*t[1]*gradb*t[9]*t[23]+0.00161944515068*t[1]*t[25]*t[21]*t[22]);
    d3fdrbgaga = 0.0;
    d3fdrbgbga = 0.0;
    d3fdrbgbgb = 0.5*(6.892903261182013E-6*t[1]*t[20]*t[28]/pow(rhob,7.666666666666667)-0.00102280746359*t[1]*t[7]*t[21]*t[26]+0.0142284263421*t[1]*t[9]*t[10]);
    d3fdgagaga = 2.876645991343355E-4*t[1]*grada*t[15]*t[24]-2.584838722943255E-6*t[1]*t[18]*t[27]/pow(rhoa,6.666666666666667);
    ds->df1000 += factor*dfdra;
    ds->df0100 += factor*dfdrb;
    ds->df0010 += factor*dfdga;
    ds->df0001 += factor*dfdgb;
    ds->df00001 += factor*dfdab;
    ds->df2000 += factor*d2fdrara;
    ds->df1100 += factor*d2fdrarb;
    ds->df1010 += factor*d2fdraga;
    ds->df1001 += factor*d2fdragb;
    ds->df10001 += factor*d2fdraab;
    ds->df0200 += factor*d2fdrbrb;
    ds->df0110 += factor*d2fdrbga;
    ds->df0101 += factor*d2fdrbgb;
    ds->df01001 += factor*d2fdrbab;
    ds->df0020 += factor*d2fdgaga;
    ds->df0002 += factor*d2fdgbgb;
    ds->df2010 += factor*d3fdraraga;
    ds->df2001 += factor*d3fdraragb;
    ds->df1101 += factor*d3fdrarbgb;
    ds->df11001 += factor*d3fdrarbab;
    ds->df1020 += factor*d3fdragaga;
    ds->df1002 += factor*d3fdragbgb;
    ds->df3000 += factor*d3fdrarara;
    ds->df2100 += factor*d3fdrararb;
    ds->df20001 += factor*d3fdraraab;
    ds->df02001 += factor*d3fdrbrbab;
    ds->df1200 += factor*d3fdrarbrb;
    ds->df1110 += factor*d3fdrarbga;
    ds->df0300 += factor*d3fdrbrbrb;
    ds->df0210 += factor*d3fdrbrbga;
    ds->df0201 += factor*d3fdrbrbgb;
    ds->df0120 += factor*d3fdrbgaga;
    ds->df0102 += factor*d3fdrbgbgb;
    ds->df0030 += factor*d3fdgagaga;
}
