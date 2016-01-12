/* Automatically generated functional code: pbex
   Maxima input:
    >> pi:3.14159265358979312;
    >> 
    >> xa:sqrt(grada*grada)/rhoa^(4/3);
    >> xb:sqrt(gradb*gradb)/rhob^(4/3);
    >> 
    >> R:0.804;
    >> d:0.066725;
    >> mu:d*pi^2/3;
    >> Sa:xa/(2*(6*pi^2)^(1/3)); 
    >> Sb:xb/(2*(6*pi^2)^(1/3)); 
    >> 
    >> F(S):=1+R-R/(1+mu*S^2/R);
    >> Ea(n):=-3/(4*pi)*(3*pi^2)^(1/3)*n^(4/3)*F(Sa);
    >> Eb(n):=-3/(4*pi)*(3*pi^2)^(1/3)*n^(4/3)*F(Sb);
    >> 
    >> K(rhoa,grada,rhob,gradb,gradab):=0.5*(Ea(2*rhoa)+Eb(2*rhob));
*/

// add "extern Functional pbexFunctional;" to 'functionals.h'
// add "&pbexFunctional," to 'functionals.c'
// add "fun-pbex.c" to 'Makefile.in'

#include <math.h>
#include <stddef.h>
#include "general.h"

#define __CVERSION__

#include "functionals.h"
#define LOG log
#define ABS fabs
#define ASINH asinh
#define SQRT sqrt

/* INTERFACE PART */
static integer pbex_isgga(void) {return 1;}
static integer pbex_read(const char* conf_line);
static real pbex_energy(const FunDensProp* dp);
static void pbex_first(FunFirstFuncDrv *ds, real factor, 
                           const FunDensProp* dp);
static void pbex_second(FunSecondFuncDrv *ds, real factor,
                            const FunDensProp* dp);
static void pbex_third(FunThirdFuncDrv *ds, real factor,
                           const FunDensProp* dp);

static void pbex_fourth(FunFourthFuncDrv *ds, real factor,
                           const FunDensProp* dp);

//static integer fun_true(void) { return 1; }
Functional PBExFunctional = {
  "PBEx",
  pbex_isgga,
  3,
  pbex_read,
  NULL,
  pbex_energy,
  pbex_first,
  pbex_second,
  pbex_third,
  pbex_fourth
};

/* IMPLEMENTATION PART */
static integer
pbex_read(const char* conf_line)
{
    fun_set_hf_weight(0);
    return 1;
}


static real
pbex_energy(const FunDensProp* dp)
{
    real t[2],zk;
    real rhoa = dp->rhoa;
    real rhob = dp->rhob;
    real grada = dp->grada;
    real gradb = dp->gradb;
    real gradab = dp->gradab;

    t[1] = pow(2.0,1.333333333333333);
    zk = 0.5*(-0.73855876638202*t[1]*pow(rhob,1.333333333333333)*(1.804-0.804/(0.0044927994802311*pow(gradb,2.0)/pow(rhob,2.666666666666667)+1.0))-0.73855876638202*t[1]*pow(rhoa,1.333333333333333)*(1.804-0.804/(0.0044927994802311*pow(grada,2.0)/pow(rhoa,2.666666666666667)+1.0)));
    return zk;
}

static void
pbex_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real t[9];
    real dfdra, dfdrb, dfdga, dfdgb, dfdab;
    real rhoa = dp->rhoa;
    real rhob = dp->rhob;
    real grada = dp->grada;
    real gradb = dp->gradb;
    real gradab = dp->gradab;

    t[1] = pow(2.0,1.333333333333333);
    t[2] = pow(grada,2.0);
    t[3] = 0.0044927994802311*t[2]/pow(rhoa,2.666666666666667)+1.0;
    t[4] = 1/pow(t[3],2.0);
    t[5] = pow(2.0,3.333333333333334);
    t[6] = pow(gradb,2.0);
    t[7] = 0.0044927994802311*t[6]/pow(rhob,2.666666666666667)+1.0;
    t[8] = 1/pow(t[7],2.0);
    dfdra = 0.5*(0.0071142131710504*t[1]*t[2]*t[4]/pow(rhoa,2.333333333333334)-0.24618625546067*(1.804-0.804/t[3])*t[5]*pow(rhoa,0.33333333333333));
    dfdrb = 0.5*(0.0071142131710504*t[1]*t[6]*t[8]/pow(rhob,2.333333333333334)-0.24618625546067*t[5]*(1.804-0.804/t[7])*pow(rhob,0.33333333333333));
    dfdga = -0.0026678299391439*t[1]*t[4]*grada/pow(rhoa,1.333333333333333);
    dfdgb = -0.0026678299391439*t[1]*t[8]*gradb/pow(rhob,1.333333333333333);
    dfdab = 0.0;
    ds->df1000 += factor*dfdra;
    ds->df0100 += factor*dfdrb;
    ds->df0010 += factor*dfdga;
    ds->df0001 += factor*dfdgb;
    ds->df00001 += factor*dfdab;
}

static void
pbex_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real t[20];
    real dfdra, dfdrb, dfdga, dfdgb, dfdab;
    real d2fdraga, d2fdrara, d2fdrarb, d2fdragb, d2fdrbrb;
    real d2fdrbgb, d2fdgaga, d2fdgbgb, d2fdrbga;
    real d2fdraab, d2fdrbab;
    real d2fdgaab, d2fdgbab, d2fdabab, d2fdgagb;
    real rhoa = dp->rhoa;
    real rhob = dp->rhob;
    real grada = dp->grada;
    real gradb = dp->gradb;
    real gradab = dp->gradab;

    t[1] = pow(2.0,1.333333333333333);
    t[2] = pow(grada,2.0);
    t[3] = 0.0044927994802311*t[2]/pow(rhoa,2.666666666666667)+1.0;
    t[4] = 1/pow(t[3],2.0);
    t[5] = 1/pow(rhoa,2.333333333333334);
    t[6] = pow(2.0,3.333333333333334);
    t[7] = 1.804-0.804/t[3];
    t[8] = pow(gradb,2.0);
    t[9] = 0.0044927994802311*t[8]/pow(rhob,2.666666666666667)+1.0;
    t[10] = 1/pow(t[9],2.0);
    t[11] = 1/pow(rhob,2.333333333333334);
    t[12] = 1.804-0.804/t[9];
    t[13] = 1/pow(rhoa,1.333333333333333);
    t[14] = 1/pow(rhob,1.333333333333333);
    t[15] = 1/pow(t[3],3.0);
    t[16] = 1/pow(rhoa,3.333333333333334);
    t[17] = pow(2.0,2.333333333333334);
    t[18] = 1/pow(t[9],3.0);
    t[19] = 1/pow(rhob,3.333333333333334);
    dfdra = 0.5*(0.0071142131710504*t[1]*t[2]*t[4]*t[5]-0.24618625546067*t[6]*t[7]*pow(rhoa,0.33333333333333));
    dfdrb = 0.5*(0.0071142131710504*t[1]*t[8]*t[10]*t[11]-0.24618625546067*t[12]*t[6]*pow(rhob,0.33333333333333));
    dfdga = -0.0026678299391439*t[1]*grada*t[4]*t[13];
    dfdgb = -0.0026678299391439*t[1]*gradb*t[10]*t[14];
    dfdab = 0.0;
    d2fdrara = 0.5*(1.704679105981247E-4*t[1]*t[15]*pow(grada,4.0)/pow(rhoa,6.0)-0.082062085153558*t[6]*t[7]/pow(rhoa,0.66666666666667)+0.0023714043903501*t[6]*t[2]*t[4]*t[16]-0.016599830732451*t[1]*t[2]*t[4]*t[16]);
    d2fdrarb = 0.0;
    d2fdraga = 0.5*(-1.2785093294859353E-4*t[1]*t[15]*pow(grada,3.0)/pow(rhoa,5.0)-0.0017785532927626*t[6]*grada*t[4]*t[5]+0.0071142131710504*t[17]*grada*t[4]*t[5]);
    d2fdragb = 0.0;
    d2fdraab = 0.0;
    d2fdrbrb = 0.5*(1.704679105981247E-4*t[1]*t[18]*pow(gradb,4.0)/pow(rhob,6.0)-0.082062085153558*t[12]*t[6]/pow(rhob,0.66666666666667)+0.0023714043903501*t[6]*t[8]*t[10]*t[19]-0.016599830732451*t[1]*t[8]*t[10]*t[19]);
    d2fdrbga = 0.0;
    d2fdrbgb = 0.5*(-1.2785093294859353E-4*t[1]*t[18]*pow(gradb,3.0)/pow(rhob,5.0)-0.0017785532927626*t[6]*gradb*t[10]*t[11]+0.0071142131710504*t[17]*gradb*t[10]*t[11]);
    d2fdrbab = 0.0;
    d2fdgaga = 4.7944099855722582E-5*t[1]*t[15]*t[2]/pow(rhoa,4.0)-0.0026678299391439*t[1]*t[4]*t[13];
    d2fdgagb = 0.0;
    d2fdgaab = 0.0;
    d2fdgbgb = 4.7944099855722582E-5*t[1]*t[18]*t[8]/pow(rhob,4.0)-0.0026678299391439*t[1]*t[10]*t[14];
    d2fdgbab = 0.0;
    d2fdabab = 0.0;
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
    ds->df0011 += factor*d2fdgagb;
    ds->df00101+= factor*d2fdgaab;
    ds->df0002 += factor*d2fdgbgb;
    ds->df00011+= factor*d2fdgbab;
    ds->df00002+= factor*d2fdabab;
}

static void
pbex_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real t[37];
    real dfdra, dfdrb, dfdga, dfdgb, dfdab;
    real d2fdrara , d2fdrarb , d2fdraga , d2fdragb ;
    real d2fdraab , d2fdrbrb , d2fdrbga , d2fdrbgb ;
    real d2fdrbab , d2fdgaga , d2fdgagb , d2fdgaab ;
    real d2fdgbgb , d2fdgbab , d2fdabab ;
    real d3fdrarara , d3fdrararb , d3fdraraga , d3fdraragb ;
    real d3fdraraab , d3fdrarbrb , d3fdrarbga , d3fdrarbgb ;
    real d3fdrarbab , d3fdragaga , d3fdragagb , d3fdragaab ;
    real d3fdragbgb , d3fdragbab , d3fdraabab , d3fdrbrbrb ;
    real d3fdrbrbga , d3fdrbrbgb , d3fdrbrbab , d3fdrbgaga ;
    real d3fdrbgagb , d3fdrbgaab , d3fdrbgbgb , d3fdrbgbab ;
    real d3fdrbabab , d3fdgagaga , d3fdgagagb , d3fdgagaab ;
    real d3fdgagbgb , d3fdgagbab , d3fdgaabab , d3fdgbgbgb ;
    real d3fdgbgbab , d3fdgbabab , d3fdababab ;
    real rhoa = dp->rhoa;
    real rhob = dp->rhob;
    real grada = dp->grada;
    real gradb = dp->gradb;
    real gradab = dp->gradab;

    t[1] = pow(2.0,1.333333333333333);
    t[2] = pow(grada,2.0);
    t[3] = 0.0044927994802311*t[2]/pow(rhoa,2.666666666666667)+1.0;
    t[4] = 1/pow(t[3],2.0);
    t[5] = 1/pow(rhoa,2.333333333333334);
    t[6] = pow(2.0,3.333333333333334);
    t[7] = 1.804-0.804/t[3];
    t[8] = pow(gradb,2.0);
    t[9] = 0.0044927994802311*t[8]/pow(rhob,2.666666666666667)+1.0;
    t[10] = 1/pow(t[9],2.0);
    t[11] = 1/pow(rhob,2.333333333333334);
    t[12] = 1.804-0.804/t[9];
    t[13] = 1/pow(rhoa,1.333333333333333);
    t[14] = 1/pow(rhob,1.333333333333333);
    t[15] = pow(grada,4.0);
    t[16] = 1/pow(t[3],3.0);
    t[17] = 1/pow(rhoa,6.0);
    t[18] = 1/pow(rhoa,3.333333333333334);
    t[19] = pow(grada,3.0);
    t[20] = 1/pow(rhoa,5.0);
    t[21] = pow(2.0,2.333333333333334);
    t[22] = pow(gradb,4.0);
    t[23] = 1/pow(t[9],3.0);
    t[24] = 1/pow(rhob,6.0);
    t[25] = 1/pow(rhob,3.333333333333334);
    t[26] = pow(gradb,3.0);
    t[27] = 1/pow(rhob,5.0);
    t[28] = 1/pow(rhoa,4.0);
    t[29] = 1/pow(rhob,4.0);
    t[30] = 1/pow(t[3],4.0);
    t[31] = 1/pow(rhoa,7.0);
    t[32] = 1/pow(rhoa,4.333333333333333);
    t[33] = pow(2.0,4.333333333333333);
    t[34] = 1/pow(t[9],4.0);
    t[35] = 1/pow(rhob,7.0);
    t[36] = 1/pow(rhob,4.333333333333333);
    dfdra = 0.5*(0.0071142131710504*t[1]*t[2]*t[4]*t[5]-0.24618625546067*t[6]*t[7]*pow(rhoa,0.33333333333333));
    dfdrb = 0.5*(0.0071142131710504*t[1]*t[8]*t[10]*t[11]-0.24618625546067*t[12]*t[6]*pow(rhob,0.33333333333333));
    dfdga = -0.0026678299391439*t[1]*grada*t[4]*t[13];
    dfdgb = -0.0026678299391439*t[1]*gradb*t[10]*t[14];
    dfdab = 0.0;
    d2fdrara = 0.5*(-0.082062085153558*t[6]*t[7]/pow(rhoa,0.66666666666667)+0.0023714043903501*t[6]*t[2]*t[4]*t[18]-0.016599830732451*t[1]*t[2]*t[4]*t[18]+1.704679105981247E-4*t[1]*t[15]*t[16]*t[17]);
    d2fdrarb = 0.0;
    d2fdraga = 0.5*(-0.0017785532927626*t[6]*grada*t[4]*t[5]+0.0071142131710504*t[21]*grada*t[4]*t[5]-1.2785093294859353E-4*t[1]*t[19]*t[16]*t[20]);
    d2fdragb = 0.0;
    d2fdraab = 0.0;
    d2fdrbrb = 0.5*(-0.082062085153558*t[12]*t[6]/pow(rhob,0.66666666666667)+0.0023714043903501*t[6]*t[8]*t[10]*t[25]-0.016599830732451*t[1]*t[8]*t[10]*t[25]+1.704679105981247E-4*t[1]*t[22]*t[23]*t[24]);
    d2fdrbga = 0.0;
    d2fdrbgb = 0.5*(-0.0017785532927626*t[6]*gradb*t[10]*t[11]+0.0071142131710504*t[21]*gradb*t[10]*t[11]-1.2785093294859353E-4*t[1]*t[26]*t[23]*t[27]);
    d2fdrbab = 0.0;
    d2fdgaga = 4.7944099855722582E-5*t[1]*t[2]*t[16]*t[28]-0.0026678299391439*t[1]*t[4]*t[13];
    d2fdgagb = 0.0;
    d2fdgaab = 0.0;
    d2fdgbgb = 4.7944099855722582E-5*t[1]*t[8]*t[23]*t[29]-0.0026678299391439*t[1]*t[10]*t[14];
    d2fdgbab = 0.0;
    d2fdabab = 0.0;
    d3fdrarara = 0.5*(6.1270251210506771E-6*t[1]*t[30]*pow(grada,6.0)/pow(rhoa,9.666666666666666)+0.027354028384519*t[33]*t[7]/pow(rhoa,1.666666666666667)+7.9046813011671027E-4*t[6]*t[2]*t[4]*t[32]-0.0039523406505836*t[33]*t[2]*t[4]*t[32]+0.027666384554085*t[21]*t[2]*t[4]*t[32]+5.6822636866041565E-5*t[6]*t[15]*t[16]*t[31]-5.1140373179437413E-4*t[21]*t[15]*t[16]*t[31]-3.9775845806229101E-4*t[1]*t[15]*t[16]*t[31]);
    d3fdrararb = 0.0;
    d3fdraraga = 0.5*(-4.5952688407880086E-6*t[1]*t[30]*pow(grada,5.0)/pow(rhoa,8.666666666666666)-5.9285109758753281E-4*t[6]*grada*t[4]*t[18]+0.0023714043903501*t[33]*grada*t[4]*t[18]-0.016599830732451*t[21]*grada*t[4]*t[18]+1.2785093294859353E-4*t[6]*t[19]*t[16]*t[17]+2.9831884354671826E-4*t[1]*t[19]*t[16]*t[17]);
    d3fdraragb = 0.0;
    d3fdraraab = 0.0;
    d3fdrarbrb = 0.0;
    d3fdrarbga = 0.0;
    d3fdrarbgb = 0.0;
    d3fdrarbab = 0.0;
    d3fdragaga = 0.5*(3.4464516305910063E-6*t[1]*t[15]*t[30]/pow(rhoa,7.666666666666667)-0.0017785532927626*t[6]*t[4]*t[5]+0.0071142131710504*t[21]*t[4]*t[5]+3.1962733237148383E-5*t[6]*t[2]*t[16]*t[20]-1.2785093294859353E-4*t[21]*t[2]*t[16]*t[20]-3.835527988457806E-4*t[1]*t[2]*t[16]*t[20]);
    d3fdragagb = 0.0;
    d3fdragaab = 0.0;
    d3fdragbgb = 0.0;
    d3fdragbab = 0.0;
    d3fdraabab = 0.0;
    d3fdrbrbrb = 0.5*(6.1270251210506771E-6*t[1]*t[34]*pow(gradb,6.0)/pow(rhob,9.666666666666666)+0.027354028384519*t[12]*t[33]/pow(rhob,1.666666666666667)+7.9046813011671027E-4*t[6]*t[8]*t[10]*t[36]-0.0039523406505836*t[33]*t[8]*t[10]*t[36]+0.027666384554085*t[21]*t[8]*t[10]*t[36]+5.6822636866041565E-5*t[6]*t[22]*t[23]*t[35]-5.1140373179437413E-4*t[21]*t[22]*t[23]*t[35]-3.9775845806229101E-4*t[1]*t[22]*t[23]*t[35]);
    d3fdrbrbga = 0.0;
    d3fdrbrbgb = 0.5*(-4.5952688407880086E-6*t[1]*t[34]*pow(gradb,5.0)/pow(rhob,8.666666666666666)-5.9285109758753281E-4*t[6]*gradb*t[10]*t[25]+0.0023714043903501*t[33]*gradb*t[10]*t[25]-0.016599830732451*t[21]*gradb*t[10]*t[25]+1.2785093294859353E-4*t[6]*t[26]*t[23]*t[24]+2.9831884354671826E-4*t[1]*t[26]*t[23]*t[24]);
    d3fdrbrbab = 0.0;
    d3fdrbgaga = 0.0;
    d3fdrbgagb = 0.0;
    d3fdrbgaab = 0.0;
    d3fdrbgbgb = 0.5*(3.4464516305910063E-6*t[1]*t[22]*t[34]/pow(rhob,7.666666666666667)+3.1962733237148383E-5*t[6]*t[8]*t[23]*t[27]-1.2785093294859353E-4*t[21]*t[8]*t[23]*t[27]-3.835527988457806E-4*t[1]*t[8]*t[23]*t[27]-0.0017785532927626*t[6]*t[10]*t[11]+0.0071142131710504*t[21]*t[10]*t[11]);
    d3fdrbgbab = 0.0;
    d3fdrbabab = 0.0;
    d3fdgagaga = -1.2924193614716277E-6*t[1]*t[19]*t[30]/pow(rhoa,6.666666666666667)+4.7944099855722582E-5*t[21]*grada*t[16]*t[28]+4.7944099855722582E-5*t[1]*grada*t[16]*t[28];
    d3fdgagagb = 0.0;
    d3fdgagaab = 0.0;
    d3fdgagbgb = 0.0;
    d3fdgagbab = 0.0;
    d3fdgaabab = 0.0;
    d3fdgbgbgb = -1.2924193614716277E-6*t[1]*t[26]*t[34]/pow(rhob,6.666666666666667)+4.7944099855722582E-5*t[21]*gradb*t[23]*t[29]+4.7944099855722582E-5*t[1]*gradb*t[23]*t[29];
    d3fdgbgbab = 0.0;
    d3fdgbabab = 0.0;
    d3fdababab = 0.0;
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
    ds->df0011 += factor*d2fdgagb;
    ds->df00101+= factor*d2fdgaab;
    ds->df0002 += factor*d2fdgbgb;
    ds->df00011+= factor*d2fdgbab;
    ds->df00002+= factor*d2fdabab;
    ds->df3000 += factor*d3fdrarara;
    ds->df2100 += factor*d3fdrararb;
    ds->df2010 += factor*d3fdraraga;
    ds->df2001 += factor*d3fdraragb;
    ds->df20001 += factor*d3fdraraab;
    ds->df1200 += factor*d3fdrarbrb;
    ds->df1110 += factor*d3fdrarbga;
    ds->df1101 += factor*d3fdrarbgb;
    ds->df11001 += factor*d3fdrarbab;
    ds->df1020 += factor*d3fdragaga;
    ds->df1011 += factor*d3fdragagb;
    ds->df10101+= factor*d3fdragaab;
    ds->df1002 += factor*d3fdragbgb;
    ds->df10011+= factor*d3fdragbab;
    ds->df10002+= factor*d3fdraabab;
    ds->df0300 += factor*d3fdrbrbrb;
    ds->df0210 += factor*d3fdrbrbga;
    ds->df0201 += factor*d3fdrbrbgb;
    ds->df02001 += factor*d3fdrbrbab;
    ds->df0120 += factor*d3fdrbgaga;
    ds->df0111 += factor*d3fdrbgagb;
    ds->df01101+= factor*d3fdrbgaab;
    ds->df0102 += factor*d3fdrbgbgb;
    ds->df01011+= factor*d3fdrbgbab;
    ds->df01002+= factor*d3fdrbgbab;
    ds->df0030 += factor*d3fdgagaga;
    ds->df0021 += factor*d3fdgagagb;
    ds->df00201+= factor*d3fdgagaab;
    ds->df0012 += factor*d3fdgagbgb;
    ds->df00111+= factor*d3fdgagbab;
    ds->df00102+= factor*d3fdgaabab;
    ds->df0003 += factor*d3fdgbgbgb;
    ds->df00021+= factor*d3fdgbgbab;
    ds->df00012+= factor*d3fdgbabab;
    ds->df00003+= factor*d3fdababab;
}

static void
pbex_fourth(FunFourthFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real t[60];
    real dfdra, dfdrb, dfdga, dfdgb, dfdab;
    real d2fdrara , d2fdrarb , d2fdraga , d2fdragb ;
    real d2fdraab , d2fdrbrb , d2fdrbga , d2fdrbgb ;
    real d2fdrbab , d2fdgaga , d2fdgagb , d2fdgaab ;
    real d2fdgbgb , d2fdgbab , d2fdabab ;
    real d3fdrarara , d3fdrararb , d3fdraraga , d3fdraragb ;
    real d3fdraraab , d3fdrarbrb , d3fdrarbga , d3fdrarbgb ;
    real d3fdrarbab , d3fdragaga , d3fdragagb , d3fdragaab ;
    real d3fdragbgb , d3fdragbab , d3fdraabab , d3fdrbrbrb ;
    real d3fdrbrbga , d3fdrbrbgb , d3fdrbrbab , d3fdrbgaga ;
    real d3fdrbgagb , d3fdrbgaab , d3fdrbgbgb , d3fdrbgbab ;
    real d3fdrbabab , d3fdgagaga , d3fdgagagb , d3fdgagaab ;
    real d3fdgagbgb , d3fdgagbab , d3fdgaabab , d3fdgbgbgb ;
    real d3fdgbgbab , d3fdgbabab , d3fdababab ;
    real d4fdrararara , d4fdrarararb , d4fdrararaga , d4fdrararagb ;
    real d4fdrararaab , d4fdrararbrb , d4fdrararbga , d4fdrararbgb ;
    real d4fdrararbab , d4fdraragaga , d4fdraragagb , d4fdraragaab ;
    real d4fdraragbgb , d4fdraragbab , d4fdraraabab , d4fdrarbrbrb ;
    real d4fdrarbrbga , d4fdrarbrbgb , d4fdrarbrbab , d4fdrarbgaga ;
    real d4fdrarbgagb , d4fdrarbgaab , d4fdrarbgbgb , d4fdrarbgbab ;
    real d4fdrarbabab , d4fdragagaga , d4fdragagagb , d4fdragagaab ;
    real d4fdragagbgb , d4fdragagbab , d4fdragaabab , d4fdragbgbgb ;
    real d4fdragbgbab , d4fdragbabab , d4fdraababab , d4fdrbrbrbrb ;
    real d4fdrbrbrbga , d4fdrbrbrbgb , d4fdrbrbrbab , d4fdrbrbgaga ;
    real d4fdrbrbgagb , d4fdrbrbgaab , d4fdrbrbgbgb , d4fdrbrbgbab ;
    real d4fdrbrbabab , d4fdrbgagaga , d4fdrbgagagb , d4fdrbgagaab ;
    real d4fdrbgagbgb , d4fdrbgagbab , d4fdrbgaabab , d4fdrbgbgbgb ;
    real d4fdrbgbgbab , d4fdrbgbabab , d4fdrbababab , d4fdgagagaga ;
    real d4fdgagagagb , d4fdgagagaab , d4fdgagagbgb , d4fdgagagbab ;
    real d4fdgagaabab , d4fdgagbgbgb , d4fdgagbgbab , d4fdgagbabab ;
    real d4fdgaababab , d4fdgbgbgbgb , d4fdgbgbgbab , d4fdgbgbabab ;
    real d4fdgbababab , d4fdabababab ;
    real rhoa = dp->rhoa;
    real rhob = dp->rhob;
    real grada = dp->grada;
    real gradb = dp->gradb;
    real gradab = dp->gradab;

    t[1] = pow(2.0,1.333333333333333);
    t[2] = pow(grada,2.0);
    t[3] = 1/pow(rhoa,2.666666666666667);
    t[4] = 0.0044927994802311*t[2]*t[3]+1.0;
    t[5] = 1/pow(t[4],2.0);
    t[6] = 1/pow(rhoa,2.333333333333334);
    t[7] = pow(2.0,3.333333333333334);
    t[8] = 1.804-0.804/t[4];
    t[9] = pow(gradb,2.0);
    t[10] = 1/pow(rhob,2.666666666666667);
    t[11] = 0.0044927994802311*t[9]*t[10]+1.0;
    t[12] = 1/pow(t[11],2.0);
    t[13] = 1/pow(rhob,2.333333333333334);
    t[14] = 1.804-0.804/t[11];
    t[15] = 1/pow(rhoa,1.333333333333333);
    t[16] = 1/pow(rhob,1.333333333333333);
    t[17] = pow(grada,4.0);
    t[18] = 1/pow(t[4],3.0);
    t[19] = 1/pow(rhoa,6.0);
    t[20] = 1/pow(rhoa,3.333333333333334);
    t[21] = pow(grada,3.0);
    t[22] = 1/pow(rhoa,5.0);
    t[23] = pow(2.0,2.333333333333334);
    t[24] = pow(gradb,4.0);
    t[25] = 1/pow(t[11],3.0);
    t[26] = 1/pow(rhob,6.0);
    t[27] = 1/pow(rhob,3.333333333333334);
    t[28] = pow(gradb,3.0);
    t[29] = 1/pow(rhob,5.0);
    t[30] = 1/pow(rhoa,4.0);
    t[31] = 1/pow(rhob,4.0);
    t[32] = pow(grada,6.0);
    t[33] = 1/pow(t[4],4.0);
    t[34] = 1/pow(rhoa,9.666666666666666);
    t[35] = 1/pow(rhoa,7.0);
    t[36] = 1/pow(rhoa,4.333333333333333);
    t[37] = pow(2.0,4.333333333333333);
    t[38] = pow(grada,5.0);
    t[39] = 1/pow(rhoa,8.666666666666666);
    t[40] = 1/pow(rhoa,7.666666666666667);
    t[41] = pow(gradb,6.0);
    t[42] = 1/pow(t[11],4.0);
    t[43] = 1/pow(rhob,9.666666666666666);
    t[44] = 1/pow(rhob,7.0);
    t[45] = 1/pow(rhob,4.333333333333333);
    t[46] = pow(gradb,5.0);
    t[47] = 1/pow(rhob,8.666666666666666);
    t[48] = 1/pow(rhob,7.666666666666667);
    t[49] = 1/pow(rhoa,6.666666666666667);
    t[50] = 1/pow(rhob,6.666666666666667);
    t[51] = 1/pow(t[4],5.0);
    t[52] = 1/pow(rhoa,10.66666666666667);
    t[53] = 1/pow(rhoa,8.0);
    t[54] = 1/pow(rhoa,5.333333333333333);
    t[55] = pow(2.0,5.333333333333333);
    t[56] = 1/pow(t[11],5.0);
    t[57] = 1/pow(rhob,10.66666666666667);
    t[58] = 1/pow(rhob,8.0);
    t[59] = 1/pow(rhob,5.333333333333333);
    dfdra = 0.5*(0.0071142131710504*t[1]*t[2]*t[5]*t[6]-0.24618625546067*t[7]*t[8]*pow(rhoa,0.33333333333333));
    dfdrb = 0.5*(0.0071142131710504*t[1]*t[9]*t[12]*t[13]-0.24618625546067*t[14]*t[7]*pow(rhob,0.33333333333333));
    dfdga = -0.0026678299391439*t[1]*grada*t[5]*t[15];
    dfdgb = -0.0026678299391439*t[1]*gradb*t[12]*t[16];
    dfdab = 0.0;
    d2fdrara = 0.5*(-0.082062085153558*t[7]*t[8]/pow(rhoa,0.66666666666667)+0.0023714043903501*t[7]*t[2]*t[5]*t[20]-0.016599830732451*t[1]*t[2]*t[5]*t[20]+1.704679105981247E-4*t[1]*t[17]*t[18]*t[19]);
    d2fdrarb = 0.0;
    d2fdraga = 0.5*(-0.0017785532927626*t[7]*grada*t[5]*t[6]+0.0071142131710504*t[23]*grada*t[5]*t[6]-1.2785093294859353E-4*t[1]*t[21]*t[18]*t[22]);
    d2fdragb = 0.0;
    d2fdraab = 0.0;
    d2fdrbrb = 0.5*(-0.082062085153558*t[14]*t[7]/pow(rhob,0.66666666666667)+0.0023714043903501*t[7]*t[9]*t[12]*t[27]-0.016599830732451*t[1]*t[9]*t[12]*t[27]+1.704679105981247E-4*t[1]*t[24]*t[25]*t[26]);
    d2fdrbga = 0.0;
    d2fdrbgb = 0.5*(-0.0017785532927626*t[7]*gradb*t[12]*t[13]+0.0071142131710504*t[23]*gradb*t[12]*t[13]-1.2785093294859353E-4*t[1]*t[28]*t[25]*t[29]);
    d2fdrbab = 0.0;
    d2fdgaga = 4.7944099855722582E-5*t[1]*t[2]*t[18]*t[30]-0.0026678299391439*t[1]*t[5]*t[15];
    d2fdgagb = 0.0;
    d2fdgaab = 0.0;
    d2fdgbgb = 4.7944099855722582E-5*t[1]*t[9]*t[25]*t[31]-0.0026678299391439*t[1]*t[12]*t[16];
    d2fdgbab = 0.0;
    d2fdabab = 0.0;
    d3fdrarara= 0.5*(0.027354028384519*t[37]*t[8]/pow(rhoa,1.666666666666667)+7.9046813011671027E-4*t[7]*t[2]*t[5]*t[36]-0.0039523406505836*t[37]*t[2]*t[5]*t[36]+0.027666384554085*t[23]*t[2]*t[5]*t[36]+5.6822636866041565E-5*t[7]*t[17]*t[18]*t[35]-5.1140373179437413E-4*t[23]*t[17]*t[18]*t[35]-3.9775845806229101E-4*t[1]*t[17]*t[18]*t[35]+6.1270251210506771E-6*t[1]*t[32]*t[33]*t[34]);
    d3fdrararb = 0.0;
    d3fdraraga = 0.5*(0.0023714043903501*t[37]*grada*t[5]*t[20]-5.9285109758753281E-4*t[7]*grada*t[5]*t[20]-0.016599830732451*t[23]*grada*t[5]*t[20]+1.2785093294859353E-4*t[7]*t[21]*t[18]*t[19]+2.9831884354671826E-4*t[1]*t[21]*t[18]*t[19]-4.5952688407880086E-6*t[1]*t[38]*t[33]*t[39]);
    d3fdraragb = 0.0;
    d3fdraraab = 0.0;
    d3fdrarbrb = 0.0;
    d3fdrarbga = 0.0;
    d3fdrarbgb = 0.0;
    d3fdrarbab = 0.0;
    d3fdragaga = 0.5*(-0.0017785532927626*t[7]*t[5]*t[6]+0.0071142131710504*t[23]*t[5]*t[6]+3.1962733237148383E-5*t[7]*t[2]*t[18]*t[22]-1.2785093294859353E-4*t[23]*t[2]*t[18]*t[22]-3.835527988457806E-4*t[1]*t[2]*t[18]*t[22]+3.4464516305910063E-6*t[1]*t[17]*t[33]*t[40]);
    d3fdragagb = 0.0;
    d3fdragaab = 0.0;
    d3fdragbgb = 0.0;
    d3fdragbab = 0.0;
    d3fdraabab = 0.0;
    d3fdrbrbrb = 0.5*(0.027354028384519*t[14]*t[37]/pow(rhob,1.666666666666667)+7.9046813011671027E-4*t[7]*t[9]*t[12]*t[45]-0.0039523406505836*t[37]*t[9]*t[12]*t[45]+0.027666384554085*t[23]*t[9]*t[12]*t[45]+5.6822636866041565E-5*t[7]*t[24]*t[25]*t[44]-5.1140373179437413E-4*t[23]*t[24]*t[25]*t[44]-3.9775845806229101E-4*t[1]*t[24]*t[25]*t[44]+6.1270251210506771E-6*t[1]*t[41]*t[42]*t[43]);
    d3fdrbrbga = 0.0;
    d3fdrbrbgb = 0.5*(0.0023714043903501*t[37]*gradb*t[12]*t[27]-5.9285109758753281E-4*t[7]*gradb*t[12]*t[27]-0.016599830732451*t[23]*gradb*t[12]*t[27]+1.2785093294859353E-4*t[7]*t[28]*t[25]*t[26]+2.9831884354671826E-4*t[1]*t[28]*t[25]*t[26]-4.5952688407880086E-6*t[1]*t[46]*t[42]*t[47]);
    d3fdrbrbab = 0.0;
    d3fdrbgaga = 0.0;
    d3fdrbgagb = 0.0;
    d3fdrbgaab = 0.0;
    d3fdrbgbgb= 0.5*(-0.0017785532927626*t[7]*t[12]*t[13]+0.0071142131710504*t[23]*t[12]*t[13]+3.1962733237148383E-5*t[7]*t[9]*t[25]*t[29]-1.2785093294859353E-4*t[23]*t[9]*t[25]*t[29]-3.835527988457806E-4*t[1]*t[9]*t[25]*t[29]+3.4464516305910063E-6*t[1]*t[24]*t[42]*t[48]);
    d3fdrbgbab = 0.0;
    d3fdrbabab = 0.0;
    d3fdgagaga = 4.7944099855722582E-5*t[23]*grada*t[18]*t[30]+4.7944099855722582E-5*t[1]*grada*t[18]*t[30]-1.2924193614716277E-6*t[1]*t[21]*t[33]*t[49];
    d3fdgagagb = 0.0;
    d3fdgagaab = 0.0;
    d3fdgagbgb = 0.0;
    d3fdgagbab = 0.0;
    d3fdgaabab = 0.0;
    d3fdgbgbgb = 4.7944099855722582E-5*t[23]*gradb*t[25]*t[31]+4.7944099855722582E-5*t[1]*gradb*t[25]*t[31]-1.2924193614716277E-6*t[1]*t[28]*t[42]*t[50];
    d3fdgbgbab = 0.0;
    d3fdgbabab = 0.0;
    d3fdababab = 0.0;
    d4fdrararara = 0.5*(2.9362661631167268E-7*t[1]*t[51]*pow(grada,8.0)/pow(rhoa,13.33333333333333)-0.0034253618971724*t[7]*t[2]*t[5]*t[54]+0.016863320109156*t[37]*t[2]*t[5]*t[54]-0.11988766640103*t[23]*t[2]*t[5]*t[54]-3.7881757910694376E-4*t[7]*t[17]*t[18]*t[53]-9.4704394776735953E-5*t[37]*t[17]*t[18]*t[53]+0.0042427568859978*t[23]*t[17]*t[18]*t[53]+0.002784309206436*t[1]*t[17]*t[18]*t[53]+2.0423417070168925E-6*t[7]*t[32]*t[33]*t[52]-1.8381075363152035E-5*t[23]*t[32]*t[33]*t[52]-7.3524301452608125E-5*t[1]*t[32]*t[33]*t[52]-0.045590047307532*t[37]*t[8]*t[3]);
    d4fdrarararb = 0.0;
    d4fdrararaga = 0.5*(-2.2021996223375453E-7*t[1]*t[51]*pow(grada,7.0)/pow(rhoa,12.33333333333333)+0.027666384554085*t[7]*grada*t[5]*t[36]-0.0039523406505836*t[55]*grada*t[5]*t[36]+9.8808516264588774E-4*t[37]*grada*t[5]*t[36]-4.1196411727880143E-4*t[7]*t[21]*t[18]*t[35]+5.6822636866041565E-5*t[55]*t[21]*t[18]*t[35]-4.4037543571182216E-4*t[37]*t[21]*t[18]*t[35]-4.9719807257786388E-4*t[23]*t[21]*t[18]*t[35]-1.5317562802626695E-6*t[7]*t[38]*t[33]*t[34]+3.216688188551606E-5*t[23]*t[38]*t[33]*t[34]+1.0722293961838687E-5*t[1]*t[38]*t[33]*t[34]);
    d4fdrararagb = 0.0;
    d4fdrararaab = 0.0;
    d4fdrararbrb = 0.0;
    d4fdrararbga = 0.0;
    d4fdrararbgb = 0.0;
    d4fdrararbab = 0.0;
    d4fdraragaga = 0.5*(1.6516497167531593E-7*t[1]*t[32]*t[51]/pow(rhoa,11.33333333333333)-3.4464516305910063E-6*t[7]*t[17]*t[33]*t[39]-3.1018064675319062E-5*t[1]*t[17]*t[33]*t[39]-5.9285109758753281E-4*t[7]*t[5]*t[20]+0.0023714043903501*t[37]*t[5]*t[20]-0.016599830732451*t[23]*t[5]*t[20]+3.9420704325816337E-4*t[7]*t[2]*t[18]*t[19]-4.2616977649531175E-5*t[37]*t[2]*t[18]*t[19]+2.9831884354671826E-4*t[23]*t[2]*t[18]*t[19]+8.9495653064015478E-4*t[1]*t[2]*t[18]*t[19]);
    d4fdraragagb = 0.0;
    d4fdraragaab = 0.0;
    d4fdraragbgb = 0.0;
    d4fdraragbab = 0.0;
    d4fdraraabab = 0.0;
    d4fdrarbrbrb = 0.0;
    d4fdrarbrbga = 0.0;
    d4fdrarbrbgb = 0.0;
    d4fdrarbrbab = 0.0;
    d4fdrarbgaga = 0.0;
    d4fdrarbgagb = 0.0;
    d4fdrarbgaab = 0.0;
    d4fdrarbgbgb = 0.0;
    d4fdrarbgbab = 0.0;
    d4fdrarbabab = 0.0;
    d4fdragagaga = 0.5*(-1.2387372875648693E-7*t[1]*t[38]*t[51]/pow(rhoa,10.33333333333333)+2.5848387229432549E-6*t[7]*t[21]*t[33]*t[40]+3.4464516305910063E-6*t[23]*t[21]*t[33]*t[40]+1.033935489177302E-5*t[1]*t[21]*t[33]*t[40]-9.588819971144515E-5*t[7]*grada*t[18]*t[22]+3.1962733237148383E-5*t[37]*grada*t[18]*t[22]-5.1140373179437413E-4*t[23]*grada*t[18]*t[22]);
    d4fdragagagb = 0.0;
    d4fdragagaab = 0.0;
    d4fdragagbgb = 0.0;
    d4fdragagbab = 0.0;
    d4fdragaabab = 0.0;
    d4fdragbgbgb = 0.0;
    d4fdragbgbab = 0.0;
    d4fdragbabab = 0.0;
    d4fdraababab = 0.0;
    d4fdrbrbrbrb= 0.5*(2.9362661631167268E-7*t[1]*t[56]*pow(gradb,8.0)/pow(rhob,13.33333333333333)-0.0034253618971724*t[7]*t[9]*t[12]*t[59]+0.016863320109156*t[37]*t[9]*t[12]*t[59]-0.11988766640103*t[23]*t[9]*t[12]*t[59]-3.7881757910694376E-4*t[7]*t[24]*t[25]*t[58]-9.4704394776735953E-5*t[37]*t[24]*t[25]*t[58]+0.0042427568859978*t[23]*t[24]*t[25]*t[58]+0.002784309206436*t[1]*t[24]*t[25]*t[58]+2.0423417070168925E-6*t[7]*t[41]*t[42]*t[57]-1.8381075363152035E-5*t[23]*t[41]*t[42]*t[57]-7.3524301452608125E-5*t[1]*t[41]*t[42]*t[57]-0.045590047307532*t[37]*t[14]*t[10]);
    d4fdrbrbrbga = 0.0;
    d4fdrbrbrbgb = 0.5*(-2.2021996223375453E-7*t[1]*t[56]*pow(gradb,7.0)/pow(rhob,12.33333333333333)+0.027666384554085*t[7]*gradb*t[12]*t[45]-0.0039523406505836*t[55]*gradb*t[12]*t[45]+9.8808516264588774E-4*t[37]*gradb*t[12]*t[45]-4.1196411727880143E-4*t[7]*t[28]*t[25]*t[44]+5.6822636866041565E-5*t[55]*t[28]*t[25]*t[44]-4.4037543571182216E-4*t[37]*t[28]*t[25]*t[44]-4.9719807257786388E-4*t[23]*t[28]*t[25]*t[44]-1.5317562802626695E-6*t[7]*t[46]*t[42]*t[43]+3.216688188551606E-5*t[23]*t[46]*t[42]*t[43]+1.0722293961838687E-5*t[1]*t[46]*t[42]*t[43]);
    d4fdrbrbrbab = 0.0;
    d4fdrbrbgaga = 0.0;
    d4fdrbrbgagb = 0.0;
    d4fdrbrbgaab = 0.0;
    d4fdrbrbgbgb = 0.5*(1.6516497167531593E-7*t[1]*t[41]*t[56]/pow(rhob,11.33333333333333)-3.4464516305910063E-6*t[7]*t[24]*t[42]*t[47]-3.1018064675319062E-5*t[1]*t[24]*t[42]*t[47]-5.9285109758753281E-4*t[7]*t[12]*t[27]+0.0023714043903501*t[37]*t[12]*t[27]-0.016599830732451*t[23]*t[12]*t[27]+3.9420704325816337E-4*t[7]*t[9]*t[25]*t[26]-4.2616977649531175E-5*t[37]*t[9]*t[25]*t[26]+2.9831884354671826E-4*t[23]*t[9]*t[25]*t[26]+8.9495653064015478E-4*t[1]*t[9]*t[25]*t[26]);
    d4fdrbrbgbab = 0.0;
    d4fdrbrbabab = 0.0;
    d4fdrbgagaga = 0.0;
    d4fdrbgagagb = 0.0;
    d4fdrbgagaab = 0.0;
    d4fdrbgagbgb = 0.0;
    d4fdrbgagbab = 0.0;
    d4fdrbgaabab = 0.0;
    d4fdrbgbgbgb = 0.5*(-1.2387372875648693E-7*t[1]*t[46]*t[56]/pow(rhob,10.33333333333333)+2.5848387229432549E-6*t[7]*t[28]*t[42]*t[48]+3.4464516305910063E-6*t[23]*t[28]*t[42]*t[48]+1.033935489177302E-5*t[1]*t[28]*t[42]*t[48]-9.588819971144515E-5*t[7]*gradb*t[25]*t[29]+3.1962733237148383E-5*t[37]*gradb*t[25]*t[29]-5.1140373179437413E-4*t[23]*gradb*t[25]*t[29]);
    d4fdrbgbgbab = 0.0;
    d4fdrbgbabab = 0.0;
    d4fdrbababab = 0.0;
    d4fdgagagaga = 4.6452648283682616E-8*t[1]*t[17]*t[51]/pow(rhoa,9.333333333333334)-1.2924193614716277E-6*t[23]*t[2]*t[33]*t[49]-5.1696774458865107E-6*t[1]*t[2]*t[33]*t[49]+4.7944099855722582E-5*t[23]*t[18]*t[30]+4.7944099855722582E-5*t[1]*t[18]*t[30];
    d4fdgagagagb = 0.0;
    d4fdgagagaab = 0.0;
    d4fdgagagbgb = 0.0;
    d4fdgagagbab = 0.0;
    d4fdgagaabab = 0.0;
    d4fdgagbgbgb = 0.0;
    d4fdgagbgbab = 0.0;
    d4fdgagbabab = 0.0;
    d4fdgaababab = 0.0;
    d4fdgbgbgbgb = 4.6452648283682616E-8*t[1]*t[24]*t[56]/pow(rhob,9.333333333333334)-1.2924193614716277E-6*t[23]*t[9]*t[42]*t[50]-5.1696774458865107E-6*t[1]*t[9]*t[42]*t[50]+4.7944099855722582E-5*t[23]*t[25]*t[31]+4.7944099855722582E-5*t[1]*t[25]*t[31];
    d4fdgbgbgbab = 0.0;
    d4fdgbgbabab = 0.0;
    d4fdgbababab = 0.0;
    d4fdabababab = 0.0;
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
    ds->df0011 += factor*d2fdgagb;
    ds->df00101 += factor*d2fdgaab;
    ds->df0002 += factor*d2fdgbgb;
    ds->df00011 += factor*d2fdgbab;
    ds->df00002 += factor*d2fdabab;
    ds->df3000 += factor*d3fdrarara;
    ds->df2100 += factor*d3fdrararb;
    ds->df2010 += factor*d3fdraraga;
    ds->df2001 += factor*d3fdraragb;
    ds->df20001 += factor*d3fdraraab;
    ds->df1200 += factor*d3fdrarbrb;
    ds->df1110 += factor*d3fdrarbga;
    ds->df1101 += factor*d3fdrarbgb;
    ds->df11001 += factor*d3fdrarbab;
    ds->df1020 += factor*d3fdragaga;
    ds->df1011 += factor*d3fdragagb;
    ds->df10101 += factor*d3fdragaab;
    ds->df1002 += factor*d3fdragbgb;
    ds->df10011 += factor*d3fdragbab;
    ds->df10002 += factor*d3fdraabab;
    ds->df0300 += factor*d3fdrbrbrb;
    ds->df0210 += factor*d3fdrbrbga;
    ds->df0201 += factor*d3fdrbrbgb;
    ds->df02001 += factor*d3fdrbrbab;
    ds->df0120 += factor*d3fdrbgaga;
    ds->df0111 += factor*d3fdrbgagb;
    ds->df01101 += factor*d3fdrbgaab;
    ds->df0102 += factor*d3fdrbgbgb;
    ds->df01011 += factor*d3fdrbgbab;
    ds->df01002 += factor*d3fdrbabab;
    ds->df0030 += factor*d3fdgagaga;
    ds->df0021 += factor*d3fdgagagb;
    ds->df00201 += factor*d3fdgagaab;
    ds->df0012 += factor*d3fdgagbgb;
    ds->df00111 += factor*d3fdgagbab;
    ds->df00102 += factor*d3fdgaabab;
    ds->df0003 += factor*d3fdgbgbgb;
    ds->df00021 += factor*d3fdgbgbab;
    ds->df00012 += factor*d3fdgbabab;
    ds->df00003 += factor*d3fdababab;
    ds->df4000 += factor*d4fdrararara;
    ds->df3100 += factor*d4fdrarararb;
    ds->df3010 += factor*d4fdrararaga;
    ds->df3001 += factor*d4fdrararagb;
    ds->df30001 += factor*d4fdrararaab;
    ds->df2200 += factor*d4fdrararbrb;
    ds->df2110 += factor*d4fdrararbga;
    ds->df2101 += factor*d4fdrararbgb;
    ds->df21001 += factor*d4fdrararbab;
    ds->df2020 += factor*d4fdraragaga;
    ds->df2011 += factor*d4fdraragagb;
    ds->df20101 += factor*d4fdraragaab;
    ds->df2002 += factor*d4fdraragbgb;
    ds->df20011 += factor*d4fdraragbab;
    ds->df20002 += factor*d4fdraraabab;
    ds->df1300 += factor*d4fdrarbrbrb;
    ds->df1210 += factor*d4fdrarbrbga;
    ds->df1201 += factor*d4fdrarbrbgb;
    ds->df12001 += factor*d4fdrarbrbab;
    ds->df1120 += factor*d4fdrarbgaga;
    ds->df1111 += factor*d4fdrarbgagb;
    ds->df11101 += factor*d4fdrarbgaab;
    ds->df1102 += factor*d4fdrarbgbgb;
    ds->df11011 += factor*d4fdrarbgbab;
    ds->df11002 += factor*d4fdrarbabab;
    ds->df1030 += factor*d4fdragagaga;
    ds->df1021 += factor*d4fdragagagb;
    ds->df10201 += factor*d4fdragagaab;
    ds->df1012 += factor*d4fdragagbgb;
    ds->df10111 += factor*d4fdragagbab;
    ds->df10102 += factor*d4fdragaabab;
    ds->df1003 += factor*d4fdragbgbgb;
    ds->df10021 += factor*d4fdragbgbab;
    ds->df10012 += factor*d4fdragbabab;
    ds->df10003 += factor*d4fdraababab;
    ds->df0400 += factor*d4fdrbrbrbrb;
    ds->df0310 += factor*d4fdrbrbrbga;
    ds->df0301 += factor*d4fdrbrbrbgb;
    ds->df03001 += factor*d4fdrbrbrbab;
    ds->df0220 += factor*d4fdrbrbgaga;
    ds->df0211 += factor*d4fdrbrbgagb;
    ds->df02101 += factor*d4fdrbrbgaab;
    ds->df0202 += factor*d4fdrbrbgbgb;
    ds->df02011 += factor*d4fdrbrbgbab;
    ds->df02002 += factor*d4fdrbrbabab;
    ds->df0130 += factor*d4fdrbgagaga;
    ds->df0121 += factor*d4fdrbgagagb;
    ds->df01201 += factor*d4fdrbgagaab;
    ds->df0112 += factor*d4fdrbgagbgb;
    ds->df01111 += factor*d4fdrbgagbab;
    ds->df01102 += factor*d4fdrbgaabab;
    ds->df0103 += factor*d4fdrbgbgbgb;
    ds->df01021 += factor*d4fdrbgbgbab;
    ds->df01012 += factor*d4fdrbgbabab;
    ds->df01003 += factor*d4fdrbababab;
    ds->df0040 += factor*d4fdgagagaga;
    ds->df0031 += factor*d4fdgagagagb;
    ds->df00301 += factor*d4fdgagagaab;
    ds->df0022 += factor*d4fdgagagbgb;
    ds->df00211 += factor*d4fdgagagbab;
    ds->df00202 += factor*d4fdgagaabab;
    ds->df0013 += factor*d4fdgagbgbgb;
    ds->df00121 += factor*d4fdgagbgbab;
    ds->df00112 += factor*d4fdgagbabab;
    ds->df00103 += factor*d4fdgaababab;
    ds->df0004 += factor*d4fdgbgbgbgb;
    ds->df00031 += factor*d4fdgbgbgbab;
    ds->df00022 += factor*d4fdgbgbabab;
    ds->df00013 += factor*d4fdgbababab;
    ds->df00004 += factor*d4fdabababab;
}
