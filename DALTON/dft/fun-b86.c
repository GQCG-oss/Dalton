/* Automatically generated functional code: becke86
   Maxima input:
    >> BETA:0.00787;
    >> GAMMA:0.004;
    >> PI:3.14159265358979312;
    >> PREF:3/4*(6/PI)^(1/3);
    >> ra43:rhoa^(4/3);
    >> rb43:rhob^(4/3);
    >> denoma:1+GAMMA*xa*xa;
    >> denomb:1+GAMMA*xb*xb;
    >> 
    >> Exa: -PREF*ra43*(1 + BETA*xa*xa)/denoma;
    >> Exb: -PREF*rb43*(1 + BETA*xb*xb)/denomb;
    >> 
    >> K(rhoa,grada,rhob,gradb,gradab) := Exa + Exb;
*/

// add "extern Functional becke86Functional;" to 'functionals.h'
// add "&becke86Functional," to 'functionals.c'
// add "fun-becke86.c" to 'Makefile.in'

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
static integer becke86_isgga(void) {return 1;}
static integer becke86_read(const char* conf_line);
static real becke86_energy(const FunDensProp* dp);
static void becke86_first(FunFirstFuncDrv *ds, real factor, 
                           const FunDensProp* dp);
static void becke86_second(FunSecondFuncDrv *ds, real factor,
                            const FunDensProp* dp);
static void becke86_third(FunThirdFuncDrv *ds, real factor,
                           const FunDensProp* dp);

static void becke86_fourth(FunFourthFuncDrv *ds, real factor,
                           const FunDensProp* dp);

//static integer fun_true(void) { return 1; }
Functional B86xFunctional = {
  "B86x",
  becke86_isgga,
  3,
  becke86_read,
  NULL,
  becke86_energy,
  becke86_first,
  becke86_second,
  becke86_third,
  becke86_fourth
};

/* IMPLEMENTATION PART */
static integer
becke86_read(const char* conf_line)
{
    fun_set_hf_weight(0);
    return 1;
}


static real
becke86_energy(const FunDensProp* dp)
{
    real t[5],zk;
    real rhoa = dp->rhoa;
    real rhob = dp->rhob;
    real grada = dp->grada;
    real gradb = dp->gradb;
    real gradab = dp->gradab;

    t[1] = pow(grada,2.0);
    t[2] = 1/pow(rhoa,2.666666666666667);
    t[3] = pow(gradb,2.0);
    t[4] = 1/pow(rhob,2.666666666666667);
    zk = -0.9305257363491*(0.00787*t[3]*t[4]+1.0)*pow(rhob,1.333333333333333)/(0.004*t[3]*t[4]+1.0)-0.9305257363491*(0.00787*t[1]*t[2]+1.0)*pow(rhoa,1.333333333333333)/(0.004*t[1]*t[2]+1.0);
    return zk;
}

static void
becke86_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real t[17];
    real dfdra, dfdrb, dfdga, dfdgb, dfdab;
    real rhoa = dp->rhoa;
    real rhob = dp->rhob;
    real grada = dp->grada;
    real gradb = dp->gradb;
    real gradab = dp->gradab;

    t[1] = pow(grada,2.0);
    t[2] = 1/pow(rhoa,2.666666666666667);
    t[3] = 0.004*t[1]*t[2]+1.0;
    t[4] = 1/t[3];
    t[5] = 1/pow(rhoa,2.333333333333333);
    t[6] = 1/pow(t[3],2.0);
    t[7] = 0.00787*t[1]*t[2]+1.0;
    t[8] = pow(gradb,2.0);
    t[9] = 1/pow(rhob,2.666666666666667);
    t[10] = 0.004*t[8]*t[9]+1.0;
    t[11] = 1/t[10];
    t[12] = 1/pow(rhob,2.333333333333333);
    t[13] = 1/pow(t[10],2.0);
    t[14] = 0.00787*t[8]*t[9]+1.0;
    t[15] = 1/pow(rhoa,1.333333333333333);
    t[16] = 1/pow(rhob,1.333333333333333);
    dfdra = -1.2407009817988*t[4]*t[7]*pow(rhoa,.3333333333333333)-0.0099256078543904*t[1]*t[6]*t[7]*t[5]+.01952863345351311*t[1]*t[4]*t[5];
    dfdrb = -1.2407009817988*t[11]*t[14]*pow(rhob,.3333333333333333)-0.0099256078543904*t[8]*t[13]*t[14]*t[12]+.01952863345351311*t[8]*t[11]*t[12];
    dfdga = 0.0074442058907928*grada*t[6]*t[7]*t[15]-.01464647509013483*grada*t[4]*t[15];
    dfdgb = 0.0074442058907928*gradb*t[13]*t[14]*t[16]-.01464647509013483*gradb*t[11]*t[16];
    dfdab = 0.0;
    ds->df1000 += factor*dfdra;
    ds->df0100 += factor*dfdrb;
    ds->df0010 += factor*dfdga;
    ds->df0001 += factor*dfdgb;
    ds->df00001 += factor*dfdab;
}

static void
becke86_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real t[31];
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

    t[1] = pow(grada,2.0);
    t[2] = 1/pow(rhoa,2.666666666666667);
    t[3] = 0.004*t[1]*t[2]+1.0;
    t[4] = 1/t[3];
    t[5] = 1/pow(rhoa,2.333333333333333);
    t[6] = 1/pow(t[3],2.0);
    t[7] = 0.00787*t[1]*t[2]+1.0;
    t[8] = pow(gradb,2.0);
    t[9] = 1/pow(rhob,2.666666666666667);
    t[10] = 0.004*t[8]*t[9]+1.0;
    t[11] = 1/t[10];
    t[12] = 1/pow(rhob,2.333333333333333);
    t[13] = 1/pow(t[10],2.0);
    t[14] = 0.00787*t[8]*t[9]+1.0;
    t[15] = 1/pow(rhoa,1.333333333333333);
    t[16] = 1/pow(rhob,1.333333333333333);
    t[17] = pow(grada,4.0);
    t[18] = 1/pow(rhoa,6.0);
    t[19] = 1/pow(t[3],3.0);
    t[20] = 1/pow(rhoa,3.333333333333333);
    t[21] = pow(grada,3.0);
    t[22] = 1/pow(rhoa,5.0);
    t[23] = pow(gradb,4.0);
    t[24] = 1/pow(rhob,6.0);
    t[25] = 1/pow(t[10],3.0);
    t[26] = 1/pow(rhob,3.333333333333333);
    t[27] = pow(gradb,3.0);
    t[28] = 1/pow(rhob,5.0);
    t[29] = 1/pow(rhoa,4.0);
    t[30] = 1/pow(rhob,4.0);
    dfdra = -1.2407009817988*t[4]*t[7]*pow(rhoa,.3333333333333333)-0.0099256078543904*t[1]*t[6]*t[7]*t[5]+.01952863345351311*t[1]*t[4]*t[5];
    dfdrb = -1.2407009817988*t[11]*t[14]*pow(rhob,.3333333333333333)-0.0099256078543904*t[8]*t[13]*t[14]*t[12]+.01952863345351311*t[8]*t[11]*t[12];
    dfdga = 0.0074442058907928*grada*t[6]*t[7]*t[15]-.01464647509013483*grada*t[4]*t[15];
    dfdgb = 0.0074442058907928*gradb*t[13]*t[14]*t[16]-.01464647509013483*gradb*t[11]*t[16];
    dfdab = 0.0;
    d2fdrara = -.4135669939329333*t[4]*t[7]/pow(rhoa,.6666666666666666)+.009925607854390403*t[1]*t[6]*t[7]*t[20]-.01952863345351312*t[1]*t[4]*t[20]-2.1174630089366187E-4*t[17]*t[19]*t[7]*t[18]+4.166108470082797E-4*t[17]*t[6]*t[18];
    d2fdrarb = 0.0;
    d2fdraga = -0.0099256078543904*grada*t[6]*t[7]*t[5]+.01952863345351311*grada*t[4]*t[5]+1.588097256702464E-4*t[21]*t[19]*t[7]*t[22]-3.124581352562098E-4*t[21]*t[6]*t[22];
    d2fdragb = 0.0;
    d2fdraab = 0.0;
    d2fdrbrb = -.4135669939329333*t[11]*t[14]/pow(rhob,.6666666666666666)+.009925607854390403*t[8]*t[13]*t[14]*t[26]-.01952863345351312*t[8]*t[11]*t[26]-2.1174630089366187E-4*t[23]*t[25]*t[14]*t[24]+4.166108470082797E-4*t[23]*t[13]*t[24];
    d2fdrbga = 0.0;
    d2fdrbgb = -0.0099256078543904*gradb*t[13]*t[14]*t[12]+.01952863345351311*gradb*t[11]*t[12]+1.588097256702464E-4*t[27]*t[25]*t[14]*t[28]-3.124581352562098E-4*t[27]*t[13]*t[28];
    d2fdrbab = 0.0;
    d2fdgaga = 0.0074442058907928*t[6]*t[7]*t[15]-.01464647509013483*t[4]*t[15]-1.1910729425268479E-4*t[1]*t[19]*t[7]*t[29]+2.3434360144215735E-4*t[1]*t[6]*t[29];
    d2fdgagb = 0.0;
    d2fdgaab = 0.0;
    d2fdgbgb = 0.0074442058907928*t[13]*t[14]*t[16]-.01464647509013483*t[11]*t[16]-1.1910729425268479E-4*t[8]*t[25]*t[14]*t[30]+2.3434360144215735E-4*t[8]*t[13]*t[30];
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
becke86_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real t[49];
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

    t[1] = pow(grada,2.0);
    t[2] = 1/pow(rhoa,2.666666666666667);
    t[3] = 0.004*t[1]*t[2]+1.0;
    t[4] = 1/t[3];
    t[5] = 1/pow(rhoa,2.333333333333333);
    t[6] = 1/pow(t[3],2.0);
    t[7] = 0.00787*t[1]*t[2]+1.0;
    t[8] = pow(gradb,2.0);
    t[9] = 1/pow(rhob,2.666666666666667);
    t[10] = 0.004*t[8]*t[9]+1.0;
    t[11] = 1/t[10];
    t[12] = 1/pow(rhob,2.333333333333333);
    t[13] = 1/pow(t[10],2.0);
    t[14] = 0.00787*t[8]*t[9]+1.0;
    t[15] = 1/pow(rhoa,1.333333333333333);
    t[16] = 1/pow(rhob,1.333333333333333);
    t[17] = pow(grada,4.0);
    t[18] = 1/pow(rhoa,6.0);
    t[19] = 1/pow(t[3],3.0);
    t[20] = 1/pow(rhoa,3.333333333333333);
    t[21] = pow(grada,3.0);
    t[22] = 1/pow(rhoa,5.0);
    t[23] = pow(gradb,4.0);
    t[24] = 1/pow(rhob,6.0);
    t[25] = 1/pow(t[10],3.0);
    t[26] = 1/pow(rhob,3.333333333333333);
    t[27] = pow(gradb,3.0);
    t[28] = 1/pow(rhob,5.0);
    t[29] = 1/pow(rhoa,4.0);
    t[30] = 1/pow(rhob,4.0);
    t[31] = pow(grada,6.0);
    t[32] = 1/pow(rhoa,9.666666666666666);
    t[33] = 1/pow(t[3],4.0);
    t[34] = 1/pow(rhoa,7.0);
    t[35] = 1/pow(rhoa,4.333333333333333);
    t[36] = pow(grada,5.0);
    t[37] = 1/pow(rhoa,8.666666666666666);
    t[38] = 1/pow(rhoa,7.666666666666667);
    t[39] = pow(gradb,6.0);
    t[40] = 1/pow(rhob,9.666666666666666);
    t[41] = 1/pow(t[10],4.0);
    t[42] = 1/pow(rhob,7.0);
    t[43] = 1/pow(rhob,4.333333333333333);
    t[44] = pow(gradb,5.0);
    t[45] = 1/pow(rhob,8.666666666666666);
    t[46] = 1/pow(rhob,7.666666666666667);
    t[47] = 1/pow(rhoa,6.666666666666667);
    t[48] = 1/pow(rhob,6.666666666666667);
    dfdra = -1.2407009817988*t[4]*t[7]*pow(rhoa,.3333333333333333)-0.0099256078543904*t[1]*t[6]*t[7]*t[5]+.01952863345351311*t[1]*t[4]*t[5];
    dfdrb = -1.2407009817988*t[11]*t[14]*pow(rhob,.3333333333333333)-0.0099256078543904*t[8]*t[13]*t[14]*t[12]+.01952863345351311*t[8]*t[11]*t[12];
    dfdga = 0.0074442058907928*grada*t[6]*t[7]*t[15]-.01464647509013483*grada*t[4]*t[15];
    dfdgb = 0.0074442058907928*gradb*t[13]*t[14]*t[16]-.01464647509013483*gradb*t[11]*t[16];
    dfdab = 0.0;
    d2fdrara = -.4135669939329333*t[4]*t[7]/pow(rhoa,.6666666666666666)+.009925607854390403*t[1]*t[6]*t[7]*t[20]-.01952863345351312*t[1]*t[4]*t[20]-2.1174630089366187E-4*t[17]*t[19]*t[7]*t[18]+4.166108470082797E-4*t[17]*t[6]*t[18];
    d2fdrarb = 0.0;
    d2fdraga = -0.0099256078543904*grada*t[6]*t[7]*t[5]+.01952863345351311*grada*t[4]*t[5]+1.588097256702464E-4*t[21]*t[19]*t[7]*t[22]-3.124581352562098E-4*t[21]*t[6]*t[22];
    d2fdragb = 0.0;
    d2fdraab = 0.0;
    d2fdrbrb = -.4135669939329333*t[11]*t[14]/pow(rhob,.6666666666666666)+.009925607854390403*t[8]*t[13]*t[14]*t[26]-.01952863345351312*t[8]*t[11]*t[26]-2.1174630089366187E-4*t[23]*t[25]*t[14]*t[24]+4.166108470082797E-4*t[23]*t[13]*t[24];
    d2fdrbga = 0.0;
    d2fdrbgb = -0.0099256078543904*gradb*t[13]*t[14]*t[12]+.01952863345351311*gradb*t[11]*t[12]+1.588097256702464E-4*t[27]*t[25]*t[14]*t[28]-3.124581352562098E-4*t[27]*t[13]*t[28];
    d2fdrbab = 0.0;
    d2fdgaga = 0.0074442058907928*t[6]*t[7]*t[15]-.01464647509013483*t[4]*t[15]-1.1910729425268479E-4*t[1]*t[19]*t[7]*t[29]+2.3434360144215735E-4*t[1]*t[6]*t[29];
    d2fdgagb = 0.0;
    d2fdgaab = 0.0;
    d2fdgbgb = 0.0074442058907928*t[13]*t[14]*t[16]-.01464647509013483*t[11]*t[16]-1.1910729425268479E-4*t[8]*t[25]*t[14]*t[30]+2.3434360144215735E-4*t[8]*t[13]*t[30];
    d2fdgbab = 0.0;
    d2fdabab = 0.0;
    d3fdrarara= .2757113292886222*t[4]*t[7]/pow(rhoa,1.666666666666667)-.03749674078325263*t[1]*t[6]*t[7]*t[35]+.07377483749104954*t[1]*t[4]*t[35]+.001482224106255633*t[17]*t[19]*t[7]*t[34]-.002916275929057958*t[17]*t[6]*t[34]-6.77588162859718E-6*t[31]*t[33]*t[7]*t[32]+1.333154710426495E-5*t[31]*t[19]*t[32];
    d3fdrararb = 0.0;
    d3fdraraga = .02315975166024427*grada*t[6]*t[7]*t[20]-0.0455668113915306*grada*t[4]*t[20]-.001005794929244894*t[21]*t[19]*t[7]*t[18]+.001978901523289329*t[21]*t[6]*t[18]+5.081911221447886E-6*t[36]*t[33]*t[7]*t[37]-9.998660328198714E-6*t[36]*t[19]*t[37];
    d3fdraragb = 0.0;
    d3fdraraab = 0.0;
    d3fdrarbrb = 0.0;
    d3fdrarbga = 0.0;
    d3fdrarbgb = 0.0;
    d3fdrarbab = 0.0;
    d3fdragaga = -0.0099256078543904*t[6]*t[7]*t[5]+.01952863345351311*t[4]*t[5]+6.352389026809856E-4*t[1]*t[19]*t[7]*t[22]-.001249832541024839*t[1]*t[6]*t[22]-3.8114334160859137E-6*t[17]*t[33]*t[7]*t[38]+7.498995246149036E-6*t[17]*t[19]*t[38];
    d3fdragagb = 0.0;
    d3fdragaab = 0.0;
    d3fdragbgb = 0.0;
    d3fdragbab = 0.0;
    d3fdraabab = 0.0;
    d3fdrbrbrb = .2757113292886222*t[11]*t[14]/pow(rhob,1.666666666666667)-.03749674078325263*t[8]*t[13]*t[14]*t[43]+.07377483749104954*t[8]*t[11]*t[43]+.001482224106255633*t[23]*t[25]*t[14]*t[42]-.002916275929057958*t[23]*t[13]*t[42]+1.333154710426495E-5*t[39]*t[25]*t[40]-6.77588162859718E-6*t[39]*t[41]*t[14]*t[40];
    d3fdrbrbga = 0.0;
    d3fdrbrbgb = .02315975166024427*gradb*t[13]*t[14]*t[26]-0.0455668113915306*gradb*t[11]*t[26]-.001005794929244894*t[27]*t[25]*t[14]*t[24]+.001978901523289329*t[27]*t[13]*t[24]+5.081911221447886E-6*t[44]*t[41]*t[14]*t[45]-9.998660328198714E-6*t[44]*t[25]*t[45];
    d3fdrbrbab = 0.0;
    d3fdrbgaga= 0.0;
    d3fdrbgagb = 0.0;
    d3fdrbgaab = 0.0;
    d3fdrbgbgb = -0.0099256078543904*t[13]*t[14]*t[12]+.01952863345351311*t[11]*t[12]+6.352389026809856E-4*t[8]*t[25]*t[14]*t[28]-.001249832541024839*t[8]*t[13]*t[28]-3.8114334160859137E-6*t[23]*t[41]*t[14]*t[46]+7.498995246149036E-6*t[23]*t[25]*t[46];
    d3fdrbgbab = 0.0;
    d3fdrbabab = 0.0;
    d3fdgagaga = -3.573218827580544E-4*grada*t[19]*t[7]*t[29]+7.030308043264719E-4*grada*t[6]*t[29]+2.858575062064435E-6*t[21]*t[33]*t[7]*t[47]-5.624246434611776E-6*t[21]*t[19]*t[47];
    d3fdgagagb = 0.0;
    d3fdgagaab = 0.0;
    d3fdgagbgb = 0.0;
    d3fdgagbab = 0.0;
    d3fdgaabab = 0.0;
    d3fdgbgbgb = -3.573218827580544E-4*gradb*t[25]*t[14]*t[30]+7.030308043264719E-4*gradb*t[13]*t[30]+2.858575062064435E-6*t[27]*t[41]*t[14]*t[48]-5.624246434611776E-6*t[27]*t[25]*t[48];
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
becke86_fourth(FunFourthFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real t[71];
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

    t[1] = pow(grada,2.0);
    t[2] = 1/pow(rhoa,2.666666666666667);
    t[3] = 0.004*t[1]*t[2]+1.0;
    t[4] = 1/t[3];
    t[5] = 1/pow(rhoa,2.333333333333333);
    t[6] = 1/pow(t[3],2.0);
    t[7] = 0.00787*t[1]*t[2]+1.0;
    t[8] = pow(gradb,2.0);
    t[9] = 1/pow(rhob,2.666666666666667);
    t[10] = 0.004*t[8]*t[9]+1.0;
    t[11] = 1/t[10];
    t[12] = 1/pow(rhob,2.333333333333333);
    t[13] = 1/pow(t[10],2.0);
    t[14] = 0.00787*t[8]*t[9]+1.0;
    t[15] = 1/pow(rhoa,1.333333333333333);
    t[16] = 1/pow(rhob,1.333333333333333);
    t[17] = pow(grada,4.0);
    t[18] = 1/pow(rhoa,6.0);
    t[19] = 1/pow(t[3],3.0);
    t[20] = 1/pow(rhoa,3.333333333333333);
    t[21] = pow(grada,3.0);
    t[22] = 1/pow(rhoa,5.0);
    t[23] = pow(gradb,4.0);
    t[24] = 1/pow(rhob,6.0);
    t[25] = 1/pow(t[10],3.0);
    t[26] = 1/pow(rhob,3.333333333333333);
    t[27] = pow(gradb,3.0);
    t[28] = 1/pow(rhob,5.0);
    t[29] = 1/pow(rhoa,4.0);
    t[30] = 1/pow(rhob,4.0);
    t[31] = pow(grada,6.0);
    t[32] = 1/pow(rhoa,9.666666666666666);
    t[33] = 1/pow(t[3],4.0);
    t[34] = 1/pow(rhoa,7.0);
    t[35] = 1/pow(rhoa,4.333333333333333);
    t[36] = pow(grada,5.0);
    t[37] = 1/pow(rhoa,8.666666666666666);
    t[38] = 1/pow(rhoa,7.666666666666667);
    t[39] = pow(gradb,6.0);
    t[40] = 1/pow(rhob,9.666666666666666);
    t[41] = 1/pow(t[10],4.0);
    t[42] = 1/pow(rhob,7.0);
    t[43] = 1/pow(rhob,4.333333333333333);
    t[44] = pow(gradb,5.0);
    t[45] = 1/pow(rhob,8.666666666666666);
    t[46] = 1/pow(rhob,7.666666666666667);
    t[47] = 1/pow(rhoa,6.666666666666667);
    t[48] = 1/pow(rhob,6.666666666666667);
    t[49] = pow(grada,8.0);
    t[50] = 1/pow(rhoa,13.33333333333333);
    t[51] = 1/pow(t[3],5.0);
    t[52] = 1/pow(rhoa,10.66666666666667);
    t[53] = 1/pow(rhoa,8.0);
    t[54] = 1/pow(rhoa,5.333333333333333);
    t[55] = pow(grada,7.0);
    t[56] = 1/pow(rhoa,12.33333333333333);
    t[57] = 1/pow(rhoa,11.33333333333333);
    t[58] = 1/pow(rhoa,10.33333333333333);
    t[59] = pow(gradb,8.0);
    t[60] = 1/pow(rhob,13.33333333333333);
    t[61] = 1/pow(t[10],5.0);
    t[62] = 1/pow(rhob,10.66666666666667);
    t[63] = 1/pow(rhob,8.0);
    t[64] = 1/pow(rhob,5.333333333333333);
    t[65] = pow(gradb,7.0);
    t[66] = 1/pow(rhob,12.33333333333333);
    t[67] = 1/pow(rhob,11.33333333333333);
    t[68] = 1/pow(rhob,10.33333333333333);
    t[69] = 1/pow(rhoa,9.333333333333334);
    t[70] = 1/pow(rhob,9.333333333333334);
    dfdra = -1.2407009817988*t[4]*t[7]*pow(rhoa,.3333333333333333)-0.0099256078543904*t[1]*t[6]*t[7]*t[5]+.01952863345351311*t[1]*t[4]*t[5];
    dfdrb = -1.2407009817988*t[11]*t[14]*pow(rhob,.3333333333333333)-0.0099256078543904*t[8]*t[13]*t[14]*t[12]+.01952863345351311*t[8]*t[11]*t[12];
    dfdga = 0.0074442058907928*grada*t[6]*t[7]*t[15]-.01464647509013483*grada*t[4]*t[15];
    dfdgb = 0.0074442058907928*gradb*t[13]*t[14]*t[16]-.01464647509013483*gradb*t[11]*t[16];
    dfdab = 0.0;
    d2fdrara = -.4135669939329333*t[4]*t[7]/pow(rhoa,.6666666666666666)+.009925607854390403*t[1]*t[6]*t[7]*t[20]-.01952863345351312*t[1]*t[4]*t[20]-2.1174630089366187E-4*t[17]*t[19]*t[7]*t[18]+4.166108470082797E-4*t[17]*t[6]*t[18];
    d2fdrarb = 0.0;
    d2fdraga = -0.0099256078543904*grada*t[6]*t[7]*t[5]+.01952863345351311*grada*t[4]*t[5]+1.588097256702464E-4*t[21]*t[19]*t[7]*t[22]-3.124581352562098E-4*t[21]*t[6]*t[22];
    d2fdragb = 0.0;
    d2fdraab = 0.0;
    d2fdrbrb =-.4135669939329333*t[11]*t[14]/pow(rhob,.6666666666666666)+.009925607854390403*t[8]*t[13]*t[14]*t[26]-.01952863345351312*t[8]*t[11]*t[26]-2.1174630089366187E-4*t[23]*t[25]*t[14]*t[24]+4.166108470082797E-4*t[23]*t[13]*t[24];
    d2fdrbga = 0.0;
    d2fdrbgb = -0.0099256078543904*gradb*t[13]*t[14]*t[12]+.01952863345351311*gradb*t[11]*t[12]+1.588097256702464E-4*t[27]*t[25]*t[14]*t[28]-3.124581352562098E-4*t[27]*t[13]*t[28];
    d2fdrbab = 0.0;
    d2fdgaga = 0.0074442058907928*t[6]*t[7]*t[15]-.01464647509013483*t[4]*t[15]-1.1910729425268479E-4*t[1]*t[19]*t[7]*t[29]+2.3434360144215735E-4*t[1]*t[6]*t[29];
    d2fdgagb = 0.0;
    d2fdgaab = 0.0;
    d2fdgbgb= 0.0074442058907928*t[13]*t[14]*t[16]-.01464647509013483*t[11]*t[16]-1.1910729425268479E-4*t[8]*t[25]*t[14]*t[30]+2.3434360144215735E-4*t[8]*t[13]*t[30];
    d2fdgbab = 0.0;
    d2fdabab =0.0;
    d3fdrarara = .2757113292886222*t[4]*t[7]/pow(rhoa,1.666666666666667)-.03749674078325263*t[1]*t[6]*t[7]*t[35]+.07377483749104954*t[1]*t[4]*t[35]+.001482224106255633*t[17]*t[19]*t[7]*t[34]-.002916275929057958*t[17]*t[6]*t[34]-6.77588162859718E-6*t[31]*t[33]*t[7]*t[32]+1.333154710426495E-5*t[31]*t[19]*t[32];
    d3fdrararb = 0.0;
    d3fdraraga = .02315975166024427*grada*t[6]*t[7]*t[20]-0.0455668113915306*grada*t[4]*t[20]-.001005794929244894*t[21]*t[19]*t[7]*t[18]+.001978901523289329*t[21]*t[6]*t[18]+5.081911221447886E-6*t[36]*t[33]*t[7]*t[37]-9.998660328198714E-6*t[36]*t[19]*t[37];
    d3fdraragb = 0.0;
    d3fdraraab = 0.0;
    d3fdrarbrb = 0.0;
    d3fdrarbga = 0.0;
    d3fdrarbgb = 0.0;
    d3fdrarbab = 0.0;
    d3fdragaga = -0.0099256078543904*t[6]*t[7]*t[5]+.01952863345351311*t[4]*t[5]+6.352389026809856E-4*t[1]*t[19]*t[7]*t[22]-.001249832541024839*t[1]*t[6]*t[22]-3.8114334160859137E-6*t[17]*t[33]*t[7]*t[38]+7.498995246149036E-6*t[17]*t[19]*t[38];
    d3fdragagb = 0.0;
    d3fdragaab = 0.0;
    d3fdragbgb = 0.0;
    d3fdragbab = 0.0;
    d3fdraabab = 0.0;
    d3fdrbrbrb = .2757113292886222*t[11]*t[14]/pow(rhob,1.666666666666667)-.03749674078325263*t[8]*t[13]*t[14]*t[43]+.07377483749104954*t[8]*t[11]*t[43]+.001482224106255633*t[23]*t[25]*t[14]*t[42]-.002916275929057958*t[23]*t[13]*t[42]+1.333154710426495E-5*t[39]*t[25]*t[40]-6.77588162859718E-6*t[39]*t[41]*t[14]*t[40];
    d3fdrbrbga = 0.0;
    d3fdrbrbgb = .02315975166024427*gradb*t[13]*t[14]*t[26]-0.0455668113915306*gradb*t[11]*t[26]-.001005794929244894*t[27]*t[25]*t[14]*t[24]+.001978901523289329*t[27]*t[13]*t[24]+5.081911221447886E-6*t[44]*t[41]*t[14]*t[45]-9.998660328198714E-6*t[44]*t[25]*t[45];
    d3fdrbrbab =0.0;
    d3fdrbgaga = 0.0;
    d3fdrbgagb = 0.0;
    d3fdrbgaab = 0.0;
    d3fdrbgbgb = -0.0099256078543904*t[13]*t[14]*t[12]+.01952863345351311*t[11]*t[12]+6.352389026809856E-4*t[8]*t[25]*t[14]*t[28]-.001249832541024839*t[8]*t[13]*t[28]-3.8114334160859137E-6*t[23]*t[41]*t[14]*t[46]+7.498995246149036E-6*t[23]*t[25]*t[46];
    d3fdrbgbab = 0.0;
    d3fdrbabab = 0.0;
    d3fdgagaga =-3.573218827580544E-4*grada*t[19]*t[7]*t[29]+7.030308043264719E-4*grada*t[6]*t[29]+2.858575062064435E-6*t[21]*t[33]*t[7]*t[47]-5.624246434611776E-6*t[21]*t[19]*t[47];
    d3fdgagagb = 0.0;
    d3fdgagaab = 0.0;
    d3fdgagbgb = 0.0;
    d3fdgagbab = 0.0;
    d3fdgaabab = 0.0;
    d3fdgbgbgb = -3.573218827580544E-4*gradb*t[25]*t[14]*t[30]+7.030308043264719E-4*gradb*t[13]*t[30]+2.858575062064435E-6*t[27]*t[41]*t[14]*t[48]-5.624246434611776E-6*t[27]*t[25]*t[48];
    d3fdgbgbab = 0.0;
    d3fdgbabab = 0.0;
    d3fdababab = 0.0;
    d4fdrararara = -.4595188821477036*t[4]*t[7]*t[2]+.1654267975731734*t[1]*t[6]*t[7]*t[54]-.3254772242252185*t[1]*t[4]*t[54]-.01117549921383215*t[17]*t[19]*t[7]*t[53]+.02198779470321477*t[17]*t[6]*t[53]+1.1293136047661967E-4*t[31]*t[33]*t[7]*t[52]-2.2219245173774915E-4*t[31]*t[19]*t[52]-2.891042828201464E-7*t[49]*t[51]*t[7]*t[50]+5.68812676448638E-7*t[49]*t[33]*t[50];
    d4fdrarararb = 0.0;
    d4fdrararaga = -.07719917220081424*grada*t[6]*t[7]*t[35]+0.151889371305102*grada*t[4]*t[35]+.006528844277554575*t[21]*t[19]*t[7]*t[34]-.01284550111608863*t[21]*t[6]*t[34]-7.622866832171827E-5*t[36]*t[33]*t[7]*t[32]+1.499799049229807E-4*t[36]*t[19]*t[32]+2.1682821211510975E-7*t[55]*t[51]*t[7]*t[56]-4.2660950733647846E-7*t[55]*t[33]*t[56];
    d4fdrararagb = 0.0;
    d4fdrararaab = 0.0;
    d4fdrararbrb = 0.0;
    d4fdrararbga = 0.0;
    d4fdrararbgb = 0.0;
    d4fdrararbab = 0.0;
    d4fdraragaga = .02315975166024427*t[6]*t[7]*t[20]-0.0455668113915306*t[4]*t[20]-0.00338794081429859*t[1]*t[19]*t[7]*t[18]+.006665773552132477*t[1]*t[6]*t[18]+4.954863440911688E-5*t[17]*t[33]*t[7]*t[37]-9.748693819993747E-5*t[17]*t[19]*t[37]-1.6262115908633232E-7*t[31]*t[51]*t[7]*t[57]+3.1995713050235886E-7*t[31]*t[33]*t[57];
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
    d4fdrarbgagb= 0.0;
    d4fdrarbgaab = 0.0;
    d4fdrarbgbgb = 0.0;
    d4fdrarbgbab = 0.0;
    d4fdrarbabab =0.0;
    d4fdragagaga = .001429287531032218*grada*t[19]*t[7]*t[22]-.002812123217305888*grada*t[6]*t[22]-3.049146732868731E-5*t[21]*t[33]*t[7]*t[38]+5.99919619691923E-5*t[21]*t[19]*t[38]+1.2196586931474923E-7*t[36]*t[51]*t[7]*t[58]-2.399678478767691E-7*t[36]*t[33]*t[58];
    d4fdragagagb = 0.0;
    d4fdragagaab = 0.0;
    d4fdragagbgb = 0.0;
    d4fdragagbab = 0.0;
    d4fdragaabab = 0.0;
    d4fdragbgbgb = 0.0;
    d4fdragbgbab = 0.0;
    d4fdragbabab = 0.0;
    d4fdraababab = 0.0;
    d4fdrbrbrbrb = -.4595188821477036*t[11]*t[14]*t[9]+.1654267975731734*t[8]*t[13]*t[14]*t[64]-.3254772242252185*t[8]*t[11]*t[64]-.01117549921383215*t[23]*t[25]*t[14]*t[63]+.02198779470321477*t[23]*t[13]*t[63]+1.1293136047661967E-4*t[39]*t[41]*t[14]*t[62]-2.2219245173774915E-4*t[39]*t[25]*t[62]-2.891042828201464E-7*t[59]*t[61]*t[14]*t[60]+5.68812676448638E-7*t[59]*t[41]*t[60];
    d4fdrbrbrbga = 0.0;
    d4fdrbrbrbgb = -.07719917220081424*gradb*t[13]*t[14]*t[43]+0.151889371305102*gradb*t[11]*t[43]+.006528844277554575*t[27]*t[25]*t[14]*t[42]-.01284550111608863*t[27]*t[13]*t[42]-7.622866832171827E-5*t[44]*t[41]*t[14]*t[40]+1.499799049229807E-4*t[44]*t[25]*t[40]+2.1682821211510975E-7*t[65]*t[61]*t[14]*t[66]-4.2660950733647846E-7*t[65]*t[41]*t[66];
    d4fdrbrbrbab = 0.0;
    d4fdrbrbgaga = 0.0;
    d4fdrbrbgagb = 0.0;
    d4fdrbrbgaab = 0.0;
    d4fdrbrbgbgb = .02315975166024427*t[13]*t[14]*t[26]-0.0455668113915306*t[11]*t[26]-0.00338794081429859*t[8]*t[25]*t[14]*t[24]+.006665773552132477*t[8]*t[13]*t[24]+4.954863440911688E-5*t[23]*t[41]*t[14]*t[45]-9.748693819993747E-5*t[23]*t[25]*t[45]-1.6262115908633232E-7*t[39]*t[61]*t[14]*t[67]+3.1995713050235886E-7*t[39]*t[41]*t[67];
    d4fdrbrbgbab = 0.0;
    d4fdrbrbabab = 0.0;
    d4fdrbgagaga = 0.0;
    d4fdrbgagagb = 0.0;
    d4fdrbgagaab = 0.0;
    d4fdrbgagbgb = 0.0;
    d4fdrbgagbab = 0.0;
    d4fdrbgaabab = 0.0;
    d4fdrbgbgbgb = .001429287531032218*gradb*t[25]*t[14]*t[28]-.002812123217305888*gradb*t[13]*t[28]-3.049146732868731E-5*t[27]*t[41]*t[14]*t[46]+5.99919619691923E-5*t[27]*t[25]*t[46]+1.2196586931474923E-7*t[44]*t[61]*t[14]*t[68]-2.399678478767691E-7*t[44]*t[41]*t[68];
    d4fdrbgbgbab = 0.0;
    d4fdrbgbabab = 0.0;
    d4fdrbababab = 0.0;
    d4fdgagagaga = -3.573218827580544E-4*t[19]*t[7]*t[29]+7.030308043264719E-4*t[6]*t[29]+1.715145037238661E-5*t[1]*t[33]*t[7]*t[47]-3.3745478607670654E-5*t[1]*t[19]*t[47]-9.147440198606192E-8*t[17]*t[51]*t[7]*t[69]+1.7997588590757682E-7*t[17]*t[33]*t[69];
    d4fdgagagagb = 0.0;
    d4fdgagagaab = 0.0;
    d4fdgagagbgb = 0.0;
    d4fdgagagbab = 0.0;
    d4fdgagaabab = 0.0;
    d4fdgagbgbgb = 0.0;
    d4fdgagbgbab = 0.0;
    d4fdgagbabab = 0.0;
    d4fdgaababab = 0.0;
    d4fdgbgbgbgb = -3.573218827580544E-4*t[25]*t[14]*t[30]+7.030308043264719E-4*t[13]*t[30]+1.715145037238661E-5*t[8]*t[41]*t[14]*t[48]-3.3745478607670654E-5*t[8]*t[25]*t[48]-9.147440198606192E-8*t[23]*t[61]*t[14]*t[70]+1.7997588590757682E-7*t[23]*t[41]*t[70];
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
