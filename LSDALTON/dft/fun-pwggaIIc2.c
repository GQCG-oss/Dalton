/*
C...   Copyright (c) 2015 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren. 
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
*/
/*
C...   Copyright (c) 2015 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras, T. Saue, 
C...   S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren. 
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
*/
/* Automatically generated functional code: pwggaIIc2
   Maxima input:
    >> PI:    3.14159265358979312;
    >> Cc0:0.004235;
    >> 
    >> rho:  rhoa + rhob;
    >> grad: sqrt(grada*grada + gradb*gradb + 2*gradab);
    >> zeta: (rhoa-rhob)/(rhoa+rhob);
    >> 
    >> gzgr: (1-zeta)/rho*grada^2 - (1+zeta)/rho*gradb^2 - (2*zeta/rho)*gradab;
    >> t1: -0.458*zeta * gzgr/((rho*(1-zeta^2))^(1/3)*rho);
    >> gz2: (1-zeta)^2/rho^2*grada^2 + (1+zeta)^2/rho^2*gradb^2 - 2*(1-zeta^2)/rho^2*gradab;
    >> t2: ( (-0.037+0.10*zeta^2)*gz2 / rho^(1/3)*(1-zeta^2) );
    >> 
    >> 
    >> 
    >> K(rhoa,rhob,grada,gradb,gradab):=Cc0*rho*(t1+t2);
*/

#include <math.h>
#include <stddef.h>

#define __CVERSION__

#include "lsdalton_functionals.h"
#define LOG log
#define ABS fabs
#define ASINH asinh
#define SQRT sqrt

/* INTERFACE PART */
static integer pwggaIIc2_read(const char* conf_line, real *hfweight);
static real pwggaIIc2_energy(const FunDensProp* dp);
static void pwggaIIc2_first(FirstFuncDrv *ds, real factor, 
                       const FunDensProp* dp);
static void pwggaIIc2_second(SecondFuncDrv *ds, real factor,
                        const FunDensProp* dp);
static void pwggaIIc2_third(ThirdFuncDrv *ds, real factor,
                       const FunDensProp* dp);

Functional PWggaIIc2Functional = {
  "pwggaIIc2",
  fun_true,
  pwggaIIc2_read,
  NULL,
  pwggaIIc2_energy,
  pwggaIIc2_first,
  pwggaIIc2_second,
  pwggaIIc2_third
};

/* IMPLEMENTATION PART */
static integer
pwggaIIc2_read(const char* conf_line, real *hfweight)
{
    return 1;
}


static real
pwggaIIc2_energy(const FunDensProp* dp)
{
    real t[11],zk;
    real rhoa = dp->rhoa;
    real rhob = dp->rhob;
    real grada = dp->grada;
    real gradb = dp->gradb;
    real gradab = dp->gradab;

    t[1] = rhob+rhoa;
    t[2] = rhoa-1.0*rhob;
    t[3] = pow(t[2],2.0);
    t[4] = 1/pow(t[1],2.0);
    t[5] = 1.0-1.0*t[3]*t[4];
    t[6] = pow(grada,2.0);
    t[7] = 1/t[1];
    t[8] = 1.0-1.0*t[2]*t[7];
    t[9] = pow(gradb,2.0);
    t[10] = t[2]*t[7]+1.0;
    zk = 0.004235*t[1]*((0.1*t[3]*t[4]-0.037)*t[5]*(-2.0*t[4]*t[5]*gradab+pow(t[10],2.0)*t[4]*t[9]+t[4]*t[6]*pow(t[8],2.0))/pow(t[1],0.33333333333333)-0.458*t[2]*(-2.0*t[2]*t[4]*gradab-1.0*t[10]*t[7]*t[9]+t[6]*t[7]*t[8])/(pow(t[1],2.333333333333334)*pow(t[5],0.33333333333333)));
    return zk;
}

static void
pwggaIIc2_first(FirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real t[42];
    real dfdra, dfdrb, dfdga, dfdgb, dfdab;
    real rhoa = dp->rhoa;
    real rhob = dp->rhob;
    real grada = dp->grada;
    real gradb = dp->gradb;
    real gradab = dp->gradab;

    t[1] = rhob+rhoa;
    t[2] = rhoa-1.0*rhob;
    t[3] = 1/pow(t[1],2.333333333333334);
    t[4] = pow(t[2],2.0);
    t[5] = 1/pow(t[1],2.0);
    t[6] = 1.0-1.0*t[4]*t[5];
    t[7] = 1/pow(t[6],0.33333333333333);
    t[8] = 1/pow(t[1],3.0);
    t[9] = 4.0*t[2]*t[8]*gradab;
    t[10] = pow(grada,2.0);
    t[11] = 1/t[1];
    t[12] = t[2]*t[5];
    t[13] = -1.0*t[11];
    t[14] = t[13]+t[12];
    t[15] = pow(gradb,2.0);
    t[16] = -1.0*t[2]*t[5];
    t[17] = t[11]+t[16];
    t[18] = 1.0-1.0*t[11]*t[2];
    t[19] = -1.0*t[10]*t[18]*t[5];
    t[20] = t[2]*t[11]+1.0;
    t[21] = t[15]*t[5]*t[20];
    t[22] = 2.0*t[4]*t[8];
    t[23] = t[22]-2.0*t[2]*t[5];
    t[24] = 1/pow(t[6],1.333333333333333);
    t[25] = -2.0*t[2]*t[5]*gradab-1.0*t[11]*t[15]*t[20]+t[10]*t[11]*t[18];
    t[26] = 1/pow(t[1],3.333333333333334);
    t[27] = 1.068666666666667*t[2]*t[26]*t[7]*t[25];
    t[28] = 1/pow(t[1],0.33333333333333);
    t[29] = 0.1*t[4]*t[5]-0.037;
    t[30] = 4.0*t[6]*t[8]*gradab;
    t[31] = pow(t[18],2.0);
    t[32] = -2.0*t[10]*t[31]*t[8];
    t[33] = pow(t[20],2.0);
    t[34] = -2.0*t[15]*t[33]*t[8];
    t[35] = -0.2*t[4]*t[8];
    t[36] = -2.0*t[5]*t[6]*gradab+t[15]*t[5]*t[33]+t[10]*t[5]*t[31];
    t[37] = -0.33333333333333*t[29]*t[36]*t[6]/pow(t[1],1.333333333333333);
    t[38] = 0.004235*(t[28]*t[6]*t[29]*t[36]-0.458*t[2]*t[3]*t[7]*t[25]);
    t[39] = t[13]+t[16];
    t[40] = t[11]+t[12];
    t[41] = 2.0*t[2]*t[5]+t[22];
    dfdra = 0.004235*t[1]*(t[28]*t[29]*t[6]*(-2.0*t[23]*t[5]*gradab+2.0*t[15]*t[17]*t[20]*t[5]+2.0*t[10]*t[14]*t[18]*t[5]+t[34]+t[32]+t[30])-0.458*t[2]*t[3]*t[7]*(-2.0*t[5]*gradab+t[9]+t[21]+t[19]-1.0*t[11]*t[15]*t[17]+t[10]*t[11]*t[14])+t[37]+t[28]*(0.2*t[2]*t[5]+t[35])*t[6]*t[36]+t[28]*t[23]*t[29]*t[36]+t[27]-0.458*t[3]*t[7]*t[25]+0.15266666666667*t[2]*t[3]*t[23]*t[24]*t[25])+t[38];
    dfdrb = 0.004235*t[1]*(t[28]*t[29]*t[6]*(-2.0*t[41]*t[5]*gradab+2.0*t[10]*t[18]*t[40]*t[5]+2.0*t[15]*t[20]*t[39]*t[5]+t[34]+t[32]+t[30])-0.458*t[2]*t[3]*t[7]*(2.0*t[5]*gradab+t[9]+t[10]*t[11]*t[40]-1.0*t[11]*t[15]*t[39]+t[21]+t[19])+t[37]+t[28]*(t[35]-0.2*t[2]*t[5])*t[6]*t[36]+t[28]*t[41]*t[29]*t[36]+t[27]+0.458*t[3]*t[7]*t[25]+0.15266666666667*t[2]*t[3]*t[41]*t[24]*t[25])+t[38];
    dfdga = 0.004235*t[1]*(2.0*t[29]*t[3]*t[31]*t[6]*grada-0.916*grada*t[2]*t[26]*t[7]*t[18]);
    dfdgb = 0.004235*t[1]*(2.0*t[29]*t[3]*t[33]*t[6]*gradb+0.916*gradb*t[2]*t[26]*t[7]*t[20]);
    dfdab = 0.004235*t[1]*(0.916*t[4]*t[7]/pow(t[1],4.333333333333333)-2.0*t[29]*t[3]*pow(t[6],2.0));
    ds->df1000 += factor*dfdra;
    ds->df0100 += factor*dfdrb;
    ds->df0010 += factor*dfdga;
    ds->df0001 += factor*dfdgb;
    ds->df00001 += factor*dfdab;
}

static void
pwggaIIc2_second(SecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real t[96];
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

    t[1] = rhob+rhoa;
    t[2] = rhoa-1.0*rhob;
    t[3] = 1/pow(t[1],2.333333333333334);
    t[4] = pow(t[2],2.0);
    t[5] = 1/pow(t[1],2.0);
    t[6] = 1.0-1.0*t[4]*t[5];
    t[7] = 1/pow(t[6],0.33333333333333);
    t[8] = 1/pow(t[1],3.0);
    t[9] = 4.0*t[2]*t[8]*gradab;
    t[10] = pow(grada,2.0);
    t[11] = 1/t[1];
    t[12] = t[2]*t[5];
    t[13] = -1.0*t[11];
    t[14] = t[13]+t[12];
    t[15] = pow(gradb,2.0);
    t[16] = -1.0*t[2]*t[5];
    t[17] = t[11]+t[16];
    t[18] = 1.0-1.0*t[11]*t[2];
    t[19] = -1.0*t[10]*t[18]*t[5];
    t[20] = t[2]*t[11]+1.0;
    t[21] = t[15]*t[5]*t[20];
    t[22] = -2.0*t[5]*gradab+t[9]+t[21]+t[19]-1.0*t[11]*t[15]*t[17]+t[10]*t[11]*t[14];
    t[23] = 2.0*t[4]*t[8];
    t[24] = t[23]-2.0*t[2]*t[5];
    t[25] = 1/pow(t[6],1.333333333333333);
    t[26] = -2.0*t[2]*t[5]*gradab-1.0*t[11]*t[15]*t[20]+t[10]*t[11]*t[18];
    t[27] = 1/pow(t[1],3.333333333333334);
    t[28] = 1.068666666666667*t[2]*t[27]*t[7]*t[26];
    t[29] = 1/pow(t[1],0.33333333333333);
    t[30] = 0.1*t[4]*t[5]-0.037;
    t[31] = 4.0*t[6]*t[8]*gradab;
    t[32] = pow(t[18],2.0);
    t[33] = -2.0*t[10]*t[32]*t[8];
    t[34] = pow(t[20],2.0);
    t[35] = -2.0*t[15]*t[34]*t[8];
    t[36] = -2.0*t[24]*t[5]*gradab+2.0*t[15]*t[17]*t[20]*t[5]+2.0*t[10]*t[14]*t[18]*t[5]+t[35]+t[33]+t[31];
    t[37] = -0.2*t[4]*t[8];
    t[38] = 0.2*t[2]*t[5]+t[37];
    t[39] = -2.0*t[5]*t[6]*gradab+t[15]*t[5]*t[34]+t[10]*t[5]*t[32];
    t[40] = 1/pow(t[1],1.333333333333333);
    t[41] = -0.33333333333333*t[30]*t[39]*t[40]*t[6];
    t[42] = t[41]+t[29]*t[24]*t[30]*t[39]+t[29]*t[38]*t[6]*t[39]+t[29]*t[6]*t[30]*t[36]-0.458*t[3]*t[7]*t[26]+t[28]+0.15266666666667*t[2]*t[3]*t[24]*t[25]*t[26]-0.458*t[2]*t[3]*t[7]*t[22];
    t[43] = 0.004235*(t[29]*t[6]*t[30]*t[39]-0.458*t[2]*t[3]*t[7]*t[26]);
    t[44] = t[13]+t[16];
    t[45] = t[11]+t[12];
    t[46] = 2.0*t[5]*gradab+t[9]+t[10]*t[11]*t[45]-1.0*t[11]*t[15]*t[44]+t[21]+t[19];
    t[47] = 2.0*t[2]*t[5]+t[23];
    t[48] = -2.0*t[47]*t[5]*gradab+2.0*t[10]*t[18]*t[45]*t[5]+2.0*t[15]*t[20]*t[44]*t[5]+t[35]+t[33]+t[31];
    t[49] = t[37]-0.2*t[2]*t[5];
    t[50] = t[41]+t[29]*t[47]*t[30]*t[39]+t[29]*t[49]*t[6]*t[39]+t[29]*t[6]*t[30]*t[48]+0.458*t[3]*t[7]*t[26]+t[28]+0.15266666666667*t[2]*t[3]*t[47]*t[25]*t[26]-0.458*t[2]*t[3]*t[7]*t[46];
    t[51] = 2.0*t[3]*t[30]*t[32]*t[6]*grada-0.916*grada*t[2]*t[27]*t[7]*t[18];
    t[52] = 2.0*t[3]*t[30]*t[34]*t[6]*gradb+0.916*gradb*t[2]*t[27]*t[7]*t[20];
    t[53] = 1/pow(t[1],4.333333333333333);
    t[54] = pow(t[6],2.0);
    t[55] = 0.916*t[4]*t[53]*t[7]-2.0*t[3]*t[30]*t[54];
    t[56] = 1/pow(t[1],4.0);
    t[57] = -12.0*t[2]*t[56]*gradab;
    t[58] = 2.0*t[2]*t[8];
    t[59] = -2.0*t[5];
    t[60] = t[59]+t[58];
    t[61] = -2.0*t[2]*t[8];
    t[62] = 2.0*t[5];
    t[63] = t[62]+t[61];
    t[64] = 2.0*t[10]*t[18]*t[8];
    t[65] = -2.0*t[15]*t[20]*t[8];
    t[66] = 1/pow(t[6],2.333333333333334);
    t[67] = -6.0*t[4]*t[56];
    t[68] = 8.0*t[2]*t[8]+t[67]+t[59];
    t[69] = -3.562222222222223*t[2]*t[53]*t[7]*t[26];
    t[70] = -12.0*t[56]*t[6]*gradab;
    t[71] = 6.0*t[10]*t[32]*t[56];
    t[72] = 6.0*t[15]*t[34]*t[56];
    t[73] = 0.6*t[4]*t[56];
    t[74] = 0.2*t[5];
    t[75] = 0.44444444444444*t[3]*t[30]*t[39]*t[6];
    t[76] = t[62]+t[67];
    t[77] = 1/pow(t[1],5.0);
    t[78] = 0.004235*t[51];
    t[79] = 2.137333333333334*grada*t[2]*t[53]*t[7]*t[18];
    t[80] = -0.66666666666667*t[27]*t[30]*t[32]*t[6]*grada;
    t[81] = -2.0*t[18]*t[5]*grada;
    t[82] = -4.0*t[32]*t[8]*grada;
    t[83] = 0.004235*t[52];
    t[84] = -2.137333333333334*gradb*t[2]*t[53]*t[7]*t[20];
    t[85] = -0.66666666666667*t[27]*t[30]*t[34]*t[6]*gradb;
    t[86] = 2.0*t[20]*t[5]*gradb;
    t[87] = -4.0*t[34]*t[8]*gradb;
    t[88] = 0.004235*t[55];
    t[89] = -2.137333333333334*t[4]*t[7]/pow(t[1],5.333333333333333);
    t[90] = 4.0*t[2]*t[8];
    t[91] = 0.66666666666667*t[27]*t[30]*t[54];
    t[92] = 4.0*t[6]*t[8];
    t[93] = t[59]+t[61];
    t[94] = t[62]+t[58];
    t[95] = -8.0*t[2]*t[8]+t[67]+t[59];
    dfdra = t[43]+0.004235*t[1]*t[42];
    dfdrb = t[43]+0.004235*t[1]*t[50];
    dfdga = 0.004235*t[1]*t[51];
    dfdgb = 0.004235*t[1]*t[52];
    dfdab = 0.004235*t[1]*t[55];
    d2fdrara = 0.004235*t[1]*(t[29]*t[30]*t[6]
			      *(8.0*t[24]*t[8]*gradab-2.0*t[5]*t[68]*gradab-8.0*t[15]*t[17]*t[20]*t[8]-8.0*t[10]*t[14]*t[18]*t[8]+t[72]+t[71]+t[70]+2.0*t[10]*t[18]*t[5]*t[63]+2.0*t[15]*t[20]*t[5]*t[60]+2.0*t[15]*pow(t[17],2.0)*t[5]+2.0*t[10]*pow(t[14],2.0)*t[5])
			      -0.458*t[2]*t[3]*t[7]*(8.0*t[8]*gradab+t[65]+t[64]+t[10]*t[11]*t[63]-1.0*t[11]*t[15]*t[60]+t[57]+2.0*t[15]*t[17]*t[5]-2.0*t[10]*t[14]*t[5])+t[75]+t[69]-0.20355555555556*t[2]*pow(t[24],2.0)*t[26]*t[3]*t[66]
			      -0.66666666666667*t[38]*t[39]*t[40]*t[6]-0.66666666666667*t[30]*t[36]*t[40]*t[6]+2.0*t[29]*t[36]*t[38]*t[6]-0.66666666666667*t[24]*t[30]*t[39]*t[40]+t[29]*(t[74]-0.8*t[2]*t[8]+t[73])*t[6]*t[39]
			      +2.0*t[24]*t[29]*t[38]*t[39]+t[29]*t[68]*t[30]*t[39]+2.0*t[24]*t[29]*t[30]*t[36]+2.137333333333334*t[27]*t[7]*t[26]+0.15266666666667*t[2]*t[3]*t[68]*t[25]*t[26]+0.30533333333333*t[3]*t[24]*t[25]*t[26]
			      -0.71244444444444*t[2]*t[27]*t[24]*t[25]*t[26]-0.916*t[3]*t[7]*t[22]+2.137333333333334*t[2]*t[27]*t[7]*t[22]+0.30533333333333*t[2]*t[3]*t[24]*t[25]*t[22])+0.00847*t[42];
    d2fdrarb = 0.004235*t[1]
      *(t[29]*t[30]*t[6]*(4.0*t[47]*t[8]*gradab+4.0*t[24]*t[8]*gradab-2.0*t[5]*t[76]*gradab-4.0*t[10]*t[18]*t[45]*t[8]-4.0*t[15]*t[20]*t[44]*t[8]-4.0*t[15]*t[17]*t[20]*t[8]-4.0*t[10]*t[14]*t[18]*t[8]+4.0*t[15]*t[2]*t[20]*t[77]
			  -4.0*t[10]*t[18]*t[2]*t[77]+t[72]+t[71]+t[70]+2.0*t[10]*t[14]*t[45]*t[5]+2.0*t[15]*t[17]*t[44]*t[5])+t[75]-0.458*t[2]*t[3]
	*(t[65]+t[64]+t[57]-2.0*t[15]*t[2]*t[56]-2.0*t[10]*t[2]*t[56]-1.0*t[10]*t[45]*t[5]-1.0*t[10]*t[14]*t[5]+t[15]*t[5]*t[44]+t[15]*t[5]*t[17])*t[7]+t[69]-0.33333333333333*t[39]*t[40]*t[49]*t[6]
	-0.33333333333333*t[30]*t[40]*t[48]*t[6]-0.33333333333333*t[38]*t[39]*t[40]*t[6]-0.33333333333333*t[30]*t[36]*t[40]*t[6]+t[29]*t[38]*t[6]*t[48]+t[29]*t[24]*t[30]*t[48]-0.33333333333333*t[30]*t[39]*t[40]*t[47]
	-0.458*t[3]*t[7]*t[46]+1.068666666666667*t[2]*t[27]*t[7]*t[46]+0.15266666666667*t[2]*t[3]*t[24]*t[25]*t[46]-0.33333333333333*t[24]*t[30]*t[39]*t[40]+t[29]*(t[73]-0.2*t[5])*t[6]*t[39]+t[29]*t[24]*t[49]*t[39]
	+t[29]*t[38]*t[47]*t[39]+t[29]*t[76]*t[30]*t[39]+t[29]*t[49]*t[6]*t[36]+t[29]*t[47]*t[30]*t[36]-0.20355555555556*t[2]*t[3]*t[24]*t[47]*t[66]*t[26]+0.15266666666667*t[2]*t[3]*t[76]*t[25]*t[26]
	+0.15266666666667*t[3]*t[47]*t[25]*t[26]-0.35622222222222*t[2]*t[27]*t[47]*t[25]*t[26]-0.15266666666667*t[3]*t[24]*t[25]*t[26]-0.35622222222222*t[2]*t[27]*t[24]*t[25]*t[26]+0.458*t[3]*t[7]*t[22]
	+1.068666666666667*t[2]*t[27]*t[7]*t[22]+0.15266666666667*t[2]*t[3]*t[47]*t[25]*t[22])+0.004235*t[50]+0.004235*t[42];
    d2fdraga = 0.004235*t[1]*(t[29]*t[30]*t[6]*(4.0*t[14]*t[18]*t[5]*grada+t[82])-0.458*t[2]*t[3]*t[7]*(2.0*t[11]*t[14]*grada+t[81])+2.0*t[3]*t[32]*t[38]*t[6]*grada+2.0*t[24]*t[3]*t[30]*t[32]*grada+t[80]+t[79]-0.916*grada*t[27]*t[7]*t[18]+0.30533333333333*grada*t[2]*t[27]*t[24]*t[25]*t[18])+t[78];
    d2fdragb = 0.004235*t[1]*(t[29]*t[30]*t[6]*(4.0*t[17]*t[20]*t[5]*gradb+t[87])-0.458*t[2]*t[3]*t[7]*(t[86]-2.0*t[11]*t[17]*gradb)+2.0*t[3]*t[34]*t[38]*t[6]*gradb+2.0*t[24]*t[3]*t[30]*t[34]*gradb+t[85]+t[84]+0.916*gradb*t[27]*t[7]*t[20]-0.30533333333333*gradb*t[2]*t[27]*t[24]*t[25]*t[20])+t[83];
    d2fdraab = 0.004235*t[1]*(t[29]*t[30]*t[6]*(t[92]-2.0*t[24]*t[5])+t[91]+t[89]-0.458*t[2]*t[3]*(t[59]+t[90])*t[7]+0.916*t[2]*t[53]*t[7]-2.0*t[24]*t[3]*t[30]*t[6]-2.0*t[3]*t[38]*t[54]-0.30533333333333*t[4]*t[53]*t[24]*t[25])+t[88];
    d2fdrbrb = 0.004235*t[1]
      *(t[29]*t[30]*t[6]
	*(-2.0*t[5]*t[95]*gradab+8.0*t[47]*t[8]*gradab+2.0*t[15]*t[20]*t[5]*t[94]+2.0*t[10]*t[18]*t[5]*t[93]-8.0*t[10]*t[18]*t[45]*t[8]-8.0*t[15]*t[20]*t[44]*t[8]+t[72]+t[71]+t[70]+2.0*t[10]*pow(t[45],2.0)*t[5]
	  +2.0*t[15]*pow(t[44],2.0)*t[5])-0.458*t[2]*t[3]*t[7]*(-8.0*t[8]*gradab-1.0*t[11]*t[15]*t[94]+t[10]*t[11]*t[93]+t[65]+t[64]+t[57]-2.0*t[10]*t[45]*t[5]+2.0*t[15]*t[44]*t[5])+t[75]+t[69]
	-0.20355555555556*t[2]*t[26]*t[3]*pow(t[47],2.0)*t[66]+2.0*t[29]*t[48]*t[49]*t[6]-0.66666666666667*t[39]*t[40]*t[49]*t[6]-0.66666666666667*t[30]*t[40]*t[48]*t[6]+2.0*t[29]*t[39]*t[47]*t[49]
	+2.0*t[29]*t[30]*t[47]*t[48]-0.66666666666667*t[30]*t[39]*t[40]*t[47]+0.916*t[3]*t[7]*t[46]+2.137333333333334*t[2]*t[27]*t[7]*t[46]+0.30533333333333*t[2]*t[3]*t[47]*t[25]*t[46]+t[29]*(t[74]+0.8*t[2]*t[8]+t[73])*t[6]*t[39]
	+t[29]*t[95]*t[30]*t[39]-2.137333333333334*t[27]*t[7]*t[26]+0.15266666666667*t[2]*t[3]*t[95]*t[25]*t[26]-0.30533333333333*t[3]*t[47]*t[25]*t[26]-0.71244444444444*t[2]*t[27]*t[47]*t[25]*t[26])+0.00847*t[50];
    d2fdrbga = 0.004235*t[1]*(t[29]*t[30]*t[6]*(4.0*t[18]*t[45]*t[5]*grada+t[82])-0.458*t[2]*t[3]*t[7]*(2.0*t[11]*t[45]*grada+t[81])+2.0*t[3]*t[32]*t[49]*t[6]*grada+2.0*t[3]*t[30]*t[32]*t[47]*grada+t[80]+t[79]+0.916*grada*t[27]*t[7]*t[18]+0.30533333333333*grada*t[2]*t[27]*t[47]*t[25]*t[18])+t[78];
    d2fdrbgb = 0.004235*t[1]*(t[29]*t[30]*t[6]*(4.0*t[20]*t[44]*t[5]*gradb+t[87])-0.458*t[2]*t[3]*t[7]*(t[86]-2.0*t[11]*t[44]*gradb)+2.0*t[3]*t[34]*t[49]*t[6]*gradb+2.0*t[3]*t[30]*t[34]*t[47]*gradb+t[85]+t[84]-0.916*gradb*t[27]*t[7]*t[20]-0.30533333333333*gradb*t[2]*t[27]*t[47]*t[25]*t[20])+t[83];
    d2fdrbab = 0.004235*t[1]*(t[29]*t[30]*t[6]*(t[92]-2.0*t[47]*t[5])+t[91]+t[89]-0.458*t[2]*t[3]*(t[62]+t[90])*t[7]-0.916*t[2]*t[53]*t[7]-2.0*t[3]*t[30]*t[47]*t[6]-2.0*t[3]*t[49]*t[54]-0.30533333333333*t[4]*t[53]*t[47]*t[25])+t[88];
    d2fdgaga = 0.004235*t[1]*(2.0*t[3]*t[30]*t[32]*t[6]-0.916*t[2]*t[27]*t[7]*t[18]);
    d2fdgagb = 0.0;
    d2fdgaab = 0.0;
    d2fdgbgb = 0.004235*t[1]*(2.0*t[3]*t[30]*t[34]*t[6]+0.916*t[2]*t[27]*t[7]*t[20]);
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
static void pwggaIIc2_third(ThirdFuncDrv *ds, real factor,
                       const FunDensProp* dp){}
