/* -*-mode:c; c-style:bsd; c-basic-offset:4;indent-tabs-mode:nil; -*- */
/* fun-vwn.c:
   implementation of VWN functional and its derivatives 
   (c), Pawel Salek, pawsa@theochem.kth.se, sep 2001, nov 2002
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
static int  vwn_isgga(void) { return 0; }
static int  vwn_read(const char* conf_line);
static real vwn3_energy(const FunDensProp* dp);
static void vwn3_first(FunFirstFuncDrv *ds,   real factor, const FunDensProp* dp);
static void vwn3_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp);
static void vwn3_third(FunThirdFuncDrv *ds,   real factor, const FunDensProp* dp);
static real vwn_energy(const FunDensProp* dp);
static void vwn_first(FunFirstFuncDrv *ds,   real factor, const FunDensProp* dp);
static void vwn_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp);
static void vwn_third(FunThirdFuncDrv *ds,   real factor, const FunDensProp* dp);
static real vwni_energy(const FunDensProp* dp);
static void vwni_first(FunFirstFuncDrv *ds,   real factor, const FunDensProp* dp);
static void vwni_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp);
static void vwni_third(FunThirdFuncDrv *ds,   real factor, const FunDensProp* dp);

/* VWN3 is a Gaussian version of the VWN functional based on suboptimal
 * set of parameters */
Functional VWN3Functional = {
    "VWN3",      /* name */
    vwn_isgga,  /* gga-corrected */
    vwn_read,   /* no extra input expected, just set the common block */
    NULL,
    vwn3_energy, 
    vwn3_first,
    vwn3_second,
    vwn3_third
};

Functional VWN5Functional = {
    "VWN5",     /* name */
    vwn_isgga,  /* gga-corrected */
    vwn_read,   /* no extra input expected, just set the common block */
    NULL,
    vwn_energy, 
    vwn_first,
    vwn_second,
    vwn_third
};

/* VWN is used for backward compatibility only */
Functional VWNFunctional = {
    "VWN",     /* name */
    vwn_isgga,  /* gga-corrected */
    vwn_read,   /* no extra input expected, just set the common block */
    NULL,
    vwn_energy, 
    vwn_first,
    vwn_second,
    vwn_third
};

/* VWNIFunctional is a variant of VWN5 functional with another spin
   polarization dependence:

   F(r,zeta) = (E_p + f(zeta)*(E_f - E_p))*rho

   The implementation is crippled and works only for closed shell
   systems up to linear response.
 */
Functional VWNIFunctional = {
    "VWNI",      /* name */
    vwn_isgga,  /* gga-corrected */
    vwn_read,   /* no extra input expected, just set the common block */
    NULL,
    vwni_energy, 
    vwni_first,
    vwni_second,
    vwni_third
};


/* IMPLEMENTATION PART */
#define VWN_ZERO 1e-40

static int
vwn_read(const char* conf_line)
{
    fun_set_hf_weight(0);
    return 1;
}

/* vwn_params contains two sets of parameters for paramagnetic and
   ferromagnetic cases.  See Table 5 in VWN paper.
*/


static const struct vwn_params {
    real X0, A, B, C;
}   vwn_paramagnetic  = { -0.1049800, 0.0621814, 3.72744, 12.9352 },
    vwn_ferromagnetic = { -0.3250000, 0.0310907, 7.06042, 18.0578 },
    vwn_interp        = { -0.0047584,-0.0337737, 1.13107, 13.0045 },
    vwn3_paramagnetic = { -0.4092860, 0.0621814, 13.0720, 42.7198 },
    vwn3_ferromagnetic= { -0.7432940, 0.0310907, 20.1231, 101.578 };
            
        
static const real SPINPOLF    = 1.92366105093154; /* 1/(2(2^(1/3)-1)) */
static const real THREEFTHRD2 = 0.584822305543806;/* hm? 4.5/(4*SPINPOLF) */
static const real FOURTHREE   = 1.333333333333333;

/* vwn_en_pot:
   returns  "order" numbers in enpot array
   enpot[0]: energy of given type p E = rho F - THIS IS AN EXCEPTION!
                                                DO NOT BLAME IT ON ME.
   enpot[1]: E'=enpot[0] + rho*F'
   enpot[2]: E''
   enpot[3]: E'''
*/
static void
vwn_en_pot(real* enpot, real rho, int order, const struct vwn_params* p)
{
    const real
        AI   = p->A,
        BI   = p->B,
        CI   = p->C,
        X0I  = p->X0;
    const real 
        Q    = sqrt(4*CI - BI*BI),
        XF0I = X0I*X0I + BI*X0I + CI,
        YF0I = Q/(BI + 2*X0I),
        DCRS = pow(3.0/(4*M_PI),1.0/6.0),
        B    = X0I/XF0I,
        C    = XF0I*YF0I,
        ACON = B*BI - 1.0,
        BCON = 2*ACON + 2.0,
        CCON = 2*BI*(1.0/Q - X0I/C);

    real rho13 = pow(rho,1.0/3.0);
    real x     = DCRS/sqrt(rho13);
    real xrho  = -DCRS*pow(rho,-7.0/6.0)/6.0;
    real xxrho = +DCRS*pow(rho,-13.0/6.0)*7.0/36.0;
    real xf   = x*x + BI*x+CI;
    real xfx  = 2*x + BI;
    real yf   = Q/xfx;
    real e1, ex1, exx1, exxx1;
    e1  = 2*log(x) 
        + ACON*log(xf)
        - BCON*log(x - X0I)
        + CCON*atan(yf);
    enpot[0] = 0.5*AI*e1;
    if(order<1) return;

    ex1 = 2.0/x 
        + ACON*xfx/xf
        - BCON/(x - X0I)
        - CCON*(2*yf/xfx)/(1.0 + yf*yf);
    enpot[1] = 0.5*AI*(e1 + rho*ex1*xrho);
    if(order<2) return;

    exx1= -2.0/(x*x) 
        + ACON*(2.0/xf - xfx*xfx/(xf*xf))
        + BCON/((x - X0I)*(x - X0I))
        + CCON*8*Q*xfx/((Q*Q + xfx*xfx)*(Q*Q + xfx*xfx));
    enpot[2] = 0.5*AI*xrho*(2*ex1 + rho*exx1*xrho - ex1*7.0/6.0);
    if(order<3) return;
  
    exxx1= 4.0/(x*x*x)
        + ACON*(2.0*xfx/(xf*xf))*(xfx*xfx/xf-3.0)
        - BCON*2.0/pow(x - X0I,3.0)
        + CCON*16.0*Q*(Q*Q-3.0*xfx*xfx)/pow(Q*Q + xfx*xfx,3.0);
    enpot[3] = 0.5*AI*xxrho*(2*ex1 + rho*exx1*xrho - 7.0/6.0*ex1)
        + 0.5*AI*xrho*(2*exx1*xrho 
                       +exx1*xrho + rho*(exxx1*xrho*xrho+exx1*xxrho)
                       -7.0/6.0*exx1*xrho);
}

static real
par_energy(const FunDensProp* dp, const struct vwn_params* para,
           const struct vwn_params* ferro)
{
    real ep_p[2], ep_f[2], ep_i[2], zeta, zeta4, f_zeta, delta;
    real rhoa = dp->rhoa, rhob = dp->rhob, rho;

    if(rhoa<VWN_ZERO) rhoa = VWN_ZERO;
    if(rhob<VWN_ZERO) rhob = VWN_ZERO;
    rho = rhoa + rhob;
    vwn_en_pot(ep_p, rho, 0, para);

    if(fabs(dp->rhoa-dp->rhob)<VWN_ZERO) return ep_p[0]*rho;
    vwn_en_pot(ep_f, rho, 0, ferro);
    vwn_en_pot(ep_i, rho, 0, &vwn_interp);

    zeta   = (dp->rhoa-dp->rhob)/rho;
    zeta4  = pow(zeta,4.0);
    f_zeta = SPINPOLF*(pow(1+zeta,FOURTHREE)+pow(1-zeta,FOURTHREE)-2.0);
    delta  = f_zeta*((ep_f[0]-ep_p[0])*zeta4 + ep_i[0]*(1-zeta4)*THREEFTHRD2);
    
    return (ep_p[0]+ delta)*rho;
}

static void
par_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp,
          const struct vwn_params* para, const struct vwn_params* ferro)
{
    real zeta, zeta3,zeta4, f_zeta, f_zet1, e_f,acp, dacp, vcfp, g_f;
    real delta, ep_p[2], ep_f[2], ep_i[2];
    real rhoa = dp->rhoa, rhob = dp->rhob, rho;

    if(rhoa<VWN_ZERO) rhoa = VWN_ZERO;
    if(rhob<VWN_ZERO) rhob = VWN_ZERO;
    rho = rhoa + rhob;
    vwn_en_pot(ep_p, rho, 1, para);

    ds->df1000 += ep_p[1]*factor;
    ds->df0100 += ep_p[1]*factor;

    if(fabs(dp->rhoa-dp->rhob)<VWN_ZERO) return;

    /* contribution from spin-polarized case; first order */
    zeta   = (dp->rhoa-dp->rhob)/rho;
    zeta3  = pow(zeta,3.0);
    zeta4  = zeta3*zeta;
    f_zeta = SPINPOLF*(pow(1+zeta,FOURTHREE)+pow(1-zeta,FOURTHREE)-2.0);
    f_zet1 = SPINPOLF*4.0/3.0*(pow(1+zeta,1.0/3.0)-pow(1-zeta,1.0/3.0));
    vwn_en_pot(ep_f, rho, 1, ferro);
    e_f    = ep_f[0] - ep_p[0];
    g_f    = ep_f[1] - ep_p[1];
    vwn_en_pot(ep_i, rho, 1, &vwn_interp);
    acp    = ep_i[0]*THREEFTHRD2;
    dacp   = ep_i[1]*THREEFTHRD2;

    vcfp = f_zeta*(g_f*zeta4 + dacp*(1-zeta4));
    delta= (f_zet1*(e_f*zeta4 + acp*(1-zeta4)) +
            4*f_zeta*(e_f - acp)*zeta3);

    /* the final section: begin */
    ds->df1000 += (vcfp + delta*(1-zeta))*factor;
    ds->df0100 += (vcfp - delta*(1+zeta))*factor; 
    /* the final section: end */
}

/* vwn_second:
   CAUTION: may raise zeros to a negative power!
*/
static void
par_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp,
           const struct vwn_params* para, const struct vwn_params* ferro)
{
    real zeta, zeta2, zeta3,zeta4, f_zeta, f_zet1, f_zet2, vcfp;
    real delta, ep_p[3], ep_f[3], ep_i[3];
    real vcfp1, ef0, ef1, ef2, ei0, ei1, ei2, bterm, cterm, dterm;
    real spA, spB;
    real rhoa = dp->rhoa, rhob = dp->rhob, rho = dp->rhoa + dp->rhob;
    real rho2 = rho*rho;
    real dAA, dAB, dBB;

    vwn_en_pot(ep_p, rho, 2, para);

    ds->df1000 += ep_p[1]*factor;
    ds->df0100 += ep_p[1]*factor;
    ds->df2000 += ep_p[2]*factor;
    ds->df0200 += ep_p[2]*factor;
    ds->df1100 += ep_p[2]*factor;

    /* if(0&&dp->rhoa==dp->rhob) return; */
    /* contribution from spin-polarized case; second order */
    zeta   = (dp->rhoa-dp->rhob)/rho;
    zeta2  = zeta*zeta;
    zeta3  = zeta2*zeta;
    zeta4  = zeta3*zeta;
    f_zeta = SPINPOLF*(pow(1+zeta,FOURTHREE)+pow(1-zeta,FOURTHREE)-2.0);
    f_zet1 = SPINPOLF*4.0/3.0*(pow(1+zeta,1.0/3.0)-pow(1-zeta,1.0/3.0));
    /* CAUTION: may raise 0 to negative power! */
    f_zet2 = SPINPOLF*4.0/9.0*(pow(1+zeta,-2.0/3.0)+pow(1-zeta,-2.0/3.0));
    vwn_en_pot(ep_f, rho, 2, ferro);
    ef0   = ep_f[0] - ep_p[0];
    ef1   = ep_f[1] - ep_p[1];
    ef2   = ep_f[2] - ep_p[2];
    vwn_en_pot(ep_i, rho, 2, &vwn_interp);
    ei0   = ep_i[0]*THREEFTHRD2;
    ei1   = ep_i[1]*THREEFTHRD2;
    ei2   = ep_i[2]*THREEFTHRD2;

    bterm = ef1*zeta4 + ei1*(1-zeta4);
    vcfp  = f_zeta*bterm;
    delta = (f_zet1*(ef0*zeta4 + ei0*(1-zeta4)) 
             + 4*f_zeta*(ef0 - ei0)*zeta3);

    spA = 2*rhob/rho2; /* =  2(1-zeta)/rho */
    spB =-2*rhoa/rho2; /* = -2(1+zeta)/rho */
    /* contribution from spin-polarized case; second order */
    /* spin independent part of vcfp */
    vcfp1 = f_zeta*(ef2*zeta4 + ei2*(1-zeta4));
    /* spin dependent part of vcfp */
    cterm = 4*f_zeta*(ef1-ei1)*zeta3 + bterm*f_zet1 - delta;

    /* spin dependent part of delta */
    dterm = (f_zet2*(ef0*zeta4+ei0*(1-zeta4))
             +8*f_zet1*(ef0-ei0)*zeta3
             +12*f_zeta*(ef0-ei0)*zeta2)*rho;

    dAA =  dterm*spA*spA;
    dAB =  dterm*spA*spB;
    dBB =  dterm*spB*spB;

    /* the final section: begin */
    ds->df1000 += (vcfp + delta*(1-zeta))*factor;
    ds->df0100 += (vcfp - delta*(1+zeta))*factor;

    ds->df2000 += (vcfp1+ cterm*(spA+spA) + dAA)*factor;
    ds->df1100 += (vcfp1+ cterm*(spA+spB) + dAB)*factor;
    ds->df0200 += (vcfp1+ cterm*(spB+spB) + dBB)*factor;
    /* the final section: end */
}

/* third not tested for open-shell! */
static void
par_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp,
          const struct vwn_params* para, const struct vwn_params* ferro)
{
    real zeta, zeta2, zeta3,zeta4, f_zeta, f_zet1, f_zet2, f_zet3, vcfp;
    real delta, ep_p[4], ep_f[4], ep_i[4];
    real vcfp1, vcfp2, ef0, ef1, ef2, ef3, ei0, ei1, ei2, ei3;
    real bterm, cterm, dterm, eterm, ctrm1, ctrm2, dtrm1, dtrm2;
    real ef2bi, ei2bi;
    real spA, spB, spAA, spAB, spBB;
    real rhoa = dp->rhoa, rhob = dp->rhob, rho, rho2, rho3;

    if(rhoa<VWN_ZERO) rhoa = VWN_ZERO;
    if(rhob<VWN_ZERO) rhob = VWN_ZERO;
    rho = rhoa + rhob;
    rho2 = rho*rho; rho3 = rho2*rho;
    
    vwn_en_pot(ep_p, rho, 3, para);

    ds->df1000 += ep_p[1]*factor;
    ds->df0100 += ep_p[1]*factor;
    ds->df2000 += ep_p[2]*factor;
    ds->df0200 += ep_p[2]*factor;
    ds->df1100 += ep_p[2]*factor;

    ds->df3000 += ep_p[3]*factor;
    ds->df2100 += ep_p[3]*factor;
    ds->df1200 += ep_p[3]*factor;
    ds->df0300 += ep_p[3]*factor;

    /* if(0&&dp->rhoa==dp->rhob) return; */
    /* contribution from spin-polarized case; second order */
    zeta   = (dp->rhoa-dp->rhob)/rho;
    zeta2  = zeta*zeta;
    zeta3  = zeta2*zeta;
    zeta4  = zeta3*zeta;
    f_zeta = SPINPOLF*    (pow(1+zeta,FOURTHREE)+pow(1-zeta,FOURTHREE)-2.0);
    f_zet1 = SPINPOLF*4.0/3.0 *(pow(1+zeta, 1.0/3.0)-pow(1-zeta, 1.0/3.0));
    f_zet2 = SPINPOLF*4.0/9.0 *(pow(1+zeta,-2.0/3.0)+pow(1-zeta,-2.0/3.0));
    f_zet3 =-SPINPOLF*8.0/27.0*(pow(1+zeta,-5.0/3.0)-pow(1-zeta,-5.0/3.0));

    vwn_en_pot(ep_f, rho, 3, ferro);
    ef0   = ep_f[0] - ep_p[0];
    ef1   = ep_f[1] - ep_p[1];
    ef2   = ep_f[2] - ep_p[2];
    ef3   = ep_f[3] - ep_p[3];
    vwn_en_pot(ep_i, rho, 3,&vwn_interp);
    ei0   = ep_i[0]*THREEFTHRD2;
    ei1   = ep_i[1]*THREEFTHRD2;
    ei2   = ep_i[2]*THREEFTHRD2;
    ei3   = ep_i[3]*THREEFTHRD2;

    bterm = ef1*zeta4 + ei1*(1-zeta4);
    vcfp  = f_zeta*bterm;
    delta = (f_zet1*(ef0*zeta4 + ei0*(1-zeta4)) 
             + 4*f_zeta*(ef0 - ei0)*zeta3);

    spA = 2*rhob/rho2; /* =  2(1-zeta)/rho */
    spB =-2*rhoa/rho2; /* = -2(1+zeta)/rho */
    spAA = -4*rhob/rho3;
    spAB = 2*(rho-2*rhob)/rho3;
    spBB = 4*rhoa/rho3;
    /* contribution from spin-polarized case; second order */
    /* spin independent part of vcfp */
    vcfp1 = f_zeta*(ef2*zeta4 + ei2*(1-zeta4));
    /* spin dependent part of vcfp */
    cterm = 4*f_zeta*(ef1-ei1)*zeta3 + bterm*f_zet1 - delta;

    /* spin dependent part of delta */
    dterm = (f_zet2*(ef0*zeta4+ei0*(1-zeta4))
             +8*f_zet1*(ef0-ei0)*zeta3
             +12*f_zeta*(ef0-ei0)*zeta2)*rho;

    /* third order terms */
    vcfp2 = f_zeta*(ef3*zeta4+ei3*(1-zeta4));
    eterm = f_zet1*(ef2*zeta4 + ei2*(1-zeta4)) + f_zeta*(ef2-ei2)*4*zeta3;
    ef2bi = ef2-(ef1-ef0)/rho;
    ei2bi = ei2-(ei1-ei0)/rho;
    ctrm1 = 4*f_zeta*(ef2bi-ei2bi)*zeta3 +f_zet1*(ef2bi*zeta4+ei2bi*(1-zeta4));

    ctrm2 = (ef1-ei1-ef0+ei0)*(8*f_zet1*zeta3+12*f_zeta*zeta2)
        +f_zet2*(bterm-(ef0*zeta4 + ei0*(1-zeta4)));

    dtrm1 = f_zet2*((ef1-ef0)*zeta4+(ei1-ei0)*(1-zeta4))
        +(8*f_zet1*zeta3+12*f_zeta*zeta2)*(ef1-ei1-ef0+ei0)
        +dterm/rho;
    dtrm2 = ((12*f_zet2*zeta3 + 36*f_zet1*zeta2 + 24*f_zeta*zeta)*(ef0-ei0)+
             f_zet3*(ef0*zeta4+ei0*(1-zeta4)))*rho;

    /* the final section: begin */
    ds->df1000 += (vcfp + delta*(1-zeta))*factor;
    ds->df0100 += (vcfp - delta*(1+zeta))*factor;

    ds->df2000 += (vcfp1+ cterm*(spA+spA) + dterm*spA*spA)*factor;
    ds->df1100 += (vcfp1+ cterm*(spA+spB) + dterm*spA*spB)*factor;
    ds->df0200 += (vcfp1+ cterm*(spB+spB) + dterm*spB*spB)*factor;
    ds->df3000 += (vcfp2+ eterm*spA + 
                   ctrm1*(spA+spA)+ ctrm2*spA*(spA+spA) + cterm*(spAA+spAA) +
                   dtrm1*(spA*spA)+ dtrm2*spA*spA*spA   + dterm*(2*spAA*spA)
        )*factor;
    ds->df2100 += (vcfp2+ eterm*spB + 
                   ctrm1*(spA+spA)+ ctrm2*spB*(spA+spA) + cterm*(spAB+spAB) +
                   dtrm1*(spA*spA)+ dtrm2*spA*spA*spB   + dterm*(2*spAB*spA)
        )*factor;
    ds->df1200 += (vcfp2+ eterm*spA + 
                   ctrm1*(spB+spB)+ ctrm2*spA*(spB+spB) + cterm*(spAB+spAB) +
                   dtrm1*(spB*spB)+ dtrm2*spB*spB*spA   + dterm*(2*spAB*spB)
        )*factor;
    ds->df0300 += (vcfp2+ eterm*spB + 
                   ctrm1*(spB+spB)+ ctrm2*spB*(spB+spB) + cterm*(spBB+spBB) +
                   dtrm1*(spB*spB)+ dtrm2*spB*spB*spB   + dterm*(2*spBB*spB)
        )*factor;
    /* the final section: end */
}

/* The dispatch part of the functional implementation */
static real
vwn3_energy(const FunDensProp* dp)
{
    return par_energy(dp, &vwn3_paramagnetic, &vwn3_ferromagnetic);
}

static void
vwn3_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
    par_first(ds, factor, dp, &vwn3_paramagnetic, &vwn3_ferromagnetic);
}

static void
vwn3_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    par_second(ds, factor, dp, &vwn3_paramagnetic, &vwn3_ferromagnetic);
}

static void
vwn3_third(FunThirdFuncDrv *ds,   real factor, const FunDensProp* dp)
{
    par_third(ds, factor, dp, &vwn3_paramagnetic, &vwn3_ferromagnetic);
}


static real
vwn_energy(const FunDensProp* dp)
{
    return par_energy(dp, &vwn_paramagnetic, &vwn_ferromagnetic);
}

static void
vwn_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
    par_first(ds, factor, dp, &vwn_paramagnetic, &vwn_ferromagnetic);
}

static void
vwn_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    par_second(ds, factor, dp, &vwn_paramagnetic, &vwn_ferromagnetic);
}

static void
vwn_third(FunThirdFuncDrv *ds,   real factor, const FunDensProp* dp)
{
    par_third(ds, factor, dp, &vwn_paramagnetic, &vwn_ferromagnetic);
}


/* Other spin interpolation scheme */
static real
spni_energy(const FunDensProp* dp, const struct vwn_params* para,
            const struct vwn_params* ferro)
{
    real ep_p[2], ep_f[2], ep_i[2], zeta, zeta4, f_zeta, delta;
    real rhoa = dp->rhoa, rhob = dp->rhob, rho;

    if(rhoa<VWN_ZERO) rhoa = 1e-40;
    if(rhob<VWN_ZERO) rhob = 1e-40;
    rho = rhoa + rhob;
    vwn_en_pot(ep_p, rho, 0, para);

    if( fabs(dp->rhoa - dp->rhob)<VWN_ZERO) return ep_p[0]*rho;
    vwn_en_pot(ep_f, rho, 0, ferro);
    vwn_en_pot(ep_i, rho, 0, &vwn_interp);

    zeta   = (dp->rhoa-dp->rhob)/rho;
    zeta4  = pow(zeta,4.0);
    f_zeta = SPINPOLF*(pow(1+zeta,FOURTHREE)+pow(1-zeta,FOURTHREE)-2.0);
    delta  = f_zeta*(ep_f[0]-ep_p[0]);
    
    return (ep_p[0]+ delta)*rho;
}

static void
spni_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp,
           const struct vwn_params* para, const struct vwn_params* ferro)
{
    real zeta, f_zeta, f_zet1, vcfp;
    real delta, ep_p[2], ep_f[2];
    real rhoa = dp->rhoa, rhob = dp->rhob, rho;

    if(rhoa<VWN_ZERO) rhoa = 1e-40;
    if(rhob<VWN_ZERO) rhob = 1e-40;
    rho = rhoa + rhob;
    vwn_en_pot(ep_p, rho, 1, para);

    ds->df1000 += ep_p[1]*factor;
    ds->df0100 += ep_p[1]*factor;

    /* if(dp->rhoa==dp->rhob) return; */

    /* contribution from spin-polarized case; first order */
    zeta   = (dp->rhoa-dp->rhob)/rho;
    f_zeta = SPINPOLF*(pow(1+zeta,FOURTHREE)+pow(1-zeta,FOURTHREE)-2.0);
    f_zet1 = SPINPOLF*4.0/3.0*(pow(1+zeta,1.0/3.0)-pow(1-zeta,1.0/3.0));
    vwn_en_pot(ep_f, rho, 1, ferro);

    vcfp = f_zeta*(ep_f[1] - ep_p[1]);
    delta= f_zet1*(ep_f[0] - ep_p[0]);

    /* the final section: begin */
    ds->df1000 += (vcfp + delta*(1-zeta))*factor;
    ds->df0100 += (vcfp - delta*(1+zeta))*factor; 
    /* the final section: end */
}

static void
spni_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp,
            const struct vwn_params* para, const struct vwn_params* ferro)
{
    real zeta, f_zeta, f_zet1, f_zet2, vcfp;
    real delta, ep_p[3], ep_f[3];
    real rhoa = dp->rhoa, rhob = dp->rhob, rho = dp->rhoa + dp->rhob;
    real rho2 = rho*rho;
    real rho3 = rho*rho2;
    real vcf2, fac2, vap2, del2, ef0, ef1, ef2;
    real zA, zB, zAAr, zABr, zBBr;

    vwn_en_pot(ep_p, rho, 2, para);

    ds->df1000 += ep_p[1]*factor;
    ds->df0100 += ep_p[1]*factor;
    ds->df2000 += ep_p[2]*factor;
    ds->df1100 += ep_p[2]*factor;
    ds->df0200 += ep_p[2]*factor;

    /* if( fabs(rhoa - rhob)<VWN_ZERO) return; */
    /* contribution from spin-polarized case; first order */
    zeta   = (rhoa - rhob)/rho;
    f_zeta = SPINPOLF*(pow(1+zeta,FOURTHREE)+pow(1-zeta,FOURTHREE)-2.0);
    f_zet1 = SPINPOLF*4.0/3.0*(pow(1+zeta, 1.0/3.0)-pow(1-zeta, 1.0/3.0));
    f_zet2 = SPINPOLF*4.0/9.0*(pow(1+zeta,-2.0/3.0)+pow(1-zeta,-2.0/3.0));
    vwn_en_pot(ep_f, rho, 2, ferro);

    ef0   = ep_f[0] - ep_p[0];
    ef1   = ep_f[1] - ep_p[1];
    ef2   = ep_f[2] - ep_p[2];
    vcfp = f_zeta*ef1;
    delta= f_zet1*ef0;

    vcf2 = f_zeta*ef2;
    vap2 = f_zet1*ef1;
    fac2 = f_zet2*ef0*rho;
    zA   =  2*rhob/rho2;
    zB   = -2*rhoa/rho2;
    zAAr = -4*rhob/rho2;
    zABr =  2*zeta/rho;
    zBBr =  4*rhoa/rho2;
    
    /* the final section: begin */
    ds->df1000 += (vcfp + delta*rho*zA)*factor;
    ds->df0100 += (vcfp - delta*(1+zeta))*factor; 

    ds->df2000 += (vcf2 + vap2*(zA+zA) + fac2*zA*zA + delta*zAAr)*factor;
    ds->df1100 += (vcf2 + vap2*(zA+zB)  +fac2*zA*zB + delta*zABr)*factor;
    ds->df0200 += (vcf2 + vap2*(zB+zB) + fac2*zB*zB + delta*zBBr)*factor;
    /* the final section: end */
}

static real
vwni_energy(const FunDensProp* dp)
{
    return spni_energy(dp, &vwn_paramagnetic, &vwn_ferromagnetic);
}

static void
vwni_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
    spni_first(ds, factor, dp, &vwn_paramagnetic, &vwn_ferromagnetic);
}

static void
vwni_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    spni_second(ds, factor, dp, &vwn_paramagnetic, &vwn_ferromagnetic);
}

static void
vwni_third(FunThirdFuncDrv *ds,   real factor, const FunDensProp* dp)
{
    fun_printf("vwni_third not implemented."); exit(1);
}
