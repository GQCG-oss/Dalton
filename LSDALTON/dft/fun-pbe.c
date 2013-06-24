/* fun-pbe.c:
 *
 *    implementation of the PBE exchange-correlation functional by
 *    Perdew, Burke and Ernzerhof J.P.Perdew, K. Burke, M. Ernzerhof,
 *    Phys. Rew. Letters, 77, 3865 (1996) [Ref 1]
 *
 *    The pbe functional and its derivatives up to fourth order.
 *    Implementation by B. Jansik, brano@theochem.kth.se, march 2003
 *    
 *    NOTE:
 *
 *    This functional is dependent on pw92 functional. It can not be
 *    separated, pw92 is intrinsically entangled in PBE correlation
 *    energy NOTE: Improvement over PW91: Improvements over PW91
 *    include an accurate description of the linear response of the
 *    uniform electron gas, correct behavior under uniform scaling,
 *    and a smoother potential. [1]
 */

/* strictly conform to XOPEN ANSI C standard */
#if !defined(SYS_DEC)
/* XOPEN compliance is missing on old Tru64 4.0E Alphas and pow() prototype
 * is not specified. */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1
#endif

/* #define FOURTH_ORDER_DERIVATIVES */

#include <math.h>
#include <stdio.h>

#define __CVERSION__

#include "lsdalton_functionals.h"

/* INTERFACE PART */
static integer pbe_isgga(void) { return 0; }
static integer pbe_read(const char* conf_line, real *hfweight);
static real pbe_energy(const DftDensProp* dp);
static void pbe_first(FirstFuncDrv *ds,   real factor, const DftDensProp* dp);
static void pbe_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp);
static void pbe_third(ThirdFuncDrv *ds,   real factor, const DftDensProp* dp);
#ifdef FOURTH_ORDER_DERIVATIVES
static void pbe_fourth(FourthFuncDrv *ds, real factor, const DftDensProp* dp);
#endif

Functional PBEFunctional = {
    "PBE",       /* name */
    pbe_isgga,   /* gga-corrected */
    pbe_read, 
    NULL,
    pbe_energy, 
    pbe_first,
    pbe_second,
    pbe_third
#ifdef FOURTH_ORDER_DERIVATIVES
    ,pbe_fourth
#endif
};

/* Structure to store intermediate pbe exchange derivatives */
struct TwoVarDrv_ {
    real df10;
    real df01;
    real df20;
    real df11;
    real df02;
    real df30;
    real df21;
    real df12;
    real df03;
    real df40;
    real df31;
    real df22;
    real df13;
    real df04;
};

typedef struct TwoVarDrv_ TwoVarDrv;

/* Declarations of auxaliary functions
 * used for pbe functonal and its
 * derivatives */

static real pbe_cenergy(const DftDensProp* dp);
static real pbe_xenergy(const DftDensProp* dp);

static real pbe_ro2energy(const real *p_ro, const real *p_gro);
static void pbe_ro2d(TwoVarDrv *ds, real *p_ro, real *p_gro);

static void pbe_xd1(FirstFuncDrv *ds, real factor, const DftDensProp *dp);
static void pbe_cd1(FirstFuncDrv *ds, real factor, const DftDensProp *dp);
static void pbe_xd2(SecondFuncDrv *ds, real factor, const DftDensProp *dp);
static void pbe_cd2(SecondFuncDrv *ds, real factor, const DftDensProp *dp);
static void pbe_xd3(ThirdFuncDrv *ds, real factor, const DftDensProp *dp);
	
static void pbe_cd3(ThirdFuncDrv *ds, real factor, const DftDensProp *dp);
#ifdef FOURTH_ORDER_DERIVATIVES
static void pbe_xd4(FourthFuncDrv *ds, real factor, const DftDensProp *dp);
static void pbe_cd4(FourthFuncDrv *ds, real factor, const DftDensProp *dp);
#endif

/* End of declarations */


/* IMPLEMENTATION PART */
static integer
pbe_read(const char* conf_line, real *hfweight)
{
    return 1;
}

#define PBETHR 1e-14
static real
pbe_energy(const DftDensProp* dp)
{
    if ((dp->rhoa + dp->rhob) < PBETHR) return(0.0);
    return (
        pbe_xenergy(dp) +
        pbe_cenergy(dp) 
        );
}

static void
pbe_first(FirstFuncDrv *ds, real factor, const DftDensProp* dp)
{
    if ((dp->rhoa + dp->rhob) < PBETHR) return;
    pbe_xd1(ds, factor, dp);
    pbe_cd1(ds, factor, dp);
}
static void
pbe_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp)
{
    if ((dp->rhoa + dp->rhob) < PBETHR) return;
    pbe_xd2(ds, factor, dp);
    pbe_cd2(ds, factor, dp);
}

static void
pbe_third(ThirdFuncDrv *ds, real factor, const DftDensProp* dp)
{
    if ((dp->rhoa + dp->rhob) < PBETHR) return;
    pbe_xd3(ds, factor, dp);
    pbe_cd3(ds, factor, dp);
}

#ifdef FOURTH_ORDER_DERIVATIVES
static void
pbe_fourth(FourthFuncDrv *ds, real factor, const DftDensProp *dp)
{
    if ((dp->rhoa + dp->rhob) < PBETHR) return;
    pbe_xd4(ds, factor, dp);
    pbe_cd4(ds, factor, dp);
}
#endif

/* Constants for PBE, revPBE and RPBE Exchange energy functional
 *
 * PBE
 * See: J.P.Perdew, K. Burke, M. Ernzerhof, Phys. Rew. Letters, 77,
 * 3865 (1996) 
 */
/* exchange-related constants */
#define FAC_BETA 0.06672455060314922
static const real EXC_kappa   =    0.804;
static const real EXC_beta    =    FAC_BETA;
static const real EXC_micro   =    FAC_BETA*(M_PI*M_PI/3.0); 
static const real EXC_nEXunif =   -0.9305257363491;
static const real EXC_nKF     =    3.89777708972075;

/* correlation-related constants */
#define FAC_GAMMA 0.03109069086965489503494086371273
static const real COR_gamma = FAC_GAMMA;
static const real COR_beta  = 0.06672455060314922;
static const real COR_delta = FAC_BETA/FAC_GAMMA;
static const real COR_nZ    = 1.5874010519682; /* (2^(2/3)) */
static const real COR_nKF   = 3.09366772628014; /* (3*pi^2)^(1/3) */
static const real COR_nKS   = 1.98468639521986; /* sqrt(4*nKF/pi) */
static const real COR_nT    = 0.03272492347489; /* 2/(nKS^6)       */

/* PBE functional exchange energy */
static real
pbe_xenergy(const DftDensProp* dp)
{
    real E;

    /* E=(roa*pbe_ro2energy(roa,groa) +
     * rob*pbe_ro2energy(rob,grob))/(roa+rob); */
    E= (dp->rhoa*pbe_ro2energy(&dp->rhoa, &dp->grada) + 
        dp->rhob*pbe_ro2energy(&dp->rhob, &dp->gradb))/(dp->rhoa+dp->rhob);

    return (dp->rhoa+dp->rhob)*E;
}

/* PBE functional correlation energy */

static real
pbe_cenergy(const DftDensProp* dp)
{
    real  A;
    real  A_p2;
    real  E;
    real  EC;
    real  FI;
    real  FI_p3;
    real  H;
    real  H0;
    real  H1;
    real  KS;
    real  T;
    real  T_p2;
    real  T_p4;
    real  ZETA;
    real  MZETA, PZETA;
    real  groa;
    real  grob;
    real  roa;
    real  rob;
    real  roa_rob_p1f3, roa_rob_p1f6;
    real  roa_rob;

    /* Set up roa, rob, groa, grob */
    roa = dp->rhoa; 
    rob = dp->rhob;

    groa = dp->grada;
    grob = dp->gradb;

/* roa, rob and powers */
    roa_rob=roa+rob;
    roa_rob_p1f3=pow(roa_rob,(1.0/3.0));
    roa_rob_p1f6=pow(roa_rob,(1.0/6.0));

/* Auxaliary functions */

    ZETA=(roa-rob)/(roa_rob);
    KS=COR_nKS*roa_rob_p1f6;
    MZETA=pow((1.0-ZETA),(2.0/3.0));
    PZETA=pow((1.0+ZETA),(2.0/3.0));
    FI=(PZETA + MZETA)/2.0;
    T=(groa+grob)/((roa_rob)*(2.0*KS*FI));


    FI_p3=FI*FI*FI;
    T_p2=T*T;
    T_p4=T_p2*T_p2;

    EC  = PW92peFunctional.func(dp);

    H0  = COR_gamma*FI_p3;

    A   = COR_delta*( 1.0/(exp(-EC/H0) -1.0)  ) ;
    A_p2= A*A;

    H1=1.0 + COR_delta*T_p2*((1.0 + A*T_p2)/(1.0 + A*T_p2 + A_p2*T_p4));

    H=H0*log( H1 );
    E=EC+H;

    return(roa_rob*E);
}


static real
pbe_ro2energy(const real *p_ro, const real *p_gro)
{
    real  E;
    real  EXunif;
    real  Fxs;
    real  KF;
    real  S;
    real  S_p2;
    real  gro;
    real  ro;
    real  ro_p1f3;
    
/* Set up ro, gro */
    ro = *p_ro;
    gro = *p_gro;
  

/* Powers of ro */
    ro_p1f3=pow(ro,(1.0/3.0));

    KF=EXC_nKF*ro_p1f3;
    S=(1.0*gro)/(2.0*KF*ro);
    S_p2=S*S;

 /* EXunif=-3*KF/(4*pi); */
    EXunif= EXC_nEXunif*ro_p1f3;

    Fxs=1.0 + EXC_kappa - (EXC_kappa/( 1.0 + (EXC_micro*S_p2/EXC_kappa) ) );

    E=EXunif*Fxs;

    return(E);
}


static void 
pbe_ro2d(TwoVarDrv *ds, real *p_ro, real *p_gro)
{
    const real micro = EXC_micro, kappa = EXC_kappa;
    real  EXunif;
    real  EXunif_1;
    real  EXunif_2;
    real  EXunif_3;
    real  EXunif_4;
    real  Fxs;
    real  Fxs_01;
    real  Fxs_02;
    real  Fxs_03;
    real  Fxs_04;
    real  Fxs_10;
    real  Fxs_11;
    real  Fxs_12;
    real  Fxs_13;
    real  Fxs_20;
    real  Fxs_21;
    real  Fxs_22;
    real  Fxs_30;
    real  Fxs_31;
    real  Fxs_40;
    real  Fxs_term;
    real  Fxs_term_p2;
    real  Fxs_term_p3;
    real  Fxs_term_p4;
    real  Fxs_term_p5;
    real  KF;
    real  S;
    real  S_p2;
    real  gro;
    real  gro_p2;
    real  gro_p3;
    real  gro_p4;
    real  gro_p5;
    real  gro_p6;
    real  ro;
    real  ro_p10f3;
    real  ro_p11f3;
    real  ro_p13f3;
    real  ro_p16f3;
    real  ro_p1f3;
    real  ro_p2;
    real  ro_p2f3;
    real  ro_p4;
    real  ro_p4f3;
    real  ro_p5f3;
    real  ro_p6;
    real  ro_p7;
    real  ro_p7f3;
    real  ro_p8;
    real  ro_p8f3;
	
/* Powers of numbers */
    const real nKF_p2=EXC_nKF*EXC_nKF;
    const real nKF_p4=nKF_p2*nKF_p2;
    const real nKF_p6=nKF_p4*nKF_p2;

    const real kappa_p2=EXC_kappa*EXC_kappa;
    const real kappa_p3=kappa_p2*EXC_kappa;

    const real micro_p2=EXC_micro*EXC_micro;
    const real micro_p3=micro_p2*EXC_micro;

/* Set up ro, gro */

    ro = *p_ro;
    gro = *p_gro;

/* Powers of ro */
    ro_p2=ro*ro;
    ro_p4=ro_p2*ro_p2;
    ro_p6=ro_p4*ro_p2;
    ro_p7=ro_p6*ro;
    ro_p8=ro_p4*ro_p4;
    ro_p1f3=pow(ro,(1.0/3.0));
    ro_p2f3=ro_p1f3*ro_p1f3;
    ro_p4f3=ro_p1f3*ro;
    ro_p5f3=ro_p2f3*ro;
    ro_p7f3=ro_p2*ro_p1f3;
    ro_p8f3=ro_p2f3*ro_p2;
    ro_p10f3=ro_p7f3*ro;
    ro_p11f3=ro_p10f3*ro_p1f3;
    ro_p13f3=ro_p4*ro_p1f3;
    ro_p16f3=ro_p10f3*ro_p2;

/* Powers of gro */
    gro_p2=gro*gro;
    gro_p3=gro_p2*gro;
    gro_p4=gro_p2*gro_p2;
    gro_p5=gro_p4*gro;
    gro_p6=gro_p5*gro;

    KF=EXC_nKF*ro_p1f3;
    S=(1.0*gro)/(2.0*KF*ro);
    S_p2=S*S;

    EXunif=EXC_nEXunif*ro_p1f3;

/* Exunif derivatives */
    EXunif_1= ( 1.0/3.0)*EXC_nEXunif/ro_p2f3;
    EXunif_2=-( 2.0/9.0)*EXC_nEXunif/ro_p5f3;
    EXunif_3= (10.0/27.0)*EXC_nEXunif/ro_p8f3;
    EXunif_4=-(80.0/81.0)*EXC_nEXunif/ro_p11f3;

    Fxs=1.0 + kappa - (kappa/( 1.0 + (micro*S_p2/kappa) ) );

    Fxs_term=gro_p2*micro + 4.0*kappa*nKF_p2*ro_p8f3;
    Fxs_term_p2=Fxs_term*Fxs_term;
    Fxs_term_p3=Fxs_term_p2*Fxs_term;
    Fxs_term_p4=Fxs_term_p3*Fxs_term;
    Fxs_term_p5=Fxs_term_p4*Fxs_term;

/* Fxs function derivatives */

    Fxs_10=(-32.0/3.0)*gro_p2*kappa_p2*micro*nKF_p2*ro_p5f3/(Fxs_term_p2);
    Fxs_01=8.0*gro*kappa_p2*micro*nKF_p2*ro_p8f3/Fxs_term_p2;

    Fxs_20=-(32.0/9.0)*gro_p2*kappa_p2*micro*nKF_p2*ro_p2f3*
        (5.0*gro_p2*micro - 44.0*kappa*nKF_p2*ro_p8f3) / 
	(Fxs_term_p3);

    Fxs_02=8.0*kappa_p2*micro*nKF_p2*ro_p8f3*
        (-3.0*gro_p2*micro + 4.0*kappa*nKF_p2*ro_p8f3) / Fxs_term_p3;

    Fxs_11=64.0*gro*kappa_p2*micro*nKF_p2*ro_p5f3*
        (gro_p2*micro - 4.0*kappa*nKF_p2*ro_p8f3) / 
	(3.0*Fxs_term_p3);

    Fxs_30=-(64.0*gro_p2*kappa_p2*micro*nKF_p2*
             (5.0*gro_p4*micro_p2 - 440.0*gro_p2*kappa*micro*nKF_p2*ro_p8f3 +
              1232.0*kappa_p2*nKF_p4*ro_p16f3)) / 
	(27.0*ro_p1f3*Fxs_term_p4);

    Fxs_21=(64.0*gro*kappa_p2*micro*nKF_p2*ro_p2f3 *
            (5.0*gro_p4*micro_p2 - 128.0*gro_p2*kappa*micro*nKF_p2*ro_p8f3 + 
             176.0*kappa_p2*nKF_p4*ro_p16f3)) / 
	(9.0*Fxs_term_p4);

    Fxs_12=-(64.0*kappa_p2*micro*nKF_p2*ro_p5f3*
             (3.0*gro_p4*micro_p2 - 32.0*gro_p2*kappa*micro*nKF_p2*ro_p8f3 + 
              16.0*kappa_p2*nKF_p4*ro_p16f3)) / 
	(3.0*Fxs_term_p4);

    Fxs_03=96.0*gro*kappa_p2*micro_p2*nKF_p2*ro_p8f3*
        (gro_p2*micro -4.0*kappa*nKF_p2*ro_p8f3) / 
	Fxs_term_p4;

    Fxs_40=(64.0*gro_p2*kappa_p2*micro*nKF_p2*
            (5.0*gro_p6*micro_p3+3740.0*gro_p4*kappa*micro_p2*nKF_p2*ro_p8f3 - 
             62480.0*gro_p2*kappa_p2*micro*nKF_p4*ro_p16f3 + 
             83776.0*kappa_p3*nKF_p6*ro_p8)) / 
	(81.0*ro_p4f3*Fxs_term_p5);

    Fxs_31=(128.0*gro*kappa_p2*micro*nKF_p2*
            (5.0*gro_p6*micro_p3 - 940.0*gro_p4*micro_p2*kappa*nKF_p2*ro_p8f3 +
             7216.0*gro_p2*kappa_p2*micro*nKF_p4*ro_p16f3 - 
             4928.0*kappa_p3*nKF_p6*ro_p8)) / 
	(27.0*ro_p1f3*Fxs_term_p5);

    Fxs_22=(64.0*kappa_p2*micro*nKF_p2*ro_p2f3*
            (-15.0*gro_p6*micro_p3 +
             740.0*gro_p4*micro_p2*kappa*nKF_p2*ro_p8f3 - 
             2768.0*gro_p2*kappa_p2*micro*nKF_p4*ro_p16f3 + 
             704.0*kappa_p3*nKF_p6*ro_p8)) / 
	(9.0*Fxs_term_p5);

    Fxs_13=(256.0*gro*kappa_p2*micro_p2*nKF_p2*ro_p5f3*
            (gro_p4*micro_p2 - 20.0*gro_p2*kappa*micro*nKF_p2*ro_p8f3 +
             32.0*kappa_p2*nKF_p4*ro_p16f3)) / 
	Fxs_term_p5;

    Fxs_04=-(96.0*kappa_p2*micro_p2*nKF_p2*ro_p8f3*
             (5.0*gro_p4*micro_p2 - 40.0*gro_p2*kappa*micro*nKF_p2*ro_p8f3 +
              16.0*kappa_p2*nKF_p4*ro_p16f3)) / 
	Fxs_term_p5;

/* Derivatives */

    ds->df10 = Fxs*EXunif_1 + EXunif*Fxs_10;
    ds->df01 = EXunif*Fxs_01;

    ds->df20 = Fxs*EXunif_2 + 2.0*EXunif_1*Fxs_10 + EXunif*Fxs_20;
    ds->df11 = EXunif_1*Fxs_01 + EXunif*Fxs_11;
    ds->df02 = EXunif*Fxs_02;

    ds->df30 = Fxs*EXunif_3 + 3.0*Fxs_10*EXunif_2 + 
	3.0*Fxs_20*EXunif_1 + Fxs_30*EXunif;

    ds->df21 = EXunif_2*Fxs_01 + 2.0*EXunif_1*Fxs_11 + EXunif*Fxs_21;

    ds->df12 = EXunif_1*Fxs_02 + EXunif*Fxs_12;

    ds->df03 = EXunif*Fxs_03;

    ds->df40 = EXunif_4*Fxs + 4.0*EXunif_3*Fxs_10 + 6.0*EXunif_2*Fxs_20 + 
	4.0*EXunif_1*Fxs_30 + EXunif*Fxs_40;

    ds->df31 = EXunif_3*Fxs_01 + 3.0*EXunif_2*Fxs_11 + 3.0*EXunif_1*Fxs_21 + 
	EXunif*Fxs_31;

    ds->df22 = EXunif_2*Fxs_02 + 2.0*EXunif_1*Fxs_12 + EXunif*Fxs_22;

    ds->df13 = EXunif_1*Fxs_03 + EXunif*Fxs_13;

    ds->df04 = EXunif*Fxs_04;

}

static void
pbe_xd1(FirstFuncDrv *ds, real factor, const DftDensProp *dp)
{
real    roa = dp->rhoa;
real    rob = dp->rhob;

real    groa = dp->grada;
real    grob = dp->gradb;
	
real    Exa=pbe_ro2energy(&roa, &groa);
real    Exb=pbe_ro2energy(&rob, &grob);

TwoVarDrv dera;
TwoVarDrv derb;

    pbe_ro2d(&dera, &roa, &groa);
    pbe_ro2d(&derb, &rob, &grob);

real    Exa_10=dera.df10;
real    Exa_01=dera.df01;

real    Exb_10=derb.df10;
real    Exb_01=derb.df01;

/* Powers of total density */
real    roa_rob=roa+rob;
real    roa_rob_p2=roa_rob*roa_rob;

/* Energy */
real E = (roa*Exa + rob*Exb)/roa_rob;

/* Derivatives */

real  df10000  = (rob*Exa - rob*Exb + roa*roa_rob*Exa_10)/roa_rob_p2;
  ds->df1000 += factor*( E + roa_rob*df10000);

real  df01000  = (-roa*Exa + roa*Exb + rob*roa_rob*Exb_10)/roa_rob_p2;
  ds->df0100 += factor*( E + roa_rob*df01000);

real  df00100  = (roa*Exa_01)/roa_rob;
  ds->df0010 += factor*roa_rob*df00100;

real  df00010  =  (rob*Exb_01)/roa_rob;
  ds->df0001 += factor*roa_rob*df00010;

}

static void
pbe_xd2(SecondFuncDrv *ds, real factor, const DftDensProp *dp)
{
real    roa = dp->rhoa;
real    rob = dp->rhob;

real    groa = dp->grada;
real    grob = dp->gradb;
	
real    Exa=pbe_ro2energy(&roa, &groa);
real    Exb=pbe_ro2energy(&rob, &grob);

TwoVarDrv dera;
TwoVarDrv derb;

    pbe_ro2d(&dera, &roa, &groa);
    pbe_ro2d(&derb, &rob, &grob);

real    Exa_10=dera.df10;
real    Exa_01=dera.df01;
real    Exa_20=dera.df20;
real    Exa_11=dera.df11;
real    Exa_02=dera.df02;

real    Exb_10=derb.df10;
real    Exb_01=derb.df01;
real    Exb_20=derb.df20;
real    Exb_11=derb.df11;
real    Exb_02=derb.df02;

/* Powers of total density */
real    roa_rob=roa+rob;
real    roa_rob_p2=roa_rob*roa_rob;
real    roa_rob_p3=roa_rob_p2*roa_rob;

/* Energy */
real E = (roa*Exa + rob*Exb)/roa_rob;

/* Derivatives */

real  df10000  = (rob*Exa - rob*Exb + roa*roa_rob*Exa_10)/roa_rob_p2;
  ds->df1000 += factor*( E + roa_rob*df10000);

real  df01000  = (-roa*Exa + roa*Exb + rob*roa_rob*Exb_10)/roa_rob_p2;
  ds->df0100 += factor*( E + roa_rob*df01000);

real  df00100  = (roa*Exa_01)/roa_rob;
  ds->df0010 += factor*roa_rob*df00100;

real  df00010  =  (rob*Exb_01)/roa_rob;
  ds->df0001 += factor*roa_rob*df00010;

/* Second order derivatives */

real  df20000  =  (1.0/roa_rob_p3)*(-2.0*rob*Exa + 2.0*rob*Exb + 
                  roa_rob*(2.0*rob*Exa_10 + roa*roa_rob*Exa_20));
  ds->df2000 += factor*(2.0*df10000 + roa_rob*df20000); 

real  df11000  =  (1.0/roa_rob_p3)*((roa-rob)*Exa + (rob - roa)*Exb - 
                  roa_rob*(roa*Exa_10 + rob*Exb_10));
  ds->df1100 += factor*(df10000 + df01000 + roa_rob*df11000);

real  df10100  = (rob*Exa_01 + roa*roa_rob*Exa_11)/roa_rob_p2;
  ds->df1010 += factor*(df00100 + roa_rob*df10100); 

real  df10010  = -rob*Exb_01/roa_rob_p2;
  ds->df1001 += factor*(df00010 + roa_rob*df10010);

real  df02000  =  (1.0/roa_rob_p3)*(2.0*roa*Exa - 2.0*roa*Exb + 
                   roa_rob*(2.0*roa*Exb_10 + rob*roa_rob*Exb_20));
  ds->df0200 += factor*(2.0*df01000 + roa_rob*df02000);

real  df01100  =  -roa*Exa_01/roa_rob_p2;
  ds->df0110 += factor*(df00100 + roa_rob*df01100);

real  df01010  = (roa*Exb_01 + rob*roa_rob*Exb_11)/roa_rob_p2;
  ds->df0101 += factor*(df00010 + roa_rob*df01010);

real  df00200  =  roa*Exa_02/roa_rob;
  ds->df0020 += factor*roa_rob*df00200;

/* ds->df00110=0.0 */

real  df00020  = rob*Exb_02/roa_rob;
  ds->df0002 += factor*roa_rob*df00020;

}

static void
pbe_xd3(ThirdFuncDrv *ds, real factor, const DftDensProp *dp)
{
real    roa = dp->rhoa;
real    rob = dp->rhob;

real    groa = dp->grada;
real    grob = dp->gradb;
	
real    Exa=pbe_ro2energy(&roa, &groa);
real    Exb=pbe_ro2energy(&rob, &grob);

TwoVarDrv dera;
TwoVarDrv derb;

    pbe_ro2d(&dera, &roa, &groa);
    pbe_ro2d(&derb, &rob, &grob);

real    Exa_10=dera.df10;
real    Exa_01=dera.df01;
real    Exa_20=dera.df20;
real    Exa_11=dera.df11;
real    Exa_02=dera.df02;
real    Exa_30=dera.df30;
real    Exa_21=dera.df21;
real    Exa_12=dera.df12;
real    Exa_03=dera.df03;

real    Exb_10=derb.df10;
real    Exb_01=derb.df01;
real    Exb_20=derb.df20;
real    Exb_11=derb.df11;
real    Exb_02=derb.df02;
real    Exb_30=derb.df30;
real    Exb_21=derb.df21;
real    Exb_12=derb.df12;
real    Exb_03=derb.df03;

/* Powers of total density */
real    roa_rob=roa+rob;
real    roa_rob_p2=roa_rob*roa_rob;
real    roa_rob_p3=roa_rob_p2*roa_rob;
real    roa_rob_p4=roa_rob_p3*roa_rob;

/* Energy */
real E = (roa*Exa + rob*Exb)/roa_rob;

/* Derivatives */

real  df10000  = (rob*Exa - rob*Exb + roa*roa_rob*Exa_10)/roa_rob_p2;
  ds->df1000 += factor*( E + roa_rob*df10000);

real  df01000  = (-roa*Exa + roa*Exb + rob*roa_rob*Exb_10)/roa_rob_p2;
  ds->df0100 += factor*( E + roa_rob*df01000);

real  df00100  = (roa*Exa_01)/roa_rob;
  ds->df0010 += factor*roa_rob*df00100;

real  df00010  =  (rob*Exb_01)/roa_rob;
  ds->df0001 += factor*roa_rob*df00010;

/* Second order derivatives */

real  df20000  =  (1.0/roa_rob_p3)*(-2.0*rob*Exa + 2.0*rob*Exb + 
                  roa_rob*(2.0*rob*Exa_10 + roa*roa_rob*Exa_20));
  ds->df2000 += factor*(2.0*df10000 + roa_rob*df20000); 

real  df11000  =  (1.0/roa_rob_p3)*((roa-rob)*Exa + (rob - roa)*Exb - 
                  roa_rob*(roa*Exa_10 + rob*Exb_10));
  ds->df1100 += factor*(df10000 + df01000 + roa_rob*df11000);

real  df10100  = (rob*Exa_01 + roa*roa_rob*Exa_11)/roa_rob_p2;
  ds->df1010 += factor*(df00100 + roa_rob*df10100); 

real  df10010  = -rob*Exb_01/roa_rob_p2;
  ds->df1001 += factor*(df00010 + roa_rob*df10010);

real  df02000  =  (1.0/roa_rob_p3)*(2.0*roa*Exa - 2.0*roa*Exb + 
                   roa_rob*(2.0*roa*Exb_10 + rob*roa_rob*Exb_20));
  ds->df0200 += factor*(2.0*df01000 + roa_rob*df02000);

real  df01100  =  -roa*Exa_01/roa_rob_p2;
  ds->df0110 += factor*(df00100 + roa_rob*df01100);

real  df01010  = (roa*Exb_01 + rob*roa_rob*Exb_11)/roa_rob_p2;
  ds->df0101 += factor*(df00010 + roa_rob*df01010);

real  df00200  =  roa*Exa_02/roa_rob;
  ds->df0020 += factor*roa_rob*df00200;

/* ds->df00110=0.0 */

real  df00020  = rob*Exb_02/roa_rob;
  ds->df0002 += factor*roa_rob*df00020;

/* Third order derivatives */

real  df30000  =  (1.0/roa_rob_p4)*(6.0*rob*Exa - 6.0*rob*Exb + 
                   roa_rob*(-6.0*rob*Exa_10 +   roa_rob*(3.0*rob*Exa_20 + 
                   roa*roa_rob*Exa_30)));
  ds->df3000 += factor*(3.0*df20000 + roa_rob*df30000);

real  df21000  =  (1.0/roa_rob_p4)*(-2.0*(roa-2.0*rob)*Exa + 
                   2.0*(roa-2.0*rob)*Exb - 
                   roa_rob*(-2.0*(roa-rob)*Exa_10 - 
                   2.0*rob*Exb_10 + roa*roa_rob*Exa_20));
  ds->df2100 += factor*(2*df11000 + df20000 + roa_rob*df21000);

real  df20100  =  (1.0/roa_rob_p3)*(-2.0*rob*Exa_01 + roa_rob*(2.0*rob*Exa_11 + 
                   roa*roa_rob*Exa_21));
  ds->df2010 += factor*(2.0*df10100 + roa_rob*df20100);


real  df20010  =  2.0*rob*Exb_01/roa_rob_p3;
  ds->df2001 += factor*(2.0*df10010 + roa_rob*df20010);

real  df12000  = (1.0/roa_rob_p4)*(-2.0*(2.0*roa-rob)*Exa + 
                  2.0*(2.0*roa-rob)*Exb - 
                  roa_rob*(-2.0*roa*Exa_10 + 
                  2.0*(roa-rob)*Exb_10 + rob*roa_rob*Exb_20));
  ds->df1200 += factor*(2.0*df11000 + df02000 + roa_rob*df12000);


real  df11100  = ((roa-rob)*Exa_01 - roa*roa_rob*Exa_11)/roa_rob_p3;
  ds->df1110 += factor*(df10100 + df01100 + roa_rob*df11100);

real  df11010  = ((rob-roa)*Exb_01 - rob*roa_rob*Exb_11)/roa_rob_p3;
  ds->df1101 += factor*(df10010 + df01010 + roa_rob*df11010);

real  df10200   = (rob*Exa_02 + roa*roa_rob*Exa_12)/roa_rob_p2;
  ds->df1020 += factor*(df00200 + roa_rob*df10200);

/* ds->df10110=0.0 */

real  df10020  =  -rob*Exb_02/roa_rob_p2;
  ds->df1002 += factor*(df00020 + roa_rob*df10020);

real  df03000  = (1.0/roa_rob_p4)*(-6.0*roa*Exa + 6.0*roa*Exb + 
                  roa_rob*(-6.0*roa*Exb_10 + 
                  roa_rob*(3.0*roa*Exb_20 + 
                  rob*roa_rob*Exb_30)));
  ds->df0300 += factor*(3.0*df02000 + roa_rob*df03000);   

real  df02100  =       2.0*roa*Exa_01/roa_rob_p3;
  ds->df0210 += factor*(2.0*df01100 + roa_rob*df02100);

real  df02010  = (1.0/roa_rob_p3)*(-2.0*roa*Exb_01 + roa_rob*(2.0*roa*Exb_11 + 
                  rob*roa_rob*Exb_21));
  ds->df0201 += factor*(2.0*df01010 + roa_rob*df02010);

real  df01200  =  -roa*Exa_02/roa_rob_p2;
  ds->df0120 += factor*(df00200 + roa_rob*df01200);

/* ds->df01110=0.0 */

real  df01020  =  (roa*Exb_02 + rob*roa_rob*Exb_12)/roa_rob_p2;
  ds->df0102 += factor*(df00020 + roa_rob*df01020);


real  df00300  =  roa*Exa_03/roa_rob;
  ds->df0030 += factor*roa_rob*df00300;

/* ds->df00210=0.0 */

/* ds->df00120=0.0 */

real  df00030  = rob*Exb_03/roa_rob;
  ds->df0003 += factor*roa_rob*df00030;

}

#ifdef FOURTH_ORDER_DERIVATIVES
static void
pbe_xd4(FourthFuncDrv *ds, real factor, const DftDensProp *dp)
{
real    roa = dp->rhoa;
real    rob = dp->rhob;

real    groa = dp->grada;
real    grob = dp->gradb;
	
real    Exa=pbe_ro2energy(&roa, &groa);
real    Exb=pbe_ro2energy(&rob, &grob);

TwoVarDrv dera;
TwoVarDrv derb;

    pbe_ro2d(&dera, &roa, &groa);
    pbe_ro2d(&derb, &rob, &grob);

real    Exa_10=dera.df10;
real    Exa_01=dera.df01;
real    Exa_20=dera.df20;
real    Exa_11=dera.df11;
real    Exa_02=dera.df02;
real    Exa_30=dera.df30;
real    Exa_21=dera.df21;
real    Exa_12=dera.df12;
real    Exa_03=dera.df03;
real    Exa_40=dera.df40;
real    Exa_31=dera.df31;
real    Exa_22=dera.df22;
real    Exa_13=dera.df13;
real    Exa_04=dera.df04;

real    Exb_10=derb.df10;
real    Exb_01=derb.df01;
real    Exb_20=derb.df20;
real    Exb_11=derb.df11;
real    Exb_02=derb.df02;
real    Exb_30=derb.df30;
real    Exb_21=derb.df21;
real    Exb_12=derb.df12;
real    Exb_03=derb.df03;
real    Exb_40=derb.df40;
real    Exb_31=derb.df31;
real    Exb_22=derb.df22;
real    Exb_13=derb.df13;
real    Exb_04=derb.df04;

/* Powers of total density */
real    roa_rob=roa+rob;
real    roa_rob_p2=roa_rob*roa_rob;
real    roa_rob_p3=roa_rob_p2*roa_rob;
real    roa_rob_p4=roa_rob_p3*roa_rob;
real    roa_rob_p5=roa_rob_p4*roa_rob;

/* Energy */
real E = (roa*Exa + rob*Exb)/roa_rob;

/* Derivatives */

real  df10000  = (rob*Exa - rob*Exb + roa*roa_rob*Exa_10)/roa_rob_p2;
  ds->df10000 += factor*( E + roa_rob*df10000);

real  df01000  = (-roa*Exa + roa*Exb + rob*roa_rob*Exb_10)/roa_rob_p2;
  ds->df01000 += factor*( E + roa_rob*df01000);

real  df00100  = (roa*Exa_01)/roa_rob;
  ds->df00100 += factor*roa_rob*df00100;

real  df00010  =  (rob*Exb_01)/roa_rob;
  ds->df00010 += factor*roa_rob*df00010;

/* Second order derivatives */

real  df20000  =  (1.0/roa_rob_p3)*(-2.0*rob*Exa + 2.0*rob*Exb + 
                  roa_rob*(2.0*rob*Exa_10 + roa*roa_rob*Exa_20));
  ds->df20000 += factor*(2.0*df10000 + roa_rob*df20000); 

real  df11000  =  (1.0/roa_rob_p3)*((roa-rob)*Exa + (rob - roa)*Exb - 
                  roa_rob*(roa*Exa_10 + rob*Exb_10));
  ds->df11000 += factor*(df10000 + df01000 + roa_rob*df11000);

real  df10100  = (rob*Exa_01 + roa*roa_rob*Exa_11)/roa_rob_p2;
  ds->df10100 += factor*(df00100 + roa_rob*df10100); 

real  df10010  = -rob*Exb_01/roa_rob_p2;
  ds->df10010 += factor*(df00010 + roa_rob*df10010);

real  df02000  =  (1.0/roa_rob_p3)*(2.0*roa*Exa - 2.0*roa*Exb + 
                   roa_rob*(2.0*roa*Exb_10 + rob*roa_rob*Exb_20));
  ds->df02000 += factor*(2.0*df01000 + roa_rob*df02000);

real  df01100  =  -roa*Exa_01/roa_rob_p2;
  ds->df01100 += factor*(df00100 + roa_rob*df01100);

real  df01010  = (roa*Exb_01 + rob*roa_rob*Exb_11)/roa_rob_p2;
  ds->df01010 += factor*(df00010 + roa_rob*df01010);

real  df00200  =  roa*Exa_02/roa_rob;
  ds->df00200 += factor*roa_rob*df00200;

/* ds->df00110=0.0 */

real  df00020  = rob*Exb_02/roa_rob;
  ds->df00020 += factor*roa_rob*df00020;

/* Third order derivatives */

real  df30000  =  (1.0/roa_rob_p4)*(6.0*rob*Exa - 6.0*rob*Exb + 
                   roa_rob*(-6.0*rob*Exa_10 +   roa_rob*(3.0*rob*Exa_20 + 
                   roa*roa_rob*Exa_30)));
  ds->df30000 += factor*(3.0*df20000 + roa_rob*df30000);

real  df21000  =  (1.0/roa_rob_p4)*(-2.0*(roa-2.0*rob)*Exa + 
                   2.0*(roa-2.0*rob)*Exb - 
                   roa_rob*(-2.0*(roa-rob)*Exa_10 - 
                   2.0*rob*Exb_10 + roa*roa_rob*Exa_20));
  ds->df21000 += factor*(2*df11000 + df20000 + roa_rob*df21000);

real  df20100  =  (1.0/roa_rob_p3)*(-2.0*rob*Exa_01 + roa_rob*(2.0*rob*Exa_11 + 
                   roa*roa_rob*Exa_21));
  ds->df20100 += factor*(2.0*df10100 + roa_rob*df20100);


real  df20010  =  2.0*rob*Exb_01/roa_rob_p3;
  ds->df20010 += factor*(2.0*df10010 + roa_rob*df20010);

real  df12000  = (1.0/roa_rob_p4)*(-2.0*(2.0*roa-rob)*Exa + 
                  2.0*(2.0*roa-rob)*Exb - 
                  roa_rob*(-2.0*roa*Exa_10 + 
                  2.0*(roa-rob)*Exb_10 + rob*roa_rob*Exb_20));
  ds->df12000 += factor*(2.0*df11000 + df02000 + roa_rob*df12000);


real  df11100  = ((roa-rob)*Exa_01 - roa*roa_rob*Exa_11)/roa_rob_p3;
  ds->df11100 += factor*(df10100 + df01100 + roa_rob*df11100);

real  df11010  = ((rob-roa)*Exb_01 - rob*roa_rob*Exb_11)/roa_rob_p3;
  ds->df11010 += factor*(df10010 + df01010 + roa_rob*df11010);

real  df10200   = (rob*Exa_02 + roa*roa_rob*Exa_12)/roa_rob_p2;
  ds->df10200 += factor*(df00200 + roa_rob*df10200);

/* ds->df10110=0.0 */

real  df10020  =  -rob*Exb_02/roa_rob_p2;
  ds->df10020 += factor*(df00020 + roa_rob*df10020);

real  df03000  = (1.0/roa_rob_p4)*(-6.0*roa*Exa + 6.0*roa*Exb + 
                  roa_rob*(-6.0*roa*Exb_10 + 
                  roa_rob*(3.0*roa*Exb_20 + 
                  rob*roa_rob*Exb_30)));
  ds->df03000 += factor*(3.0*df02000 + roa_rob*df03000);   

real  df02100  =       2.0*roa*Exa_01/roa_rob_p3;
  ds->df02100 += factor*(2.0*df01100 + roa_rob*df02100);

real  df02010  = (1.0/roa_rob_p3)*(-2.0*roa*Exb_01 + roa_rob*(2.0*roa*Exb_11 + 
                  rob*roa_rob*Exb_21));
  ds->df02010 += factor*(2.0*df01010 + roa_rob*df02010);

real  df01200  =  -roa*Exa_02/roa_rob_p2;
  ds->df01200 += factor*(df00200 + roa_rob*df01200);

/* ds->df01110=0.0 */

real  df01020  =  (roa*Exb_02 + rob*roa_rob*Exb_12)/roa_rob_p2;
  ds->df01020 += factor*(df00020 + roa_rob*df01020);


real  df00300  =  roa*Exa_03/roa_rob;
  ds->df00300 += factor*roa_rob*df00300;

/* ds->df00210=0.0 */

/* ds->df00120=0.0 */

real  df00030  = rob*Exb_03/roa_rob;
  ds->df00030 += factor*roa_rob*df00030;

/* Fourth order derivatives */

real  df40000  = (1.0/roa_rob_p5)*(-24.0*rob*Exa + 24.0*rob*Exb + 
                  roa_rob*(24.0*rob*Exa_10 +
                  roa_rob*(-12.0*rob*Exa_20 + 
                  roa_rob*(4.0*rob*Exa_30 +
                  roa*roa_rob*Exa_40))));
  ds->df40000 += factor*(4*df30000 + roa_rob*df40000);

real  df31000  = (1.0/roa_rob_p5)*(6.0*(roa-3.0*rob)*Exa - 6.0*(roa-3.0*rob)*Exb - 
                 roa_rob*(6.0*(roa - 2.0*rob)*Exa_10 +
                 6.0*rob*Exb_10 + 
                 roa_rob*(-3.0*(roa-rob)*Exa_20 +
                 roa*roa_rob*Exa_30)));
   ds->df31000 += factor*(3.0*df21000 + df30000 + roa_rob*df31000);

real  df30100  =  (1.0/roa_rob_p4)*(6.0*rob*Exa_01 + 
                   roa_rob*(-6.0*rob*Exa_11 + 
                   roa_rob*(3.0*rob*Exa_21 +
                   roa*roa_rob*Exa_31)));
  ds->df30100 += factor*(3.0*df20100 + roa_rob*df30100);

real  df30010  =   -6.0*rob*Exb_01/roa_rob_p4;
  ds->df30010 += factor*(3.0*df20010 + roa_rob*df30010);

real  df22000  =    (1.0/roa_rob_p5)*(2.0*(6.0*(roa-rob)*Exa - 6.0*(roa-rob)*Exb + 
                     roa_rob*((-4.0*roa+2.0*rob)*Exa_10 +
                     2.0*(roa-2.0*rob)*Exb_10 + 
                     roa_rob*(roa*Exa_20 + rob*Exb_20))));
  ds->df22000 += factor*(2.0*df12000 + 2.0*df21000 + roa_rob*df22000);

real  df21100  =  -(1.0/roa_rob_p4)*(2.0*(roa-2.0*rob)*Exa_01 + 
                    roa_rob*(-2.0*(roa-rob)*Exa_11 +
                    roa*roa_rob*Exa_21));
  ds->df21100 +=factor*(2.0*df11100 + df20100 + roa_rob*df21100);

real  df21010  = (2.0*((roa-2.0*rob)*Exb_01 + rob*roa_rob*Exb_11))/roa_rob_p4;
  ds->df21010 +=factor*(2.0*df11010 + df20010 + roa_rob*df21010);

real  df20200  =  (1.0/roa_rob_p3)*(-2.0*rob*Exa_02 + 
                   roa_rob*(2.0*rob*Exa_12 + roa*roa_rob*Exa_22));
  ds->df20200 += factor*(2.0*df10200 + roa_rob*df20200);

/* ds->df20110=0.0 */

real  df20020  =  2.0*rob*Exb_02/roa_rob_p3;
  ds->df20020 += factor*(2.0*df10020 + roa_rob*df20020);

real  df13000  = -((6.0*(-3.0*roa+rob)*Exa + 6.0*(3.0*roa-rob)*Exb + 
                  roa_rob*(6.0*roa*Exa_10 + 
                  6.0*(-2.0*roa+rob)*Exb_10 + 
                  roa_rob*(3.0*(roa-rob)*Exb_20 +
                  rob*roa_rob*Exb_30)))/roa_rob_p5);
  ds->df13000 += factor*(3.0*df12000 + df03000 + roa_rob*df13000);

real  df12100  = (2.0*((-2.0*roa+rob)*Exa_01 + roa*roa_rob*Exa_11))/roa_rob_p4;
  ds->df12100 += factor*(2.0*df11100 + df02100 + roa_rob*df12100);

real  df12010  =  (1.0/roa_rob_p4)*((4.0*roa-2.0*rob)*Exb_01 -
                   roa_rob*(2.0*(roa-rob)*Exb_11 + rob*roa_rob*Exb_21));
  ds->df12010 += factor*(2.0*df11010 + df02010 + roa_rob*df12010);

real  df11200  = ((roa-rob)*Exa_02 - roa*roa_rob*Exa_12)/roa_rob_p3;
  ds->df11200 += factor*(df10200 + df01200 + roa_rob*df11200);

/* ds->df11110=0.0 */

real  df11020  = ((-roa+rob)*Exb_02 - rob*roa_rob*Exb_12)/roa_rob_p3;
  ds->df11020 += factor*(df10020 + df01020 + roa_rob*df11020);

real  df10300  = (rob*Exa_03 + roa*roa_rob*Exa_13)/roa_rob_p2;
  ds->df10300 += factor*(df00300 + roa_rob*df10300);

/* ds->df10210=0.0 */

/* ds->df10120=0.0 */

real  df10030  = -rob*Exb_03/roa_rob_p2;
  ds->df10030 += factor*(df00030 + roa_rob*df10030);

real  df04000  = (1.0/roa_rob_p5)*(24.0*roa*Exa - 24.0*roa*Exb + 
                  roa_rob*(24.0*roa*Exb_10 + 
                  roa_rob*(-12.0*roa*Exb_20 + 
                  roa_rob*(4.0*roa*Exb_30 +
                  rob*roa_rob*Exb_40))));
  ds->df04000 += factor*(4.0*df03000 + roa_rob*df04000);

real  df03100  = -6.0*roa*Exa_01/roa_rob_p4;
  ds->df03100 += factor*(3.0*df02100 + roa_rob*df03100);

real  df03010  = (1.0/roa_rob_p4)*(6.0*roa*Exb_01 + 
                  roa_rob*(-6.0*roa*Exb_11 + 
                  roa_rob*(3.0*roa*Exb_21 +
                  rob*roa_rob*Exb_31)));
  ds->df03010 += factor*(3.0*df02010 + roa_rob*df03010);

real  df02200  = 2.0*roa*Exa_02/roa_rob_p3;
  ds->df02200 += factor*(2.0*df01200 + roa_rob*df02200);

/* ds->df02110=0.0 */

real  df02020  = (1.0/roa_rob_p3)*(-2.0*roa*Exb_02 +
                  roa_rob*(2.0*roa*Exb_12 + rob*roa_rob*Exb_22));
  ds->df02020 += factor*(2.0*df01020 + roa_rob*df02020);

real  df01300  = -roa*Exa_03/roa_rob_p2;
  ds->df01300 += factor*(df00300 + roa_rob*df01300);

/*
  ds->df01210=0.0
  ds->df01120=0.0
*/

real  df01030  = (roa*Exb_03 + rob*roa_rob*Exb_13)/roa_rob_p2;
  ds->df01030 += factor*(df00030 + roa_rob*df01030);

real  df00400  = (roa*Exa_04)/roa_rob;
  ds->df00400 += factor*roa_rob*df00400;

/*
  ds->df00310=0.0
  ds->df00220=0.0
  ds->df00130=0.0
*/

real  df00040  = (rob*Exb_04)/roa_rob;
  ds->df00040 += factor*roa_rob*df00040;
}
#endif

static void
pbe_cd1(FirstFuncDrv *ds, real factor, const DftDensProp *dp) {

/* roa, rob and powers */
real  roa = dp->rhoa;
real  rob = dp->rhob;
real  roa_rob=roa+rob;
real  roa_rob_p1f6=pow(roa_rob,(1.0/6.0));

real  roa_rob_p2=roa_rob*roa_rob;

real  roa_rob_p5f6 = roa_rob/roa_rob_p1f6;


/* groa, grob gradients and powers */
real  gro = dp->grada + dp->gradb;

/* Auxaliary functions */
real  ZETA=(roa-rob)/(roa_rob);
real  KS=COR_nKS*roa_rob_p1f6;

real  KS_p2 = KS*KS;
real  KS_p3 = KS_p2*KS;
real  KS_p4 = KS_p3*KS;
real  KS_p7 = KS_p4*KS_p3;


real  MZETA=pow((1.0-ZETA),(2.0/3.0));
real  PZETA=pow((1.0+ZETA),(2.0/3.0));
real  FI=(PZETA + MZETA)/2.0;
real  T=gro/(COR_nT*KS_p7*FI);


/* Powers of Fi */
real  FI_p2 = FI*FI;
real  FI_p3 = FI_p2*FI;

/* Powers of T */
real  T_p2=T*T;
real  T_p4=T_p2*T_p2;

/* The pw92 EC(electron gas corelation) energy */
real  EC = PW92peFunctional.func(dp);

/* H0 and its powers */
real  H0=COR_gamma*FI_p3;
real  H0_p2 = H0*H0;

/* Calc of A */
real  AEXP=exp(-EC/H0);
real  ADENOM=AEXP -1.0;
real  ADENOM_p2=ADENOM*ADENOM;
real  A=COR_delta*(1.0/ADENOM);

/* Powers of A */
real  A_p2=A*A;

/* The pbe Ec (pbe correlation contribution) */
/* H is gga contribution to correlation energy */ 
real  H1=1.0 + COR_delta*T_p2*((1.0 + A*T_p2)/(1.0 + A*T_p2 + A_p2*T_p4));

real  H=H0*log( H1 );
real  E=EC+H;


/* Now start with derivatives *
 * first, common terms used in equations */

/* roafro, robfro and powers */
real  roafro=roa/roa_rob;
real  robfro=rob/roa_rob;

real  roafro_p1f3=pow(roafro,(1.0/3.0));
real  robfro_p1f3=pow(robfro,(1.0/3.0));

real  roafro_p2f3=roafro_p1f3*roafro_p1f3;
real  robfro_p2f3=robfro_p1f3*robfro_p1f3;

real  roafro_p5f3=roafro_p2f3*roafro;
real  robfro_p5f3=robfro_p2f3*robfro;


/* Derivatives of PZETA, MZETA */

real  MZETA_10=(-2.0/3.0)*COR_nZ*robfro_p5f3/rob;
real  PZETA_10=(2.0/3.0)*COR_nZ*rob/(roafro_p1f3*roa_rob_p2);

real  MZETA_01=(2.0/3.0)*COR_nZ*roa/(robfro_p1f3*roa_rob_p2);
real  PZETA_01=(-2.0/3.0)*COR_nZ*roafro_p5f3/roa;

/* FI  derivatives  and powers */
real  FI_10 = 0.5*(MZETA_10 + PZETA_10);
real  FI_01 = 0.5*(MZETA_01 + PZETA_01);

/* Derivatives of KS (with respect to rho) */

real  KS_1 = (1.0/6.0)*COR_nKS/roa_rob_p5f6;

/* Derivatives of EC */

FirstFuncDrv EC_drv;
drv1_clear(&EC_drv);
PW92peFunctional.first(&EC_drv, 1.0, dp);

real  EC_10=EC_drv.df1000;
real  EC_01=EC_drv.df0100;

/* Derivatives of H0 */
real  H0_10 = 3.0*COR_gamma*FI_p2*FI_10;

real  H0_01 = 3.0*COR_gamma*FI_p2*FI_01;

/* Derivatives of T */

/* Terms in T derivatives */
real  FIKS = FI*KS;
real  TDENOM =  FI*KS_p7*COR_nT;
real  TDENOM1 = TDENOM*FIKS;

/* Now T derivatives */

real  T_1000 = -((gro*(FI_10*KS + 7.0*FI*KS_1))/TDENOM1);

real  T_0100 = -((gro*(FI_01*KS + 7.0*FI*KS_1))/TDENOM1);

real  T_0010=1.0/(TDENOM);

real  T_0001=T_0010;

/* Derivatives of -EC/H0 */

real  p10 = (-(EC_10*H0) + EC*H0_10)/H0_p2;

real  p01 = (-(EC_01*H0) + EC*H0_01)/H0_p2;

/* Derivatives of AEXP */

real  AEXP_10 = AEXP*p10;

real  AEXP_01 = AEXP*p01;


/* Derivatives of A */

real  A_10 = -((AEXP_10*COR_delta)/ADENOM_p2);

real  A_01 = -((AEXP_01*COR_delta)/ADENOM_p2);


/* Derivatives of H with respect to H0, T, A */
/* Derivatives with respect to H0 higher then first order vanish */
/* Terms used in the derivatives */
real  At2 = A*T_p2;
real  At2_p2=At2*At2;

real  At_nom = 1.0 + At2;
real  At_denom = At_nom + At2_p2;
real  At_denom_p2 = At_denom*At_denom;

real  At_term0 = 1.0 + 2.0*At2;

real  HDENOM = At_denom_p2*H1;


/* First order derivatives */
real  H_A = ((At_denom - At_nom*At_term0)*COR_delta*H0*T_p4)/HDENOM;

real  H_T = (2.0*(At_denom - At2*At_nom)*At_term0*COR_delta*H0*T)/HDENOM;

real  H_H = log(H1);

/* Final expressions (general derivatives of F[H[roa,rob], T[roa,rob,groa,grob], A[roa,rob]]) */

/* First order derivatives */
real      df10000  = A_10*H_A+H0_10*H_H+H_T*T_1000 + EC_10;
ds->df1000 += factor*(E + roa_rob*df10000);

real      df01000  = A_01*H_A+H0_01*H_H+H_T*T_0100 + EC_01;
ds->df0100 += factor*(E + roa_rob*df01000);

real      df00100  = H_T*T_0010;
ds->df0010 += factor*roa_rob*df00100;

real      df00010  = H_T*T_0001;
ds->df0001 += factor*roa_rob*df00010;

}

static void
pbe_cd2(SecondFuncDrv *ds, real factor, const DftDensProp *dp) {

/* roa, rob and powers */
real  roa = dp->rhoa;
real  rob = dp->rhob;

real  roa_rob=roa+rob;
real  roa_rob_p1f6=pow(roa_rob,(1.0/6.0));

real  roa_rob_p2=roa_rob*roa_rob;
real  roa_rob_p3=roa_rob_p2*roa_rob;
real  roa_rob_p4=roa_rob_p2*roa_rob_p2;

real  roa_rob_p5f6 = roa_rob/roa_rob_p1f6;
real  roa_rob_p11f6 = roa_rob*roa_rob_p5f6;


/* groa, grob gradients and powers */
real  gro = dp->grada + dp->gradb;

/* Auxaliary functions */
real  ZETA=(roa-rob)/(roa_rob);
real  KS=COR_nKS*roa_rob_p1f6;

real  KS_p2 = KS*KS;
real  KS_p3 = KS_p2*KS;
real  KS_p4 = KS_p3*KS;
real  KS_p7 = KS_p4*KS_p3;


real  MZETA=pow((1.0-ZETA),(2.0/3.0));
real  PZETA=pow((1.0+ZETA),(2.0/3.0));
real  FI=(PZETA + MZETA)/2.0;
real  T=gro/(COR_nT*KS_p7*FI);


/* Powers of Fi */
real  FI_p2 = FI*FI;
real  FI_p3 = FI_p2*FI;

/* Powers of T */
real  T_p2=T*T;
real  T_p3=T_p2*T;
real  T_p4=T_p2*T_p2;
real  T_p5=T_p4*T;
real  T_p6=T_p4*T_p2;

/* The pw92 EC(electron gas corelation) energy */
real  EC = PW92peFunctional.func(dp);

/* H0 and its powers */
real  H0=COR_gamma*FI_p3;
real  H0_p2 = H0*H0;
real  H0_p3 = H0_p2*H0;

/* Calc of A */
real  AEXP=exp(-EC/H0);
real  ADENOM=AEXP -1.0;
real  ADENOM_p2=ADENOM*ADENOM;
real  ADENOM_p3=ADENOM_p2*ADENOM;
real  A=COR_delta*(1.0/ADENOM);

/* Powers of A */
real  A_p2=A*A;
real  A_p3=A_p2*A;

/* The pbe Ec (pbe correlation contribution) */
/* H is gga contribution to correlation energy */ 
real  H1=1.0 + COR_delta*T_p2*((1.0 + A*T_p2)/(1.0 + A*T_p2 + A_p2*T_p4));

real  H=H0*log( H1 );
real  E=EC+H;


/* Now start with derivatives *
 * first, common terms used in equations */

/* roafro, robfro and powers */
real  roafro=roa/roa_rob;
real  robfro=rob/roa_rob;

real  roafro_p1f3=pow(roafro,(1.0/3.0));
real  robfro_p1f3=pow(robfro,(1.0/3.0));

real  roafro_p2f3=roafro_p1f3*roafro_p1f3;
real  robfro_p2f3=robfro_p1f3*robfro_p1f3;

real  roafro_p4f3=roafro_p1f3*roafro;
real  robfro_p4f3=robfro_p1f3*robfro;

real  roafro_p5f3=roafro_p2f3*roafro;
real  robfro_p5f3=robfro_p2f3*robfro;


/* Derivatives of PZETA, MZETA */

real  MZETA_10=(-2.0/3.0)*COR_nZ*robfro_p5f3/rob;
real  PZETA_10=(2.0/3.0)*COR_nZ*rob/(roafro_p1f3*roa_rob_p2);

real  MZETA_01=(2.0/3.0)*COR_nZ*roa/(robfro_p1f3*roa_rob_p2);
real  PZETA_01=(-2.0/3.0)*COR_nZ*roafro_p5f3/roa;

real  MZETA_20=(10.0/9.0)*COR_nZ*robfro_p2f3/roa_rob_p2;
real  PZETA_20=(-2.0/9.0)*COR_nZ*rob*(6.0*roa + rob)/(roafro_p4f3*roa_rob_p4);

real  MZETA_11=(2.0/9.0)*COR_nZ*(-2.0*roa+3.0*rob)/(robfro_p1f3*roa_rob_p3);
real  PZETA_11=(2.0/9.0)*COR_nZ*(3.0*roa-2.0*rob)/(roafro_p1f3*roa_rob_p3);

real  MZETA_02=(-2.0/9.0)*COR_nZ*roa*(6.0*rob + roa)/(robfro_p4f3*roa_rob_p4);
real  PZETA_02=(10.0/9.0)*COR_nZ*roafro_p2f3/roa_rob_p2;

/* FI  derivatives  and powers */
real  FI_10 = 0.5*(MZETA_10 + PZETA_10);
real  FI_01 = 0.5*(MZETA_01 + PZETA_01);

real  FI_20 = 0.5*(MZETA_20 + PZETA_20);
real  FI_02 = 0.5*(MZETA_02 + PZETA_02);
real  FI_11 = 0.5*(MZETA_11 + PZETA_11);

real  FI_10_p2 = FI_10*FI_10;
real  FI_01_p2 = FI_01*FI_01;

/* Derivatives of KS (with respect to rho) */

real  KS_1 = (1.0/6.0)*COR_nKS/roa_rob_p5f6;
real  KS_2 = (-5.0/36.0)*COR_nKS/roa_rob_p11f6;

/* Powers of KS */
real  KS_1_p2 = KS_1*KS_1;


/* Derivatives of EC */

SecondFuncDrv EC_drv;
drv2_clear(&EC_drv);
PW92peFunctional.second(&EC_drv, 1.0, dp);

real  EC_10=EC_drv.df1000;
real  EC_01=EC_drv.df0100;
real  EC_20=EC_drv.df2000;
real  EC_02=EC_drv.df0200;
real  EC_11=EC_drv.df1100;


/* Derivatives of H0 */
real  H0_10 = 3.0*COR_gamma*FI_p2*FI_10;

real  H0_01 = 3.0*COR_gamma*FI_p2*FI_01;

real  H0_20 = 3.0*COR_gamma*FI*(2.0*FI_10_p2 + FI*FI_20);

real  H0_11 = 3.0*COR_gamma*FI*(2.0*FI_01*FI_10 + FI*FI_11);

real  H0_02 = 3.0*COR_gamma*FI*(2.0*FI_01_p2 + FI*FI_02);

/* Powers of H0 derivatives */
real  H0_10_p2=H0_10*H0_10;
real  H0_01_p2=H0_01*H0_01;



/* Derivatives of T */

/* Terms in T derivatives */
real  FIKS_2 = FI*KS_2;

real  FIKS = FI*KS;
real  TDENOM =  FI*KS_p7*COR_nT;
real  TDENOM1 = TDENOM*FIKS;
real  TDENOM2 = TDENOM1*FIKS;

/* Now T derivatives */

real  T_1000 = -((gro*(FI_10*KS + 7.0*FI*KS_1))/TDENOM1);

real  T_0100 = -((gro*(FI_01*KS + 7.0*FI*KS_1))/TDENOM1);

real  T_0010=1.0/(TDENOM);

real  T_0001=T_0010;


real  T_2000 = (gro*(2.0*FI_10_p2*KS_p2 + 56.0*FI_p2*KS_1_p2 - 
       FIKS*(7.0*FIKS_2 + FI_20*KS - 14.0*FI_10*KS_1)))/(TDENOM2);


real  T_1100 = (gro*(2.0*FI_01*FI_10*KS_p2 + 56.0*FI_p2*KS_1_p2 + 
       FI*KS*(-7.0*FIKS_2 - FI_11*KS + 7.0*(FI_01 + FI_10)*KS_1)))/(TDENOM2);

real  T_1010 = -((FI_10*KS + 7.0*FI*KS_1)/TDENOM1);


real  T_1001 = T_1010;

real  T_0200 = (gro*(2.0*FI_01_p2*KS_p2 + 56.0*FI_p2*KS_1_p2 - 
       FIKS*(7.0*FIKS_2 + FI_02*KS - 14.0*FI_01*KS_1)))/(TDENOM2);


real  T_0110 = -((FI_01*KS + 7.0*FI*KS_1)/TDENOM1);

real  T_0101 = T_0110;


/* Powers of T derivatives */
real  T_1000_p2 = T_1000*T_1000;

real  T_0100_p2 = T_0100*T_0100;

real  T_0010_p2 = T_0010*T_0010;

real  T_0001_p2 = T_0001*T_0001;


/* Derivatives of -EC/H0 */

real  p10 = (-(EC_10*H0) + EC*H0_10)/H0_p2;

real  p01 = (-(EC_01*H0) + EC*H0_01)/H0_p2;

real  p20 = (-(EC_20*H0_p2) + 2.0*EC_10*H0*H0_10 - 2.0*EC*H0_10_p2 + 
      EC*H0*H0_20)/H0_p3;

real  p11 = (-(EC_11*H0_p2) + EC_10*H0*H0_01 + EC_01*H0*H0_10 - 
      2.0*EC*H0_01*H0_10 + EC*H0*H0_11)/H0_p3;

real  p02 = (-(EC_02*H0_p2) + 2.0*EC_01*H0*H0_01 - 2.0*EC*H0_01_p2 + 
      EC*H0*H0_02)/H0_p3;


/* Powers of p=-EC/H0 */

real  p10_p2 = p10*p10;

real  p01_p2 = p01*p01;




/* Derivatives of AEXP */

real  AEXP_10 = AEXP*p10;

real  AEXP_01 = AEXP*p01;

real  AEXP_20 = AEXP*(p10_p2 + p20);

real  AEXP_11 = AEXP*(p10*p01 + p11);

real  AEXP_02 = AEXP*(p01_p2 + p02);


real  AEXP_10_p2 = AEXP_10*AEXP_10;

real  AEXP_01_p2 = AEXP_01*AEXP_01;





/* Derivatives of A */

real  A_10 = -((AEXP_10*COR_delta)/ADENOM_p2);

real  A_01 = -((AEXP_01*COR_delta)/ADENOM_p2);

real  A_20 = ((2.0*AEXP_10_p2 - ADENOM*AEXP_20)*COR_delta)/ADENOM_p3;

real  A_11 = ((2.0*AEXP_10*AEXP_01 - ADENOM*AEXP_11)*COR_delta)/ADENOM_p3;

real  A_02 = ((2.0*AEXP_01_p2 - ADENOM*AEXP_02)*COR_delta)/ADENOM_p3;




/* Derivatives of H with respect to H0, T, A */
/* Derivatives with respect to H0 higher then first order vanish */
/* Terms used in the derivatives */
real  At2 = A*T_p2;
real  At2_p2=At2*At2;

real  At_nom = 1.0 + At2;
real  At_denom = At_nom + At2_p2;

real  At_term0 = 1.0 + 2.0*At2;
real  At_term1 = T + 3.0*A*T_p3 + 2.0*A_p2*T_p5;
real  At_term2 = A + 3.0*At2*A + 2.0*A_p3*T_p4;
real  At_term3 = T + 2.0*A*T_p3;


/* Powers of terms */
real  At_denom_p2 = At_denom*At_denom;
real  At_denom_p3 = At_denom_p2*At_denom;


real  HDENOM = At_denom_p2*H1;

real  HDENOM_p2 = HDENOM*HDENOM;

real  At_term0_p2 = At_term0*At_term0;

real  At_term1_p2 = At_term1*At_term1;

real  At_term2_p2 = At_term2*At_term2;

real  At_term3_p2 = At_term3*At_term3;





/* First order derivatives */
real  H_A = ((At_denom - At_nom*At_term0)*COR_delta*H0*T_p4)/HDENOM;

real  H_T = (2.0*(At_denom - At2*At_nom)*At_term0*COR_delta*H0*T)/HDENOM;

real  H_H = log(H1);

/* Second order derivatives */
real  H_AA = (COR_delta*H0*T_p6*(-(At_term1_p2*COR_delta) + 
       2.0*At_denom*At_nom*At_term0*(H1 + (2.0*H1*A + COR_delta)*T_p2) - 
       At_denom_p2*(4.0*H1 + (6.0*H1*A + COR_delta)*T_p2)))/HDENOM_p2;

real  H_TT = (2.0*COR_delta*H0*(At_denom_p3*H1 + 
       At2*At_denom_p2*(-5.0 + 6.0*At_denom)*H1 + 
       4.0*At_denom*At_nom*At_term0_p2*A*(H1*A + COR_delta)*T_p4 - 
       2.0*At_term2_p2*COR_delta*T_p6 - 
       At_denom_p2*((23.0 + 22.0*At2)*H1*A_p2*T_p4 + 
       2.0*COR_delta*At_term3_p2)))/HDENOM_p2;

real  H_HT = (2.0*(At_denom - At2*At_nom)*At_term0*COR_delta*T)/HDENOM;


real  H_TA = (2.0*(At_denom - At2*At_nom)*COR_delta*H0*T_p3*(2.0*At_denom_p2*H1 - 
       At_denom*At_term0*(2.0*H1 + (4.0*H1*A + COR_delta)*T_p2) + 
       At_nom*COR_delta*At_term3_p2))/HDENOM_p2;

real  H_HA = ((At_denom - At_nom*At_term0)*COR_delta*T_p4)/HDENOM;


/* Final expressions (general derivatives of F[H[roa,rob], T[roa,rob,groa,grob], A[roa,rob]]) */

/* First order derivatives */
real      df10000  = A_10*H_A+H0_10*H_H+H_T*T_1000 + EC_10;
ds->df1000 += factor*(E + roa_rob*df10000);

real      df01000  = A_01*H_A+H0_01*H_H+H_T*T_0100 + EC_01;
ds->df0100 += factor*(E + roa_rob*df01000);

real      df00100  = H_T*T_0010;
ds->df0010 += factor*roa_rob*df00100;

real      df00010  = H_T*T_0001;
ds->df0001 += factor*roa_rob*df00010;

/* Second order derivatives */
real      df20000  = A_20*H_A+H0_20*H_H+A_10*(A_10*H_AA+2.0*H0_10*H_HA)+
               2.0*(A_10*H_TA+H0_10*H_HT)*T_1000+H_TT*T_1000_p2+H_T*T_2000 + EC_20;
ds->df2000 += factor*(2.0*df10000 + roa_rob*df20000);

real      df11000  = A_11*H_A+H0_11*H_H+A_01*(A_10*H_AA+H0_10*H_HA+H_TA*T_1000)+T_0100*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+H0_01*(A_10*H_HA+H_HT*T_1000)+H_T*T_1100 + EC_11;
ds->df1100 += factor*(df10000 + df01000 + roa_rob*df11000);

real      df10100  = T_0010*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+H_T*T_1010;
ds->df1010  += factor*(df00100 + roa_rob*df10100);

real      df10010  = T_0001*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+H_T*T_1001;
ds->df1001  += factor*(df00010 + roa_rob*df10010);

real      df02000  = A_02*H_A+H0_02*H_H+A_01*(A_01*H_AA+2.0*H0_01*H_HA)+2.0*(A_01*H_TA+H0_01*H_HT)*T_0100+H_TT*T_0100_p2+H_T*T_0200 + EC_02;
ds->df0200 += factor*(2.0*df01000 + roa_rob*df02000);

real      df01100  = T_0010*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)+H_T*T_0110;
ds->df0110 += factor*(df00100 + roa_rob*df01100);

real      df01010  = T_0001*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)+H_T*T_0101;
ds->df0101 += factor*(df00010 + roa_rob*df01010);

real      df00200  = H_TT*T_0010_p2;
ds->df0020 += factor*roa_rob*df00200;

real      df00110  = H_TT*T_0001*T_0010;
ds->df0011 += factor*roa_rob*df00110;
    

real      df00020  = H_TT*T_0001_p2;
ds->df0002 += factor*roa_rob*df00020;

}

static void
pbe_cd3(ThirdFuncDrv *ds, real factor, const DftDensProp *dp) {

/* Constants */
const real  COR_delta_p2 = COR_delta*COR_delta;

/* roa, rob and powers */
real  roa = dp->rhoa;
real  rob = dp->rhob;
real  roa_rob=roa+rob;
real  roa_rob_p1f6=pow(roa_rob,(1.0/6.0));

real  roa_rob_p2=roa_rob*roa_rob;
real  roa_rob_p3=roa_rob_p2*roa_rob;
real  roa_rob_p4=roa_rob_p2*roa_rob_p2;
real  roa_rob_p5=roa_rob_p4*roa_rob;
real  roa_rob_p6=roa_rob_p4*roa_rob_p2;

real  roa_rob_p5f6 = roa_rob/roa_rob_p1f6;
real  roa_rob_p11f6 = roa_rob*roa_rob_p5f6;
real  roa_rob_p17f6 = roa_rob_p2*roa_rob_p5f6;

real  roa_p2=roa*roa;
real  rob_p2=rob*rob;

/* groa, grob gradients and powers */
real  gro = dp->grada + dp->gradb;

/* Auxaliary functions */
real  ZETA=(roa-rob)/(roa_rob);
real  KS=COR_nKS*roa_rob_p1f6;

real  KS_p2 = KS*KS;
real  KS_p3 = KS_p2*KS;
real  KS_p4 = KS_p3*KS;
real  KS_p7 = KS_p4*KS_p3;


real  MZETA=pow((1.0-ZETA),(2.0/3.0));
real  PZETA=pow((1.0+ZETA),(2.0/3.0));
real  FI=(PZETA + MZETA)/2.0;
real  T=gro/(COR_nT*KS_p7*FI);


/* Powers of Fi */
real  FI_p2 = FI*FI;
real  FI_p3 = FI_p2*FI;

/* Powers of T */
real  T_p2=T*T;
real  T_p3=T_p2*T;
real  T_p4=T_p2*T_p2;
real  T_p5=T_p4*T;
real  T_p6=T_p4*T_p2;
real  T_p8=T_p6*T_p2;

/* The pw92 EC(electron gas corelation) energy */
real  EC = PW92peFunctional.func(dp);

/* H0 and its powers */
real  H0=COR_gamma*FI_p3;
real  H0_p2 = H0*H0;
real  H0_p3 = H0_p2*H0;
real  H0_p4 = H0_p3*H0;

/* Calc of A */
real  AEXP=exp(-EC/H0);
real  ADENOM=AEXP -1.0;
real  ADENOM_p2=ADENOM*ADENOM;
real  ADENOM_p3=ADENOM_p2*ADENOM;
real  ADENOM_p4=ADENOM_p3*ADENOM;
real  A=COR_delta*(1.0/ADENOM);

/* Powers of A */
real  A_p2=A*A;
real  A_p3=A_p2*A;

/* The pbe Ec (pbe correlation contribution) */
/* H is gga contribution to correlation energy */ 
real  H1=1.0 + COR_delta*T_p2*((1.0 + A*T_p2)/(1.0 + A*T_p2 + A_p2*T_p4));

real  H=H0*log( H1 );
real  E=EC+H;


/* Now start with derivatives *
 * first, common terms used in equations */

/* roafro, robfro and powers */
real  roafro=roa/roa_rob;
real  robfro=rob/roa_rob;

real  roafro_p1f3=pow(roafro,(1.0/3.0));
real  robfro_p1f3=pow(robfro,(1.0/3.0));

real  roafro_p2f3=roafro_p1f3*roafro_p1f3;
real  robfro_p2f3=robfro_p1f3*robfro_p1f3;

real  roafro_p4f3=roafro_p1f3*roafro;
real  robfro_p4f3=robfro_p1f3*robfro;

real  roafro_p5f3=roafro_p2f3*roafro;
real  robfro_p5f3=robfro_p2f3*robfro;

real  roafro_p7f3=roafro_p4f3*roafro;
real  robfro_p7f3=robfro_p4f3*robfro;

/* Derivatives of PZETA, MZETA */

real  MZETA_10=(-2.0/3.0)*COR_nZ*robfro_p5f3/rob;
real  PZETA_10=(2.0/3.0)*COR_nZ*rob/(roafro_p1f3*roa_rob_p2);

real  MZETA_01=(2.0/3.0)*COR_nZ*roa/(robfro_p1f3*roa_rob_p2);
real  PZETA_01=(-2.0/3.0)*COR_nZ*roafro_p5f3/roa;

real  MZETA_20=(10.0/9.0)*COR_nZ*robfro_p2f3/roa_rob_p2;
real  PZETA_20=(-2.0/9.0)*COR_nZ*rob*(6.0*roa + rob)/(roafro_p4f3*roa_rob_p4);

real  MZETA_11=(2.0/9.0)*COR_nZ*(-2.0*roa+3.0*rob)/(robfro_p1f3*roa_rob_p3);
real  PZETA_11=(2.0/9.0)*COR_nZ*(3.0*roa-2.0*rob)/(roafro_p1f3*roa_rob_p3);

real  MZETA_02=(-2.0/9.0)*COR_nZ*roa*(6.0*rob + roa)/(robfro_p4f3*roa_rob_p4);
real  PZETA_02=(10.0/9.0)*COR_nZ*roafro_p2f3/roa_rob_p2;

real  MZETA_30=(-80.0/27.0)*COR_nZ*robfro_p2f3/roa_rob_p3;
real  PZETA_30=(4.0/27.0)*COR_nZ*rob*(27.0*roa_p2 + 9.0*roa*rob+2.0*rob_p2)/(roafro_p7f3*roa_rob_p6);

real  MZETA_21=(20.0/27.0)*COR_nZ*(roa-3.0*rob)/(robfro_p1f3*roa_rob_p4);
real  PZETA_21=(4.0/27.0)*COR_nZ*(-9.0*roa_p2 + 12.0*roa*rob+ rob_p2)/(roafro_p4f3*roa_rob_p5);

real  MZETA_12=(4.0/27.0)*COR_nZ*(roa_p2 + 12.0*roa*rob -9.0*rob_p2)/(robfro_p4f3*roa_rob_p5);
real  PZETA_12=(20.0/27.0)*COR_nZ*(-3.0*roa+rob)/(roafro_p1f3*roa_rob_p4);

real  MZETA_03=(4.0/27.0)*COR_nZ*roa*(2.0*roa_p2 + 9.0*roa*rob +27.0*rob_p2)/(robfro_p7f3*roa_rob_p6);
real  PZETA_03=(-80.0/27.0)*COR_nZ*roafro_p2f3/roa_rob_p3;

/* FI  derivatives  and powers */
real  FI_10 = 0.5*(MZETA_10 + PZETA_10);
real  FI_01 = 0.5*(MZETA_01 + PZETA_01);

real  FI_20 = 0.5*(MZETA_20 + PZETA_20);
real  FI_02 = 0.5*(MZETA_02 + PZETA_02);
real  FI_11 = 0.5*(MZETA_11 + PZETA_11);

real  FI_30 = 0.5*(MZETA_30 + PZETA_30);
real  FI_03 = 0.5*(MZETA_03 + PZETA_03);
real  FI_21 = 0.5*(MZETA_21 + PZETA_21);
real  FI_12 = 0.5*(MZETA_12 + PZETA_12);

real  FI_10_p2 = FI_10*FI_10;
real  FI_01_p2 = FI_01*FI_01;

real  FI_10_p3 = FI_10_p2*FI_10;
real  FI_01_p3 = FI_01_p2*FI_01;


/* Derivatives of KS (with respect to rho) */

real  KS_1 = (1.0/6.0)*COR_nKS/roa_rob_p5f6;
real  KS_2 = (-5.0/36.0)*COR_nKS/roa_rob_p11f6;
real  KS_3 = (55.0/216.0)*COR_nKS/roa_rob_p17f6;

/* Powers of KS */

real  KS_1_p2 = KS_1*KS_1;
real  KS_1_p3 = KS_1_p2*KS_1;


/* Derivatives of EC */

ThirdFuncDrv EC_drv;
drv3_clear(&EC_drv);
PW92peFunctional.third(&EC_drv, 1.0, dp);

real  EC_10=EC_drv.df1000;
real  EC_01=EC_drv.df0100;
real  EC_20=EC_drv.df2000;
real  EC_02=EC_drv.df0200;
real  EC_11=EC_drv.df1100;
real  EC_30=EC_drv.df3000;
real  EC_03=EC_drv.df0300;
real  EC_21=EC_drv.df2100;
real  EC_12=EC_drv.df1200;


/* Derivatives of H0 */
real  H0_10 = 3.0*COR_gamma*FI_p2*FI_10;

real  H0_01 = 3.0*COR_gamma*FI_p2*FI_01;

real  H0_20 = 3.0*COR_gamma*FI*(2.0*FI_10_p2 + FI*FI_20);

real  H0_11 = 3.0*COR_gamma*FI*(2.0*FI_01*FI_10 + FI*FI_11);

real  H0_02 = 3.0*COR_gamma*FI*(2.0*FI_01_p2 + FI*FI_02);

real  H0_30 = 3.0*COR_gamma*(2.0*FI_10_p3 + 6.0*FI*FI_10*FI_20 + FI_p2*FI_30);

real  H0_21 = 3.0*COR_gamma*(2.0*FI_01*(FI_10_p2 + FI*FI_20) + 
	FI*(4.0*FI_10*FI_11 + FI*FI_21));

real  H0_12 = 3.0*COR_gamma*(2.0*FI_10*(FI_01_p2 + FI*FI_02) + 
        FI*(4.0*FI_01*FI_11 + FI*FI_12));

real  H0_03 = 3.0*COR_gamma*(2.0*FI_01_p3 + 6.0*FI*FI_01*FI_02 + FI_p2*FI_03);

/* Powers of H0 derivatives */
real  H0_10_p2=H0_10*H0_10;
real  H0_01_p2=H0_01*H0_01;
real  H0_10_p3=H0_10_p2*H0_10;
real  H0_01_p3=H0_01_p2*H0_01;



/* Derivatives of T */

/* Terms in T derivatives */
real  FIKS_2 = FI*KS_2;
real  FIKS_3 = FI*KS_3;

real  FIKS = FI*KS;
real  TDENOM =  FI*KS_p7*COR_nT;
real  TDENOM1 = TDENOM*FIKS;
real  TDENOM2 = TDENOM1*FIKS;
real  TDENOM3 = TDENOM2*FIKS;

/* Now T derivatives */

real  T_1000 = -((gro*(FI_10*KS + 7.0*FI*KS_1))/TDENOM1);

real  T_0100 = -((gro*(FI_01*KS + 7.0*FI*KS_1))/TDENOM1);

real  T_0010=1.0/(TDENOM);

real  T_0001=T_0010;


real  T_2000 = (gro*(2.0*FI_10_p2*KS_p2 + 56.0*FI_p2*KS_1_p2 - 
       FIKS*(7.0*FIKS_2 + FI_20*KS - 14.0*FI_10*KS_1)))/(TDENOM2);


real  T_1100 = (gro*(2.0*FI_01*FI_10*KS_p2 + 56.0*FI_p2*KS_1_p2 + 
       FI*KS*(-7.0*FIKS_2 - FI_11*KS + 7.0*(FI_01 + FI_10)*KS_1)))/(TDENOM2);

real  T_1010 = -((FI_10*KS + 7.0*FI*KS_1)/TDENOM1);


real  T_1001 = T_1010;

real  T_0200 = (gro*(2.0*FI_01_p2*KS_p2 + 56.0*FI_p2*KS_1_p2 - 
       FIKS*(7.0*FIKS_2 + FI_02*KS - 14.0*FI_01*KS_1)))/(TDENOM2);


real  T_0110 = -((FI_01*KS + 7.0*FI*KS_1)/TDENOM1);

real  T_0101 = T_0110;

real  T_3000 = -((gro*(6.0*FI_10_p3*KS_p3 + 504.0*FI_p3*KS_1_p3 + 
          3.0*FI*FI_10*KS_p2*(-7.0*FIKS_2 - 2.0*FI_20*KS + 
          14.0*FI_10*KS_1) + FI_p2*KS*(7.0*FIKS_3*KS + 
          FI_30*KS_p2 + 21.0*KS_1*(-8.0*FIKS_2 - FI_20*KS + 
          8.0*FI_10*KS_1))))/(TDENOM3));

real  T_2100 = -((gro*(6.0*FI_01*FI_10_p2*KS_p3 + 504.0*FI_p3*KS_1_p3 + 
          FI*KS_p2*(-7.0*FI_01*FIKS_2 - 4.0*FI_10*FI_11*KS - 
          2.0*FI_01*FI_20*KS + 14.0*FI_10*(2.0*FI_01 + FI_10)*KS_1) + 
          FI_p2*KS*(7.0*FIKS_3*KS + FI_21*KS_p2 - 
          7.0*KS_1*(24.0*FIKS_2 + 2.0*FI_11*KS + FI_20*KS - 
          8.0*(FI_01 + 2.0*FI_10)*KS_1) - 14.0*FI_10*KS*KS_2)))/(TDENOM3));

real  T_2010 = (2.0*FI_10_p2*KS_p2 + 56.0*FI_p2*KS_1_p2 - 
         FI*KS*(7.0*FIKS_2 + FI_20*KS - 14.0*FI_10*KS_1))/(TDENOM2);

real  T_2001=T_2010;


real  T_1110 = (2.0*FI_01*FI_10*KS_p2 + 56.0*FI_p2*KS_1_p2 + 
         FI*KS*(-7.0*FIKS_2 - FI_11*KS + 7.0*(FI_01 + FI_10)*KS_1))/(TDENOM2);

real  T_1101=T_1110;

real  T_1200 = -((gro*(6.0*FI_01_p2*FI_10*KS_p3 + 504.0*FI_p3*KS_1_p3 + 
         2.0*FI*KS_p2*(-7.0*FI_01*FIKS_2 - FI_02*FI_10*KS - 
         2.0*FI_01*FI_11*KS + 7.0*FI_01*(FI_01 + 2.0*FI_10)*KS_1) + 
         FI_p2*KS*(7.0*FIKS_3*KS + FI_12*KS_p2 - 
         7.0*KS_1*(24.0*FIKS_2 + FI_02*KS + 2.0*FI_11*KS - 
         8.0*(2.0*FI_01 + FI_10)*KS_1) - 7.0*FI_10*KS*KS_2)))/(TDENOM3));


real  T_0210 = (2.0*FI_01_p2*KS_p2 + 56.0*FI_p2*KS_1_p2 - 
         FI*KS*(7.0*FIKS_2 + FI_02*KS - 14.0*FI_01*KS_1))/TDENOM2;

real  T_0201=T_0210;

real  T_0300 = -((gro*(6.0*FI_01_p3*KS_p3 + 504.0*FI_p3*KS_1_p3 + 
         3.0*FI*FI_01*KS_p2*(-7.0*FIKS_2 - 2.0*FI_02*KS + 
         14.0*FI_01*KS_1) + FI_p2*KS*(7.0*FIKS_3*KS + 
         FI_03*KS_p2 + 21.0*KS_1*(-8.0*FIKS_2 - FI_02*KS + 
         8.0*FI_01*KS_1))))/(TDENOM3));

/* Powers of T derivatives */
real  T_1000_p2 = T_1000*T_1000;

real  T_0100_p2 = T_0100*T_0100;

real  T_0010_p2 = T_0010*T_0010;
real  T_0010_p3 = T_0010_p2*T_0010;

real  T_0001_p2 = T_0001*T_0001;
real  T_0001_p3 = T_0001_p2*T_0001;


/* Derivatives of -EC/H0 */

real  p10 = (-(EC_10*H0) + EC*H0_10)/H0_p2;

real  p01 = (-(EC_01*H0) + EC*H0_01)/H0_p2;

real  p20 = (-(EC_20*H0_p2) + 2.0*EC_10*H0*H0_10 - 2.0*EC*H0_10_p2 + 
      EC*H0*H0_20)/H0_p3;

real  p11 = (-(EC_11*H0_p2) + EC_10*H0*H0_01 + EC_01*H0*H0_10 - 
      2.0*EC*H0_01*H0_10 + EC*H0*H0_11)/H0_p3;

real  p02 = (-(EC_02*H0_p2) + 2.0*EC_01*H0*H0_01 - 2.0*EC*H0_01_p2 + 
      EC*H0*H0_02)/H0_p3;

real  p30 = (-(EC_30*H0_p3) + 3.0*EC_20*H0_p2*H0_10 - 
      6.0*EC_10*H0*H0_10_p2 + 6.0*EC*H0_10_p3 + 
      3.0*EC_10*H0_p2*H0_20 - 6.0*EC*H0*H0_10*H0_20 + 
      EC*H0_p2*H0_30)/H0_p4;

real  p21 = (-(EC_21*H0_p3) + EC_20*H0_p2*H0_01 + 
      2.0*EC_11*H0_p2*H0_10 - 4.0*EC_10*H0*H0_01*H0_10 - 
      2.0*EC_01*H0*H0_10_p2 + 6.0*EC*H0_01*H0_10_p2 + 
      2.0*EC_10*H0_p2*H0_11 - 4.0*EC*H0*H0_10*H0_11 + 
      EC_01*H0_p2*H0_20 - 2.0*EC*H0*H0_01*H0_20 + 
      EC*H0_p2*H0_21)/H0_p4;

real  p12 = (-(EC_12*H0_p3) + 2.0*EC_11*H0_p2*H0_01 - 
      2.0*EC_10*H0*H0_01_p2 + EC_10*H0_p2*H0_02 + 
      EC_02*H0_p2*H0_10 - 4.0*EC_01*H0*H0_01*H0_10 + 
      6.0*EC*H0_01_p2*H0_10 - 2.0*EC*H0*H0_02*H0_10 + 
      2.0*EC_01*H0_p2*H0_11 - 4.0*EC*H0*H0_01*H0_11 + 
      EC*H0_p2*H0_12)/H0_p4;

real  p03 = (-(EC_03*H0_p3) + 3.0*EC_02*H0_p2*H0_01 - 
      6.0*EC_01*H0*H0_01_p2 + 6.0*EC*H0_01_p3 + 
      3.0*EC_01*H0_p2*H0_02 - 6.0*EC*H0*H0_01*H0_02 + 
      EC*H0_p2*H0_03)/H0_p4;


/* Powers of p=-EC/H0 */

real  p10_p2 = p10*p10;
real  p10_p3 = p10_p2*p10;

real  p01_p2 = p01*p01;
real  p01_p3 = p01_p2*p01;




/* Derivatives of AEXP */

real  AEXP_10 = AEXP*p10;

real  AEXP_01 = AEXP*p01;

real  AEXP_20 = AEXP*(p10_p2 + p20);

real  AEXP_11 = AEXP*(p10*p01 + p11);

real  AEXP_02 = AEXP*(p01_p2 + p02);

real  AEXP_30 = AEXP*(p10_p3 + 3.0*p10*p20 + p30);

real  AEXP_21 = AEXP*(2.0*p10*p11 + p01*(p10_p2 + p20) + p21);

real  AEXP_12 = AEXP*(2.0*p01*p11 + p10*(p01_p2 + p02) + p12);

real  AEXP_03 = AEXP*(p01_p3 + 3.0*p01*p02 + p03);


real  AEXP_10_p2 = AEXP_10*AEXP_10;
real  AEXP_10_p3 =AEXP_10_p2*AEXP_10;

real  AEXP_01_p2 = AEXP_01*AEXP_01;
real  AEXP_01_p3 =AEXP_01_p2*AEXP_01;





/* Derivatives of A */

real  A_10 = -((AEXP_10*COR_delta)/ADENOM_p2);

real  A_01 = -((AEXP_01*COR_delta)/ADENOM_p2);

real  A_20 = ((2.0*AEXP_10_p2 - ADENOM*AEXP_20)*COR_delta)/ADENOM_p3;

real  A_11 = ((2.0*AEXP_10*AEXP_01 - ADENOM*AEXP_11)*COR_delta)/ADENOM_p3;

real  A_02 = ((2.0*AEXP_01_p2 - ADENOM*AEXP_02)*COR_delta)/ADENOM_p3;

real  A_30 = -(((6.0*AEXP_10_p3 - 6.0*ADENOM*AEXP_10*AEXP_20 + 
        ADENOM_p2*AEXP_30)*COR_delta)/ADENOM_p4);

real  A_21 = ((-6.0*AEXP_01*AEXP_10_p2 + 4.0*ADENOM*AEXP_10*AEXP_11 + 
       2.0*ADENOM*AEXP_01*AEXP_20 - ADENOM_p2*AEXP_21)*COR_delta)/ADENOM_p4;

real  A_12 = ((-6.0*AEXP_10*AEXP_01_p2 + 4.0*ADENOM*AEXP_01*AEXP_11 + 
       2.0*ADENOM*AEXP_10*AEXP_02 - ADENOM_p2*AEXP_12)*COR_delta)/ADENOM_p4;


real  A_03 = -(((6.0*AEXP_01_p3 - 6.0*ADENOM*AEXP_01*AEXP_02 + 
        ADENOM_p2*AEXP_03)*COR_delta)/ADENOM_p4);


/* Powers of A derivatives */
real  A_10_p2 = A_10*A_10;
real  A_10_p3 = A_10_p2*A_10;

real  A_01_p2 = A_01*A_01;
real  A_01_p3 = A_01_p2*A_01;



/* Derivatives of H with respect to H0, T, A */
/* Derivatives with respect to H0 higher then first order vanish */
/* Terms used in the derivatives */
real  At2 = A*T_p2;
real  At2_p2=At2*At2;

real  At_nom = 1.0 + At2;
real  At_denom = At_nom + At2_p2;

real  At_term0 = 1.0 + 2.0*At2;
real  At_term1 = T + 3.0*A*T_p3 + 2.0*A_p2*T_p5;
real  At_term2 = A + 3.0*At2*A + 2.0*A_p3*T_p4;
real  At_term3 = T + 2.0*A*T_p3;
real  At_term4 = COR_delta + 2.0*At2*COR_delta;
real  At_term5 = H1 + 2.0*At2*H1;
real  At_term6 = 1.0 + 3.0*At2 + 2.0*At2_p2;


/* Powers of terms */
real  At_nom_p2 = At_nom*At_nom;

real  At_denom_p2 = At_denom*At_denom;
real  At_denom_p3 = At_denom_p2*At_denom;
real  At_denom_p4 = At_denom_p3*At_denom;


real  HDENOM = At_denom_p2*H1;

real  HDENOM_p2 = HDENOM*HDENOM;
real  HDENOM_p3 = HDENOM_p2*HDENOM;

real  At_term0_p2 = At_term0*At_term0;
real  At_term0_p3 = At_term0_p2*At_term0;

real  At_term1_p2 = At_term1*At_term1;

real  At_term2_p2 = At_term2*At_term2;

real  At_term3_p2 = At_term3*At_term3;

real  At_term4_p2 = At_term4*At_term4;

real  At_term5_p2 = At_term5*At_term5;

real  At_term6_p2 = At_term6*At_term6;
real  At_term6_p3 = At_term6_p2*At_term6;



real  H1_p2 = H1*H1;

/* First order derivatives */
real  H_A = ((At_denom - At_nom*At_term0)*COR_delta*H0*T_p4)/HDENOM;

real  H_T = (2.0*(At_denom - At2*At_nom)*At_term0*COR_delta*H0*T)/HDENOM;

real  H_H = log(H1);

/* Second order derivatives */
real  H_AA = (COR_delta*H0*T_p6*(-(At_term1_p2*COR_delta) + 
       2.0*At_denom*At_nom*At_term0*(H1 + (2.0*H1*A + COR_delta)*T_p2) - 
       At_denom_p2*(4.0*H1 + (6.0*H1*A + COR_delta)*T_p2)))/HDENOM_p2;

real  H_TT = (2.0*COR_delta*H0*(At_denom_p3*H1 + 
       At2*At_denom_p2*(-5.0 + 6.0*At_denom)*H1 + 
       4.0*At_denom*At_nom*At_term0_p2*A*(H1*A + COR_delta)*T_p4 - 
       2.0*At_term2_p2*COR_delta*T_p6 - 
       At_denom_p2*((23.0 + 22.0*At2)*H1*A_p2*T_p4 + 
       2.0*COR_delta*At_term3_p2)))/HDENOM_p2;

real  H_HT = (2.0*(At_denom - At2*At_nom)*At_term0*COR_delta*T)/HDENOM;


real  H_TA = (2.0*(At_denom - At2*At_nom)*COR_delta*H0*T_p3*(2.0*At_denom_p2*H1 - 
       At_denom*At_term0*(2.0*H1 + (4.0*H1*A + COR_delta)*T_p2) + 
       At_nom*COR_delta*At_term3_p2))/HDENOM_p2;

real  H_HA = ((At_denom - At_nom*At_term0)*COR_delta*T_p4)/HDENOM;


/* Third order derivatives */
real  H_TTT = (4.0*(At_denom - At2*At_nom)*COR_delta*H0*T*(3.0*At_denom_p3*H1*(2.0*At_denom*H1*A - 
        (1.0 + 6.0*At2)*At_term0*(2.0*H1*A + COR_delta)) + 
        At_denom_p2*At_term0*(12.0*At_term0_p2*H1_p2*A_p2 + 
        3.0*(5.0 + At2*(23.0 + 22.0*At2))*H1*A*COR_delta + 4.0*At_term4_p2)*T_p2 - 
        4.0*At_denom*At_nom*At_term0_p3*A*COR_delta*(3.0*H1*A + 2.0*COR_delta)*T_p4 + 
        4.0*At_nom_p2*At_term0_p3*A_p2*COR_delta_p2*T_p6))/HDENOM_p3;

real  H_AAA = (-2.0*COR_delta*H0*T_p8*(3.0*At_denom_p4*H1_p2 - 3.0*At_denom*At_term1_p2*COR_delta*(H1 + 
        (2.0*H1*A + COR_delta)*T_p2) + COR_delta_p2*T_p4*At_term6_p3 + 
        3.0*At_denom_p2*At_nom*At_term0*(H1_p2 + H1*(4.0*H1*A + 3.0*COR_delta)*T_p2 + 
        (H1*A + COR_delta)*(4.0*H1*A + COR_delta)*T_p4) - At_denom_p3*(9.0*H1_p2 + 
        6.0*H1*(5.0*H1*A + COR_delta)*T_p2 + (24.0*H1_p2*A_p2 + 9.0*H1*A*COR_delta + 
        COR_delta_p2)*T_p4)))/HDENOM_p3;

real  H_HTT = (2.0*COR_delta*(At_denom_p3*H1 + At2*At_denom_p2*(-5.0 + 6.0*At_denom)*H1 + 
        4.0*At_denom*At_nom*At_term0_p2*A*(H1*A + COR_delta)*T_p4 - 
        2.0*At_term2_p2*COR_delta*T_p6 - At_denom_p2*((23.0 + 22.0*At2)*H1*A_p2*T_p4 + 
        2.0*COR_delta*At_term3_p2)))/HDENOM_p2;

real  H_HAA = (COR_delta*T_p6*(-(At_term1_p2*COR_delta) + 
        2.0*At_denom*At_nom*At_term0*(H1 + (2.0*H1*A + COR_delta)*T_p2) - 
        At_denom_p2*(4.0*H1 + (6.0*H1*A + COR_delta)*T_p2)))/HDENOM_p2; 

real  H_TTA = (2.0*(At_denom - At2*At_nom)*COR_delta*H0*T_p2*(6.0*At_denom_p3*(At_denom - 
        (1.0 + 6.0*At2)*At_term0)*H1_p2 - 
        At_denom_p2*H1*(-12.0*At_term0_p3*H1*A + 
        (9.0 + 22.0*At2)*At_denom*COR_delta - 3.0*(3.0 + At2*(17.0 + 18.0*At2))*At_term0*COR_delta)*T_p2 + 
        4.0*At_denom*At_term0_p2*COR_delta*(At_denom*COR_delta - 
        At_nom*(3.0*H1*A + COR_delta))*T_p4 + 
        4.0*At_term0_p2*A*COR_delta*(At_nom_p2*At_term0*COR_delta - 
        3.0*At_denom*At_nom*(2.0*H1*A + COR_delta))*T_p6))/HDENOM_p3;


real  H_TAA = (4.0*(At_denom - At2*At_nom)*COR_delta*H0*T_p5*(At_nom_p2*At_term0_p3*COR_delta_p2*T_p4 + 
        At_denom*At_nom*At_term3_p2*COR_delta*(-3.0*H1 - 2.0*(3.0*H1*A + COR_delta)*T_p2) - 
        2.0*At_denom_p3*H1*(3.0*H1 + (6.0*H1*A + COR_delta)*T_p2) + 
        At_denom_p2*At_term0*(3.0*At_term5_p2 + 3.0*(2.0 + 3.0*At2)*H1*COR_delta*T_p2 + 
        COR_delta_p2*T_p4)))/HDENOM_p3;

real  H_HTA = (2.0*(At_denom - At2*At_nom)*COR_delta*T_p3*(2.0*At_denom*(At_denom - At_term0)*H1 + 
        At_nom*At_term3_p2*COR_delta - 
        At_denom*At_term0*(4.0*H1*A + COR_delta)*T_p2))/HDENOM_p2;


/* Final expressions (general derivatives of F[H[roa,rob], T[roa,rob,groa,grob], A[roa,rob]]) */

/* First order derivatives */
real      df10000  = A_10*H_A+H0_10*H_H+H_T*T_1000 + EC_10;
ds->df1000 += factor*(E + roa_rob*df10000);

real      df01000  = A_01*H_A+H0_01*H_H+H_T*T_0100 + EC_01;
ds->df0100 += factor*(E + roa_rob*df01000);

real      df00100  = H_T*T_0010;
ds->df0010 += factor*roa_rob*df00100;

real      df00010  = H_T*T_0001;
ds->df0001 += factor*roa_rob*df00010;

/* Second order derivatives */
real      df20000  = A_20*H_A+H0_20*H_H+A_10*(A_10*H_AA+2.0*H0_10*H_HA)+
               2.0*(A_10*H_TA+H0_10*H_HT)*T_1000+H_TT*T_1000_p2+H_T*T_2000 + EC_20;
ds->df2000 += factor*(2.0*df10000 + roa_rob*df20000);

real      df11000  = A_11*H_A+H0_11*H_H+A_01*(A_10*H_AA+H0_10*H_HA+H_TA*T_1000)+T_0100*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+H0_01*(A_10*H_HA+H_HT*T_1000)+H_T*T_1100 + EC_11;
ds->df1100 += factor*(df10000 + df01000 + roa_rob*df11000);

real      df10100  = T_0010*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+H_T*T_1010;
ds->df1010  += factor*(df00100 + roa_rob*df10100);

real      df10010  = T_0001*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+H_T*T_1001;
ds->df1001  += factor*(df00010 + roa_rob*df10010);

real      df02000  = A_02*H_A+H0_02*H_H+A_01*(A_01*H_AA+2.0*H0_01*H_HA)+2.0*(A_01*H_TA+H0_01*H_HT)*T_0100+H_TT*T_0100_p2+H_T*T_0200 + EC_02;
ds->df0200 += factor*(2.0*df01000 + roa_rob*df02000);

real      df01100  = T_0010*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)+H_T*T_0110;
ds->df0110 += factor*(df00100 + roa_rob*df01100);

real      df01010  = T_0001*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)+H_T*T_0101;
ds->df0101 += factor*(df00010 + roa_rob*df01010);

real      df00200  = H_TT*T_0010_p2;
ds->df0020 += factor*roa_rob*df00200;

real      df00110  = H_TT*T_0001*T_0010;
ds->df0011 += factor*roa_rob*df00110;
    

real      df00020  = H_TT*T_0001_p2;
ds->df0002 += factor*roa_rob*df00020;

/* Third order derivatives */
real      df30000  = A_30*H_A+A_10_p3*H_AAA+H0_30*H_H+3.0*A_10_p2*(H0_10*H_HAA+H_TAA*T_1000)+3.0*A_10*(A_20*H_AA+H0_20*H_HA+2.0*H0_10*H_HTA*T_1000+H_TTA*T_1000_p2+H_TA*T_2000)+T_1000*(3.0*A_20*H_TA+3.0*H0_20*H_HT+H_TTT*T_1000_p2+3.0*H_TT*T_2000)+3.0*H0_10*(A_20*H_HA+H_HTT*T_1000_p2+H_HT*T_2000)+H_T*T_3000 + EC_30;
ds->df3000 += factor*(3.0*df20000 + roa_rob*df30000);

real      df21000  = A_21*H_A+H0_21*H_H+2.0*H0_10*A_11*H_HA+H0_01*A_20*H_HA+A_20*H_TA*T_0100+H0_20*H_HT*T_0100+A_10_p2*(A_01*H_AAA+H0_01*H_HAA+H_TAA*T_0100)+2.0*A_11*H_TA*T_1000+2.0*H0_11*H_HT*T_1000+2.0*H0_10*H_HTT*T_0100*T_1000+H0_01*H_HTT*T_1000_p2+H_TTT*T_0100*T_1000_p2+2.0*H0_10*H_HT*T_1100+2.0*H_TT*T_1000*T_1100+2.0*A_10*(A_11*H_AA+H0_11*H_HA+H0_10*(A_01*H_HAA+H_HTA*T_0100)+(A_01*H_TAA+H0_01*H_HTA+H_TTA*T_0100)*T_1000+H_TA*T_1100)+(H0_01*H_HT+H_TT*T_0100)*T_2000+A_01*(A_20*H_AA+H0_20*H_HA+2.0*H0_10*H_HTA*T_1000+H_TTA*T_1000_p2+H_TA*T_2000)+H_T*T_2100 + EC_21;
ds->df2100 += factor*(2*df11000 + df20000 + roa_rob*df21000);

real      df20100  = 2.0*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)*T_1010+T_0010*(A_20*H_TA+H0_20*H_HT+A_10*(A_10*H_TAA+2.0*H0_10*H_HTA)+2.0*(A_10*H_TTA+H0_10*H_HTT)*T_1000+H_TTT*T_1000_p2+H_TT*T_2000)+H_T*T_2010;
ds->df2010 += factor*(2.0*df10100 + roa_rob*df20100);

real      df20010  = 2.0*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)*T_1001+T_0001*(A_20*H_TA+H0_20*H_HT+A_10*(A_10*H_TAA+2.0*H0_10*H_HTA)+2.0*(A_10*H_TTA+H0_10*H_HTT)*T_1000+H_TTT*T_1000_p2+H_TT*T_2000)+H_T*T_2001;
ds->df2001 += factor*(2.0*df10010 + roa_rob*df20010);

real      df12000  = A_12*H_A+H0_12*H_H+A_11*(A_01*H_AA+2.0*H0_01*H_HA)+A_02*(A_10*H_AA+H0_10*H_HA+H_TA*T_1000)+T_0200*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+T_0100_p2*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+H0_02*(A_10*H_HA+H_HT*T_1000)+A_01*(A_11*H_AA+2.0*H0_11*H_HA+A_01*(A_10*H_AAA+H0_10*H_HAA+H_TAA*T_1000)+2.0*H0_01*(A_10*H_HAA+H_HTA*T_1000))+2.0*T_0100*(A_11*H_TA+H0_11*H_HT+A_01*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+H0_01*(A_10*H_HTA+H_HTT*T_1000))+2.0*(A_01*H_TA+H0_01*H_HT)*T_1100+2.0*H_TT*T_0100*T_1100+H_T*T_1200 + EC_12;
ds->df1200 += factor*(2.0*df11000 + df02000 + roa_rob*df12000);

real      df11100  = T_0110*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_1010+T_0010*(A_11*H_TA+H0_11*H_HT+A_01*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+T_0100*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+H0_01*(A_10*H_HTA+H_HTT*T_1000)+H_TT*T_1100)+H_T*T_1110;
ds->df1110 += factor*(df10100 + df01100 + roa_rob*df11100);

real      df11010  = T_0101*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_1001+T_0001*(A_11*H_TA+H0_11*H_HT+A_01*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+T_0100*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+H0_01*(A_10*H_HTA+H_HTT*T_1000)+H_TT*T_1100)+H_T*T_1101;
ds->df1101 += factor*(df10010 + df01010 + roa_rob*df11010);

real      df10200  = T_0010*(T_0010*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+2.0*H_TT*T_1010);
ds->df1020 += factor*(df00200 + roa_rob*df10200);

real      df10110  = T_0010*(T_0001*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+H_TT*T_1001)+H_TT*T_0001*T_1010;
ds->df1011 += factor*(df00110 + roa_rob*df10110);

real      df10020  = T_0001*(T_0001*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+2.0*H_TT*T_1001);
ds->df1002 += factor*(df00020 + roa_rob*df10020);

real      df03000  = A_03*H_A+A_01_p3*H_AAA+H0_03*H_H+3.0*A_01_p2*(H0_01*H_HAA+H_TAA*T_0100)+3.0*A_01*(A_02*H_AA+H0_02*H_HA+2.0*H0_01*H_HTA*T_0100+H_TTA*T_0100_p2+H_TA*T_0200)+T_0100*(3.0*A_02*H_TA+3.0*H0_02*H_HT+H_TTT*T_0100_p2+3.0*H_TT*T_0200)+3.0*H0_01*(A_02*H_HA+H_HTT*T_0100_p2+H_HT*T_0200)+H_T*T_0300 + EC_03;
ds->df0300 += factor*(3.0*df02000 + roa_rob*df03000);

real      df02100  = 2.0*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_0110+T_0010*(A_02*H_TA+H0_02*H_HT+A_01*(A_01*H_TAA+2.0*H0_01*H_HTA)+2.0*(A_01*H_TTA+H0_01*H_HTT)*T_0100+H_TTT*T_0100_p2+H_TT*T_0200)+H_T*T_0210;
ds->df0210 += factor*(2.0*df01100 + roa_rob*df02100);

real      df02010  = 2.0*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_0101+T_0001*(A_02*H_TA+H0_02*H_HT+A_01*(A_01*H_TAA+2.0*H0_01*H_HTA)+2.0*(A_01*H_TTA+H0_01*H_HTT)*T_0100+H_TTT*T_0100_p2+H_TT*T_0200)+H_T*T_0201;
ds->df0201 += factor*(2.0*df01010 + roa_rob*df02010);

real      df01200  = T_0010*(T_0010*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)+2.0*H_TT*T_0110);
ds->df0120 += factor*(df00200 + roa_rob*df01200);

real      df01110  = T_0010*(T_0001*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)+H_TT*T_0101)+H_TT*T_0001*T_0110;
ds->df0111 += factor*(df00110 + roa_rob*df01110);

real      df01020  = T_0001*(T_0001*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)+2.0*H_TT*T_0101);
ds->df0102 += factor*(df00020 + roa_rob*df01020);

real      df00300  = H_TTT*T_0010_p3;
ds->df0030 += factor*roa_rob*df00300;

real      df00210  = H_TTT*T_0001*T_0010_p2;
ds->df0021 += factor*roa_rob*df00210;

real      df00120  = H_TTT*T_0001_p2*T_0010;
ds->df0012 += factor*roa_rob*df00120;

real      df00030  = H_TTT*T_0001_p3;
ds->df0003 += factor*roa_rob*df00030;

}

#ifdef FOURTH_ORDER_DERIVATIVES
static void
pbe_cd4(FourthFuncDrv *ds, real factor, const DftDensProp *dp) {

/* Constants */
real  COR_delta_p2 = COR_delta*COR_delta;
real  COR_delta_p3 = COR_delta_p2*COR_delta;

/* roa, rob and powers */
real  roa = dp->rhoa;
real  rob = dp->rhob;
real  roa_rob=roa+rob;
real  roa_rob_p1f6=pow(roa_rob,(1.0/6.0));

real  roa_rob_p2=roa_rob*roa_rob;
real  roa_rob_p3=roa_rob_p2*roa_rob;
real  roa_rob_p4=roa_rob_p2*roa_rob_p2;
real  roa_rob_p5=roa_rob_p4*roa_rob;
real  roa_rob_p6=roa_rob_p4*roa_rob_p2;
real  roa_rob_p7=roa_rob_p6*roa_rob;

real  roa_rob_p5f6 = roa_rob/roa_rob_p1f6;
real  roa_rob_p11f6 = roa_rob*roa_rob_p5f6;
real  roa_rob_p17f6 = roa_rob_p2*roa_rob_p5f6;
real  roa_rob_p23f6 = roa_rob_p3*roa_rob_p5f6;

real  roa_p2=roa*roa;
real  rob_p2=rob*rob;

real  roa_p3=roa_p2*roa;
real  rob_p3=rob_p2*rob;

real  roa_p4=roa_p3*roa;
real  rob_p4=rob_p3*rob;

/* groa, grob gradients and powers */
real  gro = dp->grada + dp->gradb;

/* Auxaliary functions */
real  ZETA=(roa-rob)/(roa_rob);
real  KS=COR_nKS*roa_rob_p1f6;

real  KS_p2 = KS*KS;
real  KS_p3 = KS_p2*KS;
real  KS_p4 = KS_p3*KS;
real  KS_p7 = KS_p4*KS_p3;


real  MZETA=pow((1.0-ZETA),(2.0/3.0));
real  PZETA=pow((1.0+ZETA),(2.0/3.0));
real  FI=(PZETA + MZETA)/2.0;
real  T=gro/(COR_nT*KS_p7*FI);


/* Powers of Fi */
real  FI_p2 = FI*FI;
real  FI_p3 = FI_p2*FI;
real  FI_p4 = FI_p3*FI;

/* Powers of T */
real  T_p2=T*T;
real  T_p3=T_p2*T;
real  T_p4=T_p2*T_p2;
real  T_p5=T_p4*T;
real  T_p6=T_p4*T_p2;
real  T_p7=T_p6*T;
real  T_p8=T_p6*T_p2;
real  T_p10=T_p8*T_p2;
real  T_p12=T_p10*T_p2;

/* The pw92 EC(electron gas corelation) energy */
real  EC = PW92peFunctional.func(dp);

/* H0 and its powers */
real  H0=COR_gamma*FI_p3;
real  H0_p2 = H0*H0;
real  H0_p3 = H0_p2*H0;
real  H0_p4 = H0_p3*H0;
real  H0_p5 = H0_p4*H0;

/* Calc of A */
real  AEXP=exp(-EC/H0);
real  ADENOM=AEXP -1.0;
real  ADENOM_p2=ADENOM*ADENOM;
real  ADENOM_p3=ADENOM_p2*ADENOM;
real  ADENOM_p4=ADENOM_p3*ADENOM;
real  ADENOM_p5=ADENOM_p4*ADENOM;
real  A=COR_delta*(1.0/ADENOM);

/* Powers of A */
real  A_p2=A*A;
real  A_p3=A_p2*A;

/* The pbe Ec (pbe correlation contribution) */
/* H is gga contribution to correlation energy */ 
real  H1=1.0 + COR_delta*T_p2*((1.0 + A*T_p2)/(1.0 + A*T_p2 + A_p2*T_p4));

real  H=H0*log( H1 );
real  E=EC+H;


/* Now start with derivatives *
 * first, common terms used in equations */

/* roafro, robfro and powers */
real  roafro=roa/roa_rob;
real  robfro=rob/roa_rob;

real  roafro_p1f3=pow(roafro,(1.0/3.0));
real  robfro_p1f3=pow(robfro,(1.0/3.0));

real  roafro_p2f3=roafro_p1f3*roafro_p1f3;
real  robfro_p2f3=robfro_p1f3*robfro_p1f3;

real  roafro_p4f3=roafro_p1f3*roafro;
real  robfro_p4f3=robfro_p1f3*robfro;

real  roafro_p5f3=roafro_p2f3*roafro;
real  robfro_p5f3=robfro_p2f3*robfro;

real  roafro_p7f3=roafro_p4f3*roafro;
real  robfro_p7f3=robfro_p4f3*robfro;

/* Derivatives of PZETA, MZETA */

real  MZETA_10=(-2.0/3.0)*COR_nZ*robfro_p5f3/rob;
real  PZETA_10=(2.0/3.0)*COR_nZ*rob/(roafro_p1f3*roa_rob_p2);

real  MZETA_01=(2.0/3.0)*COR_nZ*roa/(robfro_p1f3*roa_rob_p2);
real  PZETA_01=(-2.0/3.0)*COR_nZ*roafro_p5f3/roa;

real  MZETA_20=(10.0/9.0)*COR_nZ*robfro_p2f3/roa_rob_p2;
real  PZETA_20=(-2.0/9.0)*COR_nZ*rob*(6.0*roa + rob)/(roafro_p4f3*roa_rob_p4);

real  MZETA_11=(2.0/9.0)*COR_nZ*(-2.0*roa+3.0*rob)/(robfro_p1f3*roa_rob_p3);
real  PZETA_11=(2.0/9.0)*COR_nZ*(3.0*roa-2.0*rob)/(roafro_p1f3*roa_rob_p3);

real  MZETA_02=(-2.0/9.0)*COR_nZ*roa*(6.0*rob + roa)/(robfro_p4f3*roa_rob_p4);
real  PZETA_02=(10.0/9.0)*COR_nZ*roafro_p2f3/roa_rob_p2;

real  MZETA_30=(-80.0/27.0)*COR_nZ*robfro_p2f3/roa_rob_p3;
real  PZETA_30=(4.0/27.0)*COR_nZ*rob*(27.0*roa_p2 + 9.0*roa*rob+2.0*rob_p2)/(roafro_p7f3*roa_rob_p6);

real  MZETA_21=(20.0/27.0)*COR_nZ*(roa-3.0*rob)/(robfro_p1f3*roa_rob_p4);
real  PZETA_21=(4.0/27.0)*COR_nZ*(-9.0*roa_p2 + 12.0*roa*rob+ rob_p2)/(roafro_p4f3*roa_rob_p5);

real  MZETA_12=(4.0/27.0)*COR_nZ*(roa_p2 + 12.0*roa*rob -9.0*rob_p2)/(robfro_p4f3*roa_rob_p5);
real  PZETA_12=(20.0/27.0)*COR_nZ*(-3.0*roa+rob)/(roafro_p1f3*roa_rob_p4);

real  MZETA_03=(4.0/27.0)*COR_nZ*roa*(2.0*roa_p2 + 9.0*roa*rob +27.0*rob_p2)/(robfro_p7f3*roa_rob_p6);
real  PZETA_03=(-80.0/27.0)*COR_nZ*roafro_p2f3/roa_rob_p3;

real  MZETA_40=(880.0/81.0)*COR_nZ*robfro_p2f3/roa_rob_p4;
real  PZETA_40=(-8.0/81.0)*COR_nZ*rob*roafro_p2f3*(162.0*roa_p3+81.0*roa_p2*rob + 36.0*roa*rob_p2 + 7.0*rob_p3)/(roa_p4*roa_rob_p4);

real  MZETA_31=(80.0/81.0)*COR_nZ*(-2.0*roa + 9.0*rob)/(robfro_p1f3*roa_rob_p5);
real  PZETA_31=(4.0/81.0)*COR_nZ*(81.0*roa_p3-162.0*roa_p2*rob - 27.0*roa*rob_p2 - 4.0*rob_p3)/(roafro_p7f3*roa_rob_p7);

real  MZETA_22=(-20.0/81.0)*COR_nZ*(roa_p2 + 18.0*roa*rob - 27.0*rob_p2)/(robfro_p4f3*roa_rob_p6);
real  PZETA_22=(20.0/81.0)*COR_nZ*(27.0*roa_p2 - 18.0*roa*rob -rob_p2)/(roafro_p4f3*roa_rob_p6);

real  MZETA_13=(-4.0/81.0)*COR_nZ*(4.0*roa_p3 + 27.0*roa_p2*rob + 162.0*roa*rob_p2 -81.0*rob_p3)/(robfro_p7f3*roa_rob_p7);
real  PZETA_13=(80.0/81.0)*COR_nZ*(9.0*roa - 2.0*rob)/(roafro_p1f3*roa_rob_p5);

real  MZETA_04=(-8.0/81.0)*COR_nZ*roa*robfro_p2f3*(7.0*roa_p3+36.0*roa_p2*rob + 81.0*roa*rob_p2 + 162.0*rob_p3)/(rob_p4*roa_rob_p4);
real  PZETA_04=(880.0/81.0)*COR_nZ*roafro_p2f3/roa_rob_p4;

/* FI  derivatives  and powers */
real  FI_10 = 0.5*(MZETA_10 + PZETA_10);
real  FI_01 = 0.5*(MZETA_01 + PZETA_01);

real  FI_20 = 0.5*(MZETA_20 + PZETA_20);
real  FI_02 = 0.5*(MZETA_02 + PZETA_02);
real  FI_11 = 0.5*(MZETA_11 + PZETA_11);

real  FI_30 = 0.5*(MZETA_30 + PZETA_30);
real  FI_03 = 0.5*(MZETA_03 + PZETA_03);
real  FI_21 = 0.5*(MZETA_21 + PZETA_21);
real  FI_12 = 0.5*(MZETA_12 + PZETA_12);

real  FI_40 = 0.5*(MZETA_40 + PZETA_40);
real  FI_04 = 0.5*(MZETA_04 + PZETA_04);
real  FI_31 = 0.5*(MZETA_31 + PZETA_31);
real  FI_13 = 0.5*(MZETA_13 + PZETA_13);
real  FI_22 = 0.5*(MZETA_22 + PZETA_22);


real  FI_10_p2 = FI_10*FI_10;
real  FI_01_p2 = FI_01*FI_01;

real  FI_10_p3 = FI_10_p2*FI_10;
real  FI_01_p3 = FI_01_p2*FI_01;

real  FI_10_p4 = FI_10_p3*FI_10;
real  FI_01_p4 = FI_01_p3*FI_01;

real  FI_20_p2 = FI_20*FI_20;
real  FI_02_p2 = FI_02*FI_02;

real  FI_11_p2 = FI_11*FI_11;

/* Derivatives of KS (with respect to rho) */

real  KS_1 = (1.0/6.0)*COR_nKS/roa_rob_p5f6;
real  KS_2 = (-5.0/36.0)*COR_nKS/roa_rob_p11f6;
real  KS_3 = (55.0/216.0)*COR_nKS/roa_rob_p17f6;
real  KS_4 = (-935.0/1296.0)*COR_nKS/roa_rob_p23f6;

/* Powers of KS */

real  KS_1_p2 = KS_1*KS_1;
real  KS_1_p3 = KS_1_p2*KS_1;
real  KS_1_p4 = KS_1_p3*KS_1;

real  KS_2_p2 = KS_2*KS_2;

/* Derivatives of EC */
FourthFuncDrv EC_drv;
drv4_clear(&EC_drv);
PW92peFunctional.fourth(&EC_drv, 1.0, dp);

real  EC_10=EC_drv.df10000;
real  EC_01=EC_drv.df01000;
real  EC_20=EC_drv.df20000;
real  EC_02=EC_drv.df02000;
real  EC_11=EC_drv.df11000;
real  EC_30=EC_drv.df30000;
real  EC_03=EC_drv.df03000;
real  EC_21=EC_drv.df21000;
real  EC_12=EC_drv.df12000;
real  EC_40=EC_drv.df40000;
real  EC_04=EC_drv.df04000;
real  EC_31=EC_drv.df31000;
real  EC_13=EC_drv.df13000;
real  EC_22=EC_drv.df22000;


/* Derivatives of H0 */
real  H0_10 = 3.0*COR_gamma*FI_p2*FI_10;

real  H0_01 = 3.0*COR_gamma*FI_p2*FI_01;

real  H0_20 = 3.0*COR_gamma*FI*(2.0*FI_10_p2 + FI*FI_20);

real  H0_11 = 3.0*COR_gamma*FI*(2.0*FI_01*FI_10 + FI*FI_11);

real  H0_02 = 3.0*COR_gamma*FI*(2.0*FI_01_p2 + FI*FI_02);

real  H0_30 = 3.0*COR_gamma*(2.0*FI_10_p3 + 6.0*FI*FI_10*FI_20 + FI_p2*FI_30);

real  H0_21 = 3.0*COR_gamma*(2.0*FI_01*(FI_10_p2 + FI*FI_20) + 
	FI*(4.0*FI_10*FI_11 + FI*FI_21));

real  H0_12 = 3.0*COR_gamma*(2.0*FI_10*(FI_01_p2 + FI*FI_02) + 
        FI*(4.0*FI_01*FI_11 + FI*FI_12));

real  H0_03 = 3.0*COR_gamma*(2.0*FI_01_p3 + 6.0*FI*FI_01*FI_02 + FI_p2*FI_03);

real  H0_40 = 3.0*COR_gamma*(12.0*FI_10_p2*FI_20+ 8.0*FI*FI_10*FI_30 + 
	FI*(6.0*FI_20_p2 + FI*FI_40));

real  H0_31 = 3.0*COR_gamma*(6.0*FI_10_p2*FI_11 + 6.0*FI_10*(FI_01*FI_20 + FI*FI_21) + 
      FI*(6.0*FI_11*FI_20 + 2.0*FI_01*FI_30 + FI*FI_31));

real  H0_22 = 3.0*COR_gamma*(2.0*FI_01_p2*FI_20 + 2.0*FI_02*(FI_10_p2 + FI*FI_20) + 
      4.0*FI_01*(2.0*FI_10*FI_11 + FI*FI_21) + 
      FI*(4.0*FI_11_p2 + 4.0*FI_10*FI_12 + FI*FI_22));

real  H0_13 = 3.0*COR_gamma*(6.0*FI_01_p2*FI_11 + 6.0*FI_01*(FI_02*FI_10 + FI*FI_12) + 
      FI*(2.0*FI_03*FI_10 + 6.0*FI_02*FI_11 + FI*FI_13));

real  H0_04 = 3.0*COR_gamma*(12.0*FI_01_p2*FI_02 + 8.0*FI*FI_01*FI_03 + 
      FI*(6.0*FI_02_p2 + FI*FI_04));

/* Powers of H0 derivatives */
real  H0_10_p2=H0_10*H0_10;
real  H0_01_p2=H0_01*H0_01;
real  H0_10_p3=H0_10_p2*H0_10;
real  H0_01_p3=H0_01_p2*H0_01;
real  H0_10_p4=H0_10_p3*H0_10;
real  H0_01_p4=H0_01_p3*H0_01;

real  H0_20_p2=H0_20*H0_20;
real  H0_02_p2=H0_02*H0_02;

real  H0_11_p2=H0_11*H0_11;

/* Derivatives of T */

/* Terms in T derivatives */
real  FIKS_2 = FI*KS_2;
real  FIKS_3 = FI*KS_3;
real  FIKS_4 = FI*KS_4;

real  FIKS = FI*KS;
real  TDENOM =  FI*KS_p7*COR_nT;
real  TDENOM1 = TDENOM*FIKS;
real  TDENOM2 = TDENOM1*FIKS;
real  TDENOM3 = TDENOM2*FIKS;
real  TDENOM4 = TDENOM3*FIKS;

/* Now T derivatives */

real  T_1000 = -((gro*(FI_10*KS + 7.0*FI*KS_1))/TDENOM1);

real  T_0100 = -((gro*(FI_01*KS + 7.0*FI*KS_1))/TDENOM1);

real  T_0010=1.0/(TDENOM);

real  T_0001=T_0010;


real  T_2000 = (gro*(2.0*FI_10_p2*KS_p2 + 56.0*FI_p2*KS_1_p2 - 
       FIKS*(7.0*FIKS_2 + FI_20*KS - 14.0*FI_10*KS_1)))/(TDENOM2);


real  T_1100 = (gro*(2.0*FI_01*FI_10*KS_p2 + 56.0*FI_p2*KS_1_p2 + 
       FI*KS*(-7.0*FIKS_2 - FI_11*KS + 7.0*(FI_01 + FI_10)*KS_1)))/(TDENOM2);

real  T_1010 = -((FI_10*KS + 7.0*FI*KS_1)/TDENOM1);


real  T_1001 = T_1010;

real  T_0200 = (gro*(2.0*FI_01_p2*KS_p2 + 56.0*FI_p2*KS_1_p2 - 
       FIKS*(7.0*FIKS_2 + FI_02*KS - 14.0*FI_01*KS_1)))/(TDENOM2);


real  T_0110 = -((FI_01*KS + 7.0*FI*KS_1)/TDENOM1);

real  T_0101 = T_0110;

real  T_3000 = -((gro*(6.0*FI_10_p3*KS_p3 + 504.0*FI_p3*KS_1_p3 + 
          3.0*FI*FI_10*KS_p2*(-7.0*FIKS_2 - 2.0*FI_20*KS + 
          14.0*FI_10*KS_1) + FI_p2*KS*(7.0*FIKS_3*KS + 
          FI_30*KS_p2 + 21.0*KS_1*(-8.0*FIKS_2 - FI_20*KS + 
          8.0*FI_10*KS_1))))/(TDENOM3));

real  T_2100 = -((gro*(6.0*FI_01*FI_10_p2*KS_p3 + 504.0*FI_p3*KS_1_p3 + 
          FI*KS_p2*(-7.0*FI_01*FIKS_2 - 4.0*FI_10*FI_11*KS - 
          2.0*FI_01*FI_20*KS + 14.0*FI_10*(2.0*FI_01 + FI_10)*KS_1) + 
          FI_p2*KS*(7.0*FIKS_3*KS + FI_21*KS_p2 - 
          7.0*KS_1*(24.0*FIKS_2 + 2.0*FI_11*KS + FI_20*KS - 
          8.0*(FI_01 + 2.0*FI_10)*KS_1) - 14.0*FI_10*KS*KS_2)))/(TDENOM3));

real  T_2010 = (2.0*FI_10_p2*KS_p2 + 56.0*FI_p2*KS_1_p2 - 
         FI*KS*(7.0*FIKS_2 + FI_20*KS - 14.0*FI_10*KS_1))/(TDENOM2);

real  T_2001=T_2010;


real  T_1110 = (2.0*FI_01*FI_10*KS_p2 + 56.0*FI_p2*KS_1_p2 + 
         FI*KS*(-7.0*FIKS_2 - FI_11*KS + 7.0*(FI_01 + FI_10)*KS_1))/(TDENOM2);

real  T_1101=T_1110;

real  T_1200 = -((gro*(6.0*FI_01_p2*FI_10*KS_p3 + 504.0*FI_p3*KS_1_p3 + 
         2.0*FI*KS_p2*(-7.0*FI_01*FIKS_2 - FI_02*FI_10*KS - 
         2.0*FI_01*FI_11*KS + 7.0*FI_01*(FI_01 + 2.0*FI_10)*KS_1) + 
         FI_p2*KS*(7.0*FIKS_3*KS + FI_12*KS_p2 - 
         7.0*KS_1*(24.0*FIKS_2 + FI_02*KS + 2.0*FI_11*KS - 
         8.0*(2.0*FI_01 + FI_10)*KS_1) - 7.0*FI_10*KS*KS_2)))/(TDENOM3));


real  T_0210 = (2.0*FI_01_p2*KS_p2 + 56.0*FI_p2*KS_1_p2 - 
         FI*KS*(7.0*FIKS_2 + FI_02*KS - 14.0*FI_01*KS_1))/TDENOM2;

real  T_0201=T_0210;

real  T_0300 = -((gro*(6.0*FI_01_p3*KS_p3 + 504.0*FI_p3*KS_1_p3 + 
         3.0*FI*FI_01*KS_p2*(-7.0*FIKS_2 - 2.0*FI_02*KS + 
         14.0*FI_01*KS_1) + FI_p2*KS*(7.0*FIKS_3*KS + 
         FI_03*KS_p2 + 21.0*KS_1*(-8.0*FIKS_2 - FI_02*KS + 
         8.0*FI_01*KS_1))))/(TDENOM3));


real  T_4000 = (gro*(24.0*FI_10_p4*KS_p4 + 
         12.0*FI*FI_10_p2*KS_p3*(-7.0*FIKS_2 - 3.0*FI_20*KS + 14.0*FI_10*KS_1) + 
         2.0*FI_p2*KS_p2*(3.0*FI_20_p2*KS_p2 + 336.0*FI_10_p2*KS_1_p2 + 
         2.0*FI_10*KS*(7.0*FIKS_3 + 2.0*FI_30*KS - 42.0*FI_20*KS_1)) + 
         FI_p3*KS*(-7.0*FIKS_4*KS_p2 - FI_40*KS_p3 + 
         28.0*KS_1*(FI_30*KS_p2 + 12.0*KS_1*(-9.0*FIKS_2 - FI_20*KS + 
         6.0*FI_10*KS_1)) + 42.0*KS*(FI_20*KS - 16.0*FI_10*KS_1)*KS_2) + 
         56.0*FI_p4*(90.0*KS_1_p4 + 3.0*KS_p2*KS_2_p2 + 
         4.0*KS_p2*KS_1*KS_3)))/(TDENOM4);

real  T_3100 = (gro*((6.0*FI*(-3.0*FI_10_p2*FI_11 + FI*FI_11*FI_20 + 
         FI*FI_10*FI_21) + 2.0*FI_01*(12.0*FI_10_p3 - 
         9.0*FI*FI_10*FI_20 + FI_p2*FI_30) - FI_p3*FI_31)*KS_p4 + 
         5040.0*FI_p4*KS_1_p4 + 504.0*FI_p3*KS*KS_1_p2*(-6.0*FIKS_2 + 
         (FI_01 + 3.0*FI_10)*KS_1) + 
         56.0*FI_p2*KS_p2*(KS_1*(4.0*FI*FIKS_3 - 
         3.0*(FI_10*FIKS_2 - 2.0*FI_10_p2*KS_1 + 
         FI*(FI_11 + FI_20)*KS_1 + FI_01*(FIKS_2 - 2.0*FI_10*KS_1))) - 
         6.0*FI*FI_10*KS_1*KS_2 + 
         3.0*FI_p2*KS_2_p2) + 7.0*FI*KS_p3*(6.0*FI_10_p3*KS_1 + 
         FI_01*(-6.0*FI_10*FIKS_2 + 18.0*FI_10_p2*KS_1 + 
         FI*(FIKS_3 - 6.0*FI_20*KS_1)) + 
         FI_p2*((3.0*FI_21 + FI_30)*KS_1 + 3.0*(FI_11 + FI_20)*KS_2) + 
         3.0*FI*FI_10*(FIKS_3 - 2.0*(2.0*FI_11*KS_1 + 
         FI_20*KS_1 + FI_10*KS_2)) - FI_p3*KS_4)))/(TDENOM4);

real  T_3010 = -((6.0*FI_10_p3*KS_p3 + 504.0*FI_p3*KS_1_p3 + 
         3.0*FI*FI_10*KS_p2*(-7.0*FIKS_2 - 2.0*FI_20*KS + 14.0*FI_10*KS_1) + 
         FI_p2*KS*(7.0*FIKS_3*KS + FI_30*KS_p2 + 
         21.0*KS_1*(-8.0*FIKS_2 - FI_20*KS + 8.0*FI_10*KS_1)))/(TDENOM3));

real  T_3001=T_3010;


real  T_2200 = (gro*((6.0*FI_01_p2*(4.0*FI_10_p2 - FI*FI_20) + 
         4.0*FI*FI_01*(-6.0*FI_10*FI_11 + FI*FI_21) + 
         FI*(FI_02*(-6.0*FI_10_p2 + 2.0*FI*FI_20) + 
         FI*(4.0*FI_11_p2 + 4.0*FI_10*FI_12 - FI*FI_22)))*KS_p4 + 
         5040.0*FI_p4*KS_1_p4 + 
         1008.0*FI_p3*KS*KS_1_p2*(-3.0*FIKS_2 + (FI_01 + FI_10)*KS_1) + 
         56.0*FI_p2*KS_p2*(2.0*FI_01_p2*KS_1_p2 + 2.0*FI_10_p2*KS_1_p2 - 
         2.0*FI_01*KS_1*(FIKS_2 - 4.0*FI_10*KS_1 + 2.0*FI*KS_2) - 
         FI*KS_1*(-2.0*FIKS_3 + (FI_02 + 4.0*FI_11 + FI_20)*KS_1 + 
         6.0*FI_10*KS_2) + FI_p2*(3.0*KS_2_p2 + 2.0*KS_1*KS_3)) + 
         7.0*FI*KS_p3*(-2.0*FI_01_p2*(FIKS_2 - 6.0*FI_10*KS_1) + 
         2.0*FI_01*(6.0*FI_10_p2*KS_1 + FI*(FIKS_3 - 
         2.0*(2.0*FI_11*KS_1 + FI_20*KS_1 + 2.0*FI_10*KS_2))) + 
         FI*(2.0*FI_10*FIKS_3 + FI_02*(FIKS_2 - 4.0*FI_10*KS_1) + 
         4.0*FI_11*(FIKS_2 - 2.0*FI_10*KS_1) - 2.0*FI_10_p2*KS_2 + 
         FI*(2.0*(FI_12 + FI_21)*KS_1 + FI_20*KS_2 - 
         FI*KS_4)))))/(TDENOM4);

real  T_2110 = (-6.0*FI_01*FI_10_p2*KS_p3 - 504.0*FI_p3*KS_1_p3 + 
         FI*KS_p2*(7.0*FI_01*FIKS_2 + 4.0*FI_10*FI_11*KS + 
         2.0*FI_01*FI_20*KS - 14.0*FI_10*(2.0*FI_01 + FI_10)*KS_1) + 
         FI_p2*KS*(-7.0*FIKS_3*KS - FI_21*KS_p2 + 
         7.0*KS_1*(24.0*FIKS_2 + 2.0*FI_11*KS + FI_20*KS - 
         8.0*(FI_01 + 2.0*FI_10)*KS_1) + 14.0*FI_10*KS*KS_2))/(TDENOM3);

real  T_2101=T_2110;


real  T_1300 = (gro*((24.0*FI_01_p3*FI_10 - 18.0*FI*FI_01_p2*FI_11 + 
         6.0*FI*FI_01*(-3.0*FI_02*FI_10 + FI*FI_12) + 
         FI_p2*(2.0*FI_03*FI_10 + 6.0*FI_02*FI_11 - FI*FI_13))*KS_p4 + 
         5040.0*FI_p4*KS_1_p4 + 
         504.0*FI_p3*KS*KS_1_p2*(-6.0*FIKS_2 + (3.0*FI_01 + FI_10)*KS_1) + 
         56.0*FI_p2*KS_p2*(6.0*FI_01_p2*KS_1_p2 + FI_01*(-9.0*FIKS_2*KS_1 + 
         6.0*FI_10*KS_1_p2) + FI*(3.0*FIKS_3*KS_1 - 
         3.0*(FI_02 + FI_11)*KS_1_p2 - 3.0*FI_10*KS_1*KS_2 + 
         3.0*FI*KS_2_p2 + FI*KS_1*KS_3)) + 
         7.0*FI*KS_p3*(6.0*FI_01_p3*KS_1 - 6.0*FI_01_p2*(FIKS_2 - 3.0*FI_10*KS_1) + 
         3.0*FI*FI_01*(FIKS_3 - 2.0*(FI_02*KS_1 + 2.0*FI_11*KS_1 + 
         FI_10*KS_2)) + FI*(-6.0*FI_02*FI_10*KS_1 + 
         FI*(-FIKS_4 + (FI_03 + 3.0*FI_12)*KS_1 + 
         3.0*(FI_02 + FI_11)*KS_2 + FI_10*KS_3)))))/(TDENOM4);

real  T_1210 = (-6.0*FI_01_p2*FI_10*KS_p3 - 504.0*FI_p3*KS_1_p3 + 
         2.0*FI*KS_p2*(7.0*FI_01*FIKS_2 + FI_02*FI_10*KS + 
         2.0*FI_01*FI_11*KS - 7.0*FI_01*(FI_01 + 2.0*FI_10)*KS_1) + 
         FI_p2*KS*(-7.0*FIKS_3*KS - FI_12*KS_p2 + 
         7.0*KS_1*(24.0*FIKS_2 + FI_02*KS + 2.0*FI_11*KS - 
          8.0*(2.0*FI_01 + FI_10)*KS_1) + 7.0*FI_10*KS*KS_2))/(TDENOM3);


real  T_1201=T_1210;


real  T_0400 = (gro*(24.0*FI_01_p4*KS_p4 + 
         12.0*FI*FI_01_p2*KS_p3*(-7.0*FIKS_2 - 3.0*FI_02*KS + 14.0*FI_01*KS_1) + 
         2.0*FI_p2*KS_p2*(3.0*FI_02_p2*KS_p2 + 336.0*FI_01_p2*KS_1_p2 + 
         2.0*FI_01*KS*(7.0*FIKS_3 + 2.0*FI_03*KS - 42.0*FI_02*KS_1)) + 
         FI_p3*KS*(-7.0*FIKS_4*KS_p2 - FI_04*KS_p3 + 
         28.0*KS_1*(FI_03*KS_p2 + 12.0*KS_1*(-9.0*FIKS_2 - FI_02*KS + 
         6.0*FI_01*KS_1)) + 42.0*KS*(FI_02*KS - 16.0*FI_01*KS_1)*KS_2) + 
         56.0*FI_p4*(90.0*KS_1_p4 + 3.0*KS_p2*KS_2_p2 + 
         4.0*KS_p2*KS_1*KS_3)))/(TDENOM4);

real  T_0310 = -((6.0*FI_01_p3*KS_p3 + 504.0*FI_p3*KS_1_p3 + 
         3.0*FI*FI_01*KS_p2*(-7.0*FIKS_2 - 2.0*FI_02*KS + 14.0*FI_01*KS_1) + 
         FI_p2*KS*(7.0*FIKS_3*KS + FI_03*KS_p2 + 
         21.0*KS_1*(-8.0*FIKS_2 - FI_02*KS + 8.0*FI_01*KS_1)))/(TDENOM3));

real  T_0301=T_0310;

/* Powers of T derivatives */
real  T_1000_p2 = T_1000*T_1000;
real  T_1000_p3 = T_1000_p2*T_1000;
real  T_1000_p4 = T_1000_p3*T_1000;

real  T_0100_p2 = T_0100*T_0100;
real  T_0100_p3 = T_0100_p2*T_0100;
real  T_0100_p4 = T_0100_p3*T_0100;

real  T_0010_p2 = T_0010*T_0010;
real  T_0010_p3 = T_0010_p2*T_0010;
real  T_0010_p4 = T_0010_p3*T_0010;

real  T_0001_p2 = T_0001*T_0001;
real  T_0001_p3 = T_0001_p2*T_0001;
real  T_0001_p4 = T_0001_p3*T_0001;

real  T_2000_p2 = T_2000*T_2000;

real  T_0200_p2 = T_0200*T_0200;

real  T_1100_p2 = T_1100*T_1100;

real  T_1010_p2 = T_1010*T_1010;

real  T_1001_p2 = T_1001*T_1001;

real  T_0110_p2 = T_0110*T_0110;

real  T_0101_p2 = T_0101*T_0101;

/* Derivatives of -EC/H0 */

real  p10 = (-(EC_10*H0) + EC*H0_10)/H0_p2;

real  p01 = (-(EC_01*H0) + EC*H0_01)/H0_p2;

real  p20 = (-(EC_20*H0_p2) + 2.0*EC_10*H0*H0_10 - 2.0*EC*H0_10_p2 + 
      EC*H0*H0_20)/H0_p3;

real  p11 = (-(EC_11*H0_p2) + EC_10*H0*H0_01 + EC_01*H0*H0_10 - 
      2.0*EC*H0_01*H0_10 + EC*H0*H0_11)/H0_p3;

real  p02 = (-(EC_02*H0_p2) + 2.0*EC_01*H0*H0_01 - 2.0*EC*H0_01_p2 + 
      EC*H0*H0_02)/H0_p3;

real  p30 = (-(EC_30*H0_p3) + 3.0*EC_20*H0_p2*H0_10 - 
      6.0*EC_10*H0*H0_10_p2 + 6.0*EC*H0_10_p3 + 
      3.0*EC_10*H0_p2*H0_20 - 6.0*EC*H0*H0_10*H0_20 + 
      EC*H0_p2*H0_30)/H0_p4;

real  p21 = (-(EC_21*H0_p3) + EC_20*H0_p2*H0_01 + 
      2.0*EC_11*H0_p2*H0_10 - 4.0*EC_10*H0*H0_01*H0_10 - 
      2.0*EC_01*H0*H0_10_p2 + 6.0*EC*H0_01*H0_10_p2 + 
      2.0*EC_10*H0_p2*H0_11 - 4.0*EC*H0*H0_10*H0_11 + 
      EC_01*H0_p2*H0_20 - 2.0*EC*H0*H0_01*H0_20 + 
      EC*H0_p2*H0_21)/H0_p4;

real  p12 = (-(EC_12*H0_p3) + 2.0*EC_11*H0_p2*H0_01 - 
      2.0*EC_10*H0*H0_01_p2 + EC_10*H0_p2*H0_02 + 
      EC_02*H0_p2*H0_10 - 4.0*EC_01*H0*H0_01*H0_10 + 
      6.0*EC*H0_01_p2*H0_10 - 2.0*EC*H0*H0_02*H0_10 + 
      2.0*EC_01*H0_p2*H0_11 - 4.0*EC*H0*H0_01*H0_11 + 
      EC*H0_p2*H0_12)/H0_p4;

real  p03 = (-(EC_03*H0_p3) + 3.0*EC_02*H0_p2*H0_01 - 
      6.0*EC_01*H0*H0_01_p2 + 6.0*EC*H0_01_p3 + 
      3.0*EC_01*H0_p2*H0_02 - 6.0*EC*H0*H0_01*H0_02 + 
      EC*H0_p2*H0_03)/H0_p4;

real  p40 = (-(EC_40*H0_p4) + 4.0*EC_30*H0_p3*H0_10 - 
      12.0*EC_20*H0_p2*H0_10_p2 + 24.0*EC_10*H0*H0_10_p3 - 
      24.0*EC*H0_10_p4 + 6.0*EC_20*H0_p3*H0_20 - 
      24.0*EC_10*H0_p2*H0_10*H0_20 + 36.0*EC*H0*H0_10_p2*H0_20 - 
      6.0*EC*H0_p2*H0_20_p2 + 4.0*EC_10*H0_p3*H0_30 - 
      8.0*EC*H0_p2*H0_10*H0_30 + EC*H0_p3*H0_40)/H0_p5;

real  p31 = (-(EC_31*H0_p4) + EC_30*H0_p3*H0_01 + 
      3.0*EC_21*H0_p3*H0_10 - 6.0*EC_20*H0_p2*H0_01*H0_10 - 
      6.0*EC_11*H0_p2*H0_10_p2 + 18.0*EC_10*H0*H0_01*H0_10_p2 + 
      6.0*EC_01*H0*H0_10_p3 - 24.0*EC*H0_01*H0_10_p3 + 
      3.0*EC_20*H0_p3*H0_11 - 12.0*EC_10*H0_p2*H0_10*H0_11 + 
      18.0*EC*H0*H0_10_p2*H0_11 + 3.0*EC_11*H0_p3*H0_20 - 
      6.0*EC_10*H0_p2*H0_01*H0_20 - 6.0*EC_01*H0_p2*H0_10*H0_20 + 
      18.0*EC*H0*H0_01*H0_10*H0_20 - 6.0*EC*H0_p2*H0_11*H0_20 + 
      3.0*EC_10*H0_p3*H0_21 - 6.0*EC*H0_p2*H0_10*H0_21 + 
      EC_01*H0_p3*H0_30 - 2.0*EC*H0_p2*H0_01*H0_30 + 
      EC*H0_p3*H0_31)/H0_p5;

real  p22 = (-(EC_22*H0_p4) + 2.0*EC_21*H0_p3*H0_01 - 
      2.0*EC_20*H0_p2*H0_01_p2 + EC_20*H0_p3*H0_02 + 
      2.0*EC_12*H0_p3*H0_10 - 8.0*EC_11*H0_p2*H0_01*H0_10 + 
      12.0*EC_10*H0*H0_01_p2*H0_10 - 4.0*EC_10*H0_p2*H0_02*H0_10 - 
      2.0*EC_02*H0_p2*H0_10_p2 + 12.0*EC_01*H0*H0_01*H0_10_p2 - 
      24.0*EC*H0_01_p2*H0_10_p2 + 6.0*EC*H0*H0_02*H0_10_p2 + 
      4.0*EC_11*H0_p3*H0_11 - 8.0*EC_10*H0_p2*H0_01*H0_11 - 
      8.0*EC_01*H0_p2*H0_10*H0_11 + 24.0*EC*H0*H0_01*H0_10*H0_11 - 
      4.0*EC*H0_p2*H0_11_p2 + 2.0*EC_10*H0_p3*H0_12 - 
      4.0*EC*H0_p2*H0_10*H0_12 + EC_02*H0_p3*H0_20 - 
      4.0*EC_01*H0_p2*H0_01*H0_20 + 6.0*EC*H0*H0_01_p2*H0_20 - 
      2.0*EC*H0_p2*H0_02*H0_20 + 2.0*EC_01*H0_p3*H0_21 - 
      4.0*EC*H0_p2*H0_01*H0_21 + EC*H0_p3*H0_22)/H0_p5;

real  p13 = (-(EC_13*H0_p4) + 3.0*EC_12*H0_p3*H0_01 - 
      6.0*EC_11*H0_p2*H0_01_p2 + 6.0*EC_10*H0*H0_01_p3 + 
      3.0*EC_11*H0_p3*H0_02 - 6.0*EC_10*H0_p2*H0_01*H0_02 + 
      EC_10*H0_p3*H0_03 + EC_03*H0_p3*H0_10 - 
      6.0*EC_02*H0_p2*H0_01*H0_10 + 18.0*EC_01*H0*H0_01_p2*H0_10 - 
      24.0*EC*H0_01_p3*H0_10 - 6.0*EC_01*H0_p2*H0_02*H0_10 + 
      18.0*EC*H0*H0_01*H0_02*H0_10 - 2.0*EC*H0_p2*H0_03*H0_10 + 
      3.0*EC_02*H0_p3*H0_11 - 12.0*EC_01*H0_p2*H0_01*H0_11 + 
      18.0*EC*H0*H0_01_p2*H0_11 - 6.0*EC*H0_p2*H0_02*H0_11 + 
      3.0*EC_01*H0_p3*H0_12 - 6.0*EC*H0_p2*H0_01*H0_12 + 
      EC*H0_p3*H0_13)/H0_p5;

real  p04 = (-(EC_04*H0_p4) + 4.0*EC_03*H0_p3*H0_01 - 
      12.0*EC_02*H0_p2*H0_01_p2 + 24.0*EC_01*H0*H0_01_p3 - 
      24.0*EC*H0_01_p4 + 6.0*EC_02*H0_p3*H0_02 - 
      24.0*EC_01*H0_p2*H0_01*H0_02 + 36.0*EC*H0*H0_01_p2*H0_02 - 
      6.0*EC*H0_p2*H0_02_p2 + 4.0*EC_01*H0_p3*H0_03 - 
      8.0*EC*H0_p2*H0_01*H0_03 + EC*H0_p3*H0_04)/H0_p5;

/* Powers of p=-EC/H0 */

real  p10_p2 = p10*p10;
real  p10_p3 = p10_p2*p10;
real  p10_p4 = p10_p3*p10;

real  p01_p2 = p01*p01;
real  p01_p3 = p01_p2*p01;
real  p01_p4 = p01_p3*p01;

real  p20_p2 = p20*p20;
real  p02_p2 = p02*p02;

real  p11_p2 = p11*p11;


/* Derivatives of AEXP */

real  AEXP_10 = AEXP*p10;

real  AEXP_01 = AEXP*p01;

real  AEXP_20 = AEXP*(p10_p2 + p20);

real  AEXP_11 = AEXP*(p10*p01 + p11);

real  AEXP_02 = AEXP*(p01_p2 + p02);

real  AEXP_30 = AEXP*(p10_p3 + 3.0*p10*p20 + p30);

real  AEXP_21 = AEXP*(2.0*p10*p11 + p01*(p10_p2 + p20) + p21);

real  AEXP_12 = AEXP*(2.0*p01*p11 + p10*(p01_p2 + p02) + p12);

real  AEXP_03 = AEXP*(p01_p3 + 3.0*p01*p02 + p03);

real  AEXP_40 = AEXP*(p10_p4 + 6.0*p10_p2*p20 + 3.0*p20_p2 + 4.0*p10*p30 + p40);

real  AEXP_31 = AEXP*(3.0*p10_p2*p11 + 3.0*p11*p20 + 3.0*p10*p21 + 
          p01*(p10_p3 + 3.0*p10*p20 + p30) + p31);

real  AEXP_22 = AEXP*(2.0*p11_p2 + 2.0*p10*p12 + p01_p2*(p10_p2 + p20) + 
          p02*(p10_p2 + p20) + 2.0*p01*(2.0*p10*p11 + p21) + p22);

real  AEXP_13 = AEXP*(3.0*p01_p2*p11 + 3.0*p11*p02 + 3.0*p01*p12 + 
          p10*(p01_p3 + 3.0*p01*p02 + p03) + p13);

real  AEXP_04 = AEXP*(p01_p4 + 6.0*p01_p2*p02 + 3.0*p02_p2 + 4.0*p01*p03 + p04);


real  AEXP_10_p2 = AEXP_10*AEXP_10;
real  AEXP_10_p3 =AEXP_10_p2*AEXP_10;
real  AEXP_10_p4 =AEXP_10_p3*AEXP_10;

real  AEXP_01_p2 = AEXP_01*AEXP_01;
real  AEXP_01_p3 =AEXP_01_p2*AEXP_01;
real  AEXP_01_p4 =AEXP_01_p3*AEXP_01;

real  AEXP_20_p2 = AEXP_20*AEXP_20;
real  AEXP_02_p2 = AEXP_02*AEXP_02;

real  AEXP_11_p2 = AEXP_11*AEXP_11;



/* Derivatives of A */

real  A_10 = -((AEXP_10*COR_delta)/ADENOM_p2);

real  A_01 = -((AEXP_01*COR_delta)/ADENOM_p2);

real  A_20 = ((2.0*AEXP_10_p2 - ADENOM*AEXP_20)*COR_delta)/ADENOM_p3;

real  A_11 = ((2.0*AEXP_10*AEXP_01 - ADENOM*AEXP_11)*COR_delta)/ADENOM_p3;

real  A_02 = ((2.0*AEXP_01_p2 - ADENOM*AEXP_02)*COR_delta)/ADENOM_p3;

real  A_30 = -(((6.0*AEXP_10_p3 - 6.0*ADENOM*AEXP_10*AEXP_20 + 
        ADENOM_p2*AEXP_30)*COR_delta)/ADENOM_p4);

real  A_21 = ((-6.0*AEXP_01*AEXP_10_p2 + 4.0*ADENOM*AEXP_10*AEXP_11 + 
       2.0*ADENOM*AEXP_01*AEXP_20 - ADENOM_p2*AEXP_21)*COR_delta)/ADENOM_p4;

real  A_12 = ((-6.0*AEXP_10*AEXP_01_p2 + 4.0*ADENOM*AEXP_01*AEXP_11 + 
       2.0*ADENOM*AEXP_10*AEXP_02 - ADENOM_p2*AEXP_12)*COR_delta)/ADENOM_p4;


real  A_03 = -(((6.0*AEXP_01_p3 - 6.0*ADENOM*AEXP_01*AEXP_02 + 
        ADENOM_p2*AEXP_03)*COR_delta)/ADENOM_p4);

real  A_40 = ((24.0*AEXP_10_p4 - 36.0*ADENOM*AEXP_10_p2*AEXP_20 + 
       8.0*ADENOM_p2*AEXP_10*AEXP_30 + ADENOM_p2*(6.0*AEXP_20_p2 - 
       ADENOM*AEXP_40))*COR_delta)/ADENOM_p5;

real  A_31 = ((6.0*ADENOM*(-3.0*AEXP_10_p2*AEXP_11 + ADENOM*AEXP_11*AEXP_20 + 
       ADENOM*AEXP_10*AEXP_21) + 
       2.0*AEXP_01*(12.0*AEXP_10_p3 - 9.0*ADENOM*AEXP_10*AEXP_20 + 
       ADENOM_p2*AEXP_30) - ADENOM_p3*AEXP_31)*COR_delta)/ADENOM_p5;

real  A_22 = ((6.0*AEXP_01_p2*(4.0*AEXP_10_p2 - ADENOM*AEXP_20) + 
       4.0*ADENOM*AEXP_01*(-6.0*AEXP_10*AEXP_11 + ADENOM*AEXP_21) + 
       ADENOM*(AEXP_02*(-6.0*AEXP_10_p2 + 2.0*ADENOM*AEXP_20) + 
       ADENOM*(4.0*AEXP_11_p2 + 4.0*AEXP_10*AEXP_12 - 
       ADENOM*AEXP_22)))*COR_delta)/ADENOM_p5;

real  A_13 = ((6.0*ADENOM*(-3.0*AEXP_01_p2*AEXP_11 + 
       ADENOM*AEXP_11*AEXP_02 + ADENOM*AEXP_01*AEXP_12) + 
       2.0*AEXP_10*(12.0*AEXP_01_p3 - 9.0*ADENOM*AEXP_01*AEXP_02 + 
       ADENOM_p2*AEXP_03) - ADENOM_p3*AEXP_13)*COR_delta)/ADENOM_p5;

real  A_04 = ((24.0*AEXP_01_p4 - 36.0*ADENOM*AEXP_01_p2*AEXP_02 + 
       8.0*ADENOM_p2*AEXP_01*AEXP_03 + ADENOM_p2*(6.0*AEXP_02_p2 - 
       ADENOM*AEXP_04))*COR_delta)/ADENOM_p5;

/* Powers of A derivatives */
real  A_10_p2 = A_10*A_10;
real  A_10_p3 = A_10_p2*A_10;
real  A_10_p4 = A_10_p3*A_10;

real  A_01_p2 = A_01*A_01;
real  A_01_p3 = A_01_p2*A_01;
real  A_01_p4 = A_01_p3*A_01;

real  A_20_p2 = A_20*A_20;

real  A_02_p2 =  A_02*A_02;

/* Derivatives of H with respect to H0, T, A */
/* Derivatives with respect to H0 higher then first order vanish */
/* Terms used in the derivatives */
real  At2 = A*T_p2;
real  At2_p2=At2*At2;

real  At_nom = 1.0 + At2;
real  At_denom = At_nom + At2_p2;

real  At_term0 = 1.0 + 2.0*At2;
real  At_term1 = T + 3.0*A*T_p3 + 2.0*A_p2*T_p5;
real  At_term2 = A + 3.0*At2*A + 2.0*A_p3*T_p4;
real  At_term3 = T + 2.0*A*T_p3;
real  At_term4 = COR_delta + 2.0*At2*COR_delta;
real  At_term5 = H1 + 2.0*At2*H1;
real  At_term6 = 1.0 + 3.0*At2 + 2.0*At2_p2;
real  At_term7 = COR_delta + At2*COR_delta;
real  At_term8 = At_denom + 2.0*At2*At_denom;


/* Powers of terms */
real  At_nom_p2 = At_nom*At_nom;
real  At_nom_p3 = At_nom_p2*At_nom;

real  At_denom_p2 = At_denom*At_denom;
real  At_denom_p3 = At_denom_p2*At_denom;
real  At_denom_p4 = At_denom_p3*At_denom;
real  At_denom_p5 = At_denom_p4*At_denom;
real  At_denom_p6 = At_denom_p4*At_denom_p2;
real  At_denom_p7 = At_denom_p4*At_denom_p3;


real  HDENOM = At_denom_p2*H1;

real  HDENOM_p2 = HDENOM*HDENOM;
real  HDENOM_p3 = HDENOM_p2*HDENOM;
real  HDENOM_p4 = HDENOM_p3*HDENOM;

real  At_term0_p2 = At_term0*At_term0;
real  At_term0_p3 = At_term0_p2*At_term0;
real  At_term0_p4 = At_term0_p3*At_term0;

real  At_term1_p2 = At_term1*At_term1;

real  At_term2_p2 = At_term2*At_term2;
real  At_term2_p4 = At_term2_p2*At_term2_p2;

real  At_term3_p2 = At_term3*At_term3;

real  At_term4_p2 = At_term4*At_term4;

real  At_term5_p2 = At_term5*At_term5;
real  At_term5_p3 = At_term5_p2*At_term5;

real  At_term6_p2 = At_term6*At_term6;
real  At_term6_p3 = At_term6_p2*At_term6;
real  At_term6_p4 = At_term6_p3*At_term6;

real  At_term7_p2 = At_term7*At_term7;
real  At_term7_p3 = At_term7_p2*At_term7;

real  At_term8_p2 = At_term8*At_term8;

real  H1_p2 = H1*H1;
real  H1_p3 = H1_p2*H1;

/* First order derivatives */
real  H_A = ((At_denom - At_nom*At_term0)*COR_delta*H0*T_p4)/HDENOM;

real  H_T = (2.0*(At_denom - At2*At_nom)*At_term0*COR_delta*H0*T)/HDENOM;

real  H_H = log(H1);

/* Second order derivatives */
real  H_AA = (COR_delta*H0*T_p6*(-(At_term1_p2*COR_delta) + 
       2.0*At_denom*At_nom*At_term0*(H1 + (2.0*H1*A + COR_delta)*T_p2) - 
       At_denom_p2*(4.0*H1 + (6.0*H1*A + COR_delta)*T_p2)))/HDENOM_p2;

real  H_TT = (2.0*COR_delta*H0*(At_denom_p3*H1 + 
       At2*At_denom_p2*(-5.0 + 6.0*At_denom)*H1 + 
       4.0*At_denom*At_nom*At_term0_p2*A*(H1*A + COR_delta)*T_p4 - 
       2.0*At_term2_p2*COR_delta*T_p6 - 
       At_denom_p2*((23.0 + 22.0*At2)*H1*A_p2*T_p4 + 
       2.0*COR_delta*At_term3_p2)))/HDENOM_p2;

real  H_HT = (2.0*(At_denom - At2*At_nom)*At_term0*COR_delta*T)/HDENOM;


real  H_TA = (2.0*(At_denom - At2*At_nom)*COR_delta*H0*T_p3*(2.0*At_denom_p2*H1 - 
       At_denom*At_term0*(2.0*H1 + (4.0*H1*A + COR_delta)*T_p2) + 
       At_nom*COR_delta*At_term3_p2))/HDENOM_p2;

real  H_HA = ((At_denom - At_nom*At_term0)*COR_delta*T_p4)/HDENOM;


/* Third order derivatives */
real  H_TTT = (4.0*(At_denom - At2*At_nom)*COR_delta*H0*T*(3.0*At_denom_p3*H1*(2.0*At_denom*H1*A - 
        (1.0 + 6.0*At2)*At_term0*(2.0*H1*A + COR_delta)) + 
        At_denom_p2*At_term0*(12.0*At_term0_p2*H1_p2*A_p2 + 
        3.0*(5.0 + At2*(23.0 + 22.0*At2))*H1*A*COR_delta + 4.0*At_term4_p2)*T_p2 - 
        4.0*At_denom*At_nom*At_term0_p3*A*COR_delta*(3.0*H1*A + 2.0*COR_delta)*T_p4 + 
        4.0*At_nom_p2*At_term0_p3*A_p2*COR_delta_p2*T_p6))/HDENOM_p3;

real  H_AAA = (-2.0*COR_delta*H0*T_p8*(3.0*At_denom_p4*H1_p2 - 3.0*At_denom*At_term1_p2*COR_delta*(H1 + 
        (2.0*H1*A + COR_delta)*T_p2) + COR_delta_p2*T_p4*At_term6_p3 + 
        3.0*At_denom_p2*At_nom*At_term0*(H1_p2 + H1*(4.0*H1*A + 3.0*COR_delta)*T_p2 + 
        (H1*A + COR_delta)*(4.0*H1*A + COR_delta)*T_p4) - At_denom_p3*(9.0*H1_p2 + 
        6.0*H1*(5.0*H1*A + COR_delta)*T_p2 + (24.0*H1_p2*A_p2 + 9.0*H1*A*COR_delta + 
        COR_delta_p2)*T_p4)))/HDENOM_p3;

real  H_HTT = (2.0*COR_delta*(At_denom_p3*H1 + At2*At_denom_p2*(-5.0 + 6.0*At_denom)*H1 + 
        4.0*At_denom*At_nom*At_term0_p2*A*(H1*A + COR_delta)*T_p4 - 
        2.0*At_term2_p2*COR_delta*T_p6 - At_denom_p2*((23.0 + 22.0*At2)*H1*A_p2*T_p4 + 
        2.0*COR_delta*At_term3_p2)))/HDENOM_p2;

real  H_HAA = (COR_delta*T_p6*(-(At_term1_p2*COR_delta) + 
        2.0*At_denom*At_nom*At_term0*(H1 + (2.0*H1*A + COR_delta)*T_p2) - 
        At_denom_p2*(4.0*H1 + (6.0*H1*A + COR_delta)*T_p2)))/HDENOM_p2; 

real  H_TTA = (2.0*(At_denom - At2*At_nom)*COR_delta*H0*T_p2*(6.0*At_denom_p3*(At_denom - 
        (1.0 + 6.0*At2)*At_term0)*H1_p2 - 
        At_denom_p2*H1*(-12.0*At_term0_p3*H1*A + 
        (9.0 + 22.0*At2)*At_denom*COR_delta - 3.0*(3.0 + At2*(17.0 + 18.0*At2))*At_term0*COR_delta)*T_p2 + 
        4.0*At_denom*At_term0_p2*COR_delta*(At_denom*COR_delta - 
        At_nom*(3.0*H1*A + COR_delta))*T_p4 + 
        4.0*At_term0_p2*A*COR_delta*(At_nom_p2*At_term0*COR_delta - 
        3.0*At_denom*At_nom*(2.0*H1*A + COR_delta))*T_p6))/HDENOM_p3;


real  H_TAA = (4.0*(At_denom - At2*At_nom)*COR_delta*H0*T_p5*(At_nom_p2*At_term0_p3*COR_delta_p2*T_p4 + 
        At_denom*At_nom*At_term3_p2*COR_delta*(-3.0*H1 - 2.0*(3.0*H1*A + COR_delta)*T_p2) - 
        2.0*At_denom_p3*H1*(3.0*H1 + (6.0*H1*A + COR_delta)*T_p2) + 
        At_denom_p2*At_term0*(3.0*At_term5_p2 + 3.0*(2.0 + 3.0*At2)*H1*COR_delta*T_p2 + 
        COR_delta_p2*T_p4)))/HDENOM_p3;

real  H_HTA = (2.0*(At_denom - At2*At_nom)*COR_delta*T_p3*(2.0*At_denom*(At_denom - At_term0)*H1 + 
        At_nom*At_term3_p2*COR_delta - 
        At_denom*At_term0*(4.0*H1*A + COR_delta)*T_p2))/HDENOM_p2;


/* Fourth order derivatives */
real  H_TTTT = (-12.0*COR_delta*H0*(-2.0*At_denom_p7*H1_p3*A + 
         At_denom_p4*(8.0*(7.0 + At2*(37.0 + 34.0*At2))*At_term4_p2*H1*A + 
         8.0*(7.0 + At2*(37.0 + 34.0*At2))*At_term0_p2*H1_p3*A_p3 + 
         (97.0 + At2*(942.0 + At2*(3101.0 + 4.0*At2*(1049.0 + 497.0*At2))))*H1_p2*A_p2*COR_delta + 
         8.0*At_term0_p4*COR_delta_p3)*T_p4 - 
         8.0*At_denom_p3*At_nom*At_term0_p2*A*(H1*A + COR_delta)*(4.0*At_term4_p2 + 
         H1*A*(4.0*At_term0_p2*H1*A + (11.0 + At2*(53.0 + 50.0*At2))*COR_delta))*T_p6 + 
         8.0*At_denom_p2*At_nom_p2*At_term0_p2*A_p2*COR_delta*(6.0*At_term4_p2 + 
         H1*A*(6.0*At_term0_p2*H1*A + 
         (13.0 + At2*(55.0 + 54.0*At2))*COR_delta))*T_p8 - 
         32.0*At_denom*At_nom_p3*At_term0_p4*A_p3*COR_delta_p2*(H1*A + COR_delta)*T_p10 + 
         8.0*At_term2_p4*COR_delta_p3*T_p12 + At_denom_p6*H1_p2*(COR_delta + 28.0*At2*COR_delta + 
         2.0*A*(H1 + 29.0*At2*H1 + A*(69.0*H1*A + 34.0*COR_delta)*T_p4)) - 
         2.0*At_denom_p5*H1*(13.0*At2*H1*COR_delta + T_p2*(4.0*(1.0 + 6.0*At2)*At_term4_p2 + 
         H1*A_p2*((13.0 + At2*(149.0 + 432.0*At2))*H1 + 
         (149.0*COR_delta + 432.0*At2*COR_delta + 356.0*At2_p2*(H1*A + COR_delta))*T_p2)))))/HDENOM_p4;

real  H_AAAA = (-6.0*COR_delta*H0*T_p10*(At_term6_p4*COR_delta_p3*T_p6 - 
         4.0*At_denom*At_nom_p3*At_term0_p3*COR_delta_p2*T_p4*(H1 + (2.0*H1*A + COR_delta)*T_p2) - 
         4.0*At_denom_p5*H1_p2*(3.0*H1 + (5.0*H1*A + COR_delta)*T_p2) + 
         2.0*At_denom_p2*At_term1_p2*COR_delta*(3.0*At_term5_p2 + 
         COR_delta*T_p2*(2.0*(4.0 + 7.0*At2)*H1 + 3.0*COR_delta*T_p2)) + 
         At_denom_p4*(4.0*(4.0 + 5.0*At2)*At_term0_p2*H1_p3 + 
         COR_delta*T_p2*((24.0 + 76.0*At2 + 58.0*At2_p2)*H1_p2 + 
         4.0*(2.0 + 3.0*At2)*H1*COR_delta*T_p2 + COR_delta_p2*T_p4)) - 
         4.0*At_denom_p3*At_nom*At_term0*(At_term5_p3 + 
         COR_delta*T_p2*(3.0*(2.0 + At2*(7.0 + 6.0*At2))*H1_p2 + 
         (5.0 + 8.0*At2)*H1*COR_delta*T_p2 + COR_delta_p2*T_p4))))/HDENOM_p4;

real  H_HTTT = (4.0*(At_denom - At2*At_nom)*COR_delta*T*(3.0*At_denom_p3*H1*(2.0*At_denom*H1*A - 
         (1.0 + 6.0*At2)*At_term0*(2.0*H1*A + COR_delta)) + 
         At_denom_p2*At_term0*(4.0*At_term4_p2 + 3.0*H1*A*(4.0*At_term0_p2*H1*A + 
         (5.0 + At2*(23.0 + 22.0*At2))*COR_delta))*T_p2 - 
         4.0*At_denom*At_nom*At_term0_p3*A*COR_delta*(3.0*H1*A + 2.0*COR_delta)*T_p4 + 
         4.0*At_nom_p2*At_term0_p3*A_p2*COR_delta_p2*T_p6))/HDENOM_p3;

real  H_TTTA = (24.0*(At_denom - At2*At_nom)*COR_delta*H0*T*(At_denom_p6*H1_p3 - 
         At_denom_p3*At_term0*(8.0*At_term5_p3*A_p2 + 
         COR_delta*(3.0*(5.0 + At2*(27.0 + 26.0*At2))*At_term0*H1_p2*A + 
         COR_delta*((5.0 + At2*(46.0 + 115.0*At2 + 82.0*At2_p2))*H1 + 
         2.0*At_term3_p2*COR_delta)))*T_p4 -  
         2.0*At_denom*At_term0_p3*At_term7_p2*A*(4.0*H1*(A + 2.0*At2*A) + 
         (2.0 + 5.0*At2)*COR_delta)*T_p8 + 2.0*At_term0_p4*At_term7_p3*A_p2*T_p10 -  
         At_denom_p5*H1_p2*(H1 + 2.0*(9.0*H1*A + 2.0*COR_delta)*T_p2 +  
         A*(40.0*H1*A + 13.0*COR_delta)*T_p4) + 
         At_denom_p4*H1*T_p2*((7.0 + 38.0*At2)*At_term5_p2*A + 
         COR_delta*((4.0 + At2*(50.0 + 155.0*At2 + 138.0*At2_p2))*H1 + 
         (5.0 + 14.0*At2)*At_term0*COR_delta*T_p2)) + 
         At_nom*At_term8_p2*COR_delta*T_p6*(12.0*At_term0_p2*H1_p2*A_p2 + 
         COR_delta*(2.0*COR_delta*(1.0 + 6.0*At2 + 8.0*A_p2*T_p4) + 
         H1*A*(13.0 + 63.0*At2 + 66.0*A_p2*T_p4)))))/HDENOM_p4;

real  H_TAAA = (12.0*(At_denom - At2*At_nom)*COR_delta*H0*T_p7*(-4.0*At_denom_p5*H1_p3 + 
         At_term0_p4*At_term7_p3*T_p6 - 
         At_denom*At_term0_p3*At_term7_p2*T_p4*(4.0*H1 + (8.0*H1*A + 3.0*COR_delta)*T_p2) - 
         At_denom_p3*At_term0*(2.0*H1 + (4.0*H1*A + COR_delta)*T_p2)*(2.0*At_term5_p2 + 
         COR_delta*T_p2*(8.0*H1 + 10.0*At2*H1 + COR_delta*T_p2)) + 
         At_denom_p4*H1*(12.0*At_term5_p2 + 
         COR_delta*T_p2*(11.0*H1 + 20.0*At2*H1 + 2.0*COR_delta*T_p2)) + 
         At_denom_p2*At_nom*At_term3_p2*COR_delta*(6.0*At_term5_p2 + 
         COR_delta*T_p2*(4.0*(3.0 + 5.0*At2)*H1 + 3.0*COR_delta*T_p2))))/HDENOM_p4;

real  H_HAAA = (-2.0*COR_delta*T_p8*(3.0*At_denom_p4*H1_p2 + At_term6_p3*COR_delta_p2*T_p4 - 
       3.0*At_denom*At_term1_p2*COR_delta*(H1 + (2.0*H1*A + COR_delta)*T_p2) + 
       3.0*At_denom_p2*At_nom*At_term0*(H1_p2 + H1*(4.0*H1*A + 3.0*COR_delta)*T_p2 + 
       (H1*A + COR_delta)*(4.0*H1*A + COR_delta)*T_p4) - 
       At_denom_p3*(9.0*H1_p2 + 6.0*H1*(5.0*H1*A + COR_delta)*T_p2 + 
       (24.0*H1_p2*A_p2 + 9.0*H1*A*COR_delta + COR_delta_p2)*T_p4)))/HDENOM_p3;

real  H_HATT = (2.0*(At_denom - At2*At_nom)*COR_delta*T_p2*(6.0*At_denom_p3*(At_denom - 
         (1.0 + 6.0*At2)*At_term0)*H1_p2 - 
         At_denom_p2*(-12.0*At_term0*At_term5_p2*A + 
         (9.0 + 22.0*At2)*At_denom*H1*COR_delta - 
         3.0*(3.0 + At2*(17.0 + 18.0*At2))*At_term0*H1*COR_delta)*T_p2 + 
         4.0*At_denom*At_term0_p2*COR_delta*(At_denom*COR_delta - At_nom*(3.0*H1*A + COR_delta))*T_p4 + 
         4.0*At_term0_p2*A*(At_term0*At_term7_p2 - 
         3.0*At_denom*At_nom*COR_delta*(2.0*H1*A + COR_delta))*T_p6))/HDENOM_p3;


real  H_HTAA = (4.0*(At_denom - At2*At_nom)*COR_delta*T_p5*(At_term0_p3*At_term7_p2*T_p4 + 
         At_denom*At_nom*At_term3_p2*COR_delta*(-3.0*H1 - 2.0*(3.0*H1*A + COR_delta)*T_p2) - 
         2.0*At_denom_p3*H1*(3.0*H1 + (6.0*H1*A + COR_delta)*T_p2) + 
         At_denom_p2*At_term0*(3.0*At_term5_p2 + 
         COR_delta*T_p2*(6.0*H1 + 9.0*At2*H1 + COR_delta*T_p2))))/HDENOM_p3;

real  H_TTAA = (-4.0*(At_denom - At2*At_nom)*COR_delta*H0*T_p4*(-6.0*At_term0_p4*At_term7_p3*A*T_p8 + 
          2.0*At_denom_p5*H1_p2*(15.0*H1 + 7.0*(6.0*H1*A + COR_delta)*T_p2) - 
          3.0*At_term8_p2*At_nom*COR_delta*T_p4*(12.0*At_term5_p2*A + 
          COR_delta*((7.0 + 9.0*At2*(5.0 + 6.0*At2))*H1 + 
          2.0*(2.0 + 5.0*At2)*COR_delta*T_p2)) - At_denom_p4*H1*(3.0*(5.0 + 34.0*At2)*At_term5_p2 + 
          COR_delta*T_p2*(3.0*(18.0 + At2*(83.0 + 90.0*At2))*H1 + 
          (17.0 + 38.0*At2)*COR_delta*T_p2)) + 
          6.0*At_denom*At_term0_p3*At_term7_p2*T_p6*(COR_delta + 
          4.0*A*(H1 + (2.0*H1*A + COR_delta)*T_p2)) + 
          3.0*At_denom_p3*At_term0*T_p2*(8.0*At_term5_p3*A + 
          COR_delta*(3.0*(3.0 + At2*(21.0 + 22.0*At2))*At_term0*H1_p2 + 
          2.0*(7.0 + At2*(29.0 + 26.0*At2))*H1*COR_delta*T_p2 + 
          2.0*At_term0*COR_delta_p2*T_p4))))/HDENOM_p4;


/* Final expressions (general derivatives of F[H[roa,rob], T[roa,rob,groa,grob], A[roa,rob]]) */

/* First order derivatives */
real      df10000  = A_10*H_A+H0_10*H_H+H_T*T_1000 + EC_10;
ds->df10000 += factor*(E + roa_rob*df10000);

real      df01000  = A_01*H_A+H0_01*H_H+H_T*T_0100 + EC_01;
ds->df01000 += factor*(E + roa_rob*df01000);

real      df00100  = H_T*T_0010;
ds->df00100 += factor*roa_rob*df00100;

real      df00010  = H_T*T_0001;
ds->df00010 += factor*roa_rob*df00010;

/* Second order derivatives */
real      df20000  = A_20*H_A+H0_20*H_H+A_10*(A_10*H_AA+2.0*H0_10*H_HA)+
               2.0*(A_10*H_TA+H0_10*H_HT)*T_1000+H_TT*T_1000_p2+H_T*T_2000 + EC_20;
ds->df20000 += factor*(2.0*df10000 + roa_rob*df20000);

real      df11000  = A_11*H_A+H0_11*H_H+A_01*(A_10*H_AA+H0_10*H_HA+H_TA*T_1000)+T_0100*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+H0_01*(A_10*H_HA+H_HT*T_1000)+H_T*T_1100 + EC_11;
ds->df11000 += factor*(df10000 + df01000 + roa_rob*df11000);

real      df10100  = T_0010*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+H_T*T_1010;
ds->df10100  += factor*(df00100 + roa_rob*df10100);

real      df10010  = T_0001*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+H_T*T_1001;
ds->df10010  += factor*(df00010 + roa_rob*df10010);

real      df02000  = A_02*H_A+H0_02*H_H+A_01*(A_01*H_AA+2.0*H0_01*H_HA)+2.0*(A_01*H_TA+H0_01*H_HT)*T_0100+H_TT*T_0100_p2+H_T*T_0200 + EC_02;
ds->df02000 += factor*(2.0*df01000 + roa_rob*df02000);

real      df01100  = T_0010*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)+H_T*T_0110;
ds->df01100 += factor*(df00100 + roa_rob*df01100);

real      df01010  = T_0001*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)+H_T*T_0101;
ds->df01010 += factor*(df00010 + roa_rob*df01010);

real      df00200  = H_TT*T_0010_p2;
ds->df00200 += factor*roa_rob*df00200;

real      df00110  = H_TT*T_0001*T_0010;
ds->df00110 += factor*roa_rob*df00110;
    

real      df00020  = H_TT*T_0001_p2;
ds->df00020 += factor*roa_rob*df00020;

/* Third order derivatives */
real      df30000  = A_30*H_A+A_10_p3*H_AAA+H0_30*H_H+3.0*A_10_p2*(H0_10*H_HAA+H_TAA*T_1000)+3.0*A_10*(A_20*H_AA+H0_20*H_HA+2.0*H0_10*H_HTA*T_1000+H_TTA*T_1000_p2+H_TA*T_2000)+T_1000*(3.0*A_20*H_TA+3.0*H0_20*H_HT+H_TTT*T_1000_p2+3.0*H_TT*T_2000)+3.0*H0_10*(A_20*H_HA+H_HTT*T_1000_p2+H_HT*T_2000)+H_T*T_3000 + EC_30;
ds->df30000 += factor*(3.0*df20000 + roa_rob*df30000);

real      df21000  = A_21*H_A+H0_21*H_H+2.0*H0_10*A_11*H_HA+H0_01*A_20*H_HA+A_20*H_TA*T_0100+H0_20*H_HT*T_0100+A_10_p2*(A_01*H_AAA+H0_01*H_HAA+H_TAA*T_0100)+2.0*A_11*H_TA*T_1000+2.0*H0_11*H_HT*T_1000+2.0*H0_10*H_HTT*T_0100*T_1000+H0_01*H_HTT*T_1000_p2+H_TTT*T_0100*T_1000_p2+2.0*H0_10*H_HT*T_1100+2.0*H_TT*T_1000*T_1100+2.0*A_10*(A_11*H_AA+H0_11*H_HA+H0_10*(A_01*H_HAA+H_HTA*T_0100)+(A_01*H_TAA+H0_01*H_HTA+H_TTA*T_0100)*T_1000+H_TA*T_1100)+(H0_01*H_HT+H_TT*T_0100)*T_2000+A_01*(A_20*H_AA+H0_20*H_HA+2.0*H0_10*H_HTA*T_1000+H_TTA*T_1000_p2+H_TA*T_2000)+H_T*T_2100 + EC_21;
ds->df21000 += factor*(2*df11000 + df20000 + roa_rob*df21000);

real      df20100  = 2.0*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)*T_1010+T_0010*(A_20*H_TA+H0_20*H_HT+A_10*(A_10*H_TAA+2.0*H0_10*H_HTA)+2.0*(A_10*H_TTA+H0_10*H_HTT)*T_1000+H_TTT*T_1000_p2+H_TT*T_2000)+H_T*T_2010;
ds->df20100 += factor*(2.0*df10100 + roa_rob*df20100);

real      df20010  = 2.0*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)*T_1001+T_0001*(A_20*H_TA+H0_20*H_HT+A_10*(A_10*H_TAA+2.0*H0_10*H_HTA)+2.0*(A_10*H_TTA+H0_10*H_HTT)*T_1000+H_TTT*T_1000_p2+H_TT*T_2000)+H_T*T_2001;
ds->df20010 += factor*(2.0*df10010 + roa_rob*df20010);

real      df12000  = A_12*H_A+H0_12*H_H+A_11*(A_01*H_AA+2.0*H0_01*H_HA)+A_02*(A_10*H_AA+H0_10*H_HA+H_TA*T_1000)+T_0200*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+T_0100_p2*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+H0_02*(A_10*H_HA+H_HT*T_1000)+A_01*(A_11*H_AA+2.0*H0_11*H_HA+A_01*(A_10*H_AAA+H0_10*H_HAA+H_TAA*T_1000)+2.0*H0_01*(A_10*H_HAA+H_HTA*T_1000))+2.0*T_0100*(A_11*H_TA+H0_11*H_HT+A_01*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+H0_01*(A_10*H_HTA+H_HTT*T_1000))+2.0*(A_01*H_TA+H0_01*H_HT)*T_1100+2.0*H_TT*T_0100*T_1100+H_T*T_1200 + EC_12;
ds->df12000 += factor*(2.0*df11000 + df02000 + roa_rob*df12000);

real      df11100  = T_0110*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_1010+T_0010*(A_11*H_TA+H0_11*H_HT+A_01*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+T_0100*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+H0_01*(A_10*H_HTA+H_HTT*T_1000)+H_TT*T_1100)+H_T*T_1110;
ds->df11100 += factor*(df10100 + df01100 + roa_rob*df11100);

real      df11010  = T_0101*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_1001+T_0001*(A_11*H_TA+H0_11*H_HT+A_01*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+T_0100*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+H0_01*(A_10*H_HTA+H_HTT*T_1000)+H_TT*T_1100)+H_T*T_1101;
ds->df11010 += factor*(df10010 + df01010 + roa_rob*df11010);

real      df10200  = T_0010*(T_0010*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+2.0*H_TT*T_1010);
ds->df10200 += factor*(df00200 + roa_rob*df10200);

real      df10110  = T_0010*(T_0001*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+H_TT*T_1001)+H_TT*T_0001*T_1010;
ds->df10110 += factor*(df00110 + roa_rob*df10110);

real      df10020  = T_0001*(T_0001*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+2.0*H_TT*T_1001);
ds->df10020 += factor*(df00020 + roa_rob*df10020);

real      df03000  = A_03*H_A+A_01_p3*H_AAA+H0_03*H_H+3.0*A_01_p2*(H0_01*H_HAA+H_TAA*T_0100)+3.0*A_01*(A_02*H_AA+H0_02*H_HA+2.0*H0_01*H_HTA*T_0100+H_TTA*T_0100_p2+H_TA*T_0200)+T_0100*(3.0*A_02*H_TA+3.0*H0_02*H_HT+H_TTT*T_0100_p2+3.0*H_TT*T_0200)+3.0*H0_01*(A_02*H_HA+H_HTT*T_0100_p2+H_HT*T_0200)+H_T*T_0300 + EC_03;
ds->df03000 += factor*(3.0*df02000 + roa_rob*df03000);

real      df02100  = 2.0*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_0110+T_0010*(A_02*H_TA+H0_02*H_HT+A_01*(A_01*H_TAA+2.0*H0_01*H_HTA)+2.0*(A_01*H_TTA+H0_01*H_HTT)*T_0100+H_TTT*T_0100_p2+H_TT*T_0200)+H_T*T_0210;
ds->df02100 += factor*(2.0*df01100 + roa_rob*df02100);

real      df02010  = 2.0*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_0101+T_0001*(A_02*H_TA+H0_02*H_HT+A_01*(A_01*H_TAA+2.0*H0_01*H_HTA)+2.0*(A_01*H_TTA+H0_01*H_HTT)*T_0100+H_TTT*T_0100_p2+H_TT*T_0200)+H_T*T_0201;
ds->df02010 += factor*(2.0*df01010 + roa_rob*df02010);

real      df01200  = T_0010*(T_0010*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)+2.0*H_TT*T_0110);
ds->df01200 += factor*(df00200 + roa_rob*df01200);

real      df01110  = T_0010*(T_0001*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)+H_TT*T_0101)+H_TT*T_0001*T_0110;
ds->df01110 += factor*(df00110 + roa_rob*df01110);

real      df01020  = T_0001*(T_0001*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)+2.0*H_TT*T_0101);
ds->df01020 += factor*(df00020 + roa_rob*df01020);

real      df00300  = H_TTT*T_0010_p3;
ds->df00300 += factor*roa_rob*df00300;

real      df00210  = H_TTT*T_0001*T_0010_p2;
ds->df00210 += factor*roa_rob*df00210;

real      df00120  = H_TTT*T_0001_p2*T_0010;
ds->df00120 += factor*roa_rob*df00120;

real      df00030  = H_TTT*T_0001_p3;
ds->df00030 += factor*roa_rob*df00030;

/* Fourth order derivatives */

real      df40000  = A_40*H_A+3.0*A_20_p2*H_AA+A_10_p4*H_AAAA+H0_40*H_H+4.0*H0_10*A_30*H_HA+4.0*A_30*H_TA*T_1000+4.0*H0_30*H_HT*T_1000+6.0*H0_20*H_HTT*T_1000_p2+4.0*H0_10*H_HTTT*T_1000_p3+H_TTTT*T_1000_p4+4.0*A_10_p3*(H0_10*H_HAAA+H_TAAA*T_1000)+6.0*H0_20*H_HT*T_2000+12.0*H0_10*H_HTT*T_1000*T_2000+6.0*H_TTT*T_1000_p2*T_2000+3.0*H_TT*T_2000_p2+6.0*A_20*(H0_20*H_HA+A_10*(A_10*H_AAA+2.0*H0_10*H_HAA)+2.0*(A_10*H_TAA+H0_10*H_HTA)*T_1000+H_TTA*T_1000_p2+H_TA*T_2000)+6.0*A_10_p2*(H0_20*H_HAA+2.0*H0_10*H_HTAA*T_1000+H_TTAA*T_1000_p2+H_TAA*T_2000)+4.0*(H0_10*H_HT+H_TT*T_1000)*T_3000+4.0*A_10*(A_30*H_AA+H0_30*H_HA+T_1000*(3.0*H0_20*H_HTA+3.0*H0_10*H_HATT*T_1000+H_TTTA*T_1000_p2)+3.0*(H0_10*H_HTA+H_TTA*T_1000)*T_2000+H_TA*T_3000)+H_T*T_4000 + EC_40;
ds->df40000 += factor*(4*df30000 + roa_rob*df40000);

real      df31000  = A_31*H_A+3.0*A_10_p2*A_11*H_AAA+H0_31*H_H+A_30*(A_01*H_AA+H0_01*H_HA+H_TA*T_0100)+A_10_p3*(A_01*H_AAAA+H0_01*H_HAAA+H_TAAA*T_0100)+H0_30*(A_01*H_HA+H_HT*T_0100)+6.0*A_10*A_11*(H0_10*H_HAA+H_TAA*T_1000)+3.0*A_10_p2*(H0_11*H_HAA+H0_10*(A_01*H_HAAA+H_HTAA*T_0100)+(A_01*H_TAAA+H0_01*H_HTAA+H_TTAA*T_0100)*T_1000+H_TAA*T_1100)+3.0*A_11*(A_20*H_AA+H0_20*H_HA+2.0*H0_10*H_HTA*T_1000+H_TTA*T_1000_p2+H_TA*T_2000)+T_1100*(3.0*A_20*H_TA+3.0*H0_20*H_HT+H_TTT*T_1000_p2+3.0*H_TT*T_2000)+3.0*H0_11*(A_20*H_HA+H_HTT*T_1000_p2+H_HT*T_2000)+3.0*A_10*(A_21*H_AA+H0_21*H_HA+A_20*(A_01*H_AAA+H0_01*H_HAA+H_TAA*T_0100)+H0_20*(A_01*H_HAA+H_HTA*T_0100)+2.0*H0_11*H_HTA*T_1000+2.0*H0_10*(A_01*H_HTAA+H_HATT*T_0100)*T_1000+(A_01*H_TTAA+H0_01*H_HATT+H_TTTA*T_0100)*T_1000_p2+2.0*H0_10*H_HTA*T_1100+2.0*H_TTA*T_1000*T_1100+(A_01*H_TAA+H0_01*H_HTA+H_TTA*T_0100)*T_2000+H_TA*T_2100)+T_1000*(3.0*A_21*H_TA+3.0*H0_21*H_HT+3.0*A_20*(A_01*H_TAA+H0_01*H_HTA+H_TTA*T_0100)+3.0*H0_20*(A_01*H_HTA+H_HTT*T_0100)+(A_01*H_TTTA+H0_01*H_HTTT+H_TTTT*T_0100)*T_1000_p2+2.0*H_TTT*T_1000*T_1100+3.0*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)*T_2000+3.0*H_TT*T_2100)+3.0*H0_10*(A_21*H_HA+A_20*(A_01*H_HAA+H_HTA*T_0100)+(A_01*H_HATT+H_HTTT*T_0100)*T_1000_p2+2.0*H_HTT*T_1000*T_1100+(A_01*H_HTA+H_HTT*T_0100)*T_2000+H_HT*T_2100)+(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_3000+H_T*T_3100 + EC_31;
ds->df31000 += factor*(3.0*df21000 + df30000 + roa_rob*df31000);

real      df30100  = A_30*H_TA*T_0010+A_10_p3*H_TAAA*T_0010+H0_30*H_HT*T_0010+3.0*A_10_p2*(T_0010*(H0_10*H_HTAA+H_TTAA*T_1000)+H_TAA*T_1010)+T_1010*(3.0*A_20*H_TA+3.0*H0_20*H_HT+H_TTT*T_1000_p2+3.0*H_TT*T_2000)+3.0*A_10*(A_20*H_TAA*T_0010+H0_20*H_HTA*T_0010+2.0*H0_10*H_HATT*T_0010*T_1000+H_TTTA*T_0010*T_1000_p2+2.0*H0_10*H_HTA*T_1010+2.0*H_TTA*T_1000*T_1010+H_TTA*T_0010*T_2000+H_TA*T_2010)+T_1000*(3.0*A_20*H_TTA*T_0010+3.0*H0_20*H_HTT*T_0010+H_TTTT*T_0010*T_1000_p2+2.0*H_TTT*T_1000*T_1010+3.0*H_TTT*T_0010*T_2000+3.0*H_TT*T_2010)+3.0*H0_10*(A_20*H_HTA*T_0010+H_HTTT*T_0010*T_1000_p2+H_HTT*(2.0*T_1000*T_1010+T_0010*T_2000)+H_HT*T_2010)+H_TT*T_0010*T_3000+H_T*T_3010;
ds->df30100 += factor*(3.0*df20100 + roa_rob*df30100);

real      df30010  = A_30*H_TA*T_0001+A_10_p3*H_TAAA*T_0001+H0_30*H_HT*T_0001+3.0*A_10_p2*(T_0001*(H0_10*H_HTAA+H_TTAA*T_1000)+H_TAA*T_1001)+T_1001*(3.0*A_20*H_TA+3.0*H0_20*H_HT+H_TTT*T_1000_p2+3.0*H_TT*T_2000)+3.0*A_10*(A_20*H_TAA*T_0001+H0_20*H_HTA*T_0001+2.0*H0_10*H_HATT*T_0001*T_1000+H_TTTA*T_0001*T_1000_p2+2.0*H0_10*H_HTA*T_1001+2.0*H_TTA*T_1000*T_1001+H_TTA*T_0001*T_2000+H_TA*T_2001)+T_1000*(3.0*A_20*H_TTA*T_0001+3.0*H0_20*H_HTT*T_0001+H_TTTT*T_0001*T_1000_p2+2.0*H_TTT*T_1000*T_1001+3.0*H_TTT*T_0001*T_2000+3.0*H_TT*T_2001)+3.0*H0_10*(A_20*H_HTA*T_0001+H_HTTT*T_0001*T_1000_p2+H_HTT*(2.0*T_1000*T_1001+T_0001*T_2000)+H_HT*T_2001)+H_TT*T_0001*T_3000+H_T*T_3001;
ds->df30010 += factor*(3.0*df20010 + roa_rob*df30010);

real      df22000  = A_22*H_A+H0_22*H_H+A_21*(A_01*H_AA+2.0*H0_01*H_HA)+2.0*A_12*(A_10*H_AA+H0_10*H_HA+H_TA*T_1000)+2.0*H0_12*(A_10*H_HA+H_HT*T_1000)+2.0*A_11*(A_11*H_AA+2.0*H0_11*H_HA+A_01*(A_10*H_AAA+H0_10*H_HAA+H_TAA*T_1000)+2.0*H0_01*(A_10*H_HAA+H_HTA*T_1000))+4.0*T_0100*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)*T_1100+4.0*(A_11*H_TA+H0_11*H_HT+A_01*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+H0_01*(A_10*H_HTA+H_HTT*T_1000))*T_1100+2.0*H_TT*T_1100_p2+2.0*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)*T_1200+A_02*(A_20*H_AA+H0_20*H_HA+A_10*(A_10*H_AAA+2.0*H0_10*H_HAA)+2.0*(A_10*H_TAA+H0_10*H_HTA)*T_1000+H_TTA*T_1000_p2+H_TA*T_2000)+T_0200*(A_20*H_TA+H0_20*H_HT+A_10*(A_10*H_TAA+2.0*H0_10*H_HTA)+2.0*(A_10*H_TTA+H0_10*H_HTT)*T_1000+H_TTT*T_1000_p2+H_TT*T_2000)+T_0100_p2*(A_20*H_TTA+H0_20*H_HTT+A_10*(A_10*H_TTAA+2.0*H0_10*H_HATT)+2.0*(A_10*H_TTTA+H0_10*H_HTTT)*T_1000+H_TTTT*T_1000_p2+H_TTT*T_2000)+H0_02*(A_20*H_HA+A_10_p2*H_HAA+2.0*A_10*H_HTA*T_1000+H_HTT*T_1000_p2+H_HT*T_2000)+2.0*T_0100*(A_21*H_TA+H0_21*H_HT+2.0*H0_10*A_11*H_HTA+A_10_p2*(A_01*H_TAAA+H0_01*H_HTAA)+2.0*(A_11*H_TTA+H0_11*H_HTT)*T_1000+2.0*A_10*(A_11*H_TAA+H0_11*H_HTA+A_01*H0_10*H_HTAA+(A_01*H_TTAA+H0_01*H_HATT)*T_1000)+A_01*(A_20*H_TAA+H0_20*H_HTA+2.0*H0_10*H_HATT*T_1000+H_TTTA*T_1000_p2+H_TTA*T_2000)+H0_01*(A_20*H_HTA+H_HTTT*T_1000_p2+H_HTT*T_2000))+A_01*(A_21*H_AA+A_10_p2*(A_01*H_AAAA+2.0*H0_01*H_HAAA)+2.0*A_10*(A_11*H_AAA+2.0*H0_11*H_HAA+A_01*H0_10*H_HAAA+(A_01*H_TAAA+2.0*H0_01*H_HTAA)*T_1000)+A_01*(A_20*H_AAA+H0_20*H_HAA+2.0*H0_10*H_HTAA*T_1000+H_TTAA*T_1000_p2+H_TAA*T_2000)+2.0*(H0_21*H_HA+H0_10*A_11*H_HAA+(A_11*H_TAA+2.0*H0_11*H_HTA)*T_1000+H0_01*(A_20*H_HAA+H_HATT*T_1000_p2+H_HTA*T_2000)))+2.0*(A_01*H_TA+H0_01*H_HT)*T_2100+2.0*H_TT*T_0100*T_2100+H_T*T_2200 + EC_22;
ds->df22000 += factor*(2.0*df12000 + 2.0*df21000 + roa_rob*df22000);

real      df21100  = A_21*H_TA*T_0010+H0_21*H_HT*T_0010+2.0*H0_10*A_11*H_HTA*T_0010+H0_01*A_20*H_HTA*T_0010+A_20*H_TTA*T_0010*T_0100+H0_20*H_HTT*T_0010*T_0100+A_20*H_TA*T_0110+H0_20*H_HT*T_0110+A_10_p2*(T_0010*(A_01*H_TAAA+H0_01*H_HTAA+H_TTAA*T_0100)+H_TAA*T_0110)+2.0*A_11*H_TTA*T_0010*T_1000+2.0*H0_11*H_HTT*T_0010*T_1000+2.0*H0_10*H_HTTT*T_0010*T_0100*T_1000+2.0*H0_10*H_HTT*T_0110*T_1000+H0_01*H_HTTT*T_0010*T_1000_p2+H_TTTT*T_0010*T_0100*T_1000_p2+H_TTT*T_0110*T_1000_p2+2.0*A_11*H_TA*T_1010+2.0*H0_11*H_HT*T_1010+2.0*H0_10*H_HTT*T_0100*T_1010+2.0*H0_01*H_HTT*T_1000*T_1010+2.0*H_TTT*T_0100*T_1000*T_1010+2.0*H0_10*H_HTT*T_0010*T_1100+2.0*H_TTT*T_0010*T_1000*T_1100+2.0*H_TT*T_1010*T_1100+2.0*H0_10*H_HT*T_1110+2.0*H_TT*T_1000*T_1110+2.0*A_10*(A_11*H_TAA*T_0010+H0_11*H_HTA*T_0010+H0_10*(T_0010*(A_01*H_HTAA+H_HATT*T_0100)+H_HTA*T_0110)+(T_0010*(A_01*H_TTAA+H0_01*H_HATT+H_TTTA*T_0100)+H_TTA*T_0110)*T_1000+(A_01*H_TAA+H0_01*H_HTA+H_TTA*T_0100)*T_1010+H_TTA*T_0010*T_1100+H_TA*T_1110)+(T_0010*(H0_01*H_HTT+H_TTT*T_0100)+H_TT*T_0110)*T_2000+(H0_01*H_HT+H_TT*T_0100)*T_2010+A_01*(A_20*H_TAA*T_0010+H0_20*H_HTA*T_0010+2.0*H0_10*H_HATT*T_0010*T_1000+H_TTTA*T_0010*T_1000_p2+2.0*H0_10*H_HTA*T_1010+2.0*H_TTA*T_1000*T_1010+H_TTA*T_0010*T_2000+H_TA*T_2010)+H_TT*T_0010*T_2100+H_T*T_2110;
ds->df21100 +=factor*(2.0*df11100 + df20100 + roa_rob*df21100);

real      df21010  = A_21*H_TA*T_0001+H0_21*H_HT*T_0001+2.0*H0_10*A_11*H_HTA*T_0001+H0_01*A_20*H_HTA*T_0001+A_20*H_TTA*T_0001*T_0100+H0_20*H_HTT*T_0001*T_0100+A_20*H_TA*T_0101+H0_20*H_HT*T_0101+A_10_p2*(T_0001*(A_01*H_TAAA+H0_01*H_HTAA+H_TTAA*T_0100)+H_TAA*T_0101)+2.0*A_11*H_TTA*T_0001*T_1000+2.0*H0_11*H_HTT*T_0001*T_1000+2.0*H0_10*H_HTTT*T_0001*T_0100*T_1000+2.0*H0_10*H_HTT*T_0101*T_1000+H0_01*H_HTTT*T_0001*T_1000_p2+H_TTTT*T_0001*T_0100*T_1000_p2+H_TTT*T_0101*T_1000_p2+2.0*A_11*H_TA*T_1001+2.0*H0_11*H_HT*T_1001+2.0*H0_10*H_HTT*T_0100*T_1001+2.0*H0_01*H_HTT*T_1000*T_1001+2.0*H_TTT*T_0100*T_1000*T_1001+2.0*H0_10*H_HTT*T_0001*T_1100+2.0*H_TTT*T_0001*T_1000*T_1100+2.0*H_TT*T_1001*T_1100+2.0*H0_10*H_HT*T_1101+2.0*H_TT*T_1000*T_1101+2.0*A_10*(A_11*H_TAA*T_0001+H0_11*H_HTA*T_0001+H0_10*(T_0001*(A_01*H_HTAA+H_HATT*T_0100)+H_HTA*T_0101)+(T_0001*(A_01*H_TTAA+H0_01*H_HATT+H_TTTA*T_0100)+H_TTA*T_0101)*T_1000+(A_01*H_TAA+H0_01*H_HTA+H_TTA*T_0100)*T_1001+H_TTA*T_0001*T_1100+H_TA*T_1101)+(T_0001*(H0_01*H_HTT+H_TTT*T_0100)+H_TT*T_0101)*T_2000+(H0_01*H_HT+H_TT*T_0100)*T_2001+A_01*(A_20*H_TAA*T_0001+H0_20*H_HTA*T_0001+2.0*H0_10*H_HATT*T_0001*T_1000+H_TTTA*T_0001*T_1000_p2+2.0*H0_10*H_HTA*T_1001+2.0*H_TTA*T_1000*T_1001+H_TTA*T_0001*T_2000+H_TA*T_2001)+H_TT*T_0001*T_2100+H_T*T_2101;
ds->df21010 +=factor*(2.0*df11010 + df20010 + roa_rob*df21010);

real      df20200  = A_20*H_TTA*T_0010_p2+A_10_p2*H_TTAA*T_0010_p2+H0_20*H_HTT*T_0010_p2+2.0*H0_10*H_HTTT*T_0010_p2*T_1000+H_TTTT*T_0010_p2*T_1000_p2+4.0*H0_10*H_HTT*T_0010*T_1010+4.0*H_TTT*T_0010*T_1000*T_1010+2.0*H_TT*T_1010_p2+2.0*A_10*T_0010*(T_0010*(H0_10*H_HATT+H_TTTA*T_1000)+2.0*H_TTA*T_1010)+H_TTT*T_0010_p2*T_2000+2.0*H_TT*T_0010*T_2010;
ds->df20200 += factor*(2.0*df10200 + roa_rob*df20200);

real      df20110  = 2.0*(T_0001*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+H_TT*T_1001)*T_1010+T_0010*(A_20*H_TTA*T_0001+H0_20*H_HTT*T_0001+A_10*(A_10*H_TTAA+2.0*H0_10*H_HATT)*T_0001+2.0*(A_10*H_TTTA+H0_10*H_HTTT)*T_0001*T_1000+H_TTTT*T_0001*T_1000_p2+2.0*(A_10*H_TTA+H0_10*H_HTT)*T_1001+2.0*H_TTT*T_1000*T_1001+H_TTT*T_0001*T_2000+H_TT*T_2001)+H_TT*T_0001*T_2010;
ds->df20110 += factor*(2.0*df10110 + roa_rob*df20110);

real      df20020  = A_20*H_TTA*T_0001_p2+A_10_p2*H_TTAA*T_0001_p2+H0_20*H_HTT*T_0001_p2+2.0*H0_10*H_HTTT*T_0001_p2*T_1000+H_TTTT*T_0001_p2*T_1000_p2+4.0*H0_10*H_HTT*T_0001*T_1001+4.0*H_TTT*T_0001*T_1000*T_1001+2.0*H_TT*T_1001_p2+2.0*A_10*T_0001*(T_0001*(H0_10*H_HATT+H_TTTA*T_1000)+2.0*H_TTA*T_1001)+H_TTT*T_0001_p2*T_2000+2.0*H_TT*T_0001*T_2001;
ds->df20020 += factor*(2.0*df10020 + roa_rob*df20020);

real      df13000  = A_13*H_A+3.0*A_01_p2*A_11*H_AAA+H0_13*H_H+6.0*A_01*A_11*(H0_01*H_HAA+H_TAA*T_0100)+3.0*A_11*(A_02*H_AA+H0_02*H_HA+2.0*H0_01*H_HTA*T_0100+H_TTA*T_0100_p2+H_TA*T_0200)+3.0*H0_11*(A_02*H_HA+H_HTT*T_0100_p2+H_HT*T_0200)+A_03*(A_10*H_AA+H0_10*H_HA+H_TA*T_1000)+A_01_p3*(A_10*H_AAAA+H0_10*H_HAAA+H_TAAA*T_1000)+T_0300*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+H0_03*(A_10*H_HA+H_HT*T_1000)+(3.0*A_02*H_TA+3.0*H0_02*H_HT+H_TTT*T_0100_p2+3.0*H_TT*T_0200)*T_1100+3.0*A_01_p2*(H0_11*H_HAA+T_0100*(A_10*H_TAAA+H0_10*H_HTAA+H_TTAA*T_1000)+H0_01*(A_10*H_HAAA+H_HTAA*T_1000)+H_TAA*T_1100)+3.0*A_01*(A_12*H_AA+H0_12*H_HA+2.0*H0_11*H_HTA*T_0100+A_02*(A_10*H_AAA+H0_10*H_HAA+H_TAA*T_1000)+T_0200*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+T_0100_p2*(A_10*H_TTAA+H0_10*H_HATT+H_TTTA*T_1000)+H0_02*(A_10*H_HAA+H_HTA*T_1000)+2.0*H0_01*T_0100*(A_10*H_HTAA+H_HATT*T_1000)+2.0*H0_01*H_HTA*T_1100+2.0*H_TTA*T_0100*T_1100+H_TA*T_1200)+T_0100*(3.0*A_12*H_TA+3.0*H0_12*H_HT+3.0*A_02*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+3.0*T_0200*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+T_0100_p2*(A_10*H_TTTA+H0_10*H_HTTT+H_TTTT*T_1000)+3.0*H0_02*(A_10*H_HTA+H_HTT*T_1000)+2.0*H_TTT*T_0100*T_1100+3.0*H_TT*T_1200)+3.0*H0_01*(A_12*H_HA+A_02*(A_10*H_HAA+H_HTA*T_1000)+T_0200*(A_10*H_HTA+H_HTT*T_1000)+T_0100_p2*(A_10*H_HATT+H_HTTT*T_1000)+2.0*H_HTT*T_0100*T_1100+H_HT*T_1200)+H_T*T_1300 + EC_13;
ds->df13000 += factor*(3.0*df12000 + df03000 + roa_rob*df13000);

real      df12100  = T_0210*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+(A_02*H_TA+H0_02*H_HT+A_01*(A_01*H_TAA+2.0*H0_01*H_HTA)+2.0*(A_01*H_TTA+H0_01*H_HTT)*T_0100+H_TTT*T_0100_p2+H_TT*T_0200)*T_1010+2.0*T_0110*(A_11*H_TA+H0_11*H_HT+A_01*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+T_0100*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+H0_01*(A_10*H_HTA+H_HTT*T_1000)+H_TT*T_1100)+2.0*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_1110+T_0010*(A_12*H_TA+H0_12*H_HT+A_11*(A_01*H_TAA+2.0*H0_01*H_HTA)+A_02*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+T_0200*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+T_0100_p2*(A_10*H_TTTA+H0_10*H_HTTT+H_TTTT*T_1000)+H0_02*(A_10*H_HTA+H_HTT*T_1000)+A_01*(A_11*H_TAA+2.0*H0_11*H_HTA+A_01*(A_10*H_TAAA+H0_10*H_HTAA+H_TTAA*T_1000)+2.0*H0_01*(A_10*H_HTAA+H_HATT*T_1000))+2.0*T_0100*(A_11*H_TTA+H0_11*H_HTT+A_01*(A_10*H_TTAA+H0_10*H_HATT+H_TTTA*T_1000)+H0_01*(A_10*H_HATT+H_HTTT*T_1000))+2.0*(A_01*H_TTA+H0_01*H_HTT)*T_1100+2.0*H_TTT*T_0100*T_1100+H_TT*T_1200)+H_T*T_1210;
ds->df12100 += factor*(2.0*df11100 + df02100 + roa_rob*df12100);


real      df12010  = T_0201*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+(A_02*H_TA+H0_02*H_HT+A_01*(A_01*H_TAA+2.0*H0_01*H_HTA)+2.0*(A_01*H_TTA+H0_01*H_HTT)*T_0100+H_TTT*T_0100_p2+H_TT*T_0200)*T_1001+2.0*T_0101*(A_11*H_TA+H0_11*H_HT+A_01*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+T_0100*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+H0_01*(A_10*H_HTA+H_HTT*T_1000)+H_TT*T_1100)+2.0*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_1101+T_0001*(A_12*H_TA+H0_12*H_HT+A_11*(A_01*H_TAA+2.0*H0_01*H_HTA)+A_02*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+T_0200*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+T_0100_p2*(A_10*H_TTTA+H0_10*H_HTTT+H_TTTT*T_1000)+H0_02*(A_10*H_HTA+H_HTT*T_1000)+A_01*(A_11*H_TAA+2.0*H0_11*H_HTA+A_01*(A_10*H_TAAA+H0_10*H_HTAA+H_TTAA*T_1000)+2.0*H0_01*(A_10*H_HTAA+H_HATT*T_1000))+2.0*T_0100*(A_11*H_TTA+H0_11*H_HTT+A_01*(A_10*H_TTAA+H0_10*H_HATT+H_TTTA*T_1000)+H0_01*(A_10*H_HATT+H_HTTT*T_1000))+2.0*(A_01*H_TTA+H0_01*H_HTT)*T_1100+2.0*H_TTT*T_0100*T_1100+H_TT*T_1200)+H_T*T_1201;
ds->df12010 += factor*(2.0*df11010 + df02010 + roa_rob*df12010);

real      df11200  = (T_0010*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)+2.0*H_TT*T_0110)*T_1010+T_0010*(2.0*T_0110*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)*T_1010+T_0010*(A_11*H_TTA+H0_11*H_HTT+A_01*(A_10*H_TTAA+H0_10*H_HATT+H_TTTA*T_1000)+T_0100*(A_10*H_TTTA+H0_10*H_HTTT+H_TTTT*T_1000)+H0_01*(A_10*H_HATT+H_HTTT*T_1000)+H_TTT*T_1100)+2.0*H_TT*T_1110);
ds->df11200 += factor*(df10200 + df01200 + roa_rob*df11200);

real      df11110  = T_0001*T_0110*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+H_TT*T_0110*T_1001+(T_0001*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)+H_TT*T_0101)*T_1010+T_0010*(T_0101*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)*T_1001+T_0001*(A_11*H_TTA+H0_11*H_HTT+A_01*(A_10*H_TTAA+H0_10*H_HATT+H_TTTA*T_1000)+T_0100*(A_10*H_TTTA+H0_10*H_HTTT+H_TTTT*T_1000)+H0_01*(A_10*H_HATT+H_HTTT*T_1000)+H_TTT*T_1100)+H_TT*T_1101)+H_TT*T_0001*T_1110;
ds->df11110 += factor*(df10110 + df01110 + roa_rob*df11110);

real      df11020  = (T_0001*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)+2.0*H_TT*T_0101)*T_1001+T_0001*(2.0*T_0101*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)*T_1001+T_0001*(A_11*H_TTA+H0_11*H_HTT+A_01*(A_10*H_TTAA+H0_10*H_HATT+H_TTTA*T_1000)+T_0100*(A_10*H_TTTA+H0_10*H_HTTT+H_TTTT*T_1000)+H0_01*(A_10*H_HATT+H_HTTT*T_1000)+H_TTT*T_1100)+2.0*H_TT*T_1101);
ds->df11020 += factor*(df10020 + df01020 + roa_rob*df11020);

real      df10300  = T_0010_p2*(T_0010*(A_10*H_TTTA+H0_10*H_HTTT+H_TTTT*T_1000)+3.0*H_TTT*T_1010);
ds->df10300 += factor*(df00300 + roa_rob*df10300);

real      df10210  = T_0010*(T_0010*(T_0001*(A_10*H_TTTA+H0_10*H_HTTT+H_TTTT*T_1000)+H_TTT*T_1001)+2.0*H_TTT*T_0001*T_1010);
ds->df10210 += factor*(df00210 + roa_rob*df10210);

real      df10120  = T_0001*(T_0010*(T_0001*(A_10*H_TTTA+H0_10*H_HTTT+H_TTTT*T_1000)+2.0*H_TTT*T_1001)+H_TTT*T_0001*T_1010);
ds->df10120 += factor*(df00120 + roa_rob*df10120);

real      df10030  = T_0001_p2*(T_0001*(A_10*H_TTTA+H0_10*H_HTTT+H_TTTT*T_1000)+3.0*H_TTT*T_1001);
ds->df10030 += factor*(df00030 + roa_rob*df10030);

real      df04000  = A_04*H_A+3.0*A_02_p2*H_AA+A_01_p4*H_AAAA+H0_04*H_H+4.0*H0_01*A_03*H_HA+4.0*A_03*H_TA*T_0100+4.0*H0_03*H_HT*T_0100+6.0*H0_02*H_HTT*T_0100_p2+4.0*H0_01*H_HTTT*T_0100_p3+H_TTTT*T_0100_p4+4.0*A_01_p3*(H0_01*H_HAAA+H_TAAA*T_0100)+6.0*H0_02*H_HT*T_0200+12.0*H0_01*H_HTT*T_0100*T_0200+6.0*H_TTT*T_0100_p2*T_0200+3.0*H_TT*T_0200_p2+6.0*A_02*(H0_02*H_HA+A_01*(A_01*H_AAA+2.0*H0_01*H_HAA)+2.0*(A_01*H_TAA+H0_01*H_HTA)*T_0100+H_TTA*T_0100_p2+H_TA*T_0200)+6.0*A_01_p2*(H0_02*H_HAA+2.0*H0_01*H_HTAA*T_0100+H_TTAA*T_0100_p2+H_TAA*T_0200)+4.0*(H0_01*H_HT+H_TT*T_0100)*T_0300+4.0*A_01*(A_03*H_AA+H0_03*H_HA+T_0100*(3.0*H0_02*H_HTA+3.0*H0_01*H_HATT*T_0100+H_TTTA*T_0100_p2)+3.0*(H0_01*H_HTA+H_TTA*T_0100)*T_0200+H_TA*T_0300)+H_T*T_0400 + EC_04;
ds->df04000 += factor*(4.0*df03000 + roa_rob*df04000);

real      df03100  = A_03*H_TA*T_0010+A_01_p3*H_TAAA*T_0010+H0_03*H_HT*T_0010+3.0*A_01_p2*(T_0010*(H0_01*H_HTAA+H_TTAA*T_0100)+H_TAA*T_0110)+T_0110*(3.0*A_02*H_TA+3.0*H0_02*H_HT+H_TTT*T_0100_p2+3.0*H_TT*T_0200)+3.0*A_01*(A_02*H_TAA*T_0010+H0_02*H_HTA*T_0010+2.0*H0_01*H_HATT*T_0010*T_0100+H_TTTA*T_0010*T_0100_p2+2.0*H0_01*H_HTA*T_0110+2.0*H_TTA*T_0100*T_0110+H_TTA*T_0010*T_0200+H_TA*T_0210)+T_0100*(3.0*A_02*H_TTA*T_0010+3.0*H0_02*H_HTT*T_0010+H_TTTT*T_0010*T_0100_p2+2.0*H_TTT*T_0100*T_0110+3.0*H_TTT*T_0010*T_0200+3.0*H_TT*T_0210)+3.0*H0_01*(A_02*H_HTA*T_0010+H_HTTT*T_0010*T_0100_p2+H_HTT*(2.0*T_0100*T_0110+T_0010*T_0200)+H_HT*T_0210)+H_TT*T_0010*T_0300+H_T*T_0310;
ds->df03100 += factor*(3.0*df02100 + roa_rob*df03100);

real      df03010  = A_03*H_TA*T_0001+A_01_p3*H_TAAA*T_0001+H0_03*H_HT*T_0001+3.0*A_01_p2*(T_0001*(H0_01*H_HTAA+H_TTAA*T_0100)+H_TAA*T_0101)+T_0101*(3.0*A_02*H_TA+3.0*H0_02*H_HT+H_TTT*T_0100_p2+3.0*H_TT*T_0200)+3.0*A_01*(A_02*H_TAA*T_0001+H0_02*H_HTA*T_0001+2.0*H0_01*H_HATT*T_0001*T_0100+H_TTTA*T_0001*T_0100_p2+2.0*H0_01*H_HTA*T_0101+2.0*H_TTA*T_0100*T_0101+H_TTA*T_0001*T_0200+H_TA*T_0201)+T_0100*(3.0*A_02*H_TTA*T_0001+3.0*H0_02*H_HTT*T_0001+H_TTTT*T_0001*T_0100_p2+2.0*H_TTT*T_0100*T_0101+3.0*H_TTT*T_0001*T_0200+3.0*H_TT*T_0201)+3.0*H0_01*(A_02*H_HTA*T_0001+H_HTTT*T_0001*T_0100_p2+H_HTT*(2.0*T_0100*T_0101+T_0001*T_0200)+H_HT*T_0201)+H_TT*T_0001*T_0300+H_T*T_0301;
ds->df03010 += factor*(3.0*df02010 + roa_rob*df03010);

real      df02200  = A_02*H_TTA*T_0010_p2+A_01_p2*H_TTAA*T_0010_p2+H0_02*H_HTT*T_0010_p2+2.0*H0_01*H_HTTT*T_0010_p2*T_0100+H_TTTT*T_0010_p2*T_0100_p2+4.0*H0_01*H_HTT*T_0010*T_0110+4.0*H_TTT*T_0010*T_0100*T_0110+2.0*H_TT*T_0110_p2+2.0*A_01*T_0010*(T_0010*(H0_01*H_HATT+H_TTTA*T_0100)+2.0*H_TTA*T_0110)+H_TTT*T_0010_p2*T_0200+2.0*H_TT*T_0010*T_0210;
ds->df02200 += factor*(2.0*df01200 + roa_rob*df02200);

real      df02110  = 2.0*(T_0001*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)+H_TT*T_0101)*T_0110+T_0010*(A_02*H_TTA*T_0001+H0_02*H_HTT*T_0001+A_01*(A_01*H_TTAA+2.0*H0_01*H_HATT)*T_0001+2.0*(A_01*H_TTTA+H0_01*H_HTTT)*T_0001*T_0100+H_TTTT*T_0001*T_0100_p2+2.0*(A_01*H_TTA+H0_01*H_HTT)*T_0101+2.0*H_TTT*T_0100*T_0101+H_TTT*T_0001*T_0200+H_TT*T_0201)+H_TT*T_0001*T_0210;
ds->df02110 += factor*(2.0*df01110 + roa_rob*df02110);

real      df02020  = A_02*H_TTA*T_0001_p2+A_01_p2*H_TTAA*T_0001_p2+H0_02*H_HTT*T_0001_p2+2.0*H0_01*H_HTTT*T_0001_p2*T_0100+H_TTTT*T_0001_p2*T_0100_p2+4.0*H0_01*H_HTT*T_0001*T_0101+4.0*H_TTT*T_0001*T_0100*T_0101+2.0*H_TT*T_0101_p2+2.0*A_01*T_0001*(T_0001*(H0_01*H_HATT+H_TTTA*T_0100)+2.0*H_TTA*T_0101)+H_TTT*T_0001_p2*T_0200+2.0*H_TT*T_0001*T_0201;
ds->df02020 += factor*(2.0*df01020 + roa_rob*df02020);

real      df01300  = T_0010_p2*(T_0010*(A_01*H_TTTA+H0_01*H_HTTT+H_TTTT*T_0100)+3.0*H_TTT*T_0110);
ds->df01300 += factor*(df00300 + roa_rob*df01300);

real      df01210  = T_0010*(T_0010*(T_0001*(A_01*H_TTTA+H0_01*H_HTTT+H_TTTT*T_0100)+H_TTT*T_0101)+2.0*H_TTT*T_0001*T_0110);
ds->df01210 += factor*(df00210 + roa_rob*df01210);

real      df01120  = T_0001*(T_0010*(T_0001*(A_01*H_TTTA+H0_01*H_HTTT+H_TTTT*T_0100)+2.0*H_TTT*T_0101)+H_TTT*T_0001*T_0110);
ds->df01120 += factor*(df00120 + roa_rob*df01120);

real      df01030  = T_0001_p2*(T_0001*(A_01*H_TTTA+H0_01*H_HTTT+H_TTTT*T_0100)+3.0*H_TTT*T_0101);
ds->df01030 += factor*(df00030 + roa_rob*df01030);

real      df00400  = H_TTTT*T_0010_p4;
ds->df00400 += factor*roa_rob*df00400;

real      df00310  = H_TTTT*T_0001*T_0010_p3;
ds->df00310 += factor*roa_rob*df00310;

real      df00220  = H_TTTT*T_0001_p2*T_0010_p2;
ds->df00220 += factor*roa_rob*df00220;

real      df00130  = H_TTTT*T_0001_p3*T_0010;
ds->df00130 += factor*roa_rob*df00130;

real      df00040  = H_TTTT*T_0001_p4;
ds->df00040 += factor*roa_rob*df00040;

/* End (Finally!) */
}
#endif

