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

#include <math.h>
#include <stdio.h>

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int pbe_isgga(void) { return 0; }
static int pbe_read(const char* conf_line);
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
static int
pbe_read(const char* conf_line)
{
    dft_set_hf_weight(0);
    return 1;
}

static real
pbe_energy(const DftDensProp* dp)
{
    return (
        pbe_xenergy(dp) + 
        pbe_cenergy(dp) 
        );
}

static void
pbe_first(FirstFuncDrv *ds, real factor, const DftDensProp* dp)
{
    pbe_xd1(ds, factor, dp);
    pbe_cd1(ds, factor, dp);
}
static void
pbe_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp)
{
    pbe_xd2(ds, factor, dp);
    pbe_cd2(ds, factor, dp);
}

static void
pbe_third(ThirdFuncDrv *ds, real factor, const DftDensProp* dp)
{
    pbe_xd3(ds, factor, dp);
    pbe_cd3(ds, factor, dp);
}

#ifdef FOURTH_ORDER_DERIVATIVES
static void
pbe_fourth(FourthFuncDrv *ds, real factor, const DftDensProp *dp)
{
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

/* PBE functional exchange energy */
static real
pbe_xenergy(const DftDensProp* dp)
{
    real E;

    /* E=(roa*pbe_ro2energy(roa,groa) +
     * rob*pbe_ro2energy(rob,grob))/(roa+rob); */
    E= (dp->rhoa*pbe_ro2energy(&dp->rhoa, &dp->grada) + 
        dp->rhob*pbe_ro2energy(&dp->rhob, &dp->gradb))/(dp->rhoa+dp->rhob);

    return E;
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
    real  KF;
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
    KF=COR_nKF*roa_rob_p1f3;
    KS=COR_nKS*roa_rob_p1f6;
    MZETA=pow((1.0-ZETA),(2.0/3.0));
    PZETA=pow((1.0+ZETA),(2.0/3.0));
    FI=(PZETA + MZETA)/2.0;
    T=(groa+grob)/((roa_rob)*(2.0*KS*FI));


    FI_p3=FI*FI*FI;
    T_p2=T*T;
    T_p4=T_p2*T_p2;

    EC  = PW92Functional.func(dp);

    H0  = COR_gamma*FI_p3;

    A   = COR_delta*( 1.0/(exp(-EC/H0) -1.0)  ) ;
    A_p2= A*A;

    H1=1.0 + COR_delta*T_p2*((1.0 + A*T_p2)/(1.0 + A*T_p2 + A_p2*T_p4));

    H=H0*log( H1 );
    E=EC+H;

    return(E);
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

    EXunif=-3.0*KF/(4.0*M_PI);
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
    TwoVarDrv  dera;
    TwoVarDrv  derb;
    real  Exa;
    real  Exb;
    real  Exa_01;
    real  Exa_10;
    real  Exb_01;
    real  Exb_10;
    real  groa;
    real  grob;
    real  roa;
    real  roa_rob;
    real  roa_rob_p2;
    real  rob;

  
    roa = dp->rhoa;
    rob = dp->rhob;

    groa = dp->grada;
    grob = dp->gradb;
	
    Exa=pbe_ro2energy(&roa, &groa);
    Exb=pbe_ro2energy(&rob, &grob);

    pbe_ro2d(&dera, &roa, &groa);
    pbe_ro2d(&derb, &rob, &grob);

    Exa_10=dera.df10;
    Exa_01=dera.df01;

    Exb_10=derb.df10;
    Exb_01=derb.df01;

/* Powers of total density */
    roa_rob=roa+rob;
    roa_rob_p2=roa_rob*roa_rob;

/* Derivatives */

    ds->df1000 += factor*(
        (rob*Exa - rob*Exb + roa*roa_rob*Exa_10)/roa_rob_p2
        );

    ds->df0100 += factor*(
        (-roa*Exa + roa*Exb + rob*roa_rob*Exb_10)/roa_rob_p2
        );

    ds->df0010 += factor*(
        (roa*Exa_01)/roa_rob
        );

    ds->df0001 += factor*(
        (rob*Exb_01)/roa_rob
        );
}

static void
pbe_xd2(SecondFuncDrv *ds, real factor, const DftDensProp *dp)
{
    TwoVarDrv  dera;
    TwoVarDrv  derb;
    real  Exa;
    real  Exb;
    real  Exa_01;
    real  Exa_02;
    real  Exa_10;
    real  Exa_11;
    real  Exa_20;
    real  Exb_01;
    real  Exb_02;
    real  Exb_10;
    real  Exb_11;
    real  Exb_20;
    real  groa;
    real  grob;
    real  roa;
    real  roa_rob;
    real  roa_rob_p2;
    real  roa_rob_p3;
    real  rob;

  
    roa = dp->rhoa;
    rob = dp->rhob;

    groa = dp->grada;
    grob = dp->gradb;
	
    Exa=pbe_ro2energy(&roa, &groa);
    Exb=pbe_ro2energy(&rob, &grob);

    pbe_ro2d(&dera, &roa, &groa);
    pbe_ro2d(&derb, &rob, &grob);

    Exa_10=dera.df10;
    Exa_01=dera.df01;
    Exa_20=dera.df20;
    Exa_11=dera.df11;
    Exa_02=dera.df02;

    Exb_10=derb.df10;
    Exb_01=derb.df01;
    Exb_20=derb.df20;
    Exb_11=derb.df11;
    Exb_02=derb.df02;

/* Powers of total density */
    roa_rob=roa+rob;
    roa_rob_p2=roa_rob*roa_rob;
    roa_rob_p3=roa_rob_p2*roa_rob;

/* Derivatives */

    ds->df1000 += factor*(
        (rob*Exa - rob*Exb + roa*roa_rob*Exa_10)/roa_rob_p2
        );

    ds->df0100 += factor*(
        (-roa*Exa + roa*Exb + rob*roa_rob*Exb_10)/roa_rob_p2
        );

    ds->df0010 += factor*(
        (roa*Exa_01)/roa_rob
        );

    ds->df0001 += factor*(
        (rob*Exb_01)/roa_rob
        );

/* Second order derivatives */

    ds->df2000 += factor*(
        (1.0/roa_rob_p3)*(-2.0*rob*Exa + 2.0*rob*Exb + 
                          roa_rob*(2.0*rob*Exa_10 + roa*roa_rob*Exa_20))
        );

    ds->df1100 += factor*(
        (1.0/roa_rob_p3)*((roa-rob)*Exa + (rob - roa)*Exb - 
                          roa_rob*(roa*Exa_10 + rob*Exb_10))
        );

    ds->df1010 += factor*(
        (rob*Exa_01 + roa*roa_rob*Exa_11)/roa_rob_p2
        );

    ds->df1001 += factor*(
        -rob*Exb_01/roa_rob_p2
        );

    ds->df0200 += factor*(
        (1.0/roa_rob_p3)*(2.0*roa*Exa - 2.0*roa*Exb + 
                          roa_rob*(2.0*roa*Exb_10 + rob*roa_rob*Exb_20))
        );

    ds->df0110 += factor*(
        -roa*Exa_01/roa_rob_p2
        );

    ds->df0101 += factor*(
        (roa*Exb_01 + rob*roa_rob*Exb_11)/roa_rob_p2
        );

    ds->df0020 += factor*(
        roa*Exa_02/roa_rob
        );

/* ds->df00110=0.0 */

    ds->df0002 += factor*(
        rob*Exb_02/roa_rob
        );
}

static void
pbe_xd3(ThirdFuncDrv *ds, real factor, const DftDensProp *dp)
{
    TwoVarDrv  dera;
    TwoVarDrv  derb;
    real  Exa;
    real  Exb;
    real  Exa_01;
    real  Exa_02;
    real  Exa_03;
    real  Exa_10;
    real  Exa_11;
    real  Exa_12;
    real  Exa_20;
    real  Exa_21;
    real  Exa_30;
    real  Exb_01;
    real  Exb_02;
    real  Exb_03;
    real  Exb_10;
    real  Exb_11;
    real  Exb_12;
    real  Exb_20;
    real  Exb_21;
    real  Exb_30;
    real  groa;
    real  grob;
    real  roa;
    real  roa_rob;
    real  roa_rob_p2;
    real  roa_rob_p3;
    real  roa_rob_p4;
    real  rob;

  
    roa = dp->rhoa;
    rob = dp->rhob;

    groa = dp->grada;
    grob = dp->gradb;
	
    Exa=pbe_ro2energy(&roa, &groa);
    Exb=pbe_ro2energy(&rob, &grob);

    pbe_ro2d(&dera, &roa, &groa);
    pbe_ro2d(&derb, &rob, &grob);

    Exa_10=dera.df10;
    Exa_01=dera.df01;
    Exa_20=dera.df20;
    Exa_11=dera.df11;
    Exa_02=dera.df02;
    Exa_30=dera.df30;
    Exa_21=dera.df21;
    Exa_12=dera.df12;
    Exa_03=dera.df03;

    Exb_10=derb.df10;
    Exb_01=derb.df01;
    Exb_20=derb.df20;
    Exb_11=derb.df11;
    Exb_02=derb.df02;
    Exb_30=derb.df30;
    Exb_21=derb.df21;
    Exb_12=derb.df12;
    Exb_03=derb.df03;

/* Powers of total density */
    roa_rob=roa+rob;
    roa_rob_p2=roa_rob*roa_rob;
    roa_rob_p3=roa_rob_p2*roa_rob;
    roa_rob_p4=roa_rob_p3*roa_rob;

/* Derivatives */

    ds->df1000 += factor*(
        (rob*Exa - rob*Exb + roa*roa_rob*Exa_10)/roa_rob_p2
        );

    ds->df0100 += factor*(
        (-roa*Exa + roa*Exb + rob*roa_rob*Exb_10)/roa_rob_p2
        );

    ds->df0010 += factor*(
        (roa*Exa_01)/roa_rob
        );

    ds->df0001 += factor*(
        (rob*Exb_01)/roa_rob
        );

/* Second order derivatives */

    ds->df2000 += factor*(
        (1.0/roa_rob_p3)*(-2.0*rob*Exa + 2.0*rob*Exb + 
                          roa_rob*(2.0*rob*Exa_10 + roa*roa_rob*Exa_20))
        );

    ds->df1100 += factor*(
        (1.0/roa_rob_p3)*((roa-rob)*Exa + (rob - roa)*Exb - 
                          roa_rob*(roa*Exa_10 + rob*Exb_10))
        );

    ds->df1010 += factor*(
        (rob*Exa_01 + roa*roa_rob*Exa_11)/roa_rob_p2
        );

    ds->df1001 += factor*(
        -rob*Exb_01/roa_rob_p2
        );

    ds->df0200 += factor*(
        (1.0/roa_rob_p3)*(2.0*roa*Exa - 2.0*roa*Exb + 
                          roa_rob*(2.0*roa*Exb_10 + rob*roa_rob*Exb_20))
        );

    ds->df0110 += factor*(
        -roa*Exa_01/roa_rob_p2
        );

    ds->df0101 += factor*(
        (roa*Exb_01 + rob*roa_rob*Exb_11)/roa_rob_p2
        );

    ds->df0020 += factor*(
        roa*Exa_02/roa_rob
        );

/* ds->df00110=0.0 */

    ds->df0002 += factor*(
        rob*Exb_02/roa_rob
        );

/* Third order derivatives */

    ds->df3000 += factor*(
        (1.0/roa_rob_p4)*(6.0*rob*Exa - 6.0*rob*Exb + 
                          roa_rob*(-6.0*rob*Exa_10 + 
                                   roa_rob*(3.0*rob*Exa_20 + 
                                            roa*roa_rob*Exa_30)))
        );

    ds->df2100 += factor*(
        (1.0/roa_rob_p4)*(-2.0*(roa-2.0*rob)*Exa + 
                          2.0*(roa-2.0*rob)*Exb - 
                          roa_rob*(-2.0*(roa-rob)*Exa_10 - 
                                   2.0*rob*Exb_10 + roa*roa_rob*Exa_20))
        );

    ds->df2010 += factor*(
        (1.0/roa_rob_p3)*(-2.0*rob*Exa_01 + roa_rob*(2.0*rob*Exa_11 + 
                                                     roa*roa_rob*Exa_21))
        );

    ds->df2001 += factor*(
        2.0*rob*Exb_01/roa_rob_p3
        );

    ds->df1200 += factor*(
        (1.0/roa_rob_p4)*(-2.0*(2.0*roa-rob)*Exa + 
                          2.0*(2.0*roa-rob)*Exb - 
                          roa_rob*(-2.0*roa*Exa_10 + 
                                   2.0*(roa-rob)*Exb_10 + rob*roa_rob*Exb_20))
        );


    ds->df1110 += factor*(
        ((roa-rob)*Exa_01 - roa*roa_rob*Exa_11)/roa_rob_p3
        );

    ds->df1101 += factor*(
        ((rob-roa)*Exb_01 - rob*roa_rob*Exb_11)/roa_rob_p3
        );

    ds->df1020 += factor*(
        (rob*Exa_02 + roa*roa_rob*Exa_12)/roa_rob_p2
        );

/* ds->df10110=0.0 */

    ds->df1002 += factor*(
        -rob*Exb_02/roa_rob_p2
        );

    ds->df0300 += factor*(
        (1.0/roa_rob_p4)*(-6.0*roa*Exa + 6.0*roa*Exb + 
                          roa_rob*(-6.0*roa*Exb_10 + 
                                   roa_rob*(3.0*roa*Exb_20 + 
                                            rob*roa_rob*Exb_30)))
        );


    ds->df0210 += factor*(
        2.0*roa*Exa_01/roa_rob_p3
        );

    ds->df0201 += factor*(
        (1.0/roa_rob_p3)*(-2.0*roa*Exb_01 + roa_rob*(2.0*roa*Exb_11 + 
                                                     rob*roa_rob*Exb_21))
        );

    ds->df0120 += factor*(
        -roa*Exa_02/roa_rob_p2
        );

/* ds->df01110=0.0 */

    ds->df0102 += factor*(
        (roa*Exb_02 + rob*roa_rob*Exb_12)/roa_rob_p2
        );

    ds->df0030 += factor*(
        roa*Exa_03/roa_rob
        );

/* ds->df00210=0.0 */

/* ds->df00120=0.0 */

    ds->df0003 += factor*(
        rob*Exb_03/roa_rob
        );
}

#ifdef FOURTH_ORDER_DERIVATIVES
static void
pbe_xd4(FourthFuncDrv *ds, real factor, const DftDensProp *dp)
{
    TwoVarDrv  dera;
    TwoVarDrv  derb;
    real  Exa;
    real  Exb;
    real  Exa_01;
    real  Exa_02;
    real  Exa_03;
    real  Exa_04;
    real  Exa_10;
    real  Exa_11;
    real  Exa_12;
    real  Exa_13;
    real  Exa_20;
    real  Exa_21;
    real  Exa_22;
    real  Exa_30;
    real  Exa_31;
    real  Exa_40;
    real  Exb_01;
    real  Exb_02;
    real  Exb_03;
    real  Exb_04;
    real  Exb_10;
    real  Exb_11;
    real  Exb_12;
    real  Exb_13;
    real  Exb_20;
    real  Exb_21;
    real  Exb_22;
    real  Exb_30;
    real  Exb_31;
    real  Exb_40;
    real  groa;
    real  grob;
    real  roa;
    real  roa_rob;
    real  roa_rob_p2;
    real  roa_rob_p3;
    real  roa_rob_p4;
    real  roa_rob_p5;
    real  rob;

  
    roa = dp->rhoa;
    rob = dp->rhob;

    groa = dp->grada;
    grob = dp->gradb;
	
    Exa=pbe_ro2energy(&roa, &groa);
    Exb=pbe_ro2energy(&rob, &grob);

    pbe_ro2d(&dera, &roa, &groa);
    pbe_ro2d(&derb, &rob, &grob);

    Exa_10=dera.df10;
    Exa_01=dera.df01;
    Exa_20=dera.df20;
    Exa_11=dera.df11;
    Exa_02=dera.df02;
    Exa_30=dera.df30;
    Exa_21=dera.df21;
    Exa_12=dera.df12;
    Exa_03=dera.df03;
    Exa_40=dera.df40;
    Exa_31=dera.df31;
    Exa_22=dera.df22;
    Exa_13=dera.df13;
    Exa_04=dera.df04;

    Exb_10=derb.df10;
    Exb_01=derb.df01;
    Exb_20=derb.df20;
    Exb_11=derb.df11;
    Exb_02=derb.df02;
    Exb_30=derb.df30;
    Exb_21=derb.df21;
    Exb_12=derb.df12;
    Exb_03=derb.df03;
    Exb_40=derb.df40;
    Exb_31=derb.df31;
    Exb_22=derb.df22;
    Exb_13=derb.df13;
    Exb_04=derb.df04;

/* Powers of total density */
    roa_rob=roa+rob;
    roa_rob_p2=roa_rob*roa_rob;
    roa_rob_p3=roa_rob_p2*roa_rob;
    roa_rob_p4=roa_rob_p3*roa_rob;
    roa_rob_p5=roa_rob_p4*roa_rob;

/* Derivatives */

    ds->df10000 += factor*(
        (rob*Exa - rob*Exb + roa*roa_rob*Exa_10)/roa_rob_p2
        );

    ds->df01000 += factor*(
        (-roa*Exa + roa*Exb + rob*roa_rob*Exb_10)/roa_rob_p2
        );

    ds->df00100 += factor*(
        (roa*Exa_01)/roa_rob
        );

    ds->df00010 += factor*(
        (rob*Exb_01)/roa_rob
        );

/* Second order derivatives */

    ds->df20000 += factor*(
        (1.0/roa_rob_p3)*(-2.0*rob*Exa + 2.0*rob*Exb + 
                          roa_rob*(2.0*rob*Exa_10 + roa*roa_rob*Exa_20))
        );

    ds->df11000 += factor*(
        (1.0/roa_rob_p3)*((roa-rob)*Exa + (rob - roa)*Exb - 
                          roa_rob*(roa*Exa_10 + rob*Exb_10))
        );

    ds->df10100 += factor*(
        (rob*Exa_01 + roa*roa_rob*Exa_11)/roa_rob_p2
        );

    ds->df10010 += factor*(
        -rob*Exb_01/roa_rob_p2
        );

    ds->df02000 += factor*(
        (1.0/roa_rob_p3)*(2.0*roa*Exa - 2.0*roa*Exb + 
                          roa_rob*(2.0*roa*Exb_10 + rob*roa_rob*Exb_20))
        );

    ds->df01100 += factor*(
        -roa*Exa_01/roa_rob_p2
        );

    ds->df01010 += factor*(
        (roa*Exb_01 + rob*roa_rob*Exb_11)/roa_rob_p2
        );

    ds->df00200 += factor*(
        roa*Exa_02/roa_rob
        );

/* ds->df00110=0.0 */

    ds->df00020 += factor*(
        rob*Exb_02/roa_rob
        );

/* Third order derivatives */

    ds->df30000 += factor*(
        (1.0/roa_rob_p4)*(6.0*rob*Exa - 6.0*rob*Exb + 
                          roa_rob*(-6.0*rob*Exa_10 + 
                                   roa_rob*(3.0*rob*Exa_20 + 
                                            roa*roa_rob*Exa_30)))
        );

    ds->df21000 += factor*(
        (1.0/roa_rob_p4)*(-2.0*(roa-2.0*rob)*Exa + 
                          2.0*(roa-2.0*rob)*Exb - 
                          roa_rob*(-2.0*(roa-rob)*Exa_10 - 
                                   2.0*rob*Exb_10 + roa*roa_rob*Exa_20))
        );

    ds->df20100 += factor*(
        (1.0/roa_rob_p3)*(-2.0*rob*Exa_01 + roa_rob*(2.0*rob*Exa_11 + 
                                                     roa*roa_rob*Exa_21))
        );

    ds->df20010 += factor*(
        2.0*rob*Exb_01/roa_rob_p3
        );

    ds->df12000 += factor*(
        (1.0/roa_rob_p4)*(-2.0*(2.0*roa-rob)*Exa + 
                          2.0*(2.0*roa-rob)*Exb - 
                          roa_rob*(-2.0*roa*Exa_10 + 
                                   2.0*(roa-rob)*Exb_10 + rob*roa_rob*Exb_20))
        );


    ds->df11100 += factor*(
        ((roa-rob)*Exa_01 - roa*roa_rob*Exa_11)/roa_rob_p3
        );

    ds->df11010 += factor*(
        ((rob-roa)*Exb_01 - rob*roa_rob*Exb_11)/roa_rob_p3
        );

    ds->df10200 += factor*(
        (rob*Exa_02 + roa*roa_rob*Exa_12)/roa_rob_p2
        );

/* ds->df10110=0.0 */

    ds->df10020 += factor*(
        -rob*Exb_02/roa_rob_p2
        );

    ds->df03000 += factor*(
        (1.0/roa_rob_p4)*(-6.0*roa*Exa + 6.0*roa*Exb + 
                          roa_rob*(-6.0*roa*Exb_10 + 
                                   roa_rob*(3.0*roa*Exb_20 + 
                                            rob*roa_rob*Exb_30)))
        );


    ds->df02100 += factor*(
        2.0*roa*Exa_01/roa_rob_p3
        );

    ds->df02010 += factor*(
        (1.0/roa_rob_p3)*(-2.0*roa*Exb_01 + roa_rob*(2.0*roa*Exb_11 + 
                                                     rob*roa_rob*Exb_21))
        );

    ds->df01200 += factor*(
        -roa*Exa_02/roa_rob_p2
        );

/* ds->df01110=0.0 */

    ds->df01020 += factor*(
        (roa*Exb_02 + rob*roa_rob*Exb_12)/roa_rob_p2
        );

    ds->df00300 += factor*(
        roa*Exa_03/roa_rob
        );

/* ds->df00210=0.0 */

/* ds->df00120=0.0 */

    ds->df00030 += factor*(
        rob*Exb_03/roa_rob
        );

/* Fourth order derivatives */

    ds->df40000 += factor*(
        (1.0/roa_rob_p5)*(-24.0*rob*Exa + 24.0*rob*Exb + 
                          roa_rob*(24.0*rob*Exa_10 +
                                   roa_rob*(-12.0*rob*Exa_20 + 
                                            roa_rob*(4.0*rob*Exa_30 +
                                                     roa*roa_rob*Exa_40))))
        );

    ds->df31000 += factor*(
        (1.0/roa_rob_p5)*(6.0*(roa-3.0*rob)*Exa - 6.0*(roa-3.0*rob)*Exb - 
                          roa_rob*(6.0*(roa - 2.0*rob)*Exa_10 +
                                   6.0*rob*Exb_10 + 
                                   roa_rob*(-3.0*(roa-rob)*Exa_20 +
                                            roa*roa_rob*Exa_30)))
        );

    ds->df30100 += factor*(
        (1.0/roa_rob_p4)*(6.0*rob*Exa_01 + 
                          roa_rob*(-6.0*rob*Exa_11 + 
                                   roa_rob*(3.0*rob*Exa_21 +
                                            roa*roa_rob*Exa_31)))
        );

    ds->df30010 += factor*(
        -6.0*rob*Exb_01/roa_rob_p4
        );

    ds->df22000 += factor*(
        (1.0/roa_rob_p5)*(2.0*(6.0*(roa-rob)*Exa - 6.0*(roa-rob)*Exb + 
                               roa_rob*((-4.0*roa+2.0*rob)*Exa_10 +
                                        2.0*(roa-2.0*rob)*Exb_10 + 
                                        roa_rob*(roa*Exa_20 + rob*Exb_20))))
        );

    ds->df21100 += factor*(
        -(1.0/roa_rob_p4)*(2.0*(roa-2.0*rob)*Exa_01 + 
                           roa_rob*(-2.0*(roa-rob)*Exa_11 +
                                    roa*roa_rob*Exa_21))
        );

    ds->df21010 += factor*(
        (2.0*((roa-2.0*rob)*Exb_01 + rob*roa_rob*Exb_11))/roa_rob_p4
        );

    ds->df20200 += factor*(
        (1.0/roa_rob_p3)*(-2.0*rob*Exa_02 + 
                          roa_rob*(2.0*rob*Exa_12 + roa*roa_rob*Exa_22))
        );

/* ds->df20110=0.0 */

    ds->df20020 += factor*(
        2.0*rob*Exb_02/roa_rob_p3
        );

    ds->df13000 += factor*(
        -(1.0/roa_rob_p5)*(6.0*(-3.0*roa+rob)*Exa + 6.0*(3.0*roa-rob)*Exb + 
                           roa_rob*(6.0*roa*Exa_10 + 
                                    6.0*(-2.0*roa+rob)*Exb_10 + 
                                    roa_rob*(3.0*(roa-rob)*Exb_20 +
                                             rob*roa_rob*Exb_30)))
        );

    ds->df12100 += factor*(
        (2.0*((-2.0*roa+rob)*Exa_01 + roa*roa_rob*Exa_11))/roa_rob_p4
        );

    ds->df12010 += factor*(
        (1.0/roa_rob_p4)*((4.0*roa-2.0*rob)*Exb_01 -
                          roa_rob*(2.0*(roa-rob)*Exb_11 + rob*roa_rob*Exb_21))
        );

    ds->df11200 += factor*(
        ((roa-rob)*Exa_02 - roa*roa_rob*Exa_12)/roa_rob_p3
        );

/* ds->df11110=0.0 */

    ds->df11020 += factor*(
        ((-roa+rob)*Exb_02 - rob*roa_rob*Exb_12)/roa_rob_p3
        );

    ds->df10300 += factor*(
        (rob*Exa_03 + roa*roa_rob*Exa_13)/roa_rob_p2
        );

/* ds->df10210=0.0 */

/* ds->df10120=0.0 */

    ds->df10030 += factor*(
        rob*Exb_03/roa_rob_p2
        );

    ds->df04000 += factor*(
        (1.0/roa_rob_p5)*(24.0*roa*Exa - 24.0*roa*Exb + 
                          roa_rob*(24.0*roa*Exb_10 + 
                                   roa_rob*(-12.0*roa*Exb_20 + 
                                            roa_rob*(4.0*roa*Exb_30 +
                                                     rob*roa_rob*Exb_40))))
        );

    ds->df03100 += factor*(
        -6.0*roa*Exa_01/roa_rob_p4
        );

    ds->df03010 += factor*(
        (1.0/roa_rob_p4)*(6.0*roa*Exb_01 + 
                          roa_rob*(-6.0*roa*Exb_11 + 
                                   roa_rob*(3.0*roa*Exb_21 +
                                            rob*roa_rob*Exb_31)))
        );

    ds->df02200 += factor*(
        2.0*roa*Exa_02/roa_rob_p3
        );

/* ds->df02110=0.0 */

    ds->df02020 += factor*(
        (1.0/roa_rob_p3)*(-2.0*roa*Exb_02 +
                          roa_rob*(2.0*roa*Exb_12 + rob*roa_rob*Exb_22))
        );

    ds->df01300 += factor*(
        -roa*Exa_03/roa_rob_p2
        );

/*
  ds->df01210=0.0
  ds->df01120=0.0
*/

    ds->df01030 += factor*(
        (roa*Exb_03 + rob*roa_rob*Exb_13)/roa_rob_p2
        );

    ds->df00400 += factor*(
        (roa*Exa_04)/roa_rob
        );

/*
  ds->df00310=0.0
  ds->df00220=0.0
  ds->df00130=0.0
*/

    ds->df00040 += factor*(
        (rob*Exb_04)/roa_rob
        );
}
#endif

static void
pbe_cd1(FirstFuncDrv *ds, real factor, const DftDensProp *dp)
{
    const real nKF = COR_nKF, nKS = COR_nKS, delta = COR_delta;
    const real nZ = COR_nZ,   gamma = COR_gamma;
    real A;	real ADENOM;	real ADENOM_p2;	real AEXP;	real A_01;
    real A_10;	real A_p2;	real At2;	real At2_p2;
    real At_denom;	real At_nom;	real At_term;	real EC;
    real EC_01;	real EC_10;	real FI;	real FI_p3;
    real H;		real H0;	real H0_01;	real H0_10;
    real H0_p2;	real H1;	real H_A;	real H_H;
    real H_T;	real KF;	real KS;	real MZETA;
    real MZETA_01;	real MZETA_01_PZETA_01;	real MZETA_10;	real MZETA_10_PZETA_10;
    real MZETA_PZETA;	real MZETA_PZETA_p2;	real PZETA;	real PZETA_01;
    real PZETA_10;	real T;	real T_0001;	real T_0010;
    real T_0100;	real T_1000;	real T_p2;	real T_p3;
    real T_p4;	real T_p5;	real T_p6;	real ZETA;
    real delta_A;	real groa;	real groa_grob;
    real grob;
    real roa;	real roa_rob;	real roa_rob_p13f6;	real roa_rob_p1f3;
    real roa_rob_p1f6;	real roa_rob_p2;	real roa_rob_p7f6;	real roafro;
    real roafro_p1f3;	real roafro_p2f3;	real roafro_p5f3;	real rob;
    real robfro;	real robfro_p1f3;	real robfro_p2f3;	real robfro_p5f3;

    FirstFuncDrv EC_drv;

/* Set up roa, rob, groa, grob */
    roa = dp->rhoa;
    rob = dp->rhob;

    groa = dp->grada;
    grob = dp->gradb;

/* roa, rob and powers */
    roa_rob=roa+rob;
    roa_rob_p1f3=pow(roa_rob,(1.0/3.0));
    roa_rob_p1f6=pow(roa_rob,(1.0/6.0));

    roa_rob_p2=roa_rob*roa_rob;

    roa_rob_p7f6 = roa_rob*roa_rob_p1f6;
    roa_rob_p13f6 = roa_rob_p2*roa_rob_p1f6;

/* groa, grob gradients and powers */
    groa_grob = groa+ grob;

/* Auxaliary functions */
    ZETA=(roa-rob)/(roa+rob);
    KF=nKF*roa_rob_p1f3;
    KS=nKS*roa_rob_p1f6;
    MZETA=pow((1.0-ZETA),(2.0/3.0));
    PZETA=pow((1.0+ZETA),(2.0/3.0));
    FI=(PZETA + MZETA)/2.0;
    T=(groa+grob)/((roa+rob)*(2.0*KS*FI));


/* Powers of Fi */
    FI_p3=FI*FI*FI;

/* Powers of T */
    T_p2=T*T;
    T_p3=T_p2*T;
    T_p4=T_p2*T_p2;
    T_p5=T_p4*T;
    T_p6=T_p4*T_p2;

/* The pw92 EC(electron gas corelation) energy */
    EC = PW92Functional.func(dp);

/* H0 and its powers */
    H0=gamma*FI_p3;
    H0_p2 = H0*H0;

/* Calc of A */
    AEXP=exp(EC/H0);
    ADENOM=AEXP -1.0;
    ADENOM_p2=ADENOM*ADENOM;
    A=delta*( 1.0/((1.0/AEXP) -1.0)  ) ;

/* Powers of A */
    A_p2=A*A;

/* The pbe Ec (pbe correlation contribution) */
/* H is gga contribution to correlation energy */ 
    H1=1.0 + delta*T_p2*((1.0 + A*T_p2)/(1.0 + A*T_p2 + A_p2*T_p4));

    H=H0*log( H1 );
/* E=EC+H; */


/* Now start with derivatives *
 * first, common terms used in equations */

/* roafro, robfro and powers */
    roafro=roa/roa_rob;
    robfro=rob/roa_rob;

    roafro_p1f3=pow(roafro,(1.0/3.0));
    robfro_p1f3=pow(robfro,(1.0/3.0));

    roafro_p2f3=roafro_p1f3*roafro_p1f3;
    robfro_p2f3=robfro_p1f3*robfro_p1f3;

    roafro_p5f3=roafro_p2f3*roafro;
    robfro_p5f3=robfro_p2f3*robfro;


/* Derivatives of PZETA, MZETA */

    MZETA_10=(-2.0/3.0)*nZ*robfro_p5f3/rob;
    PZETA_10=(2.0/3.0)*nZ*rob/(roafro_p1f3*roa_rob_p2);

    MZETA_01=(2.0/3.0)*nZ*roa/(robfro_p1f3*roa_rob_p2);
    PZETA_01=(-2.0/3.0)*nZ*roafro_p5f3/roa;

/* MZETA, PZETA terms and derivatives */
    MZETA_PZETA = MZETA + PZETA;
    MZETA_10_PZETA_10 = MZETA_10 + PZETA_10;
    MZETA_01_PZETA_01 = MZETA_01 + PZETA_01;

    MZETA_PZETA_p2 = MZETA_PZETA*MZETA_PZETA;


/* Derivatives of EC */

/* EC_drv = pw92_d(roa,rob); */
    drv1_clear(&EC_drv);
    PW92Functional.first(&EC_drv, 1.0, dp);

    EC_10=EC_drv.df1000;
    EC_01=EC_drv.df0100;

/* Derivatives of H0 */

    H0_10=(3.0/8.0)*gamma*MZETA_PZETA_p2*(MZETA_10_PZETA_10);

    H0_01=(3.0/8.0)*gamma*MZETA_PZETA_p2*(MZETA_01_PZETA_01);

/* All terms OK up to here ! */

/* Derivatives of T */

    T_1000=-(groa_grob*(7.0*MZETA+7.0*PZETA+6.0*roa_rob*(MZETA_10_PZETA_10)))/(6.0*nKS*roa_rob_p13f6*MZETA_PZETA_p2);

    T_0100=-(groa_grob*(7.0*MZETA+7.0*PZETA+6.0*roa_rob*(MZETA_01_PZETA_01)))/(6.0*nKS*roa_rob_p13f6*MZETA_PZETA_p2);

    T_0010=1.0/(nKS*roa_rob_p7f6*MZETA_PZETA);

    T_0001=T_0010;

/* All terms ok up to here! */

/* Derivatives of A */

    A_10=(delta*AEXP*(H0*EC_10-EC*H0_10))/(ADENOM_p2*H0_p2);

    A_01=(delta*AEXP*(H0*EC_01-EC*H0_01))/(ADENOM_p2*H0_p2);


/* Derivatives of H with respect to H0, T, A */
/* Derivatives with respect to H0 higher then first order vanish */
/* Terms used in the derivatives */
    At2 = A*T_p2;
    At2_p2 = At2*At2;

    At_nom = 1.0 + At2;
    At_denom = At_nom + At2_p2;
    delta_A = delta + A;
    At_term = (1.0+delta_A*T_p2*At_nom);

/* First order derivatives */
    H_A =-((delta*A*H0*T_p6*(2.0+At2))/(At_denom*At_term));

    H_T =(2.0*delta*H0*T*(1.0+2.0*At2))/(At_denom*At_term);

    H_H = log(1.0+(delta-delta/At_denom)/A);

/* Final expressions (general derivatives of F[H[roa,rob], T[roa,rob,groa,grob], A[roa,rob]]) */

/* First order derivatives */
    ds->df1000 += factor*(
        A_10*H_A+H0_10*H_H+H_T*T_1000 + EC_10
        );

    ds->df0100 += factor*(
        A_01*H_A+H0_01*H_H+H_T*T_0100 + EC_01
        );

    ds->df0010 += factor*(
        H_T*T_0010
        );

    ds->df0001 += factor*(
        H_T*T_0001
        );
/* End (Finally!) */
}

static void
pbe_cd2(SecondFuncDrv *ds, real factor, const DftDensProp *dp)
{
    const real nKF = COR_nKF, nKS = COR_nKS, delta = COR_delta;
    const real nZ = COR_nZ,   gamma = COR_gamma;

    real A;	real ADENOM;	real ADENOM_p2;	real ADENOM_p3;	real AEXP;
    real AEXP_p2;	real AEXP_p3;	real A_01;	real A_01_p2;
    real A_02;	real A_10;	real A_10_p2;	real A_11;
    real A_20;	real A_p2;	real A_p3;	real A_p4;
    real At2;	real At2_p2;	real At_denom;	real At_denom_p2;
    real At_nom;	real At_term;	real At_term_p2;	real EC;
    real EC_01;	real EC_01_p2;	real EC_02;	real EC_10;
    real EC_10_p2;	real EC_11;	real EC_20;	real EC_p2;
    real FI;	real FI_p3;	real H;		real H0;
    real H0_01;	real H0_01_p2;	real H0_02;	real H0_10;
    real H0_10_p2;	real H0_11;	real H0_20;	real H0_p2;
    real H0_p3;	real H0_p4;	real H1;	real H_A;
    real H_AA;	real H_H;	real H_HA;	real H_HT;
    real H_T;	real H_TA;	real H_TT;	real KF;
    real KS;	real MZETA;	real MZETA_01;	real MZETA_01_PZETA_01;
    real MZETA_01_PZETA_01_p2;	real MZETA_02;	real MZETA_02_PZETA_02;	real MZETA_10;
    real MZETA_10_PZETA_10;	real MZETA_10_PZETA_10_p2;	real MZETA_11;	real MZETA_11_PZETA_11;
    real MZETA_20;	real MZETA_20_PZETA_20;	real MZETA_PZETA;	real MZETA_PZETA_p2;
    real MZETA_PZETA_p3;	real PZETA;	real PZETA_01;	real PZETA_02;
    real PZETA_10;	real PZETA_11;	real PZETA_20;	real T;
    real T_0001;	real T_0001_p2;	real T_0010;	real T_0010_p2;
    real T_0100;	real T_0100_p2;	real T_0101;	real T_0110;
    real T_0200;	real T_1000;	real T_1000_p2;	real T_1001;
    real T_1010;	real T_1100;	real T_2000;	real T_p10;
    real T_p2;	real T_p3;	real T_p4;	real T_p5;
    real T_p6;	real T_p7;	real T_p8;	real ZETA;
    real delta_A;	real groa;
    real groa_grob;	real grob;
    real roa;	real roa_rob;	real roa_rob_p13f6;
    real roa_rob_p19f6;	real roa_rob_p1f3;	real roa_rob_p1f6;	real roa_rob_p2;
    real roa_rob_p3;	real roa_rob_p4;	real roa_rob_p7f6;	real roafro;
    real roafro_p1f3;	real roafro_p2f3;	real roafro_p4f3;	real roafro_p5f3;
    real rob;	real robfro;	real robfro_p1f3;	real robfro_p2f3;
    real robfro_p4f3;	real robfro_p5f3;

    SecondFuncDrv EC_drv;
/* End of Declarations */

	
/* Set up roa, rob, groa, grob */
    roa = dp->rhoa;
    rob = dp->rhob;

    groa = dp->grada;
    grob = dp->gradb;

/* roa, rob and powers */
    roa_rob=roa+rob;
    roa_rob_p1f3=pow(roa_rob,(1.0/3.0));
    roa_rob_p1f6=pow(roa_rob,(1.0/6.0));

    roa_rob_p2=roa_rob*roa_rob;
    roa_rob_p3=roa_rob_p2*roa_rob;
    roa_rob_p4=roa_rob_p2*roa_rob_p2;

    roa_rob_p7f6 = roa_rob*roa_rob_p1f6;
    roa_rob_p13f6 = roa_rob_p2*roa_rob_p1f6;
    roa_rob_p19f6 = roa_rob_p3*roa_rob_p1f6;

/* groa, grob gradients and powers */
    groa_grob = groa+ grob;

/* Auxaliary functions */
    ZETA=(roa-rob)/(roa+rob);
    KF=nKF*roa_rob_p1f3;
    KS=nKS*roa_rob_p1f6;
    MZETA=pow((1.0-ZETA),(2.0/3.0));
    PZETA=pow((1.0+ZETA),(2.0/3.0));
    FI=(PZETA + MZETA)/2.0;
    T=(groa+grob)/((roa+rob)*(2.0*KS*FI));


/* Powers of Fi */
    FI_p3=FI*FI*FI;

/* Powers of T */
    T_p2=T*T;
    T_p3=T_p2*T;
    T_p4=T_p2*T_p2;
    T_p5=T_p4*T;
    T_p6=T_p4*T_p2;
    T_p7=T_p6*T;
    T_p8=T_p6*T_p2;
    T_p10=T_p8*T_p2;

/* The pw92 EC(electron gas corelation) energy */
    EC = PW92Functional.func(dp);
    EC_p2 = EC*EC;

/* H0 and its powers */
    H0=gamma*FI_p3;
    H0_p2 = H0*H0;
    H0_p3 = H0_p2*H0;
    H0_p4 = H0_p3*H0;

/* Calc of A */
    AEXP=exp(EC/H0);
    AEXP_p2=AEXP*AEXP;
    AEXP_p3=AEXP_p2*AEXP;
    ADENOM=AEXP -1.0;
    ADENOM_p2=ADENOM*ADENOM;
    ADENOM_p3=ADENOM_p2*ADENOM;
    A=delta*( 1.0/((1.0/AEXP) -1.0)  ) ;

/* Powers of A */
    A_p2=A*A;
    A_p3=A_p2*A;
    A_p4=A_p3*A;

/* The pbe Ec (pbe correlation contribution) */
/* H is gga contribution to correlation energy */ 
    H1=1.0 + delta*T_p2*((1.0 + A*T_p2)/(1.0 + A*T_p2 + A_p2*T_p4));

    H=H0*log( H1 );
/* E=EC+H; */


/* Now start with derivatives *
 * first, common terms used in equations */

/* roafro, robfro and powers */
    roafro=roa/roa_rob;
    robfro=rob/roa_rob;

    roafro_p1f3=pow(roafro,(1.0/3.0));
    robfro_p1f3=pow(robfro,(1.0/3.0));

    roafro_p2f3=roafro_p1f3*roafro_p1f3;
    robfro_p2f3=robfro_p1f3*robfro_p1f3;

    roafro_p4f3=roafro_p1f3*roafro;
    robfro_p4f3=robfro_p1f3*robfro;

    roafro_p5f3=roafro_p2f3*roafro;
    robfro_p5f3=robfro_p2f3*robfro;


/* Derivatives of PZETA, MZETA */

    MZETA_10=(-2.0/3.0)*nZ*robfro_p5f3/rob;
    PZETA_10=(2.0/3.0)*nZ*rob/(roafro_p1f3*roa_rob_p2);

    MZETA_01=(2.0/3.0)*nZ*roa/(robfro_p1f3*roa_rob_p2);
    PZETA_01=(-2.0/3.0)*nZ*roafro_p5f3/roa;

    MZETA_20=(10.0/9.0)*nZ*robfro_p2f3/roa_rob_p2;
    PZETA_20=(-2.0/9.0)*nZ*rob*(6.0*roa + rob)/(roafro_p4f3*roa_rob_p4);

    MZETA_11=(2.0/9.0)*nZ*(-2.0*roa+3.0*rob)/(robfro_p1f3*roa_rob_p3);
    PZETA_11=(2.0/9.0)*nZ*(3.0*roa-2.0*rob)/(roafro_p1f3*roa_rob_p3);

    MZETA_02=(-2.0/9.0)*nZ*roa*(6.0*rob + roa)/(robfro_p4f3*roa_rob_p4);
    PZETA_02=(10.0/9.0)*nZ*roafro_p2f3/roa_rob_p2;

/* MZETA, PZETA terms and derivatives */
    MZETA_PZETA = MZETA + PZETA;
    MZETA_10_PZETA_10 = MZETA_10 + PZETA_10;
    MZETA_01_PZETA_01 = MZETA_01 + PZETA_01;

    MZETA_20_PZETA_20 = MZETA_20 + PZETA_20;
    MZETA_02_PZETA_02 = MZETA_02 + PZETA_02;
    MZETA_11_PZETA_11 = MZETA_11 + PZETA_11;


    MZETA_PZETA_p2 = MZETA_PZETA*MZETA_PZETA;
    MZETA_10_PZETA_10_p2 = MZETA_10_PZETA_10*MZETA_10_PZETA_10;
    MZETA_01_PZETA_01_p2 = MZETA_01_PZETA_01*MZETA_01_PZETA_01;

    MZETA_PZETA_p3 = MZETA_PZETA_p2*MZETA_PZETA;


/* Derivatives of EC */

/* EC_drv = pw92_d(roa,rob); */
    drv2_clear(&EC_drv);
    PW92Functional.second(&EC_drv, 1.0, dp);

    EC_10=EC_drv.df1000;
    EC_01=EC_drv.df0100;
    EC_20=EC_drv.df2000;
    EC_02=EC_drv.df0200;
    EC_11=EC_drv.df1100;

/* Powers of EC derivatives */

    EC_10_p2=EC_10*EC_10;
    EC_01_p2=EC_01*EC_01;

/* Derivatives of H0 */

    H0_10=(3.0/8.0)*gamma*MZETA_PZETA_p2*(MZETA_10_PZETA_10);

    H0_01=(3.0/8.0)*gamma*MZETA_PZETA_p2*(MZETA_01_PZETA_01);

    H0_20=(3.0/8.0)*gamma*MZETA_PZETA*(2.0*MZETA_10_PZETA_10_p2 + 
                                       MZETA_PZETA*MZETA_20_PZETA_20);

    H0_11=(3.0/8.0)*gamma*MZETA_PZETA*(2.0*MZETA_01_PZETA_01*MZETA_10_PZETA_10 + 
                                       MZETA_PZETA*MZETA_11_PZETA_11);

    H0_02=(3.0/8.0)*gamma*MZETA_PZETA*(2.0*MZETA_01_PZETA_01_p2 + 
                                       MZETA_PZETA*MZETA_02_PZETA_02);

/* All terms OK up to here ! */

/* Powers of H0 derivatives */
    H0_10_p2=H0_10*H0_10;
    H0_01_p2=H0_01*H0_01;

/* Derivatives of T */

    T_1000=-(groa_grob*(7.0*MZETA+7.0*PZETA+6.0*roa_rob*(MZETA_10_PZETA_10)))/(6.0*nKS*roa_rob_p13f6*MZETA_PZETA_p2);

    T_0100=-(groa_grob*(7.0*MZETA+7.0*PZETA+6.0*roa_rob*(MZETA_01_PZETA_01)))/(6.0*nKS*roa_rob_p13f6*MZETA_PZETA_p2);

    T_0010=1.0/(nKS*roa_rob_p7f6*MZETA_PZETA);

    T_0001=T_0010;

    T_2000=(groa_grob*(91.0*MZETA_PZETA_p2+
                       84.0*roa_rob*(MZETA_PZETA)*(MZETA_10_PZETA_10)+
                       72.0*roa_rob_p2*MZETA_10_PZETA_10_p2-
                       36.0*roa_rob_p2*(MZETA_PZETA)*(MZETA_20_PZETA_20)))
        /(36.0*nKS*roa_rob_p19f6*MZETA_PZETA_p3);

    T_1100=(groa_grob*(91.0*MZETA_PZETA_p2+
                       42.0*roa_rob*(MZETA_PZETA)*(MZETA_01_PZETA_01)+
                       42.0*roa_rob*(MZETA_PZETA)*(MZETA_10_PZETA_10)+
                       72.0*roa_rob_p2*(MZETA_01_PZETA_01)*(MZETA_10_PZETA_10)-
                       36.0*roa_rob_p2*(MZETA_PZETA)*(MZETA_11_PZETA_11)))
        /(36.0*nKS*roa_rob_p19f6*MZETA_PZETA_p3);

    T_1010=(-7.0*MZETA-7.0*PZETA-6.0*roa_rob*(MZETA_10_PZETA_10))
        /(6.0*nKS*roa_rob_p13f6*MZETA_PZETA_p2);

    T_1001=T_1010;

    T_0200=(groa_grob*(91.0*MZETA_PZETA_p2+
                       84.0*roa_rob*(MZETA_PZETA)*(MZETA_01_PZETA_01)+
                       72.0*roa_rob_p2*MZETA_01_PZETA_01_p2-
                       36.0*roa_rob_p2*(MZETA_PZETA)*(MZETA_02_PZETA_02)))
        /(36.0*nKS*roa_rob_p19f6*MZETA_PZETA_p3);

    T_0110=(-7.0*MZETA-7.0*PZETA-6.0*roa_rob*(MZETA_01_PZETA_01))
        /(6.0*nKS*roa_rob_p13f6*MZETA_PZETA_p2);

    T_0101=T_0110;

/* All terms ok up to here! */
/* Powers of T derivatives */
    T_1000_p2 = T_1000*T_1000;

    T_0100_p2 = T_0100*T_0100;

    T_0010_p2 = T_0010*T_0010;

    T_0001_p2 = T_0001*T_0001;

/* Derivatives of A */

    A_10=(delta*AEXP*(H0*EC_10-EC*H0_10))/(ADENOM_p2*H0_p2);

    A_01=(delta*AEXP*(H0*EC_01-EC*H0_01))/(ADENOM_p2*H0_p2);

    A_20=(delta*AEXP*(-((1.0+AEXP)*EC_p2*H0_10_p2)+
                      2.0*EC*H0*H0_10*((1.0+AEXP)*EC_10+ADENOM*H0_10)+
                      ADENOM*H0_p3*EC_20-
                      H0_p2*((1.0+AEXP)*EC_10_p2+2.0*ADENOM*EC_10*H0_10+
                             ADENOM*EC*H0_20)))/(ADENOM_p3*H0_p4);

    A_02=(delta*AEXP*(-((1.0+AEXP)*EC_p2*H0_01_p2)+
                      2.0*EC*H0*H0_01*((1.0+AEXP)*EC_01+ADENOM*H0_01)+
                      ADENOM*H0_p3*EC_02-
                      H0_p2*((1.0+AEXP)*EC_01_p2+2.0*ADENOM*EC_01*H0_01+
                             ADENOM*EC*H0_02)))/(ADENOM_p3*H0_p4);

    A_11=(delta*AEXP*(-((1.0+AEXP)*EC_p2*H0_01*H0_10)+
                      EC*H0*((1.0+AEXP)*EC_01*H0_10+
                             H0_01*((1.0+AEXP)*EC_10+2.0*ADENOM*H0_10))+
                      ADENOM*H0_p3*EC_11-
                      H0_p2*(EC_01*((1.0+AEXP)*EC_10+ADENOM*H0_10)+
                             ADENOM*(H0_01*EC_10+EC*H0_11))))
        /(ADENOM_p3*H0_p4);

/* Powers of A derivatives */
    A_10_p2 = A_10*A_10;

    A_01_p2 = A_01*A_01;


/* Derivatives of H with respect to H0, T, A */
/* Derivatives with respect to H0 higher then first order vanish */
/* Terms used in the derivatives */
    At2 = A*T_p2;
    At2_p2=At2*At2;

    At_nom = 1.0 + At2;
    At_denom = At_nom + At2_p2;
    delta_A = delta + A;
    At_term = (1.0+delta_A*T_p2*At_nom);

/* Powers of terms */
    At_denom_p2 = At_denom*At_denom;

    At_term_p2 = At_term*At_term;

/* First order derivatives */
    H_A =-((delta*A*H0*T_p6*(2.0+At2))/(At_denom*At_term));

    H_T =(2.0*delta*H0*T*(1.0+2.0*At2))/(At_denom*At_term);

    H_H = log(1.0+(delta-delta/At_denom)/A);

/* Second order derivatives */
    H_AA =(delta*H0*T_p6*(-2.0-2.0*delta_A*T_p2-
                          2.0*(delta-2.0*A)*A*T_p4+
                          2.0*A_p2*(delta+4.0*A)*T_p6+
                          4.0*A_p3*(delta+2.0*A)*T_p8+
                          A_p4*(delta+2.0*A)*T_p10))/(At_denom_p2*At_term_p2);

    H_TT =2.0*delta*H0*
        (1.0+T_p2*(-delta+A*(4.0+T_p2*(-4.0*delta+
                                       A*(-5.0-2.0*(7.0*delta+8.0*A)*T_p2-
                                          19.0*A*delta_A*T_p4-
                                          10.0*A_p2*delta_A*T_p6)))))
        /(At_denom_p2*At_term_p2);

    H_HT =(2.0*delta*T*(1.0+2.0*At2))/(At_denom*At_term);

    H_TA =(-2.0*delta*A*H0*T_p5*(6.0+4.0*(delta+3.0*A)*T_p2+
                                 A*(7.0*delta+12.0*A)*T_p4+
                                 2.0*A_p2*(2.0*delta+3.0*A)*T_p6))
        /(At_denom_p2*At_term_p2);

    H_HA =-((delta*A*T_p6*(2.0+At2))/(At_denom*At_term));

/* Final expressions (general derivatives of F[H[roa,rob],
 * T[roa,rob,groa,grob], A[roa,rob]]) */

/* First order derivatives */
    ds->df1000 += factor*(
        A_10*H_A+H0_10*H_H+H_T*T_1000 + EC_10
        );

    ds->df0100 += factor*(
        A_01*H_A+H0_01*H_H+H_T*T_0100 + EC_01
        );

    ds->df0010 += factor*(
        H_T*T_0010
        );

    ds->df0001 += factor*(
        H_T*T_0001
        );

/* Second order derivatives */
    ds->df2000 += factor*(
        A_20*H_A+H0_20*H_H+A_10*(A_10*H_AA+2.0*H0_10*H_HA)+
        2.0*(A_10*H_TA+H0_10*H_HT)*T_1000+H_TT*T_1000_p2+H_T*T_2000 + EC_20
        );

    ds->df1100 += factor*(
        A_11*H_A+H0_11*H_H+A_01*(A_10*H_AA+H0_10*H_HA+H_TA*T_1000)+
        T_0100*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+
        H0_01*(A_10*H_HA+H_HT*T_1000)+H_T*T_1100 + EC_11
        );

    ds->df1010 += factor*(
        T_0010*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+H_T*T_1010
        );

    ds->df1001 += factor*(
        T_0001*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+H_T*T_1001
        );

    ds->df0200 += factor*(
        A_02*H_A+H0_02*H_H+A_01*(A_01*H_AA+2.0*H0_01*H_HA)+
        2.0*(A_01*H_TA+H0_01*H_HT)*T_0100+H_TT*T_0100_p2+H_T*T_0200 + EC_02
        );

    ds->df0110 += factor*(
        T_0010*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)+H_T*T_0110
        );

    ds->df0101 += factor*(
        T_0001*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)+H_T*T_0101
        );

    ds->df0020 += factor*(
        H_TT*T_0010_p2
        );

    ds->df0011 += factor*(
        H_TT*T_0001*T_0010
        );

    ds->df0002 += factor*(
        H_TT*T_0001_p2
        );
/* End (Finally!) */
}

static void
pbe_cd3(ThirdFuncDrv *ds, real factor, const DftDensProp *dp)
{
    const real nKF = COR_nKF, nKS = COR_nKS, delta = COR_delta;
    const real nZ = COR_nZ,   gamma = COR_gamma;

    real A;	real ADENOM;	real ADENOM_p2;	real ADENOM_p3;	real ADENOM_p4;
    real AEXP;	real AEXP_p2;	real AEXP_p3;	real A_01;
    real A_01_p2;	real A_01_p3;	real A_02;	real A_03;
    real A_10;	real A_10_p2;	real A_10_p3;	real A_11;
    real A_12;	real A_20;	real A_21;	real A_30;
    real A_p2;	real A_p3;	real A_p4;	real A_p5;
    real A_p6;	real At2;	real At2_p2;	real At_denom;
    real At_denom_p2;	real At_denom_p3;	real At_nom;	real At_term;
    real At_term_p2;	real At_term_p3;	real EC;	real EC_01;
    real EC_01_p2;	real EC_01_p3;	real EC_02;	real EC_03;
    real EC_10;	real EC_10_p2;	real EC_10_p3;	real EC_11;
    real EC_12;	real EC_20;	real EC_21;	real EC_30;
    real EC_p2;	real EC_p3;	real FI;	real FI_p3;
    real H;	real H0;	real H0_01;	real H0_01_p2;
    real H0_01_p3;	real H0_02;	real H0_03;	real H0_10;
    real H0_10_p2;	real H0_10_p3;	real H0_11;	real H0_12;
    real H0_20;	real H0_21;	real H0_30;	real H0_p2;
    real H0_p3;	real H0_p4;	real H0_p5;	real H0_p6;
    real H1;	real H_A;	real H_AA;	real H_AAA;
    real H_H;	real H_HA;	real H_HAA;	real H_HT;
    real H_HTA;	real H_HTT;	real H_T;	real H_TA;
    real H_TAA;	real H_TT;	real H_TTA;	real H_TTT;
    real KF;	real KS;	real MZETA;	real MZETA_01;
    real MZETA_01_PZETA_01;	real MZETA_01_PZETA_01_p2;
    real MZETA_01_PZETA_01_p3;	real MZETA_02;
    real MZETA_02_PZETA_02;	real MZETA_03;
    real MZETA_03_PZETA_03;	real MZETA_10;
    real MZETA_10_PZETA_10;	real MZETA_10_PZETA_10_p2;
    real MZETA_10_PZETA_10_p3;	real MZETA_11;
    real MZETA_11_PZETA_11;	real MZETA_12;
    real MZETA_12_PZETA_12;	real MZETA_20;
    real MZETA_20_PZETA_20;	real MZETA_21;
    real MZETA_21_PZETA_21;	real MZETA_30;
    real MZETA_30_PZETA_30;	real MZETA_PZETA;
    real MZETA_PZETA_p2;	real MZETA_PZETA_p3;
    real MZETA_PZETA_p4;	real PZETA;	real PZETA_01;	real PZETA_02;
    real PZETA_03;	real PZETA_10;	real PZETA_11;	real PZETA_12;
    real PZETA_20;	real PZETA_21;	real PZETA_30;	real T;
    real T_0001;	real T_0001_p2;	real T_0001_p3;	real T_0010;
    real T_0010_p2;	real T_0010_p3;	real T_0100;	real T_0100_p2;
    real T_0100_p3;	real T_0101;	real T_0110;	real T_0200;
    real T_0201;	real T_0210;	real T_0300;	real T_1000;
    real T_1000_p2;	real T_1000_p3;	real T_1001;	real T_1010;
    real T_1100;	real T_1101;	real T_1110;	real T_1200;
    real T_2000;	real T_2001;	real T_2010;	real T_2100;
    real T_3000;	real T_p10;	real T_p12;	real T_p14;
    real T_p16;	real T_p2;	real T_p3;	real T_p4;
    real T_p5;	real T_p6;	real T_p7;	real T_p8;
    real ZETA;	real delta_A;	real delta_A_p2;
    real groa;	real groa_grob;	real grob;
    real roa;
    real roa_p2,        roa_p3,       roa_rob,      roa_rob_p13f6;
    real roa_rob_p19f6, roa_rob_p1f3, roa_rob_p1f6, roa_rob_p2;
    real roa_rob_p25f6, roa_rob_p3,   roa_rob_p4,   roa_rob_p5;
    real roa_rob_p6,    roa_rob_p7f6, roafro,       roafro_p1f3;
    real roafro_p2f3,   roafro_p4f3,  roafro_p5f3,  roafro_p7f3;
    real rob,  rob_p2,  rob_p3,       robfro;
    real robfro_p1f3,   robfro_p2f3, robfro_p4f3,   robfro_p5f3;
    real robfro_p7f3;

    ThirdFuncDrv EC_drv;

    const real delta_p2 = delta*delta;

/* Set up roa, rob, groa, grob */
    roa = dp->rhoa;
    rob = dp->rhob;

    groa = dp->grada;
    grob = dp->gradb;

/* roa, rob and powers */
    roa_rob=roa+rob;
    roa_rob_p1f3=pow(roa_rob,(1.0/3.0));
    roa_rob_p1f6=pow(roa_rob,(1.0/6.0));

    roa_rob_p2=roa_rob*roa_rob;
    roa_rob_p3=roa_rob_p2*roa_rob;
    roa_rob_p4=roa_rob_p2*roa_rob_p2;
    roa_rob_p5=roa_rob_p4*roa_rob;
    roa_rob_p6=roa_rob_p4*roa_rob_p2;

    roa_rob_p7f6 = roa_rob*roa_rob_p1f6;
    roa_rob_p13f6 = roa_rob_p2*roa_rob_p1f6;
    roa_rob_p19f6 = roa_rob_p3*roa_rob_p1f6;
    roa_rob_p25f6 = roa_rob_p4*roa_rob_p1f6;

    roa_p2=roa*roa;
    rob_p2=rob*rob;

    roa_p3=roa_p2*roa;
    rob_p3=rob_p2*rob;

/* groa, grob gradients and powers */
    groa_grob = groa+ grob;

/* Auxaliary functions */
    ZETA=(roa-rob)/(roa+rob);
    KF=nKF*roa_rob_p1f3;
    KS=nKS*roa_rob_p1f6;
    MZETA=pow((1.0-ZETA),(2.0/3.0));
    PZETA=pow((1.0+ZETA),(2.0/3.0));
    FI=(PZETA + MZETA)/2.0;
    T=(groa+grob)/((roa+rob)*(2.0*KS*FI));


/* Powers of Fi */
    FI_p3=FI*FI*FI;

/* Powers of T */
    T_p2=T*T;
    T_p3=T_p2*T;
    T_p4=T_p2*T_p2;
    T_p5=T_p4*T;
    T_p6=T_p4*T_p2;
    T_p7=T_p6*T;
    T_p8=T_p6*T_p2;
    T_p10=T_p8*T_p2;
    T_p12=T_p10*T_p2;
    T_p14=T_p12*T_p2;
    T_p16=T_p14*T_p2;

/* The pw92 EC(electron gas corelation) energy */
    EC = PW92Functional.func(dp);
    EC_p2 = EC*EC;
    EC_p3 = EC_p2*EC;

/* H0 and its powers */
    H0=gamma*FI_p3;
    H0_p2 = H0*H0;
    H0_p3 = H0_p2*H0;
    H0_p4 = H0_p3*H0;
    H0_p5 = H0_p4*H0;
    H0_p6 = H0_p5*H0;

/* Calc of A */
    AEXP=exp(EC/H0);
    AEXP_p2=AEXP*AEXP;
    AEXP_p3=AEXP_p2*AEXP;
    ADENOM=AEXP -1.0;
    ADENOM_p2=ADENOM*ADENOM;
    ADENOM_p3=ADENOM_p2*ADENOM;
    ADENOM_p4=ADENOM_p3*ADENOM;
    A=delta*( 1.0/((1.0/AEXP) -1.0)  ) ;

/* Powers of A */
    A_p2=A*A;
    A_p3=A_p2*A;
    A_p4=A_p3*A;
    A_p5=A_p4*A;
    A_p6=A_p5*A;

/* The pbe Ec (pbe correlation contribution) */
/* H is gga contribution to correlation energy */ 
    H1=1.0 + delta*T_p2*((1.0 + A*T_p2)/(1.0 + A*T_p2 + A_p2*T_p4));

    H=H0*log( H1 );
/* E=EC+H; */


/* Now start with derivatives *
 * first, common terms used in equations */

/* roafro, robfro and powers */
    roafro=roa/roa_rob;
    robfro=rob/roa_rob;

    roafro_p1f3=pow(roafro,(1.0/3.0));
    robfro_p1f3=pow(robfro,(1.0/3.0));

    roafro_p2f3=roafro_p1f3*roafro_p1f3;
    robfro_p2f3=robfro_p1f3*robfro_p1f3;

    roafro_p4f3=roafro_p1f3*roafro;
    robfro_p4f3=robfro_p1f3*robfro;

    roafro_p5f3=roafro_p2f3*roafro;
    robfro_p5f3=robfro_p2f3*robfro;

    roafro_p7f3=roafro_p4f3*roafro;
    robfro_p7f3=robfro_p4f3*robfro;

/* Derivatives of PZETA, MZETA */

    MZETA_10=(-2.0/3.0)*nZ*robfro_p5f3/rob;
    PZETA_10=(2.0/3.0)*nZ*rob/(roafro_p1f3*roa_rob_p2);

    MZETA_01=(2.0/3.0)*nZ*roa/(robfro_p1f3*roa_rob_p2);
    PZETA_01=(-2.0/3.0)*nZ*roafro_p5f3/roa;

    MZETA_20=(10.0/9.0)*nZ*robfro_p2f3/roa_rob_p2;
    PZETA_20=(-2.0/9.0)*nZ*rob*(6.0*roa + rob)/(roafro_p4f3*roa_rob_p4);

    MZETA_11=(2.0/9.0)*nZ*(-2.0*roa+3.0*rob)/(robfro_p1f3*roa_rob_p3);
    PZETA_11=(2.0/9.0)*nZ*(3.0*roa-2.0*rob)/(roafro_p1f3*roa_rob_p3);

    MZETA_02=(-2.0/9.0)*nZ*roa*(6.0*rob + roa)/(robfro_p4f3*roa_rob_p4);
    PZETA_02=(10.0/9.0)*nZ*roafro_p2f3/roa_rob_p2;

    MZETA_30=(-80.0/27.0)*nZ*robfro_p2f3/roa_rob_p3;
    PZETA_30=(4.0/27.0)*nZ*rob*(27.0*roa_p2 + 9.0*roa*rob+2.0*rob_p2)
        /(roafro_p7f3*roa_rob_p6);

    MZETA_21=(20.0/27.0)*nZ*(roa-3.0*rob)/(robfro_p1f3*roa_rob_p4);
    PZETA_21=(4.0/27.0)*nZ*(-9.0*roa_p2 + 12.0*roa*rob+ rob_p2)
        /(roafro_p4f3*roa_rob_p5);

    MZETA_12=(4.0/27.0)*nZ*(roa_p2 + 12.0*roa*rob -9.0*rob_p2)
        /(robfro_p4f3*roa_rob_p5);
    PZETA_12=(20.0/27.0)*nZ*(-3.0*roa+rob)/(roafro_p1f3*roa_rob_p4);

    MZETA_03=(4.0/27.0)*nZ*roa*(2.0*roa_p2 + 9.0*roa*rob +27.0*rob_p2)
        /(robfro_p7f3*roa_rob_p6);
    PZETA_03=(-80.0/27.0)*nZ*roafro_p2f3/roa_rob_p3;

/* MZETA, PZETA terms and derivatives */
    MZETA_PZETA = MZETA + PZETA;
    MZETA_10_PZETA_10 = MZETA_10 + PZETA_10;
    MZETA_01_PZETA_01 = MZETA_01 + PZETA_01;

    MZETA_20_PZETA_20 = MZETA_20 + PZETA_20;
    MZETA_02_PZETA_02 = MZETA_02 + PZETA_02;
    MZETA_11_PZETA_11 = MZETA_11 + PZETA_11;

    MZETA_30_PZETA_30 = MZETA_30 + PZETA_30;
    MZETA_03_PZETA_03 = MZETA_03 + PZETA_03;
    MZETA_21_PZETA_21 = MZETA_21 + PZETA_21;
    MZETA_12_PZETA_12 = MZETA_12 + PZETA_12;

    MZETA_PZETA_p2 = MZETA_PZETA*MZETA_PZETA;
    MZETA_10_PZETA_10_p2 = MZETA_10_PZETA_10*MZETA_10_PZETA_10;
    MZETA_01_PZETA_01_p2 = MZETA_01_PZETA_01*MZETA_01_PZETA_01;

    MZETA_PZETA_p3 = MZETA_PZETA_p2*MZETA_PZETA;
    MZETA_10_PZETA_10_p3 = MZETA_10_PZETA_10_p2*MZETA_10_PZETA_10;
    MZETA_01_PZETA_01_p3 = MZETA_01_PZETA_01_p2*MZETA_01_PZETA_01;

    MZETA_PZETA_p4 = MZETA_PZETA_p3*MZETA_PZETA;

/* Derivatives of EC */

/* EC_drv = pw92_d(roa,rob); */
    drv3_clear(&EC_drv);
    PW92Functional.third(&EC_drv, 1.0, dp);

    EC_10=EC_drv.df1000;
    EC_01=EC_drv.df0100;
    EC_20=EC_drv.df2000;
    EC_02=EC_drv.df0200;
    EC_11=EC_drv.df1100;
    EC_30=EC_drv.df3000;
    EC_03=EC_drv.df0300;
    EC_21=EC_drv.df2100;
    EC_12=EC_drv.df1200;

/* Powers of EC derivatives */

    EC_10_p2=EC_10*EC_10;
    EC_01_p2=EC_01*EC_01;
    EC_10_p3=EC_10_p2*EC_10;
    EC_01_p3=EC_01_p2*EC_01;

/* Derivatives of H0 */

    H0_10=(3.0/8.0)*gamma*MZETA_PZETA_p2*(MZETA_10_PZETA_10);

    H0_01=(3.0/8.0)*gamma*MZETA_PZETA_p2*(MZETA_01_PZETA_01);

    H0_20=(3.0/8.0)*gamma*MZETA_PZETA*(2.0*MZETA_10_PZETA_10_p2 + 
                                       MZETA_PZETA*MZETA_20_PZETA_20);

    H0_11=(3.0/8.0)*gamma*MZETA_PZETA*(2.0*MZETA_01_PZETA_01*MZETA_10_PZETA_10
                                       +   MZETA_PZETA*MZETA_11_PZETA_11);

    H0_02=(3.0/8.0)*gamma*MZETA_PZETA*(2.0*MZETA_01_PZETA_01_p2 + 
                                       MZETA_PZETA*MZETA_02_PZETA_02);

    H0_30=(3.0/8.0)*gamma*(2.0*MZETA_10_PZETA_10_p3 + 
                           6.0*MZETA_PZETA*MZETA_10_PZETA_10*MZETA_20_PZETA_20
                           +   MZETA_PZETA_p2*MZETA_30_PZETA_30);

    H0_21=(3.0/8.0)*gamma*(2*MZETA_01_PZETA_01*MZETA_10_PZETA_10_p2 + 
                           4*MZETA_PZETA*MZETA_10_PZETA_10*MZETA_11_PZETA_11+
                           2*MZETA_PZETA*MZETA_01_PZETA_01*MZETA_20_PZETA_20+
                           MZETA_PZETA_p2*MZETA_21_PZETA_21);

    H0_12=(3.0/8.0)*gamma*(2*MZETA_01_PZETA_01_p2*MZETA_10_PZETA_10 + 
                           4*MZETA_PZETA*MZETA_01_PZETA_01*MZETA_11_PZETA_11 + 
                           2*MZETA_PZETA*MZETA_02_PZETA_02*MZETA_10_PZETA_10 + 
                           MZETA_PZETA_p2*MZETA_12_PZETA_12);

    H0_03=(3.0/8.0)*gamma*(2*MZETA_01_PZETA_01_p3 +
                           6*MZETA_PZETA*MZETA_01_PZETA_01*MZETA_02_PZETA_02 + 
                           MZETA_PZETA_p2*MZETA_03_PZETA_03);

/* All terms OK up to here ! */

/* Powers of H0 derivatives */
    H0_10_p2=H0_10*H0_10;
    H0_01_p2=H0_01*H0_01;
    H0_10_p3=H0_10_p2*H0_10;
    H0_01_p3=H0_01_p2*H0_01;

/* Derivatives of T */

    T_1000=-(groa_grob*(7.0*MZETA+7.0*PZETA+6.0*roa_rob*(MZETA_10_PZETA_10)))
        /(6.0*nKS*roa_rob_p13f6*MZETA_PZETA_p2);

    T_0100=-(groa_grob*(7.0*MZETA+7.0*PZETA+6.0*roa_rob*(MZETA_01_PZETA_01)))
        /(6.0*nKS*roa_rob_p13f6*MZETA_PZETA_p2);

    T_0010=1.0/(nKS*roa_rob_p7f6*MZETA_PZETA);

    T_0001=T_0010;

    T_2000=(groa_grob*(91.0*MZETA_PZETA_p2+
                       84.0*roa_rob*(MZETA_PZETA)*(MZETA_10_PZETA_10)+
                       72.0*roa_rob_p2*MZETA_10_PZETA_10_p2-
                       36.0*roa_rob_p2*(MZETA_PZETA)*(MZETA_20_PZETA_20)))
        /(36.0*nKS*roa_rob_p19f6*MZETA_PZETA_p3);

    T_1100=(groa_grob*(91.0*MZETA_PZETA_p2+
                       42.0*roa_rob*(MZETA_PZETA)*(MZETA_01_PZETA_01)+
                       42.0*roa_rob*(MZETA_PZETA)*(MZETA_10_PZETA_10)+
                       72.0*roa_rob_p2*(MZETA_01_PZETA_01)*(MZETA_10_PZETA_10)-
                       36.0*roa_rob_p2*(MZETA_PZETA)*(MZETA_11_PZETA_11)))
        /(36.0*nKS*roa_rob_p19f6*MZETA_PZETA_p3);

    T_1010=(-7.0*MZETA-7.0*PZETA-6.0*roa_rob*(MZETA_10_PZETA_10))
        /(6.0*nKS*roa_rob_p13f6*MZETA_PZETA_p2);

    T_1001=T_1010;

    T_0200=(groa_grob*(91.0*MZETA_PZETA_p2+
                       84.0*roa_rob*(MZETA_PZETA)*(MZETA_01_PZETA_01)+
                       72.0*roa_rob_p2*MZETA_01_PZETA_01_p2-
                       36.0*roa_rob_p2*(MZETA_PZETA)*(MZETA_02_PZETA_02)))
        /(36.0*nKS*roa_rob_p19f6*MZETA_PZETA_p3);

    T_0110=(-7.0*MZETA-7.0*PZETA-6.0*roa_rob*(MZETA_01_PZETA_01))
        /(6.0*nKS*roa_rob_p13f6*MZETA_PZETA_p2);

    T_0101=T_0110;

    T_3000=(groa_grob*(-1729.0*MZETA_PZETA_p3
                       -1638.0*roa_rob*MZETA_PZETA_p2*(MZETA_10_PZETA_10)
                       -1512.0*roa_rob_p2*(MZETA_PZETA)*MZETA_10_PZETA_10_p2
                       -1296.0*roa_rob_p3*MZETA_10_PZETA_10_p3
                       +756.0*roa_rob_p2*MZETA_PZETA_p2*(MZETA_20_PZETA_20)
                       +1296.0*roa_rob_p3*(MZETA_PZETA)*(MZETA_10_PZETA_10)
                       *(MZETA_20_PZETA_20)
                       -216.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_30_PZETA_30)))
        /(216.0*nKS*roa_rob_p25f6*MZETA_PZETA_p4);

    T_2100=(groa_grob*(-1729.0*MZETA_PZETA_p3
                       -546.0*roa_rob*MZETA_PZETA_p2*(MZETA_01_PZETA_01)
                       -1092.0*roa_rob*MZETA_PZETA_p2*(MZETA_10_PZETA_10)
                       -1008.0*roa_rob_p2*(MZETA_PZETA)*(MZETA_01_PZETA_01)
                       *(MZETA_10_PZETA_10)
                       -504.0*roa_rob_p2*(MZETA_PZETA)*MZETA_10_PZETA_10_p2
                       -1296.0*roa_rob_p3*(MZETA_01_PZETA_01)
                       *MZETA_10_PZETA_10_p2
                       +504.0*roa_rob_p2*MZETA_PZETA_p2*(MZETA_11_PZETA_11)
                       +864.0*roa_rob_p3*(MZETA_PZETA)*(MZETA_10_PZETA_10)
                       *(MZETA_11_PZETA_11)
                       +252.0*roa_rob_p2*MZETA_PZETA_p2*(MZETA_20_PZETA_20)
                       +432.0*roa_rob_p3*(MZETA_PZETA)*(MZETA_01_PZETA_01)
                       *(MZETA_20_PZETA_20)
                       -216.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_21_PZETA_21)))
        /(216.0*nKS*roa_rob_p25f6*MZETA_PZETA_p4);

    T_2010=(91.0*MZETA_PZETA_p2+84.0*roa_rob*(MZETA_PZETA)*(MZETA_10_PZETA_10)+
            72.0*roa_rob_p2*MZETA_10_PZETA_10_p2-
            36.0*roa_rob_p2*(MZETA_PZETA)*(MZETA_20_PZETA_20))
        /(36.0*nKS*roa_rob_p19f6*MZETA_PZETA_p3);

    T_2001=T_2010;

    T_1110=(91.0*MZETA_PZETA_p2+42.0*roa_rob*(MZETA_PZETA)*(MZETA_01_PZETA_01)+
            42.0*roa_rob*(MZETA_PZETA)*(MZETA_10_PZETA_10)+
            72.0*roa_rob_p2*(MZETA_01_PZETA_01)*(MZETA_10_PZETA_10)-
            36.0*roa_rob_p2*(MZETA_PZETA)*(MZETA_11_PZETA_11))
        /(36.0*nKS*roa_rob_p19f6*MZETA_PZETA_p3);

    T_1101=T_1110;

    T_1200=(groa_grob*(-1729.0*MZETA_PZETA_p3-
                       1092.0*roa_rob*MZETA_PZETA_p2*(MZETA_01_PZETA_01)-
                       504.0*roa_rob_p2*(MZETA_PZETA)*MZETA_01_PZETA_01_p2+
                       252.0*roa_rob_p2*MZETA_PZETA_p2*(MZETA_02_PZETA_02)-
                       546.0*roa_rob*MZETA_PZETA_p2*(MZETA_10_PZETA_10)-
                       1008.0*roa_rob_p2*(MZETA_PZETA)*(MZETA_01_PZETA_01)*
                       (MZETA_10_PZETA_10)-
                       1296.0*roa_rob_p3*MZETA_01_PZETA_01_p2*
                       (MZETA_10_PZETA_10)+
                       432.0*roa_rob_p3*(MZETA_PZETA)*(MZETA_02_PZETA_02)*
                       (MZETA_10_PZETA_10)+
                       504.0*roa_rob_p2*MZETA_PZETA_p2*(MZETA_11_PZETA_11)+
                       864.0*roa_rob_p3*(MZETA_PZETA)*(MZETA_01_PZETA_01)*
                       (MZETA_11_PZETA_11)-
                       216.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_12_PZETA_12)))
        /(216.0*nKS*roa_rob_p25f6*MZETA_PZETA_p4);

    T_0210=(91.0*MZETA_PZETA_p2+84.0*roa_rob*(MZETA_PZETA)*(MZETA_01_PZETA_01)+
            72.0*roa_rob_p2*MZETA_01_PZETA_01_p2-
            36.0*roa_rob_p2*(MZETA_PZETA)*(MZETA_02_PZETA_02))
        /(36.0*nKS*roa_rob_p19f6*MZETA_PZETA_p3);

    T_0201=T_0210;

    T_0300=(groa_grob*(-1729.0*MZETA_PZETA_p3-
                       1638.0*roa_rob*MZETA_PZETA_p2*(MZETA_01_PZETA_01)-
                       1512.0*roa_rob_p2*(MZETA_PZETA)*MZETA_01_PZETA_01_p2-
                       1296.0*roa_rob_p3*MZETA_01_PZETA_01_p3+
                       756.0*roa_rob_p2*MZETA_PZETA_p2*(MZETA_02_PZETA_02)+
                       1296.0*roa_rob_p3*(MZETA_PZETA)*(MZETA_01_PZETA_01)*
                       (MZETA_02_PZETA_02)-
                       216.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_03_PZETA_03)))
        /(216.0*nKS*roa_rob_p25f6*MZETA_PZETA_p4);

/* All terms ok up to here! */
/* Powers of T derivatives */
    T_1000_p2 = T_1000*T_1000;
    T_1000_p3 = T_1000_p2*T_1000;

    T_0100_p2 = T_0100*T_0100;
    T_0100_p3 = T_0100_p2*T_0100;

    T_0010_p2 = T_0010*T_0010;
    T_0010_p3 = T_0010_p2*T_0010;

    T_0001_p2 = T_0001*T_0001;
    T_0001_p3 = T_0001_p2*T_0001;

/* Derivatives of A */

    A_10=(delta*AEXP*(H0*EC_10-EC*H0_10))/(ADENOM_p2*H0_p2);

    A_01=(delta*AEXP*(H0*EC_01-EC*H0_01))/(ADENOM_p2*H0_p2);

    A_20=(delta*AEXP*(-((1.0+AEXP)*EC_p2*H0_10_p2)
                      +2.0*EC*H0*H0_10*((1.0+AEXP)*EC_10+ADENOM*H0_10)
                      +ADENOM*H0_p3*EC_20
                      -H0_p2*((1.0+AEXP)*EC_10_p2+2.0*ADENOM*EC_10*H0_10
                              +ADENOM*EC*H0_20)))/(ADENOM_p3*H0_p4);

    A_02=(delta*AEXP*(-((1.0+AEXP)*EC_p2*H0_01_p2)
                      +2.0*EC*H0*H0_01*((1.0+AEXP)*EC_01+ADENOM*H0_01)
                      +ADENOM*H0_p3*EC_02
                      -H0_p2*((1.0+AEXP)*EC_01_p2+2.0*ADENOM*EC_01*H0_01
                              +ADENOM*EC*H0_02)))/(ADENOM_p3*H0_p4);

    A_11=(delta*AEXP*(-((1.0+AEXP)*EC_p2*H0_01*H0_10)
                      +EC*H0*((1.0+AEXP)*EC_01*H0_10
                              +H0_01*((1.0+AEXP)*EC_10+2.0*ADENOM*H0_10))
                      +ADENOM*H0_p3*EC_11
                      -H0_p2*(EC_01*((1.0+AEXP)*EC_10+ADENOM*H0_10)
                              +ADENOM*(H0_01*EC_10+EC*H0_11))))
        /(ADENOM_p3*H0_p4);

    A_30=(delta*AEXP*(-((1.0+4.0*AEXP+AEXP_p2)*EC_p3*H0_10_p3)
                      +3.0*EC_p2*H0*H0_10_p2*
                      ((1.0+4.0*AEXP+AEXP_p2)*EC_10+2.0*(-1.0+AEXP_p2)*H0_10)
                      -3.0*EC*H0_p2*H0_10*((1.0+4.0*AEXP+AEXP_p2)*EC_10_p2
                                           +4.0*(-1.0+AEXP_p2)*EC_10*H0_10
                                           +(ADENOM)*(2.0*(ADENOM)*H0_10_p2
                                                      +(1.0+AEXP)*EC*H0_20))
                      +H0_p3*((1.0+4.0*AEXP+AEXP_p2)*EC_10_p3+
                              6.0*(-1.0+AEXP_p2)*EC_10_p2*H0_10+
                              3.0*(ADENOM)*EC*H0_10*((1.0+AEXP)*EC_20+
                                                     2.0*(ADENOM)*H0_20)+
                              3.0*(ADENOM)*EC_10*(2.0*(ADENOM)*H0_10_p2+
                                                  (1.0+AEXP)*EC*H0_20))
                      +ADENOM_p2*H0_p5*EC_30
                      -ADENOM*H0_p4*(3.0*EC_10*
                                     ((1.0+AEXP)*EC_20+(ADENOM)*H0_20)
                                     +ADENOM*(3.0*H0_10*EC_20+EC*H0_30))))
        /(ADENOM_p4*H0_p6);

    A_21=(delta*AEXP*(-((1.0+4.0*AEXP+AEXP_p2)*EC_p3*H0_01*H0_10_p2)
                      +EC_p2*H0*H0_10*((1.0+4.0*AEXP+AEXP_p2)*EC_01*H0_10
                                       +2.0*H0_01*((1.0+4.0*AEXP+AEXP_p2)*EC_10
                                                   +3.0*(-1.0+AEXP_p2)*H0_10))
                      -EC*H0_p2*(2.0*H0_10*(EC_01*((1.0+4.0*AEXP+AEXP_p2)*EC_10
                                                   +2.0*(-1.0+AEXP_p2)*H0_10)
                                            +(-1.0+AEXP_p2)*EC*H0_11)
                                 +H0_01*((1.0+4.0*AEXP+AEXP_p2)*EC_10_p2
                                         +8.0*(-1.0+AEXP_p2)*EC_10*H0_10
                                         +ADENOM*(6.0*ADENOM*H0_10_p2
                                                  +(1.0+AEXP)*EC*H0_20)))
                      +H0_p3*(EC_01*((1.0+4.0*AEXP+AEXP_p2)*EC_10_p2
                                     +4.0*(-1.0+AEXP_p2)*EC_10*H0_10
                                     +ADENOM*(2.0*ADENOM*H0_10_p2+
                                              (1.0+AEXP)*EC*H0_20))+
                              ADENOM*(2.0*EC*((1.0+AEXP)*EC_10*H0_11
                                              +H0_10*((1.0+AEXP)*EC_11
                                                      +2.0*(ADENOM)*H0_11))
                                      +H0_01*(2.0*(1.0+AEXP)*EC_10_p2
                                              +4.0*(ADENOM)*EC_10*H0_10
                                              +EC*((1.0+AEXP)*EC_20
                                                   +2.0*(ADENOM)*H0_20))))
                      +ADENOM_p2*H0_p5*EC_21
                      -ADENOM*H0_p4*(2.0*(ADENOM)*H0_10*EC_11
                                     +2.0*EC_10*((1.0+AEXP)*EC_11
                                                 +ADENOM*H0_11)+
                                     EC_01*EC_20+AEXP*EC_01*EC_20
                                     -H0_01*EC_20+AEXP*H0_01*EC_20
                                     -EC_01*H0_20+AEXP*EC_01*H0_20
                                     -EC*H0_21+AEXP*EC*H0_21)))
        /(ADENOM_p4*H0_p6);

    A_12=(delta*AEXP*(-((1.0+4.0*AEXP+AEXP_p2)*EC_p3*H0_01_p2*H0_10)
                      +EC_p2*H0*H0_01*(2.0*(1.0+4.0*AEXP+AEXP_p2)*EC_01*H0_10
                                       +H0_01*((1.0+4.0*AEXP+AEXP_p2)*EC_10
                                               +6.0*(-1.0+AEXP_p2)*H0_10))
                      -EC*H0_p2*((1.0+4.0*AEXP+AEXP_p2)*EC_01_p2*H0_10
                                 +2.0*EC_01*H0_01*((1.0+4.0*AEXP+AEXP_p2)*EC_10
                                                   +4.0*(-1.0+AEXP_p2)*H0_10)
                                 +ADENOM*((1.0+AEXP)*EC*H0_02*H0_10
                                          +H0_01_p2*(4.0*(1.0+AEXP)*EC_10
                                                     +6.0*ADENOM*H0_10)
                                          +2.0*(1.0+AEXP)*EC*H0_01*H0_11))
                      +H0_p3*(EC_01_p2*((1.0+4.0*AEXP+AEXP_p2)*EC_10
                                        +2.0*(-1.0+AEXP_p2)*H0_10)
                              +2.0*ADENOM*EC_01*(2.0*H0_01*((1.0+AEXP)*EC_10
                                                            +ADENOM*H0_10)
                                                 +(1.0+AEXP)*EC*H0_11)
                              +ADENOM*(2.0*ADENOM*H0_01_p2*EC_10
                                       +EC*((1.0+AEXP)*EC_02*H0_10
                                            +H0_02*((1.0+AEXP)*EC_10
                                                    +2.0*ADENOM*H0_10))
                                       +2.0*EC*H0_01*((1.0+AEXP)*EC_11
                                                      +2.0*ADENOM*H0_11)))
                      +ADENOM_p2*H0_p5*EC_12
                      -ADENOM*H0_p4*(ADENOM*H0_02*EC_10
                                     +EC_02*((1.0+AEXP)*EC_10+ADENOM*H0_10)
                                     +2.0*EC_01*EC_11+2.0*AEXP*EC_01*EC_11
                                     -2.0*H0_01*EC_11+2.0*AEXP*H0_01*EC_11
                                     -2.0*EC_01*H0_11+2.0*AEXP*EC_01*H0_11
                                     -EC*H0_12+AEXP*EC*H0_12)))
        /(ADENOM_p4*H0_p6);


    A_03=(delta*AEXP*(-((1.0+4.0*AEXP+AEXP_p2)*EC_p3*H0_01_p3)
                      +3.0*EC_p2*H0*H0_01_p2*((1.0+4.0*AEXP+AEXP_p2)*EC_01
                                              +2.0*(-1.0+AEXP_p2)*H0_01)
                      -3.0*EC*H0_p2*H0_01*((1.0+4.0*AEXP+AEXP_p2)*EC_01_p2
                                           +4.0*(-1.0+AEXP_p2)*EC_01*H0_01
                                           +(ADENOM)*(2.0*(ADENOM)*H0_01_p2
                                                      +(1.0+AEXP)*EC*H0_02))
                      +H0_p3*((1.0+4.0*AEXP+AEXP_p2)*EC_01_p3
                              +6.0*(-1.0+AEXP_p2)*EC_01_p2*H0_01
                              +3.0*ADENOM*EC*H0_01*((1.0+AEXP)*EC_02
                                                    +2.0*(ADENOM)*H0_02)
                              +3.0*ADENOM*EC_01*(2.0*(ADENOM)*H0_01_p2
                                                   +(1.0+AEXP)*EC*H0_02))
                      +ADENOM_p2*H0_p5*EC_03
                      -ADENOM*H0_p4*(3.0*EC_01*((1.0+AEXP)*EC_02
                                                +(ADENOM)*H0_02)
                                     +(ADENOM)*(3.0*H0_01*EC_02+EC*H0_03))))
        /(ADENOM_p4*H0_p6);

/* Powers of A derivatives */
    A_10_p2 = A_10*A_10;
    A_10_p3 = A_10_p2*A_10;

    A_01_p2 = A_01*A_01;
    A_01_p3 = A_01_p2*A_01;


/* Derivatives of H with respect to H0, T, A */
/* Derivatives with respect to H0 higher then first order vanish */
/* Terms used in the derivatives */
    At2 = A*T_p2;
    At2_p2=At2*At2;

    At_nom = 1.0 + At2;
    At_denom = At_nom + At2_p2;
    delta_A = delta + A;
    At_term = (1.0+delta_A*T_p2*At_nom);

/* Powers of terms */
    At_denom_p2 = At_denom*At_denom;
    At_denom_p3 = At_denom_p2*At_denom;

    At_term_p2 = At_term*At_term;
    At_term_p3 = At_term_p2*At_term;

    delta_A_p2 = delta_A*delta_A;

/* First order derivatives */
    H_A =-((delta*A*H0*T_p6*(2.0+At2))/(At_denom*At_term));

    H_T =(2.0*delta*H0*T*(1.0+2.0*At2))/(At_denom*At_term);

    H_H = log(1.0+(delta-delta/At_denom)/A);

/* Second order derivatives */
    H_AA =(delta*H0*T_p6*(-2.0-2.0*delta_A*T_p2-2.0*(delta-2.0*A)*A*T_p4
                          +2.0*A_p2*(delta+4.0*A)*T_p6
                          +4.0*A_p3*(delta+2.0*A)*T_p8
                          +A_p4*(delta+2.0*A)*T_p10))/(At_denom_p2*At_term_p2);

    H_TT =(2.0*delta*H0*(1.0+T_p2*(-delta+A*(4.0+T_p2*
                                             (-4.0*delta+A*(-5.0-2.0*(7.0*delta+8.0*A)*T_p2-19.0*A*delta_A*T_p4-10.0*A_p2*delta_A*T_p6))))))
        /(At_denom_p2*At_term_p2);

    H_HT =(2.0*delta*T*(1.0+2.0*At2))/(At_denom*At_term);

    H_TA =(-2.0*delta*A*H0*T_p5*(6.0+4.0*(delta+3.0*A)*T_p2
                                 +A*(7.0*delta+12.0*A)*T_p4+
                                 2.0*A_p2*(2.0*delta+3.0*A)*T_p6))
        /(At_denom_p2*At_term_p2);

    H_HA =-((delta*A*T_p6*(2.0+At2))/(At_denom*At_term));

/* Third order derivatives */
    H_TTT =(4.0*delta*H0*T*
            (-3.0*delta+T_p2*(
                delta_p2+A*
                (6.0*delta*(-3.0+delta*T_p2)+
                 A*(-30.0-3.0*(23.0*delta+32.0*A)*T_p2+
                    6.0*(2.0*delta-21.0*A)*delta_A*T_p4
                    +A*(29.0*delta_p2-39.0*delta*A-66.0*A_p2)*T_p6
                    +6.0*A_p2*delta_A*(13.0*delta+11.0*A)*T_p8+
                    126.0*A_p3*delta_A_p2*T_p10
                    +96.0*A_p4*delta_A_p2*T_p12
                    +30.0*A_p5*delta_A_p2*T_p14)))))
        /(At_denom_p3*At_term_p3);

    H_AAA =(2.0*delta*H0*T_p8*
            (3.0+
             6.0*(delta+3.0*A)*T_p2+
             3.0*(delta_p2+10.0*delta*A+11.0*A_p2)*T_p4
             +3.0*A*(4.0*delta_p2+17.0*delta*A+10.0*A_p2)*T_p6
             +3.0*delta*A_p2*(6.0*delta+13.0*A)*T_p8
             +(7.0*delta_p2*A_p3-30.0*A_p5)*T_p10
             -3.0*A_p4*(2.0*delta_p2+8.0*delta*A+11.0*A_p2)*T_p12
             -6.0*A_p5*(delta_p2+3.0*A*delta_A)*T_p14
             -A_p6*(delta_p2+3.0*A*delta_A)*T_p16))
        /(At_denom_p3*At_term_p3);

    H_HTT =(2.0*delta*
            (1.0+T_p2*
             (-delta+A*(4.0+T_p2*(-4.0*delta+A*
                                  (-5.0-2.0*(7.0*delta+8.0*A)*T_p2
                                   -19.0*A*delta_A*T_p4
                                   -10.0*A_p2*delta_A*T_p6))))))
        /(At_denom_p2*At_term_p2);

    H_HAA =(delta*T_p6*(-2.0-2.0*delta_A*T_p2
                        -2.0*(delta-2.0*A)*A*T_p4
                        +2.0*A_p2*(delta+4.0*A)*T_p6
                        +4.0*A_p3*(delta+2.0*A)*T_p8
                        +A_p4*(delta+2.0*A)*T_p10))/(At_denom_p2*At_term_p2);

    H_TTA =(2.0*delta*A*H0*T_p4*
            (-30.0
             -2.0*(17.0*delta+48.0*A)*T_p2
             -3.0*(delta+2.0*A)*(4.0*delta+21.0*A)*T_p4
             -3.0*A*(9.0*delta_p2+18.0*delta*A+22.0*A_p2)*T_p6
             +A_p2*(-2.0*delta_p2+83.0*delta*A+66.0*A_p2)*T_p8
             +6.0*A_p3*(9.0*delta_p2+31.0*delta*A+21.0*A_p2)*T_p10
             +3.0*A_p4*delta_A*(19.0*delta+32.0*A)*T_p12
             +10.0*A_p5*delta_A*(2.0*delta+3.0*A)*T_p14))
        /(At_denom_p3*At_term_p3);
    
    H_TAA =(4.0*delta*H0*T_p5*
            (-3.0+T_p2*
             (-(delta*(5.0+2.0*delta*T_p2))+
              A*(-6.0+9.0*(-delta+A)*T_p2
                 -3.0*(delta_p2-2.0*A*(2.0*delta+7.0*A))*T_p4
                 +A*(6.0*delta_p2+52.0*delta*A+75.0*A_p2)*T_p6
                 +4.0*A_p2*(5.0*delta_p2+18.0*A*delta_A)*T_p8
                 +3.0*A_p3*(6.0*delta_p2+17.0*delta*A+14.0*A_p2)*T_p10
                 +2.0*A_p4*(3.0*delta_p2+8.0*delta*A+6.0*A_p2)*T_p12))))
        /(At_denom_p3*At_term_p3);

    H_HTA =(-2.0*delta*A*T_p5*(6.0+4.0*(delta+3.0*A)*T_p2
                               +A*(7.0*delta+12.0*A)*T_p4
                               +2.0*A_p2*(2.0*delta+3.0*A)*T_p6))
        /(At_denom_p2*At_term_p2);

/* Final expressions (general derivatives of F[H[roa,rob],
 * T[roa,rob,groa,grob], A[roa,rob]]) */

/* First order derivatives */
    ds->df1000 += factor*(
        A_10*H_A+H0_10*H_H+H_T*T_1000 + EC_10
        );

    ds->df0100 += factor*(
        A_01*H_A+H0_01*H_H+H_T*T_0100 + EC_01
        );

    ds->df0010 += factor*(
        H_T*T_0010
        );

    ds->df0001 += factor*(
        H_T*T_0001
        );

/* Second order derivatives */
    ds->df2000 += factor*(
        A_20*H_A+H0_20*H_H+A_10*(A_10*H_AA+2.0*H0_10*H_HA)
        +2.0*(A_10*H_TA+H0_10*H_HT)*T_1000+H_TT*T_1000_p2+H_T*T_2000 + EC_20
        );

    ds->df1100 += factor*(
        A_11*H_A+H0_11*H_H+A_01*(A_10*H_AA+H0_10*H_HA+H_TA*T_1000)
        +T_0100*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)
        +H0_01*(A_10*H_HA+H_HT*T_1000)+H_T*T_1100 + EC_11
        );

    ds->df1010 += factor*(
        T_0010*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+H_T*T_1010
        );

    ds->df1001 += factor*(
        T_0001*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+H_T*T_1001
        );

    ds->df0200 += factor*(
        A_02*H_A+H0_02*H_H+A_01*(A_01*H_AA+2.0*H0_01*H_HA)
        +2.0*(A_01*H_TA+H0_01*H_HT)*T_0100+H_TT*T_0100_p2+H_T*T_0200 + EC_02
        );

    ds->df0110 += factor*(
        T_0010*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)+H_T*T_0110
        );

    ds->df0101 += factor*(
        T_0001*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)+H_T*T_0101
        );

    ds->df0020 += factor*(
        H_TT*T_0010_p2
        );

    ds->df0011 += factor*(
        H_TT*T_0001*T_0010
        );

    ds->df0002 += factor*(
        H_TT*T_0001_p2
        );

/* Third order derivatives */
    ds->df3000 += factor*(
        A_30*H_A+A_10_p3*H_AAA+H0_30*H_H
        +3.0*A_10_p2*(H0_10*H_HAA+H_TAA*T_1000)
        +3.0*A_10*(A_20*H_AA+H0_20*H_HA+2.0*H0_10*H_HTA*T_1000+
                  H_TTA*T_1000_p2+H_TA*T_2000)
        +T_1000*(3.0*A_20*H_TA+3.0*H0_20*H_HT+H_TTT*T_1000_p2+3.0*H_TT*T_2000)
        +3.0*H0_10*(A_20*H_HA+H_HTT*T_1000_p2+H_HT*T_2000)+H_T*T_3000 + EC_30
        );

    ds->df2100 += factor*(
        A_21*H_A+H0_21*H_H+2.0*H0_10*A_11*H_HA+H0_01*A_20*H_HA
        +A_20*H_TA*T_0100+H0_20*H_HT*T_0100
        +A_10_p2*(A_01*H_AAA+H0_01*H_HAA+H_TAA*T_0100)
        +2.0*A_11*H_TA*T_1000+2.0*H0_11*H_HT*T_1000
        +2.0*H0_10*H_HTT*T_0100*T_1000+H0_01*H_HTT*T_1000_p2
        +H_TTT*T_0100*T_1000_p2+2.0*H0_10*H_HT*T_1100+2.0*H_TT*T_1000*T_1100
        +2.0*A_10*(A_11*H_AA+H0_11*H_HA+H0_10*(A_01*H_HAA+H_HTA*T_0100)
                   +(A_01*H_TAA+H0_01*H_HTA+H_TTA*T_0100)*T_1000+H_TA*T_1100)
        +(H0_01*H_HT+H_TT*T_0100)*T_2000
        +A_01*(A_20*H_AA+H0_20*H_HA+2.0*H0_10*H_HTA*T_1000+H_TTA*T_1000_p2
               +H_TA*T_2000)+H_T*T_2100 + EC_21
        );

    ds->df2010 += factor*(
        2.0*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)*T_1010
        +T_0010*(A_20*H_TA+H0_20*H_HT+A_10*(A_10*H_TAA+2.0*H0_10*H_HTA)
                 +2.0*(A_10*H_TTA+H0_10*H_HTT)*T_1000+H_TTT*T_1000_p2
                 +H_TT*T_2000)+H_T*T_2010
        );

    ds->df2001 += factor*(
        2.0*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)*T_1001
        +T_0001*(A_20*H_TA+H0_20*H_HT+A_10*(A_10*H_TAA+2.0*H0_10*H_HTA)
                 +2.0*(A_10*H_TTA+H0_10*H_HTT)*T_1000+H_TTT*T_1000_p2
                 +H_TT*T_2000)+H_T*T_2001
        );

    ds->df1200 += factor*(
        A_12*H_A+H0_12*H_H+A_11*(A_01*H_AA+2.0*H0_01*H_HA)
        +A_02*(A_10*H_AA+H0_10*H_HA+H_TA*T_1000)
        +T_0200*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)
        +T_0100_p2*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)
        +H0_02*(A_10*H_HA+H_HT*T_1000)
        +A_01*(A_11*H_AA+2.0*H0_11*H_HA
               +A_01*(A_10*H_AAA+H0_10*H_HAA+H_TAA*T_1000)
               +2.0*H0_01*(A_10*H_HAA+H_HTA*T_1000))
        +2.0*T_0100*(A_11*H_TA+H0_11*H_HT
                     +A_01*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)
                     +H0_01*(A_10*H_HTA+H_HTT*T_1000))
        +2.0*(A_01*H_TA+H0_01*H_HT)*T_1100+2.0*H_TT*T_0100*T_1100
        +H_T*T_1200 + EC_12
        );

    ds->df1110 += factor*(
        T_0110*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)
        +(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_1010
        +T_0010*(A_11*H_TA+H0_11*H_HT
                 +A_01*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)
                 +T_0100*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)
                 +H0_01*(A_10*H_HTA+H_HTT*T_1000)+H_TT*T_1100)+H_T*T_1110
        );

    ds->df1101 += factor*(
        T_0101*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)
        +(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_1001
        +T_0001*(A_11*H_TA+H0_11*H_HT
                 +A_01*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)
                 +T_0100*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)
                 +H0_01*(A_10*H_HTA+H_HTT*T_1000)+H_TT*T_1100)+H_T*T_1101
        );

    ds->df1020 += factor*(
        T_0010*(T_0010*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+2.0*H_TT*T_1010)
        );

    ds->df1011 += factor*(
        T_0010*(T_0001*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+H_TT*T_1001)
        +H_TT*T_0001*T_1010
        );

    ds->df1002 += factor*(
        T_0001*(T_0001*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+2.0*H_TT*T_1001)
        );

    ds->df0300 += factor*(
        A_03*H_A+A_01_p3*H_AAA+H0_03*H_H
        +3.0*A_01_p2*(H0_01*H_HAA+H_TAA*T_0100)
        +3.0*A_01*(A_02*H_AA+H0_02*H_HA+2.0*H0_01*H_HTA*T_0100
                   +H_TTA*T_0100_p2+H_TA*T_0200)
        +T_0100*(3.0*A_02*H_TA+3.0*H0_02*H_HT+H_TTT*T_0100_p2+3.0*H_TT*T_0200)
        +3.0*H0_01*(A_02*H_HA+H_HTT*T_0100_p2+H_HT*T_0200)+H_T*T_0300 + EC_03
        );

    ds->df0210 += factor*(
        2.0*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_0110
        +T_0010*(A_02*H_TA+H0_02*H_HT+A_01*(A_01*H_TAA+2.0*H0_01*H_HTA)
                 +2.0*(A_01*H_TTA+H0_01*H_HTT)*T_0100+H_TTT*T_0100_p2
                 +H_TT*T_0200)+H_T*T_0210
        );

    ds->df0201 += factor*(
        2.0*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_0101
        +T_0001*(A_02*H_TA+H0_02*H_HT+A_01*(A_01*H_TAA+2.0*H0_01*H_HTA)
                 +2.0*(A_01*H_TTA+H0_01*H_HTT)*T_0100+H_TTT*T_0100_p2
                 +H_TT*T_0200)+H_T*T_0201
        );

    ds->df0120 += factor*(
        T_0010*(T_0010*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)+2.0*H_TT*T_0110)
        );

    ds->df0111 += factor*(
        T_0010*(T_0001*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)+H_TT*T_0101)
        +H_TT*T_0001*T_0110
        );

    ds->df0102 += factor*(
        T_0001*(T_0001*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)+2.0*H_TT*T_0101)
        );

    ds->df0030 += factor*(
        H_TTT*T_0010_p3
        );

    ds->df0021 += factor*(
        H_TTT*T_0001*T_0010_p2
        );

    ds->df0012 += factor*(
        H_TTT*T_0001_p2*T_0010
        );

    ds->df0003 += factor*(
        H_TTT*T_0001_p3
        );

/* End (Finally!) */
}

#ifdef FOURTH_ORDER_DERIVATIVES
static void
pbe_cd4(FourthFuncDrv *ds, real factor, const DftDensProp *dp)
{

/* Declarations */
    real A;	real ADENOM;	real ADENOM_p2;	real ADENOM_p3;	real ADENOM_p4;
    real ADENOM_p5;	real AEXP;	real AEXP_p2;	real AEXP_p3;
    real A_01;	real A_01_p2;	real A_01_p3;	real A_01_p4;
    real A_02;	real A_02_p2;	real A_03;	real A_04;
    real A_10;	real A_10_p2;	real A_10_p3;	real A_10_p4;
    real A_11;	real A_12;	real A_13;	real A_20;
    real A_20_p2;	real A_21;	real A_22;	real A_30;
    real A_31;	real A_40;	real A_p10;	real A_p2;
    real A_p3;	real A_p4;	real A_p5;	real A_p6;
    real A_p7;	real A_p8;	real A_p9;	real At2;
    real At2_p2;	real At_denom;	real At_denom_p2;	real At_denom_p3;
    real At_denom_p4;	real At_nom;	real At_term;	real At_term_p2;
    real At_term_p3;	real At_term_p4;	real EC;	real EC_01;
    real EC_01_p2;	real EC_01_p3;	real EC_01_p4;	real EC_02;
    real EC_02_p2;	real EC_03;	real EC_04;	real EC_10;
    real EC_10_p2;	real EC_10_p3;	real EC_10_p4;	real EC_11;
    real EC_11_p2;	real EC_12;	real EC_13;	real EC_20;
    real EC_20_p2;	real EC_21;	real EC_22;	real EC_30;
    real EC_31;	real EC_40;	FourthFuncDrv EC_drv;	real EC_p2;
    real EC_p3;	real EC_p4;	real EC_p5;	real FI;
    real FI_p3;	real H;	real H0;	real H0_01;
    real H0_01_p2;	real H0_01_p3;	real H0_01_p4;	real H0_02;
    real H0_02_p2;	real H0_03;	real H0_04;	real H0_10;
    real H0_10_p2;	real H0_10_p3;	real H0_10_p4;	real H0_11;
    real H0_11_p2;	real H0_12;	real H0_13;	real H0_20;
    real H0_20_p2;	real H0_21;	real H0_22;	real H0_30;
    real H0_31;	real H0_40;	real H0_p2;	real H0_p3;
    real H0_p4;	real H0_p5;	real H0_p6;	real H0_p7;
    real H0_p8;	real H1;	real H_A;	real H_AA;
    real H_AAA;	real H_AAAA;	real H_H;	real H_HA;
    real H_HAA;	real H_HAAA;	real H_HATT;	real H_HT;
    real H_HTA;	real H_HTAA;	real H_HTT;	real H_HTTT;
    real H_T;	real H_TA;	real H_TAA;	real H_TAAA;
    real H_TT;	real H_TTA;	real H_TTAA;	real H_TTT;
    real H_TTTA;	real H_TTTT;	real KF;	real KS;
    real MZETA;	real MZETA_01;	real MZETA_01_PZETA_01;	real MZETA_01_PZETA_01_p2;
    real MZETA_01_PZETA_01_p3;	real MZETA_01_PZETA_01_p4;	real MZETA_02;	real MZETA_02_PZETA_02;
    real MZETA_02_PZETA_02_p2;	real MZETA_03;	real MZETA_03_PZETA_03;	real MZETA_04;
    real MZETA_04_PZETA_04;	real MZETA_10;	real MZETA_10_PZETA_10;	real MZETA_10_PZETA_10_p2;
    real MZETA_10_PZETA_10_p3;	real MZETA_10_PZETA_10_p4;	real MZETA_11;	real MZETA_11_PZETA_11;
    real MZETA_11_PZETA_11_p2;	real MZETA_12;	real MZETA_12_PZETA_12;	real MZETA_13;
    real MZETA_13_PZETA_13;	real MZETA_20;	real MZETA_20_PZETA_20;	real MZETA_20_PZETA_20_p2;
    real MZETA_21;	real MZETA_21_PZETA_21;	real MZETA_22;	real MZETA_22_PZETA_22;
    real MZETA_30;	real MZETA_30_PZETA_30;	real MZETA_31;	real MZETA_31_PZETA_31;
    real MZETA_40;	real MZETA_40_PZETA_13;	real MZETA_40_PZETA_31;	real MZETA_40_PZETA_40;
    real MZETA_PZETA;	real MZETA_PZETA_p2;	real MZETA_PZETA_p3;	real MZETA_PZETA_p4;
    real MZETA_PZETA_p5;	real MZETA_PZETA_p6;	real PZETA;	real PZETA_01;
    real PZETA_02;	real PZETA_03;	real PZETA_04;	real PZETA_10;
    real PZETA_11;	real PZETA_12;	real PZETA_13;	real PZETA_20;
    real PZETA_21;	real PZETA_22;	real PZETA_30;	real PZETA_31;
    real PZETA_40;	real T;	real T_0001;	real T_0001_p2;
    real T_0001_p3;	real T_0001_p4;	real T_0010;	real T_0010_p2;
    real T_0010_p3;	real T_0010_p4;	real T_0100;	real T_0100_p2;
    real T_0100_p3;	real T_0100_p4;	real T_0101;	real T_0101_p2;
    real T_0110;	real T_0110_p2;	real T_0200;	real T_0200_p2;
    real T_0201;	real T_0210;	real T_0300;	real T_0301;
    real T_0310;	real T_0400;	real T_1000;	real T_1000_p2;
    real T_1000_p3;	real T_1000_p4;	real T_1001;	real T_1001_p2;
    real T_1010;	real T_1010_p2;	real T_1100;	real T_1100_p2;
    real T_1101;	real T_1110;	real T_1200;	real T_1201;
    real T_1210;	real T_1300;	real T_2000;	real T_2000_p2;
    real T_2001;	real T_2010;	real T_2100;	real T_2101;
    real T_2110;	real T_2200;	real T_3000;	real T_3001;
    real T_3010;	real T_3100;	real T_4000;	real T_p10;
    real T_p12;	real T_p14;	real T_p16;	real T_p18;
    real T_p2;	real T_p20;	real T_p22;	real T_p24;
    real T_p3;	real T_p4;	real T_p5;	real T_p6;
    real T_p7;	real T_p8;	real ZETA;
    real delta_A;	real delta_A_p2;	real delta_A_p3;
    real groa;	real groa_grob;	real grob;
    real roa;	real roa_p2;
    real roa_p3;	real roa_p4;	real roa_rob;	real roa_rob_p13f6;
    real roa_rob_p19f6;	real roa_rob_p1f3;	real roa_rob_p1f6;	real roa_rob_p2;
    real roa_rob_p25f6;	real roa_rob_p3;	real roa_rob_p31f6;	real roa_rob_p4;
    real roa_rob_p5;	real roa_rob_p6;	real roa_rob_p7;	real roa_rob_p7f6;
    real roafro;		real roafro_p1f3;	real roafro_p2f3;	real roafro_p4f3;
    real roafro_p5f3;	real roafro_p7f3;	real rob;		real rob_p2;
    real rob_p3;	real rob_p4;	real robfro;	real robfro_p1f3;
    real robfro_p2f3;	real robfro_p4f3;	real robfro_p5f3;	real robfro_p7f3;
/* End of Declarations */

	
/* Constants & Numbers */
#include "pbe_c.h"
    const real delta_p2 = delta*delta;
    const real delta_p3 = delta_p2*delta;

/* Set up roa, rob, groa, grob */
    roa = dp->rhoa;
    rob = dp->rhob;

    groa = dp->grada;
    grob = dp->gradb;

/* roa, rob and powers */
    roa_rob=roa+rob;
    roa_rob_p1f3=pow(roa_rob,(1.0/3.0));
    roa_rob_p1f6=pow(roa_rob,(1.0/6.0));

    roa_rob_p2=roa_rob*roa_rob;
    roa_rob_p3=roa_rob_p2*roa_rob;
    roa_rob_p4=roa_rob_p2*roa_rob_p2;
    roa_rob_p5=roa_rob_p4*roa_rob;
    roa_rob_p6=roa_rob_p4*roa_rob_p2;
    roa_rob_p7=roa_rob_p6*roa_rob;

    roa_rob_p7f6 = roa_rob*roa_rob_p1f6;
    roa_rob_p13f6 = roa_rob_p2*roa_rob_p1f6;
    roa_rob_p19f6 = roa_rob_p3*roa_rob_p1f6;
    roa_rob_p25f6 = roa_rob_p4*roa_rob_p1f6;
    roa_rob_p31f6 = roa_rob_p5*roa_rob_p1f6;

    roa_p2=roa*roa;
    rob_p2=rob*rob;

    roa_p3=roa_p2*roa;
    rob_p3=rob_p2*rob;

    roa_p4=roa_p3*roa;
    rob_p4=rob_p3*rob;

/* groa, grob gradients and powers */
    groa_grob = groa+ grob;

/* Auxaliary functions */
    ZETA=(roa-rob)/(roa+rob);
    KF=nKF*roa_rob_p1f3;
    KS=nKS*roa_rob_p1f6;
    MZETA=pow((1.0-ZETA),(2.0/3.0));
    PZETA=pow((1.0+ZETA),(2.0/3.0));
    FI=(PZETA + MZETA)/2.0;
    T=(groa+grob)/((roa+rob)*(2.0*KS*FI));


/* Powers of Fi */
    FI_p3=FI*FI*FI;

/* Powers of T */
    T_p2=T*T;
    T_p3=T_p2*T;
    T_p4=T_p2*T_p2;
    T_p5=T_p4*T;
    T_p6=T_p4*T_p2;
    T_p7=T_p6*T;
    T_p8=T_p6*T_p2;
    T_p10=T_p8*T_p2;
    T_p12=T_p10*T_p2;
    T_p14=T_p12*T_p2;
    T_p16=T_p14*T_p2;
    T_p18=T_p16*T_p2;
    T_p20=T_p18*T_p2;
    T_p22=T_p20*T_p2;
    T_p24=T_p22*T_p2;

/* The pw92 EC(electron gas corelation) energy */
    EC = PW92Functional.func(dp);
    EC_p2 = EC*EC;
    EC_p3 = EC_p2*EC;
    EC_p4 = EC_p3*EC;
    EC_p5 = EC_p4*EC;

/* H0 and its powers */
    H0=gamma*FI_p3;
    H0_p2 = H0*H0;
    H0_p3 = H0_p2*H0;
    H0_p4 = H0_p3*H0;
    H0_p5 = H0_p4*H0;
    H0_p6 = H0_p5*H0;
    H0_p7 = H0_p6*H0;
    H0_p8 = H0_p7*H0;

/* Calc of A */
    AEXP=exp(EC/H0);
    AEXP_p2=AEXP*AEXP;
    AEXP_p3=AEXP_p2*AEXP;
    ADENOM=AEXP -1.0;
    ADENOM_p2=ADENOM*ADENOM;
    ADENOM_p3=ADENOM_p2*ADENOM;
    ADENOM_p4=ADENOM_p3*ADENOM;
    ADENOM_p5=ADENOM_p4*ADENOM;
    A=delta*( 1.0/((1.0/AEXP) -1.0)  ) ;

/* Powers of A */
    A_p2=A*A;
    A_p3=A_p2*A;
    A_p4=A_p3*A;
    A_p5=A_p4*A;
    A_p6=A_p5*A;
    A_p7=A_p6*A;
    A_p8=A_p7*A;
    A_p9=A_p8*A;
    A_p10=A_p9*A;

/* The pbe Ec (pbe correlation contribution) */
/* H is gga contribution to correlation energy */ 
    H1=1.0 + delta*T_p2*((1.0 + A*T_p2)/(1.0 + A*T_p2 + A_p2*T_p4));

    H=H0*log( H1 );
/* E=EC+H; */


/* Now start with derivatives *
 * first, common terms used in equations */

/* roafro, robfro and powers */
    roafro=roa/roa_rob;
    robfro=rob/roa_rob;

    roafro_p1f3=pow(roafro,(1.0/3.0));
    robfro_p1f3=pow(robfro,(1.0/3.0));

    roafro_p2f3=roafro_p1f3*roafro_p1f3;
    robfro_p2f3=robfro_p1f3*robfro_p1f3;

    roafro_p4f3=roafro_p1f3*roafro;
    robfro_p4f3=robfro_p1f3*robfro;

    roafro_p5f3=roafro_p2f3*roafro;
    robfro_p5f3=robfro_p2f3*robfro;

    roafro_p7f3=roafro_p4f3*roafro;
    robfro_p7f3=robfro_p4f3*robfro;

/* Derivatives of PZETA, MZETA */

    MZETA_10=(-2.0/3.0)*nZ*robfro_p5f3/rob;
    PZETA_10=(2.0/3.0)*nZ*rob/(roafro_p1f3*roa_rob_p2);

    MZETA_01=(2.0/3.0)*nZ*roa/(robfro_p1f3*roa_rob_p2);
    PZETA_01=(-2.0/3.0)*nZ*roafro_p5f3/roa;

    MZETA_20=(10.0/9.0)*nZ*robfro_p2f3/roa_rob_p2;
    PZETA_20=(-2.0/9.0)*nZ*rob*(6.0*roa + rob)/(roafro_p4f3*roa_rob_p4);

    MZETA_11=(2.0/9.0)*nZ*(-2.0*roa+3.0*rob)/(robfro_p1f3*roa_rob_p3);
    PZETA_11=(2.0/9.0)*nZ*(3.0*roa-2.0*rob)/(roafro_p1f3*roa_rob_p3);

    MZETA_02=(-2.0/9.0)*nZ*roa*(6.0*rob + roa)/(robfro_p4f3*roa_rob_p4);
    PZETA_02=(10.0/9.0)*nZ*roafro_p2f3/roa_rob_p2;

    MZETA_30=(-80.0/27.0)*nZ*robfro_p2f3/roa_rob_p3;
    PZETA_30=(4.0/27.0)*nZ*rob*(27.0*roa_p2 + 9.0*roa*rob+2.0*rob_p2)/(roafro_p7f3*roa_rob_p6);

    MZETA_21=(20.0/27.0)*nZ*(roa-3.0*rob)/(robfro_p1f3*roa_rob_p4);
    PZETA_21=(4.0/27.0)*nZ*(-9.0*roa_p2 + 12.0*roa*rob+ rob_p2)/(roafro_p4f3*roa_rob_p5);

    MZETA_12=(4.0/27.0)*nZ*(roa_p2 + 12.0*roa*rob -9.0*rob_p2)/(robfro_p4f3*roa_rob_p5);
    PZETA_12=(20.0/27.0)*nZ*(-3.0*roa+rob)/(roafro_p1f3*roa_rob_p4);

    MZETA_03=(4.0/27.0)*nZ*roa*(2.0*roa_p2 + 9.0*roa*rob +27.0*rob_p2)/(robfro_p7f3*roa_rob_p6);
    PZETA_03=(-80.0/27.0)*nZ*roafro_p2f3/roa_rob_p3;

    MZETA_40=(880.0/81.0)*nZ*robfro_p2f3/roa_rob_p4;
    PZETA_40=(-8.0/81.0)*nZ*rob*roafro_p2f3*(162.0*roa_p3+81.0*roa_p2*rob + 36.0*roa*rob_p2 + 7.0*rob_p3)/(roa_p4*roa_rob_p4);

    MZETA_31=(80.0/81.0)*nZ*(-2.0*roa + 9.0*rob)/(robfro_p1f3*roa_rob_p5);
    PZETA_31=(4.0/81.0)*nZ*(81.0*roa_p3-162.0*roa_p2*rob - 27.0*roa*rob_p2 - 4.0*rob_p3)/(roafro_p7f3*roa_rob_p7);

    MZETA_22=(-20.0/81.0)*nZ*(roa_p2 + 18.0*roa*rob - 27.0*rob_p2)/(robfro_p4f3*roa_rob_p6);
    PZETA_22=(20.0/81.0)*nZ*(27.0*roa_p2 - 18.0*roa*rob -rob_p2)/(roafro_p4f3*roa_rob_p6);

    MZETA_13=(-4.0/81.0)*nZ*(4.0*roa_p3 + 27.0*roa_p2*rob + 162.0*roa*rob_p2 -81.0*rob_p3)/(robfro_p7f3*roa_rob_p7);
    PZETA_13=(80.0/81.0)*nZ*(9.0*roa - 2.0*rob)/(roafro_p1f3*roa_rob_p5);

    MZETA_04=(-8.0/81.0)*nZ*roa*robfro_p2f3*(7.0*roa_p3+36.0*roa_p2*rob + 81.0*roa*rob_p2 + 162.0*rob_p3)/(rob_p4*roa_rob_p4);
    PZETA_04=(880.0/81.0)*nZ*roafro_p2f3/roa_rob_p4;

/* MZETA, PZETA terms and derivatives */
    MZETA_PZETA = MZETA + PZETA;
    MZETA_10_PZETA_10 = MZETA_10 + PZETA_10;
    MZETA_01_PZETA_01 = MZETA_01 + PZETA_01;

    MZETA_20_PZETA_20 = MZETA_20 + PZETA_20;
    MZETA_02_PZETA_02 = MZETA_02 + PZETA_02;
    MZETA_11_PZETA_11 = MZETA_11 + PZETA_11;

    MZETA_30_PZETA_30 = MZETA_30 + PZETA_30;
    MZETA_03_PZETA_03 = MZETA_03 + PZETA_03;
    MZETA_21_PZETA_21 = MZETA_21 + PZETA_21;
    MZETA_12_PZETA_12 = MZETA_12 + PZETA_12;

    MZETA_40_PZETA_40 = MZETA_40 + PZETA_40;
    MZETA_04_PZETA_04 = MZETA_04 + PZETA_04;
    MZETA_40_PZETA_31 = MZETA_31 + PZETA_31;
    MZETA_40_PZETA_13 = MZETA_13 + PZETA_13;
    MZETA_31_PZETA_31 = MZETA_31 + PZETA_31;
    MZETA_13_PZETA_13 = MZETA_13 + PZETA_13;
    MZETA_22_PZETA_22 = MZETA_22 + PZETA_22;

    MZETA_PZETA_p2 = MZETA_PZETA*MZETA_PZETA;
    MZETA_10_PZETA_10_p2 = MZETA_10_PZETA_10*MZETA_10_PZETA_10;
    MZETA_01_PZETA_01_p2 = MZETA_01_PZETA_01*MZETA_01_PZETA_01;

    MZETA_PZETA_p3 = MZETA_PZETA_p2*MZETA_PZETA;
    MZETA_10_PZETA_10_p3 = MZETA_10_PZETA_10_p2*MZETA_10_PZETA_10;
    MZETA_01_PZETA_01_p3 = MZETA_01_PZETA_01_p2*MZETA_01_PZETA_01;

    MZETA_PZETA_p4 = MZETA_PZETA_p3*MZETA_PZETA;
    MZETA_10_PZETA_10_p4 = MZETA_10_PZETA_10_p3*MZETA_10_PZETA_10;
    MZETA_01_PZETA_01_p4 = MZETA_01_PZETA_01_p3*MZETA_01_PZETA_01;

    MZETA_20_PZETA_20_p2 = MZETA_20_PZETA_20*MZETA_20_PZETA_20;
    MZETA_02_PZETA_02_p2 = MZETA_02_PZETA_02*MZETA_02_PZETA_02;
    MZETA_11_PZETA_11_p2 = MZETA_11_PZETA_11*MZETA_11_PZETA_11;

    MZETA_PZETA_p5 = MZETA_PZETA_p4*MZETA_PZETA;
    MZETA_PZETA_p6 = MZETA_PZETA_p5*MZETA_PZETA;

/* Derivatives of EC */

/* EC_drv = pw92_d(roa,rob); */
    drv4_clear(&EC_drv);
    PW92Functional.fourth(&EC_drv, 1.0, dp);

    EC_10=EC_drv.df10000;
    EC_01=EC_drv.df01000;
    EC_20=EC_drv.df20000;
    EC_02=EC_drv.df02000;
    EC_11=EC_drv.df11000;
    EC_30=EC_drv.df30000;
    EC_03=EC_drv.df03000;
    EC_21=EC_drv.df21000;
    EC_12=EC_drv.df12000;
    EC_40=EC_drv.df40000;
    EC_04=EC_drv.df04000;
    EC_31=EC_drv.df31000;
    EC_13=EC_drv.df13000;
    EC_22=EC_drv.df22000;

/* Powers of EC derivatives */

    EC_10_p2=EC_10*EC_10;
    EC_01_p2=EC_01*EC_01;
    EC_10_p3=EC_10_p2*EC_10;
    EC_01_p3=EC_01_p2*EC_01;
    EC_10_p4=EC_10_p3*EC_10;
    EC_01_p4=EC_01_p3*EC_01;

    EC_20_p2=EC_20*EC_20;
    EC_02_p2=EC_02*EC_02;

    EC_11_p2=EC_11*EC_11;

/* Derivatives of H0 */

    H0_10=(3.0/8.0)*gamma*MZETA_PZETA_p2*(MZETA_10_PZETA_10);

    H0_01=(3.0/8.0)*gamma*MZETA_PZETA_p2*(MZETA_01_PZETA_01);

    H0_20=(3.0/8.0)*gamma*MZETA_PZETA*(2.0*MZETA_10_PZETA_10_p2 + 
                                       MZETA_PZETA*MZETA_20_PZETA_20);

    H0_11=(3.0/8.0)*gamma*MZETA_PZETA*(2.0*MZETA_01_PZETA_01*MZETA_10_PZETA_10 + 
                                       MZETA_PZETA*MZETA_11_PZETA_11);

    H0_02=(3.0/8.0)*gamma*MZETA_PZETA*(2.0*MZETA_01_PZETA_01_p2 + 
                                       MZETA_PZETA*MZETA_02_PZETA_02);

    H0_30=(3.0/8.0)*gamma*(2.0*MZETA_10_PZETA_10_p3 + 6.0*MZETA_PZETA*MZETA_10_PZETA_10*MZETA_20_PZETA_20 + 
                           MZETA_PZETA_p2*MZETA_30_PZETA_30);

    H0_21=(3.0/8.0)*gamma*( 2.0*MZETA_01_PZETA_01*MZETA_10_PZETA_10_p2 + 
                            4.0*MZETA_PZETA*MZETA_10_PZETA_10*MZETA_11_PZETA_11 + 
                            2.0*MZETA_PZETA*MZETA_01_PZETA_01*MZETA_20_PZETA_20 + 
                            MZETA_PZETA_p2*MZETA_21_PZETA_21);

    H0_12=(3.0/8.0)*gamma*( 2.0*MZETA_01_PZETA_01_p2*MZETA_10_PZETA_10 + 
                            4.0*MZETA_PZETA*MZETA_01_PZETA_01*MZETA_11_PZETA_11 + 
                            2.0*MZETA_PZETA*MZETA_02_PZETA_02*MZETA_10_PZETA_10 + 
                            MZETA_PZETA_p2*MZETA_12_PZETA_12);

    H0_03=(3.0/8.0)*gamma*(2.0*MZETA_01_PZETA_01_p3 + 6.0*MZETA_PZETA*MZETA_01_PZETA_01*MZETA_02_PZETA_02 + 
                           MZETA_PZETA_p2*MZETA_03_PZETA_03);

    H0_40=(3.0/8.0)*gamma*(12.0*MZETA_10_PZETA_10_p2*MZETA_20_PZETA_20 + 
                           6.0*MZETA_PZETA*MZETA_20_PZETA_20_p2 + 
                           8.0*MZETA_PZETA*MZETA_10_PZETA_10*MZETA_30_PZETA_30 + 
                           MZETA_PZETA_p2*MZETA_40_PZETA_40 );

    H0_31=(3.0/8.0)*gamma*(6.0*MZETA_10_PZETA_10_p2*MZETA_11_PZETA_11 + 
                           6.0*MZETA_01_PZETA_01*MZETA_10_PZETA_10*MZETA_20_PZETA_20 + 
                           6.0*MZETA_PZETA*MZETA_11_PZETA_11*MZETA_20_PZETA_20 + 
                           6.0*MZETA_PZETA*MZETA_10_PZETA_10*MZETA_21_PZETA_21 + 
                           2.0*MZETA_PZETA*MZETA_01_PZETA_01*MZETA_30_PZETA_30 + 
                           MZETA_PZETA_p2*MZETA_31_PZETA_31);

    H0_22=(3.0/8.0)*gamma*(2.0*MZETA_02_PZETA_02*MZETA_10_PZETA_10_p2 + 
                           8.0*MZETA_01_PZETA_01*MZETA_10_PZETA_10*MZETA_11_PZETA_11 + 
                           4.0*MZETA_PZETA*MZETA_11_PZETA_11_p2 + 
                           4.0*MZETA_PZETA*MZETA_10_PZETA_10*MZETA_12_PZETA_12 + 
                           2.0*MZETA_01_PZETA_01_p2*MZETA_20_PZETA_20 + 
                           2.0*MZETA_PZETA*MZETA_02_PZETA_02*MZETA_20_PZETA_20 + 
                           4.0*MZETA_PZETA*MZETA_01_PZETA_01*MZETA_21_PZETA_21 + 
                           MZETA_PZETA_p2*MZETA_22_PZETA_22);

    H0_13=(3.0/8.0)*gamma*(6.0*MZETA_01_PZETA_01*MZETA_02_PZETA_02*MZETA_10_PZETA_10 + 
                           2.0*MZETA_PZETA*MZETA_03_PZETA_03*MZETA_10_PZETA_10 + 
                           6.0*MZETA_01_PZETA_01_p2*MZETA_11_PZETA_11 + 
                           6.0*MZETA_PZETA*MZETA_02_PZETA_02*MZETA_11_PZETA_11 + 
                           6.0*MZETA_PZETA*MZETA_01_PZETA_01*MZETA_12_PZETA_12 + 
                           MZETA_PZETA_p2*MZETA_13_PZETA_13);

    H0_04=(3.0/8.0)*gamma*(12.0*MZETA_01_PZETA_01_p2*MZETA_02_PZETA_02 + 
                           6.0*MZETA_PZETA*MZETA_02_PZETA_02_p2 + 
                           8.0*MZETA_PZETA*MZETA_01_PZETA_01*MZETA_03_PZETA_03 + 
                           MZETA_PZETA_p2*MZETA_04_PZETA_04 );

/* All terms OK up to here ! */

/* Powers of H0 derivatives */
    H0_10_p2=H0_10*H0_10;
    H0_01_p2=H0_01*H0_01;
    H0_10_p3=H0_10_p2*H0_10;
    H0_01_p3=H0_01_p2*H0_01;
    H0_10_p4=H0_10_p3*H0_10;
    H0_01_p4=H0_01_p3*H0_01;

    H0_20_p2=H0_20*H0_20;
    H0_02_p2=H0_02*H0_02;

    H0_11_p2=H0_11*H0_11;

/* Derivatives of T */

    T_1000=-(groa_grob*(7.0*MZETA+7.0*PZETA+6.0*roa_rob*(MZETA_10_PZETA_10)))/(6.0*nKS*roa_rob_p13f6*MZETA_PZETA_p2);

    T_0100=-(groa_grob*(7.0*MZETA+7.0*PZETA+6.0*roa_rob*(MZETA_01_PZETA_01)))/(6.0*nKS*roa_rob_p13f6*MZETA_PZETA_p2);

    T_0010=1.0/(nKS*roa_rob_p7f6*MZETA_PZETA);

    T_0001=T_0010;

    T_2000=(groa_grob*(91.0*MZETA_PZETA_p2+84.0*roa_rob*(MZETA_PZETA)*(MZETA_10_PZETA_10)+72.0*roa_rob_p2*MZETA_10_PZETA_10_p2-36.0*roa_rob_p2*(MZETA_PZETA)*(MZETA_20_PZETA_20)))/(36.0*nKS*roa_rob_p19f6*MZETA_PZETA_p3);

    T_1100=(groa_grob*(91.0*MZETA_PZETA_p2+42.0*roa_rob*(MZETA_PZETA)*(MZETA_01_PZETA_01)+42.0*roa_rob*(MZETA_PZETA)*(MZETA_10_PZETA_10)+72.0*roa_rob_p2*(MZETA_01_PZETA_01)*(MZETA_10_PZETA_10)-36.0*roa_rob_p2*(MZETA_PZETA)*(MZETA_11_PZETA_11)))/(36.0*nKS*roa_rob_p19f6*MZETA_PZETA_p3);

    T_1010=(-7.0*MZETA-7.0*PZETA-6.0*roa_rob*(MZETA_10_PZETA_10))/(6.0*nKS*roa_rob_p13f6*MZETA_PZETA_p2);

    T_1001=T_1010;

    T_0200=(groa_grob*(91.0*MZETA_PZETA_p2+84.0*roa_rob*(MZETA_PZETA)*(MZETA_01_PZETA_01)+72.0*roa_rob_p2*MZETA_01_PZETA_01_p2-36.0*roa_rob_p2*(MZETA_PZETA)*(MZETA_02_PZETA_02)))/(36.0*nKS*roa_rob_p19f6*MZETA_PZETA_p3);

    T_0110=(-7.0*MZETA-7.0*PZETA-6.0*roa_rob*(MZETA_01_PZETA_01))/(6.0*nKS*roa_rob_p13f6*MZETA_PZETA_p2);

    T_0101=T_0110;

    T_3000=(groa_grob*(-1729.0*MZETA_PZETA_p3-1638.0*roa_rob*MZETA_PZETA_p2*(MZETA_10_PZETA_10)-1512.0*roa_rob_p2*(MZETA_PZETA)*MZETA_10_PZETA_10_p2-1296.0*roa_rob_p3*MZETA_10_PZETA_10_p3+756.0*roa_rob_p2*MZETA_PZETA_p2*(MZETA_20_PZETA_20)+1296.0*roa_rob_p3*(MZETA_PZETA)*(MZETA_10_PZETA_10)*(MZETA_20_PZETA_20)-216.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_30_PZETA_30)))/(216.0*nKS*roa_rob_p25f6*MZETA_PZETA_p4);

    T_2100=(groa_grob*(-1729.0*MZETA_PZETA_p3-546.0*roa_rob*MZETA_PZETA_p2*(MZETA_01_PZETA_01)-1092.0*roa_rob*MZETA_PZETA_p2*(MZETA_10_PZETA_10)-1008.0*roa_rob_p2*(MZETA_PZETA)*(MZETA_01_PZETA_01)*(MZETA_10_PZETA_10)-504.0*roa_rob_p2*(MZETA_PZETA)*MZETA_10_PZETA_10_p2-1296.0*roa_rob_p3*(MZETA_01_PZETA_01)*MZETA_10_PZETA_10_p2+504.0*roa_rob_p2*MZETA_PZETA_p2*(MZETA_11_PZETA_11)+864.0*roa_rob_p3*(MZETA_PZETA)*(MZETA_10_PZETA_10)*(MZETA_11_PZETA_11)+252.0*roa_rob_p2*MZETA_PZETA_p2*(MZETA_20_PZETA_20)+432.0*roa_rob_p3*(MZETA_PZETA)*(MZETA_01_PZETA_01)*(MZETA_20_PZETA_20)-216.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_21_PZETA_21)))/(216.0*nKS*roa_rob_p25f6*MZETA_PZETA_p4);

    T_2010=(91.0*MZETA_PZETA_p2+84.0*roa_rob*(MZETA_PZETA)*(MZETA_10_PZETA_10)+72.0*roa_rob_p2*MZETA_10_PZETA_10_p2-36.0*roa_rob_p2*(MZETA_PZETA)*(MZETA_20_PZETA_20))/(36.0*nKS*roa_rob_p19f6*MZETA_PZETA_p3);

    T_2001=T_2010;

    T_1110=(91.0*MZETA_PZETA_p2+42.0*roa_rob*(MZETA_PZETA)*(MZETA_01_PZETA_01)+42.0*roa_rob*(MZETA_PZETA)*(MZETA_10_PZETA_10)+72.0*roa_rob_p2*(MZETA_01_PZETA_01)*(MZETA_10_PZETA_10)-36.0*roa_rob_p2*(MZETA_PZETA)*(MZETA_11_PZETA_11))/(36.0*nKS*roa_rob_p19f6*MZETA_PZETA_p3);

    T_1101=T_1110;

    T_1200=(groa_grob*(-1729.0*MZETA_PZETA_p3-1092.0*roa_rob*MZETA_PZETA_p2*(MZETA_01_PZETA_01)-504.0*roa_rob_p2*(MZETA_PZETA)*MZETA_01_PZETA_01_p2+252.0*roa_rob_p2*MZETA_PZETA_p2*(MZETA_02_PZETA_02)-546.0*roa_rob*MZETA_PZETA_p2*(MZETA_10_PZETA_10)-1008.0*roa_rob_p2*(MZETA_PZETA)*(MZETA_01_PZETA_01)*(MZETA_10_PZETA_10)-1296.0*roa_rob_p3*MZETA_01_PZETA_01_p2*(MZETA_10_PZETA_10)+432.0*roa_rob_p3*(MZETA_PZETA)*(MZETA_02_PZETA_02)*(MZETA_10_PZETA_10)+504.0*roa_rob_p2*MZETA_PZETA_p2*(MZETA_11_PZETA_11)+864.0*roa_rob_p3*(MZETA_PZETA)*(MZETA_01_PZETA_01)*(MZETA_11_PZETA_11)-216.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_12_PZETA_12)))/(216.0*nKS*roa_rob_p25f6*MZETA_PZETA_p4);

    T_0210=(91.0*MZETA_PZETA_p2+84.0*roa_rob*(MZETA_PZETA)*(MZETA_01_PZETA_01)+72.0*roa_rob_p2*MZETA_01_PZETA_01_p2-36.0*roa_rob_p2*(MZETA_PZETA)*(MZETA_02_PZETA_02))/(36.0*nKS*roa_rob_p19f6*MZETA_PZETA_p3);

    T_0201=T_0210;

    T_0300=(groa_grob*(-1729.0*MZETA_PZETA_p3-1638.0*roa_rob*MZETA_PZETA_p2*(MZETA_01_PZETA_01)-1512.0*roa_rob_p2*(MZETA_PZETA)*MZETA_01_PZETA_01_p2-1296.0*roa_rob_p3*MZETA_01_PZETA_01_p3+756.0*roa_rob_p2*MZETA_PZETA_p2*(MZETA_02_PZETA_02)+1296.0*roa_rob_p3*(MZETA_PZETA)*(MZETA_01_PZETA_01)*(MZETA_02_PZETA_02)-216.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_03_PZETA_03)))/(216.0*nKS*roa_rob_p25f6*MZETA_PZETA_p4);

    T_4000=(groa_grob*(43225.0*MZETA_PZETA_p4+41496.0*roa_rob*MZETA_PZETA_p3*(MZETA_10_PZETA_10)+39312.0*roa_rob_p2*MZETA_PZETA_p2*MZETA_10_PZETA_10_p2+36288.0*roa_rob_p3*(MZETA_PZETA)*MZETA_10_PZETA_10_p3+31104.0*roa_rob_p4*MZETA_10_PZETA_10_p4-19656.0*roa_rob_p2*MZETA_PZETA_p3*(MZETA_20_PZETA_20)-36288.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_10_PZETA_10)*(MZETA_20_PZETA_20)-46656.0*roa_rob_p4*(MZETA_PZETA)*MZETA_10_PZETA_10_p2*(MZETA_20_PZETA_20)+7776.0*roa_rob_p4*MZETA_PZETA_p2*MZETA_20_PZETA_20_p2+6048.0*roa_rob_p3*MZETA_PZETA_p3*(MZETA_30_PZETA_30)+10368.0*roa_rob_p4*MZETA_PZETA_p2*(MZETA_10_PZETA_10)*(MZETA_30_PZETA_30)-1296.0*roa_rob_p4*MZETA_PZETA_p3*(MZETA_40_PZETA_40)))/(1296.0*nKS*roa_rob_p31f6*MZETA_PZETA_p5);

    T_3100=(groa_grob*(43225.0*MZETA_PZETA_p4+10374.0*roa_rob*MZETA_PZETA_p3*(MZETA_01_PZETA_01)+31122.0*roa_rob*MZETA_PZETA_p3*(MZETA_10_PZETA_10)+19656.0*roa_rob_p2*MZETA_PZETA_p2*(MZETA_01_PZETA_01)*(MZETA_10_PZETA_10)+19656.0*roa_rob_p2*MZETA_PZETA_p2*MZETA_10_PZETA_10_p2+27216.0*roa_rob_p3*(MZETA_PZETA)*(MZETA_01_PZETA_01)*MZETA_10_PZETA_10_p2+9072.0*roa_rob_p3*(MZETA_PZETA)*MZETA_10_PZETA_10_p3+31104.0*roa_rob_p4*(MZETA_01_PZETA_01)*MZETA_10_PZETA_10_p3-9828.0*roa_rob_p2*MZETA_PZETA_p3*(MZETA_11_PZETA_11)-18144.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_10_PZETA_10)*(MZETA_11_PZETA_11)-23328.0*roa_rob_p4*(MZETA_PZETA)*MZETA_10_PZETA_10_p2*(MZETA_11_PZETA_11)-9828.0*roa_rob_p2*MZETA_PZETA_p3*(MZETA_20_PZETA_20)-9072.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_01_PZETA_01)*(MZETA_20_PZETA_20)-9072.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_10_PZETA_10)*(MZETA_20_PZETA_20)-23328.0*roa_rob_p4*(MZETA_PZETA)*(MZETA_01_PZETA_01)*(MZETA_10_PZETA_10)*(MZETA_20_PZETA_20)+7776.0*roa_rob_p4*MZETA_PZETA_p2*(MZETA_11_PZETA_11)*(MZETA_20_PZETA_20)+4536.0*roa_rob_p3*MZETA_PZETA_p3*(MZETA_21_PZETA_21)+7776.0*roa_rob_p4*MZETA_PZETA_p2*(MZETA_10_PZETA_10)*(MZETA_21_PZETA_21)+1512.0*roa_rob_p3*MZETA_PZETA_p3*(MZETA_30_PZETA_30)+2592.0*roa_rob_p4*MZETA_PZETA_p2*(MZETA_01_PZETA_01)*(MZETA_30_PZETA_30)-1296.0*roa_rob_p4*MZETA_PZETA_p3*(MZETA_31_PZETA_31)))/(1296.0*nKS*roa_rob_p31f6*MZETA_PZETA_p5);

    T_3010=(-1729.0*MZETA_PZETA_p3-1638.0*roa_rob*MZETA_PZETA_p2*(MZETA_10_PZETA_10)-1512.0*roa_rob_p2*(MZETA_PZETA)*MZETA_10_PZETA_10_p2-1296.0*roa_rob_p3*MZETA_10_PZETA_10_p3+756.0*roa_rob_p2*MZETA_PZETA_p2*(MZETA_20_PZETA_20)+1296.0*roa_rob_p3*(MZETA_PZETA)*(MZETA_10_PZETA_10)*(MZETA_20_PZETA_20)-216.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_30_PZETA_30))/(216.0*nKS*roa_rob_p25f6*MZETA_PZETA_p4);

    T_3001=T_3010;

    T_2200=(groa_grob*(43225.0*MZETA_PZETA_p4+20748.0*roa_rob*MZETA_PZETA_p3*(MZETA_01_PZETA_01)+6552.0*roa_rob_p2*MZETA_PZETA_p2*MZETA_01_PZETA_01_p2-3276.0*roa_rob_p2*MZETA_PZETA_p3*(MZETA_02_PZETA_02)+20748.0*roa_rob*MZETA_PZETA_p3*(MZETA_10_PZETA_10)+26208.0*roa_rob_p2*MZETA_PZETA_p2*(MZETA_01_PZETA_01)*(MZETA_10_PZETA_10)+18144.0*roa_rob_p3*(MZETA_PZETA)*MZETA_01_PZETA_01_p2*(MZETA_10_PZETA_10)-6048.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_02_PZETA_02)*(MZETA_10_PZETA_10)+6552.0*roa_rob_p2*MZETA_PZETA_p2*MZETA_10_PZETA_10_p2+18144.0*roa_rob_p3*(MZETA_PZETA)*(MZETA_01_PZETA_01)*MZETA_10_PZETA_10_p2+31104.0*roa_rob_p4*MZETA_01_PZETA_01_p2*MZETA_10_PZETA_10_p2-7776.0*roa_rob_p4*(MZETA_PZETA)*(MZETA_02_PZETA_02)*MZETA_10_PZETA_10_p2-13104.0*roa_rob_p2*MZETA_PZETA_p3*(MZETA_11_PZETA_11)-12096.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_01_PZETA_01)*(MZETA_11_PZETA_11)-12096.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_10_PZETA_10)*(MZETA_11_PZETA_11)-31104.0*roa_rob_p4*(MZETA_PZETA)*(MZETA_01_PZETA_01)*(MZETA_10_PZETA_10)*(MZETA_11_PZETA_11)+5184.0*roa_rob_p4*MZETA_PZETA_p2*MZETA_11_PZETA_11_p2+3024.0*roa_rob_p3*MZETA_PZETA_p3*(MZETA_12_PZETA_12)+5184.0*roa_rob_p4*MZETA_PZETA_p2*(MZETA_10_PZETA_10)*(MZETA_12_PZETA_12)-3276.0*roa_rob_p2*MZETA_PZETA_p3*(MZETA_20_PZETA_20)-6048.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_01_PZETA_01)*(MZETA_20_PZETA_20)-7776.0*roa_rob_p4*(MZETA_PZETA)*MZETA_01_PZETA_01_p2*(MZETA_20_PZETA_20)+2592.0*roa_rob_p4*MZETA_PZETA_p2*(MZETA_02_PZETA_02)*(MZETA_20_PZETA_20)+3024.0*roa_rob_p3*MZETA_PZETA_p3*(MZETA_21_PZETA_21)+5184.0*roa_rob_p4*MZETA_PZETA_p2*(MZETA_01_PZETA_01)*(MZETA_21_PZETA_21)-1296.0*roa_rob_p4*MZETA_PZETA_p3*(MZETA_22_PZETA_22)))/(1296.0*nKS*roa_rob_p31f6*MZETA_PZETA_p5);

    T_2110=(-1729.0*MZETA_PZETA_p3-546.0*roa_rob*MZETA_PZETA_p2*(MZETA_01_PZETA_01)-1092.0*roa_rob*MZETA_PZETA_p2*(MZETA_10_PZETA_10)-1008.0*roa_rob_p2*(MZETA_PZETA)*(MZETA_01_PZETA_01)*(MZETA_10_PZETA_10)-504.0*roa_rob_p2*(MZETA_PZETA)*MZETA_10_PZETA_10_p2-1296.0*roa_rob_p3*(MZETA_01_PZETA_01)*MZETA_10_PZETA_10_p2+504.0*roa_rob_p2*MZETA_PZETA_p2*(MZETA_11_PZETA_11)+864.0*roa_rob_p3*(MZETA_PZETA)*(MZETA_10_PZETA_10)*(MZETA_11_PZETA_11)+252.0*roa_rob_p2*MZETA_PZETA_p2*(MZETA_20_PZETA_20)+432.0*roa_rob_p3*(MZETA_PZETA)*(MZETA_01_PZETA_01)*(MZETA_20_PZETA_20)-216.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_21_PZETA_21))/(216.0*nKS*roa_rob_p25f6*MZETA_PZETA_p4);

    T_2101=T_2110;

    T_1300=(groa_grob*(43225.0*MZETA_PZETA_p4+31122.0*roa_rob*MZETA_PZETA_p3*(MZETA_01_PZETA_01)+19656.0*roa_rob_p2*MZETA_PZETA_p2*MZETA_01_PZETA_01_p2+9072.0*roa_rob_p3*(MZETA_PZETA)*MZETA_01_PZETA_01_p3-9828.0*roa_rob_p2*MZETA_PZETA_p3*(MZETA_02_PZETA_02)-9072.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_01_PZETA_01)*(MZETA_02_PZETA_02)+1512.0*roa_rob_p3*MZETA_PZETA_p3*(MZETA_03_PZETA_03)+10374.0*roa_rob*MZETA_PZETA_p3*(MZETA_10_PZETA_10)+19656.0*roa_rob_p2*MZETA_PZETA_p2*(MZETA_01_PZETA_01)*(MZETA_10_PZETA_10)+27216.0*roa_rob_p3*(MZETA_PZETA)*MZETA_01_PZETA_01_p2*(MZETA_10_PZETA_10)+31104.0*roa_rob_p4*MZETA_01_PZETA_01_p3*(MZETA_10_PZETA_10)-9072.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_02_PZETA_02)*(MZETA_10_PZETA_10)-23328.0*roa_rob_p4*(MZETA_PZETA)*(MZETA_01_PZETA_01)*(MZETA_02_PZETA_02)*(MZETA_10_PZETA_10)+2592.0*roa_rob_p4*MZETA_PZETA_p2*(MZETA_03_PZETA_03)*(MZETA_10_PZETA_10)-9828.0*roa_rob_p2*MZETA_PZETA_p3*(MZETA_11_PZETA_11)-18144.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_01_PZETA_01)*(MZETA_11_PZETA_11)-23328.0*roa_rob_p4*(MZETA_PZETA)*MZETA_01_PZETA_01_p2*(MZETA_11_PZETA_11)+7776.0*roa_rob_p4*MZETA_PZETA_p2*(MZETA_02_PZETA_02)*(MZETA_11_PZETA_11)+4536.0*roa_rob_p3*MZETA_PZETA_p3*(MZETA_12_PZETA_12)+7776.0*roa_rob_p4*MZETA_PZETA_p2*(MZETA_01_PZETA_01)*(MZETA_12_PZETA_12)-1296.0*roa_rob_p4*MZETA_PZETA_p3*(MZETA_13_PZETA_13)))/(1296.0*nKS*roa_rob_p31f6*MZETA_PZETA_p5);

    T_1210=(-1729.0*MZETA_PZETA_p3-1092.0*roa_rob*MZETA_PZETA_p2*(MZETA_01_PZETA_01)-504.0*roa_rob_p2*(MZETA_PZETA)*MZETA_01_PZETA_01_p2+252.0*roa_rob_p2*MZETA_PZETA_p2*(MZETA_02_PZETA_02)-546.0*roa_rob*MZETA_PZETA_p2*(MZETA_10_PZETA_10)-1008.0*roa_rob_p2*(MZETA_PZETA)*(MZETA_01_PZETA_01)*(MZETA_10_PZETA_10)-1296.0*roa_rob_p3*MZETA_01_PZETA_01_p2*(MZETA_10_PZETA_10)+432.0*roa_rob_p3*(MZETA_PZETA)*(MZETA_02_PZETA_02)*(MZETA_10_PZETA_10)+504.0*roa_rob_p2*MZETA_PZETA_p2*(MZETA_11_PZETA_11)+864.0*roa_rob_p3*(MZETA_PZETA)*(MZETA_01_PZETA_01)*(MZETA_11_PZETA_11)-216.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_12_PZETA_12))/(216.0*nKS*roa_rob_p25f6*MZETA_PZETA_p4);

    T_1201=T_1210;

    T_0400=(groa_grob*(43225.0*MZETA_PZETA_p4+41496.0*roa_rob*MZETA_PZETA_p3*(MZETA_01_PZETA_01)+39312.0*roa_rob_p2*MZETA_PZETA_p2*MZETA_01_PZETA_01_p2+36288.0*roa_rob_p3*(MZETA_PZETA)*MZETA_01_PZETA_01_p3+31104.0*roa_rob_p4*MZETA_01_PZETA_01_p4-19656.0*roa_rob_p2*MZETA_PZETA_p3*(MZETA_02_PZETA_02)-36288.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_01_PZETA_01)*(MZETA_02_PZETA_02)-46656.0*roa_rob_p4*(MZETA_PZETA)*MZETA_01_PZETA_01_p2*(MZETA_02_PZETA_02)+7776.0*roa_rob_p4*MZETA_PZETA_p2*MZETA_02_PZETA_02_p2+6048.0*roa_rob_p3*MZETA_PZETA_p3*(MZETA_03_PZETA_03)+10368.0*roa_rob_p4*MZETA_PZETA_p2*(MZETA_01_PZETA_01)*(MZETA_03_PZETA_03)-1296.0*roa_rob_p4*MZETA_PZETA_p3*(MZETA_04_PZETA_04)))/(1296.0*nKS*roa_rob_p31f6*MZETA_PZETA_p5);

    T_0310=(-1729.0*MZETA_PZETA_p3-1638.0*roa_rob*MZETA_PZETA_p2*(MZETA_01_PZETA_01)-1512.0*roa_rob_p2*(MZETA_PZETA)*MZETA_01_PZETA_01_p2-1296.0*roa_rob_p3*MZETA_01_PZETA_01_p3+756.0*roa_rob_p2*MZETA_PZETA_p2*(MZETA_02_PZETA_02)+1296.0*roa_rob_p3*(MZETA_PZETA)*(MZETA_01_PZETA_01)*(MZETA_02_PZETA_02)-216.0*roa_rob_p3*MZETA_PZETA_p2*(MZETA_03_PZETA_03))/(216.0*nKS*roa_rob_p25f6*MZETA_PZETA_p4);

    T_0301=T_0310;

/* All terms ok up to here! */
/* Powers of T derivatives */
    T_1000_p2 = T_1000*T_1000;
    T_1000_p3 = T_1000_p2*T_1000;
    T_1000_p4 = T_1000_p3*T_1000;

    T_0100_p2 = T_0100*T_0100;
    T_0100_p3 = T_0100_p2*T_0100;
    T_0100_p4 = T_0100_p3*T_0100;

    T_0010_p2 = T_0010*T_0010;
    T_0010_p3 = T_0010_p2*T_0010;
    T_0010_p4 = T_0010_p3*T_0010;

    T_0001_p2 = T_0001*T_0001;
    T_0001_p3 = T_0001_p2*T_0001;
    T_0001_p4 = T_0001_p3*T_0001;

    T_2000_p2 = T_2000*T_2000;

    T_0200_p2 = T_0200*T_0200;

    T_1100_p2 = T_1100*T_1100;

    T_1010_p2 = T_1010*T_1010;

    T_1001_p2 = T_1001*T_1001;

    T_0110_p2 = T_0110*T_0110;

    T_0101_p2 = T_0101*T_0101;

/* Derivatives of A */

    A_10=(delta*AEXP*(H0*EC_10-EC*H0_10))/(ADENOM_p2*H0_p2);

    A_01=(delta*AEXP*(H0*EC_01-EC*H0_01))/(ADENOM_p2*H0_p2);

    A_20=(delta*AEXP*(-((1.0+AEXP)*EC_p2*H0_10_p2)+2.0*EC*H0*H0_10*((1.0+AEXP)*EC_10+ADENOM*H0_10)+ADENOM*H0_p3*EC_20-H0_p2*((1.0+AEXP)*EC_10_p2+2.0*ADENOM*EC_10*H0_10+ADENOM*EC*H0_20)))/(ADENOM_p3*H0_p4);

    A_02=(delta*AEXP*(-((1.0+AEXP)*EC_p2*H0_01_p2)+2.0*EC*H0*H0_01*((1.0+AEXP)*EC_01+ADENOM*H0_01)+ADENOM*H0_p3*EC_02-H0_p2*((1.0+AEXP)*EC_01_p2+2.0*ADENOM*EC_01*H0_01+ADENOM*EC*H0_02)))/(ADENOM_p3*H0_p4);

    A_11=(delta*AEXP*(-((1.0+AEXP)*EC_p2*H0_01*H0_10)+EC*H0*((1.0+AEXP)*EC_01*H0_10+H0_01*((1.0+AEXP)*EC_10+2.0*ADENOM*H0_10))+ADENOM*H0_p3*EC_11-H0_p2*(EC_01*((1.0+AEXP)*EC_10+ADENOM*H0_10)+ADENOM*(H0_01*EC_10+EC*H0_11))))/(ADENOM_p3*H0_p4);

    A_30=(delta*AEXP*(-((1.0+4.0*AEXP+AEXP_p2)*EC_p3*H0_10_p3)+3.0*EC_p2*H0*H0_10_p2*((1.0+4.0*AEXP+AEXP_p2)*EC_10+2.0*(-1.0+AEXP_p2)*H0_10)-3.0*EC*H0_p2*H0_10*((1.0+4.0*AEXP+AEXP_p2)*EC_10_p2+4.0*(-1.0+AEXP_p2)*EC_10*H0_10+(ADENOM)*(2.0*(ADENOM)*H0_10_p2+(1.0+AEXP)*EC*H0_20))+H0_p3*((1.0+4.0*AEXP+AEXP_p2)*EC_10_p3+6.0*(-1.0+AEXP_p2)*EC_10_p2*H0_10+3.0*(ADENOM)*EC*H0_10*((1.0+AEXP)*EC_20+2.0*(ADENOM)*H0_20)+3.0*(ADENOM)*EC_10*(2.0*(ADENOM)*H0_10_p2+(1.0+AEXP)*EC*H0_20))+ADENOM_p2*H0_p5*EC_30-(ADENOM)*H0_p4*(3.0*EC_10*((1.0+AEXP)*EC_20+(ADENOM)*H0_20)+(ADENOM)*(3.0*H0_10*EC_20+EC*H0_30))))/(ADENOM_p4*H0_p6);

    A_21=(delta*AEXP*(-((1.0+4.0*AEXP+AEXP_p2)*EC_p3*H0_01*H0_10_p2)+EC_p2*H0*H0_10*((1.0+4.0*AEXP+AEXP_p2)*EC_01*H0_10+2.0*H0_01*((1.0+4.0*AEXP+AEXP_p2)*EC_10+3.0*(-1.0+AEXP_p2)*H0_10))-EC*H0_p2*(2.0*H0_10*(EC_01*((1.0+4.0*AEXP+AEXP_p2)*EC_10+2.0*(-1.0+AEXP_p2)*H0_10)+(-1.0+AEXP_p2)*EC*H0_11)+H0_01*((1.0+4.0*AEXP+AEXP_p2)*EC_10_p2+8.0*(-1.0+AEXP_p2)*EC_10*H0_10+(ADENOM)*(6.0*(ADENOM)*H0_10_p2+(1.0+AEXP)*EC*H0_20)))+H0_p3*(EC_01*((1.0+4.0*AEXP+AEXP_p2)*EC_10_p2+4.0*(-1.0+AEXP_p2)*EC_10*H0_10+(ADENOM)*(2.0*(ADENOM)*H0_10_p2+(1.0+AEXP)*EC*H0_20))+(ADENOM)*(2.0*EC*((1.0+AEXP)*EC_10*H0_11+H0_10*((1.0+AEXP)*EC_11+2.0*(ADENOM)*H0_11))+H0_01*(2.0*(1.0+AEXP)*EC_10_p2+4.0*(ADENOM)*EC_10*H0_10+EC*((1.0+AEXP)*EC_20+2.0*(ADENOM)*H0_20))))+ADENOM_p2*H0_p5*EC_21-(ADENOM)*H0_p4*(2.0*(ADENOM)*H0_10*EC_11+2.0*EC_10*((1.0+AEXP)*EC_11+(ADENOM)*H0_11)+EC_01*EC_20+AEXP*EC_01*EC_20-H0_01*EC_20+AEXP*H0_01*EC_20-EC_01*H0_20+AEXP*EC_01*H0_20-EC*H0_21+AEXP*EC*H0_21)))/(ADENOM_p4*H0_p6);

    A_12=(delta*AEXP*(-((1.0+4.0*AEXP+AEXP_p2)*EC_p3*H0_01_p2*H0_10)+EC_p2*H0*H0_01*(2.0*(1.0+4.0*AEXP+AEXP_p2)*EC_01*H0_10+H0_01*((1.0+4.0*AEXP+AEXP_p2)*EC_10+6.0*(-1.0+AEXP_p2)*H0_10))-EC*H0_p2*((1.0+4.0*AEXP+AEXP_p2)*EC_01_p2*H0_10+2.0*EC_01*H0_01*((1.0+4.0*AEXP+AEXP_p2)*EC_10+4.0*(-1.0+AEXP_p2)*H0_10)+ADENOM*((1.0+AEXP)*EC*H0_02*H0_10+H0_01_p2*(4.0*(1.0+AEXP)*EC_10+6.0*ADENOM*H0_10)+2.0*(1.0+AEXP)*EC*H0_01*H0_11))+H0_p3*(EC_01_p2*((1.0+4.0*AEXP+AEXP_p2)*EC_10+2.0*(-1.0+AEXP_p2)*H0_10)+2.0*ADENOM*EC_01*(2.0*H0_01*((1.0+AEXP)*EC_10+ADENOM*H0_10)+(1.0+AEXP)*EC*H0_11)+ADENOM*(2.0*ADENOM*H0_01_p2*EC_10+EC*((1.0+AEXP)*EC_02*H0_10+H0_02*((1.0+AEXP)*EC_10+2.0*ADENOM*H0_10))+2.0*EC*H0_01*((1.0+AEXP)*EC_11+2.0*ADENOM*H0_11)))+ADENOM_p2*H0_p5*EC_12-ADENOM*H0_p4*(ADENOM*H0_02*EC_10+EC_02*((1.0+AEXP)*EC_10+ADENOM*H0_10)+2.0*EC_01*EC_11+2.0*AEXP*EC_01*EC_11-2.0*H0_01*EC_11+2.0*AEXP*H0_01*EC_11-2.0*EC_01*H0_11+2.0*AEXP*EC_01*H0_11-EC*H0_12+AEXP*EC*H0_12)))/(ADENOM_p4*H0_p6);


    A_03=(delta*AEXP*(-((1.0+4.0*AEXP+AEXP_p2)*EC_p3*H0_01_p3)+3.0*EC_p2*H0*H0_01_p2*((1.0+4.0*AEXP+AEXP_p2)*EC_01+2.0*(-1.0+AEXP_p2)*H0_01)-3.0*EC*H0_p2*H0_01*((1.0+4.0*AEXP+AEXP_p2)*EC_01_p2+4.0*(-1.0+AEXP_p2)*EC_01*H0_01+(ADENOM)*(2.0*(ADENOM)*H0_01_p2+(1.0+AEXP)*EC*H0_02))+H0_p3*((1.0+4.0*AEXP+AEXP_p2)*EC_01_p3+6.0*(-1.0+AEXP_p2)*EC_01_p2*H0_01+3.0*(ADENOM)*EC*H0_01*((1.0+AEXP)*EC_02+2.0*(ADENOM)*H0_02)+3.0*(ADENOM)*EC_01*(2.0*(ADENOM)*H0_01_p2+(1.0+AEXP)*EC*H0_02))+ADENOM_p2*H0_p5*EC_03-(ADENOM)*H0_p4*(3.0*EC_01*((1.0+AEXP)*EC_02+(ADENOM)*H0_02)+(ADENOM)*(3.0*H0_01*EC_02+EC*H0_03))))/(ADENOM_p4*H0_p6);

    A_40=(delta*AEXP*(-((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_p4*H0_10_p4)+4.0*EC_p3*H0*H0_10_p3*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_10+3.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*H0_10)-6.0*EC_p2*H0_p2*H0_10_p2*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_10_p2+6.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*EC_10*H0_10+ADENOM*(6.0*(-1.0+AEXP_p2)*H0_10_p2+(1.0+4.0*AEXP+AEXP_p2)*EC*H0_20))+2.0*EC*H0_p3*H0_10*(2.0*(1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_10_p3+18.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*EC_10_p2*H0_10+6.0*ADENOM*EC_10*(6.0*(-1.0+AEXP_p2)*H0_10_p2+(1.0+4.0*AEXP+AEXP_p2)*EC*H0_20)+3.0*ADENOM*H0_10*(4.0*ADENOM_p2*H0_10_p2+EC*((1.0+4.0*AEXP+AEXP_p2)*EC_20+6.0*(-1.0+AEXP_p2)*H0_20)))-H0_p4*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_10_p4+12.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*EC_10_p3*H0_10+6.0*ADENOM*EC_10_p2*(6.0*(-1.0+AEXP_p2)*H0_10_p2+(1.0+4.0*AEXP+AEXP_p2)*EC*H0_20)+12.0*ADENOM*EC_10*H0_10*(2.0*ADENOM_p2*H0_10_p2+EC*((1.0+4.0*AEXP+AEXP_p2)*EC_20+4.0*(-1.0+AEXP_p2)*H0_20))+ADENOM_p2*EC*(3.0*(1.0+AEXP)*EC*H0_20_p2+12.0*H0_10_p2*(2.0*(1.0+AEXP)*EC_20+3.0*ADENOM*H0_20)+4.0*(1.0+AEXP)*EC*H0_10*H0_30))+2.0*ADENOM*H0_p5*(3.0*EC_10_p2*((1.0+4.0*AEXP+AEXP_p2)*EC_20+2.0*(-1.0+AEXP_p2)*H0_20)+2.0*ADENOM*EC_10*(6.0*H0_10*((1.0+AEXP)*EC_20+ADENOM*H0_20)+(1.0+AEXP)*EC*H0_30)+ADENOM*(6.0*ADENOM*H0_10_p2*EC_20+3.0*EC*H0_20*((1.0+AEXP)*EC_20+ADENOM*H0_20)+2.0*EC*H0_10*((1.0+AEXP)*EC_30+2.0*ADENOM*H0_30)))+ADENOM_p3*H0_p7*EC_40-ADENOM_p2*H0_p6*(3.0*(1.0+AEXP)*EC_20_p2+6.0*ADENOM*EC_20*H0_20+4.0*EC_10*((1.0+AEXP)*EC_30+ADENOM*H0_30)+ADENOM*(4.0*H0_10*EC_30+EC*H0_40))))/(ADENOM_p5*H0_p8);

    A_31=(delta*AEXP*(-((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_p4*H0_01*H0_10_p3)+EC_p3*H0*H0_10_p2*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_01*H0_10+3.0*H0_01*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_10+4.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*H0_10))-3.0*EC_p2*H0_p2*H0_10*(H0_10*(EC_01*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_10+3.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*H0_10)+(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*EC*H0_11)+H0_01*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_10_p2+9.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*EC_10*H0_10+ADENOM*(12.0*(-1.0+AEXP_p2)*H0_10_p2+(1.0+4.0*AEXP+AEXP_p2)*EC*H0_20)))+EC*H0_p3*(H0_01*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_10_p3+18.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*EC_10_p2*H0_10+3.0*ADENOM*EC_10*(18.0*(-1.0+AEXP_p2)*H0_10_p2+(1.0+4.0*AEXP+AEXP_p2)*EC*H0_20)+3.0*ADENOM*H0_10*(8.0*ADENOM_p2*H0_10_p2+EC*((1.0+4.0*AEXP+AEXP_p2)*EC_20+6.0*(-1.0+AEXP_p2)*H0_20)))+3.0*H0_10*(ADENOM*EC*(2.0*(1.0+4.0*AEXP+AEXP_p2)*EC_10*H0_11+H0_10*((1.0+4.0*AEXP+AEXP_p2)*EC_11+6.0*(-1.0+AEXP_p2)*H0_11))+EC_01*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_10_p2+6.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*EC_10*H0_10+ADENOM*(6.0*(-1.0+AEXP_p2)*H0_10_p2+(1.0+4.0*AEXP+AEXP_p2)*EC*H0_20))))+ADENOM*H0_p5*(3.0*EC_10_p2*((1.0+4.0*AEXP+AEXP_p2)*EC_11+2.0*(-1.0+AEXP_p2)*H0_11)+3.0*EC_10*(4.0*ADENOM*H0_10*((1.0+AEXP)*EC_11+ADENOM*H0_11)+EC_01*((1.0+4.0*AEXP+AEXP_p2)*EC_20+2.0*(-1.0+AEXP_p2)*H0_20)+ADENOM*(2.0*H0_01*((1.0+AEXP)*EC_20+ADENOM*H0_20)+(1.0+AEXP)*EC*H0_21))+ADENOM*(6.0*ADENOM*H0_10_p2*EC_11+3.0*H0_10*(2.0*ADENOM*H0_01*EC_20+2.0*EC_01*((1.0+AEXP)*EC_20+ADENOM*H0_20)+EC*((1.0+AEXP)*EC_21+2.0*ADENOM*H0_21))+EC*(3.0*(1.0+AEXP)*EC_11*H0_20+3.0*H0_11*((1.0+AEXP)*EC_20+2.0*ADENOM*H0_20)+H0_01*EC_30+AEXP*H0_01*EC_30+EC_01*H0_30+AEXP*EC_01*H0_30-2.0*H0_01*H0_30+2.0*AEXP*H0_01*H0_30)))-H0_p4*(EC_01*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_10_p3+9.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*EC_10_p2*H0_10+3.0*ADENOM*EC_10*(6.0*(-1.0+AEXP_p2)*H0_10_p2+(1.0+4.0*AEXP+AEXP_p2)*EC*H0_20)+3.0*ADENOM*H0_10*(2.0*ADENOM_p2*H0_10_p2+EC*((1.0+4.0*AEXP+AEXP_p2)*EC_20+4.0*(-1.0+AEXP_p2)*H0_20)))+ADENOM*(3.0*EC*((1.0+4.0*AEXP+AEXP_p2)*EC_10_p2*H0_11+2.0*EC_10*H0_10*((1.0+4.0*AEXP+AEXP_p2)*EC_11+4.0*(-1.0+AEXP_p2)*H0_11)+ADENOM*(H0_10_p2*(4.0*(1.0+AEXP)*EC_11+6.0*ADENOM*H0_11)+(1.0+AEXP)*EC*H0_11*H0_20+(1.0+AEXP)*EC*H0_10*H0_21))+H0_01*(3.0*(1.0+4.0*AEXP+AEXP_p2)*EC_10_p3+18.0*(-1.0+AEXP_p2)*EC_10_p2*H0_10+3.0*EC_10*(6.0*ADENOM_p2*H0_10_p2+EC*((1.0+4.0*AEXP+AEXP_p2)*EC_20+4.0*(-1.0+AEXP_p2)*H0_20))+ADENOM*EC*(6.0*H0_10*(2.0*(1.0+AEXP)*EC_20+3.0*ADENOM*H0_20)+(1.0+AEXP)*EC*H0_30))))+ADENOM_p3*H0_p7*EC_31-ADENOM_p2*H0_p6*(3.0*ADENOM*H0_11*EC_20+3.0*EC_11*((1.0+AEXP)*EC_20+ADENOM*H0_20)+3.0*EC_10*EC_21+3.0*AEXP*EC_10*EC_21-3.0*H0_10*EC_21+3.0*AEXP*H0_10*EC_21-3.0*EC_10*H0_21+3.0*AEXP*EC_10*H0_21+EC_01*EC_30+AEXP*EC_01*EC_30-H0_01*EC_30+AEXP*H0_01*EC_30-EC_01*H0_30+AEXP*EC_01*H0_30-EC*H0_31+AEXP*EC*H0_31)))/(ADENOM_p5*H0_p8);

    A_22 = (delta*AEXP*(-((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_p4*H0_01_p2*H0_10_p2)+2.0*EC_p3*H0*H0_01*H0_10*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_01*H0_10+H0_01*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_10+6.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*H0_10))-EC_p2*H0_p2*(((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_01_p2+(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*EC*H0_02)*H0_10_p2+2.0*H0_01*H0_10*(EC_01*(2.0*(1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_10+9.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*H0_10)+2.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*EC*H0_11)+H0_01_p2*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_10_p2+18.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*EC_10*H0_10+(ADENOM)*(36.0*(-1.0+AEXP_p2)*H0_10_p2+(1.0+4.0*AEXP+AEXP_p2)*EC*H0_20)))+EC*H0_p3*(2.0*EC_01_p2*H0_10*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_10+3.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*H0_10)+(ADENOM)*(EC*H0_10*((1.0+4.0*AEXP+AEXP_p2)*EC_02*H0_10+2.0*H0_02*((1.0+4.0*AEXP+AEXP_p2)*EC_10+3.0*(-1.0+AEXP_p2)*H0_10))+4.0*EC*H0_01*((1.0+4.0*AEXP+AEXP_p2)*EC_10*H0_11+H0_10*((1.0+4.0*AEXP+AEXP_p2)*EC_11+6.0*(-1.0+AEXP_p2)*H0_11))+H0_01_p2*(6.0*(1.0+4.0*AEXP+AEXP_p2)*EC_10_p2+36.0*(-1.0+AEXP_p2)*EC_10*H0_10+24.0*ADENOM_p2*H0_10_p2+EC*((1.0+4.0*AEXP+AEXP_p2)*EC_20+6.0*(-1.0+AEXP_p2)*H0_20)))+2.0*EC_01*(2.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*EC*H0_10*H0_11+H0_01*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_10_p2+12.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*EC_10*H0_10+(ADENOM)*(18.0*(-1.0+AEXP_p2)*H0_10_p2+(1.0+4.0*AEXP+AEXP_p2)*EC*H0_20))))+(ADENOM)*H0_p5*(4.0*EC_01*EC_10*EC_11+16.0*AEXP*EC_01*EC_10*EC_11+4.0*AEXP_p2*EC_01*EC_10*EC_11-8.0*H0_01*EC_10*EC_11+8.0*AEXP_p2*H0_01*EC_10*EC_11-8.0*EC_01*H0_10*EC_11+8.0*AEXP_p2*EC_01*H0_10*EC_11+8.0*H0_01*H0_10*EC_11-16.0*AEXP*H0_01*H0_10*EC_11+8.0*AEXP_p2*H0_01*H0_10*EC_11-8.0*EC_01*EC_10*H0_11+8.0*AEXP_p2*EC_01*EC_10*H0_11+8.0*H0_01*EC_10*H0_11-16.0*AEXP*H0_01*EC_10*H0_11+8.0*AEXP_p2*H0_01*EC_10*H0_11+8.0*EC_01*H0_10*H0_11-16.0*AEXP*EC_01*H0_10*H0_11+8.0*AEXP_p2*EC_01*H0_10*H0_11-4.0*EC*EC_11*H0_11+4.0*AEXP_p2*EC*EC_11*H0_11+4.0*EC*H0_11_p2-8.0*AEXP*EC*H0_11_p2+4.0*AEXP_p2*EC*H0_11_p2-2.0*EC*H0_10*EC_12+2.0*AEXP_p2*EC*H0_10*EC_12-2.0*EC*EC_10*H0_12+2.0*AEXP_p2*EC*EC_10*H0_12+4.0*EC*H0_10*H0_12-8.0*AEXP*EC*H0_10*H0_12+4.0*AEXP_p2*EC*H0_10*H0_12+EC_01_p2*EC_20+4.0*AEXP*EC_01_p2*EC_20+AEXP_p2*EC_01_p2*EC_20-4.0*EC_01*H0_01*EC_20+4.0*AEXP_p2*EC_01*H0_01*EC_20+2.0*H0_01_p2*EC_20-4.0*AEXP*H0_01_p2*EC_20+2.0*AEXP_p2*H0_01_p2*EC_20-2.0*EC_01_p2*H0_20+2.0*AEXP_p2*EC_01_p2*H0_20+4.0*EC_01*H0_01*H0_20-8.0*AEXP*EC_01*H0_01*H0_20+4.0*AEXP_p2*EC_01*H0_01*H0_20+(ADENOM)*H0_02*(2.0*(1.0+AEXP)*EC_10_p2+4.0*(ADENOM)*EC_10*H0_10+EC*((1.0+AEXP)*EC_20+2.0*(ADENOM)*H0_20))+EC_02*((1.0+4.0*AEXP+AEXP_p2)*EC_10_p2+4.0*(-1.0+AEXP_p2)*EC_10*H0_10+(ADENOM)*(2.0*(ADENOM)*H0_10_p2+(1.0+AEXP)*EC*H0_20))-2.0*EC*H0_01*EC_21+2.0*AEXP_p2*EC*H0_01*EC_21-2.0*EC*EC_01*H0_21+2.0*AEXP_p2*EC*EC_01*H0_21+4.0*EC*H0_01*H0_21-8.0*AEXP*EC*H0_01*H0_21+4.0*AEXP_p2*EC*H0_01*H0_21)-H0_p4*(EC_01_p2*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_10_p2+6.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*EC_10*H0_10+(ADENOM)*(6.0*(-1.0+AEXP_p2)*H0_10_p2+(1.0+4.0*AEXP+AEXP_p2)*EC*H0_20))+2.0*(ADENOM)*EC_01*(2.0*EC*((1.0+4.0*AEXP+AEXP_p2)*EC_10*H0_11+H0_10*((1.0+4.0*AEXP+AEXP_p2)*EC_11+4.0*(-1.0+AEXP_p2)*H0_11))+H0_01*(3.0*(1.0+4.0*AEXP+AEXP_p2)*EC_10_p2+12.0*(-1.0+AEXP_p2)*EC_10*H0_10+6.0*ADENOM_p2*H0_10_p2+EC*((1.0+4.0*AEXP+AEXP_p2)*EC_20+4.0*(-1.0+AEXP_p2)*H0_20)))+(ADENOM)*(2.0*(ADENOM)*H0_01_p2*(3.0*(1.0+AEXP)*EC_10_p2+6.0*(ADENOM)*EC_10*H0_10+EC*(2.0*(1.0+AEXP)*EC_20+3.0*(ADENOM)*H0_20))+EC*(2.0*(EC_02*H0_10*((1.0+4.0*AEXP+AEXP_p2)*EC_10+2.0*(-1.0+AEXP_p2)*H0_10)+(-1.0+AEXP_p2)*EC*(H0_11_p2+H0_10*H0_12))+H0_02*((1.0+4.0*AEXP+AEXP_p2)*EC_10_p2+8.0*(-1.0+AEXP_p2)*EC_10*H0_10+(ADENOM)*(6.0*(ADENOM)*H0_10_p2+(1.0+AEXP)*EC*H0_20)))+2.0*EC*H0_01*(2.0*EC_10*((1.0+4.0*AEXP+AEXP_p2)*EC_11+4.0*(-1.0+AEXP_p2)*H0_11)+(ADENOM)*(4.0*H0_10*(2.0*(1.0+AEXP)*EC_11+3.0*(ADENOM)*H0_11)+(1.0+AEXP)*EC*H0_21))))+ADENOM_p3*H0_p7*EC_22-ADENOM_p2*H0_p6*(2.0*(1.0+AEXP)*EC_11_p2+4.0*(ADENOM)*EC_11*H0_11-2.0*H0_10*EC_12+2.0*AEXP*H0_10*EC_12+2.0*EC_10*((1.0+AEXP)*EC_12+(ADENOM)*H0_12)+EC_02*EC_20+AEXP*EC_02*EC_20-H0_02*EC_20+AEXP*H0_02*EC_20-EC_02*H0_20+AEXP*EC_02*H0_20+2.0*EC_01*EC_21+2.0*AEXP*EC_01*EC_21-2.0*H0_01*EC_21+2.0*AEXP*H0_01*EC_21-2.0*EC_01*H0_21+2.0*AEXP*EC_01*H0_21-EC*H0_22+AEXP*EC*H0_22)))/(ADENOM_p5*H0_p8);

    A_13=(delta*AEXP*(-((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_p4*H0_01_p3*H0_10)+EC_p3*H0*H0_01_p2*(3.0*(1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_01*H0_10+H0_01*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_10+12.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*H0_10))-3.0*EC_p2*H0_p2*H0_01*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_01_p2*H0_10+EC_01*H0_01*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_10+9.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*H0_10)+ADENOM*((1.0+4.0*AEXP+AEXP_p2)*EC*H0_02*H0_10+3.0*H0_01_p2*((1.0+4.0*AEXP+AEXP_p2)*EC_10+4.0*(-1.0+AEXP_p2)*H0_10)+(1.0+4.0*AEXP+AEXP_p2)*EC*H0_01*H0_11))+EC*H0_p3*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_01_p3*H0_10+3.0*EC_01_p2*H0_01*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_10+6.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*H0_10)+3.0*ADENOM*EC_01*((1.0+4.0*AEXP+AEXP_p2)*EC*H0_02*H0_10+6.0*H0_01_p2*((1.0+4.0*AEXP+AEXP_p2)*EC_10+3.0*(-1.0+AEXP_p2)*H0_10)+2.0*(1.0+4.0*AEXP+AEXP_p2)*EC*H0_01*H0_11)+3.0*ADENOM*H0_01*(2.0*ADENOM*H0_01_p2*(3.0*(1.0+AEXP)*EC_10+4.0*ADENOM*H0_10)+EC*((1.0+4.0*AEXP+AEXP_p2)*EC_02*H0_10+H0_02*((1.0+4.0*AEXP+AEXP_p2)*EC_10+6.0*(-1.0+AEXP_p2)*H0_10))+EC*H0_01*((1.0+4.0*AEXP+AEXP_p2)*EC_11+6.0*(-1.0+AEXP_p2)*H0_11)))+ADENOM*H0_p5*(3.0*EC_01_p2*((1.0+4.0*AEXP+AEXP_p2)*EC_11+2.0*(-1.0+AEXP_p2)*H0_11)+3.0*EC_01*(EC_02*((1.0+4.0*AEXP+AEXP_p2)*EC_10+2.0*(-1.0+AEXP_p2)*H0_10)+ADENOM*(2.0*H0_02*((1.0+AEXP)*EC_10+ADENOM*H0_10)+4.0*H0_01*((1.0+AEXP)*EC_11+ADENOM*H0_11)+(1.0+AEXP)*EC*H0_12))+ADENOM*(6.0*ADENOM*H0_01_p2*EC_11+EC*((1.0+AEXP)*EC_03*H0_10+H0_03*((1.0+AEXP)*EC_10+2.0*ADENOM*H0_10)+3.0*((1.0+AEXP)*EC_02*H0_11+H0_02*((1.0+AEXP)*EC_11+2.0*ADENOM*H0_11)))+3.0*H0_01*(2.0*ADENOM*H0_02*EC_10+2.0*EC_02*((1.0+AEXP)*EC_10+ADENOM*H0_10)+EC*((1.0+AEXP)*EC_12+2.0*ADENOM*H0_12))))-H0_p4*(EC_01_p3*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_10+3.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*H0_10)+3.0*ADENOM*EC_01_p2*(3.0*H0_01*((1.0+4.0*AEXP+AEXP_p2)*EC_10+2.0*(-1.0+AEXP_p2)*H0_10)+(1.0+4.0*AEXP+AEXP_p2)*EC*H0_11)+3.0*ADENOM*EC_01*(6.0*ADENOM*H0_01_p2*((1.0+AEXP)*EC_10+ADENOM*H0_10)+EC*((1.0+4.0*AEXP+AEXP_p2)*EC_02*H0_10+H0_02*((1.0+4.0*AEXP+AEXP_p2)*EC_10+4.0*(-1.0+AEXP_p2)*H0_10))+2.0*EC*H0_01*((1.0+4.0*AEXP+AEXP_p2)*EC_11+4.0*(-1.0+AEXP_p2)*H0_11))+ADENOM*(6.0*ADENOM_p2*H0_01_p3*EC_10+6.0*ADENOM*EC*H0_01_p2*(2.0*(1.0+AEXP)*EC_11+3.0*ADENOM*H0_11)+(-1.0+AEXP_p2)*EC_p2*(H0_03*H0_10+3.0*H0_02*H0_11)+3.0*EC*H0_01*(EC_02*((1.0+4.0*AEXP+AEXP_p2)*EC_10+4.0*(-1.0+AEXP_p2)*H0_10)+ADENOM*(H0_02*(4.0*(1.0+AEXP)*EC_10+6.0*ADENOM*H0_10)+(1.0+AEXP)*EC*H0_12))))+ADENOM_p3*H0_p7*EC_13-ADENOM_p2*H0_p6*(ADENOM*H0_03*EC_10+EC_03*((1.0+AEXP)*EC_10+ADENOM*H0_10)+3.0*EC_02*EC_11+3.0*AEXP*EC_02*EC_11-3.0*H0_02*EC_11+3.0*AEXP*H0_02*EC_11-3.0*EC_02*H0_11+3.0*AEXP*EC_02*H0_11+3.0*EC_01*EC_12+3.0*AEXP*EC_01*EC_12-3.0*H0_01*EC_12+3.0*AEXP*H0_01*EC_12-3.0*EC_01*H0_12+3.0*AEXP*EC_01*H0_12-EC*H0_13+AEXP*EC*H0_13)))/(ADENOM_p5*H0_p8);


    A_04=(delta*AEXP*(-((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_p4*H0_01_p4)+4.0*EC_p3*H0*H0_01_p3*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_01+3.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*H0_01)-6.0*EC_p2*H0_p2*H0_01_p2*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_01_p2+6.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*EC_01*H0_01+ADENOM*(6.0*(-1.0+AEXP_p2)*H0_01_p2+(1.0+4.0*AEXP+AEXP_p2)*EC*H0_02))+2.0*EC*H0_p3*H0_01*(2.0*(1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_01_p3+18.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*EC_01_p2*H0_01+6.0*ADENOM*EC_01*(6.0*(-1.0+AEXP_p2)*H0_01_p2+(1.0+4.0*AEXP+AEXP_p2)*EC*H0_02)+3.0*ADENOM*H0_01*(4.0*ADENOM_p2*H0_01_p2+EC*((1.0+4.0*AEXP+AEXP_p2)*EC_02+6.0*(-1.0+AEXP_p2)*H0_02)))-H0_p4*((1.0+11.0*AEXP+11.0*AEXP_p2+AEXP_p3)*EC_01_p4+12.0*(-1.0-3.0*AEXP+3.0*AEXP_p2+AEXP_p3)*EC_01_p3*H0_01+6.0*ADENOM*EC_01_p2*(6.0*(-1.0+AEXP_p2)*H0_01_p2+(1.0+4.0*AEXP+AEXP_p2)*EC*H0_02)+12.0*ADENOM*EC_01*H0_01*(2.0*ADENOM_p2*H0_01_p2+EC*((1.0+4.0*AEXP+AEXP_p2)*EC_02+4.0*(-1.0+AEXP_p2)*H0_02))+ADENOM_p2*EC*(3.0*(1.0+AEXP)*EC*H0_02_p2+12.0*H0_01_p2*(2.0*(1.0+AEXP)*EC_02+3.0*ADENOM*H0_02)+4.0*(1.0+AEXP)*EC*H0_01*H0_03))+2.0*ADENOM*H0_p5*(3.0*EC_01_p2*((1.0+4.0*AEXP+AEXP_p2)*EC_02+2.0*(-1.0+AEXP_p2)*H0_02)+2.0*ADENOM*EC_01*(6.0*H0_01*((1.0+AEXP)*EC_02+ADENOM*H0_02)+(1.0+AEXP)*EC*H0_03)+ADENOM*(6.0*ADENOM*H0_01_p2*EC_02+3.0*EC*H0_02*((1.0+AEXP)*EC_02+ADENOM*H0_02)+2.0*EC*H0_01*((1.0+AEXP)*EC_03+2.0*ADENOM*H0_03)))+ADENOM_p3*H0_p7*EC_04-ADENOM_p2*H0_p6*(3.0*(1.0+AEXP)*EC_02_p2+6.0*ADENOM*EC_02*H0_02+4.0*EC_01*((1.0+AEXP)*EC_03+ADENOM*H0_03)+ADENOM*(4.0*H0_01*EC_03+EC*H0_04))))/(ADENOM_p5*H0_p8);


/* Powers of A derivatives */
    A_10_p2 = A_10*A_10;
    A_10_p3 = A_10_p2*A_10;
    A_10_p4 = A_10_p3*A_10;

    A_01_p2 = A_01*A_01;
    A_01_p3 = A_01_p2*A_01;
    A_01_p4 = A_01_p3*A_01;

    A_20_p2 = A_20*A_20;

    A_02_p2 =  A_02*A_02;

/* Derivatives of H with respect to H0, T, A */
/* Derivatives with respect to H0 higher then first order vanish */
/* Terms used in the derivatives */
    At2 = A*T_p2;
    At2_p2=At2*At2;

    At_nom = 1.0 + At2;
    At_denom = At_nom + At2_p2;
    delta_A = delta + A;
    At_term = (1.0+delta_A*T_p2*At_nom);

/* Powers of terms */
    At_denom_p2 = At_denom*At_denom;
    At_denom_p3 = At_denom_p2*At_denom;
    At_denom_p4 = At_denom_p3*At_denom;

    At_term_p2 = At_term*At_term;
    At_term_p3 = At_term_p2*At_term;
    At_term_p4 = At_term_p3*At_term;

    delta_A_p2 = delta_A*delta_A;
    delta_A_p3 = delta_A_p2*delta_A;

/* First order derivatives */
    H_A =-((delta*A*H0*T_p6*(2.0+At2))/(At_denom*At_term));

    H_T =(2.0*delta*H0*T*(1.0+2.0*At2))/(At_denom*At_term);

    H_H = log(1.0+(delta-delta/At_denom)/A);

/* Second order derivatives */
    H_AA =(delta*H0*T_p6*(-2.0-2.0*delta_A*T_p2-2.0*(delta-2.0*A)*A*T_p4+2.0*A_p2*(delta+4.0*A)*T_p6+4.0*A_p3*(delta+2.0*A)*T_p8+A_p4*(delta+2.0*A)*T_p10))/(At_denom_p2*At_term_p2);

    H_TT =(2.0*delta*H0*(1.0+T_p2*(-delta+A*(4.0+T_p2*(-4.0*delta+A*(-5.0-2.0*(7.0*delta+8.0*A)*T_p2-19.0*A*delta_A*T_p4-10.0*A_p2*delta_A*T_p6))))))/(At_denom_p2*At_term_p2);

    H_HT =(2.0*delta*T*(1.0+2.0*At2))/(At_denom*At_term);

    H_TA =(-2.0*delta*A*H0*T_p5*(6.0+4.0*(delta+3.0*A)*T_p2+A*(7.0*delta+12.0*A)*T_p4+2.0*A_p2*(2.0*delta+3.0*A)*T_p6))/(At_denom_p2*At_term_p2);

    H_HA =-((delta*A*T_p6*(2.0+At2))/(At_denom*At_term));

/* Third order derivatives */
    H_TTT =(4.0*delta*H0*T*(-3.0*delta+T_p2*(delta_p2+A*(6.0*delta*(-3.0+delta*T_p2)+A*(-30.0-3.0*(23.0*delta+32.0*A)*T_p2+6.0*(2.0*delta-21.0*A)*delta_A*T_p4+A*(29.0*delta_p2-39.0*delta*A-66.0*A_p2)*T_p6+6.0*A_p2*delta_A*(13.0*delta+11.0*A)*T_p8+126.0*A_p3*delta_A_p2*T_p10+96.0*A_p4*delta_A_p2*T_p12+30.0*A_p5*delta_A_p2*T_p14)))))/(At_denom_p3*At_term_p3);

    H_AAA =(2.0*delta*H0*T_p8*(3.0+6.0*(delta+3.0*A)*T_p2+3.0*(delta_p2+10.0*delta*A+11.0*A_p2)*T_p4+3.0*A*(4.0*delta_p2+17.0*delta*A+10.0*A_p2)*T_p6+3.0*delta*A_p2*(6.0*delta+13.0*A)*T_p8+(7.0*delta_p2*A_p3-30.0*A_p5)*T_p10-3.0*A_p4*(2.0*delta_p2+8.0*delta*A+11.0*A_p2)*T_p12-6.0*A_p5*(delta_p2+3.0*A*delta_A)*T_p14-A_p6*(delta_p2+3.0*A*delta_A)*T_p16))/(At_denom_p3*At_term_p3);

    H_HTT =(2.0*delta*(1.0+T_p2*(-delta+A*(4.0+T_p2*(-4.0*delta+A*(-5.0-2.0*(7.0*delta+8.0*A)*T_p2-19.0*A*delta_A*T_p4-10.0*A_p2*delta_A*T_p6))))))/(At_denom_p2*At_term_p2);

    H_HAA =(delta*T_p6*(-2.0-2.0*delta_A*T_p2-2.0*(delta-2.0*A)*A*T_p4+2.0*A_p2*(delta+4.0*A)*T_p6+4.0*A_p3*(delta+2.0*A)*T_p8+A_p4*(delta+2.0*A)*T_p10))/(At_denom_p2*At_term_p2);

    H_TTA =(2.0*delta*A*H0*T_p4*(-30.0-2.0*(17.0*delta+48.0*A)*T_p2-3.0*(delta+2.0*A)*(4.0*delta+21.0*A)*T_p4-3.0*A*(9.0*delta_p2+18.0*delta*A+22.0*A_p2)*T_p6+A_p2*(-2.0*delta_p2+83.0*delta*A+66.0*A_p2)*T_p8+6.0*A_p3*(9.0*delta_p2+31.0*delta*A+21.0*A_p2)*T_p10+3.0*A_p4*delta_A*(19.0*delta+32.0*A)*T_p12+10.0*A_p5*delta_A*(2.0*delta+3.0*A)*T_p14))/(At_denom_p3*At_term_p3);

    H_TAA =(4.0*delta*H0*T_p5*(-3.0+T_p2*(-(delta*(5.0+2.0*delta*T_p2))+A*(-6.0+9.0*(-delta+A)*T_p2-3.0*(delta_p2-2.0*A*(2.0*delta+7.0*A))*T_p4+A*(6.0*delta_p2+52.0*delta*A+75.0*A_p2)*T_p6+4.0*A_p2*(5.0*delta_p2+18.0*A*delta_A)*T_p8+3.0*A_p3*(6.0*delta_p2+17.0*delta*A+14.0*A_p2)*T_p10+2.0*A_p4*(3.0*delta_p2+8.0*delta*A+6.0*A_p2)*T_p12))))/(At_denom_p3*At_term_p3);

    H_HTA =(-2.0*delta*A*T_p5*(6.0+4.0*(delta+3.0*A)*T_p2+A*(7.0*delta+12.0*A)*T_p4+2.0*A_p2*(2.0*delta+3.0*A)*T_p6))/(At_denom_p2*At_term_p2);

/* Fourth order derivatives */
    H_TTTT =(12.0*delta*H0*(-delta+2.0*(delta-3.0*A)*(3.0*delta+5.0*A)*T_p2-(delta_p3+4.0*A*(-12.0*delta_p2+4.0*delta*A+25.0*A_p2))*T_p4+4.0*A*(-2.0*delta_p3+45.0*delta_p2*A+47.0*delta*A_p2+10.0*A_p3)*T_p6+A_p2*(-34.0*delta_p3+530.0*delta_p2*A+1199.0*delta*A_p2+690.0*A_p3)*T_p8+2.0*A_p3*(-35.0*delta_p3+645.0*delta_p2*A+1570.0*delta*A_p2+891.0*A_p3)*T_p10+2.0*A_p4*(-27.0*delta_p3+4.0*A*(273.0*delta_p2+8.0*A*(74.0*delta+39.0*A)))*T_p12+8.0*A_p6*delta_A*(265.0*delta+264.0*A)*T_p14+A_p6*delta_A*(-137.0*delta_p2+767.0*delta*A+894.0*A_p2)*T_p16-6.0*A_p7*delta_A_p2*(76.0*delta+33.0*A)*T_p18-2.0*A_p8*delta_A_p2*(277.0*delta+252.0*A)*T_p20-308.0*A_p9*delta_A_p3*T_p22-70.0*A_p10*delta_A_p3*T_p24))/(At_denom_p4*At_term_p4);

    H_AAAA =(6.0*delta*H0*T_p12*(-2.0*(delta+10.0*A)-4.0*(delta_p2+14.0*delta*A+25.0*A_p2)*T_p2-2.0*(delta_p3+26.0*delta_p2*A+121.0*delta*A_p2+120.0*A_p3)*T_p4-4.0*A*(4.0*delta_p3+49.0*delta_p2*A+131.0*delta*A_p2+90.0*A_p3)*T_p6-2.0*A_p2*(27.0*delta_p3+2.0*A*(93.0*delta_p2+173.0*delta*A+84.0*A_p2))*T_p8-8.0*A_p3*(11.0*delta_p3+50.0*delta_p2*A+70.0*delta*A_p2+21.0*A_p3)*T_p10+2.0*A_p4*(-35.0*delta_p3+4.0*A*(-29.0*delta_p2-28.0*delta*A+3.0*A_p2))*T_p12+4.0*A_p5*(-4.0*delta_p3-6.0*delta_p2*A+10.0*delta*A_p2+27.0*A_p3)*T_p14+4.0*A_p6*(3.0*delta_p3+14.0*delta_p2*A+25.0*delta*A_p2+21.0*A_p3)*T_p16+8.0*A_p7*(delta+2.0*A)*(delta_p2+2.0*A*delta_A)*T_p18+A_p8*(delta+2.0*A)*(delta_p2+2.0*A*delta_A)*T_p20))/(At_denom_p4*At_term_p4);

    H_HTTT =(4.0*delta*T*(-3.0*delta+T_p2*(delta_p2+A*(6.0*delta*(-3.0+delta*T_p2)+A*(-30.0-3.0*(23.0*delta+32.0*A)*T_p2+6.0*(2.0*delta-21.0*A)*delta_A*T_p4+A*(29.0*delta_p2-39.0*delta*A-66.0*A_p2)*T_p6+6.0*A_p2*delta_A*(13.0*delta+11.0*A)*T_p8+126.0*A_p3*delta_A_p2*T_p10+96.0*A_p4*delta_A_p2*T_p12+30.0*A_p5*delta_A_p2*T_p14)))))/(At_denom_p3*At_term_p3);

    H_TTTA =(24.0*delta*A*H0*T_p3*(-10.0-2.0*(6.0*delta+19.0*A)*T_p2-2.0*(4.0*delta_p2+9.0*delta*A+12.0*A_p2)*T_p4+(-2.0*delta_p3-7.0*delta_p2*A+122.0*delta*A_p2+117.0*A_p3)*T_p6+A*(-delta_p3+100.0*delta_p2*A+506.0*delta*A_p2+390.0*A_p3)*T_p8+2.0*A_p2*(14.0*delta_p3+175.0*delta_p2*A+466.0*delta*A_p2+300.0*A_p3)*T_p10+4.0*A_p3*(22.0*delta_p3+128.0*delta_p2*A+245.0*delta*A_p2+138.0*A_p3)*T_p12+A_p4*(99.0*delta_p3+355.0*delta_p2*A+536.0*delta*A_p2+279.0*A_p3)*T_p14+2.0*A_p5*delta_A*(12.0*delta_p2-8.0*delta*A+3.0*A_p2)*T_p16-2.0*A_p6*delta_A*(22.0*delta_p2+73.0*delta*A+46.0*A_p2)*T_p18-2.0*A_p7*delta_A_p2*(19.0*delta+32.0*A)*T_p20-5.0*A_p8*delta_A_p2*(2.0*delta+3.0*A)*T_p22))/(At_denom_p4*At_term_p4);

    H_TAAA =(12.0*delta*H0*T_p7*(4.0+T_p2*(11.0*delta+32.0*A+2.0*(5.0*delta_p2+38.0*delta*A+42.0*A_p2)*T_p2+(3.0*delta_p3+4.0*A*(15.0*delta_p2+46.0*delta*A+27.0*A_p2))*T_p4+4.0*A*(4.0*delta_p3+33.0*delta_p2*A+52.0*delta*A_p2+6.0*A_p3)*T_p6+A_p2*(32.0*delta_p3+(7.0*delta-8.0*A)*A*(16.0*delta+21.0*A))*T_p8-2.0*A_p3*(-8.0*delta_p3+A*(27.0*delta_p2+4.0*A*(37.0*delta+42.0*A)))*T_p10-4.0*A_p4*(8.0*delta_p3+57.0*delta_p2*A+119.0*delta*A_p2+90.0*A_p3)*T_p12-8.0*A_p5*(delta+2.0*A)*(7.0*delta_p2+17.0*delta*A+15.0*A_p2)*T_p14-A_p6*(34.0*delta_p3+A*(132.0*delta_p2+5.0*A*(37.0*delta+20.0*A)))*T_p16-2.0*A_p7*(4.0*delta_p3+5.0*A*(3.0*delta_p2+2.0*A*(2.0*delta+A)))*T_p18)))/(At_denom_p4*At_term_p4);

    H_HAAA =(2.0*delta*T_p8*(3.0+6.0*(delta+3.0*A)*T_p2+3.0*(delta_p2+10.0*delta*A+11.0*A_p2)*T_p4+3.0*A*(4.0*delta_p2+17.0*delta*A+10.0*A_p2)*T_p6+3.0*delta*A_p2*(6.0*delta+13.0*A)*T_p8+(7.0*delta_p2*A_p3-30.0*A_p5)*T_p10-3.0*A_p4*(2.0*delta_p2+8.0*delta*A+11.0*A_p2)*T_p12-6.0*A_p5*(delta_p2+3.0*A*delta_A)*T_p14-A_p6*(delta_p2+3.0*A*delta_A)*T_p16))/(At_denom_p3*At_term_p3);

    H_HATT =(2.0*delta*A*T_p4*(-30.0-2.0*(17.0*delta+48.0*A)*T_p2-3.0*(delta+2.0*A)*(4.0*delta+21.0*A)*T_p4-3.0*A*(9.0*delta_p2+18.0*delta*A+22.0*A_p2)*T_p6+A_p2*(-2.0*delta_p2+83.0*delta*A+66.0*A_p2)*T_p8+6.0*A_p3*(9.0*delta_p2+31.0*delta*A+21.0*A_p2)*T_p10+3.0*A_p4*delta_A*(19.0*delta+32.0*A)*T_p12+10.0*A_p5*delta_A*(2.0*delta+3.0*A)*T_p14))/(At_denom_p3*At_term_p3);

    H_HTAA =(4.0*delta*T_p5*(-3.0+T_p2*(-(delta*(5.0+2.0*delta*T_p2))+A*(-6.0+9.0*(-delta+A)*T_p2-3.0*(delta_p2-2.0*A*(2.0*delta+7.0*A))*T_p4+A*(6.0*delta_p2+52.0*delta*A+75.0*A_p2)*T_p6+4.0*A_p2*(5.0*delta_p2+18.0*A*delta_A)*T_p8+3.0*A_p3*(6.0*delta_p2+17.0*delta*A+14.0*A_p2)*T_p10+2.0*A_p4*(3.0*delta_p2+8.0*delta*A+6.0*A_p2)*T_p12))))/(At_denom_p3*At_term_p3);

    H_TTAA =(4.0*delta*H0*T_p4*(-15.0-4.0*(8.0*delta+9.0*A)*T_p2+(-23.0*delta_p2-55.0*delta*A+132.0*A_p2)*T_p4+(-6.0*delta_p3-22.0*delta_p2*A+318.0*delta*A_p2+684.0*A_p3)*T_p6+A*(-3.0*delta_p3+A*(290.0*delta_p2+9.0*A*(148.0*delta+167.0*A)))*T_p8+4.0*A_p2*(21.0*delta_p3+254.0*delta_p2*A+600.0*delta*A_p2+486.0*A_p3)*T_p10+A_p3*(264.0*delta_p3+A*(1475.0*delta_p2+3.0*A*(787.0*delta+504.0*A)))*T_p12+A_p4*(297.0*delta_p3+986.0*delta_p2*A+996.0*delta*A_p2+504.0*A_p3)*T_p14-A_p5*(-72.0*delta_p3+A*(46.0*delta_p2+63.0*A*(8.0*delta+5.0*A)))*T_p16-2.0*A_p6*(66.0*delta_p3+308.0*delta_p2*A+489.0*delta*A_p2+240.0*A_p3)*T_p18-A_p7*delta_A*(114.0*delta_p2+5.0*A*(65.0*delta+54.0*A))*T_p20-10.0*A_p8*delta_A*(3.0*delta_p2+8.0*delta*A+6.0*A_p2)*T_p22))/(At_denom_p4*At_term_p4);

/* Final expressions (general derivatives of F[H[roa,rob], T[roa,rob,groa,grob], A[roa,rob]]) */

/* First order derivatives */
    ds->df10000 += factor*(
        A_10*H_A+H0_10*H_H+H_T*T_1000 + EC_10
        );

    ds->df01000 += factor*(
        A_01*H_A+H0_01*H_H+H_T*T_0100 + EC_01
        );

    ds->df00100 += factor*(
        H_T*T_0010
        );

    ds->df00010 += factor*(
        H_T*T_0001
        );

/* Second order derivatives */
    ds->df20000 += factor*(
        A_20*H_A+H0_20*H_H+A_10*(A_10*H_AA+2.0*H0_10*H_HA)+2.0*(A_10*H_TA+H0_10*H_HT)*T_1000+H_TT*T_1000_p2+H_T*T_2000 + EC_20
        );

    ds->df11000 += factor*(
        A_11*H_A+H0_11*H_H+A_01*(A_10*H_AA+H0_10*H_HA+H_TA*T_1000)+T_0100*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+H0_01*(A_10*H_HA+H_HT*T_1000)+H_T*T_1100 + EC_11
        );

    ds->df10100 += factor*(
        T_0010*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+H_T*T_1010
        );

    ds->df10010 += factor*(
        T_0001*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+H_T*T_1001
        );

    ds->df02000 += factor*(
        A_02*H_A+H0_02*H_H+A_01*(A_01*H_AA+2.0*H0_01*H_HA)+2.0*(A_01*H_TA+H0_01*H_HT)*T_0100+H_TT*T_0100_p2+H_T*T_0200 + EC_02
        );

    ds->df01100 += factor*(
        T_0010*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)+H_T*T_0110
        );

    ds->df01010 += factor*(
        T_0001*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)+H_T*T_0101
        );

    ds->df00200 += factor*(
        H_TT*T_0010_p2
        );

    ds->df00110 += factor*(
        H_TT*T_0001*T_0010
        );

    ds->df00020 += factor*(
        H_TT*T_0001_p2
        );

/* Third order derivatives */
    ds->df30000 += factor*(
        A_30*H_A+A_10_p3*H_AAA+H0_30*H_H+3.0*A_10_p2*(H0_10*H_HAA+H_TAA*T_1000)+3.0*A_10*(A_20*H_AA+H0_20*H_HA+2.0*H0_10*H_HTA*T_1000+H_TTA*T_1000_p2+H_TA*T_2000)+T_1000*(3.0*A_20*H_TA+3.0*H0_20*H_HT+H_TTT*T_1000_p2+3.0*H_TT*T_2000)+3.0*H0_10*(A_20*H_HA+H_HTT*T_1000_p2+H_HT*T_2000)+H_T*T_3000 + EC_30
        );

    ds->df21000 += factor*(
        A_21*H_A+H0_21*H_H+2.0*H0_10*A_11*H_HA+H0_01*A_20*H_HA+A_20*H_TA*T_0100+H0_20*H_HT*T_0100+A_10_p2*(A_01*H_AAA+H0_01*H_HAA+H_TAA*T_0100)+2.0*A_11*H_TA*T_1000+2.0*H0_11*H_HT*T_1000+2.0*H0_10*H_HTT*T_0100*T_1000+H0_01*H_HTT*T_1000_p2+H_TTT*T_0100*T_1000_p2+2.0*H0_10*H_HT*T_1100+2.0*H_TT*T_1000*T_1100+2.0*A_10*(A_11*H_AA+H0_11*H_HA+H0_10*(A_01*H_HAA+H_HTA*T_0100)+(A_01*H_TAA+H0_01*H_HTA+H_TTA*T_0100)*T_1000+H_TA*T_1100)+(H0_01*H_HT+H_TT*T_0100)*T_2000+A_01*(A_20*H_AA+H0_20*H_HA+2.0*H0_10*H_HTA*T_1000+H_TTA*T_1000_p2+H_TA*T_2000)+H_T*T_2100 + EC_21
        );

    ds->df20100 += factor*(
        2.0*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)*T_1010+T_0010*(A_20*H_TA+H0_20*H_HT+A_10*(A_10*H_TAA+2.0*H0_10*H_HTA)+2.0*(A_10*H_TTA+H0_10*H_HTT)*T_1000+H_TTT*T_1000_p2+H_TT*T_2000)+H_T*T_2010
        );

    ds->df20010 += factor*(
        2.0*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)*T_1001+T_0001*(A_20*H_TA+H0_20*H_HT+A_10*(A_10*H_TAA+2.0*H0_10*H_HTA)+2.0*(A_10*H_TTA+H0_10*H_HTT)*T_1000+H_TTT*T_1000_p2+H_TT*T_2000)+H_T*T_2001
        );

    ds->df12000 += factor*(
        A_12*H_A+H0_12*H_H+A_11*(A_01*H_AA+2.0*H0_01*H_HA)+A_02*(A_10*H_AA+H0_10*H_HA+H_TA*T_1000)+T_0200*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+T_0100_p2*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+H0_02*(A_10*H_HA+H_HT*T_1000)+A_01*(A_11*H_AA+2.0*H0_11*H_HA+A_01*(A_10*H_AAA+H0_10*H_HAA+H_TAA*T_1000)+2.0*H0_01*(A_10*H_HAA+H_HTA*T_1000))+2.0*T_0100*(A_11*H_TA+H0_11*H_HT+A_01*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+H0_01*(A_10*H_HTA+H_HTT*T_1000))+2.0*(A_01*H_TA+H0_01*H_HT)*T_1100+2.0*H_TT*T_0100*T_1100+H_T*T_1200 + EC_12
        );

    ds->df11100 += factor*(
        T_0110*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_1010+T_0010*(A_11*H_TA+H0_11*H_HT+A_01*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+T_0100*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+H0_01*(A_10*H_HTA+H_HTT*T_1000)+H_TT*T_1100)+H_T*T_1110
        );

    ds->df11010 += factor*(
        T_0101*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_1001+T_0001*(A_11*H_TA+H0_11*H_HT+A_01*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+T_0100*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+H0_01*(A_10*H_HTA+H_HTT*T_1000)+H_TT*T_1100)+H_T*T_1101
        );

    ds->df10200 += factor*(
        T_0010*(T_0010*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+2.0*H_TT*T_1010)
        );

    ds->df10110 += factor*(
        T_0010*(T_0001*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+H_TT*T_1001)+H_TT*T_0001*T_1010
        );

    ds->df10020 += factor*(
        T_0001*(T_0001*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+2.0*H_TT*T_1001)
        );

    ds->df03000 += factor*(
        A_03*H_A+A_01_p3*H_AAA+H0_03*H_H+3.0*A_01_p2*(H0_01*H_HAA+H_TAA*T_0100)+3.0*A_01*(A_02*H_AA+H0_02*H_HA+2.0*H0_01*H_HTA*T_0100+H_TTA*T_0100_p2+H_TA*T_0200)+T_0100*(3.0*A_02*H_TA+3.0*H0_02*H_HT+H_TTT*T_0100_p2+3.0*H_TT*T_0200)+3.0*H0_01*(A_02*H_HA+H_HTT*T_0100_p2+H_HT*T_0200)+H_T*T_0300 + EC_03
        );

    ds->df02100 += factor*(
        2.0*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_0110+T_0010*(A_02*H_TA+H0_02*H_HT+A_01*(A_01*H_TAA+2.0*H0_01*H_HTA)+2.0*(A_01*H_TTA+H0_01*H_HTT)*T_0100+H_TTT*T_0100_p2+H_TT*T_0200)+H_T*T_0210
        );

    ds->df02010 += factor*(
        2.0*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_0101+T_0001*(A_02*H_TA+H0_02*H_HT+A_01*(A_01*H_TAA+2.0*H0_01*H_HTA)+2.0*(A_01*H_TTA+H0_01*H_HTT)*T_0100+H_TTT*T_0100_p2+H_TT*T_0200)+H_T*T_0201
        );

    ds->df01200 += factor*(
        T_0010*(T_0010*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)+2.0*H_TT*T_0110)
        );

    ds->df01110 += factor*(
        T_0010*(T_0001*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)+H_TT*T_0101)+H_TT*T_0001*T_0110
        );

    ds->df01020 += factor*(
        T_0001*(T_0001*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)+2.0*H_TT*T_0101)
        );

    ds->df00300 += factor*(
        H_TTT*T_0010_p3
        );

    ds->df00210 += factor*(
        H_TTT*T_0001*T_0010_p2
        );

    ds->df00120 += factor*(
        H_TTT*T_0001_p2*T_0010
        );

    ds->df00030 += factor*(
        H_TTT*T_0001_p3
        );

/* Fourth order derivatives */

    ds->df40000 += factor*(
        A_40*H_A+3.0*A_20_p2*H_AA+A_10_p4*H_AAAA+H0_40*H_H+4.0*H0_10*A_30*H_HA+4.0*A_30*H_TA*T_1000+4.0*H0_30*H_HT*T_1000+6.0*H0_20*H_HTT*T_1000_p2+4.0*H0_10*H_HTTT*T_1000_p3+H_TTTT*T_1000_p4+4.0*A_10_p3*(H0_10*H_HAAA+H_TAAA*T_1000)+6.0*H0_20*H_HT*T_2000+12.0*H0_10*H_HTT*T_1000*T_2000+6.0*H_TTT*T_1000_p2*T_2000+3.0*H_TT*T_2000_p2+6.0*A_20*(H0_20*H_HA+A_10*(A_10*H_AAA+2.0*H0_10*H_HAA)+2.0*(A_10*H_TAA+H0_10*H_HTA)*T_1000+H_TTA*T_1000_p2+H_TA*T_2000)+6.0*A_10_p2*(H0_20*H_HAA+2.0*H0_10*H_HTAA*T_1000+H_TTAA*T_1000_p2+H_TAA*T_2000)+4.0*(H0_10*H_HT+H_TT*T_1000)*T_3000+4.0*A_10*(A_30*H_AA+H0_30*H_HA+T_1000*(3.0*H0_20*H_HTA+3.0*H0_10*H_HATT*T_1000+H_TTTA*T_1000_p2)+3.0*(H0_10*H_HTA+H_TTA*T_1000)*T_2000+H_TA*T_3000)+H_T*T_4000 + EC_40
        );

    ds->df31000 += factor*(
        A_31*H_A+3.0*A_10_p2*A_11*H_AAA+H0_31*H_H+A_30*(A_01*H_AA+H0_01*H_HA+H_TA*T_0100)+A_10_p3*(A_01*H_AAAA+H0_01*H_HAAA+H_TAAA*T_0100)+H0_30*(A_01*H_HA+H_HT*T_0100)+6.0*A_10*A_11*(H0_10*H_HAA+H_TAA*T_1000)+3.0*A_10_p2*(H0_11*H_HAA+H0_10*(A_01*H_HAAA+H_HTAA*T_0100)+(A_01*H_TAAA+H0_01*H_HTAA+H_TTAA*T_0100)*T_1000+H_TAA*T_1100)+3.0*A_11*(A_20*H_AA+H0_20*H_HA+2.0*H0_10*H_HTA*T_1000+H_TTA*T_1000_p2+H_TA*T_2000)+T_1100*(3.0*A_20*H_TA+3.0*H0_20*H_HT+H_TTT*T_1000_p2+3.0*H_TT*T_2000)+3.0*H0_11*(A_20*H_HA+H_HTT*T_1000_p2+H_HT*T_2000)+3.0*A_10*(A_21*H_AA+H0_21*H_HA+A_20*(A_01*H_AAA+H0_01*H_HAA+H_TAA*T_0100)+H0_20*(A_01*H_HAA+H_HTA*T_0100)+2.0*H0_11*H_HTA*T_1000+2.0*H0_10*(A_01*H_HTAA+H_HATT*T_0100)*T_1000+(A_01*H_TTAA+H0_01*H_HATT+H_TTTA*T_0100)*T_1000_p2+2.0*H0_10*H_HTA*T_1100+2.0*H_TTA*T_1000*T_1100+(A_01*H_TAA+H0_01*H_HTA+H_TTA*T_0100)*T_2000+H_TA*T_2100)+T_1000*(3.0*A_21*H_TA+3.0*H0_21*H_HT+3.0*A_20*(A_01*H_TAA+H0_01*H_HTA+H_TTA*T_0100)+3.0*H0_20*(A_01*H_HTA+H_HTT*T_0100)+(A_01*H_TTTA+H0_01*H_HTTT+H_TTTT*T_0100)*T_1000_p2+2.0*H_TTT*T_1000*T_1100+3.0*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)*T_2000+3.0*H_TT*T_2100)+3.0*H0_10*(A_21*H_HA+A_20*(A_01*H_HAA+H_HTA*T_0100)+(A_01*H_HATT+H_HTTT*T_0100)*T_1000_p2+2.0*H_HTT*T_1000*T_1100+(A_01*H_HTA+H_HTT*T_0100)*T_2000+H_HT*T_2100)+(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_3000+H_T*T_3100 + EC_31
        );

    ds->df30100 += factor*(
        A_30*H_TA*T_0010+A_10_p3*H_TAAA*T_0010+H0_30*H_HT*T_0010+3.0*A_10_p2*(T_0010*(H0_10*H_HTAA+H_TTAA*T_1000)+H_TAA*T_1010)+T_1010*(3.0*A_20*H_TA+3.0*H0_20*H_HT+H_TTT*T_1000_p2+3.0*H_TT*T_2000)+3.0*A_10*(A_20*H_TAA*T_0010+H0_20*H_HTA*T_0010+2.0*H0_10*H_HATT*T_0010*T_1000+H_TTTA*T_0010*T_1000_p2+2.0*H0_10*H_HTA*T_1010+2.0*H_TTA*T_1000*T_1010+H_TTA*T_0010*T_2000+H_TA*T_2010)+T_1000*(3.0*A_20*H_TTA*T_0010+3.0*H0_20*H_HTT*T_0010+H_TTTT*T_0010*T_1000_p2+2.0*H_TTT*T_1000*T_1010+3.0*H_TTT*T_0010*T_2000+3.0*H_TT*T_2010)+3.0*H0_10*(A_20*H_HTA*T_0010+H_HTTT*T_0010*T_1000_p2+H_HTT*(2.0*T_1000*T_1010+T_0010*T_2000)+H_HT*T_2010)+H_TT*T_0010*T_3000+H_T*T_3010
        );

    ds->df30010 += factor*(
        A_30*H_TA*T_0001+A_10_p3*H_TAAA*T_0001+H0_30*H_HT*T_0001+3.0*A_10_p2*(T_0001*(H0_10*H_HTAA+H_TTAA*T_1000)+H_TAA*T_1001)+T_1001*(3.0*A_20*H_TA+3.0*H0_20*H_HT+H_TTT*T_1000_p2+3.0*H_TT*T_2000)+3.0*A_10*(A_20*H_TAA*T_0001+H0_20*H_HTA*T_0001+2.0*H0_10*H_HATT*T_0001*T_1000+H_TTTA*T_0001*T_1000_p2+2.0*H0_10*H_HTA*T_1001+2.0*H_TTA*T_1000*T_1001+H_TTA*T_0001*T_2000+H_TA*T_2001)+T_1000*(3.0*A_20*H_TTA*T_0001+3.0*H0_20*H_HTT*T_0001+H_TTTT*T_0001*T_1000_p2+2.0*H_TTT*T_1000*T_1001+3.0*H_TTT*T_0001*T_2000+3.0*H_TT*T_2001)+3.0*H0_10*(A_20*H_HTA*T_0001+H_HTTT*T_0001*T_1000_p2+H_HTT*(2.0*T_1000*T_1001+T_0001*T_2000)+H_HT*T_2001)+H_TT*T_0001*T_3000+H_T*T_3001
        );

    ds->df22000 += factor*(
        A_22*H_A+H0_22*H_H+A_21*(A_01*H_AA+2.0*H0_01*H_HA)+2.0*A_12*(A_10*H_AA+H0_10*H_HA+H_TA*T_1000)+2.0*H0_12*(A_10*H_HA+H_HT*T_1000)+2.0*A_11*(A_11*H_AA+2.0*H0_11*H_HA+A_01*(A_10*H_AAA+H0_10*H_HAA+H_TAA*T_1000)+2.0*H0_01*(A_10*H_HAA+H_HTA*T_1000))+4.0*T_0100*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)*T_1100+4.0*(A_11*H_TA+H0_11*H_HT+A_01*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+H0_01*(A_10*H_HTA+H_HTT*T_1000))*T_1100+2.0*H_TT*T_1100_p2+2.0*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)*T_1200+A_02*(A_20*H_AA+H0_20*H_HA+A_10*(A_10*H_AAA+2.0*H0_10*H_HAA)+2.0*(A_10*H_TAA+H0_10*H_HTA)*T_1000+H_TTA*T_1000_p2+H_TA*T_2000)+T_0200*(A_20*H_TA+H0_20*H_HT+A_10*(A_10*H_TAA+2.0*H0_10*H_HTA)+2.0*(A_10*H_TTA+H0_10*H_HTT)*T_1000+H_TTT*T_1000_p2+H_TT*T_2000)+T_0100_p2*(A_20*H_TTA+H0_20*H_HTT+A_10*(A_10*H_TTAA+2.0*H0_10*H_HATT)+2.0*(A_10*H_TTTA+H0_10*H_HTTT)*T_1000+H_TTTT*T_1000_p2+H_TTT*T_2000)+H0_02*(A_20*H_HA+A_10_p2*H_HAA+2.0*A_10*H_HTA*T_1000+H_HTT*T_1000_p2+H_HT*T_2000)+2.0*T_0100*(A_21*H_TA+H0_21*H_HT+2.0*H0_10*A_11*H_HTA+A_10_p2*(A_01*H_TAAA+H0_01*H_HTAA)+2.0*(A_11*H_TTA+H0_11*H_HTT)*T_1000+2.0*A_10*(A_11*H_TAA+H0_11*H_HTA+A_01*H0_10*H_HTAA+(A_01*H_TTAA+H0_01*H_HATT)*T_1000)+A_01*(A_20*H_TAA+H0_20*H_HTA+2.0*H0_10*H_HATT*T_1000+H_TTTA*T_1000_p2+H_TTA*T_2000)+H0_01*(A_20*H_HTA+H_HTTT*T_1000_p2+H_HTT*T_2000))+A_01*(A_21*H_AA+A_10_p2*(A_01*H_AAAA+2.0*H0_01*H_HAAA)+2.0*A_10*(A_11*H_AAA+2.0*H0_11*H_HAA+A_01*H0_10*H_HAAA+(A_01*H_TAAA+2.0*H0_01*H_HTAA)*T_1000)+A_01*(A_20*H_AAA+H0_20*H_HAA+2.0*H0_10*H_HTAA*T_1000+H_TTAA*T_1000_p2+H_TAA*T_2000)+2.0*(H0_21*H_HA+H0_10*A_11*H_HAA+(A_11*H_TAA+2.0*H0_11*H_HTA)*T_1000+H0_01*(A_20*H_HAA+H_HATT*T_1000_p2+H_HTA*T_2000)))+2.0*(A_01*H_TA+H0_01*H_HT)*T_2100+2.0*H_TT*T_0100*T_2100+H_T*T_2200 + EC_22
        );

    ds->df21100  += factor*(
        A_21*H_TA*T_0010+H0_21*H_HT*T_0010+2.0*H0_10*A_11*H_HTA*T_0010+H0_01*A_20*H_HTA*T_0010+A_20*H_TTA*T_0010*T_0100+H0_20*H_HTT*T_0010*T_0100+A_20*H_TA*T_0110+H0_20*H_HT*T_0110+A_10_p2*(T_0010*(A_01*H_TAAA+H0_01*H_HTAA+H_TTAA*T_0100)+H_TAA*T_0110)+2.0*A_11*H_TTA*T_0010*T_1000+2.0*H0_11*H_HTT*T_0010*T_1000+2.0*H0_10*H_HTTT*T_0010*T_0100*T_1000+2.0*H0_10*H_HTT*T_0110*T_1000+H0_01*H_HTTT*T_0010*T_1000_p2+H_TTTT*T_0010*T_0100*T_1000_p2+H_TTT*T_0110*T_1000_p2+2.0*A_11*H_TA*T_1010+2.0*H0_11*H_HT*T_1010+2.0*H0_10*H_HTT*T_0100*T_1010+2.0*H0_01*H_HTT*T_1000*T_1010+2.0*H_TTT*T_0100*T_1000*T_1010+2.0*H0_10*H_HTT*T_0010*T_1100+2.0*H_TTT*T_0010*T_1000*T_1100+2.0*H_TT*T_1010*T_1100+2.0*H0_10*H_HT*T_1110+2.0*H_TT*T_1000*T_1110+2.0*A_10*(A_11*H_TAA*T_0010+H0_11*H_HTA*T_0010+H0_10*(T_0010*(A_01*H_HTAA+H_HATT*T_0100)+H_HTA*T_0110)+(T_0010*(A_01*H_TTAA+H0_01*H_HATT+H_TTTA*T_0100)+H_TTA*T_0110)*T_1000+(A_01*H_TAA+H0_01*H_HTA+H_TTA*T_0100)*T_1010+H_TTA*T_0010*T_1100+H_TA*T_1110)+(T_0010*(H0_01*H_HTT+H_TTT*T_0100)+H_TT*T_0110)*T_2000+(H0_01*H_HT+H_TT*T_0100)*T_2010+A_01*(A_20*H_TAA*T_0010+H0_20*H_HTA*T_0010+2.0*H0_10*H_HATT*T_0010*T_1000+H_TTTA*T_0010*T_1000_p2+2.0*H0_10*H_HTA*T_1010+2.0*H_TTA*T_1000*T_1010+H_TTA*T_0010*T_2000+H_TA*T_2010)+H_TT*T_0010*T_2100+H_T*T_2110
        );

    ds->df21010 += factor*(
        A_21*H_TA*T_0001+H0_21*H_HT*T_0001+2.0*H0_10*A_11*H_HTA*T_0001+H0_01*A_20*H_HTA*T_0001+A_20*H_TTA*T_0001*T_0100+H0_20*H_HTT*T_0001*T_0100+A_20*H_TA*T_0101+H0_20*H_HT*T_0101+A_10_p2*(T_0001*(A_01*H_TAAA+H0_01*H_HTAA+H_TTAA*T_0100)+H_TAA*T_0101)+2.0*A_11*H_TTA*T_0001*T_1000+2.0*H0_11*H_HTT*T_0001*T_1000+2.0*H0_10*H_HTTT*T_0001*T_0100*T_1000+2.0*H0_10*H_HTT*T_0101*T_1000+H0_01*H_HTTT*T_0001*T_1000_p2+H_TTTT*T_0001*T_0100*T_1000_p2+H_TTT*T_0101*T_1000_p2+2.0*A_11*H_TA*T_1001+2.0*H0_11*H_HT*T_1001+2.0*H0_10*H_HTT*T_0100*T_1001+2.0*H0_01*H_HTT*T_1000*T_1001+2.0*H_TTT*T_0100*T_1000*T_1001+2.0*H0_10*H_HTT*T_0001*T_1100+2.0*H_TTT*T_0001*T_1000*T_1100+2.0*H_TT*T_1001*T_1100+2.0*H0_10*H_HT*T_1101+2.0*H_TT*T_1000*T_1101+2.0*A_10*(A_11*H_TAA*T_0001+H0_11*H_HTA*T_0001+H0_10*(T_0001*(A_01*H_HTAA+H_HATT*T_0100)+H_HTA*T_0101)+(T_0001*(A_01*H_TTAA+H0_01*H_HATT+H_TTTA*T_0100)+H_TTA*T_0101)*T_1000+(A_01*H_TAA+H0_01*H_HTA+H_TTA*T_0100)*T_1001+H_TTA*T_0001*T_1100+H_TA*T_1101)+(T_0001*(H0_01*H_HTT+H_TTT*T_0100)+H_TT*T_0101)*T_2000+(H0_01*H_HT+H_TT*T_0100)*T_2001+A_01*(A_20*H_TAA*T_0001+H0_20*H_HTA*T_0001+2.0*H0_10*H_HATT*T_0001*T_1000+H_TTTA*T_0001*T_1000_p2+2.0*H0_10*H_HTA*T_1001+2.0*H_TTA*T_1000*T_1001+H_TTA*T_0001*T_2000+H_TA*T_2001)+H_TT*T_0001*T_2100+H_T*T_2101
        );

    ds->df20200 += factor*(
        A_20*H_TTA*T_0010_p2+A_10_p2*H_TTAA*T_0010_p2+H0_20*H_HTT*T_0010_p2+2.0*H0_10*H_HTTT*T_0010_p2*T_1000+H_TTTT*T_0010_p2*T_1000_p2+4.0*H0_10*H_HTT*T_0010*T_1010+4.0*H_TTT*T_0010*T_1000*T_1010+2.0*H_TT*T_1010_p2+2.0*A_10*T_0010*(T_0010*(H0_10*H_HATT+H_TTTA*T_1000)+2.0*H_TTA*T_1010)+H_TTT*T_0010_p2*T_2000+2.0*H_TT*T_0010*T_2010
        );

    ds->df20110 += factor*(
        2.0*(T_0001*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+H_TT*T_1001)*T_1010+T_0010*(A_20*H_TTA*T_0001+H0_20*H_HTT*T_0001+A_10*(A_10*H_TTAA+2.0*H0_10*H_HATT)*T_0001+2.0*(A_10*H_TTTA+H0_10*H_HTTT)*T_0001*T_1000+H_TTTT*T_0001*T_1000_p2+2.0*(A_10*H_TTA+H0_10*H_HTT)*T_1001+2.0*H_TTT*T_1000*T_1001+H_TTT*T_0001*T_2000+H_TT*T_2001)+H_TT*T_0001*T_2010
        );

    ds->df20020 += factor*(
        A_20*H_TTA*T_0001_p2+A_10_p2*H_TTAA*T_0001_p2+H0_20*H_HTT*T_0001_p2+2.0*H0_10*H_HTTT*T_0001_p2*T_1000+H_TTTT*T_0001_p2*T_1000_p2+4.0*H0_10*H_HTT*T_0001*T_1001+4.0*H_TTT*T_0001*T_1000*T_1001+2.0*H_TT*T_1001_p2+2.0*A_10*T_0001*(T_0001*(H0_10*H_HATT+H_TTTA*T_1000)+2.0*H_TTA*T_1001)+H_TTT*T_0001_p2*T_2000+2.0*H_TT*T_0001*T_2001
        );

    ds->df13000 += factor*(
        A_13*H_A+3.0*A_01_p2*A_11*H_AAA+H0_13*H_H+6.0*A_01*A_11*(H0_01*H_HAA+H_TAA*T_0100)+3.0*A_11*(A_02*H_AA+H0_02*H_HA+2.0*H0_01*H_HTA*T_0100+H_TTA*T_0100_p2+H_TA*T_0200)+3.0*H0_11*(A_02*H_HA+H_HTT*T_0100_p2+H_HT*T_0200)+A_03*(A_10*H_AA+H0_10*H_HA+H_TA*T_1000)+A_01_p3*(A_10*H_AAAA+H0_10*H_HAAA+H_TAAA*T_1000)+T_0300*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+H0_03*(A_10*H_HA+H_HT*T_1000)+(3.0*A_02*H_TA+3.0*H0_02*H_HT+H_TTT*T_0100_p2+3.0*H_TT*T_0200)*T_1100+3.0*A_01_p2*(H0_11*H_HAA+T_0100*(A_10*H_TAAA+H0_10*H_HTAA+H_TTAA*T_1000)+H0_01*(A_10*H_HAAA+H_HTAA*T_1000)+H_TAA*T_1100)+3.0*A_01*(A_12*H_AA+H0_12*H_HA+2.0*H0_11*H_HTA*T_0100+A_02*(A_10*H_AAA+H0_10*H_HAA+H_TAA*T_1000)+T_0200*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+T_0100_p2*(A_10*H_TTAA+H0_10*H_HATT+H_TTTA*T_1000)+H0_02*(A_10*H_HAA+H_HTA*T_1000)+2.0*H0_01*T_0100*(A_10*H_HTAA+H_HATT*T_1000)+2.0*H0_01*H_HTA*T_1100+2.0*H_TTA*T_0100*T_1100+H_TA*T_1200)+T_0100*(3.0*A_12*H_TA+3.0*H0_12*H_HT+3.0*A_02*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+3.0*T_0200*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+T_0100_p2*(A_10*H_TTTA+H0_10*H_HTTT+H_TTTT*T_1000)+3.0*H0_02*(A_10*H_HTA+H_HTT*T_1000)+2.0*H_TTT*T_0100*T_1100+3.0*H_TT*T_1200)+3.0*H0_01*(A_12*H_HA+A_02*(A_10*H_HAA+H_HTA*T_1000)+T_0200*(A_10*H_HTA+H_HTT*T_1000)+T_0100_p2*(A_10*H_HATT+H_HTTT*T_1000)+2.0*H_HTT*T_0100*T_1100+H_HT*T_1200)+H_T*T_1300 + EC_13
        );

    ds->df12100 += factor*(
        T_0210*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+(A_02*H_TA+H0_02*H_HT+A_01*(A_01*H_TAA+2.0*H0_01*H_HTA)+2.0*(A_01*H_TTA+H0_01*H_HTT)*T_0100+H_TTT*T_0100_p2+H_TT*T_0200)*T_1010+2.0*T_0110*(A_11*H_TA+H0_11*H_HT+A_01*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+T_0100*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+H0_01*(A_10*H_HTA+H_HTT*T_1000)+H_TT*T_1100)+2.0*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_1110+T_0010*(A_12*H_TA+H0_12*H_HT+A_11*(A_01*H_TAA+2.0*H0_01*H_HTA)+A_02*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+T_0200*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+T_0100_p2*(A_10*H_TTTA+H0_10*H_HTTT+H_TTTT*T_1000)+H0_02*(A_10*H_HTA+H_HTT*T_1000)+A_01*(A_11*H_TAA+2.0*H0_11*H_HTA+A_01*(A_10*H_TAAA+H0_10*H_HTAA+H_TTAA*T_1000)+2.0*H0_01*(A_10*H_HTAA+H_HATT*T_1000))+2.0*T_0100*(A_11*H_TTA+H0_11*H_HTT+A_01*(A_10*H_TTAA+H0_10*H_HATT+H_TTTA*T_1000)+H0_01*(A_10*H_HATT+H_HTTT*T_1000))+2.0*(A_01*H_TTA+H0_01*H_HTT)*T_1100+2.0*H_TTT*T_0100*T_1100+H_TT*T_1200)+H_T*T_1210
        );

    ds->df12010 += factor*(
        T_0201*(A_10*H_TA+H0_10*H_HT+H_TT*T_1000)+(A_02*H_TA+H0_02*H_HT+A_01*(A_01*H_TAA+2.0*H0_01*H_HTA)+2.0*(A_01*H_TTA+H0_01*H_HTT)*T_0100+H_TTT*T_0100_p2+H_TT*T_0200)*T_1001+2.0*T_0101*(A_11*H_TA+H0_11*H_HT+A_01*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+T_0100*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+H0_01*(A_10*H_HTA+H_HTT*T_1000)+H_TT*T_1100)+2.0*(A_01*H_TA+H0_01*H_HT+H_TT*T_0100)*T_1101+T_0001*(A_12*H_TA+H0_12*H_HT+A_11*(A_01*H_TAA+2.0*H0_01*H_HTA)+A_02*(A_10*H_TAA+H0_10*H_HTA+H_TTA*T_1000)+T_0200*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+T_0100_p2*(A_10*H_TTTA+H0_10*H_HTTT+H_TTTT*T_1000)+H0_02*(A_10*H_HTA+H_HTT*T_1000)+A_01*(A_11*H_TAA+2.0*H0_11*H_HTA+A_01*(A_10*H_TAAA+H0_10*H_HTAA+H_TTAA*T_1000)+2.0*H0_01*(A_10*H_HTAA+H_HATT*T_1000))+2.0*T_0100*(A_11*H_TTA+H0_11*H_HTT+A_01*(A_10*H_TTAA+H0_10*H_HATT+H_TTTA*T_1000)+H0_01*(A_10*H_HATT+H_HTTT*T_1000))+2.0*(A_01*H_TTA+H0_01*H_HTT)*T_1100+2.0*H_TTT*T_0100*T_1100+H_TT*T_1200)+H_T*T_1201
        );

    ds->df11200 += factor*(
        (T_0010*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)+2.0*H_TT*T_0110)*T_1010+T_0010*(2.0*T_0110*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)*T_1010+T_0010*(A_11*H_TTA+H0_11*H_HTT+A_01*(A_10*H_TTAA+H0_10*H_HATT+H_TTTA*T_1000)+T_0100*(A_10*H_TTTA+H0_10*H_HTTT+H_TTTT*T_1000)+H0_01*(A_10*H_HATT+H_HTTT*T_1000)+H_TTT*T_1100)+2.0*H_TT*T_1110)
        );

    ds->df11110 += factor*(
        T_0001*T_0110*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+H_TT*T_0110*T_1001+(T_0001*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)+H_TT*T_0101)*T_1010+T_0010*(T_0101*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)*T_1001+T_0001*(A_11*H_TTA+H0_11*H_HTT+A_01*(A_10*H_TTAA+H0_10*H_HATT+H_TTTA*T_1000)+T_0100*(A_10*H_TTTA+H0_10*H_HTTT+H_TTTT*T_1000)+H0_01*(A_10*H_HATT+H_HTTT*T_1000)+H_TTT*T_1100)+H_TT*T_1101)+H_TT*T_0001*T_1110
        );

    ds->df11020 += factor*(
        (T_0001*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)+2.0*H_TT*T_0101)*T_1001+T_0001*(2.0*T_0101*(A_10*H_TTA+H0_10*H_HTT+H_TTT*T_1000)+(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)*T_1001+T_0001*(A_11*H_TTA+H0_11*H_HTT+A_01*(A_10*H_TTAA+H0_10*H_HATT+H_TTTA*T_1000)+T_0100*(A_10*H_TTTA+H0_10*H_HTTT+H_TTTT*T_1000)+H0_01*(A_10*H_HATT+H_HTTT*T_1000)+H_TTT*T_1100)+2.0*H_TT*T_1101)
        );

    ds->df10300 += factor*(
        T_0010_p2*(T_0010*(A_10*H_TTTA+H0_10*H_HTTT+H_TTTT*T_1000)+3.0*H_TTT*T_1010)
        );

    ds->df10210 += factor*(
        T_0010*(T_0010*(T_0001*(A_10*H_TTTA+H0_10*H_HTTT+H_TTTT*T_1000)+H_TTT*T_1001)+2.0*H_TTT*T_0001*T_1010)
        );

    ds->df10120 += factor*(
        T_0001*(T_0010*(T_0001*(A_10*H_TTTA+H0_10*H_HTTT+H_TTTT*T_1000)+2.0*H_TTT*T_1001)+H_TTT*T_0001*T_1010)
        );

    ds->df10030 += factor*(
        T_0001_p2*(T_0001*(A_10*H_TTTA+H0_10*H_HTTT+H_TTTT*T_1000)+3.0*H_TTT*T_1001)
        );

    ds->df04000 += factor*(
        A_04*H_A+3.0*A_02_p2*H_AA+A_01_p4*H_AAAA+H0_04*H_H+4.0*H0_01*A_03*H_HA+4.0*A_03*H_TA*T_0100+4.0*H0_03*H_HT*T_0100+6.0*H0_02*H_HTT*T_0100_p2+4.0*H0_01*H_HTTT*T_0100_p3+H_TTTT*T_0100_p4+4.0*A_01_p3*(H0_01*H_HAAA+H_TAAA*T_0100)+6.0*H0_02*H_HT*T_0200+12.0*H0_01*H_HTT*T_0100*T_0200+6.0*H_TTT*T_0100_p2*T_0200+3.0*H_TT*T_0200_p2+6.0*A_02*(H0_02*H_HA+A_01*(A_01*H_AAA+2.0*H0_01*H_HAA)+2.0*(A_01*H_TAA+H0_01*H_HTA)*T_0100+H_TTA*T_0100_p2+H_TA*T_0200)+6.0*A_01_p2*(H0_02*H_HAA+2.0*H0_01*H_HTAA*T_0100+H_TTAA*T_0100_p2+H_TAA*T_0200)+4.0*(H0_01*H_HT+H_TT*T_0100)*T_0300+4.0*A_01*(A_03*H_AA+H0_03*H_HA+T_0100*(3.0*H0_02*H_HTA+3.0*H0_01*H_HATT*T_0100+H_TTTA*T_0100_p2)+3.0*(H0_01*H_HTA+H_TTA*T_0100)*T_0200+H_TA*T_0300)+H_T*T_0400 + EC_04
        );

    ds->df03100 += factor*(
        A_03*H_TA*T_0010+A_01_p3*H_TAAA*T_0010+H0_03*H_HT*T_0010+3.0*A_01_p2*(T_0010*(H0_01*H_HTAA+H_TTAA*T_0100)+H_TAA*T_0110)+T_0110*(3.0*A_02*H_TA+3.0*H0_02*H_HT+H_TTT*T_0100_p2+3.0*H_TT*T_0200)+3.0*A_01*(A_02*H_TAA*T_0010+H0_02*H_HTA*T_0010+2.0*H0_01*H_HATT*T_0010*T_0100+H_TTTA*T_0010*T_0100_p2+2.0*H0_01*H_HTA*T_0110+2.0*H_TTA*T_0100*T_0110+H_TTA*T_0010*T_0200+H_TA*T_0210)+T_0100*(3.0*A_02*H_TTA*T_0010+3.0*H0_02*H_HTT*T_0010+H_TTTT*T_0010*T_0100_p2+2.0*H_TTT*T_0100*T_0110+3.0*H_TTT*T_0010*T_0200+3.0*H_TT*T_0210)+3.0*H0_01*(A_02*H_HTA*T_0010+H_HTTT*T_0010*T_0100_p2+H_HTT*(2.0*T_0100*T_0110+T_0010*T_0200)+H_HT*T_0210)+H_TT*T_0010*T_0300+H_T*T_0310
        );

    ds->df03010 += factor*(
        A_03*H_TA*T_0001+A_01_p3*H_TAAA*T_0001+H0_03*H_HT*T_0001+3.0*A_01_p2*(T_0001*(H0_01*H_HTAA+H_TTAA*T_0100)+H_TAA*T_0101)+T_0101*(3.0*A_02*H_TA+3.0*H0_02*H_HT+H_TTT*T_0100_p2+3.0*H_TT*T_0200)+3.0*A_01*(A_02*H_TAA*T_0001+H0_02*H_HTA*T_0001+2.0*H0_01*H_HATT*T_0001*T_0100+H_TTTA*T_0001*T_0100_p2+2.0*H0_01*H_HTA*T_0101+2.0*H_TTA*T_0100*T_0101+H_TTA*T_0001*T_0200+H_TA*T_0201)+T_0100*(3.0*A_02*H_TTA*T_0001+3.0*H0_02*H_HTT*T_0001+H_TTTT*T_0001*T_0100_p2+2.0*H_TTT*T_0100*T_0101+3.0*H_TTT*T_0001*T_0200+3.0*H_TT*T_0201)+3.0*H0_01*(A_02*H_HTA*T_0001+H_HTTT*T_0001*T_0100_p2+H_HTT*(2.0*T_0100*T_0101+T_0001*T_0200)+H_HT*T_0201)+H_TT*T_0001*T_0300+H_T*T_0301
        );

    ds->df02200 += factor*(
        A_02*H_TTA*T_0010_p2+A_01_p2*H_TTAA*T_0010_p2+H0_02*H_HTT*T_0010_p2+2.0*H0_01*H_HTTT*T_0010_p2*T_0100+H_TTTT*T_0010_p2*T_0100_p2+4.0*H0_01*H_HTT*T_0010*T_0110+4.0*H_TTT*T_0010*T_0100*T_0110+2.0*H_TT*T_0110_p2+2.0*A_01*T_0010*(T_0010*(H0_01*H_HATT+H_TTTA*T_0100)+2.0*H_TTA*T_0110)+H_TTT*T_0010_p2*T_0200+2.0*H_TT*T_0010*T_0210
        );

    ds->df02110 += factor*(
        2.0*(T_0001*(A_01*H_TTA+H0_01*H_HTT+H_TTT*T_0100)+H_TT*T_0101)*T_0110+T_0010*(A_02*H_TTA*T_0001+H0_02*H_HTT*T_0001+A_01*(A_01*H_TTAA+2.0*H0_01*H_HATT)*T_0001+2.0*(A_01*H_TTTA+H0_01*H_HTTT)*T_0001*T_0100+H_TTTT*T_0001*T_0100_p2+2.0*(A_01*H_TTA+H0_01*H_HTT)*T_0101+2.0*H_TTT*T_0100*T_0101+H_TTT*T_0001*T_0200+H_TT*T_0201)+H_TT*T_0001*T_0210
        );

    ds->df02020 += factor*(
        A_02*H_TTA*T_0001_p2+A_01_p2*H_TTAA*T_0001_p2+H0_02*H_HTT*T_0001_p2+2.0*H0_01*H_HTTT*T_0001_p2*T_0100+H_TTTT*T_0001_p2*T_0100_p2+4.0*H0_01*H_HTT*T_0001*T_0101+4.0*H_TTT*T_0001*T_0100*T_0101+2.0*H_TT*T_0101_p2+2.0*A_01*T_0001*(T_0001*(H0_01*H_HATT+H_TTTA*T_0100)+2.0*H_TTA*T_0101)+H_TTT*T_0001_p2*T_0200+2.0*H_TT*T_0001*T_0201
        );

    ds->df01300 += factor*(
        T_0010_p2*(T_0010*(A_01*H_TTTA+H0_01*H_HTTT+H_TTTT*T_0100)+3.0*H_TTT*T_0110)
        );

    ds->df01210 += factor*(
        T_0010*(T_0010*(T_0001*(A_01*H_TTTA+H0_01*H_HTTT+H_TTTT*T_0100)+H_TTT*T_0101)+2.0*H_TTT*T_0001*T_0110)
        );

    ds->df01120 += factor*(
        T_0001*(T_0010*(T_0001*(A_01*H_TTTA+H0_01*H_HTTT+H_TTTT*T_0100)+2.0*H_TTT*T_0101)+H_TTT*T_0001*T_0110)
        );

    ds->df01030 += factor*(
        T_0001_p2*(T_0001*(A_01*H_TTTA+H0_01*H_HTTT+H_TTTT*T_0100)+3.0*H_TTT*T_0101)
        );

    ds->df00400 += factor*(
        H_TTTT*T_0010_p4
        );

    ds->df00310 += factor*(
        H_TTTT*T_0001*T_0010_p3
        );

    ds->df00220 += factor*(
        H_TTTT*T_0001_p2*T_0010_p2
        );

    ds->df00130 += factor*(
        H_TTTT*T_0001_p3*T_0010
        );

    ds->df00040 += factor*(
        H_TTTT*T_0001_p4
        );

/* End (Finally!) */
}
#endif
