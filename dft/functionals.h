#if defined(__CVERSION__) /* THE C VERSION OF THE HEADERS */
/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* The variables, structures and functions related to computation
   of functional and their derivatives.
   (c) Pawel Salek, pawsa@theochem.kth.se. 2001.07.13

   NOTE1: the derivatives are computed with respect to the density,
   and SQUARE of the density gradient. This is a choice. It stems from
   the fact that the factors involved in the derivative vector
   distribution depend on the square of the density gradient.

   NOTE2: C version is included once per file, Fortran version -
   multiple times.
*/
#ifndef _FUNCTIONALS_H_
#define _FUNCTIONALS_H_

#include "general.h"

/* FirstDrv: matrix of first order derivatives with respect to two
 * parameters: density rho and SQUARE of the gradient of density grho.  
 * zeta_i = |\nabla\rho_i|²
 * mu     = |\nabla\rho_\alpha||\nabla\rho_\beta|
 */

struct FirstDrv_ {
    real fR;  /* d/drho F     */
    real fZ;  /* d/zeta F     */
};

/* FirstDrv: matrix of first order derivatives with respect to two
 * parameters: density rho and SQUARE of the gradient of density grho.  
 * zeta_i = |\nabla\rho_i|²
 * mu     = |\nabla\rho_\alpha||\nabla\rho_\beta|
 */
struct FirstFuncDrv_ {
    real df1000;  /* d/drho F     */
    real df0100;
    real df0010;  /* d/zeta F     */
    real df0001;
    real df00001;
};

/* SecondDrv:  matrix  of  second  order functional  derivatives  with
 * respect  to two  parameters: density  rho and  SQUARE  of the
 * density gradient zeta.  The derivatives are computed for alpha-alpha
 * or beta-beta spin-orbital block (i.e. include triplet flag).
 */
struct SecondDrv_ {
    real fR; /* d/drho  F */
    real fZ; /* d/dzeta F */ 
    real fRR; /* d/drho² F */
    real fRZ; /* d/(drho dzeta) F */ 
    real fZZ; /* d/dzeta² F */ 
    /* additional derivatives required by  */ 
    /* general linear response scheme     */ 
    real fRG; /* d/(drho dgamma) F */ 
    real fZG; /* d/(dzeta dgamma) F */
    real fGG; /* d/dzgamma² F */  
    real fG;  /* d/dgamma F */ 
};

/* SecondFuncDrv: this structure is used by functional derivative
 * evaluation procedures. Do not include "triplet" transformation.
 */
struct SecondFuncDrv_ {
    real df1000;  /* d/drho_alpha F               */
    real df0100;  /* d/drho_beta F                */
    real df0010;  /* d/|zeta_alpha| F             */
    real df0001;  /* d/|zeta_beta| F              */
    real df00001;
    real df2000;  /* d/drho_alpha^2 F             */
    real df1100;  /* d/(drho_alpha drho_beta) F   */
    real df1010;  /* d/drho_alpha d/dzeta_alpha F */
    real df1001;  /* d/drho_alpha d/dzeta_beta F  */
    real df10001;
    real df0200;  /* d/drho_beta^2 F              */
    real df0110;  /* d/drho_beta d/dzeta_alpha F  */ 
    real df0101;  /* d/drho_beta d/dzeta_beta F   */
    real df01001;
    real df0020;  /* d/dzeta_alpha^2 F          */
    real df0011;  /* d2/dzeta_a zeta_b F          */
    real df00101;
    real df0002;  /* d/dzeta_beta^2 F             */
    real df00011;
    real df00002;
};

/* ThirdDrv: matrix of third derivatives with respect to two parameters:
   density rho and SQUARE of the density gradient zeta.
*/
struct ThirdDrv_ {
    real fR;   /* d/drho  F */
    real fZ;   /* d/dzeta F */
    real fG;   /* d/dgamma F */
    real fRR[2];  /* d/drho² F */
    real fRZ[2];  /* d/(drho dzeta) F */
    real fZZ[2];  /* d/dzeta² F */
    real fRG[2];  /* d/(drho dgamma) F */
    real fRRR[2]; /* d/drho³ F */
    real fRRZ[2][2]; /* d/(drho² dzeta) F */
    /* two forms of fRRG needed as the formulae is non symmetric */
    real fRRG[2];     /* d/(drho? dgamma) F */
    real fRRGX[2][2]; /* d/(drho? dgamma) F */   
    real fRZZ[2][2]; /* d/(drho dzeta²) F */
    real fZZZ[2]; /* d/dzeta³ F */
};

/* ThirdFuncDrv: matrix of third derivatives with respect to five
   parameters: density rho_alpha and SQUARE of the density gradient
   zeta.  and mu.
*/

struct ThirdFuncDrv_ {
    real df1000;   /* d/drho F          */
    real df0100; 
    real df0010;   /* d/|zeta| F        */
    real df0001;
    real df00001;
     
    real df2000;  /* d/drho_alpha^2 F  */
    real df0200; 
    real df1100;  /* d/(drho_alpha drho_beta) F   */
    real df1010;  /* d/drho_alpha d/dzeta_alpha F */
    real df1001;  /* d/drho_alpha d/dzeta_beta F  */
    real df0101;
    real df0110;
    real df10001;
    real df01001;   
    real df0020;  /* d/dzeta^2 F                  */
    real df0002;
    real df0011;  /* d2/dzeta_a zeta_b F          */

    real df3000; /* d/drho_alpha^3 F              */
    real df0300; 
    real df2100; /* d/d(rho_alpha^2rho_beta) F    */
    real df1200;
    real df2010; /* d/drhoa^2 d/dzetaa F */
    real df2001;
    real df0201; 
    real df0210;     
    real df1110; /* d/drho^2 d/dzeta F */
    real df1101;
    real df1020; /* d/drho d/dzeta^2 F */
    real df1002;
    real df0102;
    real df0120; 
    real df0030; /* d/dzeta^3 F        */
    real df0003;
    real df1011;
    real df0111;
    real df0021;
    real df0012;
    real df20001;
    real df02001;
    real df11001;  
  /* most derivatives for fith variable not included *
   * this makes functionals valid only to square grad rho 
   */
};

/* EnergyFunc: the function returning the energy for given densities
   and gradients. Note that some functionals(like LYP) depend explicitely
   on separately alpha and beta densities
*/
typedef int (*IsGGAFunc)(void);
typedef int (*ReadInputFunc)(const char* conf_string);
typedef void (*ReportFunc)(void);
typedef real (*EnergyFunc)(const DftDensProp* dens_prop);
typedef void (*FirstOrderFun)(FirstFuncDrv *ds, real factor,
                              const DftDensProp* dns_prp);

typedef void (*SecondOrderFun)(SecondFuncDrv *ds, real factor,
                               const DftDensProp* dens_prop);

typedef void (*ThirdOrderFun)(ThirdFuncDrv *ds, real factor,
                              const DftDensProp* dens_prop);

struct Functional_ {
    const char* name; /* descriptive functional name (usually 5 characters) */
    IsGGAFunc       is_gga;
    ReadInputFunc   read;
    ReportFunc      report;
   /* Only unrestricted implementations are needed. A benchmark for
     * a CO molecule with 28 basis function reveals a 4% time difference.
     * This difference will only decrease for larger systems. */
    EnergyFunc      func;
    FirstOrderFun   first;
    SecondOrderFun  second;
    ThirdOrderFun   third;
};

int  dft_isgga(void);
void dftreport_(void);
void dftlistfuncs_(void);
real dftenergy_(const real* rho, const real* grad);
void dftpot0_(FirstDrv *ds, const real* w, const DftDensProp* dp);
void dftpot1_(SecondDrv*ds, const real* w, const DftDensProp* dp,
              const int* triplet); 
void dftpot2_(ThirdDrv *ds, real factor, const DftDensProp* dp, int isgga,
              int triplet);

void drv1_clear(FirstFuncDrv* gga);  /* set all components to 0 */
void drv2_clear(SecondFuncDrv* gga); /* set all components to 0 */
void drv3_clear(ThirdFuncDrv* gga);  /* set all components to 0 */
void dft_set_hf_weight(real w);

/* The list of functionals */
/* generic functionals */
extern Functional BeckeFunctional;
extern Functional Example2Functional;
extern Functional Example3Functional;
extern Functional Example4Functional;
extern Functional Example5Functional;
extern Functional Example6Functional;
extern Functional Example7Functional;
extern Functional Example8Functional;
extern Functional Example9Functional;
extern Functional ExampleFunctional;
extern Functional KTFunctional;
extern Functional LB94Functional;
extern Functional LYPFunctional;
extern Functional OPTXFunctional;
extern Functional P86cFunctional;
extern Functional PbecFunctional;
extern Functional PbexFunctional;
extern Functional PW86xFunctional;
extern Functional PW91cFunctional;
extern Functional PW92Functional;
extern Functional PW92peFunctional;
extern Functional PWggaIIc2Functional;
extern Functional PWggaIIcFunctional;
extern Functional PWggaIIxFunctional;
extern Functional PZ81Functional;
extern Functional SlaterFunctional;
extern Functional VWN3Functional;
extern Functional VWN5Functional;
extern Functional VWNFunctional;
extern Functional XAlphaFunctional;

/* mixed functionals */
extern Functional LDAFunctional;
extern Functional LDAGaussFunctional;
extern Functional GGAKeyFunctional;
extern Functional BLYPFunctional;
extern Functional B3LYPFunctional;
extern Functional B3LYPGaussFunctional;
extern Functional BP86Functional;
extern Functional B3P86Functional;
extern Functional Camb3lypFunctional;
extern Functional KT1Functional;
extern Functional KT2Functional;
extern Functional KT3Functional;
extern Functional OLYPFunctional;
extern Functional PBE0Functional;
extern Functional PBEFunctional;
extern Functional PW91Functional;
extern Functional SVWN5Functional;
extern Functional SVWN3Functional;

/* the list of the functionals */
extern Functional* available_functionals[];

int fun_true(void);
int fun_false(void);
extern void dftsethf_(real *);
extern real dftgethf_(void);
extern void dftsetcam_(real *, real *);
#define dft_get_hf_weight() dftgethf_()
#define dft_set_hf_weight(w) do{real x=w;dftsethf_(&x);}while(0)
#define dft_set_cam_param(w,b) do{real x=w,be=b;dftsetcam_(&x, &be);}while(0)
#endif /* _FUNCTIONALS_H_ */
#else /* THE FORTRAN VERSION OF THE HEADERS; out of date as of 20021120 */

      LOGICAL DFT_ISGGA
      EXTERNAL DFT_ISGGA
c     The derivatives with respect to rho and zeta
      INTEGER FR0, FZ0, FRR, FRZ, FZZ, FRRR, FRRZ, FRZZ, FZZZ
      INTEGER DERVS1, DERVS2, DERVS3
      PARAMETER(FR0=1,  FZ0=2,   FRR=3,  FRZ=4)
      PARAMETER(FZZ=5,  FRRR=6,  FRRZ=7, FRZZ=8, FZZZ=9)
      PARAMETER(DERVS1=3, DERVS2=5, DERVS3=9)

c     The derivatives are a vector that can be declared as
c     DIMENSION DRVS(DERVS3)
c     where 3 stands for the functional derivative order.
c     Second order derivatives can be declared as
c     DIMENSION DRVS(DERVS2)
c     etc.
c     The derivatives can be accessed as DRVS(F1000), DRVS(F0010) etc.
c
      INTEGER F1000, F0010, F2000, F1100, F1010, F1001, F0020
      INTEGER F3000, F2100, F2010, F1110, F2001
      INTEGER F1020, F0120, F0030, F1011, F0021
      INTEGER FUN_DERVS1, FUN_DERVS2, FUN_DERVS3
      PARAMETER(F1000=1, F0010=2, F2000=3, F1100=4, F1010=5, F1001=6)
      PARAMETER(F0020=7)
      PARAMETER(F3000=8, F2100=9, F2010=10,F1110=11,F2001=12)
      PARAMETER(F1020=13,F0120=14,F0030=15,F1011=16,F0021=17)
      PARAMETER(FUN_DERVS1=2, FUN_DERVS2=7, FUN_DERVS3=17)
#endif /* __CVERSION__ */
