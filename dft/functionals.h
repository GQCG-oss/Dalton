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

typedef double real;

/* FirstDrv: matrix of first order derivatives with respect to two
 * parameters: density rho and SQUARE of the gradient of density grho.  
 * zeta_i = |\nabla\rho_i|²
 * mu     = |\nabla\rho_\alpha||\nabla\rho_\beta|
 */
typedef struct FunFirstFuncDrv_ {
    real df1000;  /* d/drho F     */
    real df0100;
    real df0010;  /* d/zeta F     */
    real df0001;
    real df00001;
} FunFirstFuncDrv;

/* SecondDrv:  matrix  of  second  order functional  derivatives  with
 * respect  to two  parameters: density  rho and  SQUARE  of the
 * density gradient zeta.  The derivatives are computed for alpha-alpha
 * or beta-beta spin-orbital block (i.e. include triplet flag).
 */
/* SecondFuncDrv: this structure is used by functional derivative
 * evaluation procedures. Do not include "triplet" transformation.
 */
typedef struct FunSecondFuncDrv_ {
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
} FunSecondFuncDrv;


/* ThirdFuncDrv: matrix of third derivatives with respect to five
   parameters: density rho_alpha and SQUARE of the density gradient
   zeta.  and mu.
*/

typedef struct FunThirdFuncDrv_ {
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
} FunThirdFuncDrv;
typedef struct Functional_ Functional;

enum FunError { FUN_OK, FUN_UNKNOWN, FUN_CONF_ERROR };
enum FunError fun_select_by_name(const char *conf_string);
extern Functional *selected_func;
extern int (*fun_printf)(const char *fmt, ...);
extern void (*fun_set_hf_weight)(real w);
extern real (*fun_get_hf_weight)(void);
extern void (*fun_set_cam_param)(real w, real b);

/* FunDensProp structure contains properties of the density that are
   needed for functional evaluation and possibly other purposes.
*/
typedef struct FunDensProp_ {
    real rhoa,  rhob;
    real grada, gradb; /* norms of the density gradient, not squares */
    real gradab;       /* scalar product of grada and gradb */
    /* real current[3] or something may come in the future :-) */
} FunDensProp;

/* EnergyFunc: the function returning the energy for given densities
   and gradients. Note that some functionals(like LYP) depend explicitely
   on separately alpha and beta densities
*/
typedef int (*IsGGAFunc)(void);
typedef int (*ReadInputFunc)(const char* conf_string);
typedef void (*ReportFunc)(void);
typedef real (*EnergyFunc)(const FunDensProp* dens_prop);
typedef void (*FirstOrderFun)(FunFirstFuncDrv *ds, real factor,
                              const FunDensProp* dns_prp);

typedef void (*SecondOrderFun)(FunSecondFuncDrv *ds, real factor,
                               const FunDensProp* dens_prop);

typedef void (*ThirdOrderFun)(FunThirdFuncDrv *ds, real factor,
                              const FunDensProp* dens_prop);

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

void drv1_clear(FunFirstFuncDrv* gga);  /* set all components to 0 */
void drv2_clear(FunSecondFuncDrv* gga); /* set all components to 0 */
void drv3_clear(FunThirdFuncDrv* gga);  /* set all components to 0 */

/* The list of functionals */
/* sorted list of generic functionals */
extern Functional BeckeFunctional;
extern Functional Example2Functional;
extern Functional ExampleFunctional;
extern Functional KTFunctional;
extern Functional LB94Functional;
extern Functional LYPFunctional;
extern Functional OPTXFunctional;
extern Functional P86cFunctional;
extern Functional PW86xFunctional;
extern Functional Pw91cFunctional;
extern Functional PZ81Functional;
extern Functional PbecFunctional;
extern Functional PbexFunctional;
extern Functional SlaterFunctional;
extern Functional VWN3Functional;
extern Functional VWN5Functional;
extern Functional VWNIFunctional;
extern Functional VWNFunctional;
extern Functional XAlphaFunctional;

/* sorted list of mixed functionals */
extern Functional B3LYPFunctional;
extern Functional B3LYPGaussFunctional;
extern Functional B3P86Functional;
extern Functional B3P86GFunctional;
extern Functional BLYPFunctional;
extern Functional BP86Functional;
extern Functional BPW91Functional;
extern Functional Camb3lypFunctional;
extern Functional GGAKeyFunctional;
extern Functional KT1Functional;
extern Functional KT2Functional;
extern Functional KT3Functional;
extern Functional LDAFunctional;
extern Functional OLYPFunctional;
extern Functional PBE0Functional;
extern Functional PBEFunctional;
extern Functional SVWN3Functional;
extern Functional SVWN5Functional;

/* the list of the functionals */
extern Functional* available_functionals[];

extern int fun_true(void);
extern int fun_false(void);
/* fortran (and not only) functional stub routines */
void dftlistfuncs_(void);
int dft_isgga_(void);
int dft_isgga__(void);

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
