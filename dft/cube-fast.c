/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* quadra-strict.c:
   (c) Brano Jansik, brano@theochem.kth.se, Nov 2003

   The DFT cubic response subroutines that follow string matrix
   commutator expressions. 

   This file also demonstrates usage of DFT integrator with a single
   callback only. While multiple callback can be easily emulated by a
   single one with a help of a wrapper routine, the extra flexibility
   is a nice feature.

   ONE THING TO REMEMBER: stick to separate treatment of alpha and
   beta densities. Really do.

   NOTE: 
   The code can add the fock-like contribution [k,[k,omega]] if needed.

   For example QFOCK returns complete FOCK/Kohn-Sham matrix, including
   DFT contribution, as opposed to RSPFXD/RSPFX which return
   electronic part only (w/o DFT contribution). If this is the case
   (i.e. .NOT.DIRFCK) we add the DFT contribution to KS matrix later
   in DFTQRC. field addfock determines, if the fock-type contribution
   should be added or not.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define __CVERSION__

/* temporary code */
#undef  TIMING
#undef  TEST

#include "integrator.h"
#include "functionals.h"
#include "general.h"

#include "inforb.h"
#include "priunit.h" 

#if !defined __inline__
/* inline some stuff whenever possible */
#define __inline__
#endif

/* derivatives prefactors */
#define RRRR 8.0
#define RRRZ 16.0
#define RRZZ 32.0
#define RZZZ 64.0
#define ZZZZ 128.0
#define RRR 4.0
#define RRZ 8.0
#define RZZ 16.0
#define ZZZ 32.0
#define RR 2.0
#define RZ 4.0
#define ZZ 8.0
#define R 1.0
#define Z 2.0

static void
matn_times_vec_full(real PM, int sym, const real* mat, const real* vec, real a, real* res);
static void
matt_times_vec_full(real PM, int sym, const real* mat, const real* vec, real a, real* res);
static void
matn_times_vec_symm(real PM, int sym, const real* mat, const real* vec, real a, real* res);
static void
matt_times_vec_symm(real PM, int sym, const real* mat, const real* vec, real a, real* res);


typedef struct {
    /* Commutator structure
     (out of these data we construct correponding commutator matrix) */

    /* Each commutator matrix is consturcted as COM = fac*A*B *
     * using add_commutator(CommData *X, real pref, real *Y)  *
     * where Y = Y + pref*COM(X)                              */
    
     real **A;     /* Start of A map */
     real **B;     /* Start of B map */
     real *coeAB;  /* Coefficient of particular AB term */ 
     real fac;     /* Whole Commutator matrix will be scaled by this factor fac */
     real pref;    /* Whole Commutator matrix will be added with this weigth    */
     int  vec_num; /* Number of elements (vectors or columns) in maps A, B */
     int  sym;     /* Commutator symmetry */

} CommData;

static void
add_commutators_full(CommData **d, int maxlst, real *dftcontr);
static void
add_commutators_ij(CommData **d, int maxlst, real *dftcontr);

#if 0
static void
add_commutators_symm(CommData **d, int maxlst, real *dftcontr);
#endif

static const real PLUS = 1.0, MINUS = -1.0, NEW = 0.0, ADD = 1.0;
static const real HALFR = 0.5, SIXTHR = 1.0/6.0;

/* the temporary data used for integration */
struct CubeFastData_ {
    /* pointers to external data (data is not owned) */
    const real *kappaB, *kappaC, *kappaD;
    int           symB,    symC,    symD;

    /* dft contribution matrix */
    real* dftcontr;

    /* Omega and its commutators (density dependent part) */
    CommData omega; /* omega = phi_p*phi_q */
    CommData omB;   /* [kappaB, omega] (commutator data) */
    CommData omC;   /* [kappaC, omega] (commutator data) */
    CommData omD;   /* [kappaD, omega] (commutator data) */
    CommData omBC;  /* [kappaB,[kappaC,omega]] + [kappaC,[kappaB,omega]] */
    CommData omCD;  /* [kappaC,[kappaD,omega]] + [kappaD,[kappaC,omega]] */
    CommData omBD;  /* [kappaB,[kappaD,omega]] + [kappaD,[kappaB,omega]]*/
    CommData omBCD; /* sum of all six permutations [kappaB,[kappaC,[kappaD,omega]]] */
  
    /* Omegax and its commutators */
    /* omega[x] = a*b[x] + b[x]*a; first index is a*b[x], b[x]*a contributions,
       second index are x component                                          */
    CommData omegax[2][3]; /* omega[x] = phi_p*dphi_q/dx + dphi_p/dx*phi_q */
    CommData omBx[2][3];   /* [kappaB, omega(x)] (commutator data) */
    CommData omCx[2][3];   /* [kappaC, omega(x)] (commutator data) */
    CommData omDx[2][3];   /* [kappaD, omega(x)] (commutator data) */
    CommData omBCx[2][3];  /* [kappaB,[kappaC,omega(x)]] + [kappaC,[kappaB,omega(x)]] */
    CommData omCDx[2][3];  /* [kappaC,[kappaD,omega(x)]] + [kappaD,[kappaC,omega(x)]] */
    CommData omBDx[2][3];  /* [kappaB,[kappaD,omega(x)]] + [kappaD,[kappaB,omega(x)]]*/
    CommData omBCDx[2][3]; /* sum of all six permutations [kappaB,[kappaC,[kappaD,omega]]] */

    /* Auxaliary variables used in GGA prefactors */
    real sB;   real sC;   real sD;
    real sBC;  real sCD;  real sBD;
    real sBCD;
    real pBC;  real pCD;  real pBD;
    real pBCD;

    /* Expectation values */
    real roB;
    real roC;
    real roD;

    real roBC;
    real roCD;
    real roBD;

    real roBCD;

    /* rox = rox[0], roy=rox[1], roz = rox[2] */
    real rox[3];

    real roBx[3];
    real roCx[3];
    real roDx[3];

    real roBCx[3];
    real roCDx[3];
    real roBDx[3];

    real roBCDx[3];

    /* Matrix vector multiplication scheme */
    void (*matn_times_vec)();
    void (*matt_times_vec)();
    char mm_code_version[24];

    /* Matrix distribution scheme */
    void (*add_commutators)();
    
    /* temporary space */
    real *mapA[27];
    real *mapB[27];
    real *mapAx[162];
    real *mapBx[162];
    real *vecA;
    real *vecB;
    real *vecAx[3];
    real *vecBx[3];

#ifdef TEST
    real *tmp_f;
    real *tmp_s;
#endif
    

#ifdef TIMING
    /* timing */
    real timedata[3];
#endif

};

typedef struct CubeFastData_ CubeFastData;

static CubeFastData*
cubefast_data_new(  const real *kB, 
                    const real *kC,    
		    const real *kD,
                    const int symB,
                    const int symC,
                    const int symD)  
{
/* This routine allocates space for CubeFastData structure and
   initialize its fields */

    /* Local vars */
    int x; int ab;
    int ksymBC, ksymCD, ksymBD, ksymBCD;
    int vec_size;
    int offset = 0;
    real **mapA;
    real **mapB;
    real *buf;

    /* AB coefficients */
    static real coeX[2] = {+1.0, -1.0};
    static real coeXY[4] = {-2.0,+1.0,-2.0,+1.0};
    static real coeXYZ[8] = {-1.0,+3.0,+3.0,+3.0,-3.0,-3.0,-3.0,+1.0}; 
	
    /* Create empty structure */
    CubeFastData *res = calloc(1,sizeof(CubeFastData));

    /* Assign kappa response vectors and symmetries */
    res->kappaB = kB;     res->kappaC = kC;	res->kappaD = kD;
    res->symB = symB-1;     res->symC = symC-1;     res->symD = symD-1;


   /* perhaps this should be allocated by DALTOM GETMEM */ 
    res->dftcontr = alloc_mat_MO(1);

#ifdef TEST
    res->tmp_f = alloc_mat_MO(1);
    res->tmp_s = alloc_mat_MO(1);
#endif

    /* set matrix vector multiplicaton scheme */
//    if (inforb_.nsym==1 || inforb_.norbt<31) {
    if (inforb_.nsym==1) {
     res->matn_times_vec = &matn_times_vec_full;
     res->matt_times_vec = &matt_times_vec_full;
     res->add_commutators = &add_commutators_full;
     sprintf(res->mm_code_version,"full");
    } else {
     res->matn_times_vec = &matn_times_vec_symm;
     res->matt_times_vec = &matt_times_vec_symm;
     res->add_commutators = &add_commutators_full;
//     res->add_commutators = &add_commutators_ij;
     sprintf(res->mm_code_version,"symmetry adapted");
    }

    /* allocating temporary space */ 
    /* there is (7*2)*4 vectors for omega commutators */
    /* This should be done with only one malloc */
    vec_size = sizeof(real)*inforb_.norbt;
    buf  = dal_malloc(vec_size*56);
    
    res->vecA = buf;

    offset += 7;
    res->vecB = &buf[offset*inforb_.norbt];

    /* there is 2*(2*3*16) vectors for omega(x) commutators */
      for (x=0; x<3; x++) {
        offset +=7;
        res->vecAx[x] = &buf[offset*inforb_.norbt];
        offset +=7;
        res->vecBx[x] = &buf[offset*inforb_.norbt];
      }


    /* assigning space for omega and commutators */
    /* symmetries (kommutator)sym ksymXXX */
    ksymBC  = inforb_.muld2h[res->symB][res->symC]-1;
    ksymCD  = inforb_.muld2h[res->symC][res->symD]-1;
    ksymBD  = inforb_.muld2h[res->symB][res->symD]-1;
    ksymBCD = inforb_.muld2h[res->symB][ksymCD]-1;

    /* SetComm macro */
#define SetComm(x,AA,BB,cAB,ssym,ffac,vvec_num) x.A = AA; x.B = BB; \
    x.fac = ffac; x.vec_num = vvec_num; x.sym = ssym; x.coeAB = cAB

    /* just s shortcut */
    mapA = res->mapA;
    mapB = res->mapB;

    SetComm(res->omega, &mapA[0], &mapB[0], coeX, 0,  1.0, 1);

    SetComm(res->omB,  &mapA[1],  &mapB[1], coeX, res->symB,  1.0, 2);
    SetComm(res->omC,  &mapA[3],  &mapB[3], coeX, res->symC,  1.0, 2);  
    SetComm(res->omD,  &mapA[5],  &mapB[5], coeX, res->symD,  1.0, 2); 

    SetComm(res->omBC, &mapA[7],  &mapB[7], coeXY, ksymBC, 0.5, 4);
    SetComm(res->omCD, &mapA[11], &mapB[11],coeXY, ksymCD, 0.5, 4);
    SetComm(res->omBD, &mapA[15], &mapB[15],coeXY, ksymBD, 0.5, 4);

    SetComm(res->omBCD, &mapA[19], &mapB[19], coeXYZ, ksymBCD, SIXTHR, 8);

    /* assigning space for omegax */
    for (ab=0; ab<2; ab++) {
      for (x=0; x<3; x++) {

       mapA=&res->mapAx[ab*81 + x*27];
       mapB=&res->mapBx[ab*81 + x*27];

       SetComm(res->omegax[ab][x],&mapA[0],&mapB[0],coeX,0, 1.0, 1);
       
       SetComm(res->omBx[ab][x],  &mapA[1],   &mapB[1], coeX, res->symB, 1.0, 2);
       SetComm(res->omCx[ab][x],  &mapA[3],   &mapB[3], coeX, res->symC, 1.0, 2);  
       SetComm(res->omDx[ab][x],  &mapA[5],   &mapB[5], coeX, res->symD, 1.0, 2); 

       SetComm(res->omBCx[ab][x], &mapA[7],   &mapB[7],  coeXY, ksymBC, 0.5, 4);
       SetComm(res->omCDx[ab][x], &mapA[11],  &mapB[11], coeXY, ksymCD, 0.5, 4);
       SetComm(res->omBDx[ab][x], &mapA[15],  &mapB[15], coeXY, ksymBD, 0.5, 4);

       SetComm(res->omBCDx[ab][x], &mapA[19], &mapB[19], coeXYZ, ksymBCD, SIXTHR, 8);

      }
    }

    return res;
}

static void 
cubefast_data_free(CubeFastData* tmp)
{
 /* This routine will free the CubeFastData structure */
    free(tmp->dftcontr);
    
    free(tmp->vecA);

#ifdef TEST
    free(tmp->tmp_f);
    free(tmp->tmp_s);
#endif

    free(tmp);
}

static void
matn_times_vec_full(real PM, int sym, const real* mat, const real* vec, real a, real* res)
{
    dgemv_("N", &inforb_.norbt, &inforb_.norbt, &PM, mat,
           &inforb_.norbt, vec, &ONEI, &a, res, &ONEI);
}

static void
matt_times_vec_full(real PM, int sym, const real* mat, const real* vec, real a, real* res)
{
    dgemv_("T", &inforb_.norbt, &inforb_.norbt, &PM, mat,
           &inforb_.norbt, vec, &ONEI, &a, res, &ONEI);
}

static void
matn_times_vec_symm(real PM, int ops,const real* mat, const real* vec, real a, real* res)
{
    int isym;
                                                                                                              
    for(isym=0; isym<inforb_.nsym; isym++) {
        int iorbs = inforb_.iorb[isym];
        integer noccs = inforb_.nocc[isym];
        integer nvirs = inforb_.nvir[isym];
        int i     = inforb_.muld2h[ops][isym]-1;
        int iorbi = inforb_.iorb[i];
        integer nocci = inforb_.nocc[i];
        integer nviri = inforb_.nvir[i];
        int begll = (iorbi+nocci)+iorbs        *inforb_.norbt;
        int begur = iorbi        +(iorbs+noccs)*inforb_.norbt;
        if(nviri>0) {
            if(noccs>0)
                dgemv_("N", &nviri,&noccs, &PM, &mat[begll],&inforb_.norbt,
                       &vec[iorbs], &ONEI, &a, &res[iorbi+nocci], &ONEI);
            else if(a==0) dzero_(&res[iorbi+nocci], &nviri);
        }
        if(nocci>0) {
            if(nvirs>0)
                dgemv_("N", &nocci,&nvirs, &PM, &mat[begur],&inforb_.norbt,
                       &vec[iorbs+noccs], &ONEI, &a, &res[iorbi], &ONEI);
            else if(a==0) dzero_(&res[iorbi], &nocci);
        }
    }
}

static void
matt_times_vec_symm(real PM, int ops,const real* mat, const real* vec, real a, real* res)
{
    int isym;
                                                                                                              
    for(isym=0; isym<inforb_.nsym; isym++) {
        int iorbs = inforb_.iorb[isym];
        integer noccs = inforb_.nocc[isym];
        integer nvirs = inforb_.nvir[isym];
        int i     = inforb_.muld2h[ops][isym]-1;
        int iorbi = inforb_.iorb[i];
        integer nocci = inforb_.nocc[i];
        integer nviri = inforb_.nvir[i];
        int begll = (iorbi+nocci)+iorbs        *inforb_.norbt;
        int begur = iorbi        +(iorbs+noccs)*inforb_.norbt;
        if(noccs>0) {
            if(nviri>0)
                dgemv_("T", &nviri,&noccs, &PM, &mat[begll],&inforb_.norbt,
                       &vec[iorbi+nocci], &ONEI, &a, &res[iorbs], &ONEI);
            else if(a==0) dzero_(&res[iorbs], &noccs);
        }
        if(nvirs>0) {
            if(nocci>0)
                dgemv_("T", &nocci,&nvirs, &PM, &mat[begur],&inforb_.norbt,
                       &vec[iorbi], &ONEI, &a, &res[iorbs+noccs], &ONEI);
            else if(a==0) dzero_(&res[iorbs+noccs], &nvirs);
        }
    }
}

static real inactive_trace_AB(CommData *d)
/* This function returns inactive trace of commutator d */
{
   real **mata = d->A;
   real **matb = d->B;
   real *coef   = d->coeAB;
   int vec_num = d->vec_num;
   int vec;

   /* inactive trace part */
   int isym, iorb;
   real result = 0.0;

   for (isym=0; isym<inforb_.nsym; isym++) {
      if (inforb_.nocc[isym] > 0) {
        iorb = inforb_.iorb[isym];
        for (vec=0; vec<vec_num; vec++) result +=
          coef[vec]*ddot_(&inforb_.nocc[isym], &mata[vec][iorb], &ONEI, &matb[vec][iorb], &ONEI);
      }
   }
   return d->fac*result;
}

static void
eval_omega_vectors(CubeFastData *t,
	       	const real *KB, const real *KC, const real *KD,
                real *moa, real *mob,
		real *vecA, real *vecB)
{
/* This routine returns commutator vectors A and B, see comments below */
    int  symB = t->symB;
    int  symC = t->symC;
    int  symD = t->symD;
    
    void (*matn_times_vec)() = t->matn_times_vec;
    void (*matt_times_vec)() = t->matt_times_vec;

    real *kb_mo, *kc_mo, *kd_mo;  
    real *kbT_mo, *kcT_mo, *kdT_mo;

    real *kb_kc_mo, *kc_kd_mo, *kb_kd_mo;
    real *kbT_kcT_mo, *kcT_kdT_mo, *kbT_kdT_mo;
 
    real *kb_kc_kd_mo;

    real *kbT_kcT_kdT_mo;

    /* assign vectors to tmpA memory */
    kb_mo           = &(vecA[0*inforb_.norbt]);
    kc_mo           = &(vecA[1*inforb_.norbt]);
    kd_mo           = &(vecA[2*inforb_.norbt]);
  
    kb_kc_mo        = &(vecA[3*inforb_.norbt]); 
    kc_kd_mo        = &(vecA[4*inforb_.norbt]); 
    kb_kd_mo        = &(vecA[5*inforb_.norbt]); 

    kb_kc_kd_mo     = &(vecA[6*inforb_.norbt]);


    /* assign vectors to tmpB memory (same way as tmpA) */
    kbT_mo           = &(vecB[ 0*inforb_.norbt]);
    kcT_mo           = &(vecB[ 1*inforb_.norbt]);
    kdT_mo           = &(vecB[ 2*inforb_.norbt]);

    kbT_kcT_mo       = &(vecB[ 3*inforb_.norbt]); 
    kcT_kdT_mo       = &(vecB[ 4*inforb_.norbt]);
    kbT_kdT_mo       = &(vecB[ 5*inforb_.norbt]);

    kbT_kcT_kdT_mo   = &(vecB[ 6*inforb_.norbt]);
    
    /* Computing vectors (matrix vector product) */

    /* kx_mo vectors */
    matn_times_vec(PLUS,symB,KB,moa,NEW,kb_mo);
    matn_times_vec(PLUS,symC,KC,moa,NEW,kc_mo);
    matn_times_vec(PLUS,symD,KD,moa,NEW,kd_mo);
   
    /* kxT_mo vectors */
    matt_times_vec(PLUS,symB,KB,mob,NEW,kbT_mo); 
    matt_times_vec(PLUS,symC,KC,mob,NEW,kcT_mo); 
    matt_times_vec(PLUS,symD,KD,mob,NEW,kdT_mo); 

    /* kx_kx_mo vectors */
    matn_times_vec(PLUS,symB,KB,kc_mo,NEW,kb_kc_mo);
    matn_times_vec(PLUS,symC,KC,kb_mo,ADD,kb_kc_mo);

    matn_times_vec(PLUS,symB,KB,kd_mo,NEW,kb_kd_mo);
    matn_times_vec(PLUS,symD,KD,kb_mo,ADD,kb_kd_mo);

    matn_times_vec(PLUS,symC,KC,kd_mo,NEW,kc_kd_mo);
    matn_times_vec(PLUS,symD,KD,kc_mo,ADD,kc_kd_mo);
    
    /* kxT_kxT_mo vectors */
    matt_times_vec(PLUS,symB,KB,kcT_mo,NEW,kbT_kcT_mo);
    matt_times_vec(PLUS,symC,KC,kbT_mo,ADD,kbT_kcT_mo);

    matt_times_vec(PLUS,symB,KB,kdT_mo,NEW,kbT_kdT_mo);
    matt_times_vec(PLUS,symD,KD,kbT_mo,ADD,kbT_kdT_mo);

    matt_times_vec(PLUS,symC,KC,kdT_mo,NEW,kcT_kdT_mo);
    matt_times_vec(PLUS,symD,KD,kcT_mo,ADD,kcT_kdT_mo);
    
    /* kx_kx_kx_mo vectors */
    matn_times_vec(PLUS,symB,KB,kc_kd_mo,NEW,kb_kc_kd_mo);
    matn_times_vec(PLUS,symC,KC,kb_kd_mo,ADD,kb_kc_kd_mo);
    matn_times_vec(PLUS,symD,KD,kb_kc_mo,ADD,kb_kc_kd_mo);
    
    /* kxT_kxT_kxT_mo vectors */
    matt_times_vec(PLUS,symB,KB,kcT_kdT_mo,NEW,kbT_kcT_kdT_mo);
    matt_times_vec(PLUS,symC,KC,kbT_kdT_mo,ADD,kbT_kcT_kdT_mo);
    matt_times_vec(PLUS,symD,KD,kbT_kcT_mo,ADD,kbT_kcT_kdT_mo);

     /* That's it! */

/* COMMENTS ON COMMUTATOR ROUTINES */
/*
   We need to construct sums of commutators such as
    [kappaB,[kappaC,omega(x)]] + [kappaC,[kappaB,omega(x)]],
    [kappaC,[kappaD,omega(x)]] + [kappaD,[kappaC,omega(x)]],
    [kappaB,[kappaD,omega(x)]] + [kappaD,[kappaB,omega(x)]],
    and sum of  [kappaB,[kappaC,[kappaD,omega]]] permutations,
  
    where kappaB,kappaC,kappaD and omega(x) are square matrices of size
    inforb_.norbt*inforb_.norbt.

    The matrix (Matlab) notation: a*b is always matrix matrix multiplication, a' is
    transpose of a.

    Let's call omega matirx O and use matrix notation. The omega matrix has
    special structure: 
                   O=mo*mo' and O(x)=mo*mog(x)' + mog(x)*mo'
    where mo is column vector.

    The commutator [KB,O] is then expressed as:
                   [KB,O] = KB*O - O*KB
    using O=mo*mo' we get
                   [KB,O] = KB*O - O*KB
                          = KB*mo*mo' - mo*mo'*KB
                          = KB*mo*mo' - mo*(KB'*mo)'
                          = kb_mo*mo' - mo*kbT_mo'

                    kb_mo = KB*mo
                    kbT_mo= KB'*mo
    
    In this manner we find all single, double and triple commutators:
    
    // commutator [KB,O] 
    Ko1 = kb_mo*mo' - mo*kbT_mo';
    // commutator [KC,O] 
    Ko1 = kc_mo*mo' - mo*kcT_mo';
    // commutator [KD,O] 
    Ko1 = kd_mo*mo' - mo*kdT_mo';

    // commutator [KC,[KB,O]] 
    KKo1 = kc_kb_mo*mo' - kc_mo*kbT_mo' - kb_mo*kcT_mo' + mo*kcT_kbT_mo'
    // commutator [KB,[KC,O]] 
    KKo1 = kb_kc_mo*mo' - kb_mo*kcT_mo' - kc_mo*kbT_mo' + mo*kbT_kcT_mo'


    // commutator [KB,[KD,O]] 
    KKo1 = kb_kd_mo*mo' - kb_mo*kdT_mo' - kd_mo*kbT_mo' + mo*kbT_kdT_mo'
    // commutator [KD,[KB,O]] 
    KKo1 = kd_kb_mo*mo' - kd_mo*kbT_mo' - kb_mo*kdT_mo' + mo*kdT_kbT_mo'


    // commutator [KC,[KD,O]] 
    KKo1 = kc_kd_mo*mo' - kc_mo*kdT_mo' - kd_mo*kcT_mo' + mo*kcT_kdT_mo'
    // commutator [KD,[KC,O]] 
    KKo1 = kd_kc_mo*mo' - kd_mo*kcT_mo' - kc_mo*kdT_mo' + mo*kdT_kcT_mo'

    // commutator [KD,[KC,[KB,O]]] 

            24               25                26                  27
    KKKo = kd_kc_kb_mo*mo' - kd_kc_mo*kbT_mo' - kd_kb_mo*kcT_mo' + kd_mo*kcT_kbT_mo' - ...
         - kc_kb_mo*kdT_mo' + kc_mo*kdT_kbT_mo' + kb_mo*kdT_kcT_mo' - mo*kdT_kcT_kbT_mo'
              0                1      CB           2                      3            

    // commutator [KC,[KD,[KB,O]]] 

           36                 37                 38                 39
    KKKo = kc_kd_kb_mo*mo' - kc_kd_mo*kbT_mo' - kc_kb_mo*kdT_mo' + kc_mo*kdT_kbT_mo' - ...
         - kd_kb_mo*kcT_mo' + kd_mo*kcT_kbT_mo' + kb_mo*kcT_kdT_mo' - mo*kcT_kdT_kbT_mo'
             12                13     DB            14                  15              

    // commutator [KD,[KB,[KC,O]]] 

           28                 29                 30                  31
    KKKo = kd_kb_kc_mo*mo' - kd_kb_mo*kcT_mo' - kd_kc_mo*kbT_mo' + kd_mo*kbT_kcT_mo' - ...
         - kb_kc_mo*kdT_mo' + kb_mo*kdT_kcT_mo' + kc_mo*kdT_kbT_mo' - mo*kdT_kbT_kcT_mo'
             4                  5     BC            6                   7               

    // commutator [KB,[KD,[KC,O]]] 

            44                45                 46                  47
    KKKo = kb_kd_kc_mo*mo' - kb_kd_mo*kcT_mo' - kb_kc_mo*kdT_mo' + kb_mo*kdT_kcT_mo' - ...
         - kd_kc_mo*kbT_mo' + kd_mo*kbT_kcT_mo' + kc_mo*kbT_kdT_mo' - mo*kbT_kdT_kcT_mo'
              20                21    DC           22                  23               

    // commutator [KC,[KB,[KD,O]]] 

             32                33                 34                  35
    KKKo = kc_kb_kd_mo*mo' - kc_kb_mo*kdT_mo' - kc_kd_mo*kbT_mo' + kc_mo*kbT_kdT_mo' - ...
         - kb_kd_mo*kcT_mo' + kb_mo*kcT_kdT_mo' + kd_mo*kcT_kbT_mo' - mo*kcT_kbT_kdT_mo'
               8                  9     BD           10                  11             

    // commutator [KB,[KC,[KD,O]]] 

             40                 41                42                 43
    KKKo = kb_kc_kd_mo*mo' - kb_kc_mo*kdT_mo' - kb_kd_mo*kcT_mo' + kb_mo*kcT_kdT_mo' - ...
         - kc_kd_mo*kbT_mo' + kc_mo*kbT_kdT_mo' + kd_mo*kbT_kcT_mo' - mo*kbT_kcT_kdT_mo'
                16                 17   CD           18                  19             

After simplification
BCD = - mo*(kbT_kcT_kdT_mo + kbT_kdT_kcT_mo + kcT_kbT_kdT_mo + kcT_kdT_kbT_mo + kdT_kbT_kcT_mo + kdT_kcT_kbT_mo)
      + 3*kb_mo*(kcT_kdT_mo + kdT_kcT_mo)
      + 3*kc_mo*(kbT_kdT_mo + kdT_kbT_mo)
      + 3*kd_mo*(kbT_kcT_mo + kcT_kbT_mo)
      - 3*(kb_kc_mo + kc_kb_mo)*kdT_mo
      - 3*(kc_kd_mo + kd_kc_mo)*kbT_mo
      - 3*(kb_kd_mo + kd_kb_mo)*kcT_mo
      + (kb_kc_kd_mo + kb_kd_kc_mo + kc_kb_kd_mo + kc_kd_kb_mo + kd_kb_kc_mo + kd_kc_kb_mo)*mo

BC = - 2*kc_mo*kbT_mo
     + mo*(kbT_kcT_mo + kcT_kbT_mo)
     - 2*kb_mo*kcT_mo
     + (kb_kc_mo + kc_kb_mo)*mo

CD = - 2*kd_mo*kcT_mo
     + mo*(kcT_kdT_mo + kdT_kcT_mo)
     - 2*kc_mo*kdT_mo
     + (kc_kd_mo + kd_kc_mo)*mo

BD = - 2*kd_mo*kbT_mo
     + mo*(kbT_kdT_mo + kdT_kbT_mo)
     - 2*kb_mo*kdT_mo
     + (kb_kd_mo + kd_kb_mo)*mo

     
*/  
} /* End of eval_omega_matrices() */

static void
assign_omega_vectors(
                real *moa, real *mob,
		real *vecA, real *vecB,
                CommData *omB, CommData *omC, CommData *omD,
                CommData *omBC, CommData *omCD, CommData *omBD,
                CommData *omBCD)
{
    real *kb_mo, *kc_mo, *kd_mo;  
    real *kbT_mo, *kcT_mo, *kdT_mo;

    real *kb_kc_mo, *kc_kd_mo, *kb_kd_mo;
    real *kbT_kcT_mo, *kcT_kdT_mo, *kbT_kdT_mo;
 
    real *kb_kc_kd_mo;

    real *kbT_kcT_kdT_mo;

    /* assign vectors to tmpA memory */
    kb_mo           = &(vecA[0*inforb_.norbt]);
    kc_mo           = &(vecA[1*inforb_.norbt]);
    kd_mo           = &(vecA[2*inforb_.norbt]);
  
    kb_kc_mo        = &(vecA[3*inforb_.norbt]); 
    kc_kd_mo        = &(vecA[4*inforb_.norbt]); 
    kb_kd_mo        = &(vecA[5*inforb_.norbt]); 

    kb_kc_kd_mo     = &(vecA[6*inforb_.norbt]);


    /* assign vectors to vecB memory (same way as vecA) */
    kbT_mo           = &(vecB[ 0*inforb_.norbt]);
    kcT_mo           = &(vecB[ 1*inforb_.norbt]);
    kdT_mo           = &(vecB[ 2*inforb_.norbt]);

    kbT_kcT_mo       = &(vecB[ 3*inforb_.norbt]); 
    kcT_kdT_mo       = &(vecB[ 4*inforb_.norbt]);
    kbT_kdT_mo       = &(vecB[ 5*inforb_.norbt]);

    kbT_kcT_kdT_mo   = &(vecB[ 6*inforb_.norbt]);

    /* Form A and B vector maps  */
    omB->A[0]= kb_mo; omB->B[0]= mob;
    omB->A[1]= moa  ; omB->B[1]= kbT_mo;          // -1
    
    omC->A[0]= kc_mo; omC->B[0]= mob;
    omC->A[1]= moa  ; omC->B[1]= kcT_mo;          // -1

    omD->A[0]= kd_mo; omD->B[0]= mob;
    omD->A[1]= moa  ; omD->B[1]= kdT_mo;          // -1

    omBC->A[0]= kc_mo   ; omBC->B[0]= kbT_mo;     /* coeAB =-2 */
    omBC->A[1]= moa     ; omBC->B[1]= kbT_kcT_mo; 
    omBC->A[2]= kb_mo   ; omBC->B[2]= kcT_mo;     /* coeAB =-2 */
    omBC->A[3]= kb_kc_mo; omBC->B[3]= mob;

    omCD->A[0]= kd_mo   ; omCD->B[0]= kcT_mo;     /* coeAB =-2 */
    omCD->A[1]= moa     ; omCD->B[1]= kcT_kdT_mo;
    omCD->A[2]= kc_mo   ; omCD->B[2]= kdT_mo;     /* coeAB =-2 */
    omCD->A[3]= kc_kd_mo; omCD->B[3]= mob;

    omBD->A[0]= kd_mo   ; omBD->B[0]= kbT_mo;     /* coeAB =-2 */
    omBD->A[1]= moa     ; omBD->B[1]= kbT_kdT_mo;
    omBD->A[2]= kb_mo   ; omBD->B[2]= kdT_mo;     /* coeAB =-2 */
    omBD->A[3]= kb_kd_mo; omBD->B[3]= mob;

    omBCD->A[0]= moa        ; omBCD->B[0]= kbT_kcT_kdT_mo; // -1
    omBCD->A[1]= kb_mo      ; omBCD->B[1]= kcT_kdT_mo; /* coeAB =3 */
    omBCD->A[2]= kc_mo      ; omBCD->B[2]= kbT_kdT_mo; /* coeAB =3 */
    omBCD->A[3]= kd_mo      ; omBCD->B[3]= kbT_kcT_mo; /* coeAB =3 */
    omBCD->A[4]= kb_kc_mo   ; omBCD->B[4]= kdT_mo;         // -3          
    omBCD->A[5]= kc_kd_mo   ; omBCD->B[5]= kbT_mo;         // -3            
    omBCD->A[6]= kb_kd_mo   ; omBCD->B[6]= kcT_mo;         // -3          
    omBCD->A[7]= kb_kc_kd_mo; omBCD->B[7]= mob;            // +1

} /* end assign_omega_vectors */

static void
add_commutators_full(CommData **CommList, int maxlst, real *dftcontr) {
     register int vec;
     int i,lst;
     real factor;

     int vec_num;
     real **mapA; 
     real **mapB;
     real *coef;
     real pref;
     real fac;
     real part_fac;
 
     for (i=0;i<inforb_.norbt; i++ ) {
      for (lst=0; lst<=maxlst; lst++)  {
       mapA = CommList[lst]->A;
       mapB = CommList[lst]->B;
       pref = CommList[lst]->pref;
       vec_num = CommList[lst]->vec_num;
       coef = CommList[lst]->coeAB;
       fac = CommList[lst]->fac;
     
       if(pref*pref < 1e-32) continue; 
    
       part_fac = fac*pref;
 
       for (vec=0; vec<vec_num; vec++) {
          factor = coef[vec]*mapB[vec][i]*part_fac;
          if (factor*factor > 1e-32)
          daxpy_(&inforb_.norbt, &factor, mapA[vec], &ONEI, &dftcontr[i*inforb_.norbt], &ONEI);
       } /* end vec */
     }   /* end lst */
    }    /* end i   */
}

static void
add_commutators_ij(CommData **CommList, int maxlst, real *dftcontr) {
     register int vec;
     int i,j,lst;
     real factor;
                                                                                                                                       
     int vec_num;
     real **mapA;
     real **mapB;
     real *coef;
     real pref;
     real fac;
     real part_fac;
                                                                                                                                       
     for (i=0;i<inforb_.norbt; i++ ) {
     for (j=0;j<inforb_.norbt; j++ ) {
      for (lst=0; lst<=maxlst; lst++)  {
       mapA = CommList[lst]->A;
       mapB = CommList[lst]->B;
       pref = CommList[lst]->pref;
       vec_num = CommList[lst]->vec_num;
       coef = CommList[lst]->coeAB;
       fac = CommList[lst]->fac;
                                                                                                                                       
       part_fac = fac*pref;
                                                                                                                                       
       for (vec=0; vec<vec_num; vec++) {
          factor = coef[vec]*mapB[vec][i]*part_fac;
          dftcontr[i*inforb_.norbt + j] += factor*mapA[vec][j];
       } /* end vec */
     }   /* end lst */
    }    /* end j   */
    }    /* end i   */
}

#if 0
static void
add_commutators_symm(CommData **commlist, int maxlst, real *dftcontr) {
     register int ii;
     int i, isym, vec;
     int vec_num = d->vec_num;
     real **mapA = d->A;
     real **mapB = d->B;
     int  sym = d->sym;
     real factor;

     int iorbs, noccs, nvirs, nuorb;
     int iorbi, nocci, nviri, nuorbi;

     for (vec=0; vec<vec_num; vec++) {
       for (isym=0; isym<inforb_.nsym; isym++) {
         iorbs = inforb_.iorb[isym];
         noccs = inforb_.nocc[isym];
         nvirs = inforb_.nvir[isym]; 
         nuorb  = noccs+nvirs;
         i     = inforb_.muld2h[sym][isym]-1;
         iorbi = inforb_.iorb[i];
         nocci = inforb_.nocc[i];
         nviri = inforb_.nvir[i];
         nuorbi = nocci + nviri;
         
         for(ii=iorbs; ii<iorbs+nuorb; ii++) {
           factor = mapB[vec][ii]*d->fac*pref;
           daxpy_(&nuorbi, &factor, &mapA[vec][iorbi], &ONEI, &dftcontr[ii*inforb_.norbt + iorbi], &ONEI);
         }
       }
     } /* end of vec */


} /* end of add_commutator_symm */
#endif


/* eval_rho_vars: rho-dependent temporary vars */
static void
eval_rho_vars(DftGrid *grid, CubeFastData* d)
{
    /* Set Omega */
    d->omega.A[0] = d->omega.B[0] = grid->mov;

    /* Evaluate unscaled omega commutators  vectors */
    eval_omega_vectors(
         d, d->kappaB, d->kappaC, d->kappaD,
         grid->mov, grid->mov,
         d->vecA, d->vecB);

    assign_omega_vectors(
         grid->mov, grid->mov,
         d->vecA, d->vecB,
         &d->omB, &d->omC, &d->omD,
         &d->omBC, &d->omCD, &d->omBD,
         &d->omBCD
         );

    
    /* Properly scaled expectation values (roB, roC, roD) */
    d->roB   = inactive_trace_AB(&d->omB);
    d->roC   = inactive_trace_AB(&d->omC);
    d->roD   = inactive_trace_AB(&d->omD);

    d->roBC  = inactive_trace_AB(&d->omBC);
    d->roCD  = inactive_trace_AB(&d->omCD);
    d->roBD  = inactive_trace_AB(&d->omBD);

    d->roBCD = inactive_trace_AB(&d->omBCD);

}

static void
eval_grad_vars(DftGrid *grid, CubeFastData* d)
{
    /* temporary variables */
    int x;

    /* set aux variables to zero */
    d->sB = d->sC =  d->sD = d->sBC = d->sCD = d->sBD = d->sBCD =
    d->pBC = d->pCD = d->pBD = d->pBCD = 0.0;

    for(x=0; x<3; x++) {

    /* Set omegax */
    d->omegax[0][x].A[0] = &grid->mog[x*inforb_.norbt];
    d->omegax[0][x].B[0] = grid->mov;

    d->omegax[1][x].A[0] = grid->mov;
    d->omegax[1][x].B[0] = &grid->mog[x*inforb_.norbt];

    eval_omega_vectors(
         d, d->kappaB, d->kappaC, d->kappaD,
         &grid->mog[x*inforb_.norbt], &grid->mog[x*inforb_.norbt],
         d->vecAx[x], d->vecBx[x]
         );

    assign_omega_vectors(
         grid->mov, &grid->mog[x*inforb_.norbt],
         d->vecA, d->vecBx[x],
         &d->omBx[0][x], &d->omCx[0][x], &d->omDx[0][x],
         &d->omBCx[0][x], &d->omCDx[0][x], &d->omBDx[0][x],
         &d->omBCDx[0][x]
         );

    assign_omega_vectors(
         &grid->mog[x*inforb_.norbt], grid->mov,
         d->vecAx[x], d->vecB,
         &d->omBx[1][x], &d->omCx[1][x], &d->omDx[1][x],
         &d->omBCx[1][x], &d->omCDx[1][x], &d->omBDx[1][x],
         &d->omBCDx[1][x]
         );

    /* auxaliary variables for prefactors */
    /* traces first */
    d->rox[x] = inactive_trace_AB(&d->omegax[0][x]) + inactive_trace_AB(&d->omegax[1][x]);

    d->roBx[x] = inactive_trace_AB(&d->omBx[0][x]) + inactive_trace_AB(&d->omBx[1][x]);
    d->roCx[x] = inactive_trace_AB(&d->omCx[0][x]) + inactive_trace_AB(&d->omCx[1][x]);
    d->roDx[x] = inactive_trace_AB(&d->omDx[0][x]) + inactive_trace_AB(&d->omDx[1][x]);
                 
    d->roBCx[x] = inactive_trace_AB(&d->omBCx[0][x]) + inactive_trace_AB(&d->omBCx[1][x]);
    d->roCDx[x] = inactive_trace_AB(&d->omCDx[0][x]) + inactive_trace_AB(&d->omCDx[1][x]);
    d->roBDx[x] = inactive_trace_AB(&d->omBDx[0][x]) + inactive_trace_AB(&d->omBDx[1][x]);

    d->roBCDx[x] = inactive_trace_AB(&d->omBCDx[0][x]) + inactive_trace_AB(&d->omBCDx[1][x]);


    /* and now auxaliary variables sB, sBC, sBCD, pBC, pBCD */
    d->sB   += d->roBx[x]*d->rox[x];
    d->sC   += d->roCx[x]*d->rox[x];
    d->sD   += d->roDx[x]*d->rox[x];
    d->sBC  += d->roBCx[x]*d->rox[x];
    d->sCD  += d->roCDx[x]*d->rox[x];
    d->sBD  += d->roBDx[x]*d->rox[x];
    d->sBCD += d->roBCDx[x]*d->rox[x];

    d->pBC  += d->roBx[x]*d->roCx[x];
    d->pCD  += d->roCx[x]*d->roDx[x];
    d->pBD  += d->roBx[x]*d->roDx[x];

    d->pBCD += (d->roBx[x]*d->roCDx[x] +
               d->roCx[x]*d->roBDx[x] + 
               d->roDx[x]*d->roBCx[x]); 
    }

}

static void
add_dft_contribution(DftGrid* grid, CubeFastData* d)
{
    int x;
    real pref;
    real* dftcontr = d->dftcontr;
    void (*add_commutators)() = d->add_commutators;
    FourthDrv drvs; /* the functional derivatives */
    int lst = 0;
    CommData *commlist[49];
    

    /* for closed shell ngradab = ngrada = ngradb */
    dftpot3ab_(&drvs, &grid->curr_weight, &grid->dp,
               &grid->dogga);

#if 0
    /* distribute three-index transformed densities (prefactor r) */
    if(d->addfock) {
        pref = 0.5*drvs.fR;
        //pref = drvs.fR;
	//pref = 0.0;
        daxpy_(&norbt2, &pref, d->omBCD, &ONEI, dftcontr, &ONEI);
    }
    //return;
#endif

#ifdef TEST
    add_commutator_full(&d->omega, grid->curr_weight, d->tmp_f);
    add_commutator_symm(&d->omega, grid->curr_weight, d->tmp_s);
#endif

    /* distribute two-index transformed densities (prefactors rB, rC, rD */
    /* rB */
    pref = RR*d->roB*drvs.fRR + RZ*2*d->sB*drvs.fRZ;
    d->omCD.pref = pref; commlist[lst] = &d->omCD;
    //daxpy_(&norbt2, &pref, d->omCD, &ONEI, dftcontr, &ONEI);
    
    /* rC */
    pref = RR*d->roC*drvs.fRR + RZ*2*d->sC*drvs.fRZ;
    d->omBD.pref = pref; commlist[++lst] = &d->omBD;
    //daxpy_(&norbt2, &pref, d->omBD, &ONEI, dftcontr, &ONEI);

    /* rD */
    pref = RR*d->roD*drvs.fRR + RZ*2*d->sD*drvs.fRZ;
    d->omBC.pref = pref; commlist[++lst] = &d->omBC;
    //daxpy_(&norbt2, &pref, d->omBC, &ONEI, dftcontr, &ONEI);

    /* distribute one-index transformed densities (prefactors rBC, rBD, rCD */
    /* rBC */
    pref  = RRR*d->roB*d->roC*drvs.fRRR;
    pref += RRZ*2*(d->roB*d->sC + d->roC*d->sB)*drvs.fRRZ;
    pref += RZZ*4*d->sB*d->sC*drvs.fRZZ; 
    pref += RR*d->roBC*drvs.fRR; 
    pref += RZ*2*(d->sBC + d->pBC)*drvs.fRZ;
    d->omD.pref = pref; commlist[++lst] = &d->omD;
    //daxpy_(&norbt2, &pref, d->omD, &ONEI, dftcontr, &ONEI);

    /* rCD */
    pref  = RRR*d->roC*d->roD*drvs.fRRR;
    pref += RRZ*2*(d->roC*d->sD + d->roD*d->sC)*drvs.fRRZ;
    pref += RZZ*4*d->sC*d->sD*drvs.fRZZ;
    pref += RR*d->roCD*drvs.fRR;
    pref += RZ*2*(d->sCD + d->pCD)*drvs.fRZ;
    d->omB.pref = pref; commlist[++lst] = &d->omB;
    //daxpy_(&norbt2, &pref, d->omB, &ONEI, dftcontr, &ONEI);
    
    /* rBD */
    pref  = RRR*d->roB*d->roD*drvs.fRRR;
    pref += RRZ*2*(d->roB*d->sD + d->roD*d->sB)*drvs.fRRZ;
    pref += RZZ*4*d->sB*d->sD*drvs.fRZZ;
    pref += RR*d->roBD*drvs.fRR;
    pref += RZ*2*(d->sBD + d->pBD)*drvs.fRZ;
    d->omC.pref = pref; commlist[++lst] = &d->omC;
    //daxpy_(&norbt2, &pref, d->omC, &ONEI, dftcontr, &ONEI);


    /*  distribute rho (prefactor rBCD) */
    pref  = RRRR*d->roB*d->roC*d->roD*drvs.fRRRR;
    pref += RRRZ*2*(d->roB*d->roC*d->sD + d->roB*d->roD*d->sC +
                    d->roC*d->roD*d->sB)*drvs.fRRRZ;
    pref += RRZZ*4*(d->roB*d->sC*d->sD + d->roC*d->sB*d->sD +
                    d->roD*d->sB*d->sC)*drvs.fRRZZ;
    pref += RZZZ*8*d->sB*d->sC*d->sD*drvs.fRZZZ;
    pref += RRR*(d->roB*d->roCD + d->roC*d->roBD + d->roD*d->roBC)*drvs.fRRR;
    pref += RRZ*2*(d->roB*d->sCD + d->roC*d->sBD + d->roD*d->sBC +
                   d->roB*d->pCD + d->roC*d->pBD + d->roD*d->pBC +
                   d->sB*d->roCD + d->sC*d->roBD + d->sD*d->roBC )*drvs.fRRZ;
    pref += RZZ*4*(d->sB*d->sCD + d->sC*d->sBD + d->sD*d->sBC +
                   d->sB*d->pCD + d->sC*d->pBD + d->sD*d->pBC)*drvs.fRZZ;
    pref += RR*d->roBCD*drvs.fRR;
    pref += RZ*2*(d->sBCD + d->pBCD)*drvs.fRZ;
    d->omega.pref = pref; commlist[++lst] = &d->omega;
    //daxpy_(&norbt2, &pref, d->omega, &ONEI, dftcontr, &ONEI);
    
    if(!grid->dogga) {add_commutators(commlist,lst,dftcontr); return;} /* happilly home! */

    /* now the same thing has to be done for the grho matrices...
     * "It's a long way to the top if you want to rock'n'roll!" */

    for (x=0; x<3; x++) {

    /* prefactors qBx, qBy, qBz */
    /* qBx */
    pref = (RZ*2*d->roB*drvs.fRZ + ZZ*4*d->sB*drvs.fZZ)*d->rox[x] +
             Z*2*d->roBx[x]*drvs.fZ;
    d->omCDx[0][x].pref = pref;  commlist[++lst] = &d->omCDx[0][x];
    d->omCDx[1][x].pref = pref;  commlist[++lst] = &d->omCDx[1][x];
    //daxpy_(&norbt2, &pref, &d->omCDx[x*inforb_.n2orbx], &ONEI, dftcontr, &ONEI);
	    
    /* prefactors qCx, qCy, qCz */
    /* qCx */
    pref = (RZ*2*d->roC*drvs.fRZ + ZZ*4*d->sC*drvs.fZZ)*d->rox[x] +
             Z*2*d->roCx[x]*drvs.fZ;
    d->omBDx[0][x].pref = pref;  commlist[++lst] = &d->omBDx[0][x];
    d->omBDx[1][x].pref = pref;  commlist[++lst] = &d->omBDx[1][x];
    //daxpy_(&norbt2, &pref, &d->omBDx[x*inforb_.n2orbx], &ONEI, dftcontr, &ONEI);
	    
	    
    /* prefactors qDx, qDy, qDz */
    /* qDx */
    pref = (RZ*2*d->roD*drvs.fRZ + ZZ*4*d->sD*drvs.fZZ)*d->rox[x] +
             Z*2*d->roDx[x]*drvs.fZ;
    d->omBCx[0][x].pref = pref;  commlist[++lst] = &d->omBCx[0][x];
    d->omBCx[1][x].pref = pref;  commlist[++lst] = &d->omBCx[1][x];
    //daxpy_(&norbt2, &pref, &d->omBCx[x*inforb_.n2orbx], &ONEI, dftcontr, &ONEI);
	    

    /* prefactors qBCx, qBCy, qBCz */
    /* qBCx */
    pref  = RRZ*2*d->roB*d->roC*drvs.fRRZ*d->rox[x];
    pref += RZZ*4*(d->roB*d->sC + d->roC*d->sB)*drvs.fRZZ*d->rox[x];
    pref += ZZZ*8*d->sB*d->sC*drvs.fZZZ*d->rox[x];
    pref += ZZ*4*(d->roBx[x]*d->sC + d->sB*d->roCx[x] + d->sBC*d->rox[x] + 
               d->pBC*d->rox[x])*drvs.fZZ;
    pref += RZ*2*(d->roB*d->roCx[x] + d->roC*d->roBx[x] +
	       d->roBC*d->rox[x])*drvs.fRZ;
    pref += Z*2*d->roBCx[x]*drvs.fZ;
    d->omDx[0][x].pref = pref;  commlist[++lst] = &d->omDx[0][x];
    d->omDx[1][x].pref = pref;  commlist[++lst] = &d->omDx[1][x];
    //daxpy_(&norbt2, &pref, &d->omDx[x*inforb_.n2orbx], &ONEI, dftcontr, &ONEI);


    /* prefactors qCDx, qCDy, qCDz */
    /* qCDx */
    pref  = RRZ*2*d->roC*d->roD*drvs.fRRZ*d->rox[x];
    pref += RZZ*4*(d->roC*d->sD + d->roD*d->sC)*drvs.fRZZ*d->rox[x];
    pref += ZZZ*8*d->sC*d->sD*drvs.fZZZ*d->rox[x];
    pref += ZZ*4*(d->roCx[x]*d->sD + d->sC*d->roDx[x] + d->sCD*d->rox[x] +
               d->pCD*d->rox[x])*drvs.fZZ;
    pref += RZ*2*(d->roC*d->roDx[x] + d->roD*d->roCx[x] +
	       d->roCD*d->rox[x])*drvs.fRZ;
    pref += Z*2*d->roCDx[x]*drvs.fZ;
    d->omBx[0][x].pref = pref;  commlist[++lst] = &d->omBx[0][x];
    d->omBx[1][x].pref = pref;  commlist[++lst] = &d->omBx[1][x];
    //daxpy_(&norbt2, &pref, &d->omBx[x*inforb_.n2orbx], &ONEI, dftcontr, &ONEI);

	    
    /* prefactors qBDx, qBDy, qBDz */
    /* qBDx */
    pref  = RRZ*2*d->roB*d->roD*drvs.fRRZ*d->rox[x];
    pref += RZZ*4*(d->roB*d->sD + d->roD*d->sB)*drvs.fRZZ*d->rox[x];
    pref += ZZZ*8*d->sB*d->sD*drvs.fZZZ*d->rox[x];
    pref += ZZ*4*(d->roBx[x]*d->sD + d->sB*d->roDx[x] + d->sBD*d->rox[x] +
               d->pBD*d->rox[x])*drvs.fZZ;
    pref += RZ*2*(d->roB*d->roDx[x] + d->roD*d->roBx[x] +
               d->roBD*d->rox[x])*drvs.fRZ;
    pref += Z*2*d->roBDx[x]*drvs.fZ;
    d->omCx[0][x].pref = pref;  commlist[++lst] = &d->omCx[0][x];
    d->omCx[1][x].pref = pref;  commlist[++lst] = &d->omCx[1][x];
    //daxpy_(&norbt2, &pref, &d->omCx[x*inforb_.n2orbx], &ONEI, dftcontr, &ONEI);

    /* Prefactors qBCDx, qBCDy, qBCDz */
    pref  =  ZZZZ*16*d->rox[x]*d->sB*d->sC*d->sD*drvs.fZZZZ;
    pref +=  RZZZ*8*d->rox[x]*(d->roB*d->sC*d->sD + d->roC*d->sB*d->sD +
               d->roD*d->sB*d->sC)*drvs.fRZZZ;
    pref +=  RRZZ*4*d->rox[x]*(d->roB*d->roC*d->sD + d->roB*d->roD*d->sC +
               d->roD*d->roC*d->sB)*drvs.fRRZZ;
    pref +=  RRRZ*2*d->rox[x]*d->roB*d->roC*d->roD*drvs.fRRRZ;
    pref +=  ZZZ*8*d->rox[x]*(d->sB*d->sCD + d->sC*d->sBD + d->sD*d->sBC +
               d->sB*d->pCD + d->sC*d->pBD + d->sD*d->pBC)*drvs.fZZZ;
    pref +=  ZZZ*8*(d->roBx[x]*d->sC*d->sD + d->roCx[x]*d->sB*d->sD +
               d->roDx[x]*d->sB*d->sC)*drvs.fZZZ;
    pref +=  RZZ*4*d->rox[x]*(d->roB*d->sCD + d->roC*d->sBD + d->roD*d->sBC +
               d->roB*d->pCD + d->roC*d->pBD + d->roD*d->pBC)*drvs.fRZZ;
    pref +=  RZZ*4*d->rox[x]*(d->sB*d->roCD + d->sC*d->roBD +
               d->sD*d->roBC)*drvs.fRZZ;
    pref +=  RZZ*4*(d->roBx[x]*(d->roC*d->sD + d->roD*d->sC) + 
               d->roCx[x]*(d->roD*d->sB + d->roB*d->sD) +
               d->roDx[x]*(d->roC*d->sB + d->roB*d->sC))*drvs.fRZZ;
    pref +=  RRZ*2*(d->roBx[x]*d->roC*d->roD + d->roCx[x]*d->roB*d->roD +
               d->roDx[x]*d->roB*d->roC)*drvs.fRRZ;
    pref +=  RRZ*2*d->rox[x]*(d->roB*d->roCD + d->roC*d->roBD +
               d->roD*d->roBC)*drvs.fRRZ;
    pref +=  ZZ*4*(d->roBx[x]*d->sCD + d->roCx[x]*d->sBD + d->roDx[x]*d->sBC +
               d->roBx[x]*d->pCD + d->roCx[x]*d->pBD + d->roDx[x]*d->pBC +
               d->rox[x]*d->sBCD + d->rox[x]*d->pBCD)*drvs.fZZ;
    pref +=  ZZ*4*(d->roBCx[x]*d->sD + d->roCDx[x]*d->sB +
               d->roBDx[x]*d->sC)*drvs.fZZ;
    pref +=  RZ*2*(d->roBx[x]*d->roCD + d->roCx[x]*d->roBD +
               d->roDx[x]*d->roBC + d->roB*d->roCDx[x] + d->roC*d->roBDx[x] +
               d->roD*d->roBCx[x] + d->rox[x]*d->roBCD)*drvs.fRZ;
    pref +=  Z*2*d->roBCDx[x]*drvs.fZ;

    d->omegax[0][x].pref = pref;  commlist[++lst] = &d->omegax[0][x];
    d->omegax[1][x].pref = pref;  commlist[++lst] = &d->omegax[1][x];
    //daxpy_(&norbt2, &pref, &d->omegax[x*inforb_.n2orbx], &ONEI, dftcontr, &ONEI);

   } /* end of x loop */

    add_commutators(commlist,lst,dftcontr);

} /* end of add_dft_contribution */


/* MAIN: */

void
fast_callback(DftGrid* grid, CubeFastData* data)
{
#ifdef TIMING
    real sec, tmpsec, dummy;
    void gettim_(real *a, real *b);
    gettim_(&sec,&dummy);
#endif

    eval_rho_vars(grid, data);
    if(grid->dogga) {
	eval_grad_vars(grid, data);     
    }

#ifdef TIMING
    gettim_(&tmpsec,&dummy);
    data->timedata[0] += tmpsec-sec;
    sec=tmpsec;
#endif

    add_dft_contribution(grid, data);

#ifdef TIMING
    gettim_(&tmpsec,&dummy);
    data->timedata[1] += tmpsec-sec;
    data->timedata[2] += 1.0;
#endif

}

#ifdef VAR_MPI
#include <mpi.h>
#include <our_extra_mpi.h>
#include "infpar.h"

void
dft_cr_resp_slave(real* work, integer* lwork, const integer* iprint)
{
    real* fi    = calloc(inforb_.n2orbx, sizeof(real));              /* OUT */
    real *cmo = dal_malloc(inforb_.norbt*inforb_.nbast*sizeof(real)); /*IN */
    real *kappaB= dal_malloc(inforb_.n2orbx*sizeof(real));            /*IN */
    real *kappaC= dal_malloc(inforb_.n2orbx*sizeof(real));            /*IN */
    real *kappaD= dal_malloc(inforb_.n2orbx*sizeof(real));            /*IN */
    integer  symB, symC, symD;                                        /*IN */
    dftcrcf_(fi,     cmo,
             kappaB, &symB,
             kappaC, &symC,
             kappaD, &symD,
             work,   lwork);
    free(kappaB);
    free(kappaC);
    free(kappaD);
    free(cmo);
    free(fi);
}


static void
dft_cr_resp_sync_slaves(real* cmo, real* kappaB, real* kappaC, real *kappaD,
                        integer* symB, integer* symC, integer* symD)
{
    static const SyncData sync_data[] = {
        { inforb_.nocc,   8, fortran_MPI_INT },
        { inforb_.nvir,   8, fortran_MPI_INT },
        { &inforb_.nocct, 1, fortran_MPI_INT },
        { &inforb_.nvirt, 1, fortran_MPI_INT },
    };
#ifdef C99_COMPILER
    const SyncData data2[] = {
        { cmo,     inforb_.norbt*inforb_.nbast,MPI_DOUBLE },
        { kappaB,  inforb_.n2orbx,             MPI_DOUBLE },
        { kappaC,  inforb_.n2orbx,             MPI_DOUBLE },
        { kappaD,  inforb_.n2orbx,             MPI_DOUBLE },
        { symB,    1,                          fortran_MPI_INT    },
        { symC,    1,                          fortran_MPI_INT    },
        { symD,    1,                          fortran_MPI_INT    },
    };
#else /* C99_COMPILER */
    /* this is more error-prone but some compilers (HP/UX)... */
    static SyncData data2[] =
    { {NULL, 0, MPI_DOUBLE}, {NULL, 0, MPI_DOUBLE}, {NULL, 0, MPI_DOUBLE},
      {NULL, 0, MPI_DOUBLE},
      {NULL, 0, fortran_MPI_INT   }, {NULL, 0, fortran_MPI_INT   }, {NULL, 0, fortran_MPI_INT   }};
    data2[0].data = cmo;     data2[0].count = inforb_.norbt*inforb_.nbast;
    data2[1].data = kappaB;  data2[1].count = inforb_.n2orbx;
    data2[2].data = kappaC;  data2[2].count = inforb_.n2orbx;
    data2[3].data = kappaD;  data2[3].count = inforb_.n2orbx;
    data2[4].data = symB;    data2[4].count = 1;
    data2[5].data = symC;    data2[5].count = 1;
    data2[6].data = symD;    data2[6].count = 1;
#endif /* C99_COMPILER */
                                                                                       
    mpi_sync_data(sync_data, ELEMENTS(sync_data));
    mpi_sync_data(data2,     ELEMENTS(data2));
}

static __inline__ void
dft_cr_resp_collect_info(real* fi, real*work, int lwork)
{
    int sz = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &sz);
    if(sz<=1) return;
    CHECK_WRKMEM(inforb_.n2orbx, lwork);
    dcopy_(&inforb_.n2orbx, fi, &ONEI, work, &ONEI);
    MPI_Reduce(work, fi, inforb_.n2orbx, MPI_DOUBLE, MPI_SUM,
               daltoninfpar_.master, MPI_COMM_WORLD);
}

#else /* VAR_MPI */
#define dft_cr_resp_sync_slaves(cmo,kappaB,kappaC,kappaD,symB,symC,symD)
#define dft_cr_resp_collect_info(fi,work,lwork)
#endif /* VAR_MPI */

/* this is the routine for computing the DFT exchange-correlation 
   contribution to cubic response.
*/
void
FSYM(dftcrcf)(real* fi, real* cmo,
	      real* kappaB, integer* symB,
	      real* kappaC, integer* symC,
	      real* kappaD, integer* symD, 
	      real* work, integer* lwork)
{
    static int msg_printed = 0;
    integer norbt2 = inforb_.norbt*inforb_.norbt;
    void gettim_(real *a, real *b);
    real sec, tmpsec, dummy;
    DftCallbackData cbdata[1];
    CubeFastData* data;
    void dump_mat(char *name, real* mat, int dimm, int dimn);

    /* WARNING: NO work MAY BE done before syncing slaves! */
    dft_wake_slaves((DFTPropEvalMaster)dftcrcf_);       /* NO-OP in serial */
    dft_cr_resp_sync_slaves(cmo,kappaB,kappaC,kappaD,
                                symB,  symC,  symD );

    gettim_(&sec,&dummy);

    data = cubefast_data_new(kappaB, kappaC, kappaD,
                              *symB,  *symC,  *symD );
    

    if(!msg_printed) {
      fort_print("\n DFT-CR uses %s matrix multipication code.",
                 data->mm_code_version);
      msg_printed = 1;
    }
    
    cbdata[0].callback = (DftCallback)fast_callback;
    cbdata[0].cb_data = data;

    dft_integrate(cmo, work, lwork, cbdata, ELEMENTS(cbdata));

    dft_cr_resp_collect_info(data->dftcontr,work,*lwork); /* NO-OP in serial */
    daxpy_(&norbt2, &ONER, data->dftcontr, &ONEI, fi, &ONEI);

    cubefast_data_free(data);
 
    gettim_(&tmpsec,&dummy);
    fort_print("DFT-CR integration time: %10.2f s",tmpsec-sec);

}

#ifdef TEST
void dump_mat(char *name, real* mat, int dimm, int dimn)
{
    printf("%s\n",name); 
    int i, j;
    for(j=0; j<dimm; j++) {
        for(i=0; i<dimn; i++)
            printf("%5.3f ", mat[i*dimm + j]);
        puts("");
    }
        puts("");
}
#endif
