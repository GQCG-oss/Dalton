/*


!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!

!

*/
/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* quad-fast.c:

   (c) Pawel Salek, pawsa@theochem.kth.se, feb 2002

   It duplicates the functionality of quad-strict except that the
   matrix-matrix operations are replaced with matrix-vector.  

   NOTE: 
   The code can add the fock-like contribution [k,[k,omega]] if needed.

   For example QFOCK returns complete FOCK/Kohn-Sham matrix, including
   DFT contribution, as opposed to RSPFXD/RSPFX which return
   electronic part only (w/o DFT contribution). If this is the case
   (i.e. .NOT.DIRFCK) we add the DFT contribution to KS matrix later
   in DFTQRC. field addfock determines, if the fock-type contribution
   should be added or not.

   We allocate single temporary memory block to be able to use
   GEMV methods to generate distribution vectors.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>

#if defined(VAR_MPI)
#include <mpi.h>
#include <our_extra_mpi.h>
#endif

#ifdef FAST_TEST
#include <math.h>
#endif

#define __CVERSION__

#include "integrator.h"
#include "inforb.h"

enum {
    OFF_MOV = 0,
    OFF_YN  = 1,
    OFF_ZN  = 2,
    OFF_YNZN= 3,
    OFF_YT  = 4,
    OFF_ZT  = 5,
    OFF_YTZT= 6
};
/* Singleton objects */
struct QuadFastData_ {
    /* pointers to external data (data is not owned) */
    integer symY, symZ;
    const real* kappaY, *kappaZ; 
    integer ispinY, ispinZ, addfock;

    /* temporary variables */
    real* dftcontr;
    real pref2sum; /* this is used for checksumming, debugging etc. */
    real grada[3];
    real yy, zz, yzzy;
    real trY[3], trZ[3]; /* trY[x] = [kappa, grho_x])*/
    real trYsum, trZsum; /* = sum trY[x]*grad_x/|grho| */
    real trYZZYsum; 
    real trYZ[3], trZY[3]; /* trYZ[x] = [k_Z,[k_Y,grho_x]])*grad_x/|grho| */
    real trYZZY[3];  /* trYZZY=<Tr([kappa,[kappa, grho_x]])*grad_x/|grho| */
    real trYtimesZ; /* = sum [kY, grho_x]*[kappaZ, grho_x]*(grho_x/|grho|)² */

    real *memblock;
    real *mov;
    real *u_YN, *u_YT, *u_ZN, *u_ZT;
    real *v_YN_ZNp;  /* == v_YN_ZN + v_ZN_YN */
    real *v_YT_ZTp;  /* == v_YT_ZT + v_ZT_YT */

    /* GGA vectors */
    real *g_YN, *g_YT, *g_ZN, *g_ZT;
    real *x_YN_ZNp;  /* == v_YN_ZN + v_ZN_YN */
    real *x_YT_ZTp;  /* == v_YT_ZT + v_ZT_YT */
};
typedef struct QuadFastData_ QuadFastData;

/* struct QuadFastData:

   structure holding temporary variables for the computing the
   quadratic response contribution at given point. 

   The data inside should be only accessed only from the creation and
   release routines, and from the routine computing the contribution
   itself.
*/
static QuadFastData*
quadfast_data_new(const real* kY, integer symY, integer spinY, 
                  const real *kZ, integer symZ, integer spinZ, integer addfock)
{
    QuadFastData* res = malloc(sizeof(QuadFastData));

    res->kappaY = kY; res->symY = symY-1; res->ispinY = spinY;
    res->kappaZ = kZ; res->symZ = symZ-1; res->ispinZ = spinZ;
    res->addfock = addfock;
    

    res->dftcontr = alloc_mat_MO(1);
    res->pref2sum = 0;
    res->trYsum = res->trZsum = 0; /* initialize them for dogga == FALSE */
    res->trYZZYsum = res->trYtimesZ = 0;

    res->memblock = calloc(inforb_.norbt*7, sizeof(real));
    /* aliases below */
    res->mov      = res->memblock+inforb_.norbt*OFF_MOV;
    res->u_YN     = res->memblock+inforb_.norbt*OFF_YN;
    res->u_ZN     = res->memblock+inforb_.norbt*OFF_ZN;
    res->v_YN_ZNp = res->memblock+inforb_.norbt*OFF_YNZN;
    res->u_YT     = res->memblock+inforb_.norbt*OFF_YT;
    res->u_ZT     = res->memblock+inforb_.norbt*OFF_ZT;
    res->v_YT_ZTp = res->memblock+inforb_.norbt*OFF_YTZT;

    /* GGA data */
    res->g_YN = calloc(3*inforb_.norbt, sizeof(real));
    res->g_YT = calloc(3*inforb_.norbt, sizeof(real));
    res->g_ZN = calloc(3*inforb_.norbt, sizeof(real));
    res->g_ZT = calloc(3*inforb_.norbt, sizeof(real));
    res->x_YN_ZNp = calloc(3*inforb_.norbt, sizeof(real));
    res->x_YT_ZTp = calloc(3*inforb_.norbt, sizeof(real));

    return res;
}

static void 
quadfast_data_free(QuadFastData* tmp)
{
    free(tmp->dftcontr);

    free(tmp->memblock);

    free(tmp->g_YN);
    free(tmp->g_YT);
    free(tmp->g_ZN);
    free(tmp->g_ZT);
    free(tmp->x_YN_ZNp);
    free(tmp->x_YT_ZTp);
    free(tmp);
}

#if 0
// USE_SIMPLE_BUT_USUALLY_SLOWER_MM_CODE
/* this version of the routines multiplies vectors by a vaste number
   of zeros.  But it is simpler, easier to debug AND faster for
   systems with, say, less than 30 basis functions so it should
   definetely stay here just in case.
*/
static const char* mm_code_version = "full (slower)";
static __inline__ void
matn_times_vec_MO(int sym, const real* mat, const real* vec, real a, real* res)
{
    dgemv_("N", &inforb_.norbt, &inforb_.norbt, &ONER, mat,
	   &inforb_.norbt, vec, &ONEI, &a, res, &ONEI);
}
static __inline__ void
matt_times_vec_MO(int sym, const real* mat, const real* vec, real a, real* res)
{
    dgemv_("T", &inforb_.norbt, &inforb_.norbt, &ONER, mat,
	   &inforb_.norbt, vec, &ONEI, &a, res, &ONEI);
}
static __inline__ void
matn_times_vec3_MO(int sym, const real* mat, const real* vec3, real a,
                   real* res3)
{
    dgemm_("N", "N", &inforb_.norbt, &THREEI, &inforb_.norbt, &ONER, mat,
	   &inforb_.norbt, vec3, &inforb_.norbt, &a, res3, &inforb_.norbt);
}
static __inline__ void
matt_times_vec3_MO(int sym, const real* mat, const real* vec3, real a,
                   real* res3)
{
    dgemm_("T", "N", &inforb_.norbt, &THREEI, &inforb_.norbt, &ONER, mat,
	   &inforb_.norbt, vec3, &inforb_.norbt, &a, res3, &inforb_.norbt);
}
#else
/* mat{n,t}_times_vec_MO: routines for matrix-vector multiplication.
   We use the fact that kappa matrices have plenty of zeros for blocks
   corresponding to transformation between different symmetries not coupled
   by the operator.
   This code is a big win for large symmetric molecules.
   For CO2 with 90 orbitals (33,27,17,13) time per iteration goes
   down from 23.65s to 17.42s.
*/
static const char* mm_code_version = "symmetry adapted";
static __inline__ void
matn_times_vec_MO(int ops,const real* mat, const real* vec, real a, real* res)
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
		dgemv_("N", &nviri,&noccs, &ONER, &mat[begll],&inforb_.norbt, 
		       &vec[iorbs], &ONEI, &a, &res[iorbi+nocci], &ONEI);
	    else if(a==0) dzero_(&res[iorbi+nocci], &nviri);
	}
	if(nocci>0) {
	    if(nvirs>0)
		dgemv_("N", &nocci,&nvirs, &ONER, &mat[begur],&inforb_.norbt, 
		       &vec[iorbs+noccs], &ONEI, &a, &res[iorbi], &ONEI);
	    else if(a==0) dzero_(&res[iorbi], &nocci);
	}
    }
}

static __inline__ void
matt_times_vec_MO(int ops,const real* mat, const real* vec, real a, real* res)
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
                dgemv_("T", &nviri,&noccs, &ONER, &mat[begll],&inforb_.norbt,
                       &vec[iorbi+nocci], &ONEI, &a, &res[iorbs], &ONEI);
	    else if(a==0) dzero_(&res[iorbs], &noccs);
	}
        if(nvirs>0) {
            if(nocci>0)
                dgemv_("T", &nocci,&nvirs, &ONER, &mat[begur],&inforb_.norbt, 
                       &vec[iorbi], &ONEI, &a, &res[iorbs+noccs], &ONEI);
	    else if(a==0) dzero_(&res[iorbs+noccs], &nvirs);
	}
    }
}

static __inline__ void
matn_times_vec3_MO(int ops, const real* mat, const real* vec3, real a, 
                   real* res3)
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
                dgemm_("N","N", &nviri, &THREEI, &noccs, &ONER, 
                       &mat[begll],&inforb_.norbt, &vec3[iorbs],&inforb_.norbt,
                       &a, &res3[iorbi+nocci], &inforb_.norbt);
	    else if(a==0) { 
                dzero_(&res3[iorbi+nocci], &nviri);
                dzero_(&res3[iorbi+nocci+inforb_.norbt], &nviri);
                dzero_(&res3[iorbi+nocci+inforb_.norbt*2], &nviri);
            }
	}
	if(nocci>0) {
	    if(nvirs>0)
                dgemm_("N","N", &nocci, &THREEI, &nvirs, &ONER, 
                       &mat[begur],&inforb_.norbt, 
                       &vec3[iorbs+noccs],&inforb_.norbt,
                       &a, &res3[iorbi], &inforb_.norbt);
	    else if(a==0) {
                dzero_(&res3[iorbi], &nocci);
                dzero_(&res3[iorbi+inforb_.norbt], &nocci);
                dzero_(&res3[iorbi+inforb_.norbt*2], &nocci);
            }
	}
    }

}

static __inline__ void
matt_times_vec3_MO(int ops, const real* mat, const real* vec3, real a,
                   real* res3)
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
                dgemm_("T", "N", &noccs, &THREEI, &nviri, &ONER, 
                       &mat[begll], &inforb_.norbt, 
                       &vec3[iorbi+nocci],&inforb_.norbt,
                       &a, &res3[iorbs], &inforb_.norbt);
	    else if(a==0) {
                dzero_(&res3[iorbs], &noccs);
                dzero_(&res3[iorbs+inforb_.norbt], &noccs);
                dzero_(&res3[iorbs+inforb_.norbt*2], &noccs);
            }
	}
        if(nvirs>0) {
            if(nocci>0)
                dgemm_("T", "N", &nvirs, &THREEI, &nocci, &ONER, 
                       &mat[begur], &inforb_.norbt, 
                       &vec3[iorbi], &inforb_.norbt, 
                       &a, &res3[iorbs+noccs], &inforb_.norbt);
	    else if(a==0) {
                dzero_(&res3[iorbs+noccs], &nvirs);
                dzero_(&res3[iorbs+noccs+inforb_.norbt], &nvirs);
                dzero_(&res3[iorbs+noccs+inforb_.norbt*2], &nvirs);
            }
	}
    }
}
#endif

static __inline__ void
eval_rho_vars(DftGrid* grid, QuadFastData* d)
{
    int isym,iorbi;
    integer nocci;

    /* u_YN = kappaY*mov; u_YT = kappaY'*mov */
    matn_times_vec_MO(d->symY, d->kappaY, grid->mov, 0.0, d->u_YN);
    matt_times_vec_MO(d->symY, d->kappaY, grid->mov, 0.0, d->u_YT);
    matn_times_vec_MO(d->symZ, d->kappaZ, grid->mov, 0.0, d->u_ZN);
    matt_times_vec_MO(d->symZ, d->kappaZ, grid->mov, 0.0, d->u_ZT);

    /* v_YN_ZN = kappaY*kappaZ*mov; v_YT_ZT = kappaY'*(kappaZ'*mov) */
    matn_times_vec_MO(d->symY, d->kappaY, d->u_ZN, 0.0, d->v_YN_ZNp);
    matt_times_vec_MO(d->symY, d->kappaY, d->u_ZT, 0.0, d->v_YT_ZTp);
    matn_times_vec_MO(d->symZ, d->kappaZ, d->u_YN, 1.0, d->v_YN_ZNp);
    matt_times_vec_MO(d->symZ, d->kappaZ, d->u_YT, 1.0, d->v_YT_ZTp);
    d->yy=d->zz=d->yzzy=0;
    for (isym=0; isym<inforb_.nsym; isym++) {
       nocci=inforb_.nocc[isym];
       if (nocci > 0) {
          iorbi = inforb_.iorb[isym];
    d->yy += (ddot_(&nocci, d->u_YN+iorbi,   &ONEI, grid->mov+iorbi, &ONEI) -
              ddot_(&nocci, grid->mov+iorbi, &ONEI, d->u_YT+iorbi,   &ONEI));
    d->zz += (ddot_(&nocci, d->u_ZN+iorbi,   &ONEI, grid->mov+iorbi, &ONEI) -
              ddot_(&nocci, grid->mov+iorbi, &ONEI, d->u_ZT+iorbi,   &ONEI));
    d->yzzy += 
	0.5*(+  ddot_(&nocci, d->v_YN_ZNp+iorbi,&ONEI, grid->mov+iorbi,&ONEI)
             +  ddot_(&nocci, d->v_YT_ZTp+iorbi,&ONEI, grid->mov+iorbi,&ONEI)
             -2*ddot_(&nocci, d->u_YN+iorbi,    &ONEI, d->u_ZT+iorbi,  &ONEI)
             -2*ddot_(&nocci, d->u_ZN+iorbi,    &ONEI, d->u_YT+iorbi,  &ONEI));
       }
    }
}

static __inline__ void
eval_grad_vars(DftGrid* grid, QuadFastData* d)
{
    int x, i, isym;
    
    matn_times_vec3_MO(d->symY, d->kappaY, grid->mog, 0.0, d->g_YN);
    matt_times_vec3_MO(d->symY, d->kappaY, grid->mog, 0.0, d->g_YT);
    matn_times_vec3_MO(d->symZ, d->kappaZ, grid->mog, 0.0, d->g_ZN);
    matt_times_vec3_MO(d->symZ, d->kappaZ, grid->mog, 0.0, d->g_ZT);
    matn_times_vec3_MO(d->symY, d->kappaY, d->g_ZN,   0.0, d->x_YN_ZNp);
    matt_times_vec3_MO(d->symY, d->kappaY, d->g_ZT,   0.0, d->x_YT_ZTp);
    matn_times_vec3_MO(d->symZ, d->kappaZ, d->g_YN,   1.0, d->x_YN_ZNp);
    matt_times_vec3_MO(d->symZ, d->kappaZ, d->g_YT,   1.0, d->x_YT_ZTp);
    d->trYsum = d->trZsum = d->trYZZYsum = 0;
    d->trYtimesZ = 0;
    for(x=0; x<3; x++) {
        int ix = x*inforb_.norbt;
        d->trY[x] = d->trZ[x] = d->trYZZY[x] = 0;
        for (isym=0; isym<inforb_.nsym; isym++) {
           for(i=inforb_.iorb[isym]; i<inforb_.iorb[isym]+inforb_.nocc[isym];i++) {
               d->trY[x]    += (grid->mog[ix+i]*(d->u_YN[i]-d->u_YT[i])+
                                grid->mov[i]*(d->g_YN[ix+i]-d->g_YT[ix+i]));
               d->trZ[x]    += (grid->mog[ix+i]*(d->u_ZN[i]-d->u_ZT[i])+
                                grid->mov[i]*(d->g_ZN[ix+i]-d->g_ZT[ix+i]));
               d->trYZZY[x] += 0.5*(grid->mog[ix+i]*(d->v_YN_ZNp[i]+
                                                     d->v_YT_ZTp[i])+
                                    grid->mov[i]*(d->x_YN_ZNp[ix+i]+
                                                  d->x_YT_ZTp[ix+i])-
                                    2*d->u_YN[i]*d->g_ZT[ix+i]-
                                    2*d->g_YN[ix+i]*d->u_ZT[i]-
                                    2*d->u_ZN[i]*d->g_YT[ix+i]-
                                    2*d->g_ZN[ix+i]*d->u_YT[i]);
           }
        }
        
	d->trYsum    += d->trY[x]*d->grada[x];
	d->trZsum    += d->trZ[x]*d->grada[x];
	d->trYZZYsum += d->trYZZY[x]*d->grada[x]; 
	d->trYtimesZ += d->trY[x]*d->trZ[x];
    }
}

/* add_dft_contribution:
   computes dft contribution to the quadratic response using precomputed
   traces of the operators (trY, trZ, etc). 
   When computing the contributions to the dftcontr matrix, we use
   intentionally aligned transformed vectors in order to utilize
   DGEMV to contract four or three vectors into one column of the outcome
   matrix. The gain is rather for small vectors fitting in cache.
   FIXME: contracting with GGA could use the same approach - it would provide
   more gain then.
*/
static void
add_dft_contribution(DftGrid* grid, QuadFastData* d)
{
    static const real sgn[2]={1.0,-1.0};
    int j, x, p, q;
    real pref, pref2b, pref2c, pref3;
    real *mov = grid->mov; /* just a shortcut */
    real *dftcontr = d->dftcontr;
    ThirdDrv drvs; /* the functional derivatives */
    int sY = d->ispinY, sZ = d->ispinZ;

    dftpot2_(&drvs, grid->curr_weight, &grid->dp,
             grid->dogga, d->ispinY != d->ispinZ);
    /* and first order contribution */
    pref = 0.5*drvs.fR;
    pref2b = d->yy*drvs.fRR[sY] + 2*d->trYsum*drvs.fRZ[sY] +
	d->trYsum*drvs.fRG[sY];
    pref2c = d->zz*drvs.fRR[sZ] + 2*d->trZsum*drvs.fRZ[sZ] +
	d->trZsum*drvs.fRG[sZ];

    /* third order, and part of second order contribution */
     pref3 = d->yy*d->zz*drvs.fRRR[sY|sZ] +
        2*(drvs.fRRZ[sZ][sY]*d->yy*d->trZsum +
           drvs.fRRZ[sY][sZ]*d->zz*d->trYsum) +
	4*d->trZsum*d->trYsum*drvs.fRZZ[sY^sZ][sY] +
	d->yzzy*drvs.fRR[sY^sZ] +
        2*(d->trYtimesZ+d->trYZZYsum)*drvs.fRZ[sY^sZ];
     pref3 += drvs.fRRG[sY]*d->yy*d->trZsum*(1+sgn[sZ])+
             drvs.fRRG[sZ]*d->zz*d->trYsum*(1+sgn[sY])+
             0.5*drvs.fRG[0]*(1+sgn[sY]*sgn[sZ])*d->trYZZYsum +
  	     0.5*drvs.fRG[0]*(sgn[sY]+sgn[sZ])*d->trYtimesZ;

    /* #define OLD_CODE */
#ifndef OLD_CODE
    dcopy_(&inforb_.norbt, grid->mov, &ONEI, d->mov, &ONEI);
#endif

    for(j=0; j<inforb_.norbt; j++) {
#ifdef OLD_CODE
        real f_ZN = mov[j]*pref2b, f_YN = mov[j]*pref2c;
        real f_mov = mov[j]*pref3 - (d->u_ZT[j]*pref2b + d->u_YT[j]*pref2c);
        int i;
        if(d->addfock) {
            f_YN  += -2*d->u_ZT[j]*pref;
            f_ZN  += -2*d->u_YT[j]*pref;
            f_mov += d->v_YT_ZTp[j]*pref;
        }
	for(i=0; i<inforb_.norbt; i++) {
            if(d->addfock)
                dftcontr[i+j*inforb_.norbt] += 
                    d->v_YN_ZNp[i]*mov[j]*pref;
	    /* second order contribution */
	    dftcontr[i+j*inforb_.norbt] += 
                f_ZN*d->u_ZN[i] + f_YN*d->u_YN[i] + f_mov*mov[i];
	}
#else
        real fac[4];
        fac[OFF_MOV] = mov[j]*pref3 - (d->u_ZT[j]*pref2b + d->u_YT[j]*pref2c);
        fac[OFF_YN]  = mov[j]*pref2c;
        fac[OFF_ZN]  = mov[j]*pref2b;
        if(d->addfock) {
            fac[OFF_YN]  += -2*d->u_ZT[j]*pref;
            fac[OFF_ZN]  += -2*d->u_YT[j]*pref;
            fac[OFF_MOV] += d->v_YT_ZTp[j]*pref;
            fac[OFF_YNZN] = mov[j]*pref;
            dgemv_("N", &inforb_.norbt, &FOURI, &ONER, d->mov,
                   &inforb_.norbt, fac, &ONEI, &ONER, 
                   &dftcontr[j*inforb_.norbt], &ONEI);
        } else 
            dgemv_("N", &inforb_.norbt, &THREEI,&ONER, d->mov,
                   &inforb_.norbt, fac, &ONEI, &ONER, 
                   &dftcontr[j*inforb_.norbt], &ONEI);
#endif
    }

    if(!grid->dogga) return;
    /* now the same thing has to be done for the grho matrices...
     * "It's a long way to the top if you want to rock'n'roll!" */
    for(x=0; x<3; x++) {
        int ix = x*inforb_.norbt;
	/* distribute two-index transformed densities */
        if(d->addfock) {
            pref = (drvs.fZ+0.5*drvs.fG)*d->grada[x];
            for(q=0; q<inforb_.norbt; q++) {
                int dq=q*inforb_.norbt;
                for(p=0; p<inforb_.norbt; p++) {
                    /* FIXME: cache conflicts? */
                    real ompq =
                        (+d->v_YN_ZNp[p]*grid->mog[ix+q]
                         +grid->mog[ix+p]*d->v_YT_ZTp[q]
                         +d->x_YN_ZNp[ix+p]*grid->mov[q]
                         +grid->mov[p]*d->x_YT_ZTp[ix+q]
                         -2*d->u_YN[p]*d->g_ZT[ix+q]
                         -2*d->g_YN[ix+p]*d->u_ZT[q]
                         -2*d->u_ZN[p]*d->g_YT[ix+q]
                         -2*d->g_ZN[ix+p]*d->u_YT[q]);
                    dftcontr[p+dq] += ompq*pref;
                }
            }
        }

       	/* distribute one-index transformed densities [k, rho] */
	pref2b = (2*d->yy*drvs.fRZ[sY] +4*d->trYsum*drvs.fZZ[sY])*d->grada[x]
	       +2*drvs.fZ*d->trY[x];
        pref2b += d->yy*drvs.fRG[sY]*d->grada[x] + sgn[sY]*d->trY[x]*drvs.fG;  

	pref2c  = (2*d->zz*drvs.fRZ[sZ] + 4*d->trZsum*drvs.fZZ[sZ])*d->grada[x]
	        +2*drvs.fZ*d->trZ[x];
        pref2c += d->zz*drvs.fRG[sZ]*d->grada[x] + sgn[sZ]*d->trZ[x]*drvs.fG;

	/*  distribute grho_x */
	pref3 =
	    (8*drvs.fZZZ[sY|sZ]*d->trYsum*d->trZsum +
             4*(drvs.fRZZ[sY][sZ]*d->yy*d->trZsum + 
                drvs.fRZZ[sZ][sY]*d->zz*d->trYsum)+
             2*drvs.fRRZ[sY^sZ][sY]*d->yy*d->zz)*d->grada[x];
        /* gamma 1st term */
        pref3 += drvs.fRRGX[sY][sZ]*d->yy*d->zz*d->grada[x];
        pref3 +=
            4*drvs.fZZ[sY^sZ]*d->trYZZYsum*d->grada[x] + 
            4*( drvs.fZZ[sY]*d->trZ[x]*d->trYsum + 
                drvs.fZZ[sZ]*d->trY[x]*d->trZsum + 
                drvs.fZZ[sY^sZ]*d->trYtimesZ*d->grada[x]
                  ) + 
            2*(drvs.fRZ[sY]*d->yy*d->trZ[x]
              +drvs.fRZ[sZ]*d->zz*d->trY[x])+
	      2*drvs.fRZ[sY^sZ]*d->yzzy*d->grada[x];
        /*gamma 2st term*/ 
        pref3 += drvs.fRG[sY^sZ]*d->yzzy*d->grada[x]+
                 drvs.fRG[sY]*d->yy*d->trZ[x]*sgn[sZ]+
  	         drvs.fRG[sZ]*d->zz*d->trY[x]*sgn[sY]; 
        pref3 += 2*drvs.fZ*d->trYZZY[x];
	/* gamma 3st term */
	pref3 += sgn[sY]*sgn[sZ]*drvs.fG*d->trYZZY[x]; 

        for(q=0; q<inforb_.norbt; q++) {
            int dq=q*inforb_.norbt;
            for(p=0; p<inforb_.norbt; p++) {
		/* distribute one-index transformed densities [k, rho] */
                dftcontr[p+dq] += 
                    (+d->u_ZN[p]*grid->mog[ix+q]
		     +d->g_ZN[ix+p]*grid->mov[q]
		     -grid->mog[ix+p]*d->u_ZT[q]
		     -grid->mov[p]*d->g_ZT[ix+q])*pref2b;

                dftcontr[p+dq] += 
                    (+d->u_YN[p]*grid->mog[ix+q]
		     +d->g_YN[ix+p]*grid->mov[q]
		     -grid->mog[ix+p]*d->u_YT[q]
		     -grid->mov[p]*d->g_YT[ix+q])*pref2c;

		/*  distribute grho_x */
                dftcontr[p+dq] += 
                    (+grid->mog[ix+p]*grid->mov[q]
		     +grid->mov[p]*grid->mog[ix+q])*pref3;
            }
        }
    }
}

static void
fast_callback(DftGrid* grid, QuadFastData* data)
{
    eval_rho_vars(grid, data);
    if(grid->dogga) {
	data->grada[0] = grid->grada[0];
	data->grada[1] = grid->grada[1];
	data->grada[2] = grid->grada[2];
        eval_grad_vars(grid, data); 
    }
    add_dft_contribution(grid, data);
}


#if defined(VAR_MPI)
#include <mpi.h>
#include <our_extra_mpi.h>
#define MASTER_NO 0
/* dft_qr_resp_slave:
   this is a slave driver. It's task is to allocate memory needed by
   the main property evaluator (dftqrcf_ in this case) and call it.
*/
void
dft_qr_resp_slave(real* work, integer* lwork, const integer* iprint)
{
    real* fi    = malloc(inforb_.n2basx*sizeof(real));              /* OUT */
    real *cmo   = malloc(inforb_.norbt*inforb_.nbast*sizeof(real)); /* IN  */
    real *kappaY= malloc(inforb_.n2orbx*sizeof(real));              /* IN  */
    real *kappaZ= malloc(inforb_.n2orbx*sizeof(real));              /* IN  */
    integer addfock, symY, symZ, spinY, spinZ;                      /* IN  */
    dftqrcf_(fi, cmo, kappaY, &symY, &spinY, kappaZ, &symZ, &spinZ, 
	     &addfock, work, lwork);
    free(kappaZ);
    free(kappaY);
    free(cmo);
    free(fi);
}

static void
dft_qr_resp_sync_slaves(real* cmo, real* kappaY, real* kappaZ, integer* addfock,
                        integer* symY, integer* symZ, integer* spinY, integer* spinZ)
{
    static const SyncData sync_data[] = {
 	{ inforb_.nocc,   8, fortran_MPI_INT },
 	{ &inforb_.nocct, 1, fortran_MPI_INT },
 	{ &inforb_.nvirt, 1, fortran_MPI_INT },
    };
#ifdef C99_COMPILER
    const SyncData data2[] = {
	{ cmo,     inforb_.norbt*inforb_.nbast,MPI_DOUBLE },
	{ kappaY,  inforb_.n2orbx,             MPI_DOUBLE },
	{ kappaZ,  inforb_.n2orbx,             MPI_DOUBLE },
	{ addfock, 1,                          fortran_MPI_INT    },
	{ symY,    1,                          fortran_MPI_INT    },
	{ symZ,    1,                          fortran_MPI_INT    },
	{ spinY,   1,                          fortran_MPI_INT    },
	{ spinZ,   1,                          fortran_MPI_INT    }
    };
#else /* C99_COMPILER */
    /* this is more error-prone but some compilers (HP/UX)... */
    static SyncData data2[] = 
    { {NULL, 0, MPI_DOUBLE}, {NULL, 0, MPI_DOUBLE}, {NULL, 0, MPI_DOUBLE}, 
      {NULL, 0, fortran_MPI_INT   }, {NULL, 0, fortran_MPI_INT   }, {NULL, 0, fortran_MPI_INT   },
      {NULL, 0, fortran_MPI_INT   }, {NULL, 0, fortran_MPI_INT   } };
    data2[0].data = cmo;     data2[0].count = inforb_.norbt*inforb_.nbast;
    data2[1].data = kappaY;  data2[1].count = inforb_.n2orbx;
    data2[2].data = kappaZ;  data2[2].count = inforb_.n2orbx;
    data2[3].data = addfock; data2[3].count = 1;
    data2[4].data = symY;    data2[4].count = 1;
    data2[5].data = symZ;    data2[5].count = 1;
    data2[6].data = spinY;   data2[6].count = 1;
    data2[7].data = spinZ;   data2[7].count = 1;
#endif /* C99_COMPILER */

    mpi_sync_data(sync_data, ELEMENTS(sync_data));
    mpi_sync_data(data2,     ELEMENTS(data2));
}
static __inline__ void
dft_qr_resp_collect_info(real* fi, real*work, integer lwork)
{
    int sz = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &sz);
    if(sz<=1) return;
    CHECK_WRKMEM(inforb_.n2orbx,lwork);
    dcopy_(&inforb_.n2orbx, fi, &ONEI, work, &ONEI);
    MPI_Reduce(work, fi, inforb_.n2orbx, MPI_DOUBLE, MPI_SUM, 
	       MASTER_NO, MPI_COMM_WORLD);
}
#else /* VAR_MPI */
#define dft_qr_resp_sync_slaves(cmo,kappaY,kappaZ,addfck,symY,symZ,spinY,spinZ)
#define dft_qr_resp_collect_info(fi,work,lwork)
#endif /* VAR_MPI */

/* this is the routine for computing the DFT exchange-correlation 
   contribution to quadratic response.
   NOTES:
   addfock - should Fock-like contribution [k,[k,o]] be added? 
*/
void
dftqrcf_(real* fi, real* cmo, real* kappaY, integer* symY, integer* spinY, 
         real* kappaZ, integer* symZ, integer* spinZ, integer* addfock, 
         real* work, integer* lwork)
{
    static int msg_printed = 0;
    struct tms starttm, endtm; clock_t utm;
    DftCallbackData cbdata[1];
    real electrons;
    QuadFastData* data;

    /* WARNING: NO work MAY BE done before syncing slaves! */
    dft_wake_slaves((DFTPropEvalMaster)dftqrcf_);       /* NO-OP in serial */
    dft_qr_resp_sync_slaves(cmo,kappaY,kappaZ,addfock,
                            symY,symZ, spinY,spinZ); 

    if(!msg_printed) {
        fort_print("DFT-QR uses %s matrix multipication code.\n",
                   mm_code_version);
        msg_printed = 1;
    }

    data = quadfast_data_new(kappaY, *symY, *spinY, kappaZ, *symZ, *spinZ,
                             *addfock); 

    cbdata[0].callback = (DftCallback)fast_callback;
    cbdata[0].cb_data  = data;
    times(&starttm);

    electrons = dft_integrate(cmo, work, lwork, cbdata, ELEMENTS(cbdata));
    
    dft_qr_resp_collect_info(data->dftcontr, work,*lwork); /* NO-OP in serial */
    daxpy_(&inforb_.n2orbx, &ONER, data->dftcontr, &ONEI, fi, &ONEI);
    if(inforb_.norbt<=4) {
        fort_print("Dumping DFT quadratic contribution");
        outmat_(data->dftcontr, &ONEI, &inforb_.norbt, &ONEI, &inforb_.norbt,
                &inforb_.norbt, &inforb_.norbt);
    }
    quadfast_data_free(data);
    times(&endtm);
    utm = endtm.tms_utime-starttm.tms_utime;
    fort_print("      Electrons: %15.7f. Quadratic response time: %10.2f s\n", 
               electrons, utm/(double)sysconf(_SC_CLK_TCK));
}

#ifdef FAST_TEST
#include "integrator.h"

struct common_inforb inforb_;

static void
dump_grY(char* name, DftGrid*grid, QuadFastData* d, int ix, int dim)
{
    int p, q;
    printf("%s:\n", name);
    for(p=0; p<dim; p++) {
	for(q=0; q<dim; q++)
	    printf("%10.7f ", 
		   +d->u_YN[p]*grid->mog[ix+q]
		   -grid->mog[ix+p]*d->u_YT[q]
		   +d->g_YN[ix+p]*grid->mov[q]
		   -grid->mov[p]*d->g_YT[ix+q]);
	puts("");
    }
    puts("");
}    
int main(int argc, char* argv[])
{
    static const real kappaY[] = { 0, 1, 1, 0};
    static const real kappaZ[] = { 0, 0, 0, 0};
    DftGrid *grid;
    QuadFastData* data;
    
    inforb_.norbt = 2;
    inforb_.nocct = 1;

    grid = dft_grid_new(&Example3Functional); 
    printf("main - grid->mov: %p\n", grid->mov);
    data = quadfast_data_new(kappaY, kappaZ);
    grid->grad[0] = grid->grad[1] = grid->grad[2] = 0;
    grid->ngrad = sqrt(2);
    grid->mov[0] = 1; grid->mov[1] = 0.5;
    grid->mog[0+0] = 1; grid->mog[0+1] = 0.5;
    grid->mog[2+0] = 0; grid->mog[2+1] = 0;
    grid->mog[4+0] = 0; grid->mog[4+1] = 0;

    grid->curr_weight = 1;

    fast_callback(grid, data);

    dump_grY("[Y,grho_x]", grid, data, 0, inforb_.norbt);
    dump_grY("[Y,grho_y]", grid, data, 4, inforb_.norbt);
    dump_grY("[Y,grho_z]", grid, data, 8, inforb_.norbt);
    dft_grid_free(grid);
    quadfast_data_free(data);
    return 0;
}
#endif
