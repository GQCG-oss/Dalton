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
/* quadra-strict.c:
   (c) Pawel Salek, pawsa@theochem.kth.se, feb 2002

   The DFT quadratic response subroutines that follow string matrix
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

#define __CVERSION__

#include "integrator.h"
#include "functionals.h"

#include "inforb.h"

#if !defined __inline__
/* inline some stuff whenever possible */
#define __inline__
#endif

/* the temporary data used for integration */
struct QuadStrictData_ {
    /* pointers to external data (data is not owned) */
    const real* kappaY, *kappaZ; 
    integer ispinY,ispinZ; /* rank of kappa vectors */
    integer addfock; /* add fock contribution */

    /* temporary variables */
    real grada[3]; /* alpha density gradients */
    real pref2sum; /* this is used for checksumming, debugging etc. */
    real yy, zz, yzzy;
    real trY[3], trZ[3]; /* trY[x] = [kappa, grho_x])*/
    real trYsum, trZsum; /* = sum trY[x]*grad_x/|grho| */
    real trYZZYsum; 
    real trYZ[3], trZY[3]; /* trYZ[x] = [k_Z,[k_Y,grho_x]])*grad_x/|grho| */
    real trYZZY[3];  /* trYZZY=<Tr([kappa,[kappa, grho_x]])*grad_x/|grho| */
    real trYtimesZ; /* = sum [kY, grho_x]*[kappaZ, grho_x]*(grho_x/|grho|)^2 */

    /* GGA temporary variables */
    real *grho; /* grho_x = phi_p^x*phi_q^x */
    real *grY;   /* grho(kappaY) = [kappaY, grho] (commutator) */
    real *grZ;   /* grho(kappaZ) = [kappaZ, grho] (commutator) */
    real *grYZZY; /* grho(kappaY,kappaZ)+grho(kappaZ,kappaY) */
    real trYtZ;   /* =sum_xy Tr([k, grho_y])Tr([k, grho_z])
		   *  *grho_x*grho_y/|grho|^2 */
    real trYZg;   /* sum_x Tr([k, grho_y])Tr([k, grho_y])/|grho| */

    /* Allocatable arrays */
    real* dftcontr;
    real *omega; /* omega = phi_p*phi_q */
    real *omY;   /* omega(kappaY) = [kappaY, omega] (commutator) */
    real *omZ;   /* omega(kappaZ) = [kappaZ, omega] (commutator) */
    real *omYZ; /* omega(kappaY,kappaZ) */
    real *omZY; /* omega(kappaZ,kappaY) */

    real* d1, *d2, *d3, *d4;
};

typedef struct QuadStrictData_ QuadStrictData;

static QuadStrictData*
quadstrict_data_new(const real* kY, const real* kZ, integer ispinY, integer ispinZ,
                    integer addfock)
{
    QuadStrictData* res = malloc(sizeof(QuadStrictData));

    res->kappaY = kY;     res->kappaZ = kZ;
    res->ispinY = ispinY; res->ispinZ = ispinZ;
    res->addfock = addfock;
    res->pref2sum = 0;

    res->trYsum = res->trZsum = 0; /* initialize them for dogga == FALSE */
    res->trYtimesZ = 0;

    res->dftcontr = alloc_mat_MO(1);
    res->omega = alloc_mat_MO(1);
    res->omY  = alloc_mat_MO(1);
    res->omZ  = alloc_mat_MO(1);
    res->omYZ = alloc_mat_MO(1);
    res->omZY = alloc_mat_MO(1);
    res->grho = alloc_mat_MO(3);
    res->grY  = alloc_mat_MO(3);
    res->grZ  = alloc_mat_MO(3);
    res->grYZZY = alloc_mat_MO(3);

    res->d1 = alloc_mat_MO(1);
    res->d2 = alloc_mat_MO(1);
    res->d3 = alloc_mat_MO(1);
    res->d4 = alloc_mat_MO(1);
    return res;
}

static void 
quadstrict_data_free(QuadStrictData* tmp)

{
    free(tmp->dftcontr);
    free(tmp->omega);
    free(tmp->omY);  free(tmp->omZ);
    free(tmp->omYZ); free(tmp->omZY); 
    free(tmp->grho);
    free(tmp->grY);  free(tmp->grZ);
    free(tmp->grYZZY);

    free(tmp->d1);
    free(tmp->d2);
    free(tmp->d3);
    free(tmp->d4);
    free(tmp);
}

static __inline__ void
commute_matrices(real alpha, const real* a, const real* b, real* c, int addp)
{
    const real* firstpref = addp ? &ONER : &ZEROR;
    dgemm_("N", "N", &inforb_.norbt, &inforb_.norbt, &inforb_.norbt, &alpha, 
	   a, &inforb_.norbt, b, &inforb_.norbt, firstpref, c, &inforb_.norbt);
    alpha = -alpha; 
    dgemm_("N", "N", &inforb_.norbt, &inforb_.norbt, &inforb_.norbt, &alpha, 
           b, &inforb_.norbt, a, &inforb_.norbt, &ONER, c, &inforb_.norbt);
}

static real inactive_trace(real *mat)
{
   int isym,symoff;
   integer stride;
   real result = 0;
   symoff =0;
   stride=inforb_.norbt+1;
   for (isym=0; isym<inforb_.nsym; isym++) {
      if (inforb_.nocc[isym] > 0) {
         result += dsum_(&inforb_.nocc[isym], mat+symoff, &stride);
      }
      symoff += inforb_.norbt*inforb_.norb[isym] + inforb_.norb[isym];
   }
   return result;
}
         

/* eval_rho_vars: rho-dependent temporary vars */
static void
eval_rho_vars(DftGrid *grid, QuadStrictData* d)
{
    dgemm_("N", "N", &inforb_.norbt, &inforb_.norbt, &ONEI, &ONER, 
	   grid->mov, &inforb_.norbt, grid->mov, &ONEI, 
	   &ZEROR, d->omega, &inforb_.norbt); 
    /* dger won't do: it does not reset the matrix first */

    commute_matrices(1.0, d->kappaY, d->omega, d->omY, 0);
    commute_matrices(1.0, d->kappaZ, d->omega, d->omZ, 0);
    commute_matrices(1.0, d->kappaY, d->omZ, d->omYZ, 0);
    commute_matrices(1.0, d->kappaZ, d->omY, d->omZY, 0);
    d->yy = inactive_trace(d->omY);
    d->zz = inactive_trace(d->omZ);
    d->yzzy = 0.5*(inactive_trace(d->omYZ) + 
                   inactive_trace(d->omZY));
}

static void
eval_grad_vars(DftGrid *grid, QuadStrictData* d)
{
    int x, stride, norbt2;

    norbt2 = inforb_.norbt*inforb_.norbt;
    /* compute grho, one-index grY and grZ and double-index 
     * grYZZY... */
    for(x=0; x<3; x++) {
	int idx = x*norbt2;
	dgemm_("N", "N", &inforb_.norbt, &inforb_.norbt, &ONEI, &ONER, 
	       &grid->mog[x*inforb_.norbt], &inforb_.norbt, 
	       grid->mov, &ONEI, 
	       &ZEROR, &d->grho[idx], 
	       &inforb_.norbt); 
	dgemm_("N", "N", &inforb_.norbt, &inforb_.norbt, &ONEI, &ONER, 
	       grid->mov, &inforb_.norbt, 
	       &grid->mog[x*inforb_.norbt], &ONEI, 
	       &ONER, &d->grho[idx], 
	       &inforb_.norbt); 
	commute_matrices(1.0, d->kappaY, &d->grho[idx], &d->grY[idx], 0);
	commute_matrices(1.0, d->kappaZ, &d->grho[idx], &d->grZ[idx], 0);
	commute_matrices(1.0, d->kappaY, &d->grZ[idx],  &d->grYZZY[idx], 0);
	commute_matrices(1.0, d->kappaZ, &d->grY[idx],  &d->grYZZY[idx], 1);
    }
    /* ... and some traces... */
    stride = inforb_.norbt+1;

    d->trYsum = d->trZsum = d->trYZZYsum = 0;
    d->trYtimesZ = 0;
    for(x=0; x<3; x++) {
	d->trY[x]    = inactive_trace(&d->grY[x*norbt2]);
	d->trZ[x]    = inactive_trace(&d->grZ[x*norbt2]); 
	d->trYZZY[x] = 
            0.5*inactive_trace(&d->grYZZY[x*norbt2]);

	d->trYsum    += d->trY[x]*d->grada[x];
	d->trZsum    += d->trZ[x]*d->grada[x];
	d->trYZZYsum += d->trYZZY[x]*d->grada[x]; 
	d->trYtimesZ += d->trY[x]*d->trZ[x];
    }
    /* fort_print("%10f %12.8g %12.8g %12.8g %12.8g", grid->ngrad,
       d->trYsum, d->trZsum, d->trYZZYsum, d->trYtimesZ); */
}

static void
add_dft_contribution(DftGrid* grid, QuadStrictData* d)
{
    static const real sgn[2]={1.0,-1.0};
    int stride, x;
    integer norbt2 = inforb_.norbt*inforb_.norbt;
    real pref;
    real* dftcontr = d->dftcontr;
    ThirdDrv drvs; /* the functional derivatives */
    int sY = d->ispinY, sZ = d->ispinZ;

    dftpot2_(&drvs, grid->curr_weight, &grid->dp,
	     grid->dogga, d->ispinY != d->ispinZ);

    stride = inforb_.norbt+1;
    /* distribute two-index transformed densities */
    if(d->addfock) {
        pref = 0.5*drvs.fR;
        daxpy_(&norbt2, &pref, d->omYZ, &ONEI, dftcontr, &ONEI);
        daxpy_(&norbt2, &pref, d->omZY, &ONEI, dftcontr, &ONEI);
    }

    /* distribute one-index transformed densities [k, rho] */
    pref = d->yy*drvs.fRR[sY] + 2*d->trYsum*drvs.fRZ[sY] +
           d->trYsum*drvs.fRG[sY];
    daxpy_(&norbt2, &pref, d->omZ, &ONEI, dftcontr, &ONEI);
    pref = d->zz*drvs.fRR[sZ] + 2*d->trZsum*drvs.fRZ[sZ] +
           d->trZsum*drvs.fRG[sZ];
    daxpy_(&norbt2, &pref, d->omY, &ONEI, dftcontr, &ONEI);

    /*  distribute rho */
     pref = d->yy*d->zz*drvs.fRRR[sY|sZ] +
        2*(drvs.fRRZ[sZ][sY]*d->yy*d->trZsum +
           drvs.fRRZ[sY][sZ]*d->zz*d->trYsum) +
	4*d->trZsum*d->trYsum*drvs.fRZZ[sY^sZ][sY] +
	d->yzzy*drvs.fRR[sY^sZ] +
        2*(d->trYtimesZ+d->trYZZYsum)*drvs.fRZ[sY^sZ];
     pref += drvs.fRRG[sY]*d->yy*d->trZsum*(1+sgn[sZ])+
             drvs.fRRG[sZ]*d->zz*d->trYsum*(1+sgn[sY])+
             0.5*drvs.fRG[0]*(1+sgn[sY]*sgn[sZ])*d->trYZZYsum +
  	     0.5*drvs.fRG[0]*(sgn[sY]+sgn[sZ])*d->trYtimesZ;
     daxpy_(&norbt2, &pref, d->omega, &ONEI, dftcontr, &ONEI); 

    if(!grid->dogga) return; /* happilly home! */

    /* now the same thing has to be done for the grho matrices...
     * "It's a long way to the top if you want to rock'n'roll!" */
    for(x=0; x<3; x++) {
	/* distribute two-index transformed densities */
        if(d->addfock) {
            pref = (drvs.fZ+0.5*drvs.fG)*d->grada[x];
            daxpy_(&norbt2, &pref,&d->grYZZY[x*norbt2], &ONEI, 
                   dftcontr, &ONEI);
        }

       	/* distribute one-index transformed densities [k, rho]*/
	pref = (2*d->yy*drvs.fRZ[sY] +4*d->trYsum*drvs.fZZ[sY])*d->grada[x]
	       +2*drvs.fZ*d->trY[x];
        pref += d->yy*drvs.fRG[sY]*d->grada[x] + sgn[sY]*d->trY[x]*drvs.fG;  
	daxpy_(&norbt2, &pref, &d->grZ[x*norbt2], &ONEI, dftcontr, &ONEI);

	pref  = (2*d->zz*drvs.fRZ[sZ] + 4*d->trZsum*drvs.fZZ[sZ])*d->grada[x]
	        +2*drvs.fZ*d->trZ[x];
        pref += d->zz*drvs.fRG[sZ]*d->grada[x] + sgn[sZ]*d->trZ[x]*drvs.fG;
	daxpy_(&norbt2, &pref, &d->grY[x*norbt2], &ONEI, dftcontr, &ONEI); 
	/*  distribute grho_x  */
	/* the prefix is split into 3 subterms ordered by derivatives */
	pref =
	    (8*drvs.fZZZ[sY|sZ]*d->trYsum*d->trZsum +
             4*(drvs.fRZZ[sY][sZ]*d->yy*d->trZsum + 
                drvs.fRZZ[sZ][sY]*d->zz*d->trYsum)+
             2*drvs.fRRZ[sY^sZ][sY]*d->yy*d->zz)*d->grada[x];
        /* gamma 1st term */
        pref += drvs.fRRGX[sY][sZ]*d->yy*d->zz*d->grada[x];
        pref +=
            4*drvs.fZZ[sY^sZ]*d->trYZZYsum*d->grada[x] + 
            4*( drvs.fZZ[sY]*d->trZ[x]*d->trYsum + 
                drvs.fZZ[sZ]*d->trY[x]*d->trZsum + 
                drvs.fZZ[sY^sZ]*d->trYtimesZ*d->grada[x]
                  ) + 
            2*(drvs.fRZ[sY]*d->yy*d->trZ[x]
              +drvs.fRZ[sZ]*d->zz*d->trY[x])+
	      2*drvs.fRZ[sY^sZ]*d->yzzy*d->grada[x];
        /*gamma 2st term*/ 
        pref += drvs.fRG[sY^sZ]*d->yzzy*d->grada[x]+
                 drvs.fRG[sY]*d->yy*d->trZ[x]*sgn[sZ]+
  	         drvs.fRG[sZ]*d->zz*d->trY[x]*sgn[sY]; 
        pref += 2*drvs.fZ*d->trYZZY[x];
	/* gamma 3st term */
	pref += sgn[sY]*sgn[sZ]*drvs.fG*d->trYZZY[x]; 
	daxpy_(&norbt2, &pref, &d->grho[x*norbt2], &ONEI, dftcontr, &ONEI); 
    }
}

static void
dump_input_matrices(real* fi, real* kappaY, real* kappaZ, real* cmo)
{
    fort_print("Dumping FI");
    outmat_(fi, &ONEI, &inforb_.norbt, &ONEI, &inforb_.norbt,
    &inforb_.norbt, &inforb_.norbt);
    fort_print("Dumping kappaY");
    outmat_(kappaY, &ONEI, &inforb_.norbt, &ONEI, &inforb_.norbt,
    &inforb_.norbt, &inforb_.norbt);
    fort_print("Dumping kappaZ");
    outmat_(kappaZ, &ONEI, &inforb_.norbt, &ONEI, &inforb_.norbt,
    &inforb_.norbt, &inforb_.norbt);
    fort_print("Dumping CMO");
    outmat_(cmo, &ONEI, &inforb_.nbast, &ONEI, &inforb_.norbt,
    &inforb_.nbast, &inforb_.norbt);
}

/* MAIN: */

void
strict_callback(DftGrid* grid, QuadStrictData* data)
{
    eval_rho_vars(grid, data);
    if(grid->dogga) {
	/* We use here alpha and beta density gradients */
	data->grada[0] = grid->grada[0];
	data->grada[1] = grid->grada[1];
	data->grada[2] = grid->grada[2];
	eval_grad_vars(grid, data);     
    }
    add_dft_contribution(grid, data);
}


/* this is the routine for computing the DFT exchange-correlation 
   contribution to quadratic response.
*/
void
dftqrcs_(real* fi, real* cmo, real* kappaY, integer* symY, integer* ispinY,
         real* kappaZ, integer* symZ, integer* ispinZ, integer* addfock, 
         real* work, integer* lwork)
{
    integer norbt2 = inforb_.norbt*inforb_.norbt;
    DftCallbackData cbdata[1];
    QuadStrictData* data;

    if(0) dump_input_matrices(fi, kappaY, kappaZ, cmo);
    data = quadstrict_data_new(kappaY, kappaZ, *ispinY, *ispinZ, *addfock);
    
    cbdata[0].callback = (DftCallback)strict_callback;
    cbdata[0].cb_data = data;

    dft_integrate(cmo, work, lwork, cbdata, ELEMENTS(cbdata));

    /* data->dftcontr[1] = data->dftcontr[2] = 0.016642368353978; */
    /* Example2: LiH, 1, STO-2G, Example:
       data->dftcontr[1] = data->dftcontr[2] = -0.02926486447825;  */

    daxpy_(&norbt2, &ONER, data->dftcontr, &ONEI, fi, &ONEI);

    fort_print("Total DFT quadratic contribution (Exp:-0.3972284925 )");
    /* -0.09690870714785695978: LiH, 1, STO-2G, Example*/
    outmat_(data->dftcontr, &ONEI, &inforb_.norbt, &ONEI, &inforb_.norbt,
	    &inforb_.norbt, &inforb_.norbt);
    fort_print("First term:");
    outmat_(data->d1, &ONEI, &inforb_.norbt, &ONEI, &inforb_.norbt,
	    &inforb_.norbt, &inforb_.norbt);
    fort_print("Second term:");
    outmat_(data->d2, &ONEI, &inforb_.norbt, &ONEI, &inforb_.norbt,
	    &inforb_.norbt, &inforb_.norbt);
    fort_print("Third term:");
    outmat_(data->d3, &ONEI, &inforb_.norbt, &ONEI, &inforb_.norbt,
	    &inforb_.norbt, &inforb_.norbt);
    fort_print("Fourth term:");
    outmat_(data->d4, &ONEI, &inforb_.norbt, &ONEI, &inforb_.norbt,
	    &inforb_.norbt, &inforb_.norbt);

    quadstrict_data_free(data);
}

#ifdef STRICT_TEST
#include "integrator.h"

struct common_inforb inforb_;

static void
dump_mat(char* name, real* mat, int dim)
{
    int i, j;
    printf("%s:\n", name);
    for(j=0; j<dim; j++) {
	int dj = dim*j;
	for(i=0; i<dim; i++)
	    printf("%10.7f ", mat[i+dj]);
	puts("");
    }
	puts("");
}    
int main(int argc, char* argv[])
{
    static const real kappaY[] = { 0, 1, 1, 0};
    static const real kappaZ[] = { 0, 0, 0, 0};
    ThirdDrv drvs; /* the functional derivatives */
    DftGrid *grid;
    QuadStrictData* data;
    
    inforb_.norbt = 2;
    inforb_.nocct = 1;

    grid = dft_grid_new(&Example3Functional); 
    printf("main - grid->mov: %p\n", grid->mov);
    data = quadstrict_data_new(kappaY, kappaZ);
    grid->grad[0] = grid->grad[1] = grid->grad[2] = 0;
    grid->ngrad = sqrt(2);
    grid->mov[0] = 1; grid->mov[1] = 0.5;
    grid->mog[0+0] = 1; grid->mog[0+1] = 0.5;
    grid->mog[2+0] = 0; grid->mog[2+1] = 0;
    grid->mog[4+0] = 0; grid->mog[4+1] = 0;

    drvs->fRRR = 1;
    grid->curr_weight = 1;

    strict_callback(grid, data);

    dump_mat("Omega",     data->omega,   inforb_.norbt);
    dump_mat("[Y,grho_x]", &data->grY[0], inforb_.norbt);
    dump_mat("[Y,grho_y]", &data->grY[4], inforb_.norbt);
    dump_mat("[Y,grho_z]", &data->grY[8], inforb_.norbt);
    dft_grid_free(grid);
    quadstrict_data_free(data);
    return 0;
}
#endif
