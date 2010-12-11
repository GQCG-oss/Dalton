/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* Quadratic response DFT contribution.
   Pawel Salek, pawsa@theochem.kth.se.
   2001.07.13
   This file could be merged with other DFT related files later on.
*/

#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#if defined(__CVERSION__) /* THE C VERSION OF THE HEADERS */

#include "general.h"
#include "functionals.h"

struct DftGrid_ {
    /* private to integrator */
    real (*coor)[3];
    real* weight;
    integer* ncnt;
    real* atv; /* the orbital values and their derivatives at given 
		* grid point. The vector is indexed by dftinf_.kso1, etc
		*/

    real dfthri; /* threshold on orbital values */
    real dfthr0; /* threshold on density values */
    /* public, read only */
    integer ntypso; /* how many different vectors are computed for each
		     * (point,orbital) pair. i.e whether only orbital values
		     * are computed (1), orbital values and first derivatives 
		     * (4), etc. */

    integer london_off; /* offset of the "london" orbital derivatives */
    /* 1 - only values; 4 - values + (x,y,z) derivatives, etc */
    FunDensProp dp;
    real grada[3];/* alpha, also used in closed-shell code */
    real gradb[3];/* beta, used only in open-shell calculations. */ 
    real *mov;
    real *mog;
    int  curr_point;  /* index of the current point */
    real curr_weight; /* the weight at current grid point */
    integer dogga;    /* whether the functional requires gradient correction*/
    int needgrad;     /* the property evaluator requires orbital derivatives*/
    int needlap;      /* whether second order orbital derivatives are needed*/
    integer needgb;
};

typedef void (*DftCallback)(DftGrid* grid, void* cb_data);
struct DftCallbackData_ {
    DftCallback callback;
    void* cb_data;
};

typedef struct DftCallbackData_ DftCallbackData;
/* dft_grid_new: initialize grid/integrator. inform whether
   the callback will need orbital derivatives on its own.
*/
DftGrid* dft_grid_new(int needgrad, int needlap, int needgb);
void dft_grid_free(DftGrid* res);
real dft_integrate(real* cmo, real* work, integer* lwork, 
		   const DftCallbackData* cbarr, int cbcount);
real dft_integrate_ao(DftDensity* dens, real* work, integer* lwork, 
                      int needgrad, int needlap, int needgb,
                      const DftCallbackData* cbarr, int cbcount);

/* =================================================================== */
/*                     BLOCKED INTEGRATORS                             */
/* =================================================================== */
/* Max block length. Should not be too short because the loop overhead
 * will grow too large, nor too long because we run out of cache
 * then. */
#define DFT_BLLEN 140

typedef struct DftIntegratorBl_ {
    /* private to integrator */
    real (*RESTRICT coor)[3];
    real* RESTRICT weight;
    real* RESTRICT atv; /* the orbital values and their derivatives at given 
		* grid point. The vector is indexed by dftinf_.kso1, etc
                * the dimensioning is (C syntax) [ntypso][nbast][bllen].
		*/
    real dfthri; /* threshold on orbital values */
    int nsym;
    integer shl_bl_cnt, bas_bl_cnt[8];
    integer (*RESTRICT shlblocks)[2]; /* shell blocks   */
    integer (*RESTRICT basblocks)[2]; /* basis function blocks */
#define BASBLOCK(grid,isym) ((grid)->basblocks + (isym)*(grid)->shl_bl_cnt)

    integer ntypso; /* how many different vectors are computed for each
                 * (point,orbital) pair. i.e whether only orbital values
                 * are computed (1), orbital values and first derivatives 
                 * (4), etc. */

    integer london_off; /* offset of the "london" orbital derivatives */
    /* 1 - only values; 4 - values + (x,y,z) derivatives, etc */

    int ndmat; /* 1 for closed shell, 2 for open shell */
    /* for closed shell, only rho is set. For open shell, only rhoa and rhob
     * is set. */
    union {
        real *rho; /* total density vector; used in closed shell code. */
        struct {  /* used in open-shell code.    */
            real *a, *b;
        }ho;
    }r;
    union {
        real (*grad)[3]; /*total density gradient; used in closed shell code.*/
        struct {
            real (*a)[3], (*b)[3];
        }rad;
    }g;
    /* public, read only */
    real tgrad[3];/* alpha, also used in closed-shell code */
    int  curr_point;  /* index of the current point */
    real curr_weight; /* the weight at current grid point */
    integer dogga, needlap, needgb;
} DftIntegratorBl;

/* dft_integrate_ao_bl:
   numerical integration in atomic orbitals, blocked scheme.
*/
typedef void (*DftBlockCallback)(DftIntegratorBl* grid, real *tmp, 
                                 int bllen, int blstart, int blend,
                                 void* cb_data);

real dft_integrate_ao_bl(int ndmat, real *dmat, real* work, integer* lwork,
                         int needlnd, DftBlockCallback cb, void *cb_data);

#else /* THE FORTRAN VERSION OF THE HEADERS */

#endif /* _C_VERSION_ */
#endif /* _INTEGRATOR_H_ */
