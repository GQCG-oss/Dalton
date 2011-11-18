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
/* The DFT integrator.
   (c) Pawel Salek, pawsa@theochem.kth.se.
   2001.07.13

   The WRKMEM memory block is not used since it should be deprecated.
   It might be therefore useful to enable memory overcommiting. On linux-2.4.x
   it can be done with echo 1 > /proc/sys/vm/overcommit_memory or a 
   sysctl call. We use it only to pass it to other Fortran routines we call.

   OPTIMIZATIONS: ordinary calculation uses approximately only 4%
   total CPU time in this code. Most likely, the optimizations should
   be sought somewhere else. The simple optimization path is though to
   use block structure of kappa matrices to reduce time by 2 for
   larger matrices.  

   integrator.c provides dft_integrator() routine. It is passed some
   standard DALTON parameters (CMO, WORK, LWORK), and a table of
   callbacks and associated closures (callback data). The callback
   gets the grid data for current point as well as its own closure.

   The grid file DALTON.QUAD is assumed to be available on call to
   dft_integrator(). OTherwise, it is a black-box implementation.

*/

/* strictly conform to XOPEN ANSI C standard */
#if !defined(SYS_DEC)
/* XOPEN compliance is missing on old Tru64 4.0E Alphas */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1
#endif

#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define __CVERSION__
#include "integrator.h"
#include "functionals.h"
#include "grid-gen.h"
#include "inforb.h"

#if 0 && defined(VAR_MPI)
#include <mpi.h>
#include "infpar.h"
#endif /* VAR_MPI */

#if !defined __inline__
/* inline some stuff whenever possible */
#define __inline__
#endif

/* blocksz_t is a variable that matches the one used by the 
 * fortran runtime library to store the block size. This type
 * is compilator dependent but is usually int or long int.
 */
#if defined(__gnu_linux__)
/* gnu compilers use this */
typedef long blocksz_t;
#else
/* safe default for all the other compilers */
typedef int blocksz_t;
#endif

#define max(a,b) ((a)>(b)? (a):(b))

/* Fortran function prototypes. They should probably go to a DALTON
 * header interface file.
 */

extern void FSYM(setupsos)(const integer* geodrv, const integer* dolnd,
			   integer *ntypso, integer *lndoff);
extern void FSYM(getsos)(real* atv, integer* ncnt, real* coorx,
			 real* work, integer* lwork, const integer* nbast, 
			 const integer* dolnd, const integer* dogga,
			 const real* dfthri, const integer* iprint);
extern void FSYM2(dft_get_thresholds)(real* dfthri, real *dfthr0);


/* alloc_mat_MO: 
 * allocate space for a matrix or a set of matrices
 * transforming from one set of molecular orbitals to others
 * (basically, norbt*norbt with a twist for multiple symmetries).  
 */
real*
alloc_mat_MO(int cnt)
{
  real* res = calloc(cnt*inforb_.n2orbx, sizeof(real));
  if(res == NULL) abort();
  return res;
}
  

/* dft_grid_new:
   initialize grid data.
*/
static const int GRID_BUFF_SZ = 100000;  /* the DFT grid buffer length */
DftGrid*
dft_grid_new(int needgrad, int needlap, int needgb)
{
    real dfthri, dfthr0;
    integer geodrv;
    DftGrid* grid = malloc(sizeof(DftGrid));
    
    grid->coor   = calloc(3*GRID_BUFF_SZ, sizeof(real));
    grid->weight = calloc(GRID_BUFF_SZ, sizeof(real));
    grid->ncnt   = calloc(inforb_.nbast, sizeof(int));
    grid->needlap= needlap;
    grid->needgb = needgb;
    grid->dogga  = selected_func->is_gga();
    grid->needgrad= needgrad || grid->dogga;
    FSYM2(dft_get_thresholds)(&dfthri, &dfthr0);
    grid->dfthri = dfthri;
    grid->dfthr0 = dfthr0;

    geodrv = grid->dogga ? 1 : 0;
    if(needlap) geodrv++;
    FSYM(setupsos)(&geodrv, &grid->needgb, &grid->ntypso, &grid->london_off);
    grid->atv   = calloc(grid->ntypso*inforb_.nbast, sizeof(real));
    grid->mov   = calloc(inforb_.norbt, sizeof(real));
    if(grid->dogga) grid->mog = calloc(3*inforb_.norbt, sizeof(real));
    else grid->mog = NULL;
    grid->grada[0] =  grid->grada[1] =  grid->grada[2] = 0;
    grid->gradb[0] =  grid->gradb[1] =  grid->gradb[2] = 0;
    grid->dp.rhoa =  grid->dp.grada = 0;
    grid->dp.rhob =  grid->dp.gradb = 0;

    return grid;
}

void
dft_grid_free(DftGrid* res)
{
  free(res->coor);
  free(res->weight);
  free(res->ncnt);
  free(res->atv);
  free(res->mov);
  if(res->mog) free(res->mog);
  free(res);
}

static __inline__ void
dft_grid_getval(DftGrid* grid, int ipnt, real* work, integer* lwork)
{
    real thrint = grid->dfthri/grid->weight[ipnt];
    FSYM(getsos)(grid->atv, grid->ncnt, &grid->coor[ipnt][0], 
		 work, lwork, &inforb_.nbast, &grid->needgb,
		 &grid->dogga, &thrint, &ONEI);
}

void 
printptr_(const char* str, void* ptr)
{
    fputs(str, stderr); fprintf(stderr, "%p\n", ptr);
}


#if 0 && defined(VAR_MPI)
static void
dft_integrate_collect_info(real *electrons)
{
    real tmp = *electrons;
    int sz = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &sz);
    if(sz<=1) return;
    MPI_Reduce(&tmp, electrons, 1, MPI_DOUBLE, MPI_SUM, 
	       daltoninfpar_.master, MPI_COMM_WORLD);
}
#else /* VAR_MPI */
#define dft_integrate_collect_info(electrons)
#endif /* VAR_MPI */

/* dft_integrate:
   Integrate KS wave function in MO basis using provided set of callbacks.
   Returns integrated charge.
   FIXME: adapt to open-shell scheme.
 */
real
dft_integrate(real* cmo, real* work, integer* lwork, 
	      const DftCallbackData* cbarr, int cbcount)
{
    int npoints, ipnt, cbno;
    int blocksz;
    real electrons;
    DftGrid* grid;
    int isym;
    integer nbasi, norbi, nocci;
    real *atvi, *movi, *cmoi, *mogi;
    DftGridReader* rawgrid;

    grid = dft_grid_new(0, 0, 0);

    /* start integration */
    electrons  = 0.0;

    rawgrid = grid_open_cmo(inforb_.nbast, cmo, work, lwork);
    npoints = 0;
    while( (blocksz=grid_getchunk_plain(rawgrid, GRID_BUFF_SZ,
                                        grid->coor, grid->weight)) >=0) {
	npoints += blocksz;
	
	for(ipnt= 0; ipnt<blocksz; ipnt++) {
	    grid->curr_point  = ipnt;
	    grid->curr_weight = grid->weight[ipnt];

	    dft_grid_getval(grid, ipnt, work, lwork);
            grid->dp.rhoa = 0;
            for (isym=0; isym<inforb_.nsym; isym++) {
               norbi=inforb_.norb[isym];
               if (norbi > 0) {
                  nbasi=inforb_.nbas[isym];
                  atvi=grid->atv + inforb_.ibas[isym];
                  movi=grid->mov + inforb_.iorb[isym];
                  cmoi=cmo + inforb_.icmo[isym];
                  FSYM(dgemv)("T",&nbasi, &norbi, &ONER,
                         cmoi, &nbasi, atvi, &ONEI, &ZEROR, movi, 
                         &ONEI);
                  nocci = inforb_.nocc[isym];
                  if (nocci > 0) {
		    grid->dp.rhoa +=
		      FSYM(ddot)(&nocci,movi, &ONEI, movi, &ONEI);
                  }
               }
            }
            grid->dp.rhob = grid->dp.rhoa;
            electrons += grid->weight[ipnt]*grid->dp.rhoa*2;
	    if(grid->dp.rhoa*2<=grid->dfthr0) continue;

	    if(grid->dogga) {
                real ngrad;
		/* transform grad vectors to molecular orbitals */
               grid->grada[0]=grid->grada[1]=grid->grada[2]=0;
               for (isym=0; isym<inforb_.nsym; isym++) {
                  norbi=inforb_.norb[isym];
                  if (norbi > 0) {
                     nbasi=inforb_.nbas[isym];
                     atvi=grid->atv + inforb_.nbast + inforb_.ibas[isym];
                     movi=grid->mov + inforb_.iorb[isym];
                     mogi=grid->mog + inforb_.iorb[isym];
                     cmoi=cmo + inforb_.icmo[isym];
                     FSYM(dgemm)("T", "N", &norbi, &THREEI, &nbasi,
                             &ONER, cmoi, &nbasi, 
                             atvi, &inforb_.nbast, &ZEROR, 
                             mogi, &inforb_.norbt);
                     nocci=inforb_.nocc[isym];
                     if (nocci > 0) {
		       FSYM(dgemv)("T", &nocci, &THREEI, &TWOR,
				   mogi, &inforb_.norbt, 
				   movi, &ONEI,  &ONER, grid->grada, &ONEI);
                     }
                  }
               }
               ngrad = sqrt(grid->grada[0]*grid->grada[0]+
                            grid->grada[1]*grid->grada[1]+
                            grid->grada[2]*grid->grada[2]); 
               grid->dp.gradb  = grid->dp.grada = ngrad;
               grid->dp.gradab = grid->dp.grada*grid->dp.gradb;
	    }
            for(cbno=0; cbno< cbcount; cbno++)
		(cbarr[cbno].callback)(grid, cbarr[cbno].cb_data);
	}
    };
    grid_close(rawgrid); 
    dft_grid_free(grid);
    dft_integrate_collect_info(&electrons); /* NO-OP for serial code */
    return electrons;
}

/* dft_integrate_ao:
   Integrate KS wave function in atomic orbital basis,
   using provided set of callbacks.
   needgrad - logical, whether callback will need first order derivatives
   (even if the functionals need only orbital values).
   needlap - logical, whether callback will need second order derivatives
   of orbitals wrt to coordinates.
   Returns integrated charge.
 */
real
dft_integrate_ao(DftDensity* dens, real* work, integer* lwork,
                 int needgrad, int needlap, int needgb, 
                 const DftCallbackData* cbarr, int cbcount)
{
    int npoints, ipnt, cbno;
    int blocksz;
    real electrons;
    DftGrid* grid;
    real *dmagao, rho;
    DftGridReader* rawgrid;
    dmagao = calloc(inforb_.nbast*2, sizeof(real));
    grid = dft_grid_new(needgrad, needlap, needgb);

    /* start integration */
    electrons  = 0.0;
    /*printf("CALLING grid_open from dft_integrate_ao\n");*/
    rawgrid = grid_open(inforb_.nbast, dens->dmata, work, lwork);
    npoints = 0;
    while( (blocksz=grid_getchunk_plain(rawgrid, GRID_BUFF_SZ,
                                        grid->coor, grid->weight)) >=0) {
	for(ipnt=0; ipnt<blocksz; ipnt++, npoints++) {
	    grid->curr_point  = ipnt;
	    grid->curr_weight = grid->weight[ipnt];
	    dft_grid_getval(grid, ipnt, work, lwork);
            dens->evaluate(dens, &grid->dp, grid, dmagao);
            rho = grid->dp.rhoa+grid->dp.rhob;
            electrons += grid->weight[ipnt]*rho;
	    if(rho<=grid->dfthr0) continue;
            for(cbno=0; cbno<cbcount; cbno++)
		(cbarr[cbno].callback)(grid, cbarr[cbno].cb_data);
	}
    };
    grid_close(rawgrid);
    free(dmagao);
    dft_grid_free(grid);
    dft_integrate_collect_info(&electrons); /* no-op for serial code */

    return electrons; 
}

/* =================================================================== */
/*                     BLOCKED INTEGRATORS                             */
/* =================================================================== */
/* Blocked integrator(s) have altered the block structure to enhance
 * data locality and increase length of internal loops. This should
 * increase performance even for small molecules and reach linear
 * scaling for large ones by enabling vector-like optimization and
 * enhancing data locality.
 */
/* dft_grid_blocked_new:
   initialize grid data.
   ndmat - number of density matrices handled at the same time 
           needed for temporary array.
   bllen - grid point batch length.
*/
DftIntegratorBl*
dft_integrator_bl_new(Functional* f, int ndmat,
                      int bllen, int needlondon)
{
    real dfthri, dfthr0;
    integer geodrv, kmax;
    DftIntegratorBl* grid = dal_malloc(sizeof(DftIntegratorBl));

    grid->coor   = calloc(3*GRID_BUFF_SZ, sizeof(real));
    grid->weight = calloc(GRID_BUFF_SZ, sizeof(real));
    grid->dogga  = f->is_gga();
    FSYM2(dft_get_thresholds)(&dfthri, &dfthr0);
    grid->dfthri = dfthri;
    grid->needlap= 0;
    grid->needgb = needlondon;
    grid->nsym   = inforb_.nsym;
    grid->ndmat  = ndmat;

    kmax = FSYM2(ishell_cnt)();
    geodrv = grid->dogga ? 1 : 0;
    FSYM(setupsos)(&geodrv, &grid->needgb, &grid->ntypso, &grid->london_off);
    grid->london_off--; /* convert from fortran offset type */
    grid->atv   = calloc(bllen*inforb_.nbast*grid->ntypso, sizeof(real));
    grid->shlblocks = dal_malloc(2*kmax*sizeof(integer));
    grid->basblocks = dal_malloc(2*kmax*8*sizeof(integer));

    /* Allocate memory for rho, taking advantage of the union. */
    grid->r.rho    = dal_malloc(ndmat*bllen*sizeof(real));
    grid->g.grad   = dal_malloc(ndmat*3*bllen*sizeof(real));
    
    /* and set some aliases in case somebody needed them for open-shell. *
     * Observe that rho aliases with rhoa. */
    if(ndmat == 2) {
        grid->r.ho.b    = grid->r.ho.a  + bllen;
        grid->g.rad.b   = grid->g.rad.a + bllen;
    }
    return grid;
}

void
dft_integrator_bl_free(DftIntegratorBl *res)
{
  free(res->coor);
  free(res->weight);
  free(res->atv);
  free(res->r.rho);
  free(res->g.grad);
  free(res->shlblocks);
  free(res->basblocks);
  free(res);
}

/* grid_blocked_getval:
   evaluates a block of bllen orbitals
*/
void FSYM(blgetsos)(integer *nvclen, real GSO[], real COOR[],
		    integer *NBLCNT, integer IBLCKS[],
		    real WORK[], integer *LWORK,
		    integer *NBAST, integer *DOLND, integer *DOGGA,
		    real *DFTHRI, const integer *IPRINT);

void __inline__
grid_blocked_getval(DftIntegratorBl* grid, int ipnt, integer bllen,
                    real* work, integer *lwork)
{
    real thrint = grid->dfthri/grid->weight[ipnt];
    FSYM(blgetsos)(&bllen, grid->atv, &grid->coor[ipnt][0],
		   &grid->shl_bl_cnt, &grid->shlblocks[0][0],
		   work, lwork, &inforb_.nbast, &grid->needgb,
		   &grid->dogga, &thrint, &ONEI);
}
extern void
FSYM(shltoorb)(const integer* shlcnt, integer (*shlblocks)[2],
	       integer aocnt[], integer (*basblocks)[2],
	       const integer *iroidx);
extern void FSYM2(construct_ioridx)(integer *ioridx);

extern void
FSYM2(getrho_blocked_lda)(real*dmat, const real* atv, 
			  const integer *bas_bl_cnt, integer (*basblocks)[2],
			  const integer *shl,
			  real *tmp, const integer *bllen, real *rho);

extern void
FSYM2(getrho_blocked_gga)(real*dmat, const real* atv, 
			  const integer *bas_bl_cnt,
			  integer (*basblocks)[2], const integer *shl,
			  real *tmp, const integer * bllen,
			  real *rho, real (*grad)[3]);

real
dft_integrate_ao_bl(int ndmat, real *dmat, real *work, integer *lwork,
                    int needlnd, DftBlockCallback cb, void *cb_data)
{
    int npoints, ipnt, i, j;
    int blocksz, lo, hi;
    integer *ioridx = dal_malloc(FSYM2(ishell_cnt)()*2*8*sizeof(real));
    real electrons; /* alpha electrons only most of the time */
    DftIntegratorBl* grid;
    real *dmagao;
    DftGridReader* rawgrid;
    dmagao = dal_malloc(inforb_.nbast*DFT_BLLEN*sizeof(real));
    grid = dft_integrator_bl_new(selected_func, ndmat,
                                 DFT_BLLEN, needlnd);
    FSYM2(construct_ioridx)(ioridx);

    /* start integration */
    electrons  = 0.0;
    /*printf("CALLING grid_open from dft_integrate_ao_bl\n");*/
    rawgrid = grid_open(inforb_.nbast, dmat, work, lwork);
    npoints = 0;
    while( (blocksz=grid_getchunk_blocked(rawgrid, GRID_BUFF_SZ,
                                          &grid->shl_bl_cnt, 
                                          &grid->shlblocks[0],
                                          &grid->coor[0],
                                          grid->weight)) >=0) {
      FSYM(shltoorb)(&grid->shl_bl_cnt, grid->shlblocks, 
		     grid->bas_bl_cnt, grid->basblocks, ioridx);

        for(ipnt=0; ipnt<blocksz; ipnt+=DFT_BLLEN) {
            integer len = ipnt+DFT_BLLEN<blocksz ? DFT_BLLEN : blocksz-ipnt;
            grid->curr_point  = ipnt;
            grid_blocked_getval(grid, ipnt, len, work, lwork);
            for(i=0; i<ndmat; i++) {
                int doff = i*inforb_.n2basx;
                int roff = i*DFT_BLLEN;
                if(grid->dogga)
                    FSYM2(getrho_blocked_gga)(dmat+doff, grid->atv, grid->bas_bl_cnt,
                                        grid->basblocks, &grid->shl_bl_cnt, 
                                        dmagao, &len,
                                        grid->r.rho  + roff, 
                                        grid->g.grad + roff);
                else
                    FSYM2(getrho_blocked_lda)(dmat+doff, grid->atv, grid->bas_bl_cnt,
                                        grid->basblocks, &grid->shl_bl_cnt,
                                        dmagao, &len, grid->r.rho+roff);
                for(j=0; j<len; j++)
                    electrons += grid->weight[ipnt+j]*grid->r.rho[j+roff];
            }
/*#if 0
            for(lo=0; lo<len && grid->rhoa[lo]<1e-10; lo++)
                ;
            for(hi=len-1; hi>lo && grid->rhoa[hi]<1e-10; hi--)
                ;
#else */
            lo = 0; hi = len-1;
/*#endif */
            npoints += len;
            if(lo<=hi)
                cb(grid, dmagao, len, lo, hi+1, cb_data);
        }
    }
    grid_close(rawgrid);
#ifdef VAR_MPI
/*
    printf("collecting e- \n");
*/
    FSYM(dftintcollect)(&electrons);
#endif
    free(dmagao);
    free(ioridx);
    dft_integrator_bl_free(grid);
    return electrons;
}
