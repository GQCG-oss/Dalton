/* this is true spare-enabled code (well, mostly). 
 * This file contains property evaluators, that is, modules that
 * setup the stage for numerical integration, call the integrator,
 * postprocess the integration results and return home happily.
 * 
 * (c) Pawel Salek, pawsa@theochem.kth.se
 */

#if !defined(SYS_DEC)
/* XOPEN compliance is missing on old Tru64 4.0E Alphas and pow() prototype
 * is not specified. */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1
#endif

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>

#define __CVERSION__
#include "general.h"
#include "functionals.h"
#include "grid-gen.h"
#include "inforb.h"
#include "integrator.h"

/* =================================================================== */
/*          DSO spin-spin coupling contribution evaluation             */
/* =================================================================== */
/*
 * the trace is given by
 * 2 (r_K . r_L)/(r_K^3 r_L^3) * alpha^4
 * where K and L are atom indices. The matrix elements are generally:
 * ((r_K . r_L)I -r_K r_L^T)/(r_K^3 r_L^3)
 * and we only compute the (xx,yy,zz) components.
 * We also leave the dirty work to fortran part since we do not want 
 * to mess with f77 common blocks. Not yet, at least.
 */

struct dso_data {
    real * dso;      /** dso(natoms,natoms) */
    real (*rvec)[3]; /** rvec(3,natoms)  - temp space */
    real * r3i;      /** r3(nvec) = 1/|rvec(:,natoms)|^3 - temp space */
};

static void
dso_cb(DftIntegratorBl* grid, real * RESTRICT tmp,
       int bllen, int blstart, int blend, struct dso_data* dso)
{
    extern void dsocb_(real*, const int*, real (*coor)[3],
                       const real*rho,const real* wght,
                       real(*rvec)[3], real*r3i);
    int off = grid->curr_point;
    FSYM(dsocb)(dso->dso, &bllen,
           grid->coor+off, grid->r.rho,
           grid->weight+off, dso->rvec, dso->r3i);
}

/** NUMDSO:
 * computes the Dipole Spin-Orbit contribution to spin-spin couplings
 * averaged with given density matrix DMAT. The result is placed in 
 * provided matrix SPNDSO of size (MXCOOR,MXCOOR).
 *
 * NOTE: provided matrix is in folded, packed format. In dft-qr methodology,
 * we could just provide a packed density evaluator. Here, we do not have
 * this option and we unfold the matrix first. Not that it matters much
 * except for memory usage.
 */
void
FSYM(numdso)(real* dmat_folded, real* spndso, int *nucind,
             real* work, int* lwork)
{
    extern void dunfld_(const int* n, const real* dsp, real* dge);
    extern void numdso_finish_(real* spndso);
    struct tms starttm, endtm; clock_t utm;
    struct dso_data dso;
    real electrons, *dmat;

    dso.dso  = spndso;
    dso.rvec = malloc(DFT_BLLEN*(*nucind)*3*sizeof(real));
    dso.r3i  = malloc(DFT_BLLEN*(*nucind)*sizeof(real));
    dmat = malloc(inforb_.n2basx*sizeof(real));
    if(!dso.rvec || !dso.r3i || !dmat) dalton_quit("ABORT!!!");

    dunfld_(&inforb_.nbast, dmat_folded, dmat);

    times(&starttm);
    electrons = dft_integrate_ao_bl(1, dmat, work, lwork, 1,
                                    (DftBlockCallback)dso_cb,&dso);

    numdso_finish_(spndso);
    times(&endtm);
    free(dso.rvec);
    free(dso.r3i);
    free(dmat);
    utm = endtm.tms_utime-starttm.tms_utime;
    fort_print("Electrons: %11.7f(%9.3g): DFTDSO time: %9.1f s", 
               electrons, (double)(electrons-(int)(electrons+0.5)),
	       (double)(utm/(double)sysconf(_SC_CLK_TCK)));
}

