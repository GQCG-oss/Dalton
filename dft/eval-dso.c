/*
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren. 
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
*/
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

#ifdef VAR_MPI
#include <mpi.h>

void FSYM(getdsosz)(int *dsodim);
void
numdso_slave(real* work, int* lwork, const int* iprint)
{
    int dsodim,sz, new_lwork;
    int nucind;            /* IN: will be synced from master */
    real *spndso;
    
    FSYM(getdsosz)(&dsodim); sz = dsodim*dsodim;
    if(*lwork < sz)
	dalton_quit("Slave has not enough memory %d for DSO matrix %d",
		    *lwork, sz);
    new_lwork = *lwork-sz;
    FSYM(numdso)(work, &nucind, work+sz, &new_lwork);
}
#define numdso_sync_slaves(dmat,nucind,work,lwork) \
        FSYM(dsosyncslaves)(dmat,nucind,work,lwork)
static void
numdso_collect_info(real *spndso, real *work, int lwork)
{
    int dsodim,sz;

    MPI_Comm_size(MPI_COMM_WORLD, &sz);
    if(sz <=1) return;

    FSYM(getdsosz)(&dsodim); sz = dsodim*dsodim;
    dcopy_(&sz, spndso, &ONEI, work, &ONEI);
    MPI_Reduce(work, spndso, sz, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}
#else /* VAR_MPI */
#define numdso_sync_slaves(dmat,nucind,work,lwork) \
        FSYM(dftdns)(dmat, work, lwork, &ZEROI)

#define numdso_collect_info(a,work,lwork)
#endif /* VAR_MPI */

static void
dso_cb(DftIntegratorBl* grid, real * RESTRICT tmp,
       int bllen, int blstart, int blend, struct dso_data* dso)
{
    extern void FSYM(dsocb)(real*, const int*, real (*coor)[3],
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
 * provided matrix SPNDSO of size (MXCOOR,MXCOOR). Nucind is a number
 * of atoms in the molecule.
 *
 * NOTE: provided matrix is in folded, packed format. In dft-qr methodology,
 * we could just provide a packed density evaluator. Here, we do not have
 * this option and we unfold the matrix first. Not that it matters much
 * except for memory usage.
 */
extern void FSYM2(numdso_finish)(real* spndso);
void
FSYM(numdso)(real* spndso, int *nucind, real* work, int* lwork)
{
    extern void dunfld_(const int* n, const real* dsp, real* dge);
    struct tms starttm, endtm; clock_t utm;
    struct dso_data dso;
    real electrons, *dmat;

    dft_wake_slaves((DFTPropEvalMaster)FSYM(numdso));
    dmat = dal_malloc(inforb_.n2basx*sizeof(real));
    numdso_sync_slaves(dmat, nucind,work,lwork);

    dso.dso  = spndso;
    dso.rvec = malloc(DFT_BLLEN*(*nucind)*3*sizeof(real));
    dso.r3i  = malloc(DFT_BLLEN*(*nucind)*sizeof(real));
    if(!dso.rvec || !dso.r3i) dalton_quit("ABORT!!!");
    times(&starttm);
    electrons = dft_integrate_ao_bl(1, dmat, work, lwork, 0,
                                    (DftBlockCallback)dso_cb, &dso);
    numdso_collect_info(spndso, work, *lwork);

    FSYM2(numdso_finish)(spndso);
    times(&endtm);
    free(dso.rvec);
    free(dso.r3i);
    free(dmat);
    utm = endtm.tms_utime-starttm.tms_utime;
    fort_print("Electrons: %11.7f(%9.3g): DFTDSO time: %9.1f s", 
               electrons, (double)(electrons-(int)(electrons+0.5)),
	       (double)(utm/(double)sysconf(_SC_CLK_TCK)));
}

