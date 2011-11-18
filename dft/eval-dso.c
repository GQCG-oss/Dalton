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

void FSYM(getdsosz)(integer *dsodim);
void
numdso_slave(real* work, integer* lwork, const integer* iprint)
{
    integer sz, new_lwork;
    integer dsodim;
    integer nucind;            /* IN: will be synced from master */
    real *spndso;
    
    FSYM(getdsosz)(&dsodim); sz = dsodim*dsodim;
    if(*lwork < sz)
	dalton_quit("Slave has not enough memory %d for DSO matrix %d",
		    *lwork, sz);
    new_lwork = *lwork-sz;
    FSYM(dzero)(work, &sz);
    FSYM(numdso)(work, &nucind, work+sz, &new_lwork);
}
void FSYM(dsosyncslaves)(real *dmat,integer *nucind,real *work,integer *lwork);
#define numdso_sync_slaves(dmat,nucind,work,lwork) \
        FSYM(dsosyncslaves)(dmat,nucind,work,lwork)
static void
numdso_collect_info(real *spndso, real *work, integer lwork)
{
    integer dsodim, sz;
    int     sz_mpi;

    MPI_Comm_size(MPI_COMM_WORLD, &sz_mpi);
    if(sz_mpi <=1) return;

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
    extern void FSYM(dsocb)(real*, const integer*, real (*coor)[3],
                            const real*rho,const real* wght,
                            real(*rvec)[3], real*r3i);
    int off = grid->curr_point;
    integer bllen_copy = bllen; /* conversion needed for VAR_INT64 */
    FSYM(dsocb)(dso->dso, &bllen_copy,
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
FSYM(numdso)(real* spndso, integer *nucind, real* work, integer* lwork)
{
    extern void FSYM(dunfld)(const integer* n, const real* dsp, real* dge);
    struct tms starttm, endtm; clock_t utm;
    struct dso_data dso;
    real electrons, *dmat;

    dft_wake_slaves((DFTPropEvalMaster)FSYM(numdso));
    dmat = dal_malloc(inforb_.n2basx*sizeof(real));
    numdso_sync_slaves(dmat, nucind,work,lwork);

    dso.dso  = spndso;
    dso.rvec = malloc(DFT_BLLEN*(*nucind)*3*sizeof(real));
    dso.r3i  = malloc(DFT_BLLEN*(*nucind)*sizeof(real));
    if(!dso.rvec || !dso.r3i) dalton_quit("Not enough memory in numdso.");
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

