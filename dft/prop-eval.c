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
/* The DFT Propery evaluators.
   (c) Pawel Salek, pawsa@theochem.kth.se.
   2002.04.05

   This file is written in C because it is easier to write portable
   programs in C than in Fortran due to more strict syntax checking and
   existing language standard that is actually obeyed by compiler writers
   (not to mention language features).

   This module evaluates DFT contribution to different properties, 
   in particular:
   a). DFT contribution to the Fock/KS matrix (dft_kohn_sham)
   b). DFT contribution to the linear response (dft_lin_resp)
   c). DFT contribution to the molecular gradient (dft_mol_grad)
*/
/* strictly conform to XOPEN ANSI C standard */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>

#define __CVERSION__
#include "general.h"
#include "integrator.h"
#include "functionals.h"

#include "inforb.h"

const static int KOHNSH_DEBUG = 0;
const static int DFTLR_DEBUG  = 0;
const static int DFTMAG_DEBUG = 0;

#if defined(VAR_MPI)
#include <mpi.h>
#include <our_extra_mpi.h>
#define MASTER_NO 0
#endif
#if 0 && defined(VAR_MPI)

/* dft_kohn_sham_slave:
   this is a slave driver. It's task is to allocate memory needed by
   the main property evaluator (dft_kohn_sham in this case) and call it.
*/
void
dft_kohn_sham_slave(real* work, integer* lwork, const integer* iprint)
{
    real* dmat = malloc(inforb_.n2basx*sizeof(real));
    real* ksm  = calloc(inforb_.n2basx,sizeof(real));
    integer iprfck = 0;
    dft_kohn_sham_(dmat, ksm, work, lwork, &iprfck);
    free(dmat);
    free(ksm);
}

static __inline__ void
dft_kohn_sham_sync_slaves(real* dmat)
{
    MPI_Bcast(dmat, inforb_.n2basx,MPI_DOUBLE,
	      MASTER_NO, MPI_COMM_WORLD);
}

static __inline__ void
dft_kohn_sham_collect_info(real*ksm, real* energy, real* work, integer lwork)
{
    real tmp = *energy;
    int sz = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &sz);
    if(sz <=1) return;
    CHECK_WRKMEM(inforb_.n2basx, lwork);
    dcopy_(&inforb_.n2basx, ksm,&ONEI, work, &ONEI);
    MPI_Reduce(work, ksm, inforb_.n2basx, MPI_DOUBLE, MPI_SUM, 
	       MASTER_NO, MPI_COMM_WORLD);
    MPI_Reduce(&tmp, energy, 1, MPI_DOUBLE, MPI_SUM, 
	       MASTER_NO, MPI_COMM_WORLD);
}

#else /* VAR_MPI */
#define dft_kohn_sham_sync_slaves(dmat)
#define dft_kohn_sham_collect_info(myksm, ksm, energy,lwork)
#endif /* VAR_MPI */

static real energy = 0.0;
static void
kohn_sham_cb(DftGrid* grid, real* excmat)
{
    FirstDrv drvs;
    int isym, i, j;
    real de = selected_func->func(&grid->dp);

    energy += de*grid->curr_weight;
    dftpot0_(&drvs, &grid->curr_weight, &grid->dp);

    if(grid->dogga) {
	real *atvX = &grid->atv[inforb_.nbast];
	real *atvY = &grid->atv[inforb_.nbast*2];
	real *atvZ = &grid->atv[inforb_.nbast*3];
	real fx, fy, fz;

	drvs.fZ *= 1.0/grid->dp.grada;
	fx = 2*drvs.fZ*grid->grada[0];
	fy = 2*drvs.fZ*grid->grada[1];
	fz = 2*drvs.fZ*grid->grada[2];

	for(isym=0; isym<inforb_.nsym; isym++) {
	    int istr = inforb_.ibas[isym];
	    int iend = inforb_.ibas[isym]+inforb_.nbas[isym];
	    for(j=istr; j<iend; j++) { 
		real fc = drvs.fR*grid->atv[j]
		    +fx*atvX[j] + fy*atvY[j] + fz*atvZ[j];
		if(fabs(fc)>grid->dfthri) {
		    int joff = j*inforb_.nbast;
		    for(i=istr; i<iend; i++)  
			excmat[i+joff] += fc*grid->atv[i];
		}
	    }
	}
    } else {
	for(isym=0; isym<inforb_.nsym; isym++) {
	    int istr = inforb_.ibas[isym];
	    int iend = inforb_.ibas[isym]+inforb_.nbas[isym];
	    for(j=istr; j<iend; j++) { 
		real gvxc = 2.0*drvs.fR*grid->atv[j];
		if(fabs(gvxc)>grid->dfthri) {
		    int joff = j*inforb_.nbast;
		    for(i=istr; i<j; i++)  
			excmat[i+joff] += gvxc*grid->atv[i];
		    excmat[j+joff] += 0.5*gvxc*grid->atv[j];
		}
	    }
	}
    }
}

/* dft_kohn_sham:
   compute Fock matrix ksm corresponding to given density matrix dmat.
*/
void
FSYM2(dft_kohn_sham)(real* dmat, real* ksm, real *edfty, 
		     real* work, integer *lwork, integer* iprfck)
{
    int nbast2, i, j;
    DftCallbackData cbdata[1];
    DftDensity dens = { dft_dens_restricted, NULL, NULL };
    struct tms starttm, endtm; clock_t utm;
    real electrons, *ksm_exch;

    /* WARNING: NO work MAY BE done before syncing slaves! */
    dft_wake_slaves((DFTPropEvalMaster)dft_kohn_sham_); /* NO-OP in serial */
    dft_kohn_sham_sync_slaves(dmat);                    /* NO-OP in serial */

    dens.dmata = dmat;
    nbast2   = inforb_.nbast*inforb_.nbast;
    ksm_exch =  calloc(nbast2, sizeof(real));
    cbdata[0].callback = (DftCallback)kohn_sham_cb;
    cbdata[0].cb_data  = ksm_exch;

    times(&starttm);
    energy = 0.0;

    electrons = dft_integrate_ao(&dens, work, lwork, 0, 0,0, 
				 cbdata, ELEMENTS(cbdata));

    for(i=0; i<inforb_.nbast; i++) {
	int ioff = i*inforb_.nbast;
	for(j=0; j<i; j++) {
	    int joff = j*inforb_.nbast;
	    real averag = 0.5*(ksm_exch[i+joff] + ksm_exch[j+ioff]);
	    ksm_exch[i+joff] = ksm_exch[j+ioff] = averag;
	}
    }
    dft_kohn_sham_collect_info(ksm_exch, &energy, work, *lwork);
    *edfty = energy;
    daxpy_(&inforb_.n2basx, &ONER, ksm_exch, &ONEI, ksm, &ONEI);
    
    free(ksm_exch);
    times(&endtm);
    utm = endtm.tms_utime-starttm.tms_utime;
    fort_print("      Electrons: %11.7f %8.2g: Energy %12.6f KS/B time: %9.1f s", 
               electrons, (electrons-2.0*inforb_.nrhft)/(2.0*inforb_.nrhft), 
               energy, utm/(double)sysconf(_SC_CLK_TCK));
}

/* ------------------------------------------------------------------- */
/* ---------- DFT LINEAR RESPONSE CONTRIBUTION EVALUATOR ------------- */
/* ------------------------------------------------------------------- */
void FSYM(deq27)(const real* cmo, const real* ubo, const real* dv, 
                 real* dxcao, real* dxvao, real* wrk, integer* lfrsav);
void FSYM(autpv)(const integer*isym,const integer*jsym,
		 const real*u,const real* v, 
                 real*prpao, const integer*nbas,const integer*nbast, real*prpmo,
                 const integer*norb, const integer*norbt,
		 real* wrk,integer* lwrk);

typedef struct {
    real* dmat, *kappa, *res;
    real* dtgao;
    integer   trplet, ksymop;
} LinRespData;


static __inline__ real
min(real a, real b)
{ return a>b ? b : a; }
static __inline__ real
max(real a, real b)
{ return a>b ? a : b; }

/* dft_lin_resp_slave:
   this is a slave driver. It's task is to allocate memory needed by
   the main property evaluator (dft_lin_resp in this case) and call
   it.  The computed values are collected in the evaluater using
   dft_lin_resp_collect info and we can just drop the arrays on the
   floor when we are done.
   FIXME: we can use work for the temp data....
*/
#if defined(VAR_MPI)

extern void FSYM(dftlrsync)(void);
static void
dft_lin_resp_sync_slaves(real* cmo, integer *nvec, real**zymat,
                         integer* trplet, integer* ksymop,
                         real **work, integer lwork)
{
    static SyncData data2[] = 
    { {NULL, 0, MPI_DOUBLE}, {NULL, 0, MPI_DOUBLE}, {NULL, 0, fortran_MPI_INT},
      {NULL, 0, fortran_MPI_INT} };
    MPI_Bcast(nvec, 1, fortran_MPI_INT, MASTER_NO, MPI_COMM_WORLD);
    FSYM(dftlrsync)();
    if(*zymat == NULL) { /* we are a slave */
        if(lwork<*nvec*inforb_.n2orbx)
            dalton_quit("%s:  slave needs %d words for ZYMAT, available %d",
                        __FUNCTION__, *nvec*inforb_.n2orbx, lwork);
	*zymat = *work;
        *work += *nvec*inforb_.n2orbx;
    }
    data2[0].data = cmo;    data2[0].count = inforb_.ncmot;
    data2[1].data = *zymat; data2[1].count = *nvec*inforb_.n2orbx;
    data2[2].data = trplet; data2[2].count = 1;
    data2[3].data = ksymop; data2[3].count = 1;
    mpi_sync_data(data2, ELEMENTS(data2));
}

static __inline__ void
dft_lin_resp_collect_info(integer nvec, real* fmat, real*work, integer lwork)
{
    int sz_mpi = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &sz_mpi);
    if(sz_mpi<=1) return;

    integer sz = nvec*inforb_.n2basx;
    sz_mpi = sz;
    if(sz>lwork) {
        CHECK_WRKMEM(inforb_.n2basx,lwork);
        for(sz=0; sz<nvec; sz++) {
            dcopy_(&inforb_.n2basx, fmat + sz*inforb_.n2basx,
                   &ONEI, work, &ONEI);
            int len_for_reduce = inforb_.n2basx;
            MPI_Reduce(work, fmat+sz*inforb_.n2basx, len_for_reduce,
                       MPI_DOUBLE, MPI_SUM, MASTER_NO, MPI_COMM_WORLD);
        }
    } else {
        CHECK_WRKMEM(sz, lwork);
        dcopy_(&sz, fmat, &ONEI, work, &ONEI);
        MPI_Reduce(work, fmat, sz_mpi, MPI_DOUBLE, MPI_SUM, 
                   MASTER_NO, MPI_COMM_WORLD);
    }
}

#else  /* VAR_MPI */
#define dft_lin_resp_sync_slaves(cmo,nvec,zymat,trplet,ksymop,work,lwork)
#define dft_lin_resp_collect_info(nvec,fmat,work,lwork)
#endif /* VAR_MPI */

static void
lin_resp_cb(DftGrid* grid, LinRespData* data)
{
    real b0, b3[3];
    int isym, i, j;
    SecondDrv vxc;

    dgemv_("N",&inforb_.nbast,&inforb_.nbast,&ONER,
	   data->kappa,&inforb_.nbast,grid->atv,&ONEI,&ZEROR,
	   data->dtgao,&ONEI);
    b0 = ddot_(&inforb_.nbast,data->dtgao,&ONEI,grid->atv,&ONEI);
    if(grid->dogga) {
	real brg, brz, facr, facg, bmax, vt[4];
	real *atvX = &grid->atv[inforb_.nbast  ];
	real *atvY = &grid->atv[inforb_.nbast*2];
	real *atvZ = &grid->atv[inforb_.nbast*3];
        real ngrad = (grid->dp.grada + grid->dp.gradb);
	/* B3 = GAO1'*DTGAO; */
	dgemv_("T",&inforb_.nbast,&THREEI,&ONER,atvX, &inforb_.nbast,
	       data->dtgao,&ONEI, &ZEROR, b3, &ONEI);
	/* DTGAO= DTRMAT'*GAO */
	dgemv_("T",&inforb_.nbast, &inforb_.nbast, &ONER, data->kappa,
	       &inforb_.nbast,grid->atv,&ONEI,&ZEROR,data->dtgao,&ONEI);
	/*  B3 = B3 + GAO1'*DTGAO */
	dgemv_("T",&inforb_.nbast,&THREEI,&ONER,atvX,
	       &inforb_.nbast,data->dtgao,&ONEI,&ONER,b3,&ONEI);
	bmax = max(fabs(b0),max(fabs(b3[0]),max(fabs(b3[1]),fabs(b3[2]))));
	if(bmax<=grid->dfthri) return;
	brg = (b3[0]*grid->grada[0] +
               b3[1]*grid->grada[1] +
               b3[2]*grid->grada[2])*2;
        brz = brg/ngrad;

	dftpot1_(&vxc, &grid->curr_weight, &grid->dp, &data->trplet);
        facr = vxc.fRZ*b0 + (vxc.fZZ-vxc.fZ/ngrad)*brz + vxc.fZG*brg;
        facr = 2*(facr/ngrad + (vxc.fRG*b0+vxc.fZG*brz +vxc.fGG*brg));
        facg = vxc.fZ/ngrad + vxc.fG;
        vt[0] = vxc.fRR*b0 + vxc.fRZ*brz+ vxc.fRG*brg;
        vt[1] = facr*grid->grada[0] + facg*b3[0];
        vt[2] = facr*grid->grada[1] + facg*b3[1];
        vt[3] = facr*grid->grada[2] + facg*b3[2];
	for(isym=0; isym<inforb_.nsym; isym++) {
	    int istr = inforb_.ibas[isym];
	    int iend = inforb_.ibas[isym] + inforb_.nbas[isym];
	    int jsym = inforb_.muld2h[data->ksymop-1][isym]-1;
	    if (isym>=jsym) {
		int jstr = inforb_.ibas[jsym];
		int jend = inforb_.ibas[jsym] + inforb_.nbas[jsym];
		for(i=istr; i<iend; i++) {
		    real g0 = grid->atv[i];
		    real gx = atvX[i];
		    real gy = atvY[i];
		    real gz = atvZ[i];
		    int ioff = i*inforb_.nbast;
                    int jtop = min(jend,i+1);
		    for(j=jstr; j<jtop; j++) {
			real a0 = g0*grid->atv[j];
			real ax = gx*grid->atv[j] + g0*atvX[j];
			real ay = gy*grid->atv[j] + g0*atvY[j];
			real az = gz*grid->atv[j] + g0*atvZ[j];
			data->res[j+ioff] +=  
                            vt[0]*a0 + vt[1]*ax +
                            vt[2]*ay + vt[3]*az;
		    }
		}
	    }
	}
    } else { /* dogga == FALSE */
	real vt;
	if(fabs(b0)<=grid->dfthri) return;
	dftpot1_(&vxc, &grid->curr_weight, &grid->dp, &data->trplet);
	vt = vxc.fRR*b0;
	for(isym=0; isym<inforb_.nsym; isym++) {
	    int jsym = inforb_.muld2h[data->ksymop-1][isym]-1;
	    if(isym>=jsym) {
		int istr = inforb_.ibas[isym];
		int iend = inforb_.ibas[isym] + inforb_.nbas[isym];
		int jstr = inforb_.ibas[jsym];
		int jend = inforb_.ibas[jsym] + inforb_.nbas[jsym]-1;
		for(i=istr; i<iend; i++) {
		    int ioff = i*inforb_.nbast;
		    real gvi = vt*grid->atv[i];
		    for(j=min(i,jend); j>=jstr; j--)
			data->res[j+ioff] += gvi*grid->atv[j];
		}
	    }
	}
    }
}

/* dft_lin_resp_:
   Main Linear Response evaluator.
   FMAT(NORBT,NORBT) - result added to FMAT.
   cmo -const
   zymat(NORBT,NORBT) - const response vector.
   trplet - triplet excitation? (bool)
*/
void
FSYM2(dft_lin_resp)(real* fmat, real *cmo, real *zymat, integer *trplet,
		    integer *ksymop, real* work, integer* lwork)
{
    struct tms starttm, endtm; clock_t utm;
    real electrons;
    LinRespData lr_data; /* linear response data */
    DftCallbackData cbdata[1];
    DftDensity dens = { dft_dens_restricted, NULL, NULL };
    real dummy;
    integer nvec1=1;
    int i, j;
    
    /* WARNING: NO work MAY BE done before syncing slaves! */
    dft_wake_slaves((DFTPropEvalMaster)dft_lin_resp_);    /* NO-OP in serial */
    dft_lin_resp_sync_slaves(cmo,&nvec1,&zymat,trplet,ksymop,
                             &work, *lwork);              /* NO-OP in serial */

    times(&starttm);
    lr_data.dmat  = malloc(inforb_.n2basx*sizeof(real));
    lr_data.res   = calloc(inforb_.n2basx,sizeof(real));
    lr_data.kappa = calloc(inforb_.n2basx,sizeof(real));
    lr_data.dtgao = malloc(inforb_.nbast *sizeof(real));
    lr_data.trplet= *trplet;
    lr_data.ksymop = *ksymop;
    dens.dmata = lr_data.dmat;

    FSYM2(dft_get_ao_dens_mat)(cmo, lr_data.dmat, work, lwork);
    deq27_(cmo,zymat,&dummy,lr_data.kappa,&dummy,work,lwork);
    if(DFTLR_DEBUG) {
        fort_print("kappa matrix in dft_lin_resp");
        outmat_(lr_data.kappa,&ONEI,&inforb_.nbast,&ONEI,&inforb_.nbast,
                &inforb_.nbast, &inforb_.nbast);
    }

    cbdata[0].callback = (DftCallback)lin_resp_cb;
    cbdata[0].cb_data  = &lr_data;    
    electrons = dft_integrate_ao(&dens, work, lwork, 0, 0, 0,
				 cbdata, ELEMENTS(cbdata));

    if(DFTLR_DEBUG) {
        fort_print("AO Fock matrix contribution in dft_lin_resp");
        outmat_(lr_data.res,&ONEI,&inforb_.nbast,&ONEI,&inforb_.nbast,
                &inforb_.nbast, &inforb_.nbast);
    }

    for(i=0; i<inforb_.nbast; i++) {
	int ioff = i*inforb_.nbast;
	for(j=0; j<i; j++) {
	    int joff = j*inforb_.nbast;
	    real averag = lr_data.res[i+joff] + lr_data.res[j+ioff];
	    lr_data.res[i+joff] = lr_data.res[j+ioff] = averag;
	}
    }

    if(DFTLR_DEBUG) {
        fort_print("MO Fock matrix contribution in dft_lin_resp");
        outmat_(work,&ONEI,&inforb_.norbt,&ONEI,&inforb_.norbt,
                &inforb_.norbt, &inforb_.norbt);
    }
    
    /* transform to MO Fock matrix contribution  */
    FSYM(lrao2mo)(cmo, ksymop, lr_data.res, fmat, work, lwork);

    /* FIXME: parallelization here! */
    free(lr_data.dmat);
    free(lr_data.kappa);
    free(lr_data.res);
    free(lr_data.dtgao);
    times(&endtm);
    utm = endtm.tms_utime-starttm.tms_utime;
    fort_print("      Electrons: %f(%9.1g): LR-DFT eval. time: %9.1f s", 
               electrons, (electrons-2.0*inforb_.nrhft)/(2.0*inforb_.nrhft), 
               utm/(double)sysconf(_SC_CLK_TCK));
}

/* ------------------------------------------------------------------- */
/* ---------- DFT LONDON ORBITAL TRANSFORMATION EVALUATOR ------------ */
/* ------------------------------------------------------------------- */
void
FSYM(dftmag)(real* excmat,const real* x,const real* y,const real* z,
	     const real* gao, const real* gao1,const real* gab1,
	     const real* gab2,const real* vxc, const real* vxb,
	     const real* rh, const integer* dogga,const integer* fromvx);

static void
london_cb(DftGrid* grid, real* d)
{
    FirstDrv drvs;
    int lnd_off = grid->london_off*inforb_.nbast;

    dftpot0_(&drvs, &grid->curr_weight, &grid->dp);
    if(grid->dogga) drvs.fZ *= 1.0/grid->dp.grada;
    /*
    fort_print("DFTMAG: (%9.4f,%9.4f,%9.4f): RHO=%g %g %g", 
               grid->corx[grid->curr_point],
               grid->cory[grid->curr_point],
               grid->corz[grid->curr_point],
               grid->rho, 
               drvs.df1000, drvs.df0010);
    */
    dftmag_(d,&grid->coor[grid->curr_point][0],
	    &grid->coor[grid->curr_point][1],
	    &grid->coor[grid->curr_point][2],
	    grid->atv, &grid->atv[inforb_.nbast],
	    &grid->atv[lnd_off],&grid->atv[lnd_off+inforb_.nbast],
	    &drvs.fR, &drvs.fZ, grid->grada, &grid->dogga,
	    &ZEROI);
}

void
FSYM(stopit)(const char* sub,const char* place,
	     const integer* int1, const integer* int2,
	     int sub_len, int place_len);

void
FSYM2(dft_london)(real* fx, real* fy, real* fz, real* work,
		  integer* lwork, integer* iprint)
{
    DftCallbackData cbdata[1];
    DftDensity dens = { dft_dens_restricted, NULL, NULL };
    struct tms starttm, endtm; clock_t utm;
    real electrons;
    int i, j;
    real* new_work, *xcmat, *dmat;
    integer new_lwork, xc_sz=3*inforb_.n2basx;

    xcmat    = work;
    dmat     = work + 3*inforb_.n2basx;
    new_work = work + 4*inforb_.n2basx;
    new_lwork = *lwork-4*inforb_.n2basx;
    if (new_lwork<0) stopit_("DFTLND"," ", lwork, lwork, 6, 1);
    cbdata[0].callback = (DftCallback)london_cb;
    cbdata[0].cb_data  = xcmat;
    dens.dmata = dmat;

    times(&starttm);
    FSYM(dftdns)(dmat, work, &new_lwork, iprint);
    dzero_(xcmat, &xc_sz);
    electrons = dft_integrate_ao(&dens, new_work, &new_lwork, 0, 0, 1,
				 cbdata, ELEMENTS(cbdata));
    times(&endtm);
    /* dft_mol_grad_collect_info(work);                  NO-OP in serial */
    utm = endtm.tms_utime-starttm.tms_utime;
    fort_print("      Electrons: %f (%9.1g): LONDON evaluation time: %10.2f s", 
               electrons, (double)(electrons-(int)(electrons+0.5)),
	       (double)(utm/(double)sysconf(_SC_CLK_TCK)));

    /* post-transformation: dftdrv stage */
    /*output_(xcmat,&ONEI,&inforb_.nbast,&ONEI,&inforb_.nbast,
      &inforb_.nbast, &inforb_.nbast, &ONEI, &priunt_.lupri); */

    for(i=0; i<inforb_.nbast; i++)
	for(j=0; j<i; j++) {
            real averag = 0.5*(xcmat[i+j*inforb_.nbast]-
			       xcmat[j+i*inforb_.nbast]);
            xcmat[i+j*inforb_.nbast] =  averag;
            xcmat[j+i*inforb_.nbast] = -averag;
	}

    if(DFTMAG_DEBUG) {
        fort_print("DFT LDN contribution (x direction)");
        outmat_(xcmat,&ONEI,&inforb_.nbast,&ONEI,&inforb_.nbast,
                &inforb_.nbast, &inforb_.nbast);
        fort_print("DFT LDN contribution (y direction)");
        outmat_(xcmat+inforb_.n2basx,&ONEI,&inforb_.nbast,&ONEI,&inforb_.nbast,
                &inforb_.nbast, &inforb_.nbast);
        fort_print("DFT LDN contribution (z direction)");
        outmat_(xcmat+2*inforb_.n2basx,&ONEI,&inforb_.nbast,&ONEI,
                &inforb_.nbast, &inforb_.nbast, &inforb_.nbast);
    }

    /* post-transformation: dftlnd stage */
    daxpy_(&inforb_.n2basx, &ONER,xcmat,                 &ONEI,fx,&ONEI);
    daxpy_(&inforb_.n2basx, &ONER,xcmat+inforb_.n2basx,  &ONEI,fy,&ONEI);
    daxpy_(&inforb_.n2basx, &ONER,xcmat+inforb_.n2basx*2,&ONEI,fz,&ONEI);
    if(DFTMAG_DEBUG) {
        fort_print("FX (x direction)");
        outmat_(fx,&ONEI,&inforb_.nbast,&ONEI,&inforb_.nbast,
                &inforb_.nbast, &inforb_.nbast);
        fort_print("FY (y direction)");
        outmat_(fy,&ONEI,&inforb_.nbast,&ONEI,&inforb_.nbast,
                &inforb_.nbast, &inforb_.nbast);
        fort_print("FZ (z direction)");
        outmat_(fz,&ONEI,&inforb_.nbast,&ONEI,&inforb_.nbast,
                &inforb_.nbast, &inforb_.nbast);
    }

}


/* ================================================================== */
/* Open shell versions of the dft property evaluators                 */
/* ================================================================== */
/* ZR alpha/beta version */
typedef struct {
    real *ksma, *ksmb;
    real energy;
} DftKohnShamU;

#if defined(VAR_MPI)
/* dft_kohn_sham_ab_slave:
   this is a slave driver. It's task is to allocate memory needed by
   the main property evaluator (dft_kohn_shamab in this case) and call it.
*/
void
dft_kohn_shamab_slave(real* work, integer* lwork, const integer* iprint)
{
    real* dmat = malloc(2*inforb_.n2basx*sizeof(real));
    real* ksm  = calloc(2*inforb_.n2basx,sizeof(real));
    integer iprfck = 0;
    real edfty;
    FSYM2(dft_kohn_shamab)(dmat, ksm, &edfty, work, lwork, &iprfck);
    free(dmat);
    free(ksm);
}

static __inline__ void
dft_kohn_shamab_sync_slaves(real* dmat)
{
    MPI_Bcast(dmat,2*inforb_.n2basx,MPI_DOUBLE,
	      MASTER_NO, MPI_COMM_WORLD);
}

static __inline__ void
dft_kohn_shamab_collect_info(real*ksm, real* energy, real* work, integer lwork)
{
    real tmp = *energy;
    int sz = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &sz);
    if(sz<=1) return;

    sz = 2*inforb_.n2basx;
    integer sz_ftn = sz;
    CHECK_WRKMEM(sz,lwork);
    dcopy_(&sz_ftn, ksm, &ONEI, work, &ONEI);
    MPI_Reduce(work, ksm, sz, MPI_DOUBLE, MPI_SUM, 
	       MASTER_NO, MPI_COMM_WORLD);
    MPI_Reduce(&tmp, energy, 1, MPI_DOUBLE, MPI_SUM, 
	       MASTER_NO, MPI_COMM_WORLD);
}

#else /* VAR_MPI */
#define dft_kohn_shamab_sync_slaves(dmat)
#define dft_kohn_shamab_collect_info(myksm, ksm, energy,lwork)
#endif /* VAR_MPI */

static void
kohn_shamab_cb(DftGrid* grid, DftKohnShamU* exc)
{
    FunFirstFuncDrv drvs;
    int isym, i, j;

    drv1_clear(&drvs);
    exc->energy += selected_func->func(&grid->dp)*grid->curr_weight;
    if(grid->dp.rhob<1e-20) grid->dp.rhob = 1e-20;
    if(grid->dp.gradb<1e-20) grid->dp.gradb = 1e-20;
    if(grid->dp.gradab<1e-20) /* so small that we can ignore direction */
      grid->dp.gradab = 1e-20;
    selected_func->first(&drvs, grid->curr_weight, &grid->dp);
    if(isnan(drvs.df1000) || isnan(drvs.df0100) || isinf(drvs.df0100) ||
       isnan(drvs.df0010) || isnan(drvs.df0001) ||
       isnan(drvs.df00001)) {
      fort_print("AAA: %g %g %g %g %g => %g %g %g %g %g",
                 grid->dp.rhoa, grid->dp.grada,
                 grid->dp.rhob, grid->dp.gradb, grid->dp.gradab,
                 drvs.df1000, drvs.df0100, drvs.df0010, drvs.df0001,
                 drvs.df00001);
    }
    if(grid->dogga) {
	real *atvX = &grid->atv[inforb_.nbast];
	real *atvY = &grid->atv[inforb_.nbast*2];
	real *atvZ = &grid->atv[inforb_.nbast*3];
	real fxa, fya, fza, fxb, fyb, fzb;
        /* add epsilon to avoid division by zero */
        real grada = fabs(grid->dp.grada)>1e-40 ? grid->dp.grada : 1e-40;
        real gradb = fabs(grid->dp.gradb)>1e-40 ? grid->dp.gradb : 1e-40;
        /* alpha  */ 
	drvs.df0010 *= 2.0/grada;
	fxa = drvs.df0010*grid->grada[0]+2.0*drvs.df00001*grid->gradb[0];
	fya = drvs.df0010*grid->grada[1]+2.0*drvs.df00001*grid->gradb[1];
	fza = drvs.df0010*grid->grada[2]+2.0*drvs.df00001*grid->gradb[2];
        /*beta */
	drvs.df0001 *= 2.0/gradb;
	fxb = drvs.df0001*grid->gradb[0]+2.0*drvs.df00001*grid->grada[0];
	fyb = drvs.df0001*grid->gradb[1]+2.0*drvs.df00001*grid->grada[1];
	fzb = drvs.df0001*grid->gradb[2]+2.0*drvs.df00001*grid->grada[2];

	for(isym=0; isym<inforb_.nsym; isym++) {
	    int istr = inforb_.ibas[isym];
	    int iend = inforb_.ibas[isym]+inforb_.nbas[isym];
	    for(j=istr; j<iend; j++) { 
		real fca = drvs.df1000*grid->atv[j]
		    +fxa*atvX[j] + fya*atvY[j] + fza*atvZ[j];
                real fcb = drvs.df0100*grid->atv[j]
		    +fxb*atvX[j] + fyb*atvY[j] + fzb*atvZ[j];
                int joff = j*inforb_.nbast;
                for(i=istr; i<iend; i++)  
                    exc->ksma[i+joff] += fca*grid->atv[i];
                for(i=istr; i<iend; i++)  
                    exc->ksmb[i+joff] += fcb*grid->atv[i];
            }
	}
    } else {   
        for(isym=0; isym<inforb_.nsym; isym++) {
            int istr = inforb_.ibas[isym];
            int iend = inforb_.ibas[isym]+inforb_.nbas[isym];
            for(j=istr; j<iend; j++) { 
                real gvxca = 2.0*drvs.df1000*grid->atv[j];
                real gvxcb = 2.0*drvs.df0100*grid->atv[j];
                int joff = j*inforb_.nbast;
                for(i=istr; i<j; i++)  
                    exc->ksma[i+joff] += gvxca*grid->atv[i];
                exc->ksma[j+joff] += 0.5*gvxca*grid->atv[j];
                for(i=istr; i<j; i++)  
                    exc->ksmb[i+joff] += gvxcb*grid->atv[i];
                exc->ksmb[j+joff] += 0.5*gvxcb*grid->atv[j];
            }
        }
    } 
}


/* ZR alpha/beta version:
 * we require that dmata and dmatb are aligned consecutively in memory
 * in order to reduce the MPI latency penalty.
 */
void
FSYM2(dft_kohn_shamab)(real* dmat, real* ksm, real *edfty,
		       real* work, integer *lwork, integer* iprfck)
{
    int nbast2, i, j;
    DftCallbackData cbdata[1];
    DftKohnShamU res; /* res like result */
    DftDensity dens = { dft_dens_unrestricted, NULL, NULL };
    struct tms starttm, endtm; clock_t utm;
    real electrons, exp_el = 2.0*inforb_.nrhft+inforb_.nasht;
    integer sz;

    /* WARNING: NO work MAY BE done before syncing slaves! */
    dft_wake_slaves((DFTPropEvalMaster)FSYM2(dft_kohn_shamab));/* NO-OP in serial */
    dft_kohn_shamab_sync_slaves(dmat);                   /* NO-OP in serial */

    dens.dmata = dmat; dens.dmatb = dmat + inforb_.n2basx;
    nbast2   = inforb_.nbast*inforb_.nbast;
    res.ksma = calloc(2*nbast2, sizeof(real));
    res.ksmb = res.ksma + nbast2;
    cbdata[0].callback = (DftCallback)kohn_shamab_cb;
    cbdata[0].cb_data  = &res;

    times(&starttm);
    res.energy = 0.0;
    electrons = dft_integrate_ao(&dens, work, lwork, 0, 0,0,
                                 cbdata, ELEMENTS(cbdata));
    dft_kohn_shamab_collect_info(res.ksma, &res.energy, work, *lwork);

    for(i=0; i<inforb_.nbast; i++) {
	int ioff = i*inforb_.nbast;
	for(j=0; j<i; j++) {
	    int joff = j*inforb_.nbast;
	    real averag = 0.5*(res.ksma[i+joff] + res.ksma[j+ioff]);
	    res.ksma[i+joff] = res.ksma[j+ioff] = averag;
            averag = 0.5*(res.ksmb[i+joff] + res.ksmb[j+ioff]);    
            res.ksmb[i+joff] = res.ksmb[j+ioff] = averag; 
	}
    }
    *edfty=res.energy;
    sz = 2*inforb_.n2basx;
    daxpy_(&sz, &ONER, res.ksma, &ONEI, ksm, &ONEI);

    free(res.ksma);
    times(&endtm);
    utm = endtm.tms_utime-starttm.tms_utime;
    fort_print("      Electrons: %11.7f %9.1g: Energy %f KS time: %9.1f s", 
               electrons, (electrons-exp_el)/exp_el,
               res.energy, utm/(double)sysconf(_SC_CLK_TCK));
}

/* =================================================================== */
typedef struct {
    real* dmata, *dmatb, *kappaa, *kappab, *resa, *resb;
    real* dtgaoa, *dtgaob;
    integer trplet, ksymop;
} LinRespDataab;


#if defined(VAR_MPI)
void
dft_lin_respab_slave(real* work, integer* lwork, const integer* iprint)
{
    real *fmat = calloc(2*inforb_.n2orbx,sizeof(real));            /* OUT */
    real *cmo  = malloc(inforb_.norbt*inforb_.nbast*sizeof(real)); /* IN  */
    real *zymat= malloc(2*inforb_.n2orbx*sizeof(real));            /* IN  */
    integer trplet;                         /* IN: will be synced from master */
    integer ksymop;                         /* IN: will be synced from master */
    FSYM2(dft_lin_respab)(fmat, fmat+inforb_.n2orbx, cmo, zymat, &trplet,
			  &ksymop, work, lwork);
    free(fmat);
    free(cmo);
    free(zymat); 
}

static __inline__ void
dft_lin_respab_collect_info(real* fmat, real*work, integer lwork)
{
    int sz = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &sz);
    if(sz<=1) return;

    sz = 2*inforb_.n2orbx; 
    integer sz_ftn =sz;
    CHECK_WRKMEM(sz,lwork);
    dcopy_(&sz_ftn, fmat, &ONEI, work, &ONEI);
    MPI_Reduce(work, fmat, sz, MPI_DOUBLE, MPI_SUM, 
	       MASTER_NO, MPI_COMM_WORLD);
}

#else  /* VAR_MPI */
#define dft_lin_respab_collect_info(fmat,work, lwork)
#endif /* VAR_MPI */


static void
compute_trans_rho(DftGrid* grid, LinRespDataab* data, real* rhowa, real* rhowb)
{
    dgemv_("N",&inforb_.nbast,&inforb_.nbast,&ONER,
	   data->kappaa,&inforb_.nbast,grid->atv,&ONEI,&ZEROR,
	   data->dtgaoa,&ONEI);
    *rhowa = ddot_(&inforb_.nbast,data->dtgaoa,&ONEI,grid->atv,&ONEI);
    dgemv_("N",&inforb_.nbast,&inforb_.nbast,&ONER,
	   data->kappab,&inforb_.nbast,grid->atv,&ONEI,&ZEROR,
	   data->dtgaob,&ONEI);
    *rhowb = ddot_(&inforb_.nbast,data->dtgaob,&ONEI,grid->atv,&ONEI);
}


/* lin_resp_cbab:
 * ZR alpa/beta version of the linear response evaluation routine.
 * triplet: hypotetical case of triplet excitations
 *          should work for general case....       */   
static void
lin_resp_cbab_gga(DftGrid* grid, LinRespDataab* data)
{
    /*    real b0a, b0b, b3a[3],b3b[3],gradab; */
    real rhowa,rhowb, gradwa[3], gradwb[3];  
    int isym, i, j;
    int jsym, istr, iend, jstr, jend, ioff;
    FunSecondFuncDrv vxc;    
    real zetaa, zetab, zetac; 
    real znva, rxa, rya, rza; 
    real znvb, rxb, ryb, rzb;
    real fac0, facr, facz;
    real ar, ab, arb, abw;
    real g0, gx, gy, gz;
    real a0, ax, ay, az;
    real *atvX = &grid->atv[inforb_.nbast];
    real *atvY = &grid->atv[inforb_.nbast*2];
    real *atvZ = &grid->atv[inforb_.nbast*3];

    compute_trans_rho(grid, data, &rhowa, &rhowb);
     if (data->trplet) rhowb = -rhowb; 
    /* gradwa evaluation */
    dgemv_("T",&inforb_.nbast,&THREEI,&ONER,atvX, &inforb_.nbast,
           data->dtgaoa,&ONEI, &ZEROR, gradwa, &ONEI);
    dgemv_("T",&inforb_.nbast, &inforb_.nbast, &ONER, data->kappaa,
           &inforb_.nbast,grid->atv,&ONEI,&ZEROR,data->dtgaoa,&ONEI);
    dgemv_("T",&inforb_.nbast,&THREEI,&ONER,atvX,
           &inforb_.nbast,data->dtgaoa,&ONEI,&ONER,gradwa,&ONEI);
    /* gradwb evaluation */ 
    dgemv_("T",&inforb_.nbast,&THREEI,&ONER,atvX, &inforb_.nbast,
           data->dtgaob,&ONEI, &ZEROR, gradwb, &ONEI);
    dgemv_("T",&inforb_.nbast, &inforb_.nbast, &ONER, data->kappab,
           &inforb_.nbast,grid->atv,&ONEI,&ZEROR,data->dtgaob,&ONEI);
    dgemv_("T",&inforb_.nbast,&THREEI,&ONER,atvX,
           &inforb_.nbast,data->dtgaob,&ONEI,&ONER,gradwb,&ONEI); 
    if  (data->trplet) {
        gradwb[0]= -gradwb[0];
        gradwb[1]= -gradwb[1];
        gradwb[2]= -gradwb[2];
    }
    /* second derivatives calculation */
    drv2_clear(&vxc);
    selected_func->second(&vxc, grid->curr_weight, &grid->dp);
    /* alpha coeficients */
    /* add epsilon to avoid division by zero */
    znva = 1.0/(fabs(grid->dp.grada)>1e-40 ? grid->dp.grada : 1e-40);
    vxc.df0010 = znva*vxc.df0010; 
    rxa = znva*grid->grada[0];
    rya = znva*grid->grada[1];
    rza = znva*grid->grada[2];
    /* beta coeficients */
    znvb = 1.0/(fabs(grid->dp.gradb)>1e-40 ? grid->dp.gradb : 1e-40);
    vxc.df0001 = znvb*vxc.df0001; 
    rxb = znvb*grid->gradb[0];
    ryb = znvb*grid->gradb[1];
    rzb = znvb*grid->gradb[2]; 
    /* variations of functionals variables */
    zetaa = gradwa[0]*rxa + gradwa[1]*rya + gradwa[2]*rza;
    zetab = gradwb[0]*rxb + gradwb[1]*ryb + gradwb[2]*rzb;
    zetac = gradwa[0]*grid->gradb[0]+gradwa[1]*grid->gradb[1] + 
        gradwa[2]*grid->gradb[2]+gradwb[0]*grid->grada[0] +
        gradwb[1]*grid->grada[1]+gradwb[2]*grid->grada[2];  

    /* Derivatives can go in principle to infinity if evaluated at
     * inapriopriate points. Take precautions. */
    fac0 = vxc.df1010*zetaa+vxc.df1001*zetab +vxc.df10001*zetac;
    if(fabs(rhowa)>0) fac0 += vxc.df2000*rhowa;
    if(fabs(rhowa)>0) fac0 += vxc.df1100*rhowb;
    facr = vxc.df1010*rhowa+vxc.df0110*rhowb+vxc.df0020*zetaa
        +vxc.df0011*zetab
        +vxc.df00101*zetac;
    /* note: the mixed zetac derivatives not included in facr */ 
    facz = vxc.df10001*rhowa+vxc.df01001*rhowb
         + vxc.df00101*zetaa+vxc.df00011*zetab+vxc.df00002*zetac;
    /* note: the mixed zetac derivatives not included in facz */     
    for(isym=0; isym<inforb_.nsym; isym++) {
        istr = inforb_.ibas[isym];
        iend = inforb_.ibas[isym] + inforb_.nbas[isym];
        jsym = inforb_.muld2h[data->ksymop-1][isym]-1;
        if (isym>=jsym) {
            jstr = inforb_.ibas[jsym];
            jend = inforb_.ibas[jsym] + inforb_.nbas[jsym];
            for(i=istr; i<iend; i++) {
                g0 = grid->atv[i];
                gx = atvX[i];
                gy = atvY[i];
                gz = atvZ[i];
                ioff = i*inforb_.nbast;
                for(j=min(i,jend); j>=jstr; j--) {
                    a0 = g0*grid->atv[j];
                    ax = gx*grid->atv[j] + g0*atvX[j];
                    ay = gy*grid->atv[j] + g0*atvY[j];
                    az = gz*grid->atv[j] + g0*atvZ[j];
                    ar = ax*rxa + ay*rya + az*rza;
                    ab = ax*gradwa[0]+ay*gradwa[1]+az*gradwa[2]-ar*zetaa;
                    arb = ax*grid->gradb[0]+ay*grid->gradb[1]+az*grid->gradb[2];  
                    abw = ax*gradwb[0]+ay*gradwb[1]+az*gradwb[2];
                    /* triplet: 
                       abw = -(ax*gradwb[0]+ay*gradwb[1]+az*gradwb[2]);  
                    */
                    data->resa[j+ioff] += 0.5*(fac0*a0+facr*ar+arb*facz+
                                               vxc.df0010*ab+abw*vxc.df00001);
                }
	    }
	}
    }           
    /* beta Fock contribution evaluation */
    fac0 = vxc.df0101*zetab+vxc.df0110*zetaa +vxc.df01001*zetac; 
    if(fabs(rhowb)>0) fac0 += vxc.df0200*rhowb;
    if(fabs(rhowa)>0) fac0 += vxc.df1100*rhowa;
    facr = vxc.df0101*rhowb+vxc.df1001*rhowa
        +vxc.df0002*zetab+vxc.df0011*zetaa+vxc.df00011*zetac;
    /* note: the mixed zetac derivatives not included in facr */
    facz = vxc.df10001*rhowa+vxc.df01001*rhowb
         + vxc.df00011*zetab+vxc.df00101*zetaa+vxc.df00002*zetac;
    /* note: the mixed zetac derivatives not included in facr */
    for(isym=0; isym<inforb_.nsym; isym++) {
        istr = inforb_.ibas[isym];
        iend = inforb_.ibas[isym] + inforb_.nbas[isym];
        jsym = inforb_.muld2h[data->ksymop-1][isym]-1;
        if (isym>=jsym) {
            jstr = inforb_.ibas[jsym];
            jend = inforb_.ibas[jsym] + inforb_.nbas[jsym];
            for(i=istr; i<iend; i++) {
                g0 = grid->atv[i];
                gx = atvX[i];
                gy = atvY[i];
                gz = atvZ[i];
                ioff = i*inforb_.nbast;
                for(j=min(i,jend); j>=jstr; j--) {
                    a0 = g0*grid->atv[j];
                    ax = gx*grid->atv[j] + g0*atvX[j];
                    ay = gy*grid->atv[j] + g0*atvY[j];
                    az = gz*grid->atv[j] + g0*atvZ[j];
                    ar = ax*rxb + ay*ryb + az*rzb;
                    ab = ax*gradwb[0]+ay*gradwb[1]+az*gradwb[2]
                        -ar*zetab;
                    arb = ax*grid->grada[0]+ay*grid->grada[1]
                        +az*grid->grada[2];  
                    abw = ax*gradwa[0]+ay*gradwa[1]+az*gradwa[2];
                    data->resb[j+ioff] += 
                        0.5*(fac0*a0+facr*ar+arb*facz+
                             vxc.df0001*ab+abw*vxc.df00001); 
		}
	    }
	}
    }
}

/* lin_resp_cbab:
 * ZR alpa/beta version of the linear response evaluation routine.
 * triplet: hypotetical case of triplet excitations
 *          should work for general case....       */   
static void
lin_resp_cbab_nogga(DftGrid* grid, LinRespDataab* data)
{
    /*    real b0a, b0b, b3a[3],b3b[3],gradab; */
    real rhowa,rhowb;  
    int isym, i, j;
    int jsym, istr, iend, jstr, jend, ioff;
    FunSecondFuncDrv vxc;    
    real vt, gvi;

    compute_trans_rho(grid, data, &rhowa, &rhowb); 
    if (data->trplet) rhowb = -rhowb; 
    drv2_clear(&vxc);
    selected_func->second(&vxc, grid->curr_weight, &grid->dp);
    /* alpha Fock contribution */
    /* derivatives can in principle diverge, be cautious! */
    vt = 0;
    if(fabs(rhowa)>0) vt += 0.5*vxc.df2000*rhowa;
    if(fabs(rhowb)>0) vt += 0.5*vxc.df1100*rhowb;
    if(fabs(vt)>1e-15) {
        for(isym=0; isym<inforb_.nsym; isym++) {
            jsym = inforb_.muld2h[data->ksymop-1][isym]-1;
            if(isym>=jsym) {
                istr = inforb_.ibas[isym];
                iend = inforb_.ibas[isym] + inforb_.nbas[isym];
                jstr = inforb_.ibas[jsym];
                jend = inforb_.ibas[jsym] + inforb_.nbas[jsym];
                for(i=istr; i<iend; i++) {
                    ioff = i*inforb_.nbast;
                    gvi = vt*grid->atv[i];
                    for(j=min(i,jend-1); j>=jstr; j--)
                        data->resa[j+ioff] += gvi*grid->atv[j];
                }
            }
        } 
    }
    /* beta Fock contribution */ 
    /* derivatives can in principle diverge, be cautious! */
    vt = 0;
    if(fabs(rhowb)>0) vt += 0.5*vxc.df0200*rhowb;
    if(fabs(rhowa)>0) vt += 0.5*vxc.df1100*rhowa;
    if(fabs(vt)>1e-15) {
        for(isym=0; isym<inforb_.nsym; isym++) {
            jsym = inforb_.muld2h[data->ksymop-1][isym]-1;
            if(isym>=jsym) {
                istr = inforb_.ibas[isym];
                iend = inforb_.ibas[isym] + inforb_.nbas[isym];
                jstr = inforb_.ibas[jsym];
                jend = inforb_.ibas[jsym] + inforb_.nbas[jsym];
                for(i=istr; i<iend; i++) {
                    ioff = i*inforb_.nbast;
                    gvi = vt*grid->atv[i];
                    for(j=min(i,jend-1); j>=jstr; j--)
                        data->resb[j+ioff] += gvi*grid->atv[j];
                }
            }    
        }
    }
}

/* ZR linear response for open shell system:
   Evaluates DFT contributions:
   - core fock  matrix      i.e. FCONE
   - open shell fock matix  i.e  FVONE
 */

void
FSYM2(dft_lin_respab)(real* fmatc, real* fmato,  real *cmo, real *zymat, 
		      integer *trplet, integer *ksymop,
		      real* work, integer* lwork)
{
    const real DP5R = 0.5;
    const real MONER = -1.0;
    real electrons = 13.0; 
    struct tms starttm, endtm; clock_t utm;
    LinRespDataab lr_data;
    DftCallbackData cbdata[1];
    int i=1, j, ioff, joff, norbi, norbj;
    integer nvec1=1;
    real averag;
    real *fmata, *fmatb; 
    real *runit;
    DftDensity dens = { dft_dens_unrestricted, NULL, NULL };    
    integer isym, jsym;

    /* WARNING: NO work MAY BE done before syncing slaves! */
    dft_wake_slaves((DFTPropEvalMaster)FSYM2(dft_lin_respab)); /* NO-OP in serial */
    dft_lin_resp_sync_slaves(cmo,&nvec1,&zymat,trplet, ksymop,
                             &work, *lwork);              /* NO-OP in serial */

    times(&starttm);
    /* linear reponse data */
    dens.dmata = lr_data.dmata = malloc(inforb_.n2basx*sizeof(real));
    dens.dmatb = lr_data.dmatb = malloc(inforb_.n2basx*sizeof(real));
    lr_data.resa    = calloc(2*inforb_.n2basx,sizeof(real));
    lr_data.resb    = lr_data.resa + inforb_.n2basx; /* it's an alias only */
    lr_data.kappaa  = calloc(inforb_.n2basx,sizeof(real));
    lr_data.kappab  = calloc(inforb_.n2basx,sizeof(real));
    lr_data.dtgaoa  = malloc(inforb_.nbast *sizeof(real));  
    lr_data.dtgaob  = malloc(inforb_.nbast *sizeof(real)); 
    lr_data.trplet  = *trplet;
    lr_data.ksymop  = *ksymop;

    /* get alpha/beta densities and corresponding kappas */

    FSYM2(dft_get_ao_dens_matab)(cmo,lr_data.dmatb,lr_data.dmata,work,lwork);
    daxpy_(&inforb_.n2basx,&DP5R,lr_data.dmatb,&ONEI,lr_data.dmata,&ONEI);
    dscal_(&inforb_.n2basx,&DP5R,lr_data.dmatb,&ONEI);    
    runit=calloc(inforb_.n2ashx,sizeof(real));
    dunit_(runit,&inforb_.nasht);
    deq27_(cmo,zymat,runit,lr_data.kappab,lr_data.kappaa,work,lwork); 
    free(runit);
    daxpy_(&inforb_.n2basx,&ONER,lr_data.kappab,&ONEI,lr_data.kappaa,&ONEI);
    
    cbdata[0].callback = 
        (DftCallback)(selected_func->is_gga() ? 
                      lin_resp_cbab_gga : lin_resp_cbab_nogga);
    cbdata[0].cb_data  = &lr_data;
    electrons = dft_integrate_ao(&dens,work, lwork,0, 0,0, 
                                 cbdata,ELEMENTS(cbdata));

    dft_lin_respab_collect_info(lr_data.resa, work, *lwork);/*serial:NO-OP*/

    for(i=0; i<inforb_.nbast; i++) {
	ioff = i*inforb_.nbast;
	for(j=0; j<i; j++) {
	    joff = j*inforb_.nbast;
	    averag = lr_data.resa[i+joff] + lr_data.resa[j+ioff];
	    lr_data.resa[i+joff] = lr_data.resa[j+ioff] = averag;
            averag = lr_data.resb[i+joff] + lr_data.resb[j+ioff];
	    lr_data.resb[i+joff] = lr_data.resb[j+ioff] = averag;
	}
    }
    /*    fort_print("Ressa: ");
    output_(lr_data.resa,&ONEI,&inforb_.norbt,&ONEI,&inforb_.norbt,
    &inforb_.norbt, &inforb_.norbt, &ONEI, &priunt_.lupri);
    fort_print("Ressb");
    output_(lr_data.resb,&ONEI,&inforb_.norbt,&ONEI,&inforb_.norbt,
    &inforb_.norbt, &inforb_.norbt, &ONEI, &priunt_.lupri); */

    fmata=calloc(inforb_.n2orbx,sizeof(real));
    fmatb=calloc(inforb_.n2orbx,sizeof(real));   
    for(isym=1; isym<=inforb_.nsym; isym++) {
	jsym  = inforb_.muld2h[lr_data.ksymop-1][isym-1];
        norbi = inforb_.norb[isym-1];
	norbj = inforb_.norb[jsym-1];
	if(norbi>0 && norbj>0)
	    FSYM(autpv)(&isym,&jsym,&cmo[inforb_.icmo[isym-1]],
		   &cmo[inforb_.icmo[jsym-1]], lr_data.resa,
		   inforb_.nbas,&inforb_.nbast,fmata,
		   inforb_.norb,&inforb_.norbt,
		   work, lwork);
    }   
    for(isym=1; isym<=inforb_.nsym; isym++) {
	jsym  = inforb_.muld2h[lr_data.ksymop-1][isym-1];
        norbi = inforb_.norb[isym-1];
	norbj = inforb_.norb[jsym-1];
	if(norbi>0 && norbj>0)
	    FSYM(autpv)(&isym,&jsym,&cmo[inforb_.icmo[isym-1]],
		   &cmo[inforb_.icmo[jsym-1]], lr_data.resb,
		   inforb_.nbas,&inforb_.nbast,fmatb,
		   inforb_.norb,&inforb_.norbt,
		   work, lwork);
    }    
    daxpy_(&inforb_.n2orbx,&TWOR,fmata,&ONEI,fmatc,&ONEI);
    if (inforb_.nasht > 0) {
        if (*trplet) {
             daxpy_(&inforb_.n2orbx,&MONER,fmatb,&ONEI,fmato,&ONEI);
        } else {
             daxpy_(&inforb_.n2orbx,&ONER,fmatb,&ONEI,fmato,&ONEI);
        }
       daxpy_(&inforb_.n2orbx,&MONER,fmata,&ONEI,fmato,&ONEI);
    }
    free(fmata);
    free(fmatb);
    free(lr_data.dmata);
    free(lr_data.dmatb);
    free(lr_data.resa);
    free(lr_data.kappaa);
    free(lr_data.kappab);
    free(lr_data.dtgaoa);
    free(lr_data.dtgaob);
    times(&endtm);
    utm = endtm.tms_utime-starttm.tms_utime;  
    fort_print("      Electrons: %f(%9.1g): LR-DFT evaluation time: %9.1f s", 
               electrons, (double)(electrons-(int)(electrons+0.5)),
               utm/(double)sysconf(_SC_CLK_TCK));
}

/* =================================================================== */
/*                    BLOCKED PROPERTY EVALUATORS                      */
/* =================================================================== */
static void
distribute_lda_bl(DftIntegratorBl *grid, int bllen, int blstart, int blend,
                  real * RESTRICT tmp, real *RESTRICT dR,
                  real * RESTRICT excmat)
{
    int isym, jbl, j, ibl, i, k;
    real * RESTRICT aos = grid->atv;

    for(isym=0; isym<grid->nsym; isym++) {
        integer (*RESTRICT blocks)[2] = BASBLOCK(grid,isym);
        int bl_cnt = grid->bas_bl_cnt[isym];

        for(jbl=0; jbl<bl_cnt; jbl++)
            for(j=blocks[jbl][0]-1; j<blocks[jbl][1]; j++) { 
                int joff = j*bllen;
                for(k=blstart; k<blend; k++)
                    tmp[k+joff] = aos[k+joff]*dR[k];
            }
    
        for(jbl=0; jbl<bl_cnt; jbl++) {
            for(j=blocks[jbl][0]-1; j<blocks[jbl][1]; j++) { 
                real *RESTRICT tmpj = tmp + j*bllen; /* jth orbital */
                real *RESTRICT e_jcol = excmat + j*inforb_.nbast;
                for(ibl=0; ibl<bl_cnt; ibl++) {
                    int top = blocks[ibl][1] < j
                        ? blocks[ibl][1] : j; 
                    for(i=blocks[ibl][0]-1; i<top; i++) { 
                        real *RESTRICT aosi = aos + i*bllen; /* ith orbital */
                        for(k=blstart; k<blend; k++)
                            e_jcol[i] += aosi[k]*tmpj[k];
                    }
                }
                for(k=blstart; k<blend; k++)
                    e_jcol[j] += 0.5*aos[k+j*bllen]*tmpj[k];
            }
        }
    }
}

static void
distribute_gga_bl(DftIntegratorBl *grid, int bllen, int blstart, int blend,
                  real * RESTRICT tmp, real *RESTRICT dR, real *RESTRICT dZ,
                  real * RESTRICT excmat)
{
    int isym, jbl, j, ibl, i, k;
    real * RESTRICT aox = grid->atv+bllen*inforb_.nbast;
    real * RESTRICT aoy = grid->atv+bllen*inforb_.nbast*2;
    real * RESTRICT aoz = grid->atv+bllen*inforb_.nbast*3;
    real * RESTRICT aos = grid->atv;

    for(isym=0; isym<grid->nsym; isym++) {
        integer (*RESTRICT blocks)[2] = BASBLOCK(grid,isym);
        int nblocks = grid->bas_bl_cnt[isym];
        for(jbl=0; jbl<nblocks; jbl++)
            for(j=blocks[jbl][0]-1; j<blocks[jbl][1]; j++) { 
                int joff = j*bllen;
                for(k=0; k<bllen; k++)
                    tmp[k+joff] = 
                        dR[k]* aos[k+j*bllen] +
                        dZ[k]*(aox[k+j*bllen]*grid->g.rad.a[k][0]+
                               aoy[k+j*bllen]*grid->g.rad.a[k][1]+
                               aoz[k+j*bllen]*grid->g.rad.a[k][2]);
        }
        
        for(jbl=0; jbl<nblocks; jbl++) {
            for(j=blocks[jbl][0]-1; j<blocks[jbl][1]; j++) { 
                real *RESTRICT tmpj = tmp + j*bllen; /* jth orbital */
                real *RESTRICT e_jcol = excmat + j*inforb_.nbast;
                for(ibl=0; ibl<nblocks; ibl++) {
                    for(i=blocks[ibl][0]-1; i<blocks[ibl][1]; i++) { 
                        for(k=0; k<bllen; k++)
                            e_jcol[i] += aos[k+i*bllen]*tmpj[k];
                    }
                }
            }
        }
    }
}

/* =================================================================== */
/*                 blocked density and KS evaluation                   */
/* =================================================================== */

struct ks_data {
  real* excmat;
  real* dR;
  real* dZ;
  real energy;
};
static void
kohn_sham_cb_b_lda(DftIntegratorBl *grid, real * RESTRICT tmp,
                   int bllen, int blstart, int blend,
                   struct ks_data* data)
{
    FirstDrv drvs; 
    int i;
    real * RESTRICT excmat = data->excmat;
    real * RESTRICT dR     = data->dR;
    FunDensProp dp = { 0 };
    
    assert(grid->ntypso >0);
    for(i=blstart; i<bllen; i++) {
        real weight = grid->weight[grid->curr_point+i];
        dp.rhoa = dp. rhob = 0.5*grid->r.rho[i];
        data->energy += selected_func->func(&dp)*weight;
        dftpot0_(&drvs, &weight, &dp);
        dR[i] = 2*drvs.fR;
    }

    distribute_lda_bl(grid, bllen, blstart, blend, tmp, dR, excmat);
}

static void
kohn_sham_cb_b_gga(DftIntegratorBl* grid, real * RESTRICT tmp, 
                   int bllen, int blstart, int blend,
                   struct ks_data* data)
{
    FirstDrv drvs;
    int i;
    real * RESTRICT excmat = data->excmat;
    real * RESTRICT dR = data->dR;
    real * RESTRICT dZ = data->dZ;
    FunDensProp dp = { 0 };

    assert(grid->ntypso >0);
    for(i=0; i<bllen; i++) {
        real weight = grid->weight[grid->curr_point+i];
        dp.grada = 0.5*sqrt(grid->g.grad[i][0]*grid->g.grad[i][0]+
                            grid->g.grad[i][1]*grid->g.grad[i][1]+
                            grid->g.grad[i][2]*grid->g.grad[i][2]);
        dp. rhoa = dp. rhob = 0.5*grid->r.rho[i];
        dp.gradb  = dp.grada;
        dp.gradab = dp.grada*dp.gradb;
        data->energy += selected_func->func(&dp)*weight;
        dftpot0_(&drvs, &weight, &dp);
        dR[i] = drvs.fR;
        dZ[i] = drvs.fZ/dp.grada;
    }

    distribute_gga_bl(grid, bllen, blstart, blend, tmp, dR, dZ, excmat);
}

/* dft_kohn_sham:
   compute Fock matrix ksm corresponding to given density matrix dmat.
   fast version - uses memory bandwidth-efficient algorithm.
*/
void
FSYM2(dft_kohn_shamf)(real* dmat, real* ksm, real* edfty,
		      real* work, integer *lwork, integer* iprfck)
{
    int nbast2, i, j;
    struct tms starttm, endtm; clock_t utm;
    real electrons;
    struct ks_data ds;

    nbast2   = inforb_.nbast*inforb_.nbast;
    ds.excmat = calloc(nbast2, sizeof(real));
    ds.dR     = dal_malloc(DFT_BLLEN*sizeof(real));
    ds.dZ     = dal_malloc(DFT_BLLEN*sizeof(real));
        
    times(&starttm);
    ds.energy = 0.0;

    electrons = dft_integrate_ao_bl(1, dmat, work, lwork, 0, 
                                    (DftBlockCallback)
                                    (selected_func->is_gga() ?
                                    kohn_sham_cb_b_gga : kohn_sham_cb_b_lda),
                                    &ds);

    for(i=0; i<inforb_.nbast; i++) {
	int ioff = i*inforb_.nbast;
	for(j=0; j<i; j++) {
	    int joff = j*inforb_.nbast;
	    real averag = 0.5*(ds.excmat[i+joff] + ds.excmat[j+ioff]);
	    ds.excmat[i+joff] = ds.excmat[j+ioff] = averag;
	}
    }

    *edfty=ds.energy;
    daxpy_(&inforb_.n2basx, &ONER, ds.excmat, &ONEI, ksm, &ONEI);
    if(KOHNSH_DEBUG) {
        fort_print("kohn sham matrix");
        outmat_(ds.excmat,&ONEI,&inforb_.nbast,&ONEI,&inforb_.nbast,
                &inforb_.nbast, &inforb_.nbast);
    }

    free(ds.excmat);
    free(ds.dR);
    free(ds.dZ);
    times(&endtm);
    utm = endtm.tms_utime-starttm.tms_utime;
    fort_print("      Electrons: %11.7f %7.1g: Energy %f KS/B time: %9.1f s", 
               electrons, (electrons-2.0*inforb_.nrhft)/(2.0*inforb_.nrhft), 
               ds.energy, utm/(double)sysconf(_SC_CLK_TCK));
}

/* ------------------------------------------------------------------- */
/* ------ blocked DFT LINEAR RESPONSE CONTRIBUTION EVALUATOR --------- */
/* ------------------------------------------------------------------- */

typedef struct {
    real* dmat, *kappa, *res;
    real* dtgao;
    real* vt; /* dimensioned [bllen] for LDA, [bllen][4] for GGA */
    integer trplet, ksymop;
    integer vecs_in_batch;
} LinRespBlData;

static void
lin_resp_cb_b_lda(DftIntegratorBl* grid, real * RESTRICT tmp,
                  int bllen, int blstart, int blend,
                  LinRespBlData* data)
{
    real * RESTRICT aos = grid->atv;
    real * RESTRICT excmat = data->res;
    real (* RESTRICT vt) = data->vt; /* [bllen][4] */
    int ibl, i, jbl, j, k, isym, ivec;
    FunDensProp dp = { 0 };
    integer blen = bllen;

    for(ivec=0; ivec<data->vecs_in_batch; ivec++) {
        /* compute vector of transformed densities vt */
        FSYM2(getexp_blocked_lda)(&data->ksymop, data->kappa + ivec*inforb_.n2basx,
                                  grid->atv, grid->bas_bl_cnt, grid->basblocks,
                                  &grid->shl_bl_cnt, tmp, &blen, vt);

        for(i=blstart; i<blend; i++) {
            SecondDrv vxc;
            real weight = grid->weight[grid->curr_point+i];
            dp.rhoa = dp.rhob = 0.5*grid->r.rho[i];
            dftpot1_(&vxc, &weight, &dp, &data->trplet);
            vt[i] = vxc.fRR*vt[i]*2;
        }

        for(isym=0; isym<grid->nsym; isym++) {
            integer (*RESTRICT iblocks)[2] = BASBLOCK(grid,isym);
            int ibl_cnt = grid->bas_bl_cnt[isym];
            
            for(ibl=0; ibl<ibl_cnt; ibl++)
                for(i=iblocks[ibl][0]-1; i<iblocks[ibl][1]; i++) { 
                    int ioff = i*bllen;
                    for(k=blstart; k<blend; k++)
                    tmp[k+ioff] = aos[k+ioff]*vt[k];
                }
            
            for(ibl=0; ibl<ibl_cnt; ibl++) {
                for(i=iblocks[ibl][0]-1; i<iblocks[ibl][1]; i++) { 
                    int ioff = i*inforb_.nbast + ivec*inforb_.n2basx;
                    int jsym = inforb_.muld2h[data->ksymop-1][isym]-1;
                    integer (*RESTRICT jblocks)[2] = BASBLOCK(grid,jsym);
                    int jbl_cnt = grid->bas_bl_cnt[jsym];
                    real *RESTRICT tmpi = &tmp[i*bllen];
                    if (isym<jsym) continue;
                    for(jbl=0; jbl<jbl_cnt; jbl++) {
                        int jtop = min(jblocks[jbl][1],i);
                        for(j=jblocks[jbl][0]-1; j<jtop; j++) { 
                            for(k=blstart; k<blend; k++)
                                excmat[j+ioff] +=
                                    aos[k+j*bllen]*tmpi[k];
                        }
                    }
                    for(k=blstart; k<blend; k++)
                        excmat[i+ioff] += aos[k+i*bllen]*tmpi[k]*0.5;
                }
            }
        }
    }
}

static void
lin_resp_cb_b_gga(DftIntegratorBl* grid, real * RESTRICT tmp,
                  int bllen, int blstart, int blend,
                  LinRespBlData* data)
{
    int ibl, i, jbl, j, k, isym, ivec;
    real (* RESTRICT vt3)[4] = (real(*)[4])data->vt; /* [bllen][4] */
    real * RESTRICT aos = grid->atv;
    real * RESTRICT aox = grid->atv+bllen*inforb_.nbast;
    real * RESTRICT aoy = grid->atv+bllen*inforb_.nbast*2;
    real * RESTRICT aoz = grid->atv+bllen*inforb_.nbast*3;
    real * RESTRICT excmat = data->res;
    FunDensProp dp = { 0 };
    integer blen = bllen;

    for(ivec=0; ivec<data->vecs_in_batch; ivec++) {
        /* compute vector of transformed densities and dens. gradients vt3 */
        FSYM2(getexp_blocked_gga)(&data->ksymop,data->kappa + ivec*inforb_.n2basx,
                                  grid->atv, grid->bas_bl_cnt,
                                  grid->basblocks, &grid->shl_bl_cnt, tmp,
				  &blen, vt3);
        for(i=blstart; i<blend; i++) {
            SecondDrv vxc;
            real facr, facg;
            real weight = grid->weight[grid->curr_point+i];
            real ngrad  = sqrt(grid->g.grad[i][0]*grid->g.grad[i][0]+
                               grid->g.grad[i][1]*grid->g.grad[i][1]+
                               grid->g.grad[i][2]*grid->g.grad[i][2]);
            real brg, brz, b0 = vt3[i][0];
            if(ngrad<1e-15|| grid->r.rho[i]<1e-15) {
                vt3[i][0] = vt3[i][1] = vt3[i][2] = vt3[i][3] = 0;
                continue;
            }
            brg = (vt3[i][1]*grid->g.grad[i][0] +
                   vt3[i][2]*grid->g.grad[i][1] +
                   vt3[i][3]*grid->g.grad[i][2]);
            brz = brg/ngrad;
            dp. rhoa = dp. rhob = 0.5*grid->r.rho[i];
            dp.grada = dp.gradb = 0.5*ngrad;
            dp.gradab = dp.grada*dp.gradb;
            dftpot1_(&vxc, &weight, &dp, &data->trplet);
            facr = vxc.fRZ*b0 + (vxc.fZZ-vxc.fZ/ngrad)*brz + vxc.fZG*brg;
            facr = facr/ngrad + (vxc.fRG*b0+vxc.fZG*brz +vxc.fGG*brg);
            facg = vxc.fZ/ngrad + vxc.fG;
            vt3[i][0] = vxc.fRR*b0 + vxc.fRZ*brz+ vxc.fRG*brg;
            vt3[i][1] = (grid->g.grad[i][0]*facr + facg*vt3[i][1])*2;
            vt3[i][2] = (grid->g.grad[i][1]*facr + facg*vt3[i][2])*2;
            vt3[i][3] = (grid->g.grad[i][2]*facr + facg*vt3[i][3])*2;
        }

        for(isym=0; isym<grid->nsym; isym++) {
            integer (*RESTRICT iblocks)[2] = BASBLOCK(grid,isym);
            int ibl_cnt = grid->bas_bl_cnt[isym];

            for(ibl=0; ibl<ibl_cnt; ibl++) {
                for(i=iblocks[ibl][0]-1; i<iblocks[ibl][1]; i++) { 
                    real *RESTRICT g0i = &aos[i*bllen];
                    real *RESTRICT gxi = &aox[i*bllen];
                    real *RESTRICT gyi = &aoy[i*bllen];
                    real *RESTRICT gzi = &aoz[i*bllen];
                    int ioff = i*inforb_.nbast + ivec*inforb_.n2basx;
                    int jsym = inforb_.muld2h[data->ksymop-1][isym]-1;
                    integer (*RESTRICT jblocks)[2] = BASBLOCK(grid,jsym);
                    int jbl_cnt = grid->bas_bl_cnt[jsym];
                    for(k=blstart; k<blend; k++)
                        tmp[k] = (vt3[k][0]*g0i[k] + 
                                  vt3[k][1]*gxi[k] +
                                  vt3[k][2]*gyi[k] +
                                  vt3[k][3]*gzi[k]);
                    for(jbl=0; jbl<jbl_cnt; jbl++) {
                        int jtop = jblocks[jbl][1];
                        for(j=jblocks[jbl][0]-1; j<jtop; j++) { 
                            real *RESTRICT g0j = &aos[j*bllen];
                            real s = 0;
                            for(k=blstart; k<blend; k++)
                                s += g0j[k]*tmp[k];
                            excmat[j+ioff] += s;
                        }
                    }
                }
            }
        }
    }
}

/* dft_lin_respf_:
   Main Linear Response evaluator.
   FMAT(NORBT,NORBT,NOSIM) - result added to FMAT. Must not be referenced on slaves.
   cmo -const
   zymat(NORBT,NORBT,NOSIM) - const response vector.
   trplet - triplet excitation? (bool)
   NOSIM - number of simultaneously transformed response vectors.
*/
void
FSYM2(dft_lin_respf)(integer *nosim, real* fmat, real *cmo, real *zymat,
		     integer *trplet, integer *ksymop,
		     real* work, integer* lwork)
{
    /* MAX_VEC determines the number of simultaneusly transformed
       vectors.  Set it to a small mumber (eg. 3) - the only operation
       that is saved is evaluation of orbitals and density and density
       gradients which is usually small compared to the time spent in
       actual vector transformation. Larger values will only increate
       memory utilization without positive impact on performance.
       FIXME: consider using work for this purpose. */
    static const int MAX_VEC = 5;
    struct tms starttm, endtm; clock_t utm;
    real electrons = 0;
    LinRespBlData lr_data; /* linear response data */
    real dummy;
    int i, j, ivec, jvec;
    int max_vecs;

    /* WARNING: NO work MAY BE done before syncing slaves! */
    dft_wake_slaves((DFTPropEvalMaster)FSYM2(dft_lin_respf)); /* NO-OP in serial */
    dft_lin_resp_sync_slaves(cmo,nosim,&zymat,trplet,ksymop,
                             &work, *lwork);              /* NO-OP in serial */

    times(&starttm);
    max_vecs = *nosim>MAX_VEC ? MAX_VEC : *nosim;
    lr_data.dmat  = dal_malloc(inforb_.n2basx*sizeof(real));
    lr_data.res   = dal_malloc(inforb_.n2basx*sizeof(real)*max_vecs);
    lr_data.kappa = dal_malloc(inforb_.n2basx*sizeof(real)*max_vecs);
    lr_data.vt    = dal_malloc(DFT_BLLEN*4   *sizeof(real));
    lr_data.dtgao = dal_malloc(inforb_.nbast *sizeof(real));
    lr_data.trplet= *trplet;
    lr_data.ksymop= *ksymop;
    FSYM2(dft_get_ao_dens_mat)(cmo, lr_data.dmat, work, lwork);
    
    for(ivec=0; ivec<*nosim; ivec+=max_vecs) {
        integer sz;
        lr_data.vecs_in_batch = ivec + max_vecs > *nosim ? *nosim - ivec : max_vecs;
        sz = lr_data.vecs_in_batch * inforb_.n2basx;
        FSYM(dzero)(lr_data.kappa, &sz);
        for(jvec=0; jvec<lr_data.vecs_in_batch; jvec++)
            FSYM(deq27)(cmo,zymat+(ivec+jvec)*inforb_.n2orbx,&dummy,
                        lr_data.kappa+jvec*inforb_.n2basx, &dummy,work,lwork);
        FSYM(dzero)(lr_data.res, &sz);
        electrons = dft_integrate_ao_bl(1, lr_data.dmat, work, lwork, 0, 
                                        (DftBlockCallback)
                                        (selected_func->is_gga() ?
                                         lin_resp_cb_b_gga : lin_resp_cb_b_lda),
                                        &lr_data);
#ifdef VAR_MPI
        dft_lin_resp_collect_info(lr_data.vecs_in_batch,lr_data.res, work, *lwork);
	if(fmat == NULL)  /* The transformations below are done only by master. */
	    continue;
#endif
        if(DFTLR_DEBUG) {
            fort_print("AO Fock matrix contribution in dft_lin_respf");
            outmat_(lr_data.res,&ONEI,&inforb_.nbast,&ONEI,&inforb_.nbast,
                    &inforb_.nbast, &inforb_.nbast);
        }
        for(jvec=0; jvec<lr_data.vecs_in_batch; jvec++){
            for(i=0; i<inforb_.nbast; i++) {
                int ioff = i*inforb_.nbast + jvec*inforb_.n2basx;
                for(j=0; j<i; j++) {
                    int joff = j*inforb_.nbast + jvec*inforb_.n2basx;
                    real averag = 0.5*(lr_data.res[i+joff] + 
                                       lr_data.res[j+ioff]);
                    lr_data.res[i+joff] = lr_data.res[j+ioff] = averag;
                }
            }
            /* transform to MO Fock matrix contribution  */
            FSYM(lrao2mo)(cmo, &lr_data.ksymop, lr_data.res+jvec*inforb_.n2basx,
                          fmat+(ivec+jvec)*inforb_.n2orbx, work, lwork);
        }
        if(DFTLR_DEBUG) {
            fort_print("MO Fock matrix contribution in dft_lin_resp");
            outmat_(fmat,&ONEI,&inforb_.norbt,&ONEI,&inforb_.norbt,
                    &inforb_.norbt, &inforb_.norbt);
        }
    }

    free(lr_data.dmat);
    free(lr_data.res);
    free(lr_data.kappa);
    free(lr_data.vt);
    free(lr_data.dtgao);
    times(&endtm);
    utm = endtm.tms_utime-starttm.tms_utime;
    fort_print("      Electrons: %f(%9.3g): LR-DFT*%d evaluation time: %9.1f s", 
               electrons, (double)(electrons-(int)(electrons+0.5)), *nosim,
               utm/(double)sysconf(_SC_CLK_TCK));
}

void
dft_lin_respf_slave(real* work, integer* lwork, const integer* iprint)
{
    real *cmo  = malloc(inforb_.norbt*inforb_.nbast*sizeof(real)); /* IN  */
    integer trplet;                     /* IN: will be synced from master */
    integer ksymop, nosim;
    FSYM2(dft_lin_respf)(&nosim, NULL, cmo, NULL,
                         &trplet, &ksymop, work, lwork);
    free(cmo);
}


/*====================================================================*/
/*                     BLOCKED OPEN-SHELL CODE                        */
/*                   Z. Rinkevicius, June 2008                        */
/*====================================================================*/

/*                    Kohn-Sham evaluator                             */ 

struct ks_data_ab {
  real *ksma, *ksmb; 
  real energy;
  real *dRa; 
  real *dRb; 
  real (*dZa);  
  real (*dZb);  
  real (*dZab); 
  real *tmpa, *tmpb; 
};


#if defined(VAR_MPI)
void
dft_kohn_shamab_b_slave(real* work, integer* lwork, const integer* iprint)
{
  real* dmat = malloc(2*inforb_.n2basx*sizeof(real));
  integer iprfck = 0;
  FSYM2(dft_kohn_shamab_b)(dmat, NULL, NULL, work, lwork, &iprfck);
  free(dmat);
}

void
dft_kohn_shamab_b_sync_slaves(real* dmat)
{
  FSYM(dftintbcast)();
  MPI_Bcast(dmat, 2*inforb_.n2basx,MPI_DOUBLE,
            MASTER_NO, MPI_COMM_WORLD);
}

static void
dft_kohn_shamab_b_collect(real *energy, real *ksm, real* work, integer lwork)
{
  real tmp = *energy;
  int sz = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &sz);
  if(sz <=1) return; 
  sz = 2*inforb_.n2basx;
  integer sz_ftn = sz;
  CHECK_WRKMEM(sz, lwork);
  dcopy_(&sz_ftn, ksm,&ONEI, work, &ONEI);
  MPI_Reduce(work, ksm, sz, MPI_DOUBLE, MPI_SUM, 
             MASTER_NO, MPI_COMM_WORLD);
  MPI_Reduce(&tmp, energy, 1, MPI_DOUBLE, MPI_SUM, 
             MASTER_NO, MPI_COMM_WORLD);
}
#endif

static void
distribute_lda_ab_bl(DftIntegratorBl *grid, int bllen, int blstart, int blend,
                     real * RESTRICT tmpa,real * RESTRICT tmpb,
                     real *RESTRICT dRa, real *RESTRICT dRb,
                     real * RESTRICT exc_a, real * RESTRICT exc_b)
{
  int isym, jbl, j, ibl, i, k;
  real * RESTRICT aos = grid->atv;

  for(isym=0; isym<grid->nsym; isym++) {
    integer (*RESTRICT blocks)[2] = BASBLOCK(grid,isym);
    int bl_cnt = grid->bas_bl_cnt[isym];
      
    for(jbl=0; jbl<bl_cnt; jbl++)
      for(j=blocks[jbl][0]-1; j<blocks[jbl][1]; j++) {
        int joff = j*bllen;
        real *RESTRICT e_jcola = exc_a+ j*inforb_.nbast;
        real *RESTRICT e_jcolb = exc_b+ j*inforb_.nbast;
        for(k=blstart; k<blend; k++) {
          tmpa[k] = aos[k+joff]*dRa[k];
          tmpb[k] = aos[k+joff]*dRb[k];
        }      
        for(ibl=0; ibl<bl_cnt; ibl++) {
                    int top = blocks[ibl][1] < j
                      ? blocks[ibl][1] : j;
                    for(i=blocks[ibl][0]-1; i<top; i++) {
                      real *RESTRICT aosi = aos + i*bllen; /* ith orbital */
                      for(k=blstart; k<blend; k++) {
                        e_jcola[i] += aosi[k]*tmpa[k];
                        e_jcolb[i] += aosi[k]*tmpb[k];
                      }
                    }
        }
        for(k=blstart; k<blend; k++) {
          e_jcola[j] += 0.5*aos[k+j*bllen]*tmpa[k];
          e_jcolb[j] += 0.5*aos[k+j*bllen]*tmpb[k];
        }
      }
  }
}

static void
distribute_gga_ab_bl(DftIntegratorBl *grid, int bllen, int blstart, int blend,
                     real * RESTRICT tmpa, real * RESTRICT tmpb,
                     real *RESTRICT dRa, real *RESTRICT dRb,
                     real *RESTRICT dZa, real *RESTRICT dZb,
                     real *RESTRICT dZab,
                     real * RESTRICT exc_a, real * RESTRICT exc_b)
{
  int isym, jbl, j, ibl, i, k;
  real * RESTRICT aox = grid->atv+bllen*inforb_.nbast;
  real * RESTRICT aoy = grid->atv+bllen*inforb_.nbast*2;
  real * RESTRICT aoz = grid->atv+bllen*inforb_.nbast*3;
  real * RESTRICT aos = grid->atv;

  for(isym=0; isym<grid->nsym; isym++) {
    integer (*RESTRICT blocks)[2] = BASBLOCK(grid,isym);
    int nblocks = grid->bas_bl_cnt[isym];
    for(jbl=0; jbl<nblocks; jbl++)
      for(j=blocks[jbl][0]-1; j<blocks[jbl][1]; j++) {
        int joff = j*bllen;
        real *RESTRICT e_jcola = exc_a + j*inforb_.nbast;
        real *RESTRICT e_jcolb = exc_b + j*inforb_.nbast;
        for(k=blstart; k<blend; k++) {
          tmpa[k] =
            dRa[k]* aos[k+joff]+
            dZa[k]*(aox[k+joff]*grid->g.rad.a[k][0]+
                    aoy[k+joff]*grid->g.rad.a[k][1]+
                    aoz[k+joff]*grid->g.rad.a[k][2])+
            dZab[k]*(aox[k+joff]*grid->g.rad.b[k][0]+
                     aoy[k+joff]*grid->g.rad.b[k][1]+
                     aoz[k+joff]*grid->g.rad.b[k][2]);
          tmpb[k] =
            dRb[k]* aos[k+joff]+
            dZb[k]*(aox[k+joff]*grid->g.rad.b[k][0]+
                    aoy[k+joff]*grid->g.rad.b[k][1]+
                    aoz[k+joff]*grid->g.rad.b[k][2])+
            dZab[k]*(aox[k+joff]*grid->g.rad.a[k][0]+
                     aoy[k+joff]*grid->g.rad.a[k][1]+
                     aoz[k+joff]*grid->g.rad.a[k][2]);
        }
        for(ibl=0; ibl<nblocks; ibl++) {
          for(i=blocks[ibl][0]-1; i<blocks[ibl][1]; i++) {
            for(k=blstart; k<blend; k++) {
              e_jcola[i] += aos[k+i*bllen]*tmpa[k];
              e_jcolb[i] += aos[k+i*bllen]*tmpb[k];
            }
          }
        }
      }
  }
}

static void
kohn_shamab_cb_b_lda(DftIntegratorBl *grid, real * RESTRICT tmp,
                     int bllen, int blstart, int blend,
                     struct ks_data_ab* data)
{

  FunFirstFuncDrv drvs;
  int i;
  real * RESTRICT exc_a  = data->ksma;
  real * RESTRICT exc_b  = data->ksmb;
  real * RESTRICT dRa    = data->dRa;
  real * RESTRICT dRb    = data->dRb;
  FunDensProp dp = { 0 };

  assert(grid->ntypso >0);
  for(i=blstart; i<blend; i++) {
    real weight = grid->weight[grid->curr_point+i];
    dp.rhoa = grid->r.ho.a[i];
    dp.rhob = grid->r.ho.b[i];
    data->energy += selected_func->func(&dp)*weight;
    drv1_clear(&drvs);
    selected_func->first(&drvs, weight, &dp);
    dRa[i] = 2.0*drvs.df1000;
    dRb[i] = 2.0*drvs.df0100;
  }
  distribute_lda_ab_bl(grid, bllen, blstart, blend, data->tmpa, data->tmpb,
                       dRa, dRb,  exc_a, exc_b);
}

static void
kohn_shamab_cb_b_gga(DftIntegratorBl* grid, real * RESTRICT tmp,
                     int bllen, int blstart, int blend,
                     struct ks_data_ab* data)
{
  FunFirstFuncDrv drvs;
  int i;
  real * RESTRICT exc_a   = data->ksma;
  real * RESTRICT exc_b   = data->ksmb;
  real * RESTRICT dRa     = data->dRa;
  real * RESTRICT dRb     = data->dRb;
  real * RESTRICT dZa     = data->dZa;
  real * RESTRICT dZb     = data->dZb;
  real * RESTRICT dZab    = data->dZab;
  FunDensProp dp = { 0 };

  assert(grid->ntypso >0);
  for(i=blstart; i<blend; i++) {
    real weight = grid->weight[grid->curr_point+i];
    dp.rhoa = grid->r.ho.a[i];
    dp.rhob = grid->r.ho.b[i];
    dp.grada = sqrt(grid->g.rad.a [i][0]*grid->g.rad.a [i][0]+
                    grid->g.rad.a [i][1]*grid->g.rad.a [i][1]+
                    grid->g.rad.a [i][2]*grid->g.rad.a [i][2]);
    dp.grada = fabs(dp.grada)>1e-40 ? dp.grada : 1e-40;
    dp.gradb = sqrt(grid->g.rad.b [i][0]*grid->g.rad.b [i][0]+
                    grid->g.rad.b [i][1]*grid->g.rad.b [i][1]+
                    grid->g.rad.b [i][2]*grid->g.rad.b [i][2]);
    dp.gradb  = fabs(dp.gradb)>1e-40 ? dp.gradb : 1e-40;
    dp.gradab = grid->g.rad.a [i][0]*grid->g.rad.b [i][0]+
                grid->g.rad.a [i][1]*grid->g.rad.b [i][1]+
      grid->g.rad.a [i][2]*grid->g.rad.b [i][2];
    data->energy += selected_func->func(&dp)*weight;
    drv1_clear(&drvs);
    selected_func->first(&drvs, weight, &dp);
    dRa[i]    = drvs.df1000;
    dRb[i]    = drvs.df0100;
    dZa[i]    = 2.0*drvs.df0010/dp.grada;
    dZb[i]    = 2.0*drvs.df0001/dp.gradb;
    dZab[i]   = 2.0*drvs.df00001;
  }
  distribute_gga_ab_bl(grid, bllen, blstart, blend, data->tmpa, data->tmpb,
                       dRa, dRb, dZa, dZb, dZab, exc_a, exc_b);
}

void
FSYM2(dft_kohn_shamab_b)(real* dmat, real* ksm, real *edfty,
                         real* work, integer *lwork, integer* iprfck)
{
  int i, j;
  integer sz;
  struct ks_data_ab result; /* res like result */
  struct tms starttm, endtm; clock_t utm;
  real electrons, exp_el = 2.0*inforb_.nrhft+inforb_.nasht;
  
  /* parallel code insert */
#if defined(VAR_MPI)
  dft_wake_slaves((DFTPropEvalMaster)FSYM2(dft_kohn_shamab_b)); 
  dft_kohn_shamab_b_sync_slaves(dmat);
#endif  
  /* memory allocation */
  result.ksma = calloc(2*inforb_.n2basx, sizeof(real));
  result.ksmb = result.ksma + inforb_.n2basx;
  result.dRa     = dal_malloc(DFT_BLLEN*sizeof(real));
  result.dRb     = dal_malloc(DFT_BLLEN*sizeof(real));
  result.dZa     = dal_malloc(DFT_BLLEN*sizeof(real));
  result.dZb     = dal_malloc(DFT_BLLEN*sizeof(real));
  result.dZab    = dal_malloc(DFT_BLLEN*sizeof(real));
  result.tmpa    = dal_malloc(DFT_BLLEN*sizeof(real));
  result.tmpb    = dal_malloc(DFT_BLLEN*sizeof(real));
  result.energy = 0.0;
  /* start timer*/
  times(&starttm);
  /* integeration */
  electrons = dft_integrate_ao_bl(2, dmat, work, lwork, 0,
                                  (DftBlockCallback)(selected_func->is_gga() ?
						     kohn_shamab_cb_b_gga:kohn_shamab_cb_b_lda),
                                  &result);
  /* parallel code insert */
#if defined(VAR_MPI)
  /* collect data */
  dft_kohn_shamab_b_collect(&result.energy,result.ksma, work, *lwork);
  /* free memory and leave if slave*/
  if (edfty == NULL) {
    free(result.ksma);
    free(result.dRa);
    free(result.dRb);
    free(result.dZa);
    free(result.dZb);
    free(result.dZab);
    free(result.tmpa);
    free(result.tmpb);
    return;
  }; 
#endif 
  /* average Kohn-Sham matrices */
  for(i=0; i<inforb_.nbast; i++) {
    int ioff = i*inforb_.nbast;
    for(j=0; j<i; j++) {
      int joff = j*inforb_.nbast;
      real averag = 0.5*(result.ksma[i+joff] + result.ksma[j+ioff]);
      result.ksma[i+joff] = result.ksma[j+ioff] = averag;
      averag = 0.5*(result.ksmb[i+joff] + result.ksmb[j+ioff]);
      result.ksmb[i+joff] = result.ksmb[j+ioff] = averag;
    }
  }
  /* results */
  *edfty=result.energy;
  sz = 2*inforb_.n2basx;
  FSYM(daxpy)(&sz, &ONER, result.ksma, &ONEI, ksm, &ONEI);
  /* free memory */
  free(result.ksma);
  free(result.dRa);
  free(result.dRb);
  free(result.dZa);
  free(result.dZb);
  free(result.dZab);
  free(result.tmpa);
  free(result.tmpb);
  /* timings & print out*/
  times(&endtm);
  utm = endtm.tms_utime-starttm.tms_utime;
  fort_print("      Electrons: %11.7f %9.1g: Energy %f KS-AB/B time: %9.1f s",
             electrons, (electrons-exp_el)/exp_el,
             result.energy, utm/(double)sysconf(_SC_CLK_TCK));
}

/* Linear response contribution evaluator */

typedef struct {
  real* dmata, *dmatb;
  real* kappa_a, *kappa_b;
  real* res_a, *res_b;
  real* dtgaoa, *dtgaob;
  real* rhowa, *rhowb;
  real* vta, *vtb;
  real* tmpa, *tmpb;
  integer trplet, ksymop, vecs_in_batch;
} LinResp_ab_BlData;


#if defined(VAR_MPI)
void
dft_lin_respab_b_slave(real* work, integer* lwork, const integer* iprint)
{
  real *cmo  = malloc(inforb_.norbt*inforb_.nbast*sizeof(real)); 
  integer trplet;                    
  integer ksymop, nosim;
  FSYM2(dft_lin_respab_b)(&nosim, NULL, NULL, cmo, NULL,
			  &trplet, &ksymop, work, lwork);
  free(cmo);
}

static  void
dft_lin_respab_b_collect(integer nvec, real* fmata, real* fmatb,  real *work, integer lwork)
{
  int sz = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &sz);
  if(sz<=1) return;

  sz = nvec*inforb_.n2basx;
  integer sz_ftn =sz;
  if(sz>lwork) {
    CHECK_WRKMEM(inforb_.n2basx,lwork);
    for(sz=0; sz<nvec; sz++) {
      dcopy_(&inforb_.n2basx, fmata + sz*inforb_.n2basx,
             &ONEI, work, &ONEI);
      MPI_Reduce(work, fmata+sz*inforb_.n2basx, inforb_.n2basx,
                 MPI_DOUBLE, MPI_SUM, MASTER_NO, MPI_COMM_WORLD);
    }
    for(sz=0; sz<nvec; sz++) {
      dcopy_(&inforb_.n2basx, fmatb + sz*inforb_.n2basx,
             &ONEI, work, &ONEI);
      MPI_Reduce(work, fmatb+sz*inforb_.n2basx, inforb_.n2basx,
                 MPI_DOUBLE, MPI_SUM, MASTER_NO, MPI_COMM_WORLD);
    }
  } else {
    CHECK_WRKMEM(sz, lwork);
    dcopy_(&sz_ftn, fmata, &ONEI, work, &ONEI);
    MPI_Reduce(work, fmata, sz, MPI_DOUBLE, MPI_SUM, 
               MASTER_NO, MPI_COMM_WORLD);
    dcopy_(&sz_ftn, fmatb, &ONEI, work, &ONEI);
    MPI_Reduce(work, fmatb, sz, MPI_DOUBLE, MPI_SUM,
               MASTER_NO, MPI_COMM_WORLD); 
  }
}
#endif 

static void
lin_resp_cb_b_lda_ab(DftIntegratorBl* grid, real * RESTRICT tmp,
		     int bllen, int blstart, int blend,
		     LinResp_ab_BlData* data)
{
  static const real MONER = -1.0;
  real * RESTRICT aos = grid->atv;
  real * RESTRICT rhowa = data->rhowa;
  real * RESTRICT rhowb = data->rhowb;
  real * RESTRICT exc_a = data->res_a;
  real * RESTRICT exc_b = data->res_b;
  real * RESTRICT vta   = data->vta;
  real * RESTRICT vtb   = data->vtb;
  real * RESTRICT tmpa  = data->tmpa;
  real * RESTRICT tmpb  = data->tmpb;
  int ibl, i, jbl, j, k, isym, ivec;
  integer bl_size = DFT_BLLEN;
  FunSecondFuncDrv vxc;
  FunDensProp dp = { 0 };
  integer blen = bllen;

    
  for(ivec=0; ivec<data->vecs_in_batch; ivec++) {
    /* evaluate transformed densities */
    FSYM2(getexp_blocked_lda)(&data->ksymop, data->kappa_a + ivec*inforb_.n2basx,
			      grid->atv, grid->bas_bl_cnt, grid->basblocks,
			      &grid->shl_bl_cnt, tmp, &blen, rhowa);
    FSYM2(getexp_blocked_lda)(&data->ksymop, data->kappa_b + ivec*inforb_.n2basx,
			      grid->atv, grid->bas_bl_cnt, grid->basblocks,
			      &grid->shl_bl_cnt, tmp, &blen, rhowb);
    /* change sign if triplet response */
    if (data->trplet) FSYM(dscal)(&bl_size,&MONER,rhowb,&ONEI);
    for(i=blstart; i<blend; i++) {
      real weight = grid->weight[grid->curr_point+i];
      dp.rhoa = grid->r.ho.a[i];
      dp.rhob = grid->r.ho.b[i];
      drv2_clear(&vxc);
      selected_func->second(&vxc,weight,&dp);
      vta[i] = vxc.df2000*rhowa[i]+vxc.df1100*rhowb[i];
      vtb[i] = vxc.df0200*rhowb[i]+vxc.df1100*rhowa[i];
    }
    for(isym=0; isym<grid->nsym; isym++) {
      integer (*RESTRICT iblocks)[2] = BASBLOCK(grid,isym);
      int ibl_cnt = grid->bas_bl_cnt[isym];
      for(ibl=0; ibl<ibl_cnt; ibl++)
	for(i=iblocks[ibl][0]-1; i<iblocks[ibl][1]; i++) {
	  int ioff = i*bllen;
	  for(k=blstart; k<blend; k++){
	    tmpa[k+ioff] = aos[k+ioff]*vta[k];
	    tmpb[k+ioff] = aos[k+ioff]*vtb[k];
	  }
	}
      for(ibl=0; ibl<ibl_cnt; ibl++) {
	for(i=iblocks[ibl][0]-1; i<iblocks[ibl][1]; i++) {
	  int ioff = i*inforb_.nbast + ivec*inforb_.n2basx;
	  int jsym = inforb_.muld2h[data->ksymop-1][isym]-1;
	  integer (*RESTRICT jblocks)[2] = BASBLOCK(grid,jsym);
	  int jbl_cnt = grid->bas_bl_cnt[jsym];
	  real *RESTRICT tmpai = &tmpa[i*bllen];
	  real *RESTRICT tmpbi = &tmpb[i*bllen];
	  if (isym<jsym) continue;
	  for(jbl=0; jbl<jbl_cnt; jbl++) {
	    int jtop = min(jblocks[jbl][1],i);
	    for(j=jblocks[jbl][0]-1; j<jtop; j++) {
	      for(k=blstart; k<blend; k++) {
		exc_a[j+ioff] += aos[k+j*bllen]*tmpai[k];
		exc_b[j+ioff] += aos[k+j*bllen]*tmpbi[k];
	      }
	    }
	  }
	  for(k=blstart; k<blend; k++) {
	    exc_a[i+ioff] += aos[k+i*bllen]*tmpai[k]*0.5;
	    exc_b[i+ioff] += aos[k+i*bllen]*tmpbi[k]*0.5;
	  }
	}
      }
    }
  }
}

static void
lin_resp_cb_b_gga_ab(DftIntegratorBl* grid, real * RESTRICT tmp,
		     int bllen, int blstart, int blend,
		     LinResp_ab_BlData* data)
{
  static const real MONER = -1.0;
  int ibl, i, jbl, j, k, isym, ivec;
  integer bl_size = 4*DFT_BLLEN;
  real * RESTRICT aos        = grid->atv;
  real * RESTRICT aox        = grid->atv+bllen*inforb_.nbast;
  real * RESTRICT aoy        = grid->atv+bllen*inforb_.nbast*2;
  real * RESTRICT aoz        = grid->atv+bllen*inforb_.nbast*3;
  real (* RESTRICT rhowa)[4] = (real(*)[4])data->rhowa;
  real (* RESTRICT rhowb)[4] = (real(*)[4])data->rhowb;
  real (* RESTRICT vt3a)[4]  = (real(*)[4])data->vta; 
  real (* RESTRICT vt3b)[4]  = (real(*)[4])data->vtb; 
  real * RESTRICT exc_a      = data->res_a;
  real * RESTRICT exc_b      = data->res_b;
  real * RESTRICT tmpa       = data->tmpa;
  real * RESTRICT tmpb       = data->tmpb;
  FunDensProp dp = { 0 };
  FunSecondFuncDrv vxc;
  real znva, znvb;
  real zetaa, zetab, zetac;
  real fac0a, fac0b;
  real facra, facrb, facg;
  integer blen = bllen;

  for(ivec=0; ivec<data->vecs_in_batch; ivec++) {
    /* compute vector of transformed densities and dens. gradients */
    FSYM2(getexp_blocked_gga)(&data->ksymop,data->kappa_a + ivec*inforb_.n2basx,
			      grid->atv, grid->bas_bl_cnt,
			      grid->basblocks, &grid->shl_bl_cnt, tmp, &blen, rhowa);
    FSYM2(getexp_blocked_gga)(&data->ksymop,data->kappa_b + ivec*inforb_.n2basx,
			      grid->atv, grid->bas_bl_cnt,
			      grid->basblocks, &grid->shl_bl_cnt, tmp, &blen, rhowb); 
    if (data->trplet) FSYM(dscal)(&bl_size,&MONER,(real*)rhowb,&ONEI);
    for(i=blstart; i<blend; i++) {          
      /* Compute density parameters at grid point */
      real weight = grid->weight[grid->curr_point+i];
      dp.rhoa = grid->r.ho.a[i];
      dp.rhob = grid->r.ho.b[i];
      dp.grada = sqrt(grid->g.rad.a [i][0]*grid->g.rad.a [i][0]+
                            grid->g.rad.a [i][1]*grid->g.rad.a [i][1]+
		      grid->g.rad.a [i][2]*grid->g.rad.a [i][2]);
      dp.gradb = sqrt(grid->g.rad.b [i][0]*grid->g.rad.b [i][0]+
                            grid->g.rad.b [i][1]*grid->g.rad.b [i][1]+
		      grid->g.rad.b [i][2]*grid->g.rad.b [i][2]);
            dp.gradab = grid->g.rad.a [i][0]*grid->g.rad.b [i][0]+
                        grid->g.rad.a [i][1]*grid->g.rad.b [i][1]+
	      grid->g.rad.a [i][2]*grid->g.rad.b [i][2];
            if(dp.grada<1e-20 || dp.rhoa<1e-20 || dp.gradb<1e-20 || dp.rhob<1e-20) {
	      vt3a[i][0] = vt3a[i][1] = vt3a[i][2] = vt3a[i][3] = 0;
	      vt3b[i][0] = vt3b[i][1] = vt3b[i][2] = vt3b[i][3] = 0;
	      continue; 
            }  
            if(dp.grada<1e-20) dp.grada = 1e-20;
            if(dp.gradb<1e-20) dp.gradb = 1e-20;
            /* second derivatives calculation */
            drv2_clear(&vxc);
            selected_func->second(&vxc,weight,&dp);
            /* compute zetas */
            znva = 1.0/dp.grada;
            znvb = 1.0/dp.gradb;
            zetaa = rhowa[i][1]*grid->g.rad.a [i][0]+
                    rhowa[i][2]*grid->g.rad.a [i][1]+
	            rhowa[i][3]*grid->g.rad.a [i][2];
            zetaa = zetaa*znva;
            zetab = rhowb[i][1]*grid->g.rad.b [i][0]+
                    rhowb[i][2]*grid->g.rad.b [i][1]+
	            rhowb[i][3]*grid->g.rad.b [i][2];
            zetab = zetab*znvb;
            zetac = rhowa[i][1]*grid->g.rad.b [i][0]+
                    rhowa[i][2]*grid->g.rad.b [i][1]+
                    rhowa[i][3]*grid->g.rad.b [i][2]+
                    rhowb[i][1]*grid->g.rad.a [i][0]+
                    rhowb[i][2]*grid->g.rad.a [i][1]+
	            rhowb[i][3]*grid->g.rad.a [i][2];
	    /* rhoa and rhob dependent prefactors */
            fac0a = vxc.df2000*rhowa[i][0]+vxc.df1100*rhowb[i][0]+ 
	            vxc.df1010*zetaa+vxc.df1001*zetab+vxc.df10001*zetac;
            fac0b = vxc.df0200*rhowb[i][0]+vxc.df1100*rhowa[i][0]+
	            vxc.df0101*zetab+vxc.df0110*zetaa+vxc.df01001*zetac;
            vt3a[i][0] = 0.5*fac0a;      
            vt3b[i][0] = 0.5*fac0b;
            /* grada and gradb dependent prefactors*/ 
            facra = vxc.df1010*rhowa[i][0]+vxc.df0110*rhowb[i][0]+ 
	            vxc.df0020*zetaa+vxc.df0011*zetab+vxc.df00101*zetac;
            facrb = vxc.df0101*rhowb[i][0]+vxc.df1001*rhowa[i][0]+ 
	            vxc.df0002*zetab+vxc.df0011*zetaa+vxc.df00011*zetac;
            facra = facra*znva; 
            facrb = facrb*znvb;
            facg  = vxc.df10001*rhowa[i][0]+vxc.df01001*rhowb[i][0]+
	            vxc.df00101*zetaa+vxc.df00011*zetab+vxc.df00002*zetac;
            vt3a[i][1] = grid->g.rad.a[i][0]*facra+ 
                         grid->g.rad.b[i][0]*facg+ 
                         znva*vxc.df0010*rhowa[i][1]-
                         znva*znva*zetaa*vxc.df0010*grid->g.rad.a[i][0]+ 
	                 rhowb[i][1]*vxc.df00001; 
            vt3b[i][1] = grid->g.rad.b[i][0]*facrb+
                         grid->g.rad.a[i][0]*facg+ 
                         znvb*vxc.df0001*rhowb[i][1]-
                         znvb*znvb*zetab*vxc.df0001*grid->g.rad.b[i][0]+ 
	                 rhowa[i][1]*vxc.df00001; 
            vt3a[i][2] = grid->g.rad.a[i][1]*facra+
                         grid->g.rad.b[i][1]*facg+
                         znva*vxc.df0010*rhowa[i][2]- 
                         znva*znva*zetaa*vxc.df0010*grid->g.rad.a[i][1]+         
	                 rhowb[i][2]*vxc.df00001; 
            vt3b[i][2] = grid->g.rad.b[i][1]*facrb+
                         grid->g.rad.a[i][1]*facg+
                         znvb*vxc.df0001*rhowb[i][2]-
                         znvb*znvb*zetab*vxc.df0001*grid->g.rad.b[i][1]+
	                 rhowa[i][2]*vxc.df00001; 
            vt3a[i][3] = grid->g.rad.a[i][2]*facra+ 
                         grid->g.rad.b[i][2]*facg+ 
                         znva*vxc.df0010*rhowa[i][3]- 
                         znva*znva*zetaa*vxc.df0010*grid->g.rad.a[i][2]+
	                 rhowb[i][3]*vxc.df00001; 
            vt3b[i][3] = grid->g.rad.b[i][2]*facrb+
                         grid->g.rad.a[i][2]*facg+
                         znvb*vxc.df0001*rhowb[i][3]- 
                         znvb*znvb*zetab*vxc.df0001*grid->g.rad.b[i][2]+
	                 rhowa[i][3]*vxc.df00001;  
    }

 
    for(isym=0; isym<grid->nsym; isym++) {
      integer (*RESTRICT iblocks)[2] = BASBLOCK(grid,isym);
      int ibl_cnt = grid->bas_bl_cnt[isym];

      for(ibl=0; ibl<ibl_cnt; ibl++) {
	for(i=iblocks[ibl][0]-1; i<iblocks[ibl][1]; i++) { 
	  real *RESTRICT g0i = &aos[i*bllen];
	  real *RESTRICT gxi = &aox[i*bllen];
	  real *RESTRICT gyi = &aoy[i*bllen];
	  real *RESTRICT gzi = &aoz[i*bllen];
	  int ioff = i*inforb_.nbast + ivec*inforb_.n2basx;
	  int jsym = inforb_.muld2h[data->ksymop-1][isym]-1;
	  integer (*RESTRICT jblocks)[2] = BASBLOCK(grid,jsym);
	  int jbl_cnt = grid->bas_bl_cnt[jsym];
	  for(k=blstart; k<blend; k++) {
                      tmpa[k] = vt3a[k][0]*g0i[k]+
                                vt3a[k][1]*gxi[k]+
                                vt3a[k][2]*gyi[k]+
			        vt3a[k][3]*gzi[k]; 
                      tmpb[k] = vt3b[k][0]*g0i[k]+
                                vt3b[k][1]*gxi[k]+
                                vt3b[k][2]*gyi[k]+
		        	vt3b[k][3]*gzi[k];     
	  } 
	  for(jbl=0; jbl<jbl_cnt; jbl++) {
	    int jtop = jblocks[jbl][1];
	    for(j=jblocks[jbl][0]-1; j<jtop; j++) { 
	      real *RESTRICT g0j = &aos[j*bllen];
	      real sa = 0;
	      real sb = 0;
	      for(k=blstart; k<blend; k++) {
		sa += g0j[k]*tmpa[k];
		sb += g0j[k]*tmpb[k];
	      }
	      exc_a[j+ioff] += sa;
	      exc_b[j+ioff] += sb; 
	    }
	  }
	}
      }
    }
  }
}


void
FSYM2(dft_lin_respab_b)(integer *nosim, real* fmatc, real* fmato, real *cmo,
			real *zymat, integer *trplet, integer *ksymop, real* work,
			integer* lwork)
{
  static const int MAX_VEC = 5;
  static const real DP5R  = 0.5;
  static const real MDP5R = -0.5;

  struct tms starttm, endtm; clock_t utm;
  real electrons = 0;
  LinResp_ab_BlData lr_data;
  real *runit;
  int max_vecs;
  int i, j, jvec, ivec;
  real averag;
  real * fmata, * fmatb;

  /* parallel code insert */
#if defined(VAR_MPI)
  dft_wake_slaves((DFTPropEvalMaster)FSYM2(dft_lin_respab_b)); 
  dft_lin_resp_sync_slaves(cmo,nosim,&zymat,trplet,ksymop,
			   &work, *lwork);             
#endif     
  /* start timer and allocate memory */ 
  times(&starttm);
  max_vecs = *nosim>MAX_VEC ? MAX_VEC : *nosim;
  lr_data.dmata   = dal_malloc(2*inforb_.n2basx*sizeof(real));
  lr_data.dmatb   = lr_data.dmata+inforb_.n2basx;
  lr_data.res_a   = dal_malloc(inforb_.n2basx*sizeof(real)*max_vecs);
  lr_data.res_b   = dal_malloc(inforb_.n2basx*sizeof(real)*max_vecs);
  lr_data.kappa_a = dal_malloc(inforb_.n2basx*sizeof(real)*max_vecs);
  lr_data.kappa_b = dal_malloc(inforb_.n2basx*sizeof(real)*max_vecs);
  lr_data.tmpa    = dal_malloc(inforb_.nbast*DFT_BLLEN*sizeof(real));
  lr_data.tmpb    = dal_malloc(inforb_.nbast*DFT_BLLEN*sizeof(real));
  if (selected_func->is_gga()) {
    lr_data.rhowa    = dal_malloc(DFT_BLLEN*4*sizeof(real));
    lr_data.rhowb    = dal_malloc(DFT_BLLEN*4*sizeof(real));
    lr_data.vta      = dal_malloc(4*DFT_BLLEN*sizeof(real));
    lr_data.vtb      = dal_malloc(4*DFT_BLLEN*sizeof(real));  
  } else {
    lr_data.rhowa    = dal_malloc(DFT_BLLEN*sizeof(real));
    lr_data.rhowb    = dal_malloc(DFT_BLLEN*sizeof(real));
    lr_data.vta      = dal_malloc(DFT_BLLEN*sizeof(real));
    lr_data.vtb      = dal_malloc(DFT_BLLEN*sizeof(real)); 
  }
  lr_data.dtgaoa  = dal_malloc(inforb_.nbast*sizeof(real));
  lr_data.dtgaob  = dal_malloc(inforb_.nbast*sizeof(real));
  lr_data.trplet= *trplet;
  lr_data.ksymop= *ksymop;
  /* get density in AO basis*/
  FSYM2(dft_get_ao_dens_matab)(cmo,lr_data.dmatb,lr_data.dmata,work,lwork);
  FSYM(daxpy)(&inforb_.n2basx,&DP5R,lr_data.dmatb,&ONEI,lr_data.dmata,&ONEI);
  FSYM(dscal)(&inforb_.n2basx,&DP5R,lr_data.dmatb,&ONEI);
  /* get active density in MO basis */
  runit=dal_malloc(inforb_.n2ashx*sizeof(real));
  FSYM(dunit)(runit,&inforb_.nasht);
  /* MO alpha, beta Vxc contributions */
  integer sz1 = max_vecs*inforb_.n2orbx;
  int sz1_int= sz1;
  fmata=dal_malloc(sz1_int*sizeof(real));
  fmatb=dal_malloc(sz1_int*sizeof(real));
  for(ivec=0; ivec<*nosim; ivec+=max_vecs) {
    integer sz;
    lr_data.vecs_in_batch = ivec + max_vecs > *nosim ? *nosim - ivec : max_vecs;
    sz = lr_data.vecs_in_batch * inforb_.n2basx;
    FSYM(dzero)(lr_data.kappa_a, &sz);
    FSYM(dzero)(lr_data.kappa_b, &sz);
    for(jvec=0; jvec<lr_data.vecs_in_batch; jvec++)
      FSYM(deq27)(cmo,zymat+(ivec+jvec)*inforb_.n2orbx,runit,
		  lr_data.kappa_b+jvec*inforb_.n2basx,
		  lr_data.kappa_a+jvec*inforb_.n2basx,
		  work,lwork);
    FSYM(daxpy)(&sz,&ONER,lr_data.kappa_b,&ONEI,lr_data.kappa_a,&ONEI);
    FSYM(dzero)(lr_data.res_a, &sz);
    FSYM(dzero)(lr_data.res_b, &sz);
    electrons = dft_integrate_ao_bl(2, lr_data.dmata, work, lwork, 0,
				    (DftBlockCallback)(selected_func->is_gga() ?
				    lin_resp_cb_b_gga_ab : lin_resp_cb_b_lda_ab),
				    &lr_data);
#ifdef VAR_MPI
    dft_lin_respab_b_collect(lr_data.vecs_in_batch, lr_data.res_a, lr_data.res_b, work, *lwork);
    if(fmatc == NULL)  continue;
#endif
    FSYM(dzero)(fmata, &sz1);
    FSYM(dzero)(fmatb, &sz1);
    for(jvec=0; jvec<lr_data.vecs_in_batch; jvec++){
      for(i=0; i<inforb_.nbast; i++) {
	int ioff = i*inforb_.nbast + jvec*inforb_.n2basx;
	for(j=0; j<i; j++) {
	  int joff = j*inforb_.nbast + jvec*inforb_.n2basx;
	  averag = 0.5*(lr_data.res_a[i+joff] +
			lr_data.res_a[j+ioff]);
	  lr_data.res_a[i+joff] = lr_data.res_a[j+ioff] = averag;
	  averag = 0.5*(lr_data.res_b[i+joff] +
			lr_data.res_b[j+ioff]);
	  lr_data.res_b[i+joff] = lr_data.res_b[j+ioff] = averag;
	}
      }
      /* transform to MO Fock matrix contribution  */
      FSYM(lrao2mo)(cmo, &lr_data.ksymop, lr_data.res_a+jvec*inforb_.n2basx,
		    fmata+jvec*inforb_.n2orbx, work, lwork);
      FSYM(lrao2mo)(cmo, &lr_data.ksymop, lr_data.res_b+jvec*inforb_.n2basx,
		    fmatb+jvec*inforb_.n2orbx, work, lwork);
    }
    integer sz2 = inforb_.n2orbx*lr_data.vecs_in_batch;
    FSYM(daxpy)(&sz2,&ONER,fmata,&ONEI,fmatc+ivec*inforb_.n2orbx,&ONEI);
    if (inforb_.nasht > 0) {
      if (*trplet) {
	FSYM(daxpy)(&sz2,&MDP5R,fmatb,&ONEI,fmato+ivec*inforb_.n2orbx,&ONEI);
      } else {
	FSYM(daxpy)(&sz2,&DP5R,fmatb,&ONEI,fmato+ivec*inforb_.n2orbx,&ONEI);
      }
      FSYM(daxpy)(&sz2,&MDP5R,fmata,&ONEI,fmato+ivec*inforb_.n2orbx,&ONEI);
    }
  }
  free(runit);
  free(fmata);
  free(fmatb);
  free(lr_data.dmata);
  free(lr_data.res_a);
  free(lr_data.res_b);
  free(lr_data.kappa_a);
  free(lr_data.kappa_b);
  free(lr_data.rhowa);
  free(lr_data.rhowb);
  free(lr_data.vta);
  free(lr_data.vtb);
  free(lr_data.tmpa);
  free(lr_data.tmpb);
  free(lr_data.dtgaoa);
  free(lr_data.dtgaob);

  times(&endtm);
  utm = endtm.tms_utime-starttm.tms_utime;
  fort_print("      Electrons: %f(%9.3g): LR-AB-DFT*%d evaluation time: %9.1f s",
	     electrons, (double)(electrons-(int)(electrons+0.5)), *nosim,
	     utm/(double)sysconf(_SC_CLK_TCK));
}

