/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* general.c:
   (c) Pawel Salek, pawsa@theochem.kth.se, 2001-08-02
   NOTES: Adding new functionals:
   a. use fun-slater.c as template.
   b. add 'extern Functional MyFunctional;' to functionals.h
   c. add '&MyFunctional' to available_functionals below.
   d. have a beer. Or some crackers, if you prefer.
*/

/* strictly conform to XOPEN ANSI C standard */
#if !defined(SYS_DEC)
/* XOPEN compliance is missing on old Tru64 4.0E Alphas */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1
#endif

/* Use BSD's strncasecmp(); if there is a platform that has no strncasecmp()
 * ask pawsa@theochem.kth.se for replacement */
#define _BSD_SOURCE 1

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#define __CVERSION__

#include "general.h"
#include "integrator.h"
#include "functionals.h"
#include "inforb.h"

/* C-wide constants */
const int  ZEROI = 0,   ONEI = 1, THREEI = 3, FOURI = 4;
const real ZEROR = 0.0, ONER = 1.0, TWOR = 2.0, FOURR = 4.0;


/* stub subroutines for the functional code */
extern real dftgethf_(void);
extern void dftsethf_(real *w);
extern void dftsetcam_(real *w, real *b);

static real
dal_get_hf_weight(void) { return dftgethf_(); }
static void
dal_set_hf_weight(real w) { dftsethf_(&w); }
static void
dal_set_cam_param(real w, real be) { dftsetcam_(&w, &be);}

/* =================================================================== */
/* dftinput:

   read DFT functional from given line. The calling convention assumes
   Sun-style passing parameters. ATLAS linear algebra package
   http://www.netlib.org/atlas/ or http://math-atlas.sourceforge.net/
   contains an elaborate discuttion of character type variable passing
   conventions, in particular little bit below
   http://math-atlas.sourceforge.net/errata.html#RH7.0
*/
static char* DftConfString = NULL;
int
FSYM(dftsetfunc)(const char* line, int * inperr, int len)
{
    int i, off;

    /* set the functional code printf function and HF weight setting
       functions to the dalton version that appends the output to the
       DALTON.OUT file. */
    fun_printf        = fort_print;
    fun_set_hf_weight = dal_set_hf_weight;
    fun_get_hf_weight = dal_get_hf_weight;
    fun_set_cam_param = dal_set_cam_param;

    for(i=len-1; i>=0 && isspace((int)line[i]); i--)
        ;
    if(DftConfString) free(DftConfString);
    i++;
    for(off=0; line[off] && isspace((int)line[off]); off++)
        ;
    DftConfString = malloc(i+1-off);
    strncpy(DftConfString, line+off, i-off); 
    DftConfString[i-off] = '\0';
    
    switch(fun_select_by_name(DftConfString)) {
    case FUN_OK: return 1; /* SUCCESS! */
    case FUN_UNKNOWN:
        fort_print("Unknown functional '%s'. Aborting.\n", DftConfString);
        dftlistfuncs_();
        break;
    case FUN_CONF_ERROR:
        fort_print("Functional configuration '%s' is not understood. "
                   "Aborting.\n", DftConfString);
        break;
    }
    (*inperr)++;
    return 0; /* failed to find */
}

/* =================================================================== */
/* transformations of functional derivatives from from unrestricted
 * ones to closed shell case. */
real
dftene_(const real *rho, const real *grad)
{
    FunDensProp dp;
    dp.rhoa  = dp.rhob  = *rho  *0.5;
    dp.grada = dp.gradb = *grad *0.5;
    dp.gradab = dp.grada*dp.gradb;
    return selected_func->func(&dp);
}

void
dftptf0_(real *rho, real *grad, real *wght, real *vx)
{
    FunFirstFuncDrv drvs;
    FunDensProp dp;
    dp.rhoa  = dp.rhob  = *rho *0.5;
    dp.grada = dp.gradb = *grad*0.5;
    if(dp.rhoa<1e-13) dp.rhoa = dp.rhob = 1e-13;
    if(dp.grada<1e-13) dp.grada = dp.gradb = 1e-13;
    dp.gradab = dp.grada*dp.gradb;
    drv1_clear(&drvs);
    selected_func->first(&drvs, *wght, &dp);
    vx[0] = drvs.df1000;
    vx[1] = drvs.df0010 + 0.5*drvs.df00001* (*grad);
}

void
dftpot0_(FirstDrv *ds, const real* weight, const FunDensProp* dp)
{
    FunFirstFuncDrv drvs;
    drv1_clear(&drvs);
    selected_func->first(&drvs, *weight, dp);
    ds->fR = drvs.df1000;
    ds->fZ = drvs.df0010 + 0.5*drvs.df00001* (dp->grada+dp->gradb);
}

void
dftpot1_(SecondDrv *ds, const real* w, const FunDensProp* dp, 
         const int* triplet)
{
    
    FunSecondFuncDrv drvs;

    drv2_clear(&drvs);
    if(dp->rhoa + dp->rhob>1e-14)
        selected_func->second(&drvs, *w, dp);
    if (*triplet) { /* triplet */  
        ds->fZ  = drvs.df0010;
        ds->fG  = -0.5*drvs.df00001;
        ds->fRR = 0.5*(drvs.df2000 - drvs.df1100);
        ds->fRZ = 0.5*(drvs.df1010 - drvs.df1001); 
        ds->fRG = 0.0;  
        ds->fZZ = 0.5*(drvs.df0020 - drvs.df0011); 
        ds->fZG = 0.0; 
        ds->fGG = 0.0; 
    } else { /* singlet */
        ds->fR  = 0.5*(drvs.df1000 + drvs.df0100);
        ds->fZ  = drvs.df0010;
        ds->fRR = 0.5*(drvs.df2000 + drvs.df1100);
        ds->fRZ = 0.5*(drvs.df1010 + drvs.df1001);
        ds->fZZ = 0.5*(drvs.df0020 + drvs.df0011); 
        ds->fRG = 0.5*drvs.df10001;   
        ds->fZG = 0.5*drvs.df00101;   
        ds->fGG = 0.25*drvs.df00002; 
        ds->fG  = 0.5*drvs.df00001;  
    }
}



/* dftpot2_:
   computes third order derivatives of selected functional with respect
   to rho and zeta=|\nabla\rho|
*/
void
dftpot2_(ThirdDrv *ds, real factor, const FunDensProp* dp, int isgga,
         int triplet)
{
    FunThirdFuncDrv drvs;

    drv3_clear(&drvs);
    selected_func->third(&drvs, factor, dp);
    /* transform to density from density_alpha derivatives below */
    /* This could be a separate function. */
    /* we treat singlet here */
    ds->fR   = (drvs.df1000);
    ds->fRR[0] = (drvs.df2000 + drvs.df1100);
    ds->fRR[1] = (drvs.df2000 - drvs.df1100);
    ds->fRRR[0] = (drvs.df3000 + 3*drvs.df2100);
    ds->fRRR[1] = (drvs.df3000 - drvs.df2100);

    if(isgga) { /* FORMULAE seem ok but were never really tested */
        real grada = dp->grada;
        real grada2= grada*grada;
        real grada3= grada2*grada;
	/* transform to closed shell, second time. Oh, I love this mess. */
        ds->fZ  = drvs.df0010/(2*grada);
        ds->fG = drvs.df00001;
        ds->fZZ[0]  = (drvs.df0020 + drvs.df0011)/(4*grada2) 
	    -drvs.df0010/(4*grada3);
        ds->fZZ[1]  = (drvs.df0020 - drvs.df0011)/(4*grada2) 
	    -drvs.df0010/(4*grada3);
        ds->fRZ[0]  = (drvs.df1010 + drvs.df1001)/(2*grada);
        ds->fRZ[1]  = (drvs.df1010 - drvs.df1001)/(2*grada);
        ds->fRG[0]  = 2*drvs.df10001;
        ds->fRG[1]  = 0.0;  
        ds->fRRZ[0][0] = (drvs.df2010+drvs.df2001+2*drvs.df1110)/(2*grada);
        ds->fRRZ[0][1] = (drvs.df2010+drvs.df2001-2*drvs.df1110)/(2*grada);
        ds->fRRZ[1][0] = (drvs.df2010-drvs.df2001)/(2*grada);
        ds->fRRZ[1][1] = ds->fRRZ[1][0];
        ds->fRRG[0] = drvs.df20001+drvs.df11001;
        ds->fRRG[1] = drvs.df20001-drvs.df11001;
        ds->fRRGX[0][0] = 2*(drvs.df20001+drvs.df11001);
        ds->fRRGX[1][1] = 2*(drvs.df20001-drvs.df11001);
        ds->fRRGX[1][0] = ds->fRRGX[0][1] = 0;
         
        ds->fRZZ[0][0] = (drvs.df1020+drvs.df0120+2*drvs.df1011)/(4*grada2) 
                      - (drvs.df1010+drvs.df1001)/(4*grada3);
        ds->fRZZ[0][1] = (drvs.df1020+drvs.df0120-2*drvs.df1011)/(4*grada2) 
                      - (drvs.df1010+drvs.df1001)/(4*grada3);
        ds->fRZZ[1][0] = (drvs.df1020-drvs.df0120)/(4*grada2) 
                      - (drvs.df1010-drvs.df1001)/(4*grada3);
        ds->fRZZ[1][1] = ds->fRZZ[1][0];
        ds->fZZZ[0] = ((drvs.df0030 + 3*drvs.df0021)/grada3 
                   -3*(drvs.df0020 + drvs.df0011)/(grada2*grada2)
                   +3*drvs.df0010/(grada3*grada2))/8.0; 
        ds->fZZZ[1] = ((drvs.df0030 - drvs.df0021)/grada3 
                   -(3*drvs.df0020 - drvs.df0011)/(grada2*grada2)
                   +3*drvs.df0010/(grada3*grada2))/8.0; 
    } else {
        ds->fZ = ds->fZZ[0] = ds->fZZ[1] = ds->fRZ[0] = ds->fRZ[1] = 0;
        ds->fRRZ[0][0] = ds->fRRZ[0][1]  = ds->fRRZ[1][0] = 0;
        ds->fRRZ[1][1] = ds->fRZZ[0][0] = ds->fRZZ[0][1] = 0;
        ds->fRZZ[1][0] = ds->fRZZ[1][1] = ds->fZZZ[0] = ds->fZZZ[1] = 0; 
        ds->fG = ds->fRG[0] = ds->fRG[1]= ds->fRRG[0] = ds->fRRG[1] = 0;
        ds->fRRGX[0][0] = ds->fRRGX[1][0] = 
	    ds->fRRGX[0][1] = ds->fRRGX[1][1] = 0;  
    }
}


/* =================================================================== */
/*    DFT density evaluators for restricted and unrestricted cases.    */
/*    evaluate density properties                                      */
/* =================================================================== */
void
dft_dens_restricted(DftDensity* dens, FunDensProp* dp, DftGrid* grid,
                    real* tmp_vec)
{
    /* Note that dens->dmata is really the total density, i.e rho_a+rho_b */
    /* since this is the convention in the rest of dalton */
    real rho, ngrad;
    getrho_(dens->dmata, grid->atv, &rho, tmp_vec, &grid->dfthri);
    dp->rhoa = dp->rhob = 0.5*rho; 
    /* transform grad vectors to molecular orbitals */
    if(rho>grid->dfthr0 && grid->dogga) {
        /* compute only half density gradient, i.e only grad_alpha. */
        dgemv_("T", &inforb_.nbast, &THREEI, &ONER,
               &grid->atv[inforb_.nbast], &inforb_.nbast, tmp_vec,
               &ONEI, &ZEROR, grid->grada, &ONEI);
        
        ngrad = sqrt(grid->grada[0]*grid->grada[0]+
                     grid->grada[1]*grid->grada[1]+
                     grid->grada[2]*grid->grada[2]);
        dp->grada = dp->gradb = ngrad;
        dp->gradab= dp->grada*dp->gradb;
    }
}

void
dft_dens_unrestricted(DftDensity* dens, FunDensProp* dp, DftGrid* grid,
                      real* tmp_vec)
{
    real *tmpa = tmp_vec, *tmpb = tmp_vec + inforb_.nbast;
    getrho_(dens->dmata, grid->atv, &dp->rhoa, tmpa, &grid->dfthri);
    getrho_(dens->dmatb, grid->atv, &dp->rhob, tmpb, &grid->dfthri);

    if(dp->rhoa<1e-40)          dp->rhoa = 1e-40;
    if(dp->rhob<dp->rhoa*1e-15) dp->rhob = dp->rhoa*1e-15;

    if( (dp->rhoa+dp->rhob)>grid->dfthr0 && grid->dogga) {
        /* transform grad vectors to molecular orbitals */
        dgemv_("T", &inforb_.nbast, &THREEI, &TWOR,
               &grid->atv[inforb_.nbast], &inforb_.nbast, tmpa,
               &ONEI, &ZEROR, grid->grada, &ONEI);
        dp->grada = sqrt(grid->grada[0]*grid->grada[0]+
                         grid->grada[1]*grid->grada[1]+
                         grid->grada[2]*grid->grada[2]);
        dgemv_("T", &inforb_.nbast, &THREEI, &TWOR,
               &grid->atv[inforb_.nbast], &inforb_.nbast, tmpb,
               &ONEI, &ZEROR, grid->gradb, &ONEI);
        dp->gradb = sqrt(grid->gradb[0]*grid->gradb[0]+
                         grid->gradb[1]*grid->gradb[1]+
                         grid->gradb[2]*grid->gradb[2]);

        dp->gradab = grid->grada[0]*grid->gradb[0]+
            grid->grada[1]*grid->gradb[1]+grid->grada[2]*grid->gradb[2];
  }
}

#if defined(VAR_MPI)
#include <mpi.h>
#define MASTER_NO 0
/* =================================================================== */
/* General parallel routines.
 * dft_lifesupport_sync - one-time sync of basic common block data
 *                   that is crucial for initialization.
 * dft_wake_slaves - wakes slaves for a specified DFT evaluator (perhaps this
 *                   can be generalized?). It gets a registered property
 *                   evaluator address as an argument and broadcasts a
 *                   message that will call corresponding slave code
 *                   via dft_cslave().
 * FIXME: use rather a register approach that does not require modification
 * of this file every time new property evaluator is added.
 */

struct {
    DFTPropEvalMaster master_func;
    DFTPropEvalSlave  slave_func;
} PropEvaluatorList[] = {
#if 0
    { (DFTPropEvalMaster)dft_kohn_sham_,   dft_kohn_sham_slave   },
    { (DFTPropEvalMaster)dft_kohn_shamab_, dft_kohn_shamab_slave },
    { (DFTPropEvalMaster)dft_lin_resp_,    dft_lin_resp_slave    },
#endif
    { (DFTPropEvalMaster)FSYM2(dft_lin_respf),   dft_lin_respf_slave   },
#if 0
    { (DFTPropEvalMaster)dft_lin_respab_,  dft_lin_respab_slave  },
    { (DFTPropEvalMaster)dft_mol_grad_,    dft_mol_grad_slave    },
#endif
    { (DFTPropEvalMaster)dftqrcf_,         dft_qr_resp_slave     }
};

/* mpi_sync_data:
   sync given set of data with slaves.
*/
void
mpi_sync_data(const SyncData* data, int count)
{
    int i;
    for(i=0; i<count; i++) {
        MPI_Bcast(data[i].data, data[i].count, data[i].type,
                  0, MPI_COMM_WORLD);
    }
}


void
dft_wake_slaves(DFTPropEvalMaster evaluator)
{
    static int iprtyp = 5; /* magic DFT/C number */
    static int iprint = 0;
    int id, mynum;

    MPI_Comm_rank(MPI_COMM_WORLD, &mynum);
    if(mynum != 0)
        return; /* slaves do not wake up other slaves */

    for(id=0; 
        id<ELEMENTS(PropEvaluatorList) && 
            PropEvaluatorList[id].master_func != evaluator;
        id++)
        ;
        
    if(id>=ELEMENTS(PropEvaluatorList)) {
        /* this would really be an programming error.... */
        fprintf(stderr, "Evaluator not registered. No slaves activated.\n");
        return;
    }
    /* ignore MPI errors */
    MPI_Bcast(&iprtyp,1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&iprint,1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&id,    1, MPI_INT, 0, MPI_COMM_WORLD);
    FSYM(dftintbcast)();
}

/* dft_cslave:
   a slave task handler. Receives the task ID and calls apropriate
   property evaluator.
   We broadcast also some basic properties and common block data
   that slaves should absolutely know but for some reason do not.
*/
void
FSYM2(dft_cslave)(real* work,int*lwork,int*iprint)
{
    int rank, size;
    if(MPI_Comm_rank(MPI_COMM_WORLD, &rank) ||
       MPI_Comm_size(MPI_COMM_WORLD, &size)) printf("MPI error\n");
    else {
        int id;
        MPI_Bcast(&id,1,MPI_INT, MASTER_NO, MPI_COMM_WORLD);
        FSYM(dftintbcast)();
        (PropEvaluatorList[id].slave_func)(work, lwork, iprint);
    }
}
#else
void
FSYM2(dft_cslave)(real* work,int*lwork,int*iprint)
{
   fort_print("DFT slave called but does nothing now.");
}
#endif /* VAR_MPI */


extern void quit_(const char* str, int len);
void
dalton_quit(const char* format, ...)
{
    char line[128];
    int len;
    va_list a;
 
    va_start(a, format);
    vsnprintf(line, sizeof(line), format, a);
    va_end(a);
    len = strlen(line);
    quit_(line, len);
}

/* Helper functions. Could be bracketed with #ifdef DEBUG or something */
extern void FSYM2(fort_wrt)(const char* str, const int* len, int ln);
int
fort_print(const char* format, ...)
{
    char line[128];
    int len;
    va_list a;

    va_start(a, format);
    vsnprintf(line, sizeof(line), format, a);
    va_end(a);
    len = strlen(line);
    FSYM2(fort_wrt)(line, &len, len);
    return len;
}

/* ===================================================================
 * Parallel support.
 * =================================================================== */
#ifdef VAR_MPI
#include <mpi.h>

/* dftfuncsync synchronises the selected functional between all
 * nodes. */
void
FSYM(dftfuncsync)(int *mynum, int *nodes)
{
    static int done = 0;
    int len = DftConfString ? strlen(DftConfString) + 1: 0;
    if(done) return;
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(len>0) {
        int res;
        char *line = malloc(len);
        if(*mynum == 0)
            strcpy(line, DftConfString);
        MPI_Bcast(line, len, MPI_CHAR, 0, MPI_COMM_WORLD);
        if(*mynum != 0) {
            int res;
            FSYM(dftsetfunc)(line, &res, len);
        }
        free(line);
    }
    done = 1;
}

#endif
 
