/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* general.c:
   (c) Pawel Salek, pawsa@theochem.kth.se, 2001-08-02
   NOTES: Adding new functionals:
   a. use fun-dirac.c as template.
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
    char func_name[20];
    int i;
    for(i=len-1; i>=0 && isspace((int)line[i]); i--)
        ;
    if(DftConfString) free(DftConfString);
    i++;
    DftConfString = malloc(i+1);
    strncpy(DftConfString, line, i); 
    DftConfString[i] = '\0';
    sscanf(DftConfString,"%20s", func_name);

    for(i=0; available_functionals[i]; i++)
        if(strcasecmp(available_functionals[i]->name, func_name)==0) {
            int ok;
            selected_func = available_functionals[i];
            ok = selected_func->read ?
                selected_func->read(DftConfString+strlen(func_name)) : 1;
            if(!ok) (*inperr)++;
            return ok;
        }
    fort_print("Unknown functional '%s'. Aborting.\n", func_name);
    dftlistfuncs_();
    (*inperr)++;
    return 0; /* failed to find */
}

/* =================================================================== */
/*    DFT density evaluators for restricted and unrestricted cases.    */
/*    evaluate density properties                                      */
/* =================================================================== */
void
dft_dens_restricted(DftDensity* dens, DftDensProp* dp, DftGrid* grid,
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
               &grid->atv[inforb_.nbast-1], &inforb_.nbast, tmp_vec,
               &ONEI, &ZEROR, grid->grada, &ONEI);
        
        ngrad = sqrt(grid->grada[0]*grid->grada[0]+
                     grid->grada[1]*grid->grada[1]+
                     grid->grada[2]*grid->grada[2]);
        dp->grada = dp->gradb = ngrad;
        dp->gradab= dp->grada*dp->gradb;
    }
}

void
dft_dens_unrestricted(DftDensity* dens, DftDensProp* dp, DftGrid* grid,
                      real* tmp_vec)
{
    real *tmpa = tmp_vec, *tmpb = tmp_vec + inforb_.nbast;
    getrho_(dens->dmata, grid->atv, &dp->rhoa, tmpa, &grid->dfthri);
    getrho_(dens->dmatb, grid->atv, &dp->rhob, tmpb, &grid->dfthri);

    if(dp->rhoa==0) dp->rhoa = 1e-40;
    if(dp->rhob==0) dp->rhob = dp->rhoa*1e-15;

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
    { (DFTPropEvalMaster)dft_lin_respf_,   dft_lin_respf_slave   },
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
FSYM(dft_cslave)(real* work,int*lwork,int*iprint)
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
FSYM(dft_cslave)(real* work,int*lwork,int*iprint)
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
extern void fort_wrt_(const char* str, const int* len, int ln);
void
fort_print(const char* format, ...)
{
    char line[128];
    int len;
    va_list a;

    va_start(a, format);
    vsnprintf(line, sizeof(line), format, a);
    va_end(a);
    len = strlen(line);
    fort_wrt_(line, &len, len);
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
 
