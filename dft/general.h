/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* general.h: general definitions needed by C code:
   function prototypes, often used constants, etc.

   (c) Pawel Salek, pawsa@theochem.kth.se, feb 2002
*/
#ifndef _GENERAL_H_
#define _GENERAL_H_

#include <stdlib.h>

#if !defined(RESTRICT)
#define RESTRICT
#endif

#define ELEMENTS(arr) (sizeof(arr)/sizeof(arr[0]))

/* define the basic floating-point variable type used by the Fortran code */
typedef double real;

typedef struct Functional_ Functional;
extern Functional* selected_func;
typedef struct FirstDrv_  FirstDrv;
typedef struct SecondDrv_ SecondDrv;
typedef struct ThirdDrv_  ThirdDrv;

typedef struct FirstFuncDrv_  FirstFuncDrv;
typedef struct SecondFuncDrv_ SecondFuncDrv;
typedef struct ThirdFuncDrv_  ThirdFuncDrv;

/* Density evaluators */
typedef struct DftDensity_  DftDensity;
typedef struct DftDensProp_ DftDensProp;
typedef struct DftGrid_     DftGrid;

typedef void (*DftDensEvaluator)(DftDensity* dens, DftDensProp* dp,
                                 DftGrid* grid, real* tmp_vec);
struct DftDensity_ {
    DftDensEvaluator evaluate;
    real *dmata, *dmatb;
};

/* DftDensProp structure contains properties of the density that are
   needed for functional evaluation and possibly other purposes.
*/
struct DftDensProp_ {
    real rhoa,  rhob;
    real grada, gradb; /* norms of the density gradient, not squares */
    real gradab;       /* scalar product of grada and gradb */
    /* real current[3] or something may come in the future :-) */
};


void dft_dens_restricted  (DftDensity* dens, DftDensProp* dp, DftGrid* grid,
                           real* tmp_vec);
void dft_dens_unrestricted(DftDensity* dens, DftDensProp* dp, DftGrid* grid,
                           real* tmp_vec);


/* Property evaluators */
void dft_kohn_sham_(real* dmat, real* ksm, real *edfty,
                    real* work, int *lwork, int* ipr);
void dft_lin_resp_(real* fmat, real *cmo, real *zymat, int *trplet, 
		   int *ksymop, real* work,int* lwork);
void dft_mol_grad_(real* dmat, real* work, int* lwork, int* iprint);
void dftqrcf_(real* fi, real* cmo, real* kappaY, int* symY, int* spinY,
              real* kappaZ, int* symZ, int* spinZ, int* addfock,
              real* work, int* lwork);


void dft_kohn_shamab_(real* dmat, real* ksm, real *edfty, 
                      real* work, int *lwork, int* ipr);
void dft_lin_respab_(real* fmatc, real* fmato,  real *cmo, real *zymat, 
                     int *trplet, real* work,int* lwork);

typedef void (*DFTPropEvalMaster)(void);
typedef void (*DFTPropEvalSlave)(real* work, int* lwork, const int* iprint);
#if 0 && defined(VAR_MPI)
#include <mpi.h>
void dft_kohn_sham_slave(real* work, int* lwork, const int* iprint);
void dft_lin_resp_slave (real* work, int* lwork, const int* iprint);
void dft_kohn_shamab_slave(real* work, int* lwork, const int* iprint);
void dft_lin_respab_slave (real* work, int* lwork, const int* iprint);
void dft_mol_grad_slave (real* work, int* lwork, const int* iprint);
void dft_qr_resp_slave  (real* work, int* lwork, const int* iprint);
void dft_wake_slaves(DFTPropEvalMaster);
typedef struct {
    void*        data;
    int          count;
    MPI_Datatype type;
} SyncData;
void mpi_sync_data(const SyncData* data, int count);
#else
#define dft_wake_slaves(a)
#endif

void* dal_malloc(size_t sz);

void fort_print(const char* format, ...);
/* FORTRAN FUNCTION PROTOTYPES */
void dzero_(real* arr, const int* len);
void dunit_(real* arr, const int* len);
void outmat_(const real* mat, const int* rowlow, const int* rowhi,
	     const int* collow, const int* colhi,
             const int* rowdim, const int* coldim);
void getrho_(const real*dmat, const real* atv, real* rho, real* dmagao, 
	     const real* densthr);
void dftgrd_(real* work, int* lwork, const int* d1, const int* log1);
void dftdns_(real* dmat, real* work,int *lwork,int* iprint);
void gtdmso_(real* udv, real* cmo, real* di, real* dv, real* work);
int ishell_cnt_(void);
void dalton_quit(const char* format, ...);

/* BLAS and other linear algebra routines */
real dsum_(const int* cnt, const real* v, const int* stride);
void dscal_(const int* cnt, const real* fac, real* v, const int* stride);
real dnorm2_(const int* cnt, const real* v, const int* stride);
real ddot_(const int* cnt, const real* v1, const int* stride1,
	   const real* v2, const int* stride2);
real dcopy_(const int* cnt, const real* v1, const int* stride1,
	   const real* v2, const int* stride2);
void daxpy_(const int* cnt, const real* alpha, const real* v1, 
	    const int* stride1, const real* v2, const int* stride2);

void dger_(const int* m, const int* n, const real* alpha, 
	   const real* x, const int* incx, const real* y, const int* incy,
	   real* a, const int* lda);

void dgemv_(const char* tr, const int* nrows, const int* ncols, 
	    const real* alpha, const real* a, const int* lda, const real* vb,
	    const int* strideb, const real* beta, 
	    real* c, const int* stridec);

void dgemm_(const char* transa, const char* transb, const int* nrowc, 
	    const int* ncolc, const int* ncolopa, const real* alpha,
	    const real* a, const int* lda, const real* b, const int* ldb,
	    const real* beta, real* c, const int* ldc);

real* alloc_mat_MO(int cnt);

void dft_get_ao_dens_mat_(real* cmo, real* dmat, real* work, int* lwork);
void dft_get_ao_dens_matab_(real* cmo, real* dmata, real* dmatb,
                            real* work, int* lwork);

/* useful  constants for fortran interfacing */
extern const int  ZEROI, ONEI, THREEI, FOURI;
extern const real ZEROR, ONER, TWOR, FOURR;

#if !defined __inline__
/* inline some stuff whenever possible */
#define __inline__
#endif

#endif
