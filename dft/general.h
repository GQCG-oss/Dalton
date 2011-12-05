/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* general.h: general definitions needed by C code:
   function prototypes, often used constants, etc.

   (c) Pawel Salek, pawsa@theochem.kth.se, feb 2002
*/
#ifndef _GENERAL_H_
#define _GENERAL_H_

#include <stdlib.h>

#if defined(VAR_INT64)
#include <stdint.h>
typedef int64_t integer;
#else
typedef int integer;
#endif

#if !defined(RESTRICT)
#define RESTRICT
#endif

/* define the basic floating-point variable type used by the Fortran code */

#if !defined(__CVERSION)
#define __CVERSION__
#endif

#include "functionals.h"

#include "FSYMdef.h"

#if defined(VAR_PGF77)
#define __FUNCTION__ "PGI_does_not_define__FUNCTION__"
#endif
#if defined(SYS_SUN)
#define __FUNCTION__ "SUNs CC compiler_does_not_define__FUNCTION__"
#endif
#if defined(SYS_IRIX)
#define __FUNCTION__ "SGIs CC compiler_does_not_define__FUNCTION__"
#endif
#if defined(SYS_DEC)
#define __FUNCTION__ "DEC CC compiler does not define __FUNCTION__"
#endif

#define ELEMENTS(arr) (sizeof(arr)/sizeof(arr[0]))



/* Density evaluators */
typedef struct DftDensity_  DftDensity;
typedef struct DftGrid_     DftGrid;

typedef void (*DftDensEvaluator)(DftDensity* dens, FunDensProp* dp,
                                 DftGrid* grid, real* tmp_vec);
struct DftDensity_ {
    DftDensEvaluator evaluate;
    real *dmata, *dmatb;
};

/* FirstDrv: matrix of first order derivatives with respect to two
 * parameters: density rho and SQUARE of the gradient of density grho.
 * zeta_i = |\nabla\rho_i|
 * mu     = |\nabla\rho_\alpha||\nabla\rho_\beta|
 */

typedef struct {
    real fR;  /* d/drho F     */
    real fZ;  /* d/zeta F     */
} FirstDrv;

/* SecondDrv:  matrix  of  second  order functional  derivatives  with
 * respect  to two  parameters: density  rho and  SQUARE  of the
 * density gradient zeta.  The derivatives are computed for alpha-alpha
 * or beta-beta spin-orbital block (i.e. include triplet flag).
 */
typedef struct {
    real fR; /* d/drho  F */
    real fZ; /* d/dzeta F */
    real fRR; /* d/drho^2 F */
    real fRZ; /* d/(drho dzeta) F */
    real fZZ; /* d/dzeta^2 F */
    /* additional derivatives required by  */
    /* general linear response scheme     */
    real fRG; /* d/(drho dgamma) F */
    real fZG; /* d/(dzeta dgamma) F */
    real fGG; /* d/dzgamma^2 F */
    real fG;  /* d/dgamma F */
} SecondDrv;

/* ThirdDrv: matrix of third derivatives with respect to two parameters:
   density rho and SQUARE of the density gradient zeta.
*/
typedef struct {
    real fR;   /* d/drho  F */
    real fZ;   /* d/dzeta F */
    real fG;   /* d/dgamma F */
    real fRR[2];  /* d/drho^2 F */
    real fRZ[2];  /* d/(drho dzeta) F */
    real fZZ[2];  /* d/dzeta^2 F */
    real fRG[2];  /* d/(drho dgamma) F */
    real fRRR[2]; /* d/drho^2 F */
    real fRRZ[2][2]; /* d/(drho^2 dzeta) F */
    /* two forms of fRRG needed as the formulae is non symmetric */
    real fRRG[2];     /* d/(drho^2 dgamma) F */
    real fRRGX[2][2]; /* d/(drho^2 dgamma) F */
    real fRZZ[2][2]; /* d/(drho dzeta^2 F */
    real fZZZ[2]; /* d/dzeta^3 F */
} ThirdDrv;


typedef struct {
    real fR;   /* d/drho  F */
    real fZ;   /* d/dzeta F */

    real fRR;  /* d/drho^2 F */
    real fRZ;  /* d/(drho dzeta) F */
    real fZZ;  /* d/dzeta^2 F */

    real fRRR; /* d/drho^3 F */
    real fRRZ; /* d/(drho^2 dzeta) F */
    real fRZZ; /* d/(drho dzeta^2 F */
    real fZZZ;

    real fRRRR;
    real fRRRZ;
    real fRRZZ;
    real fRZZZ;
    real fZZZZ;
} FourthDrv;

void dftpot0_(FirstDrv *ds, const real* weight, const FunDensProp* dp);
void dftpot1_(SecondDrv *ds, const real* w, const FunDensProp* dp,
              const integer* triplet);
void dftpot2_(ThirdDrv *ds, real factor, const FunDensProp* dp, integer isgga,
              integer triplet);
void dftpot3ab_(FourthDrv *ds, const real *factor, const FunDensProp* dp,
                const integer *isgga);


void dft_dens_restricted  (DftDensity* dens, FunDensProp* dp, DftGrid* grid,
                           real* tmp_vec);
void dft_dens_unrestricted(DftDensity* dens, FunDensProp* dp, DftGrid* grid,
                           real* tmp_vec);
extern void
FSYM2(getexp_blocked_lda)(const integer *idsym, real*dmat, const real* atv, 
			  const integer *nblocks, integer (*orbblocks)[2],
			  const integer *lda,
			  real *tmp, const integer * bllen, real *rho);

extern void
FSYM2(getexp_blocked_gga)(const integer *idsym, real*dmat, const real* atv, 
			  const integer *nblocks, integer (*orbblocks)[2],
			  const integer *lda,
			  real *tmp, const integer * bllen, real (*grad)[4]);



/* Property evaluators */
void dft_kohn_sham_(real* dmat, real* ksm, real *edfty,
                    real* work, integer *lwork, integer* ipr);
void dft_lin_resp_(real* fmat, real *cmo, real *zymat, integer *trplet, 
		   integer *ksymop, real* work, integer* lwork);
void FSYM2(dft_lin_respf)(integer *nosim, real* fmat, real *cmo, real *zymat,
                          integer *trplet, integer *ksymop, real* work,
			  integer* lwork);
void dft_mol_grad_(real* dmat, real* work, integer* lwork, integer* iprint);
void dftqrcf_(real* fi, real* cmo, real* kappaY, integer* symY, integer* spinY,
              real* kappaZ, integer* symZ, integer* spinZ, integer* addfock,
              real* work, integer* lwork);
void FSYM2(dft_qr_respons)(real* fi, real* cmo,
                           real* kappaY, integer* symY, integer* spinY,
			   real* kappaZ, integer* symZ, integer* spinZ,
                           integer* addfock, real* work, integer* lwork);
void FSYM(dftcrcf)(real* fi, real* cmo,
		   real* kappaB, integer* symB,
		   real* kappaC, integer* symC,
		   real* kappaD, integer* symD,
		   real* work, integer* lwork);
void FSYM(numdso)(real* spndso, integer *nucind, real* work, integer* lwork);


void FSYM2(dft_kohn_shamab)(real* dmat, real* ksm, real *edfty,
                            real* work, integer *lwork, integer* iprfck);
void FSYM2(dft_lin_respab)(real* fmatc, real* fmato,  real *cmo, real *zymat, 
                           integer *trplet, integer *ksymop,
			   real* work, integer* lwork);
void FSYM(dftmolgradab)(real* work, integer* lwork, integer* iprint);
void FSYM2(dft_kohn_shamab_b)(real* dmat, real* ksm, real *edfty,
			      real* work, integer *lwork, integer* iprfck);
void FSYM2(dft_lin_respab_b)(integer *nosim, real* fmatc, real* fmato, real *cmo,
			     real *zymat, integer *trplet, integer *ksymop, real* work,
			     integer* lwork);

typedef void (*DFTPropEvalMaster)(void);
typedef void (*DFTPropEvalSlave)(real* work, integer* lwork, const integer* iprint);
#if defined(VAR_MPI)
#include <mpi.h>
void dft_kohn_sham_slave(real* work, integer* lwork, const integer* iprint);
void dft_lin_resp_slave (real* work, integer* lwork, const integer* iprint);
void dft_lin_respf_slave (real* work, integer* lwork, const integer* iprint);
void dft_kohn_shamab_slave(real* work, integer* lwork, const integer* iprint);
void dft_lin_respab_slave (real* work, integer* lwork, const integer* iprint);
void dft_mol_grad_slave (real* work, integer* lwork, const integer* iprint);
void dft_qr_resp_slave  (real* work, integer* lwork, const integer* iprint);
void dft_qrbl_slave     (real* work, integer* lwork, const integer* iprint);
void dft_cr_resp_slave(real* work, integer* lwork, const integer* iprint);
void numdso_slave(real* work, integer* lwork, const integer* iprint);
void dft_kohn_shamab_b_slave(real* work, integer* lwork, const integer* iprint);
void dft_lin_respab_b_slave(real* work, integer* lwork, const integer* iprint);
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

void* dal_malloc_(size_t sz, const char *func, unsigned line);
#define dal_malloc(sz) dal_malloc_((sz),__FUNCTION__, __LINE__)
#define dal_new(sz,tp) (tp*)dal_malloc_((sz)*sizeof(tp),__FUNCTION__, __LINE__)

int fort_print(const char* format, ...);
/* FORTRAN FUNCTION PROTOTYPES */
void dzero_(real* arr, const integer* len);
void dunit_(real* arr, const integer* len);
void outmat_(const real* mat, const integer* rowlow, const integer* rowhi,
	     const integer* collow, const integer* colhi,
             const integer* rowdim, const integer* coldim);
void getrho_(const real*dmat, const real* atv, real* rho, real* dmagao, 
	     const real* densthr);
void dftgrd_(real* work, integer* lwork, const integer* d1, const integer* log1);
void dftdns_(real* dmat, real* work,integer *lwork, const integer* iprint);
void gtdmso_(real* udv, real* cmo, real* di, real* dv, real* work);
void dftdnsab_(real* dmata,real* dmatb, real* work, integer* lwork, integer* iprint);
void udftmolgrdab_(real* gao, real* damta, real* dmatb, real* rha, real* rhb, 
                   real* vra, real* vrb, real* vza, real* vzb, real* vzg); 
integer FSYM2(ishell_cnt)(void);
void dalton_quit(const char* format, ...);

/* BLAS and other linear algebra routines */
real dsum_(const integer* cnt, const real* v, const integer* stride);
void dscal_(const integer* cnt, const real* fac, real* v, const integer* stride);
real dnorm2_(const integer* cnt, const real* v, const integer* stride);
real ddot_(const integer* cnt, const real* v1, const integer* stride1,
	   const real* v2, const integer* stride2);
real dcopy_(const integer* cnt, const real* v1, const integer* stride1,
	   const real* v2, const integer* stride2);
void daxpy_(const integer* cnt, const real* alpha, const real* v1, 
	    const integer* stride1, const real* v2, const integer* stride2);

void dger_(const integer* m, const integer* n, const real* alpha, 
	   const real* x, const integer* incx, const real* y, const integer* incy,
	   real* a, const integer* lda);

void dgemv_(const char* tr, const integer* nrows, const integer* ncols, 
	    const real* alpha, const real* a, const integer* lda,
	    const real* vb, const integer* strideb, const real* beta, 
	    real* c, const integer* stridec);

void dgemm_(const char* transa, const char* transb, const integer* nrowc, 
	    const integer* ncolc, const integer* ncolopa, const real* alpha,
	    const real* a, const integer* lda, const real* b, const integer* ldb,
	    const real* beta, real* c, const integer* ldc);

real* alloc_mat_MO(int cnt);

void FSYM2(dft_get_ao_dens_mat)(const real* cmo, real* dmat,
                                real* work, integer* lwork);
void FSYM2(dft_get_ao_dens_matab)(real* cmo, real* dmata, real* dmatb,
                                  real* work, integer* lwork);
void FSYM(lrao2mo)(const real* cmo, const integer *ksymop, 
                   const real*res, real* fmat, real* work, integer*lw);

/* useful  constants for fortran interfacing */
extern const integer ZEROI, ONEI, THREEI, FOURI;
extern const real ZEROR, ONER, TWOR, FOURR;

#if !defined __inline__
/* inline some stuff whenever possible */
#define __inline__
#endif
#define CHECK_WRKMEM(req,lwork)\
 if((req)>(lwork)){dalton_quit("%s requires %u words but only %u available\n",\
                   __FUNCTION__, (unsigned)(req), (unsigned)(lwork));}

#endif
