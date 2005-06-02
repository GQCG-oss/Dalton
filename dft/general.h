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

/* define the basic floating-point variable type used by the Fortran code */

#if !defined(__CVERSION)
#define __CVERSION__
#endif

#include "functionals.h"

/* Match Fortran name mangling. If the Fortran compiler does not
 * mangle names, define NO_UNDERSCORE in CFLAGS.  g77 and compaq fort
 * (cryptically referred to with HAVE_GCPP below) for linux-alpha both
 * insert a second underscore if routine name contains at least one
 * underscore /hjaaj Oct04 */
#ifdef NO_UNDERSCORE
#define FSYM(a) a
#define FSYM2(a) a
#else
#define FSYM(a) a ## _
#if defined(VAR_G77) || defined(HAVE_GCPP)
#define FSYM2(a) a ## __
#else
#define FSYM2(a) a ## _
#endif
#endif

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
 * zeta_i = |\nabla\rho_i|%Gï¿¿%@
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
    real fRR; /* d/drho%Gï¿¿%@ F */
    real fRZ; /* d/(drho dzeta) F */
    real fZZ; /* d/dzeta%Gï¿¿%@ F */
    /* additional derivatives required by  */
    /* general linear response scheme     */
    real fRG; /* d/(drho dgamma) F */
    real fZG; /* d/(dzeta dgamma) F */
    real fGG; /* d/dzgamma%Gï¿¿%@ F */
    real fG;  /* d/dgamma F */
} SecondDrv;

/* ThirdDrv: matrix of third derivatives with respect to two parameters:
   density rho and SQUARE of the density gradient zeta.
*/
typedef struct {
    real fR;   /* d/drho  F */
    real fZ;   /* d/dzeta F */
    real fG;   /* d/dgamma F */
    real fRR[2];  /* d/drho%Gï¿¿%@ F */
    real fRZ[2];  /* d/(drho dzeta) F */
    real fZZ[2];  /* d/dzeta%Gï¿¿%@ F */
    real fRG[2];  /* d/(drho dgamma) F */
    real fRRR[2]; /* d/drho%Gï¿¿%@ F */
    real fRRZ[2][2]; /* d/(drho%Gï¿¿%@ dzeta) F */
    /* two forms of fRRG needed as the formulae is non symmetric */
    real fRRG[2];     /* d/(drho? dgamma) F */
    real fRRGX[2][2]; /* d/(drho? dgamma) F */
    real fRZZ[2][2]; /* d/(drho dzeta%Gï¿¿%@) F */
    real fZZZ[2]; /* d/dzeta%Gï¿¿%@ F */
} ThirdDrv;


typedef struct {
    real fR;   /* d/drho  F */
    real fZ;   /* d/dzeta F */

    real fRR;  /* d/drho%G￿%@ F */
    real fRZ;  /* d/(drho dzeta) F */
    real fZZ;  /* d/dzeta%G￿%@ F */

    real fRRR; /* d/drho%G￿%@ F */
    real fRRZ; /* d/(drho%G￿%@ dzeta) F */
    real fRZZ; /* d/(drho dzeta%G￿%@) F */
    real fZZZ;

    real fRRRR;
    real fRRRZ;
    real fRRZZ;
    real fRZZZ;
    real fZZZZ;
} FourthDrv;

void dftpot0_(FirstDrv *ds, const real* weight, const FunDensProp* dp);
void dftpot1_(SecondDrv *ds, const real* w, const FunDensProp* dp,
              const int* triplet);
void dftpot2_(ThirdDrv *ds, real factor, const FunDensProp* dp, int isgga,
              int triplet);
void dftpot3ab_(FourthDrv *ds, real *factor, const FunDensProp* dp,
                int *isgga);


void dft_dens_restricted  (DftDensity* dens, FunDensProp* dp, DftGrid* grid,
                           real* tmp_vec);
void dft_dens_unrestricted(DftDensity* dens, FunDensProp* dp, DftGrid* grid,
                           real* tmp_vec);
extern void
FSYM2(getexp_blocked_lda)(const int *idsym, real*dmat, const real* atv, 
                    const int *nblocks, int (*orbblocks)[2], const int *lda,
                    real *tmp, const int * bllen, real *rho);

extern void
FSYM2(getexp_blocked_gga)(const int *idsym, real*dmat, const real* atv, 
                    const int *nblocks, int (*orbblocks)[2], const int *lda,
                    real *tmp, const int * bllen, real (*grad)[4]);



/* Property evaluators */
void dft_kohn_sham_(real* dmat, real* ksm, real *edfty,
                    real* work, int *lwork, int* ipr);
void dft_lin_resp_(real* fmat, real *cmo, real *zymat, int *trplet, 
		   int *ksymop, real* work,int* lwork);
void FSYM2(dft_lin_respf)(int *nosim, real* fmat, real *cmo, real *zymat,
                          int *trplet, int *ksymop, real* work,int* lwork);
void dft_mol_grad_(real* dmat, real* work, int* lwork, int* iprint);
void dftqrcf_(real* fi, real* cmo, real* kappaY, int* symY, int* spinY,
              real* kappaZ, int* symZ, int* spinZ, int* addfock,
              real* work, int* lwork);
void FSYM2(dft_qr_respons)(real* fi, real* cmo,
                           real* kappaY, int* symY, int* spinY,
			   real* kappaZ, int* symZ, int* spinZ,
                           int* addfock, real* work, int* lwork);
void dftcrcf_(real* fi, real* cmo,
              real* kappaB, int* symB,
              real* kappaC, int* symC,
              real* kappaD, int* symD,
              real* work, int* lwork);


void FSYM2(dft_kohn_shamab)(real* dmat, real* ksm, real *edfty,
                            real* work, int *lwork, int* iprfck);
void FSYM2(dft_lin_respab)(real* fmatc, real* fmato,  real *cmo, real *zymat, 
                           int *trplet, int *ksymop, real* work,int* lwork);
void dftmolgradab_(real* work, int* lwork, int* iprint);
typedef void (*DFTPropEvalMaster)(void);
typedef void (*DFTPropEvalSlave)(real* work, int* lwork, const int* iprint);
#if defined(VAR_MPI)
#include <mpi.h>
void dft_kohn_sham_slave(real* work, int* lwork, const int* iprint);
void dft_lin_resp_slave (real* work, int* lwork, const int* iprint);
void dft_lin_respf_slave (real* work, int* lwork, const int* iprint);
void dft_kohn_shamab_slave(real* work, int* lwork, const int* iprint);
void dft_lin_respab_slave (real* work, int* lwork, const int* iprint);
void dft_mol_grad_slave (real* work, int* lwork, const int* iprint);
void dft_qr_resp_slave  (real* work, int* lwork, const int* iprint);
void dft_qrbl_slave     (real* work, int* lwork, const int* iprint);
void dft_cr_resp_slave(real* work, int* lwork, const int* iprint);

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

int fort_print(const char* format, ...);
/* FORTRAN FUNCTION PROTOTYPES */
void dzero_(real* arr, const int* len);
void dunit_(real* arr, const int* len);
void outmat_(const real* mat, const int* rowlow, const int* rowhi,
	     const int* collow, const int* colhi,
             const int* rowdim, const int* coldim);
void getrho_(const real*dmat, const real* atv, real* rho, real* dmagao, 
	     const real* densthr);
void dftgrd_(real* work, int* lwork, const int* d1, const int* log1);
void dftdns_(real* dmat, real* work,int *lwork, const int* iprint);
void gtdmso_(real* udv, real* cmo, real* di, real* dv, real* work);
void dftdnsab_(real* dmata,real* dmatb, real* work, int* lwork, int* iprint);
void udftmolgrdab_(real* gao, real* damta, real* dmatb, real* rha, real* rhb, 
                   real* vra, real* vrb, real* vza, real* vzb, real* vzg); 
int FSYM2(ishell_cnt)(void);
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

void FSYM2(dft_get_ao_dens_mat)(const real* cmo, real* dmat,
                                real* work, int* lwork);
void FSYM2(dft_get_ao_dens_matab)(real* cmo, real* dmata, real* dmatb,
                                  real* work, int* lwork);

/* useful  constants for fortran interfacing */
extern const int  ZEROI, ONEI, THREEI, FOURI;
extern const real ZEROR, ONER, TWOR, FOURR;

#if !defined __inline__
/* inline some stuff whenever possible */
#define __inline__
#endif
#define CHECK_WRKMEM(req,lwork)\
 if((req)>(lwork)){dalton_quit("%s requires %u words but only %u available\n",\
                   __FUNCTION__, (unsigned)(req), (unsigned)(lwork));}

#endif
