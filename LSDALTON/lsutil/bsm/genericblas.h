/* general.h should be included for Failure.h */

#if defined(HAVE_LIBESSL)
/* this is just a hack, but AIX is the only system with non-standard
 * symbol mangling */
#define F(x) x
#else
#define F(x) x ## _
#endif

#if defined(VAR_INT64)
/* #include <stdint.h> */
/* typedef int64_t integer; stinne comment*/
typedef long integer;
#else
typedef int integer;
#endif

extern "C" void F(dgemm)(const char *ta,const char *tb,
		   const integer *n, const integer *k, const integer *l,
		   const double *alpha,const double *A,const integer *lda,
		   const double *B, const integer *ldb,
		   const double *beta,const double *C, const integer *ldc);
extern "C" void F(dpptrf)(const char *uplo,const integer *n, double* ap,
			  integer *info);
extern "C" void F(dspgst)(const integer *itype, const char *uplo,
			  const integer *n,
			  double* ap,const double *bp,integer *info);
extern "C" void F(dtptri)(const char *uplo,const char *diag,const integer *n,
		       double* ap,integer *info);
  /* unit triangular means that a value of 1.0 is assumed  */
  /* for the diagonal elements (hence diagonal not stored in packed format) */
extern "C" void F(dtrmm)(const char *side,const char *uplo,const char *transa,
		       const char *diag,const integer *m,const integer *n,
		       const double *alpha,const double *A,const integer *lda,
		       double *B,const integer *ldb);
extern "C" void F(dsygv)(const integer *itype,const char *jobz,
		       const char *uplo,const integer *n,
		       double *A,const integer *lda,double *B,
			 const integer *ldb,
			 double* w,double* work,const integer *lwork,
			 integer *info);
extern "C" void F(sgemm)(const char *ta,const char *tb,
		   const integer *n, const integer *k, const integer *l,
		   const float *alpha,const float *A,const integer *lda,
		   const float *B, const integer *ldb,
		   const float *beta,const float *C, const integer *ldc);
extern "C" void F(spptrf)(const char *uplo,const integer *n, float* ap,
			  integer *info);
extern "C" void F(sspgst)(const integer *itype, const char *uplo,
			  const integer *n,
			  float* ap,const float *bp,integer *info);
extern "C" void F(stptri)(const char *uplo,const char *diag,const integer *n,
		       float* ap,integer *info);
  /* unit triangular means that a value of 1.0 is assumed  */
  /* for the diagonal elements (hence diagonal not stored in packed format) */
extern "C" void F(strmm)(const char *side,const char *uplo,const char *transa,
		       const char *diag,const integer *m,const integer *n,
		       const float *alpha,const float *A,const integer *lda,
		       float *B,const integer *ldb);
extern "C" void F(ssygv)(const integer *itype,const char *jobz,
		       const char *uplo,const integer *n,
		       float *A,const integer *lda,float *B,const integer *ldb,
		       float* w,float* work,const integer *lwork,integer *info);

namespace bsm
{
  /*************** Default version throws exception */
  template<class T>
  inline void gemm(const char *ta,const char *tb,
		   const integer n, const integer k, const integer l,
		   const T *alpha,const T *A,const integer lda,
		   const T *B, const integer ldb,
		   const T *beta,T *C, const integer ldc)
  {
    throw Failure("genericblas does not support used type");
  }
  
  template<class T>
  inline void pptrf(const char *uplo,const integer *n, T* ap, integer *info)
  {
    throw Failure("genericblas does not support used type");
  }
  
  template<class T>
  inline void spgst(const integer *itype, const char *uplo,const integer *n,
		    T* ap,const T *bp,integer *info)
  {
    throw Failure("genericblas does not support used type");
  }
  
  template<class T>
  inline static void tptri(const char *uplo,const char *diag,const integer *n,
		    T* ap,integer *info)
  {
    throw Failure("genericblas does not support used type");
  }
  
  template<class T>
  inline void trmm(const char *side,const char *uplo,const char *transa,
		   const char *diag,const integer *m,const integer *n,
		   const T *alpha,const T *A,const integer *lda,
		   T *B,const integer *ldb)
  {
    throw Failure("genericblas does not support used type");
  }
  
  template<class T>
  inline void sygv(const integer *itype,const char *jobz,
		   const char *uplo,const integer *n,
		   T *A,const integer *lda,T *B,const integer *ldb,
		   T* w,T* work,const integer *lwork,integer *info)
  {
    throw Failure("genericblas does not support used type");
  }
  
  /*************** Double specialization */
  template<>
  inline void gemm<double>(const char *ta,const char *tb,
			   const integer n, const integer k, const integer l,
		   const double *alpha,const double *A, const integer lda,
		   const double *B, const integer ldb,
		   const double *beta,double *C, integer ldc)
  {
    F(dgemm)(ta,tb,&n,&k,&l,alpha,A,&lda,B,&ldb,beta,C,&ldc);
  }
  
  template<>
  inline void pptrf<double>(const char *uplo,const integer *n, double* ap, integer *info)
  {
    F(dpptrf)(uplo,n,ap,info);
  }
  
  template<>
  inline void spgst<double>(const integer *itype, const char *uplo,const integer *n,
		    double* ap,const double *bp,integer *info)
  {
    F(dspgst)(itype,uplo,n,ap,bp,info);
  }
  
  template<>
  inline void tptri<double>(const char *uplo,const char *diag,const integer *n,
		    double* ap,integer *info)
  {
    F(dtptri)(uplo,diag,n,ap,info);
  }
  
  template<>
  inline void trmm<double>(const char *side,const char *uplo,const char *transa,
		   const char *diag,const integer *m,const integer *n,
		   const double *alpha,const double *A,const integer *lda,
		   double *B,const integer *ldb)
  {
    F(dtrmm)(side,uplo,transa,diag,m,n,alpha,A,lda,B,ldb);
  }
  
  template<>
  inline void sygv<double>(const integer *itype,const char *jobz,
		   const char *uplo,const integer *n,
		   double *A,const integer *lda,double *B,const integer *ldb,
		   double* w,double* work,const integer *lwork,integer *info)
  {
    F(dsygv)(itype,jobz,uplo,n,A,lda,B,ldb,w,work,lwork,info);
  }
  
  #if 1
  /*************** Single specialization */
  template<>
  inline void gemm<float>(const char *ta,const char *tb,
			  const integer n, const integer k, const integer l,
			  const float *alpha,const float *A,const integer lda,
			  const float *B, const integer ldb,
			  const float *beta,float *C, const integer ldc)
  {
    F(sgemm)(ta,tb,&n,&k,&l,alpha,A,&lda,B,&ldb,beta,C,&ldc);
  }
template<>
  inline void pptrf<float>(const char *uplo,const integer *n, float* ap, integer *info)
  {
    F(spptrf)(uplo,n,ap,info);
  }
  
  template<>
  inline void spgst<float>(const integer *itype, const char *uplo,const integer *n,
		    float* ap,const float *bp,integer *info)
  {
    F(sspgst)(itype,uplo,n,ap,bp,info);
  }
  
  template<>
  inline void tptri<float>(const char *uplo,const char *diag,const integer *n,
		    float* ap,integer *info)
  {
    F(stptri)(uplo,diag,n,ap,info);
  }
  
  template<>
  inline void trmm<float>(const char *side,const char *uplo,const char *transa,
		   const char *diag,const integer *m,const integer *n,
		   const float *alpha,const float *A,const integer *lda,
		   float *B,const integer *ldb)
  {
    F(strmm)(side,uplo,transa,diag,m,n,alpha,A,lda,B,ldb);
  }
  
  template<>
  inline void sygv<float>(const integer *itype,const char *jobz,
		   const char *uplo,const integer *n,
		   float *A,const integer *lda,float *B,const integer *ldb,
		   float* w,float* work,const integer *lwork,integer *info)
  {
    F(ssygv)(itype,jobz,uplo,n,A,lda,B,ldb,w,work,lwork,info);
  }
  
  #endif
}
