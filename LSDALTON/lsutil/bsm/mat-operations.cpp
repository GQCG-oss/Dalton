
#include <cmath>
#include <iostream>
#include <exception>

#include "DirectDensity.h"
#include "DivideSpace.h"
#include "Failure.h"
#include "genericblas.h"

#include "mat-operations.h"

#if defined(HAVE_LIBESSL)
/* this is just a hack, but AIX is the only system with non-standard
 * symbol mangling */
#define F(x) x
#else
#define F(x) x ## _
#endif

#define EXTERNC extern "C"

EXTERNC void F(dpptrf)(const char *uplo,const integer *n, double* ap,
		       integer *info);
EXTERNC void F(dtptri)(const char *uplo,const char *diag,const integer *n,
                     double* ap,integer *info);
  /* unit triangular means that a value of 1.0 is assumed  */
  /* for the diagonal elements (hence diagonal not stored in packed format) */
EXTERNC void F(dtrmm)(const char *side,const char *uplo,const char *transa,
                    const char *diag,const integer *m,const integer *n,
                    const double *alpha,const double *A,const integer *lda,
                    double *B, const integer *ldb);

static bsm::Permutation* bsm_lib_divide = NULL;
static bsm::Permutation* bsm_lib_divide_save = NULL;
static real bsm_lib_epsilon        = 1e-6;
static real bsm_lib_tracethreshold = 0.1;
static integer  bsm_lib_ntight         = 6;

EXTERNC void
bsm_lib_init_(double* x, double* y, double* z, integer* atomstart, 
	      integer* natoms, integer* maxblocksize, real* epsilon,
	      real *tracethr, integer *ntight) {
  try {
#if 0
    for(integer i=0; i<*natoms; i++)
      printf("%4d: [%10.5f %10.5f %10.5f]: %5i\n", i, x[i], y[i], z[i],
             atomstart[i]);
#endif
    bsm_lib_divide = 
    new bsm::DivideSpace(x,y,z,atomstart,*natoms,*maxblocksize);

    bsm_lib_divide_save = bsm_lib_divide;
    //new bsm::Permutation(x,y,z,atomstart,*natoms,*maxblocksize);
    bsm_lib_epsilon = *epsilon;
    bsm_lib_tracethreshold = *tracethr;
    bsm_lib_ntight         = *ntight;
    bsm::Block::setcutoffvalue(*epsilon);
    printf("trace threshold: %g ntight: %d\n",
	   bsm_lib_tracethreshold, bsm_lib_ntight);
  }
  catch (bsm::Failure f) {
    std::cout<<f.what()<<std::endl;
    std::exit(1);
  } 
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }

}

EXTERNC void
bsm_lib_free_(long int *perm) { delete (bsm::Permutation*) *perm; *perm = 0; }


EXTERNC void
bsm_lib_setpermutation_(long int *perm) {
  bsm_lib_divide_save = (bsm::Permutation*) *perm; 
}

EXTERNC void
bsm_lib_getpermutation_(long int *perm) {
  *perm = (long int) bsm_lib_divide_save;
}

EXTERNC void
bsm_init_(DalMatrix *a, const integer *size, const integer *size2) {
  
  try {
    if(!bsm_lib_divide) throw bsm::Failure("bsmlib not initialized!");
    a->size = *size; a->size2 = *size2;
    a->permutation = (void *) bsm_lib_divide_save;
    bsm_lib_divide = bsm_lib_divide_save;
    if(*size == *size2) {
      bsm::BlockSparse* matrix = 
        new bsm::BlockSparse(*bsm_lib_divide, *size);
      a->private_data = matrix;
    } else {
      printf("BSM used to store full data %d x %d?\n",
	     integer(a->size), integer(a->size2));
      a->private_data = new real[a->size * a->size2];
    }
  }  
  catch (bsm::Failure f) {
    std::cout<<f.what()<<std::endl;
    std::exit(1);
  } 
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }
}

EXTERNC void
bsm_free_(DalMatrix *a) { 
  if(a->size == a->size2) {
    bsm_lib_divide = (bsm::Permutation*)a->permutation;
    bsm::BlockSparse *priva = (bsm::BlockSparse*)(a->private_data);
    delete priva; 
  } else {
    real *m = (real*)(a->private_data);
    delete m;
  }
}

EXTERNC void
bsm_init_from_full_(DalMatrix* a, const integer* s1, const integer *s2,
                    real* afull, real *sparsity) {
  try {
    bsm_lib_divide = (bsm::Permutation*)a->permutation;
    a->size = *s1; a->size2 = *s2;
    if( *s1 == *s2) {
      bsm::BlockSparse* matrix = 
        new bsm::BlockSparse(*bsm_lib_divide, *s1, afull, bsm_lib_epsilon);
      a->private_data = matrix;
      *sparsity = matrix->elsparsity();
    } else {
      real *v = new real[a->size * a->size2];
      a->private_data = v;
      memcpy(v, afull, a->size * a->size2 * sizeof(real));
      *sparsity = 1.0;
    }
  }  
  catch (bsm::Failure f) {
    std::cout<<f.what()<<std::endl;
    std::exit(1);
  } 
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }
}

EXTERNC void
bsm_to_full_(DalMatrix *a, real *afull, real *sparsity) {
  try {
    bsm_lib_divide = (bsm::Permutation*)a->permutation;
    if(a->size == a->size2) {
      bsm::BlockSparse *priva = (bsm::BlockSparse*)(a->private_data);
      priva->fullmatrix(afull);
      *sparsity = priva->elsparsity();
    } else {
      real *v = (real*)(a->private_data);
      memcpy(afull, v, a->size * a->size2 * sizeof(real));
      *sparsity = 1;
    }
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }
}

EXTERNC void
bsm_extract_diag_(DalMatrix *a, real *adiag) {
  try {
    bsm_lib_divide = (bsm::Permutation*)a->permutation;
    if(a->size != a->size2)
      throw bsm::DimensionF("bsm_extract_diag: matrix not square");
    bsm::BlockSparse *priva = (bsm::BlockSparse*)(a->private_data);
    priva->extract_diag(adiag);
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }
}


/* assignment operator a = b (a should be initialized) */
EXTERNC void
bsm_assign_(DalMatrix *a, const DalMatrix *b) {
  try {
    bsm_lib_divide = (bsm::Permutation*)a->permutation;
    if(b->size == b->size2) {
      bsm::BlockSparse *priva = (bsm::BlockSparse*)(a->private_data);
      bsm::BlockSparse *privb = (bsm::BlockSparse*)(b->private_data);
      *priva = *privb; 
      a->size = b->size;
      a->size2 = b->size2;
    } else {
      if(a->size != b->size || a->size2 != b->size2)
        throw bsm::DimensionF("bsm_assign: sizes do not match.");
      memcpy(a->private_data, b->private_data,
             a->size * a->size2 * sizeof(real));
    }
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }
}

EXTERNC void
bsm_copy_(real *alpha, const DalMatrix *a, DalMatrix *b) {
  try {
    bsm_lib_divide = (bsm::Permutation*)a->permutation;
    if(a->size == a->size2 && b->size == b->size2) {
      bsm::BlockSparse *priva = (bsm::BlockSparse*)(a->private_data);
      bsm::BlockSparse *privb = (bsm::BlockSparse*)(b->private_data);
      *privb = *priva; 
      *privb *= *alpha;
      b->size = a->size;
      b->size2 = a->size2;
    } else {
      if(a->size != b->size || a->size2 != b->size2)
        throw bsm::DimensionF("bsm_copy: sizes do not match.");
      real *av = (real*)(a->private_data);
      real *bv = (real*)(b->private_data);
      for(integer i=a->size * a->size2 -1; i>=0; i--)
        bv[i] = *alpha * av[i];
    }
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }
}

EXTERNC real
bsm_tr_(DalMatrix *a) {
  try {
    bsm_lib_divide = (bsm::Permutation*)a->permutation;
    if(a->size != a->size2)
      throw bsm::DimensionF("bsm_tr: matrix not square");
    bsm::BlockSparse *priva = (bsm::BlockSparse*)(a->private_data);
    return priva->trace();
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }
}

EXTERNC real
bsm_tratransb_(const DalMatrix *a, const DalMatrix *b) {
  try {
    bsm_lib_divide = (bsm::Permutation*)a->permutation;
    if(a->size == a->size2 &&
       b->size == b->size2) {
      bsm::BlockSparse *priva = (bsm::BlockSparse*)(a->private_data);
      bsm::BlockSparse *privb = (bsm::BlockSparse*)(b->private_data);
      return priva->trace_atransb(*privb);
    } else {
      if(a->size != b->size ||
         a->size2 != 1 || b->size2 != 1) {
        printf("A: %d x %d B: %d x %d\n", integer(a->size), integer(a->size2),
               integer(b->size), integer(b->size2));
        throw bsm::DimensionF("bsm_trtransb: that's sick!");
      }
      real res = 0;
      real *av = (real*)(a->private_data);
      real *bv = (real*)(b->private_data);
      for(integer i = a->size -1; i>=0; i--)
        res += av[i]*bv[i];
      return res;
    }
  }
  catch(bsm::Failure f) {
    std::cout<<f.what()<<std::endl;
    std::exit(1);
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }
} 
EXTERNC real
bsm_trab_(const DalMatrix *a, const DalMatrix *b) {
  try {
    bsm_lib_divide = (bsm::Permutation*)a->permutation;
    if(a->size != a->size2)
      throw bsm::DimensionF("bsm_trab: matrix A not square");
    if(b->size != b->size2)
      throw bsm::DimensionF("bsm_trab: matrix B not square");
    bsm::BlockSparse *priva = (bsm::BlockSparse*)(a->private_data);
    bsm::BlockSparse *privb = (bsm::BlockSparse*)(b->private_data);
    return priva->trace_ab(*privb);
  }
  catch(bsm::Failure f) {
    std::cout<<f.what()<<std::endl;
    std::exit(1);
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }
} 


EXTERNC void
bsm_add_(real *alpha, const DalMatrix *a, real *beta, const DalMatrix *b,
         DalMatrix* c) {
  try {
    bsm_lib_divide = (bsm::Permutation*)a->permutation;
    if(a->size == a->size2) {
      if(b->size != b->size2 || c->size != c->size2)
        throw bsm::DimensionF("bsm_add abused");
      bsm::BlockSparse *priva = (bsm::BlockSparse*)(a->private_data);
      bsm::BlockSparse *privb = (bsm::BlockSparse*)(b->private_data);
      bsm::BlockSparse *privc = (bsm::BlockSparse*)(c->private_data);
      bsm::BlockSparse* matrix = new bsm::BlockSparse(*priva,
                                                      *privb,
                                                      *alpha, *beta);
      delete privc;
      c->private_data = matrix;
      c->size = a->size;
      c->size2 = a->size2;
    } else {
      if(b->size == b->size2 || c->size == c->size2)
        throw bsm::DimensionF("bsm_add terribly abused");
      if(a->size != b->size || a->size2 != b->size2 ||
         a->size != c->size || a->size2 != c->size2)
        throw bsm::DimensionF("bsm_add 'orribly abused");
      real *av = (real*) (a->private_data);
      real *bv = (real*) (b->private_data);
      real *cv = (real*) (c->private_data);
      for(integer i= a->size * a->size2 -1; i>= 0; i--)
        cv[i] = *alpha *av[i] + *beta * bv[i];
    }
  }
  catch(bsm::Failure f) {
    std::cout<<f.what()<<std::endl;
    std::exit(1);
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }
 
}

EXTERNC void
bsm_daxpy_(real *alpha, const DalMatrix *a, DalMatrix* c) {
  try {
    bsm_lib_divide = (bsm::Permutation*)a->permutation;
    if(a->size == a->size2) {
      if(c->size != c->size2)
        throw bsm::DimensionF("bsm_daxpy abused");
      bsm::BlockSparse *priva = (bsm::BlockSparse*)(a->private_data);
      bsm::BlockSparse *privc = (bsm::BlockSparse*)(c->private_data);
      bsm::BlockSparse* matrix = new bsm::BlockSparse(*priva,
                                                      *privc,
                                                      *alpha, 1.0);
      delete privc;
      c->private_data = matrix;
    } else { /* VECTOR VERSION */
      if(c->size == c->size2)
        throw bsm::DimensionF("bsm_daxpy shockingly abused");
      if(a->size != c->size || a->size2 != c->size2)
        throw bsm::DimensionF("bsm_daxpy violently abused");
      real *av = (real*) (a->private_data);
      real *cv = (real*) (c->private_data);

      for(integer i=a->size * a->size2-1; i>=0; i--)
        cv[i] += *alpha * av[i];
    }
  }
  catch(bsm::Failure f) {
    std::cout<<f.what()<<std::endl;
    std::exit(1);
  } 
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }
}

EXTERNC void
bsm_precond_(const DalMatrix *d, const DalMatrix *v,  DalMatrix *vp)
{
  // xprec(i) := x(i)/M(i,i)         
  try {
    bsm_lib_divide = (bsm::Permutation*)d->permutation;
    if(d->size != d->size2)
      throw bsm::DimensionF("bsm_precond: precond matrix not square");
    if(v->size != d->size || v->size2 != 1)
      throw bsm::DimensionF("bsm_precond: first vector is not a vector");
    if(vp->size != d->size || vp->size2 != 1)
      throw bsm::DimensionF("bsm_precond: second vector is not a vector");
    bsm::BlockSparse *privd  = (bsm::BlockSparse*)(d->private_data);
    real *privv  = (real*)(v->private_data);
    real *privvp = (real*)(vp->private_data);
#if 0
    real *diag = new real[d->size];
    privd->extract_diag(diag);
    for(integer i=0; i<d->size; i++) {
      if(fabs(diag[i])>1e-9)
        privvp[i] = privv[i]/diag[i];
    }
    delete []diag;
#endif
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }
}

EXTERNC void
bsm_ao_precond_(const integer *symm, const real * omega,
                DalMatrix *FUP, DalMatrix *FUQ, 
                DalMatrix *DU,  DalMatrix *X_AO)
{
  try {
    bsm_lib_divide = (bsm::Permutation*)FUP->permutation;
    if(FUP->size != FUP->size2 ||
       FUQ->size != FUQ->size2 ||
       DU->size  != DU->size2  ||
       X_AO->size != X_AO->size2)
      throw bsm::DimensionF("bsm_ao_precond severely abused");
      
    integer n = ((bsm::BlockSparse*)(FUP->private_data))->gettotalcols();
    real *fup = new real[n];
    real *fuq = new real[n];
    real *du  = new real[n];
    ((bsm::BlockSparse*)(FUP->private_data))->extract_diag_internal(fup);
    ((bsm::BlockSparse*)(FUQ->private_data))->extract_diag_internal(fuq);
    ((bsm::BlockSparse*)(DU->private_data))->extract_diag_internal(du);
    ((bsm::BlockSparse*)(X_AO->private_data))->precond_ao(*symm, fup, fuq,
                                                          du, *omega);
    delete []fup;
    delete []fuq;
    delete []du;
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }
}

EXTERNC void
bsm_scal_(real *alpha, DalMatrix *a) {
  try {
    bsm_lib_divide = (bsm::Permutation*)a->permutation;
    if(a->size == a->size2) {
      bsm::BlockSparse *priva = (bsm::BlockSparse*)(a->private_data);
      *priva *= *alpha;
    } else {
      real *av = (real*) (a->private_data);
      for(integer i = a->size *a->size2 -1; i>=0; i--)
        av[i] *= *alpha;
    }
  }
  catch(bsm::Failure f) {
    std::cout<<f.what()<<std::endl;
    std::exit(1);
  } 
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }
}

EXTERNC void
bsm_scal_dia_(real *alpha, DalMatrix *a) {
  try {
    bsm_lib_divide = (bsm::Permutation*)a->permutation;
    if(a->size != a->size2)
      throw bsm::DimensionF("bsm_scal_dia: matrix is not square");
    bsm::BlockSparse *priva = (bsm::BlockSparse*)(a->private_data);
    priva->scal_dia(*alpha);
  }
  catch(bsm::Failure f) {
    std::cout<<f.what()<<std::endl;
    std::exit(1);
  } 
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }
}

EXTERNC void
bsm_zerohalf_(const char *triangle, DalMatrix *A) {
  try {
    bsm_lib_divide = (bsm::Permutation*)A->permutation;
    if(A->size != A->size2)
      throw bsm::DimensionF("bsm_zerohalf: matrix is not square");
    bsm::BlockSparse *priva = (bsm::BlockSparse*)(A->private_data);
    priva->zerohalf(triangle);
  }
  catch(bsm::Failure f) {
    std::cout<<f.what()<<std::endl;
    std::exit(1);
  } 
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }
}

static real*
getVec(DalMatrix *a, const char *s) {
  real *af;
  bsm_lib_divide = (bsm::Permutation*)a->permutation;
  if(a->size == a->size2) {
    bsm::BlockSparse *priva = (bsm::BlockSparse*)(a->private_data);
    af = new real[a->size * a->size2];
    priva->fullmatrix(af);
  } else af = (real*) (a->private_data);
  return af;
}

EXTERNC void
bsm_mul_(DalMatrix *a, DalMatrix *b, 
         char *transa, char *transb,
         real *alpha, real *beta, DalMatrix *c) {
  try {
  bsm_lib_divide = (bsm::Permutation*)a->permutation;
    if(a->size == a->size2 &&
       b->size == b->size2 &&
       c->size == c->size2) {
      bsm::BlockSparse *priva = (bsm::BlockSparse*)(a->private_data);
      bsm::BlockSparse *privb = (bsm::BlockSparse*)(b->private_data);
      bsm::BlockSparse *privc = (bsm::BlockSparse*)(c->private_data);
      bsm::BlockSparse *result, *privat = NULL, *privbt = NULL;
      if( *transa == 'T' || *transa == 't') {
        privat = new bsm::BlockSparse(*priva, true);
      } else privat = priva;
      if( *transb == 'T' || *transb == 't') {
        privbt = new bsm::BlockSparse(*privb, true);
      } else privbt = privb;

      result = new bsm::BlockSparse(*bsm_lib_divide, a->size);
      bsm::BlockSparse::multiply(*privat, *privbt, *result, *alpha, *beta);
      if( *transa == 'T' || *transa == 't') delete privat;
      if( *transb == 'T' || *transb == 't') delete privbt;
      delete privc;
      c->private_data = result;
    } else {
      real *af, *bf, *cf;
      printf("Multiplication of sparse matrices by dense vectors "
             "requires conversions %c %c...\n", *transa, *transb);
      af = getVec(a, "A");
      bf = getVec(b, "B");
      cf = getVec(c, "C");
      integer k = (*transa == 'n' || *transa == 'N') ? b->size : b->size2;
      F(dgemm)(transa, transb, &c->size, &c->size2, &k, alpha, af, &a->size,
               bf, &b->size, beta, cf, &c->size);
      if(a->size == a->size2) delete af;
      if(b->size == b->size2) delete bf;
      if(c->size == c->size2) {
        printf("Multiple conversions...\n");
        bsm::BlockSparse *pc = (bsm::BlockSparse*)(c->private_data);
        delete pc;
        pc = new bsm::BlockSparse(*bsm_lib_divide, c->size, cf,
                                  bsm_lib_epsilon);
        c->private_data = pc;
        delete cf;
      }
    }
  }
  catch(bsm::Failure f) {
    std::cout<<f.what()<<std::endl;
    std::exit(1);
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }
 
}

EXTERNC real
bsm_frob_(DalMatrix *a) {
  try {
  bsm_lib_divide = (bsm::Permutation*)a->permutation;
    if(a->size == a->size2) {
      bsm::BlockSparse *priva = (bsm::BlockSparse*)(a->private_data);
      return priva->frob();
    } else {
      real *av = (real*)(a->private_data);
      real r = 0;
      for(integer i=a->size *a->size2 - 1; i>=0; i--)
        r += av[i]*av[i];
      return sqrt(r);
    }
  }
  catch(bsm::Failure f) {
    std::cout<<f.what()<<std::endl;
    std::exit(1);
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }

}    

EXTERNC real
bsm_max_(DalMatrix *a) {
  try {
  bsm_lib_divide = (bsm::Permutation*)a->permutation;
    if(a->size != a->size2)
      throw bsm::DimensionF("bsm_max: matrix A is not square");
    bsm::BlockSparse *priva = (bsm::BlockSparse*)(a->private_data);
    return priva->max();
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }

}    

EXTERNC void
bsm_max_diag_(DalMatrix *a, integer *pos, real *val) {
  try {
  bsm_lib_divide = (bsm::Permutation*)a->permutation;
    if(a->size != a->size2)
      throw bsm::DimensionF("bsm_max_diag: matrix A is not square");
    bsm::BlockSparse *priva = (bsm::BlockSparse*)(a->private_data);
    integer pproxy;
    priva->max_abs_diag(&pproxy, val);
    *pos = pproxy;
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }

}    
EXTERNC real
bsm_outdia_sqnorm2_(DalMatrix *a)
{
  bsm_lib_divide = (bsm::Permutation*)a->permutation;
  if(a->size != a->size2)
    throw bsm::DimensionF("bsm_outdia_sqnorm2: matrix A is not square");
  bsm::BlockSparse *privF = (bsm::BlockSparse*)(a->private_data);
  return privF->sum_outdia_sqnorm2();
}

EXTERNC void
bsm_print_struct_(DalMatrix *a)
{
  bsm_lib_divide = (bsm::Permutation*)a->permutation;
  if(a->size != a->size2)
    throw bsm::DimensionF("bsm_print_struct: matrix A is not square");
  bsm::BlockSparse *privF = (bsm::BlockSparse*)(a->private_data);
  return privF->print_struct();
}

EXTERNC void
bsm_transpose_(DalMatrix *a, DalMatrix *b)
{
  try {
  bsm_lib_divide = (bsm::Permutation*)a->permutation;
    if(a->size != a->size2)
      throw bsm::DimensionF("bsm_transpose: matrix A is not square");
    if(b->size != b->size2)
      throw bsm::DimensionF("bsm_transpose: matrix A is not square");
    bsm::BlockSparse *priva = (bsm::BlockSparse*)(a->private_data);
    bsm::BlockSparse *privb = (bsm::BlockSparse*)(b->private_data);
    if(privb)
      delete privb;
    b->private_data = new bsm::BlockSparse(*priva, true);
  } catch(bsm::Failure f) {
    std::cout<<f.what()<<std::endl;
    std::exit(1);
  }
}


#if 0
EXTERNC void
bsm_d_from_f(DalMatrix *F, DalMatrix *Dnew, 
             real* lmin, real* lmax, integer* nshells) {
  try {
    bsm::BlockSparse *privF = (bsm::BlockSparse*)(F->private_data);
    bsm::BlockSparse *privDnew = (bsm::BlockSparse*)(Dnew->private_data);


    bsm::BlockSparse Fock(*privF);
    bsm::DirectDensity DD(Fock, *privDnew, *lmin, *lmax, *nshells,
			  bsm_lib_tracethreshold, bsm_lib_ntight);
    DD.finddensity();
  }
  catch(bsm::Failure f) {
    std::cout<<f.what()<<std::endl;
    std::exit(1);
  }
}
#endif

EXTERNC void
bsm_identity_(DalMatrix* I) {
  try {
  bsm_lib_divide = (bsm::Permutation*)I->permutation;
    if(I->size != I->size2)
      throw bsm::DimensionF("bsm_identity: matrix I is not square");
    bsm::BlockSparse *privI = (bsm::BlockSparse*)(I->private_data);  
    privI->identity();
  }
  catch(bsm::Failure f) {
    std::cout<<f.what()<<std::endl;
    std::exit(1);
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }
}

EXTERNC void
bsm_add_identity_(DalMatrix* A, real* alpha) {
  try {
  bsm_lib_divide = (bsm::Permutation*)A->permutation;
    if(A->size != A->size2)
      throw bsm::DimensionF("bsm_add_identity: matrix A is not square");
    bsm::BlockSparse *privA = (bsm::BlockSparse*)(A->private_data);  
    privA->addscaleye(*alpha);
  }
  catch(bsm::Failure f) {
    std::cout<<f.what()<<std::endl;
    std::exit(1);
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }

}

#if 1
EXTERNC void
bsm_d_from_f_(DalMatrix *F, DalMatrix* S, integer* nshells,
              DalMatrix *Dnew, integer *cycles) {
  try {
  bsm_lib_divide = (bsm::Permutation*)F->permutation;
    if(F->size != F->size2 ||
       S->size != S->size  ||
       Dnew->size != Dnew->size2)
      throw bsm::DimensionF("bsm_d_from_f: abused");
    bsm::BlockSparse *privF = (bsm::BlockSparse*)(F->private_data);
    bsm::BlockSparse *privDnew = (bsm::BlockSparse*)(Dnew->private_data);
    bsm::BlockSparse *privS = (bsm::BlockSparse*)(S->private_data);
    integer size = F->size;
    double* fockfull=new double [size*size];
    privF->fullmatrix(fockfull);

    
    // fockpacked never used?
    // double* fockpacked=new double [(size+1)*size/2];
    // bsm::fulltopacked(fockfull,fockpacked,size);

    double* overfull=new double [size*size];
    privS->fullmatrix(overfull);
    double* overpacked=new double [(size+1)*size/2];
    bsm::fulltopacked(overfull,overpacked,size);
    /**************  Reduction of eigenvalue problem to standard form */
    /**************  Step 1 S=U'U*/
    integer info;
    static real ONE = 1.0, TWO = 2.0;
    real tracediff;

    F(dpptrf)("U",&size,overpacked,&info); 
    if (info)
      {std::cout<<"Error: dpptrf_ info="<<info<<std::endl;
      std::exit(1);}
     
    /**************  Step 2 Compute inverse of U */
    F(dtptri)("U","N",&size,overpacked,&info);
    if (info)
      {std::cout<<"Error: dtptri_ info="<<info<<std::endl;
      std::exit(1);}  
  
    bsm::tripackedtofull(overpacked,overfull,size);
    delete[] overpacked;
     
    /************* Step 3 Compute 'new' Fock matrix = inv(U')*F*inv(U) */
    F(dtrmm)("L","U","T","N",&size,&size,&ONE,overfull,&size,fockfull,&size);
    F(dtrmm)("R","U","N","N",&size,&size,&ONE,overfull,&size,fockfull,&size);


    /**************** Find upper and lower bounds for eigenvalues */
    real lmin, lmax;
    bsm::gerschgorin(fockfull,size,lmin,lmax);
    /* std::cout<<"Eigenvalues between "<<lmin<<" and "<<lmax<<std::endl;*/
  
    /***************  Initialization of blocksparse matrices */   

    bsm::BlockSparse Fock(*bsm_lib_divide,size,fockfull, bsm_lib_epsilon);
    bsm::BlockSparse Dens(*bsm_lib_divide,size,fockfull, bsm_lib_epsilon);
    delete[] fockfull;

    //printf("eigenvalue estimated to be in range [%g %g]\n", lmin, lmax);
    /**************  Density computation*/
    bsm::DirectDensity DD(Fock, Dens, lmin, lmax, *nshells,
			  bsm_lib_tracethreshold, bsm_lib_ntight);
    integer cyclesProxy;
    DD.finddensity(&cyclesProxy, &tracediff);
    *cycles =cyclesProxy;
    if(tracediff>bsm_lib_tracethreshold)
      printf("purification converged in %d cycles, tracediff: %g thr=%g\n",
             integer(*cycles), tracediff, bsm_lib_tracethreshold);
     
    /**************  Backsubstitution to original basis of Fock matrix */  
    double* densappfull=new double [size*size];
    Dens.fullmatrix(densappfull);
    F(dtrmm)("L","U","N","N",&size,&size,&ONE,overfull,&size,densappfull,&size);
    F(dtrmm)("R","U","T","N",&size,&size,&ONE,overfull,&size,densappfull,&size);
    delete[] overfull;
    delete privDnew;

    bsm::BlockSparse* result =
      new bsm::BlockSparse(*bsm_lib_divide, size, densappfull,
                           bsm_lib_epsilon);
    //delete privc; ????
    delete[] densappfull;
    Dnew->private_data = result;
  }
  catch(bsm::Failure f) {
    std::cout<<f.what()<<std::endl;
    std::exit(1);
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }

}
#endif

/** Returns the ratio of stored elements to the total number of
    elements in the matrix. */
EXTERNC void
bsm_get_sparsity_(DalMatrix* A, real *sparsity) {
  bsm_lib_divide = (bsm::Permutation*)A->permutation;
  bsm::BlockSparse *priva = (bsm::BlockSparse*)(A->private_data);
  *sparsity = priva->elsparsity();
}

EXTERNC void
bsm_write_to_unit_(integer *iunit, DalMatrix* A,
                   void (*write_int) (void *unit, integer *cnt, integer *where),
                   void (*write_real)(void *unit, integer *cnt, real *where)) {
  try {
  bsm_lib_divide = (bsm::Permutation*)A->permutation;
    if(A->size == A->size2) {
      bsm::BlockSparse *privA = (bsm::BlockSparse*)(A->private_data);
      privA->write(iunit, write_int, write_real);
    } else {
      integer sz = A->size * A->size2;
      write_real(iunit, &sz, (real*)(A->private_data));
    }
  } catch(bsm::Failure f) {
    std::cout<<f.what()<<std::endl;
    std::exit(1);
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }
}

EXTERNC void
bsm_read_from_unit_(integer *iunit, DalMatrix* A,
                    void (*read_int) (void *data, integer *cnt, integer *where),
                    void (*read_real)(void *data, integer *cnt, real *where)) {

  try {
  bsm_lib_divide = (bsm::Permutation*)A->permutation;
    if(A->size == A->size2) {
      bsm::BlockSparse *privA = (bsm::BlockSparse*)(A->private_data);
      privA->read(iunit, read_int, read_real);
    } else {
      integer sz = A->size * A->size2;
      read_real(iunit, &sz, (real*) (A->private_data));
    }
  } catch(bsm::Failure f) {
    std::cout<<f.what()<<std::endl;
    std::exit(1);
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what();
    std::exit(1);
  }
}
