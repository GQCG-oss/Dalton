/* This file contains, together with "test3.cpp", the source code   *
 * to test the blocked sparse matrix library qcsparsemat            *
 * The test programs should be run without arguments and will       *
 * print an "OK" to standard output if the tests were successful    *
 * test3.cpp: Tests basic matrix manipulations like multiplication  *
 * testDD.cpp: Tests the density matrix purification scheme         *
 *                                                                  *
 *                                                                  *
 *     \\\|||///  \    Emanuel Rubensson, 2005                      *
 *     \ ~   ~ /   \                                                *
 *     | @   @ |    \  mail:  emanuel@theochem.kth.se               *
 * oOo---(_)---oOo---\----------------------------------------------*
 *                                                                  */

#include "config.h"

#include <fstream>
#include <iomanip> /* For setprecision in fstream */
#include <iostream>
#include <ctime>
#include <stdio.h>

#include "DivideSpace.h"
#include "DoubleSparse.h"
#include "DirectDensity.h"

using namespace bsm;

#if defined(HAVE_LIBESSL)
/* this is just a hack, but AIX is the only system with non-standard
 * symbol mangling */
#define F(x) x
#else
#define F(x) x ## _
#endif

extern "C" void F(dgemm)(const char *ta,const char *tb,
			 const integer *n, const integer *k, const integer *l,
			 const double *alpha,const double *A,const integer *lda,
			 const double *B, const integer *ldb,
			 const double *beta,const double *C, const integer *ldc);
extern "C" void F(dpptrf)(const char *uplo,const integer *n, double* ap, integer *info);
extern "C" void F(dspgst)(const integer *itype, const char *uplo,const integer *n,
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
		       double *A,const integer *lda,double *B,const integer *ldb,
		       double* w,double* work,const integer *lwork,integer *info);

static double maxdiff(const real* f1,const real* f2,integer size);

const integer one=1;
const double zero=0;
const double two=2;
integer main(integer argc,char* argv[])
{
  using namespace bsm;
  
  try
   {
     integer nshells=155;
     double lmin;
     double lmax;
     clock_t start, end,totalstart,totalend;
     integer info;
     float tm,tmnaive;
     integer size=237;
     integer natoms=73;
     double epsilon=pow(10,-7);
     Block::setcutoffvalue(epsilon);
     
     char atomsPath[200] = "atoms.txt";
     char fockPath[200] = "matrix_F";
     char densPath[200] = "matrix_D";
     char overPath[200] = "matrix_S";
     
     integer m=size;
     integer n=size;
     integer k=size;
     double* xpos=new double[natoms];  
     double* ypos=new double[natoms];
     double* zpos=new double[natoms];
     integer* atomstart=new integer[natoms+1];
     
     
     char* atomfile(atomsPath);
     std::ifstream input(atomfile);
     if (!input)
       {
	 std::cout << "Cannot open inputfile " << atomfile << std::endl;
	 std::exit(1);
       }
     input >> std::setprecision(10);
     for(integer i = 0 ; i < natoms ; i++)
       {
	 input >> xpos[i];
	 input >> ypos[i];
	 input >> zpos[i];
	 input >> atomstart[i];
       }
     atomstart[natoms]=size;
     
     FILE* fockfile=fopen(fockPath,"rb");
     FILE* densfile=fopen(densPath,"rb");
     FILE* overfile=fopen(overPath,"rb");
     
     if (!fockfile || !densfile || !overfile)
       {
	 std::cout<<"Cannot open inputfile "<<std::endl;
	 std::exit(1);
       }
     double* fockfull=new double [size*size];
     fread(fockfull,sizeof(double),size*size,fockfile);
     double* fockpacked=new double [(size+1)*size/2];
     fulltopacked(fockfull,fockpacked,size);
     
     double* overfull=new double [size*size];
     fread(overfull,sizeof(double),size*size,overfile);
     double* overpacked=new double [(size+1)*size/2];
     fulltopacked(overfull,overpacked,size);
     
     
     /************ Get info of how large working array that is needed  */
     double* w=new double[size];
     integer lwork=-1;
     double* work=new double;
     F(dsygv)(&one,"V","U",&size,fockfull,&size,overfull,&size,w,work,&lwork,&info);
     if (info)
       {std::cout<<"Error: dsygv_1 info="<<info<<std::endl;
       std::exit(1);}
     lwork=(integer)work[0];
     delete work;
     work=new double[lwork];
     /************** Diagonalization  */
     F(dsygv)(&one,"V","U",&size,fockfull,&size,overfull,&size,w,work,&lwork,&info);
     if (info)
       {std::cout<<"Error: dsygv_2 info="<<info<<std::endl;
       std::exit(1);}
     double* densdiag=new double [size*size];
     
     /************* Multiplication of matrices with eigenvectors : */
     F(dgemm)("N","T",&size,&size,&nshells,&two,fockfull,&size,fockfull,&size,
	    &zero,densdiag,&size);
     
     
     packedtofull(fockpacked,fockfull,size); /* Overwrites diag result*/
     delete[] overfull;
     
     /**************  Reduction of eigenvalue problem to standard form */
     /**************  Step 1 S=U'U*/
     F(dpptrf)("U",&size,overpacked,&info); 
     if (info)
       {std::cout<<"Error: dpptrf_ info="<<info<<std::endl;
       std::exit(1);}
     
     /**************  Step 2 Compute inverse of U */
     F(dtptri)("U","N",&size,overpacked,&info);
     if (info)
       {std::cout<<"Error: dtptri_ info="<<info<<std::endl;
       std::exit(1);}  
     overfull=new double [size*size];
     tripackedtofull(overpacked,overfull,size);
     delete[] overpacked;
     
     /************* Step 3 Compute 'new' Fock matrix = inv(U')*F*inv(U) */
     F(dtrmm)("L","U","T","N",&size,&size,&ONE,overfull,&size,fockfull,&size);
     F(dtrmm)("R","U","N","N",&size,&size,&ONE,overfull,&size,fockfull,&size);
     
     
     /**************** Find upper and lower bounds for eigenvalues */
     gerschgorin(fockfull,size,lmin,lmax);
     /* std::cout<<"Eigenvalues between "<<lmin<<" and "<<lmax<<std::endl;*/
     
     /***************  Initialization of blocksparse matrices */
     
     DivideSpace divide(xpos,ypos,zpos,atomstart,natoms,67);
     BlockSparse Fock(divide,size,fockfull,epsilon);
     BlockSparse Dens(divide,size,fockfull,epsilon);
     delete[] fockfull;
     
     /**************  Density computation*/
     DirectDensity DD(Fock,Dens,lmin,lmax,nshells);
     DD.finddensity();
     
     /**************  Backsubstitution to original basis of Fock matrix */  
     double* densappfull=new double [size*size];
     Dens.fullmatrix(densappfull);
     F(dtrmm)("L","U","N","N",&size,&size,&two,overfull,&size,densappfull,&size);
     F(dtrmm)("R","U","T","N",&size,&size,&ONE,overfull,&size,densappfull,&size);
     delete[] overfull;
     double* densfull=new double [size*size];
     fread(densfull,sizeof(double),size*size,densfile);
     bool ok=true;  
     if (maxdiff(densappfull,densfull,size) > epsilon*100)
       {
	 ok=false;
	 std::cout<<"The error in the trace purification scheme is larger than expected"<<std::endl;
       }
     if (maxdiff(densdiag,densfull,size) > epsilon*100)
       {
	 ok=false;
	 std::cout<<"The error in the diagonalization is larger than expected"<<std::endl;
       }
     if (ok)
       {
	 std::cout<<"OK"<<std::endl;
       }
     
     delete[] densappfull;
     delete[] densfull;
     
   }
 catch(Failure f)
   {
     std::cout<<f.what()<<std::endl;
     std::exit(1);
   }
}

static double maxdiff(const real* f1,const real* f2,integer size)
{
  double diff=0;
  double tmpdiff;
  for(integer i=0;i<size*size;i++)
    {
      tmpdiff=fabs(f1[i]-f2[i]);
      if (tmpdiff>0)
	{
	  diff=(diff>tmpdiff ? diff : tmpdiff);
	}
    }
  return diff;
}
