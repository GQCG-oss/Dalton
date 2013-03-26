/* This file contains, together with "testDD.cpp", the source code  *
 * to test the blocked sparse matrix library qcsparsemat            *
 * The test programs should be run without arguments and will       *
 * print an "OK" to standard output if the tests were successful    *
 * test3.cpp: Tests basic matrix manipulations like multiplication  *
 *            See more information below                            *
 * testDD.cpp: Tests the density matrix purification scheme         *
 *                                                                  *
 *                                                                  *
 *     \\\|||///  \    Emanuel Rubensson, 2005                      *
 *     \ ~   ~ /   \                                                *
 *     | @   @ |    \  mail:  emanuel@theochem.kth.se               *
 * oOo---(_)---oOo---\----------------------------------------------*
 *                                                                  */
#include "DivideSpace.h"
#include "DoubleSparse.h"
#include "BlockSparse.h"
#include <fstream>
#include <iomanip> /* For setprecision in fstream */
#include <iostream>
#include <ctime>
#include <stdio.h>
using namespace bsm;
extern "C" void dgemm_(const char *ta,const char *tb,
		   const integer *n, const integer *k, const integer *l,
		   const double *alpha,const double *A,const integer *lda,
		   const double *B, const integer *ldb,
		   const double *beta,const double *C, const integer *ldc);


static double maxdiff(const real* f1,const real* f2,integer size);
integer main(integer argc,char* argv[])
{
using namespace bsm;
 try
   {
  clock_t start, end;
  float tm,tmnaive;
  integer size = 237;//arg 3
  integer natoms = 73;//arg 2
  double epsilon = pow(10,-7);
  Block::setcutoffvalue(epsilon);
  integer tmpexp;
  

  char atomsPath[200] = "atoms.txt";
  char fockPath[200] = "matrix_F";
  char densPath[200] = "matrix_D";
  char overPath[200] = "matrix_S";

  integer m = size;
  integer n = size;
  integer k = size;
  double* xpos = new double[natoms];  
  double* ypos = new double[natoms];
  double* zpos = new double[natoms];
  integer* atomstart = new integer[natoms+1];
  
  
  char* atomfile(atomsPath);
  std::ifstream input(atomfile);
  if (!input)
    {
      std::cout<<"Cannot open inputfile "<<atomfile<<std::endl;
      std::exit(1);
    }
  input>>std::setprecision(10);
  for(integer i = 0 ; i < natoms ; i++)
    {
      input >> xpos[i];
      input >> ypos[i];
      input >> zpos[i];
      input >> atomstart[i];    
    }
  atomstart[natoms] = size;
  
  double fockcoradius;
  double denscoradius;
  
  integer s;
 
  FILE* fockfile=fopen(fockPath,"rb");
  FILE* densfile=fopen(densPath,"rb");
  FILE* overfile=fopen(overPath,"rb");
    
  if (!fockfile || !densfile || !overfile)
    {
      std::cout<<"Cannot open inputfile "<<std::endl;
      std::exit(1);
    }
  double* fockfull=new double [size*size];
  double* densfull=new double [size*size];
  double* overfull=new double [size*size];
  double* prodfull=new double [size*size];
  double* sumfull=new double [size*size];
  double* tmpfull=new double [size*size];
  fread(fockfull,sizeof(double),size*size,fockfile);
  fread(densfull,sizeof(double),size*size,densfile);
  fread(overfull,sizeof(double),size*size,overfile);
  double prodmaxerror;
  double fockmaxerror;
  double densmaxerror;
  double summaxerror;
  double diff;
  /* Multiplication */
  dgemm_("N","N",&m,&n,&k,&ONE,fockfull,&m,densfull,&k,&ZERO,prodfull,&m);
  integer maxbf=67;
  DivideSpace divide(xpos,ypos,zpos,atomstart,natoms,maxbf);
  BlockSparse fock(divide,size,fockfull,epsilon);//fockcoradius);
  BlockSparse dens(divide,size,densfull,epsilon);//denscoradius);
  BlockSparse prod(divide,size);
  Block::zeroseconds();
  BlockSparse::multiply(fock,dens,prod);
  
  /* Addition */
  double alpha = 0.432;
  double beta = 72;
  for (integer i = 0; i < size * size; i++) {
    sumfull[i] = alpha * fockfull[i] + beta * densfull[i]; 
  }
  BlockSparse sum(fock, dens, alpha, beta);
  
  /* Trace */
  double tr=0;
  for (integer i = 0; i < size*size; i += size + 1) {
    tr += fockfull[i];
  }
  double tracediff = fabs(tr - fock.trace());

  /* Frobenius norm */
  double fr = 0;
  for (integer i = 0; i < size*size; i++) {
    fr += fockfull[i] * fockfull[i];
  }
  fr = sqrt(fr);
  double frobdiff = fabs(fr - fock.frob());
  
  /* Trace(A * B) */
  double tr_ab = 0;
  for (integer i = 0; i < size; i++) {
    for (integer j = 0; j < size; j++) {
      tr_ab += fockfull[i + j * size] * densfull[i * size + j];
    }
  }
  double tr_abdiff = fabs(tr_ab - fock.trace_ab(dens));

  /* Trace(A' * B) */
  double tr_atransb = 0;
  for (integer i = 0; i < size; i++) {
    for (integer j = 0; j < size; j++) {
      tr_atransb += fockfull[j + i * size] * densfull[i * size + j];
    }
  }
  double tr_atransbdiff = fabs(tr_atransb - fock.trace_atransb(dens));
  

      prod.fullmatrix(tmpfull);
      prodmaxerror=maxdiff(tmpfull,prodfull,size);
      fock.fullmatrix(tmpfull);
      fockmaxerror=maxdiff(tmpfull,fockfull,size);
      dens.fullmatrix(tmpfull);
      densmaxerror=maxdiff(tmpfull,densfull,size);
      sum.fullmatrix(tmpfull); 
      summaxerror=maxdiff(tmpfull,sumfull,size);

      fock=fockfull;
      dens=densfull;

  DoubleSparse fockds(fockfull,size,size,epsilon/size);
  DoubleSparse densds(densfull,size,size,epsilon/size);
  DoubleSparse overds(overfull,size,size,epsilon/size);

  DoubleSparse prodds(size,size);

  
  integer noperations = DoubleSparse::multiply(fockds,densds,prodds);
  prodds.truncate(epsilon / size);
  bool ok=true;
  if (prodmaxerror > epsilon)
   {
     ok=false;
     std::cout<<"The error in the product is larger than expected"<<std::endl;
   }
  if (fockmaxerror > epsilon)
    {
      ok=false;
      std::cout<<"The error in the Fock matrix is larger than expected"<<std::endl;
    }
  if (densmaxerror > epsilon)
    {
      ok=false;
      std::cout<<"The error in the density matrix is larger than expected"<<std::endl;
    }
  if (summaxerror > epsilon)
    {
      ok=false;
      std::cout<<"The error in the sum is larger than expected"<<std::endl;
    }
  if (tracediff > epsilon)
    {
      ok=false;
      std::cout<<"The error in the trace is larger than expected"<<std::endl;
    }
  if (frobdiff > epsilon)
    {
      ok=false;
      std::cout<<"The error in the frobenius norm is larger than expected"<<std::endl;
    }
  if (tr_abdiff > epsilon)
    {
      ok=false;
      std::cout<<"The error in trace(A * B) is larger than expected"<<std::endl;
    }
  if (tr_atransbdiff > epsilon)
    {
      ok=false;
      std::cout<<"The error in trace(A' * B) is larger than expected"<<std::endl;
    }
  if (ok)
    {
      std::cout<<"OK"<<std::endl;
    }

  

  delete[] xpos;
  delete[] ypos;
  delete[] zpos;
  delete[] atomstart;
  delete[] fockfull;
  delete[] densfull;
  delete[] prodfull;
  delete[] tmpfull;
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


