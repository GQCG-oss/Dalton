#include "Block.h"
#include "IndexF.h"
#include "InputF.h"
#include "DimensionF.h"
#include "IOF.h"
#include "genericblas.h"
#include <iostream>
#include <fstream>
#include <iomanip> /* For setprecision in fstream */
#include <ctime>
namespace bsm
{

  real Block::cutoffvalue(pow(10.0,-11));
  float Block::seconds(0);
#if 0  
  extern "C" void dgemm_(const char *ta,const char *tb,
			 const integer *n, const integer *k, const integer *l,
			 const double *alpha,const double *A,const integer *lda,
			 const double *B, const integer *ldb,
			 const double *beta,double *C, const integer *ldc);

  extern "C" void sgemm_(const char *ta,const char *tb,
			 const integer *n, const integer *k, const integer *l,
			 const float *alpha,const float *A,const integer *lda,
			 const float *B, const integer *ldb,
			 const float *beta,float *C, const integer *ldc);

  /* dense matrix multiply */
template<class T>
static void dmm(const char *ta,const char *tb,
			 const integer *n, const integer *k, const integer *l,
			 const T *alpha,const T *A,const integer *lda,
			 const T *B, const integer *ldb,
			 const T *beta,T *C, const integer *ldc)
{
  dgemm_(ta,tb,n,k,l,alpha,A,lda,B,ldb,beta,C,ldc);
}
template<>
static void dmm<float>(const char *ta,const char *tb,
			 const integer *n, const integer *k, const integer *l,
			 const float *alpha,const float *A,const integer *lda,
			 const float *B, const integer *ldb,
			 const float *beta,float *C, const integer *ldc)
{
  sgemm_(ta,tb,n,k,l,alpha,A,lda,B,ldb,beta,C,ldc);
}
#endif

/* colat,rowat: which combination of atoms to collect   */
/* atomstart: which row/col each atoms elements start   */
Block::Block(const integer* colat,integer nrofcolat,const integer* rowat,integer nrofrowat,
	     const integer* atomstart,const integer* bfperatom,const integer natoms,
	     const real* full,integer fullnrows,integer fullncols)
  :Matrix(0,0),maxvalue(0),transposed('N')
{
  try
    {
      if (fullnrows<1 || fullncols<1 || natoms<1)
	{throw InputF("Block:fullnrows<1 || fullncols<1 || natoms<1");}
      if (nrofrowat<1 || nrofcolat<1 || nrofrowat>natoms || nrofcolat>natoms)
	{throw InputF("Block:Invalid number of row- or columnatoms");}
      
      for (integer i=0;i<nrofcolat;i++)
	{
	  if (colat[i]<0 || colat[i]>=natoms)
	    {throw IndexF("Block:columnatoms not within atom range");}
	  nrofcols+=bfperatom[colat[i]];
	  //+=atomstart[colat[i]+1]-atomstart[colat[i]];
	}
      for(integer i=0;i<nrofrowat;i++)
	{
	  if (rowat[i]<0 || rowat[i]>=natoms)
	    {throw IndexF("Block:rowatoms not within atom range");}
	  nrofrows+=bfperatom[rowat[i]];
	  //+=atomstart[rowat[i]+1]-atomstart[rowat[i]];
	}
      elements=new real[nrofrows*nrofcols];
      integer rowindex=0;
      integer colindex=0;
      integer ncolel;
      integer nrowel;
      real tmp;
      /* Loop over columnatoms */
      for (integer col=0;col<nrofcolat;col++)
	{
	  ncolel=bfperatom[colat[col]];
	  //=atomstart[colat[col]+1]-atomstart[colat[col]];
	  /* Loop over rowatoms */
	  for (integer row=0; row<nrofrowat; row++)
	    {
	      nrowel=bfperatom[rowat[row]];
	      //=atomstart[rowat[row]+1]-atomstart[rowat[row]];
	      /* Loops over single atom elements */
              if(full) {
                for(integer i=0;i<ncolel;i++) {
                  for(integer j=0;j<nrowel;j++) {
                    tmp=full[(atomstart[colat[col]]+i)*fullnrows+atomstart[rowat[row]]+j];
                    elements[(i+colindex)*nrofrows+j+rowindex]=tmp;
                    maxvalue=(maxvalue>fabs(tmp) ? maxvalue : fabs(tmp));
                  }
		}
              }
	      rowindex+=nrowel;
	    }
	  rowindex=0;
	  colindex+=ncolel;
	}
      part=nrofrows/(real)fullnrows;
      onenorm=this->norm1();
    }
  catch (...)
    {
      throw;
    }
  /*
    std::cout<<"nrofrows:"<<nrofrows<<" nrofcols:"<<nrofcols<<std::endl;
  for(integer i=0;i<nrofrows*nrofcols;i++)
    {
      std::cout<<elements[i]<<' ';
    }
    std::cout<<std::endl;*/
}
  void Block::replace(const integer* colat,integer nrofcolat,const integer* rowat,integer nrofrowat,
		      const integer* atomstart,const integer* bfperatom,const integer natoms,
		      const real* full,integer fullnrows,integer fullncols)
  {
    try
    {
      if (fullnrows<1 || fullncols<1 || natoms<1)
	{throw InputF("Block::replace: fullnrows<1 || fullncols<1 || natoms<1");}
      if (nrofrowat<1 || nrofcolat<1 || nrofrowat>natoms || nrofcolat>natoms)
	{throw InputF("Block::replace: Invalid number of row- or columnatoms");}
      integer testnrows(0);
      integer testncols(0);
      for (integer i=0;i<nrofcolat;i++)
	{
	  if (colat[i]<0 || colat[i]>=natoms)
	    {throw IndexF("Block::replace: columnatoms not within atom range");}
	  testncols+=bfperatom[colat[i]];
	}
      for(integer i=0;i<nrofrowat;i++)
	{
	  if (rowat[i]<0 || rowat[i]>=natoms)
	    {throw IndexF("Block::replace: rowatoms not within atom range");}
	  testnrows+=bfperatom[rowat[i]];
	}
      if(testnrows!=nrofrows || testncols!=nrofcols)
	{throw DimensionF("Block::replace: New matrix structure does not match old one");}
      maxvalue=0;
      transposed='N';
      integer rowindex=0;
      integer colindex=0;
      integer ncolel;
      integer nrowel;
      real tmp;
      /* Loop over columnatoms */
      for (integer col=0;col<nrofcolat;col++)
	{
	  ncolel=bfperatom[colat[col]];
	  /* Loop over rowatoms */
	  for (integer row=0;row<nrofrowat;row++)
	    {
	      nrowel=bfperatom[rowat[row]];
	      /* Loops over single atom elements */
	      for(integer i=0;i<ncolel;i++)
		{
		  for(integer j=0;j<nrowel;j++)
		    {
		      tmp=full[(atomstart[colat[col]]+i)*fullnrows+atomstart[rowat[row]]+j];
		      elements[(i+colindex)*nrofrows+j+rowindex]=tmp;
		      maxvalue=(maxvalue>fabs(tmp) ? maxvalue : fabs(tmp));
		    }
		}
	      rowindex+=nrowel;
	    }
	  rowindex=0;
	  colindex+=ncolel;
	}
      part=nrofrows/(real)fullnrows;
      onenorm=this->norm1();
    }
  catch (...)
    {
      throw;
    }
  }
Block::Block(integer rows, integer cols)
  :Matrix(rows,cols),maxvalue(0),onenorm(0),transposed('N')
{
  try
    {
      elements=new real[nrofrows*nrofcols];
    }
  catch(...)
    {throw;}
}


Block::Block(integer rows, integer cols,real* el)
  :Matrix(rows,cols),maxvalue(0),transposed('N')
{
  try
    {
      elements=new real[nrofrows*nrofcols];  
      for (integer i=0;i<nrofrows*nrofcols;i++)
	{
	  elements[i]=el[i];
	  maxvalue=(maxvalue>fabs(el[i]) ? maxvalue : fabs(el[i]));
	}
      onenorm=this->norm1();
    }
  catch(...)
    {throw;}
}


Block::Block(const Block& b)
  :Matrix(b.nrofrows,b.nrofcols),maxvalue(b.maxvalue),
   onenorm(b.onenorm),part(b.part),transposed(b.transposed)
{
try
  {
    elements=new real[nrofrows*nrofcols];
    for (integer i=0;i<nrofrows*nrofcols;i++)
      {
	elements[i]=b.elements[i];
      }
  }
catch(...)
  {throw;}
}


Block::Block(const SM<Block,real>& sm)
  :Matrix(sm.A.nrofrows,sm.A.nrofcols),maxvalue(fabs(sm.alpha)*sm.A.maxvalue),
     onenorm(fabs(sm.alpha)*sm.A.onenorm),part(sm.A.part),
   transposed('N')
{
try
  {
    elements=new real[nrofrows*nrofcols];
    if(sm.A.transposed == 'N') {
      for(integer i=0;i<nrofrows*nrofcols;i++)
	elements[i]=sm.alpha*sm.A.elements[i];
    } else {
      for(integer i=0; i<nrofcols; i++)
        for(integer j=0; j<nrofrows; j++)
          elements[i + j*nrofrows] = 
            sm.alpha*sm.A.elements[j + i*nrofrows];
    }
  }
catch(...)
  {throw;}
}




Block::Block(const SMM<Block,real>& smm)
  :Matrix(smm.A.nrofrows,smm.B.nrofcols),
   /* Upper bound estimate used for maxvalue */
     maxvalue(fabs(smm.alpha)*smm.A.maxvalue*smm.B.maxvalue*smm.A.nrofcols),
     part(smm.A.part),transposed('N')
{
  try
    {
      if (smm.A.nrofcols==smm.B.nrofrows)
	{
	  elements=new real[nrofrows*nrofcols];
	  integer LDA=(smm.A.transposed=='n' || smm.A.transposed=='N' ?
		   smm.A.nrofrows : smm.A.nrofcols);
	  integer LDB=(smm.B.transposed=='n' || smm.B.transposed=='N' ?
		   smm.B.nrofrows : smm.B.nrofcols);
	  clock_t start=clock();
	  gemm(&smm.A.transposed,&smm.B.transposed,nrofrows,nrofcols,
		 smm.A.nrofcols,&smm.alpha,
		 smm.A.elements,LDA,
		 smm.B.elements,LDB,
		 &ZERO,elements,nrofrows);
	  seconds+=((float)(clock()-start))/CLOCKS_PER_SEC;
	  onenorm=this->norm1();
	}
      else
	{
	  throw DimensionF("Block:Incorrect matrixdimensions for multiplication");
	}
    }
  catch(...)
    {throw;}
}

Block::Block(const SMpSM<Block, real>& smpsm) 
  :Matrix(smpsm.A.nrofrows,smpsm.A.nrofcols),part(smpsm.A.part),transposed('N') {
  try {
    if (smpsm.A.nrofrows == smpsm.B.nrofrows &&
        smpsm.A.nrofcols == smpsm.B.nrofcols) {
      elements=new real[nrofrows*nrofcols];
      if ((smpsm.A.transposed=='N' || smpsm.A.transposed=='n') &&
          (smpsm.B.transposed=='N' || smpsm.B.transposed=='n')) {
        maxvalue=0;
        for (integer i=0; i<nrofrows*nrofcols; i++) {
          elements[i] = smpsm.alpha * smpsm.A.elements[i] + 
            smpsm.beta * smpsm.B.elements[i];
          maxvalue=(maxvalue>elements[i] ? maxvalue : elements[i]);
        }
      } else {
        maxvalue=0;
        for (integer i=0; i<nrofcols; i++) {
          for (integer j=0; j<nrofrows; j++) {
            real a = smpsm.alpha *
              smpsm.A.elements[smpsm.A.transposed == 'N'
                               ?  j + i*nrofrows : i + j*nrofcols ];
            real b = smpsm.beta *
              smpsm.B.elements[smpsm.B.transposed == 'N'
                               ? j + i*nrofrows : i + j*nrofcols ];
            elements[j+i*nrofrows] = a + b;
            maxvalue=(maxvalue>elements[i] ? maxvalue : elements[i]);
          }
        }
	//throw Failure("Block: Addition not implemented for transposed matrices");
      }
      onenorm=this->norm1();
    }
    else {
      throw DimensionF("Block:Incorrect matrixdimensions for addition");
    }
  }
  catch(...)
    {throw;}
}



const Block& Block::operator=(const Block& b)
{
  try
    {
      if (nrofrows*nrofcols!=b.nrofrows*b.nrofcols)
	{
	  delete[] elements;
	  elements=new real[b.nrofrows*b.nrofcols];
	}
      nrofrows=b.nrofrows;
      nrofcols=b.nrofcols;
      for (integer i=0;i<nrofrows*nrofcols;i++)
	{
	  elements[i]=b.elements[i];
	}
      maxvalue=b.maxvalue;
      onenorm=b.onenorm;
      part=b.part;
      transposed=b.transposed;
      return *this;
    }
  catch(...)
    {throw;}
}

const Block& Block::operator=(const SM<Block,real>& sm)
{
  if (nrofrows*nrofcols!=sm.A.nrofrows*sm.A.nrofcols) {
    delete[] elements;
    elements=new real[sm.A.nrofrows*sm.A.nrofcols];
  }
  nrofrows = sm.A.nrofrows;
  nrofcols = sm.A.nrofcols;
  transposed = 'N';
  if(sm.A.transposed == 'N') {
    for (integer i=0; i<nrofrows*nrofcols; i++)
      elements[i] = sm.alpha*sm.A.elements[i];
  } else {
    for(integer i=0; i<nrofcols; i++)
      for(integer j=0; j<nrofrows; j++)
	elements[j + i*nrofrows] = sm.alpha*sm.A.elements[i + j*nrofrows];
  }
  maxvalue = sm.A.maxvalue;
  onenorm  = sm.A.onenorm;
  part     = sm.A.part;
  return *this;
}


bool Block::operator=(const SMM<Block,real>& smm)
{
  try
    {
      /* Upper bound estimate used for maxvalue */
      maxvalue=fabs(smm.alpha)*smm.A.maxvalue*smm.B.maxvalue*smm.A.nrofcols;
      
      /* Perform matrixoperation only if the result has significantly large values */
      if (fabs(smm.alpha)*smm.A.onenorm*smm.B.onenorm >
	  cutoffvalue*smm.A.part*smm.B.part)
      //if (maxvalue>cutoffvalue)
	{
	  part=smm.A.part;
	  transposed='N';
	  integer cap=nrofrows*nrofcols;
	  nrofrows=smm.A.nrofrows;
	  nrofcols=smm.B.nrofcols;
	  
	  if (smm.A.nrofcols==smm.B.nrofrows)
	    {
	      if(nrofrows*nrofcols!=cap)
		{
		  delete[] elements;
		  elements=new real[nrofrows*nrofcols];
		}
	      integer LDA=(smm.A.transposed=='n' || smm.A.transposed=='N' ?
		       smm.A.nrofrows : smm.A.nrofcols);
	      integer LDB=(smm.B.transposed=='n' || smm.B.transposed=='N' ?
		       smm.B.nrofrows : smm.B.nrofcols);
	      clock_t start=clock();
	      gemm(&smm.A.transposed,&smm.B.transposed, nrofrows, nrofcols,
		     smm.A.nrofcols,&smm.alpha,
		     smm.A.elements, LDA,
		     smm.B.elements, LDB,
		     &ZERO,elements, nrofrows);
	      seconds+=((float)(clock()-start))/CLOCKS_PER_SEC;
	      onenorm=this->norm1();
	    }
	  else
	    {
	      throw DimensionF("Block::operator=: Incorrect matrixdimensions for multiplication");
	    }
	  return true;
	}
      else
	{
	  return false;
	}
    }
  catch(...)
    {throw;}
}


bool Block::operator+=(const SMM<Block,real>& smm)
{
  try
    {
      /* Upper bound estimate used for tmpmax */
      real tmpmax=fabs(smm.alpha)*smm.A.maxvalue*smm.B.maxvalue*smm.A.nrofcols;
      if (fabs(smm.alpha)*smm.A.onenorm*smm.B.onenorm >
	  cutoffvalue*smm.A.part*smm.B.part)
	//if (tmpmax>cutoffvalue)
	{
	  part=smm.A.part;
	  maxvalue+=tmpmax;
	  
	  if (smm.A.nrofcols==smm.B.nrofrows && smm.A.nrofrows==nrofrows &&
	      smm.B.nrofcols==nrofcols)
	    {
	      if(transposed=='N' || transposed=='n')
		{
		  integer LDA = (smm.A.transposed=='n' || smm.A.transposed=='N') 
                    ? smm.A.nrofrows : smm.A.nrofcols;
		  integer LDB = (smm.B.transposed=='n' || smm.B.transposed=='N') 
                    ? smm.B.nrofrows : smm.B.nrofcols;

		  clock_t start=clock();
                  //printf("   bsize %d %d\n", smm.B.nrofrows, smm.B.nrofcols);
		  gemm(&smm.A.transposed,&smm.B.transposed,
			 nrofrows,nrofcols,
			 smm.A.nrofcols,&smm.alpha,
			 smm.A.elements, LDA,
			 smm.B.elements, LDB,
			 &ONE,elements, nrofrows);
		  seconds+=((float)(clock()-start))/CLOCKS_PER_SEC;
		}
	      else
		{
		  /* Instead of solving something like C'=AB+C' */
		  /* we solve C=B'A'+C                          */
		  char tA,tB;
		  integer LDA,LDB;
		  if (smm.A.transposed=='n' || smm.A.transposed=='N')
		    {
		      tA='T';
		      LDA=smm.A.nrofrows;
		    }
		  else
		    {
		      tA='N';
		      LDA=smm.A.nrofcols;
		    }
		  if (smm.B.transposed=='n' || smm.B.transposed=='N')
		    {
		      tB='T';
		      LDB=smm.B.nrofrows;
		    }
		  else
		    {
		      tB='N';
		      LDB=smm.B.nrofcols;
		    }
		  clock_t start=clock();
		  gemm(&tB,&tA, nrofcols, nrofrows,
		       smm.A.nrofcols,&smm.alpha,
		       smm.B.elements, LDB,
			 smm.A.elements, LDA,
			 &ONE,elements, nrofrows);
		  seconds+=((float)(clock()-start))/CLOCKS_PER_SEC;
		}

	    }
	  else
	    {
	      throw DimensionF("Block:operator+=: Incorrect matrixdimensions for multiplication");
	    }
	  onenorm=this->norm1();
	  return true;
	}
      else
	{
	  return false;
	}
    }
  catch(...)
    {throw;}
}

const Block& Block::operator*=(real scal)
{
      maxvalue*=scal;
      onenorm*=fabs(scal);
      for (integer ind=0;ind<nrofrows*nrofcols;ind++)
	{
	  elements[ind]*=scal;
	}
      return *this;
}

Block::~Block()
{
  delete[] elements;
}

void Block::zero(void) 
{
  for(integer el=0; el<nrofrows*nrofcols; el++)
    elements[el] = 0.0;
}

void Block::identity(void) 
{
  for(integer col=0; col<nrofcols; col++) {
    for(integer row=0; row<col; row++)
      elements[row + col*nrofrows] = 0.0;
    elements[col + col*nrofrows] = 1.0;
    for(integer row=col+1; row<nrofrows; row++)
      elements[row + col*nrofrows] = 0.0;
  }
  onenorm = maxvalue = 1.0;
}

void Block::addscaleye(real scal)
{
try
    {
      if(nrofcols==nrofrows)
	{
	  for (integer i=0;i<nrofcols*nrofrows;i=i+nrofrows+1)
	    {
	      elements[i]+=scal;
	    }
	  onenorm=this->norm1();
	}
      else
	{
	  throw Failure("Block:addscaleye: Not implemented for rectangular matrices");
	}
    }
  catch(...)
    {throw;}
}

void Block::add(const Block& b, const real alpha, const real beta)
{
  try
    {
      if (nrofcols==b.nrofcols && nrofrows==b.nrofrows)
	{
	  if ((transposed=='N' || transposed=='n') && (b.transposed=='N' || b.transposed=='n'))
	    {
	      maxvalue=0;
	      for (integer i=0;i<nrofrows*nrofcols;i++)
		{
		  elements[i] = alpha * elements[i] + beta * b.elements[i];
		  maxvalue=(maxvalue>elements[i] ? maxvalue : elements[i]);
		}
	    }
	  else
	    {
	      throw Failure("Block:add: Addition not implemented for transposed matrices");
	    }
	  onenorm=this->norm1();
	}
	else
	  {
	    throw DimensionF("Block:add: Incorrect matrixdimensions for addition");
	  }
    }
  catch(...)
    {throw;}
}

real Block::trace() const
{
  try
    {
      if(nrofcols==nrofrows)
	{
	  real tr=0;
	  for (integer i=0;i<nrofcols*nrofrows;i=i+nrofrows+1)
	    {
	      tr+=elements[i];
	    }
	  return tr;
	}
      else
	{
	  throw Failure("Block:trace: Trace is not implemented for rectangular matrices");
	}
    }
  catch(...)
    {throw;}
}


real Block::trace_ab(const Block& b) const {
  /* First equality must be valid for multiplication */
  if (nrofcols != b.nrofrows) {
    printf("trace_ab:A %d %d\n", nrofcols, b.nrofrows);
    throw DimensionF("Block:trace_ab: Incorrect matrix dimensions "
                     "for multiplication");
  }
  if (nrofrows != b.nrofcols) {
    printf("trace_ab:B: %d %d\n", nrofrows, b.nrofcols);
    throw Failure("Block:trace_ab: Trace is not implemented for "
                  "rectangular matrices");
  }
  if (transposed == b.transposed) {
    /* Trace for rectangular matrices undefined */
    real tr = 0;
    for (integer i = 0 ; i < b.nrofcols; i++) {
      for (integer j = 0 ; j < b.nrofrows; j++) {
        tr += b.elements[i * b.nrofrows + j] * elements[j * nrofrows + i];
      }
    }
    return tr;
  } else {
    real tr = 0;
    for (integer i = 0 ; i < b.nrofcols*b.nrofrows; i++)
      tr += b.elements[i] * elements[i];
    return tr;
  }
}


real Block::trace_atransb(const Block& b) const {
    /* First equality must be valid for multiplication */
  if (nrofrows != b.nrofrows)
    throw DimensionF("Block:trace_atransb: Incorrect matrix "
                     "dimensions for multiplication");
  if (nrofcols != b.nrofcols)
    throw Failure("Block:trace_atransb: Trace is not implemented "
                  "for rectangular matrices");
  
  if ((transposed=='N'   || transposed=='n') &&
      (b.transposed=='N' || b.transposed=='n')) {
    /* Trace for rectangular matrices undefined */
    real tr = 0;
    for (integer i = 0 ; i < b.nrofcols*b.nrofrows; i++)
      tr += b.elements[i] * elements[i];
    return tr;
  } else {
    if((transposed=='N'   || transposed=='n') &&
       (b.transposed=='T' || b.transposed=='t')) {
      real tr = 0;
      for (integer i = 0 ; i < nrofcols; i++)
        for (integer j = 0 ; j < nrofrows; j++)
          tr += elements[j + i*nrofrows] * b.elements[i + j*nrofrows];
      return tr;
    } else {
      printf("%c %c\n", transposed, b.transposed);
      throw Failure("Block:add: trace_atransb not implemented for "
                    "transposed matrices");
    }
  }
}


real Block::norm1()
{
  try
    {
      real res=0;
      real tmp;
      if (transposed=='N' || transposed=='n')
	{
	  for(integer col=0;col<nrofcols;col++)
	    {
	      tmp=0;
	       /* BLAS routine dasum can be used instead */
	      for(integer row=0;row<nrofrows;row++)
		{
		 
		  tmp+=fabs(elements[col*nrofrows+row]);
		}
	      res=(res>tmp ? res : tmp);
	    }
	}
      else
	{
	  for(integer col=0;col<nrofcols;col++)
	    {
	      tmp=0;
	       /* BLAS routine dasum can be used instead */
	      for(integer row=0;row<nrofrows;row++)
		{
		  tmp+=fabs(elements[col+row*nrofcols]);
		}
	      res=(res>tmp ? res : tmp);
	    }
	}
      return res;
    }
  catch(...)
    {throw;}
}

real Block::frob_squared() const {
  real res = 0;
  for (integer i = 0; i < nrofrows * nrofcols; i++) {
    res += elements[i] * elements[i];
  }
  return res;
}

real Block::max() const {
  real res = elements[0];
  for (integer i = 1; i < nrofrows * nrofcols; i++) {
    if(elements[i]>res)
      res = elements[i];
  }
  return res;
}

real Block::sum_outdia() const {
  real res = 0;
  for (integer i = 0; i < nrofrows * nrofcols; i++) {
    if(i != i+nrofcols*i)
      res += elements[i] * elements[i];
  }
  return res;
}

void Block::storefullonfile(char* file)
{
  try
    {
      std::ofstream output(file);
      if (!output)
	{throw IOF("Block:Cannot open outputfile"+std::string(file));}
      output<<std::setprecision(10);
      for(integer i=0;i<nrofcols;i++)
	{
	  for(integer j=0;j<nrofrows;j++)
	    {
	      output<<elements[i*nrofrows+j]<<' ';
	    }
	  output<<std::endl;
	}
    }
  catch(...)
    {throw;}
}


void Block::fullmatrix(const integer* colat,integer nrofcolat,const integer* rowat,
		       integer nrofrowat,const integer* atomstart,
		       const integer* bfperatom,const integer natoms,
		       real* full,integer fullnrows,integer fullncols)
{
  try
    {
      if (fullnrows<1 || fullncols<1 || natoms<1)
	{throw InputF("Block::fullmatrix:fullnrows<1 || fullncols<1 || natoms<1");}
      if (nrofrowat<1 || nrofcolat<1 || nrofrowat>natoms || nrofcolat>natoms)
	{throw InputF("Block::fullmatrix:Invalid number of row- or columnatoms");}
      
      for (integer i=0;i<nrofcolat;i++)
	{
	  if (colat[i]<0 || colat[i]>=natoms)
	    {throw IndexF("Block::fullmatrix:columnatoms not within atom range");}
	}
      for(integer i=0;i<nrofrowat;i++)
	{
	  if (rowat[i]<0 || rowat[i]>=natoms)
	    {throw IndexF("Block::fullmatrix:rowatoms not within atom range");}
	}
      
      integer rowindex=0;
      integer colindex=0;
      integer ncolel;
      integer nrowel;
      /* Loop over columnatoms */
      for (integer col=0;col<nrofcolat;col++)
	{
	  ncolel=atomstart[colat[col]+1]-atomstart[colat[col]];
	  /* Loop over rowatoms */
	  for (integer row=0;row<nrofrowat;row++)
	    {
	      nrowel=atomstart[rowat[row]+1]-atomstart[rowat[row]];
	      /* Loops over single atom elements */
	      for(integer i=0;i<ncolel;i++)
		{
		  for(integer j=0;j<nrowel;j++)
		    {
                      if(transposed == 'T') {
                        printf("Transposed!?\n");
                        full[(atomstart[colat[col]]+i)*fullnrows+atomstart[rowat[row]]+j]=
                          elements[(j+colindex)*nrofrows+i+rowindex];
                      } else
                        full[(atomstart[colat[col]]+i)*fullnrows+atomstart[rowat[row]]+j]=
                          elements[(i+colindex)*nrofrows+j+rowindex];
		    }
		}
	      rowindex+=nrowel;
	    }
	  rowindex=0;
	  colindex+=ncolel;
	}
    }
  catch(...)
    {throw;}
}


void Block::extract_diag(const integer* colat, integer nrofcolat, const integer* rowat,
                         integer nrofrowat, const integer* atomstart,
                         const integer* bfperatom, const integer natoms,
                         real* diag, integer fullnrows, integer fullncols) const
{
  if (fullnrows<1 || fullncols<1 || natoms<1)
    {throw InputF("Block::extract_diag:fullnrows<1 || fullncols<1 || natoms<1");}
  if (nrofrowat<1 || nrofcolat<1 || nrofrowat>natoms || nrofcolat>natoms)
    {throw InputF("Block::extract_diag:Invalid number of row- or columnatoms");}

  for (integer i=0; i<nrofcolat; i++)	{
    if (colat[i]<0 || colat[i]>=natoms)
      {throw IndexF("Block::extract_diag:columnatoms not within atom range");}
  }
  for(integer i=0;i<nrofrowat;i++) {
    if (rowat[i]<0 || rowat[i]>=natoms)
      {throw IndexF("Block::extract_diag:rowatoms not within atom range");}
  }

  integer rowindex=0;
  integer colindex=0;
  integer ncolel;
  integer nrowel;

  /* Loop over columnatoms */
  for (integer col=0; col<nrofcolat; col++) {
    ncolel = atomstart[colat[col]+1]-atomstart[colat[col]];
    /* Loop over rowatoms */
    for (integer row=0; row<nrofrowat; row++) {
      nrowel = atomstart[rowat[row]+1]-atomstart[rowat[row]];
      /* Loops over single atom elements */
      if(colat[col] == rowat[row]) {
        for(integer i=0; i<ncolel; i++)
          diag[atomstart[colat[col]]+i]=
            elements[(i+colindex)*nrofrows+i+rowindex];
      }
      rowindex += nrowel;
    }
    rowindex=0;
    colindex+=ncolel;
  }
}

void Block::extract_diag_intern(real *diag)
{
  if(nrofcols != nrofrows)
    throw Failure("extract_diag_intern called for a non-square block");

  for(integer i=0; i<nrofcols; i++)
    diag[i] = elements[i + i*nrofrows];
}

void Block::updateMaxAndNorm(void)
{
  integer elcnt = nrofrows * nrofcols;
  maxvalue = fabs(elements[0]);
  for(integer i=1; i<elcnt; i++) {
    real e = fabs(elements[i]);
    if(e > maxvalue) maxvalue = e;
  }
  onenorm=this->norm1();
}

void Block::precond_ao_12(integer rowoff, integer coloff,
			  const real *fup, const real *fuq,
			  real omega)
{
  for(integer c=0; c<nrofcols; c++)
    for(integer r=0; r<nrofrows; r++) {
      real denom
	= fuq[coloff+c] + fuq[rowoff+r]
	- fup[rowoff+r] - fup[coloff+c]
	- omega;
      if(fabs(denom)>1e-10)
        elements[r + c*nrofrows] /= denom;
    }
  updateMaxAndNorm();
}

void Block::precond_ao_other(integer rowoff, integer coloff,
			     const real *fup, const real *fuq,
			     const real *du, real omega)
{
  for(integer c=0; c<nrofcols; c++)
    for(integer r=0; r<nrofrows; r++) {
      real denom
        = fuq[coloff+c] + fuq[rowoff+r]
        - fup[rowoff+r] - fup[coloff+c]
        - omega*(du[coloff+c] - du[rowoff+r]);
      if(fabs(denom)>1e-10) {
        if(omega<0.5 || r != c)
          elements[r + c*nrofrows] /= denom;
      }
    }
  updateMaxAndNorm();
}

/** Zero lower half and diagonal. */
void Block::zero_lower_half() {
  for(integer c=0; c<nrofcols; c++)
    for(integer r=c; r<nrofrows; r++) 
      elements[r + c*nrofrows] = 0;
  updateMaxAndNorm();
}
/** Zero upper half of the block. */
void Block::zero_upper_half() {
  for(integer c=1; c<nrofcols; c++)
    for(integer r=0; r<c; r++) 
      elements[r + c*nrofrows] = 0;
  updateMaxAndNorm();
}

/** Multiply the diagonal by a */
void Block::scal_dia(real a) {
  for(integer c=0; c<nrofcols; c++)
      elements[c + c*nrofrows] *= a;
  updateMaxAndNorm();
}

void Block::transpose() {
    transposed = (transposed == 'n' || transposed == 'N')
      ? 'T' : 'N';
    integer tmp = nrofrows;
    nrofrows = nrofcols;
    nrofcols = tmp;
    onenorm=this->norm1();
  }


void Block::print() const {
    integer nrow, ncol;
    if(transposed == 'T') {
      nrow = nrofcols; ncol = nrofrows;
    } else {
      nrow = nrofrows; ncol = nrofcols;
    }
    for(integer i=0; i<nrow; i++) {
      for(integer j=0; j<ncol; j++) printf("%10.5f ", elements[j*nrow + i]);
      puts("");
    }                                  
  }

bool Block::write(void *unit,
                  void (*write_int) (void *unit, integer *cnt, integer *where),
                  void (*write_real)(void *unit, integer *cnt, real *where))
{
  integer dt[2];
  static integer TWO = 2, THREE = 3;

  if(transposed != 'N')
    throw Failure("saving transposed blocks is not implemented");
  dt[0] = nrofrows; dt[1] = nrofcols;
  write_int (unit, &TWO,   dt);
  write_real(unit, &THREE, &maxvalue);
  dt[0] = nrofrows*nrofcols; write_real(unit, &dt[0], elements);
  return true;
}

bool Block::read(void *unit,
                 void (*read_int) (void *unit, integer *cnt, integer *where),
                 void (*read_real)(void *unit, integer *cnt, real *where))
{
  integer dt[2];
  static integer TWO = 2, THREE = 3;
  read_int (unit, &TWO,   dt);
  read_real(unit, &THREE, &maxvalue);
  transposed = 'N';
  if(dt[0]*dt[1] != nrofrows*nrofcols) {
    delete[] elements;
    elements = new real[ dt[0]*dt[1] ];
  }
  nrofrows = dt[0];
  nrofcols = dt[1];
  dt[0] = nrofrows*nrofcols; read_real(unit, &dt[0], elements);

  updateMaxAndNorm();
  return true;
}

} /* end of bsm namespace */
