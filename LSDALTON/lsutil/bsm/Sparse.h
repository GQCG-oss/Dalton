#ifndef SPARSE
#define SPARSE
#include "Symbolic.h"
#include <iostream>
#include <cstdlib>
/*include "omp.h"*/

 
namespace bsm
{
template<class T>
class Sparse:public Symbolic
{
 protected:
  T* elements;
    
  void requeststorage(integer newcapacity);
 public:
  Sparse(const Sparse<T> &A, bool transposed = false);
  Sparse(integer rows=0,integer cols=0,integer cap=0);
  const Sparse<T>& operator=(const Sparse<T>& A);
  //virtual void full_to_sparse(double* full, integer row, integer col, bool sym=false);
  Sparse<T> operator*(const Sparse<T> &m) const;

  /* sparse_mul computes alpha*a*b + beta*c */
  /* So far only for beta=0 */
    static integer multiply(const Sparse<T> &A, const Sparse<T> &B,Sparse<T> &C,
			  real alpha=1, real beta=0);
     Sparse<T>& operator*=(real scal);
  ~Sparse();
};

template<class T>
void Sparse<T>::requeststorage(integer newcapacity)
{
  try
    {
      if (newcapacity<0)
	{
	  throw InputF("Sparse::requeststorage:negative capacity requested");
	}
  if (newcapacity!=capacity && newcapacity>=nrofel)
    {
      integer* ritmp=new integer[newcapacity];
      T* eltmp=new T[newcapacity];

      /* Copying rowindices and elements */
      for (integer i=0;i<nrofel;i++)
	{
	  ritmp[i]=rowind[i];
	  eltmp[i]=elements[i];
	}
      delete[] rowind;
      delete[] elements;
      rowind=ritmp;
      elements=eltmp;
      capacity=newcapacity;
      
    }
    }
  catch(...)
    {throw;}
}

template<class T>
Sparse<T>::Sparse(const Sparse<T> &A, bool transposed)
  :Symbolic(transposed ? A.nrofcols : A.nrofrows,
	    transposed ? A.nrofrows : A.nrofcols,
	    A.nrofel, A.nrofel)
{
  elements=new T[capacity];
  if(transposed) {
    for(integer i=0; i<=nrofcols; i++) colpoint[i]=0;
    for(integer i=0; i<A.nrofel;  i++) colpoint[A.rowind[i]+1]++;
    
    integer tmp, tmp2 = colpoint[0];
    for(integer i=0; i<nrofcols; i++) {
      tmp = colpoint[i+1];
      colpoint[i+1] = tmp2 + colpoint[i];
      tmp2=tmp;
    }
    for(integer col=0; col<A.nrofcols; col++)  {
      for(integer ind=A.colpoint[col]; ind<A.colpoint[col+1]; ind++) {
	integer row = A.rowind[ind];
	integer indT = colpoint[row+1];
	rowind[indT] = col;
	colpoint[row+1] = indT+1;
	transpose(A.elements[ind], elements[indT]);
      }
    }
  } else {
    for(integer i=0;i<nrofcols+1;i++)
      colpoint[i]=A.colpoint[i];
  
    for (integer i=0;i<nrofel;i++) {
      rowind[i]=A.rowind[i];
      elements[i]=A.elements[i];
    }
  }
}

template<class T>
Sparse<T>::Sparse(integer rows,integer cols,integer cap)
  :Symbolic(rows,cols,cap,0)
{
  try
    {
      elements=new T[capacity];
    }
  catch(...)
    {throw;}
}

template<class T>
const Sparse<T>& Sparse<T>::operator=(const Sparse<T>& A)
{ 
  try
    {
      nrofrows=A.nrofrows;
      if(nrofcols!=A.nrofcols)
	{
	  delete[] colpoint;
	  nrofcols=A.nrofcols;
	  colpoint=new integer[nrofcols+1];
	}  
      if (capacity<A.nrofel)
	{
	  delete[] rowind;
	  delete[] elements;
	  capacity=A.nrofel;
	  rowind=new integer[capacity];
	  elements=new T[capacity];
	}
      for(integer i=0;i<nrofcols+1;i++)
	{
	  colpoint[i]=A.colpoint[i];
	}
      nrofel=A.nrofel;
      for (integer i=0;i<nrofel;i++)
	{
	  rowind[i]=A.rowind[i];
	  elements[i]=A.elements[i];
	}
    }
  catch(...)
    {throw;}
}


/*
#define LOOP_MACRO(xxx) \
for(integer j=0;j<C.nrofcols;j++) \
{ \
  tmpcol[rowA]+= xxx this->elements[indA]*B.elements[indB]; \
} 


if(var == 1)
     LOOP_MACRO( )
else
     LOOP_MACRO(var*)


if(var == 1)
{
for(integer j=0;j<C.nrofcols;j++)
{
  tmpcol[rowA]+=this->elements[indA]*B.elements[indB];
}
}
else
{
for(integer j=0;j<C.nrofcols;j++)
{
  tmpcol[rowA]+= var * this->elements[indA]*B.elements[indB];
}
}
*/

template<class T>
Sparse<T> Sparse<T>::operator*(const Sparse<T> &B) const
{
  if (this->nrofcols==B.nrofrows)
    {
      Sparse<T> C(this->nrofrows,B.nrofcols);
      integer* tmpint=new integer[C.nrofrows];
      T* tmpcol=new T[C.nrofrows];
      for(integer v=0;v<C.nrofrows;v++)
	{tmpint[v]=-1;}
      integer rowB;
      integer rowA;  
      integer rowC;
      
      for(integer j=0;j<C.nrofcols;j++)
	{
	  C.colpoint[j]=C.nrofel;
	  for(integer indB=B.colpoint[j];indB<B.colpoint[j+1];indB++)
	    {
	      rowB=B.rowind[indB];
	      for(integer indA=this->colpoint[rowB];indA<this->colpoint[rowB+1];indA++)
		{
		  rowA=this->rowind[indA];
		  if(tmpint[rowA]!=j)
		    {
		      C.rowind[C.nrofel]=rowA;
		      C.nrofel++;
		      tmpint[rowA]=j;
		      tmpcol[rowA]=this->elements[indA]*B.elements[indB];
		    }
		  else
		    {
		      tmpcol[rowA]+=this->elements[indA]*B.elements[indB];
		    }
		}
	    }
	  if(C.nrofel>C.capacity-C.nrofrows)
	    {
	      integer newcapacity=C.nrofel*C.nrofcols/(j+1)+C.nrofrows;
	      C.requeststorage(newcapacity);
	    }
	  
	  for (integer indC=C.colpoint[j];indC<C.nrofel;indC++)
	    {
	      rowC=C.rowind[indC];
	      C.elements[indC]=tmpcol[rowC];
	    }
	}
      C.colpoint[C.nrofcols]=C.nrofel;
      delete[] tmpint;
      delete[] tmpcol;
      return C;
    }
  else
    {
      std::cerr<<"Incorrect matrixdimensions for multiplication"<<std::endl;
      std::exit(1);
    }
}

/* When matrices reaches this function they should have the correct */
/* number of rows and columns and permutation                       */
template<class T>
integer Sparse<T>::multiply(const Sparse<T> &A, const Sparse<T> &B,Sparse<T> &C,
			 real alpha, real beta)
{
  try
    {
      integer count(0);
      integer total(0);
      
      if (A.nrofcols==B.nrofrows && C.nrofrows==A.nrofrows && C.nrofcols==B.nrofcols)
	{
	  if (C.capacity<C.nrofrows)
	    {C.requeststorage(C.nrofrows);}
	  C.nrofel=0;
	  if (alpha==1 && beta==0)
	    {
	      integer* tmpint=new integer[C.nrofrows];
	      T* tmpcol=new T[C.nrofrows];
	      for(integer v=0;v<C.nrofrows;v++)
		{tmpint[v]=-1;}
	      integer rowB;
	      integer rowA;  
	      integer rowC;
	      
	      for(integer j=0;j<C.nrofcols;j++)
		{
		  C.colpoint[j]=C.nrofel;
		  for(integer indB=B.colpoint[j];indB<B.colpoint[j+1];indB++)
		    {
		      rowB=B.rowind[indB];
//#pragma omp parallel for shared(A.rowind, tmpint, j, tmpcol, A.elements, B.elements, indB, C.rowind, C.nrofel, count, total)
		      for(integer indA=A.colpoint[rowB];indA<A.colpoint[rowB+1];indA++)
			{
			  rowA=A.rowind[indA];
			  if(tmpint[rowA]!=j)
			    {
			      if(tmpcol[rowA]=A.elements[indA]*B.elements[indB])
				{
//#pragma omp critical (addnrofel)
				  C.rowind[C.nrofel]=rowA;
//#pragma omp critical (addnrofel)
				  C.nrofel++;
				  tmpint[rowA]=j;
				}
			      else{count++;}/////////////
			      
			    }
			  else
			    {
			      if(!(tmpcol[rowA]+=A.elements[indA]*B.elements[indB]))
				{count++;}////////////
			    }
			  total++;///////////
			}
		    }
		  if(C.nrofel>C.capacity-C.nrofrows)
		    {
		      integer newcapacity=(integer)(C.nrofel*((double)C.nrofcols/(j+1))+C.nrofrows);
		      C.requeststorage(newcapacity);
		    }
		  
		  for (integer indC=C.colpoint[j];indC<C.nrofel;indC++)
		    {
		      rowC=C.rowind[indC];
		      C.elements[indC]=tmpcol[rowC];
		    }
		}
	      C.colpoint[C.nrofcols]=C.nrofel;
	      delete[] tmpint;
	      delete[] tmpcol;
	    }
	  else if (beta==0)
	    {
	      integer* tmpint=new integer[C.nrofrows];
	      T* tmpcol=new T[C.nrofrows];
	      for(integer v=0;v<C.nrofrows;v++)
		{tmpint[v]=-1;}
	      integer rowB;
	      integer rowA;  
	      integer rowC;
	      
	      for(integer j=0;j<C.nrofcols;j++)
		{
		  C.colpoint[j]=C.nrofel;
		  for(integer indB=B.colpoint[j];indB<B.colpoint[j+1];indB++)
		    {
		      rowB=B.rowind[indB];
		      for(integer indA=A.colpoint[rowB];indA<A.colpoint[rowB+1];indA++)
			{
			  rowA=A.rowind[indA];
			  if(tmpint[rowA]!=j)
			    {
			      if(tmpcol[rowA]=alpha*A.elements[indA]*B.elements[indB])
				{
				  C.rowind[C.nrofel]=rowA;
				  C.nrofel++;
				  tmpint[rowA]=j;
				}
			      else{count++;}/////////////
			      
			    }
			  else
			    {
			      if(!(tmpcol[rowA]+=alpha*A.elements[indA]*B.elements[indB]))
				{count++;}////////////
			    }
			  total++;///////////
			}
		    }
		  if(C.nrofel>C.capacity-C.nrofrows)
		    {
		      integer newcapacity=(integer)(C.nrofel*((double)C.nrofcols/(j+1))+C.nrofrows);
		      C.requeststorage(newcapacity);
		    }
		  
		  for (integer indC=C.colpoint[j];indC<C.nrofel;indC++)
		    {
		      rowC=C.rowind[indC];
		      C.elements[indC]=tmpcol[rowC];
		    }
		}
	      C.colpoint[C.nrofcols]=C.nrofel;
	      delete[] tmpint;
	      delete[] tmpcol;
	    }
	  else
	    {
	      throw Failure("Sparse::multiply: Not implemented for beta!=0");
	    }
	}
      else
	{
	  throw DimensionF("Sparse::multiply:Incorrect matrixdimensions for multiplication");
	} 
      return total;
      std::cout<<"count: "<<count<<" total: "<<total<<std::endl;
    }
  catch(...)
    {throw;}
}

#if 0
/* AT has the right number of columns, rows, elements and everything.
 * This of course is useless and we keep this code here for reference
 * only. */
template<class T>
void Sparse<T>::transpose(const Sparse<T> &A, Sparse<T> &AT)
{
#if 0
  printf("will transpose nrofcols=%d to %d\n", A.nrofcols, AT.nrofcols); 
  for(integer i=0;i<A.nrofel;i++) A.elements[i].report();
  puts("AT");
  for(integer i=0;i<AT.nrofel;i++) AT.elements[i].report();
#endif
  for(integer i=0; i<=AT.nrofcols; i++) AT.colpoint[i]=0;
  for(integer i=0; i<A.nrofel; i++)     AT.colpoint[A.rowind[i]+1]++;

   integer tmp;
   integer tmp2=AT.colpoint[0];
   for(integer i=0; i<AT.nrofcols;i++) {
       tmp=AT.colpoint[i+1];
       AT.colpoint[i+1]=tmp2+AT.colpoint[i];
       tmp2=tmp;
   }

   integer row;
   integer indT;
   for(integer col=0;col<A.nrofcols;col++)  {
     for(integer ind=A.colpoint[col];ind<A.colpoint[col+1];ind++) {
       row=A.rowind[ind];
       indT=AT.colpoint[row+1];
       AT.rowind[indT]=col;
       transpose(A.elements[ind],AT.elements[indT]);
       AT.colpoint[row+1]=indT+1;
     }
   }
}
#endif

template<class T>
Sparse<T>& Sparse<T>::operator*=(real scal)
{
  for (integer ind=0;ind<nrofel;ind++)
    {
      elements[ind]*=scal;
    }
  return *this;
}


template<class T>
Sparse<T>::~Sparse()
{
  delete[] elements;
}

}
#endif
