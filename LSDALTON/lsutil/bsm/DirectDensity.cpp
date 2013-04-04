#include "DirectDensity.h"
#include <iostream>
namespace bsm
{
  DirectDensity::DirectDensity(BlockSparse& fock,BlockSparse& dens,
			       real lmi,real lma,
			       integer nsh,real elim, integer ntightening)
    :A(fock),B(dens),tmp(dens),lmin(lmi),lmax(lma),nshells(nsh),
     errorlimit(elim),ntight(ntightening),
     present(&fock),next(&dens),iterations(0)
  {}

  void DirectDensity::finddensity(integer *cycles, real *trdiff)
  {
    init();
    integer steps=0;
    while (steps<ntight)
      {
	if (step())
	  { steps++; iterations+=2;}
	else
	  {iterations++;}
	tracediff=present->trace()-nshells;
      }
    if(cycles) *cycles = iterations;
    if(trdiff) *trdiff = tracediff;
  }

  void DirectDensity::init()
  {
    A.addscaleye(-lmax);
    A*=-1/(lmax-lmin);
    tracediff=A.trace()-nshells;
    iterations=0;
  }

  integer DirectDensity::step()
  {
    if (fabs(tracediff)>errorlimit)
      {
	if (tracediff>0)
	  {
	    BlockSparse::multiply(*present,*present,*next);
	  }
	else
	  {
	    tmp=*present;
	    tmp.addscaleye(-2);
	    BlockSparse::multiply(tmp,*present,*next,-1);
	  }
	swappointers();
	return 0;
      }
    else
      {
	if (tracediff>0)
	  {
	    BlockSparse::multiply(*present,*present,*next);
	    swappointers();
	    tmp=*present;
	    tmp.addscaleye(-2);
	    BlockSparse::multiply(tmp,*present,*next,-1);
	  }
	else
	  {
	    tmp=*present;
	    tmp.addscaleye(-2);
	    BlockSparse::multiply(tmp,*present,*next,-1);
	    swappointers();
	    BlockSparse::multiply(*present,*present,*next);
	  }
	swappointers();
	return 1;
      }
  }


void fulltopacked(const real* full, real* packed, const integer size)
{
  integer pind=0;
  for (integer col=0;col<size;col++)
      for(integer row=0;row<=col;row++)
	  packed[pind++]=full[col*size+row];
}

void packedtofull(const real* packed, real* full, const integer size)
{
  integer psize=(size+1)*size/2;
  integer col=0;
  integer row=0;
  for(integer pind=0;pind<psize;pind++)
    {
      if (col==row)
	{
	  full[col*size+row]=packed[pind];
	  col++;
	  row=0;
	}
      else
	{
	  full[col*size+row]=packed[pind];
	  full[row*size+col]=packed[pind];
	  row++;
	}
    }
}

void tripackedtofull(const real* packed,real* full,const integer size)
{
integer psize=(size+1)*size/2;
  integer col=0;
  integer row=0;
  for(integer pind=0;pind<psize;pind++)
    {
      if (col==row)
	{
	  full[col*size+row]=packed[pind];
	  col++;
	  row=0;
	}
      else
	{
	  full[col*size+row]=packed[pind];
	  full[row*size+col]=0;
	  row++;
	}
    }
}

void gerschgorin(const real* matris,const integer size,real& low, real& high)
{
  real abssum=0;
  for (integer row=1;row<size;row++) /* All rows except diagonal element */
    {
      abssum+=fabs(matris[row]);
    }
  low=matris[0]-abssum;
  high=matris[0]+abssum;
  real tmp;
  
  for (integer col=1;col<size;col++)
    {
      abssum=0;
      for (integer row=0;row<col;row++)
	{
	  abssum+=fabs(matris[col*size+row]);
	}
      for (integer row=col+1;row<size;row++)
	{
	  abssum+=fabs(matris[col*size+row]);
	}
      tmp=matris[col*(size+1)]-abssum;
      low=(tmp<low ? tmp : low);
      tmp=matris[col*(size+1)]+abssum;
      high=(tmp>high ? tmp : high);
    }
}
}
