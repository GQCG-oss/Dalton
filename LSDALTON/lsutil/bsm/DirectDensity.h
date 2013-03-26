#ifndef DIRECTDENSITY
#define DIRECTDENSITY

/* The implementation of the trace correcting algorithm.
   (C) Emanuel Rubensson, 2004.
*/
#include "BlockSparse.h"
namespace bsm
{
  class DirectDensity
    {
      BlockSparse& A;
      BlockSparse& B;
      BlockSparse tmp; /* tmp not necessary when multiply fully implemented */
      
      BlockSparse* present;
      BlockSparse* next;
      real lmin;
      real lmax;
      integer nshells;
      real errorlimit;
      real tracediff;
      integer iterations;
      integer ntight; /* number of tightening double-iterations */
      
      void init();
      inline void swappointers()
	{
	  BlockSparse* tmp=present;
	  present=next;
	  next=tmp;
	}
      integer step();
    public:
      DirectDensity(BlockSparse& fock,BlockSparse& dens,real lmin,real lmax,
		    integer nshells, real errorlimit=0.1, integer ntight=6);
      void finddensity(integer *iterations, real *trdiff);
      integer iter()
	{return iterations;}
    };


void fulltopacked(const real* full,real* packed,const integer size);
void packedtofull(const real* packed,real* full,const integer size);
void tripackedtofull(const real* packed,real* full,const integer size);
void gerschgorin(const real* matris,const integer size,real& low, real& high);
}
#endif
