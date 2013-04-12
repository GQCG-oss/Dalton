#ifndef DIVIDESPACE
#define DIVIDESPACE
#include "Permutation.h"
namespace bsm
{
class DivideSpace : public Permutation
{
 public:
  DivideSpace(const double* xp,const double* yp,
	      const double* zp,const integer* atstart,
	      integer natoms,integer maxbsize);

 private:
  /* Sort initializes atomorder, natomsperblock and nblocks */
  void sort(integer begin,integer end);
  void quicksort(const double* pos,integer low,integer high);
};
}
#endif
