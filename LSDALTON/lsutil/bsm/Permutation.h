#ifndef PERMUTATION
#define PERMUTATION
#include "general.h"
namespace bsm
{
class Permutation
{
  friend class BlockSparse;
  
 protected:
  double* xpos;          /* Coordinates for atoms [natoms]              */
  double* ypos;          /* Ordered in original order                   */
  double* zpos;  
  integer* atomorder;        /* The ordering of the atoms [natoms]          */
  integer natoms;            /* Number of atoms                             */
  integer* natomsperblock;   /* Number of atoms per block [nblocks]         */
  integer nblocks;           /* Number of blocks in row and/or column       */
  integer* blockstart;       /* Start position in atomorder for each block  */
                         /* [nblocks+1]                                 */
  integer* atomstart;        /* Start position in full matrix for each atom */
                         /* [natoms+1]                                  */
  integer* bfperatom;        /* Number of basis functions per atom [natoms] */
  integer maxblocksize;      /* Largest allowed blocksize                   */

    
 public:
  Permutation(const double* xp,const double* yp,const double* zp,
	      const integer* atomstart,integer natoms,integer maxblocksize);
  const Permutation& operator=(const Permutation& A);
 
  /* Returns  the shortest distance between two sets of atoms */
  double distance(const integer* atomset1, const integer* atomset2,
		  const integer nset1,const integer nset2);

  /* Writes matrix with distances between blocks to file */
  void distmatrixtofile(char* file);

  /* Writes permutation to file to be read by matlab      */ 
  /* Run perm.m in matlab                                 */
  void permutationtofile(char* file);

  float blocksize()
    {
      return atomstart[natoms]/(float)nblocks;
    }

  virtual ~Permutation();
};
}
#endif
