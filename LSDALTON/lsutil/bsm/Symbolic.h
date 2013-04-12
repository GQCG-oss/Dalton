#ifndef SYMBOLIC
#define SYMBOLIC
#include "Matrix.h"
namespace bsm
{
class Symbolic:public Matrix
{
 protected:
  integer* colpoint;
  integer* rowind;
  integer nrofel;
  integer capacity;
 public:
  Symbolic(integer rows=0,integer cols=0,integer cap=0,integer el=0);

  /* Writes blockstructure to file to be read by matlab   */
  /* write spy(sparse(load('file'))) in matlab            */
  void spy(char* file);
  ~Symbolic();   
};
}
#endif
