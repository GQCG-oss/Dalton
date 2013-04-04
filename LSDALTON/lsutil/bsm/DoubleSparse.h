#ifndef DOUBLESPARSE
#define DOUBLESPARSE
#include "Sparse.h"
namespace bsm
{
  inline void transpose(const double& A,double& AT)
    { AT = A; }
class DoubleSparse:public Sparse<double>
{
 public:
  DoubleSparse(const DoubleSparse &A)
    :Sparse<double>(A)
    {}
  DoubleSparse(const Sparse<double> &A)
    :Sparse<double>(A)
    {}
  DoubleSparse(integer rows,integer cols,integer cap=0)
    :Sparse<double>(rows,cols,cap)
    {}
  DoubleSparse(char* file,integer rows,integer cols,double thres=0,integer cap=0); 
  DoubleSparse(double* full,integer rows,integer cols,double thres=0, integer cap=0);
  void storefullonfile(char* file);
  //  Sparse<T> operator*(const Sparse<T> &m) const;
  void truncate(double eps);
};
}
#endif /* DOUBLESPARSE */
