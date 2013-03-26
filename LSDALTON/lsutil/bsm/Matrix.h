//struct SM;
//struct MM;
//struct SMM;
//struct SMMpSM;
//struct MMpM;
#ifndef MATRIX
#define MATRIX
#include "general.h"
namespace bsm
{
class Matrix
  {
  protected:
    integer nrofrows; /* nrofrows */
    integer nrofcols; /* nrofcols */
  public:
    Matrix(integer rows, integer cols)
      :nrofrows(rows),nrofcols(cols)
      {
	try{if (rows<0 || cols<0){throw InputF("Matrix:Dimensions must be positive");}}
	catch(...){throw;}
      }
    integer getnrows() const
      {return nrofrows;}
    integer getncols() const
      {return nrofcols;}

    virtual ~Matrix()
      {}
  };


template<class MAT, class SCAL>
struct SM
{
  const MAT& A;
  const SCAL alpha;
  
  SM(const MAT& AA, const SCAL scalar=1)
    :A(AA),alpha(scalar)
  {}

}; 

template<class MAT, class SCAL>
  inline SM<MAT,SCAL> operator*(const SCAL scalar, const MAT& A)
{
  return SM<MAT,SCAL>(A,scalar);
}

template<class MAT, class SCAL>
struct SMpSM {
  const SCAL alpha, beta;
  const MAT& A, B;
  SMpSM(const MAT& AA, const MAT& BB,
        const SCAL scalar_a=1, const SCAL scalar_b=1)
    :A(AA), B(BB), alpha(scalar_a), beta(scalar_b)
  {}
};

template<class MAT, class SCAL>
  inline SMpSM<MAT,SCAL> operator+(const SM<MAT,SCAL> sm1, const SM<MAT,SCAL> sm2 ) {
  return SMpSM<MAT,SCAL>(sm1.A, sm2.A, sm1.alpha, sm2.alpha);
}

template<class MAT>
struct MM
  {
  const MAT& A;
  const MAT& B;
  
  MM(const MAT& AA,const MAT& BB)
  :A(AA),B(BB)
    {}
  };

template<class MAT>
inline MM<MAT> operator*(const MAT& A, const MAT& B)              
{ 
  return MM<MAT>(A,B);
}

template<class MAT>
struct MpM
  {
  const MAT& A;
  const MAT& B;
  
  MpM(const MAT& AA,const MAT& BB)
  :A(AA),B(BB)
    {}
  };

template<class MAT>
inline MpM<MAT> operator+(MAT& A, MAT& B)              
{ 
  return MpM<MAT>(A,B);
}

template<class MAT>
struct MmM
  {
  const MAT& A;
  const MAT& B;
  
  MmM(const MAT& AA,const MAT& BB)
  :A(AA),B(BB)
    {}
  };

template<class MAT>
inline MmM<MAT> operator-(const MAT& A, const MAT& B)              
{ 
  return MmM<MAT>(A,B);
}

template<class MAT, class SCAL>
struct SMM
  {
  const SCAL alpha;
  const MAT& A;
  const MAT& B;
  SMM(const MM<MAT>& mm)
      :A(mm.A),B(mm.B),alpha(1)
    {}
  SMM(const MAT& AA,const MAT& BB, const SCAL a=1)
  :A(AA),B(BB),alpha(a)
    {}
  };

template<class MAT, class SCAL>
inline SMM<MAT,SCAL> operator*(const SM<MAT,SCAL>& SM,const MAT& B)
{
  return SMM<MAT,SCAL>(SM.A,B,SM.alpha); 
}

template<class MAT>
struct MMpM
{
  const MAT& A;
  const MAT& B;
  const MAT& C;

  MMpM(const MAT& AA, const MAT& BB, const MAT& CC)
  :A(AA),B(BB),C(CC)
  {} 
};

template<class MAT, class SCAL>
struct SMMpSM
{
  const SCAL alpha;
  const MAT& A;
  const MAT& B;
  const SCAL beta;
  const MAT& C;
  SMMpSM(const MMpM<MAT>& mmpm)
    :A(mmpm.A),B(mmpm.B),C(mmpm.C),alpha(1),beta(1)
  {}
  SMMpSM(const MAT& AA, const MAT& BB, const MAT& CC, const SCAL a=1, const SCAL b=1)
  :A(AA),B(BB),C(CC),alpha(a),beta(b)
  {} 
  };

template<class MAT, class SCAL>
inline SMMpSM<MAT,SCAL> operator+(const SMM<MAT,SCAL>& smm, const SM<MAT,SCAL>& sm)
{
  return SMMpSM<MAT,SCAL>(smm.A,smm.B,sm.A,smm.alpha,sm.alpha);
}

template<class MAT, class SCAL>
inline SMMpSM<MAT,SCAL> operator+(const SMM<MAT,SCAL>& smm, const MAT& CC)
{
  return SMMpSM<MAT,SCAL>(smm.A,smm.B,CC,smm.alpha,1);
}

template<class MAT, class SCAL>
inline SMMpSM<MAT,SCAL> operator+(const MM<MAT>& mm, const SM<MAT,SCAL>& sm)
{
  return SMMpSM<MAT,SCAL>(mm.A,mm.B,sm.A,1,sm.alpha);
}




/*onodig finns redan för SMM
template<class MAT>
inline MMpM<MAT> operator+(const MM<MAT>& mm, const MAT& CC)
{
  return MMpM<MAT>(mm.A,mm.B,CC);
}*/


/*Maste ligga i arvda klassen!!*/
/*
Matrix::Matrix(const sMMmul& mm)
  :nrofrows(mm.A.nrofrows),nrofcols(mm.B.nrofcols)
{
  this.multiply(mm.A,mm.B,*this,mm.tA,mm.tB,mm.alpha,0);
}
Matrix::Matrix(const sMMmulsMadd& mm)
  :nrofrows(mm.A.nrofrows),nrofcols(mm.B.nrofcols)
{
  this->multiply(mm.A,mm.B,mm.C,mm.tA,mm.tB,mm.alpha,mm.beta);
}

*/

}
#endif
