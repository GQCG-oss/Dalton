#ifndef BLOCK
#define BLOCK
#include "Matrix.h"

#include <cmath>
#include <stdio.h>
#include <string.h>

namespace bsm
{
  /* Mycket mindre funktionalitet bor ligga i Block. */
  /* Det blir problem när viss information finns i BlockSparse */
  /* eller Sparse och viss info i Block. Battre att en skickar info */
  /* till den andra? */

class Block:public Matrix
{
 private:
  real* elements;          /* Elements sorted in column order if the matrix  */
                           /* is not transposed, row order if transposed     */
  static real cutoffvalue; /* Cutoff value used for skipping multiplication  */
                           /* of small blocks                                */
  static float seconds;    /* Used only for benchmarking                     */
  real maxvalue;           /* Largest absolute element in block              */
  real onenorm;            /* The 1-norm of the block computed by norm1      */
  real part; /*delta*/     /* Part of the entire matrix Example: entire      */
                           /* matrix(10x10) thisBlock(2x3) part=2/10=1/5     */
  char transposed;         /* 'T' for transposed and 'N'for not transposed   */

  void updateMaxAndNorm(void);
 public:
  Block(const integer* colat,integer nrofcolat,const integer* rowat,integer nrofrowat,
	const integer* atomstart,const integer* bfperatom,const integer natoms,
	const real* full,integer fullnrows,integer fullncols);
  void replace(const integer* colat,integer nrofcolat,const integer* rowat,integer nrofrowat,
	       const integer* atomstart,const integer* bfperatom,const integer natoms,
	       const real* full,integer fullnrows,integer fullncols);
  Block(integer rows=0, integer cols=0);
  Block(integer rows,integer cols,real* el);
  Block(const Block& b);
  Block(const SM<Block,real >& sm);
  Block(const SMM<Block,real >& smm);
  Block(const SMMpSM<Block,real >& smmpsm);
  Block(const SMpSM<Block, real>& smpsm);
  const Block& operator=(const Block& b);
  const Block& operator=(const SM<Block,real>& sm);
  bool operator=(const SMM<Block,real >& smm);
  bool operator+=(const SMM<Block,real >& smm);
  const Block& operator*=(real scal);
  ~Block();

  void zero(void);
  void identity(void);
  /* Adds a real number to all diagonal elements, only square matrices */
  void addscaleye(real scal);

  /* Temporary addition to be implemented in same style as mult */
  /* *this = a -> a = alpha * a + beta * b */
  void add(const Block& b, const real alpha = 1, const real beta = 1);
 

  /* Trace is at the moment only implemented for quadratic matrices */
  real trace() const;

  /* Returns trace(*this * b) */
  real trace_ab(const Block& b) const;

  /* Returns trace((*this)' * b) */
  real trace_atransb(const Block& b) const;

  /* Computes the 1-norm of the block   */
  real norm1();

  /* Computes the squared frobenius norm of the block */
  real frob_squared() const;
  /* finds the max element of the block */
  real max() const;
  /* Computes the squared frobenius norm of off-diagonal elements */
  real sum_outdia() const;
  /* Stores this block on file in txt-format  */
  void storefullonfile(char* file);

  /* Called by BlockSparse::fullmatrix to convert matrix to full matrix */
  void fullmatrix(const integer* colat,integer nrofcolat,const integer* rowat,
		  integer nrofrowat,const integer* atomstart,
		  const integer* bfperatom,const integer natoms,
		  real* full,integer fullnrows,integer fullncols);
  void extract_diag(const integer* colat, integer nrofcolat, const integer* rowat,
                    integer nrofrowat, const integer* atomstart,
                    const integer* bfperatom, const integer natoms,
                    real* diag, integer fullnrows, integer fullncols) const;
  void extract_diag_intern(real *diag);
  void precond_ao_12(integer rowoff, integer coloff,
		     const real *fup, const real *fuq,
		     real omega);
  void precond_ao_other(integer rowoff, integer coloff,
			const real *fup, const real *fuq,
			const real *du, real omega);
  void zero_lower_half(); /**< lower half and diag */
  void zero_upper_half(); /**< upper half */
  void scal_dia(real a); /**< multiply the diagonal by a */
  static void setcutoffvalue(real newvalue)
    {cutoffvalue=newvalue;}
  static void zeroseconds()
    {seconds=0;}
  static float secondstaken()
    {return seconds;}
  real maxabs()
    {return maxvalue;} 

/* Does not compute the norm!! Only returns stored value computed by norm1() */
  real norm() 
    {return onenorm;}
  real getpart()
    {return part;}

  friend void transpose(const Block& A,Block& AT);
  void transpose();
  void print() const;
  char is_transposed() const { return transposed; }
  bool write(void *unit,
             void (*write_int) (void *unit, integer *cnt, integer *where),
             void (*write_real)(void *unit, integer *cnt, real *where));
  bool read(void *unit,
            void (*write_int) (void *unit, integer *cnt, integer *where),
            void (*write_real)(void *unit, integer *cnt, real *where));
};
 
/* FIXME: is the API broken? it seems to require copying... */
inline void transpose(const Block& A,Block& AT) {
#if 0
  printf("Block::transpose. sizes: A: [%d %d]%c AT: [%d %d]%c\n",
	 A.nrofcols, A.nrofrows, A.transposed,
	 AT.nrofcols, AT.nrofrows, AT.transposed);
#endif
  if(&A != &AT) {} // AT.elements = A.elements; POINTER COPY!?
  else puts("passed same array to transpose");
  AT.maxvalue    = A.maxvalue;
  AT.onenorm     = A.onenorm;
  AT.part        = A.part;
  if(AT.nrofcols*AT.nrofrows != A.nrofrows*A.nrofcols) {
    delete[] AT.elements;
    AT.elements = new real[A.nrofrows*A.nrofcols];
  }
  AT.nrofcols = A.nrofrows;
  AT.nrofrows = A.nrofcols;
  AT.transposed = 'N';
  if(A.transposed == 'N') {
    for(integer i=0; i<A.nrofcols; i++)
      for(integer j=0; j<A.nrofrows; j++)
	AT.elements[i + j*A.nrofcols] = A.elements[j + i*A.nrofrows];
  } else { /* double transpose is noop */
    memcpy(AT.elements, A.elements,
	   A.nrofrows*A.nrofcols*sizeof(real));
  }
}
 /*
std::ostream& operator<<(std::ostream& s,const Block& b);
std::istream& operator>>(std::istream& s,Block& b);
*/

}
#endif
