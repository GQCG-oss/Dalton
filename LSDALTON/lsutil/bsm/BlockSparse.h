#ifndef BLOCKSPARSE
#define BLOCKSPARSE
#include "Sparse.h"
#include "Block.h"
#include "Permutation.h"
namespace bsm
{
  class BlockSparse:public Sparse<Block>
    {
    protected:
      Permutation& perm;
      integer totalnrows;
      integer totalncols;
      void setup_rcdim_arrs(integer *rdim, integer *cdim);
    public:
      /* Constructors for matrices with structure as fock and dens matrices */
      /* Cutoff is assumed to be a cutoff radius if cutoff>=1               */
      /* Otherwise (0<cutoff<1) truncation using 1-norm is used             */
      BlockSparse(Permutation& perm,integer fullsize, const real* full, 
		  const double cutoff);
      BlockSparse(Permutation& perm, integer fullsize); /* INTE BRA, DENNA KONSTRUKTOR SKA TAS BORT.*/
      BlockSparse(const BlockSparse &A, bool transposed = false);

      /* addition:  *this = C = alpha * A + beta * B */
      BlockSparse(const BlockSparse& A, const BlockSparse& B,
                  const real alpha = 1, const real beta = 1);

      const BlockSparse& operator=(const BlockSparse& A);
      const BlockSparse& operator=(const real* full);

      /* Converts matrix to full matrix ie columnoriented array */
      void fullmatrix(real* full); /* Only works for fock, dens structure */
      void extract_diag(real* diag) const; /* Extract the diagonal elements */
      void extract_diag_internal(real *diag);
      void precond_ao(integer symm, const real *fup, const real *fuq,
                      const real *du, real omega);
      void zerohalf(const char *triangle);
      void scal_dia(real alpha); /**< multiply diagonal by given number */
      
      /* Returns the fraction of stored elements (not blocks) */
      float elsparsity();
      void maxelementstofile(char* file);
      void norm1ofblockstofile(char* file);
      /* The following functions should perhaps lie in Sparse?... */
      void identity(void);
      void addscaleye(real scal);
      real trace() const;
      real trace_ab(const BlockSparse& B) const;
      real trace_atransb(const BlockSparse& B) const;
      void max_abs_diag(integer *pos, real *val) const;
      /* Computes the frobenius norm of the matrix */
      real frob() const;
      /* max() finds the max element */
      real max() const;
      /* Computes the sum of nondiagonal elements of the matrix */
      real sum_outdia_sqnorm2() const;
      void print_struct() const;

      bool write(void *unit,
                 void (*write_int) (void *unit, integer *cnt, integer *where),
                 void (*write_real)(void *unit, integer *cnt, real *where));
      bool read(void *data,
                void (*read_int) (void *data, integer *cnt, integer *where),
                void (*read_real)(void *data, integer *cnt, real *where));
      integer gettotalrows() const { return totalnrows; };
      integer gettotalcols() const { return totalncols; };

    };
}


#endif
