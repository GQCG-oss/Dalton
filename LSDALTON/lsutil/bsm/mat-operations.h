
#ifdef __cplusplus
extern "C" {
#endif
  typedef double real;
  /* anything that has size != size2 is stored as plain full matrices. */
  typedef struct { integer size; integer size2; void *private_data; void *permutation;} DalMatrix;
    


  void bsm_lib_init_(double* x, double* y, double* z, integer* atomstart, 
		    integer* natoms, integer* maxblocksize, real* epsilon,
		     real *tracethr, integer *ntight);

  void bsm_lib_free_(long int* perm);
  void bsm_lib_setpermutation_(long int* perm);
  void bsm_lib_getpermutation_(long int* perm);

  void bsm_init_(DalMatrix *a, const integer *size, const integer *size2); // VKT
  void bsm_free_(DalMatrix *a);       // VKT
  //void bsm_set_from_full(real *afull, real *alpha, DalMatrix *a); // VKT

  /* bsm_init and bsm_set_from_full above are replaced with             *
   * cutoff is assumed to be a cutoff radius if cutoff>=1               *
   * Otherwise (0 < cutoff < 1) truncation using 1-norm is used         *
   * so that the error in the 1-norm is smaller than cutoff             *
   * OBSERVE: truncation using the 1-norm is recommended!               */
  void bsm_init_from_full_(DalMatrix* a, const integer* size,
                           const integer *size2, real* afull,
                           real *sparsity);

  //  void bsm_to_full(DalMatrix *a, real *alpha, real *afull);       // VKT
  /* changed to (can easily be unchanged if proven necessary)  */
  void bsm_to_full_(DalMatrix *a, real *afull, real *sparsity);
  void bsm_extract_diag_(DalMatrix *a, real *adiag);

  void bsm_print_(DalMatrix *a, integer *i_row1, integer *i_rown,
		 integer *j_col1, integer *j_coln);

  void bsm_transpose_(DalMatrix *a, DalMatrix *b);

  /* assignment operator a = b (a should be initialized) */
    void bsm_assign_(DalMatrix *a, const DalMatrix *b);              // VKT
    void bsm_copy_(real *alpha,const DalMatrix *a, DalMatrix *b);
    real bsm_tr(DalMatrix *a);
    real bsm_tratransb_(const DalMatrix *a, const DalMatrix *b);     // VKT 
    real bsm_trab_(const DalMatrix *a, const DalMatrix *b);          // VKT 
    //c = alpha*ab + beta*c
    //transa = 'T'/'t' - transposed, 'N'/'n' - normal
    void bsm_mul_(DalMatrix *a, DalMatrix *b, char *transa, char *transb,
                  real *alpha, real *beta, DalMatrix *c);            // VKT 
    
    /* C = alpha*A + beta*B                                           *
     * c should be initialized                                        */
    void bsm_add_(real *alpha, const DalMatrix *a,
                  real *beta, const DalMatrix *b,
                  DalMatrix* c);   // VKT 
    void bsm_daxpy_(real *alpha, const DalMatrix *a, DalMatrix* c);
    void bsm_identity_(DalMatrix *I); /* replaced with A := A + alpha*I */
    void bsm_add_identity_(DalMatrix* A, real* alpha);

    // Y += alpha*X;
    void bsm_daxpy(real *alpha, DalMatrix *X, DalMatrix *Y);
    //frobenius norm
    real bsm_frob(DalMatrix *a);                                 // VKT 
    //max element()
    real bsm_max_(DalMatrix *a);
    void bsm_max_diag_(DalMatrix *a, integer *pos, real *val);
    //Sqares all the outer diagonal elements in a, and return the sum
    real bsm_outdia_sqnorm2(DalMatrix *a);
  void bsm_scal_(real *alpha, DalMatrix *A); /* useful for eg normalization */
  void bsm_scal_dia_(real *alpha, DalMatrix *A); /* used in response */
  void bsm_zerohalf_(const char *triangle, DalMatrix *A);
#if 0
  /* constructs new density (Dnew) from Fort                              *
   * The eigenvalue problem F * C = S * C * e must be reduced to          *
   * standard form Fort * C = C * e before calling this function          *
   * This can be done with the lapack routines pptrf, tptri and  trmm     *
   * Upper and lower bounds for the eigenvalues (lmin and lmax)           *
   * can be found using gerschgorin's theorem or lanczos for eigenvalues  *
   * Dnew should be initialized                                           */
  void bsm_d_from_f(DalMatrix *Fort, DalMatrix *Dnew, 
		    real* lmin, real* lmax, integer* nshells);// VKT
#endif
  
  void bsm_d_from_f_(DalMatrix *F, DalMatrix* S, integer* nshells,
                     DalMatrix *Dnew, integer *cycles);



 //Computes dE_SCF/dmu where mu is the damping in damped roothan Do
    //not implement.  Find dE/dmu = sum_nu,I dE/dC_nu,I dC_nu,I/dmu
    real bsm_dE_dmu(DalMatrix *Fnew, DalMatrix *Fdamp, DalMatrix *SDS,
                    DalMatrix *Cmo, integer *nocc);
    //Returns the sum of the elements Mat(from_row:to_row,ncol)
    //squared Do not implement.
    real bsm_column_norm(DalMatrix *mat, integer *ncol, integer *from_row,
			 integer*to_row);
    void bsm_precond_(const DalMatrix *M,  const DalMatrix *v,  DalMatrix *vp);
    //divide the occ-virt part of the matrix X_MO with the orbital
    //energy difference E_a - E_i. Do not implement.
    void bsm_mo_precond(integer *nocc,integer *symm, real *omega,
                        real *Eorb_final, real *X_MO);
    void bsm_ao_precond_(const integer *symm, const real * omega,
                         DalMatrix *FUP, DalMatrix *FUQ,
                         DalMatrix *DU, DalMatrix *X_AO);

    //Create or overwrite element in matrix
    void bsm_create_elm(integer *i, integer *j,real *val, DalMatrix *A);
    //scale the matrix A by the scalar alpha
    void bsm_scal(real *alpha, DalMatrix *A);
    void bsm_scal_dia(real *alpha, DalMatrix *A);
    void bsm_zero(DalMatrix *A);
    //set upper or lower triangle of A to zero
    //diagonal belongs to lower triangle
    void bsm_zerohalf(char *part, DalMatrix *A);
    void bsm_write_to_disk(integer *iunit, DalMatrix *A);             // VKT?
    void bsm_read_from_disk(integer *iunit, DalMatrix *A);            // VKT?
    //Change a vector to a matrix. The matrix is symmetric or antisymmetric.
    //PS: do not ask me...
    void bsm_VEC_TO_MAT(integer *symmetry, DalMatrix *VEC, DalMatrix *MAT);
    //Change a matrix to a vector. The matrix must be symmetric or
    //antisymmetric. However, be careful wether it is sym. or
    //antisym. when transforming back to a matrix//
    void BSM_TO_VEC(DalMatrix *MAT, DalMatrix *VEC);

    void bsm_get_sparsity_(DalMatrix* A, real *sparsity);

void bsm_read_from_unit(integer *iunit, DalMatrix* A,
                        void (*read_int) (void *data, integer *cnt, integer *where),
                        void (*read_real)(void *data, integer *cnt, real *where));
void bsm_write_to_unit(integer *iunit, DalMatrix* A,
                       void (*write_int) (void *unit, integer *cnt, integer *where),
                       void (*write_real)(void *unit, integer *cnt, real *where));


#ifdef __cplusplus
};
#endif
