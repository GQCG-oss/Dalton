#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>

/*
This file contains auxiliary functions for the sparse matrix
functionality implemented in matop_csr.f90.

Rasmus Andersen <raand@chem.au.dk>, June 2010
*/

#if defined(VAR_INT64)
#include <stdint.h>
   typedef long integer;
   /* typedef long int; */
#else
   typedef int integer;
#endif


/*
Print a Compressed Sparse Row matrix in a standard matrix fashion
*/
integer mat_csr_pretty_print_(double *val, integer *col, integer *row, integer *nrow){
  integer i,j,k,cur_col, cur_val;
  cur_val=0;
  /*printf("will loop nrows, from %i to %i\n", 0, *nrow);*/
  for (i=0; i<*nrow; i++){
    k=0;
    cur_col=0;
    printf("\n");
    /*printf("will loop row, from %i to %i\n",row[i]-1 ,row[i+1]-1);*/
    for (j=row[i]-1; j<row[i+1]-1; j++){
      k = col[j]-1;
      if (k == cur_col){
	printf("%1.3f\t", val[cur_val]);
	cur_val++;
      }
      else{
	printf("  -   \t");
	j--;
      }
      cur_col++;
    }
    /*printf("will loop remaining, from %i to %i\n",row[i+1]-1, *nrow );*/
    for (j=row[i+1]-1; j<*nrow; j++){
      printf("  -   \t");
    }
  }
  printf("\n");
  return 0;
}

/*
Print a Compressed Sparse Row matrix in sets of 4 columns
*/

/*
integer mat_csr_column_print_(integer *fd, double *val, integer *col, integer *row, integer *nrow){
  integer i,j,k,cur_col, cur_val;
  
  printf("fd is %i\n", *fd);
  write(*fd, (const void *)fd, 4);
  return 0;
}
*/

/*
Clean a CSR matrix, i.e. remove all elements in its values array that 
are below a given threshold. Adjust col and row arrays accordingly.
*/
integer mat_csr_cleanup_(double *val, integer *col, integer *row, integer *dim, integer *nnz, double *tol){
  /*printf("In dcsr_cleanup, nnz is %i!!\n", *nnz);*/
  /*printf("Last element in c is : %1.14f\n", val[*nnz-1]);*/
  double *d1, *d2;
  integer *c1, *c2, *r1, *r2;
  integer nc, ncnew, rr, ic, ir, n;
  size_t realloc_size;
  double *new_val;
  integer *new_col;

  if (*nnz == 0){
    return 0;
  }

  n = *dim;
  r1 = row;
  r2 = r1+1;
  rr = *r1;
  d1 = d2 = val;
  c1 = c2 = col;
  for (ir=0; ir<n; ir++) {
    nc = ncnew = *r2 - rr;
    for (ic=0; ic<nc; ic++) {
      if (fabs(*d2) < *tol) {
	d2++;
	c2++;
	ncnew--;
      } 
      else {
	if (d1 != d2) {
	  *d1 = *d2;
	  *c1 = *c2;
	}
	/* printf("before: %e, addr: %p\n",*d1,d1); */
	d1++;
	/* printf("after: %e, addr: %p\n",*d1,d1); */
	d2++;
	c1++;
	c2++;
      }
    }
    rr = *r2;
    *r2 = *r1 + ncnew;
    r1++;
    r2++;
  }
  *nnz = row[*dim]-1;
  /* resize the allocated memory to the new size after 
     removal of near-zeroes. If nnz is zero, just free 
     val and col arrays.
     realloc(3): if size is equal to zero, and ptr is not 
     NULL, then the  call is equivalent  to  free(ptr)
  */
  realloc_size = *nnz * sizeof(double);
  new_val = (double *) realloc(val, realloc_size);
  realloc_size = *nnz * sizeof(integer);
  new_col = (integer *) realloc(col, realloc_size);
  if (*nnz == 0){
    val = new_val;
    col = new_col;
    return 0;
  }
  assert((new_val != NULL) && (new_col != NULL));
  /*realloc(3): If the area pointed to was moved, a free(ptr) is done.*/
  val = new_val;
  col = new_col;  
  return 0;
}

