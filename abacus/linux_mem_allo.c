/******************************************************
 * linux_mem_allo: dynamically allocates memory       *
 *                 for Dalton calculation for         *
 *                 linux systems.                     *
 *                                                    *
 * Input:          waddr  - Pointer to memory address *
 *                 nbytes - Number of bytes to be     *
 *                          allocated                 *
 *                 wrk    - First element in WORK     *
 ******************************************************/
#include <stdio.h>
#include <stdlib.h>
  
#define INT long long
  
int linux_mem_allo__(INT * waddr,int * nbytes,double *wrk)
  {
   int test=0;
   double *where;
  
   /* return 0, fortran code will print error message */
   if( (where=(double *)malloc(*nbytes)) == NULL) 
      return 0;
   
   if(test) 
      fprintf(stderr,"Allocating block at %i, work: %i, offset: %ld dwords\n",
            (int)where, (int)wrk, (long)*waddr);

   *waddr = ((INT)where-(INT)wrk)/sizeof(double) +1;
  
  /*------> Finished <-------*/
  return 1;
  }
