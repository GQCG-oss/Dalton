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
#include <errno.h>
  
#define INT long long
  
int linux_mem_allo__(INT * waddr,int * nbytes,double *wrk)
  {
   int test=0;
   int success=0;
   int indx;
   int indy;
   int indz;
   double *where;

  where=(double *)malloc(*nbytes);
 
  if(where==NULL)
    {
      (void)fprintf(stderr,"Memory allocation failure");
      if( (errno>0) && (errno<sys_nerr) )
      {
        (void)fprintf(stderr," (%s) ",sys_errlist[errno]);
      }
      (void)fprintf(stderr," \n");
      return(0);
    }

   if(test) 
      fprintf(stderr,"Allocating block at %i, work: %i, offset: %ld dwords\n",
            (int)where, (int)wrk, (long)*waddr);
   indx=(int)((void *)where);
   indy=(int)((void *)wrk);
   indz=(indx-indy)/(sizeof(double)/sizeof(char))+1;
   *waddr = indz;
  
  /*------> Finished <-------*/
  return 1;
  }
