/*
       Copyright (c) 2001 by the authors of Dalton (see below).
       All Rights Reserved.
    
       The source code in this file is part of
       "Dalton, a molecular electronic structure program, Release 1.2
       (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
       J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
       O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
       T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
       H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
       H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
       T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
       B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
       and O. Vahtras"
    
       This source code is provided under a written licence and may be
       used, copied, transmitted, or stored only in accord with that
       written licence.
    
       In particular, no part of the source code or compiled modules may
       be distributed outside the research group of the licence holder.
       This means also that persons (e.g. post-docs) leaving the research
       group of the licence holder may not take any part of Dalton,
       including modified files, with him/her, unless that person has
       obtained his/her own licence.
    
       For questions concerning this copyright write to:
          dalton-admin@kjemi.uio.no
    
       For information on how to get a licence see:
          http://www.kjemi.uio.no/software/dalton/dalton.html

 ******************************************************
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
  
typedef void (*DaltonDriver)(double* work, int* lmwork, double* wrkdlm,
                             int* master, int* mynum);
extern DaltonDriver nodstr_;
extern DaltonDriver exedrv_;
void cexe_(DaltonDriver drv, int* nwords, double* wrkdlm, 
           int* master, int* mynum)
{
    double* mem_block;
    static const int debug = 0;

    mem_block = (double*)malloc((*nwords+2)*sizeof(double));
    if(!mem_block) {
        fprintf(stderr,"CEXE: Cannot allocate  work: %i dwords\n", *nwords);
        dalton_quit("CEXE: Cannot allocate work: %i dwords\n", *nwords);
    }
    mem_block[0] = *wrkdlm;
    mem_block[1+*nwords] = *wrkdlm;
    if(debug) fprintf(stderr,"CEXE: Allocating  work: %i dwords\n", *nwords);
    drv(mem_block, nwords, wrkdlm,master,mynum);
    if(debug) fprintf(stderr,"CEXE finished.\n");
}
/* ioff=iallor8(work,nwords)  work(ioff) is the first position */
int  iallor8_(double *work,int *nwords){
  double *iadd=(double *)calloc((*nwords),8);
  if(iadd==NULL){fprintf(stderr,"iallor8: cannot allocate\n");exit(2);}
  return (iadd-work +1);
}
int  ialloi4_(int *work,int *nwords){
  int *iadd=(int *)calloc((*nwords),4);
  if(iadd==NULL){fprintf(stderr,"ialloi4: cannot allocate\n");exit(2);}
  return (iadd-work +1);
}
int  ialloi2_(short int *work,int *nwords){
  short int *iadd=(short int *)calloc((*nwords),2);
  if(iadd==NULL){fprintf(stderr,"ialloi2: cannot allocate\n");exit(2);}
  return (iadd-work +1);
}
int  ialloi1_(char *work,int *nwords){
  char *iadd=(char *)calloc((*nwords),1);
  if(iadd==NULL){fprintf(stderr,"ialloi1: cannot allocate\n");exit(2);}
  return (iadd-work +1);
}
void memrelr8_(double *work){  /*call memrel8(work(ioff))  */
free(work);
}
void memreli4_(int *work){
free(work);
}
void memreli2_(short int *work){
free(work);
}
void memreli1_(char *work){
free(work);
}
