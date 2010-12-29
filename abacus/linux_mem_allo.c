/*
!
!...   Copyright (c) 2010 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2010), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!

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

#if defined(VAR_INT64)
#include <stdint.h>
typedef int64_t integer;
#else
typedef int integer;
#endif

void dalton_quit(const char* format, ...);
typedef void (*DaltonDriver)(double* work, integer* lmwork, double* wrkdlm);
extern DaltonDriver nodedriver_;
extern DaltonDriver exedrv_;
void cexe_(DaltonDriver drv, integer* nwords, double* wrkdlm)
{
    double* mem_block;
    static const int debug = 0;

    mem_block = (double*)malloc((*nwords+2)*sizeof(double));
    if(!mem_block) {
        fprintf(stderr,"CEXE: Cannot allocate  work: %ld dwords\n", (long)*nwords);
        dalton_quit("CEXE: Cannot allocate work: %ld dwords\n", (long)*nwords);
    }
    mem_block[0] = *wrkdlm;
    mem_block[1+*nwords] = *wrkdlm;
    if(debug) fprintf(stderr,"CEXE: Allocating  work: %ld dwords\n", (long)*nwords);
    drv(mem_block, nwords, wrkdlm);
    if(debug) fprintf(stderr,"CEXE finished.\n");
}
#ifdef VAR_G77
/* ioff=iallor8(work,nwords)  work(ioff) is the first position */
integer  iallor8_(double *work,integer *nwords){
  double *iadd=(double *)calloc((*nwords),8);
  if(iadd==NULL){fprintf(stderr,"iallor8: cannot allocate\n");exit(2);}
  return (iadd-work +1);
}
integer ialloi4_(int *work,integer *nwords){
  int *iadd=(int *)calloc((*nwords),4);
  if(iadd==NULL){fprintf(stderr,"ialloi4: cannot allocate\n");exit(2);}
  return (iadd-work +1);
}
integer ialloi2_(short int *work,integer *nwords){
  short int *iadd=(short int *)calloc((*nwords),2);
  if(iadd==NULL){fprintf(stderr,"ialloi2: cannot allocate\n");exit(2);}
  return (iadd-work +1);
}
integer  ialloi1_(char *work,integer *nwords){
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
#endif /* #ifdef VAR_G77 */
