/*


!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
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

!

*/
#define __CVERSION__
/* Written by Elias Rudberg, KTH, Stockholm */
#define _BSD_SOURCE 1

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <memory.h>
#include <time.h>
#include <stdarg.h>
#include "general.h"
#include "basisinfo.h"
#include "pi.h"


static int    global_outputLevel       = 1;

static void
do_output(char* s, int prio)
{
  if(prio > global_outputLevel)
    return;

  printf("%s\n", s);

#if defined USE_TEST_VERSION
  printf("%s\n", s);
#else
  fort_print(s);
  /*  printf("%s\n", s); */
#endif
}


static void 
do_output_2(int prio, const char* format, ...)
{
  char s[888];
  va_list a;

  if(prio > global_outputLevel)
    return;

  va_start(a, format);
  vsnprintf(s, sizeof(s), format, a);
  va_end(a);

#if defined USE_TEST_VERSION
  printf("%s\n", s);
#else
  fort_print(s);
  /*  printf("%s\n", s); */
#endif

  printf("%s\n", s);
}




/* FIXME: is it a way to make this routine shorter, cleaner? */
static int 
get_simple_primitives(
		      BasisFuncStruct* currBasisFunc,
		      DistributionSpecStruct* list,
		      int nInput,
		      int nListMax)
{
  int spd, contr, kk, n, j, nTerms, ii;
  real scaleFactor;

  /* make sure there is enough space left in list */
  if((nListMax - nInput) < MAX_NO_OF_PRIMITIVES_PER_BASIS_FUNC)
    {
      do_output_2(0, "error in get_simple_primitives: "
		  "not enough space left in list");
      return -1;
    }

  n = nInput;
  spd = currBasisFunc->shellType;
  contr = currBasisFunc->noOfContr;
  switch(spd)
    {
    case 0:
        /* 's' type shell */
      for(kk = 0; kk < contr; kk++)
	{
	  list[n].coeff = currBasisFunc->coeffList[kk];
	  list[n].exponent = currBasisFunc->exponentList[kk];
	  for(j = 0; j < 3; j++)
	    list[n].centerCoords[j] = currBasisFunc->centerCoords[j];
	  for(j = 0; j < 3; j++)
	    list[n].monomialInts[j] = 0;
	  n++;
	} /* END FOR kk */
      break;
    case 1:
        /* 'p' type shell */
      for(kk = 0; kk < contr; kk++)
	{
	  list[n].coeff = currBasisFunc->coeffList[kk];
	  list[n].exponent = currBasisFunc->exponentList[kk];
	  for(j = 0; j < 3; j++)
	    list[n].centerCoords[j] = currBasisFunc->centerCoords[j];
	  for(j = 0; j < 3; j++)
	    list[n].monomialInts[j] = 0;
	  switch(currBasisFunc->functionNumber)
	    {
	    case 0:
                /* function x */
	      list[n].monomialInts[0] = 1;
	      break;
	    case 1:
                /* function y */
	      list[n].monomialInts[1] = 1;
	      break;
	    case 2:
                /* function z */
	      list[n].monomialInts[2] = 1;
	      break;
	    default:
	      do_output("error: default reached", 0);
	      return -1;
	    } /* END SWITCH */
	  n++;
	} /* END FOR kk */
      break;
    case 2:
        /* 'd' type shell */
      for(kk = 0; kk < contr; kk++)
	{
	  switch(currBasisFunc->functionNumber)
	    {
	    case 0:
                /* function sqrt(3) * x * y */
	      nTerms = 1;
	      list[n].coeff = currBasisFunc->coeffList[kk];
	      /*list[n].coeff *= sqrt(3); */
	      list[n].monomialInts[0] = 1;
	      list[n].monomialInts[1] = 1;
	      list[n].monomialInts[2] = 0;
	      break;
	    case 1:
                /* function sqrt(3) * y * z */
	      nTerms = 1;
	      list[n].coeff = currBasisFunc->coeffList[kk];
	      /*list[n].coeff *= sqrt(3); */
	      list[n].monomialInts[0] = 0;
	      list[n].monomialInts[1] = 1;
	      list[n].monomialInts[2] = 1;
	      break;
	    case 2:
	      /* function (1/2) * (3*z*z - r*r) = (1/2) * (2*z*z - x*x - y*y) */
	      nTerms = 3;
	      /* term (1/2) * 2 * z * z */
	      list[n+0].coeff = currBasisFunc->coeffList[kk];
	      list[n+0].coeff *= 0.5 * 2.0 / sqrt(3);
	      list[n+0].monomialInts[0] = 0;
	      list[n+0].monomialInts[1] = 0;
	      list[n+0].monomialInts[2] = 2;
	      /* term -0.5 * x * x */
	      list[n+1].coeff = currBasisFunc->coeffList[kk];
	      list[n+1].coeff *= -0.5 / sqrt(3);
	      list[n+1].monomialInts[0] = 2;
	      list[n+1].monomialInts[1] = 0;
	      list[n+1].monomialInts[2] = 0;
	      /* term -0.5 * y * y */
	      list[n+2].coeff = currBasisFunc->coeffList[kk];
	      list[n+2].coeff *= -0.5 / sqrt(3);
	      list[n+2].monomialInts[0] = 0;
	      list[n+2].monomialInts[1] = 2;
	      list[n+2].monomialInts[2] = 0;
	      break;
	    case 3:
                /* function sqrt(3) * x * z */
	      nTerms = 1;
	      list[n].coeff = currBasisFunc->coeffList[kk];
	      /*list[n].coeff *= sqrt(3); */
	      list[n].monomialInts[0] = 1;
	      list[n].monomialInts[1] = 0;
	      list[n].monomialInts[2] = 1;
	      break;
	    case 4:
                /* function (1/2) * sqrt(3) * (x*x - y*y) */
	      nTerms = 2;
	      /* term 0.5 * sqrt(3) * x * x */
	      list[n+0].coeff = currBasisFunc->coeffList[kk];
	      list[n+0].coeff *= 0.5;/* * sqrt(3); */
	      list[n+0].monomialInts[0] = 2;
	      list[n+0].monomialInts[1] = 0;
	      list[n+0].monomialInts[2] = 0;
	      /* term -0.5 * sqrt(3) * y * y */
	      list[n+1].coeff = currBasisFunc->coeffList[kk];
	      list[n+1].coeff *= -0.5;/* * sqrt(3); */
	      list[n+1].monomialInts[0] = 0;
	      list[n+1].monomialInts[1] = 2;
	      list[n+1].monomialInts[2] = 0;
	      break;
	    default:
	      do_output_2(0, "error: default reached");
	      return -1;
	    } /* END SWITCH functionNumber */

	  for(ii = 0; ii < nTerms; ii++)
	    {
	      list[n].exponent = currBasisFunc->exponentList[kk];
	      for(j = 0; j < 3; j++)
		list[n].centerCoords[j] = currBasisFunc->centerCoords[j];
	      n++;
	    } /* END FOR ii */

	} /* END FOR kk */
      break;
    case 3:
        /* 'f' type shell */
      scaleFactor = 1.0/sqrt(15);
      for(kk = 0; kk < contr; kk++)
	{
	  switch(currBasisFunc->functionNumber)
	    {
	    case 0:
                /* function 0.5*sqrt(2.5)*(3*x*x-y*y)*y =  */
                /*    0.5*sqrt(2.5)*3*x*x*y - 0.5*sqrt(2.5)*y*y*y */
	      nTerms = 2;
	      /* term 0.5*sqrt(2.5)*3*x*x*y */
	      list[n+0].coeff = currBasisFunc->coeffList[kk];
	      list[n+0].coeff *= 0.5 * sqrt(2.5) * 3;
	      list[n+0].monomialInts[0] = 2;
	      list[n+0].monomialInts[1] = 1;
	      list[n+0].monomialInts[2] = 0;
	      /* term  - 0.5*sqrt(2.5)*y*y*y */
	      list[n+1].coeff = currBasisFunc->coeffList[kk];
	      list[n+1].coeff *= -0.5 * sqrt(2.5);
	      list[n+1].monomialInts[0] = 0;
	      list[n+1].monomialInts[1] = 3;
	      list[n+1].monomialInts[2] = 0;
	      break;
	    case 1:
                /* function sqrt(15)*x*y*z */
	      nTerms = 1;
	      list[n].coeff = currBasisFunc->coeffList[kk];
	      list[n].coeff *= sqrt(15);
	      list[n].monomialInts[0] = 1;
	      list[n].monomialInts[1] = 1;
	      list[n].monomialInts[2] = 1;
	      break;
	    case 2:
                /* function 0.5*sqrt(1.5)*(5*z*z - r*r)*y =  */
                /*    0.5*sqrt(1.5)*(4*z*z*y-x*x*y-y*y*y) */
	      nTerms = 3;
	      /* term 0.5*sqrt(1.5)*4*z*z*y */
	      list[n+0].coeff = currBasisFunc->coeffList[kk];
	      list[n+0].coeff *= 0.5*sqrt(1.5)*4.0;
	      list[n+0].monomialInts[0] = 0;
	      list[n+0].monomialInts[1] = 1;
	      list[n+0].monomialInts[2] = 2;
	      /* term -0.5*sqrt(1.5)*x*x*y */
	      list[n+1].coeff = currBasisFunc->coeffList[kk];
	      list[n+1].coeff *= -0.5*sqrt(1.5);
	      list[n+1].monomialInts[0] = 2;
	      list[n+1].monomialInts[1] = 1;
	      list[n+1].monomialInts[2] = 0;
	      /* term -0.5*sqrt(1.5)*y*y*y */
	      list[n+2].coeff = currBasisFunc->coeffList[kk];
	      list[n+2].coeff *= -0.5*sqrt(1.5);
	      list[n+2].monomialInts[0] = 0;
	      list[n+2].monomialInts[1] = 3;
	      list[n+2].monomialInts[2] = 0;
	      break;
	    case 3:
                /* function 0.5*(5*z*z-3*r*r)*z =  */
                /*    z*z*z - 1.5*x*x*z - 1.5*y*y*z */
	      nTerms = 3;
	      /* term z*z*z */
	      list[n+0].coeff = currBasisFunc->coeffList[kk];
	      list[n+0].coeff *= 1;
	      list[n+0].monomialInts[0] = 0;
	      list[n+0].monomialInts[1] = 0;
	      list[n+0].monomialInts[2] = 3;
	      /* term -1.5*x*x*z */
	      list[n+1].coeff = currBasisFunc->coeffList[kk];
	      list[n+1].coeff *= -1.5;
	      list[n+1].monomialInts[0] = 2;
	      list[n+1].monomialInts[1] = 0;
	      list[n+1].monomialInts[2] = 1;
	      /* term -1.5*y*y*z */
	      list[n+2].coeff = currBasisFunc->coeffList[kk];
	      list[n+2].coeff *= -1.5;
	      list[n+2].monomialInts[0] = 0;
	      list[n+2].monomialInts[1] = 2;
	      list[n+2].monomialInts[2] = 1;
	      break;
	    case 4:
                /* function 0.5*sqrt(1.5)*(5*z*z-r*r)*x =  */
                /*    0.5*sqrt(1.5)*(4*z*z*x - x*x*x - y*y*x) */
	      nTerms = 3;
	      /* term 0.5*sqrt(1.5)*4*z*z*x */
	      list[n+0].coeff = currBasisFunc->coeffList[kk];
	      list[n+0].coeff *= 0.5*sqrt(1.5)*4;
	      list[n+0].monomialInts[0] = 1;
	      list[n+0].monomialInts[1] = 0;
	      list[n+0].monomialInts[2] = 2;
	      /* term -0.5*sqrt(1.5)*x*x*x */
	      list[n+1].coeff = currBasisFunc->coeffList[kk];
	      list[n+1].coeff *= -0.5 * sqrt(1.5);
	      list[n+1].monomialInts[0] = 3;
	      list[n+1].monomialInts[1] = 0;
	      list[n+1].monomialInts[2] = 0;
	      /* term -0.5*sqrt(1.5)*y*y*x */
	      list[n+2].coeff = currBasisFunc->coeffList[kk];
	      list[n+2].coeff *= -0.5 * sqrt(1.5);
	      list[n+2].monomialInts[0] = 1;
	      list[n+2].monomialInts[1] = 2;
	      list[n+2].monomialInts[2] = 0;
	      break;
	    case 5:
                /* function 0.5*sqrt(15)*(x*x-y*y)*z =  */
                /*      0.5*sqrt(15)*(x*x*z - y*y*z) */
	      nTerms = 2;
	      /* term 0.5*sqrt(15)*x*x*z */
	      list[n+0].coeff = currBasisFunc->coeffList[kk];
	      list[n+0].coeff *= 0.5*sqrt(15);
	      list[n+0].monomialInts[0] = 2;
	      list[n+0].monomialInts[1] = 0;
	      list[n+0].monomialInts[2] = 1;
	      /* term 0.5*sqrt(15)*y*y*z */
	      list[n+1].coeff = currBasisFunc->coeffList[kk];
	      list[n+1].coeff *= -0.5*sqrt(15);
	      list[n+1].monomialInts[0] = 0;
	      list[n+1].monomialInts[1] = 2;
	      list[n+1].monomialInts[2] = 1;
	      break;
	    case 6:
                /* function 0.5*sqrt(2.5)*(x*x*x - 3*y*y*x) */
	      nTerms = 2;
	      /* term 0.5*sqrt(2.5)*x*x*x */
	      list[n+0].coeff = currBasisFunc->coeffList[kk];
	      list[n+0].coeff *= 0.5*sqrt(2.5);
	      list[n+0].monomialInts[0] = 3;
	      list[n+0].monomialInts[1] = 0;
	      list[n+0].monomialInts[2] = 0;
	      /* term -0.5*sqrt(2.5)*3*y*y*x */
	      list[n+1].coeff = currBasisFunc->coeffList[kk];
	      list[n+1].coeff *= -0.5*sqrt(2.5)*3;
	      list[n+1].monomialInts[0] = 1;
	      list[n+1].monomialInts[1] = 2;
	      list[n+1].monomialInts[2] = 0;
	      break;
	    default:
	      do_output_2(0, "error: default reached");
	      return -1;
	    } /* END SWITCH functionNumber */

	  for(ii = 0; ii < nTerms; ii++)
	    {
	      list[n].coeff *= scaleFactor;
	      list[n].exponent = currBasisFunc->exponentList[kk];
	      for(j = 0; j < 3; j++)
		list[n].centerCoords[j] = currBasisFunc->centerCoords[j];
	      n++;
	    } /* END FOR ii */
	} /* END FOR kk */
      break;
    default:
      do_output_2(0, "error in get_simple_primitives: "
		  "only spdf shells are currently implemented");
      return -1;
    } /* END SWITCH spd */

  if(n >= nListMax)
    {
      do_output_2(0, "error in get_simple_primitives: "
		  "(n >= nListMax)");
      return -1;
    }

  return n - nInput;
  
} /* END get_simple_primitives */

void getshellscnt_(int *);
void getshellno_(const int *no, int *contr, int *L, real *x, real *y, real *z,
                 const int *mxcontr, real *coefs, real *exps);

int 
get_shells(BasisInfoStruct* basisInfo)
{
  int i, nShells = 0;
  ShellSpecStruct* shellList;

  getshellscnt_(&nShells);
  if(nShells <= 0)
    {
      do_output_2(0, "error in getshellscnt");
      return -1;
    }
  basisInfo->noOfShells = nShells;
  basisInfo->shellList = 
    (ShellSpecStruct*)malloc(nShells * sizeof(ShellSpecStruct));
  shellList = basisInfo->shellList;
  for(i = 0; i < nShells; i++)
    {
      const int MaxContr = MAX_NO_OF_CONTR_GAUSSIANS;
      int fShell = i+1; /* fortran shell number */
      int contr, spd, kk;
      real x, y, z;
      real coeffList[MAX_NO_OF_CONTR_GAUSSIANS];
      real exponentList[MAX_NO_OF_CONTR_GAUSSIANS];
      getshellno_(&fShell, &contr, &spd, &x, &y, &z,
                  &MaxContr, coeffList, exponentList);
      if(contr < 1)
	{
	  do_output_2(0, "error reading shell number %i %i", i+1, contr);
	  return -1;
	}
      if(contr > MAX_NO_OF_CONTR_GAUSSIANS)
	{
	  do_output_2(0, "error: too many contracted gaussians, %i > %i", 
		      contr, (int)MAX_NO_OF_CONTR_GAUSSIANS);
	  do_output_2(0, "change constant MAX_NO_OF_CONTR_GAUSSIANS"
		      " in file grid-gen2.c");
	  return -1;
	}

      for(kk = 0; kk < contr; kk++)
	{
	  shellList[i].coeffList[kk] = coeffList[kk];
	  shellList[i].exponentList[kk] = exponentList[kk];
	}
      shellList[i].centerCoords[0] = x;
      shellList[i].centerCoords[1] = y;
      shellList[i].centerCoords[2] = z;
      shellList[i].extent = 0;
      shellList[i].shellType = spd;
      shellList[i].noOfContr = contr;
    } /* END FOR i */
  
  return 0;
}

int 
get_basis_funcs(BasisInfoStruct* basisInfo)
{
  /* create list of 'basis functions',  */
  /* and set startIndexInMatrix for each shell */
  int nShells = basisInfo->noOfShells;
  int count = 0;
  int i, j, kk, nFunctions;
  ShellSpecStruct* currShell;
  BasisFuncStruct* basisFuncList;

  for(i = 0; i < nShells; i++)
    {
      currShell = &basisInfo->shellList[i];
      currShell->startIndexInMatrix = count;
      nFunctions = 1 + currShell->shellType * 2;
      count += nFunctions;
      currShell->noOfBasisFuncs = nFunctions;
    }
  basisInfo->noOfBasisFuncs = count;
  basisInfo->basisFuncList = 
    (BasisFuncStruct*)malloc(count * sizeof(BasisFuncStruct));
  basisFuncList = basisInfo->basisFuncList;
  count = 0;
  for(i = 0; i < nShells; i++)
    {
      currShell = &basisInfo->shellList[i];
      nFunctions = currShell->noOfBasisFuncs;
      for(j = 0; j < nFunctions; j++)
	{
	  basisFuncList[count].noOfContr = currShell->noOfContr;
	  for(kk = 0; kk < currShell->noOfContr; kk++)
	    {
	      basisFuncList[count].coeffList[kk] = currShell->coeffList[kk];
	      basisFuncList[count].exponentList[kk] = 
		currShell->exponentList[kk];
	    } /* END FOR kk */
	  for(kk = 0; kk < 3; kk++)
	    basisFuncList[count].centerCoords[kk] = 
	      currShell->centerCoords[kk];
	  basisFuncList[count].extent = currShell->extent;
	  basisFuncList[count].shellType = currShell->shellType;
	  basisFuncList[count].functionNumber = j;
	  count++;
	} /* END FOR j */
    } /* END FOR i each shell */
  if(count != basisInfo->noOfBasisFuncs)
    {
      do_output_2(0, "error in get_basis_funcs: "
		  "(count != basisInfo->noOfBasisFuncs)");
      return -1;
    }
  return 0;
}


int
get_simple_primitives_all(BasisInfoStruct* basisInfo)
{
  int nbast = basisInfo->noOfBasisFuncs;
  int maxNoOfSimplePrimsTot = nbast * MAX_NO_OF_PRIMITIVES_PER_BASIS_FUNC;
  DistributionSpecStruct* list = 
    malloc(maxNoOfSimplePrimsTot * sizeof(DistributionSpecStruct));
  
  BasisFuncStruct* basisFuncList = basisInfo->basisFuncList;

  /* create list of 'simple primitives' */
  int n = 0;
  int i, nBytes;
  for(i = 0; i < nbast; i++)
    {
      BasisFuncStruct* currBasisFunc = &basisFuncList[i];
      int noOfPrimitives = get_simple_primitives(currBasisFunc,
						 list,
						 n,
						 maxNoOfSimplePrimsTot);
      if(noOfPrimitives <= 0)
	{
	  do_output("error in get_simple_primitives", 0);
	  return -1;
	}
      currBasisFunc->noOfSimplePrimitives = noOfPrimitives;
      currBasisFunc->simplePrimitiveIndex = n;
      n += noOfPrimitives;
    } /* END FOR i */
  nBytes = n * sizeof(DistributionSpecStruct);
  basisInfo->simplePrimitiveList = malloc(nBytes);
  memcpy(basisInfo->simplePrimitiveList, list, nBytes);
  free(list);
  basisInfo->noOfSimplePrimitives = n;
  return 0;
}


#define K_MAX_DIM 44

typedef struct{
  real a0;
  real a1;
} polydeg1struct;

static int 
multiply_polynomials(real result[], 
		     polydeg1struct* polydeg1, 
		     int dim, 
		     real a[])
{
  int i;
  real p1[K_MAX_DIM + 1];
  real p2[K_MAX_DIM + 1];
  if(dim >= (K_MAX_DIM-1))
    return -1;
  for(i = 0; i <= dim; i++)
    p1[i] = a[i]*polydeg1->a0;
  p1[dim+1] = 0;
  p2[0] = 0;
  for(i = 0; i <= dim; i++)
    p2[i+1] = a[i]*polydeg1->a1;
  for(i = 0; i <= (dim+1); i++)
    result[i] = p1[i] + p2[i];
  return 0;
} /* END multiply_polynomials */



int
get_product_simple_prims(DistributionSpecStruct* primA,
			 DistributionSpecStruct* primB,
			 DistributionSpecStruct resultList[],
			 int maxCount)
{
  /* use the Gaussian product rule */
  real sum = 0;
  real newCenter[3];
  real CxCyCz, AiAj, alphaNew;
  int k, l, m, nn;
  real poly0[K_MAX_DIM];
  real poly1[K_MAX_DIM];
  real poly2[K_MAX_DIM];
  real tempPoly[K_MAX_DIM];
  real tempPoly2[K_MAX_DIM];
  real tempPoly3[K_MAX_DIM];
  int tempPolyDegree, tempPoly2Degree;
  int poly0degree, poly1degree, poly2degree;
  polydeg1struct polyDeg1;
  real* poly;
  int* degreePtr;

  for(k = 0; k < 3; k++)
    {
      real temp = primA->centerCoords[k] - primB->centerCoords[k];
      sum += temp * temp;
    } /* END FOR k */
  CxCyCz = exp(-primA->exponent * primB->exponent * 
		    sum / (primA->exponent + primB->exponent));
  AiAj = primA->coeff * primB->coeff;
  alphaNew = primA->exponent + primB->exponent;
  for(k = 0; k < 3; k++)
    {
      newCenter[k] = 
	(primA->exponent * primA->centerCoords[k] +
	 primB->exponent * primB->centerCoords[k]) /
	(primA->exponent + primB->exponent);
    } /* END FOR k */

  /* do product of polynomials */
  /* one coordinate at a time */
  for(k = 0; k < 3; k++)
    {
      switch(k)
	{
	case 0: poly = poly0; degreePtr = &poly0degree; break;
	case 1: poly = poly1; degreePtr = &poly1degree; break;
	case 2: poly = poly2; degreePtr = &poly2degree; break;
	default: return -1;
	} /* END SWITCH k */
      tempPoly[0] = 1;
      tempPolyDegree = 0;
      for(m = 0; m < primA->monomialInts[k]; m++)
	{
	  polyDeg1.a0 = -primA->centerCoords[k];
	  polyDeg1.a1 = 1;
	  if(multiply_polynomials(tempPoly2, &polyDeg1, 
				  tempPolyDegree, tempPoly) != 0)
	    return -1;
	  tempPolyDegree++;
	  memcpy(tempPoly, 
		 tempPoly2, 
		 (tempPolyDegree+1)*sizeof(real));
	} /* END FOR m */
      for(m = 0; m < primB->monomialInts[k]; m++)
	{
	  polyDeg1.a0 = -primB->centerCoords[k];
	  polyDeg1.a1 = 1;
	  if(multiply_polynomials(tempPoly2, &polyDeg1, 
				  tempPolyDegree, tempPoly) != 0)
	    return -1;
	  tempPolyDegree++;
	  memcpy(tempPoly, 
		 tempPoly2, 
		 (tempPolyDegree+1)*sizeof(real));
	} /* END FOR m */

      /* now do variable change */
      for(m = 0; m < K_MAX_DIM; m++)
	poly[m] = 0;
      tempPoly2Degree = 0;
      for(m = 0; m <= tempPolyDegree; m++)
	{
	  tempPoly2[0] = tempPoly[m];
	  tempPoly2Degree = 0;
	  for(l = 0; l < m; l++)
	    {
	      polyDeg1.a0 = newCenter[k];
	      polyDeg1.a1 = 1;
	      if(multiply_polynomials(tempPoly3, 
				      &polyDeg1, 
				      tempPoly2Degree, 
				      tempPoly2) != 0)
		return -1;
	      tempPoly2Degree++;
	      memcpy(tempPoly2, 
		     tempPoly3, 
		     (tempPoly2Degree+1)*sizeof(real));
	    } /* END FOR l */
	  for(l = 0; l <= tempPoly2Degree; l++)
	    {
	      poly[l] += tempPoly2[l];
	    } /* END FOR l */
	} /* END FOR m */
      *degreePtr = tempPoly2Degree;
    } /* END FOR k */

  nn = 0;
  for(k = 0; k <= poly0degree; k++)
    {
      int l;
      for(l = 0; l <= poly1degree; l++)
	{
	  int m;
	  for(m = 0; m <= poly2degree; m++)
	    {
	      real newCoeff = AiAj * CxCyCz * poly0[k] * poly1[l] * poly2[m];

	      real sqrtValue = sqrt(pi / alphaNew);
	      real absvalue = newCoeff * sqrtValue * sqrtValue * sqrtValue;
	      if(absvalue < 0) absvalue *= -1;

	      /* add one function to final list */
	      resultList[nn].coeff = newCoeff;
	      resultList[nn].exponent = alphaNew;

	      memcpy(resultList[nn].centerCoords, 
		     newCenter, 
		     3 * sizeof(real));
	      resultList[nn].monomialInts[0] = k;
	      resultList[nn].monomialInts[1] = l;
	      resultList[nn].monomialInts[2] = m;
	      
	      nn++;
	      if(nn >= maxCount)
		{
		  do_output_2(0, "error in read_density_file: "
			      "maxCount exceeded");
		  do_output_2(0, "nn = %i, maxCount = %i", 
			      nn, maxCount);
		  return -1;
		}
	    } /* END FOR m */
	} /* END FOR l */
    } /* END FOR k */

  return nn;
}


int
get_product_simple_primitives(BasisInfoStruct* basisInfoA, int iA,
			      BasisInfoStruct* basisInfoB, int iB,
			      DistributionSpecStruct resultList[],
			      int maxCount)
{
    /* printf("entering get_product_simple_primitives\n"); */

  BasisFuncStruct* basisFuncA = &basisInfoA->basisFuncList[iA];
  int nPrimsA = basisFuncA->noOfSimplePrimitives;
  int Astart = basisFuncA->simplePrimitiveIndex;
  BasisFuncStruct* basisFuncB = &basisInfoB->basisFuncList[iB];
  int nPrimsB = basisFuncB->noOfSimplePrimitives;
  int Bstart = basisFuncB->simplePrimitiveIndex;
  int n = 0;
  int i;
  /* printf("nPrimsA = %i, nPrimsB = %i\n", nPrimsA, nPrimsB); */
  /* printf("Astart = %i, Bstart = %i\n", Astart, Bstart);     */
  for(i = 0; i < nPrimsA; i++)
    {
      DistributionSpecStruct* primA = 
	&basisInfoA->simplePrimitiveList[Astart + i];
      int j;
      for(j = 0; j < nPrimsB; j++)
	{
	  DistributionSpecStruct* primB = 
	    &basisInfoB->simplePrimitiveList[Bstart + j];
	  int nNewPrims = get_product_simple_prims(primA, 
						   primB, 
						   &resultList[n],
						   maxCount - n);
	  if(nNewPrims <= 0)
	    {
	      do_output_2(0, "error in get_product_simple_prims");
	      return -1;
	    }
	  n += nNewPrims;
	}
    }
  return n;
}


 
