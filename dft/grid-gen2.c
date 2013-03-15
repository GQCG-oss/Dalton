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
/* The multi-scale cartesian cubature grid generator. The original
 * implementation described in M. Challacombe, JCP 113(22),
 * p.10037. This one is modified wrt to the paper.
 *
 *  Elias Rudberg, 2004-04
*/

#define __CVERSION__
#if !defined(SYS_DEC)
/* XOPEN compliance is missing on old Tru64 4.0E Alphas and pow() prototype
 * is not specified. */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1
#endif

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <memory.h>
#include <time.h>
#include <stdarg.h>
#include "general.h"
#include "grid-gen.h"
#include "basisinfo.h"

/* #define USE_PTHREADS */
#ifdef USE_PTHREADS
#include <pthread.h>
#endif

static const int CUBATURE_RULE = 3;
static const int CUBATURE_RULE_2 = 6;
static const int NO_OF_DIMENSIONS = 3;  /* 1, 2 or 3 */
static const real CONSTPI = M_PI;


/* Whether hierarchical cubature should be used or not. Modified by
 * a keyword in the input file. */
int global_useHiCu = 0;

static int    global_maxNoOfShells     = 44444;

/* This variable is used to keep track of when it is time to 
   create a new grid, together with parameter global_nFreeze. */
static int    global_gridCount = 0;


/* ---------- CONFIGURABLE PARAMETERS ------------------------------
   ---------- CONFIGURABLE PARAMETERS ------------------------------
   ---------- CONFIGURABLE PARAMETERS ------------------------------ */
/* number of threads */
static int    global_nThreads          = 1;

/* flag for test integration. If turned on, the grid file is
reopened after it has been created, and the density is integrated
using the newly created grid. Just to check that it gives the 
same result as reported by the grid generation and by 
the dft integrator. */
static int    global_doTestIntegration = 0; 

/* Output level. 0 means minimum output, 1 means little output,
   2 means a lot of output. */
static int    global_outputLevel       = 1;

/* Number of iterations to use the same grid, before creating
a new one. nFreeze=1 gives a new grid for each iteration.
nFreeze=1000 means that the first grid is used throughout
the whole calculation. */
static int    global_nFreeze           = 1000;

/* Threshold value for distributions. A gaussian is ignored in areas
where its value is below this threshold. A low value is 
computationally expensive. */
static real   global_targetRhoError = 1.0e-10;

/* Cutoff value used to decide which gaussian products should
be ignored. A product with a lower coefficient than this 
will be thrown away. */
static real   global_distrCutoff    = 1.0e-11;

/* 3d space is partitioned into boxes of this size. */
static real   global_boxdist        = 1.5;

/* Main threshold error value for grid generation. The difference of
analytical and numerical integrals is forced below this value.
This is the most important parameter, and probably the only one
that a typical user should worry about. */
static real   global_maxerror       = 1.0e-7;






#define USE_EXP_STD
#define USE_ERF_STD
#define DO_EXTRA_ERROR_CHECKING
#define FILE_BATCH_N 200000
#define MAX_NO_OF_POINTS_PER_BATCH 100
#define MAX_NO_OF_SHLBLOCKS 44444
#define EXPONENT_DIFF_LIMIT 1e-22
#define DISTR_CENTER_DIST_LIMIT 1e-22
#define N_BATCH_JOBS 22
#define MAX_NO_OF_POINTS_PER_WRITE 50000

#if !defined(INSTALL_WRKMEM)
#define USE_TEST_VERSION
#endif

/*#define USE_TEST_VERSION */

/*#define SKIP_THREADS */


#ifdef USE_PTHREADS
pthread_mutex_t globalOutputMutex = PTHREAD_MUTEX_INITIALIZER;
#endif


/*//////////////////////////////////////////////////////////////////////// */
/*/////////////////  typedef section  //////////////////////////////////// */
/*//////////////////////////////////////////////////////////////////////// */

#if 0
struct DistributionSpecContr_{
  int noOfContr;
  real coeffList[MAX_NO_OF_CONTR_GAUSSIANS];
  real exponentList[MAX_NO_OF_CONTR_GAUSSIANS];
  real centerCoords[3]; /* x0, y0, z0 */
  int monomialInts[3];  /* nx, ny, nz */
};
typedef struct DistributionSpecContr_ DistributionSpecContr;
#endif

typedef struct
{
  int noOfShells;
  ShellSpecStruct* shellList;
  int nbast;
  const real* dmat;
  BasisFuncStruct* basisFuncList;
  int noOfDistributions;
  DistributionSpecStruct* distrList;
} DensitySpecStruct;



struct BoxStruct_{
    real min[3]; /* xmin, ymin, zmin */
    real max[3]; /* xmax, ymax, zmax */
};
typedef struct BoxStruct_ BoxStruct;

struct rhoTreeNode_{
  BoxStruct box;
    struct rhoTreeNode_* child1; /* NULL for leaf node */
    struct rhoTreeNode_* child2; /* NULL for leaf node */
    int distrIndex;      /* -1 for non-leaf node */
};
typedef struct rhoTreeNode_ rhoTreeNode;

typedef struct
{
  DensitySpecStruct density;
  int noOfNonzeroDistributions;
  int* nonZeroDistrIndexList;
  int noOfNonzeroShells;
  int* nonZeroShellsIndexList;
  real maxerrorPerBox;
} compute_grid_for_box_params_struct;

  typedef struct
  {
    DensitySpecStruct* density;
    rhoTreeNode* rhoTreeRootNode;
    rhoTreeNode* rhoTreeRootNodeShells;
    real maxerror;
    FILE* gridFile;
    BoxStruct* startBox;
    int Nx;
    int Ny;
    int Nz;
#ifdef USE_PTHREADS
    pthread_mutex_t* fileMutex;
    pthread_mutex_t* jobMutex;
    pthread_t thread;
#endif
    int* currJobNumber;
      int noOfPoints;         /* OUTPUT */
      int noOfWrittenBatches; /* OUTPUT */
      real integralResult;    /* OUTPUT */
    int threadNo;
  } compute_grid_thread_func_struct;


/*//////////////////////////////////////////////////////////////////////// */
/*/////////////////  end of typedef section  ///////////////////////////// */
/*//////////////////////////////////////////////////////////////////////// */


/* Solid harmonics based on the table 6.3 of Molecular
 * Electronic-Structure Theory by Helgaker, Jørgensen and Olsen. */

#define solid_harmonic_s_0(x, y, z, x2, y2, z2, r2) 1

#define solid_harmonic_p_0(x, y, z, x2, y2, z2, r2) x
#define solid_harmonic_p_1(x, y, z, x2, y2, z2, r2) y
#define solid_harmonic_p_2(x, y, z, x2, y2, z2, r2) z

#define solid_harmonic_d_0(x, y, z, x2, y2, z2, r2) (x * y)
#define solid_harmonic_d_1(x, y, z, x2, y2, z2, r2) (y * z)
#define solid_harmonic_d_2(x, y, z, x2, y2, z2, r2) ((2 * z2 - x2 - y2) / (2 * sqrt(3)))
#define solid_harmonic_d_3(x, y, z, x2, y2, z2, r2) (x * z)
#define solid_harmonic_d_4(x, y, z, x2, y2, z2, r2) (0.5 * (x2 - y2))

#define solid_harmonic_f_0(x, y, z, x2, y2, z2, r2) ((0.5 * sqrt(2.5) * (3 * x2 - y2) * y) / sqrt(15))
#define solid_harmonic_f_1(x, y, z, x2, y2, z2, r2) (x * y * z)
#define solid_harmonic_f_2(x, y, z, x2, y2, z2, r2) (0.5 * sqrt(1.5) * (5 * z2 - r2) * y / sqrt(15))
#define solid_harmonic_f_3(x, y, z, x2, y2, z2, r2) (0.5 * (5 * z2 - 3 * r2) * z / sqrt(15))
#define solid_harmonic_f_4(x, y, z, x2, y2, z2, r2) (0.5 * sqrt(1.5) * (5 * z2 - r2) * x / sqrt(15))
#define solid_harmonic_f_5(x, y, z, x2, y2, z2, r2) (0.5 * (x2 - y2) * z)
#define solid_harmonic_f_6(x, y, z, x2, y2, z2, r2) (0.5 * sqrt(2.5) * (x2 - 3 * y2) * x / sqrt(15))



void 
do_output_2(int prio, const char* format, ...)
{
  char s[888];
  va_list a;

  if(prio > global_outputLevel)
    return;

#ifdef USE_PTHREADS
  pthread_mutex_lock(&globalOutputMutex);
#endif

  va_start(a, format);
  vsnprintf(s, sizeof(s), format, a);
  va_end(a);

#if defined USE_TEST_VERSION
  printf("USE_TEST_VERSION defined\n");
  printf("%s\n", s);
#else
  fort_print(s);
#endif

#ifdef USE_PTHREADS
  pthread_mutex_unlock(&globalOutputMutex);
#endif
}



static void
do_error_exit(const char* s)
{
  char ss[888];
  sprintf(ss, "error_exit: %s\n", s);
  printf("%s\n",ss);
  do_output_2(0, ss);
  fprintf(stderr,"%s\n", ss);
  exit(1);
}

void*
dal_malloc_safe_(size_t sz, const char *place, int line)
{
  void* res = malloc(sz);
  if(!res) {
    do_output_2(0, "error in dal_malloc_safe, '%s', sz = %i, line %i\n", 
		place, sz, line);
    do_error_exit("error in dal_malloc_safe_");
  }
  return res;
}

void
dal_free(void *a)
{
  free(a);
}

#define dal_malloc_safe(sz) dal_malloc_safe_((sz),__FUNCTION__, __LINE__)





#if 0
void
dftgridparams_(char *line, int len)
{
    int i = 0, tokenlen;
    char *sp;
    char* endPtr;
    printf("dftgridparams_ len = %i\n", len);
    printf("parsing '%s'\n", line);
    endPtr = line + len;
    for(i=0; i<len; i+=tokenlen) {
        sp = strchr(line+i, ' ');
        /*printf("parsing '%s'\n", line+i); */
        if(sp)
            tokenlen = sp-(line+i)+1;
        else tokenlen=100000;
        /*printf("tokenlen = %i\n", tokenlen); */
        if(strncmp(line+i, "box=", 4) ==0) {
            global_boxdist = atof(line+i+4);
            printf("New box distance %f\n", global_boxdist);
            continue;
        }
        if(strncmp(line+i, "maxerror=", 9) ==0) {
            global_maxerror = atof(line+i+9);
            printf("New maxerror %g\n", global_maxerror);
            continue;
        }
        if(strncmp(line+i, "nthreads=", 9) ==0) {
            global_nThreads = atoi(line+i+9);
            printf("New global_nThreads %i\n", global_nThreads);
            continue;
        }
        if(strncmp(line+i, "dotest=", 7) ==0) {
            global_doTestIntegration = atoi(line+i+7);
            printf("New global_doTestIntegration %i\n", 
                   global_doTestIntegration);
            continue;
        }
        if(strncmp(line+i, "usehicu=", 8) ==0) {
            global_useHiCu = atoi(line+i+8);
            printf("New global_useHiCu %i\n", 
                   global_useHiCu);
            continue;
        }
        if(strncmp(line+i, "nfreeze=", 8) ==0) {
            global_nFreeze = atoi(line+i+8);
            printf("New global_nFreeze %i\n", 
                   global_nFreeze);
            continue;
        }
        if(strncmp(line+i, "output=", 7) ==0) {
            global_outputLevel = atoi(line+i+7);
            printf("New global_outputLevel %i\n", global_outputLevel);
            continue;
        }
    } /* END FOR i           */
}
#endif


static void make_float_string(char* s, real x)
{
  real temp;
  int power;
  power = 0;
  temp = x;
  while(fabs(temp) > 9.999)
    {
      temp /= 10;
      power++;
    }
  while(fabs(temp) < 0.999)
    {
      temp *= 10;
      power--;
    }
  sprintf(s, "%.1fe%i", (double)temp, power);
}




static int 
parseParam(char* s)
{
/* use #define instead of more modern static const int to please old
 * compilers. */
#define MAX_BYTES 222
  char* p = s;
  char* endPtr = s + strlen(s);
  char paramName[MAX_BYTES];
  char paramValueString[MAX_BYTES];
  /* look for = */
  char* q = p;

  while((*q != '=') && (*q != 0))
    q++;
  if(*q != '=')
    {
      do_output_2(0, "error parsing string '%s': '=' not found", s);
      return -1;
    }
  /* now q points to '=' */
  memcpy(paramName, p, q-p);
  paramName[q-p] = 0;
  /*printf("paramName = '%s'\n", paramName); */
  memcpy(paramValueString, q+1, endPtr-q-1);
  paramValueString[endPtr-q-1] = 0;
  /*printf("paramValueString = '%s'\n", paramValueString); */
  if(strlen(paramName) == 0)
    {
      do_output_2(0, "error parsing string '%s': nothing found before '='", s);
      return -1;
    }
  if(strlen(paramValueString) == 0)
    {
      do_output_2(0, "error parsing string '%s': nothing found after '='", s);
      return -1;
    }

  if(strcmp(paramName, "box") == 0)
    {
      real newBoxDist = atof(paramValueString);
      if(newBoxDist <= 0)
	{
	  do_output_2(0, "error: grid param box = %f", newBoxDist);
	  return -1;
	}
      global_boxdist = newBoxDist;
      do_output_2(2, "grid parameter box      = %f", global_boxdist);
      return 0;
    }
  if(strcmp(paramName, "maxerror") == 0)
    {
      real new_maxerror;
      char ss[MAX_BYTES];
      new_maxerror = atof(paramValueString);
      if(new_maxerror <= 0)
	{
	  do_output_2(0, "error: grid param maxerror = %f", new_maxerror);
	  return -1;
	}
      global_maxerror = new_maxerror;
      make_float_string(ss, global_maxerror);
      do_output_2(2, "grid parameter maxerror = %s", ss);
      return 0;
    }
  if(strcmp(paramName, "output") == 0)
    {
      global_outputLevel = atoi(paramValueString);
      do_output_2(2, "grid parameter output   = %i", global_outputLevel);
      return 0;
    }
  if(strcmp(paramName, "nfreeze") == 0)
    {
      int new_nFreeze = atoi(paramValueString);
      if(new_nFreeze <= 0)
	{
	  do_output_2(0, "error: grid param nfreeze = %f", new_nFreeze);
	  return -1;
	}
      global_nFreeze = new_nFreeze;
      do_output_2(2, "grid parameter nfreeze  = %i", global_nFreeze);
      return 0;
    }
  if(strcmp(paramName, "nthreads") == 0)
    {
      int new_nThreads = atoi(paramValueString);
      if(new_nThreads <= 0)
	{
	  do_output_2(0, "error: grid param nthreads = %f", new_nThreads);
	  return -1;
	}
      global_nThreads = new_nThreads;
      do_output_2(2, "grid parameter nthreads = %i", global_nThreads);
      return 0;
    }
  if(strcmp(paramName, "dotest") == 0)
    {
      global_doTestIntegration = atoi(paramValueString);
      do_output_2(2, "grid parameter dotest   = %i\n", 
		  global_doTestIntegration);
      return 0;
    }
  do_output_2(0, "error in grid input: unknown parameter '%s'", paramName);
  return -1;
}

/* dftcartesianinput: called from Fortran code to allow setting different
 * parameters of the cubature code. */
void
dftcartesianinput_(const char *line, int line_len)
{
#define MAX_BYTES 222
  int inperr;
  int* inperrPtr = &inperr;
  char *endPtr, *p;
  char line2[MAX_BYTES];

  do_output_2(1, "dftcartesianinput, line_len = %i", line_len);
  if(line_len < 0)
    return;
  memset(line2, 0, MAX_BYTES);
  memcpy(line2, line, line_len);
  line2[line_len] = 0;
  /*printf("dftgridparams_ len = %i\n", line_len); */
  /*printf("line after cutting at line_len:\n'%s'\n", line2); */
  endPtr = line2 + line_len;
  p = line2;
  while(p < endPtr)
    {
      char paramBuf[MAX_BYTES];
      char* q;
        /* skip spaces */
      while(*p == ' ')
	p++;
      if(*p == 0)
	break;
      /* now we are at the beginning of some string */
      q = p;
      while((*q != ' ') && (*q != 0))
	q++;
      memcpy(paramBuf, p, q-p);
      paramBuf[q-p] = 0;
      if(parseParam(paramBuf) != 0)
	*inperrPtr++;
      p = q;
    }
}


static void 
print_box(BoxStruct* box, int prio)
{
  int i;
  do_output_2(prio, "print_box:");
  for(i = 0; i < NO_OF_DIMENSIONS; i++)
    {
      do_output_2(prio, "min = %.11f   max = %.11f", 
		  (double)box->min[i], (double)box->max[i]);
    } /* END FOR i */
}


static int 
get_distribution_box(BoxStruct* box, 
		     DistributionSpecStruct* distr, 
		     real targetRhoError)
{
  real targetError, r1, extent, arg;
  int i;
  targetError = targetRhoError;
  arg = distr->coeff / targetError;
  if(arg < 0) arg *= -1;
  if(arg < 1e-30)
    {
      do_output_2(0, "error in get_distribution_box: (arg == 0)");
      return -1;
    }
  r1 = log(arg);
  if(r1 < 0) r1 *= -1;
  extent = sqrt(r1 / distr->exponent);
  for(i = 0; i < NO_OF_DIMENSIONS; i++)
    {
      box->min[i] = distr->centerCoords[i] - extent;
      box->max[i] = distr->centerCoords[i] + extent;
    } /* END FOR i */
  return 0;
} /* END get_distribution_box */

static int 
get_shell_box(BoxStruct* box, ShellSpecStruct* shell)
{
  int i;
  for(i = 0; i < NO_OF_DIMENSIONS; i++)
    {
      box->min[i] = shell->centerCoords[i] - shell->extent;
      box->max[i] = shell->centerCoords[i] + shell->extent;
    } /* END FOR i */
  return 0;
} /* END get_shell_box */




static real 
compute_value_at_point(
		       DensitySpecStruct* density,
		       int noOfNonzeroShells,
		       int* nonZeroShellsIndexList,
		       int noOfNonzeroBasFuncs,
		       int* nonZeroBasFuncsIndexList,
		       real (*coor)[3],
		       real* workList)
{
  ShellSpecStruct* currShell;
  int i, j, iIndex, jIndex, symmetryFactor, count;
  real expFactor, result, currivalue;
  real xdiff, ydiff, zdiff;
  real x0, y0, z0;
  real x2, y2, z2, r2;
  int nbast;
  const real* dmat;

  nbast = density->nbast;
  dmat = density->dmat;

  if(noOfNonzeroBasFuncs > nbast)
    {
      do_error_exit("error in compute_integral_from_points: "
		    "(noOfNonzeroBasFuncs > nbast)\n");
      return 0;
    }

  /* compute values of contracted distributions at given point */
  count = 0;
  for(i = 0; i < noOfNonzeroShells; i++)
    {
      currShell = &density->shellList[nonZeroShellsIndexList[i]];
      x0 = currShell->centerCoords[0];
      y0 = currShell->centerCoords[1];
      z0 = currShell->centerCoords[2];

      xdiff = coor[0][0] - x0;
      ydiff = coor[0][1] - y0;
      zdiff = coor[0][2] - z0;
      x2 = xdiff * xdiff;
      y2 = ydiff * ydiff;
      z2 = zdiff * zdiff;
      r2 = x2 + y2 + z2;

      /* compute expFactor (this is the same procedure for all shell types) */
      expFactor = 0;
      for(j = 0; j < currShell->noOfContr; j++)
	  expFactor += currShell->coeffList[j] * 
	    exp(-currShell->exponentList[j] * r2);
      /* OK, expFactor computed */

      /* now there will be a different number of entries  */
      /* depending on shell type */
      switch(currShell->shellType)
	{
	case 0:
            /* 's' type shell, 1 function */
	  workList[count] = expFactor * 
	    solid_harmonic_s_0(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
	  break;
	case 1:
            /* 'p' type shell, 3 functions */
	  workList[count] = expFactor * 
	    solid_harmonic_p_0(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
	  workList[count] = expFactor * 
	    solid_harmonic_p_1(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
	  workList[count] = expFactor * 
	    solid_harmonic_p_2(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
	  break;
	case 2:
            /* 'd' type shell, 5 functions */
	  workList[count] = expFactor * 
	    solid_harmonic_d_0(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
	  workList[count] = expFactor * 
	    solid_harmonic_d_1(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
	  workList[count] = expFactor * 
	    solid_harmonic_d_2(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
	  workList[count] = expFactor * 
	    solid_harmonic_d_3(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
	  workList[count] = expFactor * 
	    solid_harmonic_d_4(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
	  break;
	case 3:
            /* 'f' type shell, 7 functions */
	  workList[count] = expFactor * 
	    solid_harmonic_f_0(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
	  workList[count] = expFactor * 
	    solid_harmonic_f_1(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
	  workList[count] = expFactor * 
	    solid_harmonic_f_2(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
	  workList[count] = expFactor * 
	    solid_harmonic_f_3(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
	  workList[count] = expFactor * 
	    solid_harmonic_f_4(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
	  workList[count] = expFactor * 
	    solid_harmonic_f_5(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
	  workList[count] = expFactor * 
	    solid_harmonic_f_6(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
	  break;
	default:
	  do_output_2(0, "error in compute_value_at_point_2: "
		      "only spdf type shells implemented");
	  do_output_2(0, "currShell->shellType = %i\n", currShell->shellType);
	  do_error_exit("error in compute_value_at_point_2");
	  return -1;
	} /* END SWITCH shellType       */
    } /* END FOR i (for each shell) */
    
  if(count > nbast)
    {
      do_error_exit("error in compute_value_at_point: (count > nbast)");
      return -1;
    }
  
  /* now use density matrix to obtain final result */
  result = 0;
  for(i = 0; i < noOfNonzeroBasFuncs; i++)
    {
      currivalue = workList[i];
      iIndex = nonZeroBasFuncsIndexList[i];
      for(j = i; j < noOfNonzeroBasFuncs; j++)
	{
	  if(j == i)
	    symmetryFactor = 1;
	  else
	    symmetryFactor = 2;
	  jIndex = nonZeroBasFuncsIndexList[j];
	  result += symmetryFactor * dmat[iIndex*nbast+jIndex] * 
	    currivalue * workList[j];
	} /* END FOR j */
    } /* END FOR i */
  
  return result;
} /* END compute_value_at_point */



static real 
compute_integral_from_points(
			     DensitySpecStruct* density,
			     int noOfNonzeroShells,
			     int* nonZeroShellsIndexList,
			     int noOfNonzeroBasFuncs,
			     int* nonZeroBasFuncsIndexList,
			     int nPoints,
			     real (*coor)[3],
			     real* weight,
			     real* workList)
{
#if 1
  int i;
  real sum;
  sum = 0;
  for(i = 0; i < nPoints; i++)
    {
      sum += compute_value_at_point(density,
				    noOfNonzeroShells,
				    nonZeroShellsIndexList,
				    noOfNonzeroBasFuncs,
				    nonZeroBasFuncsIndexList,
				    &coor[i],
				    workList) * weight[i];
    } /* END FOR i */
  return sum;
#else
  ShellSpecStruct* currShell;
  int i, j, iIndex, jIndex, symmetryFactor, count, pointNo;
  real expFactor, result;
  real xdiff, ydiff, zdiff;
  real x2, y2, z2, r2;
  real x0, y0, z0;
  real sum;
  int savedCount;
  int nbast;
  const real* dmat;

  nbast = density->nbast;
  dmat = density->dmat;

  if(noOfNonzeroBasFuncs > nbast)
    {
      do_error_exit("error in compute_integral_from_points: "
		    "(noOfNonzeroBasFuncs > nbast)\n");
      return 0;
    }

  if(nPoints > MAX_NO_OF_POINTS_PER_BATCH)
    {
      do_error_exit("error in compute_integral_from_points: "
		    "(nPoints > MAX_NO_OF_POINTS_PER_BATCH)\n");
      return 0;
    }

  memset(workList, 0x00 ,noOfNonzeroBasFuncs * sizeof(real));

  /* compute values of contracted distributions at given point */
  count = 0;
  
  for(i = 0; i < noOfNonzeroShells; i++)
    {
      currShell = &density->shellList[nonZeroShellsIndexList[i]];
      x0 = currShell->centerCoords[0];
      y0 = currShell->centerCoords[1];
      z0 = currShell->centerCoords[2];

      savedCount = count;;
      for(pointNo = 0; pointNo < nPoints; pointNo++)
	{
	  count = savedCount;

	  xdiff = coor[pointNo][0] - x0;
	  ydiff = coor[pointNo][1] - y0;
	  zdiff = coor[pointNo][2] - z0;
	  x2 = xdiff * xdiff;
	  y2 = ydiff * ydiff;
	  z2 = zdiff * zdiff;
	  r2 = x2 + y2 + z2;

	  /* compute expFactor (this is the same procedure for all
           * shell types) */
	  expFactor = 0;
	  for(j = 0; j < currShell->noOfContr; j++)
	      expFactor += currShell->coeffList[j] * 
		exp(-currShell->exponentList[j] * r2);
	  /* OK, expFactor computed */
	  /*printf("expFactor = %.22f\n", expFactor); */

	  /* now there will be a different number of entries  */
	  /* depending on shell type */
	  switch(currShell->shellType)
	    {
	    case 0:
                /* 's' type shell, 1 function */
	      workList[count*nPoints+pointNo] = expFactor * 
		solid_harmonic_s_0(xdiff, ydiff, zdiff, x2, y2, z2, r2); 
	      count++;
	      break;
	    case 1:
                /* 'p' type shell, 3 functions */
	      workList[count*nPoints+pointNo] = expFactor * 
		solid_harmonic_p_0(xdiff, ydiff, zdiff, x2, y2, z2, r2); 
	      count++;
	      workList[count*nPoints+pointNo] = expFactor * 
		solid_harmonic_p_1(xdiff, ydiff, zdiff, x2, y2, z2, r2); 
	      count++;
	      workList[count*nPoints+pointNo] = expFactor * 
		solid_harmonic_p_2(xdiff, ydiff, zdiff, x2, y2, z2, r2); 
	      count++;
	      break;
	    case 2:
                /* 'd' type shell, 5 functions */
	      workList[count*nPoints+pointNo] = expFactor * 
		solid_harmonic_d_0(xdiff, ydiff, zdiff, x2, y2, z2, r2); 
	      count++;
	      workList[count*nPoints+pointNo] = expFactor * 
		solid_harmonic_d_1(xdiff, ydiff, zdiff, x2, y2, z2, r2); 
	      count++;
	      workList[count*nPoints+pointNo] = expFactor * 
		solid_harmonic_d_2(xdiff, ydiff, zdiff, x2, y2, z2, r2); 
	      count++;
	      workList[count*nPoints+pointNo] = expFactor * 
		solid_harmonic_d_3(xdiff, ydiff, zdiff, x2, y2, z2, r2); 
	      count++;
	      workList[count*nPoints+pointNo] = expFactor * 
		solid_harmonic_d_4(xdiff, ydiff, zdiff, x2, y2, z2, r2); 
	      count++;
	      break;
	    case 3:
                /* 'f' type shell, 7 functions */
	      workList[count*nPoints+pointNo] = expFactor * 
		solid_harmonic_f_0(xdiff, ydiff, zdiff, x2, y2, z2, r2); 
	      count++;
	      workList[count*nPoints+pointNo] = expFactor * 
		solid_harmonic_f_1(xdiff, ydiff, zdiff, x2, y2, z2, r2); 
	      count++;
	      workList[count*nPoints+pointNo] = expFactor * 
		solid_harmonic_f_2(xdiff, ydiff, zdiff, x2, y2, z2, r2); 
	      count++;
	      workList[count*nPoints+pointNo] = expFactor * 
		solid_harmonic_f_3(xdiff, ydiff, zdiff, x2, y2, z2, r2); 
	      count++;
	      workList[count*nPoints+pointNo] = expFactor * 
		solid_harmonic_f_4(xdiff, ydiff, zdiff, x2, y2, z2, r2); 
	      count++;
	      workList[count*nPoints+pointNo] = expFactor * 
		solid_harmonic_f_5(xdiff, ydiff, zdiff, x2, y2, z2, r2); 
	      count++;
	      workList[count*nPoints+pointNo] = expFactor * 
		solid_harmonic_f_6(xdiff, ydiff, zdiff, x2, y2, z2, r2); 
	      count++;
	      break;
	    default:
	      do_output_2(0, "error in compute_integral_from_points: "
			  "only spdf type shells implemented");
	      do_output_2(0, "currShell->shellType = %i\n", 
			  currShell->shellType);
	      do_error_exit("error in compute_integral_from_points");
	      return -1;
	    } /* END SWITCH shellType */

	} /* END FOR pointNo */

    } /* END FOR i (for each shell) */
    
  if(count != noOfNonzeroBasFuncs)
    {
      do_error_exit("error in compute_integral_from_points: "
		    "(count != noOfNonzeroBasFuncs)");
      return -1;
    }
  
  /* now use density matrix to obtain final result */
  result = 0;
  for(i = 0; i < noOfNonzeroBasFuncs; i++)
    {
      iIndex = nonZeroBasFuncsIndexList[i];
      for(j = i; j < noOfNonzeroBasFuncs; j++)
	{
	  if(j == i)
	    symmetryFactor = 1;
	  else
	    symmetryFactor = 2;
	  jIndex = nonZeroBasFuncsIndexList[j];
	  /*real *restrict ci = &workList[i*nPoints]; */
	  /*real *restrict cj = &workList[j*nPoints]; */
	  /*real *restrict wt = weight; */
	  sum = 0;
	  for(pointNo = 0; pointNo < nPoints; pointNo++)
	    {
	      sum += workList[i*nPoints+pointNo] * 
		workList[j*nPoints+pointNo] * weight[pointNo];
	      /*sum += ci[pointNo] * cj[pointNo] * wt[pointNo]; */
	    } /* END FOR pointNo */
      
	  result += symmetryFactor * dmat[iIndex*nbast+jIndex] * sum;
	} /* END FOR j */
    } /* END FOR i */

  return result;
#endif
} /* END compute_integral_from_points */



static real 
to_power(real x, int n)
{
  real result;
  int i;
  result = 1;
  for(i = 0; i < n; i++)
    result *= x;
  return result;
}


static real 
compute_1d_gaussian_integral_recursive(real a, real b, int n, real alpha)
{
  real result, sqrtalpha, term1, term2;
  real aToPowerNminus1, bToPowerNminus1;
  if(n == 0)
    {
      sqrtalpha = sqrt(alpha);
      result = sqrt(CONSTPI/(4*alpha)) * (erf(sqrtalpha*b) - erf(sqrtalpha*a));
      return result;
    }
  if(n == 1)
    {
      result = -(1 / (2*alpha)) * (exp(-alpha*b*b) - exp(-alpha*a*a));
      return result;
    }
  if(n < 0)
    {
      do_output_2(0, "error in 1dintegral: n < 0");
      return 0;
    }
  /* now we know that n >= 2 */
  term1 = (n - 1) * compute_1d_gaussian_integral_recursive(a, b, n-2, alpha);
  aToPowerNminus1 = to_power(a, n-1);
  bToPowerNminus1 = to_power(b, n-1);
  term2  = 
    bToPowerNminus1 * exp(-alpha*b*b) - 
    aToPowerNminus1 * exp(-alpha*a*a);
  result = (term1 - term2) / (2 * alpha);
  /*  return 0; */
  return result;
} /* END compute_1d_gaussian_integral_recursive */



static real
compute_1d_gaussian_integral(real a, real b, int n, real alpha)
{
  real result, sqrtalpha, term1, term2;
  return compute_1d_gaussian_integral_recursive(a, b, n, alpha);
  result = 0;
  switch(n)
    {
    case 0:
      sqrtalpha = sqrt(alpha);
      result = sqrt(CONSTPI/(4*alpha)) * (erf(sqrtalpha*b) - erf(sqrtalpha*a));
      break;
    case 1:
      result = -(1 / (2*alpha)) * (exp(-alpha*b*b) - exp(-alpha*a*a));
      break;
    case 2:
      sqrtalpha = sqrt(alpha);
      term1 = 
	sqrt(CONSTPI/(16*alpha*alpha*alpha)) * 
	(erf(sqrtalpha*b) - erf(sqrtalpha*a));
      term2 = -(1 / (2 * alpha)) * (b*exp(-alpha*b*b) - a*exp(-alpha*a*a));
      result = term1 + term2;
      break;
    case 3:
      result = -(1 / (2*alpha*alpha)) * ((1+alpha*b*b)*exp(-alpha*b*b) - 
					 (1+alpha*a*a)*exp(-alpha*a*a));
      break;
    default:
      compute_1d_gaussian_integral_recursive(a, b, n, alpha);
      break;
    } /* END SWITCH n */
  return result;
} /* END compute_1d_gaussian_integral */


static real 
compute_integral_over_box(DistributionSpecStruct* distr, BoxStruct* box)
{
  real result, a, b, alpha;
  int i, n;
  result = distr->coeff;
  alpha = distr->exponent;
  for(i = 0; i < NO_OF_DIMENSIONS; i++)
    {
      n = distr->monomialInts[i];
      a = box->min[i] - distr->centerCoords[i];
      b = box->max[i] - distr->centerCoords[i];
      result *= compute_1d_gaussian_integral(a, b, n, alpha);
    } /* END FOR i */
  return result;
} /* END compute_integral_over_box */


static int 
get_distrs_for_box(int* resultList, rhoTreeNode* node, BoxStruct* inputBoxPtr)
{
#define MAX_DEPTH 888
  int n, i, overlap, currDepth;
  rhoTreeNode* nodeList[MAX_DEPTH];
  int statusList[MAX_DEPTH];
  rhoTreeNode* currNode;
  BoxStruct box;
  BoxStruct* currBox;

  memcpy(&box, inputBoxPtr, sizeof(BoxStruct));

  n = 0;
  currDepth = 0;
  nodeList[0] = node;
  statusList[0] = 0;
  while(currDepth >= 0)
    {
      if(statusList[currDepth] == 2)
	currDepth--;
      else
	{

	  currNode = nodeList[currDepth];
	  currBox = &currNode->box;

	  /* check for box overlap */
	  overlap = 1;
	  for(i = 0; i < NO_OF_DIMENSIONS; i++)
	    {
	      if(currBox->min[i] > box.max[i])
		overlap = 0;
	      if(currBox->max[i] < box.min[i])
		overlap = 0;
	    } /* END FOR i */
	  if(overlap == 0)
	    currDepth--;
	  else
	    {

	      if(statusList[currDepth] == 0)
		{
		  if(currNode->distrIndex >= 0)
		    {
		      resultList[n] = currNode->distrIndex;
		      n++;
		      currDepth--;
		    }
		  else
		    {
		      statusList[currDepth] = 1;
		      currDepth++;
		      statusList[currDepth] = 0;
		      nodeList[currDepth] = currNode->child1;
		    }
		} /* END IF status 0 */
	      else
		{
                    /* status is 1 */
		  statusList[currDepth] = 2;
		  currDepth++;
		  statusList[currDepth] = 0;
		  nodeList[currDepth] = currNode->child2;	      
		} /* END ELSE status 1 */
	    }
	}
    } /* END WHILE (currDepth >= 0) */

  return n;
} /* END get_distrs_for_box */

static int
use_cubature_rule(int maxlen,
		  real (*coor)[3],
		  real *weight,
		  BoxStruct* box,
		  int ruleNumber)
{
  real volume, diff0, diff1, diff2;
  real c0, c1, c2, a, b;
  real currCoords[3];
  int Ngrid, currIndex;
  int i, j, k, ii;
  real a0, a1, a2;

  volume = 1;
  for(i = 0; i < NO_OF_DIMENSIONS; i++)
    volume *= (box->max[i] - box->min[i]);
  
  switch(ruleNumber)
    {
    case 1: /* single point in center of box */
      Ngrid = 1;
      if(Ngrid >= maxlen)
	{
	  do_output_2(0, "error in use_cubature_rule: (Ngrid >= maxlen)");
	  return -1;
	}
      for(i = 0; i < NO_OF_DIMENSIONS; i++)
	  coor[0][i] = (box->max[i] + box->min[i]) / 2;
      weight[0] = volume;
      break;
    case 2: /* eight points towards corners of box */
      Ngrid = 8;
      if(Ngrid >= maxlen)
	{
	  do_output_2(0, "error in use_cubature_rule: (Ngrid >= maxlen)");
	  return -1;
	}
      for(i = 0; i < Ngrid; i++)
	weight[i] = volume / 8;
      diff0 = box->max[0] - box->min[0];
      diff1 = box->max[1] - box->min[1];
      diff2 = box->max[2] - box->min[2];
      currIndex = 0;
      for(i = 0; i < 2; i++)
	{
	  currCoords[0] = box->min[0] + 0.25*diff0 + 0.5*diff0*i;
	  for(j = 0; j < 2; j++)
	    {
	      currCoords[1] = box->min[1] + 0.25*diff1 + 0.5*diff1*i;
	      for(k = 0; k < 2; k++)
		{
		  currCoords[2] = box->min[2] + 0.25*diff2 + 0.5*diff2*i;
		  for(ii = 0; ii < 3; ii++)
		    {
		      coor[currIndex][ii] = currCoords[ii];
		    } /* END FOR ii */
		  currIndex++;
		} /* END FOR k */
	    } /* END FOR j */
	} /* END FOR i */
      break;
    case 3: /* 14 point, degree 5 rule (Stroud 1971) */
      Ngrid = 14;
      if(Ngrid >= maxlen)
	{
	  do_output_2(0, "error in use_cubature_rule: (Ngrid >= maxlen)");
	  return -1;
	}
      for(i = 0; i < 6; i++)
	weight[i] = volume *  0.88642659279778393 / 8;
      for(i = 6; i < 14; i++)
	weight[i] = volume * 0.33518005540166204 / 8;
      diff0 = box->max[0] - box->min[0];
      diff1 = box->max[1] - box->min[1];
      diff2 = box->max[2] - box->min[2];
      c0 = (box->max[0] + box->min[0]) / 2;
      c1 = (box->max[1] + box->min[1]) / 2;
      c2 = (box->max[2] + box->min[2]) / 2;
      a = 0.79582242575422146 * 0.5;
      b = 0.75878691063932814 * 0.5;

#define MACRO_3VECT(v,x,y,z) v[0]=x; v[1]=y; v[2]=z;

      MACRO_3VECT(coor[0], c0-a*diff0, c1,         c2        );
      MACRO_3VECT(coor[1], c0+a*diff0, c1,         c2        );
      MACRO_3VECT(coor[2], c0        , c1-a*diff1, c2        );
      MACRO_3VECT(coor[3], c0        , c1+a*diff1, c2        );
      MACRO_3VECT(coor[4], c0        , c1        , c2-a*diff2);
      MACRO_3VECT(coor[5], c0        , c1        , c2+a*diff2);

      MACRO_3VECT(coor[ 6], c0-b*diff0, c1-b*diff1, c2-b*diff2);
      MACRO_3VECT(coor[ 7], c0-b*diff0, c1-b*diff1, c2+b*diff2);
      MACRO_3VECT(coor[ 8], c0-b*diff0, c1+b*diff1, c2-b*diff2);
      MACRO_3VECT(coor[ 9], c0-b*diff0, c1+b*diff1, c2+b*diff2);
      MACRO_3VECT(coor[10], c0+b*diff0, c1-b*diff1, c2-b*diff2);
      MACRO_3VECT(coor[11], c0+b*diff0, c1-b*diff1, c2+b*diff2);
      MACRO_3VECT(coor[12], c0+b*diff0, c1+b*diff1, c2-b*diff2);
      MACRO_3VECT(coor[13], c0+b*diff0, c1+b*diff1, c2+b*diff2);

      break;

    case 4: /* 25 point, degree 5 rule (Stroud 1971) */
      Ngrid = 25;
      if(Ngrid >= maxlen)
	{
	  do_output_2(0, "error in use_cubature_rule: (Ngrid >= maxlen)");
	  return -1;
	}
      weight[0] = volume * 1.6842105263157894 / 8;
      for(i = 1; i < 25; i++)
	weight[i] = volume *  0.26315789473684210 / 8;

      diff0 = box->max[0] - box->min[0];
      diff1 = box->max[1] - box->min[1];
      diff2 = box->max[2] - box->min[2];
      c0 = (box->max[0] + box->min[0]) / 2;
      c1 = (box->max[1] + box->min[1]) / 2;
      c2 = (box->max[2] + box->min[2]) / 2;
      a = 0.47800981191507060 * 0.5;
      b = 0.89982215247931316 * 0.5;

      MACRO_3VECT(coor[0], c0, c1, c2);

      ii = 1;
      a0 = a;
      a1 = a;
      a2 = b;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1-a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1-a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1+a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1+a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1-a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1-a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1+a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1+a1*diff1, c2+a2*diff2); ii++;
      a0 = a;
      a1 = b;
      a2 = a;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1-a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1-a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1+a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1+a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1-a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1-a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1+a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1+a1*diff1, c2+a2*diff2); ii++;
      a0 = b;
      a1 = a;
      a2 = a;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1-a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1-a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1+a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1+a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1-a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1-a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1+a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1+a1*diff1, c2+a2*diff2); ii++;

      break;

    case 5: /* 27 point, degree 7 rule (Stroud 1971) */
      Ngrid = 27;
      if(Ngrid >= maxlen)
	{
	  do_output_2(0, "error in use_cubature_rule: (Ngrid >= maxlen)");
	  return -1;
	}
      weight[0] = volume * 0.78807348274421057 / 8;
      for(i = 1; i < 7; i++)
	weight[i] = volume *  0.49936900230772032 / 8;
      for(i = 7; i < 19; i++)
	weight[i] = volume *  0.032303742334037395 / 8;
      for(i = 19; i < 27; i++)
	weight[i] = volume *  0.47850844942512734 / 8;

      diff0 = box->max[0] - box->min[0];
      diff1 = box->max[1] - box->min[1];
      diff2 = box->max[2] - box->min[2];
      c0 = (box->max[0] + box->min[0]) / 2;
      c1 = (box->max[1] + box->min[1]) / 2;
      c2 = (box->max[2] + box->min[2]) / 2;

      MACRO_3VECT(coor[0], c0, c1, c2);
      a = 0.84841801147225245 * 0.5;
      ii = 1;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1        , c2        ); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1        , c2        ); ii++;
      MACRO_3VECT(coor[ii], c0        , c1-a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0        , c1+a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0        , c1        , c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1        , c2+a*diff2); ii++;
      a = 1.1064128986267175 * 0.5;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1-a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1+a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1-a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1+a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1        , c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1        , c2+a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1        , c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1        , c2+a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1-a*diff1, c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1-a*diff1, c2+a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1+a*diff1, c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1+a*diff1, c2+a*diff2); ii++;
      a0 = a1 = a2 = 0.65281647210169120 * 0.5;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1-a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1-a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1+a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1+a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1-a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1-a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1+a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1+a1*diff1, c2+a2*diff2); ii++;

      break;

    case 6: /* 32 point, degree 7 rule (Beckers 1992) */
      Ngrid = 32;
      if(Ngrid >= maxlen)
	{
	  do_output_2(0, "error in use_cubature_rule: (Ngrid >= maxlen)");
	  return -1;
	}
      for(i = 0; i < 6; i++)
	weight[i] = volume * 0.14098806933910446 / 8;
      for(i = 6; i < 12; i++)
	weight[i] = volume * 0.53332245896607639 / 8;
      for(i = 12; i < 24; i++)
	weight[i] = volume *  0.049451452995044458 / 8;
      for(i = 24; i < 32; i++)
	weight[i] = volume *  0.42008992427854766 / 8;

      diff0 = box->max[0] - box->min[0];
      diff1 = box->max[1] - box->min[1];
      diff2 = box->max[2] - box->min[2];
      c0 = (box->max[0] + box->min[0]) / 2;
      c1 = (box->max[1] + box->min[1]) / 2;
      c2 = (box->max[2] + box->min[2]) / 2;

      ii = 0;

      a = 1.0 * 0.5;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1        , c2        ); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1        , c2        ); ii++;
      MACRO_3VECT(coor[ii], c0        , c1-a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0        , c1+a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0        , c1        , c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1        , c2+a*diff2); ii++;
      a = 0.66289786904352112 * 0.5;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1        , c2        ); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1        , c2        ); ii++;
      MACRO_3VECT(coor[ii], c0        , c1-a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0        , c1+a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0        , c1        , c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1        , c2+a*diff2); ii++;
      a = 1.0306143700994171 * 0.5;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1-a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1+a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1-a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1+a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1        , c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1        , c2+a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1        , c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1        , c2+a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1-a*diff1, c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1-a*diff1, c2+a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1+a*diff1, c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1+a*diff1, c2+a*diff2); ii++;
      a0 = a1 = a2 = 0.66713797405746656 * 0.5;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1-a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1-a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1+a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1+a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1-a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1-a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1+a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1+a1*diff1, c2+a2*diff2); ii++;

      break;

    default:
      do_output_2(0, "error: unknown cubature rule");
      return -1;
    } /* END SWITCH */

  /*  
  real testSum = 0;
  for(i = 0; i < Ngrid; i++)
    testSum += weight[i];
  printf("testSum = %22.11f\n", testSum);
  printf("volume  = %22.11f\n\n", volume);
  */

  return Ngrid;
} /* END use_cubature_rule */


static int 
compute_grid_for_box(compute_grid_for_box_params_struct* params,
		     int maxlen,
		     real (*coor)[3],
		     real *weight,
		     BoxStruct* box,
		     real analyticalIntegralValue,
		     real* workList,
		     real* totalIntegralResult)
{
#define MAX_NO_OF_TEST_POINTS 88

  int noOfGridPoints;
  int Ngrid;
  int i;
  real Iapprox, Iexact;
  BoxStruct box1;
  BoxStruct box2;
  int bestcoord, nPoints1, nPoints2;
  real abserror, dist, maxdist, halfway;
  real IexactAbs;
  real analyticalIntegralBox1, analyticalIntegralBox2;
  int splitBox;

  /*return 0; */

  /* Define Ngrid points inside box, with corresponding weights */
  /* this is where the 'cubature rule' is used */
  Ngrid = use_cubature_rule(maxlen, coor, weight, box, CUBATURE_RULE);
  if(Ngrid <= 0)
    {
      do_output_2(0, "error in use_cubature_rule");
      return -1;
    }
  noOfGridPoints = Ngrid;

  /* compute approximate and exact integrals for box */

  Iapprox = 0;

  Iapprox = compute_integral_from_points(
					 &params->density,
					 params->noOfNonzeroShells,
					 params->nonZeroShellsIndexList,
					 params->noOfNonzeroDistributions,
					 params->nonZeroDistrIndexList,
					 Ngrid,
					 &coor[0],
					 weight,
					 workList);

  Iexact = analyticalIntegralValue;

  IexactAbs = Iexact;
  if(IexactAbs < 0) IexactAbs *= -1;

  /* compute absolute error */
  abserror = Iexact - Iapprox;
  if(abserror < 0) abserror *= -1;

  /* check if error is too large */
  splitBox = 1;
  
  if((abserror*100) < params->maxerrorPerBox)
    splitBox = 0;

  if((splitBox == 1) && (abserror < params->maxerrorPerBox))
    {
        /* it seems that the error is small enough. */
        /* however, this could be a coincidence. */
        /* to check, compare with denser grid */
      real testCoor[MAX_NO_OF_TEST_POINTS][3];
      real testWeight[MAX_NO_OF_TEST_POINTS];
      real testIapprox;
      int Ngrid2, testAbsError;

      Ngrid2 = use_cubature_rule(MAX_NO_OF_TEST_POINTS, 
				 testCoor, testWeight, box, CUBATURE_RULE_2);
      if(Ngrid2 <= 0)
	{
	  do_output_2(0, "error in use_cubature_rule");
	  return -1;
	}
      
      testIapprox = 
	compute_integral_from_points(
				     &params->density,
				     params->noOfNonzeroShells,
				     params->nonZeroShellsIndexList,
				     params->noOfNonzeroDistributions,
				     params->nonZeroDistrIndexList,
				     Ngrid2,
				     &testCoor[0],
				     testWeight,
				     workList);
      testAbsError = fabs(Iexact - testIapprox);
      /* we demand that the denser grid should give better result */
      /*printf("abserror     = %66.55f\n", abserror); */
      /*printf("testAbsError = %66.55f\n\n", testAbsError); */
      if(testAbsError <= abserror)
	splitBox = 0;
      
      /*splitBox = 0; */
    }
  if(splitBox == 1)
    {
        /* error too large, split box into box1 and box2 */

        /* first determine in which coordinate direction to do the split */
      maxdist = 0;
      bestcoord = -1;
      for(i = 0; i < NO_OF_DIMENSIONS; i++)
	{
	  dist = box->max[i] - box->min[i];
	  if(dist > maxdist)
	    {
	      maxdist = dist;
	      bestcoord = i;
	    }
	} /* END FOR i */
      if(bestcoord < 0)
	return -1;
      /* now create new boxes box1 and box2 */
      for(i = 0; i < NO_OF_DIMENSIONS; i++)
	{
	  if(i == bestcoord)
	    {
                /* direction of split */
	      halfway = (box->max[i] + box->min[i]) / 2;
	      box1.min[i] = box->min[i];
	      box1.max[i] = halfway;
	      box2.min[i] = halfway;
	      box2.max[i] = box->max[i];
	    }
	  else
	    {
                /* other direction, simply copy bounds */
	      box1.min[i] = box->min[i];
	      box1.max[i] = box->max[i];
	      box2.min[i] = box->min[i];
	      box2.max[i] = box->max[i];
	    }
	} /* END FOR i */
      /* now boxes box1 and box2 are now created */
      
      analyticalIntegralBox1 = 0;
      for(i = 0; i < params->density.noOfDistributions; i++)
	analyticalIntegralBox1 += 
	  compute_integral_over_box(&params->density.distrList[i], &box1);

      analyticalIntegralBox2 = 
	analyticalIntegralValue - analyticalIntegralBox1;

      /* create grid points for box1 */
      nPoints1 = compute_grid_for_box(params,
				      maxlen,
				      coor,
				      weight,
				      &box1,
				      analyticalIntegralBox1,
				      workList,
				      totalIntegralResult);
      if(nPoints1 < 0)
	return -1;
      /* create grid points for box2 */
      nPoints2 = compute_grid_for_box(params,
				      maxlen-nPoints1,
				      &coor[nPoints1],
				      &weight[nPoints1],
				      &box2,
				      analyticalIntegralBox2,
				      workList,
				      totalIntegralResult);
      if(nPoints2 < 0)
	return -1;
      noOfGridPoints = nPoints1 + nPoints2;
    } /* END IF error too large */
  else
    {
        /* error acceptable,  */
        /* the computed grid points for this box are good enough. */
        /* do nothing more, just return the number of points */
      *totalIntegralResult += Iapprox;
    }

  return noOfGridPoints;
} /* END compute_grid_for_box */



static rhoTreeNode* 
BuildRhoTreeBranch(int noOfDistributionsTot,
		   DistributionSpecStruct* rho_alt_1,
		   ShellSpecStruct* rho_alt_2,
		   int distrIndexListN, 
		   int* distrIndexList,
		   real targetRhoError)
{
  rhoTreeNode* newNode;
  int i, j, samePoint, n1, n2, bestCoord;
  BoxStruct tempBox;
  real currCoord, currDiff, maxDiff, limit, extent1, extent2, testCoord;
  rhoTreeNode* child1;
  rhoTreeNode* child2;
  int* tempList;
  int tempInt;

  if(distrIndexListN < 1)
    {
      printf("error in BuildRhoTreeBranch: (distrIndexListN < 1), "
	     "distrIndexListN = %i\n", distrIndexListN);
      return NULL;
    }

  newNode = dal_malloc_safe(sizeof(rhoTreeNode));

  /* compute bounding box for this node */
  if(rho_alt_1 != NULL)
    {
      if(get_distribution_box(&newNode->box, 
			      &rho_alt_1[distrIndexList[0]], 
			      targetRhoError) != 0)
	return NULL;
    }
  else
    {
      if(get_shell_box(&newNode->box, &rho_alt_2[distrIndexList[0]]) != 0)
	return NULL;
    }
  for(i = 1; i < distrIndexListN; i++)
    {
      if(rho_alt_1 != NULL)
	{
	  if(get_distribution_box(&tempBox, 
				  &rho_alt_1[distrIndexList[i]], 
				  targetRhoError) != 0)
	    return NULL;
	}
      else
	{
	  if(get_shell_box(&tempBox, &rho_alt_2[distrIndexList[i]]) != 0)
	    return NULL;
	}
      for(j = 0; j < NO_OF_DIMENSIONS; j++)
	{
	  if(tempBox.min[j] < newNode->box.min[j]) 
	    newNode->box.min[j] = tempBox.min[j];
	  if(tempBox.max[j] > newNode->box.max[j]) 
	    newNode->box.max[j] = tempBox.max[j];
	} /* END FOR j */
    } /* END FOR i */

  /* check if only one distr */
  if(distrIndexListN == 1)
    {
        /* OK, this becomes a leaf node */
      newNode->child1 = NULL;
      newNode->child2 = NULL;
      newNode->distrIndex = distrIndexList[0];
      return newNode;
    }

  /* There is more than one distribution */
  /* Get box that encloses all distributions */
  for(i = 0; i < NO_OF_DIMENSIONS; i++)
    {
      if(rho_alt_1 != NULL)
	{
	  tempBox.min[i] = rho_alt_1[distrIndexList[0]].centerCoords[i];
	  tempBox.max[i] = rho_alt_1[distrIndexList[0]].centerCoords[i];
	}
      else
	{
	  tempBox.min[i] = rho_alt_2[distrIndexList[0]].centerCoords[i];
	  tempBox.max[i] = rho_alt_2[distrIndexList[0]].centerCoords[i];
	}
    } /* END FOR i */
  for(i = 1; i < distrIndexListN; i++)
    {
      for(j = 0; j < NO_OF_DIMENSIONS; j++)
	{
	  if(rho_alt_1 != NULL)
	    currCoord = rho_alt_1[distrIndexList[i]].centerCoords[j];
	  else
	    currCoord = rho_alt_2[distrIndexList[i]].centerCoords[j];
	  if(tempBox.min[j] > currCoord) tempBox.min[j] = currCoord;
	  if(tempBox.max[j] < currCoord) tempBox.max[j] = currCoord;
	} /* END FOR j */
    } /* END FOR i */
  
  /* check if all distrs are at the same point */

  bestCoord = -1;
  maxDiff = 0;
  for(i = 0; i < NO_OF_DIMENSIONS; i++)
    {
      currDiff = tempBox.max[i] - tempBox.min[i];
      if(currDiff > maxDiff)
	{
	  bestCoord = i;
	  maxDiff = currDiff;
	}
    } /* END FOR i */
  if(bestCoord < 0)
    samePoint = 1;
  else
    {
      if(maxDiff > 1.0e-13)
	{
	  limit = (tempBox.max[bestCoord] + tempBox.min[bestCoord]) / 2;
	  samePoint = 0;
	}
      else
	samePoint = 1;
    }

  if(samePoint == 1)
    {
        /* all distrs are at the same point */
        /* sort by extent */
        /* bubble sort (this could be optimized) */
      for(i = 0; i < (distrIndexListN-1); i++)
	{
	  for(j = 0; j < (distrIndexListN-1-i); j++)
	    {
	      if(rho_alt_1 != NULL)
		{
		  extent1 = rho_alt_1[distrIndexList[j]].extent;
		  extent2 = rho_alt_1[distrIndexList[j+1]].extent;
		}
	      else
		{
		  extent1 = rho_alt_2[distrIndexList[j]].extent;
		  extent2 = rho_alt_2[distrIndexList[j+1]].extent;
		}
	      if(extent1 > extent2)
		{
                    /* do switch */
		  tempInt = distrIndexList[j];
		  distrIndexList[j] = distrIndexList[j+1];
		  distrIndexList[j+1] = tempInt;
		} /* END IF SWITCH */
	    } /* END FOR j bubble sort */
	} /* END FOR i bubble sort       */
      /* check sort */
      for(i = 0; i < (distrIndexListN-1); i++)
	{
	  if(rho_alt_1 != NULL)
	    {
	      extent1 = rho_alt_1[distrIndexList[i]].extent;
	      extent2 = rho_alt_1[distrIndexList[i+1]].extent;
	    }
	  else
	    {
	      extent1 = rho_alt_2[distrIndexList[i]].extent;
	      extent2 = rho_alt_2[distrIndexList[i+1]].extent;
	    }
	  if(extent1 > extent2)
	    {
	      do_output_2(0, 
			  "error in BuildRhoTreeBranch: list not sorted");
	      return NULL;
	    }
	} /* END FOR i check sort */

      /* create 2 new boxes: small extent and large extent */
      n1 = distrIndexListN / 2;
      n2 = distrIndexListN - n1;
    }
  else
    {
        /* all distrs are NOT at the same point */
      limit = (tempBox.max[bestCoord] + tempBox.min[bestCoord]) / 2;
      tempList = dal_malloc_safe(distrIndexListN * sizeof(int));
      n1 = 0;
      n2 = 0;
      for(i = 0; i < distrIndexListN; i++)
	{
	  if(rho_alt_1 != NULL)
	    testCoord = rho_alt_1[distrIndexList[i]].centerCoords[bestCoord];
	  else
	    testCoord = rho_alt_2[distrIndexList[i]].centerCoords[bestCoord];
	  if(testCoord > limit)
	    {
	      tempList[n1] = distrIndexList[i];
	      n1++;
	    }
	  else
	    {
	      tempList[distrIndexListN-1-n2] = distrIndexList[i];
	      n2++;
	    }
	} /* END FOR i */
      if((n1 == 0) || (n2 == 0))
	{
	  printf("error in BuildRhoTreeBranch (after split): "
		 "n1 = %i, n2 = %i\n", n1, n2);
	  printf("maxDiff = %.55f\n", maxDiff);
	  return NULL;
	}

      memcpy(distrIndexList, tempList, distrIndexListN * sizeof(int));
      dal_free(tempList);
    }
  if((n1 == 0) || (n2 == 0))
    {
      printf("error in BuildRhoTreeBranch: n1 = %i, n2 = %i\n", n1, n2);
      return NULL;
    }
  child1 = BuildRhoTreeBranch(noOfDistributionsTot, rho_alt_1, rho_alt_2,
			      n1, distrIndexList, targetRhoError);
  if(child1 == NULL)
    return NULL;
  child2 = BuildRhoTreeBranch(noOfDistributionsTot, rho_alt_1, rho_alt_2,
			      n2, distrIndexList + n1, targetRhoError);
  if(child2 == NULL)
    return NULL;
  newNode->child1 = child1;
  newNode->child2 = child2;
  newNode->distrIndex = -1;

  return newNode;
} /* END  */


static rhoTreeNode* 
BuildRhoTree(int noOfDistributions,
	     DistributionSpecStruct* rho_alt_1,
	     ShellSpecStruct* rho_alt_2,
	     real targetRhoError)
{
  rhoTreeNode* rootNode;
  int* distrIndexList;
  int i;
  real targetError, arg, r1;
  DistributionSpecStruct* distr;

  if(rho_alt_1 != NULL)
    {
        /* compute extent for each distribution in list */
      for(i = 0; i < noOfDistributions; i++)
	{
	  distr = &rho_alt_1[i];
	  targetError = distr->coeff / 1e20;
	  arg = distr->coeff / targetError;
	  r1 = log(arg);
	  if(r1 < 0) r1 *= -1;
	  distr->extent = sqrt(r1 / distr->exponent);
	} /* END FOR i */
    }

  /* set up initial index list: all distributions included */
  distrIndexList = dal_malloc_safe(noOfDistributions * sizeof(int));
  for(i = 0; i < noOfDistributions; i++)
    distrIndexList[i] = i;

  rootNode = BuildRhoTreeBranch(noOfDistributions, rho_alt_1, rho_alt_2, 
				noOfDistributions, distrIndexList, 
				targetRhoError);

  if(rootNode == NULL)
    {
      do_output_2(0, "error in BuildRhoTreeBranch");
      return NULL;
    }

  dal_free(distrIndexList);

  return rootNode;
} /* END BuildRhoTree */






static void free_rho_tree_memory(rhoTreeNode* rootNode)
{
  rhoTreeNode* child1;
  rhoTreeNode* child2;
  child1 = rootNode->child1;
  child2 = rootNode->child2;
  if(child1 != NULL)
    free_rho_tree_memory(child1);
  if(child2 != NULL)
    free_rho_tree_memory(child2);
  dal_free(rootNode);
} /* END free_rho_tree_memory */


static int round_real(real x)
{
  int x1, x2;
  real err1, err2;

  x1 = (int)x;
  x2 = x1 + 1;
  err1 = x - (real)x1;
  err2 = (real)x2 - x;
  if(err1 <= err2)
    return x1;
  else
    return x2;
}



void*
compute_grid_thread_func(void* arg)
{
    /*  char s[888]; */
  int maxNoOfPoints;
  real (*coor)[3];
  real* weight;
  real* coorx;
  real* coory;
  real* coorz;
  int noOfShells, noOfDistributions;
  int noOfNonzeroShells;
  int* nonZeroShellIndexList;
  int noOfNonzeroDistributions;
  int* nonZeroDistrIndexList;
  int currShellNo, prevShellNo, tempInt;
  int* listShlblocks;
  compute_grid_for_box_params_struct paramsStruct;
  real* workList;
  real totalIntegralResult;
  int writeResultsToFile;
  DistributionSpecStruct* rhoForSubBox;
  int noOfWrittenBatches, noOfGridPoints;
  BoxStruct startBox;
  BoxStruct subBox;
  DensitySpecStruct* density;
  compute_grid_thread_func_struct* inputParams;
  rhoTreeNode* rhoTreeRootNode;
  rhoTreeNode* rhoTreeRootNodeShells;
  int* tempList;
  int i, j, k, m, ii, jj, kk;
  int Nx, Ny, Nz;
  int nFunctions, count, nPoints, nblocks, blockStarted;
  int startShellNo, NthisWrite;
  real Iexact, maxerror;
  FILE* gridFile;
  int jobCount, assignedJobNumber;
  ShellSpecStruct* currShell;
  
  /* get hold of input params */
  inputParams = (compute_grid_thread_func_struct*)arg;
  density = inputParams->density;
  memcpy(&startBox, inputParams->startBox, sizeof(BoxStruct));
  rhoTreeRootNode = inputParams->rhoTreeRootNode;
  rhoTreeRootNodeShells = inputParams->rhoTreeRootNodeShells;
  noOfShells = density->noOfShells;
  noOfDistributions = density->noOfDistributions;
  Nx = inputParams->Nx;
  Ny = inputParams->Ny;
  Nz = inputParams->Nz;
  maxerror = inputParams->maxerror;
  gridFile = inputParams->gridFile;

  writeResultsToFile = 1;

  do_output_2(2, "thread %i entering compute_grid_thread_func..", 
	      inputParams->threadNo);

  /* allocate memory */
  maxNoOfPoints = FILE_BATCH_N;
  coor   = dal_malloc_safe(3 * maxNoOfPoints * sizeof(real));
  weight = dal_malloc_safe(maxNoOfPoints * sizeof(real));
  coorx  = dal_malloc_safe(maxNoOfPoints * sizeof(real));
  coory  = dal_malloc_safe(maxNoOfPoints * sizeof(real));
  coorz  = dal_malloc_safe(maxNoOfPoints * sizeof(real));
  nonZeroShellIndexList = dal_malloc_safe(noOfShells * sizeof(int));
  nonZeroDistrIndexList = dal_malloc_safe(density->nbast * sizeof(int));
  workList = dal_malloc_safe(density->nbast * 
			     MAX_NO_OF_POINTS_PER_BATCH * sizeof(real));
  rhoForSubBox = dal_malloc_safe(noOfDistributions * 
				 sizeof(DistributionSpecStruct));
  tempList = dal_malloc_safe(noOfDistributions * sizeof(int));
  listShlblocks = dal_malloc_safe(MAX_NO_OF_SHLBLOCKS * 2 * sizeof(int));


  /* get initial assignedJobNumber */
#ifdef USE_PTHREADS
  pthread_mutex_lock(inputParams->jobMutex);
#endif
  assignedJobNumber = *inputParams->currJobNumber;
  *inputParams->currJobNumber += N_BATCH_JOBS;
#ifdef USE_PTHREADS
  pthread_mutex_unlock(inputParams->jobMutex);
#endif
  jobCount = 0;
  noOfWrittenBatches = 0;
  totalIntegralResult = 0;
  noOfGridPoints = 0;
  i = j = k = 0;
  for(i = 0; i < Nx; i++)
    {
      for(j = 0; j < Ny; j++)
	{
	  for(k = 0; k < Nz; k++)
	    {
	      jobCount++;
	      if(jobCount >= (assignedJobNumber + N_BATCH_JOBS))
		{
                    /* get new assignedJobNumber */
#ifdef USE_PTHREADS
		  pthread_mutex_lock(inputParams->jobMutex);
#endif
		  assignedJobNumber = *inputParams->currJobNumber;
		  *inputParams->currJobNumber += N_BATCH_JOBS;
#ifdef USE_PTHREADS
		  pthread_mutex_unlock(inputParams->jobMutex);
#endif
		}
	      if(jobCount < assignedJobNumber)
		continue;
	      if(assignedJobNumber > (Nx*Ny*Nz))
		continue;
	      /*printf("i = %i, j = %i, k = %i\n", i, j, k); */

	      /* determine current sub-box */
	      subBox.min[0] = startBox.min[0] + (real)(i + 0) * 
		(startBox.max[0] - startBox.min[0]) / Nx;
	      subBox.max[0] = startBox.min[0] + (real)(i + 1) * 
		(startBox.max[0] - startBox.min[0]) / Nx;
	      subBox.min[1] = startBox.min[1] + (real)(j + 0) * 
		(startBox.max[1] - startBox.min[1]) / Ny;
	      subBox.max[1] = startBox.min[1] + (real)(j + 1) * 
		(startBox.max[1] - startBox.min[1]) / Ny;
	      subBox.min[2] = startBox.min[2] + (real)(k + 0) * 
		(startBox.max[2] - startBox.min[2]) / Nz;
	      subBox.max[2] = startBox.min[2] + (real)(k + 1) * 
		(startBox.max[2] - startBox.min[2]) / Nz;

	      /* get list of non-zero shells for current sub-box */
	      noOfNonzeroShells = get_distrs_for_box(nonZeroShellIndexList, 
						     rhoTreeRootNodeShells, 
						     &subBox);
	      if(noOfNonzeroShells < 0)
		{
		  do_output_2(0, "error in get_distrs_for_box");
		  return NULL;
		}
	      if(noOfNonzeroShells == 0)
		continue;
	      
	      /* sort list of non-zero shells (bubble sort, could be optimized) */
	      for(kk = 0; kk < (noOfNonzeroShells - 1); kk++)
		{
		  for(jj = 0; jj < (noOfNonzeroShells - 1 - kk); jj++)
		    {
		      if(nonZeroShellIndexList[jj] > 
			 nonZeroShellIndexList[jj+1])
			{
			  tempInt = nonZeroShellIndexList[jj];
 			  nonZeroShellIndexList[jj] = 
			    nonZeroShellIndexList[jj+1];
	 		  nonZeroShellIndexList[jj+1] = tempInt;
		 	}
                    } /* END FOR jj */
	 	} /* END FOR kk */
	      
	      /* translate list of nonzero shells to list of  */
	      /* nonzero contracted distributions */
	      noOfNonzeroDistributions = 0;
	      for(kk = 0; kk < noOfNonzeroShells; kk++)
		{
		  currShell = &density->shellList[nonZeroShellIndexList[kk]];
		  nFunctions = 1 + 2 * 
		    currShell->shellType;
		  for(ii = 0; ii < nFunctions; ii++)
		    {
		      nonZeroDistrIndexList[noOfNonzeroDistributions] = 
			currShell->startIndexInMatrix + ii;
		      noOfNonzeroDistributions++;
		    } /* END FOR ii      */
		} /* END FOR kk */
	      if(noOfNonzeroDistributions > density->nbast)
		{
		  do_output_2(0, "error: (noOfNonzeroDistributions > nbast)");
		  return NULL;
		}
	      
	      /* get list of relevant distributions for sub-box */
	      count = get_distrs_for_box(tempList, rhoTreeRootNode, &subBox);
	      if(count < 0)
		{
		  do_output_2(0, "error in get_distrs_for_box");
		  return NULL;
		}
	      if(count == 0)
		continue;

	      for(m = 0; m < count; m++)
		memcpy(&rhoForSubBox[m], 
		       &density->distrList[tempList[m]], 
		       sizeof(DistributionSpecStruct));

	      Iexact = 0;
	      for(kk = 0; kk < count; kk++)
		Iexact += 
		  compute_integral_over_box(&rhoForSubBox[kk], &subBox);

	      /* create grid for sub-box */

	      memcpy(&paramsStruct.density, 
		     density, 
		     sizeof(DensitySpecStruct));
	      paramsStruct.density.noOfDistributions = count;
	      paramsStruct.density.distrList = rhoForSubBox;

	      paramsStruct.maxerrorPerBox = maxerror;
	      paramsStruct.nonZeroDistrIndexList = nonZeroDistrIndexList;
	      paramsStruct.noOfNonzeroDistributions = noOfNonzeroDistributions;
	      paramsStruct.nonZeroShellsIndexList = nonZeroShellIndexList;
	      paramsStruct.noOfNonzeroShells = noOfNonzeroShells;

	      nPoints = compute_grid_for_box(&paramsStruct,
					     maxNoOfPoints,
					     &coor[0],
					     &weight[0],
					     &subBox,
					     Iexact,
					     workList,
					     &totalIntegralResult);
	      if(nPoints < 0)
		{
		  do_output_2(0, "error in compute_grid_for_box");
		  return NULL;
		}

	      if(nPoints == 0)
		continue;

	      noOfGridPoints += nPoints;

	      if(writeResultsToFile == 1)
		{
                    int nPointsLeft;
                    /* make block-list of non-zero shells to write to file */
		  nblocks = 0;
		  blockStarted = 0;
		  startShellNo = -1;
		  prevShellNo = -1;
		  for(kk = 0; kk < noOfNonzeroShells; kk++)
		    {
		      currShellNo = nonZeroShellIndexList[kk];
		      if(blockStarted == 0)
			{
			  blockStarted = 1;
			  startShellNo = currShellNo;
			}
		      else
			{
			  if(currShellNo != (prevShellNo + 1))
			    {
                                /* register previous block */
			      listShlblocks[nblocks*2] = startShellNo+1;
			      listShlblocks[nblocks*2+1] = prevShellNo+1;
			      nblocks++;
			      startShellNo = currShellNo;
			    }
			}
		      prevShellNo = currShellNo;
		    } /* END FOR kk */
		  if(blockStarted == 1)
		    {
                        /* register previous block */
		      listShlblocks[nblocks*2] = startShellNo+1;
		      listShlblocks[nblocks*2+1] = prevShellNo+1;
		      nblocks++;
		    }

		  /* set up separate x, y, z vectors for writing to file */
		  for(kk = 0; kk < nPoints; kk++)
		    {
		      coorx[kk] = coor[kk][0];
		      coory[kk] = coor[kk][1];
		      coorz[kk] = coor[kk][2];
		    }

		  /* write grid points to file */
		   nPointsLeft = nPoints;
#ifdef USE_PTHREADS
		  pthread_mutex_lock(inputParams->fileMutex);
#endif
		  while(nPointsLeft > 0)
		    {
		      if(nPointsLeft <= MAX_NO_OF_POINTS_PER_WRITE)
			NthisWrite = nPointsLeft;
		      else
			NthisWrite = MAX_NO_OF_POINTS_PER_WRITE;
		      fwrite(&NthisWrite, sizeof(int), 1, gridFile);
		      fwrite(&nblocks, sizeof(int), 1, gridFile);
		      fwrite(listShlblocks, sizeof(int), 2*nblocks, gridFile);
#if 1
		      fwrite(&(coor[nPoints-nPointsLeft][0]),
			     sizeof(real), 3*NthisWrite, gridFile);
#else
		      fwrite(&coorx[nPoints-nPointsLeft], 
			     sizeof(real), NthisWrite, gridFile);
		      fwrite(&coory[nPoints-nPointsLeft], 
			     sizeof(real), NthisWrite, gridFile);
		      fwrite(&coorz[nPoints-nPointsLeft], 
			     sizeof(real), NthisWrite, gridFile);
#endif
		      fwrite(&weight[nPoints-nPointsLeft], 
			     sizeof(real), NthisWrite, gridFile);
		      nPointsLeft -= NthisWrite;
		      noOfWrittenBatches++;		      
		    } /* END WHILE points left */
#ifdef USE_PTHREADS
		  pthread_mutex_unlock(inputParams->fileMutex);
#endif
		} /* END IF (gridFile != NULL) */
	    } /* END FOR k */
	} /* END FOR j */
    } /* END FOR i */


  do_output_2(2, "thread %i loops done, freeing memory..", 
	      inputParams->threadNo);

  /* free memory */
  dal_free(coor);
  dal_free(weight);
  dal_free(coorx);
  dal_free(coory);
  dal_free(coorz);
  dal_free(nonZeroShellIndexList);
  dal_free(nonZeroDistrIndexList);
  dal_free(workList);
  dal_free(rhoForSubBox);
  dal_free(tempList);
  dal_free(listShlblocks);

  do_output_2(2, "thread %i mem freed OK, setting result params..", 
	      inputParams->threadNo);
  
  /* report results through input structure */
  inputParams->noOfPoints = noOfGridPoints;
  inputParams->noOfWrittenBatches = noOfWrittenBatches;
  inputParams->integralResult = totalIntegralResult;

  do_output_2(2, "thread %i exiting compute_grid_thread_func", 
	      inputParams->threadNo);

  return NULL;
} /* END compute_grid_thread_func */



int
do_test_integration(DensitySpecStruct* density, char* gridFileName)
{
  real (*coor)[3];
  real* weight;
  real* coorx;
  real* coory;
  real* coorz;
  int noOfNonzeroShells;
  int* nonZeroShellIndexList;
  int noOfNonzeroDistributions;
  int* nonZeroDistrIndexList;
  int finished, currIndex;
  int* listShlblocks;
  real* workList;
  int nRepeats, repeatNo;
  real testIntegralResult;
  int maxNoOfPoints, noOfShells, noOfDistributions;
  /*  char s[888]; */
  time_t startSeconds2, endSeconds2;
  clock_t startClock, endClock;
  FILE* gridFile;
  int nPoints, nblocks, ii, kk, nFunctions;
  int nPointsDone, nPointsThisTime, nBatches;
  real secondsTakenReal;


  noOfShells = density->noOfShells;
  if(noOfShells <= 0)
    return -1;

  noOfDistributions = density->noOfDistributions;
  if(noOfDistributions <= 0)
    return -1;

  /* allocate memory */
  maxNoOfPoints = FILE_BATCH_N;
  coor   = dal_malloc_safe(3 * maxNoOfPoints * sizeof(real));
  weight = dal_malloc_safe(maxNoOfPoints * sizeof(real));
  coorx  = dal_malloc_safe(maxNoOfPoints * sizeof(real));
  coory  = dal_malloc_safe(maxNoOfPoints * sizeof(real));
  coorz  = dal_malloc_safe(maxNoOfPoints * sizeof(real));
  nonZeroShellIndexList = dal_malloc_safe(noOfShells * sizeof(int));
  nonZeroDistrIndexList = dal_malloc_safe(density->nbast * sizeof(int));
  workList = dal_malloc_safe(density->nbast * 
			     MAX_NO_OF_POINTS_PER_BATCH * sizeof(real));
  listShlblocks = dal_malloc_safe(MAX_NO_OF_SHLBLOCKS * 2 * sizeof(int));


  nRepeats = 1;
  do_output_2(2, "grid created OK, computing integral to check time, "
	      "nRepeats = %i", nRepeats);

  time(&startSeconds2);

  for(repeatNo = 0; repeatNo < nRepeats; repeatNo++)
    {
      testIntegralResult = 0;
      startClock = clock();
      gridFile = fopen(gridFileName, "rb");
      if(gridFile == NULL)
	{
	  do_output_2(0, "error opening grid file for reading");
	  return -1;
	}
      finished = 0;
      nBatches = 0;
      while(finished == 0)
	{
	  if(fread(&nPoints, sizeof(int), 1, gridFile) != 1)
	    finished = 1;
	  else
	    {
	      if(fread(&nblocks, sizeof(int), 1, gridFile) != 1) 
		{
		  do_output_2(0, "error reading nblocks");
		  return -1;
		}
	      if(fread(listShlblocks, sizeof(int), 2*nblocks, gridFile) 
		 != (2*nblocks))
		{
		  do_output_2(0, "error reading blocks");
		  return -1;
		}
	      if(fread(coor, sizeof(real), 3*nPoints, gridFile) != 3*nPoints)
		return -1;
	      /*
	      if(fread(coorx, sizeof(real), nPoints, gridFile) != nPoints)
		return -1;
	      if(fread(coory, sizeof(real), nPoints, gridFile) != nPoints)
		return -1;
	      if(fread(coorz, sizeof(real), nPoints, gridFile) != nPoints)
		return -1;
	      for(kk = 0; kk < nPoints; kk++)
		{
		  coor[kk][0] = coorx[kk];
		  coor[kk][1] = coory[kk];
		  coor[kk][2] = coorz[kk];
		}
	      */
	      if(fread(&weight[0], sizeof(real), nPoints, gridFile) != nPoints)
		return -1;
	  

	      noOfNonzeroShells = 0;
	      for(kk = 0; kk < nblocks; kk++)
		{
		  currIndex = listShlblocks[kk*2];
		  while(currIndex <= listShlblocks[kk*2+1])
		    {
		      nonZeroShellIndexList[noOfNonzeroShells] = currIndex-1;
		      noOfNonzeroShells++;
		      currIndex++;
		    } /* END WHILE */
		} /* END FOR kk */
	  
	      /* translate list of nonzero shells to list of nonzero
               * contracted distributions */
	      noOfNonzeroDistributions = 0;
	      for(kk = 0; kk < noOfNonzeroShells; kk++)
		{
		  nFunctions = 1 + 2 * density->shellList[nonZeroShellIndexList[kk]].shellType;
		  for(ii = 0; ii < nFunctions; ii++)
		    {
		      nonZeroDistrIndexList[noOfNonzeroDistributions] = 
			density->shellList[nonZeroShellIndexList[kk]].startIndexInMatrix + ii;
		      noOfNonzeroDistributions++;
		    } /* END FOR ii */
		} /* END FOR kk */

	      nPointsDone = 0;
	      while(nPointsDone < nPoints)
		{
		  if((nPoints - nPointsDone) > MAX_NO_OF_POINTS_PER_BATCH)
		    nPointsThisTime = MAX_NO_OF_POINTS_PER_BATCH;
		  else
		    nPointsThisTime = nPoints - nPointsDone;

		  testIntegralResult += 
		    compute_integral_from_points(density,
						 noOfNonzeroShells,
						 nonZeroShellIndexList,
						 noOfNonzeroDistributions,
						 nonZeroDistrIndexList,
						 nPointsThisTime,
						 &coor[nPointsDone],
						 &weight[nPointsDone],
						 workList);
		  nPointsDone += nPointsThisTime;
		} /* END WHILE (nPointsDone < nPoints) */

#if 0
		  for(kk = 0; kk < nPoints; kk++)
		    {
		      testIntegralResult += 
			compute_value_at_point_2(nonZeroShellIndexList,
						 noOfNonzeroShells,
						 shellList2,
						 noOfShells,
						 nbast,
						 dmat,
						 &coor[kk],
						 nonZeroDistrIndexList,
						 noOfNonzeroDistributions)
			* weight[kk];
		    } /* END FOR kk */
#endif

		  nBatches++;
	    } /* END ELSE nPoints read OK */
	} /* END WHILE more in file */
      endClock = clock();
      secondsTakenReal = (real)(endClock - startClock) / CLOCKS_PER_SEC;
      do_output_2(2, "integral computed, result: %.11f, "
		  "took %.3f s (time value may overflow for > 1000 s)", 
		  testIntegralResult, secondsTakenReal);
      do_output_2(2, "nBatches from test integration: %i", nBatches);
      fclose(gridFile);
    } /* END FOR repeatNo */

  time(&endSeconds2);
  do_output_2(2, "time for test integration: %i s",
	      (int)(endSeconds2 - startSeconds2));

  /* free memory */
  dal_free(coor);
  dal_free(weight);
  dal_free(workList);
  dal_free(coorx);
  dal_free(coory);
  dal_free(coorz);
  dal_free(nonZeroShellIndexList);
  dal_free(nonZeroDistrIndexList);
  dal_free(listShlblocks);  

  return 0;
} /* END do_test_integration */


int compute_grid(
		 DensitySpecStruct* density,
		 real maxerror,
		 real boxdist,
		 real targetRhoError,
		 char* gridFileName,
		 int noOfThreads,
		 int doTestIntegration
		 )
{
  FILE* gridFile;
  char ss[888];
  char sss[888];
  int noOfGridPoints;
  BoxStruct startBox;
  BoxStruct tempBox;
  int i, j;
  rhoTreeNode* rhoTreeRootNode;
  rhoTreeNode* rhoTreeRootNodeShells;
  clock_t startTime, timeTaken;
  time_t startSeconds, endSeconds;
  real Iexact, absRelError;
  int Nxyz[3]; /* Nx Ny Nz */
  int Nx, Ny, Nz;
  int IexactInteger;
  int correctValueInt;
  real correctValue, abserrorRel;
  int noOfWrittenBatches;
  real totalIntegralResult;
  int noOfDistributions, writeResultsToFile;
  int currJobNumber, noOfShells;
  real megaBytes;
  compute_grid_thread_func_struct* threadParamsList;

  do_output_2(2, "entering compute_grid..\n");

  noOfShells = density->noOfShells;
  if(noOfShells <= 0)
    return -1;

  noOfDistributions = density->noOfDistributions;
  if(noOfDistributions <= 0)
    return -1;


  make_float_string(ss, maxerror);
  make_float_string(sss, targetRhoError);
  do_output_2(2, "Entering compute_grid, noOfDistributions = %i, "
	      "maxerror = %s, targetRhoError = %s",
	      noOfDistributions, ss, sss);

  time(&startSeconds);
  do_output_2(2, "compute_grid start time: %s", ctime(&startSeconds));

  megaBytes = (real)(density->nbast * density->nbast * sizeof(real)) / 1000000;
  do_output_2(2, "nbast = %i, dmat uses %.1f Mb", density->nbast, megaBytes);  

  /* set up starting box */
  if(get_distribution_box(&startBox, 
			  &density->distrList[0], 
			  targetRhoError) != 0)
    {
      do_output_2(0, "error in get_distribution_box");
      return -1;
    }
  for(i = 1; i < noOfDistributions; i++)
    {
      if(get_distribution_box(&tempBox, 
			      &density->distrList[i], 
			      targetRhoError) != 0)
	{
	  do_output_2(0, "error in get_distribution_box");
	  return -1;
	}
      for(j = 0; j < NO_OF_DIMENSIONS; j++)
	{
	  if(tempBox.min[j] < startBox.min[j]) 
	    startBox.min[j] = tempBox.min[j];
	  if(tempBox.max[j] > startBox.max[j]) 
	    startBox.max[j] = tempBox.max[j];
	} /* END FOR j */
    } /* END FOR i */

  do_output_2(2, "compute_grid starting box:");
  print_box(&startBox, 2);

  Iexact = 0;
  for(i = 0; i < noOfDistributions; i++)
    Iexact += compute_integral_over_box(&density->distrList[i], &startBox);
  do_output_2(2, "Analytical integral over starting box: %.22f", Iexact);
  IexactInteger = round_real(Iexact);
  absRelError = fabs((double)IexactInteger - Iexact) / (double)IexactInteger;
  make_float_string(ss, absRelError);
  do_output_2(2, "Assuming that the correct value is %i, "
	      "the relative error is %s", IexactInteger, ss);

  do_output_2(2, "calling BuildRhoTree...");
  startTime = clock();
  rhoTreeRootNode = BuildRhoTree(noOfDistributions, 
				 density->distrList, 
				 NULL, 
				 targetRhoError);
  if(rhoTreeRootNode == NULL)
    {
      do_output_2(0, "error in BuildRhoTree\n");
      return -1;
    }
  timeTaken = clock() - startTime;
  do_output_2(2, "BuildRhoTree returned OK, time taken: %.2f s", 
	      ((double)timeTaken)/CLOCKS_PER_SEC);

  do_output_2(2, "calling BuildRhoTree for shells...");
  startTime = clock();
  rhoTreeRootNodeShells = BuildRhoTree(noOfShells, 
				       NULL, 
				       density->shellList, 
				       targetRhoError);
  if(rhoTreeRootNodeShells == NULL)
    {
      do_output_2(0, "error in BuildRhoTree\n");
      return -1;
    }
  timeTaken = clock() - startTime;
  do_output_2(2, "BuildRhoTree for shells returned OK, time taken: %.2f s", 
	      ((double)timeTaken)/CLOCKS_PER_SEC);

  /* compute Nx Ny Nz */
  for(i = 0; i < 3; i++)
      Nxyz[i] = 1 + (int)((startBox.max[i] - startBox.min[i]) / boxdist);
  Nx = Nxyz[0];
  Ny = Nxyz[1];
  Nz = Nxyz[2];
  do_output_2(2, "boxdist = %f, Nx = %i, Ny = %i, Nz = %i, Ntot = %i", 
	      boxdist, Nx, Ny, Nz, Nx*Ny*Nz);

  /* create grid file */
  gridFile = fopen(gridFileName, "wb");
  if(gridFile == NULL)
    {
      do_output_2(0, "error opening grid file '%s' for writing", gridFileName);
      return -1;
    }

  writeResultsToFile = 1;


  /* up to this point there is no parallellization */
  /* this is where we start to think about threading */

#ifdef USE_PTHREADS
  pthread_mutex_t fileMutex = PTHREAD_MUTEX_INITIALIZER;
  pthread_mutex_t jobMutex = PTHREAD_MUTEX_INITIALIZER;
#endif

  threadParamsList = dal_malloc_safe(noOfThreads * 
				     sizeof(compute_grid_thread_func_struct));

  currJobNumber = 1;

  for(i = 0; i < noOfThreads; i++)
    {
      threadParamsList[i].density = density;
      threadParamsList[i].rhoTreeRootNode = rhoTreeRootNode;
      threadParamsList[i].rhoTreeRootNodeShells = rhoTreeRootNodeShells;
      threadParamsList[i].maxerror = maxerror;
      threadParamsList[i].gridFile = gridFile;
      threadParamsList[i].startBox = &startBox;
      threadParamsList[i].Nx = Nx;
      threadParamsList[i].Ny = Ny;
      threadParamsList[i].Nz = Nz;
#ifdef USE_PTHREADS
      threadParamsList[i].fileMutex = &fileMutex;
      threadParamsList[i].jobMutex = &jobMutex;
#endif
      threadParamsList[i].currJobNumber = &currJobNumber;
      threadParamsList[i].noOfPoints = -1;
      threadParamsList[i].noOfWrittenBatches = 0;
      threadParamsList[i].integralResult = 0;
      threadParamsList[i].threadNo = i;
    } /* END FOR i */


#ifndef USE_PTHREADS
  do_output_2(2, "USE_PTHREADS not set, no threads created");
  if(noOfThreads != 1)
    {
      do_output_2(0, "error: cannot skip threads when (noOfThreads != 1)");
      return -1;
    }
  compute_grid_thread_func(&threadParamsList[0]);
#else
  do_output_2(2, "Starting %i threads.", noOfThreads);

  /* start threads */
  for(i = 0; i < noOfThreads; i++)
    {
      if(pthread_create(&threadParamsList[i].thread, 
			NULL, 
			compute_grid_thread_func, 
			&threadParamsList[i]) != 0)
	{
	  do_output_2(0, "Error in pthread_create for thread %i", i);
	  do_output_2(0, "waiting for already created threads..");
	  for(j = 0; j < i; j++)
	    {
	      if(pthread_join(threadParamsList[j].thread, NULL) != 0)
		{
		  do_output_2(0, "Error in pthread_join for thread %i", j);
		}	      
	    } /* END FOR j */
	  do_output_2(0, "all threads finished, returning error code");
	  return -1;
	}
    } /* END FOR i */

  do_output_2(2, "%i threads started OK.", noOfThreads);

  /* wait for threads to finish */
  for(i = 0; i < noOfThreads; i++)
    {
      if(pthread_join(threadParamsList[i].thread, NULL) != 0)
	{
	  do_output_2(0, "Error in pthread_join for thread %i", i);
	}
    } /* END FOR i */

  do_output_2(2, "all %i threads have finished:", noOfThreads);
  for(i = 0; i < noOfThreads; i++)
    do_output_2(2, "thread %2i noOfWrittenBatches = %6i", 
		i, threadParamsList[i].noOfWrittenBatches);
#endif
  
  /* now all threads have finished, check for errors */
  for(i = 0; i < noOfThreads; i++)
    {
      if(threadParamsList[i].noOfPoints < 0)
	{
	  do_output_2(0, "error in compute_grid_thread_func"
		      " for thread %i\n", i);
	  return -1;
	}
    } /* END FOR i */



  noOfGridPoints = 0;
  noOfWrittenBatches = 0;
  totalIntegralResult = 0;
  for(i = 0; i < noOfThreads; i++)
    {
      noOfGridPoints += threadParamsList[i].noOfPoints;
      noOfWrittenBatches += threadParamsList[i].noOfWrittenBatches;
      totalIntegralResult += threadParamsList[i].integralResult;
    } /* END FOR i */

  fclose(gridFile);

  do_output_2(2, "noOfWrittenBatches = %i", noOfWrittenBatches);

  time(&endSeconds);
  do_output_2(1, "compute_grid ending OK, noOfGridPoints = %i, took %i s",
	      noOfGridPoints, (int)(endSeconds - startSeconds));
  do_output_2(2, "compute_grid totalIntegralResult = %.11f", 
	      totalIntegralResult);
  correctValueInt = round_real(totalIntegralResult);
  correctValue = correctValueInt;
  abserrorRel = fabs(totalIntegralResult - correctValue) / correctValue;
  make_float_string(ss, abserrorRel);
  do_output_2(2, "Assuming that the correct value is %i,\n"
	      "the absolute relative error is %s", correctValueInt, ss);

  if(doTestIntegration == 1)
    {
      if(do_test_integration(density, gridFileName) != 0)
	{
	  do_output_2(0, "error in do_test_integration");
	  return -1;
	}
    } /* END IF */


  free_rho_tree_memory(rhoTreeRootNode);
  free_rho_tree_memory(rhoTreeRootNodeShells);

  dal_free(threadParamsList);

  time(&endSeconds);
  do_output_2(2, "compute_grid finish time: %s", ctime(&endSeconds));

  return noOfGridPoints;
} /* END compute_grid */


static int 
do_merge_sort_distrs(int n, 
		     DistributionSpecStruct* list, 
		     DistributionSpecStruct* workList)
{
    /* merge sort:  */
    /* first sort the first half, */
    /* then sort the second half, */
    /* then merge results to form final sorted list. */
  int n1, n2, nn, decision, i1, i2, i;
  DistributionSpecStruct* d1;
  DistributionSpecStruct* d2;

  if(n < 1)
    {
      do_output_2(0, "(n < 1)");
      return -1;
    }
  if(n == 1)
    return 0;
  
  n1 = n / 2;
  n2 = n - n1;

  /* sort first half */
  if(do_merge_sort_distrs(n1, list, workList) != 0)
    return -1;

  /* sort second half */
  if(do_merge_sort_distrs(n2, &list[n1], workList) != 0)
    return -1;

  /* merge results */
  nn = 0;
  i1 = 0;
  i2 = 0;
  while(nn < n)
    {
      if((i1 < n1) && (i2 < n2))
	{
            /* compare */
	  d1 = &list[i1];
	  d2 = &list[n1+i2];
	  decision = 0;
	  for(i = 0; i < 3; i++)
	    {
	      if(decision == 0)
		{
		  if(d1->monomialInts[i] != d2->monomialInts[i])
		    {
		      if(d1->monomialInts[i] > d2->monomialInts[i])
			decision = 1;
		      else
			decision = 2;
		    }
		} /* END IF (decision == 0) */
	    } /* END FOR i */
	  if(decision == 0)
	    {
                /* check exponents */
	      if(d1->exponent > d2->exponent)
		decision = 1;
	      else
		decision = 2;
	    }
	}
      else
	{
	  if(i1 == n1)
	      decision = 2;
	  else
	      decision = 1;
	}
      if(decision <= 0)
	{
	  do_output_2(0, "(decision <= 0)");
	  return -1;
	}
      if(decision == 1)
	{
	  memcpy(&workList[nn], &list[i1], sizeof(DistributionSpecStruct));
	  i1++;
	}
      else
	{
	  memcpy(&workList[nn], &list[n1+i2], sizeof(DistributionSpecStruct));
	  i2++;
	}
      nn++;
    } /* END WHILE (nn < n) */
  if(i1 != n1)
    {
      do_output_2(0, "(i1 != n1)");
      return -1;
    }
  if(i2 != n2)
    {
      do_output_2(0, "(i2 != n2)");
      return -1;
    }
  if(nn != n)
    {
      do_output_2(0, "(nn != n)");
      return -1;
    }
  memcpy(list, workList, n * sizeof(DistributionSpecStruct));
  return 0;
} /* END do_merge_sort_distrs */



static int
compute_extent_for_shells(BasisInfoStruct* basisInfo, real targetRhoError)
{
  int i;
  for(i = 0; i < basisInfo->noOfShells; i++)
    {
      ShellSpecStruct* currShell = &basisInfo->shellList[i];
      int contr = currShell->noOfContr;
      int worstIndex = -1;
      real largestExtent = 0;
      int kk;
      for(kk = 0; kk < contr; kk++)
	{
	  DistributionSpecStruct testDistr;
	  BoxStruct testBox;
	  int j;
	  real currExtent;

	  testDistr.coeff = currShell->coeffList[kk];
	  testDistr.exponent = currShell->exponentList[kk];
	  for(j = 0; j < 3; j++)
	    testDistr.centerCoords[j] = 0;
	  for(j = 0; j < 3; j++)
	    testDistr.monomialInts[j] = 0;
	  get_distribution_box(&testBox, &testDistr, targetRhoError);
          currExtent = (testBox.max[0] - testBox.min[0]) / 2;
	  if(currExtent > largestExtent)
	    {
	      largestExtent = currExtent;
	      worstIndex = kk;
	    }
	} /* END FOR kk */
      if(worstIndex < 0)
	return -1;
      currShell->extent = largestExtent;
    } /* END FOR i */
  return 0;
}





static int
get_no_of_primitives_for_density(real cutoff,
				 const real *dmat,
				 BasisInfoStruct* basisInfo)
{
#define MAX_DISTR_IN_TEMP_LIST 888

  int i, j;
  int symmetryFactor;
  int nBasisFuncs, nn;
  
  do_output_2(2, "entering function get_no_of_primitives_for_density, cutoff = %22.15f", cutoff);

  nBasisFuncs = basisInfo->noOfBasisFuncs;
  nn = 0;
  for(i = 0; i < nBasisFuncs; i++)
    {
      for(j = 0; j < nBasisFuncs; j++)
	{
	  DistributionSpecStruct tempList[MAX_DISTR_IN_TEMP_LIST];
	  int nPrimitives, k;
	  /* the matrix M is symmetric: include diagonal terms once, */
	  /* and include upper off-diagonal terms multiplied by 2 */
	  if(i == j)
              symmetryFactor = 1;
	  else
	    symmetryFactor = 2;
	  if(i > j)
	    continue;
          nPrimitives = 
	    get_product_simple_primitives(basisInfo, i,
					  basisInfo, j,
					  tempList,
					  MAX_DISTR_IN_TEMP_LIST);
	  do_output_2(3, "get_product_simple_primitives returned %i",
		      nPrimitives);
	  if(nPrimitives <= 0)
	    {
	      do_output_2(0, "error in get_product_simple_primitives");
	      return -1;
	    }
	  for(k = 0; k < nPrimitives; k++)
	    {
	      DistributionSpecStruct* currDistr = &tempList[k];
	      real Mij = dmat[i*nBasisFuncs+j];
	      real newCoeff = currDistr->coeff * Mij * symmetryFactor;
	      if(fabs(newCoeff) > cutoff)
		nn++;
	    }
	}
    }
  return nn;
}




static int
get_density(DistributionSpecStruct** rhoPtr,
	    int maxCountShellList,
	    int* noOfShellsReturn,
	    real cutoffInp, 
	    real targetRhoError,
	    int nbast, 
	    const real *dmat,
	    ShellSpecStruct* shellList,
	    BasisFuncStruct* basisFuncList)
{
#define MAX_DISTR_IN_TEMP_LIST 888
  real cutoff = cutoffInp;

  /*char s[888]; */
  int i, j, k, kk;
  DistributionSpecStruct* workList;
  DistributionSpecStruct* rhoSaved;
  real absvalue;
  real absdiff;
  real sqrtValue;
  int sameYesNo, firstIndex, count, withinLimit, resultCount;
  real coeffSum;
  int* markList;
  int symmetryFactor;
  int nBasisFuncs, nn, nNeededForRho;
  BasisInfoStruct basisInfo;
  DistributionSpecStruct* rho;

  do_output_2(2, "entering function get_density, cutoff = %22.15f", cutoff);

  if(get_shells(&basisInfo) != 0)
    {
      do_output_2(0, "error in get_shells");
      return -1;
    }
  do_output_2(2, "get_shells returned OK, number of shells: %i",
	      basisInfo.noOfShells);
  if(compute_extent_for_shells(&basisInfo, targetRhoError) != 0)
    {
      do_output_2(0, "error in compute_extent_for_shells");
      return -1;
    }
  do_output_2(2, "compute_extent_for_shells returned OK");
  if(get_basis_funcs(&basisInfo) != 0)
    {
      do_output_2(0, "error in get_basis_funcs");
      return -1;
    }
  do_output_2(2, "get_basis_funcs returned OK, number of basis funcs: %i",
	      basisInfo.noOfBasisFuncs);
  if(get_simple_primitives_all(&basisInfo) != 0)
    {
      do_output_2(0, "error in get_simple_primitives_all");
      return -1;
    }
  do_output_2(2, "get_simple_primitives_all returned OK, n = %i",
	      basisInfo.noOfSimplePrimitives);


  /* find out how much space is needed for rho */
  nNeededForRho = get_no_of_primitives_for_density(cutoff,
                                                   dmat,
                                                   &basisInfo);
  if(nNeededForRho <= 0)
    {
      do_output_2(0, "error in get_no_of_primitives_for_density");
      return -1;
    }
  
  /* allocate rho */
  rho = dal_malloc_safe(nNeededForRho * sizeof(DistributionSpecStruct));
  *rhoPtr = rho;
  
  nBasisFuncs = basisInfo.noOfBasisFuncs;
  nn = 0;
  for(i = 0; i < nBasisFuncs; i++)
    {
      for(j = 0; j < nBasisFuncs; j++)
	{
	  DistributionSpecStruct tempList[MAX_DISTR_IN_TEMP_LIST];
	  int nPrimitives, k;
            /*printf("i = %i, j = %i\n", i, j); */
	  /* the matrix M is symmetric: include diagonal terms once, */
	  /* and include upper off-diagonal terms multiplied by 2 */
	  if(i == j)
              symmetryFactor = 1;
	  else
	    symmetryFactor = 2;
	  if(i > j)
	    continue;
	  /*printf("calling get_product_simple_primitives\n"); */
          nPrimitives = 
	    get_product_simple_primitives(&basisInfo, i,
					  &basisInfo, j,
					  tempList,
					  MAX_DISTR_IN_TEMP_LIST);
	  do_output_2(3, "get_product_simple_primitives returned %i",
		      nPrimitives);
	  if(nPrimitives <= 0)
	    {
	      do_output_2(0, "error in get_product_simple_primitives");
	      return -1;
	    }
	  for(k = 0; k < nPrimitives; k++)
	    {
	      DistributionSpecStruct* currDistr = &tempList[k];
	      real Mij = dmat[i*nBasisFuncs+j];
	      /*printf("symmetryFactor = %i\n", symmetryFactor); */
	      real newCoeff = currDistr->coeff * Mij * symmetryFactor;
	      do_output_2(4, "Mij = %33.22f", Mij);
	      do_output_2(4, "currDistr->coeff = %33.22f", currDistr->coeff);
	      do_output_2(4, "newCoeff = %33.22f", newCoeff);
	      if(fabs(newCoeff) > cutoff)
		{
		  /* add to final list */
		  if(nn > nNeededForRho)
		    {
		      do_output_2(0, "error: (nn > nNeededForRho)");
		      return -1;
		    }
		  memcpy(&rho[nn], currDistr, 
			 sizeof(DistributionSpecStruct));
		  rho[nn].coeff = newCoeff;
		  nn++;  
		}
	    }
	}
    }

  *noOfShellsReturn = basisInfo.noOfShells;

  memcpy(shellList, basisInfo.shellList, 
	 basisInfo.noOfShells * sizeof(ShellSpecStruct));

  memcpy(basisFuncList, basisInfo.basisFuncList,
	 basisInfo.noOfShells * sizeof(BasisFuncStruct));	 


  do_output_2(2, "loop ended OK; list 'rho' created, nn = %i", nn);

  /* Now all distributions are stored in the list 'rho'. */
  /* The number of entries in the list is nn. */
  /* It could happen that all entries are not unique. */
  /* We want to join distributions that have the same center  */
  /* and the same exponent. */
  /* To do this, start with sorting the list by nx, ny, nz, exponent. */
  workList = dal_malloc_safe(nn * sizeof(DistributionSpecStruct));
  rhoSaved = dal_malloc_safe(nn * sizeof(DistributionSpecStruct));
  memcpy(rhoSaved, rho, nn * sizeof(DistributionSpecStruct));
  
  do_output_2(2, "calling do_merge_sort_distrs, nn = %i", nn);
  if(do_merge_sort_distrs(nn, rho, workList) != 0)
    {
      do_output_2(0, "error in do_merge_sort_distrs");
      return -1;
    }
  do_output_2(2, "do_merge_sort_distrs returned OK");
  

  do_output_2(2, "checking sort..");
  /* check that list is sorted */
  for(i = 0; i < (nn-1); i++)
    {
      if(rho[i].exponent < rho[i+1].exponent)
	{
	  sameYesNo = 1;
	  for(j = 0; j < 3; j++)
	    {
	      if(rho[i].monomialInts[j] != rho[i+1].monomialInts[j])
		sameYesNo = 0;
	    } /* END FOR j */
	  if(sameYesNo == 1)
	    {
	      printf("error: distr list NOT properly sorted\n");
	      return -1;
	    }
	}
    } /* END FOR i */
  do_output_2(2, "sort checked OK");


  markList = dal_malloc_safe(nn * sizeof(int));
  for(i = 0; i < nn; i++)
    markList[i] = 0;

  /* now go through sorted list, joining distributions where possible */
  i = 0;
  count = 0;
  firstIndex = 0;
  while(i < nn)
    {
        /* check if this entry has the same nx ny nz as current 'firstIndex' */
      sameYesNo = 1;
      for(j = 0; j < 3; j++)
	{
	  if(rho[i].monomialInts[j] != rho[firstIndex].monomialInts[j])
	    sameYesNo = 0;
	} /* END FOR j */
      /* check exponent */
      absdiff = fabs(rho[i].exponent - rho[firstIndex].exponent);
      if(absdiff > EXPONENT_DIFF_LIMIT)
	sameYesNo = 0;
      if(sameYesNo == 0)
	{
	  for(j = firstIndex; j < i; j++)
	    {
	      if(markList[j] == 0)
		{
		  markList[j] = 1;
		  /* join distrs that have centers within  */
		  /* DISTR_CENTER_DIST_LIMIT of this one */
		  coeffSum = rho[j].coeff;
		  for(k = j+1; k < i; k++)
		    {
		      withinLimit = 1;
		      for(kk = 0; kk < 3; kk++)
			{
			  absdiff = fabs(rho[j].centerCoords[kk] - 
					 rho[k].centerCoords[kk]);
			  if(absdiff > DISTR_CENTER_DIST_LIMIT)
			    withinLimit = 0;
			} /* END FOR kk */
		      if(withinLimit == 1)
			{
			  coeffSum += rho[k].coeff;
			  markList[k] = 1;
			  /*printf("found, k = %i\n", k); */
			}
		    } /* END FOR k */
		  memcpy(&workList[count], 
			 &rho[j], 
			 sizeof(DistributionSpecStruct));
		  workList[count].coeff = coeffSum;
		  count++;
		} /* END IF (markList[j] == 0) */
	    } /* END FOR j */
	  /*printf("setting firstIndex = %i\n", i); */
	  firstIndex = i;
	}
      else
	{
            /*printf("sameYesNo = 1\n"); */
	}
      i++;
    } /* END WHILE (i < nn) */
  /* take care of last part */
  for(j = firstIndex; j < nn; j++)
    {
      if(markList[j] == 0)
	{
	  markList[j] = 1;
	  /* join distrs that have centers within  */
	  /* DISTR_CENTER_DIST_LIMIT of this one */
	  coeffSum = rho[j].coeff;
	  for(k = j+1; k < nn; k++)
	    {
	      withinLimit = 1;
	      for(kk = 0; kk < 3; kk++)
		{
		  absdiff = fabs(rho[j].centerCoords[kk] - 
				 rho[k].centerCoords[kk]);
		  if(absdiff > DISTR_CENTER_DIST_LIMIT)
		    withinLimit = 0;
		} /* END FOR kk */
	      if(withinLimit == 1)
		{
		  coeffSum += rho[k].coeff;
		  markList[k] = 1;
		}
	    } /* END FOR k */
	  memcpy(&workList[count], &rho[j], sizeof(DistributionSpecStruct));
	  workList[count].coeff = coeffSum;
	  count++;
	} /* END IF (markList[j] == 0) */
    } /* END FOR j */

  for(j = 0; j < nn; j++)
    {
      if(markList[j] != 1)
	{
	  printf("error: (markList[%i] != 1)\n", j);
	  return -1;
	}
    } /* END FOR j */


  /* now move results back to list 'rho',  */
  /* skipping those that have too small coeff */
  resultCount = 0;
  for(i = 0; i < count; i++)
    {
      sqrtValue = sqrt(CONSTPI / workList[i].exponent);
      absvalue = workList[i].coeff * sqrtValue * sqrtValue * sqrtValue;
      if(absvalue < 0) absvalue *= -1;      
      if(absvalue > cutoff)
	{
	  memcpy(&rho[resultCount], 
		 &workList[i], 
		 sizeof(DistributionSpecStruct));
	  resultCount++;
	}
    } /* END FOR i */
  /*memcpy(rho, workList, count * sizeof(DistributionSpecStruct)); */

  do_output_2(2, "nn          = %9i", nn);
  do_output_2(2, "count       = %9i", count);
  do_output_2(2, "resultCount = %9i", resultCount);

  /*dal_free(list); */
  /*dal_free(indexList); */
  dal_free(workList);
  dal_free(markList);
  dal_free(rhoSaved);

  return resultCount;
} /* END read_density_file */



void
do_cartesian_grid(int nbast, const real* dmat, DftGridReader* res)
{
  DistributionSpecStruct* rho;
  ShellSpecStruct* shellList;
  BasisFuncStruct* basisFuncList;
  int noOfShells, nGridPoints, noOfDistributions;
  DensitySpecStruct density;

  /* check if new grid must be created */
  if((global_gridCount == 0) || 
     ((global_gridCount+1) % global_nFreeze) == 0)
    {

        /* check dmat */
      real maxabs = 0;
      int i;
      for(i = 0; i < nbast*nbast; i++)
	{
	  real temp = fabs(dmat[i]);
	  if(temp > maxabs)
	    maxabs = temp;
	}
      do_output_2(2, "do_cartesian_grid checking dmat: maxabs = %22.15f",
		  maxabs);

        /* allocate memory */
        /*printf("allocating memory..\n"); */
      shellList = 
	dal_malloc(global_maxNoOfShells * 
		   sizeof(ShellSpecStruct));
      basisFuncList = 
	dal_malloc(nbast * sizeof(BasisFuncStruct));
                    
      noOfDistributions = 
	get_density(&rho, 
		    global_maxNoOfShells, 
		    &noOfShells, 
		    global_distrCutoff, 
		    global_targetRhoError,
		    nbast, 
		    dmat, 
		    shellList, 
		    basisFuncList);
      if(noOfDistributions <= 0)
	{
	  fprintf(stderr, "error in read_density_file\n");
	  fort_print("error in read_density_file");
	  abort();
	  return;
	}
  
      density.noOfShells = noOfShells;
      density.shellList = shellList;
      density.nbast = nbast;
      density.dmat = dmat;
      density.basisFuncList = basisFuncList;
      density.noOfDistributions = noOfDistributions;
      density.distrList = rho;

      /* get grid */
      nGridPoints = compute_grid(&density,
				 global_maxerror,
				 global_boxdist,
				 global_targetRhoError,
				 "DALTON.QUAD",
				 global_nThreads,
				 global_doTestIntegration);
      if(nGridPoints <= 0)
	{
	  fprintf(stderr, "error in compute_grid");
	  fort_print("error in compute_grid");
	  abort();
	  return;
	}
                    
      /* free memory */
      free(shellList);
      free(rho);
      free(basisFuncList);
                    
    } /* END IF create new grid */

  global_gridCount++;
  return;
}


int output_energy_counter = 0;
FILE* energy_file = NULL;

void
output_energy_(double* energyPtr)
{
  double energy = *energyPtr;
  char s[888];

  output_energy_counter++;
  sprintf(s, "Energy %2i: %20.12f", output_energy_counter, energy);
  do_output_2(1, s);
  if(energy_file == NULL)
    {
      energy_file = fopen("energy_file.txt", "wt");
    }
  fprintf(energy_file, "%s\n", s);
  fflush(energy_file);
}
 
