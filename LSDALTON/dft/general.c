/*
C...   Copyright (c) 2015 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren. 
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
*/
/*
C...   Copyright (c) 2015 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras, T. Saue, 
C...   S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren. 
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
*/
/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* general.c:
   (c) Pawel Salek, pawsa@theochem.kth.se, 2001-08-02
   NOTES: Adding new functionals:
   a. use fun-slater.c as template.
   b. add 'extern Functional MyFunctional;' to lsdalton_functionals.h
   c. add '&MyFunctional' to available_functionals below.
   d. have a beer. Or some crackers, if you prefer.
*/

/* strictly conform to XOPEN ANSI C standard */
#if !defined(SYS_DEC)
/* XOPEN compliance is missing on old Tru64 4.0E Alphas */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1
#endif

/* Use BSD's strncasecmp(); if there is a platform that has no strncasecmp()
 * ask pawsa@theochem.kth.se for replacement */
#define _BSD_SOURCE 1

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#define __CVERSION__

#include "lsdalton_general.h"
/*#include "integrator.h"*/
#include "lsdalton_functionals.h"

/* C-wide constants */
const integer ZEROI = 0,   ONEI = 1, THREEI = 3, FOURI = 4;
const real ZEROR = 0.0, ONER = 1.0, TWOR = 2.0, FOURR = 4.0;


/* =================================================================== */
/* dftinput:

   read DFT functional from given line. The calling convention assumes
   Sun-style passing parameters. ATLAS linear algebra package
   http://www.netlib.org/atlas/ or http://math-atlas.sourceforge.net/
   contains an elaborate discuttion of character type variable passing
   conventions, in particular little bit below
   http://math-atlas.sourceforge.net/errata.html#RH7.0
*/
static char* DftConfString = NULL;
integer
FSYM(dftsetfunc)(const char* line, real *hfweight, integer *ierror)
{
  integer i, off, len;
  len=80;
  /* set the functional code printf function and HF weight setting
     functions to the dalton version that appends the output to the
     DALTON.OUT file. */
  for(i=len-1; i>=0 && isspace((integer)line[i]); i--)
    ;
  if(DftConfString) free(DftConfString);
  i++;
  for(off=0; line[off] && isspace((integer)line[off]); off++)
    ;
  DftConfString = malloc(i+1-off);
  strncpy(DftConfString, line+off, i-off); 
  DftConfString[i-off] = '\0';
  *ierror = 0;
  switch(fun_select_by_name(DftConfString, hfweight)) {    
  case FUN_OK: free(DftConfString); DftConfString = NULL; return 1; /* SUCCESS! */
  case FUN_UNKNOWN:
    printf("Unknown functional '%s'. Aborting.\n", DftConfString);
    dftlistfuncs_();
    *ierror = 1;
    break;
  case FUN_CONF_ERROR:
    printf("Functional configuration '%s' is not understood. "
	   "Aborting.\n", DftConfString);
    *ierror = 1;
    break;
  }
  free(DftConfString);
  return 0; /* failed to find */
}

/* =================================================================== */
/* Add DFT functional from given line.                                 */
static char* DftConfAddString = NULL;
integer
FSYM(dftaddfunc)(const char* line, real *hfweight)
{
  integer i, off, len;
  len=80;
  /* set the functional code printf function and HF weight setting
     functions to the dalton version that appends the output to the
     DALTON.OUT file. */
  for(i=len-1; i>=0 && isspace((integer)line[i]); i--)
    ;
  if(DftConfAddString) free(DftConfAddString);
  i++;
  for(off=0; line[off] && isspace((integer)line[off]); off++)
    ;
  DftConfAddString = malloc(i+1-off);
  strncpy(DftConfAddString, line+off, i-off); 
  DftConfAddString[i-off] = '\0';
  
  switch(fun_add_by_name(DftConfAddString, hfweight)) {    
  case FUN_OK: free(DftConfAddString); DftConfAddString = NULL; return 1; /* SUCCESS! */
  case FUN_UNKNOWN:
    printf("Unknown functional '%s'. Aborting.\n", DftConfString);
    dftlistfuncs_();
    break;
  case FUN_CONF_ERROR:
    printf("Functional configuration '%s' is not understood. "
	   "Aborting.\n", DftConfAddString);
    break;
  }
  free(DftConfAddString);
  return 0; /* failed to find */
}

/* =================================================================== */
/* transformations of functional derivatives from from unrestricted
 * ones to closed shell case. */
real
dftene_(const real *rho, const real *grad)
{
    FunDensProp dp;
    dp.rhoa  = dp.rhob  = *rho  *0.5;
    dp.grada = dp.gradb = *grad *0.5;
    dp.gradab = dp.grada*dp.gradb;
    return selected_func->func(&dp);
}

real
dfteneunres_(const real *rhoa, const real *rhob,
	     const real *grada,const real *gradb)
{
    FunDensProp dp;
    dp.rhoa  = *rhoa;
    dp.rhob  = *rhob;
    dp.grada = *grada;
    dp.gradb = *gradb;
    dp.gradab = dp.grada*dp.gradb;
    return selected_func->func(&dp);
}

void
dft_funcderiv1_(real *rho, real *grad, real *wght, real *vx)
{
    FunFirstFuncDrv drvs;
    FunDensProp dp;
    dp.rhoa  = dp.rhob  = *rho *0.5;
    dp.grada = dp.gradb = *grad*0.5;
    if(dp.rhoa<1e-13) dp.rhoa = dp.rhob = 1e-13;
    if(dp.grada<1e-13) dp.grada = dp.gradb = 1e-13;
    dp.gradab = dp.grada*dp.gradb;
    if(dp.gradab<1e-13) dp.gradab = 1e-13;
    drv1_clear(&drvs);
    selected_func->first(&drvs, *wght, &dp);
    vx[0] = drvs.df1000;
    vx[1] = drvs.df0010; 
    vx[2] = drvs.df00001;
}

void
dft_funcderiv1unres_(real *rhoa, real *rhob, real *grada, real *gradb, real *wght, real *vx)
{
    FunFirstFuncDrv drvs;
    FunDensProp dp;
    dp.rhoa  = *rhoa;
    dp.rhob  = *rhob;
    dp.grada = *grada;
    dp.gradb = *gradb;
    if(dp.rhoa<1e-13) dp.rhoa = 1e-13;
    if(dp.rhob<1e-13) dp.rhob = 1e-13;
    if(dp.grada<1e-13) dp.grada = 1e-13;
    if(dp.gradb<1e-13) dp.gradb = 1e-13;
    dp.gradab = dp.grada*dp.gradb;
    drv1_clear(&drvs);
    selected_func->first(&drvs, *wght, &dp);
    vx[0] = drvs.df1000;
    vx[1] = drvs.df0100;
    vx[2] = drvs.df0010; 
    vx[3] = drvs.df0001;
    vx[4] = drvs.df00001;
}

void
dft_funcderiv2_(real *rho, real *grad, real *wght, real *vx)
{
    FunSecondFuncDrv drvs;
    FunDensProp dp;
    dp.rhoa  = dp.rhob  = *rho *0.5;
    dp.grada = dp.gradb = *grad*0.5;
    if(dp.rhoa<1e-13) dp.rhoa = dp.rhob = 1e-13;
    if(dp.grada<1e-13) dp.grada = dp.gradb = 1e-13;
    dp.gradab = dp.grada*dp.gradb;
    drv2_clear(&drvs);
    selected_func->second(&drvs, *wght, &dp);
    vx[0] = drvs.df1000;
    vx[1] = drvs.df0100; 
    vx[2] = drvs.df0010;
    vx[3] = drvs.df0001;
    vx[4] = drvs.df00001; 
    vx[5] = drvs.df2000;
    vx[6] = drvs.df1100;
    vx[7] = drvs.df1010;
    vx[8] = drvs.df1001;
    vx[9] = drvs.df10001;
    vx[10] = drvs.df0020;
    vx[11] = drvs.df0011;
    vx[12] = drvs.df00101;
    vx[13] = drvs.df00002;
}

void
dft_funcderiv3_(real *rho, real *grad, real *wght, real *vx)
{
    FunThirdFuncDrv drvs;
    FunDensProp dp;
    dp.rhoa  = dp.rhob  = *rho *0.5;
    dp.grada = dp.gradb = *grad*0.5;
    if(dp.rhoa<1e-13) dp.rhoa = dp.rhob = 1e-13;
    if(dp.grada<1e-13) dp.grada = dp.gradb = 1e-13;
    dp.gradab = dp.grada*dp.gradb;
    drv3_clear(&drvs);
    selected_func->third(&drvs, *wght, &dp);
    vx[0] = drvs.df1000;  /* VX(1)  */
    vx[1] = drvs.df0100; 
    vx[2] = drvs.df0010;
    vx[3] = drvs.df0001;
    vx[4] = drvs.df00001; 
    vx[5] = drvs.df2000;
    vx[6] = drvs.df1100;
    vx[7] = drvs.df1010;
    vx[8] = drvs.df1001;
    vx[9] = drvs.df10001;
    vx[10] = drvs.df0020;  /* VX(11)  */
    vx[11] = drvs.df0011;
    vx[12] = drvs.df00101;
    vx[13] = drvs.df00002;
    vx[14] = drvs.df3000;
    vx[15] = drvs.df2100;
    vx[16] = drvs.df2010;
    vx[17] = drvs.df2001;
    vx[18] = drvs.df20001;
    vx[19] =  drvs.df1110;
    vx[20] =  drvs.df11001; /* VX(21)  */
    vx[21] =  drvs.df1020;
    vx[22] =  drvs.df1011;
    vx[23] =  drvs.df0120;
    vx[24] =  drvs.df0030;
    vx[25] =  drvs.df0021;
    vx[26] =  drvs.df00003;
}

extern void lsquit_(const char* str, integer len);
void
dalton_quit(const char* format, ...)
{
    char line[128];
    integer len;
    va_list a;
 
    va_start(a, format);
    vsnprintf(line, sizeof(line), format, a); 
    va_end(a);
    len = strlen(line);
    lsquit_(line, -1);
}

/* Helper functions. Could be bracketed with #ifdef DEBUG or something */
extern void FSYM2(lsfort_wrt)(integer* lupri, const char* str, const integer* len, integer ln);
integer
lsfort_print(integer lupri, const char* format, ...)
{
    char line[128];
    integer len;
    integer l;
    va_list a;

    va_start(a, format);
    vsnprintf(line, sizeof(line), format, a); 
    va_end(a);
    len = l = strlen(line);
    /*    printf(line,&len,l);*/
    FSYM2(lsfort_wrt)(&lupri,line,&len, l);
    return len;
}

 
