/*
C...   Copyright (c) 2005 by the authors of Dalton (see below).
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
C...   Copyright (c) 2005 by the authors of Dalton (see below).
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
/* functionals.c:
   Program for computing and testing functional routines in the DFT module.
   (c) Pawel Salek, pawsa@theochem.kth.se, 2001-08-02

 */

#define _BSD_SOURCE 1

#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "lsdalton_general.h"
#define __CVERSION__

#include "lsdalton_functionals.h"

Functional* available_functionals[] = {
    /* generic functionals */
    &BeckeFunctional,
    &B97_1Functional,
    &ExampleFunctional,
    &Hcth120Functional,
    &Hcth147Functional,
    &Hcth407Functional,
    &Hcth93Functional,
    &KTFunctional,
    &LB94Functional,
    &LYPFunctional,
    &OPTXFunctional,
    &P86cFunctional,
    &PW86xFunctional,
    &Pw91cFunctional,
    &PZ81Functional,
    &PBEcFunctional,
    &PbexFunctional,
    &revPBExFunctional,
    &rPBExFunctional,
    &mPBExFunctional,
    &PW91xFunctional,
    &G96xFunctional,
    &LG93xFunctional,
    &B88X_Functional,
    &LDAX_Functional,
    &PBEX_Functional,
    &REVPBEX_Functional,
    &RPBEX_Functional,
    &MPBEX_Functional,
    &PW91X_Functional,
    &G96X_Functional,
    &LG93X_Functional,
    &KT1X_Functional,
    &KT2X_Functional,
    &KT3X_Functional,
    &OPTX_Functional,
    &SlaterFunctional,
    &VWN3Functional,
    &VWN5Functional,
    &VWNIFunctional,
    &VWNFunctional,
    &XAlphaFunctional,
    /* mixed functionals */
    &B3LYPFunctional,
    &B3LYPGaussFunctional,
    &B3P86Functional,
    &B3P86GFunctional,
    &BLYPFunctional,
    &BP86Functional,
    &BPW91Functional,
    &Camb3lypFunctional,
    &CamxFunctional,
    &CamcompxFunctional,
    &GGAKeyFunctional,
    &CombineFunctional,
    &KT1Functional,
    &KT2Functional,
    &KT3Functional,
    &LDAFunctional,
    &OLYPFunctional,
    &PBE0Functional,
    &PBEFunctional,
    &SVWN3Functional,
    &SVWN5Functional,
    NULL
};

integer 
clear_funclist();

Functional* selected_func = &LDAFunctional;
/* =================================================================== */
enum FunError
fun_select_by_name(const char *conf_string, real *hfweight)
{
    integer ok, i;
    char func_name[20];
    clear_funclist();
    sscanf(conf_string,"%20s", func_name);
    for(i=0; available_functionals[i]; i++)
        if(strcasecmp(available_functionals[i]->name, func_name)==0) {
            selected_func = available_functionals[i];
            ok = selected_func->read ?
                selected_func->read(conf_string+strlen(func_name),hfweight) : 1;
            return ok ? FUN_OK : FUN_CONF_ERROR;
        }
    return FUN_UNKNOWN;
}

Functional* added_func = &LDAFunctional;
/* =================================================================== */
enum FunError
fun_add_by_name(const char *conf_string, real *hfweight)
{
    integer ok, i;
    char func_name[20];
    sscanf(conf_string,"%20s", func_name);
    for(i=0; available_functionals[i]; i++)
        if(strcasecmp(available_functionals[i]->name, func_name)==0) {
            added_func = available_functionals[i];
            ok = added_func->read ?
                added_func->read(conf_string+strlen(func_name),hfweight) : 1;
            return ok ? FUN_OK : FUN_CONF_ERROR;
        }
    return FUN_UNKNOWN;
}

void
drv1_clear(FunFirstFuncDrv* gga)
{
    memset(gga, 0, sizeof(*gga));
}

void
drv2_clear(FunSecondFuncDrv* gga)
{
    memset(gga, 0, sizeof(*gga));
}

void
drv3_clear(FunThirdFuncDrv* gga)
{
    memset(gga, 0, sizeof(*gga));
}
void
drv4_clear(FunFourthFuncDrv* gga)
{
    memset(gga, 0, sizeof(*gga));
}

integer fun_true(void)  { return 1; }
integer fun_false(void) { return 0; }

/* Fortran interface. We specify different names suffixes so that
 * library can be linked with code compiled with different compilers
 * or different compilation options. */

void
funset(const char *str, integer *info, integer len)
{
  real dummy; 
  dummy = 0.000;
  switch(fun_select_by_name(str,&dummy)) {
  case FUN_OK:         *info = 0; break;
  case FUN_UNKNOWN   : *info = 1; break;
  case FUN_CONF_ERROR: *info = 2; break;
  }
}
void
funset_(const char *str, integer *info, integer len)
{
  real dummy ;
  dummy = 0.000;
  switch(fun_select_by_name(str,&dummy)) {
  case FUN_OK:         *info = 0; break;
  case FUN_UNKNOWN   : *info = 1; break;
  case FUN_CONF_ERROR: *info = 2; break;
  }
}

integer
funisgga(void)
{ return selected_func->is_gga(); }
integer
funisgga_(void)
{ return selected_func->is_gga(); }

real
funenergy(const FunDensProp *dp)
{ return selected_func->func(dp); }
real
funenergy_(const FunDensProp *dp)
{ return selected_func->func(dp); }

void
funfirst(FunFirstFuncDrv *df, real *factor, const FunDensProp *dp)
{ selected_func->first(df, *factor, dp); }
void
funfirst_(FunFirstFuncDrv *df, real *factor, const FunDensProp *dp)
{ selected_func->first(df, *factor, dp); }

/* =================================================================== */
/*           fortran (and not only) functional stub routines           */
/* =================================================================== */

/* dftreport:
   report the selected functional and its configuration.
*/
void
dftreport_(integer *lupri)
{
  lsfort_print(*lupri,"\n     This is a DFT calculation of type: %s",
               selected_func->name);
    if(selected_func->report)
        selected_func->report(*lupri);
}

/* check if for the functional used the empirical disp. corr. is defined */
void   
dftdispcheck_(void)
{
    int res = 0;
    const char *bp86  = "BP86" ;
    const char *blyp  = "BLYP" ;
    const char *pbe   = "PBE"  ;
    const char *b3lyp = "B3LYP";
    const char *tpss  = "TPSS" ;
    if (strcasecmp(selected_func->name, bp86 )==0) res = 1;
    if (strcasecmp(selected_func->name, blyp )==0) res = 1;
    if (strcasecmp(selected_func->name, pbe  )==0) res = 1;
    if (strcasecmp(selected_func->name, b3lyp)==0) res = 1;
    if (strcasecmp(selected_func->name, tpss )==0) res = 1;

    if (res == 0) {
       printf("###########################################################################################################\n");
       printf("The empirical dispersion correction is only defined for the following functionals:\n\n");
       printf("   BP86\n   BLYP\n   PBE\n   B3LYP\n   TPSS\n\n");
       printf("At the moment these functionals have to be specified in the DALTON.INP with these names explicitly,\n"); 
       printf("defining the functionals by specifying the exchange and correlation parts explicitly is not yet possible.\n");
       printf("TPSS is not implemented yet.\n");
       printf("###########################################################################################################\n");
       dalton_quit("Empirical dispersion correction not defined for the choosen functional");
    }
}

/* give the s6 factor for the functional, needed in the empirical disp. corr. */
void
disp_funcfac_(real *res)
{
    *res = 0.0;
    const char *bp86  = "BP86" ;
    const char *blyp  = "BLYP" ;
    const char *pbe   = "PBE"  ;
    const char *b3lyp = "B3LYP";
    const char *tpss  = "TPSS" ;
    if (strcasecmp(selected_func->name, bp86 )==0) *res = 1.05;
    if (strcasecmp(selected_func->name, blyp )==0) *res = 1.20;
    if (strcasecmp(selected_func->name, pbe  )==0) *res = 0.75;
    if (strcasecmp(selected_func->name, b3lyp)==0) *res = 1.05;
    if (strcasecmp(selected_func->name, tpss )==0) *res = 1.00;
}

void
dftlistfuncs_(void)
{
    integer i;
    printf("\nAvailable functionals:");
    for(i=0; available_functionals[i]; i++)
      printf(" %s ",available_functionals[i]->name);
}

/* declare both known fortran name-mangled variants */
integer
dft_isgga_(void)
{ return selected_func->is_gga(); }

integer
dft_isgga__(void)
{ return selected_func->is_gga(); }






