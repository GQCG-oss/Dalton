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


/* DFT-D2 Functional Dependent Parameter Setup */

/* check if for the functional used the DFT-D2 corr. is defined */
integer
dft_d2_check_(void)
{
    integer res = 0;
    /*fun_printf("\n     Check Disp: This is a DFT calculation of type: %s",
 *  *  *                selected_func->name);*/
    if (strcasecmp(selected_func->name, "BP86")==0) res = 1;
    if (strcasecmp(selected_func->name, "BLYP")==0) res = 1;
    if (strcasecmp(selected_func->name, "PBE")==0) res = 1;
    if (strcasecmp(selected_func->name, "B3LYP")==0) res = 1;
    if (strcasecmp(selected_func->name, "B97-D")==0) res = 1;
    if (strcasecmp(selected_func->name, "B2PLYP")==0) res = 1;
    if (strcasecmp(selected_func->name, "B2GPPLYP")==0) res = 1;
    if (strcasecmp(selected_func->name, "PBE0")==0) res = 1;
    /* Functionals parameterized but not setup in Dalton */
    if (strcasecmp(selected_func->name, "REVPBE")==0) res = 1;
    if (strcasecmp(selected_func->name, "PW6B95")==0) res = 1;
    if (strcasecmp(selected_func->name, "TPSS")==0) res = 1;
    if (strcasecmp(selected_func->name, "TPSS0")==0) res = 1;
    if (strcasecmp(selected_func->name, "DSD-BLYP")==0) res = 1;
    /*Now quit outside based on results*/
    /*if (res == 0) { 
 *  *  *        dalton_quit("DFT-D2 correction not defined for the chosen functional"); 
 *   *   *            }*/
    return res;
}


/* give the s6 factor for the functional, needed in the DFT-D2 corr. */
real
dft_d2_s6_(void)
{
    real res = 0.0;

    if (strcasecmp(selected_func->name, "BP86")==0) res = 1.05;
    if (strcasecmp(selected_func->name, "BLYP")==0) res = 1.20;
    if (strcasecmp(selected_func->name, "PBE")==0) res = 0.75;
    if (strcasecmp(selected_func->name, "B3LYP")==0) res = 1.05;
    if (strcasecmp(selected_func->name, "B97-D")==0) res = 1.25;
    if (strcasecmp(selected_func->name, "B2PLYP")==0) res = 0.55;
    if (strcasecmp(selected_func->name, "B2GPPLYP")==0) res = 0.4;
    if (strcasecmp(selected_func->name, "PBE0")==0) res = 0.6;
    /* Functionals parameterized but not setup in Dalton */
    if (strcasecmp(selected_func->name, "REVPBE")==0) res = 1.25;
    if (strcasecmp(selected_func->name, "PW6B95")==0) res = 0.5;
    if (strcasecmp(selected_func->name, "TPSS")==0) res = 1.0;
    if (strcasecmp(selected_func->name, "TPSS0")==0) res = 0.85;
    if (strcasecmp(selected_func->name, "DSD-BLYP")==0) res = 0.41;

    /* Now quit outside based on results*/
    /*if (res == 0) {
 *  *       dalton_quit("DFT-D2 s_6 not defined for the chosen functional");
 *   *        *           }*/
    return res;
}


/* give the alpha factor for the functional, needed in the DFT-D2 corr. */
real
dft_d2_alp_(void)
{
    real res = 0.0;

    if (strcasecmp(selected_func->name, "BP86")==0) res = 20.0;
    if (strcasecmp(selected_func->name, "BLYP")==0) res = 20.0;
    if (strcasecmp(selected_func->name, "PBE")==0) res = 20.0;
    if (strcasecmp(selected_func->name, "B3LYP")==0) res = 20.0;
    if (strcasecmp(selected_func->name, "B97-D")==0) res = 20.0;
    if (strcasecmp(selected_func->name, "B2PLYP")==0) res = 20.0;
    if (strcasecmp(selected_func->name, "B2GPPLYP")==0) res = 20.0;
    if (strcasecmp(selected_func->name, "PBE0")==0) res = 20.0;
    /* Functionals parameterized but not setup in Dalton */
    if (strcasecmp(selected_func->name, "REVPBE")==0) res = 20.0;
    if (strcasecmp(selected_func->name, "PW6B95")==0) res = 20.0;
    if (strcasecmp(selected_func->name, "TPSS")==0) res = 20.0;
    if (strcasecmp(selected_func->name, "TPSS0")==0) res = 20.0;
    if (strcasecmp(selected_func->name, "DSD-BLYP")==0) res = 60.0;

    if (res == 0) {
      dalton_quit("DFT-D2 alpha not defined for the chosen functional");
    }
    return res;
}


/* give the alpha factor for the functional, needed in the DFT-D2 corr. */
real
dft_d2_rs6_(void)
{
    real res = 0.0;

    if (strcasecmp(selected_func->name, "BP86")==0) res = 1.10;
    if (strcasecmp(selected_func->name, "BLYP")==0) res = 1.10;
    if (strcasecmp(selected_func->name, "PBE")==0) res = 1.10;
    if (strcasecmp(selected_func->name, "B3LYP")==0) res = 1.10;
    if (strcasecmp(selected_func->name, "B97-D")==0) res = 1.10;
    if (strcasecmp(selected_func->name, "B2PLYP")==0) res = 1.10;
    if (strcasecmp(selected_func->name, "B2GPPLYP")==0) res = 1.10;
    if (strcasecmp(selected_func->name, "PBE0")==0) res = 1.10;
    /* Functionals parameterized but not setup in Dalton */
    if (strcasecmp(selected_func->name, "REVPBE")==0) res = 1.10;
    if (strcasecmp(selected_func->name, "PW6B95")==0) res = 1.10;
    if (strcasecmp(selected_func->name, "TPSS")==0) res = 1.10;
    if (strcasecmp(selected_func->name, "TPSS0")==0) res = 1.10;
    if (strcasecmp(selected_func->name, "DSD-BLYP")==0) res = 1.10;

    if (res == 0) {
      dalton_quit("DFT-D2 rs_6 not defined for the chosen functional");
    }
    return res;
}
/* DFT-D2 Functional Dependent Parameter End*/



/* DFT-D3 Functional Dependent Parameter Setup */
integer
dft_d3_check_(void)
{
    integer res = 0;
    /*fun_printf("\n     Check Disp: This is a DFT calculation of type: %s",
 *  *  *                selected_func->name);*/
    if (strcasecmp(selected_func->name, "Slater")==0)    res = 1;
    if (strcasecmp(selected_func->name, "BLYP")==0)      res = 1;
    if (strcasecmp(selected_func->name, "BP86")==0)      res = 1;
    if (strcasecmp(selected_func->name, "B97-D")==0)     res = 1;
    if (strcasecmp(selected_func->name, "rev-PBE")==0)   res = 1;
    if (strcasecmp(selected_func->name, "PBE")==0)       res = 1;
    if (strcasecmp(selected_func->name, "PBEsol")==0)    res = 1;
    if (strcasecmp(selected_func->name, "rpw86-pbe")==0) res = 1;
    if (strcasecmp(selected_func->name, "rPBE")==0)      res = 1;
    if (strcasecmp(selected_func->name, "TPSS")==0)      res = 1;
    if (strcasecmp(selected_func->name, "B3LYP")==0)     res = 1;
    if (strcasecmp(selected_func->name, "PBE0")==0)      res = 1;
    if (strcasecmp(selected_func->name, "REBPBE38")==0)  res = 1;
    if (strcasecmp(selected_func->name, "PW6B95")==0)    res = 1;
    if (strcasecmp(selected_func->name, "TPSS0")==0)     res = 1;
    if (strcasecmp(selected_func->name, "B2-PLYP")==0)   res = 1;
    if (strcasecmp(selected_func->name, "PWPB95")==0)    res = 1;
    if (strcasecmp(selected_func->name, "B2GP-PLYP")==0) res = 1;
    if (strcasecmp(selected_func->name, "PTPSS")==0)     res = 1;
    if (strcasecmp(selected_func->name, "HF")==0)        res = 1;
    if (strcasecmp(selected_func->name, "MPWLYP")==0)    res = 1;
    if (strcasecmp(selected_func->name, "BPBE")==0)      res = 1;
    if (strcasecmp(selected_func->name, "BHLYP")==0)     res = 1;
    if (strcasecmp(selected_func->name, "TPSSh")==0)     res = 1;
    if (strcasecmp(selected_func->name, "PWB6K")==0)     res = 1;
    if (strcasecmp(selected_func->name, "B1B95")==0)     res = 1;
    if (strcasecmp(selected_func->name, "BOP")==0)       res = 1;
    if (strcasecmp(selected_func->name, "OLYP")==0)      res = 1;
    if (strcasecmp(selected_func->name, "OPBE")==0)      res = 1;
    if (strcasecmp(selected_func->name, "SSB")==0)       res = 1;
    if (strcasecmp(selected_func->name, "revSSB")==0)    res = 1;
    if (strcasecmp(selected_func->name, "OTPSS")==0)     res = 1;
    if (strcasecmp(selected_func->name, "B3PW91")==0)    res = 1;
    if (strcasecmp(selected_func->name, "revPBE0")==0)   res = 1;
    if (strcasecmp(selected_func->name, "PBE38")==0)     res = 1;
    if (strcasecmp(selected_func->name, "MPW1B95")==0)   res = 1;
    if (strcasecmp(selected_func->name, "MPWB1K")==0)    res = 1;
    if (strcasecmp(selected_func->name, "BMK")==0)       res = 1;
    if (strcasecmp(selected_func->name, "CAM-B3LYP")==0) res = 1;
    if (strcasecmp(selected_func->name, "LC-wPBE")==0)   res = 1;
    if (strcasecmp(selected_func->name, "M05")==0)       res = 1;
    if (strcasecmp(selected_func->name, "M052X")==0)     res = 1;
    if (strcasecmp(selected_func->name, "M06-L")==0)     res = 1;
    if (strcasecmp(selected_func->name, "M06")==0)       res = 1;
    if (strcasecmp(selected_func->name, "M062X")==0)     res = 1;
    if (strcasecmp(selected_func->name, "M06HF")==0)     res = 1;
    if (strcasecmp(selected_func->name, "DFTB3")==0)     res = 1;
    if (strcasecmp(selected_func->name, "HCTH120")==0)   res = 1;
    /* Now handle this quit outside*/
    /*if (res == 0) {
 *  *  *       dalton_quit("DFT-D3 alpha not defined for the chosen functional");
 *   *   *           }*/
    return res;
}


/* give the s6 factor for the functional, needed in the DFT-D3 corr. */
real
dft_d3_s6_(void)
{
    real res = 0.0;

    if (strcasecmp(selected_func->name, "Slater")==0)    res = 1.0;
    if (strcasecmp(selected_func->name, "BLYP")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "BP86")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "B97-D")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "rev-PBE")==0)   res = 1.0;
    if (strcasecmp(selected_func->name, "PBE")==0)       res = 1.0;
    if (strcasecmp(selected_func->name, "PBEsol")==0)    res = 1.0;
    if (strcasecmp(selected_func->name, "rpw86-pbe")==0) res = 1.0;
    if (strcasecmp(selected_func->name, "rPBE")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "TPSS")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "B3LYP")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "PBE0")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "REBPBE38")==0)  res = 1.0;
    if (strcasecmp(selected_func->name, "PW6B95")==0)    res = 1.0;
    if (strcasecmp(selected_func->name, "TPSS0")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "B2-PLYP")==0)   res = 0.64;
    if (strcasecmp(selected_func->name, "PWPB95")==0)    res = 0.82;
    if (strcasecmp(selected_func->name, "B2GP-PLYP")==0) res = 0.56;
    if (strcasecmp(selected_func->name, "PTPSS")==0)     res = 0.75;
    if (strcasecmp(selected_func->name, "HF")==0)        res = 1.0;
    if (strcasecmp(selected_func->name, "MPWLYP")==0)    res = 1.0;
    if (strcasecmp(selected_func->name, "BPBE")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "BHLYP")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "TPSSh")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "PWB6K")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "B1B95")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "BOP")==0)       res = 1.0;
    if (strcasecmp(selected_func->name, "OLYP")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "OPBE")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "SSB")==0)       res = 1.0;
    if (strcasecmp(selected_func->name, "revSSB")==0)    res = 1.0;
    if (strcasecmp(selected_func->name, "OTPSS")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "B3PW91")==0)    res = 1.0;
    if (strcasecmp(selected_func->name, "revPBE0")==0)   res = 1.0;
    if (strcasecmp(selected_func->name, "PBE38")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "MPW1B95")==0)   res = 1.0;
    if (strcasecmp(selected_func->name, "MPWB1K")==0)    res = 1.0;
    if (strcasecmp(selected_func->name, "BMK")==0)       res = 1.0;
    if (strcasecmp(selected_func->name, "CAM-B3LYP")==0) res = 1.0;
    if (strcasecmp(selected_func->name, "LC-wPBE")==0)   res = 1.0;
    if (strcasecmp(selected_func->name, "M05")==0)       res = 1.0;
    if (strcasecmp(selected_func->name, "M052X")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "M06-L")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "M06")==0)       res = 1.0;
    if (strcasecmp(selected_func->name, "M062X")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "M06HF")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "DFTB3")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "HCTH120")==0)   res = 1.0;

    if (res == 0.0) {
      dalton_quit("DFT-D3 s_6 not defined for the chosen functional");
    }
    return res;
}



/* give the alpha factor for the functional, needed in the DFT-D3 corr. */
real
dft_d3_alp_(void)
{
    real res = 0.0;

    if (strcasecmp(selected_func->name, "Slater")==0)    res = 14.0;
    if (strcasecmp(selected_func->name, "BLYP")==0)      res = 14.0;
    if (strcasecmp(selected_func->name, "BP86")==0)      res = 14.0;
    if (strcasecmp(selected_func->name, "B97-D")==0)     res = 14.0;
    if (strcasecmp(selected_func->name, "rev-PBE")==0)   res = 14.0;
    if (strcasecmp(selected_func->name, "PBE")==0)       res = 14.0;
    if (strcasecmp(selected_func->name, "PBEsol")==0)    res = 14.0;
    if (strcasecmp(selected_func->name, "rpw86-pbe")==0) res = 14.0;
    if (strcasecmp(selected_func->name, "rPBE")==0)      res = 14.0;
    if (strcasecmp(selected_func->name, "TPSS")==0)      res = 14.0;
    if (strcasecmp(selected_func->name, "B3LYP")==0)     res = 14.0;
    if (strcasecmp(selected_func->name, "PBE0")==0)      res = 14.0;
    if (strcasecmp(selected_func->name, "REBPBE38")==0)  res = 14.0;
    if (strcasecmp(selected_func->name, "PW6B95")==0)    res = 14.0;
    if (strcasecmp(selected_func->name, "TPSS0")==0)     res = 14.0;
    if (strcasecmp(selected_func->name, "B2-PLYP")==0)   res = 14.0;
    if (strcasecmp(selected_func->name, "PWPB95")==0)    res = 14.0;
    if (strcasecmp(selected_func->name, "B2GP-PLYP")==0) res = 14.0;
    if (strcasecmp(selected_func->name, "PTPSS")==0)     res = 14.0;
    if (strcasecmp(selected_func->name, "HF")==0)        res = 14.0;
    if (strcasecmp(selected_func->name, "MPWLYP")==0)    res = 14.0;
    if (strcasecmp(selected_func->name, "BPBE")==0)      res = 14.0;
    if (strcasecmp(selected_func->name, "BHLYP")==0)     res = 14.0;
    if (strcasecmp(selected_func->name, "TPSSh")==0)     res = 14.0;
    if (strcasecmp(selected_func->name, "PWB6K")==0)     res = 14.0;
    if (strcasecmp(selected_func->name, "B1B95")==0)     res = 14.0;
    if (strcasecmp(selected_func->name, "BOP")==0)       res = 14.0;
    if (strcasecmp(selected_func->name, "OLYP")==0)      res = 14.0;
    if (strcasecmp(selected_func->name, "OPBE")==0)      res = 14.0;
    if (strcasecmp(selected_func->name, "SSB")==0)       res = 14.0;
    if (strcasecmp(selected_func->name, "revSSB")==0)    res = 14.0;
    if (strcasecmp(selected_func->name, "OTPSS")==0)     res = 14.0;
    if (strcasecmp(selected_func->name, "B3PW91")==0)    res = 14.0;
    if (strcasecmp(selected_func->name, "revPBE0")==0)   res = 14.0;
    if (strcasecmp(selected_func->name, "PBE38")==0)     res = 14.0;
    if (strcasecmp(selected_func->name, "MPW1B95")==0)   res = 14.0;
    if (strcasecmp(selected_func->name, "MPWB1K")==0)    res = 14.0;
    if (strcasecmp(selected_func->name, "BMK")==0)       res = 14.0;
    if (strcasecmp(selected_func->name, "CAM-B3LYP")==0) res = 14.0;
    if (strcasecmp(selected_func->name, "LC-wPBE")==0)   res = 14.0;
    if (strcasecmp(selected_func->name, "M05")==0)       res = 14.0;
    if (strcasecmp(selected_func->name, "M052X")==0)     res = 14.0;
    if (strcasecmp(selected_func->name, "M06-L")==0)     res = 14.0;
    if (strcasecmp(selected_func->name, "M06")==0)       res = 14.0;
    if (strcasecmp(selected_func->name, "M062X")==0)     res = 14.0;
    if (strcasecmp(selected_func->name, "M06HF")==0)     res = 14.0;
    if (strcasecmp(selected_func->name, "DFTB3")==0)     res = 14.0;
    if (strcasecmp(selected_func->name, "HCTH120")==0)   res = 14.0;

    if (res == 0) {
      dalton_quit("DFT-D3 alpha not defined for the chosen functional");
    }
    return res;
}


/* give the rs18 factor for the functional, needed in the DFT-D3 corr. */
real
dft_d3_rs18_(void)
{
    real res = 0.0;

    if (strcasecmp(selected_func->name, "Slater")==0)    res = 0.687;
    if (strcasecmp(selected_func->name, "BLYP")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "BP86")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "B97-D")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "rev-PBE")==0)   res = 1.0;
    if (strcasecmp(selected_func->name, "PBE")==0)       res = 1.0;
    if (strcasecmp(selected_func->name, "PBEsol")==0)    res = 1.0;
    if (strcasecmp(selected_func->name, "rpw86-pbe")==0) res = 1.0;
    if (strcasecmp(selected_func->name, "rPBE")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "TPSS")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "B3LYP")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "PBE0")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "REBPBE38")==0)  res = 1.0;
    if (strcasecmp(selected_func->name, "PW6B95")==0)    res = 1.0;
    if (strcasecmp(selected_func->name, "TPSS0")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "B2-PLYP")==0)   res = 1.0;
    if (strcasecmp(selected_func->name, "PWPB95")==0)    res = 1.0;
    if (strcasecmp(selected_func->name, "B2GP-PLYP")==0) res = 1.0;
    if (strcasecmp(selected_func->name, "PTPSS")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "HF")==0)        res = 1.0;
    if (strcasecmp(selected_func->name, "MPWLYP")==0)    res = 1.0;
    if (strcasecmp(selected_func->name, "BPBE")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "BHLYP")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "TPSSh")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "PWB6K")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "B1B95")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "BOP")==0)       res = 1.0;
    if (strcasecmp(selected_func->name, "OLYP")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "OPBE")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "SSB")==0)       res = 1.0;
    if (strcasecmp(selected_func->name, "revSSB")==0)    res = 1.0;
    if (strcasecmp(selected_func->name, "OTPSS")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "B3PW91")==0)    res = 1.0;
    if (strcasecmp(selected_func->name, "revPBE0")==0)   res = 1.0;
    if (strcasecmp(selected_func->name, "PBE38")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "MPW1B95")==0)   res = 1.0;
    if (strcasecmp(selected_func->name, "MPWB1K")==0)    res = 1.0;
    if (strcasecmp(selected_func->name, "BMK")==0)       res = 1.0;
    if (strcasecmp(selected_func->name, "CAM-B3LYP")==0) res = 1.0;
    if (strcasecmp(selected_func->name, "LC-wPBE")==0)   res = 1.0;
    if (strcasecmp(selected_func->name, "M05")==0)       res = 1.0;
    if (strcasecmp(selected_func->name, "M052X")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "M06-L")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "M06")==0)       res = 1.0;
    if (strcasecmp(selected_func->name, "M062X")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "M06HF")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "DFTB3")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "HCTH120")==0)   res = 1.0;

    if (res == 0.0) {
      dalton_quit("DFT-D3 rs_18 not defined for the chosen functional");
    }
    return res;
}


/* give the rs6 factor for the functional, needed in the DFT-D3 corr. */
real
dft_d3_rs6_(void)
{
    real res = 0.0;

    if (strcasecmp(selected_func->name, "Slater")==0)    res = 0.999;
    if (strcasecmp(selected_func->name, "BLYP")==0)      res = 1.094;
    if (strcasecmp(selected_func->name, "BP86")==0)      res = 1.139;
    if (strcasecmp(selected_func->name, "B97-D")==0)     res = 0.892;
    if (strcasecmp(selected_func->name, "rev-PBE")==0)   res = 0.923;
    if (strcasecmp(selected_func->name, "PBE")==0)       res = 1.217;
    if (strcasecmp(selected_func->name, "PBEsol")==0)    res = 1.345;
    if (strcasecmp(selected_func->name, "rpw86-pbe")==0) res = 1.224;
    if (strcasecmp(selected_func->name, "rPBE")==0)      res = 0.872;
    if (strcasecmp(selected_func->name, "TPSS")==0)      res = 1.166;
    if (strcasecmp(selected_func->name, "B3LYP")==0)     res = 1.261;
    if (strcasecmp(selected_func->name, "PBE0")==0)      res = 1.287;
    if (strcasecmp(selected_func->name, "REBPBE38")==0)  res = 1.021;
    if (strcasecmp(selected_func->name, "PW6B95")==0)    res = 1.532;
    if (strcasecmp(selected_func->name, "TPSS0")==0)     res = 1.252;
    if (strcasecmp(selected_func->name, "B2-PLYP")==0)   res = 1.427;
    if (strcasecmp(selected_func->name, "PWPB95")==0)    res = 1.557;
    if (strcasecmp(selected_func->name, "B2GP-PLYP")==0) res = 1.586;
    if (strcasecmp(selected_func->name, "PTPSS")==0)     res = 1.541;
    if (strcasecmp(selected_func->name, "HF")==0)        res = 1.158;
    if (strcasecmp(selected_func->name, "MPWLYP")==0)    res = 1.239;
    if (strcasecmp(selected_func->name, "BPBE")==0)      res = 1.087;
    if (strcasecmp(selected_func->name, "BHLYP")==0)     res = 1.370;
    if (strcasecmp(selected_func->name, "TPSSh")==0)     res = 1.223;
    if (strcasecmp(selected_func->name, "PWB6K")==0)     res = 1.660;
    if (strcasecmp(selected_func->name, "B1B95")==0)     res = 1.613;
    if (strcasecmp(selected_func->name, "BOP")==0)       res = 0.929;
    if (strcasecmp(selected_func->name, "OLYP")==0)      res = 0.806;
    if (strcasecmp(selected_func->name, "OPBE")==0)      res = 0.837;
    if (strcasecmp(selected_func->name, "SSB")==0)       res = 1.215;
    if (strcasecmp(selected_func->name, "revSSB")==0)    res = 1.221;
    if (strcasecmp(selected_func->name, "OTPSS")==0)     res = 1.128;
    if (strcasecmp(selected_func->name, "B3PW91")==0)    res = 1.176;
    if (strcasecmp(selected_func->name, "revPBE0")==0)   res = 0.949;
    if (strcasecmp(selected_func->name, "PBE38")==0)     res = 1.333;
    if (strcasecmp(selected_func->name, "MPW1B95")==0)   res = 1.605;
    if (strcasecmp(selected_func->name, "MPWB1K")==0)    res = 1.671;
    if (strcasecmp(selected_func->name, "BMK")==0)       res = 1.931;
    if (strcasecmp(selected_func->name, "CAM-B3LYP")==0) res = 1.378;
    if (strcasecmp(selected_func->name, "LC-wPBE")==0)   res = 1.355;
    if (strcasecmp(selected_func->name, "M05")==0)       res = 1.373;
    if (strcasecmp(selected_func->name, "M052X")==0)     res = 1.417;
    if (strcasecmp(selected_func->name, "M06-L")==0)     res = 1.581;
    if (strcasecmp(selected_func->name, "M06")==0)       res = 1.325;
    if (strcasecmp(selected_func->name, "M062X")==0)     res = 1.619;
    if (strcasecmp(selected_func->name, "M06HF")==0)     res = 1.446;
    if (strcasecmp(selected_func->name, "DFTB3")==0)     res = 1.235;
    if (strcasecmp(selected_func->name, "HCTH120")==0)   res = 1.221;

    if (res == 0) {
      dalton_quit("DFT-D3 rs_6 not defined for the chosen functional");
    }
    return res;
}



/* give the s18 factor for the functional, needed in the DFT-D3 corr. */
real
dft_d3_s18_(void)
{
    real res = 0.0;

    if (strcasecmp(selected_func->name, "Slater")==0)    res = -1.957;
    if (strcasecmp(selected_func->name, "BLYP")==0)      res = 1.682;
    if (strcasecmp(selected_func->name, "BP86")==0)      res = 1.683;
    if (strcasecmp(selected_func->name, "B97-D")==0)     res = 0.909;
    if (strcasecmp(selected_func->name, "rev-PBE")==0)   res = 1.010;
    if (strcasecmp(selected_func->name, "PBE")==0)       res = 0.722;
    if (strcasecmp(selected_func->name, "PBEsol")==0)    res = 0.612;
    if (strcasecmp(selected_func->name, "rpw86-pbe")==0) res = 0.901;
    if (strcasecmp(selected_func->name, "rPBE")==0)      res = 0.514;
    if (strcasecmp(selected_func->name, "TPSS")==0)      res = 1.105;
    if (strcasecmp(selected_func->name, "B3LYP")==0)     res = 1.703;
    if (strcasecmp(selected_func->name, "PBE0")==0)      res = 0.928;
    if (strcasecmp(selected_func->name, "REBPBE38")==0)  res = 0.862;
    if (strcasecmp(selected_func->name, "PW6B95")==0)    res = 0.862;
    if (strcasecmp(selected_func->name, "TPSS0")==0)     res = 1.242;
    if (strcasecmp(selected_func->name, "B2-PLYP")==0)   res = 1.022;
    if (strcasecmp(selected_func->name, "PWPB95")==0)    res = 0.705;
    if (strcasecmp(selected_func->name, "B2GP-PLYP")==0) res = 0.760;
    if (strcasecmp(selected_func->name, "PTPSS")==0)     res = 0.879;
    if (strcasecmp(selected_func->name, "HF")==0)        res = 1.746;
    if (strcasecmp(selected_func->name, "MPWLYP")==0)    res = 1.098;
    if (strcasecmp(selected_func->name, "BPBE")==0)      res = 2.033;
    if (strcasecmp(selected_func->name, "BHLYP")==0)     res = 1.442;
    if (strcasecmp(selected_func->name, "TPSSh")==0)     res = 1.219;
    if (strcasecmp(selected_func->name, "PWB6K")==0)     res = 0.550;
    if (strcasecmp(selected_func->name, "B1B95")==0)     res = 1.868;
    if (strcasecmp(selected_func->name, "BOP")==0)       res = 1.975;
    if (strcasecmp(selected_func->name, "OLYP")==0)      res = 1.764;
    if (strcasecmp(selected_func->name, "OPBE")==0)      res = 2.055;
    if (strcasecmp(selected_func->name, "SSB")==0)       res = 0.663;
    if (strcasecmp(selected_func->name, "revSSB")==0)    res = 0.560;
    if (strcasecmp(selected_func->name, "OTPSS")==0)     res = 1.494;
    if (strcasecmp(selected_func->name, "B3PW91")==0)    res = 1.775;
    if (strcasecmp(selected_func->name, "revPBE0")==0)   res = 0.792;
    if (strcasecmp(selected_func->name, "PBE38")==0)     res = 0.998;
    if (strcasecmp(selected_func->name, "MPW1B95")==0)   res = 1.118;
    if (strcasecmp(selected_func->name, "MPWB1K")==0)    res = 1.061;
    if (strcasecmp(selected_func->name, "BMK")==0)       res = 2.168;
    if (strcasecmp(selected_func->name, "CAM-B3LYP")==0) res = 1.217;
    if (strcasecmp(selected_func->name, "LC-wPBE")==0)   res = 1.279;
    if (strcasecmp(selected_func->name, "M05")==0)       res = 0.595;
    if (strcasecmp(selected_func->name, "M052X")==0)     res = 0.000;
    if (strcasecmp(selected_func->name, "M06-L")==0)     res = 0.000;
    if (strcasecmp(selected_func->name, "M06")==0)       res = 0.000;
    if (strcasecmp(selected_func->name, "M062X")==0)     res = 0.000;
    if (strcasecmp(selected_func->name, "M06HF")==0)     res = 0.000;
    if (strcasecmp(selected_func->name, "DFTB3")==0)     res = 0.673;
    if (strcasecmp(selected_func->name, "HCTH120")==0)   res = 1.206;

    if (res == 0) {
      dalton_quit("DFT-D3 s_18 not defined for the chosen functional");
    }
    return res;
}
/* DFT-D3 Functional Dependent Parameter End */



/* DFT-D3-BJ Functional Dependent Parameter Setup */
integer
dft_d3bj_check_(void)
{
    integer res = 0;
    /*fun_printf("\n     Check Disp: This is a DFT calculation of type: %s",
 *  *  *                selected_func->name);*/
    if (strcasecmp(selected_func->name, "BP86")==0)        res = 1;
    if (strcasecmp(selected_func->name, "BLYP")==0)        res = 1;
    if (strcasecmp(selected_func->name, "revPBE")==0)      res = 1;
    if (strcasecmp(selected_func->name, "B97-D")==0)       res = 1;
    if (strcasecmp(selected_func->name, "PBE")==0)         res = 1;
    if (strcasecmp(selected_func->name, "rpw86-pbe")==0)   res = 1;
    if (strcasecmp(selected_func->name, "B3LYP")==0)       res = 1;
    if (strcasecmp(selected_func->name, "TPSS")==0)        res = 1;
    if (strcasecmp(selected_func->name, "HF")==0)          res = 1;
    if (strcasecmp(selected_func->name, "TPSS0")==0)       res = 1;
    if (strcasecmp(selected_func->name, "PBE0")==0)        res = 1;
    if (strcasecmp(selected_func->name, "revPBE38")==0)    res = 1;
    if (strcasecmp(selected_func->name, "PW6B95")==0)      res = 1;
    if (strcasecmp(selected_func->name, "B2PLYP")==0)      res = 1;
    if (strcasecmp(selected_func->name, "DSD-BLYP")==0)    res = 1;
    if (strcasecmp(selected_func->name, "DSD-BLYP-FC")==0) res = 1;
    if (strcasecmp(selected_func->name, "BOP")==0)         res = 1;
    if (strcasecmp(selected_func->name, "MPWLYP")==0)      res = 1;
    if (strcasecmp(selected_func->name, "OLYP")==0)        res = 1;
    if (strcasecmp(selected_func->name, "PBEsol")==0)      res = 1;
    if (strcasecmp(selected_func->name, "BPBE")==0)        res = 1;
    if (strcasecmp(selected_func->name, "OPBE")==0)        res = 1;
    if (strcasecmp(selected_func->name, "SSB")==0)         res = 1;
    if (strcasecmp(selected_func->name, "revSSB")==0)      res = 1;
    if (strcasecmp(selected_func->name, "OTPSS")==0)       res = 1;
    if (strcasecmp(selected_func->name, "B3PW91")==0)      res = 1;
    if (strcasecmp(selected_func->name, "BHLYP")==0)       res = 1;
    if (strcasecmp(selected_func->name, "revPBE0")==0)     res = 1;
    if (strcasecmp(selected_func->name, "TPSSh")==0)       res = 1;
    if (strcasecmp(selected_func->name, "MPW1B95")==0)     res = 1;
    if (strcasecmp(selected_func->name, "PWB6K")==0)       res = 1;
    if (strcasecmp(selected_func->name, "B1B95")==0)       res = 1;
    if (strcasecmp(selected_func->name, "BMK")==0)         res = 1;
    if (strcasecmp(selected_func->name, "CAM-B3LYP")==0)   res = 1;
    if (strcasecmp(selected_func->name, "LC-wPBE")==0)     res = 1;
    if (strcasecmp(selected_func->name, "B2GP-PLYP")==0)   res = 1;
    if (strcasecmp(selected_func->name, "PTPSS")==0)       res = 1;
    if (strcasecmp(selected_func->name, "PWPB95")==0)      res = 1;
    if (strcasecmp(selected_func->name, "HCTH120")==0)     res = 1;
    if (strcasecmp(selected_func->name, "DFTB3")==0)       res = 1;
    /* Now handle this quit outside */
    /*if (res == 0) {
 *  *  *       dalton_quit("DFT-D3 alpha not defined for the chosen functional");
 *   *   *           }*/
    return res;
}


/* give the s6 factor for the functional, needed in the DFT-D3-BJ corr. */
real
dft_d3bj_s6_(void)
{
    real res = 0.0;

    if (strcasecmp(selected_func->name, "BP86")==0)        res = 1.0;
    if (strcasecmp(selected_func->name, "BLYP")==0)        res = 1.0;
    if (strcasecmp(selected_func->name, "revPBE")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "B97-D")==0)       res = 1.0;
    if (strcasecmp(selected_func->name, "PBE")==0)         res = 1.0;
    if (strcasecmp(selected_func->name, "rpw86-pbe")==0)   res = 1.0;
    if (strcasecmp(selected_func->name, "B3LYP")==0)       res = 1.0;
    if (strcasecmp(selected_func->name, "TPSS")==0)        res = 1.0;
    if (strcasecmp(selected_func->name, "HF")==0)          res = 1.0;
    if (strcasecmp(selected_func->name, "TPSS0")==0)       res = 1.0;
    if (strcasecmp(selected_func->name, "PBE0")==0)        res = 1.0;
    if (strcasecmp(selected_func->name, "revPBE38")==0)    res = 1.0;
    if (strcasecmp(selected_func->name, "PW6B95")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "B2PLYP")==0)      res = 0.64;
    if (strcasecmp(selected_func->name, "DSD-BLYP")==0)    res = 0.50;
    if (strcasecmp(selected_func->name, "DSD-BLYP-FC")==0) res = 0.50;
    if (strcasecmp(selected_func->name, "BOP")==0)         res = 1.0;
    if (strcasecmp(selected_func->name, "MPWLYP")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "OLYP")==0)        res = 1.0;
    if (strcasecmp(selected_func->name, "PBEsol")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "BPBE")==0)        res = 1.0;
    if (strcasecmp(selected_func->name, "OPBE")==0)        res = 1.0;
    if (strcasecmp(selected_func->name, "SSB")==0)         res = 1.0;
    if (strcasecmp(selected_func->name, "revSSB")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "OTPSS")==0)       res = 1.0;
    if (strcasecmp(selected_func->name, "B3PW91")==0)      res = 1.0;
    if (strcasecmp(selected_func->name, "BHLYP")==0)       res = 1.0;
    if (strcasecmp(selected_func->name, "revPBE0")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "TPSSh")==0)       res = 1.0;
    if (strcasecmp(selected_func->name, "MPW1B95")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "PWB6K")==0)       res = 1.0;
    if (strcasecmp(selected_func->name, "B1B95")==0)       res = 1.0;
    if (strcasecmp(selected_func->name, "BMK")==0)         res = 1.0;
    if (strcasecmp(selected_func->name, "CAM-B3LYP")==0)   res = 1.0;
    if (strcasecmp(selected_func->name, "LC-wPBE")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "B2GP-PLYP")==0)   res = 0.560;
    if (strcasecmp(selected_func->name, "PTPSS")==0)       res = 0.750;
    if (strcasecmp(selected_func->name, "PWPB95")==0)      res = 0.820;
    if (strcasecmp(selected_func->name, "HCTH120")==0)     res = 1.0;
    if (strcasecmp(selected_func->name, "DFTB3")==0)       res = 1.0;

    if (res == 0) {
      dalton_quit("DFT-D3-BJ s_6 not defined for the chosen functional");
    }
    return res;
}



/* give the alpha factor for the functional, needed in the DFT-D3-BJ corr. */
real
dft_d3bj_alp_(void)
{
    real res = 0.0;

    if (strcasecmp(selected_func->name, "BP86")==0)        res = 14.0;
    if (strcasecmp(selected_func->name, "BLYP")==0)        res = 14.0;
    if (strcasecmp(selected_func->name, "revPBE")==0)      res = 14.0;
    if (strcasecmp(selected_func->name, "B97-D")==0)       res = 14.0;
    if (strcasecmp(selected_func->name, "PBE")==0)         res = 14.0;
    if (strcasecmp(selected_func->name, "rpw86-pbe")==0)   res = 14.0;
    if (strcasecmp(selected_func->name, "B3LYP")==0)       res = 14.0;
    if (strcasecmp(selected_func->name, "TPSS")==0)        res = 14.0;
    if (strcasecmp(selected_func->name, "HF")==0)          res = 14.0;
    if (strcasecmp(selected_func->name, "TPSS0")==0)       res = 14.0;
    if (strcasecmp(selected_func->name, "PBE0")==0)        res = 14.0;
    if (strcasecmp(selected_func->name, "revPBE38")==0)    res = 14.0;
    if (strcasecmp(selected_func->name, "PW6B95")==0)      res = 14.0;
    if (strcasecmp(selected_func->name, "B2PLYP")==0)      res = 14.0;
    if (strcasecmp(selected_func->name, "DSD-BLYP")==0)    res = 14.0;
    if (strcasecmp(selected_func->name, "DSD-BLYP-FC")==0) res = 14.0;
    if (strcasecmp(selected_func->name, "BOP")==0)         res = 14.0;
    if (strcasecmp(selected_func->name, "MPWLYP")==0)      res = 14.0;
    if (strcasecmp(selected_func->name, "OLYP")==0)        res = 14.0;
    if (strcasecmp(selected_func->name, "PBEsol")==0)      res = 14.0;
    if (strcasecmp(selected_func->name, "BPBE")==0)        res = 14.0;
    if (strcasecmp(selected_func->name, "OPBE")==0)        res = 14.0;
    if (strcasecmp(selected_func->name, "SSB")==0)         res = 14.0;
    if (strcasecmp(selected_func->name, "revSSB")==0)      res = 14.0;
    if (strcasecmp(selected_func->name, "OTPSS")==0)       res = 14.0;
    if (strcasecmp(selected_func->name, "B3PW91")==0)      res = 14.0;
    if (strcasecmp(selected_func->name, "BHLYP")==0)       res = 14.0;
    if (strcasecmp(selected_func->name, "revPBE0")==0)     res = 14.0;
    if (strcasecmp(selected_func->name, "TPSSh")==0)       res = 14.0;
    if (strcasecmp(selected_func->name, "MPW1B95")==0)     res = 14.0;
    if (strcasecmp(selected_func->name, "PWB6K")==0)       res = 14.0;
    if (strcasecmp(selected_func->name, "B1B95")==0)       res = 14.0;
    if (strcasecmp(selected_func->name, "BMK")==0)         res = 14.0;
    if (strcasecmp(selected_func->name, "CAM-B3LYP")==0)   res = 14.0;
    if (strcasecmp(selected_func->name, "LC-wPBE")==0)     res = 14.0;
    if (strcasecmp(selected_func->name, "B2GP-PLYP")==0)   res = 14.0;
    if (strcasecmp(selected_func->name, "PTPSS")==0)       res = 14.0;
    if (strcasecmp(selected_func->name, "PWPB95")==0)      res = 14.0;
    if (strcasecmp(selected_func->name, "HCTH120")==0)     res = 14.0;
    if (strcasecmp(selected_func->name, "DFTB3")==0)       res = 14.0;

    if (res == 0) {
      dalton_quit("DFT-D3-BJ alpha not defined for the chosen functional");
    }
    return res;
}


/* give the rs18 factor for the functional, needed in the DFT-D3-BJ corr. */
real
dft_d3bj_rs18_(void)
{
    real res = 0.0;

    if (strcasecmp(selected_func->name, "BP86")==0)        res = 4.8516;
    if (strcasecmp(selected_func->name, "BLYP")==0)        res = 4.2359;
    if (strcasecmp(selected_func->name, "revPBE")==0)      res = 3.5016;
    if (strcasecmp(selected_func->name, "B97-D")==0)       res = 3.2297;
    if (strcasecmp(selected_func->name, "PBE")==0)         res = 4.4407;
    if (strcasecmp(selected_func->name, "rpw86-pbe")==0)   res = 4.5062;
    if (strcasecmp(selected_func->name, "B3LYP")==0)       res = 4.4211;
    if (strcasecmp(selected_func->name, "TPSS")==0)        res = 4.4752;
    if (strcasecmp(selected_func->name, "HF")==0)          res = 2.8830;
    if (strcasecmp(selected_func->name, "TPSS0")==0)       res = 4.5865;
    if (strcasecmp(selected_func->name, "PBE0")==0)        res = 4.8593;
    if (strcasecmp(selected_func->name, "revPBE38")==0)    res = 3.9446;
    if (strcasecmp(selected_func->name, "PW6B95")==0)      res = 6.3750;
    if (strcasecmp(selected_func->name, "B2PLYP")==0)      res = 5.0570;
    if (strcasecmp(selected_func->name, "DSD-BLYP")==0)    res = 6.0519;
    if (strcasecmp(selected_func->name, "DSD-BLYP-FC")==0) res = 5.9807;
    if (strcasecmp(selected_func->name, "BOP")==0)         res = 3.5043;
    if (strcasecmp(selected_func->name, "MPWLYP")==0)      res = 4.5323;
    if (strcasecmp(selected_func->name, "OLYP")==0)        res = 2.8065;
    if (strcasecmp(selected_func->name, "PBEsol")==0)      res = 6.1742;
    if (strcasecmp(selected_func->name, "BPBE")==0)        res = 4.3908;
    if (strcasecmp(selected_func->name, "OPBE")==0)        res = 2.9444;
    if (strcasecmp(selected_func->name, "SSB")==0)         res = 5.2170;
    if (strcasecmp(selected_func->name, "revSSB")==0)      res = 4.0986;
    if (strcasecmp(selected_func->name, "OTPSS")==0)       res = 4.3153;
    if (strcasecmp(selected_func->name, "B3PW91")==0)      res = 4.4693;
    if (strcasecmp(selected_func->name, "BHLYP")==0)       res = 4.9615;
    if (strcasecmp(selected_func->name, "revPBE0")==0)     res = 3.7619;
    if (strcasecmp(selected_func->name, "TPSSh")==0)       res = 4.6550;
    if (strcasecmp(selected_func->name, "MPW1B95")==0)     res = 6.4177;
    if (strcasecmp(selected_func->name, "PWB6K")==0)       res = 7.7627;
    if (strcasecmp(selected_func->name, "B1B95")==0)       res = 5.5545;
    if (strcasecmp(selected_func->name, "BMK")==0)         res = 5.9197;
    if (strcasecmp(selected_func->name, "CAM-B3LYP")==0)   res = 5.4743;
    if (strcasecmp(selected_func->name, "LC-wPBE")==0)     res = 5.0987;
    if (strcasecmp(selected_func->name, "B2GP-PLYP")==0)   res = 6.3332;
    if (strcasecmp(selected_func->name, "PTPSS")==0)       res = 6.5745;
    if (strcasecmp(selected_func->name, "PWPB95")==0)      res = 7.3141;
    if (strcasecmp(selected_func->name, "HCTH120")==0)     res = 4.3359;
    if (strcasecmp(selected_func->name, "DFTB3")==0)       res = 4.1906;

    if (res == 0) {
      dalton_quit("DFT-D3-BJ rs_18 not defined for the chosen functional");
    }
    return res;
}


/* give the rs6 factor for the functional, needed in the DFT-D3-BJ corr. */
real
dft_d3bj_rs6_(void)
{
    real res = 0.0;

    if (strcasecmp(selected_func->name, "BP86")==0)        res = 0.3946;
    if (strcasecmp(selected_func->name, "BLYP")==0)        res = 0.4298;
    if (strcasecmp(selected_func->name, "revPBE")==0)      res = 0.5238;
    if (strcasecmp(selected_func->name, "B97-D")==0)       res = 0.5545;
    if (strcasecmp(selected_func->name, "PBE")==0)         res = 0.4289;
    if (strcasecmp(selected_func->name, "rpw86-pbe")==0)   res = 0.4613;
    if (strcasecmp(selected_func->name, "B3LYP")==0)       res = 0.3981;
    if (strcasecmp(selected_func->name, "TPSS")==0)        res = 0.4535;
    if (strcasecmp(selected_func->name, "HF")==0)          res = 0.3385;
    if (strcasecmp(selected_func->name, "TPSS0")==0)       res = 0.3768;
    if (strcasecmp(selected_func->name, "PBE0")==0)        res = 0.4145;
    if (strcasecmp(selected_func->name, "revPBE38")==0)    res = 0.4309;
    if (strcasecmp(selected_func->name, "PW6B95")==0)      res = 0.2076;
    if (strcasecmp(selected_func->name, "B2PLYP")==0)      res = 0.3065;
    if (strcasecmp(selected_func->name, "DSD-BLYP")==0)    res = 0.0000;
    if (strcasecmp(selected_func->name, "DSD-BLYP-FC")==0) res = 0.0009;
    if (strcasecmp(selected_func->name, "BOP")==0)         res = 0.4870;
    if (strcasecmp(selected_func->name, "MPWLYP")==0)      res = 0.4831;
    if (strcasecmp(selected_func->name, "OLYP")==0)        res = 0.5299;
    if (strcasecmp(selected_func->name, "PBEsol")==0)      res = 0.4466;
    if (strcasecmp(selected_func->name, "BPBE")==0)        res = 0.4567;
    if (strcasecmp(selected_func->name, "OPBE")==0)        res = 0.5512;
    if (strcasecmp(selected_func->name, "SSB")==0)         res = -0.0952;
    if (strcasecmp(selected_func->name, "revSSB")==0)      res = 0.4720;
    if (strcasecmp(selected_func->name, "OTPSS")==0)       res = 0.4634;
    if (strcasecmp(selected_func->name, "B3PW91")==0)      res = 0.4312;
    if (strcasecmp(selected_func->name, "BHLYP")==0)       res = 0.2793;
    if (strcasecmp(selected_func->name, "revPBE0")==0)     res = 0.4679;
    if (strcasecmp(selected_func->name, "TPSSh")==0)       res = 0.4529;
    if (strcasecmp(selected_func->name, "MPW1B95")==0)     res = 0.1955;
    if (strcasecmp(selected_func->name, "PWB6K")==0)       res = 0.1805;
    if (strcasecmp(selected_func->name, "B1B95")==0)       res = 0.2092;
    if (strcasecmp(selected_func->name, "BMK")==0)         res = 0.1940;
    if (strcasecmp(selected_func->name, "CAM-B3LYP")==0)   res = 0.3708;
    if (strcasecmp(selected_func->name, "LC-wPBE")==0)     res = 0.3919;
    if (strcasecmp(selected_func->name, "B2GP-PLYP")==0)   res = 0.0000;
    if (strcasecmp(selected_func->name, "PTPSS")==0)       res = 0.0000;
    if (strcasecmp(selected_func->name, "PWPB95")==0)      res = 0.0000;
    if (strcasecmp(selected_func->name, "HCTH120")==0)     res = 0.3563;
    if (strcasecmp(selected_func->name, "DFTB3")==0)       res = 0.7461;

    if (res == 0) {
      dalton_quit("DFT-D3-BJ rs_6 not defined for the chosen functional");
    }
    return res;
}


/* give the s18 factor for the functional, needed in the DFT-D3-BJ corr. */
real
dft_d3bj_s18_(void)
{
    real res = 0.0;

    if (strcasecmp(selected_func->name, "BP86")==0)        res = 3.2822;
    if (strcasecmp(selected_func->name, "BLYP")==0)        res = 2.6996;
    if (strcasecmp(selected_func->name, "revPBE")==0)      res = 2.3550;
    if (strcasecmp(selected_func->name, "B97-D")==0)       res = 2.2609;
    if (strcasecmp(selected_func->name, "PBE")==0)         res = 0.7875;
    if (strcasecmp(selected_func->name, "rpw86-pbe")==0)   res = 1.3845;
    if (strcasecmp(selected_func->name, "B3LYP")==0)       res = 1.9889;
    if (strcasecmp(selected_func->name, "TPSS")==0)        res = 1.9435;
    if (strcasecmp(selected_func->name, "HF")==0)          res = 0.9171;
    if (strcasecmp(selected_func->name, "TPSS0")==0)       res = 1.2576;
    if (strcasecmp(selected_func->name, "PBE0")==0)        res = 1.2177;
    if (strcasecmp(selected_func->name, "revPBE38")==0)    res = 1.4760;
    if (strcasecmp(selected_func->name, "PW6B95")==0)      res = 0.7257;
    if (strcasecmp(selected_func->name, "B2PLYP")==0)      res = 0.9147;
    if (strcasecmp(selected_func->name, "DSD-BLYP")==0)    res = 0.2130;
    if (strcasecmp(selected_func->name, "DSD-BLYP-FC")==0) res = 0.2112;
    if (strcasecmp(selected_func->name, "BOP")==0)         res = 3.2950;
    if (strcasecmp(selected_func->name, "MPWLYP")==0)      res = 2.0077;
    if (strcasecmp(selected_func->name, "OLYP")==0)        res = 2.6205;
    if (strcasecmp(selected_func->name, "PBEsol")==0)      res = 2.9491;
    if (strcasecmp(selected_func->name, "BPBE")==0)        res = 4.0728;
    if (strcasecmp(selected_func->name, "OPBE")==0)        res = 3.3816;
    if (strcasecmp(selected_func->name, "SSB")==0)         res = -0.1744;
    if (strcasecmp(selected_func->name, "revSSB")==0)      res = 0.4389;
    if (strcasecmp(selected_func->name, "OTPSS")==0)       res = 2.7465;
    if (strcasecmp(selected_func->name, "B3PW91")==0)      res = 2.8524;
    if (strcasecmp(selected_func->name, "BHLYP")==0)       res = 1.0354;
    if (strcasecmp(selected_func->name, "revPBE0")==0)     res = 1.7588;
    if (strcasecmp(selected_func->name, "TPSSh")==0)       res = 2.2382;
    if (strcasecmp(selected_func->name, "MPW1B95")==0)     res = 1.0508;
    if (strcasecmp(selected_func->name, "PWB6K")==0)       res = 0.9383;
    if (strcasecmp(selected_func->name, "B1B95")==0)       res = 1.4507;
    if (strcasecmp(selected_func->name, "BMK")==0)         res = 2.0860;
    if (strcasecmp(selected_func->name, "CAM-B3LYP")==0)   res = 2.0674;
    if (strcasecmp(selected_func->name, "LC-wPBE")==0)     res = 1.8541;
    if (strcasecmp(selected_func->name, "B2GP-PLYP")==0)   res = 0.2597;
    if (strcasecmp(selected_func->name, "PTPSS")==0)       res = 0.2804;
    if (strcasecmp(selected_func->name, "PWPB95")==0)      res = 0.2904;
    if (strcasecmp(selected_func->name, "HCTH120")==0)     res = 1.0821;
    if (strcasecmp(selected_func->name, "DFTB3")==0)       res = 3.2090;

    if (res == 0) {
      dalton_quit("DFT-D3-BJ s_18 not defined for the chosen functional");
    }
    return res;
}
/* DFT-D3-BJ Functional Dependent Parameter End */







