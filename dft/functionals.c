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
#define __CVERSION__

#include "functionals.h"

Functional* available_functionals[] = {
    /* generic functionals */
    &BeckeFunctional,
    &mBeckeFunctional,
    &B86xFunctional,
    &B86mxFunctional,
    &B97Functional,
    &B97_1Functional,
    &B97_2Functional,
    &B97_3Functional,
    &B97_KFunctional,
    &Example2Functional,
    &ExampleFunctional,
    &DK87xFunctional,
    &G96xFunctional,
    &KTxFunctional,
    &LB94Functional,
    &LG93xFunctional,
    &LRC95xFunctional,
    &LYPFunctional,
    &LYPrFunctional,
    &HCTHFunctional,
    &HCTH93Functional,
    &HCTH93mFunctional,
    &HCTH120Functional,
    &HCTH147Functional,
    &HCTH407Functional,
    &HCTH407pFunctional,
    &OPTXFunctional,
    &mPWxFunctional,
    &P86cFunctional,
    &PW86xFunctional,
    &PW91cFunctional,
    &PW91ncFunctional,
    &PW91xFunctional,
    &PW91x2Functional,
    &PW92cFunctional,
    &PW92acFunctional,
    &PZ81Functional,
    &PBEcFunctional,
    &PBExFunctional,
    &RPBExFunctional,
    &mPBExFunctional,
    &revPBExFunctional,
    &SlaterFunctional,
    &VWN3Functional,
    &VWN5Functional,
    &VWNIFunctional,
    &VWN3IFunctional,
    &VWNFunctional,
    &XAlphaFunctional,
    &WignerFunctional,
    &WL90cFunctional,
    /* mixed functionals */
    &B2PLYPFunctional,
    &B3LYPFunctional,
    &B3LYPgFunctional,
    &B3LYPGaussFunctional,
    &B3P86Functional,
    &B3P86gFunctional,
    &B3PW91Functional,
    &B1LYPFunctional,
    &B1PW91Functional,
    &BHandHFunctional,
    &BHandHLYPFunctional,
    &B86VWNFunctional,
    &B86LYPFunctional,
    &B86P86Functional,
    &B86PW91Functional,
    &BVWNFunctional,
    &BLYPFunctional,
    &BP86Functional,
    &BPW91Functional,
    &BWFunctional,
    &BFWFunctional,
    &Camb3lypFunctional,
    &CombineFunctional,
    &DBLYPFunctional,
    &DBP86Functional,
    &DBPW91Functional,
    &EDF1Functional,
    &EDF2Functional,
    &GGAKeyFunctional,
    &G96VWNFunctional,
    &G96LYPFunctional,
    &G96P86Functional,
    &G96PW91Functional,
    &G961LYPFunctional,
    &KMLYPFunctional,
    &KT1Functional,
    &KT2Functional,
    &KT3Functional,
    &LDAFunctional,
    &LG1LYPFunctional,
    &mPWVWNFunctional,
    &mPWLYPFunctional,
    &mPWP86Functional,
    &mPWPW91Functional,
    &mPW91Functional,
    &mPW1PW91Functional,
    &mPW3PW91Functional,
    &mPW1KFunctional,
    &mPW1NFunctional,
    &mPW1SFunctional,
    &OVWNFunctional,
    &OLYPFunctional,
    &OP86Functional,
    &OPW91Functional,
    &PBE0Functional,
    &PBE0PBEFunctional,
    &PBE1PBEFunctional,
    &PBEFunctional,
    &PBEPBEFunctional,
    &revPBEFunctional,
    &RPBEFunctional,
    &mPBEFunctional,
    &PW91Functional,
    &revPBEFunctional,
    &PW91VWNFunctional,
    &PW91LYPFunctional,
    &PW91P86Functional,
    &PW91PW91Functional,
    &SVWN3Functional,
    &SVWN5Functional,
    &XLYPFunctional,
    &X3LYPFunctional,
    NULL
};

static int my_printf(const char *fmt, ...)
{
    int i;va_list ap; va_start(ap, fmt); i= vprintf(fmt, ap); va_end(ap);
    puts("");
    return i;
}
 
static void set_hf_weight(real w)         {}
static real get_hf_weight(void)           {return 0;}
static void set_cam_param(int cnt, const real *w, const real *b) {}

Functional* selected_func = &LDAFunctional;
int (*fun_printf)(const char *fmt, ...) = my_printf;
void (*fun_set_hf_weight)(real w)         = set_hf_weight;
real (*fun_get_hf_weight)(void)           = get_hf_weight;
void (*fun_set_mp2_weight)(real w)        = set_hf_weight;
real (*fun_get_mp2_weight)(void)          = get_hf_weight;
void (*fun_set_cam_param)(int cnt, const real *mu, const real *b)
     = set_cam_param;

/* =================================================================== */
enum FunError
fun_select_by_name(const char *conf_string)
{
    int ok, i;
    char func_name[20];

    sscanf(conf_string,"%20s", func_name);
    for(i=0; available_functionals[i]; i++)
        if(strcasecmp(available_functionals[i]->name, func_name)==0) {
            selected_func = available_functionals[i];
            ok = selected_func->read ?
                selected_func->read(conf_string+strlen(func_name)) : 1;
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

int fun_true(void)  { return 1; }
int fun_false(void) { return 0; }

/* Fortran interface. We specify different names suffixes so that
 * library can be linked with code compiled with different compilers
 * or different compilation options. */

void
funset(const char *str, int *info, int len)
{
    switch(fun_select_by_name(str)) {
    case FUN_OK:         *info = 0; break;
    case FUN_UNKNOWN   : *info = 1; break;
    case FUN_CONF_ERROR: *info = 2; break;
    }
}
void
funset_(const char *str, int *info, int len)
{
    switch(fun_select_by_name(str)) {
    case FUN_OK:         *info = 0; break;
    case FUN_UNKNOWN   : *info = 1; break;
    case FUN_CONF_ERROR: *info = 2; break;
    }
}

int
funisgga(void)
{ return selected_func->is_gga(); }
int
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
dftreport_(void)
{
    fun_printf("\n     This is a DFT calculation of type: %s",
               selected_func->name);
    if(selected_func->report)
        selected_func->report();
}

void
dftlistfuncs_(void)
{
    int i;
    fun_printf("\nAvailable functionals:");
    for(i=0; available_functionals[i]; i++)
        fun_printf(available_functionals[i]->name);
}

/* declare both known fortran name-mangled variants */
int
dft_isgga_(void)
{ return selected_func->is_gga(); }

int
dft_isgga__(void)
{ return selected_func->is_gga(); }
