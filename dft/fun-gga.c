/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* fun-gga.c:
   implementation of a functional being a linear combination of 
   Slater, VWN, Becke, LYP functionals.
   (c) Pawel Salek, pawsa@theochem.kth.se, sep 2001
   NOTE:
   this file may seem unnecessarily complex but the structure really pays off
   when implementing multiple functionals depending on different parameters.
*/

/* strictly conform to XOPEN ANSI C standard */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1

/* Use BSD's strncasecmp() */
#define _BSD_SOURCE 1

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int  fun_false(void) { return 0; }
static int  fun_true (void) { return 1; }
static int  lda_read(const char* conf_line);
static real lda_energy(const DftDensProp* dp);
static void lda_first(FirstFuncDrv *ds,   real factor, const DftDensProp* dp);
static void lda_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp);
static void lda_third(ThirdFuncDrv *ds,   real factor, const DftDensProp* dp);
static int  ldagauss_read(const char* conf_line);
static int  blyp_read(const char* conf_line);
static int  b3lyp_read(const char* conf_line);
static int  b3lypgauss_read(const char* conf_line);
static int  bp86_read(const char* conf_line);
static int  b3p86_read(const char* conf_line);
static int  kt1_read(const char* conf_line);
static int  kt2_read(const char* conf_line);
static int  kt3_read(const char* conf_line);
static int  olyp_read(const char* conf_line);
static int  pw91_read(const char* conf_line);

static int  gga_isgga(void);
static int  xalpha_read(const char* conf_line);
static int  gga_key_read(const char* conf_line);
static void gga_report(void);
static real gga_energy(const DftDensProp* dp);
static void gga_first(FirstFuncDrv *ds,   real factor, const DftDensProp* dp);
static void gga_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp);
static void gga_third(ThirdFuncDrv *ds,   real factor, const DftDensProp* dp);


Functional XAlphaFunctional = {
    "XAlpha",       /* name */
    fun_false,   /* gga-corrected */
    xalpha_read, 
    NULL,
    gga_energy, 
    gga_first,
    gga_second,
    gga_third
};

Functional LDAFunctional = {
    "LDA",         /* name */
    fun_false,     /* not gga-corrected */
    lda_read,
    NULL,
    lda_energy, 
    lda_first,
    lda_second,
    lda_third
};

Functional LDAGaussFunctional = {
    "LDAGauss",    /* name */
    fun_false,     /* not gga-corrected */
    ldagauss_read,
    NULL,
    gga_energy, 
    gga_first,
    gga_second,
    gga_third
};

Functional GGAKeyFunctional = {
    "GGAKey",      /* name */
    gga_isgga,     /* gga-corrected */
    gga_key_read,
    gga_report,
    gga_energy, 
    gga_first,
    gga_second,
    gga_third
};

Functional BLYPFunctional = {
    "BLYP",      /* name */
    gga_isgga,     /* gga-corrected */
    blyp_read,
    gga_report,
    gga_energy, 
    gga_first,
    gga_second,
    gga_third
};


Functional B3LYPFunctional = {
    "B3LYP",      /* name */
    gga_isgga,     /* gga-corrected */
    b3lyp_read,
    gga_report,
    gga_energy, 
    gga_first,
    gga_second,
    gga_third
};

Functional B3LYPGaussFunctional = {
    "B3LYPGauss",      /* name */
    gga_isgga,     /* gga-corrected */
    b3lypgauss_read,
    gga_report,
    gga_energy, 
    gga_first,
    gga_second,
    gga_third
};


Functional BP86Functional = {
    "BP86",      /* name */
    gga_isgga,     /* gga-corrected */
    bp86_read,
    gga_report,
    gga_energy, 
    gga_first,
    gga_second,
    gga_third
};

Functional B3P86Functional = {
    "B3P86",      /* name */
    gga_isgga,     /* gga-corrected */
    b3p86_read,
    gga_report,
    gga_energy, 
    gga_first,
    gga_second,
    gga_third
};

Functional KT1Functional = {
    "KT1",      /* name */
    gga_isgga,     /* gga-corrected */
    kt1_read,
    gga_report,
    gga_energy,
    gga_first,
    gga_second,
    gga_third
};

Functional KT2Functional = {
    "KT2",      /* name */
    gga_isgga,     /* gga-corrected */
    kt2_read,
    gga_report,
    gga_energy,
    gga_first,
    gga_second,
    gga_third
};

Functional KT3Functional = {
    "KT3",      /* name */
     gga_isgga,     /* gga-corrected */
     kt3_read,
     gga_report,
     gga_energy,
     gga_first,
     gga_second,
     gga_third
};

Functional OLYPFunctional = {
    "OLYP",      /* name */
    gga_isgga,     /* gga-corrected */
    olyp_read,
    gga_report,
    gga_energy,
    gga_first,
    gga_second,
    gga_third
};

Functional PW91Functional = {
    "PW91",      /* name */
    gga_isgga,     /* gga-corrected */
    pw91_read,
    gga_report,
    gga_energy,
    gga_first,
    gga_second,
    gga_third
};


/* MIXED FUNCTIONALS */
typedef struct FuncList_ FuncList;

struct FuncList_ {
    Functional* func;
    real        weight;
    FuncList*   next;
};

FuncList* gga_fun_list = NULL;

static int
gga_isgga(void)
{ 
    int res = 0;

    FuncList* lst;
    for(lst=gga_fun_list; lst && !res; lst=lst->next) 
        res |= lst->func->is_gga();
    return res;
}
static FuncList*
add_functional(FuncList* lst, Functional* f, float weight)
{
    FuncList* n = malloc(sizeof(FuncList));
    n->func = f; n->weight = weight; n->next = lst;
    return n;
}


static int
xalpha_read(const char* conf_line)
{
    float weight;
    int res = (sscanf(conf_line, "%g", &weight)==1);
    if(res) 
        gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 
                                      1.5*weight);
    dft_set_hf_weight(0);
    return res;
}

static int
lda_read(const char* conf_line)
{
    dft_set_hf_weight(0);
    return 1;
}

static real
lda_energy(const DftDensProp* dp)
{
    return SlaterFunctional.func(dp) + VWNFunctional.func(dp);
}

static void
lda_first(FirstFuncDrv *ds, real factor,  const DftDensProp* dp)
{
    SlaterFunctional.first(ds, factor, dp);
    VWNFunctional  .first(ds, factor, dp);
}

static void
lda_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp)
{
    SlaterFunctional.second(ds, factor, dp);
    VWNFunctional  .second(ds, factor, dp);
}

static void
lda_third(ThirdFuncDrv *ds, real factor, const DftDensProp* dp)
{
    SlaterFunctional.third(ds, factor, dp);
    VWNFunctional  .third(ds, factor, dp);
}

static int
ldagauss_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &VWN3Functional,  1.0);
    dft_set_hf_weight(0);
    return 1;
}

static int
blyp_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,   1.0);
    dft_set_hf_weight(0);
    return 1;
}

static int
b3lyp_read(const char* conf_line)
{
    static const real lypw = 0.81, dirw = 0.8;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,  0.72);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,    lypw);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,    1-lypw);
    dft_set_hf_weight(1-dirw);
    return 1;
}

static int
b3lypgauss_read(const char* conf_line)
{
    static const real lypw = 0.81, dirw = 0.8;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,  0.72);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,    lypw);
    gga_fun_list = add_functional(gga_fun_list, &VWN3Functional,   1-lypw);
    dft_set_hf_weight(1-dirw);
    return 1;
}

static int
bp86_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,  1.0);
    dft_set_hf_weight(0);
    return 1;
}

static int
b3p86_read(const char* conf_line)
{
    static const real lypw = 0.81, dirw = 0.8;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,  0.72);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,   lypw);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,    1-lypw);
    dft_set_hf_weight(1-dirw);
    return 1;
}

static int
kt1_read(const char* conf_line)
{
    static const real ktgam = -0.006;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &KTFunctional,    ktgam);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,   1.0);
    dft_set_hf_weight(0);
    return 1;
}

static int
kt2_read(const char* conf_line)
{   
    static const real dirw = 1.07173, vwnw = 0.576727; 
    static const real ktgam = -0.006;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, dirw);
    gga_fun_list = add_functional(gga_fun_list, &KTFunctional,    ktgam);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,   vwnw);
    dft_set_hf_weight(0);
    return 1;
}

static int
kt3_read(const char* conf_line)
{
    static const real dirw = 1.092, lypw = 0.864409, optw = -0.925452;
    static const real ktgam = -0.004;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, dirw);
    gga_fun_list = add_functional(gga_fun_list, &KTFunctional,    ktgam);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,   lypw);
    gga_fun_list = add_functional(gga_fun_list, &OPTXFunctional, optw);
    dft_set_hf_weight(0);
    return 1;
}

static int
olyp_read(const char* conf_line)
{
    static const real optkw = -1.43169, dirw = 1.05151;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, dirw);
    gga_fun_list = add_functional(gga_fun_list, &OPTXFunctional, optkw);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,   1.0);
    dft_set_hf_weight(0);
    return 1;
}

static int
pw91_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &PWggaIIcFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &PWggaIIxFunctional, 1.0);
    dft_set_hf_weight(0);
    return 1;
}


/* gga_key_read:
   read general GGA input: 
*/
static int
gga_key_read(const char* conf_line)
{
    int res = 1, i;
    float f;
    const char* str = conf_line;

    dft_set_hf_weight(0);

    while(*str) {
        while(*str && isspace((int)*str)) str++; /* skip whitespace */
        if(*str =='\0') break; /* line ended by whitespace */
        if(strncasecmp("HF=", str, 3)==0) {
            if(sscanf(str+3,"%g", &f) != 1) {
                fort_print("GGAKey: HF not followed by the weight: ", conf_line);
                res = 0;
            } else dft_set_hf_weight(f);
        } else {
            for(i=0; available_functionals[i]; i++) {
                int len = strlen(available_functionals[i]->name);
                if(strncasecmp(available_functionals[i]->name, str, len)==0 &&
                   str[len] == '=') {
                    if(sscanf(str+len+1,"%g", &f) != 1) {
                        fort_print("GGAKey: keyword '%s' not followed by "
                                   "weight: %s", available_functionals[i]->name, 
                                   conf_line);
                        res = 0;
                    } else {
                        gga_fun_list = add_functional(gga_fun_list, 
                                                      available_functionals[i], f);
                        break; /* weight properly read, break the 'for' loop */
                    }
                }
            }  
            if(available_functionals[i] == NULL) {
                fort_print("GGAKey: functional '%s' not recognised: ", str);
                res = 0;
            }
        }
        while(*str && !isspace((int)*str)) str++; /* skip nonws */
    } return res;
}

static void
gga_report(void)
{
    FuncList* lst;
    fort_print("Weighted mixed functional:");
    if(dft_get_hf_weight()>0)
        fort_print("%30s: %10.5f", "HF exchange", dft_get_hf_weight());
    for(lst=gga_fun_list; lst; lst=lst->next) 
        fort_print("%30s: %10.5f", lst->func->name, lst->weight);
}

static real
gga_energy(const DftDensProp* dp)
{
    real res = 0;
    FuncList* lst;
    for(lst=gga_fun_list; lst; lst=lst->next) {
        real contr = lst->weight*lst->func->func(dp);
/*        fort_print("[%g,%g,w=%g] %s contributes with %g", dp->rhoa,
                   dp->grada, lst->weight, lst->func->name, contr); */
        res += contr;
    }
    return res;
}

static void
gga_first(FirstFuncDrv *ds, real factor,  const DftDensProp* dp)
{
    FuncList* lst;
    for(lst=gga_fun_list; lst; lst=lst->next) {
	real df10 = ds->df1000, df01 = ds->df0010;
        lst->func->first(ds, factor*lst->weight, dp);
/*      fort_print("[%g,%g,w=%g] %s f: %g deriv (%g,%g)", dp->rhoa,
                   dp->grada, lst->weight, lst->func->name, factor,
		   ds->df1000-df10, ds->df0010-df01); */
    }
}

static void
gga_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp)
{
    FuncList* lst;
    for(lst=gga_fun_list; lst; lst=lst->next) 
        lst->func->second(ds, factor*lst->weight, dp);
}

static void
gga_third(ThirdFuncDrv *ds, real factor, const DftDensProp* dp)
{
    FuncList* lst;
    for(lst=gga_fun_list; lst; lst=lst->next) 
        lst->func->third(ds, factor*lst->weight, dp);
}


