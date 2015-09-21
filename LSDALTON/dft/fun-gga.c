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
#include <strings.h>
#include <stdlib.h>
#include "lsdalton_general.h"

#define __CVERSION__

#include "lsdalton_functionals.h"

/* INTERFACE PART */
static integer  lda_read(const char* conf_line, real *hfweight);
static real lda_energy(const FunDensProp* dp);
static void lda_first(FunFirstFuncDrv *ds,   real fac, const FunDensProp* dp);
static void lda_second(FunSecondFuncDrv *ds, real fac, const FunDensProp* dp);
static void lda_third(FunThirdFuncDrv *ds,   real fac, const FunDensProp* dp);
static void lda_fourth(FunFourthFuncDrv *ds, real fac, const FunDensProp* dp);
static real ldax_energy(const FunDensProp* dp);
static void ldax_first(FunFirstFuncDrv *ds,   real fac, const FunDensProp* dp);
static void ldax_second(FunSecondFuncDrv *ds, real fac, const FunDensProp* dp);
static void ldax_third(FunThirdFuncDrv *ds,   real fac, const FunDensProp* dp);
static void ldax_fourth(FunFourthFuncDrv *ds, real fac, const FunDensProp* dp);
static integer  ldagauss_read(const char* conf_line, real *hfweight);
static integer  blyp_read(const char* conf_line, real *hfweight);
static integer  b88x_read(const char* conf_line, real *hfweight);
static integer  ldax_read(const char* conf_line, real *hfweight);
static integer  pbex_read(const char* conf_line, real *hfweight);
static integer  revpbex_read(const char* conf_line, real *hfweight);
static integer  rpbex_read(const char* conf_line, real *hfweight);
static integer  mpbex_read(const char* conf_line, real *hfweight);
static integer  pw91x_read(const char* conf_line, real *hfweight);
static integer  kt1x_read(const char* conf_line, real *hfweight);
static integer  kt2x_read(const char* conf_line, real *hfweight);
static integer  kt3x_read(const char* conf_line, real *hfweight);
static integer  g96x_read(const char* conf_line, real *hfweight);
static integer  lg93x_read(const char* conf_line, real *hfweight);
static integer  optx_read(const char* conf_line, real *hfweight);
static integer  b3lyp_read(const char* conf_line, real *hfweight);
static integer  b3lypgauss_read(const char* conf_line, real *hfweight);
static integer  bp86_read(const char* conf_line, real *hfweight);
static integer  bpw91_read(const char* conf_line, real *hfweight);
static integer  b3p86_read(const char* conf_line, real *hfweight);
static integer  b3p86g_read(const char* conf_line, real *hfweight);
static integer  kt1_read(const char* conf_line, real *hfweight);
static integer  kt2_read(const char* conf_line, real *hfweight);
static integer  kt3_read(const char* conf_line, real *hfweight);
static integer  olyp_read(const char* conf_line, real *hfweight);
static integer  pbe_read(const char* conf_line, real *hfweight);
static integer  pbe0_read(const char* conf_line, real *hfweight);

static integer  gga_isgga(void);
static integer  xalpha_read(const char* conf_line, real *hfweight);
static integer  gga_key_read(const char* conf_line, real *hfweight);
static void gga_report(integer lupri);
static real gga_energy(const FunDensProp* dp);
static void gga_first(FunFirstFuncDrv *ds,   real fac, const FunDensProp* dp);
static void gga_second(FunSecondFuncDrv *ds, real fac, const FunDensProp* dp);
static void gga_third(FunThirdFuncDrv *ds,   real fac, const FunDensProp* dp);
static void gga_fourth(FunFourthFuncDrv *ds, real fac, const FunDensProp* dp);

#define LDA_FUNCTIONAL(name,read) { (name), \
    fun_false, (read), NULL, lda_energy, lda_first, lda_second, \
    lda_third, lda_fourth }

#define LDAX_FUNCTIONAL(name,read) { (name), \
    fun_false, (read), NULL, ldax_energy, ldax_first, ldax_second, \
    ldax_third, ldax_fourth }

#define GGA_FUNCTIONAL(name,read) { (name), \
    gga_isgga, (read), gga_report, gga_energy, gga_first, gga_second, \
    gga_third, gga_fourth }

Functional XAlphaFunctional     = GGA_FUNCTIONAL("XAlpha",  xalpha_read);
Functional LDAFunctional        = LDA_FUNCTIONAL("LDA",     lda_read);
/* SVWN5 aliases LDA */
Functional SVWN5Functional      = LDA_FUNCTIONAL("SVWN5",   lda_read);
Functional SVWN3Functional      = GGA_FUNCTIONAL("SVWN3",   ldagauss_read);
Functional B3LYPFunctional      = GGA_FUNCTIONAL("B3LYP",   b3lyp_read);
Functional B3LYPGaussFunctional = GGA_FUNCTIONAL("B3LYP-G", b3lypgauss_read);
Functional B3P86Functional      = GGA_FUNCTIONAL("B3P86",   b3p86_read);
Functional B3P86GFunctional     = GGA_FUNCTIONAL("B3P86-G", b3p86g_read);
Functional BLYPFunctional       = GGA_FUNCTIONAL("BLYP",    blyp_read);
Functional B88X_Functional      = GGA_FUNCTIONAL("B88X",    b88x_read);
Functional LDAX_Functional      = LDAX_FUNCTIONAL("LDAX",   ldax_read);
Functional PBEX_Functional      = GGA_FUNCTIONAL("PBEX",    pbex_read);
Functional RPBEX_Functional     = GGA_FUNCTIONAL("RPBEX",   rpbex_read);
Functional REVPBEX_Functional   = GGA_FUNCTIONAL("REVPBEX", revpbex_read);
Functional MPBEX_Functional     = GGA_FUNCTIONAL("MPBEX",   mpbex_read);
Functional PW91X_Functional     = GGA_FUNCTIONAL("PW91X",   pw91x_read);
Functional KT1X_Functional      = GGA_FUNCTIONAL("KT1X",    kt1x_read);
Functional KT2X_Functional      = GGA_FUNCTIONAL("KT2X",    kt2x_read);
Functional KT3X_Functional      = GGA_FUNCTIONAL("KT3X",    kt3x_read);
Functional G96X_Functional      = GGA_FUNCTIONAL("G96X",    g96x_read);
Functional LG93X_Functional     = GGA_FUNCTIONAL("LG93X",   lg93x_read);
Functional OPTX_Functional      = GGA_FUNCTIONAL("OPTX",    optx_read);
Functional BP86Functional       = GGA_FUNCTIONAL("BP86",    bp86_read);
Functional BPW91Functional      = GGA_FUNCTIONAL("BPW91",   bpw91_read);
Functional GGAKeyFunctional     = GGA_FUNCTIONAL("GGAKey",  gga_key_read);
Functional KT1Functional        = GGA_FUNCTIONAL("KT1",     kt1_read);
Functional KT2Functional        = GGA_FUNCTIONAL("KT2",     kt2_read);
Functional KT3Functional        = GGA_FUNCTIONAL("KT3",     kt3_read);
Functional OLYPFunctional       = GGA_FUNCTIONAL("OLYP" ,   olyp_read);
Functional PBE0Functional       = GGA_FUNCTIONAL("PBE0",    pbe0_read);
Functional PBEFunctional        = GGA_FUNCTIONAL("PBE",     pbe_read); 

/* MIXED FUNCTIONALS */
typedef struct FuncList_ FuncList;

struct FuncList_ {
    Functional* func;
    real        weight;
    FuncList*   next;
};

FuncList* gga_fun_list = NULL;

static integer
gga_isgga(void)
{ 
    integer res = 0;

    FuncList* lst;
    for(lst=gga_fun_list; lst && !res; lst=lst->next) 
        res |= lst->func->is_gga();
    return res;
}
static FuncList*
add_functional(FuncList* lst, Functional* f, real weight)
{
    FuncList* n = malloc(sizeof(FuncList));
    n->func = f; n->weight = weight; n->next = lst;
    return n;
}
integer 
clear_funclist()
{
    FuncList* lst;
    FuncList* previous;
    previous = NULL;
    for(lst=gga_fun_list; lst; lst=lst->next)  {
      if (previous) free(previous);
      previous = lst;
    }
    if (previous) free(previous);
    gga_fun_list = NULL;
    return 1;
}

static integer
xalpha_read(const char* conf_line, real *hfweight)
{
    real weight;
    integer res = (sscanf(conf_line, "%lf", &weight)==1);
    if(res) 
        gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 
                                      1.5*weight);
    return res;
}

static integer
lda_read(const char* conf_line, real *hfweight)
{
    return 1;
}

static real
lda_energy(const FunDensProp* dp)
{
    return SlaterFunctional.func(dp) + VWNFunctional.func(dp);
}

static void
lda_first(FunFirstFuncDrv *ds, real factor,  const FunDensProp* dp)
{
    SlaterFunctional.first(ds, factor, dp);
    VWNFunctional  .first(ds, factor, dp);
}

static void
lda_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    SlaterFunctional.second(ds, factor, dp);
    VWNFunctional  .second(ds, factor, dp);
}

static void
lda_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    SlaterFunctional.third(ds, factor, dp);
    VWNFunctional  .third(ds, factor, dp);
}

static void
lda_fourth(FunFourthFuncDrv *ds, real factor, const FunDensProp* dp)
{
    SlaterFunctional.fourth(ds, factor, dp);
    VWNFunctional   .fourth(ds, factor, dp);
}

static integer
ldax_read(const char* conf_line, real *hfweight)
{
    return 1;
}

static real
ldax_energy(const FunDensProp* dp)
{
    return SlaterFunctional.func(dp);
}

static void
ldax_first(FunFirstFuncDrv *ds, real factor,  const FunDensProp* dp)
{
    SlaterFunctional.first(ds, factor, dp);
}

static void
ldax_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    SlaterFunctional.second(ds, factor, dp);
}

static void
ldax_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    SlaterFunctional.third(ds, factor, dp);
}

static void
ldax_fourth(FunFourthFuncDrv *ds, real factor, const FunDensProp* dp)
{
    SlaterFunctional.fourth(ds, factor, dp);
}

static integer
ldagauss_read(const char* conf_line, real *hfweight)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &VWN3Functional,  1.0);
    return 1;
}

static integer
blyp_read(const char* conf_line, real *hfweight)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,   1.0);
    return 1;
}

static integer
b88x_read(const char* conf_line, real *hfweight)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, *hfweight);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional, *hfweight);
    return 1;
}

static integer
pbex_read(const char* conf_line, real *hfweight)
{
    gga_fun_list = add_functional(gga_fun_list, &PbexFunctional, *hfweight);
    return 1;
}

static integer
revpbex_read(const char* conf_line, real *hfweight)
{
    gga_fun_list = add_functional(gga_fun_list, &revPBExFunctional, *hfweight);
    return 1;
}

static integer
rpbex_read(const char* conf_line, real *hfweight)
{
    gga_fun_list = add_functional(gga_fun_list, &rPBExFunctional, *hfweight);
    return 1;
}

static integer
mpbex_read(const char* conf_line, real *hfweight)
{
    gga_fun_list = add_functional(gga_fun_list, &mPBExFunctional, *hfweight);
    return 1;
}

static integer
pw91x_read(const char* conf_line, real *hfweight)
{
    gga_fun_list = add_functional(gga_fun_list, &PW91xFunctional, *hfweight);
    return 1;
}

static integer
kt1x_read(const char* conf_line, real *hfweight)
{
    static const real dirw = 1.0;
    static const real ktgam = -0.006;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, dirw*(*hfweight));
    gga_fun_list = add_functional(gga_fun_list, &KTFunctional,    ktgam*(*hfweight));
    return 1;
}

static integer
kt2x_read(const char* conf_line, real *hfweight)
{
    static const real dirw = 1.07173;
    static const real ktgam = -0.006;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, dirw*(*hfweight));
    gga_fun_list = add_functional(gga_fun_list, &KTFunctional,    ktgam*(*hfweight));
    return 1;
}

static integer
kt3x_read(const char* conf_line, real *hfweight)
{
    static const real dirw = 1.092, optw = -0.925452;
    static const real ktgam = -0.004;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, dirw*(*hfweight));
    gga_fun_list = add_functional(gga_fun_list, &KTFunctional,    ktgam*(*hfweight));
    gga_fun_list = add_functional(gga_fun_list, &OPTXFunctional,   optw*(*hfweight));
    return 1;
}

static integer
g96x_read(const char* conf_line, real *hfweight)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, *hfweight);
    gga_fun_list = add_functional(gga_fun_list, &G96xFunctional,   *hfweight);
    return 1;
}

static integer
lg93x_read(const char* conf_line, real *hfweight)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, *hfweight);
    gga_fun_list = add_functional(gga_fun_list, &LG93xFunctional,  *hfweight);
    return 1;
}

static integer
optx_read(const char* conf_line, real *hfweight)
{
    static const real optkw = -1.43169, dirw = 1.05151;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, dirw*(*hfweight));
    gga_fun_list = add_functional(gga_fun_list, &OPTXFunctional,  optkw*(*hfweight));
    return 1;
}

static integer
b3lyp_read(const char* conf_line, real *hfweight)
{
    static const real lypw = 0.81, dirw = 0.8;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,  0.72);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,    lypw);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,    1-lypw);
    *hfweight = 1-dirw;
    return 1;
}

static integer
b3lypgauss_read(const char* conf_line, real *hfweight)
{
    static const real lypw = 0.81, dirw = 0.8;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,  0.72);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,    lypw);
    gga_fun_list = add_functional(gga_fun_list, &VWN3Functional,   1-lypw);
    *hfweight = 1-dirw;
    return 1;
}

static integer
bp86_read(const char* conf_line, real *hfweight)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &PZ81Functional,  1.0);
    return 1;
}

static integer
b3p86_read(const char* conf_line, real *hfweight)
{
    static const real lypw = 0.81, dirw = 0.8;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,  0.72);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,   lypw);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,    1);
    *hfweight = 1-dirw;
    return 1;
}

static integer
b3p86g_read(const char* conf_line, real *hfweight)
{
    static const real lypw = 0.81, dirw = 0.8;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,  0.72);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,   lypw);
    gga_fun_list = add_functional(gga_fun_list, &VWN3Functional,   1);
    *hfweight = 1-dirw;
    return 1;
}

static integer
bpw91_read(const char* conf_line, real *hfweight)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,  1);
    gga_fun_list = add_functional(gga_fun_list, &Pw91cFunctional,  1);
    return 1;
}

static integer
kt1_read(const char* conf_line, real *hfweight)
{
    static const real ktgam = -0.006;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &KTFunctional,    ktgam);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,   1.0);
    return 1;
}

static integer
kt2_read(const char* conf_line, real *hfweight)
{   
    static const real dirw = 1.07173, vwnw = 0.576727; 
    static const real ktgam = -0.006;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, dirw);
    gga_fun_list = add_functional(gga_fun_list, &KTFunctional,    ktgam);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,   vwnw);
    return 1;
}

static integer
kt3_read(const char* conf_line, real *hfweight)
{
    static const real dirw = 1.092, lypw = 0.864409, optw = -0.925452;
    static const real ktgam = -0.004;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, dirw);
    gga_fun_list = add_functional(gga_fun_list, &KTFunctional,    ktgam);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,   lypw);
    gga_fun_list = add_functional(gga_fun_list, &OPTXFunctional, optw);
    return 1;
}

static integer
olyp_read(const char* conf_line, real *hfweight)
{
    static const real optkw = -1.43169, dirw = 1.05151;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, dirw);
    gga_fun_list = add_functional(gga_fun_list, &OPTXFunctional, optkw);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,   1.0);
    return 1;
}

static integer
pbe_read(const char* conf_line, real *hfweight)
{
    gga_fun_list = add_functional(gga_fun_list, &PBEcFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &PbexFunctional, 1.0);
    return 1;
}
 
static integer
pbe0_read(const char* conf_line, real *hfweight)
{
    gga_fun_list = add_functional(gga_fun_list, &PBEcFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &PbexFunctional, 0.75);
    *hfweight = 0.25;
    return 1;
}


/* gga_key_read:
   read general GGA input: 
*/
static integer
gga_key_read(const char* conf_line, real *hfweight)
{
    integer res = 1, i;
    real f;
    const char* str = conf_line;

    *hfweight = 0.0000;

    while(*str) {
        while(*str && isspace((integer)*str)) str++; /* skip whitespace */
        if(*str =='\0') break; /* line ended by whitespace */
        if(strncasecmp("HF=", str, 3)==0) {
            if(sscanf(str+3,"%lf", &f) != 1) {
                printf("GGAKey: HF not followed by the weight: %s",conf_line);
                res = 0;
            } else {
	      *hfweight = f;
		    }
        } else {
            for(i=0; available_functionals[i]; i++) {
                integer len = strlen(available_functionals[i]->name);
                if(strncasecmp(available_functionals[i]->name, str, len)==0 &&
                   str[len] == '=') {
                    if(sscanf(str+len+1,"%lf", &f) != 1) {
                        printf("GGAKey: keyword '%s' not followed by "
                                   "weight: %s",
                                   available_functionals[i]->name, 
                                   conf_line);
                        res = 0;
                    } else {
                        gga_fun_list = 
                            add_functional(gga_fun_list, 
                                           available_functionals[i], f);
                        break; /* weight properly read, break the 'for' loop */
                    }
                }
            }  
            if(available_functionals[i] == NULL) {
                printf("GGAKey: functional '%s' not recognised: ", str);
                res = 0;
            }
        }
        while(*str && !isspace((integer)*str)) str++; /* skip nonws */
    } return res;
}

static void
gga_report(integer lupri)
{
    FuncList* lst;
    lsfort_print(lupri,"Weighted mixed functional:");
    for(lst=gga_fun_list; lst; lst=lst->next) 
      lsfort_print(lupri,"%25s: %10.5f", lst->func->name, lst->weight);
}

static real
gga_energy(const FunDensProp* dp)
{
    real res = 0;
    FuncList* lst;
    for(lst=gga_fun_list; lst; lst=lst->next) {
        real contr = lst->weight*lst->func->func(dp);
/*        printf("[%g,%g,w=%g] %s contributes with %g", dp->rhoa,
                   dp->grada, lst->weight, lst->func->name, contr); */
        res += contr;
    }
    return res;
}

static void
gga_first(FunFirstFuncDrv *ds, real factor,  const FunDensProp* dp)
{
    FuncList* lst;
    for(lst=gga_fun_list; lst; lst=lst->next) {
        lst->func->first(ds, factor*lst->weight, dp);
/*      printf("[%g,%g,w=%g] %s f: %g deriv (%g,%g)", dp->rhoa,
                   dp->grada, lst->weight, lst->func->name, factor,
		   ds->df1000-df10, ds->df0010-df01); */
    }
}

static void
gga_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    FuncList* lst;
    for(lst=gga_fun_list; lst; lst=lst->next) 
        lst->func->second(ds, factor*lst->weight, dp);
}

static void
gga_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    FuncList* lst;
    for(lst=gga_fun_list; lst; lst=lst->next) 
        lst->func->third(ds, factor*lst->weight, dp);
}


static void
gga_fourth(FunFourthFuncDrv *ds, real factor, const FunDensProp* dp)
{
    FuncList* lst;
    for(lst=gga_fun_list; lst; lst=lst->next) 
        lst->func->fourth(ds, factor*lst->weight, dp);
}


