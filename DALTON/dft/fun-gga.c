/*-*-mode: 

!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2016 (2015), see http://daltonprogram.org"
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
/* fun-gga.c:
   implementation of a functional being a linear combination of 
   exchange and correlation functionals
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
#include "general.h"

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static integer  lda_read(const char* conf_line);
static real lda_energy(const FunDensProp* dp);
static void lda_first(FunFirstFuncDrv *ds,   real fac, const FunDensProp* dp);
static void lda_second(FunSecondFuncDrv *ds, real fac, const FunDensProp* dp);
static void lda_third(FunThirdFuncDrv *ds,   real fac, const FunDensProp* dp);
static void lda_fourth(FunFourthFuncDrv *ds, real fac, const FunDensProp* dp);
static integer  ldagauss_read(const char* conf_line);
static integer  b86vwn_read(const char* conf_line);
static integer  b86lyp_read(const char* conf_line);
static integer  b86p86_read(const char* conf_line);
static integer  b86pw91_read(const char* conf_line);
static integer  bvwn_read(const char* conf_line);
static integer  blyp_read(const char* conf_line);
static integer  b1lyp_read(const char* conf_line);
static integer  b2plyp_read(const char* conf_line);
static integer  b2tplyp_read(const char* conf_line);
static integer  mpw2plyp_read(const char* conf_line);
static integer  mpw2kplyp_read(const char* conf_line);
static integer  b2gpplyp_read(const char* conf_line);
static integer  b2piplyp_read(const char* conf_line);
static integer  pbe0dh_read(const char* conf_line);
static integer  b3lyp_read(const char* conf_line);
static integer  b3lypg_read(const char* conf_line);
static integer  bp86_read(const char* conf_line);
static integer  b3p86_read(const char* conf_line);
static integer  b3p86g_read(const char* conf_line);
static integer  bpw91_read(const char* conf_line);
static integer  b1pw91_read(const char* conf_line);
static integer  b3pw91_read(const char* conf_line);
static integer  bhandh_read(const char* conf_line);
static integer  bhandhlyp_read(const char* conf_line);
static integer  bw_read(const char* conf_line);
static integer  bfw_read(const char* conf_line);
static integer  dblyp_read(const char* conf_line);
static integer  dbp86_read(const char* conf_line);
static integer  dbpw91_read(const char* conf_line);
static integer  edf1_read(const char* conf_line);
static integer  edf2_read(const char* conf_line);
static integer  g96vwn_read(const char* conf_line);
static integer  g96lyp_read(const char* conf_line);
static integer  g96p86_read(const char* conf_line);
static integer  g96pw91_read(const char* conf_line);
static integer  g961lyp_read(const char* conf_line);
static integer  hcth_read(const char* conf_line);
static integer  kmlyp_read(const char* conf_line);
static integer  kt1_read(const char* conf_line);
static integer  kt2_read(const char* conf_line);
static integer  kt3_read(const char* conf_line);
static integer  lg1lyp_read(const char* conf_line);
static integer  mpwvwn_read(const char* conf_line);
static integer  mpwlyp_read(const char* conf_line);
static integer  mpwp86_read(const char* conf_line);
static integer  mpwpw91_read(const char* conf_line);
static integer  mpw3pw91_read(const char* conf_line);
static integer  mpw1pw91_read(const char* conf_line);
static integer  mpw1k_read(const char* conf_line);
static integer  mpw1n_read(const char* conf_line);
static integer  mpw1s_read(const char* conf_line);
static integer  ovwn_read(const char* conf_line);
static integer  olyp_read(const char* conf_line);
static integer  op86_read(const char* conf_line);
static integer  opw91_read(const char* conf_line);
static integer  pbe_read(const char* conf_line);
static integer  pbe0_read(const char* conf_line);
static integer  revpbe_read(const char* conf_line);
static integer  rpbe_read(const char* conf_line);
static integer  mpbe_read(const char* conf_line);
static integer  pw91_read(const char* conf_line);
static integer  pw91vwn_read(const char* conf_line);
static integer  pw91lyp_read(const char* conf_line);
static integer  pw91p86_read(const char* conf_line);
static integer  xlyp_read(const char* conf_line);
static integer  x3lyp_read(const char* conf_line);

static integer  gga_isgga(void);
static integer  xalpha_read(const char* conf_line);
static integer  gga_key_read(const char* conf_line);
static void gga_report(void);
static real gga_energy(const FunDensProp* dp);
static void gga_first(FunFirstFuncDrv *ds,   real fac, const FunDensProp* dp);
static void gga_second(FunSecondFuncDrv *ds, real fac, const FunDensProp* dp);
static void gga_third(FunThirdFuncDrv *ds,   real fac, const FunDensProp* dp);
static void gga_fourth(FunFourthFuncDrv *ds, real fac, const FunDensProp* dp);

#define LDA_FUNCTIONAL(name,resp_order,read) { (name), \
    fun_false, (resp_order), (read), NULL, lda_energy, lda_first, lda_second, \
    lda_third, lda_fourth }

#define GGA_FUNCTIONAL(name,resp_order,read) { (name), \
    gga_isgga, (resp_order), (read), gga_report, gga_energy, gga_first, gga_second, \
    gga_third, gga_fourth }

Functional B2TPLYPFunctional    = GGA_FUNCTIONAL("B2TPLYP",    1, b2tplyp_read);
Functional MPW2PLYPFunctional   = GGA_FUNCTIONAL("mPW2PLYP",   1, mpw2plyp_read);
Functional MPW2KPLYPFunctional  = GGA_FUNCTIONAL("mPW2KPLYP",  1, mpw2kplyp_read);
Functional B2GPPLYPFunctional   = GGA_FUNCTIONAL("B2GPPLYP",   1, b2gpplyp_read);
Functional B2PIPLYPFunctional   = GGA_FUNCTIONAL("B2PIPLYP",   1, b2piplyp_read);
Functional PBE0DHFunctional     = GGA_FUNCTIONAL("PBE0DH",     1, pbe0dh_read);
Functional XAlphaFunctional     = GGA_FUNCTIONAL("XAlpha",     3, xalpha_read);
Functional LDAFunctional        = LDA_FUNCTIONAL("LDA",        3, lda_read);
Functional SVWN5Functional      = LDA_FUNCTIONAL("SVWN5",      3, lda_read);       // SVWN5 aliases LDA
Functional SVWN3Functional      = GGA_FUNCTIONAL("SVWN3",      1, ldagauss_read);
Functional B2PLYPFunctional     = GGA_FUNCTIONAL("B2PLYP",     3, b2plyp_read);
Functional B3LYPFunctional      = GGA_FUNCTIONAL("B3LYP",      3, b3lyp_read);
Functional B3LYPgFunctional     = GGA_FUNCTIONAL("B3LYPg",     1, b3lypg_read);
Functional B3LYPGaussFunctional = GGA_FUNCTIONAL("B3LYPGauss", 1, b3lypg_read);
Functional B3P86Functional      = GGA_FUNCTIONAL("B3P86",      1, b3p86_read);
Functional B3P86gFunctional     = GGA_FUNCTIONAL("B3P86g",     1, b3p86g_read);
Functional B3PW91Functional     = GGA_FUNCTIONAL("B3PW91",     1, b3pw91_read);
Functional BHandHFunctional     = GGA_FUNCTIONAL("BHandH",     3, bhandh_read);
Functional BHandHLYPFunctional  = GGA_FUNCTIONAL("BHandHLYP",  3, bhandhlyp_read);
Functional B86VWNFunctional     = GGA_FUNCTIONAL("B86VWN",     1, b86vwn_read);
Functional B86LYPFunctional     = GGA_FUNCTIONAL("B86LYP",     1, b86lyp_read);
Functional B86P86Functional     = GGA_FUNCTIONAL("B86P86",     1, b86p86_read);
Functional B86PW91Functional    = GGA_FUNCTIONAL("B86PW91",    1, b86pw91_read);
Functional BVWNFunctional       = GGA_FUNCTIONAL("BVWN",       3, bvwn_read);
Functional BLYPFunctional       = GGA_FUNCTIONAL("BLYP",       3, blyp_read);
Functional B1LYPFunctional      = GGA_FUNCTIONAL("B1LYP",      3, b1lyp_read);
Functional BP86Functional       = GGA_FUNCTIONAL("BP86",       1, bp86_read);
Functional BPW91Functional      = GGA_FUNCTIONAL("BPW91",      1, bpw91_read);
Functional B1PW91Functional     = GGA_FUNCTIONAL("B1PW91",     1, b1pw91_read);
Functional BWFunctional         = GGA_FUNCTIONAL("BW",         1, bw_read);
Functional BFWFunctional        = GGA_FUNCTIONAL("BFW",        1, bfw_read);
Functional CombineFunctional    = GGA_FUNCTIONAL("Combine",    1, gga_key_read);
Functional DBLYPFunctional      = GGA_FUNCTIONAL("DBLYP",      1, dblyp_read);
Functional DBP86Functional      = GGA_FUNCTIONAL("DBP86",      1, dbp86_read);
Functional DBPW91Functional     = GGA_FUNCTIONAL("DBPW91",     1, dbpw91_read);
Functional EDF1Functional       = GGA_FUNCTIONAL("EDF1",       1, edf1_read);
Functional EDF2Functional       = GGA_FUNCTIONAL("EDF2",       1, edf2_read);
Functional GGAKeyFunctional     = GGA_FUNCTIONAL("GGAKey",     1, gga_key_read);
Functional G96VWNFunctional     = GGA_FUNCTIONAL("G96VWN",     1, g96vwn_read);
Functional G96LYPFunctional     = GGA_FUNCTIONAL("G96LYP",     1, g96lyp_read);
Functional G96P86Functional     = GGA_FUNCTIONAL("G96P86",     1, g96p86_read);
Functional G96PW91Functional    = GGA_FUNCTIONAL("G96PW91",    1, g96pw91_read);
Functional G961LYPFunctional    = GGA_FUNCTIONAL("G961LYP",    1, g961lyp_read);
Functional HCTHFunctional       = GGA_FUNCTIONAL("HCTH",       1, hcth_read);
Functional KMLYPFunctional      = GGA_FUNCTIONAL("KMLYP",      3, kmlyp_read);
Functional KT1Functional        = GGA_FUNCTIONAL("KT1",        1, kt1_read);
Functional KT2Functional        = GGA_FUNCTIONAL("KT2",        1, kt2_read);
Functional KT3Functional        = GGA_FUNCTIONAL("KT3",        1, kt3_read);
Functional LG1LYPFunctional     = GGA_FUNCTIONAL("LG1LYP",     1, lg1lyp_read);
Functional mPWVWNFunctional     = GGA_FUNCTIONAL("mPWVWN",     1, mpwvwn_read);
Functional mPWLYPFunctional     = GGA_FUNCTIONAL("mPWLYP",     1, mpwlyp_read);
Functional mPWP86Functional     = GGA_FUNCTIONAL("mPWP86",     1, mpwp86_read);
Functional mPWPW91Functional    = GGA_FUNCTIONAL("mPWPW91",    1, mpwpw91_read);
Functional mPW91Functional      = GGA_FUNCTIONAL("mPW91",      1, mpwpw91_read);   // mPW91 aliases mPWPW91
Functional mPW3PW91Functional   = GGA_FUNCTIONAL("mPW3PW91",   1, mpw3pw91_read);
Functional mPW1PW91Functional   = GGA_FUNCTIONAL("mPW1PW91",   1, mpw1pw91_read);
Functional mPW1KFunctional      = GGA_FUNCTIONAL("mPW1K",      1, mpw1k_read);
Functional mPW1NFunctional      = GGA_FUNCTIONAL("mPW1N",      1, mpw1n_read);
Functional mPW1SFunctional      = GGA_FUNCTIONAL("mPW1S",      1, mpw1s_read);
Functional OVWNFunctional       = GGA_FUNCTIONAL("OVWN",       1, ovwn_read);
Functional OLYPFunctional       = GGA_FUNCTIONAL("OLYP",       1, olyp_read);
Functional OP86Functional       = GGA_FUNCTIONAL("OP86",       1, op86_read);
Functional OPW91Functional      = GGA_FUNCTIONAL("OPW91",      1, opw91_read);
Functional PBE0Functional       = GGA_FUNCTIONAL("PBE0",       1, pbe0_read);
Functional PBE1PBEFunctional    = GGA_FUNCTIONAL("PBE1PBE",    1, pbe0_read);
Functional PBE0PBEFunctional    = GGA_FUNCTIONAL("PBE0PBE",    1, pbe0_read);
Functional PBEFunctional        = GGA_FUNCTIONAL("PBE",        1, pbe_read);
Functional PBEPBEFunctional     = GGA_FUNCTIONAL("PBEPBE",     1, pbe_read);
Functional revPBEFunctional     = GGA_FUNCTIONAL("revPBE",     1, revpbe_read);
Functional RPBEFunctional       = GGA_FUNCTIONAL("RPBE",       1, rpbe_read);
Functional mPBEFunctional       = GGA_FUNCTIONAL("mPBE",       1, mpbe_read);
Functional PW91Functional       = GGA_FUNCTIONAL("PW91",       1, pw91_read);
Functional PW91VWNFunctional    = GGA_FUNCTIONAL("PW91VWN",    1, pw91vwn_read);
Functional PW91LYPFunctional    = GGA_FUNCTIONAL("PW91LYP",    1, pw91lyp_read);
Functional PW91P86Functional    = GGA_FUNCTIONAL("PW91P86",    1, pw91p86_read);
Functional PW91PW91Functional   = GGA_FUNCTIONAL("PW91PW91",   1, pw91_read);
Functional XLYPFunctional       = GGA_FUNCTIONAL("XLYP",       1, xlyp_read);
Functional X3LYPFunctional      = GGA_FUNCTIONAL("X3LYP",      1, x3lyp_read);

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

int
clear_funclist()
{
    gga_fun_list = NULL;
    return 1;
}

static integer
xalpha_read(const char* conf_line)
{
    real weight;
    integer res = (sscanf(conf_line, "%lg", &weight)==1);
    if(res) 
        gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 
                                      1.5*weight);
    fun_set_hf_weight(0);
    return res;
}

static integer
lda_read(const char* conf_line)
{
    fun_set_hf_weight(0);
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
ldagauss_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &VWN3Functional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
b86vwn_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &B86xFunctional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,     1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
b86lyp_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &B86xFunctional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
b86p86_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &B86xFunctional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &PZ81Functional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
b86pw91_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &B86xFunctional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional,   1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
bvwn_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
blyp_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
b1lyp_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  0.75);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   0.75);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     1.0);
    fun_set_hf_weight(0.25);
    return 1;
}


static integer
b2plyp_read(const char* conf_line)
{
    static const real lypw = 0.73, dirw = 0.47;
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   0.47);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     lypw);
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    fun_set_hf_weight(1-dirw);
    fun_set_mp2_weight(0.27);
    return 1;
}


static integer
b2tplyp_read(const char* conf_line)
{
    static const real lypw = 0.69, dirw = 0.40;
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   0.40);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     lypw);
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    fun_set_hf_weight(1-dirw);
    fun_set_mp2_weight(0.31);
    return 1;
}

static integer
mpw2plyp_read(const char* conf_line)
{
    static const real lypw = 0.75, mpwxw = 0.45;
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     lypw);
    gga_fun_list = add_functional(gga_fun_list, &mPWxFunctional,  mpwxw);
    fun_set_hf_weight(1-mpwxw);
    fun_set_mp2_weight(0.25);
    return 1;
}

static integer
mpw2kplyp_read(const char* conf_line)
{
    static const real lypw = 0.58, mpwxw = 0.28;
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     lypw);
    gga_fun_list = add_functional(gga_fun_list, &mPWxFunctional,  mpwxw);
    fun_set_hf_weight(1-mpwxw);
    fun_set_mp2_weight(0.42);
    return 1;
}

static integer
b2gpplyp_read(const char* conf_line)
{
    static const real lypw = 0.64, dirw = 0.35;
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   0.35);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     lypw);
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    fun_set_hf_weight(1-dirw);
    fun_set_mp2_weight(0.36);
    return 1;
}

static integer
b2piplyp_read(const char* conf_line)
{
    static const real lypw = 0.727, dirw = 0.398;
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   0.398);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     lypw);
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    fun_set_hf_weight(1-dirw);
    fun_set_mp2_weight(0.273);
    return 1;
}

static integer
pbe0dh_read(const char* conf_line)
{
    static const real pbecw = 0.875, pbexw = 0.50;
    gga_fun_list = add_functional(gga_fun_list, &PBEcFunctional,    pbecw);
    gga_fun_list = add_functional(gga_fun_list, &PBExFunctional,    pbexw);
    fun_set_hf_weight(1-pbexw);
    fun_set_mp2_weight(0.125);
    return 1;
}

static integer
b3lyp_read(const char* conf_line)
{
    static const real lypw = 0.81, dirw = 0.8;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   0.72);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     lypw);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,     1-lypw);
    fun_set_hf_weight(1-dirw);
    return 1;
}

static integer
b3lypg_read(const char* conf_line)
{
    static const real lypw = 0.81, dirw = 0.8;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,   dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,    0.72);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,      lypw);
    gga_fun_list = add_functional(gga_fun_list, &VWN3IFunctional,   1-lypw);
    fun_set_hf_weight(1-dirw);
    return 1;
}

static integer
b3pw91_read(const char* conf_line)
{
    static const real corw = 0.81, dirw = 0.8;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,   dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,    0.72);
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional,    corw);
    gga_fun_list = add_functional(gga_fun_list, &PW92cFunctional,    1-corw);
    fun_set_hf_weight(1-dirw);
    return 1;
}

static integer
b1pw91_read(const char* conf_line)
{
    static const real hfw=0.25;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,   1-hfw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,    1-hfw);
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional,    1.0);
    fun_set_hf_weight(hfw);
    return 1;
}

static integer
bhandh_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 0.5);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,    1.0);
    fun_set_hf_weight(0.5);
    return 1;
}

static integer
bhandhlyp_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 0.5);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,  0.5);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,    1.0);
    fun_set_hf_weight(0.5);
    return 1;
}

static integer
bp86_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &PZ81Functional,   1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
b3p86_read(const char* conf_line)
{
    static const real p86w = 0.81, dirw = 0.8;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   0.72);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,    p86w);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,     1.0);
    fun_set_hf_weight(1-dirw);
    return 1;
}

static integer
b3p86g_read(const char* conf_line)
{
    static const real p86w = 0.81, dirw = 0.8;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   0.72);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,    p86w);
    gga_fun_list = add_functional(gga_fun_list, &VWN3Functional,    1.0);
    fun_set_hf_weight(1-dirw);
    return 1;
}

static integer
bpw91_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional,   1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
bw_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &WignerFunctional,  1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
bfw_read(const char* conf_line)
{
    static const real xw = 0.736, cw = 1.178, hfw = 0.286;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  xw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   xw);
    gga_fun_list = add_functional(gga_fun_list, &WignerFunctional,  cw);
    fun_set_hf_weight(hfw);
    return 1;
}

static integer
dblyp_read(const char* conf_line)
{
    static const real dirw = 1.030952, mbw = 10.4017, bw = -8.44793;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   bw);
    gga_fun_list = add_functional(gga_fun_list, &mBeckeFunctional,  mbw);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
dbpw91_read(const char* conf_line)
{
    static const real dirw = 1.030952, mbw = 10.4017, bw = -8.44793;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   bw);
    gga_fun_list = add_functional(gga_fun_list, &mBeckeFunctional,  mbw);
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional,   1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
dbp86_read(const char* conf_line)
{
    static const real dirw = 1.030952, mbw = 10.4017, bw = -8.44793;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   bw);
    gga_fun_list = add_functional(gga_fun_list, &mBeckeFunctional,  mbw);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &PZ81Functional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
edf1_read(const char* conf_line)
{
    static const real dirw = 1.030952, mbw = 10.4017, bw = -8.44793;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   bw);
    gga_fun_list = add_functional(gga_fun_list, &mBeckeFunctional,  mbw);
    gga_fun_list = add_functional(gga_fun_list, &LYPrFunctional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
edf2_read(const char* conf_line)
{
    static const real dirw = 0.2811, bw = 0.6227, mbw = -0.0551;
    static const real vwnw = 0.3029, lypw = 0.5998, lyprw = -0.0053;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   bw);
    gga_fun_list = add_functional(gga_fun_list, &mBeckeFunctional,  mbw);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,     vwnw);
    gga_fun_list = add_functional(gga_fun_list, &LYPrFunctional,    lyprw);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     lypw);
    fun_set_hf_weight(0.1695);
    return 1;
}

static integer
g96vwn_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &G96xFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
g96lyp_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &G96xFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
g96p86_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &G96xFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &PZ81Functional,   1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
g96pw91_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &G96xFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional,  1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
g961lyp_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  0.75);
    gga_fun_list = add_functional(gga_fun_list, &G96xFunctional,    0.75);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     1.0);
    fun_set_hf_weight(0.25);
    return 1;
}

static integer
hcth_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &HCTH407Functional, 1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
kmlyp_read(const char* conf_line)
{
    static const real hfw = 0.557, lypw = 0.448;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1-hfw);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,     1-lypw);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     lypw);
    fun_set_hf_weight(hfw);
    return 1;
}

static integer
kt1_read(const char* conf_line)
{
    static const real ktgam = -0.006;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &KTxFunctional,     ktgam);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,     1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
kt2_read(const char* conf_line)
{   
    static const real dirw = 1.07173, vwnw = 0.576727; 
    static const real ktgam = -0.006;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, dirw);
    gga_fun_list = add_functional(gga_fun_list, &KTxFunctional,    ktgam);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,    vwnw);
    fun_set_hf_weight(0);
    return 1;
}

static integer
kt3_read(const char* conf_line)
{
    static const real dirw = 1.092, lypw = 0.864409, optw = -0.925452;
    static const real ktgam = -0.004;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, dirw);
    gga_fun_list = add_functional(gga_fun_list, &KTxFunctional,    ktgam);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,    lypw);
    gga_fun_list = add_functional(gga_fun_list, &OPTXFunctional,   optw);
    fun_set_hf_weight(0);
    return 1;
}

static integer
lg1lyp_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  0.75);
    gga_fun_list = add_functional(gga_fun_list, &LG93xFunctional,   0.75);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     1.0);
    fun_set_hf_weight(0.25);
    return 1;
}

static integer
mpwvwn_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &mPWxFunctional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,     1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
mpwlyp_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &mPWxFunctional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
mpwp86_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &mPWxFunctional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &PZ81Functional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
mpwpw91_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &mPWxFunctional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional,   1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
mpw3pw91_read(const char* conf_line)
{
    static const real hfw=0.2, wmp=0.72, lypw=0.81;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1-hfw);
    gga_fun_list = add_functional(gga_fun_list, &mPWxFunctional,    wmp);
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional,   lypw);
    gga_fun_list = add_functional(gga_fun_list, &PW92cFunctional,   1-lypw);
    fun_set_hf_weight(hfw);
    return 1;
}

static integer
mpw1pw91_read(const char* conf_line)
{
    static const real hfw=0.25;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1-hfw);
    gga_fun_list = add_functional(gga_fun_list, &mPWxFunctional,    1-hfw);
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional,   1.0);
    fun_set_hf_weight(hfw);
    return 1;
}

static integer
mpw1k_read(const char* conf_line)
{
    static const real hfw=0.428;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1-hfw);
    gga_fun_list = add_functional(gga_fun_list, &mPWxFunctional,    1-hfw);
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional,   1.0);
    fun_set_hf_weight(hfw);
    return 1;
}

static integer
mpw1n_read(const char* conf_line)
{
    static const real hfw=0.406;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1-hfw);
    gga_fun_list = add_functional(gga_fun_list, &mPWxFunctional,    1-hfw);
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional,   1.0);
    fun_set_hf_weight(hfw);
    return 1;
}

static integer
mpw1s_read(const char* conf_line)
{
    static const real hfw=0.060;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1-hfw);
    gga_fun_list = add_functional(gga_fun_list, &mPWxFunctional,    1-hfw);
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional,   1.0);
    fun_set_hf_weight(hfw);
    return 1;
}

static integer
ovwn_read(const char* conf_line)
{
    static const real optw = -1.43169, dirw = 1.05151;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &OPTXFunctional,    optw);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,     1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
olyp_read(const char* conf_line)
{
    static const real optw = -1.43169, dirw = 1.05151;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &OPTXFunctional,    optw);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
op86_read(const char* conf_line)
{
    static const real optkw = -1.43169, dirw = 1.05151;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &OPTXFunctional,    optkw);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &PZ81Functional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
opw91_read(const char* conf_line)
{
    static const real optkw = -1.43169, dirw = 1.05151;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &OPTXFunctional,    optkw);
    gga_fun_list = add_functional(gga_fun_list, &PW91Functional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
pbe_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &PBEcFunctional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &PBExFunctional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}
 
static integer
pbe0_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &PBEcFunctional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &PBExFunctional,    0.75);
    fun_set_hf_weight(0.25);
    return 1;
}

static integer
revpbe_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &PBEcFunctional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &revPBExFunctional, 1.0);
    fun_set_hf_weight(0);
    return 1;
}
 
static integer
rpbe_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &PBEcFunctional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &RPBExFunctional,   1.0);
    fun_set_hf_weight(0);
    return 1;
}
 
static integer
mpbe_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &PBEcFunctional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &mPBExFunctional,   1.0);
    fun_set_hf_weight(0);
    return 1;
}
 
static integer
pw1pw_read(const char* conf_line)
{
    static const real hfw=0.25;
    gga_fun_list = add_functional(gga_fun_list, &PW91xFunctional,   1-hfw);
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional,   1.0);
    fun_set_hf_weight(hfw);
    return 1;
}

static integer
pw91_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &PW91xFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional,   1.0);
    fun_set_hf_weight(0);
    return 1;
}

static integer
pw91vwn_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &PW91xFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,     1.0);
    fun_set_hf_weight(0);
    return 1;
} 

static integer
pw91lyp_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     1.0);
    gga_fun_list = add_functional(gga_fun_list, &PW91xFunctional,   1.0);
    fun_set_hf_weight(0);
    return 1;
} 

static integer
pw91p86_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &PZ81Functional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &PW91xFunctional,   1.0);
    fun_set_hf_weight(0);
    return 1;
} 

static integer
xlyp_read(const char* conf_line)
{
    static const real b88w = 0.722, pw91w = 0.347;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1-pw91w);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   b88w);
    gga_fun_list = add_functional(gga_fun_list, &PW91xFunctional,   pw91w);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     1.0);
    fun_set_hf_weight(0);
    return 1;
} 

static integer
x3lyp_read(const char* conf_line)
/* PW91x contains Slater exchange, so subtract from dirw */
{
    static const real dirw = 0.782, b88w = 0.542, pw91w = 0.167;
    static const real vwnw = 0.129, lypw = 0.871;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw-pw91w);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   b88w);
    gga_fun_list = add_functional(gga_fun_list, &PW91xFunctional,   pw91w);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,     vwnw);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     lypw);
    fun_set_hf_weight(0.218);
    return 1;
} 


/* gga_key_read:
   read general GGA input: 
*/
static integer
gga_key_read(const char* conf_line)
{
    integer res = 1, i;
    real f;
    const char* str = conf_line;

    fun_set_hf_weight(0);

    while(*str) {
        while(*str && isspace((integer)*str)) str++; /* skip whitespace */
        if(*str =='\0') break; /* line ended by whitespace */
        if(strncasecmp("HF=", str, 3)==0) {
            if(sscanf(str+3,"%lg", &f) != 1) {
                fun_printf("GGAKey: HF not followed by the weight: ",
                           conf_line);
                res = 0;
            } else fun_set_hf_weight(f);
        } else if(strncasecmp("MP2=", str, 3)==0) {
            if(sscanf(str+3,"%lg", &f) != 1) {
                fun_printf("GGAKey: MP2 not followed by the weight: ",
                           conf_line);
                res = 0;
            } else fun_set_mp2_weight(f);
        } else {
            for(i=0; available_functionals[i]; i++) {
                integer len = strlen(available_functionals[i]->name);
                if(strncasecmp(available_functionals[i]->name, str, len)==0 &&
                   str[len] == '=') {
                    if(sscanf(str+len+1,"%lg", &f) != 1) {
                        fun_printf("GGAKey: keyword '%s' not followed by "
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
                fun_printf("GGAKey: functional '%s' not recognised: ", str);
                res = 0;
            }
        }
        while(*str && !isspace((integer)*str)) str++; /* skip nonws */
    } return res;
}

static void
gga_report(void)
{
    FuncList* lst;
    fun_printf("Weighted mixed functional:");
    if(fun_get_hf_weight()>0)
        fun_printf("%25s: %10.5f", "HF exchange", fun_get_hf_weight());
    for(lst=gga_fun_list; lst; lst=lst->next) 
        fun_printf("%25s: %10.5f", lst->func->name, lst->weight);
}

static real
gga_energy(const FunDensProp* dp)
{
    real res = 0;
    FuncList* lst;
    for(lst=gga_fun_list; lst; lst=lst->next) {
        real contr = lst->weight*lst->func->func(dp);
/*        fun_printf("[%g,%g,w=%g] %s contributes with %g", dp->rhoa,
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
/*      fun_printf("[%g,%g,w=%g] %s f: %g deriv (%g,%g)", dp->rhoa,
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


