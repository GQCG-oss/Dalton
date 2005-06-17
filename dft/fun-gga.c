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
#include <stdlib.h>

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int  lda_read(const char* conf_line);
static real lda_energy(const FunDensProp* dp);
static void lda_first(FunFirstFuncDrv *ds,   real fac, const FunDensProp* dp);
static void lda_second(FunSecondFuncDrv *ds, real fac, const FunDensProp* dp);
static void lda_third(FunThirdFuncDrv *ds,   real fac, const FunDensProp* dp);
static void lda_fourth(FunFourthFuncDrv *ds, real fac, const FunDensProp* dp);
static int  ldagauss_read(const char* conf_line);
static int  b86vwn_read(const char* conf_line);
static int  b86lyp_read(const char* conf_line);
static int  b86p86_read(const char* conf_line);
static int  b86pw91_read(const char* conf_line);
static int  bvwn_read(const char* conf_line);
static int  blyp_read(const char* conf_line);
static int  b3lyp_read(const char* conf_line);
static int  b3lypg_read(const char* conf_line);
static int  bp86_read(const char* conf_line);
static int  b3p86_read(const char* conf_line);
static int  b3p86g_read(const char* conf_line);
static int  bpw91_read(const char* conf_line);
static int  b3pw91_read(const char* conf_line);
static int  b3pw91g_read(const char* conf_line);
static int  bhandh_read(const char* conf_line);
static int  bhandhlyp_read(const char* conf_line);
static int  bw_read(const char* conf_line);
static int  dblyp_read(const char* conf_line);
static int  dbp86_read(const char* conf_line);
static int  dbpw91_read(const char* conf_line);
static int  edf1_read(const char* conf_line);
static int  edf2_read(const char* conf_line);
static int  g96vwn_read(const char* conf_line);
static int  g96lyp_read(const char* conf_line);
static int  g96p86_read(const char* conf_line);
static int  g96pw91_read(const char* conf_line);
static int  kt1_read(const char* conf_line);
static int  kt2_read(const char* conf_line);
static int  kt3_read(const char* conf_line);
static int  mpwvwn_read(const char* conf_line);
static int  mpwlyp_read(const char* conf_line);
static int  mpwp86_read(const char* conf_line);
static int  mpwpw91_read(const char* conf_line);
static int  ovwn_read(const char* conf_line);
static int  olyp_read(const char* conf_line);
static int  op86_read(const char* conf_line);
static int  opw91_read(const char* conf_line);
static int  pbe_read(const char* conf_line);
static int  pbe0_read(const char* conf_line);
static int  pw91_read(const char* conf_line);
static int  pw91vwn_read(const char* conf_line);
static int  pw91lyp_read(const char* conf_line);
static int  pw91p86_read(const char* conf_line);

static int  gga_isgga(void);
static int  xalpha_read(const char* conf_line);
static int  gga_key_read(const char* conf_line);
static void gga_report(void);
static real gga_energy(const FunDensProp* dp);
static void gga_first(FunFirstFuncDrv *ds,   real fac, const FunDensProp* dp);
static void gga_second(FunSecondFuncDrv *ds, real fac, const FunDensProp* dp);
static void gga_third(FunThirdFuncDrv *ds,   real fac, const FunDensProp* dp);
static void gga_fourth(FunFourthFuncDrv *ds, real fac, const FunDensProp* dp);

#define LDA_FUNCTIONAL(name,read) { (name), \
    fun_false, (read), NULL, lda_energy, lda_first, lda_second, \
    lda_third, lda_fourth }

#define GGA_FUNCTIONAL(name,read) { (name), \
    gga_isgga, (read), gga_report, gga_energy, gga_first, gga_second, \
    gga_third, gga_fourth }

Functional XAlphaFunctional = GGA_FUNCTIONAL("XAlpha", xalpha_read);
Functional LDAFunctional =    LDA_FUNCTIONAL("LDA",     lda_read);
/* SVWN5 aliases LDA */
Functional SVWN5Functional =  LDA_FUNCTIONAL("SVWN5",   lda_read);
Functional SVWN3Functional =  GGA_FUNCTIONAL("SVWN3",   ldagauss_read);
Functional B3LYPFunctional =  GGA_FUNCTIONAL("B3LYP",   b3lyp_read);
Functional B3LYPGaussFunctional = GGA_FUNCTIONAL("B3LYP-G", b3lypg_read);
Functional B3P86Functional =  GGA_FUNCTIONAL("B3P86",   b3p86_read);
Functional B3P86GFunctional = GGA_FUNCTIONAL("B3P86-G", b3p86g_read);
Functional B3PW91Functional =  GGA_FUNCTIONAL("B3PW91", b3pw91_read);
Functional B3PW91GFunctional = GGA_FUNCTIONAL("B3PW91G", b3pw91g_read);
Functional BHandHFunctional = GGA_FUNCTIONAL("BHandH",  bhandh_read);
Functional BHandHLYPFunctional = GGA_FUNCTIONAL("BHandHLYP", bhandhlyp_read);
Functional B86VWNFunctional =   GGA_FUNCTIONAL("B86VWN",  b86vwn_read);
Functional B86LYPFunctional =   GGA_FUNCTIONAL("B86LYP",  b86lyp_read);
Functional B86P86Functional =   GGA_FUNCTIONAL("B86P86",  b86p86_read);
Functional B86PW91Functional =  GGA_FUNCTIONAL("B86PW91", b86pw91_read);
Functional BVWNFunctional =   GGA_FUNCTIONAL("BVWN",    bvwn_read);
Functional BLYPFunctional =   GGA_FUNCTIONAL("BLYP",    blyp_read);
Functional BP86Functional =   GGA_FUNCTIONAL("BP86",    bp86_read);
Functional BPW91Functional =  GGA_FUNCTIONAL("BPW91",   bpw91_read);
Functional BWFunctional =     GGA_FUNCTIONAL("BW",      bw_read);
Functional CombineFunctional =GGA_FUNCTIONAL("Combine", gga_key_read);
Functional DBLYPFunctional =   GGA_FUNCTIONAL("DBLYP",  dblyp_read);
Functional DBP86Functional =   GGA_FUNCTIONAL("DBP86",  dbp86_read);
Functional DBPW91Functional =  GGA_FUNCTIONAL("DBPW91", dbpw91_read);
Functional EDF1Functional =   GGA_FUNCTIONAL("EDF1",    edf1_read);
Functional EDF2Functional =   GGA_FUNCTIONAL("EDF2",    edf2_read);
Functional GGAKeyFunctional = GGA_FUNCTIONAL("GGAKey",  gga_key_read);
Functional G96VWNFunctional = GGA_FUNCTIONAL("G96VWN",  g96vwn_read);
Functional G96LYPFunctional = GGA_FUNCTIONAL("G96LYP",  g96lyp_read);
Functional G96P86Functional = GGA_FUNCTIONAL("G96P86",  g96p86_read);
Functional G96PW91Functional = GGA_FUNCTIONAL("G96PW91", g96pw91_read);
Functional KT1Functional =    GGA_FUNCTIONAL("KT1",     kt1_read);
Functional KT2Functional =    GGA_FUNCTIONAL("KT2",     kt2_read);
Functional KT3Functional =    GGA_FUNCTIONAL("KT3",     kt3_read);
Functional mPWVWNFunctional = GGA_FUNCTIONAL("mPWVWN" , mpwvwn_read);
Functional mPWLYPFunctional = GGA_FUNCTIONAL("mPWLYP" , mpwlyp_read);
Functional mPWP86Functional = GGA_FUNCTIONAL("mPWP86" , mpwp86_read);
Functional mPWPW91Functional = GGA_FUNCTIONAL("mPWPW91" , mpwpw91_read);
Functional OVWNFunctional =   GGA_FUNCTIONAL("OVWN" ,   ovwn_read);
Functional OLYPFunctional =   GGA_FUNCTIONAL("OLYP" ,   olyp_read);
Functional OP86Functional =   GGA_FUNCTIONAL("OP86" ,   op86_read);
Functional OPW91Functional =   GGA_FUNCTIONAL("OPW91" ,  opw91_read);
Functional PBE0Functional =   GGA_FUNCTIONAL("PBE0",    pbe0_read);
Functional PBEFunctional =    GGA_FUNCTIONAL("PBE",     pbe_read); 
Functional PW91Functional =   GGA_FUNCTIONAL("PW91",    pw91_read); 
Functional PW91VWNFunctional =  GGA_FUNCTIONAL("PW91VWN", pw91vwn_read); 
Functional PW91LYPFunctional =  GGA_FUNCTIONAL("PW91LYP", pw91lyp_read); 
Functional PW91P86Functional =  GGA_FUNCTIONAL("PW91P86", pw91p86_read); 
Functional PW91PW91Functional =  GGA_FUNCTIONAL("PW91PW91", pw91_read); 

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
    fun_set_hf_weight(0);
    return res;
}

static int
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

static int
ldagauss_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &VWN3Functional,  1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
b86vwn_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &B86xFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
b86lyp_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &B86xFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
b86p86_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &B86xFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &PZ81Functional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
b86pw91_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &B86xFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
bvwn_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
blyp_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,   1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
b3lyp_read(const char* conf_line)
{
    static const real lypw = 0.81, dirw = 0.8;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,   dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,    0.72);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,      lypw);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,    1-lypw);
    fun_set_hf_weight(1-dirw);
    return 1;
}

static int
b3lypg_read(const char* conf_line)
{
    static const real lypw = 0.81, dirw = 0.8;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,   dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,    0.72);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,      lypw);
    gga_fun_list = add_functional(gga_fun_list, &VWN3Functional,   1-lypw);
    fun_set_hf_weight(1-dirw);
    return 1;
}

static int
b3pw91_read(const char* conf_line)
{
    static const real corw = 0.81, dirw = 0.8;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,   dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,    0.72);
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional,    corw);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,    1-corw);
    fun_set_hf_weight(1-dirw);
    return 1;
}

static int
b3pw91g_read(const char* conf_line)
{
    static const real corw = 0.81, dirw = 0.8;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,   dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,    0.72);
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional,    corw);
    gga_fun_list = add_functional(gga_fun_list, &VWN3Functional,   1-corw);
    fun_set_hf_weight(1-dirw);
    return 1;
}

static int
bhandh_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 0.5);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,    1.0);
    fun_set_hf_weight(0.5);
    return 1;
}

static int
bhandhlyp_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 0.5);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,  0.5);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,    1.0);
    fun_set_hf_weight(0.5);
    return 1;
}

static int
bp86_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &PZ81Functional,  1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
b3p86_read(const char* conf_line)
{
    static const real lypw = 0.81, dirw = 0.8;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   0.72);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,    lypw);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,        1);
    fun_set_hf_weight(1-dirw);
    return 1;
}

static int
b3p86g_read(const char* conf_line)
{
    static const real lypw = 0.81, dirw = 0.8;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,  0.72);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,   lypw);
    gga_fun_list = add_functional(gga_fun_list, &VWN3Functional,   1);
    fun_set_hf_weight(1-dirw);
    return 1;
}

static int
bpw91_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,  1);
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional,  1);
    fun_set_hf_weight(0);
    return 1;
}

static int
bw_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &WignerFunctional, 1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
dblyp_read(const char* conf_line)
{
    static const real dirw = 1.030952, b35w = 10.4017, b42w = -8.44793;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   b42w);
    gga_fun_list = add_functional(gga_fun_list, &Becke35Functional, b35w);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
dbpw91_read(const char* conf_line)
{
    static const real dirw = 1.030952, b35w = 10.4017, b42w = -8.44793;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   b42w);
    gga_fun_list = add_functional(gga_fun_list, &Becke35Functional, b35w);
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional,   1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
dbp86_read(const char* conf_line)
{
    static const real dirw = 1.030952, b35w = 10.4017, b42w = -8.44793;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   b42w);
    gga_fun_list = add_functional(gga_fun_list, &Becke35Functional, b35w);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &PZ81Functional,     1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
edf1_read(const char* conf_line)
{
    static const real dirw = 1.030952, b35w = 10.4017, b42w = -8.44793;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   b42w);
    gga_fun_list = add_functional(gga_fun_list, &Becke35Functional, b35w);
    gga_fun_list = add_functional(gga_fun_list, &LYPrFunctional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
edf2_read(const char* conf_line)
{
    static const real dirw = 0.2811, b42w = 0.6227, b35w = -0.0551;
    static const real vwnw = 0.3029, lypw = 0.5998, lyprw = -0.0053;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,   b42w);
    gga_fun_list = add_functional(gga_fun_list, &Becke35Functional, b35w);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,     vwnw);
    gga_fun_list = add_functional(gga_fun_list, &LYPrFunctional,    lyprw);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,     lypw);
    fun_set_hf_weight(0.1695);
    return 1;
}

static int
g96vwn_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &G96xFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
g96lyp_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &G96xFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
g96p86_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &G96xFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &PZ81Functional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
g96pw91_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &G96xFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional,  1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
kt1_read(const char* conf_line)
{
    static const real ktgam = -0.006;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &KTFunctional,    ktgam);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,     1.0);
    fun_set_hf_weight(0);
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
    fun_set_hf_weight(0);
    return 1;
}

static int
kt3_read(const char* conf_line)
{
    static const real dirw = 1.092, lypw = 0.864409, optw = -0.925452;
    static const real ktgam = -0.004;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, dirw);
    gga_fun_list = add_functional(gga_fun_list, &KTFunctional,    ktgam);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,    lypw);
    gga_fun_list = add_functional(gga_fun_list, &OPTXFunctional,   optw);
    fun_set_hf_weight(0);
    return 1;
}

static int
mpwvwn_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &mPWFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,   1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
mpwlyp_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &mPWFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,   1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
mpwp86_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &mPWFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &PZ81Functional,   1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
mpwpw91_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &mPWFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional,   1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
ovwn_read(const char* conf_line)
{
    static const real optw = -1.43169, dirw = 1.05151;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, dirw);
    gga_fun_list = add_functional(gga_fun_list, &OPTXFunctional,   optw);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
olyp_read(const char* conf_line)
{
    static const real optw = -1.43169, dirw = 1.05151;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, dirw);
    gga_fun_list = add_functional(gga_fun_list, &OPTXFunctional,   optw);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
op86_read(const char* conf_line)
{
    static const real optkw = -1.43169, dirw = 1.05151;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, dirw);
    gga_fun_list = add_functional(gga_fun_list, &OPTXFunctional,   optkw);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,    1.0);
    gga_fun_list = add_functional(gga_fun_list, &PZ81Functional,    1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
opw91_read(const char* conf_line)
{
    static const real optkw = -1.43169, dirw = 1.05151;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, dirw);
    gga_fun_list = add_functional(gga_fun_list, &OPTXFunctional,   optkw);
    gga_fun_list = add_functional(gga_fun_list, &PW91Functional,   1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
pbe_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &PBEcFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &PBExFunctional, 1.0);
    fun_set_hf_weight(0);
    return 1;
}
 
static int
pbe0_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &PBEcFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &PBExFunctional, 0.75);
    fun_set_hf_weight(0.25);
    return 1;
}

static int
pw91_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &PW91cFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &PW91xFunctional, 1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
pw91vwn_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &PW91xFunctional, 1.0);
    fun_set_hf_weight(0);
    return 1;
} 

static int
pw91lyp_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &PW91xFunctional, 1.0);
    fun_set_hf_weight(0);
    return 1;
} 

static int
pw91p86_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &PZ81Functional,   1.0);
    gga_fun_list = add_functional(gga_fun_list, &PW91xFunctional, 1.0);
    fun_set_hf_weight(0);
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

    fun_set_hf_weight(0);

    while(*str) {
        while(*str && isspace((int)*str)) str++; /* skip whitespace */
        if(*str =='\0') break; /* line ended by whitespace */
        if(strncasecmp("HF=", str, 3)==0) {
            if(sscanf(str+3,"%g", &f) != 1) {
                fun_printf("GGAKey: HF not followed by the weight: ",
                           conf_line);
                res = 0;
            } else fun_set_hf_weight(f);
        } else {
            for(i=0; available_functionals[i]; i++) {
                int len = strlen(available_functionals[i]->name);
                if(strncasecmp(available_functionals[i]->name, str, len)==0 &&
                   str[len] == '=') {
                    if(sscanf(str+len+1,"%g", &f) != 1) {
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
        while(*str && !isspace((int)*str)) str++; /* skip nonws */
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


