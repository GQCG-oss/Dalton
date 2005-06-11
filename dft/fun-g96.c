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
/* fun-g96x.c:
   implementation of G96 (Gill 1996) exchange functional and its derivatives.
   #### this is just the gradient corrected term  for G96 functional#### 
   The full G96 exchange functional is given by G96x + Slater
   Reference: P.M.W. Gill, Mol. Phys. 
   implemented by Dave Wilson (davidwi@kjemi.uio.no)
   NOTE:
   this file may seem unnecessarily complex but the structure does pay off
   when implementing multiple functionals depending on different parameters.
*/

/* strictly conform to XOPEN ANSI C standard */
#if !defined(SYS_DEC)
/* XOPEN compliance is missing on old Tru64 4.0E Alphas */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1
#endif
#include <math.h>
#include <stddef.h>
#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int g96x_isgga(void) { return 1; }
static int g96x_read(const char* conf_line);
static real g96x_energy(const FunDensProp* dp);
static void g96x_first(FunFirstFuncDrv *ds,   real factor, 
                        const FunDensProp* dp);
static void g96x_second(FunSecondFuncDrv *ds, real factor, 
                        const FunDensProp* dp);
static void g96x_third(FunThirdFuncDrv *ds,   real factor, 
                        const FunDensProp* dp);
static void g96x_fourth(FunFourthFuncDrv *ds, real factor,
                       const FunDensProp* dp);

Functional G96xFunctional = {
    "G96x",         /* name */
    g96x_isgga,    /* gga-corrected */
    g96x_read,     /* set bloody common blocks */
    NULL,         /* reporter */
    g96x_energy, 
    g96x_first,
    g96x_second,
    g96x_third,
    g96x_fourth
};

/* IMPLEMENTATION PART */
static int
g96x_read(const char *conf_line)
{
    fun_set_hf_weight(0);
    return 1;
}

/* g96x_energy:
   note that in reality E_G96 = E_G96,alpha + E_G96,beta
   i.e the energy is linear in alpha and beta densities.

   G96 threshold is needed to avoid numerical problems on 0/0
   divisions.  The problems are small but it is better to be on the
   safe side.
 
   GAMMA = 1/137
*/
static const real G96_THRESHOLD = 1e-14;
static const real GAMMA = 0.007299270072993;
static real
g96x_energy(const FunDensProp* dp)
{
   real ea,eb;
   if (dp->rhoa<G96_THRESHOLD)
     ea=0.0;
   else {
     real ra43 = pow(dp->rhoa,4.0/3.0);
     real gra = dp->grada;
     real xa = gra/ra43;
     real xa12 = pow(dp->grada,1.0/2.0)*pow(dp->rhoa,-2.0/3.0);
     ea = ra43*GAMMA*xa12*xa;
   }
   if (dp->rhob<G96_THRESHOLD)
     eb = 0.0;
   else {
     real rb43 = pow(dp->rhob,4.0/3.0);
     real grb = dp->gradb;
     real xb = grb/rb43;
     real xb12 = pow(dp->gradb,1.0/2.0)*pow(dp->rhob,-2.0/3.0);
     eb = rb43*GAMMA*xb12*xb; 
   } 
   return -(ea+eb);
}


static void
g96x_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real ra13, ra43, gra, xa, xa12;
    real rb13, rb43, grb, xb, xb12;
    real faR, faZ;
    real fbR, fbZ;

    if (dp->rhoa >G96_THRESHOLD) {
        ra13 = pow(dp->rhoa,1.0/3.0);
        ra43 = pow(dp->rhoa,4.0/3.0);
        gra = dp->grada;
        xa = gra/ra43;
        xa12 = pow(dp->grada,1.0/2.0)*pow(dp->rhoa,-2.0/3.0);
        faR = 2.0/3.0*GAMMA*ra13*xa12*xa;
        faZ = -3.0/2.0*GAMMA*xa12;
        ds->df1000 += factor*faR;
        ds->df0010 += factor*faZ;
    }
    if (dp->rhob >G96_THRESHOLD) {
        rb13 = pow(dp->rhob,1.0/3.0);
        rb43 = pow(dp->rhob,4.0/3.0);
        grb = dp->gradb;
        xb = grb/rb43;
        xb12 = pow(dp->gradb,1.0/2.0)*pow(dp->rhob,-2.0/3.0);
        fbR = 2.0/3.0*GAMMA*rb13*xb12*xb;
        fbZ = -3.0/2.0*GAMMA*xb12;
        ds->df0100 += factor*fbR;
        ds->df0001 += factor*fbZ;
    } 
}


static void
g96x_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real ra13, ram23, ra43, ra, gra, xa, xa12;
    real rb13, rbm23, rb43, rb, grb, xb, xb12;
    real faR, faZ, faRR, faRZ, faZZ;
    real fbR, fbZ, fbRR, fbRZ, fbZZ;
    
    if (dp->rhoa >G96_THRESHOLD) {
        ra = dp->rhoa;
        gra = dp->grada;
        ra13 = pow(dp->rhoa,1.0/3.0);
        ra43 = ra13*ra;
        ram23 = ra13/ra;
        xa = gra/ra43;
        xa12 = pow(dp->grada,1.0/2.0)*pow(dp->rhoa,-2.0/3.0);
        faR = 2.0/3.0*GAMMA*ra13*xa12*xa;
        faZ = -3.0/2.0*GAMMA*xa12;
        faRR = -10.0/9.0*GAMMA*xa12*xa*ram23;
        faRZ = GAMMA*xa12/ra;
        faZZ = -3.0/4.0*GAMMA*xa12/(xa*ra43);
        ds->df1000 += factor*faR;
        ds->df0010 += factor*faZ;
        ds->df1010 += factor*faRZ;
        ds->df2000 += factor*faRR;
        ds->df0020 += factor*faZZ;
    }   
    if (dp->rhob >G96_THRESHOLD) {
        rb = dp->rhob;
        grb = dp->gradb;
        rb13 = pow(dp->rhob,1.0/3.0);
        rb43 = rb13*rb;
        rbm23 = rb13/rb;
        xb = grb/rb43;
        xb12 = pow(dp->gradb,1.0/2.0)*pow(dp->rhob,-2.0/3.0);
        fbR = 2.0/3.0*GAMMA*rb13*xb12*xb;
        fbZ = -3.0/2.0*GAMMA*xb12;
        fbRR = -10.0/9.0*GAMMA*xb12*xb*rbm23;
        fbRZ = GAMMA*xb12/rb;
        fbZZ = -3.0/4.0*GAMMA*xb12/(xb*rb43);
        ds->df0100 += factor*fbR;
        ds->df0001 += factor*fbZ;
        ds->df0101 += factor*fbRZ;
        ds->df0200 += factor*fbRR;
        ds->df0002 += factor*fbZZ;
    } 
}

static void
g96x_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real ra13, ram23, ra43, ram43, ra, gra, xa, xa12;
    real rb13, rbm23, rb43, rbm43, rb, grb, xb, xb12;
    real faR, faZ, faRR, faRZ, faZZ, faRRR, faRRZ, faRZZ, faZZZ;
    real fbR, fbZ, fbRR, fbRZ, fbZZ, fbRRR, fbRRZ, fbRZZ, fbZZZ;

    if (dp->rhoa >G96_THRESHOLD) {
        ra = dp->rhoa;
        gra = dp->grada;
        ra13 = pow(dp->rhoa,1.0/3.0);
        ra43 = ra13*ra;
        ram43 = 1.0/ra43;
        ram23 = ra13/ra;
        xa = gra/ra43;
        xa12 = pow(dp->grada,1.0/2.0)*pow(dp->rhoa,-2.0/3.0);
        faR = 2.0/3.0*GAMMA*ra13*xa12*xa;
        faZ = -3.0/2.0*GAMMA*xa12;
        faRR = -10.0/9.0*GAMMA*xa12*xa*ram23;
        faRZ = GAMMA*xa12/ra;
        faZZ = -3.0/4.0*GAMMA*xa12/(xa*ra43);
        faRRR = 80.0/27.0*GAMMA*xa12*xa*ram23/ra;
        faRRZ = -5.0/3.0*GAMMA*xa12/(ra*ra);
        faRZZ = 1.0/2.0*GAMMA*xa12/xa*ram43/ra;
        faZZZ = 3.0/8.0*GAMMA*xa12/xa/xa*ram43*ram43;
        ds->df1000 += factor*faR;
        ds->df0010 += factor*faZ;
        ds->df1010 += factor*faRZ;
        ds->df2000 += factor*faRR;
        ds->df0020 += factor*faZZ;
        ds->df3000 += factor*faRRR;
        ds->df2010 += factor*faRRZ;
        ds->df1020 += factor*faRZZ;
        ds->df0030 += factor*faZZZ;
    }
    if (dp->rhob >G96_THRESHOLD) {
        rb = dp->rhob;
        grb = dp->gradb;
        rb13 = pow(dp->rhob,1.0/3.0);
        rb43 = rb13*rb;
        rbm23 = rb13/rb;
        rbm43 = 1.0/rb43;
        xb = grb/rb43;
        xb12 = pow(dp->gradb,1.0/2.0)*pow(dp->rhob,-2.0/3.0);
        fbR = 2.0/3.0*rb13*GAMMA*xb12*xb;
        fbZ = -3.0/2.0*GAMMA*xb12;
        fbRR = -10.0/9.0*GAMMA*xb12*xb*rbm23;
        fbRZ = GAMMA*xb12/rb;
        fbZZ = -3.0/4.0*GAMMA*xb12/(xb*rb43);
        fbRRR = 80.0/27.0*GAMMA*xb12*xb*rbm23/rb;
        fbRRZ = -5.0/3.0*GAMMA*xb12/(rb*rb);
        fbRZZ = 1.0/2.0*GAMMA*xb12/xb*rbm43/rb;
        fbZZZ = 3.0/8.0*GAMMA*xb12/xb/xb*rbm43*rbm43;
        ds->df0100 += factor*fbR;
        ds->df0001 += factor*fbZ;
        ds->df0101 += factor*fbRZ;
        ds->df0200 += factor*fbRR;
        ds->df0002 += factor*fbZZ;
        ds->df0300 += factor*fbRRR;
        ds->df0201 += factor*fbRRZ;
        ds->df0102 += factor*fbRZZ;
        ds->df0003 += factor*fbZZZ;
    } 
}

static void
g96x_fourth(FunFourthFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real ra, ra2, ra3, ra13, ram23, ra43, ram43;
    real rb, rb2, rb3, rb13, rbm23, rb43, rbm43;
    real gra, xa, xa12, xa2, xa3;
    real grb, xb, xb12, xb2, xb3;
    real faR, faZ, faRR, faRZ, faZZ, faRRR, faRRZ, faRZZ, faZZZ;
    real fbR, fbZ, fbRR, fbRZ, fbZZ, fbRRR, fbRRZ, fbRZZ, fbZZZ;
    real faRRRR, faRRRZ, faRRZZ, faRZZZ, faZZZZ; 
    real fbRRRR, fbRRRZ, fbRRZZ, fbRZZZ, fbZZZZ; 

    if (dp->rhoa >G96_THRESHOLD) {
        ra = dp->rhoa;
        gra = dp->grada;
        ra13 = pow(dp->rhoa,1.0/3.0);
        ra43 = ra13*ra;
        ra2 = ra*ra;
        ra3 = ra2*ra;
        ram43 = pow(dp->rhoa,-4.0/3.0);
        ram23 = pow(dp->rhoa,-2.0/3.0);
        xa = gra*ram43;
        xa2 = xa*xa;
        xa3 = xa2*xa;
        xa12 = ram23*pow(dp->grada,1.0/2.0);
        faR = 2.0/3.0*GAMMA*ra13*xa12*xa;
        faZ = -3.0/2.0*GAMMA*xa12;
        faRR = -10.0/9.0*GAMMA*xa12*xa*ram23;
        faRZ = GAMMA*xa12/ra;
        faZZ = -3.0/4.0*GAMMA*xa12*ram43/xa;
        faRRR = 80.0/27.0*GAMMA*xa12*xa*ram23/ra;
        faRRZ = -5.0/3.0*GAMMA*xa12/(ra*ra);
        faRZZ = 1.0/2.0*GAMMA*xa12/xa*ram43/ra;
        faZZZ = 3.0/8.0*GAMMA*ram43*ram43*xa12/(xa*xa);
        faRRRR = -880.0/81.0*GAMMA*ram43*ram43*xa12*xa;
        faRRRZ = 40.0/9.0*GAMMA/ra3*xa12;
        faRRZZ = -5.0/6.0*GAMMA*xa12/xa*ram43/ra2;
        faRZZZ = -1.0/4.0*GAMMA*ram23/ra3*xa12/xa2;
        faZZZZ = -9.0/16.0*GAMMA/(ra3*ra)*xa12/xa3;
        ds->df1000 += factor*faR;
        ds->df0010 += factor*faZ;
        ds->df1010 += factor*faRZ;
        ds->df2000 += factor*faRR;
        ds->df0020 += factor*faZZ;
        ds->df3000 += factor*faRRR;
        ds->df2010 += factor*faRRZ;
        ds->df1020 += factor*faRZZ;
        ds->df0030 += factor*faZZZ;
        ds->df4000 += factor*faRRRR;
        ds->df3010 += factor*faRRRZ;
        ds->df2020 += factor*faRRZZ;
        ds->df1030 += factor*faRZZZ;
        ds->df0040 += factor*faZZZZ;
    }
    if (dp->rhob >G96_THRESHOLD) {
        rb = dp->rhob;
        grb = dp->gradb;
        rb13 = pow(dp->rhob,1.0/3.0);
        rb43 = rb13*rb;
        rbm23 = rb13/rb;
        rbm43 = 1.0/rb43;
        rb2 = rb*rb;
        rb3 = rb*rb*rb;
        xb = grb/rb43;
        xb2 = xb*xb;
        xb3 = xb*xb*xb;
        xb12 = rbm23*pow(dp->gradb,1.0/2.0);
        fbR = 2.0/3.0*rb13*GAMMA*xb12*xb;
        fbZ = -3.0/2.0*GAMMA*xb12;
        fbRR = -10.0/9.0*GAMMA*xb12*xb*rbm23;
        fbRZ = GAMMA*xb12/rb;
        fbZZ = -3.0/4.0*GAMMA*rbm43*xb12/xb;
        fbRRR = 80.0/27.0*GAMMA*xb12*xb*rbm23/rb;
        fbRRZ = -5.0/3.0*GAMMA*xb12/rb2;
        fbRZZ = 1.0/2.0*GAMMA*xb12/xb*rbm43/rb;
        fbZZZ = 3.0/8.0*GAMMA*xb12/xb2*rbm43*rbm43;
        fbRRRR = -880.0/81.0*GAMMA*rbm43*rbm43*xb12*xb;
        fbRRRZ = 40.0/9.0*GAMMA/rb3*xb12;
        fbRRZZ = -5.0/6.0*GAMMA*xb12/xb*rbm43/rb2;
        fbRZZZ = -1.0/4.0*GAMMA*rbm23/rb3*xb12/xb2;
        fbZZZZ = -9.0/16.0*GAMMA/(rb3*rb)*xb12/xb3;
        ds->df0100 += factor*fbR;
        ds->df0001 += factor*fbZ;
        ds->df0101 += factor*fbRZ;
        ds->df0200 += factor*fbRR;
        ds->df0002 += factor*fbZZ;
        ds->df0300 += factor*fbRRR;
        ds->df0201 += factor*fbRRZ;
        ds->df0102 += factor*fbRZZ;
        ds->df0003 += factor*fbZZZ;
        ds->df0400 += factor*fbRRRR;
        ds->df0301 += factor*fbRRRZ;
        ds->df0202 += factor*fbRRZZ;
        ds->df0103 += factor*fbRZZZ;
        ds->df0004 += factor*fbZZZZ;
    }
}
