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
/* fun-becke86.c:

   A. D. Becke, J. Chem. Phys. 107, 8554 (1997) 
   A. D. Becke,J. Chem. Phys. 84, 4524 (1986)

   implementation of Becke(86) functional and its derivatives.
   or exactly: Becke(86) GGA correction to the functional. (total Becke(86)
   energy is E_LDA+E_BCK).
   Implemented by David Wilson (davidwi@kjemi.uio.no), 2005.

   Equates with the B86R functional of MOLPRO. 

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
static int  becke86_isgga(void) { return 1; }
static int  becke86_read(const char* conf_line);
static real becke86_energy(const FunDensProp* dens_prop);
static void becke86_first(FunFirstFuncDrv *ds, real factor, 
                        const FunDensProp* dens_prop);
static void becke86_second(FunSecondFuncDrv *ds, real factor,
                         const FunDensProp* dens_prop);
static void becke86_third(FunThirdFuncDrv *ds, real factor,
                        const FunDensProp* dens_prop);
#ifdef FOURTH_ORDER_DERIVATIVE
static void becke86_fourth(FunFourthFuncDrv *ds, real factor,
                        const FunDensProp* dens_prop);
#endif

Functional B86xFunctional = {
    "B86x",      /* name */
    becke86_isgga,  /* gga-corrected */
    becke86_read,   /* set bloody common blocks */
    NULL,         /* reporter */
    becke86_energy, 
    becke86_first,
    becke86_second,
    becke86_third
#ifdef FOURTH_ORDER_DERIVATIVE
    ,becke86_fourth
#endif
};

/* IMPLEMENTATION PART */

static int
becke86_read(const char* conf_line)
{
    fun_set_hf_weight(0.0);
    return 1;
}

/* becke86_energy:
   note that in reality E_BCK86 = E_BCK86,alpha + E_BCK86,beta
   i.e the is linear in alpha and beta densities.

   Becke86 threshold is needed to avoid numerical problems on 0/0
   divisions.  The problems are small but it is better to be on the
   safe side.
*/
static const real BECKE_THRESHOLD = 1e-14;
static const real BETA = 0.00787;
static const real GAMMA = 0.004;
static real
becke86_energy(const FunDensProp* dp)
{
   real ea,eb;
   if (dp->rhoa<BECKE_THRESHOLD) 
     ea=0;
   else {
       real PREF = 3.0/4.0*pow(6/M_PI, 1.0/3.0); 
       real xa = dp->grada*pow(dp->rhoa,-4.0/3.0);
       real ra43 = pow(dp->rhoa,4.0/3.0);
       real denoma = 1.0 + GAMMA*xa*xa;
       ea = PREF*ra43*(1.0 + BETA*xa*xa)/denoma;
   }
   if (dp->rhob<BECKE_THRESHOLD)
     eb = 0.0;
   else {
       real PREF = 3.0/4.0*pow(6/M_PI, 1.0/3.0); 
       real xb = dp->gradb*pow(dp->rhob,-4.0/3.0);
       real rb43 = pow(dp->rhob,4.0/3.0);
       real denomb = 1.0 +GAMMA*xb*xb;
       eb = PREF*rb43*(1.0 + BETA*xb*xb)/denomb;
   } 
   return -(ea+eb);
}

static void
becke86_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real ra13, ra43, xa, xa2, denoma, denoma2;
    real rb13, rb43, xb, xb2, denomb, denomb2;
    real faR, faZ, t1a, t2a;
    real fbR, fbZ, t1b, t2b;
    real faR1, faR2, faZ1, faZ2, PREF;
    real fbR1, fbR2, fbZ1, fbZ2;

    if (dp->rhoa >BECKE_THRESHOLD) {
       PREF = 3.0/4.0*pow(6/M_PI, 1.0/3.0);
       xa = dp->grada*pow(dp->rhoa,-4.0/3.0);
       xa2 = xa*xa;
       ra43 = pow(dp->rhoa,4.0/3.0);
       ra13 = pow(dp->rhoa,1.0/3.0);
       denoma = 1.0 +GAMMA*xa2;
       denoma2 = denoma*denoma;
       t1a = 4.0/3.0*PREF*BETA*ra13*xa2/denoma;
       faR1 = 1.0 - 2.0*GAMMA*xa2/denoma; 
       faR2 = -4.0/3.0*PREF*ra13/denoma 
              - 8.0/3.0*PREF*GAMMA*ra13*xa*xa/denoma2;
       t2a = -2.0*PREF*BETA*xa/denoma;
       faZ1 = 1.0-GAMMA*xa2/denoma;
       faZ2 = 2.0*PREF*GAMMA*xa/denoma2;
       faR = faR1*t1a + faR2;
       faZ = faZ1*t2a + faZ2;
       ds->df1000 += factor*faR;
       ds->df0010 += factor*faZ;
    }
    if (dp->rhob >BECKE_THRESHOLD) {
       PREF = 3.0/4.0*pow(6/M_PI, 1.0/3.0);
       xb = dp->gradb*pow(dp->rhob,-4.0/3.0);
       xb2 = xb*xb;
       rb43 = pow(dp->rhob,4.0/3.0);
       rb13 = pow(dp->rhob,1.0/3.0);
       denomb = 1.0 +GAMMA*xb2;
       denomb2 = denomb*denomb;
       t1b = 4.0/3.0*PREF*BETA*rb13*xb2/denomb;
       fbR1 = 1.0 - 2.0*GAMMA*xb2/denomb; 
       fbR2 = -4.0/3.0*PREF*rb13/denomb 
              - 8.0/3.0*PREF*GAMMA*rb13*xb*xb/denomb2;
       t2b = -2.0*PREF*BETA*xb/denomb;
       fbZ1 = 1.0-GAMMA*xb2/denomb;
       fbZ2 = 2.0*PREF*GAMMA*xb/denomb2;
       fbR = t1b*fbR1 + fbR2;
       fbZ = t2b*fbZ1 + fbZ2;
       ds->df0100 += factor*fbR;
       ds->df0001 += factor*fbZ;
    } 
}

static void
becke86_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real ra, ra13, ra43, ram23, ram43, xa, xa2, denoma;
    real rb, rb13, rb43, rbm23, rbm43, xb, xb2, denomb;
    real gdxa, gdxa2, denoma2, PREF;
    real gdxb, gdxb2, denomb2;
    real faR, faR1, faR2, faZ, faZ1, faZ2, faRR, faRZ, faZZ;
    real fbR, fbR1, fbR2, fbZ, fbZ1, fbZ2, fbRR, fbRZ, fbZZ;
    real faRR1, faRR2, faRZ1, faRZ2, faZZ1, faZZ2;
    real fbRR1, fbRR2, fbRZ1, fbRZ2, fbZZ1, fbZZ2;
    real t1a, t2a, t3a, t4a;
    real t1b, t2b, t3b, t4b;
    
    if (dp->rhoa >BECKE_THRESHOLD) {    
       PREF = 3.0/4.0*pow(6/M_PI, 1.0/3.0);
       xa = dp->grada*pow(dp->rhoa,-4.0/3.0);
       xa2 = xa*xa;
       ra = dp->rhoa;
       ra43 = pow(dp->rhoa,4.0/3.0);
       ra13 = pow(dp->rhoa,1.0/3.0);
       ram23 = pow(dp->rhoa,-2.0/3.0);
       ram43 = 1.0/ra43;
       denoma = 1.0 +GAMMA*xa2;
       denoma2 = denoma*denoma;
       gdxa = GAMMA*xa2/denoma;
       gdxa2 = gdxa*gdxa;
       t1a = 4.0/3.0*PREF*BETA*ra13*xa2/denoma;
       faR1 = 1.0 - 2.0*gdxa;
       faR2 = -4.0/3.0*PREF*ra13/denoma - 8.0/3.0*PREF*ra13*gdxa/denoma;
       t2a = -2.0*PREF*BETA*xa/denoma;
       faZ1 = 1.0 - gdxa;
       faZ2 = 2.0*PREF*GAMMA*xa/denoma2;
       faR = faR1*t1a + faR2;
       faZ = faZ1*t2a + faZ2;

       t3a = 1.0 - 5.0*gdxa + 4.0*gdxa2;
       t4a = 7.0 - 38.0*gdxa + 32.0*gdxa2;
       faZZ1 = -2.0*PREF*BETA*ram43/denoma;
       faZZ2 = 2.0*PREF*GAMMA*ram43/denoma2*(1.0-4.0*gdxa);
       faRR1 = -4.0/9.0*PREF*BETA*ram23*xa2/denoma;
       faRR2 = -4.0/9.0*PREF*ram23/denoma*
                  (1.0 - 6.0*gdxa + 32.0*gdxa2);
       faRZ1 = 8.0/3.0*PREF*BETA*xa/(denoma*ra);
       faRZ2 = -8.0/3.0*PREF*GAMMA*xa/ra/denoma2*(1.0 - 4.0*gdxa);
       faRR = t4a*faRR1 + faRR2;
       faRZ = t3a*faRZ1 + faRZ2;
       faZZ = t3a*faZZ1 + faZZ2;
       ds->df1000 += factor*faR;
       ds->df0010 += factor*faZ;
       ds->df1010 += factor*faRZ;
       ds->df0020 += factor*faZZ;
       ds->df2000 += factor*faRR;

  }
  if (dp->rhob >BECKE_THRESHOLD) {    
       PREF = 3.0/4.0*pow(6/M_PI, 1.0/3.0); 
       xb = dp->gradb*pow(dp->rhob,-4.0/3.0);
       xb2 = xb*xb;
       rb = dp->rhob;
       rb43 = pow(dp->rhob,4.0/3.0);
       rb13 = pow(dp->rhob,1.0/3.0); 
       rbm43 = 1.0/rb43;
       rbm23 = rb13/rb;
       denomb = 1.0 +GAMMA*xb2;
       denomb2 = denomb*denomb;
       gdxb = GAMMA*xb2/denomb;
       gdxb2 = gdxb*gdxb;

       t1b = 4.0/3.0*PREF*BETA*rb13*xb2/denomb; 
       fbR1 = 1.0 - 2.0*gdxb;
       fbR2 = -4.0/3.0*PREF*rb13/denomb
              - 8.0/3.0*PREF*rb13*gdxb/denomb;
       t2b = -2.0*PREF*BETA*xb/denomb;
       fbZ1 = 1.0 - gdxb;
       fbZ2 = 2.0*PREF*GAMMA*xb/denomb2;
       fbR = t1b*fbR1 + fbR2; 
       fbZ = t2b*fbZ1 + fbZ2;

       t3b = 1.0 - 5.0*gdxb + 4.0*gdxb2;
       t4b = 7.0 - 38.0*gdxb + 32.0*gdxb2;
       fbZZ1 = -2.0*PREF*BETA*rbm43/denomb;
       fbRR1 = -4.0/9.0*PREF*BETA*rbm23*xb2/denomb;
       fbRZ1 = 8.0/3.0*PREF*BETA*xb/(denomb*rb);
       fbZZ2 = 2.0*PREF*GAMMA*rbm43/denomb2*(1.0-4.0*gdxb);
       fbRR2 = -4.0/9.0*PREF*rbm23/denomb*
                  (1.0 - 6.0*gdxb + 32.0*gdxb2);
       fbRZ2 = -8.0/3.0*PREF*GAMMA*xb/rb/denomb2*(1.0 - 4.0*gdxb);
       fbRR = t4b*fbRR1 + fbRR2;
       fbRZ = t3b*fbRZ1 + fbRZ2;
       fbZZ = t3b*fbZZ1 + fbZZ2;
       ds->df0100 += factor*fbR;
       ds->df0001 += factor*fbZ;
       ds->df0101 += factor*fbRZ;
       ds->df0002 += factor*fbZZ;
       ds->df0200 += factor*fbRR;
    }  
}


static void
becke86_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real ra, ra2, ra13, ra43, ram23, ram43, ram53, ram73, xa, xa2;
    real rb, rb2, rb13, rb43, rbm23, rbm43, rbm53, rbm73, xb, xb2;
    real gdxa, gdxa2, gdxa3, denoma, denoma2, denoma3, PREF;
    real gdxb, gdxb2, gdxb3, denomb, denomb2, denomb3;
    real faR, faZ, faRR, faRZ, faZZ, faRRR, faRRZ, faRZZ, faZZZ;
    real fbR, fbZ, fbRR, fbRZ, fbZZ, fbRRR, fbRRZ, fbRZZ, fbZZZ;
    real faR1, faR2, faZ1, faZ2, faRR1, faRR2, faRZ1, faRZ2, faZZ1, faZZ2;
    real fbR1, fbR2, fbZ1, fbZ2, fbRR1, fbRR2, fbRZ1, fbRZ2, fbZZ1, fbZZ2;
    real faRRR1, faRRR2, faRRZ1, faRRZ2, faRZZ1, faRZZ2, faZZZ1, faZZZ2;
    real fbRRR1, fbRRR2, fbRRZ1, fbRRZ2, fbRZZ1, fbRZZ2, fbZZZ1, fbZZZ2;
    real t1a, t2a, t3a, t4a, t5a, t6a, t7a, t8a;
    real t1b, t2b, t3b, t4b, t5b, t6b, t7b, t8b;
   
    if (dp->rhoa >BECKE_THRESHOLD) {
       real PREF = 3.0/4.0*pow(6/M_PI, 1.0/3.0);
       xa = dp->grada*pow(dp->rhoa,-4.0/3.0);
       xa2 = xa*xa;
       ra = dp->rhoa;
       ra2 = ra*ra;
       ra13 = pow(dp->rhoa,1.0/3.0);
       ra43 = pow(dp->rhoa,4.0/3.0);
       ram23 = ra13/ra;
       ram43 = 1.0/ra43;
       ram53 = ram23/ra;
       ram73 = ram43/ra;
       denoma = 1.0 + GAMMA*xa2;
       denoma2 = denoma*denoma;
       denoma3 = denoma2*denoma;
       gdxa = GAMMA*xa2/denoma;
       gdxa2 = gdxa*gdxa;
       gdxa3 = gdxa2*gdxa;

       t1a = 4.0/3.0*PREF*BETA*ra13*xa2/denoma;
       faR1 = 1.0 - 2.0*gdxa;
       faR2 = -4.0/3.0*PREF*ra13/denoma - 8.0/3.0*PREF*ra13*gdxa/denoma;
       t2a = -2.0*PREF*BETA*xa/denoma;
       faZ1 = 1.0 - gdxa; 
       faZ2 = 2.0*PREF*GAMMA*xa/denoma2;
       faR = faR1*t1a + faR2;
       faZ = faZ1*t2a + faZ2;
       
       t3a = 1.0 - 5.0*gdxa + 4.0*gdxa2;
       t4a = 7.0 - 38.0*gdxa + 32.0*gdxa2;
       faZZ1 = -2.0*PREF*BETA*ram43/denoma;
       faZZ2 = 2.0*PREF*GAMMA*ram43/denoma2*(1.0-4.0*gdxa);
       faRR1 = -4.0/9.0*PREF*BETA*ram23*xa2/denoma;
       faRR2 = -4.0/9.0*PREF*ram23/denoma*
                  (1.0 - 6.0*gdxa + 32.0*gdxa2);
       faRZ1 = 8.0/3.0*PREF*BETA*xa/(denoma*ra);
       faRZ2 = -8.0/3.0*PREF*GAMMA*xa/ra/denoma2*(1.0 - 4.0*gdxa);
       faRR = t4a*faRR1 + faRR2;
       faRZ = t3a*faRZ1 + faRZ2;
       faZZ = t3a*faZZ1 + faZZ2;

       t5a = 8.0/27.0*PREF*BETA*ram53*xa2/denoma;
       t6a = 8.0/9.0*PREF*BETA*xa/(denoma*ra2);
       t7a = 8.0/3.0*PREF*BETA*ram73/denoma;
       t8a = 24.0*PREF*BETA*ram53/ra*gdxa/(denoma*xa);
       faRRR1 = 35.0 - 370.0*gdxa + 720.0*gdxa2 - 384.0*gdxa3;
       faRRR2 = 8.0/27.0*PREF/(ra43*ra13*denoma)*
               (1.0 - 34.0*gdxa + 336.0*gdxa2 - 384.0*gdxa3);
       faRRZ1 = -7.0 + 83.0*gdxa - 172.0*gdxa2 + 96.0*gdxa3;
       faRRZ2 = 8.0/9.0*PREF*GAMMA/ra/ra*xa/denoma2*
                (7.0 - 76.0*gdxa + 96.0*gdxa2);
       faRZZ1 = 1.0 - 17.0*gdxa + 40.0*gdxa2 - 24.0*gdxa3;
       faRZZ2 = -8.0/3.0*PREF*GAMMA*ram43/ra/denoma2*
                 (1.0 - 16.0*gdxa + 24.0*gdxa2);
       faZZZ1 = 1.0 - 3.0*gdxa + 2.0*gdxa2;
       faZZZ2 = -24.0*PREF*GAMMA*GAMMA*xa/(denoma3*ra43*ra43)*
                 (1.0 - 2.0*gdxa);

       faRRR = faRRR1*t5a + faRRR2;
       faRRZ = faRRZ1*t6a + faRRZ2;
       faRZZ = faRZZ1*t7a + faRZZ2;
       faZZZ = faZZZ1*t8a + faZZZ2;
       ds->df1000 += factor*faR;
       ds->df0010 += factor*faZ;
       ds->df1010 += factor*faRZ;
       ds->df0020 += factor*faZZ; 
       ds->df2000 += factor*faRR; 
       ds->df3000 += factor*faRRR;
       ds->df2010 += factor*faRRZ;
       ds->df1020 += factor*faRZZ;
       ds->df0030 += factor*faZZZ;
  }
  if (dp->rhob >BECKE_THRESHOLD) {
       real PREF = 3.0/4.0*pow(6/M_PI, 1.0/3.0);
       xb = dp->gradb*pow(dp->rhob,-4.0/3.0);
       xb2 = xb*xb;
       rb = dp->rhob;
       rb2 = rb*rb;
       rb13 = pow(dp->rhob,1.0/3.0);
       rb43 = pow(dp->rhob,4.0/3.0);
       rbm23 = rb13/rb;
       rbm43 = 1.0/rb43;
       rbm53 = rbm23/rb;
       rbm73 = rbm43/rb;
       denomb = 1.0 +GAMMA*xb2;
       denomb2 = denomb*denomb;
       denomb3 = denomb2*denomb;
       gdxb = GAMMA*xb2/denomb;
       gdxb2 = gdxb*gdxb;
       gdxb3 = gdxb2*gdxb;

       t1b = 4.0/3.0*PREF*BETA*rb13*xb2/denomb;
       fbR1 = 1.0 - 2.0*gdxb;
       fbR2 = -4.0/3.0*PREF*rb13/denomb
              - 8.0/3.0*PREF*rb13*gdxb/denomb;
       t2b = -2.0*PREF*BETA*xb/denomb;
       fbZ1 = 1.0 - gdxb;
       fbZ2 = 2.0*PREF*GAMMA*xb/denomb2;
       fbR = t1b*fbR1 + fbR2;
       fbZ = t2b*fbZ1 + fbZ2;

       t3b = 1.0 - 5.0*gdxb + 4.0*gdxb2;
       t4b = 7.0 - 38.0*gdxb + 32.0*gdxb2;
       fbZZ1 = -2.0*PREF*BETA*rbm43/denomb;
       fbRR1 = -4.0/9.0*PREF*BETA*rbm23*xb2/denomb;
       fbRZ1 = 8.0/3.0*PREF*BETA*xb/(denomb*rb);
       fbZZ2 = 2.0*PREF*GAMMA*rbm43/denomb2*(1.0-4.0*gdxb);
       fbRR2 = -4.0/9.0*PREF*rbm23/denomb*
                  (1.0 - 6.0*gdxb + 32.0*gdxb2);
       fbRZ2 = -8.0/3.0*PREF*GAMMA*xb/rb/denomb2*(1.0 - 4.0*gdxb);
       fbRR = t4b*fbRR1 + fbRR2;
       fbRZ = t3b*fbRZ1 + fbRZ2;
       fbZZ = t3b*fbZZ1 + fbZZ2;

       t5b = 8.0/27.0*PREF*BETA*rbm53*xb2/denomb;
       t6b = 8.0/9.0*PREF*BETA*xb/(denomb*rb2);
       t7b = 8.0/3.0*PREF*BETA*rbm73/denomb;
       t8b = 24.0*PREF*BETA*rbm53/rb*gdxb/(denomb*xb);
       fbRRR1 = 35.0 - 370.0*gdxb + 720.0*gdxb2 - 384.0*gdxb3;
       fbRRR2 = 8.0/27.0*PREF*rbm43/rb13/denomb*
               (1.0 - 34.0*gdxb + 336.0*gdxb2 - 384.0*gdxb3);
       fbRRZ1 = -7.0 + 83.0*gdxb - 172.0*gdxb2 + 96.0*gdxb3;
       fbRRZ2 = 8.0/9.0*PREF*GAMMA/rb/rb*xb/denomb2*
                (7.0 - 76.0*gdxb + 96.0*gdxb2);
       fbRZZ1 = 1.0 - 17.0*gdxb + 40.0*gdxb2 - 24.0*gdxb3;
       fbRZZ2 = -8.0/3.0*PREF*GAMMA*rbm43/rb/denomb2*
                 (1.0 - 16.0*gdxb + 24.0*gdxb2);
       fbZZZ1 = 1.0 - 3.0*gdxb + 2.0*gdxb2; 
       fbZZZ2 = -24.0*PREF*GAMMA*GAMMA*xb/(denomb3*rb43*rb43)*
                 (1.0 - 2.0*gdxb);
       fbRRR = fbRRR1*t5b + fbRRR2;
       fbRRZ = fbRRZ1*t6b + fbRRZ2;
       fbRZZ = fbRZZ1*t7b + fbRZZ2;
       fbZZZ = fbZZZ1*t8b + fbZZZ2;

       ds->df0100 += factor*fbR;
       ds->df0001 += factor*fbZ;
       ds->df0101 += factor*fbRZ;
       ds->df0002 += factor*fbZZ;
       ds->df0200 += factor*fbRR;
       ds->df0300 += factor*fbRRR;
       ds->df0201 += factor*fbRRZ;
       ds->df0102 += factor*fbRZZ;
       ds->df0003 += factor*fbZZZ;
    }
}


