/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* fun-kt.c:
   implementation of KT functional and its derivatives.
   or exactly: KT GGA correction to the functional for KT1,KT2 total functional
   energy is E_LDA+E_KT).
   Reference:  Keal, Tozer, J. Chem. Phys., 119, 3015 (2003).
   GAMMA is included in the KTx definition in fun-gga.c 
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
static int  kt_isgga(void) { return 1; }
static int  kt_read(const char* conf_line);
static real kt_energy(const DftDensProp* dens_prop);
static void kt_first(FirstFuncDrv *ds, real factor, 
                        const DftDensProp* dens_prop);
static void kt_second(SecondFuncDrv *ds, real factor, 
                        const DftDensProp* dens_prop);
static void kt_third(ThirdFuncDrv *ds, real factor, 
                        const DftDensProp* dens_prop);

Functional KTFunctional = {
    "KT",      /* name */
    kt_isgga,  /* gga-corrected */
    kt_read,   /* set bloody common blocks */
    NULL,         /* reporter */
    kt_energy, 
    kt_first, 
    kt_second,
    kt_third
};

/* IMPLEMENTATION PART */

static int
kt_read(const char* conf_line)
{
    dft_set_hf_weight(0.0);
    return 1;
}

/* kt_energy:
   note that in reality E_KT = E_KT,alpha + E_KT,beta
   i.e the energy is linear in alpha and beta densities.

   KT threshold is needed to avoid numerical problems on 0/0
   divisions.  The problems are small but it is better to be on the
   safe side.
*/
static const real KT_THRESHOLD = 1e-14;
static const real DELTA = 0.1;
static real
kt_energy(const DftDensProp* dp)
{
   real ea,eb;
   if (dp->rhob<KT_THRESHOLD)
     eb = 0.0;
   else {
     real xb = dp->gradb;
     real rb = pow(dp->rhob,4.0/3.0);
     real denomb = rb + DELTA;
     eb = xb*xb/denomb; 
   } 
   if (dp->rhoa<KT_THRESHOLD) 
     ea=0;
   else {
     real xa = dp->grada;
     real ra = pow(dp->rhoa,4.0/3.0);
     real denoma = ra + DELTA;
     ea = xa*xa/denoma;   
   }
   return (ea+eb);
}

static void
kt_first(FirstFuncDrv *ds, real factor, const DftDensProp* dp)
{
    real xa, ra43, ra13;
    real xb, rb43, rb13;
    real denoma, denoma2;
    real denomb, denomb2; 
    real faZ, faR;
    real fbZ, fbR;

    if (dp->rhoa >KT_THRESHOLD) {
        xa = dp->grada;
        ra43 = pow(dp->rhoa,4.0/3.0);
        ra13 = pow(dp->rhoa,1.0/3.0);
        denoma = ra43 + DELTA;
        denoma2= denoma*denoma;
        faR = -4.0/3.0*xa*xa*ra13/denoma2;
        faZ = 2.0*xa/denoma;
        ds->df1000 += factor*faR;
        ds->df0010 += factor*faZ;
    }
    if (dp->rhob >KT_THRESHOLD) {
        xb = dp->gradb;
        rb43 = pow(dp->rhob,4.0/3.0);
        rb13 = pow(dp->rhob,1.0/3.0);
        denomb = rb43 + DELTA;
        denomb2= denomb*denomb;
        fbR = -4.0/3.0*xb*xb*rb13/denomb2; 
        fbZ = 2.0*xb/denomb;
        ds->df0100 += factor*fbR;
        ds->df0001 += factor*fbZ;
    } 
}


static void
kt_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp)
{
    real xa, ra43, ra13, ram2;
    real xb, rb43, rb13, rbm2;
    real denoma, denoma2, t1a;
    real denomb, denomb2, t1b; 
    real faR, faZ, faRZ, faZZ, faRR;
    real fbR, fbZ, fbRZ, fbZZ, fbRR;

    if (dp->rhoa >KT_THRESHOLD) {
        xa = dp->grada;
        ra43 = pow(dp->rhoa,4.0/3.0);
        ra13 = pow(dp->rhoa,1.0/3.0);
        ram2 = pow(dp->rhoa,-2.0/3.0);
        denoma = ra43 + DELTA;
        denoma2= denoma*denoma;
        faR = -(4.0/3.0)*xa*xa*ra13/denoma2;
        faZ = 2.0*xa/denoma; 
	t1a = 4.0/9.0*xa*xa/denoma2;
	faRR = -ram2 + 8.0*ra13*ra13/denoma;
        faRZ  = -8.0/3.0*xa*ra13/denoma2;
        faZZ  = 2.0/denoma;
	ds->df1000 += factor*faR;
	ds->df0010 += factor*faZ; 
	ds->df1010 += factor*faRZ;
	ds->df2000 += factor*faRR*t1a;
	ds->df0020 += factor*faZZ;
    }
    if (dp->rhob >KT_THRESHOLD) {
        xb = dp->gradb;
        rb43 = pow(dp->rhob,4.0/3.0);
        rb13 = pow(dp->rhob,1.0/3.0);
        rbm2 = pow(dp->rhob,-2.0/3.0);
        denomb = rb43 + DELTA;
        denomb2= denomb*denomb;
        fbR = -4.0/3.0*xb*xb*rb13/denomb2; 
        fbZ = 2.0*xb/denomb; 
	t1b = 4.0/9.0*xb*xb/denomb2;
	fbRR = -rbm2 + 8.0*rb13*rb13/denomb;
	fbRZ  = -8.0/3.0*xb*rb13/denomb2;
	fbZZ  = 2.0/denomb;
	ds->df0100 += factor*fbR;
	ds->df0001 += factor*fbZ;
	ds->df0101 += factor*fbRZ;
	ds->df0200 += factor*fbRR*t1b;
	ds->df0002 += factor*fbZZ;
    } 
}


static void
kt_third(ThirdFuncDrv *ds, real factor, const DftDensProp* dp)
{
    real ra, xa, ra43, ra13, ra23, ram1, ram2;
    real rb, xb, rb43, rb13, rb23, rbm1, rbm2;
    real xad2, t1a, t2a;
    real xbd2, t1b, t2b;
    real denoma, denoma2;
    real denomb, denomb2; 
    real faR, faZ, faRR, faZZ, faRZ; 
    real fbR, fbZ, fbRR, fbZZ, fbRZ;
    real faRRR, faRRZ, faRZZ, faZZZ;
    real fbRRR, fbRRZ, fbRZZ, fbZZZ;

    if (dp->rhoa >KT_THRESHOLD) {
        xa = dp->grada;
        ra = dp->rhoa;
        ra13 = pow(dp->rhoa,1.0/3.0);
        ra43 = ra13*ra;
        ra23 = ra13*ra13;
        ram2 = ra13/ra;
        ram1 = ram2*ra13;
        denoma = ra43 + DELTA;
        denoma2= denoma*denoma;
	xad2 = xa/denoma2;
        t1a = 8.0/9.0*xad2*xa*ram1;
        t2a = 1.0 - 8.0*ra43/denoma;
        faR  = -4.0/3.0*xad2*xa*ra13;
        faZ  = 2.0*xa/denoma; 
        faRR = -4.0/9.0*xad2*xa*ram2;
        faRZ = -8.0/3.0*xad2*ra13;
        faZZ = 2.0/denoma;
	faZZZ = 0.0;
	faRZZ = -8.0/3.0*ra13/denoma2;
	faRRZ =  -8.0/9.0*xad2*ram2;
	faRZZ = -8.0/3.0*ra13/denoma2;
        faRRR = 1.0/(3.0*ra43) + 4.0/denoma -
                16.0*ra43/denoma2;
	ds->df1000 += factor*faR;
	ds->df0010 += factor*faZ;
	ds->df1010 += factor*faRZ;
	ds->df2000 += factor*faRR*t2a;
	ds->df0020 += factor*faZZ;
	ds->df3000 += factor*faRRR*t1a;
	ds->df0030 += factor*faZZZ;
	ds->df2010 += factor*faRRZ*t2a;
	ds->df1020 += factor*faRZZ;
    }
    if (dp->rhob >KT_THRESHOLD) {
	xb = dp->gradb;
	rb = dp->rhob;
	rb13 = pow(dp->rhob,1.0/3.0);
        rb43 = rb13*rb;
	rb23 = rb13*rb13;
	rbm2 = rb13/rb;
	rbm1 = rbm2*rb13;
	denomb = rb43 + DELTA;
	denomb2= denomb*denomb;
	xbd2 = xb/denomb2;
        t1b = 8.0/9.0*xbd2*xb*rbm1;
        t2b = 1.0 - 8.0*rb43/denomb;
	fbR  = -4.0/3.0*xbd2*xb*rb13;
	fbZ  = 2.0*xb/denomb; 
        fbRR = -4.0/9.0*xbd2*xb*rbm2;
	fbRZ = -8.0/3.0*xbd2*rb13;
	fbZZ = 2.0/denomb;
	fbZZZ = 0.0;
	fbRRZ = -8.0/9.0*xbd2*rbm2;
	fbRZZ = -8.0/3.0*rb13/denomb2;
	fbRRR = 1.0/(3.0*rb43) +  4.0/denomb -
		16.0*rb43/denomb2;
	ds->df0100 += factor*fbR;
	ds->df0001 += factor*fbZ;
	ds->df0101 += factor*fbRZ;
	ds->df0200 += factor*fbRR*t2b;
	ds->df0002 += factor*fbZZ;
	ds->df0300 += factor*fbRRR*t1b;
	ds->df0003 += factor*fbZZZ;
	ds->df0201 += factor*fbRRZ*t2b;
	ds->df0102 += factor*fbRZZ;
    } 
}




