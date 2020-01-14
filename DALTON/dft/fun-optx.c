/*


!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!

!

*/
/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* fun-optx.c:
   implementation of OPTX exchange functional and its derivatives.
   #### this is just the gradient corrected term  for KT3 functional#### 
   Reference: N.C. Handy and A.J. Cohen, Mol. Phys., 99, 403 (2001).
              Keal, Tozer, in press (2004).
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
#include "general.h"

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static integer  optx_isgga(void) { return 1; }
static integer  optxread(const char* conf_line);
static real optx_energy(const FunDensProp* dens_prop);
static void optx_first(FunFirstFuncDrv *ds, real factor, 
                        const FunDensProp* dens_prop);
static void optx_second(FunSecondFuncDrv *ds, real factor, 
                        const FunDensProp* dens_prop);
static void optx_third(FunThirdFuncDrv *ds, real factor, 
                        const FunDensProp* dens_prop);

Functional OPTXFunctional = {
    "OPTX_",      /* name */
    optx_isgga,  /* gga-corrected */
   1,
    optxread,   /* set bloody common blocks */
    NULL,         /* reporter */
    optx_energy, 
    optx_first,
    optx_second,
    optx_third

};

/* IMPLEMENTATION PART */

static integer
optxread(const char* conf_line)
{
    fun_set_hf_weight(0.0);
    return 1;
}

/* optx_energy:
   note that in reality E_OPTX = E_OPTX,alpha + E_OPTX,beta
   i.e the energy is linear in alpha and beta densities.

   OPTX threshold is needed to avoid numerical problems on 0/0
   divisions.  The problems are small but it is better to be on the
   safe side.
*/
static const real OPTX_THRESHOLD = 1e-14;
static const real GAMMA  = 0.006;
static real
optx_energy(const FunDensProp* dp)
{
   real ea,eb;
   if (dp->rhob<OPTX_THRESHOLD)
     eb = 0.0;
   else {
     real rb43 = pow(dp->rhob,4.0/3.0);
     real grb = dp->gradb;
     real xb = grb/rb43;
     real gxb2 = GAMMA*xb*xb;
     real ub = gxb2/(1.0 + gxb2); 
     eb = rb43*ub*ub; 
   } 
   if (dp->rhoa<OPTX_THRESHOLD) 
     ea=0;
   else {
     real ra43 = pow(dp->rhoa,4.0/3.0);
     real xa = (dp->grada)/ra43;
     real gxa2 = GAMMA*xa*xa;
     real ua = gxa2/(1.0 + gxa2);
     ea = ra43*ua*ua;
   }
   return (ea+eb);
}

static void
optx_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real ra, gra, xa, ra43, ra13, ua, ua2, gxa2;
    real rb, grb, xb, rb43, rb13, ub, ub2, gxb2;
    real faR, faZ;
    real fbR, fbZ;

    if (dp->rhoa >OPTX_THRESHOLD) {
        ra = dp->rhoa;
	gra = dp->grada;
        ra43 = pow(dp->rhoa,4.0/3.0);
        ra13 = pow(dp->rhoa,1.0/3.0);
        xa = gra/ra43;
        gxa2 = GAMMA*xa*xa;
        ua = gxa2/(1.0 + gxa2);
        ua2 = ua*ua;
        faZ = 4.0*ua2*(1.0 - ua)/xa;
        faR = 4.0*ua2*ra13*(4.0/3.0*ua - 1.0);
        ds->df1000 += factor*faR;
        ds->df0010 += factor*faZ;
    }
    if (dp->rhob >OPTX_THRESHOLD) {
        rb = dp->rhob;
        grb = dp->gradb;
        rb43 = pow(dp->rhob,4.0/3.0);
        rb13 = pow(dp->rhob,1.0/3.0);
        xb = grb/rb43;
        gxb2 = GAMMA*xb*xb;
        ub = gxb2/(1.0 + gxb2);
        ub2 = ub*ub;
        fbZ = 4.0*ub2*(1.0 - ub)/xb;
        fbR = 4.0*ub2*rb13*(4.0/3.0*ub - 1.0); 
        ds->df0100 += factor*fbR;
        ds->df0001 += factor*fbZ;
    } 
}


static void
optx_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real ra, gra, xa, ra43, ra13, ua, ua2, gxa2;
    real rb, grb, xb, rb43, rb13, ub, ub2, gxb2;
    real faR, faZ, faZZ, faRZ, faRR;
    real fbR, fbZ, fbZZ, fbRZ, fbRR;
    real faca, fac2a, facb, fac2b;

    if (dp->rhoa >OPTX_THRESHOLD) {
        ra = dp->rhoa;
        ra43 = pow(dp->rhoa,4.0/3.0);
        ra13 = pow(dp->rhoa,1.0/3.0);
        gra = dp->grada;
        xa = (dp->grada)/ra43;
        gxa2 = GAMMA*xa*xa;
        ua = gxa2/(1.0 + gxa2);
        ua2 = ua*ua;
        faZ = 4.0*ua2*(1.0 - ua)/xa;
	faR = 4.0/3.0*ra13*ua2*(4.0*ua - 3.0);
	faca = (1.0 - 2.0*ua)*(1.0 - ua);
        faZZ = 12.0*ua2/(xa*xa*ra43);
        faRZ = -16.0*ua2/(xa*ra);
	fac2a = 4.0*ua2*ra13/(3.0*ra);
	faRR = 4.0/3.0*ua - 1.0 + 16.0*faca;
        ds->df1000 += factor*faR;
        ds->df0010 += factor*faZ;
        ds->df2000 += factor*faRR*fac2a;
        ds->df0020 += factor*faZZ*faca;
        ds->df1010 += factor*faRZ*faca;
    }
    if (dp->rhob >OPTX_THRESHOLD) {
        rb = dp->rhob;
        rb43 = pow(dp->rhob,4.0/3.0);
        rb13 = pow(dp->rhob,1.0/3.0);
        grb = dp->gradb;
        xb = (dp->gradb)/rb43;
        gxb2 = GAMMA*xb*xb;
        ub = gxb2/(1.0 + gxb2);
        ub2 = ub*ub;
        fbZ = 4.0*ub2*(1.0 - ub)/xb;
	fbR = 4.0/3.0*rb13*ub2*(4.0*ub - 3.0);
	facb = (1.0 - 2.0*ub)*(1.0 - ub);
        fbZZ = 12.0*ub2/(xb*xb*rb43);
        fbRZ = -16.0*ub2/(xb*rb);
	fac2b = 4.0*ub2*rb13/(3.0*rb);
	fbRR = 4.0/3.0*ub - 1.0 + 16.0*facb;
        ds->df0100 += factor*fbR;
        ds->df0001 += factor*fbZ;
        ds->df0200 += factor*fbRR*fac2b;
        ds->df0002 += factor*fbZZ*facb;
        ds->df0101 += factor*fbRZ*facb;

    } 
}

static void
optx_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real ra, gra, xa, xa2, ra43, ra13, ua, ua2, gxa2;
    real rb, grb, xb, xb2, rb43, rb13, ub, ub2, gxb2;
    real faR, faZ, faZZ, faRZ, faRR;
    real fbR, fbZ, fbZZ, fbRZ, fbRR;
    real faRRR, faRRZ, faRZZ, faZZZ;
    real fbRRR, fbRRZ, fbRZZ, fbZZZ;
    real faca, fac2a, t3a, t4a, t5a;
    real facb, fac2b, t3b, t4b, t5b;

    if (dp->rhoa >OPTX_THRESHOLD) {
        ra = dp->rhoa;
        ra43 = pow(dp->rhoa,4.0/3.0);
        ra13 = pow(dp->rhoa,1.0/3.0);
        gra = dp->grada;
        xa = (dp->grada)/ra43;
	xa2 = xa*xa;
        gxa2 = GAMMA*xa2;
        ua = gxa2/(1.0 + gxa2);
        ua2 = ua*ua;
        faZ = 4.0*ua2*(1.0 - ua)/xa;
	faR = 4.0/3.0*ra13*ua2*(4.0*ua - 3.0);
	faca = (1.0 - 2.0*ua)*(1.0 - ua);
        faZZ = 12.0*ua2/(xa*xa*ra43);
        faRZ = -16.0*ua2/(xa*ra);
	fac2a = 4.0*ua2*ra13/(3.0*ra);
	faRR = 4.0/3.0*ua - 1.0 + 16.0*faca;
	t3a = 4.0*ua2*(1.0 - ua)/(xa2*ra43);
	t4a = 16.0*ua2*(1.0 - ua)/(3.0*xa*ra*ra);
	t5a = -8.0*ua2*ra13/(9.0*ra*ra);
	faZZZ = 6.0/(xa*ra43)*(8.0*ua2 - 7.0*ua + 1.0);
	faRZZ = -(4.0/ra)*(4.0*ua - 3.0)*(4.0*ua - 1.0);
	faRRZ = 64.0*ua2 - 70.0*ua + 15.0;
	faRRR = 4.0*(1.0-ua)*(128.0*ua2-140.0*ua+30.0) +
	        32.0*ua2 - 140.0/3.0*ua + 15;
        ds->df1000 += factor*faR;
        ds->df0010 += factor*faZ;
        ds->df2000 += factor*faRR*fac2a;
        ds->df0020 += factor*faZZ*faca;
        ds->df1010 += factor*faRZ*faca;
        ds->df0030 += factor*faZZZ*t3a;
        ds->df1020 += factor*faRZZ*t3a;
        ds->df2010 += factor*faRRZ*t4a;
        ds->df3000 += factor*faRRR*t5a;
    }
    if (dp->rhob >OPTX_THRESHOLD) {
        rb = dp->rhob;
        rb43 = pow(dp->rhob,4.0/3.0);
        rb13 = pow(dp->rhob,1.0/3.0);
        grb = dp->gradb;
        xb = (dp->gradb)/rb43;
	xb2 = xb*xb;
        gxb2 = GAMMA*xb2;
        ub = gxb2/(1.0 + gxb2);
        ub2 = ub*ub;
        fbZ = 4.0*ub2*(1.0 - ub)/xb;
	fbR = 4.0/3.0*rb13*ub2*(4.0*ub - 3.0);
	facb = (1.0 - 2.0*ub)*(1.0 - ub);
        fbZZ = 12.0*ub2/(xb*xb*rb43);
        fbRZ = -16.0*ub2/(xb*rb);
	fac2b = 4.0*ub2*rb13/(3.0*rb);
	fbRR = 4.0/3.0*ub - 1.0 + 16.0*facb;
	t3b = 4.0*ub2*(1.0 - ub)/(xb2*rb43);
	t4b = 16.0*ub2*(1.0 - ub)/(3.0*xb*rb*rb);
	t5b = -8.0*ub2*rb13/(9.0*rb*rb);
	fbZZZ = 6.0/(xb*rb43)*(8.0*ub2 - 7.0*ub + 1.0);
	fbRZZ = -(4.0/rb)*(4.0*ub - 3.0)*(4.0*ub - 1.0);
	fbRRZ = 64.0*ub2 - 70.0*ub + 15.0;
	fbRRR = 4.0*(1.0-ub)*(128.0*ub2-140.0*ub+30.0) +
		32.0*ub2 - 140.0/3.0*ub + 15;
        ds->df0100 += factor*fbR;
        ds->df0001 += factor*fbZ;
        ds->df0200 += factor*fbRR*fac2b;
        ds->df0002 += factor*fbZZ*facb;
        ds->df0101 += factor*fbRZ*facb;
        ds->df0003 += factor*fbZZZ*t3b;
        ds->df0102 += factor*fbRZZ*t3b;
        ds->df0201 += factor*fbRRZ*t4b;
        ds->df0300 += factor*fbRRR*t5b;

    } 
}


