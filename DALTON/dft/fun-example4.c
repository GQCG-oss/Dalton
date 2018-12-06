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
/* fun-example4.c:
   implementation of Example4 functional and its derivatives 
   (c) Pawel Salek, pawsa@theochem.kth.se, aug 2001
   NOTE:
   this file may seem unnecessarily complex but the structure really pays off
   when implementing multiple functionals depending on different parameters.
*/
/* strictly conform to XOPEN ANSI C standard */
#define __USE_XOPEN 

#include <math.h>
#include <stdio.h>
#include "general.h"

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static integer example4_isgga(void) { return 1; }
static integer example4_read(const char* conf_line);
static real example4_energy(const FunDensProp* dp);
static void example4_first(FunFirstFuncDrv *ds,   real factor, const FunDensProp* dp);
static void example4_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp);
static void example4_third(FunThirdFuncDrv *ds,   real factor, const FunDensProp* dp);

Functional Example4Functional = {
  "Example4",       /* name */
  example4_isgga,   /* gga-corrected */
   1,
  example4_read, 
  NULL,
  example4_energy, 
  example4_first,
  example4_second,
  example4_third
};

/* IMPLEMENTATION PART */
static integer
example4_read(const char* conf_line)
{
    fun_set_hf_weight(0);
    return 1;
}

static const real PREF= -5e-5;

static real
example4_energy(const FunDensProp* dp)
{
  real grad2 = dp->grada*dp->grada+dp->gradb*dp->gradb+2.0*dp->gradab;
  return PREF*dp->rhoa*dp->rhob*grad2;
}

static void
example4_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
  real grada2 = dp->grada*dp->grada;
  real gradb2 = dp->gradb*dp->gradb;    
  real grad2 = grada2+gradb2+2.0*dp->gradab; 
  real rhoab = dp->rhoa*dp->rhob;
  ds->df1000  += PREF*factor*dp->rhob*grad2;
  ds->df0100  += PREF*factor*dp->rhoa*grad2;
  ds->df0010  += 2.0*PREF*factor*rhoab*dp->grada;
  ds->df0001  += 2.0*PREF*factor*rhoab*dp->gradb;
  ds->df00001 += 2.0*PREF*factor*rhoab; 
}
static void
example4_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
  real grada2 = dp->grada*dp->grada;
  real gradb2 = dp->gradb*dp->gradb;    
  real grad2 = grada2+gradb2+2.0*dp->gradab; 
  real rhoab = dp->rhoa*dp->rhob;
  /* first derivatives */
  ds->df1000  += PREF*factor*dp->rhob*grad2;
  ds->df0100  += PREF*factor*dp->rhoa*grad2;
  ds->df0010  += 2.0*PREF*factor*rhoab*dp->grada;
  ds->df0001  += 2.0*PREF*factor*rhoab*dp->gradb;
  ds->df00001 += 2.0*PREF*factor*rhoab; 
  /* second derivatives */
  ds->df0020  += 2.0*PREF*factor*rhoab;
  ds->df0002  += 2.0*PREF*factor*rhoab; 
  /* mixed derivatives */
  ds->df1100  += PREF*factor*grad2;
  ds->df1010  += 2.0*PREF*factor*dp->rhob*dp->grada;
  ds->df1001  += 2.0*PREF*factor*dp->rhob*dp->gradb;
  ds->df0101  += 2.0*PREF*factor*dp->rhoa*dp->gradb;
  ds->df0110  += 2.0*PREF*factor*dp->rhoa*dp->grada;
  ds->df10001 += 2.0*PREF*factor*dp->rhob;
  ds->df01001 += 2.0*PREF*factor*dp->rhoa;
}

/* example4_third:
   Example4 functional derivatives.
*/
static void
example4_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
  real r2 = dp->rhoa*dp->rhoa;
  real denom = 0.3+r2;
  real d2    = denom*denom;

  ds->df1000 += PREF*(0.3-r2)/d2*dp->grada*dp->grada*factor;
  ds->df0010 += PREF*dp->rhoa/denom*2*dp->grada*factor;
  ds->df2000 += PREF*2*(r2*dp->rhoa-0.9*dp->rhoa)/(denom*d2)*dp->grada*dp->grada*factor;
  ds->df1010 += PREF*(0.3-r2)/d2*2*dp->grada*factor;
  ds->df0020 += PREF*dp->rhoa/denom*2*factor;

  ds->df3000 += PREF*(-6*r2*r2+10.8*r2-0.54)/(d2*d2)*dp->grada*dp->grada*factor;
  ds->df2010 += PREF*2*(r2*dp->rhoa-0.9*dp->rhoa)/(denom*d2)*2*dp->grada*factor;
  ds->df1020 += PREF*(0.3-r2)/d2*2*factor;
}
