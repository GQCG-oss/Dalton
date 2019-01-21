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
/* fun-example5.c:
   implementation of Example5 functional and its derivatives 
   (c) Pawel Salek, pawsa@theochem.kth.se, aug 2001
   NOTE:
   this file may seem unnecessarily complex but the structure really pays off
   when implementing multiple functionals depending on different parameters.
*/
/* strictly conform to XOPEN ANSI C standard */
#define __USE_XOPEN 

#include <math.h>
#include <stddef.h>
#include "general.h"

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static integer example5_isgga(void) { return 1; }
static integer example5_read(const char* conf_line);
static real example5_energy(const FunDensProp* dp);
static void example5_first(FunFirstFuncDrv *ds,   real factor, const FunDensProp* dp);
static void example5_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp);
static void example5_third(FunThirdFuncDrv *ds,   real factor, const FunDensProp* dp);

Functional Example5Functional = {
  "Example5",       /* name */
  example5_isgga,   /* gga-corrected */
   1,
  example5_read, 
  NULL,
  example5_energy, 
  example5_first,
  example5_second,
  example5_third
};

/* IMPLEMENTATION PART */
static integer
example5_read(const char* conf_line)
{
  fun_set_hf_weight(0.0);
  return 1;
}

static const real PREF= -1.5e-2;

static real
example5_energy(const FunDensProp* dp)
{
  return PREF*(pow(dp->rhoa,1.3)*dp->grada*dp->grada +
               pow(dp->rhob,1.3)*dp->gradb*dp->gradb);
}

static void
example5_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
  ds->df1000 += PREF*1.3*pow(dp->rhoa,0.3)*dp->grada*dp->grada*factor;
  ds->df0010 += PREF*pow(dp->rhoa,1.3)*2*dp->grada*factor;
}
static void
example5_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
  ds->df1000 += PREF*1.3*pow(dp->rhoa,0.3)*dp->grada*dp->grada*factor;
  ds->df0010 += PREF*pow(dp->rhoa,1.3)*2*dp->grada*factor;
  ds->df2000 += PREF*1.3*0.3*pow(dp->rhoa,-0.7)*dp->grada*dp->grada*factor;
  ds->df1010 += PREF*1.3*pow(dp->rhoa,0.3)*2*dp->grada*factor;
  ds->df0020 += PREF*pow(dp->rhoa,1.3)*2*factor;
}

/* example5_third:
   Example5 functional derivatives.
*/
static void
example5_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
  ds->df1000 += PREF*1.3*pow(dp->rhoa,0.3)*dp->grada*dp->grada*factor;
  ds->df0010 += PREF*pow(dp->rhoa,1.3)*2*dp->grada*factor;
  ds->df2000 += PREF*1.3*0.3*pow(dp->rhoa,-0.7)*dp->grada*dp->grada*factor;
  ds->df1010 += PREF*1.3*pow(dp->rhoa,0.3)*2*dp->grada*factor;
  ds->df0020 += PREF*pow(dp->rhoa,1.3)*2*factor;
  ds->df3000 +=-PREF*1.3*0.3*0.7*pow(dp->rhoa,-1.7)*dp->grada*dp->grada*factor;
  ds->df2010 += PREF*1.3*0.3*pow(dp->rhoa,-0.7)*2*dp->grada*factor;
  ds->df1020 += PREF*1.3*pow(dp->rhoa,0.3)*2*factor;
}
