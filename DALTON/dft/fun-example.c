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
/* fun-example.c:
   implementation of a test GGA-class functional
   (c) Pawel Salek, pawsa@theochem.kth.se, aug 2001
   NOTE:
   this file may seem unnecessarily complex but the structure really pays off
   when implementing multiple functionals depending on different parameters.
*/

#include <math.h>
#include <stdio.h>
#include "general.h"

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static integer example_isgga(void) { return 1; }
static integer example_read(const char* conf_line);
static real example_energy(const FunDensProp* dp);
static void example_first(FunFirstFuncDrv *ds,   real factor, const FunDensProp* dp);
static void example_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp);
static void example_third(FunThirdFuncDrv *ds,   real factor, const FunDensProp* dp);

Functional ExampleFunctional = {
  "Example",         /* name */
  example_isgga,     /* gga-corrected */
   1,
  example_read, 
  NULL,              /* reporter */
  example_energy, 
  example_first,
  example_second,
  example_third
};

/* IMPLEMENTATION PART */
static integer
example_read(const char* conf_line)
{
  fun_set_hf_weight(0.0);
  return 1;
}

#if LEVEL==2
static const real EPREF= -1e0;
#else
static const real EPREF= -2e-3;
#endif

static real
example_energy(const FunDensProp* dp)
{
#if LEVEL==2
  return EPREF*(dp->rhoa*dp->rhoa*dp->rhob);
#else
  return EPREF*(pow(dp->rhoa,4.0/3.0)*dp->grada*dp->grada +
                pow(dp->rhob,4.0/3.0)*dp->gradb*dp->gradb);
#endif
}
/* example_first:
   derivatives with respect to rho_alpha, and zeta=|grad_alpha|^2
 */
static void
example_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
#if LEVEL==2
  ds->df1000 +=  EPREF*2*dp->rhoa*dp->rhob*factor;
#else
  ds->df1000 += EPREF*4.0/3.0*pow(dp->rhoa,1.0/3.0)*dp->grada*dp->grada*factor;
  ds->df0010 += EPREF*pow(dp->rhoa,4.0/3.0)*2*dp->grada*factor;
#endif
}

static void
example_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
#if LEVEL==2
  ds->df1000 +=  EPREF*2*dp->rhoa*dp->rhob*factor;
  ds->df2000 +=  EPREF*2*dp->rhob*factor;
  ds->df1100 +=  EPREF*2*dp->rhoa*factor;
#else
  const real q = dp->grada*dp->grada;
  ds->df1000 +=  EPREF*4.0/3.0*pow(dp->rhoa,1.0/3.0)*q*factor;
  ds->df0010 +=  EPREF*pow(dp->rhoa,4.0/3.0)*2*dp->grada*factor;
  ds->df2000 +=  EPREF*4.0/9.0 *pow(dp->rhoa,-2.0/3.0)*q*factor;
  ds->df1010 +=  EPREF*4.0/3.0 *pow(dp->rhoa,1.0/3.0)*2*dp->grada*factor;
  ds->df0020 +=  EPREF*pow(dp->rhoa,4.0/3.0)*2*factor;
#endif
}

/* example_third:
   Test functional derivatives.
*/
static void
example_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
#if LEVEL==2
  ds->df1000 +=  EPREF*2*dp->rhoa*dp->rhob*factor;
  ds->df2000 +=  EPREF*2*dp->rhob*factor;
  ds->df1100 +=  EPREF*2*dp->rhoa*factor;
  ds->df2100 +=  EPREF*2*factor;
#else
  const real q = dp->grada*dp->grada;
  ds->df1000 +=  EPREF*4.0/3.0*pow(dp->rhoa,1.0/3.0)*q*factor;
  ds->df0010 +=  EPREF*pow(dp->rhoa,4.0/3.0)*2*dp->grada*factor;
  ds->df2000 +=  EPREF*4.0/9.0 *pow(dp->rhoa,-2.0/3.0)*q*factor;
  ds->df1010 +=  EPREF*4.0/3.0 *pow(dp->rhoa, 1.0/3.0)*2*dp->grada*factor;
  ds->df0020 +=  EPREF*pow(dp->rhoa,4.0/3.0)*2*factor;

  ds->df3000 += -EPREF*8.0/27.0*pow(dp->rhoa,-5.0/3.0)*q*factor;
  ds->df2010 +=  EPREF*4.0/9.0 *pow(dp->rhoa,-2.0/3.0)*2*dp->grada*factor;
  ds->df1020 +=  EPREF*4.0/3.0 *pow(dp->rhoa, 1.0/3.0)*2*factor;
#endif
}
