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
/* fun-example.c:
   implementation of a test GGA-class functional
   (c) Pawel Salek, pawsa@theochem.kth.se, aug 2001
   NOTE:
   this file may seem unnecessarily complex but the structure really pays off
   when implementing multiple functionals depending on different parameters.
*/

#include <math.h>
#include <stdio.h>

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int example_isgga(void) { return 1; }
static int example_read(const char* conf_line);
static real example_energy(const FunDensProp* dp);
static void example_first(FunFirstFuncDrv *ds,   real factor, const FunDensProp* dp);
static void example_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp);
static void example_third(FunThirdFuncDrv *ds,   real factor, const FunDensProp* dp);

Functional ExampleFunctional = {
  "Example",         /* name */
  example_isgga,     /* gga-corrected */
  example_read, 
  NULL,              /* reporter */
  example_energy, 
  example_first,
  example_second,
  example_third
};

/* IMPLEMENTATION PART */
static int
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
