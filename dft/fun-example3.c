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
/* fun-test.c:
   implementation of a test GGA-class functional.
   This is a third example functional.
   (c) Pawel Salek, pawsa@theochem.kth.se, aug 2001
   NOTE:
   this file may seem unnecessarily complex but the structure really pays off
   when implementing multiple functionals depending on different parameters.
*/

#include <math.h>
#include <stddef.h>

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int example3_isgga(void) { return 1; }
static int example3_read(const char* conf_line);
static real example3_energy(const FunDensProp* dp);
static void example3_first(FunFirstFuncDrv *ds,   real factor, const FunDensProp* dp);
static void example3_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp);
static void example3_third(FunThirdFuncDrv *ds,   real factor, const FunDensProp* dp);

Functional Example3Functional = {
  "Example3",         /* name */
  example3_isgga,     /* gga-corrected */
  example3_read, 
  NULL,              /* reporter */
  example3_energy, 
  example3_first,
  example3_second,
  example3_third
};

/* IMPLEMENTATION PART */
static int
example3_read(const char* conf_line)
{
  fun_set_hf_weight(0.0);
  return 1;
}

#define LEVEL 3
#if LEVEL==2
static const real EPREF= -1e-2;
#else
static const real EPREF= -5e-2;
#endif

static real
example3_energy(const FunDensProp* dp)
{
#if LEVEL==2
  return EPREF*(pow(dp->grada,4.0)+pow(dp->gradb,4.0)); 
#else
  return EPREF*(pow(dp->grada,1.7)+pow(dp->gradb,1.7));
#endif
}
/* example_first:
   derivatives with respect to rho_alpha, and zeta=|\nabla\grad_alpha|^2
 */
static void
example3_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
#if LEVEL==2
  ds->df0010 += EPREF*3*pow(dp->grada,3.0)*factor;
#else
  ds->df0010 += EPREF*1.7*pow(dp->grada,0.7)*factor;
#endif
}

/* example3_second:
   derivatives with respect to rho_alpha, and zeta=|\nabla\grad_alpha|
   (OBS: inconsistency!)
*/
static void
example3_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
#if LEVEL==2
  ds->df0010 += EPREF* 4*pow(dp->grada,3.0)*factor;
  ds->df0020 += EPREF*12*pow(dp->grada,2.0)*factor;
#else
  ds->df0010 +=  EPREF*1.7*pow(dp->grada,0.7)*factor;
  ds->df0020 +=  EPREF*1.7*0.7*pow(dp->grada,-0.3)*factor;
#endif
}

/* example_third:
   Test functional derivatives.
   derivatives with respect to rho_alpha, and zeta=|\nabla\grad_alpha|^2
*/
static void
example3_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
#if LEVEL==2
  ds->df0010 += EPREF* 4*pow(dp->grada,3.0)*factor;
  ds->df0020 += EPREF*12*pow(dp->grada,2.0)*factor;
#else
  ds->df0010 +=  EPREF*1.7*pow(dp->grada,0.7)*factor;
  ds->df0020 +=  EPREF*1.7*0.7*pow(dp->grada,-0.3)*factor;
  ds->df0030 += -EPREF*1.7*0.7*0.3*pow(dp->grada,-1.3)*factor;
#endif
}
