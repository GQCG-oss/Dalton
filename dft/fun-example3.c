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
static real example3_energy(const DftDensProp* dp);
static void example3_first(FirstFuncDrv *ds,   real factor, const DftDensProp* dp);
static void example3_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp);
static void example3_third(ThirdFuncDrv *ds,   real factor, const DftDensProp* dp);

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
example3_energy(const DftDensProp* dp)
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
example3_first(FirstFuncDrv *ds, real factor, const DftDensProp* dp)
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
example3_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp)
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
example3_third(ThirdFuncDrv *ds, real factor, const DftDensProp* dp)
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
