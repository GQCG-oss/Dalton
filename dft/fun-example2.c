/* fun-test.c:
   implementation of a test GGA-class functional.
   This is a second example functional.
   (c) Pawel Salek, pawsa@theochem.kth.se, aug 2001
   NOTE:
   this file may seem unnecessarily complex but the structure really pays off
   when implementing multiple functionals depending on different parameters.
*/

#include <stddef.h>

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int example2_isgga(void) { return 1; }
static int example2_read(const char* conf_line);
static real example2_energy(const DftDensProp* dp);
static void example2_first(FirstFuncDrv *ds,   real factor, const DftDensProp* dp);
static void example2_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp);
static void example2_third(ThirdFuncDrv *ds,   real factor, const DftDensProp* dp);

Functional Example2Functional = {
  "Example2",         /* name */
  example2_isgga,     /* gga-corrected */
  example2_read, 
  NULL,              /* reporter */
  example2_energy, 
  example2_first,
  example2_second,
  example2_third
};

/* IMPLEMENTATION PART */
static int
example2_read(const char* conf_line)
{
  dft_set_hf_weight(0.0);
  return 1;
}

static const real EPREF= -5e-5;

static real
example2_energy(const DftDensProp* dp)
{
  return EPREF*(dp->rhoa*dp->grada*dp->grada+dp->rhob*dp->gradb*dp->gradb);
}
/* example_first:
   derivatives with respect to dp->rho_alpha, and dp->grad_alpha
 */
static void
example2_first(FirstFuncDrv *ds, real factor, const DftDensProp* dp)
{
  ds->df1000 +=  EPREF*dp->grada*dp->grada*factor;
  ds->df0010 +=  EPREF*dp->rhoa*2*dp->grada*factor;
}

static void
example2_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp)
{
  ds->df1000 +=  EPREF*dp->grada*dp->grada*factor;
  ds->df0010 +=  EPREF*dp->rhoa*2*dp->grada*factor;
  ds->df1010 +=  EPREF*2*dp->grada*factor;
  ds->df0020 +=  EPREF*dp->rhoa*2*factor;
}

/* example_third:
   Test functional derivatives.
*/
static void
example2_third(ThirdFuncDrv *ds, real factor, const DftDensProp* dp)
{
  ds->df1000 +=  EPREF*dp->grada*dp->grada*factor;
  ds->df0010 +=  EPREF*dp->rhoa*2*dp->grada*factor;
  ds->df1010 +=  EPREF*2*dp->grada*factor;
  ds->df0020 +=  EPREF*dp->rhoa*2*factor;

  ds->df1020 +=  EPREF*2*factor;
}
