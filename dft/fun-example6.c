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
static int example6_isgga(void) { return 1; }
static int example6_read(const char* conf_line);
static real example6_energy(const DftDensProp* dp);
static void example6_first(FirstFuncDrv *ds,   real factor, const DftDensProp* dp);
static void example6_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp);
static void example6_third(ThirdFuncDrv *ds,   real factor, const DftDensProp* dp);

Functional Example6Functional = {
  "Example6",         /* name */
  example6_isgga,     /* gga-corrected */
  example6_read, 
  NULL,              /* reporter */
  example6_energy, 
  example6_first,
  example6_second,
  example6_third
};

/* IMPLEMENTATION PART */
static int
example6_read(const char* conf_line)
{
  fun_set_hf_weight(0.0);
  return 1;
}

static const real EPREF= -5e-2;

static real
example6_energy(const DftDensProp* dp)
{
  return EPREF*(dp->grada*dp->gradb);
}
/* example_first:
   derivatives with respect to rho_alpha, and dp->grad_alpha
 */
static void
example6_first(FirstFuncDrv *ds, real factor, const DftDensProp* dp)
{
  ds->df0010 +=  EPREF*dp->gradb*factor;
}

static void
example6_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp)
{
  ds->df0010 +=  EPREF*dp->gradb*factor;
  ds->df0011 +=  EPREF*factor;
}

/* example_third:
   Test functional derivatives.
*/
static void
example6_third(ThirdFuncDrv *ds, real factor, const DftDensProp* dp)
{
  ds->df0010 +=  EPREF*dp->gradb*factor;
  ds->df0011 +=  EPREF*factor;
}
