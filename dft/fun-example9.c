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
static int example9_isgga(void) { return 1; }
static int example9_read(const char* conf_line);
static real example9_energy(const DftDensProp* dp);
static void example9_first(FirstFuncDrv *ds,   real factor, const DftDensProp* dp);
static void example9_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp);
static void example9_third(ThirdFuncDrv *ds,   real factor, const DftDensProp* dp);

Functional Example9Functional = {
  "Example9",         /* name */
  example9_isgga,     /* gga-corrected */
  example9_read, 
  NULL,              /* reporter */
  example9_energy, 
  example9_first,
  example9_second,
  example9_third
};

/* IMPLEMENTATION PART */
static int
example9_read(const char* conf_line)
{
  fun_set_hf_weight(0.0);
  return 1;
}

static const real EPREF= -5e-2;

static real
example9_energy(const DftDensProp* dp)
{
  return EPREF*(dp->gradab);
}
/* example_first:
   derivatives with respect to rho_alpha, and grad_alpha
 */
static void
example9_first(FirstFuncDrv *ds, real factor, const DftDensProp* dp)
{
  ds->df00001 +=  EPREF*factor;
}

static void
example9_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp)
{
  ds->df00001 +=  EPREF*factor;
}

/* example_third:
   Test functional derivatives.
*/
static void
example9_third(ThirdFuncDrv *ds, real factor, const DftDensProp* dp)
{
  ds->df00001 +=  EPREF*factor;
}
