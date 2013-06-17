/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* fun-Dirac.c:
   implementation of Dirac functional and its derivatives 
   (c), Pawel Salek, pawsa@theochem.kth.se, aug 2001
   Z. Rinkevicius adapted for open shell systems: energy, first derivatives.
   NOTE:
   this file may seem unnecessarily complex but the structure really pays off
   when implementing multiple functionals depending on different parameters.
*/

/* strictly conform to XOPEN ANSI C standard */
#if !defined(SYS_DEC)
/* XOPEN compliance is missing on old Tru64 4.0E Alphas and pow() prototype
 * is not specified. */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1
#endif
#include <math.h>
#include <stdio.h>

#define __CVERSION__

#include "lsdalton_functionals.h"

/* INTERFACE PART */
static integer dirac_isgga(void) { return 0; }
static integer dirac_read(const char* conf_line, real *hfweight);
static real dirac_energy(const DftDensProp* dp);
static void dirac_first(FirstFuncDrv *ds,   real factor, const DftDensProp* dp);
static void dirac_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp);
static void dirac_third(ThirdFuncDrv *ds,   real factor, const DftDensProp* dp);

Functional DiracFunctional = {
  "Dirac",       /* name */
  dirac_isgga,   /* gga-corrected */
  dirac_read, 
  NULL,
  dirac_energy, 
  dirac_first,
  dirac_second,
  dirac_third
};

/* IMPLEMENTATION PART */
static integer
dirac_read(const char* conf_line, real *hfweight)
{
    return 1;
}

/* DIRAC_THRESHOLD Only to avoid numerical problems due to raising 0
 * to a fractional power. */
static const real DIRAC_THRESHOLD = 1e-14;
static real
dirac_energy(const DftDensProp* dp)
{
  real ea = 0.0, eb = 0.0; 
#if TEST_FUNCTIONAL==1
  return -0.05*(rhoa + rhob);
#elif TEST_FUNCTIONAL==2
  return -0.05*(rhoa*rhoa + rhob*rhob);
#else 
  const real PREF= -3.0/4.0*pow(6/M_PI, 1.0/3.0);
  if (dp->rhoa >DIRAC_THRESHOLD)
      ea= PREF*pow(dp->rhoa,4.0/3.0);
  if (dp->rhob >DIRAC_THRESHOLD)
      eb= PREF*pow(dp->rhob,4.0/3.0);   
  return ea+eb;  
#endif
}

static void
dirac_first(FirstFuncDrv *ds, real factor, const DftDensProp* dp)
{
#if TEST_FUNCTIONAL==1
  ds->df1000 += -0.05*factor;
#elif TEST_FUNCTIONAL==2
  ds->df1000 += -0.1*rhoa*factor;
#else
  if (dp->rhoa>DIRAC_THRESHOLD)
     ds->df1000 += -pow(6.0/M_PI*dp->rhoa, 1.0/3.0)*factor;
  if (dp->rhob>DIRAC_THRESHOLD)
     ds->df0100 += -pow(6.0/M_PI*dp->rhob, 1.0/3.0)*factor;
#endif
}
static void
dirac_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp)
{
#if TEST_FUNCTIONAL==1
  ds->df1000 += -0.05*factor;
#elif TEST_FUNCTIONAL==2
  ds->df1000 += -0.1*rhoa*factor;
  ds->df2000 += -0.1*factor;
#else 
  const real PREF = pow(6.0/M_PI, 1.0/3.0);
  if (dp->rhoa>DIRAC_THRESHOLD) {
    ds->df1000 += -PREF*pow(dp->rhoa,  1.0/3.0)*factor;
    ds->df2000 += -PREF*pow(dp->rhoa, -2.0/3.0)/3*factor;
  }
  if (dp->rhob>DIRAC_THRESHOLD) {
    ds->df0100 += -PREF*pow(dp->rhob,  1.0/3.0)*factor;
    ds->df0200 += -PREF*pow(dp->rhob, -2.0/3.0)/3*factor;
  }
#endif
}

/* dirac_third:
   Dirac functional derivatives.
*/
static void
dirac_third(ThirdFuncDrv *ds, real factor, const DftDensProp* dp)
{
#if TEST_FUNCTIONAL==1
  ds->df1000 += -0.05*factor;
#elif TEST_FUNCTIONAL==2
  ds->df1000 += -0.1*rhoa*factor;
  ds->df2000 += -0.1*factor;
#else
  const real PREF = pow(6.0/M_PI, 1.0/3.0);
  if (dp->rhoa>DIRAC_THRESHOLD) {
    ds->df1000 += -PREF*pow(dp->rhoa,  1.0/3.0)*factor;
    ds->df2000 += -PREF*pow(dp->rhoa, -2.0/3.0)/3*factor;
    ds->df3000 +=  PREF*pow(dp->rhoa, -5.0/3.0)*2.0/9.0*factor;
  }
  if (dp->rhob>DIRAC_THRESHOLD) {
    ds->df0100 += -PREF*pow(dp->rhob,  1.0/3.0)*factor;
    ds->df0200 += -PREF*pow(dp->rhob, -2.0/3.0)/3*factor;
    ds->df0300 +=  PREF*pow(dp->rhob, -5.0/3.0)*2.0/9.0*factor;
  }
#endif
}
