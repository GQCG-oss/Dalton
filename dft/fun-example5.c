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

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int example5_isgga(void) { return 1; }
static int example5_read(const char* conf_line);
static real example5_energy(const DftDensProp* dp);
static void example5_first(FirstFuncDrv *ds,   real factor, const DftDensProp* dp);
static void example5_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp);
static void example5_third(ThirdFuncDrv *ds,   real factor, const DftDensProp* dp);

Functional Example5Functional = {
  "Example5",       /* name */
  example5_isgga,   /* gga-corrected */
  example5_read, 
  NULL,
  example5_energy, 
  example5_first,
  example5_second,
  example5_third
};

/* IMPLEMENTATION PART */
static int
example5_read(const char* conf_line)
{
  dft_set_hf_weight(0.0);
  return 1;
}

static const real PREF= -1.5e-2;

static real
example5_energy(const DftDensProp* dp)
{
  return PREF*(pow(dp->rhoa,1.3)*dp->grada*dp->grada +
               pow(dp->rhob,1.3)*dp->gradb*dp->gradb);
}

static void
example5_first(FirstFuncDrv *ds, real factor, const DftDensProp* dp)
{
  ds->df1000 += PREF*1.3*pow(dp->rhoa,0.3)*dp->grada*dp->grada*factor;
  ds->df0010 += PREF*pow(dp->rhoa,1.3)*2*dp->grada*factor;
}
static void
example5_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp)
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
example5_third(ThirdFuncDrv *ds, real factor, const DftDensProp* dp)
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
