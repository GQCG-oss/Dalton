/* fun-pw92.c:
   implementation of the electron-gas correlation energy
   by J. P. Perdew and Y. Wang;
   J.P.Perdew,Y. Wang; Phys. Rew. B; 40, 13244, (1992) [Ref 1]
   
   The PW92 functional and its derivatives up to fourth order.
   Implementation (c)by B. Jansik, brano@theochem.kth.se, jan 2003

   NOTE: Improvement over VWN:
   While we confirm the practical accuracy of the VWN and PZ reprezentations,
   we eliminate some minor problems with these forms. [1]
   See ref [1] for details.

   NOTE: This functional is the standalone PW92 functional
         constructed as rho*PW92pe where rho is density and
         PW92pe is the 'per electron' pw92 functional as given in
         ref [1].

*/

#if !defined(SYS_DEC)
/* XOPEN compliance is missing on old Tru64 4.0E Alphas and pow() prototype
 * is not specified. */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1
#endif

#include <math.h>
#include <stddef.h>

#define __CVERSION__

#include "functionals.h"
/* #define FOURTH_ORDER_DERIVATIVES */


/* INTERFACE PART */
static integer  pw92_isgga(void) { return 0; }
static integer  pw92_read(const char* conf_line, real *hfweight);
static real pw92_energy(const DftDensProp* dens_prop);
static void pw92_first(FirstFuncDrv *ds, real factor,
                       const DftDensProp* dens_prop);
static void pw92_second(SecondFuncDrv *ds, real factor,
                        const DftDensProp* dens_prop);

static void pw92_third(ThirdFuncDrv *ds, real factor,
                       const DftDensProp* dens_prop);

#ifdef FOURTH_ORDER_DERIVATIVES
static void pw92_fourth(FourthFuncDrv *ds, real factor,
                        const DftDensProp* dens_prop);
#endif
Functional PW92Functional = {
    "PW92",      /* name */
    pw92_isgga,  /* gga-corrected */
    pw92_read,   /* set bloody common blocks */
    NULL,         /* reporter */
    pw92_energy, 
    pw92_first,
                             pw92_second,
    pw92_third
#ifdef FOURTH_ORDER_DERIVATIVES
    ,pw92_fourth
#endif
};


/* IMPLEMENTATION PART */

static integer
pw92_read(const char* conf_line, real *hfweight)
{
/*  dft_set_hf_weight(0.0); */
    return 1;
}

static real
pw92_energy(const DftDensProp* dp)
{
  real E = (dp->rhoa+dp->rhob)*PW92peFunctional.func(dp);
  return (E);
}

static void
pw92_first(FirstFuncDrv *ds, real factor, const DftDensProp* dp)
{
  /* Structure of "PW92 per electron" derivatives*/
  FirstFuncDrv drv;
  drv1_clear(&drv);
  PW92peFunctional.first(&drv, 1.0, dp);

  /* Energy */
  real roa_rob = dp->rhoa + dp->rhob;
  real E = PW92peFunctional.func(dp);

  /* First order derivatives */
  ds->df1000 += factor*(E + roa_rob*drv.df1000);
  ds->df0100 += factor*(E + roa_rob*drv.df0100);
 
}

static void
pw92_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp)
{

  /* Structure of "PW92 per electron" derivatives*/
  SecondFuncDrv drv;
  drv2_clear(&drv);
  PW92peFunctional.second(&drv, 1.0, dp);

  /* Energy */
  real roa_rob = dp->rhoa + dp->rhob;
  real E = PW92peFunctional.func(dp);

  /* First order derivatives */
  ds->df1000 += factor*(E + roa_rob*drv.df1000);
  ds->df0100 += factor*(E + roa_rob*drv.df0100);

  /* Second order derivatives */
  ds->df2000 += factor*(2.0*drv.df1000 + roa_rob*drv.df2000);
  ds->df1100 += factor*(drv.df1000 + drv.df0100 + roa_rob*drv.df1100);
  ds->df0200 += factor*(2.0*drv.df0100 + roa_rob*drv.df0200);

  /* End of derivatives */
}


static void
pw92_third(ThirdFuncDrv *ds, real factor, const DftDensProp* dp)
{

  /* Structure of "PW92 per electron" derivatives*/
  ThirdFuncDrv drv;
  drv3_clear(&drv);
  PW92peFunctional.third(&drv, 1.0, dp);

  /* Energy */
  real roa_rob = dp->rhoa + dp->rhob;
  real E = PW92peFunctional.func(dp);

  /* First order derivatives */
  ds->df1000 += factor*(E + roa_rob*drv.df1000);
  ds->df0100 += factor*(E + roa_rob*drv.df0100);

  /* Second order derivatives */
  ds->df2000 += factor*(2.0*drv.df1000 + roa_rob*drv.df2000);
  ds->df1100 += factor*(drv.df1000 + drv.df0100 + roa_rob*drv.df1100);
  ds->df0200 += factor*(2.0*drv.df0100 + roa_rob*drv.df0200);
  
  /* Third order derivatives */
  ds->df3000 += factor*(3.0*drv.df2000 + roa_rob*drv.df3000);
  ds->df2100 += factor*(2.0*drv.df1100 + drv.df2000 + roa_rob*drv.df2100);
  ds->df1200 += factor*(2.0*drv.df1100 + drv.df0200 + roa_rob*drv.df1200);
  ds->df0300 += factor*(3.0*drv.df0200 + roa_rob*drv.df0300);

  /* End of derivatives */
}

#ifdef FOURTH_ORDER_DERIVATIVES
static void
pw92_fourth(FourthFuncDrv *ds, real factor, const DftDensProp* dp)
{

  /* Structure of "PW92 per electron" derivatives*/
  FourthFuncDrv drv;
  drv4_clear(&drv);
  PW92peFunctional.fourth(&drv, 1.0, dp);

  /* Energy */
  real roa_rob = dp->rhoa + dp->rhob;
  real E = PW92peFunctional.func(dp);

  /* First order derivatives */
  ds->df10000 += factor*(E + roa_rob*drv.df10000);
  ds->df01000 += factor*(E + roa_rob*drv.df01000);

  /* Second order derivatives */
  ds->df20000 += factor*(2.0*drv.df10000 + roa_rob*drv.df20000);
  ds->df11000 += factor*(drv.df10000 + drv.df01000 + roa_rob*drv.df11000);
  ds->df02000 += factor*(2.0*drv.df01000 + roa_rob*drv.df02000);
  
  /* Third order derivatives */
  ds->df30000 += factor*(3.0*drv.df20000 + roa_rob*drv.df30000);
  ds->df21000 += factor*(2.0*drv.df11000 + drv.df20000 + roa_rob*drv.df21000);
  ds->df12000 += factor*(2.0*drv.df11000 + drv.df02000 + roa_rob*drv.df12000);
  ds->df03000 += factor*(3.0*drv.df02000 + roa_rob*drv.df03000);

  /* Fourth order derivatives */
  ds->df40000 += factor*(4.0*drv.df30000 + roa_rob*drv.df40000);
  ds->df31000 += factor*(3.0*drv.df21000 + drv.df30000 + roa_rob*drv.df31000);
  ds->df22000 += factor*(2.0*drv.df12000 + 2.0*drv.df21000 + roa_rob*drv.df22000);
  ds->df13000 += factor*(3.0*drv.df12000 + drv.df03000 + roa_rob*drv.df13000);
  ds->df04000 += factor*(4.0*drv.df03000 + roa_rob*drv.df04000);

  /* End of derivatives */
}
#endif
