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
static int example7_isgga(void) { return 1; }
static int example7_read(const char* conf_line);
static real example7_energy(const FunDensProp* dp);
static void example7_first(FunFirstFuncDrv *ds,   real factor,
                           const FunDensProp* dp);
static void example7_second(FunSecondFuncDrv *ds, real factor,
                            const FunDensProp* dp);
static void example7_third(FunThirdFuncDrv *ds,   real factor,
                           const FunDensProp* dp);

Functional Example7Functional = {
  "Example7",         /* name */
  example7_isgga,     /* gga-corrected */
  example7_read, 
  NULL,              /* reporter */
  example7_energy, 
  example7_first,
  example7_second,
  example7_third
};

/* IMPLEMENTATION PART */
static int
example7_read(const char* conf_line)
{
  fun_set_hf_weight(0.0);
  return 1;
}

static const real EPREF= -5e-2;

static real
example7_energy(const FunDensProp* dp)
{
  return EPREF*(dp->rhoa*dp->gradb+dp->rhob*dp->grada);
}
/* example_first:
   derivatives with respect to dp->rho_alpha, and dp->grad_alpha
 */
static void
example7_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
  ds->df1000 +=  EPREF*dp->gradb*factor;
  ds->df0010 +=  EPREF*dp->rhob*factor;
}

static void
example7_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
  ds->df1000 +=  EPREF*dp->gradb*factor;
  ds->df0010 +=  EPREF*dp->rhob*factor;
  ds->df1001 +=  EPREF*factor;
}

/* example_third:
   Test functional derivatives.
*/
static void
example7_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
  ds->df1000 +=  EPREF*dp->gradb*factor;
  ds->df0010 +=  EPREF*dp->rhob*factor;
  ds->df1001 +=  EPREF*factor;
}
