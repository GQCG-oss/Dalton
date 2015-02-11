/*


!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2015 (2015), see http://daltonprogram.org"
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
#include "general.h"

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static integer example6_isgga(void) { return 1; }
static integer example6_read(const char* conf_line);
static real example6_energy(const FunDensProp* dp);
static void example6_first(FunFirstFuncDrv *ds,   real factor, const FunDensProp* dp);
static void example6_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp);
static void example6_third(FunThirdFuncDrv *ds,   real factor, const FunDensProp* dp);

Functional Example6Functional = {
  "Example6",         /* name */
  example6_isgga,     /* gga-corrected */
   1,
  example6_read, 
  NULL,              /* reporter */
  example6_energy, 
  example6_first,
  example6_second,
  example6_third
};

/* IMPLEMENTATION PART */
static integer
example6_read(const char* conf_line)
{
  fun_set_hf_weight(0.0);
  return 1;
}

static const real EPREF= -5e-2;

static real
example6_energy(const FunDensProp* dp)
{
  return EPREF*(dp->grada*dp->gradb);
}
/* example_first:
   derivatives with respect to rho_alpha, and dp->grad_alpha
 */
static void
example6_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
  ds->df0010 +=  EPREF*dp->gradb*factor;
}

static void
example6_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
  ds->df0010 +=  EPREF*dp->gradb*factor;
  ds->df0011 +=  EPREF*factor;
}

/* example_third:
   Test functional derivatives.
*/
static void
example6_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
  ds->df0010 +=  EPREF*dp->gradb*factor;
  ds->df0011 +=  EPREF*factor;
}
