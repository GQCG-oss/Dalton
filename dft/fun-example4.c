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
/* fun-example4.c:
   implementation of Example4 functional and its derivatives 
   (c) Pawel Salek, pawsa@theochem.kth.se, aug 2001
   NOTE:
   this file may seem unnecessarily complex but the structure really pays off
   when implementing multiple functionals depending on different parameters.
*/
/* strictly conform to XOPEN ANSI C standard */
#define __USE_XOPEN 

#include <math.h>
#include <stdio.h>

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int example4_isgga(void) { return 1; }
static int example4_read(const char* conf_line);
static real example4_energy(const FunDensProp* dp);
static void example4_first(FunFirstFuncDrv *ds,   real factor, const FunDensProp* dp);
static void example4_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp);
static void example4_third(FunThirdFuncDrv *ds,   real factor, const FunDensProp* dp);

Functional Example4Functional = {
  "Example4",       /* name */
  example4_isgga,   /* gga-corrected */
  example4_read, 
  NULL,
  example4_energy, 
  example4_first,
  example4_second,
  example4_third
};

/* IMPLEMENTATION PART */
static int
example4_read(const char* conf_line)
{
    fun_set_hf_weight(0);
    return 1;
}

static const real PREF= -5e-5;

static real
example4_energy(const FunDensProp* dp)
{
  real grad2 = dp->grada*dp->grada+dp->gradb*dp->gradb+2.0*dp->gradab;
  return PREF*dp->rhoa*dp->rhob*grad2;
}

static void
example4_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
  real grada2 = dp->grada*dp->grada;
  real gradb2 = dp->gradb*dp->gradb;    
  real grad2 = grada2+gradb2+2.0*dp->gradab; 
  real rhoab = dp->rhoa*dp->rhob;
  ds->df1000  += PREF*factor*dp->rhob*grad2;
  ds->df0100  += PREF*factor*dp->rhoa*grad2;
  ds->df0010  += 2.0*PREF*factor*rhoab*dp->grada;
  ds->df0001  += 2.0*PREF*factor*rhoab*dp->gradb;
  ds->df00001 += 2.0*PREF*factor*rhoab; 
}
static void
example4_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
  real grada2 = dp->grada*dp->grada;
  real gradb2 = dp->gradb*dp->gradb;    
  real grad2 = grada2+gradb2+2.0*dp->gradab; 
  real rhoab = dp->rhoa*dp->rhob;
  /* first derivatives */
  ds->df1000  += PREF*factor*dp->rhob*grad2;
  ds->df0100  += PREF*factor*dp->rhoa*grad2;
  ds->df0010  += 2.0*PREF*factor*rhoab*dp->grada;
  ds->df0001  += 2.0*PREF*factor*rhoab*dp->gradb;
  ds->df00001 += 2.0*PREF*factor*rhoab; 
  /* second derivatives */
  ds->df0020  += 2.0*PREF*factor*rhoab;
  ds->df0002  += 2.0*PREF*factor*rhoab; 
  /* mixed derivatives */
  ds->df1100  += PREF*factor*grad2;
  ds->df1010  += 2.0*PREF*factor*dp->rhob*dp->grada;
  ds->df1001  += 2.0*PREF*factor*dp->rhob*dp->gradb;
  ds->df0101  += 2.0*PREF*factor*dp->rhoa*dp->gradb;
  ds->df0110  += 2.0*PREF*factor*dp->rhoa*dp->grada;
  ds->df10001 += 2.0*PREF*factor*dp->rhob;
  ds->df01001 += 2.0*PREF*factor*dp->rhoa;
}

/* example4_third:
   Example4 functional derivatives.
*/
static void
example4_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
  real r2 = dp->rhoa*dp->rhoa;
  real denom = 0.3+r2;
  real d2    = denom*denom;

  ds->df1000 += PREF*(0.3-r2)/d2*dp->grada*dp->grada*factor;
  ds->df0010 += PREF*dp->rhoa/denom*2*dp->grada*factor;
  ds->df2000 += PREF*2*(r2*dp->rhoa-0.9*dp->rhoa)/(denom*d2)*dp->grada*dp->grada*factor;
  ds->df1010 += PREF*(0.3-r2)/d2*2*dp->grada*factor;
  ds->df0020 += PREF*dp->rhoa/denom*2*factor;

  ds->df3000 += PREF*(-6*r2*r2+10.8*r2-0.54)/(d2*d2)*dp->grada*dp->grada*factor;
  ds->df2010 += PREF*2*(r2*dp->rhoa-0.9*dp->rhoa)/(denom*d2)*2*dp->grada*factor;
  ds->df1020 += PREF*(0.3-r2)/d2*2*factor;
}
