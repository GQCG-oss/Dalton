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
/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* fun-Slater.c:
   implementation of Slater functional and its derivatives 
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

#include "functionals.h"

/* INTERFACE PART */
static int slater_isgga(void) { return 0; }
static int slater_read(const char* conf_line);
static real slater_energy(const FunDensProp* dp);
static void slater_first(FunFirstFuncDrv *ds,   real fac, const FunDensProp*);
static void slater_second(FunSecondFuncDrv *ds, real fac, const FunDensProp*);
static void slater_third(FunThirdFuncDrv *ds,   real fac, const FunDensProp*);
static void slater_fourth(FunFourthFuncDrv *ds, real fac, const FunDensProp*);

Functional SlaterFunctional = {
  "Slater",       /* name */
  slater_isgga,   /* gga-corrected */
  slater_read, 
  NULL,
  slater_energy, 
  slater_first,
  slater_second,
  slater_third,
  slater_fourth
};

/* IMPLEMENTATION PART */
static int
slater_read(const char* conf_line)
{
    fun_set_hf_weight(0);
    return 1;
}

/* SLATER_THRESHOLD Only to avoid numerical problems due to raising 0
 * to a fractional power. */
static const real SLATER_THRESHOLD = 1e-20;
static real
slater_energy(const FunDensProp* dp)
{
  real ea = 0.0, eb = 0.0; 
  const real PREF= -3.0/4.0*pow(6/M_PI, 1.0/3.0);
  if (dp->rhoa >SLATER_THRESHOLD)
      ea= PREF*pow(dp->rhoa,4.0/3.0);
  if (dp->rhob >SLATER_THRESHOLD)
      eb= PREF*pow(dp->rhob,4.0/3.0);   
  return ea+eb;  
}

static void
slater_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
  if (dp->rhoa>SLATER_THRESHOLD)
     ds->df1000 += -pow(6.0/M_PI*dp->rhoa, 1.0/3.0)*factor;
  if (dp->rhob>SLATER_THRESHOLD)
     ds->df0100 += -pow(6.0/M_PI*dp->rhob, 1.0/3.0)*factor;
}
static void
slater_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
  const real PREF = pow(6.0/M_PI, 1.0/3.0);
  if (dp->rhoa>SLATER_THRESHOLD) {
    ds->df1000 += -PREF*pow(dp->rhoa,  1.0/3.0)*factor;
    ds->df2000 += -PREF*pow(dp->rhoa, -2.0/3.0)/3*factor;
  }
  if (dp->rhob>SLATER_THRESHOLD) {
    ds->df0100 += -PREF*pow(dp->rhob,  1.0/3.0)*factor;
    ds->df0200 += -PREF*pow(dp->rhob, -2.0/3.0)/3*factor;
  }
}

/* slater_third:
   Slater functional derivatives.
*/
static void
slater_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
  const real PREF = pow(6.0/M_PI, 1.0/3.0);
  if (dp->rhoa>SLATER_THRESHOLD) {
    ds->df1000 += -PREF*pow(dp->rhoa,  1.0/3.0)*factor;
    ds->df2000 += -PREF*pow(dp->rhoa, -2.0/3.0)/3*factor;
    ds->df3000 +=  PREF*pow(dp->rhoa, -5.0/3.0)*2.0/9.0*factor;
  }
  if (dp->rhob>SLATER_THRESHOLD) {
    ds->df0100 += -PREF*pow(dp->rhob,  1.0/3.0)*factor;
    ds->df0200 += -PREF*pow(dp->rhob, -2.0/3.0)/3*factor;
    ds->df0300 +=  PREF*pow(dp->rhob, -5.0/3.0)*2.0/9.0*factor;
  }
}

/* slater_fourth:
   Dirac functional fourth derivatives.
   by B. Jansik
*/
static void
slater_fourth(FunFourthFuncDrv *ds, real factor, const FunDensProp *dp)
{
    const real PREF = (-3.0/4.0)*pow(6.0/M_PI, 1.0/3.0);/* Dirac G prefactor */
    const real DPREF = 40.0/81.0;	  /* Prefactor from 4th derivative */
    const real JPREF = DPREF*PREF*factor; /* Joined prefactor */
    const real ROEXP = -8.0/3.0;	  /* Exponent on density (from 4th deriv.) */
    FunThirdFuncDrv ds_third;
   
   
   // set up lower order derivatives
   // dirac_third contain third and also lower order derivatives	

   drv3_clear(&ds_third);
   slater_third(&ds_third, factor, dp);
   
   ds->df1000 += ds_third.df1000;
   ds->df2000 += ds_third.df2000;
   ds->df3000 += ds_third.df3000;

   ds->df0100 += ds_third.df0100;
   ds->df0200 += ds_third.df0200;
   ds->df0300 += ds_third.df0300;

   if (dp->rhoa > SLATER_THRESHOLD)
     ds->df4000 += JPREF*pow(dp->rhoa, ROEXP);
   
   if (dp->rhob > SLATER_THRESHOLD)
     ds->df0400 += JPREF*pow(dp->rhob, ROEXP);
}
