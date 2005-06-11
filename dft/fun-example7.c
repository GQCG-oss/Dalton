/*
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras, T. Saue, 
C...   S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren. 
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
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
