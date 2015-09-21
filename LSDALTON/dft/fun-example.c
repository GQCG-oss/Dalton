/*
C...   Copyright (c) 2015 by the authors of Dalton (see below).
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
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
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
/*
C...   Copyright (c) 2015 by the authors of Dalton (see below).
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
/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* fun-example.c:
   implementation of a test GGA-class functional
   (c) Pawel Salek, pawsa@theochem.kth.se, aug 2001
   NOTE:
   this file may seem unnecessarily complex but the structure really pays off
   when implementing multiple functionals depending on different parameters.
*/

#include <math.h>
#include <stdio.h>
#include "lsdalton_general.h"

#define __CVERSION__

#include "lsdalton_functionals.h"

/* INTERFACE PART */
static integer example_isgga(void) { return 1; }
static integer example_read(const char* conf_line, real *hfweight);
static real example_energy(const FunDensProp* dp);
static void example_first(FunFirstFuncDrv *ds,   real factor, const FunDensProp* dp);
static void example_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp);
static void example_third(FunThirdFuncDrv *ds,   real factor, const FunDensProp* dp);

Functional ExampleFunctional = {
  "Example",         /* name */
  example_isgga,     /* gga-corrected */
  example_read, 
  NULL,              /* reporter */
  example_energy, 
  example_first,
  example_second,
  example_third
};

/* IMPLEMENTATION PART */
static integer
example_read(const char* conf_line, real *hfweight)
{
  return 1;
}

#if LEVEL==2
static const real EPREF= -1e0;
#else
static const real EPREF= -2e-3;
#endif

static real
example_energy(const FunDensProp* dp)
{
#if LEVEL==2
  return EPREF*(dp->rhoa*dp->rhoa*dp->rhob);
#else
  return EPREF*(pow(dp->rhoa,4.0/3.0)*dp->grada*dp->grada +
                pow(dp->rhob,4.0/3.0)*dp->gradb*dp->gradb);
#endif
}
/* example_first:
   derivatives with respect to rho_alpha, and zeta=|grad_alpha|^2
 */
static void
example_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
#if LEVEL==2
  ds->df1000 +=  EPREF*2*dp->rhoa*dp->rhob*factor;
#else
  ds->df1000 += EPREF*4.0/3.0*pow(dp->rhoa,1.0/3.0)*dp->grada*dp->grada*factor;
  ds->df0010 += EPREF*pow(dp->rhoa,4.0/3.0)*2*dp->grada*factor;
#endif
}

static void
example_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
#if LEVEL==2
  ds->df1000 +=  EPREF*2*dp->rhoa*dp->rhob*factor;
  ds->df2000 +=  EPREF*2*dp->rhob*factor;
  ds->df1100 +=  EPREF*2*dp->rhoa*factor;
#else
  const real q = dp->grada*dp->grada;
  ds->df1000 +=  EPREF*4.0/3.0*pow(dp->rhoa,1.0/3.0)*q*factor;
  ds->df0010 +=  EPREF*pow(dp->rhoa,4.0/3.0)*2*dp->grada*factor;
  ds->df2000 +=  EPREF*4.0/9.0 *pow(dp->rhoa,-2.0/3.0)*q*factor;
  ds->df1010 +=  EPREF*4.0/3.0 *pow(dp->rhoa,1.0/3.0)*2*dp->grada*factor;
  ds->df0020 +=  EPREF*pow(dp->rhoa,4.0/3.0)*2*factor;
#endif
}

/* example_third:
   Test functional derivatives.
*/
static void
example_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
#if LEVEL==2
  ds->df1000 +=  EPREF*2*dp->rhoa*dp->rhob*factor;
  ds->df2000 +=  EPREF*2*dp->rhob*factor;
  ds->df1100 +=  EPREF*2*dp->rhoa*factor;
  ds->df2100 +=  EPREF*2*factor;
#else
  const real q = dp->grada*dp->grada;
  ds->df1000 +=  EPREF*4.0/3.0*pow(dp->rhoa,1.0/3.0)*q*factor;
  ds->df0010 +=  EPREF*pow(dp->rhoa,4.0/3.0)*2*dp->grada*factor;
  ds->df2000 +=  EPREF*4.0/9.0 *pow(dp->rhoa,-2.0/3.0)*q*factor;
  ds->df1010 +=  EPREF*4.0/3.0 *pow(dp->rhoa, 1.0/3.0)*2*dp->grada*factor;
  ds->df0020 +=  EPREF*pow(dp->rhoa,4.0/3.0)*2*factor;

  ds->df3000 += -EPREF*8.0/27.0*pow(dp->rhoa,-5.0/3.0)*q*factor;
  ds->df2010 +=  EPREF*4.0/9.0 *pow(dp->rhoa,-2.0/3.0)*2*dp->grada*factor;
  ds->df1020 +=  EPREF*4.0/3.0 *pow(dp->rhoa, 1.0/3.0)*2*factor;
#endif
}
