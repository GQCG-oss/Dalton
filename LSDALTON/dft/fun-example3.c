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
/* fun-test.c:
   implementation of a test GGA-class functional.
   This is a third example functional.
   (c) Pawel Salek, pawsa@theochem.kth.se, aug 2001
   NOTE:
   this file may seem unnecessarily complex but the structure really pays off
   when implementing multiple functionals depending on different parameters.
*/

#include <math.h>
#include <stddef.h>

#define __CVERSION__

#include "lsdalton_functionals.h"

/* INTERFACE PART */
static integer example3_isgga(void) { return 1; }
static integer example3_read(const char* conf_line, real *hfweight);
static real example3_energy(const FunDensProp* dp);
static void example3_first(FunFirstFuncDrv *ds,   real factor, const FunDensProp* dp);
static void example3_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp);
static void example3_third(FunThirdFuncDrv *ds,   real factor, const FunDensProp* dp);

Functional Example3Functional = {
  "Example3",         /* name */
  example3_isgga,     /* gga-corrected */
  example3_read, 
  NULL,              /* reporter */
  example3_energy, 
  example3_first,
  example3_second,
  example3_third
};

/* IMPLEMENTATION PART */
static integer
example3_read(const char* conf_line, real *hfweight)
{
  fun_set_hf_weight(0.0);
  return 1;
}

#define LEVEL 3
#if LEVEL==2
static const real EPREF= -1e-2;
#else
static const real EPREF= -5e-2;
#endif

static real
example3_energy(const FunDensProp* dp)
{
#if LEVEL==2
  return EPREF*(pow(dp->grada,4.0)+pow(dp->gradb,4.0)); 
#else
  return EPREF*(pow(dp->grada,1.7)+pow(dp->gradb,1.7));
#endif
}
/* example_first:
   derivatives with respect to rho_alpha, and zeta=|\nabla\grad_alpha|^2
 */
static void
example3_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
#if LEVEL==2
  ds->df0010 += EPREF*3*pow(dp->grada,3.0)*factor;
#else
  ds->df0010 += EPREF*1.7*pow(dp->grada,0.7)*factor;
#endif
}

/* example3_second:
   derivatives with respect to rho_alpha, and zeta=|\nabla\grad_alpha|
   (OBS: inconsistency!)
*/
static void
example3_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
#if LEVEL==2
  ds->df0010 += EPREF* 4*pow(dp->grada,3.0)*factor;
  ds->df0020 += EPREF*12*pow(dp->grada,2.0)*factor;
#else
  ds->df0010 +=  EPREF*1.7*pow(dp->grada,0.7)*factor;
  ds->df0020 +=  EPREF*1.7*0.7*pow(dp->grada,-0.3)*factor;
#endif
}

/* example_third:
   Test functional derivatives.
   derivatives with respect to rho_alpha, and zeta=|\nabla\grad_alpha|^2
*/
static void
example3_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
#if LEVEL==2
  ds->df0010 += EPREF* 4*pow(dp->grada,3.0)*factor;
  ds->df0020 += EPREF*12*pow(dp->grada,2.0)*factor;
#else
  ds->df0010 +=  EPREF*1.7*pow(dp->grada,0.7)*factor;
  ds->df0020 +=  EPREF*1.7*0.7*pow(dp->grada,-0.3)*factor;
  ds->df0030 += -EPREF*1.7*0.7*0.3*pow(dp->grada,-1.3)*factor;
#endif
}
