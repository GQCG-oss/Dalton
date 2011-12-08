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
/* DFT callback test program.
   (c) Pawel Salek, pawsa@theochem.kth.se, feb 2002
*/
#define __CVERSION__

#include "integrator.h"

static void
test_callback(DftGrid* grid, real rho, real* res)
{
    sum[0] = sum[0] + rho*grid->curr_weight;
    sum[1] = sum[1] + dftenergy_(rho, grid->ngrad);
}

real
dft_test_(real* cmo, real* work, int* lwork)
{
    int norbt2 = inforb_.norbt*inforb_.norbt;
    real res[2] = {0, 0};

    cbdata[0].callback = (DftCallback)test_callback;
    cbdata[0].cb_data  = res;
    dft_integrate(cmo, work, lwork, cbdata, ELEMENTS(cbdata));
    fort_print("electrons: %20.14f energy: %20.10g", res[0], res[1]);
}

