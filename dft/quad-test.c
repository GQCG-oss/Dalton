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

