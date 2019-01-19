/*


!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!

!

*/
/* DFT callback test program.
   (c) Pawel Salek, pawsa@theochem.kth.se, feb 2002
*/
#include "general.h"
#define __CVERSION__

#include "integrator.h"

static void
test_callback(DftGrid* grid, real rho, real* res)
{
    sum[0] = sum[0] + rho*grid->curr_weight;
    sum[1] = sum[1] + dftenergy_(rho, grid->ngrad);
}

real
dft_test_(real* cmo, real* work, integer* lwork, integer* iprint)
{
    real res[2] = {0, 0};

    cbdata[0].callback = (DftCallback)test_callback;
    cbdata[0].cb_data  = res;
    dft_integrate(cmo, work, lwork, iprint, cbdata, ELEMENTS(cbdata));
    fort_print("electrons: %20.14f energy: %20.10g", res[0], res[1]);
}

