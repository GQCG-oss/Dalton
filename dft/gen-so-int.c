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
/* gen-so-int.c:
   generate SO integrals in form of three NBAST^2 matrices.

   (c) Pawel Salek, 2003, pawsa@theochem.kth.se


   For LDA case, the formulas are:
   r = dphi/dx, s = dphi/dy, t = dphi/dz: orbital derivatives
   x component: (s_p t_q - t_p s_q) dF/drho
   y component: (t_p r_q - r_p t_q) dF/drho
   z component: (r_p s_q - s_p r_q) dF/drho
   Tested against functional F = rhoa^2+rhob^2
   On atom:
INTGRL
He atom
--------
    1    0  X  Y  Z
        2.    1    2    1    1
H     .0000     .0000    0.0000
    1    1    0
 0.7000000  1.000000
    1    1    0
 0.7000000  1.000000
For this basis set and functional, the expression for the energy is:
E=2 (alpha/pi)^1.5
and the (p_x p_y dE/drho) matrix element
V=-2 alpha(alpha/pi)^1.5
(2 comes from the fact that (dF/drho=2dF/drhoa) for closed shell).
*/

/* strictly conform to XOPEN ANSI C standard */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>

#define __CVERSION__
#include "integrator.h"
#include "functionals.h"

#include "inforb.h"

#if !defined(RESTRICT)
#define RESTRICT restrict
#endif

static const int DFTSO_DEBUG = 1;

/* the computed expressions are:
   iX[i+joff] += (aY[i]*aZ[j] - aZ[i]*aY[j])*drvs.fR;
   iY[i+joff] += (aZ[i]*aX[j] - aX[i]*aZ[j])*drvs.fR;
   iZ[i+joff] += (aX[i]*aY[j] - aY[i]*aX[j])*drvs.fR;

   but we can use symmetry of the problem and reduce the number of
   operations inside and only antisymmetrize later on.
*/

static void
so_ints_lda_cb(DftGrid* grid, real* excmat)
{
    FirstDrv drvs;
    int i, j, nbast = inforb_.nbast;
    real *RESTRICT aX = &grid->atv[inforb_.nbast];
    real *RESTRICT aY = &grid->atv[inforb_.nbast*2];
    real *RESTRICT aZ = &grid->atv[inforb_.nbast*3];

    dftpot0_(&drvs, &grid->curr_weight, &grid->dp);

    for(j=0; j<nbast; j++) {
        real *RESTRICT iX = excmat + j*nbast;
        real *RESTRICT iY = excmat + j*nbast+inforb_.n2basx;
        real *RESTRICT iZ = excmat + j*nbast+inforb_.n2basx*2;
        for(i=0; i<nbast; i++) {
            iX[i] += (aY[i]*aZ[j])*drvs.fR;
            iY[i] += (aZ[i]*aX[j])*drvs.fR;
            iZ[i] += (aX[i]*aY[j])*drvs.fR;
        }
    }
}

/* dftsoi_:
   computes DFT contribution to SO integrals for given reference density
   and saves them to usual file.
*/
void FSYM2(write_soi)(const real* soints, real* wrk, integer* lwrk);
void
FSYM(dftsoi)(real* cmo, real* work, integer *lwork, integer* iprfck)
{
    int nbast2, m, i, j;
    DftCallbackData cbdata[1];
    DftDensity dens = { dft_dens_restricted, NULL, NULL };
    struct tms starttm, endtm; clock_t utm;
    real electrons, *so_ints;


    nbast2     = inforb_.nbast*inforb_.nbast;
    dens.dmata = malloc(nbast2*sizeof(real));
    so_ints    = calloc(nbast2*3, sizeof(real));
    cbdata[0].callback = (DftCallback)so_ints_lda_cb;
    cbdata[0].cb_data  = so_ints;

    fort_print("Computing DFT_SO integrals....");
    times(&starttm);
    FSYM2(dft_get_ao_dens_mat)(cmo, dens.dmata, work, lwork);
    if(DFTSO_DEBUG) {
    outmat_(dens.dmata,&ONEI,&inforb_.nbast,&ONEI,&inforb_.nbast,
            &inforb_.nbast,&inforb_.nbast);
    }
    electrons = dft_integrate_ao(&dens, work, lwork, 1,0,0, 
				 cbdata, ELEMENTS(cbdata));
    for(m=0; m<3; m++) {
        for(i=0; i<inforb_.nbast; i++)
            for(j=0; j<i; j++) {
                real averag = (so_ints[i+j*inforb_.nbast+m*nbast2]-
                               so_ints[j+i*inforb_.nbast+m*nbast2]);
                so_ints[i+j*inforb_.nbast+m*nbast2] =  averag;
                so_ints[j+i*inforb_.nbast+m*nbast2] = -averag;
            }
    }
    if(DFTSO_DEBUG) {
    fort_print("THE DFT-SO matrix: X");
    outmat_(so_ints,&ONEI,&inforb_.nbast,&ONEI,&inforb_.nbast,
            &inforb_.nbast,&inforb_.nbast);
    fort_print("THE DFT-SO matrix: Y");
    outmat_(so_ints+inforb_.n2basx,&ONEI,&inforb_.nbast,&ONEI,&inforb_.nbast,
            &inforb_.nbast,&inforb_.nbast);
    fort_print("THE DFT-SO matrix: Z");
    outmat_(so_ints+inforb_.n2basx*2,&ONEI,&inforb_.nbast,&ONEI,&inforb_.nbast,
            &inforb_.nbast,&inforb_.nbast);
    }
    FSYM2(write_soi)(so_ints, work, lwork);
    free(so_ints);
    free(dens.dmata);
    times(&endtm);
    utm = endtm.tms_utime-starttm.tms_utime;
    fort_print("Electrons: %11.7f %9.1g: SO-int time: %9.1f s\n", 
               electrons, (electrons-2.0*inforb_.nrhft)/(2.0*inforb_.nrhft), 
               utm/(double)sysconf(_SC_CLK_TCK));
}
