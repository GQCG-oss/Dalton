/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* functionals.c:
   Program for computing and testing functional routines in the DFT module.
   (c) Pawel Salek, pawsa@theochem.kth.se, 2001-08-02

   The test build can be done by:
   g77 -O -DTEST_BUILD functionals.c -o functionals -L. -ldft -lm
   or 
   cc -O -DTEST_BUILD functionals.c -o functionals -L. -ldft -lm -lg2c

   NOTES: this file is short but in a separate file to reduce the
   number of dependences and be able to easily compile the code for
   the TEST_BUILD. 

   NOTES:
   the closed shell assumptions make following derivatives equivalent:

   3000 = 0300
   2100 = 1200
   2010 = 0201
   1110 = 1101
   0210 = 2001
   1020 = 0102
   0120 = 1002
   0030 = 0003
   1011 = 0111
   0021 = 0012
 */

#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define __CVERSION__

#include "general.h"
#include "functionals.h"
#include "inforb.h"

Functional* available_functionals[] = {
    /* generic functionals */
    &BeckeFunctional,
    &Example2Functional,
    &Example3Functional,
    &Example4Functional,
    &Example5Functional,
    &Example6Functional,
    &Example7Functional,
    &Example8Functional,
    &Example9Functional,
    &ExampleFunctional,
    &KTFunctional,
    &LB94Functional,
    &LYPFunctional,
    &OPTXFunctional,
    &P86cFunctional,
    &PW86xFunctional,
    &PW91cFunctional,
    &PW92Functional,
    &PWggaIIc2Functional,
    &PWggaIIcFunctional,
    &PWggaIIxFunctional,
    &PZ81Functional,
    &PbecFunctional,
    &PbexFunctional,
    &SlaterFunctional,
    &VWN3Functional,
    &VWN5Functional,
    &VWNFunctional,
    &XAlphaFunctional,
    /* mixed functionals */
    &B3LYPFunctional,
    &B3LYPGaussFunctional,
    &B3P86Functional,
    &B3P86GFunctional,
    &BLYPFunctional,
    &BP86Functional,
    &BPW91Functional,
    &Camb3lypFunctional,
    &GGAKeyFunctional,
    &KT1Functional,
    &KT2Functional,
    &KT3Functional,
    &LDAFunctional,
    &OLYPFunctional,
    &PBE0Functional,
    &PBEFunctional,
    &SVWN3Functional,
    &SVWN5Functional,
    NULL
};
Functional* selected_func = &LDAFunctional;
int (*fun_printf)(const char *fmt, ...) = printf;

static void set_hf_weight(real w)         {}
static real get_hf_weight(void)           {return 0;}
static void set_cam_param(real w, real b) {}

void (*fun_set_hf_weight)(real w)         = set_hf_weight;
real (*fun_get_hf_weight)(void)           = get_hf_weight;
void (*fun_set_cam_param)(real w, real b) = set_cam_param;

/* =================================================================== */
enum FunError
fun_select_by_name(const char *conf_string)
{
    int ok, i;
    char func_name[20];

    sscanf(conf_string,"%20s", func_name);
    for(i=0; available_functionals[i]; i++)
        if(strcasecmp(available_functionals[i]->name, func_name)==0) {
            selected_func = available_functionals[i];
            ok = selected_func->read ?
                selected_func->read(conf_string+strlen(func_name)) : 1;
            return ok ? FUN_OK : FUN_CONF_ERROR;
        }
    return FUN_UNKNOWN;
}

void
drv1_clear(FirstFuncDrv* gga)
{
    gga->df1000 = gga->df0100 = gga->df0010 = gga->df0001 = gga->df00001 = 0;
}

void
drv2_clear(SecondFuncDrv* gga)
{
    gga->df1000 = gga->df0100 = gga->df0010 = gga->df0001 = gga->df00001 = 0;
    gga->df2000 = gga->df0200 = gga->df0020 = gga->df0002 = gga->df00002 = 0;
    gga->df1010 = gga->df0101 = gga->df1001 = gga->df0110 = 0;
    gga->df1100 = gga->df0011 = gga->df10001 = gga->df01001 = 0;   
    gga->df00101= gga->df00011 = 0;
}

void
drv3_clear(ThirdFuncDrv* gga)
{
    gga->df1000 = gga->df0100 = gga->df0010 = gga->df0001 = gga->df00001 = 0;
    gga->df2000 = gga->df0200 = gga->df0020 = gga->df0002 = 0;
    gga->df1100 = gga->df0011 = 0;
    gga->df1010 = gga->df1001 = gga->df0101 = gga->df0110 = 0;
    gga->df10001 = gga->df01001 = 0;    

    gga->df3000 = gga->df0300 = gga->df0030 = gga->df0003 = 0;
    gga->df2100 = gga->df1200 = gga->df0012 = gga->df0021 = 0;  
    gga->df2010 = gga->df2001 = gga->df0201 = gga->df0210 = 0;
    gga->df1020 = gga->df1002 = gga->df0102 = gga->df0120 = 0; 
    gga->df1110 = gga->df1101 = gga->df1011 = gga->df0111 = 0;
    gga->df20001 = gga->df02001 = gga->df11001 = 0; 
}

/* fortran (and not only) functional stub routines */
int
FSYM2(dft_isgga)(void)
{
    return selected_func->is_gga();
}

/* dftreport:
   report the data and check consistency of the input (FIXME: this should
   be done by a separate routine.  
*/
void
dftreport_(void)
{
    fun_printf("\n     This is a DFT calculation of type: %s",
               selected_func->name);
    if(selected_func->report)
        selected_func->report();

    /* dft_check_consistency(); */
}

void
dftlistfuncs_(void)
{
    int i;
    fun_printf("\nAvailable functionals:");
    for(i=0; available_functionals[i]; i++)
        fun_printf(available_functionals[i]->name);
}

real
dftene_(const real *rho, const real *grad)
{
    DftDensProp dp;
    dp.rhoa  = dp.rhob  = *rho  *0.5;
    dp.grada = dp.gradb = *grad *0.5;
    dp.gradab = dp.grada*dp.gradb;
    return selected_func->func(&dp); 
}

void
dftptf0_(real *rho, real *grad, real *wght, real *vx)
{
    FirstFuncDrv drvs;
    DftDensProp dp;
    dp.rhoa  = dp.rhob  = *rho *0.5;
    dp.grada = dp.gradb = *grad*0.5;
    if(dp.rhoa<1e-13) dp.rhoa = dp.rhob = 1e-13;
    if(dp.grada<1e-13) dp.grada = dp.gradb = 1e-13;
    dp.gradab = dp.grada*dp.gradb;
    drv1_clear(&drvs);
    selected_func->first(&drvs, *wght, &dp);
    vx[0] = drvs.df1000;
    vx[1] = drvs.df0010 + 0.5*drvs.df00001* (*grad);
}

void
dftpot0_(FirstDrv *ds, const real* weight, const DftDensProp* dp)
{
    FirstFuncDrv drvs;
    drv1_clear(&drvs);
    selected_func->first(&drvs, *weight, dp);
    ds->fR = drvs.df1000;
    ds->fZ = drvs.df0010 + 0.5*drvs.df00001* (dp->grada+dp->gradb);
}

#if NO_BACKWARD_COMP
#define DFTPOT1 dftpot1
#else
#define DFTPOT1 dftpot1_
#endif
void
DFTPOT1(SecondDrv *ds, const real* w, const DftDensProp* dp, 
        const int* triplet)
{
    
    SecondFuncDrv drvs;

    drv2_clear(&drvs);
    if(dp->rhoa + dp->rhob>1e-14)
        selected_func->second(&drvs, *w, dp);
    if (*triplet) { /* triplet */  
        ds->fZ  = drvs.df0010;
        ds->fG  = -0.5*drvs.df00001;
        ds->fRR = 0.5*(drvs.df2000 - drvs.df1100);
        ds->fRZ = 0.5*(drvs.df1010 - drvs.df1001); 
        ds->fRG = 0.0;  
        ds->fZZ = 0.5*(drvs.df0020 - drvs.df0011); 
        ds->fZG = 0.0; 
        ds->fGG = 0.0; 
    } else { /* singlet */
        ds->fR  = 0.5*(drvs.df1000 + drvs.df0100);
        ds->fZ  = drvs.df0010;
        ds->fRR = 0.5*(drvs.df2000 + drvs.df1100);
        ds->fRZ = 0.5*(drvs.df1010 + drvs.df1001);
        ds->fZZ = 0.5*(drvs.df0020 + drvs.df0011); 
        ds->fRG = 0.5*drvs.df10001;   
        ds->fZG = 0.5*drvs.df00101;   
        ds->fGG = 0.25*drvs.df00002; 
        ds->fG  = 0.5*drvs.df00001;  
    }
}



/* dftpot2_:
   computes third order derivatives of selected functional with respect
   to rho and zeta=|\nabla\rho|
*/
void
dftpot2_(ThirdDrv *ds, real factor, const DftDensProp* dp, int isgga,
         int triplet)
{
    ThirdFuncDrv drvs;

    drv3_clear(&drvs);
    selected_func->third(&drvs, factor, dp);
    /* transform to density from density_alpha derivatives below */
    /* This could be a separate function. */
    /* we treat singlet here */
    ds->fR   = (drvs.df1000);
    ds->fRR[0] = (drvs.df2000 + drvs.df1100);
    ds->fRR[1] = (drvs.df2000 - drvs.df1100);
    ds->fRRR[0] = (drvs.df3000 + 3*drvs.df2100);
    ds->fRRR[1] = (drvs.df3000 - drvs.df2100);

    if(isgga) { /* FORMULAE seem ok but were never really tested */
        real grada = dp->grada;
        real grada2= grada*grada;
        real grada3= grada2*grada;
	/* transform to closed shell, second time. Oh, I love this mess. */
        ds->fZ  = drvs.df0010/(2*grada);
        ds->fG = drvs.df00001;
        ds->fZZ[0]  = (drvs.df0020 + drvs.df0011)/(4*grada2) 
	    -drvs.df0010/(4*grada3);
        ds->fZZ[1]  = (drvs.df0020 - drvs.df0011)/(4*grada2) 
	    -drvs.df0010/(4*grada3);
        ds->fRZ[0]  = (drvs.df1010 + drvs.df1001)/(2*grada);
        ds->fRZ[1]  = (drvs.df1010 - drvs.df1001)/(2*grada);
        ds->fRG[0]  = 2*drvs.df10001;
        ds->fRG[1]  = 0.0;  
        ds->fRRZ[0][0] = (drvs.df2010+drvs.df2001+2*drvs.df1110)/(2*grada);
        ds->fRRZ[0][1] = (drvs.df2010+drvs.df2001-2*drvs.df1110)/(2*grada);
        ds->fRRZ[1][0] = (drvs.df2010-drvs.df2001)/(2*grada);
        ds->fRRZ[1][1] = ds->fRRZ[1][0];
        ds->fRRG[0] = drvs.df20001+drvs.df11001;
        ds->fRRG[1] = drvs.df20001-drvs.df11001;
        ds->fRRGX[0][0] = 2*(drvs.df20001+drvs.df11001);
        ds->fRRGX[1][1] = 2*(drvs.df20001-drvs.df11001);
        ds->fRRGX[1][0] = ds->fRRGX[0][1] = 0;
         
        ds->fRZZ[0][0] = (drvs.df1020+drvs.df0120+2*drvs.df1011)/(4*grada2) 
                      - (drvs.df1010+drvs.df1001)/(4*grada3);
        ds->fRZZ[0][1] = (drvs.df1020+drvs.df0120-2*drvs.df1011)/(4*grada2) 
                      - (drvs.df1010+drvs.df1001)/(4*grada3);
        ds->fRZZ[1][0] = (drvs.df1020-drvs.df0120)/(4*grada2) 
                      - (drvs.df1010-drvs.df1001)/(4*grada3);
        ds->fRZZ[1][1] = ds->fRZZ[1][0];
        ds->fZZZ[0] = ((drvs.df0030 + 3*drvs.df0021)/grada3 
                   -3*(drvs.df0020 + drvs.df0011)/(grada2*grada2)
                   +3*drvs.df0010/(grada3*grada2))/8.0; 
        ds->fZZZ[1] = ((drvs.df0030 - drvs.df0021)/grada3 
                   -(3*drvs.df0020 - drvs.df0011)/(grada2*grada2)
                   +3*drvs.df0010/(grada3*grada2))/8.0; 
    } else {
        ds->fZ = ds->fZZ[0] = ds->fZZ[1] = ds->fRZ[0] = ds->fRZ[1] = 0;
        ds->fRRZ[0][0] = ds->fRRZ[0][1]  = ds->fRRZ[1][0] = 0;
        ds->fRRZ[1][1] = ds->fRZZ[0][0] = ds->fRZZ[0][1] = 0;
        ds->fRZZ[1][0] = ds->fRZZ[1][1] = ds->fZZZ[0] = ds->fZZZ[1] = 0; 
        ds->fG = ds->fRG[0] = ds->fRG[1]= ds->fRRG[0] = ds->fRRG[1] = 0;
        ds->fRRGX[0][0] = ds->fRRGX[1][0] = 
	    ds->fRRGX[0][1] = ds->fRRGX[1][1] = 0;  
    }
}

int fun_true(void)  { return 1; }
int fun_false(void) { return 0; }
