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
/* general.c:
   (c) Pawel Salek, pawsa@theochem.kth.se, 2001-08-02
   NOTES: Adding new functionals:
   a. use fun-slater.c as template.
   b. add 'extern Functional MyFunctional;' to functionals.h
   c. add '&MyFunctional' to available_functionals below.
   d. have a beer. Or some crackers, if you prefer.
*/

/* strictly conform to XOPEN ANSI C standard */
#if !defined(SYS_DEC)
/* XOPEN compliance is missing on old Tru64 4.0E Alphas */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1
#endif

/* Use BSD's strncasecmp(); if there is a platform that has no strncasecmp()
 * ask pawsa@theochem.kth.se for replacement */
#define _BSD_SOURCE 1

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#define __CVERSION__

#include "general.h"
#include "integrator.h"
#include "functionals.h"
#include "inforb.h"

/* C-wide constants */
const integer ZEROI = 0,   ONEI = 1, THREEI = 3, FOURI = 4;
const real ZEROR = 0.0, ONER = 1.0, TWOR = 2.0, FOURR = 4.0;


/* stub subroutines for the functional code */
extern real dftgethf_(void);
extern void dftsethf_(real *w);
extern real dftgetmp2_(void);
extern void dftsetmp2_(real *w);
extern void dftsetcam_(const real *w, const real *b);

static real
dal_get_hf_weight(void) { return dftgethf_(); }
static void
dal_set_hf_weight(real w) { dftsethf_(&w); }

static real
dal_get_mp2_weight(void) { return dftgetmp2_(); }
static void
dal_set_mp2_weight(real w) { dftsetmp2_(&w); }

static void
dal_set_cam_param(int cnt, const real *w, const real *be) {
  dftsetcam_(w, be);
}

/* =================================================================== */
/* dftinput:

   read DFT functional from given line. The calling convention assumes
   Sun-style passing parameters. ATLAS linear algebra package
   http://www.netlib.org/atlas/ or http://math-atlas.sourceforge.net/
   contains an elaborate discuttion of character type variable passing
   conventions, in particular little bit below
   http://math-atlas.sourceforge.net/errata.html#RH7.0
*/
static char* DftConfString = NULL;
int
FSYM(dftsetfunc)(const char* line, integer * inperr, int len)
{
    int i, off;

    /* set the functional code printf function and HF weight setting
       functions to the dalton version that appends the output to the
       DALTON.OUT file. */
    fun_printf        = fort_print;
    fun_set_hf_weight = dal_set_hf_weight;
    fun_get_hf_weight = dal_get_hf_weight;
    fun_set_mp2_weight = dal_set_mp2_weight;
    fun_get_mp2_weight = dal_get_mp2_weight;
    fun_set_cam_param = dal_set_cam_param;

    for(i=len-1; i>=0 && isspace((int)line[i]); i--)
        ;
    if(DftConfString) free(DftConfString);
    i++;
    for(off=0; line[off] && isspace((int)line[off]); off++)
        ;
    DftConfString = malloc(i+1-off);
    strncpy(DftConfString, line+off, i-off); 
    DftConfString[i-off] = '\0';
    
    switch(fun_select_by_name(DftConfString)) {
    case FUN_OK: return 1; /* SUCCESS! */
    case FUN_UNKNOWN:
        fort_print("Unknown functional '%s'. Aborting.\n", DftConfString);
        dftlistfuncs_();
        break;
    case FUN_CONF_ERROR:
        fort_print("Functional configuration '%s' is not understood. "
                   "Aborting.\n", DftConfString);
        break;
    }
    (*inperr)++;
    return 0; /* failed to find */
}

/* =================================================================== */
/* transformations of functional derivatives from from unrestricted
 * ones to closed shell case. */
real
dftene_(const real *rho, const real *grad)
{
    FunDensProp dp;
    dp.rhoa  = dp.rhob  = *rho  *0.5;
    dp.grada = dp.gradb = *grad *0.5;
    dp.gradab = dp.grada*dp.gradb;
    return selected_func->func(&dp);
}

void
dftptf0_(real *rho, real *grad, real *wght, real *vx)
{
    FunFirstFuncDrv drvs;
    FunDensProp dp;
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
dftpot0_(FirstDrv *ds, const real* weight, const FunDensProp* dp)
{
    FunFirstFuncDrv drvs;
    drv1_clear(&drvs);
    selected_func->first(&drvs, *weight, dp);
    ds->fR = drvs.df1000;
    ds->fZ = drvs.df0010 + 0.5*drvs.df00001* (dp->grada+dp->gradb);
}

void
dftpot1_(SecondDrv *ds, const real* w, const FunDensProp* dp, 
         const integer* triplet)
{
    
    FunSecondFuncDrv drvs;

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
dftpot2_(ThirdDrv *ds, real factor, const FunDensProp* dp, integer isgga,
         integer triplet)
{
    FunThirdFuncDrv drvs;

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

void
dftpot3ab_(FourthDrv *ds, const real *factor, const FunDensProp* dp,
	   const integer *isgga)
{
     FunFourthFuncDrv drvs;
     
     /* Initialize drvs */ 
     drv4_clear(&drvs);
     selected_func->fourth(&drvs, *factor, dp);

     /* Transform derivatives (drvs -> ds) */

     ds->fR = 0.5*(drvs.df1000 + drvs.df0100);
     ds->fRR = 0.25*(drvs.df2000 + 2.0*drvs.df1100 + drvs.df0200);
     ds->fRRR = 0.125*(drvs.df3000 + 3.0*drvs.df2100 + 
                       3.0*drvs.df1200 + drvs.df0300);
     ds->fRRRR = 0.0625*(drvs.df4000 + 4.0*drvs.df3100 + 
                         6.0*drvs.df2200 + 4.0*drvs.df1300 + drvs.df0400);
     
     if (*isgga) {
         real igroa, igrob;
         real igroa_p2, igrob_p2;
         real igroa_p3, igrob_p3;
         real igroa_p4, igrob_p4;
         real igroa_p5, igrob_p5;
         real igroa_p6, igrob_p6;
         real igroa_p7, igrob_p7;
         
         /* 1/groa, 1/grob and its powers */
         igroa = 1.0/(dp->grada);
         igrob = 1.0/(dp->gradb);
         
         igroa_p2=igroa*igroa;
         igroa_p3=igroa_p2*igroa;
         igroa_p4=igroa_p3*igroa;
         igroa_p5=igroa_p4*igroa;
         igroa_p6=igroa_p5*igroa;
         igroa_p7=igroa_p6*igroa;
         
         igrob_p2=igrob*igrob;
         igrob_p3=igrob_p2*igrob;
         igrob_p4=igrob_p3*igrob;
         igrob_p5=igrob_p4*igrob;
         igrob_p6=igrob_p5*igrob;
         igrob_p7=igrob_p6*igrob;
         
         
         /* 1st order */
         
         ds->fZ = 0.125*(igroa*drvs.df0010+2.0*drvs.df00001+igrob*drvs.df0001);
         
         /* 2nd order */
         
         ds->fRZ = 0.0625*(2*drvs.df01001 + 2*drvs.df10001 + drvs.df0110*igroa
                           +drvs.df1010*igroa + drvs.df0101*igrob
                           + drvs.df1001*igrob);
     
         ds->fZZ = 0.015625*(4*drvs.df00002 + 4*drvs.df00101*igroa + 
                     drvs.df0020*igroa_p2 - drvs.df0010*igroa_p3 + 
                     4*drvs.df00011*igrob + 2*drvs.df0011*igroa*igrob + 
                     drvs.df0002*igrob_p2 - drvs.df0001*igrob_p3);

      /* 3rd order */
     
     ds->fRRZ = 0.03125*(2*drvs.df02001 + 4*drvs.df11001 + 2*drvs.df20001 +
		     drvs.df0210*igroa + 2*drvs.df1110*igroa +
                     drvs.df2010*igroa + drvs.df0201*igrob +
                     2*drvs.df1101*igrob + drvs.df2001*igrob);

     ds->fRZZ = 0.0078125*(4*drvs.df01002 + 4*drvs.df10002 + 4*drvs.df01101*igroa + 
	             4*drvs.df10101*igroa + drvs.df0120*igroa_p2 +       
	             drvs.df1020*igroa_p2 - drvs.df0110*igroa_p3 -       
	             drvs.df1010*igroa_p3 + 4*drvs.df01011*igrob +       
	             4*drvs.df10011*igrob + 2*drvs.df0111*igroa*igrob +       
	             2*drvs.df1011*igroa*igrob + drvs.df0102*igrob_p2 +       
	             drvs.df1002*igrob_p2 - drvs.df0101*igrob_p3 -       
	             drvs.df1001*igrob_p3);

     ds->fZZZ = 0.001953125*(8*drvs.df00003 + 12*drvs.df00102*igroa + 
		     6*drvs.df00201*igroa_p2 - 6*drvs.df00101*igroa_p3 +       
		     drvs.df0030*igroa_p3 - 3*drvs.df0020*igroa_p4 +       
		     3*drvs.df0010*igroa_p5 + 12*drvs.df00012*igrob + 
		     12*drvs.df00111*igroa*igrob +3*drvs.df0021*igroa_p2*igrob-
		     3*drvs.df0011*igroa_p3*igrob + 
		     6*drvs.df00021*igrob_p2 + 3*drvs.df0012*igroa*igrob_p2 - 
		     6*drvs.df00011*igrob_p3 + drvs.df0003*igrob_p3 - 
		     3*drvs.df0011*igroa*igrob_p3 - 3*drvs.df0002*igrob_p4 + 
	  	     3*drvs.df0001*igrob_p5);

     /* 4th order */

     ds->fRRRZ = 0.015625*(2*drvs.df03001 + 6*drvs.df12001 + 6*drvs.df21001 + 
		     2*drvs.df30001 + drvs.df0310*igroa + 3*drvs.df1210*igroa +
		     3*drvs.df2110*igroa + drvs.df3010*igroa +
		     drvs.df0301*igrob + 3*drvs.df1201*igrob + 
		     3*drvs.df2101*igrob + drvs.df3001*igrob);
     
     ds->fRRZZ = 0.00390625*(4*drvs.df02002 + 8*drvs.df11002 + 4*drvs.df20002 +
	             4*drvs.df02101*igroa + 8*drvs.df11101*igroa + 
	             4*drvs.df20101*igroa + drvs.df0220*igroa_p2 + 
 	             2*drvs.df1120*igroa_p2 + drvs.df2020*igroa_p2 - 
	             drvs.df0210*igroa_p3 - 2*drvs.df1110*igroa_p3 - 
	             drvs.df2010*igroa_p3 + 4*drvs.df02011*igrob + 
	             8*drvs.df11011*igrob + 4*drvs.df20011*igrob + 
	             2*drvs.df0211*igroa*igrob + 4*drvs.df1111*igroa*igrob + 
	             2*drvs.df2011*igroa*igrob + drvs.df0202*igrob_p2 + 
	             2*drvs.df1102*igrob_p2 + drvs.df2002*igrob_p2 - 
	             drvs.df0201*igrob_p3 - 2*drvs.df1101*igrob_p3 - 
	             drvs.df2001*igrob_p3);

     ds->fRZZZ = 0.0009765625*(8*drvs.df01003 + 8*drvs.df10003 + 12*drvs.df01102*igroa + 
		     12*drvs.df10102*igroa + 6*drvs.df01201*igroa_p2 + 
		     6*drvs.df10201*igroa_p2 - 6*drvs.df01101*igroa_p3 + 
		     drvs.df0130*igroa_p3 - 6*drvs.df10101*igroa_p3 + 
		     drvs.df1030*igroa_p3 - 3*drvs.df0120*igroa_p4 - 
		     3*drvs.df1020*igroa_p4 + 3*drvs.df0110*igroa_p5 + 
		     3*drvs.df1010*igroa_p5 + 12*drvs.df01012*igrob + 
		     12*drvs.df10012*igrob + 12*drvs.df01111*igroa*igrob + 
		     12*drvs.df10111*igroa*igrob +3*drvs.df0121*igroa_p2*igrob+
		     3*drvs.df1021*igroa_p2*igrob - 
		     3*drvs.df0111*igroa_p3*igrob-3*drvs.df1011*igroa_p3*igrob+
		     6*drvs.df01021*igrob_p2 + 6*drvs.df10021*igrob_p2 + 
		     3*drvs.df0112*igroa*igrob_p2 + 3*drvs.df1012*igroa*igrob_p2 -
		     6*drvs.df01011*igrob_p3 + drvs.df0103*igrob_p3 - 
		     6*drvs.df10011*igrob_p3 + drvs.df1003*igrob_p3 - 
		     3*drvs.df0111*igroa*igrob_p3-3*drvs.df1011*igroa*igrob_p3-
		     3*drvs.df0102*igrob_p4 - 3*drvs.df1002*igrob_p4 + 
		     3*drvs.df0101*igrob_p5 + 3*drvs.df1001*igrob_p5);

     ds->fZZZZ = (1.0/4096.0)*(16*drvs.df00004 + 32*drvs.df00103*igroa + 
		     24*drvs.df00202*igroa_p2 - 24*drvs.df00102*igroa_p3 + 
		     8*drvs.df00301*igroa_p3 - 24*drvs.df00201*igroa_p4 + 
		     drvs.df0040*igroa_p4 + 24*drvs.df00101*igroa_p5 - 
		     6*drvs.df0030*igroa_p5 + 15*drvs.df0020*igroa_p6 - 
		     15*drvs.df0010*igroa_p7 + 32*drvs.df00013*igrob + 
		     48*drvs.df00112*igroa*igrob + 24*drvs.df00211*igroa_p2*igrob - 
		     24*drvs.df00111*igroa_p3*igrob +
		     4*drvs.df0031*igroa_p3*igrob - 12*drvs.df0021*igroa_p4*igrob + 
		     12*drvs.df0011*igroa_p5*igrob +
		     24*drvs.df00022*igrob_p2 + 24*drvs.df00121*igroa*igrob_p2+
		     6*drvs.df0022*igroa_p2*igrob_p2 - 6*drvs.df0012*igroa_p3*igrob_p2 -
		     24*drvs.df00012*igrob_p3 + 8*drvs.df00031*igrob_p3 - 
		     24*drvs.df00111*igroa*igrob_p3 + 
		     4*drvs.df0013*igroa*igrob_p3 - 6*drvs.df0021*igroa_p2*igrob_p3 +
		     6*drvs.df0011*igroa_p3*igrob_p3 - 
		     24*drvs.df00021*igrob_p4 + drvs.df0004*igrob_p4 - 
		     12*drvs.df0012*igroa*igrob_p4 + 24*drvs.df00011*igrob_p5 -
		     6*drvs.df0003*igrob_p5 + 12*drvs.df0011*igroa*igrob_p5 + 
		     15*drvs.df0002*igrob_p6 - 15*drvs.df0001*igrob_p7);
     }
}


void
vxcfab_(real *rhoa, real *rhob, real *grada, real *gradb, 
        real *gradab, real *wght, real *vx)
{
  FunFirstFuncDrv drvs;
  FunDensProp dp;
 
  dp.rhoa   = *rhoa;
  dp.rhob   = *rhob; 
  dp.grada  = *grada;
  dp.gradb  = *gradb;
  dp.gradab = *gradab;

  if(dp.rhoa<1e-13)   dp.rhoa   = 1e-13;
  if(dp.rhob<1e-13)   dp.rhob   = 1e-13;
  if(dp.grada<1e-13)  dp.grada  = 1e-13;
  if(dp.gradb<1e-13)  dp.gradb  = 1e-13;
  if(dp.gradab<1e-13) dp.gradab = 1e-13;

  drv1_clear(&drvs);
  selected_func->first(&drvs, *wght, &dp);

  vx[0] = drvs.df1000;
  vx[1] = drvs.df0100;
  vx[2] = drvs.df0010; 
  vx[3] = drvs.df0001;
  vx[4] = drvs.df00001; 
}



/* =================================================================== */
/*    DFT density evaluators for restricted and unrestricted cases.    */
/*    evaluate density properties                                      */
/* =================================================================== */
void
dft_dens_restricted(DftDensity* dens, FunDensProp* dp, DftGrid* grid,
                    real* tmp_vec)
{
    /* Note that dens->dmata is really the total density, i.e rho_a+rho_b */
    /* since this is the convention in the rest of dalton */
    real rho, ngrad;
    getrho_(dens->dmata, grid->atv, &rho, tmp_vec, &grid->dfthri);
    dp->rhoa = dp->rhob = 0.5*rho; 
    /* transform grad vectors to molecular orbitals */
    if(rho>grid->dfthr0 && grid->dogga) {
        /* compute only half density gradient, i.e only grad_alpha. */
      FSYM(dgemv)("T", &inforb_.nbast, &THREEI, &ONER,
		  &grid->atv[inforb_.nbast], &inforb_.nbast, tmp_vec,
		  &ONEI, &ZEROR, grid->grada, &ONEI);
        
        ngrad = sqrt(grid->grada[0]*grid->grada[0]+
                     grid->grada[1]*grid->grada[1]+
                     grid->grada[2]*grid->grada[2]);
        dp->grada = dp->gradb = ngrad;
        dp->gradab= dp->grada*dp->gradb;
    }
}

void
dft_dens_unrestricted(DftDensity* dens, FunDensProp* dp, DftGrid* grid,
                      real* tmp_vec)
{
    real *tmpa = tmp_vec, *tmpb = tmp_vec + inforb_.nbast;
    getrho_(dens->dmata, grid->atv, &dp->rhoa, tmpa, &grid->dfthri);
    getrho_(dens->dmatb, grid->atv, &dp->rhob, tmpb, &grid->dfthri);

    if(dp->rhoa<1e-40)          dp->rhoa = 1e-40;
    if(dp->rhob<dp->rhoa*1e-15) dp->rhob = dp->rhoa*1e-15;

    if( (dp->rhoa+dp->rhob)>grid->dfthr0 && grid->dogga) {
        /* transform grad vectors to molecular orbitals */
        dgemv_("T", &inforb_.nbast, &THREEI, &TWOR,
               &grid->atv[inforb_.nbast], &inforb_.nbast, tmpa,
               &ONEI, &ZEROR, grid->grada, &ONEI);
        dp->grada = sqrt(grid->grada[0]*grid->grada[0]+
                         grid->grada[1]*grid->grada[1]+
                         grid->grada[2]*grid->grada[2]);
        dgemv_("T", &inforb_.nbast, &THREEI, &TWOR,
               &grid->atv[inforb_.nbast], &inforb_.nbast, tmpb,
               &ONEI, &ZEROR, grid->gradb, &ONEI);
        dp->gradb = sqrt(grid->gradb[0]*grid->gradb[0]+
                         grid->gradb[1]*grid->gradb[1]+
                         grid->gradb[2]*grid->gradb[2]);

        dp->gradab = grid->grada[0]*grid->gradb[0]+
            grid->grada[1]*grid->gradb[1]+grid->grada[2]*grid->gradb[2];
  }
}

#if defined(VAR_MPI)
#include <mpi.h>
#include <our_extra_mpi.h>
#define MASTER_NO 0
/* =================================================================== */
/* General parallel routines.
 * dft_lifesupport_sync - one-time sync of basic common block data
 *                   that is crucial for initialization.
 * dft_wake_slaves - wakes slaves for a specified DFT evaluator (perhaps this
 *                   can be generalized?). It gets a registered property
 *                   evaluator address as an argument and broadcasts a
 *                   message that will call corresponding slave code
 *                   via dft_cslave().
 * FIXME: use rather a register approach that does not require modification
 * of this file every time new property evaluator is added.
 */

struct {
    DFTPropEvalMaster master_func;
    DFTPropEvalSlave  slave_func;
} PropEvaluatorList[] = {
#if 0
    { (DFTPropEvalMaster)dft_kohn_sham_,   dft_kohn_sham_slave   },
    { (DFTPropEvalMaster)dft_kohn_shamab_, dft_kohn_shamab_slave },
    { (DFTPropEvalMaster)dft_lin_resp_,    dft_lin_resp_slave    },
    { (DFTPropEvalMaster)dft_lin_respab_,  dft_lin_respab_slave  },
    { (DFTPropEvalMaster)dft_mol_grad_,    dft_mol_grad_slave    },
#endif
    { (DFTPropEvalMaster)FSYM2(dft_kohn_shamab),dft_kohn_shamab_slave },
    { (DFTPropEvalMaster)FSYM2(dft_lin_respab), dft_lin_respab_slave  },
    { (DFTPropEvalMaster)FSYM2(dft_lin_respf),  dft_lin_respf_slave },
    { (DFTPropEvalMaster)dftqrcf_,              dft_qr_resp_slave   },
    { (DFTPropEvalMaster)FSYM2(dft_qr_respons), dft_qrbl_slave      },
    { (DFTPropEvalMaster)FSYM(dftcrcf),         dft_cr_resp_slave   },
    { (DFTPropEvalMaster)FSYM(numdso),          numdso_slave        }, 
    { (DFTPropEvalMaster)FSYM2(dft_kohn_shamab_b), dft_kohn_shamab_b_slave}, 
    { (DFTPropEvalMaster)FSYM2(dft_lin_respab_b),  dft_lin_respab_b_slave}
};

/* mpi_sync_data:
   sync given set of data with slaves.
*/
void
mpi_sync_data(const SyncData* data, int count)
{
    int i;
    for(i=0; i<count; i++) {
        MPI_Bcast(data[i].data, data[i].count, data[i].type,
                  0, MPI_COMM_WORLD);
    }
}


void
dft_wake_slaves(DFTPropEvalMaster evaluator)
{
    static integer iprtyp = 5; /* magic DFT/C number */
    static integer iprint = 0;
    int id, mynum;

    MPI_Comm_rank(MPI_COMM_WORLD, &mynum);
    if(mynum != 0)
        return; /* slaves do not wake up other slaves */

    for(id=0; 
        id<ELEMENTS(PropEvaluatorList) && 
            PropEvaluatorList[id].master_func != evaluator;
        id++)
        ;
        
    if(id>=ELEMENTS(PropEvaluatorList)) {
        /* this would really be an programming error.... */
        fprintf(stderr, "Evaluator not registered. No slaves activated.\n");
        return;
    }
    /* ignore MPI errors */
#ifdef VAR_INT64      
/*        
    printf("waking up slaves with iprtyp %lld and iprint %lld \n",iprtyp, iprint);
    printf(" id is %lld \n",id);
*/
#else
/*        
    printf("waking up slaves with iprtyp %d and iprint %d \n",iprtyp, iprint);
    printf(" id is %d \n",id);
*/
#endif
      
    MPI_Bcast(&iprtyp,1, fortran_MPI_INT, MASTER_NO, MPI_COMM_WORLD);
    MPI_Bcast(&iprint,1, fortran_MPI_INT, MASTER_NO, MPI_COMM_WORLD);
    MPI_Bcast(&id,    1, MPI_INT, MASTER_NO, MPI_COMM_WORLD);
    /*
    printf("id %d bcast",id);
    */
    FSYM(dftintbcast)();
}

/* dft_cslave:
   a slave task handler. Receives the task ID and calls apropriate
   property evaluator.
   We broadcast also some basic properties and common block data
   that slaves should absolutely know but for some reason do not.
*/
void
FSYM2(dft_cslave)(real* work, integer*lwork,integer*iprint)
{
    int rank, size;
    if(MPI_Comm_rank(MPI_COMM_WORLD, &rank) ||
       MPI_Comm_size(MPI_COMM_WORLD, &size)) printf("MPI error\n");
    else {
        int id;
        MPI_Bcast(&id,1,MPI_INT, MASTER_NO, MPI_COMM_WORLD);
/*        
        printf("done id bcast: id = %i mynum = %i; size = %i .\n", id, rank, size);
*/        
        FSYM(dftintbcast)();
        (PropEvaluatorList[id].slave_func)(work, lwork, iprint);
    }
}
#else
void
FSYM2(dft_cslave)(real* work,integer*lwork,integer*iprint)
{
   fort_print("DFT slave called but does nothing now.");
}
#endif /* VAR_MPI */


extern void quit_(const char* str, int len);
void
dalton_quit(const char* format, ...)
{
    char line[128];
    int len;
    va_list a;
 
    va_start(a, format);
    vsnprintf(line, sizeof(line), format, a);
    va_end(a);
    len = strlen(line);
    quit_(line, len);
}

/* Helper functions. Could be bracketed with #ifdef DEBUG or something */
extern void FSYM2(fort_wrt)(const char* str, const integer* len, int ln);
int
fort_print(const char* format, ...)
{
    char line[128];
    integer len;
    integer l;
    va_list a;

    va_start(a, format);
    vsnprintf(line, sizeof(line), format, a);
    va_end(a);
    len = l = strlen(line);
    FSYM2(fort_wrt)(line, &len, l);
    return len;
}

/* ===================================================================
 * Parallel support.
 * =================================================================== */
#ifdef VAR_MPI
#include <mpi.h>

/* dftfuncsync synchronises the selected functional between all
 * nodes. */
void
FSYM(dftfuncsync)(integer *mynum, integer *nodes)
{
    static int done = 0;
    int len = DftConfString ? strlen(DftConfString) + 1: 0;
    if(done) return;
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(len>0) {
        integer res;
        char *line = malloc(len);
        if(*mynum == 0)
            strcpy(line, DftConfString);
        MPI_Bcast(line, len, MPI_CHAR, 0, MPI_COMM_WORLD);
/*        
        printf("my line is %s its length is %i",line, len);
*/
        if(*mynum != 0) {
            integer res;
            FSYM(dftsetfunc)(line, &res, len);
        }
        free(line);
    }
    done = 1;
}

#endif
 
