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
/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* fun-tester.c:
   Program for testing functional routines in the DFT module.
   (c) Pawel Salek, pawsa@theochem.kth.se, 2001-10-15

   The test build can be done by:
   g77 -O fun-tester.c -o fun-tester -L. -ldft -lm
   or 
   cc -O  fun-tester.c -o fun-tester -L. -ldft -lm -lg2c

   NOTES: this file is short but in a separate file to reduce the
   number of dependences and be able to easily compile the code for
   the TEST_BUILD. 
 */
#define __CVERSION__
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "general.h"
#include "functionals.h"

#if !defined __inline__
/* inline some stuff whenever possible */
#define __inline__
#endif

typedef void (*DaltonEnFunc)(real* res, const real* rho, 
                             const real* rho13, const real*grad);
void edrc_(real* drc, const real* rho, const real* rho13, const real* rhogrd);
void evwn_(real* vwn, const real* rho, const real* rho13, const real* rhogrd);
void ebck_(real* bck, const real* rho, const real* rho13, const real* rhogrd);
void elyp_(real* lyp, const real* rho, const real* rho13, const real* rhogrd);

void dftpot1(SecondDrv *ds, const real *w, const real* rho, const real* grad,
             const integer* triplet);
void condft_(void);

static __inline__
void test_var(real comp, real refer, const char* fun, const char* drv, 
              integer* counter)
{
      if(fabs(comp-refer)>2e-7+5e-5*(fabs(comp)+fabs(refer))) { 
        if(*counter<95)
            printf("%s %s: fin.diff: %12g found: %12g, diff=%g\n", 
                   fun, drv, refer, comp, fabs(comp-refer));
	++*counter;
	  }/*  else printf("Test '%s:%s' passed (expected: %g found: %g).\n",
               fun,drv, refer, comp); */
}

const integer GRID_STEP = 1;

/* test_first: test first order derivatives of given functional. 
   Note that really the restricted case is only tested...
*/
static integer
test_first(const char* fun_name, EnergyFunc func, FirstOrderFun first_func)
{ 
    integer i, j, k, failed = 0;
    real drho, dgra, resp, resm, num;
    FunFirstFuncDrv gga;
    for(i=1; i<=40; i+=GRID_STEP) {
	for(j=1; j<=40; j+=GRID_STEP) {
	    for(k=-19; k<=19; k+=GRID_STEP) {
		real rho   = i/40.0;
		real ngrad = j/40.0;
		real gracos= k/20.0;
		FunDensProp dt, dp = { 0.5*rho,0.2*rho, 0.2*ngrad,0.2*ngrad };
		dp.gradab = dp.grada*dp.gradb*gracos;
		/* TEST df1000 */
		drho = rho*1e-4;
		dt = dp; dt.rhoa -= drho; resm = func(&dt);
		dt = dp; dt.rhoa += drho; resp = func(&dt);
		drv1_clear(&gga);
		first_func(&gga, 1, &dp);
		num = (resp-resm)/(2*drho);
		test_var(gga.df1000,num, fun_name, "df1000", &failed);
		
		/* TEST df0010 */
		dgra = ngrad*1e-7;
		dt = dp; dt.grada -= dgra; resm = func(&dt);
		dt = dp; dt.grada += dgra; resp = func(&dt);
		num = (resp-resm)/(2*dgra);
		test_var(gga.df0010,num, fun_name, "df0010", &failed);
		
		/* TEST df00001 */
		if(fabs(gracos)<1e-5) continue;
		dgra = gracos*1e-7;
		dt = dp; dt.gradab -= dgra; resm = func(&dt);
		dt = dp; dt.gradab += dgra; resp = func(&dt);
		num = (resp-resm)/(2*dgra);
		test_var(gga.df00001,num, fun_name, "df00001", &failed);
	    }
	}
    }
    if(failed==0) printf("%-5s (first order derivatives): OK\n", fun_name);
    return failed;
}

/* test_second:
   test second order derivatives of given functional.
   It is assumed that the first order derivatives are OK.
*/
#define COMP_DER(c,delta,eps,field, fun, ord) \
   eps = c*delta; \
   drv##ord##_clear(&m); dt = dp; dt.field -= eps; fun ## _fun(&m, 1, &dt); \
   drv##ord##_clear(&p); dt = dp; dt.field += eps; fun ## _fun(&p, 1, &dt)

#define T2(der,derdif,eps,label) \
   num = (p.derdif-m.derdif)/(2*eps);\
   test_var(d.der, num, fname, label ":" #der, &fail)

static integer
test_second(const char* fname,
            FirstOrderFun first_fun, SecondOrderFun second_fun)
{ 
    integer i, j, k, fail = 0;
    real drho, dgra, num;
    FunFirstFuncDrv m, p;
    FunSecondFuncDrv d;
    for(i=1; i<=40; i+=GRID_STEP) {
	for(j=1; j<=40; j+=GRID_STEP) {
	    for(k=-19; k<=19; k+=GRID_STEP) {
		real rho   = i/40.0;
		real ngrad = j/40.0;
		real gracos= k/20.0;
		FunDensProp dt, dp = {0.5*rho, 0.2*rho, 0.5*ngrad, 0.3*ngrad};
		dp.gradab = dp.grada*dp.gradb*gracos;
		drv2_clear(&d);
		second_fun(&d, 1, &dp);
		drv1_clear(&m);
		first_fun(&m, 1, &dp);
		test_var(d.df1000, m.df1000, fname, "df1000X", &fail);
		test_var(d.df0010, m.df0010, fname, "df0010X", &fail);
		test_var(d.df00001,m.df00001,fname, "df00001X", &fail);
		test_var(d.df0100, m.df0100, fname, "df0100X", &fail);
		test_var(d.df0001, m.df0001, fname, "df0001X", &fail);
		test_var(d.df00001,m.df00001,fname, "df00001X", &fail);

		/* TEST df2000, df1010, df1001 and df10001  */
                COMP_DER(rho,1e-5,drho,rhoa, first,1);
                T2(df2000,  df1000,  drho, "A");
                T2(df1010,  df0010,  drho, "A");
                T2(df1001,  df0001,  drho, "A");
                T2(df10001, df00001, drho, "A");
		
		/* TEST df0200, df0101, df0110 and df01001 */
                COMP_DER(rho,1e-7,drho,rhob, first,1);
                T2(df0200,  df0100,  drho, "A");
                T2(df0110,  df0010,  drho, "A");
                T2(df0101,  df0001,  drho, "A");
                T2(df01001, df00001, drho, "A");
                T2(df1100,  df1000,  drho, "A");
		
		/* TEST df1010, df0110, df0020 */
                COMP_DER(ngrad,1e-5,dgra, grada, first,1);
                T2(df1010,  df1000,  dgra, "B");
                T2(df0020,  df0010,  dgra, "B");
                T2(df0110,  df0100,  dgra, "B");

               	/* TEST df1001, df0101, df0002 */
                COMP_DER(ngrad,1e-5,dgra, gradb, first,1);
                T2(df1001,  df1000,  dgra, "B");
                T2(df0101,  df0100,  dgra, "B");
                T2(df0002,  df0001,  dgra, "B");
                T2(df00011, df00001, dgra, "B");
            }
	}
    }
    if(fail==0) printf("%-5s (second order derivatives): OK\n", fname);
    return fail;
}

/* test_third:
   test third order derivatives of given functional.
   It is assumed that the second order derivatives are OK.
*/
static integer
test_third(const char* fname,
           SecondOrderFun second_fun, ThirdOrderFun third_fun)
{ 
    integer i, j, k, fail = 0;
    real eps, num;
    FunSecondFuncDrv m, p;
    FunThirdFuncDrv d;
    for(i=1; i<=40; i+=GRID_STEP) {
	for(j=1; j<=40; j+=GRID_STEP) {
	    for(k=-19; k<=19; k+=GRID_STEP) {
		real rho   = i/40.0;
		real ngrad = j/40.0;
		real gracos= k/20.0;
		FunDensProp dt, dp = {0.5*rho, 0.3*rho, 0.5*ngrad, 1*ngrad};
		dp.gradab = dp.grada*dp.gradb*gracos;
		drv3_clear(&d);
		third_fun(&d, 1, &dp);
		drv2_clear(&m);
		second_fun(&m, 1, &dp);
                test_var(d.df1000, m.df1000,  fname, "df1000X",  &fail);
                test_var(d.df0100, m.df0100,  fname, "df0100X",  &fail);
                test_var(d.df0010, m.df0010,  fname, "df0010X",  &fail);
                test_var(d.df0001, m.df0001,  fname, "df0001X",  &fail);
                test_var(d.df00001,m.df00001, fname, "df00001X", &fail);     
                test_var(d.df2000, m.df2000,  fname, "df2000X",  &fail);
                test_var(d.df1100, m.df1100,  fname, "df1100X",  &fail);
                test_var(d.df1010, m.df1010,  fname, "df1010X",  &fail);
                test_var(d.df1001, m.df1001,  fname, "df1001X",  &fail);
                test_var(d.df10001,m.df10001, fname, "df10001X", &fail);
                test_var(d.df0200, m.df0200,  fname, "df0200X",  &fail);
                test_var(d.df0110, m.df0110,  fname, "df0110X",  &fail);
                test_var(d.df0101, m.df0101,  fname, "df0101X",  &fail);
                test_var(d.df01001,m.df01001, fname, "df01001X", &fail);
                test_var(d.df0020, m.df0020,  fname, "df0020X",  &fail);
                test_var(d.df0011, m.df0011,  fname, "df0011X",  &fail);
                test_var(d.df00101,m.df00101, fname, "df00101X", &fail);
                test_var(d.df0002, m.df0002,  fname, "df0002X",  &fail);
                test_var(d.df00011,m.df00011, fname, "df00011X", &fail);
                test_var(d.df00002,m.df00002, fname, "df00002X", &fail);

		/* drhoa: test  */
                COMP_DER(rho,1e-7,eps,rhoa, second,2);
                T2(df3000,  df2000,  eps, "A");
                T2(df2100,  df1100,  eps, "A");
                T2(df2010,  df1010,  eps, "A");
                T2(df2001,  df1001,  eps, "A");
                T2(df1200,  df0200,  eps, "A");
                T2(df1110,  df0110,  eps, "A");
                T2(df1101,  df0101,  eps, "A");
                T2(df11001, df01001, eps, "A");
                T2(df1020,  df0020,  eps, "A");
                T2(df1011,  df0011,  eps, "A");
                T2(df1002,  df0002,  eps, "A");

                /*drhob: test */
                COMP_DER(rho,1e-7,eps,rhob, second,2);
		T2(df0300,  df0200,  eps, "A");
		T2(df0201,  df0101,  eps, "A");
                T2(df0210,  df0110,  eps, "A");
		T2(df0102,  df0002,  eps, "A");
               	T2(df0120,  df0020,  eps, "A");
		T2(df1200,  df1100,  eps, "A");
		T2(df0111,  df0011,  eps, "A");
               	T2(df1101,  df1001,  eps, "A");
              	T2(df1110,  df1010,  eps, "A");
		T2(df02001, df01001, eps, "A");
		T2(df11001, df10001, eps, "A");

		/* dgrada: test */
                COMP_DER(ngrad,1e-7,eps,grada, second,2);
                T2(df2010,  df2000,  eps, "A");
                T2(df0030,  df0020,  eps, "A");
	        T2(df0021,  df0011,  eps, "A");

               	/* dgradb: test */  
                COMP_DER(ngrad,1e-7,eps,gradb, second,2);
                T2(df2001,  df2000,  eps, "B");
                T2(df0003,  df0002,  eps, "B");
	        T2(df0012,  df0011,  eps, "B");
	    }
	}
    }
    if(fail==0) printf("%-5s (third order derivatives): OK\n", fname);
    return fail;
}

/* test_fourth:
   test fourth order derivatives of given functional.
   It is assumed that the second order derivatives are OK.
*/
static integer
test_fourth(const char* fname,
            ThirdOrderFun third_fun, FourthOrderFun fourth_fun)
{ 
    integer i, j, k, fail = 0;
    real eps, num;
    FunThirdFuncDrv m, p;
    FunFourthFuncDrv d;
    for(i=1; i<=40; i+=GRID_STEP) {
	for(j=1; j<=40; j+=GRID_STEP) {
	    for(k=-19; k<=19; k+=GRID_STEP) {
		real rho   = i/40.0;
		real ngrad = j/40.0;
		real gracos= k/20.0;
		FunDensProp dt, dp = {0.5*rho, 0.3*rho, 0.5*ngrad, 1*ngrad};
		dp.gradab = dp.grada*dp.gradb*gracos;
		drv4_clear(&d);	fourth_fun(&d, 1, &dp);
		drv3_clear(&m); third_fun(&m, 1, &dp);
                test_var(d.df1000,  m.df1000,  fname, "df1000X", &fail);
                test_var(d.df0100,  m.df0100,  fname, "df0100X", &fail);
                test_var(d.df0010,  m.df0010,  fname, "df0010X", &fail);
                test_var(d.df0001,  m.df0001,  fname, "df0001X", &fail);
                test_var(d.df00001, m.df00001, fname, "df00001X", &fail);
                test_var(d.df2000,  m.df2000,  fname, "df2000X", &fail);
                test_var(d.df1100,  m.df1100,  fname, "df1100X", &fail);
                test_var(d.df1010,  m.df1010,  fname, "df1010X", &fail);
                test_var(d.df1001,  m.df1001,  fname, "df1001X", &fail);
                test_var(d.df10001, m.df10001, fname, "df10001X", &fail);
                test_var(d.df0200,  m.df0200,  fname, "df0200X", &fail);
                test_var(d.df0110,  m.df0110,  fname, "df0110X", &fail);
                test_var(d.df0101,  m.df0101,  fname, "df0101X", &fail);
                test_var(d.df01001, m.df01001, fname, "df01001X", &fail);
                test_var(d.df0020,  m.df0020,  fname, "df0020X", &fail);
                test_var(d.df0011,  m.df0011,  fname, "df0011X", &fail);
                test_var(d.df00101, m.df00101, fname, "df00101X", &fail);
                test_var(d.df0002,  m.df0002,  fname, "df0002X", &fail);
                test_var(d.df00011, m.df00011, fname, "df00011X", &fail);
                test_var(d.df00002, m.df00002, fname, "df00002X", &fail);
                test_var(d.df3000,  m.df3000,  fname, "df3000X", &fail);
                test_var(d.df2100,  m.df2100,  fname, "df2100X", &fail);
                test_var(d.df2010,  m.df2010,  fname, "df2010X", &fail);
                test_var(d.df2001,  m.df2001,  fname, "df2001X", &fail);
                test_var(d.df20001, m.df20001, fname, "df20001X", &fail);
                test_var(d.df1200,  m.df1200,  fname, "df1200X", &fail);
                test_var(d.df1110,  m.df1110,  fname, "df1110X", &fail);
                test_var(d.df1101,  m.df1101,  fname, "df1101X", &fail);
                test_var(d.df11001, m.df11001, fname, "df11001X", &fail);
                test_var(d.df1020,  m.df1020,  fname, "df1020X", &fail);
                test_var(d.df1011,  m.df1011,  fname, "df1011X", &fail);
                test_var(d.df10101, m.df10101, fname, "df10101X", &fail);
                test_var(d.df1002,  m.df1002,  fname, "df1002X", &fail);
                test_var(d.df10011, m.df10011, fname, "df10011X", &fail);
                test_var(d.df10002, m.df10002, fname, "df10002X", &fail);
                test_var(d.df0300,  m.df0300,  fname, "df0300X", &fail);
                test_var(d.df0210,  m.df0210,  fname, "df0210X", &fail);
                test_var(d.df0201,  m.df0201,  fname, "df0201X", &fail);
                test_var(d.df02001, m.df02001, fname, "df02001X", &fail);
                test_var(d.df0120,  m.df0120,  fname, "df0120X", &fail);
                test_var(d.df0111,  m.df0111,  fname, "df0111X", &fail);
                test_var(d.df01101, m.df01101, fname, "df01101X", &fail);
                test_var(d.df0102,  m.df0102,  fname, "df0102X", &fail);
                test_var(d.df01011, m.df01011, fname, "df01011X", &fail);
                test_var(d.df01002, m.df01002, fname, "df01002X", &fail);
                test_var(d.df0030,  m.df0030,  fname, "df0030X", &fail);
                test_var(d.df0021,  m.df0021,  fname, "df0021X", &fail);
                test_var(d.df00201, m.df00201, fname, "df00201X", &fail);
                test_var(d.df0012,  m.df0012,  fname, "df0012X", &fail);
                test_var(d.df00111, m.df00111, fname, "df00111X", &fail);
                test_var(d.df00102, m.df00102, fname, "df00102X", &fail);
                test_var(d.df0003,  m.df0003,  fname, "df0003X", &fail);
                test_var(d.df00021, m.df00021, fname, "df00021X", &fail);
                test_var(d.df00012, m.df00012, fname, "df00012X", &fail);
                test_var(d.df00003, m.df00003, fname, "df00003X", &fail);

		/* drhoa: test  */
                COMP_DER(rho,1e-6,eps,rhoa, third,3);
                T2(df4000,  df3000, eps, "A");
                T2(df3100,  df2100, eps, "A");
                T2(df3010,  df2010, eps, "A");
                T2(df3001,  df2001, eps, "A");
                T2(df30001, df20001,eps, "A");
                T2(df2200,  df1200, eps, "A");
                T2(df2110,  df1110, eps, "A");
                T2(df2101,  df1101, eps, "A");
                T2(df21001, df11001,eps, "A");
                T2(df2020,  df1020, eps, "A");
                T2(df2011,  df1011, eps, "A");
                T2(df20101, df10101,eps, "A");
                T2(df2002,  df1002, eps, "A");
                T2(df20011, df10011,eps, "A");
                T2(df20002, df10002,eps, "A");
                T2(df1300,  df0300, eps, "A");
                T2(df1210,  df0210, eps, "A");
                T2(df1201,  df0201, eps, "A");
                T2(df12001, df02001,eps, "A");
                T2(df1120,  df0120, eps, "A");
                T2(df1111,  df0111, eps, "A");
                T2(df11101, df01101,eps, "A");
                T2(df1102,  df0102, eps, "A");
                T2(df11011, df01011,eps, "A");
                T2(df11002, df01002,eps, "A");
                T2(df1030,  df0030, eps, "A");
                T2(df1021,  df0021, eps, "A");
                T2(df10201, df00201,eps, "A");
                T2(df1012,  df0012, eps, "A");
                T2(df10111, df00111,eps, "A");
                T2(df10102, df00102,eps, "A");
                T2(df1003,  df0003, eps, "A");
                T2(df10021, df00021,eps, "A");
                T2(df10012, df00012,eps, "A");
                T2(df10003, df00003,eps, "A");

                /*drhob: test */
                COMP_DER(rho,1e-7,eps,rhob, third,3);
                T2(df3100,  df3000, eps, "B");
                T2(df2200,  df2100, eps, "B");
                T2(df2110,  df2010, eps, "B");
                T2(df2101,  df2001, eps, "B");
                T2(df21001, df20001,eps, "B");
                T2(df1300,  df1200, eps, "B");
                T2(df1210,  df1110, eps, "B");
                T2(df1201,  df1101, eps, "B");
                T2(df12001, df11001,eps, "B");
                T2(df1120,  df1020, eps, "B");
                T2(df1111,  df1011, eps, "B");
                T2(df11101, df10101,eps, "B");
                T2(df1102,  df1002, eps, "B");
                T2(df11011, df10011,eps, "B");
                T2(df11002, df10002,eps, "B");
                T2(df0400,  df0300, eps, "B");
                T2(df0310,  df0210, eps, "B");
                T2(df0301,  df0201, eps, "B");
                T2(df03001, df02001,eps, "B");
                T2(df0220,  df0120, eps, "B");
                T2(df0211,  df0111, eps, "B");
                T2(df02101, df01101,eps, "B");
                T2(df0202,  df0102, eps, "B");
                T2(df02011, df01011,eps, "B");
                T2(df02002, df01002,eps, "B");
                T2(df0130,  df0030, eps, "B");
                T2(df0121,  df0021, eps, "B");
                T2(df01201, df00201,eps, "B");
                T2(df0112,  df0012, eps, "B");
                T2(df01111, df00111,eps, "B");
                T2(df01102, df00102,eps, "B");
                T2(df0103,  df0003, eps, "B");
                T2(df01021, df00021,eps, "B");
                T2(df01012, df00012,eps, "B");
                T2(df01003, df00003,eps, "B");
		/* dgrada: test */
                COMP_DER(ngrad,1e-7,eps,grada, third,3);
                T2(df3010,  df3000, eps, "C");
                T2(df2110,  df2100, eps, "C");
                T2(df2020,  df2010, eps, "C");
                T2(df2011,  df2001, eps, "C");
                T2(df20101, df20001,eps, "C");
                T2(df20011, df20001,eps, "C");
                T2(df1210,  df1200, eps, "C");
                T2(df1120,  df1110, eps, "C");
                T2(df1111,  df1101, eps, "C");
                T2(df11101, df11001,eps, "C");
                T2(df1030,  df1020, eps, "C");
                T2(df1021,  df1011, eps, "C");
                T2(df10201, df10101,eps, "C");
                T2(df1012,  df1002, eps, "C");
                T2(df10111, df10011,eps, "C");
                T2(df10102, df10002,eps, "C");
                T2(df1003,  df1002, eps, "C");
                T2(df0310,  df0300, eps, "C");
                T2(df0220,  df0210, eps, "C");
                T2(df0211,  df0201, eps, "C");
                T2(df02101, df02001,eps, "C");
                T2(df0130,  df0120, eps, "C");
                T2(df0121,  df0111, eps, "C");
                T2(df01201, df01101,eps, "C");
                T2(df0112,  df0102, eps, "C");
                T2(df01111, df01011,eps, "C");
                T2(df01102, df01002,eps, "C");
                T2(df0040,  df0030, eps, "C");
                T2(df0031,  df0021, eps, "C");
                T2(df00301, df00201,eps, "C");
                T2(df0022,  df0012, eps, "C");
                T2(df00211, df00111,eps, "C");
                T2(df00202, df00102,eps, "C");
                T2(df0013,  df0003, eps, "C");
                T2(df00121, df00021,eps, "C");
                T2(df00112, df00012,eps, "C");
                T2(df00103, df00003,eps, "C");

               	/* dgradb: test */  
                COMP_DER(ngrad,1e-7,eps,gradb, third,3);
                T2(df3001,  df3000, eps, "C");
                T2(df2101,  df2100, eps, "C");
                T2(df2011,  df2010, eps, "C");
                T2(df2002,  df2001, eps, "C");
                T2(df20011, df20001,eps, "C");
                T2(df1201,  df1200, eps, "C");
                T2(df1111,  df1110, eps, "C");
                T2(df1102,  df1101, eps, "C");
                T2(df11011, df11001,eps, "C");
                T2(df1021,  df1020, eps, "C");
                T2(df1012,  df1011, eps, "C");
                T2(df10111, df10101,eps, "C");
                T2(df1003,  df1002, eps, "C");
                T2(df10021, df10011,eps, "C");
                T2(df10012, df10002,eps, "C");
                T2(df0301,  df0300, eps, "C");
                T2(df0211,  df0210, eps, "C");
                T2(df0202,  df0201, eps, "C");
                T2(df02011, df02001,eps, "C");
                T2(df0121,  df0120, eps, "C");
                T2(df0112,  df0111, eps, "C");
                T2(df01111, df01101,eps, "C");
                T2(df0103,  df0102, eps, "C");
                T2(df01021, df01011,eps, "C");
                T2(df01012, df01002,eps, "C");
                T2(df0031,  df0030, eps, "C");
                T2(df0022,  df0021, eps, "C");
                T2(df00211, df00201,eps, "C");
                T2(df0013,  df0012, eps, "C");
                T2(df00121, df00101,eps, "C");
                T2(df00112, df00102,eps, "C");
                T2(df0004,  df0003, eps, "C");
                T2(df00031, df00021,eps, "C");
                T2(df00022, df00012,eps, "C");
                T2(df00013, df00003,eps, "C");
            }
	}
    }
    if(fail==0) printf("%-5s (fourth order derivatives): OK\n", fname);
    return fail;
}

static integer
test_derivatives(Functional* f, integer *orders, DaltonEnFunc dal_fun)
{
    integer res = 0;
    /* if(dal_fun) res = test_energy(f->name, f->func, dal_fun); */
    if(!res && (!orders || orders[0]) )
       res = test_first(f->name,  f->func,   f->first);
    if(!res && (!orders || orders[1]) )
        res = test_second(f->name, f->first,  f->second);
    if(!res && (!orders || orders[2]) )
        res = test_third(f->name,  f->second, f->third);
    if(!res && (!orders || orders[3]) )
        res = test_fourth(f->name,  f->third, f->fourth);
    return res;
}


/* main:
   this is the main test program.
*/
integer
main(integer argc, char* argv[])
{
    integer res = 0, i, length, argidx, funcsel = 0;
    static integer funco[] = { 0, 0, 0, 0 };
    char* arg;
    Functional* func;
    
    if(argc<=1) {
	fprintf(stderr,
                "Functional derivative tester:\n"
                "usage: fun-tester [-n] <functional> <options>\n"
                "-n - test only selected order of derivatives\n"
                "example: fun-tester GGAKey becke=1 lyp=1\n");
        return 1;
    } 
    for(argidx=1; argidx<argc && argv[argidx][0] == '-'; argidx++)
        switch(argv[argidx][1]) {
        case '1': funcsel = 1; funco[0]=1; break;
        case '2': funcsel = 1; funco[1]=1; break;
        case '3': funcsel = 1; funco[2]=1; break;
        case '4': funcsel = 1; funco[3]=1; break;
        default: fprintf(stderr, "option %s is unknown.\n", argv[argidx]);
        }

    for(i=0; available_functionals[i]; i++)
        if(strcasecmp(argv[argidx], available_functionals[i]->name)==0)
            break;
    if(available_functionals[i]==NULL) {
        fprintf(stderr, "Functional '%s' not found.\n\n"
                "Available functionals:\n", argv[1]);
        for(i=0; available_functionals[i]; i++)
            fprintf(stderr, "    %s\n", available_functionals[i]->name);
        return 2;
    } 
    func = available_functionals[i];
    argidx++;
    for(length=1, i=argidx; i<argc; i++)
        length += strlen(argv[i])+1;

    arg = malloc(length);
    if(argc>argidx+1)
        strcpy(arg, argv[argidx]);
    else
        *arg = '\0';

    for(i=argidx; i<argc; i++) {
        strcat(arg, " ");
        strcat(arg, argv[i]);
    }

    if(!func->read(arg)) {
        fprintf(stderr, "Reading configuration for %s from '%s' failed.\n",
                argv[1], arg);
        return 3;
    }
    free(arg);

    if(func->report)
        func->report();
    res += test_derivatives(func, funcsel ? funco: NULL, NULL);

    if(res>0) 
	printf("%i tests failed.\n", res);
    else printf("OK\n");
    return res;
}
