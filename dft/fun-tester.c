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
             const int* triplet);
void condft_(void);

static __inline__
void test_var(real comp, real refer, const char* fun, const char* drv, 
              int* counter)
{
      if(fabs(comp-refer)>2e-8+5e-5*(fabs(comp)+fabs(refer))) { 
        if(*counter<95)
            printf("%s %s: fin.diff: %12g found: %12g, diff=%g\n", 
                   fun, drv, refer, comp, fabs(comp-refer));
	++*counter;
	  }/*  else printf("Test '%s:%s' passed (expected: %g found: %g).\n",
               fun,drv, refer, comp); */
}

const int GRID_STEP = 1;

/* test_first: test first order derivatives of given functional. 
   Note that really the restricted case is only tested...
*/
static int
test_first(const char* fun_name, EnergyFunc func, FirstOrderFun first_func)
{ 
    int i, j, k, failed = 0;
    real drho, dgra, resp, resm, num;
    FirstFuncDrv gga;
    for(i=1; i<=40; i+=GRID_STEP) {
	for(j=1; j<=40; j+=GRID_STEP) {
	    for(k=-19; k<=19; k+=GRID_STEP) {
		real rho   = i/40.0;
		real ngrad = j/40.0;
		real gracos= k/20.0;
		DftDensProp dt, dp = { 0.5*rho,0.5*rho, 0.5*ngrad,0.5*ngrad };
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
static int
test_second(const char* fname,
            FirstOrderFun first_fun, SecondOrderFun second_fun)
{ 
    int i, j, k, fail = 0;
    real drho, dgra, num;
    FirstFuncDrv m, p;
    SecondFuncDrv d;
    for(i=1; i<=40; i+=GRID_STEP) {
	for(j=1; j<=40; j+=GRID_STEP) {
	    for(k=-19; k<=19; k+=GRID_STEP) {
		real rho   = i/40.0;
		real ngrad = j/40.0;
		real gracos= k/20.0;
		DftDensProp dt, dp = { 0.5*rho, 0.2*rho, 0.5*ngrad, 0.2*ngrad };
		dp.gradab = dp.grada*dp.gradb*gracos;
		drv2_clear(&d);
		second_fun(&d, 1, &dp);
		drv1_clear(&m);
		first_fun(&m, 1, &dp);
		test_var(d.df1000, m.df1000, fname, "df1000X", &fail);
		test_var(d.df0010, m.df0010, fname, "df0010X", &fail);
		test_var(d.df00001,m.df00001,fname, "df00001X", &fail);

		drv2_clear(&d);
		second_fun(&d, 1, &dp);
		drv1_clear(&m);
		first_fun(&m, 1, &dp);
		test_var(d.df0100, m.df0100, fname, "df0100X", &fail);
		test_var(d.df0001, m.df0001, fname, "df0001X", &fail);
		test_var(d.df00001,m.df00001,fname, "df00001X", &fail);

		/* TEST df2000, df1010, df1001 and df10001  */
		drho = rho*1e-5;
		drv1_clear(&m);
		dt = dp; dt.rhoa -= drho; first_fun(&m, 1, &dt);
		drv1_clear(&p);
		dt = dp; dt.rhoa += drho; first_fun(&p, 1, &dt);
		num = (p.df1000-m.df1000)/(2*drho);
		test_var(d.df2000, num, fname, "df2000", &fail);
		num = (p.df0010-m.df0010)/(2*drho);
		test_var(d.df1010, num, fname, "df1010a", &fail);
		num = (p.df0001-m.df0001)/(2*drho);
		test_var(d.df1001, num, fname, "df1001a", &fail);
		num = (p.df00001-m.df00001)/(2*drho);
		test_var(d.df10001, num, fname, "df10001", &fail);
		
		/* TEST df0200, df0101, df0110 and df01001 */
	       	drho = rho*1e-7;
		drv1_clear(&m);
		dt = dp; dt.rhob -= drho; first_fun(&m, 1, &dt);
		drv1_clear(&p);
		dt = dp; dt.rhob += drho; first_fun(&p, 1, &dt);
		num = (p.df0100-m.df0100)/(2*drho);
		test_var(d.df0200, num, fname, "df0200", &fail);
		// printf("%g %g %g %g\n", dp.rhoa, dp.rhob, dp.grada, dp.gradb);
		num = (p.df0010-m.df0010)/(2*drho);
		test_var(d.df0110, num, fname, "df0110a", &fail);        
		num = (p.df0001-m.df0001)/(2*drho);
		test_var(d.df0101, num, fname, "df0101a", &fail);
		num = (p.df00001-m.df00001)/(2*drho);
		test_var(d.df01001, num, fname, "df01001", &fail);
		
		/* TEST df1010, df0110, df0020 */
		 dgra = ngrad*1e-5;
	         drv1_clear(&m);
		 dt = dp; dt.grada -= dgra; first_fun(&m, 1, &dt);
		 drv1_clear(&p);
		 dt = dp; dt.grada += dgra; first_fun(&p, 1, &dt);
		 num = (p.df1000-m.df1000)/(2*dgra);
		 test_var(d.df1010, num, fname, "df1010b", &fail);
		 num = (p.df0010-m.df0010)/(2*dgra);
		 test_var(d.df0020, num, fname, "df0020", &fail);
                 num = (p.df0100-m.df0100)/(2*dgra);
		 test_var(d.df0110, num, fname, "df0110b", &fail); 
               	/* TEST df1001, df0101, df0002 */
		 dgra = ngrad*1e-5;
	         drv1_clear(&m);
		 dt = dp; dt.gradb -= dgra; first_fun(&m, 1, &dt);
		 drv1_clear(&p);
		 dt = dp; dt.gradb += dgra; first_fun(&p, 1, &dt);
		 num = (p.df1000-m.df1000)/(2*dgra);
		 test_var(d.df1001, num, fname, "df1001b", &fail);
		 num = (p.df0001-m.df0001)/(2*dgra);
		 test_var(d.df0002, num, fname, "df0002", &fail);
                 num = (p.df0100-m.df0100)/(2*dgra);
		 test_var(d.df0101, num, fname, "df0101b", &fail);   
		
		/* TEST REMAINING     df1100, df0200, df1001 */
		 drho = rho*1e-7;
		 drv1_clear(&m);
		 dt = dp; dt.rhob -= drho; first_fun(&m, 1, &dt);
		 drv1_clear(&p);
		 dt = dp; dt.rhob += drho; first_fun(&p, 1, &dt);
		 num = (p.df0100-m.df0100)/(2*drho);
		 test_var(d.df0200, num, fname, "df0200", &fail);
		 num = (p.df1000-m.df1000)/(2*drho);
		 test_var(d.df1100, num, fname, "df1100", &fail);
		 dgra = ngrad*1e-7;
		 drv1_clear(&m);
		 dt = dp; dt.gradb -= dgra; first_fun(&m, 1, &dt);
		 drv1_clear(&p);
		 dt = dp; dt.gradb += dgra; first_fun(&p, 1, &dt);
		 num = (p.df1000-m.df1000)/(2*dgra);
		 test_var(d.df1001, num, fname, "df1001b", &fail);
		 num = (p.df0010-m.df0010)/(2*dgra);
		 test_var(d.df0011, num, fname, "df0011", &fail);
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
static int
test_third(const char* fun_name,
            SecondOrderFun second_fun, ThirdOrderFun third_fun)
{ 
    int i, j, k, failed = 0;
    real drho, dgra, num;
    SecondFuncDrv m, p;
    ThirdFuncDrv d;
    for(i=1; i<=40; i+=GRID_STEP) {
	for(j=1; j<=40; j+=GRID_STEP) {
	    for(k=-19; k<=19; k+=GRID_STEP) {
		real rho   = i/40.0;
		real ngrad = j/40.0;
		real gracos= k/20.0;
		DftDensProp dt, dp = {0.5*rho, 0.5*rho, 0.5*ngrad, 0.5*ngrad};
		dp.gradab = dp.grada*dp.gradb*gracos;
		drv3_clear(&d);
		third_fun(&d, 1, &dp);
		drv2_clear(&m);
		second_fun(&m, 1, &dp);
		test_var(d.df1000, m.df1000, fun_name, "df1000X",&failed);
               	test_var(d.df0100, m.df0100, fun_name, "df0100X",&failed);
		test_var(d.df0010, m.df0010, fun_name, "df0010X",&failed);
                test_var(d.df0001, m.df0001, fun_name, "df0001X",&failed);
		test_var(d.df2000, m.df2000, fun_name, "df2000X",&failed);
               	test_var(d.df0200, m.df0200, fun_name, "df0200X",&failed);
		test_var(d.df1100, m.df1100, fun_name, "df1100X",&failed);
                test_var(d.df0011, m.df0011, fun_name, "df0011X",&failed);
		test_var(d.df1010, m.df1010, fun_name, "df1010X",&failed);
		test_var(d.df1001, m.df1001, fun_name, "df1001X",&failed);
               	test_var(d.df0101, m.df0101, fun_name, "df0101X",&failed);
		test_var(d.df0110, m.df0110, fun_name, "df0110X",&failed); 
		test_var(d.df0020, m.df0020, fun_name, "df0020X",&failed);
                test_var(d.df0002, m.df0002, fun_name, "df0002X",&failed); 
                test_var(d.df00001, m.df00001, fun_name, "df00001X",&failed);
                test_var(d.df10001, m.df10001, fun_name, "df10001X",&failed);  
                test_var(d.df01001, m.df01001, fun_name, "df01001X",&failed);

		/* drhoa: test  */
		drho = rho*1e-7;
		drv2_clear(&m);
		dt = dp; dt.rhoa -= drho; second_fun(&m, 1, &dt);
		drv2_clear(&p);
		dt = dp; dt.rhoa += drho; second_fun(&p, 1, &dt);
		num = (p.df2000-m.df2000)/(2*drho);
		test_var(d.df3000, num, fun_name, "df3000", &failed); 
		num = (p.df1100-m.df1100)/(2*drho);
		test_var(d.df2100, num, fun_name, "df2100a", &failed);
		num = (p.df0200-m.df0200)/(2*drho);
		test_var(d.df1200, num, fun_name, "df1200a", &failed);
		num = (p.df1010-m.df1010)/(2*drho);
		test_var(d.df2010, num, fun_name, "df2010a", &failed);
               	num = (p.df1001-m.df1001)/(2*drho);
		test_var(d.df2001, num, fun_name, "df2001a", &failed);
		num = (p.df0020-m.df0020)/(2*drho);
		test_var(d.df1020, num, fun_name, "df1020a", &failed);
		num = (p.df0002-m.df0002)/(2*drho);
		test_var(d.df1002, num, fun_name, "df1002a", &failed);
		num = (p.df0011-m.df0011)/(2*drho);
		test_var(d.df1011, num, fun_name, "df1011a", &failed);
		num = (p.df0101-m.df0101)/(2*drho);
		test_var(d.df1101, num, fun_name, "df1101a", &failed);
                num = (p.df0110-m.df0110)/(2*drho);
		test_var(d.df1110, num, fun_name, "df1110a", &failed); 
                num = (p.df10001-m.df10001)/(2*drho);
		test_var(d.df20001, num, fun_name, "df20001a", &failed);   
		num = (p.df01001-m.df01001)/(2*drho);
		test_var(d.df11001, num, fun_name, "df11001a", &failed);  

                /*drhob: test */
               	drho = rho*1e-7;
		drv2_clear(&m);
		dt = dp; dt.rhob -= drho; second_fun(&m, 1, &dt);
		drv2_clear(&p);
		dt = dp; dt.rhob += drho; second_fun(&p, 1, &dt);
		num = (p.df0200-m.df0200)/(2*drho);
		test_var(d.df0300, num, fun_name, "df0300", &failed); 
		num = (p.df0101-m.df0101)/(2*drho);
		test_var(d.df0201, num, fun_name, "df0201a", &failed);
                num = (p.df0110-m.df0110)/(2*drho);
		test_var(d.df0210, num, fun_name, "df0210a", &failed);  
		num = (p.df0002-m.df0002)/(2*drho);
		test_var(d.df0102, num, fun_name, "df0102a", &failed);
               	num = (p.df0020-m.df0020)/(2*drho);
		test_var(d.df0120, num, fun_name, "df0120a", &failed);  
		num = (p.df1100-m.df1100)/(2*drho);
		test_var(d.df1200, num, fun_name, "df1200a", &failed);
		num = (p.df0011-m.df0011)/(2*drho);
		test_var(d.df0111, num, fun_name, "df0111a", &failed);
               	num = (p.df1001-m.df1001)/(2*drho);
		test_var(d.df1101, num, fun_name, "df1101b", &failed);
              	num = (p.df1010-m.df1010)/(2*drho);
		test_var(d.df1110, num, fun_name, "df1110b", &failed);    
		num = (p.df01001-m.df01001)/(2*drho);
		test_var(d.df02001, num, fun_name, "df02001a", &failed);   
		num = (p.df10001-m.df10001)/(2*drho);
		test_var(d.df11001, num, fun_name, "df11001b", &failed);
		/* dgrada: test */
		dgra = ngrad*1e-7;
		drv2_clear(&m);
		dt = dp; dt.grada -= dgra; second_fun(&m, 1, &dt);
		drv2_clear(&p);
		dt = dp; dt.grada += dgra; second_fun(&p, 1, &dt);
                num = (p.df0020-m.df0020)/(2*dgra);
	        test_var(d.df0030, num, fun_name, "df0030", &failed);
	        num = (p.df0011-m.df0011)/(2*dgra);
		test_var(d.df0021, num, fun_name, "df0021", &failed);

               	/* dgradb: test */  
               	dgra = ngrad*1e-7;
		drv2_clear(&m);
		dt = dp; dt.gradb -= dgra; second_fun(&m, 1, &dt);
		drv2_clear(&p);
		dt = dp; dt.gradb += dgra; second_fun(&p, 1, &dt);
                num = (p.df0002-m.df0002)/(2*dgra);
	        test_var(d.df0003, num, fun_name, "df0003", &failed);
	        num = (p.df0011-m.df0011)/(2*dgra);
		test_var(d.df0012, num, fun_name, "df0012", &failed);
	    }
	}
    }
    if(failed==0) printf("%-5s (third order derivatives): OK\n", fun_name);
    return failed;
}

static int
test_derivatives(Functional* f, int *orders, DaltonEnFunc dal_fun)
{
    int res = 0;
    /* if(dal_fun) res = test_energy(f->name, f->func, dal_fun); */
    if(!res && (!orders || orders[0]) )
       res = test_first(f->name,  f->func,   f->first);
    if(!res && (!orders || orders[1]) )
        res = test_second(f->name, f->first,  f->second);
    if(!res && (!orders || orders[2]) )
        res = test_third(f->name,  f->second, f->third);
    return res;
}


/* main:
   this is the main test program.
*/
int
main(int argc, char* argv[])
{
    int res = 0, i, length, argidx, funcsel = 0;
    static int funco[] = { 0, 0, 0 };
    char* arg;
    Functional* func;
    DftDensProp  dp;
    FirstFuncDrv ds;

    dp.rhoa = dp.rhob = 1.0000;
    dp.rhoa -= 0.0001;
    dp.grada = dp.gradb = sqrt(3);
    dp.gradab = dp.grada*dp.gradab; 
    
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

void FSYM(dftsethf)(real *w){}
real FSYM(dftgethf)(void){return 0;}
void FSYM(dftsetcam)(real *b, real *mu) {}
void fort_print(const char *fmt, ...)
{
    va_list va_args;

    va_start(va_args, fmt);
    vprintf(fmt, va_args);puts("");
    va_end(va_args);
}
