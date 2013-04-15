#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "xcfun.h"


int main(void)
{
  double d[] = {0.39E+02,
		0.38E+02,
		0.81E+06,
		0.82E+06,
		0.82E+06};
  double out[256];
  xc_functional fun;
  int res, i;
  clock_t t0,t1;
  printf("%s",xcfun_splash());
  printf("XCFun version: %g\n",xcfun_version());
  fun = xc_new_functional();
  xc_set(fun,XC_LYPC,1.0);
  //xc_set(fun,XC_PBEX,1.0);
  if ((res = xc_eval_setup(fun,
			   XC_A_B_GAA_GAB_GBB,
			   XC_PARTIAL_DERIVATIVES,
			   2)) == 0)
    {
      t0 = clock();
      for (i=0;i<1000000;i++)
	xc_eval(fun,d,out);
      t1 = clock();
      printf("time: %g\n",(t1-t0)/(double)CLOCKS_PER_SEC);
      for (i=0;i<xc_output_length(fun);i++)
	printf("E = %.16le\n",out[i]);
    }
  else
    {
      printf("Could not set up, error %i\n",res);
    }
  return EXIT_SUCCESS;
}
