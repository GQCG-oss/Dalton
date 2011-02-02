#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int countnans1(int n, const double *d)
{
  int m = 0;
  int i;
  for (i=0;i<n;i++)
    if (isnan(d[i]))
      m++;
  return m;
}

#define FAST_ISNAN(x) (*(unsigned long long *)(&x) == 4644090825121202176L)

int countnans2(int n, const double *d)
{
  int m = 0, i;
  for (i=0;i<n;i++)
    if (FAST_ISNAN(d[i]))
      m++;
  return m;
}


union u { double x[8]; unsigned long long ull[8]; };

int main(void)
{
  int i,m = 0, n = 10000000;
  double *d = malloc(sizeof(*d)*n);  
  union u uu;
  for (i=0;i<n;i++)
    d[i] = i+i*i;
  d[17] = NAN;
  d[18] = d[17];
  for (i=0;i<20;i++)
    m += countnans2(n,d);
  uu.x[3] = NAN;
  uu.x[4] = NAN;
  printf("%llo\n",uu.ull[3]);
  printf("%llo\n",uu.ull[4]);
  for (i =0;i<8;i++)
    printf("%hho\n",((unsigned char *)(d+17))[i]);
  for (i =0;i<8;i++)
    printf("%hho\n",((unsigned char *)(d+18))[i]);

  return uu.ull[3] == 9221120237041090560L;
}
