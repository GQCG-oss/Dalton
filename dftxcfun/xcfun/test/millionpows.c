#include <stdio.h>
#include <math.h>

int main(void)
{
  double res = 0, x=0.1;
  int i;
  for (i = 0; i<56e6;i++)
    {
      res += pow(x,4.1/3.0);
      x+=1e-5;
    }
  return (int)res;
}
