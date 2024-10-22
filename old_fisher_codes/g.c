#include <stdlib.h>
#include <math.h>
#include <stdio.h>

int main(void) {
  double x, cc, sum;

  for(x=-0.1; x>-1.01; x-=.1) {
    sum = 0.;
    for(cc=0.55; cc<4000; cc+=1) sum += -(1. + exp(cc*x)*(cc*x-1) )/cc/cc;
    sum/=3.;

    printf("%9.6lf %9.6lf %9.6lf %9.6lf\n",
      x,
      x/3. - x*x/12. + x*x*x/108. - x*x*x*x*x/10800. + x*x*x*x*x*x*x/635040.,
      x/3. +0.00833333*x*x -0.00449074*x*x*x -0.000171875*x*x*x*x +0.0000775637*x*x*x*x*x,
      sum);
  }
  return(0);
}
