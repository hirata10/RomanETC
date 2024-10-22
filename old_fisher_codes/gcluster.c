#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "nr_utils.h"
#include "utils.h"
#include "ps.h"
#include "fom.h"

#include "mnemonics.c"

int main(void) {
  COSMOPARAM p;
  double b, Pk, k, n, V, nmode, dlnk=0.02;
  double Flnb, nmodetot, d2tot;
  double L, CL;

  set_default_cosmology(&p);

  V = 6.0e7;
  n = 3.0e-4;
  b = 1.13;
b=1;
  Flnb = 0.;
  nmodetot = 0.;
  d2tot = 0.;

#if 1
  for(L=2; L<=1000.01; L+=1) {get_wl_ps(&p,1,5*L,1,1,1,&CL); printf("%4d %12.5le\n", (int)L, 25*CL);}
  return(0);
#endif

  for(k=0.001; k<0.2*p.h; k*=exp(dlnk)) {
    nmode = k*k*k*dlnk/(2.*M_PI*M_PI)*V;

    get_Delta2kNL(&p, 0.07, k, dlnk, 1, &Pk);
    Pk /= k*k*k/(2.*M_PI*M_PI);
    Pk *= b*b;

    Flnb += nmode * 2. / pow(1. + 1./(n*Pk),2.);
    nmodetot += nmode;
    d2tot += k*k*k/(2.*M_PI*M_PI) * Pk / (b*b) * dlnk;

    printf("k=%12.5le Pk=%12.5le nmode=%12.5le (tot=%12.5le) d2tot=%12.5le sigma(lnb)=%12.5le\n", k, Pk, nmode, nmodetot, d2tot, 1./sqrt(Flnb));
    
  }

  return(0);
}
