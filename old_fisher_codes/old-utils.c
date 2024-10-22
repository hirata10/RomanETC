#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "utils.h"

/* Resets cosmology to default values */
void set_default_cosmology(COSMOPARAM *p) {
  int i;
  p->T_cmb = 2.728;
  p->obh2 = 0.022;
  p->omh2 = 0.147;
  p->delta_zeta = 2e-5;
  p->ns = 0.95;
  p->alphas = 0.;
  p->h = 0.7;
  p->Omega_K = 0.;
  for(i=0; i<NWPARAM; i++) p->w[i] = 0.;
}

double get_H(double z, COSMOPARAM *p) {
  int i, imax;
  double a;
  double Omega_L, Omega_M, E, w1, H;

  E=1.;
  a = 1./(1.+z);
  Omega_M = p->omh2/(p->h*p->h);
  Omega_L = 1. - p->Omega_K - Omega_M;
  imax = (int)floor((1.-a)/DELTA_A);
  if (imax<0) imax=0;
  if (imax>=NWPARAM) imax=NWPARAM;
  for(i=0; i<imax; i++)
    E *= pow(1.-DELTA_A/(1.-DELTA_A*i), -3.*p->w[i]);
  if (imax<NWPARAM)
    E *= pow(a/(1.-DELTA_A*imax), -3.*p->w[imax]);

  H = p->h * sqrt(Omega_L*E + (Omega_M/a + p->Omega_K)/a/a) / HL_HIMPC;
  return(H);
}

double get_DAC(double z, COSMOPARAM *p) {
  double a, a_current, z_current, K_chi2, r;
  int i, imax, j, jj;
  double chi = 0.;
  double astep = DELTA_A/8.; /* Want an integer fraction of DELTA_A so  */
                             /* integrand is a polynomial in each slice */
  double x[]= {0.148874338981631, 0.433395934129247, 0.679409568299024, 0.865063366688985, 0.973906528517172};
  double w[]= {0.295524224714753, 0.269266719309996, 0.219086362515982, 0.149451349150581, 0.066671344308688};

  a = 1./(1.+z);
  imax = (int)floor((1.-a)/astep);

  for(i=0; i<imax; i++) {
    for(j=0; j<10; j++) {
      jj = j%5;
      a_current = 1. - (i+0.5)*astep + 0.5*astep*(j>=5? x[jj]: -x[jj]);
      z_current = 1./a_current - 1.;
      chi += w[jj]/get_H(z_current,p)/a_current/a_current;
    }
  }
  chi *= astep/2.;
  astep = 1.-imax*astep-a;
  for(j=0; j<10; j++) {
    jj = j%5;
    a_current = a + 0.5*astep + 0.5*astep*(j>=5? x[jj]: -x[jj]);
    z_current = 1./a_current - 1.;
    chi += w[jj]/get_H(z_current,p)/a_current/a_current*astep/2.;
  }

  r=0.;
  K_chi2 = -p->Omega_K * p->h * p->h * chi * chi / (HL_HIMPC*HL_HIMPC);
  if (fabs(K_chi2)<1e-8) {
    r = chi * (1. - K_chi2/6.*(1.-K_chi2/20.));
  } else if (K_chi2>0) {
    r = chi * sin(sqrt(K_chi2))/sqrt(K_chi2);
  } else if (K_chi2<0) {
    r = chi * sinh(sqrt(-K_chi2))/sqrt(-K_chi2);
  }

  return(r);
}

int main(void) {
  COSMOPARAM p;
  double z;

  set_default_cosmology(&p);
  p.w[3] = 0.5;

  for(z=0.05; z<2.1; z+=0.1) {
    if (z>2) z=1100;
    printf("%12.6lf %12.6lf %9.3lf\n", z, HL_HIMPC*get_H(z,&p), get_DAC(z,&p));
  }

  return(0);
}
