#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "nr_utils.h"
#include "utils.h"
#include "ps.h"

void DNgrowth_derivative(COSMOPARAM *p, double lny, double *X, double *dXdy, int flag) {
  static double z, H, y, om;

  if (flag) {
    y = exp(lny);
    z = y-1;
    H = HL_HIMPC * get_H(z,p) / p->h;
    om = p->omh2*y*y*y/(H*H*p->h*p->h);
  }
  dXdy[0] = X[1]/H+X[0];
  dXdy[1] = 3*X[1]+1.5*om*H*X[0];
}

void DNgrowth_step(COSMOPARAM *p, double lny, double *Xin, double *Xout, double dlny, int flag) {
  double k1[2], k2[2], k3[2], k4[2], X[2];
  int i;
  double lnyh = lny+0.5*dlny;

  /* Get derivatives */
  DNgrowth_derivative(p, lny, Xin, k1, flag);

  for(i=0;i<2;i++) X[i] = Xin[i]+0.5*dlny*k1[i];
  DNgrowth_derivative(p, lnyh, X, k2, 1);

  for(i=0;i<2;i++) X[i] = Xin[i]+0.5*dlny*k2[i];
  DNgrowth_derivative(p, lnyh, X, k3, 0);

  for(i=0;i<2;i++) X[i] = Xin[i]+dlny*k3[i];
  DNgrowth_derivative(p, lny+dlny, X, k4, 1);

  for(i=0;i<2;i++) Xout[i] = Xin[i] + dlny/6.*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
}

/* Extrapolated growth function */
double DNget_growth_normhiz(COSMOPARAM *p, double z) {
  double Xold[2], Xnew[2];
  double lny, dlny=-0.025;
  double om_m, x, a__3, H;
  long N=240;
  long i;

  lny = log(1.+z) - N*dlny;
  Xold[0] = 1.;
  Xold[1] = -get_H(exp(lny)-1,p) * HL_HIMPC / p->h;

  for(i=0; i<N; i++) {
    DNgrowth_step(p, lny, Xold, Xnew, dlny, i? 0: 1);
    lny += dlny;
    Xold[0] = Xnew[0]; Xold[1]=Xnew[1];
  }

#ifdef GROWTH_Z1
  return(Xnew[0]/(1.+z) * (z<1? p->growthAmp: exp(p->growthSlope)) * pow(1.+z, -p->growthCurve));
#endif

  /* gamma parameterization.  x = ln Omega_m(z).
   * In Lambda CDM we find
   * Delta ln D = -1/3 * int xe^{0.55x} dx/(1-e^x).
   * Expand integrand using Bernoulli polynomials:
   * Delta ln D = Delta gamma * sum_n B_n(0.55) x^{n+1} / ( 3 * (n+1)! )
   * Here use first 5 terms.
   *
   * Variable identification: growthAmp = G_0, growthSlope = Delta gamma, growthCurve = Delta f.
   */
#ifdef GROWTH_GAMMA
  a__3 = pow(1.+z,3.);
  H = HL_HIMPC * get_H(z,p);
  om_m = p->omh2*a__3 / H / H;
  x = log(om_m);
  return(Xnew[0]/(1.+z) * p->growthAmp * pow(1.+z, -p->growthCurve)
    * exp(p->growthSlope*(
        x/3. +0.00833333*x*x -0.00449074*x*x*x -0.000171875*x*x*x*x +0.0000775637*x*x*x*x*x
    )));
#endif

  return(Xnew[0]/(1.+z) * p->growthAmp * pow(1.+z,-p->growthSlope +0.5*log(1.+z)*p->growthCurve));
}

/* Sound horizon s in Mpc -- from Eisenstein/Hu TF */
double getSoundHorizon(COSMOPARAM *p) {
  double th27,omh2,obh2;
  double b1,b2,zd,zeq,keq,rd,req,s;

  /* Numbers to re-compute */
  th27 = p->T_cmb / 2.7;
  omh2 = p->omh2;
  obh2 = p->obh2;

  /* redshift at decoupling, equality, sound horizon */
  b1 = 0.313 * pow(omh2, -0.419) * ( 1 + 0.607*pow(omh2, 0.674) );
  b2 = 0.238 * pow(omh2, 0.223);
  zd = 1291. * pow(omh2, 0.251) / ( 1 + 0.659*pow(omh2, 0.828) ) * ( 1 + b1*pow(obh2,b2) );
  zeq = 25000. * omh2 / pow(th27,4.);
  keq = 0.0746 * omh2 / th27 / th27; /* in Mpc^-1 */
  rd = 31500.*obh2 / pow(th27,4.) / zd;
  req = zd*rd/zeq;
  s = 1.632993161855/keq/sqrt(req) * log((sqrt(1+rd)+sqrt(rd+req))/(1+sqrt(req)));
  return(s);
}

/* Transfer fcn normalized to 1 at k=0 (Eisenstein/Hu)
 * Returns: T(k) in matter era
 * Note: k in Mpc^-1. (no h)
 *
 * If k<0, returns no-wiggles 
 *
 * Contains extremely crude hack for neutrino masses (C.M.H.)
 */
double DNeh_trf_mpc(COSMOPARAM *p, double k) {

  int nowflag = 0;
  double th27, kmpc, omh2, obh2;
  double q, alphac, betac, a1, a2, b1, b2, bm, zd, zeq, keq, s, req, rd;
  double f, C0, T0_k1b, T0_kab, Tc;
  double ksilk, y, alphab, betab, betanode__ks, betab__ks, Tb, stilde, T;
  double xsupp, xnumer, kmpc0;
  double alphaGamma, Gammaeffh, ks043, L0;

  /* No-wiggles flag */
  if (k<0) {
    nowflag=1;
    k=fabs(k);
  }

  /* Numbers to re-compute */
  th27 = p->T_cmb / 2.7;
  omh2 = p->omh2;
  obh2 = p->obh2;

  /* redshift at decoupling, equality, sound horizon */
  b1 = 0.313 * pow(omh2, -0.419) * ( 1 + 0.607*pow(omh2, 0.674) );
  b2 = 0.238 * pow(omh2, 0.223);
  zd = 1291. * pow(omh2, 0.251) / ( 1 + 0.659*pow(omh2, 0.828) ) * ( 1 + b1*pow(obh2,b2) );
  zeq = 25000. * omh2 / pow(th27,4.);
  keq = 0.0746 * omh2 / th27 / th27; /* in Mpc^-1 */
  rd = 31500.*obh2 / pow(th27,4.) / zd;
  req = zd*rd/zeq;
  s = 1.632993161855/keq/sqrt(req) * log((sqrt(1+rd)+sqrt(rd+req))/(1+sqrt(req)));

  /* EH parameters */
  a1 = pow(46.9*omh2, 0.670) * ( 1 + pow(32.1*omh2, -0.532) );
  a2 = pow(12.0*omh2, 0.424) * ( 1 + pow(45.0*omh2, -0.582) );
  b1 = 0.944 / ( 1 + pow(458*omh2, -0.708) );
  b2 = pow( 0.395*omh2, -0.0266);
  bm = obh2/omh2;
  alphac = pow(a1, -bm) * pow(a2, -bm*bm*bm);
  betac = 1./( 1 + b1*(pow(1-bm, b2)-1) );

  /* k-independent baryon parameters */
  ksilk = 1.6 * pow(obh2, 0.52) * pow(omh2, 0.73) * ( 1 + pow(10.4*omh2, -0.95) );
  y = (1.+zeq)/(1.+zd);
  alphab = 2.07*keq*s*pow(1+rd,-0.75) * y * ( -6.*sqrt(1+y) + (2.+3.*y)*log((sqrt(1.+y)+1.)/(sqrt(1.+y)-1.)) );
  betab = 0.5 + bm + (3.-2.*bm)*sqrt(295.84*omh2*omh2+1.);

  /* More parameters */
  xnumer = 8.41*pow(omh2,0.435)/s;
  kmpc0 = omh2/(th27*th27);

  /* wavenumber in Mpc^-1 comoving, no h */
  kmpc = k;

  /* If no-wiggles */
  if (nowflag) {

    alphaGamma = 1. - 0.318*log(431.*omh2)*bm + 0.38*log(22.3*omh2)*bm*bm;
    ks043 = 0.43*kmpc*s;
    Gammaeffh = omh2 * ( alphaGamma + (1.-alphaGamma)/(1+ks043*ks043*ks043*ks043) );
    q = kmpc/Gammaeffh * th27*th27;
    L0 = log(5.436563657+1.8*q);
    C0 = 14.2 + 731./(1. + 62.5*q);
    T = L0/(L0+C0*q*q);    
  }
  else {

    /* CDM piece */
    q = kmpc/kmpc0;
    f = 1./(1. + pow(kmpc*s/5.4, 4.));
    C0 = 386./(1+69.9*pow(q,1.08));
    xsupp = q*q/log( M_E + 1.8*betac*q );
    T0_kab = 1./(1. + (C0+14.2/alphac)*xsupp);
    T0_k1b = 1./(1. + (C0+14.2)*xsupp);
    Tc = f*T0_k1b + (1-f)*T0_kab;

    /* Baryonic piece */
    betanode__ks = xnumer/kmpc;
    betab__ks = betab/kmpc/s;
    stilde = s*pow(1.+betanode__ks*betanode__ks*betanode__ks, -0.33333333);
    Tb = 1./(1. + (C0+14.2)*q*q/log( M_E + 1.8*q ))/(1+kmpc*kmpc*s*s/27.04)
         + alphab/(1+betab__ks*betab__ks*betab__ks)*exp(-pow(kmpc/ksilk,1.4));
    Tb *= sin(kmpc*stilde)/(kmpc*stilde);

    T = bm*Tb + (1.-bm)*Tc;
  }

  /* Hack for neutrino masses; based on scaling of 0.1 eV result.
   * Intended for Fisher matrix use only.
   */
  T *= exp(-33.8*p->onuh2/(1.+1.e-4/k/k));

  return(T);
}

/* BAO amplitude parameter */
double get_A0_BAO(COSMOPARAM *p) {

  double th27, kmpc, omh2, obh2;
  double q, alphac, betac, a1, a2, b1, b2, bm, zd, zeq, keq, s, req, rd;
  double f, C0, T0_k1b, T0_kab, Tc;
  double ksilk, y, alphab, betab, betanode__ks, betab__ks, Tb, stilde, T;
  double xsupp, xnumer, kmpc0;
  double alphaGamma, Gammaeffh, ks043, L0;

  /* Numbers to re-compute */
  th27 = p->T_cmb / 2.7;
  omh2 = p->omh2;
  obh2 = p->obh2;

  /* redshift at decoupling, equality, sound horizon */
  b1 = 0.313 * pow(omh2, -0.419) * ( 1 + 0.607*pow(omh2, 0.674) );
  b2 = 0.238 * pow(omh2, 0.223);
  zd = 1291. * pow(omh2, 0.251) / ( 1 + 0.659*pow(omh2, 0.828) ) * ( 1 + b1*pow(obh2,b2) );
  zeq = 25000. * omh2 / pow(th27,4.);
  keq = 0.0746 * omh2 / th27 / th27; /* in Mpc^-1 */
  rd = 31500.*obh2 / pow(th27,4.) / zd;
  req = zd*rd/zeq;
  s = 1.632993161855/keq/sqrt(req) * log((sqrt(1+rd)+sqrt(rd+req))/(1+sqrt(req)));

  /* EH parameters */
  a1 = pow(46.9*omh2, 0.670) * ( 1 + pow(32.1*omh2, -0.532) );
  a2 = pow(12.0*omh2, 0.424) * ( 1 + pow(45.0*omh2, -0.582) );
  b1 = 0.944 / ( 1 + pow(458*omh2, -0.708) );
  b2 = pow( 0.395*omh2, -0.0266);
  bm = obh2/omh2;
  alphac = pow(a1, -bm) * pow(a2, -bm*bm*bm);
  betac = 1./( 1 + b1*(pow(1-bm, b2)-1) );

  /* k-independent baryon parameters */
  ksilk = 1.6 * pow(obh2, 0.52) * pow(omh2, 0.73) * ( 1 + pow(10.4*omh2, -0.95) );
  y = (1.+zeq)/(1.+zd);
  alphab = 2.07*keq*s*pow(1+rd,-0.75) * y * ( -6.*sqrt(1+y) + (2.+3.*y)*log((sqrt(1.+y)+1.)/(sqrt(1.+y)-1.)) );
  betab = 0.5 + bm + (3.-2.*bm)*sqrt(295.84*omh2*omh2+1.);

  /* More parameters */
  xnumer = 8.41*pow(omh2,0.435)/s;
  kmpc0 = omh2/(th27*th27);

  /* wavenumber in Mpc^-1 comoving, no h */
  kmpc = 0.2*p->h;

  /* CDM piece */
  q = kmpc/kmpc0;
  f = 1./(1. + pow(kmpc*s/5.4, 4.));
  C0 = 386./(1+69.9*pow(q,1.08));
  xsupp = q*q/log( M_E + 1.8*betac*q );
  T0_kab = 1./(1. + (C0+14.2/alphac)*xsupp);
  T0_k1b = 1./(1. + (C0+14.2)*xsupp);
  Tc = f*T0_k1b + (1-f)*T0_kab;

  /* Baryonic piece -- part with no Silk damping or oscillation */
  betanode__ks = xnumer/kmpc;
  betab__ks = betab/kmpc/s;
  stilde = s*pow(1.+betanode__ks*betanode__ks*betanode__ks, -0.33333333);
  Tb = 1./(1. + (C0+14.2)*q*q/log( M_E + 1.8*q ))/(1+kmpc*kmpc*s*s/27.04)
       + alphab/(1+betab__ks*betab__ks*betab__ks);

  return(bm/(1.-bm) * Tb/Tc / sqrt(2.) / M_PI);
}

/* Gets matter density fraction at that redshift */
double get_Om(COSMOPARAM *p, double z) {
  double H,om;
  H = HL_HIMPC * get_H(z,p);
  om = p->omh2*(1.+z)*(1.+z)*(1.+z)/(H*H);
  return(om);
}

/* Linear power spectrum at specified z, range of k: D2[0..nk-1]
 * Returns no-wiggles PS if k<0
 */
void get_Delta2k(COSMOPARAM *p, double z, double k0, double dlnk, long nk, double *D2) {
  double G;
  long ik;
  double k, T, zeta2delta, val_4pga2r;
  double kpivot = 0.05;

  /* Poisson eqn coefficient */
  val_4pga2r = 1.5*p->omh2*(1.+z)/HL_HIMPC/HL_HIMPC;

  /* Potential growth function */
  G = (1.+z)*DNget_growth_normhiz(p,z);

  /* Consider each wavenumber */
  for(ik=0; ik<nk; ik++) {
    k = fabs(k0)*exp(ik*dlnk);
    T = DNeh_trf_mpc(p,k0<0? -k: k);
    zeta2delta = -0.6 * G * T / val_4pga2r * k*k;
    D2[ik] = zeta2delta * zeta2delta * exp(2.*p->lndelta_zeta)
             * pow(k/kpivot,p->ns-1.+0.5*log(k/kpivot)*p->alphas);
  }
}

/* Gets linear fluctuation amplitude in ball of size R Mpc
 * set z=0, R=8/h to get sigma_8
 */
double get_linsigma(COSMOPARAM *p, double z, double R) {
  double k0, dlnk;
  long ik;
  double D2[750];
  double T, var, kR;

  k0 = (1e-5)/R;
  dlnk = 0.01*log(10.);
  get_Delta2k(p,z,k0,dlnk,750,D2);

  /* Integrate linear fluctuations */
  var = 0.;
  for(ik=0; ik<750; ik++) {
    kR = k0*exp(ik*dlnk)*R;
    T = kR>1e-3? 3./(kR*kR)*(-cos(kR)+sin(kR)/kR): 1.; /* tophat transfer fcn */
    var += T*T*D2[ik];
  }
  var *= dlnk;
  return(sqrt(var));
}

/* Mass based formula */
double get_linsigma_M(COSMOPARAM *p, double z, double M) {
  return(get_linsigma(p,z,pow(3.*M/(4.*M_PI),1./3.)));
}

/* Linear fluctuation amplitude, Gaussian filtering.
 * Recomputes power spectrum if recomp!=0.
 */
double get_linsigma_gauss(COSMOPARAM *p, double z, double R, int recomp) {
  double k0, dlnk;
  long ik;
  static double D2[601];
  double T, var, kR;

  k0 = 1e-4;
  dlnk = 0.01*log(10.);
  if (recomp)
    get_Delta2k(p,z,k0,dlnk,601,D2);

  /* Integrate linear fluctuations */
  var = 0.;
  for(ik=0; ik<601; ik++) {
    kR = k0*exp(ik*dlnk)*R;
    T = exp(-kR*kR/2.);
    var += T*T*D2[ik];
  }
  var *= dlnk;
  return(sqrt(var));
}

/* Nonlinear power spectrum at specified z, range of k: D2[0..nk-1]
 * Uses Smith et al
 * BAO smoothing implemented *unless* k<0.
 */
void get_Delta2kNL(COSMOPARAM *p, double z, double k0, double dlnk, long nk, double *D2) {
  int sflag = 1;
  double R_NL, factor_NL, sigma, n, C;
  double k,y;
  long ik,ikv;
  double an,bn,cn,alpha,beta,gamma,mu,nu;
  double f1,f2,f3,Om;
  double kv, vSigma, vD2[78], x, s;
  double D2Q,D2H;
  double *D2lin, *D2lin_nowiggle;
  COSMOPARAM p2;

  /* BAO wiggle smoothing flag */
  if (k0<0) sflag=0;
  k0=fabs(k0);

  /* Use linear at z>6 */
  if (z>6) {
    get_Delta2k(p,z,k0,dlnk,nk,D2);
    return;
  }

  /* Linear and no-wiggle power spectra */
  D2lin = dvector(0,nk-1);
  D2lin_nowiggle = dvector(0,nk-1);
  get_Delta2k(p,z,k0,dlnk,nk,D2lin);
  if (sflag) {
    get_Delta2k(p,z,-k0,dlnk,nk,D2lin_nowiggle);
  } else {
    for(ik=0;ik<nk;ik++)
       D2lin_nowiggle[ik] = D2lin[ik];
  }

  /* Get nonlinear scale */
  R_NL = 10.;
  get_linsigma_gauss(p,z,R_NL,1);
  while(get_linsigma_gauss(p,z,R_NL,0)>1.) R_NL *= 2.;
  while(get_linsigma_gauss(p,z,R_NL,0)<1.) R_NL /= 2.;
  factor_NL = sqrt(2.);
  while(factor_NL>1.000001) {
    sigma = get_linsigma_gauss(p,z,R_NL,0);
    R_NL *= sigma>1? factor_NL: 1./factor_NL;
    factor_NL = sqrt(factor_NL);
  }

  /* Slope at nonlinear scale */
  n = -3. + 2.*(get_linsigma_gauss(p,z,0.9*R_NL,0)-get_linsigma_gauss(p,z,1.1*R_NL,0))/0.2;

  /* Curvature at nonlinear scale */
  C = -log(get_linsigma_gauss(p,z,exp(0.2)*R_NL,0)*get_linsigma_gauss(p,z,exp(-0.2)*R_NL,0))/0.02;

#if 0
  printf("z=%5.3lf R_NL=%12.5le sigma=%12.5le n=%9.6lf C=%9.6lf\n", z, R_NL, get_linsigma_gauss(p,z,R_NL,0), n, C);
#endif

  /* Smith et al parameters */
  an = pow(10., 1.4861+n*( 1.8369+n*( 1.6762+n*( 0.7940+n*0.1670)))-0.6206*C);
  bn = pow(10., 0.9463+n*( 0.9466+n*0.3084)-0.9400*C);
  cn = pow(10.,-0.2807+n*( 0.6669+n*0.3214)-0.0793*C);
  alpha = 1.3884+n*( 0.3700-n*0.1452);
  beta  = 0.8291+n*( 0.9854+n*0.3401);
  gamma = 0.8649 + 0.2989*n + 0.1631*C;
  mu = pow(10., -3.5442 + 0.1908*n);
  nu = pow(10.,  0.9589 + 1.2857*n);

  /* Matter density dependent functions */
  Om = get_Om(p,z);
#if 0
  set_default_cosmology(&p2);
  Om = get_Om(&p2,z);
#endif
  f1 = pow(Om, -0.0307);
  f2 = pow(Om, -0.0585);
  f3 = pow(Om,  0.0743);

  /* RMS displacement per axis (for BAO smoothing) */
  vSigma = 0.;
  s = getSoundHorizon(p);
  get_Delta2k(p,z,-2e-4,0.14,78,vD2);
  for(ikv=0;ikv<78;ikv++) {
    kv = 2e-4 * exp(0.14*ikv);
    x = kv * s;
    vSigma += vD2[ikv]/(kv*kv) * 2. * (1./3. +(- sin(x)/x - 2*cos(x)/x/x + 2*sin(x)/x/x/x)*exp(-0.005*x*x));
  }
  vSigma = sqrt(vSigma*0.14);

  /* Now loop over wavenumbers */
  for(ik=0; ik<nk; ik++) {
    k = k0*exp(ik*dlnk);
    y = k*R_NL;

    /* Quasilinear contribution */
    D2Q = D2lin_nowiggle[ik] * pow(1.+D2lin_nowiggle[ik], beta) / (1.+alpha*D2lin_nowiggle[ik])
          * exp(-y/4.*(1.+y/2.));

    /* Halo-halo contribution */
    D2H = an*pow(y,3.*f1) / (1.+bn*pow(y,f2)+pow(cn*f3*y,3.-gamma))
          / (1. + (mu+nu/y)/y);

    D2[ik] = D2Q+D2H + (D2lin[ik]-D2lin_nowiggle[ik])*exp(-k*k*vSigma*vSigma/2.);
  }

#if 0
  fprintf(stderr, "D2(%12.5le)=%12.5le .. D2(%12.5le)=%12.5le\n", k0, D2[0], k0*exp((nk-1)*dlnk), D2[nk-1]);
#endif

  free_dvector(D2lin,0,nk-1);    
  free_dvector(D2lin_nowiggle,0,nk-1);    
}

/* mass -> nu relation. */
double get_nu(COSMOPARAM *p, double M, double z) {
  double dc, ratio, sigma;
  dc = 1.686;
  sigma = get_linsigma_M(p,z,M);
  ratio = dc/sigma;
  return(ratio*ratio);
}

/* mass function f(M) = d[prob]/dM */
double get_fM(COSMOPARAM *p, double M, double z) {
  double Anorm, a, pp, f_nu, nu, dnudM, epsilonX, sigma, dsigmadM;
  a = 0.707;
  pp = 0.3;
  Anorm = 0.322;
  epsilonX = 0.02;
#ifdef SHETH_TORMEN
  nu = get_nu(p,M,z);
  f_nu = Anorm*(1.+pow(a*nu,-pp))*sqrt(a/2./M_PI/nu)*exp(-a*nu/2.);
  dnudM = (get_nu(p,M*(1.+epsilonX),z)-get_nu(p,M*(1.-epsilonX),z))/(2.*M*epsilonX);
  return(f_nu*dnudM);
#endif
  sigma = get_linsigma_M(p,z,M);
  dsigmadM = (get_linsigma_M(p,z,M*(1.+epsilonX))-get_linsigma_M(p,z,M*(1.-epsilonX)))/(2.*M*epsilonX);
  return(0.301/sigma*fabs(dsigmadM)*exp(-pow(fabs(log(1./sigma)+0.64),3.82)));
}

/* bias b(M) */
double get_bM(COSMOPARAM *p, double M, double z) {
  double Anorm, a, pp, dc, nu, b;
  a = 0.707;
  pp = 0.3;
  Anorm = 0.322;

  dc = 1.686;
  nu = get_nu(p,M,z);
  b = 1. + (a*nu-1.)/dc + 2.*pp/dc/(1.+pow(a*nu,pp));
  return(b);
}

/* Concentration at delta_vir */
double getconc(COSMOPARAM *p, double M, double z) {
  double M15 = M*0.0002775*p->omh2/(p->h*p->h);
  return(8.5*pow(M15,-0.086)*pow(1.+z,-0.65));
}

/* delta_virial */
double get_delta_virial(COSMOPARAM *p, double M, double z) {
  double om_then,x;

  om_then = get_Om(p,z);
  x = om_then-1.;
  return((18*M_PI*M_PI+x*(82-39*x))/om_then);
}

/* Returns mass at x times the virial radius */
double mass_other_radius(COSMOPARAM *p, double Mv, double z, double x) {
  double c, m1, m2;

  c = getconc(p,Mv,z);
  m1 = log(1+c) - c/(1+c);
  m2 = log(1+c*x) - c*x/(1+c*x);
  return(Mv * m2/m1);
}

/* Returns mass of halo with virial mass Mv at mean SO threshold Delta.
 * Note can only go a factor of 16 from virial radius
 */
double mvir2mso(COSMOPARAM *p, double Mv, double z, double Delta) {
  double d, dv, x, M, factor;

  M = Mv;
  dv = get_delta_virial(p,M,z);
  x = 1;
  for(factor=4.; factor>1.0000001; factor=sqrt(factor)) {
    M = mass_other_radius(p,Mv,z,x);
    d = M/Mv * dv / (x*x*x);
    x *= d>Delta? factor: 1./factor;
  }
  return(M);
}

double mso2mvir(COSMOPARAM *p, double M, double z, double Delta) {
  double Mv, Mtest, factor=0;

  Mv = M;
  while(fabs(factor-1)>1e-6) {
    Mtest = mvir2mso(p,Mv,z,Delta);
    factor = Mtest/M;
    Mv /= factor;
  }
  return(Mv);
}

/* Fourier transform of halo profile @ given k.
 * Normalize to 1 at k=0.  Note "M" = SO(178) mass.
 */
double tkhalo(COSMOPARAM *p, double M, double z, double k) {
  double F0, F1, Rv, Mv, R180, ksc, x, Rsc;
  double rhox2;
  int i;

  Mv = mso2mvir(p,M,z,180);
  Rv = pow(3.*Mv/(4.*M_PI) / get_delta_virial(p,M,z), 1./3.);
  ksc = k * Rv / getconc(p,Mv,z);
  R180 = pow(3.*M/(4.*M_PI) / 179., 1./3.);
  Rsc = R180 * getconc(p,Mv,z) / Rv;

  F0=F1=0.;
  for(i=0; i<100; i++) {
    x = (double)i/100. * Rsc;
    rhox2 = x/((1.+x)*(1.+x));
    F0 += rhox2;
    F1 += rhox2 * (i>0? sin(ksc*x)/(ksc*x): 1.);
  }
  return(F1/F0);
}

/* Get mass threshold of haloes with specified number density cutoff
 * (in comoving Mpc^3 times mean density)
 */
double get_mcut(double n, double z, COSMOPARAM *p) {
  double M, n_actual_old, n_actual;
  double fM;
  double lnMstep = 0.04;

  /* Start integrating downward until we reach target density */
  M = 2.5e6; n_actual = 0;
  while(n_actual<n) {
    fM = get_fM(p, M/exp(0.5*lnMstep), z);
    n_actual_old = n_actual;
    n_actual += fM * lnMstep;
    M/=exp(lnMstep);
  }

  /* Now that we have enough haloes, figure out where the threshold was */
  return(M*exp(lnMstep*(n-n_actual)/(n_actual_old-n_actual)));
}

/* Get bias and stochasticity of central galaxies at M>Mcut in halo model picture.
 * Note definition includes galaxy shot noise.
 * Assumes occupation probability po per halo.
 */
void getbr_halo(double Mcut, double z, double k, double po, COSMOPARAM *p, double *b, double *r) {
  double Pgg, Pgm, Pmm, Plin;
  double M, dn, n, nMT, nb;
  int i;
  double lnMstep = 0.04;
  int Nstep = 300; /* 12 e-folds */

  /* Mass PS: use NL spectrum */
  get_Delta2kNL(p, z, k, 0, 1, &Pmm);  
  Pmm *= 2.*M_PI*M_PI/(k*k*k);

  /* Linear PS */
  get_Delta2kNL(p, z, k, 0, 1, &Plin);
  Plin *= 2.*M_PI*M_PI/(k*k*k);

  /* Integrate key quantities over the mass fcn: n, nMT(k), nb */
  n = nMT = nb = 0.;
  for(i=0;i<Nstep;i++) {
    M = Mcut*exp((i+0.5)*lnMstep);
    dn = get_fM(p,M,z) * lnMstep;
    n += dn;
    nMT += dn * M * tkhalo(p,M,z,k);
    nb += dn * get_bM(p,M,z);
  }

  /* Galaxy PS */
  Pgg = nb*nb/n/n * Plin + 1./(n*po);
  Pgm = nb/n * Plin + nMT/n;

  if (b!=NULL) *b = sqrt(Pgg/Pmm);
  if (r!=NULL) *r = Pgm/sqrt(Pgg*Pmm);
}

/* Lensing strength g(z1=lens,z2=source), in Mpc.
 */
double get_lensing_strength(COSMOPARAM *p, double z1, double z2) {
  double K_chi2,r,chi;

  if (z1>z2) return(0.);

  r = 0.;
  chi = get_DRC(z2,p)-get_DRC(z1,p);
  K_chi2 = -p->Omega_K * p->h * p->h * chi * chi / (HL_HIMPC*HL_HIMPC);
  if (fabs(K_chi2)<1e-8) {
    r = chi * (1. - K_chi2/6.*(1.-K_chi2/20.));
  } else if (K_chi2>0) {
    r = chi * sin(sqrt(K_chi2))/sqrt(K_chi2);
  } else if (K_chi2<0) {
    r = chi * sinh(sqrt(-K_chi2))/sqrt(-K_chi2);
  }

  return(2.*r*get_DAC(z1,p)/get_DAC(z2,p));
}

/* Weak lensing power spectrum, for nl multipoles,
 * l=l0 ... l0*exp((nl-1)*dlnl), redshifts z1,z2 into
 * Cl[0..nl-1].
 */
void get_wl_ps(COSMOPARAM *p, long nl, double l0, double dlnl,
  double z1, double z2, double *Cl) {
#define DLNALENS 0.05
  double ztemp,z,r, val_4pga2r, g1, g2, l, k, dchi;
  long iz,izmax,il;
  double *Delta2;

  Delta2 = dvector(0,nl-1);

  /* Get maximum lens plane number */
  if (z1>z2) {
    ztemp = z1;
    z1 = z2;
    z2 = ztemp;
  }
  izmax = (long)floor(log(1.+z1)/DLNALENS);

  /* Initialize lensing power spectra */
  for(il=0; il<nl; il++)
    Cl[il] = 0.;

  /* Loop over lens planes */
  for(iz=1; iz<=izmax; iz++) {
    z = exp(iz*DLNALENS)-1.;
    r = get_DAC(z,p);
    get_Delta2kNL(p, z, l0/r, dlnl, nl, Delta2);

    /* Geometric factors */
    val_4pga2r = 1.5*p->omh2*(1.+z)/HL_HIMPC/HL_HIMPC;
    g1 = get_lensing_strength(p,z,z1);
    g2 = get_lensing_strength(p,z,z2);
    dchi = DLNALENS*(1.+z) / get_H(z,p);

    for(il=0; il<nl; il++) {
      l = l0*exp(il*dlnl);
      k = l/r;
      Cl[il] += Delta2[il] * 2.*M_PI*M_PI/pow(k,3.) * val_4pga2r * val_4pga2r * g1 * g2
                * dchi / (r*r) /4.;
    }
  }
  free_dvector(Delta2,0,nl-1);
#undef DLNALENS
}

/* Gets tomographic lensing cross-spectra.
 *
 * Cl[il + TINDEX(b1,b2)*nl] = convergence power spectrum between bins b1,b2
 * multipole bin il.
 *
 * Cl[start_gg + il + TINDEX(b1,b2)*nl] = galaxy power spectrum between bins b1,b2
 * multipole bin il.
 *
 * Cl[start_gk + il + nl*(b1+b2*nz)] = galaxies(b1) * convergence(b2).
 *
 * Covariance matrix ClCov[][] is included if not null.  Based on
 * gamma^2/n for each slice, and fsky.
 *
 * Non-Gaussianity included for matter only out to z=3.
 */
void get_wl_ps_many(COSMOPARAM *p, long nl, double l0, double dlnl,
  long nz, double *z, double *noise, double fsky, double *Cl,
  double **ClCov, double *auxparam) {
#define DZLENS 0.05
#define NGZMAX 3.0
  double nmodes,l,k;
  double *Clplus;
  long b1,b2,ti,nzl,i,j,b3,b4,ti2,il,il2;
  double M, dlnM, zz, r, Nh, *Cl1h, val_4pga2r, factor;
  double *T, *P, *dr, *dz;
  long start_gg, start_gk, nxtot;
  double P02,kref;
  double corr=0.5;
  long izmax;
  double *Delta2, *g, dchi;
  double nP;

  /* Parameters for lensing of the CMB.
   * Select column from file:
   * nc[0] = Planck TT
   * nc[1] = 2' beam, 18   uK' pol
   * nc[2] = 2' beam,  9   uK' pol
   * nc[3] = 2' beam,  1.4 uK' pol
   */
#define L_CMBLENS_MAX 2048
  static double noiseCMB[L_CMBLENS_MAX+1];
  double nc[4];
  long Lc, thisL;
  FILE *fp;
  noiseCMB[0] = noiseCMB[1] = 1e6;
  fp = fopen("cmb/lensing.dat", "r");
  for(Lc=2; Lc<=L_CMBLENS_MAX; Lc++) {
    fscanf(fp, "%ld %lg %lg %lg %lg", &thisL, nc, nc+1, nc+2, nc+3);
    noiseCMB[Lc] = nc[2];
  }
  fclose(fp);

  nzl = NCROSS(nz)*nl;
  start_gg = nzl;
  start_gk = start_gg + nzl;
  nxtot = start_gk + nl*nz*nz;
  dr = dvector(0, nz-1);

  /* Clear power spectrum */
  for(i=0; i<nxtot; i++)
    Cl[i] = 0.;

  /* Thicknesses of bins -- needed for galaxies */
  dz = dvector(0, nz-1);
  dz[0] = z[1]-z[0];
  dz[nz-1] = z[nz-1]-z[nz-2];
  for(b1=1; b1<nz-1; b1++)
    dz[b1] = (z[b1+1]-z[b1-1])/2.;
  if (z[nz-1]>1000)
    dz[nz-2] = z[nz-2]-z[nz-3];
  for(b1=0; b1<nz; b1++)
    dr[b1] = dz[b1]/get_H(z[b1],p);
  free_dvector(dz, 0, nz-1);

  /* Get cross-spectra */
#define DLNALENS 0.05
  for(b1=0; b1<nz; b1++)
    for(b2=0; b2<=b1; b2++)
      for(il=0;il<nl;il++)
        Cl[il+nl*TINDEX(b1,b2)]=0;
  izmax = (long)floor(log(1.+z[nz-1])/DLNALENS);
  Delta2 = dvector(0,nl-1);
  g = dvector(0,nz-1);
  for(i=1;i<=izmax;i++) {
    zz = exp(i*DLNALENS)-1.;
    r = get_DAC(zz,p);
    get_Delta2kNL(p, zz, l0/r, dlnl, nl, Delta2);

    /* Geometric factors */
    val_4pga2r = 1.5*p->omh2*(1.+zz)/HL_HIMPC/HL_HIMPC;
    for(b1=0; b1<nz; b1++)
      g[b1] = get_lensing_strength(p,zz,z[b1]);
    dchi = DLNALENS*(1.+zz) / get_H(zz,p);

    for(b1=0; b1<nz; b1++)
      for(b2=b1; b2<nz; b2++)
        for(il=0; il<nl; il++) {
          l = l0*exp(il*dlnl);
          k = l/r;
          Cl[il+nl*TINDEX(b1,b2)] += Delta2[il] * 2.*M_PI*M_PI/(k*k*k) * val_4pga2r * val_4pga2r * g[b1] * g[b2]
                                     * dchi / (r*r) /4.;
        }
  }
  free_dvector(Delta2,0,nl-1);
  free_dvector(g,0,nz-1);
#undef DLNALENS

  for(b1=0; b1<nz; b1++)
    for(b2=0; b2<=b1; b2++) {
      get_wl_ps(p,nl,l0,dlnl,z[b1],z[b2],Cl+nl*TINDEX(b1,b2));
    }

  /* Galaxy spectra, assume bias=1 -- though this won't matter. */
  P = dvector(0,nl-1);
  for(b1=0; b1<nz; b1++) {

    val_4pga2r = 1.5*p->omh2*(1.+z[b1])/HL_HIMPC/HL_HIMPC;
    r = get_DAC(z[b1],p);
    get_Delta2kNL(p,z[b1],l0/r,dlnl,nl,P);
    for(il=0; il<nl; il++) {
      l = l0*exp(il*dlnl);
      k = l/r;
#if 1
      corr = 0.75-0.3/M_PI*atan(log(P[il]));
#if 0
      kref = 0.2/p->h;
      get_Delta2kNL(p,z[b1],kref,0,1,&P02);
      corr = sqrt(P[il]/k/k/k/(P[il]/k/k/k + P02/kref/kref/kref/(3.0))); /* (3.0) = nP @ kref */
      if (corr<0.6) corr=0.6;
      if (corr>0.96) corr=0.96;
#endif
      /* Shot noise correction -- crude.  2.53e7 is # galaxies per sr per slice */
      nP = 2.53e7 / (r*r*dr[b1]) * P[il] * 2.*M_PI*M_PI/(k*k*k);
      corr *= sqrt(nP/(1.+nP));
#endif
#if 1
      corr = 1e-4;
#endif
      /* Knock out correlations at surface of last scatter */
      if (z[b1]>1000) corr = 1e-3;

      /* Compute galaxy auto-spectra -- xspec=0. */
      Cl[start_gg + il + TINDEX(b1,b1)*nl] = P[il] * 2.*M_PI*M_PI/(k*l*l*dr[b1]) / corr / corr;

      /* Galaxy convergence x-spec */
      for(b2=b1+1; b2<nz; b2++)
        Cl[start_gk + il + nl*(b1+b2*nz)] = P[il] * M_PI*M_PI/(k*l*l) * val_4pga2r * get_lensing_strength(p,z[b1],z[b2]);
    }
  }
  free_dvector(P,0,nl-1);

  /* Covariance, if requested */
  if (ClCov!=NULL) {
    /* Covariance matrix augmented with noise */
    Clplus = dvector(0, NCROSS(nz)-1);

    /* Initialize covariance to zero */
    for(i=0;i<nzl;i++)
      for(j=0;j<nzl;j++)
        ClCov[i][j] = 0.;

    for(il=0; il<nl; il++) {
      l = l0*exp(il*dlnl);
      nmodes = fsky * l * l * 2. * dlnl; /* # modes/bin */
      for(b1=0; b1<nz; b1++)
        for(b2=0; b2<=b1; b2++) {
          ti = TINDEX(b1,b2);
          Clplus[ti] = Cl[il+nl*ti] + (b1==b2? noise[b1]: 0);

          /* Replace noise with CMB lensing curve if at z>1000 */
          if (b1==b2 && z[b1]>1000) {
            Lc = (long)floor(l);
            Clplus[ti] = Cl[il+nl*ti] + (Lc<L_CMBLENS_MAX? noiseCMB[Lc]: 1e6*noiseCMB[L_CMBLENS_MAX]);
fprintf(stderr, "bins %2ld,%2ld: l=%6.1lf, Cl=%12.5le with noise=%12.5le\n", b1, b2, l, Cl[il+nl*ti], Clplus[ti]);
          }
        }

      /* Weak lensing covariance */
      for(b1=0; b1<nz; b1++)
        for(b2=0; b2<=b1; b2++) {
          ti = TINDEX(b1,b2);
          for(b3=0; b3<nz; b3++)
            for(b4=0; b4<=b3; b4++) {
              ti2 = TINDEX(b3,b4);
              ClCov[il+nl*ti][il+nl*ti2] =(Clplus[TINDEX(b1,b3)]*Clplus[TINDEX(b2,b4)]
                                          +Clplus[TINDEX(b1,b4)]*Clplus[TINDEX(b2,b3)])
                                          /nmodes;
            }
        } /* end for b2 */

      /* Galaxy covariance */
      for(b1=0; b1<nz; b1++)
        for(b2=0; b2<=b1; b2++) {
          ti = TINDEX(b1,b2);
          ClCov[start_gg+il+nl*ti][start_gg+il+nl*ti] = (b1==b2? 2: 1)/nmodes
            * Cl[start_gg+il+nl*TINDEX(b1,b1)]*Cl[start_gg+il+nl*TINDEX(b2,b2)];
        }

      /* Cross-correlation covariance */
      for(b1=0; b1<nz; b1++)
        for(b2=0; b2<nz; b2++)
          for(b3=0; b3<nz; b3++)
            for(b4=0; b4<nz; b4++)
              ClCov[start_gk+il+nl*(b1+b2*nz)][start_gk+il+nl*(b3+b4*nz)]
                = (Cl[start_gg+il+nl*TINDEX(b1,b3)]*Clplus[TINDEX(b2,b4)]
                  +Cl[start_gk+il+nl*(b1+b4*nz)]*Cl[start_gk+il+nl*(b3+b2*nz)])/nmodes;

      /* gg x kk */
      for(b1=0; b1<nz; b1++)
        for(b2=0; b2<=b1; b2++) {
          ti = TINDEX(b1,b2);
          for(b3=0; b3<nz; b3++)
            for(b4=0; b4<=b3; b4++) {
              ti2 = TINDEX(b3,b4);
              ClCov[start_gg+il+nl*ti][il+nl*ti2] = ClCov[il+nl*ti2][start_gg+il+nl*ti]
                = (Cl[start_gk+il+nl*(b1+b3*nz)]*Cl[start_gk+il+nl*(b2+b4*nz)]
                  +Cl[start_gk+il+nl*(b1+b4*nz)]*Cl[start_gk+il+nl*(b2+b3*nz)])/nmodes;
              }
          }

      /* kk x gk */
      for(b1=0; b1<nz; b1++)
        for(b2=0; b2<=b1; b2++) {
          ti = TINDEX(b1,b2);
          for(b3=0; b3<nz; b3++)
            for(b4=0; b4<nz; b4++)
              ClCov[il+nl*ti][start_gk+il+nl*(b3+b4*nz)] = ClCov[start_gk+il+nl*(b3+b4*nz)][il+nl*ti]
                = (Clplus[TINDEX(b1,b4)]*Cl[start_gk+il+nl*(b3+b2*nz)]
                  +Clplus[TINDEX(b2,b4)]*Cl[start_gk+il+nl*(b3+b1*nz)])/nmodes;
        }

      /* gg x gk */
      for(b1=0; b1<nz; b1++)
        for(b2=0; b2<=b1; b2++) {
          ti = TINDEX(b1,b2);
          for(b3=0; b3<nz; b3++)
            for(b4=0; b4<nz; b4++)
              ClCov[start_gg+il+nl*ti][start_gk+il+nl*(b3+b4*nz)] = ClCov[start_gk+il+nl*(b3+b4*nz)][start_gg+il+nl*ti]
                = (Cl[start_gg+il+nl*TINDEX(b1,b3)]*Cl[start_gk+il+nl*(b2+b4*nz)]
                  +Cl[start_gg+il+nl*TINDEX(b2,b3)]*Cl[start_gk+il+nl*(b1+b4*nz)])/nmodes;
        }
    } /* end for il */

    /* Halo shot noise */
    dlnM = 0.1*log(10.);
    T = dvector(0,nl-1);
    Cl1h = dvector(0,nz-1);
    izmax = floor(z[nz-1]/DZLENS);
    izmax = izmax<floor(NGZMAX/DZLENS)? izmax: floor(NGZMAX/DZLENS);
    for(i=1; i<izmax; i++) {
      zz = i*DZLENS;
      r = get_DAC(zz,p);
      val_4pga2r = 1.5*p->omh2*(1.+zz)/HL_HIMPC/HL_HIMPC;

      /* Halo transfer functions */
      for(M=1e-3; M<1e6; M*=exp(dlnM)) {
        for(il=0;il<nl;il++)
          T[il] = tkhalo(p,M,zz,l0*exp(il*dlnl)/r);

        /* Number of haloes, full sky, this bin */
        Nh = 4.*M_PI*r*r*DZLENS/get_H(zz,p)*get_fM(p,M,zz)*dlnM;

        /* int(kappa d Omega) from a single halo */
        for(b1=0; b1<nz; b1++)
          Cl1h[b1] = val_4pga2r*get_lensing_strength(p,zz,z[b1])*M/(r*r);

        factor = 1./(fsky * 16. * M_PI * M_PI);
#if 0
        if (i==2) printf("z=%7.4lf dz=%7.4lf M=%12.5le dlnM=%7.4lf Nh=%12.5le intkappa=%12.5le T5=%9.6lf\n", zz, DZLENS, M, dlnM, Nh, Cl1h[nz-1], T[5]);
#endif

#if 0
        if ((zz<0.7001) && (M>2e14/2.775e11/p->omh2)) Nh /= 5.0;
#endif

#ifndef GAUSS_CONVERGENCE
        /* 4pt fcn */
        for(il=0;il<nl;il++)
          for(b1=0; b1<nz; b1++)
            for(b2=0; b2<=b1; b2++) {
              ti = TINDEX(b1,b2);
              for(il2=0;il2<nl;il2++)
                for(b3=0; b3<nz; b3++)
                  for(b4=0; b4<=b3; b4++) {
                    ti2 = TINDEX(b3,b4);
                    ClCov[il+nl*ti][il2+nl*ti2] += Nh * Cl1h[b1] * Cl1h[b2] * Cl1h[b3] * Cl1h[b4] * T[il] * T[il] * T[il2] * T[il2] * factor;
                  }
            }
#endif
      }
    }
    free_dvector(T,0,nl-1);
    free_dvector(Cl1h,0,nz-1);

    free_dvector(Clplus, 0, NCROSS(nz)-1);
  }
  free_dvector(dr,0,nz-1);
#undef DZLENS
#undef NGZMAX
}

/* Gets CCC Fisher matrix for lensing strengths g_ab.
 * Requires input Cl, and done for galaxies in bin b1.
 *
 * Output -> Fgg[0..nz-1][0..nz-1].
 */
void get_fisher_ccc1(COSMOPARAM *p, long nl, double l0, double dlnl,
  long nz, double *z, double *noise, double fsky, double *Cl, double rcorr,
  long b1, double **Fgg) {

  double **C;
  long il, nsp, db2, db3, b2, b3;
  double l, k, Nmodes, D2__Pdchi, D, D2__dchi, dz, val_4pga2r, denom;
  double *mpower, *gref, *fgref;

  /* Clear output */
  for(b2=0;b2<nz;b2++)
    for(b3=0;b3<nz;b3++)
      Fgg[b2][b3] = 0.;

  /* Need 2 bins for ccc */
  nsp = nz-b1-1;
  if (nsp<=1) return;

  /* Distances */
  dz = z[b1+1]-z[b1];
  D = get_DAC(z[b1],p);
  val_4pga2r = 1.5*p->omh2*(1.+z[b1])/HL_HIMPC/HL_HIMPC;
  D2__dchi = D*D/dz*get_H(z[b1],p);
  mpower = dvector(0,nl-1);
  get_Delta2kNL(p,z[b1],l0/D,dlnl,nl,mpower);

  for(il=0; il<nl; il++) {
    l = l0*exp(il*dlnl);
    k = l/D;
    D2__Pdchi = D2__dchi * k*k*k/(2.*M_PI*M_PI*mpower[il]);
    Nmodes = l*l*2.*fsky*dlnl;
    C = dmatrix(0,nsp-1,0,nsp-1);

    /* Get the C-inverse for the lensing strengths */
    for(db2=0; db2<nsp; db2++) {
      b2 = b1+db2+1;
      for(db3=0; db3<nsp; db3++) {
        b3 = b1+db3+1;
        C[db2][db3] = 4.*D2__Pdchi/(val_4pga2r*val_4pga2r*rcorr*rcorr*Nmodes)
                      * (Cl[il+nl*TINDEX(b2,b3)] + (b2==b3? noise[b2]: 0.));
      }
    }
    gaussjinv(C,nsp);

    /* Insert marginalization over rcorr here */
    gref = dvector(0,nsp-1);
    fgref = dvector(0,nsp-1);
    for(db2=0; db2<nsp; db2++) {
      b2 = b1+db2+1;
      gref[db2] = get_lensing_strength(p,z[b1],z[b2]);
    }
    denom = 0.;
    for(db2=0; db2<nsp; db2++)
      for(db3=0; db3<nsp; db3++) {
        fgref[db2] += C[db2][db3]*gref[db3];
        denom += C[db2][db3]*gref[db2]*gref[db3];
      }
    for(db2=0; db2<nsp; db2++)
      for(db3=0; db3<nsp; db3++)
        C[db2][db3] -= fgref[db2]*fgref[db3]/denom;
    free_dvector(gref,0,nsp-1);
    free_dvector(fgref,0,nsp-1);

    for(db2=0; db2<nsp; db2++) {
      b2 = b1+db2+1;
      for(db3=0; db3<nsp; db3++) {
        b3 = b1+db3+1;
        Fgg[b2][b3] += C[db2][db3];
      }
    }

    free_dmatrix(C,0,nsp-1,0,nsp-1);
  }
  free_dvector(mpower,0,nl-1);
}

/* ISW effect temperature perturbation.
 * Requires as input b^2*n, in Mpc^-3.
 * Output is RMS temperature per lnl from ISW and its uncertainty.
 */
#define LMAX 2000
void get_isw(COSMOPARAM *p, double l, double dlnl,
  double z, double dz, double b2n, double fsky, double *dT,
  double *sigT) {

  static int is_init=0;
  static double *Delta2TT_ref;
  long LL;
  double ctt,cte,cee;
  double Pk,k,dlnGdlna;
  double zp, zm;
  FILE *fp;

  /* Get CMB power spectrum */
  if (!is_init) {
    Delta2TT_ref = dvector(0,LMAX);

    fp = fopen("cmb/cmbout.ref", "r");
    while(fscanf(fp, "%ld %lg %lg %lg", &LL, &ctt, &cee, &cte)!=EOF) {
      if (LL<=LMAX) {
        Delta2TT_ref[LL] = ctt*LL*(LL+1.0)/(2.*M_PI);
#if 0
        if (LL>30) Delta2TT_ref[LL] *= 1.-cte*cte/ctt/(cee+4.2e-17);
#endif
      }
    }
    fclose(fp);
    Delta2TT_ref[0] = Delta2TT_ref[1] = Delta2TT_ref[2];
    is_init = 1;
  }

  /* Get wavenumber, power spectrum, growth rate */
  k = l/get_DAC(z,p);
  get_Delta2k(p,z,k,0,1,&Pk);
  Pk *= 2.*M_PI*M_PI/(k*k*k);

  zp = z+dz/2.;
  zm = z-dz/2.;
  dlnGdlna = log(DNget_growth_normhiz(p,zm)/DNget_growth_normhiz(p,zp))
             / log((1.+zp)/(1.+zm));

  *dT = 3./sqrt(2.*M_PI) * p->omh2 / (HL_HIMPC*HL_HIMPC)
        * sqrt(Pk)/k * sqrt(get_H(z,p)*dz) * (1.-dlnGdlna);

  *sigT = sqrt(
            (l<LMAX?Delta2TT_ref[(long)floor(l)]:1e49)
            * (1. + 1./(b2n*Pk))
            / (2*l*l*dlnl*fsky)
          );

}
#undef LMAX

/* CMB lensing delta kappa.
 * Requires as input b^2*n, in Mpc^-3.
 * Output is RMS kappa per lnl from this redshift range and its uncertainty.
 */
#define LMAX 2000
void get_kappa_cmblens(COSMOPARAM *p, double l, double dlnl,
  double z, double dz, double b2n, double fsky, double *dkappa,
  double *sigk) {

  FILE *fp;
  static int is_init2 = 0;
  static double *Delta2KK_ref;
  long LL;
  double k,Pk;
  static double info[7];

  if (!is_init2) {
    Delta2KK_ref = dvector(0,LMAX);

    fp = fopen("cmb/pact-lensing.dat", "r");
    while(fscanf(fp, "%ld %lg %lg %lg %lg %lg %lg %lg", &LL, info, info+1, info+2, info+3, info+4, info+5, info+6)!=EOF) {
      if (LL<=LMAX) {
//        Delta2KK_ref[LL] = info[6]*LL*(LL+1.0)/(2.*M_PI); /* [1]=TT; [2]=TE; [3]=EE; [4]=TB; [5]=EB; [6]=all */
        /* combined */
        Delta2KK_ref[LL] = info[3]*info[5]/(info[3]+info[5])*LL*(LL+1.0)/(2.*M_PI);
      }
    }
    fclose(fp);
    Delta2KK_ref[0] = Delta2KK_ref[1] = Delta2KK_ref[2];
    is_init2 = 1;
  }

  /* Get wavenumber, power spectrum, growth rate */
  k = l/get_DAC(z,p);
  get_Delta2kNL(p,z,k,0,1,&Pk);
  Pk *= 2.*M_PI*M_PI/(k*k*k);

  *dkappa = 3./sqrt(32.*M_PI) * p->omh2 / (HL_HIMPC*HL_HIMPC)
            * get_lensing_strength(p, z, 1100.)
            * (1.+z) * k*sqrt(Pk*dz/get_H(z,p));

  *sigk = sqrt(
            (l<LMAX?Delta2KK_ref[(long)floor(l)]:1e49)
            * (1. + 1./(b2n*Pk))
            / (2*l*l*dlnl*fsky)
          );

  /* Throw out nonlinear regime */
  if (k>0.1) *sigk=1e49;

}
#undef LMAX

/* velocity power spectrum and its uncertainty; uses no-wiggles */
double get_VPS(COSMOPARAM *p, double z, double k) {

  double dz;
  double zp, zm, dlnGdlna;
  double D2m, D2v;

  dz = z/10.;

  /* Value */
  zp = z+dz/2.;
  zm = z-dz/2.;
  dlnGdlna = log(DNget_growth_normhiz(p,zm)/DNget_growth_normhiz(p,zp))
             / log((1.+zp)/(1.+zm));
  get_Delta2k(p, z, -k, 0.01, 1, &D2m);
  D2v = D2m * dlnGdlna * dlnGdlna;

  return(D2v);
}
