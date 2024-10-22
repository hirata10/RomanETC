#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "nr_utils.h"
#include "utils.h"
#include "ps.h"
#include "fom.h"

#include "mnemonics.c"

/* Choice of parameterized HOD model space */
#define NHOD 3

/* Halo parameter space -- right now 1e9 .. 2.5e15 Msun */
#define M_MIN 0.027176
#define DLNM  0.115129254649702
#define NMASS 141

/* Linear parameters: central occupation probability & satellite number.
 * central[0 .. NMASS-1], satellite[0 .. NMASS-1], trunc[0 .. NMASS-1]
 * Trunc counts truncated NFW profiles of that mass (relative to number
 * of unstripped haloes)
 */
#define NPARAM (3*NMASS)

#define XKMIN 0.0316227766
#define XDLNK 0.230258509299405
#define XNK   36

/* Abundance of haloes in each bin */
void GetHaloAbundanceVector(double *n, double z, COSMOPARAM *p) {
  int i;
  double M;

  for(i=0;i<NMASS;i++) {
    M = M_MIN*exp(i*DLNM);
    n[i] = get_fM(p,M,z) * DLNM;
  }
}

/* Bias of haloes in each bin */
void GetHaloBiasVector(double *b, double z, COSMOPARAM *p) {
  int i;
  double M;

  for(i=0;i<NMASS;i++) {
    M = M_MIN*exp(i*DLNM);
    b[i] = get_bM(p,M,z);
  }
}

/* Fourier transform of halo profile @ given k.
 * Normalize to 1 at k=0.  Note "M" = SO(178) mass.
 * Truncates @ 0.4r_vir.
 */
double tkhaloTrunc(COSMOPARAM *p, double M, double z, double k) {
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
    if (i<40) F1 += rhox2 * (i>0? sin(ksc*x)/(ksc*x): 1.);
  }
  return(F1/F0);
}

/* Lensing signal, P(halo density, delta, k): dimensionless!
 * Returns vector for centrals & satellites, length NPARAM
 */
void GetHaloMassPk(double *Pk, double k, double z, COSMOPARAM *p) {
  int i;
  double M, n, b;
  double Plin, trf, trfsub;

  /* Linear power spectrum for 2 halo term */
  get_Delta2k(p,z,k,0,1,&Plin);
  Plin *= 2.*M_PI*M_PI/(k*k*k);

  for(i=0;i<NMASS;i++) {
    M = M_MIN*exp(i*DLNM);
    n = get_fM(p,M,z) * DLNM;
    b = get_bM(p,M,z);
    trf = tkhalo(p,M,z,k);
    trfsub = tkhaloTrunc(p,M,z,k);

    /* Central term: 2 halo + 1 halo */
    Pk[i] = n*( b*Plin + M*trf );

    /* Satellite term: additional halo transfer fcn */
    Pk[NMASS+i] = Pk[i]*trf;

    /* Truncated term: subhalo only, no 2-halo contribution */
    Pk[2*NMASS+i] = n*M*trfsub;
  }
}

/* Get inverse-variance on lensing signal.  Neglect 4pt contributions.
 * Formula:  Cl(cross) = dr/(r*r) * Pkcross * 0.5 * val_4pga2r * g
 * Variance[Cl(cross)] = Cl(galaxies) * Cl(shear) / sqrt(Nmode)
 * Cl(galaxies) = nL*nL * dr/(r*r) * [ 1/nL + b^2*Pkmatter ]
 * Cl(shear) = gamma^2/nS * [ 1 + (l/lc)^-1.3 ]
 * Nmode = Area * dlnk * k^2/(2pi)                     [ Area in Mpc^2 ]
 *
 * The S/N on Pkcross depends only on the total volume, so we will set
 * dr=1 and Area=Volume.
 * Also assume gamma_rms for the sources is 0.25.
 */
double GetIvarLensingPk(double k, double zl, double zs,
  double lc, double Volume, double nLens, double nSrc, double bLens,
  COSMOPARAM *p) {
#define GAMMARMS_SRC 0.25

  double r, Pm, factor;
  double Nmode, Cl_shear, Cl_gal, ivar_clcross;

  /* Physical parameters */
  r = get_DAC(zl,p);
  get_Delta2kNL(p,zl,k,0,1,&Pm);
  Pm *= 2.*M_PI*M_PI/(k*k*k);

  /* Fisher calculation */
  Nmode = Volume * XDLNK * k*k/(2.*M_PI);

  get_wl_ps(p,1,k*r,0,zs,zs,&Cl_shear);
  Cl_shear += GAMMARMS_SRC*GAMMARMS_SRC/nSrc;
#if 0
  Cl_shear = GAMMARMS_SRC*GAMMARMS_SRC/nSrc * (1. + pow(k*r/lc,-1.3));
#endif
  Cl_gal = nLens*nLens * (1./nLens + bLens*bLens*Pm)/(r*r);
  ivar_clcross = Nmode / (Cl_shear*Cl_gal);
  factor = 0.5 * 1.5*p->omh2*(1.+zl)/HL_HIMPC/HL_HIMPC * get_lensing_strength(p,zl,zs) / (r*r);

  return(ivar_clcross * factor * factor);
}

/* Adds Fisher matrix for gg-lensing to halo model.
 * Uses fiducial lens model x[].  Returns S/N ratio.
 */
double AddGGLensFisher(double **F, double zl, double zs, double lc, double Volume,
  double *x, double nSrc, COSMOPARAM *p) {

  int i,j,ik;
  double k, ivar, snr;
  double nLens, bLens;
  double dPdx[NPARAM];
  double abund[NMASS], bias[NMASS];
  double **FL;

  /* Lensing Fisher matrix */
  FL = dmatrix(0,NPARAM-1,0,NPARAM-1);

  /* Get lens density & bias */
  GetHaloAbundanceVector(abund,zl,p);
  GetHaloBiasVector(bias,zl,p);
  nLens = bLens = 0.;
  for(i=0;i<NMASS;i++) {
    nLens += abund[i]*(x[i]+x[NMASS+i]);
    bLens += abund[i]*(x[i]+x[NMASS+i])*bias[i];
  }
  bLens /= nLens;

  for(ik=0; ik<XNK; ik++) {
#if 0
    fprintf(stderr, "ik = %3d / %3d\n", ik, XNK);
#endif
    k = XKMIN * exp(XDLNK*ik);

    /* For each wavenumber, obtain P(k) for each halo mass ... */
    GetHaloMassPk(dPdx,k,zl,p);

    /* Get uncertainties & add to Fisher */
    ivar = GetIvarLensingPk(k,zl,zs,lc,Volume,nLens,nSrc,bLens,p);
    for(i=0;i<NPARAM;i++) for(j=0;j<NPARAM;j++) FL[i][j] += ivar*dPdx[i]*dPdx[j];
  }

  /* Get signal/noise ratio */
  snr = 0.;
  for(i=0;i<NPARAM;i++) for(j=0;j<NPARAM;j++) snr += FL[i][j]*x[i]*x[j];
  snr = sqrt(snr);

#if 1
  fprintf(stderr, "n=%12.5le, b=%9.6lf, S/N=%12.5le\n", nLens, bLens, snr);
#endif

#if 0
#define SNRMAX 17.44
#ifdef SNRMAX
  if (snr>SNRMAX) for(i=0;i<NPARAM;i++) for(j=0;j<NPARAM;j++) FL[i][j] *= SNRMAX/snr*SNRMAX/snr;
#endif
#endif

  /* Add to total Fisher */
  for(i=0;i<NPARAM;i++) for(j=0;j<NPARAM;j++) F[i][j] += FL[i][j];
  free_dmatrix(FL,0,NPARAM-1,0,NPARAM-1);

  return(snr);
}

/* Adds Fisher matrix contribution from knowing the number density of galaxies
 * to accuracy ivar_ln_n and the bias to ivar_ln_b.
 */
void AddNFisher(double *x, double **F, double ivar_ln_n, double ivar_ln_b,
  double z, COSMOPARAM *p) {
  int i,j;
  double abund[NMASS], bias[NMASS];
  double N,Nb;
  double dNdx[NPARAM], dNbdx[NPARAM];
  double dlnNdx[NPARAM], dlnbdx[NPARAM];

  GetHaloAbundanceVector(abund,z,p);
  GetHaloBiasVector(bias,z,p);

  /* Derivative of N and Nb wrt each parameter */
  for(i=0;i<NMASS;i++) {
    dNdx[i]  = dNdx[i+NMASS]  = abund[i];
    dNbdx[i] = dNbdx[i+NMASS] = abund[i]*bias[i];
    dNdx[i+2*NMASS] = dNbdx[i+2*NMASS] = 0.;
  }
  N = Nb = 0.;
  for(i=0;i<NPARAM;i++) {
    N  += dNdx[i]  * x[i];
    Nb += dNbdx[i] * x[i];
  }

  /* Derivatives of ln N and ln b */
  for(i=0;i<NPARAM;i++) {
    dlnNdx[i] = dNdx[i]/N;
    dlnbdx[i] = dNbdx[i]/Nb - dlnNdx[i];
  }

  /* Fisher */
  for(i=0;i<NPARAM;i++)
    for(j=0;j<NPARAM;j++)
      F[i][j] += ivar_ln_n*dlnNdx[i]*dlnNdx[j] + ivar_ln_b*dlnbdx[i]*dlnbdx[j];
}

/* Fisher for cluster analysis.  Uses:
 * z = typical redshift
 * Vol = volume surveyed
 * clMmin = minimum estimated mass (in Msun)
 * sigma_disp = dispersion in estimated mass in dex (1 sigma)
 * sProb = selection probability for clusters (a Fisher penalty is applied!)
 * sysFloor = systematic error floor in sigma(Nbar)/Nbar per log10 M
 */
void AddClusterFisher(double *x, double **F, double z, double Vol, double clMmin,
  double sigma_disp, double sProb, double sysFloor, COSMOPARAM *p) {

  int i, j, nbi, nmax;
  double nbar;
  double bin_eff, sig_bin, Ncl, Nsat, Nsat2, IvarNsat, IvarSys, Ncl_unnorm, N_interloper;
  double P[NMASS], abund[NMASS];

  /* Dispersion in bin numbers */
  sig_bin = sigma_disp * log(10.)/DLNM;
  if (sig_bin<1.) {
    fprintf(stderr, "Warning: cluster mass estimator dispersion too low [%6.4lf], refine grid.\n", sig_bin);
  }

  GetHaloAbundanceVector(abund,z,p);
  nbar = 0.;
  for(i=0;i<NMASS;i++) nbar += abund[i]*(x[i]+x[i+NMASS]);

  /* Sum over cluster estimated mass bins */
  nmax = NMASS - log(clMmin/M_MIN/2.775e11/p->omh2)/DLNM + 3*sig_bin;
  for(nbi=0;nbi<nmax;nbi++) {
    bin_eff = log(clMmin/M_MIN/2.775e11/p->omh2)/DLNM + 0.5 + nbi;

    N_interloper = nbar * 2*M_PI*1.0*1.0*0.005*(1.+z)/get_H(z,p);

    /* Probability of cluster in each true mass bin making it into this estimated bin */
    Ncl = Nsat = Nsat2 = 0.;
    for(i=0;i<NMASS;i++) {
      P[i] = exp(-(i-bin_eff)*(i-bin_eff)/2./sig_bin/sig_bin)/sqrt(2*M_PI)/sig_bin;
      Ncl += P[i]*abund[i];
      Nsat += P[i]*abund[i] * x[i+NMASS];
      Nsat2 += P[i]*abund[i] * x[i+NMASS] * (x[i+NMASS]+1.);
    }
    Ncl_unnorm = Ncl;
    Nsat /= Ncl;
    Nsat2 /= Ncl;
    Ncl *= sProb*Vol;

    IvarNsat = Ncl/(Nsat2-Nsat*Nsat + N_interloper); /* assumes Poisson errors */
    IvarSys = DLNM/log(10.)/(sysFloor*sysFloor*Nsat*Nsat);
    IvarNsat = IvarNsat*IvarSys/(IvarNsat+IvarSys);

#if 1
    fprintf(stderr, "%8.5lf %12.5le %12.5le %12.5le %12.5le %12.5le\n", (log(clMmin)+(0.5+nbi)*DLNM)/log(10.), Ncl, Nsat, sqrt(Nsat2-Nsat*Nsat), N_interloper, 1./sqrt(IvarNsat));
#endif

    /* Add to Fisher */
    for(i=0;i<NMASS;i++) for(j=0;j<NMASS;j++) F[i+NMASS][j+NMASS] += IvarNsat * P[i]*abund[i] * P[j]*abund[j] / (Ncl_unnorm*Ncl_unnorm);
  }

}

/* Galaxy HOD model space.  
 * [Cosmological parameters are fixed.]
 */
void getmodel(double *hp, double *x, double z, COSMOPARAM *p) {

  int i;
  double chi, peakwidth, ncent, nsat, epsilon, Rmin, y;
  double Mmin_Solar;
  double abund[NMASS];

  GetHaloAbundanceVector(abund,z,p);  

  /* Minimum halo mass in range, solar masses */
  Mmin_Solar = 2.775e11 * p->omh2 * M_MIN;

  /* Satellite power law */
  epsilon = NHOD>=4? hp[3]: 1.0;

  /* Minimum ratio for satellites (vs Mcentral) */
  Rmin = NHOD>=5? hp[4]: 3.0;

  /* Width of peak in dex, 1 sigma */
  peakwidth = NHOD>=6? hp[5]: 0.2;

  /* First 3 parameters in model space = {nbar, fsat, log10Mcent} */

  /* Build parts of HOD */
  ncent = nsat = 0.;
  for(i=0;i<NMASS;i++) {

    /* Central */
    chi = (log(Mmin_Solar) + i*DLNM - hp[2]*log(10.)) / (peakwidth*log(10.));
    x[i] = exp(-0.5*chi*chi)/abund[i];
    ncent += x[i]*abund[i];

    /* Satellite */
    y = Mmin_Solar*exp(i*DLNM)/pow(10.,hp[2]);
#define SATMIN_CUTOFF_STRENGTH 10.0
    x[i+NMASS] = y>Rmin? 1./(1.+pow(Rmin/y,SATMIN_CUTOFF_STRENGTH)): pow(y/Rmin,SATMIN_CUTOFF_STRENGTH)/(1.+pow(y/Rmin,SATMIN_CUTOFF_STRENGTH));
#undef SATMIN_CUTOFF_STRENGTH
    x[i+NMASS] *= pow(y,epsilon);
    nsat += x[i+NMASS]*abund[i];
  }

  /* Truncated (sub)haloes */
  for(i=0;i<NMASS;i++)
    x[i+2*NMASS] = x[i]*hp[0]*hp[1]/ncent;

  /* Normalize */
  for(i=0;i<NMASS;i++) {
    x[i] *= hp[0]*(1.-hp[1]) / ncent;
    x[i+NMASS] *= hp[0]*hp[1] / nsat;
  }

}

/* Project Fisher in Big Space down to HOD parameter space. */
void ProjectFisher(double **F, double **Fp, double *hp, double z, COSMOPARAM *p) {
  int i,j,ih,jh;
  double x[NPARAM];
  double hpPlus[NHOD], xPlus[NPARAM];
  double hpMinus[NHOD], xMinus[NPARAM];
  double **Jac;

  /* Fiducial step sizes, must change if we change the HOD model beyond # parameters */
  double stepSize[] = {0.01,0.01,0.01,0.01,0.01,0.01};

  /* Transformation matrix: HOD --> Big Space */
  Jac = dmatrix(0,NPARAM-1,0,NHOD-1);

  /* Fiducial model & Jacobian calculation */
  getmodel(hp, x, z, p);
  for(ih=0;ih<NHOD;ih++) {
    for(jh=0;jh<NHOD;jh++) hpPlus[jh] = hpMinus[jh] = hp[jh];
    hpPlus[ih] += stepSize[ih];
    hpMinus[ih] -= stepSize[ih];

    getmodel(hpPlus, xPlus, z, p);
    getmodel(hpMinus, xMinus, z, p);

    for(i=0;i<NPARAM;i++) Jac[i][ih] = (xPlus[i]-xMinus[i])/(2.*stepSize[ih]);
  }

  /* Matrix multiplication */
  for(ih=0;ih<NHOD;ih++) for(jh=0;jh<NHOD;jh++) {
    Fp[ih][jh] = 0.;
    for(i=0;i<NPARAM;i++) for(j=0;j<NPARAM;j++) Fp[ih][jh] += F[i][j]*Jac[i][ih]*Jac[j][jh];
  }

  free_dmatrix(Jac,0,NPARAM-1,0,NHOD-1);
}

int main(int argc, char **argv) {
  FILE *fp;
  int i,j,ik;
  COSMOPARAM p;
  double z = 0.07, k;
  double abund[NMASS], bias[NMASS], hmpk[NPARAM], hp[NHOD];
  double *x;
  double **F, **Fp, **biasHessian, **bHH;
  double varb;
  double nl, sigma2, delta2;

  set_default_cosmology(&p);
  F = dmatrix(0,NPARAM-1,0,NPARAM-1);
  bHH = dmatrix(0,NPARAM-1,0,NPARAM-1);
  Fp = dmatrix(0,NHOD-1,0,NHOD-1);
  biasHessian = dmatrix(0,NHOD-1,0,NHOD-1);

  GetHaloAbundanceVector(abund,z,&p);
  GetHaloBiasVector(bias,z,&p);

  k = 5.;
  GetHaloMassPk(hmpk,k,z,&p);

#if 0
  sigma2 = 0.;
  for(k=1e-5*p.h; k<0.1*p.h; k*= pow(10.,0.05)) {
    get_Delta2k(&p, z, k*pow(10.,-0.025), 1, 1, &delta2);
    sigma2 += delta2 * 0.05 * log(10.);
    printf("%12.5le %12.5le\n", k, sigma2);
  }
  return(0.);
#endif

#if 0
  for(i=NMASS-1;i>=0;i--) printf("%5.2lf %12.5le %9.6lf %12.5le %12.5le\n", 9+0.05*i, abund[i], bias[i], hmpk[i], hmpk[NMASS+i]);
#endif

  /* Make vector that takes the top 4e-5/Mpc^3 haloes, all centrals */
  x = dvector(0,NPARAM-1);

  /* HOD construction */
#ifdef YEONG_SAMPLE1
  printf("Using Yeong sample\n");
  hp[0] = 0.0030*0.10;     /* Number density */
  hp[1] = 0.14;            /* Satellite fraction */
sscanf(argv[1], "%lg", hp+1);
  hp[2] = 11.78;           /* Central mass */
  if (NHOD>3) hp[3] = 1.;  /* Exponent */
  if (NHOD>4) hp[4] = 3.;  /* Min ratio */
#else
  hp[0] = 1.39e-5*0.92;    /* Number density */
  hp[1] = 0.15;            /* Satellite fraction */
  hp[2] = 13.92;           /* Central mass */
  if (NHOD>3) hp[3] = 1.;  /* Exponent */
  if (NHOD>4) hp[4] = 3.;  /* Min ratio */
#endif
  getmodel(hp, x, z, &p);

//  for(i=0;i<NMASS;i++) printf("%5.2lf %12.5le %12.5le %12.5le %12.5le\n", 9+0.05*i, abund[i], x[i], x[i+NMASS], x[i+2*NMASS]);

//  AddClusterFisher(x, F, z, 8.0e7, 3e14, 0.1, 1.00, 0.1, &p);

  /* *.5 in volume to cut to 5000 deg^2 */
//  AddGGLensFisher(F, z, 0.38, 1., 6.0e7*.5, x, 1.43e7*(1.-0.17/0.76), &p);

  /* "LSST style" -- assumes spectroscopic redshifts available */
  AddGGLensFisher(F, z, 0.7, 1., 6.0e7, x, 1.18e7*20., &p);

//  AddNFisher(x, F, 1e8, pow(1.13/0.04,2.), z, &p);
  AddNFisher(x, F, 1e8, 0, z, &p);

  fprintf(stderr, "Fisher computed.\n");

  ProjectFisher(F, Fp, hp, z, &p);

Fp[2][2] += .01;

  gaussjinv(Fp,NHOD);

  for(i=0;i<NHOD;i++) {
    printf("%13.6le %13.6le  ", hp[i], sqrt(Fp[i][i]));
    for(j=0;j<NHOD;j++) printf(" %13.6le", Fp[i][j]/sqrt(Fp[i][i]*Fp[j][j]));
    printf("\n");
  }

  /* Uncertainty on bias */
  AddNFisher(x, bHH, 0, 1, z, &p);
  ProjectFisher(bHH, biasHessian, hp, z, &p);
  varb = 0.;
  for(i=0;i<NHOD;i++) for(j=0;j<NHOD;j++) varb += biasHessian[i][j] * Fp[i][j];
  printf("sigma(bias from HOD fit) = %6.4lf\n", sqrt(varb));

  free_dvector(x,0,NPARAM-1);
  free_dmatrix(F,0,NPARAM-1,0,NPARAM-1);
  free_dmatrix(bHH,0,NPARAM-1,0,NPARAM-1);
  free_dmatrix(Fp,0,NHOD-1,0,NHOD-1);
  free_dmatrix(biasHessian,0,NHOD-1,0,NHOD-1);
  return(0);
}
