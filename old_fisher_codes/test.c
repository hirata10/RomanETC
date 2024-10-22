#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "nr_utils.h"
#include "utils.h"
#include "ps.h"
#include "fom.h"

#include "mnemonics.c"

/* Adds FIRAS Fisher matrix (flag&0x1==0) or overwrites (flag&0x1==1). */
void get_firas_fisher_matrix(double **F, int flag) {
  COSMOPARAM p;
  int ip,jp;

  if (flag&0x1)
    for(ip=0; ip<NPARAMTOT; ip++)
      for(jp=0; jp<NPARAMTOT; jp++)
        F[ip][jp] = 0.;

  ip=0;
  while (&(p.T_cmb) != get_ptr_param(&p,ip)) ip++;
  F[ip][ip] += 1.e6; /* 2.5e5; -- old estimate */
}

/* Adds CMB Fisher matrix.
 *
 * flag&0x1 = clears matrix first
 * flag&0x2 = ON (WMAP9), OFF (Planck)
 */
#define LMAX 2000
#define PLANCK_BANDSX
void get_cmb_fisher_matrix(double **F, COSMOPARAM *p, int flag) {
  int ip,jp,ip2,k1,k2;
  long L, LL;
  char FileName[1024];
  FILE *fp;
  double **dCTT, **dCTE, **dCEE;
  double CTT, CEE, CTE;
  double Lc, fsky, NTT, NEE;
  double **Cov_L, **icov;
  double deriv_w, deriv_h;

  dCTT = dmatrix(0, NPARAMTOT-1, 0, LMAX);
  dCTE = dmatrix(0, NPARAMTOT-1, 0, LMAX);
  dCEE = dmatrix(0, NPARAMTOT-1, 0, LMAX);
  Cov_L = dmatrix(0, LMAX, 0, 2);
  icov = dmatrix(0,2,0,2);

  /* Set noise parameters */
  if (flag&0x2) {
    fsky = 0.7;
    NTT = 2.8e-15 * 4./9.;
    NEE = 5.6e-15 * 4./9.;
    Lc = 620.;
  } else {
    fsky = 0.7;
    NTT = 3.7e-17;
    NEE = 7.5e-17;
    Lc = 1012.;
  }

  /* Clear Fisher matrix if necessary */
  if (flag&0x1)
    for(ip=0; ip<NPARAMTOT; ip++)
      for(jp=0; jp<NPARAMTOT; jp++)
        F[ip][jp] = 0.;

  /* Reads reference file */
  fp = fopen("cmb/cmbout.ref", "r");
  for(L=2; L<=LMAX; L++) {
    fscanf(fp, "%ld %lg %lg %lg", &LL, &CTT, &CEE, &CTE);
#ifdef PLANCK_BANDS
    CTT += 1./(
              2.73e15*exp(-2.99e-6*L*(L+1.)) /* 70GHz */
           + 2.10e16*exp(-1.38e-6*L*(L+1.)) /* 100GHz */
           + 4.84e16*exp(-7.69e-7*L*(L+1.)) /* 143GHz */
           );
    CEE += 1./(
              1.34e15*exp(-2.99e-6*L*(L+1.)) /* 70GHz */
           + 8.19e15*exp(-1.38e-6*L*(L+1.)) /* 100GHz */
           + 1.33e16*exp(-7.69e-7*L*(L+1.)) /* 143GHz */
//         + 4.92e15*exp(-3.81e-7*L*(L+1.)) /* 217GHz */
           )/10000.;

/* Angular Resolution (arcmin) 33 24 14 9.5 7.1 5.0 5.0 5.0 5.0
 * Bandwidth (GHz)  6 8.8 14 33 47 72 116 180 283
 * Average T per pixel 2.0 2.7 4.7 2.5 2.2 4.8 14.7 147 6700
 * Average P per pixel 2.8 3.9 6.7 4.0 4.2 9.8 29.8
 */
#else
    CTT += NTT*exp(L*L/Lc/Lc);
    CEE += NEE*exp(L*L/Lc/Lc);
#endif
#ifdef KILL_LOWL_POL
    if (L<30) CEE *= 1e4;
#endif
    Cov_L[L][0] = CTT;
    Cov_L[L][1] = CEE;
    Cov_L[L][2] = CTE;
  }
  fclose(fp);

  /* Now get parameter derivatives for non-w parameters */
  for(ip=NWPARAM; ip<NPARAMTOT; ip++) {
    sprintf(FileName, "cmb/cmbout.%02d.d", ip-NWPARAM);
    fp = fopen(FileName, "r");
    for(L=2; L<=LMAX; L++) {
      fscanf(fp, "%ld %lg %lg %lg", &LL, dCTT[ip]+L, dCEE[ip]+L, dCTE[ip]+L);
    }
    fclose(fp);
  }

  /* w, Omega_K, onuh2 needs to be made just like h to avoid accidental
   * breaking of degeneracy
   */
  jp = NWPARAM;
  while (&(p->h)!=get_ptr_param(p,jp)) jp++;
  get_derivatives(1100., p, jp, NULL, &deriv_h, NULL);  /* Now jp -> h */
  for(ip=0; ip<NWPARAM; ip++) {
    get_derivatives(1100., p, ip, NULL, &deriv_w, NULL);
    for(L=2; L<=LMAX; L++) {
      dCTT[ip][L] = dCTT[jp][L] * deriv_w/deriv_h;
      dCTE[ip][L] = dCTE[jp][L] * deriv_w/deriv_h;
      dCEE[ip][L] = dCEE[jp][L] * deriv_w/deriv_h;
    }
  }
  ip = NWPARAM;
  while (&(p->Omega_K)!=get_ptr_param(p,ip)) ip++;
  get_derivatives(1100., p, ip, NULL, &deriv_w, NULL);
  for(L=2; L<=LMAX; L++) {
    dCTT[ip][L] = dCTT[jp][L] * deriv_w/deriv_h;
    dCTE[ip][L] = dCTE[jp][L] * deriv_w/deriv_h;
    dCEE[ip][L] = dCEE[jp][L] * deriv_w/deriv_h;
  }

  /* Fix neutrino masses at fixed omh2 to look like simultaneously
   * changing omh2 and also fixing the distance.
   */
  ip = NWPARAM;
  while (&(p->omh2)!=get_ptr_param(p,ip)) ip++;
  get_derivatives(1100., p, ip, NULL, &deriv_w, NULL);
  ip2 = NWPARAM;
  while (&(p->onuh2)!=get_ptr_param(p,ip2)) ip2++;
  for(L=2; L<=LMAX; L++) {
    dCTT[ip2][L] = -dCTT[ip][L] + dCTT[jp][L] * deriv_w/deriv_h;
    dCTE[ip2][L] = -dCTE[ip][L] + dCTE[jp][L] * deriv_w/deriv_h;
    dCEE[ip2][L] = -dCEE[ip][L] + dCEE[jp][L] * deriv_w/deriv_h;
  }

  /* Add contributions to the Fisher matrix */
  for(L=2; L<=LMAX; L++) {
    /* Inverse-cov of TT, EE, TE */
    CTT = Cov_L[L][0];
    CEE = Cov_L[L][1];
    CTE = Cov_L[L][2];
    icov[0][0] = 2*CTT*CTT;
    icov[1][1] = 2*CEE*CEE;
    icov[2][2] = CTT*CEE+CTE*CTE;
    icov[0][1] = icov[1][0] = 2*CTE*CTE;
    icov[0][2] = icov[2][0] = 2*CTT*CTE;
    icov[1][2] = icov[2][1] = 2*CEE*CTE;
    gaussjinv(icov,3);
    for(k1=0; k1<3; k1++)
      for(k2=0; k2<3; k2++)
        icov[k1][k2] *= (2*L+1)*fsky;

    for(ip=0; ip<NPARAMTOT; ip++)
      for(jp=0; jp<NPARAMTOT; jp++) {
        F[ip][jp] += icov[0][0] * dCTT[ip][L] * dCTT[jp][L];
        F[ip][jp] += icov[0][1] * dCTT[ip][L] * dCEE[jp][L];
        F[ip][jp] += icov[0][2] * dCTT[ip][L] * dCTE[jp][L];
        F[ip][jp] += icov[1][0] * dCEE[ip][L] * dCTT[jp][L];
        F[ip][jp] += icov[1][1] * dCEE[ip][L] * dCEE[jp][L];
        F[ip][jp] += icov[1][2] * dCEE[ip][L] * dCTE[jp][L];
        F[ip][jp] += icov[2][0] * dCTE[ip][L] * dCTT[jp][L];
        F[ip][jp] += icov[2][1] * dCTE[ip][L] * dCEE[jp][L];
        F[ip][jp] += icov[2][2] * dCTE[ip][L] * dCTE[jp][L];
      }
  }

  free_dmatrix(dCTT, 0, NPARAMTOT-1, 0, LMAX);
  free_dmatrix(dCTE, 0, NPARAMTOT-1, 0, LMAX);
  free_dmatrix(dCEE, 0, NPARAMTOT-1, 0, LMAX);
  free_dmatrix(Cov_L, 0, LMAX, 0, 2);
  free_dmatrix(icov,0,2,0,2);
}
#undef LMAX

/* Adds SN Fisher matrix
 *
 * flag&0x1 = clears matrix first
 * flag&0x2 = DETF covariance matrix
 * flag&0x4 = implements 0.5% tilt per octave
 * flag&0x8 = kills SN flux uncertainty (for GW sources)
 *
 * file format: first line is Nz, fluxcal, corrz
 * each following line is z, sigma(ln flux)
 */
void get_sn_fisher_matrix(double **F, char FileName[], COSMOPARAM *p, int flag) {
  FILE *fp;
  int ip, jp, iz, jz, Nz;
  double fluxcal, corrz;
  double *z, *sigma;
  double **lnDLcov;
  double **dlnDLdp;

  /* Clear Fisher matrix if necessary */
  if (flag&0x1)
    for(ip=0; ip<NPARAMTOT; ip++)
      for(jp=0; jp<NPARAMTOT; jp++)
        F[ip][jp] = 0.;

  /* Get SN */
  fp = fopen(FileName, "r");
  fscanf(fp, "%d %lg %lg", &Nz, &fluxcal, &corrz);
  z = dvector(0, Nz-1);
  sigma = dvector(0, Nz-1);
  lnDLcov = dmatrix(0, Nz-1, 0, Nz-1);
  dlnDLdp = dmatrix(0, Nz-1, 0, NPARAMTOT-1);
  for(iz=0; iz<Nz; iz++)
    fscanf(fp, "%lg %lg", z+iz, sigma+iz);
  fclose(fp);

  /* Construct covariance matrix for distances */
  for(iz=0; iz<Nz; iz++)
    for(jz=0; jz<Nz; jz++) {
      lnDLcov[iz][jz] = fluxcal*fluxcal*exp(-fabs(z[iz]-z[jz])/corrz);
      if (flag & 0x2) {
        lnDLcov[iz][jz] = 0.0000424*z[iz]*z[jz]*(1.+z[iz]*z[jz]);
      }
      if (iz==jz) lnDLcov[iz][jz] += sigma[iz]*sigma[iz];

      /* Implement tilt per octave if desired */
      if (flag & 0x4) lnDLcov[iz][jz] += 5.203e-5 * log(1.+z[iz]) * log(1.+z[jz]);

      lnDLcov[iz][jz] /= 4.; /* distance errors, not flux */
    }
  gaussjinv(lnDLcov, Nz); /* get inverse-cov */

  /* Build Jacobian from parameters --> lnDL */
  for(iz=0; iz<Nz; iz++)
    for(ip=0; ip<NPARAMTOT; ip++) {
      get_derivatives(z[iz], p, ip, NULL, dlnDLdp[iz]+ip, NULL);
      if (!(flag & 0x8)) {if (&(p->SNflux)==get_ptr_param(p, ip)) dlnDLdp[iz][ip] += 0.5;}
    }

  /* Fisher matrix for parameters is Jac^T * invCov * Jac */
  for(ip=0; ip<NPARAMTOT; ip++)
    for(jp=0; jp<NPARAMTOT; jp++)
      for(iz=0; iz<Nz; iz++)
        for(jz=0; jz<Nz; jz++)
          F[ip][jp] += lnDLcov[iz][jz] * dlnDLdp[iz][ip] * dlnDLdp[jz][jp];

  /* Cleanup */
  free_dmatrix(lnDLcov, 0, Nz-1, 0, Nz-1);
  free_dmatrix(dlnDLdp, 0, Nz-1, 0, NPARAMTOT-1);
  free_dvector(z, 0, Nz-1);
  free_dvector(sigma, 0, Nz-1);
}

/* Constructs BAO Fisher matrix from specified input file.
 *
 * Flags:
 * flag&0x1 causes the Fisher matrix to be cleared first (0 simply adds)
 * flag&0x2 invokes reconstruction
 * flag&0x4 kills correlation of D and H (useful for comparison)
 * flag&0x8 multiplies uncertainties by 100 (i.e. fractional, not percentage)
 * flag&0x10 removes H constraints
 * flag&0x20 reads a final column in the file (z,D,H,junk)
 * flag&0x40 removes D constraints
 *
 * File format:
 * first line gives number of lines in FM
 * other lines are:
 * 1.05 7.80 9.51 3.28 0.492 0.920 0.726 1.736
 * (z, Ngal, Vtot, Veff, D.recon, H.recon, D.base, H.base)
 *
 * Last 4 columns are accuracy in % (1sigma)
 *
 * If BAO_NEW_FORMAT is set, the other lines are
 * (z,D,H).
 */
void get_bao_fisher_matrix(double **F, char FileName[], COSMOPARAM *p, int flag) {
  FILE *fp;
  int i, N, ip, jp;
  double gradD[NPARAMTOT], gradH[NPARAMTOT], dD, dH, ds;
  double z, Ngal, Vtot, Veff, D_recon, H_recon, D_base, H_base, D, H;
  double **covDH;
  double corr = 0.4; /* correlation of D and H measurements */
  double junk;

  /* Allocate D, H covariance matrix */
  covDH = dmatrix(0,1,0,1);

  /* Clear Fisher matrix if necessary */
  if (flag&0x1)
    for(ip=0; ip<NPARAMTOT; ip++)
      for(jp=0; jp<NPARAMTOT; jp++)
        F[ip][jp] = 0.;

  /* Open file, get # redshift slices */
  fp = fopen(FileName, "r");
  fscanf(fp, "%d", &N);

  /* Go over each redshift slice */
  for(i=0; i<N; i++) {
    /* Get redshift and D,H accuracies */
#ifdef BAO_NEW_FORMAT
    if (flag & 0x20) {
      fscanf(fp, "%lg %lg %lg %lg", &z, &D, &H, &junk);
    } else {
      fscanf(fp, "%lg %lg %lg", &z, &D, &H);
   }
#else
    fscanf(fp, "%lg %lg %lg %lg %lg %lg %lg %lg", &z, &Ngal, &Vtot,
      &Veff, &D_recon, &H_recon, &D_base, &H_base);
    D = flag&0x2? D_recon: D_base;
    H = flag&0x2? H_recon: H_base;
#endif

#if 0
    /* Kill some area */
    D /= sqrt(0.75);
    H /= sqrt(0.75);
#endif
    if (flag&0x8) {
      D *= 100.;
      H *= 100.;
    }

    covDH[0][0] = 1e-4*D*D; /* 1e-4 since uncertainty is in percents */
    covDH[1][1] = 1e-4*H*H;
    covDH[0][1] = covDH[1][0] = 1e-4*D*H*corr;
    if (flag & 0x4)
      covDH[0][1] = covDH[1][0] = 0;

// if (i==0) fprintf(stderr, "cov(D,H): %12.5le %12.5le %12.5le\n", covDH[0][0], covDH[0][1], covDH[1][1]);

    /* Removal of Hubble information */
    if (flag & 0x10) {
      covDH[0][1] = covDH[1][0] = 0.;
      covDH[1][1] *= 1e12;
    }

    /* Removal of distance information */
    if (flag & 0x40) {
      covDH[0][1] = covDH[1][0] = 0.;
      covDH[0][0] *= 1e12;
    }

    /* inverse-covariance */
    gaussjinv(covDH,2);

    /* Get sensitivities of D/s and H*s to cosmological parameters */
    for(ip=0; ip<NPARAMTOT; ip++) {
      get_derivatives(z, p, ip, &dH, &dD, &ds);
      gradD[ip] = dD-ds;
      gradH[ip] = dH+ds;
      if (fabs(gradD[ip])<1e-10) gradD[ip] = 0.;
      if (fabs(gradH[ip])<1e-10) gradH[ip] = 0.;
//    if(i==0) fprintf(stderr, "z=%5.2lf param#%2d gradD=%12.5le gradH=%12.5le\n", z, ip, gradD[ip], gradH[ip]);
    }

    /* Add to Fisher matrix */
    for(ip=0; ip<NPARAMTOT; ip++)
      for(jp=0; jp<NPARAMTOT; jp++)
        F[ip][jp] += covDH[0][0]*gradD[ip]*gradD[jp]
                    +covDH[1][1]*gradH[ip]*gradH[jp]
                    +covDH[0][1]*gradD[ip]*gradH[jp]
                    +covDH[1][0]*gradH[ip]*gradD[jp];
  }
  fclose(fp);
  free_dmatrix(covDH,0,1,0,1);
}

/* Constructs BAO Fisher matrix from specified input file.
 *
 * Flags:
 * flag&0x1 causes the Fisher matrix to be cleared first (0 simply adds)
 * flag&0x2 implements BAO systematics floor
 * flag&0x4 uses nP to compute reconstruction (floor at frac_reconstruction)
 *
 * Input file:
 * #redshift fsky sigma_v
 * z dz b ngal frac_reconstruction
 *
 * ngal is in Mpc^-3; frac_reconstruction is between 0 and 1 (reconstruction
 * suppression of NL erasure in amplitude)
 * sigma_v is in units of c.
 *
 * put - sign in front of number of redshifts for Ly-a forest.
 */
void get_bao_fisher_matrix2(double **F, char FileName[], COSMOPARAM *p, int flag) {
  FILE *fp;
  int i, N, ip, jp;
  double z, fsky, dz, V;
  double gradD[NPARAMTOT], gradH[NPARAMTOT], dD, dH, ds;
  double **covDH;
  double b, ngal, b2n, recons;
  double SigT, SigR, SigS, SigZ;
  double P02, Pk, k, mu, sigv, R, integrand, beta, A0;
  double f[2];
  double rcmin, nP;

  double k0, dlnk, s;
  long ik;
  double D2[750];
  double T, var, var1;
  double SNR2, Pkw;
  int is_lya = 0;

  /* Allocate D, H covariance matrix */
  covDH = dmatrix(0,1,0,1);

  /* Clear Fisher matrix if necessary */
  if (flag&0x1)
    for(ip=0; ip<NPARAMTOT; ip++)
      for(jp=0; jp<NPARAMTOT; jp++)
        F[ip][jp] = 0.;

  /* Open file, get # redshift slices */
  fp = fopen(FileName, "r");
  fscanf(fp, "%d %lg %lg", &N, &fsky, &sigv);
  if (N<0) {N=-N; is_lya=1;}

  A0 = get_A0_BAO(p);
  fprintf(stderr, "A0=%6.4lf\n", A0);

  /* Go over each redshift slice */
  for(i=0; i<N; i++) {
    /* Get redshift and D,H accuracies */
    fscanf(fp, "%lg %lg %lg %lg %lg", &z, &dz, &b, &ngal, &recons);
    b2n = b*b*ngal;
    V = 4.*M_PI*fsky * pow(get_DAC(z,p),2.) * dz/get_H(z,p);

    /* Initialize BAO SNR -- for information only */
    SNR2 = 0;

    /* Determine the fractional reconstruction capability based on linear PT */
    if (flag & 0x4) {

      s = getSoundHorizon(p);
      k0 = 1e-5;
      dlnk = 0.01*log(10.);
      get_Delta2k(p,z,k0,dlnk,750,D2);

      /* Integrate linear displacement field */
      var = var1 = 0.;
      for(ik=0; ik<750; ik++) {
        k = k0*exp(ik*dlnk);
        nP = b2n * D2[ik] * 2.*M_PI*M_PI/(k*k*k);
        T = sin(k*s)/(k*s); /* tophat transfer fcn */
        var += T*T*D2[ik] / (k*k);
        var1 += T*T*D2[ik] / (k*k) / (1. + nP);
      }
      rcmin = sqrt(var1/var);

    } else {
      rcmin = 0.;
    }

    /* Seo & Eisenstein method */
    recons = sqrt(rcmin*rcmin*(1-recons*recons) + recons*recons);
    SigT = 14.9 * get_linsigma(p,z,8./p->h) * recons;
    SigR = SigT*(1. + pow(get_Om(p,z),0.6));
    SigS = 0.625*pow(p->obh2,-0.52)*pow(p->omh2,-0.73)/(1.+pow(10.4*p->omh2,-0.95));
#ifdef LONGFINGER
    SigZ = sqrt(sigv*sigv + 0.00097*0.00097) * (1.+z) / get_H(z,p);
#else
    SigZ = sigv * (1.+z) / get_H(z,p);
#endif

    /* For Lya, use Slosar et al model */
    beta = is_lya? 0.336/b-1: pow(get_Om(p,z),0.6)/b;

    get_Delta2k(p,z,k=0.2*p->h,0,1,D2);   
    nP = b2n * D2[0] * 2.*M_PI*M_PI/(k*k*k);
#if 0
    fprintf(stderr, "nP = %6.3lf recons = %5.3lf (%5.3lf min) Sigmas: T=%6.3lf R=%6.3lf S=%6.3lf Z=%6.3lf ", nP, recons, rcmin, SigT, SigR, SigS, SigZ);
#endif

    covDH[0][0] = covDH[0][1] = covDH[1][0] = covDH[1][1] = 0.;
    get_Delta2k(p,z,-0.2*p->h,0,1,&P02);
    P02 *= 2.*M_PI*M_PI/pow(0.2*p->h,3.);

    /* Integration -- Seo & Eisenstein Eq 26 */
    for(k=0.005*p->h; k<0.5*p->h; k+=0.01*p->h) {
      get_Delta2k(p,z,-k,0,1,&Pk);
      Pk *= 2.*M_PI*M_PI/(k*k*k);

      get_Delta2k(p,z,k,0,1,&Pkw);
      Pkw *= 2.*M_PI*M_PI/(k*k*k);

      for(mu=0.025; mu<1; mu+=0.05) {

        R = (1.+beta*mu*mu)*(1.+beta*mu*mu)*exp(-k*k*mu*mu*SigZ*SigZ);
        integrand = k*k*exp(-2.*pow(k*SigS,1.4))*R*R/pow(R*Pk/P02 + 1./(b2n*P02),2.)
                    *exp(-k*k*(1.-mu*mu)*SigT*SigT - k*k*mu*mu*SigR*SigR);

        SNR2 += pow((Pk-Pkw)*R*b2n/(Pk*R*b2n+1.),2.) *exp(-k*k*(1.-mu*mu)*SigT*SigT/recons/recons - k*k*mu*mu*SigR*SigR/recons/recons) / 2.
                  * 0.05 * 0.01*p->h * V * k*k/(2.*M_PI*M_PI);

        f[0] = mu*mu-1;
        f[1] = mu*mu;
        for(ip=0;ip<2;ip++)
          for(jp=0;jp<2;jp++)
            covDH[ip][jp] += f[ip]*f[jp]*integrand;
      }
    }
    for(ip=0;ip<2;ip++)
      for(jp=0;jp<2;jp++)
        covDH[ip][jp] *= 0.05 * 0.01*p->h * A0*A0 * V;

    gaussjinv(covDH,2);
#if 0
    fprintf(stderr, " SNR = %6.2lf\n", sqrt(SNR2));
#endif
    fprintf(stderr, "%5.3lf %12.5le %12.5le %8.6lf %8.6lf %8.5lf %8.6lf\n",
      z, V, b2n*P02, sqrt(covDH[0][0]), sqrt(covDH[1][1]), covDH[0][1]/sqrt(covDH[0][0]*covDH[1][1]),
      sqrt((covDH[0][0]*covDH[1][1]-covDH[0][1]*covDH[0][1])/(covDH[0][0]+2*covDH[0][1]+covDH[1][1])) );

    /* Implement BAO systematics floor */
    if (flag & 0x2) {
      covDH[0][0] += 1e-6 * N;
      covDH[1][1] += 1e-6 * N;
    }
    gaussjinv(covDH,2);

    /* Get sensitivities of D/s and H*s to cosmological parameters */
    for(ip=0; ip<NPARAMTOT; ip++) {
      get_derivatives(z, p, ip, &dH, &dD, &ds);
      gradD[ip] = dD-ds;
      gradH[ip] = dH+ds;
      if (fabs(gradD[ip])<1e-10) gradD[ip] = 0.;
      if (fabs(gradH[ip])<1e-10) gradH[ip] = 0.;
//    if(i==0) fprintf(stderr, "z=%5.2lf param#%2d gradD=%12.5le gradH=%12.5le\n", z, ip, gradD[ip], gradH[ip]);
    }

    /* Add to Fisher matrix */
    for(ip=0; ip<NPARAMTOT; ip++)
      for(jp=0; jp<NPARAMTOT; jp++)
        F[ip][jp] += covDH[0][0]*gradD[ip]*gradD[jp]
                    +covDH[1][1]*gradH[ip]*gradH[jp]
                    +covDH[0][1]*gradD[ip]*gradH[jp]
                    +covDH[1][0]*gradH[ip]*gradD[jp];
  }
  fclose(fp);
  free_dmatrix(covDH,0,1,0,1);
}

/* Adds WL Fisher matrix.
 * Flags:
 * flag& 0x1 causes the Fisher matrix to be cleared first (0 simply adds)
 * flag& 0x2 rejects WL auto-power spectra (kills II)
 * flag& 0x4 rejects GI directions in cross-spectra.
 * flag& 0x8 rejects galaxy auto-power spectra (i.e. LSS info).
 * flag& 0x10 marginalizes over galaxy-mass x-spectrum in each bin (keeps only CCC directions)
 * flag& 0x20 marginalizes over galaxy-ellipticity correlation.
 * flag& 0x40 marginalizes over DETF GI direction.
 * flag& 0x80 marginalizes over galaxy-galaxy x-spectra between different slices.
 * flag&0x100 includes overlapping BAO survey
 * flag&0x200 kills last lensing auto-slice (useful for CMB)
 * flag&0x400 rejects WL cross-power spectra (i.e. all shear shears)
 * flag&0x800 turns *ON* g-g lensing
 *
 * File format:
 * first line is #zbins #lbins first_l delta_ln_l f_sky
 * Remaining lines are:
 * 0.5 1e-9 0.01 0.01 (redshift, gamma^2/n, sigma(z), sigma(cal))
 */
void get_wl_fisher_matrix(double **F, char FileName[], COSMOPARAM *p, int flag) {
  FILE *fp;
  double l0,dlnl,fsky;
  long nl, nz, ip, jp, nzl;
  long start_gg, start_gk, nxtot;
  double *z, *noise, *Cl, **ClCov, *Cl_inc, *Cl_dec, **dCldp, *sigz, *sigcal, *giproj, *marg;
  double dp[] = DEFAULT_VARIATION;
  double ldp[] = LARGE_VARIATION;
  double ref,g1,g2;
  long i,j,iz;
  COSMOPARAM p2;
  long N,il,b1,b2,b3;
  double **Fbig, **Fgg, **dgdp, *aux=NULL;
  int nrot;
  double *eval, **evec;
  double cc, prod, k, l, dr, ngmodes1d, ngmodes3d, dz;
  long ilmax;
  double SN2;

#define KMAX_GGL 0.15
#define RFLOOR_PER_ZBIN 0.001

  /* Clear Fisher matrix if necessary */
  if (flag&0x1)
    for(ip=0; ip<NPARAMTOT; ip++)
      for(jp=0; jp<NPARAMTOT; jp++)
        F[ip][jp] = 0.;

  /* Get experiment information, allocate arrays */
  fp = fopen(FileName, "r");
  fscanf(fp, "%ld %ld %lg %lg %lg", &nz, &nl, &l0, &dlnl, &fsky);
  nzl = nl*NCROSS(nz);
  start_gg = nzl;
  start_gk = start_gg + nzl;
  nxtot = start_gk + nl*nz*nz;
  N = NPARAMTOT + 2*nz; /* total number of parameters including nuisance */

  marg = dvector(0,nxtot-1);
  giproj = dvector(0,nz-1);
  z = dvector(0,nz-1);
  sigz = dvector(0,nz-1);
  sigcal = dvector(0,nz-1);
  noise = dvector(0,nz-1);
  Cl = dvector(0,nxtot-1);
  dCldp = dmatrix(0,nxtot-1,0,N-1);
  dgdp = dmatrix(0,nz-1,0,N-1);
  Fbig = dmatrix(0,N-1,0,N-1);
  Fgg = dmatrix(0,nz-1,0,nz-1);
  Cl_inc = dvector(0,nxtot-1);
  Cl_dec = dvector(0,nxtot-1);
  ClCov = dmatrix(0,nxtot-1,0,nxtot-1);
  for(i=0;i<nz;i++) {
    fscanf(fp, "%lg %lg %lg %lg", z+i, noise+i, sigz+i, sigcal+i);
#if 0
    /* Kill systematics */
    sigz[i] = 1e-5;
    sigcal[i] = 1e-5;
#endif

  }
  fclose(fp);

  get_wl_ps_many(p,nl,l0,dlnl,nz,z,noise,fsky,Cl,ClCov,aux);

#if 1
  for(il=0; il<nl; il++) {
    printf("%6.1lf ", l0*exp(dlnl*il));
    for(b1=0;b1<nz;b1++) {
      b2=b1;
      printf("%10.4le(%10.4le) ", Cl[il+nl*TINDEX(b1,b2)], sqrt(ClCov[il+nl*TINDEX(b1,b2)][il+nl*TINDEX(b1,b2)]));
    }
    printf("\n");
  }
#endif

  /* Inclusion of BAO survey */
  if (flag & 0x100) {
#define GAL_KMAX  100.0
#define GAL_MUMAX 0.5
#define GAL_ZMIN  0.0
#define GAL_ZMAX  2.0

#ifndef FASTDATA
    /* Get WL+gal Fisher */
    if (sym_matrix_invert(ClCov,ClCov,nxtot)) {
      fprintf(stderr, "Failed: C-inverse #1.\n");
      exit(1);
    }
    fprintf(stderr, "Done: C^-1 #1\n");
#endif

    /* For each slice, add 3D galaxy information */
    for(il=0;il<nl;il++)
      for(b1=0; b1<nz; b1++) {
        i = start_gg + il + nl*TINDEX(b1,b1);
        l = l0*exp(dlnl*il);
        k = l/get_DAC(z[b1],p);

#define PERFECTGAL
#ifdef PERFECTGAL
        ClCov[i][i] *= 1e6;
#else
        /* Only use available (k,z) range */
        if (k<GAL_KMAX && z[b1]>GAL_ZMIN && z[b1]<GAL_ZMAX) {
          /* Get radial width and number of modes */
          dz = b1==0? z[1]-z[0]: b1==nz-1? z[nz-1]-z[nz-2]: (z[b1+1]-z[b1-1])/2.;
          if (z[nz-1]>1000 && b1==nz-2) dz = z[nz-2]-z[nz-3];
          dr = dz/get_H(z[b1],p);
          ngmodes1d = k*GAL_MUMAX*dr/M_PI - 1.;
          if (ngmodes1d>0) {
            ngmodes3d = ngmodes1d * 2. * l * l * dlnl * fsky;
            ClCov[i][i] += ngmodes3d/(2.*Cl[i]*Cl[i]);
#if 1
            fprintf(stderr, "lbin #%2ld zbin #%2ld k=%12.5le dr=%12.5le ngmodes: 1D %12.5le 3D %12.5le\n", il, b1, k, dr, ngmodes1d, ngmodes3d);
#endif
          }
        }
#endif
      }
#ifdef FASTDATA
    exit(1);
#endif

    /* Go back to covariance */
    if (sym_matrix_invert(ClCov,ClCov,nxtot)) {
      fprintf(stderr, "Failed: C-inverse #2.\n");
      exit(1);
    }
    fprintf(stderr, "Done: C^-1 #2\n");
  }

#if 0
  fp = fopen("/newman/user2/chirata/fisher-adept/WFIRST-LensPlot2/drm0ps.dat", "w");
  for(il=0; il<nl; il++) {
    fprintf(fp, "%6.1lf ", l0*exp(dlnl*il));
    for(b1=0;b1<nz;b1++) {
      b2=b1;
      fprintf(fp, "%10.4le %10.4le ", Cl[il+nl*TINDEX(b1,b2)], sqrt(ClCov[il+nl*TINDEX(b1,b2)][il+nl*TINDEX(b1,b2)]));
    }
    fprintf(fp, "\n");
  }  
  fclose(fp);
  exit(0);
#endif
#if 0
  fp = fopen("/newman/user2/chirata/fisher-adept/WFIRST-LensPlot2/afta_Cl.dat", "w");
  for(i=0; i<start_gg; i++)
    fprintf(fp, "%19.12le\n", Cl[i]);
  fclose(fp);
  fp = fopen("/newman/user2/chirata/fisher-adept/WFIRST-LensPlot2/afta_Cov.dat", "w");
  for(i=0; i<start_gg; i++)
    for(j=0; j<start_gg; j++)
      fprintf(fp, "%19.12le\n", ClCov[i][j]);
  fclose(fp);
  sym_matrix_invert(ClCov,ClCov,start_gg);
  fp = fopen("/newman/user2/chirata/fisher-adept/WFIRST-LensPlot2/afta_invCov.dat", "w");
  for(i=0; i<start_gg; i++)
    for(j=0; j<start_gg; j++)
      fprintf(fp, "%19.12le\n", ClCov[i][j]);
  fclose(fp);
  exit(0);
#endif

#if 0
  /* Eigenvectors.  NOTE: If you turn this section on, it ruins the Fisher matrix.
   * ... Therefore use only for testing.
   */
  eval = dvector(0,nxtot-1);
  evec = dmatrix(0,nxtot-1,0,nxtot-1);
  nrot = 0;
  jacobi(ClCov, nxtot, eval, evec, &nrot);
  printf("[%ld] %d rotations\n", nxtot, nrot);
  for(i=0; i<nxtot; i++) {
    printf("%10.3le ", eval[i]);
    if (i%12==11)
      printf("\n");
  }
  printf("\n");
  for(i=0; i<nxtot; i++) {
    printf("%10.3le ", evec[i][0]);
    if (i%12==11)
      printf("\n");
  }
  printf("\n");
  free_dvector(eval,0,nxtot-1);
  free_dmatrix(evec,0,nxtot-1,0,nxtot-1);
  exit(1);
#endif

  /* Get derivative of power spectrum wrt cosmological parameters.
   * First copy parameters to p2, then one by one vary them.
   */
  for(ip=0; ip<NPARAMTOT; ip++)
    *get_ptr_param(&p2,ip) = *get_ptr_param(p,ip);
  for(ip=0; ip<NPARAMTOT; ip++)
    if (get_ptr_param(p,ip)!=&(p->tau) && get_ptr_param(p,ip)!=&(p->SNflux)) {
      *get_ptr_param(&p2,ip) += dp[ip];
      get_wl_ps_many(&p2,nl,l0,dlnl,nz,z,noise,fsky,Cl_inc,NULL,aux);
      *get_ptr_param(&p2,ip) -= 2*dp[ip];
      get_wl_ps_many(&p2,nl,l0,dlnl,nz,z,noise,fsky,Cl_dec,NULL,aux);
      *get_ptr_param(&p2,ip) += dp[ip];
      for(i=0;i<nxtot;i++) {
        dCldp[i][ip] = (Cl_inc[i]-Cl_dec[i])/(2.*dp[ip]);
      }

      fprintf(stderr, "Got dPS param #%2ld --> %s\n", ip, dataLabel[ip]);
    }

  /* Derivatives of ps wrt nuisance parameters.
   * NPARAMTOT + iz => shear calibration, amplitude
   * NPARAMTOT + nz + iz => redshift bias
   */
  for(iz=0; iz<nz; iz++)
    for(i=0;i<nxtot;i++)
      dCldp[i][NPARAMTOT+iz] = dCldp[i][NPARAMTOT+nz+iz] = 0.;
  /* Calibration */
  for(il=0;il<nl;il++) {
    for(b1=0; b1<nz; b1++)
      for(b2=0; b2<=b1; b2++) {
        i = il + nl*TINDEX(b1,b2);
        dCldp[i][NPARAMTOT+b1] += Cl[i];
        dCldp[i][NPARAMTOT+b2] += Cl[i];
      }
    for(b1=0; b1<nz; b1++)
      for(b2=0; b2<nz; b2++) {
        i = start_gk + il + nl*(b1+b2*nz);
        dCldp[i][NPARAMTOT+b2] += Cl[i];
      }
  }
  /* Redshift */
  for(iz=0;iz<nz;iz++) {
    z[iz] += 0.01;
    get_wl_ps_many(&p2,nl,l0,dlnl,nz,z,noise,fsky,Cl_inc,NULL,aux);
    z[iz] -= 0.02;
    get_wl_ps_many(&p2,nl,l0,dlnl,nz,z,noise,fsky,Cl_dec,NULL,aux);
    z[iz] += 0.01;
    for(i=0;i<nxtot;i++)
      dCldp[i][NPARAMTOT+nz+iz] = (Cl_inc[i]-Cl_dec[i])/0.02;

    fprintf(stderr, "Got dPS photo-z err #%2ld\n", iz);
  }

  /* Reject convergence auto-powers, if requested */
  if (flag & 0x2)
    for(il=0;il<nl;il++)
      for(b1=0; b1<nz; b1++) {
        i = il + nl*TINDEX(b1,b1);
        ClCov[i][i] += 1e6 * Cl[i]*Cl[i];
      }

  /* Reject GI cross-powers, if requested */
  if (flag & 0x4) {
    for(b1=0; b1<nz-1; b1++) {
      for(b2=b1+1; b2<nz; b2++)
        giproj[b2] = get_lensing_strength(p, z[b1], z[b2])/get_lensing_strength(p, z[b1], z[nz-1]);
      for(il=0;il<nl;il++) {
        ref = Cl[il + nl*TINDEX(b1,b1)];
        ref = 1e6*ref*ref;
#if 1
        ref = 0.003*sqrt(nl+0.0)*sqrt(nz+0.0)*Cl[start_gk + il + nl*(b1+nz*(nz-1))];
        ref = 0.003/sqrt(dlnl)*Cl[start_gk + il + nl*(b1+nz*(nz-1))];

        /* FoMSWG prior */
        ref = il<nl/2? 0: 0.003*sqrt(nl/2.0)*sqrt(nz-1.0)*Cl[start_gk + il + nl*(b1+nz*(nz-1))];
#if 0
        /* Gary prior */
        dz = b1==0? z[1]-z[0]: b1==nz-1? z[nz-1]-z[nz-2]: (z[b1+1]-z[b1-1])/2.;
        if (z[nz-1]>1000 && b1==nz-2) dz = z[nz-2]-z[nz-3];
        ref = 0.007 / sqrt(3.) / sqrt(dlnl) / sqrt(10.*dz/(1.+z[b1])) * Cl[start_gk + il + nl*(b1+nz*(nz-1))] * pow(1.+z[b1],3.);
#endif
        ref = ref*ref;
#endif
        for(b2=b1+1; b2<nz; b2++)
          for(b3=b1+1; b3<nz; b3++)
            ClCov[il+nl*TINDEX(b1,b2)][il+nl*TINDEX(b1,b3)] += ref*giproj[b2]*giproj[b3];
      }
    }
  }

  /* Reject galaxy auto-powers, if requested */
  if (flag & 0x8)
    for(il=0;il<nl;il++)
      for(b1=0; b1<nz; b1++) {
        if (flag & 0x800) {
          l = l0*exp(dlnl*il);
          k = l/get_DAC(z[b1],p);
          if (k<KMAX_GGL) continue;
        }
        i = start_gg + il + nl*TINDEX(b1,b1);
#ifdef GARY_GALPRIOR
        ClCov[i][i] += 4. * pow((il<nl/2? 0.05: 0.1)/1.5 * Cl[i], 2.) / sqrt(0.34);
#else
        ClCov[i][i] += 1e6 * Cl[i]*Cl[i];
#endif
      }

  /* Marginalize over galaxy-mass cross correlations */
  if (flag & 0x10) {
    for(b1=0; b1<nz-1; b1++)
      for(il=0;il<nl;il++) {
        if (flag & 0x800) {
          l = l0*exp(dlnl*il);
          k = l/get_DAC(z[b1],p);
          if (k<KMAX_GGL) continue;
        }
        for(b2=b1+1; b2<nz; b2++)
          for(b3=b1+1; b3<nz; b3++) {
#ifdef GARY_GALPRIOR
            ClCov[start_gk+il+nl*(b1+b2*nz)][start_gk+il+nl*(b1+b3*nz)]
              += (il<nl/2? 0.05*0.05/1.5/1.5 + 0.05*0.05/0.9/0.9: 0.1*0.1/1.5/1.5 + 0.1*0.1/0.6/0.6) 
                *Cl[start_gk+il+nl*(b1+b2*nz)]*Cl[start_gk+il+nl*(b1+b3*nz)]  / sqrt(0.34);
#else
            ClCov[start_gk+il+nl*(b1+b2*nz)][start_gk+il+nl*(b1+b3*nz)]
              += 1e6*Cl[start_gk+il+nl*(b1+b2*nz)]*Cl[start_gk+il+nl*(b1+b3*nz)];
#endif
          }
      }
  }

#ifdef GARY_GALPRIOR
  /* Cross-correlation of galaxy-mass and galaxy-galaxy priors */
  for(b1=0; b1<nz-1; b1++)
    for(il=0;il<nl;il++)
      for(b2=b1+1; b2<nz; b2++) {
        ClCov[start_gg+il+nl*TINDEX(b1,b1)][start_gk+il+nl*(b1+b2*nz)]
          += 2. * (il<nl/2? 0.05*0.05/1.5/1.5: 0.1*0.1/1.5/1.5) * Cl[start_gg+il+nl*TINDEX(b1,b1)] * Cl[start_gk+il+nl*(b1+b2*nz)] / sqrt(0.34);
        ClCov[start_gk+il+nl*(b1+b2*nz)][start_gg+il+nl*TINDEX(b1,b1)]
          += 2. * (il<nl/2? 0.05*0.05/1.5/1.5: 0.1*0.1/1.5/1.5) * Cl[start_gg+il+nl*TINDEX(b1,b1)] * Cl[start_gk+il+nl*(b1+b2*nz)] / sqrt(0.34);
      }
#endif

  /* If requested, turn on galaxy-galaxy lensing.
   */
  if (flag & 0x800)
    for(b1=0; b1<nz-1; b1++) {
      fprintf(stderr, "GGL slice#%2ld\n", b1);

      /* Marginalize bias at each scale, but fixing r.  Need marginalization vector */
      for(il=0;il<nl;il++) {
        /* Direction to marginalize */
        for(i=0;i<nxtot;i++) marg[i] = 0;
        i = start_gg + il + nl*TINDEX(b1,b1);
        marg[i] = 2.*Cl[i];
        for(b2=b1+1;b2<nz;b2++) {
          i = start_gk+il+nl*(b1+b2*nz);
          marg[i] = Cl[i];
        }

        /* Add to covariance */
        for(i=0;i<nxtot;i++)
          if (marg[i]!=0)
            for(j=0;j<nxtot;j++)
              ClCov[i][j] += 1e5*marg[i]*marg[j];
      }

      /* Error floor due to stochasticity error */
      for(i=0;i<nxtot;i++) marg[i] = 0;
      for(il=0;il<nl;il++) {
        for(b2=b1+1;b2<nz;b2++) {
          i = start_gk+il+nl*(b1+b2*nz);
          marg[i] = RFLOOR_PER_ZBIN * Cl[i];
        }
      }
      for(i=0;i<nxtot;i++)
        if (marg[i]!=0)
          for(j=0;j<nxtot;j++)
            ClCov[i][j] += marg[i]*marg[j];

      /* Marginalize galaxy shot noise */
      for(i=0;i<nxtot;i++) marg[i] = 0;
      j = start_gg + nl/2 + nl*TINDEX(b1,b1);
      for(il=0;il<nl;il++) {
        i = start_gg + il + nl*TINDEX(b1,b1);
        marg[i] = Cl[j];  /* This is just to make sure the order of magnitude
                           * of the shot noise is right -- nothing special
                           * about the jth bin.
                           */
      }
      /* Add to covariance */
      for(i=0;i<nxtot;i++)
        if (marg[i]!=0)
          for(j=0;j<nxtot;j++)
            ClCov[i][j] += 1e5*marg[i]*marg[j];
    
    }

  /* Reject galaxy ellipticity correlation */
  if (flag & 0x20)
    for(il=0;il<nl;il++)
      for(b1=0; b1<nz; b1++) {
        i = start_gk + il + nl*(b1+b1*nz);
        ClCov[i][i] *= 1e6;
      }

  /* Marginalize over DETF GI direction */
  if (flag & 0x40) {
    for(b1=0; b1<nz-1; b1++)
      for(il=0;il<nl;il++) {
        for(b2=b1+1; b2<nz; b2++)
          giproj[b2] = Cl[start_gk+il+nl*(b1+b2*nz)]/Cl[start_gg+il+TINDEX(b1,b1)]/(0.5*0.5);
        ref = 1e6 * ClCov[start_gk+il+nl*(b1+b1*nz)][start_gk+il+nl*(b1+b1*nz)];
        for(b2=b1+1; b2<nz; b2++)
          for(b3=b1+1; b3<nz; b3++)
            ClCov[il+nl*TINDEX(b1,b2)][il+nl*TINDEX(b1,b3)] += ref*giproj[b2]*giproj[b3];
        for(b2=b1+1; b2<nz; b2++) {
          ClCov[start_gk+il+nl*(b1+b1*nz)][il+nl*TINDEX(b1,b2)] += ref*giproj[b2];
          ClCov[il+nl*TINDEX(b1,b2)][start_gk+il+nl*(b1+b1*nz)] += ref*giproj[b2];
        }
        ClCov[start_gk+il+nl*(b1+b1*nz)][start_gk+il+nl*(b1+b1*nz)] += ref;
      }
  }

  /* Reject galaxy auto-powers, if requested */
  if (flag & 0x80)
    for(il=0;il<nl;il++)
      for(b1=0; b1<nz-1; b1++)
        for(b2=b1+1; b2<nz; b2++) {
          i = start_gg + il + nl*TINDEX(b1,b2);
          ClCov[i][i] += 1e6 * Cl[start_gg + il + nl*TINDEX(b1,b1)]*Cl[start_gg + il + nl*TINDEX(b2,b2)];
        }

  /* Reject last bin convergence auto-power, if requested */
  if (flag & 0x200)
    for(il=0;il<nl;il++) {
      i = il + nl*TINDEX(nz-1,nz-1);
      ClCov[i][i] += 1e6 * Cl[i]*Cl[i];
    }

  /* Reject convergence cross-powers, if requested */
  if (flag & 0x400)
    for(il=0;il<nl;il++)
      for(b1=0; b1<nz; b1++)
        for(b2=b1; b2<nz; b2++) {
          i = il + nl*TINDEX(b1,b2);
          ClCov[i][i] += 1e6 * Cl[il + nl*TINDEX(b1,b1)]*Cl[il + nl*TINDEX(b2,b2)];
        }

  fprintf(stderr, "Inverting covariance matrix ...\n");

  /* Build WL Fisher matrix */
  if (sym_matrix_invert(ClCov,ClCov,nxtot)) {
    fprintf(stderr, "Failed: C-inverse.\n");
    exit(1);
  }

  fprintf(stderr, "Done: C-inverse\n");

  /* Compute S/N */
  SN2 = 0.;
  for(i=0;i<start_gg;i++)
    for(j=0;j<start_gg;j++)
      SN2 += ClCov[i][j] * Cl[i] * Cl[j];
  fprintf(stderr, "Cosmic shear S/N = %12.5le\n", sqrt(SN2));

  for(jp=0; jp<N; jp++) {
    for(i=0;i<nxtot;i++) {
      prod = 0.;
      for(j=0;j<nxtot;j++)
        prod += ClCov[i][j] * dCldp[j][jp];
      for(ip=0; ip<N; ip++)
          Fbig[ip][jp] += prod * dCldp[i][ip];
    }
  }
  /* z, cal priors */
  for(iz=0;iz<nz;iz++) {
    Fbig[NPARAMTOT+iz][NPARAMTOT+iz] += 1./sigcal[iz]/sigcal[iz];
    Fbig[NPARAMTOT+nz+iz][NPARAMTOT+nz+iz] += 1./sigz[iz]/sigz[iz];
  }

  /* Reduce to NPARAMTOTxNPARAMTOT cosmological Fisher matrix.
   * Increment/decrement of diagonal parameters prevents singularity
   * during inversion process.
   */
  for(ip=0; ip<NPARAMTOT; ip++)
    Fbig[ip][ip] += 1./ldp[ip]/ldp[ip];
  sym_matrix_invert(Fbig,Fbig,N);
  sym_matrix_invert(Fbig,Fbig,NPARAMTOT);
  for(ip=0; ip<NPARAMTOT; ip++)
    Fbig[ip][ip] -= 1./ldp[ip]/ldp[ip];
  for(ip=0; ip<NPARAMTOT; ip++)
    for(jp=0; jp<NPARAMTOT; jp++)
      F[ip][jp] += Fbig[ip][jp];

  free_dvector(marg,0,nxtot-1);
  free_dvector(giproj,0,nz-1);
  free_dvector(z,0,nz-1);
  free_dvector(noise,0,nz-1);
  free_dvector(sigz,0,nz-1);
  free_dvector(sigcal,0,nz-1);
  free_dvector(Cl,0,nxtot-1);
  free_dmatrix(dCldp,0,nxtot-1,0,N-1);
  free_dmatrix(dgdp,0,nxtot-1,0,N-1);
  free_dmatrix(Fbig,0,N-1,0,N-1);
  free_dmatrix(Fgg,0,nz-1,0,nz-1);
  free_dvector(Cl_inc,0,nxtot-1);
  free_dvector(Cl_dec,0,nxtot-1);
  free_dmatrix(ClCov,0,nxtot-1,0,nxtot-1);

  fprintf(stderr, "Done: WL Fisher\n");
#undef KMAX_GGL
#undef RFLOOR_PER_ZBIN
}

/* Adds the ISW effect Fisher matrix.
 * Flags:
 * flag& 0x1 causes the Fisher matrix to be cleared first (0 simply adds)
 *
 * File format:
 * #zbins l0 dlnl nl (first line)
 * z dz b nbar fsky (b=bias; nbar = galaxies/sr)
 */
void get_isw_fisher_matrix(double **F, char FileName[], COSMOPARAM *p, int flag) {
  FILE *fp;
  long ip,jp,iz,nz;
  double z, dz, b, nbar, b2n_3d, fsky;
  double l0, dlnl, l;
  long nl,il;
  double T, sigma, ivar, SN2;
  double *dTdp;
  double dp[] = DEFAULT_VARIATION;
  double T1, T2, junk;
  COSMOPARAM p2;

  dTdp = dvector(0,NPARAMTOT-1);

  /* Clear Fisher matrix if necessary */
  if (flag&0x1)
    for(ip=0; ip<NPARAMTOT; ip++)
      for(jp=0; jp<NPARAMTOT; jp++)
        F[ip][jp] = 0.;

  fp = fopen(FileName, "r");
  fscanf(fp, "%ld %lg %lg %ld", &nz, &l0, &dlnl, &nl);
  for(iz=0; iz<nz; iz++) {
    /* Get basic info for each sample */
    fscanf(fp, "%lg %lg %lg %lg %lg", &z, &dz, &b, &nbar, &fsky);
    b2n_3d = b*b*nbar/get_DAC(z,p)/get_DAC(z,p)/dz*get_H(z,p);

    SN2 = 0.;
    for(il=0;il<nl;il++) {
      /* Get prediction & uncertainty */
      get_isw(p,l0*exp(dlnl*il),dlnl,z,dz,b2n_3d,fsky,&T,&sigma);
      ivar = 1./(sigma*sigma);
      SN2 += T*T*ivar;

      /* Get variations */
      for(ip=0; ip<NPARAMTOT; ip++)
        *get_ptr_param(&p2,ip) = *get_ptr_param(p,ip);
      for(ip=0; ip<NPARAMTOT; ip++) {
        *get_ptr_param(&p2,ip) += dp[ip];
        get_isw(&p2,l0*exp(dlnl*il),dlnl,z,dz,b2n_3d,fsky,&T1,&junk);
        *get_ptr_param(&p2,ip) -= 2*dp[ip];
        get_isw(&p2,l0*exp(dlnl*il),dlnl,z,dz,b2n_3d,fsky,&T2,&junk);
        *get_ptr_param(&p2,ip) += dp[ip];
        dTdp[ip] = (T1-T2)/(2*dp[ip]);
      }

      /* Add to Fisher matrix */
      for(ip=0; ip<NPARAMTOT; ip++)
        for(jp=0; jp<NPARAMTOT; jp++)
          F[ip][jp] += dTdp[ip]*dTdp[jp]*ivar;
    }

    printf("ISW sample %2ld S/N=%8.5lf\n", iz, sqrt(SN2));
  }
  fclose(fp);

  free_dvector(dTdp,0,NPARAMTOT-1);
}

/* Adds the CMB lensing Fisher matrix.
 * Flags:
 * flag& 0x1 causes the Fisher matrix to be cleared first (0 simply adds)
 *
 * File format:
 * #zbins l0 dlnl nl (first line)
 * z dz b nbar fsky (b=bias; nbar = galaxies/sr)
 */
void get_cmblens_fisher_matrix(double **F, char FileName[], COSMOPARAM *p, int flag) {
  FILE *fp;
  long ip,jp,iz,nz;
  double z, dz, b, nbar, b2n_3d, fsky;
  double l0, dlnl, l;
  long nl,il;
  double T, sigma, ivar, SN2;
  double *dTdp;
  double dp[] = DEFAULT_VARIATION;
  double T1, T2, junk;
  COSMOPARAM p2;

  dTdp = dvector(0,NPARAMTOT-1);

  /* Clear Fisher matrix if necessary */
  if (flag&0x1)
    for(ip=0; ip<NPARAMTOT; ip++)
      for(jp=0; jp<NPARAMTOT; jp++)
        F[ip][jp] = 0.;

  fp = fopen(FileName, "r");
  fscanf(fp, "%ld %lg %lg %ld", &nz, &l0, &dlnl, &nl);
  for(iz=0; iz<nz; iz++) {
    /* Get basic info for each sample */
    fscanf(fp, "%lg %lg %lg %lg %lg", &z, &dz, &b, &nbar, &fsky);
    b2n_3d = b*b*nbar/get_DAC(z,p)/get_DAC(z,p)/dz*get_H(z,p);

    SN2 = 0.;
    for(il=0;il<nl;il++) {
      /* Get prediction & uncertainty */
      get_kappa_cmblens(p,l0*exp(dlnl*il),dlnl,z,dz,b2n_3d,fsky,&T,&sigma);
      ivar = 1./(sigma*sigma);
      SN2 += T*T*ivar;

      /* Get variations */
      for(ip=0; ip<NPARAMTOT; ip++)
        *get_ptr_param(&p2,ip) = *get_ptr_param(p,ip);
      for(ip=0; ip<NPARAMTOT; ip++) {
        *get_ptr_param(&p2,ip) += dp[ip];
        get_kappa_cmblens(&p2,l0*exp(dlnl*il),dlnl,z,dz,b2n_3d,fsky,&T1,&junk);
        *get_ptr_param(&p2,ip) -= 2*dp[ip];
        get_kappa_cmblens(&p2,l0*exp(dlnl*il),dlnl,z,dz,b2n_3d,fsky,&T2,&junk);
        *get_ptr_param(&p2,ip) += dp[ip];
        dTdp[ip] = (T1-T2)/(2*dp[ip]);
      }

      /* Add to Fisher matrix */
      for(ip=0; ip<NPARAMTOT; ip++)
        for(jp=0; jp<NPARAMTOT; jp++)
          F[ip][jp] += dTdp[ip]*dTdp[jp]*ivar;
    }

    printf("CMBLENS sample %2ld S/N=%8.5lf\n", iz, sqrt(SN2));
  }
  fclose(fp);

  free_dvector(dTdp,0,NPARAMTOT-1);
}

/* Adds a contribution to the Fisher matrix assuming a measurement at
 * redshift z of the velocity power spectrum, with fractional uncertainty
 * fracsig on the *amplitude*.
 */
void add_zspace_fisher_matrix(double **F, COSMOPARAM *p, double z, double kref, double fracsig) {
  COSMOPARAM p2;
  double k,dk;
  double V2ref, V2p, V2m;
  double dlnVdp[NPARAMTOT];
  int ip,jp;
  double dp[] = DEFAULT_VARIATION;

  /* Get reference velocity */
  k = kref;
  V2ref = get_VPS(p,z,k);

fprintf(stderr, "%9.5lf %9.5lf %12.5le\n", z, k, sqrt(V2ref));

  /* Get variations */
  for(ip=0; ip<NPARAMTOT; ip++)
    *get_ptr_param(&p2,ip) = *get_ptr_param(p,ip);
  for(ip=0; ip<NPARAMTOT; ip++) {
    dk = 0.;
    *get_ptr_param(&p2,ip) += dp[ip];
    V2p = get_VPS(&p2,z,k+dk);
    *get_ptr_param(&p2,ip) -= 2*dp[ip];
    V2m = get_VPS(&p2,z,k-dk);
    *get_ptr_param(&p2,ip) += dp[ip];
    dlnVdp[ip] = (V2p-V2m)/(4.*dp[ip]*V2ref);
  }

  for(ip=0; ip<NPARAMTOT; ip++)
    for(jp=0; jp<NPARAMTOT; jp++)
      F[ip][jp] += dlnVdp[ip] * dlnVdp[jp] / fracsig / fracsig;
}

void add_prior(double **F, char FileName[]) {
  FILE *fp;
  int i,j,m;
  double eigenval;
  double eigenvec[NWPARAM];

  fp = fopen(FileName,"r");
  for(m=0; m<NWPARAM; m++) {
    fscanf(fp, "%lg", &eigenval);
    for(i=0; i<NWPARAM; i++) {
      fscanf(fp, "%lg", eigenvec+i);
    }
    for(i=0; i<NWPARAM; i++)
      for(j=0; j<NWPARAM; j++)
        F[i][j] += eigenvec[i]*eigenvec[j]/eigenval/eigenval;
  }
  fclose(fp);
}

int main(int argc, char **argv) {
  FILE *fp;
  int i,j;
  COSMOPARAM p;
  double z, dH, D, dD, ds, R;
  double **F, **F_small;
  double M,n,k,Pk,gb,gr;
  double fomab,fomdetf;
  double Vec[3000];
  double dphi, r, xi, xinl, D2[12000], D2NL[12000], kr, T;
  double var,kx,ky,kz,dlnk,Del2,Lx,Ly,Lz,J0;
  double nc,nbc,alpha;
  double vd, vdhat,nP,b;
  double integ;
int ik;


  F = dmatrix(0, NPARAMTOT-1, 0, NPARAMTOT-1);
  F_small = dmatrix(0, NWPARAM-1, 0, NWPARAM-1);
  set_default_cosmology(&p);

#if 0
  integ=0;
  z = 19;
  get_Delta2k(&p, z, 1e-4, 0.00115, 12000, D2);
  for(ik=0;ik<12000;ik++) {
    k = 1e-4 * exp(0.00115*ik);
    integ += pow(D2[ik]*2.*M_PI*M_PI/k/k/k, 2)/k * .00115;
    if (ik%500==0) {printf("  %12.5le %12.5le\n", k, D2[ik]*2.*M_PI*M_PI/k/k/k);}
  }
  for(i=0;i<4;i++) integ *= get_H(z,&p)/(1.+z);
  integ *= 2./15./M_PI;
  printf("%12.5le\n", integ);
  printf("%12.5le\n", get_H(z,&p)/(1.+z));
  return(0);
#endif
#if 0
  M = 2.551;
  M = M*M*M;
  R = pow(3.*M/4./M_PI, 1./3.);
  printf("%8.5lf %8.5lf %8.5lf\n", R, get_linsigma(&p, 2.5, R), get_linsigma(&p, 2.5, R)*get_linsigma(&p, 2.5, R));
  return(0);
#endif

#if 0
  z=1.35;
  for(i=0;i<=40;i++) {
    M = pow(10., 10.+i*.125);
    R = pow(M/1.64e8/p.omh2*.141, 1./3.) * .1;
    printf("%8.5lf %9.6lf %9.7lf %9.6lf %9.6lf\n", log(M)/log(10.), R, get_linsigma(&p, z, R),
      get_fM(&p, 4./3.*M_PI*R*R*R, z),
      get_bM(&p, 4./3.*M_PI*R*R*R, z));
  }
  return(0);
#endif
#if 0
  z=2.; b=1.7;
  vd = vdhat = 0.;
  for(i=0;i<=5000;i++) {
    k = pow(10., -5.+0.001*i);
    get_Delta2k(&p, z, k, 0.1, 1, D2);
    nP = 1e-4 * 2.*M_PI*M_PI/(k*k*k)*(*D2)*b*b;
    vd += 0.001*log(10.)*(*D2);
    vdhat += 0.001*log(10.)*(*D2)/(1.+nP);
    if (i%100==0) printf("%11.5E %13.7lE %11.5E %11.5E %8.6lf\n", k, *D2, vd, vdhat, sqrt(vdhat)*b);
  }
  return(0);
#endif

#if 1
  sscanf(argv[1], "%lg", &z);
  get_Delta2k(&p, z, 1e-4, 0.00115, 12000, D2);
  get_Delta2kNL(&p, z, 1e-4, 0.00115, 12000, D2NL);
  for(r=1.5/p.h; r<5.65/p.h; r+=0.1/p.h) {
    xi = xinl = 0;
    for(i=0; i<12000; i++) {
      k = 1e-4 * exp(0.00115*i);
      dphi = 0.00115*k*r/2.;
      xi += D2[i] * sin(k*r)/(k*r) * sin(dphi)/dphi;
      xinl += D2NL[i] * sin(k*r)/(k*r) * sin(dphi)/dphi;
    }
    xi *= 0.00115;
    xinl *= 0.00115;
    printf("%12.5le %12.5le %12.5le\n", r*p.h, xi, xinl);
  }
  return(0);
#endif

  get_firas_fisher_matrix(F, 0x1);

// get_cmb_fisher_matrix(F, &p, 0x2);
//  F[NWPARAM+8][NWPARAM+8] += pow(0.015,-2.);

//   get_cmblens_fisher_matrix(F, "data/isw.s4", &p, 0x0);
//   get_isw_fisher_matrix(F, "data/isw.s4", &p, 0x0);
//  get_wl_fisher_matrix(F, "data/wl.JDEM_MIN", &p, 0x13c);
//  get_wl_fisher_matrix(F, argv[1], &p, 0x3c); // should be 0x3c
  get_wl_fisher_matrix(F, argv[1], &p, 0x18); // stat only

//  get_bao_fisher_matrix(F, "dan/bao.np2.z1.5_2", &p, 0x20);
//  get_bao_fisher_matrix(F, "dan/fisher.11_20.zodi2.npA.dat", &p, 0x60);
//  get_bao_fisher_matrix(F, argv[1], &p, 0x8);
//  get_bao_fisher_matrix(F, "dan/karl.110608.dat", &p, 0x0);

#if 0
//  get_bao_fisher_matrix2(F, "BAOREV/bao.s3_1", &p, 0x6);
//  get_bao_fisher_matrix2(F, "BAOREV/bao.s3_2", &p, 0x6);
//  get_bao_fisher_matrix2(F, "BAOREV/bao.s3_3", &p, 0x6);
//  get_bao_fisher_matrix2(F, "BAOREV/bao.s3_4", &p, 0x6);
//  get_bao_fisher_matrix2(F, "BAOREV/bao.s3_5", &p, 0x6);
//  get_bao_fisher_matrix2(F, "BAOREV/bao.s3_6", &p, 0x6);
#endif
//  get_bao_fisher_matrix2(F, "BAOREV/bao.PFS", &p, 0x6);
//  get_bao_fisher_matrix2(F, "EUCLID/bao.BigBoss", &p, 0x6);
//  get_bao_fisher_matrix2(F, "BAOREV/bao.Euclid.11k", &p, 0x6);
//  get_bao_fisher_matrix2(F, argv[1], &p, 0x4);

//  get_bao_fisher_matrix2(F, "BAOREV/bao.ultimate", &p, 0x6);

// add_zspace_fisher_matrix(F, &p, 0.10, 0.035, 0.16/sqrt(2.) ); /* 2dF   Verde et al x2 */

//  write_swg_fisher(F, &p, "matrices/F_BAO_Karl_1e4_1000.swg");
//  write_swg_fisher(F, &p, "anze/bao.all.24k.swg");

//  add_fisher_yw(F,&p);
  write_swg_fisher(F, &p,argv[2]);
  return(0);

  for(i=0; i<NPARAMTOT; i++)
    F[i][i] += 1e-4;

  /* OTHER PRIORS HERE */
  F[NWPARAM+5][NWPARAM+5] += 1e12; /* alphas */
// F[NWPARAM+6][NWPARAM+6] += 1./0.03/0.03/p.h/p.h; /* h */
// F[NWPARAM+7][NWPARAM+7] += 1e12; /* Omega_K */
  F[NWPARAM+13][NWPARAM+13] += 1e12; /* Omega_nu h^2 */

  /* Remove w parameters */
#if 0
  for(i=0; i<NWPARAM; i++)
  for(j=0; j<NWPARAM; j++)
    F[i][j] += 1e8;
#endif
#if 0
  for(i=0; i<NWPARAM; i++)
    F[i][i] += 1e12;
#endif

  /* Fix growth function to GR values */
#if 1
  F[NWPARAM+10][NWPARAM+10] += 1e12;
  F[NWPARAM+11][NWPARAM+11] += 1e12;
  F[NWPARAM+12][NWPARAM+12] += 1e12; /* growth curvature off */
#endif

#ifdef DAN_SPLINE
  F[0][0] += 1e12;
  F[8][8] += 1e12;
#endif

  sym_matrix_invert(F, F, NPARAMTOT);

#if 1
  for(i=0; i<NPARAMTOT; i++)
    printf("%2d: %s: %12.5le\n", i, dataLabel[i], sqrt(F[i][i]));
  printf("\n");
  for(i=0; i<NPARAMTOT; i++) {
    for(j=0; j<NPARAMTOT; j++)
      printf(" %5d", (int)floor(0.5+1000*F[i][j]/sqrt(F[i][i]*F[j][j])));
    printf("\n");
  }
  printf("\n");
#endif

#ifdef DAN_SPLINE
  return(0);
#endif

#if 1
  /* Compute FOM */
  submatrix(F, F_small, NWPARAM);
  fomab = get_fom_ab(F_small, NWPARAM);
  printf("FOM/AB = %12.5le\n", fomab);
  fomdetf = get_fom_detf(F_small, NWPARAM);
  printf("FOM/DETF = %12.5le\n", fomdetf);
#endif

#if 0
  /* Growth function derivative FOM */
  printf("sigma(Slope)=%12.5le sigma(Curvature)=%12.5le corr=%12.5le\n", sqrt(F[20][20]), sqrt(F[21][21]),
    F[20][21]/sqrt(F[20][20]*F[21][21]));
  printf("pivot a=%12.5le sigma(Slope_pivot)=%12.5le\n", exp(-F[21][20]/F[21][21]),
    sqrt(F[20][20]-F[21][20]*F[21][20]/F[21][21]));
  return(0);
#endif

#if 0
  /* Growth function FOM */
  printf("sigma(Amp)=%12.5le sigma(Slope)=%12.5le corr=%12.5le\n", sqrt(F[19][19]), sqrt(F[20][20]),
    F[20][19]/sqrt(F[19][19]*F[20][20]));
  printf("pivot a=%12.5le sigma(Amp_pivot)=%12.5le\n", exp(-F[20][19]/F[20][20]),
    sqrt(F[19][19]-F[20][19]*F[20][19]/F[20][20]));
#endif

  return(0);
}
