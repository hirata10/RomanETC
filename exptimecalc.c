/* BAO+WL Exposure Time Calculator.
 * v19 by C. Hirata, March 30, 2022
 *
 * Compilation instructions are:
 * gcc exptimecalc.c -o pzcaletc.exe -lm -Wall -O3 -DPZCAL_MODE
 * gcc exptimecalc.c -o baoetc.exe -lm -Wall -O3 -DBAO_MODE
 * gcc exptimecalc.c -o wletc.exe -lm -Wall -O3 -DWL_MODE
 * gcc exptimecalc.c -o spcontetc.exe -lm -Wall -O3 -DSPCONT_MODE
 *
 * Additional options:
 * -DKEEP_WINDOW_OPEN = requires user to hit <return> before program exits
 * -DTILT_BACK = tilts back from Sun (up to 115 degrees)
 * -DIN_DLON = requires input for longitude offset relative to Sun
 * -DREAD5 = sets read noise at 5e and dark current at 0.05 e/s/pix (instead of taking user input)
 * -DWLCUT_DEFAULT = sets default WL cuts (Res > 0.4 and sigma_e < 0.2)
 * -DOUT_WL_CAT = dumps weak lensing source catalog to a file
 * -DLOGICAL_READ_FLOOR = changes the behavior of the read noise floor from quadrature to logical
 * -DWFE_OVERRIDE = overrides control on WFE in WL mode (useful to extract detection sensitivity
 *    estimates from the imager blueward of 14*RMSWFE)
 * -DWFE_HOPKINS_EXP = put wavefront error in exponential Hopkins distribution
 * -DNO_OBJ_CONTINUUM = ignores source continuum Poisson fluctuations as a noise source
 * -DUSE_SNRE = use the expected rather than observed SNR for the BAO cuts
 * -DOUT_EXTRA_PSF_PROPERTIES = write additional properties of the PSF (for imaging modes)
 * -DOIII_GAL = predict [O III] properties instead of H alpha
 * -DOII_GAL = predict [O II] properties instead of H alpha
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/* The photo-z calibration mode uses many of the BAO mode functions,
 * with a few exceptions, so go ahead and set the BAO_MODE flag here.
 */
#ifdef PZCAL_MODE
#define BAO_MODE
#endif

/* This is useful to avoid repeated blocks of code -- defines BAO_OR_SPCONT_MODE if
 * we are in either the BAO or SPCONT modes.
 */
#ifdef BAO_MODE
#define BAO_OR_SPCONT_MODE
#endif
#ifdef SPCONT_MODE
#define BAO_OR_SPCONT_MODE
#endif

/* Default case is to input a wave front error */
#define WFE_INPUT

/* Thermal background margin factor -- should be >=1 */
#define THERMAL_MARGIN_FACTOR 1.00

/* Minimum S/N per band for a shape measurement */
#define WL_SNR_MIN 18.0

/* --- CONSTANTS ASSOCIATED WITH OUTPUT THAT WE DO NOT REPEAT FOR EACH USER INPUT --- */

/* Reference galaxy size for BAO outputs */
#define R12_REF 0.3

/* Spacing of output redshifts */
#define DZ_OUT 0.05

/* Detector parameters -- NEW in v9.
 *
 * Since these are not anticipated to change during a run of the ETC, we will
 * store them in a block here (analogous to a Fortran common block).
 */
/* == Applicable to all detectors == */
double PIXSIZE_PHYS;    /* Physical pixel size in microns                                */
double PIX_CD;          /* Charge diffusion length (microns) - 1 sigma, assumed Gaussian */
double DET_IS_CCD;      /* Detector is a CCD? (Yes=1/No=0)                               */
/* == Applicable to NIR detectors only == */
double T_FRAME;         /* Time to read out each frame, in seconds                       */
double VAR_READ_FRAME;  /* Read noise variance per frame (e-2) -- roughly cdsnoise^2/2   */

/* H alpha rest wavelength in vacuum (microns) */
#define LAMBDA0_HA 0.6565
/* [O III] rest wavelength in vacuum (microns) */
#define LAMBDA0_OIII 0.5009
/* [O III] doublet intensity ratio */
#define OIII_RATIO 3.0
/* [O II] rest wavelength in vacuum (microns) */
#define LAMBDA0_OII 0.3728

/* Cosmological parameters (this code assumes Lambda CDM).
 * H0 in km/s/Mpc
 */
#define COSMO_H0 71.903
#define COSMO_OMEGA_M 0.25648

/* Ratio of half-light to scale radius for exponential profile galaxy */
#define RAT_HL_SL_EXP 1.67834

/* Conversions */
#define METERS_PER_MPC 3.08567802e22
#define SQDEG_PER_SR 3.282806350011744293e3

/* Setp sizes for galaxy profile integrals: need more accuracy for WL than BAO */
#ifdef BAO_OR_SPCONT_MODE
#define DELTA_STEP 0.0025
#endif
#ifdef WL_MODE
#define DELTA_STEP 0.0005
#endif

/* Smearing of galaxy size due to spectral dispersion -- sigma on 1 axis in arcsec */
#define SIGMA_PROF 0.040

#ifdef NO_OBJ_CONTINUUM
#define EWD_REF 1e12
#else
/* Reference equivalent width in arcsec -- 100A EW rest frame @ D_theta = 240 arcsec. */
#define EWD_REF 3.656
#endif

/* PSF attributes structure */
typedef struct {
  double pixscale;   /* pixel scale in arcsec */
  double sigma;      /* Gaussian smearing, in arcsec */
  double lD;         /* central lambda/D * 206265 */
  double centobs;    /* central obscuration as a fraction of primary diameter */
  int is_broadband;  /* 0=monochromatic, 1=broadband */
  double dll;        /* fractional width of band, i.e. ln(lambda) = central +/- 0.5*dll;
                      * only used if broadband */
  double rmswfe;     /* rms wave front error in waves */
  double linDisp;    /* linear dispersion (arcsec/um): only used for BAO mode */
}
PSF_DATA;

/* Throughput table structure.
 * Definition of throughput curve is via linear interpolation, constants if
 * off the end of the curve.
 */
#define N_THROUGHPUT_MAX 512
typedef struct {
  int N;                               /* Number of points in interpolation grid */
  double lambda[N_THROUGHPUT_MAX];     /* Wavelengths (in microns) of nodes */
  double throughput[N_THROUGHPUT_MAX]; /* Throughput values */
}
THROUGHPUT_DATA;

/* Parameters of the thermal model for a telescope.
 */
typedef struct {
  double T_tel;                /* Telescope temperature [K] -- PM/SM + assoc struct */
  double pmsm_emissivity;      /* Emissivity of primary & secondary mirror surfaces (epsilon) */
  double outer_ann_ratio;      /* Ratio of outer diameter of beam to that of primary (rho) */
  int ismask;                  /* Cold pupil mask? 1 = yes, 0 = no */
  double T_aft;                /* Temperature of aft optics [K] */
  double aft_net_emissivity;   /* Net emissivity of aft optics (total of all surfaces in 'series') */
  double post_filter_rate;     /* Count rate in e-/pix/s from structures downstream of the filter */
}
THERMAL_DATA;

/* If applicable, wavefront error table, for use in the spectroscopy mode.
 */
typedef struct {
  double lambdabest;           /* Wavelength [microns] of best wavefront */
  double rmswfebest;           /* Best RMS wavefront error [microns] */
  double slope1, slope2;       /* d(RMS WFE)/d(lambda) [dimensionless] at blue (slope1) or red (slope2) wavelengths
                                * -- normally slope1<0 and slope2>0.
                                */
}
WAVEFRONT_ERROR_FUNCTION;

/* --- INTERFACE FUNCTIONS --- */

/* Function to wait for user -- useful on some platforms */
void wait_for_user(void) {
#ifdef KEEP_WINDOW_OPEN
  char temp;
  printf("\nEnter <return> to continue --> ");
  scanf("%c", &temp);
  scanf("%c", &temp);
#endif
}

/* --- SPECIAL FUNCTIONS --- */

/* Bessel functions, J0 and J1:
 * These are needed to convert between real and Fourier space.
 * The |x|>3 expansions are formulae 9.4.3,9.4.6 of Abramowitz &
 * Stegun, good to several parts in 10^8 (i.e. good enough for
 * PSF prediction work). At |x|<3 the integral definition
 * is used.
 */

double getJ0(double x) {
  double f,t,u;
  int i;
  
  x = fabs(x);
  
  /* Large value - Abramowitz & Stegun */
  if (x>3) {
    u=3./x;
    f = 0.79788456 - 0.00000077*u - 0.00552740*u*u - 0.00009512*u*u*u
       + u*u*u*u*(0.00137234-0.00072805*u+0.00014476*u*u);
    t = x - 0.78539816 - 0.04166397*u - 0.00003954*u*u + 0.00262573*u*u*u
       + u*u*u*u*(-0.00054125 -0.00029333*u +0.00013558*u*u);  
    return(f/sqrt(x)*cos(t));
  }
    
  /* Small value - Abramowitz & Stegun */
  f = 0.;
  for(i=0;i<50;i++) f+=cos(x*cos((i+0.5)*M_PI/50.));
  return(0.02*f);
}

double getJ1(double x) {
  double f,t,u,s;
  int i;
  
  /* Flip sign for negative values */
  s=1; if (x<0) {x=-x; s=-1;} 
  
  /* Large value - Abramowitz & Stegun */
  if (x>3) {
    u=3./x;
    f = 0.79788456 + 0.00000156*u + 0.01659667*u*u + 0.00017105*u*u*u
       + u*u*u*u*(-0.00249511 + 0.00113653*u - 0.00020033*u*u);
    t = x - 2.35619449 + 0.12499612*u + 0.0000565*u*u - 0.00637879*u*u*u
       + u*u*u*u*(0.00074348 + 0.00079824*u - 0.00029166*u*u);
    return(s*f/sqrt(x)*cos(t));
  }
  
  /* Small value - Abramowitz & Stegun */
  f = 0.;
  for(i=0;i<50;i++) f+=cos(x*sin((i+0.5)*M_PI/50.)-((i+0.5)*M_PI/50.));
  return(0.02*s*f);  
}

/* Error function */
double geterf(double x) {

  int j;
  double s, term, erfabsx, erfcx, y, dy, u;

  s=1;
  if (x<0) {
    s=-1; x=-x;
  }

  /* For values greater than 6, the result for erf will be unity to machine precision */
  if (x>6) return(s);

  /* Taylor expansion for smallest values */
  if (x<1.5) {
    erfabsx = 0.;
    term = 2./sqrt(M_PI)*x;
    for(j=0;j<=24;j++) {
      erfabsx += term/(2*j+1);
      term *= -x*x/(j+1);
    }
    return(s*erfabsx);
  }

  /* Compute erfc(x) by transforming the complementary error function integral:
   * erfc(x) = 2/sqrtpi * int_x^infty e^(-z^2) dz
   * transform: z = x + exp(y)/x, dz = exp(y)/x * dy
   * erfc(x) = 2/(sqrtpi*x) * int_-infty^infty e^{y-[x+exp(y)/x]^2} dx
   * The integrand rises rapidly, peaks near y~-0.7, and then falls.
   * Convergence is exponential or faster in both directions.
   * Current integration range -30<y<3 chosen for accuracy at 1.5<x<6.
   */
  erfcx = 0.;
  dy = 0.01;
  for(j=-3000;j<300;j++) {
    y=j*dy;
    u = x+exp(y)/x;
    erfcx += exp(y-u*u);
  }
  erfcx *= 2./sqrt(M_PI)/x*dy;
  erfabsx=1-erfcx;

  return(s*erfabsx);
}

/* --- ROUTINES TO SET UP DETECTOR PARAMETERS --- */

/* Configures detector parameters, given a type of detector.
 * (NEW in v9.)
 */
void ConfigDetector(int det_type) {
  int is_ccd = 0;

  switch(det_type) {
    case 0: /* H2RG - 18 um, 32 channel readout */
      PIXSIZE_PHYS       = 18.0;
      PIX_CD             =  2.94;
         /* This value is the variance from Barron et al 2007 PASP 119, 466.
          * They fit a sech-profile, but the MTF calculator here assumes a Gaussian
          * so matching the second moment is most appropriate. The second moment of
          * P(x) propto sech (x/l) distribution is (pi*l/2)^2.
          */
      T_FRAME            = 1.3;
      VAR_READ_FRAME     = 200.;
      break;

    case 1: /* e2v CCD (Euclid) */
      PIXSIZE_PHYS       = 12.0;
      PIX_CD             =  5.0;
      is_ccd = 1;
      break;

    case 2: /* H4RG - 10 um, 32 channel readout */
      PIXSIZE_PHYS       = 10.0;
      PIX_CD             =  2.94;  /* Assumed unchanged from H2RG until better data available */
      T_FRAME            = 2.75;  /* 200 kHz readout */
      VAR_READ_FRAME     = 200.;
#ifdef RN_EIGHTEEN
      VAR_READ_FRAME     = 162.;
#endif
      break;

    case 3: /* H4RG - 10 um, 64 channel readout */
      /* Note -- equivalent to assuming unchanged read noise for 32 channel readout @ 200 kHz */
      PIXSIZE_PHYS       = 10.0;
      PIX_CD             =  2.94;  /* Assumed unchanged from H2RG until better data available */
      T_FRAME            = 2.62;   /* 100 kHz readout */
      VAR_READ_FRAME     = 200.;
      break;

    case 4: /* a 'really good' H1RG (for the IFU) - 18 um, 16 channel readout */
      PIXSIZE_PHYS       = 18.0;
      PIX_CD             =  2.94;
      T_FRAME            = 0.66;
      VAR_READ_FRAME     = 112.5;
      break;

    case 5: /* a 'really good' H1RG (for the IFU) - 18 um, 16 channel readout, 2x2 binned in post-processing */
      PIXSIZE_PHYS       = 36.0;
      PIX_CD             =  2.94;
      T_FRAME            = 0.66;
      VAR_READ_FRAME     = 450.0;
      break;

    case 6: /* H4RG - 10 um, 32 channel readout, 15 e CDS */
      PIXSIZE_PHYS       = 10.0;
      PIX_CD             =  2.94;  /* Assumed unchanged from H2RG until better data available */
      T_FRAME            = 3.04;  /* 200 kHz readout */
#ifdef BAO_MODE
      T_FRAME            = 4.03;  /* with spectroscopic overheads */
#endif
      VAR_READ_FRAME     = 112.5;
      break;

    default:
      fprintf(stderr, "Error: ConfigDetector: Unrecognized detector type = %d\n", det_type);
      exit(1);
      break;
  }

  DET_IS_CCD=is_ccd; /* Permanently set the detector as a CCD or not */
}

/* --- ROUTINES TO COMPUTE THE THROUGHPUT --- */

/* Compute the throughput from a table and a wavelength.
 * Definition: If given T=NULL, returns 1.
 */
double get_throughput(double lambda, THROUGHPUT_DATA *T) {
  int i, N;

  if (T==NULL) return(1);

  N = T->N;

  /* If off the end of the table, give endpoint values.
   * If N=1 this automatically exits.
   */
  if (lambda<=T->lambda[0]) return(T->throughput[0]);
  if (lambda>=T->lambda[N-1]) return(T->throughput[N-1]);

  /* Interpolation */
  i=0;
  while (lambda>T->lambda[i+1] && i<N-2) i++;
  return(T->throughput[i] + (T->throughput[i+1]-T->throughput[i])*(lambda-T->lambda[i])/(T->lambda[i+1]-T->lambda[i]));
}

/* --- ROUTINES TO COMPUTE THE MODULATION TRANSFER FUNCTION (MTF) --- */

/* Compute the wavefront error from a wavefront error function.
 * The wavefront error function currently is a broken linear function.
 */
double get_rmswfe_from_fcn(double lambda, WAVEFRONT_ERROR_FUNCTION *wavefront) {
  double x;
  x = lambda - wavefront->lambdabest;
  return(wavefront->rmswfebest + (x>0?wavefront->slope2:wavefront->slope1)*x);
}

/* Monochromatic MTF at center of filter is computed for diffraction from
 * an unaberrated annular aperture, plus pixel tophat and Gaussian smoothing
 * effects. The input Fourier mode is (u,v) cycles per arcsec. Output is
 * between -1 and +1 inclusive.
 *
 * This routine now has an approximate method to account for low order
 * aberrations.
 */
double get_MTF_mono(double u, double v, PSF_DATA *psf) {
  double upabs, vpabs, kxabs, kyabs, k1, k2;
  double gt_out, gt_in, gt_x, Dc, Ac;
  double T_pix, T_sig, T_diff, T_wfe, T_wfe2;
  double cosbeta, beta, area_triangle, alp1, alp2;

  /* --- Part 1: Computation of the pixel tophat (T_pix) --- */

  /* Compute |u| and |v| in cycles per pixel (upabs,vpabs) and in
   * radians per half-pixel (kxabs,kyabs).
   */
  upabs = fabs(u*psf->pixscale);
  vpabs = fabs(v*psf->pixscale);
  kxabs = M_PI*upabs;
  kyabs = M_PI*vpabs;
  
  /* Pixel transfer: a product of two sincs */
  T_pix = (kxabs<1e-6? 1.-kxabs*kxabs/6.: sin(kxabs)/kxabs)
         *(kyabs<1e-6? 1.-kyabs*kyabs/6.: sin(kyabs)/kyabs);

  /* --- Part 2: Computation of the Gaussian smoothing (T_sig) --- */

  T_sig = exp(-2.*M_PI*M_PI*(u*u+v*v)*psf->sigma*psf->sigma);

  /* --- Part 3: Computation of the diffraction MTF (T_diff) --- */

  /* First get linear and area central obscurations.
   * The central obscuration, if smaller than 1e-14, is set to 1e-14 to
   * avoid division by zero errors.
   */
  Dc = psf->centobs;
  if (Dc<1e-14) Dc=1e-14;
  Ac = psf->centobs*psf->centobs;

  /* Wave number as a function of the cutoff at outer diameter */
  k1 = sqrt(u*u+v*v)*psf->lD;
 
  /* The method is that we want the fraction of the area of an annulus of
   * inner radius Dc and outer radius 1 that overlaps with a version of
   * itself that is translated by 2*k1. This can be thought of as having
   * three parts: we can consider the overlap of the outer circles,
   * that of the inner circles, and cross terms (overlap of one outer
   * and one inner circle).
   *
   * The area of the outer circle is pi, that of the inner is pi*Ac,
   * and that of the annulus is pi*(1-Ac).
   */

  /* Area of outer circle that overlaps itself after translation
   * by 2*k1.
   */
  gt_out = k1<1? 2.*( acos(k1) - k1*sqrt(1-k1*k1)): 0;
  
  /* Area of inner circle that overlaps itself after translation by 2*k1.
   */
  k2 = k1/Dc;
  gt_in = Ac * ( k2<1? 2.*( acos(k2) - k2*sqrt(1-k2*k2)): 0 );
  
  /* Cross-term between inner and outer circles.
   * Need area of sector bounded by circles of radii 1 (center: O) and Dc
   * (center: O') with separation O-O' equal to 2*k1.
   *
   * First consider that the two circles have centers O and O' and intersect
   * at two points P and P'. Then we can find the angle beta = OPO' (unless
   * the inner circle is entirely inside the outer, or both circles are outside
   * each other entirely, in which case P and P' do not exist).
   */
  cosbeta = (-4*k1*k1 + 1 + Ac) / (2.*Dc);
  if (cosbeta>=1) {  
    /* The inner circle is inside the outer circle: overlap area is pi*Ac. */
    gt_x = M_PI*Ac;
  } else
  if (cosbeta<=-1) {
    /* No overlap */
    gt_x = 0.;
  } else {
    /* Partial overlap: area by sum of parts of circles minus triangles */
    beta = acos(cosbeta);                           /* angle beta = OPO' */
    area_triangle = 0.5*Dc*sin(beta);               /* area of OPO' */
    alp1 = asin(area_triangle / k1);                /* angle POO' */
    alp2 = M_PI - alp1 - beta;                      /* angle PO'O */
    gt_x = alp1 + alp2*Ac - area_triangle * 2.;
  }

  /* Diffraction transfer is sum of the sectors, divided by area of
   * annulus. Note cross-term has a - sign.
   */
  T_diff = (gt_out + gt_in - 2.*gt_x)/M_PI/(1.-Ac);

  /* The effect of each type of aberration on the MTF is different. However one
   * can use the following approximate method, which has been tested against
   * the low-order aberrations (defocus, astigmatism, coma, and trefoil)
   * and is conservative in these cases.
   */
  if (psf->rmswfe>0.225) {
    fprintf(stderr, "Error: RMS WFE out of range of validity of approximations.\n");
    wait_for_user();
    exit(1);
  }

  /* Aberration contribution comes from the default, which is worst case of 2nd/3rd order
   * aberrations; or the exponential Hopkins distribution, a Gaussian with assumed
   * correlation length.
   */
  T_wfe2 = 1.; /* this line is to prevent an unused variable warning; it has no effect */
#ifdef WFE_HOPKINS_EXP
  T_wfe = exp(-4.*M_PI*M_PI*psf->rmswfe*psf->rmswfe*(1.-exp(-18.*k1*k1)));
#else
  T_wfe = exp(-72.*psf->rmswfe*psf->rmswfe*(1.-0.5*exp(-32.*k1*k1)-0.5*exp(-8.*k1*k1)));
  T_wfe2 = 0.5*exp(-720.*psf->rmswfe*psf->rmswfe*k1*k1) + 0.5*exp(-2160.*psf->rmswfe*psf->rmswfe*k1*k1);
  if (psf->rmswfe>0.175) {
    T_wfe = T_wfe2;
  } else if (psf->rmswfe>0.125) {
    T_wfe = T_wfe + (psf->rmswfe-0.125)/0.05*(T_wfe2-T_wfe);
  }
#endif
  T_diff *= T_wfe;

#ifdef NIIFRAC
  /* Ascribes a fraction NIIFRAC of the flux to [NII] instead of H alpha.
   * Treats as a smearing of the effective PSF; assumes 240 arcsec dispersion.
   */
  T_sig *= 1.-NIIFRAC + NIIFRAC*0.247*cos(2.*M_PI*u*0.538) + NIIFRAC*0.753*cos(2.*M_PI*u*0.756);
#endif

  /* The three effects are convolved in real space so multiplied in Fourier space */
  return(T_pix*T_sig*T_diff);
}

/* Full MTF is computed for diffraction from an unaberrated annular aperture,
 * plus pixel tophat and Gaussian smoothing effects. The input Fourier
 * mode is (u,v) cycles per arcsec. Output is between -1 and +1 inclusive.
 * For monochromatic PSFs, use get_MTF_mono; for general case, use 12 point
 * midpoint rule across bandpass.
 */
double get_MTF(double u, double v, PSF_DATA *psf) {

  int j;
  double mtftot = 0;
  double loffset;
  PSF_DATA psf2;

  /* Monochromatic PSFs can simply be returned using get_MTF_mono. */
  if (psf->is_broadband==0) return(get_MTF_mono(u,v,psf));

#define N_LAMBDA_INTEG 12
  psf2.pixscale = psf->pixscale;
  psf2.sigma = psf->sigma;
  psf2.centobs = psf->centobs;
  psf2.is_broadband = 0;
  for(j=0;j<N_LAMBDA_INTEG;j++) {
    loffset = exp(psf->dll*((j+0.5)/(double)N_LAMBDA_INTEG-0.5));
    psf2.lD = psf->lD * loffset;
    psf2.rmswfe = psf->rmswfe / loffset;
    mtftot += get_MTF_mono(u,v,&psf2);
  }
  return(mtftot/N_LAMBDA_INTEG);
#undef N_LAMBDA_INTEG
}

/* --- ROUTINES TO COMPUTE GALAXY IMAGE PROPERTIES --- */

/* Fraction of encircled energy in radius r (arcsec) for galaxy with stated
 * effective radius (assumed exponential) and given PSF.  The result is given by:
 *
 * = r * int dk <G~(k) f~(k)> J_1(kr).
 *
 * where f~(k) is the FT of the galaxy profile, and <> is an angular average.
 * In its present version, the galaxy is assumed circular and the angular average
 * is done with angular steps of 10 degrees.
 *
 * The function is well behaved if reff=0.
 */
double get_gal_fracE(double r, double reff, PSF_DATA *psf) {
  double k, k2pi, kmax, dk, I, rs, dphi;
  long j,ndk,jphi,nphi;
  double u,v,tildef;

  /* Get exponential scale length */
  rs = reff/RAT_HL_SL_EXP;

  /* Get maximum wavevector and integration step */
  kmax = 1./psf->lD;
  if (psf->sigma*kmax>6) kmax=6./psf->sigma;
  if (rs*kmax>25) kmax=25/rs;
  kmax *= 2.*M_PI;
  if (kmax*r<1) {
    dk = 0.005*kmax;
  } else {
    dk = 0.005/r;  
  }
  ndk = (long)ceil(kmax/dk);
  dk = kmax/ndk;
  nphi = 9;
  dphi = M_PI/(double)nphi;

  /* The actual integration: outer loop over k, inner loop over angle */
  I=0;
  for(j=0;j<ndk;j++) {
    k = (j+0.5)*dk;
    tildef = pow(1.+k*k*rs*rs, -1.5);
    k2pi = k/(2.*M_PI);
    for(jphi=0;jphi<nphi;jphi++) {
      u = k2pi*cos(dphi*jphi);
      v = k2pi*sin(dphi*jphi);
      I += get_MTF(u,v,psf) * getJ1(k*r) * tildef;
    }
  }
  I *= dk*r/(double)nphi;
      
  return(I);
}

/* Gets the apparent galaxy size (measured by encircled energy radius
 * with the specified energy fraction) for given reff (circular exponential
 * profile) and PSF. Both reff and output are in arcsec.
 *
 * This function essentially inverts get_gal_fracE.
 */
double get_gal_size(double eFrac, double reff, PSF_DATA *psf) {
  double r, ratio;
     
  r=0.2; ratio=0.5;
  while(get_gal_fracE(r,reff,psf)>eFrac && r<206265) r/=exp(ratio);
  while(get_gal_fracE(r,reff,psf)<eFrac && r>1e-6) r*=exp(ratio);
  
  for(ratio=0.5; ratio>1e-7; ratio/=2.)
    r *= get_gal_fracE(r,reff,psf)>eFrac? exp(-ratio): exp(ratio);
 
  return(r);
}

/* Gets the effective solid angle of the galaxy in arcsec^2 for noise
 * purposes, i.e.
 *
 * Om_eff^{-1} = int [G~(u) f~(u)]^2 d^2u = (2pi)^-2 int [G~(k) f~(k)]^2 d^2k,
 *
 * where G~ = MTF and f~ = galaxy Fourier transform.
 * This is defined so that if one had a perfect template for the location and
 * shape of the galaxy, and if one were background limited, the variance of the
 * galaxy counts would be the number of background photons in an area Om_eff.
 */
double get_gal_Omeff(double reff, PSF_DATA *psf) {
  double k, k2pi, kmax, dk, I, rs, dphi;
  long j,ndk,jphi,nphi;
  double u,v,tildef,prod;

  /* Get exponential scale length */
  rs = reff/RAT_HL_SL_EXP;

  /* Get maximum wavevector and integration step */
  kmax = 1./psf->lD;
  if (psf->sigma*kmax>6) kmax=6./psf->sigma;
  if (rs*kmax>25) kmax=25/rs;
  kmax *= 2.*M_PI;
  dk = DELTA_STEP*kmax;
  ndk = (long)ceil(kmax/dk);
  dk = kmax/ndk;
  nphi = 9;
  dphi = M_PI/(double)nphi;

  /* The actual integration: outer loop over k, inner loop over angle */
  I=0;
  for(j=0;j<ndk;j++) {
    k = (j+0.5)*dk;
    k2pi = k/(2.*M_PI);
    tildef = pow(1.+k*k*rs*rs, -1.5);
    for(jphi=0;jphi<nphi;jphi++) {
      u = k2pi*cos(jphi*dphi);
      v = k2pi*sin(jphi*dphi);
#ifdef BAO_OR_SPCONT_MODE
      if (fabs(u)<0.5/psf->pixscale && fabs(v)<0.5/psf->pixscale) {
#endif
      prod = get_MTF(u,v,psf) * tildef;
      I += k*prod*prod;
#ifdef BAO_OR_SPCONT_MODE
      }
#endif
    }
  }
  I *= dk/(double)nphi/2./M_PI;

  return(1./I);
}

/* Gets the peak surface intensity of the image in arcsec^-2, i.e.
 *
 * PeakInt = int [G~(u) f~(u)] d^2u = (2pi)^-2 int [G~(k) f~(k)] d^2k,
 *
 * where G~ = MTF and f~ = galaxy Fourier transform.
 * The main intent is to use this to study saturation of stars.
 */
double get_gal_PeakInt(double reff, PSF_DATA *psf) {
  double k, k2pi, kmax, dk, I, rs, dphi;
  long j,ndk,jphi,nphi;
  double u,v,tildef,prod;

  /* Get exponential scale length */
  rs = reff/RAT_HL_SL_EXP;

  /* Get maximum wavevector and integration step */
  kmax = 1./psf->lD;
  if (psf->sigma*kmax>6) kmax=6./psf->sigma;
  if (rs*kmax>25) kmax=25/rs;
  kmax *= 2.*M_PI;
  dk = DELTA_STEP*kmax;
  ndk = (long)ceil(kmax/dk);
  dk = kmax/ndk;
  nphi = 9;
  dphi = M_PI/(double)nphi;

  /* The actual integration: outer loop over k, inner loop over angle */
  I=0;
  for(j=0;j<ndk;j++) {
    k = (j+0.5)*dk;
    k2pi = k/(2.*M_PI);
    tildef = pow(1.+k*k*rs*rs, -1.5);
    for(jphi=0;jphi<nphi;jphi++) {
      u = k2pi*cos(jphi*dphi);
      v = k2pi*sin(jphi*dphi);
      prod = get_MTF(u,v,psf) * tildef;
      I += k*prod;
    }
  }
  I *= dk/(double)nphi/2./M_PI;

  return(I);
}

/* Gets the effective trace width of the galaxy (units: arcsec), in terms
 * of the PSF and the effective radius (units: arcsec).
 *
 * theta_eff^{-1} = int [G~(u) f~(u)]^2 du_2 (eval. at u_1=0).
 */
double get_gal_thetaeff(double reff, PSF_DATA *psf) {

  double rs, kmax, dk, k, u, I, prod, tildef;
  long ndk, j;

  /* Get exponential scale length */
  rs = reff/RAT_HL_SL_EXP;

  /* Get maximum wavevector and integration step */  
  kmax = 1./psf->lD;
  if (psf->sigma*kmax>6) kmax=6./psf->sigma;
  if (rs*kmax>25) kmax=25/rs;
  kmax *= 2.*M_PI;
  dk = DELTA_STEP*kmax;
  ndk = (long)ceil(kmax/dk);
  dk = kmax/ndk;

  /* The actual integration: loop over k only, no angle */
  I=0;
  for(j=0;j<ndk;j++) {
    k = (j+0.5)*dk;
    u = k/(2.*M_PI);
    tildef = pow(1.+k*k*rs*rs, -1.5);
#ifdef BAO_OR_SPCONT_MODE
    if (fabs(u)<0.5/psf->pixscale) {
#endif
    prod = get_MTF(0,u,psf) * tildef;
    I += prod*prod;
#ifdef BAO_OR_SPCONT_MODE
    }
#endif
  }
  I *= dk/M_PI; /* div by 2pi to get k->u, multiply by 2 for negative wavenumbers */
  return(1./I);
}

/* Gets the coupling integral for the galaxy (units: arcsec^-1),
 *
 * J = int [G~(u) f~(u)]^2 du_2 = (2pi)^-1 int [G~(k) f~(k)]^2 dk_2,
 *
 * where G~ = MTF and f~ = galaxy Fourier transform.
 * This is defined so that when extracting a line, the effective contribution to
 * the background counts from the source continuum is J times the source counts
 * per arcsec.
 */
double get_J_couplingIntegral(double reff, PSF_DATA *psf) {
  double k, k2pi, kmax, dk, I, rs;
  long j,ndk;
  double tildef,prod;

  /* Get exponential scale length */
  rs = reff/RAT_HL_SL_EXP;

  /* Get maximum wavevector and integration step */
  kmax = 1./psf->lD;
  if (psf->sigma*kmax>6) kmax=6./psf->sigma;
  if (rs*kmax>25) kmax=25/rs;
  kmax *= 2.*M_PI;
  dk = DELTA_STEP*kmax;
  ndk = (long)ceil(kmax/dk);
  dk = kmax/ndk;

  /* The actual integration: loop over k; only do k>0, double for k<0 */
  I=0;
  for(j=0;j<ndk;j++) {
    k = (j+0.5)*dk;
    k2pi = k/(2.*M_PI);
    if (fabs(k2pi)<0.5/psf->pixscale) {
      tildef = pow(1.+k*k*rs*rs, -1.5);
      prod = get_MTF(0,k2pi,psf) * tildef;
      I += prod*prod;
    }
  }
  I *= dk/M_PI;

  return(I);
}

/* Gets the centroiding integral of the galaxy in arcsec^-4 for centroiding
 * purposes, i.e.
 *
 * c_int = int [G~(u) f~(u) 2 pi u_x]^2 d^2u = (2pi)^-2 int [G~(k) f~(k) k_x]^2 d^2k,
 *
 * where G~ = MTF and f~ = galaxy Fourier transform.
 * This is defined so that if one had a perfect template for the location and
 * shape of the galaxy, and if one were background limited, the inverse-variance
 * of the galaxy centroid (in arcsec^-2) would be
 *
 * ivar(centroid) = c_int * [galaxy counts]^2/[background counts per solid angle]
 *
 * This is just the statistical centroid error and does not account for aliasing --
 * so should be used with caution for very small sources or if it gives errors of a
 * tiny fraction of a pixel.
 */
double get_gal_centroid_integral(double reff, PSF_DATA *psf) {
  double k, k2pi, kmax, dk, I, rs, dphi;
  long j,ndk,jphi,nphi;
  double u,v,tildef,prod;

  /* Get exponential scale length */
  rs = reff/RAT_HL_SL_EXP;

  /* Get maximum wavevector and integration step */
  kmax = 1./psf->lD;
  if (psf->sigma*kmax>6) kmax=6./psf->sigma;
  if (rs*kmax>25) kmax=25/rs;
  kmax *= 2.*M_PI;
  dk = DELTA_STEP*kmax;
  ndk = (long)ceil(kmax/dk);
  dk = kmax/ndk;
  nphi = 9;
  dphi = M_PI/(double)nphi;

  /* The actual integration: outer loop over k, inner loop over angle */
  I=0;
  for(j=0;j<ndk;j++) {
    k = (j+0.5)*dk;
    k2pi = k/(2.*M_PI);
    tildef = pow(1.+k*k*rs*rs, -1.5);
    for(jphi=0;jphi<nphi;jphi++) {
      u = k2pi*cos(jphi*dphi);
      v = k2pi*sin(jphi*dphi);
#ifdef BAO_OR_SPCONT_MODE
      if (fabs(u)<0.5/psf->pixscale && fabs(v)<0.5/psf->pixscale) {
#endif
      prod = get_MTF(u,v,psf) * tildef;
      I += k*prod*prod*u*u;
#ifdef BAO_OR_SPCONT_MODE
      }
#endif
    }
  }
  I *= dk/(double)nphi*2.*M_PI;

  return(I);
}

/* Gets the shape penalty factor for a PSF, i.e. the uncertainty per ellipticity
 * component is
 *
 * sigma(e) = 2/nu * sqrt(penalty_factor)
 *
 * where nu is the detection S/N. This is 1 for a well-resolved Gaussian,
 * and equal to the (resolution factor)^-2 for general Gaussian PSF+galaxy. This
 * routine does the full integral for an exponential profile galaxy and the true PSF.
 *
 * This can be computed from the effective area and the shape effective area,
 *
 * A_eff^{-1} = int [G~(u) f~(u)]^2 d^2u = (2pi)^-2 int [G~(k) f~(k)]^2 d^2k,
 * A_{eff,s}^{-1} = int [G~(u) d/dlnu f~(u)]^2 d^2u = (2pi)^-2 int [G~(k) d/dlnk f~(k)]^2 d^2k,
 *
 * where G~ = MTF and f~ = galaxy Fourier transform.
 * The penalty factor is A_{eff,s}/A_eff.
 * See Bernstein & Jarvis (2002) for a derivation of the 2nd equation.
 */
double get_shape_penalty_factor(double reff, PSF_DATA *psf) {
  double k, k2pi, kmax, dk, I, rs, dphi;
  long j,ndk,jphi,nphi;
  double u,v,dtildefdlnk,prod;

  /* Get exponential scale length */
  rs = reff/RAT_HL_SL_EXP;

  /* Get maximum wavevector and integration step */
  kmax = 1./psf->lD;
  if (psf->sigma*kmax>6) kmax=6./psf->sigma;
  if (rs*kmax>30) kmax=30/rs;
  kmax *= 2.*M_PI;
  dk = DELTA_STEP*kmax;
  ndk = (long)ceil(kmax/dk);
  dk = kmax/ndk;
  nphi = 9;
  dphi = M_PI/(double)nphi;

  /* The actual integration: outer loop over k, inner loop over angle */
  I=0;
  for(j=0;j<ndk;j++) {
    k = (j+0.5)*dk;
    k2pi = k/(2.*M_PI);
    dtildefdlnk = -pow(1.+k*k*rs*rs, -2.5) * 2*k*k*rs*rs;
    for(jphi=0;jphi<nphi;jphi++) {
      u = k2pi*cos(jphi*dphi);
      v = k2pi*sin(jphi*dphi);
      prod = get_MTF(u,v,psf) * dtildefdlnk;
      I += k*prod*prod;
    }
  }
  I *= dk/(double)nphi/2./M_PI;

  return(1./I/get_gal_Omeff(reff,psf));
}

/* Gets the 1 sigma flux uncertainty for a line in a single exposure for a galaxy with:
 * effective radius reff (arcsec)
 * specified psf
 * variance per unit area (e-^2/arcsec^2)
 * calib_1exp (e- per (W/m2))
 *
 * Output in W/m2.
 */
double get_1sigma_flux_1exp(double reff, PSF_DATA *psf, double var_1exp, double calib_1exp) {
  double Omeff;

  /* The effective area of the galaxy in arcsec^2 */
  Omeff = get_gal_Omeff(reff,psf);

  return(sqrt(var_1exp*Omeff)/calib_1exp);
}

/* Gets the Z sigma flux uncertainty for a line in a single exposure for a galaxy with:
 * effective radius reff (arcsec)
 * specified psf
 * variance per unit area (e-^2/arcsec^2)
 * calib_1exp (e- per (W/m2))
 * equivalent width*dispersion EWD (arcsec)
 *
 * Output in W/m2.
 */
double get_limflux(double reff, PSF_DATA *psf, double var_1exp, double calib_1exp, double EWD,
  double Nexp, double Z) {

  double Omeff, J, fd;
  double alpha, beta;

  /* Galaxy effective solid angle and coupling integral J */
  Omeff = get_gal_Omeff(reff,psf);
  J = get_J_couplingIntegral(reff,psf);
  fd = DET_IS_CCD? 1.: 1.2;

  /* Set up equation for limiting flux:
   * F = Z * sqrt((var_1exp + fd*J/EWD*F*calib_1exp)*Omeff)/calib_1exp / sqrt(Nexp)
   *
   * The solution to this equation is to first find the coefficients by squaring:
   * F^2 = alpha + 2*beta*F
   * and then
   * F = beta + sqrt(beta^2 + alpha)
   */
  alpha = Z*Z/calib_1exp/calib_1exp*var_1exp*Omeff/Nexp;
  beta = Z*Z*fd*J/EWD*Omeff/calib_1exp/2./Nexp;
  return(beta+sqrt(alpha+beta*beta));
}

/* Gets limiting AB magnitude for successful shape measurement on a source of
 * given r_eff, number of exposures, noise (var_1exp), and calibration (calib_1exp),
 * and the specified maximum ellipticity error.
 *
 * If WL_SNR_MIN is defined, we also require a detection at this number of sigmas.
 */
double get_maglim_shape(double r_eff, PSF_DATA *psf, int N_exp, double var_1exp,
  double calib_1exp, double max_ellip_err) {

  double SNR0, sigma_e, Fmin, Fmin_shape, Fmin_SNR;

  /* First get shape measurement error for AB mag = 0 galaxy */
  SNR0 = calib_1exp * sqrt(N_exp/var_1exp/get_gal_Omeff(r_eff,psf));
  sigma_e = 2./SNR0 * sqrt(get_shape_penalty_factor(r_eff,psf));
  Fmin_shape = sigma_e/max_ellip_err;

  Fmin_SNR = 0;
#ifdef WL_SNR_MIN
  Fmin_SNR = WL_SNR_MIN / SNR0;
#endif
  Fmin = Fmin_shape>Fmin_SNR? Fmin_shape: Fmin_SNR;

  /* Re-scale to specified error. */
  return(-2.5*log(Fmin)/log(10.));
}

/* Gets limiting AB magnitude for successful shape measurement on a source of
 * given r_eff, number of exposures, noise (var_1exp), and calibration (calib_1exp),
 * and the specified maximum ellipticity error. NO S/N cut is imposed in this case;
 * it is useful for estimating the ellipticity error by comparing the magnitude of a
 * particular galaxy to the "limiting" magnitude (with known sigma_e).
 */
double get_maglim_shape_noSNRcut(double r_eff, PSF_DATA *psf, int N_exp, double var_1exp,
  double calib_1exp, double max_ellip_err) {

  double SNR0, sigma_e, Fmin, Fmin_shape, Fmin_SNR;

  /* First get shape measurement error for AB mag = 0 galaxy */
  SNR0 = calib_1exp * sqrt(N_exp/var_1exp/get_gal_Omeff(r_eff,psf));
  sigma_e = 2./SNR0 * sqrt(get_shape_penalty_factor(r_eff,psf));
  Fmin_shape = sigma_e/max_ellip_err;

  Fmin_SNR = 0;
  Fmin = Fmin_shape>Fmin_SNR? Fmin_shape: Fmin_SNR;

  /* Re-scale to specified error. */
  return(-2.5*log(Fmin)/log(10.));
}

/* --- ROUTINES TO COMPUTE COSMOLOGICAL INFORMATION --- */

/* Comoving distance in Mpc */
double computeDistance(double z) {
  double r, ztemp, y;
  
  r=0;
  for(ztemp=0.0005*z; ztemp<z; ztemp+=0.001*z) {
    y = 1.+ztemp;
    r += 1./sqrt(1.-COSMO_OMEGA_M+COSMO_OMEGA_M*y*y*y);
  }
  r *= 0.001*z;
  return(299792.458*r/COSMO_H0);
}

/* Hubble rate in c/Mpc */
double computeHubble(double z) {
  double y;
  y = 1.+z;
  return(COSMO_H0*sqrt(1.-COSMO_OMEGA_M+COSMO_OMEGA_M*y*y*y)/299792.458);
}

/* --- ROUTINES TO COMPUTE PROPERTIES OF THE BAO GALAXY DISTRIBUTION --- */

/* Obtains the lognormal distribution of physical galaxy sizes for the specified
 * population model.
 *
 * Input is redshift and Ha flux (in W/m^2), output is the median reff (arcsec)
 * and dispersion (in e-folds). There is plumbing for several models.
 */
void get_galsizes(double z, double FHa, double *med_reff, double *dispersion, int model) {

  double xi;

#ifdef OIII_GAL
  /* == [O III] models == */
  /* Model from Sangeeta Malhotra, 08/15/16 */
  if (model%10==2) {
    *med_reff = 0.263;
    *dispersion = 0.502;
    return;
  }

  xi=0.;

#elif OII_GAL
  /* == [O II] models == */

  /* H alpha model 2, based on [O II] -> H alpha conversion
   * log10 reff = -0.66 + 0.27*log10 (FHa/1e-16cgs) +/- 0.20
   * uses [O II] -> H alpha conversion of 2.27
   * intrinsic ratio of 1.77 from Kennicutt (1998)
   * dust reddening of 0.27 mag between the two lines following Khostovan et al. (2015)
   */
  if (model%10==2) {
    xi = (1.75-LAMBDA0_HA*(1+z))/0.45;
    *med_reff = 0.219*pow(10.,0.06*xi)*pow(2.27*FHa/1e-19,0.27+0.01*xi);
    *dispersion = 0.46;
    return;
  }

  xi=0.;

#else
  /* == H alpha models == */

  /* Default model:
   * Fit to COSMOS Mock catalog in range of 1.5<lambda<2.0 um observed H alpha
   * Lognormal reff distributions fit in 5 bins centered at log10 FHa(cgs) = -16.0(0.2)-15.2
   * Can be approximated to ~0.01dex (except for last bin) by:
   * log10 reff = -0.62 + 0.25*log10 (FHa/1e-16cgs) +/- 0.20
   */
  if (model%10==0) {
    *med_reff = 0.240*pow(FHa/1e-19,0.25);
    *dispersion = 0.46;
    return;
  }
  /* Updated CMC model: (August 15, 2011 version)
   * log10 reff = -0.66 + 0.27*log10 (FHa/1e-16cgs) +/- 0.20
   */
  if (model%10==2) {
    xi = (1.75-LAMBDA0_HA*(1+z))/0.45;
    *med_reff = 0.219*pow(10.,0.06*xi)*pow(FHa/1e-19,0.27+0.01*xi);
    *dispersion = 0.46;
    return;
  }
#endif

  /* All galaxies at reff=0.3" (for test purposes only) */
  if (model%10==1) {
    *med_reff = 0.3;
    *dispersion = 0;
    return;
  }

  fprintf(stderr, "Error: get_galsizes: illegal model number %d\n", model);
  wait_for_user();
  exit(1);
}

/* Returns the H alpha luminosity function in observer frame flux.
 * Input fluxes in W/m2.
 *
 * Output is d[Number objects]/d[Comoving volume]/d[ln FHa]
 * Units are Mpc^-3.
 *
 * The associated function, print_HaLF_model, prints the model name
 * to the specified file. If additional models are incorporated, it
 * should be updated at the same time.
 *
 * The [NII]/Ha ratio assumed by each model is given by get_HaLF2 to the pointer NII_Ha.
 *
 * Models that are 1000-1999 switch to [O III].
 * Models that are 2000-2999 switch to [O II].
 */
double get_HaLF(double z, double FHa, int model) {
  double get_HaLF2(double,double,int,double*);
  return(get_HaLF2(z,FHa,model,NULL));
}
double get_HaLF2(double z, double FHa, int model, double *NII_Ha) {

  double L, D, Lstar, phistar, alpha;
#ifdef OIII_GAL
  double xi;
#endif
#ifdef OII_GAL
  double xi;
#endif

  /* Average of H alpha */
  if (model/10==99)
    return((get_HaLF2(z,FHa,60+model%10,NII_Ha) + get_HaLF2(z,FHa,70+model%10,NII_Ha) + get_HaLF2(z,FHa,80+model%10,NII_Ha))/3.);
  /* Average of [OIII] */
  if (model/10==199)
    return((get_HaLF2(z,FHa,1000+model%10,NII_Ha) + get_HaLF2(z,FHa,1010+model%10,NII_Ha) + get_HaLF2(z,FHa,1020+model%10,NII_Ha))/3.);

  /* Get distance (in meters) and L (in Watts) */
  D = computeDistance(z) * METERS_PER_MPC;
  L = 4.*M_PI*(1+z)*(1+z)*FHa*D*D;

#ifdef OIII_GAL
  /* == [O III] models here == */

  if (NII_Ha!=NULL)
    *NII_Ha=0.; /* turn off the "[N II]" feature since [O III] is not blended with another line */

  /* model used for FSWG Phase A updates
   * Mehta et al. (2015) -- not cosmology corrected
   */
  if (model/10==100) {
    xi = z>2.2? 1.2: z<0.8? -0.2: z-1;
    alpha = -1.42-0.15*xi;
    Lstar = pow(10., 35.21+0.35*xi);
    phistar = pow(10., -3.18+0.46*xi);
    return(phistar*pow(L/Lstar,alpha+1)*exp(-L/Lstar));
  }

  /* Colbert et al. (2013) model */
  if (model/10==101) {
    xi = z>2.3? 1.5: z<0.7? -0.5: (z-1.1)/.8;
    alpha = -1.40-0.27*xi;
    Lstar = pow(10., 35.34+0.58*xi);
    phistar = pow(10., -3.20-0.57*xi);
    return(phistar*pow(L/Lstar,alpha+1)*exp(-L/Lstar));
  }

  /* Khostovan et al. (2015) model */
  if (model/10==102) {
    alpha = -1.60;
    Lstar = pow(10., 35.85);
    phistar = pow(10., -3.35);
    if (z<3.24) {
      xi = (z-2.23)/1.01;
      Lstar = pow(10., 35.67+0.18*xi);
      phistar = pow(10., -3.06-0.29*xi);
    }
    if (z<2.23) {
      xi = (z-1.42)/0.81;
      Lstar = pow(10., 35.07+0.60*xi);
      phistar = pow(10., -2.63-0.43*xi);
    }
    if (z<1.42) {
      xi = (z-0.84)/0.58;
      Lstar = pow(10., 34.79+0.28*xi);
      phistar = pow(10., -2.56-0.07*xi);
    }
    Lstar *= 1.+1./OIII_RATIO;
    return(phistar*pow(L/Lstar,alpha+1)*exp(-L/Lstar));
  }

  fprintf(stderr, "Error: get_HaLF: illegal model number %d\n", model);
  wait_for_user();
  exit(1);
#endif

#ifdef OII_GAL
  /* == [O II] models here == */

  if (NII_Ha!=NULL)
    *NII_Ha=0.; /* turn off the "[N II]" feature since [O III] is not blended with another line */

  if (model/10==200) {
    alpha=-1.30;
    Lstar = pow(10., 35.95);
    phistar = pow(10., -3.73);
    if (z<4.69) {
      xi = (z-3.34)/1.35;
      Lstar = pow(10., 35.71+0.24*xi);
      phistar = pow(10., -3.11-0.62*xi);
    }
    if (z<3.34) {
      xi = (z-2.25)/1.09;
      Lstar = pow(10., 35.35+0.36*xi);
      phistar = pow(10., -2.51-0.60*xi);
    }
    if (z<2.25) {
      xi = (z-1.47)/0.78;
      Lstar = pow(10., 34.87+0.48*xi);
      phistar = pow(10., -2.27-0.24*xi);
    }
    return(phistar*pow(L/Lstar,alpha+1)*exp(-L/Lstar));
  }

  fprintf(stderr, "Error: get_HaLF: illegal model number %d\n", model);
  wait_for_user();
  exit(1);
#endif

  /* == H alpha models below here == */

  /* Default model: Geach et al (2009) LF
   * Monthly Notices of the Royal Astronomical Society, Volume 402, Issue 2, pp. 1330-1338.
   * Schechter function with phi* = 1.35e-3 Mpc^-3, alpha=-1.35, L*=5.1e34(1+z)^3.1 W
   * (except above z=1.3 fix at L*=6.8e35 W)
   */
  if (model/10==0) {
    Lstar = z>1.3? 6.8e35: 5.1e34*pow(1+z,3.1);
    if (NII_Ha!=NULL) *NII_Ha = 0.43;
    return(1.35e-3*pow(L/Lstar,-0.35)*exp(-L/Lstar));
  }

  /* 90% confidence lower limit: Geach et al LF
   * (-1.28155 sigma)
   * Modified Geach et al but with L*=5.1e34(z>1.3?2.3:1+z)^2.59 W
   */
  if (model/10==1) {
    Lstar = 5.1e34*pow(z>1.3?2.3:1+z,2.59);
    if (NII_Ha!=NULL) *NII_Ha = 0.43;
    return(1.35e-3*pow(L/Lstar,-0.35)*exp(-L/Lstar));
  }

  /* Geach et al / 1.257 -- as used for Euclid.
   * See Yun Wang email.
   */
  if (model/10==2) {
    Lstar = z>1.3? 6.8e35: 5.1e34*pow(1+z,3.1);
    if (NII_Ha!=NULL) *NII_Ha = 0.43;
    return(1.35e-3*pow(L/Lstar,-0.35)*exp(-L/Lstar) / 1.257);
  }

  /* Sobral et al 2012 with cosmology corrections (to WMAP5) and aperture corrections
   * assuming exponential profile with median r_eff at 1e-19 W/m^2.
   */
  if (model/10==3) {
    Lstar = z<0.84? 34.70176145+0.500472054*(z-0.40): z<1.47? 34.92196916+0.495788251*(z-0.84): z<2.23? 35.23431576+0.394312857*(z-1.47): 35.53399353;
    Lstar = pow(10,Lstar);
    phistar = z<0.84? -3.107120757+1.434831821*(z-0.40): z<1.47? -2.475794755-0.247531478*(z-0.84): z<2.23? -2.631739587-0.170782235*(z-1.47): -2.761534085;
    phistar = pow(10,phistar);
    if (NII_Ha!=NULL) *NII_Ha = 0.33;
    return(phistar*pow(L/Lstar,-0.6)*exp(-L/Lstar));
  }

  /* Colbert et al 2013 fit -- private communication, 2013-02-13 */
  if (model/10==4) {
    Lstar = 35.17 + 0.32*(z<1.5? (z-1.2)/0.6: 0.5);
    Lstar = pow(10,Lstar);
    phistar = -2.33 -0.13*(z<1.5? (z-1.2)/0.6: 0.5);
    phistar = pow(10,phistar);
    alpha = -1.27 -0.07*(z<1.5? (z-1.2)/0.6: 0.5);
    if (NII_Ha!=NULL) *NII_Ha = 0.41;
    return(phistar*pow(L/Lstar,alpha+1.)*exp(-L/Lstar));
  }

  /* Colbert et al 2013 fit -- as presented in the paper, and corrected to WMAP5 cosmology.
   * Parameters are:
   * z=0.6: logPhistar = -2.51+0.00, logLstar = 34.72-0.01, alpha = -1.27
   * z=1.2: logPhistar = -2.70-0.01, logLstar = 35.18+0.00, alpha = -1.43
   */
  if (model/10==5) {
    Lstar = 35.18 + 0.47*(z<1.5? (z-1.2)/0.6: 0.5);
    Lstar = pow(10,Lstar);
    phistar = -2.71 -0.20*(z<1.5? (z-1.2)/0.6: 0.5);
    phistar = pow(10,phistar);
    alpha = -1.43 -0.16*(z<1.5? (z-1.2)/0.6: 0.5);
    if (NII_Ha!=NULL) *NII_Ha = 0.41;
    return(phistar*pow(L/Lstar,alpha+1.)*exp(-L/Lstar));
  }

  /* Global fit to HIZELS + WISP + NICMOS (by Hirata: Euclid Model 3) */
  if (model/10==6) {
    Lstar = pow(10., 35.956 -1.223*pow(1.5/(1.+z), 1.615) );
    if (NII_Ha!=NULL) *NII_Ha = 0.37;
    return(1.202e-3*pow(L/Lstar,-0.587)/(1.+1.718281828*pow(L/Lstar,2.288)));
  }

  /* Euclid Model 1 */
  if (model/10==7) {
    Lstar = 3.1622776e34 * (1.+z) * (1.+z);
    phistar = 1.585e-3 * (1.+z) * (z>1.3? 2.3*2.3/(1+z)/(1+z): 1);
    alpha = -1.35;
    if (NII_Ha!=NULL) *NII_Ha = 0.37;
    return(phistar*pow(L/Lstar,alpha+1.)*exp(-L/Lstar));
  }

  /* Euclid Model 2 */
  if (model/10==8) {
    Lstar = pow(10., 35.59 - 0.22*(z-2.23)*(z-2.23));
    phistar = 1.995e-3;
    alpha = -1.4;
    if (NII_Ha!=NULL) *NII_Ha = 0.37;
    return(phistar*pow(L/Lstar,alpha+1.)*exp(-L/Lstar));
  }

  fprintf(stderr, "Error: get_HaLF: illegal model number %d\n", model);
  wait_for_user();
  exit(1);
}

void print_HaLF_model(FILE *fp, int model) {

#ifdef OIII_GAL
  if (model/10==100) { fprintf(fp, "Mehta et al. (2015) [O III]"); return; }
  if (model/10==101) { fprintf(fp, "Colbert et al. (2013) [O III]"); return; }
  if (model/10==102) { fprintf(fp, "Khostovan et al. (2015) H beta + [O III]"); return; }
  if (model/10==199) { fprintf(fp, "[O III] 3-model ave."); return; }
#elif OII_GAL
  if (model/10==200) { fprintf(fp, "Khostovan et al. (2015) [O II]"); return; }
#else
  if (model/10==0) { fprintf(fp, "Geach et al (2009), 50%c confidence", '%'); return; }
  if (model/10==1) { fprintf(fp, "Geach et al (2009), 90%c confidence", '%'); return; }
  if (model/10==2) { fprintf(fp, "Geach et al (2009), reduced /1.257"); return; }
  if (model/10==3) { fprintf(fp, "Sobral et al (2012), corrected, best fit"); return; }
  if (model/10==4) { fprintf(fp, "Colbert et al (2013), prepublication, best fit"); return; }
  if (model/10==5) { fprintf(fp, "Colbert et al (2013), best fit"); return; }
  if (model/10==6) { fprintf(fp, "Global fit, HIZELS+WISP+NICMOS - Model 3"); return; }
  if (model/10==7) { fprintf(fp, "Global fit, Euclid - Model 1 (LDE)"); return; }
  if (model/10==8) { fprintf(fp, "Global fit, Euclid - Model 2 (LE)"); return; }
  if (model/10==99) { fprintf(fp, "H alpha 3-model ave."); return; }
#endif

  fprintf(stderr, "Error: print_HaLF_model: illegal model number %d\n", model);
  wait_for_user();
  exit(1);
}

/* --- ROUTINES TO COMPUTE DENSITIES OF OBSERVABLE GALAXIES: BAO --- */

/* S/N boost factor from having [NII] as well as H alpha.
 * Takes as input the ratio [NII]/Ha (where [NII] is the *total* of the [NII] lines).
 * Based on methodology from the effective area determination.
 */
double get_SNR_boost_factor(double z, double reff, double NII_Ha, PSF_DATA *psf) {

  double A[3], x[3];
  double k, k2pi, kmax, dk, rs, dphi;
  long j,ndk,jphi,nphi;
  double u,v,tildef,prod,triplet_trans_sq;
  double I1, I2;  /* The integral I1 will represent the total SNR^2 integral, whereas
                   * I2 will include only the H alpha. The SNR boost is the sqrt of the ratio.
                   */

  /* Strengths of the two [NII] lines, and their offsets in arcsec. */
  A[0] = 1.;
  A[1] = NII_Ha/4.05;
  A[2] = NII_Ha-A[1];
  x[0] = 0.;
  x[1] = -0.00147*(1+z)*psf->linDisp;
  x[2] =  0.00206*(1+z)*psf->linDisp;

  /* Get exponential scale length */
  rs = reff/RAT_HL_SL_EXP;

  /* Get maximum wavevector and integration step */
  kmax = 1./psf->lD;
  if (psf->sigma*kmax>6) kmax=6./psf->sigma;
  if (rs*kmax>25) kmax=25/rs;
  kmax *= 2.*M_PI;
  dk = DELTA_STEP*kmax;
  ndk = (long)ceil(kmax/dk);
  dk = kmax/ndk;
  nphi = 9;
  dphi = M_PI/(double)nphi;
   
  /* The actual integration: outer loop over k, inner loop over angle */
  I1 = I2 = 0;
  for(j=0;j<ndk;j++) {
    k = (j+0.5)*dk;
    k2pi = k/(2.*M_PI);
    tildef = pow(1.+k*k*rs*rs, -1.5);
    for(jphi=0;jphi<nphi;jphi++) {
      u = k2pi*cos(jphi*dphi);
      v = k2pi*sin(jphi*dphi);
#ifdef BAO_OR_SPCONT_MODE
      if (fabs(u)<0.5/psf->pixscale && fabs(v)<0.5/psf->pixscale) {
#endif
      prod = get_MTF(u,v,psf) * tildef;
      triplet_trans_sq = A[0]*A[0] + A[1]*A[1] + A[2]*A[2]
                         + 2.*A[0]*A[1]*cos(2.*M_PI*u*(x[0]-x[1]))
                         + 2.*A[0]*A[2]*cos(2.*M_PI*u*(x[0]-x[2]))
                         + 2.*A[1]*A[2]*cos(2.*M_PI*u*(x[1]-x[2]));
        /* Note only u appeared in the squared transfer function for the 'triplet' since
         * the lines are separated only in the x-direction.
         */
      I1 += k*prod*prod*triplet_trans_sq;
      I2 += k*prod*prod;
#ifdef BAO_OR_SPCONT_MODE
      }
#endif
    }
  }
  I1 *= dk/(double)nphi/2./M_PI;
  I2 *= dk/(double)nphi/2./M_PI;

  return(sqrt(I1/I2));
}

/* Probability of a given galaxy being above the detection threshold at specified
 * redshift and H alpha flux. This does not include photometric scatter.
 *
 * Input units: FHa (W/m2), var_1exp(e-^2/pix), calib_1exp(e- per W/m^2)
 */
double get_P_Detectable(double z, double FHa, PSF_DATA *psf, double var_1exp, double calib_1exp,
  int N_exp, double significance_cut, int model) {

  double P, x, dx, med_reff, dispersion, reff, sf, NII_Ha, boost;

  /* Get the median effective radius and its log-dispersion */
  get_galsizes(z, FHa, &med_reff, &dispersion, model);

  /* Need to determine the critical number of sigmas of the lognormal distribution in size at
   * which the galaxies are detected. Use x = ln(reff/med_reff)/dispersion.
   * If the critical threshold is outside 8 sigma, this bisection routine will fail. However if
   * this happens then the galaxies are essentially all detected or not anyway.
   */
  x=0;
  for(dx=4;dx>1e-6;dx/=2.) {
    reff = med_reff * exp(x*dispersion);
    boost = 1.; NII_Ha = 0.;
#ifdef OIII_GAL
    boost = sqrt(1.+OIII_RATIO*OIII_RATIO)/(1.+OIII_RATIO);
#endif
#ifdef USE_NII
    get_HaLF2(z,FHa,model,&NII_Ha);
    boost = get_SNR_boost_factor(z,reff,NII_Ha,psf);

    /* This shouldn't happen. This really shouldn't be <1, but the limit is set at 0.98
     * for numerical reasons.
     */
    if (boost<0.98) {
      fprintf(stderr, "Error: boost=%8.6lf\n", boost);
      fprintf(stderr, "Diagnostics: z=%6.4lf reff=%7.4lf NII_Ha=%7.4lf disp=%7.2lf as/um\n",
        z, reff, NII_Ha, psf->linDisp);
      exit(1);
    }
#endif
    sf = get_limflux(reff,psf,var_1exp,calib_1exp,EWD_REF,N_exp,significance_cut) / boost;
    x += sf > FHa? -dx: dx;
  }

  P = 0.5*(1.+geterf(x));
  if (x<=-4.0) P=0; /* If we can only see objects with -4 sigma size fluctuations,
                     * let's just call them nondetections.
                     */
  return(P);
}

/* Obtain the comoving density of detectable galaxies per ln Flux, at a given redshift and
 * flux level and ancillary parameters.
 *
 * Output units are Mpc^-3.
 * Input units: FHa (W/m2), var_1exp(e-^2/pix), calib_1exp(e- per W/m^2)
 *
 * Here "detectable" means above threshold -- does not count losses due to bright stars etc.
 *
 * If USE_SNRE is set, then we simply take the probability of SNRe being greater than the
 * threshold. Otherwise (DEFAULT) we want the probability of SNRo being greater than the
 * threshold -- this means a 2-point integration over the distribution of SNRe.
 */
double get_n_Detectable(double z, double FHa, PSF_DATA *psf, double var_1exp, double calib_1exp,
  int N_exp, double significance_cut, int model) {

#ifdef USE_SNRE
  return( get_P_Detectable(z,FHa,psf,var_1exp,calib_1exp,N_exp,significance_cut,model) * get_HaLF(z,FHa,model) );
#else
  /* DEFAULT */
  double P1,P2;
  P1 = get_P_Detectable(z,FHa,psf,var_1exp,calib_1exp,N_exp,significance_cut+1.-0.5/significance_cut,model);
  P2 = get_P_Detectable(z,FHa,psf,var_1exp,calib_1exp,N_exp,significance_cut-1.-0.5/significance_cut,model);
  return( 0.5*(P1+P2)*get_HaLF(z,FHa,model) );
#endif
}

/* Obtain the total comoving density of detectable galaxies, and (if requested) the
 * statistics of their flux distribution:
 *
 * stats[0] = geometric mean flux (W/m2)
 * stats[1] = sigma(ln F)
 * stats[2] = dimensionless skewness of ln F (=0 for lognormal)
 * stats[3] = dimensionless kurtosis of ln F (=0 for lognormal)
 *
 * Setting stats=NULL does not return the output.
 *
 * Output units are Mpc^-3.
 * Input units: FHa (W/m2), var_1exp(e-^2/pix), calib_1exp(e- per W/m^2)
 *  
 * Here "detectable" means above threshold -- does not count losses due to bright stars etc.
 */
double get_n_galaxies(double z, PSF_DATA *psf, double var_1exp, double calib_1exp,
  int N_exp, double significance_cut, int model, double *stats) {

#define NBIN 215
  int k;
  double ngal, thisn[NBIN];
  double Fmin, dlnF, r, dr;
  double rmean, rvar, rskew, rkurt;

  /* Get minimum detectable flux for a point source in W/m2, divided by 2 to account for
   * objects that up-scatter due to statistical fluctuations.
   */
  Fmin = get_1sigma_flux_1exp(0,psf,var_1exp,calib_1exp)
         * significance_cut / sqrt((double)N_exp) / 2.;

  /* Compute comoving density per flux bin up to 4 orders of magnitude brighter than min flux.
   * Track r = ln(F/Fmin) and compute the integral over all galaxies of r^j
   * (j=0..4).
   */
  dlnF = 0.02*log(10);
  for(k=0;k<NBIN;k++) {
    r = (k+0.5)*dlnF;
    thisn[k] = get_n_Detectable(z,Fmin*exp(r),psf,var_1exp,calib_1exp,N_exp,significance_cut,model)*dlnF;
  }

  /* Get comoving number of galaxies */
  ngal = 0.;
  for(k=0;k<NBIN;k++) ngal += thisn[k];

  /* If the stats are desired */
  if (stats!=NULL) {

    /* First get the mean value of r */
    rmean = 0.;
    for(k=0;k<NBIN;k++) {
      r = (k+0.5)*dlnF;
      rmean += r*thisn[k];
    }
    rmean/=ngal;
    stats[0] = Fmin*exp(rmean);

    /* Now go for the variance of r */
    rvar = 0.;
    for(k=0;k<NBIN;k++) {
      dr = (k+0.5)*dlnF - rmean;
      rvar += dr*dr*thisn[k];
    }
    rvar/=ngal;
    stats[1] = sqrt(rvar);

    /* The dimensionless skewness */
    rskew = 0.;
    for(k=0;k<NBIN;k++) {
      dr = (k+0.5)*dlnF - rmean;
      rskew += dr*dr*dr*thisn[k];
    }
    rskew/=ngal;
    stats[2] = rskew/(stats[1]*stats[1]*stats[1]);

    /* The dimensionless kurtosis */
    rkurt = 0.;
    for(k=0;k<NBIN;k++) {
      dr = (k+0.5)*dlnF - rmean;
      rkurt += dr*dr*dr*dr*thisn[k];
    }
    rkurt/=ngal;
    stats[3] = rkurt/(stats[1]*stats[1]*stats[1]*stats[1]) - 3.;
  }

  return(ngal);
#undef NBIN
}

/* --- ROUTINES TO COMPUTE DENSITIES OF OBSERVABLE GALAXIES: WL --- */

/* Get the number density (in gal/deg^2) of observable galaxies in the specified
 * effective radius and redshift ranges. These are [reffmin,reffmax) and [zmin,zmax)
 * so as to ensure that overlapping ranges are added together consistently.
 *
 * If dNeffdA!=NULL assigns the value of effective number density (measured relative
 * to intrinsic ellipticity dispersion 0.4 rms per axis).
 *
 * If OutFile!=NULL, sends a table of galaxies to the specified file in write or
 * append mode depending on whether outmode is 'w' or 'a'.
 */
double get_dNdA_WL_gal2(double reffmin, double reffmax, double zmin, double zmax, char GalaxyCat[],
  double lambda_min, double lambda_max, PSF_DATA *psf, int N_exp, double var_1exp,
  double calib_1exp, double max_ellip_err, double *dNeffdA, char OutFile[], char outmode) {

  /* Maximum number of wavelength points allowed in galaxy SEDs */
#define NWMAX 256
  /* Number of effective radius sub-bins */
#define NRBINS 100

  char InfoLine[512];
  FILE *fp, *fq;
  int i, iw, Nw;
  long id, igal, Ngal, count;
  double z, r_eff, area, cw[NWMAX], abmag[NWMAX];
  double lambda, galflux, galmag, frac;
  double count_eff, sig_e;
  double maglim[NRBINS], maglim_noSNRcut[NRBINS];

  /* Return 0 if ranges inconsistent */
  if (reffmin>=reffmax || zmin>=zmax) return(0);

  /* Generate limiting magnitudes in each range */
  for(i=0;i<NRBINS;i++) {
    maglim[i] = get_maglim_shape( reffmin*pow(reffmax/reffmin,(i+0.5)/(double)NRBINS),
                  psf, N_exp, var_1exp, calib_1exp, max_ellip_err);
    maglim_noSNRcut[i] = get_maglim_shape_noSNRcut( reffmin*pow(reffmax/reffmin,(i+0.5)/(double)NRBINS),
                  psf, N_exp, var_1exp, calib_1exp, max_ellip_err);
  }

  /* Open galaxy catalog, skip comment lines */
  fp = fopen(GalaxyCat, "r");
  if (fp==NULL) {
    fprintf(stderr, "Error: Can't read file: %s\n", GalaxyCat);
    wait_for_user();
    exit(1);
  }
  i=0;
  do {
    i++;
    if (fgets(InfoLine, 511, fp)==NULL) {
      fprintf(stderr, "Error while reading file: %s::%d\n", GalaxyCat, i);
      wait_for_user();
      exit(1);
    }
  } while (InfoLine[0]=='#');

  /* Basic data: number of wavelength points, number of galaxies, area of catalog in deg^2 */
  if (sscanf(InfoLine, "%d %ld %lg", &Nw, &Ngal, &area)!=3) {
      fprintf(stderr, "Error while reading file: %s::%d\n", GalaxyCat, i);
      wait_for_user();
      exit(1);
  }
  for(iw=0;iw<Nw;iw++) {
    if (fscanf(fp, "%lg", cw+iw)==EOF) {
      fprintf(stderr, "Error: unexpected EOF: %s in reference filter definition\n", GalaxyCat);
      wait_for_user();
      exit(1);
    }
  }

#if 0
  /* Test circuit: galaxy file header */
  fprintf(stderr, "Galaxy info: %d wavelengths, %ld galaxies, %12.5lE deg^2\n", Nw, Ngal, area);
  for(iw=0;iw<Nw;iw++) fprintf(stderr, "Wavelength #%2d = %12.5lE um\n", iw, cw[iw]);
  exit(1);
#endif

  /* Make sure we aren't out of range */
  if (lambda_min<cw[0] || lambda_max>cw[Nw-1]) {
    fprintf(stderr, "Error: bandpass outside of specified SED set.\n");
    wait_for_user();
    exit(1);
  }

  /* Open output file */
  fq = NULL; /* This line does nothing except suppress a spurious compiler warning */
  if (OutFile!=NULL) {
    if (outmode!='w' && outmode!='a') {
      fprintf(stderr, "Error: get_dNdA_WL_gal2: outmode=%c invalid\n", outmode);
      wait_for_user();
      exit(1);
    }
    fq = fopen(OutFile, &outmode);
    if (fq==NULL) {
      fprintf(stderr, "Error: can't write to %s\n", OutFile);
      wait_for_user();
      exit(1);
    }
  }

  /* Read galaxy table */
  count=0;
  count_eff=0;
  for(igal=0; igal<Ngal; igal++) {
    if (fscanf(fp, "%ld %lg %lg", &id, &z, &r_eff)==EOF) {
      fprintf(stderr, "Error: unexpected EOF: %s at galaxy %ld of %ld\n", GalaxyCat, igal, Ngal);
      wait_for_user();
      exit(1);
    }
    for(iw=0;iw<Nw;iw++) {
      if (fscanf(fp, "%lg", abmag+iw)==EOF) {
        fprintf(stderr, "Error: unexpected EOF: %s at galaxy %ld of %ld\n", GalaxyCat, igal, Ngal);
        wait_for_user();
        exit(1);
      }
    }

    /* Determine whether this galaxy passes cuts */
    if (z>=zmin && z<zmax && r_eff>=reffmin && r_eff<reffmax) {

      /* Get AB magnitude of this galaxy */
      galflux = 0.;
      for(i=0;i<10;i++) {
        lambda = lambda_min * pow(lambda_max/lambda_min, (i+0.5)/10.0);
        iw=0;
        while (iw<Nw-2 && cw[iw+1]<lambda) iw++;
        frac = log(lambda/cw[iw])/log(cw[iw+1]/cw[iw]);
        galflux += exp(-(abmag[iw]*(1-frac)+abmag[iw+1]*frac)/1.085736204758129);
      }
      galflux /= 10;
      galmag = -1.085736204758129*log(galflux);

      /* Get size bin */
      i = (int)floor(NRBINS*log(r_eff/reffmin)/log(reffmax/reffmin));
      if (i<0) i=0;
      if (i>=NRBINS-1) i=NRBINS-1;

#if 0
      /* Test circuit: magnitude interpolation */
      if (igal<64) {
        fprintf(stderr, "%5.3lf %6.3lf %6.3lf    ", z, r_eff, galmag);
        for(iw=0;iw<Nw;iw++) fprintf(stderr, " %6.3lf", abmag[iw]);
        fprintf(stderr, "\n");
      }
#endif

      /* Count the galaxy if it is bright enough */
      if (galmag < maglim[i]) {
        count++;
        sig_e = max_ellip_err * pow(10, -0.4*(maglim_noSNRcut[i]-galmag));
        count_eff += 1./(1. + sig_e*sig_e/0.16);
        if (OutFile!=NULL) {
#ifdef MY_WL_EXTRAFIELDS
          fprintf(fq, "%7ld %6.4lf %6.4lf %7.5lf %6.4lf %11.5lE\n", id, z, r_eff, sig_e, galmag, calib_1exp*pow(10,-.4*galmag));
#else
          fprintf(fq, "%7ld %6.4lf %6.4lf %7.5lf\n", id, z, r_eff, sig_e);
#endif
        }
      }

    }
  }
  fclose(fp);
  if (OutFile!=NULL) fclose(fq);

  /* Return effective number of galaxies per square degree */
  if (dNeffdA!=NULL) *dNeffdA = count_eff/area;

  /* Number of galaxies in this range per square degree */
  return(count/area);
#undef NWMAX
#undef NRBINS
}

/* Same as get_dNdA_WL_gal2 but no output file -- for backward compatibility */
double get_dNdA_WL_gal(double reffmin, double reffmax, double zmin, double zmax, char GalaxyCat[],
  double lambda_min, double lambda_max, PSF_DATA *psf, int N_exp, double var_1exp,
  double calib_1exp, double max_ellip_err, double *dNeffdA) {

  return (
    get_dNdA_WL_gal2(reffmin,reffmax,zmin,zmax,GalaxyCat,lambda_min,lambda_max,psf,N_exp,
    var_1exp,calib_1exp,max_ellip_err,dNeffdA,NULL,'a')
  );
}

/* --- ROUTINES TO COMPUTE THE FOREGROUND RADIATION AND ABSORPTION --- */

/* Ratio of extinction A_lambda to reddening E(B-V) for Milky Way dust.
 * Input lambda in microns.
 *
 * Valid range: 0.1--10.0 um
 *
 * Uses lookup table based on:
 * ftp://ftp.astro.princeton.edu/draine/dust/mix/kext_albedo_WD_MW_3.1_60_D03.all
 * Retrieved 03/18/2011
 */
double Galactic_Alambda__EBV(double lambda) {

  int il;
  double lset, lf;
  double norm[] = {
    0.24174,  0.25504,  0.26708,  0.26392,  0.23976,  0.20540,  0.17650,  0.15484,  0.13199,  0.10599,
    0.08604,  0.07340,  0.06302,  0.05489,  0.05107,  0.05005,  0.05041,  0.05151,  0.05313,  0.05554,
    0.06009,  0.06392,  0.06130,  0.06117,  0.06239,  0.06417,  0.06629,  0.06860,  0.07123,  0.07400,
    0.07702,  0.08018,  0.08361,  0.08723,  0.09105,  0.09506,  0.09934,  0.10388,  0.10862,  0.11356,
    0.11876,  0.12416,  0.12989,  0.13581,  0.14207,  0.14858,  0.15543,  0.16261,  0.17156,  0.17781,
    0.18571,  0.19401,  0.20276,  0.21192,  0.22140,  0.23127,  0.24147,  0.25214,  0.26313,  0.27459,
    0.28650,  0.29875,  0.31152,  0.32469,  0.33831,  0.35240,  0.36708,  0.38216,  0.39776,  0.41389,
    0.43055,  0.44779,  0.46557,  0.48394,  0.50296,  0.52265,  0.54292,  0.56386,  0.58552,  0.60790,
    0.63094,  0.65477,  0.67939,  0.70441,  0.73074,  0.75839,  0.78670,  0.81567,  0.84529,  0.87689,
    0.90915,  0.94273,  0.97762,  1.01382,  1.05201,  1.09151,  1.13232,  1.17512,  1.21922,  1.26531,
    1.31336,  1.36340,  1.41475,  1.46741,  1.52205,  1.57867,  1.63594,  1.69585,  1.75708,  1.82028,
    1.88545,  1.95260,  2.02107,  2.09217,  2.16524,  2.24029,  2.31731,  2.39631,  2.47795,  2.56155,
    2.64714,  2.73469,  2.82423,  2.91573,  3.00987,  3.10533,  3.20342,  3.30349,  3.40487,  3.50823,
    3.61290,  3.71955,  3.82752,  3.93680,  4.04608,  4.15668,  4.26728,  4.37854,  4.49111,  4.60369,
    4.71692,  4.83081,  4.94470,  5.05991,  5.17643,  5.29361,  5.41145,  5.53061,  5.65109,  5.77354,
    5.89862,  6.02567,  6.15339,  6.26596,  6.38907,  6.52732,  6.68203,  6.87294,  7.09019,  7.33377,
    7.64319,  7.96577,  8.45293,  8.95326,  9.47334,  9.94075, 10.22383, 10.19750,  9.92100,  9.52600,
    9.11126,  8.75576,  8.47268,  8.26201,  8.11718,  8.01843,  7.96577,  7.95260,  7.96577,  8.02502,
    8.11718,  8.24226,  8.39368,  8.55826,  8.74259,  8.92034,  9.10467,  9.31534,  9.54575,  9.82883,
   10.11192, 10.41475, 10.72416, 11.05332, 11.40224, 11.77749, 12.17248, 12.61356, 13.06781, 13.52205,
   13.92363};

  if (lambda<.1 || lambda>10) {
    fprintf(stderr, "Error: Wavelength lambda = %12.5lE out of range for dust extinction law.\n", lambda);
    wait_for_user();
    exit(1);
  }

  /* Find place to linearly interpolate */
  lset = 100./log(10.) * log(10./lambda);
  il = (int)floor(lset);
  if (il<0) il=0;
  if (il>199) il=199;
  lf = lset-il;

  return(norm[il] + lf*(norm[il+1]-norm[il]));
}

/* Computes the zodiacal foreground radiation in units of photons/m2/s/arcsec^2.
 *
 * Inputs:
 *   ecl_lat = ecliptic latitude (input in *degrees*)
 *   ecl_dlon = ecliptic longitude relative to Sun (input in *degrees*)
 *   lambda_min = min wavelength (input in *microns*)
 *   lambda_max = max wavlenegth (input in *microns*)
 *   T = throughput table (NULL for unit throughput)
 *
 * Caveats:
 *   No allowance is made for annual variations (due to tilt of Ecliptic relative to dust midplane)
 *   Range of valid wavelengths = 0.22--2.50 microns (in particlar: neglected thermal emission)
 */
double get_zodi_bkgnd(double ecl_lat, double ecl_dlon, double lambda_min, double lambda_max,
  THROUGHPUT_DATA *T) {

  double deg = M_PI/180.; /* degrees */
  double zodi_tot;
  double elon, z, sky05, fco;
  int ilat, ilon, ilam;
  double frlat, frlon, frlam;
  double lambda, dlambda, index_lambda, sun_spec;

  int ilambda, Nlambda = 8; /* number of integration points in wavelength */

  /* Sky brightness table: rows are varying ecliptic latitude, cols are varying longitude,
   * at the values shown.
   *
   * This is at 0.5um in units of 1e-8 W/m^2/sr/um.
   * Electronic version of Table 17 of Leinert (1997), except for placeholders (1's) at
   * elongation <15 degrees (where we should not use this routine anyway!).
   */
  int nlat = 11;
  int nlon = 19;
  double betaTable[] = {0,5,10,15,20,25,30,45,60,75,90};
  double dlonTable[] = {0,5,10,15,20,25,30,35,40,45,60,75,90,105,120,135,150,165,180};
  double skyTable[] = {
       1,    1,    1, 3140, 1610, 985, 640, 275, 150, 100, 77,
       1,    1,    1, 2940, 1540, 945, 625, 271, 150, 100, 77,
       1,    1, 4740, 2470, 1370, 865, 590, 264, 148, 100, 77,
   11500, 6780, 3440, 1860, 1110, 755, 525, 251, 146, 100, 77,
    6400, 4480, 2410, 1410,  910, 635, 454, 237, 141,  99, 77,
    3840, 2830, 1730, 1100,  749, 545, 410, 223, 136,  97, 77,
    2480, 1870, 1220,  845,  615, 467, 365, 207, 131,  95, 77,
    1650, 1270,  910,  680,  510, 397, 320, 193, 125,  93, 77,
    1180,  940,  700,  530,  416, 338, 282, 179, 120,  92, 77,
     910,  730,  555,  442,  356, 292, 250, 166, 116,  90, 77,
     505,  442,  352,  292,  243, 209, 183, 134, 104,  86, 77,
     338,  317,  269,  227,  196, 172, 151, 116,  93,  82, 77,
     259,  251,  225,  193,  166, 147, 132, 104,  86,  79, 77,
     212,  210,  197,  170,  150, 133, 119,  96,  82,  77, 77,
     188,  186,  177,  154,  138, 125, 113,  90,  77,  74, 77,
     179,  178,  166,  147,  134, 122, 110,  90,  77,  73, 77,
     179,  178,  165,  148,  137, 127, 116,  96,  79,  72, 77,
     196,  192,  179,  165,  151, 141, 131, 104,  82,  72, 77,
     230,  212,  195,  178,  163, 148, 134, 105,  83,  72, 77
  };

  /* Solar spectrum: in units of W/m^2/sr/um at log10(lambda/um) = -0.80(0.01)+0.40
   * Ref: Colina, Bohlin, Castelli 1996 AJ 112, 307
   *
   * V band (550 nm) is SolarSpec[54]
   */
  double SolarSpec[] = {
    1.87138e-01, 2.61360e-01, 4.08020e-01, 6.22197e-01, 9.02552e-01, 1.51036e+00, 2.25890e+00, 2.75901e+00, 4.03384e+00, 5.42817e+00, 
    7.26182e+00, 1.01910e+01, 2.01114e+01, 3.62121e+01, 4.31893e+01, 5.43904e+01, 4.91581e+01, 4.95091e+01, 4.95980e+01, 5.93722e+01, 
    5.27380e+01, 1.02502e+02, 1.62682e+02, 2.53618e+02, 2.01084e+02, 2.08273e+02, 4.05163e+02, 5.39830e+02, 5.31917e+02, 6.31200e+02, 
    7.06134e+02, 8.13653e+02, 1.00508e+03, 9.56536e+02, 9.50568e+02, 9.82400e+02, 1.06093e+03, 1.12669e+03, 1.09922e+03, 1.10224e+03, 
    1.36831e+03, 1.72189e+03, 1.74884e+03, 1.59871e+03, 1.74414e+03, 1.98823e+03, 2.02743e+03, 2.00367e+03, 2.03584e+03, 1.90296e+03, 
    1.93097e+03, 1.86594e+03, 1.86655e+03, 1.87957e+03, 1.87978e+03, 1.83915e+03, 1.84447e+03, 1.80371e+03, 1.76779e+03, 1.70796e+03, 
    1.66589e+03, 1.61456e+03, 1.53581e+03, 1.51269e+03, 1.44957e+03, 1.39215e+03, 1.34031e+03, 1.28981e+03, 1.24501e+03, 1.19548e+03, 
    1.15483e+03, 1.10546e+03, 1.06171e+03, 9.94579e+02, 9.54006e+02, 9.15287e+02, 8.63891e+02, 8.31183e+02, 7.95761e+02, 7.62568e+02, 
    7.27589e+02, 6.94643e+02, 6.60883e+02, 6.21830e+02, 5.83846e+02, 5.59624e+02, 5.34124e+02, 5.06171e+02, 4.80985e+02, 4.63139e+02, 
    4.39482e+02, 4.13122e+02, 3.94543e+02, 3.75591e+02, 3.56069e+02, 3.35294e+02, 3.16374e+02, 2.98712e+02, 2.82737e+02, 2.69581e+02, 
    2.49433e+02, 2.36936e+02, 2.21403e+02, 2.04770e+02, 1.87379e+02, 1.75880e+02, 1.60408e+02, 1.46210e+02, 1.36438e+02, 1.24412e+02, 
    1.16500e+02, 1.07324e+02, 9.89669e+01, 9.12134e+01, 8.28880e+01, 7.71064e+01, 7.06245e+01, 6.42367e+01, 5.87697e+01, 5.39387e+01, 
    4.98208e+01
  };

  /* Put longitude between 0 and 180 */
  ecl_dlon = fabs(ecl_dlon);
  ecl_dlon -= 360*floor(ecl_dlon/360);
  if (ecl_dlon>180) ecl_dlon = 360-ecl_dlon;
  if (ecl_dlon>180) ecl_dlon = 180;

  /* Set latitude to be positive */
  ecl_lat = fabs(ecl_lat);
  if (ecl_lat>90) ecl_lat = 90;

  /* Check wavelength ranges */
  if (lambda_min<0.22 || lambda_max>2.50) {
    fprintf(stderr, "Error: range lambda = %12.5lE .. %12.5lE microns out of range.\n", lambda_min, lambda_max);
    wait_for_user();
    exit(1);
  }

  /* Compute elongation (Sun angle). Complain if <15 degrees. */
  z = cos(ecl_lat*deg)*cos(ecl_dlon*deg);
  elon = z>=1? 0: z<=-1? 180: acos(z)/deg;
  if (elon<15) {
    fprintf(stderr, "Error: get_zodi_bkgnd: elongation = %12.5lE degrees out of valid range.\n", elon);
    wait_for_user();
    exit(1);
  }

  /* Compute sky brightness at 0.5um in units of 1e-8 W/m^2/sr/um.
   * Fit to Table 17 of Leinert (1997).
   */
  ilat=0; while(betaTable[ilat+1]<ecl_lat && ilat<nlat-2) ilat++;
  ilon=0; while(dlonTable[ilon+1]<ecl_dlon && ilon<nlon-2) ilon++;
  frlat = (ecl_lat-betaTable[ilat])/(betaTable[ilat+1]-betaTable[ilat]);
  frlon = (ecl_dlon-dlonTable[ilon])/(dlonTable[ilon+1]-dlonTable[ilon]);
  sky05 = exp(
           log(skyTable[ilat  +(ilon  )*nlat]) * (1.-frlat) * (1.-frlon)
          +log(skyTable[ilat  +(ilon+1)*nlat]) * (1.-frlat) * (   frlon)
          +log(skyTable[ilat+1+(ilon  )*nlat]) * (   frlat) * (1.-frlon)
          +log(skyTable[ilat+1+(ilon+1)*nlat]) * (   frlat) * (   frlon)
          );


  /* Integrate over wavelengths */
  zodi_tot = 0.;
  dlambda = (lambda_max-lambda_min)/(double)Nlambda;
  for(ilambda=0; ilambda<Nlambda; ilambda++) {

    lambda = lambda_min + (ilambda+0.5)/(double)Nlambda * (lambda_max-lambda_min);

    /* Solar spectrum at this wavelength: F_lambda/F_{0.5um} */
    index_lambda = 100*log(lambda)/log(10.) + 80;
    ilam = (int)floor(index_lambda);
    frlam = index_lambda - ilam;
    sun_spec = (SolarSpec[ilam] + frlam*(SolarSpec[ilam+1]-SolarSpec[ilam]))/SolarSpec[50];

    /* Color correction relative to solar */
    if (lambda>0.5) {
      fco = 1. + (elon>90? 0.6: elon<30? 0.8: 0.8-0.2*(elon-30)/60.)*log(lambda/0.5)/log(10.);
    } else {
      fco = 1. + (elon>90? 0.9: elon<30? 1.2: 1.2-0.3*(elon-30)/60.)*log(lambda/0.5)/log(10.);
    }

    /* The integral for the zodiacal foreground.
     * Here sky05*fco*sun_spec*dlambda is the power per unit area per unit solid angle in this
     * wavelength range (Units: 1e-8 W/m^2/sr).
     *
     * The conversion from 1e-8 W --> photons/sec is 5.03411747e10*lambda(um).
     */
    zodi_tot += sky05 * fco * sun_spec * dlambda * 5.03411747e10*lambda
                * get_throughput(lambda,T);

  }

  /* We now have the zodi level in photons/m^2/sr/sec. Convert to photons/m^2/arcsec^2/sec. */
  return(zodi_tot / 4.2545170296152206e10);
}

/* --- ROUTINES TO COMPUTE NOISE AND RELATED QUANTITIES --- */

/* Obtains the blackbody rate in units of photons per micron^2 per second per sr through the
 * filter defined by Thr in the wavelength range lambda_min .. lambda_max (in um) at
 * temperature T (in K).
 *
 * Currently uses 8192 linearly spaced steps across the waveband.
 */
double get_BB_flux_filter(THROUGHPUT_DATA *Thr, double lambda_min, double lambda_max,
  double T) {

  long i, N = 8192;
  double dl, lambda, x, Nlambda, f, totflux=0;

  dl = (lambda_max-lambda_min)/(double)N;
  for(i=0;i<N;i++) {
    lambda = lambda_min+dl*(i+0.5);
    x = 14387.75/lambda/T;

    /* blackbody photons/um^2/s/um/sr; N_lambda = 2*c/lambda^4/(e^x-1) */
    f = x>1e-4? exp(-x)/(1.-exp(-x)): 1./x/(1.+x/2.*(1.+x/3.*(1.+x/4.)));
    Nlambda = 5.99584916e14/lambda/lambda/lambda/lambda*f;

    totflux += Nlambda*dl*get_throughput(lambda,Thr);
  }
  return(totflux);
}

/* Computes thermal emission count rate in the detector due to the telescope only, as a
 * function of:
 *   Thr                throughput table for telescope, excluding filter
 *   lambda_min         minimum wavelength of filter (microns)
 *   lambda_max         maximum wavelength of filter (microns)
 *   filter_throughput  filter throughput in band
 *   Ttel               telescope temperature (PM + SM + associated structures)
 *   fratio             f/ratio of the telescope (open aperture only -- increased by outer mask)
 *   ismask             flag: 0 = no cold mask, 1 = with cold mask
 *   centobs            linear central obscuration (only used w/o mask)
 *   th_out_rad         outer radius of thermal emission annulus divided by outer radius
 *                        of primary mirror (>1; only used w/o mask)
 *   pmsm_emissivity    emissivity of PM + SM mirror surfaces
 *
 * Output is in e/pix/s, without margin.
 */
double getThermalBkgndTel(THROUGHPUT_DATA *Thr, double lambda_min, double lambda_max,
  double filter_throughput, double Ttel, double fratio, int ismask, double centobs,
  double th_out_rad, double pmsm_emissivity) {

  double bbrate, Omega;

  /* Blackbody count rate in photons/sr/pix/s at telescope temperature -- this
   * is key to normalizing the rest of the computation.
   *
   * The reflectivities of PM/SM are divided out since some structures are seen "directly"
   * after the TM/folds, and do not have to bounce off the PM/SM on their way to the
   * detector -- so we divide by the reflection coefficients.
   */
  bbrate = filter_throughput * PIXSIZE_PHYS * PIXSIZE_PHYS
           * get_BB_flux_filter(Thr, lambda_min, lambda_max, Ttel)
           / (1.-pmsm_emissivity) / (1.-pmsm_emissivity);

  /* Now we need to figure out the emissivity-weighted solid angle that the
   * detector actually sees. This is the sum of the open solid angle that sees
   * the sky (times the effective emissivity), plus the black surface contributions.
   */
  Omega = M_PI/4. / (fratio*fratio) * (
              (1.-centobs*centobs) * (2.-pmsm_emissivity) * pmsm_emissivity
            + (ismask? 0: centobs*centobs + th_out_rad*th_out_rad-1.)
          );

  return(bbrate*Omega); /* in photons/pix/sec */
}

/* Computes the total thermal background rate in the detector from all sources (telescope,
 * aft optics, and sources downstream of the filter). Output in e-/pix/s.
 *
 * If component==0 returns sum of all sources; for component==1, 2, 3, returns just the
 * specified source.
 */
double getThermalBkgndTot(THROUGHPUT_DATA *Thr, double lambda_min, double lambda_max,
  double filter_throughput, double fratio, double centobs, THERMAL_DATA *thermal,
  int component) {

  /* The three sources */
  double bkgnd1, bkgnd2, bkgnd3;

  /* Telescope */
  bkgnd1 = getThermalBkgndTel(Thr, lambda_min, lambda_max, filter_throughput,
           thermal->T_tel, fratio, thermal->ismask, centobs, thermal->outer_ann_ratio,
           thermal->pmsm_emissivity)
           * THERMAL_MARGIN_FACTOR;

  /* Aft optics */
  bkgnd2 = filter_throughput * PIXSIZE_PHYS * PIXSIZE_PHYS
           * get_BB_flux_filter(Thr, lambda_min, lambda_max, thermal->T_aft)
           * thermal->aft_net_emissivity/(1.-thermal->aft_net_emissivity)
           / (1.-thermal->pmsm_emissivity) / (1.-thermal->pmsm_emissivity)
           * (thermal->ismask? 1: thermal->outer_ann_ratio*thermal->outer_ann_ratio)
           * M_PI/4. / (fratio*fratio)
           * THERMAL_MARGIN_FACTOR;

  /* Everything else */
  bkgnd3 = thermal->post_filter_rate * THERMAL_MARGIN_FACTOR;

  /* Choose which source to return */
  if (component==1) return(bkgnd1);
  if (component==2) return(bkgnd2);
  if (component==3) return(bkgnd3);
  if (component==0) return(bkgnd1 + bkgnd2 + bkgnd3);

  fprintf(stderr, "Error: getThermalBkgndTot: invalid component number.\n");
  wait_for_user();
  exit(1);
}

/* Returns the noise variance (in counts) for an exposure with time
 * t_exp (s), with total number of photon counts ct (including dark current),
 * and read noise floor rnfloor (counts rms).
 */

double getNoiseTotal(double t_exp, double ct, double rnfloor, int mode) {
  double varPoisson, varRead, n, stretch;
  double m;

  switch(mode) {

    /* Sum of Poisson noise and read -- appropriate for e.g. CCDs */
    case 0:
      return(rnfloor*rnfloor+ct);
      break;

    /* SUTR, unweighted fit.
     * Rauscher et al 2007 PASP 119, 768 Eq. 1.
     * Floor is added in quadrature (mode 1) or as a replacement if the formula goes lower (mode 2).
     */
    case 1:
    case 2:
      n = floor(t_exp/T_FRAME);
      if (n<2) {
        fprintf(stderr, "Error: getNoiseTotal: SUTR illegal for exposure of <2 frames.\n");
        wait_for_user();
        exit(1);
      }
      stretch = t_exp/T_FRAME * t_exp/T_FRAME;
      varRead = 12./n/(n*n-1.)*VAR_READ_FRAME * stretch;
      varPoisson = 1.2/n*(n*n+1.)/(n*n-1.)*T_FRAME*ct/t_exp * stretch;
      if (mode==1)
        varRead += rnfloor*rnfloor;
      if (mode==2)
        varRead = varRead>rnfloor*rnfloor? varRead: rnfloor*rnfloor;     
      return(varRead + varPoisson);
      break;

    /* Fowler 4 (case 3) or Fowler 8 (case 4) */
    case 3:
      m=4;
    case 4:
      if (mode==4) m=8;
      n = floor(t_exp/T_FRAME);
      if (n<2*m) {
        fprintf(stderr, "Error: getNoiseTotal: Fowler-%d illegal for exposure of <%d frames.\n", (int)m, 2*(int)m);
        wait_for_user();
        exit(1);
      }
      stretch = t_exp/(n-m)/T_FRAME * t_exp/(n-m)/T_FRAME;
      varRead = 2./m*VAR_READ_FRAME * stretch;
      varPoisson = (n-m)*T_FRAME*ct/t_exp * stretch;
      varPoisson -= 1./3.*(m-1./m)*T_FRAME*ct/t_exp * stretch;
      varRead += rnfloor*rnfloor;
      return(varRead + varPoisson);
      break;

    /* Unrecognized mode */
    default:
      fprintf(stderr, "Error: getNoiseTotal: invalid read mode #%d.\n", mode);
      wait_for_user();
      exit(1);
  }

  /* Can't actually get here */
  return(0);
}

/* --- I/O AND ASSOCIATED ROUTINES --- */

/* Checks a thermal data structure for unusual or illegal values. */
void error_check_thermal(THERMAL_DATA *thermal) {
  /* Error checks */
  if (thermal->T_tel<=0) {
    fprintf(stderr, "Error: telescope temperature %12.5lE K is illegal.\n", thermal->T_tel);
    wait_for_user();
    exit(1);
  }
  if (thermal->T_tel<160 || thermal->T_tel>300) {
    fprintf(stderr, "Warning: telescope temperature %12.5lE K is unusual.\n", thermal->T_tel);
  }
  if (thermal->pmsm_emissivity<0 || thermal->pmsm_emissivity>=1) {
    fprintf(stderr, "Error: PM/SM emissivity %12.5lE is illegal.\n", thermal->pmsm_emissivity);
    wait_for_user();
    exit(1);
  }
  if (thermal->pmsm_emissivity<0.001 || thermal->pmsm_emissivity>0.1) {
    fprintf(stderr, "Warning: PM/SM emissivity %12.5lE is unusual.\n", thermal->pmsm_emissivity);
  }
  if (thermal->outer_ann_ratio<0.999999) {
    fprintf(stderr, "Error: outer annulus ratio rho = %12.5lE is illegal.\n", thermal->outer_ann_ratio);
    wait_for_user();
    exit(1);
  }
  if (thermal->outer_ann_ratio>1.25) {
    fprintf(stderr, "Warning: outer annulus ratio rho = %12.5lE is unusual.\n", thermal->outer_ann_ratio);
  }
  if (thermal->ismask!=0 && thermal->ismask!=1) {
    fprintf(stderr, "Error: cold mask is binary (0/1); %d is illegal.\n", thermal->ismask);
    wait_for_user();
    exit(1);
  }
  if (thermal->T_aft<=0) {
    fprintf(stderr, "Error: aft optics temperature %12.5lE K is illegal.\n", thermal->T_aft);
    wait_for_user();
    exit(1);
  }
  if (thermal->T_aft<140 || thermal->T_aft>250) {
    fprintf(stderr, "Warning: aft optics temperature %12.5lE K is unusual.\n", thermal->T_aft);
  }
  if (thermal->aft_net_emissivity<0 || thermal->aft_net_emissivity>=1) {
    fprintf(stderr, "Error: PM/SM emissivity %12.5lE is illegal.\n", thermal->pmsm_emissivity);
    wait_for_user();
    exit(1);
  }
  if (thermal->post_filter_rate<=0) {
    fprintf(stderr, "Error: post-filter thermal background %12.5lE e-/pix/s is illegal.\n",
      thermal->post_filter_rate);
    wait_for_user();
    exit(1);
  }
  if (thermal->post_filter_rate>0.05) {
    fprintf(stderr, "Warning: post-filter thermal background %12.5lE e-/pix/s is unusual.\n",
      thermal->post_filter_rate);
  }
}

/* Reads thermal data from a character string.
 * Format: 7 numbers, including:
 *   T_tel pmsm_emissivity outer_ann_ratio ismask T_aft aft_net_emissivity post_filter_rate
 * Performs basic error checks.
 */
void read_thermal_data(char Data[], THERMAL_DATA *thermal) {
  if (sscanf(Data, "%lg %lg %lg %d %lg %lg %lg",
    &(thermal->T_tel),
    &(thermal->pmsm_emissivity),
    &(thermal->outer_ann_ratio),
    &(thermal->ismask),
    &(thermal->T_aft),
    &(thermal->aft_net_emissivity),
    &(thermal->post_filter_rate)
  )!=7) {
    fprintf(stderr, "Error: read_thermal_data: failed call to sscanf\n");
    wait_for_user();
    exit(1);
  }
  error_check_thermal(thermal);
}

/* === MAIN PROGRAM: INCLUDES DRIVER FOR OTHER ROUTINES === */

int main(void) {
  FILE *fp;
  char FileName[256], InfoLine[256];
  int N_exp, i, N, det_type;
  long tel_config;
  THROUGHPUT_DATA TAll, T1;
  double lambda_min = 0.22, lambda_max = 2.5;
  double lambda, D, centobs, jitter, pixscale, fratio;
  double ecl_lat, ecl_dlon, zodi_bkgnd, tot_bkgnd, gal_ebv;
  double throughput, t_exp, readnoise, darkcurrent, thermalbkgnd[4];
  double var_1exp, var_1exp_skyonly, var_1exp_bkgnd, calib_1exp, rmswfe_um;
  PSF_DATA psf;
  THERMAL_DATA thermal;
  int sample_mode;
  WAVEFRONT_ERROR_FUNCTION wavefront;
  int is_use_wavefront = 0;
#ifdef BAO_MODE
  double significance_cut, completeness, linDisp;
  int gal_pop_model;
#ifndef PZCAL_MODE
  double z;
  double nbar, stats[4], dNdA, dz;
#endif
#endif
#ifdef TILT_BACK
  double epsilon;
#endif
  double filterthroughput = 1.;
#ifdef WL_MODE
  int iz;
  char GalaxyCat[256];
  double cpp, ee50;
  double min_res_factor, max_ellip_err;
  double min_reff, r_eff;
  double dn, ntot, dneff, nefftot;
#endif
#ifdef OUT_WL_CAT
  char OutGalCat[256];
#endif
#ifdef PZCAL_MODE
  int nstep;
#endif
#ifdef SPCONT_MODE
  double spDispersion, thetaeff, r_eff, flux1, mag1;
  int nstep;
  double npix_2d, dll_1pix, l1, l2, varperpixel_1exp; /* IFU parameters */
#endif
  int is_an_ifu = 0;

  /* --- REPORT CODE CONFIGURATION --- */
  printf("Exposure time calculator v19\n");
#ifdef WL_MODE
  printf("Mode: WL");
#endif
#ifdef SPCONT_MODE
  printf("Mode: SPCONT");
#endif
#ifdef BAO_MODE
#ifdef PZCAL_MODE
  printf("Mode: PZCAL");
#else
  printf("Mode: BAO");
#endif
#endif
  printf("\nOptions: ");
#ifdef IN_DLON
  printf("-DIN_DLON ");
#endif
#ifdef TILT_BACK
  printf("-DTILT_BACK ");
#endif
#ifdef KEEP_WINDOW_OPEN
  printf("-DKEEP_WINDOW_OPEN ");
#endif
#ifdef READ5
  printf("-DREAD5 ");
#endif
#ifdef OUT_WL_CAT
  printf("-DOUT_WL_CAT ");
#endif
#ifdef LOGICAL_READ_FLOOR
  printf("-DLOGICAL_READ_FLOOR ");
#endif
#ifdef WLCUT_DEFAULT
  printf("-DWLCUT_DEFAULT ");
#endif
#ifdef WFE_OVERRIDE
  printf("-DWFE_OVERRIDE ");
#endif
#ifdef WFE_HOPKINS_EXP
  printf("-DWFE_HOPKINS_EXP ");
#endif
#ifdef NO_OBJ_CONTINUUM
  printf("-DNO_OBJ_CONTINUUM ");
#endif
#ifdef USE_SNRE
  printf("-DUSE_SNRE ");
#endif
#ifdef OUT_EXTRA_PSF_PROPERTIES
  printf("-DOUT_EXTRA_PSF_PROPERTIES ");
#endif
#ifdef USE_NII
  printf("-DUSE_NII ");
#endif
#ifdef OIII_GAL
  printf("-DOIII_GAL ");
#endif
#ifdef OII_GAL
  printf("-DOII_GAL ");
#endif
#ifdef FOWLER4
  printf("-DFOWLER4 ");
#endif
#ifdef FOWLER8
  printf("-DFOWLER8 ");
#endif
  printf("\n");

  /* ... and complain about incompatible options */
#ifdef IN_DLON
#ifdef TILT_BACK
  fprintf(stderr, "Error: IN_DLON and TILT_BACK not compatible.\n");
  wait_for_user();
  exit(1);
#endif
#endif

  /* --- SET DEFAULT VALUES OF THERMAL PARAMETERS (CAN BE CHANGED LATER) --- */

  thermal.T_tel               = 230;   /* K */
  thermal.pmsm_emissivity     = 0.02;
  thermal.outer_ann_ratio     = 1.05;
  thermal.ismask              = 0;
  thermal.T_aft               = 230;   /* K */
  thermal.aft_net_emissivity  = 0.0776;
#ifdef WL_MODE
  thermal.post_filter_rate    = 0.005;  /* counts/pix/s */
#endif
#ifdef BAO_OR_SPCONT_MODE
  thermal.post_filter_rate    = 0.03;  /* counts/pix/s */
#endif

  /* --- CONFIGURATION INPUT --- */

  printf("Enter telescope configuration [0=generic, 1=from file]: ");
  scanf("%ld", &tel_config);
  if (tel_config==0) {

    /* Generic case */
    printf("Enter aperture outer diameter (meters): ");
    scanf("%lg", &D);
    if (D<0) {
      fprintf(stderr, "Error: aperture outer diameter D=%12.5lE m is illegal.\n", D);
      wait_for_user();
      exit(1);
    }
    if (D<0.5 || D>2.5) {
      fprintf(stderr, "Warning: aperture outer diameter D=%12.5lE m is unusual.\n", D);
    }
    printf("Enter central obscuration (fractional - linear): ");
    scanf("%lg", &centobs);
    if (centobs<0 || centobs>=1) {
      fprintf(stderr, "Error: central obscuration = %12.5lE is illegal.\n", centobs);
      wait_for_user();
      exit(1);
    }
    printf("Enter pixel scale (arcsec): ");
    scanf("%lg", &pixscale);
    if (pixscale<=0) {
      fprintf(stderr, "Error: pixel scale = %12.5lE arcsec is illegal.\n", pixscale);
      wait_for_user();
      exit(1);
    }
    if (pixscale<0.03 || pixscale>1) {
      fprintf(stderr, "Warning: pixel scale = %12.5lE arcsec is unusual.\n", pixscale);
    }
#ifdef BAO_OR_SPCONT_MODE
    printf("Enter throughput (total, all orders): ");
#endif
#ifdef WL_MODE
    printf("Enter throughput: ");
#endif
    scanf("%lg", &throughput);
    if (throughput<0 || throughput>1) {
      fprintf(stderr, "Error: throughput = %12.5lE is illegal.\n", throughput);
      wait_for_user();
      exit(1);
    }
    if (throughput<.1 || throughput>.9) {
      fprintf(stderr, "Warning: throughput = %12.5lE is unusual.\n", throughput);
    }
    TAll.N=1;
    TAll.lambda[0] = lambda_min;
    TAll.throughput[0] = throughput;
#ifdef BAO_OR_SPCONT_MODE
    printf("Enter throughput (1st order): ");
    scanf("%lg", &throughput);
    if (throughput<0 || throughput>1) {
      fprintf(stderr, "Error: throughput = %12.5lE is illegal.\n", throughput);
      wait_for_user();
      exit(1);
    }
    if (throughput<.1 || throughput>.9) {
      fprintf(stderr, "Warning: throughput = %12.5lE is unusual.\n", throughput);
    }
#endif
    T1.N=1;
    T1.lambda[0] = lambda_min;
    T1.throughput[0] = throughput;
#ifdef BAO_OR_SPCONT_MODE
    if (T1.throughput[0]>TAll.throughput[0]) {
      fprintf(stderr, "Error: 1st order throughput cannot exceed all orders combined.\n");
      wait_for_user();
      exit(1);
    }
#endif
  }
  else {

    /* Read configuration from a file. First find the file. */
    printf("Input file name: ");
    if (!scanf("%255s", FileName)) {
      fprintf(stderr, "Error: Can't get input file.\n");
      wait_for_user();
      exit(1);
    }
    fp = fopen(FileName, "r");
    if (fp==NULL) {
      fprintf(stderr, "Error: Can't read file: %s\n", FileName);
      wait_for_user();
      exit(1);
    }

    /* File format:
     * For BAO the desired order is 1st, and we want the light in the central
     * spot; for WL all light is 0th. However scattered light (which contributes
     * to background but not shape measurement) is included in the 2nd column
     * but not the 3rd column.
     *
     * # Any comment lines start with # and are at the beginning.
     * configuration_name
     * D centobs pixscale N
     * lambda[0] throughput_all[0] throughput_desired_order[0]
     * ...
     * lambda[N-1] throughput_all[N-1] throughput_desired_order[N-1]
     *
     * The preamble may contain directives that start with !
     * See the manual for explicit information. The currently implemented ones are:
     * !THERMAL = specifies thermal data
     * !WAVEFRONT = wavefront error table in configuration file
     * !IMAGING = configuration file for imaging mode
     * !SLITLESS = configuration file for slitless spectroscopy
     * !IFU = configuration file for an IFU
     */

    /* Strip comments from the file; exit at first line that is not a comment. */
    i=0;
    do {
      i++;
      if (fgets(InfoLine, 511, fp)==NULL) {
        fprintf(stderr, "Error while reading file: %s::%d\n", FileName, i);
        wait_for_user();
        exit(1);
      }
      if (memcmp(InfoLine, "!THERMAL", (size_t)8)==0) {
        read_thermal_data(InfoLine+8, &thermal);
      }
      if (memcmp(InfoLine, "!WAVEFRONT", (size_t)10)==0) {
        is_use_wavefront = 1;
        sscanf(InfoLine+10, "%lg %lg %lg %lg", &(wavefront.lambdabest), &(wavefront.rmswfebest), &(wavefront.slope1), &(wavefront.slope2));
#ifndef BAO_OR_SPCONT_MODE
        fprintf(stderr, "Error: wavelength dependent wavefront error enabled only in spectroscopy modes at present.\n");
        exit(1);
#endif
      }
#ifdef BAO_OR_SPCONT_MODE
      if (memcmp(InfoLine, "!IMAGING", (size_t)8)==0) {
        fprintf(stderr, "Error: This is an imaging configuration file: not valid in a spectroscopic mode.\n");
        exit(1);
      }
#endif
#ifdef WL_MODE
      if (memcmp(InfoLine, "!SLITLESS", (size_t)9)==0) {
        fprintf(stderr, "Error: This is a slitless spectroscopy configuration file: not valid in an imaging mode.\n");
        exit(1);
      }
#endif
      if (memcmp(InfoLine, "!IFU", (size_t)4)==0) {
        is_an_ifu = 1;
#ifndef SPCONT_MODE
        fprintf(stderr, "Error: IFU currently supported only in SPCONT mode.\n");
        exit(1);
#endif
      }
    } while (InfoLine[0]=='#' || InfoLine[0]=='!');

    /* Output the configuration name */
    printf("Using configuration: %s", InfoLine);

    /* Read basic data on this configuration */
    if (fscanf(fp, "%lg %lg %lg %d", &D, &centobs, &pixscale, &N)==EOF) {
      fprintf(stderr, "Error: Unexpected end of file: %s\n", FileName);
      wait_for_user();
      exit(1);
    }
    TAll.N = T1.N = N;
    if (N>N_THROUGHPUT_MAX) {
      fprintf(stderr, "Error: %s: number of throughput grid points %d is too many (max=%d).\n", FileName, N, N_THROUGHPUT_MAX);
      wait_for_user();
      exit(1);
    }

    /* Get throughput table */
    for(i=0;i<N;i++) {
      if (fscanf(fp, "%lg %lg %lg", TAll.lambda+i, TAll.throughput+i, T1.throughput+i)==EOF) {
        fprintf(stderr, "Error: Unexpected end of file: %s\n", FileName);
        wait_for_user();
        exit(1);
      }
      T1.lambda[i] = TAll.lambda[i];
      if (i>0)
        if (T1.lambda[i]<=T1.lambda[i-1]) {
          fprintf(stderr, "Error: wavelengths %d,%d out of order in %s.\n", i-1, i, FileName);
          wait_for_user();
          exit(1);
        }
#ifdef MULT11
      TAll.lambda[i] *= 1.1;
#endif
    }

    /* End read file */
    fclose(fp);
  }

  /* Wave front error -- if the user needs to input it! */
#ifdef WFE_INPUT
  rmswfe_um = 0.;
  if (!is_use_wavefront) {
    printf("Enter rms wavefront error (microns): ");
    scanf("%lg", &rmswfe_um);
    if (rmswfe_um<0) {
      fprintf(stderr, "Error: rms wavefront error = %12.5lE microns is illegal.\n", rmswfe_um);
      wait_for_user();
      exit(1);
    }
  }
#ifdef WL_MODE
  if (rmswfe_um>0.12) {
#endif
#ifdef BAO_OR_SPCONT_MODE
  if (rmswfe_um>0.24) {
#endif
    fprintf(stderr, "Warning: rms wavefront error = %12.5lE microns is unusual.\n", rmswfe_um);
  }
#else
  rmswfe_um = 0.;
#endif

  /* Detector type */
  printf("Enter detector type [e.g. 0=H2RG]: ");
  scanf("%d", &det_type);
  ConfigDetector(det_type);

  printf("Enter pointing jitter (arcsec rms per axis): ");
  scanf("%lg", &jitter);
  if (jitter<0) {
    fprintf(stderr, "Error: jitter = %12.5lE is illegal.\n", jitter);
    wait_for_user();
    exit(1);
  }
  if (jitter>1) {
    fprintf(stderr, "Warning: jitter = %12.5lE arcsec rms is unusual.\n", jitter);
  }
#ifdef BAO_MODE
  jitter = sqrt(jitter*jitter + SIGMA_PROF*SIGMA_PROF);
#endif
  printf("Enter minimum wavelength (microns): ");
  scanf("%lg", &lambda_min);
  if (lambda_min<=0.22) {
    fprintf(stderr, "Error: lambda_min = %12.5lE micron is outside the range of validity of background model.\n", lambda_min);
    wait_for_user();
    exit(1);
  }
#ifdef BAO_MODE
#ifndef PZCAL_MODE
  if (lambda_min<LAMBDA0_HA) {
    fprintf(stderr, "Warning: lambda_min = %12.5lE micron is unusual.\n", lambda_min);
  }
#endif
#endif
  printf("Enter maximum wavelength (microns): ");
  scanf("%lg", &lambda_max);
  if (lambda_max<=lambda_min) {
    fprintf(stderr, "Error: lambda_max must be greater than lambda_min: (%12.5lE,%12.5lE) micron is illegal.\n", lambda_min, lambda_max);
    wait_for_user();
    exit(1);
  }
  if (lambda_max>2.5) {
    fprintf(stderr, "Error: lambda_max = %12.5lE micron is outside the range of validity of background model.\n", lambda_max);
    wait_for_user();
    exit(1);
  }
#ifdef WL_MODE
  printf("Enter filter throughput: ");
  scanf("%lg", &filterthroughput);
  if (filterthroughput<=0 || filterthroughput>1.00000000001) {
    fprintf(stderr, "Error: filter throughput %12.5lE is illegal.\n", filterthroughput);
    wait_for_user();
    exit(1);   
  }
#endif
#ifdef SPCONT_MODE
  printf("Enter filter throughput: ");
  scanf("%lg", &filterthroughput);
  if (filterthroughput<=0 || filterthroughput>1.00000000001) {
    fprintf(stderr, "Error: filter throughput %12.5lE is illegal.\n", filterthroughput);
    wait_for_user();
    exit(1);   
  }
  if (filterthroughput<0.9 || filterthroughput>0.99) {
    fprintf(stderr, "Warning: filter throughput %12.5lE is unusual.\n", filterthroughput);
  }
#endif
#ifdef BAO_MODE
  printf("Enter completeness: ");
  scanf("%lg", &completeness);
  if (completeness<=0 || completeness>1.00000000001) {
    fprintf(stderr, "Error: completeness %12.5lE is illegal.\n", completeness);
    wait_for_user();
    exit(1);
  }
  printf("Enter linear spectral dispersion d theta/d lambda (arcsec/um): ");
  scanf("%lg", &linDisp);
  if (linDisp<=0) {
    fprintf(stderr, "Error: dispersion %12.5lE is illegal.\n", linDisp);
    wait_for_user();
    exit(1);
  }
  if (linDisp<1 || linDisp>1000) {
    fprintf(stderr, "Warning: dispersion %12.5lE is unusual.\n", linDisp);
  }
#endif
#ifdef SPCONT_MODE
  printf("Enter spectral dispersion, Dtheta (arcsec): ");
  scanf("%lg", &spDispersion);
  if (spDispersion<=0) {
    fprintf(stderr, "Error: dispersion %12.5lE is illegal.\n", spDispersion);
    wait_for_user();
    exit(1);
  }
  if (spDispersion<1 || spDispersion>1000) {
    fprintf(stderr, "Warning: dispersion %12.5lE is unusual.\n", spDispersion);
  }
  printf("Enter source effective radius (arcsec): ");
  scanf("%lg", &r_eff);
  if (r_eff<-1e-9) {
    fprintf(stderr, "Error: source effective radius %12.5lE is illegal.\n", r_eff);
    wait_for_user();
    exit(1);
  }
  if (r_eff>10) {
    fprintf(stderr, "Error: source effective radius %12.5lE is unusual.\n", r_eff);
  }
#endif
  printf("Enter single exposure time (s): ");
  scanf("%lg", &t_exp);
  if (t_exp<=0) {
    fprintf(stderr, "Error: exposure time = %12.5lE sec is illegal.\n", t_exp);
    wait_for_user();
    exit(1);
  }
  if (t_exp<20 || t_exp>2000) {
    fprintf(stderr, "Warning: exposure time = %12.5lE sec is unusual.\n", t_exp);
  }
#ifdef READ5
  readnoise = 5.;
  darkcurrent = 0.05;
#else
  printf("Enter read noise floor (effective e- rms per pixel): ");
  scanf("%lg", &readnoise);
  if (readnoise<0) {
    fprintf(stderr, "Error: read noise = %12.5lE e- rms is illegal.\n", readnoise);
    wait_for_user();
    exit(1);
  }
  if (readnoise<2 || readnoise>30) {
    fprintf(stderr, "Warning: read noise = %12.5lE e- rms is unusual.\n", readnoise);
  }
  printf("Enter dark current (e-/pix/sec): ");
  scanf("%lg", &darkcurrent);
  if (darkcurrent<0) {
    fprintf(stderr, "Error: dark current = %12.5lE e-/pix/s is illegal.\n", darkcurrent);
    wait_for_user();
    exit(1);
  }
  if (darkcurrent<1e-4 || darkcurrent>.1) {
    fprintf(stderr, "Warning: dark current = %12.5lE e-/pix/s is unusual.\n", darkcurrent);
  }
#endif
  printf("Enter ecliptic latitude (degrees): ");
  scanf("%lg", &ecl_lat);
  if (fabs(ecl_lat)>90.000000001) {
    fprintf(stderr, "Error: ecliptic latitude = %12.5lE degrees is illegal.\n", ecl_lat);
    wait_for_user();
    exit(1);
  }
#ifdef IN_DLON
  printf("Enter ecliptic longitude relative to the Sun (degrees): ");
  scanf("%lg", &ecl_dlon);
  ecl_dlon = fabs(ecl_dlon);
  ecl_dlon -= 360.*floor(ecl_dlon/360. + 0.5);
  ecl_dlon = fabs(ecl_dlon);
#endif
  printf("Enter Galactic reddening, E(B-V) (magnitudes): ");
  scanf("%lg", &gal_ebv);
  if (gal_ebv<0) {
    fprintf(stderr, "Error: reddening = %12.5lE magnitudes is illegal.\n", gal_ebv);
    wait_for_user();
    exit(1);
  }
  if (gal_ebv>1.) {
    fprintf(stderr, "Warning: reddening = %12.5lE magnitudes is unusual.\n", gal_ebv);
    wait_for_user();
    exit(1);
  }
  printf("Enter number of exposures: ");
  scanf("%d", &N_exp);
  if (N_exp<=0) {
    fprintf(stderr, "Error: number of exposures = %6d is illegal.\n", N_exp);
    wait_for_user();
    exit(1);
  }
  if (N_exp>100) {
    fprintf(stderr, "Warning: number of exposures = %6d is unusual.\n", N_exp);
  }

#ifdef BAO_MODE
  printf("Enter significance cut (number of sigmas): ");
  scanf("%lg", &significance_cut);
  if (significance_cut<=0) {
    fprintf(stderr, "Error: significance cut at %12.5lE sigma is illegal.\n", significance_cut);
    wait_for_user();
    exit(1);
  }
#ifndef USE_SNRE
  if (significance_cut<3.5) {
    fprintf(stderr, "Error: significance cut at %12.5lE sigma is outside range of validity of", significance_cut);
    fprintf(stderr, " approximations used for SNRo vs SNRe scatter in this code.\n");
    wait_for_user();
    exit(1);
  }
#endif
  if (significance_cut<5) {
    fprintf(stderr, "Warning: using galaxies with %12.5lE sigma detections not recommended.\n", significance_cut);
  }
#ifdef PZCAL_MODE
  gal_pop_model = 0.; /* NOT actually used */
#else
  printf("Enter galaxy population model [SDT2013 report = 42]: ");
  scanf("%d", &gal_pop_model);
  printf("Using luminosity function model: ");
  print_HaLF_model(stdout, gal_pop_model);
#endif
  printf("\n");
#endif

#ifdef WL_MODE
#ifdef WLCUT_DEFAULT
  min_res_factor = 0.4;
  max_ellip_err = 0.2;
#else
#ifdef WFE_OVERRIDE
  min_res_factor = 0.4;
  max_ellip_err = 0.2;
#else
  printf("Enter minimum resolution factor: ");
  scanf("%lg", &min_res_factor);
  if (min_res_factor<=0 || min_res_factor>=1) {
    fprintf(stderr, "Error: minimum resolution factor %12.5lE is illegal.\n", min_res_factor);
    wait_for_user();
    exit(1);
  }
  if (min_res_factor<.25 || min_res_factor>.6) {
    fprintf(stderr, "Warning: minimum resolution factor %12.5lE is unusual.\n", min_res_factor);
  }
  printf("Enter maximum ellipticity error: ");
  scanf("%lg", &max_ellip_err);
  if (max_ellip_err<0 || max_ellip_err>1) {
    fprintf(stderr, "Error: maximum ellipticity error %12.5lE is illegal.\n", max_ellip_err);
    wait_for_user();
    exit(1);
  }
  if (max_ellip_err<.1 || max_ellip_err>.25) {
    fprintf(stderr, "Warning: maximum ellipticity error %12.5lE is unusual.\n", max_ellip_err);
  }
#endif
#endif
#ifndef WFE_OVERRIDE
  printf("Enter input galaxy catalog file: ");
  if (!scanf("%255s", GalaxyCat)) {
    fprintf(stderr, "Error: Can't get input file.\n");
    wait_for_user();
    exit(1);
  }
#endif
#ifdef OUT_WL_CAT
  printf("Enter output galaxy catalog file: ");
  if (!scanf("%255s", OutGalCat)) {
    fprintf(stderr, "Error: Can't get output file.\n");
    wait_for_user();
    exit(1);
  }
#endif
#endif

  /* --- (CURRENTLY) FIXED PARAMETERS --- */

#ifdef PZCAL_MODE
  /* Number of output steps in wavelength */
  nstep = 28;
#endif

  /* Longitude relative to the Sun -- right now assume a Sun angle of 90 degrees */
#ifndef IN_DLON
  ecl_dlon = 90.;
#endif

#ifdef TILT_BACK
  ecl_lat *= 1-1e-10;
  epsilon = 180. - atan(2*tan(ecl_lat*M_PI/180.))*180./M_PI;
  if (epsilon>115) epsilon=115;
  ecl_dlon = 180./M_PI*acos(cos(epsilon*M_PI/180.)/cos(ecl_lat*M_PI/180.));
  printf("Tilting back from Sun: |lambda-lambda(Sun)| = %7.3lf deg, elongation = %7.3lf deg\n", ecl_dlon, epsilon);
#endif

  /* --- RECOMPUTE THROUGHPUTS INCLUDING FILTER FOR IMAGING --- */
  for(i=0;i<TAll.N;i++) TAll.throughput[i] *= filterthroughput;
  for(i=0;i<T1.N;i++) T1.throughput[i] *= filterthroughput;

  /* --- GENERIC ETC FUNCTIONS --- */

  /* The background is computed here. It is final for all ETC modes *except* the IFU, which has to be
   * re-computed inside the lambda loop.
   */

  /* Determine zodiacal background within band, in photons/m^2/arcsec^2/sec */
  zodi_bkgnd = get_zodi_bkgnd(ecl_lat, ecl_dlon, lambda_min, lambda_max, &TAll);

  /* The zodi is the dominant sky background, and the only one considered here. */
  tot_bkgnd = zodi_bkgnd;

  /* Sampling mode */
  sample_mode = 1;      /* Default */
#ifdef LOGICAL_READ_FLOOR
  sample_mode = 2;
#endif
#ifdef FOWLER4
  sample_mode = 3;
#endif
#ifdef FOWLER8
  sample_mode = 4;
#endif
  if (DET_IS_CCD) sample_mode = 0;

  /* Thermal background:
   * f/ratio requires conversions: () = pixel scale in microradians
   */
  fratio = PIXSIZE_PHYS/(pixscale*M_PI/0.648)/D;
  for(i=0; i<4; i++)
    thermalbkgnd[i] = getThermalBkgndTot(&TAll, lambda_min, lambda_max, filterthroughput,
                      fratio, centobs, &thermal, i);

  /* Noise variance per unit area in 1 exposure in e-^2/arcsec^2 */
  var_1exp_skyonly = tot_bkgnd * M_PI/4.* D * D * (1.-centobs*centobs) * t_exp;
  var_1exp_bkgnd = var_1exp_skyonly
             + (darkcurrent+thermalbkgnd[0])*t_exp/(pixscale*pixscale);
  var_1exp = getNoiseTotal(t_exp,var_1exp_bkgnd*pixscale*pixscale,readnoise,sample_mode) / (pixscale*pixscale);

  /* --- BAO ETC --- */

#ifdef BAO_MODE
#ifdef PZCAL_MODE

  /* Photo-z calibration mode simply displays line fluxes */

  printf("Limiting fluxes in W/m2 vs wavelength and galaxy size\n");
  printf("lambda|");
  for(i=0;i<=10;i++) printf("   %4.1lf\"   |", 0.1*i);
  printf("\n");
  printf("  um  |");
  for(i=0;i<=10;i++) printf("           |");
  printf("\n");
  for(lambda=lambda_min; lambda<1.0000001*lambda_max; lambda+=(lambda_max-lambda_min)/(double)nstep) {

    /* Setup PSF parameters: */
    /* Linear dispersion */
    psf.linDisp = 0.;
#ifdef BAO_MODE
    psf.linDisp = linDisp;
#endif
    /* Pixel scale, in arcsec */
    psf.pixscale = pixscale;
    /* 1 sigma jitter + diffusion, in arcsec */
    psf.sigma = sqrt((PIX_CD/PIXSIZE_PHYS*pixscale)*(PIX_CD/PIXSIZE_PHYS*pixscale) + jitter*jitter);
    /* lambda/D in arcsec -- hence conversion factor of 0.648/pi from micron/meter=microradians */
    psf.lD = lambda/D * 0.648/M_PI;
    /* fractional central obscuration */
    psf.centobs = centobs;
    /* monochromatic */
    psf.is_broadband = 0;
    /* wave front error */
    psf.rmswfe = rmswfe_um / lambda;
    if (is_use_wavefront) psf.rmswfe = get_rmswfe_from_fcn(lambda,&wavefront)/lambda;

    /* Set up calibration parameters: e- per exposure per W/m^2
     * Note 1 joule = 5.03411747e18*lambda(um) photons
     * This accounts for Galactic extinction.
     */
    calib_1exp = get_throughput(lambda,&T1) * M_PI/4. * D * D * (1.-centobs*centobs) * 5.03411747e18*lambda * t_exp
                 * pow(10., -0.4*gal_ebv*Galactic_Alambda__EBV(lambda));

    printf("%6.4lf", lambda);
    for(i=0;i<=10;i++)
      printf(" %11.5lE", get_1sigma_flux_1exp(0.1*i,&psf,var_1exp,calib_1exp)*significance_cut/sqrt((double)N_exp));
    printf("\n");
  }

#else

  /* Write output headers */
  printf("  z  |lambda| EE50 | dV/dz/dA  | Flim@%4.2lf\"| n targets | dN/dz/dA  |geom mean F| siglnF |skew lnF|kurt lnF|Num ph 1exp|\n", R12_REF);
  printf("     |  um  |arcsec| Mpc3/deg2 |   W/m2    |  Mpc^-3   |   deg^-2  |   W/m2    |        |        |        | 1e-19 W/m2|\n");

  /* --- FINISHED WITH INPUTS: NOW PROCEED TO COMPUTATION --- */

  /* Initializations of integrals over redshift */
  dNdA = 0.;

  /* Loop over redshifts within wavelength range */
#ifdef OIII_GAL
  for(z=DZ_OUT; z<lambda_max/LAMBDA0_OIII-1; z+=DZ_OUT)
  if (z>lambda_min/LAMBDA0_OIII-1)
#elif OII_GAL
  for(z=DZ_OUT; z<lambda_max/LAMBDA0_OII-1; z+=DZ_OUT)
  if (z>lambda_min/LAMBDA0_OII-1)
#else
  for(z=DZ_OUT; z<lambda_max/LAMBDA0_HA-1; z+=DZ_OUT)
  if (z>lambda_min/LAMBDA0_HA-1)
#endif
  {
    /* Get the wavelength of H alpha at this redshift */
#ifdef OIII_GAL
    lambda = (1.+z) * LAMBDA0_OIII;
#elif OII_GAL
    lambda = (1.+z) * LAMBDA0_OII;
#else
    lambda = (1.+z) * LAMBDA0_HA;
#endif

    /* Setup PSF parameters: */
    /* Linear dispersion */
    psf.linDisp = linDisp;
    /* Pixel scale, in arcsec */
    psf.pixscale = pixscale;
    /* 1 sigma jitter + diffusion, in arcsec */
    psf.sigma = sqrt((PIX_CD/PIXSIZE_PHYS*pixscale)*(PIX_CD/PIXSIZE_PHYS*pixscale) + jitter*jitter);
    /* lambda/D in arcsec -- hence conversion factor of 0.648/pi from micron/meter=microradians */
    psf.lD = lambda/D * 0.648/M_PI;
    /* fractional central obscuration */
    psf.centobs = centobs;
    /* monochromatic */
    psf.is_broadband = 0;
    /* wave front error */
    psf.rmswfe = rmswfe_um / lambda;
    if (is_use_wavefront) psf.rmswfe = get_rmswfe_from_fcn(lambda,&wavefront)/lambda;

    /* Set up calibration parameters: e- per exposure per W/m^2
     * Note 1 joule = 5.03411747e18*lambda(um) photons
     * This accounts for Galactic extinction.
     */
    calib_1exp = get_throughput(lambda,&T1) * M_PI/4. * D * D * (1.-centobs*centobs) * 5.03411747e18*lambda * t_exp
                 * pow(10., -0.4*gal_ebv*Galactic_Alambda__EBV(lambda));

    nbar = get_n_galaxies(z,&psf,var_1exp,calib_1exp,N_exp,significance_cut,gal_pop_model,stats)
           * completeness;

    printf("%5.3lf %6.4lf %6.4lf %11.5lE %11.5lE %11.5lE %11.5lE %11.5lE %8.4lf %8.4lf %8.4lf %11.5lE\n", z,
      lambda, get_gal_size(0.5, 0.0, &psf),
      computeDistance(z)*computeDistance(z)/computeHubble(z)/SQDEG_PER_SR,
      get_limflux(R12_REF,&psf,var_1exp,calib_1exp,EWD_REF,(double)N_exp,significance_cut),
      nbar,
      nbar*computeDistance(z)*computeDistance(z)/computeHubble(z)/SQDEG_PER_SR,
      stats[0], stats[1], stats[2], stats[3], calib_1exp*1e-19
    );
    /* Sum up integrals over redshift */
#ifdef OIII_GAL
    dz =  (z+0.5*DZ_OUT>lambda_max/LAMBDA0_OIII-1? lambda_max/LAMBDA0_OIII-1: z+0.5*DZ_OUT)
        - (z-0.5*DZ_OUT<lambda_min/LAMBDA0_OIII-1? lambda_min/LAMBDA0_OIII-1: z-0.5*DZ_OUT);
#elif OII_GAL
    dz =  (z+0.5*DZ_OUT>lambda_max/LAMBDA0_OII-1? lambda_max/LAMBDA0_OII-1: z+0.5*DZ_OUT)
        - (z-0.5*DZ_OUT<lambda_min/LAMBDA0_OII-1? lambda_min/LAMBDA0_OII-1: z-0.5*DZ_OUT);
#else
    dz =  (z+0.5*DZ_OUT>lambda_max/LAMBDA0_HA-1? lambda_max/LAMBDA0_HA-1: z+0.5*DZ_OUT)
        - (z-0.5*DZ_OUT<lambda_min/LAMBDA0_HA-1? lambda_min/LAMBDA0_HA-1: z-0.5*DZ_OUT);
#endif

    /* Galaxies per square degree */
    dNdA += dz * nbar*computeDistance(z)*computeDistance(z)/computeHubble(z)/SQDEG_PER_SR;
  }
#endif

  /* --- Print redshift independent summary statistics --- */

  printf("\nSummary statistics:\n");
  printf("Sky background flux:        %12.5lE e-/pix/s\n", tot_bkgnd * M_PI/4.* D * D * (1.-centobs*centobs) * pixscale * pixscale);
  printf("Thermal background flux:    %12.5lE e-/pix/s\n", thermalbkgnd[0]);
  printf("  telescope:                %12.5lE e-/pix/s\n", thermalbkgnd[1]);
  printf("  upstream:                 %12.5lE e-/pix/s\n", thermalbkgnd[2]);
  printf("  downstream:               %12.5lE e-/pix/s\n", thermalbkgnd[3]);
  printf("Noise variance per unit solid angle in one exposure:\n");
  printf("                sky only:   %12.5lE e-^2/arcsec^2\n", var_1exp_skyonly);
  printf("                   total:   %12.5lE e-^2/arcsec^2\n", var_1exp);
#ifndef PZCAL_MODE
  printf("Available galaxy density:   %12.5lE gal/deg^2\n", dNdA);
#endif
  printf("\n");
#endif

  /* --- SPCONT ETC --- */

#ifdef SPCONT_MODE

  /* Report general properties */
  if (!is_an_ifu) {
    printf("\nGeneral properties:\n");
    printf("Trace length:               %12.5lE pixels\n", spDispersion/pixscale*log(lambda_max/lambda_min));
    printf("Trace length:               %12.5lE arcsec\n", spDispersion*log(lambda_max/lambda_min));
    printf("Sky background flux:        %12.5lE e-/pix/s\n", tot_bkgnd * M_PI/4.* D * D * (1.-centobs*centobs) * pixscale * pixscale);
    printf("Thermal background flux:    %12.5lE e-/pix/s\n", thermalbkgnd[0]);
    printf("  telescope:                %12.5lE e-/pix/s\n", thermalbkgnd[1]);
    printf("  upstream:                 %12.5lE e-/pix/s\n", thermalbkgnd[2]);
    printf("  downstream:               %12.5lE e-/pix/s\n", thermalbkgnd[3]);
    printf("Noise variance per unit solid angle in one exposure:\n");
    printf("                sky only:   %12.5lE e-^2/arcsec^2\n", var_1exp_skyonly);
    printf("                   total:   %12.5lE e-^2/arcsec^2\n", var_1exp);
    printf("\n");
  }

  /* Number of wavelength steps */
  nstep = 16;

  /* Print table header */
  if (!is_an_ifu) {
    printf("|lambda|thetaeff|Resoln|counts per|ABmag at|ABmag at|ABmag at|\n");
    printf("|      |        | l/dl |1D pix per| S/N=1  | S/N=5  |S/N=5per|\n");
    printf("|  um  | arcsec |      |exp per Jy|per pix |per resl|1e4 km/s|\n");
  } else {
    printf("|lambda|eff numb|Pixels|counts per|ABmag at|ABmag at|Zodi bkgnd|Tot bkgnd |\n");
    printf("|      | of 2D  |per ln|1D pix per| S/N=1  |S/N=5per|per pixel |per pixel |\n");
    printf("|  um  | pixels |lambda|exp per Jy|per pix |1e4 km/s| e^2/1exp | e^2/1exp |\n");
  }

  /* Now loop over wavelengths */
  for(lambda=lambda_min; lambda<1.0000001*lambda_max; lambda+=(lambda_max-lambda_min)/(double)nstep) {

    /* Setup PSF parameters: */
    /* Linear dispersion */
    psf.linDisp = 0.;
    /* Pixel scale, in arcsec */
    psf.pixscale = pixscale;
    /* 1 sigma jitter + diffusion, in arcsec. For the IFU case, charge 'diffusion' only smears the trace
     * in 1 of the spatial dimensions, since it is a detector effect and does not move light onto a
     * different slice. (Of course it also smears the spectral direction but we are using 'resolution'
     * in this mode to denote the 2-pixel width.)
     */
    if (is_an_ifu) {
      psf.sigma = sqrt((PIX_CD/PIXSIZE_PHYS*pixscale)*(PIX_CD/PIXSIZE_PHYS*pixscale)/2. + jitter*jitter);
    } else {
      psf.sigma = sqrt((PIX_CD/PIXSIZE_PHYS*pixscale)*(PIX_CD/PIXSIZE_PHYS*pixscale) + jitter*jitter);
    }
    /* lambda/D in arcsec -- hence conversion factor of 0.648/pi from micron/meter=microradians */
    psf.lD = lambda/D * 0.648/M_PI;
    /* fractional central obscuration */
    psf.centobs = centobs;
    /* monochromatic */
    psf.is_broadband = 0;
    /* wave front error */
    psf.rmswfe = rmswfe_um / lambda;
    if (is_use_wavefront) psf.rmswfe = get_rmswfe_from_fcn(lambda,&wavefront)/lambda;

    /* Set up calibration parameters: e- per exposure per W/m^2   
     * Note 1 joule = 5.03411747e18*lambda(um) photons
     * This accounts for Galactic extinction.
     */
    calib_1exp = get_throughput(lambda,&T1) * M_PI/4. * D * D * (1.-centobs*centobs) * 5.03411747e18*lambda * t_exp
                 * pow(10., -0.4*gal_ebv*Galactic_Alambda__EBV(lambda));
    /* Then convert calibration to e- per exposure per 1D-pixel per Jy */
    calib_1exp *= 2.99792458e-12 * pixscale/spDispersion / lambda;

    /* Separately calculate sensitivity for the case of [i] an IFU or [ii] slitless spectroscopy */
    if (is_an_ifu) {
      /* IFU SPECTROSCOPY */
      /* In this case, backgrounds are re-computed at each wavelength since they vary! */
      /* Effective number of 2D spatial pixels used */
      npix_2d = get_gal_Omeff(r_eff, &psf) / (pixscale*pixscale);
      /* Get the background flux per focal plane pixel */
      dll_1pix = pixscale/spDispersion;
      l1 = lambda * exp(-0.5*dll_1pix);
      l2 = lambda * exp( 0.5*dll_1pix);

      /* Noise level in photons/m^2/arcsec^2/sec */
      zodi_bkgnd = get_zodi_bkgnd(ecl_lat, ecl_dlon, l1, l2, &TAll);
      /* ... converted to e-^2/arcsec^2 */
      var_1exp_skyonly = zodi_bkgnd * M_PI/4.* D * D * (1.-centobs*centobs) * t_exp;
      /* Thermal background, e/p/s */
      for(i=0; i<4; i++)
        thermalbkgnd[i] = getThermalBkgndTot(&TAll, l1, l2, filterthroughput,
                          fratio, centobs, &thermal, i);
      /* Noise in 1 exposure per pixel */
      var_1exp_bkgnd = var_1exp_skyonly
                       + (darkcurrent+thermalbkgnd[0])*t_exp/(pixscale*pixscale);
      varperpixel_1exp = getNoiseTotal(t_exp,var_1exp_bkgnd*pixscale*pixscale,readnoise,sample_mode);
      /* The "final" varperpixel_1exp is e^2 per pixel in a single exposure */

      /* Now the flux for S/N=1 per spectral pixel */
      flux1 = sqrt(varperpixel_1exp*npix_2d)/calib_1exp/sqrt((double)N_exp); /* in Jy */
      mag1 = 2.5*log(3631./flux1)/log(10.); /* in mag AB */
      printf(" %6.4lf %8.4lf %6.1f %10.4lE %8.4lf %8.4lf %10.4lE %10.4lE\n",
        lambda, npix_2d, spDispersion/pixscale, calib_1exp, mag1,
        mag1 - 2.5*log(5.)/log(10.) + 1.25*log(spDispersion/pixscale/29.972458)/log(10.),
        var_1exp_skyonly*pixscale*pixscale, varperpixel_1exp
      );

    } else {
      /* SLITLESS SPECTROSCOPY */
      /* Trace width */
      thetaeff = get_gal_thetaeff(r_eff, &psf);
      /* Flux for S/N=1 per pixel */
      flux1 = sqrt(var_1exp*pixscale*thetaeff)/calib_1exp/sqrt((double)N_exp); /* in Jy */
      mag1 = 2.5*log(3631./flux1)/log(10.); /* in mag AB */
      printf(" %6.4lf %8.4lf %6.1f %10.4lE %8.4lf %8.4lf %8.4lf\n",
        lambda, thetaeff, spDispersion/thetaeff, calib_1exp, mag1,
        mag1 - 2.5*log(5.)/log(10.) + 1.25*log(thetaeff/pixscale)/log(10.),
        mag1 - 2.5*log(5.)/log(10.) + 1.25*log(spDispersion/pixscale/29.972458)/log(10.));
    }

  }

#endif

  /* --- WL ETC --- */

#ifdef WL_MODE

  /* Generate the PSF */
  /* Linear dispersion */
  psf.linDisp = 0.;
  /* Pixel scale, in arcsec */
  psf.pixscale = pixscale;
  /* 1 sigma jitter + diffusion, in arcsec */
  psf.sigma = sqrt((PIX_CD/PIXSIZE_PHYS*pixscale)*(PIX_CD/PIXSIZE_PHYS*pixscale) + jitter*jitter);
  /* lambda/D in arcsec -- hence conversion factor of 0.648/pi from micron/meter=microradians */
  psf.lD = sqrt(lambda_min*lambda_max)/D * 0.648/M_PI;
  /* fractional central obscuration */
  psf.centobs = centobs;
  /* filter width */
  psf.is_broadband = 1;
  psf.dll = log(lambda_max/lambda_min);
  /* wave front error */
  psf.rmswfe = rmswfe_um / sqrt(lambda_min*lambda_max);

  /* Theoretical cycles/pixel */
  cpp = D/lambda_min*pixscale/(0.648/M_PI);

  /* PSF size & minimum usable size */
  ee50 = get_gal_size(0.5,0.0,&psf);
  min_reff = ee50 * sqrt(min_res_factor/(1-min_res_factor));

  /* Set up calibration parameters: e- per exposure for a 0 mag AB source
   * Note 0 mag AB = 3.631e-23 W/m2/Hz = 5.479865e10 photons/m2 per ln(lambda)
   * This accounts for Galactic extinction.
   */
  throughput = 0;
  for(i=0;i<1000;i++) {
    lambda = lambda_min * pow(lambda_max/lambda_min, (i+0.5)/1000.0);
    throughput += get_throughput(lambda,&T1);
  }
  throughput /= 1000;
  calib_1exp = M_PI/4. * D * D * (1.-centobs*centobs) * 5.479865e10 * t_exp
               * pow(10., -0.4*gal_ebv*Galactic_Alambda__EBV(sqrt(lambda_min*lambda_max)))
               * throughput * log(lambda_max/lambda_min);

#ifdef OUT_MTF_TABLE
  printf("\n=== Table of modulation transfer functions ===\n");
  for(i=0;i<=240;i++) {
    printf("%6.3lf %7.5lf\n", i*.05, get_MTF(i*.05,0,&psf));
  }
#endif

  printf("\nGeneral properties:\n");
  printf("PSF EE50:                   %12.5lE arcsec\n", ee50);
  printf("PSF effective area:         %12.5lE arcsec^2\n", get_gal_Omeff(0,&psf));
  printf("PSF encircled energy:       %8.6lf, %8.6lf, %8.6lf (at 0.1\", 0.25\", 0.5\")\n", get_gal_fracE(0.1,0,&psf), get_gal_fracE(0.25,0,&psf), get_gal_fracE(0.5,0,&psf));
#ifdef OUT_EXTRA_PSF_PROPERTIES
#if 0
  printf("EE at 50.0/52.7/58.3: %8.6lf, %8.6lf, %8.6lf as\n", get_gal_size(0.500,0.0,&psf), get_gal_size(0.527,0.0,&psf), get_gal_size(0.583,0.0,&psf));
  printf("EE at 70.0/73.8/81.6: %8.6lf, %8.6lf, %8.6lf as\n", get_gal_size(0.700,0.0,&psf), get_gal_size(0.738,0.0,&psf), get_gal_size(0.816,0.0,&psf));
  printf("EE at 80.0/84.3/93.3: %8.6lf, %8.6lf, %8.6lf as\n", get_gal_size(0.800,0.0,&psf), get_gal_size(0.843,0.0,&psf), get_gal_size(0.933,0.0,&psf));
#endif
  printf("PSF centroid error coeff:   %12.5lE arcsec\n", 1./sqrt(get_gal_Omeff(0,&psf)*get_gal_centroid_integral(0.,&psf)));
  printf("PSF peak surface brightness:%12.5lE /pix\n", get_gal_PeakInt(0,&psf)*pixscale*pixscale);
#endif
  printf("Spatial frequency cutoff:   %12.5lE cycles/pix\n", cpp);
  if (cpp<0.5) {
    printf("Sampling case:              Oversampled\n");
  } else if (cpp<1) {
    printf("Sampling case:              Weakly undersampled\n");
  } else {
    printf("Sampling case:              Strongly undersampled\n");
  }
  printf("Min usable galaxy r_eff:    %12.5lE arcsec\n", min_reff);
  printf("Sky background flux:        %12.5lE e-/pix/s\n", tot_bkgnd * M_PI/4.* D * D * (1.-centobs*centobs) * pixscale * pixscale);
  printf("Thermal background flux:    %12.5lE e-/pix/s\n", thermalbkgnd[0]);
  printf("  telescope:                %12.5lE e-/pix/s\n", thermalbkgnd[1]);
  printf("  upstream:                 %12.5lE e-/pix/s\n", thermalbkgnd[2]);
  printf("  downstream:               %12.5lE e-/pix/s\n", thermalbkgnd[3]);
  printf("Noise variance per unit solid angle in one exposure:\n");
  printf("                sky only:   %12.5lE e-^2/arcsec^2\n", var_1exp_skyonly);
  printf("                   total:   %12.5lE e-^2/arcsec^2\n", var_1exp);
  printf("Source counts per exposure:\n");
  printf("            at AB mag 20:   %12.5lE e-\n", calib_1exp*pow(10, -8.0));
  printf("            at AB mag 21:   %12.5lE e-\n", calib_1exp*pow(10, -8.4));
  printf("            at AB mag 22:   %12.5lE e-\n", calib_1exp*pow(10, -8.8));
  printf("            at AB mag 23:   %12.5lE e-\n", calib_1exp*pow(10, -9.2));
  printf("            at AB mag 24:   %12.5lE e-\n", calib_1exp*pow(10, -9.6));
  printf("            at AB mag 25:   %12.5lE e-\n", calib_1exp*pow(10,-10.0));
  printf("5 sigma pt src threshold:   %7.3lf mag AB\n", 2.5*log(calib_1exp*sqrt(N_exp)/sqrt(var_1exp*get_gal_Omeff(0.,&psf))/5.0)/log(10));
  printf("5 sigma ext src threshold:  %7.3lf mag AB [r_eff=0.3\"]\n", 2.5*log(calib_1exp*sqrt(N_exp)/sqrt(var_1exp*get_gal_Omeff(0.3,&psf))/5.0)/log(10));
  printf("\n");

#ifdef WFE_OVERRIDE
  return(0);
#endif

  /* Return shape limiting magnitude table */
  printf("| r_eff |Om_eff |penalty|resolution|lim mag |S/N at |\n");
  printf("|       |       |factor |  factor  |(shapes)|shape  |\n");
  printf("|arcsec |arcsec2|       |          |   AB   |lim mag|\n");
  for(r_eff=.01;r_eff<=1.001;r_eff*=pow(10,.05)) if (r_eff>min_reff) {
    printf(" %7.5lf %7.4lf %7.4lf  %7.5lf   %8.5lf %7.3lf\n", r_eff, get_gal_Omeff(r_eff,&psf),
      get_shape_penalty_factor(r_eff,&psf),
      r_eff*r_eff/(r_eff*r_eff+ee50*ee50),
      get_maglim_shape(r_eff,&psf,N_exp,var_1exp,calib_1exp,max_ellip_err),
      calib_1exp*sqrt(N_exp)*pow(10,-0.4*get_maglim_shape(r_eff,&psf,N_exp,var_1exp,calib_1exp,max_ellip_err))/sqrt(var_1exp*get_gal_Omeff(r_eff,&psf))
    );
  }

  /* Output galaxy yields */
  printf("\n");
  printf("|  z  |  dN/dz/dA |dNeff/dz/dA|\n");
  printf("|     |   deg^-2  |   deg^-2  |\n");
  ntot = nefftot = 0;
  for(iz=0;iz<20;iz++) {
    dn = get_dNdA_WL_gal2(min_reff, 100*min_reff, 0.199999*iz, 0.199999*(iz+1), GalaxyCat, lambda_min, lambda_max, &psf, N_exp, var_1exp, calib_1exp, max_ellip_err, &dneff,
#ifdef OUT_WL_CAT
         OutGalCat
#else
         NULL
#endif
         , iz>0?'a':'w');

    printf(" %5.3lf %11.5lE %11.5lE\n", 0.2*iz+0.1, dn/0.2, dneff/0.2);
    ntot += dn;
    nefftot += dneff;
  }
  printf("\n");
  printf("Weak lensing n:             %12.5lE gal/deg^2\n", ntot);
  printf("                          = %12.5lE gal/arcmin^2\n", ntot/3600.);
  printf("Weak lensing n_eff:         %12.5lE gal/deg^2\n", nefftot);
  printf("                          = %12.5lE gal/arcmin^2\n", nefftot/3600.);
  
#endif

  wait_for_user();
  return(0);
}
