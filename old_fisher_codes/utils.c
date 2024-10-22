#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define DEF_PARAMETER_MAPPINGS_HERE
#include "nr_utils.h"
#include "utils.h"

/* gaussjinv
 * *** INVERTS A MATRIX ***
 *
 * Takes matrix A[0..n-1][0..n-1] and inverts it.  The original matrix A is destroyed.
 *
 * Arguments:
 * > A: matrix to be inverted
 *   n: dimension of matrix
 */

void gaussjinv(double **A, int n) {
   /* This routine inverts a matrix A and stores the result back to
    * A.  The matrices are taken to be n by n.  This routine is a
    * slight modification to Numerical Recipes gaussj, see Section
    * 2.1 p.39 in NR C 2ed.  Note we've replaced float with double.
    */

   int *indxc, *indxr, *ipiv;
   int i,icol,irow,j,k,l,ll;
   double big,dum,pivinv,temp;

   /* The integer arrays ipiv, indxr, and indxc are used for    */
   /* bookkeeping on the pivoting.                              */
   indxc = ivector(0,n-1);
   indxr = ivector(0,n-1);
   ipiv = ivector(0,n-1);

   for(j=0;j<n;j++) ipiv[j]=0;
   for(i=0;i<n;i++) { /* This is the main loop over the columns */
      big=0.0;                              /* to be reduced. */
      for(j=0;j<n;j++)       /* This is the outer loop of the */
         if (ipiv[j] != 1)   /* search for a pivot element. */
            for (k=0;k<n;k++) {
               if (ipiv[k] == 0) {
                  if (fabs(A[j][k]) >= big) {
                     big=fabs(A[j][k]);
                     irow=j;
                     icol=k;
                  }
               } else if (ipiv[k] > 1)
                  nrerror("gaussj: Singular Matrix-1");
               } /* end for(k) loop */
      /* also end for(j) loop */
      ++(ipiv[icol]);

      /* We now have the pivot element, so we interchange
       * rows, if needed, to put the pivot element on the
       * diagonal.  The columns are not physically
       * interchanged, only relabeled: indxc[i], the column
       * of the ith pivot element, is the ith column that is
       * reduced, while indxr[i] is the row in which that
       * pivot element was originally located.  If indxr[i]
       * != indxc[i] there is an implied column interchange.
       * With this form of bookkeeping, the solution Bs (of
       * which we don't have any!) will end up in the correct
       * order, and the inverse matrix (i.e. that remains in
       * A) will be scrambled by columns.
       */

      if (irow != icol) {
         for (l=0;l<n;l++) {
            /* Swap A[irow][l] and A[icol][l] */
            temp = A[irow][l];
            A[irow][l] = A[icol][l];
            A[icol][l] = temp;
         }
      } /* end if */
      indxr[i]=irow; /* We are now ready to divide the pivot  */
      indxc[i]=icol; /* row by the pivot element, located at  */
                                            /* irow and icol. */
      if (A[icol][icol] == 0.0)
         nrerror("gaussj: Singular Matrix-2");
      pivinv=1.0/A[icol][icol];
      A[icol][icol]=1.0;
      for (l=0;l<n;l++) A[icol][l] *= pivinv;
      for (ll=0;ll<n;ll++) /* Now we reduce the rows, except  */
         if (ll != icol) { /* for the pivot one, of course. */
            dum=A[ll][icol];
            A[ll][icol]=0.0;
            for(l=0;l<n;l++) A[ll][l] -= A[icol][l]*dum;
         }
   }

   /* This is the end of the main loop over columns of the
    * reduction. It only remains to unscramble the solution in
    * view of the column interchanges.  We do this by
    * interchanging pairs of columns in the reverse order that
    * the permutation was built up.
    */
   for (l=n-1;l>=0;l--) {
      if (indxr[l] != indxc[l])
         for (k=0;k<n;k++) {
            temp = A[k][indxr[l]];
            A[k][indxr[l]] = A[k][indxc[l]];
            A[k][indxc[l]] = temp;
         }
   } /* And we are done. */

   /* Clean up memory */
   free_ivector(ipiv,0,n-1);
   free_ivector(indxr,0,n-1);
   free_ivector(indxc,0,n-1);
}

/* Resets cosmology to default values */
void set_default_cosmology(COSMOPARAM *p) {
  int i;
  p->SNflux = 1.;
  p->tau = 0.07;
  p->T_cmb = 2.728;
  p->obh2 = 0.022;
  p->omh2 = 0.13;
  p->lndelta_zeta = -10.0;
  p->ns = 0.94;
  p->alphas = 0.;
  p->h = 0.71;
  p->Omega_K = 0.;
  for(i=0; i<NWPARAM; i++) p->w[i] = 0.;
  p->growthAmp = 1.;
  p->growthSlope = 0.;
  p->growthCurve = 0.;
  p->onuh2 = 0.;
#if 1
  /* DETF parameters */
  p->SNflux = 1.;
  p->tau = 0.17;
  p->T_cmb = 2.728;
  p->obh2 = 0.024;
  p->omh2 = 0.146;
  p->lndelta_zeta = -9.89;
  p->ns = 1.0;
  p->alphas = 0.;
  p->h = 0.725;
  p->Omega_K = 0.;
  for(i=0; i<NWPARAM; i++) p->w[i] = 0.;
#endif
#if 0
  /* Song/Knox parameters */
  p->SNflux = 1.;
  p->tau = 0.07;
  p->T_cmb = 2.728;
  p->obh2 = 0.021;
  p->omh2 = 0.146;
  p->lndelta_zeta = -10.127;
  p->ns = 1.0;
  p->alphas = 0.;
  p->h = 0.655;
  p->Omega_K = 0.;
  for(i=0; i<NWPARAM; i++) p->w[i] = 0.;
#endif
#if 0
  /* Pullen parameters */
  p->SNflux = 1.;
  p->tau = 0.07;
  p->T_cmb = 2.728;
  p->obh2 = 0.0245;
  p->omh2 = 0.147;
  p->lndelta_zeta = -9.927191525;
  p->ns = 1.0;
  p->alphas = 0.;
  p->h = 0.7;
  p->Omega_K = 0.;
  for(i=0; i<NWPARAM; i++) p->w[i] = 0.;
#endif
#if 1
  /* WMAP5 parameters */
  p->SNflux = 1.;
  p->tau = 0.087;
  p->T_cmb = 2.725;
  p->obh2 = 0.02273;
  p->omh2 = 0.1326;
  p->lndelta_zeta = -9.9814;
  p->ns = 0.963;
  p->alphas = 0.;
  p->h = 0.719;
  p->Omega_K = 0.;
  for(i=0; i<NWPARAM; i++) p->w[i] = 0.;
#endif
#if 0
  /* WMAP5 + SN + BAO parameters */
  p->SNflux = 1.;
  p->tau = 0.084;
  p->T_cmb = 2.725;
  p->obh2 = 0.02265;
  p->omh2 = 0.1369;
  p->lndelta_zeta = -9.9122;
  p->ns = 0.960;
  p->alphas = 0.;
  p->h = 0.701;
  p->Omega_K = 0.;
  for(i=0; i<NWPARAM; i++) p->w[i] = 0.;
#endif
//p->obh2 = 0.021125; p->omh2 = 0.147875; p->h = 0.65; p->lndelta_zeta = -9.958475212282153; p->ns = 1;
#if 0
  /* 737 parameters */
  p->SNflux = 1.;
  p->tau = 0.17;
  p->T_cmb = 2.728;
  p->obh2 = 0.05*0.7*0.7;
  p->omh2 = 0.3*0.7*0.7;
  p->lndelta_zeta = -9.89;
  p->ns = 0.95;
  p->alphas = 0.;
  p->h = 0.7;
  p->Omega_K = 0.;
  for(i=0; i<NWPARAM; i++) p->w[i] = 0.;
#endif
#ifdef PLANCK_PAR
  /* Planck 2015 + everything parameters */
  p->SNflux = 1.;
  p->tau = 0.066;
  p->T_cmb = 2.7255;
  p->obh2 = 0.02230;
  p->omh2 = 0.14170;
  p->lndelta_zeta = -9.9809;
  p->ns = 0.9667;
  p->alphas = 0.;
  p->h = 0.6774;
  p->Omega_K = 0.;
  for(i=0; i<NWPARAM; i++) p->w[i] = 0.;
#endif
#ifdef PEEPLES_PAR
  /* Planck 2015 + everything parameters */
  p->SNflux = 1.;
  p->tau = 0.066;
  p->T_cmb = 2.7255;
  p->obh2 = 0.044*.7*.7;
  p->omh2 = 0.25*.7*.7;
  p->lndelta_zeta = -9.91876960814208;
  p->ns = 0.95;
  p->alphas = 0.;
  p->h = 0.7;
  p->Omega_K = 0.;
  for(i=0; i<NWPARAM; i++) p->w[i] = 0.;
#endif
}

/* Hubble rate, units of c/Mpc */
double get_H(double z, COSMOPARAM *p) {
  int i, imax;
  double a;
  double Omega_L, Omega_M, E, w1, H;
#ifdef DAN_SPLINE
  double Omega_X, i_fl;
#endif

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

#ifdef DAN_SPLINE
  i_fl = log((1.+z)/1.05)/0.14;
  i = (long)floor(i_fl);
  i_fl -= i;
  Omega_X = i<0 || i>=NWPARAM-1? 0: (1.-i_fl)*p->w[i] + i_fl*p->w[i+1];
  H = p->h * sqrt(Omega_L + (Omega_M/a + p->Omega_K)/a/a) / HL_HIMPC / sqrt(1.-Omega_X);
#endif

  return(H);
}

/* Distance, units of Mpc */
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

/* Distance, radial, units of Mpc */
double get_DRC(double z, COSMOPARAM *p) {
  double a, a_current, z_current;
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

  return(chi);
}

/* Logarithmic derivatives of Hubble rate, comoving angular
 * diameter distance, and sound horizon s with respect to
 * i_th parameter, at redshift z.
 */
#include "eisenstein.c"
void get_derivatives(double z, COSMOPARAM *p, int i, double *dlnHdpi, double *dlnDACdpi, double *dlnsdpi) {
  COSMOPARAM pnew;
  double *pinew;
  int j;
  double xp, xm, step=0.002;

  /* Copy cosmological parameters */
  for(j=0; j<NPARAMTOT; j++) {
    pinew = get_ptr_param(&pnew, j);
    *pinew = *get_ptr_param(p,j);
  }
  pinew = get_ptr_param(&pnew, i);

  /* Hubble rate, if wanted */
  if (dlnHdpi!=NULL) {
    *pinew += step;
    xp = get_H(z,&pnew);
    *pinew -= 2*step;
    xm = get_H(z,&pnew);
    *pinew += step;
    *dlnHdpi = (xp-xm)/(xp+xm)/step;
  }

  /* Distance, if wanted */
  if (dlnDACdpi!=NULL) {
    *pinew += step;
    xp = get_DAC(z,&pnew);
    *pinew -= 2*step;
    xm = get_DAC(z,&pnew);
    *pinew += step;
    *dlnDACdpi = (xp-xm)/(xp+xm)/step;
  }

  /* Sound horizon, if wanted */
  if (dlnsdpi!=NULL) {
    *pinew += step;
    xp = TFsound_horizon(pnew.omh2, pnew.obh2/pnew.omh2, pnew.T_cmb);
    *pinew -= 2*step;
    xm = TFsound_horizon(pnew.omh2, pnew.obh2/pnew.omh2, pnew.T_cmb);
    *pinew += step;
    *dlnsdpi = (xp-xm)/(xp+xm)/step;
  }
}

/* Inversion of symmetric matrix, C-->Cinv, of dimension N.
 * Uses Cholesky algorithm.
 *
 * Returns 0 if successful, 1 if failed (non positive definite).
 */
int sym_matrix_invert(double **C, double **Cinv, int N) {
#define S(a,b) (a+b*(long)N)
  int i,j,k;
  double arg;
  double *L, *M;

  L = (double*)malloc((size_t)(N*(long)N*sizeof(double)));
  M = (double*)malloc((size_t)(N*(long)N*sizeof(double)));

  /* Construct the L-matrix: C = L L^T */
  for(j=N-1; j>=0; j--) {
    /* Diagonal element */
    arg = C[j][j];
    for(k=j+1; k<N; k++)
      arg -= L[S(j,k)]*L[S(j,k)];
    if (arg<=0)
      return(1);
    L[S(j,j)] = sqrt(arg);

    /* Off-diagonal elements */
    for(i=j-1; i>=0; i--) {
      arg = C[i][j];
      for(k=j+1; k<N; k++)
        arg -= L[S(i,k)]*L[S(j,k)];
      L[S(i,j)] = arg/L[S(j,j)];
    }
  }

  /* Now the M-matrix */
  for(i=0; i<N; i++) {
    /* Diagonal element */
    M[S(i,i)] = 1./L[S(i,i)];

    /* Off-diagonal elements */
    for(j=i+1; j<N; j++) {
      arg = 0.;
      for(k=i; k<j; k++)
        arg += M[S(i,k)]*L[S(k,j)];
      M[S(i,j)] = -arg/L[S(j,j)];
    }
  }

  /* Now the C-invese */
  for(i=0; i<N; i++)
    for(j=0; j<N; j++) {
      arg = 0.;
      for(k=0; k<=i && k<=j; k++)
        arg += M[S(k,i)]*M[S(k,j)];
      Cinv[i][j] = arg;
    }

  free((char*)L);
  free((char*)M);
#undef S
  return(0);
}

