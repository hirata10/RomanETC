#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "nr_utils.h"
#include "utils.h"
#include "mnemonics.h"

/* Gets Albrecht-Bernstein figure of merit for a covariance matrix of
 * w in redshift bins
 */
double get_fom_ab(double **cov, int nw) {
  int i,j,nrot;
  double **new_cov;
  double **eigenvecs, *eigenvals;
  double fom;

  new_cov=dmatrix(0, nw-1, 0, nw-1);
  eigenvecs=dmatrix(0, nw-1, 0, nw-1);
  eigenvals=dvector(0, nw-1);

  /* Copy covariance matrix since NR jacobi destroys input matrix */
  for(i=0; i<nw; i++)
    for(j=0; j<nw; j++)
      new_cov[i][j] = cov[i][j];

  jacobi(new_cov, nw, eigenvals, eigenvecs, &nrot);

  for(i=0; i<nw; i++) {
    printf("sigma=%12.5le (", sqrt(eigenvals[i]));
    for(j=0; j<nw; j++) printf(" %6.3lf", eigenvecs[j][i]);
    printf(")\n");
  }

  fom = 1.;
  for(i=0; i<nw; i++)
    if (eigenvals[i]<1.)
      fom /= sqrt(eigenvals[i]);

  free_dmatrix(new_cov, 0, nw-1, 0, nw-1);  
  free_dmatrix(eigenvecs, 0, nw-1, 0, nw-1);  
  free_dvector(eigenvals, 0, nw-1);
  return(fom);
}

/* DETF figure-of-merit */
double get_fom_detf(double **cov, int nw) {
  int i,j;
  double **new_cov;
  double *vec1, *vec2;
  double a, area;
  double **F_wwa;

  F_wwa = dmatrix(0,1,0,1);
  new_cov=dmatrix(0, nw-1, 0, nw-1);
  vec1 = dvector(0, nw-1);
  vec2 = dvector(0, nw-1);

  /* Construct vectors d(p^i)/dw, d(p^1)/dwa */
  for(i=0; i<nw; i++) {
    a = 1. - DELTA_A*(i+0.5);
    vec1[i] = 1.;
    vec2[i] = 1.-a;
  }

  /* Copy covariance matrix since inversion destroys it */
  for(i=0; i<nw; i++)
    for(j=0; j<nw; j++)
      new_cov[i][j] = cov[i][j];
  gaussjinv(new_cov, nw);

  /* Get inner products with new_cov */
  F_wwa[0][0] = F_wwa[0][1] = F_wwa[1][1] = 0.;
  for(i=0; i<nw; i++)
    for(j=0; j<nw; j++) {
      F_wwa[0][0] += new_cov[i][j] * vec1[i] * vec1[j];
      F_wwa[0][1] += new_cov[i][j] * vec1[i] * vec2[j];
      F_wwa[1][1] += new_cov[i][j] * vec2[i] * vec2[j];
    }
  F_wwa[1][0]=F_wwa[0][1];

  /* Get covariance of w,w_a */
  gaussjinv(F_wwa,2);
  area = sqrt(F_wwa[0][0]*F_wwa[1][1]-F_wwa[0][1]*F_wwa[0][1]);
  printf("sig(w)=%9.6lf; sig(wa)=%9.6lf; corr=%9.6lf\n", sqrt(F_wwa[0][0]), sqrt(F_wwa[1][1]), F_wwa[0][1]/sqrt(F_wwa[0][0]*F_wwa[1][1]));
  printf("a_pivot=%9.6lf, sig(wp)=%9.6lf\n", 1.+F_wwa[0][1]/F_wwa[1][1],
    sqrt(F_wwa[0][0]-F_wwa[0][1]*F_wwa[0][1]/F_wwa[1][1]));

  free_dmatrix(new_cov, 0, nw-1, 0, nw-1);
  free_dvector(vec1, 0, nw-1);
  free_dvector(vec2, 0, nw-1);
  free_dmatrix(F_wwa, 0,1,0,1);
  return(1./(area));
}

/* Construct Jacobian for transformation to SWG parameters.
 * Jac[i][j] = d(our parameter #i)/d(SWG parameter #j)
 *
 * SWG parameter order:
 * n_s, omega_m, omega_b, omega_k, omega_de
 * gamma, Mag_SN, G0, ln A_s
 * w_0 .. w_35
 */
void get_jacobian_swg(double **Jac, COSMOPARAM *p) {
  int i,j;

  /* Clear */
  for(i=0;i<NPARAMTOT;i++) for(j=0;j<NSWG;j++) Jac[i][j] = 0.;

  /* Matter composition */
  Jac[NWPARAM+1][2] = 1.; /* baryons */
  Jac[NWPARAM+2][1] = 1.; /* matter */
  Jac[NWPARAM+7][3] = 1./p->h/p->h; /* curvature */

  /* We use Hubble constant instead of omega_de.
   * h = sqrt(omega_m + omega_k + omega_de)
   */
  Jac[NWPARAM+6][1] = 1./(2.*p->h);
  Jac[NWPARAM+6][3] = 1./(2.*p->h);
  Jac[NWPARAM+6][4] = 1./(2.*p->h);

  /* Primordial spectrum */
  Jac[NWPARAM+4][0] = 1.; /* ns */
  Jac[NWPARAM+3][8] = 0.5; /* ln delta^2 -> ln delta */

  /* SN magnitude */
  Jac[NWPARAM+9][6] = -0.4*log(10.);

  /* Growth */
  Jac[NWPARAM+10][7] = 1.;
  Jac[NWPARAM+11][5] = 1.;

  /* w's */
  for(i=0; i<NWPARAM; i++) Jac[i][9+i] = 1.;
}

/* Write SWG space Fisher matrix. */
void write_swg_fisher(double **F, COSMOPARAM *p, char FileName[]) {
  FILE *fp;
  double **Jac, **myF, *tc;
  double **F_swg;
  int i,j,i_s,j_s;

  Jac = dmatrix(0,NPARAMTOT-1,0,NSWG-1);
  myF = dmatrix(0,NPARAMTOT-1,0,NPARAMTOT-1);
  F_swg = dmatrix(0,NSWG-1,0,NSWG-1);
  tc = dvector(0,NPARAMTOT-1);

  /* Get Fisher matrix after marginalization over tau. */
  for(i=0; i<NPARAMTOT; i++)
    tc[i] = F[NWPARAM+8][NWPARAM+8]>1e-16? F[i][NWPARAM+8]/sqrt(F[NWPARAM+8][NWPARAM+8]): 0;
  for(i=0; i<NPARAMTOT; i++)
    for(j=0; j<NPARAMTOT; j++)
      myF[i][j] = F[i][j] - tc[i]*tc[j];

  /* Projection */
  get_jacobian_swg(Jac,p);
  for(i_s=0; i_s<NSWG; i_s++)
    for(j_s=0; j_s<NSWG; j_s++) {
      F_swg[i_s][j_s] = 0.;
      for(i=0; i<NPARAMTOT; i++)
        for(j=0; j<NPARAMTOT; j++)
          F_swg[i_s][j_s] += myF[i][j] * Jac[i][i_s] * Jac[j][j_s];
    }
  fp = fopen(FileName, "w");
  for(i_s=0; i_s<NSWG; i_s++)
    for(j_s=i_s; j_s<NSWG; j_s++) {
#if 0
      if (F_swg[i_s][i_s]>0 && F_swg[j_s][j_s]>0) {
        fprintf(fp, "%19.12le # %13.10lf %s %s\n", F_swg[i_s][j_s],
         F_swg[i_s][j_s]/sqrt(F_swg[i_s][i_s]*F_swg[j_s][j_s]),
         swgParam[i_s], swgParam[j_s]);
      } else {
        fprintf(fp, "%19.12le #      nan      %s %s\n", F_swg[i_s][j_s], swgParam[i_s], swgParam[j_s]);
      }
#endif
      if (F_swg[i_s][j_s]!=0) fprintf(fp, "%2d %2d %19.12le\n", i_s, j_s, F_swg[i_s][j_s]);
    }
  fclose(fp);

#if 0
  printf("%19.12le %19.12le %19.12le\n", F_swg[5][5], F_swg[5][7], F_swg[7][7]);
  exit(1);
#endif

  free_dmatrix(Jac,0,NPARAMTOT-1,0,NSWG-1);
  free_dmatrix(myF,0,NPARAMTOT-1,0,NPARAMTOT-1);
  free_dvector(tc,0,NPARAMTOT-1);
  free_dmatrix(F_swg,0,NSWG-1,0,NSWG-1);
}
