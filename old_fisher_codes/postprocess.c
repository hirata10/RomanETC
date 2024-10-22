#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define hFIDUCIAL 0.719

const double sigmaS3[] = {
0.02117745582484,
0.04406683490357,
0.06885951269987,
0.09425648689388,
0.12232051602279,
0.16042890163756,
0.19156225771009,
0.24580719841290,
0.28685730339407,
0.34222661959731,
0.41812669087086,
0.42413856796260,
0.50576965954295,
0.61019752392088,
0.63229944816836,
0.74089113331915,
0.89800071314793,
1.05482206655566,
1.18234759518891,
1.30913761640752,
1.45674099381325,
1.69671608877407,
2.12207345319778,
2.51235528574334,
3.20457641071040,
3.69895654073885,
4.24037957867629,
5.03135976815914,
10.09112051129622,
36.36702910468918,
59.09210691328457,
72.05824848624925,
138.84770774411041,
195.27937539529722,
440.36940468675442,
2370.82604641045873};


void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
     printf("Numerical Recipes run-time error...\n");
     printf("%s\n",error_text);
     printf("...now exiting to system...\n");
     exit(1);
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
/* replaced macros, as with dmatrix etc.                   */
{
   double *v;
   long i;

   v=(double *)malloc((size_t) ((nh-nl+2)*sizeof(double)));
   if (!v) nrerror("allocation failure in dvector()");

   /* Sets the newly created vector to zero */
   for(i=0;i<nh-nl+2;i++) v[i] = 0.;

   return(v-nl+1);
}
/* End dvector */

/* free_dvector
 * *** DE-ALLOCATES DOUBLE PRECISION VECTORS ***
 *
 * the vector has range m[nl..nh]
 *
 * This code is in the public domain.
 */

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
   free((char*) (v+nl-1));
}
/* End free_dvector */

/* dmatrix
 * *** ALLOCATES DOUBLE PRECISION MATRICES ***
 *
 * the matrix has range m[nrl..nrh][ncl..nch]
 *
 * This code is in the public domain.
 */

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range               */
/* m[nrl..nrh][ncl..nch]                                       */
/* NR_END has been replaced with its value, 1.                 */
{
   long i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1;
   double **m;

   /* allocate pointers to rows */
   m=(double **) malloc((size_t)((nrow+1)*sizeof(double*)));
   if (!m) nrerror("allocation failure 1 in matrix()");
   m += 1;
   m -= nrl;

   /* allocate rows and set pointers to them */
   m[nrl]=(double *)malloc((size_t)((nrow*ncol+1)*sizeof(double)));
   if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
   m[nrl] += 1;
   m[nrl] -= ncl;

   for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

   /* Sets the newly created matrix to zero */
   for(i=nrl;i<=nrh;i++) for(j=ncl;j<=nch;j++) m[i][j] = 0.;

   /* return pointer to array of pointers to rows */
   return m;
}

/* free_dmatrix
 * *** DE-ALLOCATES DOUBLE PRECISION MATRICES ***
 *
 * the matrix has range m[nrl..nrh][ncl..nch]
 *
 * This code is in the public domain.
 */

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free an double matrix allocated by dmatrix() */
/* replaced NR_END => 1, FREE_ARG => (char *)   */
{
   free((char *) (m[nrl]+ncl-1));
   free((char *) (m+nrl-1));
}

/* ivector
 * *** ALLOCATES LONG INTEGER VECTORS ***
 *
 * the vector has range m[nl..nh]
 *
 * This code is in the public domain.
 */

int *ivector(long nl, long nh)
/* allocate an integer vector with subscript range v[nl..nh] */
{
   int *v;

   v=(int *)malloc((size_t) ((nh-nl+2)*sizeof(int)));
   if (!v) nrerror("allocation failure in ivector()");
   return(v-nl+1);
}
/* End ivector */

/* free_ivector
 * *** DE-ALLOCATES INTEGER VECTORS ***
 *
 * the vector has range m[nl..nh]
 *
 * This code is in the public domain.
 */

void free_ivector(int *v, long nl, long nh)
/* free an integer vector allocated with ivector() */
{
   free((char*) (v+nl-1));
}
/* End free_ivector */

/* Integrator */
void rk4_dbl(double y[], int n, double x, double h, double yout[],
void (*derivs)(double, double [], double []))
{
int i;
double xh,hh,h6,*dym,*dyt,*yt,*dydx;
dydx=dvector(0,n-1);
dym=dvector(0,n-1);
dyt=dvector(0,n-1);
yt=dvector(0,n-1);
hh=h*0.5;
h6=h/6.0;
xh=x+hh;
(*derivs)(x,y,dydx);
for (i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i];
(*derivs)(xh,yt,dyt);
for (i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i];
(*derivs)(xh,yt,dym);
for (i=0;i<n;i++) {
yt[i]=y[i]+h*dym[i];
dym[i] += dyt[i];
}
(*derivs)(x+h,yt,dyt);
for (i=0;i<n;i++)
yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
free_dvector(dydx,0,n-1);
free_dvector(yt,0,n-1);
free_dvector(dyt,0,n-1);
free_dvector(dym,0,n-1);
}

/* Submatrix */
void submatrix(double **in, double **out, int dim) {
  int i,j;
  for(i=0;i<dim;i++)
    for(j=0;j<dim;j++)
      out[i][j] = in[i][j];
}

/* Submatrix, but skips first nskip rows/cols */
void submatrix2(double **in, double **out, int nskip, int dim) {
  int i,j;
  for(i=0;i<dim;i++)
    for(j=0;j<dim;j++)
      out[i][j] = in[nskip+i][nskip+j];
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
     
/* Eigenvalue routines */
void eigsrt(double d[], double **v, int n)
{
        int k,j,i;
        double p;

        for (i=0;i<n-1;i++) {
                p=d[k=i];
                for (j=i+1;j<n;j++)
                        if (d[j] <= p) p=d[k=j];
                if (k != i) {
                        d[k]=d[i];
                        d[i]=p;
                        for (j=0;j<n;j++) {
                                p=v[j][i];
                                v[j][i]=v[j][k];
                                v[j][k]=p;
                        }
                }
        }
}
/* (C) Copr. 1986-92 Numerical Recipes Software vO(1Z712,$. */

#define NRANSI
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);
void jacobi1(double **a, int n, double d[], double **v, int *nrot)
{
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

	b=dvector(0,n-1);
	z=dvector(0,n-1);
	for (ip=0;ip<n;ip++) {
		for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=0;ip<n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=1;i<=75;i++) {
		sm=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			free_dvector(z,0,n-1);
			free_dvector(b,0,n-1);
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (double)(fabs(d[ip])+g) == fabs(d[ip])
					&& (double)(fabs(d[iq])+g) == fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((double)(fabs(h)+g) == (double)fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=0;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=0;j<n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	nrerror("Too many iterations in routine jacobi");
}
#undef ROTATE
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software vO(1Z712,$. */

/* Obtain eigenvalues & eigenvectors */
void jacobi(double **a, int n, double d[], double **v, int *nrot) {
  int i,j,vSign;

  jacobi1(a,n,d,v,nrot);
  eigsrt(d,v,n);

  /* Point eigenvectors to have first element positive */
  for(i=0;i<n;i++) {
    j=0;
    while (j<n-1 && v[j][i]==0) j++;
    vSign = v[j][i]>0? 1: -1;
    for(j=0;j<n;j++) v[j][i] *= vSign;
  }
}

/* Start our routines */

#define DELTA_A 0.025

/* DETF figure-of-merit & projection onto w0-wa:
 * projCov = {Var(w0), Cov, Var(wa)}.
 */
double get_fom_detf2(double **cov, int nw, double *projCov) {
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
  sym_matrix_invert(new_cov, new_cov, nw);
  for(i=0; i<nw; i++)
    new_cov[i][i] -= WBIN_PRIOR;
    
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
  sym_matrix_invert(F_wwa,F_wwa,2);
  area = sqrt(F_wwa[0][0]*F_wwa[1][1]-F_wwa[0][1]*F_wwa[0][1]);
#if 0
  printf("sig(w)=%9.6lf; sig(wa)=%9.6lf; corr=%9.6lf\n", sqrt(F_wwa[0][0]), sqrt(F_wwa[1][1]), F_wwa[0][1]/sqrt(F_wwa[0][0]*F_wwa[1][1]));
  printf("a_pivot=%9.6lf, sig(wp)=%9.6lf\n", 1.+F_wwa[0][1]/F_wwa[1][1],
    sqrt(F_wwa[0][0]-F_wwa[0][1]*F_wwa[0][1]/F_wwa[1][1]));
#endif
  projCov[0] = F_wwa[0][0];
  projCov[1] = F_wwa[1][0];
  projCov[2] = F_wwa[1][1];
  
  free_dmatrix(new_cov, 0, nw-1, 0, nw-1);
  free_dvector(vec1, 0, nw-1);
  free_dvector(vec2, 0, nw-1);
  free_dmatrix(F_wwa, 0,1,0,1);
  return(1./(area));
}

/* Gets Albrecht-Bernstein figure of merit for a covariance matrix of
 * w in redshift bins
 */
double get_fom_ab(double **cov, int nw) {
#ifdef EVOUTPUT
  FILE *fp;
#endif
  int i,j,nrot;
  double **new_cov;
  double **eigenvecs, *eigenvals;
  double fom,temp;
  double prod, sratio;

  new_cov=dmatrix(0, nw-1, 0, nw-1);
  eigenvecs=dmatrix(0, nw-1, 0, nw-1);
  eigenvals=dvector(0, nw-1);

  prod=sratio=1;

  /* Copy covariance matrix since NR jacobi destroys input matrix */
  for(i=0; i<nw; i++)
    for(j=0; j<nw; j++)
      new_cov[i][j] = cov[i][j];

  jacobi(new_cov, nw, eigenvals, eigenvecs, &nrot);
  printf("Eigenvalues:\n");
  for(i=0; i<nw && eigenvals[i]<MAX_EIGEN_PRINT; i++) {
    temp = 1./eigenvals[i]-WBIN_PRIOR;
    temp = 1./sqrt(40.*temp + 1e-12);
    sratio = sigmaS3[i]/temp;
    prod *= sratio;
    if (i<5) {
      printf("%2d %11.6lf ", i+1, temp);
      printf("%9.3lf %9.3le ", sratio, prod);
      for(j=0; j<nw; j++) printf("%c", eigenvecs[j][i]>0? '+': eigenvecs[j][i]<0? '-': 'o');
      printf("\n");
    }
  }

#ifdef EVOUTPUT
  fp = fopen("pc.dat", "w");
  for(j=0; j<nw; j++) {
    for(i=0; i<5; i++) fprintf(fp, " %10.7lf", eigenvecs[j][i]);
    fprintf(fp, "\n");
  }
  fclose(fp);
#endif

  fom = 1.;
  for(i=0; i<nw; i++)
    if (eigenvals[i]<36.)
      fom /= sqrt(eigenvals[i]/36.);

  free_dmatrix(new_cov, 0, nw-1, 0, nw-1);
  free_dmatrix(eigenvecs, 0, nw-1, 0, nw-1);
  free_dvector(eigenvals, 0, nw-1);
  return(fom);
}

/* Gets FoM for step function w(z) parameterization */
double get_fom_lohi(double **cov, int nw) {
  int i,j;
  double **new_cov;
  double Fll, Flh, Fhh, det;

  new_cov = dmatrix(0, nw-1, 0, nw-1);

  /* Copy covariance matrix since inversion destroys it */
  for(i=0; i<nw; i++)
    for(j=0; j<nw; j++)
      new_cov[i][j] = cov[i][j];
  sym_matrix_invert(new_cov, new_cov, nw);
  for(i=0; i<nw; i++)
    new_cov[i][i] -= WBIN_PRIOR;

  /* Get Fisher for lo-z (z<1) and hi-z (z>1) */
  Fll=Flh=Fhh=0.;
#define ICUT 20
  for(i=0;i<ICUT;i++) for(j=0;j<ICUT;j++) Fll+=new_cov[i][j];
  for(i=ICUT;i<nw;i++) for(j=0;j<ICUT;j++) Flh+=new_cov[i][j];
  for(i=ICUT;i<nw;i++) for(j=ICUT;j<nw;j++) Fhh+=new_cov[i][j];
#undef ICUT  

  det = Fll*Fhh-Flh*Flh;

  printf("LoHi FoM info: sig(w_lo) = %9.5lf sig(w_hi) = %9.5lf corr = %9.5lf\n",
    sqrt(Fhh/det), sqrt(Fll/det), -Flh/sqrt(Fll*Fhh));

  free_dmatrix(new_cov, 0, nw-1, 0, nw-1);

  /* Return inverse-area */
  return(sqrt(det));
}

/* Gets uncertainty on w in redshift range, bins nskip .. nskip+nuse-1.
 * Forces w to be constant in this range.
 */
double get_sigw(double **C, int nskip, int nuse) {
  double **subC;
  double ivar;
  int i,j;

  /* Construct sub-covariance matrix */
  subC = dmatrix(0,nuse-1,0,nuse-1);
  submatrix2(C,subC,nskip,nuse);

  /* Get Fisher via inversion */
  if (sym_matrix_invert(subC,subC,nuse)) {
    fprintf(stderr, "Failed: matrix inversion for w subrange.\n");
    exit(1);
  }

  /* Inverse-variance */
  ivar = -nuse * WBIN_PRIOR;
  for(i=0;i<nuse;i++) for(j=0;j<nuse;j++) ivar += subC[i][j];

  free_dmatrix(subC,0,nuse-1,0,nuse-1);
  return(1./sqrt(ivar));
}

/* Gets uncertainty on w in redshift range, bins nskip .. nskip+nuse-1.
 * Forces w to be constant in this range, and fixes w=-1 elsewhere.
 */
double get_sigw2(double **C, int nskip, int nuse) {
  double **F;
  double ivar;
  int i,j;

  /* Construct sub-covariance matrix */
  F = dmatrix(0,35,0,35);
  submatrix2(C,F,0,36);

  /* Get Fisher via inversion */
  if (sym_matrix_invert(F,F,36)) {
    fprintf(stderr, "Failed: matrix inversion for w.\n");
    exit(1);
  }

  /* Inverse-variance */
  ivar = -nuse * WBIN_PRIOR;
  for(i=0;i<nuse;i++) for(j=0;j<nuse;j++) ivar += F[nskip+i][nskip+j];

  free_dmatrix(F,0,35,0,35);
  return(1./sqrt(ivar));
}

/* MAIN PROGRAM.
 * Reads in first a number of Fisher matrices to read, then a list
 * of those Fishers.
 *
 * Preceding a Fisher by e.g. "\*2.0" multiplies it by 2.0.
 */
int main(int argc, char **argv) {
  int kF, NF;
  double **Fbig, **Fsmall, **Ftemp, **Cbig;
  FILE *fp;
  int i,j;
  double f, sigGamma;
  double fomDETF, fomGamma, fomAB, fomLH;
  double w0waCov[3];
  char InLine[1024];
  double mFactor;
  double temp;

  Cbig = dmatrix(0,44,0,44);
  Fbig = dmatrix(0,44,0,44);
  Fsmall = dmatrix(0,35,0,35);

  NF = argc-1;
  kF=0;
  while(kF<NF) {
    mFactor = 1.;

    /* See if the Fisher is supposed to be re-scaled */
    if (argv[kF+1][0] == '*') {
      sscanf(argv[kF+1]+1, "%lg", &mFactor);
      kF++;
    }

    Ftemp = dmatrix(0,44,0,44);
    fp = fopen(argv[kF+1], "r");
    while(fgets(InLine,1023,fp)!=NULL) if (InLine[0]!='#') {
      if(sscanf(InLine, "%d %d %lg", &i, &j, &f)==3)
        Ftemp[i][j] = Ftemp[j][i] = f * mFactor;
    }
    fclose(fp);
    for(i=0;i<45;i++) for(j=0;j<45;j++) Fbig[i][j] += Ftemp[i][j];
    free_dmatrix(Ftemp,0,44,0,44);
    kF++;
  }

  printf("sigma(G0|everything) = %8.5lf\n", 1./sqrt(Fbig[7][7]));

  /* Get giant Covariance matrix */
  for(i=0;i<45;i++) Fbig[i][i] += 1.e-10*(1.+Fbig[i][i]);
  for(i=9;i<45;i++) Fbig[i][i] += WBIN_PRIOR;
  if (sym_matrix_invert(Fbig,Cbig,45)) nrerror("Error: non semipositive definite Cov [#1]\n");
  submatrix2(Cbig,Fsmall,9,36);

#if 1
  /* Uncertainties in other parameters */
  /* Hubble */
  temp = 0.; for(i=1;i<=4;i++) for(j=1;j<=4;j++) if (i!=2 && j!=2) temp+=Cbig[i][j]; temp=sqrt(temp); /* sigma(h^2) */
  printf("sigma(ln H0) = %8.6lf sigma(omega_m) = %8.6lf sigma(omega_de) = %8.6lf\n",
    temp/(2.*0.719*0.719), sqrt(Cbig[1][1]), sqrt(Cbig[4][4]));
  printf("COV:\n%12.5le %12.5le %12.5le\n%12.5le %12.5le %12.5le\n%12.5le %12.5le %12.5le\n",
    Cbig[1][1], Cbig[1][3], Cbig[1][4], Cbig[3][1], Cbig[3][3], Cbig[3][4], Cbig[4][1], Cbig[4][3], Cbig[4][4]);
#endif

  /* growth of structure */
  fomGamma = 1./Cbig[5][5];

  /* More parameters */
  printf("sigma(Omega_K) = %9.6lf; sigma(ns) = %9.6lf\n", sqrt(Cbig[3][3])/.719/.719, sqrt(Cbig[0][0]));

  /* high-z growth */
  submatrix(Fbig,Fsmall,9);
  if (sym_matrix_invert(Fsmall,Cbig,9)) nrerror("Error: non semipositive definite Cov [#2]\n");
  sigGamma = sqrt(Cbig[5][5]);

  /* project growth of structure out of the Fisher matrix */
  Fbig[5][5] += 1.e24;
  Fbig[7][7] += 1.e24;
  if (sym_matrix_invert(Fbig,Cbig,45)) nrerror("Error: non semipositive definite Cov [#4]\n");
  submatrix2(Cbig,Fsmall,9,36);

  /* SN calibration */
  printf("sigma(SN mag | assume GR) = %7.5lf\n", sqrt(Cbig[6][6]));

  /* w0,wa */
  fomDETF = get_fom_detf2(Fsmall,36,w0waCov);
  printf("sigma(wp) = %8.6lf; sigma(wa) = %8.5lf; zp = %8.5lf\n",
    sqrt(w0waCov[0]-w0waCov[1]*w0waCov[1]/w0waCov[2]), sqrt(w0waCov[2]), 1./(1.+w0waCov[1]/w0waCov[2])-1.);
#if 0
  printf("ssc: %8.6lf %8.6lf %9.6lf\n", sqrt(w0waCov[0]), sqrt(w0waCov[2]), w0waCov[1]/sqrt(w0waCov[0]*w0waCov[2]));
#endif

  fomLH = get_fom_lohi(Fsmall,36);

#if 0
  /* w in redshift bins */
  printf("sig w: [0.00-0.48] %9.6lf [0.48-1.00] %9.6lf [1.00-1.50] %9.6lf [1.50-2.08] %9.6lf [2.08-3.00] %9.6lf [3.00+] %9.6lf\n",
    get_sigw(Fsmall,0,13), get_sigw(Fsmall,13,7), get_sigw(Fsmall,20,4), get_sigw(Fsmall,24,3), get_sigw(Fsmall,27,3), get_sigw(Fsmall,30,6));
  printf("sig w: [0.00-1.00] %9.6lf [1.00-2.08] %9.6lf [2.08+] %9.6lf\nsig w: [1.00+] %9.6lf [1.50+] %9.6lf\n",
    get_sigw(Fsmall,0,20), get_sigw(Fsmall,20,7), get_sigw(Fsmall,27,9), get_sigw(Fsmall,20,16), get_sigw(Fsmall,24,12));
  printf("sig w: [0.00-0.18] %9.6lf [0.18-0.43] %9.6lf [0.43-0.82] %9.6lf [0.82-1.50] %9.6lf [1.50-3.00] %9.6lf [3.00+] %9.6lf\n",
    get_sigw(Fsmall,0,6), get_sigw(Fsmall,6,6), get_sigw(Fsmall,12,6), get_sigw(Fsmall,18,6), get_sigw(Fsmall,24,6), get_sigw(Fsmall,30,6));
#endif
  printf("sig(w@z>0.48) = %9.6lf [%9.6lf]\n", get_sigw2(Fsmall,13,23), get_sigw(Fsmall,13,23));
  printf("sig(w@z>1.00) = %9.6lf [%9.6lf]\n", get_sigw2(Fsmall,20,16), get_sigw(Fsmall,20,16));
  printf("sig(w@z>2.08) = %9.6lf [%9.6lf]\n", get_sigw2(Fsmall,27,9) , get_sigw(Fsmall,27,9) );

  /* AB */
  fomAB = get_fom_ab(Fsmall,36);

  /* Report Figures of Merit */
  printf("FOM = %7.2lf [DETF] %7.2lf [LH] %7.2lf [gamma] %6ld [w0-wa-gamma]\n", fomDETF, fomLH, fomGamma, (long)floor(0.5 + fomDETF/sigGamma));

  free_dmatrix(Cbig,0,44,0,44);
  free_dmatrix(Fbig,0,44,0,44);
  free_dmatrix(Fsmall,0,35,0,35);
  return(0);
}
