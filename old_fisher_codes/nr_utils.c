#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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

/* Eigenvalue routines */
#include "jacobi.c"
