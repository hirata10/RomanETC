void nrerror(char error_text[]);

double *dvector(long nl, long nh);

void free_dvector(double *v, long nl, long nh);

double **dmatrix(long nrl, long nrh, long ncl, long nch);

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

int *ivector(long nl, long nh);

void free_ivector(int *v, long nl, long nh);

void rk4_dbl(double y[], int n, double x, double h, double yout[],
  void (*derivs)(double, double [], double []));

void submatrix(double **in, double **out, int dim);
