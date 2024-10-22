/* DE parameterization */
#ifdef SWG36
#define NWPARAM 36
#define DELTA_A 0.025
#define NPARAMTOT 50
#else
#define NWPARAM 9
#define DELTA_A 0.08888888889
#define NPARAMTOT 23
#endif

#define NSWG 45

/* Speed of light / 100 km/s */
#define HL_HIMPC 2997.92458

typedef struct {
  double T_cmb, SNflux, tau;
  double obh2, omh2;
  double lndelta_zeta, ns, alphas;
  double h, Omega_K;
  double growthAmp, growthSlope, growthCurve;
  double onuh2;
  double w[NWPARAM];
} COSMOPARAM;

/* Tomography bin pair indices */
#define TINDEX(b1,b2)  ((b1)<=(b2)? (((b2)+1)*(b2))/2+(b1): (((b1)+1)*(b1))/2+(b2))
#define NCROSS(zb)     (((zb)*((zb)+1))/2)

/* Parameter mappings */

/* Parameter mapping: returns pointer to i_th parameter of p.
 * Change this if you want to change the parameter mappings.
 */
#ifdef DEF_PARAMETER_MAPPINGS_HERE
double *get_ptr_param(COSMOPARAM *p, int i) {
  if (i<NWPARAM) return(p->w+i);
  i-=NWPARAM;
  switch(i) {
    case 0:
      return(&(p->T_cmb       ));
      break;
    case 1:
      return(&(p->obh2        ));
      break;
    case 2:
      return(&(p->omh2        ));
      break;
    case 3:
      return(&(p->lndelta_zeta));
      break;
    case 4:
      return(&(p->ns          ));
      break;
    case 5:
      return(&(p->alphas      ));
      break;
    case 6:
      return(&(p->h           ));
      break;
    case 7:
      return(&(p->Omega_K     ));
      break;
    case 8:
      return(&(p->tau         ));
      break;
    case 9:
      return(&(p->SNflux      ));
      break;
    case 10:
      return(&(p->growthAmp   ));
      break;
    case 11:
      return(&(p->growthSlope ));
      break;
    case 12:
      return(&(p->growthCurve ));
      break;
    case 13:
      return(&(p->onuh2       ));
      break;
    default:
      fprintf(stderr, "get_ptr_param: Invalid parameter\n");
      exit(1);
      break;
  }
}
#else
double *get_ptr_param(COSMOPARAM *p, int i);
#endif

#ifdef SWG36

#define DEFAULT_VARIATION {.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1, .02,.0002,.001,.02,.005,.001,.01,.005,.01,.01,.005,.005,.005,.002}
#define LARGE_VARIATION {10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10, 1,1,1,1,1,1,1,1,1,1,1,1,1,.1}

#else

#define DEFAULT_VARIATION {.1,.1,.1,.1,.1,.1,.1,.1,.1, .02,.0002,.001,.02,.005,.001,.01,.005,.01,.01,.01,.01,.01,.002}
#define LARGE_VARIATION {10,10,10,10,10,10,10,10,10, 1,1,1,1,1,1,1,1,1,1,1,1,1,.1}

#endif

/* Other functions */

void set_default_cosmology(COSMOPARAM *p);

double get_H(double z, COSMOPARAM *p);

double get_DAC(double z, COSMOPARAM *p);

double get_DRC(double z, COSMOPARAM *p);

void get_derivatives(double z, COSMOPARAM *p, int i, double *dlnHdpi, double *dlnDACdpi, double *dlnsdpi);

void gaussjinv(double **A, int n);

#define float double
void jacobi(float **a, int n, float d[], float **v, int *nrot);
#undef float

int sym_matrix_invert(double **C, double **Cinv, int N);
