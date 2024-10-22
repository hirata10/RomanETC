double DNget_growth_normhiz(COSMOPARAM *p, double z);

double getSoundHorizon(COSMOPARAM *p);

double DNeh_trf_mpc(COSMOPARAM *p, double k);

double get_A0_BAO(COSMOPARAM *p);

double get_Om(COSMOPARAM *p, double z);

void get_Delta2k(COSMOPARAM *p, double z, double k0, double dlnk, long nk, double *D2);

double get_linsigma(COSMOPARAM *p, double z, double R);

void get_Delta2kNL(COSMOPARAM *p, double z, double k0, double dlnk, long nk, double *D2);

double get_nu(COSMOPARAM *p, double M, double z);

double get_fM(COSMOPARAM *p, double M, double z);

double get_bM(COSMOPARAM *p, double M, double z);

double getconc(COSMOPARAM *p, double M, double z);

double get_delta_virial(COSMOPARAM *p, double M, double z);

double mso2mvir(COSMOPARAM *p, double M, double z, double Delta);

double tkhalo(COSMOPARAM *p, double M, double z, double k);

double get_mcut(double n, double z, COSMOPARAM *p);

void getbr_halo(double Mcut, double z, double k, double po, COSMOPARAM *p, double *b, double *r);

double get_lensing_strength(COSMOPARAM *p, double z1, double z2);

void get_wl_ps(COSMOPARAM *p, long nl, double l0, double dlnl,
  double z1, double z2, double *Cl);

void get_wl_ps_many(COSMOPARAM *p, long nl, double l0, double dlnl,
  long nz, double *z, double *noise, double fsky, double *Cl,
  double **ClCov, double *auxparam);

void get_fisher_ccc1(COSMOPARAM *p, long nl, double l0, double dlnl,
  long nz, double *z, double *noise, double fsky, double *Cl, double rcorr,
  long b1, double **Fgg);

void get_isw(COSMOPARAM *p, double l, double dlnl,
  double z, double dz, double b2n, double fsky, double *dT,
  double *sigT);

void get_kappa_cmblens(COSMOPARAM *p, double l, double dlnl,
  double z, double dz, double b2n, double fsky, double *dkappa,
  double *sigk);

double get_VPS(COSMOPARAM *p, double z, double k);
