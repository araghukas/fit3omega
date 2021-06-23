// =================================================================================================
// Olson, Graham, and Chen model
// =================================================================================================

// OGC model for surface impedance, Z
double CHIS[N_XPTS];
double complex Phi_[MAX_N_LAYERS];
double complex *fPhi(double chi, double omega);

double complex zs_[MAX_N_LAYERS];
double complex *fzs(double chi, double omega);

double complex ogc_integrand(double lambda, double omega);
double complex ogc_integral_result_[MAX_N_OMEGAS];
double complex *ogc_integral(void);


// OGC model for dZ/dX_k, for sample parameter X_k
double complex Xis_[MAX_N_LAYERS];
double complex *fXis(double chi, double omega, double complex *Phis, double complex *zs);

double complex dz0_dky_[MAX_N_LAYERS];
double complex dz0_dCv_[MAX_N_LAYERS];
double complex dz0_dpsi_[MAX_N_LAYERS];
double complex dz0_dRc_[MAX_N_LAYERS];
