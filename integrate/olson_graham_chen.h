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
double complex dz0_dkx_[MAX_N_LAYERS];
double complex dz0_dRc_[MAX_N_LAYERS];

/*
	{{ dZ/dky_0, ... , dZ/dky_N-1 },
	 { dZ/dCv_0, ... , dZ/dCv_N-1 },
	 { dZ/dkx_0, ... , dZ/dkx_N-1 },
	 { dZ/dRc_0, ... , dZ/dRc_N-1 }}
*/
