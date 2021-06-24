// =================================================================================================
// Olson, Graham, and Chen model
// =================================================================================================

// OGC model for surface impedance, Z
double CHIS[N_XPTS];
double complex Phis_[MAX_n_LAYERS];
void fPhis(double chi, double omega);

double complex zs_[MAX_n_LAYERS];
void fzs(double chi, double omega);

double complex ogc_integrand(double lambda, double omega);
double complex ogc_integral_result_[MAX_n_OMEGAS];
double complex *ogc_integral(void);


// OGC model for dZ/dX_k, for sample parameter X_k
double complex Xis_[MAX_n_LAYERS];
void fXis(double chi, double omega);

double complex dz0_dky_[MAX_n_LAYERS];
double complex dz0_dCv_[MAX_n_LAYERS];
double complex dz0_dpsi_[MAX_n_LAYERS];
double complex dz0_dRc_[MAX_n_LAYERS];

int n_PARAMS;
/*
indices identifying the degree of freedom, from

	{{   k[0], ... , k[n_LAYERS-1]   },
	 {  Cv[0], ... , Cv[n_LAYERS-1]  },
	 { psi[0], ... , psi[n_LAYERS-1] },
	 {  Rc[0], ... , Rc[n_LAYERS-1]  }}

*/
int param_ids_[MAX_n_PARAMS][2];

double complex jac_Z_result_[MAX_n_PARAMS][MAX_n_OMEGAS];
double complex jac_Z_fs_buff_[MAX_n_PARAMS][N_XPTS];
double complex (*jac_Z(void))[MAX_n_OMEGAS];