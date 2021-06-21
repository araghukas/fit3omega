// =================================================================================================
// Borca-Tasciuc model (no contact resistance)
// =================================================================================================
char boundary_type_;
double LAMBDAS[N_XPTS];

double complex bt_integrand_x_[N_XPTS];
double complex bt_integrand(double lambda, double omega);


double complex bt_integral_result_[MAX_N_OMEGAS];
double complex *bt_integral(void);
