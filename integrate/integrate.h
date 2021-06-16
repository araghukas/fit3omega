#include <math.h>
#include <complex.h>

#define MAX_N_LAYERS 10
#define MAX_N_OMEGAS 150
#define N_XPTS    200

// measurement domain; array of angular frequencies
const double *OMEGAS;
int n_OMEGAS;


// utility functions
void make_logspace(double *arr, double min, double max, int size);
double sinc_sq(double x);


// sample parameters (set from Python side)
int N_LAYERS;  // number of layers
double HALF_WIDTH;      // heater half-width

double ds_[MAX_N_LAYERS];   // heights
double kxs_[MAX_N_LAYERS];  // x-thermal conductivity
double kys_[MAX_N_LAYERS];  // y-thermal conductivity
double Cvs_[MAX_N_LAYERS];  // heat capacities
double Rcs_[MAX_N_LAYERS];  // thermal contact resistances  


// =================================================================================================
// Borca-Tasciuc model (no contact resistance)
// =================================================================================================
char boundary_type_;
double complex bt_integral_result_[MAX_N_OMEGAS];
double complex bt_integrand(double lambda, double omega, double *ds, double *kxs, double *kys, double *Cvs);
double complex *bt_integral(double *ds, double *kxs, double *kys, double *Cvs);
double LAMBDAS[N_XPTS];

// =================================================================================================
// Olson, Graham, and Chen model
// =================================================================================================
double complex fz_[MAX_N_OMEGAS];
double complex fPhi_[MAX_N_OMEGAS];
double complex ogc_integral_result_[MAX_N_OMEGAS];
double complex ogc_integrand(double lambda, double omega,
	                         double *ds, double *kxs, double *kys, double *Cvs, double *Rcs);
double complex *ogc_integral(double *ds, double *kxs, double *kys, double *Cvs, double *Rcs);
double CHIS[N_XPTS];