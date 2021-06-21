#include <math.h>
#include <complex.h>

#define MAX_N_LAYERS 10
#define MAX_N_OMEGAS 150
#define N_XPTS       200  // number of x sample points for integrations


// measurement domain; array of angular frequencies
const double *OMEGAS;
int n_OMEGAS;


// utility functions
void make_logspace(double *arr, double min, double max, int size);
double sinc_sq(double x);
double complex *trapz(double complex (*fp)(double,double), double *xs, double complex *Fs);


// sample parameters (set from Python side)
int N_LAYERS;  // number of layers
double HALF_WIDTH;      // heater half-width

// 5 layer parameters in general
double ds_[MAX_N_LAYERS];   // heights
double kxs_[MAX_N_LAYERS];  // x-thermal conductivity
double kys_[MAX_N_LAYERS];  // y-thermal conductivity
double Cvs_[MAX_N_LAYERS];  // heat capacities
double Rcs_[MAX_N_LAYERS];  // thermal contact resistances
