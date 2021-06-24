#include <math.h>
#include <complex.h>

#define MAX_n_LAYERS 10
#define MAX_n_OMEGAS 150
#define N_XPTS       200  // number of x sample points for integrations


// measurement domain; array of angular frequencies
const double *OMEGAS;
int n_OMEGAS;

// maximum number of fitting parameters (for either method)
const int MAX_n_PARAMS = 4 * MAX_n_LAYERS;

// utility functions
void make_logspace(double *arr, double min, double max, int size);
double sinc_sq(double x);
double complex *trapz(double complex (*fp)(double,double), double *xs, double complex *Fs);
double complex *val_trapz(double complex *fs, double *xs, double complex *Fs);


// sample parameters (set from Python side)
int n_LAYERS;       // number of layers
double HALF_WIDTH;  // heater half-width

// 5 layer parameters in general
double ds_[MAX_n_LAYERS];   // heights
double psis_[MAX_n_LAYERS];  // (x/y)-thermal conductivity ratio
double kys_[MAX_n_LAYERS];  // y-thermal conductivity
double Cvs_[MAX_n_LAYERS];  // heat capacities
double Rcs_[MAX_n_LAYERS];  // thermal contact resistances
