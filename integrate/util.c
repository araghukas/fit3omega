#include "integrate.h"


// writes into array points linear in log-space, between min and max
void make_logspace(double *arr, double min, double max, int size)
{
	for (int k = 0; k < size; k++)
    	arr[k] = min * pow(max / min, ((double) k) / (size - 1));
}


// a function that computes (sinc(x))^2
double sinc_sq(double x)
{
	double sinc = sin(x) / x;
	return sinc*sinc;
}


// a general trapezoidal-rule integrator of f(x,ω_i)dx, for each ω_i 
double complex *trapz(double complex (*fp)(double,double),
                 	  double *xs,
               	      double complex *Fs)
{
	for (int i = 0; i < n_OMEGAS; i++) {
		Fs[i] = 0.0*I;
		double complex f0 = fp(xs[0],OMEGAS[i]);
		double complex f_prev = f0;
		for (int k = 1; k < N_XPTS; k++) {
			double complex fk = fp(xs[k],OMEGAS[i]);
			double dx = xs[k] - xs[k-1];
			Fs[i] += (dx / 2.0) * (fk + f_prev);
			f_prev = fk;
		}
	}

	return Fs;
}
