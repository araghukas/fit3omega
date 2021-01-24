#include <complex>
#include <iostream>


const char SEMI_INFINITE = 's';
const char ADIABATIC = 'a';
const char ISOTHERMAL = 'i';


struct Sample {
  const int nf; // number of frequencies in `omegas`
  const int nl; // number of layers on sample
  
  const double* omegas;
  
  const double* ds;   // layer thicknesses
  const double* kxs;  // x-thermal conductivity
  const double* kys;  // y-thermal conductivity
  const double* Cvs;  // heat capacities
  
  const double b; // heater half-width
};


class Integral {
  // the integral from Borca-Tasciuc et al. Eq (1-4)  
  const double lambda_i; // lower integration limit
  const double lambda_f; // upper integration limit
  const int N;           // number of sample points
  const Sample& s;       // sample properties
  const char b_type;           // boundary type

  std::complex<double> fB(int i_layer, double lambda, double omega);
  std::complex<double> phi(int i_layer, double lambda, double omega);
  std::complex<double> fA(int i_layer, double lambda, double omega);
  std::complex<double> integrand(double lambda, double omega);

  std::complex<double>* result;

public:
  Integral(double lambda_i, double lambda_f, int N, Sample& s, char b_type)
    : lambda_i{lambda_i}, lambda_f{lambda_f}, N{N}, s{s}, b_type{b_type}
  {
    result = new std::complex<double>[s.nf];
  }

  ~Integral()
  {
    delete result;
  }

  std::complex<double>* integrate();
};


double sinc(double x); // the sinc(x) function returns sin(x) / x
