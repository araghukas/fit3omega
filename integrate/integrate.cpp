#include "integrate.h"


std::complex<double> Integral::fB(int i_layer, double lambda, double omega)
{
  int i = i_layer - 1;

  double x = s.kxs[i] / s.kys[i] * lambda * lambda;
  double y = 2.0 * omega * s.Cvs[i] / s.kys[i];
  std::complex<double> z {x, y};
  return std::sqrt(z);
}

std::complex<double> Integral::phi(int i_layer, double lambda, double omega)
{
  return fB(i_layer, lambda, omega) * s.ds[i_layer - 1];
}

std::complex<double> Integral::fA(int i_layer, double lambda, double omega)
{
  if (i_layer == s.nl)
    switch (b_type) {
      case SEMI_INFINITE:
        return std::complex<double> {-1, 0};
      case ADIABATIC:
        return std::tanh(phi(s.nl, lambda, omega));
      case ISOTHERMAL:
        return  -1.0 / std::tanh(phi(s.nl, lambda, omega));
      default:
        throw std::invalid_argument("unrecognized boundary type");
    }
  int i = i_layer - 1;
  std::complex<double> A_next = fA(i_layer + 1, lambda, omega);
  std::complex<double> B_next = fB(i_layer + 1, lambda, omega);
  std::complex<double> B = fB(i_layer, lambda, omega);
  std::complex<double> kBA_term = A_next * s.kys[i + 1]* B_next / (s.kys[i] * B);
  std::complex<double> tanh_term = std::tanh(phi(i_layer, lambda, omega));
  return (kBA_term - tanh_term) / (1.0 - kBA_term * tanh_term);
}

std::complex<double> Integral::integrand(double lambda, double omega)
{
  std::complex<double> A1 = fA(1, lambda, omega);
  std::complex<double> B1 = fB(1, lambda, omega);
  double sinc_term = sinc(s.b * lambda);

  return 1.0 / (A1 * B1) * sinc_term * sinc_term;
}

std::complex<double>* Integral::integrate()
{
  double h = (lambda_f - lambda_i) / N;
  for (int k = 0; k < s.nf; k++) {
    std::complex<double> fi = integrand(lambda_i, s.omegas[k]);
    std::complex<double> ff = integrand(lambda_f, s.omegas[k]);
    result[k] = 0.5 * h * (fi + ff);
  }
  for (int k = 0; k < s.nf; k++)
    for (int j = 0; j < N; j++)
      result[k] += h * integrand(lambda_i + j * h, s.omegas[k]);

  return result;
}


double sinc(double x)
{
  return std::sin(x) / x;
}

/*
int main()
{
  const double omegas[] {1, 10, 100, 1000, 10000};
  const double ds[] {3e-6, 300e-6};
  const double kxs[] {0.25, 150};
  const double kys[] {0.25, 150};
  const double Cvs[] {100000, 2000000};
  const double b = 60e-6;

  const int nf = 5;
  const int nl = 2;

  const double lambda_i = 1e-3;
  const double lambda_f = 1e7;
  const int N = 1500;
  const char b_type = 'a';


  Sample s {nf, nl, omegas, ds, kxs, kys, Cvs, b};
  Integral g {lambda_i, lambda_f, N, s, b_type};
  
  std::complex<double>* r = g.integrate();
  for (int i = 0; i < nf; i++)
    std::cout << r[i].real() << ", " << r[i].imag() << std::endl;
}
*/
