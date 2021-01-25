#include "intg.h"


complex<double>* IntegralTermBT_EQ1::fB(int i_layer, double lambda)
{
  double ky = kys[i_layer - 1];
  double kx = kxs[i_layer - 1];
  double Cv = Cvs[i_layer - 1];

  for (int i = 0; i < nf; i++) {
    double re = kx / ky * lambda * lambda;
    double im = 2.0 * omegas[i] * Cv / ky;
    fB_omegas_[i] = sqrt(complex<double> {re, im});
  }

  return fB_omegas_;
}


complex<double>* IntegralTermBT_EQ1::phi(int i_layer, double lambda)
{
  double d = ds[i_layer - 1];
  complex<double>* fB_omegas_ = fB(i_layer, lambda);

  for (int i = 0; i < nf; i++)
    phi_omegas_[i] = fB_omegas_[i] * d;

  return phi_omegas_;
}


complex<double>* IntegralTermBT_EQ1::fA(int i_layer, double lambda)
{
  if (i_layer == nl) {
    switch (b_type) {
    case 's':
      {
        for (int i = 0; i < nf; i++)
          fA_omegas_[i] = complex<double>{-1.0, 0.0};
        break;
      }
    case 'i':
      {
        complex<double>* phi_omegas_ = phi(nl, lambda);
        for (int i = 0; i < nf; i++)
          fA_omegas_[i] = complex<double>{-1.0, 0.0} * tanh(phi_omegas_[i]);
        break;
      }
    case 'a':
      {
        complex<double>* phi_omegas_ = phi(nl, lambda);
        for (int i = 0; i < nf; i++)
          fA_omegas_[i] = complex<double>{-1.0, 0.0} / tanh(phi_omegas_[i]);
        break;
      }
    default:
      throw invalid_argument("unknown boundary specifier");
    }

    return fA_omegas_;
  }

  complex<double>* A_next = fA(i_layer + 1, lambda);
  complex<double>* B_next = fB(i_layer + 1, lambda);
  for (int i = 0; i < nf; i++)
    _AB_next[i] = A_next[i] * B_next[i];

  complex<double>* B = fB(i_layer, lambda);
  double k_next = kys[i_layer + 1];
  double k = kys[i_layer];
  for (int i = 0; i < nf; i++)
    _kkB[i] = k_next / (k  * B[i]);

  complex<double>* phi_ = phi(i_layer, lambda);
  for (int i = 0; i < nf; i++)
    _tanh_term[i] = tanh(phi_[i]);

  for (int i = 0; i < nf; i++) {
    complex<double> ABkkB = _AB_next[i] * _kkB[i];
    complex<double> tanh_term = _tanh_term[i];
    complex<double> z {1, 0};
    fA_omegas_[i] = (ABkkB - tanh_term) / (z - ABkkB * tanh_term);
  }

  return fA_omegas_;
}


complex<double> IntegralTermBT_EQ1::sinc_sq(double x)
{
  double sinc = sin(x) / x;
  return complex<double>{sinc * sinc, 0.0};
}


complex<double>* IntegralTermBT_EQ1::integrand(double lambda)
{
  complex<double>* A1 = fA(1, lambda);
  complex<double>* B1 = fB(1, lambda);
  complex<double> sinc_sq_ = sinc_sq(b * lambda);
  for (int i = 0; i < nf; i++)
    _integrand[i] = sinc_sq_ / (A1[i] * B1[i]);
  return _integrand;
}


complex<double>* IntegralTermBT_EQ1::integrand(int nl,
                                               double lambda,
                                               double* d_,
                                               double* kx_,
                                               double* ky_,
                                               double* Cv_,
                                               char b_type)
{
  this->nl = nl;
  ds  = d_;
  kxs = kx_;
  kys = ky_;
  Cvs = Cv_;
  this->b_type = b_type;

  complex<double>* A1 = fA(1, lambda);
  complex<double>* B1 = fB(1, lambda);
  complex<double> sinc_sq_ = sinc_sq(b * lambda);
  for (int i = 0; i < nf; i++)
    _integrand[i] = complex<double>{1.0, 0.0} * sinc_sq_ / (A1[i] * B1[i]);
  return _integrand;
}


complex<double>* IntegralTermBT_EQ1::integral(int nl,
                                              double* d_,
                                              double* kx_,
                                              double* ky_,
                                              double* Cv_,
                                              char b_type)
{
  this->nl  = nl;
  ds  = d_;
  kxs = kx_;
  kys = ky_;
  Cvs = Cv_;
  this->b_type = b_type;


  // get 1D logarithmic grid of `N` sample points
  double* lambdas = new double[N];
  for (int k = 0; k < N; k++)
    lambdas[k] = lambda_i * pow((lambda_f / lambda_i), double(k) / (N - 1));

  // prime function values (for each omega) at lambdas[0]
  complex<double>* f_prev = new complex<double>[nf];
  complex<double>* f0 = integrand(lambdas[0]);
  for (int i = 0; i < nf; i++) {
    result[i] = complex<double>{0.0, 0.0};
    f_prev[i] = f0[i];
  }

  // accumulate sums (for each omega)
  for (int k = 1; k < N; k++) {
    complex<double>* fk = integrand(lambdas[k]);
    double dx = lambdas[k] - lambdas[k - 1];
    for (int i = 0; i < nf; i++) {
      result[i] += complex<double>{dx / 2.0, 0.0} * (fk[i] + f_prev[i]);
      f_prev[i] = fk[i];
    }
  }

  delete[] lambdas;
  delete[] f_prev;

  return result;
}
