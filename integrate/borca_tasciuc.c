#include "integrate.h"


// =================================================================================================
//
// IMPLEMENTS INTEGRAL AND INTEGRAND FROM BORCA-TASCIUC Eq. (1) for complex ΔT
//
// [REF: Rev. Sci. Instrum., Vol 72, No. 4, April 2001]
// =================================================================================================


double complex fB(int i_layer, double lambda, double omega, double *kxs, double *kys, double *Cvs)
{
  /* Borca-Tasciuc Eq. (3) */
  return csqrt(kxs[i_layer]/kys[i_layer]*lambda*lambda + I*2.0*omega*Cvs[i_layer]/kys[i_layer]);
}

double complex fA(int i_layer, double lambda, double omega, double *ds, double *kxs, double *kys, double *Cvs)
{
  /* Borca-Tasciuc Eq. (2); with i -> i+1 */
  // base case - depends on boundary type
  if (i_layer == N_LAYERS - 1) {
    switch (boundary_type_)
    {
      case 's':
        return -1.0;
      case 'a':
      {
        complex double Bn = fB(i_layer,lambda,omega,kxs,kys,Cvs);
        return -1.0 * ctanh(Bn * ds[i_layer]);
      }
      case 'i':
      {
        complex double Bn = fB(i_layer,lambda,omega,kxs,kys,Cvs);
        return -1.0 / ctanh(Bn * ds[i_layer]);
      }
      default:
        return -1.0;
    }
  }

  // recursive call until base case
  double complex A_ii = fA(i_layer+1,lambda,omega,ds,kxs,kys,Cvs);

  double complex B_i = fB(i_layer,lambda,omega,kxs,kys,Cvs);
  double complex B_ii = fB(i_layer+1,lambda,omega,kxs,kys,Cvs);
  double k_i = kys[i_layer];
  double k_ii = kys[i_layer+1];
  double d_i = ds[i_layer];
  
  double complex AkB = A_ii * k_ii * B_ii / (k_i * B_i);
  double complex tanh_term = ctanh(B_i * d_i);
  return (AkB - tanh_term) / (1 - AkB * tanh_term);
}


// =================================================================================================


double complex bt_integrand(double lambda, double omega, double *ds, double *kxs, double *kys, double *Cvs)
{
  /* Borca-Tasciuc Eq.(1) integrand */
  double complex A_top = fA(0,lambda,omega,ds,kxs,kys,Cvs);
  double complex B_top = fB(0,lambda,omega,kxs,kys,Cvs);
  return 1.0 / (A_top * B_top) * sinc_sq(HALF_WIDTH * lambda);
}


double complex *bt_integral(double *ds, double *kxs, double *kys, double *Cvs)
{
  /* Borca-Tasciuc Eq.(1) integral - here, trapezoidal integration in log-space */
  for (int i = 0; i < n_OMEGAS; i++) {

    bt_integral_result_[i] = 0.0*I;
    double complex f0 = bt_integrand(LAMBDAS[0],OMEGAS[i],ds,kxs,kys,Cvs);
    double complex f_prev = f0;

    for (int k = 1; k < N_XPTS; k++) {

      double complex fk = bt_integrand(LAMBDAS[k],OMEGAS[i],ds,kxs,kys,Cvs);
      double dx = LAMBDAS[k] - LAMBDAS[k-1];
      bt_integral_result_[i] += (dx / 2.0) * (fk + f_prev);
      f_prev = fk;

    }
  }
  return bt_integral_result_;  // pointer to results for each ω
}
