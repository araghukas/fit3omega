#include "integrate.h"
#include "borca_tasciuc.h"


// =================================================================================================
//
// IMPLEMENTS INTEGRAL AND INTEGRAND FROM BORCA-TASCIUC Eq. (1) for complex ΔT
//
// [REF: Rev. Sci. Instrum., Vol 72, No. 4, April 2001]
// =================================================================================================


double complex fB(int i_layer, double lambda, double omega)
{
  /* Borca-Tasciuc Eq. (3) */
  return csqrt(kxs_[i_layer]/kys_[i_layer]*lambda*lambda
  	           + I*2.0*omega*Cvs_[i_layer]/kys_[i_layer]);
}

double complex fA(int i_layer, double lambda, double omega)
{
  /* Borca-Tasciuc Eq. (2); with i -> i+1 */
  // base case - depends_ on boundary type
  if (i_layer == N_LAYERS - 1) {
    switch (boundary_type_)
    {
      case 's':
        return -1.0;
      case 'a':
      {
        complex double Bn = fB(i_layer,lambda,omega);
        return -1.0 * ctanh(Bn * ds_[i_layer]);
      }
      case 'i':
      {
        complex double Bn = fB(i_layer,lambda,omega);
        return -1.0 / ctanh(Bn * ds_[i_layer]);
      }
      default:
        return -1.0;
    }
  }

  // recursive call until base case
  double complex A_ii = fA(i_layer+1,lambda,omega);

  double complex B_i = fB(i_layer,lambda,omega);
  double complex B_ii = fB(i_layer+1,lambda,omega);
  double k_i = kys_[i_layer];
  double k_ii = kys_[i_layer+1];
  double d_i = ds_[i_layer];
  
  double complex AkB = A_ii * k_ii * B_ii / (k_i * B_i);
  double complex tanh_term = ctanh(B_i * d_i);
  return (AkB - tanh_term) / (1 - AkB * tanh_term);
}


// =================================================================================================


double complex bt_integrand(double lambda, double omega)
{
  /* Borca-Tasciuc Eq. (1) integrand */
  double complex A_top = fA(0,lambda,omega);
  double complex B_top = fB(0,lambda,omega);
  return sinc_sq(HALF_WIDTH * lambda) / (A_top * B_top);
}


double complex *bt_integral()
{
	/* Borca-Tasciuc Eq. (1) integral */
	return trapz(bt_integrand,LAMBDAS,bt_integral_result_);
}
