#include "integrate.h"


// =================================================================================================
//
// IMPLEMENTS INTEGRAL AND INTEGRAND FROM OLSON, GRAHAM, AND CHEN Eq. (4) for complex Z
//
// [REF: Rev. Sci. Instrum. 76, 053901 (2005)]
// =================================================================================================


double complex *fPhi(int i_layer, double chi, double *kxs, double *kys, double *Cvs)
{
  /* OGC Eq. (6) */
  for (int i = 0; i < n_OMEGAS; i++) {
    double re = kxs[i_layer] / kys[i_layer] * chi * chi;
    double im = 2.0 * OMEGAS[i] * b_ * b_ * Cvs[i_layer] / kys[i_layer];
    fPhi_[i] = csqrt(re + I*im);
  }
  return fPhi_;
}


double complex *fz(int i_layer, double chi, double *kxs, double *kys, double *Cvs, double *Rcs)
{

}