#include "integrate.h"
#include "olson_graham_chen.h"


// =================================================================================================
//
// IMPLEMENTS INTEGRAL AND INTEGRAND FROM OLSON, GRAHAM, AND CHEN Eq. (4) for complex Z
//
// [REF: Rev. Sci. Instrum. 76, 053901 (2005)]
// =================================================================================================


double complex Phi(int i_layer, double chi, double omega)
{
    /* OGC Eq. (6) */
    double b = HALF_WIDTH;
    return csqrt(psis_[i_layer]*chi*chi + I*b*b*2.0*omega*Cvs_[i_layer]/kys_[i_layer]);
}


void fPhis(double chi, double omega)
{
	/* OGC Eq. (6);  for all layers at (χ,ω) */
	int i_layer = n_LAYERS - 1;
	Phis_[i_layer] = Phi(i_layer,chi,omega);
	i_layer--;
	while (i_layer >= 0) {
		Phis_[i_layer] = Phi(i_layer,chi,omega);
		i_layer--;
	}
}


void fzs(double chi, double omega)
{
	/* OGC Eq. (5); (the z w/ no tilde) for all layers at (χ,ω)  */
	int i_layer = n_LAYERS - 1;
	zs_[i_layer] = -HALF_WIDTH / (kys_[i_layer] * Phis_[i_layer]);
	i_layer--;
	while (i_layer >= 0) {
		double complex P = Phis_[i_layer];
		double complex kPhi_b = kys_[i_layer] * P / HALF_WIDTH;
		double complex tanh_term = ctanh(P * ds_[i_layer] / HALF_WIDTH);
		double complex z_tilde = zs_[i_layer+1] - Rcs_[i_layer+1];
		zs_[i_layer] = (kPhi_b * z_tilde - tanh_term)
		                / (kPhi_b - kPhi_b * kPhi_b * z_tilde * tanh_term);
		i_layer--;
	}
}


// =================================================================================================


double complex ogc_integrand(double chi, double omega)
{
	/* OGC Eq. (4) integrand */
	static const double A = 2.0 / M_PI; // 2x because integrand is symmetric in chi [-MAX,MAX]

	fPhis(chi,omega);
	fzs(chi,omega);

	return A * (zs_[0] - Rcs_[0]) * sinc_sq(chi);
}


double complex *ogc_integral(void)
{
	/* OGC Eq. (4) integral */
	return trapz(ogc_integrand,CHIS,ogc_integral_result_);
}
