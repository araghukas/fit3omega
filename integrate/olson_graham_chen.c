#include "integrate.h"


// =================================================================================================
//
// IMPLEMENTS INTEGRAL AND INTEGRAND FROM OLSON, GRAHAM, AND CHEN Eq. (4) for complex Z
//
// [REF: Rev. Sci. Instrum. 76, 053901 (2005)]
// =================================================================================================


double complex Phi(int i_layer, double chi, double omega, double *kxs, double *kys, double *Cvs)
{
  /* OGC Eq. (6) */
  double b = HALF_WIDTH;
	return csqrt(kxs[i_layer]/kys[i_layer]*chi*chi + I*b*b*2.0*omega*Cvs[i_layer]/kys[i_layer]);
}


double complex fz(int i_layer, double chi, double omega,
	                double *ds, double *kxs, double *kys, double *Cvs, double *Rcs)
{
	/* OGC Eq. (5); rearrange for z_n = f(z_n+1) and number top-to-bottom */
	if (i_layer == N_LAYERS - 1)  // only semi-infinite B.C.
		return -HALF_WIDTH / (kys[i_layer] * Phi(i_layer,chi,omega,kxs,kys,Cvs));

	// recursive call until base case
	double complex z_ii = fz(i_layer+1,chi,omega,ds,kxs,kys,Cvs,Rcs);
	double complex Phi_ii = Phi(i_layer+1,chi,omega,kxs,kys,Cvs);	

	double complex C1 = kys[i_layer+1] * Phi_ii / HALF_WIDTH;
	double complex C2 = ctanh(Phi_ii * ds[i_layer+1] / HALF_WIDTH);
	return Rcs[i_layer] + (C1*z_ii + C2) / (C1 + C1*C1*C2*z_ii);
}


// =================================================================================================


double complex ogc_integrand(double chi, double omega,
	                           double *ds, double *kxs, double *kys, double *Cvs, double *Rcs)
{
	/* OGC Eq. (4) integrand */
	static const double A = 2.0 / M_PI; // 2x because integrand is symmetric in chi [-MAX,MAX]
	return A * (fz(0,chi,omega,ds,kxs,kys,Cvs,Rcs) + Rcs[0]) * sinc_sq(chi);
}


double complex *ogc_integral(double *ds, double *kxs, double *kys, double *Cvs, double *Rcs)
{
	/* OGC Eq. (4) integral */

	for (int i = 0; i < n_OMEGAS; i++) {

		ogc_integral_result_[i] = 0.0*I;
		double complex f0 = ogc_integrand(CHIS[0],OMEGAS[i],ds,kxs,kys,Cvs,Rcs);
		double complex f_prev = f0;

		for (int k = 1; k < N_XPTS; k++) {

			double complex fk = ogc_integrand(CHIS[k],OMEGAS[i],ds,kxs,kys,Cvs,Rcs);
			double complex dx = CHIS[k] - CHIS[k-1];
			ogc_integral_result_[i] += (dx / 2.0) * (fk + f_prev);
			f_prev = fk;

		}
	}
	return ogc_integral_result_;  // pointer to results for each ω
}