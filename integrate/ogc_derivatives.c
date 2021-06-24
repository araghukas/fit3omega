#include "integrate.h"
#include "olson_graham_chen.h"
#include "exceptions.h"


// =================================================================================================
//
// IMPLEMENTS OPTIMIZATION ALGORITHM FROM OLSON, GRAHAM, AND CHEN Eqs. (10-23)
//
// [REF: Rev. Sci. Instrum. 76, 053901 (2005)]
// ================================================================================================


double complex *fXis(double chi, double omega, double complex *Phis, double complex *zs)
{
	/* OGC Eq. (11); chain rule terms for all layers at (χ,ω) */
	int i_layer = n_LAYERS - 1;
	Xis_[i_layer] = 0.0*I;

	i_layer--;
	while (i_layer >= 0) {
		double complex kPhi_b = kys_[i_layer] * Phis[i_layer] / HALF_WIDTH;
		double complex kPhi_b_sq = kPhi_b * kPhi_b;
		double complex z_i = zs[i_layer];
		double complex z_tilde = zs[i_layer+1] - Rcs_[i_layer+1];
		Xis_[i_layer] = (1.0 - kPhi_b_sq * z_i * z_i)
						 /(1.0 - kPhi_b_sq * z_tilde * z_tilde);
		i_layer--;
	}

	return Xis_;
}


// ================================================================================================


double complex *fdz0_dky(double complex *zs, double complex *Xis)
{
	/* OGC Eq. (12); for all layers at (χ,ω) */
	int i_layer = n_LAYERS - 1;

	dz0_dky_[i_layer] = 1.0 / kys_[i_layer] + 0.0*I;

	for (int j = i_layer-1; j >= 0; j--)
		dz0_dky_[i_layer] *= Xis[j];
	dz0_dky_[i_layer] *= -zs[i_layer];  // z_tilde is 0 for substrate

	i_layer--;
	while (i_layer >= 0) {
		dz0_dky_[i_layer] = 1.0 / kys_[i_layer] + 0.0*I;
		for (int j = i_layer-1; j >= 0; j--)
			dz0_dky_[i_layer] *= Xis[j];

		double complex z_tilde = zs[i_layer+1] - Rcs_[i_layer+1];
		dz0_dky_[i_layer] *= (Xis[i_layer]*z_tilde - zs[i_layer]);
		i_layer--;
	}

	return dz0_dky_;
}


double complex *fdz0_dCv(double chi, double omega, double complex *Phis, double complex *zs, double complex *Xis)
{
	/*
	OGC Eq. (13); for all layers at (χ,ω)

	Note: dz/dCv = dz/da * da/dCv = (-ky / Cv^2) * dz/da
	*/

	int i_layer = n_LAYERS - 1;

	const double b = HALF_WIDTH;
	const double H = ds_[n_LAYERS-1];
	double ky;
	double Cv;
	double complex z;
	double complex P;
	double complex alphaPhi;

	ky = kys_[i_layer];
	Cv = Cvs_[i_layer];
	dz0_dCv_[i_layer] = -ky / (Cv * Cv);  // chain rule factor

	for (int j = i_layer-1; j >= 0; j--)
		dz0_dCv_[i_layer] *= Xis[j];

	P = Phis[i_layer];
	alphaPhi = ky * P / Cv;
	dz0_dCv_[i_layer] *= -I*omega*b*b / (alphaPhi*alphaPhi);

	z = zs[i_layer];
	dz0_dCv_[i_layer] *= H/ky * (z*z*ky*ky*P*P/(b*b) - 1.0) - z;

	i_layer--;
	while (i_layer >= 0) {
		ky = kys_[i_layer];
		Cv = Cvs_[i_layer];
		dz0_dCv_[i_layer] = -ky / (Cv * Cv);

		for (int j = i_layer-1; j >= 0; j--)
			dz0_dCv_[i_layer] *= Xis[j];

		P = Phis[i_layer];
		alphaPhi = ky * P / Cv;
		dz0_dCv_[i_layer] *= -I*omega*b*b / (alphaPhi*alphaPhi);

		z = zs[i_layer];
		double complex z_tilde = zs[i_layer+1] - Rcs_[i_layer+1];
		dz0_dCv_[i_layer] *= H/ky * (z*z*ky*ky*P*P/(b*b) - 1.0) + Xis[i_layer]*z_tilde - z;
		i_layer--;
	}

	return dz0_dCv_;
}



double complex *fdz0_dpsi_(double chi, double omega, double complex *Phis, double complex *zs, double complex *Xis)
{
	/* OGC Eq. (14); for all layers at (χ,ω) */

	int i_layer = n_LAYERS - 1;

	const double b = HALF_WIDTH;
	const double H = ds_[n_LAYERS-1];
	double ky;
	double Cv;
	double complex z;
	double complex P;

	ky = kys_[i_layer];
	Cv = Cvs_[i_layer];
	dz0_dpsi_[i_layer] = 1.0;  // chain rule factor

	for (int j = i_layer-1; j >= 0; j--)
		dz0_dpsi_[i_layer] *= Xis[j];

	P = Phis[i_layer];
	dz0_dpsi_[i_layer] *=  chi*chi/(2.0*P*P);

	z = zs[i_layer];
	dz0_dpsi_[i_layer] *= H/ky * (z*z*ky*ky*P*P/(b*b) - 1.0) - z;

	i_layer--;
	while (i_layer >= 0) {
		ky = kys_[i_layer];
		Cv = Cvs_[i_layer];
		dz0_dpsi_[i_layer] = 1.0;

		for (int j = i_layer-1; j >= 0; j--)
			dz0_dpsi_[i_layer] *= Xis[j];

		P = Phis[i_layer];
		dz0_dpsi_[i_layer] *=  chi*chi/(2.0*P*P);

		z = zs[i_layer];
		double complex z_tilde = zs[i_layer+1] - Rcs_[i_layer+1];
		dz0_dpsi_[i_layer] *= H/ky * (z*z*ky*ky*P*P/(b*b) - 1.0) + Xis[i_layer]*z_tilde - z;
		i_layer--;
	}

	return dz0_dpsi_;
}


double complex *fdz0_dRc(double complex *Xis)
{
	/* OGC Eq. (15); for all layers at (χ,ω) */
	int i_layer = n_LAYERS - 1;

	while (i_layer >= 0) {
		dz0_dRc_[i_layer] = -Xis[i_layer];
		for (int j = i_layer-1; j >= 0; j--)
			dz0_dRc_[i_layer] *= Xis[j];
		i_layer--;
	}

	return dz0_dRc_;
};


// ================================================================================================

double complex **jac_Z(void)
{

	for (int i = 0; i < n_OMEGAS; i++) {
		for (int n = 0; n < n_PARAMS; n++) {

			double complex *row = jac_J_result_[n];
			
			for (int k = 0; k < N_XPTS; k++) {
				double chi = CHIS[k];
				double omega = OMEGAS[i];
				Phis_ = fPhis(chi,omega);
				zs_ = fzs(chi,omega);
				Xis_ = fXis(chi,omega,Phis_,zs_);
			}

		}
	}

	return jac_Z_result_;
}


