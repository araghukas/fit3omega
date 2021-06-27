#include "integrate.h"
#include "olson_graham_chen.h"
#include "exceptions.h"


// =================================================================================================
//
// IMPLEMENTS OPTIMIZATION ALGORITHM FROM OLSON, GRAHAM, AND CHEN Eqs. (10-23)
//
// [REF: Rev. Sci. Instrum. 76, 053901 (2005)]
// ================================================================================================


void fXis(double chi, double omega)
{
	/* OGC Eq. (11); for all layers at (χ,ω) */
	int i_layer = n_LAYERS - 1;
	Xis_[i_layer] = 0.0*I;

	i_layer--;
	while (i_layer >= 0) {
		double complex kPhi_b = kys_[i_layer] * Phis_[i_layer] / HALF_WIDTH;
		double complex kPhi_b_sq = kPhi_b * kPhi_b;
		double complex z_i = zs_[i_layer];
		double complex z_tilde = zs_[i_layer+1] - Rcs_[i_layer+1];
		Xis_[i_layer] = (1.0 - kPhi_b_sq * z_i * z_i) / (1.0 - kPhi_b_sq * z_tilde * z_tilde);
		i_layer--;
	}
}


// ================================================================================================


double complex *fdz0_dky(double c, double o)
{
	/* OGC Eq. (12); for all layers at (χ,ω) */
	int i_layer = n_LAYERS - 1;

	dz0_dky_[i_layer] = 1.0 / kys_[i_layer] + 0.0*I;

	for (int j = i_layer-1; j >= 0; j--)
		dz0_dky_[i_layer] *= Xis_[j];
	dz0_dky_[i_layer] *= -zs_[i_layer];  // z_tilde is 0 for substrate

	i_layer--;
	while (i_layer >= 0) {
		dz0_dky_[i_layer] = 1.0 / kys_[i_layer] + 0.0*I;
		for (int j = i_layer-1; j >= 0; j--)
			dz0_dky_[i_layer] *= Xis_[j];

		double complex z_tilde = zs_[i_layer+1] - Rcs_[i_layer+1];
		dz0_dky_[i_layer] *= (Xis_[i_layer]*z_tilde - zs_[i_layer]);
		i_layer--;
	}
	return dz0_dky_;
}


double complex *fdz0_dCv(double chi, double omega)
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
		dz0_dCv_[i_layer] *= Xis_[j];

	P = Phis_[i_layer];
	alphaPhi = ky * P / Cv;
	dz0_dCv_[i_layer] *= -I*omega*b*b / (alphaPhi*alphaPhi);

	z = zs_[i_layer];
	dz0_dCv_[i_layer] *= H/ky * (z*z*ky*ky*P*P/(b*b) - 1.0) - z;

	i_layer--;
	while (i_layer >= 0) {
		ky = kys_[i_layer];
		Cv = Cvs_[i_layer];
		dz0_dCv_[i_layer] = -ky / (Cv * Cv);

		for (int j = i_layer-1; j >= 0; j--)
			dz0_dCv_[i_layer] *= Xis_[j];

		P = Phis_[i_layer];
		alphaPhi = ky * P / Cv;
		dz0_dCv_[i_layer] *= -I*omega*b*b / (alphaPhi*alphaPhi);

		z = zs_[i_layer];
		double complex z_tilde = zs_[i_layer+1] - Rcs_[i_layer+1];
		dz0_dCv_[i_layer] *= H/ky * (z*z*ky*ky*P*P/(b*b) - 1.0) + Xis_[i_layer]*z_tilde - z;
		i_layer--;
	}
	return dz0_dCv_;
}



double complex *fdz0_dpsi(double chi, double omega)
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
		dz0_dpsi_[i_layer] *= Xis_[j];

	P = Phis_[i_layer];
	dz0_dpsi_[i_layer] *=  chi*chi/(2.0*P*P);

	z = zs_[i_layer];
	dz0_dpsi_[i_layer] *= H/ky * (z*z*ky*ky*P*P/(b*b) - 1.0) - z;

	i_layer--;
	while (i_layer >= 0) {
		ky = kys_[i_layer];
		Cv = Cvs_[i_layer];
		dz0_dpsi_[i_layer] = 1.0;

		for (int j = i_layer-1; j >= 0; j--)
			dz0_dpsi_[i_layer] *= Xis_[j];

		P = Phis_[i_layer];
		dz0_dpsi_[i_layer] *=  chi*chi/(2.0*P*P);

		z = zs_[i_layer];
		double complex z_tilde = zs_[i_layer+1] - Rcs_[i_layer+1];
		dz0_dpsi_[i_layer] *= H/ky * (z*z*ky*ky*P*P/(b*b) - 1.0) + Xis_[i_layer]*z_tilde - z;
		i_layer--;
	}
	return dz0_dpsi_;
}


double complex *fdz0_dRc(double c, double o)
{
	/* OGC Eq. (15); for all layers at (χ,ω) */
	int i_layer = n_LAYERS - 1;

	while (i_layer >= 0) {
		dz0_dRc_[i_layer] = -Xis_[i_layer];
		for (int j = i_layer-1; j >= 0; j--)
			dz0_dRc_[i_layer] *= Xis_[j];
		i_layer--;
	}
	return dz0_dRc_;
};


// ================================================================================================


double complex (*jac_Z(void))[MAX_n_OMEGAS]
{
	static const double A = 2.0 / M_PI;  // 2x because integrand is symmetric in chi [-MAX,MAX]

	for (int i = 0; i < n_OMEGAS; i++) {
		double omega = OMEGAS[i];

		for (int k = 0; k < N_XPTS; k++) {
			double chi = CHIS[k];
			double sinq_sq_ = sinc_sq(chi);
			fPhis(chi,omega);
			fzs(chi,omega);
			fXis(chi,omega);

			for (int n = 0; n < n_PARAMS; n++) {
				int i_param = param_ids_[n][0];
				int i_layer = param_ids_[n][1];

				double complex *f_values = NULL;
				switch(i_param)
				{
					case 0:
						f_values = fdz0_dky(chi,omega);
						break;
					case 1:
						f_values = fdz0_dpsi(chi,omega);
						break;
					case 2:
						f_values = fdz0_dCv(chi,omega);
						break;
					case 3:
						f_values = fdz0_dRc(chi,omega);
						break;
					default:
						PyErr_SetString(ParameterIDError, ParameterIDError_MSG);
						return jac_Z_result_;
				}
				jac_Z_fs_buff_[n][k] = A * sinq_sq_ * f_values[i_layer];
			}
		}

		for (int n = 0; n < n_PARAMS; n++)
			jac_Z_result_[n][i] = val_trapz(jac_Z_fs_buff_[n],CHIS);
	}

	return jac_Z_result_;
}
