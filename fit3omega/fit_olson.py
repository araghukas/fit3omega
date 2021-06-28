import numpy as np

from fit3omega.general import BasinhoppingOptimizer
import fit3omega.plot as plot


class FitOlson(BasinhoppingOptimizer):
    FIT_PARAMS = ["kys", "ratio_xys", "Cvs", "Rcs"]

    def __init__(self, sample, data):
        super().__init__(sample, data)

        # clean up any heights targeted for fit (hence read into a str)
        heights = []
        for height in self.sample.heights:
            if type(height) is str:
                heights.append(float(height.rstrip('*')))
            else:
                heights.append(height)

        self._config_values = (
            heights,
            self.sample.kys,
            self.sample.ratio_xys,
            self.sample.Cvs,
            self.sample.Rcs
        )

    def Z2_func(self, kys, ratio_xys, Cvs, Rcs) -> np.ndarray:
        """
        Synthetic Z2 prediction from model and sample properties.
        This is the fit curve function.
        """
        if not self._integrators_ready or self._refresh_dependents:
            self._init_integrators()

        area = self.heater.width * self.heater.length
        return self.integrators.ogc_integral(self.config_values[0], kys, ratio_xys, Cvs, Rcs) / area

    def dZ2_func(self, kys, ratio_xys, Cvs, Rcs) -> np.ndarray:
        if not self._integrators_ready or self._refresh_dependents:
            self._init_integrators()

        area = self.heater.width * self.heater.length
        return self.integrators.ogc_jacobian(
            self.config_values[0], kys, ratio_xys, Cvs, Rcs) / area

    def T2_func(self, kys, ratio_xys, Cvs, Rcs):
        """synthetic average T2 at each ω, from sample parameters"""
        self._T2_func_latest = -self.power.norm * self.Z2_func(kys, ratio_xys, Cvs, Rcs)
        return self._T2_func_latest

    def dT2_func(self, kys, ratio_xys, Cvs, Rcs) -> np.ndarray:
        return -self.power.norm * self.dZ2_func(kys, ratio_xys, Cvs, Rcs)

    def error_func_jac(self, args) -> np.ndarray:
        """
        Jacobian matrix for the error function: sum( |T - T_measured|^2 / |T_measured|^2 )
        """
        args_T2 = self._sub_args_into_complete_params(args)
        jac_T2 = self.dT2_func(*args_T2)

        T2_func_ = self._T2_func_latest
        tx = (T2_func_.real - self.T2.x) / self.T2.norm_sq
        ty = (T2_func_.imag - self.T2.y) / self.T2.norm_sq
        error_func_jacobian = np.dot(jac_T2.real, tx) + np.dot(jac_T2.imag, ty)

        return (2. / len(T2_func_)) * error_func_jacobian

    def plot_fit(self, show=False):
        return plot.plot_fitted_Z2(self, show=show)

    def _init_integrators(self):
        self.integrators.ogc_set(self.omegas,
                                 self._ids,
                                 self.heater.width / 2.0,
                                 1e-6,
                                 15.,
                                 len(self.sample.layers))
        self._integrators_ready = True

    def get_current_T2(self):
        kys = self.sample.kys
        ratio_xys = self.sample.ratio_xys
        Cvs = self.sample.Cvs
        Rcs = self.sample.Rcs

        T2_complex = self.T2_func(kys, ratio_xys, Cvs, Rcs)
        T2_x = T2_complex.real
        T2_y = T2_complex.imag
        err = sum(((T2_x - self.T2.x)**2 + (T2_y - self.T2.y)**2) / self.T2.norm_sq)
        return T2_x, T2_y, err / (len(T2_complex))

    def _insert_extra_minimizer_kwargs(self, kwargs: dict) -> dict:
        kwargs['jac'] = self.error_func_jac
        return kwargs
