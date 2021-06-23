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

    def Z2_func_jac(self, ids) -> np.ndarray:
        """
        Jacobian matrix for the above function
        """
        return np.zeros(1)

    def T2_func(self, kys, ratio_xys, Cvs, Rcs):
        """synthetic average T2 at each ω, from sample parameters"""
        return -self.power.norm * self.Z2_func(kys, ratio_xys, Cvs, Rcs)

    def plot_fit(self, show=False):
        return plot.plot_fitted_Z2(self, show=show)

    def _init_integrators(self):
        self.integrators.ogc_set(self.omegas,
                                 self.heater.width / 2.0,
                                 1e-6,
                                 15.,
                                 len(self.sample.layers))
        self._integrators_ready = True

    def get_current_T2(self):
        heights = self.sample.heights
        kys = self.sample.kys
        ratio_xys = self.sample.ratio_xys
        Cvs = self.sample.Cvs
        Rcs = self.sample.Rcs

        T2_complex = self.T2_func(heights, kys, ratio_xys, Cvs, Rcs)
        T2_x = T2_complex.real
        T2_y = T2_complex.imag
        err = (sum(abs(self.T2.x - T2_complex.real))
               + sum(abs(self.T2.y - T2_complex.imag)))
        return T2_x, T2_y, err / (2 * len(T2_x))
