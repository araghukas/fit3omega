import numpy as np

from fit3omega.general import BasinhoppingOptimizer
import fit3omega.plot as plot


class FitOlson(BasinhoppingOptimizer):
    FIT_ARG_NAMES = ["heights", "kys", "ratio_xys", "Cvs", "Rcs"]

    def __init__(self, sample, data):
        super().__init__(sample, data)

        self._defaults = (
            self.sample.heights,
            self.sample.kys,
            self.sample.ratio_xys,
            self.sample.Cvs,
            self.sample.Rcs
        )

    def Z2_func(self, heights, kys, ratio_xys, Cvs, Rcs) -> np.ndarray:
        """
        Synthetic Z2 prediction from model and sample properties.
        This is the fit curve function.
        """
        if not self._integrators_ready or self._refresh_dependents:
            self._init_integrators()

        area = self.heater.width * self.heater.length
        return self.integrators.ogc_integral(heights, kys, ratio_xys, Cvs, Rcs) / area

    def T2_func(self, heights, kys, ratio_xys, Cvs, Rcs):
        """synthetic average T2 at each ω, from sample parameters"""
        return -self.power.norm * self.Z2_func(heights, kys, ratio_xys, Cvs, Rcs)

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
