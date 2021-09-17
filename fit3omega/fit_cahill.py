import numpy as np

from fit3omega.base_optimizer import BasinhoppingOptimizer
import fit3omega.plot as plot


class FitCahill(BasinhoppingOptimizer):
    BOUNDARY_TYPES = ['s', 'i', 'a']
    FIT_PARAMS = ["heights", "kys", "ratio_xys", "Cvs"]

    def __init__(self, sample, data, b_type):
        super().__init__(sample, data)

        if b_type not in FitCahill.BOUNDARY_TYPES:
            raise ValueError("boundary type {} is not one of {}"
                             .format(b_type, FitCahill.BOUNDARY_TYPES))
        self.b_type = b_type

    def T2_func(self, heights, kys, ratio_xys, Cvs) -> np.ndarray:
        """
        Synthetic average T2 at each ω, from sample parameters.
        This is the fit curve function.
        """
        if not self._integrators_ready or self._refresh_dependents:
            self._init_integrators()

        # Using average power to get average temperature
        p_ave = -self.power.x / (np.pi * self.heater.length * kys[0])
        return p_ave * self.integrators.bt_integral(heights, kys, ratio_xys, Cvs)

    def plot_fit(self, show=False):
        """plot fit result"""
        return plot.plot_fitted_T2(self, show=show)

    def _init_integrators(self):
        self.integrators.bt_set(self.omegas,
                                self.heater.width / 2,
                                1e-3,
                                1e7,
                                len(self.sample.layers),
                                self.b_type.encode('utf-8'))
        self._integrators_ready = True

    def get_current_T2(self):
        heights = self.sample.heights
        kys = self.sample.kys
        ratio_xys = self.sample.ratio_xys
        Cvs = self.sample.Cvs

        T2_complex = self.T2_func(heights, kys, ratio_xys, Cvs)
        T2_x = T2_complex.real
        T2_y = T2_complex.imag
        err = sum(((T2_x - self.T2.x)**2 + (T2_y - self.T2.y)**2) / self.T2.norm_sq)
        return T2_x, T2_y, err / (len(T2_complex))

    def _insert_extra_minimizer_kwargs(self, kwargs: dict) -> dict:
        return kwargs
