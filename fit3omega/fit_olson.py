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
        return self.power.norm * self.Z2_func(heights, kys, ratio_xys, Cvs, Rcs)

    def plot_fit(self, show=False):
        return plot.plot_fitted_Z2(self, show=show)

    def error_func(self, args):
        """objective function for the fit method"""
        err_args = self._full_args.copy()
        for j, index in enumerate(self._fit_indices):
            err_args[index] = args[j]

        # reconstruct Z2_func args
        args_Z2 = tuple()
        for k in range(len(self._full_args) // self.n_layers):
            i_min = k * self.n_layers
            i_max = i_min + self.n_layers
            args_Z2 += (err_args[i_min:i_max],)

        Z2_func_ = self.Z2_func(*args_Z2)

        err = sum((self.Z2.x - Z2_func_.real)**2)
        err += sum((self.Z2.y - Z2_func_.imag)**2)
        return err / len(Z2_func_)

    def _init_integrators(self):
        self.integrators.ogc_set(self.omegas,
                                 self.heater.width / 2.0,
                                 0.,
                                 15.,
                                 len(self.sample.layers))
        self._integrators_ready = True
