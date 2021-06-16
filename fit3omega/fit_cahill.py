import numpy as np

from fit3omega.general import BasinhoppingOptimizer
from fit3omega.model import ROOT2
import fit3omega.plot as plot


class FitCahill(BasinhoppingOptimizer):

    BOUNDARY_TYPES = ['s', 'i', 'a']
    FIT_ARG_NAMES = ["heights", "kys", "ratio_xys", "Cvs"]

    def __init__(self, sample, data, b_type):
        super().__init__(sample, data)

        self._defaults = (
            self.sample.heights, self.sample.kys, self.sample.ratio_xys, self.sample.Cvs
        )

        for lst in self._defaults:
            self._full_args += lst

        if b_type not in FitCahill.BOUNDARY_TYPES:
            raise ValueError("boundary type {} is not one of {}"
                             .format(b_type, FitCahill.BOUNDARY_TYPES))
        self.b_type = b_type

    @property
    def defaults(self):
        return self._defaults

    def T2_func(self, heights, kys, ratio_xys, Cvs) -> np.ndarray:
        """T2 prediction from physical model and provided properties"""
        if not self._integrators_ready or self._refresh_dependents:
            self._init_integrators()

        # extra divisor of ROOT2 since measured T2 amplitudes are RMS
        P = -1. / (np.pi * self.heater.length * kys[0] * ROOT2) * self.power.x
        return P * self.integrators.bt_integral(heights, kys, ratio_xys, Cvs)

    def plot_fit(self, show=False):
        """plot fit result"""
        return plot.plot_fitted_T2(self, show=show)

    def error_func(self, args) -> float:
        """objective function for the fit method"""
        err_args = self._full_args.copy()
        for j, index in enumerate(self._fit_indices):
            err_args[index] = args[j]

        # reconstruct T2_func args
        args_T2 = tuple()
        for k in range(len(self._full_args) // self.n_layers):
            i_min = k * self.n_layers
            i_max = i_min + self.n_layers
            args_T2 += (err_args[i_min:i_max],)

        T2_func_ = self.T2_func(*args_T2)

        err = sum(np.abs(self.T2.x - T2_func_.real))
        err += sum(np.abs(self.T2.y - T2_func_.imag))
        return err / len(T2_func_)

    def _init_integrators(self):
        self.integrators.bt_set(self.omegas,
                                self.heater.width / 2,
                                1e-3,
                                1e7,
                                len(self.sample.layers),
                                self.b_type.encode('utf-8'))
        self._integrators_ready = True

    def _get_initial_values(self):
        return [
            ("heights", self.sample.heights.copy()),
            ("kys", self.sample.kys.copy()),
            ("ratio_xys", self.sample.ratio_xys.copy()),
            ("Cvs", self.sample.Cvs.copy())
        ]