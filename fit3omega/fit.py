"""
A class for fitting the measured data with a given the sample configuration.
"""
from dataclasses import dataclass
from typing import Union, List
from scipy.optimize import minimize, OptimizeResult
import numpy as np

from fit3omega.model import Model
from fit3omega.sample import SampleParameters, load_sample_parameters
from fit3omega.data import Data, ACReading
import fit3omega.utils as utils


class Fit3omega(Model):
    """fits sample parameters to the measured voltage data"""

    @property
    def T2(self) -> ACReading:
        """
        Average 2ω temperature oscillations at each ω
        CORRECTED for heater properties using Borca-Tasciuc (Ref. 1) Eq. (20)

        NOTE: The correction has no effect if Rth and Cv or d are 0.
        """
        if self._T2 is None or self._refresh_dependents:
            T2_measured_raw = super().T2
            cT2_raw = T2_measured_raw.as_complex()
            Rth = self.sample.heater.Rc
            Cv = self.sample.heater.Cv
            d = self.sample.heater.height
            area = self._heater_area
            T2 = ((cT2_raw + Rth * self.power.norm / area)
                  / (1.0 + 2.0j * self.data.omegas * Cv * d * (
                            Rth + cT2_raw * area / self.power.norm)))
            # NOTE: error is not scaled
            self._T2 = ACReading(T2.real, T2.imag, T2_measured_raw.xerr, T2_measured_raw.yerr)

        return self._T2

    @property
    def result(self) -> 'FitResult':
        """result of the latest fit"""
        return self._result

    @property
    def fitted_T2(self) -> np.ndarray:
        """complex array for fitted temperature rise at each frequency"""
        return self.T2_function(*self.sample.substitute(self.sample.x))

    def __init__(self,
                 sample: Union[str, SampleParameters],
                 data: Union[str, Data]):

        if type(sample) is str:
            sample = load_sample_parameters(sample)
        if type(data) is str:
            data = Data(data)
        super().__init__(sample, data)

        # store initial configuration to track changes
        self._previous_sample = sample.copy()
        self._original_sample = sample.copy()

        # C-extension for computing integrals
        self._integrator_module = __import__('integrate')
        self._integrators_ready = False

        # some constants
        self._layer_heights = [layer.height for layer in self.sample.layers]
        self._heater_area = self.sample.heater.width * self.sample.heater.length
        self._n_omegas = len(self.data.omegas)

        # fit result (assigned by calling the `fit` method)
        self._result = None

    def fit(self, tol: float = 1e-6, x0: np.ndarray = None) -> None:
        """
        Run the fitting algorithm to estimate parameters.

        :param tol: termination tolerance
        :param x0: initial fit arguments vector
        """
        if x0 is None:
            x0 = self.sample.x

        result = minimize(fun=self.objective_func,
                          x0=x0,
                          args=None,
                          method='TNC',
                          tol=tol,
                          bounds=utils.positive_bounds(x0, min_frac=1e-6, max_frac=1e3),
                          options={'disp': False, 'maxiter': 200, 'stepmx': 100})

        self._record_result(result)

    def objective_func(self, *args) -> float:
        """returns the value of the objective function (MSE)."""
        args_T2 = self.sample.substitute(args[0])
        T2_func_values = self.T2_function(*args_T2)
        dx = T2_func_values.real - self.T2.x
        dy = T2_func_values.imag - self.T2.y
        return sum((dx**2 + dy**2) / self.T2.norm_sq) / self._n_omegas

    def T2_function(self,
                    kys: List[float],
                    ratio_xys: List[float],
                    Cvs: List[float],
                    Rcs: List[float]) -> np.ndarray:
        """
        Computes a prediction of the 2ω temperature rise based on
        arbitrary layer parameters.
        """
        if not self._integrators_ready:
            self._init_integrators()

        integral = self._integrator_module.ogc_integral(
            self._layer_heights,
            kys,
            ratio_xys,
            Cvs,
            Rcs
        )
        return -self.power.norm / self._heater_area * integral

    def _init_integrators(self) -> None:
        """initialize the integrator module"""
        self._integrator_module.ogc_set(self.data.omegas,
                                        self.sample.fit_indices,
                                        self.sample.heater.width / 2.0,
                                        1e-6,
                                        15.,
                                        len(self.sample.layers))
        self._integrators_ready = True

    def _record_result(self, result: OptimizeResult):
        self._result = FitResult(result, self._previous_sample)
        self._previous_sample = self.sample.copy()


@dataclass(frozen=True)
class FitResult:
    """a container for results of a data fit"""
    result: OptimizeResult
    previous_sample: SampleParameters

    @property
    def x(self) -> np.ndarray:
        """fitted argument vector"""
        return self.result.x

    @property
    def error(self) -> float:
        """residual value of the objective function"""
        return self.result.fun

    @property
    def summary(self) -> str:
        """return a string summarized the result"""
        x0 = self.previous_sample.x
        diff_percents = [1e2 * (xf - xi) / xi for xi, xf in zip(x0, self.x)]
        diff_signs = ['+' if p >= 0 else '-' for p in diff_percents]

        props_by_layer = {}
        for i, idx in enumerate(self.previous_sample.fit_indices):
            layer_name = self.previous_sample.layers[idx[1]].name
            param_name = self.previous_sample.FIELDS[idx[0]].rstrip('s')
            change_str = (
                "{:>8} --> {:.2e} ({} %)".format(param_name,
                                                 self.x[i],
                                                 diff_signs[i] + "%.2f" % abs(diff_percents[i]))
            )
            if layer_name in props_by_layer:
                props_by_layer[layer_name].append(change_str)
            else:
                props_by_layer[layer_name] = [change_str]

        lines = []
        for layer_name, change_str_list in props_by_layer.items():
            lines.append(layer_name + ":")
            for change_str in change_str_list:
                lines.append("    " + change_str)

        return "\n".join(lines)

    def __repr__(self):
        return self.summary
