"""objects that perform fits on measured data are defined here"""
from typing import List, Union, Tuple
import numpy as np
from scipy.optimize import minimize, OptimizeResult
from scipy.interpolate import interp1d

from fit3omega.base import BaseFitter, _set_mpl_defaults
from fit3omega.sample import Sample
from fit3omega.varsample import VarSample
from fit3omega.model import Model
from fit3omega.data import ACReading, Data
from fit3omega.result import FitResult
from fit3omega.utils import positive_bounds
from fit3omega.params import ParameterMap


class FitLinear(BaseFitter):
    """
    Determines film and substrate thermal conductivity by fitting
    the 3-omega voltage in a linear region.
    """
    eta = 3 / 2 - np.euler_gamma

    def __init__(self,
                 sample: Union[str, Sample],
                 data: Union[str, Data]):
        model = Model(sample, data)
        super().__init__(model)

        # default values for fitline options
        self._thresh = 0.99
        self._min_fitline_length = 5

        # linear subset of the total data
        self._omegas_linear = None
        self._T2_linear = None
        self._coeffs_re, self._Rsq_re = None, None
        self._coeffs_im, self._Rsq_im = None, None
        self._f_linear = None

    @property
    def T2_fit(self) -> ACReading:
        if self._result is None:
            raise ValueError("fit line has not been computed.")
        x = np.log(self._omegas_linear)
        m_re, b_re = self._coeffs_re[self._i_linear]
        m_im, b_im = self._coeffs_im[self._i_linear]
        return ACReading(m_re * x + b_re,
                         m_im * x + b_im,
                         0.0,  # TODO: error estimate
                         0.0)

    @property
    def result(self) -> FitResult:
        return self._result

    @property
    def error(self) -> float:
        return self._error

    def T2_func(self, omega) -> np.ndarray:
        return self._f_linear(omega)

    def fit(self, *args, **kwargs) -> FitResult:
        self._calculate_linearity()

        p = self.model.power.x[:self._i_linear]
        l = self.model.heater.length
        b = self.model.heater.width / 2.
        a1, a0 = self._coeffs_re[self._i_linear]
        eta = FitLinear.eta

        fitted_kwargs = {}

        if len(self.model.sample.layers) == 1:
            # heater-on-substrate sample
            subs = self.model.sample.layers[0]
            rS = subs.ratio_xy
            kSy = -p / (2. * np.pi * l * a1 * np.sqrt(rS))  # substrate conductivity
            fitted_kwargs['kys'] = [np.average(kSy)]
            if type(subs.Cv) is str:
                Cv = kSy / b**2 * np.exp((2. * p * eta - 2 * np.pi * l * kSy * a0) / p)
                Cv = np.average(Cv)
                fitted_kwargs['Cvs'] = [Cv]
        else:
            # heater-on-film-on-substrate sample
            film, subs = self.model.sample.layers[:2]
            rS = subs.ratio_xy
            root_rS = np.sqrt(rS)

            fitted_substrate = False
            if type(subs.ky) is str:
                kSy = -p / (2. * np.pi * l * a1 * root_rS)
                fitted_substrate = True
            else:
                kSy = subs.ky

            CvS = subs.Cv
            TS = p / (np.pi * l * kSy * root_rS) * (
                    .5 * np.log((kSy * root_rS) / (CvS * b**2))
                    - .5 * np.log(self._omegas_linear)
                    + eta
            )

            dF = film.height
            kF = np.average(p * dF / (2. * b * l) / (self._T2_linear.x - TS))
            kSy = np.average(kSy)
            fitted_kwargs['kys'] = [kF, np.kSy] if fitted_substrate else [None, kSy]

        self._result = FitResult(self.model.sample, fitted_kwargs)
        self._error = sum((self.T2_fit.x - self._T2_linear.x)**2) / len(self.T2_fit.x)
        return self._result

    @property
    def Rsq(self) -> float:
        """R-squared value of the in-phase voltage fit line"""
        return self._Rsq_re

    @staticmethod
    def Rsq_func(line_y, points_y):
        """data points Y; fit values y"""
        mean_points_y = np.average(points_y)
        SS_res = sum((line_y - points_y)**2)
        SS_tot = sum((points_y - mean_points_y)**2)
        if SS_tot == 0.0 and SS_res < np.finfo(float).eps:
            return 0.
        return 1. - SS_res / SS_tot

    def _calculate_linearity(self, thresh=None, min_fitline_length=None):
        if thresh is None:
            thresh = self._thresh
        if min_fitline_length is None:
            min_fitline_length = self._min_fitline_length

        coeffs_re, Rsq_re = [None, ], [None, ]
        coeffs_im, Rsq_im = [None, ], [None, ]

        N = len(self.model.omegas)
        x = np.log(self.model.omegas)
        for i in range(N - 1):
            i_ = i + 2
            xi = x[:i_]
            yi_re = self.model.T2.x[:i_]
            yi_im = self.model.T2.y[:i_]
            mb_re = np.polyfit(xi, yi_re, 1)
            mb_im = np.polyfit(xi, yi_im, 1)
            y_line_re = mb_re[0] * xi + mb_re[1]
            y_line_im = mb_im[0] * xi + mb_im[1]
            Rsq_re.append(FitLinear.Rsq_func(y_line_re, yi_re))
            Rsq_im.append(FitLinear.Rsq_func(y_line_im, yi_im))
            coeffs_re.append(mb_re)
            coeffs_im.append(mb_im)

        k = min_fitline_length - 1
        while k < N and Rsq_re[k] > thresh:
            k += 1

        self._i_linear = k - 1
        self._coeffs_re, self._Rsq_re = coeffs_re, Rsq_re
        self._coeffs_im, self._Rsq_im = coeffs_im, Rsq_im

        self._omegas_linear = self.model.omegas[:self._i_linear]
        self._T2_linear = ACReading(self.model.T2.x[:self._i_linear],
                                    self.model.T2.y[:self._i_linear],
                                    self.model.T2.xerr[:self._i_linear],
                                    self.model.T2.yerr[:self._i_linear])

        m_re, b_re = self._coeffs_re[self._i_linear]
        m_im, b_im = self._coeffs_im[self._i_linear]
        x = np.log(self._omegas_linear)
        y = (m_re * x + b_re) + 1.j * (m_im * x + b_im)
        self._f_linear = interp1d(x, y)


class FitGeneral(BaseFitter):
    """
    Determines film and substrate thermal conductivity by fitting
    the entire 3-omega voltage curve.
    """

    def __init__(self,
                 sample: Union[str, Sample],
                 data: Union[str, Data]):
        model = Model(sample, data)
        super().__init__(model)

        self._n_omegas = len(model.data.omegas)

        # for fitting of an arbitrary parameter set
        self._map = ParameterMap(model.sample)

        # C-extension integrator module
        self._integrator_module = __import__('integrate')
        self._integrators_ready = False

        # convenience attributes
        self._heights = []
        for height in self.model.sample.heights:  # can't fit height, clean up any targets
            if type(height) is str:
                self._heights.append(float(height.rstrip('*')))
            else:
                self._heights.append(height)

        # set after the fit method is called
        self._T2_fit = None

    @property
    def T2_fit(self) -> ACReading:
        """synthetic values along the optimized fit line"""
        return self._T2_fit

    @property
    def result(self) -> FitResult:
        return self._result

    @property
    def error(self) -> float:
        return self._error

    def fit(self, tol: float = 1e-6) -> FitResult:
        """
        Run the fitting algorithm to estimate parameters.

        :param tol: termination tolerance
        """
        x0 = self._map.initial_values
        result = minimize(fun=self.objective_func,
                          x0=x0,
                          args=None,
                          method='TNC',
                          tol=tol,
                          bounds=positive_bounds(x0, min_frac=1e-6, max_frac=1e3),
                          options={'disp': True, 'maxiter': 200, 'stepmx': 100})

        self._record_result(result)
        return self._result

    def objective_func(self, *args) -> float:
        """returns the value of the objective function (MSE)."""
        args_T2 = self._map.substitute(args[0])
        T2_func_values = self.T2_func(*args_T2)
        dx = T2_func_values.real - self.model.T2.x
        dy = T2_func_values.imag - self.model.T2.y
        return sum((dx**2 + dy**2) / self.model.T2.norm_sq) / self._n_omegas

    def Z2_func(self,
                kys: List[float],
                ratio_xys: List[float],
                Cvs: List[float],
                Rcs: List[float]) -> np.ndarray:
        """
        Synthetic Z2 prediction from model and sample properties.
        This is the fit curve function.
        """
        if not self._integrators_ready or self.model.refresh:
            self._init_integrators()

        area = self.model.heater.width * self.model.heater.length
        return self._integrator_module.ogc_integral(
            self._heights, kys, ratio_xys, Cvs, Rcs) / area

    def T2_func(self,
                kys: List[float],
                ratio_xys: List[float],
                Cvs: List[float],
                Rcs: List[float]) -> np.ndarray:
        return -self.model.power.norm * self.Z2_func(kys, ratio_xys, Cvs, Rcs)

    def get_T2_kwargs(self):
        return self._map.get_complete_kwargs(self._map.initial_values)

    def _init_integrators(self) -> None:
        """initialize the integrator module"""
        self._integrator_module.ogc_set(self.model.omegas,
                                        self._map.integrator_ids,
                                        self.model.heater.width / 2.0,
                                        1e-6,
                                        15.,
                                        len(self.model.sample.layers))
        self._integrators_ready = True

    def _record_result(self, fit_result: OptimizeResult) -> None:
        """set post-fit parameters"""

        fitted_kwargs = self._map.get_fitted_kwargs(fit_result.x)
        self._result = FitResult(self.model.sample, fitted_kwargs)
        self._error = fit_result.fun

        T2_kwargs = self._map.get_complete_kwargs(fit_result.x)
        T2_vals = self.T2_func(**T2_kwargs)
        self._T2_fit = ACReading(T2_vals.real, T2_vals.imag, 0.0, 0.0)
