import numpy as np
from scipy.optimize import basinhopping

import fit3omega.plot as plot
from fit3omega.model import Model
from fit3omega.data import ACReading

ROOT2 = np.sqrt(2)


class Stepper:
    def __init__(self, step_sizes):
        self.step_sizes = step_sizes
        self._bounds_func = None

    def __call__(self, x):
        for i, s in enumerate(self.step_sizes):
            x[i] += np.random.uniform(-s, s)
        return x

    def get_bounds(self, *args, **kwargs):
        return self._bounds_func(*args, **kwargs)

    @classmethod
    def by_fraction(cls, guesses: list, frac: float):
        step_sizes = [guess * frac for guess in guesses]
        stp = cls(step_sizes)
        stp._bounds_func = Stepper.bounds_by_fraction
        return stp

    @staticmethod
    def bounds_by_fraction(guesses: list, frac: float):
        f = frac / 2
        a1 = 1 - f
        a2 = 1 + f
        if a1 > a2:
            a1, a2 = a2, a1

        a2 = 0.0 if frac < 0 < a2 else a2  # neg. fraction clamp to [1+f  ,   0]
        a1 = 0.0 if a1 < 0 < frac else a1  # pos. fraction clamp to [    0, 1+f]

        return [(a1 * guess, a2 * guess) for guess in guesses]


class BasicPrinterCallBack:
    def __init__(self, i_max):
        self.counter = 1

        if i_max <= 0:
            raise ValueError("non-positive value for i_max")
        self.i_max = i_max

        self._index_field_width = int(np.log(i_max) / np.log(10.)) + 2
        self._min_f = None

    def __call__(self, x, f, a):
        if self._min_f is None:
            self._min_f = f

        s = (("{:>%d,d} | {:.6e} ({:.4f}) | " % self._index_field_width)
             .format(self.counter, f, f - self._min_f))
        for arg_val in x:
            s += "{:>16,.5e}".format(arg_val)
        print(s.replace(",", ""))

        self.counter += 1
        if f < self._min_f:
            self._min_f = f


class FitGeneralResult:
    W = 10
    SEP = "\t"

    def __init__(self, sample, _guesses, _result, _ids):
        self.sample = sample
        self.guesses = _guesses
        self.result = _result
        self.ids = _ids

        self._row = ["{:>%d}" % self.W] + (  # layer name
                ["{:>%d,.2e}{:>%d}" % (self.W, self.W - 3)] * len(FitGeneral.FIT_ARG_NAMES)
        )
        self._header = [" " * self.W] + (  # layer name
                ["{:>%d}{:>%d}" % (self.W, self.W - 3)] * len(FitGeneral.FIT_ARG_NAMES)
        )

    def __str__(self):
        vals = []
        errs = []
        for layer in self.sample.layers:
            vals.append([layer.__getattribute__(name[:-1]) for name in FitGeneral.FIT_ARG_NAMES])
            errs.append([0] * len(FitGeneral.FIT_ARG_NAMES))

        for j, id_pair in enumerate(self.ids):
            i_arg, i_layer = id_pair
            vals[i_layer][i_arg] = self.result[j]
            errs[i_layer][i_arg] = (self.result[j] - self.guesses[j]) / self.guesses[j]

        blank = self.SEP.join(self._row)
        h_blank = self.SEP.join(self._header)

        h_format_args = tuple()
        for name in FitGeneral.FIT_ARG_NAMES:
            h_format_args += (name, " ")

        ROWS = [h_blank.format(*h_format_args)]

        for i_layer, layer in enumerate(self.sample.layers):
            format_args = [layer.name]

            for i_arg in range(len(FitGeneral.FIT_ARG_NAMES)):
                format_args.append(vals[i_layer][i_arg])

                if errs[i_layer][i_arg] != 0:
                    p = errs[i_layer][i_arg] * 100
                    s = "({}{:.2f}%)".format('+' if p >= 0 else '-', abs(p))
                    format_args.append(s)
                else:
                    format_args.append("")
            ROWS.append(blank.format(*format_args))

        return "\n".join(ROWS)


class FitGeneral(Model):
    DEFAULT_GUESS = 100.0
    BOUNDARY_TYPES = ['s', 'i', 'a']
    FIT_ARG_NAMES = ["heights", "kys", "ratio_xys", "Cvs"]

    def __init__(self, sample, data, b_type):
        super().__init__(sample, data)
        self.n_layers = len(self.sample.layers)
        self._full_args = None

        self._fit_indices = []
        self._guesses = []
        self._ids = []

        self._bias = None
        self._defaults = (
            self.sample.heights, self.sample.kys, self.sample.ratio_xys, self.sample.Cvs
        )

        if b_type not in FitGeneral.BOUNDARY_TYPES:
            raise ValueError("boundary type {} is not one of {}"
                             .format(b_type, FitGeneral.BOUNDARY_TYPES))
        self.b_type = b_type

        self._full_args = []
        for lst in self._defaults:
            self._full_args += lst

        for i, arg in enumerate(self._full_args):
            if type(arg) is str and arg.endswith('*'):
                self._fit_indices.append(i)
                self._guesses.append(float(arg.rstrip('*')))
                self._ids.append(self._identify_fit_index(i))

        self._intg_set = False
        self._result = None
        self._fitted_kwargs = None
        self._error = None

        # C-extension integrator module
        self.intg = __import__('intg')

    @property
    def T2_fit(self):
        if self._result is None:
            return None
        T2 = self.T2_func(**self._fitted_kwargs)
        x = T2.real
        y = T2.imag
        xerr = None
        yerr = None  # TODO: an error estimate
        return ACReading(x, y, xerr, yerr)

    @property
    def result(self):
        """fitted kwargs for `T2_func` method"""
        return self._result

    @property
    def error(self):
        """mean-squared error for fit"""
        return self._error

    @property
    def defaults(self):
        return self._defaults

    def fit(self, niter=30, **kwargs):
        kwargs["guesses"] = self._guesses
        stepper = Stepper.by_fraction(**kwargs)
        callback = BasicPrinterCallBack(niter)
        fit_result = basinhopping(self.T2_err_func, self._guesses, niter=niter,
                                  minimizer_kwargs={
                                      'method': 'L-BFGS-B',
                                      'bounds': stepper.get_bounds(**kwargs)
                                  },
                                  take_step=stepper,
                                  callback=callback)
        self._record_result(fit_result)

    def T2_func(self, heights, kys, ratio_xys, Cvs) -> np.ndarray:
        """T2 prediction from physical model and provided properties"""
        if not self._intg_set:
            self._set_intg()

        # extra divisor of ROOT2 since measured T2 amplitudes are RMS
        P = -1 / (np.pi * self.heater.length * kys[0] * ROOT2) * self.power.x
        return P * self.intg.integral(heights, kys, ratio_xys, Cvs)

    def T2_err_func(self, args) -> float:
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
        # err += sum(np.abs(self.T2.y - T2_func_.imag))
        return err / len(T2_func_)

    def plot_fit(self, show=False):
        """plot fit result"""
        return plot.plot_fitted_T2(self, show=show)

    def set_data_limits(self, start, end):
        super().set_data_limits(start, end)

    def _identify_fit_index(self, index) -> tuple:
        i_arg = index // len(self.sample.heights)
        i_layer = index - i_arg * len(self.sample.heights)
        return i_arg, i_layer

    def _record_result(self, fit_result):
        fitted_argv = fit_result.x
        result = [
            ("heights", self.sample.heights.copy()),
            ("kys", self.sample.kys.copy()),
            ("ratio_xys", self.sample.ratio_xys.copy()),
            ("Cvs", self.sample.Cvs.copy())
        ]
        for i, index in enumerate(self._fit_indices):
            i_source = index // len(self.sample.heights)
            i_field = index - i_source * len(self.sample.heights)
            result[i_source][1][i_field] = fitted_argv[i]

        self._fitted_kwargs = {r[0]: r[1] for r in result}
        self._error = fit_result.fun
        self._result = FitGeneralResult(self.sample, self._guesses, fit_result.x, self._ids)

    def _set_intg(self):
        self.intg.set(self.omegas,
                      self.heater.width / 2,
                      1e-3,
                      1e7,
                      len(self.sample.layers),
                      self.b_type.encode('utf-8'))
        self._intg_set = True
