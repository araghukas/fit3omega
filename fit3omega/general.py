import numpy as np
from scipy.optimize import basinhopping

import fit3omega.plot as plot
from fit3omega.model import Model, ROOT2
from fit3omega.data import ACReading


class Stepper:
    """
    Controls random walk used in "basinhopping" fit engine
    """

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
        # TODO: not quite right here
        f = frac / 2
        a1 = 1 - f
        a2 = 1 + f
        if a1 > a2:
            a1, a2 = a2, a1

        a2 = 0.0 if frac < 0 < a2 else a2  # neg. fraction clamp to [1+f  ,   0]
        a1 = 0.0 if a1 < 0 < frac else a1  # pos. fraction clamp to [    0, 1+f]

        return [(a1 * guess, a2 * guess) for guess in guesses]


class BasicPrinterCallBack:
    """
    Handles neatly printing results during fitting iterations
    """
    MIN_F_THRESH = 1e-6

    def __init__(self, i_max, silent=False):
        if i_max <= 0:
            raise ValueError("non-positive value for i_max")
        self.i_max = i_max

        self._counter = 1
        self._idx_col_width = int(np.log(i_max) / np.log(10.)) + 2
        self._min_f = None
        self._min_counter = 1

    def __call__(self, x, f, a):
        if self._min_f is None:
            self._min_f = f
            is_new_min = False
        else:
            is_new_min = (self._min_f - f) > self.MIN_F_THRESH

        g = "-" if is_new_min else "+"
        s = " | ".join([
            ("{:>%d}" % self._idx_col_width).format(self._counter),
            ("{:.3e} (%s{:.5f})" % g).format(f, abs(f - self._min_f)),
            "".join(["{:>16,.3e}".format(arg_val) for arg_val in x])
        ])

        if is_new_min:
            h_line = "-" * len(s)
            s = "\n".join([h_line, s + " min %d" % self._min_counter, h_line])
            self._min_counter += 1
            self._min_f = f

        print(s)
        self._counter += 1


class BasinhoppingOptimizer(Model):
    FIT_ARG_NAMES = []

    def __init__(self, sample, data):  # , b_type):
        super().__init__(sample, data)
        self.n_layers = len(self.sample.layers)
        self._full_args = None

        self._fit_indices = []
        self._guesses = []
        self._ids = []

        self._full_args = []

        for i, arg in enumerate(self._full_args):
            if type(arg) is str and arg.endswith('*'):
                self._fit_indices.append(i)
                self._guesses.append(float(arg.rstrip('*')))
                self._ids.append(self._identify_fit_index(i))

        self._integrators_ready = False
        self._result = None
        self._fitted_kwargs = {}
        self._error = None

        # C-extension integrator module
        self.integrators = __import__('integrate')

    def _identify_fit_index(self, index) -> tuple:
        i_arg = index // len(self.sample.heights)
        i_layer = index - i_arg * len(self.sample.heights)
        return i_arg, i_layer

    @property
    def fitted_kwargs(self):
        return self._fitted_kwargs.copy()

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
        return self._result

    @property
    def error(self):
        """mean-squared error for fit"""
        return self._error

    def set_data_limits(self, start, end):
        super().set_data_limits(start, end)

    # override methods below this line
    @property
    def defaults(self):
        raise NotImplementedError

    def T2_func(self, **kwargs):
        raise NotImplementedError

    def plot_fit(self, **kwargs):
        raise NotImplementedError

    def error_func(self, **kwargs):
        raise NotImplementedError

    def _init_integrators(self):
        raise NotImplementedError

    def _record_result(self, fit_result):
        raise NotImplementedError


class FitCBT(BasinhoppingOptimizer):
    BOUNDARY_TYPES = ['s', 'i', 'a']

    def __init__(self, sample, data, b_type):
        super().__init__(sample, data)

        self._defaults = (
            self.sample.heights, self.sample.kys, self.sample.ratio_xys, self.sample.Cvs
        )

        for lst in self._defaults:
            self._full_args += lst

        if b_type not in FitCBT.BOUNDARY_TYPES:
            raise ValueError("boundary type {} is not one of {}"
                             .format(b_type, FitCBT.BOUNDARY_TYPES))
        self.b_type = b_type

    def defaults(self):
        return self._defaults

    def T2_func(self, heights, kys, ratio_xys, Cvs) -> np.ndarray:
        """T2 prediction from physical model and provided properties"""
        if not self._intg_set or self._refresh_dependents:
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
        self._intg_set = True


class FitCBTResult:
    """
    A probably-unnecessary container for the fit result.
    This class mostly just enables consistently formatted output.
    """
    W = 10
    SEP = "\t"

    def __init__(self, sample, _guesses, x, _ids):
        self.sample = sample
        self.guesses = _guesses
        self.x = x
        self.ids = _ids

        self._row = ["{:>%d}" % self.W] + (  # layer name
                ["{:>%d,.2e}{:>%d}" % (self.W, self.W - 3)] * len(FitCBT.FIT_ARG_NAMES)
        )
        self._header = [" " * self.W] + (  # layer name
                ["{:>%d}{:>%d}" % (self.W, self.W - 3)] * len(FitCBT.FIT_ARG_NAMES)
        )

    def __str__(self):
        vals = []
        errs = []
        for layer in self.sample.layers:
            vals.append([layer.__getattribute__(name[:-1]) for name in FitCBT.FIT_ARG_NAMES])
            errs.append([0] * len(FitCBT.FIT_ARG_NAMES))

        for j, id_pair in enumerate(self.ids):
            i_arg, i_layer = id_pair
            vals[i_layer][i_arg] = self.x[j]
            errs[i_layer][i_arg] = (self.x[j] - self.guesses[j]) / self.guesses[j]

        blank = self.SEP.join(self._row)
        h_blank = self.SEP.join(self._header)

        h_format_args = tuple()
        for name in FitCBT.FIT_ARG_NAMES:
            h_format_args += (name, " ")

        ROWS = [h_blank.format(*h_format_args)]

        for i_layer, layer in enumerate(self.sample.layers):
            format_args = [layer.name]

            for i_arg in range(len(FitCBT.FIT_ARG_NAMES)):
                format_args.append(vals[i_layer][i_arg])

                if errs[i_layer][i_arg] != 0:
                    p = errs[i_layer][i_arg] * 100
                    s = "({}{:.2f}%)".format('+' if p >= 0 else '-', abs(p))
                    format_args.append(s)
                else:
                    format_args.append("")
            ROWS.append(blank.format(*format_args))

        return "\n".join(ROWS)