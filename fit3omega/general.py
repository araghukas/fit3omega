import numpy as np
from scipy.optimize import basinhopping, Bounds

from fit3omega.model import Model
from fit3omega.data import ACReading


def positive_bounds(guesses: list) -> Bounds:
    lb, ub = [], []
    for g in guesses:
        lb.append(0.1 * g)
        ub.append(1.9 * g)
    return Bounds(lb, ub, keep_feasible=True)


class TakeStep:
    def __init__(self, guesses, frac):
        self.guesses = guesses
        self.frac = frac if frac < 1.0 else 0.99
        self.step_sizes = [frac * guess for guess in guesses]

    def __call__(self, x):
        for i in range(len(x)):
            s = self.step_sizes[i]
            x[i] += np.random.uniform(-s, s)
        return x


class BasicPrinterCallBack:
    """
    Handles neatly printing results during fitting iterations
    """
    MIN_F_THRESH = 1e-6

    def __init__(self, niter, quiet=True, silent=False):
        if niter <= 0:
            raise ValueError("non-positive value for i_max")
        self.i_max = niter

        self._silent = silent
        self._quiet = quiet

        self._counter = 1
        self._idx_col_width = int(np.log(niter) / np.log(10.)) + 2
        self._min_f = None
        self._min_counter = 1

    def __call__(self, x, f, a):
        if self._silent:
            return

        if self._min_f is None:
            self._min_f = f
            is_new_min = False
        else:
            is_new_min = (self._min_f - f) > self.MIN_F_THRESH

        g = "-" if is_new_min else "+"
        s = " | ".join([
            ("{:>%d}" % self._idx_col_width).format(self._counter),
            ("{:.3e} (%s{:.5e})" % g).format(f, abs(f - self._min_f)),
            "".join(["{:>16,.3e}".format(arg_val) for arg_val in x])
        ])

        if is_new_min:
            s += " min %d" % self._min_counter
            self._min_counter += 1
            self._min_f = f

        if self._quiet:
            if is_new_min:
                print(s)
        else:
            print(s)

        self._counter += 1


class BasinhoppingOptimizer(Model):
    FIT_PARAMS = []

    def __init__(self, sample, data):
        super().__init__(sample, data)
        self.n_layers = len(self.sample.layers)
        self._full_args = None
        self._config_values = None

        self._fit_indices = []
        self._guesses = []
        self._ids = []

        self._full_args = []
        for lst in self.config_values:
            self._full_args += lst

        for i, arg in enumerate(self._full_args):
            if type(arg) is str:
                if arg.endswith('*'):
                    self._fit_indices.append(i)
                    self._guesses.append(float(arg.rstrip('*')))
                    self._ids.append(self._identify_fit_index(i))
                else:
                    raise ValueError("invalid sample parameter '%s'" % arg)

        self._integrators_ready = False
        self._result = None
        self._fitted_kwargs = {}
        self._error = None

        self._T2_func_latest = None
        self._min_err = None

        # C-extension integrator module
        self.integrators = __import__('integrate')

    def _identify_fit_index(self, index) -> tuple:
        n_layers = len(self.sample.heights)
        i_arg = index // n_layers
        i_layer = index - i_arg * n_layers
        return i_arg, i_layer

    @property
    def config_values(self):
        # lists of starting values for applicable fit parameters
        if self._config_values is None:
            self._config_values = tuple(self.sample.__getattribute__(k) for k in self.FIT_PARAMS)
        return self._config_values

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
        xerr = 0.0
        yerr = 0.0  # TODO: an error estimate
        return ACReading(x, y, xerr, yerr)

    @property
    def result(self):
        return self._result

    @property
    def error(self):
        """mean-squared error for fit"""
        return self._error

    def error_func(self, args) -> float:
        """objective function for the fit method"""
        args_T2 = self._sub_args_into_complete_params(args)
        T2_func_ = self.T2_func(*args_T2)

        err = sum(
            ((T2_func_.real - self.T2.x)**2 + (T2_func_.imag - self.T2.y)**2)
            / self.T2.norm_sq
        )

        return err / (len(T2_func_))

    def fit(self, niter=30, niter_success=1000, frac=0.9, tol=1e-1,
            quiet=True, silent=False, min_err=None, T=None):

        self._min_err = min_err
        callback = BasicPrinterCallBack(niter, quiet=quiet, silent=silent)
        minimizer_kwargs = self._insert_extra_minimizer_kwargs({
            'method': 'L-BFGS-B',
            'bounds': positive_bounds(self._guesses),
            'tol': tol
        })
        take_step = TakeStep(self._guesses, frac)
        fit_result = basinhopping(self.error_func, self._guesses, niter=niter,
                                  minimizer_kwargs=minimizer_kwargs,
                                  callback=callback,
                                  take_step=take_step,
                                  niter_success=niter_success,
                                  accept_test=None if min_err is None else self._force_improvement,
                                  T=1e-4 if T is None else 0)
        self._record_result(fit_result)

    def _force_improvement(self, f_new, **ignored_kwargs):
        if f_new < self._min_err:
            return True
        return False

    def _sub_args_into_complete_params(self, args) -> tuple:
        """produces a tuple of all arguments for `T2_func` from a partial args list"""
        err_args = self._full_args.copy()
        for j, index in enumerate(self._fit_indices):  # `self._fit_indices` provides mapping
            err_args[index] = args[j]

        # reconstruct T2_func args
        args_T2 = tuple()
        for k in range(len(self._full_args) // self.n_layers):
            i_min = k * self.n_layers
            i_max = i_min + self.n_layers
            args_T2 += (err_args[i_min:i_max],)

        return args_T2

    def _record_result(self, fit_result):
        fitted_argv = fit_result.x
        result = self._get_initial_result()
        for i, index in enumerate(self._fit_indices):
            i_source = index // len(self.sample.heights)
            i_field = index - i_source * len(self.sample.heights)
            result[i_source][1][i_field] = fitted_argv[i]

        self._fitted_kwargs = {r[0]: r[1] for r in result}
        self._error = fit_result.fun
        self._result = OptimizerResult(self.sample, self._guesses, fit_result.x, self._ids,
                                       self.FIT_PARAMS)

    def _get_initial_result(self):
        return [(name, self.sample.__getattribute__(name).copy()) for name in self.FIT_PARAMS]

    # override methods below this line
    def T2_func(self, *args, **kwargs):
        raise NotImplementedError

    def plot_fit(self, **kwargs):
        raise NotImplementedError

    def _init_integrators(self):
        raise NotImplementedError

    def get_current_T2(self):
        raise NotImplementedError

    def _insert_extra_minimizer_kwargs(self, kwargs: dict) -> dict:
        raise NotImplementedError


class OptimizerResult:
    """
    A probably-unnecessary container for the fit result.
    This class mostly just enables consistently formatted output.
    """
    W = 10
    SEP = "\t"

    def __init__(self, sample, _guesses, x, _ids, fit_arg_names):
        self.sample = sample
        self.guesses = _guesses
        self.x = x
        self.ids = _ids
        self.fit_arg_names = fit_arg_names

        self._row = ["{:>%d}" % self.W] + (  # layer name
                ["{:>%d,.2e}{:>%d}" % (self.W, self.W - 3)] * len(self.fit_arg_names)
        )
        self._header = [" " * self.W] + (  # layer name
                ["{:>%d}{:>%d}" % (self.W, self.W - 3)] * len(self.fit_arg_names)
        )

    def __str__(self):
        vals = []
        errs = []
        for layer in self.sample.layers:
            vals.append([layer.__getattribute__(name[:-1]) for name in self.fit_arg_names])
            errs.append([0] * len(self.fit_arg_names))

        for j, id_pair in enumerate(self.ids):
            i_arg, i_layer = id_pair
            vals[i_layer][i_arg] = self.x[j]
            errs[i_layer][i_arg] = (self.x[j] - self.guesses[j]) / self.guesses[j]

        blank = self.SEP.join(self._row)
        h_blank = self.SEP.join(self._header)

        h_format_args = tuple()
        for name in self.fit_arg_names:
            h_format_args += (name, " ")

        ROWS = [h_blank.format(*h_format_args)]

        for i_layer, layer in enumerate(self.sample.layers):
            format_args = [layer.name]

            for i_arg in range(len(self.fit_arg_names)):
                format_args.append(vals[i_layer][i_arg])

                if errs[i_layer][i_arg] != 0:
                    p = errs[i_layer][i_arg] * 100
                    s = "({}{:.2f}%)".format('+' if p >= 0 else '-', abs(p))
                    format_args.append(s)
                else:
                    format_args.append("")
            ROWS.append(blank.format(*format_args))

        return "\n".join(ROWS)
