import numpy as np
from scipy.optimize import basinhopping

from fit3omega.model import Model


class Stepper:
    def __init__(self, step_sizes):
        self.step_sizes = step_sizes

    def __call__(self, x):
        for i, s in enumerate(self.step_sizes):
            x[i] += np.random.uniform(-s, s)
        return x

    @classmethod
    def by_fraction(cls, guesses: list, frac: float):
        step_sizes = [guess * frac for guess in guesses]
        return cls(step_sizes)


def bounds_by_fraction(guesses: list, frac: float):
    f = frac / 2
    return [((1 - f) * guess, (1 + f) * guess) for guess in guesses]


class FitGeneral(Model):
    DEFAULT_GUESS = 100.0
    BOUNDARY_TYPES = ['s', 'i', 'a']

    METHODS = {
        "fraction": (Stepper.by_fraction, bounds_by_fraction)
    }

    def __init__(self, sample, data, b_type):
        super().__init__(sample, data)
        self.n_layers = len(self.sample.layers)
        self._full_args = None
        self._fit_indices = []
        self._guesses = []
        self._bias = None
        self._defaults = (self.sample.heights, self.sample.kxs, self.sample.kys, self.sample.Cvs)

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

        self._intg_set = False
        self._result = None
        self._error = None

        # C-extension integrator module
        self.intg = __import__('intg')

    @property
    def result(self):
        """fitted kwargs for `T2_func` method"""
        if self._result is None:
            raise ValueError("fit result has not been calculated")
        return self._result

    @property
    def error(self):
        """mean-squared error for fit"""
        if self._error is None:
            raise ValueError("fit result has not been calculated")
        return self._error

    @property
    def defaults(self):
        return self._defaults

    @property
    def bias(self):
        if self._bias is None:
            self._bias = self.omegas / np.max(self.omegas)
        return self._bias

    def T2_func(self, heights, kxs, kys, Cvs) -> np.ndarray:
        """T2 prediction from physical model and provided properties"""
        if not self._intg_set:
            self._set_intg()

        P = -1 / (np.pi * self.heater.length * kys[0]) * self.power.phasor()
        return P * self.intg.integral(heights, kxs, kys, Cvs)

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

        err = np.dot(np.abs(self.T2.phasor() - T2_func_), self.bias)
        return err / len(T2_func_)

    def fit(self, method: str = None, **kwargs):
        if method is None:
            method = "fraction"
            kwargs = dict(frac=0.1)
        elif method not in FitGeneral.METHODS:
            raise ValueError("unknown method '%s'" % method)

        """MAIN FIT FUNCTION"""
        stepper, bound_func = FitGeneral.METHODS[method]
        kwargs["guesses"] = self._guesses
        fit_result = basinhopping(self.T2_err_func, self._guesses, niter=10,
                                  minimizer_kwargs={
                                      'method': 'L-BFGS-B',
                                      'bounds': bound_func(**kwargs)
                                  },
                                  take_step=stepper(**kwargs))
        self._record_result(fit_result)

    def plot(self):
        """plot fit result"""
        # result = self.result
        # T2_func_ = self.T2_func(**self.result)
        raise NotImplementedError

    def _identify_fit_index(self, index) -> str:
        prop_names = ["height", "kx", "ky", "Cv"]
        i_prop = index // len(self.defaults[0])
        i_name = index - i_prop * len(self.defaults[0])
        name = self.sample.layers[i_name].name
        prop = prop_names[i_prop]
        return '.'.join([name, prop])

    def _record_result(self, fit_result):
        fitted_argv = fit_result.x
        result = [
            ("heights", self.sample.heights.copy()),
            ("kxs", self.sample.kxs.copy()),
            ("kys", self.sample.kys.copy()),
            ("Cvs", self.sample.Cvs.copy())
        ]
        for i, index in enumerate(self._fit_indices):
            i_source = index // len(self.sample.heights)
            i_field = index - i_source * len(self.sample.heights)
            result[i_source][1][i_field] = fitted_argv[i]

        self._result = {r[0]: r[1] for r in result}
        self._error = fit_result.fun

    def _set_intg(self):
        self.intg.set(self.omegas,
                      self.heater.width / 2,
                      1e-3,
                      1e7,
                      len(self.sample.layers),
                      self.b_type.encode('utf-8'))
        self._intg_set = True


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    sample_ = "../sample/2232_2.f3oc"
    data_ = "../sample/tc3omega_data_3.0_V.csv"
    g = FitGeneral(sample_, data_, 'i')
    g.data.drop_row(42)
    g.data.drop_row(48)

    g.fit("fraction", frac=1.0)
    for k1, v1 in g.result.items():
        print(k1, v1)
    print(g.error)

    Xm = g.T2.x
    Ym = g.T2.y
    Rm = np.sqrt(Xm**2 + Ym**2)

    Tf = g.T2_func(**g.result)
    Xf = Tf.real
    Yf = Tf.imag
    Rf = np.abs(Tf)

    fig, ax = plt.subplots()
    ax.plot(g.omegas, Xm, marker='o', markersize=5, markerfacecolor='white', color="red",
            linewidth=0)
    ax.plot(g.omegas, Ym, marker='o', markersize=5, markerfacecolor='white', color="blue",
            linewidth=0)
    ax.plot(g.omegas, Rm, marker='o', markersize=5, color="black",
            linewidth=0)

    ax.plot(g.omegas, Xf, color="red", linestyle=':')
    ax.plot(g.omegas, Yf, color="blue", linestyle=':')
    ax.plot(g.omegas, Rf, color='black', linewidth=1.5)
    ax.set_xscale('log')

    plt.show()
