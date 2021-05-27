import numpy as np
from warnings import warn

import fit3omega.plot as plot
from fit3omega.data import ACReading
from fit3omega.model import Model


class FitLinearResult:
    W = 10
    SEP = "\t"

    def __init__(self, sample, fitted_kwargs):
        self.sample = sample
        self.kwargs = fitted_kwargs

        self._row = ["{:>%d}" % self.W] + (  # layer name
                ["{:>%d,.2e}{:>%d}" % (self.W, self.W - 3)] * len(sample.layers[:2])
        )
        self._header = [" " * self.W] + (  # layer name
                ["{:>%d}{:>%d}" % (self.W, self.W - 3)] * len(sample.layers[:2])
        )

    def __str__(self):
        blank = self.SEP.join(self._row)
        h_blank = self.SEP.join(self._header)
        COL_NAMES = self.kwargs.keys()

        h_format_args = tuple()
        for name in COL_NAMES:
            h_format_args += (name, " ")
        ROWS = [h_blank.format(*h_format_args)]

        for i_layer, layer in enumerate(self.sample.layers[:2]):
            format_args = [layer.name]
            for k, v in self.kwargs.items():
                guess = self._parse_param(layer.__getattribute__(k.rstrip('s')))
                val = self._parse_param(v[i_layer])
                err = (val - guess) / guess

                format_args.append(val)
                if err != 0:
                    p = err * 100.
                    s = "({}{:.2f}%)".format('+' if p >= 0 else '-', abs(p))
                    format_args.append(s)
                else:
                    format_args.append("")
            ROWS.append(blank.format(*format_args))

        return "\n".join(ROWS)

    @staticmethod
    def _parse_param(x):
        if type(x) is str:
            return float(x.rstrip('*'))
        return float(x)


class FitLinear(Model):
    eta = 3 / 2 - np.euler_gamma

    def __init__(self, sample, data, warnings=True, thresh=0.99, min_length=5):
        self.warnings = warnings
        self._min_length = min_length
        self._thresh = thresh
        super().__init__(sample, data)

        # count layers
        self.n_layers = len(self.sample.layers)
        if self.n_layers > 2:
            self.warn("only top 2/%d layers will be considered in linear model"
                      % self.n_layers)

        self._fitted_kwargs = {}
        self._result = None
        self._error = None
        self._i_linear = None
        self._T2_linear = None

        # linear coefficients
        self._coeffs_re, self._Rsq_re = None, None
        self._coeffs_im, self._Rsq_im = None, None

    @property
    def fitted_kwargs(self):
        return self._fitted_kwargs.copy()

    @property
    def thresh(self):
        return self._thresh

    @thresh.setter
    def thresh(self, t: float):
        if 0 < t < 1:
            self._thresh = t
        else:
            raise ValueError("Rsq threshold t must satisfy 0 < t < 1")

    @property
    def min_length(self):
        return self._min_length

    @min_length.setter
    def min_length(self, ml):
        if ml > 0:
            self._min_length = ml
        else:
            raise ValueError("minimum length must be a positive integer")

    @property
    def Rsq(self):
        return self._Rsq_re[self._i_linear]

    @property
    def T2_linear(self) -> ACReading:
        return self._T2_linear

    @property
    def omegas_linear(self) -> np.ndarray:
        return self.omegas[:self._i_linear]

    @property
    def T2_fit(self) -> ACReading:
        x = np.log(self.omegas_linear)
        m_re, b_re = self._coeffs_re[self._i_linear]
        m_im, b_im = self._coeffs_im[self._i_linear]
        return ACReading(m_re * x + b_re,
                         m_im * x + b_im,
                         None, None)

    @property
    def result(self):
        return self._result

    @property
    def error(self):
        return self._error

    def fit(self):
        self._calculate_linearity()

        p = self.power.x[:self._i_linear]
        l = self.heater.length
        b = self.heater.width / 2.
        a1, a0 = self._coeffs_re[self._i_linear]
        eta = FitLinear.eta

        # heater-on-substrate sample
        if len(self.sample.layers) == 1:
            subs = self.sample.layers[0]
            rS = subs.ratio_xy
            kSy = -p / (2. * np.pi * l * a1 * np.sqrt(rS))  # substrate conductivity
            if type(subs.Cv) is str:
                Cv = kSy / b**2 * np.exp((2. * p * eta - 2 * np.pi * l * kSy * a0) / p)
                Cv = np.average(Cv)
            else:
                Cv = subs.Cv
            self._fitted_kwargs = {
                "heights": [subs.height],
                "kys": [np.average(kSy)],
                "ratio_xys": [rS],
                "Cvs": [Cv],
            }
        # heater-on-film-on-substrate sample
        else:
            film, subs = self.sample.layers[:2]
            rS = subs.ratio_xy
            root_rS = np.sqrt(rS)
            kSy = -p / (2. * np.pi * l * a1 * root_rS) if type(subs.ky) is str else subs.ky
            CvS = subs.Cv
            TS = p / (np.pi * l * kSy * root_rS) * (
                    .5 * np.log((kSy * root_rS) / (CvS * b**2))
                    - .5 * np.log(self.omegas_linear)
                    + eta
            )

            dF = film.height
            kF = np.average(p * dF / (2. * b * l) / (self.T2_linear.x - TS))
            CvF = film.Cv

            self._fitted_kwargs = {
                "heights": [layer.height for layer in self.sample.layers[:2]],
                "kys": [kF, np.average(kSy)],
                "ratio_xys": [1.0, rS],
                "Cvs": [CvF, CvS]
            }

        self._result = FitLinearResult(self.sample, self._fitted_kwargs)
        self._error = sum((self.T2_fit.x - self.T2_linear.x)**2) / len(self.T2_fit.x)
        return

    def warn(self, msg):
        if self.warnings:
            warn(msg)

    def plot_fit(self, show=False):
        """plot fit result"""
        return plot.plot_fitted_T2_linear(self, show=show)

    def set_data_limits(self, start, end):
        super().set_data_limits(start, end)

    @staticmethod
    def Rsq_func(line_y, points_y):
        """data points Y; fit values y"""
        mean_points_y = np.average(points_y)
        SS_res = sum((line_y - points_y)**2)
        SS_tot = sum((points_y - mean_points_y)**2)
        if SS_tot == 0.0 and SS_res < np.finfo(float).eps:
            return 0.
        return 1. - SS_res / SS_tot

    def _calculate_linearity(self):
        coeffs_re, Rsq_re = [None, ], [None, ]
        coeffs_im, Rsq_im = [None, ], [None, ]

        N = len(self.omegas)
        x = np.log(self.omegas)
        for i in range(N - 1):
            i_ = i + 2
            xi = x[:i_]
            yi_re = self.T2.x[:i_]
            yi_im = self.T2.y[:i_]
            mb_re = np.polyfit(xi, yi_re, 1)
            mb_im = np.polyfit(xi, yi_im, 1)
            y_line_re = mb_re[0] * xi + mb_re[1]
            y_line_im = mb_im[0] * xi + mb_im[1]
            Rsq_re.append(FitLinear.Rsq_func(y_line_re, yi_re))
            Rsq_im.append(FitLinear.Rsq_func(y_line_im, yi_im))
            coeffs_re.append(mb_re)
            coeffs_im.append(mb_im)

        k = self._min_length - 1
        while k < N and Rsq_re[k] > self.thresh:
            k += 1
        self._i_linear = k - 1
        self._coeffs_re, self._Rsq_re = coeffs_re, Rsq_re
        self._coeffs_im, self._Rsq_im = coeffs_im, Rsq_im
        self._T2_linear = ACReading(self.T2.x[:self._i_linear],
                                    self.T2.y[:self._i_linear],
                                    self.T2.xerr[:self._i_linear],
                                    self.T2.yerr[:self._i_linear])


if __name__ == "__main__":
    SMPL = "/Users/araghukasyan/3omegaPi/control_4/control_4.f3oc"
    DATA = "/Users/araghukasyan/3omegaPi/control_4/May26_2021_m3/tc3omega_data_2.8_V.csv"
    fl = FitLinear(SMPL, DATA, thresh=0.9)
    fl.set_data_limits(4, 50)
    fl.fit()
    print(fl.result)
    fl.plot_fit(show=True)

    # import matplotlib.pyplot as plt
    #
    # plt.plot(fl.omegas, fl.T2.x)
    # plt.plot(fl.omegas, fl.T2.y)
    #
    # plt.plot(fl.omegas_linear, fl.T2_linear.x, linestyle='None', marker='s')
    # plt.plot(fl.omegas_linear, fl.T2_linear.y, linestyle='None', marker='s')
    #
    # plt.plot(fl.omegas_linear, fl.T2_fit.x, color='b')
    # plt.plot(fl.omegas_linear, fl.T2_fit.y, color='g')
    # plt.xscale("log")
    # plt.show()
