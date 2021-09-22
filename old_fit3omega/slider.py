import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider
from typing import Tuple
from math import pi as PI
import numpy as np

from old_fit3omega.base import _set_mpl_defaults
from old_fit3omega.varsample import VarSample
from old_fit3omega.fit import FitGeneral


class SliderPlot(object):
    """
    An interactive plot providing a real-time visual for manual data fitting.
    """
    meas_plot_kw = dict(
        markersize=5,
        marker='o',
        markerfacecolor='w',
        zorder=0,
        linestyle='None'
    )

    fit_plot_kw = dict(
        linestyle='-',
        linewidth=1.5,
        markersize=0
    )

    # left, bottom, width, height; dimensions normalized (0,1)
    button_hovercolor = "0.98"

    slider_start_dims = [0.2, 0.38, 0.6, 0.01]
    slider_delta = [0.0, -0.035, 0.0, 0.0]
    slider_valfmt = "%.2e"

    error_fmt = "error: {:<10,.6e}"
    error_green_thresh = 0.005

    def __init__(self, fitter: FitGeneral, frac: float = 0.99, niter: int = None,
                 enable_heater_params: bool = False):
        mpl.rc("font", family="monospace")
        self.frac = frac
        self.niter = niter
        self.enable_heater_params = enable_heater_params
        fitter.model.sample = VarSample.from_sample(fitter.model.sample)

        self.fitter = fitter
        self.fitter.model.refresh = True
        self.model_type = type(self.fitter)
        self.fig, self.ax = plt.subplots(figsize=(8, 8))
        self.fig.subplots_adjust(bottom=.5, top=0.95)
        self.param_sliders = {}
        self.sample_sliders = {}
        self.buttons = {}

        self._x_scale = "log(x)"
        self._n_sliders = 0
        self._sample_values_contain_strings = True

    def plot_initial_state(self):
        _set_mpl_defaults()

        self.ax.set_xscale("log")
        self.ax.set_ylabel(r"Sample $\widebar{T}_{2\omega}$")
        self.ax.set_xlabel(r"$2\pi\omega$ [Hz]")

        X = self.fitter.model.omegas / 2. / PI
        cx = 'blue'
        cy = 'red'
        fcx = u'#0B554C'
        fcy = u'#550B14'

        # measured data
        self.ax.plot(X, self.fitter.model.T2.x, color=cx, **self.meas_plot_kw)
        self.ax.plot(X, self.fitter.model.T2.y, color=cy, **self.meas_plot_kw)

        # fitted data (initial values)
        fit = self.get_fitline_data()
        self.ax.plot(X, fit[0], color=fcx, **self.fit_plot_kw)
        self.ax.plot(X, fit[1], color=fcy, **self.fit_plot_kw)

        # prepare the fit parameters sliders
        for i, layer in enumerate(self.fitter.model.sample.layers):
            d = layer.as_dict()
            for k, v in d.items():
                k_ = k + 's'
                if type(v) is str and v.endswith('*') and k_ in VarSample.FIT_PARAMS:
                    name = layer.name
                    label = name + '.' + k
                    guess_val = float(v.rstrip('*'))
                    axes = plt.axes(self._get_slider_dims())
                    self.param_sliders[label] = Slider(
                        ax=axes,
                        label=label.replace('.', ' '),
                        valmin=(1 - self.frac) * guess_val,
                        valmax=(1 + self.frac) * guess_val,
                        valinit=guess_val,
                        valfmt=self.slider_valfmt
                    )
                    self.param_sliders[label].on_changed(self._apply_sliders)
                    self._n_sliders += 1

        # prepare setup parameter sliders (Rsh, dRdT, length, width)
        if self.enable_heater_params:
            heater_params = dict(
                dRdT=self.fitter.model.sample.heater.dRdT,
                # width=self.model.sample.heater.width,
                # length=self.model.sample.heater.length,
                Cv=self.fitter.model.sample.heater.Cv,
                height=self.fitter.model.sample.heater.height,
                Rc=self.fitter.model.sample.heater.Rc
            )
            for k, v in heater_params.items():
                if v == 0:
                    continue
                self.sample_sliders[k] = Slider(
                    ax=plt.axes(self._get_slider_dims()),
                    label=k,
                    valmin=(1 - self.frac) * v,
                    valmax=(1 + self.frac) * v,
                    valinit=v,
                    valfmt=self.slider_valfmt
                )
                self.sample_sliders[k].on_changed(self._apply_sliders)
                self._n_sliders += 1

        # TODO: clean up hard-coded nonsense below
        # sloppy reset button creation
        reset_button_dims = self._get_slider_dims()
        reset_button_dims[2:] = [0.1, 0.03]
        reset_button_dims[1] -= 0.035
        self.buttons["Reset"] = Button(plt.axes(reset_button_dims), "Reset",
                                       hovercolor=self.button_hovercolor)
        self.buttons["Reset"].on_clicked(self._reset_sliders)

        # sloppy x-scale toggle creation
        xscale_button_dims = reset_button_dims
        xscale_button_dims[0] += 0.12
        self.buttons["x-scale"] = Button(plt.axes(xscale_button_dims), self._x_scale,
                                         hovercolor=self.button_hovercolor)
        self.buttons["x-scale"].on_clicked(self._toggle_xscale)

        # fit button
        fit_button_dims = xscale_button_dims
        fit_button_dims[0] += 0.12
        fit_button_dims[2] += 0.005
        self.buttons["Fit"] = Button(plt.axes(fit_button_dims), "Fit",
                                     hovercolor=self.button_hovercolor)
        self.buttons["Fit"].on_clicked(self._run_fit_and_update)

        self.ax.text(0.0, 1.025, self.error_fmt.format(fit[2]), transform=self.ax.transAxes,
                     fontsize=10, color=self._get_error_color(fit[2]), fontweight="bold")
        plt.show()

    def get_fitline_data(self) -> Tuple[np.ndarray, np.ndarray, float]:
        """return a tuple: fitline_real, fitline_imag, fitline_err"""
        T2_kwargs = self.fitter.get_T2_kwargs()
        T2_func_values = self.fitter.T2_func(**T2_kwargs)
        dx = T2_func_values.real - self.fitter.model.T2.x
        dy = T2_func_values.imag - self.fitter.model.T2.y
        err = sum((dx**2 + dy**2) / self.fitter.model.T2.norm_sq) / len(dx)
        return T2_func_values.real, T2_func_values.imag, err

    def _apply_sliders(self, _):
        for label, slider in self.param_sliders.items():
            layer_name, attr_name = label.split('.')
            self.fitter.model.sample.param_modify(layer_name, attr_name + 's', slider.val)
        for label, slider in self.sample_sliders.items():
            self.fitter.model.sample.heater.modify(label, slider.val)
        self._update_graph()

    def _update_graph(self):
        fit = self.get_fitline_data()
        self.ax.lines[0].set_ydata(self.fitter.model.T2.x)
        self.ax.lines[1].set_ydata(self.fitter.model.T2.y)
        self.ax.lines[2].set_ydata(fit[0])
        self.ax.lines[3].set_ydata(fit[1])
        self.ax.texts[0].set_text(self.error_fmt.format(fit[2]))
        self.ax.texts[0].set_color(self._get_error_color(fit[2]))
        self.fig.canvas.draw_idle()

    def _reset_sliders(self, _):
        for label in self.param_sliders:
            self.param_sliders[label].reset()
        for label in self.sample_sliders:
            self.sample_sliders[label].reset()
        self.fitter.model.sample.reset_params()
        self.fitter.model.heater.reset()

    def _toggle_xscale(self, _):
        if self._x_scale == "log(x)":
            self.ax.set_xscale("linear")
            self._x_scale = "x"
        else:
            self.ax.set_xscale("log")
            self._x_scale = "log(x)"

        self.buttons["x-scale"].label.set_text(self._x_scale)

    def _get_slider_dims(self):
        dims = []
        n = self._n_sliders
        for x, d in zip(self.slider_start_dims, self.slider_delta):
            dims.append(x + n * d)
        return dims

    def _get_error_color(self, err):
        t = self.error_green_thresh
        if err > 1.0:
            return 1.0, 0.0, 0.0
        elif t < err < 1.0:
            return (err - t), 0.0, 0.0
        elif 0.0 <= err <= t:
            return 0.0, (t - err) / t, 0.0
        else:
            return 0.0, 0.0, 0.0

    def _run_fit_and_update(self, _):
        self.fitter.fit()
        for attr_name, values in self.fitter.result.fitted_kwargs.items():
            self.fitter.model.sample.param_modify(None, attr_name, values)
        self._update_graph()
        print("--> Result:\n", self.fitter.result.summary(), "\n")

        for label, slider in self.param_sliders.items():
            layer_name, attr_name = label.split('.')
            for i, layer in enumerate(self.fitter.model.sample.layers):
                if layer.name == layer_name:
                    slider.set_val(self.fitter.result.fitted_kwargs[attr_name + "s"][i])
                    break
