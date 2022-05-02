"""the matplotlib-based real time fitting GUI is defined here"""
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from typing import Tuple, Union, List
from matplotlib.widgets import Slider, Button
from math import pi as PI

from fit3omega.sample import SampleParameters
from fit3omega.data import Data
from fit3omega.fit import Fit3omega
from fit3omega.plots import _set_mpl_defaults


class SliderFit(Fit3omega):
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
    error_green_thresh = 0.01

    def __init__(self,
                 sample: Union[str, SampleParameters],
                 data: Union[str, Data],
                 frac: float = 0.99,
                 enable_heater_params: bool = True):
        super().__init__(sample, data)

        mpl.rc("font", family="monospace")
        self.frac = frac
        self.sample_sliders = {}
        self.heater_sliders = {}
        self.buttons = {}

        self._x_scale = "log(x)"
        self._n_sliders = 0

        # PLOT INITIAL STATE
        self._refresh_dependents = True
        self._enable_heater_params = enable_heater_params

    def show(self):
        """display the slider plot"""
        self._plot_initial_state(self._enable_heater_params)

    def _create_axes(self):
        self.fig, self.ax = plt.subplots(figsize=(8, 8))
        self.fig.subplots_adjust(bottom=.5, top=0.95)
        self.ax.set_xscale("log")
        self.ax.set_ylabel(r"$\widebar{T}_{2\omega}$")
        self.ax.set_xlabel(r"$\omega$ [Hz]")

    def _plot_markers_and_curves(self):
        X = self.data.omegas / 2. / PI
        cx = 'blue'
        cy = 'red'
        fcx = u'#0B554C'
        fcy = u'#550B14'

        # measured data
        self.ax.plot(X, self.T2.x, color=cx, **self.meas_plot_kw)
        self.ax.plot(X, self.T2.y, color=cy, **self.meas_plot_kw)

        # shade to show error ranges
        self.ax.fill_between(X,
                             self.T2.x * (1. - self.T2.xerr),
                             self.T2.x * (1. + self.T2.xerr),
                             alpha=0.2, color=cx)
        self.ax.fill_between(X,
                             self.T2.y * (1. - self.T2.yerr),
                             self.T2.y * (1. + self.T2.yerr),
                             alpha=0.2, color=cy)

        # fitted data curve (initial values)
        fit = self.T2_function(*self.sample.argv)
        self.ax.plot(X, fit.real, color=fcx, **self.fit_plot_kw)
        self.ax.plot(X, fit.imag, color=fcy, **self.fit_plot_kw)

    def _create_sliders(self, enable_heater_params):
        # prepare the fit parameters sliders
        for name_string, value in self.sample.parameters.items():
            axes = plt.axes(self._get_slider_dims())
            parameter_slider = Slider(
                ax=axes,
                label=name_string,
                valmin=(1 - self.frac) * value,
                valmax=(1 + self.frac) * value,
                valinit=value,
                valfmt=self.slider_valfmt
            )
            # register the new slider
            self.sample_sliders[name_string] = parameter_slider
            self.sample_sliders[name_string].on_changed(self._apply_sample_sliders)
            self._n_sliders += 1

        # prepare setup parameter sliders (Rsh, dRdT, length, width)
        if enable_heater_params:
            heater_params = dict(
                # dRdT=self.sample.heater.dRdT,
                # width=self.model.sample.heater.width,
                # length=self.model.sample.heater.length,
                Cv=self.sample.heater.Cv,
                height=self.sample.heater.height,
                Rc=self.sample.heater.Rc
            )
            for k, v in heater_params.items():
                if v == 0:
                    continue
                self.heater_sliders[k] = Slider(
                    ax=plt.axes(self._get_slider_dims()),
                    label=k,
                    valmin=(1 - self.frac) * v,
                    valmax=(1 + self.frac) * v,
                    valinit=v,
                    valfmt=self.slider_valfmt
                )
                self.heater_sliders[k].on_changed(self._apply_heater_sliders)
                self._n_sliders += 1

    def _create_buttons(self):
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

        # save state button
        save_button_dims = fit_button_dims
        save_button_dims[0] += 0.12
        self.buttons["Save"] = Button(plt.axes(save_button_dims), "Save",
                                      hovercolor=self.button_hovercolor)
        self.buttons["Save"].on_clicked(self._save_sample_state)

        err0 = self.objective_func(self.sample.x)
        self.ax.text(0.0, 1.025, self.error_fmt.format(err0), transform=self.ax.transAxes,
                     fontsize=10, color=self._get_error_color(err0), fontweight="bold")
        plt.show()

    def _plot_initial_state(self, enable_heater_params: bool) -> None:
        """prepare the tool and show the window"""
        _set_mpl_defaults()
        self._create_axes()
        self._plot_markers_and_curves()
        self._create_sliders(enable_heater_params)
        self._create_buttons()

    def _run_fit_and_update(self, _) -> None:
        """fit from current position and update the fit information with the result"""
        self.fit()
        x = self.result.x

        # set sample parameters to fitted values
        i = 0
        for k in self.sample.parameters.keys():
            param_name, layer_name = k.split('.')
            self.sample.modify_layer(layer_name, param_name, x[i])
            self.sample_sliders[k].set_val(x[i])
            i += 1

        # update the graph to show fitted curve
        print(self.result)
        self._update_graph()

    def _save_sample_state(self, _) -> None:
        """save the sample state into a new config file"""
        default_name = "./saved_state{}.txt"
        save_name = default_name.format("")
        i = 1
        while os.path.isfile(save_name):
            save_name = default_name.format("(%d)" % i)
            i += 1

        self.sample.write_state(save_name)
        print("==> fit3omega: saved sample config\n%s" % save_name)

    def _apply_sample_sliders(self, _) -> None:
        """change sample parameters based on slider values"""
        for label, slider in self.sample_sliders.items():
            param_name, layer_name = label.split('.')
            self.sample.modify_layer(layer_name, param_name, slider.val)
        self._update_graph()

    def _apply_heater_sliders(self, _) -> None:
        """change heater parameters based on slider values"""
        for label, slider in self.heater_sliders.items():
            self.sample.modify_heater(label, slider.val)
        self._update_graph()
        self._update_error_shading()

    def _reset_sliders(self, _) -> None:
        """reset sliders and calculation parameters to starting positions"""
        for label in self.sample_sliders:
            self.sample_sliders[label].reset()
        for label in self.heater_sliders:
            self.heater_sliders[label].reset()
        self.sample = self._original_sample.copy()

    def _toggle_xscale(self, _) -> None:
        """toggle the scale of the x-axis on the graph between log and linear"""
        if self._x_scale == "log(x)":
            self.ax.set_xscale("linear")
            self._x_scale = "x"
        else:
            self.ax.set_xscale("log")
            self._x_scale = "log(x)"

        self.buttons["x-scale"].label.set_text(self._x_scale)

    def _get_slider_dims(self) -> List[float]:
        """returns a location for the slider on the canvas"""
        dims = []
        n = self._n_sliders
        for x, d in zip(self.slider_start_dims, self.slider_delta):
            dims.append(x + n * d)
        return dims

    def _get_error_color(self, err) -> Tuple[float, float, float]:
        """returns a color value for the error (greener = better)"""
        t = self.error_green_thresh
        if err > 1.0:
            return 1.0, 0.0, 0.0
        elif t < err < 1.0:
            return (err - t), 0.0, 0.0
        elif 0.0 <= err <= t:
            return 0.0, (t - err) / t, 0.0
        else:
            return 0.0, 0.0, 0.0

    def _get_fitline_data(self) -> Tuple[np.ndarray, np.ndarray, float]:
        """returns current T_real, T_imag, and the error vs. measured data"""
        T_func_vals = self.fitted_T2
        err = self.objective_func(self.sample.x)
        return T_func_vals.real, T_func_vals.imag, err

    def _update_graph(self):
        fit = self._get_fitline_data()
        self.ax.lines[0].set_ydata(self.T2.x)
        self.ax.lines[1].set_ydata(self.T2.y)
        self.ax.lines[2].set_ydata(fit[0])
        self.ax.lines[3].set_ydata(fit[1])
        self.ax.texts[0].set_text(self.error_fmt.format(fit[2]))
        self.ax.texts[0].set_color(self._get_error_color(fit[2]))
        self.fig.canvas.draw_idle()

    def _update_error_shading(self):
        # re-shade to show error ranges
        self.ax.collections = []
        X = self.data.omegas / 2. / PI
        cx = 'blue'
        cy = 'red'
        self.ax.fill_between(X,
                             self.T2.x * (1. - self.T2.xerr),
                             self.T2.x * (1. + self.T2.xerr),
                             alpha=0.2, color=cx)
        self.ax.fill_between(X,
                             self.T2.y * (1. - self.T2.yerr),
                             self.T2.y * (1. + self.T2.yerr),
                             alpha=0.2, color=cy)
        self.fig.canvas.draw_idle()
