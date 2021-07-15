import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.widgets import Slider, Button
from math import pi as PI

from fit3omega.varsample import VarSample

"""
Standard plots and visualization elements.
"""


def _set_mpl_defaults():
    mpl.rc('xtick', direction='in')
    mpl.rc('ytick', direction='in')
    mpl.rc('axes', labelsize=15)
    mpl.rc('lines', marker='o')
    mpl.rc('lines', markerfacecolor='white')
    mpl.rc('lines', markersize=4)
    mpl.rc('errorbar', capsize=2)
    mpl.rc('legend', fontsize=13)


def plot_measured_data(m, show: bool = False) -> plt.Figure:
    """basic plot of sample voltages, shunt current, and temperature rise"""
    _set_mpl_defaults()
    fig = plt.figure(tight_layout=True, figsize=(10, 8))

    ax_V = fig.add_subplot(221)
    ax_Ish = fig.add_subplot(222)
    ax_V3 = fig.add_subplot(223)
    ax_T2 = fig.add_subplot(224)

    ax_V.set_xscale('log')
    ax_Ish.set_xscale('log')
    ax_V3.set_xscale('log')
    ax_T2.set_xscale('log')

    ax_V.set_ylabel(r"Sample V$_{1\omega}$")
    ax_V.set_xlabel(r"Source Frequency [Hz]")

    ax_Ish.set_ylabel(r"Shunt Current")
    ax_Ish.set_xlabel(r"Source Frequency [Hz]")

    ax_V3.set_ylabel(r"Sample V$_{3\omega}$")
    ax_V3.set_xlabel(r"Source Frequency [Hz]")

    ax_T2.set_ylabel(r"Sample $\widebar{T}_{2\omega}$")
    ax_T2.set_xlabel(r"Source Frequency [Hz]")

    cx = 'blue'
    cy = 'red'

    X = m.omegas / 2 / PI

    ax_V.errorbar(X, m.V.x, m.V.xerr * m.V.x, color=cx, elinewidth=.5)
    ax_V.errorbar(X, m.V.y, m.V.yerr * m.V.y, color=cy, elinewidth=.5)
    ax_V.grid(which="both")

    ax_Ish.errorbar(X, m.Ish.x, m.Ish.xerr * m.Ish.x, color=cx, label='X',
                    elinewidth=.5)
    ax_Ish.errorbar(X, m.Ish.y, m.Ish.yerr * m.Ish.y, color=cy, label='Y',
                    elinewidth=.5)
    ax_Ish.legend(frameon=False, fontsize=15)
    ax_Ish.grid(which="both")

    ax_V3.errorbar(X, m.V3.x, m.V3.xerr * m.V3.x, color=cx, elinewidth=.5)
    ax_V3.errorbar(X, m.V3.y, m.V3.yerr * m.V3.y, color=cy, elinewidth=.5)
    ax_V3.grid(which="both")

    ax_T2.errorbar(X, m.T2.x, m.T2.xerr * m.T2.x, color=cx, markerfacecolor=cx, elinewidth=.5)
    ax_T2.errorbar(X, m.T2.y, m.T2.yerr * m.T2.y, color=cy, markerfacecolor=cy, elinewidth=.5)
    ax_T2.grid(which="both")

    if show:
        plt.show()
    return fig


def plot_compare_measured_data(m1, m2, show: bool = False) -> plt.Figure:
    """basic plot of sample voltages, shunt current, and temperature rise"""
    _set_mpl_defaults()
    mpl.rc('lines', markersize=0)

    fig = plt.figure(tight_layout=True, figsize=(10, 8))

    ax_V = fig.add_subplot(221)
    ax_Ish = fig.add_subplot(222)
    ax_V3 = fig.add_subplot(223)
    ax_T2 = fig.add_subplot(224)

    ax_V.set_xscale('log')
    ax_Ish.set_xscale('log')
    ax_V3.set_xscale('log')
    ax_T2.set_xscale('log')

    ax_V.set_ylabel(r"Sample V$_{1\omega}$")
    ax_V.set_xlabel(r"Source Frequency [Hz]")

    ax_Ish.set_ylabel(r"Shunt Current")
    ax_Ish.set_xlabel(r"Source Frequency [Hz]")

    ax_V3.set_ylabel(r"Sample V$_{3\omega}$")
    ax_V3.set_xlabel(r"Source Frequency [Hz]")

    ax_T2.set_ylabel(r"Sample $\widebar{T}_{2\omega}$")
    ax_T2.set_xlabel(r"Source Frequency [Hz]")

    X2 = m1.omegas / 2 / PI

    cx1 = 'blue'
    cy1 = 'red'

    ax_V.plot(X2, m1.V.x, color=cx1)
    ax_V.plot(X2, m1.V.y, color=cy1)
    ax_V.text(0.1, 0.5, "dark: %s" % m1.data.data_file, fontsize=7,
              transform=ax_V.transAxes)
    ax_V.text(0.1, 0.45, "light: %s" % m2.data.data_file, fontsize=7,
              transform=ax_V.transAxes)
    ax_V.grid(which="both")

    ax_Ish.plot(X2, m1.Ish.x, color=cx1, label='X')
    ax_Ish.plot(X2, m1.Ish.y, color=cy1, label='Y')
    ax_Ish.legend(frameon=False, fontsize=15)
    ax_Ish.grid(which="both")

    ax_V3.plot(X2, m1.V3.x, color=cx1)
    ax_V3.plot(X2, m1.V3.y, color=cy1)
    ax_V3.grid(which="both")

    ax_T2.plot(X2, m1.T2.x, color=cx1)
    ax_T2.plot(X2, m1.T2.y, color=cy1)
    ax_T2.grid(which="both")

    X2 = m2.omegas / 2 / PI

    cx2 = 'turquoise'
    cy2 = 'salmon'

    ax_V.plot(X2, m2.V.x, color=cx2)
    ax_V.plot(X2, m2.V.y, color=cy2)
    ax_V.grid(which="both")

    ax_Ish.plot(X2, m2.Ish.x, color=cx2)
    ax_Ish.plot(X2, m2.Ish.y, color=cy2)
    ax_Ish.legend(frameon=False, fontsize=15)
    ax_Ish.grid(which="both")

    ax_V3.plot(X2, m2.V3.x, color=cx2)
    ax_V3.plot(X2, m2.V3.y, color=cy2)
    ax_V3.grid(which="both")

    ax_T2.plot(X2, m2.T2.x, color=cx2)
    ax_T2.plot(X2, m2.T2.y, color=cy2)
    ax_T2.grid(which="both")

    if show:
        plt.show()
    return fig


def plot_fitted_T2(m, show: bool = False) -> tuple:
    """plot fitted curves against measured T2 components"""
    _set_mpl_defaults()
    fig, ax = plt.subplots(tight_layout=True)

    ax.set_xlabel(r"Source Frequency [Hz]")
    ax.set_ylabel(r"$T_{2\omega,rms}$ [K]")

    cs = ["blue", "red"]
    ls = ["X", "Y"]

    X = m.omegas / 2. / PI
    ax.errorbar(X, m.T2.x, m.T2.xerr * m.T2.x,
                linewidth=0, elinewidth=.5, color=cs[0], label=ls[0])
    ax.errorbar(X, m.T2.y, m.T2.yerr * m.T2.y,
                linewidth=0, elinewidth=.5, color=cs[1], label=ls[1])

    ax.plot(X, m.T2_fit.x, markersize=0, color=cs[0])
    ax.plot(X, m.T2_fit.y, markersize=0, color=cs[1])

    ax.legend(frameon=False)

    ax.set_xscale('log')
    if show:
        plt.show()
    return fig, ax


def plot_fitted_Z2(m, show: bool = False) -> tuple:
    """plot fitted curves against measured Z2 components"""
    # TODO: this won't work
    _set_mpl_defaults()
    fig, ax = plt.subplots(tight_layout=True)

    ax.set_xlabel(r"Source Frequency [Hz]")
    ax.set_ylabel(r"$Z_{2\omega,rms}$ [K]")

    cs = ["purple", "crimson"]
    ls = ["X", "Y"]

    X = m.omegas / 2 / PI
    ax.errorbar(X, m.Z2.x, m.Z2.xerr * m.Z2.x,
                linewidth=0, elinewidth=.5, color=cs[0], label=ls[0])
    ax.errorbar(X, m.Z2.y, m.Z2.yerr * m.Z2.y,
                linewidth=0, elinewidth=.5, color=cs[1], label=ls[1])

    ax.plot(X, m.Z2_fit.x, markersize=0, color=cs[0])
    ax.plot(X, m.Z2_fit.y, markersize=0, color=cs[1])

    ax.legend(frameon=False)

    ax.set_xscale('log')
    if show:
        plt.show()
    return fig, ax


def plot_fitted_T2_linear(m, show: bool = False) -> tuple:
    """plot fitted curves against measured T2 components"""
    _set_mpl_defaults()
    fig, ax = plt.subplots(tight_layout=True)

    ax.set_xlabel(r"Source Frequency [Hz]")
    ax.set_ylabel(r"$T_{2\omega,rms}$ [K]")

    cs = ["blue", "red", "black"]
    ls = ["X", "Y", "R"]

    X = m.omegas / 2. / PI
    ax.errorbar(X, m.T2.x, m.T2.xerr * m.T2.x,
                linewidth=0, elinewidth=.5, color=cs[0], label=ls[0])
    ax.errorbar(X, m.T2.y, m.T2.yerr * m.T2.y,
                linewidth=0, elinewidth=.5, color=cs[1], label=ls[1])

    X_linear = m.omegas_linear / 2. / PI
    ax.plot(X_linear, m.T2_fit.x, markersize=0, linewidth=2, color="cyan")
    ax.plot(X_linear, m.T2_fit.y, markersize=0, linewidth=2, color="fuchsia")

    ax.legend(frameon=False)

    ax.set_xscale('log')
    if show:
        plt.show()
    return fig, ax


def plot_diagnostics(m, show: bool = False) -> plt.Figure:
    _set_mpl_defaults()
    """plot V3/I^3 and power from measured data"""
    from numpy import sqrt
    fig = plt.figure(figsize=(10, 5))

    ax_power = fig.add_subplot(121)
    ax_VI = fig.add_subplot(122)

    ax_power.set_xscale("log")
    ax_VI.set_xscale("log")

    fig.text(0.5, 0.05, "Source Frequency [Hz]",
             ha="center", va="top", fontsize=15)

    ax_power.set_ylabel("Real Power [mW]")
    ax_VI.set_ylabel(r"$V_{3\omega,x} (I_{1\omega})^{-3}$")

    ax_power.grid("both")
    ax_VI.grid("both")

    X = m.omegas / 2. / PI

    power = m.power.x * 1e3  # mW
    dpower = m.power.abserr()

    ax_power.errorbar(X, power, dpower, elinewidth=.5, color='k')
    min_power = min(power)
    max_power = max(power)
    power_range = max_power - min_power
    ax_power.set_ylim(min(power) - power_range, max(power) + power_range)

    V3x = m.V3.x
    dV3x = V3x * m.V3.xerr
    I_cubed = m.Ish.x**3
    dI_cubed = sqrt(3) * m.Ish.xerr * I_cubed

    VI = V3x / I_cubed
    dVI = sqrt(dV3x**2 + dI_cubed**2)

    ax_VI.errorbar(X, VI, dVI, elinewidth=.5, color='coral')

    if show:
        plt.show()
    return fig


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

    def __init__(self, m, frac: float = 0.99, niter: int = None,
                 enable_heater_params: bool = False):
        mpl.rc("font", family="monospace")
        self.frac = frac
        self.niter = niter
        self.enable_heater_params = enable_heater_params
        m.sample = VarSample.from_sample(m.sample)

        self.model = m
        self.model.set_refresh(True)
        self.model_type = type(self.model)
        self.fig, self.ax = plt.subplots(figsize=(8, 8))
        self.fig.subplots_adjust(bottom=.5, top=0.95)
        self.param_sliders = {}
        self.sample_sliders = {}
        self.buttons = {}

        self._x_scale = "log(x)"
        self._n_sliders = 0

    def plot_initial_state(self):
        _set_mpl_defaults()

        self.ax.set_xscale("log")
        self.ax.set_ylabel(r"Sample $\widebar{T}_{2\omega}$")
        self.ax.set_xlabel(r"$2\pi\omega$ [Hz]")

        X = self.model.omegas / 2. / PI
        cx = 'blue'
        cy = 'red'
        fcx = u'#0B554C'
        fcy = u'#550B14'

        # measured data
        self.ax.plot(X, self.model.T2.x, color=cx, **self.meas_plot_kw)
        self.ax.plot(X, self.model.T2.y, color=cy, **self.meas_plot_kw)

        # fitted data (initial values)
        fit = self.get_fitline_data()
        self.ax.plot(X, fit[0], color=fcx, **self.fit_plot_kw)
        self.ax.plot(X, fit[1], color=fcy, **self.fit_plot_kw)

        # prepare the fit parameters sliders
        for i, layer in enumerate(self.model.sample.layers):
            d = layer.as_dict()
            for k, v in d.items():
                k_ = k + 's'
                if type(v) is str and v.endswith('*') and k_ in self.model.FIT_PARAMS:
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
                dRdT=self.model.sample.heater.dRdT,
                # width=self.model.sample.heater.width,
                # length=self.model.sample.heater.length,
                Cv=self.model.sample.heater.Cv,
                height=self.model.sample.heater.height,
                Rc=self.model.sample.heater.Rc
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

    def get_fitline_data(self):
        return self.model.get_current_T2()

    def _apply_sliders(self, _):
        for label, slider in self.param_sliders.items():
            layer_name, attr_name = label.split('.')
            self.model.sample.param_modify(layer_name, attr_name + 's', slider.val)
        for label, slider in self.sample_sliders.items():
            self.model.sample.heater.modify(label, slider.val)
        self._update_graph()

    def _update_graph(self):
        fit = self.get_fitline_data()
        self.ax.lines[0].set_ydata(self.model.T2.x)
        self.ax.lines[1].set_ydata(self.model.T2.y)
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
        self.model.sample.reset_params()
        self.model.heater.reset()

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
        niter = self._get_niter_estimate() if self.niter is None else self.niter
        err = self.get_fitline_data()[-1]
        self.model.fit(niter=niter, min_err=err)
        for attr_name, values in self.model.fitted_kwargs.items():
            self.model.sample.param_modify(None, attr_name, values)
        self._update_graph()
        print("--> Result:\n", self.model.result, "\n")

        for label, slider in self.param_sliders.items():
            layer_name, attr_name = label.split('.')
            for i, layer in enumerate(self.model.sample.layers):
                if layer.name == layer_name:
                    slider.set_val(self.model.fitted_kwargs[attr_name + "s"][i])
                    break

    def _get_niter_estimate(self):
        N = len(self.model.omegas)
        if 0 < N < 10:
            return 100
        elif 10 <= N <= 30:
            return 50
        else:
            return 30
