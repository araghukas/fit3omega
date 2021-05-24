import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.widgets import Slider, Button
from math import pi as PI


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

    ax_T2.set_ylabel(r"Sample T$_{2\omega}$")
    ax_T2.set_xlabel(r"Source Frequency [Hz]")

    cx = 'blue'
    cy = 'red'
    cz = 'black'

    X = m.omegas / 2 / PI

    ax_V.errorbar(X, m.V.x, m.V.xerr * m.V.x, color=cx, elinewidth=.5)
    ax_V.errorbar(X, m.V.y, m.V.yerr * m.V.y, color=cy, elinewidth=.5)
    ax_V.errorbar(X, m.V.norm(), m.V.abserr(), color=cz, elinewidth=.5)
    ax_V.grid(which="both")

    ax_Ish.errorbar(X, m.Ish.x, m.Ish.xerr * m.Ish.x, color=cx, label='X',
                    elinewidth=.5)
    ax_Ish.errorbar(X, m.Ish.y, m.Ish.yerr * m.Ish.y, color=cy, label='Y',
                    elinewidth=.5)
    ax_Ish.errorbar(X, m.Ish.norm(), m.Ish.abserr(), color=cz, label='R', elinewidth=.5)
    ax_Ish.legend(frameon=False, fontsize=15)
    ax_Ish.grid(which="both")

    ax_V3.errorbar(X, m.V3.x, m.V3.xerr * m.V3.x, color=cx, elinewidth=.5)
    ax_V3.errorbar(X, m.V3.y, m.V3.yerr * m.V3.y, color=cy, elinewidth=.5)
    ax_V3.errorbar(X, m.V3.norm(), m.V3.abserr(), color=cz, elinewidth=.5)
    ax_V3.grid(which="both")

    ax_T2.errorbar(X, m.T2.x, m.T2.xerr * m.T2.x, color=cx, markerfacecolor=cx, elinewidth=.5)
    ax_T2.errorbar(X, m.T2.y, m.T2.yerr * m.T2.y, color=cy, markerfacecolor=cy, elinewidth=.5)
    ax_T2.errorbar(X, m.T2.norm(), m.T2.abserr(), color=cz, markerfacecolor=cz, elinewidth=.5)
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

    ax_T2.set_ylabel(r"Sample T$_{2\omega}$")
    ax_T2.set_xlabel(r"Source Frequency [Hz]")

    X2 = m1.omegas / 2 / PI

    cx1 = 'blue'
    cy1 = 'red'
    cz1 = 'black'

    ax_V.plot(X2, m1.V.x, color=cx1)
    ax_V.plot(X2, m1.V.y, color=cy1)
    ax_V.plot(X2, m1.V.norm(), color=cz1)
    ax_V.text(0.1, 0.5, "dark: %s" % m1.data.data_file, fontsize=7,
              transform=ax_V.transAxes)
    ax_V.text(0.1, 0.45, "light: %s" % m2.data.data_file, fontsize=7,
              transform=ax_V.transAxes)
    ax_V.grid(which="both")

    ax_Ish.plot(X2, m1.Ish.x, color=cx1, label='X')
    ax_Ish.plot(X2, m1.Ish.y, color=cy1, label='Y')
    ax_Ish.plot(X2, m1.Ish.norm(), color=cz1, label='R')
    ax_Ish.legend(frameon=False, fontsize=15)
    ax_Ish.grid(which="both")

    ax_V3.plot(X2, m1.V3.x, color=cx1)
    ax_V3.plot(X2, m1.V3.y, color=cy1)
    ax_V3.plot(X2, m1.V3.norm(), color=cz1)
    ax_V3.grid(which="both")

    ax_T2.plot(X2, m1.T2.x, color=cx1)
    ax_T2.plot(X2, m1.T2.y, color=cy1)
    ax_T2.plot(X2, m1.T2.norm(), color=cz1)
    ax_T2.grid(which="both")

    X2 = m2.omegas / 2 / PI

    cx2 = 'turquoise'
    cy2 = 'salmon'
    cz2 = 'grey'

    ax_V.plot(X2, m2.V.x, color=cx2)
    ax_V.plot(X2, m2.V.y, color=cy2)
    ax_V.plot(X2, m2.V.norm(), color=cz2)
    ax_V.grid(which="both")

    ax_Ish.plot(X2, m2.Ish.x, color=cx2)
    ax_Ish.plot(X2, m2.Ish.y, color=cy2)
    ax_Ish.plot(X2, m2.Ish.norm(), color=cz2)
    ax_Ish.legend(frameon=False, fontsize=15)
    ax_Ish.grid(which="both")

    ax_V3.plot(X2, m2.V3.x, color=cx2)
    ax_V3.plot(X2, m2.V3.y, color=cy2)
    ax_V3.plot(X2, m2.V3.norm(), color=cz2)
    ax_V3.grid(which="both")

    ax_T2.plot(X2, m2.T2.x, color=cx2)
    ax_T2.plot(X2, m2.T2.y, color=cy2)
    ax_T2.plot(X2, m2.T2.norm(), color=cz2)
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

    cs = ["blue", "red", "black"]
    ls = ["X", "Y", "R"]

    X = m.omegas / 2 / PI
    ax.errorbar(X, m.T2.x, m.T2.xerr * m.T2.x,
                linewidth=0, elinewidth=.5, color=cs[0], label=ls[0])
    ax.errorbar(X, m.T2.y, m.T2.yerr * m.T2.y,
                linewidth=0, elinewidth=.5, color=cs[1], label=ls[1])
    ax.errorbar(X, m.T2.norm(), m.T2.relerr() * m.T2.norm(),
                linewidth=0, elinewidth=1, color=cs[2], label=ls[2])

    ax.plot(X, m.T2_fit.x, markersize=0, color=cs[0])
    ax.plot(X, m.T2_fit.y, markersize=0, color=cs[1])
    ax.plot(X, m.T2_fit.norm(), markersize=0, color=cs[2])

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

    ax_power.set_ylabel("Power [mW]")
    ax_VI.set_ylabel(r"$V_{3\omega,x} (I_{1\omega})^{-3}$")

    ax_power.grid("both")
    ax_VI.grid("both")

    X = m.omegas / 2. / PI

    power = m.power.norm() * 1e3  # mW
    dpower = m.power.abserr()

    ax_power.errorbar(X, power, dpower, elinewidth=.5, color='k')
    min_power = min(power)
    max_power = max(power)
    power_range = max_power - min_power
    ax_power.set_ylim(min(power) - power_range, max(power) + power_range)

    V3x = m.V3.x
    dV3x = V3x * m.V3.xerr
    I_cubed = m.Ish.norm()**3
    dI_cubed = sqrt(3) * m.Ish.abserr()

    VI = V3x / I_cubed
    dVI = sqrt(dV3x**2 + dI_cubed**2)

    ax_VI.errorbar(X, VI, dVI, elinewidth=.5, color='coral')

    if show:
        plt.show()
    return fig


class SliderPlot:
    # TODO: modify non-array params
    # TODO: some slider customization

    meas_plot_kw = dict(
        markersize=6,
        marker='x',
        markerfacecolor='w',
        zorder=0,
        linestyle='None'
    )

    fit_plot_kw = dict(
        linestyle='--',
        linewidth=1,
        markersize=0
    )

    # left, bottom, width, height; dimensions normalized (0,1)
    button_dims = [0.6, 0.1, 0.05, 0.05]
    button_hovercolor = "0.98"

    text_box_dims = [0.6, 0.15, 0.05, 0.05]

    slider_start_dims = [0.6, 0.9, 0.3, 0.03]

    slider_delta = [0.0, -0.05, 0.0, 0.0]

    error_fmt = "error: {:<10,.6e}"

    def __init__(self, m):
        mpl.rc("font", family="monospace")
        self.frac = 0.99
        self._x_scale = "log(x)"

        m.sample = m.sample.as_var_sample()

        self.m = m
        self.fig, self.ax = plt.subplots(figsize=(12, 5))
        self.fig.subplots_adjust(right=.5)
        self.sliders = {}
        self.buttons = {}

    def get_fitline_data(self):
        heights = self.m.sample.heights
        kys = self.m.sample.kys
        ratio_xys = self.m.sample.ratio_xys
        Cvs = self.m.sample.Cvs

        T2_complex = self.m.T2_func(heights, kys, ratio_xys, Cvs)
        T2_x = T2_complex.real
        T2_y = T2_complex.imag
        T2_norm = abs(T2_complex)
        err = sum(abs(self.m.T2.x - T2_complex.real))
        err += sum(abs(self.m.T2.y - T2_complex.imag))
        return T2_x, T2_y, T2_norm, err / len(T2_x)

    def plot_initial_state(self):
        _set_mpl_defaults()

        self.ax.set_xscale("log")
        self.ax.set_ylabel(r"Sample T$_{2\omega}$")
        self.ax.set_xlabel(r"$2\pi\omega$ [Hz]")

        X = self.m.omegas / 2. / PI
        cx = 'blue'
        cy = 'red'
        cz = 'black'

        # measured data
        self.ax.plot(X, self.m.T2.x, color=cx, **self.meas_plot_kw)
        self.ax.plot(X, self.m.T2.y, color=cy, **self.meas_plot_kw)
        self.ax.plot(X, self.m.T2.norm(), color=cz, **self.meas_plot_kw)

        # fitted data (initial values)
        fit = self.get_fitline_data()
        self.ax.plot(X, fit[0], color=cx, **self.fit_plot_kw)
        self.ax.plot(X, fit[1], color=cy, **self.fit_plot_kw)
        self.ax.plot(X, fit[2], color=cz, **self.fit_plot_kw)
        self.ax.text(0.7, 1.05, self.error_fmt.format(fit[3]), transform=self.ax.transAxes,
                     fontsize=8, color=self._get_error_color(fit[3]), fontweight="bold")

        # prepare the sliders
        for i, layer in enumerate(self.m.sample.layers):
            d = layer.as_dict()
            for k, v in d.items():
                if type(v) is str and v.endswith('*'):
                    name = layer.name
                    label = name + '.' + k
                    guess_val = float(v.rstrip('*'))
                    axes = plt.axes(self._get_slider_dims())
                    self.sliders[label] = Slider(
                        ax=axes,
                        label=label.replace('.', ' '),
                        valmin=(1 - self.frac) * guess_val,
                        valmax=(1 + self.frac) * guess_val,
                        valinit=guess_val,
                        valfmt="%.2e"
                    )
                    self.sliders[label].on_changed(self._update_sliders)

        # sloppy reset button creation
        reset_button_dims = self._get_slider_dims()
        reset_button_dims[2:] = [0.05, 0.05]
        reset_button_dims[1] -= 0.025
        self.buttons["Reset"] = Button(plt.axes(reset_button_dims), "Reset",
                                       hovercolor=self.button_hovercolor)
        self.buttons["Reset"].on_clicked(self._reset_sliders)

        # sloppy x-scale toggle creation
        xscale_button_dims = reset_button_dims
        xscale_button_dims[0] += 0.055
        self.buttons["x-scale"] = Button(plt.axes(xscale_button_dims), self._x_scale,
                                         hovercolor=self.button_hovercolor)
        self.buttons["x-scale"].on_clicked(self._toggle_xscale)

    def _update_sliders(self, _):
        for label, slider in self.sliders.items():
            layer_name, attr_name = label.split('.')
            self.m.sample.param_modify(layer_name, attr_name + 's', slider.val)
        fit = self.get_fitline_data()
        self.ax.lines[3].set_ydata(fit[0])
        self.ax.lines[4].set_ydata(fit[1])
        self.ax.lines[5].set_ydata(fit[2])
        self.ax.texts[0].set_text(self.error_fmt.format(fit[3]))
        self.ax.texts[0].set_color(self._get_error_color(fit[3]))
        self.fig.canvas.draw_idle()

    def _reset_sliders(self, _):
        for label in self.sliders:
            self.sliders[label].reset()
        self.m.sample.reset_params()

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
        n = len(self.sliders)
        for x, d in zip(self.slider_start_dims, self.slider_delta):
            dims.append(x + n * d)
        return dims

    @staticmethod
    def _get_error_color(err):
        if err > 1.0:
            return 1.0, 0.0, 0.0
        elif 0.1 < err < 1.0:
            return (err - 0.1), 0.0, 0.0
        elif 0.0 <= err <= 0.1:
            return 0.0, (0.1 - err) / 0.1, 0.0
        else:
            return 0.0, 0.0, 0.0

