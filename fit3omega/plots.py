"""a module of plotting function to visualize measured and fitted data"""
import matplotlib as mpl
import matplotlib.pyplot as plt

from math import pi as PI

from fit3omega.fit import Fit3omega


def _set_mpl_defaults():
    """set default parameters for plots"""
    mpl.rc('xtick', direction='in')
    mpl.rc('ytick', direction='in')
    mpl.rc('axes', labelsize=15)
    mpl.rc('lines', marker='o')
    mpl.rc('lines', markerfacecolor='white')
    mpl.rc('lines', markersize=4)
    mpl.rc('errorbar', capsize=2)
    mpl.rc('legend', fontsize=13)


def plot_measured_data(fitter: Fit3omega,
                       show: bool = True) -> plt.Figure:
    """plot the measured data with frequency on the x-axis"""
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

    # alias for brevity below
    V = fitter.data.V
    Ish = fitter.Ish
    V3 = fitter.data.V3
    T2 = fitter.T2
    # -----------------------

    X = fitter.data.omegas / (2.0 * PI)

    ax_V.errorbar(X, V.x, V.xerr * V.x, color=cx, elinewidth=.5)
    ax_V.errorbar(X, V.y, V.yerr * V.y, color=cy, elinewidth=.5)
    ax_V.grid(which="both")

    ax_Ish.errorbar(X, Ish.x, Ish.xerr * Ish.x, color=cx, label='X', elinewidth=.5)
    ax_Ish.errorbar(X, Ish.y, Ish.yerr * Ish.y, color=cy, label='Y', elinewidth=.5)
    ax_Ish.legend(frameon=False, fontsize=15)
    ax_Ish.grid(which="both")

    ax_V3.errorbar(X, V3.x, V3.xerr * V3.x, color=cx, elinewidth=.5)
    ax_V3.errorbar(X, V3.y, V3.yerr * V3.y, color=cy, elinewidth=.5)
    ax_V3.grid(which="both")

    ax_T2.errorbar(X, T2.x, T2.xerr * T2.x, color=cx, markerfacecolor=cx, elinewidth=.5)
    ax_T2.errorbar(X, T2.y, T2.yerr * T2.y, color=cy, markerfacecolor=cy, elinewidth=.5)
    ax_T2.grid(which="both")

    if show:
        plt.show()
    return fig


def plot_fitted_data(fitter: Fit3omega,
                     show: bool = True) -> plt.Figure:
    """plot the fitted temperature curves over the measured temperature data"""
    _set_mpl_defaults()
    fig, ax = plt.subplots(tight_layout=True)

    ax.set_xlabel(r"Source Frequency [Hz]")
    ax.set_ylabel(r"$T_{2\omega,rms}$ [K]")

    cs = ["blue", "red"]
    ls = ["X", "Y"]

    X = fitter.data.omegas / 2. / PI
    ax.errorbar(X, fitter.T2.x, fitter.T2.xerr * fitter.T2.x,
                linewidth=0, elinewidth=.5, color=cs[0], label=ls[0])
    ax.errorbar(X, fitter.T2.y, fitter.T2.yerr * fitter.T2.y,
                linewidth=0, elinewidth=.5, color=cs[1], label=ls[1])

    fitted_T2 = fitter.T2_function(*fitter.sample.substitute(fitter.result.x))
    ax.plot(X, fitted_T2.real, markersize=0, color=cs[0])
    ax.plot(X, fitted_T2.imag, markersize=0, color=cs[1])
    ax.text(0.1, 0.37, fitter.result.summary.replace('\t', "    ").replace("\n", "\n\n"),
            transform=ax.transAxes, fontsize=8)

    ax.legend(frameon=False)

    ax.set_xscale('log')
    if show:
        plt.show()
    return fig
