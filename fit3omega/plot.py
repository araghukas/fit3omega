import matplotlib.pyplot as plt
import matplotlib as mpl
from math import pi as PI

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

    ax_T2.errorbar(X, m.T2.x, m.T2.xerr * m.T2.x, color=cx, elinewidth=.5)
    ax_T2.errorbar(X, m.T2.y, m.T2.yerr * m.T2.y, color=cy, elinewidth=.5)
    ax_T2.errorbar(X, m.T2.norm(), m.T2.abserr(), color=cz, elinewidth=.5)
    ax_T2.grid(which="both")

    if show:
        plt.show()
    return fig


def plot_fitted_T2(fg, show: bool = False) -> plt.Figure:
    """plot fitted curves against measured T2 components"""
    fig, ax = plt.subplots(tight_layout=True)

    ax.set_xlabel(r"Source Frequency [Hz]")
    ax.set_ylabel(r"$T_{2\omega,rms}$ [K]")

    cs = ["red", "blue", "black"]
    ls = ["cos", "sin", "abs"]

    X = fg.omegas / 2 / PI
    ax.errorbar(X, fg.T2.x, fg.T2.xerr * fg.T2.x, linewidth=0, elinewidth=.5,
                color=cs[0], label=ls[0])
    ax.errorbar(X, fg.T2.y, fg.T2.yerr * fg.T2.y, linewidth=0, elinewidth=.5,
                color=cs[1], label=ls[1])
    ax.errorbar(X, fg.T2.norm(), fg.T2.relerr() * fg.T2.norm(), linewidth=0, elinewidth=1,
                color=cs[2], label=ls[2], markersize=mpl.rcParams["lines.markersize"] * 1.5)

    ax.plot(X, fg.T2_fit.x, markersize=0, color=cs[0])
    ax.plot(X, fg.T2_fit.y, markersize=0, color=cs[1])
    ax.plot(X, fg.T2_fit.norm(), markersize=0, color=cs[2],
            linewidth=mpl.rcParams["lines.linewidth"] + 1)

    ax.legend(frameon=False)

    ax.set_xscale('log')
    if show:
        plt.show()
    return fig
