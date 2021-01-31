import matplotlib.pyplot as plt


def plot_data(m, show: bool = False) -> plt.Figure:
    """basic plot of sample voltages, shunt current, and temperature rise"""
    import matplotlib as mpl
    mpl.rc('xtick', direction='in')
    mpl.rc('ytick', direction='in')
    mpl.rc('axes', labelsize=15)
    mpl.rc('lines', marker='o')
    mpl.rc('lines', markerfacecolor='white')
    mpl.rc('lines', markersize=5)
    mpl.rc('errorbar', capsize=(0 if m.data.no_error else 3))

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
    ax_V.set_xlabel(r"$\omega$ [Hz]")

    ax_Ish.set_ylabel(r"Shunt Current")
    ax_Ish.set_xlabel(r"$\omega$ [Hz]")

    ax_V3.set_ylabel(r"Sample V$_{3\omega}$")
    ax_V3.set_xlabel(r"$\omega$ [Hz]")

    ax_T2.set_ylabel(r"Sample T$_{2\omega}$")
    ax_T2.set_xlabel(r"$\omega$ [Hz]")

    cx = 'blue'
    cy = 'red'
    cz = 'black'

    ax_V.errorbar(m.omegas, m.V.x, yerr=m.V.xerr * m.V.x, color=cx, elinewidth=.8)
    ax_V.errorbar(m.omegas, m.V.y, yerr=m.V.yerr * m.V.y, color=cy, elinewidth=.8)
    ax_V.errorbar(m.omegas, m.V.norm(), yerr=m.V.abserr(), color=cz, elinewidth=.8)
    ax_V.grid(which="both")

    ax_Ish.errorbar(m.omegas, m.Ish.x, yerr=m.Ish.xerr * m.Ish.x, color=cx, label='X', elinewidth=.8)
    ax_Ish.errorbar(m.omegas, m.Ish.y, yerr=m.Ish.yerr * m.Ish.y, color=cy, label='Y', elinewidth=.8)
    ax_Ish.errorbar(m.omegas, m.Ish.norm(), yerr=m.Ish.abserr(), color=cz, label='R', elinewidth=.8)
    ax_Ish.legend(frameon=False, fontsize=15)
    ax_Ish.grid(which="both")

    ax_V3.errorbar(m.omegas, m.V3.x, yerr=m.V3.xerr * m.V3.x, color=cx, elinewidth=.8)
    ax_V3.errorbar(m.omegas, m.V3.y, yerr=m.V3.yerr * m.V3.y, color=cy, elinewidth=.8)
    ax_V3.errorbar(m.omegas, m.V3.norm(), yerr=m.V3.abserr(), color=cz, elinewidth=.8)
    ax_V3.grid(which="both")

    ax_T2.errorbar(m.omegas, m.T2.x, yerr=m.T2.xerr * m.T2.x, color=cx, elinewidth=.8)
    ax_T2.errorbar(m.omegas, m.T2.y, yerr=m.T2.yerr * m.T2.y, color=cy, elinewidth=.8)
    ax_T2.errorbar(m.omegas, m.T2.norm(), yerr=m.T2.abserr(), color=cz, elinewidth=.8)
    ax_T2.grid(which="both")

    if show:
        plt.show()
    return fig
