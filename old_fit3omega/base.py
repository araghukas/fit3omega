"""module defining base class for object performing data fits"""
from abc import ABC, abstractmethod
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import pi as PI
from typing import Tuple
import numpy as np

from old_fit3omega.data import ACReading
from old_fit3omega.model import Model
from old_fit3omega.result import FitResult


class BaseFitter(ABC):
    """base class for data fitting objects"""

    @property
    @abstractmethod
    def T2_fit(self) -> ACReading:
        """discrete T2 values of the fit line at the measured frequencies"""
        raise NotImplementedError

    @property
    @abstractmethod
    def result(self) -> FitResult:
        """resulting sample parameters extracted from fit line"""
        raise NotImplementedError

    @property
    @abstractmethod
    def error(self) -> float:
        """MSE of the fit line versus measured data"""
        raise NotImplementedError

    @abstractmethod
    def T2_func(self, *args, **kwargs) -> np.ndarray:
        """synthetic average T2 at each ω, from sample parameters"""
        raise NotImplementedError

    @abstractmethod
    def fit(self, *args, **kwargs) -> FitResult:
        """calculate the fit line and extract sample parameters"""
        raise NotImplementedError

    def __init__(self, m: Model):
        self.model = m
        self._result = None
        self._error = None

    @property
    def model(self) -> Model:
        """the data model that transfers voltages to temperatures"""
        return self._model

    @model.setter
    def model(self, m):
        if type(m) is not Model:
            raise ValueError("model must be a old_fit3omega.model.Model instance.")
        self._model = m

    def set_data_limits(self, start: int, end: int) -> None:
        """truncate the measured data by taking an interior set of points"""
        self.model.set_data_limits(start, end)

    def drop_row(self, row_index: int) -> None:
        """delete the data row at the given index"""
        self.model.data.drop_row(row_index)

    # PLOTTING METHODS -----------------------------------------------------------------------------
    def plot_measured_data(self, show: bool = True) -> plt.Figure:
        """plot the experimental data with frequency on the x-axis"""
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
        V = self.model.V
        Ish = self.model.Ish
        V3 = self.model.V3
        T2 = self.model.T2
        # -----------------------

        X = self.model.omegas / (2.0 * PI)

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

    def plot_compare_measured_data(self,
                                   fitter2: 'BaseFitter',
                                   show: bool = True) -> plt.Figure:
        """
        Plot sample voltages, shunt current, and temperature rise
        between this and another fitter's data.
        """
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

        X2 = self.model.omegas / 2 / PI

        cx1 = 'blue'
        cy1 = 'red'

        ax_V.plot(X2, self.model.V.x, color=cx1)
        ax_V.plot(X2, self.model.V.y, color=cy1)
        ax_V.text(0.1, 0.5, "dark: %s" % self.model.data.data_file, fontsize=7,
                  transform=ax_V.transAxes)
        ax_V.text(0.1, 0.45, "light: %s" % fitter2.model.data.data_file, fontsize=7,
                  transform=ax_V.transAxes)
        ax_V.grid(which="both")

        ax_Ish.plot(X2, self.model.Ish.x, color=cx1, label='X')
        ax_Ish.plot(X2, self.model.Ish.y, color=cy1, label='Y')
        ax_Ish.legend(frameon=False, fontsize=15)
        ax_Ish.grid(which="both")

        ax_V3.plot(X2, self.model.V3.x, color=cx1)
        ax_V3.plot(X2, self.model.V3.y, color=cy1)
        ax_V3.grid(which="both")

        ax_T2.plot(X2, self.model.T2.x, color=cx1)
        ax_T2.plot(X2, self.model.T2.y, color=cy1)
        ax_T2.grid(which="both")

        X2 = fitter2.model.omegas / 2 / PI

        cx2 = 'turquoise'
        cy2 = 'salmon'

        ax_V.plot(X2, fitter2.model.V.x, color=cx2)
        ax_V.plot(X2, fitter2.model.V.y, color=cy2)
        ax_V.grid(which="both")

        ax_Ish.plot(X2, fitter2.model.Ish.x, color=cx2)
        ax_Ish.plot(X2, fitter2.model.Ish.y, color=cy2)
        ax_Ish.legend(frameon=False, fontsize=15)
        ax_Ish.grid(which="both")

        ax_V3.plot(X2, fitter2.model.V3.x, color=cx2)
        ax_V3.plot(X2, fitter2.model.V3.y, color=cy2)
        ax_V3.grid(which="both")

        ax_T2.plot(X2, fitter2.model.T2.x, color=cx2)
        ax_T2.plot(X2, fitter2.model.T2.y, color=cy2)
        ax_T2.grid(which="both")

        if show:
            plt.show()
        return fig

    def plot_fitted_T2(self, show: bool = True) -> Tuple[plt.Figure, plt.Axes]:
        """plot the fitted temperature curves over the measured temperature data"""
        _set_mpl_defaults()
        fig, ax = plt.subplots(tight_layout=True)

        ax.set_xlabel(r"Source Frequency [Hz]")
        ax.set_ylabel(r"$T_{2\omega,rms}$ [K]")

        cs = ["blue", "red"]
        ls = ["X", "Y"]

        X = self.model.omegas / 2. / PI
        ax.errorbar(X, self.model.T2.x, self.model.T2.xerr * self.model.T2.x,
                    linewidth=0, elinewidth=.5, color=cs[0], label=ls[0])
        ax.errorbar(X, self.model.T2.y, self.model.T2.yerr * self.model.T2.y,
                    linewidth=0, elinewidth=.5, color=cs[1], label=ls[1])

        ax.plot(X, self.T2_fit.x, markersize=0, color=cs[0])
        ax.plot(X, self.T2_fit.y, markersize=0, color=cs[1])
        ax.text(0.1, 0.5, self.result.summary().replace('\t', "    "),
                transform=ax.transAxes, fontsize=8)

        ax.legend(frameon=False)

        ax.set_xscale('log')
        if show:
            plt.show()
        return fig, ax


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
