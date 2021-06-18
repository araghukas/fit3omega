import numpy as np

from .sample import *
from .data import Data, ACReading
import fit3omega.plot as plot

ROOT2 = np.sqrt(2)


class Model:
    """
    Thermal model for 2ω temperature rise, calculated from AC signals.

    All measurement data are RMS voltages, V_rms = Vx_rms + jVy_rms.

    The lock-in amplifier (SRS830) displays and outputs (X, Y) or (R, θ), where

        X^2 + Y^2 = R^2 = (V_signal_rms)^2
        X = R*cos(θ)
        Y = R*sin(θ)
    """

    def __init__(self, sample, data):
        if type(sample) is str:
            self.sample = Sample(sample)
        else:
            self.sample = sample

        if type(data) is str:
            self.data = Data(data)
        elif type(data) is Data:
            self.data = data
        else:
            raise ValueError("invalid `data` type; need str or fit3omega.data.Data")

        # if true, always recalculate data-derived values (below) when accessing props.
        self._refresh_dependents = False

        self._Ish = None
        self._T2 = None
        self._Z2 = None
        self._dT2 = None
        self._power = None

    def set_refresh(self, b: bool) -> None:
        if b is True:
            self._refresh_dependents = True
        elif b is False:
            self._refresh_dependents = False
        else:
            raise ValueError("argument of set_refresh is not a boolean")

    @property
    def shunt(self) -> Shunt:
        return self.sample.shunt

    @property
    def heater(self) -> Heater:
        return self.sample.heater

    @property
    def omegas(self) -> np.array:
        # array of sample points ω = 2πf
        return self.data.omegas

    @property
    def V(self) -> ACReading:
        # sample V[1ω] (measured)
        return self.data.V

    @property
    def V3(self) -> ACReading:
        # sample V[3ω] RMS
        return self.data.V3

    @property
    def Vsh(self) -> ACReading:
        # resistor V[1ω] RMS
        return self.data.Vsh

    @property
    def Ish(self) -> ACReading:
        """RMS series current at each ω"""
        if self._Ish is None or self._refresh_dependents:
            x = self.Vsh.x / self.shunt.R
            y = self.Vsh.y / self.shunt.R
            xerr = np.sqrt(self.Vsh.xerr**2 + self.shunt.err**2)
            yerr = np.sqrt(self.Vsh.yerr**2 + self.shunt.err**2)
            self._Ish = ACReading(x, y, xerr, yerr)
        return self._Ish

    @property
    def T2(self) -> ACReading:
        """average 2ω temperature oscillations at each ω"""
        if self._T2 is None or self._refresh_dependents:
            x = 2. * np.abs(self.V3.x) / (self.heater.dRdT * self.Ish.norm)
            y = -2. * np.abs(self.V3.y) / (self.heater.dRdT * self.Ish.norm)
            xerr = np.sqrt(self.V3.xerr**2 + self.heater.dRdT_err**2 + self.Ish.norm_err**2)
            yerr = np.sqrt(self.V3.yerr**2 + self.heater.dRdT_err**2 + self.Ish.norm_err**2)
            self._T2 = ACReading(x, y, xerr, yerr)

        return self._T2

    @property
    def Z2(self) -> ACReading:
        """
        Average 2ω surface thermal impedance at each ω
        (Z2 = T/Q, Joule heat input Q) [K / W ]
        """
        if self._Z2 is None or self._refresh_dependents:
            x = self.T2.x / self.power.norm
            y = self.T2.y / self.power.norm
            xerr = np.sqrt(self.T2.xerr**2 + self.power.norm_err**2)
            yerr = np.sqrt(self.T2.yerr**2 + self.power.norm_err**2)
            self._Z2 = ACReading(x, y, xerr, yerr)

        return self._Z2

    @property
    def power(self) -> ACReading:
        """
        Average active (x) and reactive (y) power (IEEE Std 1459-2010): S = VI*
        P_ave = I_rms V_rms = |power|
        """
        if self._power is None or self._refresh_dependents:
            x = self.V.x * self.Ish.x + self.V.y * self.Ish.y
            xerr = np.sqrt(
                (self.V.x * self.Ish.x)**2 * (self.V.xerr**2 + self.Ish.xerr**2)
                + (self.V.y * self.Ish.y)**2 * (self.V.yerr**2 + self.Ish.yerr**2)
            ) / x

            y = self.V.y * self.Ish.x - self.V.x * self.Ish.y
            yerr = np.sqrt(
                (self.V.y * self.Ish.x)**2 * (self.V.yerr**2 + self.Ish.xerr**2)
                + (self.V.x * self.Ish.y)**2 * (self.V.xerr**2 + self.Ish.yerr**2)
            ) / y

            self._power = ACReading(x, y, xerr, yerr)
        return self._power

    def set_data_limits(self, start, end):
        self.data.set_limits(start, end)
        self._Ish = None
        self._T2 = None
        self._Z2 = None
        self._dT2 = None
        self._power = None

    def plot_data(self, show=False):
        return plot.plot_measured_data(self, show)


class ApproximationWarning(Warning):
    pass
