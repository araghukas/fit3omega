"""
This module contains the class that converts between voltages
and derived quantities. Namely current, power, and temperature amplitude.
"""
import numpy as np

from fit3omega.data import Data, ACReading
from fit3omega.sample import SampleParameters

_ROOT2 = 1.4142135623730951


class Model:
    """
    Thermal model for 2ω temperature rise calculated from AC voltages.

    All measurement data are RMS voltages, V_rms = Vx_rms + jVy_rms.

    The lock-in amplifier (e.g. SRS830) displays and outputs (X, Y) or (R, θ), where

        X^2 + Y^2 = R^2 = (V_signal_rms)^2
        X = R*cos(θ)
        Y = R*sin(θ)
    """

    def __init__(self,
                 sample: SampleParameters,
                 data: Data):

        self.sample = sample
        self.data = data

        self._refresh_dependents = False
        self._Ish = None
        self._T2 = None
        self._Z2 = None
        self._dT2 = None
        self._power = None

    @property
    def refresh(self) -> bool:
        """if true, always recalculate data-derived values (below) when accessing props"""
        return self._refresh_dependents

    @refresh.setter
    def refresh(self, b: bool) -> None:
        if b is True:
            self._refresh_dependents = True
        elif b is False:
            self._refresh_dependents = False
        else:
            raise ValueError("argument is not a boolean")

    @property
    def Ish(self) -> ACReading:
        """RMS series current at each ω"""
        if self._Ish is None or self.refresh:
            x = self.data.Vsh.x / self.sample.shunt.R
            y = self.data.Vsh.y / self.sample.shunt.R
            xerr = np.sqrt(self.data.Vsh.xerr**2 + self.sample.shunt.err**2)
            yerr = np.sqrt(self.data.Vsh.yerr**2 + self.sample.shunt.err**2)
            self._Ish = ACReading(x, y, xerr, yerr)
        return self._Ish

    @property
    def T2(self) -> ACReading:
        """peak temperature oscillations at each ω"""
        if self._T2 is None or self.refresh:
            x = 2. * np.abs(self.data.V3.x) / (self.sample.heater.dRdT * self.Ish.norm)
            y = -2. * np.abs(self.data.V3.y) / (self.sample.heater.dRdT * self.Ish.norm)
            xerr = np.sqrt(
                self.data.V3.xerr**2 + self.sample.heater.dRdT_err**2 + self.Ish.norm_err**2)
            yerr = np.sqrt(
                self.data.V3.yerr**2 + self.sample.heater.dRdT_err**2 + self.Ish.norm_err**2)
            self._T2 = ACReading(x, y, xerr, yerr)

        return self._T2

    @property
    def Z2(self) -> ACReading:
        """
        Average 2ω surface thermal impedance at each ω
        (Z2 = T/Q, Joule heat input Q) [K / W ]
        """
        if self._Z2 is None or self.refresh:
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
            x = self.data.V.x * self.Ish.x + self.data.V.y * self.Ish.y
            xerr = np.sqrt(
                (self.data.V.x * self.Ish.x)**2 * (self.data.V.xerr**2 + self.Ish.xerr**2)
                + (self.data.V.y * self.Ish.y)**2 * (self.data.V.yerr**2 + self.Ish.yerr**2)
            ) / x

            y = self.data.V.y * self.Ish.x - self.data.V.x * self.Ish.y
            yerr = np.sqrt(
                (self.data.V.y * self.Ish.x)**2 * (self.data.V.yerr**2 + self.Ish.xerr**2)
                + (self.data.V.x * self.Ish.y)**2 * (self.data.V.xerr**2 + self.Ish.yerr**2)
            ) / y

            self._power = ACReading(x, y, xerr, yerr)
        return self._power
