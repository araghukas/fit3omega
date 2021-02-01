import numpy as np

from .sample import *
from .data import Data, ACReading
import fit3omega.plot as plot

ROOT2 = np.sqrt(2)


class Model:
    """manages sample configuration and measured data"""

    def __init__(self, sample, data):
        if type(sample) is str:
            self.sample = Sample(sample)
        elif type(sample) is Sample:
            self.sample = sample
        else:
            raise ValueError("invalid `sample` type; need str or fit3omega.sample.Sample")

        if type(data) is str:
            self.data = Data(data)
        elif type(data) is Data:
            self.data = data
        else:
            raise ValueError("invalid `data` type; need str or fit3omega.data.Data")

        self._V = None
        self._V3 = None
        self._Vsh = None
        self._Ish = None
        self._T2 = None
        self._dT2 = None
        self._power = None

    @property
    def shunt(self) -> Shunt:
        return self.sample.shunt

    @property
    def heater(self) -> Heater:
        return self.sample.heater

    @property
    def omegas(self) -> np.array:
        return self.data.omegas

    @property
    def V(self) -> ACReading:
        return self.data.V

    @property
    def V3(self) -> ACReading:
        return self.data.V3

    @property
    def Vsh(self) -> ACReading:
        return self.data.Vsh

    @property
    def Ish(self) -> ACReading:
        """RMS series current through heater and shunt"""
        if self._Ish is None:
            x = self.Vsh.x / self.shunt.R
            y = self.Vsh.y / self.shunt.R
            xerr = np.sqrt(self.Vsh.xerr**2 + self.shunt.err**2)
            yerr = np.sqrt(self.Vsh.yerr**2 + self.shunt.err**2)
            self._Ish = ACReading(x, y, xerr, yerr)
        return self._Ish

    @property
    def T2(self) -> ACReading:
        """Temperature wave component at twice driving frequency"""
        if self._T2 is None:
            x = ROOT2 * self.V3.x / self.heater.dRdT / self.Ish.norm()
            y = -ROOT2 * self.V3.y / self.heater.dRdT / self.Ish.norm()
            xerr = np.sqrt(self.V3.xerr**2 + self.heater.dRdT_err**2 + self.Ish.relerr()**2)
            yerr = np.sqrt(self.V3.yerr**2 + self.heater.dRdT_err**2 + self.Ish.relerr()**2)
            self._T2 = ACReading(x, y, xerr, yerr)
        return self._T2

    @property
    def power(self) -> ACReading:
        """
        Average complex power absorbed by heater line

        Great reference
            "electrical-engineering-portal.com": tinyurl.com/3yp8c5qt
        """
        if self._power is None:
            A = np.cos(self.V.phi() - self.Ish.phi())
            B = np.sin(self.V.phi() - self.Ish.phi())
            phi = self.V.phi() - self.Ish.phi()  # power factor angle

            dphi = np.sqrt(self.V.abserr_phi()**2 + self.Ish.abserr_phi()**2)
            dA = np.sin(phi) * np.sin(dphi)
            dB = np.sin(dphi) * np.cos(phi)

            x = self.V.norm() * self.Ish.norm() * A  # actual average dissipated power
            y = self.V.norm() * self.Ish.norm() * B
            xerr = np.sqrt(self.V.relerr()**2 + self.Ish.relerr()**2 + (dA / A)**2)
            yerr = np.sqrt(self.V.relerr()**2 + self.Ish.relerr()**2 + (dB / B)**2)
            self._power = ACReading(x, y, xerr, yerr)
        return self._power

    def plot_data(self, show=False):
        return plot.plot_data(self, show)
