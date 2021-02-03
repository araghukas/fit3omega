import pandas as pd
import numpy as np
from typing import NamedTuple


class ACReading(NamedTuple):
    x: np.array  # cos voltage
    y: np.array  # sin voltage
    xerr: np.array  # (abserr x) / abs x
    yerr: np.array  # (abserr y) / abs y

    def norm(self):
        """pythagorean norm"""
        return np.sqrt(self.x**2 + self.y**2)

    def phi(self):
        """phase"""
        return np.arctan(self.y / self.x)

    def abserr(self):
        """absolute error of norm"""
        return np.sqrt((self.x * self.xerr)**2 + (self.y * self.yerr)**2)

    def relerr(self):
        """relative error of norm"""
        return self.abserr() / self.norm()

    def abserr_phi(self):
        """absolute error (in radians) of phase angle"""
        r = self.y / self.x
        dr = r * np.sqrt(self.xerr**2 + self.yerr**2)
        # d[ arctan(r) ] = 1 / (1 + r**2) * dr
        return dr / (1 + r**2)

    def relerr_phi(self):
        """relative error of phase angle"""
        return self.abserr_phi() / self.phi()

    def phasor(self):
        """complex vector x + jy"""
        return self.x + 1j * self.y


class Data:
    CSV_COLS = {
        "V": ['Vs_1w', 'Vs_1w_o'],
        "V3": ['Vs_3w', 'Vs_3w_o'],
        "Vsh": ['Vsh_1w', 'Vsh_1w_o'],

        "dV": ['dVs_1w', 'dVs_1w_o'],
        "dV3": ['dVs_3w', 'dVs_3w_o'],
        "dVsh": ['dVsh_1w', 'dVsh_1w_o']
    }

    def __init__(self, data_csv: str, error_csv: str = None):
        self._data = pd.read_csv(data_csv, header="infer")
        self._data_file = data_csv

        if error_csv:
            self._error = pd.read_csv(error_csv, header="infer")
            self._error_file = error_csv
        else:
            # try substituting .error.csv at the end of the `data_csv`
            error_csv = '.'.join(data_csv.split('.')[:-1]) + ".error.csv"
            try:
                self._error = pd.read_csv(error_csv, header="infer")
                self._error_file = error_csv
            except FileNotFoundError:
                self._error = zero_error_data(self.data)
                self._error_file = None

        if len(self._data) != len(self._error):
            raise ValueError("data-error length mismatch")

        # Data limits and voltage readings
        self._start = None
        self._end = None
        self._V = None
        self._V3 = None
        self._Vsh = None

    def set_limits(self, start: int, end: int):
        self._V = None
        self._V3 = None
        self._Vsh = None
        self._start = int(start)
        self._end = int(end)

    def reset(self):
        self._start = None
        self._end = None
        self._V = None
        self._V3 = None
        self._Vsh = None
        self._data = pd.read_csv(self._data_file, header="infer")
        if self._error_file is not None:
            self._error = pd.read_csv(self._error_file, header="infer")
        else:
            self._error = zero_error_data(self._data)

    def drop_row(self, row_index):
        self._data = self._data.drop(row_index, axis=0)
        if self._error is not None:
            self._error = self._error.drop(row_index, axis=0)

    @property
    def data(self) -> pd.DataFrame:
        return self._data[self._start:self._end]

    @property
    def data_file(self) -> str:
        return self._data_file

    @property
    def error(self) -> pd.DataFrame:
        if self._error is None:
            raise ValueError("no error data has been initialized")
        return self._error[self._start:self._end]

    @error.setter
    def error(self, error_csv):
        if self._error is not None:
            raise ValueError("error data already set")

        e = pd.read_csv(error_csv, header="infer")

        for i in range(len(self.data['freq'].values)):
            if e['freq'].values[i] != self.data['freq'].values[i]:
                raise ValueError("frequency mismatch")

        if len(e['freq']) != len(self.data['freq']):
            raise ValueError("data length mismatch")
        self._error = e
        self._error_file = error_csv

    @property
    def no_error(self):
        return self._error_file is None

    @property
    def omegas(self) -> np.array:
        omegas_ = 2 * np.pi * self.data['freq'].values
        return np.ascontiguousarray(omegas_)

    @property
    def V(self) -> ACReading:
        if self._V is None:
            self._V = self._get_reading("V")
        return self._V

    @property
    def V3(self) -> ACReading:
        if self._V3 is None:
            self._V3 = self._get_reading("V3")
        return self._V3

    @property
    def Vsh(self) -> ACReading:
        if self._Vsh is None:
            self._Vsh = self._get_reading("Vsh")
        return self._Vsh

    def _get_reading(self, key) -> ACReading:
        args = tuple()
        for k in self.CSV_COLS[key]:
            # average voltages
            args += (self.data[k].values,)
        for k in self.CSV_COLS['d' + key]:
            # standard deviations
            data_values = self.data[k[1:]].values
            data_values[data_values == 0] = 1e-12  # avoid division errors
            args += (self.error[k].values / data_values,)
        return ACReading(*args)


def zero_error_data(df: pd.DataFrame) -> pd.DataFrame:
    cols = ['d' + c for c in df.columns]
    return pd.DataFrame(np.zeros((len(df), len(cols))), columns=cols)
