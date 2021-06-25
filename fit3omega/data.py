import pandas as pd
import numpy as np


class ACReading:
    def __init__(self, x, y, xerr, yerr):
        self.x = x  # "real" part
        self.y = y  # "imag" part
        self.xerr = xerr  # (abserr x) / abs x
        self.yerr = yerr  # (abserr y) / abs y
        self.norm_sq = x**2 + y**2
        self.norm = np.sqrt(self.norm_sq)  # R
        self.norm_err = np.sqrt(xerr**2 + yerr**2)  # (abserr R) / R


class Data:
    """
    Handles the measurement data, given the CSV file containing it
    """
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

        # Data limits and voltage readings
        self._start = 0
        self._end = None
        self._V = None
        self._V3 = None
        self._Vsh = None

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
            raise ValueError("data vs. error-data length mismatch")

    def __len__(self):
        return len(self.data)

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
        """2*PI*f for each measurement frequency f"""
        omegas_ = 2 * np.pi * self.data['freq'].values
        return np.ascontiguousarray(omegas_)

    @property
    def V(self) -> ACReading:
        """(V_x,RMS, V_y,RMS, d(V_x,RMS), d(V_y,RMS)"""
        if self._V is None:
            self._V = self._get_reading("V")
        return self._V

    @property
    def V3(self) -> ACReading:
        """(V3_x,RMS, V3_y,RMS, d(V3_x,RMS), d(V3_y,RMS)"""
        if self._V3 is None:
            self._V3 = self._get_reading("V3")
        return self._V3

    @property
    def Vsh(self) -> ACReading:
        """(Vsh_x,RMS, Vsh_y,RMS, d(Vsh_x,RMS), d(Vsh_y,RMS)"""
        if self._Vsh is None:
            self._Vsh = self._get_reading("Vsh")
        return self._Vsh

    def _get_reading(self, key) -> ACReading:
        args = tuple()
        for k in self.CSV_COLS[key]:
            # average voltages (x, y)
            args += (self.data[k].values,)
        for k in self.CSV_COLS['d' + key]:
            # standard deviations (xerr, yerr)
            data_values = self.data[k[1:]].values
            data_values[data_values == 0] = 1e-12  # avoid division errors
            args += (self.error[k].values / data_values,)
        return ACReading(*args)


def zero_error_data(df: pd.DataFrame) -> pd.DataFrame:
    cols = ['d' + c for c in df.columns]
    return pd.DataFrame(np.zeros((len(df), len(cols))), columns=cols)
