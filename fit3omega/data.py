import pandas as pd
import numpy as np
from typing import NamedTuple


class Voltage(NamedTuple):
    x: np.array
    y: np.array


class Data:
    def __init__(self, data_csv: str, error_csv: str = None):
        self._data = pd.read_csv(data_csv, header="infer")

        if error_csv:
            self._error = pd.read_csv(error_csv, header="infer")
        else:
            error_csv = data_csv.split(".")[:-1] + ".error.csv"
            try:
                self._error = pd.read_csv(error_csv, header="infer")
                self.error = self._error.copy()
            except FileNotFoundError:
                self._error = None

        # Data limits
        self._start = None
        self._end = None

    def set_limits(self, start: int, end: int):
        self._start = int(start)
        self._end = int(end)

    def reset_limits(self):
        self._start = None
        self._end = None

    def drop(self, row_index):
        self._data.drop(row_index, axis=0)
        if self._error:
            self._error.drop(row_index, axis=0)

    @property
    def data(self):
        return self._data[self._start:self._end]

    @property
    def error(self):
        if not self._error:
            raise ValueError("no error data has been initialized")
        return self._error[self._start:self._end]

    @error.setter
    def error(self, error_csv):
        if self._error:
            raise ValueError("error data already set")

        e = pd.read_csv(error_csv, header="infer")

        for i in range(len(self.data['freq'].values)):
            if e['freq'].values[i] != self.data['freq'].values[i]:
                raise ValueError("frequency mismatch")

        if len(e['freq']) != len(self.data['freq']):
            raise ValueError("data length mismatch")
        self._error = e

    @property
    def V(self):
        return Voltage(self.data['Vs_1w'].values, self.data['Vs_1w_o'].values)

    @property
    def V3w(self):
        return Voltage(self.data['Vs_3w'].values, self.data['Vs_3w_o'].values)

    @property
    def Vsh(self):
        return Voltage(self.data['Vsh_1w'].values, self.data['Vsh_1w_o'].values)