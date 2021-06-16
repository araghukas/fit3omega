import yaml
from typing import NamedTuple

"""
Containers for crucial sample data.

ALL quantities should be specified in S.I. units:
    
    m, kg, J, Ohm, K
    
Use exponentiation instead of unit prefixes. (ex: 2e-9 [m] for 2 nm)
"""


class Sample:
    def __init__(self, config_file):
        self.shunt = None
        self.heater = None
        self.layers = []
        self.config_file = None
        self.load_config(config_file)

    def load_config(self, config_file):
        with open(config_file) as f:
            d = yaml.safe_load(f)
            if type(d) is str:
                raise SampleError(f"invalid sample configuration file:\n {config_file}")
        for k, v in d.items():
            if k == "shunt":
                self.shunt = Shunt(**v)
            elif k == "heater":
                self.heater = Heater(**v)
            elif k == "layers":
                for _, layer_dict in v.items():
                    try:
                        self.layers.append(Layer(**layer_dict))
                    except TypeError as te:
                        raise ValueError("unknown sample property '%s'"
                                         % str(te).split(' ')[-1].strip("'"))
            else:
                raise ValueError("unknown configuration key '%s'" % k)
        self.config_file = config_file

    @property
    def heights(self):
        return [layer.height for layer in self.layers]

    @property
    def kys(self):
        return [layer.ky for layer in self.layers]

    @property
    def ratio_xys(self):
        return [layer.ratio_xy for layer in self.layers]

    @property
    def Cvs(self):
        return [layer.Cv for layer in self.layers]


class Heater:
    """Metal heater/transducer/thermometer line properties"""

    def __init__(self, length: float, width: float, dRdT: float, dRdT_err: float):
        self.length = self._length = length
        self.width = self._width = width
        self.dRdT = self._dRdT = dRdT
        self.dRdT_err = dRdT_err

    def reset(self):
        self.dRdT = self._dRdT
        self.width = self._width
        self.length = self._length

    def modify(self, attr_name, new_value):
        if not attr_name.startswith("_"):
            try:
                self.__setattr__(attr_name, new_value)
            except KeyError:
                pass


class SampleError(Exception):
    """just a new name for it"""
    pass


class Shunt(NamedTuple):
    R: float
    err: float


class Layer(NamedTuple):
    name: str
    height: float
    ky: float
    ratio_xy: float
    Cv: float  # heat capacity [J/m^3/K]

    def as_dict(self):
        return {
            "name": self.name,
            "height": self.height,
            "ky": self.ky,
            "ratio_xy": self.ratio_xy,
            "Cv": self.Cv
        }
