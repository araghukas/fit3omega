import yaml
from typing import NamedTuple

"""
Containers for crucial data.

ALL quantities should be specified in S.I. units:
    
    m, kg, J, Ohm, K
    
Use exponentiation instead of unit prefixes. (ex: 2e-9 [m] for 2 nm)
"""


class Sample:
    def __init__(self, config_file):
        self.shunt = None
        self.heater = None
        self.layers = []
        self.load_config(config_file)

    def load_config(self, config_file):
        with open(config_file, 'r') as f:
            d = yaml.safe_load(f)
        for k, v in d.items():
            if k == "shunt":
                self.shunt = Shunt(**v)
            elif k == "heater":
                self.heater = Heater(**v)
            elif k == "layers":
                for k2, layer_dict in v.items():
                    self.layers.append(Layer(**layer_dict))
            else:
                raise ValueError("unrecognized configuration key '%s'" % k)

    @property
    def heights(self):
        return [layer.height for layer in self.layers]

    @property
    def Cvs(self):
        return [layer.Cv for layer in self.layers]


class Shunt(NamedTuple):
    """shunt resistor for current determination"""
    R: float
    err: float


class Heater(NamedTuple):
    """metal heater/transducer/thermometer line properties"""
    length: float
    width: float
    dRdT: float
    dRdT_err: float


class Layer(NamedTuple):
    name: str
    height: float
    Cv: float  # heat capacity [J/m^3/K]
    ratio_xy: float  # thermal conductivity anisotropy ratio_xy = k_x / k_y
