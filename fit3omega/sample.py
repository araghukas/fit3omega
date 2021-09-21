import yaml
from dataclasses import dataclass
from typing import Dict, Union
from os.path import expanduser

"""
Containers for crucial sample data.

ALL quantities should be specified in S.I. units:
    
    m, kg, J, Ohm, K
    
Use exponentiation instead of unit prefixes. (ex: 2e-9 [m] for 2 nm)
"""


class Sample:
    """The sample data reader and container"""

    FIT_PARAMS = ["kys", "ratio_xys", "Cvs", "Rcs"]

    def __init__(self, config_file):
        self.shunt = None
        self.heater = None
        self.layers = []
        self.config_file = None
        self._load_config(config_file)

    def as_dict(self):
        """reconstruct the current configuration as a YAML readable dictionary"""
        d = {
            "heater": {
                "Cv": self.heater.Cv,
                "Rc": self.heater.Rc,
                "dRdT": self.heater.dRdT,
                "dRdT_err": self.heater.dRdT_err,
                "height": self.heater.height,
                "length": self.heater.length,
                "width": self.heater.width
            },
            "layers": {},
            "shunt": {
                "R": self.shunt.R,
                "err": self.shunt.err
            }
        }

        for i, layer in enumerate(self.layers):
            d["layers"][str(i)] = {
                "name": layer.name,
                "height": layer.height,
                "ky": layer.ky,
                "ratio_xy": layer.ratio_xy,
                "Cv": layer.Cv,
                "Rc": layer.Rc
            }

        return d

    def __str__(self):
        """reconstructs and returns the config file as a string"""
        return yaml.safe_dump(self.as_dict())

    def _load_config(self, config_file):
        """read the config file into instance attributes"""
        config_file = expanduser(config_file)
        with open(config_file) as f:
            d = yaml.safe_load(f)
            if type(d) is str:
                raise SampleError(f"invalid sample configuration file:\n {config_file}")
        for k, v in d.items():
            if k == "shunt":
                self.shunt = Shunt(**{_: float(x) for _, x in v.items()})
            elif k == "heater":
                self.heater = Heater(**{_: float(x) for _, x in v.items()})
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

    @property
    def Rcs(self):
        return [layer.Rc for layer in self.layers]


class Heater:
    """properties of metal heater/transducer/thermometer line"""

    def __init__(self, length: float, width: float, dRdT: float, dRdT_err: float,
                 Cv: float = 0.0, height: float = 0.0, Rc: float = 0.0):
        self.length = self._length = length
        self.width = self._width = width
        self.height = self._height = height
        self.dRdT = self._dRdT = dRdT
        self.Cv = self._Cv = Cv
        self.Rc = self._Rc = Rc
        self.dRdT_err = dRdT_err

    def reset(self):
        self.dRdT = self._dRdT
        self.width = self._width
        self.height = self._height
        self.length = self._length
        self.Cv = self._Cv
        self.Rc = self._Rc

    def modify(self, attr_name, new_value):
        if not attr_name.startswith("_"):
            try:
                self.__setattr__(attr_name, new_value)
            except KeyError:
                pass


class SampleError(Exception):
    pass


@dataclass
class Shunt:
    """properties of the shunt resistor (current measurement)"""
    R: float
    err: float


@dataclass
class Layer:
    """properties of a layer in the thermal model"""
    name: str
    height: str or float
    ky: str or float
    ratio_xy: str or float
    Cv: str or float  # heat capacity [J/m^3/K]
    Rc: str or float = 0.0  # Olson model only, thermal contact resistance to layer below

    def as_dict(self) -> Dict[str, Union[str, float]]:
        return {
            "name": self.name,
            "height": self.height,
            "ky": self.ky,
            "ratio_xy": self.ratio_xy,
            "Cv": self.Cv,
            "Rc": self.Rc,
        }
