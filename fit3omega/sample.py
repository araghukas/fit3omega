"""
Containers for crucial sample data.

ALL quantities should be specified in S.I. units:

    m, kg, J, Ohm, K

Use exponentiation instead of unit prefixes. (ex: 2e-9 [m] for 2 nm)
"""
import copy
import os
import yaml
import numpy as np
from dataclasses import dataclass
from typing import List, Tuple, Sequence, Dict
from collections import OrderedDict


@dataclass
class Heater:
    """container for the heater parameters"""
    length: float
    width: float
    dRdT: float
    dRdT_err: float
    height: float = 0.0  # height of the heater line [m]
    Cv: float = 0.0  # heat capacity [J/m^3/K]
    Rc: float = 0.0  # thermal contact resistance [K/W]


@dataclass
class Layer:
    """container for the layer parameters"""
    name: str
    height: float
    ky: float
    ratio_xy: float  # in-plane to cross-plane thermal conductivity ratio
    Cv: float  # heat capacity [J/m^3/K]
    Rc: float = 0.0  # thermal contact resistance [K/W]


@dataclass
class ShuntResistor:
    """container for the shunt resistor parameters"""
    R: float
    err: float


@dataclass
class SampleParameters:
    """
    A sample configuration. Contains the relevant parameters connecting
    measured voltages to temperature oscillation amplitudes.
    """
    heater: Heater
    layers: List[Layer]
    shunt: ShuntResistor
    fit_indices: List[Tuple[int, int]]

    FIELDS = ("kys", "ratio_xys", "Cvs", "Rcs")

    def __post_init__(self):
        self._validate_layer_names()

    def _validate_layer_names(self) -> None:
        """make no two layers have the same name"""
        layer_names = [layer.name for layer in self.layers]
        names_set = set(layer_names)
        if len(names_set) != len(layer_names) or len(names_set) == 0:
            raise ConfigFileError("layer must have unique names.")

    @property
    def parameters(self) -> Dict[str, float]:
        """a descriptive dict of variable fitting parameters"""
        d = OrderedDict()
        args = self.x
        i = 0
        for idx in self.fit_indices:
            p, q = idx
            param_name = self.FIELDS[p].rstrip('s')
            layer_name = self.layers[q].name
            d[f"{param_name}.{layer_name}"] = args[i]
            i += 1
        return d

    @property
    def argv(self) -> Tuple[List[float]]:
        """
        Get the complete set of present parameters:
        i.e. lists of each fit parameter for every layer
        """
        args = tuple()
        for param in self.FIELDS:
            values = []
            k = param.rstrip("s")
            for layer in self.layers:
                values.append(layer.__dict__[k])
            args += (values,)
        return args

    @property
    def x(self) -> np.ndarray:
        """initial values of the selected fitting parameters; optimizer starting point"""
        _x = np.zeros(len(self.fit_indices))
        argv = self.argv
        for i, idx in enumerate(self.fit_indices):
            _x[i] = argv[idx[0]][idx[1]]
        return _x

    @property
    def state(self) -> Dict:
        """a dictionary of the present sample parameters"""
        return {
            "heater": {
                "Cv": float(self.heater.Cv),
                "Rc": float(self.heater.Rc),
                "dRdT": float(self.heater.dRdT),
                "dRdT_err": float(self.heater.dRdT_err),
                "height": float(self.heater.height),
                "length": float(self.heater.length),
                "width": float(self.heater.width)
            },
            "layers": {
                str(i): {
                    "name": layer.name,
                    "height": float(layer.height),
                    "ky": float(layer.ky),
                    "ratio_xy": float(layer.ratio_xy),
                    "Cv": float(layer.Cv),
                    "Rc": float(layer.Rc)
                } for i, layer in enumerate(self.layers)
            },
            "shunt": {
                "R": self.shunt.R,
                "err": self.shunt.err
            }
        }

    def modify_heater(self,
                      field_name: str,
                      new_value: float) -> None:
        """modify the value of a field in the `heater` attribute"""
        self.heater.__setattr__(field_name, new_value)

    def modify_layer(self,
                     layer_name: str,
                     field_name: str,
                     new_value: float) -> None:
        """modify the value of a field in a specific layer"""
        self.get_layer(layer_name).__setattr__(field_name, new_value)

    def get_layer(self, layer_name: str) -> Layer:
        """return (a ref to) the layer with the specified name"""
        for layer in self.layers:
            if layer.name == layer_name:
                return layer

        raise ValueError(f"no layer named '{layer_name}'")

    def get_value(self, layer_name: str, param_name: str) -> float:
        """return the value of the parameter from the specified layer"""
        return self.get_layer(layer_name).__getattribute__(param_name)

    def substitute(self, partial_argv: Sequence[float]) -> Tuple[List[float]]:
        """substitute the partial argument vector into complete arguments at fit indices"""
        complete_argv = self.argv
        for arg, idx in zip(partial_argv, self.fit_indices):
            i_param, i_layer = idx
            complete_argv[i_param][i_layer] = arg
        return complete_argv

    def write_state(self, filename: str) -> None:
        """write the present state as a new configuration file"""
        filename = os.path.expanduser(filename)
        with open(filename, 'w') as f:
            yaml.safe_dump(self.state, f)

    def copy(self) -> 'SampleParameters':
        """return a copy of this instance"""
        return copy.deepcopy(self)


class ConfigFileError(Exception):
    """for error encountered while reading sample configuration files"""


def load_sample_parameters(config_filename: str) -> SampleParameters:
    """Read a sample configuration in YAML format into  SampleParameters instance"""
    filename = os.path.expanduser(config_filename)
    with open(filename) as f:
        d = yaml.safe_load(f)
        if type(d) is str:
            raise ConfigFileError(f"invalid sample configuration file:\n {filename}")

    d = convert_parameters_to_float(d)

    if "heater" not in d:
        raise ConfigFileError("missing 'heater' section.")
    if "layers" not in d:
        raise ConfigFileError("missing 'layers' section.")
    if "shunt" not in d:
        raise ConfigFileError("missing 'shunt' section.")

    try:
        heater = Heater(**d["heater"])
    except TypeError as te:
        kw_name = str(te).split(' ')[-1]
        raise ConfigFileError(f"invalid parameter in 'heater' section: {kw_name}.")

    try:
        shunt = ShuntResistor(**d["shunt"])
    except TypeError as te:
        kw_name = str(te).split(' ')[-1]
        raise ConfigFileError(f"invalid parameter in 'shunt' section: {kw_name}.")

    for layer_id, layer_params in d["layers"].items():
        if type(layer_params) is not dict:
            raise ConfigFileError("layers section is not double nested.")

    fit_indices = []
    layers = []
    for i_layer, label in enumerate(d["layers"]):
        layer_params = d["layers"][label]
        layer_kwargs = layer_params.copy()
        for param_name in layer_params:
            value = layer_params[param_name]
            if type(value) is str and param_name != 'name':
                i_param = SampleParameters.FIELDS.index(param_name + 's')
                fit_indices.append((i_param, i_layer))
                try:
                    layer_kwargs[param_name] = float(value.rstrip("*"))
                except ValueError as ve:
                    val_string = str(ve).split(' ')[-1]
                    raise ConfigFileError(f"invalid parameter value in layer {label}: {val_string}")
        layers.append(Layer(**layer_kwargs))

    return SampleParameters(heater, layers, shunt, fit_indices)


def convert_parameters_to_float(d: dict) -> dict:
    """
    Try to convert all numerical parameters to float, because the
    YAML reader may fail to do so.
    """

    # convert everything to float if possible
    for k, v in d["heater"].items():
        try:
            d["heater"][k] = float(v)
        except ValueError:
            raise ConfigFileError(f"failed to convert heater parameter '{k}' to float.")
    for k, v in d["layers"].items():
        for k1, v1 in v.items():
            typ_v1 = type(v1)
            if typ_v1 is str and k1 != 'name' and not v1.endswith('*'):
                try:
                    v[k1] = float(v1)
                except ValueError:
                    raise ConfigFileError(
                        f"failed to convert layer '{k}' parameter "
                        f"'{k1}' value {v1} (type {typ_v1}) to float."
                    )
    for k, v in d["shunt"].items():
        try:
            d["shunt"][k] = float(v)
        except ValueError:
            raise ConfigFileError(
                f"failed to convert shunt parameter "
                f"'{k}' value '{v}' (type {type(v)}) to float."
            )

    return d
