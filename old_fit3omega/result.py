"""the universal fit result object is defined here"""
from dataclasses import dataclass
from typing import List, Dict
import yaml

from old_fit3omega.sample import Sample


@dataclass(frozen=True)
class FitResult:
    """a container for results of a data fit"""
    sample: Sample
    fitted_kwargs: Dict[str, List[float]]  # values for parameters derived from fit

    def as_dict(self) -> dict:
        """convert results to a new sample config dictionary, but with fitted values"""
        d = self.sample.as_dict()
        for param, values in self.fitted_kwargs.items():
            for i_layer, value in enumerate(values):
                if value is None:
                    continue
                d["layers"][str(i_layer)][param.rstrip('s')] = f"{value}*"
        return d

    def summary(self) -> str:
        """return a string summarizing the result"""
        total_lines = []

        layers_initial = self.sample.as_dict()["layers"]
        layers_fitted = self.as_dict()["layers"]

        for layer_number, layer_params in layers_fitted.items():
            total_lines.append(
                "[" + layer_params["name"] + "]"
            )
            layer_lines = []
            for param, value in layer_params.items():
                if type(value) is not str or param == 'name':
                    continue
                v_fitted = self._get_numerical(value)
                v_guess = self._get_numerical(layers_initial[layer_number][param])
                diff_percent = (v_fitted - v_guess) / v_guess * 100.0
                sign = "+" if diff_percent >= 0 else "-"

                layer_lines.append(
                    "\t{}: {:.2e} ({}{:.2f}%)".format(param, v_fitted, sign, abs(diff_percent))
                )
            if len(layer_lines) == 0:
                total_lines.pop()
            else:
                for line in layer_lines:
                    total_lines.append(line)

        return "\n".join(total_lines)

    def __str__(self):
        return yaml.safe_dump(self.as_dict())

    @staticmethod
    def _get_numerical(config_val_string: str) -> float:
        """strip the trailing asterisk and return the float value"""
        return float(config_val_string.rstrip("*"))
