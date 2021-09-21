import numpy as np
from typing import List, Tuple, Dict, Union

from fit3omega.sample import Sample


class ParameterMap:
    """
    Maps all fitting parameters to a flat list, which can be sampled to fit
    an arbitrary subset using a generic n-parameter optimizer.
    """

    @property
    def initial_values(self) -> np.ndarray:
        """initial values for for the fitting parameters"""
        return np.asarray(self._initial_values)

    @property
    def integrator_ids(self) -> List[Tuple[int, int]]:
        """list of integer IDs passed to the C-extension integrators"""
        return self._ids

    @property
    def sample(self) -> Sample:
        """the associated Sample object"""
        return self._sample

    @sample.setter
    def sample(self, _sample):
        if type(_sample) is not Sample:
            raise ValueError("sample must be a fit3omega.sample.Sample instance.")
        self._sample = _sample

    def __init__(self, sample: Sample):
        self.sample = sample

        self._n_layers = len(sample.heights)
        self._full_args = []  # a flat list of fitting parameters across all layers
        for lst in [sample.__getattribute__(k) for k in Sample.FIT_PARAMS]:
            self._full_args += lst

        self._fit_indices = []  # indices identifying fit parameters
        self._initial_values = []  # initial values for fit parameters
        self._ids: List[tuple] = []  # layer and argument indices for each fit parameter

        for i, arg in enumerate(self._full_args):
            if type(arg) is str:
                if arg.endswith('*'):
                    self._fit_indices.append(i)
                    self._initial_values.append(float(arg.rstrip('*')))
                    self._ids.append(self._identify_fit_index(i))
                else:
                    raise ValueError("invalid sample parameter '%s'" % arg)

    def substitute(self, fitting_args) -> Tuple[List[float]]:
        """return a complete parameter set with the fitting values substituted in"""
        objective_args = self._full_args.copy()
        for j, index in enumerate(self._fit_indices):
            objective_args[index] = fitting_args[j]

        # reconstruct complete args
        args_complete = tuple()
        for k in range(len(self._full_args) // self._n_layers):
            i_min = k * self._n_layers
            i_max = i_min + self._n_layers
            args_complete += (objective_args[i_min:i_max],)

        return args_complete

    def get_fitted_kwargs(self, fitted_argv: np.ndarray) -> Dict[str, List[Union[None, float]]]:
        """return a partial parameter set where only fitted values are not None"""
        fitted_kwargs = {param: ([None] * self._n_layers) for param in Sample.FIT_PARAMS}
        for i, value in enumerate(fitted_argv):
            arg_index, layer_index = self._ids[i]
            param = Sample.FIT_PARAMS[arg_index]
            fitted_kwargs[param][layer_index] = value
        return fitted_kwargs

    def get_complete_kwargs(self, fitted_argv: np.ndarray) -> Dict[str, List[float]]:
        """construct keyword arguments for the temperature function from the optimized output vector"""
        fitted_args = self.substitute(fitted_argv)
        fitted_kwargs = {}
        for i, argv in enumerate(fitted_args):
            param = Sample.FIT_PARAMS[i]
            fitted_kwargs[param] = argv
        return fitted_kwargs

    def _identify_fit_index(self, index) -> tuple:
        """return the parameter and layer indices based on index in total args"""
        i_arg = index // self._n_layers
        i_layer = index - i_arg * self._n_layers
        return i_arg, i_layer

    def _get_initial_kwargs(self) -> Dict[str, List[float]]:
        """starting parameter names and value sets"""
        return {name: self.sample.__getattribute__(name).copy() for name in Sample.FIT_PARAMS}
