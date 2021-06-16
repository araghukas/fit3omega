from fit3omega.sample import Sample


class VarSample(Sample):
    """
    A variable-parameters version of the above Sample class.
    """

    def __init__(self, config_file):
        super().__init__(config_file)
        self._heights = None
        self._kys = None
        self._ratio_xys = None
        self._Cvs = None
        self.reset_params()

    def reset_params(self):
        self._heights = [self._parse_item(layer.height) for layer in self.layers]
        self._kys = [self._parse_item(layer.ky) for layer in self.layers]
        self._ratio_xys = [self._parse_item(layer.ratio_xy) for layer in self.layers]
        self._Cvs = [self._parse_item(layer.Cv) for layer in self.layers]

    @staticmethod
    def _parse_item(item):
        try:
            return float(item)
        except ValueError:
            return float(item.rstrip("*"))

    @classmethod
    def from_sample(cls, sample):
        return cls(sample.config_file)

    @property
    def heights(self):
        return self._heights

    @heights.setter
    def heights(self, x):
        if len(x) == len(self.layers):
            self._heights = x
        else:
            raise ValueError("array length incompatible with sample config")

    @property
    def kys(self):
        return self._kys

    @kys.setter
    def kys(self, x):
        if len(x) == len(self.layers):
            self._kys = x
        else:
            raise ValueError("array length incompatible with sample config")

    @property
    def ratio_xys(self):
        return self._ratio_xys

    @ratio_xys.setter
    def ratio_xys(self, x):
        if len(x) == len(self.layers):
            self._ratio_xys = x
        else:
            raise ValueError("array length incompatible with sample config")

    @property
    def Cvs(self):
        return self._Cvs

    @Cvs.setter
    def Cvs(self, x):
        if len(x) == len(self.layers):
            self._Cvs = x
        else:
            raise ValueError("array length incompatible with sample config")

    def param_modify(self, layer_name: str, attr_name: str, new_value: float) -> None:
        if layer_name is None:
            if attr_name == "heights":
                self.heights = new_value
            elif attr_name == "kys":
                self.kys = new_value
            elif attr_name == "ratio_xys":
                self.ratio_xys = new_value
            elif attr_name == "Cvs":
                self.Cvs = new_value
            else:
                raise ValueError("unrecognized attribute name `%s`" % attr_name)
            return

        for i, layer in enumerate(self.layers):
            if layer.name == layer_name:
                if attr_name == "heights":
                    x = self.heights.copy()
                    x[i] = new_value
                    self.heights = x
                elif attr_name == "kys":
                    x = self.kys.copy()
                    x[i] = new_value
                    self.kys = x
                elif attr_name == "ratio_xys":
                    x = self.ratio_xys.copy()
                    x[i] = new_value
                    self.ratio_xys = x
                elif attr_name == "Cvs":
                    x = self.Cvs.copy()
                    x[i] = new_value
                    self.Cvs = x
                else:
                    raise ValueError("unrecognized attribute name `%s`" % attr_name)
                return
