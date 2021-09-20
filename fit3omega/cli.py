"""
A CLI for creating sample files and running various fitting operations.
"""
import yaml
import os


class CLI:
    HEADER = (
            "=================================\n"
            "    fit3omega config generator   \n"
            "=================================\n"
            "%s\n"
            "\n"
            "Enter sample parameters.\n"
            % os.getcwd()
    )

    SECTION_TITLES = {
        "shunt": "Shunt Resistor",
        "heater": "Heater/Transducer",
        "layers": ("Layer Configuration (top to bottom)\n"
                   "\t   Enter a guess for unknown values followed by '*' to fit")
    }

    PROMPTS = {
        "shunt": {
            "R": "resistance [Ohm]",
            "err": "resistance error (0.01 = 1%)"
        },
        "heater": {
            "Cv": "heat capacity [J/m^3/K]",
            "Rc": "thermal contact resistance [K/W]",
            "length": "line length [m]",
            "width": "line width[m]",
            "dRdT": "line dR/dT [Ohm/m]",
            "dRdT_err": "line dR/dT error (0.01 = 1%)"
        },
        "layers": {
            "name": "name",
            "height": "height [m]",
            "ky": "cross-plane thermal conductivity [W/m/K]",
            "ratio_xy": "in/cross thermal conductivity ratio",
            "Cv": "heat capacity [J/m^3/K]",
            "Rc": "thermal contact resistance [K/W]",
        },
    }

    TYPES = {
        "shunt": {
            "R": float,
            "err": float
        },
        "heater": {
            "Cv": float,
            "Rc": float,
            "length": float,
            "width": float,
            "dRdT": float,
            "dRdT_err": float
        },
        "layers": {
            "name": str,
            "height": float,
            "ky": float,
            "ratio_xy": float,
            "Cv": float,
            "Rc": float,
        },
    }

    CAN_FIT = {
        "shunt": {
            "R": False,
            "err": False,
        },
        "heater": {
            "length": False,
            "width": False,
            "dRdT": False,
            "dRdT_err": False
        },
        "layers": {
            "name": False,
            "height": False,
            "ky": True,
            "ratio_xy": True,
            "Cv": True,
            "Rc": True
        }
    }

    MULTIPLES = ["layers"]

    Q_LINE_BLANK = "\t   --> {:>40}: "
    T_LINE_BLANK = "\t:: {}: "

    QQ_LINE_BLANK = "\t" + Q_LINE_BLANK
    TT_LINE_BLANK = "\t    :: [{}] #{:d}: "

    INCOMPLETE = "_incomplete.f3oc"
    BLANK = "blank.f3oc"

    def __init__(self, data=None):
        if data is None:
            self.data = {}
            for k, v in CLI.PROMPTS.items():
                self.data[k] = {k_: None for k_ in v}
        else:
            self.data = data
        self.cleanup_incomplete = False

    @classmethod
    def from_incomplete(cls, file_):
        with open(file_) as f:
            data = yaml.safe_load(f)
        for k in CLI.MULTIPLES:
            data[k] = {k_: None for k_ in cls.PROMPTS[k]}
        cli = cls(data)
        cli.cleanup_incomplete = True
        cli.INCOMPLETE = file_
        return cli

    def run(self):
        print(CLI.HEADER)
        for k1, v in CLI.PROMPTS.items():
            print(CLI.T_LINE_BLANK.format(CLI.SECTION_TITLES[k1]))
            if k1 in self.MULTIPLES:
                self._read_multiple(k1)
            else:
                for k2 in v:
                    self._read(k1, k2)
            print()
        self.save_data()
        exit(0)

    def save_data(self):
        x = ""
        while not x:
            x = input("==> fit3omega: save configuration as? ")
        file_ = os.path.join(os.getcwd(), x)
        with open(file_, 'w') as f:
            yaml.safe_dump(self.data, f, sort_keys=False)
        if self.cleanup_incomplete:
            os.remove(self.INCOMPLETE)

    def _read_multiple(self, k1):
        N = 0
        while N == 0:
            try:
                N = int(input(self.Q_LINE_BLANK.format("How many")))
            except ValueError:
                N = 0

        keys = self.data[k1].keys()
        self.data[k1] = {}
        for i in range(N):
            print(self.TT_LINE_BLANK.format(k1, i + 1))
            d = {k: None for k in keys}
            for k2 in self.PROMPTS[k1]:
                x = None
                while not x:
                    x = input(self.QQ_LINE_BLANK.format(self.PROMPTS[k1][k2]))
                try:
                    t = self.TYPES[k1][k2]
                    d[k2] = t(x)
                except ValueError:
                    if self.CAN_FIT[k1][k2] and x.endswith('*'):
                        try:
                            d[k2] = self.TYPES[k1][k2](x.rstrip('*'))
                            d[k2] = str(d[k2]) + '*'
                            continue
                        except ValueError:
                            pass

                    with open(self.INCOMPLETE, 'w') as f:
                        yaml.safe_dump(self.data, f, sort_keys=False)
                    print("==> fit3omega: input error '%s'; dumped incomplete config..." % x)
                    exit(1)
            self.data[k1][str(i + 1)] = d

    def _read(self, k1, k2):
        if self.data[k1][k2]:
            return

        x = None
        while not x:
            x = input(self.Q_LINE_BLANK.format(self.PROMPTS[k1][k2]))
        try:
            t = self.TYPES[k1][k2]
            self.data[k1][k2] = t(x)
        except ValueError:
            # delete null entries
            for k1, v1 in self.data.items():
                for k2, v2 in v1.items():
                    if v2 is None:
                        del v2

            if self.CAN_FIT[k1][k2] and x.endswith('*'):
                try:
                    self.data[k1][k2] = self.TYPES[k1][k2](x.rstrip('*'))
                    self.data[k1][k2] = str(self.data[k1][k2]) + '*'
                    return
                except ValueError:
                    pass

            with open(self.INCOMPLETE, 'w') as f:
                yaml.safe_dump(self.data, f, sort_keys=False)
            print("==> fit3omega: dumped incomplete config...")
            exit(1)