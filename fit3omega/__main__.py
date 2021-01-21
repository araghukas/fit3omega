import sys
import os
import yaml


class CLI(object):
    HEADER = (
            "=================================\n"
            "    fit3omega config generator   \n"
            "=================================\n"
            "%s\n"
            "\n"
            "Enter sample parameters."
            % os.getcwd()
    )

    SECTION_TITLES = {
        "shunt": "Shunt Resistor",
        "heater": "Heater/Transducer",
        "nanowires": "Nanowire Sample",
        "layers": "Layer Configuration"
    }

    PROMPTS = {
        "shunt": {
            "R_shunt": "shunt resistance [Ohm]"
        },
        "heater": {
            "line_length": "line length [m]",
            "line_width": "line width[m]",
            "line_dRdT": "line dR/dT [Ohm/m]",
        },
        "nanowires": {
            "height": "nanowire height [m]",
            "width": "nanowire diameter [m]",
            "pitch": "nanowire array pitch [m]",
        },
        "layers": {
            "name": "layer name",
            "thickness": "layer thickness",
            "Cv": "layer heat capacity [J/m^3/K]",
            "ratio_xy": "conductivity ratio (k_x/k_y)"
        },
    }

    TYPES = {
        "shunt": {
            "R_shunt": float
        },
        "heater": {
            "line_length": float,
            "line_width": float,
            "line_dRdT": float,
        },
        "nanowires": {
            "height": float,
            "width": float,
            "pitch": float
        },
        "layers": {
            "name": str,
            "thickness": float,
            "Cv": float,
            "ratio_xy": float,
        },
    }

    Q_LINE_BLANK = "\t   --> {:>33}: "
    T_LINE_BLANK = "\t:: {}"

    def __init__(self, data=None):
        if data is None:
            self.data = {}
            for k, v in self.PROMPTS.items():
                if type(v) is dict:
                    self.data[k] = {k_: None for k_ in v}
                else:
                    self.data[k] = None
        else:
            self.data = data

    @classmethod
    def from_incomplete(cls, file_):
        with open(file_, 'r') as f:
            data = yaml.safe_load(f)
        return cls(data)

    def start(self):
        print(CLI.HEADER)
        for k1, v in CLI.PROMPTS.items():
            print(CLI.T_LINE_BLANK.format(CLI.SECTION_TITLES[k1]))
            for k2 in v:
                self.get_input_value(k1, k2)
            print()

        self.save_data()

    def get_input_value(self, k1, k2):
        if self.data[k1][k2] is None:
            x = None
            while not x:
                x = input(self.Q_LINE_BLANK.format(self.PROMPTS[k1][k2]))
            try:
                t = self.TYPES[k1][k2]
                self.data[k1][k2] = t(x)
            except ValueError:
                with open(os.path.join(os.getcwd(), "_incomplete.f3oc")) as f:
                    yaml.safe_dump(self.data, f)
                print("==> fit3omega: dumped incomplete config...")

    def save_data(self):
        x = ""
        while not x:
            x = input("==> fit3omega: save configuration as: ")
        file_ = os.path.join(os.getcwd(), x)
        with open(file_, 'w') as f:
            yaml.safe_dump(self.data, f)


def _write_blank_config(path_):
    if os.path.isdir(path_):
        data = {}
        for k, v in CLI.PROMPTS.items():
            if type(v) is dict:
                data[k] = {k_: None for k_ in v}
            else:
                data[k] = None
        file_ = os.path.join(path_, "blank.f3oc")
        with open(file_, 'w') as f:
            yaml.safe_dump(data, f)
        print("wrote file: %s" % file_)
    else:
        raise NotADirectoryError("'%s'" % path_)


def _launch_cli():
    for file_ in os.listdir(os.getcwd()):
        if file_ == "_incomplete.f3oc":
            print("==> fit3omega: detected incomplete file")

            x = ""
            while not x:
                x = input("Recover? [Y/n]")

            if x in ['y', 'Y', 'yes', 'Yes']:
                CLI.from_incomplete("./incomplete.f3oc").start()
                return
            else:
                break
    CLI().start()


if __name__ == "__main__":
    # TODO: what about the data file? error file? ... looks for files with '3omega' in the name in cwd
    if len(sys.argv) >= 2:
        d = os.path.expanduser(sys.argv[1])
        if os.path.isdir(d):
            _write_blank_config(d)
            print("==> fit3omega: wrote blank config file to '%s'" % d)
    else:
        _launch_cli()

    # TODO: test dump incomplete
    # TODO: test read incomplete
