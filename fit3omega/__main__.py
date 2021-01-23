import os
import yaml


class CLI:
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
        "layers": "Layer Configuration (top to bottom)"
    }

    PROMPTS = {
        "shunt": {
            "R": "resistance [Ohm]",
            "err": "resistance error (0.01 = 1%)"
        },
        "heater": {
            "length": "line length [m]",
            "width": "line width[m]",
            "dRdT": "line dR/dT [Ohm/m]",
            "dRdT_err": "line dR/dT error (0.01 = 1%)"
        },
        "layers": {
            "name": "name",
            "height": "height [m]",
            "Cv": "heat capacity [J/m^3/K]",
            "ratio_xy": "conductivity ratio (k_x/k_y)"
        },
    }

    TYPES = {
        "shunt": {
            "R": float,
            "err": float
        },
        "heater": {
            "length": float,
            "width": float,
            "dRdT": float,
            "dRdT_err": float
        },
        "layers": {
            "name": str,
            "height": float,
            "Cv": float,
            "ratio_xy": float,
        },
    }

    MULTIPLES = ["layers"]

    Q_LINE_BLANK = "\t   --> {:>33}: "
    T_LINE_BLANK = "\t:: {}: "

    QQ_LINE_BLANK = "\t" + Q_LINE_BLANK
    TT_LINE_BLANK = "\t\t:: [{}] #{:d}: "

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
        with open(file_, 'r') as f:
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
        exit()

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
                    with open(self.INCOMPLETE, 'w') as f:
                        yaml.safe_dump(self.data, f, sort_keys=False)
                    print("==> fit3omega: dumped incomplete config...")
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

            with open(self.INCOMPLETE, 'w') as f:
                yaml.safe_dump(self.data, f, sort_keys=False)
            print("==> fit3omega: dumped incomplete config...")
            exit(1)


def _write_blank_config(path_):
    if os.path.isdir(path_):
        data = {}
        for k, v in CLI.PROMPTS.items():
            data[k] = {k_: None for k_ in v}

        for k in CLI.MULTIPLES:
            keys = data[k].keys()
            data[k] = {"1": {k: None for k in keys}}

        file_ = os.path.join(path_, CLI.BLANK)
        with open(file_, 'w') as f:
            yaml.safe_dump(data, f, sort_keys=False)
        print("wrote file: %s" % file_)
    else:
        raise NotADirectoryError("'%s'" % path_)


def _launch_cli():
    i_max = 1000
    for i, file_ in enumerate(os.listdir(os.getcwd())):
        if i > i_max:
            break
        if file_ == CLI.INCOMPLETE:
            print("==> fit3omega: detected incomplete file")

            x = ""
            while not x:
                x = input("Recover? [y/N]: ")

            if x in ['y', 'Y', 'yes', 'Yes']:
                CLI.from_incomplete(file_).run()
            else:
                break
    CLI().run()


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="create a 3-omega sample file")
    parser.add_argument("--blank",
                        help="write blank config file template",
                        action="store_true", default=True)
    parser.add_argument("--new",
                        help="start CLI to create a new config file",
                        action="store_true", default=False)
    parser.add_argument("-d",
                        help="directory that config file is written to",
                        default=None)

    args = parser.parse_args()

    d_ = os.path.expanduser(args.d) if args.d else os.getcwd()
    if args.new:
        _launch_cli()
    elif args.blank:
        _write_blank_config(d_)
