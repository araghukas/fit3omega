import os
import yaml


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
            "ky": float,
            "ratio_xy": float,
            "Cv": float
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
            "height": True,
            "ky": True,
            "ratio_xy": True,
            "Cv": True
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


def _write_blank_config():
    data = {}
    for k, v in CLI.PROMPTS.items():
        data[k] = {k_: None for k_ in v}

    for k in CLI.MULTIPLES:
        keys = data[k].keys()
        data[k] = {"1": {k: None for k in keys}}

    file_ = os.path.join(os.getcwd(), CLI.BLANK)
    with open(file_, 'w') as f:
        yaml.safe_dump(data, f, sort_keys=False)
    print("wrote file: %s" % file_)
    exit(0)


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


def _plot_data(sample, data, label, show=True):
    from .model import Model
    m = Model(sample, data)
    fig = m.plot_data(show=show)
    save_name = os.path.abspath(data).strip(".csv") + label + "_plot.pdf"
    fig.savefig(save_name)
    print("==> fit3omega: saved plot as %s" % save_name)
    exit(0)


def _plot_diagnostics(sample, data, label, show=True):
    from .model import Model
    from .plot import plot_diagnostics
    m = Model(sample, data)
    fig = plot_diagnostics(m, show)
    save_name = os.path.abspath(data).strip(".csv") + label + "_diagnostic.pdf"
    fig.savefig(save_name)
    print("==> fig3omega: saved plot as %s" % save_name)
    exit(0)


def _plot_compare_data(sample, data1, data2, label, show=True):
    from .model import Model
    from .plot import plot_compare_measured_data
    m1 = Model(sample, data1)
    m2 = Model(sample, data2)
    fig = plot_compare_measured_data(m1, m2, show)
    save_name = os.path.abspath(data1).strip(".csv") + label + "_compare_plot.pdf"
    fig.savefig(save_name)
    print("==> fit3omega: saved plot as %s" % save_name)
    exit(0)


def _fit_data(sample, data, frac, niter, data_lims, b_type, label, show=True):
    from .fit_general import FitGeneral
    niter = 30 if niter is None else niter
    fg = FitGeneral(sample, data, b_type)
    start, end = (0, len(fg.data)) if data_lims is None else data_lims
    fg.set_data_limits(start, end)
    fg.fit(niter=niter, frac=frac)

    print()
    print(fg.result)
    print()

    fig, ax = fg.plot_fit(show=show)
    save_name = os.path.abspath(data).rstrip(".csv") + label + "_fit.pdf"
    fig.savefig(save_name)
    print("==> fit3omega: saved plot as %s" % save_name)
    exit(0)


def _launch_slider_plot(sample, data, frac, niter, data_lims, b_type, label):
    from .fit_general import FitGeneral
    from .plot import SliderPlot
    fg = FitGeneral(sample, data, b_type)
    start, end = (0, len(fg.data)) if data_lims is None else data_lims
    fg.set_data_limits(start, end)

    sp = SliderPlot(fg, frac=frac, niter=niter)
    sp.plot_initial_state()
    save_name = os.path.abspath(data).rstrip(".csv") + label + "_slider_fit.pdf"
    sp.fig.savefig(save_name)
    print("==> fit3omega: saved plot as %s" % save_name)
    exit(0)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="create a 3-omega sample file or plot data")
    parser.add_argument("--blank",
                        help="write blank config file template",
                        action="store_true", default=False)
    parser.add_argument("--new",
                        help="start CLI to create a new config file",
                        action="store_true", default=False)
    parser.add_argument("--show",
                        help="toggle showing plot window after plotting",
                        action="store_true", default=False)

    parser.add_argument("-plot_data",
                        help="plot 3omega voltage data from sample config and data csv",
                        nargs=2, type=str, default=None)

    parser.add_argument("-plot_diagnostics",
                        help="plot power and 3omega voltage over current cubed "
                             "from sample config and data csv",
                        nargs=2, type=str, default=None)

    parser.add_argument("-compare_data",
                        help="plot 3omega voltage from two different runs on the same sample",
                        nargs=3, type=str, default=None)

    parser.add_argument("-fit_data",
                        help="fit and plot 3omega voltage data from sample config and data csv",
                        nargs=2, type=str, default=None)

    parser.add_argument("-slider_plot",
                        help="create a responsive plot where fit parameters are adjustable",
                        nargs=2, type=str, default=None)

    parser.add_argument("-frac",
                        help="perturbation fraction in interval (0,1)",
                        type=float, default=1.0)

    parser.add_argument("-niter",
                        help="number of basinhopping iterations",
                        type=int, default=None)

    parser.add_argument("-data_lims",
                        help="indices of first and last data points",
                        nargs=2, type=int, default=None)

    parser.add_argument("-b_type",
                        help="boundary type for general heat transfer model",
                        type=str, default='i')

    parser.add_argument("-label",
                        help="label to append to saved files",
                        type=str, default="")

    args = parser.parse_args()

    if args.new:
        _launch_cli()
    elif args.blank:
        _write_blank_config()
    elif args.plot_data:
        _plot_data(*args.plot_data, show=args.show, label=args.label)
    elif args.plot_diagnostics:
        _plot_diagnostics(*args.plot_diagnostics, show=args.show, label=args.label)
    elif args.compare_data:
        _plot_compare_data(*args.compare_data, show=args.show, label=args.label)
    elif args.fit_data:
        _fit_data(*args.fit_data,
                  frac=args.frac,
                  niter=args.niter,
                  data_lims=args.data_lims,
                  show=args.show,
                  b_type=args.b_type,
                  label=args.label)
    elif args.slider_plot:
        _launch_slider_plot(*args.slider_plot,
                            frac=args.frac,
                            niter=args.niter,
                            data_lims=args.data_lims,
                            b_type=args.b_type,
                            label=args.label)
    else:
        print("==> fit3omega: no arguments detected")
        exit(0)
