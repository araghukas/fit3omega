import os
import yaml

from fit3omega.cli import CLI


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


def _plot_data(sample, data, save_name, show=True):
    from .fit import FitGeneral
    fitter = FitGeneral(sample, data)
    fig = fitter.plot_measured_data(show=show)
    if save_name is None:
        save_name = os.path.abspath(data).strip(".csv") + "_plot.pdf"
    fig.savefig(save_name)
    print("==> fit3omega: saved plot as %s" % save_name)
    exit(0)


def _plot_compare_data(sample, data1, data2, save_name, show=True):
    from .fit import FitGeneral
    fitter1 = FitGeneral(sample, data1)
    fitter2 = FitGeneral(sample, data2)
    fig = fitter1.plot_compare_measured_data(fitter2, show)
    if save_name is None:
        save_name = os.path.abspath(data1).strip(".csv") + "_compare_plot.pdf"
    fig.savefig(save_name)
    print("==> fit3omega: saved plot as %s" % save_name)
    exit(0)


def _fit_data_general_model(sample,
                            data,
                            tol,
                            data_lims,
                            save_name,
                            show=True):
    from .fit import FitGeneral
    fg = FitGeneral(sample, data)
    start, end = (0, len(fg.model.data)) if data_lims is None else data_lims
    fg.set_data_limits(start, end)
    fg.fit(tol=tol)
    print(fg.result)

    fig, ax = fg.plot_measured_data(show=show)
    if save_name is None:
        save_name = os.path.abspath(data).rstrip(".csv") + "_fit.pdf"
    fig.savefig(save_name)
    print("==> fit3omega: saved plot as %s" % save_name)
    exit(0)


def _fit_data_linear_model(sample, data, data_lims, save_name, thresh, min_length, show=True):
    from .fit import FitLinear
    fl = FitLinear(sample, data)
    start, end = (0, len(fl.model.data)) if data_lims is None else data_lims
    fl.set_data_limits(start, end)
    fl.fit(thresh=thresh, min_length=min_length)

    print()
    print(fl.result)
    print()
    print("error: %.2e" % fl.error)
    print("Rsq: %.6f" % fl.Rsq)
    print()

    fig, ax = fl.plot_fitted_T2(show)
    if save_name is None:
        save_name = os.path.abspath(data).rstrip(".csv") + "_fit.pdf"
    fig.savefig(save_name)
    print("==> fit3omega: saved plot as %s" % save_name)
    exit(0)


def _launch_slider_plot(sample, data, data_lims, save_name, ehp):
    from .fit import FitGeneral
    from .slider import SliderPlot
    fg = FitGeneral(sample, data)
    start, end = (0, len(fg.model.data)) if data_lims is None else data_lims
    fg.set_data_limits(start, end)

    sp = SliderPlot(fg, enable_heater_params=ehp)
    sp.plot_initial_state()
    if save_name is not None:
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

    parser.add_argument("-plot_data", metavar=("<config file>", "<data file>"),
                        help="plot 3omega voltage data from sample config and data csv",
                        nargs=2, type=str, default=None)

    parser.add_argument("-compare_data",
                        metavar=("<config file>", "<data file 1>", "<data file 2>"),
                        help="plot 3omega voltage from two different runs on the same sample",
                        nargs=3, type=str, default=None)

    parser.add_argument("-fit", metavar=("<config file>", "<data file>"),
                        help="fit and plot 3omega voltage data from sample config and data csv",
                        nargs=2, type=str, default=None)

    parser.add_argument("-tol",
                        help="tolerance for fit termination",
                        type=float,
                        default=1e-8)

    parser.add_argument("-slider_plot", metavar=("<config file>", "<data file>"),
                        help="create a responsive plot where fit parameters are adjustable",
                        nargs=2, type=str, default=None)

    parser.add_argument("--enable_heater_params",
                        help="allow slider plot to adjust heater dRdT, width, and length",
                        action="store_true", default=False)

    parser.add_argument("-data_lims", metavar=("<lower>", "<upper>"),
                        help="indices of first and last data points",
                        nargs=2, type=int, default=None)

    parser.add_argument("-save_name",
                        help="output file name",
                        type=str, default=None)

    parser.add_argument("-fit_linear", metavar=("<config file>", "<data file>"),
                        help="fit and plot 3omega voltage data from sample config and data csv",
                        nargs=2, type=str, default=None)

    parser.add_argument("-thresh",
                        help="Rsq threshold for determining continuous linear subset data",
                        type=float, default=0.99)

    parser.add_argument("-min_length",
                        help="minimum length of continuous linear subset data",
                        type=int, default=5)

    args = parser.parse_args()

    _save_name = os.path.expanduser(args.save_name) if args.save_name else args.save_name

    if args.new:
        _launch_cli()
    elif args.blank:
        _write_blank_config()
    elif args.plot_data:
        _plot_data(*args.plot_data, show=args.show, save_name=args.save_name)
    elif args.compare_data:
        _plot_compare_data(*args.compare_data, show=args.show, save_name=args.save_name)
    elif args.fit:
        _fit_data_general_model(*args.fit,
                                tol=args.tol,
                                data_lims=args.data_lims,
                                show=args.show,
                                save_name=args.save_name)
    elif args.fit_linear:
        _fit_data_linear_model(*args.fit_linear,
                               data_lims=args.data_lims,
                               save_name=args.save_name,
                               thresh=args.thresh,
                               min_length=args.min_length,
                               show=args.show)
    elif args.slider_plot:
        _launch_slider_plot(*args.slider_plot,
                            data_lims=args.data_lims,
                            save_name=args.save_name,
                            ehp=args.enable_heater_params)
    else:
        print("==> fit3omega: no arguments, exiting.")
        exit(0)

# TODO: double check this module still works
# TODO: final edits and version bump
