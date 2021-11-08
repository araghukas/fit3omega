"""a command line interface for fit3omega"""
import os
import argparse

from .fit import Fit3omega
from .plots import plot_fitted_data, plot_measured_data
from .slider_gui import SliderFit


def main(args: argparse.Namespace) -> None:
    ft = Fit3omega(args.sample_file, args.data_file)
    if args.data_lims:
        a, b = args.data_lims
        ft.data.set_limits(a, b)
    if args.fit:
        _run_fit(args, ft)
        exit()
    if args.plot:
        _plot_measured_data(args, ft)
        exit()

    _launch_slider_plot(args)


def _run_fit(args: argparse.Namespace, ft: Fit3omega) -> None:
    """create a fitter instance, run a fit, and display the results"""
    ft.fit()
    print(ft.result)

    if args.plot:
        save_name = os.path.abspath(args.data_file).strip(".csv") + "_fit_plot.pdf"
        fig = plot_fitted_data(ft, show=(not args.hide))
        fig.savefig(save_name)
        print("==> fit3omega: saved plot\n%s" % save_name)


def _plot_measured_data(args: argparse.Namespace, ft: Fit3omega) -> None:
    """create a plot of the measured data"""
    save_name = os.path.abspath(args.data_file).strip(".csv") + "_measured_plot.pdf"
    fig = plot_measured_data(ft, show=(not args.hide))
    fig.savefig(save_name)
    print("==> fit3omega: saved plot\n%s" % save_name)


def _launch_slider_plot(args: argparse.Namespace) -> None:
    """create and display a slider plot"""
    sf = SliderFit(args.sample_file, args.data_file)
    if args.data_lims:
        sf.data.set_limits(*args.data_lims)
    sf.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="The fit3omega command line interface.",
                                     epilog="""
                                     If the 'fit' option is selected, the data will be fitted.
                                     
                                     If the 'plot' option is selected, the data will be plotted
                                     with any fit that was done.
                                     
                                     If neither 'fit' nor 'plot' are selected, a slider plot
                                     will be displayed.
                                     """)

    parser.add_argument("sample_file",
                        help="path to YAML formatted sample configuration file.",
                        type=str)
    parser.add_argument("data_file",
                        help="path to CSV file containing experimental voltage data.",
                        type=str)

    parser.add_argument("-fit",
                        help="fit the sample parameters to the 3-omega voltage data.",
                        action='store_true',
                        default=False)

    parser.add_argument("-plot",
                        help="plot the measurement data.",
                        action='store_true',
                        default=False)

    parser.add_argument("-hide", "-d",
                        help="suppress showing plots",
                        action='store_true',
                        default=False)

    parser.add_argument("-data_lims",
                        help="limit the data range by taking data[a:b]",
                        nargs=2,
                        type=int,
                        default=None)

    main(parser.parse_args())
