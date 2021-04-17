import numpy as np

import matplotlib
import matplotlib.pyplot as plt

from tiwaz.post_process import extract_solver_data

import tinga
from tinga.io import write_pickle


def linestyle(solver):
    wb = solver["wb"]
    order = np.max(solver["order"])

    colors_ = {1: "green", 2: "blue", 3: "purple", 4: "red", 5: "orange"}

    kwargs_ = {
        "constant": {
            "marker": "P",
            "linestyle": "-",
            "linewidth": 3,
            "color": colors_[order],
        },
        "isentropic": {
            "marker": "o",
            "linestyle": ":",
            "linewidth": 3,
            "color": colors_[order],
        },
    }

    return kwargs_[wb]


def convergence_plot_kwargs(solver):
    return linestyle(solver)


class SplitConvergencePlot:
    """Facilitates converence plots.

    Typical usage:
       plot = SplitCovergencePlot()
       for method in methods:
           nxs, err, rate = load_data(...)
           plot.add(nxs, err, rate, method)

       plot.finalize(title_params, ylabel_params)
       plot.save(filename)

    Customization through inheritance and overwriding.
    """

    def __init__(self):
        self.fig, axs = plt.subplots(
            2, sharex=True, gridspec_kw={"hspace": 0, "height_ratios": [3, 1]}
        )

        self.ax1, self.ax2 = axs[0], axs[1]
        self.err_min, self.err_max = np.inf, 0.0

    def __del__(self):
        plt.close(self.fig)

    def add(self, resolutions, errors, rates, method):
        """Add the convergence data of 'method' to the plot."""
        plt.figure(self.fig.number)

        self.err_min = min(self.err_min, np.min(errors))
        self.err_max = max(self.err_max, np.max(errors))

        self.ax1.loglog(resolutions, errors, **self.error_style(method))

        rates = np.array([rates[0]] + [r for r in rates])
        self.ax2.step(resolutions, rates, **self.rate_style(method))

        self.tick_values = resolutions

    def finalize(self, title_params, ylabel_params, tick_values=None):
        plt.figure(self.fig.number)

        if tick_values is None:
            tick_values = self.tick_values

        # plt.minorticks_off()
        self.xticks(tick_values)
        self.xlim(tick_values)
        self.yticks()
        self.ylim()
        self.title(title_params)
        self.ylabel(ylabel_params)
        self.xlabel()

        self.fig.tight_layout()

    def save(self, filename, dpi=300):
        plt.figure(self.fig.number)
        plt.savefig(filename, dpi=dpi)

    # -- Implementation ----------------------------------------------------------
    def error_style(self, method):
        return convergence_plot_kwargs(method)

    def rate_style(self, method):
        return convergence_plot_kwargs(method)

    def format_title(self, title_params):
        return title_params

    def format_ylabel(self, ylabel_params):
        return ylabel_params

    def format_xlabel(self):
        return "Resolution $dx$"

    def xlim(self, tick_values):
        pass
        # k, K = log10_bounds(np.min(tick_values), np.max(tick_values))
        # self.ax1.set_xlim((10 ** k, 10 ** K))

    def xticks(self, tick_values):
        self.ax2.tick_params(axis="x", which="minor", bottom=True)

    def yticks(self):
        fig = self.fig

        logy_min, logy_max = log10_bounds(self.err_min, self.err_max)
        ax1, ax2 = self.ax1, self.ax2

        ax1.set_yticks(
            [10 ** k for k in range(int(logy_min), int(logy_max) + 1)], minor=True
        )
        ax1.set_yticklabels(
            ["" for k in range(int(logy_min), int(logy_max) + 1)],
            minor=True,
        )
        ax1.grid(which="major", axis="y", color="k", linewidth=0.8, alpha=1.0)
        ax1.grid(which="minor", axis="y", color="k", linewidth=0.5, alpha=0.6)

        ax2.grid(which="major", axis="y", color="k", linewidth=0.8, alpha=1.0)
        ax2.grid(which="minor", axis="y", color="k", linewidth=0.5, alpha=0.6)

        ax2.set_yticks(self.ytick_values_rate())
        ax2.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        ax2.yaxis.set_label_position("right")
        ax2.yaxis.tick_right()

        ax2.yaxis.set_tick_params(which="major", length=5, width=2)
        rtick_offset = matplotlib.transforms.ScaledTranslation(
            0.0, 5.0 / 72.0, fig.dpi_scale_trans
        )

        for l in ax2.yaxis.get_majorticklabels():
            l.set_transform(l.get_transform() + rtick_offset)

    def ytick_values_rate(self):
        return [1, 2, 3, 4, 5]

    def ylim(self):
        self.ylim_error()
        self.ylim_rate()

    def ylim_error(self):
        logy_min, logy_max = log10_bounds(self.err_min, self.err_max)
        self.ax1.set_ylim((10 ** logy_min, 10 ** logy_max))

    def ylim_rate(self):
        ticks = self.ytick_values_rate()

        self.ax2.set_ylim((min(ticks) - 0.2, max(ticks) + 0.2))

    def title(self, title_params):
        self.ax1.set_title(self.format_title(title_params))

    def xlabel(self):
        self.ax2.set_xlabel(self.format_xlabel())

    def ylabel(self, ylabel_params):
        self.ax1.set_ylabel(self.format_ylabel(ylabel_params))
        self.ax2.set_ylabel(r"Rate")


def log10_bounds(m, M, lower_clip=-200):
    if lower_clip is not None:
        m, M = max(m, 10 ** lower_clip), max(M, 10 ** lower_clip)

    k, K = int(np.floor(np.log10(m))), int(np.ceil(np.log10(M)))
    return k, K


def compute_rates(dx, err):
    log_err = np.log(err)
    log_dx = np.log(dx)

    return (log_err[1:] - log_err[:-1]) / (log_dx[1:] - log_dx[:-1])


def write_convergence_plots(results, solvers, labels, filename):
    for key in ["l1_error", "l1_eq_error"]:
        solver_keys = [labels(s) for s in solvers]

        plot = SplitConvergencePlot()
        plot_data = []

        for solver in solvers:
            result = extract_solver_data(solver, results)
            resolutions = np.array([r["dx_max"] for r in result])

            l1_err = np.array([r[key] for r in result])
            l1_rate = compute_rates(resolutions, l1_err)

            plot_data.append(
                {
                    "solver": solver,
                    "l1_err": l1_err,
                    "l1_rate": l1_rate,
                    "resolutions": resolutions,
                }
            )

            plot.add(resolutions, l1_err, l1_rate, method=solver)

        plot.finalize(key, "")
        plot.save(filename + "_" + key + ".png")
        write_pickle(filename + "_" + key + ".pkl", plot_data)
