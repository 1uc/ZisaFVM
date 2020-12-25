#!/usr/bin/env python3

import os
import sys

import h5py
import numpy as np
import argparse

import matplotlib.pyplot as plt

import tiwaz
from tiwaz.post_process import find_grid
from tiwaz.post_process import find_steady_state_file
from tiwaz.post_process import find_data_files
from tiwaz.post_process import load_grid

from tiwaz.scatter_plot import ScatterPlot
from tiwaz.scatter_plot import ScatterSemilogY
from tiwaz.tri_plot import TriPlot
from tinga.io import write_pickle


def load(data_file):
    steady_state_file = find_steady_state_file(os.path.dirname(data_file))
    if not os.path.isfile(steady_state_file):
        steady_state_file = None

    return tiwaz.post_process.load_data(data_file, steady_state_file)


def load_array(h5_file, key):
    with h5py.File(h5_file, "r") as h5:
        return np.array(h5[key])


def tri_plot(grid, q):
    plot = TriPlot()
    plot.color_plot(grid, q)
    plt.show()


def stencil_indicator(grid, I):
    data = np.zeros(grid.n_cells)
    data[I] = 1.0
    data[I[0]] = 2.0

    return data


def plot_all(grid, data_files, keys, delta_keys, log_keys, with_ghost_cells=False):
    n_dims = 2 if grid.vertex_indices.shape[1] == 3 else 3

    for f in data_files:
        u = load(f)

        for key in keys:
            if n_dims == 2:
                plot = TriPlot()
                plot.color_plot(grid, u[key])
                plot.save(plot.filename(f, key))

            plot = ScatterPlot(with_ghost_cells)
            plot(grid, u[key])
            plot.save(plot.filename(f, key))

        for key in delta_keys:
            if n_dims == 2:
                plot = TriPlot("PRGn", symmetric_limits=True)
                plot.color_plot(grid, u[key])
                plot.save(plot.filename(f, key))

            plot = ScatterPlot(with_ghost_cells)
            plot(grid, u[key])
            plot.save(plot.filename(f, key))

        for key in log_keys:
            if n_dims == 2:
                plot = TriPlot()
                plot.color_plot(grid, np.log10(u[key]))
                plot.save(plot.filename(f, "log" + key))

            plot = ScatterSemilogY(with_ghost_cells)
            plot(grid, u[key])
            plot.save(plot.filename(f, "log" + key))


def plot_radial_average(grid, data_files, lin_keys, log_keys):
    keys = list(set(lin_keys + log_keys))
    for key in keys:
        n_bins = 50
        n_steps = len(data_files)

        r = np.linalg.norm(grid.cell_centers, axis=1)
        r_min, r_max = np.min(r), np.max(r)
        r_bins = np.linspace(r_min, r_max, n_bins + 1)
        indices = [
            np.argwhere(np.logical_and(a < r, r <= b))
            for a, b in zip(r_bins[:-1], r_bins[1:])
        ]

        t = np.empty(n_steps)
        avg_values = np.empty((n_bins, n_steps))
        for l, f in enumerate(data_files):
            u = load(f)
            values = u[key]
            t[l] = u.time

            for k, I in enumerate(indices):
                avg_values[k, l] = np.average(values[I])

        T, R = np.meshgrid(t, 0.5 * (r_bins[1:] + r_bins[:-1]))

        if key in lin_keys:
            fig = plt.figure()
            plt.contourf(T, R, avg_values, 100, cmap="jet")
            figname = os.path.dirname(data_files[0]) + f"radial_average-{key}.png"
            plt.colorbar()
            plt.savefig(figname, dpi=300)
            plt.close(fig)

            write_pickle(
                os.path.dirname(data_files[0]) + f"radial_average-{key}.pkl",
                {"R": R, "T": T, "Z": avg_values},
            )

        if key in log_keys:
            fig = plt.figure()
            levels = np.linspace(14.0, 19.0, 101)
            plt.contourf(T, R, np.log10(avg_values), levels, cmap="jet", extend="min")
            figname = os.path.dirname(data_files[0]) + f"radial_average-log10{key}.png"
            plt.colorbar()
            plt.savefig(figname, dpi=300)
            plt.close(fig)

            write_pickle(
                os.path.dirname(data_files[0]) + f"radial_average-log10{key}.pkl",
                {"R": R, "T": T, "Z": np.log10(avg_values)},
            )


if __name__ == "__main__":
    parser_help = "Visualize the output Zisa."
    parser = argparse.ArgumentParser(description=parser_help)

    parser.add_argument(
        "--grid",
        help="The filename of the grid.",
    )

    parser.add_argument(
        "--vars",
        nargs="*",
        type=str,
        default=["rho"],
        help="Plot these variables, e.g. 'rho', 'mv1', 'E', 'p', 'h', etc.",
    )

    parser.add_argument(
        "--delta-vars",
        nargs="*",
        type=str,
        default=["drho"],
        help="Plot these variables with symmetric colorbar, e.g. 'drho', 'dp', etc.",
    )

    parser.add_argument(
        "--log-vars",
        nargs="*",
        type=str,
        default=[],
        help="Plot these variables in a semilog-y plot, e.g. 'rho', 'E', 'p', 'h', etc.",
    )

    parser.add_argument(
        "--with-ghost-cells",
        action="store_true",
        help="The ghost-cells are only included in the plots if this flag is passed.",
    )

    parser.add_argument(
        "--radial-average",
        action="store_true",
        help="Plot the radial average of `vars` over time.",
    )

    parser.add_argument(
        "data_files",
        nargs="*",
        type=str,
        default=["."],
        help="Either the directory to visualize, or the datafiles.",
    )

    args = parser.parse_args()

    if os.path.isdir(args.data_files[0]):
        base_directory = args.data_files[0]
        files = find_data_files(base_directory)
    else:
        base_directory = os.path.dirname(args.data_files[0])
        files = args.data_files

    grid = load_grid(args.grid or find_grid(base_directory))
    if not getattr(sys, "ps1", sys.flags.interactive):
        if args.radial_average:
            plot_radial_average(grid, files, args.vars, args.log_vars)

        else:
            plot_all(
                grid,
                files,
                args.vars,
                args.delta_vars,
                args.log_vars,
                args.with_ghost_cells,
            )
