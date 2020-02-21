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
from tiwaz.tri_plot import TriPlot


def load(data_file):
    steady_state_file = find_steady_state_file(".")
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


def plot_all(grid, data_files, keys):
    for f in data_files:
        u = load(f)

        for key in keys:
            plot = TriPlot()

            plot.color_plot(grid, u[key])
            plot.save(plot.filename(f, key))

            plot = ScatterPlot()
            plot(grid, u[key])
            plot.save(plot.filename(f, key))


if __name__ == "__main__":
    parser_help = "Visualize the output Zisa."
    parser = argparse.ArgumentParser(description=parser_help)

    parser.add_argument(
        "--grid", help="The filename of the grid.",
    )

    parser.add_argument(
        "--vars",
        nargs="*",
        type=str,
        default=["rho", "drho"],
        help="Plot these variables, e.g. 'rho', 'drho', 'mv1', 'E', 'p', 'h', etc.",
    )

    parser.add_argument(
        "data_files",
        nargs="+",
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
        plot_all(grid, files, args.vars)
