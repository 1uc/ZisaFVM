#! /usr/bin/env python3

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


def load(data_file):
    steady_state_file = find_steady_state_file(os.path.dirname(data_file))
    if not os.path.isfile(steady_state_file):
        steady_state_file = None

    return tiwaz.post_process.load_data(data_file, steady_state_file)


def load_array(h5_file, key):
    with h5py.File(h5_file, "r") as h5:
        return np.array(h5[key])


def transform(key, grid, u):
    x = grid.cell_centers
    xhat = x / np.linalg.norm(x, axis=1).reshape((-1, 1))

    rho = u["rho"]

    n_cells = x.shape[0]
    v = np.empty((n_cells, 3))
    v[:, 0] = u["mv1"] / rho
    v[:, 1] = u["mv2"] / rho
    v[:, 2] = u["mv3"] / rho

    vr = np.sum(v * xhat, axis=1)
    vh = np.linalg.norm(v - vr.reshape((-1, 1)) * xhat, axis=1)

    if key == "vr":
        return vr

    elif key == "vh":
        return vh

    elif key == "v":
        return np.linalg.norm(v, axis=1)

    else:
        raise Exception("Invalid key.")


def plot_all(grid, data_files, with_ghost_cells=False):
    n_dims = 2 if grid.vertex_indices.shape[1] == 3 else 3

    keys = ["v"]

    for f in data_files:
        u = load(f)

        vh = transform("vh", grid, u)
        with h5py.File(f[:-3] + "--velocity.h5", "w") as h5:
            h5["vh"] = vh

        # for key in keys:
        #     q = transform(key, grid, u)

        #     plot = ScatterPlot(with_ghost_cells)
        #     plot(grid, q)
        #     plot.save(plot.filename(f, key))


if __name__ == "__main__":
    parser_help = "Visualize the output Zisa."
    parser = argparse.ArgumentParser(description=parser_help)

    parser.add_argument(
        "--grid", help="The filename of the grid.",
    )

    parser.add_argument(
        "--with-ghost-cells",
        action="store_true",
        help="The ghost-cells are only included in the plots if this flag is passed.",
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
        plot_all(grid, files, args.with_ghost_cells)
