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


def compute_velocity(grid, u):
    v = np.empty(grid.cell_centers.shape)
    rho = u["rho"]
    v[:, 0] = u["mv1"] / rho
    v[:, 1] = u["mv2"] / rho
    v[:, 2] = u["mv3"] / rho

    return v


def compute_radial_velocity(grid, u):
    x = grid.cell_centers
    rho = u["rho"]
    v = compute_velocity(grid, u)

    vr = np.sum(x * v, axis=1) / np.linalg.norm(x, axis=1)
    return vr


def plot_radial_velocity(grid, u, with_ghost_cells):
    vr = compute_radial_velocity(grid, u)

    plot = ScatterPlot(with_ghost_cells)
    plot(grid, vr)
    plot(grid, u.xvars["cs"])
    return plot


def plot_horizontal_velocity(grid, u, with_ghost_cells):
    x = grid.cell_centers
    vr = compute_radial_velocity(grid, u)
    v = compute_velocity(grid, u)
    vh = np.linalg.norm(
        v - vr.reshape((-1, 1)) * x / np.linalg.norm(x, axis=1, keepdims=True), axis=1
    )

    plot = ScatterPlot(with_ghost_cells)
    plot(grid, vh)
    plot(grid, u.xvars["cs"])
    return plot


def plot_massfractions(grid, u, with_ghost_cells):
    rho = u["rho"]
    mu = [u[f"mq{k}"] / rho for k in range(6)]

    plot = ScatterSemilogY(with_ghost_cells)
    for k in range(6):
        plot(grid, mu[k])
    return plot


def plot_all(grid, data_files, with_ghost_cells=False):
    n_dims = 2 if grid.vertex_indices.shape[1] == 3 else 3

    for f in data_files:
        print(f)
        u = load(f)

        # vr_plot = plot_radial_velocity(grid, u, with_ghost_cells)
        # vr_plot.save(vr_plot.filename(f, "vr"))

        # vh_plot = plot_horizontal_velocity(grid, u, with_ghost_cells)
        # vh_plot.save(vh_plot.filename(f, "vh"))

        mu_plot = plot_massfractions(grid, u, with_ghost_cells)
        mu_plot.save(mu_plot.filename(f, "mu"))


if __name__ == "__main__":
    parser_help = "Visualize the output Zisa."
    parser = argparse.ArgumentParser(description=parser_help)

    parser.add_argument(
        "--grid",
        help="The filename of the grid.",
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
