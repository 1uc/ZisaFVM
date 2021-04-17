#!/usr/bin/env python3

import os
import sys
import glob
import itertools

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
        tex.append(f)
        u = load(f)

        # vr_plot = plot_radial_velocity(grid, u, with_ghost_cells)
        # vr_plot.save(vr_plot.filename(f, "vr"))

        # vh_plot = plot_horizontal_velocity(grid, u, with_ghost_cells)
        # vh_plot.save(vh_plot.filename(f, "vh"))

        mu_plot = plot_massfractions(grid, u, with_ghost_cells)
        mu_plot.save(mu_plot.filename(f, "mu"))


def plot_error_growth_rate():
    levels = range(2, 5)
    dirs = [
        "rayleigh_taylor_atol1e-13_0.00e+00_G2.500000_o3222_SSP3_isentropic_L2",
        "rayleigh_taylor_atol1e-13_0.00e+00_G2.500000_o3222_SSP3_isentropic_L3",
        "rayleigh_taylor_atol1e-13_0.00e+00_G2.500000_o3222_SSP3_isentropic_L4",
    ]

    r_min, r_max = 0.35, 0.45

    linf_errors = dict()

    for l, d in zip(levels, dirs):
        files = sorted(glob.glob(f"{d}/*_data-????.h5"))
        gridname = f"{d}/grids/rayleigh_taylor_with_halo-{l}/grid.h5"
        grid = load_grid(gridname)

        rho_eq = load_array(f"{d}/steady_state.h5", "rho")
        r = np.linalg.norm(grid.cell_centers, axis=1)
        I = np.logical_and(r_min < r, r < r_max)

        linf_errors[l] = np.empty(len(files))
        for k, f in enumerate(files):
            rho = load_array(f, "rho")
            drho = rho - rho_eq

            linf_errors[l][k] = np.max(np.abs(drho[I]))

    tex.append(linf_errors)

    for l in linf_errors:
        plt.semilogy(linf_errors[l], label=f"L = {l}")

    plt.legend()
    plt.savefig(dirs[0][:-3] + "_linf.png")
    plt.show()


def plot_rayleigh_comparison():

    amps = [1e-6, 1e-4]
    drhos = [1e-4]

    levels = []
    dirs = []

    for l in range(4, 5):
        for a, d in itertools.product(amps, drhos):
            levels.append(l)
            dirs.append(
                (
                    f"rayleigh_taylor_drho{d:.2e}_vamp{a:.2e}_PT_G4.000000_o3222_SSP3_constant_L{l}",
                    f"rayleigh_taylor_drho{d:.2e}_vamp{a:.2e}_PT_G4.000000_o3222_SSP3_isentropic_L{l}",
                )
            )

    for l, (dc, di) in zip(levels, dirs):
        const_files = sorted(glob.glob(f"{dc}/*_data-????.h5"))
        isen_files = sorted(glob.glob(f"{di}/*_data-????.h5"))
        subfiles = slice(0, None, 10)
        # subfiles = slice(0, 1, 10)

        gridname = f"{dc}/grids/rayleigh_taylor_with_halo-{l}/grid.h5"
        grid = load_grid(gridname)

        rho_eq = load_array(f"{di}/steady_state.h5", "rho")
        x = grid.cell_centers
        cI = x[:, 0] < 0
        iI = np.logical_not(cI)

        for (fc, fi) in zip(const_files[subfiles], isen_files[subfiles]):
            crho = load_array(fc, "rho")
            cdrho = crho - rho_eq

            irho = load_array(fi, "rho")
            idrho = irho - rho_eq

            drho = np.empty_like(idrho)
            drho[cI] = cdrho[cI]
            drho[iI] = 3.0 * idrho[iI]

            mask = np.linalg.norm(x, axis=1) > 0.6
            plot = TriPlot(cmap="PRGn", symmetric_limits=True)
            plot.color_plot(grid, drho, mask=mask)

            ax = plt.gca()

            t = load_array(fc, "time")[()]
            plt.title(f"t = {t:.1f}")
            plt.axvline(0.0, color="black", linewidth=2.0)
            plt.text(0.05, 0.95, "unbalanced", transform=ax.transAxes)
            plt.text(0.7, 0.95, "well-balanced", transform=ax.transAxes)
            plt.text(0.85, 0.05, "x3", transform=ax.transAxes)
            plot.save(plot.filename(fi, "comparison_drho"))


def generate_rayleigh_comparison_latex():
    amps = [1e-6, 1e-4, 1e-2]
    drhos = [1e-4, 1e-3, 1e-2]

    plot = TriPlot()

    tex = []
    tex.append(r"\documentclass[11pt,a4paper]{article}")
    tex.append(r"\usepackage{graphicx}")
    tex.append(r"\usepackage{grffile}")
    tex.append(r"\usepackage{verbatim}")
    tex.append(r"\begin{document}")

    for l in range(4, 5):
        subfiles = range(0, 100, 10)
        for k in subfiles:
            tex.append(r"\begin{figure}[htpb]")
            tex.append(r"\begin{center}")
            for a in amps:
                for d in drhos:
                    di = f"rayleigh_taylor_drho{d:.2e}_vamp{a:.2e}_PT_G4.000000_o3222_SSP3_isentropic_L{l}"
                    fi = f"{di}/rayleigh_taylor_data-{k:04d}.h5"
                    if not os.path.exists(fi):
                        fi = f"{di}/rayleigh_taylor_data-0000.h5"
                    fig_name = plot.filename(fi, "comparison_drho")

                    tex.append(
                        r"\includegraphics[width=0.3\textwidth]{" + fig_name + "}"
                    )

                tex.append(r"\\")

            tex.append(r"\end{center}")
            tex.append(r"\caption{" + f"L = {l}" + r" $k = $" + f"{k}" + "}")
            tex.append(r"\end{figure}")

    tex.append(r"\end{document}")

    tex = "\n".join(tex)
    with open("report.tex", "w") as f:
        f.write(tex)


if __name__ == "__main__":
    plot_rayleigh_comparison()
    # generate_rayleigh_comparison_latex()
