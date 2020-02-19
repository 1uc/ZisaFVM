#! /usr/bin/env python3

import os
import shutil
import glob

from datetime import timedelta

import matplotlib.pyplot as plt

import tiwaz
import tiwaz.scheme as sc

from tiwaz.launch_params import all_combinations, pointwise_combinations
from tiwaz.launch_params import build_zisa

from tiwaz.post_process import load_results

from tiwaz.cli_parser import default_cli_parser
from tiwaz.launch_job import launch_all
from tiwaz.latex_tables import write_convergence_table
from tiwaz.scatter_plot import plot_visual_convergence
from tiwaz.gmsh import generate_cube_grids, generate_circular_grids

from tiwaz.launch_params import folder_name
from tiwaz.post_process import load_data, load_grid
from tiwaz.post_process import find_data_files, find_last_data_file
from tiwaz.post_process import find_steady_state_file
from tiwaz.post_process import write_xdmf
from tiwaz.tri_plot import TriPlot

from tiwaz.queue_args import MPIQueueArgs
from tiwaz.work_estimate import ZisaFixedMemoryWorkEstimate


class ShockBubbleExperiment(sc.Subsection):
    def __init__(self, amplitude, width):
        super().__init__(
            {
                "name": "shock_bubble",
                "initial_conditions": {"amplitude": amplitude, "width": width},
            }
        )

    def short_id(self):
        amp = self["initial_conditions"]["amplitude"]
        return self["name"] + "_amp{:.2e}".format(amp)


amplitudes = [1e-2]
bump_width = 0.05

eos = sc.IdealGasEOS(gamma=2.0, r_gas=1.0)
gravity = sc.ConstantGravity()
euler = sc.Euler(eos, gravity)

t_end = 0.09
time = sc.Time(t_end=t_end)
io = sc.IO("hdf5", "shock_bubble", n_snapshots=1)


def grid_name_stem(l):
    return "grids/shock_bubble-{}".format(l)


def grid_name_geo(l):
    return grid_name_stem(l) + ".geo"


def grid_name_hdf5(l):
    return grid_name_stem(l) + ".msh.h5"


mesh_width = 1.0
mesh_levels = list(range(0, 6))
lc_rel = {l: 0.1 * 0.5 ** l for l in mesh_levels}

coarse_grid_levels = list(range(5, 6))
coarse_grid_names = [grid_name_hdf5(level) for level in coarse_grid_levels]

coarse_grid_choices = {
    "grid": [sc.Grid(grid_name_hdf5(l), l) for l in coarse_grid_levels]
}
coarse_grids = all_combinations(coarse_grid_choices)
reference_grid = sc.Grid(grid_name_hdf5(4), 4)

independent_choices = {
    "euler": [euler],
    "flux-bc": [sc.FluxBC("constant")],
    "well-balancing": [sc.WellBalancing("constant")],
    "io": [io],
    "time": [time],
    "parallelization": [{"mode": "mpi"}],
}

dependent_choices = {
    "reconstruction": [
        sc.Reconstruction(
            "CWENO-AO", [3, 2, 2, 2, 2], overfit_factors=[3.0, 2.0, 2.0, 2.0, 2.0]
        )
    ],
    "ode": [sc.ODE("SSP2")],
    "quadrature": [sc.Quadrature(2)],
}

reference_choices = {
    "euler": [euler],
    "flux-bc": [sc.FluxBC("constant")],
    "well-balancing": [sc.WellBalancing("constant")],
    "time": [time],
    "io": [io],
    "reconstruction": [sc.Reconstruction("CWENO-AO", [5, 2, 2, 2])],
    "ode": [sc.ODE("SSP3")],
    "quadrature": [sc.Quadrature(4)],
    "grid": [reference_grid],
    "reference": [sc.Reference("constant", coarse_grid_names)],
}

base_choices = all_combinations(independent_choices)
model_choices = base_choices.product(pointwise_combinations(dependent_choices))

coarse_runs_ = model_choices.product(coarse_grids)
reference_runs_ = all_combinations(reference_choices)


def make_runs(amplitude):
    coarse_runs = coarse_runs_.product(
        [{"experiment": ShockBubbleExperiment(amplitude, bump_width)}]
    )

    reference_runs = reference_runs_.product(
        [{"experiment": ShockBubbleExperiment(amplitude, bump_width)}]
    )

    coarse_runs = [sc.Scheme(choice) for choice in coarse_runs]
    reference_runs = [sc.Scheme(choice) for choice in reference_runs]

    return coarse_runs, reference_runs


all_runs = [make_runs(amp) for amp in amplitudes]


def visualize_snapshot(grid, data_files):
    for data_file in data_files[::10]:
        u = load_data(data_file, find_steady_state_file(coarse_dir))

        rho = u.cvars["rho"]
        vx = u.cvars["mv1"] / rho
        vy = u.cvars["mv2"] / rho

        trip = TriPlot()
        trip.color_plot(coarse_grid, rho)
        trip.quiver(coarse_grid, vx, vy)

        plt.show()


def post_process(coarse_runs):
    for coarse_run in coarse_runs:
        write_xdmf(coarse_run)


class TableLabels:
    def __call__(self, col):
        return " ".join(str(v) for v in col.values())


def generate_grids():
    generate_cube_grids(grid_name_geo, mesh_width, lc_rel, mesh_levels)


def make_work_estimate():
    n0 = 200_000
    t0 = 2 * timedelta(seconds=t_end / 1e-1 * 30 * 96)
    b0 = 0.0 * 1e9
    o0 = 2.0 * 1e9

    return ZisaFixedMemoryWorkEstimate(n0=n0, t0=t0, b0=b0, o0=o0)


def main():
    parser = default_cli_parser("'shock_bubble' numerical experiment.")
    args = parser.parse_args()

    if args.generate_grids:
        generate_grids()

    if args.run:
        build_zisa()

        t_min = timedelta(minutes=10)
        t_max = timedelta(hours=4)
        work_estimate = make_work_estimate()

        queue_args = MPIQueueArgs(work_estimate, t_min=t_min, t_max=t_max)
        # queue_args = dict()

        for c, r in all_runs:
            launch_all(c, force=args.force, queue_args=queue_args)

            if args.reference:
                launch_all(r, force=args.force, queue_args=queue_args)

    if args.post_process:
        for c, r in all_runs:
            post_process(c)


if __name__ == "__main__":
    main()
