#! /usr/bin/env python3

import os
import shutil
import glob
import h5py

from datetime import timedelta

import numpy as np

import matplotlib.pyplot as plt

import tiwaz
import tiwaz.scheme as sc

from tiwaz.launch_params import all_combinations, pointwise_combinations
from tiwaz.launch_params import build_zisa

from tiwaz.cli_parser import default_cli_parser
from tiwaz.launch_job import launch_all
from tiwaz.launch_params import folder_name
from tiwaz.post_process import load_data, load_grid
from tiwaz.post_process import find_data_files, find_last_data_file
from tiwaz.post_process import find_steady_state_file
from tiwaz.tri_plot import tri_plot
from tiwaz.scatter_plot import scatter_plot
from tiwaz.gmsh import generate_circular_grids
from tiwaz.queue_args import MPIQueueArgs

from tiwaz.work_estimate import ZisaWorkEstimate


class RayleighTaylorExperiment(sc.Subsection):
    def __init__(self):
        super().__init__(
            {
                "name": "rayleigh_taylor",
                "initial_conditions": {
                    "drho": 0.1,
                    "amplitude": 0.0001,
                    "width": 0.1,
                    "n_bumps": 6,
                },
            }
        )

    def short_id(self):
        return self["name"]


# This is just in case we want to run the same experiment with many parameters.
experiment_params = [None]

eos = sc.IdealGasEOS(gamma=2.0, r_gas=1.0)
gravity = sc.PolytropeGravityWithJump(rhoC=1.0, K_inner=1.0, K_outer=1.0, G=3.0)
euler = sc.Euler(eos, gravity)

t_end = 5.0
time = sc.Time(t_end=t_end)
io = sc.IO("hdf5", "rayleigh_taylor", n_snapshots=20)
# io = sc.IO("opengl", "rayleigh_taylor", steps_per_frame=2)

radius = 0.6
mesh_levels = list(range(1, 6))
lc_rel = {l: 0.1 * 0.5 ** l for l in mesh_levels}


def grid_name_stem(l):
    return "grids/rayleigh_taylor-{}".format(l)


def grid_name_geo(l):
    return grid_name_stem(l) + ".geo"


def grid_name_msh(l):
    return grid_name_stem(l) + ".msh.h5"


def grid_name_hdf5(l):
    return grid_name_stem(l) + ".msh.h5"


def generate_grids():
    generate_circular_grids(grid_name_geo, radius, lc_rel, mesh_levels)


def make_work_estimate():
    n0 = sc.read_n_cells(grid_name_hdf5(4))

    # measured on Euler on L=4 with 96 cores.
    t0 = 2 * timedelta(seconds=t_end / 1e-1 * 30 * 96)

    # measured on Euler on L=4 with 2 and 96 cores.
    b0 = 0.0 * 1e9
    o0 = 2.0 * 1e9

    return ZisaWorkEstimate(n0=n0, t0=t0, b0=b0, o0=o0)


coarse_grid_levels = list(range(2, 5))
coarse_grid_names = [grid_name_msh(level) for level in coarse_grid_levels]

coarse_grid_choices = {
    "grid": [sc.Grid(grid_name_msh(l), l) for l in coarse_grid_levels]
}
coarse_grids = all_combinations(coarse_grid_choices)
reference_grid = sc.Grid(grid_name_msh(2), 2)

independent_choices = {
    "euler": [euler],
    "flux-bc": [sc.FluxBC("isentropic")],
    "well-balancing": [sc.WellBalancing("isentropic")],
    "io": [io],
    "time": [time],
    "parallelization": [{"mode": "mpi"}],
}

dependent_choices = {
    "reconstruction": [
        # sc.Reconstruction("CWENO-AO", [1]),
        sc.Reconstruction(
            "CWENO-AO", [3, 2, 2, 2], overfit_factors=[3.0, 2.0, 2.0, 2.0]
        ),
    ],
    "ode": [
        # sc.ODE("ForwardEuler"),
        # sc.ODE("SSP3")
        sc.ODE("SSP3", cfl_number=0.8)
    ],
    "quadrature": [
        # sc.Quadrature(1),
        sc.Quadrature(2)
        # sc.Quadrature(3)
    ],
}

# dependent_choices = {
#     "reconstruction": [
#         sc.Reconstruction("CWENO-AO", [1]),
#         sc.Reconstruction("CWENO-AO", [2, 2, 2, 2], overfit_factors=[3.0, 2.0, 2.0, 2.0]),
#         sc.Reconstruction("CWENO-AO", [3, 2, 2, 2]),
#         sc.Reconstruction("CWENO-AO", [4, 2, 2, 2])
#     ],
#
#     "ode": [
#         sc.ODE("ForwardEuler"),
#         sc.ODE("SSP2"),
#         sc.ODE("SSP3"),
#         sc.ODE("SSP3")
#     ],
#
#     "quadrature": [
#         sc.Quadrature(1),
#         sc.Quadrature(2),
#         sc.Quadrature(3),
#         sc.Quadrature(4)
#     ]
# }

reference_choices = {
    "euler": [euler],
    "flux-bc": [sc.FluxBC("isentropic")],
    "well-balancing": [sc.WellBalancing("isentropic")],
    "time": [time],
    "io": [io],
    "reconstruction": [sc.Reconstruction("CWENO-AO", [5, 2, 2, 2])],
    "ode": [sc.ODE("SSP3")],
    "quadrature": [sc.Quadrature(4)],
    "grid": [reference_grid],
    "reference": [sc.Reference("isentropic", coarse_grid_names)],
    "parallelization": [{"mode": "mpi"}],
}

base_choices = all_combinations(independent_choices)
model_choices = base_choices.product(pointwise_combinations(dependent_choices))

coarse_runs_ = model_choices.product(coarse_grids)
reference_runs_ = all_combinations(reference_choices)


def make_runs(params):
    coarse_runs = coarse_runs_.product([{"experiment": RayleighTaylorExperiment()}])

    reference_runs = reference_runs_.product(
        [{"experiment": RayleighTaylorExperiment()}]
    )

    coarse_runs = [sc.Scheme(choice) for choice in coarse_runs]
    reference_runs = [sc.Scheme(choice) for choice in reference_runs]

    return coarse_runs, reference_runs


all_runs = [make_runs(param) for param in experiment_params]


def post_process(coarse_runs, reference_run):
    for coarse_run in coarse_runs:
        coarse_dir = folder_name(coarse_run)
        coarse_grid = load_grid(coarse_dir)
        data_files = find_data_files(coarse_dir)

        for filename in data_files:
            u_coarse = load_data(filename, find_steady_state_file(coarse_dir))

            rho = u_coarse.cvars["rho"]
            drho = u_coarse.dvars["rho"]

            plot = tri_plot(coarse_grid, rho)
            plot.save(filename[:-3] + "-rho.png")

            plot = tri_plot(coarse_grid, drho)
            plot.save(filename[:-3] + "-drho.png")

            plot = scatter_plot(coarse_grid, rho)
            plot.save(filename[:-3] + "-scatter-rho.png")

            plot = scatter_plot(coarse_grid, drho)
            plot.save(filename[:-3] + "-scatter-drho.png")


class TableLabels:
    def __call__(self, col):
        return " ".join(str(v) for v in col.values())


def main():
    parser = default_cli_parser("'rayleigh_taylor' numerical experiment.")
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
            post_process(c, r)

    if args.copy_to_paper:
        dir = "${HOME}/git/papers/LucGrosheintz/papers/unstructured_well_balancing/img/rayleigh_taylor"
        dir = os.path.expandvars(dir)

        for c, _ in all_runs:
            stem = c[0]["experiment"].short_id()
            files = glob.glob(stem + "*.tex")

            for f in files:
                shutil.copy(f, os.path.join(dir, os.path.basename(f)))


if __name__ == "__main__":
    main()
