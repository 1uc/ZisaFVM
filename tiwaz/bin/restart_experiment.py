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
from tiwaz.gmsh import decompose_grids
from tiwaz.gmsh import renumber_grids
from tiwaz.gmsh import GridNamingScheme
from tiwaz.site_details import MPIHeuristics
from tiwaz.queue_args import MPIQueueArgs

from tiwaz.work_estimate import ZisaWorkEstimate


class RayleighTaylorExperiment(sc.Subsection):
    def __init__(self, drho):
        super().__init__(
            {
                "name": "rayleigh_taylor",
                "initial_conditions": {
                    "drho": drho,
                    "amplitude": 0.4,
                    "width": 0.1,
                    "n_bumps": 6,
                },
            }
        )

        self.drho = drho

    def short_id(self):
        return "restart_experiment"


G = 2.5
mangle_str = f"G{G:f}"

# This is just in case we want to run the same experiment with many parameters.
experiment_params = [0.1]

eos = sc.IdealGasEOS(gamma=2.0, r_gas=1.0)
gravity = sc.PolytropeGravityWithJump(rhoC=1.0, K_inner=1.0, K_outer=1.0, G=G)
euler = sc.Euler(eos, gravity)

t_end = 0.4
time = sc.Time(t_end=t_end)
io = sc.IO(
    "hdf5",
    "restart_experiment",
    n_snapshots=100,
    parallel_strategy="gathered",
    n_writers=4,
)

parallelization = {"mode": "mpi"}
radius = 0.6
mesh_levels = [2]
lc_rel = {l: 0.1 * 0.5 ** l for l in mesh_levels}
grid_name = GridNamingScheme("restart_experiment")

coarse_grid_levels = [2]
coarse_grid_names = [grid_name.msh_h5(level) for level in coarse_grid_levels]

coarse_grid_choices = {
    "grid": [
        sc.Grid(grid_name.config_string(l, parallelization), l)
        for l in coarse_grid_levels
    ]
}
coarse_grids = all_combinations(coarse_grid_choices)
reference_grid = sc.Grid(grid_name.config_string(2, parallelization), 2)


independent_choices = {
    "euler": [euler],
    "flux-bc": [sc.FluxBC("isentropic")],
    "boundary-condition": [sc.BoundaryCondition("frozen")],
    "well-balancing": [sc.WellBalancing("constant")],
    "io": [io],
    "time": [time],
    "parallelization": [parallelization],
    "debug": [{"global_indices": False, "stencils": False}],
}

dependent_choices = {
    "reconstruction": [
        sc.Reconstruction(
            "CWENO-AO", [3, 2, 2, 2], overfit_factors=[3.0, 3.0, 3.0, 3.0]
        ),
    ],
    "ode": [
        sc.ODE("SSP3"),
    ],
    "quadrature": [sc.Quadrature(2)],
}

base_choices = all_combinations(independent_choices)
model_choices = base_choices.product(pointwise_combinations(dependent_choices))

coarse_runs_ = model_choices.product(coarse_grids)


def make_runs(params):
    coarse_runs = coarse_runs_.product(
        [
            {
                "experiment": RayleighTaylorExperiment(params),
                "mangle": sc.Mangle(mangle_str),
            }
        ]
    )

    coarse_runs = [sc.Scheme(choice) for choice in coarse_runs]

    return coarse_runs, []


all_runs = [make_runs(param) for param in experiment_params]


def compute_parts(mesh_levels, host):
    work_estimate = make_work_estimate()

    heuristics = MPIHeuristics(host=host)
    queue_args = MPIQueueArgs(
        work_estimate, t_min=None, t_max=None, heuristics=heuristics
    )

    parts_ = dict()
    for l in mesh_levels:
        parts_[l] = [queue_args.n_mpi_tasks({"grid": {"file": grid_name.msh_h5(l)}})]

        if l <= 5:
            parts_[l].append(2)

    return parts_


def generate_grids(cluster, must_generate, must_decompose):
    geo_name = lambda l: grid_name.geo(l)
    msh_h5_name = lambda l: grid_name.msh_h5(l)

    if must_generate:
        for l in mesh_levels:
            shutil.rmtree(grid_name.dir(l), ignore_errors=True)

        generate_circular_grids(geo_name, radius, lc_rel, mesh_levels, with_halo=True)
        renumber_grids(msh_h5_name, mesh_levels)

    if must_decompose:
        decompose_grids(msh_h5_name, mesh_levels, compute_parts(mesh_levels, cluster))


def make_work_estimate():
    n0 = 200_000

    # 'measured' on Euler on L=4 with 96 cores.
    t0 = 1.5 * timedelta(seconds=t_end / 1e-1 * 60 * 96)

    # measured on Euler on L=4 with 2 and 96 cores.
    b0 = 0.0

    MB = 1e6
    o0 = 100 * MB

    return ZisaWorkEstimate(n0=n0, t0=t0, b0=b0, o0=o0)


def main():
    parser = default_cli_parser("'restart_experiment' numerical experiment.")
    args = parser.parse_args()

    generate_grids(args.cluster, args.generate_grids, args.decompose_grids)

    if args.run:
        build_zisa()

        t_min = timedelta(minutes=10)
        t_max = timedelta(hours=24)
        work_estimate = make_work_estimate()

        queue_args = MPIQueueArgs(work_estimate, t_min=t_min, t_max=t_max)

        for c, r in all_runs:
            launch_all(c, force=args.force, queue_args=queue_args)

            if args.reference:
                launch_all(r, force=args.force, queue_args=queue_args)


if __name__ == "__main__":
    main()
