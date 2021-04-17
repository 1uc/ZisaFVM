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
    def __init__(self, param):
        drho, v_amp = param
        super().__init__(
            {
                "name": "rayleigh_taylor",
                "initial_conditions": {
                    "r_crit": 0.25,
                    "drho": drho,
                    "amplitude": v_amp,
                    "A_noise": 0.0,
                    "A_shape": 0.0,
                    "r_noise": 0.4,
                    "width": 0.1,
                    "n_bumps": 6,
                },
            }
        )

        self.drho = drho
        self.v_amp = v_amp

    def short_id(self):
        return self["name"] + f"_drho{self.drho:.2e}_vamp{self.v_amp:.2e}"


G = 2.5
mangle_str = f"bcPT_G{G:f}"

experiment_params = [(drho, v_amp) for drho in [1e-1] for v_amp in [1e-4, 1e-2]]

eos = sc.IdealGasEOS(gamma=2.0, r_gas=1.0)
gravity = sc.PolytropeGravity(rhoC=1.0, K=1.0, G=G)
euler = sc.Euler(eos, gravity)

t_end = 10.0
time = sc.Time(t_end=t_end)
io = sc.IO(
    "hdf5",
    "rayleigh_taylor",
    n_snapshots=10,
    parallel_strategy="gathered",
    n_writers=4,
)

parallelization = {"mode": "mpi"}
radius = 0.6
mesh_levels = [6, 7, 8]
lc_rel = {l: 0.1 * 0.5 ** l for l in mesh_levels}
grid_name = GridNamingScheme("rayleigh_taylor_with_halo")

coarse_grid_levels = [6, 7, 8]
coarse_grid_names = [grid_name.msh_h5(level) for level in coarse_grid_levels]

coarse_grid_choices = {
    "grid": [
        sc.Grid(grid_name.config_string(l, parallelization), l)
        for l in coarse_grid_levels
    ]
}
coarse_grids = all_combinations(coarse_grid_choices)


independent_choices = {
    "euler": [euler],
    "flux-bc": [sc.FluxBC("none")],
    "boundary-condition": [sc.BoundaryCondition("frozen")],
    "well-balancing": [sc.WellBalancing("isentropic")],
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
    return coarse_runs, None


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


def generate_grids(cluster, must_generate, must_decompose):
    geo_name = lambda l: grid_name.geo(l)
    msh_h5_name = lambda l: grid_name.msh_h5(l)

    if must_generate:
        for l in mesh_levels:
            shutil.rmtree(grid_name.dir(l), ignore_errors=True)

        generate_circular_grids(geo_name, radius, lc_rel, mesh_levels, with_halo=True)
        renumber_grids(msh_h5_name, mesh_levels)

    # if must_decompose:
    #     decompose_grids(msh_h5_name, mesh_levels, compute_parts(mesh_levels, cluster))


def main():
    parser = default_cli_parser("'rayleigh_taylor' numerical experiment.")
    args = parser.parse_args()

    generate_grids(args.cluster, args.generate_grids, args.decompose_grids)

    if args.run:
        build_zisa()

        t_min = timedelta(minutes=10)
        t_max = timedelta(hours=24)
        work_estimate = make_work_estimate()

        queue_args = MPIQueueArgs(work_estimate, t_min=t_min, t_max=t_max)
        # queue_args = dict()

        for c, r in all_runs:
            launch_all(c, force=args.force, queue_args=queue_args)

    if args.post_process:
        for c, r in all_runs:
            post_process(c, r)


if __name__ == "__main__":
    main()
