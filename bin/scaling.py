#! /usr/bin/env python3

import os
import shutil
import glob
from datetime import timedelta

import numpy as np
import matplotlib.pyplot as plt

import tiwaz
import tiwaz.scheme as sc

from tiwaz.launch_params import all_combinations, pointwise_combinations
from tiwaz.launch_params import build_zisa

from tiwaz.cli_parser import default_cli_parser
from tiwaz.launch_job import launch_all
from tiwaz.gmsh import generate_spherical_grids
from tiwaz.gmsh import decompose_grids
from tiwaz.gmsh import renumber_grids
from tiwaz.gmsh import GridNamingScheme
from tiwaz.site_details import MPIHeuristics
from tiwaz.queue_args import TabularMPIQueueArgs
from tiwaz.launch_params import folder_name


class ScalingExperiment(sc.Subsection):
    def __init__(self, n_proc):
        super().__init__(
            {
                "name": "scaling_experiment",
                "initial_conditions": {"amplitude": 0.01, "width": 0.05},
                "n_proc": n_proc,
            }
        )

    def short_id(self):
        return f"{self['name']}_{self['n_proc']}"


mangle_str = f"noIO"

steps_per_frame = 1000
n_steps = 100
amplitude = 0.1
width = 0.05

eos = sc.IdealGasEOS(gamma=2.0, r_gas=1.0)
gravity = sc.PolytropeGravity()
euler = sc.Euler(eos, gravity)

grid_name = GridNamingScheme("scaling")
parallelization = {"mode": "mpi"}

radius = 0.5
# mesh_levels = range(5)
# mesh_levels = [0, 1, 2, 3]
mesh_levels = [3]
# coarse_grid_levels = [0, 1, 2]
coarse_grid_levels = [3]

lc_rel = {0: 0.1, 1: 0.042, 2: 0.042 * 0.44, 3: 0.042 * 0.44 * 0.48, 4: 0.00625}

processors = {l: [2 ** (3 * l) * 2 * 2 ** k for k in range(3)] for l in mesh_levels}
# processors[3] = [16 * 2 ** k for k in range(10)]
# processors[mesh_levels[-1]] = processors[mesh_levels[-1]][:1]


no_io_mpi_table = {
    "wall-clock": {
        l: {n: timedelta(minutes=30) for n in processors[l]} for l in processors
    },
    "mem-per-core": {l: {n: 1.0 * 1e9 for n in processors[l]} for l in processors},
}

no_io_mpi_table["wall-clock"][3][16] = timedelta(minutes=120)
no_io_mpi_table["wall-clock"][3][32] = timedelta(minutes=70)
no_io_mpi_table["wall-clock"][3][64] = timedelta(minutes=55)
no_io_mpi_table["wall-clock"][3][128] = timedelta(minutes=40)

no_io_mpi_table["mem-per-core"][3][16] = 32.0 * 1e9
no_io_mpi_table["mem-per-core"][3][32] = 16.0 * 1e9
no_io_mpi_table["mem-per-core"][3][64] = 8.0 * 1e9
no_io_mpi_table["mem-per-core"][3][128] = 4.0 * 1e9
no_io_mpi_table["mem-per-core"][3][256] = 2.0 * 1e9

mpi_table = no_io_mpi_table

independent_choices = {
    "euler": [euler],
    "mangle": [sc.Mangle(mangle_str)],
    "io": [
        sc.IO(
            "hdf5",
            f"scaling",
            steps_per_frame=steps_per_frame,
            parallel_strategy="gathered",
            n_writers=4,
        )
    ],
    "time": [sc.Time(n_steps=n_steps)],
    "parallelization": [parallelization],
    "boundary-condition": [sc.BoundaryCondition("frozen")],
    "debug": [{"global_indices": False, "stencils": False}],
}

dependent_choices_a = {
    # "flux-bc": [sc.FluxBC("none"), sc.FluxBC("none")],
    # "well-balancing": [sc.WellBalancing("constant"), sc.WellBalancing("isentropic")],
    "flux-bc": [sc.FluxBC("none")],
    "well-balancing": [sc.WellBalancing("isentropic")],
}

dependent_choices_b = {
    "reconstruction": [sc.Reconstruction("CWENO-AO", [3, 2, 2, 2, 2])],
    "ode": [sc.ODE("SSP3")],
    "quadrature": [sc.Quadrature(2)],
}

base_choices = all_combinations(independent_choices)
choices_a = pointwise_combinations(dependent_choices_a)
choices_b = pointwise_combinations(dependent_choices_b)
model_choices = base_choices.product(choices_a).product(choices_b)


def make_runs():
    scaling_choices = []
    for l in coarse_grid_levels:
        experiment_choices = all_combinations(
            {
                "experiment": [
                    ScalingExperiment(n_proc=n_proc) for n_proc in processors[l]
                ]
            }
        )

        grid_choices = [
            {"grid": sc.Grid(grid_name.config_string(l, parallelization), l)}
        ]

        scaling_choices += experiment_choices.product(grid_choices)

    coarse_runs = model_choices.product(scaling_choices)
    coarse_runs = [sc.Scheme(choice) for choice in coarse_runs]

    return coarse_runs, []


all_runs = make_runs()


def read_float(fn):
    if os.path.exists(fn):
        return float(tiwaz.io.read_txt(fn))
    else:
        print(fn)
        return -1.0


def post_process(coarse_runs, reference_run):
    results = []

    for cr in coarse_runs:
        runtime_file = os.path.join(cr.folder_name(), "run_time.txt")
        runtime = read_float(runtime_file)

        halo_files = glob.glob(
            os.path.join(cr.folder_name(), "mpi_exchange_runtime-*.txt")
        )
        halo_cost = np.array([read_float(fn) for fn in halo_files])

        local_cell_files = glob.glob(
            os.path.join(cr.folder_name(), "local_cells-*.txt")
        )
        local_cells = np.array([read_float(fn) for fn in local_cell_files])

        non_local_cell_files = glob.glob(
            os.path.join(cr.folder_name(), "non_local_cells-*.txt")
        )
        non_local_cells = np.array([read_float(fn) for fn in non_local_cell_files])

        results.append(
            {
                "solver": {"wb": cr.well_balancing(), "order": cr.order()},
                "grid_level": cr.grid_level(),
                "n_cells": sc.read_n_cells(
                    "../zisa-grids/" + os.path.basename(cr["grid"]["file"]) + "/grid.h5"
                ),
                "n_proc": cr["experiment"]["n_proc"],
                "runtime": runtime,
                "halo_cost": list(halo_cost),
                "local_cells": list(local_cells),
                "non_local_cells": list(non_local_cells),
                "avg_halo_cost": np.mean(halo_cost),
                "min_halo_cost": np.min(halo_cost),
                "max_halo_cost": np.max(halo_cost),
            }
        )

    tiwaz.io.write_json(f"scaling_{mangle_str}.json", results)


def generate_grids(cluster, must_generate, must_decompose):
    geo_name = lambda l: grid_name.geo(l)
    msh_h5_name = lambda l: grid_name.msh_h5(l)

    print(processors)

    if must_generate:
        for l in mesh_levels:
            shutil.rmtree(grid_name.dir(l), ignore_errors=True)

        generate_spherical_grids(geo_name, radius, lc_rel, mesh_levels, with_halo=True)
        renumber_grids(msh_h5_name, mesh_levels)

    # if must_decompose:
    #     decompose_grids(msh_h5_name, mesh_levels, processors)


def main():
    parser = default_cli_parser("'scaling' numerical experiment.")
    args = parser.parse_args()

    if args.generate_grids or args.decompose_grids:
        generate_grids(args.cluster, args.generate_grids, args.decompose_grids)

    if args.post_process:
        post_process(*all_runs)

    if args.run:
        build_zisa()

        queue_args = TabularMPIQueueArgs(mpi_table)
        queue_args.lsf_args = ["-R", '"model == EPYC_7742"', "-R", "span[ptile=128]"]

        c, _ = all_runs
        launch_all(c, force=args.force, queue_args=queue_args)


if __name__ == "__main__":
    main()
