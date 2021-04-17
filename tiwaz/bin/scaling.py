#! /usr/bin/env python3

import os
import shutil
import glob
import itertools
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
from tiwaz.gmsh import generate_spherical_grids
from tiwaz.gmsh import decompose_grids
from tiwaz.gmsh import renumber_grids
from tiwaz.gmsh import GridNamingScheme
from tiwaz.site_details import MPIHeuristics
from tiwaz.queue_args import TabularMPIQueueArgs
from tiwaz.launch_params import folder_name


class ScalingExperiment(sc.Subsection):
    def __init__(self, n_proc, n_threads):
        super().__init__(
            {
                "name": "scaling_experiment",
                "initial_conditions": {"amplitude": 0.01, "width": 0.05},
                "n_proc": n_proc,
                "n_threads": n_threads,
            }
        )

    def short_id(self):
        return f"{self['name']}_{self['n_proc']}_{self['n_threads']}"


mangle_str = f"noIO"

steps_per_frame = 1000
n_steps = 20
amplitude = 0.1
width = 0.05

eos = sc.IdealGasEOS(gamma=2.0, r_gas=1.0)
gravity = sc.PolytropeGravity()
euler = sc.Euler(eos, gravity)

grid_name = GridNamingScheme("scaling")
parallelization = {"mode": "mpi"}

radius = 0.5
# mesh_levels = range(5)
mesh_levels = [0, 1, 2, 3, 4]
# mesh_levels = [3]
coarse_grid_levels = [0, 1, 2, 3, 4]
# coarse_grid_levels = [1]
# coarse_grid_levels = [3]

lc_rel = {
    0: 0.1,
    1: 0.1 * 0.5 ** 1,
    2: 0.1 * 0.5 ** 2,
    3: 0.1 * 0.5 ** 3,
    4: 0.1 * 0.5 ** 4,
}

# fmt: off
pes = {
    # 2, 4, 8
    0: [
        (2, 1),
        (2, 2), (4, 1),
        (2, 4), (4, 2), (8, 1),
    ],
    # 16, 32, 64
    1: [
        (2,  8), (4, 4), ( 8, 2), (16, 1),
        (2, 16), (4, 8), ( 8, 4), (16, 2), (32, 1),
        (4, 16), (8, 8), (16, 4), (32, 2), (64, 1),
    ],
    # 128, 256, 512
    2: [
        ( 8, 16), (16, 8), ( 32, 4), ( 64, 2), (128, 1),
        (16, 16), (32, 8), ( 64, 4), (128, 2), (256, 1),
        (32, 16), (64, 8), (128, 4), (256, 2), (512, 1),
    ],
    # 1024, 2048, 4096
    3: [
        ( 32, 32), ( 64, 16), (128, 8), ( 256, 4), ( 512, 2),
        ( 64, 32), (128, 16), (256, 8), ( 512, 4), (1024, 2),
        (128, 32), (256, 16), (512, 8), (1024, 4), (2048, 2),
    ],

    # 8192
    4: [
        (512, 16), (1024, 8), (2048, 4)
    ],
}
# fmt: on
# processors[3] = [16 * 2 ** k for k in range(10)]
# processors[mesh_levels[-1]] = processors[mesh_levels[-1]][:1]


no_io_mpi_table = {
    "wall-clock": {
        l: {n: timedelta(minutes=30) for n in [2 ** k for k in range(1, 12)]}
        for l in mesh_levels
    },
    "mem-per-core": {l: dict() for l in mesh_levels},
}

# no_io_mpi_table["wall-clock"][3][16] = timedelta(minutes=120)
# no_io_mpi_table["wall-clock"][3][32] = timedelta(minutes=70)
# no_io_mpi_table["wall-clock"][3][64] = timedelta(minutes=55)
# no_io_mpi_table["wall-clock"][3][128] = timedelta(minutes=40)


for l in mesh_levels:
    no_io_mpi_table["mem-per-core"][l][2] = 200.0 * 1e9
    no_io_mpi_table["mem-per-core"][l][4] = 100.0 * 1e9
    no_io_mpi_table["mem-per-core"][l][8] = 62.0 * 1e9
    no_io_mpi_table["mem-per-core"][l][16] = 31.0 * 1e9
    no_io_mpi_table["mem-per-core"][l][32] = 15.0 * 1e9
    no_io_mpi_table["mem-per-core"][l][64] = 7.6 * 1e9
    no_io_mpi_table["mem-per-core"][l][128] = 3.8 * 1e9
    no_io_mpi_table["mem-per-core"][l][256] = 1.9 * 1e9
    no_io_mpi_table["mem-per-core"][l][512] = 1.9 * 1e9
    no_io_mpi_table["mem-per-core"][l][1024] = 1.9 * 1e9
    no_io_mpi_table["mem-per-core"][l][2048] = 1.9 * 1e9
    no_io_mpi_table["mem-per-core"][l][4096] = 1.9 * 1e9
    no_io_mpi_table["mem-per-core"][l][8192] = 1.9 * 1e9

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
                    ScalingExperiment(n_proc=n_proc, n_threads=n_threads)
                    for n_proc, n_threads in pes[l]
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


def load_array(f, key):
    with h5py.File(f, "r") as h5:
        return np.array(h5[key])


def post_process(coarse_runs, reference_run):
    results = []

    for cr in coarse_runs:
        runtime_file = os.path.join(cr.folder_name(), "run_time.txt")
        runtime = read_float(runtime_file)

        halo_files = glob.glob(
            os.path.join(cr.folder_name(), "mpi_exchange_runtime_wait-*.txt")
        )

        if not halo_files:
            continue

        halo_wait_cost = np.array([read_float(fn) for fn in halo_files])

        halo_files = glob.glob(
            os.path.join(cr.folder_name(), "mpi_exchange_runtime_post-*.txt")
        )
        halo_post_cost = np.array([read_float(fn) for fn in halo_files])

        local_cell_files = glob.glob(
            os.path.join(cr.folder_name(), "local_cells-*.txt")
        )
        local_cells = np.array([read_float(fn) for fn in local_cell_files])

        non_local_cell_files = glob.glob(
            os.path.join(cr.folder_name(), "non_local_cells-*.txt")
        )
        non_local_cells = np.array([read_float(fn) for fn in non_local_cell_files])

        gn = "../zisa-grids/" + os.path.basename(cr["grid"]["file"]) + "/grid.msh.h5"
        vi = load_array(gn, "vertex_indices")
        v = load_array(gn, "vertices")

        x = 1.0 / 3.0 * (v[vi[:, 0], :] + v[vi[:, 1], :] + v[vi[:, 2], :])
        r = np.linalg.norm(x, axis=1)
        n_cells = np.count_nonzero(r < 0.5)

        results.append(
            {
                "solver": {"wb": cr.well_balancing(), "order": cr.order()},
                "grid_level": cr.grid_level(),
                "n_cells": n_cells,
                "n_tasks": cr["experiment"]["n_proc"],
                "n_threads": cr["experiment"]["n_threads"],
                "n_cores": cr["experiment"]["n_proc"] * cr["experiment"]["n_threads"],
                "runtime": runtime,
                "halo_wait_cost": list(halo_wait_cost),
                "halo_post_cost": list(halo_post_cost),
                "local_cells": list(local_cells),
                "non_local_cells": list(non_local_cells),
            }
        )

    tiwaz.io.write_json(f"scaling_{mangle_str}.json", results)


def generate_grids(cluster, must_generate, must_decompose):
    geo_name = lambda l: grid_name.geo(l)
    msh_h5_name = lambda l: grid_name.msh_h5(l)

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
        queue_args.lsf_args = ["-R", '"model == EPYC_7742"']
        queue_args.is_hybrid = True

        c, _ = all_runs
        launch_all(c, force=args.force, queue_args=queue_args)


if __name__ == "__main__":
    main()
