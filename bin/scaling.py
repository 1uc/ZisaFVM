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
from tiwaz.convergence_plots import write_convergence_plots
from tiwaz.scatter_plot import plot_visual_convergence
from tiwaz.gmsh import generate_circular_grids
from tiwaz.gmsh import decompose_grids
from tiwaz.gmsh import renumber_grids
from tiwaz.site_details import MPIHeuristics
from tiwaz.work_estimate import ZisaWorkEstimate
from tiwaz.queue_args import FixedMPIQueueArgs
from tiwaz.launch_params import folder_name
from tiwaz.post_process import load_data, load_grid
from tiwaz.post_process import find_data_files, find_last_data_file
from tiwaz.post_process import find_steady_state_file
from tiwaz.tri_plot import TriPlot


class ScalingExperiment(sc.Subsection):
    def __init__(self, amplitude, width, n_proc):
        super().__init__(
            {
                "name": "scaling_experiment",
                "initial_conditions": {"amplitude": amplitude, "width": width},
                "n_proc": n_proc,
            }
        )

    def short_id(self):
        return f"{self['name']}_{self['n_proc']}"


processors = [36, 72, 120]
amplitude = 0.1
width = 0.05

eos = sc.IdealGasEOS(gamma=2.0, r_gas=1.0)
gravity = sc.PolytropeGravity()
euler = sc.Euler(eos, gravity)

t_end = 0.01
time = sc.Time(t_end=t_end)
io = sc.IO(
    "hdf5", "gaussian_bump", n_snapshots=0, parallel_strategy="gathered", n_writers=4
)


def make_work_estimate():
    n0 = sc.read_n_cells(grid_name_hdf5(4))

    t0 = 2 * timedelta(seconds=t_end / 1e-1 * 60 * 96)
    b0 = 0.0
    o0 = 100 * 1e6

    return ZisaWorkEstimate(n0=n0, t0=t0, b0=b0, o0=o0)


def grid_name_stem(l):
    return "grids/scaling-{}".format(l)


def grid_name_geo(l):
    return grid_name_stem(l) + ".geo"


def grid_name_hdf5(l):
    return grid_name_stem(l) + ".msh.h5"


def grid_config_string(l):
    if parallelization["mode"] == "mpi":
        return grid_name_stem(l)

    return grid_name_hdf5(l)


parallelization = {"mode": "mpi"}

radius = 0.5
mesh_levels = [4]
lc_rel = {l: 0.1 * 0.5 ** l for l in mesh_levels}

coarse_grid_levels = [4]

coarse_grid_choices = {
    "grid": [sc.Grid(grid_config_string(l), l) for l in coarse_grid_levels]
}
coarse_grids = all_combinations(coarse_grid_choices)
reference_grid = sc.Grid(grid_config_string(7), 7)

independent_choices = {
    "euler": [euler],
    "io": [io],
    "time": [time],
    "parallelization": [parallelization],
    "boundary-condition": [sc.BoundaryCondition("frozen")],
    "debug": [{"global_indices": False, "stencils": False}],
}

dependent_choices_a = {
    "flux-bc": [sc.FluxBC("constant"), sc.FluxBC("isentropic")],
    "well-balancing": [sc.WellBalancing("constant"), sc.WellBalancing("isentropic")],
}

dependent_choices_b = {
    "reconstruction": [
        sc.Reconstruction("CWENO-AO", [1]),
        sc.Reconstruction(
            "CWENO-AO", [2, 2, 2, 2], overfit_factors=[3.0, 2.0, 2.0, 2.0]
        ),
        sc.Reconstruction("CWENO-AO", [3, 2, 2, 2]),
        sc.Reconstruction("CWENO-AO", [4, 2, 2, 2]),
        sc.Reconstruction("CWENO-AO", [5, 2, 2, 2]),
    ],
    "ode": [
        sc.ODE("ForwardEuler"),
        sc.ODE("SSP2"),
        sc.ODE("SSP3"),
        sc.ODE("RK4"),
        sc.ODE("Fehlberg"),
    ],
    "quadrature": [
        sc.Quadrature(1),
        sc.Quadrature(1),
        sc.Quadrature(2),
        sc.Quadrature(3),
        sc.Quadrature(4),
    ],
}

base_choices = all_combinations(independent_choices)
choices_a = pointwise_combinations(dependent_choices_a)
choices_b = pointwise_combinations(dependent_choices_b)
model_choices = base_choices.product(choices_a).product(choices_b)

coarse_runs_ = model_choices.product(coarse_grids)


def make_runs(n_proc):
    coarse_runs = coarse_runs_.product(
        [{"experiment": ScalingExperiment(amplitude, width, n_proc)}]
    )

    coarse_runs = [sc.Scheme(choice) for choice in coarse_runs]

    return coarse_runs, []


all_runs = [make_runs(n_proc) for n_proc in processors]


def post_process(coarse_runs, reference_run):
    pass


def generate_grids(cluster):
    generate_circular_grids(grid_name_geo, radius, lc_rel, mesh_levels, with_halo=True)
    renumber_grids(grid_name_hdf5, mesh_levels)
    decompose_grids(grid_name_hdf5, mesh_levels, {k: processors for k in mesh_levels})


def main():
    parser = default_cli_parser("'scaling' numerical experiment.")
    args = parser.parse_args()

    if args.generate_grids:
        generate_grids(args.cluster)

    if args.run:
        build_zisa()

        queue_args = FixedMPIQueueArgs(
            mem_per_core=2.0 * 1e9, wall_clock=timedelta(minutes=30)
        )

        for c, r in all_runs:
            launch_all(c, force=args.force, queue_args=queue_args)


if __name__ == "__main__":
    main()
