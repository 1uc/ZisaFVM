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


G = 4.0
mangle_str = f"PT_G{G:f}"

# This is just in case we want to run the same experiment with many parameters.
# experiment_params = [0.0, 1e-4]
# experiment_params = [1e-4 * 3 ** k for k in range(1, 5)]
# experiment_params = [0.0]
experiment_params = [
    (drho, v_amp) for drho in [1e-4, 1e-3, 1e-2] for v_amp in [1e-6, 1e-4, 1e-2]
]

eos = sc.IdealGasEOS(gamma=2.0, r_gas=1.0)
gravity = sc.PolytropeGravity(rhoC=1.0, K=1.0, G=G)
# gravity = sc.ConstantGravity(g=0.0)
euler = sc.Euler(eos, gravity)

t_end = 100.0
time = sc.Time(t_end=t_end)
# io = sc.IO(
#     "hdf5", "rayleigh_taylor", n_snapshots=20, parallel_strategy="gathered", n_writers=4
# )
io = sc.IO(
    "hdf5",
    "rayleigh_taylor",
    n_snapshots=100,
    parallel_strategy="gathered",
    n_writers=4,
)
# io = sc.IO("opengl", "rayleigh_taylor", steps_per_frame=2)

parallelization = {"mode": "mpi"}
radius = 0.6
# mesh_levels = list(range(1, 8))
# mesh_levels = [2, 3, 4]
mesh_levels = [1, 2, 3, 4]
lc_rel = {l: 0.1 * 0.5 ** l for l in mesh_levels}
grid_name = GridNamingScheme("rayleigh_taylor_with_halo")

coarse_grid_levels = [4]
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
    "well-balancing": [sc.WellBalancing("constant"), sc.WellBalancing("isentropic")],
    "io": [io],
    "time": [time],
    "parallelization": [parallelization],
    "debug": [{"global_indices": False, "stencils": False}],
}

dependent_choices = {
    "reconstruction": [
        # sc.Reconstruction("CWENO-AO", [1]),
        sc.Reconstruction(
            "CWENO-AO", [3, 2, 2, 2], overfit_factors=[3.0, 3.0, 3.0, 3.0]
        ),
        # sc.Reconstruction("CWENO-AO", [1], overfit_factors=[1.0],),
    ],
    "ode": [
        # sc.ODE("ForwardEuler"),
        sc.ODE("SSP3"),
        # sc.ODE("SSP3", cfl_number=0.4)
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
    "parallelization": [parallelization],
}

base_choices = all_combinations(independent_choices)
model_choices = base_choices.product(pointwise_combinations(dependent_choices))

coarse_runs_ = model_choices.product(coarse_grids)
reference_runs_ = all_combinations(reference_choices)


def make_runs(params):
    coarse_runs = coarse_runs_.product(
        [
            {
                "experiment": RayleighTaylorExperiment(params),
                "mangle": sc.Mangle(mangle_str),
            }
        ]
    )

    reference_runs = reference_runs_.product(
        [
            {
                "experiment": RayleighTaylorExperiment(params),
                "mangle": sc.Mangle(mangle_str),
            }
        ]
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
    n0 = sc.read_n_cells(grid_name.msh_h5(4))

    # 'measured' on Euler on L=4 with 96 cores.
    t0 = 1.5 * timedelta(seconds=t_end / 1e-1 * 60 * 96)

    # measured on Euler on L=4 with 2 and 96 cores.
    b0 = 0.0

    MB = 1e6
    o0 = 100 * MB

    return ZisaWorkEstimate(n0=n0, t0=t0, b0=b0, o0=o0)


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
