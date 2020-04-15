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
from tiwaz.gmsh import generate_spherical_grids
from tiwaz.work_estimate import ZisaWorkEstimate
from tiwaz.queue_args import MPIQueueArgs
from tiwaz.launch_params import folder_name
from tiwaz.post_process import load_data, load_grid
from tiwaz.post_process import find_data_files, find_last_data_file
from tiwaz.post_process import find_steady_state_file
from tiwaz.scatter_plot import ScatterPlot


class GaussianBumpExperiment(sc.Subsection):
    def __init__(self, amplitude, width):
        super().__init__(
            {
                "name": "gaussian_bump_3d",
                "initial_conditions": {"amplitude": amplitude, "width": width},
            }
        )

    def short_id(self):
        amp = self["initial_conditions"]["amplitude"]
        return self["name"] + "_amp{:.2e}".format(amp)


amplitudes = [0.1]
width = 0.05

eos = sc.IdealGasEOS(gamma=2.0, r_gas=1.0)
gravity = sc.PolytropeGravity()
# gravity = sc.NoGravity()
# gravity = sc.ConstantGravity(g=1.0)
euler = sc.Euler(eos, gravity)

t_end = 0.09
time = sc.Time(t_end=t_end)
io = sc.IO("hdf5", "gaussian_bump", n_snapshots=1)
# io = sc.IO("opengl", "gaussian_bump", steps_per_frame=1)


def make_work_estimate():
    n0 = sc.read_n_cells(grid_name_hdf5(1))

    t0 = 2 * timedelta(seconds=t_end / 0.09 * 90 * 24)

    b0 = 0.5 * 1e9
    o0 = 1.0 * 1e9

    return ZisaWorkEstimate(n0=n0, t0=t0, b0=b0, o0=o0, n_dims=3)


def grid_name_stem(l):
    return "grids/gaussian_bump_3d-{}".format(l)


def grid_name_geo(l):
    return grid_name_stem(l) + ".geo"


def grid_name_hdf5(l):
    return grid_name_stem(l) + ".msh.h5"


radius = 0.5
mesh_levels = list(range(0, 5))
lc_rel = {l: 0.1 * 0.5 ** l for l in mesh_levels}

coarse_grid_levels = list(range(0, 3))
coarse_grid_names = [grid_name_hdf5(level) for level in coarse_grid_levels]

coarse_grid_choices = {
    "grid": [sc.Grid(grid_name_hdf5(l), l) for l in coarse_grid_levels]
}
coarse_grids = all_combinations(coarse_grid_choices)
reference_grid = sc.Grid(grid_name_hdf5(3), 3)

independent_choices = {
    "euler": [euler],
    "io": [io],
    "boundary-condition": [sc.BoundaryCondition("frozen")],
    "time": [time],
    "parallelization": [{"mode": "mpi"}],
    "debug": [{"global_indices": False, "stencils": False}],
}

dependent_choices_a = {
    "flux-bc": [
        sc.FluxBC("constant"),
        sc.FluxBC("isentropic")
    ],
    "well-balancing": [
        sc.WellBalancing("constant"),
        sc.WellBalancing("isentropic")
    ],
}

dependent_choices_b = {
    "reconstruction": [
        sc.Reconstruction("CWENO-AO", [1]),
        sc.Reconstruction(
            "CWENO-AO", [2, 2, 2, 2, 2], overfit_factors=[3.0, 2.5, 2.5, 2.5, 2.5]
        ),
        sc.Reconstruction(
            "CWENO-AO", [3, 2, 2, 2, 2], overfit_factors=[3.0, 2.5, 2.5, 2.5, 2.5]
        ),
        sc.Reconstruction(
            "CWENO-AO", [4, 2, 2, 2, 2], overfit_factors=[3.0, 2.5, 2.5, 2.5, 2.5]
        ),
    ],
    "ode": [
        sc.ODE("ForwardEuler"),
        sc.ODE("SSP2"),
        sc.ODE("SSP3"),
        sc.ODE("Fehlberg")
    ],
    "quadrature": [
        sc.Quadrature(1),
        sc.Quadrature(1),
        sc.Quadrature(2),
        sc.Quadrature(3),
    ],
}

reference_choices = {
    "euler": [euler],
    "flux-bc": [sc.FluxBC("isentropic")],
    "well-balancing": [sc.WellBalancing("isentropic")],
    "time": [time],
    "io": [io],
    "reconstruction": [sc.Reconstruction("CWENO-AO", [4, 2, 2, 2])],
    "ode": [sc.ODE("Fehlberg")],
    "quadrature": [sc.Quadrature(4)],
    "grid": [reference_grid],
    "reference": [sc.Reference("isentropic", coarse_grid_names)],
    "parallelization": [{"mode": "mpi"}],
    "debug": [{"global_indices": False, "stencils": False}],
}

base_choices = all_combinations(independent_choices)
choices_a = pointwise_combinations(dependent_choices_a)
choices_b = pointwise_combinations(dependent_choices_b)
model_choices = base_choices.product(choices_a).product(choices_b)

coarse_runs_ = model_choices.product(coarse_grids)
reference_runs_ = all_combinations(reference_choices)


def make_runs(amplitude):
    coarse_runs = coarse_runs_.product(
        [{"experiment": GaussianBumpExperiment(amplitude, width)}]
    )

    reference_runs = reference_runs_.product(
        [{"experiment": GaussianBumpExperiment(amplitude, width)}]
    )

    coarse_runs = [sc.Scheme(choice) for choice in coarse_runs]
    reference_runs = [sc.Scheme(choice) for choice in reference_runs]

    return coarse_runs, reference_runs


all_runs = [make_runs(amp) for amp in amplitudes]


def post_process(coarse_runs, reference_run):
    # results, columns = load_results(coarse_runs, reference_run)
    # labels = TableLabels()

    # filename = coarse_runs[0]["experiment"].short_id()
    # write_convergence_table(results, columns, labels, filename)
    # write_convergence_plots(results, columns, labels, filename)
    # plot_visual_convergence(results, columns, labels, filename)

    for coarse_run in coarse_runs:
        coarse_dir = folder_name(coarse_run)
        coarse_grid = load_grid(coarse_dir)
        data_files = find_data_files(coarse_dir)

        key = "rho"
        for l, data_file in enumerate(data_files):
            u_coarse = load_data(data_file, find_steady_state_file(coarse_dir))

            q = u_coarse.cvars[key]
            dq = u_coarse.dvars[key]

            scatter = ScatterPlot()
            scatter(coarse_grid, q)
            scatter.save(f"{coarse_dir}/scatter_{key}_{l:04d}.png")

            scatter = ScatterPlot()
            scatter(coarse_grid, dq)
            scatter.save(f"{coarse_dir}/triplot_d{key}_{l:04d}.png")


class TableLabels:
    def __call__(self, col):
        return " ".join(str(v) for v in col.values())


def generate_grids():
    generate_spherical_grids(grid_name_geo, radius, lc_rel, mesh_levels)


def main():
    parser = default_cli_parser("'gaussian_bump' numerical experiment.")
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
            if not args.reference_only:
                launch_all(c, force=args.force, queue_args=queue_args)

            if args.reference or args.reference_only:
                launch_all(r, force=args.force, queue_args=queue_args)

    if args.post_process:
        for c, r in all_runs:
            post_process(c, r[0])

    if args.copy_to_paper:
        d = "${HOME}/git/papers/LucGrosheintz/papers/unstructured_well_balancing/img/gaussian_bump"
        d = os.path.expandvars(d)

        for c, _ in all_runs:
            stem = c[0]["experiment"].short_id()
            files = glob.glob(stem + "*.tex")

            for f in files:
                shutil.copy(f, os.path.join(d, os.path.basename(f)))


if __name__ == "__main__":
    main()
