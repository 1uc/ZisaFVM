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
from tiwaz.work_estimate import ZisaWorkEstimate
from tiwaz.queue_args import MPIQueueArgs
from tiwaz.launch_params import folder_name
from tiwaz.post_process import load_data, load_grid
from tiwaz.post_process import find_data_files, find_last_data_file
from tiwaz.post_process import find_steady_state_file
from tiwaz.tri_plot import TriPlot


class GaussianBumpExperiment(sc.Subsection):
    def __init__(self, amplitude, width):
        super().__init__(
            {
                "name": "gaussian_bump",
                "initial_conditions": {"amplitude": amplitude, "width": width},
            }
        )

    def short_id(self):
        amp = self["initial_conditions"]["amplitude"]
        return self["name"] + "_amp{:.2e}".format(amp)


amplitudes = [0.0, 1e-6, 1e-2]
width = 0.05

eos = sc.IdealGasEOS(gamma=2.0, r_gas=1.0)
gravity = sc.PolytropeGravity()
euler = sc.Euler(eos, gravity)

t_end = 0.09
time = sc.Time(t_end=t_end)
io = sc.IO("hdf5", "gaussian_bump", n_snapshots=1)
# io = sc.IO("opengl", "gaussian_bump", steps_per_frame=1)


def make_work_estimate():
    n0 = sc.read_n_cells(grid_name_hdf5(4))

    # measured on Euler on L=4 with 96 cores.
    t0 = 2 * timedelta(seconds=t_end / 1e-1 * 30 * 96)

    b0 = 0.5 * 1e9
    o0 = 1.0 * 1e9

    return ZisaWorkEstimate(n0=n0, t0=t0, b0=b0, o0=o0)


def grid_name_stem(l):
    return "grids/gaussian_bump-{}".format(l)


def grid_name_geo(l):
    return grid_name_stem(l) + ".geo"


def grid_name_hdf5(l):
    return grid_name_stem(l) + ".msh.h5"


radius = 0.5
mesh_levels = list(range(0, 6))
lc_rel = {l: 0.1 * 0.5 ** l for l in mesh_levels}

coarse_grid_levels = list(range(0, 4))
coarse_grid_names = [grid_name_hdf5(level) for level in coarse_grid_levels]

coarse_grid_choices = {
    "grid": [sc.Grid(grid_name_hdf5(l), l) for l in coarse_grid_levels]
}
coarse_grids = all_combinations(coarse_grid_choices)
reference_grid = sc.Grid(grid_name_hdf5(4), 4)

independent_choices = {
    "euler": [euler],
    "flux-bc": [sc.FluxBC("isentropic")],
    "well-balancing": [sc.WellBalancing("constant"), sc.WellBalancing("isentropic")],
    "io": [io],
    "time": [time],
    "parallelization": [{"mode": "mpi"}],
}

dependent_choices = {
    "reconstruction": [
        sc.Reconstruction("CWENO-AO", [1]),
        sc.Reconstruction(
            "CWENO-AO", [2, 2, 2, 2], overfit_factors=[3.0, 2.0, 2.0, 2.0]
        ),
        sc.Reconstruction("CWENO-AO", [3, 2, 2, 2]),
        sc.Reconstruction("CWENO-AO", [4, 2, 2, 2]),
    ],
    "ode": [sc.ODE("ForwardEuler"), sc.ODE("SSP2"), sc.ODE("SSP3"), sc.ODE("SSP3")],
    "quadrature": [
        sc.Quadrature(1),
        sc.Quadrature(2),
        sc.Quadrature(3),
        sc.Quadrature(4),
    ],
}

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
    results, columns = load_results(coarse_runs, reference_run)
    labels = TableLabels()

    filename = coarse_runs[0]["experiment"].short_id()
    write_convergence_table(results, columns, labels, filename)
    write_convergence_plots(results, columns, labels, filename)
    plot_visual_convergence(results, columns, labels, filename)

    for coarse_run in coarse_runs:
        coarse_dir = folder_name(coarse_run)
        coarse_grid = load_grid(coarse_dir)

        data_files = find_data_files(coarse_dir)

        key = "rho"
        for l, data_file in enumerate(data_files):
            u_coarse = load_data(data_file, find_steady_state_file(coarse_dir))

            rho = u_coarse.cvars[key]
            drho = u_coarse.dvars[key]

            # vx = u_coarse.cvars["mv1"] / rho
            # vy = u_coarse.cvars["mv2"] / rho

            trip = TriPlot()
            trip.color_plot(coarse_grid, rho)
            # trip.quiver(coarse_grid, vx, vy)

            raise Exception("next line is broken. fix first.")
            plt.savefig(f"gaussian_bump_scatter_{key}.png", dpi=300)

            trip = TriPlot()
            trip.color_plot(coarse_grid, drho)
            # trip.quiver(coarse_grid, vx, vy)

            plt.savefig(f"{coarse_dir}_d{key}_{l:04d}.png", dpi=300)


class TableLabels:
    def __call__(self, col):
        return " ".join(str(v) for v in col.values())


def generate_grids():
    generate_circular_grids(grid_name_geo, radius, lc_rel, mesh_levels)


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
