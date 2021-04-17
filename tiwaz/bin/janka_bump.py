#! /usr/bin/env python3

import os
import shutil
import glob

import numpy as np
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
from tiwaz.gmsh import generate_circular_grids

from tiwaz.launch_params import folder_name
from tiwaz.post_process import load_data, load_grid
from tiwaz.post_process import find_data_files, find_last_data_file
from tiwaz.post_process import find_steady_state_file
from tiwaz.tri_plot import TriPlot
from tiwaz.scatter_plot import ScatterPlot


class JankaBumpExperiment(sc.Subsection):
    def __init__(self, amplitude, width):
        super().__init__(
            {
                "name": "janka_bump",
                "initial_conditions": {"amplitude": amplitude, "width": width},
            }
        )

    def short_id(self):
        amp = self["initial_conditions"]["amplitude"]
        return self["name"] + "_amp{:.2e}".format(amp)


rho_center = 1e10
K_center = 4.897e14
G = 6.674e-8
radius = (2.0 * K_center / (4.0 * np.pi * G)) ** 0.5 * np.pi
amplitudes = [0.0, 1e-8, 1e-4, 1e-1, 1e2]
# amplitudes = [0.0]
# amplitudes = [1e2]
width = 0.10 * radius


eos = sc.IdealGasEOS(gamma=2.0, r_gas=1.0)
gravity = sc.PolytropeGravity(rhoC=rho_center, K=K_center, G=G)
euler = sc.Euler(eos, gravity)

time = sc.Time(t_end=0.01)
io = sc.IO("hdf5", "janka_bump", n_snapshots=1)


def grid_name_stem(l):
    return "grids/janka_bump-{}".format(l)


def grid_name_geo(l):
    return grid_name_stem(l) + ".geo"


def grid_name_msh(l):
    return grid_name_stem(l) + ".msh"


mesh_levels = list(range(0, 6))
lc_rel = {l: 0.1 * 0.5 ** l for l in mesh_levels}

coarse_grid_levels = list(range(0, 4))
coarse_grid_names = [grid_name_msh(level) for level in coarse_grid_levels]

coarse_grid_choices = {
    "grid": [sc.Grid(grid_name_msh(l), l) for l in coarse_grid_levels]
}
coarse_grids = all_combinations(coarse_grid_choices)
reference_grid = sc.Grid(grid_name_msh(5), 5)

independent_choices = {
    "euler": [euler],
    "io": [io],
    "time": [time],
}

dependent_choices_a = {
    "flux-bc": [sc.FluxBC("constant"), sc.FluxBC("isentropic")],
    "well-balancing": [sc.WellBalancing("constant"), sc.WellBalancing("isentropic")],
}


# dependent_choices = {
#     "reconstruction": [
#         # sc.Reconstruction("CWENO-AO", [1]),
#         sc.Reconstruction("CWENO-AO", [2, 2, 2, 2])
#         #sc.Reconstruction("CWENO-AO", [3, 2, 2, 2])
#     ],
#
#     "ode": [
#         # sc.ODE("ForwardEuler"),
#         sc.ODE("SSP3")
#         # sc.ODE("SSP3")
#     ],
#
#     "quadrature": [
#         # sc.Quadrature(1),
#         sc.Quadrature(1)
#         # sc.Quadrature(3)
#     ]
# }

dependent_choices_b = {
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
}

base_choices = pointwise_combinations(independent_choices)
base_choices = base_choices.product(pointwise_combinations(dependent_choices_a))
model_choices = base_choices.product(pointwise_combinations(dependent_choices_b))

coarse_runs_ = model_choices.product(coarse_grids)
reference_runs_ = all_combinations(reference_choices)


def make_runs(amplitude):
    coarse_runs = coarse_runs_.product(
        [{"experiment": JankaBumpExperiment(amplitude, width)}]
    )

    reference_runs = reference_runs_.product(
        [{"experiment": JankaBumpExperiment(amplitude, width)}]
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
    plot_visual_convergence(results, columns, labels, filename)

    for coarse_run in coarse_runs:
        coarse_dir = folder_name(coarse_run)
        coarse_grid = load_grid(coarse_dir)
        data_files = find_data_files(coarse_dir)

        # for data_file in data_files[-1:]:
        #     print(data_file)
        #     u_coarse = load_data(data_file, find_steady_state_file(coarse_dir))
        #
        #     drho = u_coarse.dvars["rho"]
        #     rho = u_coarse.cvars["rho"]
        #     vx = u_coarse.cvars["mv1"] / rho
        #     vy = u_coarse.cvars["mv2"] / rho
        #     p = u_coarse.xvars["p"]
        #     dE = u_coarse.dvars["E"]
        #
        #     # trip = TriPlot()
        #     # trip.color_plot(coarse_grid, rho)
        #     # trip.quiver(coarse_grid, vx, vy)
        #
        #     plot = ScatterPlot()
        #     plot(coarse_grid, np.log10(rho))
        #     plt.xlabel("Radius [cm]")
        #     plt.ylabel(r"$log_{10}(\rho)$ [g/cm^3]")
        #     plt.savefig(coarse_dir + "/img/rho_scatter.png")
        #
        #     # plot = ScatterPlot()
        #     # plot(coarse_grid, dE)
        #     # plt.xlabel("Radius [cm]")
        #     # plt.ylabel("log_10(\\rho) [g/cm^3]")
        #     # plt.savefig(coarse_dir + "/img/rho_scatter.png")
        #
        #     plot = ScatterPlot()
        #     plot(coarse_grid, np.log10(p))
        #     plt.xlabel("Radius [cm]")
        #     plt.ylabel(r"$log_{10}(p)$ [barye]")
        #     plt.savefig(coarse_dir + "/img/p_scatter.png")


class TableLabels:
    def __call__(self, col):
        return " ".join(str(v) for v in col.values())


def generate_grids():
    generate_circular_grids(grid_name_geo, radius, lc_rel, mesh_levels)


def main():
    parser = default_cli_parser("'janka_bump' numerical experiment.")
    args = parser.parse_args()

    if args.generate_grids:
        generate_grids()

    if args.run:
        build_zisa()
        queue_args = None

        for c, r in all_runs:
            launch_all(c, force=args.force, queue_args=queue_args)

            if args.reference:
                launch_all(r, force=args.force, queue_args=queue_args)

    if args.post_process:
        for c, r in all_runs:
            post_process(c, r[0])

    if args.copy_to_paper:
        dir = "${HOME}/git/papers/LucGrosheintz/papers/unstructured_well_balancing/img/janka_bump"
        dir = os.path.expandvars(dir)

        for c, _ in all_runs:
            stem = c[0]["experiment"].short_id()
            files = glob.glob(stem + "*.tex")

            for f in files:
                shutil.copy(f, os.path.join(dir, os.path.basename(f)))


if __name__ == "__main__":
    main()
