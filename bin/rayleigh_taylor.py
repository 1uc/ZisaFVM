#! /usr/bin/env python3

import os
import shutil
import glob

import matplotlib.pyplot as plt

import tiwaz
import tiwaz.scheme as sc

from tiwaz.launch_params import all_combinations, pointwise_combinations
from tiwaz.launch_params import build_zisa

from tiwaz.post_process import load_results

from tiwaz.cli_parser import default_cli_parser
from tiwaz.launch_job import launch_all
from tiwaz.latex_tables import write_convergence_table
from tiwaz.scatter_plot import ScatterPlot, plot_visual_convergence
from tiwaz.launch_params import folder_name
from tiwaz.post_process import load_data, load_grid
from tiwaz.post_process import find_last_data_file, find_steady_state_file
from tiwaz.tri_plot import tri_plot

class RayleighTaylorExperiment(sc.Subsection):
    def __init__(self):
        super().__init__({
            "name": "rayleigh_taylor",
            "initial_conditions": {
                "drho" : 0.02,
                "amplitude": 0.01,
                "width": 0.1,
                "n_bumps": 2
            }
        })

    def short_id(self):
        return self["name"]

# This is just in case we want to run the same experiment with many parameters.
experiment_params = [None]

eos = sc.IdealGasEOS(gamma=2.0, r_gas=1.0)
gravity = sc.PolytropeGravityWithJump(rhoC=1.0, K_inner=1.0, K_outer=1.0)
euler = sc.Euler(eos, gravity)

time = sc.Time(t_end=0.01)
io = sc.IO("hdf5", "rayleigh_taylor", n_snapshots=1)

def grid_name(level):
    return "grids/polytrope-{:}.msh".format(level)

coarse_grid_levels = list(range(2, 3))
coarse_grid_names = [grid_name(level) for level in coarse_grid_levels]

coarse_grid_choices = {
    "grid": [
        sc.Grid(grid_name(l), l) for l in coarse_grid_levels
    ]
}
coarse_grids = all_combinations(coarse_grid_choices)
reference_grid = sc.Grid("grids/polytrope-4.msh", 4)

independent_choices = {
    "euler": [euler],

    "flux-bc": [sc.FluxBC("isentropic")],

    "well-balancing": [ # sc.WellBalancing("constant"),
                       sc.WellBalancing("isentropic")],

    "io": [io],

    "time": [time],
}

dependent_choices = {
    "reconstruction": [
        # sc.Reconstruction("CWENO-AO", [1]),
        # sc.Reconstruction("CWENO-AO", [2, 2, 2, 2])
        sc.Reconstruction("CWENO-AO", [3, 2, 2, 2])
    ],

    "ode": [
        # sc.ODE("ForwardEuler"),
        # sc.ODE("SSP3")
        sc.ODE("SSP3")
    ],

    "quadrature": [
        # sc.Quadrature(1),
        # sc.Quadrature(1)
        sc.Quadrature(3)
    ]
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
    "reference": [sc.Reference("isentropic", coarse_grid_names)]
}

base_choices = all_combinations(independent_choices)
model_choices = base_choices.product(pointwise_combinations(dependent_choices))

coarse_runs_ = model_choices.product(coarse_grids)
reference_runs_ = all_combinations(reference_choices)

def make_runs(params):
    coarse_runs = coarse_runs_.product([{
        "experiment": RayleighTaylorExperiment()
    }])

    reference_runs = reference_runs_.product([{
        "experiment": RayleighTaylorExperiment()
    }])

    coarse_runs = [sc.Scheme(choice) for choice in coarse_runs]
    reference_runs = [sc.Scheme(choice) for choice in reference_runs]

    return coarse_runs, reference_runs

all_runs = [make_runs(param) for param in experiment_params]

def post_process(coarse_runs, reference_run):
    coarse_run = coarse_runs[0]
    coarse_dir = folder_name(coarse_run)
    coarse_grid = load_grid(coarse_dir)
    u_coarse = load_data(find_last_data_file(coarse_dir),
                         find_steady_state_file(coarse_dir))

    tri_plot(coarse_grid, u_coarse.cvars["E"])

    plot = ScatterPlot()
    plot(coarse_grid, u_coarse.cvars["rho"])
    plt.show()

class TableLabels:
    def __call__(self, col):
        return " ".join(str(v) for v in col.values())

def main():
    parser = default_cli_parser("'rayleigh_taylor' numerical experiment.")
    args = parser.parse_args()

    if args.run:
        build_zisa()

        for c, r in all_runs:
            launch_all(c, force=args.force)

            if args.reference:
                launch_all(r, force=args.force)

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







