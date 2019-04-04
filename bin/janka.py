#!/usr/bin/env python3

import os
import shutil
import glob
import subprocess

import numpy as np
import matplotlib.pyplot as plt

import tiwaz
import tiwaz.scheme as sc
import tiwaz.gmsh as gmsh

from tiwaz.launch_params import all_combinations, pointwise_combinations
from tiwaz.launch_params import build_zisa

from tiwaz.cli_parser import default_cli_parser
from tiwaz.launch_job import launch_all
from tiwaz.scatter_plot import ScatterPlot, plot_visual_convergence
from tiwaz.launch_params import folder_name
from tiwaz.post_process import load_data, load_grid
from tiwaz.post_process import find_data_files, find_last_data_file
from tiwaz.post_process import find_steady_state_file
from tiwaz.utils import read_txt, write_txt
from tiwaz.tri_plot import tri_plot

janka_params = dict(rho_bounce = 2e14,
                    gamma1 = 1.33,
                    gamma2 = 2.5,
                    gamma_thermal = 1.5,
                    E1 = 1.46925e15)

polytropic_index_n = 3.0
class JankaExperiment(sc.Subsection):
    def __init__(self):
        super().__init__({
            "name": "janka",
            "initial_conditions": {
                "gamma": (polytropic_index_n + 1) / polytropic_index_n
            }
        })

    def short_id(self):
        return self["name"]

# This is just in case we want to run the same experiment with many parameters.
experiment_params = [None]

eos = sc.JankaEOS(**janka_params)
gravity = sc.GeneralPolytropeGravity(
    rho_center=1e10,
    polytropic_index_n=polytropic_index_n,
    radius=1.2e8,
    G = 6.674e-8
)
euler = sc.Euler(eos, gravity)

time = sc.Time(t_end=0.1)
# io = sc.IO("hdf5", "janka", n_snapshots=1)
io = sc.IO("opengl", "janka", steps_per_frame=1)

def grid_name(l):
    return "grids/janka_{}.msh".format(l)


mesh_levels = list(range(0, 5))
mesh_sizes = [1e-1 * 0.5**l for l in mesh_levels]

coarse_grid_levels = list(range(0, 1))
coarse_grid_names = [grid_name(level) for level in coarse_grid_levels]

coarse_grid_choices = {
    "grid": [
        sc.Grid(grid_name(l), l) for l in coarse_grid_levels
    ]
}
coarse_grids = all_combinations(coarse_grid_choices)
reference_grid = sc.Grid(grid_name(mesh_levels[-1]), mesh_levels[-1])

independent_choices = {
    "euler": [euler],

    "flux-bc": [sc.FluxBC("constant")],

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
        sc.ODE("ForwardEuler")
        # sc.ODE("SSP3")
        # sc.ODE("SSP3")
    ],

    "quadrature": [
        # sc.Quadrature(1),
        # sc.Quadrature(1)
        sc.Quadrature(1)
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
        "experiment": JankaExperiment()
    }])

    reference_runs = reference_runs_.product([{
        "experiment": JankaExperiment()
    }])

    coarse_runs = [sc.Scheme(choice) for choice in coarse_runs]
    reference_runs = [sc.Scheme(choice) for choice in reference_runs]

    return coarse_runs, reference_runs

all_runs = [make_runs(param) for param in experiment_params]

def post_process(coarse_runs, reference_run):
    coarse_run = coarse_runs[0]
    coarse_dir = folder_name(coarse_run)
    coarse_grid = load_grid(coarse_dir)


    data_file = find_data_files(coarse_dir)[1]
    u_coarse = load_data(data_file,
                         find_steady_state_file(coarse_dir))


    rho = u_coarse.cvars["rho"]
    vx = u_coarse.cvars["mv1"] / rho
    vy = u_coarse.cvars["mv2"] / rho

    plt.figure()
    tri_plot(coarse_grid, np.log10(rho))

    plt.figure()
    tri_plot(coarse_grid, vx / rho)

    plt.figure()
    tri_plot(coarse_grid, vy / rho)

    plot = ScatterPlot()
    plot(coarse_grid, np.log10(u_coarse.cvars["rho"]))
    plt.show()

class TableLabels:
    def __call__(self, col):
        return " ".join(str(v) for v in col.values())

def generate_grids():
    gmsh_template = read_txt("grids/janka.tmpl")
    filename = "grids/janka-{}.geo"

    for l in mesh_levels:
        lc_rel = mesh_sizes[l]
        geo = gmsh_template.replace("LC_REL", str(lc_rel))
        write_txt(filename.format(l), geo)

    gmsh.generate_grids([filename.format(l) for l in mesh_levels])

def main():
    parser = default_cli_parser("'janka' numerical experiment.")

    parser.add_argument(
        "--config-only",
        action='store_true',
        help="Write the config file to the working directory."
    )

    args = parser.parse_args()

    if args.generate_grids:
        generate_grids()

    if args.run:
        build_zisa()

        for c, r in all_runs:
            launch_all(c, force=args.force)

            if args.reference:
                launch_all(r, force=args.force)

    if args.post_process:
        for c, r in all_runs:
            post_process(c, r)

    if args.config_only:
        for c, r in all_runs:
            for cc in c:
                cc.save("config--" + cc.folder_name() + ".json")

            for rr in r:
                rr.save("config--" + rr.folder_name() + ".json")

    if args.copy_to_paper:
        dir = "${HOME}/git/papers/LucGrosheintz/papers/unstructured_well_balancing/img/janka"
        dir = os.path.expandvars(dir)

        for c, _ in all_runs:
            stem = c[0]["experiment"].short_id()
            files = glob.glob(stem + "*.tex")

            for f in files:
                shutil.copy(f, os.path.join(dir, os.path.basename(f)))

if __name__ == "__main__":
    main()







