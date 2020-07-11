#! /usr/bin/env python3

import os
import shutil
import glob
import h5py
import numpy as np
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
from tiwaz.gmsh import generate_spherical_shell_grids
from tiwaz.gmsh import generate_spherical_grids
from tiwaz.gmsh import decompose_grids
from tiwaz.gmsh import renumber_grids
from tiwaz.gmsh import GridNamingScheme
from tiwaz.site_details import MPIHeuristics
from tiwaz.work_estimate import ZisaWorkEstimate
from tiwaz.queue_args import MPIQueueArgs
from tiwaz.launch_params import folder_name
from tiwaz.post_process import load_data, load_grid
from tiwaz.post_process import find_data_files, find_last_data_file
from tiwaz.post_process import find_steady_state_file
from tiwaz.tri_plot import TriPlot

cm = 1
m = 1e2 * cm
km = 1e3 * m

one_dimensional_profile = "data/stellar_convection/non_ideal_gas_profile.h5"


def mass_mixing_ratio(filename):
    with h5py.File(filename, "r") as h5:
        r = np.array(h5["x1"])
        rho = np.array(h5["conserved/density"])
        C12 = np.array(h5["advected/C12"])
        H1 = np.array(h5["advected/H1"])
        He4 = np.array(h5["advected/He4"])
        Ne20 = np.array(h5["advected/Ne20"])
        O16 = np.array(h5["advected/O16"])
        Si28 = np.array(h5["advected/Si28"])

    km = 1e5
    I = np.logical_and(5000 * km < r, r < 40000 * km)
    # print(C12[I][0])
    # print(H1[I][0])
    # print(He4[I][0])
    # print(Ne20[I][0])
    # print(O16[I][0])
    # print(Si28[I][0])

    print(np.min(C12 + H1 + He4 + Ne20 + O16 + Si28))
    print(np.max(C12 + H1 + He4 + Ne20 + O16 + Si28))
    # plt.plot(r[I] / km, rho[I])
    # plt.plot(C12)
    # plt.plot(H1)
    # plt.plot(He4)
    # plt.plot(Ne20)
    # plt.plot(O16)
    # plt.plot(Si28)


class StellarConvectionExperiment(sc.Subsection):
    def __init__(self):
        mass_mixing_ratio(one_dimensional_profile)

        super().__init__(
            {
                "name": "stellar_convection",
                "initial_conditions": {"profile": one_dimensional_profile,},
            }
        )

    def short_id(self):
        return self["name"]


eos = sc.HelmholtzEOS(
    "data/stellar_convection/helm_table.dat",
    element_keys=["C12", "H1", "He4", "Ne20", "O16", "Si28"],
    mass_number=[12, 1, 4, 20, 16, 28],
    charge_number=[6, 1, 2, 10, 8, 14],
)
gravity = sc.RadialGravity(one_dimensional_profile)
euler = sc.Euler(eos, gravity)

t_end = 0.09
time = sc.Time(t_end=t_end)
io = sc.IO(
    "hdf5",
    "stellar_convection",
    n_snapshots=1,
    parallel_strategy="gathered",
    n_writers=4,
)


def make_work_estimate():
    n0 = 10_000

    t0 = 2 * timedelta(seconds=t_end / 1e-1 * 60 * 96)
    b0 = 0.0
    o0 = 100 * 1e6
    unit_work = 256

    return ZisaWorkEstimate(n0=n0, t0=t0, b0=b0, o0=o0, unit_work=512)


grid_name = GridNamingScheme("stellar_convection")

parallelization = {"mode": "mpi"}

radii = [5000 * km, 40_000 * km]
# mesh_levels = list(range(0, 6)) + [7]
mesh_levels = list(range(1, 2))
lc_rel = {l: 0.1 * 0.5 ** l for l in mesh_levels}

coarse_grid_levels = [1]

coarse_grid_choices = {
    "grid": [
        sc.Grid(grid_name.config_string(l, parallelization), l)
        for l in coarse_grid_levels
    ]
}
coarse_grids = all_combinations(coarse_grid_choices)
reference_grid = sc.Grid(grid_name.config_string(7, parallelization), 7)

independent_choices = {
    "euler": [euler],
    "io": [io],
    "time": [time],
    "parallelization": [parallelization],
    "boundary-condition": [sc.BoundaryCondition("frozen")],
    "debug": [{"global_indices": False, "stencils": False}],
}

dependent_choices_a = {
    # "flux-bc": [sc.FluxBC("constant")],
    # "well-balancing": [sc.WellBalancing("constant")],
    "flux-bc": [sc.FluxBC(k) for k in ["constant", "isentropic"]],
    "well-balancing": [sc.WellBalancing(k) for k in ["constant", "isentropic"]],
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

coarse_runs_ = model_choices.product(coarse_grids)


def make_runs():
    coarse_runs = coarse_runs_.product([{"experiment": StellarConvectionExperiment()}])
    coarse_runs = [sc.Scheme(choice) for choice in coarse_runs]

    return coarse_runs


all_runs = make_runs()


def compute_parts(mesh_levels, host):
    work_estimate = make_work_estimate()

    heuristics = MPIHeuristics(host=host)
    queue_args = MPIQueueArgs(
        work_estimate, t_min=None, t_max=None, heuristics=heuristics
    )

    parts_ = dict()
    for l in mesh_levels:
        parts_[l] = [2] + [
            queue_args.n_mpi_tasks({"grid": {"file": grid_name.msh_h5(l)}})
        ]

    return parts_


def generate_grids(cluster, must_generate, must_decompose):
    geo_name = lambda l: grid_name.geo(l)
    msh_h5_name = lambda l: grid_name.msh_h5(l)

    if must_generate:
        for l in mesh_levels:
            shutil.rmtree(grid_name.dir(l), ignore_errors=True)

        generate_spherical_grids(
            geo_name, radii[1], lc_rel, mesh_levels, with_halo=True
        )
        renumber_grids(msh_h5_name, mesh_levels)

    if must_decompose:
        decompose_grids(msh_h5_name, mesh_levels, compute_parts(mesh_levels, cluster))


def main():
    parser = default_cli_parser("'stellar_convection' numerical experiment.")
    args = parser.parse_args()

    generate_grids(args.cluster, args.generate_grids, args.decompose_grids)

    if args.run:
        build_zisa()

        t_min = timedelta(minutes=30)
        t_max = timedelta(hours=24)
        work_estimate = make_work_estimate()

        if parallelization["mode"] == "mpi":
            queue_args = MPIQueueArgs(work_estimate, t_min=t_min, t_max=t_max)
        else:
            queue_args = dict()

        launch_all(all_runs, force=args.force, queue_args=queue_args)

    if args.post_process:
        for c in all_runs:
            post_process(c)

    if args.copy_to_paper:
        d = "${HOME}/git/papers/LucGrosheintz/papers/unstructured_well_balancing/img/stellar_convection"
        d = os.path.expandvars(d)

        for c in all_runs:
            stem = c[0]["experiment"].short_id()
            files = glob.glob(stem + "*.tex")

            for f in files:
                shutil.copy(f, os.path.join(d, os.path.basename(f)))


if __name__ == "__main__":
    main()
