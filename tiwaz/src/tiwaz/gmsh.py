import subprocess
import multiprocessing
import meshio
import h5py
import os
import datetime
import argparse
import scibs

from tiwaz.utils import read_txt, write_txt
from tiwaz.site_details import zisa_home_directory
from tiwaz.site_details import make_batch_system
from tiwaz.site_details import zisa_build_directory
from tiwaz.site_details import max_cores_per_node
from tiwaz.launch_params import build_target


class GridNamingScheme:
    def __init__(self, root_dir, basename, version):
        self.root_dir = root_dir
        self.basename = basename

        if isinstance(version, int):
            self.version = f"{version:04d}"
        else:
            self.version = version

    def dir(self, l):
        return os.path.join(
            self.root_dir, f"{self.basename}", f"l{l:02d}", self.version
        )

    def geo(self, l):
        return self._stem(l) + ".geo"

    def msh_h5(self, l):
        return self._stem(l) + ".msh.h5"

    def config_string(self, l, parallelization):
        if parallelization["mode"] == "mpi":
            return self.dir(l)
        else:
            return self.msh_h5(l)

    def _stem(self, l):
        return os.path.join(self.dir(l), "grid")


def generate_geo_from_template(template_name, filename, substitutions, levels):
    template = read_txt(template_name)

    for l, s in zip(levels, substitutions):
        geo = template
        for key, value in s.items():
            geo = geo.replace(key, value)

        write_txt(filename(l), geo)

    return [filename(l) for l in levels]


def generate_spherical_shell_grids(filename, radii, lc_rel, levels, with_halo):
    template_name = "grids.light/spherical_shell" + (
        "_with_halo.tmpl" if with_halo else ".tmpl"
    )

    with open("grids.light/sphere.macro", "r") as f:
        sphere_macro = f.read()

    substitutions = [
        {
            "INNER_RADIUS": str(radii[0]),
            "OUTER_RADIUS": str(radii[1]),
            "LC_REL": str(lc_rel[l]),
            "SPHERE_MACRO": sphere_macro,
        }
        for l in levels
    ]

    return generate_geo_from_template(template_name, filename, substitutions, levels)


def generate_geo_spherical_grids(filename, radius, lc_rel, levels, with_halo):
    template_name = (
        "grids.light/sphere_with_halo.tmpl" if with_halo else "grids.light/sphere.tmpl"
    )
    with open("grids.light/sphere.macro", "r") as f:
        sphere_macro = f.read()

    substitutions = [
        {"RADIUS": str(radius), "LC_REL": str(lc_rel[l]), "SPHERE_MACRO": sphere_macro}
        for l in levels
    ]

    return generate_geo_from_template(template_name, filename, substitutions, levels)


def generate_geo_circular_grids(filename, radius, lc_rel, levels, with_halo):
    template_name = (
        "grids.light/circle_with_halo.tmpl" if with_halo else "grids.light/circle.tmpl"
    )
    substitutions = [{"RADIUS": str(radius), "LC_REL": str(lc_rel[l])} for l in levels]

    return generate_geo_from_template(template_name, filename, substitutions, levels)


def generate_geo_cube_grids(filename, domain, lc_rel, levels, with_halo):
    template_name = (
        "grids.light/cube_with_halo.tmpl" if with_halo else "grids.light/cube.tmpl"
    )

    with open("grids.light/cube.macro", "r") as f:
        cube_macro = f.read()

    substitutions = [
        {
            "X0": str(domain[0]),
            "X1": str(domain[1]),
            "CUBE_MACRO": cube_macro,
            "LC": str(lc_rel[l] * (domain[1] - domain[0])),
        }
        for l in levels
    ]

    return generate_geo_from_template(template_name, filename, substitutions, levels)


def generate_geo_square_grids(filename, domain, lc_rel, levels, with_halo):
    template_name = (
        "grids.light/square_with_halo.tmpl" if with_halo else "grids.light/square.tmpl"
    )

    with open("grids.light/square.macro", "r") as f:
        square_macro = f.read()

    substitutions = [
        {
            "X0": str(domain[0]),
            "X1": str(domain[1]),
            "SQUARE_MACRO": square_macro,
            "LC": str(lc_rel[l] * (domain[1] - domain[0])),
        }
        for l in levels
    ]

    return generate_geo_from_template(template_name, filename, substitutions, levels)


def convert_msh_to_hdf5(msh):
    grid = meshio.read(msh, file_format="gmsh")

    h5_file = msh[:-4] + ".msh.h5"
    with h5py.File(h5_file, "w") as h5:
        h5["vertices"] = grid.points

        cells = grid.cells_dict
        if "tetra" in cells:
            h5["n_dims"] = 3
            h5["vertex_indices"] = cells["tetra"]
        else:
            h5["n_dims"] = 2
            h5["vertex_indices"] = cells["triangle"]


def renumber_grid_command(grid_name):
    assert grid_name.endswith(".msh.h5")

    zisa_home = zisa_home_directory()
    binary = os.path.join(zisa_home, "build-release/renumber-grid")
    return [binary, "--grid", grid_name]


def renumber_grid(grid_name):
    subprocess.run(renumber_grid_command(grid_name), check=True)


def renumber_grids(grid_name_generator, mesh_levels):
    zisa_home = zisa_home_directory()
    build_target("renumber-grid")

    for l in mesh_levels:
        renumber_grid(grid_name_generator(l))


def decompose_grid_command(msh_h5, n_proc, n_parts):
    outdir = os.path.join(os.path.dirname(msh_h5), "partitioned")

    build_dir = zisa_build_directory()
    bin = os.path.join(build_dir, "domain-decomposition")

    # fmt: off
    return [
        bin,
        "--partitions", str(n_parts),
        "--workers", str(n_proc),
        "--grid", msh_h5,
        "-o", outdir,
    ]
    # fmt: on


def submit_grid_decomposition(batch_system, msh_h5, n_proc, n_parts, resources):
    cmd = decompose_grid_command(msh_h5, n_proc=n_proc, n_parts=n_parts)
    batch_system.submit(scibs.Job(cmd=cmd, resources=resources))


def generate_grid(geo):
    # fmt: off
    cmd = [
        "gmsh",
        "-optimize_netgen",
        "-optimize_threshold", "0.9",
        "-3",
        "-format", "msh22",
        geo,
    ]
    # fmt: on

    subprocess.run(cmd, check=True)
    convert_msh_to_hdf5(geo[:-4] + ".msh")
    renumber_grid(geo[:-4] + ".msh.h5")

    return geo


def generate_grid_command(geo):
    return ["python3", __file__, "--geo", geo]


def submit_grid_generation(batch_system, geo, resources):
    cmd = generate_grid_command(geo)
    batch_system.submit(scibs.Job(cmd=cmd, resources=resources))


def generate_smallish_grids(geo_files, hours=24, mem=16e9):
    build_target("renumber-grid")

    batch_system = make_batch_system()

    wall_clock = datetime.timedelta(hours=hours)
    resources = scibs.OMPResource(
        n_omp_threads=1, wall_clock=wall_clock, total_memory=mem
    )

    for geo in geo_files:
        submit_grid_generation(batch_system, geo, resources)


def decompose_smallish_grids(msh_h5_files, parts, hours=48, mem=16e9):
    build_target("domain-decomposition")

    batch_system = make_batch_system()
    max_cores = max_cores_per_node()

    wall_clock = datetime.timedelta(hours=hours)
    resources = scibs.OMPResource(
        n_threads=max_cores, wall_clock=wall_clock, total_memory=mem
    )

    for msh_h5, n_parts in zip(msh_h5_files, parts):
        n_proc = min(n_parts, max_cores)
        submit_grid_decomposition(
            batch_system, msh_h5, n_proc=n_proc, n_parts=n_parts, resources=resources
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run GMSH to generate grids")
    parser.add_argument("--geo", type=str, help="The `*.geo` on which to run `gmsh`.")

    args = parser.parse_args()
    generate_grid(args.geo)
