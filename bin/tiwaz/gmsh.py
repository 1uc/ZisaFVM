import subprocess
import multiprocessing
import meshio
import h5py
import os

from .utils import read_txt, write_txt
from .site_details import zisa_home_directory
from .launch_params import build_target


class GridNamingScheme:
    def __init__(self, basename):
        self.basename = basename

    def dir(self, l):
        return f"grids/{self.basename}-{l}"

    def geo(self, l):
        return self._stem(l) + ".geo"

    def msh_h5(self, l):
        return self._stem(l) + ".msh.h5"

    def config_string(self, l, parallelization):
        return self.dir(l)

    def _stem(self, l):
        return f"{self.dir(l)}/grid"


def generate_grids_from_template(template_name, filename, substitutions, levels):
    template = read_txt(template_name)

    for l, s in zip(levels, substitutions):
        geo = template
        for key, value in s.items():
            geo = geo.replace(key, value)

        print(filename(l))
        write_txt(filename(l), geo)

    generate_grids([filename(l) for l in levels])


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

    generate_grids_from_template(template_name, filename, substitutions, levels)


def generate_spherical_grids(filename, radius, lc_rel, levels, with_halo):
    template_name = (
        "grids.light/sphere_with_halo.tmpl" if with_halo else "grids.light/sphere.tmpl"
    )
    with open("grids.light/sphere.macro", "r") as f:
        sphere_macro = f.read()

    substitutions = [
        {"RADIUS": str(radius), "LC_REL": str(lc_rel[l]), "SPHERE_MACRO": sphere_macro}
        for l in levels
    ]

    generate_grids_from_template(template_name, filename, substitutions, levels)


def generate_circular_grids(filename, radius, lc_rel, levels, with_halo):
    template_name = (
        "grids.light/circle_with_halo.tmpl" if with_halo else "grids.light/circle.tmpl"
    )
    substitutions = [{"RADIUS": str(radius), "LC_REL": str(lc_rel[l])} for l in levels]

    generate_grids_from_template(template_name, filename, substitutions, levels)


def generate_cube_grids(filename, domain, lc_rel, levels, with_halo):
    template_name = (
        "grids.light/cube_with_halo.tmpl" if with_halo else "grids.light/cube.tmpl"
    )

    substitutions = [
        {
            "X0": str(domain[0]),
            "X1": str(domain[1]),
            "CUBE_MACRO": cube_macro,
            "LC": str(lc_rel[l] * (domain[1] - domain[0])),
        }
        for l in levels
    ]

    generate_grids_from_template(template_name, filename, substitutions, levels)


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


def renumber_grids(grid_name_generator, mesh_levels):
    zisa_home = zisa_home_directory()
    build_target("renumber-grid")

    binary = zisa_home + "/build-release/renumber-grid"
    for l in mesh_levels:
        grid_name = grid_name_generator(l)
        assert grid_name.endswith(".msh.h5")

        subprocess.run([binary, "--grid", grid_name])


def decompose_grids(grid_name_generator, mesh_levels, parts):
    zisa_home = zisa_home_directory()
    build_target("domain-decomposition")

    dd_binary = zisa_home + "/build-release/domain-decomposition"
    for l in mesh_levels:
        grid_name = grid_name_generator(l)
        outdir = os.path.dirname(grid_name) + "/partitioned"

        for n in parts[l]:
            output = outdir + f"/{n}"
            os.makedirs(output, exist_ok=True)
            subprocess.run([dd_binary, "-n", str(n), "--grid", grid_name, "-o", output])


def run_gmesh(geo):
    gmsh = "bin/gmsh"
    cmd = [
        gmsh,
        "-optimize_netgen",
        "-optimize_threshold",
        "0.9",
        "-3",
        "-format",
        "msh22",
        geo,
    ]
    subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    convert_msh_to_hdf5(geo[:-4] + ".msh")

    return geo


def generate_grids(geos):
    with multiprocessing.Pool() as pool:
        for i, geo in enumerate(pool.imap_unordered(run_gmesh, geos)):
            print("[{:2d}/{}] {}".format(i + 1, len(geos), geo))
