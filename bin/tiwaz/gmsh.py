import subprocess
import multiprocessing
import meshio
import h5py
import os

from .utils import read_txt, write_txt
from .site_details import zisa_home_directory
from .launch_params import build_target


def generate_grids_from_template(template_name, filename, substitutions, levels):
    template = read_txt(template_name)

    for l, s in zip(levels, substitutions):
        geo = template
        for key, value in s.items():
            geo = geo.replace(key, value)

        write_txt(filename(l), geo)

    generate_grids([filename(l) for l in levels])


def generate_spherical_grids(filename, radius, lc_rel, levels):
    template_name = "grids/sphere.tmpl"
    substitutions = [{"RADIUS": str(radius), "LC_REL": str(lc_rel[l])} for l in levels]

    generate_grids_from_template(template_name, filename, substitutions, levels)


def generate_circular_grids(filename, radius, lc_rel, levels):
    template_name = "grids/circle.tmpl"
    substitutions = [{"RADIUS": str(radius), "LC_REL": str(lc_rel[l])} for l in levels]

    generate_grids_from_template(template_name, filename, substitutions, levels)


def generate_cube_grids(filename, width, lc_rel, levels):
    template_name = "grids/cube.tmpl"
    substitutions = [{"WIDTH": str(radius), "LC_REL": str(lc_rel[l])} for l in levels]

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


def decompose_grids(grid_name_generator, mesh_levels, parts):
    zisa_home = zisa_home_directory()
    build_target("domain-decomposition")

    dd_binary = zisa_home + "/build-release/domain-decomposition"
    for l in mesh_levels:
        grid_name = grid_name_generator(l)
        if grid_name.endswith(".msh.h5"):
            outdir = grid_name[:-7]
        elif grid_name.endswith(".h5"):
            outdir = grid_name[:-3]
        else:
            raise Exception(f"Implement the missing logic. [{grid_name}]")

        for n in parts[l]:
            output = outdir + f"/{n}"
            os.makedirs(output, exist_ok=True)
            subprocess.run([dd_binary, "-n", str(n), "--grid", grid_name, "-o", output])


def run_gmesh(geo):
    gmsh = "bin/gmsh"
    cmd = [gmsh, "-3", "-format", "msh22", geo]
    subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    convert_msh_to_hdf5(geo[:-4] + ".msh")

    return geo


def generate_grids(geos):
    with multiprocessing.Pool() as pool:
        for i, geo in enumerate(pool.imap_unordered(run_gmesh, geos)):
            print("[{:2d}/{}] {}".format(i + 1, len(geos), geo))
