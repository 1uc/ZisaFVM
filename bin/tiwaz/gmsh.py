import subprocess
import multiprocessing
import meshio
import h5py

from .utils import read_txt, write_txt


def generate_circular_grids(filename, radius, lc_rel, levels):
    template = read_txt("grids/circle.tmpl")

    for l in levels:
        geo = template.replace("RADIUS", str(radius))
        geo = geo.replace("LC_REL", str(lc_rel[l]))

        write_txt(filename(l), geo)

    generate_grids(([filename(l) for l in levels]))


def generate_cube_grids(filename, width, lc_rel, levels):
    template = read_txt("grids/cube.tmpl")

    for l in levels:
        geo = template.replace("WIDTH", str(width))
        geo = geo.replace("LC_REL", str(lc_rel[l]))

        write_txt(filename(l), geo)

    generate_grids(([filename(l) for l in levels]))


def convert_msh_to_hdf5(msh):
    grid = meshio.read(msh, file_format="gmsh")

    h5_file = msh[:-4] + ".msh.h5"
    with h5py.File(h5_file, "w") as h5:
        h5["vertices"] = grid.points

        if "tetra" in grid.cells:
            h5["n_dims"] = 3
            h5["vertex_indices"] = grid.cells["tetra"]
        else:
            h5["n_dims"] = 2
            h5["vertex_indices"] = grid.cells["triangle"]


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
