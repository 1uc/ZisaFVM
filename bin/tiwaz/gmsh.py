import subprocess
import multiprocessing

from . utils import read_txt, write_txt

def generate_circular_grids(filename, radius, lc_rel, levels):
    template = read_txt("grids/circle.tmpl")

    for l in levels:
        geo = template.replace("RADIUS", str(radius))
        geo = geo.replace("LC_REL", str(lc_rel[l]))

        write_txt(filename(l), geo)

    generate_grids(([filename(l) for l in levels]))


def run_gmesh(geo):
    gmsh = "bin/gmsh"
    cmd = [gmsh, "-2", geo]
    output = subprocess.check_call(cmd,
                                   stdout=subprocess.DEVNULL,
                                   stderr=subprocess.DEVNULL)

    return geo


def generate_grids(geos):
    with multiprocessing.Pool() as pool:
        for i, geo in enumerate(pool.imap_unordered(run_gmesh, geos)):
            print("[{:2d}/{}] {}".format(i+1, len(geos), geo))
