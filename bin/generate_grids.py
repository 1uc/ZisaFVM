#! /usr/bin/env python3
# encoding: utf-8

import argparse
import glob

import tiwaz
import tiwaz.gmsh as gmsh


def minimal_geo_files():
    path = "grids"
    return [path + "/dbg.geo", path + "/small.geo"] + [
        path + f"/convergence/unit_square_{d}.geo" for d in range(4)
    ]


def matching_geo_files(pattern):
    return glob.glob("{}*.geo".format(pattern))


def all_geo_files():
    path = "grids"
    return glob.glob("{}/*.geo".format(path)) + glob.glob(f"{path}/**/*.geo")


def generate_geo_files():
    generate_unit_square_geo_files()
    generate_unit_cube_geo_files()


def generate_convergence_geo_files(basename, n_grids):
    template = tiwaz.utils.read_txt(f"grids/convergence/{basename}.tmpl")

    for k in range(n_grids):
        geo = template.replace("LC", str(0.1 * 2 ** (-k)))
        tiwaz.utils.write_txt(f"grids/convergence/{basename}_{k}.geo", geo)


def generate_unit_square_geo_files():
    generate_convergence_geo_files("unit_square", 5)
    generate_convergence_geo_files("unit_square_with_halo", 5)


def generate_unit_cube_geo_files():
    generate_convergence_geo_files("unit_cube", 4)
    generate_convergence_geo_files("unit_cube_with_halo", 4)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate .msh files.")

    parser.add_argument(
        "--minimal",
        action="store_true",
        help="Generate only those grids required for testing.",
    )

    parser.add_argument(
        "pattern", nargs="?", default=None, help="Generate only for the matching grids."
    )

    args = parser.parse_args()

    generate_geo_files()
    if args.minimal:
        geos = minimal_geo_files()
    elif args.pattern:
        geos = matching_geo_files(args.pattern)
    else:
        geos = all_geo_files()

    gmsh.generate_grids(geos)
