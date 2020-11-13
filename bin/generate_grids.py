#! /usr/bin/env python3
# encoding: utf-8

import argparse
import glob

import tiwaz
import tiwaz.gmsh as gmsh


def minimal_geo_files():
    path = "grids"
    return [path + "/dbg.geo", path + "/small.geo"]


def matching_geo_files(pattern):
    return glob.glob(f"{path}*/*.geo")


def all_geo_files():
    path = "grids"
    return glob.glob(f"{path}/*.geo") + glob.glob(f"{path}/**/*.geo")


def generate_geo_files():
    generate_unit_square_geo_files()
    generate_unit_cube_geo_files()


def generate_convergence_geo_files(basename, n_grids):
    template = tiwaz.utils.read_txt(f"grids.light/convergence/{basename}.tmpl")

    for k in range(n_grids):
        geo = template.replace("LC", str(0.1 * 2 ** (-k)))
        tiwaz.utils.write_txt(f"grids/convergence/{basename}_{k}/grid.geo", geo)


def generate_unit_square_geo_files():
    n_grids = 4
    domain = [0.0, 1.0]
    levels = range(n_grids)
    lc_rel = [0.1 * 2 ** (-k) for k in levels]

    def no_halo_name(l):
        return f"grids/convergence/unit_square_{l}/grid.geo"

    def with_halo_name(l):
        return f"grids/convergence/unit_square_with_halo_{l}/grid.geo"

    gmsh.generate_square_grids(no_halo_name, domain, lc_rel, levels, with_halo=False)
    gmsh.generate_square_grids(with_halo_name, domain, lc_rel, levels, with_halo=True)


def generate_unit_cube_geo_files():
    n_grids = 4
    domain = [0.0, 1.0]
    levels = range(n_grids)
    lc_rel = [0.1 * 2 ** (-k) for k in levels]

    def no_halo_name(l):
        return f"grids/convergence/unit_cube_{l}/grid.geo"

    def with_halo_name(l):
        return f"grids/convergence/unit_cube_with_halo_{l}/grid.geo"

    gmsh.generate_cube_grids(no_halo_name, domain, lc_rel, levels, with_halo=False)
    gmsh.generate_cube_grids(with_halo_name, domain, lc_rel, levels, with_halo=True)


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
        geos = []

    gmsh.generate_grids(geos)
