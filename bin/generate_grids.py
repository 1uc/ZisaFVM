#! /usr/bin/env python3
# encoding: utf-8

import argparse
import glob

import tiwaz
import tiwaz.gmsh as gmsh

def minimal_geo_files(path):
    return [
        path + "/dbg.geo",
        path + "/small.geo"
    ] + [
        path + "/convergence/unit_square_{:d}.geo".format(d) for d in range(4)
    ]

def matching_geo_files(pattern):
    return glob.glob("{}*.geo".format(pattern))

def all_geo_files(path):
    return glob.glob("{}/*.geo".format(path)) + glob.glob("{}/**/*.geo".format(path))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate .msh files.")

    parser.add_argument('--minimal',
                        action='store_true',
                        help="Generate only those grids required for testing.")

    parser.add_argument("pattern",
                        nargs="?",
                        default=None,
                        help="Generate only for the matching grids.")

    args = parser.parse_args()

    if args.minimal:
        geos = minimal_geo_files("grids")
    elif args.pattern:
        geos = matching_geo_files(args.pattern)
    else:
        geos = all_geo_files("grids")

    gmsh.generate_grids(geos)
