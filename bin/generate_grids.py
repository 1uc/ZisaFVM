#! /usr/bin/env python3
# encoding: utf-8

import argparse
import subprocess
import glob
import multiprocessing


def run_gmesh(geo):
    gmsh = "bin/gmsh"
    cmd = [gmsh, "-2", geo]
    output = subprocess.check_call(cmd,
                                   stdout=subprocess.DEVNULL,
                                   stderr=subprocess.DEVNULL)

    return geo


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
                        default=None,
                        help="Generate only for the matching grids.")

    args = parser.parse_args()

    if args.minimal:
        geos = minimal_geo_files("grids")
    elif args.pattern:
        geos = matching_geo_files(args.pattern)
    else:
        geos = all_geo_files("grids")

    with multiprocessing.Pool() as pool:
        for i, geo in enumerate(pool.imap_unordered(run_gmesh, geos)):
            print("[{:2d}/{}] {}".format(i+1, len(geos), geo))
