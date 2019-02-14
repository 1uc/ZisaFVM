#! /usr/bin/env python3
# encoding: utf-8

import argparse
import subprocess
import glob
import multiprocessing

gmsh = "bin/gmsh"

def run_gmesh(geo):
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

def all_geo_files(path):
    return glob.glob("{}/*.geo".format(path)) + glob.glob("{}/**/*.geo".format(path))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generate .msh files.")

    parser.add_argument('--minimal',
                        action='store_true',
                        help="Generate only those grids required for testing.")

    args = parser.parse_args()

    if args.minimal:
        geos = minimal_geo_files("grids")
    else:
        geos = all_geo_files("grids")

    with multiprocessing.Pool() as pool:
        for i, geo in enumerate(pool.imap_unordered(run_gmesh, geos)):
            print("[{:2d}/{}] {}".format(i, len(geos), geo))
