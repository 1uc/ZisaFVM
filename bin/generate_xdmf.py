#! /usr/bin/env python3
import numpy as np
import h5py
import argparse

import matplotlib.pyplot as plt

import tiwaz.xdmf as xdmf
from tiwaz.tri_plot import TriPlot
from tiwaz.post_process import load_grid


if __name__ == "__main__":
    parser_help = "Generate XDMF file."
    parser = argparse.ArgumentParser(description=parser_help)

    parser.add_argument(
        "-o", "--output",
        help="Name of XDMF file.",
    )

    parser.add_argument(
        "--grid",
        help="Name of the HDF5 grid.",
    )

    parser.add_argument(
        "data_files",
        nargs="+")


    args = parser.parse_args()

    data_files = args.data_files
    grid_file = args.grid
    xdmf_file = args.output

    components = ["rho", "E"]

    with open(xdmf_file, "w") as f:
        f.write(xdmf.xml_to_string(xdmf.generate_xdmf(grid_file, data_files, components)))



