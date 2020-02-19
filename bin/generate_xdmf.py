#! /usr/bin/env python3

import os
import numpy as np
import h5py
import argparse

import matplotlib.pyplot as plt

import tiwaz.xdmf as xdmf
from tiwaz.tri_plot import TriPlot
from tiwaz.post_process import load_grid
from tiwaz.post_process import find_grid
from tiwaz.post_process import find_data_files


if __name__ == "__main__":
    parser_help = "Generate XDMF file."
    parser = argparse.ArgumentParser(description=parser_help)

    parser.add_argument(
        "-o", "--output", help="Name of XDMF file.",
    )

    parser.add_argument(
        "--grid", help="Name of the HDF5 grid.",
    )

    parser.add_argument(
        "data_files",
        nargs="+",
        help="Either a list of filenames to visualize or a Zisa output directory.",
    )

    args = parser.parse_args()

    if len(args.data_files) == 1 and os.path.isdir(args.data_files[0]):
        base_directory = args.data_files[0]
        data_files = find_data_files(base_directory)
    else:
        data_files = args.data_files

    grid_file = args.grid or find_grid(base_directory)

    default_xdmf_file = data_files[0][: -len("_data-????.h5")] + ".xdmf"
    xdmf_file = args.output or default_xdmf_file

    components = ["rho", "E"]

    with open(xdmf_file, "w") as f:
        f.write(
            xdmf.xml_to_string(xdmf.generate_xdmf(grid_file, data_files, components))
        )
