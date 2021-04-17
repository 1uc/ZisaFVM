#!/usr/bin/env python3

import os
import sys

import h5py
import numpy as np
import argparse

import tiwaz
from tiwaz.post_process import load_grid


if __name__ == "__main__":
    parser_help = "Visualize the output Zisa."
    parser = argparse.ArgumentParser(description=parser_help)

    parser.add_argument(
        "--grid", help="The filename of the grid.",
    )

    args = parser.parse_args()

    grid = load_grid(args.grid)

    vol = grid.volumes
    i_small = sorted(np.arange(vol.shape[0]), key=lambda j: vol[j])[:10]

    I = np.zeros(vol.shape[0])
    I[i_small] = 1.0

    with h5py.File("indices.h5", "w") as h5:
        h5["indices"] = I
