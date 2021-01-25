#! /usr/bin/env python3

import os
import numpy as np
import h5py


def load_array(fn, key):
    with h5py.File(fn, "r") as h5:
        return np.array(h5[key])


def parse_part(base_dir, n_parts, k_part):
    gn = os.path.join(base_dir, str(n_parts), f"subgrid-{k_part:04d}.msh.h5")
    partition = load_array(gn, "partition")

    n_cells_local = np.count_nonzero(partition == k_part)
    return n_cells_local


base_dir = "/mnt/euler/zisa-grids/scaling-2/partitioned"
n_parts = 128 * 4

local_cells = [parse_part(base_dir, n_parts, k) for k in range(n_parts)]
print(local_cells)
print(np.min(local_cells) / np.max(local_cells))
