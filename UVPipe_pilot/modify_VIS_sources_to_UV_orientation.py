#!/usr/bin/env python3

import numpy as np

from glob import glob


def UV_orient(coo_file, channel):
    VIS_to_FUV = np.array(
        [
            [0.83988, -0.57725, 146.13500],
            [0.57725, 0.83988, -109.14374],
            [0.0, 0.0, 1.0],
        ]
    )
    VIS_to_NUV = np.array(
        [
            [-1.01474, 0.04521, -17.28794],
            [-0.04521, -1.01474, -12.51485],
            [0.0, 0.0, 1.0],
        ]
    )

    VIS_to_UV_matrix = {"FUV": VIS_to_FUV, "NUV": VIS_to_NUV}
    VIS_to_UV = VIS_to_UV_matrix[channel]

    positions = np.genfromtxt(coo_file, delimiter=",")
    xcentroid = positions[:, 0]
    ycentroid = positions[:, 1]

    xcentroid = xcentroid * (8 / 3)
    ycentroid = ycentroid * (8 / 3)

    xcentroid = xcentroid - 2400
    ycentroid = ycentroid - 2400

    xcentroid_UV_orientation = (
        xcentroid * VIS_to_UV[0, 0] + ycentroid * VIS_to_UV[0, 1] + VIS_to_UV[0, 2]
    )
    ycentroid_UV_orientation = (
        xcentroid * VIS_to_UV[1, 0] + ycentroid * VIS_to_UV[1, 1] + VIS_to_UV[1, 2]
    )

    xcentroid_UV_orientation = xcentroid_UV_orientation + 2400
    ycentroid_UV_orientation = ycentroid_UV_orientation + 2400

    new_positions = np.transpose((xcentroid_UV_orientation, ycentroid_UV_orientation))

    UV_oriented_coo_file = "UV_oriented_" + coo_file
    np.savetxt(UV_oriented_coo_file, new_positions, fmt="%s", delimiter=",")


coo_files = glob("combined_VIS_for_*_coo.dat")
for coo_file in coo_files:
    channel = coo_file.split("_")[3]
    print(f"\nModifying {coo_file} to {channel} orientation.")
    UV_orient(coo_file, channel)
