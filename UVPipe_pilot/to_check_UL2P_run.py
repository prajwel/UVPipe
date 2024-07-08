#!/usr/bin/env python3

import re
import os
import numpy as np
import matplotlib.pyplot as plt

from glob import glob
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from matplotlib.colors import LogNorm


def rebin(arr, new_shape):
    shape = (
        new_shape[0],
        arr.shape[0] // new_shape[0],
        new_shape[1],
        arr.shape[1] // new_shape[1],
    )
    return arr.reshape(shape).mean(-1).mean(1)


def tile_array(arr, factor):
    m, n = arr.shape
    return np.repeat(np.repeat(arr, factor, axis=1), factor, axis=0).reshape(
        m * factor, n * factor
    )


params_check_list = ["framesDiscardpc,i,h,2", "framesDiscardpcfuv,i,h,2"]


def check_pipeline_parameters():
    param_file = glob("pipeline/*")

    if len(param_file) == 1:
        param_file = param_file[0]
    elif len(param_file) == 0:
        raise RuntimeError("No L2 parameter file!")
    else:
        raise RuntimeError(f"Multiple parameter files should not exist: {param_file}")

    for d in params_check_list:
        with open(param_file) as f:
            if d not in f.read():
                raise RuntimeError(f"Unexpected L2 parameter file settings: {d}")

    print("L2 parameters are OK.")


check_pipeline_parameters()

# original resolution is 4800. change it to factors of 2.
resolution = 600

root = "."
for dirpath, dirnames, files in os.walk(root):
    for s in files:
        fnam = dirpath + "/" + s
        if fnam[-12:] in ["I_l2img.fits"]:
            hdu = fits.open(fnam)
            path = os.path.normpath(fnam)
            parent = path.split(os.sep)[-2]
            plt.figure(figsize=(10, 10))
            if parent[0] == "V":
                data = hdu[0].data
                mean, median, std = sigma_clipped_stats(data, sigma=5.0, maxiters=1)
                plt.imshow(data, cmap="hot", origin="lower", vmax=median + 5 * std)
            else:
                data = rebin(hdu[0].data, [resolution, resolution])
                plt.imshow(data, cmap="hot", origin="lower", norm=LogNorm())
            filename = path.split(os.sep)[-1]
            figure_name = filename.replace(".fits", "_quick_look.png")
            figure_name = parent + "_" + figure_name
            plt.savefig(figure_name, format="png", bbox_inches="tight", dpi=150)
            plt.close()

        if fnam[-9:] in ["l2dr.fits"]:
            hdu = fits.open(fnam)
            history = hdu[0].header["history"]
            re_result = re.search(r"V/uvtV\.\d{2}/", str(history))
            ras_num = re_result.group()[5:9]
            fig, ax = plt.subplots()
            ax.plot(hdu[2].data["Time"], hdu[2].data["X_shift"], label="X_shift")
            ax.plot(hdu[2].data["Time"], hdu[2].data["Y_shift"], label="Y_shift")
            ax.axhline(44)
            ax.axhline(-44)
            ax.set_xlabel("Time")
            ax.set_ylabel("Pixel")
            ax.set_title(ras_num)
            ax.legend()
            path = os.path.normpath(fnam)
            parent = path.split(os.sep)[-2]
            filename = path.split(os.sep)[-1]
            figure_name = filename.replace(".fits", "_drift.png")
            figure_name = parent + "_" + figure_name
            plt.savefig(figure_name, format="png", bbox_inches="tight", dpi=150)
            plt.close()
