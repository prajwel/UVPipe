#!/usr/bin/env python3

import os
import shutil
import numpy as np

from astropy.io import fits

# Window Vs Frame rate dictionary.
window_rate_dict = {
    "511": 28.719,
    "349": 59.82,
    "299": 80.29,
    "249": 113.42,
    "199": 172.30,
    "149": 292.76,
    "99": 605.08,
}

# The threshold is defined in frames.
exp_in_frames_threshold = 60 * 28.7185

r = []
root = "."
for dirpath, dirnames, files in os.walk(root):
    for s in files:
        fnam = dirpath + "/" + s
        if fnam[-12:] == "I_l2exp.fits" and os.path.dirname(fnam)[-2:] != "ci":
            hdulist = fits.open(fnam)
            window = str(hdulist[0].header["win_x_sz"])
            framecount_per_sec = float(window_rate_dict[window])
            exp_map = hdulist[0].data
            exp_in_frames = np.round(np.max(exp_map))
            if exp_in_frames < exp_in_frames_threshold:
                dir_to_remove = os.path.dirname(fnam)
                exp_time = np.round(exp_in_frames / framecount_per_sec).astype(int)
                shutil.rmtree(dir_to_remove)
                os.mkdir(dir_to_remove)
                print(
                    f"\nRemoving {dir_to_remove} due to low exposure time ({exp_time} seconds)"
                )
