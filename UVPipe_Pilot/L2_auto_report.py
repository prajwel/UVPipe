#!/usr/bin/env python3


import os
import numpy as np
from astropy.io import fits
from itertools import zip_longest


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

r = []
root = "."
for dirpath, dirnames, files in os.walk(root):
    for s in files:
        fnam = dirpath + "/" + s
        if fnam[-12:] == "I_l2exp.fits" and os.path.dirname(fnam)[-2:] != "ci":
            hdulist = fits.open(fnam)
            filterr = hdulist[0].header["filter"]
            exp_time2 = hdulist[0].header["exp_time"]
            exp_time2 = np.around(exp_time2)
            window = str(hdulist[0].header["win_x_sz"])
            framecount_per_sec = float(window_rate_dict[window])
            exp_map = hdulist[0].data
            shape = exp_map.shape[0]
            x, y = np.ogrid[:shape, :shape]
            cx, cy = shape / 2.0, shape / 2.0
            radius = shape / 4
            r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy)
            circmask = r2 <= radius * radius
            exp_time = np.max(exp_map[circmask]) / framecount_per_sec
            exp_time = np.around(exp_time)
            r.append([fnam[-53:-49], exp_time, exp_time2, filterr])

uvtN = [d[0:4] for d in r if d[0][:1] == "N"]
uvtN = sorted(uvtN)
uvtF = [d[0:4] for d in r if d[0][:1] == "F"]
uvtF = sorted(uvtF)

obsID = s[3:24]
L2_exptime_filename = "L2_stats_" + obsID + ".csv"
stats_header = "#NUV_dir, exptime, exptime (header), filter, FUV_dir, exptime, exptime (header), filter\n"

with open(L2_exptime_filename, "wt") as l2file:
    l2file.write(stats_header)
    for combined_row in zip_longest(uvtN, uvtF, fillvalue=np.array(["", "", "", ""])):
        l2file.write(
            "{z[0]},{z[1]},{z[2]},{z[3]},{z[4]},{z[5]},{z[6]},{z[7]}\n".format(
                z=np.concatenate(combined_row)
            )
        )

print("\nDone!\n")
