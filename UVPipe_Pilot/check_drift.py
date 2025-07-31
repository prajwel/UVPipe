import re
import os
import numpy as np

from astropy.io import fits


def find_drift_anomalies(drift_series, T=1, T_rev=1):
    # Compute all jumps
    jump = np.diff(drift_series)

    # Define jump segments relative to index n:
    j0 = jump[:-3]
    j1 = jump[1:-2]
    j2 = jump[2:-1]
    j3 = jump[3:]

    # Candidate middle indices n, offset by +2 due to slicing
    n_indices = np.arange(2, len(drift_series) - 2)

    # Conditions:
    jump_cond = (np.abs(j1) > T) & (np.abs(j2) > T_rev) & (np.sign(j1) != np.sign(j2))

    local_smooth = (0.5 * j1 <= j0) & (j0 <= j1) & (0.5 * j2 <= j3) & (j3 <= j2)

    # Final condition: jump–reverse-jump not explained by neighbourhood
    suspect_mask = jump_cond & ~local_smooth
    suspect_indices = n_indices[suspect_mask]
    return suspect_indices


root = "."
for dirpath, dirnames, files in os.walk(root):
    for s in files:
        fnam = dirpath + "/" + s
        if fnam[-9:] in ["l2dr.fits"]:
            print(f"Checking {fnam}")
            hdu = fits.open(fnam)
            history = hdu[0].header["history"]
            re_result = re.search(r"V/uvtV\.\d{2}/", str(history))
            ras_num = re_result.group()[5:9]
            time = hdu[2].data["Time"]
            X_shift = hdu[2].data["X_shift"]
            Y_shift = hdu[2].data["Y_shift"]
            suspect_X_shift_indices = find_drift_anomalies(X_shift)
            suspect_Y_shift_indices = find_drift_anomalies(Y_shift)
            if len(suspect_X_shift_indices) != 0:
                raise RuntimeError(f"Check {suspect_X_shift_indices} in X_shift")
            if len(suspect_Y_shift_indices) != 0:
                raise RuntimeError(f"Check {suspect_Y_shift_indices} in Y_shift")
