#!/usr/bin/env python3

from glob import glob
from astropy.io import fits

tolerance_in_percentage = 0.1

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


def check_framerate(filename):
    hdu = fits.open(filename)
    calculated_framerate = 1 / hdu[0].header["INT_TIME"]
    window_size = str(hdu[0].header["WIN_X_SZ"])
    expected_framerate = window_rate_dict[window_size]

    diff_percentage = (
        (calculated_framerate - expected_framerate) * 100 / expected_framerate
    )
    abs_diff_percentage = abs(diff_percentage)

    if abs_diff_percentage > tolerance_in_percentage:
        raise RuntimeError(
            f"""Detected a file having a outlier framerate value:
        file = {filename},
        Tolerance (%) = {tolerance_in_percentage},
        Expected framerate = {expected_framerate},
        Calculated framrate = {calculated_framerate}"""
        )


files_list = glob("uvt_*/*/*A_l2img.fits")
for filename in files_list:
    check_framerate(filename)
