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
    raise_error_flag = False
    hdu = fits.open(filename)
    exp_time = hdu[0].header["EXP_TIME"]
    calculated_framerate = 1 / hdu[0].header["INT_TIME"]
    window_size = str(hdu[0].header["WIN_X_SZ"])
    expected_framerate = window_rate_dict[window_size]

    diff_percentage = (
        (calculated_framerate - expected_framerate) * 100 / expected_framerate
    )
    abs_diff_percentage = abs(diff_percentage)

    if abs_diff_percentage > tolerance_in_percentage:
        print(
            f"""\nDetected a file with an outlier framerate value:
        file = {filename},
        Exposure time (s) = {exp_time},
        Tolerance (%) = {tolerance_in_percentage},
        Expected framerate = {expected_framerate},
        Calculated framerate = {calculated_framerate}
        Difference (%) = {abs_diff_percentage}\n"""
        )
        raise_error_flag = True
    return raise_error_flag


files_list = glob("uvt_*/[F,N]*/*I_l2img.fits")
flags = []
print("The following files have been checked:")
for filename in files_list:
    print(filename)
    flags.append(check_framerate(filename))

if True in flags:
    raise RuntimeError("Please check episodes with outlier framerate values.")
