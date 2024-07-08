#!/usr/bin/env python3

import os
import warnings
import numpy as np
import subprocess

from glob import glob
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs.utils import proj_plane_pixel_scales


def get_average_plate_scale(channel, window_filter):
    XP_FUV = 0.4168705
    YP_FUV = 0.416848964

    XP_NUV = 0.417310
    YP_NUV = 0.417218

    XP_NUV_F2 = 0.419686
    YP_NUV_F2 = 0.419184

    XP_NUV_F1 = 0.4173677
    YP_NUV_F1 = 0.417299

    if channel == "FUV":
        plate_scales = [XP_FUV, YP_FUV]
    elif channel == "NUV":
        if window_filter[-2:] == "F2":
            plate_scales = [XP_NUV_F2, YP_NUV_F2]
        elif window_filter[-2:] == "F1":
            plate_scales = [XP_NUV_F1, YP_NUV_F1]
        else:
            plate_scales = [XP_NUV, YP_NUV]
    else:
        raise ValueError

    average_plate_scale = np.average(plate_scales)
    return average_plate_scale


def remove_warnings(text):
    text = text.split(sep="\n")
    modified_text = []
    for line in text:
        if "warning" not in line.lower():
            modified_text.append(line)
    result = "\n".join(modified_text)
    return result


def remove_useless_lines(text):
    text = text.split(sep="\n")
    modified_text = []
    for line in text:
        if "did not solve." not in line.lower():
            modified_text.append(line)
    result = "\n".join(modified_text)
    return result


def run_solve_field(args):
    output = subprocess.run(args, capture_output=True)
    modified_stdout = remove_warnings(output.stdout.decode())
    if output.returncode == 0:
        # print('\nsolve-field run output:')
        # print(modified_stdout)
        if "no WCS file was written" in modified_stdout:
            print("\nsolve-field run output:")
            print(remove_useless_lines(modified_stdout))
            raise RuntimeError("solve-field failed.")
    else:
        print(modified_stdout)
        raise RuntimeError("solve-field failed.")


def call_solve_field(
    coo_file,
    corr_file,
    wcs_file,
    RA_pointing,
    DEC_pointing,
    plate_scale_low,
    plate_scale_high,
):
    xy_file = coo_file.replace(".dat", ".xy")

    args = ["text2fits", "-f", "dd", "-s", ",", "-H", "x,y", coo_file, xy_file]

    text2fits_output = subprocess.run(args, capture_output=True)

    if text2fits_output.returncode == 0:
        pass
    else:
        raise RuntimeError("text2fits failed")

    intermediate_wcs_file = "intermediate.wcs"

    args = [
        "solve-field",
        "--config",
        "astrometry.cfg",
        "--no-plots",
        "--no-remove-lines",
        "--overwrite",
        "--scale-units",
        "arcsecperpix",
        "--scale-low",
        str(plate_scale_low),
        "--scale-high",
        str(plate_scale_high),
        "--ra",
        str(RA_pointing),
        "--dec",
        str(DEC_pointing),
        "--radius",
        "1.0",
        "--width",
        "4800",
        "--height",
        "4800",
        "--cpulimit",
        "1200",
        "--crpix-center",
        "--no-tweak",
        "--uniformize",
        "0",
        "--wcs",
        intermediate_wcs_file,
        "--rdls",
        "none",
        "--match",
        "none",
        "--index-xyls",
        "none",
        "--corr",
        corr_file,
        "--temp-axy",
        "--solved",
        "none",
        xy_file,
        "--odds-to-tune-up",
        "1e4",
        "--odds-to-solve",
        "1e6",
    ]

    run_solve_field(args)

    args = [
        "solve-field",
        "--config",
        "astrometry.cfg",
        "--no-plots",
        "--no-remove-lines",
        "--overwrite",
        "--scale-units",
        "arcsecperpix",
        "--scale-low",
        str(plate_scale_low),
        "--scale-high",
        str(plate_scale_high),
        "--ra",
        str(RA_pointing),
        "--dec",
        str(DEC_pointing),
        "--radius",
        "1.0",
        "--width",
        "4800",
        "--height",
        "4800",
        "--cpulimit",
        "1200",
        "--crpix-center",
        "--no-tweak",
        "--uniformize",
        "0",
        "--verify",
        intermediate_wcs_file,
        "--wcs",
        wcs_file,
        "--rdls",
        "none",
        "--match",
        "none",
        "--index-xyls",
        "none",
        "--corr",
        corr_file,
        "--temp-axy",
        "--solved",
        "none",
        "--pixel-error",
        "1",
        "--no-verify-uniformize",
        xy_file,
    ]

    run_solve_field(args)

    os.remove(intermediate_wcs_file)
    os.remove(xy_file)

    print("solve-field run is completed.")
    return wcs_file, corr_file


def uv_astrometry(VIS_WCS_file=None, error_factor=None):
    channel, window_filter = VIS_WCS_file.split("_")[-3:-1]
    coo_file = f"combined_{channel}_{window_filter}_CPS_coo.dat"

    VIS_WCS_hdu = fits.open(VIS_WCS_file)
    RA_pointing = VIS_WCS_hdu[0].header["CRVAL1"]
    DEC_pointing = VIS_WCS_hdu[0].header["CRVAL2"]

    plate_scale = get_average_plate_scale(channel, window_filter)
    plate_scale_error = error_factor * plate_scale
    plate_scale_low = plate_scale - plate_scale_error
    plate_scale_high = plate_scale + plate_scale_error

    print("\nInputs are the following:")
    print(f"input coo file = {coo_file}")
    print(f"RA = {RA_pointing}")
    print(f"DEC = {DEC_pointing}")
    print(f"plate scale = {plate_scale:.6f}''")
    print(f"plate scale (low) = {plate_scale_low:.6f}''")
    print(f"plate scale (high) = {plate_scale_high:.6f}''")

    corr_file = coo_file.replace(".dat", ".corr")
    wcs_file = coo_file.replace(".dat", ".wcs")

    wcs_filename, corr_filename = call_solve_field(
        coo_file,
        corr_file,
        wcs_file,
        RA_pointing,
        DEC_pointing,
        plate_scale_low,
        plate_scale_high,
    )

    corr_hdu = fits.open(corr_filename, mode="update")

    coo_data = np.genfromtxt(coo_file, delimiter=",")
    coo_x, coo_y = coo_data.T

    # While the corr file field values have double precision,
    # the values from the coo_file will be kept in the corr file.
    matched_indices = corr_hdu[1].data["field_id"]
    field_x = coo_x[matched_indices]
    field_y = coo_y[matched_indices]
    corr_hdu[1].data["field_x"] = field_x
    corr_hdu[1].data["field_y"] = field_y
    corr_hdu.close()

    w_astrometry_net = fits.open(wcs_file)
    wcs = WCS(w_astrometry_net[0].header)
    pixel_scales = proj_plane_pixel_scales(wcs) * 3600
    pixel_scale = pixel_scales[0]
    print(f"WCS plate scale = {pixel_scale:.6f}''")


warnings.simplefilter("ignore", category=AstropyWarning)

VIS_WCS_files = glob("UV_oriented_combined_VIS_for_*_coo.wcs")
VIS_WCS_files.sort(reverse=True)
for VIS_WCS_file in VIS_WCS_files:
    try:
        uv_astrometry(VIS_WCS_file, 10 / 4000)
    except RuntimeError:
        try:
            uv_astrometry(VIS_WCS_file, 20 / 4000)
        except RuntimeError:
            channel, window_filter = VIS_WCS_file.split("_")[-3:-1]
            print(f"UV astrometry failed for {channel} {window_filter}")
