#!/usr/bin/env python3

import os
import warnings
import subprocess

from glob import glob
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs.utils import proj_plane_pixel_scales


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


def add_radec_to_args(args, RA_pointing, DEC_pointing):
    if RA_pointing is not None and DEC_pointing is not None:
        ra_dec_args = [
            "--ra",
            str(RA_pointing),
            "--dec",
            str(DEC_pointing),
            "--radius",
            "1.0",
        ]
        args = args + ra_dec_args
    return args


def call_solve_field(
    coo_file,
    corr_file,
    wcs_file,
    RA_pointing,
    DEC_pointing,
    plate_scale_low,
    plate_scale_high,
    size_x,
    size_y,
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
        "--width",
        str(size_x),
        "--height",
        str(size_y),
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
    ]

    args = add_radec_to_args(args, RA_pointing, DEC_pointing)
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
        "--width",
        str(size_x),
        "--height",
        str(size_y),
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

    args = add_radec_to_args(args, RA_pointing, DEC_pointing)
    run_solve_field(args)

    os.remove(intermediate_wcs_file)
    os.remove(xy_file)

    print("solve-field run is completed.")
    return wcs_file, corr_file


def vis_astrometry(coo_file, plate_scale_low, plate_scale_high, size_x, size_y):
    VIS_image = coo_file.replace("_coo.dat", ".fits")
    VIS_image = VIS_image.replace("UV_oriented_", "")
    vis_hdu = fits.open(VIS_image)
    RA_pointing = vis_hdu[0].header["RA_PNT"]
    DEC_pointing = vis_hdu[0].header["DEC_PNT"]

    print("\nInputs are the following:")
    print(f"input coo file = {coo_file}")
    print(f"RA = {RA_pointing}")
    print(f"DEC = {DEC_pointing}")
    print(f"plate scale (low) = {plate_scale_low:.6f}''")
    print(f"plate scale (high) = {plate_scale_high:.6f}''")

    corr_file = coo_file.replace(".dat", ".corr")
    wcs_file = coo_file.replace(".dat", ".wcs")

    try:
        wcs_filename, corr_filename = call_solve_field(
            coo_file,
            corr_file,
            wcs_file,
            RA_pointing,
            DEC_pointing,
            plate_scale_low,
            plate_scale_high,
            size_x,
            size_y,
        )
    except RuntimeError:
        print(f"solve-field run with {RA_pointing}, {DEC_pointing} failed.")
        print("Trying full sky now.")
        wcs_filename, corr_filename = call_solve_field(
            coo_file,
            corr_file,
            wcs_file,
            None,
            None,
            plate_scale_low,
            plate_scale_high,
            size_x,
            size_y,
        )

    w_astrometry_net = fits.open(wcs_file)
    wcs = WCS(w_astrometry_net[0].header)
    pixel_scales = proj_plane_pixel_scales(wcs) * 3600
    pixel_scale = pixel_scales[0]
    print(f"WCS plate scale = {pixel_scale:.6f}''")


warnings.simplefilter("ignore", category=AstropyWarning)

coo_files = glob("UV_oriented_combined_VIS_for*_coo.dat")
for coo_file in coo_files:
    vis_astrometry(coo_file, 0.40, 0.43, 4800, 4800)

coo_files = glob("combined_VIS_IM*_coo.dat")
for coo_file in coo_files:
    vis_astrometry(coo_file, 1.07, 1.15, 1800, 1800)
