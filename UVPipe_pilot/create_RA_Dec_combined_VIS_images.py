#!/usr/bin/env python3

import numpy as np

from glob import glob
from astropy.io import fits
from astropy.nddata import Cutout2D
from reproject import reproject_exact
from reproject.mosaicking import find_optimal_celestial_wcs


def reproject_image(filename):
    hdu = fits.open(filename)
    wcs_out, shape_out = find_optimal_celestial_wcs([hdu[0]], frame="icrs")
    array, footprint = reproject_exact(hdu[0], wcs_out, shape_out, parallel=True)
    array[footprint < 0.5] = 0
    position = np.array(shape_out) / 2
    cutout = Cutout2D(
        array, position=position.astype("int"), size=(1800, 1800), wcs=wcs_out
    )
    data = cutout.data.astype(np.single)

    if np.isfinite(data).all() is False:
        raise ValueError("The array contains non-finite values.")

    new_hdu = fits.PrimaryHDU(data, header=cutout.wcs.to_header())
    new_hdu.header["N_EPISDS"] = hdu[0].header["N_EPISDS"]
    new_hdu.header["EXPSTART"] = hdu[0].header["EXPSTART"]
    new_hdu.header["EXPEND"] = hdu[0].header["EXPEND"]
    new_hdu.header["EXP_TIME"] = hdu[0].header["EXP_TIME"]
    new_hdu.writeto("RA_Dec_" + filename, overwrite=True)


combined_VIS_images = glob("combined_VIS_IM*.fits")

if __name__ == "__main__":
    for combined_VIS_image in combined_VIS_images:
        channel, window_filter = combined_VIS_image.split("_")[1:3]
        window_filter = window_filter[:6]
        print(f"Reprojecting {channel} {window_filter} combined VIS image.")
        reproject_image(combined_VIS_image)
