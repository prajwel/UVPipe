#!/usr/bin/env python3

from glob import glob
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS


def apply_wcs(fits_file, wcs_filename, method):
    wcs_hdu = fits.open(wcs_filename)
    wcs_hdu[0].header["NAXIS"] = 2
    wcs = WCS(wcs_hdu[0].header)
    hdu = fits.open(fits_file, mode="update")
    hdu[0].header.update(wcs.to_header(), "update")
    if method:
        hdu[0].header["ASTMETH"] = method
    hdu.close()


CPS_images = glob("combined_*CPS.fits")
for CPS_image in CPS_images:
    channel, window_filter = CPS_image.split("_")[1:3]
    print(f"Applying WCS solution to {channel} {window_filter} combined products.")
    wcs_filename = f"combined_{channel}_{window_filter}_CPS_coo_fixed_platescale.wcs"

    if Path(wcs_filename).is_file():
        method = "Direct"
    else:
        print(
            f"No UV astrometry, applying VIS astrometry to {channel} {window_filter}."
        )
        method = "Indirect"
        wcs_filename = f"UV_oriented_combined_VIS_for_{channel}_{window_filter}_coo_fixed_platescale.wcs"

    CPS_error_image = CPS_image[:-5] + "_error.fits"
    EXP_MAP_image = CPS_image[:-9] + "_exp_map_in_seconds.fits"
    apply_wcs(CPS_image, wcs_filename, method)
    apply_wcs(CPS_error_image, wcs_filename, method)
    apply_wcs(EXP_MAP_image, wcs_filename, method)

combined_VIS_images = glob("combined_VIS_IM*.fits")
for combined_VIS_image in combined_VIS_images:
    channel, window_filter = combined_VIS_image.split("_")[1:3]
    window_filter = window_filter[:6]
    print(f"Applying WCS solution to {channel} {window_filter} combined VIS image.")
    wcs_filename = combined_VIS_image[:-5] + "_coo.wcs"
    apply_wcs(combined_VIS_image, wcs_filename, None)
