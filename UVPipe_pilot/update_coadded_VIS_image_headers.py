#!/usr/bin/env python3

import numpy as np

from glob import glob
from astropy.io import fits


# Check if the input list contains the same values or not.
def check_for_consistency(value_list, input_name):
    if len(set(value_list)) == 1:
        value = value_list[0]
        return value
    else:
        print(f"\nDifferent {input_name} values found. Taking the mean.")
        value_list = np.array(value_list)
        unique, unique_counts = np.unique(value_list, return_counts=True)
        print(f"unique values = {unique}")
        print(f"unique value counts = {unique_counts}")
        max_to_min_difference = max(value_list) - min(value_list)
        print(f"Maximum - Minimum = {max_to_min_difference:.4}")
        value = np.mean(value_list)
        return value


print("\nUpdating the coadded VIS image headers.")

BZERO_MJD = []
BSCALE_MJD = []
UV_image_files = glob("uvt_*/[F,N]_*/*A_l2img.fits")
for UV_image_file in UV_image_files:
    hdu = fits.open(UV_image_file)
    BZERO_MJD.append(hdu[0].header["BZERO_MJD"])
    BSCALE_MJD.append(hdu[0].header["BSCALE_MJD"])

BZERO_MJD = check_for_consistency(BZERO_MJD, "BZERO_MJD")
BSCALE_MJD = check_for_consistency(BSCALE_MJD, "BSCALE_MJD")


def get_MJD_time(LOCAL_TIMEV=None, BZERO_MJD=BZERO_MJD, BSCALE_MJD=BSCALE_MJD):
    met = BZERO_MJD + BSCALE_MJD * LOCAL_TIMEV
    jan2010 = 55197.0  # 2010.0(UTC) expressed with MJD format and scale UTC.
    mjd = (met / 86400.0) + jan2010  # 1 julian day = 86400 seconds.
    return mjd


VIS_coadded_images = glob("uvt_*/V_*/*I_l2img.fits")

for coadded_image in VIS_coadded_images:
    drift_series_file = coadded_image.replace("I_l2img", "_l2dr")
    drift_series_hdu = fits.open(drift_series_file)
    time_data = drift_series_hdu[2].data["Time"]
    time_data = time_data[1:]
    exp_time = time_data[-1] - time_data[0]
    exp_start_time = get_MJD_time(time_data[0])
    exp_end_time = get_MJD_time(time_data[-1])
    coadded_image_hdu = fits.open(coadded_image, mode="update")
    header = coadded_image_hdu[0].header
    header.set("EXP_TIME", exp_time, "Total exposure time")
    header.set("EXPSTART", exp_start_time, "Exposure start time in MJD")
    header.set("EXPEND", exp_end_time, "Exposure end time in MJD")
    coadded_image_hdu.close()

print("Done!")
