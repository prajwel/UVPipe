#!/usr/bin/env python

import warnings
import numpy as np

from glob import glob
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.exceptions import AstropyWarning


def get_header_values(combined_eventslist_hdu):
    framecount_per_sec = combined_eventslist_hdu[0].header["AVGFRMRT"]
    framecount_per_sec_stdev = combined_eventslist_hdu[0].header["STDFRMRT"]
    RA_pointing = combined_eventslist_hdu[0].header["RA_PNT"]
    DEC_pointing = combined_eventslist_hdu[0].header["DEC_PNT"]
    return framecount_per_sec, framecount_per_sec_stdev, RA_pointing, DEC_pointing


def read_columns(events_list_hdu):
    # Reading few columns.
    time = events_list_hdu[1].data["MJD_L2"]
    fx = events_list_hdu[1].data["Fx"]
    fy = events_list_hdu[1].data["Fy"]
    photons = events_list_hdu[1].data["EFFECTIVE_NUM_PHOTONS"]
    bad_flag = events_list_hdu[1].data["BAD FLAG"]
    mask = photons > 0
    mask = np.logical_and(mask, bad_flag)
    time = time[mask]
    fx = fx[mask]
    fy = fy[mask]
    photons = photons[mask]
    return time, fx, fy, photons


def instrument_to_astronomical(fx, fy, PC1_1, PC1_2, PC2_1, PC2_2):
    fx_prime = (PC1_1 * (fx - 2400)) + (PC1_2 * (fy - 2400)) + 2400
    fy_prime = (PC2_1 * (fx - 2400)) + (PC2_2 * (fy - 2400)) + 2400
    return fx_prime, fy_prime


def create_exp_map_in_seconds(
    exp_map_hdu, combined_eventslist_hdu, channel, window_filter, method
):
    (
        framecount_per_sec,
        framecount_per_sec_stdev,
        RA_pointing,
        DEC_pointing,
    ) = get_header_values(combined_eventslist_hdu)
    ten_percent_limit = 0.1 * np.max(exp_map_hdu[0].data)
    exp_map_hdu[0].data[exp_map_hdu[0].data < ten_percent_limit] = 0
    data = exp_map_hdu[0].data / framecount_per_sec
    data = data.astype(np.single)
    hdu = fits.PrimaryHDU(data=data, header=exp_map_hdu[0].header)
    hdu.header["RA_PNT"] = RA_pointing
    hdu.header["DEC_PNT"] = DEC_pointing
    hdu.header["AVGFRMRT"] = framecount_per_sec
    hdu.header["STDFRMRT"] = framecount_per_sec_stdev
    hdu.header["ASTMETH"] = method
    expmap_in_sec_filename = (
        f"RA_Dec_combined_{channel}_{window_filter}_exp_map_in_seconds.fits"
    )
    hdu.writeto(expmap_in_sec_filename, overwrite=True)
    print(f"\nCombined exposure map in seconds: {expmap_in_sec_filename}")


def create_CPS_and_error_images(
    exp_map_hdu, combined_eventslist_hdu, channel, window_filter, wcs_filename, method
):
    (
        framecount_per_sec,
        framecount_per_sec_stdev,
        RA_pointing,
        DEC_pointing,
    ) = get_header_values(combined_eventslist_hdu)

    combined_eventslist_hdu[1].data.sort(order="MJD_L2")

    if "DEC-RADEC" in combined_eventslist_hdu[1].data.columns.names:
        combined_eventslist_hdu[1].data.columns["RA"].name = "Sky_RA"
        combined_eventslist_hdu[1].data.columns["DEC"].name = "Sky_DEC"
        combined_eventslist_hdu[1].data.columns["RA-RADEC"].name = "Fx_astronomical"
        combined_eventslist_hdu[1].data.columns["DEC-RADEC"].name = "Fy_astronomical"

    wcs_hdu = fits.open(wcs_filename)
    PC1_1 = wcs_hdu[0].header["PC1_1"]
    PC1_2 = wcs_hdu[0].header["PC1_2"]
    PC2_1 = wcs_hdu[0].header["PC2_1"]
    PC2_2 = wcs_hdu[0].header["PC2_2"]

    CRPIX1 = wcs_hdu[0].header["CRPIX1"]
    CRPIX2 = wcs_hdu[0].header["CRPIX2"]

    CDELT1 = wcs_hdu[0].header["CDELT1"]
    CDELT2 = wcs_hdu[0].header["CDELT2"]

    if CRPIX1 == CRPIX2 == 2400.5:
        pass
    else:
        raise ValueError("Check WCS CRPIX values!")

    if CDELT1 > 0:
        PC1_1 = -1 * PC1_1
        PC1_2 = -1 * PC1_2
        CDELT1 = -1 * CDELT1

    if CDELT2 < 0:
        PC2_1 = -1 * PC2_1
        PC2_2 = -1 * PC2_2
        CDELT2 = -1 * CDELT2

    ten_percent_limit = 0.1 * np.max(exp_map_hdu[0].data)
    exp_map_hdu[0].data[exp_map_hdu[0].data < ten_percent_limit] = 0
    time, fx, fy, photons = read_columns(combined_eventslist_hdu)
    fx_prime, fy_prime = instrument_to_astronomical(fx, fy, PC1_1, PC1_2, PC2_1, PC2_2)

    warnings.simplefilter("ignore", category=AstropyWarning)
    w = WCS(wcs_hdu[0].header)
    sky = w.pixel_to_world(
        combined_eventslist_hdu[1].data["Fx"], combined_eventslist_hdu[1].data["Fy"]
    )

    combined_eventslist_hdu[1].data["Sky_RA"] = sky.ra.degree
    combined_eventslist_hdu[1].data["Sky_DEC"] = sky.dec.degree
    (
        combined_eventslist_hdu[1].data["Fx_astronomical"],
        combined_eventslist_hdu[1].data["Fy_astronomical"],
    ) = instrument_to_astronomical(
        combined_eventslist_hdu[1].data["Fx"],
        combined_eventslist_hdu[1].data["Fy"],
        PC1_1,
        PC1_2,
        PC2_1,
        PC2_2,
    )
    combined_eventslist_hdu.close()

    bins_range = (-0.5, 4800.5, 1)
    bins = np.arange(*bins_range)
    events, _, _ = np.histogram2d(
        fy_prime, fx_prime, bins=(bins, bins), weights=photons
    )
    counts, _, _ = np.histogram2d(fy_prime, fx_prime, bins=(bins, bins))

    zero_weights_number = np.sum(photons == 0)
    if zero_weights_number != 0:
        print(f"\nWarning! {zero_weights_number} event weights have value zero.")

    CPS = events / exp_map_hdu[0].data
    CPS_error = CPS / np.sqrt(counts)

    # To remove Nan values from CPS and CPS_error maps using exposure map.
    CPS[exp_map_hdu[0].data == 0] = 0
    CPS_error[exp_map_hdu[0].data == 0] = 0

    # To remove Nan values due to the np.sqrt(counts) denominator.
    CPS_error[counts == 0] = 0

    CPS = CPS.astype(np.single)
    CPS_error = CPS_error.astype(np.single)

    wcs_hdu[0].header["PC1_1"] = 1
    wcs_hdu[0].header["PC1_2"] = 0
    wcs_hdu[0].header["PC2_1"] = 0
    wcs_hdu[0].header["PC2_2"] = 1
    wcs_hdu[0].header["CDELT1"] = CDELT1
    wcs_hdu[0].header["CDELT2"] = CDELT2
    if channel == "NUV":
        wcs_hdu[0].header["CRPIX1"] = 2401.5
        wcs_hdu[0].header["CRPIX2"] = 2401.5

    CPS_hdu = fits.PrimaryHDU(CPS)
    CPS_hdu.header.update(wcs_hdu[0].header)
    CPS_hdu.header["RA_PNT"] = RA_pointing
    CPS_hdu.header["DEC_PNT"] = DEC_pointing
    CPS_hdu.header["AVGFRMRT"] = framecount_per_sec
    CPS_hdu.header["STDFRMRT"] = framecount_per_sec_stdev
    CPS_hdu.header["ASTMETH"] = method
    CPS_filename = f"RA_Dec_combined_{channel}_{window_filter}_CPS.fits"
    CPS_hdu.writeto(CPS_filename, overwrite=True)
    print(f"Combined CPS image: {CPS_filename}")

    CPS_error_hdu = fits.PrimaryHDU(CPS_error)
    CPS_error_hdu.header.update(wcs_hdu[0].header)
    CPS_error_hdu.header["RA_PNT"] = RA_pointing
    CPS_error_hdu.header["DEC_PNT"] = DEC_pointing
    CPS_error_hdu.header["AVGFRMRT"] = framecount_per_sec
    CPS_error_hdu.header["STDFRMRT"] = framecount_per_sec_stdev
    CPS_error_hdu.header["ASTMETH"] = method
    CPS_error_filename = f"RA_Dec_combined_{channel}_{window_filter}_CPS_error.fits"
    CPS_error_hdu.writeto(CPS_error_filename, overwrite=True)
    print(f"Combined CPS error image: {CPS_error_filename}")


def make_combined_final_products():
    filelist = glob("[F,N]UV*_transformations.txt")
    if len(filelist) > 0:
        for file in filelist:
            channel = file[:3]
            window_filter = file[4:10]
            print(f"\nWorking on {channel} {window_filter}:")
            combined_eventslist_filename = (
                f"combined_{channel}_{window_filter}_events_list.fits"
            )
            exp_map_filename = f"RA_Dec_combined_{channel}_{window_filter}_exp_map.fits"
            wcs_filename = (
                f"combined_{channel}_{window_filter}_CPS_coo_fixed_platescale.wcs"
            )
            if Path(wcs_filename).is_file():
                method = "Direct"
            else:
                print(
                    f"No UV astrometry, using VIS astrometry for {channel} {window_filter}."
                )
                method = "Indirect"
                wcs_filename = f"UV_oriented_combined_VIS_for_{channel}_{window_filter}_coo_fixed_platescale.wcs"
            combined_eventslist_hdu = fits.open(
                combined_eventslist_filename, mode="update"
            )
            exp_map_hdu = fits.open(exp_map_filename)
            create_exp_map_in_seconds(
                exp_map_hdu, combined_eventslist_hdu, channel, window_filter, method
            )
            create_CPS_and_error_images(
                exp_map_hdu,
                combined_eventslist_hdu,
                channel,
                window_filter,
                wcs_filename,
                method,
            )
    print("\nProcessing completed")


np.seterr(divide="ignore", invalid="ignore")
make_combined_final_products()
