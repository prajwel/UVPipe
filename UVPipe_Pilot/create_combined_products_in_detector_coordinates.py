#!/usr/bin/env python

import numpy as np

from glob import glob
from astropy.io import fits


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


def create_exp_map_in_seconds(
    exp_map_hdu, combined_eventslist_hdu, channel, window_filter
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
    hdu = fits.PrimaryHDU(data)
    hdu.header["RA_PNT"] = RA_pointing
    hdu.header["DEC_PNT"] = DEC_pointing
    hdu.header["AVGFRMRT"] = framecount_per_sec
    hdu.header["STDFRMRT"] = framecount_per_sec_stdev
    expmap_in_sec_filename = (
        f"combined_{channel}_{window_filter}_exp_map_in_seconds.fits"
    )
    hdu.writeto(expmap_in_sec_filename, overwrite=True)
    print(f"\nCombined exposure map in seconds: {expmap_in_sec_filename}")


def create_CPS_and_error_images(
    exp_map_hdu, combined_eventslist_hdu, channel, window_filter
):
    (
        framecount_per_sec,
        framecount_per_sec_stdev,
        RA_pointing,
        DEC_pointing,
    ) = get_header_values(combined_eventslist_hdu)
    ten_percent_limit = 0.1 * np.max(exp_map_hdu[0].data)
    exp_map_hdu[0].data[exp_map_hdu[0].data < ten_percent_limit] = 0
    time, fx, fy, photons = read_columns(combined_eventslist_hdu)
    bins_range = (-0.5, 4800.5, 1)
    bins = np.arange(*bins_range)
    events, _, _ = np.histogram2d(fy, fx, bins=(bins, bins), weights=photons)
    counts, _, _ = np.histogram2d(fy, fx, bins=(bins, bins))

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

    CPS_hdu = fits.PrimaryHDU(CPS)
    CPS_hdu.header["RA_PNT"] = RA_pointing
    CPS_hdu.header["DEC_PNT"] = DEC_pointing
    CPS_hdu.header["AVGFRMRT"] = framecount_per_sec
    CPS_hdu.header["STDFRMRT"] = framecount_per_sec_stdev
    CPS_filename = f"combined_{channel}_{window_filter}_CPS.fits"
    CPS_hdu.writeto(CPS_filename, overwrite=True)
    print(f"Combined CPS image: {CPS_filename}")

    CPS_error_hdu = fits.PrimaryHDU(CPS_error)
    CPS_error_hdu.header["RA_PNT"] = RA_pointing
    CPS_error_hdu.header["DEC_PNT"] = DEC_pointing
    CPS_error_hdu.header["AVGFRMRT"] = framecount_per_sec
    CPS_error_hdu.header["STDFRMRT"] = framecount_per_sec_stdev
    CPS_error_filename = f"combined_{channel}_{window_filter}_CPS_error.fits"
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
            exp_map_filename = f"combined_{channel}_{window_filter}_exp_map.fits"
            combined_eventslist_hdu = fits.open(combined_eventslist_filename)
            exp_map_hdu = fits.open(exp_map_filename)
            create_exp_map_in_seconds(
                exp_map_hdu, combined_eventslist_hdu, channel, window_filter
            )
            create_CPS_and_error_images(
                exp_map_hdu, combined_eventslist_hdu, channel, window_filter
            )
    print("\nProcessing completed")


np.seterr(divide="ignore", invalid="ignore")
make_combined_final_products()
