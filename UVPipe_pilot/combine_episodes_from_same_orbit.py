#!/usr/bin/env python3

import os
import sys
import shutil
import numpy as np

from glob import glob
from astropy.io import fits
from pathlib import Path


def get_first_frametime(VIS_RAS_file):
    hdu = fits.open(VIS_RAS_file)
    reference_frame_time = hdu[2].data["Time"][0]
    return reference_frame_time


def combine_episodes_from_same_orbit(to_combine_paths, channel):
    INT_TIME_from_header = []
    i = 0
    for path in to_combine_paths:
        parent_path = Path(path).parents[1]
        img_path = str(parent_path) + f"/{channel[0]}*/*I_l2img*"

        if len(glob(img_path)) == 1:
            img_path = glob(img_path)[0]
        else:
            print(
                "Requires exactly one file matching *I_l2img* in directory {}!".format(
                    os.path.split(img_path)[0]
                )
            )
            sys.exit()

        eventslist_path = str(parent_path) + f"/{channel[0]}*/**l2ce.fits"
        eventslist_path = glob(eventslist_path)[0]

        exp_path = str(parent_path) + f"/{channel[0]}*/*I_l2exp*"
        exp_path = glob(exp_path)[0]

        img_hdu = fits.open(img_path)
        INT_TIME_from_header.append(img_hdu[0].header["INT_TIME"])

        if i == 0:
            first_eventslist_path = eventslist_path
            first_img_path = img_path
            first_exp_path = exp_path
        else:
            hdu_add = fits.open(eventslist_path)
            hdu_base = fits.open(first_eventslist_path)

            nrows_base = hdu_base[1].data.shape[0]
            nrows_add = hdu_add[1].data.shape[0]
            nrows = nrows_base + nrows_add

            hdu = fits.BinTableHDU.from_columns(
                hdu_base[1].columns, nrows=nrows, fill=True
            )
            for colname in hdu_base[1].columns.names:
                hdu.data[colname][:nrows_base] = hdu_base[1].data[colname]
                hdu.data[colname][nrows_base:] = hdu_add[1].data[colname]

            hduP = fits.PrimaryHDU()
            hduP.header = hdu_base[0].header

            hdu.name = hdu_base[1].name
            hdu_base = fits.HDUList([hduP, hdu])
            hdu_base.writeto(first_eventslist_path, overwrite=True)

            exp_hdu_add = fits.open(exp_path)
            exp_hdu_base = fits.open(first_exp_path)
            exp_hdu_base[0].header["EXP_TIME"] = (
                exp_hdu_base[0].header["EXP_TIME"] + exp_hdu_add[0].header["EXP_TIME"]
            )
            exp_hdu_base[0].data = exp_hdu_base[0].data + exp_hdu_add[0].data
            exp_hdu_base.writeto(first_exp_path, overwrite=True)
            dir_to_remove = os.path.dirname(exp_path)
            shutil.rmtree(dir_to_remove)
            os.mkdir(dir_to_remove)
        i = i + 1

    AVG_INT_TIME = np.mean(INT_TIME_from_header)
    first_img_hdu = fits.open(first_img_path, mode="update")
    first_img_hdu[0].header["INT_TIME"] = AVG_INT_TIME
    first_img_hdu.close()


def combine_episodes_in_channel(channel):
    all_events_lists = glob(f"uvt_*/{channel[0]}*/*l2ce.fits")
    if len(all_events_lists) > 0:
        window_filter_list = [file[42:48] for file in all_events_lists]
        window_filter_list = list(set(window_filter_list))
        for window_filter in window_filter_list:
            chosen_window_filter_events_lists = glob(
                f"uvt_*/{channel[0]}*/*{window_filter}*l2ce.fits"
            )
            VIS_RAS_files = [
                glob(os.path.split(d)[0].replace(f"{channel[0]}", "V") + "/*dr.fits")[0]
                for d in chosen_window_filter_events_lists
            ]
            reference_frame_times = [
                get_first_frametime(VIS_RAS) for VIS_RAS in VIS_RAS_files
            ]
            reference_frame_times = np.array(reference_frame_times)
            unique_frame_times, unique_counts = np.unique(
                reference_frame_times, return_counts=True
            )
            unique_frame_times = unique_frame_times[unique_counts > 1]
            for unique_frame_time in unique_frame_times:
                mask = reference_frame_times == unique_frame_time
                to_combine_paths = np.array(VIS_RAS_files)[mask]
                print(
                    (
                        f"\nCombining the {channel} {window_filter} files corresponding "
                        f"to the following episodes:\n{to_combine_paths}"
                    )
                )
                combine_episodes_from_same_orbit(to_combine_paths, channel)


def is_dir_filled(dir_name):
    if os.path.isdir(dir_name):
        if not os.listdir(dir_name):
            return False
        else:
            return True
    else:
        return False


def clear_VIS_dirs():
    VIS_RAS_files = glob("uvt_*/V_*/*dr.fits")
    for VIS_RAS_file in np.sort(VIS_RAS_files):
        VIS_path = os.path.dirname(VIS_RAS_file)
        NUV_dir = VIS_path.replace("V_", "N_")
        FUV_dir = VIS_path.replace("V_", "F_")
        filled_UV_dir = np.logical_or(is_dir_filled(NUV_dir), is_dir_filled(FUV_dir))
        if not filled_UV_dir:
            print(f"Removing contents of {VIS_path} due to empty UV directories.")
            shutil.rmtree(VIS_path)
            os.mkdir(VIS_path)


combine_episodes_in_channel("FUV")
combine_episodes_in_channel("NUV")
clear_VIS_dirs()
print("\nDone!")
