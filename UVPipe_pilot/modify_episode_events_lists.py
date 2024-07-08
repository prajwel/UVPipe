#!/usr/bin/env python3

import numpy as np

from glob import glob
from pathlib import Path
from astropy.io import fits


def get_plate_scales(channel, filter_id):
    XP_FUV = 0.4168705
    YP_FUV = 0.416848964

    XP_NUV = 0.417310
    YP_NUV = 0.417218

    XP_NUV_F2 = 0.419686
    YP_NUV_F2 = 0.419184

    XP_NUV_F1 = 0.4173677
    YP_NUV_F1 = 0.417299

    if channel == "FUV":
        return XP_FUV, YP_FUV
    elif channel == "NUV":
        if filter_id[-2:] == "F2":
            return XP_NUV_F2, YP_NUV_F2
        elif filter_id[-2:] == "F1":
            return XP_NUV_F1, YP_NUV_F1
        else:
            return XP_NUV, YP_NUV
    else:
        raise ValueError


def modify_events_list(events_list_filename):
    channel_filterid = Path(events_list_filename).name[24:36]
    filter_id = channel_filterid[-6:]

    if channel_filterid[:4] == "uvtF":
        channel = "FUV"
    elif channel_filterid[:4] == "uvtN":
        channel = "NUV"
    else:
        raise ValueError

    print(f"\nModifying plate scale of {events_list_filename}")
    print(channel, filter_id)

    hdu = fits.open(events_list_filename, mode="update")

    if "Fx_original" in hdu[1].data.columns.names:
        print("No changes, event list has already been modified.")
        return

    X_platescale, Y_platescale = get_plate_scales(channel, filter_id)
    average_platescale = (X_platescale + Y_platescale) / 2

    print(f"Plate scale along X = {X_platescale}")
    print(f"Plate scale along Y = {Y_platescale}")

    modified_Fx = fits.Column(name="Fx_original", format="D", array=hdu[1].data["Fx"])
    modified_Fy = fits.Column(name="Fy_original", format="D", array=hdu[1].data["Fy"])

    new_cols = fits.ColDefs([modified_Fx, modified_Fy])
    new_table = fits.FITS_rec.from_columns(hdu[1].data.columns + new_cols)
    hdu[1].data = new_table

    hdu[1].data["Fx"] = (
        hdu[1].data["Fx"] - 2400
    ) * X_platescale / average_platescale + 2400
    hdu[1].data["Fy"] = (
        hdu[1].data["Fy"] - 2400
    ) * Y_platescale / average_platescale + 2400

    hdu.flush()


def revert_modified_events_list(events_list_filename):
    print(f"\nReverting {events_list_filename}")

    hdu = fits.open(events_list_filename, mode="update")

    if "Fx_original" not in hdu[1].data.columns.names:
        print("No changes done, event list is unmodified or already reverted.")
        return

    hdu[1].data["Fx"] = hdu[1].data["Fx_original"]
    hdu[1].data["Fy"] = hdu[1].data["Fy_original"]

    hdu[1].data.columns.del_col("Fx_original")
    hdu[1].data.columns.del_col("Fy_original")

    hdu.flush()


def reject_frames_from_events_list(events_list_filename):
    print(f"\nRejecting frames of {events_list_filename}")

    hdu = fits.open(events_list_filename, mode="update")

    if "FRMRJCTD" in hdu[0].header:
        print("No changes, event list has already been modified.")
        return

    first_framecount = hdu[1].data["FrameCount"][0]
    frames_to_remove = np.array([0, 1]) + first_framecount
    frames_to_remove_mask = ~np.isin(hdu[1].data["FrameCount"], frames_to_remove)
    hdu[1].data = hdu[1].data[frames_to_remove_mask]
    hdu[0].header["FRMRJCTD"] = 2

    hdu.flush()


# events_list_files = glob("uvt_*/[N,F]*/*l2ce.fits")
# for events_list_file in events_list_files:
#     revert_modified_events_list(events_list_file)

events_list_files = glob("uvt_*/[N,F]*/*l2ce.fits")
for events_list_file in events_list_files:
    modify_events_list(events_list_file)
    reject_frames_from_events_list(events_list_file)
