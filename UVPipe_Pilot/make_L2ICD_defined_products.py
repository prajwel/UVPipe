#!/usr/bin/env python3

import os
import shutil
import numpy as np

from glob import glob
from astropy.io import fits
from astropy.time import Time
from datetime import datetime


allowed_CALDB_versions = ["2.1.0", "2.1.0-frm05274"]

# Dictionary containing all channels
all_filters = {
    "VIS": {
        "F1": {"old name": "VIS3", "new name": "V461W"},
        "F2": {"old name": "VIS2", "new name": "V391M"},
        "F3": {"old name": "VIS1", "new name": "V347M"},
        "F4": {"old name": "ND1", "new name": "V435ND"},
        "F5": {"old name": "BK7", "new name": "V420W"},
    },
    "NUV": {
        "F1": {"old name": "Silica - 1", "new name": "N242W"},
        "F2": {"old name": "NUVB15", "new name": "N219M"},
        "F3": {"old name": "NUVB13", "new name": "N245M"},
        "F4": {"old name": "Grating", "new name": ""},
        "F5": {"old name": "NUVB4", "new name": "N263M"},
        "F6": {"old name": "NUVN2", "new name": "N279N"},
        "F7": {"old name": "Silica - 2", "new name": "N242Wa"},
    },
    "FUV": {
        "F1": {"old name": "CaF2 - 1", "new name": "F148W"},
        "F2": {"old name": "BaF2", "new name": "F154W"},
        "F3": {"old name": "Sapphire", "new name": "F169M"},
        "F4": {"old name": "Grating - 1", "new name": ""},
        "F5": {"old name": "Silica", "new name": "F172M"},
        "F6": {"old name": "Grating - 2", "new name": ""},
        "F7": {"old name": "CaF2 - 2", "new name": "F148Wa"},
    },
}

imcoosys_dict = {"I": "Instrument", "A": "Astronomical"}

datatype_dict = {
    "img": "count-rate image",
    "err": "count-rate error image",
    "exp": "exposure map",
}


# To rename curvit generated files.
def new_lightcurve_name(filename):
    _, _, xp, yp, _, channel, window_filter, _, _ = filename.split("_")
    new_name = f"{channel}_{window_filter}_X_{xp}_Y_{yp}_lightcurve.png"
    return new_name


def new_fullfield_name(filename):
    _, xp, yp, _, channel, window_filter, _, _ = filename.split("_")
    new_name = f"{channel}_{window_filter}_X_{xp}_Y_{yp}_fullfield.png"
    return new_name


def new_source_zoomed_name(filename):
    _, _, xp, yp, _, channel, window_filter, _, _ = filename.split("_")
    new_name = f"{channel}_{window_filter}_X_{xp}_Y_{yp}_source.png"
    return new_name


# To rename astrometry error vector plots.
def new_error_vector_name(filename):
    _, channel, window_filter, _, _, _, _ = filename.split("_")
    new_name = f"{channel}_{window_filter}_astrometry_error_vectors.png"
    return new_name


# To keep only the header of a FITS file.
def strip_data(episode_file, write_location):
    hdu = fits.open(episode_file)
    skeleton_hdu = fits.PrimaryHDU(header=hdu[0].header)
    skeleton_hdu.writeto(write_location, overwrite=True)


# To change mission elapsed time in seconds to modified julian date.
def met_to_mjd(met):
    jan2010 = 55197.0  # 2010.0(UTC) expressed with MJD format and scale UTC.
    mjd = (met / 86400.0) + jan2010  # 1 julian day = 86400 seconds.
    return mjd


# Check if the input list contains the same values or not.
def check_for_consistency(value_list):
    if len(set(value_list)) == 1:
        value = value_list[0]
        return value
    else:
        raise RuntimeError(f"Please check {value_list}")


def update_header(
    filename,
    combined_events_list,
    combined_exp_map_in_frames,
    combined_exp_map_in_seconds,
    channel,
    window_filter,
    level1_dataset,
    CALDB_version,
    SEER_1,
    SEER_2,
    header_template_file="header_template.txt",
):
    inputfile_hdu = fits.open(filename, mode="update")

    events_hdu = fits.open(combined_events_list)
    exp_start_time = met_to_mjd(min(events_hdu[1].data["MJD_L2"]))
    exp_end_time = met_to_mjd(max(events_hdu[1].data["MJD_L2"]))
    start_time = Time(exp_start_time, format="mjd")
    end_time = Time(exp_end_time, format="mjd")
    AVGFRMRT = events_hdu[0].header["AVGFRMRT"]
    AVGFRMRT = np.round(AVGFRMRT, 6)
    STDFRMRT = events_hdu[0].header["STDFRMRT"]
    STDFRMRT = f"{STDFRMRT:.4E}"

    transformations_filename = f"{channel}_{window_filter}_transformations.txt"
    transformations = np.genfromtxt(transformations_filename, dtype=None, encoding=None)
    transformations = np.atleast_1d(transformations)
    N_EPISDS = len(transformations)

    exp_map_in_frames_hdu = fits.open(combined_exp_map_in_frames)
    exp_map_in_seconds_hdu = fits.open(combined_exp_map_in_seconds)
    MEDFRAME = np.median(
        exp_map_in_frames_hdu[0].data[exp_map_in_frames_hdu[0].data > 0]
    )
    EXP_TIME = np.median(
        exp_map_in_seconds_hdu[0].data[exp_map_in_seconds_hdu[0].data > 0]
    )
    RA_PNT = exp_map_in_seconds_hdu[0].header["CRVAL1"]
    DEC_PNT = exp_map_in_seconds_hdu[0].header["CRVAL2"]

    chosen_window_filter_image_paths = glob(
        f"uvt_*/{channel[0]}*/*{window_filter}*A_l2img.fits*"
    )
    OBS_MODE = []
    WIN_XOFF = []
    WIN_YOFF = []
    WIN_X_SZ = []
    WIN_Y_SZ = []
    CENTROID = []
    OBSERVER = []
    OBJECT = []

    for path in chosen_window_filter_image_paths:
        image_hdu = fits.open(path)
        OBS_MODE.append(image_hdu[0].header["OBS_MODE"])
        WIN_XOFF.append(image_hdu[0].header["WIN_XOFF"])
        WIN_YOFF.append(image_hdu[0].header["WIN_YOFF"])
        WIN_X_SZ.append(image_hdu[0].header["WIN_X_SZ"])
        WIN_Y_SZ.append(image_hdu[0].header["WIN_Y_SZ"])
        CENTROID.append(image_hdu[0].header["CENTROID"])
        OBSERVER.append(image_hdu[0].header["OBSERVER"])
        OBJECT.append(image_hdu[0].header["OBJECT"])

    OBS_MODE = check_for_consistency(OBS_MODE)
    WIN_XOFF = check_for_consistency(WIN_XOFF)
    WIN_YOFF = check_for_consistency(WIN_YOFF)
    WIN_X_SZ = check_for_consistency(WIN_X_SZ)
    WIN_Y_SZ = check_for_consistency(WIN_Y_SZ)
    CENTROID = check_for_consistency(CENTROID)
    OBSERVER = check_for_consistency(OBSERVER)
    OBJECT = check_for_consistency(OBJECT)

    header_template = fits.Header.fromtextfile(header_template_file)

    # Basic information
    header_template["FILENAME"] = filename
    header_template["FILEDATE"] = now_time
    header_template["CALDBVER"] = CALDB_version
    header_template["IMCOOSYS"] = imcoosys_dict[header_template["FILENAME"][36]]
    header_template["DATATYPE"] = datatype_dict[header_template["FILENAME"][40:43]]

    # Instrument configuration information
    filterid = window_filter[-2:]
    header_template["DETECTOR"] = channel
    header_template["FILTERID"] = filterid
    header_template["FILNAMEN"] = all_filters[channel][filterid]["new name"]
    header_template["FILNAMEO"] = all_filters[channel][filterid]["old name"]
    header_template["OBS_MODE"] = OBS_MODE
    header_template["WIN_XOFF"] = WIN_XOFF
    header_template["WIN_YOFF"] = WIN_YOFF
    header_template["WIN_X_SZ"] = WIN_X_SZ
    header_template["WIN_Y_SZ"] = WIN_Y_SZ
    header_template["CENTROID"] = CENTROID
    header_template["RA_PNT"] = RA_PNT
    header_template["DEC_PNT"] = DEC_PNT

    # Observation identifiers
    header_template["DATE-OBS"] = start_time.isot.split("T")[0]
    header_template["TIME-OBS"] = start_time.isot.split("T")[1]
    header_template["DATE-BEG"] = start_time.isot
    header_template["DATE-END"] = end_time.isot
    header_template["OBS_ID"] = observation_ID
    header_template["PROPOSID"] = proposal_ID
    header_template["TARGETID"] = target_ID
    header_template["OBSERVER"] = OBSERVER
    header_template["OBSLABEL"] = OBJECT

    # Exposure parameters
    header_template["EXP_TIME"] = EXP_TIME
    header_template["MEDFRAME"] = MEDFRAME
    header_template["AVGFRMRT"] = AVGFRMRT
    header_template["STDFRMRT"] = STDFRMRT
    header_template["EXPSTART"] = exp_start_time
    header_template["EXPEND"] = exp_end_time

    # WCS information
    header_template["ASTMETH"] = inputfile_hdu[0].header["ASTMETH"]
    header_template["WCSAXES"] = inputfile_hdu[0].header["WCSAXES"]
    header_template["CRPIX1"] = inputfile_hdu[0].header["CRPIX1"]
    header_template["CRPIX2"] = inputfile_hdu[0].header["CRPIX2"]
    header_template["CDELT1"] = inputfile_hdu[0].header["CDELT1"]
    header_template["CDELT2"] = inputfile_hdu[0].header["CDELT2"]
    header_template["CRVAL1"] = inputfile_hdu[0].header["CRVAL1"]
    header_template["CRVAL2"] = inputfile_hdu[0].header["CRVAL2"]
    header_template["LATPOLE"] = inputfile_hdu[0].header["LATPOLE"]
    header_template["MJDREF"] = inputfile_hdu[0].header["MJDREF"]
    header_template["PC1_1"] = inputfile_hdu[0].header["PC1_1"]
    header_template["PC1_2"] = inputfile_hdu[0].header["PC1_2"]
    header_template["PC2_1"] = inputfile_hdu[0].header["PC2_1"]
    header_template["PC2_2"] = inputfile_hdu[0].header["PC2_2"]
    header_template["CUNIT1"] = inputfile_hdu[0].header["CUNIT1"]
    header_template["CUNIT2"] = inputfile_hdu[0].header["CUNIT2"]
    header_template["CTYPE1"] = inputfile_hdu[0].header["CTYPE1"]
    header_template["CTYPE2"] = inputfile_hdu[0].header["CTYPE2"]
    header_template["LONPOLE"] = inputfile_hdu[0].header["LONPOLE"]
    header_template["RADESYS"] = inputfile_hdu[0].header["RADESYS"]

    # Level1 to Level2 processing
    header_template["LVL1FILE"] = level1_dataset

    # Combining information
    header_template["N_EPISDS"] = N_EPISDS
    header_template["COMBMETH"] = events_hdu[0].header["COMBMETH"]

    # Image quality information
    header_template["SEER_1"] = SEER_1
    header_template["SEER_2"] = SEER_2

    inputfile_hdu[0].header = header_template
    inputfile_hdu.close()


def update_events_list_header(
    filename,
    combined_exp_map_in_frames,
    combined_exp_map_in_seconds,
    channel,
    window_filter,
    level1_dataset,
    CALDB_version,
    SEER_1,
    SEER_2,
    header_template_file="header_template_events_list.txt",
):
    inputfile_hdu = fits.open(filename, mode="update")

    exp_start_time = met_to_mjd(min(inputfile_hdu[1].data["MJD_L2"]))
    exp_end_time = met_to_mjd(max(inputfile_hdu[1].data["MJD_L2"]))
    start_time = Time(exp_start_time, format="mjd")
    end_time = Time(exp_end_time, format="mjd")
    AVGFRMRT = inputfile_hdu[0].header["AVGFRMRT"]
    AVGFRMRT = np.round(AVGFRMRT, 6)
    STDFRMRT = inputfile_hdu[0].header["STDFRMRT"]
    STDFRMRT = f"{STDFRMRT:.4E}"

    transformations_filename = f"{channel}_{window_filter}_transformations.txt"
    transformations = np.genfromtxt(transformations_filename, dtype=None, encoding=None)
    transformations = np.atleast_1d(transformations)
    N_EPISDS = len(transformations)

    exp_map_in_frames_hdu = fits.open(combined_exp_map_in_frames)
    exp_map_in_seconds_hdu = fits.open(combined_exp_map_in_seconds)
    MEDFRAME = np.median(
        exp_map_in_frames_hdu[0].data[exp_map_in_frames_hdu[0].data > 0]
    )
    EXP_TIME = np.median(
        exp_map_in_seconds_hdu[0].data[exp_map_in_seconds_hdu[0].data > 0]
    )
    RA_PNT = exp_map_in_seconds_hdu[0].header["CRVAL1"]
    DEC_PNT = exp_map_in_seconds_hdu[0].header["CRVAL2"]

    chosen_window_filter_image_paths = glob(
        f"uvt_*/{channel[0]}*/*{window_filter}*A_l2img.fits*"
    )
    OBS_MODE = []
    WIN_XOFF = []
    WIN_YOFF = []
    WIN_X_SZ = []
    WIN_Y_SZ = []
    CENTROID = []
    OBSERVER = []
    OBJECT = []

    for path in chosen_window_filter_image_paths:
        image_hdu = fits.open(path)
        OBS_MODE.append(image_hdu[0].header["OBS_MODE"])
        WIN_XOFF.append(image_hdu[0].header["WIN_XOFF"])
        WIN_YOFF.append(image_hdu[0].header["WIN_YOFF"])
        WIN_X_SZ.append(image_hdu[0].header["WIN_X_SZ"])
        WIN_Y_SZ.append(image_hdu[0].header["WIN_Y_SZ"])
        CENTROID.append(image_hdu[0].header["CENTROID"])
        OBSERVER.append(image_hdu[0].header["OBSERVER"])
        OBJECT.append(image_hdu[0].header["OBJECT"])

    OBS_MODE = check_for_consistency(OBS_MODE)
    WIN_XOFF = check_for_consistency(WIN_XOFF)
    WIN_YOFF = check_for_consistency(WIN_YOFF)
    WIN_X_SZ = check_for_consistency(WIN_X_SZ)
    WIN_Y_SZ = check_for_consistency(WIN_Y_SZ)
    CENTROID = check_for_consistency(CENTROID)
    OBSERVER = check_for_consistency(OBSERVER)
    OBJECT = check_for_consistency(OBJECT)

    header_template = fits.Header.fromtextfile(header_template_file)

    # Basic information
    header_template["FILENAME"] = filename
    header_template["FILEDATE"] = now_time
    header_template["CALDBVER"] = CALDB_version

    # Instrument configuration information
    filterid = window_filter[-2:]
    header_template["DETECTOR"] = channel
    header_template["FILTERID"] = filterid
    header_template["FILNAMEN"] = all_filters[channel][filterid]["new name"]
    header_template["FILNAMEO"] = all_filters[channel][filterid]["old name"]
    header_template["OBS_MODE"] = OBS_MODE
    header_template["WIN_XOFF"] = WIN_XOFF
    header_template["WIN_YOFF"] = WIN_YOFF
    header_template["WIN_X_SZ"] = WIN_X_SZ
    header_template["WIN_Y_SZ"] = WIN_Y_SZ
    header_template["CENTROID"] = CENTROID
    header_template["RA_PNT"] = RA_PNT
    header_template["DEC_PNT"] = DEC_PNT

    # Observation identifiers
    header_template["DATE-OBS"] = start_time.isot.split("T")[0]
    header_template["TIME-OBS"] = start_time.isot.split("T")[1]
    header_template["DATE-BEG"] = start_time.isot
    header_template["DATE-END"] = end_time.isot
    header_template["OBS_ID"] = observation_ID
    header_template["PROPOSID"] = proposal_ID
    header_template["TARGETID"] = target_ID
    header_template["OBSERVER"] = OBSERVER
    header_template["OBSLABEL"] = OBJECT

    # Exposure parameters
    header_template["EXP_TIME"] = EXP_TIME
    header_template["MEDFRAME"] = MEDFRAME
    header_template["AVGFRMRT"] = AVGFRMRT
    header_template["STDFRMRT"] = STDFRMRT
    header_template["EXPSTART"] = exp_start_time
    header_template["EXPEND"] = exp_end_time

    # Level1 to Level2 processing
    header_template["LVL1FILE"] = level1_dataset

    # Combining information
    header_template["N_EPISDS"] = N_EPISDS
    header_template["COMBMETH"] = inputfile_hdu[0].header["COMBMETH"]

    # Image quality information
    header_template["SEER_1"] = SEER_1
    header_template["SEER_2"] = SEER_2

    inputfile_hdu[0].header = header_template
    inputfile_hdu.close()


def generate_ICD_UV_products(channel, window_filter, CALDB_version):
    combined_events_list = f"combined_{channel}_{window_filter}_events_list.fits"

    combined_CPS_map = f"combined_{channel}_{window_filter}_CPS.fits"
    RA_Dec_combined_CPS_map = f"RA_Dec_combined_{channel}_{window_filter}_CPS.fits"

    combined_CPS_error_map = f"combined_{channel}_{window_filter}_CPS_error.fits"
    RA_Dec_combined_CPS_error_map = (
        f"RA_Dec_combined_{channel}_{window_filter}_CPS_error.fits"
    )

    combined_exp_map_in_frames = f"combined_{channel}_{window_filter}_exp_map.fits"
    combined_exp_map_in_seconds = (
        f"combined_{channel}_{window_filter}_exp_map_in_seconds.fits"
    )
    RA_Dec_combined_exp_map_in_seconds = (
        f"RA_Dec_combined_{channel}_{window_filter}_exp_map_in_seconds.fits"
    )

    ICD_combined_events_list = (
        f"AS1{observation_ID}uvt{channel[0]}II{window_filter}_l2ce.fits"
    )

    ICD_combined_CPS_map = (
        f"AS1{observation_ID}uvt{channel[0]}II{window_filter}I_l2img.fits"
    )
    ICD_RA_Dec_combined_CPS_map = (
        f"AS1{observation_ID}uvt{channel[0]}II{window_filter}A_l2img.fits"
    )

    ICD_combined_CPS_error_map = (
        f"AS1{observation_ID}uvt{channel[0]}II{window_filter}I_l2err.fits"
    )
    ICD_RA_Dec_combined_CPS_error_map = (
        f"AS1{observation_ID}uvt{channel[0]}II{window_filter}A_l2err.fits"
    )

    ICD_combined_exp_map_in_seconds = (
        f"AS1{observation_ID}uvt{channel[0]}II{window_filter}I_l2exp.fits"
    )
    ICD_RA_Dec_combined_exp_map_in_seconds = (
        f"AS1{observation_ID}uvt{channel[0]}II{window_filter}A_l2exp.fits"
    )

    shutil.copy(combined_events_list, ICD_combined_events_list)

    shutil.copy(combined_CPS_map, ICD_combined_CPS_map)
    shutil.copy(RA_Dec_combined_CPS_map, ICD_RA_Dec_combined_CPS_map)

    shutil.copy(combined_CPS_error_map, ICD_combined_CPS_error_map)
    shutil.copy(RA_Dec_combined_CPS_error_map, ICD_RA_Dec_combined_CPS_error_map)

    shutil.copy(combined_exp_map_in_seconds, ICD_combined_exp_map_in_seconds)
    shutil.copy(
        RA_Dec_combined_exp_map_in_seconds, ICD_RA_Dec_combined_exp_map_in_seconds
    )

    CPS_hdu = fits.open(combined_CPS_map)
    SEER_1 = CPS_hdu[0].header["SEER_1"]
    SEER_2 = CPS_hdu[0].header["SEER_2"]

    common_args = dict(
        combined_events_list=combined_events_list,
        combined_exp_map_in_frames=combined_exp_map_in_frames,
        combined_exp_map_in_seconds=combined_exp_map_in_seconds,
        channel=channel,
        window_filter=window_filter,
        level1_dataset=level1_dataset,
        CALDB_version=CALDB_version,
        SEER_1=SEER_1,
        SEER_2=SEER_2,
    )

    update_header(ICD_combined_CPS_map, **common_args)
    update_header(ICD_RA_Dec_combined_CPS_map, **common_args)
    update_header(ICD_combined_CPS_error_map, **common_args)
    update_header(ICD_RA_Dec_combined_CPS_error_map, **common_args)
    update_header(ICD_combined_exp_map_in_seconds, **common_args)
    update_header(ICD_RA_Dec_combined_exp_map_in_seconds, **common_args)
    update_events_list_header(
        ICD_combined_events_list,
        combined_exp_map_in_frames,
        combined_exp_map_in_seconds,
        channel,
        window_filter,
        level1_dataset,
        CALDB_version,
        SEER_1,
        SEER_2,
    )


def update_VIS_header(
    filename,
    channel,
    window_filter,
    level1_dataset,
    CALDB_version,
    header_template_file="header_template_VIS.txt",
):
    inputfile_hdu = fits.open(filename, mode="update")

    exp_start_time = inputfile_hdu[0].header["EXPSTART"]
    exp_end_time = inputfile_hdu[0].header["EXPEND"]
    start_time = Time(exp_start_time, format="mjd")
    end_time = Time(exp_end_time, format="mjd")
    RA_PNT = inputfile_hdu[0].header["CRVAL1"]
    DEC_PNT = inputfile_hdu[0].header["CRVAL2"]

    chosen_window_filter_image_paths = glob(
        f"uvt_*/{channel[0]}*/*{window_filter}*_l2dr.fits*"
    )
    OBSERVER = []
    OBJECT = []

    for path in chosen_window_filter_image_paths:
        image_hdu = fits.open(path)
        OBSERVER.append(image_hdu[0].header["OBSERVER"])
        OBJECT.append(image_hdu[0].header["OBJECT"])

    OBSERVER = check_for_consistency(OBSERVER)
    OBJECT = check_for_consistency(OBJECT)

    header_template = fits.Header.fromtextfile(header_template_file)

    # Basic information
    header_template["FILENAME"] = filename
    header_template["FILEDATE"] = now_time
    header_template["CALDBVER"] = CALDB_version
    header_template["IMCOOSYS"] = imcoosys_dict[header_template["FILENAME"][36]]

    # Instrument configuration information
    filterid = window_filter[-2:]
    header_template["DETECTOR"] = channel
    header_template["FILTERID"] = filterid
    header_template["FILNAMEN"] = all_filters[channel][filterid]["new name"]
    header_template["FILNAMEO"] = all_filters[channel][filterid]["old name"]
    header_template["RA_PNT"] = RA_PNT
    header_template["DEC_PNT"] = DEC_PNT

    # Observation identifiers
    header_template["DATE-OBS"] = start_time.isot.split("T")[0]
    header_template["TIME-OBS"] = start_time.isot.split("T")[1]
    header_template["DATE-BEG"] = start_time.isot
    header_template["DATE-END"] = end_time.isot
    header_template["OBS_ID"] = observation_ID
    header_template["PROPOSID"] = proposal_ID
    header_template["TARGETID"] = target_ID
    header_template["OBSERVER"] = OBSERVER
    header_template["OBSLABEL"] = OBJECT

    # Exposure parameters
    header_template["EXP_TIME"] = inputfile_hdu[0].header["EXP_TIME"]
    header_template["EXPSTART"] = exp_start_time
    header_template["EXPEND"] = exp_end_time
    header_template["N_EPISDS"] = inputfile_hdu[0].header["N_EPISDS"]

    # WCS information
    header_template["WCSAXES"] = inputfile_hdu[0].header["WCSAXES"]
    header_template["CRPIX1"] = inputfile_hdu[0].header["CRPIX1"]
    header_template["CRPIX2"] = inputfile_hdu[0].header["CRPIX2"]
    header_template["CDELT1"] = inputfile_hdu[0].header["CDELT1"]
    header_template["CDELT2"] = inputfile_hdu[0].header["CDELT2"]
    header_template["CRVAL1"] = inputfile_hdu[0].header["CRVAL1"]
    header_template["CRVAL2"] = inputfile_hdu[0].header["CRVAL2"]
    header_template["LATPOLE"] = inputfile_hdu[0].header["LATPOLE"]
    header_template["MJDREF"] = inputfile_hdu[0].header["MJDREF"]

    if filename[36] == "I":
        header_template["PC1_1"] = inputfile_hdu[0].header["PC1_1"]
        header_template["PC1_2"] = inputfile_hdu[0].header["PC1_2"]
        header_template["PC2_1"] = inputfile_hdu[0].header["PC2_1"]
        header_template["PC2_2"] = inputfile_hdu[0].header["PC2_2"]
    else:
        header_template["PC1_1"] = 1
        header_template["PC1_2"] = 0
        header_template["PC2_1"] = 0
        header_template["PC2_2"] = 1

    header_template["CUNIT1"] = inputfile_hdu[0].header["CUNIT1"]
    header_template["CUNIT2"] = inputfile_hdu[0].header["CUNIT2"]
    header_template["CTYPE1"] = inputfile_hdu[0].header["CTYPE1"]
    header_template["CTYPE2"] = inputfile_hdu[0].header["CTYPE2"]
    header_template["LONPOLE"] = inputfile_hdu[0].header["LONPOLE"]
    header_template["RADESYS"] = inputfile_hdu[0].header["RADESYS"]

    # Level1 to Level2 processing
    header_template["LVL1FILE"] = level1_dataset

    inputfile_hdu[0].header = header_template
    inputfile_hdu.close()


def generate_ICD_VIS_products(channel, window_filter, CALDB_version):
    combined_VIS_image = f"combined_{channel}_{window_filter}.fits"
    RA_Dec_combined_VIS_image = f"RA_Dec_combined_{channel}_{window_filter}.fits"

    ICD_combined_VIS_image = (
        f"AS1{observation_ID}uvt{channel[0]}II{window_filter}I_l2ql.fits"
    )
    ICD_RA_Dec_combined_VIS_image = (
        f"AS1{observation_ID}uvt{channel[0]}II{window_filter}A_l2ql.fits"
    )

    shutil.copy(combined_VIS_image, ICD_combined_VIS_image)
    shutil.copy(RA_Dec_combined_VIS_image, ICD_RA_Dec_combined_VIS_image)

    common_args = dict(
        channel=channel,
        window_filter=window_filter,
        level1_dataset=level1_dataset,
        CALDB_version=CALDB_version,
    )

    update_VIS_header(ICD_combined_VIS_image, **common_args)
    update_VIS_header(ICD_RA_Dec_combined_VIS_image, **common_args)


def get_caldb_version(
    driver_module_param_file, allowed_CALDB_versions=allowed_CALDB_versions
):
    with open(driver_module_param_file, "r") as file:
        lines = file.readlines()
        CALDB_line = lines[8]
        CALDB_suffix = CALDB_line.split("/")[-2]
        CALDB_version = CALDB_suffix.split("CALDB_L2_v")[1]
        if CALDB_version not in allowed_CALDB_versions:
            raise ValueError(f"Incorrect CALDB version = {CALDB_version}")
        else:
            return CALDB_version


def make_L2ICD_defined_products(driver_module_param_file):
    CALDB_version = get_caldb_version(driver_module_param_file)

    filelist = glob("[F,N]UV*_transformations.txt")
    if len(filelist) > 0:
        for file in filelist:
            channel = file[:3]
            window_filter = file[4:10]
            print(f"\nWorking on {channel} {window_filter}")
            generate_ICD_UV_products(channel, window_filter, CALDB_version)

    combined_VIS_images = glob("combined_VIS_IM*.fits")
    for combined_VIS_image in combined_VIS_images:
        channel, window_filter = combined_VIS_image.split("_")[1:3]
        window_filter = window_filter[:6]
        print(f"\nWorking on {channel} {window_filter}")
        generate_ICD_VIS_products(channel, window_filter, CALDB_version)


now_time = datetime.utcnow().isoformat(sep="T", timespec="seconds")
level1_dataset = glob("LEVL1AS1UVT*.tar_V*2.gz")
if len(level1_dataset) == 1:
    level1_dataset = level1_dataset[0]
    level1_dataset = level1_dataset[:-3]
else:
    raise RuntimeError(
        f"Just a single Level1 dataset expected, please check {level1_dataset}"
    )

date_from_level1 = level1_dataset[11:19]
observation_ID = level1_dataset[19:40]
proposal_ID = observation_ID[:7]
target_ID = observation_ID[7:10]

if len(glob("pipeline/*txt")) == 1:
    driver_module_param_file = glob("pipeline/*txt")[0]
else:
    raise RuntimeError("Missing UVIT driver parameter file!")

parent_dir = f"{date_from_level1}_{observation_ID}_level2"
data_dir = f"{parent_dir}/uvit/data_products"
aux_data_dir = f"{parent_dir}/uvit/subsidiary_data_and_information"
combine_info_dir = f"{aux_data_dir}/combining_process_information"
episode_dir = f"{aux_data_dir}/episode_data"

shutil.rmtree(parent_dir, ignore_errors=True)
os.mkdir(parent_dir)
os.mkdir(f"{parent_dir}/uvit")
os.mkdir(data_dir)
os.mkdir(f"{data_dir}/FUV_data")
os.mkdir(f"{data_dir}/NUV_data")
os.mkdir(f"{data_dir}/VIS_data")
os.mkdir(aux_data_dir)
os.mkdir(combine_info_dir)
os.mkdir(f"{combine_info_dir}/astrometry_error_vector_plots")
os.mkdir(f"{combine_info_dir}/curvit_plots")
os.mkdir(f"{combine_info_dir}/logs")
os.mkdir(episode_dir)
os.mkdir(f"{episode_dir}/driver_module")

make_L2ICD_defined_products(driver_module_param_file)

for FUV_file in glob("AS1*uvtF*fits"):
    shutil.move(FUV_file, f"{data_dir}/FUV_data")

for NUV_file in glob("AS1*uvtN*fits"):
    shutil.move(NUV_file, f"{data_dir}/NUV_data")

for VIS_file in glob("AS1*uvtV*fits"):
    shutil.move(VIS_file, f"{data_dir}/VIS_data")

for file in glob("combined_*UV*CPS_coo_error_vectors.png"):
    shutil.copy(
        file,
        f"{combine_info_dir}/astrometry_error_vector_plots/{new_error_vector_name(file)}",
    )

for file in glob("curve_orbitwise_*png"):
    shutil.copy(
        file,
        f"{combine_info_dir}/curvit_plots/{new_lightcurve_name(file)}",
    )

for file in glob("source_[0-9]*png"):
    shutil.copy(
        file,
        f"{combine_info_dir}/curvit_plots/{new_fullfield_name(file)}",
    )

for file in glob("source_zoomed_*png"):
    shutil.copy(
        file,
        f"{combine_info_dir}/curvit_plots/{new_source_zoomed_name(file)}",
    )

shutil.copy("combining.log", f"{combine_info_dir}/logs/")
shutil.copy(driver_module_param_file, f"{episode_dir}/driver_module/")

for episode_file in glob("uvt_*[0-9]/*/*fits"):
    avoid_strings = ["artefacts", "A_l2err"]
    if not any(x in episode_file for x in avoid_strings):
        dst_dir = f"{episode_dir}/{os.path.dirname(episode_file)}"
        if not os.path.exists(dst_dir):
            os.makedirs(dst_dir)
        dst_path = f"{dst_dir}/{os.path.basename(episode_file)}"
        shutil.copy(episode_file, dst_path)

for uvt_dir in glob(f"{episode_dir}/uvt*"):
    dir_num = uvt_dir.split("uvt_")[-1]
    fuv_dir = f"{uvt_dir}/F_{dir_num}"
    if not os.path.exists(fuv_dir):
        os.makedirs(fuv_dir)
    nuv_dir = f"{uvt_dir}/N_{dir_num}"
    if not os.path.exists(nuv_dir):
        os.makedirs(nuv_dir)
    vis_dir = f"{uvt_dir}/V_{dir_num}"
    if not os.path.exists(vis_dir):
        os.makedirs(vis_dir)

print("\nDone")
