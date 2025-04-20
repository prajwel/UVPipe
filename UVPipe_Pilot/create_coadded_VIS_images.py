#!/usr/bin/env python3

import gc
import re
import numpy as np

from glob import glob
from scipy import stats
from joblib import Parallel, delayed
from astropy.io import fits
from pathlib import Path
from photutils.aperture import CircularAperture


zoom_factor = 3
all_VIS_RAS = glob("uvt_[0-9][0-9]/*/*dr.fits")
dir_num = 0
pad_width = 44


def get_framedata(file):
    hdu = fits.open(file, memmap=False)
    data = hdu[0].data
    hdu.close()
    return data


def get_framedata_joblib_loop(VIS_frames):
    all_frames_list = Parallel(n_jobs=-1, prefer="threads", max_nbytes=None)(
        delayed(get_framedata)(i) for i in VIS_frames
    )
    all_frames = np.array(all_frames_list).astype(np.single)
    del all_frames_list
    gc.collect()
    return all_frames


def tile_array(arr, factor):
    m, n = arr.shape
    return np.repeat(np.repeat(arr, factor, axis=1), factor, axis=0).reshape(
        m * factor, n * factor
    )


dimension = 512 + (2 * pad_width)

size = (dimension * zoom_factor) + 1
xx, yy = np.meshgrid(np.arange(1, size), np.arange(1, size))

detector_distortion_hdu = fits.open("detector_distortion_VIS.fits.gz")
detector_distortion_offsetx = -1 * detector_distortion_hdu[1].data
detector_distortion_offsety = -1 * detector_distortion_hdu[2].data
detector_distortion_offsetx = np.pad(
    detector_distortion_offsetx, pad_width=44, mode="constant", constant_values=0
)
detector_distortion_offsety = np.pad(
    detector_distortion_offsety, pad_width=44, mode="constant", constant_values=0
)
detector_distortion_offsetx_zoomed_image = tile_array(
    detector_distortion_offsetx, zoom_factor
)
detector_distortion_offsety_zoomed_image = tile_array(
    detector_distortion_offsety, zoom_factor
)
detector_distortion_offsetx_zoomed_image = (
    detector_distortion_offsetx_zoomed_image * zoom_factor
)
detector_distortion_offsety_zoomed_image = (
    detector_distortion_offsety_zoomed_image * zoom_factor
)

optical_distortion_hdu = fits.open("optical_distortion_VIS.fits.gz")
optical_distortion_offsetx = -1 * optical_distortion_hdu[1].data
optical_distortion_offsety = -1 * optical_distortion_hdu[2].data
optical_distortion_offsetx = np.pad(
    optical_distortion_offsetx, pad_width=44, mode="constant", constant_values=0
)
optical_distortion_offsety = np.pad(
    optical_distortion_offsety, pad_width=44, mode="constant", constant_values=0
)
optical_distortion_offsetx_zoomed_image = tile_array(
    optical_distortion_offsetx, zoom_factor
)
optical_distortion_offsety_zoomed_image = tile_array(
    optical_distortion_offsety, zoom_factor
)
optical_distortion_offsetx_zoomed_image = (
    optical_distortion_offsetx_zoomed_image * zoom_factor
)
optical_distortion_offsety_zoomed_image = (
    optical_distortion_offsety_zoomed_image * zoom_factor
)


def correct_single_frame(
    frame=None,
    xshift=None,
    yshift=None,
    frame_median=None,
    xx=xx,
    yy=yy,
    dimension=dimension,
    zoom_factor=zoom_factor,
    detector_distortion_offsetx_zoomed_image=detector_distortion_offsetx_zoomed_image,
    detector_distortion_offsety_zoomed_image=detector_distortion_offsety_zoomed_image,
    optical_distortion_offsetx_zoomed_image=optical_distortion_offsetx_zoomed_image,
    optical_distortion_offsety_zoomed_image=optical_distortion_offsety_zoomed_image,
):
    zoomed_frame = tile_array(frame, zoom_factor)
    frame_data = np.array(
        [xx.flatten() - 1, yy.flatten() - 1, zoomed_frame.flatten()]
    ).T
    frame_data = frame_data.astype(int)
    frame_data = frame_data[frame_data[:, 2] != 0]
    ddx = detector_distortion_offsetx_zoomed_image[frame_data[:, 1], frame_data[:, 0]]
    ddy = detector_distortion_offsety_zoomed_image[frame_data[:, 1], frame_data[:, 0]]
    x_plus_ddx = frame_data[:, 0] + ddx
    y_plux_ddy = frame_data[:, 1] + ddy
    frame_data = np.array([x_plus_ddx, y_plux_ddy, frame_data[:, 2]]).T
    detector_distortion_corrected_coo = np.round(frame_data[:, :2], 0).astype(int)
    odx = optical_distortion_offsetx_zoomed_image[
        detector_distortion_corrected_coo[:, 1], detector_distortion_corrected_coo[:, 0]
    ]
    ody = optical_distortion_offsety_zoomed_image[
        detector_distortion_corrected_coo[:, 1], detector_distortion_corrected_coo[:, 0]
    ]
    x_plus_ddx_plus_odx = detector_distortion_corrected_coo[:, 0] + np.array(odx)
    y_plux_ddy_plus_ody = detector_distortion_corrected_coo[:, 1] + np.array(ody)
    half_image_size = (dimension * zoom_factor) / 2
    x_plus_ddx_plus_odx = (
        x_plus_ddx_plus_odx + (x_plus_ddx_plus_odx - half_image_size) * 0.002
    )
    y_plux_ddy_plus_ody = (
        y_plux_ddy_plus_ody + (y_plux_ddy_plus_ody - half_image_size) * 0.0002
    )
    x_plus_ddx_plus_odx = x_plus_ddx_plus_odx - (xshift * zoom_factor)
    y_plux_ddy_plus_ody = y_plux_ddy_plus_ody - (yshift * zoom_factor)

    size = (dimension * zoom_factor) + 1
    bins = np.arange(0, size)
    dist_corrected_frame, _, _ = np.histogram2d(
        y_plux_ddy_plus_ody + 0.5,
        x_plus_ddx_plus_odx + 0.5,
        bins=(bins, bins),
        weights=frame_data[:, 2],
    )
    dist_corrected_frame = dist_corrected_frame + frame_median
    dist_corrected_frame = dist_corrected_frame + 50
    dist_corrected_frame[dist_corrected_frame < 0] = 50
    dist_corrected_frame[dist_corrected_frame > 65000] = 65535
    return np.round(dist_corrected_frame).astype(np.uint16)


def correct_single_frames_joblib_loop(all_vis_frames, xshifts, yshifts, frame_medians):
    transformed_frames = Parallel(n_jobs=-1, max_nbytes=None)(
        delayed(correct_single_frame)(frame, xshift, yshift, frame_median)
        for frame, xshift, yshift, frame_median in zip(
            all_vis_frames, xshifts, yshifts, frame_medians.flatten()
        )
    )

    return np.array(transformed_frames).astype(np.uint16)


def coadd_vis_image(VIS_RAS_file=None):
    VIS_RAS = fits.open(VIS_RAS_file)

    RA_pointing = VIS_RAS[0].header["RA_PNT"]
    DEC_pointing = VIS_RAS[0].header["DEC_PNT"]

    if "history" in VIS_RAS[0].header:
        history = VIS_RAS[0].header["history"]
        if len(history) > 0:
            re_result = re.search(r"V/uvtV\.\d{2}/", str(history))
            if type(re_result.group()) is str:
                ras_num = re_result.group()[2:9]
    else:
        print("History missing for events list in {}".format(VIS_RAS_file))

    VIS_frames = glob("output/*/*/*/" + ras_num + "/*/SignalFrames/*fits")

    # Extract time values
    VIS_frame_time_values = [
        float(re.search(r"_t(\d+\.\d+)_f", VIS_frame_file).group(1))
        for VIS_frame_file in VIS_frames
    ]

    VIS_frame_time_values = np.array(VIS_frame_time_values)
    VIS_frames = np.array(VIS_frames)
    VIS_frames = VIS_frames[np.argsort(VIS_frame_time_values)]
    VIS_frames = VIS_frames.tolist()

    all_vis_frames = get_framedata_joblib_loop(VIS_frames[1:])
    gc.collect()

    #    image_artefacts = np.median(all_vis_frames, axis = 0)
    image_artefacts, _ = stats.mode(all_vis_frames, axis=0, keepdims=False)
    hdu = fits.PrimaryHDU(image_artefacts)
    hdu.header["RA_PNT"] = RA_pointing
    hdu.header["DEC_PNT"] = DEC_pointing
    filename = Path(VIS_RAS_file).name.replace("_l2dr.fits", "I_l2img_artefacts.fits")
    hdu.writeto(Path(VIS_RAS_file).with_name(filename), overwrite=True)

    all_vis_frames = all_vis_frames - image_artefacts
    all_vis_frames = all_vis_frames * 10

    aperture = CircularAperture((256, 256), r=240)
    mask = aperture.to_mask(method="center")
    mask = mask.to_image(shape=((512, 512)))
    mask = mask.astype(bool)
    mask = np.logical_not(mask)
    mask = np.broadcast_to(mask, all_vis_frames.shape)
    all_vis_frames[mask] = 0

    frame_medians = np.median(all_vis_frames, axis=[1, 2])
    frame_medians = frame_medians[:, np.newaxis, np.newaxis]
    all_vis_frames = all_vis_frames - frame_medians

    all_vis_frames[mask] = 0

    pad_width = 44
    padded_all_vis_frames = np.pad(
        all_vis_frames, pad_width=pad_width, mode="constant", constant_values=0
    )
    padded_all_vis_frames = padded_all_vis_frames[
        pad_width:-pad_width
    ]  # remember that np.pad pads along all axes.
    all_vis_frames = None

    xshifts = VIS_RAS[2].data["X_shift"]
    yshifts = VIS_RAS[2].data["Y_shift"]

    transformed_frames = correct_single_frames_joblib_loop(
        padded_all_vis_frames, xshifts, yshifts, frame_medians
    )
    padded_all_vis_frames = None

    transformed_frames_median = np.median(transformed_frames, axis=0)
    lower_limit = transformed_frames_median - (50 + (0.1 * transformed_frames_median))
    upper_limit = transformed_frames_median + (50 + (0.1 * transformed_frames_median))
    lower_limit[lower_limit < 0] = 0
    upper_limit[upper_limit > 65535] = 65535
    lower_limit = lower_limit.astype(np.uint16)
    upper_limit = upper_limit.astype(np.uint16)
    mask = np.logical_and(
        transformed_frames >= lower_limit, transformed_frames <= upper_limit
    )

    mask_sum = np.sum(mask, axis=0)
    if np.sum(mask_sum < 1) > 0:
        print("Problem, one or more pixels masked across all frames")

    masked_transformed_frames = np.ma.masked_array(transformed_frames, ~mask)
    coadded_frame = np.ma.average(masked_transformed_frames, axis=0)
    coadded_frame = coadded_frame.data
    coadded_frame = coadded_frame.astype(np.float32)

    hdu = fits.PrimaryHDU(coadded_frame)
    hdu.header["RA_PNT"] = RA_pointing
    hdu.header["DEC_PNT"] = DEC_pointing
    filename = Path(VIS_RAS_file).name.replace("_l2dr.fits", "I_l2img.fits")
    hdu.writeto(Path(VIS_RAS_file).with_name(filename), overwrite=True)


for VIS_RAS_file in all_VIS_RAS:
    print(VIS_RAS_file)
    coadd_vis_image(VIS_RAS_file)

print("\nProcessing completed")
