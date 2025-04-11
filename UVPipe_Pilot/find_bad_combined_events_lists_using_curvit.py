#!/usr/bin/env python3

import curvit
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from glob import glob
from functools import partial
from astropy.io import fits
from astropy.stats import SigmaClip
from photutils.aperture import CircularAperture, aperture_photometry
from photutils.background import Background2D, MedianBackground
from matplotlib.colors import LogNorm


time_separation = 1 * 60  # seconds
image_size = 30
aperture_radius = 2.5


def rebin(arr, bin_factor):
    shape = (
        int(arr.shape[0] / bin_factor),
        bin_factor,
        int(arr.shape[1] / bin_factor),
        bin_factor,
    )
    binned_arr = arr.reshape(shape).mean(-1).mean(1)
    return binned_arr


def tile_array(arr, factor):
    m, n = arr.shape
    return np.repeat(np.repeat(arr, factor, axis=1), factor, axis=0).reshape(
        m * factor, n * factor
    )


def get_SEERs(coo_file):
    _, channel, window_filter, _, _ = coo_file.split("_")
    brightest_source_coos = np.genfromtxt(coo_file, delimiter=",")[:2]
    brightest_source_coos = brightest_source_coos - 1
    brightest_source_coos = np.round(brightest_source_coos, 2)

    CPS_image = coo_file[:-12] + "_CPS.fits"
    exp_map = coo_file[:-12] + "_exp_map_in_seconds.fits"
    hdu = fits.open(CPS_image, mode="update")
    exp_map_hdu = fits.open(exp_map)

    exp_map_mask = exp_map_hdu[0].data < 5
    sigma_clip = SigmaClip(sigma=3.0)
    bkg_estimator = MedianBackground()
    for_bkg_binned = rebin(hdu[0].data, 64)
    for_bkg_zoomed = tile_array(for_bkg_binned, 64)
    bkg = Background2D(
        for_bkg_zoomed,
        (384, 384),
        filter_size=(3, 3),
        sigma_clip=sigma_clip,
        bkg_estimator=bkg_estimator,
        coverage_mask=exp_map_mask,
    )

    small_aperture = CircularAperture(brightest_source_coos, r=2.5)
    big_aperture = CircularAperture(brightest_source_coos, r=5)
    phot_table_small_aperture = aperture_photometry(hdu[0].data, small_aperture)
    phot_table_big_aperture = aperture_photometry(
        hdu[0].data - bkg.background, big_aperture
    )
    SEERs = (
        phot_table_small_aperture["aperture_sum"]
        / phot_table_big_aperture["aperture_sum"]
    )
    SEER_1, SEER_2 = np.round(SEERs, 2)
    print(f"\n{channel} {window_filter}:")
    print(f"SEER_1, SEER_2 = {SEER_1}, {SEER_2}")
    hdu[0].header["SEER_1"] = SEER_1
    hdu[0].header["SEER_2"] = SEER_2
    hdu.close()
    return SEER_1, SEER_2


def variability_measure(filename):
    data = np.genfromtxt(filename)
    if len(data.shape) == 1:
        R = 0
    else:
        max_min_drop = np.max(data[:, 1]) - np.min(data[:, 1])
        # max_min_drop = np.mean(data[:, 1]) - np.min(data[:, 1])
        mean = np.mean(data[:, 1])
        R = max_min_drop / mean
        R = np.round(R, 2)
    return R


def read_columns(events_list):
    # Reading few columns.
    hdu = fits.open(events_list)
    time = hdu[1].data["MJD_L2"]
    fx = hdu[1].data["Fx"]
    fy = hdu[1].data["Fy"]
    photons = hdu[1].data["EFFECTIVE_NUM_PHOTONS"]
    bad_flag = hdu[1].data["BAD FLAG"]
    mask = photons > 0
    mask = np.logical_and(mask, bad_flag)
    if "FrameCount" in hdu[1].data.names:
        first_frame_mask = hdu[1].data["FrameCount"] != 1
        mask = np.logical_and(mask, first_frame_mask)
    time = time[mask]
    fx = fx[mask]
    fy = fy[mask]
    photons = photons[mask]
    return time, fx, fy, photons


# To change mission elapsed time in seconds to modified julian date.
def met_to_mjd(met):
    jan2010 = 55197.0  # 2010.0(UTC) expressed with MJD format and scale UTC.
    mjd = (met / 86400.0) + jan2010  # 1 julian day = 86400 seconds.
    return mjd


def choose_events(starts, ends, time, fx, fy, weights):
    for start, end in zip(starts, ends):
        XY = [
            (t, x, y, w)
            for t, x, y, w in zip(time, fx, fy, weights)
            if start <= t <= end
        ]
        yield XY


def create_sub_image(
    zipped_data,
    x_centre=None,
    y_centre=None,
    image_size=image_size,
    aperture_radius=aperture_radius,
):
    ax.clear()
    global orbit_no
    left_x = x_centre - image_size / 2
    right_x = x_centre + image_size / 2
    bottom_y = y_centre - image_size / 2
    top_y = y_centre + image_size / 2
    time, fx, fy, weights = zip(*zipped_data)
    mean_MJD_time = np.round(met_to_mjd(np.mean(time)), 2)
    H, _, _ = np.histogram2d(
        fx,
        fy,
        bins=image_size,
        range=[[left_x, right_x], [bottom_y, top_y]],
        weights=weights,
    )
    if np.max(H) != 0:
        norm = LogNorm(vmin=np.mean(H) / 2, vmax=np.max(H))
        ax.imshow(H, interpolation="nearest", origin="lower", norm=norm)
    else:
        ax.imshow(H, interpolation="nearest", origin="lower")
    obj_circle = plt.Circle(
        (image_size / 2 - 1.5, image_size / 2 - 1.5),
        aperture_radius,
        color="k",
        fill=False,
    )
    ax.add_artist(obj_circle)
    ax.set_title(
        f"Orbit = {orbit_no}, MJD = {mean_MJD_time}, radius = {aperture_radius} px"
    )
    if orbit_no != 0:
        fig.savefig(
            f"{channel}_{window_filter}_X_{x_centre}_Y_{y_centre}_orbit_{orbit_no}.png",
            dpi=150,
        )
    orbit_no = orbit_no + 1


warnings.filterwarnings("ignore")

coo_files = glob("combined_[F,N]*coo.dat")
for coo_file in coo_files:
    _, channel, window_filter, _, _ = coo_file.split("_")
    brightest_source_coos = np.genfromtxt(coo_file, delimiter=",")[:2]
    brightest_source_coos = np.round(brightest_source_coos, 2)
    events_list = coo_file[:-12] + "_events_list.fits"
    hdu = fits.open(events_list)
    framecount_per_sec = hdu[0].header["AVGFRMRT"]
    time, fx, fy, photons = read_columns(events_list)
    weights = photons / framecount_per_sec
    for source_coo in brightest_source_coos:
        x_centre, y_centre = source_coo
        curvit.curve_orbitwise(
            events_list=events_list,
            xp=x_centre - 1,
            yp=y_centre - 1,
            radius=aperture_radius,
            framecount_per_sec=framecount_per_sec,
            time_separation=time_separation,
        )
        mask = ((fx - x_centre) ** 2 + (fy - y_centre) ** 2) <= image_size**2
        masked_time = time[mask]
        masked_fx = fx[mask]
        masked_fy = fy[mask]
        masked_weights = weights[mask]
        unique_time = np.unique(masked_time)
        unique_time = np.sort(unique_time)
        diffs = np.diff(unique_time)
        boundaries = np.where(diffs > time_separation)[0] + 1
        starts = np.concatenate([[0], boundaries])
        ends = np.concatenate([boundaries - 1, [len(unique_time) - 1]])
        starts = unique_time[starts]
        ends = unique_time[ends]
        orbit_no = 0
        fig, ax = plt.subplots()
        ani = animation.FuncAnimation(
            fig,
            partial(create_sub_image, x_centre=x_centre, y_centre=y_centre),
            choose_events(
                starts, ends, masked_time, masked_fx, masked_fy, masked_weights
            ),
            blit=False,
            interval=400,
            repeat=True,
            cache_frame_data=False,
        )

        gif_filename = (
            f"{channel}_{window_filter}_X_{x_centre}_Y_{y_centre}_animation.gif"
        )
        ani.save(gif_filename, writer="imagemagick")


datfiles = glob("curve*dat")
sortorder = np.argsort([datfile.split("combined_")[1] for datfile in datfiles])
datfiles = np.array(datfiles)[sortorder]
for datfile in datfiles:
    R = variability_measure(datfile)
    print("{}, R = {}".format(datfile, R))

for coo_file in coo_files:
    SEER_1, SEER_2 = get_SEERs(coo_file)
