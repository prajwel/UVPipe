#!/usr/bin/env python3

import warnings
import numpy as np
import matplotlib.pyplot as plt

from glob import glob
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture, RectangularAperture
from photutils.centroids import centroid_com
from photutils.background import Background2D, MedianBackground
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.convolution import Gaussian2DKernel, convolve


exp_map_cutoff_fraction = 0.8
maximum_detections = 200
minimum_detections = 30
string_to_window_size = {
    "00": 512,
    "22": 100,
    "33": 150,
    "44": 200,
    "55": 250,
    "66": 300,
    "77": 350,
}


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


def detect_sources(data, threshold, exp_map_mask, outer_radius=2020, inner_radius=2000):
    xcentre, ycentre = centroid_com(~exp_map_mask)
    aperture = CircularAperture((xcentre, ycentre), r=outer_radius)
    outer_radius_mask = aperture.to_mask(method="center")
    outer_radius_mask = outer_radius_mask.to_image(shape=((4800, 4800)))
    outer_radius_mask = np.logical_not(outer_radius_mask)

    sigma_clip = SigmaClip(sigma=3.0)
    bkg_estimator = MedianBackground()
    for_bkg_binned = rebin(data, 64)
    for_bkg_zoomed = tile_array(for_bkg_binned, 64)
    bkg = Background2D(
        for_bkg_zoomed,
        (384, 384),
        filter_size=(3, 3),
        sigma_clip=sigma_clip,
        bkg_estimator=bkg_estimator,
        coverage_mask=outer_radius_mask,
    )

    kernel = Gaussian2DKernel(x_stddev=1.5, x_size=5)
    smoothed_data = convolve(data, kernel, nan_treatment="fill", mask=outer_radius_mask)

    threshold = threshold * np.median(bkg.background_rms[bkg.background_rms > 0])
    daofind = DAOStarFinder(fwhm=3.33, threshold=threshold, exclude_border=True)
    sources = daofind(smoothed_data - bkg.background, mask=outer_radius_mask)

    source_pixel_locations = np.round(
        np.array([sources["xcentroid"], sources["ycentroid"]])
    )
    source_pixel_locations = source_pixel_locations.astype("int")

    aperture = CircularAperture((xcentre, ycentre), r=inner_radius)
    inner_radius_mask = aperture.to_mask(method="center")
    inner_radius_mask = inner_radius_mask.to_image(shape=((4800, 4800)))
    inner_radius_mask = np.logical_not(inner_radius_mask)

    combined_mask = np.logical_or(exp_map_mask, inner_radius_mask)

    sources_to_keep_mask = combined_mask[
        source_pixel_locations[1], source_pixel_locations[0]
    ]
    sources_to_keep_mask = np.logical_not(sources_to_keep_mask)
    sources = sources[sources_to_keep_mask]

    if sources is not None:
        sources.sort("mag")
        sources = sources["xcentroid", "ycentroid"]
    else:
        sources = []
    return sources, (xcentre, ycentre)


def detect_sources_window_mode(data, threshold, exp_map_mask, window_size):
    outer_size = window_size * 8 + 20
    inner_size = window_size * 8
    xcentre, ycentre = centroid_com(~exp_map_mask)
    aperture = RectangularAperture((xcentre, ycentre), outer_size, outer_size)
    outer_mask = aperture.to_mask(method="center")
    outer_mask = outer_mask.to_image(shape=((4800, 4800)))
    outer_mask = np.logical_not(outer_mask)

    sigma_clip = SigmaClip(sigma=3.0)
    bkg_estimator = MedianBackground()
    for_bkg_binned = rebin(data, 64)
    for_bkg_zoomed = tile_array(for_bkg_binned, 64)
    bkg = Background2D(
        for_bkg_zoomed,
        (384, 384),
        filter_size=(3, 3),
        sigma_clip=sigma_clip,
        bkg_estimator=bkg_estimator,
        coverage_mask=outer_mask,
    )

    kernel = Gaussian2DKernel(x_stddev=1.5, x_size=5)
    smoothed_data = convolve(data, kernel, nan_treatment="fill", mask=outer_mask)

    threshold = threshold * np.median(bkg.background_rms[bkg.background_rms > 0])
    daofind = DAOStarFinder(fwhm=3.33, threshold=threshold, exclude_border=True)
    sources = daofind(smoothed_data - bkg.background, mask=outer_mask)

    source_pixel_locations = np.round(
        np.array([sources["xcentroid"], sources["ycentroid"]])
    )
    source_pixel_locations = source_pixel_locations.astype("int")

    aperture = RectangularAperture((xcentre, ycentre), inner_size, inner_size)
    inner_mask = aperture.to_mask(method="center")
    inner_mask = inner_mask.to_image(shape=((4800, 4800)))
    inner_mask = np.logical_not(inner_mask)

    combined_mask = np.logical_or(exp_map_mask, inner_mask)

    sources_to_keep_mask = combined_mask[
        source_pixel_locations[1], source_pixel_locations[0]
    ]
    sources_to_keep_mask = np.logical_not(sources_to_keep_mask)
    sources = sources[sources_to_keep_mask]

    if sources is not None:
        sources.sort("mag")
        sources = sources["xcentroid", "ycentroid"]
    else:
        sources = []
    return sources, (xcentre, ycentre)


def UV_detect_sources_in_range(
    data,
    exp_map_mask,
    window_size,
    minimum_detections=minimum_detections,
    maximum_detections=maximum_detections,
):
    threshold = 10
    if window_size == 512:
        sources_coo, centre = detect_sources(data, threshold, exp_map_mask)
        while len(sources_coo) <= minimum_detections:
            sources_coo, centre = detect_sources(data, threshold, exp_map_mask)
            threshold = threshold - (threshold / 3)
            if threshold <= 0.5:
                RuntimeError("Detection threshold went below 0.5!")
    else:
        sources_coo, centre = detect_sources_window_mode(
            data, threshold, exp_map_mask, window_size
        )
        while len(sources_coo) <= minimum_detections:
            sources_coo, centre = detect_sources_window_mode(
                data, threshold, exp_map_mask, window_size
            )
            threshold = threshold - (threshold / 3)
            if threshold <= 0.5:
                RuntimeError("Detection threshold went below 0.5!")
    sources_coo = sources_coo[:maximum_detections]
    print(f"{len(sources_coo)} sources detected")
    return sources_coo, centre


def get_sources(image_file, exp_map_file):
    image_hdu = fits.open(image_file)
    exp_map_hdu = fits.open(exp_map_file)

    _, channel, window_filter, _ = image_file.split("_")
    window_id = window_filter[2:4]
    window_size = string_to_window_size[window_id]

    # mask creation
    exp_threshold = exp_map_cutoff_fraction * np.median(
        exp_map_hdu[0].data[exp_map_hdu[0].data > 0]
    )
    exp_map_mask = exp_map_hdu[0].data < exp_threshold

    ra = image_hdu[0].header["RA_PNT"]
    dec = image_hdu[0].header["DEC_PNT"]
    print("ra={}".format(ra))
    print("dec={}".format(dec))

    warnings.simplefilter("ignore", category=AstropyWarning)
    coo, centre = UV_detect_sources_in_range(
        image_hdu[0].data, exp_map_mask, window_size
    )
    positions = np.transpose((coo["xcentroid"] + 1, coo["ycentroid"] + 1))

    if image_file[-2:] == "gz":
        filename = image_file.replace(".fits.gz", "_coo.dat")
    else:
        filename = image_file.replace(".fits", "_coo.dat")

    np.savetxt(filename, positions, fmt="%s", delimiter=",")

    kernel = Gaussian2DKernel(x_stddev=1.5, x_size=5)
    binned_data = rebin(image_hdu[0].data, 2)
    binned_mask = rebin(exp_map_mask, 2)
    smoothed_data = convolve(
        binned_data, kernel, nan_treatment="fill", mask=binned_mask
    )
    zoomed_data = tile_array(smoothed_data, 2)

    mean, median, std = sigma_clipped_stats(zoomed_data, sigma=5.0, maxiters=1)
    plt.imshow(zoomed_data, cmap="Greys", origin="lower", vmax=median + 5 * std)
    apertures = CircularAperture(positions - 1, r=10.0)
    apertures.plot(color="blue", lw=1.5, alpha=0.5)

    if window_size == 512:
        apertures = CircularAperture(centre, r=2000)
        apertures.plot(color="red", lw=1.5)
    else:
        apertures = RectangularAperture(centre, window_size * 8, window_size * 8)
        apertures.plot(color="red", lw=1.5)

    if image_file[-2:] == "gz":
        filename = image_file.replace(".fits.gz", "_detections.png")
    else:
        filename = image_file.replace(".fits", "_detections.png")

    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close("all")


image_files = glob("combined_*CPS.fits*")

for image_file in image_files:
    print(f"\nWorking on {image_file}:")
    exp_map_file = image_file.replace("CPS", "exp_map")
    get_sources(image_file, exp_map_file)
