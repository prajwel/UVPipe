#!/usr/bin/env python3


import warnings
import numpy as np
import matplotlib.pyplot as plt

from glob import glob
from astropy.io import fits
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.utils.exceptions import AstropyWarning
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture
from photutils.background import Background2D, MedianBackground

maximum_detections = 100
minimum_detections = 30


def detect_sources(data, threshold):
    aperture = CircularAperture((900, 900), r=650)
    mask = aperture.to_mask(method="center")
    mask = mask.to_image(shape=((1800, 1800)))
    data[~mask.astype(bool)] = np.nan

    sigma_clip = SigmaClip(sigma=3.0)
    bkg_estimator = MedianBackground()
    bkg = Background2D(
        data,
        (50, 50),
        filter_size=(3, 3),
        sigma_clip=sigma_clip,
        bkg_estimator=bkg_estimator,
    )

    daofind = DAOStarFinder(
        fwhm=1.5 * 3,
        threshold=threshold * bkg.background_rms_median,
        exclude_border=True,
    )
    sources = daofind(data - bkg.background)

    if sources is not None:
        sources.sort("mag")
        sources = sources["xcentroid", "ycentroid"]
    else:
        sources = []
    return sources


def VIS_detect_sources_in_range(
    data,
    minimum_detections=minimum_detections,
    maximum_detections=maximum_detections,
):
    threshold = 50
    sources_coo = detect_sources(data, threshold)
    while len(sources_coo) <= minimum_detections:
        sources_coo = detect_sources(data, threshold)
        threshold = threshold - (threshold / 3)
        if threshold <= 0.5:
            RuntimeError("Detection threshold went below 0.5!")
    sources_coo = sources_coo[:maximum_detections]
    print(f"{len(sources_coo)} sources detected")
    return sources_coo


def get_sources(image_file):
    image_hdu = fits.open(image_file)

    ra = image_hdu[0].header["RA_PNT"]
    dec = image_hdu[0].header["DEC_PNT"]
    print("ra={}".format(ra))
    print("dec={}".format(dec))

    warnings.simplefilter("ignore", category=AstropyWarning)
    coo = VIS_detect_sources_in_range(image_hdu[0].data)
    positions = np.transpose((coo["xcentroid"] + 1, coo["ycentroid"] + 1))

    if image_file[-2:] == "gz":
        filename = image_file.replace(".fits.gz", "_coo.dat")
    else:
        filename = image_file.replace(".fits", "_coo.dat")

    np.savetxt(filename, positions, fmt="%s", delimiter=",")

    mean, median, std = sigma_clipped_stats(image_hdu[0].data, sigma=5.0, maxiters=1)
    plt.imshow(image_hdu[0].data, cmap="Greys", origin="lower", vmax=median + 5 * std)
    apertures = CircularAperture(positions - 1, r=10.0)
    apertures.plot(color="blue", lw=1.5, alpha=0.5)

    if image_file[-2:] == "gz":
        filename = image_file.replace(".fits.gz", "_detections.png")
    else:
        filename = image_file.replace(".fits", "_detections.png")

    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close("all")


image_files = glob("combined_VIS_*.fits")
for image_file in image_files:
    print(f"\nWorking on {image_file}:")
    get_sources(image_file)
