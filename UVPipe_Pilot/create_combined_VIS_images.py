#!/usr/bin/env python3

import sys
import warnings
import numpy as np
import astroalign as aa

from glob import glob
from astropy.io import fits
from astropy.stats import SigmaClip
from astropy.convolution import convolve
from astropy.utils.exceptions import AstropyWarning
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture
from photutils.background import Background2D, MedianBackground
from joblib import Parallel, delayed
from skimage import transform
from aafitrans import find_transform, MaxIterError
from itertools import compress
from scipy.spatial import KDTree


def get_framedata(file):
    hdu = fits.open(file)
    return hdu[0].data


def get_framedata_joblib_loop(VIS_frames):
    all_frames = Parallel(n_jobs=-1)(delayed(get_framedata)(i) for i in VIS_frames)
    return all_frames


def opt_align_frames(matrix=None, frame=None, reference_frame=None):
    transf = transform.SimilarityTransform(matrix=matrix)
    transformed_frame, _ = aa.apply_transform(transf, frame, reference_frame)
    return transformed_frame


def get_rigid_transform(frame_sources, reference):
    try:
        transf, (source_list, target_list) = find_transform(
            frame_sources,
            reference,
            ttype="euclidean",
            seed=1,
            max_control_points=50,
            pixel_tolerance=1,
            min_matches=4,
            num_nearest_neighbors=12,
            kdtree_search_radius=0.02,
            n_samples=1,
            get_best_fit=True,
        )
        return transf.params
    except MaxIterError:
        return np.full((3, 3), np.nan)


def match_frames_joblib_loop(sources_all_frames):
    trans_matrices = Parallel(n_jobs=-1)(
        delayed(get_rigid_transform)(i, sources_all_frames[0])
        for i in sources_all_frames
    )
    return trans_matrices


def detect_sources(data, threshold):
    aperture = CircularAperture((900, 900), r=650)
    mask = aperture.to_mask(method="center")
    mask = mask.to_image(shape=((1800, 1800)))
    data[~mask.astype(bool)] = np.nan

    high_value_mask = data > 50000
    box_kernel = np.ones((15, 15))
    high_value_mask = convolve(high_value_mask, box_kernel, normalize_kernel=False)
    high_value_mask = high_value_mask > 0
    data[high_value_mask] = np.nan

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
        uA = np.array([sources["xcentroid"].data, sources["ycentroid"].data]).T
        uA = uA.astype(np.float64)
    else:
        uA = []
    return uA


def remove_close_pairs(detections, distance):
    tree = KDTree(detections)
    pairs = tree.query_pairs(r=distance)
    close_sources = np.unique(list(pairs))
    if len(close_sources) > 0:
        # print(f'Removed {len(close_sources)} sources.')
        detections = np.delete(detections, close_sources, axis=0)
    return detections


def VIS_detect_sources_in_range(
    VIS_images, minimum_detections=30, maximum_detections=50
):
    print("\nDetecting sources in VIS images:")
    detected_sources = []
    original_threshold = 500
    uA = []
    for path in VIS_images:
        hdu = fits.open(path)
        uA = []
        threshold = original_threshold
        while len(uA) <= minimum_detections:
            uA = detect_sources(hdu[0].data, threshold)
            uA = remove_close_pairs(uA, 5 / 1.13)
            threshold = threshold - (threshold / 3)
            if threshold <= 1:
                print("If you see this, please contact developer.")
                sys.exit()
        while len(uA) > maximum_detections:
            uA = detect_sources(hdu[0].data, threshold)
            uA = remove_close_pairs(uA, 5 / 1.13)
            threshold = 2 * threshold
            if threshold >= 1e20:
                print("If you see this, please contact developer.")
                sys.exit()

        print(f"{len(uA)} sources detected in {path}")
        detected_sources.append(uA)

    return detected_sources


def opt_make_coadded_vis_image(
    select_VIS_images_paths, trans_matrices, combined_vis_filename
):
    VIS_hdu = fits.open(select_VIS_images_paths[0])

    RA_pointing = VIS_hdu[0].header["RA_PNT"]
    DEC_pointing = VIS_hdu[0].header["DEC_PNT"]

    all_vis_frames = get_framedata_joblib_loop(select_VIS_images_paths)
    all_vis_frames = np.array(all_vis_frames)
    normalisation = np.mean(all_vis_frames, axis=(2, 1))

    def transform_frames_joblib_loop():
        transformed_frames = Parallel(n_jobs=-1, max_nbytes=None)(
            delayed(opt_align_frames)(matrix, frame, all_vis_frames[0])
            for matrix, frame in zip(trans_matrices, all_vis_frames)
        )

        return transformed_frames

    if len(all_vis_frames) > 1:
        transformed_frames = transform_frames_joblib_loop()
    else:
        transformed_frames = all_vis_frames

    all_vis_frames = None
    transformed_frames = np.array(transformed_frames)
    coadded_frame = np.average(transformed_frames, axis=0, weights=normalisation)
    coadded_frame = coadded_frame.astype(np.float32)
    hdu = fits.PrimaryHDU(coadded_frame)
    hdu.header["RA_PNT"] = RA_pointing
    hdu.header["DEC_PNT"] = DEC_pointing
    hdu.header["N_EPISDS"] = len(select_VIS_images_paths)

    expstarts = [fits.getheader(path)["EXPSTART"] for path in select_VIS_images_paths]
    expends = [fits.getheader(path)["EXPEND"] for path in select_VIS_images_paths]
    exptimes = [fits.getheader(path)["EXP_TIME"] for path in select_VIS_images_paths]
    hdu.header["EXPSTART"] = np.min(expstarts)
    hdu.header["EXPEND"] = np.max(expends)
    hdu.header["EXP_TIME"] = np.sum(exptimes)

    hdu.writeto(combined_vis_filename, overwrite=True)
    print(f"\nCombined VIS image: {combined_vis_filename}")


def make_combined_VIS_images(ref_episode=None):
    files_list = glob("uvt_*/V*/*img.fits")
    if len(files_list) > 0:
        window_filter_list = [file[42:48] for file in files_list]
        window_filter_list = list(set(window_filter_list))
        for window_filter in window_filter_list:
            print(f"\nCombining VIS {window_filter} files:")
            combined_vis_filename = f"combined_VIS_{window_filter}.fits"
            chosen_window_filter_paths = glob(f"uvt*/V*/*{window_filter}*img.fits")
            expstarts = [
                fits.getheader(path)["EXPSTART"] for path in chosen_window_filter_paths
            ]
            _, unique_index = np.unique(expstarts, return_index=True)
            chosen_window_filter_paths = np.array(chosen_window_filter_paths)[
                unique_index
            ]
            chosen_window_filter_paths.sort()

            chosen_window_filter_paths = list(chosen_window_filter_paths)
            if ref_episode:
                indices = [
                    i
                    for i, s in enumerate(chosen_window_filter_paths)
                    if ref_episode in s
                ]
                if len(indices) == 1:
                    index = indices[0]
                    chosen_window_filter_paths.insert(
                        0, chosen_window_filter_paths.pop(index)
                    )
                else:
                    raise RuntimeError(f"Please check {ref_episode}.")

            print(np.array(chosen_window_filter_paths))
            detected_sources = VIS_detect_sources_in_range(chosen_window_filter_paths)
            trans_matrices = match_frames_joblib_loop(detected_sources)

            # To mask failed cases.
            trans_matrices = np.array(trans_matrices)
            mask = np.logical_not(np.isnan(trans_matrices[:, 0, 0]))
            trans_matrices = trans_matrices[mask]

            if False in mask:
                print("\nFailed to find VIS transformations for:")
                print(np.array(chosen_window_filter_paths)[~mask])

            chosen_window_filter_paths = list(
                compress(chosen_window_filter_paths, mask)
            )

            with np.printoptions(suppress=True, formatter={"float": "{:.5f}".format}):
                print("\nVIS channel derived transformation matrices:")
                print(np.array(trans_matrices))
            opt_make_coadded_vis_image(
                chosen_window_filter_paths, trans_matrices, combined_vis_filename
            )


warnings.simplefilter("ignore", category=AstropyWarning)
make_combined_VIS_images()
