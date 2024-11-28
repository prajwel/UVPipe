#!/usr/bin/env python3

import sys
import ntpath
import warnings
import numpy as np
import astroalign as aa
import matplotlib.pyplot as plt

from glob import glob
from pathlib import Path
from astropy.io import fits
from astropy.stats import SigmaClip
from astropy.convolution import Gaussian2DKernel, convolve, Box2DKernel
from astropy.utils.exceptions import AstropyWarning
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture
from photutils.background import Background2D, MedianBackground
from scipy import interpolate
from scipy.spatial import KDTree
from joblib import Parallel, delayed
from skimage import transform
from aafitrans import find_transform, MaxIterError
from itertools import compress


# To change mission elapsed time in seconds to modified julian date.
def met_to_mjd(met):
    jan2010 = 55197.0  # 2010.0(UTC) expressed with MJD format and scale UTC.
    mjd = (met / 86400.0) + jan2010  # 1 julian day = 86400 seconds.
    return mjd


def get_MJD_start_time(path):
    events_list = glob(f"{path}/*ce.fits")[0]
    hdu = fits.open(events_list)
    median_MJD = met_to_mjd(np.median(hdu[1].data["MJD_L2"]))
    return median_MJD


def get_exptime(filename):
    hdu = fits.open(filename)
    exptime = hdu[0].header["EXP_TIME"]
    if exptime <= 0:
        print(f"Warning! EXP_TIME is {exptime}")
    return exptime


def plot_transform_table(transform_table, transformations_filename, median_MJDs):
    sort_order = np.argsort(median_MJDs)
    median_MJDs = median_MJDs[sort_order]
    transform_table = transform_table[sort_order]
    sin_theta = transform_table[:, 2].astype(np.double)
    theta_degrees = np.rad2deg(np.arcsin(sin_theta))

    episode = [int(d.split("_")[-1]) for d in transform_table[:, 0]]

    fig, ax = plt.subplots()
    axplot1 = ax.plot(
        median_MJDs, theta_degrees, label="$\\theta$", marker=".", color="tab:blue"
    )

    ax2 = ax.twinx()
    ax2plot1 = ax2.plot(
        median_MJDs,
        transform_table[:, 3].astype(np.double),
        label="X shift",
        marker=".",
        color="tab:orange",
    )
    ax2plot2 = ax2.plot(
        median_MJDs,
        transform_table[:, 4].astype(np.double),
        label="Y shift",
        marker=".",
        color="tab:green",
    )
    ax.set_xlabel("MJD")
    ax.set_ylabel("Degrees")
    ax2.set_ylabel("Sub-pixels")

    axplots = axplot1 + ax2plot1 + ax2plot2
    labels = [axplot.get_label() for axplot in axplots]
    ax.legend(axplots, labels)

    fig.savefig(transformations_filename.replace("txt", "png"), dpi=150)


def read_columns(events_list):
    # Reading few columns.
    f = fits.open(events_list)
    time = f[1].data["MJD_L2"]
    fx = f[1].data["Fx"]
    fy = f[1].data["Fy"]
    photons = f[1].data["EFFECTIVE_NUM_PHOTONS"]
    bad_flag = f[1].data["BAD FLAG"]
    mask = photons > 0
    mask = np.logical_and(mask, bad_flag)
    time = time[mask]
    fx = fx[mask]
    fy = fy[mask]
    photons = photons[mask]
    return time, fx, fy, photons


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
    hdu.writeto(combined_vis_filename, overwrite=True)
    print(f"\nCombined VIS image: {combined_vis_filename}")


class correlate_events_lists:
    def __init__(self, reference_events_list, to_match_events_list):
        self.reference_events_list = reference_events_list
        self.to_match_events_list = to_match_events_list
        self.upper_limit = 160

    def get_columns(self, hdu):
        time = hdu[1].data["MJD_L2"]
        fx = hdu[1].data["Fx"]
        fy = hdu[1].data["Fy"]
        photons = hdu[1].data["EFFECTIVE_NUM_PHOTONS"]
        bad_flag = hdu[1].data["BAD FLAG"]
        mask = photons > 0
        mask = np.logical_and(mask, bad_flag)
        time = time[mask]
        fx = fx[mask]
        fy = fy[mask]
        photons = photons[mask]
        photons = photons / np.median(photons[photons > 0])
        return time, fx, fy, photons

    def get_image(self, hdu, bins_range):
        time, fx, fy, photons = self.get_columns(hdu)
        bins = np.arange(*bins_range)
        ndarray, yedges, xedges = np.histogram2d(
            fy, fx, bins=(bins, bins), weights=photons
        )
        return ndarray

    def get_lag(self, correlations, shift_range):
        f = interpolate.interp1d(shift_range, correlations, kind="quadratic")
        new_range = np.linspace(
            shift_range[0],
            shift_range[-1],
            num=(self.upper_limit * 2 + 1) * 100,
            endpoint=True,
        )
        interpolated_correlations = f(new_range)
        max_index = np.argmax(interpolated_correlations)
        peak_position = new_range[max_index]
        return peak_position

    def get_shifts(self, bin_size=1, bg_bin_size=64, window_filter_episode=None):
        # To get shifts
        reference_image = self.get_image(
            self.reference_events_list, (0, 4801, bin_size)
        )
        to_match_image = self.get_image(self.to_match_events_list, (0, 4801, bin_size))

        box_kernel = Box2DKernel(5)
        boxed_reference_image = convolve(reference_image, box_kernel)
        reference_image_background_mask = boxed_reference_image <= (4 / 25)
        reference_image[reference_image_background_mask] = 0

        boxed_to_match_image = convolve(to_match_image, box_kernel)
        to_match_image_background_mask = boxed_to_match_image <= (4 / 25)
        to_match_image[to_match_image_background_mask] = 0

        kernel = Gaussian2DKernel(x_stddev=1)
        smoothed_reference_image = convolve(reference_image, kernel)
        smoothed_to_match_image = convolve(to_match_image, kernel)

        smoothed_reference_image = smoothed_reference_image / np.std(
            smoothed_reference_image
        )
        smoothed_to_match_image = smoothed_to_match_image / np.std(
            smoothed_to_match_image
        )

        smoothed_reference_image[smoothed_reference_image < 0] = 0
        smoothed_to_match_image[smoothed_to_match_image < 0] = 0

        smoothed_reference_X = np.sum(smoothed_reference_image, axis=0)
        smoothed_reference_Y = np.sum(smoothed_reference_image, axis=1)

        smoothed_to_match_X = np.sum(smoothed_to_match_image, axis=0)
        smoothed_to_match_Y = np.sum(smoothed_to_match_image, axis=1)

        shift_limit = int(self.upper_limit / bin_size)
        shift_range = np.arange(-shift_limit, shift_limit + 1)

        X_correlations = []
        Y_correlations = []
        for shift in shift_range:
            X_product = np.sum(
                smoothed_reference_X * np.roll(smoothed_to_match_X, shift)
            )
            Y_product = np.sum(
                smoothed_reference_Y * np.roll(smoothed_to_match_Y, shift)
            )
            X_correlations.append(X_product)
            Y_correlations.append(Y_product)

        x_shift = self.get_lag(X_correlations, shift_range) * bin_size
        y_shift = self.get_lag(Y_correlations, shift_range) * bin_size

        fig, ax = plt.subplots()
        ax.plot(shift_range, X_correlations, label="X correlations")
        ax.plot(shift_range, Y_correlations, label="Y correlations")
        ax.set_xlabel("Shifts in sub-pixels")
        ax.set_ylabel("Correlation")
        ax.legend()
        corrplot_filename = window_filter_episode + "_correlations.png"
        fig.savefig(corrplot_filename, dpi=150)
        plt.close("all")

        return x_shift, y_shift


def combine_event_lists(
    paths, matrices, channel, combined_eventslist_name, transformations_filename
):
    transform_table = []
    framerate_from_header = []
    shift_x_corr = 0
    shift_y_corr = 0
    i = 0
    for path, matrix in zip(paths, matrices):
        VIS_to_UV = VIS_to_UV_matrix[channel]

        # trans_matrices = orig_trans_matrices
        cos_rot = matrix[0, 0]
        sin_rot = matrix[1, 0]
        shift_x = matrix[0, 2] * 8 / 3
        shift_y = matrix[1, 2] * 8 / 3

        shift_x_UV_orientation = shift_x * VIS_to_UV[0, 0] + shift_y * VIS_to_UV[0, 1]
        shift_y_UV_orientation = shift_x * VIS_to_UV[1, 0] + shift_y * VIS_to_UV[1, 1]

        img_path = ntpath.split(path)[0] + "/*I_l2img*"
        eventslist_header = fits.getheader(path)
        if len(glob(img_path)) == 1:
            img_hdu = fits.open(glob(img_path)[0])
            framerate_from_header.append(1 / img_hdu[0].header["INT_TIME"])
            RA_pointing = img_hdu[0].header["RA_PNT"]
            DEC_pointing = img_hdu[0].header["DEC_PNT"]
        elif "AVGFRMRT" in eventslist_header:
            framerate_from_header.append(eventslist_header["AVGFRMRT"])
            RA_pointing = eventslist_header["RA_PNT"]
            DEC_pointing = eventslist_header["DEC_PNT"]
        else:
            print(
                "Requires exactly one file matching *I_l2img* in directory {}!".format(
                    ntpath.split(path)[0]
                )
            )
            sys.exit()

        if i == 0:
            hdu_base = fits.open(path)

            photons = hdu_base[1].data["EFFECTIVE_NUM_PHOTONS"]
            mask = photons > 0

            nrows_base = hdu_base[1].data[mask].shape[0]

            X_and_Y = (
                np.array([hdu_base[1].data["Fx"], hdu_base[1].data["Fy"]]).T - 2400
            )
            X_centroid = (
                X_and_Y[:, 0] * cos_rot
                - X_and_Y[:, 1] * sin_rot
                + shift_x_UV_orientation
            )
            Y_centroid = (
                X_and_Y[:, 0] * sin_rot
                + X_and_Y[:, 1] * cos_rot
                + shift_y_UV_orientation
            )

            hdu_base[1].data["Fx"] = X_centroid + 2400
            hdu_base[1].data["Fy"] = Y_centroid + 2400

            hdu = fits.BinTableHDU.from_columns(
                hdu_base[1].columns, nrows=nrows_base, fill=True
            )
            for colname in hdu_base[1].columns.names:
                hdu.data[colname][:nrows_base] = hdu_base[1].data[colname][mask]

            hduP = fits.PrimaryHDU()
            hduP.header = hdu_base[0].header
            hdu.name = hdu_base[1].name
            hdu_base = fits.HDUList([hduP, hdu])
            hdu_base.writeto(combined_eventslist_name, overwrite=True)

        else:
            hdu_add = fits.open(path)
            hdu_base = fits.open(combined_eventslist_name)

            photons = hdu_add[1].data["EFFECTIVE_NUM_PHOTONS"]
            mask = photons > 0

            nrows_base = hdu_base[1].data.shape[0]
            nrows_add = hdu_add[1].data[mask].shape[0]
            nrows = nrows_base + nrows_add

            X_and_Y = np.array([hdu_add[1].data["Fx"], hdu_add[1].data["Fy"]]).T - 2400
            X_centroid = (
                X_and_Y[:, 0] * cos_rot
                - X_and_Y[:, 1] * sin_rot
                + shift_x_UV_orientation
            )
            Y_centroid = (
                X_and_Y[:, 0] * sin_rot
                + X_and_Y[:, 1] * cos_rot
                + shift_y_UV_orientation
            )

            hdu_add[1].data["Fx"] = X_centroid + 2400
            hdu_add[1].data["Fy"] = Y_centroid + 2400

            window_filter_episode = (
                transformations_filename[:10] + "_" + Path(path).parts[-2]
            )

            init_correlation = correlate_events_lists(hdu_base, hdu_add)
            shift_x_corr, shift_y_corr = init_correlation.get_shifts(
                window_filter_episode=window_filter_episode
            )

            hdu_add[1].data["Fx"] = hdu_add[1].data["Fx"] + shift_x_corr
            hdu_add[1].data["Fy"] = hdu_add[1].data["Fy"] + shift_y_corr

            hdu = fits.BinTableHDU.from_columns(
                hdu_base[1].columns, nrows=nrows, fill=True
            )
            for colname in hdu_base[1].columns.names:
                hdu.data[colname][:nrows_base] = hdu_base[1].data[colname]
                hdu.data[colname][nrows_base:] = hdu_add[1].data[colname][mask]

            hdu.name = hdu_base[1].name
            hdu_base = fits.HDUList([hduP, hdu])
            hdu_base.writeto(combined_eventslist_name, overwrite=True)
            # break

        AVGFRMRT = np.mean(framerate_from_header)
        STDFRMRT = np.std(framerate_from_header)

        hdu_base = fits.open(combined_eventslist_name, mode="update")
        hdu_base[0].header["AVGFRMRT"] = AVGFRMRT
        hdu_base[0].header["STDFRMRT"] = STDFRMRT
        hdu_base[0].header["RA_PNT"] = RA_pointing
        hdu_base[0].header["DEC_PNT"] = DEC_pointing
        hdu_base[0].header["COMBMETH"] = "ALT_METHOD1"
        hdu_base.close()

        transform_table.append(
            [
                ntpath.split(path)[0],
                cos_rot,
                sin_rot,
                shift_x_UV_orientation + shift_x_corr,
                shift_y_UV_orientation + shift_y_corr,
            ]
        )
        i = i + 1

    print(f"Combined events list: {combined_eventslist_name}")
    transform_table = np.array(transform_table)
    np.savetxt(transformations_filename, transform_table, fmt="%s")
    median_MJDs = np.array([get_MJD_start_time(path) for path in transform_table[:, 0]])
    plot_transform_table(transform_table, transformations_filename, median_MJDs)


def make_combined_eventslists(channel, ref_episode=None):
    files_list = glob(f"uvt_*/{channel[0]}*/*l2ce.fits")
    if len(files_list) > 0:
        window_filter_list = [file[42:48] for file in files_list]
        window_filter_list = list(set(window_filter_list))
        for window_filter in window_filter_list:
            print(f"\nCombining {channel} {window_filter} files:")
            transformations_filename = f"{channel}_{window_filter}_transformations.txt"
            combined_vis_filename = f"combined_VIS_for_{channel}_{window_filter}.fits"
            combined_eventslist_name = (
                f"combined_{channel}_{window_filter}_events_list.fits"
            )
            chosen_window_filter_expmap_paths = glob(
                f"uvt_*/{channel[0]}*/*{window_filter}*exp.fits*"
            )
            exp_times = (
                get_exptime(expmap) for expmap in chosen_window_filter_expmap_paths
            )
            exp_times = np.array(list(exp_times))
            sort_order = np.argsort(exp_times)
            sort_order = np.flip(sort_order)
            chosen_window_filter_expmap_paths = np.array(
                chosen_window_filter_expmap_paths
            )[sort_order]
            chosen_window_filter_paths = (
                expmap.replace("I_l2exp", "_l2ce")
                for expmap in chosen_window_filter_expmap_paths
            )
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
            VIS_images = [
                glob(ntpath.split(d)[0].replace(f"{channel[0]}", "V") + "/*img.fits")[0]
                for d in chosen_window_filter_paths
            ]
            detected_sources = VIS_detect_sources_in_range(VIS_images)
            trans_matrices_for_VIS = match_frames_joblib_loop(detected_sources)

            # To mask failed cases.
            trans_matrices_for_VIS = np.array(trans_matrices_for_VIS)
            mask = np.logical_not(np.isnan(trans_matrices_for_VIS[:, 0, 0]))
            trans_matrices_for_VIS = trans_matrices_for_VIS[mask]
            detected_sources = list(compress(detected_sources, mask))

            if False in mask:
                print("\nFailed to find VIS transformations for:")
                print(np.array(chosen_window_filter_paths)[~mask])

            VIS_images = list(compress(VIS_images, mask))
            chosen_window_filter_paths = list(
                compress(chosen_window_filter_paths, mask)
            )

            detected_sources_centre_changed = [d - 900 for d in detected_sources]
            trans_matrices = match_frames_joblib_loop(detected_sources_centre_changed)

            with np.printoptions(suppress=True, formatter={"float": "{:.5f}".format}):
                print("\nVIS channel derived transformation matrices:")
                print(np.array(trans_matrices))
            opt_make_coadded_vis_image(
                VIS_images, trans_matrices_for_VIS, combined_vis_filename
            )
            combine_event_lists(
                chosen_window_filter_paths,
                trans_matrices,
                channel,
                combined_eventslist_name,
                transformations_filename,
            )
        print("\nProcessing completed")


warnings.simplefilter("ignore", category=AstropyWarning)

VIS_to_FUV = np.array([[0.85093, -0.56645, 0], [0.56645, 0.85093, 0], [0, 0, 1]])
VIS_to_NUV = np.array([[-1.01652, 0.04649, 0], [-0.04649, -1.01652, 0], [0, 0, 1]])
VIS_to_UV_matrix = {"FUV": VIS_to_FUV, "NUV": VIS_to_NUV}

make_combined_eventslists("NUV")
make_combined_eventslists("FUV")
