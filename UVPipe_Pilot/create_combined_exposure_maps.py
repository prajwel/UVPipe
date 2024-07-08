#!/usr/bin/env python3

import numpy as np

from glob import glob
from astropy.io import fits
from joblib import Parallel, delayed
from astropy.convolution import convolve


def get_framedata(file):
    hdu = fits.open(file)
    return hdu[0].data


def get_framedata_joblib_loop(fits_files):
    all_frames = Parallel(n_jobs=-1)(delayed(get_framedata)(i) for i in fits_files)
    return all_frames


def back_sampling(x_prime, y_prime, transform):
    _, cos_rot, sin_rot, shift_x, shift_y = transform
    x_prime = x_prime - 2400
    y_prime = y_prime - 2400
    x = (x_prime - shift_x) * cos_rot + (y_prime - shift_y) * sin_rot
    y = -(x_prime - shift_x) * sin_rot + (y_prime - shift_y) * cos_rot
    x = x + 2400
    y = y + 2400
    x = np.round(x).astype(int)
    y = np.round(y).astype(int)
    return x, y


xx, yy = np.meshgrid(np.arange(0, 4800), np.arange(0, 4800))


def transform_single_expmap(expmap, transform, xx=xx, yy=yy):
    frame_data = np.array([xx.flatten(), yy.flatten()]).T
    xx_back, yy_back = back_sampling(frame_data[:, 0], frame_data[:, 1], transform)

    zero_padding = np.zeros(len(xx_back)).astype(int)
    zero_padded_data = np.array([xx_back, yy_back, zero_padding]).T
    frame_data = np.hstack([frame_data, zero_padded_data])

    mask_x = np.logical_and(frame_data[:, 2] >= 0, frame_data[:, 2] <= 4799)
    mask_y = np.logical_and(frame_data[:, 3] >= 0, frame_data[:, 3] <= 4799)
    mask = np.logical_and(mask_x, mask_y)

    back_sampled_data = expmap[frame_data[:, 3][mask], frame_data[:, 2][mask]]
    frame_data[:, 4][mask] = back_sampled_data
    frame_data = frame_data[:, (0, 1, 4)]

    return frame_data


def transform_expmaps_joblib_loop(expmaps, transformations):
    all_frame_data = Parallel(n_jobs=-1, max_nbytes=None)(
        delayed(transform_single_expmap)(expmap, transformation)
        for expmap, transformation in zip(expmaps, transformations)
    )

    return all_frame_data


def combine_expmaps(transformations_filename, combined_expmap_filename):
    transformations = np.genfromtxt(transformations_filename, dtype=None, encoding=None)
    transformations = np.atleast_1d(transformations)
    all_expmap_files = [glob(d[0] + "/*I_l2exp*fits")[0] for d in transformations]
    print(np.array(all_expmap_files))

    expmaps = get_framedata_joblib_loop(all_expmap_files)
    expmaps = np.array(expmaps)

    # SNT kernel
    snt_kernel = np.ones([5, 5])
    snt_kernel[:] = 0.4
    snt_kernel[1:4, 1:4] = 0.6
    snt_kernel[2, 2] = 1

    final_expmap_frame_data = transform_expmaps_joblib_loop(expmaps, transformations)
    final_expmap_frame_data = np.vstack(final_expmap_frame_data)
    expmaps = None
    mask = final_expmap_frame_data[:, 2] != 0
    final_expmap_frame_data = final_expmap_frame_data[mask]

    bins = np.arange(-0.5, 4800.5)
    ndarray, yedges, xedges = np.histogram2d(
        final_expmap_frame_data[:, 1],
        final_expmap_frame_data[:, 0],
        bins=(bins, bins),
        weights=final_expmap_frame_data[:, 2],
    )

    final_expmap_frame_data = None

    filtered_ndarray = convolve(ndarray, snt_kernel, normalize_kernel=True)
    hdu = fits.PrimaryHDU(data=filtered_ndarray.astype(np.int64))
    hdu.writeto(combined_expmap_filename, overwrite=True)
    print(f"\nCombined exposure map: {combined_expmap_filename}")


def make_combined_expmaps():
    filelist = glob("[F,N]UV*_transformations.txt")
    if len(filelist) > 0:
        for file in filelist:
            channel = file[:3]
            window_filter = file[4:10]
            print(f"\nCombining {channel} {window_filter} exposure maps:")
            combined_expmap_filename = (
                f"combined_{channel}_{window_filter}_exp_map.fits"
            )
            combine_expmaps(file, combined_expmap_filename)
    print("\nProcessing completed")


make_combined_expmaps()
