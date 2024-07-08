#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


from glob import glob
from astropy.wcs import WCS
from astropy.io import fits
from scipy.optimize import least_squares
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist


def get_average_plate_scale(channel, window_filter):
    XP_FUV = 0.4168705
    YP_FUV = 0.416848964

    XP_NUV = 0.417310
    YP_NUV = 0.417218

    XP_NUV_F2 = 0.419686
    YP_NUV_F2 = 0.419184

    XP_NUV_F1 = 0.4173677
    YP_NUV_F1 = 0.417299

    if channel == "FUV":
        plate_scales = [XP_FUV, YP_FUV]
    elif channel == "NUV":
        if window_filter[-2:] == "F2":
            plate_scales = [XP_NUV_F2, YP_NUV_F2]
        elif window_filter[-2:] == "F1":
            plate_scales = [XP_NUV_F1, YP_NUV_F1]
        else:
            plate_scales = [XP_NUV, YP_NUV]
    else:
        raise ValueError

    average_plate_scale = np.average(plate_scales)
    return average_plate_scale


# Equal plate scales along X and Y, and only rotation are assumed.
def pixel_to_skycoo(px, py, theta, pixel_scale, CRPIX1, CRPIX2, CRVAL1, CRVAL2):
    # For conversion of (px, py) to (x, y) in projection plane coordinates.
    CDELT1 = pixel_scale
    CDELT2 = pixel_scale

    PC1_1 = np.cos(theta)
    PC1_2 = -1 * np.sin(theta)
    PC2_1 = np.sin(theta)
    PC2_2 = np.cos(theta)

    x = CDELT1 * ((PC1_1 * (px - CRPIX1)) + (PC1_2 * (py - CRPIX2)))
    y = CDELT2 * ((PC2_1 * (px - CRPIX1)) + (PC2_2 * (py - CRPIX2)))

    # For conversion of (x, y) to (phi, theta) in native spherical coordinates.
    phi = np.arctan2(x, -1 * y)
    R_theta = (x**2 + y**2) ** 0.5
    theta = np.arctan(180 / (np.pi * R_theta))

    # For conversion of (phi, theta) to (alpha, delta) in celestial spherical coordinates.
    if CRVAL2 == 90:
        phi_p = 0
    else:
        phi_p = np.deg2rad(180)

    alpha_p = np.deg2rad(CRVAL1)
    delta_p = np.deg2rad(CRVAL2)

    alpha = alpha_p + np.arctan2(
        -1 * np.cos(theta) * np.sin(phi - phi_p),
        (np.sin(theta) * np.cos(delta_p))
        - (np.cos(theta) * np.sin(delta_p) * np.cos(phi - phi_p)),
    )
    delta = np.arcsin(
        np.sin(theta) * np.sin(delta_p)
        + np.cos(theta) * np.cos(delta_p) * np.cos(phi - phi_p)
    )

    RA = np.rad2deg(alpha)
    DEC = np.rad2deg(delta)

    return RA, DEC


def get_residuals(
    params, lon, lat, px, py, pixel_scale, CRPIX1, CRPIX2, pixel_to_skycoo
):
    theta = params[0]
    CRVAL1 = params[1]
    CRVAL2 = params[2]

    lon2, lat2 = pixel_to_skycoo(
        px, py, theta, pixel_scale, CRPIX1, CRPIX2, CRVAL1, CRVAL2
    )

    lat_resids = lat - lat2
    lon_resids = lon - lon2
    # In case the longitude has wrapped around
    lon_resids = np.mod(lon_resids - 180.0, 360.0) - 180.0

    resids = np.concatenate((lon_resids * np.cos(np.radians(lat)), lat_resids))

    return resids


def make_error_vector_plot(
    index_xx, index_yy, field_xx, field_yy, num_sources, rms_sky_difference, filename
):
    fig, ax = plt.subplots(figsize=(6, 6))
    for index_x, index_y, field_x, field_y in zip(
        index_xx, index_yy, field_xx, field_yy
    ):
        dummy = ax.arrow(
            field_x,
            field_y,
            (index_x - field_x) * 20,
            (index_y - field_y) * 20,
            width=5,
            head_width=20,
        )

    # Add graphical legend indicating arrow length is multiplied by 20
    ax.legend([dummy], ["Arrow length multiplied by 20"], loc="upper left")

    ax.set_xbound(0, 4800)
    ax.set_ybound(0, 4800)

    ax.set_title(
        f"RMS difference = {rms_sky_difference:.3} arcsec, Number of sources = {num_sources}",
        wrap=True,
    )

    fig.savefig(filename, dpi=200, bbox_inches="tight")


def get_fixed_plate_scale_wcs(filename, channel, window_filter):
    hdu = fits.open(filename)

    field_x = hdu[1].data["field_x"]
    field_y = hdu[1].data["field_y"]
    index_x = hdu[1].data["index_x"]
    index_y = hdu[1].data["index_y"]
    index_ra = hdu[1].data["index_ra"]
    index_dec = hdu[1].data["index_dec"]

    threshold = 10
    if channel == "NUV" and window_filter[-2:] == "F2":
        threshold = 20

    mask = (field_x - index_x) ** 2 + (field_y - index_y) ** 2 <= threshold

    if np.sum(mask) >= 5:
        positions = np.transpose((field_x, field_y))
        positions = positions[mask]

        hull = ConvexHull(positions)
        hullpoints = positions[hull.vertices, :]
        hdist = cdist(hullpoints, hullpoints, metric="euclidean")
        convex_hull_max_distance = hdist.max()

        if convex_hull_max_distance > 1000 and np.sum(mask) != len(mask):
            print("Pruning the list of detections:")
            print(f"Threshold = {threshold}")
            print(f"Number of sources in the original matched list = {len(mask)}")
            print(f"Number of sources in the pruned matched list = {np.sum(mask)}")
        else:
            print("No pruning has been carried out.")
            mask = np.ones(len(mask)).astype(np.bool_)
    else:
        print("No pruning has been carried out.")
        mask = np.ones(len(mask)).astype(np.bool_)

    # FIXED values.
    pixel_scale = get_average_plate_scale(channel, window_filter) / 3600
    CRPIX1 = 2400.5
    CRPIX2 = 2400.5

    params = [0, np.mean(index_ra), np.mean(index_dec)]
    fit = least_squares(
        get_residuals,
        params,
        args=(
            index_ra[mask],
            index_dec[mask],
            field_x[mask],
            field_y[mask],
            pixel_scale,
            CRPIX1,
            CRPIX2,
            pixel_to_skycoo,
        ),
        bounds=[
            [-2 * np.pi, np.mean(index_ra) - 1, np.mean(index_dec) - 1],
            [2 * np.pi, np.mean(index_ra) + 1, np.mean(index_dec) + 1],
        ],
    )

    theta = fit.x[0]

    CRVAL1, CRVAL2 = fit.x[1:3]

    CDELT1 = pixel_scale
    CDELT2 = pixel_scale

    PC1_1 = np.cos(theta)
    PC1_2 = -1 * np.sin(theta)
    PC2_1 = np.sin(theta)
    PC2_2 = np.cos(theta)

    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crpix = [CRPIX1, CRPIX2]
    w.wcs.crval = [CRVAL1, CRVAL2]
    w.wcs.cdelt = [CDELT1, CDELT2]
    w.wcs.pc = np.array([[PC1_1, PC1_2], [PC2_1, PC2_2]])

    fixed_platescale_wcs_hdu = fits.PrimaryHDU(header=w.to_header())
    fixed_platescale_wcs_filename = filename.replace(".corr", "_fixed_platescale.wcs")

    index_xx, index_yy = w.wcs_world2pix(index_ra, index_dec, 1)
    field_ra, field_dec = w.wcs_pix2world(field_x, field_y, 1)

    df = pd.DataFrame(hdu[1].data[mask])
    df = df.iloc[:, 0:8]
    df["index_x"] = index_xx[mask]
    df["index_y"] = index_yy[mask]
    df["field_ra"] = field_ra[mask]
    df["field_dec"] = field_dec[mask]
    df["x_difference"] = df["index_x"] - df["field_x"]
    df["y_difference"] = df["index_y"] - df["field_y"]
    df = df.round(6)
    output_filename = filename.replace("corr", "xlsx")
    writer = pd.ExcelWriter(output_filename)
    df.to_excel(writer, sheet_name="sheet1", index=False, na_rep="NaN")

    for column in df:
        column_length = max(df[column].astype(str).map(len).max(), len(column))
        col_idx = df.columns.get_loc(column)
        writer.sheets["sheet1"].set_column(col_idx, col_idx, column_length + 1)

    writer.close()

    num_sources = len(index_xx[mask])
    difference = np.sqrt(df["x_difference"] ** 2 + df["y_difference"] ** 2)
    rms_difference = np.sqrt(np.mean(np.square(difference)))
    rms_sky_difference = rms_difference * get_average_plate_scale(
        channel, window_filter
    )

    fixed_platescale_wcs_hdu.header["N_SOURCE"] = num_sources
    fixed_platescale_wcs_hdu.header["RMSSKY_D"] = np.round(rms_sky_difference, 3)
    fixed_platescale_wcs_hdu.writeto(fixed_platescale_wcs_filename, overwrite=True)

    make_error_vector_plot(
        index_xx[mask],
        index_yy[mask],
        field_x[mask],
        field_y[mask],
        num_sources,
        rms_sky_difference,
        filename.replace(".corr", "_error_vectors.png"),
    )


corr_files = glob("combined_[F,N]UV*corr")
for corr_file in corr_files:
    print(f"\nWorking on {corr_file}")
    channel, window_filter = corr_file.split("_")[1:3]
    get_fixed_plate_scale_wcs(corr_file, channel, window_filter)

corr_files = glob("UV_oriented_combined_VIS_for_[F,N]UV*corr")
for corr_file in corr_files:
    print(f"\nWorking on {corr_file}")
    channel, window_filter = corr_file.split("_")[-3:-1]
    get_fixed_plate_scale_wcs(corr_file, channel, window_filter)
