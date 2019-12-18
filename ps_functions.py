# -*- coding: utf-8 -*-
"""
A number of functions I keep reusing throughout code.

@author: Tom Williams
"""

import copy

import numpy as np
from pafit.fit_kinematic_pa import fit_kinematic_pa as pafit
from scipy.odr import Model, ODR, RealData


def star_mask(vel):
    # Outlier velocities are almost certainly artifacts so throw away anything  with velocities higher/lower than
    # (-)300km/s

    star_mask_idx = np.where((vel > 300) | (vel < -300))

    return star_mask_idx


def pafit_wrapper(vel, vel_err, x_cen, y_cen):
    # PAfit works by providing a list of x- and y-positions, and the velocity (with systemic subtracted) at that
    # position. IMPORTANT: The *centre* must be at (0,0).

    # Create coordinates from a meshgrid

    x_pafit, y_pafit = np.meshgrid(range(vel.shape[1]), range(vel.shape[0]))

    # Shift x- and y- coords to (0,0) at the centre

    x_pafit -= int(x_cen)
    y_pafit -= int(y_cen)

    # Now flatten everything out and get rid of the NaNs

    idx = np.where(np.isnan(vel) == False)

    x_pafit = x_pafit[idx].flatten()

    # For speed, take 1% of the data.

    x_pafit = x_pafit[::10]
    y_pafit = y_pafit[idx].flatten()[::10]
    vel_pafit = vel[idx].flatten()[::10]
    vel_err_pafit = vel_err[idx].flatten()[::10]

    # Put this into pafit to get any residual systemic velocity and also the kinematic PA

    pafit_pa, pafit_pa_err, v_sys = pafit(x_pafit, y_pafit, vel_pafit, dvel=vel_err_pafit,
                                          quiet=True, plot=False)

    return pafit_pa, pafit_pa_err, v_sys


def get_sample_table_info(galaxy, galaxy_table):
    # From the sample table overview, pull out the RA/Dec, inclination and distance.

    galaxy_found = False

    galaxy_homogenised = copy.copy(galaxy)

    # If the leading 0 after NGC is missing (only actually a problem for NGC628), force it in.

    if galaxy_homogenised[:3] == 'NGC' and len(galaxy_homogenised) == 6:
        galaxy_homogenised = galaxy_homogenised[:3] + '0' + galaxy_homogenised[3:]

    for i in range(len(galaxy_table)):
        if galaxy_table['NAME'][i].strip() == galaxy_homogenised:
            ra, dec, inclination, dist, pa, pa_err = galaxy_table['ORIENT_RA'][i], galaxy_table['ORIENT_DEC'][i], \
                                                     galaxy_table['ORIENT_INCL'][i], galaxy_table['DIST'][i], \
                                                     galaxy_table['ORIENT_POSANG'][i], galaxy_table['ORIENT_POSANG_UNC'][
                                                         i]
            galaxy_found = True
            break

    if not galaxy_found:
        raise Warning(galaxy + ' not found!')

    return ra, dec, inclination, dist, pa, pa_err


def bar_mask(x_grid, y_grid, pa, bar_r):
    # Apply rotation matrix

    rot_matrix = np.matrix([[np.cos(pa * np.pi / 180), np.sin(pa * np.pi / 180)],
                            [-np.sin(pa * np.pi / 180), np.cos(pa * np.pi / 180)]])

    x_rot, y_rot = np.asarray(
        rot_matrix * np.vstack([x_grid.flatten(), y_grid.flatten()]))

    x_rot = x_rot.reshape(x_grid.shape)
    y_rot = y_rot.reshape(y_grid.shape)

    # Find the coordinates where we're greater than the bar and mask them in everything

    idx = np.where(np.abs(x_rot) > 1.2 * bar_r)  # | (np.abs(y_rot) > 1.2*bar_r) )

    return idx


def symmetrize_tw_integral(flux, x_lon, y_lon):
    flux_flat = flux.flatten()
    flux_len = np.arange(len(flux_flat))

    for y_lon_unique in np.unique(y_lon):

        # Take into account any NaNs in the image

        nan_idx = np.where(np.isnan(flux_flat) == False)

        y_lon_idx = np.where(y_lon[nan_idx] == y_lon_unique)

        if len(y_lon_idx[0]) == 0:
            continue

        # Find the most negative and positive x_lon at these positions.

        x_lon_slice = x_lon[nan_idx][y_lon_idx]

        x_lon_min = np.min(x_lon_slice)
        x_lon_max = np.max(x_lon_slice)

        # For things where we only have positive or negative data, remove everything

        if x_lon_min * x_lon_max > 0:

            x_lon_thresh = 0

        else:

            x_lon_thresh = np.min([np.abs(x_lon_max),
                                   np.abs(x_lon_min)])

        x_lon_idx = np.where(np.abs(x_lon[nan_idx][y_lon_idx]) > x_lon_thresh)

        if len(x_lon_idx[0]) == 0:
            continue

        final_idxs = flux_len[nan_idx][y_lon_idx][x_lon_idx]

        for final_idx in final_idxs:
            flux_flat[final_idx] = np.nan

    # Reshape the array and return

    flux_final = flux_flat.reshape(flux.shape)

    return flux_final


def linear_fit(theta, x):
    return theta[0] * x + theta[1]


def odr_fit(x, x_err, y, y_err):
    linear = Model(linear_fit)

    # Filter out any NaNs else ODR breaks.

    nan_idx = np.where(np.isnan(x) == False)

    # Fit using SciPy's ODR

    odr_data = RealData(x[nan_idx], y[nan_idx],
                        sx=x_err[nan_idx], sy=y_err[nan_idx])

    odr_run = ODR(odr_data, linear, beta0=[2, 0])

    odr_output = odr_run.run()

    m, c = odr_output.beta
    m_err, c_err = odr_output.sd_beta

    return m, m_err, c, c_err
