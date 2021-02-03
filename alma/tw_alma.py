# -*- coding: utf-8 -*-
"""
Apply Tremaine-Weinberg method to all of the PHANGS ALMA galaxies

@author: Tom Williams
"""

import os
import sys
import time
import warnings
from collections import defaultdict

import cmocean
import corner
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pybar.pybar as pybar
from astropy.io import fits
from astropy.wcs import WCS

sys.path.append(os.getcwd())

import ps_functions
from bootstraps import bootstrap_tw, mcmc_tw_multiple
from vars import phangs_folder, galaxy_table, alma_version, alma_plot, alma_output, mask_outside_bars, slit_widths, \
    slit_lengths, use_sim, n_pattern_speeds, alma_galaxies, overwrite_pafit, plot, use_gmc_mask, \
    gmc_removal_factors

warnings.simplefilter("ignore")

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

# Set up random seed for testing.

np.random.seed(42)

start = time.time()

os.chdir(phangs_folder)

if not os.path.exists(alma_plot + alma_version):
    os.makedirs(alma_plot + alma_version)
if not os.path.exists(alma_output + alma_version):
    os.makedirs(alma_output + alma_version)

use_strict = False
use_lowres = False
use_lang_rot = False
null_hypothesis = True
correct_for_null = False

correction_type = 'strictest'

overwrite_bootstraps = False

null_hypothesis_sigma = 3

# And test galaxies

# alma_galaxies = ['ngc1546']

for galaxy in alma_galaxies:

    galaxy = galaxy.upper()

    print('Beginning ' + galaxy)

    for mask_outside_bar in mask_outside_bars:

        for slit_width in slit_widths:

            for slit_length in slit_lengths:

                for gmc_removal_factor in gmc_removal_factors:

                    hdu_file_name = 'alma/' + alma_version + '/' + galaxy + '_mom0'

                    if use_strict:
                        hdu_file_name += '_strict'

                    if use_lowres:
                        hdu_file_name += '_lowres'

                    hdu_file_name += '.fits'

                    try:
                        flux_hdu = fits.open(hdu_file_name)[0]
                    except FileNotFoundError:
                        print(galaxy + ' not found!')
                        continue

                    # Set up file names depending on the various parameters

                    file_name = alma_version + '/'

                    if np.sum([len(slit_widths) > 1, len(slit_lengths) > 1, len(mask_outside_bars) > 1]) > 1:
                        raise Warning('Too many things varying!')

                    if use_gmc_mask:

                        file_name += 'gmc_mask/'

                    elif len(slit_widths) > 1:

                        file_name += 'slit_width/'

                    elif len(slit_lengths) > 1:

                        file_name += 'slit_length/'

                    elif mask_outside_bar:

                        file_name += 'bmask/'

                    elif not mask_outside_bar:

                        file_name += 'nobmask/'

                    else:

                        raise Warning('Not sure what is going on here!')

                    if not os.path.exists(alma_plot + file_name):
                        os.makedirs(alma_plot + file_name)
                    if not os.path.exists(alma_output + file_name):
                        os.makedirs(alma_output + file_name)

                    if n_pattern_speeds > 1:
                        file_name += '%d_pattern_speeds/' % n_pattern_speeds

                    # If we're doing a simulation then that should go in a folder of its own

                    if use_sim:
                        file_name += 'sim/'
                    if use_lang_rot:
                        file_name += 'lang/'
                    if null_hypothesis:
                        file_name += 'null/' + correction_type + '/'
                    if correct_for_null:
                        file_name += 'corr/' + correction_type + '/'

                    if not os.path.exists(alma_plot + file_name):
                        os.makedirs(alma_plot + file_name)
                    if not os.path.exists(alma_output + file_name):
                        os.makedirs(alma_output + file_name)

                    file_name += galaxy

                    # Set up file extensions for mask

                    if mask_outside_bar:
                        file_name += '_bmask_'
                    else:
                        file_name += '_nobmask_'

                    if len(slit_lengths) > 1:

                        file_name += 'sl_%.0f_' % slit_length

                    elif len(slit_widths) > 1:

                        file_name += 'sw_%.1f_' % slit_width

                    elif len(gmc_removal_factors) > 1:

                        file_name += 'gmc_%.1f_' % gmc_removal_factor

                    if use_strict:
                        file_name += 'strict_'

                    if use_lowres:
                        file_name += 'lowres_'

                    if use_lang_rot:
                        file_name += 'lang_'

                    if null_hypothesis:
                        file_name += 'null_'

                    if correct_for_null:
                        file_name += 'corr_%d_sigma_' % null_hypothesis_sigma

                    pattern_speed_plot_filename = alma_plot + file_name + 'alma'
                    pattern_speed_filename = alma_output + file_name + 'pattern_speed_alma.txt'
                    bootstrap_filename = alma_output + file_name + 'bootstrap_alma.txt'
                    max_r_filename = alma_output + file_name + 'max_r.npy'

                    # Remove any extraneous third dimensions from the ALMA data

                    del flux_hdu.header['PC3_1'], flux_hdu.header['PC3_2'], flux_hdu.header['PC1_3'], flux_hdu.header[
                        'PC2_3'], \
                        flux_hdu.header['PC3_3'], flux_hdu.header['CTYPE3'], flux_hdu.header['CRVAL3'], flux_hdu.header[
                        'CDELT3'], \
                        flux_hdu.header['CRPIX3'], flux_hdu.header['CUNIT3']

                    wcs = WCS(flux_hdu)

                    # Read in the various .fits images we'll need.

                    flux = flux_hdu.data
                    flux_err = fits.open(hdu_file_name.replace('mom0', 'emom0'))[0].data

                    vel = fits.open(hdu_file_name.replace('mom0', 'mom1'))[0].data

                    if use_lang_rot:
                        try:
                            vel = fits.open(os.path.join('vrot_losmodels', galaxy.lower() + 'LOSmodel.fits'))[0].data

                            # Make sure  the  shapes match
                            if not np.all(vel.shape == flux.shape):
                                continue
                        except FileNotFoundError:
                            continue

                    vel_err = fits.open(hdu_file_name.replace('mom0', 'emom1'))[0].data

                    # Remove any 0s in the flux image, replace with NaNs

                    idx = np.where(flux == 0)

                    flux[idx] = np.nan
                    flux_err[idx] = np.nan
                    vel[idx] = np.nan
                    vel_err[idx] = np.nan

                    # And also for places where the velocity isn't defined, mask that too

                    idx = np.where(np.isnan(vel))

                    flux[idx] = np.nan
                    flux_err[idx] = np.nan
                    vel[idx] = np.nan
                    vel_err[idx] = np.nan

                    idx = np.where(np.isnan(flux))

                    if np.all(np.isnan(flux)):
                        print('No flux detected. Skipping')
                        continue

                    # Pull out parameters we need from the sample table.

                    sample_table_parameters = ['orient_ra', 'orient_dec',
                                               'orient_ra_unc', 'orient_dec_unc',
                                               'orient_incl', 'dist',
                                               'orient_posang', 'orient_posang_unc',
                                               'morph_bar_r', 'orient_vlsr'
                                               ]

                    parameters_dict = ps_functions.get_sample_table_info(galaxy.lower(), galaxy_table,
                                                                         sample_table_parameters)

                    ra = parameters_dict['orient_ra']
                    dec = parameters_dict['orient_dec']
                    centering_err = (parameters_dict['orient_ra_unc'] + parameters_dict['orient_dec_unc']) / 2
                    inclination = parameters_dict['orient_incl']
                    pa = parameters_dict['orient_posang']
                    pa_err = parameters_dict['orient_posang_unc']
                    dist = parameters_dict['dist']
                    bar_r = parameters_dict['morph_bar_r']

                    # This won't work for an inclination of 90 so if it is, skip

                    if inclination == 90:
                        print('Edge-on: skipping')
                        continue

                    # If the galaxy doesn't have a measured PA, skip

                    if np.isnan(pa):
                        print('PA undefined: skipping')
                        continue

                    # Calculate a physical distance conversion factor from the galaxy distance.

                    kpc_per_arcsec = dist * 1e3 * np.sin(1 / 3600 * np.pi / 180)

                    # Convert the RA and Dec to a central pixel.

                    x_cen, y_cen = wcs.all_world2pix(ra, dec, 1)

                    if not os.path.exists(
                            'pattern_speeds_output/pafit/' + galaxy + '_pafit_output_alma.txt') or overwrite_pafit:

                        pafit_pa, pafit_pa_err, v_sys = ps_functions.pafit_wrapper(vel, vel_err, x_cen, y_cen)

                        np.savetxt('pattern_speeds_output/pafit/' + galaxy + '_pafit_output_alma.txt',
                                   np.c_[pafit_pa, pafit_pa_err, v_sys],
                                   header='kin_pa, kin_pa_err, v_sys')

                    else:

                        pafit_pa, pafit_pa_err, v_sys = np.loadtxt(
                            'pattern_speeds_output/pafit/' + galaxy + '_pafit_output_alma.txt',
                            unpack=True)

                    # print('pafit/Lang fit: %.2f/%.2f' % (pafit_pa, pa))

                    # Subtract any residual velocity

                    v_sys = parameters_dict['orient_vlsr']

                    if not (use_sim or use_lang_rot):
                        vel -= v_sys

                    # Set up a grid of x- and y-coordinates

                    grid_shape = flux.shape
                    grid_step = np.abs(flux_hdu.header['CDELT1'] * 3600)
                    grid_x = (np.arange(grid_shape[1]) - x_cen) * grid_step
                    grid_y = (np.arange(grid_shape[0]) - y_cen) * grid_step

                    x_coords, y_coords = np.meshgrid(grid_x, grid_y)

                    if use_sim:
                        rot_matrix = np.matrix([[np.cos(pa * np.pi / 180), np.sin(pa * np.pi / 180)],
                                                [-np.sin(pa * np.pi / 180), np.cos(pa * np.pi / 180)]])

                        x_rot, y_rot = np.asarray(
                            rot_matrix * np.vstack([x_coords.flatten(), y_coords.flatten()]))

                        x_rot = x_rot.reshape(grid_shape)
                        y_rot = y_rot.reshape(grid_shape)

                        # We want to inject a fake velocity field of 50 km/s/kpc

                        vel = 50 * y_rot * kpc_per_arcsec * np.sin(inclination * np.pi / 180)
                        vel_err = 50 * np.ones(vel.shape)

                        r = np.sqrt(x_rot ** 2 + y_rot ** 2)

                        flux = np.exp(-0.5 * r)
                        flux *= 10 / np.nanmax(flux)

                        # And mask to match the observation FOV

                        flux[idx] = np.nan
                        flux_err[idx] = np.nan
                        vel[idx] = np.nan
                        vel_err[idx] = np.nan

                    if mask_outside_bar:

                        # Mask anything Beyond the Bar in line with the kinematic PA.

                        if bar_r > 0:
                            idx = ps_functions.bar_mask(x_coords, y_coords, pa, bar_r)

                            flux[idx] = np.nan
                            flux_err[idx] = np.nan
                            vel[idx] = np.nan
                            vel_err[idx] = np.nan

                        # else:
                        #
                        #     print('Bar not defined. Skipping masking')

                    if len(slit_lengths) > 1:
                        rot_matrix = np.matrix([[np.cos(pa * np.pi / 180), np.sin(pa * np.pi / 180)],
                                                [-np.sin(pa * np.pi / 180), np.cos(pa * np.pi / 180)]])

                        x_rot, y_rot = np.asarray(
                            rot_matrix * np.vstack([x_coords.flatten(), y_coords.flatten()]))

                        x_rot = x_rot.reshape(grid_shape)
                        y_rot = y_rot.reshape(grid_shape)

                        # Find the coordinates where we're longer than specified slit length

                        idx = np.where((np.abs(x_coords) > slit_length) | (np.abs(y_coords) > slit_length))

                        flux[idx] = np.nan
                        flux_err[idx] = np.nan
                        vel[idx] = np.nan
                        vel_err[idx] = np.nan

                    if len(gmc_removal_factors) > 1:

                        # Read in the projected GMC mask

                        try:
                            antenna_combination = {'NGC0300': '7m+tp'}[galaxy]
                        except KeyError:
                            antenna_combination = '12m+7m+tp'

                        gmc_hdu_name = galaxy.lower() + '_' + antenna_combination + '_co21_native_asgn.fits.bz2'

                        gmc_mask = fits.open(phangs_folder + '/pattern_speeds_output/gmc_catalogues/' + gmc_hdu_name)[0]

                        max_value = np.nanmax(gmc_mask.data)

                        mask = defaultdict(lambda: defaultdict(list))

                        for i in range(gmc_mask.data.shape[1]):
                            for j in range(gmc_mask.data.shape[2]):
                                mask[i][j] = [x for x in gmc_mask.data[:, i, j] if x > 0]

                        masked_pixel_idx = np.random.choice(range(1, int(max_value)),
                                                            size=int(gmc_removal_factor * max_value),
                                                            replace=False)

                        mask_collapse = np.zeros(flux.shape)

                        for i in range(mask_collapse.shape[0]):
                            for j in range(mask_collapse.shape[1]):

                                if any(x in masked_pixel_idx for x in mask[i][j]):
                                    mask_collapse[i, j] = 1

                        flux[mask_collapse == 1] = np.nan
                        flux_err[mask_collapse == 1] = np.nan
                        vel[mask_collapse == 1] = np.nan
                        vel_err[mask_collapse == 1] = np.nan

                    # Finally, make sure the integrals will be symmetrical about the PA

                    bar = pybar.mybar(Flux=flux, Flux_err=flux_err,
                                      Velocity=vel, Velocity_err=vel_err,
                                      Xin=x_coords, Yin=y_coords, PAnodes=pa,
                                      inclin=inclination)

                    y_lon_round = bar.Y_lon.round(0)
                    x_lon = bar.X_lon

                    flux = ps_functions.symmetrize_tw_integral(flux, x_lon, y_lon_round)

                    flux_err[np.isnan(flux)] = np.nan
                    vel[np.isnan(flux)] = np.nan
                    vel_err[np.isnan(flux)] = np.nan

                    # Pull out the maximum slit length and save

                    x_lon_reshape = bar.X_lon.reshape(flux.shape)
                    x_lon_reshape[np.isnan(flux)] = np.nan

                    y_lon_reshape = bar.Y_lon.reshape(flux.shape)
                    y_lon_reshape[np.isnan(flux)] = np.nan

                    r_max = np.nanmax(np.abs(x_lon_reshape))

                    np.save(max_r_filename, r_max)

                    if null_hypothesis or correct_for_null:

                        r = np.sqrt(x_lon_reshape ** 2 + (y_lon_reshape / np.cos(np.radians(inclination))) ** 2)

                        r_edges = np.linspace(0, np.nanmax(r), 100)
                        spacing = (r_edges[1] - r_edges[0]) / 2

                        vel_null_hdu_name = os.path.join('fullprojVrot', galaxy + 'LOSmodel')
                        if use_lowres:
                            vel_null_hdu_name += '7m'
                        vel_null_hdu_name += 'FULL.fits'

                        if not os.path.exists(vel_null_hdu_name):

                            # Create a low-budget symmetrised velocity field

                            vel_null = np.zeros_like(vel)
                            vel_null[vel_null == 0] = np.nan

                            for r_edge in r_edges:
                                idx = np.where((r >= r_edge - spacing) & (r <= r_edge + spacing))

                                vel_null[idx] = np.nanmedian(np.abs(vel[idx])) * - np.sign(x_lon_reshape[idx])

                        else:
                            vel_null = fits.open(vel_null_hdu_name)[0].data

                        # Make sure  the shapes match
                        if not np.all(vel_null.shape == flux.shape):
                            continue

                        coverage = np.ones_like(flux)
                        coverage[(flux == 0) | (np.isnan(flux)) | (vel_null == 0)] = 0

                        # Create a low-budget symmetrised flux field

                        # flux_null = np.zeros_like(flux).flatten()
                        # flux_flat = flux.flatten()
                        #
                        # for y_lon_unique in np.unique(y_lon_null_round):
                        #
                        #     idx = np.where(y_lon_null_round == y_lon_unique)[0]
                        #     # symmetrised_val = np.nanmean(flux_flat[idx])
                        #     # flux_null[idx] = symmetrised_val
                        #     for i in range(len(idx)):
                        #         symmetrised_val = (np.abs(flux_flat[idx[i]]) + np.abs(flux_flat[idx[-i]])) / 2
                        #
                        #         flux_null[idx[i]] = symmetrised_val
                        #
                        # flux_null = flux_null.reshape(flux.shape)

                        # flux_null = np.ones_like(vel.data)

                        if correction_type == 'moderate':

                            flux_null = np.zeros_like(flux)
                            flux_null[flux_null == 0] = np.nan

                            for r_edge in r_edges:

                                idx = np.where((r >= r_edge - spacing) & (r <= r_edge + spacing))

                                flux_null[idx] = np.nanmedian(flux[idx])

                        else:

                            flux_null = flux.copy()

                        if correction_type == 'strictest':

                            vel_null = vel.copy()

                        flux_null[coverage == 0] = np.nan
                        vel_null[coverage == 0] = np.nan

                        if null_hypothesis:
                            flux = flux_null
                            flux_err = flux_err
                            vel = vel_null
                            vel_err = vel_err

                    # Pull out the average fit.

                    bar = pybar.mybar(Flux=flux, Flux_err=flux_err,
                                      Velocity=vel, Velocity_err=vel_err,
                                      Xin=x_coords, Yin=y_coords, PAnodes=pa,
                                      inclin=inclination)

                    bar.tremaine_weinberg(slit_width=slit_width)

                    x_tw = bar.dfx_tw
                    v_tw = bar.dfV_tw

                    x_tw_err = bar.dfx_tw_err
                    v_tw_err = bar.dfV_tw_err

                    if correct_for_null:

                        if correction_type is not 'strictest':

                            bar_null = pybar.mybar(Flux=flux_null, Flux_err=np.ones_like(flux_null),
                                                   Velocity=vel_null, Velocity_err=np.ones_like(flux_null),
                                                   Xin=x_coords, Yin=y_coords, PAnodes=pa,
                                                   inclin=inclination)
                            bar_null.tremaine_weinberg(slit_width=slit_width)

                            x_tw_null = bar_null.dfx_tw
                            v_tw_null = bar_null.dfV_tw

                        else:

                            x_tw_null = np.zeros_like(x_tw)
                            v_tw_null = np.zeros_like(v_tw)

                            dig = bar.dig.reshape(flux.shape)
                            x_lon = bar.X_lon.reshape(flux.shape)

                            for i in np.unique(dig):
                                idx = np.where(dig == i)

                                flux_dig = flux[idx]
                                vel_dig = vel[idx]
                                x_lon_dig = x_lon[idx]

                                nan_idx = np.where(np.isnan(flux_dig))[0]
                                anti_nan_idx = -1 - nan_idx

                                flux_anti_nan = np.nansum(flux_dig[anti_nan_idx])
                                v_tw_null[i] = np.nansum(vel_dig[anti_nan_idx] * flux_dig[anti_nan_idx]) / \
                                               flux_anti_nan
                                x_tw_null[i] = np.nansum(x_lon_dig[anti_nan_idx] * flux_dig[anti_nan_idx]) / \
                                               flux_anti_nan

                        # NaN out anything that isn't above the null

                        mask = np.where((np.abs(x_tw) - np.abs(x_tw_null) <= null_hypothesis_sigma * x_tw_err) |
                                        (np.abs(v_tw) - np.abs(v_tw_null) <= null_hypothesis_sigma * v_tw_err))

                        # NaN out anything that isn't significantly different in <x> or <v>

                        # mask = np.where((np.abs(x_tw_null - x_tw) <= null_hypothesis_sigma * x_tw_err) |
                        #                 (np.abs(v_tw_null - v_tw) <= null_hypothesis_sigma * v_tw_err))

                        x_tw[mask] = np.nan
                        v_tw[mask] = np.nan

                        x_tw_err[mask] = np.nan
                        v_tw_err[mask] = np.nan

                        # x_tw_err = x_tw_null  # np.sqrt(x_tw_err ** 2 + x_tw_null ** 2)
                        # v_tw_err = v_tw_null  # np.sqrt(v_tw_err ** 2 + v_tw_null ** 2)

                    # Because some position angles are from HYPERLEDA, and only go from 0 to 180, sometimes they have
                    # the wrong sense. Check that here, and if this is the case swap around.

                    m, m_err, c, c_err = ps_functions.odr_fit(x_tw, x_tw_err, v_tw, v_tw_err)

                    if m < 0:
                        pa += 180

                        bar = pybar.mybar(Flux=flux, Flux_err=flux_err,
                                          Velocity=vel, Velocity_err=vel_err,
                                          Xin=x_coords, Yin=y_coords, PAnodes=pa,
                                          inclin=inclination)

                        bar.tremaine_weinberg(slit_width=slit_width)

                        x_tw = bar.dfx_tw
                        v_tw = bar.dfV_tw

                        x_tw_err = bar.dfx_tw_err
                        v_tw_err = bar.dfV_tw_err

                        if correct_for_null:

                            if correction_type is not 'strictest':
                                bar_null = pybar.mybar(Flux=flux_null, Flux_err=np.ones_like(flux_null),
                                                       Velocity=vel_null, Velocity_err=np.ones_like(flux_null),
                                                       Xin=x_coords, Yin=y_coords, PAnodes=pa,
                                                       inclin=inclination)
                                bar_null.tremaine_weinberg(slit_width=slit_width)

                                x_tw_null = bar_null.dfx_tw
                                v_tw_null = bar_null.dfV_tw

                            else:

                                x_tw_null = np.zeros_like(x_tw)
                                v_tw_null = np.zeros_like(v_tw)

                                dig = bar.dig.reshape(flux.shape)
                                x_lon = bar.X_lon.reshape(flux.shape)

                                for i in np.unique(dig):
                                    idx = np.where(dig == i)

                                    flux_dig = flux[idx]
                                    vel_dig = vel[idx]
                                    x_lon_dig = x_lon[idx]

                                    nan_idx = np.where(np.isnan(flux_dig))[0]
                                    anti_nan_idx = -1 - nan_idx

                                    flux_anti_nan = np.nansum(flux_dig[anti_nan_idx])
                                    v_tw_null[i] = np.nansum(vel_dig[anti_nan_idx] * flux_dig[anti_nan_idx]) / \
                                                   flux_anti_nan
                                    x_tw_null[i] = np.nansum(x_lon_dig[anti_nan_idx] * flux_dig[anti_nan_idx]) / \
                                                   flux_anti_nan

                            # NaN out anything that isn't above the null

                            mask = np.where((np.abs(x_tw) - np.abs(x_tw_null) <= null_hypothesis_sigma * x_tw_err) |
                                            (np.abs(v_tw) - np.abs(v_tw_null) <= null_hypothesis_sigma * v_tw_err))

                            # NaN out anything that isn't significantly different in <x> or <v>

                            # mask = np.where((np.abs(x_tw_null - x_tw) <= null_hypothesis_sigma * x_tw_err) |
                            #                 (np.abs(v_tw_null - v_tw) <= null_hypothesis_sigma * v_tw_err))

                            x_tw[mask] = np.nan
                            v_tw[mask] = np.nan

                            x_tw_err[mask] = np.nan
                            v_tw_err[mask] = np.nan

                        m, m_err, c, c_err = ps_functions.odr_fit(x_tw, x_tw_err, v_tw, v_tw_err)

                    # Convert fitted numbers to something usefully physical -- km/s/kpc, and get rid of the
                    # inclination

                    km_s_kpc_conversion = (np.sin(inclination * np.pi / 180) * kpc_per_arcsec) ** -1

                    omega_bar = m.copy()
                    omega_bar_err = np.abs(m_err.copy())

                    omega_bar *= km_s_kpc_conversion
                    omega_bar_err *= km_s_kpc_conversion

                    if plot and n_pattern_speeds == 1:

                        # Plot on a best fit line and associated errors for this single bootstrap

                        x_max = 1.2 * np.nanpercentile(x_tw, 99)
                        x_min = 1.2 * np.nanpercentile(x_tw, 1)

                        # If the integrals are messed up, skip.

                        if np.any(np.isnan([x_min, x_max])):
                            print('Integral issues here')
                            continue

                        # x_min, x_max = 1.2 * np.nanmin(x_tw), 1.2 * np.nanmax(x_tw)

                        x_fit = np.linspace(x_min, x_max, 500)

                        y_max = 1.2 * np.nanpercentile(v_tw, 99)
                        y_min = 1.2 * np.nanpercentile(v_tw, 1)

                        # y_min, y_max = 1.2 * np.nanmin(v_tw), 1.7 * np.nanmax(v_tw)

                        n_lines = 200

                        y_samples = np.zeros([n_lines, len(x_fit)])

                        for i in range(n_lines):
                            y_samples[i, :] = np.random.normal(loc=m, scale=np.abs(m_err)) * x_fit + \
                                              np.random.normal(loc=c, scale=np.abs(c_err))

                        y_upper = np.percentile(y_samples, 84, axis=0)
                        y_lower = np.percentile(y_samples, 16, axis=0)
                        y_median = np.median(y_samples, axis=0)

                        # Calculate the position of the slits

                        slit_positions = bar.y_slits / grid_step

                        # Calculate the axis line

                        rad = np.sqrt((flux.shape[0] - y_cen) ** 2 + (flux.shape[1] - x_cen) ** 2)
                        ang = [0, np.pi] + np.radians(pa)

                        # Plotting

                        # Set up colourbar limits for the image

                        vmin_flux = np.nanpercentile(flux, 0.5)
                        vmax_flux = np.nanpercentile(flux, 99.75)

                        vmax_vel = 150
                        vmin_vel = -150

                        plt.figure(figsize=(8, 9))

                        plt.subplot(221, projection=wcs)

                        plt.imshow(flux,
                                   cmap=cmocean.cm.gray_r,
                                   origin='lower', interpolation='none',
                                   vmin=vmin_flux, vmax=vmax_flux)

                        colour = iter(cmocean.cm.thermal(np.linspace(0, 1, len(slit_positions))))

                        # Plot only 25% of the slits

                        n_slits = 4

                        # Pull out the LoN matrix.

                        lon_matrix = bar._mat_lon * bar._mat_NE

                        for i, slit_position in enumerate(slit_positions):
                            c = next(colour)

                            if i % n_slits == 0:
                                # Use the rotation matrix to transform into the LoN plane.

                                x0, x1 = -flux.shape[1], flux.shape[1]

                                coord_rotate = np.asarray(np.array([x0, slit_position]) * lon_matrix)
                                x0_rotate, y0_rotate = coord_rotate[0][0], coord_rotate[0][1]

                                coord_rotate = np.asarray(np.array([x1, slit_position]) * lon_matrix)
                                x1_rotate, y1_rotate = coord_rotate[0][0], coord_rotate[0][1]

                                plt.plot([x0_rotate + x_cen, x1_rotate + x_cen],
                                         [y0_rotate + y_cen, y1_rotate + y_cen],
                                         c=c, linewidth=1,
                                         )

                        plt.xlim([0, flux.shape[1]])
                        plt.ylim([0, flux.shape[0]])

                        # ax.coords[0].set_ticks(spacing=5 * u.arcsec)

                        # plt.axis('off')
                        plt.xlabel('RA (J2000)')
                        plt.ylabel('Dec (J2000)')

                        ax = plt.subplot(222, projection=wcs)

                        plt.imshow(vel,
                                   cmap=cmocean.cm.balance,
                                   origin='lower', interpolation='none',
                                   vmin=vmin_vel, vmax=vmax_vel)

                        # Also plot on the axis line

                        plt.plot(-rad * np.sin(ang) + x_cen, rad * np.cos(ang) + y_cen, 'k--', linewidth=2)

                        plt.xlim([0, flux.shape[1]])
                        plt.ylim([0, flux.shape[0]])

                        # plt.axis('off')
                        plt.xlabel('RA (J2000)')
                        plt.ylabel('Dec (J2000)')

                        ax.coords[1].set_ticklabel_position('r')
                        ax.coords[1].set_axislabel_position('r')

                        ax = plt.subplot(212)

                        colour = iter(cmocean.cm.thermal(np.linspace(0, 1, len(x_tw))))

                        for i in range(len(x_tw)):
                            c = next(colour)

                            plt.errorbar(x_tw[i], v_tw[i],
                                         xerr=x_tw_err[i], yerr=v_tw_err[i],
                                         marker='o', c=c, linestyle='none')

                        if correct_for_null:
                            plt.plot(x_tw_null, v_tw_null, 'o', c='r', mfc='none')

                        plt.plot(x_fit, y_median, c='k', lw=2)
                        plt.fill_between(x_fit, y_upper, y_lower,
                                         edgecolor='k', facecolor='k', alpha=0.25)

                        plt.xlim([x_min, x_max])
                        plt.ylim([y_min, y_max])

                        plt.text(0.05, 0.85,
                                 r'$\Omega_p = %.2f \pm %.2f\, \mathrm{km\,s}^{-1} \mathrm{kpc}^{-1}$ ' % (
                                     omega_bar, omega_bar_err),
                                 transform=ax.transAxes,
                                 bbox=dict(facecolor='white', alpha=0.7))

                        plt.xlabel(r'<$x$> $\left(^{\prime \prime}\right)$')
                        plt.ylabel(r'<$v$> $\left(\mathrm{km\,s}^{-1}\right)$')

                        # plt.tight_layout()
                        # plt.show()

                        plt.savefig(pattern_speed_plot_filename + '.png', bbox_inches='tight')
                        plt.savefig(pattern_speed_plot_filename + '.pdf', bbox_inches='tight')
                        plt.close()

                    if use_sim:

                        # We're testing the FOV here, not the position angle. Save out the results of the one fit

                        np.savetxt(pattern_speed_filename, np.c_[omega_bar, omega_bar_err, omega_bar_err],
                                   header='omega_bar, omega_barr_err_up, omega_bar_err_down (all km/s/kpc)')

                    else:

                        if n_pattern_speeds <= 0:

                            raise Warning('This is not a number of pattern speeds I can deal with')

                        elif n_pattern_speeds == 1:

                            # Put everything into the bootstrap

                            if not correct_for_null:
                                flux_null, vel_null = None, None

                            bootstrap_tw(flux, flux_err, vel, vel_err,
                                         flux_null=flux_null, vel_null=vel_null, null_sigma=null_hypothesis_sigma,
                                         grid_step=grid_step, slit_width=slit_width,
                                         x_cen=x_cen, y_cen=y_cen, centering_err=centering_err, pa=pa, pa_err=pa_err,
                                         inclination=inclination,
                                         dist=dist, bootstrap_filename=bootstrap_filename,
                                         overwrite_bootstraps=overwrite_bootstraps,
                                         n_bootstraps=1000, pattern_speed_filename=pattern_speed_filename,
                                         correction_type=correction_type
                                         )

                        elif n_pattern_speeds > 1:

                            bootstrap_filename = bootstrap_filename.replace('bootstrap', 'samples').replace('txt',
                                                                                                            'pkl')

                            n_steps, n_walkers = 1000, 50

                            sampler, pattern_speeds = mcmc_tw_multiple(flux, flux_err, vel, vel_err,
                                                                       n_pattern_speeds=n_pattern_speeds,
                                                                       grid_step=grid_step, slit_width=slit_width,
                                                                       x_cen=x_cen, y_cen=y_cen,
                                                                       centering_err=centering_err,
                                                                       pa=pa, pa_err=pa_err,
                                                                       inclination=inclination, dist=dist,
                                                                       n_steps=n_steps, n_walkers=n_walkers,
                                                                       samples_filename=bootstrap_filename,
                                                                       overwrite_samples=overwrite_bootstraps,
                                                                       pattern_speed_filename=pattern_speed_filename,
                                                                       )

                            if plot:

                                n_dim = 2 + n_pattern_speeds

                                labels = [r'$r_0$', 'c']
                                labels.extend([r'$m_%d$' % (i + 1) for i in range(n_pattern_speeds)])

                                fig, axes = plt.subplots(n_dim, figsize=(10, 7), sharex=True)
                                samples = sampler.get_chain()
                                for i in range(n_dim):
                                    ax = axes[i]
                                    ax.plot(samples[:, :, i], "k", alpha=0.3)
                                    ax.set_xlim(0, len(samples))
                                    ax.set_ylabel(labels[i])
                                    ax.yaxis.set_label_coords(-0.1, 0.5)

                                axes[-1].set_xlabel("Step number")

                                plt.savefig(pattern_speed_plot_filename + '_steps.png',
                                            bbox_inches='tight')
                                plt.savefig(pattern_speed_plot_filename + '_steps.pdf',
                                            bbox_inches='tight')
                                plt.close()

                                flat_samples = sampler.get_chain(discard=int(n_steps / 2), flat=True)

                                fig = corner.corner(flat_samples, labels=labels,
                                                    truth_color='k')

                                plt.savefig(pattern_speed_plot_filename + '_corner.png',
                                            bbox_inches='tight')
                                plt.savefig(pattern_speed_plot_filename + '_corner.pdf',
                                            bbox_inches='tight')
                                plt.close()

                                # Finally, plot on best-fit line

                                x_min, x_max = 1.2 * np.nanmin(x_tw), 1.2 * np.nanmax(x_tw)
                                x_fit = np.linspace(x_min, x_max, 500)

                                y_min, y_max = 1.2 * np.nanmin(v_tw), 1.7 * np.nanmax(v_tw)
                                y_median = np.zeros([n_pattern_speeds, len(x_fit)])
                                y_upper = y_median.copy()
                                y_lower = y_median.copy()

                                n_lines = 200

                                # Use the samples to calculate the values and uncertainties

                                for i in range(n_pattern_speeds):

                                    y_fit_single = np.zeros([n_lines, len(x_fit)])

                                    for line in range(n_lines):
                                        choice = np.random.randint(low=0, high=flat_samples.shape[0])

                                        y_fit_single[line, :] = x_fit * flat_samples[choice, 2 + i] + flat_samples[
                                            choice, 1]

                                    y_median[i, :] = np.nanmedian(y_fit_single, axis=0)
                                    y_upper[i, :] = np.percentile(y_fit_single, 84, axis=0)
                                    y_lower[i, :] = np.percentile(y_fit_single, 16, axis=0)

                                # Calculate the position of the slits

                                slit_positions = bar.y_slits / grid_step

                                # Calculate the axis line

                                rad = np.sqrt((flux.shape[0] - y_cen) ** 2 + (flux.shape[1] - x_cen) ** 2)
                                ang = [0, np.pi] + np.radians(pa)

                                # Plotting

                                # Set up colourbar limits for the image

                                vmin_flux = np.nanpercentile(flux, 0.5)
                                vmax_flux = np.nanpercentile(flux, 99.75)

                                vmax_vel = 150
                                vmin_vel = -150

                                plt.figure(figsize=(6, 6))

                                plt.subplot(221)

                                plt.imshow(flux,
                                           cmap=cmocean.cm.gray_r,
                                           origin='lower', interpolation='none',
                                           vmin=vmin_flux, vmax=vmax_flux)

                                colour = iter(cmocean.cm.thermal(np.linspace(0, 1, len(slit_positions))))

                                # Plot only 25% of the slits

                                n_slits = 4

                                # Pull out the LoN matrix.

                                lon_matrix = bar._mat_lon * bar._mat_NE

                                for i, slit_position in enumerate(slit_positions):
                                    c = next(colour)

                                    if i % n_slits == 0:
                                        # Use the rotation matrix to transform into the LoN plane.

                                        x0, x1 = -flux.shape[1], flux.shape[1]

                                        coord_rotate = np.asarray(np.array([x0, slit_position]) * lon_matrix)
                                        x0_rotate, y0_rotate = coord_rotate[0][0], coord_rotate[0][1]

                                        coord_rotate = np.asarray(np.array([x1, slit_position]) * lon_matrix)
                                        x1_rotate, y1_rotate = coord_rotate[0][0], coord_rotate[0][1]

                                        plt.plot([x0_rotate + x_cen, x1_rotate + x_cen],
                                                 [y0_rotate + y_cen, y1_rotate + y_cen],
                                                 c=c, linewidth=1,
                                                 )

                                plt.xlim([0, flux.shape[1]])
                                plt.ylim([0, flux.shape[0]])

                                plt.axis('off')

                                plt.subplot(222)

                                plt.imshow(vel,
                                           cmap=cmocean.cm.balance,
                                           origin='lower', interpolation='none',
                                           vmin=vmin_vel, vmax=vmax_vel)

                                # Also plot on the axis line

                                plt.plot(-rad * np.sin(ang) + x_cen, rad * np.cos(ang) + y_cen, 'k--', linewidth=2)

                                plt.xlim([0, flux.shape[1]])
                                plt.ylim([0, flux.shape[0]])

                                plt.axis('off')

                                ax = plt.subplot(212)

                                colour = iter(cmocean.cm.thermal(np.linspace(0, 1, len(x_tw))))

                                for i in range(len(x_tw)):
                                    c = next(colour)

                                    plt.errorbar(x_tw[i], v_tw[i],
                                                 xerr=x_tw_err[i], yerr=v_tw_err[i],
                                                 marker='o', c=c, linestyle='none')

                                for i in range(n_pattern_speeds):
                                    plt.plot(x_fit, y_median[i, :], c='k', lw=2)
                                    plt.fill_between(x_fit, y_upper[i, :], y_lower[i, :],
                                                     edgecolor='k', facecolor='k', alpha=0.25)

                                plt.xlim([x_min, x_max])
                                plt.ylim([y_min, y_max])

                                plot_text = ''

                                for i in range(n_pattern_speeds):
                                    plot_text += r'$\Omega_{p,%d} = %.2f^{+%.2f}_{-%.2f}\, \mathrm{km\,s}^{-1} ' \
                                                 r'\mathrm{''kpc}^{-1}$' \
                                                 % (i + 1, pattern_speeds[i, 0], pattern_speeds[i, 1],
                                                    pattern_speeds[i, 2])
                                    if i != n_pattern_speeds - 1:
                                        plot_text += '\n'

                                plt.text(0.05, 0.95,
                                         plot_text,
                                         transform=ax.transAxes,
                                         verticalalignment='top',
                                         bbox=dict(facecolor='white', alpha=0.7))

                                plt.xlabel(r'<$x$> $\left(^{\prime \prime}\right)$')
                                plt.ylabel(r'<$v$> $\left(\mathrm{km\,s}^{-1}\right)$')

                                plt.savefig(pattern_speed_plot_filename + '.png',
                                            bbox_inches='tight')
                                plt.savefig(pattern_speed_plot_filename + '.pdf',
                                            bbox_inches='tight')

print('Complete! Took %.2fm' % ((time.time() - start) / 60))
