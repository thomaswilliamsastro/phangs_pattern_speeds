# -*- coding: utf-8 -*-
"""
Apply Tremaine-Weinberg method to all of the PHANGS MUSE galaxies

@author: Tom Williams
"""

import os
import sys
import time
import warnings

import cmocean
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pybar.pybar as pybar
from astropy.io import fits
from astropy.wcs import WCS
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from reproject import reproject_interp

sys.path.append(os.getcwd())

import ps_functions
from bootstraps import bootstrap_tw
from vars import phangs_folder, muse_version, muse_galaxies, muse_plot, muse_output, galaxy_table, overwrite_bootstraps, plot, hdu_types, mask_outside_bars, star_masks, slit_widths, slit_lengths, n_pattern_speeds

warnings.simplefilter("ignore")

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

# Set up random seed for testing.

np.random.seed(42)

start = time.time()

os.chdir(phangs_folder)

if not os.path.exists(muse_plot):
    os.mkdir(muse_plot)
if not os.path.exists(muse_output):
    os.mkdir(muse_output)
if not os.path.exists(muse_plot + muse_version):
    os.mkdir(muse_plot + muse_version)
if not os.path.exists(muse_output + muse_version):
    os.mkdir(muse_output + muse_version)

# Pixel size (arcsec)

grid_step = 0.2

null_hypothesis = False
fov_check = False

# muse_galaxies = ['ngc3351']

for galaxy in muse_galaxies:

    galaxy = galaxy.upper()

    if muse_version == 'DR1.0' and galaxy == 'NGC0628':
        galaxy = 'NGC628'

    print('Beginning ' + galaxy)

    for hdu_type in hdu_types:

        # Sigma-masking to test covering factor effects

        if hdu_type == 'toy_sim_cov':
            cov_factor_sigmas = np.arange(1, 11)
        else:
            cov_factor_sigmas = [1]

        for cov_factor_sigma in cov_factor_sigmas:

            for mask_outside_bar in mask_outside_bars:

                for star_mask in star_masks:

                    for slit_width in slit_widths:

                        for slit_length in slit_lengths:

                            # Set up file names depending on various parameters

                            file_name = muse_version + '/'

                            if np.sum([len(slit_widths) > 1, len(slit_lengths) > 1, len(mask_outside_bars) > 1]) > 1:
                                raise Warning('Too many things varying!')

                            if len(slit_widths) > 1:

                                file_name += 'slit_width/'

                            elif len(slit_lengths) > 1:

                                file_name += 'slit_length/'

                            elif mask_outside_bar and star_mask:

                                file_name += hdu_type + '_smask_bmask/'

                            elif mask_outside_bar and not star_mask:

                                file_name += hdu_type + '_nosmask_bmask/'

                            elif not mask_outside_bar and star_mask:

                                file_name += hdu_type + '_smask_nobmask/'

                            elif not mask_outside_bar and not star_mask:

                                file_name += hdu_type + '_nosmask_nobmask/'

                            else:

                                raise Warning('Not sure what is going on here!')

                            if n_pattern_speeds > 1:
                                file_name += '%d_pattern_speeds/' % n_pattern_speeds

                            if null_hypothesis:
                                file_name += 'null/'

                            if fov_check:
                                file_name += 'fov_check/'

                            if not os.path.exists(muse_plot + file_name):
                                os.mkdir(muse_plot + file_name)
                            if not os.path.exists(muse_output + file_name):
                                os.mkdir(muse_output + file_name)

                            file_name += galaxy

                            # Set up file extensions for mask

                            file_name += '_' + hdu_type

                            if star_mask:
                                file_name += '_smask'
                            else:
                                file_name += '_nosmask'

                            if mask_outside_bar:
                                file_name += '_bmask_'
                            else:
                                file_name += '_nobmask_'

                            if len(slit_lengths) > 1:

                                file_name += 'sl_%.0f_' % slit_length

                            elif len(slit_widths) > 1:

                                file_name += 'sw_%.1f_' % slit_width

                            if len(cov_factor_sigmas) > 1:
                                file_name += 'cov_sig_%d_' % cov_factor_sigma

                            if null_hypothesis:
                                file_name += 'null_'

                            if fov_check:
                                file_name += 'fov_check_'

                            pattern_speed_plot_filename = muse_plot + file_name + 'muse'
                            pattern_speed_filename = muse_output + file_name + 'pattern_speed_muse.txt'
                            bootstrap_filename = muse_output + file_name + 'bootstrap_muse.txt'
                            max_r_filename = muse_output + file_name + 'max_r.npy'

                            # Read in the various .fits images we'll need.

                            try:
                                hdu = fits.open('muse/' + muse_version + '/' + galaxy + '_MAPS.fits')
                            except FileNotFoundError:
                                print('%s maps not found: skipping' % galaxy)
                                continue
                            wcs = WCS(hdu[1])

                            # Flux-weighted

                            if hdu_type == 'flux':

                                flux = hdu['FLUX'].data
                                snr = hdu['SNR'].data
                                flux_err = flux / snr

                                vel = hdu['V_STARS'].data
                                vel_err = hdu['FORM_ERR_V_STARS'].data

                            # Mass weighted

                            elif 'mass' in hdu_type:

                                flux = hdu['STELLAR_MASS_DENSITY'].data
                                flux_err = hdu['STELLAR_MASS_DENSITY_ERR'].data

                                if 'spitzer' in hdu_type:
                                    # Read in the Spitzer heaader and reproject to the MUSE scale

                                    hdr = hdu['STELLAR_MASS_DENSITY'].header
                                    nans = np.where(np.isnan(flux.data))

                                    flux = fits.open('spitzer_irac/' + galaxy + '.phot.2.rs.fits')[0]

                                    flux, _ = reproject_interp(flux, hdr)
                                    flux[nans] = np.nan

                                    # Throw in a 10% flux uncertainty for now
                                    flux_err = 0.1 * flux

                                vel = hdu['V_STARS'].data
                                vel_err = hdu['FORM_ERR_V_STARS'].data

                            # H-alpha

                            elif hdu_type == 'ha':

                                flux = hdu['HA6562_FLUX'].data
                                flux_err = hdu['HA6562_FLUX_ERR'].data

                                vel = hdu['HA6562_VEL'].data
                                vel_err = hdu['HA6562_VEL_ERR'].data

                            elif 'toy_sim' in hdu_type:

                                # Setup arrays for the simulation

                                flux = np.zeros(hdu['HA6562_FLUX'].data.shape)
                                flux_err = flux.copy()

                                vel = np.zeros(hdu['HA6562_VEL'].data.shape)
                                vel_err = vel.copy()

                                # We also want a mask to mimic the FOV of the observations
                                mask = hdu['HA6562_FLUX'].data
                                mask = np.isnan(mask)

                                # If we're testing the covering factor, put that in here

                                if 'cov' in hdu_type:
                                    obs_flux = hdu['HA6562_FLUX'].data
                                    obs_flux_err = hdu['HA6562_FLUX_ERR'].data

                                    mask[obs_flux / obs_flux_err <= cov_factor_sigma] = np.nan

                            else:

                                raise Warning('%s not a valid selection!' % hdu_type)

                            flux_pristine = flux.copy()
                            vel_pristine = vel.copy()

                            if star_mask:
                                # Do a simple star mask by velocities.

                                star_mask_idx = ps_functions.star_mask(vel)

                                flux[star_mask_idx] = np.nan
                                flux_err[star_mask_idx] = np.nan
                                vel[star_mask_idx] = np.nan
                                vel_err[star_mask_idx] = np.nan

                            # Pull out parameters we need from the sample table.

                            sample_table_parameters = ['orient_ra', 'orient_dec',
                                                       'orient_ra_unc', 'orient_dec_unc',
                                                       'orient_incl', 'dist',
                                                       'orient_posang', 'orient_posang_unc',
                                                       'morph_bar_r'
                                                       ]

                            galaxy_edit = galaxy
                            if galaxy_edit == 'NGC628':
                                galaxy_edit = 'NGC0628'

                            parameters_dict = ps_functions.get_sample_table_info(galaxy_edit.lower(), galaxy_table,
                                                                                 sample_table_parameters)

                            ra = parameters_dict['orient_ra']
                            dec = parameters_dict['orient_dec']
                            centering_err = (parameters_dict['orient_ra_unc'] + parameters_dict['orient_dec_unc']) / 2
                            inclination = parameters_dict['orient_incl']
                            pa = parameters_dict['orient_posang']
                            pa_err = parameters_dict['orient_posang_unc']
                            dist = parameters_dict['dist']
                            bar_r = parameters_dict['morph_bar_r']

                            # If the galaxy doesn't have a measured PA, skip

                            if np.isnan(pa):
                                print('PA undefined: skipping')
                                continue

                            # If we want to check the effects of field of view, force in the test HDU here

                            if fov_check:
                                coverage = np.zeros_like(flux)
                                coverage[~np.isnan(flux)] = 1

                                fov_hdu = fits.open('muse/' + muse_version + '/NGC0628_MAPS.fits')
                                fov_flux = fov_hdu['FLUX'].data
                                fov_snr = fov_hdu['SNR'].data
                                fov_flux_err = fov_flux / fov_snr

                                fov_vel = fov_hdu['V_STARS'].data
                                fov_vel_err = fov_hdu['FORM_ERR_V_STARS'].data

                                # Convert the original RA and Dec to a central pixel.

                                x_cen_orig, y_cen_orig = wcs.all_world2pix(ra, dec, 1)
                                x_cen_orig, y_cen_orig = np.round(x_cen_orig), np.round(y_cen_orig)

                                parameters_dict = ps_functions.get_sample_table_info('ngc0628', galaxy_table,
                                                                                     sample_table_parameters)

                                ra = parameters_dict['orient_ra']
                                dec = parameters_dict['orient_dec']
                                centering_err = (parameters_dict['orient_ra_unc'] + parameters_dict[
                                    'orient_dec_unc']) / 2
                                inclination = parameters_dict['orient_incl']
                                pa = parameters_dict['orient_posang']
                                pa_err = parameters_dict['orient_posang_unc']
                                dist = parameters_dict['dist']
                                bar_r = parameters_dict['morph_bar_r']

                                wcs = WCS(fov_hdu[1])

                                x_cen, y_cen = wcs.all_world2pix(ra, dec, 1)
                                x_cen_round, y_cen_round = np.round(x_cen), np.round(y_cen)

                                x_pad = [int(x_cen_round - x_cen_orig),
                                         int(fov_flux.shape[1] - flux.shape[1] - x_cen_round + x_cen_orig)]
                                y_pad = [int(y_cen_round - y_cen_orig),
                                         int(fov_flux.shape[0] - flux.shape[0] - y_cen_round + y_cen_orig)]

                                # If any of these are negative, then actually trim off
                                slice_x_start = 0
                                slice_x_end = None
                                slice_y_start = 0
                                slice_y_end = None
                                if x_pad[0] < 0:
                                    slice_x_start = - x_pad[0]
                                    x_pad[0] = 0
                                if x_pad[1] < 0:
                                    slice_x_end = x_pad[1]
                                    x_pad[1] = 0
                                if y_pad[0] < 0:
                                    slice_y_start = - y_pad[0]
                                    y_pad[0] = 0
                                if y_pad[1] < 0:
                                    slice_y_end = y_pad[1]
                                    y_pad[1] = 0

                                coverage_pad = np.pad(coverage, (tuple(y_pad), tuple(x_pad)),
                                                      'constant', constant_values=0)[
                                               slice_y_start:slice_y_end, slice_x_start:slice_x_end]

                                fov_flux[coverage_pad == 0] = np.nan
                                fov_flux_err[coverage_pad == 0] = np.nan
                                fov_vel[coverage_pad == 0] = np.nan
                                fov_vel_err[coverage_pad == 0] = np.nan

                                flux = fov_flux
                                flux_err = fov_flux_err
                                vel = fov_vel
                                vel_err = fov_vel_err

                            # Calculate a physical distance conversion factor from the galaxy distance.

                            kpc_per_arcsec = dist * 1e3 * np.sin(1 / 3600 * np.pi / 180)

                            # Convert the RA and Dec to a central pixel.

                            x_cen, y_cen = wcs.all_world2pix(ra, dec, 1)

                            # If we're doing the simulation, now we'll figure out the flux and velocity field

                            if 'toy_sim' in hdu_type:
                                # Set up a grid of x- and y-coordinates

                                grid_shape = flux.shape
                                grid_x = (np.arange(grid_shape[1]) - x_cen) * grid_step
                                grid_y = (np.arange(grid_shape[0]) - y_cen) * grid_step

                                x_coords, y_coords = np.meshgrid(grid_x, grid_y)

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
                                flux_err = 0.1 * flux

                                # Mask to match the observational FOV

                                flux[mask] = np.nan
                                flux_err[mask] = np.nan
                                vel[mask] = np.nan
                                vel_err[mask] = np.nan

                            # if not os.path.exists('pattern_speeds_output/pafit/' + galaxy + '_pafit_output_muse.txt') \
                            #         or overwrite_pafit:
                            #
                            #     pafit_pa, pafit_pa_err, v_sys = ps_functions.pafit_wrapper(vel, vel_err, x_cen, y_cen)
                            #
                            #     np.savetxt('pattern_speeds_output/pafit/' + galaxy + '_pafit_output_muse.txt',
                            #                np.c_[pafit_pa, pafit_pa_err, v_sys],
                            #                header='kin_pa, kin_pa_err, v_sys')
                            #
                            # else:
                            #
                            #     pafit_pa, pafit_pa_err, v_sys = np.loadtxt(
                            #         'pattern_speeds_output/pafit/' + galaxy + '_pafit_output_muse.txt',
                            #         unpack=True)

                            v_sys = 0

                            # print('pafit/Lang fit: %.2f/%.2f' % (pafit_pa, pa))

                            # Subtract any residual velocity if we're not doing a simulation

                            if not 'toy_sim' in hdu_type:
                                vel -= v_sys

                            # Set up a grid of x- and y-coordinates

                            grid_shape = flux.shape
                            grid_x = (np.arange(grid_shape[1]) - x_cen) * grid_step
                            grid_y = (np.arange(grid_shape[0]) - y_cen) * grid_step

                            x_coords, y_coords = np.meshgrid(grid_x, grid_y)

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

                            # Finally, make sure the integrals will be symmetrical about the PA

                            bar = pybar.mybar(Flux=flux, Flux_err=flux_err,
                                              Velocity=vel, Velocity_err=vel_err,
                                              Xin=x_coords, Yin=y_coords, PAnodes=pa,
                                              inclin=inclination)

                            y_lon_round = bar.Y_lon.round(1)
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

                            if null_hypothesis:

                                r = np.sqrt(x_lon_reshape ** 2 + (y_lon_reshape / np.cos(np.radians(inclination))) ** 2)

                                r_edges = np.linspace(0, np.nanmax(r), 100)
                                spacing = (r_edges[1] - r_edges[0]) / 2

                                vel_null_hdu_name = os.path.join('halpha_rot_curves',
                                                                 galaxy_edit + '_model_RC_MUSE_Ha.fits')

                                if not os.path.exists(vel_null_hdu_name):
                                    print('%s vel not found' % galaxy)
                                    continue

                                else:
                                    vel_null = fits.open(vel_null_hdu_name)[0].data

                                # Make sure  the shapes match
                                if not np.all(vel_null.shape == flux.shape):
                                    continue

                                coverage = np.ones_like(flux)
                                coverage[(flux == 0) | (np.isnan(flux)) | (vel_null == 0)] = 0

                                flux_null = np.zeros_like(flux)
                                flux_null[flux_null == 0] = np.nan

                                for r_edge in r_edges:
                                    idx = np.where((r >= r_edge - spacing) & (r <= r_edge + spacing))

                                    flux_null[idx] = np.nanmedian(flux[idx])

                                flux_null[coverage == 0] = np.nan
                                vel_null[coverage == 0] = np.nan

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

                            m, m_err, c, c_err = ps_functions.odr_fit(x_tw, x_tw_err, v_tw, v_tw_err)

                            # Convert fitted numbers to something usefully physical -- km/s/kpc, and get rid of the
                            # inclination

                            km_s_kpc_conversion = (np.sin(inclination * np.pi / 180) * kpc_per_arcsec) ** -1

                            omega_bar = m.copy()
                            omega_bar_err = m_err.copy()

                            omega_bar *= km_s_kpc_conversion
                            omega_bar_err *= km_s_kpc_conversion

                            if plot:

                                # Figure out the 1kpc scale

                                kpc_scale = 1 / (grid_step/3600 * np.pi/180 * dist * 1e3)

                                # Plot on a best fit line and associated errors for this single bootstrap

                                x_max = 1.2 * np.nanpercentile(x_tw, 99)
                                x_min = 1.2 * np.nanpercentile(x_tw, 1)
                                # x_min, x_max = 1.2 * np.nanmin(x_tw), 1.2 * np.nanmax(x_tw)
                                x_fit = np.linspace(x_min, x_max, 500)

                                y_max = 1.2 * np.nanpercentile(v_tw, 98)
                                y_min = 1.2 * np.nanpercentile(v_tw, 2)
                                # y_min, y_max = 1.2 * np.nanmin(v_tw), 1.7 * np.nanmax(v_tw)

                                n_lines = 200

                                y_samples = np.zeros([n_lines, len(x_fit)])

                                for i in range(n_lines):
                                    y_samples[i, :] = np.random.normal(loc=m, scale=m_err) * x_fit + \
                                                      np.random.normal(loc=c, scale=c_err)

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

                                vmin_flux = np.nanpercentile(np.log10(flux[flux != 0]), 0.5)
                                vmax_flux = np.nanpercentile(np.log10(flux[flux != 0]), 99.75)

                                vmax_vel = 150
                                vmin_vel = -150

                                plt.figure(figsize=(8, 9))

                                ax = plt.subplot(221, projection=wcs)

                                plt.imshow(np.log10(flux_pristine),
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

                                plt.xlabel('RA (J2000)')
                                plt.ylabel('Dec (J2000)')

                                fancy_label = {'mass': r'$\mathbf{M_\ast}}$',
                                               'ha': r'H$\mathbf{\alpha}$'}[hdu_type]

                                fancy_text = '%s, %s' % (galaxy.upper(), fancy_label)

                                plt.text(0.05, 0.95, fancy_text, ha='left', va='top',
                                         fontweight='bold', transform=ax.transAxes,
                                         bbox=dict(facecolor='white', alpha=0.7))

                                scalebar = AnchoredSizeBar(ax.transData,
                                                           kpc_scale, '1kpc', 3,
                                                           pad=0.3,
                                                           borderpad=0.3,
                                                           color='black',
                                                           frameon=True)

                                ax.add_artist(scalebar)

                                ax.coords[0].set_ticks(exclude_overlapping=True)
                                ax.coords[1].set_ticks(exclude_overlapping=True)

                                # plt.axis('off')

                                ax = plt.subplot(222, projection=wcs)

                                plt.imshow(vel_pristine,
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

                                scalebar = AnchoredSizeBar(ax.transData,
                                                           kpc_scale, '1kpc', 3,
                                                           pad=0.3,
                                                           borderpad=0.3,
                                                           color='black',
                                                           frameon=True)

                                ax.add_artist(scalebar)

                                plt.text(0.05, 0.95, 'Velocity', ha='left', va='top',
                                         fontweight='bold', transform=ax.transAxes,
                                         bbox=dict(facecolor='white', alpha=0.7))

                                ax.coords[1].set_ticklabel_position('r')
                                ax.coords[1].set_axislabel_position('r')

                                ax.coords[0].set_ticks(exclude_overlapping=True)
                                ax.coords[1].set_ticks(exclude_overlapping=True)

                                ax = plt.subplot(212)

                                colour = iter(cmocean.cm.thermal(np.linspace(0, 1, len(x_tw))))

                                for i in range(len(x_tw)):
                                    c = next(colour)

                                    plt.errorbar(x_tw[i], v_tw[i],
                                                 xerr=x_tw_err[i], yerr=v_tw_err[i],
                                                 marker='o', c=c, linestyle='none')

                                plt.plot(x_fit, y_median, c='k', lw=2)
                                plt.fill_between(x_fit, y_upper, y_lower,
                                                 edgecolor='k', facecolor='k', alpha=0.25)

                                plt.xlim([x_min, x_max])
                                plt.ylim([y_min, y_max])

                                plt.grid(which='major', color='gray', linestyle='-', alpha=0.5)

                                plt.text(0.05, 0.85,
                                         r'$\Omega_\mathrm{P} = %.2f \pm %.2f\, \mathrm{km\,s}^{-1} \mathrm{kpc}^{-1}$ '
                                         % (omega_bar, omega_bar_err),
                                         transform=ax.transAxes,
                                         bbox=dict(facecolor='white', alpha=0.7))

                                plt.xlabel(r'<$x$> $\left(^{\prime \prime}\right)$')
                                plt.ylabel(r'<$v$> $\left(\mathrm{km\,s}^{-1}\right)$')

                                # plt.tight_layout()

                                # plt.show()

                                plt.savefig(pattern_speed_plot_filename + '.png', bbox_inches='tight')
                                plt.savefig(pattern_speed_plot_filename + '.pdf', bbox_inches='tight')

                                plt.close()

                            if 'toy_sim' in hdu_type:

                                # We're testing the FOV here, not the position angle. Save out the results of the one
                                # fit

                                np.savetxt(pattern_speed_filename, np.c_[omega_bar, omega_bar_err, omega_bar_err],
                                           header='omega_bar, omega_barr_err_up, omega_bar_err_down (all km/s/kpc)')

                            else:

                                # Put everything into the bootstrap

                                bootstrap_tw(flux, flux_err, vel, vel_err, grid_step=grid_step, slit_width=slit_width,
                                             x_cen=x_cen, y_cen=y_cen, centering_err=centering_err,
                                             pa=pa, pa_err=pa_err, inclination=inclination,
                                             dist=dist, bootstrap_filename=bootstrap_filename,
                                             overwrite_bootstraps=overwrite_bootstraps,
                                             n_bootstraps=1000, pattern_speed_filename=pattern_speed_filename,
                                             )

print('Complete! Took %.2fm' % ((time.time() - start) / 60))
