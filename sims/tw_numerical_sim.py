# -*- coding: utf-8 -*-
"""
Test on Eric's numerical simulations

@author: Tom Williams
"""

import os
import time
import warnings

import cmocean
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pybar.pybar as pybar
from astropy.io import fits

import ps_functions
from bootstraps import bootstrap_tw
from vars import phangs_folder, plot_folder, output_folder, slit_lengths, plot

warnings.simplefilter("ignore")

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

# Set up random seed for testing.

np.random.seed(42)

start = time.time()

os.chdir(phangs_folder)

if not os.path.exists(plot_folder):
    os.mkdir(plot_folder)
if not os.path.exists(output_folder):
    os.mkdir(output_folder)
if not os.path.exists(plot_folder + 'num_sim'):
    os.mkdir(plot_folder + 'num_sim')
if not os.path.exists(output_folder + 'num_sim'):
    os.mkdir(output_folder + 'num_sim')

galaxies = ['n628_70', 'n3351_02', 'n3351_20']
inclinations = ['09', '60']
emissions = ['gas', 'ns', 'os']

# start_cf = 0.5
# stop_cf = 3
# step_cf = 0.25
# covering_factors = np.arange(start_cf, stop_cf+step_cf, step_cf)
covering_factors = [-99]

n_pattern_speeds = 1
pa = 0
pa_err = 1
dist = 10
pix_size = 1

slit_widths = [pix_size * 3]

overwrite_bootstraps = False

for galaxy in galaxies:

    snapshots = {'n628_70': ['0278'],
                 'n3351_02': ['0065', '0226', '0268'],
                 'n3351_20': ['0097', '0050']}[galaxy]

    for snapshot in snapshots:

        for emission in emissions:

            for inclination in inclinations:

                for slit_width in slit_widths:

                    for slit_length in slit_lengths:

                        if len(covering_factors) > 1:
                            masked_fluxes = []

                        for covering_factor in covering_factors:

                            hdu_name = 'FigFits_2gal/Fits/'

                            hdu_name += galaxy + 'hr_' + emission + '_' + snapshot + '_M_i' + inclination + '.fits'

                            try:
                                flux_hdu = fits.open(hdu_name)[0]
                            except FileNotFoundError:
                                continue

                            # Set up file names depending on the various parameters

                            file_name = 'num_sim/'

                            if np.sum([len(slit_widths) > 1, len(slit_lengths) > 1]) > 1:
                                raise Warning('Too many things varying!')

                            if len(covering_factors) > 1:
                                file_name += 'covering_factors/'

                            if len(slit_widths) > 1:

                                file_name += 'slit_width/'

                            elif len(slit_lengths) > 1:

                                file_name += 'slit_length/'

                            if not os.path.exists(plot_folder + file_name):
                                os.mkdir(plot_folder + file_name)
                            if not os.path.exists(output_folder + file_name):
                                os.mkdir(output_folder + file_name)

                            if n_pattern_speeds > 1:
                                file_name += '%d_pattern_speeds/' % n_pattern_speeds

                            if not os.path.exists(plot_folder + file_name):
                                os.mkdir(plot_folder + file_name)
                            if not os.path.exists(output_folder + file_name):
                                os.mkdir(output_folder + file_name)

                            file_name += galaxy + '_' + snapshot + '_' + emission + '_i' + inclination

                            if len(slit_lengths) > 1:

                                file_name += '_sl_%.0f' % slit_length

                            elif len(slit_widths) > 1:

                                file_name += '_sw_%.1f' % slit_width

                            if len(covering_factors) > 1:
                                file_name += '_cf_%.2f' % covering_factor

                            pattern_speed_plot_filename = plot_folder + file_name
                            pattern_speed_filename = output_folder + file_name + '_pattern_speed.txt'
                            bootstrap_filename = output_folder + file_name + '_bootstrap.txt'
                            max_r_filename = output_folder + file_name + '_max_r.npy'

                            inclination = float(inclination)

                            # Read in the various .fits images we'll need. We don't have an error map here, use 10%

                            flux = flux_hdu.data
                            vel = fits.open(hdu_name.replace('_M_', '_V_'))[0].data

                            # If we're looking at covering factors, mask here

                            if len(covering_factors) > 1:

                                n_pix_orig = flux.shape[0]*flux.shape[1]
                                # n_pix_orig = len(np.where(np.isnan(flux) == False)[0])

                                mask = np.where(flux <= covering_factor)
                                flux[mask] = np.nan
                                vel[mask] = np.nan

                                n_pix_masked = len(np.where(np.isnan(flux) == False)[0])
                                masked_fluxes.append(n_pix_masked/n_pix_orig)

                            flux_err = flux * 0.1
                            vel_err = vel * 0.1

                            # Remove any 0s in the flux image, replace with NaNs

                            idx = np.where(flux == 0)

                            flux[idx] = np.nan
                            flux_err[idx] = np.nan
                            vel[idx] = np.nan
                            vel_err[idx] = np.nan

                            idx = np.where(np.isnan(flux))

                            if np.all(np.isnan(flux)):
                                print('No flux detected. Skipping')
                                continue

                            # Calculate a physical distance conversion factor from the galaxy distance.

                            kpc_per_arcsec = dist * 1e3 * np.sin(1 / 3600 * np.pi / 180)

                            # Find the brightest pixel, use that as the centre

                            y_cen, x_cen = np.where(flux == np.nanmax(flux))
                            y_cen = y_cen[0]
                            x_cen = x_cen[0]

                            # Set up a grid of x- and y-coordinates

                            grid_shape = flux.shape
                            grid_step = np.abs(pix_size)
                            grid_x = (np.arange(grid_shape[1]) - x_cen) * grid_step
                            grid_y = (np.arange(grid_shape[0]) - y_cen) * grid_step

                            x_coords, y_coords = np.meshgrid(grid_x, grid_y)

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

                            y_lon_round = bar.Y_lon.round(0)
                            x_lon = bar.X_lon

                            flux = ps_functions.symmetrize_tw_integral(flux, x_lon, y_lon_round)

                            flux_err[np.isnan(flux)] = np.nan
                            vel[np.isnan(flux)] = np.nan
                            vel_err[np.isnan(flux)] = np.nan

                            y_lon_reshape = bar.Y_lon.reshape(flux.shape)
                            y_lon_reshape[np.isnan(flux)] = np.nan

                            r_max = np.nanmax(np.abs(y_lon_reshape))

                            np.save(max_r_filename, r_max)

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

                            # Because some position angles are from HYPERLEDA, and only go from 0 to 180, sometimes
                            # they have the wrong sense. Check that here, and if this is the case swap around.

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

                                m, m_err, c, c_err = ps_functions.odr_fit(x_tw, x_tw_err, v_tw, v_tw_err)

                            # Convert fitted numbers to something usefully physical -- km/s/kpc, and get rid of the
                            # inclination

                            km_s_kpc_conversion = (np.sin(inclination * np.pi / 180) * kpc_per_arcsec) ** -1

                            omega_bar = m.copy()
                            omega_bar_err = m_err.copy()

                            omega_bar *= km_s_kpc_conversion
                            omega_bar_err *= km_s_kpc_conversion

                            if plot and n_pattern_speeds == 1:

                                # Plot on a best fit line and associated errors for this single bootstrap

                                x_max = 1.2 * np.nanpercentile(x_tw, 99)
                                x_min = 1.2 * np.nanpercentile(x_tw, 1)

                                # x_min, x_max = 1.2 * np.nanmin(x_tw), 1.2 * np.nanmax(x_tw)

                                x_fit = np.linspace(x_min, x_max, 500)

                                y_max = 1.2 * np.nanpercentile(v_tw, 99)
                                y_min = 1.2 * np.nanpercentile(v_tw, 1)

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

                                if len(covering_factors) > 1:
                                    plt.suptitle(str(covering_factor))

                                plt.tight_layout()

                                plt.savefig(pattern_speed_plot_filename + '.png',
                                            bbox_inches='tight')
                                plt.savefig(pattern_speed_plot_filename + '.pdf',
                                            bbox_inches='tight')
                                plt.close()

                            if n_pattern_speeds <= 0:

                                raise Warning('This is not a number of pattern speeds I can deal with')

                            elif n_pattern_speeds == 1:

                                # Put everything into the bootstrap

                                bootstrap_tw(flux, flux_err, vel, vel_err, grid_step=grid_step, slit_width=slit_width,
                                             x_cen=x_cen, y_cen=y_cen, centering_err=1, pa=pa, pa_err=pa_err,
                                             inclination=inclination,
                                             dist=dist, bootstrap_filename=bootstrap_filename,
                                             overwrite_bootstraps=overwrite_bootstraps,
                                             n_bootstraps=1000, pattern_speed_filename=pattern_speed_filename,
                                             )

                        if len(covering_factors) > 1:
                            np.savetxt(output_folder + 'spiral_models/covering_factors/' +
                                       galaxy + '_fractions.txt',
                                       masked_fluxes)

print('Complete! Took %.2fm' % ((time.time() - start) / 60))
