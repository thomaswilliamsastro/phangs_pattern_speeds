# -*- coding: utf-8 -*-
"""
Apply Tremaine-Weinberg method to all of the PHANGS MUSE galaxies

@author: Tom Williams
"""

import os
import warnings

import cmocean
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pybar.pybar as pybar
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

import ps_functions
from bootstraps import bootstrap_tw
from muse.folders import phangs_folder, plot_folder, output_folder

warnings.simplefilter("ignore")

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

# Set up random seed for testing.

np.random.seed(42)

os.chdir(phangs_folder)

if not os.path.exists(plot_folder):
    os.mkdir(plot_folder)
if not os.path.exists(output_folder):
    os.mkdir(output_folder)

# Galaxies to fit

galaxies = ['NGC1087', 'NGC1512', 'NGC1672', 'NGC3351', 'NGC4254', 'NGC5068',
            'IC5332', 'NGC1365', 'NGC1566', 'NGC2835', 'NGC3627', 'NGC4535', 'NGC628']
galaxies = sorted(galaxies)

# Pixel size (arcsec)

grid_step = 0.2

# Read in the basic galaxy info

galaxy_table = fits.open('documents/phangs_sample_table_v1p1.fits')
galaxy_table = Table(galaxy_table[1].data)

# Read in bar information for later masking

bar_galaxy, bar_rs = np.loadtxt('environment/PHANGSmasks_v2.dat',
                                usecols=(0, 12),
                                unpack=True,
                                skiprows=4,
                                dtype=str)

overwrite_pafit = False
overwrite_bootstraps = True
plot = True

# Various options we can play around with for testing

hdu_types = ['mass']
mask_outside_bars = [True]
star_masks = [True]
slit_widths = [1]
slit_lengths = [0]

# For testing, here are the settings. Trying to do all of these things simultaneously will cause weird behaviour so
# they're grouped as they work. For emission:

# hdu_types = ['mass', 'flux', 'ha']
# mask_outside_bars = [True, False]
# star_masks = [True, False]

# For slit widths:

# start = 1
# stop = 10
# step = 0.5
#
# slit_widths = np.arange(start, stop + step, step)

# For slit lengths:

# start = 30
# stop = 110
# step = 5
#
# slit_lengths = np.arange(start, stop + step, step)

for galaxy in galaxies:

    print('Beginning ' + galaxy)

    for hdu_type in hdu_types:

        for mask_outside_bar in mask_outside_bars:

            for star_mask in star_masks:

                for slit_width in slit_widths:

                    for slit_length in slit_lengths:

                        # If the length of any of 'hdu_types', 'slit_widths', or 'slit_lengths is  greater than 1 we're
                        # in test mode so move into galaxy folder accordingly

                        if len(hdu_types) > 1 or len(slit_widths) > 1 or len(slit_lengths) > 1:

                            file_name = galaxy + '/'

                            if not os.path.exists(plot_folder+galaxy):
                                os.mkdir(plot_folder+galaxy)
                            if not os.path.exists(output_folder+galaxy):
                                os.mkdir(output_folder+galaxy)

                            # Let's set this so only one thing at a time can be changed. If we're changing HDU types,
                            # it's going in an emission folder

                            if len(hdu_types) > 1:

                                file_name += 'emission/'

                            # If the slit lengths are changing, it's going in a slit_length folder

                            elif len(slit_lengths) > 1:

                                file_name += 'slit_length/'

                            # If the slit widths are changing, it's going in a slit_width folder

                            elif len(slit_widths) > 1:

                                file_name += 'slit_width/'

                            if not os.path.exists(plot_folder+file_name):
                                os.mkdir(plot_folder+file_name)
                            if not os.path.exists(output_folder+file_name):
                                os.mkdir(output_folder+file_name)

                        else:

                            file_name = ''

                        # Set up file extensions for masks

                        if star_mask:
                            star_ext = '_smask'
                        else:
                            star_ext = '_nosmask'

                        if mask_outside_bar:
                            mask_ext = '_bmask'
                        else:
                            mask_ext = '_nobmask'

                        # Set up file names

                        file_name += galaxy + '_' + hdu_type + star_ext + mask_ext + '_'

                        if len(slit_lengths) > 1:

                            file_name += 'sl_%.0f_' % slit_length

                        elif len(slit_widths) > 1:

                            file_name += 'sw_%.1f_' % slit_width

                        pattern_speed_plot_filename = plot_folder + file_name + 'muse'
                        pattern_speed_filename = output_folder + file_name + 'pattern_speed_muse.txt'
                        bootstrap_filename = output_folder + file_name + 'bootstrap_muse.txt'

                        # Read in the various .fits images we'll need.

                        hdu = fits.open('muse/' + galaxy + '_MAPS.fits')
                        wcs = WCS(hdu[1])

                        # Flux-weighted

                        if hdu_type == 'flux':

                            flux = hdu[2].data
                            snr = hdu[3].data
                            flux_err = flux / snr

                            vel = hdu[6].data
                            vel_err = hdu[7].data

                        # Mass weighted

                        elif hdu_type == 'mass':

                            flux = hdu[115].data
                            flux_err = hdu[116].data

                            vel = hdu[6].data
                            vel_err = hdu[7].data

                        # H-alpha

                        elif hdu_type == 'ha':

                            flux = hdu[80].data
                            flux_err = hdu[81].data

                            vel = hdu[82].data
                            vel_err = hdu[83].data

                        else:

                            raise Warning('%s not a valid selection!' % hdu_type)

                        if star_mask:
                            # Do a simple star mask by velocities.

                            star_mask_idx = ps_functions.star_mask(vel)

                            flux[star_mask_idx] = np.nan
                            flux_err[star_mask_idx] = np.nan
                            vel[star_mask_idx] = np.nan
                            vel_err[star_mask_idx] = np.nan

                        ra, dec, inclination, dist, pa, pa_err = ps_functions.get_sample_table_info(galaxy,
                                                                                                    galaxy_table)

                        # Calculate a physical distance conversion factor from the galaxy distance.

                        kpc_per_arcsec = dist * 1e3 * np.sin(1 / 3600 * np.pi / 180)

                        # Convert the RA and Dec to a central pixel.

                        x_cen, y_cen = wcs.all_world2pix(ra, dec, 1)

                        if not os.path.exists(
                                'pattern_speeds_output/pafit/' + galaxy + '_pafit_output_muse.txt') or overwrite_pafit:

                            pafit_pa, pafit_pa_err, v_sys = ps_functions.pafit_wrapper(vel, vel_err, x_cen, y_cen)

                            np.savetxt('pattern_speeds_output/pafit/' + galaxy + '_pafit_output_muse.txt',
                                       np.c_[pafit_pa, pafit_pa_err, v_sys],
                                       header='kin_pa, kin_pa_err, v_sys')

                        else:

                            pafit_pa, pafit_pa_err, v_sys = np.loadtxt(
                                'pattern_speeds_output/pafit/' + galaxy + '_pafit_output_muse.txt',
                                unpack=True)

                        # print('pafit/Lang fit: %.2f/%.2f' % (pafit_pa, pa))

                        # Subtract any residual velocity

                        vel -= v_sys

                        # Set up a grid of x- and y-coordinates

                        grid_shape = flux.shape
                        grid_x = (np.arange(grid_shape[1]) - x_cen) * grid_step
                        grid_y = (np.arange(grid_shape[0]) - y_cen) * grid_step

                        x_coords, y_coords = np.meshgrid(grid_x, grid_y)

                        if mask_outside_bar:

                            # Mask anything Beyond the Bar in line with the kinematic PA.

                            try:

                                bar_r = np.float(bar_rs[bar_galaxy == galaxy])

                                if bar_r > 0:

                                    idx = ps_functions.bar_mask(x_coords, y_coords, pa, bar_r)

                                    flux[idx] = np.nan
                                    flux_err[idx] = np.nan
                                    vel[idx] = np.nan
                                    vel_err[idx] = np.nan

                                else:

                                    print('Bar not defined. Skipping masking')

                            except TypeError:

                                print('No bar information found. Not masking anything.')

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

                        if plot:

                            # Plot on a best fit line and associated errors for this single bootstrap

                            # Convert fitted numbers to something usefully physical -- km/s/kpc, and get rid of the
                            # inclination

                            km_s_kpc_conversion = (np.sin(inclination * np.pi / 180) * kpc_per_arcsec) ** -1

                            omega_bar = m.copy()
                            omega_bar_err = m_err.copy()

                            omega_bar *= km_s_kpc_conversion
                            omega_bar_err *= km_s_kpc_conversion

                            x_min, x_max = 1.2 * np.nanmin(x_tw), 1.2 * np.nanmax(x_tw)
                            x_fit = np.linspace(x_min, x_max, 500)

                            y_min, y_max = 1.2 * np.nanmin(v_tw), 1.7 * np.nanmax(v_tw)

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

                            plt.figure(figsize=(6, 6))

                            plt.subplot(221)

                            plt.imshow(np.log10(flux),
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

                            plt.tight_layout()

                            plt.savefig(pattern_speed_plot_filename + '.png',
                                        bbox_inches='tight')
                            plt.savefig(pattern_speed_plot_filename + '.pdf',
                                        bbox_inches='tight')
                            plt.close()

                        # Put everything into the bootstrap

                        bootstrap_tw(flux, flux_err, vel, vel_err, grid_step=grid_step, slit_width=slit_width,
                                     x_cen=x_cen, y_cen=y_cen, centering_err=1, pa=pa, pa_err=pa_err,
                                     inclination=inclination,
                                     dist=dist, bootstrap_filename=bootstrap_filename,
                                     overwrite_bootstraps=overwrite_bootstraps,
                                     n_bootstraps=1000, pattern_speed_filename=pattern_speed_filename,
                                     )

print('Complete!')