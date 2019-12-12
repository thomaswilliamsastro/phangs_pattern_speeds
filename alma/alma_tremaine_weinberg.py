# -*- coding: utf-8 -*-
"""
Apply Tremaine-Weinberg method to all of the PHANGS ALMA galaxies

@author: Tom Williams
"""

import copy
import os
import warnings

import cmocean
import corner
import emcee
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pybar.pybar as pybar
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from pafit.fit_kinematic_pa import fit_kinematic_pa as pafit
from scipy.stats import spearmanr


def ln_like(theta, x, y, y_err):
    m, c = theta
    model = m * x + c
    return -0.5 * np.nansum((y - model) ** 2 / y_err ** 2)


def ln_prob(theta, x, y, y_err):
    lp = ln_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + ln_like(theta, x, y, y_err)


def ln_prior(theta):

    # We just have totally unbounded priors here.

    return 0

# TODO THIS IS PROBABLY REDUNDANT.

warnings.simplefilter("ignore")

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

# Set up random seed for testing.

np.random.seed(42)

# Overwrite switches for slow bits of code.

overwrite_pafit = False
overwrite_emcee_samples = False

os.chdir('/Users/williams/Documents/phangs')

# Read in the basic galaxy info

galaxy_table = fits.open('documents/phangs_sample_table_v1p1.fits')
galaxy_table = Table(galaxy_table[1].data)

galaxies = galaxy_table['NAME'][galaxy_table['ALMA'] == 1]

galaxies = sorted(galaxies)

for galaxy in galaxies:

    # Remove any white space from the name

    galaxy = galaxy.strip()

    print('Beginning ' + galaxy)

    try:
        flux_hdu = fits.open('alma/' + galaxy + '_mom0.fits')[0]
    except FileNotFoundError:
        print(galaxy+' not found!')
        continue

    # Remove any extraneous third dimensions from the ALMA data

    del flux_hdu.header['PC3_1'], flux_hdu.header['PC3_2'], flux_hdu.header['PC1_3'], flux_hdu.header['PC2_3'], \
        flux_hdu.header['PC3_3'], flux_hdu.header['CTYPE3'], flux_hdu.header['CRVAL3'], flux_hdu.header['CDELT3'], \
        flux_hdu.header['CRPIX3'], flux_hdu.header['CUNIT3']

    wcs = WCS(flux_hdu)

    # Read in the various .fits images we'll need.

    flux = flux_hdu.data
    flux_err = fits.open('alma/' + galaxy + '_emom0.fits')[0].data

    vel = fits.open('alma/' + galaxy + '_mom1.fits')[0].data
    vel_err = fits.open('alma/' + galaxy + '_emom1.fits')[0].data

    # For no detections, we have a 0 so set these all to NaNs

    flux[flux == 0] = np.nan
    flux_err[flux_err == 0] = np.nan
    vel[vel == 0] = np.nan
    vel_err[vel_err == 0] = np.nan

    if len(np.where(np.isnan(flux) == False)[0]) == 0:
        print('No emission found')
        continue

    # From the sample table overview, pull out the RA/Dec, inclination and distance.

    galaxy_found = False

    galaxy_homogenised = copy.copy(galaxy)

    # If the leading 0 for the galaxy is missing, force it in.

    if galaxy_homogenised[:3] == 'NGC' and len(galaxy_homogenised) == 6:
        galaxy_homogenised = galaxy_homogenised[:3]+'0'+galaxy_homogenised[3:]

    for i in range(len(galaxy_table)):
        if galaxy_table['NAME'][i].strip() == galaxy_homogenised:
            ra, dec, inclination, dist, kinematic_pa = galaxy_table['RA_DEG'][i], galaxy_table['DEC_DEG'][i], \
                                                        galaxy_table['INCL'][i], galaxy_table['DIST'][i], \
                                                        galaxy_table['POSANG'][i]
            galaxy_found = True
            break

    # This won't work for an inclination of 90 so if it is, skip

    if inclination == 90:
        print('Edge-on: skipping')
        continue

    # If the galaxy isn't found in the table, break!

    if not galaxy_found:
        raise Warning(galaxy + ' not found!')

    # If the galaxy doesn't have a measured PA, skip

    if np.isnan(kinematic_pa):
        print('PA undefined: skipping')
        continue

    # Calculate a physical distance conversion factor from the galaxy distance.

    kpc_per_arcsec = dist * 1e3 * np.sin(1 / 3600 * np.pi / 180)

    # Convert the RA and Dec to a central pixel.

    x_cen, y_cen = wcs.all_world2pix(ra, dec, 1)

    if not os.path.exists('pattern_speeds_output/' + galaxy + '_pafit_output_alma.txt') or overwrite_pafit:

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

        # For speed, take 10% of the data.

        x_pafit = x_pafit[::10]

        y_pafit = y_pafit[idx].flatten()[::10]
        vel_pafit = vel[idx].flatten()[::10]
        vel_err_pafit = vel_err[idx].flatten()[::10]

        # Put this into pafit to get any residual systemic velocity and also the kinematic PA

        pafit_pa, pafit_pa_err, v_sys = pafit(x_pafit, y_pafit, vel_pafit, dvel=vel_err_pafit,
                                              quiet=True, plot=False)

        np.savetxt('pattern_speeds_output/' + galaxy + '_pafit_output_alma.txt',
                   np.c_[pafit_pa, pafit_pa_err, v_sys],
                   header='kin_pa, kin_pa_err, v_sys')

    else:

        pafit_pa, pafit_pa_err, v_sys = np.loadtxt('pattern_speeds_output/' + galaxy + '_pafit_output_alma.txt',
                                                   unpack=True)

    print('pafit/Lang fit: %.2f/%.2f' % (pafit_pa, kinematic_pa))

    # Subtract any residual velocity

    vel -= v_sys

    grid_shape = vel.shape
    grid_step = np.abs(flux_hdu.header['CDELT1']*3600)
    grid_x = (np.arange(grid_shape[1]) - x_cen) * grid_step
    grid_y = (np.arange(grid_shape[0]) - y_cen) * grid_step

    x_coords, y_coords = np.meshgrid(grid_x, grid_y)

    # Set up pybar. Use default slit width (for now).

    slit_width = 1

    bar = pybar.mybar(Flux=flux, Flux_err=flux_err,
                      Velocity=vel, Velocity_err=vel_err,
                      Xin=x_coords, Yin=y_coords, PAnodes=kinematic_pa,
                      inclin=inclination)

    bar.tremaine_weinberg(slit_width=slit_width)

    x_tw = bar.dfx_tw
    v_tw = bar.dfV_tw

    # As a quick test to see that the position angle is right, have a look at the correlation between these.

    correlation = spearmanr(x_tw[np.isnan(x_tw) == False], v_tw[np.isnan(x_tw) == False])[0]

    # If the correlation is negative, flip the PA round 180 degrees and repeat.

    if correlation < 0:

        kinematic_pa += 180

        bar = pybar.mybar(Flux=flux, Flux_err=flux_err,
                          Velocity=vel, Velocity_err=vel_err,
                          Xin=x_coords, Yin=y_coords, PAnodes=kinematic_pa,
                          inclin=inclination)

        bar.tremaine_weinberg(slit_width=slit_width)

        x_tw = bar.dfx_tw
        v_tw = bar.dfV_tw

        correlation = spearmanr(x_tw[np.isnan(x_tw) == False], v_tw[np.isnan(x_tw) == False])[0]

        if correlation < 0:
            print('Correlation still negative. Beware!')

    # Diagnostic plot

    # Set up colourbar limits for the image

    vmin_flux = np.nanpercentile(flux, 1)
    vmax_flux = np.nanpercentile(flux, 99)

    vmax_vel = 150
    vmin_vel = -150

    plt.figure(figsize=(12, 6))

    plt.subplot(1, 3, 1)

    plt.imshow(flux,
               cmap=cmocean.cm.thermal,
               origin='lower', interpolation='none',
               vmin=vmin_flux, vmax=vmax_flux)

    plt.axis('off')

    plt.subplot(1, 3, 2)

    plt.imshow(vel,
               cmap=cmocean.cm.balance,
               origin='lower', interpolation='none',
               vmin=vmin_vel, vmax=vmax_vel)

    # Also plot on the axis line

    rad = np.sqrt((vel.shape[0] - y_cen) ** 2 + (vel.shape[1] - x_cen) ** 2)
    ang = [0, np.pi] + np.radians(kinematic_pa)

    plt.plot(-rad * np.sin(ang) + vel.shape[1] / 2, rad * np.cos(ang) + vel.shape[0] / 2, 'k--',
             linewidth=2)

    plt.xlim([0, vel.shape[1]])
    plt.ylim([0, vel.shape[0]])

    plt.axis('off')

    plt.subplot(1, 3, 3)

    plt.imshow(bar.X_lon.reshape(vel.shape),
               cmap=cmocean.cm.balance,
               origin='lower', interpolation='none',
               vmin=vmin_vel,vmax=vmax_vel)

    plt.axis('off')

    plt.savefig('plots/pattern_speeds/' + galaxy + '_alma.png',
                bbox_inches='tight')
    plt.savefig('plots/pattern_speeds/' + galaxy + '_alma.pdf',
                bbox_inches='tight')

    plt.close()

    x_tw_err = bar.dfx_tw_err
    v_tw_err = bar.dfV_tw_err

    samples_filename = 'pattern_speeds_output/' + galaxy + '_samples_alma.hd5'

    # Initial guesses for m (2 km/s/arcsec) and c (which should be, by definition 0 km/s)

    initial_guesses = [2, 0]

    pos = initial_guesses + 1e-2 * np.random.randn(32, 2)

    n_walkers, n_dim = pos.shape
    n_steps = 5000

    if not os.path.exists(samples_filename) or overwrite_emcee_samples:

        try:
            os.remove(samples_filename)
        except FileNotFoundError:
            pass

        # Fit the pattern speed.

        backend = emcee.backends.HDFBackend(samples_filename)

        sampler = emcee.EnsembleSampler(n_walkers, n_dim, ln_prob, args=(x_tw, v_tw, v_tw_err), backend=backend)
        sampler.run_mcmc(pos, n_steps, progress=True, )

    else:

        sampler = emcee.backends.HDFBackend(samples_filename)

    # Plot the chains

    labels = ["m", "c"]

    fig, axes = plt.subplots(2, figsize=(8, 6), sharex=True)
    samples = sampler.get_chain()

    for i in range(n_dim):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel("Step Number")

    plt.savefig('plots/pattern_speeds/' + galaxy + '_steps_alma.png',
                bbox_inches='tight')
    plt.savefig('plots/pattern_speeds/' + galaxy + '_steps_alma.pdf',
                bbox_inches='tight')

    plt.close()

    # Flatten and thin samples, and plot the corner. Thin size should be about 1/3 the autocorrelation time, discard
    # the first half of samples as burn-in.

    tau = sampler.get_autocorr_time(quiet=True)
    thin = int(np.max(tau) / 3)

    flat_samples = sampler.get_chain(discard=int(n_steps / 2), thin=thin, flat=True)

    m_median = np.median(flat_samples[:, 0])
    c_median = np.median(flat_samples[:, 1])

    fig = corner.corner(flat_samples, labels=labels,
                        truths=[m_median, c_median], truth_color='k',
                        quantiles=[0.16, 0.84],
                        show_titles=True)

    plt.savefig('plots/pattern_speeds/' + galaxy + '_corner_alma.png',
                bbox_inches='tight')
    plt.savefig('plots/pattern_speeds/' + galaxy + '_corner_alma.pdf',
                bbox_inches='tight')

    plt.close()

    x_min, x_max = 1.2 * np.nanmin(x_tw), 1.2 * np.nanmax(x_tw)
    x_fit = np.linspace(x_min, x_max, 500)

    y_min, y_max = 1.2 * np.nanmin(v_tw), 1.2 * np.nanmax(v_tw)

    n_samples = 200

    y_samples = np.zeros([n_samples, len(x_fit)])

    for i in range(200):
        rand_i = np.random.randint(low=0, high=len(flat_samples[:, 0]))

        y_samples[i, :] = flat_samples[rand_i, 0] * x_fit + flat_samples[rand_i, 1]

    y_upper = np.percentile(y_samples, 84, axis=0)
    y_lower = np.percentile(y_samples, 16, axis=0)
    y_median = np.median(y_samples, axis=0)

    omega_bar = m_median.copy()
    omega_bar_err_up = np.percentile(flat_samples[:, 0], 84) - omega_bar
    omega_bar_err_down = omega_bar - np.percentile(flat_samples[:, 0], 16)

    # Convert to something usefully physical -- km/s/kpc, and get rid of the inclination

    km_s_kpc_conversion = (np.sin(inclination * np.pi / 180) * kpc_per_arcsec) ** -1

    omega_bar *= km_s_kpc_conversion
    omega_bar_err_up *= km_s_kpc_conversion
    omega_bar_err_down *= km_s_kpc_conversion

    # Write this pattern speed out

    np.savetxt('pattern_speeds_output/' + galaxy + '_pattern_speed_alma.txt',
               np.c_[omega_bar, omega_bar_err_up, omega_bar_err_down],
               header='omega_bar, omega_barr_err_up, omega_bar_err_down (all km/s/kpc)')

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_axes([0, 0, 1, 1])

    plt.errorbar(x_tw, v_tw,
                 xerr=x_tw_err, yerr=v_tw_err,
                 marker='o', c='r', linestyle='none')

    plt.plot(x_fit, y_median, c='k', lw=2)
    plt.fill_between(x_fit, y_upper, y_lower,
                     edgecolor='k', facecolor='k', alpha=0.25)

    plt.xlim([x_min, x_max])
    plt.ylim([y_min, y_max])

    plt.text(0.05, 0.95,
             r'$\Omega_p = %.2f^{+%.2f}_{-%.2f}\, \mathrm{km\,s}^{-1} \mathrm{kpc}^{-1}$ ' % (
                 omega_bar, omega_bar_err_up, omega_bar_err_down),
             transform=ax.transAxes)

    plt.xlabel(r'<$x$> $\left(^{\prime \prime}\right)$')
    plt.ylabel(r'<$v$> $\left(\mathrm{km\,s}^{-1}\right)$')

    plt.savefig('plots/pattern_speeds/' + galaxy + '_pattern_speed_tw_alma.png',
                bbox_inches='tight')
    plt.savefig('plots/pattern_speeds/' + galaxy + '_pattern_speed_tw_alma.pdf',
                bbox_inches='tight')

    plt.close()

print('Complete!')
