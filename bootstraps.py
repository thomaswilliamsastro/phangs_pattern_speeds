# Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import numpy as np
import pybar.pybar as pybar
from tqdm import tqdm

import ps_functions


def bootstrap_tw(flux, flux_err,
                 vel, vel_err,
                 grid_step=0.2, slit_width=1,
                 x_cen=None, y_cen=None,
                 centering_err=1,
                 pa=45, pa_err=1,
                 inclination=30,
                 dist=1,
                 n_bootstraps=1000,
                 bootstrap_filename='bootstraps.txt',
                 overwrite_bootstraps=False,
                 pattern_speed_filename='pattern_speeds.txt'
                 ):
    """Bootstrapped errors for the Tremaine-Weinberg pattern speed method.

    Bootstrap wrapper (bootstwrapper?) around pybar. For n_bootstraps, will perturb the position angle and centre by
    pa_err and centering_err, respectively. pybar will then do the TW integral using this setup, and fitting is done
    accounting for x- and y-errors using SciPy's ODR routines. Will finally produce a pattern speed and associated
    errors (16th and 84th percentiles), and save those to pattern_speed_filename.

    Args:
        flux (numpy.ndarray): 2D array representing the flux (e.g. stellar mass, H alpha flux).
        flux_err (numpy.ndarray): 2D array representing per-pixel flux error.
        vel (numpy.ndarray): 2D array representing velocity associated with the flux above (e.g. stellar velocity, H
            alpha velocity). Should be in units of km/s
        vel_err (numpy.ndarray): 2D array representing per-pixel velocity error. N.B. these four arrays should all be
            the same size!
        grid_step (float, optional): The pixel size of the image, in arcsec. Defaults to 0.2, which is for MUSE data.
        slit_width (float, optional): The width of each slit, in arcsec. Defaults to 1, which is the pybar default.
        x_cen (float, optional): The x-centre of the galaxy, in pixels. Defaults to None, which will use half way along
            the image.
        y_cen (float, optional): The y-centre of the galaxy, in pixels. Defaults to None, which will use half way along
            the image.
        centering_err (float, optional): The error in this centering, in pixels. Defaults to 1 pixel.
        pa (float, optional): The position angle of the galaxy (measured left of north), in degrees. Defaults to 45
            degrees.
        pa_err (float, optional): The error on this position angle, in degrees. Defaults to 1 degree
        inclination (float, optional): The galaxy inclination, in degrees. Defaults to 30 degrees.
        dist (float, optional): The distance to the galaxy, in Mpc. Defaults to 1Mpc.
        n_bootstraps (int, optional): The number of bootstraps to run. Defaults to 1000.
        bootstrap_filename (str, optional): Will save the calculated m and c values from the straight line fit to a text
            file. This can speed up later runs. Defaults to 'bootstraps.txt'
        overwrite_bootstraps (bool, optional): If False, will attempt to read in existing bootstraps and fit from there.
            Otherwise will fit for every bootstrap iteration regardless of whether fitting has been done before.
            Defaults to False.
        pattern_speed_filename (str, optional): Where to save the final pattern speed (and error) output. Defaults to
            'pattern_speeds.txt'.

    """

    # If centres not specified, use centre of image.

    if x_cen is None:
        x_cen = flux.shape[1] / 2

    if y_cen is None:
        y_cen = flux.shape[0] / 2

    # Calculate a physical distance conversion factor from the galaxy distance.

    kpc_per_arcsec = dist * 1e3 * np.sin(1 / 3600 * np.pi / 180)

    # Set up arrays for m and c. Load in any if applicable

    m_bootstrap = np.zeros(n_bootstraps)
    c_bootstrap = np.zeros_like(m_bootstrap)

    if not overwrite_bootstraps:

        try:
            m_loaded, c_loaded = np.loadtxt(bootstrap_filename,
                                            unpack=True)

            m_bootstrap[:len(m_loaded)] = m_loaded
            c_bootstrap[:len(c_loaded)] = c_loaded

        except OSError:
            pass

    for bootstrap_i in tqdm(range(n_bootstraps)):

        if m_bootstrap[bootstrap_i] != 0:
            # If we've already fitted here, just skip.

            continue

        pa_bootstrap = pa + np.random.normal(loc=0, scale=pa_err)

        # Swing this around if it passes through 0/360 degrees.

        if pa_bootstrap < 0:
            pa_bootstrap += 360
        if pa_bootstrap > 360:
            pa_bootstrap -= 360

        x_cen_bootstrap = x_cen + np.random.normal(loc=0, scale=centering_err)
        y_cen_bootstrap = y_cen + np.random.normal(loc=0, scale=centering_err)

        grid_shape = flux.shape
        grid_x = (np.arange(grid_shape[1]) - x_cen_bootstrap) * grid_step
        grid_y = (np.arange(grid_shape[0]) - y_cen_bootstrap) * grid_step

        x_coords, y_coords = np.meshgrid(grid_x, grid_y)

        # TODO: Replace with pydisc!

        bar = pybar.mybar(Flux=flux, Flux_err=flux_err,
                          Velocity=vel, Velocity_err=vel_err,
                          Xin=x_coords, Yin=y_coords, PAnodes=pa_bootstrap,
                          inclin=inclination)

        bar.tremaine_weinberg(slit_width=slit_width)

        x_tw = bar.dfx_tw
        v_tw = bar.dfV_tw

        x_tw_err = bar.dfx_tw_err
        v_tw_err = bar.dfV_tw_err

        m, m_err, c, c_err = ps_functions.odr_fit(x_tw, x_tw_err, v_tw, v_tw_err)

        m_bootstrap[bootstrap_i] = m
        c_bootstrap[bootstrap_i] = c

    # Now we've bootstrapped, pull out the pattern speed and associated errors.

    km_s_kpc_conversion = (np.sin(inclination * np.pi / 180) * kpc_per_arcsec) ** -1

    omega_bar = np.nanmedian(m_bootstrap)
    omega_bar_err_up = np.nanpercentile(m_bootstrap, 84) - omega_bar
    omega_bar_err_down = omega_bar - np.nanpercentile(m_bootstrap, 16)

    omega_bar *= km_s_kpc_conversion
    omega_bar_err_up *= km_s_kpc_conversion
    omega_bar_err_down *= km_s_kpc_conversion

    # Write this pattern speed out

    np.savetxt(pattern_speed_filename,
               np.c_[omega_bar, omega_bar_err_up, omega_bar_err_down],
               header='omega_bar, omega_barr_err_up, omega_bar_err_down (all km/s/kpc)')

    # Also save out each individual m and c

    np.savetxt(bootstrap_filename,
               np.c_[m_bootstrap, c_bootstrap],
               header='m, c')
