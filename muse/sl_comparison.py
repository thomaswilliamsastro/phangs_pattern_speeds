# -*- coding: utf-8 -*-
"""
Test to see how varying the slit lengths affects the measured pattern speed

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from vars import phangs_folder, muse_version, muse_galaxies, muse_plot, muse_output, slit_lengths, \
    mask_outside_bars, star_masks, hdu_types

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

if not os.path.exists(muse_plot):
    os.mkdir(muse_plot)
if not os.path.exists(muse_plot + muse_version):
    os.mkdir(muse_plot + muse_version)
if not os.path.exists(muse_plot + muse_version + '/slit_length/'):
    os.mkdir(muse_plot + muse_version + '/slit_length/')

bar_galaxy, bar_rs = np.loadtxt('environment/PHANGSmasks_v2.dat',
                                usecols=(0, 12),
                                unpack=True,
                                skiprows=4,
                                dtype=str)

# muse_galaxies = ['NGC1512']

for galaxy in muse_galaxies:

    for hdu_type in hdu_types:

        for star_mask in star_masks:

            for mask_outside_bar in mask_outside_bars:

                plot_filename = muse_plot + muse_version + '/slit_length/' + galaxy + '_' + hdu_type

                if star_mask:
                    plot_filename += '_smask'
                else:
                    plot_filename += '_nosmask'

                if mask_outside_bar:
                    plot_filename += '_bmask'
                else:
                    plot_filename += '_nobmask'

                plot_filename += '_muse_sl_comparison'

                omega_bars = np.zeros(len(slit_lengths))
                omega_bars_err_up = np.zeros_like(omega_bars)
                omega_bars_err_down = np.zeros_like(omega_bars)

                for i, slit_length in enumerate(slit_lengths):

                    file_name = muse_output + muse_version + '/slit_length/' + galaxy + '_' + hdu_type

                    if star_mask:
                        file_name += '_smask'
                    else:
                        file_name += '_nosmask'

                    if mask_outside_bar:
                        file_name += '_bmask'
                    else:
                        file_name += '_nobmask'

                    file_name += '_sl_' + str(slit_length) + '_pattern_speed_muse.txt'

                    try:
                        omega_bar, omega_bar_err_up, omega_bar_err_down = np.loadtxt(file_name, unpack=True)
                    except OSError:
                        omega_bar = np.nan
                        omega_bar_err_down = np.nan
                        omega_bar_err_up = np.nan

                    omega_bars[i] = omega_bar
                    omega_bars_err_up[i] = omega_bar_err_up
                    omega_bars_err_down[i] = omega_bar_err_down

                if np.all(np.isnan(omega_bars)):
                    print('Nothing found for %s: skipping' % galaxy)
                    continue

                plt.figure(figsize=(8, 6))

                plt.errorbar(slit_lengths, omega_bars,
                             yerr=[omega_bars_err_down, omega_bars_err_up],
                             fmt='o',
                             c='k')

                # Plot on the bar radius.

                idx = np.where(bar_galaxy == galaxy)
                if len(idx[0] > 0):

                    r = np.float(bar_rs[bar_galaxy == galaxy])

                    if r > 0:
                        plt.axvline(r, c='k', ls='--', lw=2)

                plt.xlabel(r'$r$ ($^{\prime \prime}$)')
                plt.ylabel(r'$\Omega_{p, \mathrm{TW}}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')

                plt.tight_layout()

                plt.savefig(plot_filename + '.png',
                            bbox_inches='tight')
                plt.savefig(plot_filename + '.pdf',
                            bbox_inches='tight')

                plt.close()

print('Complete!')
