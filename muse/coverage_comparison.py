# -*- coding: utf-8 -*-
"""
Test on simulated data with covering fraction

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from vars import phangs_folder, muse_version, muse_galaxies, muse_plot, muse_output, star_masks, mask_outside_bars

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

if not os.path.exists(muse_plot + muse_version + '/covering_factors/'):
    os.mkdir(muse_plot + muse_version + '/covering_factors/')

# muse_galaxies = ['NGC1512']

hdu_type = 'toy_sim_cov'
cov_factor_sigmas = np.arange(1, 11)

for galaxy in muse_galaxies:

    for star_mask in star_masks:

        for mask_outside_bar in mask_outside_bars:

            plot_filename = muse_plot + muse_version + '/covering_factors/' + galaxy + '_' + hdu_type

            if star_mask:
                plot_filename += '_smask'
            else:
                plot_filename += '_nosmask'

            if mask_outside_bar:
                plot_filename += '_bmask'
            else:
                plot_filename += '_nobmask'

            plot_filename += '_muse_comparison'

            omega_bars = np.zeros(len(cov_factor_sigmas))
            omega_bars_err_up = np.zeros_like(omega_bars)
            omega_bars_err_down = np.zeros_like(omega_bars)

            for i, cov_factor_sigma in enumerate(cov_factor_sigmas):

                file_name = muse_output + muse_version + '/' + hdu_type

                if star_mask:
                    file_name += '_smask'
                else:
                    file_name += '_nosmask'

                if mask_outside_bar:
                    file_name += '_bmask'
                else:
                    file_name += '_nobmask'

                file_name += '/' + galaxy + '_' + hdu_type

                if star_mask:
                    file_name += '_smask'
                else:
                    file_name += '_nosmask'

                if mask_outside_bar:
                    file_name += '_bmask'
                else:
                    file_name += '_nobmask'

                file_name += '_cov_sig_' + str(cov_factor_sigma) + '_pattern_speed_muse.txt'

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

            plt.errorbar(cov_factor_sigmas, omega_bars - 50,
                         yerr=[omega_bars_err_down, omega_bars_err_up],
                         fmt='o',
                         c='k')

            plt.axhline(0, c='k', ls='--')

            plt.ylim([-10, 10])

            plt.xlabel(r'$\sigma$')
            plt.ylabel(r'$\Omega_{p, \mathrm{obs}}-\Omega_{p, \mathrm{true}}\, '
                       r'(\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')

            plt.tight_layout()

            plt.savefig(plot_filename + '.png',
                        bbox_inches='tight')
            plt.savefig(plot_filename + '.pdf',
                        bbox_inches='tight')

            plt.close()

print('Complete!')
