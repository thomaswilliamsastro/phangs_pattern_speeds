# -*- coding: utf-8 -*-
"""
Test to see how varying the slit length affects the measured pattern speed

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from vars import phangs_folder, alma_version, alma_plot, alma_output, alma_galaxies, slit_lengths, mask_outside_bars

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

if not os.path.exists(alma_plot):
    os.mkdir(alma_plot)
if not os.path.exists(alma_plot + alma_version):
    os.mkdir(alma_plot + alma_version)
if not os.path.exists(alma_plot + alma_version + '/slit_length/'):
    os.mkdir(alma_plot + alma_version + '/slit_length/')

bar_galaxy, bar_rs = np.loadtxt('environment/PHANGSmasks_v2.dat',
                                usecols=(0, 12),
                                unpack=True,
                                skiprows=4,
                                dtype=str)

for mask_outside_bar in mask_outside_bars:

    for galaxy in alma_galaxies:

        galaxy = galaxy.strip()

        omega_bars = np.zeros(len(slit_lengths))
        omega_bars_err_up = np.zeros_like(omega_bars)
        omega_bars_err_down = np.zeros_like(omega_bars)

        for i, slit_length in enumerate(slit_lengths):

            file_name = galaxy

            if mask_outside_bar:
                file_name += '_bmask'
            else:
                file_name += '_nobmask'

            file_name += '_sl_' + str(slit_length) + '_pattern_speed_alma.txt'

            try:
                omega_bar, omega_bar_err_up, omega_bar_err_down = np.loadtxt(
                    alma_output + alma_version + '/slit_length/' + file_name, unpack=True)
            except OSError:
                # print(file_name)
                omega_bar = np.nan
                omega_bar_err_down = np.nan
                omega_bar_err_up = np.nan

            omega_bars[i] = omega_bar
            omega_bars_err_up[i] = omega_bar_err_up
            omega_bars_err_down[i] = omega_bar_err_down

        if np.all(np.isnan(omega_bars)):
            print('%s not found' % galaxy)
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

        plt.ylim([-10, 120])

        plt.xlabel(r'$r$ ($^{\prime \prime}$)')
        plt.ylabel(r'$\Omega_{p, \mathrm{TW}}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')

        plt.tight_layout()

        plot_name = alma_plot + alma_version + '/slit_length/' + galaxy

        if mask_outside_bar:
            plot_name += '_bmask'
        else:
            plot_name += '_nobmask'

        plot_name += '_alma_sl_comparison'

        plt.savefig(plot_name + '.png',
                    bbox_inches='tight')
        plt.savefig(plot_name + '.pdf',
                    bbox_inches='tight')

        plt.close()

print('Complete!')
