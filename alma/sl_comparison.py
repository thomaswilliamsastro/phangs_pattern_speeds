# -*- coding: utf-8 -*-
"""
Test to see how varying the slit length affects the measured pattern speed

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

from vars import phangs_folder, alma_version, alma_plot, alma_output, alma_galaxies, slit_lengths, mask_outside_bars, \
    galaxy_table, pattern_speed_version, output_folder

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

pattern_speed_table = Table.read(output_folder + 'pattern_speed_table_' + pattern_speed_version + '.fits')

for mask_outside_bar in mask_outside_bars:

    for galaxy in alma_galaxies:

        galaxy_row = galaxy_table[galaxy_table['name'] == galaxy]
        r_25 = galaxy_row['size_r25'][0]

        galaxy = galaxy.upper()

        pattern_speed_row = pattern_speed_table[pattern_speed_table['GALAXY'] == galaxy]

        try:
            om_p = pattern_speed_row['OM_P_ALMA'][0]
            om_p_err_up = om_p + pattern_speed_row['OM_P_ALMA_ERR_UP'][0]
            om_p_err_down = om_p - pattern_speed_row['OM_P_ALMA_ERR_DOWN'][0]
        except IndexError:
            om_p, om_p_err_up, om_p_err_down = np.nan, np.nan, np.nan

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

        # Read in the maximum slit length

        max_r = np.load(alma_output + alma_version + '/bmask/' + galaxy + '_bmask_max_r.npy')

        fig, ax = plt.subplots(figsize=(4, 3.5))

        plt.errorbar(slit_lengths, omega_bars,
                     yerr=[omega_bars_err_down, omega_bars_err_up],
                     fmt='o',
                     c='k')

        # Plot on the bar radius.

        bar_r = galaxy_table[galaxy_table['name'] == galaxy.lower()]['morph_bar_r'][0]
        plt.axvline(bar_r/2, c='k', ls='--', lw=2)

        plt.axvline(max_r, c='k', ls=':', lw=2)
        plt.axhline(om_p, c='k', ls='-')
        plt.axhline(om_p_err_up, c='k', ls='--')
        plt.axhline(om_p_err_down, c='k', ls='--')

        # Add a second axis normalised by r25

        ax_2 = ax.secondary_xaxis('top', functions=(lambda r: r / r_25, lambda r: r / r_25))
        ax_2.set_xlabel(r'$R/R_{25, \mathrm{opt}}$')

        plt.ylim([-10, 120])

        plt.xlabel(r'Slit length ($^{\prime \prime}$)')
        plt.ylabel(r'$\Omega_{p}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')

        plt.tight_layout()

        plot_name = alma_plot + alma_version + '/slit_length/' + galaxy

        if mask_outside_bar:
            plot_name += '_bmask'
        else:
            plot_name += '_nobmask'

        plot_name += '_alma_sl_comparison'

        # plt.show()

        plt.savefig(plot_name + '.png',
                    bbox_inches='tight')
        plt.savefig(plot_name + '.pdf',
                    bbox_inches='tight')

        plt.close()

print('Complete!')
