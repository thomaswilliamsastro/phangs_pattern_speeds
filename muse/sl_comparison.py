# -*- coding: utf-8 -*-
"""
Test to see how varying the slit lengths affects the measured pattern speed

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

from vars import phangs_folder, muse_version, muse_plot, muse_output, slit_lengths, \
    mask_outside_bars, star_masks, hdu_types, galaxy_table, output_folder, pattern_speed_version

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

pattern_speed_table = Table.read(output_folder + 'pattern_speed_table_' + pattern_speed_version + '.fits')

muse_galaxies = ['ngc3351']

for galaxy in muse_galaxies:

    galaxy_row = galaxy_table[galaxy_table['name'] == galaxy]
    r_25 = galaxy_row['size_r25'][0]
    bar_r = galaxy_row['morph_bar_r'][0]

    galaxy = galaxy.upper()

    pattern_speed_row = pattern_speed_table[pattern_speed_table['GALAXY'] == galaxy]

    for hdu_type in hdu_types:

        if hdu_type == 'mass':
            om_p_table_extension = 'MUSE_MASS'
        else:
            om_p_table_extension = 'MUSE_HA'

        try:
            om_p = pattern_speed_row['OM_P_'+om_p_table_extension][0]
            om_p_err_up = om_p + pattern_speed_row['OM_P_'+om_p_table_extension+'_ERR_UP'][0]
            om_p_err_down = om_p - pattern_speed_row['OM_P_'+om_p_table_extension+'_ERR_DOWN'][0]
        except IndexError:
            om_p, om_p_err_up, om_p_err_down = np.nan, np.nan, np.nan

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

                # Read in the maximum slit length

                max_r = np.load(muse_output + muse_version + '/' + hdu_type + '_smask_bmask/' +
                                galaxy + '_' + hdu_type + '_smask_bmask_max_r.npy')

                slit_lengths = np.array(slit_lengths)
                omega_bars = np.array(omega_bars)

                idx = np.where(slit_lengths <= max_r)
                slit_lengths = slit_lengths[idx]
                omega_bars = omega_bars[idx]
                omega_bars_err_down = omega_bars_err_down[idx]
                omega_bars_err_up = omega_bars_err_up[idx]

                fig, ax = plt.subplots(figsize=(4, 3.5))

                plt.errorbar(slit_lengths, omega_bars,
                             yerr=[omega_bars_err_down, omega_bars_err_up],
                             fmt='o',
                             c='k')

                # Plot on the bar radius.

                if bar_r > 0:
                    plt.axvline(bar_r, c='k', ls='--', lw=2)

                # plt.axvline(max_r, c='k', ls=':', lw=2)
                plt.axhline(om_p, c='k', ls='-')
                plt.axhline(om_p_err_up, c='k', ls='--')
                plt.axhline(om_p_err_down, c='k', ls='--')

                plt.xlim(0, max_r)

                # Add a second axis normalised by r25

                ax_2 = ax.secondary_xaxis('top', functions=(lambda r: r / r_25, lambda r: r / r_25))
                ax_2.set_xlabel(r'$R/R_{25, \mathrm{opt}}$')

                plt.xlabel(r'$r$ ($^{\prime \prime}$)')
                plt.ylabel(r'$\Omega_\mathrm{P}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')

                plt.tight_layout()

                # plt.show()
                plt.savefig(plot_filename + '.png', bbox_inches='tight')
                plt.savefig(plot_filename + '.pdf', bbox_inches='tight')

                plt.close()

print('Complete!')
