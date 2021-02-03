# -*- coding: utf-8 -*-
"""
Test to see how varying the slit width affects the measured pattern speed

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

from vars import phangs_folder, muse_version, muse_galaxies, muse_plot, muse_output, slit_widths, \
    hdu_types, star_masks, mask_outside_bars, output_folder, pattern_speed_version

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

# muse_galaxies = ['NGC1512']

if not os.path.exists(muse_plot + muse_version + '/slit_width/'):
    os.mkdir(muse_plot + muse_version + '/slit_width/')

pattern_speed_table = Table.read(output_folder + 'pattern_speed_table_' + pattern_speed_version + '.fits')

for galaxy in muse_galaxies:

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

                plot_filename = muse_plot + muse_version + '/slit_width/' + galaxy + '_' + hdu_type

                if star_mask:
                    plot_filename += '_smask'
                else:
                    plot_filename += '_nosmask'

                if mask_outside_bar:
                    plot_filename += '_bmask'
                else:
                    plot_filename += '_nobmask'

                plot_filename += '_muse_sw_comparison'

                omega_bars = np.zeros(len(slit_widths))
                omega_bars_err_up = np.zeros_like(omega_bars)
                omega_bars_err_down = np.zeros_like(omega_bars)

                for i, slit_width in enumerate(slit_widths):

                    file_name = muse_output + muse_version + '/slit_width/' + galaxy + '_' + hdu_type

                    if star_mask:
                        file_name += '_smask'
                    else:
                        file_name += '_nosmask'

                    if mask_outside_bar:
                        file_name += '_bmask'
                    else:
                        file_name += '_nobmask'

                    file_name += '_sw_' + str(slit_width) + '_pattern_speed_muse.txt'

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

                plt.figure(figsize=(4, 3))

                plt.errorbar(slit_widths, omega_bars,
                             yerr=[omega_bars_err_down, omega_bars_err_up],
                             fmt='o',
                             c='k')
                plt.axhline(om_p, c='k', ls='-')
                plt.axhline(om_p_err_up, c='k', ls='--')
                plt.axhline(om_p_err_down, c='k', ls='--')

                plt.xlabel(r'Slit Width ($^{\prime \prime}$)')
                plt.ylabel(r'$\Omega_\mathrm{P}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')

                plt.tight_layout()
                # plt.show()

                plt.savefig(plot_filename + '.png', bbox_inches='tight')
                plt.savefig(plot_filename + '.pdf', bbox_inches='tight')

                plt.close()

print('Complete!')
