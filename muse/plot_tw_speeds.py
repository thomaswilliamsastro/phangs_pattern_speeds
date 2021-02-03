# -*- coding: utf-8 -*-
"""
Plot the MUSE pattern speeds for all the PHANGS galaxies

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

from vars import phangs_folder, muse_version, muse_galaxies, muse_plot, muse_output, hdu_types, star_masks, \
    mask_outside_bars

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

for hdu_type in hdu_types:

    for star_mask in star_masks:

        for mask_outside_bar in mask_outside_bars:

            plot_filename = muse_plot + muse_version + '/' + hdu_type

            if star_mask:
                plot_filename += '_smask'
            else:
                plot_filename += '_nosmask'

            if mask_outside_bar:
                plot_filename += '_bmask'
            else:
                plot_filename += '_nobmask'

            plot_filename += '_pattern_speeds_overview_muse'

            plt.figure(figsize=(4, 0.5 * len(muse_galaxies) + 1))

            position = 0

            colours = iter(cm.rainbow(np.linspace(0, 1, len(muse_galaxies))))

            for galaxy in muse_galaxies:
                c = next(colours)

                try:

                    muse_filename = muse_output + muse_version + '/' + hdu_type

                    if star_mask:
                        muse_filename += '_smask'
                    else:
                        muse_filename += '_nosmask'

                    if mask_outside_bar:
                        muse_filename += '_bmask/'
                    else:
                        muse_filename += '_nobmask/'

                    muse_filename += galaxy + '_' + hdu_type

                    if star_mask:
                        muse_filename += '_smask'
                    else:
                        muse_filename += '_nosmask'

                    if mask_outside_bar:
                        muse_filename += '_bmask'
                    else:
                        muse_filename += '_nobmask'

                    muse_filename += '_pattern_speed_muse.txt'

                    omega_bar, omega_bar_err_up, omega_bar_err_down = np.loadtxt(muse_filename, unpack=True)

                    if omega_bar < 0:
                        print(galaxy)

                    if hdu_type == 'toy_sim':
                        omega_bar -= 50

                    plt.errorbar(omega_bar, position,
                                 xerr=np.array([[omega_bar_err_down, omega_bar_err_up]]).T,
                                 fmt='o',
                                 c=c)
                except OSError:
                    print(galaxy + ' not found!')
                    pass

                position += 0.5

            if hdu_type == 'toy_sim':
                plt.axvline(0, c='k', ls='--', lw=2)

            if hdu_type == 'toy_sim':
                x_label = r'$\Omega_{p, \mathrm{obs}}-\Omega_{p, \mathrm{true}}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$'
            else:
                x_label = r'$\Omega_{p}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$'

            plt.yticks(0.5 * np.array(range(len(muse_galaxies))), muse_galaxies)
            plt.xlabel(x_label)

            plt.ylim([-0.25, 0.5 * len(muse_galaxies) - 0.25])

            if hdu_type == 'toy_sim':
                plt.xlim([-10, 10])
            else:
                plt.xlim([0, 60])

            plt.tight_layout()

            plt.savefig(plot_filename + '.png',
                        bbox_inches='tight')
            plt.savefig(plot_filename + '.pdf',
                        bbox_inches='tight')

            plt.close()

print('Complete!')
