# -*- coding: utf-8 -*-
"""
Plot pattern speeds from numerical sim

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from vars import phangs_folder, output_folder, plot_folder

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

galaxies = ['n628_70', 'n3351_02', 'n3351_20']
inclinations = ['09', '60']
emissions = ['gas', 'ns', 'os']

for galaxy in galaxies:

    snapshots = {'n628_70': ['0278'],
                 'n3351_02': ['0065', '0226', '0268'],
                 'n3351_20': ['0097', '0050']}[galaxy]

    for snapshot in snapshots:

        plot_filename = plot_folder + 'num_sim' + '/'
        plot_filename += galaxy + '_' + snapshot + '_pattern_speeds_overview'

        labels = []
        positions = []
        omega_bars = []
        omega_bars_err_up = []
        omega_bars_err_down = []

        i = 0

        for emission in emissions:

            for inclination in inclinations:

                filename = output_folder + 'num_sim/' + galaxy + '_' + snapshot + '_' + emission + '_i' + \
                           inclination + '_pattern_speed.txt'

                try:

                    omega_bar, omega_bar_err_up, omega_bar_err_down = np.loadtxt(filename, unpack=True)

                except IOError:
                    omega_bar, omega_bar_err_up, omega_bar_err_down = np.nan, np.nan, np.nan

                omega_bars.append(omega_bar)
                omega_bars_err_up.append(omega_bar_err_up)
                omega_bars_err_down.append(omega_bar_err_down)

                i += 1

                positions.append(i)
                labels.append('%s, %s' % (emission, inclination))

        plt.figure(figsize=(8, 6))

        plt.errorbar(positions, omega_bars, yerr=[omega_bars_err_down, omega_bars_err_up],
                     c='k', ls='none', marker='o')

        plt.xticks(positions, labels=labels, rotation=45)

        plt.ylabel(r'$\Omega_{p}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')
        plt.xlabel('Model')

        plt.ylim([0, 80])

        plt.tight_layout()

        # plt.show()

        plt.savefig(plot_filename + '.png',
                    bbox_inches='tight')
        plt.savefig(plot_filename + '.pdf',
                    bbox_inches='tight')

        plt.close()

print('Complete!')
