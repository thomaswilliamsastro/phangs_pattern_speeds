# -*- coding: utf-8 -*-
"""
Compare the density/material spiral models

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

plot_filename = plot_folder + 'density_material' + '/'

plot_filename += 'density_material_spiral_pattern_speeds_comparison'

bar_lens = ['115.0', '130.0', '150.0']
bar_ellips = ['0.25', '0.5']
bar_angs = ['5.0', '15.0', '30.0', '45.0', '60.0', '90.0', '120.0', '150.0', '180.0', '210.0']

colours = ['r', 'g']

bar_angs_label = [float(bar_ang) for bar_ang in bar_angs]

plt.figure(figsize=(8, 6))

for i, bar_len in enumerate(bar_lens):

    plt.subplot(len(bar_lens), 1, i + 1)

    for j, bar_ellip in enumerate(bar_ellips):

        omega_bars_material = []
        omega_bars_material_err_up = []
        omega_bars_material_err_down = []

        omega_bars_density = []
        omega_bars_density_err_up = []
        omega_bars_density_err_down = []

        for bar_ang in bar_angs:

            file_name = output_folder + 'density_material/material_' + bar_len + bar_ellip + bar_ang + \
                        '_pattern_speed.txt'

            try:
                omega_bar, omega_bar_err_up, omega_bar_err_down = np.loadtxt(file_name)
            except OSError:
                omega_bar, omega_bar_err_up, omega_bar_err_down = np.nan, np.nan, np.nan

            omega_bars_material.append(omega_bar)
            omega_bars_material_err_up.append(omega_bar_err_up)
            omega_bars_material_err_down.append(omega_bar_err_down)

            file_name = output_folder + 'density_material/density_' + bar_len + bar_ellip + bar_ang + \
                        '_pattern_speed.txt'

            try:
                omega_bar, omega_bar_err_up, omega_bar_err_down = np.loadtxt(file_name)
            except OSError:
                omega_bar, omega_bar_err_up, omega_bar_err_down = np.nan, np.nan, np.nan

            omega_bars_density.append(omega_bar)
            omega_bars_density_err_up.append(omega_bar_err_up)
            omega_bars_density_err_down.append(omega_bar_err_down)

        plt.errorbar(bar_angs_label, omega_bars_material,
                     yerr=[omega_bars_material_err_down, omega_bars_material_err_up],
                     c=colours[j], ls='none', marker='o', label='Ellipticity: %s' % bar_ellip)

        plt.errorbar(bar_angs_label, omega_bars_density,
                     yerr=[omega_bars_density_err_down, omega_bars_density_err_up],
                     c=colours[j], ls='none', marker='x')

        if bar_len == '115.0':
            plt.axhline(73, c='k', ls='--')
        elif bar_len == '130.0':
            plt.axhline(62, c='k', ls='--')
        elif bar_len == '150.0':
            plt.axhline(48, c='k', ls='--')

        plt.title('Bar len: %s arcsec' % bar_len[1:])

        if i == 1:
            plt.legend(loc='upper left', bbox_to_anchor=(1.05, 0.9), frameon=False)

        if i == len(bar_lens) - 1:

            plt.xlabel('Bar angle (deg)')

        plt.ylabel(r'$\Omega_\mathrm{p}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')

        plt.ylim([30, 75])

plt.tight_layout()

# plt.show()

plt.savefig(plot_filename + '.png',
            bbox_inches='tight')
plt.savefig(plot_filename + '.pdf',
            bbox_inches='tight')

plt.close()

print('Complete!')
