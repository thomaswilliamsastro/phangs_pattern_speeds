# -*- coding: utf-8 -*-
"""
Plot the pattern speeds on simulated barred galaxies

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

bar_lens = ['15.0', '30.0']  #, '60.0']
bar_ellips = ['0.25', '0.5']
bar_angs = ['5.0', '15.0', '30.0', '45.0', '60.0', '90.0', '120.0', '150.0', '180.0', '210.0']
bar_angs_float = np.array(bar_angs, dtype=float)

plot_filename = plot_folder + 'spiral_models' + '/'

plot_filename += 'barred_pattern_speeds_overview'

c = ['r', 'g']

plt.figure(figsize=(6, 3 * len(bar_lens)))

for i, bar_len in enumerate(bar_lens):

    plt.subplot(len(bar_lens), 1, i + 1)

    for j, bar_ellip in enumerate(bar_ellips):

        omega_bars = []
        omega_bars_err_up = []
        omega_bars_err_down = []

        for bar_ang in bar_angs:
            filename = output_folder + 'spiral_models/' + bar_len + bar_ellip + bar_ang + '_pattern_speed.txt'

            omega_bar, omega_bar_err_up, omega_bar_err_down = np.loadtxt(filename, unpack=True)

            omega_bars.append(omega_bar)
            omega_bars_err_up.append(omega_bar_err_up)
            omega_bars_err_down.append(omega_bar_err_down)

        plt.errorbar(bar_angs, omega_bars, yerr=[omega_bars_err_down, omega_bars_err_up],
                     c=c[j], ls='none', marker='o', label='Ellipticity=' + bar_ellip)

    if i == 0:
        plt.legend(loc='upper right', frameon=False)

    plt.ylim([30, 55])

    plt.title('Bar length = %s arcsec' % bar_len)

    plt.xlabel(r'$\Delta \mathrm{PA}$ (deg)')
    plt.ylabel(r'$\Omega_{p}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')

plt.tight_layout()

plt.savefig(plot_filename + '.png',
            bbox_inches='tight')
plt.savefig(plot_filename + '.pdf',
            bbox_inches='tight')

print('Complete!')
