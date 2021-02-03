# -*- coding: utf-8 -*-
"""
Test to see how varying the slit lengths affects the measured pattern speed

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from vars import phangs_folder, plot_folder, output_folder, slit_lengths

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

sim_extensions = ['15', '130', '160', '190', '1120', '1150', '1180']
sim_extensions = ['130']
full_fields = False

plot_filename = plot_folder + 'spiral_models' + '/'

plot_filename += 'pattern_speeds_overview'

# muse_galaxies = ['NGC1512']

for sim_extension in sim_extensions:

    omega_bars = np.zeros(len(slit_lengths))
    omega_bars_err_up = np.zeros_like(omega_bars)
    omega_bars_err_down = np.zeros_like(omega_bars)

    for i, slit_length in enumerate(slit_lengths):

        file_name = output_folder + 'spiral_models/slit_length/' + sim_extension + '_sl_' + str(slit_length) + \
                    '_pattern_speed.txt'
        omega_bar, omega_bar_err_up, omega_bar_err_down = np.loadtxt(file_name, unpack=True)

        omega_bars[i] = omega_bar
        omega_bars_err_up[i] = omega_bar_err_up
        omega_bars_err_down[i] = omega_bar_err_down

    max_r = np.load(output_folder + 'spiral_models/' + sim_extension + '_max_r.npy')

    plt.figure(figsize=(8, 6))

    plt.errorbar(slit_lengths, omega_bars,
                 yerr=[omega_bars_err_down, omega_bars_err_up],
                 fmt='o',
                 c='k')

    # Plot on the bar radius.

    plt.axvline(max_r, c='k', ls='--')

    plt.xlabel(r'$r$ ($^{\prime \prime}$)')
    plt.ylabel(r'$\Omega_{p}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')

    plt.tight_layout()

    plt.savefig(plot_filename + '.png',
                bbox_inches='tight')
    plt.savefig(plot_filename + '.pdf',
                bbox_inches='tight')

    plt.close()

print('Complete!')
