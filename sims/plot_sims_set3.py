# -*- coding: utf-8 -*-
"""
Compare the set3 sims

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

plot_filename = plot_folder + 'mock_set3/'

plot_filename += 'mock_set3_pattern_speeds_comparison'

vrots = ['9.8', '10.2', '10.6']
bar_angs = ['60.0', '210.0']
bar_ellips = ['0.25', '0.5']
bar_lens = ['10.0', '30.0', '50.0']
bar_len = '50.0'
spiral_arm_types = ['DW5', 'DW15', 'MAT']

colours = ['r', 'g']

bar_angs_label = [float(bar_ang) for bar_ang in bar_angs]

plt.figure(figsize=(8, 6))

for i, vrot in enumerate(vrots):

    plt.subplot(len(vrots), 1, i + 1)

    x_labels = []
    omega_bars = []
    omega_bars_err_up = []
    omega_bars_err_down = []

    for spiral_arm_type in spiral_arm_types:

        for bar_ang in bar_angs:

            for bar_ellip in bar_ellips:

                file_name = output_folder + 'mock_set3/' + vrot + bar_ang + bar_ellip + bar_len + spiral_arm_type + \
                            '_pattern_speed.txt'

                try:
                    omega_bar, omega_bar_err_up, omega_bar_err_down = np.loadtxt(file_name)
                except OSError:
                    omega_bar, omega_bar_err_up, omega_bar_err_down = np.nan, np.nan, np.nan

                x_labels.append(", ".join([spiral_arm_type, bar_ang, bar_ellip]))

                omega_bars.append(omega_bar)
                omega_bars_err_up.append(omega_bar_err_up)
                omega_bars_err_down.append(omega_bar_err_down)

    plt.errorbar(np.arange(len(omega_bars)), omega_bars, yerr=[omega_bars_err_down, omega_bars_err_up],
                 c='k', marker='o', ls='none')
    plt.title(vrot)

    plt.xticks([])

plt.xticks(np.arange(len(omega_bars)), x_labels, rotation=90)
plt.xlabel('spiral arm type, bar angle, bar ellipticity')

# plt.ylim([0, 30])

plt.tight_layout()
# plt.subplots_adjust(hspace=0)
plt.savefig(plot_filename + '.png', bbox_inches='tight')
plt.savefig(plot_filename + '.pdf', bbox_inches='tight')
plt.close()

print('Complete!')
