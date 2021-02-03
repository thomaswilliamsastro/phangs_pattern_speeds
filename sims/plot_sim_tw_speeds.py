# -*- coding: utf-8 -*-
"""
Plot the pattern speeds on simulated galaxies

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

sim_extensions = np.arange(0, 185, 5)
full_fields = False

plot_filename = plot_folder + 'spiral_models' + '/'

plot_filename += 'pattern_speeds_overview'

omega_bars = []
omega_bars_err_up = []
omega_bars_err_down = []

if full_fields:
    omega_bars_ff = []
    omega_bars_err_up_ff = []
    omega_bars_err_down_ff = []

for sim_extension in sim_extensions:

    sim_extension = str(sim_extension)
    filename = output_folder + 'spiral_models/' + sim_extension + '_pattern_speed.txt'

    try:

        omega_bar, omega_bar_err_up, omega_bar_err_down = np.loadtxt(filename, unpack=True)

        omega_bars.append(omega_bar)
        omega_bars_err_up.append(omega_bar_err_up)
        omega_bars_err_down.append(omega_bar_err_down)

    except:

        omega_bars.append(np.nan)
        omega_bars_err_up.append(np.nan)
        omega_bars_err_down.append(np.nan)

    if full_fields:
        filename = output_folder + 'spiral_models/full_fields_' + sim_extension + '_pattern_speed.txt'

        omega_bar, omega_bar_err_up, omega_bar_err_down = np.loadtxt(filename, unpack=True)

        omega_bars_ff.append(omega_bar)
        omega_bars_err_up_ff.append(omega_bar_err_up)
        omega_bars_err_down_ff.append(omega_bar_err_down)

plt.figure(figsize=(8, 6))

plt.errorbar(sim_extensions, omega_bars, yerr=[omega_bars_err_down, omega_bars_err_up],
             c='k', ls='none', marker='o')

if full_fields:
    plt.errorbar(sim_extensions, omega_bars_ff, yerr=[omega_bars_err_down_ff, omega_bars_err_up_ff],
                 c='r', ls='none', marker='o')

plt.ylabel(r'$\Omega_{p}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')
plt.xlabel('Spiral Orientation')

# plt.ylim([30, 40])

plt.axhline(43, c='k', ls='--')

plt.tight_layout()

plt.show()

plt.savefig(plot_filename + '.png',
            bbox_inches='tight')
plt.savefig(plot_filename + '.pdf',
            bbox_inches='tight')

plt.close()

print('Complete!')
