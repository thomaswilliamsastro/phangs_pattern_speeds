# -*- coding: utf-8 -*-
"""
Test effect of covering factor on recovered pattern speeds

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

sim_extensions = ['130', '55']

start_cf = 0.5
stop_cf = 3
step_cf = 0.25
covering_factors = np.arange(start_cf, stop_cf+step_cf, step_cf)

true_om_p = 43

for sim_extension in sim_extensions:

    plot_filename = plot_folder + 'spiral_models' + '/'

    plot_filename += sim_extension + '_cf_comparison'

    omega_bars = []
    omega_bars_err_up = []
    omega_bars_err_down = []

    covering_fractions = np.loadtxt(output_folder + 'spiral_models/covering_factors/' + sim_extension +
                                    '_fractions.txt')

    for covering_factor in covering_factors:

        filename = output_folder + 'spiral_models/covering_factors/' + sim_extension + \
                   '_cf_%.2f_pattern_speed.txt' % covering_factor

        omega_bar, omega_bar_err_up, omega_bar_err_down = np.loadtxt(filename, unpack=True)

        omega_bars.append(omega_bar)
        omega_bars_err_up.append(omega_bar_err_up)
        omega_bars_err_down.append(omega_bar_err_down)

    omega_bars = np.array(omega_bars)
    omega_bars_err_up = np.array(omega_bars_err_up)
    omega_bars_err_down = np.array(omega_bars_err_down)

    plt.figure(figsize=(4, 3))

    plt.errorbar(covering_fractions * 100, omega_bars, yerr=[omega_bars_err_down, omega_bars_err_up],
                 c='k', ls='none', marker='o')

    plt.axhline(true_om_p, c='k', ls='--')

    plt.ylabel(r'$\Omega_{p}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')
    plt.xlabel('Covering factor (%)')

    plt.tight_layout()

    # plt.show()

    plt.savefig(plot_filename + '.png',
                bbox_inches='tight')
    plt.savefig(plot_filename + '.pdf',
                bbox_inches='tight')

    plt.close()

print('Complete!')
