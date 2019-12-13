# -*- coding: utf-8 -*-
"""
Test to see how varying the slit width affects the measured pattern speed

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from alma.folders import phangs_folder, plot_folder, output_folder

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

galaxy = 'NGC2090'

start = 1
stop = 10
step = 0.5

slit_widths = np.arange(start, stop + step, step)

plt.figure(figsize=(8, 6))

omega_bars = np.zeros(len(slit_widths))
omega_bars_err_up = np.zeros_like(omega_bars)
omega_bars_err_down = np.zeros_like(omega_bars)

for i, slit_width in enumerate(slit_widths):

    file_name = output_folder + galaxy + '/slit_width/' + galaxy + '_bmask_sw_' \
                + str(slit_width) + '_pattern_speed_alma.txt'

    try:
        omega_bar, omega_bar_err_up, omega_bar_err_down = np.loadtxt(file_name, unpack=True)
    except OSError:
        print(file_name)
        omega_bar = np.nan
        omega_bar_err_down = np.nan
        omega_bar_err_up = np.nan

    omega_bars[i] = omega_bar
    omega_bars_err_up[i] = omega_bar_err_up
    omega_bars_err_down[i] = omega_bar_err_down

plt.errorbar(slit_widths, omega_bars,
             yerr=[omega_bars_err_down, omega_bars_err_up],
             fmt='o',
             c='k')

plt.xlabel(r'Slit Width ($^{\prime \prime}$)')
plt.ylabel(r'$\Omega_{p, \mathrm{TW}}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')

plt.tight_layout()

plt.savefig(plot_folder+galaxy+'/'+galaxy+'_alma_sw_comparison.png',
            bbox_inches='tight')
plt.savefig(plot_folder+galaxy+'/'+galaxy+'_alma_sw_comparison.pdf',
            bbox_inches='tight')

print('Complete!')
