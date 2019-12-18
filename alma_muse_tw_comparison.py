# -*- coding: utf-8 -*-
"""
Compare MUSE/ALMA TW pattern speeds

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from muse.folders import phangs_folder

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

galaxies = ['NGC1087', 'NGC1512', 'NGC1672', 'NGC3351', 'NGC4254', 'NGC5068',
            'IC5332', 'NGC1365', 'NGC1566', 'NGC2835', 'NGC3627', 'NGC4535', 'NGC628']

galaxies = sorted(galaxies)

plt.figure(figsize=(4, 0.5 * len(galaxies) + 1))
ax1 = plt.subplot(1, 1, 1)

frame1 = plt.gca()

position = 0

# Plot on some dummy errorbars for labels

plt.errorbar(-100, 100, xerr=3,
             fmt='o', c='r', label='MUSE')
plt.errorbar(-100, 100, xerr=3,
             fmt='o', c='b', label='ALMA')

for galaxy in galaxies:

    try:

        # Load in the MUSE pattern speeds
        omega_bar, omega_bar_err_up, omega_bar_err_down = np.loadtxt(
            'pattern_speeds_output/muse/' + galaxy + '_mass_smask_bmask_pattern_speed_muse.txt',
            unpack=True)

        plt.errorbar(omega_bar, position - 0.125,
                     xerr=np.array([[omega_bar_err_down, omega_bar_err_up]]).T,
                     fmt='o',
                     c='r')

        # NGC0628 is inconsistently named >:(

        try:
            galaxy = {'NGC628': 'NGC0628',
                      }[galaxy]
        except KeyError:
            pass

        # Load in the ALMA pattern speeds
        omega_bar, omega_bar_err_up, omega_bar_err_down = np.loadtxt(
            'pattern_speeds_output/alma/' + galaxy + '_bmask_pattern_speed_alma.txt',
            unpack=True)

        plt.errorbar(omega_bar, position + 0.125,
                     xerr=np.array([[omega_bar_err_down, omega_bar_err_up]]).T,
                     fmt='o',
                     c='b')

    except OSError:

        print(galaxy)
        pass

    # Plot on horizontal lines to delineate galaxies

    plt.axhline(position + 0.25,
                c='k', ls='--')
    plt.axhline(position - 0.25,
                c='k', ls='--')

    position += 0.5

plt.yticks(0.5 * np.array(range(len(galaxies))), galaxies)
plt.xlabel(r'$\Omega_{p, \mathrm{TW}}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')

plt.legend(loc='upper right',
           frameon=False)

plt.ylim([-0.25, 0.5 * len(galaxies) + 0.5])
plt.xlim([0, 100])

ax1.yaxis.set_ticks_position('none')

plt.tight_layout()

plt.savefig('plots/pattern_speeds/pattern_speeds_comparison.png',
            bbox_inches='tight')
plt.savefig('plots/pattern_speeds/pattern_speeds_comparison.pdf',
            bbox_inches='tight')

print('Complete!')
