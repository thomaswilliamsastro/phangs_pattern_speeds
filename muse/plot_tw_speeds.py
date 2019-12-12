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

from muse.folders import phangs_folder, plot_folder, output_folder

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

galaxies = ['NGC1087', 'NGC1512', 'NGC1672', 'NGC3351', 'NGC4254', 'NGC5068',
            'IC5332', 'NGC1365', 'NGC1566', 'NGC2835', 'NGC3627', 'NGC4535', 'NGC628']

galaxies = sorted(galaxies)

plt.figure(figsize=(4, 0.5 * len(galaxies)))
ax1 = plt.subplot(1, 1, 1)

frame1 = plt.gca()

position = 0

colours = iter(cm.rainbow(np.linspace(0, 1, len(galaxies))))

for galaxy in galaxies:
    c = next(colours)

    try:
        omega_bar, omega_bar_err_up, omega_bar_err_down = np.loadtxt(
            output_folder + galaxy + '_mass_smask_bmask_pattern_speed_muse.txt',
            unpack=True)

        if omega_bar < 0:
            print(galaxy)

        plt.errorbar(omega_bar, position,
                     xerr=np.array([[omega_bar_err_down, omega_bar_err_up]]).T,
                     fmt='o',
                     c=c)
    except OSError:
        print(galaxy + ' not found!')
        pass

    position += 0.5

plt.yticks(0.5 * np.array(range(len(galaxies))), galaxies)
plt.xlabel(r'$\Omega_{p, \mathrm{TW}}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')

plt.ylim([-0.25, 0.5 * len(galaxies) - 0.25])
plt.xlim([0, 60])

plt.tight_layout()

plt.savefig(plot_folder + 'pattern_speeds_overview_muse.png',
            bbox_inches='tight')
plt.savefig(plot_folder + 'pattern_speeds_overview_muse.pdf',
            bbox_inches='tight')

print('Complete!')
