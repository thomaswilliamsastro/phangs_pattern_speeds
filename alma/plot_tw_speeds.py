# -*- coding: utf-8 -*-
"""
Plot the ALMA pattern speeds for all the PHANGS galaxies

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from astropy.io import fits
from astropy.table import Table
from matplotlib import cm

from alma.folders import phangs_folder, plot_folder, output_folder

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

# Read in the basic galaxy info

galaxy_table = fits.open('documents/phangs_sample_table_v1p4.fits')
galaxy_table = Table(galaxy_table[1].data)

galaxies = galaxy_table['NAME'][galaxy_table['HAS_ALMA'] == 1]

galaxies = sorted(galaxies)

plt.figure(figsize=(4, 0.5 * len(galaxies)))
ax1 = plt.subplot(1, 1, 1)

frame1 = plt.gca()

position = 0

colours = iter(cm.rainbow(np.linspace(0, 1, len(galaxies))))

omega_bars = []

for galaxy in galaxies:

    galaxy = galaxy.strip()

    c = next(colours)

    try:
        omega_bar, omega_bar_err_up, omega_bar_err_down = np.loadtxt(
            output_folder + galaxy + '_bmask_pattern_speed_alma.txt',
            unpack=True)

        omega_bars.append(omega_bar)

        plt.errorbar(omega_bar, position,
                     xerr=np.array([[omega_bar_err_down, omega_bar_err_up]]).T,
                     fmt='o',
                     c=c)
    except OSError:

        print(galaxy + ' not found')
        pass

    position += 0.25

plt.yticks(0.25 * np.array(range(len(galaxies))), galaxies)
plt.xlabel(r'$\Omega_{p, \mathrm{TW}}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')

plt.ylim([-0.125, 0.25 * len(galaxies) - 0.125])
plt.xlim([0, 100])

plt.tight_layout()

plt.savefig(plot_folder + 'pattern_speeds_overview_alma.png',
            bbox_inches='tight')
plt.savefig(plot_folder + 'pattern_speeds_overview_alma.pdf',
            bbox_inches='tight')

plt.close()

omega_bars = np.array(omega_bars)

# Also create a KDE plot, which is probably the best way to show all this data

# Filter out any particularly weird pattern speeds (i.e. negative, very high)

kde_mask = np.where((omega_bars < 100) & (omega_bars > 0))

plt.figure(figsize=(8, 6))

sns.kdeplot(omega_bars[kde_mask],
            bw='silverman', color='k', shade=True,
            clip=[0, 100])

plt.ylabel('Probability Density')
plt.xlabel(r'$\Omega_{p, \mathrm{TW}}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')

plt.xlim([0, 100])

plt.savefig(plot_folder + 'pattern_speeds_kde_alma.png',
            bbox_inches='tight')
plt.savefig(plot_folder + 'pattern_speeds_kde_alma.pdf',
            bbox_inches='tight')

print('Complete!')
