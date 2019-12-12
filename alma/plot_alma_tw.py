# -*- coding: utf-8 -*-
"""
Plot the ALMA pattern speeds for all the PHANGS galaxies

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from matplotlib import cm

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir('/Users/williams/Documents/phangs')

# Read in the basic galaxy info

galaxy_table = fits.open('documents/phangs_sample_table_v1p1.fits')
galaxy_table = Table(galaxy_table[1].data)

galaxies = galaxy_table['NAME'][galaxy_table['ALMA'] == 1]

galaxies = sorted(galaxies)

plt.figure(figsize=(4, 0.5 * len(galaxies)))
ax1 = plt.subplot(1, 1, 1)

frame1 = plt.gca()

position = 0

colours = iter(cm.rainbow(np.linspace(0, 1, len(galaxies))))

for galaxy in galaxies:

    galaxy = galaxy.strip()

    c = next(colours)

    try:
        omega_bar, omega_bar_err_up, omega_bar_err_down = np.loadtxt(
            'pattern_speeds_output/' + galaxy + '_pattern_speed_alma.txt',
            unpack=True)

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

plt.savefig('plots/pattern_speeds/pattern_speeds_overview_alma.png',
            bbox_inches='tight')
plt.savefig('plots/pattern_speeds/pattern_speeds_overview_alma.pdf',
            bbox_inches='tight')

print('Complete!')
