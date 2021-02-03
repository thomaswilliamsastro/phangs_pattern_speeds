# -*- coding: utf-8 -*-
"""
Plot bars onto rotation curves

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

from vars import galaxy_table, phangs_folder, plot_folder

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

rot_curve_table = Table.read('rotation_curves/RCtable_Nov2019.fits')
co_lum_table = Table.read('documents/galaxy_table.fits')

if not os.path.exists(os.path.join(plot_folder, 'r_bar_rot_curves')):
    os.makedirs(os.path.join(plot_folder, 'r_bar_rot_curves'))

for row in galaxy_table:

    galaxy = row['name']
    dist = row['dist']
    pa_galaxy = row['orient_ra']
    incl = row['orient_incl']

    r_bar = row['morph_bar_r']
    pa_bar = row['morph_bar_pa']

    m_star = row['props_mstar']
    m_hi = row['props_mhi']

    co_lum_row = co_lum_table[co_lum_table['NAME'] == galaxy.ljust(10)]
    m_h2 = co_lum_row['TDEP_PHANGS'] * co_lum_row['SFR_Z0MGS']

    if np.isnan(r_bar):
        continue

    # Deproject bar length, convert to kpc

    delta_pa = np.radians(pa_bar - pa_galaxy)
    incl = np.radians(incl)

    r_bar_deproj = r_bar * np.sqrt(np.cos(delta_pa) ** 2 + (np.sin(delta_pa) * (1 / np.cos(incl))) ** 2)
    r_bar_phys = r_bar_deproj * 1 / 3600 * np.pi / 180 * dist * 1e3

    # Read in the rotation curves

    rot_rows = rot_curve_table[rot_curve_table['Galaxy'] == galaxy]

    if len(rot_rows) == 0:
        continue

    gas_frac = (m_hi + m_h2) / (m_hi + m_h2 + m_star)

    plot_name = os.path.join(plot_folder, 'r_bar_rot_curves', galaxy)

    plt.figure()

    ax = plt.subplot(1, 1, 1)

    plt.errorbar(rot_rows['Radius'],
                 rot_rows['Vrot'], yerr=[rot_rows['Vrot_lower'], rot_rows['Vrot_upper']],
                 c='k', marker='o', ls='none')

    plt.axvline(r_bar_phys, c='k', ls='--')

    plt.xlabel('R (kpc)')
    plt.ylabel(r'v (km s$^{-1}$)')

    plt.text(0.95, 0.05, r'$\log10(M_\ast): %.2f$, $f_\mathregular{gas}$: %.2f' % (np.log10(m_star), gas_frac),
             va='bottom', ha='right',
             transform=ax.transAxes)

    plt.savefig(plot_name + '.png', bbox_inches='tight')

    plt.close()

print('Complete!')
